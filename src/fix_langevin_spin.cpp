/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "fix_langevin_spin.h"
#include <cmath>
#include <cstring>
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
// #include "random_park.h"
#include <mpi.h>
#include "random_mars.h"
#include "respa.h"
#include "update.h"
#include "utils.h"
#include "atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixLangevinSpin::FixLangevinSpin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_temp(NULL), random(NULL), ptr_T_el(NULL), emrd(NULL)//, emrd_all(NULL)
{
  if ((narg != 6) && (narg != 8))
    error->all(FLERR,"Illegal langevin/spin command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  nevery = 1;
  D = sigma = 0.0;

  temp = force->numeric(FLERR,arg[3]);
  alpha_t = force->numeric(FLERR,arg[4]);
  seed = force->inumeric(FLERR,arg[5]);

  if (alpha_t < 0.0) {
    error->all(FLERR,"Illegal langevin/spin command");
  } else if (alpha_t == 0.0) {
    tdamp_flag = 0;
  } else {
    tdamp_flag = 1;
  }

  // to be modified so that temp_flag = 0 if ttm

  if (temp < 0.0) {
    error->all(FLERR,"Illegal langevin/spin command");
  } else if (temp == 0.0) {
    temp_flag = 0;
  } else {
    temp_flag = 1;
  }

  int ttm_err = 1;
  int iarg = 6;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"lang")==0) {
      if (iarg+1 > narg-1) 
        error->all(FLERR,"Illegal fix langevin/spin/3tm command");
      strncpy(lang_name,arg[iarg+1],100);
      ttm_err = 0;
      ttm_flag = 1;
    }
    iarg++;
  }

  if (ttm_err == 1) error->all(FLERR,"Coupling to fix ttm not possible");

  // initialize Marsaglia RNG with processor-unique seed

  // random = new RanPark(lmp,seed + comm->me);
  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

FixLangevinSpin::~FixLangevinSpin()
{
  delete random;
  if (ttm_flag == 1) {
    memory->destroy(sigma_ttm);
    memory->destroy(emrd);
//         memory->destroy(emrd_all);
  }
}

/* ---------------------------------------------------------------------- */

int FixLangevinSpin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::init()
{
  // fix_langevin_spin has to be the last defined fix

  for (int whichfix = 0; whichfix < modify->nfix; whichfix++) {
    if (strcmp(lang_name,modify->fix[whichfix]->id) == 0) {
      id_lang = whichfix;
      break;
    }
    if (whichfix == (modify->nfix-1)) 
      error->universe_all(FLERR,"langevin fix ID for lattice is not defined");
  }

  int flag_force = 0;
  int flag_lang = 0;
  for (int i = 0; i < modify->nfix; i++) {
     if (strcmp("precession/spin",modify->fix[i]->style)==0) flag_force = MAX(flag_force,i);
     if (strcmp("langevin/spin",modify->fix[i]->style)==0) flag_lang = i;
  }
  if (flag_force >= flag_lang) error->all(FLERR,"Fix langevin/spin has to come after all other spin fixes");

  gil_factor = 1.0/(1.0+(alpha_t)*(alpha_t));
  dts = 0.25 * update->dt;
  double hbar = force->hplanck/MY_2PI;  // eV/(rad.THz)
  double kb = force->boltz;             // eV/K

  // compute sigma

  if (ttm_flag == 1) {
    nlocal_max = atom->nlocal;
    memory->grow(emrd,nlocal_max,"fix/langevin/spin:emrd");
    memory->grow(emrd_all,atom->natoms,"fix/langevin/spin:emrd_all");
    memory->grow(sigma_ttm,nlocal_max,"fix/langevin/spin:sigma_ttm");
    int tmp;
    ptr_T_el = (double **) modify->fix[id_lang]->extract("ptr_tforce",tmp);
    for (int i = 0; i < nlocal_max; i++) {
      D = (alpha_t*gil_factor*kb*(*ptr_T_el)[i]);
      D /= (hbar*dts);
      sigma_ttm[i] = sqrt(2.0*D); // to be checked
    }
  } else {
    D = (alpha_t*gil_factor*kb*temp);
    // D = (12.0/MY_2PI)*(MY_2PI*alpha_t*gil_factor*kb*temp);
    D /= (hbar*dts);
    sigma = sqrt(2.0*D);
  }

}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^respa")) {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  } else post_force(vflag);
  
  if (ttm_flag == 1) {
    int tmp;
    int nlocal = atom->nlocal; 
    ptr_T_el = (double **) modify->fix[id_lang]->extract("ptr_tforce",tmp);
    if (nlocal_max < nlocal) {    // grow emag lists if necessary
      nlocal_max = nlocal;
      memory->grow(emrd,nlocal_max,"fix/langevin/spin:emrd");
      // memory->grow(emrd_all,atom->natoms,"fix/langevin/spin:emrd_all");
      memory->grow(sigma_ttm,nlocal_max,"fix/langevin/spin:sigma_ttm");
    }  
  }

}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::add_tdamping(double spi[3], double fmi[3])
{
  double cpx = fmi[1]*spi[2] - fmi[2]*spi[1];
  double cpy = fmi[2]*spi[0] - fmi[0]*spi[2];
  double cpz = fmi[0]*spi[1] - fmi[1]*spi[0];

  // adding the transverse damping

  fmi[0] -= alpha_t*cpx;
  fmi[1] -= alpha_t*cpy;
  fmi[2] -= alpha_t*cpz;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::add_temperature(double fmi[3])
{

  // double rx = sigma*(2.0*random->uniform() - 1.0);
  // double ry = sigma*(2.0*random->uniform() - 1.0);
  // double rz = sigma*(2.0*random->uniform() - 1.0);
  double rx = sigma*random->gaussian();
  double ry = sigma*random->gaussian();
  double rz = sigma*random->gaussian();

  // adding the random field

  fmi[0] += rx;
  fmi[1] += ry;
  fmi[2] += rz;

  // adding gilbert's prefactor

  fmi[0] *= gil_factor;
  fmi[1] *= gil_factor;
  fmi[2] *= gil_factor;

}


/* ---------------------------------------------------------------------- */

void FixLangevinSpin::add_temperature_3tm(int i, double spi[3], double fmi[3])
{

//  double sigma[i] = f(Te[i]);
  int tmp;
//  double **ptr_T_el = (double **) modify->fix[id_lang]->extract("ptr_tforce",tmp);
  dts = 0.25 * update->dt;
  double hbar = force->hplanck/MY_2PI;  // eV/(rad.THz)
  double kb = force->boltz;  
  gil_factor = 1.0/(1.0+(alpha_t)*(alpha_t));
//  for (int j = 0; j < nlocal_max; j++){
      
  D = (alpha_t*gil_factor*kb*(*ptr_T_el)[i]);
  D /= (hbar*dts);
  sigma_ttm[i] = sqrt(2.0*D); // to be checked

//  }

  double rx = sigma_ttm[i]*random->gaussian();
  double ry = sigma_ttm[i]*random->gaussian();
  double rz = sigma_ttm[i]*random->gaussian();

  // compute random mag. energy

  emrd[i] = hbar*(rx*spi[0]+ry*spi[1]+rz*spi[2]);

  // adding the random field

  fmi[0] += rx;
  fmi[1] += ry;
  fmi[2] += rz;

  // adding gilbert's prefactor

  fmi[0] *= gil_factor;
  fmi[1] *= gil_factor;
  fmi[2] *= gil_factor;

}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
 *    extract energy of fluctuating langevin term (to be added to electrons)
 *    ------------------------------------------------------------------------- */

void *FixLangevinSpin::extract(const char *str, int &dim)
{
  if (strcmp(str,"emrd") == 0) {
    dim = 1;
  //  if (me==0) 
  //  MPI_Gather(&emrd, atom->nlocal, MPI_DOUBLE, emrd_all, atom->natoms, MPI_DOUBLE, 0, world);
    return &emrd;
  }
  if (strcmp(str,"flag") == 0) {
    dim = 0;
    return &ttm_flag;
  }
  return NULL;
}
