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

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
------------------------------------------------------------------------- */

#include "fix_ttm.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "compute.h"
#include "iostream"
#include "utils.h"
#include <math.h>

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixTTM::FixTTM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fp(NULL), fpr(NULL), nsum(NULL), nsum_all(NULL),
  T_initial_set(NULL), gfactor1(NULL), gfactor2(NULL),
  flangevin(NULL), T_electron(NULL), T_electron_old(NULL), sum_vsq(NULL),
  sum_mass_vsq(NULL), sum_vsq_all(NULL), sum_mass_vsq_all(NULL),
  net_energy_transfer(NULL), net_energy_transfer_all(NULL), u_node(NULL), v_node(NULL), w_node(NULL), nvel(NULL),
  u_node_all(NULL), v_node_all(NULL), w_node_all(NULL), nvel_all(NULL), id_lang(NULL), id_spin(NULL), electronic_specific_heat(NULL)
{
  int arg_count = 12;
  if (narg < arg_count) error->all(FLERR,"Illegal fix ttm command");
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;
  biasflag = 0;
  id_temp = NULL;
  //  temperature = NULL;
  check_temp_flag = 0;
  int conv_err = 1;
  int lang_err = 1;
  int spin_err = 1;
//  int estop_err = 1;
  maxatom = 0;
  maxatom1 = 0;
  for (int index = 0; index < (narg-1); index++){

     if (strcmp(arg[index],"lang")==0){
        if (index+1 > narg-1) error->all(FLERR,"Illegal fix ttm command");
        strncpy(lang_fix_name,arg[index+1],100); 
        lang_err = 0;
     }

     if (strcmp(arg[index],"spin")==0){
        if (index+1 > narg-1) error->all(FLERR,"Illegal fix ttm command");
        strncpy(spin_fix_name,arg[index+1],100);
        spin_err = 0;
        arg_count = arg_count + 2;
     }

     if (strcmp(arg[index],"walls")==0){
        if (index+3 > narg-1) error->all(FLERR,"Illegal fix ttm command");
        walls[0] = force->inumeric(FLERR,arg[index+1]);
        walls[1] = force->inumeric(FLERR,arg[index+2]);
        walls[2] = force->inumeric(FLERR,arg[index+3]);
        arg_count = arg_count + 4;       
     }

     if (strcmp(arg[index],"conv")==0){
         conv_err = 0;
         if (index+2 > narg-1) error->all(FLERR,"Illegal fix ttm command");
         Nlimit = force->inumeric(FLERR,arg[index+2]);
         if (strcmp(arg[index+1],"yes") == 0) convflag = 1;
         else if (strcmp(arg[index+1],"no") == 0) convflag = 0;
         else error->all(FLERR,"Illegal fix ttm command");
     }

  }

//     if (strcmp(arg[index],"e_stop")==0){
//        if (index+1 > narg-1) {error->all(FLERR,"Improper e_stop option for fix ttm");}
//        if (strcmp(arg[index+1],"yes") == 0){
//            estopflag = 1;
//        }
//        else {
//            estopflag = 0;
//        }
//        estop_err = 0;
//     }

//  }
  
  if (lang_err == 1) error->all(FLERR,"No langevin fix name provided: coupling to fix ttm not possible");
  if (conv_err == 1) error->all(FLERR,"Convective option not specified");
  if (spin_err == 1) error->all(FLERR,"No langevin/spin fix name provided: coupling to fix ttm not possible");
//  if (estop_err == 1) {error->all(FLERR,"Electron stopping option not specified");}
  const char *filename = arg[3];
  FILE *fpr_2 = force->open_potential(arg[3]);
  if (fpr_2 == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open file %s",arg[3]);
    error->all(FLERR,str);
  }

  char linee[MAXLINE];
  double tresh_d;
  int tresh_i;
  specific_heat_flag = 0;
  // electronic_specific_heat
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%lg",&tresh_d);
  electronic_specific_heat_const = tresh_d;
  if (electronic_specific_heat_const==0){
	specific_heat_flag = 1;
  }
  // electronic_density 
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%lg",&tresh_d);
  electronic_density = tresh_d;

  // electronic_thermal_conductivity
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%lg",&tresh_d);
  electronic_thermal_conductivity = tresh_d;

  // gamma_p
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%lg",&tresh_d);
  gamma_p = tresh_d;

  // gamma_s
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%lg",&tresh_d);
  gamma_s = tresh_d;

  // v_0
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%lg",&tresh_d);
  v_0 = tresh_d;

  // nxnodes
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%ld",&tresh_i);
  nxnodes = tresh_i;

  // nynodes
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%ld",&tresh_i);
  nynodes = tresh_i;

  // nznodes
  utils::sfgets(FLERR,linee,MAXLINE,fpr_2,filename,error);
  sscanf(linee,"%ld",&tresh_i);
  nznodes = tresh_i;

  nodes_xyz[0] = nxnodes;
  nodes_xyz[1] = nynodes;
  nodes_xyz[2] = nznodes;

  fpr = fopen(arg[4],"r");
  if (fpr == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open file %s",arg[4]);
    error->one(FLERR,str);
  }

  nfileevery = force->inumeric(FLERR,arg[5]);
  strcpy(fname,arg[6]);

  if (nfileevery)
    if (narg != arg_count) error->all(FLERR,"Illegal fix ttm command");

  // error check

//  if (electronic_specific_heat <= 0.0)
//    error->all(FLERR,"Fix ttm electronic_specific_heat must be > 0.0");
  if (electronic_density <= 0.0)
    error->all(FLERR,"Fix ttm electronic_density must be > 0.0");
  if (electronic_thermal_conductivity < 0.0)
    error->all(FLERR,"Fix ttm electronic_thermal_conductivity must be >= 0.0");
  if (gamma_p <= 0.0) error->all(FLERR,"Fix ttm gamma_p must be > 0.0");
  if (gamma_s < 0.0) error->all(FLERR,"Fix ttm gamma_s must be >= 0.0");
  if (v_0 < 0.0) error->all(FLERR,"Fix ttm v_0 must be >= 0.0");
  if (nxnodes <= 0 || nynodes <= 0 || nznodes <= 0)
    error->all(FLERR,"Fix ttm number of nodes must be > 0");
  v_0_sq = v_0*v_0;

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];
  temperature_lang = NULL;
  // allocate 3d grid variables

  total_nnodes = nxnodes*nynodes*nznodes;

  memory->create(nsum,nxnodes,nynodes,nznodes,"ttm:nsum");
  memory->create(nsum_all,nxnodes,nynodes,nznodes,"ttm:nsum_all");
  memory->create(T_initial_set,nxnodes,nynodes,nznodes,"ttm:T_initial_set");
  memory->create(sum_vsq,nxnodes,nynodes,nznodes,"ttm:sum_vsq");
  memory->create(sum_mass_vsq,nxnodes,nynodes,nznodes,"ttm:sum_mass_vsq");
  memory->create(sum_vsq_all,nxnodes,nynodes,nznodes,"ttm:sum_vsq_all");
  memory->create(sum_mass_vsq_all,nxnodes,nynodes,nznodes,
                 "ttm:sum_mass_vsq_all");
  memory->create(T_electron_old,nxnodes,nynodes,nznodes,"ttm:T_electron_old");
  memory->create(T_electron,nxnodes,nynodes,nznodes,"ttm:T_electron");
  memory->create(u_node,nxnodes,nynodes,nznodes,"ttm:u_node");
  memory->create(v_node,nxnodes,nynodes,nznodes,"ttm:v_node");
  memory->create(w_node,nxnodes,nynodes,nznodes,"ttm:w_node");
  memory->create(nvel,nxnodes,nynodes,nznodes,"ttm:nvel"); 
  memory->create(u_node_all,nxnodes,nynodes,nznodes,"ttm:u_node_all");
  memory->create(v_node_all,nxnodes,nynodes,nznodes,"ttm:v_node_all");
  memory->create(w_node_all,nxnodes,nynodes,nznodes,"ttm:w_node_all");
  memory->create(nvel_all,nxnodes,nynodes,nznodes,"ttm:nvel_all");
  memory->create(net_energy_transfer,nxnodes,nynodes,nznodes,
                 "TTM:net_energy_transfer");
  memory->create(net_energy_transfer_all,nxnodes,nynodes,nznodes,
                 "TTM:net_energy_transfer_all");
  memory->create(electronic_specific_heat,nxnodes,nynodes,nznodes,"TTM:electronic_specific_heat");
  atom->add_callback(0);
  atom->add_callback(1);

  // set initial electron temperatures from user input file

  if (me == 1) read_initial_electron_temperatures();
  MPI_Bcast(&T_electron[0][0][0],total_nnodes,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

FixTTM::~FixTTM()
{
  if (nfileevery && me == 0) fclose(fp);

  delete [] gfactor1;
  delete [] gfactor2;

  memory->destroy(nsum);
  memory->destroy(nsum_all);
  memory->destroy(T_initial_set);
  memory->destroy(sum_vsq);
  memory->destroy(sum_mass_vsq);
  memory->destroy(sum_vsq_all);
  memory->destroy(sum_mass_vsq_all);
  memory->destroy(T_electron_old);
  memory->destroy(T_electron);
  memory->destroy(u_node);
  memory->destroy(v_node);
  memory->destroy(w_node);
  memory->destroy(nvel);
  memory->destroy(u_node_all);
  memory->destroy(v_node_all);
  memory->destroy(w_node_all);
  memory->destroy(nvel_all);
  memory->destroy(net_energy_transfer);
  memory->destroy(net_energy_transfer_all);
  memory->destroy(electronic_specific_heat);
}

/* ---------------------------------------------------------------------- */

int FixTTM::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTTM::init()
{

  for (int whichfix = 0; whichfix < modify->nfix; whichfix++) {
     if (strcmp(lang_fix_name,modify->fix[whichfix]->id) == 0){
         id_lang = whichfix;
         break;
     }
     if (whichfix == (modify->nfix-1)) error->universe_all(FLERR,"langevin fix ID is not defined");
  }

  for (int whichfix = 0; whichfix < modify->nfix; whichfix++) {
     if (strcmp(spin_fix_name,modify->fix[whichfix]->id) == 0){
         id_spin = whichfix;
         break;
     }
     if (whichfix == (modify->nfix-1)) error->universe_all(FLERR,"langevin fix ID is not defined");
  }

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix ttm with 2d simulation");
//  if (domain->nonperiodic != 0)
//    error->all(FLERR,"Cannot use non-periodic boundaries with fix ttm");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix ttm with triclinic box");

  for (int ixnode = 0; ixnode < nxnodes; ixnode++){
    for (int iynode = 0; iynode < nynodes; iynode++){
      for (int iznode = 0; iznode < nznodes; iznode++){
         net_energy_transfer_all[ixnode][iynode][iznode] = 0;
         u_node_all[ixnode][iynode][iznode] = 0;
         v_node_all[ixnode][iynode][iznode] = 0;
         w_node_all[ixnode][iynode][iznode] = 0;
         nvel_all[ixnode][iynode][iznode] = 0;
      }
    }
  }

  int tmp;
  int nlocal = atom->nlocal;
  double ***ptr_flangevin = (double ***) modify->fix[id_lang]->extract("flangevin",tmp);
  double **x = atom->x;

  if (atom->nmax > maxatom) {
         memory->destroy(*ptr_flangevin);
         maxatom = atom->nmax;
         memory->create(*ptr_flangevin,maxatom,3,"langevin:flangevin");
  }
  
//  if (specific_heat_flag==1){
     for (int i = 0; i < nlocal; i++) {
          (*ptr_flangevin)[i][0] = 0;
          (*ptr_flangevin)[i][1] = 0;
          (*ptr_flangevin)[i][2] = 0;
     }
   
     if (specific_heat_flag==1){
        for (int ixnode = 0; ixnode < nxnodes; ixnode++)
          for (int iynode = 0; iynode < nynodes; iynode++)
            for (int iznode = 0; iznode < nznodes; iznode++){
              electronic_specific_heat[ixnode][iynode][iznode] = 3*tanh(0.0002*T_electron[ixnode][iynode][iznode])*force->boltz;
            }
     }
 
  int nxnodes = nodes_xyz[0];
  int nynodes = nodes_xyz[1];
  int nznodes = nodes_xyz[2];

  int *mask = atom->mask;

  double **ptr_tforce = (double **) modify->fix[id_lang]->extract("ptr_tforce",tmp);
  if (atom->nmax > maxatom1) {
    maxatom1 = atom->nmax;
    memory->destroy(*ptr_tforce);
    memory->create(*ptr_tforce,maxatom1,"ttm:langevin:tforce");
  }

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      domain->x2lamda_remap(x[i], lamda);
      int ixnode = static_cast<int>(lamda[0]*nxnodes);
      int iynode = static_cast<int>(lamda[1]*nynodes);
      int iznode = static_cast<int>(lamda[2]*nznodes);

      (*ptr_tforce)[i] = T_electron[ixnode][iynode][iznode];
      if ((*ptr_tforce)[i] < 0.0)
        error->one(FLERR, "Fix ttm returned negative electron temperature");
    }
  }

}
/* ---------------------------------------------------------------------- */

void FixTTM::read_initial_electron_temperatures()
{
  char line[MAXLINE];

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_initial_set[ixnode][iynode][iznode] = 0;

  // read initial electron temperature values from file

  int ixnode,iynode,iznode;
  double T_tmp;
  while (1) {
    if (fgets(line,MAXLINE,fpr) == NULL) break;
    sscanf(line,"%d %d %d %lg",&ixnode,&iynode,&iznode,&T_tmp);
    if (T_tmp < 0.0)
      error->one(FLERR,"Fix ttm electron temperatures must be > 0.0");
    T_electron[ixnode][iynode][iznode] = T_tmp;
    T_initial_set[ixnode][iynode][iznode] = 1;
  }

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        if (T_initial_set[ixnode][iynode][iznode] == 0)
          error->one(FLERR,"Initial temperatures not all set in fix ttm");

  // close file

  fclose(fpr);

}

/* ---------------------------------------------------------------------- */

void FixTTM::setup(int vflag)
{
  post_force_setup(vflag);
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_setup(int /*vflag*/)
{
  int tmp;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double ***ptr_flangevin = (double ***) modify->fix[id_lang]->extract("flangevin",tmp);
  if (atom->nmax > maxatom) {
         memory->destroy(*ptr_flangevin);
         maxatom = atom->nmax;
         memory->create(*ptr_flangevin,maxatom,3,"langevin:flangevin");
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      atom->f[i][0] += (*ptr_flangevin)[i][0];
      atom->f[i][1] += (*ptr_flangevin)[i][1];
      atom->f[i][2] += (*ptr_flangevin)[i][2];
    }
  }
}
/* ---------------------------------------------------------------------- */

void FixTTM::pre_force(int /* vflag */)
{
  int tmp;
  int nxnodes = nodes_xyz[0];
  int nynodes = nodes_xyz[1];
  int nznodes = nodes_xyz[2];
  
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;              
 
  double **ptr_gfactor1 = (double **) modify->fix[id_lang]->extract("gfactor1",tmp); 
  double **ptr_gfactor2 = (double **) modify->fix[id_lang]->extract("gfactor2",tmp);
  double **ptr_tforce = (double **) modify->fix[id_lang]->extract("ptr_tforce",tmp);
  if (atom->nmax > maxatom1) {
    maxatom1 = atom->nmax;
    memory->destroy(*ptr_tforce);
    memory->create(*ptr_tforce,maxatom1,"ttm:langevin:tforce");
  }

  for (int i = 1; i <= atom->ntypes; i++) {
      (*ptr_gfactor1)[i] = - gamma_p / force->ftm2v;
      (*ptr_gfactor2)[i] = sqrt(45.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      domain->x2lamda_remap(x[i], lamda);
      int ixnode = static_cast<int>(lamda[0]*nxnodes);
      int iynode = static_cast<int>(lamda[1]*nynodes);
      int iznode = static_cast<int>(lamda[2]*nznodes);
 
      (*ptr_tforce)[i] = T_electron[ixnode][iynode][iznode];
      if ((*ptr_tforce)[i] < 0.0)
        error->one(FLERR, "Fix ttm returned negative electron temperature");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::end_of_step()
{

  int tmp;
  temperature_lang = (LAMMPS_NS::Compute *) modify->fix[id_lang]->extract("temperature",tmp);     
  if (temperature_lang && temperature_lang->tempbias){ //checking if temperature_lang is NULL (possible if no fix_modify is given)
     temperature_lang->compute_scalar();
  }
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;              

  double ***flangevin = (double ***) modify->fix[id_lang]->extract("flangevin",tmp);
  double **ptr_emrd = (double **) modify->fix[id_spin]->extract("emrd",tmp);
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++){
        net_energy_transfer[ixnode][iynode][iznode] = 0;
        nvel[ixnode][iynode][iznode] = 0;
        if (convflag==1){
           u_node[ixnode][iynode][iznode] = 0;
           v_node[ixnode][iynode][iznode] = 0;
           w_node[ixnode][iynode][iznode] = 0;
	}
      }

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {

      domain->x2lamda_remap(x[i], lamda);
      int ixnode = static_cast<int>(lamda[0]*nxnodes);
      int iynode = static_cast<int>(lamda[1]*nynodes);
      int iznode = static_cast<int>(lamda[2]*nznodes);

      if (temperature_lang && temperature_lang->tempbias) temperature_lang->remove_bias(i,v[i]);

      if (convflag==1){
         u_node[ixnode][iynode][iznode] += v[i][0];
         v_node[ixnode][iynode][iznode] += v[i][1];
         w_node[ixnode][iynode][iznode] += v[i][2];
         nvel[ixnode][iynode][iznode] += 1;
      }

      net_energy_transfer[ixnode][iynode][iznode] +=
        ((*flangevin)[i][0]*v[i][0] + (*flangevin)[i][1]*v[i][1] +
         (*flangevin)[i][2]*v[i][2] + (*ptr_emrd)[i]);
    }
  }

  MPI_Allreduce(&net_energy_transfer[0][0][0],&net_energy_transfer_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&u_node[0][0][0],&u_node_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&v_node[0][0][0],&v_node_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&w_node[0][0][0],&w_node_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&nvel[0][0][0],&nvel_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Comm_rank(world,&me);

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve
  double el_spec_heat_min = electronic_specific_heat_const;
  if (specific_heat_flag==1){
    el_spec_heat_min = 10000;
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++){
	    electronic_specific_heat[ixnode][iynode][iznode] = 3*tanh(0.0002*T_electron[ixnode][iynode][iznode])*force->boltz;
	    if (electronic_specific_heat[ixnode][iynode][iznode] < el_spec_heat_min) el_spec_heat_min = electronic_specific_heat[ixnode][iynode][iznode];

        }
  }

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;
  double stability_criterion = 1.0 -
    2.0*inner_dt/(el_spec_heat_min*electronic_density) *
    (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
  if (stability_criterion < 0.0) {
    inner_dt = 0.5*(el_spec_heat_min*electronic_density) /
      (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm",0);
  }

  for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps;
       ith_inner_timestep++) {

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          T_electron_old[ixnode][iynode][iznode] = T_electron[ixnode][iynode][iznode];

    // compute new electron T profile

    double Tc, Txr, Txl, Tyr, Tyl, Tzr, Tzl, u_vel, v_vel, w_vel, npar;
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++){

          Tc = T_electron_old[ixnode][iynode][iznode];
          u_vel = u_node_all[ixnode][iynode][iznode];
          v_vel = v_node_all[ixnode][iynode][iznode];
          w_vel = w_node_all[ixnode][iynode][iznode];
          npar  = nvel_all[ixnode][iynode][iznode];

          if (ixnode < nxnodes-1) Txr = T_electron_old[ixnode+1][iynode][iznode];
          else {
            if (domain->periodicity[0]) Txr = T_electron_old[0][iynode][iznode];
            else if (walls[0]==1) Txr = T_electron_old[0][iynode][iznode];
            else Txr = 0.0;
          }

          if (iynode < nynodes-1) Tyr = T_electron_old[ixnode][iynode+1][iznode];
          else {
            if (domain->periodicity[1]) Tyr = T_electron_old[ixnode][0][iznode];
            else if (walls[1]==1) Tyr = T_electron_old[ixnode][0][iznode];
            else Tyr = 0.0; 
          } 

          if (iznode < nznodes-1) Tzr = T_electron_old[ixnode][iynode][iznode+1];
          else {
            if (domain->periodicity[2]) Tzr = T_electron_old[ixnode][iynode][0];
            else if (walls[2]==1) Tzr = T_electron_old[ixnode][iynode][0];
            else Tzr = 0.0;
          }

          if (ixnode > 0) Txl = T_electron_old[ixnode-1][iynode][iznode];
          else {
            if (domain->periodicity[0]) Txl = T_electron_old[nxnodes-1][iynode][iznode];
            else if (walls[0]==1) Txl = T_electron_old[nxnodes-1][iynode][iznode];
            else Txl = 0.0; 
          } 

          if (iynode > 0) Tyl = T_electron_old[ixnode][iynode-1][iznode];
          else {
            if (domain->periodicity[1]) Tyl = T_electron_old[ixnode][nynodes-1][iznode];
            else if (walls[1]==1) Tyl = T_electron_old[ixnode][nynodes-1][iznode];
            else Tyl = 0.0;
          }

          if (iznode > 0) Tzl = T_electron_old[ixnode][iynode][iznode-1];
          else {
            if (domain->periodicity[2]) Tzl = T_electron_old[ixnode][iynode][nznodes-1];
            else if (walls[2]==1) Tzl = T_electron_old[ixnode][iynode][nznodes-1];
            else Tzl = 0.0;
          }

	if (specific_heat_flag==1){
          T_electron[ixnode][iynode][iznode] =
            Tc +
            inner_dt/(electronic_specific_heat[ixnode][iynode][iznode]*electronic_density) *
            (electronic_thermal_conductivity *
             ((Txr + Txl - 2*Tc)/dx/dx +
              (Tyr + Tyl - 2*Tc)/dy/dy +
              (Tzr + Tzl - 2*Tc)/dz/dz) -
             (net_energy_transfer_all[ixnode][iynode][iznode])/del_vol);
	}
	else {
          T_electron[ixnode][iynode][iznode] =
            Tc +
            inner_dt/(electronic_specific_heat_const*electronic_density) *
            (electronic_thermal_conductivity *
             ((Txr + Txl - 2*Tc)/dx/dx +
              (Tyr + Tyl - 2*Tc)/dy/dy +
              (Tzr + Tzl - 2*Tc)/dz/dz) -
             (net_energy_transfer_all[ixnode][iynode][iznode])/del_vol);
        }

          if (convflag == 1) {
            T_electron[ixnode][iynode][iznode] -=
              u_vel*Txr*(0.5*(inner_dt/dx)/npar) + u_vel*Txl*(0.5*(inner_dt/dx)/npar) -
              v_vel*Tyr*(0.5*(inner_dt/dy)/npar) + v_vel*Tyl*(0.5*(inner_dt/dy)/npar) -
              w_vel*Tzr*(0.5*(inner_dt/dz)/npar) + w_vel*Tzl*(0.5*(inner_dt/dz)/npar);
          } 
          
    }

  }

  // output nodal temperatures for current timestep
//  std::cout << specific_heat_flag << " " << electronic_specific_heat_const << std::endl;
  if ((nfileevery) && !(update->ntimestep % nfileevery)) {

    // compute atomic Ta for each grid point

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          nsum[ixnode][iynode][iznode] = 0;
          nsum_all[ixnode][iynode][iznode] = 0;
          sum_vsq[ixnode][iynode][iznode] = 0.0;
          sum_mass_vsq[ixnode][iynode][iznode] = 0.0;
          sum_vsq_all[ixnode][iynode][iznode] = 0.0;
          sum_mass_vsq_all[ixnode][iynode][iznode] = 0.0;
        }

    double massone;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        domain->x2lamda_remap(x[i], lamda);
        int ixnode = static_cast<int>(lamda[0]*nxnodes);
        int iynode = static_cast<int>(lamda[1]*nynodes);
        int iznode = static_cast<int>(lamda[2]*nznodes);
        double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        nsum[ixnode][iynode][iznode] += 1;
        sum_vsq[ixnode][iynode][iznode] += vsq;
        sum_mass_vsq[ixnode][iynode][iznode] += massone*vsq;
        if (temperature_lang && temperature_lang->tempbias) temperature_lang->restore_bias(i,v[i]);
      }

    MPI_Allreduce(&nsum[0][0][0],&nsum_all[0][0][0],total_nnodes,
                  MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&sum_vsq[0][0][0],&sum_vsq_all[0][0][0],total_nnodes,
                  MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&sum_mass_vsq[0][0][0],&sum_mass_vsq_all[0][0][0],
                  total_nnodes,MPI_DOUBLE,MPI_SUM,world);

    if (me == 0) {
      FILE *fp_1;
      std::string file_name(fname);
      file_name = file_name +  "." + std::to_string(update->ntimestep);
      char const *fnchar = file_name.c_str();
      fp_1 = fopen(fnchar,"w");
      if (fp_1 == NULL) {
         char str[128];
         snprintf(str,128,"Cannot open fix ttm file %s",fnchar);
         error->one(FLERR,str);
         fclose(fp_1);
      }
      std::string ts_str = "Timestep: " + std::to_string(update->ntimestep);
      fprintf(fp_1, ts_str.c_str());  //BIGINT_FORMAT,update->ntimestep);      
      double T_a;
      fprintf(fp_1,"\ny = columns\nz = rows");
      fprintf(fp_1,"\nAtomic temperatures\n\nx = 0\n");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++){
          if (ixnode > 0) fprintf(fp_1,"\n\nx = %d\n", ixnode);
        for (int iynode = 0; iynode < nynodes; iynode++){
            if (iynode > 0) fprintf(fp_1,"\n");
          for (int iznode = 0; iznode < nznodes; iznode++) {
            T_a = 0;
            if (nsum_all[ixnode][iynode][iznode] > 0)
              T_a = sum_mass_vsq_all[ixnode][iynode][iznode]/
                (3.0*force->boltz*nsum_all[ixnode][iynode][iznode]/force->mvv2e);
            fprintf(fp_1," %f",T_a);
          }
        }
      }

      fprintf(fp_1,"\n\nElectron temperatures\n\nx = 0\n");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++){
        if (ixnode > 0) fprintf(fp_1,"\n\nx = %d\n", ixnode);
        for (int iynode = 0; iynode < nynodes; iynode++){
            if (iynode > 0) fprintf(fp_1,"\n");
          for (int iznode = 0; iznode < nznodes; iznode++)
            fprintf(fp_1," %f ",T_electron[ixnode][iynode][iznode]);
         }
      }
       fprintf(fp_1,"\n");
       fclose(fp_1);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of 3d grid
------------------------------------------------------------------------- */

double FixTTM::memory_usage()
{
  double bytes = 0.0;
  bytes += 5*total_nnodes * sizeof(int);
  bytes += 14*total_nnodes * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
  return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTM::compute_vector(int n)
{
  double e_energy = 0.0;
  double transfer_energy = 0.0;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;
  int    grid_cells_low_count = 0;
  int    min_atoms_in_cell = 100000000;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        if (specific_heat_flag==1)
           e_energy += T_electron[ixnode][iynode][iznode]*electronic_specific_heat[ixnode][iynode][iznode]*electronic_density*del_vol;
	else
	   e_energy += T_electron[ixnode][iynode][iznode]*electronic_specific_heat_const*electronic_density*del_vol;

        transfer_energy += net_energy_transfer_all[ixnode][iynode][iznode]*update->dt;
        if (nvel_all[ixnode][iynode][iznode] < Nlimit){
           grid_cells_low_count += 1;
        }       
        if (nvel_all[ixnode][iynode][iznode] < min_atoms_in_cell){
           min_atoms_in_cell = nvel_all[ixnode][iynode][iznode];
        }
  }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  if (n == 2) return grid_cells_low_count;
  if (n == 3) return min_atoms_in_cell; 
  return 0.0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTM::write_restart(FILE *fp)
{
  double *rlist;
  memory->create(rlist,nxnodes*nynodes*nznodes,"TTM:rlist");

  int n = 0;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        rlist[n++] =  T_electron[ixnode][iynode][iznode];

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),n,fp);
  }

  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTM::restart(char *buf)
{
  int n = 0;
  double *rlist = (double *) buf;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_electron[ixnode][iynode][iznode] = rlist[n++];

}

// /* ---------------------------------------------------------------------- */

