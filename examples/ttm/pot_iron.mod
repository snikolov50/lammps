#pair_style eam/fs
#pair_coeff * * FeVoter-ChenRecheck.fs Fe

#pair_style      hybrid/overlay snap zbl 4.0 4.8 spin/exchange/biquadratic 5.0
pair_style      hybrid/overlay eam/fs spin/exchange/biquadratic 5.0
#pair_coeff      * * snap Fe_pot_snappy.snapcoeff Fe_pot_snappy.snapparam Fe
#pair_coeff      * * zbl 26 26
pair_coeff      * * eam/fs FeVoter-ChenRecheck.fs Fe
pair_coeff      * * spin/exchange/biquadratic biquadratic 5.0 0.2827 -4.747 0.7810 -0.03619 -2.973 0.5273

# Setup neighbor style
neighbor 1.0    bin
neigh_modify    every 10 check yes delay 20

fix             AN all precession/spin cubic 0.001 0.0005 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
fix_modify      AN energy yes
fix 1z1         all nve/spin lattice yes
set             type 1 spin 2.2 0.0 0.0 1.0

# Setup output
thermo          1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
