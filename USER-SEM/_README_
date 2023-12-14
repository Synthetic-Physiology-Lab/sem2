--------------------------
HOW TO CITE:
--------------------------
Chattaraj Sandipan, Pasqualini Francesco, SEM2: A computational framework for modeling cell and tissue mechanics with coarse-grained subcellular elements.


--------------------------
To compile: (once LAMMPS compiles)
--------------------------
- put USER-SEM in src
- add the package to the makefile by appending user-sem to this line in the Makefile:

	PACKUSER = user-adios user-atc user-awpmd user-bocs user-cgdna user-cgsdk user-colvars \
	   user-diffraction user-dpd user-drude user-eff user-fep user-h5md \
	   user-intel user-lb user-manifold user-meamc user-mesodpd user-mesont \
	   user-mgpt user-misc user-mofff user-molfile \
	   user-netcdf user-omp user-phonon user-plumed user-ptm user-qmmm \
	   user-qtb user-quip user-reaction user-reaxc user-scafacos user-smd user-smtbq \
	   user-sdpd user-sph user-tally user-uef user-vtk user-yaff user-sem

copy variable.cpp from the USER-SEM/src_changed folder to src

- add the package
	$make yes-user-sem

- also add the molecule package
	$make yes-molecule

- if you change the code in the package, to recompile, you have to update the package:
	$make package-update

- recompile LAMMPS from the src folder with
  make serial
  and
  make mpi
  

- test with bench/in.*.txt
-> in.growN1000.txt: grow cell from 1 to N = 1000 SCE -> result stored as grownUp.restart
-> in.randN1000.txt: equilibrate cell with N = 1000 randomly positioned SCE -> result stored as grownUp.restart
-> in.stress15000.txt: using grownUp.restart produced above, fix top/bottom and stretch with 15 Pa
-> in.nuc_growN1000.txt: like in.growN1000.txt but with addition of nucleus (5% of cell volume)
-> in.nuc_stress15000.txt: stretching result of in.nuc_growN1000.txt with 15 Pa
=> xyz files can be visualized in VMD (see below for how to deal with growth)
-> these scripts are just for demonstration!
=> in practice, we run LAMMPS as a library and compute constants at runtime instead of hard-coding them as in these demo-scripts

Creep experiments:
1) Run the input file in.growN1000.txt: It will grow cell from 1 to Np = 1000 particles -> result stored as cell_1000p.restart
2) Copy cell_1000p.restart in the folder 'equilibration' and run the file in.equil_1000.txt: It will equilibrate the sample and store the result as cell_1000p_equil_298K.restart
3) Copy cell_1000p_equil_298K.restart into the folder 'creep_top_1'and run the input file in.creepstress_1000.txt: It will perform the creep experiment and store the trajectories of all particles, top and bottom layer particles in the files all.xyz, top.xyz and bottom.xyz These cells can be viewed in visualization software such as Ovito.

Migration experiments:
1) In src, open the file fix_sem_PMN.cpp. In lines 172-173, ensure that cpCell[i] has the value CP_Q. Uncomment line 173 and comment line 172. as follows:
//cpCell[i]=CP_P;
cpCell[i]=CP_Q;
2) If you have made a change, compile LAMMPS once again by running the following commands one by one from src:
make serial
make mpi
3) Then go to the migration folder in bench and run the input file in.migration_3D_1000.txt. Ensure that the restart file cell_1000p.restart is there in the folder from where you are running the input script. The output trajectories will be created in the .xyz files

Proliferation experiments:
1) In src, open the file fix_sem_PMN.cpp. In lines 172-173, ensure that cpCell[i] has the value CP_P. Uncomment line 172 and comment line 173. as follows:
cpCell[i]=CP_P;
//cpCell[i]=CP_Q;
2) If you have made a change, compile LAMMPS once again by running the following commands one by one from src:
make serial
make mpi
2) If you have made a change, compile LAMMPS once again by running the following commands one by one from src:
make serial
make mpi
3) Then go to the proliferation folder in bench and run the input file in.nuc_divPMN1000.txt Ensure that the restart file grownUp.restart is there in the folder from where you are running the input script. The output trajectories will be created in the .xyz files
