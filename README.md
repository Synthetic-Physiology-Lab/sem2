# HOW TO CITE

Sandipan Chattaraj, Michelle Torre, Constanze Kalcher, Alexander Stukowski, Simone Morganti, Alessandro Reali, Francesco Pasqualini, “SEM2: Introducing mechanics in cell and tissue models using
coarse-grained homogeneous particle dynamics”, APL Bioengineering, 7, 046118 (2023) 1-14, https://doi.org/10.1063/5.0166829

# Instructions for Installation

## To compile

* Dowload LAMMPS version `stable_29Oct2020`
* Unpack and `cd` into `src` folder
* build LAMMPS by running `make serial` and `make mpi`
* In this repository: Unpack `USER-SEM.tar.gz` and cp `USER-SEM` into LAMMPS `src` directory
* Add the package to the makefile by appending `user-sem` to this line in the Makefile in the LAMMPS `src` directory:
        ```
        PACKUSER = user-adios user-atc user-awpmd user-bocs user-cgdna user-cgsdk user-colvars \
           user-diffraction user-dpd user-drude user-eff user-fep user-h5md \
           user-intel user-lb user-manifold user-meamc user-mesodpd user-mesont \
           user-mgpt user-misc user-mofff user-molfile \
           user-netcdf user-omp user-phonon user-plumed user-ptm user-qmmm \
           user-qtb user-quip user-reaction user-reaxc user-scafacos user-smd user-smtbq \
           user-sdpd user-sph user-tally user-uef user-vtk user-yaff user-sem
        ```

* Add the package by `make yes-user-sem`. It will throw a lot of warnings. Enter `y` to every open point.
* also add the molecule package `make yes-molecule`
* apply `make package-update`
* copy the content of the original `USER-SEM/src_changed` folder as downloaed from GitHub 
  to `src/USER-SEM/src_changed`
* run `make yes-user-sem` (again `y` to all) and `make package-update`
* recompile LAMMPS with `make serial` and `make mpi`
* either copy LAMMPS binaries `lmp_serial` and `lmp_mpi` to a directory included in `PATH` or add the `src` directory to `PATH`

## Ovito installation

* Download the ovito binaries from their website
* Add the `bin` directory to the path

### Ovito Python scripts

TODO: how to use? Undocumented?
  
# Instructions for Running Simulations

First: unpack the `example_simulations.tar.xz` folder.
**TODO**: actually, the simulations are in `USER-SEM/bench`. Navigate there.

## Creep experiments

1. Browse to the folder `creep/without_nucleus`. 
   From here, run the input file `in.growN1000.txt` using `mpirun -np 4 lmp_mpi -in in.growN1000.txt`.
   It will grow cell from 1 to Np = 1000 particles -> result stored as `cell_1000p.restart`.
   Serial alternative: `lmp_serial -in in.growN1000.txt`
2. Copy `cell_1000p.restart` in the folder `equilibration` and run the file `in.equil_1000.txt`.
   It will equilibrate the sample and store the result as `cell_1000p_equil_298K.restart`.
3. Copy `cell_1000p_equil_298K.restart` into the folder `creep_top_1` and run the input file `in.creepstress_1000.txt`.
   It will perform the creep experiment and store the trajectories of all particles, top and bottom layer particles in the files `all.xyz`, `top.xyz` and `bottom.xyz`.
   These cells can be viewed in visualization software such as Ovito.
4. You can repeat steps 1-3 for cells with nucleus by starting from the directory `example_simulations/with_nucleus`.

## Migration experiments

1. In `src` of LAMMPS, open the file `fix_sem_PMN.cpp`. In lines 172-173, ensure that `cpCell[i]` has the value `CP_Q`.
   Uncomment line 173 and comment line 172 as follows:
   ```
   //cpCell[i]=CP_P;
   cpCell[i]=CP_Q;
   ```
   **Please note that adding "//" before a line comments it out in C++**
2. If you have made a change, compile LAMMPS once again (see above, `make serial` and `make mpi`.
3. Then go to the migration folder in `migration/migration_3D_without_bias/tm_100`.
   Ensure that the restart file `cell_1000p.restart` is there in the folder.
   Run the input file `in.migration_3D_1000.txt` with similar commands as mentioned earlier.
   Note that migration simulations only run in serial mode.
   The output trajectories will be created in the `.xyz` files.
4. You can run simulations for various cases by browsing to the appropriate folder inside `migration` and running the input script.
   Check that the restart file is already present before running. 

## Proliferation experiments

1. In `src`, open the file `fix_sem_PMN.cpp`. In lines 172-173, ensure that `cpCell[i]` has the value `CP_P`.
   Uncomment line 172 and comment line 173 as follows:
   ```
   cpCell[i]=CP_P;
   //cpCell[i]=CP_Q;
   ```
   **Please note that adding "//" before a line comments it out in C++**
2. If you have made a change, compile LAMMPS once again (see above, `make serial` and `make mpi`.
3. Then go to the proliferation folder: `proliferation/prolif_upto_2_cells`. 
   Ensure that the restart file `grownUp.restart` is there in the folder.
   From here, run the input file `in.nuc_divPMN1000.txt` with similar commands as mentioned earlier.
   Note that proliferation simulations run in both serial and parallel mode.
   The output trajectories will be created in the `.xyz` files.
4. You can run unconstrained and constrained simulations by browsing to the appropriate folder 
   (`prolif_unconstrained_8_cells` and `prolif_constrained_8_cells`) inside `migration` 
   and running the input script.
   Check that the restart file is already present before running.
