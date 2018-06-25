# LAMMPS ARTn 
The LAMMPS implication of Activation Relaxation Technique nouveau.


## Compiling
1. Put two files (min_artn.cpp min_artn.h) in the src directory.
2. Modify your Makefile.{machine} in LAMMPS to make sure that the MKL lib and include files can be linked. E.g. add the
   following lines in Makefile.{machine}
```
    PKG_INC = /opt/intel/compiler/composer_xe_2013_sp1.0.080/mkl/include
    PKG_HOM = /opt/intel/compiler/composer_xe_2013_sp1.0.080/mkl/lib/intel64
    PKG_LIB = -Wl,--start-group ${PKG_HOM}/libmkl_core.a ${PKG_HOM}/libmkl_sequential.a \
    ${PKG_HOM}/libmkl_intel_lp64.a -Wl,--end-group \
    /opt/intel/compiler/composer_xe_2013_sp1.0.080/compiler/lib/intel64/libiomp5.a \
    -lm -lpthread
```

3. Compile your LAMMPS codes.

## Benchmark 
After Compiling, you can run  prepared tests on  *./example* directory. One is Cu(110) surface diffusion events
searching,  and the other one is  thermally activated events searching  in Cu$_{64}$Zr$_{36}$ metallic glasses.

To launch a ARTn searching, two files are often needed.

  * in.lammps: use the following command to activate ARTn process.
    ```
    min_style artn
    minimize etol fotl maxiter maxeval
    ```

  * artn.control: Parameters that control ARTn simulation are defined in this file.

You can do the benchmark use the command that launch a LAMMPS simulation, like
```
 mpirun -np 4 lmp -in in.lammps 
```

## Output files
1. **log.event**: is expected to record all the important thermo info of events.
   Meaning of each column

   * Event: event index.
   * del-E: activation energy, or called energy barrier.
   * egv-sad: the lowest eigenvalue of saddle point
   * nsadl: number of atoms that are displaced more than *atom_disp_thr*(defined in **artn.control**) in saddle point
     configuration choosing the initial configuration as the reference configuration.
   * sad-dx: |x(sad) - x(minimum)|
   * sad-dy: |y(sad) - y(minimum)|
   * sad-dz: |z(sad) - z(minimum)|
   * sad-dr: |r(sad) - r(minimum)|
   * ref: id of reference configuration (initial minimum)
   * sad: id of reference configuration 
   * min: id of new minimum configuration (final minimum)
   * Center: id of center atom to do initial activation (perturbation)
   * Eref: energy of reference configuration (initial minimum)
   * Emin: energy of new minimum configuration (final minimum)
   * nMove: number of atoms that are displaced more than *atom_disp_thr*(defined in **artn.control**) in final minimum 
     configuration choose the initial configuration as the reference configuration.
   * pxx: xx component of pressure in final minimum configuration (set flag_press to 1 to show the pressure info) 
   * pyy: yy component of pressure in final minimum configuration 
   * pxx: zz component of pressure in final minimum configuration 
   * pxy: xy component of pressure in final minimum configuration 
   * pxz: xz component of pressure in final minimum configuration 
   * pyz: yz component of pressure in final minimum configuration 
   * Efinal: Energy of final minimum configuration 
   * status: 1, accept the new minimum configuration as the initial minimum configuration to start new search. 0, reject.
     Used in Metropolis condition. Set tempertature to a positive float number to activate this feature. Negative
     temperature will alway reject the new minimum configuration and start new search from the same initial minimum
     configuration.
   * dr: |r(final) - r(minimum)|

2. **log.artn** : is expected to record the  detailed info of ARTn simulation process and give a brief summary on the successful attempt rate and number of force evaluation.
   Some important parameters:

   * E-Eref: current potential energy.
   * ftot:  total force
   * fpar: force component along the eigenvector of lowest eigenvalue of Hessian matrix.
   * fperp: force component perpendicular to the eigenvector of lowest eigenvalue of Hessian matrix.
   * eigen: lowest eigenvalue of Hessian matrix.
   * evalf: number of force evaluation.
   * h.h': previous  eigenvector of lowest eigenvalue dot current eigenvector of lowest eigenvalue.

3. **sadl\_press.dat**: Recording the pressure of all the saddle point in the order of xx, yy, zz, xy, xz, yz.

4. **min.lammpstrj**: configurations of all the minimums. Dump atom format. TIMESTEP: id of the configuration (also used
   in log.event).

5. **sad.lammpstrj**: configurations of all the saddle points. Dump atom format. TIMESTEP: id of the configuration (also
   used in log.event).

## Parameters
All the parameters are defined in **artn.control**. Many parameters are self-explanatory. We review the others below.

Basic parameters (Standard ARTn process in the literature):

  * group_4_activat: The LAMMPS group ID of the atoms that will be activated.
  * cluster_radius: < 0, all atoms in the *group_4_activat* will be kicked in the initial pertubation in each ARTn loop; > 0, a cluster
    centers on a random atom  in *group_4_activat* with a radius of cluster_radius will be kicked in each ARTn loop. = 0, a random atom in *group_4_activat* will be kicked
    in each ARTn loop.
  * max_num_events: max number of events to be found.
  * force_th_saddle: Force threshold for convergence at saddle point.
  * flag_press: = 1, pressure of final minimum will be recorded in the **log.event**.
  * flag_sadl_press: = 1, pressure of saddle point configuration will be recorded in the **sadl_press.dat**

Extended parameters (Customized ARTn process):

  * events_per_atom: number of events that will be found on each atom. > 0, the max number of events will be the result of atom number in group_4_activat  times
    events_per_atom
  
  * flag_dump_direction: = 1, Define the initial kick direction from a dump file (name of this file is defined by
    fdump_direction).

## Contributor

Bin Xu, xubinrun@gmail.com

Lingti Kong, konglt@sjtu.edu.cn

## Reference

1. G. T. Barkema, et al., Phys. Rev. Lett. 77, 4358 (1996).
2. R. Malek, et al., Phys. Rev. E 62, 7723-7728 (2000).
3. E. Cances, et al., J. Chem. Phys. 130, 114711 (2009).
4. N. Mousseau, el al., J. At., Mol., Opt. Phys. 2012, 1 (2012).
    

