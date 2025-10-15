# LAMMPS ARTn 
The LAMMPS implementation of Activation Relaxation Technique nouveau.

## What it includes

This package includes a LAMMPS implementation of the Activation Relaxation Technique nouveau (ARTn) for finding saddle points and transition pathways in atomic systems. It consists of:
- Two source files (`min_artn.cpp` and `min_artn.h`) that implement ARTn as a new minimization style in LAMMPS
- Example cases demonstrating applications to different material systems (Cu surface diffusion, Cu-Zr metallic glasses, 2D systems)
- Control parameter files and input scripts for running ARTn simulations


## Why this package

ARTn is a powerful open-ended saddle point search method that can systematically explore the potential energy landscape of materials without requiring prior knowledge of final states. This LAMMPS implementation provides several advantages:
- **Efficient saddle point searching**: Finds transition states and energy barriers for thermally activated processes
- **Scalable**: Leverages LAMMPS's parallel computing capabilities for large-scale atomistic simulations
- **Versatile**: Applicable to various material systems including surfaces, bulk materials, metallic glasses, and 2D systems
- **Automated event discovery**: Can systematically sample multiple transition pathways from a given configuration
- **Integration with LAMMPS**: Seamlessly works with LAMMPS's extensive potential library and existing analysis tools

## About ARTn

The **Activation-Relaxation Technique nouveau (ARTn)** is an open-ended saddle point search method originally developed by Mousseau et al. for exploring the potential energy surface (PES) of complex systems (see https://normandmousseau.com/en/research/art-nouveau for more details). Unlike traditional transition state search methods that require knowledge of both initial and final states, ARTn only requires the initial configuration to systematically discover transition pathways and energy barriers. 

We developed this package to make the entire ARTn process parallel and benefit from the utilities provided by LAMMPS C++ modules, which facilitates the implementation of additional features.

### Core ARTn Algorithm

The standard ARTn method works through the following stages:

1. **Initial Activation**: Random displacement is applied to selected atoms to push the system out of its local minimum
2. **Basin Escape**: The system moves away from the initial basin using the lowest curvature direction computed via the Lanczos algorithm
3. **Saddle Point Convergence**: Once the eigenvalue becomes negative, the system converges to the nearest saddle point by minimizing forces perpendicular to the unstable mode
4. **Relaxation**: After crossing the saddle point, the system relaxes to a new local minimum

### This Implementation

This LAMMPS implementation follows the standard ARTn procedure while providing several extended features:

**Standard ARTn Features:**
- Lanczos algorithm for efficient computation of the lowest eigenvalue and eigenvector
- Configurable activation strategies (single atom, cluster, or all atoms)
- Force-based convergence criteria for saddle points
- Automated push-over to find connected minima

**Extended Features:**
- **Events-per-atom mode**: Systematically search for multiple events on each atom in the group
- **Custom initial directions**: Define initial kick directions from dump files, delta direction files, or deformation gradients (useful for studying specific mechanisms)
- **Validation mechanisms**: Push-back verification to confirm saddle points connect to the initial minimum
- **Pressure monitoring**: Track stress tensor components at saddle points and final minima (valuable for studying pressure-dependent processes)
- **Metropolis sampling**: Optional acceptance/rejection of new minima based on energy barriers and temperature
- **FIRE algorithm integration**: Alternative minimization method for improved efficiency in certain systems
- **Comprehensive output**: Detailed logging of eigenvalues, forces, displacements, and atomic configurations at each stage

These extensions make this implementation particularly suitable for studying complex materials phenomena such as defect migration in crystals, structural relaxations in amorphous materials, and surface diffusion processes under varying conditions.



## Installation

### Prerequisites

Intel Math Kernel Library (MKL) package.

### Compile

Command for make: (Makefile example lammps/src/MAKE/OPTIONS/Makefle.oneapi)

```bash
cd ${software}/ARTn
cp -r src/USER-ARTn ${software}/lammps/src
cd ${software}/lammps/src
make yes-USER-ARTn
make machine # need MKL
```

Command for cmake: 
```bash
cd ${software}/ARTn
cp -r src/USER-ARTn ${software}/lammps/src
cp src/USER-ARTn.cmake ${software}/lammps/cmake/Modules/Packages/

cd ${software}/lammps/cmake
sed -i '/foreach(PKG_WITH_INCL / s/)/ USER-ARTn)/' CMakeLists.txt
sed -i '/set(STANDARD_PACKAGES/,/)/ s/)/  \n  USER-ARTn)/' CMakeLists.txt

cd ..; mkdir build; cd build
# module load mkl mpi, one may need those environment
cmake -D PKG_USER-ARTn=on -C ../cmake/presets/basic.cmake ../cmake
cmake --build .
make install
```


## Usage and Examples

After compiling, you can run prepared tests in the *./example* directory. One is Cu(110) surface diffusion events
searching, and the second one is thermally activated events searching in Cu$_{64}$Zr$_{36}$ metallic glasses.
The third one is for testing new features. You don't need it.

To launch an ARTn search, two files are often needed:

  * **in.artn** (or in.lammps): Use the following commands to activate the ARTn process.
    ```
    min_style artn
    minimize etol ftol maxiter maxeval
    ```
    where `etol` is energy tolerance, `ftol` is force tolerance, `maxiter` is maximum iterations, and `maxeval` is maximum force evaluations.

  * **artn.control**: Parameters that control ARTn simulation are defined in this file.

You can run the benchmark using the command that launches a LAMMPS simulation, like:
```bash
mpirun -np 4 lmp -in in.artn
```
(where `lmp` is your LAMMPS executable, which may be named `lmp_mpi`, `lmp_serial`, etc.)

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
   * sad: id of saddle point configuration 
   * min: id of new minimum configuration (final minimum)
   * Center: id of center atom to do initial activation (perturbation)
   * Eref: energy of reference configuration (initial minimum)
   * Emin: energy of new minimum configuration (final minimum)
   * nMove: number of atoms that are displaced more than *atom_disp_thr*(defined in **artn.control**) in final minimum 
     configuration choose the initial configuration as the reference configuration.
   * pxx: xx component of pressure in final minimum configuration (set flag_press to 1 to show the pressure info) 
   * pyy: yy component of pressure in final minimum configuration 
   * pzz: zz component of pressure in final minimum configuration 
   * pxy: xy component of pressure in final minimum configuration 
   * pxz: xz component of pressure in final minimum configuration 
   * pyz: yz component of pressure in final minimum configuration 
   * Efinal: Energy of final minimum configuration 
   * status: 1, accept the new minimum configuration as the initial minimum configuration to start new search. 0, reject.
     Used in Metropolis condition. Set temperature to a positive float number to activate this feature. Negative
     temperature will always reject the new minimum configuration and start new search from the same initial minimum
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

3. **sadl_press.dat**: Recording the pressure of all the saddle points in the order of xx, yy, zz, xy, xz, yz.

4. **min.lammpstrj**: configurations of all the minimums. Dump atom format. TIMESTEP: id of the configuration (also used
   in log.event).

5. **sad.lammpstrj**: configurations of all the saddle points. Dump atom format. TIMESTEP: id of the configuration (also
   used in log.event).

## Parameters
All the parameters are defined in **artn.control**. Many parameters are self-explanatory. We review the others below.

Basic parameters (Standard ARTn process in the literature):

  * group_4_activat: The LAMMPS group ID of the atoms that will be activated.
  * cluster_radius: < 0, all atoms in the *group_4_activat* will be kicked in the initial perturbation in each ARTn loop; > 0, a cluster
    centered on a random atom in *group_4_activat* with a radius of cluster_radius will be kicked in each ARTn loop. = 0, a random atom in *group_4_activat* will be kicked
    in each ARTn loop.
  * max_num_events: max number of events to be found.
  * force_th_saddle: Force threshold for convergence at saddle point.
  * flag_press: = 1, pressure of final minimum will be recorded in the **log.event**.
  * flag_sadl_press: = 1, pressure of saddle point configuration will be recorded in the **sadl_press.dat**

Extended parameters (Customized ARTn process):

  * events_per_atom: Number of events that will be found on each atom. If > 0, the maximum number of events will be 
    the number of atoms in group_4_activat multiplied by events_per_atom.
  
  * flag_dump_direction: = 1, define the initial kick direction from a dump file (the name of this file is defined by
    fdump_direction).

## Contributors

Bin Xu, xubinrun@gmail.com, xubin@nimte.ac.cn

Lingti Kong, konglt@sjtu.edu.cn

## Reference

If you find this package useful, please cite our paper and the original ARTn method papers.

**This Implementation:**

1. B. Xu, et al., Phys. Rev. Lett. 120, 125503 (2018).

**ARTn Method Papers:**

1. G. T. Barkema, et al., Phys. Rev. Lett. 77, 4358 (1996).
2. R. Malek, et al., Phys. Rev. E 62, 7723-7728 (2000).
3. E. Cances, et al., J. Chem. Phys. 130, 114711 (2009).
4. N. Mousseau, et al., J. At., Mol., Opt. Phys. 2012, 1 (2012).
    

