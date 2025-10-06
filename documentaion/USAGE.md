# Quick reference guide

This file gives short example `input.ini` snippets for common use-cases. You may edit all the sections and veriable parameters in the given `input.ini`, however you **must not** delete any section and variable from the input file.

For running new calculations, we recommend copying the provided `input.ini` file and `input_data` directory to a separate directory, ideally, where the calculation data is intended to be stored. Then you may edit the calculation parameters in the `input.ini` file (even the file's name) and also provide different input data in the `input_data` directory, such as potential energy surfaces and transition dipole moments.

## Run the executable (after build):

Once the project executable is built using Nix, one can run the calculations using the following command. 

```
$ ./results/bin/ML-TDSE input_file          # Run with input file
```

## Input File

The test input file (`input.ini`) is thoroughly commented with each parameter and available options. However, you may refer to the follwing for a brief description of each parameter.

This section documents all namelists and parameters present in the example `input.ini`. For each parameter the type, units (if applicable), valid values and a short description are given. Do not remove namelist sections or required variables — change values only.

### &grid
- NR (integer)  
  - Example: NR=1024  
  - Description: Number of grid points in the coordinate (R) grid. Must be a power-of-two for best FFT performance but not strictly required.
- Rmin / Rmax (real)  
  - Example: Rmin=0.5, Rmax=51.2  
  - Units: 10^-10 m (Angstrom) in the input file. Converted to atomic units internally.  
  - Description: Left and right boundaries of the R grid (grid extent). Rmin < Rmax required.

Notes: The code constructs an evenly spaced grid of NR points between Rmin and Rmax.

### &nucl_masses
- m1, m2 (real)  
  - Example: m1=1.0, m2=1.0  
  - Units: atomic mass units (amu) as given in the input. Converted to internal mass units in initialization.  
  - Description: Masses of the two nuclei/particles used to form reduced and effective masses used in kinetic energy expressions.

### &time_grid
- dt (real)  
  - Example: dt=0.005  
  - Units: femtoseconds (fs) in the input. Converted to atomic units (a.u.) internally.  
  - Description: Real-time propagation time step.
- Nt (integer)  
  - Example: Nt=50000  
  - Description: Number of time steps to propagate (total simulated time = dt * Nt).

### &elec_states
- Nstates (integer)  
  - Example: Nstates=2  
  - Description: Number of electronic (Born–Oppenheimer) states / channels included in the simulation.
- Elec_pot_kind (string)  
  - Example: Elec_pot_kind="on_grid" or "Morse"  
  - Valid: "on_grid" | "Morse"  
  - Description: Selects how electronic potentials are provided. "on_grid" reads files specified in &input_files, "Morse" generates a simple Morse potential for the ground state (useful for tests).

### &vib_states
- guess_vstates (integer)  
  - Example: guess_vstates=100  
  - Description: Number of vibrational eigenstates to compute in imaginary time propagation (ITP). Determines how many bound states are found.

### &ini_guess_wf
- RI (real)  
  - Example: RI=0.7  
  - Units: Angstrom (10^-10 m) in input.  
  - Description: Center of the Gaussian used as an initial guess for ITP or for Gaussian TDSE initial distributions.
- kappa (real)  
  - Example: kappa=-5.0  
  - Description: Width parameter of the initial Gaussian guess. Often negative in ITP conventions used here.

### &laser_param
Parameters describing one or more laser fields. Fields are encoded using suffixes (e.g. lambda1, tp1, E01 for field #1).
Common keys per field:
- envelope_shape_laserX (string)  
  - Example: "cos2" | "gaussian" | "trapazoidal"  
  - Description: Envelope shape for pulse X.
- lambda1 (real)  
  - Units: nm  
  - Description: Central wavelength of pulse X.
- tp1 (real)  
  - Units: fs  
  - Description: Pulse duration parameter (interpretation depends on envelope shape).
- t_mid1 (real)  
  - Units: fs  
  - Description: Time center of the pulse within the simulation time window.
- E01 (real)  
  - Units: atomic units (a.u.) of field amplitude  
  - Description: Peak field amplitude for pulse X.
- phi1 (real)  
  - Units: multiples of pi (phase in units of π)  
- rise_time1 (real)  
  - Units: fs  
  - Description: Rise time for "trapazoidal" pulses.

Notes: The code supports at most two fields (suffix 1 and 2) in the example input. Add field #2 following the same naming convention.

### &input_files
- input_data_dir (string)  
  - Example: input_data_dir = "input_data/"  
  - Description: Directory containing external data files (potential surfaces, dipoles).
- adb_pot (string)  
  - Example: adb_pot="H2+_BO.dat"  
  - Description: Filename (inside input_data_dir) for adiabatic Born–Oppenheimer potential surfaces used with Elec_pot_kind="on_grid".
- trans_dip_prefix (string)  
  - Example: trans_dip_prefix=""  
  - Description: Optional prefix for transition dipole filenames. If empty, dipole files are expected to be named using state indices (e.g. "12.dat").

### &output_files
- output_data_dir (string)  
  - Example: output_data_dir = "output_data/"  
  - Description: Directory where simulation output files are written. Create the directory or the code will attempt to create it.

### &trans_dip_off
- total_trans_off (integer)  
  - Example: total_trans_off=0  
  - Description: Number of transition dipole entries to parse from `trans_off` string.
- trans_off (string)  
  - Example: trans_off="" or "12 23"  
  - Description: Space-separated list of transitions to switch off (pairs of state indices concatenated, e.g. "12" disables transition between states 1 and 2).

### &absorber_choice
- absorber (string)  
  - Example: absorber="mask"  
  - Valid: "mask" | "CAP"  
  - Description: Selects absorber type used in propagation (mask function or complex absorbing potential). Additional absorber parameters (start position, strength, width) are defined in code and can be exposed in input if required by read routines.

### &ini_state
Initial state for the real-time TDSE propagation.
- initial_distribution (string)  
  - Example: "single vibrational state" | "gaussian distribution" | "Boltzmann distribution"  
  - Description: Selects how to prepare the initial nuclear wavefunction (a single computed vibrational eigenstate, a Gaussian wavepacket or a thermal Boltzmann ensemble).
- N_ini (integer)  
  - Example: N_ini = 1  
  - Description: Electronic-state index in which the initial nuclear state is prepared (1-based).
- v_ini (integer)  
  - Example: v_ini = 2  
  - Description: Vibrational quantum number to use when `initial_distribution = "single vibrational state"`.
- RI_tdse (real)  
  - Example: RI_tdse = 1.5  
  - Units: Angstrom in input. Center of Gaussian initial TDSE wavepacket.
- kappa_tdse (real)  
  - Example: kappa_tdse = -5.0  
  - Description: Width parameter for Gaussian TDSE initial wavepackets.
- temperature / T_boltz (real) *(optional, if supported)*  
  - Example: T_boltz = 300.0  
  - Units: Kelvin (K)  
  - Description: Temperature used for Boltzmann vibrational populations when `initial_distribution = "Boltzmann distribution"`.

Notes: The code supports supplying an explicit vibrational-population vector as an alternative. Check `readinputmodule.f90` for supported keys in your build.

### & parallelization
- prop_par_FFTW (string)  
  - Example: prop_par_FFTW = "parallel" | ""  
  - Description: Controls FFTW threading/parallelization for the real-time propagation FFTs.
- ITP_par_FFTW (string)  
  - Example: ITP_par_FFTW = "" | "parallel"  
  - Description: Controls FFTW threading for imaginary-time propagation (vibrational eigenstates).

Notes: Set to "parallel" to enable FFTW multithreading where supported. Also check environment and FFTW build options.

### &openmp_threads
- omp_nthreads (integer)  
  - Example: omp_nthreads = 4  
  - Description: Number of OpenMP threads to request (0 = use all available threads). Affects threaded BLAS/FFTW/OpenMP regions depending on build.



**Notes:**
- The examples below only show the relevant namelists to change for each use-case.
- For full parameter descriptions consult `input.ini` and `src/readinputmodule.f90` / `src/variablesmodule.f90`.

1) Single vibrational state (start TDSE in one bound vibrational state)
```fortran
&ini_state
  initial_distribution = "single vibrational state"
  N_ini = 1        ! initial electronic state index (1-based)
  v_ini = 2        ! vibrational quantum number to start in
  RI_tdse = 2.0    ! (if using a Gaussian TDSE init; ignored for single vibrational state)
  kappa_tdse = -5.0
/
```

2) Boltzmann distribution (thermal ensemble of vibrational states)
If your build supports a temperature key, set temperature here; otherwise provide a population file as supported by the code. Example (adjust variable name if your read module uses a different key):

```fortran
&ini_state  initial_distribution = "Boltzmann distribution"  
! T_boltz (K) or other temperature key may be required by your input reader  
T_boltz = 300.0  
N_ini = 1           ! electronic state in which the Boltzmann vib populations are prepared  
! alternatively: point to a file with vibrational population (check readinputmodule)
/
```

3) Use CAP absorber (complex absorbing potential) instead of mask absorber
```fortran
&absorber_choice
  absorber = "CAP"    ! switch absorber to CAP (options: "mask" | "CAP")
/
```

Currently the absorber functions are calculated based on the (provided) grid. There are customization chances unless you change the soruce code values in `src/propagation.f90`
