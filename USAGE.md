# Quick reference — USAGE

This file gives short example `input.ini` snippets for common use-cases. Merge these namelists into your main `input.ini` (keep other required sections such as &grid, &elec_states, &time_grid, &laser_param, etc.).

## Run the executable (after build):
### Test Run:
After successfully building the package, a test calculation can be run using the following command.
```
$ ./result/bin/MLTDSE input.ini
```
The test run should create the `output_data` directory with the following directory tree.
```
$ tree output_data
output_data/
├── 12_read.out
├── H2+_BO.dat_read.out
├── nuclear_wavepacket_data
│   ├── BO_Electronic-state-g0_chi0-Evib.out
│   ├── BO_Electronic-state-g0_Evib.out
│   ├── BO_Electronic-state-g0_vibstates.out
│   ├── BO_Electronic-state-g1_chi0-Evib.out
│   ├── BO_Electronic-state-g1_Evib.out
│   ├── BO_Electronic-state-g1_vibstates.out
│   └── Bound-vibstates_in_Nthstates.out
├── pulse_data
│   ├── electric_field1_E0.10_width110.out
│   ├── electric_field2_E0.0000_width0.out
│   ├── envelope1.out
│   ├── envelope2.out
│   ├── Total_electric_field_phi0.00pi.out
│   └── Total_vector_field_phi0.00pi.out
└── time_prop
    ├── absorber_function.out
    ├── avgR_1d.out
    ├── density_1d_pm3d.out
    ├── ex_density_1d_pm3d.out
    ├── field_1d.out
    ├── KER_spectra_from_state_g0.out
    ├── KER_spectra_from_state_g0_unnormalized.out
    ├── KER_spectra_from_state_g1.out
    ├── KER_spectra_from_state_g1_unnormalized.out
    ├── momentum_1d.out
    ├── momt_spectra_from_state_g0.out
    ├── momt_spectra_from_state_g0_unnormalized.out
    ├── momt_spectra_from_state_g1.out
    ├── momt_spectra_from_state_g1_unnormalized.out
    ├── norm_1d.out
    ├── norm_pn_1d.out
    ├── psi0_1d.out
    ├── psi_outR_momt_density_1d_pm3d.out
    ├── psi_outR_norm_1d.out
    ├── Total_KER_spectra_normalized.out
    ├── Total_KER_spectra.out
    ├── Total_momt_spectra_normalized.out
    ├── Total_momt_spectra.out
    └── vibpop1D_lambda.out

```

### Custom runs:
For running new custom calculations, we recommend copying the provided `input.ini` file and `input_data` directory to a separate directory, ideally, where the calculation data is intended to be stored. Then you may edit the calculation parameters in the `input.ini` file (even the file's name) and also provide different input data in the `input_data` directory, such as potential energy surfaces and transition dipole moments. In the current version, it is necessary to include all input parameter sections in the input file, as already given in the test `input.ini` file, and to provide a path to the `input_data` directory in it. 

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
