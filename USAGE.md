# Quick reference — USAGE

This file gives short example `input.ini` snippets for common use-cases. Merge these namelists into your main `input.ini` (keep other required sections such as &grid, &elec_states, &time_grid, &laser_param, etc.).

Run the executable (after build):

Notes:
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