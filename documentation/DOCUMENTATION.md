# Multi-level 1D TDSE solver — Concise Documentation

This document summarizes the repository layout, build & run steps, the main modules/routines, and the important input/output artifacts.

## Quick start

- Build (Nix as used by the project):
  - From project root run:
    ```
    nix build
    ```
  - The executable is produced at `result/bin/ML-TDSE` (see [meson.build](../meson.build) and [src/meson.build](src/meson.build)).

- Run: `./result/bin/ML-TDSE input.ini`

- Example input template: [input.ini](../input.ini)
- Top-level README: [README.md](../README.md)

## High-level workflow

1. Parse command line and read `input.ini`:
 - Command-line handling: [`CommandLineModule.CommandLine`](../src/input_modules/commandlinemodule.f90) ([file](../src/input_modules/commandlinemodule.f90))
 - Input file reading: [`ReadInputFile.InputFilePath`](../src/input_modules/readinputmodule.f90) ([file](../src/input_modules/readinputmodule.f90))

2. Initialize grids, potentials and dipoles:
 - Main program entry: [`main.TDSE_main`](../src/main.f90) ([file](../src/main.f90))
 - Initializers and grid helpers: [`main.initializer`](../src/main.f90), [`main.p_grid`](../src/main.f90)
 - Potential read / prefix handling: [`main.pot_read`](../src/main.f90), [`main.trans_dipole_read`](../src/main.f90)
 - Optional Morse potential helper: [`main.morse_potential`](../src/main.f90)

3. Compute bound vibrational states (imaginary-time propagation):
 - Vibrational eigenstates generator: [`nuclear_wv.nuclear_wavefkt`](../src/nuclear_wv.f90) ([file](../src/nuclear_wv.f90))

4. Generate laser pulses:
 - Pulse type and generation object: [`pulse_mod.pulse_param`](../src/input_modules/pulse.f90) ([file](../src/input_modules/pulse.f90))

5. Real-time propagation & observables:
 - Propagation driver: [`propagation.propagation_1D`](../src/propagation.f90) ([file](../src/propagation.f90))
 - FFTW integration: [`FFTW3` module](../src/input_modules/fftw3.f90) ([file](../src/input_modules/fftw3.f90))

6. Outputs: many `.out` files are written into the configured `output_data_dir` (see `input.ini`).

## Important files & symbols (openable)

- Project root:
- [README.md](../README.md)
- [input.ini](../input.ini)
- [meson.build](../meson.build)

- Build config:
- [src/meson.build](../src/meson.build)

- Main program and helpers:
- [`main.TDSE_main`](../src/main.f90) — program entry ([file](../src/main.f90))
- [`main.initializer`](../src/main.f90) — creates output dirs, allocates grids ([file](../src/main.f90))
- [`main.p_grid`](../src/main.f90) — momentum grid ([file](../src/main.f90))
- [`main.pot_read`](../src/main.f90) — read adiabatic potentials ([file](../src/main.f90))
- [`main.trans_dipole_read`](../src/main.f90) — read transition dipoles ([file](../src/main.f90))
- [`main.morse_potential`](../src/main.f90) — helper potential ([file](../src/main.f90))

- Nuclear / vibrational routines:
- [`nuclear_wv.nuclear_wavefkt`](../src/nuclear_wv.f90) — computes vibrational states (ITP) ([file](../src/nuclear_wv.f90))
- Helpers in that file: `eigenvalue_R`, `integ_r`, FFTW plan setup in ITP

- Propagation & observables:
- [`propagation.propagation_1D`](../src/propagation.f90) — main real-time TDSE loop, absorbers, KER/momentum spectra ([file](../src/propagation.f90))
- Absorber and mask functions: `complex_absorber_function`, `mask_function_cos`, `mask_function_ex`
- BLAS/LAPACK interfaces and calls are in [src/propagation.f90](../src/propagation.f90) (see `blas_interfaces_module`)

- Pulse / field handling:
- [`pulse_mod.pulse_param`](../src/input_modules/pulse.f90) — read, initialize, generate and write pulses (envelopes, spectra) ([file](../src/input_modules/pulse.f90))

- Input & variables modules:
- [`VarPrecision`](../src/input_modules/variablesmodule.f90) — precision and type aliases ([file](../src/input_modules/variablesmodule.f90))
- [`InputVars` / `global_vars`](../src/input_modules/variablesmodule.f90) — main shared simulation variables and allocatables ([file](../src/input_modules/variablesmodule.f90))

- FFTW wrapper:
- [`FFTW3` module](../src/input_modules/fftw3.f90) — includes FFTW Fortran bindings ([file](../src/input_modules/fftw3.f90))

## Input (input.ini)
- The code uses Fortran namelists read by `ReadInputFile`:
- Grid: NR, R0/Rend derived from adb or config
- Time grid: dt, Nt
- Electronic states: Nstates, Elec_pot_kind
- Vibrational options: guess_vstates, initial distribution params
- Laser params: read into `pulse_mod.pulse_param` via `pulse_mod%read`
- See [input.ini](../input.ini) for a runnable example.

## Typical outputs (examples)
- Vibrational states and energies:
- BO electronic-state files: `BO_Electronic-state-g##_vibstates.out` (written by `nuclear_wv`)
- Bound-state wavefunctions: `BO_Electronic-state-g##_chi0-Evib.out`
- Propagation outputs (written by `propagation_1D`):
- `density_1d.out`, `Pdensity_1d.out`, `avgR_1d.out`, `avgPR_1d.out`, `norm_1d.out`
- KER and momentum spectra: `KER_spectra_from_state_g##_unnormalized.out`, `momt_spectra_from_state_g##_unnormalized.out`
- Pulse/pulse spectra written by `pulse_mod` into `output_data_dir/pulse_data/`

## Parallelization / performance notes
- FFTW threading is used optionally:
- ITP and TDSE propagation threads are selectable via namelist fields `ITP_par_FFTW` and `prop_par_FFTW` (see [src/input_modules/readinputmodule.f90](../src/input_modules/readinputmodule.f90)).
- FFTW init & plan code appears in both [`nuclear_wv.nuclear_wavefkt`](../src/nuclear_wv.f90) and [`propagation.propagation_1D`](../src/propagation.f90).
- BLAS/LAPACK dependencies are declared in [src/meson.build](../src/meson.build).

## Where to look for changes / extension points
- Change pulse shapes & spectra: [src/input_modules/pulse.f90](../src/input_modules/pulse.f90) — the `pulse_param` type and its methods.
- Add/modify potentials or dipoles: [src/main.f90](../src/main.f90) (`pot_read`, `trans_dipole_read`) and the `adb` / `mu_all` arrays defined in [src/input_modules/variablesmodule.f90](../src/input_modules/variablesmodule.f90).
- Modify propagation/observables: [src/propagation.f90](../src/propagation.f90) — main loop, absorbers, spectra writers.

## Tests / verification
- There is a test input for H2+ referenced in the top-level README. Run the executable with that `input.ini` and inspect the produced `output_data_dir` files.
- Useful debug prints are present across modules (normalization, energies, run-time).