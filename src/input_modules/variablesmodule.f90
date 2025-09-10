!> This file contains variable and precision definitions for the TDSE program.
!> It includes modules for variable precision, input variables, global variables,
!> atomic units, and potential parameters.
module VarPrecision
! Defines application-wide numeric kinds and standard I/O unit constants.
! Keep these short aliases so the rest of the code can use `dp` / `sp` etc.
    use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64, &
                                input_unit, output_unit, error_unit
    implicit none
    integer, parameter :: sp        = real32
    integer, parameter :: dp        = real64
    integer, parameter :: idp        = int64
    integer, parameter :: isp       = int32
    integer, parameter :: stdin     = input_unit
    integer, parameter :: stdout    = output_unit
    integer, parameter :: stderr    = error_unit
end module VarPrecision
    
module InputVars
! High-level input variables read from input file.
! These are the *declarative* variables describing the simulation setup.
! Keep names and comments in sync with readinputmodule.f90.
    use CommandLineModule
    use VarPrecision, only: dp, idp
    use, intrinsic :: iso_c_binding

    ! R-grid
    integer(C_INT):: NR                    ! number of grid points (coordinate space)
    
    ! electronic states
    integer:: Nstates                      ! number of electronic BO states
    character(200):: Elec_pot_kind         ! "on_grid" | "Morse" (select potential source)
    
    ! vibrational states 
    integer:: guess_vstates                ! number of vibrational eigenstates to compute
    integer, allocatable:: Vstates(:)      ! storage for computed vibrational energies
    
    ! time grid 
    integer:: Nt                           ! number of time steps for time propagation
    
    ! masses (input in atomic mass units or as specified; converted later)
    real(dp):: m1, m2                      ! masses of particle 1 and 2 (in code units before conversion)
    
    ! guess initial wavefunction
    real(dp):: RI, kappa                   ! RI: center of initial Gaussian (units: Angstrom unless converted)
    
    ! initial TDSE state (how to prepare the initial wavefunction for real-time propagation)
    integer:: N_ini, v_ini                 ! N_ini: electronic state index; v_ini: vibrational quantum number
    integer, allocatable:: v_dist_ini(:)   ! optional explicit vibrational-population vector
    real(dp):: temperature, kappa_tdse, RI_tdse ! parameters used for Boltzmann/Gaussian TDSE initial distributions
    character(2000):: initial_distribution ! string selecting initial distribution type ("single vibrational state", "Boltzmann distribution", etc.)
    
    ! input / output file paths and prefixes
    character(2000):: input_data_dir       ! directory with input grids, dipoles, potentials
    character(2000):: adb_pot, trans_dip_prefix ! adb_pot: filename for BO surfaces; trans_dip_prefix: optional prefix for dipole files
    character(2000):: output_data_dir      ! directory to write outputs
    
    ! transitions to be switched off (e.g. "12 23")
    integer:: total_trans_off
    character(2000):: trans_off
    
    ! Absorber choice for propagation (mask or CAP)
    character(5):: absorber                ! "mask" | "CAP"
    
    ! FFTW parallelization flags read from input
    character(10):: prop_par_FFTW
    character(10):: ITP_par_FFTW
    
end module InputVars
    
module global_vars
! Global arrays and derived variables used across modules.
    use InputVars
    use InputVars

    real(dp):: dR                           ! grid spacing in coordinate space (R)
    real(dp), allocatable:: R(:)            ! coordinate grid
    real(dp), allocatable:: en(:)           ! energy array
    real(dp), allocatable:: PR(:)           ! momentum grid
    real(dp), allocatable:: Pot(:,:), chi0(:,:,:) ! potentials & vibrational wavefunctions
    real(dp), allocatable, dimension(:,:,:):: mu_all ! transition dipole arrays
    real(dp), allocatable, dimension(:,:):: adb ! adiabatic BO potentials
    real(dp):: kap, lam                      ! derived coefficients used in dipole / kinetic expressions
    real(dp):: dt                            ! time step 
    real(dp):: dpr                           ! momentum-grid spacing
    real(dp):: mn, mn1, mn2                  ! total mass and individual mass ratios
end module global_vars
    
module data_au
! Physical constants and unit-conversion factors used throughout the code.
! All constants are in double precision (dp).
    use VarPrecision, only: dp
    real(dp),parameter:: au2a=0.52917706d0  ! length: atomic units -> Angstrom
    real(dp),parameter:: cm2au=4.5554927d-6 ! energy: wavenumbers (cm^-1) -> a.u.
    real(dp),parameter:: au2fs=0.024        ! time: atomic units -> femtoseconds
    real(dp),parameter:: j2eV = 6.242D18    ! energy: Joule -> eV (approx.)
    real(dp),parameter:: au2eV=27.2116d0    ! energy: a.u. -> eV
    real(dp),parameter:: eV2nm=1239.84      ! conversion: eV -> nm (lambda)
    real(dp),parameter:: kB = 3.167d-6      ! Boltzmann constant in hartree/K
    real(dp),parameter:: i2au=2.0997496D-9  ! intensity conversion (W/cm^2 -> a.u.)
    real(dp),parameter:: e02au=5.142206707e11 ! electric field conversion (V/m -> a.u.)
    real(dp),parameter:: pi=3.141592653589793d0
    real(dp),parameter:: mass=1836.15d0     ! electron mass unit reference (e.g. proton mass / electron mass)
    real(dp),parameter:: me =1.d0           ! electron mass 1 a.u.
    real(dp):: m_eff, m_red                  ! placeholders for effective and reduced mass
    real(dp),parameter:: c_speed =137.03604  ! speed of light in a.u. (approx.)
    complex(dp),parameter:: im=(0.d0,1.d0)   ! imaginary unit
end module

module pot_param
! Parameters related to potentials and absorbers.
    use data_au
    real(dp):: R0     ! Grid-Parameter, start..
    real(dp)::Rend   !..and end
    real(dp),parameter:: cpmR=3.2d0*2*2 !*2 !absorber position from the end of R-grid
end module pot_param
    