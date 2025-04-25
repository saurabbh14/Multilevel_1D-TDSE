module VarPrecision
! This module defines the precision and data types used in the program.
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
    use CommandLineModule
    use VarPrecision, only: dp, idp
    use, intrinsic :: iso_c_binding
    ! R-grid
    integer(C_INT):: NR 
    
    ! electronic states
    integer:: Nstates
    character(200):: Elec_pot_kind
    
    ! vibrational states
    integer:: guess_vstates
    integer, allocatable:: Vstates(:)
     
    ! time grid 
    integer:: Nt
    
    ! masses
    real(dp):: m1, m2
    
    ! guess initial wf
    real(dp):: RI, kappa
    
    ! initial TDSE state
    integer:: N_ini, v_ini
    integer, allocatable:: v_dist_ini(:)
    real(dp):: temperature, kappa_tdse, RI_tdse ! for Boltzmann distribution
    character(2000):: initial_distribution  
    
    ! laser parameters
    character(150):: envelope_shape_laser1, envelope_shape_laser2
    real(dp):: tp1, fwhm, t_mid1, rise_time1
    real(dp):: tp2, t_mid2, rise_time2
    real(dp):: e01, e02, phi1, phi2
    real(dp):: lambda1, lambda2
    
    ! input files
    character(2000):: input_data_dir
    character(2000):: adb_pot, trans_dip_prefix
    character(2000):: output_data_dir
    
    ! transitions switched off
    integer:: total_trans_off
    character(2000):: trans_off
    
    ! Absorber choice
    character(5):: absorber
    
    ! FFTW parallelization
    character(10):: prop_par_FFTW
    character(10):: ITP_par_FFTW
    
    end module InputVars
    
module global_vars
    use InputVars
    real(dp):: dR
    real(dp), allocatable:: R(:)
    real(dp), allocatable:: en(:)
    real(dp), allocatable:: PR(:)
    real(dp), allocatable:: Pot(:,:), chi0(:,:,:)
    real(dp), allocatable, dimension(:,:,:):: mu_all
    real(dp), allocatable, dimension(:,:):: adb
    real(dp):: kap, lam
    real(dp):: dt
    real(dp):: dpr 
    real(dp):: omega1, omega2
    real(dp):: mn, mn1, mn2 !all relevent mass veriables
end module global_vars
    
module data_au
    use VarPrecision, only: dp
    real(dp),parameter:: au2a=0.52917706d0  ! length in a.u. to Angstrom
    real(dp),parameter:: cm2au=4.5554927d-6 ! energy in wavenumbers to a.u.
    real(dp),parameter:: au2fs=0.024        ! time in a.u. to femtoseconds
    real(dp),parameter:: j2eV = 6.242D18    ! transforms energy in J to eV
    real(dp),parameter:: au2eV=27.2116d0    ! energy in a.u. to eV
    real(dp),parameter:: eV2nm=1239.84      ! energy from eV to wavelength(nm)
    real(dp),parameter:: kB = 3.167d-6      ! Boltzmann constant hartree/K
    real(dp),parameter:: i2au=2.0997496D-9
    real(dp),parameter:: e02au=5.142206707e11
    real(dp),parameter:: pi=3.141592653589793d0       ! that's just pi
    real(dp),parameter:: mass=1836.15d0    ! reduced mass of deuterium
    real(dp),parameter:: me =1.d0          ! mass of the electron
    real(dp):: m_eff, m_red                ! effecitve mass of el. and nucl.
    real(dp),parameter:: c_speed =137.03604 ! speed of light in a.u.
    complex(dp),parameter:: im=(0.d0,1.d0) 
end module

module pot_param
    use data_au
    real(dp):: R0     ! Grid-Parameter, start..
    real(dp)::Rend   !..and end
    real(dp),parameter:: cpmR=3.2d0*2*2 !*2 !absorber position from the end of R-grid
end module pot_param
    