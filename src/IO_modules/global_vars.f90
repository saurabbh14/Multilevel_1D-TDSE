module global_vars
! Global arrays and derived variables used across modules.
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
    real(dp):: m_eff, m_red                  ! effective and reduced mass
    real(dp):: mn, mn1, mn2                  ! total mass and individual mass ratios
    real(dp), allocatable:: time(:)         ! time grid
    character(len=2000) :: pulse_data_dir   ! output directory path for pulse data
    character(len=2000) :: nucl_wf_dir      ! output directory path for nuclear wavefunction data
    character(len=2000) :: time_prop_dir    ! output directory path for time propagation data
end module global_vars