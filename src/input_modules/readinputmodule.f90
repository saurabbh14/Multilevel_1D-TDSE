module ReadInputFile
    use global_vars
    use pulse_mod
    implicit none
    type :: InputFilePath
      character(2000) :: path
    contains
      procedure :: read => read_input_file
      procedure :: read_pulse_params => read_pulse_params
    end type InputFilePath
  contains
    subroutine read_input_file(this)
      class(InputFilePath), intent(inout) :: this
      integer :: input_tk

      namelist /grid/NR
      namelist /nucl_masses/m1,m2
      namelist /time_grid/dt,Nt
      namelist /elec_states/Nstates, Elec_pot_kind
      namelist /vib_states/guess_vstates
      namelist /ini_guess_wf/Ri, kappa
      namelist /input_files/input_data_dir,adb_pot, trans_dip_prefix
      namelist /output_files/output_data_dir
      namelist /trans_dip_off/total_trans_off, trans_off
      namelist /absorber_choice/absorber
      namelist /ini_state/v_ini,N_ini,initial_distribution,temperature,kappa_tdse, RI_tdse
      namelist /parallelization/prop_par_FFTW,ITP_par_FFTW
     
      open(newunit=input_tk, file=adjustl(trim(this%path)), status='old')
      read(input_tk, nml=grid)
      read(input_tk, nml=nucl_masses)
      read(input_tk,nml=time_grid)
      read(input_tk,nml=elec_states)
      read(input_tk,nml=vib_states)
      read(input_tk,nml=ini_guess_wf)
      read(input_tk,nml=input_files)
      read(input_tk,nml=output_files) 
      read(input_tk,nml=trans_dip_off)
      read(input_tk,nml=absorber_choice)
      read(input_tk,nml=ini_state)
      read(input_tk,nml=parallelization)
      close(input_tk)

      ! Assign values to the pulse_param components

    end subroutine read_input_file

    subroutine read_pulse_params(this, pulse)
      class(InputFilePath), intent(inout) :: this
      type(pulse_param), intent(inout) :: pulse
      integer :: input_tk

      ! Intermediate variables for pulse_param components
      character(150) :: envelope_shape_laser1, envelope_shape_laser2
      real(dp) :: lambda1, lambda2, tp1, tp2, t_mid1, t_mid2
      real(dp) :: E01, E02, phi1, phi2, rise_time1, rise_time2

      namelist /laser_param/envelope_shape_laser1, envelope_shape_laser2, &
      & lambda1, lambda2, tp1, tp2, t_mid1, t_mid2, E01, E02, & 
      & phi1, phi2, rise_time1, rise_time2

      open(newunit=input_tk, file=adjustl(trim(this%path)), status='old')
      read(input_tk,nml=laser_param)
      close(input_tk)

      ! Assign values to the pulse_param components
      pulse%envelope_shape_laser1 = envelope_shape_laser1
      pulse%envelope_shape_laser2 = envelope_shape_laser2
      print*, "test"
      print*, "Envelope shape laser 1:", trim(pulse%envelope_shape_laser1)
      print*, "Envelope shape laser 2:", trim(pulse%envelope_shape_laser2)
      print*, "test"
      pulse%lambda1 = lambda1
      pulse%lambda2 = lambda2
      pulse%tp1 = tp1
      pulse%tp2 = tp2
      pulse%t_mid1 = t_mid1
      pulse%t_mid2 = t_mid2
      pulse%E01 = E01
      pulse%E02 = E02
      pulse%phi1 = phi1
      pulse%phi2 = phi2
      pulse%rise_time1 = rise_time1
      pulse%rise_time2 = rise_time2

    end subroutine read_pulse_params
  
  end module ReadInputFile