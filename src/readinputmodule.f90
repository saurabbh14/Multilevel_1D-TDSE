module ReadInputFile
    use global_vars
    implicit none
    type :: InputFilePath
      character(2000) :: path
    contains
      procedure :: read => read_input_file
    end type InputFilePath
  contains
    subroutine read_input_file(this)
      class(InputFilePath), intent(in) :: this
      integer :: input_tk
      namelist /grid/NR
      namelist /nucl_masses/m1,m2
      namelist /time_grid/dt,Nt
      namelist /elec_states/Nstates, Elec_pot_kind
      namelist /vib_states/guess_vstates
      namelist /ini_guess_wf/Ri, kappa
      namelist /laser_param/envelope_shape_laser1, envelope_shape_laser2, &
              & lambda1,lambda2,tp1,tp2,t_mid1,t_mid2,E01,E02,phi1,phi2, &
              & rise_time1, rise_time2
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
      read(input_tk,nml=laser_param)
      read(input_tk,nml=input_files)
      read(input_tk,nml=output_files) 
      read(input_tk,nml=trans_dip_off)
      read(input_tk,nml=absorber_choice)
      read(input_tk,nml=ini_state)
      read(input_tk,nml=parallelization)
      close(input_tk)
      
    end subroutine read_input_file
  
  end module ReadInputFile