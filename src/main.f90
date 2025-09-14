!> Main driver for the Multi-level 1D TDSE solver.
!> This file contains:
!>  - program TDSE_main : entry point that reads input, initializes data,
!>                       computes bound vibrational states, generates pulses
!>                       and runs the real-time propagation.

program TDSE_main

  use CommandLineModule         ! parse command line options
  use ReadInputFile             ! read input.ini (namelists / parameters)
  use PrintInputVars            ! printing helper for chosen simulation variables
  use global_vars               ! shared allocatables and simulation variables
  use pulse_mod                 ! pulse generation and IO
  use data_au                   ! atomic units / unit conversion constants
  use initializer               ! subroutine to initialize grids, arrays
  use nuclear_wavefkt          ! compute vibrational eigenstates via ITP
  use propagation_mod
  implicit none

  ! Local types/objects
  type(CommandLine) :: cmd_line      ! command line handler object
  type(InputFilePath) :: input_path  ! input file path holder
  type(pulse_param) :: pulse         ! pulse object with methods (read, init, gen)
  type(initializer_type) :: init_obj ! initializer object
  type(nuclear_wavefkt_class) :: nwf_obj ! nuclear wavefunction object
  type(time_prop) :: prop
  
  ! Local counters and timers
  integer :: I, J, I_Emax
  integer :: scount, ecount        ! system clock ticks (start/end)
  integer :: rate                  ! clock ticks per second
  real(dp) :: st, ft, timer_elapsed_time  ! cpu_time start/end and elapsed
  real(dp) :: Emax                 ! placeholder for energy ranges
  
  ! Start overall timing
  call cpu_time(st)
  call system_clock(scount,rate)

  ! Read command line and input file
  call cmd_line%read()
  print*, "reading input:"
  print*, "General Inputs from ", trim(cmd_line%input)
  input_path%path = trim(cmd_line%input)
  call input_path%read()

  ! Read pulse parameters from input file
  call pulse%read(input_path%path)
  print*, "Done reading input"
  print*, "_________________________"

  ! Prepare working directories, grids and arrays
  call init_obj%setup

  ! Print and initialize pulse parameters
  print*, "Printing input variables"
  call pulse%initialize()
  call print_input_vars()
  call pulse%param_print()

  ! Allocate arrays for computed vibrational states
  ! Vstates: eigenvalues/energies, chi0: wavefunctions (R, state index, vib index)
  allocate(Vstates(Nstates), chi0(NR,guess_vstates,Nstates))
  chi0 = 0._dp
 
  ! Compute vibrational eigenstates via Imaginary Time Propagation (ITP)
  call nwf_obj%nuclear_wf_calc

  ! Generate and store pulse(s) then run time propagation
  call pulse%generate()
  call pulse%write_to_file()

  print*, "Time dependent calcuations"

  call prop%propagation_1D(pulse%El, pulse%Al)

  ! Clean up allocated arrays and pulse internals
  deallocate (adb, mu_all, chi0)
  call pulse%deallocate_all()

  ! Final timing and report
  call cpu_time(ft)
  print*,'Run time=', ft-st, 'seconds' 
  call system_clock(ecount)
  timer_elapsed_time = real(ecount-scount,8)/real(rate,8)
  write(*,*) "Calculated run time is ",timer_elapsed_time," seconds"
  
end program


