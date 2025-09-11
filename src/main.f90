!> Main driver for the Multi-level 1D TDSE solver.
!> This file contains:
!>  - program TDSE_main : entry point that reads input, initializes data,
!>                       computes bound vibrational states, generates pulses
!>                       and runs the real-time propagation.
!>  - helper subroutines:
!>       initializer   : prepare grids, allocate arrays, read potentials/dipoles
!>       p_grid        : construct momentum grid used by FFT-based propagate
!>       pot_read      : read or generate electronic potential surfaces
!>       trans_dipole_read : read transition dipole moments for all pairs
!>       morse_potential   : utility to create a Morse potential (if requested)
!>
!> Notes:
!>  - For detailed variable definitions, see src/input_modules/variablesmodule.f90
!>  - Time and unit conversions are done here (e.g. dt conversion to atomic units).

program TDSE_main

  use CommandLineModule         ! parse command line options
  use ReadInputFile             ! read input.ini (namelists / parameters)
  use PrintInputVars            ! printing helper for chosen simulation variables
  use global_vars               ! shared allocatables and simulation variables
  use pulse_mod                 ! pulse generation and IO
  use data_au                   ! atomic units / unit conversion constants
  implicit none

  ! Local types/objects
  type(CommandLine) :: cmd_line      ! command line handler object
  type(InputFilePath) :: input_path  ! input file path holder
  type(pulse_param) :: pulse         ! pulse object with methods (read, init, gen)
  
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
  call initializer

  ! Print and initialize pulse parameters
  print*, "Printing input variables"
  call pulse%initialize()
  call print_input_vars()
  call pulse%param_print()

  ! Construct momentum grid used by spectral operators (FFT)
  call p_grid

  ! Allocate arrays for computed vibrational states
  ! Vstates: eigenvalues/energies, chi0: wavefunctions (R, state index, vib index)
  allocate(Vstates(Nstates), chi0(NR,guess_vstates,Nstates))
  chi0 = 0.d0
 
  ! Compute vibrational eigenstates via Imaginary Time Propagation (ITP)
  call nuclear_wavefkt

  ! Generate and store pulse(s) then run time propagation
  call pulse%generate()
  call pulse%write_to_file()
  call propagation_1D(pulse%El, pulse%Al)

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


! _______________ Subroutines __________________________________________________

!> Initializes directories, grids, masses and reads potentials/dipoles.
!> Allocations and unit conversions happen here.

subroutine initializer

use global_vars
use pot_param
implicit none
  real:: dummy, dummy2, dummy3, dummy4
  character(200):: mk_out_dir, input_path
  integer:: i, j, input_tk
 
  ! Ensure output directory exists (create if missing)
  write(mk_out_dir, '(a)') adjustl(trim(output_data_dir))
  print*, "creating output directory ", trim(mk_out_dir)
  call execute_command_line("mkdir -p " // adjustl(trim(mk_out_dir)))

  ! Allocate coordinate and potential arrays now that NR is known
  allocate(R(NR), en(NR))
  allocate(pR(NR))
  allocate(adb(NR,Nstates),mu_all(Nstates,Nstates,NR))

  ! Read potentials and transition dipoles (or generate synthetic ones)
  call pot_read
  call trans_dipole_read

  ! Set derived grid parameters and unit conversions
  R0 = R(1)                  ! leftmost grid point (coordinate space)
  Rend = R(NR)               ! rightmost grid point
  dR = R(2) - R(1)           ! grid spacing (assumed uniform)
  dpr = (2._dp * pi) / (dR * NR)  ! momentum-grid spacing via FFT conventions
  

 !Grids: co-ordinate space
!  do I =1,NR
!    R(I) = R0 + (I - 1) * dR
!  enddo
!--------------------------
  ! Masses: convert to internal mass units and compute reduced/effective masses
  m1 = m1 * mass
  m2 = m2 * mass
  mn = m1 + m2
  m_red = m1*m2/(m1+m2)
  ! Effective mass used in some reduced-dimension expressions (keeps consistency)
  m_eff = (m1 + m2) / (m1 + m2 + 1.0_dp)
  mn1 = m1 / mn
  mn2 = m2 / mn
!----------------------------  

  ! diapole parameters
  kap = (mn + 2._dp) / (mn + 1._dp)
  lam = (m2 - m1) / mn

  ! Convert initial Gaussian center from Angstrom (or given units) to atomic units
  RI = RI / au2a
 
  ! Convert time step from femtoseconds to atomic units for propagation
  dt = dt / au2fs
                
end subroutine
 
!...................... Impulsgrid......................
!> Construct the momentum grid consistent with the FFT layout used in split-operator.
!> Uses the standard ordering: 0, dp, 2dp, ..., (NR/2-1)dp, -NR/2 dp, ..., -dp 

subroutine p_grid

use global_vars, only:pR, dpR, NR, R
use pot_param
 implicit none 
 integer:: I 
  
  do I = 1, NR  
    if (I.le.(NR / 2)) then    
      PR(I) = (I - 1) * dpR    
    else    
      PR(I) = - (NR + 1 - I) * dpR    
    end if
  end do
    print*, 'R0=', R(1), 'Rend=', R(NR)
    print*, 'PR0=', PR((NR/2)+1), 'PRend=', PR(NR/2)
  return
end subroutine  



!------------------------------------------------------------------------------
!%%%%%% File IO: reading potentials and dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------------------

!> Read electronic Born–Oppenheimer potential surfaces from input file or
!> construct a Morse potential when Elec_pot_kind="Morse".
subroutine pot_read
use global_vars, only:R, NR, adb, adb_pot, &
        & input_data_dir, output_data_dir, &
        & Elec_pot_kind, dp
use data_au
implicit none

  character(2000):: filepath
  integer:: I, pot_tk, pot_out_tk
  real(dp):: morse_potential
  real(dp):: dummy

  select case(Elec_pot_kind)
    case ("on_grid")
      ! Compose full path and read potential data file with NR lines.
      write(filepath,'(a,a)') adjustl(trim(input_data_dir)), adjustl(trim(adb_pot))  
      print*, "Potential surfaces in path:", trim(filepath)
      open(newunit=pot_tk,file=adjustl(trim(filepath)),status='unknown')
      do I = 1, NR
        read(pot_tk,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
            ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
      end do
      close(pot_tk)

      ! Also write a copy into output directory for verification
      write(filepath,'(a,a,a)') adjustl(trim(output_data_dir)), adjustl(trim(adb_pot)),&
            & "_read.out"  
      open(newunit=pot_out_tk,file=adjustl(trim(filepath)),status='unknown')
      do I = 1, NR
        write(pot_out_tk,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
            ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
      end do
      close(pot_out_tk)
 
    case("Morse")
      ! Fill only the ground state with a Morse potential
      do I =1, NR
        adb(I,1) = morse_potential(0.17_dp,1.85_dp,0.743_dp/au2a,R(I))
      enddo
   
      ! Write generated Morse surface to output 
      write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Morse_pot_read.out"  
      open(newunit=pot_out_tk,file=adjustl(trim(filepath)),status='unknown')
      do I = 1, NR
        write(pot_out_tk,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
            ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
      end do
      close(pot_out_tk)
  end select

end subroutine

!------------------------------------------------------------------------------

!> Read transition dipole moments for each electronic-state pair.
!> Supports optional prefixing and allows switching off specific transitions.
subroutine trans_dipole_read
use global_vars, only:R, NR, mu_all, trans_dip_prefix, Nstates, &
       & input_data_dir, output_data_dir, dp, total_trans_off, &
       & trans_off
implicit none
  integer:: I, L, M, N1, N2
  character(2000):: fn
  character(2):: trans_off_parse(total_trans_off)
  character(2):: tr
  integer:: input_tk, output_tk
  real(dp):: dummy
  
  ! Parse list of transitions to disable (e.g. "12 23")
  trans_off = trim(adjustl(trans_off))
  read(trans_off,*) trans_off_parse

  print*, "Transitions to be switched off"
  do L =1, total_trans_off
    trans_off_parse(L) = trim(adjustl(trans_off_parse(L)))
    print*, "Transition", L,": ", trans_off_parse(L) 
  enddo
   
  print*, "Transition dipoles with file prefix \", trim(trans_dip_prefix), " \."
  mu_all = 0._dp

  ! Two supported modes:
  !  - no prefix: files are named input_data_dir + "<L><M>.dat"
  !  - with prefix: files are input_data_dir + trans_dip_prefix + "<L><M>.dat"
  if (trim(trans_dip_prefix) .eq.'') then
    do L = 1, Nstates
      do M = L+1, Nstates
        write(fn,fmt='(a,i0,i0,a)') adjustl(trim(input_data_dir)),L,M,'.dat'
        print*, trim(fn)
        open(newunit=input_tk, file=adjustl(trim(fn)), form='formatted')
        do I = 1, NR
          read(input_tk,*) dummy, mu_all(L,M,I)
        enddo
        close(input_tk)
      enddo
    enddo
    ! Zero-out any transitions explicitly switched off in input (trans_off)
    do L = 1, total_trans_off
      tr = trans_off_parse(L)
      read(tr(1:1),*) N1
      read(tr(2:2),*) N2 
      mu_all(N1,N2,:) = 0._dp 
    enddo
    mu_all = abs(mu_all) ! ensure positive magnitudes
  
    ! Write read dipoles to output for inspection
    do L = 1, Nstates
      do M = L+1, Nstates
        write(fn,fmt='(a,i0,i0,a)') adjustl(trim(output_data_dir)),L,M,'_read.out'
        print*, trim(fn)
        open(newunit=output_tk, file=adjustl(trim(fn)), form='formatted')
        do I=1,NR
          write(output_tk,*) R(I), mu_all(L,M,I)
        enddo
        close(output_tk)
      enddo
    enddo

  else
    ! Prefix-mode: read files with a common prefix + state indices
    do L = 1, Nstates
      do M = L+1, Nstates
        write(fn,fmt='(a,a,i0,i0,a)') adjustl(trim(input_data_dir)), &
                 & adjustl(trim(trans_dip_prefix)),L,M,'.dat'
        print*, trim(fn)
        open(newunit=input_tk, file=adjustl(trim(fn)), form='formatted')
        do I = 1, NR
          read(input_tk,*) dummy, mu_all(L,M,I)
        enddo
        close(input_tk)
      enddo
    enddo
 
    ! Write read dipoles to output directory with prefix
    do L = 1, Nstates
      do M = L+1, Nstates
        write(fn,fmt='(a,a,i0,i0,a)') adjustl(trim(output_data_dir)), &
                 & adjustl(trim(trans_dip_prefix)),L,M,'_read.out'
        print*, trim(fn)
        open(newunit=output_tk, file=adjustl(trim(fn)), form='formatted')
        do I=1,NR
          write(output_tk,*) R(I), mu_all(L,M,I)
        enddo
        close(output_tk)
      enddo
    enddo
  endif

end subroutine

!------------------------------------------------
function morse_potential(de,a,re,r) result(pot)
! Simple Morse potential generator used when Elec_pot_kind="Morse".
! Parameters:
!  - de : dissociation energy (in same units as returned pot)
!  - a  : range parameter controlling width
!  - re : equilibrium bond length (in same units as r)
!  - r  : coordinate at which to evaluate potential
use global_vars, only: dp
implicit none
  real(dp), intent(in):: de, a, re, r
  real(dp) :: pot

  pot = de * (1._dp - exp(-a * (r-re)))**2

end function
