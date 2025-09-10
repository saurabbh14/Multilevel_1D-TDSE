program TDSE_main

  use CommandLineModule
  use ReadInputFile
  use PrintInputVars
  use global_vars
  use pulse_mod
  use data_au
  implicit none
  type(CommandLine) :: cmd_line
  type(InputFilePath) :: input_path
  type(pulse_param) :: pulse
  integer :: I, J, I_Emax
  integer :: scount,  & ! Starting "time"
             ecount ! Ending "time"
  integer  rate       ! number of clock ticks per second

  real(dp) st, ft,timer_elapsed_time
  real(dp) Emax
  
  call cpu_time(st)
  call system_clock(scount,rate)
  call cmd_line%read()
  print*, "reading input:"
  print*, "General Inputs from ", trim(cmd_line%input)
  input_path%path = trim(cmd_line%input)
  call input_path%read()
  call pulse%read(input_path%path)
  print*, "Done reading input"
  print*, "_________________________"
  call initializer
  print*, "Printing input variables"
  call pulse%initialize()
  call print_input_vars()
  call pulse%param_print()
  call p_grid
  allocate(Vstates(Nstates), chi0(NR,guess_vstates,Nstates))
  chi0 = 0.d0
 
!  call blas_check

!print*,"test"    
!  ewf = 0.d0
!  adb = 0.d0 
!
!  call adiabatic_surface(adb, mu_all) 
 
  call nuclear_wavefkt

  call pulse%generate()
  call pulse%write_to_file()
  call propagation_1D(pulse%El, pulse%Al)
  deallocate (adb, mu_all, chi0)
  call pulse%deallocate_all()

  call cpu_time(ft)
  print*,'Run time=', ft-st, 'seconds' 
  call system_clock(ecount)
  timer_elapsed_time = real(ecount-scount,8)/real(rate,8)
  write(*,*) "Calculated run time is ",timer_elapsed_time," seconds"
  
end program


! _______________ Subroutines __________________________________________________



subroutine initializer

use global_vars
use pot_param
implicit none
  real:: dummy, dummy2, dummy3, dummy4
  character(200):: mk_out_dir, input_path
  integer:: i, j, input_tk
 
  write(mk_out_dir, '(a)') adjustl(trim(output_data_dir))
  print*, "creating output directory ", trim(mk_out_dir)
  call execute_command_line("mkdir -p " // adjustl(trim(mk_out_dir)))

  allocate(R(NR), en(NR))
  allocate(pR(NR))
  allocate(adb(NR,Nstates),mu_all(Nstates,Nstates,NR))

  call pot_read

  call trans_dipole_read

!  R0=0.1/au2a
!  Rend=15.d0/au2a
!  R =R/au2a
  R0 =R(1)
  Rend = R(NR)
  dR =R(2)-R(1)!(Rend - R0) / (NR - 1)
  dpr = (2._dp * pi) / (dR * NR) 
  

 !Grids: co-ordinate space
!  do I =1,NR
!    R(I) = R0 + (I - 1) * dR
!  enddo
!--------------------------
!Masses:
  m1=m1*mass
  m2=m2*mass 
  mn=m1+m2
  m_red = m1*m2/(m1+m2)
!  m_eff= (m1*me+m2*me+m1*m2)/(m1*m2*me)
  m_eff = (m1+m2) / (m1+m2+1.0_dp)!(4.d0 * me * mass) / (2.d0 * mass + me) 
  mn1=m1/mn
  mn2=m2/mn
!----------------------------  
  kap=(mn+2._dp)/(mn+1._dp)
  lam=(m2-m1)/mn
  RI= RI / au2a   
 
  dt = dt / au2fs 
                
end subroutine
 
!...................... Impulsgrid......................
 
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
!%%%%%% File read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------------------
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
      write(filepath,'(a,a)') adjustl(trim(input_data_dir)), adjustl(trim(adb_pot))  
      print*, "Potential surfaces in path:", trim(filepath)
      open(newunit=pot_tk,file=adjustl(trim(filepath)),status='unknown')
      do I = 1, NR
        read(pot_tk,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
            ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
      end do
      close(pot_tk)
      write(filepath,'(a,a,a)') adjustl(trim(output_data_dir)), adjustl(trim(adb_pot)),&
            & "_read.out"  
      open(newunit=pot_out_tk,file=adjustl(trim(filepath)),status='unknown')
      do I = 1, NR
        write(pot_out_tk,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
            ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
      end do
      close(pot_out_tk)
 
  ! Morse potential as the ground state
    case("Morse")
      do I =1, NR
        adb(I,1) = morse_potential(0.17_dp,1.85_dp,0.743_dp/au2a,R(I))
      enddo
   
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
  
  trans_off = trim(adjustl(trans_off))
  read(trans_off,*) trans_off_parse

  print*, "Transitions to be switched off"
  do L =1, total_trans_off
    trans_off_parse(L) = trim(adjustl(trans_off_parse(L)))
    print*, "Transition", L,": ", trans_off_parse(L) 
  enddo
   
  print*, "Transition dipoles with file prefix \", trim(trans_dip_prefix), " \."
  mu_all = 0._dp
! transition dipole moments of all states
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
    do L = 1, total_trans_off
      tr = trans_off_parse(L)
      read(tr(1:1),*) N1
      read(tr(2:2),*) N2 
      mu_all(N1,N2,:) = 0._dp 
    enddo
    mu_all = abs(mu_all)
  
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
use global_vars, only: dp
implicit none
  real(dp), intent(in):: de, a, re, r
  real(dp) :: pot

  pot = de * (1._dp - exp(-a * (r-re)))**2

end function
