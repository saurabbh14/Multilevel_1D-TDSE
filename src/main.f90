module commandline_args
character(2000):: command_line
character(2000):: input

!call read_command_line
!call parse_command_line
contains

 subroutine read_command_line
 integer :: exenamelength
 integer :: io, io2

 command_line = ""
! call get_command(command = command_line,status = io)
 call get_command_argument(1,command_line,status = io)
 if (io==0) then
   input = trim(command_line)
! Uncomment if parsing of commandline is needed, i.e. multiple argments in the commandline         
!  call get_command_argument(0,length = exenamelength,status = io2)
!  if (io2==0) then
!    command_line = "&cmd "//adjustl(trim(command_line(exenamelength+1:)))//" /"
!  else
!    command_line = "&cmd "//adjustl(trim(command_line))//" /"
!  end if
   else
     write(*,*) io,"Error getting command line."
   end if
  end subroutine

! subroutine parse_command_line
!   character(256) :: msg
!   namelist /cmd/ input, adb_pot, trans_dip
!   integer :: io
!
!   if (len_trim(command_line)>0) then
!     msg = ''
!     read(command_line,nml = cmd,iostat = io,iomsg = msg)
!     if (io/=0) then
!       error stop "Error parsing the command line or cmd.conf " // msg
!     end if
!   end if
! end subroutine

end module commandline_args

module var_precision
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
end module var_precision

module input_vars
 use commandline_args
 use var_precision, only: dp, idp
 use, intrinsic :: iso_c_binding
! R-grid
 integer(C_INT):: NR 

! electronic states
 integer:: Nstates

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

end module input_vars

module global_vars
 use input_vars
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
end module

module data_au
 use var_precision, only: dp
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

module FFTW3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
!  include '/usr/include/fftw3.f03'                                        ! Desktop packet
!  include '/home/me23jok/ProjectX/FFTW3/include/fftw3.f03' ! ARA cluster
!  include '/usr/local/include/fftw3.f03'
!  include '/scratch/Saurabh/FFTW3/install/include/fftw3.f03'  ! nias
!  include '/home/me23jok/fftw-3.3.10/include/fftw3.f03' ! draco
end module



program TDSE_main

 use commandline_args
 use global_vars
 use data_au
 implicit none
 integer:: I, J, I_Emax
 integer :: scount,  & ! Starting "time"
           ecount ! Ending "time"
 integer  rate       ! number of clock ticks per second

 real(dp) st, ft,timer_elapsed_time
 real(dp) Emax
 real(dp), allocatable:: El(:), Al(:)
  
  call cpu_time(st)
  call system_clock(scount,rate)
  call read_command_line
  !call parse_command_line
  print*, "reading input:"
  print*, "General Inputs from ", trim(input)
!  call execute_command_line("pwd")
 
  call read_input
  call p_grid
  allocate(Vstates(Nstates), chi0(NR,guess_vstates,Nstates))
  chi0 = 0.d0
 allocate(El(Nt), Al(Nt)) 
  
!  call blas_check

!print*,"test"    
!  ewf = 0.d0
!  adb = 0.d0 
!
!  call adiabatic_surface(adb, mu_all) 
 
  call nuclear_wavefkt

  call pulse(El, Al)
  call propagation_1D(El, Al)
  deallocate (adb, mu_all, chi0)

 call cpu_time(ft)
 print*,'Run time=', ft-st, 'seconds' 
 call system_clock(ecount)
 timer_elapsed_time = real(ecount-scount,8)/real(rate,8)
 write(*,*) "Calculated run time is ",timer_elapsed_time," seconds"
  
end program


! _______________ Subroutines __________________________________________________


subroutine read_input

use global_vars
use pot_param
implicit none
 real:: dummy, dummy2, dummy3, dummy4
 character(200):: mk_out_dir
 integer:: i, j, input_tk
 namelist /grid/NR
 namelist /nucl_masses/m1,m2
 namelist /time_grid/dt,Nt
 namelist /elec_states/Nstates
 namelist /vib_states/guess_vstates
 namelist /ini_guess_wf/Ri, kappa
 namelist /laser_param/envelope_shape_laser1, envelope_shape_laser2, &
         & lambda1,lambda2,tp1,tp2,t_mid1,t_mid2,E01,E02,phi1,phi2, &
         & rise_time1, rise_time2
 namelist /input_files/input_data_dir,adb_pot, trans_dip_prefix, output_data_dir
 namelist /trans_dip_off/total_trans_off, trans_off
 namelist /absorber_choice/absorber
 namelist /ini_state/v_ini,N_ini,initial_distribution,temperature,kappa_tdse, RI_tdse
 namelist /parallelization/prop_par_FFTW,ITP_par_FFTW

 open(newunit=input_tk, file=input, status='old')
 read(input_tk, nml=grid)
 print*, "NR =", NR
 read(input_tk, nml=nucl_masses)
 print*, "masses:"
 print*, "m1 =", m1, "m2 =", m2
 read(input_tk,nml=time_grid)
 print*, "time grid:"
 print*, "dt =", dt, "fs"
 print*, "Nt =", Nt, "steps"
 read(input_tk,nml=elec_states)
 print*, "No. of electronic states:", Nstates
 read(input_tk,nml=vib_states)
 print*, "No. of maximum considered vibrational states:", guess_vstates
 read(input_tk,nml=ini_guess_wf)
 print*, "Guess wavefunction:" 
 print*, "Initial position (RI):", RI
 print*, "Initial width (kappa):", kappa
 read(input_tk,nml=laser_param)
 print*, "Laser parameters:"
 print*, "Laser #1:"
 print*, "Envelope shape:", trim(envelope_shape_laser1)
 print*, "Lambda:", lambda1, "nm"
 print*, "Electric field strength:", E01, "a.u."
 print*, "Pulse envelope: Cos**2"
 print*, "Pulse width (tp):", tp1, "fs"
 print*, "Pulse midpoint:", t_mid1, "fs"
 print*, "phi1:", phi1, "pi"
 print*, "Rise time:", rise_time1, "fs"
 print*, "Laser #2:"
 print*, "Envelope shape:", trim(envelope_shape_laser2)
 print*, "Lambda:", lambda2, "nm"
 print*, "Electric field strength:", E02, "a.u."
 print*, "Pulse envelope: Cos**2"
 print*, "Pulse width (tp):", tp2, "fs"
 print*, "Pulse midpoint:", t_mid2, "fs"
 print*, "phi1:", phi2, "pi"
 print*, "Rise time:", rise_time2, "fs"
 read(input_tk,nml=input_files)
 
 write(mk_out_dir, '(a)') adjustl(trim(output_data_dir))
 print*, "creating output directory ", trim(mk_out_dir)
 call execute_command_line("mkdir -p " // adjustl(trim(mk_out_dir)))

 allocate(R(NR), en(NR))
 allocate(pR(NR))
 allocate(adb(NR,Nstates),mu_all(Nstates,Nstates,NR))

 call pot_read

 read(input_tk,nml=trans_dip_off)
 print*, "Total trans dipole switched off:", total_trans_off
 call trans_dipole_read

 read(input_tk,nml=absorber_choice)
 print*, "Absorber function: ", absorber

 read(input_tk,nml=ini_state)
 print*, "TDSE Initial State:"
 print*, "Mode:", trim(initial_distribution)
 print*, "electronic state(s)", (N_ini-1)
 print*, "vibrational state(s)", (v_ini-1)
 print*, "Gaussian Distribution TDSE:" 
 print*, "centered at RI: ", RI_tdse
 print*, "standard deviation: ", kappa_tdse

 read(input_tk,nml=parallelization)
 print*, "FFTW Parallelization:"
 print*, "TDSE Propagation FFTW: ", trim(prop_par_FFTW)
 print*, "ITP FFTW: ", trim(ITP_par_FFTW)

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
  tp1 = tp1 / au2fs  
  tp2 = tp2 / au2fs  
  t_mid1 = t_mid1 / au2fs   
  t_mid2 = t_mid2 / au2fs   
  rise_time1 = rise_time1 / au2fs
  rise_time2 = rise_time2 / au2fs
  fwhm = (4._dp * log(2._dp)) / tp1**2 
  omega1=(1._dp / (lambda1 * 1.e-7_dp)) *cm2au
  omega2=(1._dp / (lambda2 * 1.e-7_dp))* cm2au
  phi1=phi1*pi
  phi2=phi2*pi
  print*,'_________________________'
  print*
  print*,'Final Parameters'
  print*,'_________________________'
  print*
  print*,'dt = ', SNGL(dt), 'a.u.'
  print*,'dR = ', SNGL(dR), 'a.u.'
  print*,'dPR = ', SNGL(dpR), 'a.u.'
  print*,'RI=', sngl(RI), 'a.u.'
  print*,'R0=', sngl(R0), 'a.u.', 'Rend=',sngl(Rend), 'a.u.'
  print*,'Wavelength 1 =', sngl(lambda1), 'nm'
  print*,'Phase 1 =', sngl(phi1)
  print*,'Field strength =', sngl(e01), 'a.u.', sngl(e01*e02au), 'V/m'
  print*,'Intensity =', sngl(e01**2*3.509e16_dp), 'W/cm2'
  print*,'Wavelength 2 =', sngl(lambda2), 'nm'
  print*,'Phase 2 =', sngl(phi2)
  print*,'Field strength =', sngl(e02), 'a.u.', sngl(e02*e02au), 'V/m'
  print*,'Intensity =', sngl(e02**2*3.509e16_dp), 'W/cm2'
  print*,'Wave duration =', sngl(tp1*au2fs), 'fs'
  print*
  print*, 'kap =', kap
  print*, 'lam =', lam
  print*
  print*,'__________________________'
  print*   
                 
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
        & input_data_dir, output_data_dir, dp
use data_au
 implicit none
 character(2000):: filepath
 integer:: I, pot_tk, pot_out_tk
 real(dp):: dummy
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
 !transition dipole moments of all states
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
