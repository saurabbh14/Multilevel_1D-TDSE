module commandline_args
character(10000):: command_line
character(2000):: input

!call read_command_line

!call parse_command_line


contains

     subroutine read_command_line
       integer :: exenamelength
       integer :: io, io2

       command_line = ""
       call get_command(command = command_line,status = io)
       if (io==0) then
         call get_command_argument(0,length = exenamelength,status = io2)
         if (io2==0) then
           command_line = "&cmd "//adjustl(trim(command_line(exenamelength+1:)))//" /"
         else
           command_line = "&cmd "//adjustl(trim(command_line))//" /"
         end if
       else
         write(*,*) io,"Error getting command line."
       end if
     end subroutine

     subroutine parse_command_line
       character(256) :: msg
       namelist /cmd/ input, adb_pot, trans_dip
       integer :: io

       if (len_trim(command_line)>0) then
         msg = ''
         read(command_line,nml = cmd,iostat = io,iomsg = msg)
         if (io/=0) then
           error stop "Error parsing the command line or cmd.conf " // msg
         end if
       end if
     end subroutine


end module commandline_args

module input_vars
 use commandline_args
! R-grid
 integer:: NR 

! electronic states
 integer:: Nstates

! vibrational states
 integer:: Vstates
 
! time grid 
 integer:: Nt

! masses
 double precision:: m1, m2

! guess initial wf
 double precision:: RI, kappa

! laser parameters
 double precision:: tp1, fwhm, t_start1
 double precision:: tp2, t_start2
 double precision:: e01, e02, phi1, phi2
 double precision:: lambda1, lambda2

! input files
 character(2000):: adb_pot, trans_dip_prefix

end module input_vars

module global_vars
 use input_vars
 integer, parameter:: Nx=2048 !, NR=1024 !changed back  !This program has been edited to check momentum grid dependance on number of grid point
 !reverted  back the changes while running calculations(Change date: 11/03/20)
 double precision:: dR
 double precision, allocatable:: R(:), x(:)
 double precision, allocatable:: en(:)
 double precision, allocatable:: Px(:),PR(:)
 double precision, allocatable:: Pot(:,:)
 double precision, allocatable, dimension(:,:,:):: mu_all
 double precision, allocatable, dimension(:,:):: adb
 double precision:: kap, lam
 double precision:: dx, dt, xeq
 double precision:: dpr, dpx 
 double precision:: omega1, omega2
 double precision:: mn, mn1, mn2 !all relevent mass veriables
end module

module data_au
 double precision,parameter:: au2a=0.52917706d0  ! length in a.u. to Angstrom
 double precision,parameter:: cm2au=4.5554927d-6 ! energy in wavenumbers to a.u.
 double precision,parameter:: au2fs=0.024        ! time in a.u. to femtoseconds
 double precision,parameter:: j2eV = 6.242D18    ! transforms energy in J to eV
 double precision,parameter:: au2eV=27.2116d0    ! energy in a.u. to eV
 double precision,parameter:: eV2nm=1239.84      ! energy from eV to wavelength(nm)
 double precision,parameter:: kB = 3.167d-6      ! Boltzmann constant hartree/K
 double precision,parameter:: i2au=2.0997496D-9
 double precision,parameter:: e02au=5.142206707e11
 double precision,parameter:: pi=3.141592653589793d0       ! that's just pi
 double precision,parameter:: mass=1836.15d0    ! reduced mass of deuterium
 double precision,parameter:: me =1.d0          ! mass of the electron
 double precision:: m_eff, m_red                ! effecitve mass of el. and nucl.
 double precision,parameter:: c_speed =137.03604 ! speed of light in a.u.
 complex*16,parameter:: im=(0.d0,1.d0) 
end module

module pot_param
use data_au
 double precision:: R0     ! Grid-Parameter, start..
 double precision,parameter:: x0 = -51.2d0*2*2!/au2a
 double precision::Rend   !..and end
 double precision,parameter:: xend = 51.2d0*2*2!/au2a
 double precision,parameter:: cpm=6.4d0*2*2!absorber position from the end of x-grid
 double precision,parameter:: cpmx=6.4d0*2*2*2!absorber position from the end of x-grid
 double precision,parameter:: cpmR=3.2d0*2 !*2 !absorber position from the end of R-grid

end module pot_param

module FFTW3
  use, intrinsic :: iso_c_binding
!   include '/usr/include/fftw3.f03'                                        ! Desktop packet
!   include '/home/me23jok/ProjectX/FFTW3/include/fftw3.f03' ! ARA cluster
  include '/usr/local/include/fftw3.f03'
!  include '/scratch/Saurabh/FFTW3/install/include/fftw3.f03'  ! nias
!  include '/home/me23jok/fftw-3.3.10/include/fftw3.f03' ! draco
end module



program TDSE_main

use commandline_args
use global_vars
use data_au
 implicit none
 integer:: I, J, I_Emax
 INTEGER :: scount,  & ! Starting "time"
           ecount ! Ending "time"
 integer  rate       ! number of clock ticks per second
 
 Real*4 st, ft,timer_elapsed_time
 double precision Emax
 double precision, allocatable, dimension(:,:,:):: ewf
 double precision, allocatable, dimension(:,:):: chi0
 double precision, allocatable:: El(:), Al(:)
  
  call cpu_time(st)
  call system_clock(scount,rate)
  call read_command_line
  call parse_command_line
  print*, "reading input:"
  print*, "General Inputs from", trim(input)
  
  call read_input
  call p_grid

  
 allocate(ewf(Nx,NR,Nstates))
 allocate(chi0(nr,vstates))
 allocate(El(Nt), Al(Nt)) 

print*,"test"
  call potential
print*,"test"    
  ewf = 0.d0
 ! adb = 0.d0 

!  call adiabatic_surface(adb, ewf, mu_all)  
  call nuclear_wavefkt(chi0)
  call pulse(El, Al)
  call propagation_1D(chi0, El, Al)
  deallocate (adb, mu_all)
 deallocate (ewf, chi0)

 call cpu_time(ft)
 print*,'Run time=', ft-st, 'seconds' 
 call system_clock(ecount)
 timer_elapsed_time = real(ecount-scount,8)/real(rate,8)
 write(*,*) "Calculated run time is ",timer_elapsed_time," seconds"

 
stop  
end program


! _______________ Subroutines __________________________________________________


subroutine read_input

use global_vars
use pot_param
implicit none
 real*4 :: dummy, dummy2, dummy3, dummy4
 integer:: i, j, input_tk
 namelist /grid/NR
 namelist /nucl_masses/m1,m2
 namelist /time_grid/dt,Nt
 namelist /elec_states/Nstates
 namelist /vib_states/Vstates
 namelist /ini_guess_wf/Ri, kappa
 namelist /laser_param/lambda1,lambda2,tp1,tp2,t_start1,t_start2,E01,E02,phi1,phi2
 namelist /input_files/adb_pot, trans_dip_prefix

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
 print*, "No. of vibrational states:", Vstates
 read(input_tk,nml=ini_guess_wf)
 print*, "Guess wavefunction:" 
 print*, "Initial position (RI):", RI
 print*, "Initial width (kappa):", kappa
 read(input_tk,nml=laser_param)
 print*, "Laser parameters:"
 print*, "Laser #1:"
 print*, "Lambda:", lambda1, "nm"
 print*, "Electric field strength:", E01, "a.u."
 print*, "Pulse envelope: Cos**2"
 print*, "Pulse width (tp):", tp1, "fs"
 print*, "Pulse midpoint (t_start):", t_start1, "fs"
 print*, "phi1:", phi1, "pi"
 print*, "Laser #2:"
 print*, "Lambda:", lambda2, "nm"
 print*, "Electric field strength:", E02, "a.u."
 print*, "Pulse envelope: Cos**2"
 print*, "Pulse width (tp):", tp2, "fs"
 print*, "Pulse midpoint (t_start):", t_start2, "fs"
 print*, "phi1:", phi2, "pi"
 read(input_tk,nml=input_files)


! open(10,file='input',status='old')
!  
!  read(10,*) NR                   ! NR = number of grid points
!  read(10,*) m1, m2               ! m1=mass1, m2=mass2
!  read(10,*) dt                   ! dt = time step in a.u.
!  read(10,*) Nt                   ! Nt = number of time steps.	
!  read(10,*) Nstates              ! Nstates = Number of calculated excited states.
!  read(10,*) Vstates              ! Vstates = Number of calculated vibrational states.
!  read(10,*) RI                   ! RI = center of initial Gaussian in 10^-10 m.
!  read(10,*) kappa                ! kappa = width of initial Gaussian 
!  read(10,*) lambda1, lambda2     ! wavelength of pulse
!  read(10,*) tp1, t_start1, E01, phi1 
!  read(10,*) tp2, t_start2, E02, phi2                  ! carrier-envelope phase
  
  allocate(R(NR), x(Nx), en(NR))
  allocate(pR(NR), px(Nx))
  allocate(adb(NR,Nstates),mu_all(Nstates,Nstates,NR))

  call pot_read
  call trans_dipole_read
!  R0=0.1/au2a
!  Rend=15.d0/au2a
!  R =R/au2a
  R0 =R(1)
  Rend = R(NR)
  dR =R(2)-R(1)!(Rend - R0) / (NR - 1)
  dx = (xend - x0) / (Nx - 1)
  dpx = (2.d0 * pi) / (dx * Nx)      
  dpr = (2.d0 * pi) / (dR * NR) 



 !Grids: co-ordinate space
!  do I =1,NR
!    R(I) = R0 + (I - 1) * dR
!  enddo

  do J=1,Nx
   x(J) = x0 + (J - 1) * dx
  enddo
!--------------------------
!Masses:
  m1=m1*mass
  m2=m2*mass 
  mn=m1+m2
  m_red = m1*m2/(m1+m2)
!  m_eff= (m1*me+m2*me+m1*m2)/(m1*m2*me)
  m_eff = (m1+m2) / (m1+m2+1.0d0)!(4.d0 * me * mass) / (2.d0 * mass + me) 
  mn1=m1/mn
  mn2=m2/mn
!----------------------------  
  kap=(mn+2.0)/(mn+1)
  lam=(m2-m1)/mn
  RI= RI / au2a   
 
  dt = dt / au2fs 
  tp1 = tp1 / au2fs  
  tp2 = tp2 / au2fs  
  t_start1 = t_start1 / au2fs   
  t_start2 = t_start2 / au2fs   
  fwhm = (4.d0 * dlog(2.d0)) / tp1**2 
  omega1=(1.d0 / (lambda1 * 1.d-7)) *cm2au
  omega2=(1.0d0/ (lambda2 *1.d-7))* cm2au
  phi1=phi1*pi
  phi2=phi2*pi
  print*,'_________________________'
  print*
  print*,'Parameters'
  print*,'_________________________'
  print*
  print*,'dt = ', SNGL(dt), 'a.u.'
  print*,'dx = ', SNGL(dx), 'a.u.'
  print*,'dR = ', SNGL(dR), 'a.u.'
  print*,'dpx = ', SNGL(dpx), 'a.u.'
  print*,'dPR = ', SNGL(dpR), 'a.u.'
  print*,'RI=', sngl(RI), 'a.u.'
  print*,'R0=', sngl(R0), 'a.u.', 'Rend=',sngl(Rend), 'a.u.'
  print*,'Wavelength 1 =', sngl(lambda1), 'nm'
  print*,'Phase 1 =', sngl(phi1)
  print*,'Field strength =', sngl(e01), 'a.u.', sngl(e01*e02au), 'V/m'
  print*,'Intensity =', sngl(e01**2*3.509e16), 'W/cm2'
  print*,'Wavelength 2 =', sngl(lambda2), 'nm'
  print*,'Phase 2 =', sngl(phi2)
  print*,'Field strength =', sngl(e02), 'a.u.', sngl(e02*e02au), 'V/m'
  print*,'Intensity =', sngl(e02**2*3.509e16), 'W/cm2'
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

use global_vars, only:px, pR, dpx, dpR, Nx, NR, R, x
use pot_param
 implicit none 
 integer:: I 
  
      
  do I = 1, Nx  
    if (I.le.(Nx / 2)) then    
    Px(I) = (I - 1) * dpx    
    else    
    Px(I) = - (Nx + 1 - I) * dpx    
    end if    
  end do
    print*, 'x0=', x0, 'xend=', xend
    print*, 'Px0=', Px((Nx/2)+1), 'Pxend=', Px(Nx/2) 
  
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
use global_vars, only:R, NR, adb, adb_pot
use data_au
 implicit none
 integer:: I, pot_tk
 double precision:: dummy
 print*, "Potential surfaces:", trim(adb_pot)
 open(newunit=pot_tk,file=adjustl(trim(adb_pot)),status='unknown')
 do I = 1, NR
   read(pot_tk,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
       ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
 end do
 close(pot_tk)
 open(1061,file='H2+_BO_read.out',status='unknown')
 do I = 1, NR
   write(1061,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
       ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
 end do
 close(1061)

end subroutine

!------------------------------------------------------------------------------

subroutine trans_dipole_read
use global_vars, only:R, NR, mu_all, trans_dip_prefix, Nstates
 implicit none
 integer:: I, L, M
 character(2000):: fn
 double precision:: dummy

 print*, "Transition dipoles with file prefix", trim(trans_dip_prefix)
 mu_all = 0.d0
 !transition dipole moments of all states
 if (trim(trans_dip_prefix) .eq.'') then
  do L = 1, Nstates
   do M = L+1, Nstates
     write(fn,fmt='(i0,i0,a)') L,M,'.dat'
     print*, trim(fn)
     open(unit=2000, file=adjustl(trim(fn)), form='formatted')
     do I = 1, NR
       read(2000,*) dummy, mu_all(L,M,I)
     enddo
     close(2000)
   enddo
  enddo
 
  do L = 1, Nstates
   do M = L+1, Nstates
     write(fn,fmt='(i0,i0,a)') L,M,'_read.out'
     print*, trim(fn)
     open(unit=2001, file=adjustl(trim(fn)), form='formatted')
     do I=1,NR
       write(2001,*) R(I), mu_all(L,M,I)
     enddo
     close(2001)
   enddo
  enddo

 else

  do L = 1, Nstates
   do M = L+1, Nstates
     write(fn,fmt='(a,i0,i0,a)') adjustl(trim(trans_dip_prefix)),L,M,'.dat'
     print*, trim(fn)
     open(unit=2000, file=adjustl(trim(fn)), form='formatted')
     do I = 1, NR
       read(2000,*) dummy, mu_all(L,M,I)
     enddo
     close(2000)
   enddo
  enddo
 
  do L = 1, Nstates
   do M = L+1, Nstates
     write(fn,fmt='(a,i0,i0,a)') adjustl(trim(trans_dip_prefix)),L,M,'_read.out'
     print*, trim(fn)
     open(unit=2001, file=adjustl(trim(fn)), form='formatted')
     do I=1,NR
       write(2001,*) R(I), mu_all(L,M,I)
     enddo
     close(2001)
   enddo
  enddo
 endif

end subroutine
