module data_grid
 integer:: Nt, Nstates, Vstates, NR
 integer, parameter:: Nx=2048 !, NR=1024 !changed back  !This program has been edited to check momentum grid dependance on number of grid point
 !reverted  back the changes while running calculations(Change date: 11/03/20)
 double precision:: RI, kappa,dR
 double precision, allocatable:: R(:), x(:)
 double precision, allocatable:: en(:)
 double precision, allocatable:: Px(:),PR(:)
 double precision, allocatable:: Pot(:,:)
 double precision:: kap, lam
 double precision:: dx, dt, xeq
 double precision:: dpr, dpx 
 double precision:: tp1, fwhm, t_start1
 double precision:: tp2, t_start2
 double precision:: e01, e02, phi1, phi2
 double precision:: omega1, omega2
 double precision:: lambda1, lambda2
 double precision:: mn, mn1, mn2, m1, m2 !all relevent mass veriables
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



program metiu

use data_grid
use data_au
 implicit none
 integer:: I, J, I_Emax
 INTEGER :: scount,  & ! Starting "time"
           ecount ! Ending "time"
 integer  rate       ! number of clock ticks per second
 
 Real*4 st, ft,timer_elapsed_time
 double precision Emax
 double precision, allocatable, dimension(:,:,:):: mu_all
 double precision, allocatable, dimension(:,:):: adb
 double precision, allocatable, dimension(:,:,:):: ewf
 double precision, allocatable, dimension(:,:):: chi0
 double precision, allocatable:: El(:), Al(:)
  
  call cpu_time(st)
  call system_clock(scount,rate)
  call input
  call p_grid

  
 allocate(adb(NR,Nstates),ewf(Nx,NR,Nstates),mu_all(Nstates,Nstates,NR))
 allocate(chi0(nr,vstates))
 allocate(El(Nt), Al(Nt)) 

print*,"test"
  call potential
print*,"test"    
  ewf = 0.d0
  adb = 0.d0 

!  call adiabatic_surface(adb, ewf, mu_all)  
  call pot_read(adb)
  call trans_dipole_read(mu_all)
  call nuclear_wavefkt(adb,chi0)
  call pulse(El, Al)
  call propagation_1D(adb, mu_all, chi0, El, Al)
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


subroutine input

use data_grid
use pot_param
implicit none
 real*4 :: dummy, dummy2, dummy3, dummy4
 integer:: i, j

 open(10,file='input',status='old')
  
  read(10,*) NR                   ! NR = number of grid points
  read(10,*) m1, m2               ! m1=mass1, m2=mass2
  read(10,*) dt                   ! dt = time step in a.u.
  read(10,*) Nt                   ! Nt = number of time steps.	
  read(10,*) Nstates              ! Nstates = Number of calculated excited states.
  read(10,*) Vstates              ! Vstates = Number of calculated vibrational states.
  read(10,*) RI                   ! RI = center of initial Gaussian in 10^-10 m.
  read(10,*) kappa                ! kappa = width of initial Gaussian 
  read(10,*) lambda1, lambda2     ! wavelength of pulse
  read(10,*) tp1, t_start1, E01, phi1 
  read(10,*) tp2, t_start2, E02, phi2                  ! carrier-envelope phase
  
  allocate(R(NR), x(Nx), en(NR))
  allocate(pR(NR), px(Nx))
! open(11,file='aeN_new_512.dat',status='old')
 open(11,file='H2+_softcore-param_smooth-interpolated_2048.dat',status='old')
 do i=1,NR
 read(11,*) R(I), en(I)
! read(11,*) dummy, dummy2 !, dummy3, dummy4
! read(11,*) dummy, dummy2 !, dummy3, dummy4
! read(11,*) dummy, dummy2 !, dummy3, dummy4
 enddo
 close(11)
!  R0=0.1/au2a
!  Rend=15.d0/au2a
  R =R/au2a
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

use data_grid
use pot_param
use data_au
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
subroutine pot_read(adb)
use data_grid
use pot_param
use data_au
 implicit none
 integer:: I
 double precision:: adb(NR, Nstates), dummy
 open(106,file='H2+_BO.dat',status='unknown')
 do I = 1, NR
   read(106,*) dummy, adb(I,:) !, sngl(adb(I,2)*au2eV), &
       ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
 end do
 close(106)
 open(1061,file='H2+_BO_read.out',status='unknown')
 do I = 1, NR
   write(1061,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
       ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
 end do
 close(1061)

end subroutine

!------------------------------------------------------------------------------

subroutine trans_dipole_read(mu_all)
use data_grid
use pot_param
use data_au
 implicit none
 integer:: I, L, M
 character(150):: fn
 double precision:: mu_all(Nstates,Nstates,NR), dummy

 mu_all = 0.d0
 !transition dipole moments of all states
 do L = 1, Nstates
  do M = L+1, Nstates
    write(fn,fmt='(i0,i0,a)') L,M,'.dat'
    print*, fn
    open(unit=2000, file=fn, form='formatted')
    do I = 1, NR
      read(2000,*) dummy, mu_all(L,M,I)
    enddo
    close(2000)
  enddo
 enddo

 do L = 1, Nstates
  do M = L+1, Nstates
    write(fn,fmt='(i0,i0,a)') L,M,'_read.out'
    print*, fn
    open(unit=2001, file=fn, form='formatted')
    do I=1,NR
      write(2001,*) R(I), mu_all(L,M,I)
    enddo
    close(2001)
  enddo
 enddo

end subroutine
