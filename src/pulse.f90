! A subroutine for defining the field 
subroutine pulse(E2,A)

use global_vars
use data_au
use FFTW3
use omp_lib

implicit none

integer k, void, Nt2
integer*8 planTF, planTB, planTF2, planTB2
character(150) filename
real(dp) time
real(dp) E2(Nt), E21(Nt), E22(Nt), E(Nt)
!double precision, allocatable:: E_dum(:), F_dum(:)
complex(dp), allocatable:: E_dum(:), F_dum(:)
complex(dp), allocatable:: E2_dum(:), F2_dum(:)
real(dp) A01, A02
real(dp) A2(Nt), A21(Nt), A22(Nt)
real(dp) A(Nt), E2_new(Nt)
real(dp) IR(Nt), Eeff(Nt), IREeff(Nt), Eeff2(Nt)
real(dp) cos2, trapazoidal, gaussian ! envelope shapes
real(dp) g1(Nt), g2(Nt), en_curve
real(dp) time_end, time_start 
real(dp) pulse_offset, rise_time
real(dp) dummy

! file tokens
integer:: elec_field_tk, vec_field_tk
integer:: field1_tk, field2_tk
integer:: envelope1_tk, envelope2_tk
integer:: IR_Eeff_field_tk
integer:: IR_field_spec_tk, IR_field_time_tk
integer:: Eeff_field_spec_tk, Eeff_field_time_tk

write(filename,fmt='(a,a,f4.2,a)') adjustl(trim(output_data_dir)), &
        & 'Total_electric_field_phi', phi2/pi,'pi.out'
open(newunit=elec_field_tk, file=filename,status="unknown")
write(filename,fmt='(a,a,f4.2,a)') adjustl(trim(output_data_dir)), &
        & 'Total_vector_field_phi', phi2/pi, 'pi.out'
open(newunit=vec_field_tk, file=filename,status="unknown")
write(filename,fmt='(a,a,f4.2,a,i0,a)') adjustl(trim(output_data_dir)), &
        & 'electric_field1_E', E01,'_width',Int(tp1*au2fs),'.out'
open(newunit=field1_tk, file=filename,status="unknown")
write(filename,fmt='(a,a,f6.4,a,i0,a)') adjustl(trim(output_data_dir)), &
        & 'electric_field2_E', E02,'_width',Int(tp2*au2fs),'.out'
open(newunit=field2_tk, file=filename,status="unknown")
write(filename,fmt='(a,a)') adjustl(trim(output_data_dir)), 'envelope1.out'
open(newunit=envelope1_tk, file=filename,status="unknown")
write(filename,fmt='(a,a)') adjustl(trim(output_data_dir)), 'envelope2.out'
open(newunit=envelope2_tk, file=filename,status="unknown")
write(filename,fmt='(a,a)') adjustl(trim(output_data_dir)), 'IR+2Eeff.dat'
open(newunit=IR_Eeff_field_tk, file=filename,status="unknown")

 print*
 print*,'Tuning FFTW...'
 void=fftw_init_threads( )
 if (void==0) then
    write(*,*) 'Error in fftw_init_threads, quitting'
    stop
 endif

Nt2 = Nt - 1623
 
allocate(E_dum(Nt2), F_dum(Nt2))
allocate(E2_dum(Nt), F2_dum(Nt))

call fftw_plan_with_nthreads(omp_get_max_threads())
call dfftw_plan_dft_1d(planTF, Nt2, E_dum, E_dum, FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_plan_dft_1d(planTB, Nt2, E_dum, E_dum, FFTW_BACKWARD, FFTW_ESTIMATE)
call dfftw_plan_dft_1d(planTF2, Nt, E2_dum, F2_dum, FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_plan_dft_1d(planTB2, Nt, F2_dum, E2_dum, FFTW_BACKWARD, FFTW_ESTIMATE)


time_end = Nt*dt
time_start = (t_mid2- tp2/2)
pulse_offset = 0 !5*pi/omega2
rise_time = 5 !fs
rise_time = rise_time/au2fs

tp1=tp1/(1-2/pi)
E21 =0._dp
E22 =0._dp
g1=0._dp
g2=0._dp
A01=E01/omega1
A02=E02/omega2
E2_New=0._dp

 envelope_shape_laser1 = trim(envelope_shape_laser1)
 envelope_shape_laser2 = trim(envelope_shape_laser2)
 select case(envelope_shape_laser1)
  case("cos2")
   do K = 1, Nt
    time = K*dt 
    g1(k) = cos2(time,tp1,t_mid1,pulse_offset)
   enddo
  case("gaussian")
   do K = 1, Nt
    time = K*dt
    g1(K) = gaussian(time, tp1, t_mid1)
   enddo
  case("trapazoidal")
    do K = 1, Nt
      time = K*dt
      g1(K) = trapazoidal(time, tp1, t_mid1, rise_time)
    enddo
  case default
   print*, "Laser1: Default pulse shape is CW."
 end select

 select case(envelope_shape_laser2)
  case("cos2")
   do K = 1, Nt
    time = k*dt 
    g2(k) = cos2(time,tp1,t_mid2,pulse_offset)
   enddo
  case("gaussian")
   do K = 1, Nt
    time = K*dt
    g1(K) = gaussian(time, tp2, t_mid2)
   enddo
  case("trapazoidal")
    do K = 1, Nt
      time = K*dt
      g2(K) = trapazoidal(time, tp2, t_mid2, rise_time)
    enddo
  case default
   print*, "Laser1: Default pulse shape is CW."
 end select
 
 

timeloop: do K = 1, Nt
  time = k*dt 
  E21(K) = E01*g1(K) * cos(omega1 * (time-t_mid1-pulse_offset)+phi1)
  A21(K) = A01*g1(K) * cos(omega1 * (time-t_mid1-pulse_offset)+phi1)
  write(field1_tk,*) time*au2fs, E21(K), A21(K)
  write(envelope1_tk,*) time*au2fs, g1(K)

  E22(K) = E02*g2(K)* sin(omega2 * (time-t_mid2-tp2/2-rise_time)+phi2)
  A22(K) = A02*g2(K)* sin(omega2 * (time-t_mid2-tp2/2-rise_time)+phi2)
  write(field2_tk,*) time*au2fs, E22(K), A22(K)
  write(envelope2_tk,*) time*au2fs, g2(K)


    E2(K)=E21(K)+E22(K)
    A2(K)=A21(K)+A22(K)
    E(K) = E2(K)
  !  E2(K)=0.0d0
    write(elec_field_tk,*) time*au2fs, E2(K), A2(K)
!! if (K.ge. 1800 .and. K.le. 4000) then
! if (K .ge. 1961 .and. K.le. 4000) then
!    read(IR_Eeff_field_tk,*) dummy, IR(K), Eeff(K), IREeff(K), Eeff2(K)
!!    E2(K) = IREeff(K)
! !   E2(K) = 2*Eeff(K)
! else
!    E2(K) = 0.d0
! endif
 !vector field
  A(k)=(-1._dp)*sum(E2(1:K))*dt
  if (K .gt. 1) then
  E2_New(k)=-(A2(k)-A2(k-1))/dt
  endif
  write(vec_field_tk,*) time*au2fs, A(k), E2_New(k), A2(K)
enddo timeloop 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%% IR field spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write(filename,fmt='(a,a,i0,a)') adjustl(trim(output_data_dir)), &
        & 'IR-field_spec_lambda',int(lambda1),'nm.out'
open(newunit=IR_field_spec_tk, file=filename,status="unknown")
write(filename,fmt='(a,a,i0,a)') adjustl(trim(output_data_dir)), &
        & 'IR-field_time_lambda',int(lambda1),'nm.out'
open(newunit=IR_field_time_tk, file=filename,status="unknown")
E2_dum = E
call dfftw_execute(planTF2,E2_dum, F2_dum)
F2_dum = F2_dum/sqrt(dble(Nt))
call dfftw_execute(planTB, F2_dum, E2_dum)
E2_dum = E2_dum/sqrt(dble(Nt))

do K = Nt/2+1, Nt
write(IR_field_spec_tk,*) -(Nt+1-K)*2*pi/(dt*Nt), real(F2_dum(K)), imag(F2_dum(K)), abs(F2_dum(K))
enddo
do K = 1, Nt/2
write(IR_field_spec_tk,*) (K-1)*2*pi/(dt*Nt), real(F2_dum(K)), imag(F2_dum(K)), abs(F2_dum(K))
enddo
close(IR_field_spec_tk)

do K=1,Nt
write(IR_field_time_tk,*) K*dt*au2fs, real(E2_dum(K)), imag(E2_dum(K)), abs(E2_dum(K))
enddo
close(IR_field_time_tk)

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%% E_eff spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!write(filename,fmt='(a,i0,a)') 'Eeff-field_spec_lambda',int(lambda1),'nm.dat'
!open(newunit=Eeff_field_spec_tk, file=filename,status="unknown")
!write(filename,fmt='(a,i0,a)') 'Eeff-field_time_lambda',int(lambda1),'nm.dat'
!open(newunit=Eeff_field_time_tk, file=filename,status="unknown")
!E_dum(:) = E2(1624:Nt)
!call dfftw_execute(planTF,E_dum, E_dum)
!F_dum = E_dum/sqrt(dble(Nt2))
!E_dum = E_dum/sqrt(dble(Nt2))
!call dfftw_execute(planTB, E_dum, E_dum)
!E_dum = E_dum/sqrt(dble(Nt2))
!
!do K = Nt2/2+1, Nt2
!write(Eeff_field_spec_tk,*) -(Nt2+1-K)*2*pi/(dt*Nt2), real(F_dum(K)), imag(F_dum(K)), abs(F_dum(K))
!enddo
!do K = 1, Nt2/2
!write(Eeff_field_spec_tk,*) (K-1)*2*pi/(dt*Nt2), real(F_dum(K)), imag(F_dum(K)), abs(F_dum(K))
!enddo
!close(Eeff_field_spec_tk)
!
!do K=1,Nt2
!write(Eeff_field_time_tk,*) K*dt*au2fs, real(E_dum(K)), imag(E_dum(K)), abs(E_dum(K))
!enddo
!close(Eeff_field_time_tk)
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E2 = E ! IR field

call dfftw_destroy_plan(planTF)
call dfftw_destroy_plan(planTB)
call dfftw_destroy_plan(planTF2)
call dfftw_destroy_plan(planTB2)

end subroutine
 
!------------------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%% pulse envelope functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cos2(time, tp, t_mid, pulse_offset)
 use global_vars, only: Nt, dp
 use data_au, only: pi
 implicit none 
 real(dp):: time, tp, t_mid, pulse_offset
 real(dp):: cos2
 if (time .gt. (t_mid+pulse_offset-tp/2) .and. time .lt. (t_mid+pulse_offset+tp/2)) then
   cos2 = cos((time - t_mid-pulse_offset)*pi/tp)**2      
 else
   cos2 = 0._dp
 endif
end function

function trapazoidal(time, tp, t_mid, rise_time)
 use global_vars, only: Nt, dp
 use data_au, only: pi
 real(dp):: time, tp, t_mid, rise_time
 real(dp):: trapazoidal, slope, yc
 
  if (time .ge. t_mid-(tp/2 + rise_time) .and. time .le. t_mid-tp/2) then
     slope = 1._dp/rise_time
     yc = (t_mid -(tp/2 + rise_time)) * slope
     trapazoidal = slope * time - yc
  elseif (time .gt. t_mid-tp/2 .and. time .le. t_mid+tp/2) then
     trapazoidal = 1._dp
  elseif (time .gt. t_mid+tp/2 .and. time .le. t_mid+(tp/2+rise_time)) then 
     slope = -1._dp/rise_time
     yc = (t_mid +tp/2 + rise_time) * slope
     trapazoidal = slope * time - yc
  else
     trapazoidal = 0._dp
  endif

end function

function gaussian(time, tp, t_mid)
 use global_vars, only: Nt, dp 
 implicit none
 real(dp):: time, tp, t_mid
 real(dp):: gaussian, fwhm

  fwhm = (4._dp * log(2._dp)) / tp**2
  gaussian = exp(-fwhm * (time - t_mid)**2)
 
end function
