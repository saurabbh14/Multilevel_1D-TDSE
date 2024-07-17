module blas_interfaces_module

  use var_precision, only: wp=>dp

  implicit none 

  interface
  subroutine zgemv (trans, m, n, alpha,a,lda,x,incx,beta,y,incy)
    import wp
    character, intent(in):: trans
    integer, intent(in):: m
    integer, intent(in):: n
    complex(wp), intent(in)::  alpha
    complex(wp), dimension(*), intent(in):: a
    integer, intent(in):: lda
    complex(wp), dimension(*), intent(in):: x
    integer, intent(in):: incx
    complex(wp), intent(in):: beta
    complex(wp), dimension(*), intent(out):: y
    integer, intent(in):: incy
  end subroutine zgemv
  subroutine dgemv (trans, m, n, alpha,a,lda,x,incx,beta,y,incy)
    import wp
    character, intent(in):: trans
    integer, intent(in):: m
    integer, intent(in):: n
    real(wp), intent(in)::  alpha
    real(wp), dimension(*), intent(in):: a
    integer, intent(in):: lda
    real(wp), dimension(*), intent(in):: x
    integer, intent(in):: incx
    real(wp), intent(in):: beta
    real(wp), dimension(*), intent(out):: y
    integer, intent(in):: incy
  end subroutine dgemv
  subroutine write_matrix(a)
    import wp
    real(wp), dimension(:,:) :: a
  end subroutine write_matrix
  end interface
  
end module blas_interfaces_module

subroutine blas_check
use var_precision, only: wp=>dp
use blas_interfaces_module, only : zgemv, dgemv, write_matrix
integer:: I, J
 real(wp):: A(2,3), x(3), y(2)
 
 A = reshape((/1._wp,0._wp,-1._wp,-3._wp,2._wp,1._wp/),(/2,3/))
 call write_matrix(A)
 
 x = (/2._wp,1._wp,0._wp/)
 write(*,*)
 y = matmul(A,x)
 write(*,*) y
 write(*,*)
 y =0._wp
 call dgemv('N', 2, 3, 1._wp, A, size(A,dim=1), x, 1, 0._wp, y, 1)
 write(*,*) y
end subroutine blas_check

subroutine write_matrix(a)
use var_precision, only: wp=>dp
   real(wp), dimension(:,:) :: a
   write(*,*)

   do i = lbound(a,1), ubound(a,1)
      write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
   end do
end subroutine write_matrix



subroutine propagation_1D(E, A)

use global_vars
use pot_param
use data_au
use FFTW3
use omp_lib
use blas_interfaces_module, only : zgemv, dgemv

 implicit none   
! include "/usr/include/fftw3.f" 
 
 integer I, J, K,II
 integer L, M, N, void
 integer*8 planF, planB
 integer eR
 real(dp) dummy, dummy2, dummy3
 character(150):: filepath

 real(dp) dt2, time
 real(dp) c, sp, deR 
 real(dp):: normpn(Nstates), spec(NR,Nstates)
 real(dp) E(Nt), E1, norm(Nstates), E21, E22, norm_out(Nstates)
 real(dp):: norm_diss, norm_bound
 real(dp) A(Nt)
 real(dp) evR(Nstates), epr(Nstates),momt(Nstates), tot_momt
 real(dp) :: H_ac(NR,Nstates), Energy_axis(NR)
! double precision :: Boltzmann_populations(Vstates)
 real(dp) :: norm_overlap(Nstates),norm_outR(Nstates), norm_outP(Nstates)
 real(dp) :: norm_gesP(Nstates), norm_gP_over(Nstates)
 real(dp) :: vib_pop(guess_vstates,Nstates) 
 real(dp), allocatable, dimension(:):: cof
 real(dp), allocatable:: psi_Nstates_real(:), psi_Nstates_imag(:)
 complex(dp):: tout(Nstates,Nstates)
 complex(dp), allocatable, dimension(:,:):: psi_ges, psi_out
 complex(dp), allocatable, dimension(:):: psi_diss, psi_bound
 complex(dp), allocatable, dimension(:):: psi_Nstates, psi_Nstates1
 complex(dp), allocatable, dimension(:,:):: psi_loc, psi_ges1
 complex(dp), allocatable, dimension(:):: psi, kprop, kprop1
 complex(dp), allocatable, dimension(:,:):: psi_out1
 complex(dp), allocatable, dimension(:,:):: psi_outR, psi_outR1
 complex(dp), allocatable, dimension(:,:):: psi_gesP
 complex(dp), allocatable, dimension(:):: psi_chi

 integer:: chi0_tk, vstates_tk
 integer:: psi_1d_tk, cof_1d_tk
 integer:: dens_1d_tk, ex_dens_1d_tk
 integer:: avgR_1d_tk, avgpR_1d_tk
 integer:: norm_1d_tk, norm_pn_1d_tk, field_1d_tk
 integer:: accumulation_1d_tk, momt_1d_tk
 integer:: vibpop_1d_tk
 
!  write(filename4,fmt='(a,f4.2,a,f4.2,a,f4.2,a)') 'dens_1D_HeH+_mH',m1/mass,'_mHe',m2/mass,'_mr',mr,'.out'
!  open(200,file=filename4,status='unknown')
!  write(filename5,fmt='(a,f4.2,a,f4.2,a,f4.2,a)') 'ex_dens_1D_HeH+_mH',m1/mass,'_mHe',m2/mass,'_mr',mr,'.out'
!  open(201,file=filename5,status='unknown')
!  write(filename6,fmt='(a,f4.2,a,f4.2,a,f4.2,a)') 'R_1D_HeH+_mH',m1/mass,'_mHe',m2/mass,'_mr',mr,'.out'
!  open(800,file=filename6,status='unknown')
!  write(filename7,fmt='(a)') 'KER_all_spctra_HeH+_mH.out'
!  open(2222,file=filename7,status='unknown')
!  write(filename8,fmt='(a)') 'KER_all_ex_spctra_HeH+_mH.out'
!  open(2223,file=filename8,status='unknown')
          
 allocate(psi(NR),kprop(NR),psi_ges(NR,Nstates),cof(NR),kprop1(NR))
 allocate(psi_loc(nr,Nstates), psi_out(nr,Nstates))
 allocate(psi_Nstates(Nstates), psi_Nstates1(Nstates))
 allocate(psi_Nstates_real(Nstates), psi_Nstates_imag(Nstates))
 allocate(psi_out1(nr,Nstates), psi_outR(nr,Nstates),psi_gesP(nr,Nstates))
 allocate(psi_outR1(NR,Nstates))
 allocate(psi_chi(guess_vstates),psi_diss(NR),psi_bound(NR))

 print*
 print*,'Tuning FFTW...'

 void=fftw_init_threads()
 if (void .eq. 0) then
    print*, 'Error in fftw_init_threads, quitting'
 endif
 
 call fftw_plan_with_nthreads(omp_get_max_threads())    
 call dfftw_plan_dft_1d(planF, NR, psi, psi, FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_1d(planB, NR, psi, psi, FFTW_BACKWARD,FFTW_MEASURE)
  
             
 print*,'Done.'
 print*  
  
 write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "Bound-vibstates_in_Nthstates.out"
 open(newunit=vstates_tk,file=filepath,status='unknown')
 chi0 = 0._dp
 do N = 1, 1 
  read(vstates_tk,*) II, Vstates(N)
  write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "BO_Elelectronic-state-g", &
         & int(N-1), "_vibstates.out"
  open(newunit=chi0_tk,file=filepath,status='unknown')
  do I=1,NR
    read(chi0_tk,*) dummy, chi0(I,1:Vstates(N),N)
  enddo 
  close(chi0_tk)
 enddo
 close(vstates_tk)

 psi_ges = (0._dp,0._dp)    
 do I = 1, NR
!   psi_ges(I,1) = (0.55d0*chi0(I,1)+0.23d0*chi0(I,2)+0.11d0*chi0(I,3) &
!                 & + 0.07d0 * chi0(I,4) + 0.04d0*chi0(I,5)+ 2.0450413213653231E-002*chi0(I,6) &
!                 & + 1.4033977605843639E-002 *chi0(I,7)+ 1.0613089628800979E-002*chi0(I,8) &
!                 & + 8.8496818475779886E-003*chi0(I,9)+ 8.0611229210941892E-003*chi0(I,10))   !Thermal distribution over vibrational levels
   psi_ges(I,1) = chi0(I,1,1) !exp(kappa*(R(I)-RI)**2) !sqrt(0.55d0)*chi0(I,1)+sqrt(0.23d0)*chi0(I,2)+sqrt(0.11d0)*chi0(I,3) &
                 !& +sqrt( 0.07d0 )* chi0(I,4) + sqrt(0.04d0)*chi0(I,5) ! + sqrt(2.0450413213653231E-002)*chi0(I,6) &
!                 & +sqrt( 1.4033977605843639E-002) *chi0(I,7)+ sqrt(1.0613089628800979E-002)*chi0(I,8) &
!                 & +sqrt( 8.8496818475779886E-003)*chi0(I,9)+ sqrt(8.0611229210941892E-003)*chi0(I,10)   !Thermal distribution over vibrational levels

!   psi_ges(I,1)=chi0(i,1)!exp(kappa*(R(I)-RI)**2)
!   psi_ges(I,1) = sum(chi0(I,:)*sqrt(Boltzmann_populations(:)))
   kprop(I) = exp(-im *dt * PR(I)**2/(4._dp*m_red))  ! pR**2 /2 * red_mass UND Half time step
   kprop1(I) =exp(-im *dt * PR(I)**2/(2._dp*m_red))
 end do 
 
 ! cpm = 3.d0/ au2a
 ! call cutoff_cos(cpm,cof)
  call cutoff_ex(cof)

 
 call integ(psi_ges, norm) 
 print*,'Initial norm:', sngl(norm) 
  
 psi_ges(:,1) = psi_ges(:,1) / sqrt(norm(1))
 
 call integ(psi_ges, norm)
 print*,'norm after normalization:', sngl(norm) 
 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "psi0_1d.out"
 open(newunit=psi_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "cof_1d.out"
 open(newunit=cof_1d_tk,file=filepath,status='unknown') 
 do I = 1, NR
  write(psi_1d_tk,*) sngl(R(I)),sngl(abs(psi_ges(I,1))**2)
  write(cof_1d_tk,*) sngl(R(I)), sngl(cof(i))
 end do
close(psi_1d_tk)
close(cof_1d_tk)

           
!______________________________________________________________________
!                                                                     
!                   Propagation Loop                         
!______________________________________________________________________   

 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "norm_1d.out"
 open(newunit=norm_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "density_1d.out"
 open(newunit=dens_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "ex_density_1d.out"
 open(newunit=ex_dens_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "avgR_1d.out"
 open(newunit=avgR_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "avgPR_1d.out"
 open(newunit=avgpR_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "norm_pn_1d.out"
 open(newunit=norm_pn_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "field_1d.out"
 open(newunit=field_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "accumulation_1d.out"
 open(newunit=accumulation_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "momentum_1d.out"
 open(newunit=momt_1d_tk,file=filepath,status='unknown') 
 write(filepath,'(a,a)') adjustl(trim(output_data_dir)), 'vibpop1D_lambda.out'
 open(newunit=vibpop_1d_tk,file=filepath,status='unknown')

 print*
 print*,'1D propagation...'
 print*
    
 psi_out=(0._dp,0._dp)
 
timeloop: do K = 1, Nt

   time = K * dt
   epr = 0._dp
   evr = 0._dp
   psi=0._dp
   momt=0._dp
   norm_out=0._dp
   
   do J = 1,Nstates
     psi(:) = psi_ges(:,J)  ! Hilfsgroesse
     call dfftw_execute(planF)
     psi = psi * kprop
     call dfftw_execute(planB)
     psi = psi / dble(NR)
     psi_ges(:,J) = psi(:)      
   end do
   
   !$OMP PARALLEL DO DEFAULT(NONE) FIRSTPRIVATE(tout, psi_Nstates, psi_Nstates1) &
   !$OMP FIRSTPRIVATE(psi_Nstates_real, psi_Nstates_imag) &
   !$OMP SHARED(mu_all, E, psi_ges, K, Nstates, NR)
   do i = 1, NR
     call pulse2(tout, mu_all(:,:,I), E(K)) 
!     psi_ges(i,1:Nstates) = matmul(tout(1:Nstates,1:Nstates),psi_ges(i,1:Nstates))  
     psi_Nstates(:) = psi_ges(i,:)
     call zgemv('N', Nstates, Nstates, (1._dp, 0._dp), tout(1:Nstates,1:Nstates), &
             & size(tout,dim=1), psi_Nstates, 1, (0._dp,0._dp), psi_Nstates1, 1)
     psi_ges(i,:) = psi_Nstates1(:)
 !    call dgemv('N', Nstates, Nstates, 1._dp, tout(1:Nstates,1:Nstates), &
 !            & size(tout,dim=1), real(psi_Nstates), 1, 0._dp, psi_Nstates_real, 1)
 !    call dgemv('N', Nstates, Nstates, 1._dp, tout(1:Nstates,1:Nstates), &
 !            & size(tout,dim=1), aimag(psi_Nstates), 1, 0._dp, psi_Nstates_imag, 1)
 !    psi_ges(i,:) = cmplx(psi_Nstates_real(:),psi_Nstates_imag(:),kind=dp)
   end do
   !$OMP END PARALLEL DO
    

   do j = 1, Nstates   
!     psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * (adb(i,j)+kap*(mu_all(I,J,J) &
!              & +0.5d0*R(I))*E(K)+((2-kap)*mr+kap-1)*R(I)*E(K)))!+H_ac(i,j))) !         
     psi_ges(:,j) = psi_ges(:,j) * exp(-im * dt * adb(:,j)) !+0.8d0*R(I)*E(K)))!+H_ac(i,j))) !         
   end do
    
   do J = 1,Nstates
     psi(:) = psi_ges(:,J)  ! Hilfsgroesse
     call dfftw_execute(planF) 
     psi = psi * kprop 
     psi = psi /sqrt(dble(nr))     
     epr(j) =  sum(abs(psi(:))**2 * pr(:))
     epr(j) = epr(j) * dr
     call dfftw_execute(planB)
     psi = psi / sqrt(dble(NR))
     psi_ges(:,J) = psi(:)      
   end do     
 

!-----------------------------------------

  
   do N = 1, Nstates
     evR(N) = evR(N) + sum(abs(psi_ges(:,N))**2 * R(:))
   enddo
   evR = evR * dR
   psi_loc(:,1) = 1._dp/sqrt(2._dp)*(psi_ges(:,1) + psi_ges(:,2))
   psi_loc(:,2) = 1._dp/sqrt(2._dp)*(psi_ges(:,1) - psi_ges(:,2))

   call integ(psi_ges, norm)    
   call integ(psi_loc, normPn)
 
   do j = 1,Nstates
     if (norm(j).ge.1.e-8_dp) then
       evR(j) = evR(j)/norm(j)
       epr(j) = epr(j)/norm(j)
     end if
   end do
    
! ------------ Popolation analysis in vibrational states ----------------------
   vib_pop = 0._dp
   do N = 1, 1 !Nstates
     do L=1,vstates(N)
       psi_chi(L)=0._dp
       psi_chi(L)=sum(chi0(:,L,N) * (psi_ges(:,N)))
       psi_chi(L)=psi_chi(L)*dR
       vib_pop(L,N)=abs(psi_chi(L)**2)
     enddo
   enddo

   write(vibpop_1d_tk,*) sngl(time*au2fs), vib_pop(1:vstates(1),1) 
     
   write(avgR_1d_tk,*) sngl(time *au2fs), sngl(evR)
   write(avgpR_1d_tk,*) sngl(time *au2fs), sngl(epr)
   write(norm_1d_tk,*) sngl(time *au2fs), norm
   write(norm_pn_1d_tk,*) sngl(time *au2fs), sngl(normPn)
   write(field_1d_tk,*) sngl(time *au2fs), sngl(E(K))      

   if(mod(K,100).eq.0) then
    do I = 1, NR, 4   
      write(dens_1d_tk,*) sngl(time *au2fs), sngl(R(I)), sngl(abs(psi_ges(I,1)**2))
      write(ex_dens_1d_tk,*) sngl(time *au2fs), sngl(R(I)), sngl(abs(psi_ges(I,2:Nstates)**2))
    end do 
    write(dens_1d_tk,*)
    write(ex_dens_1d_tk,*)
  
   end if  
 


 ! -------------Continuum treatment ---------------------
!    do j=1,Nstates
!     do i=1,nr
!      psi(i)=psi_out(i,j)
!     end do
!
!      call dfftw_execute(planB)
!      psi=psi/sqrt(dble(nr))
!     do i=1,NR
!      psi_out(i,j)=psi(i)
!     enddo
!   enddo
!
!   do i=1,NR
!    tout(1,1) = cos(-kap*mu_all(I,1,2)*E(K)*dt)
!    tout(1,2) = -im*sin(-kap*mu_all(I,1,2)*E(K)*dt)
!    tout(2,1) = tout(1,2)
!    tout(2,2) = tout(1,1)
!    ! call pulse2(tout, mu_all(i,:,:),E(K))
!     psi_out(i,1:Nstates)=matmul(tout(1:Nstates,1:Nstates), psi_out(i,1:Nstates))
!   enddo
!   do j=1,Nstates
!     do i=1,nr
!      psi(i)=psi_out(i,j)
!     end do
!
!      call dfftw_execute(planF)
!      psi=psi/sqrt(dble(nr))
!     do i=1,NR
!      psi_out(i,j)=psi(i)
!     enddo
!   enddo
!
!   do j=1,Nstates
!    do i=1,NR
!     psi_out(i,j)=psi_out(i,j)*kprop1(I)
!    enddo
!   enddo
   
   do J = 1, Nstates
     psi_outR(:,J) = psi_outR(:,J) * kprop1(:)

     psi_outR1(:,J) = psi_ges(:,J) * (1._dp-cof(:))
     psi_ges(:,J) = psi_ges(:,J) * cof(:) ! psi_ges = psi_nondiss
   enddo

!   norm_overlap(1)=2.0d0*dot_product(psi_outR(:,1),psi_ges(:,1))*dr
!   norm_overlap(2)=2.0d0*dot_product(psi_outR(:,2),psi_ges(:,2))*dr
!
!   call integ(psi_outR,norm_outR)
!   norm1=norm1+norm_overlap(1)+norm_overlap(2)+norm_outR(1)+norm_outR(2)
!
   do J=1,Nstates
     psi(:)=psi_outR1(:,j)
     call dfftw_execute(planF)
     psi=psi/sqrt(dble(nr))
     psi_outR1(:,J)=psi(:)
   enddo
   psi_outR = psi_outR + psi_outR1

!   psi_gesP=psi_ges
!   do j=1,Nstates
!     do i=1,nr
!      psi(i)=psi_gesP(i,j)
!     end do
!      call dfftw_execute(planF)
!      psi=psi/sqrt(dble(nr))
!     do i=1,NR
!      psi_gesP(i,j)=psi(i)
!     enddo
!   enddo
!   norm_gP_over(1)=2.0d0*dot_product(psi_outR(:,1),psi_gesP(:,1))*dr
!   norm_gP_over(2)=2.0d0*dot_product(psi_outR(:,2),psi_gesP(:,2))*dr
!
!!   write(1000,*) sngl(time*au2fs), sngl(norm_overlap), sngl(norm1), sngl(norm_gesP),sngl(norm_gP_over)
!
!    call integ(psi_outR, norm_outP)
!
!    do j=1,Nstates
!  !   if (norm_outP(J).ge.3.d-8) then
!        psi_out(:,j)=psi_out(:,j)+psi_outR(:,j)!/sqrt(norm_outP(J))
!  !   end if
!   enddo
!   norm_out=0.00d0
!
!   do J=1,Nstates
!     do i=1,NR
!        norm_out(J)=norm_out(J)+abs(psi_out(I,J))**2
!     enddo
!    enddo
!      norm_out= norm_out*dr
!
!   do J=1,Nstates
!      psi_out(:,J)=psi_out(:,J)!/sqrt(norm_out(J))
!   enddo
!
!   do I=1,NR
!     momt(1)=momt(1)+pr(I)*abs(psi_out(I,1))**2
!     momt(2)=momt(2)+pr(I)*abs(psi_out(I,2))**2
!   enddo
!!   momt=momt*dr
!    momt(1)=momt(1)/(sum(abs(psi_out(:,1))**2))
!    momt(2)=momt(2)/(sum(abs(psi_out(:,2))**2))
!


!   write(998,*) sngl(time* au2fs),sngl((momt(1))),sngl((momt(2))), norm_outP

! ------------   

end do timeloop
!  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), &
!          & 'Momentum_spectra_1d.out'
!  open(momt_spec_1d_tk,file=filename2,status='unknown')
!  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), &
!          & 'Momentum_spectra_total_1D.out'
!  open(1122,file=filename2,status='unknown')
!  write(filename2,fmt='(a)') 'KER_spctra_1g_1D.out'
!  open(1114,file=filename2,status='unknown')
!  write(filename3,fmt='(a)') 'KER_spctra_2u_1D.out'
!  open(1115,file=filename3,status='unknown')
!  write(filename2,fmt='(a)') 'KER_spctra_total_1D.out'
!  open(1144,file=filename2,status='unknown')
  
!Final vibrastional poppulation in ground state
  do N =1, 1
    print*, "Final vibrational poppulation in the ground state"
    psi_bound = 0._dp
    do J = 1,vstates(N)
    psi_chi(J) = 0._dp
    do I = 1, NR
     psi_chi(J) = psi_chi(J)+ chi0(I,J,N) * (psi_ges(I,1))
    enddo
     psi_chi(J) = psi_chi(J)*dR
    do I = 1,NR
     psi_bound(I) = psi_bound(I)+psi_chi(J)*chi0(I,J,N)
    enddo
    norm_bound = sum(abs(psi_bound(:))**2)*dR
    print*, 'Vibpop (',J,') =',norm_bound
    enddo
  enddo

!  psi_diss(:)=psi_ges(:,1)-psi_bound(:)
!  psi = psi_diss
!  call dfftw_execute(planF)
!  psi=psi/sqrt(dble(NR))
!  psi_diss=psi
!  norm_diss=sum(abs(psi_diss(:))**2)*dR
!  write(*,*) "Dissociation yield 1g =",norm_diss
!  write(*,*) lambda1,norm_diss
!  psi_diss=psi_diss/sqrt(norm_diss)
!   do I=NR/2 +1,NR
!     write(1112,*) pR(I),abs(psi_diss(I))**2
!   enddo
!   do I=1,NR/2
!     write(1112,*) pR(I),abs(psi_diss(I))**2
!     write(1114,*) (pR(I)**2/(2*m_red)),abs(psi_diss(I))**2
!   enddo
!
!  psi(:)=psi_ges(:,2)
!  call dfftw_execute(planF)
!  psi=psi/sqrt(dble(NR))
!  psi_ges(:,2)=psi(:)
!  norm_diss=sum(abs(psi_ges(:,2))**2)*dR
!  write(*,*) "Dissociation yield 2u =",norm_diss
!  write(*,*) lambda1, norm_diss
!  psi_ges(:,2)=psi_ges(:,2)/sqrt(norm_diss)
!   do I=NR/2 +1,NR
!     write(1113,*) pR(I),abs(psi_ges(I,2))**2
!   enddo
!   do I=1,NR/2
!     write(1113,*) pR(I),abs(psi_ges(I,2))**2
!     write(1115,*) (pR(I)**2/(2*m_red)),abs(psi_ges(I,2))**2
!   enddo
!  
!!   psi_diss(:)=psi_diss(:)+psi_ges(:,2)
!  
!  psi(:)=psi_diss(:)+ psi_ges(:,2)
!  call dfftw_execute(planF)
!  psi=psi/sqrt(dble(NR))
!  norm_diss=sum(abs(psi(:))**2)*dR
!  write(*,*) "Dissociation yield 1g =",norm_diss
!  psi=psi/sqrt(norm_diss)
!   do I=NR/2 +1,NR
!     write(1122,*) pR(I),abs(psi(I))**2
!   enddo
!   do I=1,NR/2
!     write(1122,*) pR(I),abs(psi(I))**2
!     write(1144,*) (pR(I)**2/(2*m_red)),abs(psi(I))**2
!   enddo
  

!   do I=NR/2 +1,NR
!     write(2222,*) 0.1*(II-1), pR(I), (pR(I)**2)/(2*m_red),abs(psi_diss(I))**2
!   enddo
!   do I=1,NR/2
!     write(2222,*) mr, pR(I), (pR(I)**2)/(2*m_red),abs(psi_diss(I))**2
!   enddo
!   write(2222,*)
!   do I=NR/2 +1,NR
!     write(2223,*) 0.1*(II-1), pR(I), (pR(I)**2)/(2*m_red),abs(psi_ges(I,2))**2
!   enddo
!   do I=1,NR/2
!     write(2223,*) mr, pR(I), (pR(I)**2)/(2*m_red),abs(psi_ges(I,2))**2
!   enddo
!   write(2223,*)
   
!  call integ(psi_out, norm_out)
!  psi_out(:,1)=psi_out(:,1)/sqrt(norm_out(1))
!  psi_out(:,2)=psi_out(:,2)/sqrt(norm_out(2))
!  spec=0.0d0
!  sp=0.00d0
!  tot_momt=0.0d0
!  do J=1,Nstates
!    do I=1,NR
!      tot_momt=tot_momt+pr(I)*abs(psi_out(I,J))**2
!    enddo
!  enddo
!    tot_momt=tot_momt/(sum(abs(psi_out(:,:))**2))
!
!
!
!   do I=NR/2+1,NR
!    spec(I,1)=abs(psi_out(I,1))**2
!    spec(I,2)=abs(psi_out(I,2))**2
! !   write(999,*) sngl(pr(I)), sngl(spec(I,1)), sngl(spec(I,2))
!   enddo
!
!   do I=1,NR/2
!    spec(I,1)=abs(psi_out(I,1))**2
!    spec(I,2)=abs(psi_out(I,2))**2
!    write(999,*) sngl(pr(I)), sngl(spec(I,1)), sngl(spec(I,2))
!   enddo
!  print*, 'average momentum=', sngl(tot_momt), 'a.u.'
!  print*, 'average velocity=', sngl((tot_momt*au2a)/(mass*au2fs)), 'A°/fs'


 
! ------------   



 

 close(100, status='keep')
 close(101, status='keep') 
 close(800, status='keep')
 close(801, status='keep')
 close(906, status='keep')
 close(908, status='keep')   
 close(909, status='keep')
 close(200, status='keep')
 close(201, status='keep')
 close(999, status='keep')
! close(998, status='keep')
!close(1000) 

                                      
 call dfftw_destroy_plan(planF)
 call dfftw_destroy_plan(planB)
   
 deallocate(psi, kprop, psi_ges, cof, psi_loc, psi_out,psi_outR, psi_gesP)

return
end subroutine


!_________________________________________________________


subroutine integ(psi, norm)
      
use global_vars, only:NR, Nstates, dR, dp
 implicit none
 integer I, J
      
 real(dp) norm(Nstates)
 complex(dp) psi(NR,Nstates)
      
 norm = 0.d0
 
 do J = 1, Nstates
  do I = 1, NR
    norm(J)= norm(J) + abs(psi(I,J))**2   
  end do
 end do
 
   norm = norm * dR
     
     
return 
end subroutine


!______________________________________________________________

!subroutine light_ind_pot(mu, w, d, pott)
!
!use data_grid
!
!implicit none
! double precision:: w, u, d, mu(3), pott(2)
! dimension u(2,2), d(2)
!
!
!u(1,1) = pott(1)-mu(1)*w
!u(1,2) = -mu(3)*w
!
!u(2,1) = -mu(3)*w
!u(2,2) = pott(2)-mu(2)*w
!
!
!
! call jacobi(u,2,d)
! 
! 
!return
!end subroutine
!
!!------------------------------------------

subroutine pulse2(tout,mu,E)

use global_vars, only:dt, Nstates,kap, lam, dp
use data_au, only:im

implicit none

 integer:: i, J
 real(dp):: w, u, d, mu(Nstates,Nstates), q, E
 integer Info, Lwork

 complex(dp):: tout, b, z
 dimension:: u(Nstates,Nstates), d(Nstates), b(Nstates,Nstates), z(Nstates,Nstates), tout(Nstates,Nstates)
 real(dp) work(1000), u1(Nstates, Nstates)
 character(len=1):: JOBZ

    
!u(1,1) = 0.d0
!u(1,2) = -kap*mu(1,2) * E

!u(2,1) = -kap*mu(2,1) * E
!u(2,2) = 0.d0

!Diapole matrix
u=0._dp
do I=1, Nstates-1
 do J=I+1, Nstates
       u(I,J)= -kap*mu(I,J) * E
       u(J,I)= -kap*mu(I,J) * E
 enddo
enddo
!print*, "u"
!write( * , * ) ((u(i,j),j=1,Nstates), i=1,Nstates )


Lwork=-1
JOBZ='V'
 call dsyev(JOBZ,'U', Nstates, u,Nstates, d, work, Lwork,info )
! call jacobi(u,Nstates,d)
 Lwork = min(1000, int(work(1)))
!     Solve eigenproblem.

 call dsyev('V', 'U', Nstates, u, Nstates, d, work, Lwork, info)

 if( info.GT.0 ) then
     write(*,*)'The algorithm failed to compute eigenvalues.'
     stop
 endif
!print*,"eigen vector u"
!write( * , * ) ( (u(i,j),j=1,Nstates), i=1,Nstates )
b= (0._dp,0._dp)
 
do J = 1,Nstates
  b(J,J) = exp(-im * dt * d(J))
end do
!print*, "b"
!write( * , * ) ((b(i,j),j=1,Nstates), i=1,Nstates )

z = matmul(u,b)
!print*, "z"
!write( * , * ) ((z(i,j),j=1,Nstates), i=1,Nstates )
tout = matmul(z,transpose(u))
!print*, "tout"
!write( * , * ) ((tout(i,j),j=1,Nstates), i=1,Nstates )



return
end subroutine

!!------------------------------------------------
subroutine cutoff_cos(cof)
use global_vars, only:NR, R, dR, dp
use data_au
use pot_param

implicit none

   integer :: J
   real(dp):: cof(NR)

   do j = 1, NR
    R(J) = R0 + (j - 1) * dR
    if(R(J).lt.(Rend - cpmR)) then
    cof(j) = 1.d0
    else
    cof(j) = dcos(((R(J) - Rend + cpmR) / -cpmR) * (0.5d0 * pi))
    cof(j) = cof(j)**2
    end if
   end do

 return
 end subroutine
!------------------------------------------------
 subroutine cutoff_ex(cof)
 use global_vars, only:NR, R, dR, dp
 use data_au
 use pot_param

  implicit none

  integer :: J
  real(dp):: cof(NR),c

  c=1.800d0
 do j = 1, NR
   R(J) = R0 + (j - 1) * dR
   cof(j)=1.0d0/(1.0d0+exp(c*(R(J)-Rend+cpmR)))
 end do

 return
 end subroutine

!!-------------------------
!
!
!      subroutine jacobi (mat,dim,ewerte)
!
!         implicit       none
!
!         real*8         genau
!         parameter      (genau=1.d-15)
!
!         integer        Jmax,mmax
!         parameter      (Jmax=15,mmax=18)
!         integer        matdim
!         parameter      (matdim=2)
!
!         real*8         mat(matdim,matdim)
!         integer        dim
!         real*8         ewerte(matdim)
!
!         real*8         s(matdim,matdim)
!         integer        ca,cb,p,q
!         real*8         c1,c2,t1,t2,t3,v1,v2,v3
!         real*8         tmp,l,n,t,m1,w,m
!         logical        flag
!
!         s= 0.d0
!         !!!!call fillmt(s,dim,dim,0.d0,matdim,matdim)
!
!         do 1 ca=1,dim,1
!            s(ca,ca)=1.d0
!1           continue
!
!         l=0.d0
!         do 2 ca=2,dim,1
!            do 2 cb=1,dim,1
!               tmp=mat(ca,cb)
!               l=l+2.d0*tmp*tmp
!2              continue
!
!         n=dsqrt(l)
!         m=genau*n/dim
!         t=n
!
!3        t=t/dim
!4           do 6 q=2,dim,1
!               do 6 p=1,q-1,1
!                  flag=.false.
!                  if (dabs(mat(p,q)).gt.t) then
!                     flag=.true.
!                     v1=mat(p,p)
!                     v2=mat(p,q)
!                     v3=mat(q,q)
!                     m1=(v1-v3)/2.d0
!                     if (m1.eq.0.d0) then
!                           w=-1.d0
!                        else
!                           if (m1.gt.0.d0) then
!                                 w=-v2/(dsqrt(v2*v2+m1*m1))
!                              else
!                                 w=v2/(dsqrt(v2*v2+m1*m1))
!                              endif
!                        endif
!
!                     t1=w/dsqrt(2.d0*(1+dsqrt(1.d0-w/2.d0)))
!                     t2=t1*t1
!                     c1=dsqrt(1.d0-t2)
!                     c2=c1*c1
!                     t3=t1*c1
!
!                     do 7 ca=1,dim,1
!                        l=mat(ca,p)*c1-mat(ca,q)*t1
!                        mat(ca,q)=mat(ca,p)*t1+mat(ca,q)*c1
!                        mat(ca,p)=l
!                        l=s(ca,p)*c1-s(ca,q)*t1
!                        s(ca,q)=s(ca,p)*t1+s(ca,q)*c1
!                        s(ca,p)=l
!7                       continue
!                     do 8 ca=1,dim,1
!                        mat(p,ca)=mat(ca,p)
!                        mat(q,ca)=mat(ca,q)
!8                       continue
!                     mat(p,p)=v1*c2+v3*t2-2*v2*t3
!                     mat(q,q)=v1*t2+v3*c2+2*v2*t3
!                     tmp=(v1-v3)*t3+v2*(c2-t2)
!                     mat(p,q)=tmp
!                     mat(q,p)=tmp
!                     end if
!6                 continue
!               if (flag) go to 4
!            if (m.lt.t) go to 3
!ewerte=0.d0
!         !!!call fillvc(ewerte,dim,0.d0)
!         do 9 ca=1,dim,1
!            ewerte(ca)=mat(ca,ca)
!9           continue
!         do 10 ca=1,dim,1
!            do 10 cb=1,dim,1
!               mat(ca,cb)=s(ca,cb)
!10             continue
!
!         return
!         end




