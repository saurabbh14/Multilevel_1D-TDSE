module blas_interfaces_module

  use var_precision, only: wp=>dp, idp

  implicit none 

  interface
  subroutine zgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    import wp, idp
    character*1, intent(in):: transa
    character*1, intent(in):: transb
    integer(idp), intent(in):: m
    integer(idp), intent(in):: n
    integer(idp), intent(in):: k
    complex(wp), intent(in) :: alpha
    complex(wp), dimension(*), intent(in):: a
    integer(idp), intent(in):: lda
    complex(wp), dimension(*), intent(in):: b
    integer(idp), intent(in):: ldb
    complex(wp), intent(in):: beta
    complex(wp), dimension(*), intent(out):: c
    integer(idp), intent(in):: ldc
  end subroutine zgemm
  subroutine dgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    import wp, idp
    character*1, intent(in):: transa
    character*1, intent(in):: transb
    integer(idp), intent(in):: m
    integer(idp), intent(in):: n
    integer(idp), intent(in):: k
    real(wp), intent(in) :: alpha
    real(wp), dimension(*), intent(in):: a
    integer(idp), intent(in):: lda
    real(wp), dimension(*), intent(in):: b
    integer(idp), intent(in):: ldb
    real(wp), intent(in):: beta
    real(wp), dimension(*), intent(out):: c
    integer(idp), intent(in):: ldc
  end subroutine dgemm
  subroutine zgemv (trans, m, n, alpha,a,lda,x,incx,beta,y,incy)
    import wp, idp
    character*1, intent(in):: trans
    integer(idp), intent(in):: m
    integer(idp), intent(in):: n
    complex(wp), intent(in)::  alpha
    complex(wp), dimension(*), intent(in):: a
    integer(idp), intent(in):: lda
    complex(wp), dimension(*), intent(in):: x
    integer(idp), intent(in):: incx
    complex(wp), intent(in):: beta
    complex(wp), dimension(*), intent(out):: y
    integer(idp), intent(in):: incy
  end subroutine zgemv
  subroutine dgemv (trans, m, n, alpha,a,lda,x,incx,beta,y,incy)
    import wp, idp
    character*1, intent(in):: trans
    integer(idp), intent(in):: m
    integer(idp), intent(in):: n
    real(wp), intent(in)::  alpha
    real(wp), dimension(*), intent(in):: a
    integer(idp), intent(in):: lda
    real(wp), dimension(*), intent(in):: x
    integer(idp), intent(in):: incx
    real(wp), intent(in):: beta
    real(wp), dimension(*), intent(out):: y
    integer(idp), intent(in):: incy
  end subroutine dgemv
  subroutine dsyev ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO) 
    import wp, idp
    character*1:: JOBZ
    character*1:: UPLO
    integer(idp):: N
    real(wp), dimension(*):: A
    integer(idp):: LDA
    real(wp), dimension(*):: W
    real(wp), dimension(*):: WORK
    integer(idp):: LWORK
    integer(idp):: INFO
  end subroutine dsyev
  subroutine write_matrix(a)
    import wp
    real(wp), dimension(:,:) :: a
  end subroutine write_matrix
  end interface
  
end module blas_interfaces_module

subroutine blas_check
use var_precision, only: wp=>dp, idp
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
 call dgemv('N', 2_idp, 3_idp, 1._wp, A, size(A,dim=1,kind=idp), x, 1_idp, 0._wp, y, 1_idp)
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
 
 integer I, J, K,I_cpmR, II
 integer L, M, N, void
 integer(idp) planF, planB
 integer eR
 real(dp) dummy, dummy2, dummy3
 character(150):: filepath

 real(dp) dt2, time
 real(dp) c, sp, deR 
 real(dp):: normpn(Nstates), spec(NR,Nstates)
 real(dp) E(Nt), E1, norm(Nstates), E21, E22
 real(dp):: norm_outR(Nstates), norm_outR1(Nstates), norm_SE_outR1(Nstates)
 real(dp):: norm_diss(Nstates), norm_bound(Nstates)
 real(dp) A(Nt)
 real(dp) evR(Nstates), epr(Nstates),momt(Nstates), tot_momt
 real(dp) :: H_ac(NR,Nstates), Energy_axis(NR)
! double precision :: Boltzmann_populations(Vstates)
 real(dp) :: norm_overlap(Nstates),norm_outP(Nstates)
 real(dp) :: norm_gesP(Nstates), norm_gP_over(Nstates)
 real(dp) :: vib_pop(guess_vstates,Nstates) 
 real(dp), allocatable, dimension(:):: cof, V_abs
 real(dp), allocatable:: psi_Nstates_real(:), psi_Nstates_imag(:)
 complex(dp):: tout(Nstates,Nstates)
 real(dp), allocatable:: SE_outR1(:,:)
 complex(dp), allocatable, dimension(:):: exp_abs, abs_func
 complex(dp), allocatable, dimension(:,:):: psi_ges, psi_out, psi_ges_p
 complex(dp), allocatable, dimension(:,:):: psi_diss, psi_bound
 complex(dp), allocatable, dimension(:):: psi_Nstates, psi_Nstates1
 complex(dp), allocatable, dimension(:,:):: psi_loc, psi_ges1
 complex(dp), allocatable, dimension(:):: psi, kprop, kprop1
 complex(dp), allocatable, dimension(:,:):: psi_out1
 complex(dp), allocatable, dimension(:,:):: psi_outR, psi_outR1
 complex(dp), allocatable, dimension(:,:):: psi_gesP
 complex(dp), allocatable, dimension(:):: psi_chi

 integer:: chi0_tk, vstates_tk
 integer:: psi_1d_tk, cof_1d_tk, complex_abs_tk
 integer:: dens_1d_tk, ex_dens_1d_tk, Pdens_1d_tk
 integer:: avgR_1d_tk, avgpR_1d_tk
 integer:: norm_1d_tk, norm_pn_1d_tk, field_1d_tk
 integer:: accumulation_1d_tk, momt_1d_tk
 integer:: vibpop_1d_tk
 integer:: psi_outR_norm_1d_tk, psi_outR_Pdens_1d_tk
 integer:: KER_spectra_tk, momt_spectra_tk
 integer:: KER_spectra_un_tk, momt_spectra_un_tk
 
 allocate(psi(NR),kprop(NR),psi_ges(NR,Nstates),cof(NR),kprop1(NR))
 allocate(exp_abs(NR), V_abs(NR), psi_ges_p(NR,Nstates),abs_func(NR))
 allocate(psi_loc(nr,Nstates))
 allocate(psi_Nstates(Nstates), psi_Nstates1(Nstates))
 allocate(psi_Nstates_real(Nstates), psi_Nstates_imag(Nstates))
 allocate(psi_outR(NR,Nstates),psi_gesP(NR,Nstates))
 allocate(psi_outR1(NR,Nstates))
 allocate(psi_chi(guess_vstates),psi_diss(NR,Nstates),psi_bound(NR,Nstates))

 print*
 print*,'Tuning FFTW...'

 void=fftw_init_threads()
 if (void .eq. 0) then
    print*, 'Error in fftw_init_threads, quitting'
 endif
 
! call fftw_plan_with_nthreads(1)    
 call fftw_plan_with_nthreads(omp_get_max_threads())    
 call dfftw_plan_dft_1d(planF, int(NR,kind=idp), psi, psi, FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_1d(planB, int(NR,kind=idp), psi, psi, FFTW_BACKWARD,FFTW_MEASURE)
  
             
 print*,'Done.'
 print*  
  
 write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "Bound-vibstates_in_Nthstates.out"
 open(newunit=vstates_tk,file=filepath,status='unknown')
 chi0 = 0._dp
 do N = 1, Nstates 
  read(vstates_tk,*) II, Vstates(N)
  write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "BO_Electronic-state-g", &
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

!   psi_ges(I,1)=exp(kappa*(R(I)-RI)**2)
!   psi_ges(I,1) = sum(chi0(I,:)*sqrt(Boltzmann_populations(:)))
   kprop(I) = exp(-im *dt * pR(I)*pR(I) /(4._dp*m_red))  ! pR**2 /2 * red_mass UND Half time step
   kprop1(I) =exp(-im *dt * pR(I)*pR(I) /(2._dp*m_red))
 end do 
 
 ! cpm = 3.d0/ au2a
 ! call cutoff_cos(cpm,cof)

 I_cpmR = minloc(abs(R(:)-cpmR),1) - 50
 print*, "I_cpmR =", I_cpmR, ", NR-I_cpmR", NR-I_cpmR
 print*, "R(NR-I_cpmR) =", R(NR-I_cpmR)

 select case(absorber)
 case ("CAP")
   call Complex_absorber_function(V_abs, exp_abs)
   abs_func = exp_abs
 case("mask")
   call mask_function_ex(cof)
   abs_func = cof
 end select
 
 call integ(psi_ges, norm) 
 print*,'Initial norm:', norm 
  
 psi_ges(:,1) = psi_ges(:,1) / sqrt(norm(1))
 
 call integ(psi_ges, norm)
 print*,'norm after normalization:', norm 
 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "psi0_1d.out"
 open(newunit=psi_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "absorber_function.out"
 open(newunit=cof_1d_tk,file=filepath,status='unknown') 
 do I = 1, NR
  write(psi_1d_tk,*) R(I), abs(psi_ges(I,1))**2
  write(cof_1d_tk,*) R(I), abs_func(I)
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
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "Pdensity_1d.out"
 open(newunit=Pdens_1d_tk,file=filepath,status='unknown') 
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
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "psi_outR_norm_1d.out"
 open(newunit=psi_outR_norm_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "psi_outR_Pdensity_1d.out"
 open(newunit=psi_outR_Pdens_1d_tk,file=filepath,status='unknown') 

 print*
 print*,'1D propagation...'
 print*
    
 psi_outR = (0._dp,0._dp)
 norm = 0._dp
 normPn = 0._dp
 norm_outR = 0._dp
 
timeloop: do K = 1, Nt

   time = K * dt
   epr = 0._dp
   evr = 0._dp
   psi=0._dp
   momt=0._dp
!   norm_out=0._dp
   norm = 0._dp
   normPn = 0._dp
   norm_outR = 0._dp
   
   do J = 1,Nstates
     psi = (0._dp, 0._dp)
     psi(:) = psi_ges(:,J)  ! Hilfsgroesse
     call dfftw_execute_dft(planF, psi, psi)
     psi = psi * kprop
     call dfftw_execute_dft(planB, psi, psi)
     psi = psi / dble(NR)
     psi_ges(:,J) = psi(:)      
   end do
   
   do j = 1, Nstates   
!     psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * (adb(i,j)+kap*(mu_all(I,J,J) &
!              & +0.5d0*R(I))*E(K)+((2-kap)*mr+kap-1)*R(I)*E(K)))!+H_ac(i,j))) !         
     psi_ges(:,j) = psi_ges(:,j) * exp(-im * 0.5_dp * dt * (adb(:,j))) !+0.8d0*R(I)*E(K)))!+H_ac(i,j))) !         
   end do
   
!   !$OMP PARALLEL DO DEFAULT(NONE) FIRSTPRIVATE(tout, psi_Nstates, psi_Nstates1) &
!   !$OMP FIRSTPRIVATE(psi_Nstates_real, psi_Nstates_imag) &
 !  !$OMP SHARED(mu_all, E, psi_ges, K, Nstates, NR, dt)
   do i = 1, NR
     call pulse2(tout, mu_all(:,:,I), E(K)) 
!     psi_ges(i,1:Nstates) = matmul(tout(1:Nstates,1:Nstates),psi_ges(i,1:Nstates))  
     psi_Nstates(:) = psi_ges(i,:)
     call zgemv('N', int(Nstates,kind=idp), int(Nstates, kind=idp), (1._dp, 0._dp),  &
             & tout, size(tout,dim=1,kind=idp), psi_Nstates, 1_idp, (0._dp,0._dp), &
             & psi_Nstates1, 1_idp)
     psi_ges(i,:) = psi_Nstates1(:)
   end do
!   !$OMP END PARALLEL DO
    

   do j = 1, Nstates   
!     psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * (adb(i,j)+kap*(mu_all(I,J,J) &
!              & +0.5d0*R(I))*E(K)+((2-kap)*mr+kap-1)*R(I)*E(K)))!+H_ac(i,j))) !         
     psi_ges(:,j) = psi_ges(:,j) * exp(-im * 0.5_dp * dt * (adb(:,j))) !+0.8d0*R(I)*E(K)))!+H_ac(i,j))) !         
   end do
    
   do J = 1,Nstates
     psi = (0._dp, 0._dp)
     psi(:) = psi_ges(:,J)  ! Hilfsgroesse
     call dfftw_execute_dft(planF, psi, psi) 
     psi = psi * kprop 
     psi = psi /sqrt(dble(nr))     
     epr(j) =  sum(abs(psi(:))**2 * pr(:))
     epr(j) = epr(j) * dr
     psi_ges_p(:,J) = psi(:)
     call dfftw_execute_dft(planB, psi, psi)
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
       vib_pop(L,N)=abs(psi_chi(L))**2
     enddo
   enddo

   write(vibpop_1d_tk,*) time*au2fs, vib_pop(1:vstates(1),1) 
   write(avgR_1d_tk,*) time *au2fs, evR
   write(avgpR_1d_tk,*) time *au2fs, epr
   write(norm_1d_tk,*) time *au2fs, norm
   write(norm_pn_1d_tk,*) time *au2fs, normPn
   write(field_1d_tk,*) time *au2fs, E(K)      
  
   ! Coordinate space density 
   if(mod(K,100).eq.0) then
    do I = 1, NR, 4   
      write(dens_1d_tk,*) time *au2fs, R(I), abs(psi_ges(I,1))**2
      write(ex_dens_1d_tk,*) time *au2fs, R(I), abs(psi_ges(I,2:Nstates))**2
    end do 
    write(dens_1d_tk,*)
    write(ex_dens_1d_tk,*)
   end if
   
   ! Momentum space density
   if (mod(K,100) .eq. 0) then
    do I = NR/2 +1, NR, 4   
      write(Pdens_1d_tk,*) time *au2fs, pR(I), abs(psi_ges_p(I,:))**2
    end do 
    do I = 1, NR/2, 4   
      write(Pdens_1d_tk,*) time *au2fs, pR(I), abs(psi_ges_p(I,:))**2
    end do 
    write(Pdens_1d_tk,*)
   endif

 ! -------------Continuum treatment ---------------------
  
   psi_outR1 = (0._dp, 0._dp)
   do J = 1, Nstates
     psi_outR(:,J) = psi_outR(:,J) * kprop1(:) * exp(-im*dt*adb(NR-I_cpmR,J))
     psi_outR1(:,J) = psi_ges(:,J) * (1._dp-abs_func(:))
     psi_ges(:,J) = psi_ges(:,J) * abs_func(:) ! psi_ges = psi_nondiss
   enddo
   write(psi_outR_norm_1d_tk,*) time*au2fs, norm_outR(:)
   do J = 1, Nstates
     psi = (0._dp, 0._dp)
     psi(:) = psi_outR1(:,J)
     call dfftw_execute_dft(planF, psi, psi)
     psi = psi/sqrt(dble(NR))
     psi_outR1(:,J) = psi(:)
   enddo
   psi_outR = psi_outR + psi_outR1
   ! Momentum space density
   if (mod(K,100) .eq. 0) then
    do I = NR/2 +1, NR, 4
      write(psi_outR_Pdens_1d_tk,*) time *au2fs, pR(I), abs(psi_outR(I,:))**2
    end do
    do I = 1, NR/2, 4
      write(psi_outR_Pdens_1d_tk,*) time *au2fs, pR(I), abs(psi_outR(I,:))**2
    end do 
    write(psi_outR_Pdens_1d_tk,*)
   endif

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
  
!Final vibrational population in ground state
  do N =1, Nstates
    print*
    print*, "Final vibrational poppulation in the elec. state", int(N-1)
    print*, "Number of vibrational states:", vstates(N)
    psi_bound = 0._dp
    do J = 1,vstates(N)
      psi_chi(J) = 0._dp
      do I = 1, NR
        psi_chi(J) = psi_chi(J)+ chi0(I,J,N) * (psi_ges(I,N))
      enddo
      psi_chi(J) = psi_chi(J)*dR
      do I = 1,NR
        psi_bound(I,N) = psi_bound(I,N)+psi_chi(J)*chi0(I,J,N)
      enddo
      norm_bound(N) = sum(abs(psi_bound(:,N))**2)*dR
      print*, 'Vibpop (',J,') =',norm_bound(N)
    enddo
    print*, 'Total population in state', int(N-1), ":", sum(abs(psi_ges(:,N))**2)*dR
    print*, 'Bound population in state', int(N-1), ":", sum(abs(psi_bound(:,N))**2)*dR
    psi_diss(:,N)=psi_ges(:,N)-psi_bound(:,N)
    print*, 'Unbound population in state', int(N-1), ":", sum(abs(psi_diss(:,N))**2)*dR
            !& sum(abs(psi_ges(:,N)-psi_bound(:,N))**2)*dR
    print*, 'Unbound population reached to the absorber in state', int(N-1), ":", &
            sum(abs(psi_outR(:,N))**2)*dR
    print*, 'Calculating KER spectra in state ', int(N-1)
    psi(:) = psi_diss(:,N)
    call dfftw_execute_dft(planF, psi, psi)
    psi=psi/sqrt(dble(NR))
    psi_diss(:,N)=psi(:)
    norm_diss(N)=sum(abs(psi_diss(:,N))**2)*dR
    norm_outR(N) = sum(abs(psi_outR(:,N))**2)*dR
    print*, "Writing KER spectra in state ", int(N-1)  
    write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "KER_spectra_from_state_g", &
           &  int(N-1), "_unnormalized.out"
    open(newunit=KER_spectra_un_tk,file=filepath,status='unknown')
    write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "momt_spectra_from_state_g",&
           &  int(N-1), "_unnormalized.out"
    open(newunit=momt_spectra_un_tk,file=filepath,status='unknown')
    write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "KER_spectra_from_state_g", &
           &  int(N-1), ".out"
    open(newunit=KER_spectra_tk,file=filepath,status='unknown')
    write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "momt_spectra_from_state_g",&
           &  int(N-1), ".out"
    open(newunit=momt_spectra_tk,file=filepath,status='unknown')
    do I=NR/2 +1,NR
      write(momt_spectra_un_tk,*) pR(I), abs(psi_diss(I,N))**2, abs(psi_outR(I,N))**2
      write(momt_spectra_tk,*) pR(I), abs(psi_diss(I,N)/sqrt(norm_diss(N)))**2, &
              & abs(psi_outR(I,N)/sqrt(norm_diss(N)))**2
    enddo
    do I=1,NR/2
      write(momt_spectra_un_tk,*) pR(I), abs(psi_diss(I,N))**2, abs(psi_outR(I,N))**2
      write(momt_spectra_tk,*) pR(I), abs(psi_diss(I,N)/sqrt(norm_diss(N)))**2, &
              & abs(psi_outR(I,N)/sqrt(norm_outR(N)))**2
      write(KER_spectra_un_tk,*) pR(I)**2/(2*m_red), m_red*abs(psi_diss(I,N))**2 / pR(I),&
              & m_red*abs(psi_outR(I,N))**2 / pR(I) 
      write(KER_spectra_tk,*) pR(I)**2/(2*m_red), m_red*abs(psi_diss(I,N)/sqrt(norm_diss(N)))**2 / pR(I), &
              & m_red*abs(psi_outR(I,N)/sqrt(norm_outR(N)))**2 / pR(I) 
    enddo
    close(momt_spectra_un_tk)
    close(KER_spectra_un_tk)
    close(momt_spectra_tk)
    close(KER_spectra_tk)
  enddo

  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Total_KER_spectra_unnormalized.out"
  open(newunit=KER_spectra_un_tk,file=filepath,status='unknown')
  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Total_momt_spectra_unnormalized.out"
  open(newunit=momt_spectra_un_tk,file=filepath,status='unknown')
  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Total_KER_spectra.out"
  open(newunit=KER_spectra_tk,file=filepath,status='unknown')
  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Total_momt_spectra.out"
  open(newunit=momt_spectra_tk,file=filepath,status='unknown')
  do I=NR/2 +1,NR
    write(momt_spectra_un_tk,*) pR(I), sum(abs(psi_diss(I,:))**2), sum(abs(psi_outR(I,:))**2)
    write(momt_spectra_tk,*) pR(I), sum(abs(psi_diss(I,:)/sqrt(norm_diss(:)))**2), &
              & sum(abs(psi_outR(I,:)/sqrt(norm_outR(:)))**2)
  enddo
  do I=1,NR/2
    write(momt_spectra_un_tk,*) pR(I), sum(abs(psi_diss(I,:))**2), sum(abs(psi_outR(I,:))**2)
    write(momt_spectra_tk,*) pR(I), sum(abs(psi_diss(I,:)/sqrt(norm_diss(:)))**2), &
              & sum(abs(psi_outR(I,:)/sqrt(norm_outR(:)))**2)
    write(KER_spectra_un_tk,*) pR(I)**2/(2*m_red), m_red*sum(abs(psi_diss(I,:))**2) / pR(I), &
              & m_red*sum(abs(psi_outR(I,:))**2) / pR(I)
    write(KER_spectra_tk,*) pR(I)**2/(2*m_red), m_red*sum(abs(psi_diss(I,:)/sqrt(norm_diss(:)))**2) /pR(I), &
              & m_red*sum(abs(psi_outR(I,:)/sqrt(norm_outR(:)))**2) / pR(I)
  enddo
  close(momt_spectra_un_tk)
  close(KER_spectra_un_tk)
  close(momt_spectra_tk)
  close(KER_spectra_tk)

 
! ------------   
                                      
 call dfftw_destroy_plan(planF)
 call dfftw_destroy_plan(planB)
   
 deallocate(psi, kprop, psi_ges, cof, psi_loc, psi_outR1,psi_outR, psi_gesP)

return
end subroutine


!_________________________________________________________

subroutine complex_absorber_function(v_abs, f)
use global_vars, only: NR, dp, dR, dt, R
use data_au, only: im
use pot_param, only: cpmR
 implicit none
 integer I
 real(dp):: a, eps, V_abs(NR), n, R0, p
 complex(dp):: iV_abs(NR), f(NR)
 
 eps = epsilon(a) 
 print*, "Lower limit of the precision:", eps
 n = 4 ! power of absorber function
 R0 = R(NR)- cpmR ! start of the absorber
 p = 20._dp ! optimal absorption momentum
 a = -log(eps) *(n+1) *p / (2*(R(NR)-R0)**(n+1))
 print*, "Absorber prefactor a:", a

 do I = 1, NR
  if (R(I) .gt. abs(R0)) then
    V_abs(I) = a*(R(I)-R0)**n
  else
    V_abs(I) = 0._dp
  endif
 enddo
 f(:) = exp(-dt *V_abs(:))
end subroutine
!------------------------------------------------------------------
subroutine integ(psi, norm)
      
use global_vars, only:NR, Nstates, dR, dp
 implicit none
 integer I, J
      
 real(dp) norm(Nstates)
 complex(dp) psi(NR,Nstates)
      
 norm = 0.d0
 
 do J = 1, Nstates
    norm(J)= sum(abs(psi(:,J))**2)   
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

use global_vars, only:dt, Nstates,kap, lam, dp, idp
use data_au, only:im
use blas_interfaces_module, only : zgemv, dgemv

implicit none

 integer:: i, J
 real(dp):: w, u, uv, d, mu(Nstates,Nstates), q, E
 integer(idp) Info, Lwork

 complex(dp):: tout, b, z, p, pT
 dimension:: u(Nstates,Nstates), d(Nstates), b(Nstates,Nstates), &
         & z(Nstates,Nstates), tout(Nstates,Nstates), &
         & p(Nstates,Nstates), uv(Nstates,Nstates), &
         & pT(Nstates,Nstates)
 real(dp) work(1000) 
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
 
uv = 0._dp
uv = u

Lwork=-1_idp
JOBZ='V'
 call dsyev(JOBZ,'U', int(Nstates,kind=idp), uv, int(Nstates,kind=idp), &
        & d, work, Lwork,info)
! call jacobi(u,Nstates,d)
 Lwork = min(1000_idp, int(work(1)))
!     Solve eigenproblem.
 call dsyev('V', 'U', int(Nstates,kind=idp), uv, int(Nstates,kind=idp), &
        & d, work, Lwork, info)

 if( info.gt.0 ) then
     write(*,*)'The algorithm failed to compute eigenvalues.'
     stop
 endif
 
 p = (0._dp,0._dp)
 p = uv ! transfering eigenvectors to a complex array
!print*,"eigen vector p"
!write( * , * ) ( (p(i,j),j=1,Nstates), i=1,Nstates )

 b= (0._dp,0._dp)
 do J = 1,Nstates  
   b(J,J) = exp(-im * dt * d(J)) ! exponentiate diagonal matrix i.e. eigenvalues
 end do
!print*, "b"
!write( * , * ) ((b(i,j),j=1,Nstates), i=1,Nstates )

! Calculating e^u = P (e^d) P^(-1)
!z = matmul(p,b)
 call zgemm('N', 'N', int(Nstates,kind=idp), int(Nstates,kind=idp), &
         & int(Nstates,kind=idp),(1._dp,0._dp), p, size(p,dim=1,kind=idp), &
         & b, size(b,dim=1,kind=idp), (0._dp,0._dp), z, size(z,dim=1,kind=idp))
!print*, "z"
!write( * , * ) ((z(i,j),j=1,Nstates), i=1,Nstates )
pT = transpose(p)
!tout = matmul(z,transpose(p))
 call zgemm('N', 'N', int(Nstates,kind=idp), int(Nstates,kind=idp), &
         & int(Nstates,kind=idp),(1._dp,0._dp), z, size(z,dim=1,kind=idp), &
         & pT, size(pT,dim=1,kind=idp), (0._dp,0._dp), tout, size(tout,dim=1,kind=idp))
!print*, "tout"
!write( * , * ) ((tout(i,j),j=1,Nstates), i=1,Nstates )

return
end subroutine

!!------------------------------------------------
subroutine mask_function_cos(cof)
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
    cof(j) = cos(((R(J) - Rend + cpmR) / -cpmR) * (0.5d0 * pi))
    cof(j) = cof(j)**2
    end if
   end do

 return
 end subroutine
!------------------------------------------------
 subroutine mask_function_ex(cof)
 use global_vars, only:NR, R, dR, dp
 use data_au
 use pot_param

  implicit none

  integer :: J
  real(dp):: cof(NR),c

  c=1.00d0
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




