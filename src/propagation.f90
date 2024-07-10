subroutine propagation_1D(chi0, E, A)

use global_vars
use pot_param
use data_au
use FFTW3
use omp_lib

 implicit none   
! include "/usr/include/fftw3.f" 
 
 integer I, J, K,II
 integer L, M, N, void
 integer*8 planF, planB
 integer eR
 real*4 dummy, dummy2, dummy3
 character(150):: filename, filename1, filename2
 character(150):: filename3, filename4, filename5
 character(150):: filename6, filename7, filename8
 character(150):: filename9

 double precision dt2, time
 double precision c, sp, deR 
 double precision:: normpn(Nstates), spec(NR,Nstates)
 double precision E(Nt), E1, norm(Nstates), E21, E22, norm_out(Nstates),norm1,norm_diss
 double precision A(Nt)
 double precision evR(Nstates), epr(Nstates),momt(Nstates), tot_momt
 double precision :: H_ac(NR,Nstates), Energy_axis(NR)
 double precision ::chi0(nr,vstates), Boltzmann_populations(Vstates)
 double precision :: norm_overlap(Nstates),norm_outR(Nstates), norm_outP(Nstates)
 double precision :: norm_gesP(Nstates), norm_gP_over(Nstates)
 double precision :: vib_pop(Vstates) 
 double precision, allocatable, dimension(:):: cof
 complex*16:: tout(Nstates,Nstates)
 complex*16, allocatable, dimension(:,:):: psi_ges, psi_out
 complex*16, allocatable, dimension(:):: psi_diss, psi_ex_diss
 complex*16, allocatable, dimension(:,:):: psi_loc, psi_ges1
 complex*16, allocatable, dimension(:):: psi, kprop, kprop1
 complex*16, allocatable, dimension(:,:):: psi_out1, psi_outR
 complex*16, allocatable, dimension(:,:):: psi_gesP
 complex*16, allocatable, dimension(:):: psi_chi

 
 open(100,file='psi0_1d.out',status='unknown') 
 open(101,file='cof_1d.out',status='unknown')
 open(102,file='effective_dipole.out',status='unknown')
 open(200,file='dens_1d.out',status='unknown')
 open(201,file='ex_dens_1d.out',status='unknown')
 open(800,file='R_1d.out',status='unknown')
 open(801,file='PR_1d.out',status='unknown')
 open(906,file='norm_pn_1d.out',status='unknown')
 open(909,file='field_1d.out',status='replace')
 open(999,file='accumulation.out',status='unknown')
 open(998,file='momentum_1d.out',status='unknown')
! open(1000,file='norm.out')
! open(1111,file='vib_pop.out',status='unknown')
! open(1112,file='KER_spectra.out',status='unknown')
!open(1001,file='dipcurves_3d.dat')
! open(1002,file='Dissociation_yield.out',status='unknown')
! open(1003,file='Ex_dissociation_yield.out',status='unknown')
 
  write(filename,fmt='(a,i0,a)') 'norm1D_lambda',int(lambda1),'nm.out'
  open(908,file=filename,status='unknown')
  write(filename1,fmt='(a,i0,a)') 'vibpop1D_lambda', int(lambda1),'nm.out'
  open(1111,file=filename1,status='unknown')
  write(filename2,fmt='(a)') 'Momentum_spectra_1g_1D.out'
  open(1112,file=filename2,status='unknown')
  write(filename3,fmt='(a)') 'Momentum_spctra_2u_1D.out'
  open(1113,file=filename3,status='unknown')
  write(filename2,fmt='(a)') 'Momentum_spectra_total_1D.out'
  open(1122,file=filename2,status='unknown')
  write(filename2,fmt='(a)') 'KER_spctra_1g_1D.out'
  open(1114,file=filename2,status='unknown')
  write(filename3,fmt='(a)') 'KER_spctra_2u_1D.out'
  open(1115,file=filename3,status='unknown')
  write(filename2,fmt='(a)') 'KER_spctra_total_1D.out'
  open(1144,file=filename2,status='unknown')
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
 allocate(psi_out1(nr,Nstates), psi_outR(nr,Nstates),psi_gesP(nr,Nstates))
 allocate(psi_chi(Vstates),psi_diss(NR),psi_ex_diss(NR))

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
     
 psi_ges = (0.d0,0.d0)    
 
 do I = 1, NR
!   psi_ges(I,1) = (0.55d0*chi0(I,1)+0.23d0*chi0(I,2)+0.11d0*chi0(I,3) &
!                 & + 0.07d0 * chi0(I,4) + 0.04d0*chi0(I,5)+ 2.0450413213653231E-002*chi0(I,6) &
!                 & + 1.4033977605843639E-002 *chi0(I,7)+ 1.0613089628800979E-002*chi0(I,8) &
!                 & + 8.8496818475779886E-003*chi0(I,9)+ 8.0611229210941892E-003*chi0(I,10))   !Thermal distribution over vibrational levels
   psi_ges(I,1) = chi0(I, 1) !exp(kappa*(R(I)-RI)**2) !sqrt(0.55d0)*chi0(I,1)+sqrt(0.23d0)*chi0(I,2)+sqrt(0.11d0)*chi0(I,3) &
                 !& +sqrt( 0.07d0 )* chi0(I,4) + sqrt(0.04d0)*chi0(I,5) ! + sqrt(2.0450413213653231E-002)*chi0(I,6) &
!                 & +sqrt( 1.4033977605843639E-002) *chi0(I,7)+ sqrt(1.0613089628800979E-002)*chi0(I,8) &
!                 & +sqrt( 8.8496818475779886E-003)*chi0(I,9)+ sqrt(8.0611229210941892E-003)*chi0(I,10)   !Thermal distribution over vibrational levels

!   psi_ges(I,1)=chi0(i,1)!exp(kappa*(R(I)-RI)**2)
!   psi_ges(I,1) = sum(chi0(I,:)*sqrt(Boltzmann_populations(:)))
   kprop(I) = exp(-im *dt * PR(I)**2/(4.d0*m_red))  ! pR**2 /2 * red_mass UND Half time step
   kprop1(I) =exp(-im *dt * PR(I)**2/(2.d0*m_red))
 end do 
 
 ! cpm = 3.d0/ au2a
 ! call cutoff_cos(cpm,cof)
  call cutoff_ex(cof)

 
 call integ(psi_ges, norm) 
 print*,'norm =', sngl(norm) 
  
 psi_ges(:,1) = psi_ges(:,1) / sqrt(norm(1))
 
 call integ(psi_ges, norm)
 print*,'norm =', sngl(norm) 
 

 do I = 1, NR      
  write(100,*) sngl(R(I) *au2a),sngl(abs(psi_ges(I,1))**2)         
  write(101,*) sngl(r(I)*au2a), sngl(cof(i))
 end do
close(102)

           
!______________________________________________________________________
!                                                                     
!                   Propagation Loop                         
!______________________________________________________________________                                                                     

 print*
 print*,'1D propagation...'
 print*

!   do I=1,NR
!     read(1000,*) dummy, mu(i,1), mu(i,2), mu(i,3)
!   enddo
    
  psi_out=0.00d0
  do I=1,NR
    psi_out(I,:)=psi_ges(I,:)*(1.0d0-cof(I))
  enddo


  do J = 1,1
    do I = 1,NR
      psi(I) = psi_out(I,J)  ! Hilfsgroesse
    end do
      call dfftw_execute(planF)
      psi = psi /sqrt(dble(nr))
    do I = 1,NR
      psi_out(I,J) = psi(I)
    end do
  end do
 
! kap =1.d0
! call pulse(E,A)  
!E22=0.0d0
timeloop: do K = 1, Nt

    time = K * dt
    epr = 0.d0
    evr = 0.d0
    psi=0.0d0
    momt=0.0d0
    norm_out=0.00d0


! if(K .lt. 1720) then
! if(K .lt. 1800) then  ! k = 1800, t = 9 fs
!    write(800,*) sngl(time*au2fs), RI
!    cycle
!endif
 

    
    do J = 1,Nstates
      do I = 1,NR
         psi(I) = psi_ges(I,J)  ! Hilfsgroesse
      end do 
       call dfftw_execute(planF)
       psi = psi * kprop
       call dfftw_execute(planB)
       psi = psi / dble(NR)
 
      do I = 1,NR
        psi_ges(I,J) = psi(I)
      end do      
    end do
 
 

    
   do i = 1, NR
!    mu_all(i,1,2)=0.8*R(I)
!    mu_all(i,2,1)=0.8*R(I)
    call pulse2(tout, mu_all(:,:,I), E(K)) 
!    tout(1,1) = cos(-kap*mu_all(I,1,2)*E(K)*dt)
!    tout(1,2) = -im*sin(-kap*mu_all(I,1,2)*E(K)*dt)
!    tout(2,1) = tout(1,2)
!    tout(2,2) = tout(1,1)

!     print*, "tout 2"
!     write( * , * ) ((tout(l,j),j=1,Nstates), l=1,Nstates)
!    stop
!    psi_ges(i,1:2) = matmul(tout(1:2,1:2),psi_ges(i,1:2))  
    psi_ges(i,1:Nstates) = matmul(tout(1:Nstates,1:Nstates),psi_ges(i,1:Nstates))  
   end do
    

   do j = 1, Nstates   
    do i = 1, NR  
!    psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * (adb(i,j)+kap*(mu_all(I,J,J) &
!            & +(mr-0.5)*R(I))*E(K) +lam*R(I)*E(K)))!+H_ac(i,j))) !         
!    psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * (adb(i,j)+kap*(mu_all(I,J,J) &
!            & +0.5d0*R(I))*E(K)+((2-kap)*mr+kap-1)*R(I)*E(K)))!+H_ac(i,j))) !         
    psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * (adb(i,j))) !+0.8d0*R(I)*E(K)))!+H_ac(i,j))) !         
    end do
   end do
  
 
  
    do J = 1,Nstates
      do I = 1,NR
         psi(I) = psi_ges(I,J)  ! Hilfsgroesse
      end do 
      call dfftw_execute(planF) 
      psi = psi * kprop 
      psi = psi /sqrt(dble(nr))     
      do i = 1, nr
        epr(j) = epr(j) + abs(psi(i))**2 * pr(i)
      end do
      epr(j) = epr(j) * dr
      call dfftw_execute(planB)
      psi = psi / sqrt(dble(NR))

      do I = 1,NR
         psi_ges(I,J) = psi(I)
      end do      
    end do     
 

!-----------------------------------------

  
   do I = 1, NR   
    evR(1) = evR(1) + abs(psi_ges(I,1))**2 * R(I)
    evR(2) = evR(2) + abs(psi_ges(I,2))**2 * R(I)
    psi_loc(i,1) = 1.d0/sqrt(2.d0)*(psi_ges(i,1) + psi_ges(i,2))
    psi_loc(i,2) = 1.d0/sqrt(2.d0)*(psi_ges(i,1) - psi_ges(i,2))
   end do   
    evR = evR * dR




    call integ(psi_ges, norm)    
    call integ(psi_loc, normPn)
    

    do j = 1,Nstates
       if (norm(j).ge.1.d-8) then
        evR(j) = evR(j)/norm(j)
        epr(j) = epr(j)/norm(j)
       end if
    end do
    
! ------------ Popolation analysis in vibrational states ----------------------
   vib_pop = 0.d0
   do J=1,(vstates-1)
   psi_chi(J)=0.0d0
   do I=1, NR
    psi_chi(J)=psi_chi(J)+ chi0(I,J) * (psi_ges(I,1))
   enddo
    psi_chi(J)=psi_chi(J)*dR
    vib_pop(J)=abs(psi_chi(J)**2)
   enddo


    write(1111,*) sngl(time*au2fs), (vib_pop(1)), (vib_pop(2)), (vib_pop(3)), (vib_pop(4)),&
                 &  vib_pop(5), vib_pop(6), vib_pop(7), vib_pop(8), vib_pop(9),&
                 &  vib_pop(10), vib_pop(11), vib_pop(12), sum(vib_pop(:)) 
         
     
    write(800,*) sngl(time *au2fs), sngl(evR)
    write(801,*) sngl(time *au2fs), sngl(epr)
    write(908,*) sngl(time *au2fs), norm
    write(906,*) sngl(time *au2fs), sngl(normPn)
    write(909,*) sngl(time *au2fs), sngl(E(K))      

   if(mod(K,100).eq.0) then
    do I = 1, NR/2   
    write(200,*) sngl(time *au2fs), sngl(R(I) *au2a), sngl(abs(psi_ges(I,1)**2))
    write(201,*) sngl(time *au2fs), sngl(R(I) *au2a), sngl(abs(psi_ges(I,2)**2))
    end do 
    write(200,*)
    write(201,*)
  
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

   do i = 1, NR
    psi_outR(i,:) = psi_ges(i,:) * (1.0d0-cof(i))
    psi_ges(i,:) = psi_ges(i,:) * cof(i) ! psi_ges = psi_nondiss
   enddo

!   norm_overlap(1)=2.0d0*dot_product(psi_outR(:,1),psi_ges(:,1))*dr
!   norm_overlap(2)=2.0d0*dot_product(psi_outR(:,2),psi_ges(:,2))*dr
!
!   call integ(psi_outR,norm_outR)
!   norm1=norm1+norm_overlap(1)+norm_overlap(2)+norm_outR(1)+norm_outR(2)
!
!   do j=1,Nstates
!     do i=1,nr
!      psi(i)=psi_outR(i,j)
!     end do
! 
!      call dfftw_execute(planF)
!      psi=psi/sqrt(dble(nr))
!     do i=1,NR
!      psi_outR(i,j)=psi(i)
!     enddo
!   enddo
!
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
  
!KER spectra
  psi_diss=0.0d0
  do J=1,vstates
  psi_chi(J)=0.0d0
  do I=1, NR
   psi_chi(J)=psi_chi(J)+ chi0(I,J) * (psi_ges(I,1))
  enddo
   psi_chi(J)=psi_chi(J)*dR
  do I=1,NR
   psi_diss(I)=psi_diss(I)+psi_chi(J)*chi0(I,J)
  enddo
  norm_diss=sum(abs(psi_diss(:))**2)*dR
  print*, 'Diss_yield (',J,') =',norm_diss
  enddo
  

  psi_diss(:)=psi_ges(:,1)-psi_diss(:)
!   psi_diss(:)=psi_diss(:)+psi_ges(:,2)
  
  psi(:)=psi_diss(:)+ psi_ges(:,2)
  call dfftw_execute(planF)
  psi=psi/sqrt(dble(NR))
  norm_diss=sum(abs(psi(:))**2)*dR
  write(*,*) "Dissociation yield 1g =",norm_diss
  psi=psi/sqrt(norm_diss)
   do I=NR/2 +1,NR
     write(1122,*) pR(I),abs(psi(I))**2
   enddo
   do I=1,NR/2
     write(1122,*) pR(I),abs(psi(I))**2
     write(1144,*) (pR(I)**2/(2*m_red)),abs(psi(I))**2
   enddo

  
  
  psi=psi_diss
  call dfftw_execute(planF)
  psi=psi/sqrt(dble(NR))
  psi_diss=psi
  norm_diss=sum(abs(psi_diss(:))**2)*dR
  write(*,*) "Dissociation yield 1g =",norm_diss
  write(*,*) lambda1,norm_diss
  psi_diss=psi_diss/sqrt(norm_diss)
   do I=NR/2 +1,NR
     write(1112,*) pR(I),abs(psi_diss(I))**2
   enddo
   do I=1,NR/2
     write(1112,*) pR(I),abs(psi_diss(I))**2
     write(1114,*) (pR(I)**2/(2*m_red)),abs(psi_diss(I))**2
   enddo

  psi(:)=psi_ges(:,2)
  call dfftw_execute(planF)
  psi=psi/sqrt(dble(NR))
  psi_ges(:,2)=psi(:)
  norm_diss=sum(abs(psi_ges(:,2))**2)*dR
  write(*,*) "Dissociation yield 2u =",norm_diss
  write(*,*) lambda1, norm_diss
  psi_ges(:,2)=psi_ges(:,2)/sqrt(norm_diss)
   do I=NR/2 +1,NR
     write(1113,*) pR(I),abs(psi_ges(I,2))**2
   enddo
   do I=1,NR/2
     write(1113,*) pR(I),abs(psi_ges(I,2))**2
     write(1115,*) (pR(I)**2/(2*m_red)),abs(psi_ges(I,2))**2
   enddo

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
      
use global_vars, only:NR, Nstates, dR
 implicit none
 integer I, J
      
 double precision norm(Nstates)
 complex*16 psi(NR,Nstates)
      
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

use global_vars, only:dt, Nstates,kap, lam
use data_au, only:im

implicit none

 integer:: i, J
 double precision:: w, u, d, mu(Nstates,Nstates), q, E
 integer Info, Lwork

 complex*16:: tout, b, z
 dimension u(Nstates,Nstates), d(Nstates), b(Nstates,Nstates), z(Nstates,Nstates), tout(Nstates,Nstates)
 double precision work(1000), u1(Nstates, Nstates)
 character(len=1):: JOBZ

    
!u(1,1) = 0.d0
!u(1,2) = -kap*mu(1,2) * E

!u(2,1) = -kap*mu(2,1) * E
!u(2,2) = 0.d0

!Diapole matrix
u=0.0d0
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
 LWORK = MIN( 1000, INT( WORK( 1 ) ) )
!     Solve eigenproblem.

 CALL DSYEV( 'V', 'U', Nstates, u, Nstates, d, WORK, LWORK, INFO )

 IF( INFO.GT.0 ) THEN
     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
     STOP
 END IF
!print*,"eigen vector u"
!write( * , * ) ( (u(i,j),j=1,Nstates), i=1,Nstates )
b= (0.d0,0.d0)
 
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
use global_vars, only:NR, R, dR
use data_au
use pot_param

implicit none

   integer :: J
   double precision:: cof(NR)

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
 use global_vars, only:NR, R, dR
 use data_au
 use pot_param

  implicit none

    integer :: J
    double precision:: cof(NR),c

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




