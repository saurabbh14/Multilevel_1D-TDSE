!module FFTW_check
!
!use FFTW3
!
!implicit none
!
!interface
!function fftw_plan_r2r_1d(plan,NR,pin,pout,*)
!  import C_INT, C_DOUBLE
!  type(C_PTR):: plan
!  integer(C_INT):: NR
!  real(C_DOUBLE):: pin
!  real(C_DOUBLE):: pout
!end function 
!end interface
!end module


subroutine nuclear_wavefkt

 use global_vars
 use pot_param
 use data_au
 use FFTW3
 use omp_lib

implicit none
  
  integer:: I, J, K, M, V, G   
  integer:: istep,void  
  type(C_PTR):: planF, planB, p
  character(1000):: filename, filepath
  integer :: scount, scount1, scount2,  & ! Starting "time"
            ecount ! Ending "time"
  integer  rate, rate1, rate2       ! number of clock ticks per second
  real(dp) timer_elapsed_time
  logical ext
     
  real(dp):: dt2
  real(dp):: E1, norm
  real(dp), allocatable:: E(:,:)
  real(dp):: CONS, thresh
  real(dp):: dummy, R_e, D_e, alpha, Rin
!  double precision, allocatable:: chi0(:,:)
!  double precision:: total_pop, Boltzmann_populations(guess_Vstates) 
!  double precision:: trans_dipole(guess_Vstates,guess_Vstates) 
 
  real(dp), allocatable, dimension(:):: vprop
  real(dp), allocatable, dimension(:):: psi1
  !real(C_DOUBLE), allocatable, dimension(:):: psi
  real(C_DOUBLE), pointer:: psi(:)
  real(dp), allocatable, dimension(:,:):: ref

  integer trans_dipole_tk, vib_en_tk, vstates_tk
  integer chi0_vib_en_tk, chi0_tk


! allocate(psi(NR))
 allocate(psi1(NR))
 allocate(vprop(NR), E(guess_vstates,Nstates))
 allocate(ref(NR,guess_vstates))
 p = fftw_alloc_real(int(NR, C_SIZE_T)) 
 call c_f_pointer(p,psi,[NR])


 void=fftw_init_threads( )
 if (void==0) then
    write(*,*) 'Error in fftw_init_threads, quitting'
    stop
 else
    write(*,*) 'Thread initialization =', void
 endif

 select case(ITP_par_FFTW)
 case("parallel")
   call fftw_plan_with_nthreads(omp_get_max_threads())    
   planF = fftw_plan_r2r_1d(NR, psi, psi, FFTW_R2HC, FFTW_MEASURE)
   planB = fftw_plan_r2r_1d(NR, psi, psi, FFTW_HC2R, FFTW_MEASURE)

 case default
   planF = fftw_plan_r2r_1d(NR, psi, psi, FFTW_R2HC, FFTW_MEASURE)
   planB = fftw_plan_r2r_1d(NR, psi, psi, FFTW_HC2R, FFTW_MEASURE)
 end select

  dt2 = dt!*10
  thresh = 1d-18 
  istep = 1e8

  print*
  print*, "reduced mass:", m_red
!  print*, "mass ratio:", mr
  print*, "kap:", kap
  print*, "lam:", lam
  print*

 !Morse potential
!  R_e = 0.7743/au2a
!  D_e = 2.04/au2eV
!  alpha = 2.7407*au2a
!
!  do I =1, NR
!   adb(I,1) = D_e * (exp(-alpha*(R(I)-R_e))-1)**2
!   write(90,*) R(I), adb(I,1)*au2eV
!  enddo

!............... Main Propagation ........................                   

  print*
  print*, 'Start of energy calculation...'
  print*


 E = 0._dp
  write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), &
          & "Bound-vibstates_in_Nthstates.out"
 open(newunit=vstates_tk,file=filepath,status='unknown')
 call system_clock(scount,rate)
Nloop: do J = 1, Nstates ! varying the different adiabatic states
  write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), &
          & "BO_Electronic-state-g", int(J-1), "_vibstates.out"
  inquire(file=filepath, Exist=Ext)
  print*, trim(filepath), " ", Ext
  if (Ext) then 
    print*, "File already exists." ! Exiting the ITP if the vibrational wavepackets exist.
    cycle
  endif

  call system_clock(scount1,rate1)

  do i = 1, NR ! new vprop
     vprop(i) = exp(-0.5_dp * dt2 * adb(i,J))
  end do
  Rin = R(minloc(adb(:,J),1))
  
  print*, 'Dp =', adb(NR,J)*au2eV, 'eV'
  write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), &
         & "BO_Electronic-state-g", int(J-1), "_Evib.out"
  open(newunit=vib_en_tk,file=filepath,status='unknown')
  write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), &
         & "BO_Electronic-state-g", int(J-1), "_chi0-Evib.out"
  open(newunit=chi0_vib_en_tk,file=filepath,status='unknown')
Vloop: do V = 1, guess_Vstates ! loop over the vibrational state
 call system_clock(scount2,rate2)
  
      do i = 1, NR  ! symmetry for the startup function
        psi(i) = exp(kappa * (R(I) -Rin )**2) + &
        &   (-1._dp)**(V - 1) * exp(-0.5_dp * (R(I) + Rin)**2)
      end do
      
      call integ_r(psi, psi, norm)
      psi = psi / sqrt(norm)

!.......... Imaginary Time Propagation ........................

  do K = 1, istep

      psi1 = psi        ! storing wave function of iteration step (N - 1)
      E1 = E(V,J)            ! storing eigenvalue of iteration step (N - 1)
      if (V.gt.1) then   !projecting out the vibrational ground state...
           do G = 1, (V - 1)

            call integ_r(ref(1:NR,G), psi, norm)

            do i = 1, NR
             psi(i) = psi(i) - norm * ref(i, G)
            end do

           end do
      end if

      psi = psi * vprop
      call fftw_execute_r2r(planF, psi, psi)   
      psi = psi * exp((-dt2 * pR**2) / (2._dp*m_red))   
      call fftw_execute_r2r(planB, psi, psi)   
      psi = psi / dble(Nr)         
      psi = psi * vprop

      call eigenvalue_r(psi, psi1, E(V,J), dt2)
      call integ_r(psi, psi, norm)

      psi = psi / sqrt(norm)

      if(abs(E(V,J) - E1).le.thresh) then
        print*, 'Surface', J-1, 'Vibrational state', V-1, E(V,J) * au2eV, 'eV'
        print*, 'Resonance freq for dissociation', (adb(NR,J)-E(V,J))*au2eV, 'eV'
        print*, 'Resonance freq for dissociation', eV2nm/((adb(NR,J)-E(V,J))*au2eV), 'nm'
        if (V.gt.1) then
          do G=1,(V-1)
            print*, "freq for the transition from", G-1, "to", V-1, ":", &
                    & eV2nm/(abs(E(G,J)-E(V,J))*au2eV), "nm E=", &
                    & abs(E(G,J)-E(V,J))*au2eV, "eV"
          enddo
        endif
        write(vib_en_tk,*) v-1, e(V,J)*au2ev
        do I = 1, NR
         ref(I,V) = psi(I)             ! storing as reference for the next loop
         write(chi0_vib_en_tk,*) R(I), ref(I,V)+E(V,J)*au2eV     
        end do
        write(chi0_vib_en_tk,*)
        exit
      elseif(K .eq. istep) then
        print*,'Iteration not converged!'
        print*,'Program stopped!'
        print*
        print*,'Surface:', J-1
        print*,'Vibstate:', V-1
        print*,'E =', E(V,J) *au2eV
        print*,'E1 =', E1 *au2eV
        print*,'thresh =', thresh *au2eV
        print*,'step =', K
        write(filepath, '(a,a,i0,a,i0,a)') adjustl(trim(output_data_dir)), &
                & "Electronic-state-g", int(J-1), "_vib-state", int(V-1), &
                & "_uncoverged-final-wf.out"
        open(102,file=filepath,status='unknown')
        do i = 1, nr
         write(102,*) R(I), psi(i), abs(psi(I))**2
        end do
        close(102)
        stop
      else
        cycle
      end if
  end do
  if (E(V,J) .gt. adb(NR,J)) then
     Vstates(J) = V-1
     exit
  endif
  call system_clock(ecount)
  timer_elapsed_time = real(ecount-scount2,8)/real(rate2,8)
  write(*,*) "Calculated run time for vibstate", int(V-1),"on the surface", &
          & int(J-1), "is",timer_elapsed_time," seconds"
  print*,""
 end do Vloop               ! end of vibrational states loop
 close(vib_en_tk)
 close(chi0_vib_en_tk)
 print*, "Number of bound vibrational states in the ", int(J-1), &
         & " electronic state:", Vstates(J)
 write(vstates_tk,*) J-1, Vstates(J)
 do v = 1, vstates(J)
   do I=1,NR
     chi0(I, V, J) = ref(I, V)
   enddo
 enddo
 
 write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), &
        & "BO_Electronic-state-g", int(J-1), "_vibstates.out"
 open(newunit=chi0_tk,file=filepath,status='unknown')
 do I = 1, NR
   write(chi0_tk,*) R(I), chi0(I,1:vstates(J), J)
 enddo 
 close(chi0_tk)

 call system_clock(ecount)
 timer_elapsed_time = real(ecount-scount1,8)/real(rate1,8)
 write(*,*) "Calculated run time for the surface", int(J-1), &
         & "is",timer_elapsed_time," seconds"
 print*
end do Nloop            ! end of surface loop
close(vstates_tk)
 call system_clock(ecount)
 timer_elapsed_time = real(ecount-scount,8)/real(rate,8)
 write(*,*) "Calculated ITP run time is ",timer_elapsed_time," seconds"

 call fftw_destroy_plan(planF)
 call fftw_destroy_plan(planB)
 call fftw_free(p) 

 deallocate(psi1, vprop, ref)

return
end

!_________________ Subroutines______________________________________


subroutine eigenvalue_R(A, B, E, dt2)      
      
use global_vars, only: NR, dp
 implicit none  
 real(dp):: E, e1, e2, norm
 real(dp), intent(in):: dt2, A(nr), B(nr)
 
            
  call integ_r(B, B, norm)  
  e1 = norm
  
  call integ_r(A, A, norm)  
  e2 = norm
  
  
  E = (-0.5d0/dt2) * log(e2/e1)
 

return
end subroutine    
  
! ........................................................
                                                          
subroutine integ_r(A, B, C)

use global_vars, only:NR, dR, dp  
 implicit none  
 integer I 
 real(dp),intent(in):: A(Nr), B(Nr)
 real(dp):: C
  
  C = 0.d0
  
  do I = 1, Nr  
   C = C + A(I) * B(I)   
  end do
  
  C = C * dr
  
return  
end subroutine

!!! COMMENTED Boltzmann Distribution part
! total_pop = sum(exp(-E(:)/(kB*temperature)))
! print*, "total populations =", total_pop

! do V = 1, Vstates
!   if (V .gt. 1) then
!   Boltzmann_populations(V) = exp(-(E(V)-E(1))/(kB*temperature))!/total_pop
!   else
!   Boltzmann_populations(V)=1.d0
!   endif
!   print*, "state", V-1, "population ratio =", Boltzmann_populations(V)
! enddo
!
! total_pop = sum(Boltzmann_populations(:))
! print*, "sum of ratios =", total_pop
! Boltzmann_populations = Boltzmann_populations/total_pop
! total_pop = sum(Boltzmann_populations(:))
! print*, "total population =", total_pop
!
! do V =1, Vstates
!   print*, "state", V-1, "population probability =", Boltzmann_populations(V)
! enddo
! 
! print*
! do V=1, Vstates
!   write(filename,'(a,i0,a)') "trans_dipole_v_",V,".out"
!   open(unit=trans_dipole_tk, file=filename, status="unknown")
!   !if (V .gt. 1) then
!    do G = V,Vstates
!    trans_dipole(G,V) = abs(sum(ref(:,G)*(kap*(mu_all(:,1,1)+(mr-0.5d0)*R(:))+lam*R(:))*ref(:,V))*dR) !**2
!  !   trans_dipole(G,V) = abs(sum(ref(:,G)*(kap*mu_all(:,1,1)+lam*R(:))*ref(:,V))*dR) !**2
!     print*, 'transition dipole', V, "to", G, ":", trans_dipole(G,V)
!     write(trans_dipole_tk,*) G, trans_dipole(G,V) 
!    enddo
!   !endif
!   print*
! enddo
