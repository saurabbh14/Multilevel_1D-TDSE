subroutine nuclear_wavefkt(adb,chi0)

 use data_grid
 use pot_param
 use data_au
 use FFTW3
 use omp_lib

implicit none
!  include "/usr/include/fftw3.f"

  integer:: I, J, K, M, V, G   
  integer:: istep,void  
  integer*8 planF, planB
  character(150):: filename
     
  double precision:: dt2
  double precision:: E(vstates), E1, norm
  double precision:: CONS, thresh
  double precision:: dummy, R_e, D_e, alpha, Rin
  double precision:: adb(NR, Nstates),chi0(nr,vstates)
  double precision, parameter:: temperature=3100 
  double precision:: total_pop, Boltzmann_populations(Vstates) 
  double precision:: trans_dipole(Vstates,Vstates) 
  double precision:: mu_all(NR, Nstates, Nstates) 
 
  double precision, allocatable, dimension(:):: vprop
  double precision, allocatable, dimension(:):: psi, psi1
  double precision, allocatable, dimension(:,:):: ref

  integer trans_dipole_tk

 open(99,file='chi0.out',status='unknown')
 open(98,file='Evib.out',status='unknown')
 open(97,file='chi_indivisiual.out',status='unknown')
! open(90,file='Morse_potential.out',status='unknown')

  allocate(psi(NR), psi1(NR))
  allocate(vprop(NR), ref(NR, Vstates))

 void=fftw_init_threads( )
 if (void==0) then
    write(*,*) 'Error in fftw_init_threads, quitting'
    stop
 endif

 call dfftw_plan_r2r_1d(planF, NR, psi, psi, FFTW_R2HC, FFTW_ESTIMATE)
 call dfftw_plan_r2r_1d(planB, NR, psi, psi, FFTW_HC2R, FFTW_ESTIMATE)


  dt2 = dt!*10
  thresh = 1.d-18 
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



Nloop: do J = 1,1! Nstates ! varying the different adiabatic states


    do i = 1, NR ! new vprop
      vprop(i) = dexp(-0.5d0 * dt2 * adb(i,J))
    end do
    Rin = R(minloc(adb(:,J),1))

    print*, 'Dp =', adb(NR,J)*au2eV, 'eV'
    
    E = 0.d0

Vloop:   do V = 1, Vstates ! loop over the vibrational state


      do i = 1, NR  ! symmetry for the startup function
        psi(i) = exp(kappa * (R(I) -Rin )**2) +&
     &   (-1.d0)**(V - 1) * exp(-0.5d0 * (R(I) + Rin)**2)
      end do
      
      call integ_r(psi, psi, norm)
      psi = psi / sqrt(norm)
      


!.......... Imaginary Time Propagation ........................



    do K = 1, istep

      psi1 = psi        ! storing wave function of iteration step (N - 1)
      E1 = E(V)            ! storing eigenvalue of iteration step (N - 1)


      if (V.gt.1) then   !projecting out the vibrational ground state...
           do G = 1, (V - 1)

            call integ_r(ref(1:NR,G), psi, norm)

            do i = 1, NR
             psi(i) = psi(i) - norm * ref(i, G)
            end do

           end do
      end if


      psi = psi * vprop
      call dfftw_execute(planF)   
      psi = psi * exp((-dt2 * pR**2) / (2.d0*m_red))   
      call dfftw_execute(planB)   
      psi = psi / dble(Nr)         
      psi = psi * vprop


      call eigenvalue_r(psi, psi1, E(V), dt2)
      call integ_r(psi, psi, norm)

      psi = psi / sqrt(norm)

      if(abs(E(V) - E1).le.thresh) then
        print*, 'Surface', J, 'Vibrational state', V, E(V) * au2eV, 'eV'
        print*, 'Resonance freq for dissociation', eV2nm/((adb(NR,J)-E(V))*au2eV), 'nm'
        if (V.gt.1) then
          do G=1,(V-1)
            print*, "freq for the transition from", G, "to", V, ":", &
                    & eV2nm/(abs(E(G)-E(V))*au2eV), "nm E=", abs(E(G)-E(V))*au2eV, "eV"
          enddo
        endif
        print*,""
        
        write(98,*) v, e(V)*au2ev
      
        do I = 1, NR
         ref(I,V) = psi(I)             ! storing as reference for the next loop
         write(99,*) R(I)*au2a, ref(I,v)+E(V)*au2eV     
        end do
        write(99,*)
        exit
      elseif(K .eq. istep) then
        print*,'Iteration not converged!'
        print*,'Program stopped!'
        print*
        print*,'E =', E(V) *au2eV
        print*,'E1 =', E1 *au2eV
        print*,'thresh =', thresh *au2eV
        print*,'step =', K
        do i = 1, nr
         write(102,*) i, psi(i)
        end do
        stop
      else
        cycle
      end if
  end do

        
 end do Vloop               ! end of vibrational states loop

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
         


end do Nloop            ! end of surface loop

 do I=1,NR
  write(97,*) R(I)*au2a, ref(I,:)
 enddo 

  chi0 = ref

  call dfftw_destroy_plan(planF)
  call dfftw_destroy_plan(planB)

  close(99,status='keep')
  close(98,status='keep')

 deallocate(psi, psi1, vprop, ref)

return
end

!_________________ Subroutines______________________________________


subroutine eigenvalue_R(A, B, E, dt2)      
      
use data_grid
 implicit none  
 double precision:: E, e1, e2, norm
 double precision, intent(in):: dt2, A(nr), B(nr)
 
            
  call integ_r(B, B, norm)  
  e1 = norm
  
  call integ_r(A, A, norm)  
  e2 = norm
  
  
  E = (-0.5d0/dt2) * log(e2/e1)
 

return
end subroutine    
  
! ........................................................
                                                          
subroutine integ_r(A, B, C)

use data_grid  
 implicit none  
 integer I 
 double precision,intent(in):: A(Nr), B(Nr)
 double precision C
  
  C = 0.d0
  
  do I = 1, Nr  
   C = C + A(I) * B(I)   
  end do
  
  C = C * dr
  
return  
end subroutine

