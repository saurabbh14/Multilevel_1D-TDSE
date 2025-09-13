module propagation_mod
    use varprecision, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: time_prop !, propagation_1D

    type :: time_prop
        real(dp), allocatable :: chi0(:,:,:) 
        complex(dp), allocatable :: psi_ges(:,:)
        real(dp), allocatable :: vib_en(:,:)
        integer :: i_cpmR
        complex(dp), allocatable :: abs_func(:)

        complex(dp), allocatable :: psi_outR(:,:)
        real(dp), allocatable :: psi_outR_inc(:,:)
        
        integer:: psi_1d_tk, cof_1d_tk
        integer:: norm_1d_tk, dens_1d_tk, ex_dens_1d_tk, Pdens_1d_tk
        integer:: avgR_1d_tk, avgpR_1d_tk, norm_pn_1d_tk
        integer:: field_1d_tk, accumulation_1d_tk, momt_1d_tk
        integer:: vibpop_1d_tk, psi_outR_norm_1d_tk, psi_outR_Pdens_1d_tk
    contains
        procedure :: initialize
        procedure :: read_pot_files
        procedure :: ini_dist_choice
        procedure :: time_evolution
        procedure :: localized_states_norm, expected_position
        procedure :: vib_pop_analysis
        procedure :: propagation_1D
        procedure :: absorber_gen
        procedure :: open_files_to_write
        procedure :: continuum_prop
        procedure :: post_propagation_analysis => post_prop_analysis
    end type time_prop

    type :: split_operator_type
        complex(dp), allocatable :: kprop_half(:), kprop_full(:)
        complex(dp), allocatable :: vprop(:,:)
        type(C_PTR) :: planF, planB, p
        complex(C_DOUBLE), pointer:: psi(:)
    contains
        procedure :: fft_initialize
        procedure :: kprop_gen, vprop_gen
        procedure :: split_operator
    end type split_operator_type

contains

    subroutine propagation_1D(this, E, A)
        use global_vars, only: Nt
        class(time_prop), intent(inout) :: this
        real(dp), intent(in) :: E(Nt), A(Nt)
        call this%initialize()
        call this%read_pot_files()
        call this%ini_dist_choice()
        call this%absorber_gen()
        call this%open_files_to_write()
        call this%time_evolution(E, A)


    end subroutine propagation_1D

    subroutine initialize(this)
        use global_vars, only: NR, Nstates, guess_vstates
        class(time_prop), intent(inout) :: this
        allocate(this%chi0(NR,guess_vstates,Nstates))
        allocate(this%psi_ges(NR,Nstates))
        allocate(this%psi_outR(NR,Nstates), this%psi_outR_inc(NR,Nstates))
        allocate(this%vib_en(guess_vstates,Nstates))
    end subroutine initialize

    subroutine read_pot_files(this)
        use global_vars, only: NR, Nstates, Vstates, output_data_dir
        use varprecision, only: dp, sp
        class(time_prop), intent(inout) :: this
        character(150):: filepath
        integer:: chi0_tk, vstates_tk, vib_en_tk
        integer:: i, N, V, i_dummy
        real(sp):: dummy

        ! implementation for reading potential files
        write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "Bound-vibstates_in_Nthstates.out"
        open(newunit=vstates_tk,file=filepath,status='unknown')

        this%chi0 = 0._dp
        this%vib_en = 0._dp
        do N = 1, Nstates 
            read(vstates_tk,*) i_dummy, Vstates(N)

            write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), &
                & "BO_Electronic-state-g", int(N-1), "_Evib.out"
            open(newunit=vib_en_tk,file=filepath,status='unknown')
            do V = 1, Vstates(N)
                read(vib_en_tk,*) i_dummy, this%vib_en(V,N)
            enddo
            close(vib_en_tk)

            write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "BO_Electronic-state-g", &
                & int(N-1), "_vibstates.out"
            open(newunit=chi0_tk,file=filepath,status='unknown')
            print*, "NR:", NR, "Vstates(N)", Vstates(N)

            do i = 1, NR
                read(chi0_tk,*) dummy, this%chi0(i,1:Vstates(N),N)
                print*, i, dummy, this%chi0(i,Vstates(N),N)
            enddo 
            close(chi0_tk)
        enddo
        close(vstates_tk)
    end subroutine read_pot_files

    subroutine ini_dist_choice(this)
        use global_vars, only: NR, Nstates, Vstates, v_ini, N_ini, Ri_tdse, kappa_tdse, &
            initial_distribution, R
        use data_au, only: au2a, au2eV
        class(time_prop), intent(inout) :: this
        integer :: i, N, v
        real(dp):: norm(Nstates)
        real(dp), allocatable, dimension(:) ::vib_dist
        select case(initial_distribution) 
            case("single vibrational state")
                print*, "initial wavefunction in..."
                print*, N_ini-1, "electronic state and in", v_ini-1, "vibrational state"
                do i = 1, NR
                    this%psi_ges(i,N_ini) = this%chi0(i,v_ini,N_ini) 
                enddo  

            case("gaussian distribution")
                print*, "initial wavefunction in..."
                print*, N_ini-1, "electronic state and with a Gaussian distribution centered around",&
                    & Ri_tdse/au2a, "a.u. \n with deviation of", kappa_tdse, "."
                do i = 1, NR
                    this%psi_ges(i,1)=exp(kappa_tdse*(R(i)-Ri_tdse/au2a)**2) 
                enddo

            ! case ("input dist")
                !  allocate(v_dist_ini(N_ini))
                !  v_ini_check = 0
                !  do N = 1, N_ini
                !    write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "vib_dist_", int(N-1), ".out"
                !    open(newunit=vib_dist_tk,file=filepath,status='unknown')
                !    do 
                !       read(vib_dist_tk,*,iostat=io)
                !       if (io /= 0) exit
                !       v_ini_check = v_ini_check + 1
                !    enddo
                !    if (v_ini /= v_ini_check) then
                !       write(*,'(a,a,i0,a)') "Number of vibstates in 'vib_dist_", N, ".out' not equal to input" 
                !    endif
                !  enddo
                !  allocate(vib_ini(guess_vstates,N_ini), vib_dist(guess_vstates, N_ini))
                !  do N = 1, N_ini
                !     do v = 1, v_dist_ini(N)
                !         read(vib_dist_tk,*) vib_ini(v,N), vib_dist(v,N) 
                !     enddo
                !  enddo
                !  do i = 1, NR
                !    do v = 1, v_ini
                !       do N = 1, N_ini
                !         psi_ges(i,N) = psi_ges(i,N) + vib_dist(vib_ini(v,N),N) * chi0(i,vib_ini(v,N),N)
                !       enddo
                !    enddo
                !  enddo

            case ("Boltzmann dist")
                do N = 1, N_ini
                    call Boltzmann_distribution(N, this%vib_en(:,N)/au2eV, vib_dist)
                    do v = 1, Vstates(N)
                        do i = 1, NR
                            this%psi_ges(i,N) = this%psi_ges(i,N) + vib_dist(v) * this%chi0(i,v,N)
                        enddo
                    enddo
                enddo
                deallocate(vib_dist)
              
            case default
                do i = 1, NR
                    this%psi_ges(i,1) = this%chi0(i,1,1)
                enddo
              
        end select

        ! Normalization of the initial wavefunction
        call integ(this%psi_ges, norm)
        this%psi_ges(:,1) = this%psi_ges(:,1) / sqrt(norm(1))
        this%psi_ges(:,2) = this%psi_ges(:,2) / sqrt(norm(2))
        
        ! Writing initial wavefunction to file
        do i = 1, NR
            write(this%psi_1d_tk,*) R(i), abs(this%psi_ges(i,1))**2
        end do
        close(this%psi_1d_tk)
    end subroutine ini_dist_choice

    subroutine absorber_gen(this)
        use global_vars, only: NR, R, dr, absorber
        use pot_param, only: cpmR
        use varprecision, only: dp
        class(time_prop), intent(inout) :: this
        real(dp), allocatable:: cof(:), V_abs(:)
        complex(dp), allocatable:: exp_abs(:)
        integer:: i
  
        allocate(this%abs_func(NR))

        this%i_cpmR = minloc(abs(R(:) - cpmR), 1) - 50
        print*, "i_cpmR =", this%i_cpmR, ", NR-i_cpmR", NR - this%i_cpmR
        print*, "R(NR-i_cpmR) =", R(NR - this%i_cpmR)

        select case(absorber)
            case ("CAP")
                allocate(V_abs(NR), exp_abs(NR))
                call Complex_absorber_function(V_abs, exp_abs)
                this%abs_func = exp_abs
            case("mask")
                allocate(cof(NR))
                call mask_function_ex(cof)
                this%abs_func = cof
            case default
                print*, "No absorber selected. Reflections off the grid boundary may occur!"
        end select

        ! Writing absorber function to file
        do i = 1, NR
            write(this%cof_1d_tk,*) R(i), abs(this%abs_func(i))
        end do
        close(this%cof_1d_tk)

    end subroutine absorber_gen

    subroutine open_files_to_write(this)
        use global_vars, only: output_data_dir
        class(time_prop), intent(inout) :: this
        character(150):: filepath
        ! inintial wavefunction
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "psi0_1d.out"
        open(newunit=this%psi_1d_tk,file=filepath,status='unknown') 
        ! Absorber function
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "absorber_function.out"
        open(newunit=this%cof_1d_tk,file=filepath,status='unknown') 
        
        ! Propagation outputs
        ! time dependent norm
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "norm_1d.out"
        open(newunit=this%norm_1d_tk,file=filepath,status='unknown')
        ! time dependent ground state density 
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "density_1d.out"
        open(newunit=this%dens_1d_tk,file=filepath,status='unknown')
        ! time dependent excited state density
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "ex_density_1d.out"
        open(newunit=this%ex_dens_1d_tk,file=filepath,status='unknown')
        ! time dependent momentum density 
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "Pdensity_1d.out"
        open(newunit=this%Pdens_1d_tk,file=filepath,status='unknown')
        ! time dependent average position 
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "avgR_1d.out"
        open(newunit=this%avgR_1d_tk,file=filepath,status='unknown') 
        ! time dependent average momentum
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "avgPR_1d.out"
        open(newunit=this%avgpR_1d_tk,file=filepath,status='unknown') 
        ! time dependent norm in localized states
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "norm_pn_1d.out"
        open(newunit=this%norm_pn_1d_tk,file=filepath,status='unknown')
        ! time dependent electric field 
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "field_1d.out"
        open(newunit=this%field_1d_tk,file=filepath,status='unknown')
        ! time dependent accumulated absorbed density 
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "accumulation_1d.out"
        open(newunit=this%accumulation_1d_tk,file=filepath,status='unknown') 
        ! time dependent absorbed momentum
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "momentum_1d.out"
        open(newunit=this%momt_1d_tk,file=filepath,status='unknown')
        ! time dependent vibrational populations 
        write(filepath,'(a,a)') adjustl(trim(output_data_dir)), 'vibpop1D_lambda.out'
        open(newunit=this%vibpop_1d_tk,file=filepath,status='unknown')
        ! time dependent norm of absorbed wavepacket
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "psi_outR_norm_1d.out"
        open(newunit=this%psi_outR_norm_1d_tk,file=filepath,status='unknown')
        ! time dependent momentum density of absorbed wavepacket 
        write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "psi_outR_Pdensity_1d.out"
        open(newunit=this%psi_outR_Pdens_1d_tk,file=filepath,status='unknown') 

    end subroutine open_files_to_write

    subroutine fft_initialize(this)
        use global_vars, only: NR, prop_par_FFTW
        use data_au, only: im
        use FFTW3
        class(split_operator_type), intent(inout) :: this

        ! Creating aligned memory for FFTW
        this%p = fftw_alloc_real(int(NR, C_SiZE_T)) 
        call c_f_pointer(this%p,this%psi,[NR])

        call fftw_initialize_threads
        call fftw_create_c2c_plans(this%psi, NR, this%planF, this%planB, prop_par_FFTW)

    end subroutine fft_initialize

    subroutine kprop_gen(this)
        use global_vars, only: NR, dpR, dt, m_red, PR
        use data_au, only: im
        class(split_operator_type), intent(inout) :: this

        allocate(this%kprop_half(NR), this%kprop_full(NR))
        
        this%kprop_half = exp(-im *dt * pR*pR /(4._dp*m_red))  ! pR**2 /2 * red_mass UND Half time step
        this%kprop_full = exp(-im *dt * pR*pR /(2._dp*m_red))  
         
    end subroutine kprop_gen

    subroutine vprop_gen(this)
        use global_vars, only: NR, Nstates, dt, adb, dp
        use data_au, only: im
        class(split_operator_type), intent(inout) :: this
        integer:: j
        
        allocate(this%vprop(NR,Nstates))
        do j = 1, Nstates           
            this%vprop(:,j) = exp(-im * 0.5_dp * dt * (adb(:,j))) !+0.8d0*R(i)*E(K)))!+H_ac(i,j))) !         
        end do

    end subroutine vprop_gen

    subroutine split_operator(this, psi_ges)
        use global_vars, only: NR, Nstates, m_red, pR
        use data_au, only: im
        use FFTW3
        class(split_operator_type), intent(inout) :: this
        complex(dp), intent(inout):: psi_ges(NR, Nstates)
        integer:: j

        do j = 1,Nstates
            this%psi = (0._dp, 0._dp)
            this%psi(:) = psi_ges(:,J)  ! Hilfsgroesse
            call fftw_execute_dft(this%planF, this%psi, this%psi)
            this%psi = this%psi * this%kprop_half
            call fftw_execute_dft(this%planB, this%psi, this%psi)
            this%psi = this%psi / dble(NR)
            psi_ges(:,J) = this%psi(:)      
        end do
    end subroutine split_operator

    subroutine time_evolution(this, E, A)
        use global_vars, only: NR, Nstates, time, mu_all, Nt, adb, &
            & dp, dR, guess_vstates, dt, Vstates, R, pR
        use data_au, only: im, au2fs
        use blas_interfaces_module, only: zgemv
        class(time_prop), intent(inout) :: this
        type(split_operator_type) :: split_operator
        integer :: i, j, k
        real(dp), allocatable, dimension(:) :: evr, momt
        real(dp), allocatable, dimension(:) :: norm, normPn, norm_outR
        real(dp), allocatable, dimension(:) :: norm_SE_outR
        real(dp) :: E(Nt), A(Nt) 
        complex(dp), allocatable :: tout(:,:)
        complex(dp), allocatable :: psi_Nstates(:), psi_Nstates1(:)
        real(dp), allocatable :: psi_Nstates_real(:), psi_Nstates_imag(:)
        real(dp), allocatable :: vib_pop(:,:)

        allocate(psi_Nstates(Nstates), psi_Nstates1(Nstates))
        allocate(psi_Nstates_real(Nstates), psi_Nstates_imag(Nstates))
        allocate(evr(Nstates), momt(Nstates))
        allocate(norm(Nstates), normPn(Nstates), norm_outR(Nstates))
        allocate(norm_SE_outR(Nstates))
        allocate(tout(Nstates,Nstates))
        allocate(vib_pop(guess_vstates,Nstates))

        call split_operator%fft_initialize()
        call split_operator%kprop_gen()
        call split_operator%vprop_gen()

        print*
        print*,'1D propagation...'
        print*
    
        this%psi_outR = (0._dp,0._dp)
        this%psi_outR_inc = 0._dp
        
        timeloop: do k = 1, Nt

            evr = 0._dp
            momt=0._dp
          !   norm_out=0._dp
            norm = 0._dp
            normPn = 0._dp
            norm_outR = 0._dp

            ! Kinetic propagation half step
            call split_operator%split_operator(this%psi_ges)
            ! Potential propagation full step
            ! Diagonal potential matrix half step
            do j = 1, Nstates         
                this%psi_ges(:,j) = this%psi_ges(:,j) * split_operator%vprop(:,j)
            end do
            
            ! Off-diagonal potential matrix full step i.e interstates coupling
            !$OMP PARALLEL DO DEFAULT(NONE) FiRSTPRiVATE(tout, psi_Nstates, psi_Nstates1) &
            !$OMP FiRSTPRiVATE(psi_Nstates_real, psi_Nstates_imag) &
            !$OMP SHARED(mu_all, E, this, k, Nstates, NR, dt)
            do i = 1, NR
                call pulse2(tout, mu_all(:,:,i), E(k))   
                psi_Nstates(:) = this%psi_ges(i,:)
                call zgemv('N', int(Nstates), int(Nstates), (1._dp, 0._dp),  &
                    & tout, size(tout,dim=1), psi_Nstates, 1, (0._dp,0._dp), &
                    & psi_Nstates1, 1)
                this%psi_ges(i,:) = psi_Nstates1(:)
            end do
            !$OMP END PARALLEL DO

            ! Diagonal potential matrix half step
            do j = 1, Nstates         
                this%psi_ges(:,j) = this%psi_ges(:,j) * split_operator%vprop(:,j)
            end do

            ! Kinetic propagation half step
            call split_operator%split_operator(this%psi_ges)

            ! On the fly analysis and writing to files
            ! norm
            call integ(this%psi_ges, norm)
            write(this%norm_1d_tk,*) time(k) * au2fs, norm
            ! localized states norm
            call this%localized_states_norm(normPn)
            write(this%norm_pn_1d_tk,*) time(k) * au2fs, normPn
            ! expected position
            call this%expected_position(norm, evr)
            write(this%avgR_1d_tk,*) time(k) * au2fs, evr(1:Nstates)
            ! vibrational populations (writing only for the ground state)
            call this%vib_pop_analysis(vib_pop)
            write(this%vibpop_1d_tk,*) time(k) *au2fs, vib_pop(1:Vstates(1),1)

            ! density
            ! Coordinate space density 
            if(mod(k,100).eq.0) then
                do i = 1, NR, 4   
                    write(this%dens_1d_tk,*) time(k) *au2fs, R(i), abs(this%psi_ges(i,1))**2
                    write(this%ex_dens_1d_tk,*) time(k) *au2fs, R(i), abs(this%psi_ges(i,2:Nstates))**2
                end do 
                write(this%dens_1d_tk,*)
                write(this%ex_dens_1d_tk,*)
            end if
   
            ! Momentum space density
            ! if (mod(k,100) .eq. 0) then
            !    do i = NR/2 +1, NR, 4   
            !        write(this%Pdens_1d_tk,*) time(k) *au2fs, pR(i), abs(this%psi_ges_p(i,:))**2
            !    end do 
            !    do i = 1, NR/2, 4   
            !        write(this%Pdens_1d_tk,*) time(k) *au2fs, pR(i), abs(this%psi_ges_p(i,:))**2
            !    end do 
            !    write(Pdens_1d_tk,*)
            ! endif

            !-----------------------------------
            ! Cotinuum part treatment
            call this%continuum_prop(split_operator)
            
            ! Absorbed density norm
            call integ(this%psi_outR, norm_outR)
            do J = 1, Nstates
                norm_SE_outR(J) = sum(this%psi_outR_inc(:,J))*dR
            enddo
            write(this%psi_outR_norm_1d_tk,*) time(k), norm_outR

            ! Absorbed momentum density
            ! Momentum space density
            if (mod(K,100) .eq. 0) then
                do i = NR/2 +1, NR, 4
                    write(this%psi_outR_Pdens_1d_tk,*) time(k) *au2fs, & 
                        & pR(i), abs(this%psi_outR(i,:))**2
                end do
                do i = 1, NR/2, 4
                    write(this%psi_outR_Pdens_1d_tk,*) time(k) *au2fs, &
                        & pR(i), abs(this%psi_outR(i,:))**2
                end do 
                write(this%psi_outR_Pdens_1d_tk,*)
            endif
        
        end do timeloop




    end subroutine time_evolution

    subroutine localized_states_norm(this, normPn)
        use global_vars, only: NR, Nstates
        class(time_prop), intent(inout) :: this
        real(dp):: normPn(Nstates), time
        complex(dp) :: psi_loc(NR,Nstates)
        normPn = 0._dp

        psi_loc(:,1) = 1._dp/sqrt(2._dp)*(this%psi_ges(:,1) + this%psi_ges(:,2))
        psi_loc(:,2) = 1._dp/sqrt(2._dp)*(this%psi_ges(:,1) - this%psi_ges(:,2))

        call integ(psi_loc, normPn)
        
    end subroutine localized_states_norm

    subroutine expected_position(this, norm, evr)
        use global_vars, only: NR, Nstates, dR, R
        class(time_prop), intent(inout) :: this
        real(dp):: evR(Nstates), norm(Nstates)
        integer:: N

        evr = 0._dp
        do N = 1, Nstates
            evR(N) = sum(abs(this%psi_ges(:,N))**2 * R(:))
        enddo
        evR = evR * dR

        do N = 1,Nstates
            if (norm(N).ge.1.e-8_dp) then
                evR(N) = evR(N)/norm(N)
            endif
        enddo

    end subroutine expected_position

    subroutine vib_pop_analysis(this, vib_pop)
        use global_vars, only: Vstates, Nstates, guess_vstates, dR
        class(time_prop), intent(inout) :: this
        integer:: N, L
        real(dp) :: vib_pop(guess_vstates,Nstates)
        complex(dp), allocatable :: psi_chi(:)
        allocate(psi_chi(guess_vstates))
        vib_pop = 0._dp
        do N = 1, Nstates
            do L=1,vstates(N)
                psi_chi(L)=0._dp
                psi_chi(L)=sum(this%chi0(:,L,N) * (this%psi_ges(:,N)))
                psi_chi(L)=psi_chi(L)*dR
                vib_pop(L,N)=abs(psi_chi(L))**2
            enddo
        enddo
        deallocate(psi_chi)
    end subroutine vib_pop_analysis

    subroutine continuum_prop(this, split_operator)
        use global_vars, only: NR, Nstates, dt, adb
        use data_au, only: im
        use FFTW3
        class(time_prop), intent(inout) :: this
        type(split_operator_type), intent(in) :: split_operator
        integer:: i, J
        complex(dp), allocatable :: psi_outR1(:,:)

        allocate(psi_outR1(NR,Nstates))
      
        psi_outR1 = (0._dp, 0._dp)
        do J = 1, Nstates
            this%psi_outR(:,J) = this%psi_outR(:,J) * split_operator%kprop_full(:) &
                & * exp(-im*dt*adb(NR-this%i_cpmR,J))
            psi_outR1(:,J) = this%psi_ges(:,J) * (1._dp-this%abs_func(:))
            this%psi_ges(:,J) = this%psi_ges(:,J) * this%abs_func(:) ! psi_ges = psi_nondiss
        enddo

        do J = 1, Nstates
            split_operator%psi = (0._dp, 0._dp)
            split_operator%psi(:) = psi_outR1(:,J)
            call fftw_execute_dft(split_operator%planF, split_operator%psi, split_operator%psi)
            split_operator%psi = split_operator%psi/sqrt(dble(NR))
            psi_outR1(:,J) = split_operator%psi(:)
        enddo
        this%psi_outR = this%psi_outR + psi_outR1
        this%psi_outR_inc = this%psi_outR_inc + abs(psi_outR1)**2

    end subroutine continuum_prop

    subroutine post_prop_analysis(this)
        class(time_prop), intent(inout) :: this
        print*
        print*, "Post-propagation analysis..."
        print*
        ! implementation of post-propagation analysis
    end subroutine post_prop_analysis
    !------------------------------------------------------------------
    subroutine integ(psi, norm)
         
        use global_vars, only:NR, Nstates, dR, dp
        implicit none
        integer i, J
         
        real(dp) norm(Nstates)
        complex(dp) psi(NR,Nstates)
         
        norm = 0.d0
    
        do J = 1, Nstates
            norm(J)= sum(abs(psi(:,J))**2)   
        end do
    
        norm = norm * dR
        return 
    end subroutine

    subroutine Boltzmann_distribution(N, E, Boltzmann_populations)
        use global_vars, only: temperature, dp, Vstates, guess_vstates
        use data_au, only: kB

        integer:: N, V
        real(dp), intent(in):: E(guess_vstates)
        real(dp), allocatable, intent(out):: Boltzmann_populations(:)
        real(dp):: total_pop


           
        total_pop = sum(exp(-E(:)/(kB*temperature)))
        print*, "total populations =", total_pop
   
        allocate(Boltzmann_populations(Vstates(N)))
        do V = 1, Vstates(N)
            if (V .gt. 1) then
                Boltzmann_populations(V) = exp(-(E(V)-E(1))/(kB*temperature))!/total_pop
            else
                Boltzmann_populations(V)=1._dp
            endif
            print*, "state", V-1, "population ratio =", Boltzmann_populations(V)
        enddo
   
        total_pop = sum(Boltzmann_populations(:))
        print*, "sum of ratios =", total_pop
        Boltzmann_populations = Boltzmann_populations/total_pop
        total_pop = sum(Boltzmann_populations(:))
        print*, "total population =", total_pop
   
        do V =1, Vstates(N)
            print*, "state", V-1, "population probability =", Boltzmann_populations(V)
        enddo  

    end subroutine Boltzmann_distribution

    subroutine complex_absorber_function(v_abs, f)
        use global_vars, only: NR, dp, dR, dt, R
        use data_au, only: im
        use pot_param, only: cpmR
        
        integer i
        real(dp):: a, eps, V_abs(NR), n, R0, p
        complex(dp):: iV_abs(NR), f(NR)
    
        eps = epsilon(a) 
        print*, "Lower limit of the precision:", eps
        n = 4 ! power of absorber function
        R0 = R(NR)- cpmR ! start of the absorber
        p = 20._dp ! optimal absorption momentum
        a = -log(eps) *(n+1) *p / (2*(R(NR)-R0)**(n+1))
        print*, "Absorber prefactor a:", a
   
        do i = 1, NR
            if (R(i) .gt. abs(R0)) then
                V_abs(i) = a*(R(i)-R0)**n
            else
                V_abs(i) = 0._dp
            endif
        enddo
        f(:) = exp(-dt *V_abs(:))
    end subroutine

    subroutine mask_function_cos(cof)
        use global_vars, only:NR, R, dR, dp
        use data_au
        use pot_param
   
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
    
        integer :: J
        real(dp):: cof(NR),c
  
        c=1.00d0
        do j = 1, NR
            R(J) = R0 + (j - 1) * dR
            cof(j)=1.0d0/(1.0d0+exp(c*(R(J)-Rend+cpmR)))
        end do
  
        return
    end subroutine

    subroutine pulse2(tout,mu,E)
        use global_vars, only:dt, Nstates,kap, lam, dp, idp
        use data_au, only:im
        use blas_interfaces_module, only : zgemv, dgemv
   
        integer:: i, J
        real(dp):: w, u, uv, d, mu(Nstates,Nstates), q, E
        integer info, Lwork
   
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
        do i=1, Nstates-1
            do J=i+1, Nstates
                u(i,J)= -kap*mu(i,J) * E
                u(J,i)= -kap*mu(i,J) * E
            enddo
        enddo
        !print*, "u"
        !write( * , * ) ((u(i,j),j=1,Nstates), i=1,Nstates )
    
        uv = 0._dp
        uv = u
   
        Lwork=-1
        JOBZ='V'
        call dsyev(JOBZ,'U', int(Nstates), uv, int(Nstates), &
            & d, work, Lwork,info)
        ! call jacobi(u,Nstates,d)
        Lwork = min(1000, int(work(1)))
        
        !---- Solve eigenproblem.
        call dsyev('V', 'U', int(Nstates), uv, int(Nstates), &
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
        call zgemm('N', 'N', int(Nstates), int(Nstates), &
            & int(Nstates),(1._dp,0._dp), p, size(p,dim=1), &
            & b, size(b,dim=1), (0._dp,0._dp), z, size(z,dim=1))
        !print*, "z"
        !write( * , * ) ((z(i,j),j=1,Nstates), i=1,Nstates )
   
        pT = transpose(p)
        !tout = matmul(z,transpose(p))
        call zgemm('N', 'N', int(Nstates), int(Nstates), &
            & int(Nstates),(1._dp,0._dp), z, size(z,dim=1), &
            & pT, size(pT,dim=1), (0._dp,0._dp), tout, size(tout,dim=1))
        !print*, "tout"
        !write( * , * ) ((tout(i,j),j=1,Nstates), i=1,Nstates )
   
   return
   end subroutine pulse2

end module propagation_mod





subroutine propagation_1D(E, A)
use global_vars
use pot_param
use data_au
use FFTW3
use omp_lib
use blas_interfaces_module, only : zgemv, dgemv

 implicit none   
 
 integer i, J, K,i_cpmR, ii, io
 integer L, M, N, void, v
 integer v_ini_check
 type(C_PTR) planF, planB
 integer eR
 real(dp) dummy, dummy2, dummy3
 character(150):: filepath

 real(dp) dt2 !time
 real(dp) c, sp, deR 
 real(dp):: normpn(Nstates), spec(NR,Nstates)
 real(dp) E(Nt), E1, norm(Nstates), E21, E22
 real(dp):: norm_outR(Nstates), norm_outR1(Nstates), norm_SE_outR(Nstates)
 real(dp):: norm_diss(Nstates), norm_bound(Nstates)
 real(dp) A(Nt)
 real(dp) evR(Nstates), epr(Nstates),momt(Nstates), tot_momt
 real(dp) :: H_ac(NR,Nstates), Energy_axis(NR)
! double precision :: Boltzmann_populations(Vstates)
 real(dp) :: norm_overlap(Nstates),norm_outP(Nstates)
 real(dp) :: norm_gesP(Nstates), norm_gP_over(Nstates)
 real(dp) :: vib_pop(guess_vstates,Nstates) 
 integer, allocatable, dimension(:) ::vib_ini
 real(dp), allocatable, dimension(:) ::vib_dist
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
 real(dp), allocatable, dimension(:,:):: psi_outR_inc
 complex(dp), allocatable, dimension(:,:):: psi_gesP
 complex(dp), allocatable, dimension(:):: psi_chi

 integer:: chi0_tk, vstates_tk
 integer:: vib_dist_tk, Boltzmann_dist_tk
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
 allocate(psi_outR1(NR,Nstates), psi_outR_inc(NR,Nstates))
 allocate(psi_chi(guess_vstates),psi_diss(NR,Nstates),psi_bound(NR,Nstates))

 print*
 print*,'Tuning FFTW...'
 void=fftw_init_threads( )
 if (void .eq. 0) then
    print*, 'Error in fftw_init_threads, quitting'
 else
    print*, 'number of threads found =', void
 endif

 select case(prop_par_FFTW)
 case ("parallel")
   call fftw_plan_with_nthreads(omp_get_max_threads())    
   planF = fftw_plan_dft_1d(NR, psi, psi, FFTW_FORWARD,FFTW_MEASURE)
   planB = fftw_plan_dft_1d(NR, psi, psi, FFTW_BACKWARD,FFTW_MEASURE)

 case default ! single thread FFTW
   planF = fftw_plan_dft_1d(NR, psi, psi, FFTW_FORWARD,FFTW_MEASURE)
   planB = fftw_plan_dft_1d(NR, psi, psi, FFTW_BACKWARD,FFTW_MEASURE)
 end select  
             
 print*,'Done.'
 print*  
  
 write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "Bound-vibstates_in_Nthstates.out"
 open(newunit=vstates_tk,file=filepath,status='unknown')
 chi0 = 0._dp
 do N = 1, Nstates 
  read(vstates_tk,*) ii, Vstates(N)
  write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "BO_Electronic-state-g", &
         & int(N-1), "_vibstates.out"
  open(newunit=chi0_tk,file=filepath,status='unknown')
  print*, "NR:", NR, "Vstates(N)", Vstates(N)
  do i=1,NR
    read(chi0_tk,*) dummy, chi0(i,1:Vstates(N),N)
    print*, i, dummy, chi0(i,Vstates(N),N)
  enddo 
  close(chi0_tk)
 enddo
 close(vstates_tk)

 psi_ges = (0._dp,0._dp)   

 select case(initial_distribution) 
 case("single vibrational state")
  print*, "initial wavefunction in..."
  print*, N_ini-1, "electronic state and in", v_ini-1, "vibrational state"
  do i = 1, NR
    psi_ges(i,N_ini) = chi0(i,v_ini,N_ini) 
  enddo  

 case("gaussian distribution")
  print*, "initial wavefunction in..."
  print*, N_ini-1, "electronic state and with a Gaussian distribution centered around",&
         & Ri_tdse, "a.u. \n with deviation of", kappa_tdse, "."
  do i = 1, NR
    psi_ges(i,1)=exp(kappa_tdse*(R(i)-Ri_tdse)**2) 
  enddo

! case ("input dist")
!  allocate(v_dist_ini(N_ini))
!  v_ini_check = 0
!  do N = 1, N_ini
!    write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "vib_dist_", int(N-1), ".out"
!    open(newunit=vib_dist_tk,file=filepath,status='unknown')
!    do 
!     read(vib_dist_tk,*,iostat=io)
!     if (io /= 0) exit
!      v_ini_check = v_ini_check + 1
!    enddo
!    if (v_ini /= v_ini_check) then
!      write(*,'(a,a,i0,a)') "Number of vibstates in 'vib_dist_", N, ".out' not equal to input" 
!    endif
!  enddo
!  allocate(vib_ini(guess_vstates,N_ini), vib_dist(guess_vstates, N_ini))
!  do N = 1, N_ini
!    do v = 1, v_dist_ini(N)
!      read(vib_dist_tk,*) vib_ini(v,N), vib_dist(v,N) 
!    enddo
!  enddo
!  do i = 1, NR
!   do v = 1, v_ini
!    do N = 1, N_ini
!     psi_ges(i,N) = psi_ges(i,N) + vib_dist(vib_ini(v,N),N) * chi0(i,vib_ini(v,N),N)
!    enddo
!   enddo
!  enddo

 case ("Boltzmann dist")
  do N = 1, N_ini
   call Boltzmann_distribution(N,vib_dist)
   do v = 1, Vstates(N)
    do i = 1, NR
      psi_ges(i,N) = psi_ges(i,N) + vib_dist(v) * chi0(i,v,N)
    enddo
   enddo
  enddo

 case default
  do i = 1, NR
   psi_ges(i,1) = chi0(i,1,1)
  enddo

 end select
 
 ! cpm = 3.d0/ au2a
 ! call cutoff_cos(cpm,cof)

 i_cpmR = minloc(abs(R(:)-cpmR),1) - 50
 print*, "i_cpmR =", i_cpmR, ", NR-i_cpmR", NR-i_cpmR
 print*, "R(NR-i_cpmR) =", R(NR-i_cpmR)

 select case(absorber)
 case ("CAP")
   call Complex_absorber_function(V_abs, exp_abs)
   abs_func = exp_abs
 case("mask")
   call mask_function_ex(cof)
   abs_func = cof
 end select
 
 call integ(psi_ges, norm) 
 print*,'initial norm:', norm 
  
 psi_ges(:,1) = psi_ges(:,1) / sqrt(norm(1))
 
 call integ(psi_ges, norm)
 print*,'norm after normalization:', norm 
 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "psi0_1d.out"
 open(newunit=psi_1d_tk,file=filepath,status='unknown') 
 write(filepath, '(a,a)') adjustl(trim(output_data_dir)), "absorber_function.out"
 open(newunit=cof_1d_tk,file=filepath,status='unknown') 
 do i = 1, NR
  write(psi_1d_tk,*) R(i), abs(psi_ges(i,1))**2
  write(cof_1d_tk,*) R(i), abs_func(i)
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
 psi_outR_inc = 0._dp
 norm = 0._dp
 normPn = 0._dp
 norm_outR = 0._dp
 do i = 1, NR
   kprop(i) = exp(-im *dt * pR(i)*pR(i) /(4._dp*m_red))  ! pR**2 /2 * red_mass UND Half time step
   kprop1(i) =exp(-im *dt * pR(i)*pR(i) /(2._dp*m_red))
 end do 
 
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
     call fftw_execute_dft(planF, psi, psi)
     psi = psi * kprop
     call fftw_execute_dft(planB, psi, psi)
     psi = psi / dble(NR)
     psi_ges(:,J) = psi(:)      
   end do
   
   do j = 1, Nstates   
!     psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * (adb(i,j)+kap*(mu_all(i,J,J) &
!              & +0.5d0*R(i))*E(K)+((2-kap)*mr+kap-1)*R(i)*E(K)))!+H_ac(i,j))) !         
     psi_ges(:,j) = psi_ges(:,j) * exp(-im * 0.5_dp * dt * (adb(:,j))) !+0.8d0*R(i)*E(K)))!+H_ac(i,j))) !         
   end do
   
   !$OMP PARALLEL DO DEFAULT(NONE) FiRSTPRiVATE(tout, psi_Nstates, psi_Nstates1) &
   !$OMP FiRSTPRiVATE(psi_Nstates_real, psi_Nstates_imag) &
   !$OMP SHARED(mu_all, E, psi_ges, K, Nstates, NR, dt)
   do i = 1, NR
     call pulse2(tout, mu_all(:,:,i), E(K)) 
!     psi_ges(i,1:Nstates) = matmul(tout(1:Nstates,1:Nstates),psi_ges(i,1:Nstates))  
     psi_Nstates(:) = psi_ges(i,:)
     call zgemv('N', int(Nstates), int(Nstates), (1._dp, 0._dp),  &
             & tout, size(tout,dim=1), psi_Nstates, 1, (0._dp,0._dp), &
             & psi_Nstates1, 1)
     psi_ges(i,:) = psi_Nstates1(:)
   end do
   !$OMP END PARALLEL DO
    

   do j = 1, Nstates   
!     psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * (adb(i,j)+kap*(mu_all(i,J,J) &
!              & +0.5d0*R(i))*E(K)+((2-kap)*mr+kap-1)*R(i)*E(K)))!+H_ac(i,j))) !         
     psi_ges(:,j) = psi_ges(:,j) * exp(-im * 0.5_dp * dt * (adb(:,j))) !+0.8d0*R(i)*E(K)))!+H_ac(i,j))) !         
   end do
    
   do J = 1,Nstates
     psi = (0._dp, 0._dp)
     psi(:) = psi_ges(:,J)  ! Hilfsgroesse
     call fftw_execute_dft(planF, psi, psi) 
     psi = psi * kprop 
     psi = psi /sqrt(dble(nr))     
     epr(j) =  sum(abs(psi(:))**2 * pr(:))
     epr(j) = epr(j) * dr
     psi_ges_p(:,J) = psi(:)
     call fftw_execute_dft(planB, psi, psi)
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
    do i = 1, NR, 4   
      write(dens_1d_tk,*) time *au2fs, R(i), abs(psi_ges(i,1))**2
      write(ex_dens_1d_tk,*) time *au2fs, R(i), abs(psi_ges(i,2:Nstates))**2
    end do 
    write(dens_1d_tk,*)
    write(ex_dens_1d_tk,*)
   end if
   
   ! Momentum space density
   if (mod(K,100) .eq. 0) then
    do i = NR/2 +1, NR, 4   
      write(Pdens_1d_tk,*) time *au2fs, pR(i), abs(psi_ges_p(i,:))**2
    end do 
    do i = 1, NR/2, 4   
      write(Pdens_1d_tk,*) time *au2fs, pR(i), abs(psi_ges_p(i,:))**2
    end do 
    write(Pdens_1d_tk,*)
   endif

 ! -------------Continuum treatment ---------------------
  
   psi_outR1 = (0._dp, 0._dp)
   do J = 1, Nstates
     psi_outR(:,J) = psi_outR(:,J) * kprop1(:) * exp(-im*dt*adb(NR-i_cpmR,J))
     psi_outR1(:,J) = psi_ges(:,J) * (1._dp-abs_func(:))
     psi_ges(:,J) = psi_ges(:,J) * abs_func(:) ! psi_ges = psi_nondiss
   enddo
   call integ(psi_outR, norm_outR)
   do J = 1, Nstates
     norm_SE_outR(J) = sum(psi_outR_inc(:,J))*dR
   enddo
   write(psi_outR_norm_1d_tk,*) time*au2fs, norm_outR(:), norm_SE_outR(:)
   do J = 1, Nstates
     psi = (0._dp, 0._dp)
     psi(:) = psi_outR1(:,J)
     call fftw_execute_dft(planF, psi, psi)
     psi = psi/sqrt(dble(NR))
     psi_outR1(:,J) = psi(:)
   enddo
   psi_outR = psi_outR + psi_outR1
   psi_outR_inc = psi_outR_inc + abs(psi_outR1)**2
   ! Momentum space density
   if (mod(K,100) .eq. 0) then
    do i = NR/2 +1, NR, 4
      write(psi_outR_Pdens_1d_tk,*) time *au2fs, pR(i), abs(psi_outR(i,:))**2
    end do
    do i = 1, NR/2, 4
      write(psi_outR_Pdens_1d_tk,*) time *au2fs, pR(i), abs(psi_outR(i,:))**2
    end do 
    write(psi_outR_Pdens_1d_tk,*)
   endif

! ------------   

end do timeloop
  
!Final vibrational population in ground state
  do N =1, Nstates
    print*
    print*, "Final vibrational population in the elec. state", int(N-1)
    print*, "Number of vibrational states:", vstates(N)
    psi_bound = 0._dp
    do J = 1,vstates(N)
      psi_chi(J) = 0._dp
      do i = 1, NR
        psi_chi(J) = psi_chi(J)+ chi0(i,J,N) * (psi_ges(i,N))
      enddo
      psi_chi(J) = psi_chi(J)*dR
      do i = 1,NR
        psi_bound(i,N) = psi_bound(i,N)+psi_chi(J)*chi0(i,J,N)
      enddo
      print*, 'Vibpop (',J-1,') =', abs(psi_chi(J))**2
      norm_bound(N) = sum(abs(psi_bound(:,N))**2)*dR
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
    call fftw_execute_dft(planF, psi, psi)
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
    do i=NR/2 +1,NR
      write(momt_spectra_un_tk,*) pR(i), abs(psi_diss(i,N))**2, abs(psi_outR(i,N))**2
      write(momt_spectra_tk,*) pR(i), abs(psi_diss(i,N)/sqrt(norm_diss(N)))**2, &
              & abs(psi_outR(i,N)/sqrt(norm_diss(N)))**2
    enddo
    do i=1,NR/2
      write(momt_spectra_un_tk,*) pR(i), abs(psi_diss(i,N))**2, abs(psi_outR(i,N))**2
      write(momt_spectra_tk,*) pR(i), abs(psi_diss(i,N)/sqrt(norm_diss(N)))**2, &
              & abs(psi_outR(i,N)/sqrt(norm_outR(N)))**2
      write(KER_spectra_un_tk,*) pR(i)**2/(2*m_red), m_red*abs(psi_diss(i,N))**2 / pR(i),&
              & m_red*abs(psi_outR(i,N))**2 / pR(i), m_red*abs(psi_outR_inc(i,N))/pR(i) 
      write(KER_spectra_tk,*) pR(i)**2/(2*m_red), m_red*abs(psi_diss(i,N)/sqrt(norm_diss(N)))**2 / pR(i), &
              & m_red*abs(psi_outR(i,N)/sqrt(norm_outR(N)))**2 / pR(i) 
    enddo
    close(momt_spectra_un_tk)
    close(KER_spectra_un_tk)
    close(momt_spectra_tk)
    close(KER_spectra_tk)
  enddo

  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Total_KER_spectra.out"
  open(newunit=KER_spectra_un_tk,file=filepath,status='unknown')
  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Total_momt_spectra.out"
  open(newunit=momt_spectra_un_tk,file=filepath,status='unknown')
  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Total_KER_spectra_normalized.out"
  open(newunit=KER_spectra_tk,file=filepath,status='unknown')
  write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Total_momt_spectra_normalized.out"
  open(newunit=momt_spectra_tk,file=filepath,status='unknown')
  do i=NR/2 +1,NR
    write(momt_spectra_un_tk,*) pR(i), sum(abs(psi_diss(i,:))**2), sum(abs(psi_outR(i,:))**2), sum(abs(psi_outR_inc(i,:)))
    write(momt_spectra_tk,*) pR(i), sum(abs(psi_diss(i,:)/sqrt(norm_diss(:)))**2), &
              & sum(abs(psi_outR(i,:)/sqrt(norm_outR(:)))**2)
  enddo
  do i=1,NR/2
    write(momt_spectra_un_tk,*) pR(i), sum(abs(psi_diss(i,:))**2), sum(abs(psi_outR(i,:))**2), sum(abs(psi_outR_inc(i,:)))
    write(momt_spectra_tk,*) pR(i), sum(abs(psi_diss(i,:)/sqrt(norm_diss(:)))**2), &
              & sum(abs(psi_outR(i,:)/sqrt(norm_outR(:)))**2)
    write(KER_spectra_un_tk,*) pR(i)**2/(2*m_red), m_red*sum(abs(psi_diss(i,:))**2) / pR(i), &
              & m_red*sum(abs(psi_outR(i,:))**2) / pR(i), m_red*sum(abs(psi_outR_inc(i,:)))/pR(i)
    write(KER_spectra_tk,*) pR(i)**2/(2*m_red), m_red*sum(abs(psi_diss(i,:)/sqrt(norm_diss(:)))**2) /pR(i), &
              & m_red*sum(abs(psi_outR(i,:)/sqrt(norm_outR(:)))**2) / pR(i)
  enddo
  close(momt_spectra_un_tk)
  close(KER_spectra_un_tk)
  close(momt_spectra_tk)
  close(KER_spectra_tk)

 
! ------------   
                                      
 call fftw_destroy_plan(planF)
 call fftw_destroy_plan(planB)
   
 deallocate(psi, kprop, psi_ges, cof, psi_loc, psi_outR1,psi_outR, psi_gesP)

 contains

!_________________________________________________________

   subroutine Boltzmann_distribution(N, Boltzmann_populations)
   use global_vars, only: temperature, dp
   implicit none
    integer:: N, V
    real(dp), allocatable, intent(out):: Boltzmann_populations(:)
    real(dp):: total_pop
           
    total_pop = sum(exp(-E(:)/(kB*temperature)))
    print*, "total populations =", total_pop
   
    allocate(Boltzmann_populations(Vstates(N)))
    do V = 1, Vstates(N)
      if (V .gt. 1) then
      Boltzmann_populations(V) = exp(-(E(V)-E(1))/(kB*temperature))!/total_pop
      else
      Boltzmann_populations(V)=1._dp
      endif
      print*, "state", V-1, "population ratio =", Boltzmann_populations(V)
    enddo
   
    total_pop = sum(Boltzmann_populations(:))
    print*, "sum of ratios =", total_pop
    Boltzmann_populations = Boltzmann_populations/total_pop
    total_pop = sum(Boltzmann_populations(:))
    print*, "total population =", total_pop
   
    do V =1, Vstates(N)
      print*, "state", V-1, "population probability =", Boltzmann_populations(V)
    enddo  

   end subroutine

!---------------------------------------------------------------

   subroutine complex_absorber_function(v_abs, f)
   use global_vars, only: NR, dp, dR, dt, R
   use data_au, only: im
   use pot_param, only: cpmR
    implicit none
    integer i
    real(dp):: a, eps, V_abs(NR), n, R0, p
    complex(dp):: iV_abs(NR), f(NR)
    
    eps = epsilon(a) 
    print*, "Lower limit of the precision:", eps
    n = 4 ! power of absorber function
    R0 = R(NR)- cpmR ! start of the absorber
    p = 20._dp ! optimal absorption momentum
    a = -log(eps) *(n+1) *p / (2*(R(NR)-R0)**(n+1))
    print*, "Absorber prefactor a:", a
   
    do i = 1, NR
     if (R(i) .gt. abs(R0)) then
       V_abs(i) = a*(R(i)-R0)**n
     else
       V_abs(i) = 0._dp
     endif
    enddo
    f(:) = exp(-dt *V_abs(:))
   end subroutine
!------------------------------------------------------------------
   subroutine integ(psi, norm)
         
   use global_vars, only:NR, Nstates, dR, dp
    implicit none
    integer i, J
         
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
    integer info, Lwork
   
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
   do i=1, Nstates-1
    do J=i+1, Nstates
          u(i,J)= -kap*mu(i,J) * E
          u(J,i)= -kap*mu(i,J) * E
    enddo
   enddo
   !print*, "u"
   !write( * , * ) ((u(i,j),j=1,Nstates), i=1,Nstates )
    
   uv = 0._dp
   uv = u
   
   Lwork=-1
   JOBZ='V'
    call dsyev(JOBZ,'U', int(Nstates), uv, int(Nstates), &
           & d, work, Lwork,info)
   ! call jacobi(u,Nstates,d)
    Lwork = min(1000, int(work(1)))
   !     Solve eigenproblem.
    call dsyev('V', 'U', int(Nstates), uv, int(Nstates), &
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
    call zgemm('N', 'N', int(Nstates), int(Nstates), &
            & int(Nstates),(1._dp,0._dp), p, size(p,dim=1), &
            & b, size(b,dim=1), (0._dp,0._dp), z, size(z,dim=1))
   !print*, "z"
   !write( * , * ) ((z(i,j),j=1,Nstates), i=1,Nstates )
   pT = transpose(p)
   !tout = matmul(z,transpose(p))
    call zgemm('N', 'N', int(Nstates), int(Nstates), &
            & int(Nstates),(1._dp,0._dp), z, size(z,dim=1), &
            & pT, size(pT,dim=1), (0._dp,0._dp), tout, size(tout,dim=1))
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




