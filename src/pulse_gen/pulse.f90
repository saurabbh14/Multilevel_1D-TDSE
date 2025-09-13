!> This module contains the subroutine for generating the laser pulse
!> and the envelope functions.
module pulse_mod
    use global_vars, only: dp, Nt, output_data_dir, dt
    use data_au
    use FFTW3
    use omp_lib
    implicit none

    type :: pulse_param
        character(150) :: envelope_shape_laser1, envelope_shape_laser2
        real(dp) :: tp1, fwhm, t_mid1, rise_time1
        real(dp) :: tp2, t_mid2, rise_time2
        real(dp) :: e01, e02, phi1, phi2
        real(dp) :: lambda1, lambda2
        real(dp) :: omega1, omega2
        real(dp) :: pulse_offset1, pulse_offset2
        real(dp), allocatable :: El(:), Al(:)
        real(dp), allocatable :: E21(:), E22(:)
        real(dp), allocatable :: A21(:), A22(:)
        real(dp), allocatable :: g1(:), g2(:)
  
    contains

        procedure :: read => read_pulse_params
        procedure :: initialize => initialize_pulse_param
        procedure :: param_print => print_pulse_param
        procedure :: generate => generate_pulse
        procedure :: write_to_file => write_pulse_to_file
        procedure :: deallocate_env => deallocate_envelope
        procedure :: deallocate_field => deallocate_field
        procedure :: deallocate_all => deallocate_all
        procedure :: spectra => field_spectra
    end type pulse_param

contains
  
    subroutine read_pulse_params(this, input_path)
        class(pulse_param), intent(inout) :: this
        character(2000), intent(in) :: input_path
        ! file tokens
        integer :: input_tk

        ! Intermediate variables for pulse_param components
        character(150) :: envelope_shape_laser1, envelope_shape_laser2
        real(dp) :: lambda1, lambda2, tp1, tp2, t_mid1, t_mid2
        real(dp) :: E01, E02, phi1, phi2, rise_time1, rise_time2

        namelist /laser_param/envelope_shape_laser1, envelope_shape_laser2, &
        & lambda1, lambda2, tp1, tp2, t_mid1, t_mid2, E01, E02, & 
        & phi1, phi2, rise_time1, rise_time2

        open(newunit=input_tk, file=adjustl(trim(input_path)), status='old')
        read(input_tk,nml=laser_param)
        close(input_tk)

        ! Assign values to the pulse_param components
        this%envelope_shape_laser1 = envelope_shape_laser1
        this%envelope_shape_laser2 = envelope_shape_laser2
        this%lambda1 = lambda1
        this%lambda2 = lambda2
        this%tp1 = tp1
        this%tp2 = tp2
        this%t_mid1 = t_mid1
        this%t_mid2 = t_mid2
        this%E01 = E01
        this%E02 = E02
        this%phi1 = phi1
        this%phi2 = phi2
        this%rise_time1 = rise_time1
        this%rise_time2 = rise_time2

    end subroutine read_pulse_params
    
    ! A subroutine for printing the pulse parameters
    subroutine print_pulse_param(this)
        class(pulse_param), intent(in) :: this
        print*, "Laser parameters:"
        print*, "Laser #1:"
        print*, "Envelope shape:", trim(this%envelope_shape_laser1)
        print*, "Lambda:", this%lambda1, "nm"
        print*, "Electric field strength:", this%E01, "a.u."
        print*, "Pulse width (tp):", this%tp1, "fs"
        print*, "Pulse midpoint:", this%t_mid1, "fs"
        print*, "phi1:", this%phi1, "pi"
        print*, "Rise time:", this%rise_time1, "fs"
        print*, "Laser #2:"
        print*, "Envelope shape:", trim(this%envelope_shape_laser2)
        print*, "Lambda:", this%lambda2, "nm"
        print*, "Electric field strength:", this%E02, "a.u."
        print*, "Pulse width (tp):", this%tp2, "fs"
        print*, "Pulse midpoint:", this%t_mid2, "fs"
        print*, "phi2:", this%phi2, "pi"
        print*, "Rise time:", this%rise_time2, "fs"
        print*, "------------------------------------------------------"
        print*, "Final pulse parameters:"
        print*, "Wavelength 1 =", sngl(this%lambda1), "nm"
        print*, "Phase 1 =", sngl(this%phi1)
        print*, "Field strength =", sngl(this%e01), "a.u.", sngl(this%e01*e02au), "V/m"
        print*, "Intensity =", sngl(this%e01**2*3.509e16_dp), "W/cm2"
        print*, "Wavelength 2 =", sngl(this%lambda2), "nm"
        print*, "Phase 2 =", sngl(this%phi2)
        print*, "Field strength =", sngl(this%e02), "a.u.", sngl(this%e02*e02au), "V/m"
        print*, "Intensity =", sngl(this%e02**2*3.509e16_dp), "W/cm2"
        print*, "Wave duration =", sngl(this%tp1*au2fs), "fs"
        print*, "------------------------------------------------------"
    end subroutine print_pulse_param

    ! A subroutine for initializing the pulse parameters
    subroutine initialize_pulse_param(this)
        class(pulse_param), intent(inout) :: this
        ! Initialize the pulse parameters
        this%tp1 = this%tp1 / au2fs  
        this%tp2 = this%tp2 / au2fs  
        this%t_mid1 = this%t_mid1 / au2fs   
        this%t_mid2 = this%t_mid2 / au2fs   
        this%rise_time1 = this%rise_time1 / au2fs
        this%rise_time2 = this%rise_time2 / au2fs
        this%omega1 = (1._dp / (this%lambda1 * 1.e-7_dp)) * cm2au
        this%omega2 = (1._dp / (this%lambda2 * 1.e-7_dp)) * cm2au
        this%phi1 = this%phi1 * pi
        this%phi2 = this%phi2 * pi
    end subroutine initialize_pulse_param


    ! A subroutine for defining the field 
    subroutine generate_pulse(this)
        class(pulse_param), intent(inout) :: this
        integer :: K
        real(dp) :: time
        real(dp) :: A01, A02

        print*
        print*, "Pulse generation..."
    
        ! Initialize arrays
        allocate(this%El(Nt), this%Al(Nt))
        allocate(this%E21(Nt), this%E22(Nt))
        allocate(this%A21(Nt), this%A22(Nt))
        allocate(this%g1(Nt), this%g2(Nt))
        this%El = 0.0_dp
        this%Al = 0.0_dp
        this%E21 = 0.0_dp; this%E22 = 0.0_dp
        this%A21 = 0.0_dp; this%A22 = 0.0_dp
        this%g1 = 0.0_dp; this%g2 = 0.0_dp
        this%tp1 = this%tp1/(1-2/pi) ! check this
        this%pulse_offset1 = 0.0_dp
        this%pulse_offset2 = 0.0_dp

        ! Calculate amplitudes
        A01 = this%e01 / this%omega1
        A02 = this%e02 / this%omega2
        ! Calculate the envelope shapes
        ! Envelope shape for laser 1
        select case(trim(this%envelope_shape_laser1))
            case("cos2")
                do K = 1, Nt
                    time = K*dt 
                    this%g1(k) = cos2(time, this%tp1, this%t_mid1, this%pulse_offset1)
                enddo
            case("gaussian")
                do K = 1, Nt
                    time = K*dt
                    this%g1(K) = gaussian(time, this%tp1, this%t_mid1)
                enddo
            case("trapazoidal")
                do K = 1, Nt
                    time = K*dt
                    this%g1(K) = trapazoidal(time, this%tp1, this%t_mid1, this%rise_time1)
                enddo
            case default
                print*, "Laser1: Default pulse shape is CW."
        end select
        ! Envelope shape for laser 2
        select case(trim(this%envelope_shape_laser2))
            case("cos2")
                do K = 1, Nt
                    time = k*dt 
                    this%g2(k) = cos2(time, this%tp1, this%t_mid2, this%pulse_offset2)
                enddo
            case("gaussian")
                do K = 1, Nt
                    time = K*dt
                    this%g2(K) = gaussian(time, this%tp2, this%t_mid2)
                enddo
            case("trapazoidal")
                do K = 1, Nt
                    time = K*dt
                    this%g2(K) = trapazoidal(time, this%tp2, this%t_mid2, this%rise_time2)
            enddo
            case default
                print*, "Laser2: Default pulse shape is CW."
        end select

        ! Generate the electric field
        timeloop: do K = 1, Nt
            time = k*dt 
            this%E21(K) = this%E01 * this%g1(K) * cos(this%omega1 * (time - this%t_mid1 &
                  & - this%pulse_offset1) + this%phi1)   
            this%E22(K) = this%E02 * this%g2(K) * cos(this%omega2 * (time - this%t_mid2 &
                  & - this%pulse_offset2) + this%phi2)
            this%A21(k) = (-1._dp) * sum(this%E21(1:K)) * dt
            this%A22(k) = (-1._dp) * sum(this%E22(1:K)) * dt

            this%El(K) = this%E21(K) + this%E22(K)
            this%Al(K) = (-1._dp)*sum(this%El(1:K)) * dt
        enddo timeloop
        print*, "Pulse generation complete."
    end subroutine generate_pulse

    subroutine field_spectra(this)
        class(pulse_param), intent(inout) :: this
        integer :: K, void
        type(C_PTR) planTF
        complex(dp), allocatable:: E_dum(:)
        ! file tokens
        integer:: field_spec_tk

        call fftw_initialize_threads
    
        allocate(E_dum(Nt))

        call fftw_plan_with_nthreads(omp_get_max_threads())
        planTF = fftw_plan_dft_1d(Nt, E_dum, E_dum, FFTW_FORWARD, FFTW_ESTIMATE)
        print*, "Done"

        E_dum = this%El
        call fftw_execute_dft(planTF,E_dum, E_dum)
        E_dum = E_dum/sqrt(dble(Nt))

        ! write the field spectra to file
        do K = Nt/2+1, Nt
            write(field_spec_tk,*) -(Nt + 1 - K) * 2 *pi/(dt * Nt), & 
                & real(E_dum(K)), imag(E_dum(K)), abs(E_dum(K))
        enddo
        do K = 1, Nt/2
            write(field_spec_tk,*) (K-1)*2*pi/(dt*Nt), real(E_dum(K)), & 
                & imag(E_dum(K)), abs(E_dum(K))
        enddo
        close(field_spec_tk)

        call fftw_destroy_plan(planTF)
        deallocate(E_dum)
    end subroutine field_spectra

    ! A subroutine for writing the pulse to files
    subroutine write_pulse_to_file(this)
        class(pulse_param), intent(in) :: this
        integer :: K
        character(150) filename
        character(2000) mk_out_dir
        real(dp) time     
        ! file tokens
        integer:: envelope1_tk, envelope2_tk
        integer:: field1_tk, field2_tk
        integer:: elec_field_tk, vec_field_tk

        write(mk_out_dir, '(a,a)') adjustl(trim(output_data_dir)), 'pulse_data/'
        print*, "creating pulse output directory ", trim(mk_out_dir)
        call execute_command_line("mkdir -p " // adjustl(trim(mk_out_dir)))
    
        write(filename,fmt='(a,a)') adjustl(trim(mk_out_dir)), 'envelope1.out'
        open(newunit=envelope1_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a)') adjustl(trim(mk_out_dir)), 'envelope2.out'
        open(newunit=envelope2_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f4.2,a,i0,a)') adjustl(trim(mk_out_dir)), &
            & 'electric_field1_E', this%E01,'_width',Int(this%tp1*au2fs),'.out'
        open(newunit=field1_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f6.4,a,i0,a)') adjustl(trim(mk_out_dir)), &
            & 'electric_field2_E', this%E02,'_width',Int(this%tp2*au2fs),'.out'
        open(newunit=field2_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f4.2,a)') adjustl(trim(mk_out_dir)), &
            & 'Total_electric_field_phi', this%phi2/pi,'pi.out'
        open(newunit=elec_field_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f4.2,a)') adjustl(trim(mk_out_dir)), &
            & 'Total_vector_field_phi', this%phi2/pi, 'pi.out'
        open(newunit=vec_field_tk, file=filename,status="unknown")

        timeloop: do K = 1, Nt
            time = k*dt 
            write(field1_tk,*) time*au2fs, this%E21(K), this%A21(K)
            write(field2_tk,*) time*au2fs, this%E22(K), this%A22(K)
      
            write(envelope1_tk,*) time*au2fs, this%g1(K)
            write(envelope2_tk,*) time*au2fs, this%g2(K)
            write(elec_field_tk,*) time*au2fs, this%El(K)
            write(vec_field_tk,*) time*au2fs, this%Al(K)
        enddo timeloop
        close(envelope1_tk)
        close(envelope2_tk)
        close(field1_tk)
        close(field2_tk)
        close(elec_field_tk)
        close(vec_field_tk)
    end subroutine write_pulse_to_file
  
    ! pulse envelope functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cos2(time, tp, t_mid, pulse_offset)
        real(dp), intent(in) :: time, tp, t_mid, pulse_offset
        real(dp) :: cos2
        if (time .gt. (t_mid+pulse_offset-tp/2) .and. time .lt. (t_mid+pulse_offset+tp/2)) then
            cos2 = cos((time - t_mid-pulse_offset)*pi/tp)**2      
        else
            cos2 = 0._dp
        endif
    end function cos2

    function trapazoidal(time, tp, t_mid, rise_time)
        real(dp) :: time, tp, t_mid, rise_time
        real(dp) :: trapazoidal, slope, yc 
        if (time .ge. t_mid - (tp/2 + rise_time) .and. time .le. t_mid - tp/2) then
            slope = 1._dp/rise_time
            yc = (t_mid - (tp/2 + rise_time)) * slope
            trapazoidal = slope * time - yc
        elseif (time .gt. t_mid - tp/2 .and. time .le. t_mid + tp/2) then
            trapazoidal = 1._dp
        elseif (time .gt. t_mid + tp/2 .and. time .le. t_mid+(tp/2 + rise_time)) then 
            slope = -1._dp/rise_time
            yc = (t_mid + tp/2 + rise_time) * slope
            trapazoidal = slope * time - yc
        else
            trapazoidal = 0._dp
        endif
    end function trapazoidal
 
    function gaussian(time, tp, t_mid)
        implicit none
        real(dp):: time, tp, t_mid
        real(dp):: gaussian, fwhm
 
        fwhm = (4._dp * log(2._dp)) / tp**2
        gaussian = exp(-fwhm * (time - t_mid)**2)
    end function gaussian
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! A subroutine for deallocating the envelope arrays
    subroutine deallocate_envelope(this)
        class(pulse_param), intent(inout) :: this
        if (allocated(this%g1)) then
            deallocate(this%g1)
        end if
        if (allocated(this%g2)) then
            deallocate(this%g2)
        end if
    end subroutine deallocate_envelope
    ! A subroutine for deallocating the field arrays
    subroutine deallocate_field(this)
        class(pulse_param), intent(inout) :: this
        if (allocated(this%E21)) then
            deallocate(this%E21)
        end if
        if (allocated(this%E22)) then
            deallocate(this%E22)
        end if
    end subroutine deallocate_field
    ! A subroutine for deallocating all arrays
    subroutine deallocate_all(this)
        class(pulse_param), intent(inout) :: this
        call deallocate_envelope(this)
        call deallocate_field(this)
        if (allocated(this%El)) then
            deallocate(this%El)
        end if
        if (allocated(this%Al)) then
            deallocate(this%Al)
        end if
    end subroutine deallocate_all
  !------------------------------------------------------------------------------

end module pulse_mod
 

