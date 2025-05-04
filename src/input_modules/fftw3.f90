module FFTW3
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
  !  include '/usr/include/fftw3.f03'                                        ! Desktop packet
  !  include '/home/me23jok/ProjectX/FFTW3/include/fftw3.f03' ! ARA cluster
  !  include '/usr/local/include/fftw3.f03'
  !  include '/scratch/Saurabh/FFTW3/install/include/fftw3.f03'  ! nias
  !  include '/home/me23jok/fftw-3.3.10/include/fftw3.f03' ! draco
  end module