module PrintInputVars
    use global_vars
    use pot_param
    implicit none
    contains
      subroutine print_input_vars()
        print*, "Number of grid points: NR =", NR
        print*, "Masses: m1 =", m1, "m2 =", m2
        print*, "Time grid: dt =", dt, "fs"
        print*, "Number of time steps: Nt =", Nt, "steps"
        print*, "Number of electronic states: Nstates =", Nstates
        print*, "Electronic potential kind: Elec_pot_kind =", trim(Elec_pot_kind)
        print*, "Number of maximum considered vibrational states: guess_vstates =", guess_vstates
        print*
        print*, "Guess vibrational wavefunction (Gaussian): Initial position (RI) =", RI
        print*, "with initial width (kappa) =", kappa
        print*
        print*, "Laser parameters:"
        print*, "Laser #1:"
        print*, "Envelope shape:", trim(envelope_shape_laser1)
        print*, "Lambda:", lambda1, "nm"
        print*, "Electric field strength:", E01, "a.u."
        print*, "Pulse width (tp):", tp1, "fs"
        print*, "Pulse midpoint:", t_mid1, "fs"
        print*, "phi1:", phi1, "pi"
        print*, "Rise time:", rise_time1, "fs"
        print*, "Laser #2:"
        print*, "Envelope shape:", trim(envelope_shape_laser2)
        print*, "Lambda:", lambda2, "nm"
        print*, "Electric field strength:", E02, "a.u."
        print*, "Pulse width (tp):", tp2, "fs"
        print*, "Pulse midpoint:", t_mid2, "fs"
        print*, "phi2:", phi2, "pi"
        print*, "Rise time:", rise_time2, "fs"
        print*
        print*, "Input and Output Directories:"
        print*, "Input data directory:", trim(input_data_dir)
        print*, "Output data directory:", trim(output_data_dir)
        print*, "Transition Dipole switched off:", total_trans_off
        print*
        print*, "Absorber function: ", absorber
        print*
        print*, "TDSE Initial State:"
        print*, "Mode:", trim(initial_distribution)
        print*, "electronic state(s)", (N_ini-1)
        print*, "vibrational state(s)", (v_ini-1)
        print*, "Gaussian Distribution TDSE:"
        print*, "centered at RI: ", RI_tdse
        print*, "standard deviation: ", kappa_tdse
        print*
        print*, "FFTW Parallelization:"
        print*, "TDSE Propagation FFTW: ", trim(prop_par_FFTW)
        print*, "ITP FFTW: ", trim(ITP_par_FFTW)
        print*
        print*, "_________________________"
        print*
        print*, "Final Parameters"
        print*, "_________________________"
        print*
        print*, "dt = ", SNGL(dt), "a.u."
        print*, "dR = ", SNGL(dR), "a.u."
        print*, "dPR = ", SNGL(dpR), "a.u."
        print*, "RI=", sngl(RI), "a.u."
        print*, "R0=", sngl(R0), "a.u.", "Rend=",sngl(Rend), "a.u."
        print*, "Wavelength 1 =", sngl(lambda1), "nm"
        print*, "Phase 1 =", sngl(phi1)
        print*, "Field strength =", sngl(e01), "a.u.", sngl(e01*e02au), "V/m"
        print*, "Intensity =", sngl(e01**2*3.509e16_dp), "W/cm2"
        print*, "Wavelength 2 =", sngl(lambda2), "nm"
        print*, "Phase 2 =", sngl(phi2)
        print*, "Field strength =", sngl(e02), "a.u.", sngl(e02*e02au), "V/m"
        print*, "Intensity =", sngl(e02**2*3.509e16_dp), "W/cm2"
        print*, "Wave duration =", sngl(tp1*au2fs), "fs"
        print*
        print*, "kap =", kap
        print*, "lam =", lam
        print*
        print*, "__________________________"
        print*
  
      end subroutine print_input_vars
  end module PrintInputVars