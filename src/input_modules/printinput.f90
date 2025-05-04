module PrintInputVars
    use global_vars
    use pot_param
    use pulse_mod
   
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
        print*, "Final grid Parameters"
        print*, "_________________________"
        print*
        print*, "dt = ", SNGL(dt), "a.u."
        print*, "dR = ", SNGL(dR), "a.u."
        print*, "dPR = ", SNGL(dpR), "a.u."
        print*, "RI=", sngl(RI), "a.u."
        print*, "R0=", sngl(R0), "a.u.", "Rend=",sngl(Rend), "a.u."
        print*
        print*, "kap =", kap
        print*, "lam =", lam
        print*
        print*, "__________________________"
        print*
  
      end subroutine print_input_vars
  end module PrintInputVars