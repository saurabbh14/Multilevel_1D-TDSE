module output_dir_mod
    use global_vars, only: output_data_dir, pulse_data_dir, &
        & nucl_wf_dir, time_prop_dir
    implicit none

contains
    subroutine output_dir_check
        use global_vars, only: output_data_dir
        character(2000):: mk_out_dir
        ! Ensure output directory exists (create if missing)
        write(mk_out_dir, '(a)') adjustl(trim(output_data_dir))
        print*, "creating output directory ", trim(mk_out_dir)
        call execute_command_line("mkdir -p " // adjustl(trim(mk_out_dir)))
    end subroutine output_dir_check

    subroutine setup_output_dir()
        ! Output directory for pulse data
        write(pulse_data_dir, '(a,a)') adjustl(trim(output_data_dir)), 'pulse_data/'
        print*, "checking/creating time propagation output directory ", trim(pulse_data_dir)
        call execute_command_line("mkdir -p " // adjustl(trim(pulse_data_dir)))

        ! Output directory for nuclear wavepacket data
        write(nucl_wf_dir, '(a,a)') adjustl(trim(output_data_dir)), 'nuclear_wavepacket_data/'
        print*, "checking/creating time propagation output directory ", trim(nucl_wf_dir)
        call execute_command_line("mkdir -p " // adjustl(trim(nucl_wf_dir)))

        ! Output directory for pulse data
        write(time_prop_dir, '(a,a)') adjustl(trim(output_data_dir)), 'time_prop/'
        print*, "checking/creating time propagation output directory ", trim(time_prop_dir)
        call execute_command_line("mkdir -p " // adjustl(trim(time_prop_dir)))        

    end subroutine setup_output_dir

end module output_dir_mod