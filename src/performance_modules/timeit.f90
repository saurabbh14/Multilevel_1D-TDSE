module timeit 
    use VarPrecision, only: dp
    implicit none
    private
    public:: timer, start, stop, elapsed

    type :: timer
        integer :: sys_scount, sys_ecount        ! system clock ticks (start/end)
        integer :: rate                     ! clock ticks per second
        real(dp) :: cpu_st, cpu_ft, timer_elapsed_time  ! cpu_time start/end and elapsed
    contains
        procedure :: start
        procedure :: stop
        procedure :: elapsed
    end type timer

contains
    subroutine start(this)
        class(timer), intent(inout) :: this
        call cpu_time(this%cpu_st)
        call system_clock(this%sys_scount, this%rate)
    end subroutine start

    subroutine stop(this)
        class(timer), intent(inout) :: this
        call cpu_time(this%cpu_ft)
        call system_clock(this%sys_ecount)
    end subroutine stop

    subroutine elapsed(this)
        class(timer), intent(inout) :: this
        this%timer_elapsed_time = this%cpu_ft - this%cpu_st
        write(*,*) "CPU time elapsed: ", this%timer_elapsed_time, " seconds"
        write(*,*) "System clock time elapsed: ", real(this%sys_ecount - this%sys_scount, dp) / real(this%rate, dp), " seconds"
    end subroutine elapsed

end module timeit