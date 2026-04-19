program benchmark_diff
    use iso_fortran_env, only: dp => real64
    use diff_solver_mod

    implicit none

    integer(8)            :: t1, t2, rate
    integer               :: i,j , n_trials = 1000
    real(dp)              :: h
    real(dp), allocatable :: x(:), y(:), yp(:), yp_converged(:)
    real(dp), allocatable :: times_math(:), times_conv(:)
    real(dp)              :: a_tol = 1.0d-4, r_tol = 1.0d-4
    integer(8)            :: np = 101

    ! Allocate timing arrays
    allocate(times_math(n_trials))
    allocate(times_conv(n_trials))

    call system_clock(count_rate=rate)

    ! =========================================================================
    ! BENCHMARK 1: Pure Math Routine (100,000 points)
    ! =========================================================================
    h = 0.001_dp
    allocate(x(100000), y(100000))
    x = [(real(i, dp) * h, i=1, 100000)]
    y = sin(x**2)

    ! Warm up
    call central_difs_4th_uniform(y, x, yp)

    do i = 1, n_trials
        if (allocated(yp)) deallocate(yp) ! Do not time the deallocation
        
        call system_clock(count=t1)
        call central_difs_4th_uniform(y, x, yp)
        call system_clock(count=t2)
        
        times_math(i) = real(t2 - t1, dp) / real(rate, dp)
    end do
    
    deallocate(x, y, yp)

    ! =========================================================================
    ! BENCHMARK 2: Full Convergence Routine
    ! =========================================================================
    ! Warm up loop
    do i = 1, n_trials
        ! 1. RESET STATE (Outside the timer!)
        h = 0.1_dp
        np =101
        allocate(x(np))
        x = [(-5.0_dp + real(j-1, dp)*h, j=1, 101)]
        y = sin(x**2)
        call central_difs_4th_uniform(y, x, yp)
        
        ! 2. START TIMING
        call system_clock(count=t1)
        


        call converge_derivative(   grid_inout           = x , &
                            function_values_inout= y ,&
                            number_grid_points   = np , &
                            derivative_in        = yp ,  &
                            absolut_tol          = a_tol , &
                            relative_tol         = r_tol , &
                            derivative_out       = yp_converged, &
                            grid_type_sub        = 1)
                                  
        call system_clock(count=t2)
        ! 3. END TIMING

        times_conv(i) = real(t2 - t1, dp) / real(rate, dp)
        
        ! Cleanup for the next trial
        deallocate(x, y, yp, yp_converged)
    end do

    ! =========================================================================
    ! PRINT RESULTS
    ! =========================================================================
    print *, "=================================================="
    print *, "              BENCHMARK RESULTS                   "
    print *, "=================================================="
    print *, "Trials per test: ", n_trials
    print *, ""
    
    call print_stats("1. Math Routine (100k points)", times_math)
    call print_stats("2. Full Convergence (-5 to 5)  ", times_conv)

contains

    ! Helper subroutine to calculate mean and sample standard deviation
    subroutine print_stats(test_name, times)
        character(len=*), intent(in) :: test_name
        real(dp), intent(in)         :: times(:)
        real(dp)                     :: mean_time, std_dev, variance
        integer                      :: n

        n = size(times)
        
        ! Calculate Mean
        mean_time = sum(times) / real(n, dp)
        
        ! Calculate Sample Standard Deviation
        variance = sum((times - mean_time)**2) / real(n - 1, dp)
        std_dev = sqrt(variance)

        ! Convert to milliseconds for readability
        mean_time = mean_time * 1000.0_dp
        std_dev   = std_dev * 1000.0_dp

        print "(A)", test_name
        print "(A, F10.4, A)", "  Mean Time: ", mean_time, " ms"
        print "(A, F10.4, A)", "  Std Dev:   ", std_dev,   " ms"
        print *, ""
    end subroutine print_stats

end program benchmark_diff