program differenzieren

use iso_fortran_env, only: dp => real64

implicit none
integer                 :: config_file
real(dp)                :: x1, x2, resolution, a_tol_exp, r_tol_exp
real(dp), allocatable   :: grid(:), function_values(:), derivative(:), derivative_converged(:)
integer(8)              :: tick_start, tick_end, tick_rate, config_io

open(newunit = config_file, file= "diff_config", status="old", iostat=config_io)

read(config_file, *) x1
read(config_file, *) x2
read(config_file, *) resolution
read(config_file, *) a_tol_exp
read(config_file, *) r_tol_exp

close(config_file)

call system_clock(count_rate=tick_rate)
call system_clock(count=tick_start)

call make_grid(x_start          = x1, &
               x_end            = x2, &
               grid_resolution  = resolution, &
               grid_out         = grid) 

function_values = vector_function(grid) ! this is ~5% faster then wraping in a subroutine. Should i unwrap the others to? no

call central_difs(  y       = function_values, &
                    h       = resolution, &
                    y_prime = derivative)
                    ! 15% faster as subroutine: Time per Grid Point:    2.1624915527673716E-008 as a function 
                    !                        vs Time per Grid Point:    1.7023371002457020E-008 as subroutine

call converge_derivative(   grid_inout           = grid , &
                            function_values_inout= function_values ,&
                            resolution_current   = resolution , &
                            derivative_in        = derivative ,  &
                            absolut_tol          = 10d0**(a_tol_exp) , &
                            relative_tol         = 10d0**(r_tol_exp) , &
                            derivative_out       = derivative_converged)



call system_clock(count=tick_end)

print*, "Time elapsed: ", (tick_end - tick_start) / real(tick_rate, 8)
print*, "Time per Grid Point: ", (tick_end - tick_start) / (((x2 - x1) / resolution) + 1d0 ) / real(tick_rate, 8)

call quick_plot(grid, function_values, derivative_converged)





















contains 

subroutine converge_derivative(grid_inout, function_values_inout, resolution_current, derivative_in, absolut_tol, relative_tol, derivative_out)

    implicit none

    real(dp), allocatable, intent(in)    :: derivative_in(:)    
    real(dp), allocatable, intent(inout) :: grid_inout(:), function_values_inout(:)
    real(dp), allocatable, intent(out)   :: derivative_out(:)
    real(dp), allocatable                :: derivative_working(:), derivative_old(:)
    real(dp), intent(inout)              :: resolution_current
    integer                              :: j
    real(dp), intent(in)                 :: absolut_tol, relative_tol
    logical                              :: converged
    

    converged = .false.
    j = 0

    derivative_working = derivative_in

    do while (.not. converged) 
        j= j+1
        if (j >= 20 .or. size(grid_inout) > 10**7) then
            !print *,  "maximum number of iterations reached"
            exit
        end if

        call move_alloc(from=derivative_working, to= derivative_old)

        call double_grid_function(grid_inout, resolution_current, function_values_inout) ! there is no noticable benefit moving the double_function back out into a function  

        call central_difs(  y       = function_values_inout, &
                            h       = resolution_current, &
                            y_prime = derivative_working)

        converged = check_convergence(derivative_old, derivative_working, absolut_tol, relative_tol)
    end do

    call move_alloc(from=derivative_working, to=derivative_out) 
end subroutine converge_derivative

subroutine double_grid_function(grid_inout, grid_resolution, f_inout)
    real(dp), allocatable, intent(inout)    :: grid_inout(:), f_inout(:)
    real(dp), allocatable                   :: grid_old(:), f_old(:)

    !real(dp), allocatable, intent(out)   :: grid_out(:)
    real(dp), intent(inout)              :: grid_resolution
    real(dp)                             :: x_start, x_end
    integer(8)                          :: n_new, n_old

    grid_resolution = grid_resolution * 0.5d0

    call move_alloc(from=grid_inout,    to=grid_old)
    call move_alloc(from=f_inout,       to=f_old)

    x_start = grid_old(1)
    x_end   = grid_old(size(grid_old))

    n_old = size(grid_old)
    n_new = size(grid_old) * 2 - 1

    print * , "grid points: ", n_new

    allocate(grid_inout(n_new))
    allocate(f_inout(n_new))

    grid_inout(1: n_new: 2) =  grid_old(1 : n_old)
    grid_inout(2: n_new: 2) = (grid_old(2 : n_old) + grid_old(1 : n_old-1) ) / 2d0

    f_inout(1:n_new:2) = f_old(:)
    f_inout(2:n_new:2) = vector_function(grid_inout(2: n_new: 2))

end subroutine


logical function check_convergence(y_prime_old, y_prime_new, abs_tol, rel_tol) result(convergence)
    real(dp), allocatable, intent(in)    :: y_prime_old(:), y_prime_new(:)
    real(dp)                             :: combined_error, epsilon_val
    integer(8)                          :: new_size, old_size
    real(dp), intent(in)                 :: abs_tol, rel_tol


    new_size = size(y_prime_new)
    old_size = size(y_prime_old)

    convergence = .false.
    combined_error = maxval( abs(y_prime_new(3:new_size-2:2) - y_prime_old(2:old_size-1) ) / (abs_tol + rel_tol * abs((y_prime_old(2:old_size-1)) )) )
    print*, "Convergence: ", combined_error

    if (combined_error < 1d0) then
        convergence = .true.
        print*, "Final Convergence: ", combined_error
    end if
end function 


elemental function vector_function(x_in) result(f)
    implicit none

    real(dp), intent(in) :: x_in!(:)
    real(dp)             :: f!(size(x_in))

    f = sin(x_in**2)
     
end function vector_function

subroutine quick_plot(x, y1, y2)
    real(dp), intent(in)     :: x(:), y1(:), y2(:)
    real(dp), allocatable    :: y3(:)
    integer                 :: j, unit_dn

    !allocate(y3(size(x)))

    ! 1. Write data to a temporary file
    open(newunit=unit_dn, file='plot_data.dat', status='replace')
    y3 = 2*x *cos (x**2)

    do j = 1, size(x)
        write(unit_dn, *) x(j), y1(j), y2(j), y3(j)
    end do
    close(unit_dn)

    call execute_command_line("gnuplot -p -e &
        & ""set title 'Function and Derivative'; &
        & set grid; &
        & plot 'plot_data.dat' using 1:2 with lines title 'f(x)', &
        & '' using 1:3 with lines title 'f''(x)', &
        & '' using 1:4 with lines title 'f''(x) analytical' "" ")
end subroutine quick_plot

subroutine make_grid(x_start, x_end, grid_resolution, grid_out)
    implicit none

    real(dp), intent(in)                 :: x_start, x_end, grid_resolution
    real(dp), intent(out), allocatable   :: grid_out(:)
    integer(8)                          :: j, number_of_points


    number_of_points = nint((x_end - x_start) / grid_resolution) + 1 !+1 includes x2
    print * , "grid points: ", number_of_points

    allocate(grid_out(number_of_points))

    do j = 1, size(grid_out)
        grid_out(j) = x_start + grid_resolution*(j-1)
    end do 

end subroutine make_grid

subroutine forward_difs(x, y, h, y_prime)
    implicit none

    real(dp), intent(in)                 :: x(:), y(:) 
    real(dp), allocatable, intent(out)   :: y_prime(:)
    real(dp)                             :: h
    integer(8)                          :: n_elements
    
    allocate(y_prime(size(x)))

    n_elements = size(x)


    y_prime(1:n_elements-1) = (y(2:n_elements) - y(1:n_elements-1)) / h

    !fix last element with backward diff
    y_prime(n_elements) = (y(n_elements) - y(n_elements-1)) / h

end subroutine forward_difs

subroutine central_difs(y, h, y_prime)
    implicit none

    real(dp), intent(in)                 :: y(:) 
    real(dp), allocatable, intent(out)   :: y_prime(:)
    real(dp)                             :: h, inv_2h
    integer(8)                           :: n_elements
    

    n_elements = size(y)
    allocate(y_prime(n_elements))


    ! first element to calculate is i=2. so we need i= 3 and i=1. Last one is i = n-1 so we need i=n and i=n-2
    y_prime(2:n_elements-1) = (y(3:n_elements) - y(:n_elements-2)) / (2d0*h)

    !fix last element with backward diff
    y_prime(n_elements) = (y(n_elements) - y(n_elements-1)) / h

    !fix first element with forward diff
    y_prime(1) = (y(2) - y(1)) / h
    

end subroutine central_difs


end program differenzieren


