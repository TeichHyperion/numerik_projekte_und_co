
program differenzieren

implicit none
integer                 :: config_file
real(8)                 :: x1, x2, resolution, t1, t2
real(8), allocatable    :: grid(:), function_values(:), derivative(:)

open(newunit = config_file, file= "diff_config", status="old")

read(config_file, *) x1
read(config_file, *) x2
read(config_file, *) resolution

close(config_file)

call cpu_time(t1)

call make_grid(x_start          = x1, &
               x_end            = x2, &
               grid_resolution  = resolution, &
               grid_out         = grid) 

function_values = vector_function(grid)

call central_difs(  y       = function_values, &
                    h       = resolution, &
                    y_prime = derivative)

call cpu_time(t2)

print*, "Time elapsed: ", (t2 - t1)
print*, "Time per Grid Point: ", (t2 - t1) / (((x2 - x1) / resolution) + 1d0 )

call quick_plot(grid, function_values, derivative)





















contains 

logical function check_convergence(y_prime_old, y_prime_new, threshhold) result(converged)
    real(8), allocatable, intent(in)    :: y_prime_old(:), y_prime_new(:)
    real(8), intent(in)                 :: threshhold
    real(8), allocatable                :: rel_error(:)


    allocate(rel_error(size(y_prime_old)))

    rel_error = (y_prime_old -y_prime_new) / (y_prime_old)
    converged = .false.
    if (maxval(rel_error) < threshhold) then
        converged = .true.
    end if

end function 


elemental function vector_function(x_in) result(f)
    implicit none

    real(8), intent(in) :: x_in!(:)
    real(8)             :: f!(size(x_in))

    f = sin(x_in)
     
end function vector_function

subroutine quick_plot(x, y1, y2)
    real(8), intent(in)     :: x(:), y1(:), y2(:)
    real(8), allocatable    :: y3(:)
    integer                 :: i, unit_dn

    allocate(y3(size(x)))

    ! 1. Write data to a temporary file
    open(newunit=unit_dn, file='plot_data.dat', status='replace')
    y3= cos(x)

    do i = 1, size(x)
        write(unit_dn, *) x(i), y1(i), y2(i), y3(i)
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

    real(8), intent(in)                 :: x_start, x_end, grid_resolution
    real(8), intent(inout), allocatable :: grid_out(:)
    integer(8)                          :: i, number_of_points


    number_of_points = nint((x_end - x_start) / grid_resolution) + 1 !+1 includes x2
    print * , "grid points: ", number_of_points

    allocate(grid_out(number_of_points))

    do i = 1, size(grid_out)
        grid_out(i) = x_start + grid_resolution*(i-1)
    end do 

end subroutine make_grid

subroutine forward_difs(x, y, h, y_prime)
    implicit none

    real(8), intent(in)                 :: x(:), y(:) 
    real(8), allocatable, intent(out)   :: y_prime(:)
    real(8)                             :: h
    integer(8)                          :: n_elements
    
    allocate(y_prime(size(x)))

    n_elements = size(x)


    y_prime(1:n_elements-1) = (y(2:n_elements) - y(1:n_elements-1)) / h

    !fix last element with backward diff
    y_prime(n_elements) = (y(n_elements) - y(n_elements-1)) / h

end subroutine forward_difs

subroutine central_difs(y, h, y_prime)
    implicit none

    real(8), intent(in)                 :: y(:) 
    real(8), allocatable, intent(out)   :: y_prime(:)
    real(8)                             :: h
    integer(8)                          :: n_elements
    
    allocate(y_prime(size(y)))

    n_elements = size(y)


    ! first element to calculate is i=2. so we need i= 3 and i=1. Last one is i = n-1 so we need i=n and i=n-2
    y_prime(2:n_elements-1) = (y(3:n_elements) - y(:n_elements-2)) / (2d0*h)

    !fix last element with backward diff
    y_prime(n_elements) = (y(n_elements) - y(n_elements-1)) / h

    !fix first element with forward diff
    y_prime(1) = (y(2) - y(1)) / h
    

end subroutine central_difs


end program differenzieren


