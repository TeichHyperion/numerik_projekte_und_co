program differenzieren

implicit none
integer                 :: config_file
real(8)                 :: x1, x2, resolution, t1, t2, a_tol_exp, r_tol_exp
real(8), allocatable    :: grid(:), function_values(:), derivative(:), derivative_converged(:)
logical                 :: converged

open(newunit = config_file, file= "diff_config", status="old")

read(config_file, *) x1
read(config_file, *) x2
read(config_file, *) resolution
read(config_file, *) a_tol_exp
read(config_file, *) r_tol_exp

close(config_file)

call cpu_time(t1)


call make_grid(x_start          = x1, &
               x_end            = x2, &
               grid_resolution  = resolution, &
               grid_out         = grid) 

function_values = vector_function(grid) ! this is ~5% faster then wraping in a subroutine. Should i unwrap the others to? 

call central_difs(  y       = function_values, &
                    h       = resolution, &
                    y_prime = derivative)
                    ! 15% faster as subroutine: Time per Grid Point:    2.1624915527673716E-008 as a function 
                    !                        vs Time per Grid Point:    1.7023371002457020E-008 as subroutine

call converge_derivative(   grid_inout          = grid , &
                            resolution_current  = resolution , &
                            derivative_in       = derivative ,  &
                            absolut_tol         = 10**(a_tol_exp) , &
                            relative_tol        = 10**(r_tol_exp) , &
                            derivative_out      = derivative_converged)



call cpu_time(t2)

print*, "Time elapsed: ", (t2 - t1)
print*, "Time per Grid Point: ", (t2 - t1) / (((x2 - x1) / resolution) + 1d0 )

call quick_plot(grid, function_values, derivative_converged)





















contains 

subroutine converge_derivative(grid_inout, resolution_current, derivative_in, absolut_tol, relative_tol, derivative_out)

    implicit none

    real(8), allocatable, intent(in)    :: derivative_in(:)    
    real(8), allocatable, intent(inout) :: grid_inout(:) 
    real(8), allocatable, intent(out)   :: derivative_out(:)
    real(8), allocatable                :: derivative_working(:), derivative_old(:)
    real(8), intent(inout)              :: resolution_current
    integer                             :: j
    real(8), intent(in)                 :: absolut_tol, relative_tol
    

    converged = .false.
    j = 0
    !call move_alloc(from= derivative_in, to=derivative_working) !'from' argument of 'move_alloc' intrinsic at (1) cannot be INTENT(IN)
    derivative_working = derivative_in

    do while (.not. converged) 
        j= j+1
        if (j >= 20 .or. size(grid_inout) > 10**7) then
            print *,  "maximum number of iterations reached"
            exit
        end if

        call move_alloc(from=derivative_working, to= derivative_old)
        call double_grid_function(grid_inout, resolution_current, function_values)

        !function_values = vector_function(grid_inout)

        call central_difs(  y       = function_values, &
                            h       = resolution_current, &
                            y_prime = derivative_working)

        converged = check_convergence(derivative_old, derivative_working, absolut_tol, relative_tol)
        print*, converged
    end do

    call move_alloc(from=derivative_working, to=derivative_out) 
end subroutine converge_derivative

subroutine double_grid_function(grid_inout, grid_resolution, f_inout)
    real(8), allocatable, intent(inout)    :: grid_inout(:), f_inout(:)
    real(8), allocatable                   :: grid_old(:), f_old(:)

    !real(8), allocatable, intent(out)   :: grid_out(:)
    real(8), intent(inout)              :: grid_resolution
    real(8)                             :: x_start, x_end
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
    real(8), allocatable, intent(in)    :: y_prime_old(:), y_prime_new(:)
    real(8)                             :: combined_error, epsilon_val
    integer(8)                          :: new_size, old_size
    real(8), intent(in)                 :: abs_tol, rel_tol


    new_size = size(y_prime_new)
    old_size = size(y_prime_old)

    convergence = .false.
    combined_error = maxval( abs(y_prime_new(3:new_size-2:2) - y_prime_old(2:old_size-1) ) / (abs_tol + rel_tol * abs((y_prime_old(2:old_size-1))+ epsilon(epsilon_val) )) )
    print*, "Convergence: ", combined_error

    if (combined_error < 1d0) then
        convergence = .true.
        print*, "Final Convergence: ", combined_error
    end if
end function 


elemental function vector_function(x_in) result(f)
    implicit none

    real(8), intent(in) :: x_in!(:)
    real(8)             :: f!(size(x_in))

    f = sin(x_in**2)
     
end function vector_function

subroutine quick_plot(x, y1, y2)
    real(8), intent(in)     :: x(:), y1(:), y2(:)
    real(8), allocatable    :: y3(:)
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

    real(8), intent(in)                 :: x_start, x_end, grid_resolution
    real(8), intent(out), allocatable   :: grid_out(:)
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


