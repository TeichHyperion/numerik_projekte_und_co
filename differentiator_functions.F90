module diff_solver_mod

    use iso_fortran_env, only: dp => real64
    implicit none

contains 

subroutine converge_derivative(grid_inout, function_values_inout, number_grid_points, derivative_in, absolut_tol, relative_tol, derivative_out, grid_type_sub)

    implicit none

    real(dp), allocatable, intent(in)   :: derivative_in(:)    
    real(dp), allocatable, intent(inout):: grid_inout(:), function_values_inout(:)
    real(dp), allocatable, intent(out)  :: derivative_out(:)
    real(dp), allocatable               :: derivative_working(:), derivative_old(:)
    integer(8), intent(inout)           :: number_grid_points
    integer                             :: j
    integer, intent(in)                 :: grid_type_sub
    real(dp), intent(in)                :: absolut_tol, relative_tol
    logical                             :: converged
    

    converged = .false.
    j = 0

    derivative_working = derivative_in


    do while (.not. converged) 
        j= j+1
        if (j >= 22 .or. size(grid_inout) > 10**6) then
            !print *,  "maximum number of iterations reached"
            exit
        end if

        call move_alloc(from    = derivative_working, &
                        to      = derivative_old)

        call double_grid_function(  grid_inout          = grid_inout, &
                                    number_grid_points  = number_grid_points, &
                                    f_inout             = function_values_inout,&
                                    grid_type_sub       = grid_type_sub) ! there is no noticable benefit moving the double_function back out into a function  
        
                                  
        call central_difs_4th_uniform(  y       = function_values_inout, &
                            grid_in = grid_inout, &
                            y_prime = derivative_working)

        converged = check_convergence(y_prime_old   = derivative_old,&
                                      y_prime_new   = derivative_working, &
                                      abs_tol       = absolut_tol,&
                                      rel_tol       = relative_tol)
    end do

    call move_alloc(from=derivative_working, to=derivative_out) 
end subroutine converge_derivative

subroutine double_grid_function(grid_inout, number_grid_points, f_inout, grid_type_sub)
    real(dp), allocatable, intent(inout)    :: grid_inout(:), f_inout(:)
    real(dp), allocatable                   :: grid_old(:), f_old(:), delta(:)
    integer, intent(in)                  :: grid_type_sub

    !real(dp), allocatable, intent(out)     :: grid_out(:)
    integer(8), intent(inout)               :: number_grid_points
    real(dp)                                :: x_start, x_end, safty
    integer(8)                              :: n_new, n_old

    number_grid_points = number_grid_points * 2 -1 ! doubles the resolution not the points

    call move_alloc(from=grid_inout,    to=grid_old)
    call move_alloc(from=f_inout,       to=f_old)

    x_start = grid_old(1)
    x_end   = grid_old(size(grid_old))

    n_old = size(grid_old) 
    n_new = number_grid_points

    !print * , "grid points new: ", n_new

    allocate(grid_inout(n_new))
    allocate(f_inout(n_new))

    grid_inout(1: n_new: 2) =  grid_old(1 : n_old)

    select case(grid_type_sub)
    case(1)
        grid_inout(2: n_new: 2) = (grid_old(2 : n_old) + grid_old(1 : n_old-1) ) / 2d0
    case(2)
        allocate(delta(n_old-1))
        call random_number(delta)
        !delta = delta + tiny(safty) 
        delta = 0.1_dp + delta * 0.8_dp
        grid_inout(2: n_new: 2) = grid_old(1:n_old-1) + ( delta * (grid_old(2:n_old) - grid_old(1:n_old-1)))
    end select


    f_inout(1:n_new:2) = f_old(:)
    f_inout(2:n_new:2) = vector_function(grid_inout(2: n_new: 2)) 

end subroutine


logical function check_convergence(y_prime_old, y_prime_new, abs_tol, rel_tol) result(convergence)
    real(dp), allocatable, intent(in)    :: y_prime_old(:), y_prime_new(:)
    real(dp)                             :: combined_error, max_delta
    integer(8)                           :: new_size, old_size
    real(dp), intent(in)                 :: abs_tol, rel_tol


    new_size = size(y_prime_new) 
    old_size = size(y_prime_old) 

    convergence = .false.

    !old_vals  =  y_prime_old(3:old_size-2)   ![1,      2,      3,       4,      5]
    !new_vals  =  y_prime_new(5:new_size-4:2) ![1, 1.5, 2, 2.5, 3, 3.5,  4, 4.5, 5]

    
    max_delta = maxval(abs(y_prime_new(5:new_size-4:2) - y_prime_old(3:old_size-2)  ))

    combined_error = maxval( abs(y_prime_new(5:new_size-4:2) - y_prime_old(3:old_size-2) )  &
                            /(abs_tol + rel_tol * abs(y_prime_old(3:old_size-2) )))

   ! print*, "Convergence: ", combined_error
   ! print*, "Max difference", max_delta
    if (combined_error < 1d0) then
        convergence = .true.
   !     print*, "Final Convergence: ", combined_error
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
    integer                  :: j, unit_dn

    !allocate(y3(size(x)))

    ! 1. Write data to a temporary file
    open(newunit=unit_dn, file='plot_data.dat', status='replace')
    y3 = 3*x**2 !* (x**2)

    do j = 1, size(x)
        write(unit_dn, *) x(j), y1(j), y2(j), y3(j)
    end do
    close(unit_dn)

    call execute_command_line("gnuplot -p -e &
        & ""set title 'Function and Derivative'; &
        & set grid; &
        & plot 'plot_data.dat' using 1:2 w l title 'f(x)', &
        & '' using 1:3 w l title 'f''(x)', &
        & '' using 1:4 w l title 'f''(x) analytical' "" ")
end subroutine quick_plot

subroutine make_grid(x_start, x_end, number_grid_points, grid_type_in, grid_out)
    real(dp), intent(in)                :: x_start, x_end
    integer(8), intent(in)              :: number_grid_points
    integer, intent(in)                 :: grid_type_in
    real(dp), allocatable, intent(out)  :: grid_out(:)

    select case((grid_type_in))
    case(1)
        call make_grid_uniform( x_start             = x_start, &
                                x_end               = x_end, &
                                number_grid_points  = number_grid_points, &
                                grid_out            = grid_out)
    case(2)
        call make_grid_random(  x_start             = x_start, &
                                x_end               = x_end, &
                                number_grid_points  = number_grid_points, &
                                grid_out            = grid_out)
    end select

end subroutine make_grid

subroutine make_grid_uniform(x_start, x_end, number_grid_points, grid_out)
    implicit none

    real(dp), intent(in)                :: x_start, x_end
    integer(8), intent(in)              :: number_grid_points
    real(dp)                            :: grid_resolution
    real(dp), allocatable, intent(out)  :: grid_out(:)
    integer(8)                          :: j


    !number_of_points = nint((x_end - x_start) / grid_resolution) + 1 ! +1 includes x2
    !print * , "grid points: ", number_grid_points
    !(number_of_points /nint((x_end - x_start))  - 1 = (1 / grid_resolution) 
    grid_resolution  = (x_end - x_start) / (number_grid_points-1)


    !print*, "resolution: ", grid_resolution
    allocate(grid_out(number_grid_points))

    do j = 1, size(grid_out)
        grid_out(j) = x_start + grid_resolution*(j-1)
      !  print*, j
    end do 

end subroutine make_grid_uniform

subroutine make_grid_random(x_start, x_end, number_grid_points, grid_out)
    implicit none

    real(dp), intent(in)                :: x_start, x_end
    integer(8), intent(in)              :: number_grid_points
    real(dp),allocatable, intent(out)   :: grid_out(:)
    real(dp), allocatable               :: grid_tmp(:)
    integer(8)                          :: j, k
    real(dp)                            :: temp

    allocate(grid_out(number_grid_points))
    allocate(grid_tmp(number_grid_points))

    !call random_number(grid_tmp)
    !grid_out = (x_end - x_start)  * grid_tmp + x_start
    grid_tmp(1) = x_start

    do j = 2, number_grid_points
        call random_number(temp)
        grid_tmp(j) = grid_tmp(j-1) + temp 
    end do

    grid_out = ((grid_tmp - grid_tmp(1) )/ (grid_tmp(number_grid_points)- grid_tmp(1)) )* (x_end- x_start) + grid_tmp(1)

    !print*, grid_out

end subroutine make_grid_random

subroutine forward_difs(x, y, h, y_prime)
    implicit none

    real(dp), intent(in)                 :: x(:), y(:) 
    real(dp), allocatable, intent(out)   :: y_prime(:)
    real(dp)                             :: h
    integer(8)                           :: n_elements
    
    allocate(y_prime(size(x)))

    n_elements = size(x)


    y_prime(1:n_elements-1) = (y(2:n_elements) - y(1:n_elements-1)) / h

    !fix last element with backward diff
    y_prime(n_elements) = (y(n_elements) - y(n_elements-1)) / h

end subroutine forward_difs

subroutine central_difs(y, grid_in, y_prime)
    implicit none

    real(dp), intent(in)                 :: y(:) 
    real(dp), allocatable, intent(out)   :: y_prime(:)
    real(dp)                             :: h
    real(dp), intent(in)                 :: grid_in(:)
    !real(dp), allocatable                :: dx(:)
    integer(8)                           :: n_elements
    

    n_elements = size(y)
    allocate(y_prime(n_elements))
    !dx = grid_in(2:) - grid_in(1:size(grid_in)-1)

    ! first element to calculate is i=2. so we need i= 3 and i=1. Last one is i = n-1 so we need i=n and i=n-2
    y_prime(2:n_elements-1) = (y(3:n_elements) - y(:n_elements-2)) /(grid_in(3:n_elements) - grid_in(:n_elements-2))

    !fix last element with backward diff
    y_prime(n_elements) = (y(n_elements) - y(n_elements-1)) / (grid_in(n_elements) - grid_in(n_elements-1))

    !fix first element with forward diff
    y_prime(1) = (y(2) - y(1)) / (grid_in(2) - grid_in(1))
    

end subroutine central_difs

subroutine central_difs_4th_uniform(y, grid_in, y_prime)
    implicit none

    real(dp), intent(in)                 :: y(:) 
    real(dp), allocatable, intent(out)   :: y_prime(:)
    real(dp)                             :: h
    real(dp), intent(in)                 :: grid_in(:)
    !real(dp), allocatable                :: y1(:), y2(:), y_1(:), y_2(:)
    integer(8)                           :: n_elements
    

    n_elements = size(y)
    allocate(y_prime(n_elements))
    !dx = grid_in(2:) - grid_in(1:size(grid_in)-1)
    h   =  grid_in(2) - grid_in(1)


   ! y2  = - y(5:n_elements)
   ! y1  = 8*y(4:n_elements-1)
   ! y_1 =-8*y(2:n_elements-3)
   ! y_2 =   y(1:n_elements-4)

    ! first element to calculate is i=3. so we need i= 4, 5 and i=1, 2. Last one is i = n-2 so we need i=n, n-1 and i= n-3, n-4
    y_prime(3:n_elements-2) =  (-     y(5:n_elements)   + &
                                  8 * y(4:n_elements-1) + &
                                - 8 * y(2:n_elements-3) + &
                                      y(1:n_elements-4))/ &
                                       (12d0*h)

    !fix last element with backward diff
    y_prime(n_elements)   = (y(n_elements) - y(n_elements-1)) / (1*h)
    y_prime(n_elements-1) = (y(n_elements) - y(n_elements-2)) / (2*h)


    !fix first element with forward diff
    y_prime(1) = (y(2) - y(1)) / (1*h)
    y_prime(2)   = (y(3) - y(1)) / (2.0_dp * h)
    

end subroutine central_difs_4th_uniform

subroutine central_difs_rand(y, x, y_prime)
    real(dp), intent(in)               :: y(:), x(:) 
    real(dp), allocatable, intent(out) :: y_prime(:)
    real(dp)                           :: h1, h2
    integer(8)                         :: n, i

    n = size(y)
    allocate(y_prime(n))

    ! Interior points: Weighted Non-Uniform Central Difference
    do i = 2, n-1
        h1 = x(i) - x(i-1)
        h2 = x(i+1) - x(i)
        
        ! 2nd order formula for non-uniform grids
        y_prime(i) = (h1**2 * y(i+1) - h2**2 * y(i-1) + (h2**2 - h1**2) * y(i)) / &
                     (h1 * h2 * (h1 + h2))
    end do

    ! Boundaries: 1st order Forward/Backward
    y_prime(1) = (y(2) - y(1)) / (x(2) - x(1))
    y_prime(n) = (y(n) - y(n-1)) / (x(n) - x(n-1))
end subroutine central_difs_rand

end module diff_solver_mod