program differenzieren

use iso_fortran_env, only: dp => real64
use diff_solver_mod

implicit none
integer                 :: config_file
real(dp)                :: x1, x2, a_tol_exp, r_tol_exp
real(dp), allocatable   :: grid(:), function_values(:), derivative(:), derivative_converged(:)
integer(8)              :: tick_start, tick_end, tick_rate, config_io, n_grid_points
integer                 :: grid_type

open(newunit = config_file, file= "diff_config", status="old", iostat=config_io)

read(config_file, *) x1
read(config_file, *) x2
read(config_file, *) grid_type
read(config_file, *) n_grid_points
read(config_file, *) a_tol_exp
read(config_file, *) r_tol_exp

close(config_file)

call system_clock(count_rate=tick_rate)
call system_clock(count=tick_start)

call make_grid( x_start             = x1, &
                x_end               = x2, &
                number_grid_points  = n_grid_points, &
                grid_type_in        = grid_type, &
                grid_out            = grid)

print*, grid(1)
print*, grid(size(grid))

function_values = vector_function(grid) ! this is ~5% faster then wraping in a subroutine. Should i unwrap the others to? no



call central_difs_4th_uniform(  y       = function_values, &
                                grid_in = grid ,&
                        ! h       = resolution, &
                                y_prime = derivative)

!call central_difs(  y       = function_values, &
!                    grid_in = grid ,&
!                   ! h       = resolution, &
!                    y_prime = derivative)
                    ! 15% faster as subroutine: Time per Grid Point:    2.1624915527673716E-008 as a function 
                    !                        vs Time per Grid Point:    1.7023371002457020E-008 as subroutine

!call quick_plot(grid, function_values, derivative)



call converge_derivative(   grid_inout           = grid , &
                            function_values_inout= function_values ,&
                            number_grid_points   = n_grid_points , &
                            derivative_in        = derivative ,  &
                            absolut_tol          = 10d0**(a_tol_exp) , &
                            relative_tol         = 10d0**(r_tol_exp) , &
                            derivative_out       = derivative_converged, &
                            grid_type_sub        = grid_type)



call system_clock(count=tick_end)

print*, "Time elapsed: ", (tick_end - tick_start) / real(tick_rate, 8)
print*, "Time per Grid Point: ", (tick_end - tick_start) / n_grid_points / real(tick_rate, 8)

call quick_plot(grid, function_values, derivative_converged)




end program differenzieren


