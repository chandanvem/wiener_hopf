PROGRAM main

  USE bessel_utils  
  USE io_utils
  USE contour_generate_utils
  USE contour_init_utils
  USE omp_lib
  USE input_params
  USE kernel_integral_utils
  USE fplus_utils
  USE IFT_integral_utils
  USE user_defined_precompute
  USE user_defined_functions 
  USE user_defined_IFT 
  USE user_defined_fplus

  IMPLICIT none

  type(input_params_t)   :: input_data
  type(contour_params_t) :: contour_data

  integer :: start_time, end_time, clock_rate
  real :: elapsed_time

  call system_clock(start_time,clock_rate)
  
  call create_req_dirs
  call define_input_params(input_data)
  call compute_contour_params(input_data,contour_data)
  call compute_contours(input_data,contour_data)
  
  call precompute(input_data,contour_data)
!
  if ((input_data%farswitch == 1) .OR. (input_data%farswitch == 2)) then

    print*,''
    print*,'solve: Near-field computation only:'
    print*,'' 

  else

    call meshgrid(input_data)

    print*,'solve: Now computing F+. This takes a while:'
    
    call compute_fplus(input_data,contour_data)

    print*,''
    print*,'solve: Now starting the IFT:'
    print*,''

    call computeift(input_data,contour_data)

    print*,'solve: IFT computed.'
    print*,''

 end if

 call system_clock(end_time)
 elapsed_time = real(end_time - start_time)/real(clock_rate)

 print*, 'Elapsed time for entire computation (seconds):', elapsed_time

END PROGRAM main   
