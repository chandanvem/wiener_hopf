PROGRAM main

  USE bessel_utils  
  USE io_utils
  USE contour_utils
  USE omp_lib
  USE input_params

  IMPLICIT none

  type(input_params_t)   :: input_data
  type(contour_params_t) :: contour_data
  
  call create_req_dirs
  call define_input_params(input_data)
  call compute_contour_params(input_data,contour_data)
 

END PROGRAM main   
