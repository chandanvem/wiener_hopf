PROGRAM main

  USE bessel_utils  
  USE io_utils
  USE omp_lib
  USE input_params

  IMPLICIT none

  type(input_params_t) :: input_params_data

  call define_input_params(input_params_data)


END PROGRAM main   
