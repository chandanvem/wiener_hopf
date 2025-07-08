PROGRAM main

  USE bessel_utils  
  USE io_utils
  USE contour_generate_utils
  USE contour_init_utils
  USE omp_lib
  USE input_params

  IMPLICIT none

  type(input_params_t)   :: input_data
  type(contour_params_t) :: contour_data
  
  call create_req_dirs
  call define_input_params(input_data)
  call compute_contour_params(input_data,contour_data)
  call compute_contours(input_data,contour_data)
  
  call precompute(input_data,contour_data)

  if ((farswitch == 1) .OR. (farswitch == 2)) then

    print*,''
    print*,'solve: Near-field computation only:'
    print*,'' 

  else

    call meshgrid

    print*,'solve: Now computing F+. This takes a while:'
    call compute_fplus(restart,0)

    print*,''
    print*,'solve: Now starting the IFT:'
    print*,''

    call computeift

  end if


END PROGRAM main   
