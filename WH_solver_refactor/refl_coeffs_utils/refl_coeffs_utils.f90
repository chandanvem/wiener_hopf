Module refl_coeffs_utils
  
  USE input_params 
  USE user_defined_fplus
  USE omp_lib

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC  :: init_refl_coeffs, compute_refl_coeffs

  CONTAINS 

  SUBROUTINE init_refl_coeffs(input_data,contour_data)

    type(input_params_t)          :: input_data
    type(contour_params_t)        :: contour_data

  END SUBROUTINE init_refl_coeffs

  SUBROUTINE compute_refl_coeffs(input_data,contour_data)

    type(input_params_t)          :: input_data
    type(contour_params_t)        :: contour_data


  END SUBROUTINE compute_refl_coeffs


 END MODULE refl_coeffs_utils

