Module user_defined_fplus

  USE input_params
  USE bessel_utils
  USE kernel_integral_utils 

  IMPLICIT NONE

  PUBLIC  :: get_fplus_value 

  CONTAINS 

  FUNCTION get_fplus_value(s_target,input_data,contour_data)  result(fplus)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Actual computation of \xi^{+}(s)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)       :: fplus, s_target, gpz, kpz, int_A1_at_s_target
    complex(dpk)       :: lambda_1_minus_mu_plus, lambda_2_minus_mu_plus, lambda_3_minus_mu_plus
    complex(dpk)       :: lambda_1_plus_mu_plus, lambda_2_plus_mu_plus, lambda_3_plus_mu_plus, lambda_prod
    integer            :: f
    real(dpk)          :: PI

    type(input_params_t) :: input_data
    type(contour_params_t) :: contour_data


    PI = 4._dpk*ATAN(1.)
    lambda_1_minus_mu_plus = sqrt(1._dpk - input_data%mu_plus*(input_data%M1 - 1._dpk))
    lambda_2_minus_mu_plus = sqrt(input_data%kapT - input_data%mu_plus*(input_data%kapT*input_data%M2 - 1._dpk))
    lambda_3_minus_mu_plus = sqrt(input_data%kapT - input_data%mu_plus*(input_data%kapT*input_data%M3 - 1._dpk))
    lambda_1_plus_mu_plus = sqrt(1._dpk - s_target*(input_data%M1 + 1._dpk))
    lambda_2_plus_mu_plus = sqrt(input_data%kapT - s_target*(input_data%kapT*input_data%M2 + 1._dpk))
    lambda_3_plus_mu_plus = sqrt(input_data%kapT - s_target*(input_data%kapT*input_data%M3 + 1._dpk))

    gpz = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M2)/ &
                  (input_data%k_minus_at_mu_plus*(input_data%mu_plus - input_data%KH_zero_2))

    call compute_eqn_A1_integral(s_target,int_A1_at_s_target,0,0,1,input_data,contour_data)

    if ((input_data%farswitch == 1) .OR. (input_data%farswitch == 2)) then
       if (REAL(s_target) >= contour_data%cont_cross_over_pt) then
          kpz = EXP(-int_A1_at_s_target/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) + & 
               LOG(compute_kernel(0,s_target,input_data)/compute_U_s_factor(0,s_target,input_data)))  !! s_target is above contour; cont_cross_over_pt is crossover pt
       else
          kpz = EXP(-int_A1_at_s_target/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! s_target is below
       end if
    else
       kpz = EXP(-int_A1_at_s_target/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! s_target is always below normally
    end if

    fplus = gpz*( (s_target - input_data%KH_zero_2)/(input_data%mu_plus - s_target) &
                           + input_data%vs_param_gamma)/(kpz*compute_U_s_factor(0,s_target,input_data))

  END FUNCTION get_fplus_value

END MODULE user_defined_fplus
