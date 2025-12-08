Module user_defined_farfield
  
  USE input_params
  USE bessel_utils
  USE kernel_integral_utils 
  USE user_defined_functions 
  USE user_defined_fplus
 
  IMPLICIT NONE
  
  PUBLIC  :: compute_directivity 

  CONTAINS 
 
  SUBROUTINE compute_directivity(input_data,contour_data,phi_list_index) 

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the asymptotic (steepest descent) directivity
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    type(input_params_t)                   :: input_data
    type(contour_params_t)                 :: contour_data
    integer                                :: phi_list_index
    real(dpk)                              :: phi_value,Theta_updated, Theta
    real(dpk)                              :: PI
    complex(dpk)                           :: saddle_point, Dmn_value, Dmn_num, Dmn_den
    complex(dpk)                           :: hank_arg, hank_arg_num, hank_arg_den

    PI = 4._dpk*ATAN(1.)
 
    phi_value = input_data%phi_list(phi_list_index)
    Theta = ATAN(SQRT(1._dpk - (input_data%kapT**2)*input_data%M2*input_data%M2)*TAN(phi_value))

    !if (phi_value == 0.) then
    !   Theta_updated = 0.
    !else if (phi_value == PI) then
    !   Theta_updated = PI
    !else if (TAN(phi_value) < 0.) then
    !   Theta_updated = Theta + PI
    !else 
    !   Theta_updated = Theta
    !END if

    saddle_point =  input_data%kapT*COS(phi_value)/SQRT(1._dpk - (input_data%kapT*input_data%M2*SIN(phi_value))**2)
    saddle_point =  saddle_point - input_data%kapT**2*input_data%M2
    saddle_point  = saddle_point/(1._dpk - (input_data%kapT**2)* (input_data%M2*input_data%M2))  ! stationary point

    Dmn_num = (input_data%omega_r)*(1._dpk - ((saddle_point*input_data%M2)**2))
    Dmn_num = Dmn_num*get_fplus_value(saddle_point,input_data,contour_data) 

    hank_arg_num = input_data%kapT*input_data%omega_r*SIN(phi_value)
    hank_arg_den = SQRT(1._dpk - ((input_data%kapT*input_data%M2*SIN(phi_value))**2))
    hank_arg = hank_arg_num/hank_arg_den
      
    Dmn_den = PI*dhank1(hank_arg,input_data%azim_mode,1)*input_data%kapT*SIN(phi_value)
    
    Dmn_value = Dmn_num/Dmn_den
    !Dmn_value = 20._dpk*LOG10(ABS(Dmn_value))
    Dmn_value = CMPLX(20._dpk*LOG10(ABS(Dmn_value)),0._dpk,kind=dpk)
    input_data%Dmn_list(phi_list_index) = Dmn_value

  END SUBROUTINE compute_directivity


END MODULE user_defined_farfield


!    n = wr/PI*(1._dpk - saddle_point*input_data%M2)**2*fplus(CMPLX(saddle_point,0._dpk,kind=dpk))

!    d1 = input_data%kapT*wr*SIN(phi_value)/SQRT(1._dpk - input_data%kapT**2*input_data%M2**2*SIN(phi_value)*SIN(phi_value))
    
!    d = dhank1(d1,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*d1)*input_data%kapT*SIN(phi_value)

!    dmn = n/d


