Module user_defined_precompute
  
  USE input_params
  USE bessel_utils
  USE kernel_integral_utils 
  USE user_defined_functions 
 
  IMPLICIT NONE
  
  PUBLIC  :: precompute 

  CONTAINS 
 
  SUBROUTINE precompute(input_data,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the various parameters needed while computing the IFT contour
!! 2. A check if the radial wavenumber is negative
!! 3. Constants computed: K^{-}(\mu_{mn}^{+}); K^{+}(s_{z1}); K^{+}(s_{z2})
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)          :: intgrl_A1_at_mu_plus, intgrl_A1_at_KH_zero_1, intgrl_A1_at_KH_zero_2, intgrl_at_sup_zero
    complex(dpk)          :: k_plus_at_mu_plus, f1, f2, f3
    integer               :: f, i1
    real(dpk)             :: PI

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data 
 

    PI = 4._dpk*ATAN(1.)
 
    if ((input_data%vortswitch == 1) .OR. (input_data%vortswitch == 2)) then

      input_data%psi = 1

    else

!! check if the radial wavenumber is negative:

       if((REAL(input_data%mu_plus) < -(input_data%kapT/(1._dpk - ((input_data%kapT)*(input_data%M1)) ))) .OR. &
                                                      (REAL(input_data%mu_plus) > (1._dpk/(1._dpk + input_data%M1)))) then 

          print*,'Radial wave number(s) is negative: Check the incident wave axial wave number!'
          STOP
      end if

      input_data%alpha1 = input_data%omega_r*SQRT((1._dpk - input_data%mu_plus*input_data%M1)**2 - (input_data%mu_plus)**2)

       print*,'precompute: The radial wave numbers:'
       write(*,'(/A12,2X,2F15.10)') ' alpha1:->', input_data%alpha1

!! the factor \Psi_{mn}(1) of (3.30) [see the JFM]:

       f1 = bessj(input_data%alpha1,input_data%azim_mode,1)
       f1 = f1*EXP(ABS(AIMAG(input_data%alpha1*input_data%h)))
      
       input_data%psi = f1
       print*, 'precompute: Evaluated psi for the incident wave'
   
    end if
    
!!  the factor Kt^{-}(\mu_{mn}^{+}):

    if (input_data%vortswitch .EQ. 0) then
     
        if (REAL(input_data%mu_plus) < contour_data%cont_cross_over_pt) then  !! mu_plus is below the contour; cont_cross_over_pt being the crossover pt
           print*, ''
           print*, 'precompute: The incident acoustic mode needs to be INSIDE the contour'
           print*, ''
          STOP
       end if

     end if
     
    print*, 'precompute: Evaluating kernel at mu_plus'

    call compute_eqn_A1_integral(input_data%mu_plus,intgrl_A1_at_mu_plus,0,0,1,input_data,contour_data)
       
    k_plus_at_mu_plus = EXP(-intgrl_A1_at_mu_plus/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) + & 
                        LOG(compute_kernel(0,input_data%mu_plus,input_data)/compute_U_s_factor(0,input_data%mu_plus,input_data)))

  !  write(*,'(/A12,2X,2F15.10)') ' integral at mu_plus:->', intgrl_A1_at_mu_plus
 
    input_data%k_minus_at_mu_plus =  compute_kernel(0,input_data%mu_plus,input_data)/ &
                          (k_plus_at_mu_plus*compute_U_s_factor(0,input_data%mu_plus,input_data)) 

!!  the factor Kt^{+}(s_{z1}):

    if (input_data%vortswitch .EQ. 0) then
  
       print*, 'precompute: Evaluating Kt^{+}(s_{z1}) at KH_zero_1'
       call compute_eqn_A1_integral(input_data%KH_zero_1,intgrl_A1_at_KH_zero_1,0,0,1,input_data,contour_data)

       input_data%k_plus_sz1 = EXP(-intgrl_A1_at_KH_zero_1/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk))) 
       !! NOTE: zero KH_zero_1 has to lie below  !! the contour

   !    write(*,'(/A22,2X,2F20.10/)') 'Integral at KH zero 1:->', intgrl_A1_at_KH_zero_1

    end if
!!  the factor K^{+}(s_{z2}):

    print*, 'precompute: Evaluating Kt^{+}(s_{z2}) at KH_zero_2'
    call compute_eqn_A1_integral(input_data%KH_zero_2,intgrl_A1_at_KH_zero_2,0,0,1,input_data,contour_data)

    input_data%k_plus_sz2 = EXP(-intgrl_A1_at_KH_zero_2/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! same for KH_zero_2

  !  write(*,'(/A22,2X,2F20.10/)') 'Integral at KH zero 2:->', intgrl_A1_at_KH_zero_2



  END SUBROUTINE precompute

END MODULE user_defined_precompute 
