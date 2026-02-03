Module user_defined_precompute
  
  USE input_params
  USE bessel_utils
  USE kernel_integral_utils 
  USE user_defined_functions 
  USE user_defined_fplus

 
  IMPLICIT NONE
  
  PUBLIC  :: precompute 
  PRIVATE :: precompute_hard_duct_mode, precompute_guided_jet_mode


  CONTAINS 
 
  SUBROUTINE precompute(input_data,contour_data)

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data

    if (input_data%solution_mode == 'guided_jet') then
       call precompute_guided_jet_mode(input_data,contour_data)
    else
       call precompute_hard_duct_mode(input_data,contour_data)
    end if


  END SUBROUTINE precompute
 
  SUBROUTINE precompute_hard_duct_mode(input_data,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the various parameters needed while computing the IFT contour
!! 2. A check if the radial wavenumber is negative
!! 3. Constants computed: K^{-}(\mu_{mn}^{+}); K^{+}(s_{z1}); K^{+}(s_{z2})
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)          :: intgrl_A1_at_mu_plus, intgrl_A1_at_KH_zero_1, intgrl_A1_at_KH_zero_2, intgrl_at_sup_zero
    complex(dpk)          :: k_plus_at_mu_plus, f1, f2, f3
    integer               :: f, i1
    real(dpk)             :: PI, res_mu_plus, res_KH_1

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data 
 
    PI = 4._dpk*ATAN(1.)
 
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
    f1 = f1*EXP(ABS(AIMAG(input_data%alpha1)))
    
    input_data%psi = f1
    write(*,'(/A12,2X,2F15.10)') ' psi:->', input_data%psi
    print*, 'precompute: Evaluated psi for the incident wave'

    res_mu_plus = ABS(compute_kernel(1,input_data%mu_plus,input_data))
    res_KH_1    = ABS(compute_kernel(1,input_data%KH_zero_1,input_data))

    write(*,'(/A12,2X,2F15.10)') ' Residue mu plus:->', res_mu_plus
    write(*,'(/A12,2X,2F15.10)') ' Residue KH zero:->', res_KH_1

!!  the factor Kt^{-}(\mu_{mn}^{+}):

     
     if (REAL(input_data%mu_plus) < contour_data%cont_cross_over_pt) then  !! mu_plus is below the contour; cont_cross_over_pt being the crossover pt
        print*, ''
        print*, 'precompute: The incident acoustic mode needs to be INSIDE the contour'
        print*, ''
    !   STOP
    end if

     
    print*, 'precompute: Evaluating kernel at mu_plus'

    call compute_eqn_A1_integral(input_data%mu_plus,intgrl_A1_at_mu_plus,0,0,1,'not_derivative',&
                                                                   'mu_plus',input_data,contour_data)
       
    k_plus_at_mu_plus = EXP(-intgrl_A1_at_mu_plus/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) + & 
                        LOG(compute_kernel(0,input_data%mu_plus,input_data)/&
                                 compute_U_s_factor(input_data%mu_plus,input_data,'mu_plus')))

    write(*,'(/A12,2X,2F15.10)') ' integral at mu_plus:->', intgrl_A1_at_mu_plus
 
    input_data%k_minus_at_mu_plus =  compute_kernel(0,input_data%mu_plus,input_data)/ &
                          (k_plus_at_mu_plus*compute_U_s_factor(input_data%mu_plus,input_data,'mu_plus')) 

    write(*,'(/A22,2X,2F20.10/)') 'K- at mu+ :->', input_data%k_minus_at_mu_plus 

!!  the factor Kt^{+}(s_{z1}):

    print*, 'precompute: Evaluating Kt^{+}(s_{z1}) at KH_zero_1'
    call compute_eqn_A1_integral(input_data%KH_zero_1,intgrl_A1_at_KH_zero_1,0,0,1,&
                                                             'not_derivative','KH_zero',input_data,contour_data)

    input_data%k_plus_sz1 = EXP(-intgrl_A1_at_KH_zero_1/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk))) 
    !! NOTE: zero KH_zero_1 has to lie below  !! the contour

    write(*,'(/A22,2X,2F20.10/)') 'K+ at KH1   :->', input_data%k_plus_sz1

  
  END SUBROUTINE precompute_hard_duct_mode


  SUBROUTINE precompute_guided_jet_mode(input_data,contour_data)

    complex(dpk)          :: intgrl_A1_at_mu_plus, intgrl_A1_at_KH_zero_1, intgrl_A1_at_KH_zero_2, intgrl_at_sup_zero
    complex(dpk)          :: intgrl_A1_at_k_d_plus, intgrl_A1_at_duct_mode
    complex(dpk)          :: intgrl_A1_at_GJ
    complex(dpk)          :: k_plus_at_mu_plus, k_plus_at_GJ
    complex(dpk)          :: f1, f2, f3
    complex(dpk)          :: lambda1_at_k_d_plus, A_mn_k_plus_num, A_mn_k_plus_den, limit_term_L_k_d_plus
    integer               :: f, i1, j
    real(dpk)             :: PI, res_mu_plus, res_KH_1, res_k_d_plus
    real(dpk)             :: B_mn

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data 
 
    PI = 4._dpk*ATAN(1.)

 
    input_data%alpha1 = input_data%omega_r*SQRT((1._dpk - input_data%mu_plus*input_data%M1)**2&
                                                        - (input_data%mu_plus)**2)

    input_data%alpha2 = input_data%omega_r*SQRT(((input_data%kapT)**2)*(1._dpk -input_data%mu_plus*input_data%M2)**2 - &
                                                                                       (input_data%mu_plus)**2)


    f1 = bessj(input_data%alpha1,input_data%azim_mode,1)

    B_mn = 1._dpk

    input_data%psi = f1

    input_data%C0 = input_data%psi*B_mn

    write(*,'(/A,2X,2F15.10)') 'psi for guided jet mode :->', input_data%psi
 
    res_mu_plus =  ABS(compute_kernel(1,input_data%mu_plus,input_data))
    res_KH_1    =  ABS(compute_kernel(1,input_data%KH_zero_1,input_data))
    res_k_d_plus = ABS(compute_kernel(1,input_data%k_d_plus,input_data))

    write(*,'(/A,2X,2F15.10)') ' Value of K(s) at s = mu plus:->', res_mu_plus
    write(*,'(/A,2X,2F15.10)') ' Value of K(s) at s = KH zero:->', res_KH_1
    write(*,'(/A,2X,2F15.10)') ' Value of K(s) at s = k_d_plus:->', res_k_d_plus

    print*, 'precompute_guided_jet_mode: Evaluating kernel integral at mu_plus'

    call compute_eqn_A1_integral(input_data%mu_plus,intgrl_A1_at_mu_plus,0,0,1,'not_derivative','guided_jet',&
                                                                     input_data,contour_data)
    k_plus_at_mu_plus = EXP(-intgrl_A1_at_mu_plus/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk))) 
    input_data%k_plus_at_mu_plus = k_plus_at_mu_plus

    write(*,'(/A,2X,2F15.10)') 'precompute_guided_jet_mode: k+ at  mu_plus:->',input_data%k_plus_at_mu_plus
 
    input_data%k_minus_at_mu_plus = dkernel_ds(input_data%mu_plus,input_data)/(input_data%mu_plus - input_data%KH_zero_1)

    write(*,'(/A,2X,2F20.10/)') 'precompute_guided_jet_mode: K- at mu+ :->', input_data%k_minus_at_mu_plus 

    print*, 'precompute_guided_jet_mode: Evaluating Kt^{+}(s_{z1}) at KH_zero_1'
    call compute_eqn_A1_integral(input_data%KH_zero_1,intgrl_A1_at_KH_zero_1,0,0,1,'not_derivative','KH_mode_1',&
                                                                                       input_data,contour_data)
    input_data%k_plus_sz1 = EXP(-intgrl_A1_at_KH_zero_1/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk))) 
    
    write(*,'(/A,2X,2F20.10/)') 'precompute_guided_jet_mode: K+ at KH1   :->', input_data%k_plus_sz1

    if ((input_data%near_far_field_mode == 'near_field')) then

         call compute_eqn_A1_integral(input_data%k_d_plus,intgrl_A1_at_k_d_plus,0,0,1,'not_derivative','k_d_plus',&
                                                                            input_data,contour_data)

          limit_term_L_k_d_plus = &
          EXP(-intgrl_A1_at_k_d_plus/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) +&
              LOG(dkernel_ds(input_data%k_d_plus,input_data)/compute_U_s_factor(input_data%k_d_plus,&
                                                           input_data,'not_k_d_plus')))

         write(*,'(/A,2X,2F15.10)') 'precompute_guided_jet_mode:  k_plus_at_k_d_plus:->',&
                                                                           input_data%k_plus_at_k_d_plus                                                                


         lambda1_at_k_d_plus = sqrt(1._dpk - input_data%k_d_plus*(input_data%M1+1._dpk)) *&
                                sqrt(1._dpk-  input_data%k_d_plus*(input_data%M1-1._dpk))
 
         A_mn_k_plus_num = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%C0*(1 - (input_data%mu_plus*input_data%M1))*&
                                   (1 - (input_data%M1*input_data%k_d_plus))**2
  
         A_mn_k_plus_den  = limit_term_L_k_d_plus*input_data%k_minus_at_mu_plus*&
                             (input_data%k_d_plus - input_data%mu_plus)
         A_mn_k_plus_den  = A_mn_k_plus_den*compute_U_s_factor(input_data%k_d_plus,input_data,'non_k_d_plus')
         A_mn_k_plus_den =  A_mn_k_plus_den*lambda1_at_k_d_plus*&
                           dbessj(lambda1_at_k_d_plus*input_data%omega_r,input_data%azim_mode,1)

         input_data%A_mn_k_plus = A_mn_k_plus_num/A_mn_k_plus_den
     
         do j = 1, input_data%num_duct_modes
               
             if ( (j .EQ. 1) .AND. (input_data%azim_mode .EQ. 0)  )then
                  input_data%A_mn_duct_modes_list(j) =  & 
                     -2._dpk*input_data%omega_r*CMPLX(0._dpk,1._dpk,kind=dpk)*&
                     ((1 - (input_data%duct_modes_list(j)*input_data%M1))**2)*&
                      get_fplus_value(input_data%duct_modes_list(j),input_data,contour_data)
             else
                  input_data% A_mn_duct_modes_list(j) =  & 
                     input_data%omega_r*CMPLX(0._dpk,1._dpk,kind=dpk)*&
                     ((1 - (input_data%duct_modes_list(j)*input_data%M1))**2)*&
                      get_fplus_value(input_data%duct_modes_list(j),input_data,contour_data)
                   input_data%A_mn_duct_modes_list(j) = input_data%A_mn_duct_modes_list(j)/&
                                 (d2hank1(input_data%duct_modes_list(j), input_data%azim_mode, 1)*&
                                 (input_data%M1 + (input_data%duct_modes_list(j)*(1._dpk - ( (input_data%M1)**2)))))
             end if
            
             
             write(*,'(/A,2X,2F15.10)') 'precompute_guided_jet_mode:  A_mn_duct_modes_list:->',&
                                                                          input_data%A_mn_duct_modes_list(j)

         end do

     end if


  END SUBROUTINE precompute_guided_jet_mode


END MODULE user_defined_precompute 
