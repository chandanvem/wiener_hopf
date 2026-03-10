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
    integer                       :: ios, freq_idx
    real(dpk)                     :: KH_real_part,&
                                     KH_imag_part
    real(dpk)                     :: d1,d2,d3,d4,d5
 
    input_data%RC_num_freq = 0
    

    open(10,file='RC_input_data.dat',status='old')

    do
        read(10,*,iostat=ios) d1, d2, d3, d4, d5
        if(ios /= 0) exit
        input_data%RC_num_freq = input_data%RC_num_freq + 1
    end do

    rewind(10)

    write(*,'(/A,2X,I10)') 'Number of frequencies:->', input_data%RC_num_freq

    allocate(input_data%RC_St_list(input_data%RC_num_freq))
    allocate(input_data%RC_KH_list(input_data%RC_num_freq))
    allocate(input_data%RC_k_plus_list(input_data%RC_num_freq))
    allocate(input_data%RC_GJM_list(input_data%RC_num_freq))

   do freq_idx = 1, input_data%RC_num_freq

       read(10,*,iostat=ios) d1, d2, d3, d4, d5

       input_data%RC_St_list(freq_idx) = d1
       KH_real_part = d2
       KH_imag_part = d3
       input_data%RC_GJM_list(freq_idx) = d4
       input_data%RC_k_plus_list = d5

       input_data%RC_KH_list(freq_idx) =  REAL(KH_real_part) &
                                        + CMPLX(0,1._dpk)*REAL(KH_imag_part)


       input_data%RC_GJM_list(freq_idx) =  REAL(input_data%RC_GJM_list(freq_idx)) &
                                           + 1E-20*CMPLX(0,1._dpk)

       input_data%RC_k_plus_list(freq_idx) =  REAL(input_data%RC_k_plus_list(freq_idx)) &
                                           + 1E-20*CMPLX(0,1._dpk)



       write(*,'(/A,2X,I10)') '=============== Frequency index =============', freq_idx

       write(*,'(4(2F14.6,3X))')  input_data%RC_St_list(freq_idx), &
                                  input_data%RC_KH_list(freq_idx), &
                                   input_data%RC_GJM_list(freq_idx), &
                                   input_data%RC_k_plus_list(freq_idx)

  end do
   

   close(10)

  END SUBROUTINE init_refl_coeffs

  SUBROUTINE compute_refl_coeffs(input_data,contour_data)

    complex(dpk)          :: intgrl_A1_at_mu_plus, intgrl_A1_at_KH_zero_1, intgrl_A1_at_KH_zero_2, intgrl_at_sup_zero
    complex(dpk)          :: intgrl_A1_at_k_d_plus, intgrl_A1_at_duct_mode
    complex(dpk)          :: intgrl_A1_at_GJ
    complex(dpk)          :: k_plus_at_mu_plus, k_plus_at_GJ
    complex(dpk)          :: f1, f2, f3, A_mn_k_plus_op
    complex(dpk)          :: lambda1_at_k_d_plus, A_mn_k_plus_num, A_mn_k_plus_den, limit_term_L_k_d_plus
    integer               :: f, i1, j, freq_idx
    real(dpk)             :: PI, res_mu_plus, res_KH_1, res_k_d_plus
    real(dpk)             :: B_mn

    complex(dpk)          :: GJ_mode_amp
    complex(dpk), allocatable, dimension(:) :: RC_coeffs

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data

    allocate(RC_coeffs(input_data%RC_num_freq))
    
    PI = 4._dpk*ATAN(1.)

    do freq_idx = 1, input_data%RC_num_freq

         write(*,'(/A,2X,I10)') '=============== Frequency index ', freq_idx, '=============='

         input_data%mu_plus = input_data%RC_GJM_list(freq_idx)
         input_data%KH_zero_1 = input_data%RC_KH_list(freq_idx)
         input_data%k_d_plus = input_data%RC_k_plus_list(freq_idx)
         input_data%omega_st = input_data%RC_St_list(freq_idx) 

         input_data%omega_r = PI*input_data%M1*input_data%omega_st
 
         
         input_data%alpha1 = input_data%omega_r*SQRT((1._dpk - input_data%mu_plus*input_data%M1)**2&
                                                             - (input_data%mu_plus)**2)

         f1 = bessj(input_data%alpha1,input_data%azim_mode,1)

         B_mn = 1._dpk

         input_data%psi = f1

         input_data%C0 = input_data%psi*B_mn


         res_mu_plus =  ABS(compute_kernel(1,input_data%mu_plus,input_data))
         res_KH_1    =  ABS(compute_kernel(1,input_data%KH_zero_1,input_data))
         res_k_d_plus = ABS(compute_kernel(1,input_data%k_d_plus,input_data))

         write(*,'(/A,2X,2F15.10)') ' Value of K(s) at s = mu plus:->', res_mu_plus
         write(*,'(/A,2X,2F15.10)') ' Value of K(s) at s = KH zero:->', res_KH_1
         write(*,'(/A,2X,2F15.10)') ' Value of K(s) at s = k_d_plus:->', res_k_d_plus

         print*, 'compute_refl_coeffs: Evaluating kernel integral at mu_plus'

         call compute_eqn_A1_integral(input_data%mu_plus,intgrl_A1_at_mu_plus,0,0,1,'not_derivative','guided_jet',&
                                                                          input_data,contour_data)
         k_plus_at_mu_plus = EXP(-intgrl_A1_at_mu_plus/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) + &
                             LOG(dkernel_ds(input_data%mu_plus,input_data)/compute_U_s_factor(input_data%mu_plus,&
                                                                              input_data,'not_k_d_plus')))


         input_data%k_plus_at_mu_plus = k_plus_at_mu_plus

         write(*,'(/A,2X,2F15.10)') 'compute_refl_coeffs: k+ at  mu_plus:->',input_data%k_plus_at_mu_plus

         input_data%k_minus_at_mu_plus = dkernel_ds(input_data%mu_plus,input_data)/(input_data%mu_plus - input_data%KH_zero_1)

         write(*,'(/A,2X,2F20.10/)') 'compute_refl_coeffs: K- at mu+ :->', input_data%k_minus_at_mu_plus

         print*, 'compute_refl_coeffs: Evaluating Kt^{+}(s_{z1}) at KH_zero_1'
         call compute_eqn_A1_integral(input_data%KH_zero_1,intgrl_A1_at_KH_zero_1,0,0,1,'not_derivative','KH_mode_1',&
                                                                                            input_data,contour_data)
         input_data%k_plus_sz1 = EXP(-intgrl_A1_at_KH_zero_1/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))

         write(*,'(/A,2X,2F20.10/)') 'compute_refl_coeffs: K+ at KH1   :->', input_data%k_plus_sz1

         call compute_eqn_A1_integral(input_data%k_d_plus,intgrl_A1_at_k_d_plus,0,0,1,'not_derivative','k_d_plus',&
                                                                            input_data,contour_data)

          limit_term_L_k_d_plus = &
          EXP(-intgrl_A1_at_k_d_plus/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))
     !+&
     !              LOG(dkernel_ds(input_data%k_d_plus,input_data)/compute_U_s_factor(input_data%k_d_plus,&
      !                                                          input_data,'not_k_d_plus')))

          write(*,'(/A,2X,2F15.10)') 'compute_refl_coeffs:  limit_term_L_k_d_plus:->',&
                                                                            limit_term_L_k_d_plus


          lambda1_at_k_d_plus = sqrt(1._dpk - input_data%k_d_plus*(input_data%M1+1._dpk)) *&
                                 sqrt(1._dpk-  input_data%k_d_plus*(input_data%M1-1._dpk))

          A_mn_k_plus_num = CMPLX(0._dpk,1._dpk,kind=dpk)*((input_data%omega_r)**2)*&
                                    input_data%C0*(1 - (input_data%mu_plus*input_data%M1))*&
                                    (1 - (input_data%M1*input_data%k_d_plus))**2

          A_mn_k_plus_den  = limit_term_L_k_d_plus*input_data%k_minus_at_mu_plus*&
                              (input_data%k_d_plus - input_data%mu_plus)
          A_mn_k_plus_den  = A_mn_k_plus_den*compute_U_s_factor(input_data%k_d_plus,input_data,'non_k_d_plus')
          A_mn_k_plus_den =  A_mn_k_plus_den*lambda1_at_k_d_plus*&
                            dbessj(lambda1_at_k_d_plus*input_data%omega_r,input_data%azim_mode,1)

          A_mn_k_plus_op = A_mn_k_plus_num/A_mn_k_plus_den

          GJ_mode_amp = input_data%omega_r*(1._dpk - (input_data%mu_plus*input_data%M1))

              write(*,'(/A,2X,2F15.10)') 'compute_refl_coeffs:  guided jet mode amplitude:->',&
                                                                             GJ_mode_amp

          RC_coeffs(freq_idx) = ABS(A_mn_k_plus_op)/ABS(GJ_mode_amp)

          write(*,'(/A,2X,2F15.10)') 'compute_refl_coeffs:  k plus reflection coeff. :=', RC_coeffs(freq_idx)
                                                                             
     end do

     open(10,file='RC_coeffs.out',form='FORMATTED')
     do freq_idx = 1, input_data%RC_num_freq
        write(10,'(F8.2,2F20.15)')  input_data%RC_St_list(freq_idx),RC_coeffs(freq_idx)
     end do
     close(10)


  END SUBROUTINE compute_refl_coeffs


 END MODULE refl_coeffs_utils

