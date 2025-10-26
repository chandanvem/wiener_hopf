Module user_defined_functions

  USE input_params
  USE bessel_utils
  USE io_utils  

  IMPLICIT NONE

  PUBLIC  :: compute_U_s_factor, compute_kernel, &
             integrand_ifT_pot, integrand_ifT_pr,&
             compute_Trs_lambda, compute_d_ds_Trs_lambda

  CONTAINS 

  FUNCTION compute_U_s_factor(s_target,input_data) result(u_s)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The factor U(s)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)  :: s_target, u_s
    type(input_params_t) :: input_data

    if (input_data%solution_mode == 'guided_jet') then
       u_s = (s_target-input_data%KH_zero_1)
        !*(s_target - input_data%s_GJ)
    else 
       u_s = (s_target-input_data%KH_zero_1)
    end if

  END FUNCTION compute_U_s_factor


  FUNCTION compute_kernel(ss,zeta,input_data)  result(kernel)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the kernel function (3.22) or (4.10)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

   complex(dpk)    :: kernel, zeta
   complex(dpk)    :: lambda1, lambda2
   complex(dpk)    :: F1n, F1d, F1, F1f, F2n, F2d, F2, F2f
   integer         :: ss
   real(dpk)       :: s1_1, s2_1, s1_2, s2_2
   type(input_params_t) :: input_data

   lambda1 =   sqrt(1._dpk - zeta*(input_data%M1+1._dpk))*sqrt(1._dpk - zeta*(input_data%M1-1._dpk))
   lambda2 =   sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M2+1._dpk))* &
               sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M2-1._dpk))


   if ((ABS(lambda1*input_data%omega_r) < input_data%asymplim)) then

      F1n =   bessj(lambda1*input_data%omega_r,input_data%azim_mode,1)
      F1d =  dbessj(lambda1*input_data%omega_r,input_data%azim_mode,1)
      F1f =  F1n/F1d

   else   ! asymtotic limit  ~ i

      F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
      
   end if

   F1 = ((input_data%kap_rho*((1._dpk - zeta*input_data%M1)**2))/lambda1)*F1f

   if ( ABS(AIMAG(lambda2*input_data%omega_r)) < input_data%asymplim1) then
      
      F2n =  hank1(lambda2*input_data%omega_r,input_data%azim_mode,1)
      F2d = dhank1(lambda2*input_data%omega_r,input_data%azim_mode,1)
      F2f = F2n/F2d
      
   else   ! asymtotic limit   

      F2f = (8._dpk*lambda2*input_data%omega_r +&
            4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%azim_mode*input_data%azim_mode - &
                   CMPLX(0._dpk,1._dpk,kind=dpk))

      F2f = F2f/ &
           (8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*input_data%omega_r - &
            4._dpk*input_data%azim_mode*input_data%azim_mode - 3._dpk)

   end if

   F2 =(((1._dpk - zeta*input_data%M2)**2)/lambda2)*F2f
      

   kernel = input_data%omega_r*(F1 - F2)

  END FUNCTION compute_kernel

  FUNCTION integrand_ifT_pr(ri,zi,IFT_contour_idx,input_data,contour_data) result(integrandIFTpr)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the ifT integrand when computing for prestream_idxure
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)       :: ri, zi
    integer         :: IFT_contour_idx, stream_idx
    complex(dpk)    :: u
    complex(dpk)    :: integrandIFTpr

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data

    u = contour_data%IFTpoints(IFT_contour_idx)

    if (ri .LE. 1) then 
       stream_idx = 1
    else
       stream_idx = 2
    end if

    if (stream_idx==1) then

       integrandIFTpr = (1._dpk - u*input_data%M1)**2
       integrandIFTpr = integrandIFTpr*compute_Trs_lambda(ri,u,stream_idx,input_data)* &
                        input_data%fplusz(IFT_contour_idx)* & 
                        EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)

    else if (stream_idx==2) then
       
       integrandIFTpr = (1._dpk - u*input_data%M2)**2
       integrandIFTpr = integrandIFTpr*compute_Trs_lambda(ri,u,stream_idx,input_data)* &
                        input_data%fplusz(IFT_contour_idx)* &
                        EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)
   end if


  END FUNCTION integrand_ifT_pr

  FUNCTION integrand_ifT_pot(ri,zi,IFT_contour_idx,input_data,contour_data) result(integrandIFTpot)
 !
 !!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
 !!1. Compute the ifT integrand when computing for potential
 !!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
 !
    real(dpk)       :: ri, zi
    integer         :: IFT_contour_idx, stream_idx
    complex(dpk)    :: u
    complex(dpk)    :: integrandIFTpot
 !
    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data
 !
 !
    if (ri .LE. 1) then 
       stream_idx = 1
    else
       stream_idx = 2
    end if
 !
    u = contour_data%IFTpoints(IFT_contour_idx)
 !
    if (stream_idx==1) then
 !
 !!potential:
 !
       integrandIFTpot = ((1._dpk - u*input_data%M1))
       integrandIFTpot = integrandIFTpot *compute_Trs_lambda(ri,u,stream_idx,input_data)* &
                         input_data%fplusz(IFT_contour_idx)* & 
                          EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)
 !
    else if (stream_idx==2) then
 !
 !!potential:
       
       integrandIFTpot = ((1._dpk - u*input_data%M2))
       integrandIFTpot = integrandIFTpot*compute_Trs_lambda(ri,u,stream_idx,input_data)* &
                         input_data%fplusz(IFT_contour_idx)* &
                         EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)
 !
    end if
 !
 !
  END FUNCTION integrand_ifT_pot
 !
 
  FUNCTION compute_Trs_lambda(ri,si,stream_idx,input_data) result(Trs)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute Trs of (3.33)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: ri
    complex(dpk)  :: sri, Trs
    complex(dpk)  :: si
    complex(dpk)  :: lambda1, lambda2
    complex(dpk)  :: F1n, F1d, F2n, F2d, F1f, F2f
    integer       :: stream_idx
    real(dpk)     :: s1_1, s2_1, s1_2, s2_2
    logical       :: NAN_flag 

    type(input_params_t)  :: input_data

    NAN_flag = .FALSE.

    lambda1 = sqrt(1._dpk - si*(input_data%M1+1._dpk))*sqrt(1._dpk - si*(input_data%M1-1._dpk))

    lambda2 = sqrt(input_data%kapT - si*(input_data%kapT*input_data%M2+1._dpk))* &
              sqrt(input_data%kapT - si*(input_data%kapT*input_data%M2-1._dpk))

    if (stream_idx==1) then
       
       if (ABS(lambda1*input_data%omega_r) < input_data%asymplim) then

          F1n =  bessj(lambda1*input_data%omega_r*ri,input_data%azim_mode,1)
          F1n = F1n*EXP(ABS(AIMAG(lambda1*input_data%omega_r*ri))) 
  
          F1d = dbessj(lambda1*input_data%omega_r,input_data%azim_mode,1)
          F1d = F1d*EXP(ABS(AIMAG(lambda1*input_data%omega_r))) 

          F1f = F1n/F1d

       else   ! asymtotic limit  ~ i

          F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
          
       end if

       Trs = F1f/lambda1

   else if (stream_idx==2) then

       if ( ABS(AIMAG(lambda2*input_data%omega_r)) < input_data%asymplim1) then
      
         F2n =  hank1(lambda2*input_data%omega_r*ri,input_data%azim_mode,1)
         F2n = F2n*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*input_data%omega_r*ri)

         F2d = dhank1(lambda2*input_data%omega_r,input_data%azim_mode,1)
         F2d = F2d*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*input_data%omega_r)

         F2f = F2n/F2d
       
       else   ! asymtotic limit   

         F2f = (8._dpk*lambda2*input_data%omega_r +&
               4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%azim_mode*input_data%azim_mode - &
              CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*input_data%omega_r - &
              4._dpk*input_data%azim_mode*input_data%azim_mode - 3._dpk)
        
       end if
      
       Trs = F2f/lambda2

   end if

   call check_NAN(Trs,NAN_flag)

   if(NAN_flag) then
      print*,'compute_Trs_lambda: NAN found at s = ', si
      STOP
   end if

  END FUNCTION compute_Trs_lambda


FUNCTION compute_d_ds_Trs_lambda(ri, si, stream_idx, input_data) result(dTrs)
  
  real(dpk)            :: ri
  complex(dpk)         :: si, dTrs
  complex(dpk)         :: lambda1, lambda2, dlambda1_ds, dlambda2_ds
  complex(dpk)         :: term_1_num, term_1_den, term_1
  complex(dpk)         :: term_2_num, term_2_den, term_2
  complex(dpk)         :: pre_factor
  real(dpk)            :: omega
  integer              :: stream_idx
  logical              :: NAN_flag

  type(input_params_t)  :: input_data

  NAN_flag = .FALSE.

  omega = input_data%omega_r

  lambda1 = sqrt(1._dpk - si*(input_data%M1+1._dpk)) * sqrt(1._dpk -si*(input_data%M1-1._dpk))
  lambda2 = sqrt(input_data%kapT - si*(input_data%kapT*input_data%M2+1._dpk)) *&
            sqrt(input_data%kapT - si*(input_data%kapT*input_data%M2-1._dpk))

  dlambda1_ds = -(input_data%M1 - si*( (input_data%M1)**2 - 1))
  dlambda1_ds = dlambda1_ds/lambda1

  dlambda2_ds = -( ((input_data%kapT)**2)*input_data%M2 - si*((input_data%kapT*input_data%M1)**2 - 1))
  dlambda2_ds = dlambda2_ds/lambda2

  if (stream_idx == 1) then
     
     pre_factor = dlambda1_ds*omega
     
     term_1_num = ri*dbessj(lambda1*omega*ri, input_data%azim_mode, 1)
     term_1_den = dbessj(lambda1*omega, input_data%azim_mode, 1)
     term_1     = term_1_num/term_1_den

     term_2_num =   bessj(lambda1*omega*ri, input_data%azim_mode, 1)
     term_2_num =   term_2_num*d2bessj(lambda1*omega*ri, input_data%azim_mode, 1)
     term_2_den =   (dbessj(lambda1*omega, input_data%azim_mode, 1))**2
     term_2     =   term_2_num/term_2_den

  else if (stream_idx == 2 ) then

     pre_factor = dlambda2_ds*omega
     
     term_1_num = ri*dhank1(lambda2*omega*ri, input_data%azim_mode, 1)
     term_1_den = dhank1(lambda2*omega, input_data%azim_mode, 1)
     term_1     = term_1_num/term_1_den

     term_2_num =   hank1(lambda2*omega*ri, input_data%azim_mode, 1)
     term_2_num =   term_2_num*d2hank1(lambda2*omega*ri, input_data%azim_mode, 1)
     term_2_den =   (dhank1(lambda2*omega, input_data%azim_mode, 1))**2
     term_2     =   term_2_num/term_2_den

  end if 

  dTrs = pre_factor*(term_1 - term_2)

  call check_NAN(dTrs,NAN_flag)

  if(NAN_flag) then
      print*,'compute_dTrs_ds: NAN found at s = ', si
      STOP
  end if


END FUNCTION compute_d_ds_Trs_lambda


END MODULE user_defined_functions


