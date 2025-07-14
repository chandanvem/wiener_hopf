Module user_defined_functions

  USE input_params
  USE bessel_utils
  

  IMPLICIT NONE

  PRIVATE :: bc, compute_Trs_lambda
  PUBLIC  :: compute_U_s_factor, compute_kernel, &
             integrand_IFT_pot, integrand_IFT_pr, &
             residuepr, residuepot

  CONTAINS 

  FUNCTION compute_U_s_factor(ss,s_target,input_data) result(u_s)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The factor U(s)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)  :: s_target, u_s
    integer       :: j, jj, ss
    type(input_params_t) :: input_data

    u_s = (s_target-input_data%KH_zero_1)

  END FUNCTION compute_U_s_factor



  FUNCTION compute_kernel(ss,zeta,input_data)  result(kernel)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the kernel function (3.22) or (4.10)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

   complex(dpk)    :: kernel, zeta
   complex(dpk)    :: lambda1, lambda2
   complex(dpk)    :: F1n, F1d, F1, F1f, F2n, F2d, F2, F2f
   integer         :: ss
   type(input_params_t) :: input_data


   lambda1 =  sqrt(1._dpk - zeta*(input_data%M1+1._dpk))*sqrt(1._dpk - zeta*(input_data%M1-1._dpk))
   lambda2 =  sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M2+1._dpk))* &
               sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M2-1._dpk))


   if ((ABS(lambda1*input_data%omega_r) < input_data%asymplim .AND.&
                     ABS(AIMAG(lambda1*input_data%omega_r)) < input_data%asymplim1)) then

      F1n =   bessj(lambda1*input_data%omega_r,input_data%azim_mode,1)
      F1d =  dbessj(lambda1*input_data%omega_r,input_data%azim_mode,1)
      F1f =  F1n/F1d

   else   ! asymtotic limit  ~ i

      F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
      
   end if

   F1 = ((input_data%kap_rho*((1._dpk - zeta*input_data%M1)**2))/lambda1)*F1f

   if (ABS(lambda2*input_data%omega_r) < input_data%asymplim .AND. &
        ABS(AIMAG(lambda2*input_data%omega_r)) < input_data%asymplim1) then
      
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

!!$    print*,'F2n:',F2n
!!$    print*,'F2d:',F2d
!!$    print*,'F2f:',F2f

   F2 =(((1._dpk - zeta*input_data%M2)**2)/lambda2)*F2f
      

   kernel = input_data%omega_r*(F1 - F2)

  END FUNCTION compute_kernel



  FUNCTION bc(z,input_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the term to be factored out from the kernel function (3.22)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)  :: z, bc
    integer       :: i
    type(input_params_t) :: input_data

    bc = 1._dpk


  END FUNCTION bc



  FUNCTION integrand_IFT_pot(ri,zi,ift_contour_idx,stream_idx,input_data,contour_data) result(integrandiftpot)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the IFT integrand when computing for potential
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)       :: ri, zi
    integer         :: ift_contour_idx, stream_idx
    complex(dpk)    :: u
    complex(dpk)    :: integrandiftpot

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data

    u = contour_data%iftpoints(ift_contour_idx)

    if (stream_idx==1) then

!! potential:

       integrandiftpot = ((1._dpk - u*input_data%M1))
       integrandiftpot = integrandiftpot *compute_Trs_lambda(ri,u,stream_idx,input_data)* &
                         input_data%fplusz(ift_contour_idx)* & 
                          EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)

    else if (stream_idx==2) then

!! potential:
       
       integrandiftpot = ((1._dpk - u*input_data%M2))
       integrandiftpot = integrandiftpot*compute_Trs_lambda(ri,u,stream_idx,input_data)* &
                         input_data%fplusz(ift_contour_idx)* &
                         EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)

    end if


  END FUNCTION integrand_IFT_pot



  FUNCTION integrand_IFT_pr(ri,zi,ift_contour_idx,stream_idx,input_data,contour_data) result(integrandiftpr)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the IFT integrand when computing for prestream_idxure
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)       :: ri, zi
    integer         :: ift_contour_idx, stream_idx
    complex(dpk)    :: u
    complex(dpk)    :: integrandiftpr

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data

    u = contour_data%iftpoints(ift_contour_idx)

    if (stream_idx==1) then

       integrandiftpr = (1._dpk - u*input_data%M2)**2
       integrandiftpr = integrandiftpr*compute_Trs_lambda(ri,u,stream_idx,input_data)* &
                        input_data%fplusz(ift_contour_idx)* & 
                        EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)

    else if (stream_idx==2) then
       
       integrandiftpr = (1._dpk - u*input_data%M2)**2
       integrandiftpr = integrandiftpr*compute_Trs_lambda(ri,u,stream_idx,input_data)*input_data%fplusz(ift_contour_idx)* &
                        EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)
   end if


  END FUNCTION integrand_IFT_pr

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
  
    type(input_params_t)  :: input_data

    lambda1 = sqrt(1._dpk - si*(input_data%M1+1._dpk))*sqrt(1._dpk - si*(input_data%M1-1._dpk))

    lambda2 = sqrt(input_data%kapT - si*(input_data%kapT*input_data%M2+1._dpk))* &
              sqrt(input_data%kapT - si*(input_data%kapT*input_data%M2-1._dpk))

    if (stream_idx==1) then
       
       if ((ABS(lambda1*input_data%omega_r) < input_data%asymplim .AND.&
           ABS(AIMAG(lambda1*input_data%omega_r)) < input_data%asymplim1)) then

          F1n =  bessj(lambda1*input_data%omega_r,input_data%azim_mode,1)
          F1d = dbessj(lambda1*input_data%omega_r,input_data%azim_mode,1)
          F1f = F1n/F1d

       else   ! asymtotic limit  ~ i

          F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
          
       end if

       Trs = F1f/lambda1

   else if (stream_idx==2) then

       if (ABS(lambda2*input_data%omega_r) < input_data%asymplim .AND. &
          ABS(AIMAG(lambda2*input_data%omega_r)) < input_data%asymplim1) then
      
         F2n = hank1(lambda2*input_data%omega_r,input_data%azim_mode,1)
         F2d = dhank1(lambda2*input_data%omega_r,input_data%azim_mode,1)
         F2f = F2n/F2d
       
       else   ! asymtotic limit   

         F2f = (8._dpk*lambda2*input_data%omega_r +&
               4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%azim_mode*input_data%azim_mode - &
              CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*input_data%omega_r - &
              4._dpk*input_data%azim_mode*input_data%azim_mode - 3._dpk)
        
       end if
      
       Trs = F2f/lambda2

   end if


  END FUNCTION compute_Trs_lambda


  FUNCTION residuepot(r,z,ss,input_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Same as above but for computing velocity potential
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)            :: r, z
    integer              :: ss, ii, jj
    complex(dpk)         :: residuepot, res, fn
   
    type(input_params_t) :: input_data


    if (input_data%vs_param_gamma .EQ. 0) then

       residuepot = 0.

    else 

        res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M1)/((input_data%mu_plus - input_data%KH_zero_1)* & 
                                                                 input_data%k_minus_at_mu_plus*input_data%k_plus_sz1)

        if (input_data%vortswitch .EQ. 0) then
         
           if (r <= 1) then
              fn = (1._dpk - (input_data%KH_zero_1*input_data%M1))*compute_Trs_lambda(r,input_data%KH_zero_1,1,input_data)
           else
              fn = (1._dpk - (input_data%KH_zero_1*input_data%M2))*compute_Trs_lambda(r,input_data%KH_zero_1,2,input_data)
           end if
        end if
         
        residuepot = input_data%omega_r*res*fn*EXP(CMPLX(0.,1._dpk,kind=dpk)*input_data%omega_r*input_data%KH_zero_1*z)
        
    end if      
    
  END FUNCTION residuepot
 
  FUNCTION residuepr(r,z,ss,input_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Same as above but for computing velocity potential
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)            :: r, z
    integer              :: ss, ii, jj
    complex(dpk)         :: residuepr, res, fn
   
    type(input_params_t) :: input_data


    if (input_data%vs_param_gamma .EQ. 0) then

       residuepr = 0.
    

    else 

        res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M1)/((input_data%mu_plus - input_data%KH_zero_1)* & 
                                                                 input_data%k_minus_at_mu_plus*input_data%k_plus_sz1)

        if (input_data%vortswitch .EQ. 0) then
         
           if (r <= 1) then
              fn = ((1._dpk -(input_data%KH_zero_1*input_data%M1))**2)*compute_Trs_lambda(r,input_data%KH_zero_1,1,input_data)
           else
              fn = ((1._dpk -(input_data%KH_zero_1*input_data%M2))**2)*compute_Trs_lambda(r,input_data%KH_zero_1,2,input_data)
           end if
        end if
         
        residuepr = ((input_data%omega_r)**2)*res*fn*EXP(CMPLX(0.,1._dpk,kind=dpk)*input_data%omega_r*input_data%KH_zero_1*z)
        
    end if      
    
  END FUNCTION residuepr

 FUNCTION compute_psi_incident(r,z,stream_idx,input_data) result(psi0)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the \Psi_{mn}(r) of (2.14) or (4.11)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)        :: psi0, psimn, resp
    complex(dpk)        :: f1, f2, f3
    real(dpk)           :: r, z
    integer             :: stream_idx

    type(input_params_t) :: input_data

    if ((input_data%vortswitch == 1) .OR. (input_data%vortswitch == 2)) then
        
          print*,'Incident wave has to be an pipe mode of the nozzle  '
          STOP
      end if

    

    else

       if (stream_idx == 1) then
         psimn = bessj(input_data%alpha1*r,input_data%azim_mode,1)*EXP(ABS(AIMAG(input_data%alpha1*r)))
       
         psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*(1._dpk - input_data%M1*input_data%mu_plus)* &
              psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus*z)


       else
          psi0 = 0.

       end if

    end if


  END FUNCTION compute_psi_incident


  FUNCTION compute_Trs_lambda_instab(ri,si,ss,input_data)  result(Trsin)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute Trs of (4.8)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: ri
    complex(dpk)  :: sri, Trsin
    complex(dpk)  :: si
    complex(dpk)  :: l1, l2
    complex(dpk)  :: Fn, Fd, F
    integer       :: ss

    type(input_params_t) :: input_data


    l1 = sqrt(1._dpk - si*(input_data%M1+1._dpk))*sqrt(1._dpk - si*(input_data%M1-1._dpk))
    l2 = sqrt(1._dpk - si*(input_data%M2+1._dpk))*sqrt(1._dpk - si*(input_data%M2-1._dpk))


    if (ss==1) then
       
       Fn = bessj(l1*input_data%omega_r*ri,input_data%azim_mode,1)*EXP(ABS(AIMAG(l1*input_data%omega_r*ri)))
       Fd = dbessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1)*EXP(ABS(AIMAG(l1*input_data%omega_r*input_data%h)))
       
       F = Fn/Fd

       Trsin = F/l1

    else

       Fn = dhank1(l2*input_data%omega_r,input_data%azim_mode,1)*bessj(l2*input_data%omega_r*ri,input_data%azim_mode,1)* &
                                                                              EXP(ABS(AIMAG(l2*input_data%omega_r*ri))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r) - dbessj(l2*input_data%omega_r,input_data%azim_mode,1)* &
                                                  hank1(l2*input_data%omega_r*ri,input_data%azim_mode,1)* &
            EXP(ABS(AIMAG(l2*input_data%omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r*ri)

       Fd = dhank1(l2*input_data%omega_r,input_data%azim_mode,1)* &
                                   dbessj(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
                                              EXP(ABS(AIMAG(l2*input_data%omega_r*input_data%h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r) - dbessj(l2*input_data%omega_r,input_data%azim_mode,1)* &
                                                   dhank1(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            EXP(ABS(AIMAG(l2*input_data%omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r*input_data%h)
       
       F = Fn/Fd

       Trsin = F/l2

    end if
       

  END FUNCTION compute_Trs_lambda_instab

end module user_defined_functions


