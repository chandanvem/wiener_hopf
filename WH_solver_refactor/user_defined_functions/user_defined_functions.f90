Module user_defined_functions

  USE input_params
  USE bessel_utils
  

  IMPLICIT NONE

  PRIVATE :: bc, compute_Trs
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

    if (ss == 1) then
       u_s = (s_target-input_data%KH_pole_1)
    else
       if (input_data%vortswitch .EQ. 0) then
          u_s = (s_target-input_data%KH_zero_1)*(s_target-input_data%KH_zero_2)/(s_target-input_data%KH_pole_1)
       else
          u_s = (s_target-input_data%KH_zero_2)
       end if
    end if

    if (ss == 1) then
       do j = 1, input_data%num_sup_poles
          u_s = u_s*(s_target - input_data%sup_poles_list(j))
       end do
    else
       if (input_data%vortswitch .EQ. 0) then
          do j = 1,input_data%num_sup_zeros
             do jj = 1,input_data%num_sup_poles
                u_s = u_s*(s_target - input_data%sup_zeros_list(j))/(s_target - input_data%sup_poles_list(jj))
             end do
          end do
       end if
    end if

  END FUNCTION compute_U_s_factor

  FUNCTION compute_kernel(ss,zeta,input_data)  result(kernel)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the kernel function (3.22) or (4.10)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: kernel, zeta
    complex(dpk)    :: lambda_1, lambda_2, lambda_3, lambda_z
    complex(dpk)    :: F1n, F1d, F1, F1f, F2n, F2d, F2, F2f, Rn, Rd, Rz
    integer         :: ss
    type(input_params_t) :: input_data


    if (ss == 1) then  !! compute (4.10): the inc vort upstream kernel

       lambda_1 = sqrt(1._dpk - zeta*(input_data%M1+1._dpk))*sqrt(1._dpk - zeta*(input_data%M1-1._dpk))
       lambda_2 = sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M2+1._dpk))* &
                   sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M2-1._dpk))


       if ((ABS(lambda_1*input_data%omega_r*input_data%h) < input_data%asymplim .AND.&
                         ABS(AIMAG(lambda_1*input_data%omega_r*input_data%h)) < input_data%asymplim1)) then

          F1n = bessj(lambda_1*input_data%omega_r*input_data%h,input_data%azim_mode,1)
          F1d = dbessj(lambda_1*input_data%omega_r*input_data%h,input_data%azim_mode,1)
          F1f = F1n/F1d

       else

          F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
          
       end if

!!$    print*,'F1n:',F1n
!!$    print*,'F1d:',F1d
!!$    print*,'F1f:',F1f
!!$    print*,'lambda_z:',lambda_z

       F1 = input_data%kap_rho*(1._dpk - zeta*input_data%M1)**2/lambda_1*F1f

       if ((ABS(lambda_2*input_data%omega_r) < input_data%asymplim .AND.&
             ABS(AIMAG(lambda_2*input_data%omega_r)) < input_data%asymplim1) .AND. &
            (ABS(lambda_2*input_data%omega_r*input_data%h) < input_data%asymplim &
                      .AND. ABS(AIMAG(lambda_2*input_data%omega_r*input_data%h)) < input_data%asymplim1)) then

          F2n = dhank1(lambda_2*input_data%omega_r,input_data%azim_mode,1)*&
               bessj(lambda_2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
                EXP(ABS(AIMAG(lambda_2*input_data%omega_r*input_data%h))+ &
                CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*input_data%omega_r) - &
                dbessj(lambda_2*input_data%omega_r,input_data%azim_mode,1)* &
                 hank1(lambda_2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* & 
                EXP(ABS(AIMAG(lambda_2*input_data%omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*input_data%omega_r*input_data%h)

          F2d = dhank1(lambda_2*input_data%omega_r,input_data%azim_mode,1)* &
                dbessj(lambda_2*input_data%omega_r*input_data%h,input_data%azim_mode,1)*&
                EXP(ABS(AIMAG(lambda_2*input_data%omega_r*input_data%h))+ &
                CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*input_data%omega_r)- &
                dbessj(lambda_2*input_data%omega_r,input_data%azim_mode,1)*&
                dhank1(lambda_2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
                EXP(ABS(AIMAG(lambda_2*input_data%omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*input_data%omega_r*input_data%h)
          F2f = F2n/F2d

       else
          
          F2f = CMPLX(0._dpk,1._dpk,kind=dpk)

       end if

       F2 = (1._dpk - zeta*input_data%M2)**2/lambda_2*F2f

       kernel = input_data%omega_r*(F1 - F2)

    else  !! compute (3.22): the kernel

       lambda_1 = sqrt(1._dpk - zeta*(input_data%M1+1._dpk))*sqrt(1._dpk - zeta*(input_data%M1-1._dpk))

       lambda_2 = sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M2+1._dpk))* &
                  sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M2-1._dpk))

       lambda_3 = sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M3+1._dpk))* &
                  sqrt(input_data%kapT - zeta*(input_data%kapT*input_data%M3-1._dpk))

       lambda_z = input_data%kap_rho*lambda_2/lambda_1*(1._dpk - zeta*input_data%M1)**2/&
                  (1._dpk - zeta*input_data%M2)**2


       if ((ABS(lambda_2*input_data%omega_r) < input_data%asymplim .AND.&
              ABS(AIMAG(lambda_2*input_data%omega_r)) < input_data%asymplim1) .AND. &
            (ABS(lambda_2*input_data%omega_r*input_data%h) < input_data%asymplim .AND.&
                ABS(AIMAG(lambda_2*input_data%omega_r*input_data%h)) < input_data%asymplim1) .AND. &
            (ABS(lambda_1*input_data%omega_r*input_data%h) < input_data%asymplim .AND. &
                ABS(AIMAG(lambda_1*input_data%omega_r*input_data%h)) < input_data%asymplim1)) then

          Rd = (lambda_z*bessj(lambda_1*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
                 dhank1(lambda_2*input_data%omega_r*input_data%h,input_data%azim_mode,1)- &
               hank1(lambda_2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
               dbessj(lambda_1*input_data%omega_r*input_data%h,input_data%azim_mode,1))* &
               EXP(ABS(AIMAG(lambda_1*input_data%omega_r*input_data%h)) &
                +CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*input_data%omega_r*input_data%h)

          Rn = (bessj(lambda_2*input_data%omega_r*input_data%h,input_data%azim_mode,1)*&
                dbessj(lambda_1*input_data%omega_r*input_data%h,input_data%azim_mode,1)- &
               lambda_z*bessj(lambda_1*input_data%omega_r*input_data%h,input_data%azim_mode,1)*&
                dbessj(lambda_2*input_data%omega_r*input_data%h,input_data%azim_mode,1))* &
               EXP(ABS(AIMAG(lambda_1*input_data%omega_r*input_data%h))+ABS(AIMAG(lambda_2*input_data%omega_r*input_data%h)))

          Rz = Rn/Rd


          F1n = Rz*hank1(lambda_2*input_data%omega_r,input_data%azim_mode,1)*&
              EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*input_data%omega_r)+ &
               bessj(lambda_2*input_data%omega_r,input_data%azim_mode,1)*EXP(ABS(AIMAG(lambda_2*input_data%omega_r)))

          F1d = Rz*dhank1(lambda_2*input_data%omega_r,input_data%azim_mode,1)* &
                  EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*input_data%omega_r)+ &
               dbessj(lambda_2*input_data%omega_r,input_data%azim_mode,1)*EXP(ABS(AIMAG(lambda_2*input_data%omega_r)))
          F1f = F1n/F1d

       else

          F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
       
       end if

!!$    print*,'F1n:',F1n
!!$    print*,'F1d:',F1d
!!$    print*,'F1f:',F1f
!!$    print*,'lambda_z:',lambda_z

       F1 = (1._dpk - zeta*input_data%M2)**2/lambda_2*F1f

       if (ABS(lambda_3*input_data%omega_r) < input_data%asymplim .AND. &
            ABS(AIMAG(lambda_3*input_data%omega_r)) < input_data%asymplim1) then
          
          F2n = hank1(lambda_3*input_data%omega_r,input_data%azim_mode,1)
          F2d = dhank1(lambda_3*input_data%omega_r,input_data%azim_mode,1)
          F2f = F2n/F2d
          
       else

          F2f = (8._dpk*lambda_3*input_data%omega_r +&
                4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%azim_mode*input_data%azim_mode - &
               CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_3*input_data%omega_r - &
               4._dpk*input_data%azim_mode*input_data%azim_mode - 3._dpk)

       end if

!!$    print*,'F2n:',F2n
!!$    print*,'F2d:',F2d
!!$    print*,'F2f:',F2f

       F2 = (1._dpk - zeta*input_data%M3)**2/lambda_3*F2f

       if (input_data%num_zeros_s1_s2 == 0 .AND. input_data%num_poles_s1_s2 == 0) then
          
          kernel = input_data%omega_r*(F1 - F2)
 
       else

          kernel = input_data%omega_r*(F1 - F2)*bc(zeta,input_data)

       end if

    end if


  END FUNCTION compute_kernel

  FUNCTION bc(z,input_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the term to be factored out from the kernel function (3.22)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)  :: z, bc
    integer       :: i
    type(input_params_t) :: input_data

    bc = 1._dpk

    do i = 1, input_data%num_zeros_s1_s2

       bc = bc/(z - input_data%zeros_list_bw_s1_s2(i))

    end do

    do i = 1, input_data%num_poles_s1_s2

       bc = bc*(z - input_data%poles_list_bw_s1_s2(i))

    end do


  END FUNCTION bc

  FUNCTION integrand_IFT_pot(ri,zi,i,ss,input_data,contour_data) result(integrandiftpot)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the IFT integrand when computing for potential
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)       :: ri, zi
    integer         :: i, ss
    complex(dpk)    :: u
    complex(dpk)    :: integrandiftpot

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data

    u = contour_data%iftpoints(i)

    if (ss==1) then

!! potential:

       integrandiftpot = ((1._dpk - u*input_data%M2)**2)/(1._dpk - u*input_data%M1)
       integrandiftpot = integrandiftpot *compute_Trs(ri,u,1,input_data)*input_data%fplusz(i)* & 
                           EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)

    else if (ss==2) then

!! potential:
       
       integrandiftpot = (1._dpk - u*input_data%M2)
       integrandiftpot = integrandiftpot*compute_Trs(ri,u,2,input_data)*input_data%fplusz(i)* &
                         EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)

    else

!! potential:

       integrandiftpot = (1._dpk - u*input_data%M3) 
       integrandiftpot = integrandiftpot*compute_Trs(ri,u,3,input_data)*input_data%fplusz(i)* &
                         EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)

    end if


  END FUNCTION integrand_IFT_pot


  FUNCTION integrand_IFT_pr(ri,zi,i,ss,input_data,contour_data) result(integrandiftpr)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the IFT integrand when computing for pressure
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)       :: ri, zi
    integer         :: i, ss
    complex(dpk)    :: u
    complex(dpk)    :: integrandiftpr

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data

    u = contour_data%iftpoints(i)


    if (ss==1) then

!! pressure:

       integrandiftpr = (1._dpk - u*input_data%M2)**2
       integrandiftpr = integrandiftpr*compute_Trs(ri,u,1,input_data)*input_data%fplusz(i)* & 
                        EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)
    else if (ss==2) then

!! pressure:
       
       integrandiftpr = (1._dpk - u*input_data%M2)**2
       integrandiftpr = integrandiftpr*compute_Trs(ri,u,2,input_data)*input_data%fplusz(i)* &
                        EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)
    else

!! pressure:

       integrandiftpr = (1._dpk - u*input_data%M3)**2
       integrandiftpr = integrandiftpr*compute_Trs(ri,u,3,input_data)*input_data%fplusz(i)* &
                        EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*u*zi)
    end if


  END FUNCTION integrand_IFT_pr

!
  FUNCTION compute_Trs(ri,si,ss,input_data) result(Trs)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute Trs of (3.33)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: ri
    complex(dpk)  :: sri, Trs
    complex(dpk)  :: si
    complex(dpk)  :: l1, l2, l3, lz
    complex(dpk)  :: Fn, Fd, F, Rn, Rd, Rz
    integer       :: ss
  
    type(input_params_t)  :: input_data

    l1 = sqrt(1._dpk - si*(input_data%M1+1._dpk))*sqrt(1._dpk - si*(input_data%M1-1._dpk))
    l2 = sqrt(input_data%kapT - si*(input_data%kapT*input_data%M2+1._dpk))* &
         sqrt(input_data%kapT - si*(input_data%kapT*input_data%M2-1._dpk))
    l3 = sqrt(input_data%kapT - si*(input_data%kapT*input_data%M3+1._dpk))* &
         sqrt(input_data%kapT - si*(input_data%kapT*input_data%M3-1._dpk))

    lz = input_data%kap_rho*l2/l1
    lz = lz*((1._dpk - si*input_data%M1)**2)/((1._dpk - si*input_data%M2)**2)


    if (ss==1) then
       
       Rd = (lz*bessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            dhank1(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)- &
            hank1(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            dbessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1))* &
            EXP(ABS(AIMAG(l1*input_data%omega_r*input_data%h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r*input_data%h)
       
       Rn = (bessj(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            dbessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1)- &
            lz*bessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1)*&
            dbessj(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1))* &
            EXP(ABS(AIMAG(l1*input_data%omega_r*input_data%h))+ABS(AIMAG(l2*input_data%omega_r*input_data%h)))
       
       Rz = Rn/Rd
          
       Fn = bessj(l1*input_data%omega_r*ri,input_data%azim_mode,1)*EXP(ABS(AIMAG(l1*input_data%omega_r*ri)))* &
            (Rz*hank1(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r*input_data%h) +&
            bessj(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            EXP(ABS(AIMAG(l2*input_data%omega_r*input_data%h))))
       
       Fd = bessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            EXP(ABS(AIMAG(l1*input_data%omega_r*input_data%h)))*(Rz*dhank1(l2*input_data%omega_r,input_data%azim_mode,1)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r) + dbessj(l2*input_data%omega_r,input_data%azim_mode,1)* &
            EXP(ABS(AIMAG(l2*input_data%omega_r))))
       
       F = Fn/Fd

       Trs = F/l2

    else if (ss==2) then

       Rd = (lz*bessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
             dhank1(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)- &
            hank1(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            dbessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1))* &
            EXP(ABS(AIMAG(l1*input_data%omega_r*input_data%h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r*input_data%h)
       
       Rn = (bessj(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            dbessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1)- &
            lz*bessj(l1*input_data%omega_r*input_data%h,input_data%azim_mode,1)* &
            dbessj(l2*input_data%omega_r*input_data%h,input_data%azim_mode,1))* &
            EXP(ABS(AIMAG(l1*input_data%omega_r*input_data%h))+ABS(AIMAG(l2*input_data%omega_r*input_data%h)))
       
       Rz = Rn/Rd
          
       Fn = Rz*hank1(l2*input_data%omega_r*ri,input_data%azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r*ri)+ &
            bessj(l2*input_data%omega_r*ri,input_data%azim_mode,1)*EXP(ABS(AIMAG(l2*input_data%omega_r*ri)))
       
       Fd = Rz*dhank1(l2*input_data%omega_r,input_data%azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*input_data%omega_r)+ &
            dbessj(l2*input_data%omega_r,input_data%azim_mode,1)*EXP(ABS(AIMAG(l2*input_data%omega_r)))
       
       F = Fn/Fd

       Trs = F/l2

    else

!!$       if ((ABS(l3*omega_r) < asymplim .AND. ABS(AIMAG(l3*omega_r)) < asymplim1) .AND. &
!!$            (ABS(l3*omega_r*ri) < asymplim .AND. ABS(AIMAG(l3*omega_r*ri)) < asymplim1)) then

          Fn = hank1(l3*input_data%omega_r*ri,input_data%azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l3*input_data%omega_r*ri)
       
          Fd = dhank1(l3*input_data%omega_r,input_data%azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l3*input_data%omega_r)
       
          F = Fn/Fd

!!$       else 
!!$       
!!$          F = (8._dpk*l3*omega_r*ri + 4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*azim_mode*azim_mode - &
!!$               CMPLX(0._dpk,1._dpk,kind=dpk))/((8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*l3*omega_r - &
!!$               4._dpk*azim_mode*azim_mode - 3._dpk)*ri**(1.5_dpk))*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
!!$               l3*omega_r*(ri - 1._dpk))
!!$
!!$       end if
  
       Trs = F/l3

    end if
       
!    print*,'Rz:',Rz


  END FUNCTION compute_Trs


  FUNCTION residuepot(r,z,ss,input_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Same as above but for computing velocity potential
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)            :: r, z
    integer              :: ss, ii, jj
    complex(dpk)         :: residuepot, res, fn
   
    type(input_params_t) :: input_data

    if (ss == 1) then

      res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M2)/((input_data%mu_plus - input_data%KH_zero_1)* & 
             input_data%k_minus_at_mu_plus*input_data%k_plus_sz1*(input_data%KH_zero_1 - input_data%KH_zero_2))

       if (input_data%vortswitch .EQ. 0) then
          do ii = 1, input_data%num_sup_zeros
             res = res/(input_data%KH_zero_1 - input_data%sup_zeros_list(ii))
          end do
          do jj = 1, input_data%num_sup_poles
             res = res*(input_data%KH_zero_1 - input_data%sup_poles_list(jj))
          end do
       end if

       if (r <= input_data%h) then
          fn = (1._dpk - input_data%KH_zero_1*input_data%M2)**2/(1._dpk - input_data%KH_zero_1*input_data%M1)*&
                                                                        compute_Trs(r,input_data%KH_zero_1,1,input_data)

       else if (r <= 1. .AND. r > input_data%h) then
          fn = (1._dpk - input_data%KH_zero_1*input_data%M2)*compute_Trs(r,input_data%KH_zero_1,2,input_data)
       
       else
          fn = (1._dpk - input_data%KH_zero_1*input_data%M3)*compute_Trs(r,input_data%KH_zero_1,3,input_data)
       end if
       
       residuepot = input_data%omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*input_data%omega_r*input_data%KH_zero_1*z)*res*fn
          
!!$       else
!!$          residuepot = 0.
!!$
!!$       end if
          
    elseif (ss == 2) then

       res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M2)/((input_data%mu_plus - input_data%KH_zero_2)* &
                            input_data%k_minus_at_mu_plus*input_data%k_plus_sz2*(input_data%KH_zero_2 - input_data%KH_zero_1))

       if (input_data%vortswitch .EQ. 0) then
          do ii = 1, input_data%num_sup_zeros
             res = res/(input_data%KH_zero_1 - input_data%sup_zeros_list(ii))
          end do
          do jj = 1, input_data%num_sup_poles
             res = res*(input_data%KH_zero_1 - input_data%sup_poles_list(jj))
          end do
       end if
          
       if (r <= input_data%h) then
          fn = (1._dpk - input_data%KH_zero_2*input_data%M2)**2/ &
               (1._dpk - input_data%KH_zero_2*input_data%M1)*compute_Trs(r,input_data%KH_zero_2,1,input_data)

       else if (r <= 1. .AND. r > input_data%h) then
          fn = (1._dpk - input_data%KH_zero_2*input_data%M2)*compute_Trs(r,input_data%KH_zero_2,2,input_data)
       else
          fn = (1._dpk - input_data%KH_zero_2*input_data%M3)*compute_Trs(r,input_data%KH_zero_2,3,input_data)
       end if

       residuepot = input_data%omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*input_data%omega_r*input_data%KH_zero_2*z)*res*fn

!!$       else
!!$          residuepot = 0.
!!$
!!$       end if

    else

       res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M2)/((input_data%mu_plus - input_data%sup_zeros_list(ss-2))* &
            input_data%k_minus_at_mu_plus*input_data%kzsp(ss-2)*(input_data%sup_zeros_list(ss-2) - input_data%KH_zero_1)* &
               (input_data%sup_zeros_list(ss-2) - input_data%KH_zero_2))
       do ii = 1, input_data%num_sup_zeros
          if (ii .NE. (ss-2)) then
             res = res/(input_data%sup_zeros_list(ss-2) - input_data%sup_zeros_list(ii))
          end if
       end do
       do jj = 1, input_data%num_sup_poles
          res = res*(input_data%sup_zeros_list(ss-2) - input_data%sup_poles_list(jj))
       end do
          
       if (r <= input_data%h) then
          fn = ((1._dpk - input_data%sup_zeros_list(ss-2)*input_data%M2)**2)/ &
                                        (1._dpk - input_data%sup_zeros_list(ss-2)*input_data%M1)* &
                                                    compute_Trs(r,input_data%sup_zeros_list(ss-2),1,input_data)
       else if (r <= 1. .AND. r > input_data%h) then
          fn = (1._dpk - input_data%sup_zeros_list(ss-2)*input_data%M2)*compute_Trs(r,input_data%sup_zeros_list(ss-2),2,input_data)
       else
          fn = (1._dpk - input_data%sup_zeros_list(ss-2)*input_data%M3)*compute_Trs(r,input_data%sup_zeros_list(ss-2),3,input_data)
       end if

       residuepot = input_data%omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*input_data%omega_r*input_data%sup_zeros_list(ss-2)*z)*res*fn

!!$       else
!!$          residuepot = 0.
!!$
!!$       end if

    end if

    
  END FUNCTION residuepot

 
 FUNCTION residuepr(r,z,ss,input_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the residue pressure term of (3.38)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)            :: r, z
    integer              :: ss, ii, jj
    complex(dpk)         :: residuepr, res, fn
   
    type(input_params_t) :: input_data

    if (ss == 1) then

       res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M2)/((input_data%mu_plus-input_data%KH_zero_1)* &
                             input_data%k_minus_at_mu_plus*input_data%k_plus_sz1*(input_data%KH_zero_1 - input_data%KH_zero_2))


      if (input_data%vortswitch .EQ. 0) then
          do ii = 1, input_data%num_sup_zeros
             res = res/(input_data%KH_zero_1 - input_data%sup_zeros_list(ii))
          end do
          do jj = 1, input_data%num_sup_poles
             res = res*(input_data%KH_zero_1 - input_data%sup_poles_list(jj))
          end do
       end if

       if (r <= input_data%h) then
          fn = (1._dpk - input_data%KH_zero_1*input_data%M2)**2*compute_Trs(r,input_data%KH_zero_1,1,input_data)
       else if (r <= 1. .AND. r > input_data%h) then
          fn = (1._dpk - input_data%KH_zero_1*input_data%M2)**2*compute_Trs(r,input_data%KH_zero_1,2,input_data)
       else
          fn = (1._dpk - input_data%KH_zero_1*input_data%M3)**2*compute_Trs(r,input_data%KH_zero_1,3,input_data)
       end if
       
       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            input_data%omega_r*input_data%KH_zero_1*z)*res*fn
          
         
    elseif (ss == 2) then


      res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M2)/((input_data%mu_plus - input_data%KH_zero_2)* &
                    input_data%k_minus_at_mu_plus*input_data%k_plus_sz2*(input_data%KH_zero_2 - input_data%KH_zero_1))

       if (input_data%vortswitch .EQ. 0) then

          do ii = 1, input_data%num_sup_zeros
             res = res/(input_data%KH_zero_2 - input_data%sup_zeros_list(ii))
          end do
          do jj = 1, input_data%num_sup_poles
             res = res*(input_data%KH_zero_2 - input_data%sup_poles_list(jj))
          end do
       end if
          
       if (r <= input_data%h) then
          fn = (1._dpk - input_data%KH_zero_2*input_data%M2)**2*compute_Trs(r,input_data%KH_zero_2,1,input_data)
       else if (r <= 1. .AND. r > input_data%h) then
          fn = (1._dpk - input_data%KH_zero_2*input_data%M2)**2*compute_Trs(r,input_data%KH_zero_2,2,input_data)
       else
          fn = (1._dpk - input_data%KH_zero_2*input_data%M3)**2*compute_Trs(r,input_data%KH_zero_2,3,input_data)
       end if

       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            input_data%omega_r*input_data%KH_zero_2*z)*res*fn

    else

      res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M2)/((input_data%mu_plus - input_data%sup_zeros_list(ss-2))* &
            input_data%k_minus_at_mu_plus*input_data%kzsp(ss-2)*(input_data%sup_zeros_list(ss-2) -input_data%KH_zero_1)* &
              (input_data%sup_zeros_list(ss-2) - input_data%KH_zero_2)) 


      do ii = 1, input_data%num_sup_zeros
          if (ii .NE. (ss-2)) then
             res = res/(input_data%sup_zeros_list(ss-2) - input_data%sup_zeros_list(ii))
          end if
       end do

       do jj = 1, input_data%num_sup_poles
          res = res*(input_data%sup_zeros_list(ss-2) - input_data%sup_poles_list(jj))
       end do
          
       if (r <= input_data%h) then
          fn =( (1._dpk - input_data%sup_zeros_list(ss-2)*input_data%M2)**2) * &
                             compute_Trs(r,input_data%sup_zeros_list(ss-2),1,input_data)
       else if (r <= 1. .AND. r > input_data%h) then
          fn = ((1._dpk - input_data%sup_zeros_list(ss-2)*input_data%M2)**2) * &
                          compute_Trs(r,input_data%sup_zeros_list(ss-2),2,input_data)
       else
          fn = ((1._dpk - input_data%sup_zeros_list(ss-2)*input_data%M3)**2)* &
                         compute_Trs(r,input_data%sup_zeros_list(ss-2),3,input_data)
       end if

       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            input_data%omega_r*input_data%sup_zeros_list(ss-2)*z)*res*fn

   end if

    
  END FUNCTION residuepr


 FUNCTION compute_psi_incident(r,z,ss,input_data) result(psi0)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the \Psi_{mn}(r) of (2.14) or (4.11)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)        :: psi0, psimn, resp
    complex(dpk)        :: f1, f2, f3
    real(dpk)           :: r, z
    integer             :: ss

    type(input_params_t) :: input_data

    if ((input_data%vortswitch == 1) .OR. (input_data%vortswitch == 2)) then
        
       if (z >= input_data%Zo) then 
          if (ss == 1) then
                 psimn = input_data%omega_r*resp* &
                         (1._dpk - input_data%M1*input_data%mu_plus)*compute_Trs_instab(r,input_data%mu_plus,1,input_data)

                 psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)* &
                        input_data%omega_r*(1._dpk - input_data%M1*input_data%mu_plus)* &
                        psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus*(z-input_data%Zo))

           else if (ss == 2) then
     
                psimn = input_data%omega_r*resp* &
                       (1._dpk - input_data%M2*input_data%mu_plus)*compute_Trs_instab(r,input_data%mu_plus,2,input_data)

                psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)* &
                       input_data%omega_r*(1._dpk - input_data%M2*input_data%mu_plus)* &
                       psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus*(z-input_data%Zo))

           else

              psi0 = 0.

           end if

       else

          psi0 = 0.

       end if

    else

       if (ss == 1) then
         psimn = bessj(input_data%alpha1*r,input_data%azim_mode,1)*EXP(ABS(AIMAG(input_data%alpha1*r)))
       
         psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*(1._dpk - input_data%M1*input_data%mu_plus)* &
              psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus*z)


       else if (ss == 2) then

          f1 = (1._dpk - input_data%mu_plus*input_data%M1)/(1._dpk - input_data%mu_plus*input_data%M2)* &
                                         bessj(input_data%alpha1*input_data%h,input_data%azim_mode,1)
       
          f2 = bessj(input_data%alpha2*r,input_data%azim_mode,1)* &
              dhank1(input_data%alpha2,input_data%azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* & 
              input_data%alpha2 + ABS(AIMAG(input_data%alpha2*r)))- & 
              hank1(input_data%alpha2*r,input_data%azim_mode,1)* &
              dbessj(input_data%alpha2,input_data%azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* & 
              input_data%alpha2*r + ABS(AIMAG(input_data%alpha2)))
  

          f3 = bessj(input_data%alpha2*input_data%h,input_data%azim_mode,1)* &
               dhank1(input_data%alpha2,input_data%azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
               input_data%alpha2 + ABS(AIMAG(input_data%alpha2*input_data%h)))- & 
               hank1(input_data%alpha2*input_data%h,input_data%azim_mode,1)* &
               dbessj(input_data%alpha2,input_data%azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
               input_data%alpha2*input_data%h + ABS(AIMAG(input_data%alpha2)))
            

          psimn = input_data%kap_rho*f1*f2/f3

          psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*(1._dpk - input_data%M2*input_data%mu_plus)* &
                psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus*z)

!!$       open(10,file='test.out',form='FORMATTED',position='APPEND')
!!$       write(10,'(2E20.10)') psimn
!!$       close(10)

       else
          psi0 = 0.

       end if

    end if


  END FUNCTION compute_psi_incident


  FUNCTION compute_Trs_instab(ri,si,ss,input_data)  result(Trsin)

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
       
!    print*,'Rz:',Rz


  END FUNCTION compute_Trs_instab

!
    
!
!
end module user_defined_functions


