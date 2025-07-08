Module user_defined_functions

  USE input_params
  USE bessel_utils
  

  IMPLICIT NONE

  PRIVATE ::
  PUBLIC  :: compute_U_s_factor, compute_kernel

  FUNCTION compute_U_s_factor(ss,s_target) result(u_s)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The factor U(s)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)  :: s_target, u_s
    integer       :: j, jj, ss

    if (ss == 1) then
       u_s = (s_target-KH_pole_1)
    else
       if (vortswitch .EQ. 0) then
          u_s = (s_target-KH_zero_1)*(s_target-KH_zero_2)/(s_target-KH_pole_1)
       else
          u_s = (s_target-KH_zero_2)
       end if
    end if

    if (ss == 1) then
       do j = 1, num_sup_poles
          u_s = u_s*(s_target - sup_poles_list(j))
       end do
    else
       if (vortswitch .EQ. 0) then
          do j = 1, num_sup_zeros
             do jj = 1, num_sup_poles
                u_s = u_s*(s_target - sup_zeros_list(j))/(s_target - sup_poles_list(jj))
             end do
          end do
       end if
    end if

  END FUNCTION compute_U_s_factor

  FUNCTION compute_kernel(ss,zeta)  result(kernel)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the kernel function (3.22) or (4.10)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: kernel, zeta
    complex(dpk)    :: lambda_1, lambda_2, lambda_3, lambda_z
    complex(dpk)    :: F1n, F1d, F1, F1f, F2n, F2d, F2, F2f, Rn, Rd, Rz
    integer         :: ss

    if (ss == 1) then  !! compute (4.10): the inc vort upstream kernel

       lambda_1 = sqrt(1._dpk - zeta*(M1+1._dpk))*sqrt(1._dpk - zeta*(M1-1._dpk))
       lambda_2 = sqrt(kapT - zeta*(kapT*M2+1._dpk))*sqrt(kapT - zeta*(kapT*M2-1._dpk))


       if ((ABS(lambda_1*omega_r*h) < asymplim .AND. ABS(AIMAG(lambda_1*omega_r*h)) < asymplim1)) then

          F1n = bessj(lambda_1*omega_r*h,azim_mode,1)
          F1d = dbessj(lambda_1*omega_r*h,azim_mode,1)
          F1f = F1n/F1d

       else

          F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
          
       end if

!!$    print*,'F1n:',F1n
!!$    print*,'F1d:',F1d
!!$    print*,'F1f:',F1f
!!$    print*,'lambda_z:',lambda_z

       F1 = kap_rho*(1._dpk - zeta*M1)**2/lambda_1*F1f

       if ((ABS(lambda_2*omega_r) < asymplim .AND. ABS(AIMAG(lambda_2*omega_r)) < asymplim1) .AND. &
            (ABS(lambda_2*omega_r*h) < asymplim .AND. ABS(AIMAG(lambda_2*omega_r*h)) < asymplim1)) then

          F2n = dhank1(lambda_2*omega_r,azim_mode,1)*bessj(lambda_2*omega_r*h,azim_mode,1)* &
                EXP(ABS(AIMAG(lambda_2*omega_r*h))+ &
                CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r) - &
                dbessj(lambda_2*omega_r,azim_mode,1)*hank1(lambda_2*omega_r*h,azim_mode,1)* & 
                EXP(ABS(AIMAG(lambda_2*omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r*h)

          F2d = dhank1(lambda_2*omega_r,azim_mode,1)*dbessj(lambda_2*omega_r*h,azim_mode,1)*&
                EXP(ABS(AIMAG(lambda_2*omega_r*h))+ &
                CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r) -&
                dbessj(lambda_2*omega_r,azim_mode,1)*dhank1(lambda_2*omega_r*h,azim_mode,1)* &
                EXP(ABS(AIMAG(lambda_2*omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r*h)
          F2f = F2n/F2d

       else
          
          F2f = CMPLX(0._dpk,1._dpk,kind=dpk)

       end if

       F2 = (1._dpk - zeta*M2)**2/lambda_2*F2f

       kernel = omega_r*(F1 - F2)

    else  !! compute (3.22): the kernel

       lambda_1 = sqrt(1._dpk - zeta*(M1+1._dpk))*sqrt(1._dpk - zeta*(M1-1._dpk))
       lambda_2 = sqrt(kapT - zeta*(kapT*M2+1._dpk))*sqrt(kapT - zeta*(kapT*M2-1._dpk))
       lambda_3 = sqrt(kapT - zeta*(kapT*M3+1._dpk))*sqrt(kapT - zeta*(kapT*M3-1._dpk))

       lambda_z = kap_rho*lambda_2/lambda_1*(1._dpk - zeta*M1)**2/(1._dpk - zeta*M2)**2


       if ((ABS(lambda_2*omega_r) < asymplim .AND. ABS(AIMAG(lambda_2*omega_r)) < asymplim1) .AND. &
            (ABS(lambda_2*omega_r*h) < asymplim .AND. ABS(AIMAG(lambda_2*omega_r*h)) < asymplim1) .AND. &
            (ABS(lambda_1*omega_r*h) < asymplim .AND. ABS(AIMAG(lambda_1*omega_r*h)) < asymplim1)) then

          Rd = (lambda_z*bessj(lambda_1*omega_r*h,azim_mode,1)*dhank1(lambda_2*omega_r*h,azim_mode,1)- &
               hank1(lambda_2*omega_r*h,azim_mode,1)*dbessj(lambda_1*omega_r*h,azim_mode,1))* &
               EXP(ABS(AIMAG(lambda_1*omega_r*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r*h)

          Rn = (bessj(lambda_2*omega_r*h,azim_mode,1)*dbessj(lambda_1*omega_r*h,azim_mode,1)- &
               lambda_z*bessj(lambda_1*omega_r*h,azim_mode,1)*dbessj(lambda_2*omega_r*h,azim_mode,1))* &
               EXP(ABS(AIMAG(lambda_1*omega_r*h))+ABS(AIMAG(lambda_2*omega_r*h)))

          Rz = Rn/Rd


          F1n = Rz*hank1(lambda_2*omega_r,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r)+ &
               bessj(lambda_2*omega_r,azim_mode,1)*EXP(ABS(AIMAG(lambda_2*omega_r)))
          F1d = Rz*dhank1(lambda_2*omega_r,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r)+ &
               dbessj(lambda_2*omega_r,azim_mode,1)*EXP(ABS(AIMAG(lambda_2*omega_r)))
          F1f = F1n/F1d

       else

          F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
       
       end if

!!$    print*,'F1n:',F1n
!!$    print*,'F1d:',F1d
!!$    print*,'F1f:',F1f
!!$    print*,'lambda_z:',lambda_z

       F1 = (1._dpk - zeta*M2)**2/lambda_2*F1f

       if (ABS(lambda_3*omega_r) < asymplim .AND. ABS(AIMAG(lambda_3*omega_r)) < asymplim1) then
          
          F2n = hank1(lambda_3*omega_r,azim_mode,1)
          F2d = dhank1(lambda_3*omega_r,azim_mode,1)
          F2f = F2n/F2d
          
       else

          F2f = (8._dpk*lambda_3*omega_r + 4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*azim_mode*azim_mode - &
               CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_3*omega_r - &
               4._dpk*azim_mode*azim_mode - 3._dpk)

       end if

!!$    print*,'F2n:',F2n
!!$    print*,'F2d:',F2d
!!$    print*,'F2f:',F2f

       F2 = (1._dpk - zeta*M3)**2/lambda_3*F2f

       if (num_zeros_s1_s2 == 0 .AND. num_poles_s1_s2 == 0) then
          
          kernel = omega_r*(F1 - F2)
 
       else

          kernel = omega_r*(F1 - F2)*bc(zeta)

       end if

    end if


  END FUNCTION compute_kernel


!  FUNCTION integrandiftpr(ri,zi,i,ss)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Compute the IFT integrand when computing for pressure
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    real(dpk)       :: ri, zi
!    integer         :: i, ss
!    complex(dpk)    :: u
!    complex(dpk)    :: integrandiftpr
!
!
!    u = iftpoints(i)
!
!    if (ss==1) then
!
!!! pressure:
!
!       integrandiftpr = (1._dpk - u*M2)**2*compute_Trs(ri,u,1)*fplusz(i)* & 
!            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!!!$       integrandiftpr = compute_Trs(ri,u,1)
!    else if (ss==2) then
!
!!! pressure:
!       
!       integrandiftpr = (1._dpk - u*M2)**2*compute_Trs(ri,u,2)*fplusz(i)* &
!            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!!!$       integrandiftpr = compute_Trs(ri,u,2)
!    else
!
!!! pressure:
!
!       integrandiftpr = (1._dpk - u*M3)**2*compute_Trs(ri,u,3)*fplusz(i)* &
!            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!!!$       integrandiftpr =  compute_Trs(ri,u,3)
!    end if
!
!
!  END FUNCTION integrandiftpr
!
!
!  FUNCTION integrand_IFT_pot(ri,zi,i,ss) result(integrandiftpot)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Compute the IFT integrand when computing for potential
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    real(dpk)       :: ri, zi
!    integer         :: i, ss
!    complex(dpk)    :: u
!    complex(dpk)    :: integrandiftpot
!
!
!    u = iftpoints(i)
!
!    if (ss==1) then
!
!!! potential:
!
!       integrandiftpot = (1._dpk - u*M2)**2/(1._dpk - u*M1)*compute_Trs(ri,u,1)*fplusz(i)* & 
!            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!
!    else if (ss==2) then
!
!!! potential:
!       
!       integrandiftpot = (1._dpk - u*M2)*compute_Trs(ri,u,2)*fplusz(i)* &
!            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!
!    else
!
!!! potential:
!
!       integrandiftpot = (1._dpk - u*M3)*compute_Trs(ri,u,3)*fplusz(i)* &
!            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!
!    end if
!
!
!  END FUNCTION integrand_IFT_pot
!
!   FUNCTION bc(z)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Compute the term to be factored out from the kernel function (3.22)
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    complex(dpk)  :: z, bc
!    integer       :: i
!
!    bc = 1._dpk
!
!    do i = 1, num_zeros_s1_s2
!
!       bc = bc/(z - zeros_list_bw_s1_s2(i))
!
!    end do
!
!    do i = 1, num_poles_s1_s2
!
!       bc = bc*(z - poles_list_bw_s1_s2(i))
!
!    end do
!
!
!  END FUNCTION bc
!
!
!  FUNCTION compute_Trs(ri,si,ss) result(Trs)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Compute Trs of (3.33)
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    real(dpk)     :: ri
!    complex(dpk)  :: sri, Trs
!    complex(dpk)  :: si
!    complex(dpk)  :: l1, l2, l3, lz
!    complex(dpk)  :: Fn, Fd, F, Rn, Rd, Rz
!    integer       :: ss
!
!
!    l1 = sqrt(1._dpk - si*(M1+1._dpk))*sqrt(1._dpk - si*(M1-1._dpk))
!    l2 = sqrt(kapT - si*(kapT*M2+1._dpk))*sqrt(kapT - si*(kapT*M2-1._dpk))
!    l3 = sqrt(kapT - si*(kapT*M3+1._dpk))*sqrt(kapT - si*(kapT*M3-1._dpk))
!
!    lz = kap_rho*l2/l1*(1._dpk - si*M1)**2/(1._dpk - si*M2)**2
!
!
!    if (ss==1) then
!       
!       Rd = (lz*bessj(l1*omega_r*h,azim_mode,1)*dhank1(l2*omega_r*h,azim_mode,1)- &
!            hank1(l2*omega_r*h,azim_mode,1)*dbessj(l1*omega_r*h,azim_mode,1))* &
!            EXP(ABS(AIMAG(l1*omega_r*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*h)
!       
!       Rn = (bessj(l2*omega_r*h,azim_mode,1)*dbessj(l1*omega_r*h,azim_mode,1)- &
!            lz*bessj(l1*omega_r*h,azim_mode,1)*dbessj(l2*omega_r*h,azim_mode,1))* &
!            EXP(ABS(AIMAG(l1*omega_r*h))+ABS(AIMAG(l2*omega_r*h)))
!       
!       Rz = Rn/Rd
!          
!       Fn = bessj(l1*omega_r*ri,azim_mode,1)*EXP(ABS(AIMAG(l1*omega_r*ri)))*(Rz*hank1(l2*omega_r*h,azim_mode,1)* & 
!            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*h) + bessj(l2*omega_r*h,azim_mode,1)* &
!            EXP(ABS(AIMAG(l2*omega_r*h))))
!       
!       Fd = bessj(l1*omega_r*h,azim_mode,1)*EXP(ABS(AIMAG(l1*omega_r*h)))*(Rz*dhank1(l2*omega_r,azim_mode,1)* & 
!            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r) + dbessj(l2*omega_r,azim_mode,1)* &
!            EXP(ABS(AIMAG(l2*omega_r))))
!       
!       F = Fn/Fd
!
!       Trs = F/l2
!
!    else if (ss==2) then
!
!       Rd = (lz*bessj(l1*omega_r*h,azim_mode,1)*dhank1(l2*omega_r*h,azim_mode,1)- &
!            hank1(l2*omega_r*h,azim_mode,1)*dbessj(l1*omega_r*h,azim_mode,1))* &
!            EXP(ABS(AIMAG(l1*omega_r*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*h)
!       
!       Rn = (bessj(l2*omega_r*h,azim_mode,1)*dbessj(l1*omega_r*h,azim_mode,1)- &
!            lz*bessj(l1*omega_r*h,azim_mode,1)*dbessj(l2*omega_r*h,azim_mode,1))* &
!            EXP(ABS(AIMAG(l1*omega_r*h))+ABS(AIMAG(l2*omega_r*h)))
!       
!       Rz = Rn/Rd
!          
!       Fn = Rz*hank1(l2*omega_r*ri,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*ri)+ &
!            bessj(l2*omega_r*ri,azim_mode,1)*EXP(ABS(AIMAG(l2*omega_r*ri)))
!       
!       Fd = Rz*dhank1(l2*omega_r,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r)+ &
!            dbessj(l2*omega_r,azim_mode,1)*EXP(ABS(AIMAG(l2*omega_r)))
!       
!       F = Fn/Fd
!
!       Trs = F/l2
!
!    else
!
!!!$       if ((ABS(l3*omega_r) < asymplim .AND. ABS(AIMAG(l3*omega_r)) < asymplim1) .AND. &
!!!$            (ABS(l3*omega_r*ri) < asymplim .AND. ABS(AIMAG(l3*omega_r*ri)) < asymplim1)) then
!
!          Fn = hank1(l3*omega_r*ri,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l3*omega_r*ri)
!       
!          Fd = dhank1(l3*omega_r,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l3*omega_r)
!       
!          F = Fn/Fd
!
!!!$       else 
!!!$       
!!!$          F = (8._dpk*l3*omega_r*ri + 4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*azim_mode*azim_mode - &
!!!$               CMPLX(0._dpk,1._dpk,kind=dpk))/((8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*l3*omega_r - &
!!!$               4._dpk*azim_mode*azim_mode - 3._dpk)*ri**(1.5_dpk))*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
!!!$               l3*omega_r*(ri - 1._dpk))
!!!$
!!!$       end if
!  
!       Trs = F/l3
!
!    end if
!       
!!    print*,'Rz:',Rz
!
!
!  END FUNCTION compute_Trs
!
!
!  FUNCTION compute_Trs_instab(ri,si,ss)  result(Trsin)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Compute Trs of (4.8)
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    real(dpk)     :: ri
!    complex(dpk)  :: sri, Trsin
!    complex(dpk)  :: si
!    complex(dpk)  :: l1, l2
!    complex(dpk)  :: Fn, Fd, F
!    integer       :: ss
!
!
!    l1 = sqrt(1._dpk - si*(M1+1._dpk))*sqrt(1._dpk - si*(M1-1._dpk))
!    l2 = sqrt(1._dpk - si*(M2+1._dpk))*sqrt(1._dpk - si*(M2-1._dpk))
!
!
!    if (ss==1) then
!       
!       Fn = bessj(l1*omega_r*ri,azim_mode,1)*EXP(ABS(AIMAG(l1*omega_r*ri)))
!       Fd = dbessj(l1*omega_r*h,azim_mode,1)*EXP(ABS(AIMAG(l1*omega_r*h)))
!       
!       F = Fn/Fd
!
!       Trsin = F/l1
!
!    else
!
!       Fn = dhank1(l2*omega_r,azim_mode,1)*bessj(l2*omega_r*ri,azim_mode,1)*EXP(ABS(AIMAG(l2*omega_r*ri))+ &
!            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r) - dbessj(l2*omega_r,azim_mode,1)*hank1(l2*omega_r*ri,azim_mode,1)* &
!            EXP(ABS(AIMAG(l2*omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*ri)
!
!       Fd = dhank1(l2*omega_r,azim_mode,1)*dbessj(l2*omega_r*h,azim_mode,1)*EXP(ABS(AIMAG(l2*omega_r*h))+ &
!            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r) - dbessj(l2*omega_r,azim_mode,1)*dhank1(l2*omega_r*h,azim_mode,1)* &
!            EXP(ABS(AIMAG(l2*omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*h)
!       
!       F = Fn/Fd
!
!       Trsin = F/l2
!
!    end if
!       
!!    print*,'Rz:',Rz
!
!
!  END FUNCTION compute_Trs_instab
!
!
!  FUNCTION compute_psi_incident(r,z,ss) result(psi0)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Compute the \Psi_{mn}(r) of (2.14) or (4.11)
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    complex(dpk)        :: psi0, psimn
!    complex(dpk)        :: f1, f2, f3
!    real(dpk)           :: r, z
!    integer             :: ss
!
!    if ((vortswitch == 1) .OR. (vortswitch == 2)) then
!        
!       if (z >= Zo) then 
!          if (ss == 1) then
!                 psimn = omega_r*resp*(1._dpk - M1*mu_plus)*compute_Trs_instab(r,mu_plus,1)
!                 psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*(1._dpk - M1*mu_plus)* &
!                        psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*(z-Zo))
!
!           else if (ss == 2) then
!     
!                psimn = omega_r*resp*(1._dpk - M2*mu_plus)*compute_Trs_instab(r,mu_plus,2)
!                psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*(1._dpk - M2*mu_plus)* &
!                       psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*(z-Zo))
!
!           else
!
!              psi0 = 0.
!
!           end if
!
!       else
!
!          psi0 = 0.
!
!       end if
!
!    else
!
!       if (ss == 1) then
!         psimn = bessj(alpha1*r,azim_mode,1)*EXP(ABS(AIMAG(alpha1*r)))
!       
!         psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*(1._dpk - M1*mu_plus)* &
!              psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*z)
!
!
!       else if (ss == 2) then
!
!          f1 = (1._dpk - mu_plus*M1)/(1._dpk - mu_plus*M2)*bessj(alpha1*h,azim_mode,1)
!       
!          f2 = bessj(alpha2*r,azim_mode,1)* &
!              dhank1(alpha2,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* & 
!              alpha2 + ABS(AIMAG(alpha2*r)))- & 
!              hank1(alpha2*r,azim_mode,1)* &
!              dbessj(alpha2,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* & 
!              alpha2*r + ABS(AIMAG(alpha2)))
!  
!
!          f3 = bessj(alpha2*h,azim_mode,1)* &
!               dhank1(alpha2,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
!               alpha2 + ABS(AIMAG(alpha2*h)))- & 
!               hank1(alpha2*h,azim_mode,1)* &
!               dbessj(alpha2,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
!               alpha2*h + ABS(AIMAG(alpha2)))
!            
!
!          psimn = kap_rho*f1*f2/f3
!
!          psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*(1._dpk - M2*mu_plus)* &
!                psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*z)
!
!!!$       open(10,file='test.out',form='FORMATTED',position='APPEND')
!!!$       write(10,'(2E20.10)') psimn
!!!$       close(10)
!
!       else
!          psi0 = 0.
!
!       end if
!
!    end if
!
!
!  END FUNCTION compute_psi_incident
!
!  FUNCTION residueprpolar(r,phi,ss)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Compute the residue pressure term of (3.38)
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    real(dpk)            :: r, phi
!    integer              :: ss, ii, jj
!    complex(dpk)         :: residueprpolar, res, fn
!   
!    res = psi*(1._dpk - mu_plus*M2)/((mu_plus - sup_zeros_list(ss))* &
!         k_minus_at_mu_plus*kzsp(ss)*(sup_zeros_list(ss) - KH_zero_1)*(sup_zeros_list(ss) - KH_zero_2))
!!!$    res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss) - KH_pole_1)/((mu_plus - sup_zeros_list(ss))* &
!!!$         k_minus_at_mu_plus*kzsp(ss)*(sup_zeros_list(ss) - KH_zero_1)*(sup_zeros_list(ss) - KH_zero_2))
!!!$    res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss)-KH_pole_1)*(sup_zeros_list(ss)-CONJG(KH_pole_1))/ &
!!!$            ((mu_plus - sup_zeros_list(ss))*k_minus_at_mu_plus*kzsp(ss)*(sup_zeros_list(ss)-CONJG(KH_zero_1))* &
!!!$            (sup_zeros_list(ss)-KH_zero_1)*(sup_zeros_list(ss) - KH_zero_2)*(sup_zeros_list(ss)-CONJG(KH_zero_2)))
!    do ii = 1, num_sup_zeros
!       if (ii .NE. (ss)) then
!          res = res/(sup_zeros_list(ss) - sup_zeros_list(ii))
!       end if
!    end do
!    do jj = 1, num_sup_poles
!       res = res*(sup_zeros_list(ss) - sup_poles_list(jj))
!    end do
!    
!    if (r*SIN(phi) <= h) then
!       fn = (1._dpk - sup_zeros_list(ss)*M2)**2*compute_Trs(r*SIN(phi),sup_zeros_list(ss),1)
!    else if (r*SIN(phi) <= 1. .AND. r*SIN(phi) > h) then
!       fn = (1._dpk - sup_zeros_list(ss)*M2)**2*compute_Trs(r*SIN(phi),sup_zeros_list(ss),2)
!    else
!       fn = (1._dpk - sup_zeros_list(ss)*M3)**2*compute_Trs(r*SIN(phi),sup_zeros_list(ss),3)  
!    end if
!
!    residueprpolar = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
!         omega_r*sup_zeros_list(ss)*r*COS(phi))*res*fn
!
!    
!  END FUNCTION residueprpolar
!    
!
!  FUNCTION residuepr(r,z,ss)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Compute the residue pressure term of (3.38)
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    real(dpk)            :: r, z
!    integer              :: ss, ii, jj
!    complex(dpk)         :: residuepr, res, fn
!   
!
!    if (ss == 1) then
!
!!!$       if (z > 0.) then
!!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_1-KH_pole_1)*(KH_zero_1-CONJG(KH_pole_1))/ &
!!!$            ((mu_plus - KH_zero_1)*k_minus_at_mu_plus*kpKH_zero_1*(KH_zero_1-CONJG(KH_zero_1))*(KH_zero_1-KH_zero_2)*(KH_zero_1-CONJG(KH_zero_2)))
!!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_1 - KH_pole_1)/((mu_plus - KH_zero_1)*kmmu_plus*k_plus_sz1*(KH_zero_1 - KH_zero_2))
!       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - KH_zero_1)*k_minus_at_mu_plus*k_plus_sz1*(KH_zero_1 - KH_zero_2))
!       if (vortswitch .EQ. 0) then
!          do ii = 1, num_sup_zeros
!             res = res/(KH_zero_1 - sup_zeros_list(ii))
!          end do
!          do jj = 1, num_sup_poles
!             res = res*(KH_zero_1 - sup_poles_list(jj))
!          end do
!       end if
!
!       if (r <= h) then
!          fn = (1._dpk - KH_zero_1*M2)**2*compute_Trs(r,KH_zero_1,1)
!       else if (r <= 1. .AND. r > h) then
!          fn = (1._dpk - KH_zero_1*M2)**2*compute_Trs(r,KH_zero_1,2)
!       else
!          fn = (1._dpk - KH_zero_1*M3)**2*compute_Trs(r,KH_zero_1,3)
!       end if
!       
!       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
!            omega_r*KH_zero_1*z)*res*fn
!          
!!!$       else
!!!$          residuepr = 0.
!!!$
!!!$       end if
!          
!    elseif (ss == 2) then
!
!!!$       if (z > 0.) then
!!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_2-KH_pole_1)*(KH_zero_2-CONJG(KH_pole_1))/ &
!!!$            ((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2-CONJG(KH_zero_1))*(KH_zero_2-KH_zero_1)*(KH_zero_2-CONJG(KH_zero_2)))
!!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_2 - KH_pole_1)/((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2 - KH_zero_1))
!       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2 - KH_zero_1))
!       if (vortswitch .EQ. 0) then
!          do ii = 1, num_sup_zeros
!             res = res/(KH_zero_2 - sup_zeros_list(ii))
!          end do
!          do jj = 1, num_sup_poles
!             res = res*(KH_zero_2 - sup_poles_list(jj))
!          end do
!       end if
!          
!       if (r <= h) then
!          fn = (1._dpk - KH_zero_2*M2)**2*compute_Trs(r,KH_zero_2,1)
!       else if (r <= 1. .AND. r > h) then
!          fn = (1._dpk - KH_zero_2*M2)**2*compute_Trs(r,KH_zero_2,2)
!       else
!          fn = (1._dpk - KH_zero_2*M3)**2*compute_Trs(r,KH_zero_2,3)
!       end if
!
!       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
!            omega_r*KH_zero_2*z)*res*fn
!
!!!$       else
!!!$          residuepr = 0.
!!!$
!!!$       end if
!
!    else
!
!!!$       if (z > 0.) then
!!!$       res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss-2)-KH_pole_1)*(sup_zeros_list(ss-2)-CONJG(KH_pole_1))/ &
!!!$            ((mu_plus - sup_zeros_list(ss-2))*k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2)-CONJG(KH_zero_1))* &
!!!$            (sup_zeros_list(ss-2)-KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2)*(sup_zeros_list(ss-2)-CONJG(KH_zero_2)))
!!!$       res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss-2) - KH_pole_1)/((mu_plus - sup_zeros_list(ss-2))* &
!!!$            k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2) - KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2))
!       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - sup_zeros_list(ss-2))* &
!            k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2) - KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2))
!       do ii = 1, num_sup_zeros
!          if (ii .NE. (ss-2)) then
!             res = res/(sup_zeros_list(ss-2) - sup_zeros_list(ii))
!          end if
!       end do
!       do jj = 1, num_sup_poles
!          res = res*(sup_zeros_list(ss-2) - sup_poles_list(jj))
!       end do
!          
!       if (r <= h) then
!          fn = (1._dpk - sup_zeros_list(ss-2)*M2)**2*compute_Trs(r,sup_zeros_list(ss-2),1)
!       else if (r <= 1. .AND. r > h) then
!          fn = (1._dpk - sup_zeros_list(ss-2)*M2)**2*compute_Trs(r,sup_zeros_list(ss-2),2)
!       else
!          fn = (1._dpk - sup_zeros_list(ss-2)*M3)**2*compute_Trs(r,sup_zeros_list(ss-2),3)
!       end if
!
!       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
!            omega_r*sup_zeros_list(ss-2)*z)*res*fn
!
!!!$       else
!!!$          residuepr = 0.
!!!$
!!!$       end if
!
!    end if
!
!    
!  END FUNCTION residuepr
!
!
!  FUNCTION residuepot(r,z,ss)
!
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!!! 1. Same as above but for computing velocity potential
!!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!
!    real(dpk)            :: r, z
!    integer              :: ss, ii, jj
!    complex(dpk)         :: residuepot, res, fn
!   
!
!    if (ss == 1) then
!
!!!$       if (z > 0.) then
!!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_1-KH_pole_1)*(KH_zero_1-CONJG(KH_pole_1))/ &
!!!$            ((mu_plus - KH_zero_1)*k_minus_at_mu_plus*k_plus_sz1*(KH_zero_1-CONJG(KH_zero_1))*(KH_zero_1-KH_zero_2)*(KH_zero_1-CONJG(KH_zero_2)))
!!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_1 - KH_pole_1)/((mu_plus - KH_zero_1)*k_minus_at_mu_plus*k_plus_sz1*(KH_zero_1 - KH_zero_2))
!       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - KH_zero_1)*k_minus_at_mu_plus*k_plus_sz1*(KH_zero_1 - KH_zero_2))
!       if (vortswitch .EQ. 0) then
!          do ii = 1, num_sup_zeros
!             res = res/(KH_zero_1 - sup_zeros_list(ii))
!          end do
!          do jj = 1, num_sup_poles
!             res = res*(KH_zero_1 - sup_poles_list(jj))
!          end do
!       end if
!
!       if (r <= h) then
!          fn = (1._dpk - KH_zero_1*M2)**2/(1._dpk - KH_zero_1*M1)*compute_Trs(r,KH_zero_1,1)
!       else if (r <= 1. .AND. r > h) then
!          fn = (1._dpk - KH_zero_1*M2)*compute_Trs(r,KH_zero_1,2)
!       else
!          fn = (1._dpk - KH_zero_1*M3)*compute_Trs(r,KH_zero_1,3)
!       end if
!       
!       residuepot = omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*omega_r*KH_zero_1*z)*res*fn
!          
!!!$       else
!!!$          residuepot = 0.
!!!$
!!!$       end if
!          
!    elseif (ss == 2) then
!
!!!$       if (z > 0.) then
!!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_2-KH_pole_1)*(KH_zero_2-CONJG(KH_pole_1))/ &
!!!$            ((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2-CONJG(KH_zero_1))*(KH_zero_2-KH_zero_1)*(KH_zero_2-CONJG(KH_zero_2)))
!!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_2 - KH_pole_1)/((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2 - KH_zero_1))
!       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2 - KH_zero_1))
!       if (vortswitch .EQ. 0) then
!          do ii = 1, num_sup_zeros
!             res = res/(KH_zero_1 - sup_zeros_list(ii))
!          end do
!          do jj = 1, num_sup_poles
!             res = res*(KH_zero_1 - sup_poles_list(jj))
!          end do
!       end if
!          
!       if (r <= h) then
!          fn = (1._dpk - KH_zero_2*M2)**2/(1._dpk - KH_zero_2*M1)*compute_Trs(r,KH_zero_2,1)
!       else if (r <= 1. .AND. r > h) then
!          fn = (1._dpk - KH_zero_2*M2)*compute_Trs(r,KH_zero_2,2)
!       else
!          fn = (1._dpk - KH_zero_2*M3)*compute_Trs(r,KH_zero_2,3)
!       end if
!
!       residuepot = omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*omega_r*KH_zero_2*z)*res*fn
!
!!!$       else
!!!$          residuepot = 0.
!!!$
!!!$       end if
!
!    else
!
!!!$       if (z > 0.) then
!!!$      res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss-2)-KH_pole_1)*(sup_zeros_list(ss-2)-CONJG(KH_pole_1))/ &
!!!$            ((mu_plus - sup_zeros_list(ss-2))*k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2)-CONJG(KH_zero_1))* &
!!!$            (sup_zeros_list(ss-2)-KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2)*(sup_zeros_list(ss-2)-CONJG(KH_zero_2)))
!!!$       res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss-2) - KH_pole_1)/((mu_plus - sup_zeros_list(ss-2))* &
!!!$            k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2) - KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2))
!       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - sup_zeros_list(ss-2))* &
!            k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2) - KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2))
!       do ii = 1, num_sup_zeros
!          if (ii .NE. (ss-2)) then
!             res = res/(sup_zeros_list(ss-2) - sup_zeros_list(ii))
!          end if
!       end do
!       do jj = 1, num_sup_poles
!          res = res*(sup_zeros_list(ss-2) - sup_poles_list(jj))
!       end do
!          
!       if (r <= h) then
!          fn = (1._dpk - sup_zeros_list(ss-2)*M2)**2/(1._dpk - sup_zeros_list(ss-2)*M1)*compute_Trs(r,sup_zeros_list(ss-2),1)
!       else if (r <= 1. .AND. r > h) then
!          fn = (1._dpk - sup_zeros_list(ss-2)*M2)*compute_Trs(r,sup_zeros_list(ss-2),2)
!       else
!          fn = (1._dpk - sup_zeros_list(ss-2)*M3)*compute_Trs(r,sup_zeros_list(ss-2),3)
!       end if
!
!       residuepot = omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*omega_r*sup_zeros_list(ss-2)*z)*res*fn
!
!!!$       else
!!!$          residuepot = 0.
!!!$
!!!$       end if
!
!    end if
!
!    
!  END FUNCTION residuepot
!
!
end module user_defined_functions


