PROGRAM check_residue

!!$ List of functions: (func_case = X)
!!$-----------------------------------------------------------
!!$ X = 0    :: Samanta & Freund, JFM, 2008 (incident wave)
!!$ X = 100  :: Samanta & Freund, JFM, 2008 (reflected wave)
!!$ X = 200  :: Munt's jet (incident wave)
!!$-----------------------------------------------------------
!!$ X = 1    :: Samanta & Freund, JFM, 2008 (zero mode) (3.22)
!!$ X = -1   :: Samanta & Freund, JFM, 2008 (pole mode) (3.22)
!!$ X = 2    :: Gabard & Astley, JFM, 2008 (zero mode)
!!$ X = -2   :: Gabard & Astley, JFM, 2008 (pole mode)
!!$ X = 3    :: Taylor et al, JSV, 1993 (zero mode)
!!$ X = -3   :: Taylor et al, JSV, 1993 (pole mode)
!!$ X = 4    :: Samanta & Freund, JFM, 2008 (zero mode) (4.10)
!!$ X = -4   :: Samanta & Freund, JFM, 2008 (pole mode) (4.10)
!!$ X = 11   :: Munt's jet (zeros)
!!$ X = -11  :: Munt's jet (poles)
!!$ X = 101  :: Test kernel 1 (zeros)
!!$ X = -101 :: Test kernel 1 (poles)

  USE bessel_utils
  USE io_utils

  IMPLICIT none
  
  integer, parameter                        :: dpk = kind(1.0d0) !! Double precision kind
  integer                                   :: Nx, Ny !! Mesh points
  real(dpk)                                 :: Xmin, Xmax, Ymin, Ymax !! Search box
  real(dpk)                                 :: Xm1, Xm2, Xmid
  real(dpk), allocatable, dimension(:)      :: X, Y !! ''
  integer                                   :: num_of_points
  complex(dpk), allocatable, dimension(:)   :: input_list !! ''
  integer                                   :: prec !! nag_bessel precision digits
  real(dpk)                                 :: M1, M2, M3, h
  real(dpk)                                 :: St_flag
  complex(dpk)                              :: omega
  real(dpk)                                 :: azim_mode
  real(dpk)                                 :: del  !! phase angle in degrees
  real(dpk)                                 :: PI
  real(dpk), parameter                      :: asymplim = 3.27E4
  real(dpk), parameter                      :: asymplim1 = 100
  character(len=200)                        :: dummy

  call readdata                          !! Reads the input file


CONTAINS

  SUBROUTINE readdata

!! the basic data file:
    real(dpk)           :: temp1, temp2

    PI = 4._dpk*ATAN(1._dpk)


    open(unit=10, file='input_residue.list', status='old', action='read')

      read(10,*) func_case, dummy
      print '(A, I5, A)', ' func_case = (', func_case, ')'


      read(10,*) M1,dummy
      read(10,*) M2,dummy
      read(10,*) M3,dummy
      print '(A, F8.4, ", ", F8.4, ", ", F8.4, A)', &
                   ' (M1, M2, M3) = (', M1, M2, M3, ')'

      read(10,*) h,dummy
      read(10,*) del,dummy
      read(10,*) azim_mode,dummy
      print '(A, F8.4, ", ", F8.4, ", ", F8.4,A)', &
                   ' (h,del, azim_mode) = (', h, del, azim_mode, ')'

      read(10,*) St_flag,dummy
      read(10,*) omega,dummy
      read(10,*) num_of_inputs,dummy

       if (num_of_inputs > 0) then
         allocate(input_list(num_of_inputs))
         do i = 1, num_of_inputs
           read(10,*) input_list(i)
           print*, input(i)
         end do
       end if

     close(10)

    END SUBROUTINE readdata
 
  FUNCTION fun(z,w)

    !***! This subroutine specifies the function.

    complex(dpk)              :: z, fun, w
    complex(dpk)              :: lambda1, lambda2, l3, lz, F1n, F1d, F2n, F2d, F1s, F2s, F1, F2
    complex(dpk)              :: Rnz, Rdz, Rz, R1nz, R1dz, R1z, R2nz, R2dz, R2z


    lambda1 = sqrt(1._dpk - z*(M1+1._dpk))*sqrt(1._dpk - z*(M1-1._dpk))
    lambda2 = sqrt(1._dpk*kap_T - z*(kap_T*M2+1._dpk))*sqrt(1._dpk*kap_T - z*(kap_T*M2-1._dpk))

    SELECT CASE (func_case)

    CASE (0)

!!$!! Samanta & Freund, JFM, 2008 (Acoustic Incident Waves)

       lambda1 = w*sqrt((1._dpk - z*M1)**2 - z**2)
       lambda2 = w*sqrt(kap_T**2*(1._dpk - z*M2)**2 - z**2)

       F1n = dbessj(h*lambda2,azim_mode,1)*dhank1(lambda2,azim_mode,1)*EXP(ABS(AIMAG(h*lambda2)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2) - dhank1(h*lambda2,azim_mode,1)*dbessj(lambda2,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2)) + CMPLX(0._dpk,1._dpk,kind=dpk)*h*lambda2)

       F1 = lambda2*(1._dpk - z*M1)**2*bessj(h*lambda1,azim_mode,1)*EXP(ABS(AIMAG(h*lambda1)))*F1n

       F2n = bessj(h*lambda2,azim_mode,1)*dhank1(lambda2,azim_mode,1)*EXP(ABS(AIMAG(h*lambda2)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2) - hank1(h*lambda2,azim_mode,1)*dbessj(lambda2,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2)) + CMPLX(0._dpk,1._dpk,kind=dpk)*h*lambda2)

       F2 = lambda1*(1._dpk - z*M2)**2*dbessj(h*lambda1,azim_mode,1)*EXP(ABS(AIMAG(h*lambda1)))*F2n

       fun = F1 - F2

    CASE (100)

!!$!! Samanta & Freund, JFM, 2008 (Acoustic Reflected Waves)

       lambda1 = w*sqrt((1._dpk - z*M1)**2 - z**2)
       lambda2 = w*sqrt(kap_T**2*(1._dpk - z*M2)**2 - z**2)

       F1n = dbessj(h*lambda2,azim_mode,1)*dhank1(lambda2,azim_mode,1)*EXP(ABS(AIMAG(h*lambda2)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2) - dhank1(h*lambda2,azim_mode,1)*dbessj(lambda2,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2)) + CMPLX(0._dpk,1._dpk,kind=dpk)*h*lambda2)

       F1 = lambda2*(1._dpk - z*M1)**2*bessj(h*lambda1,azim_mode,1)*EXP(ABS(AIMAG(h*lambda1)))*F1n

       F2n = bessj(h*lambda2,azim_mode,1)*dhank1(lambda2,azim_mode,1)*EXP(ABS(AIMAG(h*lambda2)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2) - hank1(h*lambda2,azim_mode,1)*dbessj(lambda2,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2)) + CMPLX(0._dpk,1._dpk,kind=dpk)*h*lambda2)

       F2 = lambda1*(1._dpk - z*M2)**2*dbessj(h*lambda1,azim_mode,1)*EXP(ABS(AIMAG(h*lambda1)))*F2n

       fun = F1 - F2

    CASE (200)

!!$!! Hollow duct incident waves

       lambda1 = w*sqrt((1._dpk - z*M1)**2 - z**2)

       F1n = dbessj(lambda1,azim_mode,1)*EXP(ABS(AIMAG(lambda1)))

       fun = F1n

    CASE (1)

!!$!! Samanta & Freund, JFM, 2008 (3.22) (computes the zeros)

       lz = kap_rho*lambda2/lambda1*(1._dpk-z*M1)*(1._dpk-z*M1)/((1._dpk-z*M2)*(1._dpk-z*M2))
       l3 = sqrt(1._dpk*kap_T - z*(kap_T*M3+1._dpk))*sqrt(1._dpk*kap_T - z*(kap_T*M3-1._dpk))

       if ((ABS(lambda2*w) < asymplim .AND. ABS(AIMAG(lambda2*w)) < asymplim1) .AND. &
            (ABS(lambda2*w*h) < asymplim .AND. ABS(AIMAG(lambda2*w*h)) < asymplim1) .AND. &
            (ABS(lambda1*w*h) < asymplim .AND. ABS(AIMAG(lambda1*w*h)) < asymplim1)) then
      
          Rdz = (lz*bessj(lambda1*w*h,azim_mode,1)*dhank1(lambda2*w*h,azim_mode,1) -&
                       hank1(lambda2*w*h,azim_mode,1)*dbessj(lambda1*w*h,azim_mode,1))* &
                       EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h)
          Rnz = (bessj(lambda2*w*h,azim_mode,1)*dbessj(lambda1*w*h,azim_mode,1) - &
                      lz*bessj(lambda1*w*h,azim_mode,1)*dbessj(lambda2*w*h,azim_mode,1))* &
                         EXP(ABS(AIMAG(lambda2*w*h)))

!!$              Rnz = lz*bessj(lambda1*w*h,azim_mode,0)*dhank1(lambda2*w*h,azim_mode,0) - hank1(lambda2*w*h,azim_mode,0)*dbessj(lambda1*w*h,azim_mode,0)
!!$              Rdz = bessj(lambda2*w*h,azim_mode,0)*dbessj(lambda1*w*h,azim_mode,0) - lz*bessj(lambda1*w*h,azim_mode,0)*dbessj(lambda2*w*h,azim_mode,0)

          Rz = Rnz/Rdz
          
          F1n = Rz*hank1(lambda2*w,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) + bessj(lambda2*w,azim_mode,1)* &
               EXP(ABS(AIMAG(lambda2*w)))
          F1d = Rz*dhank1(lambda2*w,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) + dbessj(lambda2*w,azim_mode,1)* &
               EXP(ABS(AIMAG(lambda2*w)))

!!$              F1n = hank1(lambda2*w,azim_mode,0) + Rz*bessj(lambda2*w,azim_mode,0)
!!$              F1d = dhank1(lambda2*w,azim_mode,0) + Rz*dbessj(lambda2*w,azim_mode,0)
          
          F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/lambda2*F1n/F1d

       else

          F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/lambda2*CMPLX(0._dpk,1._dpk,kind=dpk)

       end if

       if (ABS(l3*w) < asymplim .AND. ABS(AIMAG(l3*w)) < asymplim1) then
       
          F2n = hank1(l3*w,azim_mode,1)
          F2d = dhank1(l3*w,azim_mode,1)

!!$              F2n = hank1(l3*w,azim_mode,0)
!!$              F2d = dhank1(l3*w,azim_mode,0)

          F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*F2n/F2d

       else

          F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*(8._dpk*l3*w + &
            4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*azim_mode*azim_mode - &
            CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*l3*w - &
            4._dpk*azim_mode*azim_mode - 3._dpk)

       end if
      
       fun = w*(F1 - F2)
       
    CASE (-1)
       
!!$!! Samanta & Freund, JFM, 2008 (3.22) (computes the poles)

       lz = kap_rho*lambda2/lambda1*(1._dpk-z*M1)*(1._dpk-z*M1)/((1._dpk-z*M2)*(1._dpk-z*M2))
       l3 = sqrt(1._dpk*kap_T - z*(kap_T*M3+1._dpk))*sqrt(1._dpk*kap_T - z*(kap_T*M3-1._dpk))

       if ((ABS(lambda2*w) < asymplim .AND. ABS(AIMAG(lambda2*w)) < asymplim1) .AND. &
            (ABS(lambda2*w*h) < asymplim .AND. ABS(AIMAG(lambda2*w*h)) < asymplim1) .AND. &
            (ABS(lambda1*w*h) < asymplim .AND. ABS(AIMAG(lambda1*w*h)) < asymplim1)) then
       
          Rdz = (lz*bessj(lambda1*w*h,azim_mode,1)*dhank1(lambda2*w*h,azim_mode,1) - &
                             hank1(lambda2*w*h,azim_mode,1)*dbessj(lambda1*w*h,azim_mode,1))* &
                             EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h)
          Rnz = (bessj(lambda2*w*h,azim_mode,1)*dbessj(lambda1*w*h,azim_mode,1) -&
                             lz*bessj(lambda1*w*h,azim_mode,1)*dbessj(lambda2*w*h,azim_mode,1))* &
                              EXP(ABS(AIMAG(lambda2*w*h)))
       
          !    Rnz = lz*bessj(lambda1*w*h,azim_mode,0)*dhank1(lambda2*w*h,azim_mode,0) - hank1(lambda2*w*h,azim_mode,0)*dbessj(lambda1*w*h,azim_mode,0)
          !    Rdz = bessj(lambda2*w*h,azim_mode,0)*dbessj(lambda1*w*h,azim_mode,0) - lz*bessj(lambda1*w*h,azim_mode,0)*dbessj(lambda2*w*h,azim_mode,0)
       
          Rz = Rnz/Rdz
       
          F1n = Rz*hank1(lambda2*w,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) + &
                                                                    bessj(lambda2*w,azim_mode,1)* &
                                                                    EXP(ABS(AIMAG(lambda2*w)))
          F1d = Rz*dhank1(lambda2*w,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) + &
                                                                 dbessj(lambda2*w,azim_mode,1)* &
                                                                     EXP(ABS(AIMAG(lambda2*w)))
       
          !    F1n = hank1(lambda2*w,azim_mode,0) + Rz*bessj(lambda2*w,azim_mode,0)
          !    F1d = dhank1(lambda2*w,azim_mode,0) + Rz*dbessj(lambda2*w,azim_mode,0)
       
          F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/lambda2*F1n/F1d

       else

          F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/lambda2*CMPLX(0._dpk,1._dpk,kind=dpk)

       end if

       if (ABS(l3*w) < asymplim .AND. ABS(AIMAG(l3*w)) < asymplim1) then
       
          F2n = hank1(l3*w,azim_mode,1)
          F2d = dhank1(l3*w,azim_mode,1)
       
          !    F2n = hank1(l3*w,azim_mode,0)
          !    F2d = dhank1(l3*w,azim_mode,0)
       
          F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*F2n/F2d

       else

          F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*(8._dpk*l3*w + &
               4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*azim_mode*azim_mode - &
               CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*l3*w - &
               4._dpk*azim_mode*azim_mode - 3._dpk)

       end if
       
       fun = w*(F1 - F2)
       fun = 1._dpk/fun
       
    CASE (2)

!!$ Gabard & Astley, JFM, 2006 (compute the zeros)

       F1n = dhank1(lambda1*w*h,azim_mode,1)*bessj(lambda1*w,azim_mode,1)*EXP(ABS(AIMAG(lambda1*w))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda1*w*h) - dbessj(lambda1*w*h,azim_mode,1)*hank1(lambda1*w,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda1*w*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda1*w)
       F1d = dhank1(lambda1*w*h,azim_mode,1)*dbessj(lambda1*w,azim_mode,1)*EXP(ABS(AIMAG(lambda1*w))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda1*w*h) - dbessj(lambda1*w*h,azim_mode,1)*dhank1(lambda1*w,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda1*w*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda1*w)
       F1s = F1n/F1d

       F1 = kap_rho*(1._dpk - z*M1)*(1._dpk - z*M1)*lambda2*F1s
       
       F2n = hank1(lambda2*w,azim_mode,1)
       F2d = dhank1(lambda2*w,azim_mode,1)
       F2s = F2n/F2d

       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)*lambda1*F2s
       
       fun = F1 - F2

    CASE (-2)

!!$ Gabard & Astley, JFM, 2006 (compute the poles)

       F1n = dhank1(lambda1*w*h,azim_mode,1)*bessj(lambda1*w,azim_mode,1)*EXP(ABS(AIMAG(lambda1*w))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda1*w*h) - dbessj(lambda1*w*h,azim_mode,1)*hank1(lambda1*w,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda1*w*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda1*w)
       F1d = dhank1(lambda1*w*h,azim_mode,1)*dbessj(lambda1*w,azim_mode,1)*EXP(ABS(AIMAG(lambda1*w))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda1*w*h) - dbessj(lambda1*w*h,azim_mode,1)*dhank1(lambda1*w,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda1*w*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda1*w)
       F1s = F1n/F1d

       F1 = kap_rho*(1._dpk - z*M1)*(1._dpk - z*M1)*lambda2*F1s
       
       F2n = hank1(lambda2*w,azim_mode,1)
       F2d = dhank1(lambda2*w,azim_mode,1)
       F2s = F2n/F2d

       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)*lambda1*F2s
       
       fun = F1 - F2
       fun = 1._dpk/fun

    CASE (3)

!! Taylor et al, JSV, 1993 (zeros)

       lz = kap_rho*lambda1/lambda2*(1._dpk-z*M2)*(1._dpk-z*M2)/((1._dpk-z*M1)*(1._dpk-z*M1))
       l3 = sqrt(1._dpk*kap_T - z*(kap_T*M3+1._dpk))*sqrt(1._dpk*kap_T - z*(kap_T*M3-1._dpk))
       
       Rnz = (lz*bessj(lambda2*w*h,azim_mode,1)*dbessj(lambda1*w*h,azim_mode,1) - &
                                      bessj(lambda1*w*h,azim_mode,1)*dbessj(lambda2*w*h,azim_mode,1))* &
                                      EXP(ABS(AIMAG(lambda2*w*h)))
       Rdz = (bessj(lambda1*w*h,azim_mode,1)*dhank2(lambda2*w*h,azim_mode,1) - &
                lz*dbessj(lambda1*w*h,azim_mode,1)*hank2(lambda2*w*h,azim_mode,1))* &
                 EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h)

       Rz = Rnz/Rdz
       
       F1n = bessj(lambda2*w,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w))) + Rz*hank2(lambda2*w,azim_mode,1)* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w)
       F1d = dbessj(lambda2*w,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w))) + Rz*dhank2(lambda2*w,azim_mode,1)* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w)

       F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/lambda2*F1n/F1d
       
       F2n = hank2(l3*w,azim_mode,1)
       F2d = dhank2(l3*w,azim_mode,1)
       F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*F2n/F2d
       
       fun = w*(F1 - F2)

    CASE (-3)

!! Taylor et al, JSV, 1993 (poles)

       lz = kap_rho*lambda1/lambda2*(1._dpk-z*M2)*(1._dpk-z*M2)/((1._dpk-z*M1)*(1._dpk-z*M1))
       l3 = sqrt(1._dpk*kap_T - z*(kap_T*M3+1._dpk))*sqrt(1._dpk*kap_T - z*(kap_T*M3-1._dpk))
       
       Rnz = (lz*bessj(lambda2*w*h,azim_mode,1)*dbessj(lambda1*w*h,azim_mode,1) - &
                                 bessj(lambda1*w*h,azim_mode,1)*dbessj(lambda2*w*h,azim_mode,1))* &
                                 EXP(ABS(AIMAG(lambda2*w*h)))
       Rdz = (bessj(lambda1*w*h,azim_mode,1)*dhank2(lambda2*w*h,azim_mode,1) - &
                              lz*dbessj(lambda1*w*h,azim_mode,1)*hank2(lambda2*w*h,azim_mode,1))* &
                               EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h)

       Rz = Rnz/Rdz
       
       F1n = bessj(lambda2*w,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w))) + Rz*hank2(lambda2*w,azim_mode,1)* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w)
       F1d = dbessj(lambda2*w,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w))) + Rz*dhank2(lambda2*w,azim_mode,1)* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w)

       F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/lambda2*F1n/F1d
       
       F2n = hank2(l3*w,azim_mode,1)
       F2d = dhank2(l3*w,azim_mode,1)
       F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*F2n/F2d
       
       fun = w*(F1 - F2)
       fun = 1._dpk/fun

    CASE (4)

!!$ Samanta & Freund, JFM, 2008 (4.10) (zeros)

       F1n = bessj(lambda1*w*h,azim_mode,1)
       F1d = dbessj(lambda1*w*h,azim_mode,1)
       F1s = F1n/F1d

       F1 = kap_rho*(1._dpk - z*M1)*(1._dpk - z*M1)/lambda1*F1s

       F2n = dhank1(lambda2*w,azim_mode,1)*bessj(lambda2*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) - dbessj(lambda2*w,azim_mode,1)*hank1(lambda2*w*h,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2*w))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h)
       F2d = dhank1(lambda2*w,azim_mode,1)*dbessj(lambda2*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) - dbessj(lambda2*w,azim_mode,1)*dhank1(lambda2*w*h,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2*w))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h)
       F2s = F2n/F2d
       
       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)/lambda2*F2s
       
       fun = w*(F1 - F2)

    CASE (-4)

!!$ Samanta & Freund, JFM, 2008 (4.10) (poles)

       F1n = bessj(lambda1*w*h,azim_mode,1)
       F1d = dbessj(lambda1*w*h,azim_mode,1)
       F1s = F1n/F1d

       F1 = kap_rho*(1._dpk - z*M1)*(1._dpk - z*M1)/lambda1*F1s

       F2n = dhank1(lambda2*w,azim_mode,1)*bessj(lambda2*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) - dbessj(lambda2*w,azim_mode,1)*hank1(lambda2*w*h,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2*w))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h)
       F2d = dhank1(lambda2*w,azim_mode,1)*dbessj(lambda2*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) - dbessj(lambda2*w,azim_mode,1)*dhank1(lambda2*w*h,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2*w))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h)
       F2s = F2n/F2d
       
       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)/lambda2*F2s
       
       fun = (F1 - F2)
       fun = 1._dpk/fun

    CASE (11)

!!$ Munt's Duct (zeros)

       F1n = bessj(lambda1*w,azim_mode,1)
       F1d = dbessj(lambda1*w,azim_mode,1)
       F1s = F1n/F1d
       
       F1 = kap_rho*(1._dpk - z*M1)*(1._dpk - z*M1)*F1s/lambda1
       
       F2n = hank1(lambda2*w,azim_mode,1)
       F2d = dhank1(lambda2*w,azim_mode,1)
       F2s = F2n/F2d
       
       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)*F2s/lambda2
       
       fun = (F1 - F2)

    CASE (-11)

!!$ Munt's Duct (poles)

       F1n = bessj(lambda1*w,azim_mode,1)
       F1d = dbessj(lambda1*w,azim_mode,1)
       F1s = F1n/F1d
       
       F1 = kap_rho*(1._dpk - z*M1)*(1._dpk - z*M1)*lambda2*F1s
       
       F2n = hank1(lambda2*w,azim_mode,1)
       F2d = dhank1(lambda2*w,azim_mode,1)
       F2s = F2n/F2d
       
       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)*lambda1*F2s
       
       fun = F1 - F2
       fun = 1._dpk/fun


    CASE (101)

!! Test functions: Supersonic kernel (2nd part)  ! No density ratio

       F1s = dbessj(lambda1*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda1*w*h)))
       F2s = bessj(lambda1*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda1*w*h)))

       F1n = (dhank1(lambda2*w,azim_mode,1)*bessj(lambda2*w,azim_mode,1) - &
            dbessj(lambda2*w,azim_mode,1)*hank1(lambda2*w,azim_mode,1))*EXP(ABS(AIMAG(lambda2*w)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w)  ! this term has no zeros; do not include

       F1 = F1s

       F1d = lz*F2s*(dhank1(lambda2*w,azim_mode,1)*dbessj(lambda2*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) - dbessj(lambda2*w,azim_mode,1)*dhank1(lambda2*w*h,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2*w)) + CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h))
       
       F2d = F1s*(dhank1(lambda2*w,azim_mode,1)*bessj(lambda2*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) - dbessj(lambda2*w,azim_mode,1)*hank1(lambda2*w*h,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2*w)) + CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h))

       F2 = F1d - F2d
       
       fun = F1/F2
       
    CASE (-101)

!! Test functions: Supersonic kernel (2nd part)  ! No density ratio

       F1s = dbessj(lambda1*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda1*w*h)))
       F2s = bessj(lambda1*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda1*w*h)))

       F1n = (dhank1(lambda2*w,azim_mode,1)*bessj(lambda2*w,azim_mode,1) - &
            dbessj(lambda2*w,azim_mode,1)*hank1(lambda2*w,azim_mode,1))*EXP(ABS(AIMAG(lambda2*w)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w)  ! this term has no zeros; do not include

       F1 = F1s

       F1d = lz*F2s*(dhank1(lambda2*w,azim_mode,1)*dbessj(lambda2*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) - dbessj(lambda2*w,azim_mode,1)*dhank1(lambda2*w*h,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2*w)) + CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h))

       F2d = F1s*(dhank1(lambda2*w,azim_mode,1)*bessj(lambda2*w*h,azim_mode,1)*EXP(ABS(AIMAG(lambda2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w) - dbessj(lambda2*w,azim_mode,1)*hank1(lambda2*w*h,azim_mode,1)* &
            EXP(ABS(AIMAG(lambda2*w)) + CMPLX(0._dpk,1._dpk,kind=dpk)*lambda2*w*h))
       
       F2 = F1d - F2d
       
       fun = F1/F2
       fun = 1._dpk/fun

    CASE DEFAULT 

       STOP 'Not a valid function to find roots!! Program Terminated.'
       
    END SELECT


  END FUNCTION fun
!
!
  FUNCTION dfun(z,w,h,err,Nt)

!***! Computes the numerical derivative using the Ridders' method (Adv. Eng. Software.,
!***! Vol. 4, No. 2, 1982; also see Numerical Recipes, Section 5.7).

    complex(dpk)                          :: z, dfun, w
    real(dpk)                             :: err, fac, Nt, h, hh
    integer, parameter                    :: NTAB = 20 !! Tableau size
    complex(dpk), dimension(NTAB,NTAB)    :: a
    real(dpk), dimension(NTAB-1)          :: errt
    real(dpk), parameter                  :: DIV = 1.5_dpk, BIG = HUGE(ABS(z)), SAFE = 2._dpk
    integer, dimension(1)                 :: imin
    integer                               :: ierrmin, i, j


    hh = h
    a(1,1) = (fun(z+hh,w) - fun(z-hh,w))/(2._dpk*hh)
    err = BIG
    
    do i = 2, NTAB

       hh = hh/DIV !! Higher extrapolation orders are also of reduced step size
       a(1,i) = (fun(z+hh,w) - fun(z-hh,w))/(2._dpk*hh)
       fac = DIV*DIV

       do j = 2, i
          a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac - 1._dpk)
          fac = DIV*DIV*fac
       end do

       errt(1:i-1) = MAX(ABS(a(2:i,i) - a(1:i-1,i)), ABS(a(2:i,i) - a(1:i-1,i-1)))
       imin = MINLOC(errt(1:i-1))
       ierrmin = imin(1)

       if (errt(ierrmin) <= err) then !! Error checking
          err = errt(ierrmin)
          dfun = a(ierrmin+1,i)
       end if

       if (ABS(a(i,i) - a(i-1,i-1)) >= SAFE*err) then !! When the higher order becomes worse,
          Nt = i/NTAB                                 !! quit early
          RETURN
       end if

    end do

    Nt = 1._dpk


  END FUNCTION dfun


  SUBROUTINE mesh

!***! The computation grid.

    real(dpk)    :: dx, dy, dx1
    integer      :: i, j

    if (func_case .EQ. 200) then
       if ( M1 < 1) then
          dx = (Xmax - Xmin)/(Nx - 1)
          dy = (Ymax - Ymin)/(Ny - 1)
    
          allocate(X(Nx))
          allocate(Y(Ny))
    
          X(1) = Xmin
          X(Nx) = Xmax
          Y(1) = Ymin
          Y(Ny) = Ymax
    
          do i = 1, Nx-2
             X(i+1) = X(i) + dx
          end do
          do i = 1, Ny-2
             Y(i+1) = Y(i) + dy
          end do
       else
          dx = (Xmax - Xm2)/(Nx/2 - 1)
          dx1 = (Xm1 - Xmin)/(Nx/2 - 1)
          dy = (Ymax - Ymin)/(Ny - 1)
    
          allocate(X(Nx))
          allocate(Y(Ny))
    
          X(1) = Xmin
          X(Nx/2) = Xm1
          X(Nx/2+1) = Xm2
          X(Nx) = Xmax
          Y(1) = Ymin
          Y(Ny) = Ymax
    
          do i = 1, Nx/2-2
             X(i+1) = X(i) + dx1
          end do
          do i = Nx/2+1, Nx-2
             X(i+1) = X(i) + dx
          end do
          do i = 1, Ny-2
             Y(i+1) = Y(i) + dy
          end do
       end if

    else

       dx = (Xmax - Xmin)/(Nx - 1)
       dy = (Ymax - Ymin)/(Ny - 1)
    
       allocate(X(Nx))
       allocate(Y(Ny))
    
       X(1) = Xmin
       X(Nx) = Xmax
       Y(1) = Ymin
       Y(Ny) = Ymax
    
       do i = 1, Nx-2
          X(i+1) = X(i) + dx
       end do
       do i = 1, Ny-2
          Y(i+1) = Y(i) + dy
       end do

    end if

    open(20,file='./DataDump/mesh.x')
    do i = 1, Nx
       write(20,'(I10,2X,F20.10)') i, X(i)
    end do
    close(20)
    
    open(20,file='./DataDump/mesh.y')
    do i = 1, Ny
       write(20,'(I10,2X,F20.10)') i, Y(i)
    end do
    
  END SUBROUTINE mesh

!
  FUNCTION newt(gz,w,flg,grid_count,total_grid_points,file_ID)

!***! The Newton-Raphson routine for complex functions using derivatives. The
!***! algorithm removes zeros from the function as soon as it is found, thereby
!***! preventing duplication and speeding up the process. This also enables
!***! to find multiple zeros, if any.

    complex(dpk)       :: newt, gz, w
    complex(dpk)       :: f, df, dz, zp, zq, zw
    real(dpk)          :: newtr, newti
    real(dpk)          :: err1, err2
    character(10)      :: flg
    integer            :: i, j, k, grid_count, total_grid_points, file_ID


    newt = gz
    flg = 'red' !! 'Green' flag means zero found
    do i = 1, MAX_ITE
       f = fun(newt,w)
       df = dfun(newt,w,delta,err1,err2) !! Numerical derivative

       dz = f/df
       newt = newt - dz

       newtr = REAL(newt)
       newti = AIMAG(newt)

       if(newtr > Xmax+ZERO_ACC/10 .OR. newtr < Xmin-ZERO_ACC/10 .OR. &
            newti > Ymax+ZERO_ACC/10 .OR. newti < Ymin-ZERO_ACC/10) then
           
             write(file_ID,*) grid_count, "/", total_grid_points, " Jumped out of bounds"
           RETURN
       end if

       if (ABS(REAL(dz)) < ZERO_ACC .AND. ABS(AIMAG(dz)) < ZERO_ACC) then 
          flg = 'green'
          RETURN
       end if

    end do
             
    write(file_ID,*) grid_count, "/", total_grid_points, " WARNING: Max Iterations Exceeded!"


  END FUNCTION newt

!

END PROGRAM check_residue
