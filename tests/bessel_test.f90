PROGRAM nrcomplex

  USE bessel_utils

  IMPLICIT none
 
  integer, parameter                        :: dpk = kind(1.0d0) !! Double precision kind
  integer                                   :: Nx, Ny !! Mesh points
  real(dpk)                                 :: Xmin, Xmax, Ymin, Ymax !! Search box
  real(dpk)                                 :: Xm1, Xm2, Xmid
  real(dpk), allocatable, dimension(:)      :: X, Y !! ''
  integer                                   :: Nzero, Nzero_inc !! Computed zeros
  complex(dpk), allocatable, dimension(:,:) :: zero !! Records all the activity over (X,Y)
  complex(dpk), allocatable, dimension(:)   :: zerolist !! Lists only the zeros
  complex(dpk), allocatable, dimension(:)   :: checkzero !! Checks if it's really one
  real(dpk), allocatable, dimension(:,:)    :: alpha, alphap !! Real wave nos for the inc case
  real(dpk)                                 :: delta !! Derivative step size
  integer                                   :: MAX_ITE !! Max iterations in N-R scheme
  integer                                   :: MAX_ZERO !! Max estimated zeros
  real(dpk)                                 :: ZERO_ACC !! The accuracy in finding zeros
  real(dpk)                                 :: ZERO_PREC !! The precision in finding zeros
  integer                                   :: prec !! nag_bessel precision digits
  real(dpk)                                 :: M1, M2, M3
  real(dpk)                                 :: h, w0
  real(dpk)                                 :: m
  real(dpk)                                 :: del  !! phase angle in degrees
  real(dpk)                                 :: kap1  !! sqrt(T1/T0)
  real(dpk)                                 :: kap2  !! rho1/rho0
  complex(dpk)                              :: omega, z_input, residue
  real(dpk)                                 :: PI
  integer                                   :: step, limit, func_case
  real(dpk), parameter                      :: asymplim = 3.27E4
  real(dpk), parameter                      :: asymplim1 = 100
  character(len=30)                         :: dummy

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

  call readdata                          !! Reads the input file

  residue = fun(z_input,omega)

  print '(A, 2F8.3, A)', ' residue = (', REAL(residue), AIMAG(residue), ')'

CONTAINS


  SUBROUTINE readdata

!! the basic data file:
    real(dpk)           :: temp1, temp2

    PI = 4._dpk*ATAN(1._dpk)


    open(unit=10, file='input.list', status='old', action='read')

      read(10,*) z_input, dummy
      read(10,*) func_case, dummy
      print '(A, I5, A)', ' func_case = (', func_case, ')'


      read(10,*) M1,        dummy
      read(10,*) M2,        dummy
      read(10,*) M3,        dummy
      PRINT '(A, F8.4, ", ", F8.4, ", ", F8.4, A)', &
                   ' (M1, M2, M3) = (', M1, M2, M3, ')'

      read(10,*) h,         dummy
      read(10,*) w0,        dummy
      read(10,*) del,       dummy
      read(10,*) m,         dummy
      PRINT '(A, F8.4, ", ", F8.4, ", ", F8.4, ", ", F8.4,A)', &
                   ' (h, w0, del, m) = (', h, w0, del, m, ')'


      read(10,*) Nx,        dummy
      read(10,*) Ny,        dummy

      PRINT '(A, I5, ", ", I5,A)', &
                   ' (Nx, Ny) = (', Nx, Ny, ')'


      read(10,*) Xmin, dummy           
      read(10,*) Xmax, dummy           

      PRINT '(A, F8.4, ", ", F8.4,A)', &
                   ' (Xmin, Xmax) = (', Xmin, Xmax, ')'


      read(10,*) Ymin, dummy           
      read(10,*) Ymax, dummy           

      PRINT '(A, F8.4, ", ", F8.4,A)', &
                   ' (Ymin, Ymax) = (', Ymin, Ymax, ')'


      read(10,*) MAX_ITE,   dummy
      PRINT '(A, I5, A)', &
                   ' Max iterations = (', MAX_ITE, ')'


      read(10,*) ZERO_ACC,  dummy
      read(10,*) ZERO_PREC, dummy
      PRINT '(A, E12.4, ", ", E12.4,A)', &
                   ' (ZERO_ACC, ZERO_PREC) = (', ZERO_ACC, ZERO_PREC, ')'

      read(10,*) MAX_ZERO,  dummy
      print *, "MAX_ZERO =", MAX_ZERO

      read(10,*) prec,      dummy
      read(10,*) step,      dummy
      read(10,*) delta,     dummy
      read(10,*) limit,     dummy
      PRINT '(A, I5, ", ", I5, ", ", E12.4, ", ", I5,A)', &
                   ' (prec, step, delta, limit) = (', prec, step, delta, limit, ')'


      
     close(10)

     omega = ABS(w0)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*del*PI/180)

     IF(func_case .EQ. 0 .OR. func_case .EQ. 100) THEN
        Ny = 2
        Xmin = -kap1/(1._dpk - kap1*M3) + 1.E-6
        Xmax = 1._dpk/(1._dpk + M1) - 1.E-6
        Ymin = 0.
        Ymax = 0.
     END IF

     IF(func_case .EQ. 200) THEN
        Ny = 2
        Xmin = -1._dpk/(1._dpk - M1) + 1.E-6
        Xmax = 1._dpk/(1._dpk + M1) - 1.E-6

        IF(M1 > 1) THEN
           temp1 = Xmin
           temp2 = Xmax
           Xmid = (Xmin + Xmax)/2
           Xm1 = temp2 - 1.E-6
           Xmin = 2*temp2 - Xmid
           Xm2 = temp1 + 1.E-6
           Xmax = 2*temp1 - Xmid
        END IF

        Ymin = 0.
        Ymax = 0.
     END IF

     ALLOCATE(zero(Nx,Ny))
     ALLOCATE(zerolist(MAX_ZERO))
     ALLOCATE(checkzero(MAX_ZERO))

    END SUBROUTINE readdata


 FUNCTION fun(z,w)


    complex(dpk)              :: z, fun, w
    complex(dpk)              :: l1, l2, l3, lz, F1n, F1d, F2n, F2d, F1s, F2s, F1, F2
    complex(dpk)              :: Rnz, Rdz, Rz, R1nz, R1dz, R1z, R2nz, R2dz, R2z


    l1 = sqrt(1._dpk - z*(M1+1._dpk))*sqrt(1._dpk - z*(M1-1._dpk))
    l2 = sqrt(1._dpk*kap1 - z*(kap1*M2+1._dpk))*sqrt(1._dpk*kap1 - z*(kap1*M2-1._dpk))

    print '(A, 1F8.3, A)', ' l1 = (', ABS(l1), ')'

!!$!! Samanta & Freund, JFM, 2008 (3.22) (computes the zeros)

       lz = kap2*l2/l1*(1._dpk-z*M1)*(1._dpk-z*M1)/((1._dpk-z*M2)*(1._dpk-z*M2))
       l3 = sqrt(1._dpk*kap1 - z*(kap1*M3+1._dpk))*sqrt(1._dpk*kap1 - z*(kap1*M3-1._dpk))

       if ((ABS(l2*w) < asymplim .AND. ABS(AIMAG(l2*w)) < asymplim1) .AND. &
            (ABS(l2*w*h) < asymplim .AND. ABS(AIMAG(l2*w*h)) < asymplim1) .AND. &
            (ABS(l1*w*h) < asymplim .AND. ABS(AIMAG(l1*w*h)) < asymplim1)) then
      
          Rdz = (lz*bessj(l1*w*h,m,1)*dhank1(l2*w*h,m,1) - hank1(l2*w*h,m,1)*dbessj(l1*w*h,m,1))* &
               EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h)
          Rnz = (bessj(l2*w*h,m,1)*dbessj(l1*w*h,m,1) - lz*bessj(l1*w*h,m,1)*dbessj(l2*w*h,m,1))* &
               EXP(ABS(AIMAG(l2*w*h)))

!!$              Rnz = lz*bessj(l1*w*h,m,0)*dhank1(l2*w*h,m,0) - hank1(l2*w*h,m,0)*dbessj(l1*w*h,m,0)
!!$              Rdz = bessj(l2*w*h,m,0)*dbessj(l1*w*h,m,0) - lz*bessj(l1*w*h,m,0)*dbessj(l2*w*h,m,0)

          Rz = Rnz/Rdz
          
          F1n = Rz*hank1(l2*w,m,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) + bessj(l2*w,m,1)* &
               EXP(ABS(AIMAG(l2*w)))
          F1d = Rz*dhank1(l2*w,m,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) + dbessj(l2*w,m,1)* &
               EXP(ABS(AIMAG(l2*w)))

!!$              F1n = hank1(l2*w,m,0) + Rz*bessj(l2*w,m,0)
!!$              F1d = dhank1(l2*w,m,0) + Rz*dbessj(l2*w,m,0)
          
          F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/l2*F1n/F1d

       else

          F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/l2*CMPLX(0._dpk,1._dpk,kind=dpk)

       end if

       if (ABS(l3*w) < asymplim .AND. ABS(AIMAG(l3*w)) < asymplim1) then
       
          F2n = hank1(l3*w,m,1)
          F2d = dhank1(l3*w,m,1)

!!$              F2n = hank1(l3*w,m,0)
!!$              F2d = dhank1(l3*w,m,0)

          F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*F2n/F2d

       else

          F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*(8._dpk*l3*w + &
            4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*m*m - &
            CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*l3*w - &
            4._dpk*m*m - 3._dpk)

       end if
      
       fun = w*(F1 - F2)


  END FUNCTION fun  


END PROGRAM nrcomplex
