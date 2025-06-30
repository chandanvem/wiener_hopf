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
  integer                                   :: num_of_inputs,func_case
  complex(dpk), allocatable, dimension(:)   :: input_list !! ''
  integer                                   :: prec !! nag_bessel precision digits
  real(dpk)                                 :: M1, M2, M3, h, kap_T, kap_rho
  real(dpk)                                 :: St_flag
  real(dpk)                                 :: omega_abs
  complex(dpk)                              :: omega
  real(dpk)                                 :: azim_mode
  real(dpk)                                 :: del  !! phase angle in degrees
  real(dpk)                                 :: PI
  real(dpk), parameter                      :: asymplim = 3.27E4
  real(dpk), parameter                      :: asymplim1 = 100
  character(len=200)                        :: dummy

  call readdata                          !! Reads the input file
  call compute_residues


CONTAINS

  SUBROUTINE readdata
  

!! the basic data file:
    real(dpk)           :: temp1, temp2
    integer             :: i
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


      read(10,*)  kap_T
      read(10,*)  kap_rho

      read(10,*) St_flag,dummy
      read(10,*) omega_abs,dummy


      if (St_flag == 1) then
          omega_abs = PI*M1*omega_abs
      end if

      omega = ABS(omega_abs)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*del*PI/180)

      print '(A, F8.4, ",", F8.4,A)', &
                   ' (kap_T,kap_rho) = (', kap_T,kap_rho, ')'

      print '(A, F8.4, ",",F8.4,A)', &
                   ' ( omega St. #, omega helm. #) = (', omega_abs/(PI*M1),omega_abs, ')'


      read(10,*) num_of_inputs,dummy

       if (num_of_inputs > 0) then
         allocate(input_list(num_of_inputs))
         do i = 1, num_of_inputs
           read(10,*) input_list(i)
           print*, input_list(i)
         end do
       end if

     close(10)

    END SUBROUTINE readdata
 
 
    SUBROUTINE compute_residues
        
        integer      :: i
        complex(dpk) :: residue

        do i = 1,num_of_inputs
 
           residue = fun(input_list(i),omega)
           print*, ABS(residue)

        end do

    END SUBROUTINE compute_residues

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


        CASE DEFAULT 

       STOP 'Not a valid function to find roots!! Program Terminated.'
       
    END SELECT


  END FUNCTION fun
!

END PROGRAM check_residue
