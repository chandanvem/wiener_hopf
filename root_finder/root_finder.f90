PROGRAM root_finder

!***! This program computes the zeros of a complex function using a modified
!***! Newton-Raphson technique. It uses a numerical derivative computing
!***! scheme due to Ridders'. This scheme finds multiple zeros by checking
!***! against a database of already-found zeros.

!***! Last modified: Jan 16 2015

!!$  USE nag_bessel_fun, ONLY : nag_bessel_j, nag_bessel_y
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
  complex(dpk)                              :: omega
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

  call mesh                              !! The grid. Make it fine enough to detect all the zeros

  IF(step .NE. 0) CALL teststep(delta,omega)  !! Use it to find an optimum stepsize. Takes a lot of
                                              !! time for a fine mesh

  call limitcase                         !! Prints the limiting values for instability zeros/poles


  call findzero(omega,Nzero,zerolist,checkzero,zero)   !! The main subroutine

  IF(func_case .EQ. 0 .OR. func_case .EQ. 100) THEN

     allocate(alpha(Nzero,3),alphap(Nzero,3))
     CALL findzeroplus(omega,Nzero,Nzero_inc,zerolist,checkzero,zero,alpha,alphap)  
                                                       !! For inc waves
  END IF

  call printinfo(omega,Nzero,zerolist,checkzero,zero)  !! Print & write data

  IF(func_case .EQ. 0 .OR. func_case .EQ. 100) &
       CALL printinfoplus(omega,Nzero,Nzero_inc,alpha,alphap)   
                                                       !! For inc waves




CONTAINS


  SUBROUTINE readdata

!! the basic data file:
    real(dpk)           :: temp1, temp2

    PI = 4._dpk*ATAN(1._dpk)


    open(unit=10, file='input.list', status='old', action='read')

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

      read(10,*) kap1, dummy           
      read(10,*) kap2, dummy           

      PRINT '(A, F8.4, ", ", F8.4,A)', &
                   ' (kap1, kap2) = (', kap1, kap2, ')'


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

  SUBROUTINE limitcase

    complex(dpk)    :: sz1, sz2, sp1, beta

    IF(limit .NE. 0) THEN

       sz1 = 1._dpk/M1

       sz2 = 1._dpk/M2
       
       beta = SQRT((1._dpk - h**2)/h**2)

       sp1 = 1._dpk/(M2**2 + (beta*M1)**2)*(M2**2 + beta**2*M1 - CMPLX(0._dpk,1._dpk,kind=dpk)* &
            beta*(M1 - M2))

       print*,'Low frequency and low Mach no limits:'

       write(*,'(/A4,F15.10,1X,A1,1X,F15.10,A1)'),'sz1:',REAL(sz1),'+',AIMAG(sz1),'i'
       write(*,'(A4,F15.10,1X,A1,1X,F15.10,A1)'),'sz2:',REAL(sz2),'+',AIMAG(sz2),'i'
       write(*,'(A4,F15.10,1X,A1,1X,F15.10,A1/)'),'sp1:',REAL(sp1),'+',AIMAG(sp1),'i'

    END IF

  END SUBROUTINE limitcase


  SUBROUTINE teststep(hopt,w)

!***! This subroutine can be used to find an optimum stepsize for the numerical
!***! derivative. It works by testing step sizes in the range (hmin,hmax) over
!***! the specified grid (X,Y) and then comparing the resulting error measures.
!***! Two such error measures are used, viz., max_error and error_index. The  
!***! former is just the maximum error made over (X,Y) for a given h. The latter
!***! one is, although, supposed to be more useful since it measures how 'deep'
!***! the tableau in Ridders' scheme is traversed for a given h. Early exit from 
!***! the tableau would generally indicate poor extrapolation.

    real(dpk), parameter                  :: hmin = 0.01_dpk, hmax = 10._dpk, dh = 0.01_dpk
    integer                               :: Nh
    real(dpk), allocatable, dimension(:)  :: H
    real(dpk)                             :: hopt
    real(dpk), dimension(Nx,Ny)           :: error, ierror
    complex(dpk)                          :: Z, w
    complex(dpk), dimension(Nx,Ny)        :: testfun
    real(dpk), allocatable, dimension(:)  :: max_error
    integer, allocatable, dimension(:)    :: error_index
    integer, dimension(1)                 :: imin
    integer                               :: i, j, k, optindex

    Nh = (hmax - hmin)/dh + 1
    
    allocate(H(Nh))
    allocate(max_error(Nh))
    allocate(error_index(Nh))
    
    H(1) = hmin
    
    do i = 2, Nh
       H(i) = H(i-1) + dh
    end do
    
    do k = 1, Nh
       do i = 1, Nx
          do j = 1, Ny
             Z = CMPLX(X(i),Y(j),kind=dpk)
             testfun(i,j) = dfun(Z,w,H(k),error(i,j),ierror(i,j))
          end do
       end do
       
       max_error(k) = MAXVAL(error)
       error_index(k) = INT(SUM(ierror))
       
    end do
    
    imin = MAXLOC(error_index(1:Nh))
    optindex = imin(1)

!!$    if (max_error(optindex) > 10.) then
!!$       imin = MINLOC(max_error(1:Nh))
!!$       optindex = imin(1)
!!$    end if

    hopt = H(optindex)

    write(*,'(/A25,F6.2/)'),'Optimum step size:',hopt
    write(*,'(A25,F15.10/)'),'Maximum error value:',max_error(optindex)
    write(*,'(A25,I5/)'),'Error index:',error_index(optindex)

  END SUBROUTINE teststep
    

  SUBROUTINE findzero(w,Nz,zl,cz,zz)

!***! The main subroutine that finds the zeros. It actually calls the Newton-Raphson
!***! method function 'newt' to find the roots. It then collects the zeros and
!***! sorts them in ascending fashion.

    complex(dpk)                               :: w
    complex(dpk), dimension(Nx,Ny)             :: zz
    complex(dpk), dimension(MAX_ZERO)          :: zl, cz
    complex(dpk)                               :: Z, a, check
    real(dpk)                                  :: zeror, zeroi
    integer                                    :: Nz
    character(10)                              :: zflag
    integer                                    :: i, j, k

    zl(:) = (0.,0.)
    Nz = 0

    do i = 1, Nx
       do j = 1, Ny
          Z = CMPLX(X(i),Y(j),kind=dpk)
          zz(i,j) = newt(Z,w,zflag)
          
          if (zflag == 'green') then
             zeror = REAL(zz(i,j))
             zeroi = AIMAG(zz(i,j))

             check = fun(zz(i,j),w)

             if (ABS(REAL(check)) > ZERO_ACC .OR. & 
                  ABS(AIMAG(check)) > ZERO_ACC) zflag = 'red'

             do k = 1,Nz !! Check against the database of zeros already found
                if ((((REAL(zl(k))-ZERO_PREC) < zeror).AND. &
                     (zeror < (REAL(zl(k))+ZERO_PREC))).AND. & 
                     (((AIMAG(zl(k))-ZERO_PREC) < zeroi).AND. &
                     (zeroi < (AIMAG(zl(k))+ZERO_PREC)))) then
                   zflag = 'red' !! The zero found is not unique
                   EXIT
                end if
             end do

          end if

          if (zflag == 'green') then
             if (func_case .GE. 0) print*, "A new zero found!"
             if (func_case .LT. 0) print*, "A new pole found!"
             Nz = Nz + 1
             zl(Nz) = zz(i,j)
             cz(Nz) = check
             
             write(*,'(A4,F15.10,1X,A1,1X,F15.10,A1)'),'at:',&
                  REAL(zz(i,j)),'+',AIMAG(zz(i,j)),'i'

             if (func_case .GE. 0) then
                write(*,'(A11,E16.8,1X,A1,1X,E16.8,A1)'),'func@zero:',&
                     REAL(check),'+',AIMAG(check),'i'
             else
                write(*,'(A11,E16.8,1X,A1,1X,E16.8,A1)'),'func@pole:',&
                     REAL(check),'+',AIMAG(check),'i'
             end if
          end if
          
       end do
    end do

    do j = 2, Nz  !! Sorted in terms of ascending order of the ABS value
       a = zl(j)
       do i = j-1, 1, -1
          if (ABS(zl(i)) <= ABS(a)) EXIT
          zl(i+1) = zl(i)
          zl(i) = a
       end do
    end do


  END SUBROUTINE findzero


  SUBROUTINE findzeroplus(w,Nz,Nzi,zl,cz,zz,al,alp)

    complex(dpk)                               :: w, wdel
    real(dpk), dimension(Nz,3)                 :: al, alp
    complex(dpk), dimension(Nx,Ny)             :: zz, zz1, zz2
    complex(dpk), dimension(MAX_ZERO)          :: zl, cz, zl1, cz1, zl2, cz2
    complex(dpk)                               :: Z, a, check
    real(dpk)                                  :: zeror, zeroi, factor, dwdmu, temp, temp1, temp2
    integer                                    :: Nz, Nzi, zerocount
    character(10)                              :: zflag
    integer                                    :: i, j, k
    
    factor = 1.E-3_dpk
    wdel = ABS(factor)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*del*PI/180)  

    do i = 1, Nzero
       al(i,1) = REAL(zl(i))
       al(i,2) = w*sqrt((1._dpk - REAL(zl(i))*M1)**2 - REAL(zl(i))**2)
       al(i,3) = w*sqrt(kap1**2*(1._dpk - REAL(zl(i))*M2)**2 - REAL(zl(i))**2)
    end do

    CALL findzero(w+wdel,Nz,zl1,cz1,zz1)
    CALL findzero(w-wdel,Nz,zl2,cz2,zz2)

    zerocount = 0
    do i = 1, Nz
       dwdmu = (2._dpk*REAL(wdel))/(REAL(zl1(i)) - REAL(zl2(i)))
       if (func_case .EQ. 0) then
          if (dwdmu > 0) then
             zerocount = zerocount + 1
             alp(zerocount,1) = al(i,1)
             alp(zerocount,2) = al(i,2)
             alp(zerocount,3) = al(i,3)
          end if
       else
          if (dwdmu < 0) then
             zerocount = zerocount + 1
             alp(zerocount,1) = al(i,1)
             alp(zerocount,2) = al(i,2)
             alp(zerocount,3) = al(i,3)
          end if
       end if
    end do

    do j = 2, zerocount  !! Sorted in terms of ascending value of alpha1(or alpha2)
       temp = alp(j,2)
       temp1 = alp(j,3)
       temp2 = alp(j,1)
       do i = j-1, 1, -1
          if (alp(i,2) <= temp) EXIT
          alp(i+1,2) = alp(i,2)
          alp(i+1,3) = alp(i,3)
          alp(i+1,1) = alp(i,1)
          alp(i,2) = temp
          alp(i,3) = temp1
          alp(i,1) = temp2
       end do
    end do

    Nzi = zerocount


  END SUBROUTINE findzeroplus


  SUBROUTINE printinfo(w,Nz,zl,cz,zz)

    complex(dpk)                               :: w
    complex(dpk), dimension(Nx,Ny)             :: zz
    complex(dpk), dimension(MAX_ZERO)          :: zl, cz
    integer                                    :: Nz
    integer                                    :: i, j

    if (Nz == 0) then

       print*,'No zero/pole found in the box!!'

    else

       if (func_case .EQ. 200) then

          if (func_case .GE. 0) then

             if(M1 > 1) then
                write(*,'(/A10)') 'Limits:-->'
                write(*,'(A1,F6.3,A1,F6.3,A1,1X,A3,1X,A1,F6.3,A1,F6.3,A1)') '[',Xmin,',',Xm1,']','AND', &
                     '[',Xm2,',',Xmax,']'
             end if

             write(*,'(/A44)'),'--------------------------------------------'
             write(*,'(A44)'),'              Summary of Zeros              '
             write(*,'(A44)'),'--------------------------------------------'
             do i = 1, Nz
                write(*,'(A5,I4,A1,F15.10,1X,A1,1X,F15.10,A1)'),'ZERO',i,':',&
                     REAL(zl(i)),'+',AIMAG(zl(i)),'i'
             end do
             write(*,'(A44/)'),'--------------------------------------------'

             open(10,file='zerolist.out',form='FORMATTED')
             write(10,'(A2,5F10.5)') '#',M1, M2, M3, m, h
             write(10,'(A2,2F15.5)') '#',w
             if (M1 < 1) then
                write(10,'(A2,2F10.5)') '#',Xmin, Xmax
             else
                write(10,'(A2,4F10.5)') '#',Xmin, Xm1, Xm2, Xmax
             end if
             write(10,'(A2,2F10.5)') '#',Ymin, Ymax
             do i = 1, Nz
                write(10,'(I5,2F20.10,2X,E16.8,2X,E16.8)') i,zl(i),cz(i)
             end do
             close(10)

          else

             write(*,'(/A44)'),'--------------------------------------------'
             write(*,'(A44)'),'              Summary of Poles              '
             write(*,'(A44)'),'--------------------------------------------'
             do i = 1, Nz
                write(*,'(A5,I4,A1,F15.10,1X,A1,1X,F15.10,A1)'),'POLE',i,':',&
                     REAL(zl(i)),'+',AIMAG(zl(i)),'i'
             end do
             write(*,'(A44/)'),'--------------------------------------------'
             
             open(10,file='polelist.out',form='FORMATTED')
             write(10,'(A2,5F10.5)') '#',M1, M2, M3, m, h
             write(10,'(A2,2F15.5)') '#',w
             if (M1 < 1) then
                write(10,'(A2,2F10.5)') '#',Xmin, Xmax
             else
                write(10,'(A2,4F10.5)') '#',Xmin, Xm1, Xm2, Xmax
             end if
             write(10,'(A2,2F10.5)') '#',Ymin, Ymax
             do i = 1, Nz
                write(10,'(I5,2F20.10,2X,E16.8,2X,E16.8)') i,zl(i),cz(i)
             end do
             close(10)

          end if

       else

          if (func_case .GE. 0) then

             write(*,'(/A44)'),'--------------------------------------------'
             write(*,'(A44)'),'              Summary of Zeros              '
             write(*,'(A44)'),'--------------------------------------------'
             do i = 1, Nz
                write(*,'(A5,I4,A1,F15.10,1X,A1,1X,F15.10,A1)'),'ZERO',i,':',&
                     REAL(zl(i)),'+',AIMAG(zl(i)),'i'
             end do
             write(*,'(A44/)'),'--------------------------------------------'

             open(10,file='zerolist.out',form='FORMATTED')
             write(10,'(A2,5F10.5)') '#',M1, M2, M3, m, h
             write(10,'(A2,2F15.5)') '#',w
             write(10,'(A2,2F10.5)') '#',Xmin, Xmax
             write(10,'(A2,2F10.5)') '#',Ymin, Ymax
             do i = 1, Nz
                write(10,'(I5,2F20.10,2X,E16.8,2X,E16.8)') i,zl(i),cz(i)
             end do
             close(10)

          else

             write(*,'(/A44)'),'--------------------------------------------'
             write(*,'(A44)'),'              Summary of Poles              '
             write(*,'(A44)'),'--------------------------------------------'
             do i = 1, Nz
                write(*,'(A5,I4,A1,F15.10,1X,A1,1X,F15.10,A1)'),'POLE',i,':',&
                     REAL(zl(i)),'+',AIMAG(zl(i)),'i'
             end do
             write(*,'(A44/)'),'--------------------------------------------'
             
             open(10,file='polelist.out',form='FORMATTED')
             write(10,'(A2,5F10.5)') '#',M1, M2, M3, m, h
             write(10,'(A2,2F15.5)') '#',w
             write(10,'(A2,2F10.5)') '#',Xmin, Xmax
             write(10,'(A2,2F10.5)') '#',Ymin, Ymax
             do i = 1, Nz
                write(10,'(I5,2F20.10,2X,E16.8,2X,E16.8)') i,zl(i),cz(i)
             end do
             close(10)

          end if

       end if

    end if
    
    open(10,file='log.out',form='FORMATTED')
    do i = 1, Nx
       do j = 1, Ny
          write(10,'(4F20.10)') X(i),Y(j),zz(i,j)
       end do
    end do
    close(10)


  END SUBROUTINE printinfo


  SUBROUTINE printinfoplus(w,Nz,Nzi,al,alp)

    real(dpk), dimension(Nz,3)                 :: al, alp
    complex(dpk)                               :: w
    integer                                    :: Nz, Nzi
    integer                                    :: i, j

    if (Nz == 0) then

       print*,'No zero found in the box!!'

    else

       write(*,'(/A52)'),'----------------------------------------------------'
       write(*,'(A52)'),'   n         mu            alpha            beta    '
       write(*,'(A52)'),'----------------------------------------------------'
       do i = 1, Nz
          write(*,'(I4,1X,F15.10,1X,F15.10,1X,F15.10)'),i,al(i,1),al(i,2),al(i,3)
       end do
       write(*,'(A52/)'),'-----------------------------------------------------'

       open(10,file='zerolist.out',form='FORMATTED')
       write(10,'(A2,5F10.5)') '#',M1, M2, M3, m, h
       write(10,'(A2,2F15.5)') '#',w
       write(10,'(A2,2F10.5)') '#',Xmin, Xmax
       write(10,'(A2,2F10.5)') '#',Ymin, Ymax
       do i = 1, Nz
          write(10,'(I5,3F20.10,2X)') i,al(i,1),al(i,2),al(i,3)
       end do
       close(10)

       print*,''
       if (func_case .EQ. 0) then
          print*,'Considering only the +Z (RIGHT) moving waves.....'
       else
          print*,'Considering only the -Z (LEFT) moving waves......'
       end if

       write(*,'(/A52)'),'----------------------------------------------------'
       write(*,'(A52)'),'   n         mu            alpha            beta    '
       write(*,'(A52)'),'----------------------------------------------------'
       do i = 1, Nzi
          write(*,'(I4,1X,F15.10,1X,F15.10,1X,F15.10)'),i,alp(i,1),alp(i,2),alp(i,3)
       end do
       write(*,'(A52/)'),'-----------------------------------------------------'

       if (func_case .EQ. 0) then
          open(10,file='incident_waves.out',form='FORMATTED')
       else
          open(10,file='reflected_waves.out',form='FORMATTED')
       end if
       write(10,'(A2,5F10.5)') '#',M1, M2, M3, m, h
       write(10,'(A2,2F15.5)') '#',w
       write(10,'(A2,2F10.5)') '#',Xmin, Xmax
       write(10,'(A2,2F10.5)') '#',Ymin, Ymax
       do i = 1, Nzi
          write(10,'(I5,3F20.10,2X)') i,alp(i,1),alp(i,2),alp(i,3)
       end do
       close(10)

    end if


  END SUBROUTINE printinfoplus

  
  FUNCTION fun(z,w)

    !***! This subroutine specifies the function.

    complex(dpk)              :: z, fun, w
    complex(dpk)              :: l1, l2, l3, lz, F1n, F1d, F2n, F2d, F1s, F2s, F1, F2
    complex(dpk)              :: Rnz, Rdz, Rz, R1nz, R1dz, R1z, R2nz, R2dz, R2z


    l1 = sqrt(1._dpk - z*(M1+1._dpk))*sqrt(1._dpk - z*(M1-1._dpk))
    l2 = sqrt(1._dpk*kap1 - z*(kap1*M2+1._dpk))*sqrt(1._dpk*kap1 - z*(kap1*M2-1._dpk))

    SELECT CASE (func_case)

    CASE (0)

!!$!! Samanta & Freund, JFM, 2008 (Acoustic Incident Waves)

       l1 = w*sqrt((1._dpk - z*M1)**2 - z**2)
       l2 = w*sqrt(kap1**2*(1._dpk - z*M2)**2 - z**2)

       F1n = dbessj(h*l2,m,1)*dhank1(l2,m,1)*EXP(ABS(AIMAG(h*l2)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2) - dhank1(h*l2,m,1)*dbessj(l2,m,1)* &
            EXP(ABS(AIMAG(l2)) + CMPLX(0._dpk,1._dpk,kind=dpk)*h*l2)

       F1 = l2*(1._dpk - z*M1)**2*bessj(h*l1,m,1)*EXP(ABS(AIMAG(h*l1)))*F1n

       F2n = bessj(h*l2,m,1)*dhank1(l2,m,1)*EXP(ABS(AIMAG(h*l2)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2) - hank1(h*l2,m,1)*dbessj(l2,m,1)* &
            EXP(ABS(AIMAG(l2)) + CMPLX(0._dpk,1._dpk,kind=dpk)*h*l2)

       F2 = l1*(1._dpk - z*M2)**2*dbessj(h*l1,m,1)*EXP(ABS(AIMAG(h*l1)))*F2n

       fun = F1 - F2

    CASE (100)

!!$!! Samanta & Freund, JFM, 2008 (Acoustic Reflected Waves)

       l1 = w*sqrt((1._dpk - z*M1)**2 - z**2)
       l2 = w*sqrt(kap1**2*(1._dpk - z*M2)**2 - z**2)

       F1n = dbessj(h*l2,m,1)*dhank1(l2,m,1)*EXP(ABS(AIMAG(h*l2)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2) - dhank1(h*l2,m,1)*dbessj(l2,m,1)* &
            EXP(ABS(AIMAG(l2)) + CMPLX(0._dpk,1._dpk,kind=dpk)*h*l2)

       F1 = l2*(1._dpk - z*M1)**2*bessj(h*l1,m,1)*EXP(ABS(AIMAG(h*l1)))*F1n

       F2n = bessj(h*l2,m,1)*dhank1(l2,m,1)*EXP(ABS(AIMAG(h*l2)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2) - hank1(h*l2,m,1)*dbessj(l2,m,1)* &
            EXP(ABS(AIMAG(l2)) + CMPLX(0._dpk,1._dpk,kind=dpk)*h*l2)

       F2 = l1*(1._dpk - z*M2)**2*dbessj(h*l1,m,1)*EXP(ABS(AIMAG(h*l1)))*F2n

       fun = F1 - F2

    CASE (200)

!!$!! Hollow duct incident waves

       l1 = w*sqrt((1._dpk - z*M1)**2 - z**2)

       F1n = dbessj(l1,m,1)*EXP(ABS(AIMAG(l1)))

       fun = F1n

    CASE (1)

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
       
    CASE (-1)
       
!!$!! Samanta & Freund, JFM, 2008 (3.22) (computes the poles)

       lz = kap2*l2/l1*(1._dpk-z*M1)*(1._dpk-z*M1)/((1._dpk-z*M2)*(1._dpk-z*M2))
       l3 = sqrt(1._dpk*kap1 - z*(kap1*M3+1._dpk))*sqrt(1._dpk*kap1 - z*(kap1*M3-1._dpk))

       if ((ABS(l2*w) < asymplim .AND. ABS(AIMAG(l2*w)) < asymplim1) .AND. &
            (ABS(l2*w*h) < asymplim .AND. ABS(AIMAG(l2*w*h)) < asymplim1) .AND. &
            (ABS(l1*w*h) < asymplim .AND. ABS(AIMAG(l1*w*h)) < asymplim1)) then
       
          Rdz = (lz*bessj(l1*w*h,m,1)*dhank1(l2*w*h,m,1) - hank1(l2*w*h,m,1)*dbessj(l1*w*h,m,1))* &
               EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h)
          Rnz = (bessj(l2*w*h,m,1)*dbessj(l1*w*h,m,1) - lz*bessj(l1*w*h,m,1)*dbessj(l2*w*h,m,1))* &
               EXP(ABS(AIMAG(l2*w*h)))
       
          !    Rnz = lz*bessj(l1*w*h,m,0)*dhank1(l2*w*h,m,0) - hank1(l2*w*h,m,0)*dbessj(l1*w*h,m,0)
          !    Rdz = bessj(l2*w*h,m,0)*dbessj(l1*w*h,m,0) - lz*bessj(l1*w*h,m,0)*dbessj(l2*w*h,m,0)
       
          Rz = Rnz/Rdz
       
          F1n = Rz*hank1(l2*w,m,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) + bessj(l2*w,m,1)* &
               EXP(ABS(AIMAG(l2*w)))
          F1d = Rz*dhank1(l2*w,m,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) + dbessj(l2*w,m,1)* &
               EXP(ABS(AIMAG(l2*w)))
       
          !    F1n = hank1(l2*w,m,0) + Rz*bessj(l2*w,m,0)
          !    F1d = dhank1(l2*w,m,0) + Rz*dbessj(l2*w,m,0)
       
          F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/l2*F1n/F1d

       else

          F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/l2*CMPLX(0._dpk,1._dpk,kind=dpk)

       end if

       if (ABS(l3*w) < asymplim .AND. ABS(AIMAG(l3*w)) < asymplim1) then
       
          F2n = hank1(l3*w,m,1)
          F2d = dhank1(l3*w,m,1)
       
          !    F2n = hank1(l3*w,m,0)
          !    F2d = dhank1(l3*w,m,0)
       
          F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*F2n/F2d

       else

          F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*(8._dpk*l3*w + &
               4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*m*m - &
               CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*l3*w - &
               4._dpk*m*m - 3._dpk)

       end if
       
       fun = w*(F1 - F2)
       fun = 1._dpk/fun
       
    CASE (2)

!!$ Gabard & Astley, JFM, 2006 (compute the zeros)

       F1n = dhank1(l1*w*h,m,1)*bessj(l1*w,m,1)*EXP(ABS(AIMAG(l1*w))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l1*w*h) - dbessj(l1*w*h,m,1)*hank1(l1*w,m,1)* &
            EXP(ABS(AIMAG(l1*w*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l1*w)
       F1d = dhank1(l1*w*h,m,1)*dbessj(l1*w,m,1)*EXP(ABS(AIMAG(l1*w))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l1*w*h) - dbessj(l1*w*h,m,1)*dhank1(l1*w,m,1)* &
            EXP(ABS(AIMAG(l1*w*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l1*w)
       F1s = F1n/F1d

       F1 = kap2*(1._dpk - z*M1)*(1._dpk - z*M1)*l2*F1s
       
       F2n = hank1(l2*w,m,1)
       F2d = dhank1(l2*w,m,1)
       F2s = F2n/F2d

       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)*l1*F2s
       
       fun = F1 - F2

    CASE (-2)

!!$ Gabard & Astley, JFM, 2006 (compute the poles)

       F1n = dhank1(l1*w*h,m,1)*bessj(l1*w,m,1)*EXP(ABS(AIMAG(l1*w))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l1*w*h) - dbessj(l1*w*h,m,1)*hank1(l1*w,m,1)* &
            EXP(ABS(AIMAG(l1*w*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l1*w)
       F1d = dhank1(l1*w*h,m,1)*dbessj(l1*w,m,1)*EXP(ABS(AIMAG(l1*w))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l1*w*h) - dbessj(l1*w*h,m,1)*dhank1(l1*w,m,1)* &
            EXP(ABS(AIMAG(l1*w*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l1*w)
       F1s = F1n/F1d

       F1 = kap2*(1._dpk - z*M1)*(1._dpk - z*M1)*l2*F1s
       
       F2n = hank1(l2*w,m,1)
       F2d = dhank1(l2*w,m,1)
       F2s = F2n/F2d

       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)*l1*F2s
       
       fun = F1 - F2
       fun = 1._dpk/fun

    CASE (3)

!! Taylor et al, JSV, 1993 (zeros)

       lz = kap2*l1/l2*(1._dpk-z*M2)*(1._dpk-z*M2)/((1._dpk-z*M1)*(1._dpk-z*M1))
       l3 = sqrt(1._dpk*kap1 - z*(kap1*M3+1._dpk))*sqrt(1._dpk*kap1 - z*(kap1*M3-1._dpk))
       
       Rnz = (lz*bessj(l2*w*h,m,1)*dbessj(l1*w*h,m,1) - bessj(l1*w*h,m,1)*dbessj(l2*w*h,m,1))* &
            EXP(ABS(AIMAG(l2*w*h)))
       Rdz = (bessj(l1*w*h,m,1)*dhank2(l2*w*h,m,1) - lz*dbessj(l1*w*h,m,1)*hank2(l2*w*h,m,1))* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h)

       Rz = Rnz/Rdz
       
       F1n = bessj(l2*w,m,1)*EXP(ABS(AIMAG(l2*w))) + Rz*hank2(l2*w,m,1)* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w)
       F1d = dbessj(l2*w,m,1)*EXP(ABS(AIMAG(l2*w))) + Rz*dhank2(l2*w,m,1)* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w)

       F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/l2*F1n/F1d
       
       F2n = hank2(l3*w,m,1)
       F2d = dhank2(l3*w,m,1)
       F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*F2n/F2d
       
       fun = w*(F1 - F2)

    CASE (-3)

!! Taylor et al, JSV, 1993 (poles)

       lz = kap2*l1/l2*(1._dpk-z*M2)*(1._dpk-z*M2)/((1._dpk-z*M1)*(1._dpk-z*M1))
       l3 = sqrt(1._dpk*kap1 - z*(kap1*M3+1._dpk))*sqrt(1._dpk*kap1 - z*(kap1*M3-1._dpk))
       
       Rnz = (lz*bessj(l2*w*h,m,1)*dbessj(l1*w*h,m,1) - bessj(l1*w*h,m,1)*dbessj(l2*w*h,m,1))* &
            EXP(ABS(AIMAG(l2*w*h)))
       Rdz = (bessj(l1*w*h,m,1)*dhank2(l2*w*h,m,1) - lz*dbessj(l1*w*h,m,1)*hank2(l2*w*h,m,1))* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h)

       Rz = Rnz/Rdz
       
       F1n = bessj(l2*w,m,1)*EXP(ABS(AIMAG(l2*w))) + Rz*hank2(l2*w,m,1)* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w)
       F1d = dbessj(l2*w,m,1)*EXP(ABS(AIMAG(l2*w))) + Rz*dhank2(l2*w,m,1)* &
            EXP(-CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w)

       F1 = (1._dpk-z*M2)*(1._dpk-z*M2)/l2*F1n/F1d
       
       F2n = hank2(l3*w,m,1)
       F2d = dhank2(l3*w,m,1)
       F2 = (1._dpk-z*M3)*(1._dpk-z*M3)/l3*F2n/F2d
       
       fun = w*(F1 - F2)
       fun = 1._dpk/fun

    CASE (4)

!!$ Samanta & Freund, JFM, 2008 (4.10) (zeros)

       F1n = bessj(l1*w*h,m,1)
       F1d = dbessj(l1*w*h,m,1)
       F1s = F1n/F1d

       F1 = kap2*(1._dpk - z*M1)*(1._dpk - z*M1)/l1*F1s

       F2n = dhank1(l2*w,m,1)*bessj(l2*w*h,m,1)*EXP(ABS(AIMAG(l2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) - dbessj(l2*w,m,1)*hank1(l2*w*h,m,1)* &
            EXP(ABS(AIMAG(l2*w))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h)
       F2d = dhank1(l2*w,m,1)*dbessj(l2*w*h,m,1)*EXP(ABS(AIMAG(l2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) - dbessj(l2*w,m,1)*dhank1(l2*w*h,m,1)* &
            EXP(ABS(AIMAG(l2*w))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h)
       F2s = F2n/F2d
       
       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)/l2*F2s
       
       fun = w*(F1 - F2)

    CASE (-4)

!!$ Samanta & Freund, JFM, 2008 (4.10) (poles)

       F1n = bessj(l1*w*h,m,1)
       F1d = dbessj(l1*w*h,m,1)
       F1s = F1n/F1d

       F1 = kap2*(1._dpk - z*M1)*(1._dpk - z*M1)/l1*F1s

       F2n = dhank1(l2*w,m,1)*bessj(l2*w*h,m,1)*EXP(ABS(AIMAG(l2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) - dbessj(l2*w,m,1)*hank1(l2*w*h,m,1)* &
            EXP(ABS(AIMAG(l2*w))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h)
       F2d = dhank1(l2*w,m,1)*dbessj(l2*w*h,m,1)*EXP(ABS(AIMAG(l2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) - dbessj(l2*w,m,1)*dhank1(l2*w*h,m,1)* &
            EXP(ABS(AIMAG(l2*w))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h)
       F2s = F2n/F2d
       
       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)/l2*F2s
       
       fun = w*(F1 - F2)
       fun = 1._dpk/fun

    CASE (11)

!!$ Munt's Duct (zeros)

       F1n = bessj(l1*w,m,1)
       F1d = dbessj(l1*w,m,1)
       F1s = F1n/F1d
       
       F1 = kap2*(1._dpk - z*M1)*(1._dpk - z*M1)*l2*F1s
       
       F2n = hank1(l2*w,m,1)
       F2d = dhank1(l2*w,m,1)
       F2s = F2n/F2d
       
       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)*l1*F2s
       
       fun = F1 - F2

    CASE (-11)

!!$ Munt's Duct (poles)

       F1n = bessj(l1*w,m,1)
       F1d = dbessj(l1*w,m,1)
       F1s = F1n/F1d
       
       F1 = kap2*(1._dpk - z*M1)*(1._dpk - z*M1)*l2*F1s
       
       F2n = hank1(l2*w,m,1)
       F2d = dhank1(l2*w,m,1)
       F2s = F2n/F2d
       
       F2 = (1._dpk - z*M2)*(1._dpk - z*M2)*l1*F2s
       
       fun = F1 - F2
       fun = 1._dpk/fun


    CASE (101)

!! Test functions: Supersonic kernel (2nd part)  ! No density ratio

       F1s = dbessj(l1*w*h,m,1)*EXP(ABS(AIMAG(l1*w*h)))
       F2s = bessj(l1*w*h,m,1)*EXP(ABS(AIMAG(l1*w*h)))

       F1n = (dhank1(l2*w,m,1)*bessj(l2*w,m,1) - &
            dbessj(l2*w,m,1)*hank1(l2*w,m,1))*EXP(ABS(AIMAG(l2*w)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w)  ! this term has no zeros; do not include

       F1 = F1s

       F1d = lz*F2s*(dhank1(l2*w,m,1)*dbessj(l2*w*h,m,1)*EXP(ABS(AIMAG(l2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) - dbessj(l2*w,m,1)*dhank1(l2*w*h,m,1)* &
            EXP(ABS(AIMAG(l2*w)) + CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h))
       
       F2d = F1s*(dhank1(l2*w,m,1)*bessj(l2*w*h,m,1)*EXP(ABS(AIMAG(l2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) - dbessj(l2*w,m,1)*hank1(l2*w*h,m,1)* &
            EXP(ABS(AIMAG(l2*w)) + CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h))

       F2 = F1d - F2d
       
       fun = F1/F2
       
    CASE (-101)

!! Test functions: Supersonic kernel (2nd part)  ! No density ratio

       F1s = dbessj(l1*w*h,m,1)*EXP(ABS(AIMAG(l1*w*h)))
       F2s = bessj(l1*w*h,m,1)*EXP(ABS(AIMAG(l1*w*h)))

       F1n = (dhank1(l2*w,m,1)*bessj(l2*w,m,1) - &
            dbessj(l2*w,m,1)*hank1(l2*w,m,1))*EXP(ABS(AIMAG(l2*w)) + &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w)  ! this term has no zeros; do not include

       F1 = F1s

       F1d = lz*F2s*(dhank1(l2*w,m,1)*dbessj(l2*w*h,m,1)*EXP(ABS(AIMAG(l2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) - dbessj(l2*w,m,1)*dhank1(l2*w*h,m,1)* &
            EXP(ABS(AIMAG(l2*w)) + CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h))

       F2d = F1s*(dhank1(l2*w,m,1)*bessj(l2*w*h,m,1)*EXP(ABS(AIMAG(l2*w*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w) - dbessj(l2*w,m,1)*hank1(l2*w*h,m,1)* &
            EXP(ABS(AIMAG(l2*w)) + CMPLX(0._dpk,1._dpk,kind=dpk)*l2*w*h))
       
       F2 = F1d - F2d
       
       fun = F1/F2
       fun = 1._dpk/fun

    CASE DEFAULT 

       STOP 'Not a valid function to find roots!! Program Terminated.'
       
    END SELECT


  END FUNCTION fun


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

    open(20,file='mesh.x')
    do i = 1, Nx
       write(20,'(I10,2X,F20.10)'), i, X(i)
    end do
    close(20)
    
    open(20,file='mesh.y')
    do i = 1, Ny
       write(20,'(I10,2X,F20.10)'), i, Y(i)
    end do
    
  END SUBROUTINE mesh


  FUNCTION newt(gz,w,flg)

!***! The Newton-Raphson routine for complex functions using derivatives. The
!***! algorithm removes zeros from the function as soon as it is found, thereby
!***! preventing duplication and speeding up the process. This also enables
!***! to find multiple zeros, if any.

    complex(dpk)       :: newt, gz, w
    complex(dpk)       :: f, df, dz, zp, zq, zw
    real(dpk)          :: newtr, newti
    real(dpk)          :: err1, err2
    character(10)      :: flg
    integer            :: i, j, k


    newt = gz
    flg = 'red' !! 'Green' flag means zero found
    do i = 1, MAX_ITE
       f = fun(newt,w)
!       print*,f, newt
       df = dfun(newt,w,delta,err1,err2) !! Numerical derivative
!       print*, err1

       dz = f/df
       newt = newt - dz

       newtr = REAL(newt)
       newti = AIMAG(newt)

       if(newtr > Xmax+ZERO_ACC/10 .OR. newtr < Xmin-ZERO_ACC/10 .OR. &
            newti > Ymax+ZERO_ACC/10 .OR. newti < Ymin-ZERO_ACC/10) then
!          print*, "Jumped out of bounds: may be approaching a pole"
          RETURN
       end if

       if (ABS(REAL(dz)) < ZERO_ACC .AND. ABS(AIMAG(dz)) < ZERO_ACC) then 
          flg = 'green'
          RETURN
       end if

    end do

    print*, "WARNING: Max Iterations Exceeded!" !! Failure to achieve accuracy

  END FUNCTION newt



END PROGRAM root_finder
