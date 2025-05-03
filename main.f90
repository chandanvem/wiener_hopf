PROGRAM main

!!$ Single-piece integration contour (algebraic) used
!!$ Comments: 1. Incident vorticity: inner instability inside contour
!!$           2. Suitable for incident acoustic/instability mode
!!$           3. Restart mode added
!!$           4. Supersonic modes added
!!$           5. Option to compute asymptotic far field
!!$           6. Vortex shedding parameter
!!$ Last modified: APR-07-2015

!!$ USE nag_bessel_fun, ONLY : nag_bessel_j, nag_bessel_y

  IMPLICIT none

  integer, parameter                         :: dpk = kind(1.0d0)  !! double precision kind
  real(dpk)                                  :: tol  !! specified tolerance level
  complex(dpk), allocatable, dimension(:)    :: initpoints  !! location of the starting pts
  complex(dpk), allocatable, dimension(:)    :: intpanel  !! integral value at each panel
  integer, allocatable, dimension(:)         :: Npanel  !! number of times a panel is divided
  integer                                    :: max_points !! number of starting pts for kernel contour
  integer                                    :: N1, totinitpts, totiftpts
  integer                                    :: prswitch, restart, reflswitch, farswitch, vortswitch
  integer                                    :: Nzbc, Npbc  !! number of zeros & poles in between s_1^- & s_2^-
  integer                                    :: Nzsp, Npsp  !! number of supersonic zeros & poles
  integer                                    :: Nphi  !! polar mesh resolution for directivity computation
  complex(dpk), allocatable, dimension(:)    :: szi  !! the points defining the IFT contour
  complex(dpk), allocatable, dimension(:)    :: szk  !! the points defining the Kernel contour
  complex(dpk), allocatable, dimension(:)    :: szbc, spbc, szsp, spsp
  complex(dpk), allocatable, dimension(:)    :: fplusz, fplusz_temp
  complex(dpk)                               :: kmmup, kpsz1, kpsz2, kpsp1
  complex(dpk), allocatable, dimension(:)    :: kzsp
  real(dpk), allocatable, dimension(:)       :: phi  !! directivity stuff
  complex(dpk), allocatable, dimension(:)    :: Dmn  !! directivity stuff
  real(dpk)                                  :: phii, phif  !! start and end angles for phi
  real(dpk)                                  :: cnangle   !! largest cone angle for instability modes
  real(dpk), allocatable, dimension(:)       :: cnanglei  !! cone angles for instability modes
  complex(dpk)                               :: alpha1, alpha2
  real(dpk)                                  :: PI
  real(dpk)                                  :: M1, M2, M3
  real(dpk)                                  :: h
  real(dpk)                                  :: kap1  !! sqrt(T1/T0)
  real(dpk)                                  :: kap2  !! rho1/rho0
  real(dpk)                                  :: theta, crpt, w0, delr, deli, offset
  complex(dpk)                               :: wi, wr
  real(dpk)                                  :: Rmax, Rmin, Zmin, Zmax !! dimensions of the physical domain
  real(dpk), allocatable, dimension(:)       :: R, Z
  real(dpk)                                  :: Zo  !! start point of incident vorticity (negative)
  integer                                    :: Nmeshr, Nmeshz !! 400X400 -> Seg Fault!!
  real(dpk)                                  :: circmod
  real(dpk)                                  :: vparm  !! Vortex shedding parameter
  complex(dpk)                               :: sz1 !! instability zero 1
  complex(dpk)                               :: sz2 !! instability zero 2
  complex(dpk)                               :: sp1 !! instability pole 1
  complex(dpk)                               :: mup !! axial wave no (of incident wave)
  complex(dpk)                               :: mu0 !! incident axial wave no (for incident vort)
  complex(dpk)                               :: resp!! portion of residue (for inc vort)
  complex(dpk)                               :: psi
  integer                                    :: Nmax !! number of points for IFT contour
  complex(dpk), allocatable, dimension(:)    :: iftpoints
  complex(dpk), allocatable, dimension(:,:)  :: pressure
  complex(dpk), allocatable, dimension(:,:)  :: acoupressure,totpressure,inspressure1,inspressure2,initpressure
  real(dpk)                                  :: asymplim, asymplim1 !! asymptotic limit for the kernel functions


!! Only two calls:

  call initialize  !! this subrtn reads and initializes data
  call solve  !! this subrtn is the main solver


CONTAINS

  SUBROUTINE initialize

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Read data
!! 2. Compute the kernel and IFT integration contour points and various contour params
!! 3. Find the *first* kernel crossover point
!! 4. Check if the instability zeros and pole lie appropriately wrt the contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer             :: i1, i2, i3
    real(dpk)           :: t, ai, ai_im, af, af_im
    real(dpk)           :: aki, aki_im, akf, akf_im
    character(60)       :: pos

!! the basic data file:

    open(1,file='input.list.p')

    read(1,*) vortswitch  !! 1 = Use incident vorticity mode; 2 = First sup ins mode; Else = acoustic mode
    print *, 'vortswitch=', vortswitch

    read(1,*) M1  !! core jet Mach number
    print *, 'M1=', M1

    read(1,*) M2  !! coflow Mach number
    print *, 'M2=', M2

    read(1,*) M3  !! ambient flow Mach number
    print *, 'M3=', M3

    read(1,*) h   !! the width of the core jet, h = Ri/Ro
    print *, 'width of the core jet h=', h

    if ((vortswitch == 1) .OR. (vortswitch == 2)) then
       read(1,*) Zo  !! starting point of inc instability (-ve)
    end if

    read(1,*) w0  !! the Helmholtz number
    print *, 'Helmholtz number omega=', w0

    read(1,*) kap1  !! sqrt(temp ratio)
    print *, 'Sqrt of temperature ratio=', kap1
 
    read(1,*) kap2  !! density ratio
    print*, 'Density ratio=', kap2

    read(1,*) circmod  !! the circumferential mode no
    print*, 'Azimuthal wavenumber (circmod)=', circmod

    read(1,*) sz1
    print*, 'sz1 (First instability zero) =', sz1


    read(1,*) sz2
    print*, 'sz2 (Second instability zero) =', sz2


    read(1,*) sp1
    print*, 'sp1 (Instability pole) =', sp1


    if ((vortswitch == 1) .OR. (vortswitch == 2)) then
       read(1,*) mu0  !! upstream inc acoustic mode
       print*, 'Incident Vorticity mode activated'
    else
       read(1,*) mup
       print*, '================== Non-incident vorticity mode ====================='
       print*,''
       print*, 'mu for the incident mode =',mup
    end if

    read(1,*) offset  !! the offset between the two integration contours
    print*,'Offset between the IFT and K split contours=', offset

 
    read(1,*) tol  !! the tolerance of the adaptive contour integration routine
    print*,'Tolerance for adaptive contour integration =', tol

    read(1,*) Nzbc  !! zeros needed to be removed from the kernel
    print*,'Number of zeros to be removed from kernel =', Nzbc

    read(1,*) Npbc  !! poles needed to be removed from the kernel
    print*,'Number of poles to be removed from kernel =', Npbc

    read(1,*) Nzsp  !! supersonic zeros
    print*,'Number of supersonic zeros =', Nzsp

    read(1,*) Npsp  !! supersonic poles
    print*,'Number of supersonic poles =', Npsp

    if (vortswitch == 2) then
       if (Npsp == 0) then
          print*, "For vortswitch mode 2, at least one upstream supersonic pole needed! Exiting..."
          STOP
       end if
    end if

    print*,''
    print*,'','====== Mesh and contour parameters ======',''
    
    read(1,*) max_points
    print*,'Number of kernel points in each loop =', max_points

    read(1,*) theta  !! the stretching parameter for kernel contour meshpoints
    print*,'Stretching parameter for kernel contour =', theta

    read(1,*) Nmax
    print*,'Number of IFT points in each loop Nmax =', Nmax

    read(1,*) Rmin
    read(1,*) Rmax
    read(1,*) Zmin
    read(1,*) Zmax
   
    read(1,*) Nmeshr
    read(1,*) Nmeshz

    print*,'Dimensions of domain in R = [',Rmin,Rmax,']' 
    print*,'Dimensions of domain in Z = [',Zmin,Zmax,']' 
    print*,'Nr x Nz =',Nmeshr,'x',Nmeshz

    print*,'=============================='

    read(1,*) asymplim
    print*,'Asymptotic limit of ??? = ',asymplim

    read(1,*) asymplim1
    print*,'Asymptotic limit of ??? = ',asymplim1

    read(1,*) vparm  !! Default = 1.0 (0,1)
    print*,'Vortex shedding parameter gamma = ',vparm 
      

!!============================
    read(1,*) prswitch  !! 0 = Potential; 1 = Pressure

    if (prswitch == 0) then
       print*,'Solution in potential mode'
    elseif (prswitch == 1) then
       print*,'Solution in pressure mode '
    end if
!!============================
    read(1,*) reflswitch  !! 1 = Reflection mode: incident mode not added
    if (reflswitch == 1) then
       print*,'Reflection mode: incident mode not added'
    else 
       print*,'Reflection mode: incident mode added'
    end if
!!============================
    read(1,*) farswitch  !! 1 = Far-field mode: compute directivity; 2 = 1 + nearfield of sup zeros in polar mode

    if ((farswitch == 1) .OR. (farswitch == 2)) then
       read(1,*) Nphi
       allocate(cnanglei(Nzsp+2))
    end if

    if ((farswitch == 2) .AND. Nzsp <= 0) then
       print*, "For farswitch mode 2, non-zero number of supersonic zero needed! Exiting..."
       STOP
    end if

!!=====================================================

    read(1,*) restart  !! restart status: 0 = fresh job, i.e., no "fplus_part.out" exists

!! read the zeros & poles that need to be excluded, if any:
!! USE: when the effect of particular modes need to be studied individually
    
    if (Nzbc > 0) then
        allocate(szbc(Nzbc))
        do i1 = 1, Nzbc
           read(1,*) szbc(i1)
!          print*, szbc(i1)
        end do
    end if


    if (Npbc > 0) then
       allocate(spbc(Npbc))

       do i1 = 1, Npbc
          read(1,*) spbc(i1)
!         print*, spbc(i1)
       end do
    end if

    if (Nzsp > 0) then
        allocate(szsp(Nzsp))
        do i1 = 1, Nzsp
          read(1,*) szsp(i1)
!         print*, szsp(i1)
        end do
   end if

   
   if (Npsp > 0) then 
       allocate(spsp(Npsp))
       do i1 = 1, Npsp
         read(1,*) spsp(i1)
!        print*, spsp(i1)
       end do
  end if

    close(1)

!!================================================================


    PI = 4._dpk*ATAN(1.)
    delr = 0._dpk   !! for real omega
    deli = PI/2._dpk  !! for purely imaginary omega

!! defining the omegas (Helmholtz numbers):

    wr = ABS(w0)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*delr)  !! real
!!$    wi = ABS(w0)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*deli)  !! imaginary
!!$    wi = wr

!! whether to compute velocity potential or pressure:

    if (prswitch == 0) then
       if ((farswitch == 1) .OR. (farswitch == 2)) then
          print*,''
          print*,'Far-field is computed ONLY in pressure mode. Exiting...'
          STOP
       else
          print*,''
          print*,'Computing Potential...'
       end if
    else
       print*,''
       print*,'Computing Pressure...'
    end if

 !! now allocate for the contour points:

    allocate(szi(5))
    allocate(szk(5))

!! read the contour points, only IFT contour is read in full;
!! the end points of the kernel integration contour read while
!! the rest are obtained via the "offset" parameter
!! szi(1), szi(5) are end points; szi(2), szi(4) are where the
!! algebraic contour starts & ends; szi(3) is where it crosses
!! real axis

    open(1,file='input.iftpts')

    read(1,*) t
    ai = t
    read(1,*) szi(2)
    read(1,*) t
    szi(3) = CMPLX(t,0._dpk,kind=dpk)
    read(1,*) szi(4)
    read(1,*) t
    af = t

    call yfunc_p(ai,ai_im,REAL(szi(3)),REAL(szi(2)-szi(3)),AIMAG(szi(2)-szi(3)))
    call yfunc_p(af,af_im,REAL(szi(3)),REAL(szi(4)-szi(3)),AIMAG(szi(4)-szi(3)))

    szi(1) = CMPLX(ai,ai_im,kind=dpk)  
    szi(5) = CMPLX(af,af_im,kind=dpk)

    read(1,*) t
    aki = t
    read(1,*) t
    akf = t

    close(1)

!! now determine the kernel integration contour points:
!! NOTE: the kernel contour lies above the IFT one and thus it's further away from the real
!! axis when above, but closer to it when below

    szk(2) = szi(2) + CMPLX(0._dpk,offset,kind=dpk)  !! Kernel points
    szk(3) = szi(3) + CMPLX(5._dpk*offset,0._dpk,kind=dpk)
    szk(4) = szi(4) + CMPLX(0._dpk,offset,kind=dpk)

    call yfunc_p(aki,aki_im,REAL(szk(3)),REAL(szk(2)-szk(3)),AIMAG(szk(2)-szk(3)))
    call yfunc_p(akf,akf_im,REAL(szk(3)),REAL(szk(4)-szk(3)),AIMAG(szk(4)-szk(3)))
    
    szk(1) = CMPLX(aki,aki_im,kind=dpk)
    szk(5) = CMPLX(akf,akf_im,kind=dpk)

!! find where the contour crosses the real axis (at "crpt"):

    crpt = REAL(szk(3))

!! here both the instability zeros and the pole need to lie outside (under) the contour,
!! since we use the residue theorem to compute them anyway:

    call checkloc(sz1,i1) 
    call checkloc(sz2,i2)
    call checkloc(sp1,i3)

!! if any of them lie inside ask to redefine the contour:

    if (vortswitch .EQ. 0) then

       if(i1==0 .OR. i2==0 .OR. i3==0) then
          print*,''
          print*,'Redefine the contour: One of instability zero/pole is INSIDE!'
          print*,''
          if (i1==0) print*, 'Zero1 is inside'
          if (i2==0) print*, 'Zero2 is inside'
          if (i3==0) print*, 'Pole1 is inside'
          print*,''
          STOP
       end if

    else

       if(i1==1 .OR. i2==0 .OR. i3==1) then
          print*,''
          print*,'Redefine the contour:'
          print*,''
          if (i1==1) print*, 'Zero1 is outside'
          if (i2==0) print*, 'Zero2 is inside'
          if (i3==1) print*, 'Pole1 is outside'
          print*,''
          STOP
       end if

    end if

    do i1 = 1, Nzsp
       call checkloc(szsp(i1),i2)
       if (vortswitch .EQ. 0) then
          if(i2 == 0) then
             print*,''
             print*,'Redefine the contour: The following supersonic zero is inside:'
             print*,szsp(i1)
          end if
       else
          if(i2 == 1) then
             print*,''
             print*,'Redefine the contour: The following supersonic zero is outside:'
             print*,szsp(i1)
          end if
       end if
    end do

    do i1 = 1, Npsp
       call checkloc(spsp(i1),i2)
       if (vortswitch .EQ. 0) then
          if(i2 == 0) then
             print*,''
             print*,'Redefine the contour: The following supersonic pole is inside:'
             print*,spsp(i1)
          end if
       else
          if(i2 == 1) then
             print*,''
             print*,'Redefine the contour: The following supersonic pole is outside:'
             print*,spsp(i1)
          end if
       end if
    end do

!! compute the various contour parameters:

    N1 = max_points  !! for each of the two kernel contour loops

    totinitpts = 2*N1 + 3  !! total number of kernel contor pts

    totiftpts = 2*Nmax + 3  !! total number of IFT contor pts

!! print the contour points:

    write(*,'(/A26)'),'The Integration Contours:'
    write(*,'(/A27)'),'1. The kernel integration:'
    write(*,'(/A14)'),'key points:->'
    do i1 = 1, 5
       write(pos,"(I2,' =')"),i1
       write(*,'(/1X,A5,2F15.6)'), 'i'//pos, szk(i1)
    end do

    write(*,'(/A30)'),'2. The inv Fourier transform:'
    write(*,'(/A14)'),'key points:->' 
    do i1 = 1, 5
       write(pos,"(I2,' =')"),i1
       write(*,"(/1X,A5,2F15.6)"), 'i'//pos, szi(i1) 
    end do
    
    write(*,"(/A6,I6/)"),'N =',totiftpts

    if ((farswitch == 1) .OR. (farswitch == 2)) call initfarfield


  END SUBROUTINE initialize


  SUBROUTINE solve

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Calls the various other subroutines which do the actual job
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer     :: index

    call definecontours

    print*,''
    print*,'Starting...'
    print*,''
    print*,'Computing K+ at the zeros and K- at mu+:'

    call precompute

    if ((farswitch == 1) .OR. (farswitch == 2)) then

       print*,''
       print*,'Now computing Dmn:'
       print*,''

       call computefarfield

    else

       call meshgrid
    
       print*,''
       print*,'Now computing F+:'
       print*,''

       if (restart .NE. 0) then
          call readfplus(index)  !! assumes the fplus_part.out file to be present
          call computefplus(restart,index)
       else
          call computefplus(restart,0)
       end if

       print*,''
       print*,'Now starting the IFT:'
       print*,''

       call computeift

    end if

    print*,'Done!'


  END SUBROUTINE solve


  SUBROUTINE initfarfield

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Initializes the far-field computations
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer            :: j

    allocate(phi(Nphi))
    allocate(Dmn(Nphi))

    ! 0 < phi < 180

    phii = 0.5_dpk*PI/180._dpk
    phif = 179.5_dpk*PI/180._dpk

    ! Generate the directivity mesh

    do j = 1, Nphi
       phi(j) = phii + (j-1)*(phif-phii)/(Nphi-1)
    end do

    ! Find the instability wave cone angles
    
    cnanglei(1) = ATAN(- (AIMAG(sz1))/(AIMAG(SQRT(kap1**2*(1._dpk - sz1*M3)**2 - sz1**2))))
    cnanglei(2) = ATAN(- (AIMAG(sz2))/(AIMAG(SQRT(kap1**2*(1._dpk - sz2*M3)**2 - sz2**2))))

    do j = 3, 2+Nzsp
       cnanglei(j) = ATAN(- (AIMAG(szsp(j-2)))/(AIMAG(SQRT(kap1**2*(1._dpk - szsp(j-2)*M3)**2 -&
            szsp(j-2)**2))))
    end do

!!$    print*,'coneangle1=',cnangle1,'coneangle2=',cnangle2

    cnangle = cnanglei(1)
    do j = 2, Nzsp+1
       if (cnanglei(j) >= cnangle) then
          cnangle = cnanglei(j)
       end if
    end do

!!$    print*,'coneangle=',coneangle*180./PI

    open(10,file='coneangles.out',form='FORMATTED')
    write(10,'(A5,F20.5)') 'max:', cnangle*180._dpk/PI
    do j = 1, Nzsp+2
       write(10,'(I3,F20.5)') j, cnanglei(j)*180._dpk/PI
    end do
    close(10)

    if (farswitch == 2) call meshgrid


  END SUBROUTINE initfarfield


  SUBROUTINE computefarfield

    integer                    :: i, k, jj

    jj = 1

    if (restart .NE. 0) then

       open(10,file='directivity_temp.out',form='FORMATTED')
       do i = 1, Nphi
          read(10,'(I5,2F20.15)',end = 100) k, Dmn(i)
       end do
       close(10)

100    do i = 1, k

          write(*,'(/A10,2X,F15.10,A1)'),'angle :: ', phi(i)*180._dpk/PI,':'
          write(*,'(A10,2X,2F15.10/)'),'Dmn =', Dmn(i)
          if (mod(i,10) == 0) write(*,'(A64)'),'----------------------------------------------------------------'

          jj = k+1

       end do

    end if

    do i = jj, Nphi
       
       write(*,'(/A10,2X,F15.10,A1)'),'angle :: ', phi(i)*180._dpk/PI,':'

       call directcomp(phi(i),Dmn(i))
       
       write(*,'(A10,2X,2F15.10/)'),'Dmn =', Dmn(i)

       if (i == 1) then
          open(10,file='directivity_temp.out',form='FORMATTED',status='UNKNOWN')
       else
          open(10,file='directivity_temp.out',form='FORMATTED',position='APPEND')
       end if
       write(10,'(I5,2F20.15)') i, Dmn(i)
       close(10)

       if (mod(i,10) == 0) write(*,'(A64)'),'----------------------------------------------------------------'
    end do

    open(10,file='directivity_ins.out',form='FORMATTED')
    do i = 1, Nphi
       if (phi(i) .LE. cnangle) write(10,'(F8.2,2F20.15)') phi(i)*180._dpk/PI, Dmn(i) 
    end do
    close(10)

    open(10,file='directivity.out',form='FORMATTED')
    do i = 1, Nphi
       if (phi(i) > cnangle) write(10,'(F8.2,2F20.15)') phi(i)*180._dpk/PI, Dmn(i)
    end do
    close(10)

    if (vortswitch .EQ. 0) then
       if (farswitch == 2) call residuepolarsup
    end if


  END SUBROUTINE computefarfield


  SUBROUTINE definecontours

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Computes the panel lengths for the integration contours
!! 2. Calls subrtn "initcontour" which does the the actual contour defining
!! 3. Writes the contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                     :: panel_len1, panel_len2  !! length of each panel
    integer                       :: i

    allocate(initpoints(totinitpts))

!! length of each panel (kernel contour):

    panel_len1 = (REAL(szk(3))-REAL(szk(1)))/(N1+1)  !! panel length in the left section
    panel_len2 = (REAL(szk(5))-REAL(szk(3)))/(N1+1)  !! panel length in the right section

    call initcontour(szk,panel_len1,panel_len2,N1,1,initpoints)

    open(10,file='initialpoints.out',form='FORMATTED')
    do i = 1,totinitpts
       write(10,'(I10,2F30.20)') i, initpoints(i)
    end do
    close(10)

    allocate(iftpoints(totiftpts))

!! length of each panel (IFT contour):

    panel_len1 = (REAL(szi(3))-REAL(szi(1)))/(Nmax+1)
    panel_len2 = (REAL(szi(5))-REAL(szi(3)))/(Nmax+1)

    call initcontour(szi,panel_len1,panel_len2,Nmax,2,iftpoints)

    open(10,file='iftpoints.out',form='FORMATTED')
    do i = 1,totiftpts
       write(10,'(I10,2F30.20)') i, iftpoints(i)
    end do
    close(10)


  END SUBROUTINE definecontours


  SUBROUTINE initcontour(p,l1,l2,Ns,sw,inpts)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Initialize the integration contours
!! 2. First compute the individual segments
!! 3. Combine the individual segments to get the full contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), allocatable, dimension(:)      :: initpnts1, initpnts2
    complex(dpk), dimension(2*Ns+3)              :: inpts
    complex(dpk), dimension(5)                   :: p
    real(dpk)                                    :: l1, l2
    integer                                      :: Ns, sw


    allocate(initpnts1(Ns+2))
    allocate(initpnts2(Ns+2))

    call initsubdiv(REAL(p(1)),p(2)-REAL(p(3)),REAL(p(3)),Ns,l1,sw,1,initpnts1)
    call initsubdiv(REAL(p(3)),p(4)-REAL(p(3)),REAL(p(5)),Ns,l2,sw,2,initpnts2)

    call combine(Ns,initpnts1,initpnts2,inpts)

  END SUBROUTINE initcontour


  SUBROUTINE initsubdiv(p,q,r,N,len,sw,ss,zi)
    
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Define the different segments of the integration contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
    
    complex(dpk), dimension(N+2)    :: zi  !! each segment is N+2 long
    real(dpk), dimension(N+2)       :: xi
    real(dpk), dimension(N+2)       :: yi
    complex(dpk)                    :: q
    real(dpk)                       :: p, r, len
    integer                         :: N, i, ss, sw


    if (sw == 1) then  !! indicates kernel contour; equispaced points are a waste
                       !! use some stretching of the grids
       select case (ss)

       case(1)
          call stretchx(r,p,N+2,1,xi)  

       case(2)
          call stretchx(p,r,N+2,2,xi)  

       end select

    else
       xi(1) = REAL(p)

    end if


    do i = 1,N+1
       select case (ss)

       case(1)
          call yfunc_p(xi(i),yi(i),r,REAL(q),AIMAG(q))  !! the imaginary points
          if (sw == 2) xi(i+1) = xi(i) + len  !! for ift contour the points are equispaced
          
       case(2)
          call yfunc_p(xi(i),yi(i),p,REAL(q),AIMAG(q))
          if (sw == 2) xi(i+1) = xi(i) + len

       end select

       zi(i) = CMPLX(xi(i),yi(i),kind=dpk)
    end do


    select case (ss)  !! the right end point

    case(1)
       call yfunc_p(xi(N+2),yi(N+2),r,REAL(q),AIMAG(q))
       
    case(2)
       call yfunc_p(xi(N+2),yi(N+2),p,REAL(q),AIMAG(q))

    end select

    zi(N+2) = CMPLX(xi(N+2),yi(N+2),kind=dpk)

  END SUBROUTINE initsubdiv
  

  SUBROUTINE combine(Np,init1,init2,init)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Combine two successive segments
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), dimension(2*Np+3)                  :: init
    complex(dpk), dimension(Np+2)                    :: init1, init2
    integer                                          :: i
    integer                                          :: Np

    init = (0._dpk,0._dpk)

    do i = 1, Np+2
       init(i) = init1(i)
    end do

    do i = Np+3, 2*Np+3
       init(i) = init2(i-Np-1)
    end do

  END SUBROUTINE combine


  SUBROUTINE stretchx(a,b,N,ss,xst)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Define the stretching function needed for the kernel contour
!! 2. Points are to be clustered near the crossover location
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                      :: a, b
    real(dpk)                      :: xi, xf, dx
    real(dpk)                      :: c
    integer                        :: N, i, ss
    real(dpk), dimension(N)        :: xbar, xst

    c = 1._dpk   !! a constant

    xi = 0._dpk
    xf = 1._dpk
    dx = (xf - xi)/(N-1)

    xbar(1) = xi  !! xbar takes uniform spaced values between 0 and 1
    do i = 1, N-1
       xbar(i+1) = xbar(i) + dx
    end do

    if (ss == 2) then  !! right hand part of the contour
       do i = 1, N
          xst(i) = (b-a)*(1._dpk + TANH(theta*(xbar(i)-c))/TANH(theta*c)) + a
       end do
    else  !! left hand part of the contour
       do i = 1, N
          xst(N+1-i) = (b-a)*(1._dpk + TANH(theta*(xbar(i)-c))/TANH(theta*c)) + a
       end do
    end if


  END SUBROUTINE stretchx


  SUBROUTINE meshgrid

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Define the physical mesh
!! 2. Write the mesh in PLOT3D format
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: dsz, dsr
    integer       :: i, j

    dsr = (Rmax - Rmin)/(Nmeshr - 1)
    dsz = (Zmax - Zmin)/(Nmeshz - 1)

    allocate(R(Nmeshr))
    allocate(Z(Nmeshz))

    R(1) = Rmin
    R(Nmeshr) = Rmax
    Z(1) = Zmin
    Z(Nmeshz) = Zmax

    do i = 1, Nmeshr-2
       R(i+1) = R(i) + dsr
    end do
    do i = 1, Nmeshz-2
       Z(i+1) = Z(i) + dsz
    end do

    open(1,file='mesh.out',form='UNFORMATTED')
    write(1) Nmeshz,Nmeshr
    write(1) ((Z(i),i=1,Nmeshz),j=1,Nmeshr),((R(j),&
         i=1,Nmeshz),j=1,Nmeshr)
    close(1)   

    open(20,file='mesh.r')
    do i = 1, Nmeshr
       write(20,'(I10,2X,F20.10)'), i, R(i)
    end do
    close(20)

    open(20,file='mesh.z')
    do i = 1, Nmeshz
       write(20,'(I10,2X,F20.10)'), i, Z(i)
    end do
    close(20)

    if (farswitch == 2) then
       open(1,file='mesh_polar.out',form='UNFORMATTED')
       write(1) Nphi,Nmeshr
       write(1) ((phi(i),i=1,Nphi),j=1,Nmeshr),((R(j),&
            i=1,Nphi),j=1,Nmeshr)
       close(1)   
    end if


  END SUBROUTINE meshgrid


  SUBROUTINE directcomp(iphi,dmn)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the asymptotic (steepest descent) directivity
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                  :: iphi, idmn, so
    real(dpk)                  :: theta2, theta1
    complex(dpk)               :: n, d1, d, dmn

    theta1 = ATAN(SQRT(1._dpk - kap1**2*M3*M3)*TAN(iphi))

    if (iphi == 0.) then
       theta2 = 0.
    else if (iphi == PI) then
       theta2 = PI
    else if (TAN(iphi) < 0.) then
       theta2 = theta1 + PI
    else 
       theta2 = theta1
    end if

!!$    print*,'theta1=',theta1*180/PI,'theta=',theta*180/PI

!    so = (kap1*COS(theta2) - kap1**2*M3)/(1._dpk - kap1**2*M3*M3)  ! stationary point

    so = (kap1*COS(iphi)/SQRT(1._dpk - (kap1*M3*SIN(iphi))**2) - kap1**2*M3)/ &
         (1._dpk - kap1**2*M3*M3)  ! stationary point

!!$    print*,''
!!$    print*,'The saddle point:'
!!$    write(*,'(A12,2X,F15.10)'), 'so:->', so   
    
    n = wr/PI*(1._dpk - so*M3)**2*fplus(CMPLX(so,0._dpk,kind=dpk))

    d1 = kap1*wr*SIN(iphi)/SQRT(1._dpk - kap1**2*M3**2*SIN(iphi)*SIN(iphi))
    
    d = dhank1(d1,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*d1)*kap1*SIN(iphi)

    dmn = n/d

    idmn = 20._dpk*LOG10(ABS(dmn))
    dmn = CMPLX(idmn,0._dpk,kind=dpk)
!!$    print*,'dmn:',dmn


  END SUBROUTINE directcomp


  SUBROUTINE computeift

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the pressure/potential at each of the physical points
!! 2. Compute the constituent parts of the full pressure (potential)
!! 3. Write the data
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), allocatable, dimension(:)       :: pr
    complex(dpk), allocatable, dimension(:,:)     :: prsum
    complex(dpk), allocatable, dimension(:,:,:)   :: supinspressure
    integer                                       :: i,j,k,f
    integer                                       :: switch
    character(60)                                 :: pos, supind

!! allocate the various arrays:
    
    allocate(pr(totiftpts-1))
    allocate(prsum(Nmeshr,Nmeshz))
    allocate(pressure(Nmeshr,Nmeshz))
    allocate(acoupressure(Nmeshr,Nmeshz))
    allocate(inspressure1(Nmeshr,Nmeshz))
    allocate(inspressure2(Nmeshr,Nmeshz))
    allocate(totpressure(Nmeshr,Nmeshz))
    allocate(initpressure(Nmeshr,Nmeshz))

    if (Nzsp > 0) then
       allocate(supinspressure(Nzsp,Nmeshr,Nmeshz))
    end if

!! the basic do loop:

    do j = 1, Nmeshz  !! the axial direction
       do i = 1, Nmeshr  !! the radial (vertical in a 2D plot) direction

          print*,'Pressure at Z=',Z(j),'and R=',R(i)

!! different parts of the physical domain has separate formulae, so check:

          if (R(i) <= h) switch = 1
          if ((R(i) > h) .AND. (R(i) <= 1.)) switch = 2
          if (R(i) > 1.) switch = 3

          write(pos,"('R=',F10.5,'Z=',F10.5)"),R(i),Z(j)

          do k = 1, totiftpts-1

             call trapezoid1(R(i),Z(j),k,switch,pr(k))
             
          end do

          open(10,file='DataDump/Ift/intift.'//pos,form='FORMATTED')
          do k = 1, totiftpts-1
             write(10,'(I5,2E20.10)') k,pr(k)
          end do
          close(10)

          call evalintegral1(pr,prsum(i,j))
          
          if (prswitch == 0) then  !! compute velocity potential
             
             pressure(i,j) = wr/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk))*prsum(i,j)

!!$          pressure(i,j) = 0.
          
             acoupressure(i,j) = pressure(i,j)  !! acoustic part

             if (vortswitch .EQ. 0) then
                inspressure1(i,j) = residuepot(R(i),Z(j),1)  !! inner instability wave 
             end if

             inspressure2(i,j) = residuepot(R(i),Z(j),2)  !! outer instability wave

             if (vortswitch .EQ. 0) then
                do k = 1, Nzsp
                   supinspressure(k,i,j) = residuepot(R(i),Z(j),k+2)
                end do
             end if
             
             if (vortswitch .EQ. 0) then
                totpressure(i,j) = acoupressure(i,j) + inspressure1(i,j) + inspressure2(i,j)  !! total part 
             else
                totpressure(i,j) = acoupressure(i,j) + inspressure2(i,j)
             end if

             if (vortswitch .EQ. 0) then
                do k = 1, Nzsp
                   totpressure(i,j) = totpressure(i,j) + supinspressure(k,i,j)
                end do
             end if

             print*,'pressure:',REAL(totpressure(i,j))

          else  !! compute pressure

             if (reflswitch == 1) then

                pressure(i,j) = wr*wr/(2._dpk*PI)*prsum(i,j)   !! incident wave NOT added

             else

                pressure(i,j) = wr*wr/(2._dpk*PI)*prsum(i,j) + psi0(R(i),Z(j),switch)  !! the incident wave added

             end if
          
             acoupressure(i,j) = pressure(i,j)  !! acoustic part

             if (vortswitch .EQ. 0) then
                inspressure1(i,j) = residuepr(R(i),Z(j),1)  !! inner instability wave
             end if

             inspressure2(i,j) = residuepr(R(i),Z(j),2)  !! outer instability wave

             if (vortswitch .EQ. 0) then
                do k = 1, Nzsp
                   supinspressure(k,i,j) = residuepr(R(i),Z(j),k+2)
                end do
             end if

             if (vortswitch .EQ. 0) then
                totpressure(i,j) = acoupressure(i,j) + inspressure1(i,j) + inspressure2(i,j)  !! total part
             else
                totpressure(i,j) = acoupressure(i,j) + inspressure2(i,j)
             end if

             if (vortswitch .EQ. 0) then
                do k = 1, Nzsp
                   totpressure(i,j) = totpressure(i,j) + supinspressure(k,i,j)
                end do
             end if

             print*,'pressure:',REAL(totpressure(i,j))

          end if

          initpressure(i,j) = psi0(R(i),Z(j),switch)  !! the incident wave
          
          print*,''
          
       end do
    end do

!! dump data
    
    open(1,file='pressure.out',form='UNFORMATTED')
    write(1) Nmeshz,Nmeshr,1
    write(1) ((REAL(totpressure(i,j)),j=1,Nmeshz),i=1,Nmeshr)
    close(1)   

    open(1,file='incident.out',form='UNFORMATTED')
    write(1) Nmeshz,Nmeshr,1
    write(1) ((REAL(initpressure(i,j)),j=1,Nmeshz),i=1,Nmeshr)
    close(1)   

    open(1,file='acousticpr.out',form='UNFORMATTED')
    write(1) Nmeshz,Nmeshr,1
    write(1) ((REAL(acoupressure(i,j)),j=1,Nmeshz),i=1,Nmeshr)
    close(1)   

    if (vortswitch .EQ. 0) then
       open(1,file='instabilitypr1.out',form='UNFORMATTED')
       write(1) Nmeshz,Nmeshr,1
       write(1) ((REAL(inspressure1(i,j)),j=1,Nmeshz),i=1,Nmeshr)
       close(1)   
    end if
    
    open(1,file='instabilitypr2.out',form='UNFORMATTED')
    write(1) Nmeshz,Nmeshr,1
    write(1) ((REAL(inspressure2(i,j)),j=1,Nmeshz),i=1,Nmeshr)
    close(1)  

    if (vortswitch .EQ. 0) then
       do k = 1, Nzsp
          write(supind,"(I1)"),k
          open(1,file='supinstabilitypr.out.'//supind,form='UNFORMATTED')
          write(1) Nmeshz,Nmeshr,1
          write(1) ((REAL(supinspressure(k,i,j)),j=1,Nmeshz),i=1,Nmeshr)
          close(1)  
       end do
    end if

    open(1,file='totprdata.out',form='UNFORMATTED')
    write(1) Nmeshz,Nmeshr,totpressure
    close(1)   

    open(1,file='acousprdata.out',form='UNFORMATTED')
    write(1) Nmeshz,Nmeshr,acoupressure
    close(1)   
    
    
  END SUBROUTINE computeift


  SUBROUTINE residuepolarsup

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the near-field pressure for supersonic instability modes in polar form
!! 2. Write the data
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), allocatable, dimension(:,:,:)   :: supinspressure
    integer                                       :: i, j, k
    character(60)                                 :: pos, supind

!! allocate the various arrays:

    allocate(supinspressure(Nzsp,Nmeshr,Nphi))

    print*,''
    print*,'Also, computing polar near-field for Nzsp > 0:'
    print*,''

    do j = 1, Nphi  !! the polar sweep anle
       do i = 1, Nmeshr  !! the radial (vertical in a 2D plot) direction

          print*,'Pressure at phi=', phi(j)*180._dpk/PI,'and R=',R(i)
          
          do k = 1, Nzsp
             supinspressure(k,i,j) = residueprpolar(R(i),phi(j),k)
          end do
          
          print*,''
          
       end do
    end do

!! dump data
    
    do k = 1, Nzsp
       write(supind,"(I1)"),k
       open(1,file='supinstabilitypr_polar.out.'//supind,form='UNFORMATTED')
       write(1) Nphi,Nmeshr,1
       write(1) ((REAL(supinspressure(k,i,j)),j=1,Nphi),i=1,Nmeshr)
       close(1)  
    end do

    do k = 1, Nzsp
       write(supind,"(I1)"),k
       open(10,file='supinstabilitypr_trace.out'//supind,form='FORMATTED')
       do i = 1, Nmeshr
          if (mod(i,20) == 0) then
             do j = 1, Nphi
                if (mod(j,3) == 0) write(10,'(2F10.4,F20.10)') R(i), phi(j)*180._dpk/PI,&
                     REAL(supinspressure(k,i,j))
             end do
          end if
       end do
    close(10)
    end do
    
    
  END SUBROUTINE residuepolarsup


  SUBROUTINE precompute

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the various parameters needed while computing the IFT contour
!! 2. A check if the radial wavenumber is negative
!! 3. Constants computed: K^{-}(\mu_{mn}^{+}); K^{+}(s_{z1}); K^{+}(s_{z2})
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)          :: k, kp, f1, f2, f3
    integer               :: f, i1

    if ((vortswitch == 1) .OR. (vortswitch == 2)) then

!! the factor \Psi_{mn}(1) of (4.1) [see the JFM 2008]:

       psi = wr*resp*(1._dpk - M2*mup)*Trsin(1._dpk,mup,2)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*mup*(-Zo))

    else

!! check if the radial wavenumber is negative:

       if((REAL(mup) < -(kap1/(1._dpk - kap1*M2))) .OR. (REAL(mup) > (1._dpk/(1._dpk + M1)))) then 
          print*,'Radial wave number(s) is negative: Check the incident wave axial wave number!'
          STOP
       end if

       alpha1 = wr*SQRT((1._dpk - mup*M1)**2 - mup**2)
       alpha2 = wr*SQRT(kap1**2*(1._dpk - mup*M2)**2 - mup**2)

       print*,''
       print*,'The radial wave numbers:'
       write(*,'(/A12,2X,2F15.10)'), 'alpha1:->', alpha1
       write(*,'(A12,2X,2F15.10/)'), 'alpha2:->', alpha2

!! the factor \Psi_{mn}(1) of (3.30) [see the JFM]:

       f1 = (1._dpk - mup*M1)/(1._dpk - mup*M2)*bessj(alpha1*h,circmod,1)*EXP(ABS(AIMAG(alpha1*h)))
       
       f2 = (bessj(alpha2,circmod,1)*dhank1(alpha2,circmod,1)- & 
            hank1(alpha2,circmod,1)*dbessj(alpha2,circmod,1))* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*alpha2 + ABS(AIMAG(alpha2)))
       
       f3 = bessj(alpha2*h,circmod,1)*dhank1(alpha2,circmod,1)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*alpha2 + ABS(AIMAG(alpha2*h)))- & 
            hank1(alpha2*h,circmod,1)*dbessj(alpha2,circmod,1)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*alpha2*h + ABS(AIMAG(alpha2)))

       psi = kap2*f1*f2/f3
    
!!$    print*,'f1:',f1
!!$    print*,'f2:',f2
!!$    print*,'f3:',f3
!!$    print*,'psi:',psi

    end if
    
!!  the factor Kt^{-}(\mu_{mn}^{+}):

    if (vortswitch .EQ. 0) then
       if (REAL(mup) < crpt) then  !! muplus is below the contour; crpt being the crossover pt
          print*, ''
          print*, 'The incident acoustic mode needs to be INSIDE the contour'
          print*, ''
          STOP
       end if
    end if

!    call kernel_eval(mup,k,0,0,0)
    call kernel_eval(mup,k,0,0,1)
       
    kp = EXP(-k/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) + & 
         LOG(kernel(0,mup)/zp(0,mup)))

!    kmmup =  kernel(0,mup)/kp   !! since: Kt^{-}(s) = K^{-}(s)
    kmmup =  kernel(0,mup)/(kp*zp(0,mup)) 
!!$    print*,'kmmup:',kmmup

!!  the factor Kt^{+}(s_{z1}):

    if (vortswitch .EQ. 0) then

       call kernel_eval(sz1,k,0,0,1)

       kpsz1 = EXP(-k/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! NOTE: zero sz1 has to lie below
                                                                  !! the contour
    end if
!!  the factor K^{+}(s_{z2}):

    call kernel_eval(sz2,k,0,0,1)

    kpsz2 = EXP(-k/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! same for sz2

    if (Nzsp > 0) allocate(kzsp(Nzsp))

    if (vortswitch .EQ. 0) then
       do i1 = 1, Nzsp
          
          call kernel_eval(szsp(i1),k,0,0,1)
          
          kzsp(i1) = EXP(-k/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! NOTE: instability zeros are below
          
       end do
    end if


  END SUBROUTINE precompute


  SUBROUTINE kernel_eval(zi,integral,kswitch,ch,switch)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Evaluate (A1)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)           :: zi
    complex(dpk)           :: integral ! final value of integration
    integer                :: N, ch
    integer                :: switch  !! 1: K/U; else: K
    integer                :: kswitch  !! 1: compute (4.10); else: compute (3.22)

    allocate(intpanel(totinitpts-1))
    allocate(Npanel(totinitpts-1))

    call adaptive(zi,kswitch,ch,switch)

    call evalintegral(integral,N,zi)

    deallocate(intpanel)
    deallocate(Npanel)


  END SUBROUTINE kernel_eval


  SUBROUTINE adaptive(si,ksw,check,sw)
    
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The kernel integration is computed using an adaptive routine; which not feasible
!!    for the IFT integral
!! 2. This is integrated with the trapezoidal rule used for quadrature
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), allocatable, dimension(:)    :: zp, T_temp
    complex(dpk), allocatable, dimension(:)    :: knl, zall, zall_temp
    real(dpk), allocatable, dimension(:)       :: xp, yp
    real(dpk)                                  :: ds, len
    complex(dpk)                               :: T, si, scale
    integer                                    :: Np, Nq, Nr, panel_no, Nloopindex, check, sw, ksw
    integer                                    :: i, j, ii


    do j = 1, totinitpts-1  !! the integration points already specified

       len = initpoints(j+1) - initpoints(j)  !! the length of a mini panel

!! compute "scale" needed to ensure the tolerance check:

       if (j >= 1 .AND. j < N1+2 ) then
          scale = len/(szk(3) - szk(1))

       else
          scale = len/(szk(5) - szk(3))

       end if

       call trapezoid(si,initpoints(j),initpoints(j+1),ksw,sw,T)  !! the basis for comparison

!       print*, 'T:', T

       Np = 1

!! the adaptive loop starts here:

       do
          
          panel_no = 2**Np  !! min no of panels = 2
          ds = len/panel_no

!          print*, 'panels:', panel_no
          
          allocate(xp(panel_no+1))
          allocate(yp(panel_no+1))
          allocate(zp(panel_no+1))

          xp(1) = REAL(initpoints(j))
          xp(panel_no+1) = REAL(initpoints(j+1))
          yp(1) = AIMAG(initpoints(j))
          yp(panel_no+1) = AIMAG(initpoints(j+1))
          zp(1) = initpoints(j)
          zp(panel_no+1) = initpoints(j+1)

          do i = 1,panel_no-1
      
             if (j >= 1 .AND. j < N1+2) then
                xp(i+1) = xp(i) + ds
                call yfunc_p(xp(i+1),yp(i+1),REAL(szk(3)),REAL(szk(2)-szk(3)),AIMAG(szk(2)-szk(3)))
             else
                xp(i+1) = xp(i) + ds
                call yfunc_p(xp(i+1),yp(i+1),REAL(szk(3)),REAL(szk(4)-szk(3)),AIMAG(szk(4)-szk(3)))
                  
             end if

             zp(i+1) = CMPLX(xp(i+1),yp(i+1))
             
          end do

!          print*, 'zp(i):', zp(i), 'zp(i+1):', zp(i+1)
          
          allocate(T_temp(panel_no))
          
          do i = 1,panel_no
             call trapezoid(si,zp(i),zp(i+1),ksw,sw,T_temp(i))
          end do
          
          intpanel(j) = (0._dpk,0._dpk)

          do i = 1,panel_no
             intpanel(j) = intpanel(j) + T_temp(i)
          end do

          deallocate(xp)
          deallocate(yp)
          deallocate(T_temp)

!          print*, 'initpoints:', initpoints(j)
!          print*, 'intpanel:', intpanel(j)

!! the tolerance check:

          if (1._dpk/(4.**Np-1._dpk)*ABS(intpanel(j)-T) <= ABS(scale*tol)) then
             Npanel(j) = panel_no
!            write(*,'(/A15,I10,5X,A15,I10)'), 'Panel no: ', j, 'Adaptive div: ', Npanel(j)

!! the kernel check (see Pg 436, JFM, Vol. 612, 2008):

             if (check == 1) then
                   if (j == 1) then
                      allocate(zall(Npanel(j)+1))
                      do i = 1, Npanel(j)+1
                         zall(i) = zp(i)
                      end do
                   else 
                      Nq = 1
                      do i = 1, j-1
                         Nq = Nq + Npanel(i)
                      end do
                      Nr = Nq + Npanel(j)
                      allocate(zall_temp(Nr))
                      do ii = 1, Nq
                         zall_temp(ii) = zall(ii)
                      end do
                      do i = Nq+1, Nr
                         zall_temp(i) = zp(i+1-Nq)
                      end do
                      deallocate(zall)
                      allocate(zall(Nr))
                      zall = zall_temp
                      deallocate(zall_temp)
                   end if
             end if

             deallocate(zp)

             EXIT  !! exit to the next panel as soon the tolerance is satisfied
          end if

          deallocate(zp)
          Np = Np+1

       end do

    end do

    if (check == 1) then
       allocate(knl(Nr))
       call kernelcheck(Nr,ksw,zall,knl)
       call kernelplot(Nr,zall,knl)
    end if

  END SUBROUTINE adaptive


  SUBROUTINE computefplus(switch,index)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Evaluate (3.30)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer          :: i, index, switch

    allocate(fplusz(totiftpts))  !! fplus is same as \xi^{+}(s) of (3.30)
    fplusz = (0._dpk,0._dpk)

    if (switch == 0) then  !! switch = 0 :> fresh job; no restart file read before

       do i = 1,totiftpts
          print*, 'F+ at:', iftpoints(i)
          fplusz(i) = fplus(iftpoints(i))
          
          write(*,'(A22,2X,2F20.10/)'),'F+:->',fplusz(i)

          if (i == 1) then
             open(10,file='fplus_part.out',form='FORMATTED',status='UNKNOWN')
          else
             open(10,file='fplus_part.out',form='FORMATTED',position='APPEND')
          end if

          write(10,'(I5,4E20.10)') i,iftpoints(i),fplusz(i)
          close(10)

       end do

    else

       fplusz = fplusz_temp

       if (index .NE. totiftpts) then

          do i = index+1,totiftpts
             print*, 'F+ at:', iftpoints(i)
             fplusz(i) = fplus(iftpoints(i))
          
             write(*,'(A22,2X,2F20.10/)'),'F+:->',fplusz(i)

             open(10,file='fplus_part.out',form='FORMATTED',position='APPEND')
             write(10,'(I5,4E20.10)') i,iftpoints(i),fplusz(i)
             close(10)

          end do

       end if

    end if

    open(10,file='fplus.out',form='FORMATTED')
    do i = 1, totiftpts
       write(10,'(I5,4E20.10)') i,iftpoints(i),fplusz(i)
    end do
    close(10)
    

  END SUBROUTINE computefplus


  SUBROUTINE readfplus(j)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Use the computed F^{+} in case running a restart job
!! 2. Note that restarting is only possible during the "computing F+" portion
!! 3. If a job exits during (R,Z) computation, a "restart" just starts the job back at (R=0,Z=0)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer          :: i, j
    real             :: dummy1, dummy2

    allocate(fplusz_temp(totiftpts))

    fplusz_temp = (0._dpk,0._dpk)

    open(10,file='fplus_part.out',form='FORMATTED')
    do i = 1, totiftpts
       read(10,'(I5,4E20.10)',end = 100) j,dummy1,dummy2,fplusz_temp(i)
    end do
    close(10)

100 do i = 1,j
       print*, 'F+ at:', iftpoints(i)

       write(*,'(A22,2X,2F20.10/)'),'F+:->',fplusz_temp(i)
    end do

    close(10)

  END SUBROUTINE readfplus


  FUNCTION fplus(z)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Actual computation of \xi^{+}(s)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)       :: fplus, z, gpz, kpz, k
    complex(dpk)       :: l1mmup, l2mmup, l3mmup, l1pz, l2pz, l3pz, lprod
    integer            :: f


    l1mmup = sqrt(1._dpk - mup*(M1 - 1._dpk))
    l2mmup = sqrt(kap1 - mup*(kap1*M2 - 1._dpk))
    l3mmup = sqrt(kap1 - mup*(kap1*M3 - 1._dpk))
    l1pz = sqrt(1._dpk - z*(M1 + 1._dpk))
    l2pz = sqrt(kap1 - z*(kap1*M2 + 1._dpk))
    l3pz = sqrt(kap1 - z*(kap1*M3 + 1._dpk))

!!$    gpz = psi*(1._dpk - mup*M2)/(kmmup*(mup - z))
    gpz = psi*(1._dpk - mup*M2)/(kmmup*(mup - sz2))

!!$    gpz = psi*(1._dpk - mup*M2)/(kmmup*(mup - z)) - &
!!$         psi*(1._dpk - mup*M2)/(kmmup*(mup - sz2))

!!$    lprod = l2mmup*l3mmup*l2pz*l3pz

!    call kernel_eval(z,k,0,0,0)
    call kernel_eval(z,k,0,0,1)

    if ((farswitch == 1) .OR. (farswitch == 2)) then
       if (REAL(z) >= crpt) then
          kpz = EXP(-k/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) + & 
               LOG(kernel(0,z)/zp(0,z)))  !! z is above contour; crpt is crossover pt
       else
          kpz = EXP(-k/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! z is below
       end if
    else
       kpz = EXP(-k/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! z is always below normally
    end if

!!$    fplus = gpz/kpz*( (z - sz2)/(mup - z) + vparm)
    fplus = gpz*( (z - sz2)/(mup - z) + vparm)/(kpz*zp(0,z))

!!$    fplus = lprod*gpz/(kmmup*kpz)

!!$    fplus = lprod*gpz*(z - sp1)/(kmmup*(z - sz1)*(z - sz2)*kpz)

!!$    print*,'gpz',gpz


  END FUNCTION fplus


  FUNCTION zp(ss,z)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The factor U(s)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)  :: z, zp
    integer       :: j, jj, ss

    if (ss == 1) then
!       zp = ((z-sp1)*(z-CONJG(sp1)))
       zp = (z-sp1)
    else
       if (vortswitch .EQ. 0) then
!!$       zp = (z-sz1)*(z-CONJG(sz1))*(z-sz2)*(z-CONJG(sz2)) &
!!$            /((z-sp1)*(z-CONJG(sp1)))
          zp = (z-sz1)*(z-sz2)/(z-sp1)
       else
          zp = (z-sz2)
       end if
    end if

!!$    zp = (z-sz1)*(z-sz2)/(z-sp1)

    if (ss == 1) then
       do j = 1, Npsp
          zp = zp*(z - spsp(j))
       end do
    else
       if (vortswitch .EQ. 0) then
          do j = 1, Nzsp
             do jj = 1, Npsp
                zp = zp*(z - szsp(j))/(z - spsp(jj))
             end do
          end do
       end if
    end if
!!$    zp = (z-sz2)

!!$    zp = CMPLX(1._dpk,0._dpk,kind=dpk)

  END FUNCTION zp


  SUBROUTINE yfunc_p2(x,y,p,a,e1,e2)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Parabola used for the last loop (down) of integration contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk) :: p
    real(dpk)    :: x, y
    real(dpk)    :: a, e1, e2

    y = -(x-REAL(p))*(1._dpk-(x-REAL(p))*(a-e1-e2)/a**2) + e1


  END SUBROUTINE yfunc_p2


  SUBROUTINE yfunc_p1(x,y,p,e1,e2)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Parabola used for the other loops (up) of integration contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk) :: p
    real(dpk)    :: x, y
    real(dpk)    :: e1, e2

    y = (x-REAL(p))*(1._dpk-(x-REAL(p))/e1) + e2


  END SUBROUTINE yfunc_p1


  SUBROUTINE yfunc_p0(x,y,p,e1,e2)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Parabola used for the other loops (down) of integration contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk) :: p
    real(dpk)    :: x, y
    real(dpk)    :: e1, e2

    y = -(x-REAL(p))*(1._dpk-(x-REAL(p))/e1) + e2


  END SUBROUTINE yfunc_p0


  SUBROUTINE yfunc_p(x,y,p,a,b)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Algebraic "one piece" contour:  see Rienstra, JEM, 2007
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)    :: x, y
    real(dpk)    :: a, b, p

    y = b*(4._dpk*(x-p)/a)/(3._dpk + ((x-p)/a)**4) ! algebraic


  END SUBROUTINE yfunc_p


  SUBROUTINE trapezoid(s,z1,z2,ksw,sw,I)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The trapezoidal rule for the kernel contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)                :: z1, z2, s, dz, I
    complex(dpk)                :: f1, f2
    integer                     :: sw, ksw

    f1 = integrand(ksw,sw,s,z1)
    f2 = integrand(ksw,sw,s,z2)

    dz = (z2-z1)

    I = dz/2._dpk*(f1+f2)


  END SUBROUTINE trapezoid


  SUBROUTINE trapezoid1(r,z,i,ss,Int)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The trapezoidal rule for the IFT contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                   :: r, z
    complex(dpk)                :: ds, Int
    complex(dpk)                :: f1, f2
    integer                     :: i, ss

    if (prswitch == 0) then

       f1 = integrandiftpot(r,z,i,ss)
       f2 = integrandiftpot(r,z,i+1,ss)

    else

       f1 = integrandiftpr(r,z,i,ss)
       f2 = integrandiftpr(r,z,i+1,ss)

    end if

    ds = iftpoints(i+1) - iftpoints(i)

    Int = ds/2._dpk*(f1+f2)


  END SUBROUTINE trapezoid1


  SUBROUTINE evalintegral(Int,totpoints,z)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Sum up the individual trapezoids for the kernel contour
!! 2. Dump data
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)                               :: Int, z
    integer                                    :: totpoints
    integer                                    :: i
    character(60)                              :: arg

    Int = (0._dpk,0._dpk)
    do i = 1,totinitpts-1
       Int = Int + intpanel(i)  !! summing up the individual panels
    end do

!! dump data:

    write(arg,"('u=',E12.5,'+',E12.5,'i')"),REAL(z),AIMAG(z)

    open(10,file='DataDump/Kernel/intkernel.'//arg,form='FORMATTED')
    do i = 1, totinitpts-1
       write(10,'(I5,2E20.10)') i,intpanel(i)
    end do
    close(10)
    
    write(*,'(/A22,2X,2F20.10/)'),'Integral value:->', Int

    totpoints = 1
    do i = 1,totinitpts-1
       totpoints = totpoints + Npanel(i)  !! total quadrature points
    end do

    write(*,'(A22,2X,I20/)'),'Quadrature pts:->', totpoints
    

  END SUBROUTINE evalintegral


 SUBROUTINE evalintegral1(panel,Int)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Sum up the individual trapezoids for the IFT contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)                               :: Int
    complex(dpk),dimension(totiftpts-1)        :: panel
    integer                                    :: i

    Int = (0._dpk,0._dpk)
    do i = 1,totiftpts-1
       Int = Int + panel(i)
    end do


  END SUBROUTINE evalintegral1


  FUNCTION integrand(kswitch,switch,si,zi)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the kernel integrand
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: si, zi, integrand
    complex(dpk)    :: In, Id
    integer         :: switch, kswitch

    if (switch == 1)  then !! K/U
       In = LOG(kernel(kswitch,zi)/zp(kswitch,zi))
    else
       In = LOG(kernel(kswitch,zi))
    end if
    Id = zi-si
    integrand = In/Id
    
!!$    if (ABS(integrand) > 1.0E3) then
!!$       print*,'zi:',zi
!!$!       print*,'integrand:',integrand
!!$!       print*,'kernel:',kernel(zi)
!!$       print*,'In:',In
!!$       print*,'Id:',Id
!!$    end if

    if (kernel(kswitch,zi) == 0.) then 
       print*,''
       print*,'FATAL ERROR: The kernel has encountered a zero (at zi)!!'
       print*,'si:-->',si
       print*,'zi:-->',zi
       STOP
    end if


  END FUNCTION integrand


  FUNCTION integrandiftpr(ri,zi,i,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the IFT integrand when computing for pressure
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)       :: ri, zi
    integer         :: i, ss
    complex(dpk)    :: u
    complex(dpk)    :: integrandiftpr


    u = iftpoints(i)

    if (ss==1) then

!! pressure:

       integrandiftpr = (1._dpk - u*M2)**2*Trs(ri,u,1)*fplusz(i)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*u*zi)
!!$       integrandiftpr = Trs(ri,u,1)
    else if (ss==2) then

!! pressure:
       
       integrandiftpr = (1._dpk - u*M2)**2*Trs(ri,u,2)*fplusz(i)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*u*zi)
!!$       integrandiftpr = Trs(ri,u,2)
    else

!! pressure:

       integrandiftpr = (1._dpk - u*M3)**2*Trs(ri,u,3)*fplusz(i)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*u*zi)
!!$       integrandiftpr =  Trs(ri,u,3)
    end if


  END FUNCTION integrandiftpr


  FUNCTION integrandiftpot(ri,zi,i,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the IFT integrand when computing for potential
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)       :: ri, zi
    integer         :: i, ss
    complex(dpk)    :: u
    complex(dpk)    :: integrandiftpot


    u = iftpoints(i)

    if (ss==1) then

!! potential:

       integrandiftpot = (1._dpk - u*M2)**2/(1._dpk - u*M1)*Trs(ri,u,1)*fplusz(i)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*u*zi)

    else if (ss==2) then

!! potential:
       
       integrandiftpot = (1._dpk - u*M2)*Trs(ri,u,2)*fplusz(i)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*u*zi)

    else

!! potential:

       integrandiftpot = (1._dpk - u*M3)*Trs(ri,u,3)*fplusz(i)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*u*zi)

    end if


  END FUNCTION integrandiftpot


  FUNCTION kernel(ss,z)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the kernel function (3.22) or (4.10)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: kernel, z
    complex(dpk)    :: l1, l2, l3, lz
    complex(dpk)    :: F1n, F1d, F1, F1f, F2n, F2d, F2, F2f, Rn, Rd, Rz
    integer         :: ss

    if (ss == 1) then  !! compute (4.10): the inc vort upstream kernel

       l1 = sqrt(1._dpk - z*(M1+1._dpk))*sqrt(1._dpk - z*(M1-1._dpk))
       l2 = sqrt(kap1 - z*(kap1*M2+1._dpk))*sqrt(kap1 - z*(kap1*M2-1._dpk))


       if ((ABS(l1*wr*h) < asymplim .AND. ABS(AIMAG(l1*wr*h)) < asymplim1)) then

          F1n = bessj(l1*wr*h,circmod,1)
          F1d = dbessj(l1*wr*h,circmod,1)
          F1f = F1n/F1d

       else

          F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
          
       end if

!!$    print*,'F1n:',F1n
!!$    print*,'F1d:',F1d
!!$    print*,'F1f:',F1f
!!$    print*,'lz:',lz

       F1 = kap2*(1._dpk - z*M1)**2/l1*F1f

       if ((ABS(l2*wr) < asymplim .AND. ABS(AIMAG(l2*wr)) < asymplim1) .AND. &
            (ABS(l2*wr*h) < asymplim .AND. ABS(AIMAG(l2*wr*h)) < asymplim1)) then

          F2n = dhank1(l2*wr,circmod,1)*bessj(l2*wr*h,circmod,1)*EXP(ABS(AIMAG(l2*wr*h))+ &
               CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr) - dbessj(l2*wr,circmod,1)*hank1(l2*wr*h,circmod,1)* & 
               EXP(ABS(AIMAG(l2*wr))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*h)
          F2d = dhank1(l2*wr,circmod,1)*dbessj(l2*wr*h,circmod,1)*EXP(ABS(AIMAG(l2*wr*h))+ &
               CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr) - dbessj(l2*wr,circmod,1)*dhank1(l2*wr*h,circmod,1)* &
               EXP(ABS(AIMAG(l2*wr))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*h)
          F2f = F2n/F2d

       else
          
          F2f = CMPLX(0._dpk,1._dpk,kind=dpk)

       end if

       F2 = (1._dpk - z*M2)**2/l2*F2f

       kernel = wr*(F1 - F2)

    else  !! compute (3.22): the kernel

       l1 = sqrt(1._dpk - z*(M1+1._dpk))*sqrt(1._dpk - z*(M1-1._dpk))
       l2 = sqrt(kap1 - z*(kap1*M2+1._dpk))*sqrt(kap1 - z*(kap1*M2-1._dpk))
       l3 = sqrt(kap1 - z*(kap1*M3+1._dpk))*sqrt(kap1 - z*(kap1*M3-1._dpk))

       lz = kap2*l2/l1*(1._dpk - z*M1)**2/(1._dpk - z*M2)**2


       if ((ABS(l2*wr) < asymplim .AND. ABS(AIMAG(l2*wr)) < asymplim1) .AND. &
            (ABS(l2*wr*h) < asymplim .AND. ABS(AIMAG(l2*wr*h)) < asymplim1) .AND. &
            (ABS(l1*wr*h) < asymplim .AND. ABS(AIMAG(l1*wr*h)) < asymplim1)) then

          Rd = (lz*bessj(l1*wr*h,circmod,1)*dhank1(l2*wr*h,circmod,1)- &
               hank1(l2*wr*h,circmod,1)*dbessj(l1*wr*h,circmod,1))* &
               EXP(ABS(AIMAG(l1*wr*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*h)

          Rn = (bessj(l2*wr*h,circmod,1)*dbessj(l1*wr*h,circmod,1)- &
               lz*bessj(l1*wr*h,circmod,1)*dbessj(l2*wr*h,circmod,1))* &
               EXP(ABS(AIMAG(l1*wr*h))+ABS(AIMAG(l2*wr*h)))

          Rz = Rn/Rd


          F1n = Rz*hank1(l2*wr,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr)+ &
               bessj(l2*wr,circmod,1)*EXP(ABS(AIMAG(l2*wr)))
          F1d = Rz*dhank1(l2*wr,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr)+ &
               dbessj(l2*wr,circmod,1)*EXP(ABS(AIMAG(l2*wr)))
          F1f = F1n/F1d

       else

          F1f = CMPLX(0._dpk,1._dpk,kind=dpk)
       
       end if

!!$    print*,'F1n:',F1n
!!$    print*,'F1d:',F1d
!!$    print*,'F1f:',F1f
!!$    print*,'lz:',lz

       F1 = (1._dpk - z*M2)**2/l2*F1f

       if (ABS(l3*wr) < asymplim .AND. ABS(AIMAG(l3*wr)) < asymplim1) then
          
          F2n = hank1(l3*wr,circmod,1)
          F2d = dhank1(l3*wr,circmod,1)
          F2f = F2n/F2d
          
       else

          F2f = (8._dpk*l3*wr + 4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*circmod*circmod - &
               CMPLX(0._dpk,1._dpk,kind=dpk))/(8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*l3*wr - &
               4._dpk*circmod*circmod - 3._dpk)

       end if

!!$    print*,'F2n:',F2n
!!$    print*,'F2d:',F2d
!!$    print*,'F2f:',F2f

       F2 = (1._dpk - z*M3)**2/l3*F2f

       if (Nzbc == 0 .AND. Npbc == 0) then
          
          kernel = wr*(F1 - F2)
 
       else

          kernel = wr*(F1 - F2)*bc(z)

       end if

    end if


  END FUNCTION kernel


  SUBROUTINE kernelcheck(N,ss,z,I)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Check so that the kernel does not cross the negative real axis on C
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), dimension(N) :: z, I
    real(dpk)                  :: x1, x2, y1, y2, intercept
    integer                    :: j, N, ss

    do j = 1, N
       I(j) = kernel(ss,z(j))/zp(ss,z(j))
    end do

    do j = 1, N-1

       x1 = REAL(I(j))
       x2 = REAL(I(j+1))
       y1 = AIMAG(I(j))
       y2 = AIMAG(I(j+1))

       intercept = (x1*y2 - x2*y1)/(y2 - y1)

!! check the validity of the split integral:

       if((y1 < 0. .AND. y2 > 0.) .OR. (y1 > 0. .AND. y2 < 0.)) then

          ! if y1, y2 are of the same sign real axis is not crossed

          if (intercept < 0.) then 

             ! kernel function crosses the branch cut of the log function

             print*,''
             print*,'FATAL ERROR: The kernel split functions are invalid!!'
             print*,'FATAL ERROR: The kernel has crossed the negative Re axis between:'
             print*,'C(z1):z1-->',z(j)
             print*,'and'
             print*,'C(z2):z2-->',z(j+1)
             print*,'where:'
             print*,'K(z1)-->', I(j)
             print*,'K(z2)-->', I(j+1)
!!$             STOP

          end if

       end if

    end do

  END SUBROUTINE kernelcheck


  SUBROUTINE kernelplot(N,z,K)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Plot the kernel if requested
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer                        :: N, i
    complex(dpk), dimension(N)     :: z, K
    character(40)                  :: arg

    write(arg,"('w=',E12.5,'Np=',I7.7)"),REAL(wr), N

    open(10,file='DataDump/Kernel_Trace/kernel.'//arg,form='FORMATTED')
    do i = 1, N
       write(10,'(I8,2F20.10,2F30.10)') i,z(i),K(i)
    end do
    close(10)


  END SUBROUTINE kernelplot


  FUNCTION bc(z)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the term to be factored out from the kernel function (3.22)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)  :: z, bc
    integer       :: i

    bc = 1._dpk

    do i = 1, Nzbc

       bc = bc/(z - szbc(i))

    end do

    do i = 1, Npbc

       bc = bc*(z - spbc(i))

    end do


  END FUNCTION bc


  FUNCTION Trs(ri,si,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute Trs of (3.33)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: ri
    complex(dpk)  :: sri, szi, Trs
    complex(dpk)  :: si
    complex(dpk)  :: l1, l2, l3, lz
    complex(dpk)  :: Fn, Fd, F, Rn, Rd, Rz
    integer       :: ss


    l1 = sqrt(1._dpk - si*(M1+1._dpk))*sqrt(1._dpk - si*(M1-1._dpk))
    l2 = sqrt(kap1 - si*(kap1*M2+1._dpk))*sqrt(kap1 - si*(kap1*M2-1._dpk))
    l3 = sqrt(kap1 - si*(kap1*M3+1._dpk))*sqrt(kap1 - si*(kap1*M3-1._dpk))

    lz = kap2*l2/l1*(1._dpk - si*M1)**2/(1._dpk - si*M2)**2


    if (ss==1) then
       
       Rd = (lz*bessj(l1*wr*h,circmod,1)*dhank1(l2*wr*h,circmod,1)- &
            hank1(l2*wr*h,circmod,1)*dbessj(l1*wr*h,circmod,1))* &
            EXP(ABS(AIMAG(l1*wr*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*h)
       
       Rn = (bessj(l2*wr*h,circmod,1)*dbessj(l1*wr*h,circmod,1)- &
            lz*bessj(l1*wr*h,circmod,1)*dbessj(l2*wr*h,circmod,1))* &
            EXP(ABS(AIMAG(l1*wr*h))+ABS(AIMAG(l2*wr*h)))
       
       Rz = Rn/Rd
          
       Fn = bessj(l1*wr*ri,circmod,1)*EXP(ABS(AIMAG(l1*wr*ri)))*(Rz*hank1(l2*wr*h,circmod,1)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*h) + bessj(l2*wr*h,circmod,1)* &
            EXP(ABS(AIMAG(l2*wr*h))))
       
       Fd = bessj(l1*wr*h,circmod,1)*EXP(ABS(AIMAG(l1*wr*h)))*(Rz*dhank1(l2*wr,circmod,1)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr) + dbessj(l2*wr,circmod,1)* &
            EXP(ABS(AIMAG(l2*wr))))
       
       F = Fn/Fd

       Trs = F/l2

    else if (ss==2) then

       Rd = (lz*bessj(l1*wr*h,circmod,1)*dhank1(l2*wr*h,circmod,1)- &
            hank1(l2*wr*h,circmod,1)*dbessj(l1*wr*h,circmod,1))* &
            EXP(ABS(AIMAG(l1*wr*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*h)
       
       Rn = (bessj(l2*wr*h,circmod,1)*dbessj(l1*wr*h,circmod,1)- &
            lz*bessj(l1*wr*h,circmod,1)*dbessj(l2*wr*h,circmod,1))* &
            EXP(ABS(AIMAG(l1*wr*h))+ABS(AIMAG(l2*wr*h)))
       
       Rz = Rn/Rd
          
       Fn = Rz*hank1(l2*wr*ri,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*ri)+ &
            bessj(l2*wr*ri,circmod,1)*EXP(ABS(AIMAG(l2*wr*ri)))
       
       Fd = Rz*dhank1(l2*wr,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr)+ &
            dbessj(l2*wr,circmod,1)*EXP(ABS(AIMAG(l2*wr)))
       
       F = Fn/Fd

       Trs = F/l2

    else

!!$       if ((ABS(l3*wr) < asymplim .AND. ABS(AIMAG(l3*wr)) < asymplim1) .AND. &
!!$            (ABS(l3*wr*ri) < asymplim .AND. ABS(AIMAG(l3*wr*ri)) < asymplim1)) then

          Fn = hank1(l3*wr*ri,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l3*wr*ri)
       
          Fd = dhank1(l3*wr,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l3*wr)
       
          F = Fn/Fd

!!$       else 
!!$       
!!$          F = (8._dpk*l3*wr*ri + 4._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*circmod*circmod - &
!!$               CMPLX(0._dpk,1._dpk,kind=dpk))/((8._dpk*CMPLX(0._dpk,1._dpk,kind=dpk)*l3*wr - &
!!$               4._dpk*circmod*circmod - 3._dpk)*ri**(1.5_dpk))*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
!!$               l3*wr*(ri - 1._dpk))
!!$
!!$       end if
  
       Trs = F/l3

    end if
       
!    print*,'Rz:',Rz


  END FUNCTION Trs


  FUNCTION Trsin(ri,si,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute Trs of (4.8)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: ri
    complex(dpk)  :: sri, szi, Trsin
    complex(dpk)  :: si
    complex(dpk)  :: l1, l2
    complex(dpk)  :: Fn, Fd, F
    integer       :: ss


    l1 = sqrt(1._dpk - si*(M1+1._dpk))*sqrt(1._dpk - si*(M1-1._dpk))
    l2 = sqrt(1._dpk - si*(M2+1._dpk))*sqrt(1._dpk - si*(M2-1._dpk))


    if (ss==1) then
       
       Fn = bessj(l1*wr*ri,circmod,1)*EXP(ABS(AIMAG(l1*wr*ri)))
       Fd = dbessj(l1*wr*h,circmod,1)*EXP(ABS(AIMAG(l1*wr*h)))
       
       F = Fn/Fd

       Trsin = F/l1

    else

       Fn = dhank1(l2*wr,circmod,1)*bessj(l2*wr*ri,circmod,1)*EXP(ABS(AIMAG(l2*wr*ri))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr) - dbessj(l2*wr,circmod,1)*hank1(l2*wr*ri,circmod,1)* &
            EXP(ABS(AIMAG(l2*wr))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*ri)

       Fd = dhank1(l2*wr,circmod,1)*dbessj(l2*wr*h,circmod,1)*EXP(ABS(AIMAG(l2*wr*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr) - dbessj(l2*wr,circmod,1)*dhank1(l2*wr*h,circmod,1)* &
            EXP(ABS(AIMAG(l2*wr))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*wr*h)
       
       F = Fn/Fd

       Trsin = F/l2

    end if
       
!    print*,'Rz:',Rz


  END FUNCTION Trsin


  FUNCTION psi0(r,z,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the \Psi_{mn}(r) of (2.14) or (4.11)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)        :: psi0, psimn
    complex(dpk)        :: f1, f2, f3
    real(dpk)           :: r, z
    integer             :: ss

    if ((vortswitch == 1) .OR. (vortswitch == 2)) then

       if (z >= Zo) then 

          if (ss == 1) then

             psimn = wr*resp*(1._dpk - M1*mup)*Trsin(r,mup,1)
          
             psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*wr*(1._dpk - M1*mup)* &
                  psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*mup*(z-Zo))

          else if (ss == 2) then
       
             psimn = wr*resp*(1._dpk - M2*mup)*Trsin(r,mup,2)
          
             psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*wr*(1._dpk - M2*mup)* &
                  psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*mup*(z-Zo))

          else

             psi0 = 0.

          end if

       else

          psi0 = 0.

       end if

    else


       if (ss == 1) then

          psimn = bessj(alpha1*r,circmod,1)*EXP(ABS(AIMAG(alpha1*r)))
       
          psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*wr*(1._dpk - M1*mup)* &
               psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*mup*z)


       else if (ss == 2) then

          f1 = (1._dpk - mup*M1)/(1._dpk - mup*M2)*bessj(alpha1*h,circmod,1)
       
          f2 = bessj(alpha2*r,circmod,1)* &
               dhank1(alpha2,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* & 
               alpha2 + ABS(AIMAG(alpha2*r)))- & 
               hank1(alpha2*r,circmod,1)* &
               dbessj(alpha2,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* & 
               alpha2*r + ABS(AIMAG(alpha2)))
          

          f3 = bessj(alpha2*h,circmod,1)* &
               dhank1(alpha2,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
               alpha2 + ABS(AIMAG(alpha2*h)))- & 
               hank1(alpha2*h,circmod,1)* &
               dbessj(alpha2,circmod,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
               alpha2*h + ABS(AIMAG(alpha2)))
            

          psimn = kap2*f1*f2/f3

          psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*wr*(1._dpk - M2*mup)* &
               psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*wr*mup*z)

!!$       open(10,file='test.out',form='FORMATTED',position='APPEND')
!!$       write(10,'(2E20.10)') psimn
!!$       close(10)

       else

          psi0 = 0.

       end if

    end if


  END FUNCTION psi0


  FUNCTION residueprpolar(r,phi,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the residue pressure term of (3.38)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)            :: r, phi
    integer              :: ss, ii, jj
    complex(dpk)         :: residueprpolar, res, fn
   
    res = psi*(1._dpk - mup*M2)/((mup - szsp(ss))* &
         kmmup*kzsp(ss)*(szsp(ss) - sz1)*(szsp(ss) - sz2))
!!$    res = psi*(1._dpk - mup*M2)*(szsp(ss) - sp1)/((mup - szsp(ss))* &
!!$         kmmup*kzsp(ss)*(szsp(ss) - sz1)*(szsp(ss) - sz2))
!!$    res = psi*(1._dpk - mup*M2)*(szsp(ss)-sp1)*(szsp(ss)-CONJG(sp1))/ &
!!$            ((mup - szsp(ss))*kmmup*kzsp(ss)*(szsp(ss)-CONJG(sz1))* &
!!$            (szsp(ss)-sz1)*(szsp(ss) - sz2)*(szsp(ss)-CONJG(sz2)))
    do ii = 1, Nzsp
       if (ii .NE. (ss)) then
          res = res/(szsp(ss) - szsp(ii))
       end if
    end do
    do jj = 1, Npsp
       res = res*(szsp(ss) - spsp(jj))
    end do
    
    if (r*SIN(phi) <= h) then
       fn = (1._dpk - szsp(ss)*M2)**2*Trs(r*SIN(phi),szsp(ss),1)
    else if (r*SIN(phi) <= 1. .AND. r*SIN(phi) > h) then
       fn = (1._dpk - szsp(ss)*M2)**2*Trs(r*SIN(phi),szsp(ss),2)
    else
       fn = (1._dpk - szsp(ss)*M3)**2*Trs(r*SIN(phi),szsp(ss),3)  
    end if

    residueprpolar = CMPLX(0._dpk,1._dpk,kind=dpk)*wr*wr*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
         wr*szsp(ss)*r*COS(phi))*res*fn

    
  END FUNCTION residueprpolar
    

  FUNCTION residuepr(r,z,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the residue pressure term of (3.38)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)            :: r, z
    integer              :: ss, ii, jj
    complex(dpk)         :: residuepr, res, fn
   

    if (ss == 1) then

!!$       if (z > 0.) then
!!$       res = psi*(1._dpk - mup*M2)*(sz1-sp1)*(sz1-CONJG(sp1))/ &
!!$            ((mup - sz1)*kmmup*kpsz1*(sz1-CONJG(sz1))*(sz1-sz2)*(sz1-CONJG(sz2)))
!!$       res = psi*(1._dpk - mup*M2)*(sz1 - sp1)/((mup - sz1)*kmmup*kpsz1*(sz1 - sz2))
       res = psi*(1._dpk - mup*M2)/((mup - sz1)*kmmup*kpsz1*(sz1 - sz2))
       if (vortswitch .EQ. 0) then
          do ii = 1, Nzsp
             res = res/(sz1 - szsp(ii))
          end do
          do jj = 1, Npsp
             res = res*(sz1 - spsp(jj))
          end do
       end if

       if (r <= h) then
          fn = (1._dpk - sz1*M2)**2*Trs(r,sz1,1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - sz1*M2)**2*Trs(r,sz1,2)
       else
          fn = (1._dpk - sz1*M3)**2*Trs(r,sz1,3)
       end if
       
       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*wr*wr*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            wr*sz1*z)*res*fn
          
!!$       else
!!$          residuepr = 0.
!!$
!!$       end if
          
    elseif (ss == 2) then

!!$       if (z > 0.) then
!!$       res = psi*(1._dpk - mup*M2)*(sz2-sp1)*(sz2-CONJG(sp1))/ &
!!$            ((mup - sz2)*kmmup*kpsz2*(sz2-CONJG(sz1))*(sz2-sz1)*(sz2-CONJG(sz2)))
!!$       res = psi*(1._dpk - mup*M2)*(sz2 - sp1)/((mup - sz2)*kmmup*kpsz2*(sz2 - sz1))
       res = psi*(1._dpk - mup*M2)/((mup - sz2)*kmmup*kpsz2*(sz2 - sz1))
       if (vortswitch .EQ. 0) then
          do ii = 1, Nzsp
             res = res/(sz2 - szsp(ii))
          end do
          do jj = 1, Npsp
             res = res*(sz2 - spsp(jj))
          end do
       end if
          
       if (r <= h) then
          fn = (1._dpk - sz2*M2)**2*Trs(r,sz2,1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - sz2*M2)**2*Trs(r,sz2,2)
       else
          fn = (1._dpk - sz2*M3)**2*Trs(r,sz2,3)
       end if

       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*wr*wr*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            wr*sz2*z)*res*fn

!!$       else
!!$          residuepr = 0.
!!$
!!$       end if

    else

!!$       if (z > 0.) then
!!$       res = psi*(1._dpk - mup*M2)*(szsp(ss-2)-sp1)*(szsp(ss-2)-CONJG(sp1))/ &
!!$            ((mup - szsp(ss-2))*kmmup*kzsp(ss-2)*(szsp(ss-2)-CONJG(sz1))* &
!!$            (szsp(ss-2)-sz1)*(szsp(ss-2) - sz2)*(szsp(ss-2)-CONJG(sz2)))
!!$       res = psi*(1._dpk - mup*M2)*(szsp(ss-2) - sp1)/((mup - szsp(ss-2))* &
!!$            kmmup*kzsp(ss-2)*(szsp(ss-2) - sz1)*(szsp(ss-2) - sz2))
       res = psi*(1._dpk - mup*M2)/((mup - szsp(ss-2))* &
            kmmup*kzsp(ss-2)*(szsp(ss-2) - sz1)*(szsp(ss-2) - sz2))
       do ii = 1, Nzsp
          if (ii .NE. (ss-2)) then
             res = res/(szsp(ss-2) - szsp(ii))
          end if
       end do
       do jj = 1, Npsp
          res = res*(szsp(ss-2) - spsp(jj))
       end do
          
       if (r <= h) then
          fn = (1._dpk - szsp(ss-2)*M2)**2*Trs(r,szsp(ss-2),1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - szsp(ss-2)*M2)**2*Trs(r,szsp(ss-2),2)
       else
          fn = (1._dpk - szsp(ss-2)*M3)**2*Trs(r,szsp(ss-2),3)
       end if

       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*wr*wr*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            wr*szsp(ss-2)*z)*res*fn

!!$       else
!!$          residuepr = 0.
!!$
!!$       end if

    end if

    
  END FUNCTION residuepr


  FUNCTION residuepot(r,z,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Same as above but for computing velocity potential
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)            :: r, z
    integer              :: ss, ii, jj
    complex(dpk)         :: residuepot, res, fn
   

    if (ss == 1) then

!!$       if (z > 0.) then
!!$       res = psi*(1._dpk - mup*M2)*(sz1-sp1)*(sz1-CONJG(sp1))/ &
!!$            ((mup - sz1)*kmmup*kpsz1*(sz1-CONJG(sz1))*(sz1-sz2)*(sz1-CONJG(sz2)))
!!$       res = psi*(1._dpk - mup*M2)*(sz1 - sp1)/((mup - sz1)*kmmup*kpsz1*(sz1 - sz2))
       res = psi*(1._dpk - mup*M2)/((mup - sz1)*kmmup*kpsz1*(sz1 - sz2))
       if (vortswitch .EQ. 0) then
          do ii = 1, Nzsp
             res = res/(sz1 - szsp(ii))
          end do
          do jj = 1, Npsp
             res = res*(sz1 - spsp(jj))
          end do
       end if

       if (r <= h) then
          fn = (1._dpk - sz1*M2)**2/(1._dpk - sz1*M1)*Trs(r,sz1,1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - sz1*M2)*Trs(r,sz1,2)
       else
          fn = (1._dpk - sz1*M3)*Trs(r,sz1,3)
       end if
       
       residuepot = wr*EXP(CMPLX(0.,1._dpk,kind=dpk)*wr*sz1*z)*res*fn
          
!!$       else
!!$          residuepot = 0.
!!$
!!$       end if
          
    elseif (ss == 2) then

!!$       if (z > 0.) then
!!$       res = psi*(1._dpk - mup*M2)*(sz2-sp1)*(sz2-CONJG(sp1))/ &
!!$            ((mup - sz2)*kmmup*kpsz2*(sz2-CONJG(sz1))*(sz2-sz1)*(sz2-CONJG(sz2)))
!!$       res = psi*(1._dpk - mup*M2)*(sz2 - sp1)/((mup - sz2)*kmmup*kpsz2*(sz2 - sz1))
       res = psi*(1._dpk - mup*M2)/((mup - sz2)*kmmup*kpsz2*(sz2 - sz1))
       if (vortswitch .EQ. 0) then
          do ii = 1, Nzsp
             res = res/(sz1 - szsp(ii))
          end do
          do jj = 1, Npsp
             res = res*(sz1 - spsp(jj))
          end do
       end if
          
       if (r <= h) then
          fn = (1._dpk - sz2*M2)**2/(1._dpk - sz2*M1)*Trs(r,sz2,1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - sz2*M2)*Trs(r,sz2,2)
       else
          fn = (1._dpk - sz2*M3)*Trs(r,sz2,3)
       end if

       residuepot = wr*EXP(CMPLX(0.,1._dpk,kind=dpk)*wr*sz2*z)*res*fn

!!$       else
!!$          residuepot = 0.
!!$
!!$       end if

    else

!!$       if (z > 0.) then
!!$      res = psi*(1._dpk - mup*M2)*(szsp(ss-2)-sp1)*(szsp(ss-2)-CONJG(sp1))/ &
!!$            ((mup - szsp(ss-2))*kmmup*kzsp(ss-2)*(szsp(ss-2)-CONJG(sz1))* &
!!$            (szsp(ss-2)-sz1)*(szsp(ss-2) - sz2)*(szsp(ss-2)-CONJG(sz2)))
!!$       res = psi*(1._dpk - mup*M2)*(szsp(ss-2) - sp1)/((mup - szsp(ss-2))* &
!!$            kmmup*kzsp(ss-2)*(szsp(ss-2) - sz1)*(szsp(ss-2) - sz2))
       res = psi*(1._dpk - mup*M2)/((mup - szsp(ss-2))* &
            kmmup*kzsp(ss-2)*(szsp(ss-2) - sz1)*(szsp(ss-2) - sz2))
       do ii = 1, Nzsp
          if (ii .NE. (ss-2)) then
             res = res/(szsp(ss-2) - szsp(ii))
          end if
       end do
       do jj = 1, Npsp
          res = res*(szsp(ss-2) - spsp(jj))
       end do
          
       if (r <= h) then
          fn = (1._dpk - szsp(ss-2)*M2)**2/(1._dpk - szsp(ss-2)*M1)*Trs(r,szsp(ss-2),1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - szsp(ss-2)*M2)*Trs(r,szsp(ss-2),2)
       else
          fn = (1._dpk - szsp(ss-2)*M3)*Trs(r,szsp(ss-2),3)
       end if

       residuepot = wr*EXP(CMPLX(0.,1._dpk,kind=dpk)*wr*szsp(ss-2)*z)*res*fn

!!$       else
!!$          residuepot = 0.
!!$
!!$       end if

    end if

    
  END FUNCTION residuepot


  SUBROUTINE checkloc(z,switch)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Subroutine to check the location of the instability zeros and pole wrt to the contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: z
    real(dpk)       :: zx, zy, zfy
    integer         :: switch

    zx = REAL(z)
    zy = AIMAG(z)

    call yfunc_p(zx,zfy,REAL(szi(3)),REAL(szi(4)-szi(3)),AIMAG(szi(4)-szi(3)))
!!$    zfy = AIMAG(c1)

    if(zy > zfy) then
       switch = 0  !! inside
    else
       switch = 1  !! outside
    end if


  END SUBROUTINE checkloc

  
  FUNCTION bessj(z,nu,s)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! All the Bessel and similar functions use the NAG library routines. The present
!! versions also accept negative mode numbers, \nu. In addition, an accuracy check
!! ensures they compute the corresponding asymptotic forms whenever needed instead
!! of giving garbage values, which NAG routines give beyond a prescribed accuracy.
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
    
    complex(dpk)      :: z, bessj, psi, bj, by, cwrk
    real(dpk)         :: nu
    integer           :: s, nz, ifail
    
    if (nu < 0.) then
       if (s == 0) then
          call S17DEF(-nu,z,1,'U',bj,nz,ifail)
          call S17DCF(-nu,z,1,'U',by,nz,cwrk,ifail)
          bessj = bj*COS(-PI*nu) - by*SIN(-PI*nu)
       else
          call S17DEF(-nu,z,1,'S',bj,nz,ifail)
          call S17DCF(-nu,z,1,'S',by,nz,cwrk,ifail)
          bessj = bj*COS(-PI*nu) - by*SIN(-PI*nu)
       end if
    else
       if (s == 0) then
          call S17DEF(nu,z,1,'U',bessj,nz,ifail)
       else
          call S17DEF(nu,z,1,'S',bessj,nz,ifail)
       end if
    end if

   
  END FUNCTION bessj


  FUNCTION abessj(z,nu)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! These are asymptotic forms (starting with an "a")
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)      :: z, abessj, psi
    real(dpk)         :: nu


    psi = z - PI*(nu/2._dpk + 0.25_dpk)
    abessj = SQRT(2._dpk/(PI*z))*(COS(psi) - (4._dpk*nu*nu - 1._dpk)/ &
            (8._dpk*z)*SIN(psi))


  END FUNCTION abessj


  FUNCTION bessy(z,nu,s)

    complex(dpk)      :: z, bessy, psi, bj, by, cwrk
    real(dpk)         :: nu
    integer           :: s, nz, ifail

    if (nu < 0.) then
       if (s == 0) then
          call S17DEF(-nu,z,1,'U',bj,nz,ifail)
          call S17DCF(-nu,z,1,'U',by,nz,cwrk,ifail)
          bessy = by*COS(-PI*nu) + bj*SIN(-PI*nu)
       else
          call S17DEF(-nu,z,1,'S',bj,nz,ifail)
          call S17DCF(-nu,z,1,'S',by,nz,cwrk,ifail)
          bessy = by*COS(-PI*nu) + bj*SIN(-PI*nu)
       end if
    else
       if (s == 0) then
          call S17DCF(nu,z,1,'U',bessy,nz,cwrk,ifail)
       else
          call S17DCF(nu,z,1,'S',bessy,nz,cwrk,ifail)
       end if
    end if

    
  END FUNCTION bessy


  FUNCTION abessy(z,nu)

    complex(dpk)      :: z, abessy, psi
    real(dpk)         :: nu


    psi = z - PI*(nu/2._dpk + 0.25_dpk)
    abessy = SQRT(2._dpk/(PI*z))*(SIN(psi) + (4._dpk*nu*nu - 1._dpk)/ &
            (8._dpk*z)*COS(psi))


  END FUNCTION abessy


  FUNCTION hank1(z,nu,s)

    complex(dpk)      :: z, hank1, psi, h1
    real(dpk)         :: nu
    integer           :: s, nz, ifail

    if (nu < 0.) then
       if (s == 0) then
          call S17DLF(1,-nu,z,1,'U',h1,nz,ifail)
          hank1 = EXP(-nu*PI*CMPLX(0._dpk,1._dpk,kind=dpk))*h1
       else
          call S17DLF(1,-nu,z,1,'S',h1,nz,ifail)
          hank1 = EXP(-nu*PI*CMPLX(0._dpk,1._dpk,kind=dpk))*h1
       end if
    else
       if (s == 0) then
          call S17DLF(1,nu,z,1,'U',hank1,nz,ifail)
       else
          call S17DLF(1,nu,z,1,'S',hank1,nz,ifail)
       end if
    end if

    
  END FUNCTION hank1


  FUNCTION ahank1(z,nu)

    complex(dpk)      :: z, ahank1, psi
    real(dpk)         :: nu


    psi = z - PI*(nu/2._dpk + 0.25_dpk)
    ahank1 = SQRT(2._dpk/(PI*z))*(1._dpk + CMPLX(0._dpk,1._dpk,kind=dpk)*(4._dpk*nu*nu - 1._dpk)/ &
            (8._dpk*z))*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*psi)

    print*,'ahank1:',ahank1
    print*,'z:',z


  END FUNCTION ahank1


  FUNCTION dbessj(z,nu,s)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! These are the first derivatives (starting with a "d")
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: z, dbessj, psi
    real(dpk)       :: nu
    integer         :: s

    dbessj = 0.5_dpk*(bessj(z,nu-1,s) - bessj(z,nu+1,s))


  END FUNCTION dbessj


  FUNCTION adbessj(z,nu)

    complex(dpk)      :: z, adbessj, psi
    real(dpk)         :: nu


    psi = z - PI*(nu/2._dpk + 0.25_dpk)
    adbessj = SQRT(2._dpk/(PI*z))*(-SIN(psi) - (4._dpk*nu*nu + 3._dpk)/ &
            (8._dpk*z)*COS(psi))


  END FUNCTION adbessj


  FUNCTION dbessy(z,nu,s)

    complex(dpk)    :: z, dbessy, psi
    real(dpk)       :: nu
    integer         :: s

    dbessy = 0.5_dpk*(bessy(z,nu-1,s) - bessy(z,nu+1,s))


  END FUNCTION dbessy


  FUNCTION adbessy(z,nu)

    complex(dpk)      :: z, adbessy, psi
    real(dpk)         :: nu


    psi = z - PI*(nu/2._dpk + 0.25_dpk)
    adbessy = SQRT(2._dpk/(PI*z))*(COS(psi) - (4._dpk*nu*nu + 3._dpk)/ &
            (8._dpk*z)*SIN(psi))


  END FUNCTION adbessy


  FUNCTION dhank1(z,nu,s)

    complex(dpk)    :: z, dhank1, psi
    real(dpk)       :: nu
    integer         :: s

    dhank1 = 0.5_dpk*(hank1(z,nu-1,s) - hank1(z,nu+1,s))


  END FUNCTION dhank1


  FUNCTION adhank1(z,nu)

    complex(dpk)      :: z, adhank1, psi
    real(dpk)         :: nu


    psi = z - PI*(nu/2._dpk + 0.25_dpk)
    adhank1 = SQRT(2._dpk/(PI*z))*(CMPLX(0._dpk,1._dpk,kind=dpk) - (4._dpk*nu*nu + 3._dpk)/ &
            (8._dpk*z))*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*psi)

    print*,'adhank1:',adhank1
    print*,'z:',z


  END FUNCTION adhank1


END PROGRAM main
