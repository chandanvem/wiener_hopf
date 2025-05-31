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
  complex(dpk), allocatable, dimension(:)    :: ker_int_points  !! location of the starting pts
  complex(dpk), allocatable, dimension(:)    :: intpanel  !! integral value at each panel
  integer, allocatable, dimension(:)         :: Npanel  !! number of times a panel is divided
  integer                                    :: num_ker_pts_loop !! number of starting pts for kernel contour
  integer                                    :: tot_ker_points, tot_IFT_pts
  integer                                    :: prswitch, restart, reflswitch, farswitch, vortswitch
  integer                                    :: num_zeros_s1_s2, num_poles_s1_s2  !! num  zeros & poles in between s_1- s_2
  integer                                    :: num_sup_zeros, num_sup_poles  !! number of supersonic zeros & poles
  integer                                    :: Nphi  !! polar mesh resolution for directivity computation
  complex(dpk), allocatable, dimension(:)    :: def_pts_IFT_cntr  !! the points defining the IFT contour
  complex(dpk), allocatable, dimension(:)    :: def_pts_ker_cntr  !! the points defining the Kernel contour
  complex(dpk), allocatable, dimension(:)    :: zeros_list_bw_s1_s2, poles_list_bw_s1_s2, sup_zeros_list, sup_poles_list
  complex(dpk), allocatable, dimension(:)    :: fplusz, fplusz_temp
  complex(dpk)                               :: k_minus_at_mu_plus, k_plus_sz1, k_plus_sz2, kpsp1
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
  real(dpk)                                  :: kapT  !! sqrt(T1/T0)
  real(dpk)                                  :: kap_rho  !! rho1/rho0
  real(dpk)                                  :: theta, cont_cross_over_pt, w0, delr, deli, offset
  complex(dpk)                               :: omega_i, omega_r
  real(dpk)                                  :: Rmax, Rmin, Zmin, Zmax !! dimensions of the physical domain
  real(dpk), allocatable, dimension(:)       :: R, Z
  real(dpk)                                  :: Zo  !! start point of incident vorticity (negative)
  integer                                    :: Nmeshr, Nmeshz !! 400X400 -> Seg Fault!!
  real(dpk)                                  :: azim_mode
  real(dpk)                                  :: vs_param_gamma  !! Vortex shedding parameter
  complex(dpk)                               :: KH_zero_1 !! instability zero 1
  complex(dpk)                               :: KH_zero_2 !! instability zero 2
  complex(dpk)                               :: KH_pole_1 !! instability pole 1
  complex(dpk)                               :: mu_plus !! axial wave no (of incident wave)
  complex(dpk)                               :: mu0 !! incident axial wave no (for incident vort)
  complex(dpk)                               :: resp!! portion of residue (for inc vort)
  complex(dpk)                               :: psi
  integer                                    :: num_IFT_pts_loop !! number of points for IFT contour
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
	    print *, 'initialize:  vortswitch=', vortswitch

	    read(1,*) M1  !! core jet Mach number
	    print *, 'initialize:  M1=', M1

	    read(1,*) M2  !! coflow Mach number
	    print *, 'initialize:  M2=', M2

	    read(1,*) M3  !! ambient flow Mach number
	    print *, 'initialize:  M3=', M3

	    read(1,*) h   !! the width of the core jet, h = Ri/Ro
	    print *, 'initialize:  width of the core jet h=', h

	    if ((vortswitch == 1) .OR. (vortswitch == 2)) then
	       read(1,*) Zo  !! starting point of inc instability (-ve)
	    end if

	    read(1,*) w0  !! the Helmholtz number
	    print *, 'initialize:  Helmholtz number omega=', w0

	    read(1,*) kapT  !! sqrt(temp ratio)
	    print *, 'initialize:  Sqrt of temperature ratio=', kapT
	 
	    read(1,*) kap_rho  !! density ratio
	    print*, 'initialize:  Density ratio=', kap_rho

	    read(1,*) azim_mode  !! the circumferential mode no
	    print*, 'initialize:  Azimuthal wavenumber (azim_mode)=', azim_mode

	    read(1,*) KH_zero_1
	    print*, 'initialize:  KH_zero_1 (First instability zero) =', KH_zero_1


	    read(1,*) KH_zero_2
	    print*, 'initialize:  KH_zero_2 (Second instability zero) =', KH_zero_2


	    read(1,*) KH_pole_1
	    print*, 'initialize:  KH_pole_1 (Instability pole) =', KH_pole_1


	    if ((vortswitch == 1) .OR. (vortswitch == 2)) then
		       read(1,*) mu0  !! upstream inc acoustic mode
		       print*, 'initialize:  Incident Vorticity mode activated'
	    else
		       read(1,*) mu_plus
		       print*, '================== Non-incident vorticity mode ====================='
		       print*,''
		       print*, 'initialize:  mu for the incident mode =',mu_plus
	    end if

	    read(1,*) offset  !! the offset between the two integration contours
	    print*,'initialize:  Offset between the IFT and K split contours=', offset

	    read(1,*) tol  !! the tolerance of the adaptive contour integration routine
	    print*,'initialize:  Tolerance for adaptive contour integration =', tol

	    read(1,*) num_zeros_s1_s2  !! zeros needed to be removed from the kernel
	    print*,'initialize:  Number of zeros to be removed from kernel =', num_zeros_s1_s2

	    read(1,*) num_poles_s1_s2  !! poles needed to be removed from the kernel
	    print*,'initialize:  Number of poles to be removed from kernel =', num_poles_s1_s2

	    read(1,*) num_sup_zeros  !! supersonic zeros
	    print*,'initialize:  Number of supersonic zeros =', num_sup_zeros

	    read(1,*) num_sup_poles  !! supersonic poles
	    print*,'initialize:  Number of supersonic poles =', num_sup_poles

	    if (vortswitch == 2) then
	       if (num_sup_poles == 0) then
			  print*, "initialize:  For vortswitch mode 2, at least one upstream supersonic pole needed! Exiting..."
			  STOP
	       end if
	    end if

	    print*,''
	    print*,'','====== Mesh and contour parameters ======',''
	    
	    read(1,*) num_ker_pts_loop
	    print*,'initialize:  Number of kernel points in each loop =', num_ker_pts_loop

	    read(1,*) theta  !! the stretching parameter for kernel contour meshpoints
	    print*,'initialize:  Stretching parameter for kernel contour =', theta

	    read(1,*) num_IFT_pts_loop
	    print*,'initialize:  Number of IFT points in each loop num_IFT_pts_loop =', num_IFT_pts_loop

	    read(1,*) Rmin
	    read(1,*) Rmax
	    read(1,*) Zmin
	    read(1,*) Zmax
	   
	    read(1,*) Nmeshr
	    read(1,*) Nmeshz

	    print*,'initialize:  Dimensions of domain in R = [',Rmin,Rmax,']' 
	    print*,'initialize:  Dimensions of domain in Z = [',Zmin,Zmax,']' 
	    print*,'initialize:  Nr x Nz =',Nmeshr,'x',Nmeshz

	    print*,'=============================='

	    read(1,*) asymplim
	    print*,'initialize:  Asymptotic limit of ??? = ',asymplim

	    read(1,*) asymplim1
	    print*,'initialize:  Asymptotic limit of ??? = ',asymplim1

	    read(1,*) vs_param_gamma  !! Default = 1.0 (0,1)
	    print*,'initialize:  Vortex shedding parameter gamma = ',vs_param_gamma 
	      
	!!============================
	    read(1,*) prswitch  !! 0 = Potential; 1 = Pressure

	    if (prswitch == 0) then
	       print*,'initialize:  Solution in potential mode'
	    elseif (prswitch == 1) then
	       print*,'initialize:  Solution in pressure mode '
	    end if
	!!============================
	    read(1,*) reflswitch  !! 1 = Reflection mode: incident mode not added
	    if (reflswitch == 1) then
	       print*,'initialize:  Reflection mode: incident mode not added'
	    else 
	       print*,'initialize:  Reflection mode: incident mode added'
	    end if
	!!============================
	    read(1,*) farswitch  !! 1 = Far-field mode: compute directivity; 2 = 1 + nearfield of sup zeros in polar mode

	    if ((farswitch == 1) .OR. (farswitch == 2)) then
	       read(1,*) Nphi
	       allocate(cnanglei(num_sup_zeros+2))
	    end if

	    if ((farswitch == 2) .AND. num_sup_zeros <= 0) then
	       print*, "For farswitch mode 2, non-zero number of supersonic zero needed! Exiting..."
	       STOP
	    end if

	!!=====================================================

	    read(1,*) restart  !! restart status: 0 = fresh job, i.e., no "fplus_part.out" exists

	!! read the zeros & poles that need to be excluded, if any:
	!! USE: when the effect of particular modes need to be studied individually
	    
	    if (num_zeros_s1_s2 > 0) then
		allocate(zeros_list_bw_s1_s2(num_zeros_s1_s2))
		do i1 = 1, num_zeros_s1_s2
		   read(1,*) zeros_list_bw_s1_s2(i1)
	!          print*, zeros_list_bw_s1_s2(i1)
		end do
	    end if


	    if (num_poles_s1_s2 > 0) then
	       allocate(poles_list_bw_s1_s2(num_poles_s1_s2))

	       do i1 = 1, num_poles_s1_s2
		  read(1,*) poles_list_bw_s1_s2(i1)
	!         print*, poles_list_bw_s1_s2(i1)
	       end do
	    end if

	    if (num_sup_zeros > 0) then
		allocate(sup_zeros_list(num_sup_zeros))
		do i1 = 1, num_sup_zeros
		  read(1,*) sup_zeros_list(i1)
	!         print*, sup_zeros_list(i1)
		end do
	   end if

	   
	   if (num_sup_poles > 0) then 
	       allocate(sup_poles_list(num_sup_poles))
	       do i1 = 1, num_sup_poles
		 read(1,*) sup_poles_list(i1)
	!        print*, sup_poles_list(i1)
	       end do
	  end if

    close(1)

!!================================================================

    PI = 4._dpk*ATAN(1.)
    delr = 0._dpk   !! for real omega
    deli = PI/2._dpk  !! for purely imaginary omega

!! defining the omegas (Helmholtz numbers):

    omega_r = ABS(w0)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*delr)  !! real

!! whether to compute velocity potential or pressure:

    if (prswitch == 0) then
       if ((farswitch == 1) .OR. (farswitch == 2)) then
		  print*,''
		  print*,'initialize: Far-field is computed ONLY in pressure mode. Exiting...'
		  STOP
       else
		  print*,''
		  print*,'initialize: Computing Potential...'
       end if
    
    else
	       print*,''
	       print*,'initialize: Computing Pressure...'
    end if

 !! now allocate for the contour points:

    allocate(def_pts_IFT_cntr(5))
    allocate(def_pts_ker_cntr(5))

!! read the contour points, only IFT contour is read in full;
!! the end points of the kernel integration contour read while
!! the rest are obtained via the "offset" parameter
!! def_pts_IFT_cntr(1), def_pts_IFT_cntr(5) are end points; def_pts_IFT_cntr(2), def_pts_IFT_cntr(4) are where the
!! algebraic contour starts & ends; def_pts_IFT_cntr(3) is where it crosses
!! real axis

    open(1,file='input.iftpts')

	    read(1,*) t
	    ai = t
	    read(1,*) def_pts_IFT_cntr(2)
	    read(1,*) t
	    def_pts_IFT_cntr(3) = CMPLX(t,0._dpk,kind=dpk)
	    read(1,*) def_pts_IFT_cntr(4)
	    read(1,*) t
	    af = t

	    call get_y_alg_int_contour(ai,ai_im,REAL(def_pts_IFT_cntr(3)),REAL(def_pts_IFT_cntr(2)-def_pts_IFT_cntr(3)), &
					AIMAG(def_pts_IFT_cntr(2)-def_pts_IFT_cntr(3)))

	    call get_y_alg_int_contour(af,af_im,REAL(def_pts_IFT_cntr(3)),REAL(def_pts_IFT_cntr(4)-def_pts_IFT_cntr(3)), &
					     AIMAG(def_pts_IFT_cntr(4)-def_pts_IFT_cntr(3)))

	    def_pts_IFT_cntr(1) = CMPLX(ai,ai_im,kind=dpk)  
	    def_pts_IFT_cntr(5) = CMPLX(af,af_im,kind=dpk)

	    read(1,*) t
	    aki = t
	    read(1,*) t
	    akf = t

    close(1)

!! now determine the kernel integration contour points:
!! NOTE: the kernel contour lies above the IFT one and thus it's further away from the real
!! axis when above, but closer to it when below

    def_pts_ker_cntr(2) = def_pts_IFT_cntr(2) + CMPLX(0._dpk,offset,kind=dpk)  !! Kernel points
    def_pts_ker_cntr(3) = def_pts_IFT_cntr(3) + CMPLX(5._dpk*offset,0._dpk,kind=dpk)
    def_pts_ker_cntr(4) = def_pts_IFT_cntr(4) + CMPLX(0._dpk,offset,kind=dpk)

    call get_y_alg_int_contour(aki,aki_im,REAL(def_pts_ker_cntr(3)), &
                           REAL(def_pts_ker_cntr(2)-def_pts_ker_cntr(3)),AIMAG(def_pts_ker_cntr(2)-def_pts_ker_cntr(3)))
    call get_y_alg_int_contour(akf,akf_im,REAL(def_pts_ker_cntr(3)), &
                           REAL(def_pts_ker_cntr(4)-def_pts_ker_cntr(3)),AIMAG(def_pts_ker_cntr(4)-def_pts_ker_cntr(3)))
    
    def_pts_ker_cntr(1) = CMPLX(aki,aki_im,kind=dpk)
    def_pts_ker_cntr(5) = CMPLX(akf,akf_im,kind=dpk)

!! find where the contour crosses the real axis (at "cont_cross_over_pt"):

    cont_cross_over_pt = REAL(def_pts_ker_cntr(3))

!! here both the instability zeros and the pole need to lie outside (under) the contour,
!! since we use the residue theorem to compute them anyway:

    call checkloc(KH_zero_1,i1) 
    call checkloc(KH_zero_2,i2)
    call checkloc(KH_pole_1,i3)

!! if any of them lie inside ask to redefine the contour:

    if (vortswitch .EQ. 0) then

       if(i1==0 .OR. i2==0 .OR. i3==0) then
		  print*,''
		  print*,'initialize: Redefine the contour: One of instability zero/pole is INSIDE!'
		  print*,''
		  if (i1==0) print*, 'initialize: Zero1 is inside'
		  if (i2==0) print*, 'initialize: Zero2 is inside'
		  if (i3==0) print*, 'initialize: Pole1 is inside'
		  print*,''
		  STOP
       end if

    else

       if(i1==1 .OR. i2==0 .OR. i3==1) then
		  print*,''
		  print*,'initialize: Redefine the contour:'
		  print*,''
		  if (i1==1) print*, 'initialize: Zero1 is outside'
		  if (i2==0) print*, 'initialize: Zero2 is inside'
		  if (i3==1) print*, 'initialize: Pole1 is outside'
		  print*,''
		  STOP
       end if

    end if

    do i1 = 1, num_sup_zeros
       call checkloc(sup_zeros_list(i1),i2)
       if (vortswitch .EQ. 0) then
          if(i2 == 0) then
		     print*,''
		     print*,'initialize: Redefine the contour: The following supersonic zero is inside:'
		     print*,sup_zeros_list(i1)
		  end if
       else
          if(i2 == 1) then
		     print*,''
		     print*,'initialize: Redefine the contour: The following supersonic zero is outside:'
		     print*,sup_zeros_list(i1)
          end if
       end if
    end do

    do i1 = 1, num_sup_poles
       call checkloc(sup_poles_list(i1),i2)
       if (vortswitch .EQ. 0) then
          if(i2 == 0) then
		     print*,''
		     print*,'initialize: Redefine the contour: The following supersonic pole is inside:'
		     print*,sup_poles_list(i1)
          end if
       else
          if(i2 == 1) then
		     print*,''
		     print*,'initialize: Redefine the contour: The following supersonic pole is outside:'
		     print*,sup_poles_list(i1)
          end if
       end if
    end do

!! compute the various contour parameters:

    tot_ker_points = 2*num_ker_pts_loop + 3  !! total number of kernel contor pts

    tot_IFT_pts = 2*num_IFT_pts_loop + 3  !! total number of IFT contor pts

!! print the contour points:

    write(*,'(/A26)'),'initialize: The Integration Contours:'
    write(*,'(/A27)'),'1. The kernel integration contour:'
    write(*,'(/A14)'),'key points:->'
    do i1 = 1, 5
       write(pos,"(I2,' =')"),i1
       write(*,'(/1X,A5,2F15.6)'), 'i'//pos, def_pts_ker_cntr(i1)
    end do

    write(*,'(/A30)'),'2. The inv Fourier transform contour:'
    write(*,'(/A14)'),'key points:->' 
    do i1 = 1, 5
       write(pos,"(I2,' =')"),i1
       write(*,"(/1X,A5,2F15.6)"), 'i'//pos, def_pts_IFT_cntr(i1) 
    end do

   print '(A, I0)', 'initialize: total number of IFT points = ', tot_IFT_pts 
   

  END SUBROUTINE initialize


  SUBROUTINE solve

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Calls the various other subroutines which do the actual job
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer     :: index

    print*,'solve: Defining contours for kernel and IFT integration...'

    call definecontours

    print*,'solve: Done...' 
    print*,'solve: Starting pre compute routine...'
    print*,'solve: Computing K+ at the zeros and K- at mu+:'

    call precompute

    if ((farswitch == 1) .OR. (farswitch == 2)) then

	       print*,''
	       print*,'solve: Near-field computation only:'
	       print*,'' 

    else

	       call meshgrid

	       if (restart .NE. 0) then
		           call read_fplus(index)  !! assumes the fplus_part.out file to be present
		           call compute_fplus(restart,index)
	       else
                 print*,'solve: Now computing F+. This takes a while:'
		           call compute_fplus(restart,0)
	       end if

	       print*,''
	       print*,'solve: Now starting the IFT:'
	       print*,''

	       call computeift

    end if

    print*,'Done!'


  END SUBROUTINE solve


  SUBROUTINE definecontours

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Computes the panel lengths for the integration contours
!! 2. Calls subrtn "initialize_contour" which does the the actual contour defining
!! 3. Writes the contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                     :: panel_len_left, panel_len_right  !! length of each panel
    integer                       :: i

    allocate(ker_int_points(tot_ker_points))

!! length of each panel (kernel contour):

    panel_len_left =  (REAL(def_pts_ker_cntr(3))-REAL(def_pts_ker_cntr(1)))/(num_ker_pts_loop+1)  !! panel length in the left section
    panel_len_right = (REAL(def_pts_ker_cntr(5))-REAL(def_pts_ker_cntr(3)))/(num_ker_pts_loop+1)  !! panel length in the right section

    call initialize_contour(def_pts_ker_cntr,panel_len_left,panel_len_right,num_ker_pts_loop,1,ker_int_points)

    print*,''
    print*,'definecontours: Writing kernel integration points to file:'

    open(10,file='initialpoints.out',form='FORMATTED')
   
	    do i = 1,tot_ker_points
	       write(10,'(I10,2F30.20)') i, ker_int_points(i)
	    end do
   
    close(10)

    allocate(iftpoints(tot_IFT_pts))

!! length of each panel (IFT contour):

    panel_len_left = (REAL(def_pts_IFT_cntr(3))-REAL(def_pts_IFT_cntr(1)))/(num_IFT_pts_loop+1)
    panel_len_right = (REAL(def_pts_IFT_cntr(5))-REAL(def_pts_IFT_cntr(3)))/(num_IFT_pts_loop+1)

    call initialize_contour(def_pts_IFT_cntr,panel_len_left,panel_len_right,num_IFT_pts_loop,2,iftpoints)

    print*,'definecontours: Writing IFT integration points to file:'

    open(10,file='iftpoints.out',form='FORMATTED')

	    do i = 1,tot_IFT_pts
	       write(10,'(I10,2F30.20)') i, iftpoints(i)
	    end do

    close(10)


  END SUBROUTINE definecontours


  SUBROUTINE initialize_contour(def_pts,left_panel_length,right_panel_length,num_pts_per_loop,sw,init_pts_combined)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Initialize the integration contours
!! 2. First compute the individual segments
!! 3. Combine the individual segments to get the full contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), allocatable, dimension(:)       :: init_pts_left, init_pts_right
    complex(dpk), dimension(2*num_pts_per_loop+3) :: init_pts_combined
    complex(dpk), dimension(5)                    :: def_pts
    real(dpk)                                     :: left_panel_length, right_panel_length
    integer                                       :: num_pts_per_loop, sw


    allocate(init_pts_left(num_pts_per_loop+2))
    allocate(init_pts_right(num_pts_per_loop+2))

    call initsubdiv(REAL(def_pts(1)),def_pts(2)-REAL(def_pts(3)),REAL(def_pts(3)), &
                         num_pts_per_loop,left_panel_length,sw,1,init_pts_left)

    call initsubdiv(REAL(def_pts(3)),def_pts(4)-REAL(def_pts(3)),REAL(def_pts(5)), &
                         num_pts_per_loop,right_panel_length,sw,2,init_pts_right)

    call combine(num_pts_per_loop,init_pts_left,init_pts_right,init_pts_combined)

  END SUBROUTINE initialize_contour


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
          call create_stretched_grid(r,p,N+2,1,xi)  

       case(2)
          call create_stretched_grid(p,r,N+2,2,xi)  

       end select

    else
       xi(1) = REAL(p)

    end if


    do i = 1,N+1
       select case (ss)

       case(1)
          call get_y_alg_int_contour(xi(i),yi(i),r,REAL(q),AIMAG(q))  !! the imaginary points
          if (sw == 2) xi(i+1) = xi(i) + len  !! for ift contour the points are equispaced
          
       case(2)
          call get_y_alg_int_contour(xi(i),yi(i),p,REAL(q),AIMAG(q))
          if (sw == 2) xi(i+1) = xi(i) + len

       end select

       zi(i) = CMPLX(xi(i),yi(i),kind=dpk)
    end do


    select case (ss)  !! the right end point

        case(1)
           call get_y_alg_int_contour(xi(N+2),yi(N+2),r,REAL(q),AIMAG(q))
       
        case(2)
           call get_y_alg_int_contour(xi(N+2),yi(N+2),p,REAL(q),AIMAG(q))

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


  SUBROUTINE create_stretched_grid(a,b,N,ss,xst)

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


  END SUBROUTINE create_stretched_grid


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

    print*,'meshgrid: Writing mesh grids to file mesh.out'
    open(1,file='mesh.out',form='UNFORMATTED')
    write(1) Nmeshz,Nmeshr
    write(1) ((Z(i),i=1,Nmeshz),j=1,Nmeshr),((R(j),&
         i=1,Nmeshz),j=1,Nmeshr)
    close(1)   

   print*,'meshgrid: Writing mesh grids to file mesh.r'
    open(20,file='mesh.r')
    do i = 1, Nmeshr
       write(20,'(I10,2X,F20.10)'), i, R(i)
    end do
    close(20)


    print*,'meshgrid: Writing mesh grids to file mesh.z'
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
    
    allocate(pr(tot_IFT_pts-1))
    allocate(prsum(Nmeshr,Nmeshz))
    allocate(pressure(Nmeshr,Nmeshz))
    allocate(acoupressure(Nmeshr,Nmeshz))
    allocate(inspressure1(Nmeshr,Nmeshz))
    allocate(inspressure2(Nmeshr,Nmeshz))
    allocate(totpressure(Nmeshr,Nmeshz))
    allocate(initpressure(Nmeshr,Nmeshz))

    if (num_sup_zeros > 0) then
       allocate(supinspressure(num_sup_zeros,Nmeshr,Nmeshz))
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

          do k = 1, tot_IFT_pts-1

             call IFT_trapz_int(R(i),Z(j),k,switch,pr(k))
             
          end do

          open(10,file='DataDump/Ift/intift.'//pos,form='FORMATTED')
          do k = 1, tot_IFT_pts-1
             write(10,'(I5,2E20.10)') k,pr(k)
          end do
          close(10)

          call sum_panel_contribution_IFT(pr,prsum(i,j))
          
          if (prswitch == 0) then  !! compute velocity potential
             
             pressure(i,j) = omega_r/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk))*prsum(i,j)

             acoupressure(i,j) = pressure(i,j)  !! acoustic part

             if (vortswitch .EQ. 0) then
                inspressure1(i,j) = residuepot(R(i),Z(j),1)  !! inner instability wave 
             end if

             inspressure2(i,j) = residuepot(R(i),Z(j),2)  !! outer instability wave

             if (vortswitch .EQ. 0) then
                do k = 1, num_sup_zeros
                   supinspressure(k,i,j) = residuepot(R(i),Z(j),k+2)
                end do
             end if
             
             if (vortswitch .EQ. 0) then
                totpressure(i,j) = acoupressure(i,j) + inspressure1(i,j) + inspressure2(i,j)  !! total part 
             else
                totpressure(i,j) = acoupressure(i,j) + inspressure2(i,j)
             end if

             if (vortswitch .EQ. 0) then
                do k = 1, num_sup_zeros
                   totpressure(i,j) = totpressure(i,j) + supinspressure(k,i,j)
                end do
             end if

             print*,'pressure:',REAL(totpressure(i,j))

          else  !! compute pressure

             if (reflswitch == 1) then

                pressure(i,j) = omega_r*omega_r/(2._dpk*PI)*prsum(i,j)   !! incident wave NOT added

             else

                pressure(i,j) = omega_r*omega_r/(2._dpk*PI)*prsum(i,j) + compute_psi_incident(R(i),Z(j),switch) 
                                                !! the incident wave added

             end if
          
             acoupressure(i,j) = pressure(i,j)  !! acoustic part

             if (vortswitch .EQ. 0) then
                inspressure1(i,j) = residuepr(R(i),Z(j),1)  !! inner instability wave
             end if

             inspressure2(i,j) = residuepr(R(i),Z(j),2)  !! outer instability wave

             if (vortswitch .EQ. 0) then
                do k = 1, num_sup_zeros
                   supinspressure(k,i,j) = residuepr(R(i),Z(j),k+2)
                end do
             end if

             if (vortswitch .EQ. 0) then
                totpressure(i,j) = acoupressure(i,j) + inspressure1(i,j) + inspressure2(i,j)  !! total part
             else
                totpressure(i,j) = acoupressure(i,j) + inspressure2(i,j)
             end if

             if (vortswitch .EQ. 0) then
                do k = 1, num_sup_zeros
                   totpressure(i,j) = totpressure(i,j) + supinspressure(k,i,j)
                end do
             end if

             print*,'pressure:',REAL(totpressure(i,j))

          end if

          initpressure(i,j) = compute_psi_incident(R(i),Z(j),switch)  !! the incident wave
          
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
       do k = 1, num_sup_zeros
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

    allocate(supinspressure(num_sup_zeros,Nmeshr,Nphi))

    print*,''
    print*,'Also, computing polar near-field for num_sup_zeros > 0:'
    print*,''

    do j = 1, Nphi  !! the polar sweep anle
       do i = 1, Nmeshr  !! the radial (vertical in a 2D plot) direction

          print*,'Pressure at phi=', phi(j)*180._dpk/PI,'and R=',R(i)
          
          do k = 1, num_sup_zeros
             supinspressure(k,i,j) = residueprpolar(R(i),phi(j),k)
          end do
          
          print*,''
          
       end do
    end do

!! dump data
    
    do k = 1, num_sup_zeros
       write(supind,"(I1)"),k
       open(1,file='supinstabilitypr_polar.out.'//supind,form='UNFORMATTED')
       write(1) Nphi,Nmeshr,1
       write(1) ((REAL(supinspressure(k,i,j)),j=1,Nphi),i=1,Nmeshr)
       close(1)  
    end do

    do k = 1, num_sup_zeros
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

    complex(dpk)          :: intgrl_A1_at_mu_plus, intgrl_A1_at_KH_zero_1, intgrl_A1_at_KH_zero_2, intgrl_at_sup_zero
    complex(dpk)          :: k_plus_at_mu_plus, f1, f2, f3
    integer               :: f, i1

  
    if ((vortswitch == 1) .OR. (vortswitch == 2)) then
    !! the factor \Psi_{mn}(1) of (4.1) [see the JFM 2008]:

	       psi = omega_r*resp*(1._dpk - M2*mu_plus)*compute_Trs_instab(1._dpk,mu_plus,2)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*(-Zo))

    else

!! check if the radial wavenumber is negative:

	       if((REAL(mu_plus) < -(kapT/(1._dpk - kapT*M2))) .OR. (REAL(mu_plus) > (1._dpk/(1._dpk + M1)))) then 
		          print*,'Radial wave number(s) is negative: Check the incident wave axial wave number!'
		          STOP
	       end if

	       alpha1 = omega_r*SQRT((1._dpk - mu_plus*M1)**2 - mu_plus**2)
	       alpha2 = omega_r*SQRT(kapT**2*(1._dpk - mu_plus*M2)**2 - mu_plus**2)

	       print*,'precompute: The radial wave numbers:'
	       write(*,'(/A12,2X,2F15.10)'), ' alpha1:->', alpha1
	       write(*,'(A12,2X,2F15.10/)'), ' alpha2:->', alpha2

!! the factor \Psi_{mn}(1) of (3.30) [see the JFM]:

	       f1 = (1._dpk - mu_plus*M1)/(1._dpk - mu_plus*M2)*bessj(alpha1*h,azim_mode,1)*EXP(ABS(AIMAG(alpha1*h)))
	       
	       f2 = (bessj(alpha2,azim_mode,1)*dhank1(alpha2,azim_mode,1)- & 
		         hank1(alpha2,azim_mode,1)*dbessj(alpha2,azim_mode,1))* &
		         EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*alpha2 + ABS(AIMAG(alpha2)))
	       
	       f3 = bessj(alpha2*h,azim_mode,1)*dhank1(alpha2,azim_mode,1)* &
		         EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*alpha2 + ABS(AIMAG(alpha2*h)))- & 
		         hank1(alpha2*h,azim_mode,1)*dbessj(alpha2,azim_mode,1)* &
		         EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*alpha2*h + ABS(AIMAG(alpha2)))

	       psi = kap_rho*f1*f2/f3
	       print*, 'precompute: Evaluated psi for the incident wave'
	   
    end if
    
!!  the factor Kt^{-}(\mu_{mn}^{+}):

    if (vortswitch .EQ. 0) then
	     
              if (REAL(mu_plus) < cont_cross_over_pt) then  !! mu_plus is below the contour; cont_cross_over_pt being the crossover pt
		           print*, ''
		           print*, 'precompute: The incident acoustic mode needs to be INSIDE the contour'
		           print*, ''
		           STOP
	           end if

     end if
     
    print*, 'precompute: Evaluating kernel at mu_plus'
    call compute_eqn_A1_integral(mu_plus,intgrl_A1_at_mu_plus,0,0,1)
       
    k_plus_at_mu_plus = EXP(-intgrl_A1_at_mu_plus/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) + & 
                        LOG(compute_kernel(0,mu_plus)/compute_U_s_factor(0,mu_plus)))

    k_minus_at_mu_plus =  compute_kernel(0,mu_plus)/(k_plus_at_mu_plus*compute_U_s_factor(0,mu_plus)) 

!!  the factor Kt^{+}(s_{z1}):

    if (vortswitch .EQ. 0) then
  
	       print*, 'precompute: Evaluating Kt^{+}(s_{z1}) at KH_zero_1'
	       call compute_eqn_A1_integral(KH_zero_1,intgrl_A1_at_KH_zero_1,0,0,1)

	       k_plus_sz1 = EXP(-intgrl_A1_at_KH_zero_1/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! NOTE: zero KH_zero_1 has to lie below
									  !! the contour
    end if
!!  the factor K^{+}(s_{z2}):

    print*, 'precompute: Evaluating Kt^{+}(s_{z2}) at KH_zero_2'
    call compute_eqn_A1_integral(KH_zero_2,intgrl_A1_at_KH_zero_2,0,0,1)

    k_plus_sz2 = EXP(-intgrl_A1_at_KH_zero_2/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! same for KH_zero_2

    if (num_sup_zeros > 0) allocate(kzsp(num_sup_zeros))

    if (vortswitch .EQ. 0) then
	       do i1 = 1, num_sup_zeros
		  
		       print '(A, I0)', 'precompute: Evaluating kernel at supersonic zero ', i1
		       call compute_eqn_A1_integral(sup_zeros_list(i1),intgrl_at_sup_zero,0,0,1)
		  
		       kzsp(i1) = EXP(-intgrl_at_sup_zero/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! NOTE: instability zeros are below
		  
	       end do
    end if


  END SUBROUTINE precompute


  SUBROUTINE compute_eqn_A1_integral(s_target,integral_value,kswitch,ch,KU_K_switch)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Evaluate (A1)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)           :: s_target
    complex(dpk)           :: integral_value ! final value of integration
    integer                :: num_of_quad_points, ch
    integer                :: KU_K_switch  !! 1: K/U; else: K
    integer                :: kswitch  !! 1: compute (4.10); else: compute (3.22)

    allocate(intpanel(tot_ker_points-1))
    allocate(Npanel(tot_ker_points-1))

    call adaptive(s_target,kswitch,ch,KU_K_switch)

    call sum_panel_contributions_kernel(integral_value,num_of_quad_points,s_target)

    deallocate(intpanel)
    deallocate(Npanel)


  END SUBROUTINE compute_eqn_A1_integral


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


    do j = 1, tot_ker_points-1  !! the integration points already specified

       len = ker_int_points(j+1) - ker_int_points(j)  !! the length of a mini panel

!! compute "scale" needed to ensure the tolerance check:

       if (j >= 1 .AND. j < num_ker_pts_loop+2 ) then
          scale = len/(def_pts_ker_cntr(3) - def_pts_ker_cntr(1))

       else
          scale = len/(def_pts_ker_cntr(5) - def_pts_ker_cntr(3))

       end if

       call kernel_trapz_int(si,ker_int_points(j),ker_int_points(j+1),ksw,sw,T)  !! the basis for comparison

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

          xp(1) = REAL(ker_int_points(j))
          xp(panel_no+1) = REAL(ker_int_points(j+1))
          yp(1) = AIMAG(ker_int_points(j))
          yp(panel_no+1) = AIMAG(ker_int_points(j+1))
          zp(1) = ker_int_points(j)
          zp(panel_no+1) = ker_int_points(j+1)

          do i = 1,panel_no-1
      
             if (j >= 1 .AND. j < num_ker_pts_loop+2) then
                xp(i+1) = xp(i) + ds
                call get_y_alg_int_contour(xp(i+1),yp(i+1), & 
                           REAL(def_pts_ker_cntr(3)),REAL(def_pts_ker_cntr(2)-def_pts_ker_cntr(3)), &
                           AIMAG(def_pts_ker_cntr(2)-def_pts_ker_cntr(3)))
             else
                xp(i+1) = xp(i) + ds
                call get_y_alg_int_contour(xp(i+1),yp(i+1), & 
                           REAL(def_pts_ker_cntr(3)),REAL(def_pts_ker_cntr(4)-def_pts_ker_cntr(3)), &
                           AIMAG(def_pts_ker_cntr(4)-def_pts_ker_cntr(3)))
                  
             end if

             zp(i+1) = CMPLX(xp(i+1),yp(i+1))
             
          end do

!          print*, 'zp(i):', zp(i), 'zp(i+1):', zp(i+1)
          
          allocate(T_temp(panel_no))
          
          do i = 1,panel_no
             call kernel_trapz_int(si,zp(i),zp(i+1),ksw,sw,T_temp(i))
          end do
          
          intpanel(j) = (0._dpk,0._dpk)

          do i = 1,panel_no
             intpanel(j) = intpanel(j) + T_temp(i)
          end do

          deallocate(xp)
          deallocate(yp)
          deallocate(T_temp)

!          print*, 'ker_int_points:', ker_int_points(j)
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


  SUBROUTINE compute_fplus(switch,index)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Evaluate (3.30)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer          :: i, index, switch

    allocate(fplusz(tot_IFT_pts))  !! fplus is same as \xi^{+}(s) of (3.30)
    fplusz = (0._dpk,0._dpk)

    if (switch == 0) then  !! switch = 0 :> fresh job; no restart file read before

       do i = 1,tot_IFT_pts
          print *, 'compute_fplus: F+ at index: ', i, '  point = ', iftpoints(i)
          fplusz(i) = get_fplus_value(iftpoints(i))
          
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

       if (index .NE. tot_IFT_pts) then

          do i = index+1,tot_IFT_pts
             print*, 'F+ at:', iftpoints(i)
             fplusz(i) = get_fplus_value(iftpoints(i))
          
             write(*,'(A22,2X,2F20.10/)'),'F+:->',fplusz(i)

             open(10,file='fplus_part.out',form='FORMATTED',position='APPEND')
             write(10,'(I5,4E20.10)') i,iftpoints(i),fplusz(i)
             close(10)

          end do

       end if

    end if

    open(10,file='fplus.out',form='FORMATTED')
    do i = 1, tot_IFT_pts
       write(10,'(I5,4E20.10)') i,iftpoints(i),fplusz(i)
    end do
    close(10)
    

  END SUBROUTINE compute_fplus


  SUBROUTINE read_fplus(j)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Use the computed F^{+} in case running a restart job
!! 2. Note that restarting is only possible during the "computing F+" portion
!! 3. If a job exits during (R,Z) computation, a "restart" just starts the job back at (R=0,Z=0)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer          :: i, j
    real             :: dummy1, dummy2

    allocate(fplusz_temp(tot_IFT_pts))

    fplusz_temp = (0._dpk,0._dpk)

    open(10,file='fplus_part.out',form='FORMATTED')
    do i = 1, tot_IFT_pts
       read(10,'(I5,4E20.10)',end = 100) j,dummy1,dummy2,fplusz_temp(i)
    end do
    close(10)

100 do i = 1,j
       print*, 'F+ at:', iftpoints(i)

       write(*,'(A22,2X,2F20.10/)'),'F+:->',fplusz_temp(i)
    end do

    close(10)

  END SUBROUTINE read_fplus


  FUNCTION get_fplus_value(s_target)  result(fplus)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Actual computation of \xi^{+}(s)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)       :: fplus, s_target, gpz, kpz, int_A1_at_s_target
    complex(dpk)       :: lambda_1_minus_mu_plus, lambda_2_minus_mu_plus, lambda_3_minus_mu_plus
    complex(dpk)       :: lambda_1_plus_mu_plus, lambda_2_plus_mu_plus, lambda_3_plus_mu_plus, lambda_prod
    integer            :: f


    lambda_1_minus_mu_plus = sqrt(1._dpk - mu_plus*(M1 - 1._dpk))
    lambda_2_minus_mu_plus = sqrt(kapT - mu_plus*(kapT*M2 - 1._dpk))
    lambda_3_minus_mu_plus = sqrt(kapT - mu_plus*(kapT*M3 - 1._dpk))
    lambda_1_plus_mu_plus = sqrt(1._dpk - s_target*(M1 + 1._dpk))
    lambda_2_plus_mu_plus = sqrt(kapT - s_target*(kapT*M2 + 1._dpk))
    lambda_3_plus_mu_plus = sqrt(kapT - s_target*(kapT*M3 + 1._dpk))

    gpz = psi*(1._dpk - mu_plus*M2)/(k_minus_at_mu_plus*(mu_plus - KH_zero_2))

    call compute_eqn_A1_integral(s_target,int_A1_at_s_target,0,0,1)

    if ((farswitch == 1) .OR. (farswitch == 2)) then
       if (REAL(s_target) >= cont_cross_over_pt) then
          kpz = EXP(-int_A1_at_s_target/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)) + & 
               LOG(compute_kernel(0,s_target)/compute_U_s_factor(0,s_target)))  !! s_target is above contour; cont_cross_over_pt is crossover pt
       else
          kpz = EXP(-int_A1_at_s_target/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! s_target is below
       end if
    else
       kpz = EXP(-int_A1_at_s_target/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))  !! s_target is always below normally
    end if

    fplus = gpz*( (s_target - KH_zero_2)/(mu_plus - s_target) + vs_param_gamma)/(kpz*compute_U_s_factor(0,s_target))

  END FUNCTION get_fplus_value


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


  SUBROUTINE get_y_alg_int_contour(x,y,p,a,b)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Algebraic "one piece" contour:  see Rienstra, JEM, 2007
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)    :: x, y
    real(dpk)    :: a, b, p

    y = b*(4._dpk*(x-p)/a)/(3._dpk + ((x-p)/a)**4) ! algebraic


  END SUBROUTINE get_y_alg_int_contour


  SUBROUTINE kernel_trapz_int(s,z1,z2,ksw,sw,int_value)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The trapezoidal rule for the kernel contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)                :: z1, z2, s, dz, int_value
    complex(dpk)                :: f1, f2
    integer                     :: sw, ksw

    f1 = kernel_integrand(ksw,sw,s,z1)
    f2 = kernel_integrand(ksw,sw,s,z2)

    dz = (z2-z1)

    int_value = dz/2._dpk*(f1+f2)


  END SUBROUTINE kernel_trapz_int


  SUBROUTINE IFT_trapz_int(r,z,i,ss,Int)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The trapezoidal rule for the IFT contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                   :: r, z
    complex(dpk)                :: ds, Int
    complex(dpk)                :: f1, f2
    integer                     :: i, ss

    if (prswitch == 0) then

       f1 = integrand_IFT_pot(r,z,i,ss)
       f2 = integrand_IFT_pot(r,z,i+1,ss)

    else

       f1 = integrandiftpr(r,z,i,ss)
       f2 = integrandiftpr(r,z,i+1,ss)

    end if

    ds = iftpoints(i+1) - iftpoints(i)

    Int = ds/2._dpk*(f1+f2)


  END SUBROUTINE IFT_trapz_int


  SUBROUTINE sum_panel_contributions_kernel(Int,totpoints,z)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Sum up the individual trapezoids for the kernel contour
!! 2. Dump data
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)                               :: Int, z
    integer                                    :: totpoints
    integer                                    :: i
    character(60)                              :: arg

    Int = (0._dpk,0._dpk)
    do i = 1,tot_ker_points-1
       Int = Int + intpanel(i)  !! summing up the individual panels
    end do

!! dump data:

    write(arg,"('u=',E12.5,'+',E12.5,'i')"),REAL(z),AIMAG(z)

    open(10,file='DataDump/Kernel/intkernel.'//arg,form='FORMATTED')
    do i = 1, tot_ker_points-1
       write(10,'(I5,2E20.10)') i,intpanel(i)
    end do
    close(10)
    
    write(*,'(/A22,2X,2F20.10/)'),'Integral value:->', Int

    totpoints = 1
    do i = 1,tot_ker_points-1
       totpoints = totpoints + Npanel(i)  !! total quadrature points
    end do

    write(*,'(A22,2X,I20/)'),'Quadrature pts:->', totpoints
    

  END SUBROUTINE sum_panel_contributions_kernel


 SUBROUTINE sum_panel_contribution_IFT(panel,Int)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Sum up the individual trapezoids for the IFT contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)                               :: Int
    complex(dpk),dimension(tot_IFT_pts-1)        :: panel
    integer                                    :: i

    Int = (0._dpk,0._dpk)
    do i = 1,tot_IFT_pts-1
       Int = Int + panel(i)
    end do


  END SUBROUTINE sum_panel_contribution_IFT


  FUNCTION kernel_integrand(kswitch,switch,si,zi) result(integrand)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the kernel integrand
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: si, zi, integrand
    complex(dpk)    :: In, Id
    integer         :: switch, kswitch

    if (switch == 1)  then !! K/U
       In = LOG(compute_kernel(kswitch,zi)/compute_U_s_factor(kswitch,zi))
    else
       In = LOG(compute_kernel(kswitch,zi))
    end if
    Id = zi-si
    integrand = In/Id
    
!!$    if (ABS(integrand) > 1.0E3) then
!!$       print*,'zi:',zi
!!$!       print*,'integrand:',integrand
!!$!       print*,'kernel:',compute_kernel(zi)
!!$       print*,'In:',In
!!$       print*,'Id:',Id
!!$    end if

    if (compute_kernel(kswitch,zi) == 0.) then 
       print*,''
       print*,'FATAL ERROR: The kernel has encountered a zero (at zi)!!'
       print*,'si:-->',si
       print*,'zi:-->',zi
       STOP
    end if


  END FUNCTION kernel_integrand


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

       integrandiftpr = (1._dpk - u*M2)**2*compute_Trs(ri,u,1)*fplusz(i)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!!$       integrandiftpr = compute_Trs(ri,u,1)
    else if (ss==2) then

!! pressure:
       
       integrandiftpr = (1._dpk - u*M2)**2*compute_Trs(ri,u,2)*fplusz(i)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!!$       integrandiftpr = compute_Trs(ri,u,2)
    else

!! pressure:

       integrandiftpr = (1._dpk - u*M3)**2*compute_Trs(ri,u,3)*fplusz(i)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)
!!$       integrandiftpr =  compute_Trs(ri,u,3)
    end if


  END FUNCTION integrandiftpr


  FUNCTION integrand_IFT_pot(ri,zi,i,ss) result(integrandiftpot)

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

       integrandiftpot = (1._dpk - u*M2)**2/(1._dpk - u*M1)*compute_Trs(ri,u,1)*fplusz(i)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)

    else if (ss==2) then

!! potential:
       
       integrandiftpot = (1._dpk - u*M2)*compute_Trs(ri,u,2)*fplusz(i)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)

    else

!! potential:

       integrandiftpot = (1._dpk - u*M3)*compute_Trs(ri,u,3)*fplusz(i)* &
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*u*zi)

    end if


  END FUNCTION integrand_IFT_pot


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

          F2n = dhank1(lambda_2*omega_r,azim_mode,1)*bessj(lambda_2*omega_r*h,azim_mode,1)*EXP(ABS(AIMAG(lambda_2*omega_r*h))+ &
               CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r) - dbessj(lambda_2*omega_r,azim_mode,1)*hank1(lambda_2*omega_r*h,azim_mode,1)* & 
               EXP(ABS(AIMAG(lambda_2*omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r*h)
          F2d = dhank1(lambda_2*omega_r,azim_mode,1)*dbessj(lambda_2*omega_r*h,azim_mode,1)*EXP(ABS(AIMAG(lambda_2*omega_r*h))+ &
               CMPLX(0._dpk,1._dpk,kind=dpk)*lambda_2*omega_r) - dbessj(lambda_2*omega_r,azim_mode,1)*dhank1(lambda_2*omega_r*h,azim_mode,1)* &
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


  SUBROUTINE kernelcheck(N,ss,z,I)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Check so that the kernel does not cross the negative real axis on C
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), dimension(N) :: z, I
    real(dpk)                  :: x1, x2, y1, y2, intercept
    integer                    :: j, N, ss

    do j = 1, N
       I(j) = compute_kernel(ss,z(j))/compute_U_s_factor(ss,z(j))
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

    write(arg,"('w=',E12.5,'Np=',I7.7)"),REAL(omega_r), N

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

    do i = 1, num_zeros_s1_s2

       bc = bc/(z - zeros_list_bw_s1_s2(i))

    end do

    do i = 1, num_poles_s1_s2

       bc = bc*(z - poles_list_bw_s1_s2(i))

    end do


  END FUNCTION bc


  FUNCTION compute_Trs(ri,si,ss) result(Trs)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute Trs of (3.33)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: ri
    complex(dpk)  :: sri, Trs
    complex(dpk)  :: si
    complex(dpk)  :: l1, l2, l3, lz
    complex(dpk)  :: Fn, Fd, F, Rn, Rd, Rz
    integer       :: ss


    l1 = sqrt(1._dpk - si*(M1+1._dpk))*sqrt(1._dpk - si*(M1-1._dpk))
    l2 = sqrt(kapT - si*(kapT*M2+1._dpk))*sqrt(kapT - si*(kapT*M2-1._dpk))
    l3 = sqrt(kapT - si*(kapT*M3+1._dpk))*sqrt(kapT - si*(kapT*M3-1._dpk))

    lz = kap_rho*l2/l1*(1._dpk - si*M1)**2/(1._dpk - si*M2)**2


    if (ss==1) then
       
       Rd = (lz*bessj(l1*omega_r*h,azim_mode,1)*dhank1(l2*omega_r*h,azim_mode,1)- &
            hank1(l2*omega_r*h,azim_mode,1)*dbessj(l1*omega_r*h,azim_mode,1))* &
            EXP(ABS(AIMAG(l1*omega_r*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*h)
       
       Rn = (bessj(l2*omega_r*h,azim_mode,1)*dbessj(l1*omega_r*h,azim_mode,1)- &
            lz*bessj(l1*omega_r*h,azim_mode,1)*dbessj(l2*omega_r*h,azim_mode,1))* &
            EXP(ABS(AIMAG(l1*omega_r*h))+ABS(AIMAG(l2*omega_r*h)))
       
       Rz = Rn/Rd
          
       Fn = bessj(l1*omega_r*ri,azim_mode,1)*EXP(ABS(AIMAG(l1*omega_r*ri)))*(Rz*hank1(l2*omega_r*h,azim_mode,1)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*h) + bessj(l2*omega_r*h,azim_mode,1)* &
            EXP(ABS(AIMAG(l2*omega_r*h))))
       
       Fd = bessj(l1*omega_r*h,azim_mode,1)*EXP(ABS(AIMAG(l1*omega_r*h)))*(Rz*dhank1(l2*omega_r,azim_mode,1)* & 
            EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r) + dbessj(l2*omega_r,azim_mode,1)* &
            EXP(ABS(AIMAG(l2*omega_r))))
       
       F = Fn/Fd

       Trs = F/l2

    else if (ss==2) then

       Rd = (lz*bessj(l1*omega_r*h,azim_mode,1)*dhank1(l2*omega_r*h,azim_mode,1)- &
            hank1(l2*omega_r*h,azim_mode,1)*dbessj(l1*omega_r*h,azim_mode,1))* &
            EXP(ABS(AIMAG(l1*omega_r*h))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*h)
       
       Rn = (bessj(l2*omega_r*h,azim_mode,1)*dbessj(l1*omega_r*h,azim_mode,1)- &
            lz*bessj(l1*omega_r*h,azim_mode,1)*dbessj(l2*omega_r*h,azim_mode,1))* &
            EXP(ABS(AIMAG(l1*omega_r*h))+ABS(AIMAG(l2*omega_r*h)))
       
       Rz = Rn/Rd
          
       Fn = Rz*hank1(l2*omega_r*ri,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*ri)+ &
            bessj(l2*omega_r*ri,azim_mode,1)*EXP(ABS(AIMAG(l2*omega_r*ri)))
       
       Fd = Rz*dhank1(l2*omega_r,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r)+ &
            dbessj(l2*omega_r,azim_mode,1)*EXP(ABS(AIMAG(l2*omega_r)))
       
       F = Fn/Fd

       Trs = F/l2

    else

!!$       if ((ABS(l3*omega_r) < asymplim .AND. ABS(AIMAG(l3*omega_r)) < asymplim1) .AND. &
!!$            (ABS(l3*omega_r*ri) < asymplim .AND. ABS(AIMAG(l3*omega_r*ri)) < asymplim1)) then

          Fn = hank1(l3*omega_r*ri,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l3*omega_r*ri)
       
          Fd = dhank1(l3*omega_r,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*l3*omega_r)
       
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


  FUNCTION compute_Trs_instab(ri,si,ss)  result(Trsin)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute Trs of (4.8)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: ri
    complex(dpk)  :: sri, Trsin
    complex(dpk)  :: si
    complex(dpk)  :: l1, l2
    complex(dpk)  :: Fn, Fd, F
    integer       :: ss


    l1 = sqrt(1._dpk - si*(M1+1._dpk))*sqrt(1._dpk - si*(M1-1._dpk))
    l2 = sqrt(1._dpk - si*(M2+1._dpk))*sqrt(1._dpk - si*(M2-1._dpk))


    if (ss==1) then
       
       Fn = bessj(l1*omega_r*ri,azim_mode,1)*EXP(ABS(AIMAG(l1*omega_r*ri)))
       Fd = dbessj(l1*omega_r*h,azim_mode,1)*EXP(ABS(AIMAG(l1*omega_r*h)))
       
       F = Fn/Fd

       Trsin = F/l1

    else

       Fn = dhank1(l2*omega_r,azim_mode,1)*bessj(l2*omega_r*ri,azim_mode,1)*EXP(ABS(AIMAG(l2*omega_r*ri))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r) - dbessj(l2*omega_r,azim_mode,1)*hank1(l2*omega_r*ri,azim_mode,1)* &
            EXP(ABS(AIMAG(l2*omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*ri)

       Fd = dhank1(l2*omega_r,azim_mode,1)*dbessj(l2*omega_r*h,azim_mode,1)*EXP(ABS(AIMAG(l2*omega_r*h))+ &
            CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r) - dbessj(l2*omega_r,azim_mode,1)*dhank1(l2*omega_r*h,azim_mode,1)* &
            EXP(ABS(AIMAG(l2*omega_r))+CMPLX(0._dpk,1._dpk,kind=dpk)*l2*omega_r*h)
       
       F = Fn/Fd

       Trsin = F/l2

    end if
       
!    print*,'Rz:',Rz


  END FUNCTION compute_Trs_instab


  FUNCTION compute_psi_incident(r,z,ss) result(psi0)

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
		     psimn = omega_r*resp*(1._dpk - M1*mu_plus)*compute_Trs_instab(r,mu_plus,1)
		     psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*(1._dpk - M1*mu_plus)* &
			  psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*(z-Zo))

		  else if (ss == 2) then
	       
		     psimn = omega_r*resp*(1._dpk - M2*mu_plus)*compute_Trs_instab(r,mu_plus,2)
		     psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*(1._dpk - M2*mu_plus)* &
			  psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*(z-Zo))

		  else

		     psi0 = 0.

		  end if

	       else
		  psi0 = 0.
	       end if

    else

	       if (ss == 1) then
		  psimn = bessj(alpha1*r,azim_mode,1)*EXP(ABS(AIMAG(alpha1*r)))
	       
		  psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*(1._dpk - M1*mu_plus)* &
		       psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*z)


	       else if (ss == 2) then

		  f1 = (1._dpk - mu_plus*M1)/(1._dpk - mu_plus*M2)*bessj(alpha1*h,azim_mode,1)
	       
		  f2 = bessj(alpha2*r,azim_mode,1)* &
		       dhank1(alpha2,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* & 
		       alpha2 + ABS(AIMAG(alpha2*r)))- & 
		       hank1(alpha2*r,azim_mode,1)* &
		       dbessj(alpha2,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* & 
		       alpha2*r + ABS(AIMAG(alpha2)))
		  

		  f3 = bessj(alpha2*h,azim_mode,1)* &
		       dhank1(alpha2,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
		       alpha2 + ABS(AIMAG(alpha2*h)))- & 
		       hank1(alpha2*h,azim_mode,1)* &
		       dbessj(alpha2,azim_mode,1)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)* &
		       alpha2*h + ABS(AIMAG(alpha2)))
		    

		  psimn = kap_rho*f1*f2/f3

		  psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*(1._dpk - M2*mu_plus)* &
		       psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*mu_plus*z)

	!!$       open(10,file='test.out',form='FORMATTED',position='APPEND')
	!!$       write(10,'(2E20.10)') psimn
	!!$       close(10)

	       else

		  psi0 = 0.

	       end if

    end if


  END FUNCTION compute_psi_incident



  SUBROUTINE checkloc(z,switch)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Subroutine to check the location of the instability zeros and pole wrt to the contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: z
    real(dpk)       :: zx, zy, zfy
    integer         :: switch

    zx = REAL(z)
    zy = AIMAG(z)

    call get_y_alg_int_contour(zx,zfy,REAL(def_pts_IFT_cntr(3)),REAL(def_pts_IFT_cntr(4)-def_pts_IFT_cntr(3)),AIMAG(def_pts_IFT_cntr(4)-def_pts_IFT_cntr(3)))
!!$    zfy = AIMAG(c1)

    if(zy > zfy) then
       switch = 0  !! inside
    else
       switch = 1  !! outside
    end if


  END SUBROUTINE checkloc


  FUNCTION residueprpolar(r,phi,ss)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the residue pressure term of (3.38)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)            :: r, phi
    integer              :: ss, ii, jj
    complex(dpk)         :: residueprpolar, res, fn
   
    res = psi*(1._dpk - mu_plus*M2)/((mu_plus - sup_zeros_list(ss))* &
         k_minus_at_mu_plus*kzsp(ss)*(sup_zeros_list(ss) - KH_zero_1)*(sup_zeros_list(ss) - KH_zero_2))
!!$    res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss) - KH_pole_1)/((mu_plus - sup_zeros_list(ss))* &
!!$         k_minus_at_mu_plus*kzsp(ss)*(sup_zeros_list(ss) - KH_zero_1)*(sup_zeros_list(ss) - KH_zero_2))
!!$    res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss)-KH_pole_1)*(sup_zeros_list(ss)-CONJG(KH_pole_1))/ &
!!$            ((mu_plus - sup_zeros_list(ss))*k_minus_at_mu_plus*kzsp(ss)*(sup_zeros_list(ss)-CONJG(KH_zero_1))* &
!!$            (sup_zeros_list(ss)-KH_zero_1)*(sup_zeros_list(ss) - KH_zero_2)*(sup_zeros_list(ss)-CONJG(KH_zero_2)))
    do ii = 1, num_sup_zeros
       if (ii .NE. (ss)) then
          res = res/(sup_zeros_list(ss) - sup_zeros_list(ii))
       end if
    end do
    do jj = 1, num_sup_poles
       res = res*(sup_zeros_list(ss) - sup_poles_list(jj))
    end do
    
    if (r*SIN(phi) <= h) then
       fn = (1._dpk - sup_zeros_list(ss)*M2)**2*compute_Trs(r*SIN(phi),sup_zeros_list(ss),1)
    else if (r*SIN(phi) <= 1. .AND. r*SIN(phi) > h) then
       fn = (1._dpk - sup_zeros_list(ss)*M2)**2*compute_Trs(r*SIN(phi),sup_zeros_list(ss),2)
    else
       fn = (1._dpk - sup_zeros_list(ss)*M3)**2*compute_Trs(r*SIN(phi),sup_zeros_list(ss),3)  
    end if

    residueprpolar = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
         omega_r*sup_zeros_list(ss)*r*COS(phi))*res*fn

    
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
!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_1-KH_pole_1)*(KH_zero_1-CONJG(KH_pole_1))/ &
!!$            ((mu_plus - KH_zero_1)*k_minus_at_mu_plus*kpKH_zero_1*(KH_zero_1-CONJG(KH_zero_1))*(KH_zero_1-KH_zero_2)*(KH_zero_1-CONJG(KH_zero_2)))
!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_1 - KH_pole_1)/((mu_plus - KH_zero_1)*kmmu_plus*k_plus_sz1*(KH_zero_1 - KH_zero_2))
       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - KH_zero_1)*k_minus_at_mu_plus*k_plus_sz1*(KH_zero_1 - KH_zero_2))
       if (vortswitch .EQ. 0) then
          do ii = 1, num_sup_zeros
             res = res/(KH_zero_1 - sup_zeros_list(ii))
          end do
          do jj = 1, num_sup_poles
             res = res*(KH_zero_1 - sup_poles_list(jj))
          end do
       end if

       if (r <= h) then
          fn = (1._dpk - KH_zero_1*M2)**2*compute_Trs(r,KH_zero_1,1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - KH_zero_1*M2)**2*compute_Trs(r,KH_zero_1,2)
       else
          fn = (1._dpk - KH_zero_1*M3)**2*compute_Trs(r,KH_zero_1,3)
       end if
       
       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            omega_r*KH_zero_1*z)*res*fn
          
!!$       else
!!$          residuepr = 0.
!!$
!!$       end if
          
    elseif (ss == 2) then

!!$       if (z > 0.) then
!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_2-KH_pole_1)*(KH_zero_2-CONJG(KH_pole_1))/ &
!!$            ((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2-CONJG(KH_zero_1))*(KH_zero_2-KH_zero_1)*(KH_zero_2-CONJG(KH_zero_2)))
!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_2 - KH_pole_1)/((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2 - KH_zero_1))
       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2 - KH_zero_1))
       if (vortswitch .EQ. 0) then
          do ii = 1, num_sup_zeros
             res = res/(KH_zero_2 - sup_zeros_list(ii))
          end do
          do jj = 1, num_sup_poles
             res = res*(KH_zero_2 - sup_poles_list(jj))
          end do
       end if
          
       if (r <= h) then
          fn = (1._dpk - KH_zero_2*M2)**2*compute_Trs(r,KH_zero_2,1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - KH_zero_2*M2)**2*compute_Trs(r,KH_zero_2,2)
       else
          fn = (1._dpk - KH_zero_2*M3)**2*compute_Trs(r,KH_zero_2,3)
       end if

       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            omega_r*KH_zero_2*z)*res*fn

!!$       else
!!$          residuepr = 0.
!!$
!!$       end if

    else

!!$       if (z > 0.) then
!!$       res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss-2)-KH_pole_1)*(sup_zeros_list(ss-2)-CONJG(KH_pole_1))/ &
!!$            ((mu_plus - sup_zeros_list(ss-2))*k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2)-CONJG(KH_zero_1))* &
!!$            (sup_zeros_list(ss-2)-KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2)*(sup_zeros_list(ss-2)-CONJG(KH_zero_2)))
!!$       res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss-2) - KH_pole_1)/((mu_plus - sup_zeros_list(ss-2))* &
!!$            k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2) - KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2))
       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - sup_zeros_list(ss-2))* &
            k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2) - KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2))
       do ii = 1, num_sup_zeros
          if (ii .NE. (ss-2)) then
             res = res/(sup_zeros_list(ss-2) - sup_zeros_list(ii))
          end if
       end do
       do jj = 1, num_sup_poles
          res = res*(sup_zeros_list(ss-2) - sup_poles_list(jj))
       end do
          
       if (r <= h) then
          fn = (1._dpk - sup_zeros_list(ss-2)*M2)**2*compute_Trs(r,sup_zeros_list(ss-2),1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - sup_zeros_list(ss-2)*M2)**2*compute_Trs(r,sup_zeros_list(ss-2),2)
       else
          fn = (1._dpk - sup_zeros_list(ss-2)*M3)**2*compute_Trs(r,sup_zeros_list(ss-2),3)
       end if

       residuepr = CMPLX(0._dpk,1._dpk,kind=dpk)*omega_r*omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)* &
            omega_r*sup_zeros_list(ss-2)*z)*res*fn

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
!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_1-KH_pole_1)*(KH_zero_1-CONJG(KH_pole_1))/ &
!!$            ((mu_plus - KH_zero_1)*k_minus_at_mu_plus*k_plus_sz1*(KH_zero_1-CONJG(KH_zero_1))*(KH_zero_1-KH_zero_2)*(KH_zero_1-CONJG(KH_zero_2)))
!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_1 - KH_pole_1)/((mu_plus - KH_zero_1)*k_minus_at_mu_plus*k_plus_sz1*(KH_zero_1 - KH_zero_2))
       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - KH_zero_1)*k_minus_at_mu_plus*k_plus_sz1*(KH_zero_1 - KH_zero_2))
       if (vortswitch .EQ. 0) then
          do ii = 1, num_sup_zeros
             res = res/(KH_zero_1 - sup_zeros_list(ii))
          end do
          do jj = 1, num_sup_poles
             res = res*(KH_zero_1 - sup_poles_list(jj))
          end do
       end if

       if (r <= h) then
          fn = (1._dpk - KH_zero_1*M2)**2/(1._dpk - KH_zero_1*M1)*compute_Trs(r,KH_zero_1,1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - KH_zero_1*M2)*compute_Trs(r,KH_zero_1,2)
       else
          fn = (1._dpk - KH_zero_1*M3)*compute_Trs(r,KH_zero_1,3)
       end if
       
       residuepot = omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*omega_r*KH_zero_1*z)*res*fn
          
!!$       else
!!$          residuepot = 0.
!!$
!!$       end if
          
    elseif (ss == 2) then

!!$       if (z > 0.) then
!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_2-KH_pole_1)*(KH_zero_2-CONJG(KH_pole_1))/ &
!!$            ((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2-CONJG(KH_zero_1))*(KH_zero_2-KH_zero_1)*(KH_zero_2-CONJG(KH_zero_2)))
!!$       res = psi*(1._dpk - mu_plus*M2)*(KH_zero_2 - KH_pole_1)/((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2 - KH_zero_1))
       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - KH_zero_2)*k_minus_at_mu_plus*k_plus_sz2*(KH_zero_2 - KH_zero_1))
       if (vortswitch .EQ. 0) then
          do ii = 1, num_sup_zeros
             res = res/(KH_zero_1 - sup_zeros_list(ii))
          end do
          do jj = 1, num_sup_poles
             res = res*(KH_zero_1 - sup_poles_list(jj))
          end do
       end if
          
       if (r <= h) then
          fn = (1._dpk - KH_zero_2*M2)**2/(1._dpk - KH_zero_2*M1)*compute_Trs(r,KH_zero_2,1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - KH_zero_2*M2)*compute_Trs(r,KH_zero_2,2)
       else
          fn = (1._dpk - KH_zero_2*M3)*compute_Trs(r,KH_zero_2,3)
       end if

       residuepot = omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*omega_r*KH_zero_2*z)*res*fn

!!$       else
!!$          residuepot = 0.
!!$
!!$       end if

    else

!!$       if (z > 0.) then
!!$      res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss-2)-KH_pole_1)*(sup_zeros_list(ss-2)-CONJG(KH_pole_1))/ &
!!$            ((mu_plus - sup_zeros_list(ss-2))*k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2)-CONJG(KH_zero_1))* &
!!$            (sup_zeros_list(ss-2)-KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2)*(sup_zeros_list(ss-2)-CONJG(KH_zero_2)))
!!$       res = psi*(1._dpk - mu_plus*M2)*(sup_zeros_list(ss-2) - KH_pole_1)/((mu_plus - sup_zeros_list(ss-2))* &
!!$            k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2) - KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2))
       res = psi*(1._dpk - mu_plus*M2)/((mu_plus - sup_zeros_list(ss-2))* &
            k_minus_at_mu_plus*kzsp(ss-2)*(sup_zeros_list(ss-2) - KH_zero_1)*(sup_zeros_list(ss-2) - KH_zero_2))
       do ii = 1, num_sup_zeros
          if (ii .NE. (ss-2)) then
             res = res/(sup_zeros_list(ss-2) - sup_zeros_list(ii))
          end if
       end do
       do jj = 1, num_sup_poles
          res = res*(sup_zeros_list(ss-2) - sup_poles_list(jj))
       end do
          
       if (r <= h) then
          fn = (1._dpk - sup_zeros_list(ss-2)*M2)**2/(1._dpk - sup_zeros_list(ss-2)*M1)*compute_Trs(r,sup_zeros_list(ss-2),1)
       else if (r <= 1. .AND. r > h) then
          fn = (1._dpk - sup_zeros_list(ss-2)*M2)*compute_Trs(r,sup_zeros_list(ss-2),2)
       else
          fn = (1._dpk - sup_zeros_list(ss-2)*M3)*compute_Trs(r,sup_zeros_list(ss-2),3)
       end if

       residuepot = omega_r*EXP(CMPLX(0.,1._dpk,kind=dpk)*omega_r*sup_zeros_list(ss-2)*z)*res*fn

!!$       else
!!$          residuepot = 0.
!!$
!!$       end if

    end if

    
  END FUNCTION residuepot

  
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
