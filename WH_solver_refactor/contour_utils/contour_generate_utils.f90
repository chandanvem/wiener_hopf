Module contour_generate_utils
  
  USE input_params 
  USE contour_init_utils
  USE contour_GJ_utils
  
  IMPLICIT NONE
  
  PRIVATE :: compute_kernel_contour, compute_IFT_contour,&
             initialize_contour,initsubdiv,combine,create_stretched_grid
  PUBLIC  :: compute_contours, meshgrid 


  CONTAINS 


  SUBROUTINE compute_contours(input_data,contour_data)

    type(input_params_t)          :: input_data
    type(contour_params_t)        :: contour_data
    integer                       :: i
  
    call compute_kernel_contour(input_data,contour_data)
    call compute_IFT_contour(input_data,contour_data)
   
!    if (input_data%solution_mode == 'guided_jet') then
!       print*,'definecontours: Updating contours'
!       call update_GJ_IFT_cont(input_data,contour_data)
!       call update_GJ_kernel_cont(input_data,contour_data)
!    end if

    print*,''
    print*,'definecontours: Writing kernel integration points to file:'
    open(10,file='initialpoints.out',form='FORMATTED')
   
    do i = 1,contour_data%total_ker_points
       write(10,'(I10,2F30.20)') i,contour_data%ker_int_points(i)
    end do
 
    print*,'definecontours: Writing IFT integration points to file:'
    open(10,file='iftpoints.out',form='FORMATTED')

    do i = 1,contour_data%tot_IFT_pts
       write(10,'(I10,2F30.20)') i, contour_data%iftpoints(i)
    end do

    close(10)

  
  END SUBROUTINE compute_contours
 
  SUBROUTINE compute_kernel_contour(input_data,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Computes the panel lengths for the integration contours
!! 2. Calls subrtn "initialize_contour" which does the the actual contour defining
!! 3. Writes the contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                     :: panel_len_left, panel_len_right  !! length of each panel
    integer                       :: i
    type(input_params_t)          :: input_data
    type(contour_params_t)        :: contour_data

    allocate(contour_data%ker_int_points(contour_data%total_ker_points))

!! length of each panel (kernel contour):

    panel_len_left =  (REAL(contour_data%def_pts_ker_cntr(3))-REAL(contour_data%def_pts_ker_cntr(1))) &
                                     /(input_data%num_ker_pts_loop + 1)  !! panel length in the left section
    panel_len_right = (REAL(contour_data%def_pts_ker_cntr(5))-REAL(contour_data%def_pts_ker_cntr(3))) &
                                     /(input_data%num_ker_pts_loop + 1)  !! panel length in the right section

    call initialize_contour(contour_data%def_pts_ker_cntr,panel_len_left, &
                               panel_len_right,input_data%num_ker_pts_loop,input_data%theta,&
                               1,contour_data%ker_int_points)


 END SUBROUTINE compute_kernel_contour

 SUBROUTINE compute_IFT_contour(input_data,contour_data)


    real(dpk)                     :: panel_len_left, panel_len_right  !! length of each panel
    integer                       :: i
    type(input_params_t)          :: input_data
    type(contour_params_t)        :: contour_data
 

    allocate(contour_data%iftpoints(contour_data%tot_IFT_pts))

!! length of each panel (IFT contour):

    panel_len_left =  (REAL(contour_data%def_pts_IFT_cntr(3))-REAL(contour_data%def_pts_IFT_cntr(1))) &
                                                                     /(input_data%num_IFT_pts_loop+1)
    panel_len_right = (REAL(contour_data%def_pts_IFT_cntr(5))-REAL(contour_data%def_pts_IFT_cntr(3))) &
                                                                     /(input_data%num_IFT_pts_loop+1)

    call initialize_contour(contour_data%def_pts_IFT_cntr,panel_len_left,&
                             panel_len_right,input_data%num_IFT_pts_loop,input_data%theta, & 
                             2,contour_data%iftpoints)

  
 END SUBROUTINE compute_IFT_contour

 SUBROUTINE initialize_contour(def_pts,left_panel_length,right_panel_length,num_pts_per_loop,theta,sw,init_pts_combined)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Initialize the integration contours
!! 2. First compute the individual segments
!! 3. Combine the individual segments to get the full contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk), allocatable, dimension(:)       :: init_pts_left, init_pts_right
    complex(dpk), dimension(2*num_pts_per_loop+3) :: init_pts_combined
    complex(dpk), dimension(5)                    :: def_pts
    real(dpk)                                     :: left_panel_length, right_panel_length, theta
    integer                                       :: num_pts_per_loop, sw


    allocate(init_pts_left(num_pts_per_loop+2))
    allocate(init_pts_right(num_pts_per_loop+2))

    call initsubdiv(REAL(def_pts(1)),def_pts(2)-REAL(def_pts(3)),REAL(def_pts(3)), &
                         num_pts_per_loop,theta,left_panel_length,sw,1,init_pts_left)

    call initsubdiv(REAL(def_pts(3)),def_pts(4)-REAL(def_pts(3)),REAL(def_pts(5)), &
                         num_pts_per_loop,theta,right_panel_length,sw,2,init_pts_right)

    call combine(num_pts_per_loop,init_pts_left,init_pts_right,init_pts_combined)

  END SUBROUTINE initialize_contour


  SUBROUTINE initsubdiv(p,q,r,N,theta,len,sw,ss,zi)
    
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Define the different segments of the integration contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
    
    complex(dpk), dimension(N+2)    :: zi  !! each segment is N+2 long
    real(dpk), dimension(N+2)       :: xi
    real(dpk), dimension(N+2)       :: yi
    complex(dpk)                    :: q
    real(dpk)                       :: p, r, len, theta
    integer                         :: N, i, ss, sw


    if (sw == 1) then  !! indicates kernel contour; equispaced points are a waste
                       !! use some stretching of the grids
       select case (ss)

       case(1)
          call create_stretched_grid(r,p,N+2,theta,1,xi)  

       case(2)
          call create_stretched_grid(p,r,N+2,theta,2,xi)  

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


 SUBROUTINE create_stretched_grid(a,b,N,theta,ss,xst)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Define the stretching function needed for the kernel contour
!! 2. Points are to be clustered near the crossover location
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                      :: a, b
    real(dpk)                      :: xi, xf, dx
    real(dpk)                      :: c, theta
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

  SUBROUTINE meshgrid(input_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Define the physical mesh
!! 2. Write the mesh in PLOT3D format
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)     :: dsz, dsr
    integer       :: i, j
    
    type(input_params_t) :: input_data

    dsr = (input_data%Rmax - input_data%Rmin)/(input_data%Nmeshr - 1)
    dsz = (input_data%Zmax - input_data%Zmin)/(input_data%Nmeshz - 1)

    allocate(input_data%R(input_data%Nmeshr))
    allocate(input_data%Z(input_data%Nmeshz))

    input_data%R(1) = input_data%Rmin
    input_data%R(input_data%Nmeshr) = input_data%Rmax
    input_data%Z(1) = input_data%Zmin
    input_data%Z(input_data%Nmeshz) = input_data%Zmax

    do i = 1, input_data%Nmeshr-2
       input_data%R(i+1) = input_data%R(i) + dsr
    end do
    do i = 1, input_data%Nmeshz-2
       input_data%Z(i+1) = input_data%Z(i) + dsz
    end do

    print*,'meshgrid: Writing mesh grids to file mesh.out'
    open(1,file='mesh.out',form='UNFORMATTED')
    write(1) input_data%Nmeshz,input_data%Nmeshr
    write(1) ((input_data%Z(i),i=1,input_data%Nmeshz),j=1,input_data%Nmeshr),((input_data%R(j),&
         i=1,input_data%Nmeshz),j=1,input_data%Nmeshr)
    close(1)

   print*,'meshgrid: Writing mesh grids to file mesh.r'
    open(20,file='mesh.r')
    do i = 1, input_data%Nmeshr
       write(20,'(I10,2X,F20.10)')  i, input_data%R(i)
    end do
    close(20)


    print*,'meshgrid: Writing mesh grids to file mesh.z'
    open(20,file='mesh.z')
    do i = 1, input_data%Nmeshz
       write(20,'(I10,2X,F20.10)') i, input_data%Z(i)
    end do
    close(20)

    if (input_data%farswitch == 2) then
       open(1,file='mesh_polar.out',form='UNFORMATTED')
       write(1) input_data%Nphi,input_data%Nmeshr
       write(1) ((input_data%phi(i),i=1,input_data%Nphi),j=1,input_data%Nmeshr),((input_data%R(j),&
            i=1,input_data%Nphi),j=1,input_data%Nmeshr)
       close(1)
    end if


  END SUBROUTINE meshgrid

end module contour_generate_utils

