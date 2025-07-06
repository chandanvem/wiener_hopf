Module contour_utils
  
  USE input_params 
  
  IMPLICIT NONE
  
  PRIVATE :: check_location_wrt_contour, check_location_of_zeros_poles
  PUBLIC  :: compute_contour_params, get_y_alg_int_contour 


  CONTAINS 
 
     SUBROUTINE compute_contour_params(input_data,contour_data)

        !! read the contour points, only IFT contour is read in full;
        !! the end points of the kernel integration contour read while
        !! the rest are obtained via the "offset" parameter
        !! def_pts_IFT_cntr(1), def_pts_IFT_cntr(5) are end points; def_pts_IFT_cntr(2), def_pts_IFT_cntr(4) are where the
        !! algebraic contour starts & ends; def_pts_IFT_cntr(3) is where it crosses
        !! real axis

            type(input_params_t)   :: input_data
            type(contour_params_t) :: contour_data
            
            integer             :: i, i1, i2, i3
            real(dpk)           :: t, ai, ai_im, af, af_im
            real(dpk)           :: aki, aki_im, akf, akf_im
            character(60)       :: pos
 
            allocate(contour_data%def_pts_IFT_cntr(5))
            allocate(contour_data%def_pts_ker_cntr(5))
             
            open(10,status='old',file='input.iftpts')

            read(10,*) t
            ai = t
            read(10,*) contour_data%def_pts_IFT_cntr(2)
            read(10,*) t
            contour_data%def_pts_IFT_cntr(3) = CMPLX(t,0._dpk,kind=dpk)
            read(10,*) contour_data%def_pts_IFT_cntr(4)
            read(10,*) t
            af = t

            call get_y_alg_int_contour(ai,ai_im,REAL(contour_data%def_pts_IFT_cntr(3)), &
                                     REAL(contour_data%def_pts_IFT_cntr(2)-contour_data%def_pts_IFT_cntr(3)), &
                                        AIMAG(contour_data%def_pts_IFT_cntr(2)-contour_data%def_pts_IFT_cntr(3)))

            call get_y_alg_int_contour(af,af_im,REAL(contour_data%def_pts_IFT_cntr(3)), &
                                        REAL(contour_data%def_pts_IFT_cntr(4)-contour_data%def_pts_IFT_cntr(3)), &
                                             AIMAG(contour_data%def_pts_IFT_cntr(4)-contour_data%def_pts_IFT_cntr(3)))

            contour_data%def_pts_IFT_cntr(1) = CMPLX(ai,ai_im,kind=dpk)  
            contour_data%def_pts_IFT_cntr(5) = CMPLX(af,af_im,kind=dpk)

            read(10,*) t
            aki = t
            read(10,*) t
            akf = t

            close(1)

        !! now determine the kernel integration contour points:
        !! NOTE: the kernel contour lies above the IFT one and thus it's further away from the real
        !! axis when above, but closer to it when below

            contour_data%def_pts_ker_cntr(2) =  contour_data%def_pts_IFT_cntr(2) + CMPLX(0._dpk,input_data%offset,kind=dpk)  !! Kernel points
            contour_data%def_pts_ker_cntr(3) =  contour_data%def_pts_IFT_cntr(3) + CMPLX(5._dpk*(input_data%offset),0._dpk,kind=dpk)
            contour_data%def_pts_ker_cntr(4) =  contour_data%def_pts_IFT_cntr(4) + CMPLX(0._dpk,(input_data%offset),kind=dpk)

            call get_y_alg_int_contour(aki,aki_im,REAL(contour_data%def_pts_ker_cntr(3)), &
                                   REAL(contour_data%def_pts_ker_cntr(2)-contour_data%def_pts_ker_cntr(3)), &
                                   AIMAG(contour_data%def_pts_ker_cntr(2)-contour_data%def_pts_ker_cntr(3)))
            call get_y_alg_int_contour(akf,akf_im,REAL(contour_data%def_pts_ker_cntr(3)), &
                                   REAL(contour_data%def_pts_ker_cntr(4)-contour_data%def_pts_ker_cntr(3)), &
                                   AIMAG(contour_data%def_pts_ker_cntr(4)-contour_data%def_pts_ker_cntr(3)))
            
            contour_data%def_pts_ker_cntr(1) = CMPLX(aki,aki_im,kind=dpk)
            contour_data%def_pts_ker_cntr(5) = CMPLX(akf,akf_im,kind=dpk)

        !! find where the contour crosses the real axis (at "cont_cross_over_pt"):

            contour_data%cont_cross_over_pt = REAL(contour_data%def_pts_ker_cntr(3))

            contour_data%total_ker_points = 2*(input_data%num_ker_pts_loop) + 3  !! total number of kernel contor pts

            contour_data%tot_IFT_pts = 2*(input_data%num_IFT_pts_loop) + 3  !! total number of IFT contor pts


            call check_location_of_zeros_poles(input_data,contour_data)

        !! print the contour points:

            write(*,'(/A26)') 'initialize: The Integration Contours:'
            write(*,'(/A27)') '1. The kernel integration contour:'
            write(*,'(/A14)') 'key points:->'
            do i1 = 1, 5
               write(pos,"(I2,' =')") i1
               write(*,'(/1X,A5,2F15.6)')  'i'//pos, contour_data%def_pts_ker_cntr(i1)
            end do

            write(*,'(/A30)') '2. The inv Fourier transform contour:'
            write(*,'(/A14)') 'key points:->' 
            do i1 = 1, 5
               write(pos,"(I2,' =')") i1
               write(*,"(/1X,A5,2F15.6)")  'i'//pos, contour_data%def_pts_IFT_cntr(i1) 
            end do

           print '(A, I0)', 'initialize: total number of IFT points = ', contour_data%tot_IFT_pts 
           


     END SUBROUTINE compute_contour_params


     SUBROUTINE check_location_of_zeros_poles(input_data,contour_data)

             type(input_params_t)      :: input_data
             type(contour_params_t)    :: contour_data 
             integer                   :: i1, i2, i3

!! here both the instability zeros and the pole need to lie outside (under) the contour,
!! since we use the residue theorem to compute them anyway:

             call check_location_wrt_contour(input_data%KH_zero_1,i1,contour_data) 
             call check_location_wrt_contour(input_data%KH_zero_2,i2,contour_data)
             call check_location_wrt_contour(input_data%KH_pole_1,i3,contour_data)

!! if any of them lie inside ask to redefine the contour:

    if (input_data%vortswitch .EQ. 0) then

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

    do i1 = 1, input_data%num_sup_zeros
       call check_location_wrt_contour(input_data%sup_zeros_list(i1),i2,contour_data)
       if (input_data%vortswitch .EQ. 0) then
          if(i2 == 0) then
             print*,''
             print*,'initialize: Redefine the contour: The following supersonic zero is inside:'
             print*,input_data%sup_zeros_list(i1)
          end if
       else
          if(i2 == 1) then
             print*,''
             print*,'initialize: Redefine the contour: The following supersonic zero is outside:'
             print*,input_data%sup_zeros_list(i1)
          end if
       end if
    end do

    do i1 = 1, input_data%num_sup_poles
       call check_location_wrt_contour(input_data%sup_poles_list(i1),i2,contour_data)
       if (input_data%vortswitch .EQ. 0) then
          if(i2 == 0) then
             print*,''
             print*,'initialize: Redefine the contour: The following supersonic pole is inside:'
             print*,input_data%sup_poles_list(i1)
          end if
       else
          if(i2 == 1) then
            print*,''
            print*,'initialize: Redefine the contour: The following supersonic pole is outside:'
            print*,input_data%sup_poles_list(i1)
          end if
       end if
    end do


     END SUBROUTINE check_location_of_zeros_poles 

     SUBROUTINE check_location_wrt_contour(z,switch,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Subroutine to check the location of the instability zeros and pole wrt to
!the contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

      complex(dpk)    :: z
      real(dpk)       :: zx, zy, zfy
      integer         :: switch
      type(contour_params_t) :: contour_data

      zx = REAL(z)
      zy = AIMAG(z)

      call get_y_alg_int_contour(zx,zfy,REAL(contour_data%def_pts_IFT_cntr(3)), &
           REAL(contour_data%def_pts_IFT_cntr(4)-contour_data%def_pts_IFT_cntr(3)), &
           AIMAG(contour_data%def_pts_IFT_cntr(4)-contour_data%def_pts_IFT_cntr(3)))
!!$    zfy = AIMAG(c1)

    if(zy > zfy) then
       switch = 0  !! inside
    else
       switch = 1  !! outside
    end if


  END SUBROUTINE check_location_wrt_contour


  SUBROUTINE get_y_alg_int_contour(x,y,p,a,b)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Algebraic "one piece" contour:  see Rienstra, JEM, 2007
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)    :: x, y
    real(dpk)    :: a, b, p

    y = b*(4._dpk*(x-p)/a)/(3._dpk + ((x-p)/a)**4) ! algebraic


  END SUBROUTINE get_y_alg_int_contour



end module contour_utils

