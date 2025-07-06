Module contour_utils
  
  USE input_params 
  
  IMPLICIT NONE
  
  PRIVATE
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


  SUBROUTINE get_y_alg_int_contour(x,y,p,a,b)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Algebraic "one piece" contour:  see Rienstra, JEM, 2007
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)    :: x, y
    real(dpk)    :: a, b, p

    y = b*(4._dpk*(x-p)/a)/(3._dpk + ((x-p)/a)**4) ! algebraic


  END SUBROUTINE get_y_alg_int_contour



end module contour_utils

