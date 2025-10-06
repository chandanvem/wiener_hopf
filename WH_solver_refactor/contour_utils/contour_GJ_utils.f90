Module contour_GJ_utils
  
  USE input_params 
  
  IMPLICIT NONE
  
  PRIVATE :: solve_gaussian, compute_FD_sec_order, compute_A_B_matrices,&
             get_pos_idx, get_y_nth_degree_poly, get_dy_nth_degree_poly,&
             compute_offset_point, swap_arrays

  PUBLIC  :: update_GJ_IFT_cont, update_GJ_kernel_cont


  CONTAINS 

  SUBROUTINE update_GJ_IFT_cont(input_data,contour_data)

      type(input_params_t)   :: input_data
      type(contour_params_t) :: contour_data
  
      real(dpk) :: temp1(1), temp3(1)
      integer   :: max_IFT_idx, temp2(1), temp4(1), cross_over_pt_idx, i
      complex(dpk)  :: max_IFT_cont, cross_over_pt

      real(dpk) :: dx, slope_at_cross_over
      real(dpk) :: A(4,4), B(4)
      real(dpk) :: x_in,imag_part
 
      integer   :: pos_idx
 
      temp2 =  maxloc(AIMAG(contour_data%iftpoints))
      max_IFT_idx = temp2(1)
      
      max_IFT_cont = contour_data%iftpoints(max_IFT_idx)

      temp3 = minloc(ABS(AIMAG(contour_data%iftpoints)))
      cross_over_pt_idx = temp3(1)

      cross_over_pt = contour_data%iftpoints(cross_over_pt_idx)

      contour_data%IFT_cross_over_pt = cross_over_pt

      dx =  REAL(contour_data%iftpoints(cross_over_pt_idx + 1) - contour_data%iftpoints(cross_over_pt_idx))

      slope_at_cross_over = compute_FD_sec_order(AIMAG(contour_data%iftpoints(cross_over_pt_idx)),&
                                                AIMAG(contour_data%iftpoints(cross_over_pt_idx + 1)),&
                                                AIMAG(contour_data%iftpoints(cross_over_pt_idx + 2)),&
                                                 dx, dx)

      call compute_A_B_matrices(max_IFT_cont,contour_data%GJ_ref_pt,0._dpk,0._dpk,A,B)  
      call solve_gaussian(A, B, contour_data%poly_coeffs_cubic_1, 4)
    
      call compute_A_B_matrices(contour_data%GJ_ref_pt,contour_data%GJ_cntr_maxima,0._dpk,0._dpk,A,B)  
      call solve_gaussian(A, B, contour_data%poly_coeffs_cubic_2, 4)
   
      call compute_A_B_matrices(contour_data%GJ_cntr_maxima,cross_over_pt,0._dpk,slope_at_cross_over,A,B)  
      call solve_gaussian(A, B, contour_data%poly_coeffs_cubic_3, 4)

      print*,'update_GJ_IFT_cont: num of points = ', cross_over_pt_idx - max_IFT_idx + 1

      do i= max_IFT_idx,cross_over_pt_idx

         x_in = REAL(contour_data%iftpoints(i))
  
         call get_pos_idx(x_in,REAL(max_IFT_cont),REAL(contour_data%GJ_ref_pt),REAL(contour_data%GJ_cntr_maxima), &
                                                                               REAL(cross_over_pt),pos_idx)

         select case (pos_idx)
            case (1)
                imag_part = get_y_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_1, 4)
            case (2)
                imag_part = get_y_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_2, 4)
            case (3)
                imag_part = get_y_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_3, 4)
            case default
                print*,'update_GJ_IFT_cont: STOPPING. out of bounds for three cubic contours'
                STOP
         end select
         
         contour_data%iftpoints(i) = CMPLX(x_in,imag_part,kind=dpk)

     end do
 
  END SUBROUTINE update_GJ_IFT_cont

  SUBROUTINE update_GJ_kernel_cont(input_data,contour_data)

      type(input_params_t)   :: input_data
      type(contour_params_t) :: contour_data
  
      real(dpk) :: temp1(1), temp3(1)
      integer   :: max_IFT_idx, temp2(1), temp4(1), cross_over_pt_idx, i, max_kernel_idx
      integer   :: cross_over_pt_idx_kernel

      complex(dpk)  :: max_IFT_cont, max_kernel_cont
      complex(dpk)  :: cross_over_pt_kernel, cross_over_pt_IFT, start_new_kernel_contour
      real(dpk)     :: dx, imag_part, offset,  x_in, x_out, y_out, slope

      integer   :: pos_idx, updated_total_num_kernel_pts
      complex(dpk), allocatable, dimension(:)    :: kernel_pts_cubics_list
      complex(dpk), allocatable, dimension(:)    :: updated_ker_int_points

      integer   :: j, num_temp_points, count_temp_point, updated_cross_over_pt_index
      complex(dpk)  :: start_point, end_point


      allocate(kernel_pts_cubics_list(contour_data%num_ker_pts_in_cubics(1) + &
                                      contour_data%num_ker_pts_in_cubics(2) + &
                                      contour_data%num_ker_pts_in_cubics(3) - 2 ))
 
      temp2 =  maxloc(AIMAG(contour_data%ker_int_points))
      max_kernel_idx = temp2(1)
      max_kernel_cont = contour_data%ker_int_points(max_kernel_idx)

      temp2 =  maxloc(AIMAG(contour_data%iftpoints))
      max_IFT_idx = temp2(1)
      max_IFT_cont = contour_data%iftpoints(max_IFT_idx)

      temp3 = minloc(ABS(AIMAG(contour_data%ker_int_points)))
      cross_over_pt_idx_kernel = temp3(1)

      start_new_kernel_contour =  contour_data%ker_int_points(max_kernel_idx + 1)

      cross_over_pt_kernel = contour_data%ker_int_points(cross_over_pt_idx_kernel)
      cross_over_pt_IFT    = contour_data%IFT_cross_over_pt

  !    open(10,file='kernelpoints_temp.out',form='FORMATTED')

      count_temp_point = 1
      do j = 1,3
        select case (j)
               case (1)
                  start_point =  max_IFT_cont
                  end_point   =  contour_data%GJ_ref_pt 
                  num_temp_points = contour_data%num_ker_pts_in_cubics(1)-1
               case (2)
                  start_point =  contour_data%GJ_ref_pt 
                  end_point   =  contour_data%GJ_cntr_maxima
                  num_temp_points = contour_data%num_ker_pts_in_cubics(2)-1
               case (3)
                  start_point =  contour_data%GJ_cntr_maxima
                  end_point   =  contour_data%IFT_cross_over_pt
                  num_temp_points = contour_data%num_ker_pts_in_cubics(3)
        end select
    
        dx = REAL(end_point- start_point)/(contour_data%num_ker_pts_in_cubics(j) - 1)

        do i = 1,num_temp_points
     
           x_in  = REAL(start_point) + (i-1)*dx

           if (i == 1) then
               x_in = REAL(start_point)
           else if ( i == num_temp_points) then
               x_in = REAL(end_point)           
           end if

           call get_pos_idx(x_in,REAL(max_IFT_cont),REAL(contour_data%GJ_ref_pt),REAL(contour_data%GJ_cntr_maxima), &
                                                                                  REAL(cross_over_pt_IFT),pos_idx)

           select case (pos_idx)
              case (1)
                  imag_part =  get_y_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_1, 4)
                  slope     =  get_dy_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_1, 4)
              case (2)
                  imag_part = get_y_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_2, 4)
                  slope     = get_dy_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_2, 4)
              case (3)
                  imag_part = get_y_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_3, 4)
                  slope     = get_dy_nth_degree_poly(x_in, contour_data%poly_coeffs_cubic_3, 4)
              case default
                  print*,'update_GJ_IFT_cont: STOPPING. out of bounds for three cubic contours for cubic',j,&
                                                                                          'i = ',i,'x_in = ', x_in 
                  STOP
           end select
        
           call compute_offset_point(x_in, imag_part, slope, input_data%offset, x_out, y_out)
          
           !print*,'x_in = ',x_in, 'y_in = ',imag_part, 'slope = ', slope, 'x_out = ', x_out, 'y_out = ', y_out

           kernel_pts_cubics_list(count_temp_point) = CMPLX(x_out,y_out,kind=dpk)

           if(j == 3 .AND. i == num_temp_points) then
               kernel_pts_cubics_list(count_temp_point)  = cross_over_pt_kernel
           end if

     !      write(10,'(I10,2F30.20)') count_temp_point,kernel_pts_cubics_list(count_temp_point)
           count_temp_point = count_temp_point + 1
         end do     
       end do
    ! close(10)

     
    updated_cross_over_pt_index = max_kernel_idx + count_temp_point - 1
    updated_total_num_kernel_pts = updated_cross_over_pt_index + &
                                   contour_data%total_ker_points - cross_over_pt_idx_kernel


    allocate(updated_ker_int_points(updated_total_num_kernel_pts))
    
    open(10,file='kernelpoints_temp.out',form='FORMATTED')
  
    do i = 1, updated_total_num_kernel_pts
        if (i <= max_kernel_idx) then
            updated_ker_int_points(i) = contour_data%ker_int_points(i)
        else if (i > max_kernel_idx .AND. i <= updated_cross_over_pt_index ) then
           updated_ker_int_points(i) = kernel_pts_cubics_list(i- max_kernel_idx)
        else
           updated_ker_int_points(i) = contour_data%ker_int_points(i + cross_over_pt_idx_kernel &
                                                                 - updated_cross_over_pt_index )
        end if
       write(10,'(I10,2F30.20)') i,updated_ker_int_points(i)
     print*,'i =  ',i,'x = ', REAL(updated_ker_int_points(i)), 'y_out = ',AIMAG(updated_ker_int_points(i))

     end do
    close(10)
      
    call swap_arrays(contour_data%ker_int_points,updated_ker_int_points) 
!
    contour_data%total_ker_points = updated_total_num_kernel_pts

 END SUBROUTINE update_GJ_kernel_cont
 
   
 FUNCTION compute_FD_sec_order(f0, f1, f2, dx1, dx2) result(df)

    real(dpk), intent(in)  :: f0, f1, f2   
    real(dpk), intent(in)  :: dx1, dx2           
    real(dpk)              :: df         

    if (dx1 == dx2) then
       df = (-3._dpk*f0 + 4._dpk*f1 - f2) / (2._dpk*dx1)
    else
       df =  -(2._dpk*dx1 + dx2) / (dx1*(dx1+dx2)) * f0  &
             + ((dx1+dx2) / (dx1*dx2)) * f1              &
             - (dx1 / (dx2*(dx1+dx2))) * f2
    end if
   
  END FUNCTION compute_FD_sec_order

  SUBROUTINE compute_A_B_matrices(z1,z2,slope1,slope2,A,B)

    complex(dpk), intent(in)  :: z1,z2
    real(dpk), intent(in)     :: slope1, slope2
    real(dpk), intent(out)    :: A(4,4), B(4)

    integer :: i,n

    n = 4

     do i = 1,n
        A(1,i) = (REAL(z1))**(n-i)
     end do  
    
     do i = 1,n-1
        if (REAL(z1) == 0._dpk .AND. n-i-1 .EQ. 0) then
           A(2,i) = 1._dpk
        else
           A(2,i) = (n-i)*(REAL(z1))**(n-i-1)
       end if
     end do  
     A(2,n) = 0._dpk
     
     do i = 1,n
        A(3,i) = (REAL(z2))**(n-i)
     end do  
   
     do i = 1,n-1
      if (REAL(z2) == 0._dpk .AND. n-i-1 .EQ. 0) then
           A(4,i) = 1._dpk
        else
           A(4,i) = (n-i)*(REAL(z2))**(n-i-1)
       end if
     end do  
     A(4,n) = 0._dpk

     B = [AIMAG(z1), slope1, AIMAG(z2), slope2]
     

   END SUBROUTINE 

  SUBROUTINE solve_gaussian(A, b, x, n)

     integer, intent(in) :: n
     real(dpk), intent(inout) :: A(n,n)
     real(dpk), intent(inout) :: b(n)
     real(dpk), intent(out)   :: x(n)

     integer :: i, j, k, piv
     real(dpk) :: factor, maxval, tmp
     real(dpk) :: row(n)

     ! Forward elimination with partial pivoting
     do k = 1, n-1
        piv = k
        maxval = abs(A(k,k))
        do i = k+1, n
           if (abs(A(i,k)) > maxval) then
              piv = i
              maxval = abs(A(i,k))
           end if
        end do

        if (piv /= k) then
           row      = A(k,:)   ! swap rows in A
           A(k,:)   = A(piv,:)
           A(piv,:) = row

           tmp = b(k); b(k) = b(piv); b(piv) = tmp
        end if

        do i = k+1, n
           factor = A(i,k) / A(k,k)
           A(i,k) = 0.0_dpk
           A(i,k+1:n) = A(i,k+1:n) - factor * A(k,k+1:n)
           b(i) = b(i) - factor * b(k)
        end do
     end do

     ! Back substitution
     do i = n, 1, -1
        tmp = b(i)
        do j = i+1, n
           tmp = tmp - A(i,j) * x(j)
        end do
        x(i) = tmp / A(i,i)
     end do

  END SUBROUTINE solve_gaussian

  FUNCTION get_y_nth_degree_poly(x_in,poly_coeffs_cubic,n) result(y)

    integer                :: n, i
    real(dpk), intent(in)  :: x_in 
    real(dpk), intent(in)  :: poly_coeffs_cubic(n)          
    real(dpk)              :: y
    
    y = poly_coeffs_cubic(1)
    do i = 2, n
       y = y * x_in + poly_coeffs_cubic(i)
    end do

  END FUNCTION get_y_nth_degree_poly  


  FUNCTION get_dy_nth_degree_poly(x_in,poly_coeffs_cubic,n) result(y)

    integer                :: n, i
    real(dpk), intent(in)  :: x_in 
    real(dpk), intent(in)  :: poly_coeffs_cubic(n)          
    real(dpk)              :: y
    
    y = (n-1)*poly_coeffs_cubic(1)
    do i = 2, n-1
       y = y * x_in + (n - i)* poly_coeffs_cubic(i)
    end do

  END FUNCTION get_dy_nth_degree_poly  


  SUBROUTINE get_pos_idx(x_in, a, b, c, d, pos_idx)
    
    real(dpk), intent(in) :: x_in, a, b, c, d
    integer, intent(out) :: pos_idx

    pos_idx = 0  

    if ( (x_in >= min(a,b)) .AND. (x_in <= max(a,b)) ) then
       pos_idx = 1
    else if ( (x_in >= min(b,c)) .AND. (x_in <= max(b,c)) ) then
       pos_idx = 2
    else if ( (x_in >= min(c,d)) .AND. (x_in <= max(c,d)) ) then
       pos_idx = 3
    end if

  END SUBROUTINE get_pos_idx
    

  SUBROUTINE compute_offset_point(x_in, y_in, slope, r, x_out, y_out)

    real(dpk), intent(in)  :: x_in, y_in, slope, r
    real(dpk), intent(out) :: x_out, y_out
    real(dpk) :: dx, dy, denom

    if (ABS(slope) > 1.0d12) then
       dx = 0._dpk
       dy = r
    else
       denom = sqrt(1._dpk + slope*slope)
       dx = r*slope / denom
       dy = r  / denom
    end if

    if (slope < 0._dpk .OR. slope .EQ. 0) then
       x_out = x_in + ABS(dx)
       y_out = y_in + ABS(dy)
    else 
       x_out = x_in - ABS(dx)
       y_out = y_in + ABS(dy)
    end if

  END SUBROUTINE compute_offset_point

  SUBROUTINE swap_arrays(A, B)

    complex(dpk), allocatable, intent(inout) :: A(:), B(:)
    complex(dpk), allocatable :: tmp(:)

    call move_alloc(A, tmp)
    call move_alloc(B, A)
    call move_alloc(tmp, B)
    
  END SUBROUTINE swap_arrays


END MODULE contour_GJ_utils

 
