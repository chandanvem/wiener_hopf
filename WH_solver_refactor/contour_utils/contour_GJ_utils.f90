Module contour_GJ_utils
  
  USE input_params 
  
  IMPLICIT NONE
  
  PRIVATE :: solve_gaussian, compute_FD_sec_order, compute_A_B_matrices,&
             get_pos_idx, get_y_nth_degree_poly

  PUBLIC  :: update_GJ_IFT_cont, get_dy_dx_alg_int_contour


  CONTAINS 

  SUBROUTINE update_GJ_IFT_cont(input_data,contour_data)

      type(input_params_t)   :: input_data
      type(contour_params_t) :: contour_data
  
      real(dpk) :: temp1(1), temp3(1)
      integer   :: max_IFT_idx, temp2(1), temp4(1), cross_over_pt_idx, i
      complex(dpk)  :: max_IFT_cont, cross_over_pt

      real(dpk) :: dx, slope_at_cross_over
      real(dpk) :: A(4,4), B(4)
      real(dpk) :: poly_coeffs_1(4), poly_coeffs_2(4), poly_coeffs_3(4)
      real(dpk) :: x_in,imag_part
 
      integer   :: pos_idx
 
      temp2 =  maxloc(AIMAG(contour_data%iftpoints))
      max_IFT_idx = temp2(1)
      
      max_IFT_cont = contour_data%iftpoints(max_IFT_idx)

      temp3 = minloc(ABS(AIMAG(contour_data%iftpoints)))
      cross_over_pt_idx = temp3(1)

      cross_over_pt = contour_data%iftpoints(cross_over_pt_idx)
      dx =  REAL(contour_data%iftpoints(cross_over_pt_idx + 1) - contour_data%iftpoints(cross_over_pt_idx))

      slope_at_cross_over = compute_FD_sec_order(AIMAG(contour_data%iftpoints(cross_over_pt_idx)),&
                                                AIMAG(contour_data%iftpoints(cross_over_pt_idx + 1)),&
                                                AIMAG(contour_data%iftpoints(cross_over_pt_idx + 2)),&
                                                 dx, dx)

      call compute_A_B_matrices(max_IFT_cont,contour_data%GJ_ref_pt,0._dpk,0._dpk,A,B)  
      call solve_gaussian(A, B, poly_coeffs_1, 4)
    
      call compute_A_B_matrices(contour_data%GJ_ref_pt,contour_data%GJ_cntr_maxima,0._dpk,0._dpk,A,B)  
      call solve_gaussian(A, B, poly_coeffs_2, 4)
   
      call compute_A_B_matrices(contour_data%GJ_cntr_maxima,cross_over_pt,0._dpk,slope_at_cross_over,A,B)  
      call solve_gaussian(A, B, poly_coeffs_3, 4)


      do i= max_IFT_idx,cross_over_pt_idx

         x_in = REAL(contour_data%iftpoints(i))
  
         call get_pos_idx(x_in,REAL(max_IFT_cont),REAL(contour_data%GJ_ref_pt),REAL(contour_data%GJ_cntr_maxima), &
                                                                               REAL(cross_over_pt),pos_idx)

         select case (pos_idx)
            case (1)
                imag_part = get_y_nth_degree_poly(x_in, poly_coeffs_1, 4)
            case (2)
                imag_part = get_y_nth_degree_poly(x_in, poly_coeffs_2, 4)
            case (3)
                imag_part = get_y_nth_degree_poly(x_in, poly_coeffs_3, 4)
            case default
                print*,'update_GJ_IFT_cont: STOPPING. out of bounds for three cubic contours'
                STOP
         end select
         
         contour_data%iftpoints(i) = CMPLX(x_in,imag_part,kind=dpk)

     end do
 
  END SUBROUTINE update_GJ_IFT_cont
  
FUNCTION compute_FD_sec_order(f0, f1, f2, dx1, dx2) result(df)

  real(dpk), intent(in)  :: f0, f1, f2   
  real(dpk), intent(in)  :: dx1, dx2           
  real(dpk)              :: df         

  if (dx1 == dx2) then
     df = (-3.0*f0 + 4.0*f1 - f2) / (2.0*dx1)
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

FUNCTION get_y_nth_degree_poly(x_in,poly_coeffs,n) result(y)

  integer                :: n, i
  real(dpk), intent(in)  :: x_in 
  real(dpk), intent(in)  :: poly_coeffs(n)          
  real(dpk)              :: y
  
  y = poly_coeffs(1)
  do i = 2, n
     y = y * x_in + poly_coeffs(i)
  end do

END FUNCTION get_y_nth_degree_poly  

SUBROUTINE get_pos_idx(x_in, a, b, c, d, pos_idx)
  
  real(8), intent(in) :: x_in, a, b, c, d
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
  
SUBROUTINE get_dy_dx_alg_int_contour(x, dy_dx, p, a, b)
  !! Computes exact derivative dy/dx of the algebraic contour

    real(dpk), intent(in)  :: x, a, b, p
    real(dpk), intent(out) :: dy_dx

    real(dpk) :: u

    u  = (x - p) / a
    dy_dx = b * (12.0_dpk * (1.0_dpk - u**4)) / (a * (3.0_dpk + u**4)**2)

END SUBROUTINE get_dy_dx_alg_int_contour



END MODULE contour_GJ_utils

!program solve_axb
!  implicit none
!  integer, parameter :: dpk = kind(1.0d0)
!  real(dpk) :: A(8,8), B(8), X(8)
!  integer :: i, j
!
!  ! Example system: A * X = B
!  A = reshape([ &
!      4.0_dpk, -1.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk, &
!     -1.0_dpk,  4.0_dpk, -1.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk, &
!      0.0_dpk, -1.0_dpk,  4.0_dpk, -1.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk, &
!      0.0_dpk,  0.0_dpk, -1.0_dpk,  4.0_dpk, -1.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk, &
!      0.0_dpk,  0.0_dpk,  0.0_dpk, -1.0_dpk,  4.0_dpk, -1.0_dpk,  0.0_dpk,  0.0_dpk, &
!      0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk, -1.0_dpk,  4.0_dpk, -1.0_dpk,  0.0_dpk, &
!      0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk, -1.0_dpk,  4.0_dpk, -1.0_dpk, &
!      0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk,  0.0_dpk, -1.0_dpk,  4.0_dpk  ],
!&
!      shape(A))
!
!  B = [1.0_dpk, 2.0_dpk, 3.0_dpk, 4.0_dpk, 5.0_dpk, 6.0_dpk, 7.0_dpk, 8.0_dpk]
!
!  call gauss_solve(A, B, 8, X)
!
!  print *, "Solution X:"
!  do i = 1, 8
!     print '(F12.6)', X(i)
!  end do
!
!

!     print*,'max_IFT_idx = ', max_IFT_idx
!      print*,'cross_over_pt_idx ', cross_over_pt_idx 

!      print*, 'max IFT contour =', max_IFT_cont
!      print*, 'GJ_ref_pt = ', contour_data%GJ_ref_pt
!      print*, 'GJ_cntr_maxima = ', contour_data%GJ_cntr_maxima
!      print*, 'cross_over_pt = ', cross_over_pt
!      print*, 'slope_at_cross_over =', slope_at_cross_over

  
!       print*,'------------------------------------'
!      do i = 1, 4
!        print* , poly_coeffs_3(i)   ! 8 columns, width=8, 2 decimals
!      end do
!      print*,'-------------------------------------'
  
