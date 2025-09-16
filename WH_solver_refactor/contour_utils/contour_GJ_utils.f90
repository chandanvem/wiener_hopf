Module contour_GJ_utils
  
  USE input_params 
  
  IMPLICIT NONE
  
  PRIVATE :: gauss_solve 
  PUBLIC  :: update_GJ_IFT_cont, get_dy_dx_alg_int_contour


  CONTAINS 


  SUBROUTINE update_GJ_IFT_cont(input_data,contour_data)

      type(input_params_t)   :: input_data
      type(contour_params_t) :: contour_data
  
      real(dpk) :: max_IFT_imag, temp1(1)
      integer ::   max_IFT_idx, temp2(1)

      temp1 = maxval(AIMAG(contour_data%iftpoints))
      max_IFT_imag = temp1(1)

      temp2 =  maxloc(AIMAG(contour_data%iftpoints))
      max_IFT_idx = temp2(1)

  END SUBROUTINE update_GJ_IFT_cont

 
  SUBROUTINE get_dy_dx_alg_int_contour(x, dy_dx, p, a, b)
  !! Computes exact derivative dy/dx of the algebraic contour

    real(dpk), intent(in)  :: x, a, b, p
    real(dpk), intent(out) :: dy_dx

    real(dpk) :: u

    u  = (x - p) / a
    dy_dx = b * (12.0_dpk * (1.0_dpk - u**4)) / (a * (3.0_dpk + u**4)**2)

  END SUBROUTINE get_dy_dx_alg_int_contour


  SUBROUTINE gauss_solve(A, b, n, x)

    real(dpk), intent(inout) :: A(n,n)
    real(dpk), intent(inout) :: b(n)
    integer, intent(in)     :: n
    real(dpk), intent(out)   :: x(n)
    integer :: i, j, k, piv
    real(dpk) :: factor, tmp, maxval

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
          A([k,piv],:) = A([piv,k],:)  ! swap rows
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
  
END SUBROUTINE gauss_solve


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
!
end module contour_GJ_utils


