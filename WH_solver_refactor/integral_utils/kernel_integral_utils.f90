Module kernel_integral_utils

  USE input_params
  USE contour_generate_utils
  USE contour_init_utils


  IMPLICIT NONE

  PRIVATE :: kernel_trapz_int 
  PUBLIC  :: compute_eqn_A1_integral,sum_panel_contributions_kernel

  SUBROUTINE compute_eqn_A1_integral(s_target,integral_value,kswitch,ch,KU_K_switch, &
                                                              input_data,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Evaluate (A1)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)           :: s_target
    complex(dpk)           :: integral_value ! final value of integration
    integer                :: num_of_quad_points, ch
    integer                :: KU_K_switch  !! 1: K/U; else: K
    integer                :: kswitch  !! 1: compute (4.10); else: compute (3.22)
    complex(dpk), allocatable, dimension(:)    :: intpanel !! integral value at each panel
    integer(dpk), allocatable, dimension(:)    ::  Npanel 

    type(input_params_t)   :: input_data
    type(contour_params_t) :: contour_data

    allocate(intpanel(contour_data%total_ker_points-1))
    allocate(Npanel(contour%total_ker_points-1))

    call adaptive(s_target,kswitch,ch,KU_K_switch,intpanel,Npanel,input_data,contour_data)

    call sum_panel_contributions_kernel(integral_value,num_of_quad_points,s_target, &
                                                                     intpanel,Npanel)

    deallocate(intpanel)
    deallocate(Npanel)


  END SUBROUTINE compute_eqn_A1_integral


  SUBROUTINE adaptive(si,ksw,check,sw,intpanel,Npanel)
    
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The kernel integration is computed using an adaptive routine; which not feasible
!!    for the IFT integral
!! 2. This is integrated with the trapezoidal rule used for quadrature
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)  ::  intpanel(total_ker_points-1)
    integer(dpk)  ::  Npanel(total_ker_points-1)
    complex(dpk), allocatable, dimension(:)    :: zp, T_temp
    complex(dpk), allocatable, dimension(:)    :: knl, zall, zall_temp
    real(dpk), allocatable, dimension(:)       :: xp, yp
    real(dpk)                                  :: ds, len
    complex(dpk)                               :: T, si, scale
    integer                                    :: Np, Nq, Nr, panel_no, Nloopindex, check, sw, ksw
    integer                                    :: i, j, ii


    do j = 1, contour_data%total_ker_points-1  !! the integration points already specified

       len = contour_data%ker_int_points(j+1) - contour_data%ker_int_points(j)  !! the length of a mini panel

!! compute "scale" needed to ensure the tolerance check:

       if (j >= 1 .AND. j < input_data%num_ker_pts_loop+2 ) then
          scale = len/(input_data%def_pts_ker_cntr(3) - input_data%def_pts_ker_cntr(1))

       else
          scale = len/(input_data%def_pts_ker_cntr(5) - input_data%def_pts_ker_cntr(3))

       end if

       call kernel_trapz_int(si,contour_data%ker_int_points(j),contour_data%ker_int_points(j+1),ksw,sw,T)  !! the basis for comparison

       Np = 1

!! the adaptive loop starts here:

       do
          
          panel_no = 2**Np  !! min no of panels = 2
          ds = len/panel_no
          
          allocate(xp(panel_no+1))
          allocate(yp(panel_no+1))
          allocate(zp(panel_no+1))

          xp(1) = REAL(contour_data%ker_int_points(j))
          xp(panel_no+1) = REAL(contour_data%ker_int_points(j+1))
          yp(1) = AIMAG(contour_data%ker_int_points(j))
          yp(panel_no+1) = AIMAG(contour_data%ker_int_points(j+1))
          zp(1) = contour_data%ker_int_points(j)
          zp(panel_no+1) = contour_data%ker_int_points(j+1)

          do i = 1,panel_no-1
      
             if (j >= 1 .AND. j < input_data%num_ker_pts_loop+2) then
                xp(i+1) = xp(i) + ds
                call get_y_alg_int_contour(xp(i+1),yp(i+1), & 
                           REAL(input_data%def_pts_ker_cntr(3)),REAL(input_data%def_pts_ker_cntr(2)-input_data%def_pts_ker_cntr(3)), &
                           AIMAG(input_data%def_pts_ker_cntr(2)-input_data%def_pts_ker_cntr(3)))
             else
                xp(i+1) = xp(i) + ds
                call get_y_alg_int_contour(xp(i+1),yp(i+1), & 
                           REAL(input_data%def_pts_ker_cntr(3)),REAL(input_data%def_pts_ker_cntr(4)-input_data%def_pts_ker_cntr(3)), &
                           AIMAG(input_data%def_pts_ker_cntr(4)-input_data%def_pts_ker_cntr(3)))
                  
             end if

             zp(i+1) = CMPLX(xp(i+1),yp(i+1))
             
          end do

          
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


!! the tolerance check:

          if (1._dpk/(4.**Np-1._dpk)*ABS(intpanel(j)-T) <= ABS(scale*input_data%tol)) then
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


  SUBROUTINE sum_panel_contributions_kernel(Int,totpoints,z,intpanel,Npanel)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Sum up the individual trapezoids for the kernel contour
!! 2. Dump data
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
    
    complex(dpk)  ::  intpanel(total_ker_points-1)
    integer(dpk)  ::   Npanel(total_ker_points-1)
    complex(dpk)    :: Int, z
    integer         :: totpoints
    integer         :: i
    character(60)   :: arg

    Int = (0._dpk,0._dpk)
    do i = 1,total_ker_points-1
       Int = Int + intpanel(i)  !! summing up the individual panels
    end do

!! dump data:

  !  write(arg,"('u=',E12.5,'+',E12.5,'i')") REAL(z),AIMAG(z)

   ! open(10,file='DataDump/Kernel/intkernel.'//arg,form='FORMATTED')
  !  do i = 1, total_ker_points-1
  !     write(10,'(I5,2E20.10)') i,intpanel(i)
  !  end do
  !  close(10)
    
  !  write(*,'(/A22,2X,2F20.10/)') 'Integral value:->', Int

    totpoints = 1
    do i = 1,total_ker_points-1
       totpoints = totpoints + Npanel(i)  !! total quadrature points
    end do

   ! write(*,'(A22,2X,I20/)') 'Quadrature pts:->', totpoints
    

  END SUBROUTINE sum_panel_contributions_kernel


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

    write(arg,"('w=',E12.5,'Np=',I7.7)") REAL(omega_r), N

    open(10,file='DataDump/Kernel_Trace/kernel.'//arg,form='FORMATTED')
    do i = 1, N
       write(10,'(I8,2F20.10,2F30.10)') i,z(i),K(i)
    end do
    close(10)


  END SUBROUTINE kernelplot




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

    write(arg,"('w=',E12.5,'Np=',I7.7)") REAL(omega_r), N

    open(10,file='DataDump/Kernel_Trace/kernel.'//arg,form='FORMATTED')
    do i = 1, N
       write(10,'(I8,2F20.10,2F30.10)') i,z(i),K(i)
    end do
    close(10)


  END SUBROUTINE kernelplot


END SUBROUTINE kernel_integral_units
