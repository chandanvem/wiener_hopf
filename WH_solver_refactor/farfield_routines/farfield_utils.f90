Module farfield_utils
  
  USE input_params 
  USE user_defined_fplus

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC  :: init_farfield, compute_farfield

  CONTAINS 

  SUBROUTINE init_farfield(input_data,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Initializes the far-field computations
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    type(input_params_t)          :: input_data
    type(contour_params_t)        :: contour_data
    integer                       :: i
    real(dpk)                     :: PI
    real(dpk)                     :: phi_start, phi_end
    

    PI = 4._dpk*ATAN(1.)

    allocate(input_data%phi_list(input_data%num_phi))
    allocate(input_data%Dmn_list(input_data%num_phi))

    ! 0 < phi < 180

    phi_start =   0.5_dpk*PI/180._dpk
    phi_end   = 179.5_dpk*PI/180._dpk

    ! Generate the directivity mesh

    do j = 1,input_data%num_phi
       input_data%phi_list(j) = phi_start + (j-1)*(phi_start-phi_end)/(input_data%num_phi-1)
    end do


  END SUBROUTINE init_farfield

  SUBROUTINE compute_farfield(input_data,contour_data)

    type(input_params_t)          :: input_data
    type(contour_params_t)        :: contour_data
    integer                       :: phi_list_index 
    real(dpk)                     :: PI

    PI = 4._dpk*ATAN(1.)
    jj = 1

    if (restart .NE. 0) then

       open(10,file='directivity_temp.out',form='FORMATTED')
       do i = 1, num_phi
          read(10,'(I5,2F20.15)',END = 100) k, Dmn(i)
       END do
       close(10)

      do i = 1, k

          write(*,'(/A10,2X,F15.10,A1)'),'angle :: ', phi(i)*180._dpk/PI,':'
          write(*,'(A10,2X,2F15.10/)'),'Dmn =', Dmn(i)
          if (mod(i,10) == 0)         write(*,'(A64)'),'----------------------------------------------------------------'

          jj = k+1

       END do

    END if

    do phi_list_index = jj,input_data%num_phi
       
       write(*,'(/A10,2X,F15.10,A1)'),'angle in degrees :: ',input_data%phi_list(phi_list_index)*180._dpk/PI,':'

       call compute_directivity(input_data%phi_list(phi_list_index),input_data%Dmn_list(phi_list_index))
       
       write(*,'(A10,2X,2F15.10/)'),'Dmn =',input_data%Dmn_list(phi_list_index)

       if (phi_list_index == 1) then
          open(10,file='directivity_temp.out',form='FORMATTED',status='UNKNOWN')
       else
          open(10,file='directivity_temp.out',form='FORMATTED',position='APPEND')
       END if
       write(10,'(I5,2F20.15)') phi_list_index,input_data%Dmn_list(phi_list_index)
       close(10)

       if (mod(phi_list_index,10) == 0) write(*,'(A64)'),'----------------------------------------------------------------'
    END do

    open(10,file='directivity.out',form='FORMATTED')
    do i = 1, num_phi
       if (phi(i) > cnangle) write(10,'(F8.2,2F20.15)') input_data%phi_list(phi_index_list)*180._dpk/PI,input_data%Dmn_list(phi_index_list)
    END do
    close(10)

  END SUBROUTINE compute_farfield


 END MODULE farfield_utils



    ! Find the instability wave cone angles
    
!    cnanglei(1) = ATAN(- (AIMAG(sz1))/(AIMAG(SQRT(input_data%kapT**2*(1._dpk - sz1*input_data%M2)**2 -sz1**2))))
!    cnanglei(2) = ATAN(- (AIMAG(sz2))/(AIMAG(SQRT(input_data%kapT**2*(1._dpk - sz2*input_data%M2)**2 -sz2**2))))
!
!    do j = 3, 2+Nzsp
!       cnanglei(j) = ATAN(- (AIMAG(szsp(j-2)))/(AIMAG(SQRT(input_data%kapT**2*(1._dpk -szsp(j-2)*input_data%M2)**2 -&
!            szsp(j-2)**2))))
!    END do

!!!$    print*,'coneangle1=',cnangle1,'coneangle2=',cnangle2
!
!    cnangle = cnanglei(1)
!    do j = 2, Nzsp+1
!       if (cnanglei(j) >= cnangle) then
!          cnangle = cnanglei(j)
!       END if
!    END do
!
!!!$    print*,'coneangle=',coneangle*180./PI
!
!    open(10,file='coneangles.out',form='FORMATTED')
!    write(10,'(A5,F20.5)') 'max:', cnangle*180._dpk/PI
!    do j = 1, Nzsp+2
!       write(10,'(I3,F20.5)') j, cnanglei(j)*180._dpk/PI
!    END do
!    close(10)
!
!    if (farswitch == 2) call meshgrid


