Module IFT_integral_utils

  USE input_params
  USE contour_generate_utils
  USE contour_init_utils
  USE user_defined_functions

  IMPLICIT NONE

  PUBLIC ::  sum_panel_contribution_IFT, IFT_trapz_int 

  CONTAINS

!  SUBROUTINE computeift(input_data,contour_data)
!
!  !!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!  !! 1. Compute the pressure/potential at each of the physical points
!  !! 2. Compute the constituent parts of the full pressure (potential)
!  !! 3. Write the data
!  !!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!      real(dpk)                                     :: PI
!      complex(dpk), allocatable, dimension(:)       ::  pr
!      complex(dpk), allocatable, dimension(:,:)     ::  prsum
!      complex(dpk), allocatable, dimension(:,:,:)   ::  supinspressure
!      complex(dpk), allocatable, dimension(:,:)     ::  pressure
!      complex(dpk), allocatable, dimension(:,:)     ::  acoupressure,totpressure, &
!                                                       inspressure1,inspressure2,initpressure
!
!      integer                                       :: i,j,k,f, IFT_pt_idx
!      integer                                       :: switch
!      character(60)                                 :: pos, supind
!
!      type(input_params_t)                          :: input_data
!      type(contour_params_t)                        :: contour_data
!     
!      PI = 4._dpk*ATAN(1.)
!
!  !! allocate the various arrays:
!      
!      allocate(pr(contour_data%tot_IFT_pts-1))
!      allocate(prsum(input_data%Nmeshr,input_data%Nmeshz))
!      allocate(pressure(input_data%Nmeshr,input_data%Nmeshz))
!      allocate(acoupressure(input_data%Nmeshr,input_data%Nmeshz))
!      allocate(inspressure1(input_data%Nmeshr,input_data%Nmeshz))
!
!      if (input_data%num_of_streams == 3) then
!         allocate(inspressure2(input_data%Nmeshr,input_data%Nmeshz))
!      end if
!
!      allocate(totpressure(input_data%Nmeshr,input_data%Nmeshz))
!      allocate(initpressure(input_data%Nmeshr,input_data%Nmeshz))
!
!      if (input_data%num_sup_zeros > 0) then
!         allocate(supinspressure(input_data%num_sup_zeros,input_data%Nmeshr,input_data%Nmeshz))
!      end if
!
!  !! the basic do loop:
!
!      do j = 1, input_data%Nmeshz  !! the axial direction
!         do i = 1, input_data%Nmeshr  !! the radial (vertical in a 2D plot) direction
!
!            print*,'Pressure at Z=',input_data%Z(j),'and R=',input_data%R(i)
!
!  !! different parts of the physical domain has separate formulae, so check:
!
!            write(pos,"('R=',F10.5,'Z=',F10.5)") input_data%R(i),input_data%Z(j)
!
!            do IFT_pt_idx = 1, contour_data%tot_IFT_pts-1
!
!               call IFT_trapz_int(input_data%R(i),input_data%Z(j),IFT_pt_idx,pr(IFT_pt_idx), &
!                                                             input_data,contour_data)
!               
!            end do
!
!            open(10,file='DataDump/Ift/intift.'//pos,form='FORMATTED')
!            do k = 1,contour_data%tot_IFT_pts-1
!               write(10,'(I5,2E20.10)') k,pr(k)
!            end do
!            close(10)
!
!            call sum_panel_contribution_IFT(pr,prsum(i,j),contour_data)
!            
!            if (input_data%prswitch == 0) then  !! compute velocity potential
!               
!               pressure(i,j) = (input_data%omega_r/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))*prsum(i,j)
!
!               acoupressure(i,j) = pressure(i,j)  !! acoustic part
!
!               if (input_data%vortswitch .EQ. 0) then
!                  inspressure1(i,j) = residuepot(input_data%R(i),input_data%Z(j),1,input_data)  !! inner instability wave 
!               end if
!
!               if (input_data%num_of_streams .EQ. 3) then
!                   inspressure2(i,j) = residuepot(input_data%R(i),input_data%Z(j),2,input_data)  !! outer instability wave
!               end if
!
!               if (input_data%vortswitch .EQ. 0  .AND. input_data%num_sup_zeros .GE. 1) then
!                  do k = 1, input_data%num_sup_zeros
!                     supinspressure(k,i,j) = residuepot(input_data%R(i),input_data%Z(j),k+2,input_data)
!                  end do
!               end if
!               
!               if (input_data%vortswitch .EQ. 0) then
!                  totpressure(i,j) = acoupressure(i,j) + inspressure1(i,j)  !! total part 
!                  if (input_data%num_of_streams .EQ. 3) then
!                    totpressure(i,j) =  totpressure(i,j)  + inspressure2(i,j) 
!                  end if
!               else
!                  totpressure(i,j) = acoupressure(i,j) + inspressure2(i,j)
!               end if
!
!               if (input_data%vortswitch .EQ. 0 .AND. input_data%num_sup_zeros .GE. 1) then
!                  do k = 1, input_data%num_sup_zeros
!                     totpressure(i,j) = totpressure(i,j) + supinspressure(k,i,j)
!                  end do
!               end if
!
!               print*,'pressure:',REAL(totpressure(i,j))
!
!            else  !! compute pressure
!
!               if (input_data%reflswitch == 1) then
!
!                  pressure(i,j) = (input_data%omega_r*input_data%omega_r/(2._dpk*PI))*prsum(i,j)   !! incident wave NOT added
!
!               else
!
!                  pressure(i,j) = input_data%omega_r*input_data%omega_r/(2._dpk*PI)*prsum(i,j) + &
!                                   compute_psi_incident(input_data%R(i),input_data%Z(j),input_data) 
!                                                  !! the incident wave added
!
!               end if
!            
!               acoupressure(i,j) = pressure(i,j)  !! acoustic part
!
!               if (input_data%vortswitch .EQ. 0) then
!                  inspressure1(i,j) = residuepr(input_data%R(i),input_data%Z(j),1,input_data)  !! inner instability wave
!               end if
!
!               if (input_data%num_of_streams .EQ. 3) then
!                    inspressure2(i,j) = residuepr(input_data%R(i),input_data%Z(j),2,input_data)  !! outer instability wave
!               end if
!
!               if (input_data%vortswitch .EQ. 0 .AND. input_data%num_sup_zeros .GE. 1) then
!                  do k = 1, input_data%num_sup_zeros
!                     supinspressure(k,i,j) = residuepr(input_data%R(i),input_data%Z(j),k+2,input_data)
!                  end do
!               end if
!   
!               if (input_data%vortswitch .EQ. 0) then
!                  totpressure(i,j) = acoupressure(i,j) + inspressure1(i,j)  !! total part 
!                  if (input_data%num_of_streams .EQ. 3) then
!                    totpressure(i,j) =  totpressure(i,j)  + inspressure2(i,j) 
!                  end if
!               else
!                  totpressure(i,j) = acoupressure(i,j) + inspressure2(i,j)
!               end if
!
!               if (input_data%vortswitch .EQ. 0 .AND. input_data%num_sup_zeros .GE. 1) then
!                  do k = 1, input_data%num_sup_zeros
!                     totpressure(i,j) = totpressure(i,j) + supinspressure(k,i,j)
!                  end do
!               end if
!
!               print*,'pressure:',REAL(totpressure(i,j))
!
!            end if
!
!            initpressure(i,j) = compute_psi_incident(input_data%R(i),input_data%Z(j),input_data)  !! the incident wave
!            
!            print*,''
!            
!         end do
!      end do
!
!  !! dump data
!      
!      open(1,file='pressure.out',form='UNFORMATTED')
!      write(1) input_data%Nmeshz,input_data%Nmeshr,1
!      write(1) ((REAL(totpressure(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
!      close(1)   
!
!      open(1,file='incident.out',form='UNFORMATTED')
!      write(1) input_data%Nmeshz,input_data%Nmeshr,1
!      write(1) ((REAL(initpressure(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
!      close(1)   
!
!      open(1,file='acousticpr.out',form='UNFORMATTED')
!      write(1) input_data%Nmeshz,input_data%Nmeshr,1
!      write(1) ((REAL(acoupressure(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
!      close(1)   
!
!      if (input_data%vortswitch .EQ. 0) then
!         open(1,file='instabilitypr1.out',form='UNFORMATTED')
!         write(1) input_data%Nmeshz,input_data%Nmeshr,1
!         write(1) ((REAL(inspressure1(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
!         close(1)   
!      end if
!   
!
!      if (input_data%num_of_streams .EQ. 0) then
!          open(1,file='instabilitypr2.out',form='UNFORMATTED')
!          write(1) input_data%Nmeshz,input_data%Nmeshr,1
!          write(1) ((REAL(inspressure2(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
!          close(1)  
!      end if
!
!      if (input_data%vortswitch .EQ. 0 .AND. input_data%num_sup_zeros .GE. 1) then
!         do k = 1, input_data%num_sup_zeros
!            write(supind,"(I1)") k
!            open(1,file='supinstabilitypr.out.'//supind,form='UNFORMATTED')
!            write(1) input_data%Nmeshz,input_data%Nmeshr,1
!            write(1) ((REAL(supinspressure(k,i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
!            close(1)  
!         end do
!      end if
!
!      open(1,file='totprdata.out',form='UNFORMATTED')
!      write(1) input_data%Nmeshz,input_data%Nmeshr,totpressure
!      close(1)   
!
!      open(1,file='acousprdata.out',form='UNFORMATTED')
!      write(1) input_data%Nmeshz,input_data%Nmeshr,acoupressure
!      close(1)   
!      
!      
!    END SUBROUTINE computeift

SUBROUTINE IFT_trapz_int(r,z,IFT_pt_idx,integral,input_data,contour_data)
             
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. The trapezoidal rule for the IFT contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)                   :: r, z
    complex(dpk)                :: ds, integral
    complex(dpk)                :: f1, f2
    integer                     :: IFT_pt_idx, ss

    type(input_params_t)        :: input_data
    type(contour_params_t)     :: contour_data

    if (input_data%prswitch == 0) then

       f1 = integrand_IFT_pot(r,z,IFT_pt_idx,input_data,contour_data)
       f2 = integrand_IFT_pot(r,z,IFT_pt_idx+1,input_data,contour_data)

    else

       f1 = integrand_IFT_pr(r,z,IFT_pt_idx,input_data,contour_data)
       f2 = integrand_IFT_pr(r,z,IFT_pt_idx+1,input_data,contour_data)

    end if

    ds = contour_data%iftpoints(IFT_pt_idx+1) - contour_data%iftpoints(IFT_pt_idx)

    integral = (ds/2._dpk)*(f1+f2)

  END SUBROUTINE IFT_trapz_int

 SUBROUTINE sum_panel_contribution_IFT(panel,integral,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Sum up the individual trapezoids for the IFT contour
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    type(contour_params_t)                     :: contour_data

    complex(dpk)                               :: integral
    complex(dpk),dimension(contour_data%tot_IFT_pts-1)      :: panel
    integer                                    :: i

    integral = (0._dpk,0._dpk)
    do i = 1,contour_data%tot_IFT_pts-1
       integral = integral + panel(i)
    end do


  END SUBROUTINE sum_panel_contribution_IFT

END MODULE IFT_integral_utils

