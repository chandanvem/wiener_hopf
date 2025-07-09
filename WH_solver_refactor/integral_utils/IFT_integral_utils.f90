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

          write(pos,"('R=',F10.5,'Z=',F10.5)") R(i),Z(j)

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
          write(supind,"(I1)") k
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


 
