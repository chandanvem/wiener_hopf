Module user_defined_IFT

  USE omp_lib
  
  USE input_params
  USE user_defined_functions
  USE IFT_integral_utils
  USE user_defined_incident_modes

  IMPLICIT NONE

  PUBLIC :: computeift

  CONTAINS
 
  SUBROUTINE computeift(input_data,contour_data)

  !!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
  !! 1. Compute the pressure/potential at each of the physical points
  !! 2. Compute the constituent parts of the full pressure (potential)
  !! 3. Write the data
  !!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
 
     real(dpk)                                     ::  PI
     complex(dpk), allocatable, dimension(:)       ::  pr_at_IFT_points
     complex(dpk), allocatable, dimension(:,:)     ::  pr_integral
     complex(dpk), allocatable, dimension(:,:,:)   ::  supinstab_pressure
     complex(dpk), allocatable, dimension(:,:)     ::  pressure
     complex(dpk), allocatable, dimension(:,:)     ::  acoupressure,totpressure, &
                                                      instab_pressure1,instab_pressure2,inc_pressure

     integer                                       :: i,j,k,f, IFT_pt_idx
     integer                                       :: switch
     character(60)                                 :: pos, supind
 
     integer                                       :: thread_id, num_threads, chunk
     integer                                       :: start_idx, end_idx, file_ID
     character(len=200)                            :: filename

     type(input_params_t)                          :: input_data
     type(contour_params_t)                        :: contour_data
     
     PI = 4._dpk*ATAN(1.)

  !! allocate the various arrays:
      
     allocate(pr_integral(input_data%Nmeshr,input_data%Nmeshz))
     allocate(pressure(input_data%Nmeshr,input_data%Nmeshz))
     allocate(acoupressure(input_data%Nmeshr,input_data%Nmeshz))
     allocate(instab_pressure1(input_data%Nmeshr,input_data%Nmeshz))
     allocate(inc_pressure(input_data%Nmeshr,input_data%Nmeshz))
     allocate(totpressure(input_data%Nmeshr,input_data%Nmeshz))

    !$omp parallel private(i,j,IFT_pt_idx,thread_id,num_threads,chunk,start_idx,end_idx,filename,file_ID,pr_at_IFT_points)
     block 
        thread_id = omp_get_thread_num()
        num_threads = omp_get_num_threads()
        chunk = input_data%Nmeshz / num_threads
        start_idx = thread_id * chunk + 1
        end_idx = (thread_id + 1) * chunk
        
        if (thread_id == num_threads - 1) end_idx = input_data%Nmeshz

        allocate(pr_at_IFT_points(contour_data%tot_IFT_pts-1))

        write(filename, '("./DataDump/compute_IFT_log/log_", I0, ".out")') thread_id + 1
        file_ID = 10 + thread_id

        do j = start_idx, end_idx
           do i = 1, input_data%Nmeshr  !! the radial (vertical in a 2D plot) direction


              do IFT_pt_idx = 1, contour_data%tot_IFT_pts-1
 
                  call IFT_trapz_int(input_data%R(i),input_data%Z(j),IFT_pt_idx,pr_at_IFT_points(IFT_pt_idx), &
                                                               input_data,contour_data)
               
              end do

              call sum_panel_contribution_IFT(pr_at_IFT_points,pr_integral(i,j),contour_data)
           
              if (input_data%prswitch == 0) then  !! compute velocity potential
                 
                 pressure(i,j) = (input_data%omega_r/(2._dpk*PI*CMPLX(0._dpk,1._dpk,kind=dpk)))*pr_integral(i,j)
                 acoupressure(i,j) = pressure(i,j)  !! acoustic part
                 instab_pressure1(i,j) = residuepot(input_data%R(i),input_data%Z(j),1,input_data)  !! inner instability wave 
                 totpressure(i,j) = acoupressure(i,j) + instab_pressure2(i,j)
                

              else  !! compute pressure

               if (input_data%reflswitch == 1) then

                  pressure(i,j) = (input_data%omega_r*input_data%omega_r/(2._dpk*PI))*pr_integral(i,j)   !! incident wave NOT added

               else

                  pressure(i,j) = (input_data%omega_r*input_data%omega_r/(2._dpk*PI))*pr_integral(i,j) + &
                                   compute_psi_incident(input_data%R(i),input_data%Z(j),input_data) 
                                                  !! the incident wave added

               end if
            
               acoupressure(i,j) = pressure(i,j)  !! acoustic part
               instab_pressure1(i,j) = residuepr(input_data%R(i),input_data%Z(j),input_data)  !! instability wave
               inc_pressure(i,j) = compute_psi_incident(input_data%R(i),input_data%Z(j),input_data)  !! the incident wave
               totpressure(i,j) = acoupressure(i,j) + instab_pressure1(i,j)  !! total part 

            end if

           
            if (j == start_idx) then
                open(file_ID,file=filename,form='FORMATTED',status='UNKNOWN')
            else
                open(file_ID,file=filename,form='FORMATTED',position='APPEND')
            end if
     
             write(file_ID,'(2I5,3E20.10)') i, j, input_data%R(i), input_data%Z(j), REAL(totpressure(i,j)) 
             close(file_ID)
         end do
      end do
     end block
    !$omp end parallel

    !! dump data

    open(1,file='pressure.out',form='UNFORMATTED')
    write(1) input_data%Nmeshz,input_data%Nmeshr,1
    write(1) ((REAL(totpressure(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
    close(1)   

    open(1,file='incident.out',form='UNFORMATTED')
    write(1) input_data%Nmeshz,input_data%Nmeshr,1
    write(1) ((REAL(inc_pressure(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
    close(1)   

    open(1,file='acousticpr.out',form='UNFORMATTED')
    write(1) input_data%Nmeshz,input_data%Nmeshr,1
    write(1) ((REAL(acoupressure(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
    close(1)   

    if (input_data%vortswitch .EQ. 0) then
       open(1,file='instabilitypr1.out',form='UNFORMATTED')
       write(1) input_data%Nmeshz,input_data%Nmeshr,1
       write(1) ((REAL(instab_pressure1(i,j)),j=1,input_data%Nmeshz),i=1,input_data%Nmeshr)
       close(1)   
    end if

    open(1,file='totprdata.out',form='UNFORMATTED')
    write(1) input_data%Nmeshz,input_data%Nmeshr,totpressure
    close(1)   

    open(1,file='acousprdata.out',form='UNFORMATTED')
    write(1) input_data%Nmeshz,input_data%Nmeshr,acoupressure
    close(1)   
    
      
    END SUBROUTINE computeift

end module user_defined_IFT
