MODULE fplus_utils

  USE input_params
  USE omp_lib
  USE user_defined_fplus
 
  IMPLICIT NONE

  PUBLIC  :: compute_fplus, read_fplus

  CONTAINS

  SUBROUTINE compute_fplus(input_data,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Evaluate (3.30)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
    
    integer :: i
    integer :: thread_id, num_threads, chunk
    integer :: start_idx, end_idx, file_ID
    character(len=200) :: filename
    integer :: start_time, end_time, clock_rate
    real :: elapsed_time

    type(input_params_t) :: input_data
    type(contour_params_t) :: contour_data

    allocate(input_data%fplusz(contour_data%tot_IFT_pts))  !! fplus is same as \xi^{+}(s) of (3.30)
    input_data%fplusz = (0._dpk,0._dpk)

    !$omp parallel private(i,thread_id,num_threads,chunk,start_idx,end_idx,filename,file_ID, &
    !$omp& start_time,end_time,clock_rate,elapsed_time)

     block 
        thread_id = omp_get_thread_num()
        num_threads = omp_get_num_threads()
        chunk = contour_data%tot_IFT_pts / num_threads
        start_idx = thread_id * chunk + 1
        end_idx = (thread_id + 1) * chunk
        
        if (thread_id == num_threads - 1) end_idx = contour_data%tot_IFT_pts

        write(filename, '("./DataDump/compute_fplus_log/log_", I0, ".out")') thread_id + 1
        file_ID = 10 + thread_id

        call system_clock(start_time,clock_rate)

        do i = start_idx, end_idx
          input_data%fplusz(i) = get_fplus_value(contour_data%iftpoints(i),input_data,contour_data)
          if (i == start_idx) then
              open(file_ID,file=filename,form='FORMATTED',status='UNKNOWN')
          else
              open(file_ID,file=filename,form='FORMATTED',position='APPEND')
          end if
          write(file_ID,'(I5,4E20.10)') i,contour_data%iftpoints(i),input_data%fplusz(i)
        end do
        call system_clock(end_time)
        elapsed_time = real(end_time - start_time)/real(clock_rate)

        write(file_ID, *) 'Elapsed CPU time (seconds):', elapsed_time
        close(file_ID)

     end block
    !$omp end parallel
    
    open(10,file='fplus.out',form='FORMATTED')
    do i = 1, contour_data%tot_IFT_pts
       write(10,'(I5,4E20.10)') i,contour_data%iftpoints(i),input_data%fplusz(i)
     end do
    close(10)


  END SUBROUTINE compute_fplus


  SUBROUTINE read_fplus(input_data,contour_data)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Use the computed F^{+} in case running a restart job
!! 2. Note that restarting is only possible during the "computing F+" portion
!! 3. If a job exits during (R,Z) computation, a "restart" just starts the job
!back at (R=0,Z=0)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    integer          :: i, j, dummy1
    real(dpk)        :: dummy2, dummy3, fplus_real,fplus_imag
  
    type(input_params_t) :: input_data
    type(contour_params_t) :: contour_data


    allocate(input_data%fplusz(contour_data%tot_IFT_pts))  !! fplus is same as \xi^{+}(s) of (3.30)
    input_data%fplusz = (0._dpk,0._dpk)

    open(10,file='fplus.out',form='FORMATTED')
    do i = 1, contour_data%tot_IFT_pts
       read(10,'(I5,4E20.10)') dummy1,dummy2,dummy3,fplus_real,fplus_imag
       input_data%fplusz(i) = fplus_real + CMPLX(0._dpk,1._dpk,kind=dpk)*fplus_imag
    end do
    close(10)

    do i = 1, contour_data%tot_IFT_pts
       print*, 'F+ at:', contour_data%iftpoints(i)
       write(*,'(A22,2X,2F20.10/)') 'F+:->',input_data%fplusz(i)
    end do

    close(10)

  END SUBROUTINE read_fplus


END MODULE fplus_utils

