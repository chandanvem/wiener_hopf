Module io_utils
  
  USE input_params 
  USE, intrinsic :: ieee_arithmetic 

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC  :: create_req_dirs, define_input_params, &
             check_NAN 


  CONTAINS 
 
   SUBROUTINE create_directory(dir_name)
    
       character(len=*)   :: dir_name
       integer              :: status
       logical              :: dir_exists

       inquire(file=trim(dir_name), exist=dir_exists)

        if (dir_exists) then
           print *, 'Directory exists: ', trim(dir_name)
        else
           call execute_command_line('mkdir -p ' // trim(dir_name), exitstat=status)
           if (status /= 0) then
             print *, 'Directory creation failed with status:', status
           else
             print *, 'Directory created:', trim(dir_name)
           end if
        end if


   END SUBROUTINE create_directory

   SUBROUTINE create_req_dirs

     call create_directory('./DataDump/Ift')
     call create_directory('./DataDump/Kernel')
     call create_directory('./DataDump/compute_fplus_log')
     call create_directory('./DataDump/compute_IFT_log')
  
   END SUBROUTINE create_req_dirs

   SUBROUTINE define_input_params(input_data)

      type(input_params_t) :: input_data
      real(dpk)            :: PI, delr, deli
      integer              :: i
      character(len=40)    :: input_string      

      !! the basic data file:

      PI = 4._dpk*ATAN(1.)

      input_data%solution_mode = 'hard_duct_mode'

      open(10, file='input.list.p', status='old')

      read(10,*) input_data%near_far_field_mode
      if (trim(input_data%near_far_field_mode) .NE. 'far_field') then
          input_data%near_far_field_mode = 'near_field'
      end if 

      read(10,*) input_data%num_of_streams  !!

      read(10,*) input_data%M1  !! core jet Mach number
      read(10,*) input_data%M2  !! coflow Mach number

      if (input_data%num_of_streams .EQ. 3) then
           read(10,*) input_data%M3  !! ambient flow Mach number
           PRINT '(A, F8.2, ", ", F8.2, ", ", F8.2, A)', &
                     ' (M1, M2, M3) = (',input_data%M1, input_data%M2, input_data%M3, ')'
      else
          PRINT '(A, F8.2, ", ", F8.2, A)', &
                     ' (M1, M2) = (',input_data%M1, input_data%M2, ')'
      end if
 
      read(10,*) input_data%h
      
      read(10,*) input_data%St_flag
      read(10,*) input_data%omega_st

      if (input_data%St_flag .EQ. 1 ) then
          input_data%w0 = PI*input_data%M1*input_data%omega_st  ! Helmholtz number 
      end if

      read(10,*) input_data%kapT  !! sqrt(temp ratio)
      read(10,*) input_data%kap_rho  !! density ratio
      read(10,*) input_data%azim_mode  !! the circumferential mode no


      print*,'azim_mode = ', input_data%azim_mode

      PRINT '(A, F8.2, ", ", F8.2, ", ", F8.2, ", ", F8.2,A)', &
            ' (h, w0, kapT, kap_rho, azim_mode) = (', input_data%h, input_data%w0, input_data%kapT, &
                                                       input_data%kap_rho , input_data%azim_mode, ')'


      read(10,*) input_data%KH_zero_1
      print*, 'initialize:  KH_zero_1 (First instability zero) =', input_data%KH_zero_1

      read(10,*) input_data%mu_plus

      print*, '================== Non-incident vorticity mode ====================='
      print*,''

      print*, 'initialize:  mu for the incident mode=',input_data%mu_plus


      if (input_data%num_of_streams .EQ. 3) then

        read(10,*) input_data%KH_zero_2
        print*, 'initialize:  KH_zero_2 (Second instability zero) =', input_data%KH_zero_2


        read(10,*) input_data%KH_pole_1
        print*, 'initialize:  KH_pole_1 (Instability pole) =', input_data%KH_pole_1
       
      end if

      read(10,*) input_data%offset  !! the offset between the two integration contours
      print*,'initialize:  Offset between the IFT and K split contours=', input_data%offset

      read(10,*) input_data%tol  !! the tolerance of the adaptive contour integration routine
      print*,'initialize:  Tolerance for adaptive contour integration =', input_data%tol
      
      PRINT '(A, E12.4, ", ", E12.4,A)', &
                     ' initialize: (offset, tol) = (', input_data%offset, input_data%tol, ')'

      print*,''
      print*,'','====== Mesh and contour parameters ======',''
      
      read(10,*) input_data%num_ker_pts_loop
       PRINT '(A, I5, A)', &
         'initialize:  Number of kernel points in each loop =', input_data%num_ker_pts_loop

      read(10,*) input_data%theta  !! the stretching parameter for kernel contour meshpoints
      print*,'initialize:  Stretching parameter for kernel contour =', input_data%theta

      read(10,*) input_data%num_IFT_pts_loop
      print*,'initialize:  Number of IFT points in each loop num_IFT_pts_loop =', input_data%num_IFT_pts_loop

      if (trim(input_data%near_far_field_mode) == 'near_field') then
          read(10,*) input_data%Rmin
          read(10,*) input_data%Rmax
          read(10,*) input_data%Zmin
          read(10,*) input_data%Zmax
         
          read(10,*) input_data%Nmeshr
          read(10,*) input_data%Nmeshz

          print*,'initialize:  Dimensions of domain in R = [',input_data%Rmin,input_data%Rmax,']'  
          print*,'initialize:  Dimensions of domain in Z = [',input_data%Zmin,input_data%Zmax,']' 
          print*,'initialize:  Nr x Nz =',input_data%Nmeshr,'x',input_data%Nmeshz
     
      else if (trim(input_data%near_far_field_mode) == 'far_field') then
         read(10,*) input_data%num_phi
         print*,'initialize:  Number of farfield phi locations = ',input_data%num_phi

      end if

      print*,'=============================='

      read(10,*) input_data%asymplim
      print*,'initialize:  Asymptotic limit of bess/dbess = ',input_data%asymplim

      read(10,*) input_data%asymplim1
      print*,'initialize:  Asymptotic limit of hank/dhank = ',input_data%asymplim1

      read(10,*) input_data%vs_param_gamma  !! Default = 1.0 (0,1)
      print*,'initialize:  Vortex shedding parameter gamma = ',input_data%vs_param_gamma 
        
  !!============================
      read(10,*) input_string
      if (trim(input_string) .EQ. 'pressure_mode') then
          input_data%prswitch = 1
      else if (trim(input_string) .EQ. 'potential_mode') then
          input_data%prswitch = 0
      else
          print*, 'initialize: solution mode undefined. Exiting...'
          STOP 
      end if 

      if (input_data%prswitch == 0) then
         print*,'initialize:  Solution in potential mode'
      elseif (input_data%prswitch == 1) then
         print*,'initialize:  Solution in pressure mode '
      end if
  !!============================
     read(10,*) input_string
     if (trim(input_string) .EQ. 'add_incident_mode') then
          input_data%reflswitch = 0
     else
          input_data%reflswitch = 1
     end if 

     if (input_data%reflswitch == 1) then
         print*,'initialize:  Reflection mode: incident mode not added'
     else 
         print*,'initialize:  Reflection mode: incident mode added'
     end if
  !!============================

     read(10,*) input_data%fplus_compute_restart  

     if (len_trim(input_data%solution_mode) > 0) then
          read(10,*) input_data%solution_mode
          if ( (trim(input_data%solution_mode) == 'guided_jet') .OR. & 
                   (trim(input_data%solution_mode) == 'guided_jet_mode')) then
             read(10,*) input_data%s_GJ
             print*, 'initialize: Choosing guided jet mode...'
          else
            print*, 'initialize: Choosing hard duct mode:'
          end if
     else
          print*, 'Error: solution_mode not set or missing in input file.'
     end if

    close(10)
      
    delr = 0._dpk   !! for real omega
    deli = PI/2._dpk  !! for purely imaginary omega

  !! defining the omegas (Helmholtz numbers):

    input_data%omega_r = ABS(input_data%w0)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*delr)  !! real

  !! whether to compute velocity potential or pressure:

    if (input_data%prswitch == 0) then
       if (trim(input_data%near_far_field_mode) =='far_field') then
          print*,''
          print*,'initialize: Far-field is computed ONLY in pressure mode. Exiting...'
          STOP
       else
          print*,''
          print*,'initialize: Computing Potential...'
       end if
    
    else
       print*,''
       print*,'initialize: Computing Pressure...'
    end if

   !now allocate for the contour points:

 END SUBROUTINE define_input_params 


 SUBROUTINE check_NAN(input,is_nan_flag)

   complex(dpk)  :: input
   logical       :: is_nan_flag
   real(dpk)     :: input_real,input_imag, input_abs

   input_real = REAL(input)
   input_imag = AIMAG(input)
   input_abs  = ABS(input_abs)

   if(ieee_is_nan(input_real) .OR. ieee_is_nan(input_imag) .OR. input_abs .GT. 1E8 ) then
        is_nan_flag = .true.
   end if 

 END SUBROUTINE check_NAN

end module io_utils

