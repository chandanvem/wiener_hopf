Module io_utils
  
  USE input_params 
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC  :: create_req_dirs, define_input_params


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

     END SUBROUTINE create_req_dirs

     SUBROUTINE define_input_params(input_data)

        type(input_params_t) :: input_data
        real(dpk)            :: PI, delr, deli
        integer              :: i
        !! the basic data file:
        open(10, file='input.list.p', status='old')

        read(10,*) input_data%vortswitch  !! 1 = Use incident vorticity mode; 2 = First sup ins mode; Else = acoustic mode

        PRINT '(A, I1, A)', &
                      ' vortswitch = (', input_data%vortswitch, ')'

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
   

        if ((input_data%vortswitch == 1) .OR. (input_data%vortswitch == 2)) then
           read(10,*) input_data%Zo  !! starting point of inc instability (-ve)
        end if

        read(10,*) input_data%St_flag
        read(10,*) input_data%omega_st

        if (input_data%St_flag == 1 ) then
            input_data%w0 = PI*input_data%M1*input_data%omega_st  ! Helmholtz number 
        end if

        read(10,*) input_data%kapT  !! sqrt(temp ratio)
        read(10,*) input_data%kap_rho  !! density ratio
        read(10,*) input_data%azim_mode  !! the circumferential mode no

        PRINT '(A, F8.2, ", ", F8.2, ", ", F8.2, ", ", F8.2,A)', &
              ' (h, w0, kapT, kap_rho, azim_mode) = (', input_data%h, input_data%w0, input_data%kapT, &
                                                         input_data%kap_rho , input_data%azim_mode, ')'


        read(10,*) input_data%KH_zero_1
        print*, 'initialize:  KH_zero_1 (First instability zero) =', input_data%KH_zero_1
 
        if (input_data%num_of_streams .EQ. 3) then

          read(10,*) input_data%KH_zero_2
          print*, 'initialize:  KH_zero_2 (Second instability zero) =', input_data%KH_zero_2


          read(10,*) input_data%KH_pole_1
          print*, 'initialize:  KH_pole_1 (Instability pole) =', input_data%KH_pole_1
         
        end if

        if ((input_data%vortswitch == 1) .OR. (input_data%vortswitch == 2)) then
                   read(10,*)input_data%mu0  !! upstream inc acoustic mode
                   print*, 'initialize:  Incident Vorticity mode activated'
        else
                   read(10,*) input_data%mu_plus
                   print*, '================== Non-incident vorticity mode ====================='
                   print*,''
                   print*, 'initialize:  mu for the incident mode =',input_data%mu_plus
        end if

        read(10,*) input_data%offset  !! the offset between the two integration contours
        print*,'initialize:  Offset between the IFT and K split contours=', input_data%offset

        read(10,*) input_data%tol  !! the tolerance of the adaptive contour integration routine
        print*,'initialize:  Tolerance for adaptive contour integration =', input_data%tol
        
        PRINT '(A, E12.4, ", ", E12.4,A)', &
                       ' initialize: (offset, tol) = (', input_data%offset, input_data%tol, ')'


        read(10,*) input_data%num_zeros_s1_s2  !! zeros needed to be removed from the kernel
        PRINT '(A, I3, A)', &
                       'initialize:  Number of zeros to be removed from kernel ', input_data%num_zeros_s1_s2


        read(10,*) input_data%num_poles_s1_s2  !! poles needed to be removed from the kernel
        PRINT '(A, I3, A)', &
                       'initialize:  Number of poles to be removed from kernel ', input_data%num_poles_s1_s2


        read(10,*) input_data%num_sup_zeros  !! supersonic zeros
        PRINT '(A, I3, A)', &
                       'initialize:  Number of supersonic zeros ', input_data%num_sup_zeros


        read(10,*) input_data%num_sup_poles  !! supersonic poles
        PRINT '(A, I3, A)', &
                       'initialize:  Number of supersonic poles ', input_data%num_sup_poles

        if (input_data%num_sup_zeros > 0) then
            allocate(input_data%sup_zeros_list(input_data%num_sup_zeros))
            do i = 1, input_data%num_sup_zeros
              read(10,*) input_data%sup_zeros_list(i)
              print*, input_data%sup_zeros_list(i)
            end do
        end if

         if (input_data%num_sup_poles > 0) then
            allocate(input_data%sup_poles_list(input_data%num_sup_poles))
            do i = 1, input_data%num_sup_poles
              !print*, i1
              read(10,*) input_data%sup_poles_list(i)
              print*, input_data%sup_poles_list(i)
            end do
        end if

        if (input_data%vortswitch == 2) then
           if (input_data%num_sup_poles == 0) then
                      print*, "initialize:  For vortswitch mode 2, at least one upstream supersonic pole needed! Exiting..."
                      STOP
           end if
        end if

        print*,''
        print*,'','====== Mesh and contour parameters ======',''
        
        read(10,*) input_data%num_ker_pts_loop
         PRINT '(A, I5, A)', &
           'initialize:  Number of kernel points in each loop =', input_data%num_ker_pts_loop

        read(10,*) input_data%theta  !! the stretching parameter for kernel contour meshpoints
        print*,'initialize:  Stretching parameter for kernel contour =', input_data%theta

        read(10,*) input_data%num_IFT_pts_loop
        print*,'initialize:  Number of IFT points in each loop num_IFT_pts_loop =', input_data%num_IFT_pts_loop

        read(10,*) input_data%Rmin
        read(10,*) input_data%Rmax
        read(10,*) input_data%Zmin
        read(10,*) input_data%Zmax
       
        read(10,*) input_data%Nmeshr
        read(10,*) input_data%Nmeshz

        print*,'initialize:  Dimensions of domain in R = [',input_data%Rmin,input_data%Rmax,']'  
        print*,'initialize:  Dimensions of domain in Z = [',input_data%Zmin,input_data%Zmax,']' 
        print*,'initialize:  Nr x Nz =',input_data%Nmeshr,'x',input_data%Nmeshz

        print*,'=============================='

        read(10,*) input_data%asymplim
        print*,'initialize:  Asymptotic limit of bess/dbess = ',input_data%asymplim

        read(10,*) input_data%asymplim1
        print*,'initialize:  Asymptotic limit of hank/dhank = ',input_data%asymplim1

        read(10,*) input_data%vs_param_gamma  !! Default = 1.0 (0,1)
        print*,'initialize:  Vortex shedding parameter gamma = ',input_data%vs_param_gamma 
          
    !!============================
        read(10,*) input_data%prswitch  !! 0 = Potential; 1 = Pressure

        if (input_data%prswitch == 0) then
           print*,'initialize:  Solution in potential mode'
        elseif (input_data%prswitch == 1) then
           print*,'initialize:  Solution in pressure mode '
        end if
    !!============================
        read(10,*) input_data%reflswitch  !! 1 = Reflection mode: incident mode not added
        if (input_data%reflswitch == 1) then
           print*,'initialize:  Reflection mode: incident mode not added'
        else 
           print*,'initialize:  Reflection mode: incident mode added'
        end if
    !!============================
        read(10,*) input_data%farswitch  !! 1 = Far-field mode: compute directivity; 2 = 1 + nearfield of sup zeros in polar mode

        if ((input_data%farswitch == 1) .OR. (input_data%farswitch == 2)) then
           read(10,*) input_data%Nphi
           allocate(input_data%cnanglei(input_data%num_sup_zeros + 2))
        end if

        if ((input_data%farswitch == 2) .AND. input_data%num_sup_zeros <= 0) then
           print*, "For farswitch mode 2, non-zero number of supersonic zero needed! Exiting..."
           STOP
        end if

    !!=====================================================

        read(10,*) input_data%restart  !! restart status: 0 = fresh job, i.e., no "fplus_part.out" exists

      close(10)
        
      PI = 4._dpk*ATAN(1.)
      delr = 0._dpk   !! for real omega
      deli = PI/2._dpk  !! for purely imaginary omega

    !! defining the omegas (Helmholtz numbers):

      input_data%omega_r = ABS(input_data%w0)*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*delr)  !! real

    !! whether to compute velocity potential or pressure:

        if (input_data%prswitch == 0) then
           if ((input_data%farswitch == 1) .OR. (input_data%farswitch == 2)) then
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

     !! now allocate for the contour points:

   END SUBROUTINE define_input_params 



end module io_utils

