Module io_utils
  
  USE input_params 
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC  :: create_directory, define_input_params


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


     SUBROUTINE define_input_params(input_data)

        type(input_params_t) :: input_data

        !! the basic data file:
        open(10, file='input.list.p', status='old')

        read(10,*) input_data%vortswitch  !! 1 = Use incident vorticity mode; 2 = First sup ins mode; Else = acoustic mode

        PRINT '(A, I3, A)', &
                      ' vortswitch = (', input_data%vortswitch, ')'

        END SUBROUTINE define_input_params 



end module io_utils

