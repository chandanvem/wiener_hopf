Module io_utils
  
  IMPLICIT NONE
  ! iii INTEGER, PARAMETER, PUBLIC  :: dp = SELECTED_REAL_KIND(12, 60)

  PRIVATE
  PUBLIC  :: create_directory


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


end module io_utils

