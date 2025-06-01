PROGRAM nrcomplex

  USE bessel_utils

  IMPLICIT none
  
  integer, parameter                        :: dpk = kind(1.0d0) !! Double precision kind
  complex(dpk)                              :: input_z
  complex(dpk)                              :: output_bessj, output_bessy
  complex(dpk)                              :: output_bessh
  real(dpk)                                 :: azim_mode, PI

  call readdata                          !! Reads the input file

  output_bessj = bessj(input_z,azim_mode,1)
  output_bessy = bessy(input_z,azim_mode,1)
  output_bessh = hank1(input_z,azim_mode,1)

  print '(A, 2F8.3, A)', ' z = (', REAL(input_z), AIMAG(input_z), ')'
  print '(A, 2F8.3, A)', ' azim mode = ', azim_mode
  print '(A, 2F8.4, A)', ' output bessel j scaled = (', REAL(output_bessj), AIMAG(output_bessj), ')'  
  print '(A, 2F8.4, A)', ' output bessel y scaled = (', REAL(output_bessy), AIMAG(output_bessy), ')'  
  print '(A, 2F8.4, A)', ' output hank 1 scaled = (', REAL(output_bessh), AIMAG(output_bessh), ')'  


CONTAINS

  SUBROUTINE readdata

    open(1,file='input.list')
         read(1,*) input_z                    
         read(1,*) azim_mode                     
    close(1)

  END SUBROUTINE readdata

END PROGRAM nrcomplex
