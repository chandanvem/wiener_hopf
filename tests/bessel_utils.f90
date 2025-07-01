Module bessel_utils

  USE Complex_Bessel

  IMPLICIT NONE
  ! iii INTEGER, PARAMETER, PUBLIC  :: dp = SELECTED_REAL_KIND(12, 60)

  PRIVATE
  PUBLIC  :: bessj, bessy, hank1, hank2, dbessj, dbessy, dhank1, dhank2

  CONTAINS 

    FUNCTION bessj(z,order,s)

        complex(dp)        :: z, bessj, psi, bj, by, cwrk
        real(dp)           :: order, q, expo, pow
        integer            :: s, nz, ifail, unscaled_op_flag = 1, scaled_op_flag = 2
        COMPLEX (dp)       :: bessj_op(20), bessy_op(20)
        REAL(dp)           :: PI

        PI = 4._dp*ATAN(1._dp)

    
        if (order < 0.) then
              if (s == 0) then
              !   call S17DEF(-order,z,1,'U',bj,nz,ifail)
                  call cbesj(z, -order, unscaled_op_flag, 1,  bessj_op,  nz, ifail)

              !   call S17DCF(-order,z,1,'U',by,nz,cwrk,ifail)
                  call cbesy(z, -order, unscaled_op_flag, 1,  bessy_op,  nz, ifail)

              !   bessj = bj*COS(-PI*nu) - by*SIN(-PI*nu)              
                  bessj = bessj_op(1)*COS(-PI*order) - bessy_op(1)*SIN(-PI*order)

              else
              !   call S17DEF(-order,z,1,'S',bj,nz,ifail)
                  call cbesj(z, -order, scaled_op_flag, 1,  bessj_op,  nz, ifail)

              !   call S17DCF(-order,z,1,'S',by,nz,cwrk,ifail)
                  call cbesy(z, -order, scaled_op_flag, 1,  bessy_op,  nz, ifail)
           
              !   bessj = bj*COS(-PI*nu) - by*SIN(-PI*nu)              
                  bessj = bessj_op(1)*COS(-PI*order) - bessy_op(1)*SIN(-PI*order)
              end if
        else
           if (s == 0) then
            !  call S17DEF(order,z,1,'U',bessj,nz,ifail)
               call cbesj(z, order, unscaled_op_flag, 1,  bessj_op,  nz, ifail)
           else
            !  call S17DEF(nu,z,1,'S',bessj,nz,ifail)
               call cbesj(z, order, scaled_op_flag, 1,  bessj_op,  nz, ifail)
           end if

           bessj = bessj_op(1)

         end if


    END FUNCTION bessj


    FUNCTION bessy(z,order,s)


        complex(dp)      :: z, bessy, psi, bj, by, cwrk
        real(dp)         :: order, q, expo, pow
        integer           :: s, nz, ifail, unscaled_op_flag = 1, scaled_op_flag = 2
        COMPLEX (dp)      :: bessj_op(20), bessy_op(20)
        REAL(dp)                    :: PI

        PI = 4._dp*ATAN(1._dp)


        if (order < 0.) then
              if (s == 0) then
                       !  call S17DEF(-order,z,1,'U',bj,nz,ifail)
                       call cbesj(z, -order, unscaled_op_flag, 1,  bessj_op,  nz, ifail)

                       !  call S17DCF(-order,z,1,'U',by,nz,cwrk,ifail)
                       call cbesy(z, -order, unscaled_op_flag, 1,  bessy_op,  nz, ifail)

                       !  bessy = by*COS(-PI*order) + bj*SIN(-PI*order)
                       bessy = bessy_op(1)*COS(-PI*order) + bessj_op(1)*SIN(-PI*order)

              else
                       !  call S17DEF(-order,z,1,'S',bj,nz,ifail)
                       call cbesj(z, -order, scaled_op_flag, 1,  bessj_op,  nz, ifail)

                       !  call S17DCF(-order,z,1,'S',by,nz,cwrk,ifail)
                       call cbesy(z, -order, scaled_op_flag, 1,  bessy_op,  nz, ifail)

                       !  bessy = by*COS(-PI*order) + bj*SIN(-PI*order)
                       bessy = bessy_op(1)*COS(-PI*order) + bessj_op(1)*SIN(-PI*order)

              end if

           else
              if (s == 0) then
                       ! call S17DCF(order,z,1,'U',bessy,nz,cwrk,ifail)
                       call cbesy(z, order, unscaled_op_flag, 1,  bessy_op,  nz, ifail)
              else
                       !  call S17DCF(order,z,1,'S',bessy,nz,cwrk,ifail)
                       call cbesy(z, order, scaled_op_flag, 1,  bessy_op,  nz, ifail)
              end if

              bessy = bessy_op(1)

           end if
        
      END FUNCTION bessy

    FUNCTION hank1(z,order,s)

        complex(dp)      :: z, hank1, psi, bj, by, cwrk
        real(dp)         :: order, q, expo, pow
        integer           :: s, nz, ifail, unscaled_op_flag = 1, scaled_op_flag = 2
        integer           :: hank_kind_index = 1
        COMPLEX (dp)      :: bessj_op(20), bessy_op(20), bessh_op(20)
        REAL(dp)                    :: PI

        PI = 4._dp*ATAN(1._dp)


        if (order < 0.) then
            if (s == 0) then
                     !  call S17DLF(1,-order,z,1,'U',h1,nz,ifail)
                     !  hank1 = EXP(-order*PI*CMPLX(0._dp,1._dp,kind=dp))*h1
                     call cbesh(z, -order, unscaled_op_flag,hank_kind_index, 1,  bessh_op,  nz, ifail)
                     hank1 = EXP(-order*PI*CMPLX(0._dp,1._dp,kind=dp))*bessh_op(1)
            else
                     !  call S17DLF(1,-order,z,1,'S',h1,nz,ifail)
                     !  hank1 = EXP(-order*PI*CMPLX(0._dp,1._dp,kind=dp))*h1
                     call cbesh(z, -order, scaled_op_flag,hank_kind_index, 1,  bessh_op,  nz, ifail)
                     hank1 = EXP(-order*PI*CMPLX(0._dp,1._dp,kind=dp))*bessh_op(1)
            end if
         else
            if (s == 0) then
                     !  call S17DLF(1,order,z,1,'U',hank1,nz,ifail)
                     call cbesh(z, order, unscaled_op_flag,hank_kind_index, 1,  bessh_op,  nz, ifail)

            else
                     !  call S17DLF(1,order,z,1,'S',hank1,nz,ifail)
                     call cbesh(z, order, scaled_op_flag,hank_kind_index, 1,  bessh_op,  nz, ifail)
            end if

             hank1 = bessh_op(1)

         end if

       
        
    END FUNCTION hank1

    FUNCTION hank2(z,order,s)

        complex(dp)      :: z, hank2, psi, bj, by, cwrk
        real(dp)         :: order, q, expo, pow
        integer           :: s, nz, ifail, unscaled_op_flag = 1, scaled_op_flag = 2
        integer           :: hank_kind_index = 2
        COMPLEX (dp)      :: bessj_op(20), bessy_op(20), bessh_op(20)
        REAL(dp)                    :: PI

        PI = 4._dp*ATAN(1._dp)

        if (order < 0.) then
            if (s == 0) then
                     !  call S17DLF(1,-order,z,1,'U',h1,nz,ifail)
                     !  hank2 = EXP(-order*PI*CMPLX(0._dp,1._dp,kind=dp))*h1
                     call cbesh(z, -order, unscaled_op_flag,hank_kind_index, 1,  bessh_op,  nz, ifail)
                     hank2 = EXP(-order*PI*CMPLX(0._dp,1._dp,kind=dp))*bessh_op(1)
            else
                     !  call S17DLF(1,-order,z,1,'S',h1,nz,ifail)
                     !  hank2 = EXP(-order*PI*CMPLX(0._dp,1._dp,kind=dp))*h1
                     call cbesh(z, -order, scaled_op_flag,hank_kind_index, 1,  bessh_op,  nz, ifail)
                     hank2 = EXP(-order*PI*CMPLX(0._dp,1._dp,kind=dp))*bessh_op(1)
            end if
         else
            if (s == 0) then
                     !  call S17DLF(1,order,z,1,'U',hank2,nz,ifail)
                     call cbesh(z, order, unscaled_op_flag,hank_kind_index, 1,  bessh_op,  nz, ifail)

            else
                     !  call S17DLF(1,order,z,1,'S',hank2,nz,ifail)
                     call cbesh(z, order, scaled_op_flag,hank_kind_index, 1,  bessh_op,  nz, ifail)
            end if

             hank2 = bessh_op(1)

         end if

       
        
    END FUNCTION hank2


    FUNCTION dbessj(z,nu,s)

        complex(dp)    :: z, dbessj, psi
        real(dp)       :: nu, q, expo, pow
        integer         :: s

        dbessj = 0.5*(bessj(z,nu-1,s) - bessj(z,nu+1,s))


    END FUNCTION dbessj


    FUNCTION dbessy(z,nu,s)

        complex(dp)    :: z, dbessy, psi
        real(dp)       :: nu, q, expo, pow
        integer         :: s

        dbessy = 0.5*(bessy(z,nu-1,s) - bessy(z,nu+1,s))


    END FUNCTION dbessy


    FUNCTION dhank1(z,nu,s)

      complex(dp)    :: z, dhank1, psi
      real(dp)       :: nu, q, expo, pow
      integer         :: s

      
      dhank1 = 0.5*(hank1(z,nu-1,s) - hank1(z,nu+1,s))

    END FUNCTION dhank1


    FUNCTION dhank2(z,nu,s)

        complex(dp)    :: z, dhank2, psi
        real(dp)       :: nu, q, expo, pow
        integer         :: s

        dhank2 = 0.5*(hank2(z,nu-1,s) - hank2(z,nu+1,s))


    END FUNCTION dhank2

END Module bessel_utils
