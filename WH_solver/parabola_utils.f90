
  SUBROUTINE yfunc_p2(x,y,p,a,e1,e2)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Parabola used for the last loop (down) of integration contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk) :: p
    real(dpk)    :: x, y
    real(dpk)    :: a, e1, e2

    y = -(x-REAL(p))*(1._dpk-(x-REAL(p))*(a-e1-e2)/a**2) + e1


  END SUBROUTINE yfunc_p2


  SUBROUTINE yfunc_p1(x,y,p,e1,e2)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Parabola used for the other loops (up) of integration contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk) :: p
    real(dpk)    :: x, y
    real(dpk)    :: e1, e2

    y = (x-REAL(p))*(1._dpk-(x-REAL(p))/e1) + e2


  END SUBROUTINE yfunc_p1


  SUBROUTINE yfunc_p0(x,y,p,e1,e2)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Parabola used for the other loops (down) of integration contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk) :: p
    real(dpk)    :: x, y
    real(dpk)    :: e1, e2

    y = -(x-REAL(p))*(1._dpk-(x-REAL(p))/e1) + e2


  END SUBROUTINE yfunc_p0


 SUBROUTINE yfunc_p(x,y,p,a,b)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Algebraic "one piece" contour:  see Rienstra, JEM, 2007
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    real(dpk)    :: x, y
    real(dpk)    :: a, b, p

    y = b*(4._dpk*(x-p)/a)/(3._dpk + ((x-p)/a)**4) ! algebraic


  END SUBROUTINE yfunc_p


    SUBROUTINE checkloc(z,switch)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Subroutine to check the location of the instability zeros and pole wrt to the contours
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)    :: z
    real(dpk)       :: zx, zy, zfy
    integer         :: switch

    zx = REAL(z)
    zy = AIMAG(z)

    call yfunc_p(zx,zfy,REAL(szi(3)),REAL(szi(4)-szi(3)),AIMAG(szi(4)-szi(3)))
!!$    zfy = AIMAG(c1)

    if(zy > zfy) then
       switch = 0  !! inside
    else
       switch = 1  !! outside
    end if


  END SUBROUTINE checkloc
