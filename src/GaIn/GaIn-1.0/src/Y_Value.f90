
!> @file Y_Value.f90
!!
!! Defines utility routine for computing value
!! of a solid gaussian basis elements at point r.
!!
!! Author: I. Duchemin July 2015
!!

!> Value of a Solid real spherical harmonic centered in (0,0,0) at point r
!! 
recursive function Y_Value(a,l,m,r)
 
  implicit none
 
  real(kind=8), intent(in):: a                 !< a exponent for Ylm
  integer     , intent(in):: l                 !< l for Ylm, l must be < 6
  integer     , intent(in):: m                 !< m for Ylm
  real(kind=8), intent(in),dimension(3)   :: r !< relative position in space
  
  ! return value
  real(kind=8)                            :: Y_Value
 
  ! constants
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
 
  ! local variables
  real(kind=8) :: x
  real(kind=8) :: y
  real(kind=8) :: z
  real(kind=8) :: r2
 
  ! set coordinates
  x=r(1)
  y=r(2)
  z=r(3)
  r2=x**2+y**2+z**2
 
  ! selection on l
  select case(l)
    case (0)
      ! selection on m: l=0
      select case(m)
        case (0)
          Y_Value=(1.0d0/2.0d0/sqrt(PI))*exp(-a*r2)
        case default
          print*,'Error: Y_Value not implemented for l=',l,'m=',m
          stop
      end select
    case (1)
      ! selection on m: l=1
      select case(m)
        case (-1)
          Y_Value=(1.0d0/2.0d0*sqrt(3.0d0)/sqrt(PI)*y)*exp(-a*r2)
        case (0)
          Y_Value=(1.0d0/2.0d0*sqrt(3.0d0)/sqrt(PI)*z)*exp(-a*r2)
        case (1)
          Y_Value=(1.0d0/2.0d0*sqrt(3.0d0)/sqrt(PI)*x)*exp(-a*r2)
        case default
          print*,'Error: Y_Value not implemented for l=',l,'m=',m
          stop
      end select
    case (2)
      ! selection on m: l=2
      select case(m)
        case (-2)
          Y_Value=(1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)*x*y)*exp(-a*r2)
        case (-1)
          Y_Value=(1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)*y*z)*exp(-a*r2)
        case (0)
          Y_Value=(1.0d0/4.0d0*sqrt(5.0d0)/sqrt(PI)*(3.0d0*z**2-r2))*exp(-a*r2)
        case (1)
          Y_Value=(1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)*x*z)*exp(-a*r2)
        case (2)
          Y_Value=(-1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)*(y**2-x**2))*exp(-a*r2)
        case default
          print*,'Error: Y_Value not implemented for l=',l,'m=',m
          stop
      end select
    case (3)
      ! selection on m: l=3
      select case(m)
        case (-3)
          Y_Value=(-1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)*y*(y**2-3.0d0*x**2))*exp(-a*r2)
        case (-2)
          Y_Value=(1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)*x*y*z)*exp(-a*r2)
        case (-1)
          Y_Value=(1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)/sqrt(PI)*y*(5.0d0*z**2-r2))*exp(-a*r2)
        case (0)
          Y_Value=(1.0d0/4.0d0*sqrt(7.0d0)/sqrt(PI)*z*(5.0d0*z**2-3.0d0*r2))*exp(-a*r2)
        case (1)
          Y_Value=(1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)/sqrt(PI)*x*(5.0d0*z**2-r2))*exp(-a*r2)
        case (2)
          Y_Value=(-1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)*(y**2-x**2)*z)*exp(-a*r2)
        case (3)
          Y_Value=(-1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)*x*(3.0d0*y**2-x**2))*exp(-a*r2)
        case default
          print*,'Error: Y_Value not implemented for l=',l,'m=',m
          stop
      end select
    case (4)
      ! selection on m: l=4
      select case(m)
        case (-4)
          Y_Value=(-3.0d0/4.0d0*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)*x*y*(y**2-x**2))*exp(-a*r2)
        case (-3)
          Y_Value=(-3.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)*y*(y**2-3.0d0*x**2)*z)*exp(-a*r2)
        case (-2)
          Y_Value=(3.0d0/4.0d0*sqrt(5.0d0)/sqrt(PI)*x*y*(7.0d0*z**2-r2))*exp(-a*r2)
        case (-1)
          Y_Value=(3.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)/sqrt(PI)*y*z*(7.0d0*z**2-3.0d0*r2))*exp(-a*r2)
        case (0)
          Y_Value=(3.0d0/16.0d0/sqrt(PI)*(35.0d0*z**4+3.0d0*r2*(r2-10.0d0*z**2)))*exp(-a*r2)
        case (1)
          Y_Value=(3.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)/sqrt(PI)*x*z*(7.0d0*z**2-3.0d0*r2))*exp(-a*r2)
        case (2)
          Y_Value=(-3.0d0/8.0d0*sqrt(5.0d0)/sqrt(PI)*(7.0d0*(y**2-x**2)*z**2+r2*(x**2-y**2)))*exp(-a*r2)
        case (3)
          Y_Value=(-3.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)*x*(3.0d0*y**2-x**2)*z)*exp(-a*r2)
        case (4)
          Y_Value=(3.0d0/16.0d0*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)*(y**4+x**2*(x**2-6.0d0*y**2)))*exp(-a*r2)
        case default
          print*,'Error: Y_Value not implemented for l=',l,'m=',m
          stop
      end select
    case (5)
      ! selection on m: l=5
      select case(m)
        case (-5)
          Y_Value=(3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)*y*(y**4+5.0d0*x**2*(x**2-2.0d0*y**2)))*exp(-a*r2 &
)
        case (-4)
          Y_Value=(-3.0d0/4.0d0*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)*x*y*(y**2-x**2)*z)*exp(-a*r2)
        case (-3)
          Y_Value=(-1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)*y*(9.0d0*(y**2-3.0d0*x**2)*z**2+r2* &
(3.0d0*x**2-y**2)))*exp(-a*r2)
        case (-2)
          Y_Value=(1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)*x*y*z*(3.0d0*z**2-r2))*exp(-a*r2)
        case (-1)
          Y_Value=(1.0d0/16.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)*y*(7.0d0*3.0d0*z**4+r2*(r2-14.0d0*z**2)))*exp(-a*r2)
        case (0)
          Y_Value=(1.0d0/16.0d0*sqrt(11.0d0)/sqrt(PI)*z*(63.0d0*z**4+5.0d0*r2*(3.0d0*r2-14.0d0*z**2)))*exp(-a*r2)
        case (1)
          Y_Value=(1.0d0/16.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)*x*(7.0d0*3.0d0*z**4+r2*(r2-14.0d0*z**2)))*exp(-a*r2)
        case (2)
          Y_Value=(-1.0d0/8.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)*z*(3.0d0*(y**2-x**2)*z**2+r2*(x**2-y**2) &
))*exp(-a*r2)
        case (3)
          Y_Value=(-1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)*x*(9.0d0*(3.0d0*y**2-x**2)*z**2+r2* &
(x**2-3.0d0*y**2)))*exp(-a*r2)
        case (4)
          Y_Value=(3.0d0/16.0d0*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)*(y**4+x**2*(x**2-6.0d0*y**2))*z)*exp(-a*r2)
        case (5)
          Y_Value=(3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)*x*(5.0d0*y**4+x**2*(x**2-10.0d0*y**2)))*exp( &
-a*r2)
        case default
          print*,'Error: Y_Value not implemented for l=',l,'m=',m
          stop
      end select
    case (6)
      ! selection on m: l=6
      select case(m)
        case (-6)
          Y_Value=(1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)*x*y*(3.0d0*y**4+x**2* &
(3.0d0*x**2-10.0d0*y**2)))*exp(-a*r2)
        case (-5)
          Y_Value=(3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)*y*(y**4+5.0d0*x**2*(x**2-2.0d0*y**2 &
))*z)*exp(-a*r2)
        case (-4)
          Y_Value=(-3.0d0/8.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)*x*y*(11.0d0*(y**2-x**2)*z**2+r2*(x**2-y**2)))*exp(-a*r2)
        case (-3)
          Y_Value=(-1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)*y*z*(11.0d0*(y**2 &
-3.0d0*x**2)*z**2+3.0d0*r2*(3.0d0*x**2-y**2)))*exp(-a*r2)
        case (-2)
          Y_Value=(1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)*x*y*(11.0d0*3.0d0*z**4 &
+r2*(r2-2.0d0*3.0d0**2*z**2)))*exp(-a*r2)
        case (-1)
          Y_Value=(1.0d0/16.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)*y*z*(11.0d0*3.0d0*z**4+5.0d0*r2*(r2-2.0d0*3.0d0*z**2 &
)))*exp(-a*r2)
        case (0)
          Y_Value=(1.0d0/32.0d0*sqrt(13.0d0)/sqrt(PI)*(231.0d0*z**6+5.0d0*r2*(r2*(21.0d0*z**2-r2)-63.0d0*z**4)))*exp(-a*r2)
        case (1)
          Y_Value=(1.0d0/16.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)*x*z*(11.0d0*3.0d0*z**4+5.0d0*r2*(r2-2.0d0*3.0d0*z**2 &
)))*exp(-a*r2)
        case (2)
          Y_Value=(-1.0d0/(2.0d0**5*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)*(11.0d0*3.0d0*(y**2 &
-x**2)*z**4+r2*(2.0d0*3.0d0**2*(x**2-y**2)*z**2+r2*(y**2-x**2))))*exp(-a*r2)
        case (3)
          Y_Value=(-1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)*x*z*(11.0d0*(3.0d0*y**2 &
-x**2)*z**2+3.0d0*r2*(x**2-3.0d0*y**2)))*exp(-a*r2)
        case (4)
          Y_Value=(3.0d0/32.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)*(11.0d0*(y**4+x**2*(x**2-6.0d0*y**2))*z**2+r2*(x**2*(6.0d0*y**2 &
-x**2)-y**4)))*exp(-a*r2)
        case (5)
          Y_Value=(3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)*x*(5.0d0*y**4+x**2*(x**2 &
-10.0d0*y**2))*z)*exp(-a*r2)
        case (6)
          Y_Value=(-1.0d0/(2.0d0**5*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)*(y**6+x**2*(x**2* &
(5.0d0*3.0d0*y**2-x**2)-5.0d0*3.0d0*y**4)))*exp(-a*r2)
        case default
          print*,'Error: Y_Value not implemented for l=',l,'m=',m
          stop
      end select
    case default
      print*,'Error: Y_Value not implemented for l=',l,'m=',m
      stop
  end select
 
end function

