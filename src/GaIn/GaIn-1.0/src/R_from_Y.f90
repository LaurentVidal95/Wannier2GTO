
module mod_R_from_Y

  implicit none

contains

!> Expand a Solid real spherical harmonic in terms of x^i*y^j*z^k
!! 
recursive subroutine R_from_Y(l,m,c)
 
  implicit none
 
  integer     , intent(in):: l!< l for Ylm, l must be < 6
  integer     , intent(in):: m!< m for Ylm
  real(kind=8), intent(out),dimension(455)  :: c!< coefficient of decomposition in the basis {1,x,y,z,x^2,...}
 
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
 
  ! init c coefffs
  c=0.0d0
  ! selection on l
  select case(l)
    case (0)
      ! selection on m: l=0
      select case(m)
        case (0)
          c(1)=1.0d0/2.0d0/sqrt(PI)
        case default
          print*,'Error: R_from_Y not implemented for l=',l,'m=',m
          stop
      end select
    case (1)
      ! selection on m: l=1
      select case(m)
        case (-1)
          c(3)=1.0d0/2.0d0*sqrt(3.0d0)/sqrt(PI)
        case (0)
          c(4)=1.0d0/2.0d0*sqrt(3.0d0)/sqrt(PI)
        case (1)
          c(2)=1.0d0/2.0d0*sqrt(3.0d0)/sqrt(PI)
        case default
          print*,'Error: R_from_Y not implemented for l=',l,'m=',m
          stop
      end select
    case (2)
      ! selection on m: l=2
      select case(m)
        case (-2)
          c(6)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)
        case (-1)
          c(9)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)
        case (0)
          c(5)=-1.0d0/4.0d0*sqrt(5.0d0)/sqrt(PI)
          c(8)=-1.0d0/4.0d0*sqrt(5.0d0)/sqrt(PI)
          c(10)=1.0d0/2.0d0*sqrt(5.0d0)/sqrt(PI)
        case (1)
          c(7)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)
        case (2)
          c(5)=1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)
          c(8)=-1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)/sqrt(PI)
        case default
          print*,'Error: R_from_Y not implemented for l=',l,'m=',m
          stop
      end select
    case (3)
      ! selection on m: l=3
      select case(m)
        case (-3)
          c(12)=3.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(17)=-1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
        case (-2)
          c(15)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
        case (-1)
          c(12)=-1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(17)=-1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(19)=1.0d0/sqrt(2.0d0)*sqrt(3.0d0)*sqrt(7.0d0)/sqrt(PI)
        case (0)
          c(13)=-3.0d0/4.0d0*sqrt(7.0d0)/sqrt(PI)
          c(18)=-3.0d0/4.0d0*sqrt(7.0d0)/sqrt(PI)
          c(20)=1.0d0/2.0d0*sqrt(7.0d0)/sqrt(PI)
        case (1)
          c(11)=-1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(14)=-1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(16)=1.0d0/sqrt(2.0d0)*sqrt(3.0d0)*sqrt(7.0d0)/sqrt(PI)
        case (2)
          c(13)=1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(18)=-1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
        case (3)
          c(11)=1.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(14)=-3.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
        case default
          print*,'Error: R_from_Y not implemented for l=',l,'m=',m
          stop
      end select
    case (4)
      ! selection on m: l=4
      select case(m)
        case (-4)
          c(22)=3.0d0/4.0d0*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(27)=-3.0d0/4.0d0*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
        case (-3)
          c(25)=9.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(32)=-3.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
        case (-2)
          c(22)=-3.0d0/4.0d0*sqrt(5.0d0)/sqrt(PI)
          c(27)=-3.0d0/4.0d0*sqrt(5.0d0)/sqrt(PI)
          c(29)=9.0d0/2.0d0*sqrt(5.0d0)/sqrt(PI)
        case (-1)
          c(25)=-9.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)/sqrt(PI)
          c(32)=-9.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)/sqrt(PI)
          c(34)=3.0d0/sqrt(2.0d0)*sqrt(5.0d0)/sqrt(PI)
        case (0)
          c(21)=9.0d0/16.0d0/sqrt(PI)
          c(24)=9.0d0/8.0d0/sqrt(PI)
          c(26)=-9.0d0/2.0d0/sqrt(PI)
          c(31)=9.0d0/16.0d0/sqrt(PI)
          c(33)=-9.0d0/2.0d0/sqrt(PI)
          c(35)=3.0d0/2.0d0/sqrt(PI)
        case (1)
          c(23)=-9.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)/sqrt(PI)
          c(28)=-9.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)/sqrt(PI)
          c(30)=3.0d0/sqrt(2.0d0)*sqrt(5.0d0)/sqrt(PI)
        case (2)
          c(21)=-3.0d0/8.0d0*sqrt(5.0d0)/sqrt(PI)
          c(26)=9.0d0/4.0d0*sqrt(5.0d0)/sqrt(PI)
          c(31)=3.0d0/8.0d0*sqrt(5.0d0)/sqrt(PI)
          c(33)=-9.0d0/4.0d0*sqrt(5.0d0)/sqrt(PI)
        case (3)
          c(23)=3.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(28)=-9.0d0/(2.0d0**2*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
        case (4)
          c(21)=3.0d0/16.0d0*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(24)=-9.0d0/8.0d0*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
          c(31)=3.0d0/16.0d0*sqrt(5.0d0)*sqrt(7.0d0)/sqrt(PI)
        case default
          print*,'Error: R_from_Y not implemented for l=',l,'m=',m
          stop
      end select
    case (5)
      ! selection on m: l=5
      select case(m)
        case (-5)
          c(37)=15.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(42)=-15.0d0/(2.0d0**3*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(51)=3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (-4)
          c(40)=3.0d0/4.0d0*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(47)=-3.0d0/4.0d0*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (-3)
          c(37)=-3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(42)=-1.0d0/(2.0d0**3*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(44)=3.0d0/(2.0d0*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(51)=1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(53)=-1.0d0/(2.0d0*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (-2)
          c(40)=-1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(47)=-1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(49)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (-1)
          c(37)=1.0d0/16.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(42)=1.0d0/8.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(44)=-1.0d0/4.0d0*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(51)=1.0d0/16.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(53)=-1.0d0/4.0d0*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(55)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (0)
          c(38)=15.0d0/16.0d0*sqrt(11.0d0)/sqrt(PI)
          c(43)=15.0d0/8.0d0*sqrt(11.0d0)/sqrt(PI)
          c(45)=-5.0d0/2.0d0*sqrt(11.0d0)/sqrt(PI)
          c(52)=15.0d0/16.0d0*sqrt(11.0d0)/sqrt(PI)
          c(54)=-5.0d0/2.0d0*sqrt(11.0d0)/sqrt(PI)
          c(56)=1.0d0/2.0d0*sqrt(11.0d0)/sqrt(PI)
        case (1)
          c(36)=1.0d0/16.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(39)=1.0d0/8.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(41)=-1.0d0/4.0d0*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(46)=1.0d0/16.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(48)=-1.0d0/4.0d0*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(50)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (2)
          c(38)=-1.0d0/8.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(45)=1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(52)=1.0d0/8.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(54)=-1.0d0/4.0d0*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (3)
          c(36)=-1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(39)=1.0d0/(2.0d0**3*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(41)=1.0d0/(2.0d0*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(46)=3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(48)=-3.0d0/(2.0d0*sqrt(2.0d0))*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (4)
          c(38)=3.0d0/16.0d0*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(43)=-9.0d0/8.0d0*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(52)=3.0d0/16.0d0*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
        case (5)
          c(36)=3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(39)=-15.0d0/(2.0d0**3*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
          c(46)=15.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)/sqrt(PI)
        case default
          print*,'Error: R_from_Y not implemented for l=',l,'m=',m
          stop
      end select
    case (6)
      ! selection on m: l=6
      select case(m)
        case (-6)
          c(58)=1.0d0/(2.0d0**4*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(63)=-5.0d0/(2.0d0**3*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(72)=1.0d0/(2.0d0**4*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (-5)
          c(61)=15.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(68)=-15.0d0/(2.0d0**3*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(79)=3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (-4)
          c(58)=-3.0d0/8.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(65)=15.0d0/4.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(72)=3.0d0/8.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(74)=-15.0d0/4.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (-3)
          c(61)=-1.0d0/(2.0d0**4*sqrt(2.0d0))*3.0d0**2*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(68)=-1.0d0/(2.0d0**3*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(70)=1.0d0/(2.0d0*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(79)=1.0d0/(2.0d0**4*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(81)=-1.0d0/(2.0d0*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (-2)
          c(58)=1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(63)=1.0d0/(2.0d0**3*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(65)=-1.0d0/sqrt(2.0d0)*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(72)=1.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(74)=-1.0d0/sqrt(2.0d0)*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(76)=1.0d0/sqrt(2.0d0)*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (-1)
          c(61)=5.0d0/16.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(68)=5.0d0/8.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(70)=-5.0d0/4.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(79)=5.0d0/16.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(81)=-5.0d0/4.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(83)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (0)
          c(57)=-5.0d0/32.0d0*sqrt(13.0d0)/sqrt(PI)
          c(60)=-15.0d0/32.0d0*sqrt(13.0d0)/sqrt(PI)
          c(62)=45.0d0/16.0d0*sqrt(13.0d0)/sqrt(PI)
          c(67)=-15.0d0/32.0d0*sqrt(13.0d0)/sqrt(PI)
          c(69)=45.0d0/8.0d0*sqrt(13.0d0)/sqrt(PI)
          c(71)=-15.0d0/4.0d0*sqrt(13.0d0)/sqrt(PI)
          c(78)=-5.0d0/32.0d0*sqrt(13.0d0)/sqrt(PI)
          c(80)=45.0d0/16.0d0*sqrt(13.0d0)/sqrt(PI)
          c(82)=-15.0d0/4.0d0*sqrt(13.0d0)/sqrt(PI)
          c(84)=1.0d0/2.0d0*sqrt(13.0d0)/sqrt(PI)
        case (1)
          c(59)=5.0d0/16.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(64)=5.0d0/8.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(66)=-5.0d0/4.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(73)=5.0d0/16.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(75)=-5.0d0/4.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(77)=1.0d0/2.0d0*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (2)
          c(57)=1.0d0/(2.0d0**5*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(60)=1.0d0/(2.0d0**5*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(62)=-1.0d0/(2.0d0*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(67)=-1.0d0/(2.0d0**5*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(71)=1.0d0/(2.0d0*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(78)=-1.0d0/(2.0d0**5*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(80)=1.0d0/(2.0d0*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(82)=-1.0d0/(2.0d0*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (3)
          c(59)=-1.0d0/(2.0d0**4*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(64)=1.0d0/(2.0d0**3*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(66)=1.0d0/(2.0d0*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(73)=1.0d0/(2.0d0**4*sqrt(2.0d0))*3.0d0**2*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(75)=-1.0d0/(2.0d0*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(5.0d0)*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (4)
          c(57)=-3.0d0/32.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(60)=15.0d0/32.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(62)=15.0d0/16.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(67)=15.0d0/32.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(69)=-45.0d0/8.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(78)=-3.0d0/32.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(80)=15.0d0/16.0d0*sqrt(7.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (5)
          c(59)=3.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(64)=-15.0d0/(2.0d0**3*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(73)=15.0d0/(2.0d0**4*sqrt(2.0d0))*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
        case (6)
          c(57)=1.0d0/(2.0d0**5*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(60)=-5.0d0/(2.0d0**5*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(67)=5.0d0/(2.0d0**5*sqrt(2.0d0))*3.0d0**1*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
          c(78)=-1.0d0/(2.0d0**5*sqrt(2.0d0))*sqrt(3.0d0)*sqrt(7.0d0)*sqrt(11.0d0)*sqrt(13.0d0)/sqrt(PI)
        case default
          print*,'Error: R_from_Y not implemented for l=',l,'m=',m
          stop
      end select
    case default
      print*,'Error: R_from_Y not implemented for l=',l,'m=',m
      stop
  end select
 
end subroutine

end module
