module mod_D_X_R_shell_opt_0
  
  use mod_A_functions
  
  implicit none
  
  private
  public :: D_X_R_shell_opt_0
  
contains
  
  

!> Compute D X Y/R coulomb integral for the right hand shell 
!!
recursive function D_X_R_shell_opt_0(a12,a3,r3,c,lmax,value)
  
  use mod_A_functions
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12      !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3       !< exponent of the third solid harmonics 
  real(kind=8), intent(in)  :: r3(3)    !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)     !< coeffcient of the hermit decomposition of the product of orbital 1 and 2
  integer     , intent(in)  :: lmax     !< l of maximum non 0 coeff in hermit decomposition~%
  real(kind=8), intent(out) :: value(*) !< result value
  
  ! return value
  
  logical :: D_X_R_shell_opt_0
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u50
  real(kind=8) :: u51
  real(kind=8) :: u55
  real(kind=8) :: u61
  real(kind=8) :: u62
  real(kind=8) :: u63
  real(kind=8) :: u72
  real(kind=8) :: u74
  real(kind=8) :: u75
  real(kind=8) :: u76
  real(kind=8) :: u77
  real(kind=8) :: u78
  real(kind=8) :: u8
  real(kind=8) :: u81
  real(kind=8) :: u82
  real(kind=8) :: u83
  real(kind=8) :: u84
  real(kind=8) :: u85
  real(kind=8) :: u86
  real(kind=8) :: u87
  real(kind=8) :: u88
  real(kind=8) :: u89
  real(kind=8) :: u9
  real(kind=8) :: u90
  real(kind=8) :: u91
  real(kind=8) :: u92
  real(kind=8) :: u93
  real(kind=8) :: u94
  real(kind=8) :: u95
  real(kind=8) :: u96
  real(kind=8) :: u98
  real(kind=8) :: u99

  
  ! init computation flag
  D_X_R_shell_opt_0=.true.
  
  ! compute local quantities
  p=sqrt((a12*a3)/(a12+a3))
  d=sqrt(sum(r3**2))
  
  ! renormalize xyz
  x=0.0d0
  y=0.0d0
  z=0.0d0
  if (d>0.0d0) then
    x=r3(1)/d
    y=r3(2)/d
    z=r3(3)/d
  end if
  
  ! compute local quantities
  pd=p*d
  pd2=pd*pd
  erfpd =erf(pd)
  exppd2=exp(-pd2)
  
  ! computation section
  value(1:1)=0.0d0
  u61=A1_(p,pd,erfpd,exppd2)
  value(1)=value(1)+(c(1)*(u61))
  if ( lmax .eq. 0 ) return
  u61=1.0d0*A2_(p,pd,erfpd,exppd2)*p**2
  value(1)=value(1)+(c(2)*(u61*x))
  value(1)=value(1)+(c(3)*(u61*y))
  value(1)=value(1)+(c(4)*(u61*z))
  if ( lmax .eq. 1 ) return
  u61=A3_(p,pd,erfpd,exppd2)
  u55=p**3
  u76=-1.88063194515918762d-1*(4.0d0*exppd2+5.31736155271654808d0*u61)*&
 &u55
  u8=3.0d0*u61*u55
  u55=x**2
  value(1)=value(1)+(c(5)*(u8*u55+u76))
  u61=x*y
  value(1)=value(1)+(c(6)*(u8*u61))
  u63=x*z
  value(1)=value(1)+(c(7)*(u8*u63))
  u50=y**2
  value(1)=value(1)+(c(8)*(u8*u50+u76))
  u62=y*z
  value(1)=value(1)+(c(9)*(u8*u62))
  u51=z**2
  value(1)=value(1)+(c(10)*(u8*u51+u76))
  if ( lmax .eq. 2 ) return
  u76=A4_(p,pd,erfpd,exppd2)
  u8=p**4
  u82=u76*u8
  u85=-9.0d0*u82
  u72=-8.0d0*exppd2
  u88=5.64189583547756287d-1*(2.65868077635827404d1*u76+u72*pd)*u8
  u76=u88*u55
  value(1)=value(1)+(c(11)*(x*(u76+u85)))
  u8=-3.0d0*u82
  u82=(u76+u8)
  value(1)=value(1)+(c(12)*(u82*y))
  value(1)=value(1)+(c(13)*(u82*z))
  u76=u88*u50
  u83=(u76+u8)
  value(1)=value(1)+(c(14)*(x*u83))
  u82=u61*z
  value(1)=value(1)+(c(15)*(u88*u82))
  u90=u88*u51
  u84=(u90+u8)
  value(1)=value(1)+(c(16)*(x*u84))
  value(1)=value(1)+(c(17)*(y*(u76+u85)))
  value(1)=value(1)+(c(18)*(u83*z))
  value(1)=value(1)+(c(19)*(y*u84))
  value(1)=value(1)+(c(20)*(z*(u90+u85)))
  if ( lmax .eq. 3 ) return
  u76=A5_(p,pd,erfpd,exppd2)
  u90=p**5
  u88=(2.65868077635827404d1*u76-u72)*u90
  u8=3.38513750128653772d-1*u88
  u85=u76*u90
  u74=-9.0d1*u85
  u83=-1.6d1*exppd2
  u84=+u83*pd2
  u96=5.64189583547756287d-1*(1.86107654345079183d2*u76+u84)*u90
  u76=u96*u55
  value(1)=value(1)+(c(21)*(u55*(u74+u76)+u8))
  u90=-4.5d1*u85
  u75=x*(u76+u90)
  value(1)=value(1)+(c(22)*(u75*y))
  value(1)=value(1)+(c(23)*(u75*z))
  u9=1.12837916709551257d-1*u88
  u88=-1.5d1*u85
  u85=(u88+u76)
  u93=u88*u55+u9
  value(1)=value(1)+(c(24)*(u85*u50+u93))
  value(1)=value(1)+(c(25)*((u76+u88)*u62))
  value(1)=value(1)+(c(26)*(u85*u51+u93))
  u76=u96*u50
  u86=(u76+u90)
  value(1)=value(1)+(c(27)*(u61*u86))
  value(1)=value(1)+(c(28)*(x*(u76+u88)*z))
  u85=u96*u51
  value(1)=value(1)+(c(29)*(u61*(u85+u88)))
  u87=(u85+u90)
  value(1)=value(1)+(c(30)*(u63*u87))
  value(1)=value(1)+(c(31)*(u50*(u74+u76)+u8))
  value(1)=value(1)+(c(32)*(y*u86*z))
  value(1)=value(1)+(c(33)*((u88+u76)*u51+u88*u50+u9))
  value(1)=value(1)+(c(34)*(u62*u87))
  value(1)=value(1)+(c(35)*(u51*(u74+u85)+u8))
  if ( lmax .eq. 4 ) return
  u76=A6_(p,pd,erfpd,exppd2)
  u85=p**6
  u96=u76*u85
  u8=2.25d2*u96
  u9=+u83*pd
  u88=(1.86107654345079183d2*u76+u9)*u85
  u90=-5.64189583547756287d0*u88
  u74=2.0d0*pd2
  u81=5.64189583547756287d-1*(1.67496888910571265d3*u76+u9*(9.0d0+u74))&
 &*u85
  u76=u81*u55
  value(1)=value(1)+(c(36)*(x*(u55*(u90+u76)+u8)))
  u85=4.5d1*u96
  u9=-3.38513750128653772d0*u88
  u87=(u55*(u9+u76)+u85)
  value(1)=value(1)+(c(37)*(u87*y))
  value(1)=value(1)+(c(38)*(u87*z))
  u87=-1.69256875064326886d0*u88
  u86=-5.64189583547756287d-1*u88
  u88=(u87+u76)
  u77=u86*u55
  u94=u77+u85
  value(1)=value(1)+(c(39)*(x*(u88*u50+u94)))
  value(1)=value(1)+(c(40)*(x*(u76+u87)*u62))
  value(1)=value(1)+(c(41)*(x*(u88*u51+u94)))
  u89=(u86+u76)
  u76=u89*u50
  u95=u87*u55+u85
  value(1)=value(1)+(c(42)*(y*(u76+u95)))
  u88=1.5d1*u96
  u96=u77+u88
  value(1)=value(1)+(c(43)*((u76+u96)*z))
  u77=u89*u51
  value(1)=value(1)+(c(44)*(y*(u77+u96)))
  value(1)=value(1)+(c(45)*(z*(u77+u95)))
  u96=u81*u50
  u77=(u50*(u9+u96)+u85)
  value(1)=value(1)+(c(46)*(x*u77))
  value(1)=value(1)+(c(47)*(u61*(u96+u87)*z))
  u78=(u86+u96)*u51
  u89=u86*u50
  value(1)=value(1)+(c(48)*(x*(u78+u89+u88)))
  u76=u81*u51
  value(1)=value(1)+(c(49)*(u82*(u76+u87)))
  u82=(u51*(u9+u76)+u85)
  value(1)=value(1)+(c(50)*(x*u82))
  value(1)=value(1)+(c(51)*(y*(u50*(u90+u96)+u8)))
  value(1)=value(1)+(c(52)*(u77*z))
  value(1)=value(1)+(c(53)*(y*((u87+u96)*u51+u89+u85)))
  value(1)=value(1)+(c(54)*(z*(u78+u87*u50+u85)))
  value(1)=value(1)+(c(55)*(y*u82))
  value(1)=value(1)+(c(56)*(z*(u51*(u90+u76)+u8)))
  if ( lmax .eq. 5 ) return
  u96=A7_(p,pd,erfpd,exppd2)
  u76=p**7
  u89=(-u83+1.86107654345079183d2*u96)*u76
  u82=-1.20897767903090633d0*u89
  u90=u96*u76
  u77=4.725d3*u90
  u88=-3.2d1*exppd2*pd2
  u81=(1.67496888910571265d3*u96+u88)*u76
  u85=-8.46284375321634431d0*u81
  u9=(1.84246577801628391d4*u96+u88*(1.1d1+u74))*u76
  u88=5.64189583547756287d-1*u9
  u83=u88*u55
  value(1)=value(1)+(c(57)*(u55*(u77+u55*(u83+u85))+u82))
  u78=1.575d3*u90
  u94=-5.64189583547756287d0*u81
  u87=x*(u55*(u94+u83)+u78)
  value(1)=value(1)+(c(58)*(u87*y))
  value(1)=value(1)+(c(59)*(u87*z))
  u87=-2.41795535806181266d-1*u89
  u89=3.15d2*u90
  u86=6.3d2*u90
  u8=-3.38513750128653772d0*u81
  u95=-5.64189583547756287d-1*u81
  u93=(u89+u55*(u83+u8))
  u99=u95*u55
  u74=u55*(u86+u99)+u87
  value(1)=value(1)+(c(60)*(u93*u50+u74))
  value(1)=value(1)+(c(61)*((u55*(u8+u83)+u89)*u62))
  value(1)=value(1)+(c(62)*(u93*u51+u74))
  u93=9.45d2*u90
  u74=-1.69256875064326886d0*u81
  u90=(u74+u83)
  u91=u90*u50
  u75=u74*u55
  u92=u75+u93
  value(1)=value(1)+(c(63)*(u61*(u91+u92)))
  u98=u99+u89
  value(1)=value(1)+(c(64)*(x*(u91+u98)*z))
  u91=u90*u51
  value(1)=value(1)+(c(65)*(u61*(u91+u98)))
  value(1)=value(1)+(c(66)*(u63*(u91+u92)))
  u91=(u83+u95)
  u99=u86+u8*u55
  u90=u89*u55+u87
  value(1)=value(1)+(c(67)*(u50*(u99+u91*u50)+u90))
  u92=(u95+u83)
  u83=u75+u89
  value(1)=value(1)+(c(68)*(y*(u92*u50+u83)*z))
  u98=1.61197023870787511d-1*(5.58322963035237549d2*u96+u72)*u76
  u75=5.64189583547756287d-1*u81
  u72=2.25675833419102515d0*(5.02490666731713794d3*u96+u84*(6.0d0+pd2))&
 &*u76
  u81=-5.64189583547756287d-1*u9
  u9=u81*u55
  value(1)=value(1)+(c(69)*(u50*(u95+u55*(u9+u72)+(u9+u75)*u50)+u55*(u9&
 &5+u75*u55)+u98))
  value(1)=value(1)+(c(70)*(u62*(u92*u51+u83)))
  value(1)=value(1)+(c(71)*(u51*(u99+u91*u51)+u90))
  u55=u88*u50
  u72=(u50*(u94+u55)+u78)
  value(1)=value(1)+(c(72)*(u61*u72))
  value(1)=value(1)+(c(73)*(x*(u50*(u8+u55)+u89)*z))
  u81=(u74+u55)*u51
  u99=u95*u50
  value(1)=value(1)+(c(74)*(u61*(u81+u99+u89)))
  u90=u74*u50
  value(1)=value(1)+(c(75)*(u63*((u95+u55)*u51+u90+u89)))
  u91=u88*u51
  value(1)=value(1)+(c(76)*(u61*(u51*(u8+u91)+u89)))
  u61=(u51*(u94+u91)+u78)
  value(1)=value(1)+(c(77)*(u63*u61))
  value(1)=value(1)+(c(78)*(u50*(u77+u50*(u55+u85))+u82))
  value(1)=value(1)+(c(79)*(y*u72*z))
  value(1)=value(1)+(c(80)*((u89+u50*(u55+u8))*u51+u50*(u86+u99)+u87))
  value(1)=value(1)+(c(81)*(u62*(u81+u90+u93)))
  value(1)=value(1)+(c(82)*(u51*(u86+u8*u50+(u55+u95)*u51)+u89*u50+u87)&
 &)
  value(1)=value(1)+(c(83)*(u62*u61))
  value(1)=value(1)+(c(84)*(u51*(u77+u51*(u91+u85))+u82))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_R_shell_opt_0=.false.
end function

  
end module
