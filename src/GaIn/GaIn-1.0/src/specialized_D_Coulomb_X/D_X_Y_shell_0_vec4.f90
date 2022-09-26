module mod_D_X_Y_shell_optv4_0
  
  implicit none
  
  private
  public :: D_X_Y_shell_optv4_0
  integer, parameter :: NVEC=4
  
contains
  
  

!> Compute D X Y/R coulomb integral for the right hand shell 
!!
recursive function D_X_Y_shell_optv4_0(a12,a3,r3,c,lmax,value)
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12          !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3(NVEC)     !< exponent of the third solid harmonics 
  real(kind=8), intent(in)  :: r3(NVEC,3)   !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)         !< coeffcient of the hermit decomposition of the product of orbital 1 and 2
  integer     , intent(in)  :: lmax         !< l of maximum non 0 coeff in hermit decomposition~%
  real(kind=8), intent(out) :: value(*)     !< result value
  
  ! return value
  logical :: D_X_Y_shell_optv4_0
  
  ! local variables
  integer      :: i
  real(kind=8) :: d(NVEC)
  real(kind=8) :: p(NVEC)
  real(kind=8) :: pd(NVEC)
  real(kind=8) :: pd2(NVEC)
  real(kind=8) :: erfpd(NVEC)
  real(kind=8) :: exppd2(NVEC)
  real(kind=8) :: x(NVEC)
  real(kind=8) :: y(NVEC)
  real(kind=8) :: z(NVEC)
  
  ! C function interface
  INTERFACE
    PURE RECURSIVE SUBROUTINE my_erf_v4(pd, exppd2, erfpd) BIND(C,NAME="my_erf_v4_")
     real(kind=8), intent(in)  :: pd(4)
     real(kind=8), intent(out) :: erfpd(4)
     real(kind=8), intent(out) :: exppd2(4)
    END SUBROUTINE my_erf_v4
  END INTERFACE
  
  ! temporary variables
  real(kind=8) :: u56(NVEC)
  real(kind=8) :: u57(NVEC)
  real(kind=8) :: u60(NVEC)
  real(kind=8) :: u62(NVEC)
  real(kind=8) :: u63(NVEC)
  real(kind=8) :: u64(NVEC)
  real(kind=8) :: u66(NVEC)
  real(kind=8) :: u74(NVEC)
  real(kind=8) :: u75(NVEC)
  real(kind=8) :: u76(NVEC)
  real(kind=8) :: u77(NVEC)
  real(kind=8) :: u78(NVEC)
  real(kind=8) :: u79(NVEC)
  real(kind=8) :: u8(NVEC)
  real(kind=8) :: u80(NVEC)
  real(kind=8) :: u81(NVEC)
  real(kind=8) :: u83(NVEC)
  real(kind=8) :: u84(NVEC)
  real(kind=8) :: u85(NVEC)
  real(kind=8) :: u86(NVEC)
  real(kind=8) :: u87(NVEC)
  real(kind=8) :: u88(NVEC)
  real(kind=8) :: u89(NVEC)
  real(kind=8) :: u9(NVEC)
  real(kind=8) :: u90(NVEC)
  real(kind=8) :: u91(NVEC)
  real(kind=8) :: u92(NVEC)
  real(kind=8) :: u93(NVEC)
  real(kind=8) :: u94(NVEC)
  real(kind=8) :: u95(NVEC)
  real(kind=8) :: u96(NVEC)
  real(kind=8) :: u97(NVEC)
  real(kind=8) :: u98(NVEC)
  real(kind=8) :: u99(NVEC)
  real(kind=8) :: value_1_(NVEC)

  
  ! init values
  value_1_ = 0.0d0
  
  ! init computation flag
  D_X_Y_shell_optv4_0=.true.
  
  ! compute local quantities
  d=sqrt(sum(r3**2,dim=2))+1.0d-32
  p=sqrt((a12*a3)/(a12+a3))
  pd=p*d
  call my_erf_v4(pd, exppd2, erfpd)
  pd2=pd*pd
  
  ! renormalize xyz
  x=0.0d0
  y=0.0d0
  z=0.0d0
  x=r3(:,1)/d
  y=r3(:,2)/d
  z=r3(:,3)/d
  
  ! computation section
  u62=2.82094791773878144d-1*p*A1_v(pd,exppd2,erfpd)
  value_1_=value_1_+(c(1)*(u62))
  if ( lmax .eq. 0 ) go to 100
  u62=2.82094791773878144d-1*p**2*A2_v(pd,exppd2,erfpd)
  value_1_=value_1_+(c(2)*(u62*x))
  value_1_=value_1_+(c(3)*(u62*y))
  value_1_=value_1_+(c(4)*(u62*z))
  if ( lmax .eq. 1 ) go to 100
  u62=A3_v(pd,exppd2,erfpd)
  u60=p**3
  u99=-2.99311870208610936d-2*u60*(7.08981540362206411d0*exppd2+9.42477&
 &796076937972d0*u62)
  u75=8.46284375321634431d-1*u60*u62
  u60=x**2
  value_1_=value_1_+(c(5)*(u75*u60+u99))
  u62=x*y
  value_1_=value_1_+(c(6)*(u75*u62))
  u64=x*z
  value_1_=value_1_+(c(7)*(u75*u64))
  u56=y**2
  value_1_=value_1_+(c(8)*(u75*u56+u99))
  u63=y*z
  value_1_=value_1_+(c(9)*(u75*u63))
  u57=z**2
  value_1_=value_1_+(c(10)*(u75*u57+u99))
  if ( lmax .eq. 2 ) go to 100
  u99=A5_v(pd,exppd2,erfpd)
  u75=p**4*pd
  u66=1.41796308072441282d1*exppd2
  u85=u75*(u66+4.71238898038468986d1*u99)
  u93=-5.38761366375499685d-2*u85
  u9=4.23142187660817215d0*u75*u99
  u99=u9*u60
  value_1_=value_1_+(c(11)*(x*(u99+u93)))
  u75=-1.79587122125166562d-2*u85
  u85=(u99+u75)
  value_1_=value_1_+(c(12)*(u85*y))
  value_1_=value_1_+(c(13)*(u85*z))
  u99=u9*u56
  u86=(u99+u75)
  value_1_=value_1_+(c(14)*(x*u86))
  u85=u62*z
  value_1_=value_1_+(c(15)*(u9*u85))
  u95=u9*u57
  u87=(u95+u75)
  value_1_=value_1_+(c(16)*(x*u87))
  value_1_=value_1_+(c(17)*(y*(u99+u93)))
  value_1_=value_1_+(c(18)*(u86*z))
  value_1_=value_1_+(c(19)*(y*u87))
  value_1_=value_1_+(c(20)*(z*(u95+u93)))
  if ( lmax .eq. 3 ) go to 100
  u75=A7_v(pd,exppd2,erfpd)
  u99=pd2*u75
  u95=p**5
  u86=2.0d0*pd2
  u9=u95*(3.2986722862692829d2*u99+u66*(7.0d0+u86))
  u93=7.69659094822142407d-3*u9
  u79=3.2986722862692829d2*u75
  u87=2.83592616144882564d1*exppd2
  u88=u95*pd2*(u87+u79)
  u91=-7.69659094822142407d-2*u88
  u8=2.96199531362572051d1*u95*u99
  u99=u8*u60
  value_1_=value_1_+(c(21)*(u60*(u91+u99)+u93))
  u95=-3.84829547411071204d-2*u88
  u78=x*(u99+u95)
  value_1_=value_1_+(c(22)*(u78*y))
  value_1_=value_1_+(c(23)*(u78*z))
  u74=2.56553031607380803d-3*u9
  u9=-1.28276515803690401d-2*u88
  u88=(u9+u99)
  u96=u9*u60+u74
  value_1_=value_1_+(c(24)*(u88*u56+u96))
  value_1_=value_1_+(c(25)*((u99+u9)*u63))
  value_1_=value_1_+(c(26)*(u88*u57+u96))
  u99=u8*u56
  u89=(u99+u95)
  value_1_=value_1_+(c(27)*(u62*u89))
  value_1_=value_1_+(c(28)*(x*(u99+u9)*z))
  u80=u8*u57
  value_1_=value_1_+(c(29)*(u62*(u80+u9)))
  u90=(u80+u95)
  value_1_=value_1_+(c(30)*(u64*u90))
  value_1_=value_1_+(c(31)*(u56*(u91+u99)+u93))
  value_1_=value_1_+(c(32)*(y*u89*z))
  value_1_=value_1_+(c(33)*((u9+u99)*u57+u9*u56+u74))
  value_1_=value_1_+(c(34)*(u63*u90))
  value_1_=value_1_+(c(35)*(u57*(u91+u80)+u93))
  if ( lmax .eq. 4 ) go to 100
  u95=p**6*pd
  u99=u95*(u79+u87)
  u90=1.92414773705535602d-1*u99
  u91=u95*u75
  u9=-2.96199531362572051d2*u91
  u74=-5.67185232289765129d1*exppd2
  u8=8.97935610625832809d-2*u95*(2.96880505764235461d3*u75+u74*pd2)
  u79=u8*u60
  value_1_=value_1_+(c(36)*(x*(u60*(u9+u79)+u90)))
  u75=3.84829547411071204d-2*u99
  u95=-1.7771971881754323d2*u91
  u93=(u60*(u95+u79)+u75)
  value_1_=value_1_+(c(37)*(u93*y))
  value_1_=value_1_+(c(38)*(u93*z))
  u93=-8.88598594087716152d1*u91
  u88=-2.96199531362572051d1*u91
  u91=(u93+u79)
  u80=u88*u60
  u97=u80+u75
  value_1_=value_1_+(c(39)*(x*(u91*u56+u97)))
  value_1_=value_1_+(c(40)*(x*(u79+u93)*u63))
  value_1_=value_1_+(c(41)*(x*(u91*u57+u97)))
  u92=(u88+u79)
  u79=u92*u56
  u98=u93*u60+u75
  value_1_=value_1_+(c(42)*(y*(u79+u98)))
  u97=1.28276515803690401d-2*u99
  u99=u80+u97
  value_1_=value_1_+(c(43)*((u79+u99)*z))
  u80=u92*u57
  value_1_=value_1_+(c(44)*(y*(u80+u99)))
  value_1_=value_1_+(c(45)*(z*(u80+u98)))
  u79=u8*u56
  u80=(u56*(u95+u79)+u75)
  value_1_=value_1_+(c(46)*(x*u80))
  value_1_=value_1_+(c(47)*(u62*(u79+u93)*z))
  u81=(u88+u79)*u57
  u92=u88*u56
  value_1_=value_1_+(c(48)*(x*(u81+u92+u97)))
  u98=u8*u57
  value_1_=value_1_+(c(49)*(u85*(u98+u93)))
  u97=(u57*(u95+u98)+u75)
  value_1_=value_1_+(c(50)*(x*u97))
  value_1_=value_1_+(c(51)*(y*(u56*(u9+u79)+u90)))
  value_1_=value_1_+(c(52)*(u80*z))
  value_1_=value_1_+(c(53)*(y*((u93+u79)*u57+u92+u75)))
  value_1_=value_1_+(c(54)*(z*(u81+u93*u56+u75)))
  value_1_=value_1_+(c(55)*(y*u97))
  value_1_=value_1_+(c(56)*(z*(u57*(u9+u98)+u90)))
  if ( lmax .eq. 5 ) go to 100
  u79=A9_v(pd,exppd2,erfpd)
  u98=pd2*u79
  u92=2.96880505764235461d3*u98
  u75=p**7
  u80=u75*(u87*(u86+9.0d0)+u92)
  u86=-2.13794193006150669d-2*u80
  u97=u75*pd2
  u93=u97*(2.96880505764235461d3*u79-u74)
  u8=4.48967805312916404d-1*u93
  u95=u75*u98
  u98=-3.99869367339472268d3*u95
  u96=u97*(3.26568556340659007d4*u79-1.13437046457953026d2*exppd2*pd2)
  u9=8.97935610625832809d-2*u96
  u90=u9*u60
  value_1_=value_1_+(c(57)*(u60*(u8+u60*(u90+u98))+u86))
  u88=1.49655935104305468d-1*u93
  u85=-2.66579578226314846d3*u95
  u91=x*(u60*(u85+u90)+u88)
  value_1_=value_1_+(c(58)*(u91*y))
  value_1_=value_1_+(c(59)*(u91*z))
  u89=-4.27588386012301337d-3*u80
  u91=2.99311870208610936d-2*u93
  u80=5.98623740417221873d-2*u93
  u81=-1.59947746935788907d3*u95
  u99=-2.66579578226314846d2*u95
  u78=+u99
  u74=(u91+u60*(u90+u81))
  u83=u78*u60
  u77=u60*(u80+u83)+u89
  value_1_=value_1_+(c(60)*(u74*u56+u77))
  value_1_=value_1_+(c(61)*((u60*(u81+u90)+u91)*u63))
  value_1_=value_1_+(c(62)*(u74*u57+u77))
  u74=8.97935610625832809d-2*u93
  u77=-7.99738734678944537d2*u95
  u93=(u77+u90)
  u95=u93*u56
  u76=u77*u60
  u94=u76+u74
  value_1_=value_1_+(c(63)*(u62*(u95+u94)))
  u84=u83+u91
  value_1_=value_1_+(c(64)*(x*(u95+u84)*z))
  u83=u93*u57
  value_1_=value_1_+(c(65)*(u62*(u83+u84)))
  value_1_=value_1_+(c(66)*(u64*(u83+u94)))
  u94=(u90+u78)
  u84=u80+u81*u60
  u83=u91*u60+u89
  value_1_=value_1_+(c(67)*(u56*(u84+u94*u56)+u83))
  u95=(u78+u90)
  u93=u76+u91
  value_1_=value_1_+(c(68)*(y*(u95*u56+u93)*z))
  u76=8.55176772024602675d-3*u75*(u92+u66*(4.0d0*pd2-3.0d0))
  u66=-u99
  u75=3.59174244250333123d-1*u97*(8.90641517292706383d3*u79-u87*pd2)
  u79=-8.97935610625832809d-2*u96
  u96=u79*u60
  value_1_=value_1_+(c(69)*(u56*(u78+u60*(u96+u75)+(u96+u66)*u56)+u60*(&
 &u78+u66*u60)+u76))
  value_1_=value_1_+(c(70)*(u63*(u95*u57+u93)))
  value_1_=value_1_+(c(71)*(u57*(u84+u94*u57)+u83))
  u60=u9*u56
  u75=(u56*(u85+u60)+u88)
  value_1_=value_1_+(c(72)*(u62*u75))
  value_1_=value_1_+(c(73)*(x*(u56*(u81+u60)+u91)*z))
  u84=(u77+u60)*u57
  u95=u78*u56
  value_1_=value_1_+(c(74)*(u62*(u84+u95+u91)))
  u83=u77*u56
  value_1_=value_1_+(c(75)*(u64*((u78+u60)*u57+u83+u91)))
  u94=u9*u57
  value_1_=value_1_+(c(76)*(u62*(u57*(u81+u94)+u91)))
  u76=(u57*(u85+u94)+u88)
  value_1_=value_1_+(c(77)*(u64*u76))
  value_1_=value_1_+(c(78)*(u56*(u8+u56*(u60+u98))+u86))
  value_1_=value_1_+(c(79)*(y*u75*z))
  value_1_=value_1_+(c(80)*((u91+u56*(u60+u81))*u57+u56*(u80+u95)+u89))
  value_1_=value_1_+(c(81)*(u63*(u84+u83+u74)))
  value_1_=value_1_+(c(82)*(u57*(u80+u81*u56+(u60+u78)*u57)+u91*u56+u89&
 &))
  value_1_=value_1_+(c(83)*(u63*u76))
  value_1_=value_1_+(c(84)*(u57*(u8+u57*(u94+u98))+u86))
  if ( lmax .eq. 6 ) go to 100

  ! set computation flag
  D_X_Y_shell_optv4_0 = .false.

  ! setting results
100  continue
  value(1:1+1*(NVEC-1):1) = value_1_

  contains

    include "Av_functions.h"

end function

  
end module
