module mod_D_X_Y_opt_1
  
  use mod_A_functions
  
  implicit none
  
  private
  public :: D_X_Y_opt_1
  
contains
  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_Y_opt_1_m1(a12,a3,r3,c,lmax,value)
  
  use mod_A_functions
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12      !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3       !< exponent of the third solid harmonics 
  real(kind=8), intent(in)  :: r3(3)    !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)     !< coeffcient of the hermit decomposition of the product of orbital 1 and 2
  integer     , intent(in)  :: lmax     !< l of maximum non 0 coeff in hermit decomposition~%
  real(kind=8), intent(out) :: value    !< result value
  
  ! return value
  
  logical :: D_X_Y_opt_1_m1
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u61
  real(kind=8) :: u62
  real(kind=8) :: u63
  real(kind=8) :: u64
  real(kind=8) :: u65
  real(kind=8) :: u66
  real(kind=8) :: u67
  real(kind=8) :: u68
  real(kind=8) :: u69
  real(kind=8) :: u7
  real(kind=8) :: u70
  real(kind=8) :: u71
  real(kind=8) :: u72
  real(kind=8) :: u73
  real(kind=8) :: u74
  real(kind=8) :: u75
  real(kind=8) :: u76
  real(kind=8) :: u77
  real(kind=8) :: u78
  real(kind=8) :: u79
  real(kind=8) :: u8
  real(kind=8) :: u80
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
  real(kind=8) :: u97
  real(kind=8) :: u98
  real(kind=8) :: u99

  
  ! init computation flag
  D_X_Y_opt_1_m1=.true.
  
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
  value=0.0d0
  u71=a3**(-1)
  u65=-2.44301255951459961d-1*u71*A2_(p,pd,erfpd,exppd2)*p**2
  value=value+(c(1)*(u65*y))
  if ( lmax .eq. 0 ) return
  u65=A3_(p,pd,erfpd,exppd2)
  u76=p**3
  u66=-7.32903767854379883d-1*u71*u65*u76
  u87=x*y
  value=value+(c(2)*(u66*u87))
  u88=2.59211683254887779d-2*u71*(9.42477796076937972d0*u65+7.089815403&
 &62206411d0*exppd2)*u76
  u68=y**2
  value=value+(c(3)*(u66*u68+u88))
  u88=y*z
  value=value+(c(4)*(u66*u88))
  if ( lmax .eq. 1 ) return
  u66=A4_(p,pd,erfpd,exppd2)
  u77=p**4
  u84=u71*u66*u77
  u98=7.32903767854379883d-1*u84
  u67=-7.77635049764663336d-2*u71
  u89=exppd2*pd
  u83=-1.41796308072441282d1*u89
  u91=+u67*(4.71238898038468986d1*u66+u83)*u77
  u69=x**2
  value=value+(c(5)*((u91*u69+u98)*y))
  u74=u91*u68
  u63=(u74+u98)
  value=value+(c(6)*(x*u63))
  u81=u87*z
  value=value+(c(7)*(u91*u81))
  u70=2.19871130356313965d0*u84
  value=value+(c(8)*(y*(u74+u70)))
  value=value+(c(9)*(u63*z))
  u70=z**2
  value=value+(c(10)*(y*(u91*u70+u98)))
  if ( lmax .eq. 2 ) return
  u63=A5_(p,pd,erfpd,exppd2)
  u74=p**5
  u65=u71*u63*u74
  u98=1.09935565178156983d1*u65
  u91=-2.83592616144882564d1*exppd2
  u84=+u91*pd2
  u94=+u67*(3.2986722862692829d2*u63+u84)*u74
  u62=u94*u69
  value=value+(c(11)*(x*(u62+u98)*y))
  u92=1.41796308072441282d1*exppd2
  u75=u71*(u92+4.71238898038468986d1*u63)*u74
  u78=-1.55527009952932667d-2*u75
  u90=3.66451883927189941d0*u65
  value=value+(c(12)*((u90+u62)*u68+u90*u69+u78))
  value=value+(c(13)*((u62+u90)*u88))
  u62=u94*u68
  u73=(u62+u98)
  value=value+(c(14)*(u87*u73))
  value=value+(c(15)*(x*(u62+u90)*z))
  u85=u94*u70
  value=value+(c(16)*(u87*(u85+u90)))
  u64=-4.66581029858798002d-2*u75
  u75=2.19871130356313965d1*u65
  value=value+(c(17)*(u68*(u75+u62)+u64))
  value=value+(c(18)*(y*u73*z))
  value=value+(c(19)*((u90+u62)*u70+u90*u68+u78))
  value=value+(c(20)*(u88*(u85+u98)))
  if ( lmax .eq. 3 ) return
  u64=A6_(p,pd,erfpd,exppd2)
  u75=p**6
  u73=u71*u64*u75
  u62=-1.09935565178156983d1*u73
  u85=-2.83592616144882564d1*u89
  u94=u71*(3.2986722862692829d2*u64+u85)*u75
  u78=4.66581029858798002d-1*u94
  u90=2.0d0*pd2
  u98=+u67*(2.96880505764235461d3*u64+u85*(9.0d0+u90))*u75
  u65=u98*u69
  value=value+(c(21)*((u69*(u78+u65)+u62)*y))
  u79=2.33290514929399001d-1*u94
  u63=7.77635049764663336d-2*u94
  u86=u63*u69
  value=value+(c(22)*(x*((u79+u65)*u68+u86+u62)))
  value=value+(c(23)*(x*(u65+u79)*u88))
  u96=(u63+u65)
  u65=u96*u68
  value=value+(c(24)*(y*(u65+u79*u69+u62)))
  u8=-3.66451883927189941d0*u73
  u67=u86+u8
  value=value+(c(25)*((u65+u67)*z))
  value=value+(c(26)*(y*(u96*u70+u67)))
  u65=u98*u68
  u96=(u68*(u78+u65)+u62)
  value=value+(c(27)*(x*u96))
  value=value+(c(28)*(u87*(u65+u79)*z))
  u67=(u63+u65)*u70
  u86=u63*u68
  value=value+(c(29)*(x*(u67+u86+u8)))
  u8=u98*u70
  value=value+(c(30)*(u81*(u8+u79)))
  u61=-5.49677825890784912d1*u73
  u73=7.77635049764663336d-1*u94
  value=value+(c(31)*(y*(u68*(u73+u65)+u61)))
  value=value+(c(32)*(u96*z))
  value=value+(c(33)*(y*((u79+u65)*u70+u86+u62)))
  value=value+(c(34)*(z*(u67+u79*u68+u62)))
  value=value+(c(35)*(y*(u70*(u78+u8)+u62)))
  if ( lmax .eq. 4 ) return
  u61=A7_(p,pd,erfpd,exppd2)
  u73=p**7
  u67=u71*u61*u73
  u96=-3.84774478123549439d2*u67
  u86=-5.67185232289765129d1*exppd2*pd2
  u8=u71*(2.96880505764235461d3*u61+u86)*u73
  u63=7.77635049764663336d-1*u8
  u79=(1.1d1+u90)
  u78=u71*(3.26568556340659007d4*u61+u86*u79)*u73
  u94=-7.77635049764663336d-2*u78
  u62=u94*u69
  value=value+(c(36)*(x*(u69*(u63+u62)+u96)*y))
  u65=u71*(3.2986722862692829d2*u61-u91)*u73
  u98=3.33272164184855715d-2*u65
  u91=-7.69548956247098877d1*u67
  u74=-1.53909791249419775d2*u67
  u9=4.66581029858798002d-1*u8
  u80=7.77635049764663336d-2*u8
  u97=u80*u69
  value=value+(c(37)*((u91+u69*(u62+u9))*u68+u69*(u74+u97)+u98))
  value=value+(c(38)*((u69*(u9+u62)+u91)*u88))
  u64=-2.30864686874129663d2*u67
  u7=2.33290514929399001d-1*u8
  u72=(u7+u62)
  u75=u72*u68
  u95=u7*u69
  value=value+(c(39)*(u87*(u75+u95+u64)))
  u82=u97+u91
  value=value+(c(40)*(x*(u75+u82)*z))
  value=value+(c(41)*(u87*(u72*u70+u82)))
  value=value+(c(42)*(u68*(u74+u9*u69+(u62+u80)*u68)+u91*u69+u98))
  u75=(u80+u62)
  u72=u95+u91
  value=value+(c(43)*(y*(u75*u68+u72)*z))
  u82=-2.2218144278990381d-2*u71*(9.8960168588078487d2*u61-u92)*u73
  u62=-7.77635049764663336d-2*u8
  u97=-3.11054019905865334d-1*u71*(8.90641517292706383d3*u61+u84*(6.0d0&
 &+pd2))*u73
  u84=7.77635049764663336d-2*u78
  u78=u84*u69
  value=value+(c(44)*(u68*(u80+u69*(u78+u97)+(u78+u62)*u68)+u69*(u80+u6&
 &2*u69)+u82))
  value=value+(c(45)*(u88*(u75*u70+u72)))
  u82=u94*u68
  u97=(u68*(u63+u82)+u96)
  value=value+(c(46)*(u87*u97))
  value=value+(c(47)*(x*(u68*(u9+u82)+u91)*z))
  u75=(u7+u82)*u70
  u84=u80*u68
  value=value+(c(48)*(u87*(u75+u84+u91)))
  u78=u7*u68
  value=value+(c(49)*(x*z*((u80+u82)*u70+u78+u91)))
  u95=u94*u70
  value=value+(c(50)*(u87*(u70*(u9+u95)+u91)))
  u62=1.66636082092427858d-1*u65
  u72=-1.15432343437064832d3*u67
  u65=1.166452574646995d0*u8
  value=value+(c(51)*(u68*(u72+u68*(u82+u65))+u62))
  value=value+(c(52)*(y*u97*z))
  value=value+(c(53)*((u91+u68*(u82+u9))*u70+u68*(u74+u84)+u98))
  value=value+(c(54)*(u88*(u75+u78+u64)))
  value=value+(c(55)*(u70*(u74+u9*u68+(u82+u80)*u70)+u91*u68+u98))
  value=value+(c(56)*(u88*(u70*(u63+u95)+u96)))
  if ( lmax .eq. 5 ) return
  u62=A8_(p,pd,erfpd,exppd2)
  u72=p**8
  u65=u71*u62*u72
  u75=3.84774478123549439d2*u65
  u82=-5.67185232289765129d1*u89
  u84=u71*(2.96880505764235461d3*u62+u82)*u72
  u64=-3.49935772394098501d0*u84
  u95=3.26568556340659007d4*u62
  u74=u71*(u95+u82*u79)*u72
  u73=1.166452574646995d0*u74
  u78=pd2**2
  u80=(u90+1.1d1)
  u7=u71*(4.24539123242856709d5*u62+u82*(1.3d1*u80+4.0d0*u78))*u72
  u79=-7.77635049764663336d-2*u7
  u77=u79*u69
  value=value+(c(57)*((u69*(u64+u69*(u77+u73))+u75)*y))
  u91=-1.166452574646995d0*u84
  u96=-7.77635049764663336d-1*u84
  u86=7.77635049764663336d-1*u74
  u9=7.77635049764663336d-2*u74
  u97=u9*u69
  value=value+(c(58)*(x*((u91+u69*(u77+u86))*u68+u69*(u96+u97)+u75)))
  value=value+(c(59)*(x*(u69*(u86+u77)+u91)*u88))
  u63=2.30864686874129663d2*u65
  u61=-2.33290514929399001d-1*u84
  u67=-1.39974308957639401d0*u84
  u88=4.66581029858798002d-1*u74
  u8=2.33290514929399001d-1*u74
  u98=(u61+u69*(u77+u88))
  u94=u98*u68
  u92=u8*u69
  value=value+(c(60)*(y*(u94+u69*(u67+u92)+u63)))
  u66=7.69548956247098877d1*u65
  u85=-4.66581029858798002d-1*u84
  u99=u69*(u85+u97)+u66
  value=value+(c(61)*((u94+u99)*z))
  value=value+(c(62)*(y*(u98*u70+u99)))
  u94=u88*u69
  u99=u61*u69
  value=value+(c(63)*(x*(u68*(u67+u94+(u77+u8)*u68)+u99+u63)))
  u97=-6.99871544788197003d-1*u84
  u98=(u8+u77)
  u93=u92+u97
  value=value+(c(64)*(u87*(u98*u68+u93)*z))
  u92=-1.86632411943519201d0*u71*(3.2986722862692829d2*u62-7.0898154036&
 &2206411d0*u89)*u72
  u76=-2.33290514929399001d-1*u74
  u89=7.77635049764663336d-2*u71*(3.85944657493506099d4*u62+u82*(1.3d1+&
 &u90))*u72
  u62=-1.24421607962346134d0*u71*(u95+u83*(4.0d0*u80+u78))*u72
  u72=7.77635049764663336d-2*u7
  u7=-7.77635049764663336d-2*u74
  u78=u72*u69
  u95=u69*(u78+u62)
  value=value+(c(65)*(x*(u68*(u8+u95+(u78+u76)*u68)+u69*(u89+u7*u69)+u9&
 &2)))
  value=value+(c(66)*(u81*(u98*u70+u93)))
  u80=(u77+u9)
  u93=u80*u68
  value=value+(c(67)*(y*(u68*(u96+u86*u69+u93)+u91*u69+u75)))
  u71=u85+u94
  u94=u99+u66
  value=value+(c(68)*((u68*(u71+u93)+u94)*z))
  u93=u89+u95
  u95=u69*(u8+u76*u69)+u92
  u98=(u78+u7)
  value=value+(c(69)*(y*(u68*(u93+u98*u68)+u95)))
  value=value+(c(70)*(z*(u70*(u93+u98*u70)+u95)))
  value=value+(c(71)*(y*(u70*(u71+u80*u70)+u94)))
  u69=u79*u68
  u93=(u68*(u64+u68*(u69+u73))+u75)
  value=value+(c(72)*(x*u93))
  value=value+(c(73)*(u87*(u68*(u86+u69)+u91)*z))
  u95=(u61+u68*(u69+u88))*u70
  u87=u9*u68
  value=value+(c(74)*(x*(u95+u68*(u85+u87)+u66)))
  u80=u8*u68
  value=value+(c(75)*(u81*((u8+u69)*u70+u80+u97)))
  u97=(u69+u9)*u70
  u94=u88*u68
  u71=u61*u68
  value=value+(c(76)*(x*(u70*(u85+u94+u97)+u71+u66)))
  u66=u79*u70
  value=value+(c(77)*(u81*(u70*(u86+u66)+u91)))
  u81=2.69342134686484607d3*u65
  u65=-8.16516802252896503d0*u84
  u84=1.63303360450579301d0*u74
  value=value+(c(78)*(y*(u68*(u65+u68*(u69+u84))+u81)))
  value=value+(c(79)*(u93*z))
  value=value+(c(80)*(y*((u91+u68*(u69+u86))*u70+u68*(u96+u87)+u75)))
  value=value+(c(81)*(z*(u95+u68*(u67+u80)+u63)))
  value=value+(c(82)*(y*(u70*(u67+u94+(u69+u8)*u70)+u71+u63)))
  value=value+(c(83)*(z*(u70*(u96+u86*u68+u97)+u91*u68+u75)))
  value=value+(c(84)*(y*(u70*(u64+u70*(u66+u73))+u75)))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_opt_1_m1=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_Y_opt_1_0(a12,a3,r3,c,lmax,value)
  
  use mod_A_functions
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12      !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3       !< exponent of the third solid harmonics 
  real(kind=8), intent(in)  :: r3(3)    !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)     !< coeffcient of the hermit decomposition of the product of orbital 1 and 2
  integer     , intent(in)  :: lmax     !< l of maximum non 0 coeff in hermit decomposition~%
  real(kind=8), intent(out) :: value    !< result value
  
  ! return value
  
  logical :: D_X_Y_opt_1_0
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u63
  real(kind=8) :: u64
  real(kind=8) :: u65
  real(kind=8) :: u66
  real(kind=8) :: u67
  real(kind=8) :: u68
  real(kind=8) :: u69
  real(kind=8) :: u7
  real(kind=8) :: u70
  real(kind=8) :: u71
  real(kind=8) :: u72
  real(kind=8) :: u73
  real(kind=8) :: u74
  real(kind=8) :: u75
  real(kind=8) :: u76
  real(kind=8) :: u77
  real(kind=8) :: u78
  real(kind=8) :: u79
  real(kind=8) :: u8
  real(kind=8) :: u80
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
  real(kind=8) :: u97
  real(kind=8) :: u98
  real(kind=8) :: u99

  
  ! init computation flag
  D_X_Y_opt_1_0=.true.
  
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
  value=0.0d0
  u71=a3**(-1)
  u65=-2.44301255951459961d-1*u71*A2_(p,pd,erfpd,exppd2)*p**2
  value=value+(c(1)*(u65*z))
  if ( lmax .eq. 0 ) return
  u65=A3_(p,pd,erfpd,exppd2)
  u76=p**3
  u77=-7.32903767854379883d-1*u71*u65*u76
  u89=x*z
  value=value+(c(2)*(u77*u89))
  u87=y*z
  value=value+(c(3)*(u77*u87))
  u66=2.59211683254887779d-2*u71*(9.42477796076937972d0*u65+7.089815403&
 &62206411d0*exppd2)*u76
  u69=z**2
  value=value+(c(4)*(u77*u69+u66))
  if ( lmax .eq. 1 ) return
  u66=A4_(p,pd,erfpd,exppd2)
  u77=p**4
  u92=u71*u66*u77
  u91=7.32903767854379883d-1*u92
  u67=-7.77635049764663336d-2*u71
  u90=exppd2*pd
  u83=-1.41796308072441282d1*u90
  u86=+u67*(4.71238898038468986d1*u66+u83)*u77
  u68=x**2
  value=value+(c(5)*((u86*u68+u91)*z))
  u81=x*u87
  value=value+(c(6)*(u86*u81))
  u74=u86*u69
  u63=(u74+u91)
  value=value+(c(7)*(x*u63))
  u70=y**2
  value=value+(c(8)*((u86*u70+u91)*z))
  value=value+(c(9)*(y*u63))
  u63=2.19871130356313965d0*u92
  value=value+(c(10)*(z*(u74+u63)))
  if ( lmax .eq. 2 ) return
  u63=A5_(p,pd,erfpd,exppd2)
  u74=p**5
  u86=u71*u63*u74
  u91=1.09935565178156983d1*u86
  u92=-2.83592616144882564d1*exppd2
  u84=+u92*pd2
  u94=+u67*(3.2986722862692829d2*u63+u84)*u74
  u75=u94*u68
  value=value+(c(11)*(x*(u75+u91)*z))
  u9=3.66451883927189941d0*u86
  value=value+(c(12)*((u75+u9)*u87))
  u93=1.41796308072441282d1*exppd2
  u79=u71*(u93+4.71238898038468986d1*u63)*u74
  u85=-1.55527009952932667d-2*u79
  value=value+(c(13)*((u9+u75)*u69+u9*u68+u85))
  u75=u94*u70
  value=value+(c(14)*(x*(u75+u9)*z))
  u88=x*y
  u7=u94*u69
  value=value+(c(15)*(u88*(u7+u9)))
  u64=(u7+u91)
  value=value+(c(16)*(u89*u64))
  value=value+(c(17)*(y*(u75+u91)*z))
  value=value+(c(18)*((u9+u75)*u69+u9*u70+u85))
  value=value+(c(19)*(u87*u64))
  u64=-4.66581029858798002d-2*u79
  u75=2.19871130356313965d1*u86
  value=value+(c(20)*(u69*(u75+u7)+u64))
  if ( lmax .eq. 3 ) return
  u64=A6_(p,pd,erfpd,exppd2)
  u75=p**6
  u7=u71*u64*u75
  u63=-1.09935565178156983d1*u7
  u85=-2.83592616144882564d1*u90
  u79=u71*(3.2986722862692829d2*u64+u85)*u75
  u9=4.66581029858798002d-1*u79
  u91=2.0d0*pd2
  u94=+u67*(2.96880505764235461d3*u64+u85*(9.0d0+u91))*u75
  u74=u94*u68
  value=value+(c(21)*((u68*(u9+u74)+u63)*z))
  u67=2.33290514929399001d-1*u79
  value=value+(c(22)*(x*(u74+u67)*u87))
  u86=7.77635049764663336d-2*u79
  u73=u86*u68
  value=value+(c(23)*(x*((u67+u74)*u69+u73+u63)))
  u75=-3.66451883927189941d0*u7
  u98=(u86+u74)
  u74=u73+u75
  value=value+(c(24)*((u98*u70+u74)*z))
  u73=u98*u69
  value=value+(c(25)*(y*(u73+u74)))
  value=value+(c(26)*(z*(u73+u67*u68+u63)))
  u73=u94*u70
  value=value+(c(27)*(u88*(u73+u67)*z))
  u98=(u86+u73)*u69
  u74=u86*u70
  value=value+(c(28)*(x*(u98+u74+u75)))
  u75=u94*u69
  value=value+(c(29)*(u81*(u75+u67)))
  u97=(u69*(u9+u75)+u63)
  value=value+(c(30)*(x*u97))
  value=value+(c(31)*((u70*(u9+u73)+u63)*z))
  value=value+(c(32)*(y*((u67+u73)*u69+u74+u63)))
  value=value+(c(33)*(z*(u98+u67*u70+u63)))
  value=value+(c(34)*(y*u97))
  u98=-5.49677825890784912d1*u7
  u73=7.77635049764663336d-1*u79
  value=value+(c(35)*(z*(u69*(u73+u75)+u98)))
  if ( lmax .eq. 4 ) return
  u98=A7_(p,pd,erfpd,exppd2)
  u73=p**7
  u75=u71*u98*u73
  u74=-3.84774478123549439d2*u75
  u86=-5.67185232289765129d1*exppd2*pd2
  u67=u71*(2.96880505764235461d3*u98+u86)*u73
  u9=7.77635049764663336d-1*u67
  u79=(1.1d1+u91)
  u63=u71*(3.26568556340659007d4*u98+u86*u79)*u73
  u97=-7.77635049764663336d-2*u63
  u64=u97*u68
  value=value+(c(36)*(x*(u68*(u9+u64)+u74)*z))
  u94=-7.69548956247098877d1*u75
  u7=4.66581029858798002d-1*u67
  value=value+(c(37)*((u68*(u7+u64)+u94)*u87))
  u8=u71*(3.2986722862692829d2*u98-u92)*u73
  u99=3.33272164184855715d-2*u8
  u92=-1.53909791249419775d2*u75
  u78=7.77635049764663336d-2*u67
  u82=u78*u68
  value=value+(c(38)*((u94+u68*(u64+u7))*u69+u68*(u92+u82)+u99))
  u80=2.33290514929399001d-1*u67
  u66=(u80+u64)
  u72=u82+u94
  value=value+(c(39)*(x*(u66*u70+u72)*z))
  u82=u66*u69
  value=value+(c(40)*(u88*(u82+u72)))
  u66=-2.30864686874129663d2*u75
  u72=u80*u68
  value=value+(c(41)*(u89*(u82+u72+u66)))
  u82=(u78+u64)
  u95=u72+u94
  value=value+(c(42)*(y*(u82*u70+u95)*z))
  u72=-2.2218144278990381d-2*u71*(9.8960168588078487d2*u98-u93)*u73
  u96=-7.77635049764663336d-2*u67
  u76=-3.11054019905865334d-1*u71*(8.90641517292706383d3*u98+u84*(6.0d0&
 &+pd2))*u73
  u98=7.77635049764663336d-2*u63
  u63=u98*u68
  value=value+(c(43)*(u70*(u78+u68*(u63+u76)+(u63+u96)*u70)+u68*(u78+u9&
 &6*u68)+u72))
  value=value+(c(44)*(u87*(u82*u69+u95)))
  value=value+(c(45)*(u69*(u92+u7*u68+(u64+u78)*u69)+u94*u68+u99))
  u72=u97*u70
  value=value+(c(46)*(x*(u70*(u7+u72)+u94)*z))
  u76=(u80+u72)*u69
  u63=u78*u70
  value=value+(c(47)*(u88*(u76+u63+u94)))
  u82=u80*u70
  value=value+(c(48)*(u89*((u78+u72)*u69+u82+u94)))
  u64=u97*u69
  value=value+(c(49)*(u88*(u69*(u7+u64)+u94)))
  u98=(u69*(u9+u64)+u74)
  value=value+(c(50)*(u89*u98))
  value=value+(c(51)*(y*(u70*(u9+u72)+u74)*z))
  value=value+(c(52)*((u94+u70*(u72+u7))*u69+u70*(u92+u63)+u99))
  value=value+(c(53)*(u87*(u76+u82+u66)))
  value=value+(c(54)*(u69*(u92+u7*u70+(u72+u78)*u69)+u94*u70+u99))
  value=value+(c(55)*(u87*u98))
  u76=1.66636082092427858d-1*u8
  u72=-1.15432343437064832d3*u75
  u63=1.166452574646995d0*u67
  value=value+(c(56)*(u69*(u72+u69*(u64+u63))+u76))
  if ( lmax .eq. 5 ) return
  u76=A8_(p,pd,erfpd,exppd2)
  u72=p**8
  u63=u71*u76*u72
  u64=3.84774478123549439d2*u63
  u82=-5.67185232289765129d1*u90
  u65=u71*(2.96880505764235461d3*u76+u82)*u72
  u66=-3.49935772394098501d0*u65
  u96=3.26568556340659007d4*u76
  u74=u71*(u96+u82*u79)*u72
  u93=1.166452574646995d0*u74
  u78=pd2**2
  u80=(u91+1.1d1)
  u7=u71*(4.24539123242856709d5*u76+u82*(1.3d1*u80+4.0d0*u78))*u72
  u79=-7.77635049764663336d-2*u7
  u92=u79*u68
  value=value+(c(57)*((u68*(u66+u68*(u92+u93))+u64)*z))
  u67=-1.166452574646995d0*u65
  u73=7.77635049764663336d-1*u74
  value=value+(c(58)*(x*(u68*(u73+u92)+u67)*u87))
  u97=-7.77635049764663336d-1*u65
  u9=7.77635049764663336d-2*u74
  u95=u9*u68
  value=value+(c(59)*(x*((u67+u68*(u92+u73))*u69+u68*(u97+u95)+u64)))
  u75=7.69548956247098877d1*u63
  u89=-2.33290514929399001d-1*u65
  u98=-4.66581029858798002d-1*u65
  u87=4.66581029858798002d-1*u74
  u99=(u89+u68*(u92+u87))
  u85=u68*(u98+u95)+u75
  value=value+(c(60)*((u99*u70+u85)*z))
  u95=u99*u69
  value=value+(c(61)*(y*(u95+u85)))
  u85=2.30864686874129663d2*u63
  u99=-1.39974308957639401d0*u65
  u8=2.33290514929399001d-1*u74
  u86=u8*u68
  value=value+(c(62)*(z*(u95+u68*(u99+u86)+u85)))
  u84=-6.99871544788197003d-1*u65
  u94=(u8+u92)
  u95=u86+u84
  value=value+(c(63)*(u88*(u94*u70+u95)*z))
  u86=-1.86632411943519201d0*u71*(3.2986722862692829d2*u76-7.0898154036&
 &2206411d0*u90)*u72
  u77=-2.33290514929399001d-1*u74
  u90=7.77635049764663336d-2*u71*(3.85944657493506099d4*u76+u82*(1.3d1+&
 &u91))*u72
  u76=-1.24421607962346134d0*u71*(u96+u83*(4.0d0*u80+u78))*u72
  u71=7.77635049764663336d-2*u7
  u7=-7.77635049764663336d-2*u74
  u80=u71*u68
  u96=u68*(u80+u76)
  value=value+(c(64)*(x*(u70*(u8+u96+(u80+u77)*u70)+u68*(u90+u7*u68)+u8&
 &6)))
  value=value+(c(65)*(u81*(u94*u69+u95)))
  u94=u87*u68
  u78=u89*u68
  value=value+(c(66)*(x*(u69*(u99+u94+(u92+u8)*u69)+u78+u85)))
  u72=(u92+u9)
  u95=u98+u94
  u92=u78+u75
  value=value+(c(67)*((u70*(u95+u72*u70)+u92)*z))
  u94=u90+u96
  u96=u68*(u8+u77*u68)+u86
  u78=(u80+u7)
  value=value+(c(68)*(y*(u70*(u94+u78*u70)+u96)))
  value=value+(c(69)*(z*(u69*(u94+u78*u69)+u96)))
  u96=u72*u69
  value=value+(c(70)*(y*(u69*(u95+u96)+u92)))
  value=value+(c(71)*(z*(u69*(u97+u73*u68+u96)+u67*u68+u64)))
  u68=u79*u70
  value=value+(c(72)*(u88*(u70*(u73+u68)+u67)*z))
  u96=(u89+u70*(u68+u87))*u69
  u88=u9*u70
  value=value+(c(73)*(x*(u96+u70*(u98+u88)+u75)))
  u92=u8*u70
  value=value+(c(74)*(u81*((u8+u68)*u69+u92+u84)))
  u84=(u68+u9)*u69
  u72=u87*u70
  u95=u89*u70
  value=value+(c(75)*(x*(u69*(u98+u72+u84)+u95+u75)))
  u75=u79*u69
  value=value+(c(76)*(u81*(u69*(u73+u75)+u67)))
  u94=(u69*(u66+u69*(u75+u93))+u64)
  value=value+(c(77)*(x*u94))
  value=value+(c(78)*((u70*(u66+u70*(u68+u93))+u64)*z))
  value=value+(c(79)*(y*((u67+u70*(u68+u73))*u69+u70*(u97+u88)+u64)))
  value=value+(c(80)*(z*(u96+u70*(u99+u92)+u85)))
  value=value+(c(81)*(y*(u69*(u99+u72+(u68+u8)*u69)+u95+u85)))
  value=value+(c(82)*(z*(u69*(u97+u73*u70+u84)+u67*u70+u64)))
  value=value+(c(83)*(y*u94))
  u64=2.69342134686484607d3*u63
  u63=-8.16516802252896503d0*u65
  u65=1.63303360450579301d0*u74
  value=value+(c(84)*(z*(u69*(u63+u69*(u75+u65))+u64)))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_opt_1_0=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_Y_opt_1_1(a12,a3,r3,c,lmax,value)
  
  use mod_A_functions
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12      !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3       !< exponent of the third solid harmonics 
  real(kind=8), intent(in)  :: r3(3)    !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)     !< coeffcient of the hermit decomposition of the product of orbital 1 and 2
  integer     , intent(in)  :: lmax     !< l of maximum non 0 coeff in hermit decomposition~%
  real(kind=8), intent(out) :: value    !< result value
  
  ! return value
  
  logical :: D_X_Y_opt_1_1
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u60
  real(kind=8) :: u61
  real(kind=8) :: u62
  real(kind=8) :: u63
  real(kind=8) :: u64
  real(kind=8) :: u65
  real(kind=8) :: u66
  real(kind=8) :: u67
  real(kind=8) :: u68
  real(kind=8) :: u69
  real(kind=8) :: u7
  real(kind=8) :: u70
  real(kind=8) :: u71
  real(kind=8) :: u72
  real(kind=8) :: u73
  real(kind=8) :: u74
  real(kind=8) :: u75
  real(kind=8) :: u76
  real(kind=8) :: u77
  real(kind=8) :: u78
  real(kind=8) :: u79
  real(kind=8) :: u8
  real(kind=8) :: u80
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
  real(kind=8) :: u97
  real(kind=8) :: u98
  real(kind=8) :: u99

  
  ! init computation flag
  D_X_Y_opt_1_1=.true.
  
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
  value=0.0d0
  u71=a3**(-1)
  u65=-2.44301255951459961d-1*u71*A2_(p,pd,erfpd,exppd2)*p**2
  value=value+(c(1)*(u65*x))
  if ( lmax .eq. 0 ) return
  u65=A3_(p,pd,erfpd,exppd2)
  u76=p**3
  u87=2.59211683254887779d-2*u71*(9.42477796076937972d0*u65+7.089815403&
 &62206411d0*exppd2)*u76
  u66=-7.32903767854379883d-1*u71*u65*u76
  u68=x**2
  value=value+(c(2)*(u66*u68+u87))
  u87=x*y
  value=value+(c(3)*(u66*u87))
  u88=x*z
  value=value+(c(4)*(u66*u88))
  if ( lmax .eq. 1 ) return
  u66=A4_(p,pd,erfpd,exppd2)
  u77=p**4
  u69=u71*u66*u77
  u74=2.19871130356313965d0*u69
  u67=-7.77635049764663336d-2*u71
  u90=exppd2*pd
  u83=-1.41796308072441282d1*u90
  u63=+u67*(4.71238898038468986d1*u66+u83)*u77
  u81=u63*u68
  value=value+(c(5)*(x*(u81+u74)))
  u74=7.32903767854379883d-1*u69
  u69=(u81+u74)
  value=value+(c(6)*(u69*y))
  value=value+(c(7)*(u69*z))
  u69=y**2
  value=value+(c(8)*(x*(u63*u69+u74)))
  u81=u87*z
  value=value+(c(9)*(u63*u81))
  u70=z**2
  value=value+(c(10)*(x*(u63*u70+u74)))
  if ( lmax .eq. 2 ) return
  u63=A5_(p,pd,erfpd,exppd2)
  u74=p**5
  u92=1.41796308072441282d1*exppd2
  u7=u71*(u92+4.71238898038468986d1*u63)*u74
  u8=-4.66581029858798002d-2*u7
  u75=u71*u63*u74
  u85=2.19871130356313965d1*u75
  u93=-2.83592616144882564d1*exppd2
  u84=+u93*pd2
  u79=+u67*(3.2986722862692829d2*u63+u84)*u74
  u64=u79*u68
  value=value+(c(11)*(u68*(u85+u64)+u8))
  u8=1.09935565178156983d1*u75
  u85=x*(u64+u8)
  value=value+(c(12)*(u85*y))
  value=value+(c(13)*(u85*z))
  u85=-1.55527009952932667d-2*u7
  u7=3.66451883927189941d0*u75
  u75=(u7+u64)
  u73=u7*u68+u85
  value=value+(c(14)*(u75*u69+u73))
  u89=y*z
  value=value+(c(15)*((u64+u7)*u89))
  value=value+(c(16)*(u75*u70+u73))
  u64=u79*u69
  value=value+(c(17)*(u87*(u64+u8)))
  value=value+(c(18)*(x*(u64+u7)*z))
  u64=u79*u70
  value=value+(c(19)*(u87*(u64+u7)))
  value=value+(c(20)*(u88*(u64+u8)))
  if ( lmax .eq. 3 ) return
  u64=A6_(p,pd,erfpd,exppd2)
  u75=p**6
  u73=u71*u64*u75
  u79=-5.49677825890784912d1*u73
  u85=-2.83592616144882564d1*u90
  u7=u71*(3.2986722862692829d2*u64+u85)*u75
  u8=7.77635049764663336d-1*u7
  u91=2.0d0*pd2
  u82=+u67*(2.96880505764235461d3*u64+u85*(9.0d0+u91))*u75
  u61=u82*u68
  value=value+(c(21)*(x*(u68*(u8+u61)+u79)))
  u79=-1.09935565178156983d1*u73
  u8=4.66581029858798002d-1*u7
  u98=(u68*(u8+u61)+u79)
  value=value+(c(22)*(u98*y))
  value=value+(c(23)*(u98*z))
  u86=2.33290514929399001d-1*u7
  u97=7.77635049764663336d-2*u7
  u7=(u86+u61)
  u60=u97*u68
  u95=u60+u79
  value=value+(c(24)*(x*(u7*u69+u95)))
  value=value+(c(25)*(x*(u61+u86)*u89))
  value=value+(c(26)*(x*(u7*u70+u95)))
  u7=(u97+u61)
  u95=u7*u69
  u61=u86*u68+u79
  value=value+(c(27)*(y*(u95+u61)))
  u96=-3.66451883927189941d0*u73
  u73=u60+u96
  value=value+(c(28)*((u95+u73)*z))
  u77=u7*u70
  value=value+(c(29)*(y*(u77+u73)))
  value=value+(c(30)*(z*(u77+u61)))
  u61=u82*u69
  value=value+(c(31)*(x*(u69*(u8+u61)+u79)))
  value=value+(c(32)*(u87*(u61+u86)*z))
  value=value+(c(33)*(x*((u97+u61)*u70+u97*u69+u96)))
  u61=u82*u70
  value=value+(c(34)*(u81*(u61+u86)))
  value=value+(c(35)*(x*(u70*(u8+u61)+u79)))
  if ( lmax .eq. 4 ) return
  u61=A7_(p,pd,erfpd,exppd2)
  u73=p**7
  u7=u71*(3.2986722862692829d2*u61-u93)*u73
  u77=1.66636082092427858d-1*u7
  u9=u71*u61*u73
  u97=-1.15432343437064832d3*u9
  u86=-5.67185232289765129d1*exppd2*pd2
  u8=u71*(2.96880505764235461d3*u61+u86)*u73
  u95=1.166452574646995d0*u8
  u79=(1.1d1+u91)
  u60=u71*(3.26568556340659007d4*u61+u86*u79)*u73
  u96=-7.77635049764663336d-2*u60
  u82=u96*u68
  value=value+(c(36)*(u68*(u97+u68*(u82+u95))+u77))
  u77=-3.84774478123549439d2*u9
  u97=7.77635049764663336d-1*u8
  u95=x*(u68*(u97+u82)+u77)
  value=value+(c(37)*(u95*y))
  value=value+(c(38)*(u95*z))
  u86=3.33272164184855715d-2*u7
  u7=-7.69548956247098877d1*u9
  u98=-1.53909791249419775d2*u9
  u85=4.66581029858798002d-1*u8
  u94=7.77635049764663336d-2*u8
  u99=(u7+u68*(u82+u85))
  u80=u94*u68
  u67=u68*(u98+u80)+u86
  value=value+(c(39)*(u99*u69+u67))
  value=value+(c(40)*((u68*(u85+u82)+u7)*u89))
  value=value+(c(41)*(u99*u70+u67))
  u67=-2.30864686874129663d2*u9
  u9=2.33290514929399001d-1*u8
  u93=(u9+u82)
  u62=u93*u69
  u78=u9*u68
  u75=u78+u67
  value=value+(c(42)*(u87*(u62+u75)))
  u72=u80+u7
  value=value+(c(43)*(x*(u62+u72)*z))
  u62=u93*u70
  value=value+(c(44)*(u87*(u62+u72)))
  value=value+(c(45)*(u88*(u62+u75)))
  u62=(u82+u94)
  u93=u98+u85*u68
  u75=u7*u68+u86
  value=value+(c(46)*(u69*(u93+u62*u69)+u75))
  u72=(u94+u82)
  u82=u78+u7
  value=value+(c(47)*(y*(u72*u69+u82)*z))
  u80=-2.2218144278990381d-2*u71*(9.8960168588078487d2*u61-u92)*u73
  u64=-7.77635049764663336d-2*u8
  u78=-3.11054019905865334d-1*u71*(8.90641517292706383d3*u61+u84*(6.0d0&
 &+pd2))*u73
  u67=7.77635049764663336d-2*u60
  u60=u67*u68
  value=value+(c(48)*(u69*(u94+u68*(u60+u78)+(u60+u64)*u69)+u68*(u94+u6&
 &4*u68)+u80))
  value=value+(c(49)*(u89*(u72*u70+u82)))
  value=value+(c(50)*(u70*(u93+u62*u70)+u75))
  u62=u96*u69
  value=value+(c(51)*(u87*(u69*(u97+u62)+u77)))
  value=value+(c(52)*(x*(u69*(u85+u62)+u7)*z))
  value=value+(c(53)*(u87*((u9+u62)*u70+u94*u69+u7)))
  value=value+(c(54)*(u88*((u94+u62)*u70+u9*u69+u7)))
  u62=u96*u70
  value=value+(c(55)*(u87*(u70*(u85+u62)+u7)))
  value=value+(c(56)*(u88*(u70*(u97+u62)+u77)))
  if ( lmax .eq. 5 ) return
  u62=A8_(p,pd,erfpd,exppd2)
  u72=p**8
  u93=u71*u62*u72
  u75=2.69342134686484607d3*u93
  u82=-5.67185232289765129d1*u90
  u67=u71*(2.96880505764235461d3*u62+u82)*u72
  u64=-8.16516802252896503d0*u67
  u60=3.26568556340659007d4*u62
  u77=u71*(u60+u82*u79)*u72
  u98=1.63303360450579301d0*u77
  u78=pd2**2
  u80=(u91+1.1d1)
  u7=u71*(4.24539123242856709d5*u62+u82*(1.3d1*u80+4.0d0*u78))*u72
  u79=-7.77635049764663336d-2*u7
  u85=u79*u68
  value=value+(c(57)*(x*(u68*(u64+u68*(u85+u98))+u75)))
  u75=3.84774478123549439d2*u93
  u64=-3.49935772394098501d0*u67
  u98=1.166452574646995d0*u77
  u94=(u68*(u64+u68*(u85+u98))+u75)
  value=value+(c(58)*(u94*y))
  value=value+(c(59)*(u94*z))
  u74=-1.166452574646995d0*u67
  u86=-7.77635049764663336d-1*u67
  u94=7.77635049764663336d-1*u77
  u9=7.77635049764663336d-2*u77
  u63=(u74+u68*(u85+u94))
  u97=u9*u68
  u92=u68*(u86+u97)+u75
  value=value+(c(60)*(x*(u63*u69+u92)))
  value=value+(c(61)*(x*(u68*(u94+u85)+u74)*u89))
  value=value+(c(62)*(x*(u63*u70+u92)))
  u63=2.30864686874129663d2*u93
  u92=-2.33290514929399001d-1*u67
  u61=-1.39974308957639401d0*u67
  u89=4.66581029858798002d-1*u77
  u8=2.33290514929399001d-1*u77
  u99=(u92+u68*(u85+u89))
  u96=u99*u69
  u84=u8*u68
  u88=u68*(u61+u84)+u63
  value=value+(c(63)*(y*(u96+u88)))
  u66=7.69548956247098877d1*u93
  u93=-4.66581029858798002d-1*u67
  u65=u68*(u93+u97)+u66
  value=value+(c(64)*((u96+u65)*z))
  u97=u99*u70
  value=value+(c(65)*(y*(u97+u65)))
  value=value+(c(66)*(z*(u97+u88)))
  u97=(u85+u8)
  u88=u89*u68
  u96=u61+u88
  u99=u92*u68
  u65=u99+u63
  value=value+(c(67)*(x*(u69*(u96+u97*u69)+u65)))
  u73=-6.99871544788197003d-1*u67
  u67=(u8+u85)
  u95=u84+u73
  value=value+(c(68)*(u87*(u67*u69+u95)*z))
  u84=-1.86632411943519201d0*u71*(3.2986722862692829d2*u62-7.0898154036&
 &2206411d0*u90)*u72
  u76=-2.33290514929399001d-1*u77
  u90=7.77635049764663336d-2*u71*(3.85944657493506099d4*u62+u82*(1.3d1+&
 &u91))*u72
  u62=-1.24421607962346134d0*u71*(u60+u83*(4.0d0*u80+u78))*u72
  u60=7.77635049764663336d-2*u7
  u7=-7.77635049764663336d-2*u77
  u71=u60*u68
  u77=u68*(u71+u62)
  value=value+(c(69)*(x*(u69*(u8+u77+(u71+u76)*u69)+u68*(u90+u7*u68)+u8&
 &4)))
  value=value+(c(70)*(u81*(u67*u70+u95)))
  value=value+(c(71)*(x*(u70*(u96+u97*u70)+u65)))
  u61=(u85+u9)
  u65=u61*u69
  u85=u86+u94*u68
  u63=u74*u68+u75
  value=value+(c(72)*(y*(u69*(u85+u65)+u63)))
  u96=u93+u88
  u88=u99+u66
  value=value+(c(73)*((u69*(u96+u65)+u88)*z))
  u65=u90+u77
  u77=u68*(u8+u76*u68)+u84
  u99=(u71+u7)
  value=value+(c(74)*(y*(u69*(u65+u99*u69)+u77)))
  value=value+(c(75)*(z*(u70*(u65+u99*u70)+u77)))
  u60=u61*u70
  value=value+(c(76)*(y*(u70*(u96+u60)+u88)))
  value=value+(c(77)*(z*(u70*(u85+u60)+u63)))
  u60=u79*u69
  value=value+(c(78)*(x*(u69*(u64+u69*(u60+u98))+u75)))
  value=value+(c(79)*(u87*(u69*(u94+u60)+u74)*z))
  value=value+(c(80)*(x*((u92+u69*(u60+u89))*u70+u69*(u93+u9*u69)+u66))&
 &)
  value=value+(c(81)*(u81*((u8+u60)*u70+u8*u69+u73)))
  value=value+(c(82)*(x*(u70*(u93+u89*u69+(u60+u9)*u70)+u92*u69+u66)))
  u60=u79*u70
  value=value+(c(83)*(u81*(u70*(u94+u60)+u74)))
  value=value+(c(84)*(x*(u70*(u64+u70*(u60+u98))+u75)))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_opt_1_1=.false.
end function

  
  
!> Compute D X Y coulomb integral 
!!
recursive function D_X_Y_opt_1(a12,a3,r3,m3,c,lmax,value)
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12      !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3       !< exponent of the third solid harmonics 
  integer     , intent(in)  :: m3       !< m index of the solid harmonic
  real(kind=8), intent(in)  :: r3(3)    !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)     !< hermit coefficients of the two fisrt obrital produt
  integer     , intent(in)  :: lmax     !< l of maximum non 0 coeff in hermit decomposition
  real(kind=8), intent(out) :: value    !< result 
  logical                   :: D_X_Y_opt_1 !< true if the routine knows how to calculate the requested value
  
  ! init return value
  D_X_Y_opt_1 = .true.
  ! select case on m component
  select case (m3)
    case (-1)
      D_X_Y_opt_1 = D_X_Y_opt_1_m1(a12,a3,r3,c,lmax,value)
    case (0)
      D_X_Y_opt_1 = D_X_Y_opt_1_0(a12,a3,r3,c,lmax,value)
    case (1)
      D_X_Y_opt_1 = D_X_Y_opt_1_1(a12,a3,r3,c,lmax,value)
    case default
      D_X_Y_opt_1 = .false.
  end select
end function
  
end module
