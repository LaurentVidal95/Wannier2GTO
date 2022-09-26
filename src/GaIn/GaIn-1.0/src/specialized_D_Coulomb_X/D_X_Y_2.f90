module mod_D_X_Y_opt_2
  
  use mod_A_functions
  
  implicit none
  
  private
  public :: D_X_Y_opt_2
  
contains
  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_Y_opt_2_m2(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_Y_opt_2_m2
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u126
  real(kind=8) :: u177
  real(kind=8) :: u31
  real(kind=8) :: u52
  real(kind=8) :: u53
  real(kind=8) :: u54
  real(kind=8) :: u55
  real(kind=8) :: u56
  real(kind=8) :: u57
  real(kind=8) :: u58
  real(kind=8) :: u59
  real(kind=8) :: u6
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
  real(kind=8) :: u94
  real(kind=8) :: u95
  real(kind=8) :: u96
  real(kind=8) :: u97
  real(kind=8) :: u98
  real(kind=8) :: u99

  
  ! init computation flag
  D_X_Y_opt_2_m2=.true.
  
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
  u6=a3**(-2)
  u99=8.19411322944059304d-1*u6*A3_(p,pd,erfpd,exppd2)*p**3
  u126=x*y
  value=value+(c(1)*(u99*u126))
  if ( lmax .eq. 0 ) return
  u99=A4_(p,pd,erfpd,exppd2)
  u177=p**4
  u77=-8.19411322944059304d-1*u6*u99*u177
  u31=8.69422416480109528d-2*u6
  u71=exppd2*pd
  u69=-1.41796308072441282d1*u71
  u52=u31*(4.71238898038468986d1*u99+u69)*u177
  u177=x**2
  value=value+(c(2)*((u52*u177+u77)*y))
  u99=y**2
  value=value+(c(3)*(x*(u52*u99+u77)))
  u77=u126*z
  value=value+(c(4)*(u52*u77))
  if ( lmax .eq. 1 ) return
  u97=A5_(p,pd,erfpd,exppd2)
  u74=p**5
  u90=u6*u97*u74
  u78=-1.22911698441608896d1*u90
  u52=exppd2*pd2
  u60=-2.83592616144882564d1*u52
  u91=u31*(3.2986722862692829d2*u97+u60)*u74
  u98=u91*u177
  value=value+(c(5)*(x*(u98+u78)*y))
  u59=1.41796308072441282d1*exppd2
  u88=1.73884483296021906d-2*u6*(4.71238898038468986d1*u97+u59)*u74
  u74=-4.09705661472029652d0*u90
  value=value+(c(6)*((u74+u98)*u99+u74*u177+u88))
  u90=y*z
  value=value+(c(7)*((u98+u74)*u90))
  u98=u91*u99
  value=value+(c(8)*(u126*(u98+u78)))
  value=value+(c(9)*(x*(u98+u74)*z))
  u78=z**2
  value=value+(c(10)*(u126*(u91*u78+u74)))
  if ( lmax .eq. 2 ) return
  u98=A6_(p,pd,erfpd,exppd2)
  u88=p**6
  u91=u6*u98*u88
  u83=1.22911698441608896d1*u91
  u97=-2.83592616144882564d1*u71
  u68=u6*(3.2986722862692829d2*u98+u97)*u88
  u64=-5.21653449888065717d-1*u68
  u74=2.0d0*pd2
  u67=u31*(2.96880505764235461d3*u98+u97*(9.0d0+u74))*u88
  u88=u67*u177
  value=value+(c(11)*((u177*(u64+u88)+u83)*y))
  u97=-2.60826724944032859d-1*u68
  u65=-8.69422416480109528d-2*u68
  u85=u65*u177
  value=value+(c(12)*(x*((u97+u88)*u99+u85+u83)))
  value=value+(c(13)*(x*(u88+u97)*u90))
  u95=(u65+u88)
  u88=u95*u99
  value=value+(c(14)*(y*(u88+u97*u177+u83)))
  u84=4.09705661472029652d0*u91
  u92=u85+u84
  value=value+(c(15)*((u88+u92)*z))
  value=value+(c(16)*(y*(u95*u78+u92)))
  u95=u67*u99
  value=value+(c(17)*(x*(u99*(u64+u95)+u83)))
  value=value+(c(18)*(u126*(u95+u97)*z))
  value=value+(c(19)*(x*((u65+u95)*u78+u65*u99+u84)))
  value=value+(c(20)*(u77*(u67*u78+u97)))
  if ( lmax .eq. 3 ) return
  u95=A7_(p,pd,erfpd,exppd2)
  u92=p**7
  u88=u6*u95*u92
  u85=4.30190944545631135d2*u88
  u84=-5.67185232289765129d1*u52
  u67=u6*(2.96880505764235461d3*u95+u84)*u92
  u68=-8.69422416480109529d-1*u67
  u97=(1.1d1+u74)
  u64=u6*(3.26568556340659007d4*u95+u84*u97)*u92
  u65=8.69422416480109528d-2*u64
  u96=u65*u177
  value=value+(c(21)*(x*(u177*(u68+u96)+u85)*y))
  u58=-3.72609607062904084d-2*u6*(2.83592616144882564d1*exppd2+3.298672&
 &2862692829d2*u95)*u92
  u82=8.60381889091262269d1*u88
  u86=1.72076377818252454d2*u88
  u8=-5.21653449888065717d-1*u67
  u7=-8.69422416480109528d-2*u67
  u61=u7*u177
  value=value+(c(22)*((u82+u177*(u96+u8))*u99+u177*(u86+u61)+u58))
  value=value+(c(23)*((u177*(u8+u96)+u82)*u90))
  u75=2.58114566727378681d2*u88
  u81=-2.60826724944032859d-1*u67
  u89=(u81+u96)
  u70=u89*u99
  u66=u81*u177
  value=value+(c(24)*(u126*(u70+u66+u75)))
  u75=u61+u82
  value=value+(c(25)*(x*(u70+u75)*z))
  value=value+(c(26)*(u126*(u89*u78+u75)))
  value=value+(c(27)*(u99*(u86+u8*u177+(u96+u7)*u99)+u82*u177+u58))
  u70=(u7+u96)
  u89=u66+u82
  value=value+(c(28)*(y*(u70*u99+u89)*z))
  u75=2.48406404708602722d-2*u6*(9.8960168588078487d2*u95-u59)*u92
  u59=8.69422416480109528d-2*u67
  u67=3.47768966592043811d-1*u6
  u96=u67*(8.90641517292706383d3*u95+u60*(6.0d0+pd2))*u92
  u66=-8.69422416480109528d-2*u64
  u88=u66*u177
  value=value+(c(29)*(u99*(u7+u177*(u88+u96)+(u88+u59)*u99)+u177*(u7+u5&
 &9*u177)+u75))
  value=value+(c(30)*(u90*(u70*u78+u89)))
  u96=u65*u99
  value=value+(c(31)*(u126*(u99*(u68+u96)+u85)))
  value=value+(c(32)*(x*(u99*(u8+u96)+u82)*z))
  value=value+(c(33)*(u126*((u81+u96)*u78+u7*u99+u82)))
  u59=x*z
  value=value+(c(34)*(u59*((u7+u96)*u78+u81*u99+u82)))
  value=value+(c(35)*(u126*(u78*(u8+u65*u78)+u82)))
  if ( lmax .eq. 4 ) return
  u96=A8_(p,pd,erfpd,exppd2)
  u70=p**8
  u89=u6*u96*u70
  u75=-4.30190944545631135d2*u89
  u88=-5.67185232289765129d1*u71
  u66=u6*(2.96880505764235461d3*u96+u88)*u70
  u60=3.91240087416049288d0*u66
  u92=3.26568556340659007d4*u96
  u65=u6*(u92+u88*u97)*u70
  u97=-1.30413362472016429d0*u65
  u82=pd2**2
  u85=4.0d0*u82
  u64=(u74+1.1d1)
  u7=u6*(4.24539123242856709d5*u96+u88*(1.3d1*u64+u85))*u70
  u81=8.69422416480109528d-2*u7
  u8=u81*u177
  value=value+(c(36)*((u177*(u60+u177*(u8+u97))+u75)*y))
  u68=1.30413362472016429d0*u66
  u61=8.69422416480109529d-1*u66
  u98=-8.69422416480109529d-1*u65
  u94=-8.69422416480109528d-2*u65
  u63=u94*u177
  value=value+(c(37)*(x*((u68+u177*(u8+u98))*u99+u177*(u61+u63)+u75)))
  value=value+(c(38)*(x*(u177*(u98+u8)+u68)*u90))
  u79=-2.58114566727378681d2*u89
  u55=2.60826724944032859d-1*u66
  u62=1.56496034966419715d0*u66
  u87=-5.21653449888065717d-1*u65
  u9=-2.60826724944032859d-1*u65
  u73=(u55+u177*(u8+u87))
  u76=u73*u99
  u54=u9*u177
  value=value+(c(39)*(y*(u76+u177*(u62+u54)+u79)))
  u72=-8.60381889091262269d1*u89
  u86=5.21653449888065717d-1*u66
  u53=u177*(u86+u63)+u72
  value=value+(c(40)*((u76+u53)*z))
  value=value+(c(41)*(y*(u73*u78+u53)))
  u76=u87*u177
  u73=u55*u177
  value=value+(c(42)*(x*(u99*(u62+u76+(u8+u9)*u99)+u73+u79)))
  u53=7.82480174832098576d-1*u66
  u63=(u9+u8)
  u80=u54+u53
  value=value+(c(43)*(u126*(u63*u99+u80)*z))
  u54=2.08661379955226287d0*u6*(3.2986722862692829d2*u96-7.089815403622&
 &06411d0*u71)*u70
  u66=2.60826724944032859d-1*u65
  u58=(1.3d1+u74)
  u71=-8.69422416480109528d-2*u6*(3.85944657493506099d4*u96+u88*u58)*u7&
 &0
  u88=1.39107586636817525d0*u6*(u92+u69*(4.0d0*u64+u82))*u70
  u70=-8.69422416480109528d-2*u7
  u64=8.69422416480109528d-2*u65
  u69=u70*u177
  u65=u177*(u69+u88)
  value=value+(c(44)*(x*(u99*(u9+u65+(u69+u66)*u99)+u177*(u71+u64*u177)&
 &+u54)))
  value=value+(c(45)*(u77*(u63*u78+u80)))
  u63=(u8+u94)
  u80=u63*u99
  value=value+(c(46)*(y*(u99*(u61+u98*u177+u80)+u68*u177+u75)))
  u8=u86+u76
  u76=u73+u72
  value=value+(c(47)*((u99*(u8+u80)+u76)*z))
  u80=u71+u65
  u65=u177*(u9+u66*u177)+u54
  u92=(u69+u64)
  value=value+(c(48)*(y*(u99*(u80+u92*u99)+u65)))
  value=value+(c(49)*(z*(u78*(u80+u92*u78)+u65)))
  value=value+(c(50)*(y*(u78*(u8+u63*u78)+u76)))
  u80=u81*u99
  value=value+(c(51)*(x*(u99*(u60+u99*(u80+u97))+u75)))
  value=value+(c(52)*(u126*(u99*(u98+u80)+u68)*z))
  value=value+(c(53)*(x*((u55+u99*(u80+u87))*u78+u99*(u86+u94*u99)+u72)&
 &))
  value=value+(c(54)*(u77*((u9+u80)*u78+u9*u99+u53)))
  value=value+(c(55)*(x*(u78*(u86+u87*u99+(u80+u94)*u78)+u55*u99+u72)))
  value=value+(c(56)*(u77*(u78*(u98+u81*u78)+u68)))
  if ( lmax .eq. 5 ) return
  u94=A9_(p,pd,erfpd,exppd2)
  u77=p**9
  u87=u6*u94*u77
  u80=-2.71020295063747615d4*u87
  u92=-1.13437046457953026d2*u52
  u65=u6*(3.26568556340659007d4*u94+u92)*u77
  u63=9.12893537304115006d0*u65
  u66=u6*(4.24539123242856709d5*u94+u92*u58)*u77
  u68=-1.82578707460823001d0*u66
  u88=(u74+1.3d1)
  u8=u6*(6.36808684864285064d6*u94+u92*(1.5d1*u88+u85))*u77
  u60=8.69422416480109528d-2*u8
  u72=u60*u177
  value=value+(c(57)*(x*(u177*(u63+u177*(u72+u68))+u80)*y))
  u89=2.96880505764235461d3*u94
  u70=u6*(u89+5.67185232289765129d1*exppd2)*u77
  u69=1.44903736080018255d-1*u70
  u71=-3.87171850091068022d3*u87
  u76=-1.16151555027320407d4*u87
  u97=3.91240087416049288d0*u65
  u58=1.30413362472016429d0*u65
  u9=-1.30413362472016429d0*u66
  u98=-8.69422416480109528d-2*u66
  u91=u98*u177
  value=value+(c(58)*((u71+u177*(u177*(u9+u72)+u97))*u99+u177*(u76+u177&
 &*(u91+u58))+u69))
  value=value+(c(59)*((u177*(u97+u177*(u72+u9))+u71)*u90))
  u64=2.60826724944032859d0*u65
  u53=-8.69422416480109529d-1*u66
  u86=-2.60826724944032859d-1*u66
  u81=(u58+u177*(u72+u53))
  u73=u81*u99
  u55=u86*u177
  value=value+(c(60)*(u126*(u73+u177*(u64+u55)+u76)))
  u83=8.69422416480109529d-1*u65
  u95=u177*(u83+u91)+u71
  value=value+(c(61)*(x*(u73+u95)*z))
  value=value+(c(62)*(u126*(u81*u78+u95)))
  u73=8.69422416480109528d-2*u70
  u81=-4.64606220109281626d3*u87
  u95=2.60826724944032859d-1*u65
  u91=3.1299206993283943d0*u65
  u70=-5.21653449888065717d-1*u66
  u61=u70*u177
  u62=u95*u177
  value=value+(c(63)*(u99*(u81+u177*(u61+u91)+(u177*(u70+u72)+u95)*u99)&
 &+u177*(u81+u62)+u73))
  u73=-2.32303110054640813d3*u87
  u91=1.56496034966419715d0*u65
  u87=(u95+u177*(u72+u70))
  u57=u177*(u91+u55)+u73
  value=value+(c(64)*(y*(u87*u99+u57)*z))
  u55=-2.31845977728029207d-1*u6*(u89-7.08981540362206411d0*exppd2)*u77
  u89=-2.60826724944032859d-1*u65
  u56=2.60826724944032859d-1*u6*(6.23449062104894468d4*u94-2.2687409291&
 &5906051d2*u52)*u77
  u54=-2.60826724944032859d-1*u6
  u75=+u54*(8.81735102119779319d5*u94+u92*(2.7d1+4.0d0*pd2))*u77
  u7=5.21653449888065717d-1*u66
  u96=-1.73884483296021906d-1*u6*(2.93911700706593106d5*u94+u92*(9.0d0+&
 &pd2))*u77
  u79=u31*(8.9153215880999909d6*u94+u92*(2.1d1*u88+u85))*u77
  u85=-8.69422416480109528d-2*u8
  u31=8.69422416480109528d-2*u66
  u65=u85*u177
  u8=u177*(u79+u65)
  value=value+(c(65)*(u99*(u95+u177*(u8+u75)+(u177*(u7+u65)+u89)*u99)+u&
 &177*(u56+u177*(u31*u177+u96))+u55))
  value=value+(c(66)*(u90*(u87*u78+u57)))
  u52=(u72+u86)
  u57=u52*u99
  u87=u53*u177
  u81=u58*u177
  value=value+(c(67)*(u126*(u99*(u64+u87+u57)+u81+u76)))
  u64=u91+u61
  u61=u62+u73
  value=value+(c(68)*(x*(u99*(u64+u57)+u61)*z))
  u57=1.56496034966419715d0*u6*(1.48440252882117731d4*u94+u84)*u77
  u84=+u54*(4.89852834510988511d5*u94+u92*(1.5d1+u74))*u77
  u54=2.60826724944032859d-1*u66
  u6=u67*(2.12269561621428355d6*u94+u92*(5.0d0*u88+u82))*u77
  u82=u84+u177*(u65+u6)
  u77=u177*(u84+u54*u177)+u57
  u67=(u65+u54)
  value=value+(c(69)*(u126*(u99*(u82+u67*u99)+u77)))
  value=value+(c(70)*(u59*(u78*(u82+u67*u78)+u77)))
  value=value+(c(71)*(u126*(u78*(u64+u52*u78)+u61)))
  value=value+(c(72)*(u99*(u76+u97*u177+u99*((u98+u72)*u99+u9*u177+u58)&
 &)+u71*u177+u69))
  u69=(u72+u98)
  u76=u83+u87
  u52=u81+u71
  value=value+(c(73)*(y*(u99*(u76+u69*u99)+u52)*z))
  u81=u56+u177*(u7*u177+u75)
  u87=u8+u96
  u8=u177*(u95+u89*u177)+u55
  u72=(u31+u65)
  value=value+(c(74)*(u99*(u81+u99*(u72*u99+u87))+u8))
  u65=u85*u99
  value=value+(c(75)*(u90*(u78*(u84+u99*(u65+u6)+(u65+u54)*u78)+u99*(u8&
 &4+u54*u99)+u57)))
  value=value+(c(76)*(u78*(u81+u78*(u72*u78+u87))+u8))
  value=value+(c(77)*(u90*(u78*(u76+u69*u78)+u52)))
  u177=u60*u99
  value=value+(c(78)*(u126*(u99*(u63+u99*(u177+u68))+u80)))
  value=value+(c(79)*(x*(u99*(u97+u99*(u177+u9))+u71)*z))
  value=value+(c(80)*(u126*((u58+u99*(u177+u53))*u78+u99*(u83+u98*u99)+&
 &u71)))
  value=value+(c(81)*(u59*((u95+u99*(u177+u70))*u78+u99*(u91+u86*u99)+u&
 &73)))
  value=value+(c(82)*(u126*(u78*(u91+u70*u99+(u177+u86)*u78)+u95*u99+u7&
 &3)))
  value=value+(c(83)*(u59*(u78*(u83+u53*u99+(u177+u98)*u78)+u58*u99+u71&
 &)))
  value=value+(c(84)*(u126*(u78*(u97+u78*(u60*u78+u9))+u71)))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_opt_2_m2=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_Y_opt_2_m1(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_Y_opt_2_m1
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u114
  real(kind=8) :: u126
  real(kind=8) :: u52
  real(kind=8) :: u57
  real(kind=8) :: u58
  real(kind=8) :: u59
  real(kind=8) :: u6
  real(kind=8) :: u60
  real(kind=8) :: u61
  real(kind=8) :: u62
  real(kind=8) :: u63
  real(kind=8) :: u64
  real(kind=8) :: u65
  real(kind=8) :: u66
  real(kind=8) :: u67
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
  D_X_Y_opt_2_m1=.true.
  
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
  u69=a3**(-2)
  u99=8.19411322944059304d-1*u69*A3_(p,pd,erfpd,exppd2)*p**3
  u126=y*z
  value=value+(c(1)*(u99*u126))
  if ( lmax .eq. 0 ) return
  u99=A4_(p,pd,erfpd,exppd2)
  u63=8.69422416480109528d-2*u69
  u114=p**4
  u93=exppd2*pd
  u79=-1.41796308072441282d1*u93
  u52=u63*(4.71238898038468986d1*u99+u79)*u114
  u75=x*u126
  value=value+(c(2)*(u52*u75))
  u77=-8.19411322944059304d-1*u69*u99*u114
  u114=y**2
  value=value+(c(3)*((u52*u114+u77)*z))
  u99=z**2
  value=value+(c(4)*(y*(u52*u99+u77)))
  if ( lmax .eq. 1 ) return
  u97=A5_(p,pd,erfpd,exppd2)
  u84=p**5
  u90=u69*u97*u84
  u74=-4.09705661472029652d0*u90
  u6=exppd2*pd2
  u85=-2.83592616144882564d1*u6
  u65=u63*(3.2986722862692829d2*u97+u85)*u84
  u52=x**2
  value=value+(c(5)*((u65*u52+u74)*u126))
  u95=u65*u114
  value=value+(c(6)*(x*(u95+u74)*z))
  u77=x*y
  u92=u65*u99
  value=value+(c(7)*(u77*(u92+u74)))
  u78=-1.22911698441608896d1*u90
  value=value+(c(8)*(y*(u95+u78)*z))
  u58=1.41796308072441282d1*exppd2
  u98=1.73884483296021906d-2*u69*(4.71238898038468986d1*u97+u58)*u84
  value=value+(c(9)*((u74+u95)*u99+u74*u114+u98))
  value=value+(c(10)*(u126*(u92+u78)))
  if ( lmax .eq. 2 ) return
  u98=A6_(p,pd,erfpd,exppd2)
  u84=p**6
  u95=-2.83592616144882564d1*u93
  u92=u69*(3.2986722862692829d2*u98+u95)*u84
  u65=-2.60826724944032859d-1*u92
  u78=2.0d0*pd2
  u67=u63*(2.96880505764235461d3*u98+u95*(9.0d0+u78))*u84
  u95=u67*u52
  value=value+(c(11)*(x*(u95+u65)*u126))
  u91=u69*u98*u84
  u84=4.09705661472029652d0*u91
  u7=-8.69422416480109528d-2*u92
  u86=(u7+u95)
  u95=u7*u52+u84
  value=value+(c(12)*((u86*u114+u95)*z))
  value=value+(c(13)*(y*(u86*u99+u95)))
  u86=u67*u114
  value=value+(c(14)*(u77*(u86+u65)*z))
  u95=(u7+u86)*u99
  u88=u7*u114
  value=value+(c(15)*(x*(u95+u88+u84)))
  u82=u67*u99
  value=value+(c(16)*(u75*(u82+u65)))
  u83=1.22911698441608896d1*u91
  u89=-5.21653449888065717d-1*u92
  value=value+(c(17)*((u114*(u89+u86)+u83)*z))
  value=value+(c(18)*(y*((u65+u86)*u99+u88+u83)))
  value=value+(c(19)*(z*(u95+u65*u114+u83)))
  value=value+(c(20)*(y*(u99*(u89+u82)+u83)))
  if ( lmax .eq. 3 ) return
  u95=A7_(p,pd,erfpd,exppd2)
  u86=p**7
  u88=u69*u95*u86
  u82=8.60381889091262269d1*u88
  u81=-5.67185232289765129d1*u6
  u67=u69*(2.96880505764235461d3*u95+u81)*u86
  u7=-5.21653449888065717d-1*u67
  u65=(1.1d1+u78)
  u89=u69*(3.26568556340659007d4*u95+u81*u65)*u86
  u92=8.69422416480109528d-2*u89
  u90=u92*u52
  value=value+(c(21)*((u52*(u7+u90)+u82)*u126))
  u61=-2.60826724944032859d-1*u67
  u62=-8.69422416480109528d-2*u67
  u96=(u61+u90)
  u66=u62*u52+u82
  value=value+(c(22)*(x*(u96*u114+u66)*z))
  value=value+(c(23)*(u77*(u96*u99+u66)))
  u96=(u62+u90)
  u66=u61*u52+u82
  value=value+(c(24)*(y*(u96*u114+u66)*z))
  u90=2.48406404708602722d-2*u69*(9.8960168588078487d2*u95-u58)*u86
  u58=8.69422416480109528d-2*u67
  u83=3.47768966592043811d-1*u69
  u84=u83*(8.90641517292706383d3*u95+u85*(6.0d0+pd2))*u86
  u85=-8.69422416480109528d-2*u89
  u89=u85*u52
  value=value+(c(25)*(u114*(u62+u52*(u89+u84)+(u89+u58)*u114)+u52*(u62+&
 &u58*u52)+u90))
  value=value+(c(26)*(u126*(u96*u99+u66)))
  u90=u92*u114
  value=value+(c(27)*(x*(u114*(u7+u90)+u82)*z))
  u84=(u61+u90)*u99
  u96=u62*u114
  value=value+(c(28)*(u77*(u84+u96+u82)))
  u58=x*z
  u66=u61*u114
  value=value+(c(29)*(u58*((u62+u90)*u99+u66+u82)))
  u89=u92*u99
  value=value+(c(30)*(u77*(u99*(u7+u89)+u82)))
  u85=4.30190944545631135d2*u88
  u80=-8.69422416480109529d-1*u67
  value=value+(c(31)*(y*(u114*(u80+u90)+u85)*z))
  u59=-3.72609607062904084d-2*u69*(2.83592616144882564d1*exppd2+3.29867&
 &22862692829d2*u95)*u86
  u86=1.72076377818252454d2*u88
  value=value+(c(32)*((u82+u114*(u90+u7))*u99+u114*(u86+u96)+u59))
  u96=2.58114566727378681d2*u88
  value=value+(c(33)*(u126*(u84+u66+u96)))
  value=value+(c(34)*(u99*(u86+u7*u114+(u90+u62)*u99)+u82*u114+u59))
  value=value+(c(35)*(u126*(u99*(u80+u89)+u85)))
  if ( lmax .eq. 4 ) return
  u96=A8_(p,pd,erfpd,exppd2)
  u84=p**8
  u90=-5.67185232289765129d1*u93
  u66=u69*(2.96880505764235461d3*u96+u90)*u84
  u67=1.30413362472016429d0*u66
  u92=3.26568556340659007d4*u96
  u62=u69*(u92+u90*u65)*u84
  u65=-8.69422416480109529d-1*u62
  u59=pd2**2
  u85=4.0d0*u59
  u61=(u78+1.1d1)
  u7=u69*(4.24539123242856709d5*u96+u90*(1.3d1*u61+u85))*u84
  u80=8.69422416480109528d-2*u7
  u88=u80*u52
  value=value+(c(36)*(x*(u52*(u65+u88)+u67)*u126))
  u89=u69*u96*u84
  u72=-8.60381889091262269d1*u89
  u97=2.60826724944032859d-1*u66
  u82=5.21653449888065717d-1*u66
  u94=-5.21653449888065717d-1*u62
  u87=-8.69422416480109528d-2*u62
  u91=(u97+u52*(u88+u94))
  u9=u52*(u82+u87*u52)+u72
  value=value+(c(37)*((u91*u114+u9)*z))
  value=value+(c(38)*(y*(u91*u99+u9)))
  u91=7.82480174832098576d-1*u66
  u9=-2.60826724944032859d-1*u62
  u70=(u9+u88)
  u8=u9*u52+u91
  value=value+(c(39)*(u77*(u70*u114+u8)*z))
  u98=2.08661379955226287d0*u69*(3.2986722862692829d2*u96-7.08981540362&
 &206411d0*u93)*u84
  u93=2.60826724944032859d-1*u62
  u95=(1.3d1+u78)
  u60=-8.69422416480109528d-2*u69*(3.85944657493506099d4*u96+u90*u95)*u&
 &84
  u90=1.39107586636817525d0*u69*(u92+u79*(4.0d0*u61+u59))*u84
  u84=-8.69422416480109528d-2*u7
  u61=8.69422416480109528d-2*u62
  u96=u84*u52
  u79=u52*(u96+u90)
  value=value+(c(40)*(x*(u114*(u9+u79+(u96+u93)*u114)+u52*(u60+u61*u52)&
 &+u98)))
  value=value+(c(41)*(u75*(u70*u99+u8)))
  u70=(u88+u87)
  u8=u82+u94*u52
  u88=u97*u52+u72
  value=value+(c(42)*((u114*(u8+u70*u114)+u88)*z))
  u71=u60+u79
  u79=u52*(u9+u93*u52)+u98
  u73=(u96+u61)
  value=value+(c(43)*(y*(u114*(u71+u73*u114)+u79)))
  value=value+(c(44)*(z*(u99*(u71+u73*u99)+u79)))
  value=value+(c(45)*(y*(u99*(u8+u70*u99)+u88)))
  u84=u80*u114
  value=value+(c(46)*(u77*(u114*(u65+u84)+u67)*z))
  u71=(u97+u114*(u84+u94))*u99
  u79=u87*u114
  value=value+(c(47)*(x*(u71+u114*(u82+u79)+u72)))
  u70=u9*u114
  value=value+(c(48)*(u75*((u9+u84)*u99+u70+u91)))
  u73=(u84+u87)*u99
  u8=u94*u114
  u88=u97*u114
  value=value+(c(49)*(x*(u99*(u82+u8+u73)+u88+u72)))
  u96=u80*u99
  value=value+(c(50)*(u75*(u99*(u65+u96)+u67)))
  u75=-4.30190944545631135d2*u89
  u60=3.91240087416049288d0*u66
  u98=-1.30413362472016429d0*u62
  value=value+(c(51)*((u114*(u60+u114*(u84+u98))+u75)*z))
  u61=8.69422416480109529d-1*u66
  value=value+(c(52)*(y*((u67+u114*(u84+u65))*u99+u114*(u61+u79)+u75)))
  u79=-2.58114566727378681d2*u89
  u62=1.56496034966419715d0*u66
  value=value+(c(53)*(z*(u71+u114*(u62+u70)+u79)))
  value=value+(c(54)*(y*(u99*(u62+u8+(u84+u9)*u99)+u88+u79)))
  value=value+(c(55)*(z*(u99*(u61+u65*u114+u73)+u67*u114+u75)))
  value=value+(c(56)*(y*(u99*(u60+u99*(u96+u98))+u75)))
  if ( lmax .eq. 5 ) return
  u94=A9_(p,pd,erfpd,exppd2)
  u61=p**9
  u87=u69*u94*u61
  u71=-3.87171850091068022d3*u87
  u79=-1.13437046457953026d2*u6
  u65=u69*(3.26568556340659007d4*u94+u79)*u61
  u84=3.91240087416049288d0*u65
  u98=u69*(4.24539123242856709d5*u94+u79*u95)*u61
  u91=-1.30413362472016429d0*u98
  u66=(u78+1.3d1)
  u8=u69*(6.36808684864285064d6*u94+u79*(1.5d1*u66+u85))*u61
  u62=8.69422416480109528d-2*u8
  u82=u62*u52
  value=value+(c(57)*((u52*(u84+u52*(u82+u91))+u71)*u126))
  u96=1.30413362472016429d0*u65
  u88=8.69422416480109529d-1*u65
  u9=-8.69422416480109529d-1*u98
  u67=-8.69422416480109528d-2*u98
  u73=(u96+u52*(u82+u9))
  u92=u52*(u88+u67*u52)+u71
  value=value+(c(58)*(x*(u73*u114+u92)*z))
  value=value+(c(59)*(u77*(u73*u99+u92)))
  u73=-2.32303110054640813d3*u87
  u92=2.60826724944032859d-1*u65
  u95=1.56496034966419715d0*u65
  u60=-5.21653449888065717d-1*u98
  u75=-2.60826724944032859d-1*u98
  u80=(u92+u52*(u82+u60))
  u70=u52*(u95+u75*u52)+u73
  value=value+(c(60)*(y*(u80*u114+u70)*z))
  u93=2.96880505764235461d3*u94
  u97=-2.31845977728029207d-1*u69*(u93-7.08981540362206411d0*exppd2)*u6&
 &1
  u7=-2.60826724944032859d-1*u65
  u90=2.60826724944032859d-1*u69*(6.23449062104894468d4*u94-2.268740929&
 &15906051d2*u6)*u61
  u6=-2.60826724944032859d-1*u69
  u76=+u6*(8.81735102119779319d5*u94+u79*(2.7d1+4.0d0*pd2))*u61
  u89=5.21653449888065717d-1*u98
  u86=-1.73884483296021906d-1*u69*(2.93911700706593106d5*u94+u79*(9.0d0&
 &+pd2))*u61
  u74=u63*(8.9153215880999909d6*u94+u79*(2.1d1*u66+u85))*u61
  u85=-8.69422416480109528d-2*u8
  u72=8.69422416480109528d-2*u98
  u63=u85*u52
  u64=u52*(u74+u63)
  value=value+(c(61)*(u114*(u92+u52*(u64+u76)+(u52*(u89+u63)+u7)*u114)+&
 &u52*(u90+u52*(u72*u52+u86))+u97))
  value=value+(c(62)*(u126*(u80*u99+u70)))
  u80=(u82+u75)
  u70=u95+u60*u52
  u8=u92*u52+u73
  value=value+(c(63)*(x*(u114*(u70+u80*u114)+u8)*z))
  u57=1.56496034966419715d0*u69*(1.48440252882117731d4*u94+u81)*u61
  u81=+u6*(4.89852834510988511d5*u94+u79*(1.5d1+u78))*u61
  u78=2.60826724944032859d-1*u98
  u6=u83*(2.12269561621428355d6*u94+u79*(5.0d0*u66+u59))*u61
  u66=u81+u52*(u63+u6)
  u59=u52*(u81+u78*u52)+u57
  u79=(u63+u78)
  value=value+(c(64)*(u77*(u114*(u66+u79*u114)+u59)))
  value=value+(c(65)*(u58*(u99*(u66+u79*u99)+u59)))
  value=value+(c(66)*(u77*(u99*(u70+u80*u99)+u8)))
  u79=(u82+u67)
  u8=u88+u9*u52
  u80=u96*u52+u71
  value=value+(c(67)*(y*(u114*(u8+u79*u114)+u80)*z))
  u59=u90+u52*(u89*u52+u76)
  u70=u64+u86
  u64=u52*(u92+u7*u52)+u97
  u82=(u72+u63)
  value=value+(c(68)*(u114*(u59+u114*(u82*u114+u70))+u64))
  u63=u85*u114
  value=value+(c(69)*(u126*(u99*(u81+u114*(u63+u6)+(u63+u78)*u99)+u114*&
 &(u81+u78*u114)+u57)))
  value=value+(c(70)*(u99*(u59+u99*(u82*u99+u70))+u64))
  value=value+(c(71)*(u126*(u99*(u8+u79*u99)+u80)))
  u8=u62*u114
  value=value+(c(72)*(x*(u114*(u84+u114*(u8+u91))+u71)*z))
  u74=(u96+u114*(u8+u9))*u99
  u79=u67*u114
  value=value+(c(73)*(u77*(u74+u114*(u88+u79)+u71)))
  u72=u75*u114
  value=value+(c(74)*(u58*((u92+u114*(u8+u60))*u99+u114*(u95+u72)+u73))&
 &)
  u52=(u8+u75)*u99
  u59=u60*u114
  u7=u92*u114
  value=value+(c(75)*(u77*(u99*(u95+u59+u52)+u7+u73)))
  u95=u9*u114
  u73=u96*u114
  value=value+(c(76)*(u58*(u99*(u88+u95+(u8+u67)*u99)+u73+u71)))
  u58=u62*u99
  value=value+(c(77)*(u77*(u99*(u84+u99*(u58+u91))+u71)))
  u80=-2.71020295063747615d4*u87
  u63=9.12893537304115006d0*u65
  u77=-1.82578707460823001d0*u98
  value=value+(c(78)*(y*(u114*(u63+u114*(u8+u77))+u80)*z))
  u70=u69*(u93+5.67185232289765129d1*exppd2)*u61
  u69=1.44903736080018255d-1*u70
  u76=-1.16151555027320407d4*u87
  value=value+(c(79)*((u71+u114*(u114*(u91+u8)+u84))*u99+u114*(u76+u114&
 &*(u79+u96))+u69))
  u64=2.60826724944032859d0*u65
  value=value+(c(80)*(u126*(u74+u114*(u64+u72)+u76)))
  u72=8.69422416480109528d-2*u70
  u81=-4.64606220109281626d3*u87
  u87=3.1299206993283943d0*u65
  value=value+(c(81)*(u99*(u81+u114*(u59+u87)+(u114*(u60+u8)+u92)*u99)+&
 &u114*(u81+u7)+u72))
  value=value+(c(82)*(u126*(u99*(u64+u95+u52)+u73+u76)))
  value=value+(c(83)*(u99*(u76+u84*u114+u99*((u67+u8)*u99+u91*u114+u96)&
 &)+u71*u114+u69))
  value=value+(c(84)*(u126*(u99*(u63+u99*(u58+u77))+u80)))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_opt_2_m1=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_Y_opt_2_0(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_Y_opt_2_0
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u237
  real(kind=8) :: u238
  real(kind=8) :: u239
  real(kind=8) :: u240
  real(kind=8) :: u339
  real(kind=8) :: u34
  real(kind=8) :: u340
  real(kind=8) :: u341
  real(kind=8) :: u342
  real(kind=8) :: u343
  real(kind=8) :: u344
  real(kind=8) :: u345
  real(kind=8) :: u346
  real(kind=8) :: u347
  real(kind=8) :: u348
  real(kind=8) :: u349
  real(kind=8) :: u35
  real(kind=8) :: u350
  real(kind=8) :: u351
  real(kind=8) :: u352
  real(kind=8) :: u353
  real(kind=8) :: u354
  real(kind=8) :: u355
  real(kind=8) :: u356
  real(kind=8) :: u357
  real(kind=8) :: u358
  real(kind=8) :: u359
  real(kind=8) :: u36
  real(kind=8) :: u360
  real(kind=8) :: u361
  real(kind=8) :: u362
  real(kind=8) :: u363
  real(kind=8) :: u364
  real(kind=8) :: u365
  real(kind=8) :: u366
  real(kind=8) :: u367
  real(kind=8) :: u368
  real(kind=8) :: u369
  real(kind=8) :: u37
  real(kind=8) :: u370
  real(kind=8) :: u371
  real(kind=8) :: u372
  real(kind=8) :: u373
  real(kind=8) :: u374
  real(kind=8) :: u375
  real(kind=8) :: u376
  real(kind=8) :: u377
  real(kind=8) :: u378
  real(kind=8) :: u379
  real(kind=8) :: u38
  real(kind=8) :: u380
  real(kind=8) :: u39
  real(kind=8) :: u4
  real(kind=8) :: u40
  real(kind=8) :: u41
  real(kind=8) :: u42
  real(kind=8) :: u43
  real(kind=8) :: u44
  real(kind=8) :: u45
  real(kind=8) :: u46
  real(kind=8) :: u47
  real(kind=8) :: u48
  real(kind=8) :: u49
  real(kind=8) :: u5
  real(kind=8) :: u50
  real(kind=8) :: u51
  real(kind=8) :: u52
  real(kind=8) :: u53
  real(kind=8) :: u54
  real(kind=8) :: u55
  real(kind=8) :: u56
  real(kind=8) :: u57
  real(kind=8) :: u58
  real(kind=8) :: u59
  real(kind=8) :: u6
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
  D_X_Y_opt_2_0=.true.
  
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
  u240=a3**(-2)
  u237=u240*A3_(p,pd,erfpd,exppd2)*p**3
  u62=4.73087347878780009d-1*u237
  u56=-2.36543673939390005d-1*u237
  u237=x**2
  u238=y**2
  u239=z**2
  value=value+(c(1)*(u62*u239+u56*u238+u56*u237))
  if ( lmax .eq. 0 ) return
  u62=A4_(p,pd,erfpd,exppd2)
  u56=p**4
  u66=u240*u62*u56
  u86=4.73087347878780009d-1*u66
  u83=-1.41796308072441282d1*exppd2
  u340=+u83*pd
  u57=u240*(4.71238898038468986d1*u62+u340)*u56
  u62=5.0196126619428616d-2*u57
  u56=-2.5098063309714308d-2*u57
  u57=u62*u239+u56*u238+u56*u237
  u78=(u57+u86)
  value=value+(c(2)*(x*u78))
  value=value+(c(3)*(y*u78))
  u86=-9.46174695757560018d-1*u66
  value=value+(c(4)*(z*(u57+u86)))
  if ( lmax .eq. 1 ) return
  u86=A5_(p,pd,erfpd,exppd2)
  u56=4.71238898038468986d1*u86
  u62=p**5
  u66=-1.00392253238857232d-2*u240*(-u83+u56)*u62
  u78=u240*u86*u62
  u57=-2.36543673939390005d0*u78
  u71=1.18271836969695002d0*u78
  u346=5.91359184848475012d0*u78
  u8=exppd2*pd2
  u362=-2.83592616144882564d1*u8
  u354=u240*(3.2986722862692829d2*u86+u362)*u62
  u86=5.0196126619428616d-2*u354
  u95=-2.5098063309714308d-2*u354
  u354=u95*u237
  value=value+(c(5)*((u57+u86*u237)*u239+(u71+u354)*u238+u237*(u346+u35&
 &4)+u66))
  u380=4.73087347878780009d0*u78
  u54=x*y
  u352=u86*u239
  u353=u95*u238
  u4=u352+u353+u354
  value=value+(c(6)*(u54*(u4+u380)))
  u380=x*z
  u376=(u4+u57)
  value=value+(c(7)*(u380*u376))
  u4=u354+u353
  u354=u71*u237
  value=value+(c(8)*((u57+u86*u238)*u239+u238*(u346+u4)+u354+u66))
  u353=y*z
  value=value+(c(9)*(u353*u376))
  u57=2.00784506477714464d-2*u240*(u56-u83)*u62
  u346=-1.18271836969695002d1*u78
  value=value+(c(10)*(u239*(u346+u4+u352)+u71*u238+u354+u57))
  if ( lmax .eq. 2 ) return
  u57=A6_(p,pd,erfpd,exppd2)
  u346=p**6
  u86=u240*u57*u346
  u71=-1.41926204363634003d1*u86
  u78=-2.83592616144882564d1*exppd2
  u56=+u78*pd
  u62=u240*(3.2986722862692829d2*u57+u56)*u346
  u66=-1.50588379858285848d-1*u62
  u376=7.52941899291429239d-2*u62
  u4=2.25882569787428772d-1*u62
  u352=2.0d0*pd2
  u354=(9.0d0+u352)
  u51=exppd2*pd*u354
  u90=u240*(2.96880505764235461d3*u57-2.83592616144882564d1*u51)*u346
  u57=5.0196126619428616d-2*u90
  u346=-2.5098063309714308d-2*u90
  u45=u346*u237
  u355=u57*u237
  value=value+(c(11)*(x*((u66+u355)*u239+(u376+u45)*u238+u237*(u4+u45)+&
 &u71)))
  u359=-4.73087347878780009d0*u86
  u36=-5.0196126619428616d-2*u62
  u37=2.5098063309714308d-2*u62
  u49=1.75686443168000156d-1*u62
  u43=(u37+u45)
  u356=(u36+u355)*u239+u43*u238
  value=value+(c(12)*(y*(u356+u237*(u49+u45)+u359)))
  u85=2.36543673939390005d0*u86
  value=value+(c(13)*(z*(u356+u237*u43+u85)))
  u356=u57*u238
  u43=(u36+u356)*u239
  u357=u346*u238
  u87=u45+u357
  u358=u37*u237
  value=value+(c(14)*(x*(u43+u238*(u49+u87)+u358+u359)))
  u49=u54*z
  u359=u57*u239
  value=value+(c(15)*(u49*(u359+u357+u45)))
  u45=-2.00784506477714464d-1*u62
  u372=u87+u359
  u55=u358+u85
  u89=(u239*(u45+u372)+u37*u238+u55)
  value=value+(c(16)*(x*u89))
  u360=u376*u237
  value=value+(c(17)*(y*((u66+u356)*u239+u238*(u4+u87)+u360+u71)))
  value=value+(c(18)*(z*(u43+u238*(u37+u87)+u55)))
  value=value+(c(19)*(y*u89))
  u36=2.83852408727268006d1*u86
  u66=-4.51765139574857544d-1*u62
  value=value+(c(20)*(z*(u239*(u66+u372)+u376*u238+u360+u36)))
  if ( lmax .eq. 3 ) return
  u36=A7_(p,pd,erfpd,exppd2)
  u66=3.2986722862692829d2*u36
  u45=p**7
  u37=4.30252513880816708d-2*u240*(u66-u78)*u45
  u376=u240*u36*u45
  u4=4.9674171527271901d1*u376
  u62=-2.48370857636359505d1*u376
  u71=-3.22882114927267356d2*u376
  u85=2.96880505764235461d3*u36
  u360=-5.67185232289765129d1*u8
  u86=u240*(u85+u360)*u45
  u346=-3.01176759716571696d-1*u86
  u89=1.50588379858285848d-1*u86
  u372=3.51372886336000312d-1*u86
  u43=(1.1d1+u352)
  u87=u240*(3.26568556340659007d4*u36+u360*u43)*u45
  u55=5.0196126619428616d-2*u87
  u69=-2.5098063309714308d-2*u87
  u364=u69*u237
  u339=u55*u237
  value=value+(c(21)*((u4+u237*(u339+u346))*u239+(u62+u237*(u364+u89))*&
 &u238+u237*(u71+u237*(u364+u372))+u37))
  u9=-1.49022514581815703d2*u376
  u58=+u9
  u99=-1.50588379858285848d-1*u86
  u367=7.52941899291429239d-2*u86
  u363=2.76078696406857388d-1*u86
  u97=(u367+u364)
  u356=(u99+u339)*u239+u97*u238
  value=value+(c(22)*(u54*(u356+u237*(u363+u364)+u58)))
  u38=1.2549031654857154d-1*u86
  value=value+(c(23)*(u380*(u356+u237*(u38+u364))))
  u356=5.67185232289765129d1*exppd2
  u369=7.17087523134694514d-3*u240*(u85+u356)*u45
  u64=8.90641517292706383d3*u36
  u358=-2.5098063309714308d-2*u240
  u359=-1.13437046457953026d2*u8
  u77=+u358*(u64+u359)*u45
  u57=u8*(9.0d0+pd2)
  u76=1.00392253238857232d-1*u240*(2.67192455187811915d4*u36-5.67185232&
 &289765129d1*u57)*u45
  u36=-7.52941899291429239d-2*u87
  u73=u367*u237
  u361=u36*u237
  value=value+(c(24)*(u238*(u77+u237*(u361+u76)+(u361+u367)*u238)+u237*&
 &(u77+u73)+u369))
  u369=-5.0196126619428616d-2*u86
  u76=2.5098063309714308d-2*u86
  value=value+(c(25)*(u353*((u369+u339)*u239+(u76+u364)*u238+u237*u97))&
 &)
  u77=-3.58543761567347257d-3*u240*(u356+u85)*u45
  u85=7.45112572909078515d1*u376
  u97=-2.5098063309714308d-2*u86
  u36=-1.00392253238857232d-1*u240
  u70=+u36*(u64+u362*(6.0d0+pd2))*u45
  u64=2.5098063309714308d-2*u87
  u362=u64*u237
  u72=u362+u70
  u342=u76*u237
  value=value+(c(26)*(u239*(u85+u237*(u364+u99)+(u339+u369)*u239)+u238*&
 &(u76+u237*(u72)+(u362+u97)*u238)+u342+u77))
  u339=u55*u238
  u42=(u99+u339)*u239
  u343=u69*u238
  u348=u364+u343
  value=value+(c(27)*(u54*(u42+u238*(u363+u348)+u73+u58)))
  value=value+(c(28)*(u380*((u369+u339)*u239+u238*(u367+u348)+u342)))
  u363=u55*u239
  u58=u348+u363
  value=value+(c(29)*(u54*(u239*(u99+u58)+u76*u238+u342)))
  u347=-u9
  u9=-4.01569012955428928d-1*u86
  u96=(u239*(u9+u58)+u367*u238+u73+u347)
  value=value+(c(30)*(u380*u96))
  u58=u343+u364
  u364=u89*u237
  u365=u62*u237
  value=value+(c(31)*((u4+u238*(u339+u346))*u239+u238*(u71+u364+u238*(u&
 &58+u372))+u365+u37))
  value=value+(c(32)*(u353*(u42+u238*(u38+u348)+u73)))
  value=value+(c(33)*(u239*(u85+u238*(u343+u99)+(u339+u369)*u239)+u238*&
 &(u76+u237*(u64*u238+u72))+u237*(u76+u97*u237)+u77))
  value=value+(c(34)*(u353*u96))
  u369=-8.60505027761633417d-2*u240*(-u78+u66)*u45
  u346=6.45764229854534713d2*u376
  u97=-7.02745772672000623d-1*u86
  value=value+(c(35)*(u239*(u346+u364+u89*u238+u239*(u363+u58+u97))+u62&
 &*u238+u365+u369))
  if ( lmax .eq. 4 ) return
  u369=A8_(p,pd,erfpd,exppd2)
  u346=p**8
  u97=u240*u369*u346
  u9=7.45112572909078515d2*u97
  u367=-u356*pd
  u76=u240*(2.96880505764235461d3*u369+u367)*u346
  u89=7.5294189929142924d-1*u76
  u372=-3.7647094964571462d-1*u76
  u362=-2.38431601442285926d0*u76
  u77=3.26568556340659007d4*u369
  u86=u240*(u77+u367*u43)*u346
  u37=-5.0196126619428616d-1*u86
  u62=2.5098063309714308d-1*u86
  u71=5.0196126619428616d-1*u86
  u339=pd2**2
  u363=4.0d0*u339
  u347=(u352+1.1d1)
  u376=u240*(4.24539123242856709d5*u369+u367*(1.3d1*u347+u363))*u346
  u95=5.0196126619428616d-2*u376
  u85=-2.5098063309714308d-2*u376
  u45=u85*u237
  u4=(u372+u237*(u45+u62))
  u69=u95*u237
  value=value+(c(36)*(x*((u89+u237*(u69+u37))*u239+u4*u238+u237*(u362+u&
 &237*(u45+u71))+u9)))
  u43=1.49022514581815703d2*u97
  u78=1.50588379858285848d-1*u76
  u96=-7.52941899291429239d-2*u76
  u42=-1.28000122879542971d0*u76
  u348=-3.01176759716571696d-1*u86
  u58=1.50588379858285848d-1*u86
  u72=4.01569012955428928d-1*u86
  u66=(u78+u237*(u69+u348))*u239+(u96+u237*(u45+u58))*u238
  value=value+(c(37)*(y*(u66+u237*(u42+u237*(u45+u72))+u43)))
  value=value+(c(38)*(z*(u66+u237*u4)))
  u66=6.02353519433143391d-1*u240*(9.8960168588078487d2*u369+u340)*u346
  u4=9.79705669021977021d4*u369
  u365=-7.52941899291429239d-2*u240
  u364=4.0d0*pd2
  u38=+u365*(u4+u367*(3.3d1+u364))*u346
  u73=2.25882569787428772d-1*u86
  u39=+u358*(1.1578339724805183d5*u369+u367*(3.9d1+u364))*u346
  u91=4.01569012955428928d-1*u240*(u4+u56*(6.0d0*u347+u339))*u346
  u4=-7.52941899291429239d-2*u376
  u56=7.52941899291429239d-2*u86
  u344=u4*u237
  u41=u237*(u344+u91)
  u345=u56*u237
  value=value+(c(39)*(x*(u238*(u38+u41+(u344+u73)*u238)+u237*(u39+u345)&
 &+u66)))
  u65=-1.50588379858285848d-1*u76
  u377=-1.50588379858285848d-1*u86
  u67=1.75686443168000156d-1*u86
  value=value+(c(40)*(u49*((u377+u69)*u239+(u56+u45)*u238+u237*(u67+u45&
 &)+u65)))
  u64=5.27059329504000468d-1*u76
  u60=-7.52941899291429239d-2*u86
  u55=2.5098063309714308d-2*u240
  u74=u55*(2.67192455187811915d4*u369-5.67185232289765129d1*u51)*u346
  u51=-5.0196126619428616d-2*u86
  u6=u240*(u77+u340*(4.0d0*u347+u339))*u346
  u77=-4.01569012955428928d-1*u6
  u90=2.5098063309714308d-2*u376
  u340=u90*u237
  u98=u340+u77
  u50=u237*(u98)
  u368=u56+u50
  u34=u340+u60
  value=value+(c(41)*(x*(u239*(u64+u237*(u45+u51)+(u69+u377)*u239)+u238&
 &*(u368+(u34)*u238)+u74*u237+u96)))
  value=value+(c(42)*(y*(u238*(u39+u41+(u344+u56)*u238)+u237*(u38+u73*u&
 &237)+u66)))
  u66=(1.3d1+u352)
  u91=u240*(3.85944657493506099d4*u369+u367*u66)*u346
  u41=-5.0196126619428616d-2*u91
  u75=5.0196126619428616d-2*u86
  u40=-2.25882569787428772d-1*u76
  u59=2.5098063309714308d-2*u86
  u52=+u365*(7.42201264410588653d4*u369+u367*(2.5d1+u364))*u346
  u92=8.03138025910857856d-1*u6
  u6=-5.0196126619428616d-2*u376
  u44=3.01176759716571696d-1*u86
  u366=u6*u237
  value=value+(c(43)*(z*(u239*(u41+u237*(u366+u92)+(u366+u75)*u239)+u23&
 &8*(u40+u237*(u45+u44)+(u45+u59)*u238)+u237*(u52+u67*u237)+u78)))
  u41=2.25882569787428772d-1*u76
  u75=2.5098063309714308d-2*u91
  u40=-2.5098063309714308d-2*u86
  value=value+(c(44)*(y*(u239*(u41+u237*(u45+u377)+(u69+u51)*u239)+u238&
 &*(u75+u50+(u340+u40)*u238)+u237*(u56+u51*u237)+u96)))
  u59=+u365*(4.94800842940392435d3*u369+u367)*u346
  u52=4.45320758646353191d4*u369
  u361=5.0196126619428616d-2*u240
  u6=(1.5d1+pd2)
  u92=u361*(u52+u367*u6)*u346
  u44=7.52941899291429239d-2*u240
  u366=(1.5d1+u352)
  u73=u44*(u52+u367*u366)*u346
  u52=-5.0196126619428616d-2*u240
  u91=2.0d0*u339
  u69=+u52*(4.89852834510988511d5*u369+u367*(1.5d1*u347+u91))*u346
  u344=7.52941899291429239d-2*u376
  value=value+(c(45)*(z*(u239*(u92+u69*u237+(u344*u237+u60)*u239)+u73*u&
 &237+u59)))
  u73=u95*u238
  u69=(u78+u238*(u73+u348))*u239
  u344=u85*u238
  u376=u344+u45
  u346=u58*u237
  u347=u96*u237
  value=value+(c(46)*(x*(u69+u238*(u42+u346+u238*(u376+u72))+u347+u43))&
 &)
  u42=u45+u344
  value=value+(c(47)*(u49*((u377+u73)*u239+u238*(u67+u42)+u345+u65)))
  value=value+(c(48)*(x*(u239*(u41+u238*(u344+u377)+(u73+u51)*u239)+u23&
 &8*(u368+(u340+u51)*u238)+u237*(u75+u40*u237)+u96)))
  u65=3.01176759716571696d-1*u76
  u41=-3.51372886336000312d-1*u86
  u367=u95*u239
  u368=u56*u238
  value=value+(c(49)*(u49*(u239*(u41+u42+u367)+u368+u345+u65)))
  u49=-u43
  u43=1.65647217844114433d0*u76
  u42=-6.52549646052572008d-1*u86
  u45=u367+u376
  u345=(u239*(u43+u346+u58*u238+u239*(u45+u42))+u96*u238+u347+u49)
  value=value+(c(50)*(x*u345))
  u369=u62*u237
  u370=u372*u237
  value=value+(c(51)*(y*((u89+u238*(u73+u37))*u239+u238*(u362+u369+u238&
 &*(u376+u71))+u370+u9)))
  value=value+(c(52)*(z*(u69+u238*(u372+u346+u238*(u376+u62))+u347)))
  u367=u237*(u56+u60*u237)
  value=value+(c(53)*(y*(u239*(u64+u238*(u344+u51)+(u73+u377)*u239)+u23&
 &8*(u74+u237*(u90*u238+u98))+u367+u96)))
  value=value+(c(54)*(z*(u239*(u92+u50+u238*(u344+u41)+(u73+u34)*u239)+&
 &u238*(u65+u368)+u367+u59)))
  value=value+(c(55)*(y*u345))
  u367=-1.49022514581815703d3*u97
  u347=4.76863202884571852d0*u76
  u345=-1.00392253238857232d0*u86
  value=value+(c(56)*(z*(u239*(u347+u369+u62*u238+u239*(u45+u345))+u372&
 &*u238+u370+u367)))
  if ( lmax .eq. 5 ) return
  u367=A9_(p,pd,erfpd,exppd2)
  u347=2.96880505764235461d3*u367
  u345=p**9
  u42=-2.5098063309714308d-1*u240*(u356+u347)*u345
  u56=u240*u367*u345
  u37=-2.23533771872723555d3*u56
  u96=+u37
  u67=1.11766885936361777d3*u56
  u73=2.79417214840904443d4*u56
  u71=u240*(3.26568556340659007d4*u367+u359)*u345
  u92=2.25882569787428772d0*u71
  u340=-1.12941284893714386d0*u71
  u369=-4.89412234539429006d0*u71
  u85=u240*(4.24539123242856709d5*u367+u359*u66)*u345
  u43=-7.5294189929142924d-1*u85
  u65=3.7647094964571462d-1*u85
  u76=6.77647709362286316d-1*u85
  u49=6.36808684864285064d6*u367
  u344=(u352+1.3d1)
  u89=1.5d1*u344
  u75=u240*(u49+u359*(u89+u363))*u345
  u9=5.0196126619428616d-2*u75
  u39=-2.5098063309714308d-2*u75
  u64=u39*u237
  u95=u9*u237
  value=value+(c(57)*((u96+u237*(u237*(u43+u95)+u92))*u239+(u67+u237*(u&
 &237*(u65+u64)+u340))*u238+u237*(u73+u237*(u237*(u76+u64)+u369))+u42))
  u362=8.94135087490894218d3*u56
  u62=7.5294189929142924d-1*u71
  u97=-3.7647094964571462d-1*u71
  u370=-2.88627728061714542d0*u71
  u90=-5.0196126619428616d-1*u85
  u59=2.5098063309714308d-1*u85
  u77=5.52157392813714776d-1*u85
  u372=(u62+u237*(u95+u90))*u239+(u97+u237*(u64+u59))*u238
  value=value+(c(58)*(u54*(u372+u237*(u370+u237*(u64+u77))+u362)))
  u357=-u37
  u37=-1.38039348203428694d0*u71
  u66=4.01569012955428928d-1*u85
  u377=u237*(u64+u66)
  value=value+(c(59)*(u380*(u372+u237*(u37+u377)+u357)))
  u372=-2.00784506477714464d-1*u240*(-u83+u347)*u345
  u346=-2.26874092915906051d2*u8
  u83=u44*(9.79705669021977021d4*u367+u346)*u345
  u342=-2.25882569787428772d-1*u71
  u51=u44*(1.8703471863146834d5*u367-4.53748185831812103d2*u8)*u345
  u44=+u365*(2.64520530635933796d6*u367+u359*(8.1d1+8.0d0*pd2))*u345
  u78=4.51765139574857544d-1*u85
  u69=8.81735102119779319d5*u367
  u45=+u52*(u69+u359*(2.7d1+u352))*u345
  u52=u55*(2.67459647642999727d7*u367+u359*(6.3d1*u344+8.0d0*u339))*u34&
 &5
  u50=-7.52941899291429239d-2*u75
  u60=7.52941899291429239d-2*u85
  u355=u50*u237
  u376=u237*(u52+u355)
  u348=u60*u237
  value=value+(c(60)*(u238*(u83+u237*(u376+u44)+(u237*(u78+u355)+u342)*&
 &u238)+u237*(u51+u237*(u348+u45))+u372))
  u368=4.47067543745447109d2*u56
  u38=1.50588379858285848d-1*u71
  u98=-7.52941899291429239d-2*u71
  u74=-6.77647709362286316d-1*u71
  u34=-3.01176759716571696d-1*u85
  u61=1.50588379858285848d-1*u85
  u68=3.01176759716571696d-1*u85
  value=value+(c(61)*(u353*((u38+u237*(u95+u34))*u239+(u98+u237*(u64+u6&
 &1))*u238+u237*(u74+u237*(u64+u68))+u368)))
  u41=2.23533771872723554d2*u56
  u40=-1.56473640310906488d3*u56
  u86=7.52941899291429239d-2*u71
  u87=-1.50588379858285848d-1*u240
  u4=+u87*(2.67192455187811915d4*u367+u359)*u345
  u58=6.77647709362286316d-1*u71
  u88=u240*(u69+u359*(2.7d1+u364))*u345
  u69=7.52941899291429239d-2*u88
  u72=-1.50588379858285848d-1*u85
  u351=2.93911700706593106d5*u367
  u79=u55*(u351+u359*u354)*u345
  u84=u240*(8.9153215880999909d6*u367+u359*(2.1d1*u344+u363))*u345
  u5=-2.5098063309714308d-2*u84
  u354=2.5098063309714308d-2*u75
  u99=u354*u237
  u343=u5+u99
  u70=u237*(u343)
  value=value+(c(62)*(u239*(u40+u237*(u237*(u60+u64)+u58)+(u237*(u34+u9&
 &5)+u38)*u239)+u238*(u98+u237*(u70+u69)+(u237*(u72+u99)+u86)*u238)+u23&
 &7*(u4+u79*u237)+u41))
  u379=u240*(4.45320758646353191d4*u367+u359)*u345
  u53=4.51765139574857544d-1*u379
  u350=1.46955850353296553d6*u367
  u93=u240*(u350+u359*(4.5d1+u364))*u345
  u46=-7.52941899291429239d-2*u93
  u80=2.25882569787428772d-1*u85
  u371=u240*(u49+u359*(u89+u91))*u345
  u49=1.00392253238857232d-1*u371
  value=value+(c(63)*(u54*(u238*(u46+u237*(u355+u49)+(u355+u80)*u238)+u&
 &237*(u46+u80*u237)+u53)))
  u53=9.03530279149715087d-1*u240*(1.68232286599733428d4*u367+u360)*u34&
 &5
  u49=4.89852834510988511d5*u367
  u89=u240*(u49+u359*u366)*u345
  u47=-1.50588379858285848d-1*u89
  u375=-8.28236089220572163d-1*u71
  u48=+u358*(3.36365613030878777d6*u367+u359*(1.03d2+1.2d1*pd2))*u345
  u364=u240*(2.12269561621428355d6*u367+u359*(5.0d0*u344+u339))*u345
  u349=2.00784506477714464d-1*u364
  u378=-5.0196126619428616d-2*u75
  u81=1.75686443168000156d-1*u85
  u341=u378*u237
  value=value+(c(64)*(u380*(u239*(u47+u237*(u341+u349)+(u341+u61)*u239)&
 &+u238*(u375+u377+(u64+u60)*u238)+u237*(u48+u81*u237)+u53)))
  u373=-2.25882569787428772d-1*u240
  u339=+u373*(3.0677652262304331d4*u367+u359)*u345
  u46=5.27059329504000468d-1*u71
  u63=7.52941899291429239d-2*u89
  u91=-7.52941899291429239d-2*u85
  u82=u55*(1.33893108099670193d6*u367+u359*(4.1d1+6.0d0*pd2))*u345
  u35=-5.0196126619428616d-2*u85
  u366=-1.00392253238857232d-1*u364
  u364=u237*(u99+u366)
  value=value+(c(65)*(u54*(u239*(u46+u237*(u64+u35)+(u95+u72)*u239)+u23&
 &8*(u63+u364+(u99+u91)*u238)+u237*(u82+u35*u237)+u339)))
  u80=-2.25882569787428772d-1*u379
  u374=1.50588379858285848d-1*u240
  u94=u374*(u49+u359*u6)*u345
  u49=-2.25882569787428772d-1*u85
  u7=-5.0196126619428616d-2*u371
  u6=7.52941899291429239d-2*u75
  u371=u6*u237
  value=value+(c(66)*(u380*(u239*(u94+u7*u237+(u371+u49)*u239)+u63*u237&
 &+u80)))
  value=value+(c(67)*(u238*(u51+u237*(u78*u237+u44)+u238*((u60+u355)*u2&
 &38+u376+u45))+u237*(u83+u342*u237)+u372))
  u372=u378*u238
  value=value+(c(68)*(u353*(u239*(u47+u238*(u372+u349)+(u372+u61)*u239)&
 &+u238*(u48+u377+(u64+u81)*u238)+u237*(u375+u348)+u53)))
  u375=u374*(6.23449062104894468d4*u367+u346)*u345
  u83=+u36*(u351-1.13437046457953026d2*u57)*u345
  u51=5.0196126619428616d-2*u85
  u53=-7.2282422331977207d0*u240*(9.8960168588078487d2*u367-3.544907701&
 &81103205d0*u8)*u345
  u52=u55*(6.85793968315383915d5*u367+u359*(2.1d1+u352))*u345
  u50=-2.5098063309714308d-2*u85
  u351=1.48440252882117731d4*u367
  u349=+u87*(u351+u360)*u345
  u377=-1.50588379858285848d-1*u88
  u374=5.0196126619428616d-2*u84
  u84=1.50588379858285848d-1*u88
  u8=+u358*(1.14625563275571312d7*u367+u359*(2.7d1*u344+u363))*u345
  u376=2.5098063309714308d-2*u89
  value=value+(c(69)*(u239*(u375+u237*(u68*u237+u377)+u239*((u51+u341)*&
 &u239+u237*(u374+u341)+u83))+u238*(u53+u237*(u237*(u8+u99)+u84)+u238*(&
 &(u50+u99)*u238+u237*(u8+u95)+u52))+u237*(u349+u237*(u50*u237+u376))))
  u375=+u373*(3.66152623775890402d4*u367+u359)*u345
  u83=u361*(1.01236252465604292d6*u367+u359*(3.1d1+3.0d0*pd2))*u345
  u51=-1.2549031654857154d-1*u85
  u53=3.01176759716571696d-1*u71
  u52=-3.51372886336000312d-1*u85
  u349=u354*u238
  value=value+(c(70)*(u353*(u239*(u83+u237*(u64+u52)+u238*(u349+u366)+(&
 &u349+u95+u51)*u239)+u238*(u63+u91*u238)+u237*(u53+u348)+u375)))
  u95=u55*(u351+1.13437046457953026d2*exppd2)*u345
  u348=+u365*(1.33596227593905957d5*u367+u346)*u345
  u84=u55*(u350+u359*(4.5d1+u352))*u345
  u351=-7.52941899291429239d-2*u379
  u350=7.52941899291429239d-2*u93
  u373=+u358*(1.91042605459285519d7*u367+u359*(4.5d1*u344+u363))*u345
  value=value+(c(71)*(u239*(u348+u350*u237+u239*((u91+u371)*u239+u373*u&
 &237+u84))+u351*u237+u95))
  u358=u9*u238
  u371=(u62+u238*(u358+u90))*u239
  u359=u39*u238
  u344=u359+u64
  u350=u59*u237
  u351=u97*u237
  value=value+(c(72)*(u54*(u371+u238*(u370+u350+u238*(u344+u77))+u351+u&
 &362)))
  u373=u61*u237
  u374=u98*u237
  value=value+(c(73)*(u380*((u38+u238*(u358+u34))*u239+u238*(u74+u373+u&
 &238*(u344+u68))+u374+u368)))
  u378=u237*(u63+u91*u237)
  value=value+(c(74)*(u54*(u239*(u46+u238*(u359+u35)+(u358+u72)*u239)+u&
 &238*(u82+u364+(u99+u35)*u238)+u378+u339)))
  value=value+(c(75)*(u380*(u239*(u83+u364+u238*(u359+u52)+(u358+u99+u5&
 &1)*u239)+u238*(u53+u60*u238)+u378+u375)))
  u364=-8.94135087490894218d2*u56
  u51=1.35529541872457263d0*u71
  u53=-6.02353519433143391d-1*u85
  u375=u9*u239
  u378=u375+u344
  u376=u61*u238
  u377=u98*u238
  value=value+(c(76)*(u54*(u239*(u51+u373+u376+u239*(u378+u53))+u377+u3&
 &74+u364)))
  u373=-1.11766885936361777d4*u56
  u366=4.26667076265143236d0*u71
  u54=-9.53726405769143704d-1*u85
  u374=(u239*(u366+u350+u59*u238+u239*(u378+u54))+u97*u238+u351+u373)
  value=value+(c(77)*(u380*u374))
  u52=u64+u359
  u378=u340*u237
  u379=u65*u237
  u380=u67*u237
  value=value+(c(78)*((u96+u238*(u238*(u43+u358)+u92))*u239+u238*(u73+u&
 &378+u238*(u238*(u76+u52)+u379+u369))+u380+u42))
  value=value+(c(79)*(u353*(u371+u238*(u37+u350+u238*(u344+u66))+u351+u&
 &357)))
  u344=u237*(u72*u237+u69)
  u350=u237*(u98+u86*u237)
  value=value+(c(80)*(u239*(u40+u238*(u238*(u60+u359)+u58)+(u238*(u34+u&
 &358)+u38)*u239)+u238*(u4+u344+u238*(u237*(u343+u349)+u79))+u350+u41))
  value=value+(c(81)*(u353*(u239*(u94+u7*u238+(u6*u238+u49)*u239)+u63*u&
 &238+u80)))
  value=value+(c(82)*(u239*(u348+u344+u238*(u376+u51)+u239*((u91+u99+u3&
 &58)*u239+u238*(u53+u359)+u70+u84))+u238*(u364+u377)+u350+u95))
  value=value+(c(83)*(u353*u374))
  u353=5.0196126619428616d-1*u240*(u347+u356)*u345
  u240=-5.58834429681808887d4*u56
  u56=9.78824469078858012d0*u71
  u71=-1.35529541872457263d0*u85
  value=value+(c(84)*(u239*(u240+u378+u340*u238+u239*(u239*(u71+u52+u37&
 &5)+u65*u238+u379+u56))+u67*u238+u380+u353))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_opt_2_0=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_Y_opt_2_1(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_Y_opt_2_1
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u126
  real(kind=8) :: u127
  real(kind=8) :: u178
  real(kind=8) :: u53
  real(kind=8) :: u54
  real(kind=8) :: u55
  real(kind=8) :: u56
  real(kind=8) :: u57
  real(kind=8) :: u58
  real(kind=8) :: u59
  real(kind=8) :: u6
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
  D_X_Y_opt_2_1=.true.
  
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
  u6=a3**(-2)
  u99=8.19411322944059304d-1*u6*A3_(p,pd,erfpd,exppd2)*p**3
  u127=x*z
  value=value+(c(1)*(u99*u127))
  if ( lmax .eq. 0 ) return
  u99=A4_(p,pd,erfpd,exppd2)
  u126=p**4
  u77=-8.19411322944059304d-1*u6*u99*u126
  u57=8.69422416480109528d-2*u6
  u73=exppd2*pd
  u93=-1.41796308072441282d1*u73
  u97=u57*(4.71238898038468986d1*u99+u93)*u126
  u178=x**2
  value=value+(c(2)*((u97*u178+u77)*z))
  u126=x*y
  u87=u126*z
  value=value+(c(3)*(u97*u87))
  u99=z**2
  value=value+(c(4)*(x*(u97*u99+u77)))
  if ( lmax .eq. 1 ) return
  u97=A5_(p,pd,erfpd,exppd2)
  u98=p**5
  u90=u6*u97*u98
  u78=-1.22911698441608896d1*u90
  u77=exppd2*pd2
  u63=-2.83592616144882564d1*u77
  u83=u57*(3.2986722862692829d2*u97+u63)*u98
  u95=u83*u178
  value=value+(c(5)*(x*(u95+u78)*z))
  u74=-4.09705661472029652d0*u90
  u90=y*z
  value=value+(c(6)*((u95+u74)*u90))
  u59=1.41796308072441282d1*exppd2
  u91=1.73884483296021906d-2*u6*(4.71238898038468986d1*u97+u59)*u98
  value=value+(c(7)*((u74+u95)*u99+u74*u178+u91))
  u97=y**2
  value=value+(c(8)*(x*(u83*u97+u74)*z))
  u98=u83*u99
  value=value+(c(9)*(u126*(u98+u74)))
  value=value+(c(10)*(u127*(u98+u78)))
  if ( lmax .eq. 2 ) return
  u98=A6_(p,pd,erfpd,exppd2)
  u95=p**6
  u91=u6*u98*u95
  u83=1.22911698441608896d1*u91
  u8=-2.83592616144882564d1*u73
  u68=u6*(3.2986722862692829d2*u98+u8)*u95
  u80=-5.21653449888065717d-1*u68
  u78=2.0d0*pd2
  u85=u57*(2.96880505764235461d3*u98+u8*(9.0d0+u78))*u95
  u95=u85*u178
  value=value+(c(11)*((u178*(u80+u95)+u83)*z))
  u8=-2.60826724944032859d-1*u68
  value=value+(c(12)*(x*(u95+u8)*u90))
  u67=-8.69422416480109528d-2*u68
  u60=u67*u178
  value=value+(c(13)*(x*((u8+u95)*u99+u60+u83)))
  u84=4.09705661472029652d0*u91
  u88=(u67+u95)
  u95=u60+u84
  value=value+(c(14)*((u88*u97+u95)*z))
  u60=u88*u99
  value=value+(c(15)*(y*(u60+u95)))
  value=value+(c(16)*(z*(u60+u8*u178+u83)))
  u95=u85*u97
  value=value+(c(17)*(u126*(u95+u8)*z))
  value=value+(c(18)*(x*((u67+u95)*u99+u67*u97+u84)))
  u95=u85*u99
  value=value+(c(19)*(u87*(u95+u8)))
  value=value+(c(20)*(x*(u99*(u80+u95)+u83)))
  if ( lmax .eq. 3 ) return
  u95=A7_(p,pd,erfpd,exppd2)
  u60=p**7
  u88=u6*u95*u60
  u85=4.30190944545631135d2*u88
  u83=-5.67185232289765129d1*u77
  u67=u6*(2.96880505764235461d3*u95+u83)*u60
  u8=-8.69422416480109529d-1*u67
  u84=(1.1d1+u78)
  u80=u6*(3.26568556340659007d4*u95+u83*u84)*u60
  u92=8.69422416480109528d-2*u80
  u89=u92*u178
  value=value+(c(21)*(x*(u178*(u8+u89)+u85)*z))
  u82=8.60381889091262269d1*u88
  u55=-5.21653449888065717d-1*u67
  value=value+(c(22)*((u178*(u55+u89)+u82)*u90))
  u58=-3.72609607062904084d-2*u6*(2.83592616144882564d1*exppd2+3.298672&
 &2862692829d2*u95)*u60
  u86=1.72076377818252454d2*u88
  u64=-8.69422416480109528d-2*u67
  u65=u64*u178
  value=value+(c(23)*((u82+u178*(u89+u55))*u99+u178*(u86+u65)+u58))
  u7=-2.60826724944032859d-1*u67
  u96=(u7+u89)
  u75=u65+u82
  value=value+(c(24)*(x*(u96*u97+u75)*z))
  u65=u96*u99
  value=value+(c(25)*(u126*(u65+u75)))
  u96=2.58114566727378681d2*u88
  u75=u7*u178
  value=value+(c(26)*(u127*(u65+u75+u96)))
  u96=(u64+u89)
  u65=u75+u82
  value=value+(c(27)*(y*(u96*u97+u65)*z))
  u75=2.48406404708602722d-2*u6*(9.8960168588078487d2*u95-u59)*u60
  u59=8.69422416480109528d-2*u67
  u67=3.47768966592043811d-1*u6
  u66=u67*(8.90641517292706383d3*u95+u63*(6.0d0+pd2))*u60
  u60=-8.69422416480109528d-2*u80
  u63=u60*u178
  value=value+(c(28)*(u97*(u64+u178*(u63+u66)+(u63+u59)*u97)+u178*(u64+&
 &u59*u178)+u75))
  value=value+(c(29)*(u90*(u96*u99+u65)))
  value=value+(c(30)*(u99*(u86+u55*u178+(u89+u64)*u99)+u82*u178+u58))
  u96=u92*u97
  value=value+(c(31)*(x*(u97*(u55+u96)+u82)*z))
  value=value+(c(32)*(u126*((u7+u96)*u99+u64*u97+u82)))
  value=value+(c(33)*(u127*((u64+u96)*u99+u7*u97+u82)))
  u96=u92*u99
  value=value+(c(34)*(u126*(u99*(u55+u96)+u82)))
  value=value+(c(35)*(u127*(u99*(u8+u96)+u85)))
  if ( lmax .eq. 4 ) return
  u96=A8_(p,pd,erfpd,exppd2)
  u65=p**8
  u89=u6*u96*u65
  u75=-4.30190944545631135d2*u89
  u63=-5.67185232289765129d1*u73
  u66=u6*(2.96880505764235461d3*u96+u63)*u65
  u60=3.91240087416049288d0*u66
  u92=3.26568556340659007d4*u96
  u80=u6*(u92+u63*u84)*u65
  u59=-1.30413362472016429d0*u80
  u58=pd2**2
  u85=4.0d0*u58
  u64=(u78+1.1d1)
  u7=u6*(4.24539123242856709d5*u96+u63*(1.3d1*u64+u85))*u65
  u55=8.69422416480109528d-2*u7
  u8=u55*u178
  value=value+(c(36)*((u178*(u60+u178*(u8+u59))+u75)*z))
  u84=1.30413362472016429d0*u66
  u71=-8.69422416480109529d-1*u80
  value=value+(c(37)*(x*(u178*(u71+u8)+u84)*u90))
  u61=8.69422416480109529d-1*u66
  u69=-8.69422416480109528d-2*u80
  u53=u69*u178
  value=value+(c(38)*(x*((u84+u178*(u8+u71))*u99+u178*(u61+u53)+u75)))
  u72=-8.60381889091262269d1*u89
  u86=2.60826724944032859d-1*u66
  u95=5.21653449888065717d-1*u66
  u94=-5.21653449888065717d-1*u80
  u79=(u86+u178*(u8+u94))
  u62=u178*(u95+u53)+u72
  value=value+(c(39)*((u79*u97+u62)*z))
  u53=u79*u99
  value=value+(c(40)*(y*(u53+u62)))
  u79=-2.58114566727378681d2*u89
  u62=1.56496034966419715d0*u66
  u9=-2.60826724944032859d-1*u80
  u54=u9*u178
  value=value+(c(41)*(z*(u53+u178*(u62+u54)+u79)))
  u53=7.82480174832098576d-1*u66
  u76=(u9+u8)
  u66=u54+u53
  value=value+(c(42)*(u126*(u76*u97+u66)*z))
  u54=2.08661379955226287d0*u6*(3.2986722862692829d2*u96-7.089815403622&
 &06411d0*u73)*u65
  u88=2.60826724944032859d-1*u80
  u73=(1.3d1+u78)
  u70=-8.69422416480109528d-2*u6*(3.85944657493506099d4*u96+u63*u73)*u6&
 &5
  u96=1.39107586636817525d0*u6*(u92+u93*(4.0d0*u64+u58))*u65
  u65=-8.69422416480109528d-2*u7
  u64=8.69422416480109528d-2*u80
  u93=u65*u178
  u80=u178*(u93+u96)
  value=value+(c(43)*(x*(u97*(u9+u80+(u93+u88)*u97)+u178*(u70+u64*u178)&
 &+u54)))
  value=value+(c(44)*(u87*(u76*u99+u66)))
  u76=u94*u178
  u66=u86*u178
  value=value+(c(45)*(x*(u99*(u62+u76+(u8+u9)*u99)+u66+u79)))
  u63=(u8+u69)
  u8=u95+u76
  u76=u66+u72
  value=value+(c(46)*((u97*(u8+u63*u97)+u76)*z))
  u62=u70+u80
  u80=u178*(u9+u88*u178)+u54
  u66=(u93+u64)
  value=value+(c(47)*(y*(u97*(u62+u66*u97)+u80)))
  value=value+(c(48)*(z*(u99*(u62+u66*u99)+u80)))
  u80=u63*u99
  value=value+(c(49)*(y*(u99*(u8+u80)+u76)))
  value=value+(c(50)*(z*(u99*(u61+u71*u178+u80)+u84*u178+u75)))
  u80=u55*u97
  value=value+(c(51)*(u126*(u97*(u71+u80)+u84)*z))
  value=value+(c(52)*(x*((u86+u97*(u80+u94))*u99+u97*(u95+u69*u97)+u72)&
 &))
  value=value+(c(53)*(u87*((u9+u80)*u99+u9*u97+u53)))
  value=value+(c(54)*(x*(u99*(u95+u94*u97+(u80+u69)*u99)+u86*u97+u72)))
  u94=u55*u99
  value=value+(c(55)*(u87*(u99*(u71+u94)+u84)))
  value=value+(c(56)*(x*(u99*(u60+u99*(u94+u59))+u75)))
  if ( lmax .eq. 5 ) return
  u94=A9_(p,pd,erfpd,exppd2)
  u72=p**9
  u87=u6*u94*u72
  u80=-2.71020295063747615d4*u87
  u66=-1.13437046457953026d2*u77
  u65=u6*(3.26568556340659007d4*u94+u66)*u72
  u63=9.12893537304115006d0*u65
  u96=u6*(4.24539123242856709d5*u94+u66*u73)*u72
  u62=-1.82578707460823001d0*u96
  u54=(u78+1.3d1)
  u8=u6*(6.36808684864285064d6*u94+u66*(1.5d1*u54+u85))*u72
  u9=8.69422416480109528d-2*u8
  u88=u9*u178
  value=value+(c(57)*(x*(u178*(u63+u178*(u88+u62))+u80)*z))
  u71=-3.87171850091068022d3*u87
  u74=3.91240087416049288d0*u65
  u95=-1.30413362472016429d0*u96
  value=value+(c(58)*((u178*(u74+u178*(u88+u95))+u71)*u90))
  u93=2.96880505764235461d3*u94
  u70=u6*(u93+5.67185232289765129d1*exppd2)*u72
  u69=1.44903736080018255d-1*u70
  u76=-1.16151555027320407d4*u87
  u59=1.30413362472016429d0*u65
  u60=-8.69422416480109528d-2*u96
  u73=u60*u178
  value=value+(c(59)*((u71+u178*(u178*(u95+u88)+u74))*u99+u178*(u76+u17&
 &8*(u73+u59))+u69))
  u84=8.69422416480109529d-1*u65
  u75=-8.69422416480109529d-1*u96
  u64=(u59+u178*(u88+u75))
  u53=u178*(u84+u73)+u71
  value=value+(c(60)*(x*(u64*u97+u53)*z))
  u73=u64*u99
  value=value+(c(61)*(u126*(u73+u53)))
  u64=2.60826724944032859d0*u65
  u53=-2.60826724944032859d-1*u96
  u55=u53*u178
  value=value+(c(62)*(u127*(u73+u178*(u64+u55)+u76)))
  u73=-2.32303110054640813d3*u87
  u86=2.60826724944032859d-1*u65
  u91=1.56496034966419715d0*u65
  u98=-5.21653449888065717d-1*u96
  u61=(u86+u178*(u88+u98))
  u81=u178*(u91+u55)+u73
  value=value+(c(63)*(y*(u61*u97+u81)*z))
  u55=-2.31845977728029207d-1*u6*(u93-7.08981540362206411d0*exppd2)*u72
  u82=-2.60826724944032859d-1*u65
  u56=2.60826724944032859d-1*u6*(6.23449062104894468d4*u94-2.2687409291&
 &5906051d2*u77)*u72
  u68=-2.60826724944032859d-1*u6
  u93=+u68*(8.81735102119779319d5*u94+u66*(2.7d1+4.0d0*pd2))*u72
  u79=5.21653449888065717d-1*u96
  u7=-1.73884483296021906d-1*u6*(2.93911700706593106d5*u94+u66*(9.0d0+p&
 &d2))*u72
  u89=u57*(8.9153215880999909d6*u94+u66*(2.1d1*u54+u85))*u72
  u85=-8.69422416480109528d-2*u8
  u92=8.69422416480109528d-2*u96
  u77=u85*u178
  u8=u178*(u89+u77)
  value=value+(c(64)*(u97*(u86+u178*(u8+u93)+(u178*(u79+u77)+u82)*u97)+&
 &u178*(u56+u178*(u92*u178+u7))+u55))
  value=value+(c(65)*(u90*(u61*u99+u81)))
  u61=8.69422416480109528d-2*u70
  u81=-4.64606220109281626d3*u87
  u87=3.1299206993283943d0*u65
  u65=u98*u178
  u57=u86*u178
  value=value+(c(66)*(u99*(u81+u178*(u65+u87)+(u178*(u98+u88)+u86)*u99)&
 &+u178*(u81+u57)+u61))
  u81=(u88+u53)
  u87=u91+u65
  u61=u57+u73
  value=value+(c(67)*(x*(u97*(u87+u81*u97)+u61)*z))
  u57=1.56496034966419715d0*u6*(1.48440252882117731d4*u94+u83)*u72
  u83=+u68*(4.89852834510988511d5*u94+u66*(1.5d1+u78))*u72
  u68=2.60826724944032859d-1*u96
  u6=u67*(2.12269561621428355d6*u94+u66*(5.0d0*u54+u58))*u72
  u54=u83+u178*(u77+u6)
  u58=u178*(u83+u68*u178)+u57
  u67=(u77+u68)
  value=value+(c(68)*(u126*(u97*(u54+u67*u97)+u58)))
  value=value+(c(69)*(u127*(u99*(u54+u67*u99)+u58)))
  u67=u81*u99
  value=value+(c(70)*(u126*(u99*(u87+u67)+u61)))
  u81=u75*u178
  u61=u59*u178
  value=value+(c(71)*(u127*(u99*(u64+u81+u67)+u61+u76)))
  u64=(u88+u60)
  u87=u84+u81
  u81=u61+u71
  value=value+(c(72)*(y*(u97*(u87+u64*u97)+u81)*z))
  u61=u56+u178*(u79*u178+u93)
  u67=u8+u7
  u8=u178*(u86+u82*u178)+u55
  u54=(u92+u77)
  value=value+(c(73)*(u97*(u61+u97*(u54*u97+u67))+u8))
  u77=u85*u97
  value=value+(c(74)*(u90*(u99*(u83+u97*(u77+u6)+(u77+u68)*u99)+u97*(u8&
 &3+u68*u97)+u57)))
  value=value+(c(75)*(u99*(u61+u99*(u54*u99+u67))+u8))
  value=value+(c(76)*(u90*(u99*(u87+u64*u99)+u81)))
  value=value+(c(77)*(u99*(u76+u74*u178+u99*((u60+u88)*u99+u95*u178+u59&
 &))+u71*u178+u69))
  u178=u9*u97
  value=value+(c(78)*(x*(u97*(u74+u97*(u178+u95))+u71)*z))
  value=value+(c(79)*(u126*((u59+u97*(u178+u75))*u99+u97*(u84+u60*u97)+&
 &u71)))
  value=value+(c(80)*(u127*((u86+u97*(u178+u98))*u99+u97*(u91+u53*u97)+&
 &u73)))
  value=value+(c(81)*(u126*(u99*(u91+u98*u97+(u178+u53)*u99)+u86*u97+u7&
 &3)))
  value=value+(c(82)*(u127*(u99*(u84+u75*u97+(u178+u60)*u99)+u59*u97+u7&
 &1)))
  u178=u9*u99
  value=value+(c(83)*(u126*(u99*(u74+u99*(u178+u95))+u71)))
  value=value+(c(84)*(u127*(u99*(u63+u99*(u178+u62))+u80)))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_opt_2_1=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_Y_opt_2_2(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_Y_opt_2_2
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u151
  real(kind=8) :: u152
  real(kind=8) :: u228
  real(kind=8) :: u284
  real(kind=8) :: u285
  real(kind=8) :: u286
  real(kind=8) :: u287
  real(kind=8) :: u288
  real(kind=8) :: u289
  real(kind=8) :: u29
  real(kind=8) :: u290
  real(kind=8) :: u291
  real(kind=8) :: u292
  real(kind=8) :: u3
  real(kind=8) :: u30
  real(kind=8) :: u31
  real(kind=8) :: u32
  real(kind=8) :: u33
  real(kind=8) :: u34
  real(kind=8) :: u35
  real(kind=8) :: u36
  real(kind=8) :: u37
  real(kind=8) :: u38
  real(kind=8) :: u39
  real(kind=8) :: u4
  real(kind=8) :: u40
  real(kind=8) :: u41
  real(kind=8) :: u42
  real(kind=8) :: u43
  real(kind=8) :: u44
  real(kind=8) :: u45
  real(kind=8) :: u46
  real(kind=8) :: u47
  real(kind=8) :: u48
  real(kind=8) :: u49
  real(kind=8) :: u5
  real(kind=8) :: u50
  real(kind=8) :: u51
  real(kind=8) :: u52
  real(kind=8) :: u53
  real(kind=8) :: u54
  real(kind=8) :: u55
  real(kind=8) :: u56
  real(kind=8) :: u57
  real(kind=8) :: u58
  real(kind=8) :: u59
  real(kind=8) :: u6
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
  D_X_Y_opt_2_2=.true.
  
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
  u228=a3**(-2)
  u151=u228*A3_(p,pd,erfpd,exppd2)*p**3
  u73=-4.09705661472029652d-1*u151
  u87=4.09705661472029652d-1*u151
  u151=x**2
  u152=y**2
  value=value+(c(1)*(u73*u152+u87*u151))
  if ( lmax .eq. 0 ) return
  u73=A4_(p,pd,erfpd,exppd2)
  u87=p**4
  u3=u228*u73*u87
  u68=-8.19411322944059304d-1*u3
  u61=-1.41796308072441282d1*exppd2
  u63=+u61*pd
  u98=u228*(4.71238898038468986d1*u73+u63)*u87
  u73=-4.34711208240054764d-2*u98
  u87=4.34711208240054764d-2*u98
  u35=u73*u152+u87*u151
  value=value+(c(2)*(x*(u35+u68)))
  u68=8.19411322944059304d-1*u3
  value=value+(c(3)*(y*(u35+u68)))
  value=value+(c(4)*((u35)*z))
  if ( lmax .eq. 1 ) return
  u68=A5_(p,pd,erfpd,exppd2)
  u3=4.71238898038468986d1*u68
  u35=p**5
  u88=1.73884483296021906d-2*u228*(u3-u61)*u35
  u95=u228*u68*u35
  u92=2.04852830736014826d0*u95
  u86=-1.02426415368007413d1*u95
  u54=+u86
  u284=-2.83592616144882564d1*exppd2
  u78=+u284*pd2
  u99=u228*(3.2986722862692829d2*u68+u78)*u35
  u68=-4.34711208240054764d-2*u99
  u82=4.34711208240054764d-2*u99
  u73=u82*u151
  value=value+(c(5)*((u92+u68*u151)*u152+u151*(u54+u73)+u88))
  u99=x*y
  u54=u68*u152
  u88=u54+u73
  value=value+(c(6)*(u99*(u88)))
  u97=-4.09705661472029652d0*u95
  u49=+u97
  value=value+(c(7)*(x*(u88+u49)*z))
  u49=-1.73884483296021906d-2*u228*(-u61+u3)*u35
  u3=-u86
  u35=-u92
  u86=u73+u54
  u54=u35*u151
  value=value+(c(8)*(u152*(u3+u86)+u54+u49))
  u73=-u97
  value=value+(c(9)*(y*(u88+u73)*z))
  u73=z**2
  value=value+(c(10)*((u86)*u73+u92*u152+u54))
  if ( lmax .eq. 2 ) return
  u49=A6_(p,pd,erfpd,exppd2)
  u35=p**6
  u95=u228*u49*u35
  u3=2.45823396883217791d1*u95
  u92=+u284*pd
  u97=u228*(3.2986722862692829d2*u49+u92)*u35
  u88=1.30413362472016429d-1*u97
  u86=-3.91240087416049288d-1*u97
  u68=2.0d0*pd2
  u54=u228*(2.96880505764235461d3*u49+u92*(9.0d0+u68))*u35
  u49=-4.34711208240054764d-2*u54
  u35=4.34711208240054764d-2*u54
  u92=u35*u151
  u94=u49*u151
  value=value+(c(11)*(x*((u88+u94)*u152+u151*(u86+u92)+u3)))
  u86=4.34711208240054764d-2*u97
  u74=-1.30413362472016429d-1*u97
  u48=(u86+u94)*u152
  value=value+(c(12)*(y*(u48+u151*(u74+u92))))
  u94=4.09705661472029652d0*u95
  u72=-2.17355604120027382d-1*u97
  value=value+(c(13)*((u48+u151*(u72+u92)+u94)*z))
  u72=-4.34711208240054764d-2*u97
  u48=u49*u152
  u70=u92+u48
  u285=u72*u151
  value=value+(c(14)*(x*(u152*(u88+u70)+u285)))
  value=value+(c(15)*(u99*(u48+u92)*z))
  u92=-8.69422416480109528d-2*u97
  u48=u86*u152+u285
  value=value+(c(16)*(x*((u92+u70)*u73+u48+u94)))
  u92=-u3
  u3=3.91240087416049288d-1*u97
  u96=u74*u151
  value=value+(c(17)*(y*(u152*(u3+u70)+u96+u92)))
  u92=-u94
  u3=2.17355604120027382d-1*u97
  value=value+(c(18)*((u152*(u3+u70)+u285+u92)*z))
  u3=8.69422416480109528d-2*u97
  value=value+(c(19)*(y*((u3+u70)*u73+u48+u92)))
  value=value+(c(20)*(z*((u70)*u73+u88*u152+u96)))
  if ( lmax .eq. 3 ) return
  u3=A7_(p,pd,erfpd,exppd2)
  u92=3.2986722862692829d2*u3
  u95=p**7
  u88=-7.45219214125808167d-2*u228*(-u284+u92)*u95
  u94=u228*u3*u95
  u70=-4.30190944545631135d1*u94
  u48=+u70
  u285=5.59248227909320475d2*u94
  u49=2.96880505764235461d3*u3
  u86=-5.67185232289765129d1*exppd2
  u35=+u86*pd2
  u96=u228*(u49+u35)*u95
  u89=2.60826724944032859d-1*u96
  u84=-6.0859569153607667d-1*u96
  u43=(1.1d1+u68)
  u55=u228*(3.26568556340659007d4*u3+u35*u43)*u95
  u76=-4.34711208240054764d-2*u55
  u34=4.34711208240054764d-2*u55
  u54=u34*u151
  u9=u76*u151
  value=value+(c(21)*((u48+u151*(u9+u89))*u152+u151*(u285+u151*(u54+u84&
 &))+u88))
  u88=8.60381889091262269d1*u94
  u84=1.30413362472016429d-1*u96
  u4=-3.04297845768038335d-1*u96
  u41=(u84+u9)*u152
  value=value+(c(22)*(u99*(u41+u151*(u4+u54)+u88)))
  u4=1.72076377818252454d2*u94
  u85=-3.91240087416049288d-1*u96
  value=value+(c(23)*(x*(u41+u151*(u85+u54)+u4)*z))
  u85=4.34711208240054764d-2*u96
  u41=-u70
  u70=-4.34711208240054764d-2*u96
  u55=x**4
  u66=u70*u151
  value=value+(c(24)*(u152*(u48+u34*u55+(u9+u85)*u152)+u151*(u41+u66)))
  u71=-1.30413362472016429d-1*u96
  value=value+(c(25)*(y*((u85+u9)*u152+u151*(u71+u54))*z))
  u46=-6.21016011771506806d-3*u228*(-u86+u49)*u95
  u67=u228*(4.94800842940392435d3*u3+u35)*u95
  u74=4.34711208240054764d-2*u67
  u75=-2.60826724944032859d-1*u96
  u62=u228*(8.90641517292706383d3*u3+u78*(6.0d0+pd2))*u95
  u3=-1.73884483296021906d-1*u62
  u78=-8.69422416480109528d-2*u96
  value=value+(c(26)*((u41+u151*(u54+u75))*u73+u152*(u85+u151*(u54+u3)+&
 &(u54+u70)*u152)+u151*(u74+u78*u151)+u46))
  u46=-u88
  u74=3.04297845768038335d-1*u96
  u3=u76*u152
  u51=u54+u3
  u286=u71*u151
  value=value+(c(27)*(u99*(u152*(u74+u51)+u286+u46)))
  value=value+(c(28)*(x*(u152*(u84+u51)+u66)*z))
  value=value+(c(29)*(u99*((u51)*u73+u85*u152+u66)))
  u74=x*z
  u66=u84*u152+u286
  value=value+(c(30)*(u74*((u78+u51)*u73+u66+u88)))
  u88=7.45219214125808167d-2*u228*(u92-u284)*u95
  u92=-u285
  u284=6.0859569153607667d-1*u96
  u285=u3+u54
  u78=u75*u151
  u54=u41*u151
  value=value+(c(31)*(u152*(u92+u78+u152*(u285+u284))+u54+u88))
  u92=-u4
  u94=3.91240087416049288d-1*u96
  value=value+(c(32)*(y*(u152*(u94+u51)+u286+u92)*z))
  u94=6.21016011771506806d-3*u228*(u49-u86)*u95
  u88=-4.34711208240054764d-2*u67
  u92=8.69422416480109528d-2*u96
  u49=1.73884483296021906d-1*u62
  value=value+(c(33)*((u48+u152*(u3+u89))*u73+u152*(u88+u151*(u9+u49)+(&
 &u9+u92)*u152)+u151*(u70+u85*u151)+u94))
  u85=y*z
  value=value+(c(34)*(u85*((u92+u51)*u73+u66+u46)))
  value=value+(c(35)*(u73*(u78+u89*u152+(u285)*u73)+u48*u152+u54))
  if ( lmax .eq. 4 ) return
  u88=A8_(p,pd,erfpd,exppd2)
  u89=p**8
  u48=u228*u88*u89
  u46=-1.29057283363689341d3*u48
  u41=+u46
  u94=2.96880505764235461d3*u88
  u49=+u86*pd
  u95=u228*(u94+u49)*u89
  u76=-6.52066812360082147d-1*u95
  u284=4.12975647828052026d0*u95
  u4=3.26568556340659007d4*u88
  u51=u228*(u4+u49*u43)*u89
  u43=4.34711208240054764d-1*u51
  u78=-8.69422416480109529d-1*u51
  u96=pd2**2
  u54=4.0d0*u96
  u285=(u68+1.1d1)
  u9=u228*(4.24539123242856709d5*u88+u49*(1.3d1*u285+u54))*u89
  u67=-4.34711208240054764d-2*u9
  u3=4.34711208240054764d-2*u9
  u286=u3*u151
  u75=u67*u151
  value=value+(c(36)*(x*((u76+u151*(u75+u43))*u152+u151*(u284+u151*(u28&
 &6+u78))+u41)))
  u41=-8.60381889091262269d1*u48
  u78=+u41
  u66=-1.30413362472016429d-1*u95
  u289=1.17372026224814786d0*u95
  u39=2.60826724944032859d-1*u51
  u6=-5.21653449888065717d-1*u51
  u84=(u66+u151*(u75+u39))*u152
  value=value+(c(37)*(y*(u84+u151*(u289+u151*(u286+u6))+u78)))
  u6=-1.72076377818252454d2*u48
  u37=+u6
  u50=1.69537371213621358d0*u95
  u288=-6.0859569153607667d-1*u51
  value=value+(c(38)*((u84+u151*(u50+u151*(u286+u288))+u37)*z))
  u37=1.30413362472016429d-1*u51
  u288=3.04297845768038335d-1*u95
  u84=-1.73884483296021906d-1*u51
  u53=-4.34711208240054764d-2*u51
  u34=u53*u151
  value=value+(c(39)*(x*(u152*(u66+u151*(u286+u84)+(u75+u37)*u152)+u151&
 &*(u288+u34)+u78)))
  u288=2.60826724944032859d-1*u95
  u8=-3.04297845768038335d-1*u51
  value=value+(c(40)*(u99*((u37+u75)*u152+u151*(u8+u286)+u288)*z))
  u288=u228*(4.28827397215006777d3*u88+u49)*u89
  u8=-1.30413362472016429d-1*u288
  u90=6.52066812360082147d-1*u95
  u77=-1.30413362472016429d-1*u51
  u56=u228*(6.8282516325774156d4*u88+u49*(2.3d1+u68))*u89
  u44=4.34711208240054764d-2*u56
  u29=-4.34711208240054764d-1*u51
  u40=u228*(u4+u63*(4.0d0*u285+u96))*u89
  u4=-6.95537933184087623d-1*u40
  u62=-8.69422416480109528d-2*u51
  u63=u151*(u286+u4)
  value=value+(c(41)*(x*((u90+u151*(u286+u29))*u73+u152*(u37+u63+(u286+&
 &u77)*u152)+u151*(u44+u62*u151)+u8)))
  u8=-u41
  u44=-3.04297845768038335d-1*u95
  u41=4.34711208240054764d-2*u51
  u80=1.30413362472016429d-1*u95
  u45=1.73884483296021906d-1*u51
  u93=(u75+u41)
  u52=u93*u152
  u32=u77*u151
  value=value+(c(42)*(y*(u152*(u44+u151*(u286+u45)+u52)+u151*(u80+u32)+&
 &u8)))
  value=value+(c(43)*((u152*(u66+u3*u55+u52)+u151*(u80+u34))*z))
  u44=(1.3d1+u68)
  u52=u228*(3.85944657493506099d4*u88+u49*u44)*u89
  u34=4.34711208240054764d-2*u52
  u38=1.30413362472016429d-1*u52
  u69=-2.60826724944032859d-1*u51
  u36=(u286+u53)
  value=value+(c(44)*(y*((u80+u151*(u286+u69))*u73+u152*(u34+u63+u36*u1&
 &52)+u151*(u38+u84*u151)+u66)))
  u34=u228*(3.62853951489621119d3*u88+u49)*u89
  u38=-1.30413362472016429d-1*u34
  u4=(8.0d0+pd2)
  u63=u228*(u94-7.08981540362206411d0*exppd2*pd*u4)*u89
  u94=6.95537933184087623d-1*u63
  u87=1.30413362472016429d-1*u228
  u84=u87*(5.04696859799200284d4*u88+u49*(1.7d1+u68))*u89
  u70=2.0d0*u96
  u79=-8.69422416480109528d-2*u228*(3.59225411974724908d5*u88+u49*(1.1d&
 &1*u285+u70))*u89
  u88=8.69422416480109528d-2*u9
  u9=u69*u151
  value=value+(c(45)*(z*(u73*(u94+u151*(u88*u151+u79)+u36*u73)+u151*(u8&
 &4+u9)+u38)))
  u89=-u289
  u94=5.21653449888065717d-1*u51
  u84=u67*u152
  u79=u84+u286
  u36=u80*u151
  value=value+(c(46)*(x*(u152*(u89+u9+u152*(u79+u94))+u36+u8)))
  u89=-2.60826724944032859d-1*u95
  u94=3.04297845768038335d-1*u51
  u53=u286+u84
  value=value+(c(47)*(u99*(u152*(u94+u53)+u32+u89)*z))
  u89=-1.30413362472016429d-1*u52
  u94=-4.34711208240054764d-2*u52
  u52=6.95537933184087623d-1*u40
  u40=u151*(u75+u52)
  u289=u152*(u84+u39)
  value=value+(c(48)*(x*((u66+u289)*u73+u152*(u89+u40+(u75+u45)*u152)+u&
 &151*(u94+u41*u151)+u80)))
  u89=u37*u152
  value=value+(c(49)*(u99*z*((u53)*u73+u89+u32)))
  u94=5.21653449888065717d-1*u95
  u53=u9+u39*u152
  u88=u66*u152+u36
  value=value+(c(50)*(x*(u73*(u94+u53+(u79+u62)*u73)+u88+u78)))
  u94=-u46
  u45=-u284
  u46=8.69422416480109529d-1*u51
  u284=u29*u151
  u285=u90*u151
  value=value+(c(51)*(y*(u152*(u45+u284+u152*(u79+u46))+u285+u94)))
  u94=-u6
  u45=-u50
  u46=6.0859569153607667d-1*u51
  value=value+(c(52)*((u152*(u45+u9+u152*(u79+u46))+u36+u94)*z))
  u94=1.30413362472016429d-1*u288
  u50=-4.34711208240054764d-2*u56
  u46=8.69422416480109528d-2*u51
  u45=u151*(u77+u37*u151)
  value=value+(c(53)*(y*((u76+u152*(u84+u43))*u73+u152*(u50+u40+(u75+u4&
 &6)*u152)+u45+u94)))
  u94=1.30413362472016429d-1*u34
  u50=-6.95537933184087623d-1*u63
  u34=-7.82480174832098576d-1*u95
  value=value+(c(54)*(z*(u73*(u50+u40+u289+u93*u73)+u152*(u34+u89)+u45+&
 &u94)))
  u52=-5.21653449888065717d-1*u95
  value=value+(c(55)*(y*(u73*(u52+u53+(u79+u46)*u73)+u88+u8)))
  value=value+(c(56)*(z*(u73*(u284+u43*u152+(u79)*u73)+u76*u152+u285)))
  if ( lmax .eq. 5 ) return
  u43=A9_(p,pd,erfpd,exppd2)
  u285=2.96880505764235461d3*u43
  u76=p**9
  u67=u228*(u285-u86)*u76
  u52=4.34711208240054764d-1*u67
  u34=u228*u43*u76
  u46=1.93585925045534011d3*u34
  u95=-4.83964812613835028d4*u34
  u37=+u95
  u284=3.26568556340659007d4*u43
  u45=-1.13437046457953026d2*exppd2
  u88=+u45*pd2
  u94=u228*(u284+u88)*u76
  u78=-1.95620043708024644d0*u94
  u77=+u78
  u84=8.4768685606810679d0*u94
  u50=u228*(4.24539123242856709d5*u43+u88*u44)*u76
  u44=6.52066812360082147d-1*u50
  u62=-1.17372026224814786d0*u50
  u69=+u62
  u48=6.36808684864285064d6*u43
  u90=(u68+1.3d1)
  u6=1.5d1*u90
  u8=u228*(u48+u88*(u6+u54))*u76
  u80=-4.34711208240054764d-2*u8
  u66=4.34711208240054764d-2*u8
  u9=u66*u151
  u75=u80*u151
  value=value+(c(57)*((u46+u151*(u151*(u44+u75)+u77))*u152+u151*(u37+u1&
 &51*(u151*(u69+u9)+u84))+u52))
  u52=-7.74343700182136044d3*u34
  u37=+u52
  u69=-6.52066812360082147d-1*u94
  u98=3.26033406180041074d0*u94
  u40=4.34711208240054764d-1*u50
  u289=-7.82480174832098576d-1*u50
  u288=(u69+u151*(u75+u40))*u152
  value=value+(c(58)*(u99*(u288+u151*(u98+u151*(u9+u289))+u37)))
  u37=-1.16151555027320407d4*u34
  u289=+u37
  u41=4.12975647828052026d0*u94
  u79=-8.69422416480109529d-1*u50
  value=value+(c(59)*(x*(u288+u151*(u41+u151*(u9+u79))+u289)*z))
  u289=2.89807472160036509d-2*u67
  u79=3.87171850091068022d2*u34
  u67=-1.30413362472016429d-1*u94
  u288=-3.4845466508196122d3*u34
  u53=+u288
  u93=3.91240087416049288d-1*u94
  u36=2.60826724944032859d-1*u50
  u89=5.21653449888065717d-1*u94
  u29=-3.91240087416049288d-1*u50
  u49=-4.34711208240054764d-2*u50
  u286=u49*u151
  value=value+(c(60)*(u152*(u79+u151*(u151*(u29+u9)+u93)+(u151*(u36+u75&
 &)+u67)*u152)+u151*(u53+u151*(u286+u89))+u289))
  u289=-7.74343700182136043d2*u34
  u53=+u289
  u89=1.17372026224814786d0*u94
  u29=-5.21653449888065717d-1*u50
  value=value+(c(61)*(y*((u67+u151*(u75+u36))*u152+u151*(u89+u151*(u9+u&
 &29))+u53)*z))
  u3=3.85944657493506099d4*u43
  u32=1.44903736080018255d-2*u228
  u39=2.26874092915906051d2*exppd2
  u38=u32*(u3+u39)*u76
  u97=-u46
  u81=1.30413362472016429d-1*u94
  u92=u228*(2.67192455187811915d4*u43+u35)*u76
  u83=-5.21653449888065717d-1*u92
  u91=-u78
  u78=4.0d0*pd2
  u57=u228*(8.81735102119779319d5*u43+u88*(2.7d1+u78))*u76
  u72=1.30413362472016429d-1*u57
  u63=-2.60826724944032859d-1*u50
  u58=u228*(1.07767623592417472d6*u43+u88*(3.3d1+u68))*u76
  u292=4.34711208240054764d-2*u58
  u71=-6.52066812360082147d-1*u50
  u51=u228*(8.9153215880999909d6*u43+u88*(2.1d1*u90+u54))*u76
  u47=-4.34711208240054764d-2*u51
  u30=-8.69422416480109528d-2*u50
  value=value+(c(62)*((u97+u151*(u151*(u71+u9)+u91))*u73+u152*(u67+u151&
 &*(u151*(u47+u9)+u72)+(u151*(u63+u9)+u81)*u152)+u151*(u83+u151*(u30*u1&
 &51+u292))+u38))
  u38=1.30413362472016429d-1*u50
  u83=6.52066812360082147d-1*u94
  u72=-1.30413362472016429d-1*u50
  u292=(u75+u38)
  u47=u292*u152
  u290=u72*u151
  value=value+(c(63)*(u99*(u152*(u69+u66*u55+u47)+u151*(u83+u290))))
  u64=3.04297845768038335d-1*u94
  u31=-1.73884483296021906d-1*u50
  value=value+(c(64)*(x*(u152*(u67+u151*(u9+u31)+u47)+u151*(u64+u286)+u&
 &53)*z))
  u291=u228*(3.46360590058274705d4*u43+u88)*u76
  u64=-3.91240087416049288d-1*u291
  u53=u228*(4.89852834510988511d5*u43+u88*(1.5d1+u68))*u76
  u47=1.30413362472016429d-1*u53
  u59=u228*(1.79612705987362454d6*u43+u88*(5.5d1+6.0d0*pd2))*u76
  u42=4.34711208240054764d-2*u59
  u56=-4.34711208240054764d-1*u50
  u286=u228*(2.12269561621428355d6*u43+u88*(5.0d0*u90+u96))*u76
  u5=-1.73884483296021906d-1*u286
  u65=(u9+u72)
  value=value+(c(65)*(u99*((u83+u151*(u9+u56))*u73+u152*(u47+u151*(u9+u&
 &5)+u65*u152)+u151*(u42+u31*u151)+u64)))
  u31=u228*(4.45320758646353191d4*u43+u88)*u76
  u42=-3.91240087416049288d-1*u31
  u64=u228*(1.63284278170329504d5*u43+u35*(1.0d1+pd2))*u76
  u7=5.21653449888065717d-1*u64
  u60=u228*(8.16421390851647518d5*u43+u88*(2.5d1+u68))*u76
  u33=1.30413362472016429d-1*u60
  u82=u228*(u48+u88*(u6+u70))*u76
  u48=-8.69422416480109528d-2*u82
  u6=8.69422416480109528d-2*u8
  u70=u63*u151
  u287=u6*u151
  value=value+(c(66)*(u74*(u73*(u7+u151*(u287+u48)+u65*u73)+u151*(u33+u&
 &70)+u42)))
  u96=u228*(-u86+u285)*u76
  u7=-2.89807472160036509d-2*u96
  u33=-u288
  u48=-5.21653449888065717d-1*u94
  u42=4.34711208240054764d-2*u50
  u288=-u79
  u79=-3.91240087416049288d-1*u94
  u65=3.91240087416049288d-1*u50
  u86=(u42+u75)
  u285=u81*u151
  value=value+(c(67)*(u152*(u33+u151*(u70+u79)+u152*(u86*u152+u151*(u65&
 &+u9)+u48))+u151*(u288+u285)+u7))
  u7=-u289
  u288=-3.04297845768038335d-1*u94
  u48=1.73884483296021906d-1*u50
  value=value+(c(68)*(y*(u152*(u288+u151*(u9+u48)+(u75+u42)*u152)+u151*&
 &(u81+u290)+u7)*z))
  u288=u228*(1.48440252882117731d4*u43+u35)*u76
  u65=-2.60826724944032859d-1*u288
  u33=4.34711208240054764d-2*u53
  u289=2.60826724944032859d-1*u288
  u288=-4.34711208240054764d-2*u53
  u35=(u49+u9)
  value=value+(c(69)*(u152*(u65+u55*(u75+u66)+u152*(u35*u152+u75+u33))+&
 &u151*(u289+u151*(u42*u151+u288))))
  u65=u228*(u284+u61*pd2*u4)*u76
  u33=2.08661379955226287d0*u65
  u289=7.82480174832098576d-1*u94
  u288=u66*u152
  value=value+(c(70)*(u85*(u73*(u33+u151*(u9+u63)+u152*(u288+u5)+(u288+&
 &u72)*u73)+u152*(u47+u72*u152)+u151*(u289+u290)+u79)))
  u289=u32*(u284-u45)*u76
  u288=u228*(8.01577365563435745d4*u43-u39*pd2)*u76
  u32=-1.30413362472016429d-1*u288
  u61=u228*(6.85793968315383915d5*u43+u88*(2.1d1+u68))*u76
  u47=4.34711208240054764d-2*u61
  u68=-1.30413362472016429d-1*u228*(5.04696859799200284d4*u43+u88)*u76
  u290=u87*(1.27361736972857013d6*u43+u88*(3.9d1+u78))*u76
  u78=-4.34711208240054764d-2*u228*(1.14625563275571312d7*u43+u88*(2.7d&
 &1*u90+u54))*u76
  u43=2.60826724944032859d-1*u94
  value=value+(c(71)*(u73*(u32+u151*(u29*u151+u290)+u73*(u35*u73+u151*(&
 &u78+u287)+u47))+u151*(u68+u43*u151)+u289))
  u289=-u52
  u32=-u98
  u49=7.82480174832098576d-1*u50
  u52=u80*u152
  u98=u52+u9
  u287=u56*u151
  u29=u83*u151
  value=value+(c(72)*(u99*(u152*(u32+u287+u152*(u98+u49))+u29+u289)))
  u32=-u89
  u289=5.21653449888065717d-1*u50
  value=value+(c(73)*(x*(u152*(u32+u70+u152*(u98+u289))+u285+u7)*z))
  u89=3.91240087416049288d-1*u291
  u289=-4.34711208240054764d-2*u59
  u32=-1.30413362472016429d-1*u53
  u7=1.73884483296021906d-1*u286
  u49=u151*(u75+u7)
  u290=u151*(u32+u38*u151)
  value=value+(c(74)*(u99*((u69+u152*(u52+u40))*u73+u152*(u289+u49+(u75&
 &+u48)*u152)+u290+u89)))
  u291=-u33
  u289=-7.82480174832098576d-1*u94
  value=value+(c(75)*(u74*(u73*(u291+u49+u152*(u52+u36)+u292*u73)+u152*&
 &(u289+u38*u152)+u290+u93)))
  u32=u36*u152
  u289=u67*u152
  value=value+(c(76)*(u99*(u73*(u70+u32+(u98)*u73)+u289+u285)))
  u70=-3.87171850091068022d3*u34
  u33=+u70
  u7=8.69422416480109529d-1*u94
  u285=u287+u40*u152
  u63=u69*u152+u29
  value=value+(c(77)*(u74*(u73*(u7+u285+(u98+u30)*u73)+u63+u33)))
  u33=-4.34711208240054764d-1*u96
  u7=-u95
  u49=-u84
  u84=-u62
  u62=u9+u52
  u290=u91*u151
  u291=u71*u151
  u292=u97*u151
  value=value+(c(78)*(u152*(u7+u290+u152*(u152*(u84+u62)+u291+u49))+u29&
 &2+u33))
  u7=-u37
  u49=-u41
  u33=8.69422416480109529d-1*u50
  value=value+(c(79)*(y*(u152*(u49+u287+u152*(u98+u33))+u29+u7)*z))
  u287=-1.44903736080018255d-2*u228
  u29=+u287*(u39+u3)*u76
  u3=5.21653449888065717d-1*u92
  u92=-4.34711208240054764d-2*u58
  u49=8.69422416480109528d-2*u50
  u33=-1.30413362472016429d-1*u57
  u7=4.34711208240054764d-2*u51
  u51=u151*(u36*u151+u33)
  u228=u151*(u7+u75)
  u57=u151*(u81+u67*u151)
  value=value+(c(80)*((u46+u152*(u152*(u44+u52)+u77))*u73+u152*(u3+u51+&
 &u152*((u49+u75)*u152+u228+u92))+u57+u29))
  u29=3.91240087416049288d-1*u31
  u31=-5.21653449888065717d-1*u64
  u64=-1.30413362472016429d-1*u60
  u60=8.69422416480109528d-2*u82
  u82=-8.69422416480109528d-2*u8
  value=value+(c(81)*(u85*(u73*(u31+u152*(u82*u152+u60)+(u52+u38)*u73)+&
 &u152*(u64+u32)+u29)))
  u29=+u287*(-u45+u284)*u76
  u284=1.30413362472016429d-1*u288
  u288=-4.34711208240054764d-2*u61
  u61=2.32303110054640813d3*u34
  u287=-1.56496034966419715d0*u94
  value=value+(c(82)*(u73*(u284+u51+u152*(u32+u287)+u73*(u86*u73+u152*(&
 &u36+u52)+u228+u288))+u152*(u61+u289)+u57+u29))
  u228=-u70
  u34=-8.69422416480109529d-1*u94
  value=value+(c(83)*(u85*(u73*(u34+u285+(u98+u49)*u73)+u63+u228)))
  value=value+(c(84)*(u73*(u290+u77*u152+u73*((u62)*u73+u44*u152+u291))&
 &+u46*u152+u292))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_opt_2_2=.false.
end function

  
  
!> Compute D X Y coulomb integral 
!!
recursive function D_X_Y_opt_2(a12,a3,r3,m3,c,lmax,value)
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12      !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3       !< exponent of the third solid harmonics 
  integer     , intent(in)  :: m3       !< m index of the solid harmonic
  real(kind=8), intent(in)  :: r3(3)    !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)     !< hermit coefficients of the two fisrt obrital produt
  integer     , intent(in)  :: lmax     !< l of maximum non 0 coeff in hermit decomposition
  real(kind=8), intent(out) :: value    !< result 
  logical                   :: D_X_Y_opt_2 !< true if the routine knows how to calculate the requested value
  
  ! init return value
  D_X_Y_opt_2 = .true.
  ! select case on m component
  select case (m3)
    case (-2)
      D_X_Y_opt_2 = D_X_Y_opt_2_m2(a12,a3,r3,c,lmax,value)
    case (-1)
      D_X_Y_opt_2 = D_X_Y_opt_2_m1(a12,a3,r3,c,lmax,value)
    case (0)
      D_X_Y_opt_2 = D_X_Y_opt_2_0(a12,a3,r3,c,lmax,value)
    case (1)
      D_X_Y_opt_2 = D_X_Y_opt_2_1(a12,a3,r3,c,lmax,value)
    case (2)
      D_X_Y_opt_2 = D_X_Y_opt_2_2(a12,a3,r3,c,lmax,value)
    case default
      D_X_Y_opt_2 = .false.
  end select
end function
  
end module
