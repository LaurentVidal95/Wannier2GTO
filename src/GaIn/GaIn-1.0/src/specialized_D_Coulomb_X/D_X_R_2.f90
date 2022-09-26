module mod_D_X_R_opt_2
  
  use mod_A_functions
  
  implicit none
  
  private
  public :: D_X_R_opt_2
  
contains
  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_R_opt_2_1(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_R_opt_2_1
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u290
  real(kind=8) :: u356
  real(kind=8) :: u385
  real(kind=8) :: u404
  real(kind=8) :: u406
  real(kind=8) :: u411
  real(kind=8) :: u412
  real(kind=8) :: u416
  real(kind=8) :: u417
  real(kind=8) :: u418
  real(kind=8) :: u419
  real(kind=8) :: u42
  real(kind=8) :: u420
  real(kind=8) :: u423
  real(kind=8) :: u424
  real(kind=8) :: u425
  real(kind=8) :: u426
  real(kind=8) :: u427
  real(kind=8) :: u428
  real(kind=8) :: u429
  real(kind=8) :: u43
  real(kind=8) :: u430
  real(kind=8) :: u431
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
  D_X_R_opt_2_1=.true.
  
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
  u65=A3_(p,pd,erfpd,exppd2)
  u426=p**3
  u290=u65*u426
  u406=a3**(-2)
  u63=-4.70157986289796906d-2*u406*(5.31736155271654808d0*u290+2.0d0*(2&
 &.0d0*exppd2*u426-5.31736155271654808d0*a3*A1_(p,pd,erfpd,exppd2)))
  u431=7.5d-1*u406*u290
  u290=x**2
  value=value+(c(1)*(u431*u290+u63))
  if ( lmax .eq. 0 ) return
  u63=A4_(p,pd,erfpd,exppd2)
  u431=p**2
  u7=u63*u431
  u5=-2.0d0*a3
  u66=+u5*A2_(p,pd,erfpd,exppd2)
  u83=-2.5d-1*u406*u431
  u67=+u83*(9.0d0*u7+u66)
  u43=2.65868077635827404d1*u63
  u9=1.41047395886939072d-1*u406
  u80=p**4
  u430=exppd2*pd
  u58=-8.0d0*u430
  u76=u9*(u43+u58)*u80
  u416=u76*u290
  value=value+(c(2)*(x*(u416+u67)))
  u67=+u83*(3.0d0*u7+u66)
  u66=(u416+u67)
  value=value+(c(3)*(u66*y))
  value=value+(c(4)*(u66*z))
  if ( lmax .eq. 1 ) return
  u76=4.0d0*exppd2
  u66=-5.0d0*a3*(u76+5.31736155271654808d0*u65)
  u83=A5_(p,pd,erfpd,exppd2)
  u7=u83*u431
  u53=u406*u426
  u426=9.40315972579593812d-3*u53
  u87=exppd2*u431
  u81=u426*(2.0d0*(3.6d1*u87+u66)+2.39281269872244664d2*u7)
  u356=1.5d1*u7
  u68=-a3*u65
  u65=-1.5d0*u53*(u356+u68)
  u99=1.86107654345079183d2*u83
  u92=p**5
  u416=exppd2*pd2
  u67=-1.6d1*u416
  u88=(u99+u67)*u92
  u429=u9*u88
  u412=u429*u290
  value=value+(c(5)*(u290*(u65+u412)+u81))
  u65=-2.0d0*(-u68)
  u81=-7.5d-1*u53
  u68=+u81*(u356+u65)
  u356=x*(u412+u68)
  value=value+(c(6)*(u356*y))
  value=value+(c(7)*(u356*z))
  u53=u426*(2.0d0*(1.2d1*u87+u66)+7.97604232907482212d1*u7)
  u66=+u81*(5.0d0*u7+u65)
  u7=u83*u92
  u65=-3.75d0*u406*u7
  u356=y**2
  u81=(u66+u412)
  u426=u65*u290+u53
  value=value+(c(8)*(u81*u356+u426))
  u68=y*z
  value=value+(c(9)*((u412+u66)*u68))
  u412=z**2
  value=value+(c(10)*(u81*u412+u426))
  if ( lmax .eq. 2 ) return
  u429=A6_(p,pd,erfpd,exppd2)
  u65=u429*u431
  u81=+u5*u63
  u63=u406*u80
  u80=2.25d0*u63*(2.5d1*u65+u81)
  u5=1.86107654345079183d2*u429
  u426=-1.6d1*u430
  u53=(u5+u426)
  u78=u53*u431
  u82=a3*(-u58-u43)
  u43=-2.82094791773878144d-1*u63
  u47=+u43*(u82+5.0d0*u78)
  u423=1.67496888910571265d3*u429
  u404=2.0d0*pd2
  u93=+u426*(9.0d0+u404)
  u42=p**6
  u62=(u423+u93)*u42
  u85=u9*u62
  u411=u85*u290
  value=value+(c(11)*(x*(u290*(u47+u411)+u80)))
  u80=7.5d-1*u63*(1.5d1*u65+u81)
  u47=+u43*(u82+3.0d0*u78)
  u43=(u290*(u47+u411)+u80)
  value=value+(c(12)*(u43*y))
  value=value+(c(13)*(u43*z))
  u43=-1.41047395886939072d-1*u63
  u78=u430*u431
  u89=+u43*(5.58322963035237549d2*u65+2.0d0*(u82-2.4d1*u78))
  u70=u53*u42
  u96=u406*u70
  u74=-1.41047395886939072d-1*u96
  u385=(u89+u411)
  u417=u74*u290
  u64=u417+u80
  value=value+(c(14)*(x*(u385*u356+u64)))
  value=value+(c(15)*(x*(u411+u89)*u68))
  value=value+(c(16)*(x*(u385*u412+u64)))
  u89=u63*(5.0d0*u65+u81)
  u81=2.25d0*u89
  u63=+u43*(1.86107654345079183d2*u65+2.0d0*(u82-8.0d0*u78))
  u43=-4.23142187660817215d-1*u96
  u96=(u63+u411)
  u65=u96*u356
  u64=u43*u290+u81
  value=value+(c(17)*(y*(u65+u64)))
  u82=7.5d-1*u89
  u385=u417+u82
  value=value+(c(18)*((u65+u385)*z))
  u65=u96*u412
  value=value+(c(19)*(y*(u65+u385)))
  value=value+(c(20)*(z*(u65+u64)))
  if ( lmax .eq. 3 ) return
  u43=8.0d0*exppd2
  u74=-7.0d0*a3*(u43+2.65868077635827404d1*u83)
  u81=A7_(p,pd,erfpd,exppd2)
  u63=u81*u431
  u65=u406*u92
  u92=-1.20897767903090633d-2*u65
  u85=+u92*(4.65269135862697957d3*u63+2.0d0*(2.0d2*u87+u74))
  u96=a3*u83
  u82=1.125d1*u65
  u83=u82*(1.05d2*u63-4.0d0*u96)
  u64=a3*(-u67-u99)
  u99=-1.41047395886939072d-1*u65
  u80=u416*u431
  u52=-2.4d2*u80
  u61=+u99*(2.51245333365856897d4*u63+2.0d0*(u64+u52))
  u411=1.84246577801628391d4*u81
  u8=(1.1d1+u404)
  u66=-3.2d1*u416
  u417=+u66*u8
  u385=(u411+u417)
  u89=p**7
  u47=u385*u89
  u420=u406*u47
  u60=1.41047395886939072d-1*u420
  u55=u60*u290
  value=value+(c(21)*(u290*(u83+u290*(u55+u61))+u85))
  u85=-2.0d0*u96
  u83=u82*(3.5d1*u63+u85)
  u61=1.67496888910571265d3*u81
  u82=(u61+u66)
  u44=u82*u431
  u54=-2.82094791773878144d-1*u65
  u48=+u54*(u64+5.0d0*u44)
  u84=x*(u290*(u48+u55)+u83)
  value=value+(c(22)*(u84*y))
  value=value+(c(23)*(u84*z))
  u84=-4.02992559676968777d-3*u65*(2.79161481517618774d3*u63+2.0d0*(1.2&
 &d2*u87+u74))
  u46=2.1d1*u63
  u90=u65*(u46+u85)
  u79=3.75d0*u90
  u91=7.5d0*u65*(u46-u96)
  u46=+u54*(u64+3.0d0*u44)
  u48=u82*u89
  u44=u406*u48
  u54=-1.41047395886939072d-1*u44
  u45=(u79+u290*(u55+u46))
  u418=u54*u290
  u72=u290*(u91+u418)+u84
  value=value+(c(24)*(u45*u356+u72))
  value=value+(c(25)*((u290*(u46+u55)+u79)*u68))
  value=value+(c(26)*(u45*u412+u72))
  u84=1.125d1*u90
  u45=+u99*(5.02490666731713794d3*u63+2.0d0*(u64-4.8d1*u80))
  u72=-4.23142187660817215d-1*u44
  u90=x*y
  u77=(u45+u55)
  u71=u77*u356
  u419=u72*u290
  u50=u419+u84
  value=value+(c(27)*(u90*(u71+u50)))
  u46=u418+u79
  value=value+(c(28)*(x*(u71+u46)*z))
  u71=u77*u412
  value=value+(c(29)*(u90*(u71+u46)))
  u46=x*z
  value=value+(c(30)*(u46*(u71+u50)))
  u45=+u92*(9.30538271725395914d2*u63+2.0d0*(4.0d1*u87+u74))
  u92=u65*(7.0d0*u63+u85)
  u85=2.25d1*u92
  u418=1.67496888910571265d3*u63
  u63=-1.6d1*u80
  u71=+u99*(u418+2.0d0*(u64+u63))
  u79=u81*u89
  u64=7.875d1*u406*u79
  u77=-8.46284375321634431d-1*u44
  u50=(u55+u71)
  u98=u85+u77*u290
  u424=u64*u290+u45
  value=value+(c(31)*(u356*(u98+u50*u356)+u424))
  u86=1.125d1*u92
  u49=(u71+u55)
  u55=u419+u86
  value=value+(c(32)*(y*(u49*u356+u55)*z))
  u419=5.58322963035237549d2*u81
  u425=8.05985119353937553d-3*u65*(5.0d0*(u419-u43)*u431-u74)
  u84=a3**(-1)
  u65=-7.5d0*u84*u7
  u7=+u99*(u418+2.0d0*(2.65868077635827404d1*u96+u63))
  u92=2.82094791773878144d-1*u84
  u96=u92*u88
  u88=1.41047395886939072d-1*u44
  u91=5.02490666731713794d3*u81
  u418=5.64189583547756287d-1*u406
  u99=u418*(u91+u67*(6.0d0+pd2))*u89
  u63=-1.41047395886939072d-1*u420
  u420=u63*u290
  value=value+(c(33)*((u65+u96*u356)*u412+u356*(u7+u290*(u420+u99)+(u42&
 &0+u88)*u356)+u290*(u54+u88*u290)+u425))
  value=value+(c(34)*(u68*(u49*u412+u55)))
  value=value+(c(35)*(u412*(u98+u50*u412)+u424))
  if ( lmax .eq. 4 ) return
  u425=A8_(p,pd,erfpd,exppd2)
  u65=u425*u431
  u7=a3*u429
  u96=-2.0d0*u7
  u99=u406*u42
  u63=-5.625d1*u99*(4.9d1*u65+u96)
  u60=a3*(-u426-u5)
  u54=7.05236979434695359d-1*u99*(4.0d0*(u60-1.68d2*u78)+3.517434667121&
 &99656d4*u65)
  u72=a3*(-u426*(u404+9.0d0)-u423)
  u77=-1.41047395886939072d-1*u99
  u88=u430*u8
  u71=u88*u431
  u64=+u77*(3.86917813383419621d5*u65+2.0d0*(u72-3.36d2*u71))
  u86=pd2**2
  u85=4.0d0*u86
  u420=p**8
  u8=(u404+1.1d1)
  u83=-3.2d1*u430
  u44=u406*(2.39520551142116908d5*u425+u83*(1.3d1*u8+u85))*u420
  u74=1.41047395886939072d-1*u44
  u42=u74*u290
  value=value+(c(36)*(x*(u290*(u54+u290*(u42+u64))+u63)))
  u63=-1.125d1*u99*(3.5d1*u65+u96)
  u54=4.23142187660817215d-1*u99
  u64=u54*(4.0d0*(u60-1.2d2*u78)+2.51245333365856897d4*u65)
  u5=+u77*(2.76369866702442587d5*u65+2.0d0*(u72-2.4d2*u71))
  u50=(u290*(u64+u290*(u42+u5))+u63)
  value=value+(c(37)*(u50*y))
  value=value+(c(38)*(u50*z))
  u64=u54*(2.0d0*(u60-8.0d1*u78)+8.37484444552856323d3*u65)
  u54=(1.67496888910571265d3*u425+u83)
  u50=u54*u431
  u49=2.82094791773878144d-1*u99*(5.0d0*u50+u60)
  u98=1.84246577801628391d4*u425
  u424=(u98-3.2d1*u88)
  u88=u424*u431
  u55=a3*(-u93-u423)
  u423=-2.82094791773878144d-1*u99
  u93=+u423*(u55+5.0d0*u88)
  u45=u406*u424*u420
  u424=-1.41047395886939072d-1*u45
  u429=(u64+u290*(u42+u93))
  u56=u424*u290
  u95=u290*(u49+u56)+u63
  value=value+(c(39)*(x*(u429*u356+u95)))
  value=value+(c(40)*(x*(u290*(u93+u42)+u64)*u68))
  value=value+(c(41)*(x*(u429*u412+u95)))
  u93=u99*(2.1d1*u65+u96)
  u64=-1.125d1*u93
  u429=u99*(2.0d0*(u60-4.8d1*u78)+5.02490666731713794d3*u65)
  u95=1.41047395886939072d-1*u429
  u57=u99*(3.0d0*u50+u60)
  u50=8.46284375321634431d-1*u57
  u428=+u423*(u55+3.0d0*u88)
  u88=-4.23142187660817215d-1*u45
  u55=(u95+u290*(u42+u428))
  u423=u55*u356
  u97=u88*u290
  u73=u290*(u50+u97)+u64
  value=value+(c(42)*(y*(u423+u73)))
  u69=-3.75d0*u93
  u51=2.82094791773878144d-1*u57
  u75=u290*(u51+u56)+u69
  value=value+(c(43)*((u423+u75)*z))
  u423=u55*u412
  value=value+(c(44)*(y*(u423+u75)))
  value=value+(c(45)*(z*(u423+u73)))
  u428=8.46284375321634431d-1*u429
  u49=5.52739733404885173d4*u65
  u95=-4.8d1*u71
  u51=+u77*(u49+2.0d0*(u72+u95))
  u69=u406*u54*u420
  u54=4.23142187660817215d-1*u69
  u73=-8.46284375321634431d-1*u45
  u75=(u42+u51)
  u423=u73*u290
  u55=u428+u423
  u424=u54*u290
  u56=u424+u64
  value=value+(c(46)*(x*(u356*(u55+u75*u356)+u56)))
  u50=4.23142187660817215d-1*u429
  u429=(u51+u42)
  u427=u97+u50
  value=value+(c(47)*(u90*(u429*u356+u427)*z))
  u97=4.0d0*(1.86107654345079183d2*u425-4.0d0*u430)*u431
  u430=8.46284375321634431d-1*u99
  u94=u430*(u97+8.86226925452758014d0*u7)
  u6=u84*u70
  u70=-2.82094791773878144d-1*u6
  u59=a3*u53
  u53=+u77*(u49+2.0d0*(u59+u95))
  u49=u92*u62
  u62=4.23142187660817215d-1*u45
  u95=(1.3d1+u404)
  u5=-1.41047395886939072d-1*u406*(2.17745955583742644d4*u425+u83*u95)*&
 &u420
  u425=2.25675833419102515d0*u406*(u98+u58*(4.0d0*u8+u86))*u420
  u98=-1.41047395886939072d-1*u44
  u63=1.41047395886939072d-1*u45
  u57=u98*u290
  u58=u290*(u57+u425)
  u420=u49*u356
  value=value+(c(48)*(x*((u70+u420)*u412+u356*(u53+u58+(u57+u62)*u356)+&
 &u290*(u5+u63*u290)+u94)))
  value=value+(c(49)*(x*u68*(u429*u412+u427)))
  value=value+(c(50)*(x*(u412*(u55+u75*u412)+u56)))
  u94=u99*(7.0d0*u65+u96)
  u70=-5.625d1*u94
  u5=u99*(2.0d0*(u60-1.6d1*u78)+1.67496888910571265d3*u65)
  u93=1.41047395886939072d0*u5
  u51=+u77*(1.84246577801628391d4*u65+2.0d0*(u72-1.6d1*u71))
  u72=2.11571093830408608d0*u69
  u56=-1.41047395886939072d0*u45
  u75=(u42+u51)
  u60=u75*u356
  u429=u93+u56*u290
  u64=u72*u290+u70
  value=value+(c(51)*(y*(u356*(u429+u60)+u64)))
  u71=-1.125d1*u94
  u99=8.46284375321634431d-1*u5
  u5=u99+u423
  u50=u424+u71
  value=value+(c(52)*((u356*(u5+u60)+u50)*z))
  u60=u430*(u97+2.65868077635827404d1*u7)
  u97=-8.46284375321634431d-1*u6
  u6=u95*u431
  u7=+u77*(2.17745955583742644d4*u65+2.0d0*(u59+u426*u6))
  u430=u7+u58
  u59=u290*(u88+u62*u290)+u60
  u424=(u57+u63)
  value=value+(c(53)*(y*((u97+u420)*u412+u356*(u430+u424*u356)+u59)))
  value=value+(c(54)*(z*(u412*(u430+u420+u424*u412)+u97*u356+u59)))
  u60=u75*u412
  value=value+(c(55)*(y*(u412*(u5+u60)+u50)))
  value=value+(c(56)*(z*(u412*(u429+u60)+u64)))
  if ( lmax .eq. 5 ) return
  u64=1.86107654345079183d2*u81
  u60=1.6d1*exppd2
  u88=a3*(u60+u64)
  u51=-9.0d0*u88
  u5=A9_(p,pd,erfpd,exppd2)
  u420=u5*u431
  u50=u406*u89
  u71=3.35827133064140647d-2*u50
  u97=u71*(2.0d0*(7.84d2*u87+u51)+8.20734755661799196d4*u420)
  u89=-a3*u81
  u72=-2.3625d3*u50*(4.2d1*u420+u89)
  u63=-6.4d1*u416
  u49=(1.84246577801628391d4*u5+u63)
  u54=u49*u431
  u81=a3*(-u66-u61)
  u429=4.23142187660817215d0*u50*(7.0d0*u54+u81)
  u56=(2.39520551142116908d5*u5+u63*u95)
  u430=u56*u431
  u98=a3*(-u417-u411)
  u70=-2.82094791773878144d-1*u50
  u94=+u70*(u98+1.4d1*u430)
  u424=p**9
  u423=(u404+1.3d1)
  u93=u406*(3.59280826713175363d6*u5+u63*(1.5d1*u423+u85))*u424
  u99=1.41047395886939072d-1*u93
  u83=u99*u290
  value=value+(c(57)*(u290*(u72+u290*(u290*(u94+u83)+u429))+u97))
  u94=-2.0d0*(-u89)
  u72=-3.9375d2*u50*(6.3d1*u420+u94)
  u429=7.05236979434695359d-1*u50*(4.0d0*(u81-3.36d2*u80)+3.86917813383&
 &419621d5*u420)
  u417=a3*(-u66*u8-u411)
  u97=5.02993157398445508d6*u420
  u411=-1.41047395886939072d-1*u50
  u7=u416*u6
  u6=+u411*(u97+2.0d0*(u417-6.72d2*u7))
  u95=x*(u290*(u429+u290*(u83+u6))+u72)
  value=value+(c(58)*(u95*y))
  value=value+(c(59)*(u95*z))
  u429=6.71654266128281294d-3*u50*(2.0d0*(5.6d2*u87+u51)+5.862391111869&
 &99426d4*u420)
  u95=u50*(4.5d1*u420+u94)
  u62=-7.875d1*u95
  u73=-7.875d1*u50*(1.35d2*u420-4.0d0*(-u89))
  u6=2.76369866702442587d5*u420
  u72=4.23142187660817215d-1*u50
  u74=u72*(4.0d0*(u81+u52)+u6)
  u69=1.41047395886939072d-1*u50
  u59=u69*(2.0d0*(u81-4.8d2*u80)+u6)
  u75=+u411*(3.59280826713175363d6*u420+2.0d0*(u417-4.8d2*u7))
  u57=u406*u56*u424
  u56=-1.41047395886939072d-1*u57
  u58=(u62+u290*(u290*(u75+u83)+u74))
  u425=u56*u290
  u52=u290*(u73+u290*(u425+u59))+u429
  value=value+(c(60)*(u58*u356+u52))
  value=value+(c(61)*((u290*(u74+u290*(u83+u75))+u62)*u68))
  value=value+(c(62)*(u58*u412+u52))
  u74=-2.3625d2*u95
  u59=u72*(2.0d0*(u81-1.6d2*u80)+9.21232889008141955d4*u420)
  u58=u50*(5.0d0*u54+u81)
  u52=8.46284375321634431d-1*u58
  u75=+u70*(u98+5.0d0*u430)
  u429=-4.23142187660817215d-1*u57
  u45=(u59+u290*(u83+u75))
  u77=u45*u356
  u426=u429*u290
  u95=u290*(u52+u426)+u74
  value=value+(c(63)*(u90*(u77+u95)))
  u53=2.82094791773878144d-1*u58
  u58=u290*(u53+u425)+u62
  value=value+(c(64)*(x*(u77+u58)*z))
  u77=u45*u412
  value=value+(c(65)*(u90*(u77+u58)))
  value=value+(c(66)*(u46*(u77+u95)))
  u62=1.17247822237399885d4*u420
  u74=1.12d2*u87
  u73=2.01496279838484388d-2*u50*(2.0d0*(u74-3.0d0*u88)+u62)
  u59=2.7d1*u420
  u88=u50*(u59+u94)
  u75=-1.575d2*u88
  u77=2.0d0*(u81-9.6d1*u80)
  u45=5.52739733404885173d4*u420
  u95=u50*(u77+u45)
  u58=1.41047395886939072d-1*u95
  u87=-1.575d2*u50*(u59+u89)
  u59=u50*(3.0d0*u54+u81)
  u54=1.69256875064326886d0*u59
  u56=3.0d0*u430
  u430=+u70*(u98+u56)
  u53=u406*u49*u424
  u425=4.23142187660817215d-1*u53
  u52=-8.46284375321634431d-1*u57
  u49=(u290*(u430+u83)+u58)
  u427=u52*u290
  u96=u75+u290*(u427+u54)
  u428=u425*u290
  u78=u290*(u87+u428)+u73
  value=value+(c(67)*(u356*(u96+u49*u356)+u78))
  u65=-7.875d1*u88
  u55=8.46284375321634431d-1*u59
  u61=(u58+u290*(u83+u430))
  u8=u290*(u55+u426)+u65
  value=value+(c(68)*(y*(u61*u356+u8)*z))
  u426=(1.67496888910571265d3*u5-u76)*u431
  u76=-2.68661706451312518d-2*u50*(3.0d0*a3*(u43-u419)+1.4d1*u426)
  u43=+u411*(u45+u77)
  u77=u69*(2.0d0*(u81-1.92d2*u80)+1.05523040013659897d5*u420)
  u419=u416*(2.7d1+4.0d0*pd2)
  u45=+u411*(1.49239728019318997d6*u420+8.0d0*(a3*(-u67*(pd2+6.0d0)-u91&
 &)-2.4d1*u419*u431))
  u59=2.82094791773878144d-1*u50
  u44=u59*(u56+u417)
  u91=(1.65821920021465552d5*u5+u63*(9.0d0+pd2))*u431
  u56=+u70*(u81+u91)
  u67=u416*(2.1d1*u423+u85)
  u42=u69*(2.0d0*(u98-3.2d1*u67*u431)+u97)
  u97=-1.41047395886939072d-1*u93
  u98=1.41047395886939072d-1*u57
  u69=u97*u290
  value=value+(c(69)*(u356*(u58+u290*(u290*(u42+u69)+u45)+(u290*(u44+u6&
 &9)+u43)*u356)+u290*(u77+u290*(u98*u290+u56))+u76))
  value=value+(c(70)*(u68*(u61*u412+u8)))
  value=value+(c(71)*(u412*(u96+u49*u412)+u78))
  u76=-3.9375d2*u88
  u43=1.41047395886939072d0*u95
  u77=+u411*(7.18561653426350725d5*u420+2.0d0*(u417-9.6d1*u7))
  u58=2.11571093830408608d0*u53
  u88=-1.41047395886939072d0*u57
  u56=(u83+u77)
  u61=u56*u356
  u429=u88*u290
  u44=u43+u429
  u430=u58*u290
  u45=u430+u76
  value=value+(c(72)*(u90*(u356*(u44+u61)+u45)))
  u55=8.46284375321634431d-1*u95
  u75=u55+u427
  u87=u428+u65
  value=value+(c(73)*(x*(u356*(u75+u61)+u87)*z))
  u427=(8.37484444552856323d3*u5+u66)*u431
  u428=2.53885312596490329d0*u50
  u61=u428*(u427+6.2035884781693061d1*(-u89))
  u66=u84*u48
  u96=-8.46284375321634431d-1*u66
  u48=a3*u82
  u82=u416*(1.5d1+u404)
  u54=u82*u431
  u8=+u411*(8.2910960010732776d5*u420+2.0d0*(u48-9.6d1*u54))
  u416=u92*u47
  u47=4.23142187660817215d-1*u57
  u92=-4.23142187660817215d-1*u406
  u42=+u92*(2.76369866702442587d5*u5-6.4d1*u82)*u424
  u82=(1.19760275571058454d6*u5+u63*(5.0d0*u423+u86))
  u78=u418*u82*u424
  u418=u8+u290*(u69+u78)
  u404=u290*(u42+u47*u290)+u61
  u63=(u69+u47)
  u423=u416*u356
  value=value+(c(74)*(u90*((u96+u423)*u412+u356*(u418+u63*u356)+u404)))
  value=value+(c(75)*(u46*(u412*(u418+u423+u63*u412)+u96*u356+u404)))
  u42=u56*u412
  value=value+(c(76)*(u90*(u412*(u75+u42)+u87)))
  value=value+(c(77)*(u46*(u412*(u44+u42)+u45)))
  u42=u71*(2.0d0*(u74+u51)+u62)
  u96=u50*(9.0d0*u420+u94)
  u77=-1.18125d3*u96
  u44=-3.2d1*u80
  u43=u50*(2.0d0*(u81+u44)+1.84246577801628391d4*u420)
  u80=2.11571093830408608d0*u43
  u81=+u411*(2.39520551142116908d5*u420+2.0d0*(u417-3.2d1*u7))
  u7=-3.54375d3*u406*u5*u424
  u411=6.34713281491225823d0*u53
  u417=-2.11571093830408608d0*u57
  u94=(u81+u83)
  u51=u77+u411*u290
  u62=u417*u290+u80
  u404=u7*u290+u42
  value=value+(c(78)*(u356*(u51+u356*(u94*u356+u62))+u404))
  u78=-3.9375d2*u96
  u96=1.41047395886939072d0*u43
  u43=(u83+u81)
  u83=u96+u429
  u406=u430+u78
  value=value+(c(79)*(y*(u356*(u83+u43*u356)+u406)*z))
  u56=-1.34330853225656259d-2*u50*(9.0d0*a3*(u64+u60)+2.8d1*u426)
  u426=1.575d2*u84*u79
  u79=1.86107654345079183d2*(-u89)
  u45=u72*(4.0d0*(u79+u44)+3.51743466712199656d4*u420)
  u420=-1.69256875064326886d0*u66
  u44=+u70*(u48+u91)
  u70=+u92*(4.97465760064396656d5*u5-6.4d1*u419)*u424
  u419=u9*(5.02993157398445508d6*u5-6.4d1*u67)*u424
  u424=-4.23142187660817215d-1*u53
  u53=8.46284375321634431d-1*u57
  u57=u45+u290*(u53*u290+u70)
  u5=u290*(u419+u69)+u44
  u67=u290*(u425+u424*u290)+u56
  u92=(u98+u69)
  value=value+(c(80)*((u426+u356*(u423+u420))*u412+u356*(u57+u356*(u92*&
 &u356+u5))+u67))
  u69=u428*(u427+u79)
  u9=-4.23142187660817215d-1*u50*(u6+2.0d0*(u48-3.2d1*u54))
  u48=u59*(2.0d0*u82*u431+a3*u385)
  u431=u97*u356
  value=value+(c(81)*(u68*(u412*(u9+u356*(u431+u48)+(u431+u47)*u412)+u3&
 &56*(u9+u47*u356)+u69)))
  value=value+(c(82)*(u412*(u57+u420*u356+u412*(u92*u412+u423+u5))+u426&
 &*u356+u67))
  value=value+(c(83)*(u68*(u412*(u83+u43*u412)+u406)))
  value=value+(c(84)*(u412*(u51+u412*(u94*u412+u62))+u404))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_R_opt_2_1=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_R_opt_2_2(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_R_opt_2_2
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
  D_X_R_opt_2_2=.true.
  
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
  u99=7.5d-1*u6*A3_(p,pd,erfpd,exppd2)*p**3
  u126=x*y
  value=value+(c(1)*(u99*u126))
  if ( lmax .eq. 0 ) return
  u99=A4_(p,pd,erfpd,exppd2)
  u177=p**4
  u77=-7.5d-1*u6*u99*u177
  u31=1.41047395886939072d-1*u6
  u71=exppd2*pd
  u69=-8.0d0*u71
  u52=u31*(2.65868077635827404d1*u99+u69)*u177
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
  u78=-1.125d1*u90
  u52=exppd2*pd2
  u60=-1.6d1*u52
  u91=u31*(1.86107654345079183d2*u97+u60)*u74
  u98=u91*u177
  value=value+(c(5)*(x*(u98+u78)*y))
  u59=8.0d0*exppd2
  u88=2.82094791773878144d-2*u6*(2.65868077635827404d1*u97+u59)*u74
  u74=-3.75d0*u90
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
  u83=1.125d1*u91
  u97=-1.6d1*u71
  u68=u6*(1.86107654345079183d2*u98+u97)*u88
  u64=-8.46284375321634431d-1*u68
  u74=2.0d0*pd2
  u67=u31*(1.67496888910571265d3*u98+u97*(9.0d0+u74))*u88
  u88=u67*u177
  value=value+(c(11)*((u177*(u64+u88)+u83)*y))
  u97=-4.23142187660817215d-1*u68
  u65=-1.41047395886939072d-1*u68
  u85=u65*u177
  value=value+(c(12)*(x*((u97+u88)*u99+u85+u83)))
  value=value+(c(13)*(x*(u88+u97)*u90))
  u95=(u65+u88)
  u88=u95*u99
  value=value+(c(14)*(y*(u88+u97*u177+u83)))
  u84=3.75d0*u91
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
  u85=3.9375d2*u88
  u84=-3.2d1*u52
  u67=u6*(1.67496888910571265d3*u95+u84)*u92
  u68=-1.41047395886939072d0*u67
  u97=(1.1d1+u74)
  u64=u6*(1.84246577801628391d4*u95+u84*u97)*u92
  u65=1.41047395886939072d-1*u64
  u96=u65*u177
  value=value+(c(21)*(x*(u177*(u68+u96)+u85)*y))
  u58=-6.04488839515453165d-2*u6*(1.6d1*exppd2+1.86107654345079183d2*u9&
 &5)*u92
  u82=7.875d1*u88
  u86=1.575d2*u88
  u8=-8.46284375321634431d-1*u67
  u7=-1.41047395886939072d-1*u67
  u61=u7*u177
  value=value+(c(22)*((u82+u177*(u96+u8))*u99+u177*(u86+u61)+u58))
  value=value+(c(23)*((u177*(u8+u96)+u82)*u90))
  u75=2.3625d2*u88
  u81=-4.23142187660817215d-1*u67
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
  u75=4.02992559676968776d-2*u6*(5.58322963035237549d2*u95-u59)*u92
  u59=1.41047395886939072d-1*u67
  u67=5.64189583547756287d-1*u6
  u96=u67*(5.02490666731713794d3*u95+u60*(6.0d0+pd2))*u92
  u66=-1.41047395886939072d-1*u64
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
  u75=-3.9375d2*u89
  u88=-3.2d1*u71
  u66=u6*(1.67496888910571265d3*u96+u88)*u70
  u60=6.34713281491225823d0*u66
  u92=1.84246577801628391d4*u96
  u65=u6*(u92+u88*u97)*u70
  u97=-2.11571093830408608d0*u65
  u82=pd2**2
  u85=4.0d0*u82
  u64=(u74+1.1d1)
  u7=u6*(2.39520551142116908d5*u96+u88*(1.3d1*u64+u85))*u70
  u81=1.41047395886939072d-1*u7
  u8=u81*u177
  value=value+(c(36)*((u177*(u60+u177*(u8+u97))+u75)*y))
  u68=2.11571093830408608d0*u66
  u61=1.41047395886939072d0*u66
  u98=-1.41047395886939072d0*u65
  u94=-1.41047395886939072d-1*u65
  u63=u94*u177
  value=value+(c(37)*(x*((u68+u177*(u8+u98))*u99+u177*(u61+u63)+u75)))
  value=value+(c(38)*(x*(u177*(u98+u8)+u68)*u90))
  u79=-2.3625d2*u89
  u55=4.23142187660817215d-1*u66
  u62=2.53885312596490329d0*u66
  u87=-8.46284375321634431d-1*u65
  u9=-4.23142187660817215d-1*u65
  u73=(u55+u177*(u8+u87))
  u76=u73*u99
  u54=u9*u177
  value=value+(c(39)*(y*(u76+u177*(u62+u54)+u79)))
  u72=-7.875d1*u89
  u86=8.46284375321634431d-1*u66
  u53=u177*(u86+u63)+u72
  value=value+(c(40)*((u76+u53)*z))
  value=value+(c(41)*(y*(u73*u78+u53)))
  u76=u87*u177
  u73=u55*u177
  value=value+(c(42)*(x*(u99*(u62+u76+(u8+u9)*u99)+u73+u79)))
  u53=1.26942656298245165d0*u66
  u63=(u9+u8)
  u80=u54+u53
  value=value+(c(43)*(u126*(u63*u99+u80)*z))
  u54=3.38513750128653772d0*u6*(1.86107654345079183d2*u96-4.0d0*u71)*u7&
 &0
  u66=4.23142187660817215d-1*u65
  u58=(1.3d1+u74)
  u71=-1.41047395886939072d-1*u6*(2.17745955583742644d4*u96+u88*u58)*u7&
 &0
  u88=2.25675833419102515d0*u6*(u92+u69*(4.0d0*u64+u82))*u70
  u70=-1.41047395886939072d-1*u7
  u64=1.41047395886939072d-1*u65
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
  u80=-2.480625d4*u87
  u92=-6.4d1*u52
  u65=u6*(1.84246577801628391d4*u94+u92)*u77
  u63=1.48099765681286025d1*u65
  u66=u6*(2.39520551142116908d5*u94+u92*u58)*u77
  u68=-2.96199531362572051d0*u66
  u88=(u74+1.3d1)
  u8=u6*(3.59280826713175363d6*u94+u92*(1.5d1*u88+u85))*u77
  u60=1.41047395886939072d-1*u8
  u72=u60*u177
  value=value+(c(57)*(x*(u177*(u63+u177*(u72+u68))+u80)*y))
  u89=1.67496888910571265d3*u94
  u70=u6*(u89+3.2d1*exppd2)*u77
  u69=2.35078993144898453d-1*u70
  u71=-3.54375d3*u87
  u76=-1.063125d4*u87
  u97=6.34713281491225823d0*u65
  u58=2.11571093830408608d0*u65
  u9=-2.11571093830408608d0*u66
  u98=-1.41047395886939072d-1*u66
  u91=u98*u177
  value=value+(c(58)*((u71+u177*(u177*(u9+u72)+u97))*u99+u177*(u76+u177&
 &*(u91+u58))+u69))
  value=value+(c(59)*((u177*(u97+u177*(u72+u9))+u71)*u90))
  u64=4.23142187660817215d0*u65
  u53=-1.41047395886939072d0*u66
  u86=-4.23142187660817215d-1*u66
  u81=(u58+u177*(u72+u53))
  u73=u81*u99
  u55=u86*u177
  value=value+(c(60)*(u126*(u73+u177*(u64+u55)+u76)))
  u83=1.41047395886939072d0*u65
  u95=u177*(u83+u91)+u71
  value=value+(c(61)*(x*(u73+u95)*z))
  value=value+(c(62)*(u126*(u81*u78+u95)))
  u73=1.41047395886939072d-1*u70
  u81=-4.2525d3*u87
  u95=4.23142187660817215d-1*u65
  u91=5.07770625192980658d0*u65
  u70=-8.46284375321634431d-1*u66
  u61=u70*u177
  u62=u95*u177
  value=value+(c(63)*(u99*(u81+u177*(u61+u91)+(u177*(u70+u72)+u95)*u99)&
 &+u177*(u81+u62)+u73))
  u73=-2.12625d3*u87
  u91=2.53885312596490329d0*u65
  u87=(u95+u177*(u72+u70))
  u57=u177*(u91+u55)+u73
  value=value+(c(64)*(y*(u87*u99+u57)*z))
  u55=-3.76126389031837525d-1*u6*(u89-4.0d0*exppd2)*u77
  u89=-4.23142187660817215d-1*u65
  u56=4.23142187660817215d-1*u6*(3.51743466712199656d4*u94-1.28d2*u52)*&
 &u77
  u54=-4.23142187660817215d-1*u6
  u75=+u54*(4.97465760064396656d5*u94+u92*(2.7d1+4.0d0*pd2))*u77
  u7=8.46284375321634431d-1*u66
  u96=-2.82094791773878144d-1*u6*(1.65821920021465552d5*u94+u92*(9.0d0+&
 &pd2))*u77
  u79=u31*(5.02993157398445508d6*u94+u92*(2.1d1*u88+u85))*u77
  u85=-1.41047395886939072d-1*u8
  u31=1.41047395886939072d-1*u66
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
  u57=2.53885312596490329d0*u6*(8.37484444552856323d3*u94+u84)*u77
  u84=+u54*(2.76369866702442587d5*u94+u92*(1.5d1+u74))*u77
  u54=4.23142187660817215d-1*u66
  u6=u67*(1.19760275571058454d6*u94+u92*(5.0d0*u88+u82))*u77
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
  D_X_R_opt_2_2=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_R_opt_2_3(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_R_opt_2_3
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
  D_X_R_opt_2_3=.true.
  
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
  u99=7.5d-1*u6*A3_(p,pd,erfpd,exppd2)*p**3
  u127=x*z
  value=value+(c(1)*(u99*u127))
  if ( lmax .eq. 0 ) return
  u99=A4_(p,pd,erfpd,exppd2)
  u126=p**4
  u77=-7.5d-1*u6*u99*u126
  u57=1.41047395886939072d-1*u6
  u73=exppd2*pd
  u93=-8.0d0*u73
  u97=u57*(2.65868077635827404d1*u99+u93)*u126
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
  u78=-1.125d1*u90
  u77=exppd2*pd2
  u63=-1.6d1*u77
  u83=u57*(1.86107654345079183d2*u97+u63)*u98
  u95=u83*u178
  value=value+(c(5)*(x*(u95+u78)*z))
  u74=-3.75d0*u90
  u90=y*z
  value=value+(c(6)*((u95+u74)*u90))
  u59=8.0d0*exppd2
  u91=2.82094791773878144d-2*u6*(2.65868077635827404d1*u97+u59)*u98
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
  u83=1.125d1*u91
  u8=-1.6d1*u73
  u68=u6*(1.86107654345079183d2*u98+u8)*u95
  u80=-8.46284375321634431d-1*u68
  u78=2.0d0*pd2
  u85=u57*(1.67496888910571265d3*u98+u8*(9.0d0+u78))*u95
  u95=u85*u178
  value=value+(c(11)*((u178*(u80+u95)+u83)*z))
  u8=-4.23142187660817215d-1*u68
  value=value+(c(12)*(x*(u95+u8)*u90))
  u67=-1.41047395886939072d-1*u68
  u60=u67*u178
  value=value+(c(13)*(x*((u8+u95)*u99+u60+u83)))
  u84=3.75d0*u91
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
  u85=3.9375d2*u88
  u83=-3.2d1*u77
  u67=u6*(1.67496888910571265d3*u95+u83)*u60
  u8=-1.41047395886939072d0*u67
  u84=(1.1d1+u78)
  u80=u6*(1.84246577801628391d4*u95+u83*u84)*u60
  u92=1.41047395886939072d-1*u80
  u89=u92*u178
  value=value+(c(21)*(x*(u178*(u8+u89)+u85)*z))
  u82=7.875d1*u88
  u55=-8.46284375321634431d-1*u67
  value=value+(c(22)*((u178*(u55+u89)+u82)*u90))
  u58=-6.04488839515453165d-2*u6*(1.6d1*exppd2+1.86107654345079183d2*u9&
 &5)*u60
  u86=1.575d2*u88
  u64=-1.41047395886939072d-1*u67
  u65=u64*u178
  value=value+(c(23)*((u82+u178*(u89+u55))*u99+u178*(u86+u65)+u58))
  u7=-4.23142187660817215d-1*u67
  u96=(u7+u89)
  u75=u65+u82
  value=value+(c(24)*(x*(u96*u97+u75)*z))
  u65=u96*u99
  value=value+(c(25)*(u126*(u65+u75)))
  u96=2.3625d2*u88
  u75=u7*u178
  value=value+(c(26)*(u127*(u65+u75+u96)))
  u96=(u64+u89)
  u65=u75+u82
  value=value+(c(27)*(y*(u96*u97+u65)*z))
  u75=4.02992559676968776d-2*u6*(5.58322963035237549d2*u95-u59)*u60
  u59=1.41047395886939072d-1*u67
  u67=5.64189583547756287d-1*u6
  u66=u67*(5.02490666731713794d3*u95+u63*(6.0d0+pd2))*u60
  u60=-1.41047395886939072d-1*u80
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
  u75=-3.9375d2*u89
  u63=-3.2d1*u73
  u66=u6*(1.67496888910571265d3*u96+u63)*u65
  u60=6.34713281491225823d0*u66
  u92=1.84246577801628391d4*u96
  u80=u6*(u92+u63*u84)*u65
  u59=-2.11571093830408608d0*u80
  u58=pd2**2
  u85=4.0d0*u58
  u64=(u78+1.1d1)
  u7=u6*(2.39520551142116908d5*u96+u63*(1.3d1*u64+u85))*u65
  u55=1.41047395886939072d-1*u7
  u8=u55*u178
  value=value+(c(36)*((u178*(u60+u178*(u8+u59))+u75)*z))
  u84=2.11571093830408608d0*u66
  u71=-1.41047395886939072d0*u80
  value=value+(c(37)*(x*(u178*(u71+u8)+u84)*u90))
  u61=1.41047395886939072d0*u66
  u69=-1.41047395886939072d-1*u80
  u53=u69*u178
  value=value+(c(38)*(x*((u84+u178*(u8+u71))*u99+u178*(u61+u53)+u75)))
  u72=-7.875d1*u89
  u86=4.23142187660817215d-1*u66
  u95=8.46284375321634431d-1*u66
  u94=-8.46284375321634431d-1*u80
  u79=(u86+u178*(u8+u94))
  u62=u178*(u95+u53)+u72
  value=value+(c(39)*((u79*u97+u62)*z))
  u53=u79*u99
  value=value+(c(40)*(y*(u53+u62)))
  u79=-2.3625d2*u89
  u62=2.53885312596490329d0*u66
  u9=-4.23142187660817215d-1*u80
  u54=u9*u178
  value=value+(c(41)*(z*(u53+u178*(u62+u54)+u79)))
  u53=1.26942656298245165d0*u66
  u76=(u9+u8)
  u66=u54+u53
  value=value+(c(42)*(u126*(u76*u97+u66)*z))
  u54=3.38513750128653772d0*u6*(1.86107654345079183d2*u96-4.0d0*u73)*u6&
 &5
  u88=4.23142187660817215d-1*u80
  u73=(1.3d1+u78)
  u70=-1.41047395886939072d-1*u6*(2.17745955583742644d4*u96+u63*u73)*u6&
 &5
  u96=2.25675833419102515d0*u6*(u92+u93*(4.0d0*u64+u58))*u65
  u65=-1.41047395886939072d-1*u7
  u64=1.41047395886939072d-1*u80
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
  u80=-2.480625d4*u87
  u66=-6.4d1*u77
  u65=u6*(1.84246577801628391d4*u94+u66)*u72
  u63=1.48099765681286025d1*u65
  u96=u6*(2.39520551142116908d5*u94+u66*u73)*u72
  u62=-2.96199531362572051d0*u96
  u54=(u78+1.3d1)
  u8=u6*(3.59280826713175363d6*u94+u66*(1.5d1*u54+u85))*u72
  u9=1.41047395886939072d-1*u8
  u88=u9*u178
  value=value+(c(57)*(x*(u178*(u63+u178*(u88+u62))+u80)*z))
  u71=-3.54375d3*u87
  u74=6.34713281491225823d0*u65
  u95=-2.11571093830408608d0*u96
  value=value+(c(58)*((u178*(u74+u178*(u88+u95))+u71)*u90))
  u93=1.67496888910571265d3*u94
  u70=u6*(u93+3.2d1*exppd2)*u72
  u69=2.35078993144898453d-1*u70
  u76=-1.063125d4*u87
  u59=2.11571093830408608d0*u65
  u60=-1.41047395886939072d-1*u96
  u73=u60*u178
  value=value+(c(59)*((u71+u178*(u178*(u95+u88)+u74))*u99+u178*(u76+u17&
 &8*(u73+u59))+u69))
  u84=1.41047395886939072d0*u65
  u75=-1.41047395886939072d0*u96
  u64=(u59+u178*(u88+u75))
  u53=u178*(u84+u73)+u71
  value=value+(c(60)*(x*(u64*u97+u53)*z))
  u73=u64*u99
  value=value+(c(61)*(u126*(u73+u53)))
  u64=4.23142187660817215d0*u65
  u53=-4.23142187660817215d-1*u96
  u55=u53*u178
  value=value+(c(62)*(u127*(u73+u178*(u64+u55)+u76)))
  u73=-2.12625d3*u87
  u86=4.23142187660817215d-1*u65
  u91=2.53885312596490329d0*u65
  u98=-8.46284375321634431d-1*u96
  u61=(u86+u178*(u88+u98))
  u81=u178*(u91+u55)+u73
  value=value+(c(63)*(y*(u61*u97+u81)*z))
  u55=-3.76126389031837525d-1*u6*(u93-4.0d0*exppd2)*u72
  u82=-4.23142187660817215d-1*u65
  u56=4.23142187660817215d-1*u6*(3.51743466712199656d4*u94-1.28d2*u77)*&
 &u72
  u68=-4.23142187660817215d-1*u6
  u93=+u68*(4.97465760064396656d5*u94+u66*(2.7d1+4.0d0*pd2))*u72
  u79=8.46284375321634431d-1*u96
  u7=-2.82094791773878144d-1*u6*(1.65821920021465552d5*u94+u66*(9.0d0+p&
 &d2))*u72
  u89=u57*(5.02993157398445508d6*u94+u66*(2.1d1*u54+u85))*u72
  u85=-1.41047395886939072d-1*u8
  u92=1.41047395886939072d-1*u96
  u77=u85*u178
  u8=u178*(u89+u77)
  value=value+(c(64)*(u97*(u86+u178*(u8+u93)+(u178*(u79+u77)+u82)*u97)+&
 &u178*(u56+u178*(u92*u178+u7))+u55))
  value=value+(c(65)*(u90*(u61*u99+u81)))
  u61=1.41047395886939072d-1*u70
  u81=-4.2525d3*u87
  u87=5.07770625192980658d0*u65
  u65=u98*u178
  u57=u86*u178
  value=value+(c(66)*(u99*(u81+u178*(u65+u87)+(u178*(u98+u88)+u86)*u99)&
 &+u178*(u81+u57)+u61))
  u81=(u88+u53)
  u87=u91+u65
  u61=u57+u73
  value=value+(c(67)*(x*(u97*(u87+u81*u97)+u61)*z))
  u57=2.53885312596490329d0*u6*(8.37484444552856323d3*u94+u83)*u72
  u83=+u68*(2.76369866702442587d5*u94+u66*(1.5d1+u78))*u72
  u68=4.23142187660817215d-1*u96
  u6=u67*(1.19760275571058454d6*u94+u66*(5.0d0*u54+u58))*u72
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
  D_X_R_opt_2_3=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_R_opt_2_4(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_R_opt_2_4
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u296
  real(kind=8) :: u347
  real(kind=8) :: u348
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
  real(kind=8) :: u381
  real(kind=8) :: u382
  real(kind=8) :: u386
  real(kind=8) :: u387
  real(kind=8) :: u388
  real(kind=8) :: u389
  real(kind=8) :: u39
  real(kind=8) :: u390
  real(kind=8) :: u391
  real(kind=8) :: u392
  real(kind=8) :: u393
  real(kind=8) :: u394
  real(kind=8) :: u395
  real(kind=8) :: u4
  real(kind=8) :: u40
  real(kind=8) :: u41
  real(kind=8) :: u42
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
  D_X_R_opt_2_4=.true.
  
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
  u395=A3_(p,pd,erfpd,exppd2)
  u361=p**3
  u296=u395*u361
  u391=a3**(-2)
  u39=-4.70157986289796906d-2*u391*(5.31736155271654808d0*u296+2.0d0*(2&
 &.0d0*exppd2*u361-5.31736155271654808d0*a3*A1_(p,pd,erfpd,exppd2)))
  u78=7.5d-1*u391*u296
  u296=y**2
  value=value+(c(1)*(u78*u296+u39))
  if ( lmax .eq. 0 ) return
  u39=A4_(p,pd,erfpd,exppd2)
  u78=p**2
  u90=u39*u78
  u48=-2.0d0*a3
  u89=+u48*A2_(p,pd,erfpd,exppd2)
  u59=-2.5d-1*u391*u78
  u67=+u59*(3.0d0*u90+u89)
  u386=2.65868077635827404d1*u39
  u97=1.41047395886939072d-1*u391
  u64=p**4
  u382=exppd2*pd
  u98=-8.0d0*u382
  u99=u97*(u386+u98)*u64
  u363=u99*u296
  u347=(u363+u67)
  value=value+(c(2)*(x*u347))
  u362=+u59*(9.0d0*u90+u89)
  value=value+(c(3)*(y*(u363+u362)))
  value=value+(c(4)*(u347*z))
  if ( lmax .eq. 1 ) return
  u362=4.0d0*exppd2
  u99=-5.0d0*a3*(u362+5.31736155271654808d0*u395)
  u90=A5_(p,pd,erfpd,exppd2)
  u89=u90*u78
  u59=u391*u361
  u361=9.40315972579593812d-3*u59
  u67=exppd2*u78
  u38=u361*(2.0d0*(1.2d1*u67+u99)+7.97604232907482212d1*u89)
  u93=p**5
  u370=u90*u93
  u88=-3.75d0*u391*u370
  u81=a3*u395
  u395=-2.0d0*u81
  u71=-7.5d-1*u59
  u66=+u71*(5.0d0*u89+u395)
  u55=1.86107654345079183d2*u90
  u392=exppd2*pd2
  u363=-1.6d1*u392
  u96=(u55+u363)*u93
  u377=u97*u96
  u347=x**2
  value=value+(c(5)*((u88+u377*u347)*u296+u66*u347+u38))
  u9=1.5d1*u89
  u68=+u71*(u9+u395)
  u395=x*y
  u71=u377*u296
  u348=(u71+u68)
  value=value+(c(6)*(u395*u348))
  value=value+(c(7)*(x*(u71+u66)*z))
  u40=u361*(2.0d0*(3.6d1*u67+u99)+2.39281269872244664d2*u89)
  u99=-1.5d0*u59*(u9-u81)
  value=value+(c(8)*(u296*(u99+u71)+u40))
  value=value+(c(9)*(y*u348*z))
  u348=z**2
  value=value+(c(10)*((u66+u71)*u348+u88*u296+u38))
  if ( lmax .eq. 2 ) return
  u99=A6_(p,pd,erfpd,exppd2)
  u377=u99*u78
  u88=+u48*u39
  u9=u391*u64
  u89=u9*(5.0d0*u377+u88)
  u81=2.25d0*u89
  u39=1.86107654345079183d2*u99
  u361=-1.6d1*u382
  u59=(u39+u361)
  u40=p**6
  u64=u59*u40
  u48=u391*u64
  u71=-4.23142187660817215d-1*u48
  u57=a3*(-u98-u386)
  u386=-1.41047395886939072d-1*u9
  u68=u382*u78
  u390=+u386*(1.86107654345079183d2*u377+2.0d0*(u57-8.0d0*u68))
  u51=1.67496888910571265d3*u99
  u38=2.0d0*pd2
  u58=+u361*(9.0d0+u38)
  u6=(u51+u58)
  u380=u6*u40
  u371=u97*u380
  u364=u371*u347
  u365=u390*u347
  value=value+(c(11)*(x*((u71+u364)*u296+u365+u81)))
  u80=7.5d-1*u9*(1.5d1*u377+u88)
  u36=-1.41047395886939072d-1*u48
  u48=+u386*(5.58322963035237549d2*u377+2.0d0*(u57-2.4d1*u68))
  u386=(u36+u364)*u296
  value=value+(c(12)*(y*(u386+u48*u347+u80)))
  u82=7.5d-1*u89
  value=value+(c(13)*((u386+u365+u82)*z))
  u386=u59*u78
  u83=-2.82094791773878144d-1*u9
  u47=+u83*(u57+3.0d0*u386)
  u388=u371*u296
  u365=(u296*(u47+u388)+u80)
  value=value+(c(14)*(x*u365))
  value=value+(c(15)*(u395*(u388+u48)*z))
  u364=(u390+u388)*u348
  u366=u36*u296
  value=value+(c(16)*(x*(u364+u366+u82)))
  u66=2.25d0*u9*(2.5d1*u377+u88)
  u377=+u83*(u57+5.0d0*u386)
  value=value+(c(17)*(y*(u296*(u377+u388)+u66)))
  value=value+(c(18)*(u365*z))
  value=value+(c(19)*(y*((u48+u388)*u348+u366+u80)))
  value=value+(c(20)*(z*(u364+u71*u296+u81)))
  if ( lmax .eq. 3 ) return
  u66=8.0d0*exppd2
  u377=-7.0d0*a3*(u66+2.65868077635827404d1*u90)
  u9=A7_(p,pd,erfpd,exppd2)
  u36=u9*u78
  u71=u391*u93
  u390=-1.20897767903090633d-2*u71
  u386=+u390*(9.30538271725395914d2*u36+2.0d0*(4.0d1*u67+u377))
  u365=p**7
  u388=u9*u365
  u57=7.875d1*u391*u388
  u88=a3*u90
  u371=-2.0d0*u88
  u90=u71*(7.0d0*u36+u371)
  u83=2.25d1*u90
  u80=1.67496888910571265d3*u9
  u364=-3.2d1*u392
  u81=(u80+u364)
  u47=u81*u365
  u48=u391*u47
  u93=-8.46284375321634431d-1*u48
  u379=1.67496888910571265d3*u36
  u8=a3*(-u363-u55)
  u55=-1.41047395886939072d-1*u71
  u82=u392*u78
  u394=-1.6d1*u82
  u69=+u55*(u379+2.0d0*(u8+u394))
  u393=1.84246577801628391d4*u9
  u95=(1.1d1+u38)
  u61=+u364*u95
  u89=(u393+u61)
  u366=u89*u365
  u369=u391*u366
  u5=1.41047395886939072d-1*u369
  u360=u5*u347
  u367=u69*u347
  value=value+(c(21)*((u57+u347*(u360+u93))*u296+u347*(u83+u367)+u386))
  u52=2.1d1*u36
  u91=u71*(u52+u371)
  u84=1.125d1*u91
  u372=-4.23142187660817215d-1*u48
  u44=+u55*(5.02490666731713794d3*u36+2.0d0*(u8-4.8d1*u82))
  u87=(u372+u360)*u296
  u368=u44*u347
  value=value+(c(22)*(u395*(u87+u368+u84)))
  u85=1.125d1*u90
  value=value+(c(23)*(x*(u87+u367+u85)*z))
  u87=-4.02992559676968777d-3*u71*(2.79161481517618774d3*u36+2.0d0*(1.2&
 &d2*u67+u377))
  u92=7.5d0*u71*(u52-u88)
  u52=-1.41047395886939072d-1*u48
  u79=3.75d0*u91
  u381=u81*u78
  u389=-2.82094791773878144d-1*u71
  u46=+u389*(u8+3.0d0*u381)
  value=value+(c(24)*(u296*(u92+u46*u347+(u360+u52)*u296)+u79*u347+u87)&
 &)
  value=value+(c(25)*(y*((u52+u360)*u296+u368+u79)*z))
  u360=5.58322963035237549d2*u9
  u56=8.05985119353937553d-3*u71*(5.0d0*(u360-u66)*u78-u377)
  u367=a3**(-1)
  u86=-7.5d0*u367*u370
  u370=1.41047395886939072d-1*u48
  u48=+u55*(u379+2.0d0*(2.65868077635827404d1*u88+u394))
  u91=2.82094791773878144d-1*u367
  u394=u91*u96
  u368=5.02490666731713794d3*u9
  u90=5.64189583547756287d-1*u391
  u379=u90*(u368+u363*(6.0d0+pd2))*u365
  u96=-1.41047395886939072d-1*u369
  u369=u96*u347
  value=value+(c(26)*((u86+u394*u347)*u348+u296*(u52+u347*(u369+u379)+(&
 &u369+u370)*u296)+u347*(u48+u370*u347)+u56))
  u96=1.125d1*u71
  u86=u96*(3.5d1*u36+u371)
  u48=+u389*(u8+5.0d0*u381)
  u56=u5*u296
  u381=(u296*(u48+u56)+u86)
  value=value+(c(27)*(u395*u381))
  value=value+(c(28)*(x*(u296*(u46+u56)+u79)*z))
  u379=(u44+u56)*u348
  u370=u52*u296
  value=value+(c(29)*(u395*(u379+u370+u79)))
  u369=x*z
  u371=u372*u296
  value=value+(c(30)*(u369*((u69+u56)*u348+u371+u85)))
  u71=+u390*(4.65269135862697957d3*u36+2.0d0*(2.0d2*u67+u377))
  u377=u96*(1.05d2*u36-4.0d0*u88)
  u85=-2.4d2*u82
  u390=+u55*(2.51245333365856897d4*u36+2.0d0*(u8+u85))
  value=value+(c(31)*(u296*(u377+u296*(u56+u390))+u71))
  value=value+(c(32)*(y*u381*z))
  value=value+(c(33)*((u79+u296*(u56+u46))*u348+u296*(u92+u370)+u87))
  u370=y*z
  value=value+(c(34)*(u370*(u379+u371+u84)))
  value=value+(c(35)*(u348*(u83+u93*u296+(u56+u69)*u348)+u57*u296+u386)&
 &)
  if ( lmax .eq. 4 ) return
  u377=A8_(p,pd,erfpd,exppd2)
  u390=u377*u78
  u5=a3*u99
  u52=-2.0d0*u5
  u48=u391*u40
  u93=u48*(7.0d0*u390+u52)
  u69=-5.625d1*u93
  u44=-3.2d1*u382
  u57=(1.67496888910571265d3*u377+u44)
  u8=p**8
  u36=u391*u57*u8
  u386=2.11571093830408608d0*u36
  u87=a3*(-u361-u39)
  u39=u48*(2.0d0*(u87-1.6d1*u68)+1.67496888910571265d3*u390)
  u99=1.41047395886939072d0*u39
  u55=1.84246577801628391d4*u377
  u371=u382*u95
  u71=(u55-3.2d1*u371)
  u88=u391*u71*u8
  u40=-1.41047395886939072d0*u88
  u95=a3*(-u361*(u38+9.0d0)-u51)
  u381=1.84246577801628391d4*u390
  u86=u371*u78
  u371=-1.6d1*u86
  u84=-1.41047395886939072d-1*u48
  u96=+u84*(u381+2.0d0*(u95+u371))
  u79=pd2**2
  u83=4.0d0*u79
  u56=(u38+1.1d1)
  u379=u391*(2.39520551142116908d5*u377+u44*(1.3d1*u56+u83))*u8
  u394=1.41047395886939072d-1*u379
  u389=u394*u347
  u372=u96*u347
  value=value+(c(36)*(x*((u386+u347*(u389+u40))*u296+u347*(u99+u372)+u6&
 &9)))
  u94=u48*(2.1d1*u390+u52)
  u63=-1.125d1*u94
  u41=4.23142187660817215d-1*u36
  u36=u48*(2.0d0*(u87-4.8d1*u68)+5.02490666731713794d3*u390)
  u46=8.46284375321634431d-1*u36
  u73=-8.46284375321634431d-1*u88
  u62=5.52739733404885173d4*u390
  u65=-4.8d1*u86
  u77=+u84*(u62+2.0d0*(u95+u65))
  u42=(u41+u347*(u389+u73))*u296
  u373=u77*u347
  value=value+(c(37)*(y*(u42+u347*(u46+u373)+u63)))
  u70=-1.125d1*u93
  u76=8.46284375321634431d-1*u39
  value=value+(c(38)*((u42+u347*(u76+u372)+u70)*z))
  u42=u57*u78
  u57=u48*(3.0d0*u42+u87)
  u49=8.46284375321634431d-1*u57
  u4=-4.23142187660817215d-1*u88
  u72=1.41047395886939072d-1*u36
  u93=u71*u78
  u71=a3*(-u58-u51)
  u51=-2.82094791773878144d-1*u48
  u58=+u51*(u71+3.0d0*u93)
  u374=u58*u347
  u375=u72*u347
  value=value+(c(39)*(x*(u296*(u49+u374+(u389+u4)*u296)+u375+u63)))
  u54=4.23142187660817215d-1*u36
  value=value+(c(40)*(u395*((u4+u389)*u296+u373+u54)*z))
  u387=4.0d0*(1.86107654345079183d2*u377-4.0d0*u382)*u78
  u382=8.46284375321634431d-1*u48
  u60=u382*(u387+2.65868077635827404d1*u5)
  u36=u367*u64
  u64=-8.46284375321634431d-1*u36
  u378=4.23142187660817215d-1*u88
  u92=a3*u59
  u59=(1.3d1+u38)
  u74=u59*u78
  u7=+u84*(2.17745955583742644d4*u390+2.0d0*(u92+u361*u74))
  u39=u91*u380
  u380=(u55+u98*(4.0d0*u56+u79))
  u55=2.25675833419102515d0*u391*u380*u8
  u98=-1.41047395886939072d-1*u379
  u379=1.41047395886939072d-1*u88
  u361=u98*u347
  u373=u347*(u361+u55)
  u376=u39*u347
  value=value+(c(41)*(x*((u64+u376)*u348+u296*(u4+u373+(u361+u378)*u296&
 &)+u347*(u7+u379*u347)+u60)))
  u64=-1.125d1*u48*(3.5d1*u390+u52)
  u75=2.82094791773878144d-1*u48
  u50=u75*(5.0d0*u42+u87)
  u42=-1.41047395886939072d-1*u88
  u88=4.23142187660817215d-1*u48
  u45=u88*(2.0d0*(u87-8.0d1*u68)+8.37484444552856323d3*u390)
  u53=+u51*(u71+5.0d0*u93)
  u93=(u389+u42)*u296
  value=value+(c(42)*(y*(u296*(u50+u53*u347+u93)+u45*u347+u64)))
  u71=-3.75d0*u94
  u51=2.82094791773878144d-1*u57
  value=value+(c(43)*((u296*(u51+u374+u93)+u375+u71)*z))
  u93=u382*(u387+8.86226925452758014d0*u5)
  u372=-2.82094791773878144d-1*u36
  u382=-1.41047395886939072d-1*u391*(2.17745955583742644d4*u377+u44*u59&
 &)*u8
  u5=+u84*(u62+2.0d0*(u92+u65))
  u94=(u361+u379)
  u377=u378*u347
  value=value+(c(44)*(y*((u372+u376)*u348+u296*(u382+u373+u94*u296)+u34&
 &7*(u5+u377)+u93)))
  u62=-4.23142187660817215d-1*u48*(u381+2.0d0*(u92+u371))
  u5=u75*(8.0d0*u380*u78+a3*u6)
  value=value+(c(45)*(z*(u348*(u7+u347*(u361+u5)+u94*u348)+u347*(u62+u3&
 &77)+u60)))
  u60=u88*(4.0d0*(u87-1.2d2*u68)+2.51245333365856897d4*u390)
  u5=+u84*(2.76369866702442587d5*u390+2.0d0*(u95-2.4d2*u86))
  u62=u394*u296
  u372=(u296*(u60+u296*(u62+u5))+u64)
  value=value+(c(46)*(x*u372))
  value=value+(c(47)*(u395*(u296*(u53+u62)+u45)*z))
  u371=(u72+u296*(u62+u58))*u348
  u378=u42*u296
  value=value+(c(48)*(x*(u371+u296*(u51+u378)+u71)))
  u379=u4*u296
  value=value+(c(49)*(u395*z*((u77+u62)*u348+u379+u54)))
  u94=(u62+u96)*u348
  u380=u73*u296
  u381=u41*u296
  value=value+(c(50)*(x*(u348*(u76+u380+u94)+u381+u70)))
  u373=-5.625d1*u48*(4.9d1*u390+u52)
  u52=7.05236979434695359d-1*u48*(4.0d0*(u87-1.68d2*u68)+3.517434667121&
 &99656d4*u390)
  u87=+u84*(3.86917813383419621d5*u390+2.0d0*(u95-3.36d2*u86))
  value=value+(c(51)*(y*(u296*(u52+u296*(u62+u87))+u373)))
  value=value+(c(52)*(u372*z))
  value=value+(c(53)*(y*((u45+u296*(u62+u53))*u348+u296*(u50+u378)+u64)&
 &))
  value=value+(c(54)*(z*(u371+u296*(u49+u379)+u63)))
  value=value+(c(55)*(y*(u348*(u46+u380+(u62+u77)*u348)+u381+u63)))
  value=value+(c(56)*(z*(u348*(u99+u40*u296+u94)+u386*u296+u69)))
  if ( lmax .eq. 5 ) return
  u50=A9_(p,pd,erfpd,exppd2)
  u63=u50*u78
  u87=1.17247822237399885d4*u63
  u96=1.86107654345079183d2*u9
  u372=1.6d1*exppd2
  u58=a3*(u372+u96)
  u380=-9.0d0*u58
  u68=u391*u365
  u5=3.35827133064140647d-2*u68
  u94=1.12d2*u67
  u40=u5*(2.0d0*(u94+u380)+u87)
  u60=p**9
  u371=-3.54375d3*u391*u50*u60
  u378=a3*u9
  u49=-2.0d0*u378
  u95=u68*(9.0d0*u63+u49)
  u72=-1.18125d3*u95
  u62=-6.4d1*u392
  u52=(1.84246577801628391d4*u50+u62)
  u390=u391*u52*u60
  u373=6.34713281491225823d0*u390
  u69=1.84246577801628391d4*u63
  u84=a3*(-u364-u80)
  u45=-3.2d1*u82
  u41=u68*(2.0d0*(u84+u45)+u69)
  u377=2.11571093830408608d0*u41
  u71=(2.39520551142116908d5*u50+u62*u59)
  u80=u391*u71*u60
  u86=-2.11571093830408608d0*u80
  u93=a3*(-u364*u56-u393)
  u64=-1.41047395886939072d-1*u68
  u70=u392*u74
  u381=+u64*(2.39520551142116908d5*u63+2.0d0*(u93-3.2d1*u70))
  u77=(u38+1.3d1)
  u46=u391*(3.59280826713175363d6*u50+u62*(1.5d1*u77+u83))*u60
  u99=1.41047395886939072d-1*u46
  u389=u99*u347
  u382=u381*u347
  value=value+(c(57)*((u371+u347*(u347*(u86+u389)+u373))*u296+u347*(u72&
 &+u347*(u382+u377))+u40))
  u75=2.7d1*u63
  u88=u68*(u75+u49)
  u73=-3.9375d2*u88
  u98=2.11571093830408608d0*u390
  u365=2.0d0*(u84-9.6d1*u82)
  u57=5.52739733404885173d4*u63
  u53=u68*(u365+u57)
  u7=1.41047395886939072d0*u53
  u51=-1.41047395886939072d0*u80
  u361=+u64*(7.18561653426350725d5*u63+2.0d0*(u93-9.6d1*u70))
  u42=(u98+u347*(u389+u51))*u296
  u374=u361*u347
  value=value+(c(58)*(u395*(u42+u347*(u7+u374)+u73)))
  u74=-3.9375d2*u95
  u376=1.41047395886939072d0*u41
  value=value+(c(59)*(x*(u42+u347*(u376+u382)+u74)*z))
  u42=2.01496279838484388d-2*u68*(2.0d0*(u94-3.0d0*u58)+u87)
  u87=-1.575d2*u68*(u75-u378)
  u95=4.23142187660817215d-1*u390
  u75=-1.575d2*u88
  u94=u52*u78
  u58=u68*(3.0d0*u94+u84)
  u52=1.69256875064326886d0*u58
  u379=-8.46284375321634431d-1*u80
  u48=1.41047395886939072d-1*u53
  u41=u71*u78
  u382=3.0d0*u41
  u71=a3*(-u61-u393)
  u92=-2.82094791773878144d-1*u68
  u36=+u92*(u71+u382)
  u9=u36*u347
  u61=u48*u347
  value=value+(c(60)*(u296*(u87+u347*(u9+u52)+(u347*(u379+u389)+u95)*u2&
 &96)+u347*(u75+u61)+u42))
  u65=-7.875d1*u88
  u375=8.46284375321634431d-1*u53
  value=value+(c(61)*(y*((u95+u347*(u389+u379))*u296+u347*(u375+u374)+u&
 &65)*z))
  u393=(1.67496888910571265d3*u50-u362)*u78
  u53=-1.34330853225656259d-2*u68*(9.0d0*a3*(u96+u372)+2.8d1*u393)
  u96=1.575d2*u367*u388
  u372=-4.23142187660817215d-1*u390
  u390=1.86107654345079183d2*u378
  u374=(u390+u45)
  u39=4.23142187660817215d-1*u68
  u45=u39*(4.0d0*u374+3.51743466712199656d4*u63)
  u8=u367*u47
  u76=-1.69256875064326886d0*u8
  u56=-4.23142187660817215d-1*u391
  u391=u392*(2.7d1+4.0d0*pd2)
  u59=+u56*(4.97465760064396656d5*u50-6.4d1*u391)*u60
  u394=8.46284375321634431d-1*u80
  u47=(1.65821920021465552d5*u50+u62*(9.0d0+pd2))*u78
  u367=a3*u81
  u44=+u92*(u367+u47)
  u81=u91*u366
  u6=u392*(2.1d1*u77+u83)
  u54=u97*(5.02993157398445508d6*u50-6.4d1*u6)*u60
  u97=-1.41047395886939072d-1*u46
  u366=1.41047395886939072d-1*u80
  u46=u97*u347
  u362=u81*u347
  value=value+(c(62)*((u96+u347*(u362+u76))*u348+u296*(u95+u347*(u347*(&
 &u54+u46)+u59)+(u347*(u394+u46)+u372)*u296)+u347*(u45+u347*(u366*u347+&
 &u44))+u53))
  u96=u68*(4.5d1*u63+u49)
  u76=-2.3625d2*u96
  u59=u68*(5.0d0*u94+u84)
  u54=8.46284375321634431d-1*u59
  u88=-4.23142187660817215d-1*u80
  u91=u39*(2.0d0*(u84-1.6d2*u82)+9.21232889008141955d4*u63)
  u83=+u92*(u71+5.0d0*u41)
  u4=(u389+u88)*u296
  u386=u83*u347
  u387=u91*u347
  value=value+(c(63)*(u395*(u296*(u54+u386+u4)+u387+u76)))
  u55=8.46284375321634431d-1*u58
  value=value+(c(64)*(x*(u296*(u55+u9+u4)+u61+u65)*z))
  u4=(8.37484444552856323d3*u50+u364)*u78
  u9=2.53885312596490329d0*u68
  u61=u9*(u4+6.2035884781693061d1*u378)
  u364=-8.46284375321634431d-1*u8
  u8=u392*(1.5d1+u38)
  u58=+u56*(2.76369866702442587d5*u50-6.4d1*u8)*u60
  u38=4.23142187660817215d-1*u80
  u56=u8*u78
  u8=+u64*(8.2910960010732776d5*u63+2.0d0*(u367-9.6d1*u56))
  u392=(1.19760275571058454d6*u50+u62*(5.0d0*u77+u79))
  u50=u90*u392*u60
  u62=(u46+u38)
  u388=u38*u347
  value=value+(c(65)*(u395*((u364+u362)*u348+u296*(u58+u347*(u46+u50)+u&
 &62*u296)+u347*(u8+u388)+u61)))
  u90=u9*(u4+u390)
  u4=2.76369866702442587d5*u63
  u60=-4.23142187660817215d-1*u68
  u9=+u60*(u4+2.0d0*(u367-3.2d1*u56))
  u390=a3*u89
  u89=2.82094791773878144d-1*u68
  u77=u89*(2.0d0*u392*u78+u390)
  value=value+(c(66)*(u369*(u348*(u9+u347*(u46+u77)+u62*u348)+u347*(u9+&
 &u388)+u90)))
  u90=6.71654266128281294d-3*u68*(2.0d0*(5.6d2*u67+u380)+5.862391111869&
 &99426d4*u63)
  u77=-7.875d1*u68*(1.35d2*u63-4.0d0*u378)
  u392=1.41047395886939072d-1*u68
  u388=u392*(2.0d0*(u84-4.8d2*u82)+u4)
  u9=-1.41047395886939072d-1*u80
  u62=-7.875d1*u96
  u96=u39*(4.0d0*(u84+u85)+u4)
  u4=+u64*(3.59280826713175363d6*u63+2.0d0*(u93-4.8d2*u70))
  value=value+(c(67)*(u296*(u77+u96*u347+u296*((u9+u389)*u296+u4*u347+u&
 &388))+u62*u347+u90))
  u56=2.82094791773878144d-1*u59
  value=value+(c(68)*(y*(u296*(u56+u386+(u389+u9)*u296)+u387+u62)*z))
  u389=-2.68661706451312518d-2*u68*(3.0d0*a3*(u66-u360)+1.4d1*u393)
  u66=u392*(2.0d0*(u84-1.92d2*u82)+1.05523040013659897d5*u63)
  u360=+u92*(u84+u47)
  u393=u391*u78
  u47=+u64*(1.49239728019318997d6*u63+8.0d0*(a3*(-u363*(pd2+6.0d0)-u368&
 &)-2.4d1*u393))
  u391=5.02993157398445508d6*u63
  u363=-3.2d1*u6*u78
  u78=u392*(2.0d0*(u71+u363)+u391)
  u368=+u64*(u57+u365)
  u6=u89*(u382+u93)
  u382=(u366+u46)
  value=value+(c(69)*(u296*(u66+u347*(u6*u347+u47)+u296*(u382*u296+u347&
 &*(u78+u46)+u360))+u347*(u48+u368*u347)+u389))
  u389=u97*u296
  value=value+(c(70)*(u370*(u348*(u8+u362+u296*(u389+u50)+(u389+u38)*u3&
 &48)+u296*(u58+u38*u296)+u364*u347+u61)))
  u78=u39*(2.0d0*u374+u69)
  u374=+u60*(4.97465760064396656d5*u63+4.0d0*(u367-1.6d1*u393))
  u6=u392*(2.0d0*(u390+u363)+u391)
  value=value+(c(71)*(u348*(u45+u347*(u394*u347+u374)+u348*(u382*u348+u&
 &347*(u6+u46)+u44))+u347*(u78+u372*u347)+u53))
  u78=-3.9375d2*u68*(6.3d1*u63+u49)
  u49=7.05236979434695359d-1*u68*(4.0d0*(u84-3.36d2*u82)+3.869178133834&
 &19621d5*u63)
  u6=+u64*(u391+2.0d0*(u93-6.72d2*u70))
  u93=u99*u296
  u70=(u296*(u49+u296*(u93+u6))+u78)
  value=value+(c(72)*(u395*u70))
  value=value+(c(73)*(x*(u296*(u96+u296*(u93+u4))+u62)*z))
  u64=(u91+u296*(u93+u83))*u348
  u390=u9*u296
  value=value+(c(74)*(u395*(u64+u296*(u56+u390)+u62)))
  u391=u88*u296
  value=value+(c(75)*(u369*((u48+u296*(u93+u36))*u348+u296*(u55+u391)+u&
 &65)))
  u55=(u93+u361)*u348
  u392=u379*u296
  u393=u95*u296
  value=value+(c(76)*(u395*(u348*(u375+u392+u55)+u393+u65)))
  u394=u51*u296
  u395=u98*u296
  value=value+(c(77)*(u369*(u348*(u376+u394+(u93+u381)*u348)+u395+u74))&
 &)
  u369=u5*(2.0d0*(7.84d2*u67+u380)+8.20734755661799196d4*u63)
  u380=-2.3625d3*u68*(4.2d1*u63-u378)
  u378=4.23142187660817215d0*u68*(7.0d0*u94+u84)
  u68=+u92*(u71+1.4d1*u41)
  value=value+(c(78)*(u296*(u380+u296*(u296*(u68+u93)+u378))+u369))
  value=value+(c(79)*(y*u70*z))
  value=value+(c(80)*((u62+u296*(u296*(u4+u93)+u96))*u348+u296*(u77+u29&
 &6*(u390+u388))+u90))
  value=value+(c(81)*(u370*(u64+u296*(u54+u391)+u76)))
  value=value+(c(82)*(u348*(u75+u296*(u392+u52)+(u296*(u36+u93)+u48)*u3&
 &48)+u296*(u87+u393)+u42))
  value=value+(c(83)*(u370*(u348*(u7+u394+u55)+u395+u73)))
  value=value+(c(84)*(u348*(u72+u373*u296+u348*((u381+u93)*u348+u86*u29&
 &6+u377))+u371*u296+u40))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_R_opt_2_4=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_R_opt_2_5(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_R_opt_2_5
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
  D_X_R_opt_2_5=.true.
  
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
  u99=7.5d-1*u69*A3_(p,pd,erfpd,exppd2)*p**3
  u126=y*z
  value=value+(c(1)*(u99*u126))
  if ( lmax .eq. 0 ) return
  u99=A4_(p,pd,erfpd,exppd2)
  u63=1.41047395886939072d-1*u69
  u114=p**4
  u93=exppd2*pd
  u79=-8.0d0*u93
  u52=u63*(2.65868077635827404d1*u99+u79)*u114
  u75=x*u126
  value=value+(c(2)*(u52*u75))
  u77=-7.5d-1*u69*u99*u114
  u114=y**2
  value=value+(c(3)*((u52*u114+u77)*z))
  u99=z**2
  value=value+(c(4)*(y*(u52*u99+u77)))
  if ( lmax .eq. 1 ) return
  u97=A5_(p,pd,erfpd,exppd2)
  u84=p**5
  u90=u69*u97*u84
  u74=-3.75d0*u90
  u6=exppd2*pd2
  u85=-1.6d1*u6
  u65=u63*(1.86107654345079183d2*u97+u85)*u84
  u52=x**2
  value=value+(c(5)*((u65*u52+u74)*u126))
  u95=u65*u114
  value=value+(c(6)*(x*(u95+u74)*z))
  u77=x*y
  u92=u65*u99
  value=value+(c(7)*(u77*(u92+u74)))
  u78=-1.125d1*u90
  value=value+(c(8)*(y*(u95+u78)*z))
  u58=8.0d0*exppd2
  u98=2.82094791773878144d-2*u69*(2.65868077635827404d1*u97+u58)*u84
  value=value+(c(9)*((u74+u95)*u99+u74*u114+u98))
  value=value+(c(10)*(u126*(u92+u78)))
  if ( lmax .eq. 2 ) return
  u98=A6_(p,pd,erfpd,exppd2)
  u84=p**6
  u95=-1.6d1*u93
  u92=u69*(1.86107654345079183d2*u98+u95)*u84
  u65=-4.23142187660817215d-1*u92
  u78=2.0d0*pd2
  u67=u63*(1.67496888910571265d3*u98+u95*(9.0d0+u78))*u84
  u95=u67*u52
  value=value+(c(11)*(x*(u95+u65)*u126))
  u91=u69*u98*u84
  u84=3.75d0*u91
  u7=-1.41047395886939072d-1*u92
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
  u83=1.125d1*u91
  u89=-8.46284375321634431d-1*u92
  value=value+(c(17)*((u114*(u89+u86)+u83)*z))
  value=value+(c(18)*(y*((u65+u86)*u99+u88+u83)))
  value=value+(c(19)*(z*(u95+u65*u114+u83)))
  value=value+(c(20)*(y*(u99*(u89+u82)+u83)))
  if ( lmax .eq. 3 ) return
  u95=A7_(p,pd,erfpd,exppd2)
  u86=p**7
  u88=u69*u95*u86
  u82=7.875d1*u88
  u81=-3.2d1*u6
  u67=u69*(1.67496888910571265d3*u95+u81)*u86
  u7=-8.46284375321634431d-1*u67
  u65=(1.1d1+u78)
  u89=u69*(1.84246577801628391d4*u95+u81*u65)*u86
  u92=1.41047395886939072d-1*u89
  u90=u92*u52
  value=value+(c(21)*((u52*(u7+u90)+u82)*u126))
  u61=-4.23142187660817215d-1*u67
  u62=-1.41047395886939072d-1*u67
  u96=(u61+u90)
  u66=u62*u52+u82
  value=value+(c(22)*(x*(u96*u114+u66)*z))
  value=value+(c(23)*(u77*(u96*u99+u66)))
  u96=(u62+u90)
  u66=u61*u52+u82
  value=value+(c(24)*(y*(u96*u114+u66)*z))
  u90=4.02992559676968776d-2*u69*(5.58322963035237549d2*u95-u58)*u86
  u58=1.41047395886939072d-1*u67
  u83=5.64189583547756287d-1*u69
  u84=u83*(5.02490666731713794d3*u95+u85*(6.0d0+pd2))*u86
  u85=-1.41047395886939072d-1*u89
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
  u85=3.9375d2*u88
  u80=-1.41047395886939072d0*u67
  value=value+(c(31)*(y*(u114*(u80+u90)+u85)*z))
  u59=-6.04488839515453165d-2*u69*(1.6d1*exppd2+1.86107654345079183d2*u&
 &95)*u86
  u86=1.575d2*u88
  value=value+(c(32)*((u82+u114*(u90+u7))*u99+u114*(u86+u96)+u59))
  u96=2.3625d2*u88
  value=value+(c(33)*(u126*(u84+u66+u96)))
  value=value+(c(34)*(u99*(u86+u7*u114+(u90+u62)*u99)+u82*u114+u59))
  value=value+(c(35)*(u126*(u99*(u80+u89)+u85)))
  if ( lmax .eq. 4 ) return
  u96=A8_(p,pd,erfpd,exppd2)
  u84=p**8
  u90=-3.2d1*u93
  u66=u69*(1.67496888910571265d3*u96+u90)*u84
  u67=2.11571093830408608d0*u66
  u92=1.84246577801628391d4*u96
  u62=u69*(u92+u90*u65)*u84
  u65=-1.41047395886939072d0*u62
  u59=pd2**2
  u85=4.0d0*u59
  u61=(u78+1.1d1)
  u7=u69*(2.39520551142116908d5*u96+u90*(1.3d1*u61+u85))*u84
  u80=1.41047395886939072d-1*u7
  u88=u80*u52
  value=value+(c(36)*(x*(u52*(u65+u88)+u67)*u126))
  u89=u69*u96*u84
  u72=-7.875d1*u89
  u97=4.23142187660817215d-1*u66
  u82=8.46284375321634431d-1*u66
  u94=-8.46284375321634431d-1*u62
  u87=-1.41047395886939072d-1*u62
  u91=(u97+u52*(u88+u94))
  u9=u52*(u82+u87*u52)+u72
  value=value+(c(37)*((u91*u114+u9)*z))
  value=value+(c(38)*(y*(u91*u99+u9)))
  u91=1.26942656298245165d0*u66
  u9=-4.23142187660817215d-1*u62
  u70=(u9+u88)
  u8=u9*u52+u91
  value=value+(c(39)*(u77*(u70*u114+u8)*z))
  u98=3.38513750128653772d0*u69*(1.86107654345079183d2*u96-4.0d0*u93)*u&
 &84
  u93=4.23142187660817215d-1*u62
  u95=(1.3d1+u78)
  u60=-1.41047395886939072d-1*u69*(2.17745955583742644d4*u96+u90*u95)*u&
 &84
  u90=2.25675833419102515d0*u69*(u92+u79*(4.0d0*u61+u59))*u84
  u84=-1.41047395886939072d-1*u7
  u61=1.41047395886939072d-1*u62
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
  u75=-3.9375d2*u89
  u60=6.34713281491225823d0*u66
  u98=-2.11571093830408608d0*u62
  value=value+(c(51)*((u114*(u60+u114*(u84+u98))+u75)*z))
  u61=1.41047395886939072d0*u66
  value=value+(c(52)*(y*((u67+u114*(u84+u65))*u99+u114*(u61+u79)+u75)))
  u79=-2.3625d2*u89
  u62=2.53885312596490329d0*u66
  value=value+(c(53)*(z*(u71+u114*(u62+u70)+u79)))
  value=value+(c(54)*(y*(u99*(u62+u8+(u84+u9)*u99)+u88+u79)))
  value=value+(c(55)*(z*(u99*(u61+u65*u114+u73)+u67*u114+u75)))
  value=value+(c(56)*(y*(u99*(u60+u99*(u96+u98))+u75)))
  if ( lmax .eq. 5 ) return
  u94=A9_(p,pd,erfpd,exppd2)
  u61=p**9
  u87=u69*u94*u61
  u71=-3.54375d3*u87
  u79=-6.4d1*u6
  u65=u69*(1.84246577801628391d4*u94+u79)*u61
  u84=6.34713281491225823d0*u65
  u98=u69*(2.39520551142116908d5*u94+u79*u95)*u61
  u91=-2.11571093830408608d0*u98
  u66=(u78+1.3d1)
  u8=u69*(3.59280826713175363d6*u94+u79*(1.5d1*u66+u85))*u61
  u62=1.41047395886939072d-1*u8
  u82=u62*u52
  value=value+(c(57)*((u52*(u84+u52*(u82+u91))+u71)*u126))
  u96=2.11571093830408608d0*u65
  u88=1.41047395886939072d0*u65
  u9=-1.41047395886939072d0*u98
  u67=-1.41047395886939072d-1*u98
  u73=(u96+u52*(u82+u9))
  u92=u52*(u88+u67*u52)+u71
  value=value+(c(58)*(x*(u73*u114+u92)*z))
  value=value+(c(59)*(u77*(u73*u99+u92)))
  u73=-2.12625d3*u87
  u92=4.23142187660817215d-1*u65
  u95=2.53885312596490329d0*u65
  u60=-8.46284375321634431d-1*u98
  u75=-4.23142187660817215d-1*u98
  u80=(u92+u52*(u82+u60))
  u70=u52*(u95+u75*u52)+u73
  value=value+(c(60)*(y*(u80*u114+u70)*z))
  u93=1.67496888910571265d3*u94
  u97=-3.76126389031837525d-1*u69*(u93-4.0d0*exppd2)*u61
  u7=-4.23142187660817215d-1*u65
  u90=4.23142187660817215d-1*u69*(3.51743466712199656d4*u94-1.28d2*u6)*&
 &u61
  u6=-4.23142187660817215d-1*u69
  u76=+u6*(4.97465760064396656d5*u94+u79*(2.7d1+4.0d0*pd2))*u61
  u89=8.46284375321634431d-1*u98
  u86=-2.82094791773878144d-1*u69*(1.65821920021465552d5*u94+u79*(9.0d0&
 &+pd2))*u61
  u74=u63*(5.02993157398445508d6*u94+u79*(2.1d1*u66+u85))*u61
  u85=-1.41047395886939072d-1*u8
  u72=1.41047395886939072d-1*u98
  u63=u85*u52
  u64=u52*(u74+u63)
  value=value+(c(61)*(u114*(u92+u52*(u64+u76)+(u52*(u89+u63)+u7)*u114)+&
 &u52*(u90+u52*(u72*u52+u86))+u97))
  value=value+(c(62)*(u126*(u80*u99+u70)))
  u80=(u82+u75)
  u70=u95+u60*u52
  u8=u92*u52+u73
  value=value+(c(63)*(x*(u114*(u70+u80*u114)+u8)*z))
  u57=2.53885312596490329d0*u69*(8.37484444552856323d3*u94+u81)*u61
  u81=+u6*(2.76369866702442587d5*u94+u79*(1.5d1+u78))*u61
  u78=4.23142187660817215d-1*u98
  u6=u83*(1.19760275571058454d6*u94+u79*(5.0d0*u66+u59))*u61
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
  u80=-2.480625d4*u87
  u63=1.48099765681286025d1*u65
  u77=-2.96199531362572051d0*u98
  value=value+(c(78)*(y*(u114*(u63+u114*(u8+u77))+u80)*z))
  u70=u69*(u93+3.2d1*exppd2)*u61
  u69=2.35078993144898453d-1*u70
  u76=-1.063125d4*u87
  value=value+(c(79)*((u71+u114*(u114*(u91+u8)+u84))*u99+u114*(u76+u114&
 &*(u79+u96))+u69))
  u64=4.23142187660817215d0*u65
  value=value+(c(80)*(u126*(u74+u114*(u64+u72)+u76)))
  u72=1.41047395886939072d-1*u70
  u81=-4.2525d3*u87
  u87=5.07770625192980658d0*u65
  value=value+(c(81)*(u99*(u81+u114*(u59+u87)+(u114*(u60+u8)+u92)*u99)+&
 &u114*(u81+u7)+u72))
  value=value+(c(82)*(u126*(u99*(u64+u95+u52)+u73+u76)))
  value=value+(c(83)*(u99*(u76+u84*u114+u99*((u67+u8)*u99+u91*u114+u96)&
 &)+u71*u114+u69))
  value=value+(c(84)*(u126*(u99*(u63+u99*(u58+u77))+u80)))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_R_opt_2_5=.false.
end function

  
  

!> Compute YY X Y coulomb integral for the right hand shell 
!!
recursive function D_X_R_opt_2_6(a12,a3,r3,c,lmax,value)
  
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
  
  logical :: D_X_R_opt_2_6
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u294
  real(kind=8) :: u319
  real(kind=8) :: u346
  real(kind=8) :: u363
  real(kind=8) :: u364
  real(kind=8) :: u368
  real(kind=8) :: u369
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
  real(kind=8) :: u381
  real(kind=8) :: u382
  real(kind=8) :: u383
  real(kind=8) :: u384
  real(kind=8) :: u385
  real(kind=8) :: u386
  real(kind=8) :: u387
  real(kind=8) :: u388
  real(kind=8) :: u389
  real(kind=8) :: u39
  real(kind=8) :: u390
  real(kind=8) :: u391
  real(kind=8) :: u392
  real(kind=8) :: u393
  real(kind=8) :: u394
  real(kind=8) :: u395
  real(kind=8) :: u396
  real(kind=8) :: u397
  real(kind=8) :: u398
  real(kind=8) :: u399
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
  D_X_R_opt_2_6=.true.
  
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
  u319=A3_(p,pd,erfpd,exppd2)
  u57=p**3
  u294=u319*u57
  u363=a3**(-2)
  u43=-4.70157986289796906d-2*u363*(5.31736155271654808d0*u294+2.0d0*(2&
 &.0d0*exppd2*u57-5.31736155271654808d0*a3*A1_(p,pd,erfpd,exppd2)))
  u387=7.5d-1*u363*u294
  u294=z**2
  value=value+(c(1)*(u387*u294+u43))
  if ( lmax .eq. 0 ) return
  u43=A4_(p,pd,erfpd,exppd2)
  u387=p**2
  u54=u43*u387
  u397=-2.0d0*a3
  u95=+u397*A2_(p,pd,erfpd,exppd2)
  u9=-2.5d-1*u363*u387
  u73=+u9*(3.0d0*u54+u95)
  u38=2.65868077635827404d1*u43
  u393=1.41047395886939072d-1*u363
  u394=p**4
  u7=exppd2*pd
  u77=-8.0d0*u7
  u384=u393*(u38+u77)*u394
  u364=u384*u294
  u346=(u364+u73)
  value=value+(c(2)*(x*u346))
  value=value+(c(3)*(y*u346))
  u392=+u9*(9.0d0*u54+u95)
  value=value+(c(4)*(z*(u364+u392)))
  if ( lmax .eq. 1 ) return
  u392=4.0d0*exppd2
  u384=-5.0d0*a3*(u392+5.31736155271654808d0*u319)
  u54=A5_(p,pd,erfpd,exppd2)
  u95=u54*u387
  u9=u363*u57
  u57=9.40315972579593812d-3*u9
  u73=exppd2*u387
  u42=u57*(2.0d0*(1.2d1*u73+u384)+7.97604232907482212d1*u95)
  u398=p**5
  u80=-3.75d0*u363*u54*u398
  u87=a3*u319
  u319=-2.0d0*u87
  u68=-7.5d-1*u9
  u72=+u68*(5.0d0*u95+u319)
  u81=1.86107654345079183d2*u54
  u39=exppd2*pd2
  u52=-1.6d1*u39
  u6=(u81+u52)
  u93=u393*u6*u398
  u364=x**2
  value=value+(c(5)*((u80+u93*u364)*u294+u72*u364+u42))
  u84=x*y
  u44=u93*u294
  value=value+(c(6)*(u84*(u44+u72)))
  u89=1.5d1*u95
  u74=+u68*(u89+u319)
  u319=x*z
  u68=(u44+u74)
  value=value+(c(7)*(u319*u68))
  u346=y**2
  value=value+(c(8)*((u80+u93*u346)*u294+u72*u346+u42))
  u42=y*z
  value=value+(c(9)*(u42*u68))
  u80=u57*(2.0d0*(3.6d1*u73+u384)+2.39281269872244664d2*u95)
  u384=-1.5d0*u9*(u89-u87)
  value=value+(c(10)*(u294*(u384+u44)+u80))
  if ( lmax .eq. 2 ) return
  u80=A6_(p,pd,erfpd,exppd2)
  u384=u80*u387
  u93=+u397*u43
  u89=u363*u394
  u95=u89*(5.0d0*u384+u93)
  u87=2.25d0*u95
  u43=1.86107654345079183d2*u80
  u57=-1.6d1*u7
  u9=(u43+u57)
  u394=p**6
  u397=u9*u394
  u74=u363*u397
  u68=-4.23142187660817215d-1*u74
  u44=a3*(-u77-u38)
  u38=-1.41047395886939072d-1*u89
  u71=u7*u387
  u92=+u38*(1.86107654345079183d2*u384+2.0d0*(u44-8.0d0*u71))
  u376=1.67496888910571265d3*u80
  u47=2.0d0*pd2
  u64=+u57*(9.0d0+u47)
  u45=(u376+u64)
  u379=u45*u394
  u41=u393*u379
  u86=u41*u364
  u368=u92*u364
  value=value+(c(11)*(x*((u68+u86)*u294+u368+u87)))
  u88=7.5d-1*u95
  u95=-1.41047395886939072d-1*u74
  u74=(u95+u86)*u294
  value=value+(c(12)*(y*(u74+u368+u88)))
  u86=7.5d-1*u89*(1.5d1*u384+u93)
  u62=+u38*(5.58322963035237549d2*u384+2.0d0*(u44-2.4d1*u71))
  value=value+(c(13)*(z*(u74+u62*u364+u86)))
  u38=u41*u346
  u74=(u95+u38)*u294
  u368=u92*u346
  value=value+(c(14)*(x*(u74+u368+u88)))
  u67=u84*z
  u381=u41*u294
  value=value+(c(15)*(u67*(u381+u62)))
  u63=u9*u387
  u96=-2.82094791773878144d-1*u89
  u53=+u96*(u44+3.0d0*u63)
  u88=(u294*(u53+u381)+u86)
  value=value+(c(16)*(x*u88))
  value=value+(c(17)*(y*((u68+u38)*u294+u368+u87)))
  value=value+(c(18)*(z*(u74+u62*u346+u86)))
  value=value+(c(19)*(y*u88))
  u95=2.25d0*u89*(2.5d1*u384+u93)
  u68=+u96*(u44+5.0d0*u63)
  value=value+(c(20)*(z*(u294*(u68+u381)+u95)))
  if ( lmax .eq. 3 ) return
  u95=8.0d0*exppd2
  u68=-7.0d0*a3*(u95+2.65868077635827404d1*u54)
  u53=A7_(p,pd,erfpd,exppd2)
  u63=u53*u387
  u92=u363*u398
  u44=-1.20897767903090633d-2*u92
  u384=+u44*(9.30538271725395914d2*u63+2.0d0*(4.0d1*u73+u68))
  u88=p**7
  u93=7.875d1*u363*u53*u88
  u41=a3*u54
  u54=-2.0d0*u41
  u96=u92*(7.0d0*u63+u54)
  u89=2.25d1*u96
  u74=1.67496888910571265d3*u53
  u62=-3.2d1*u39
  u38=(u74+u62)
  u86=u38*u88
  u398=u363*u86
  u381=-8.46284375321634431d-1*u398
  u98=1.67496888910571265d3*u63
  u99=a3*(-u52-u81)
  u81=-1.41047395886939072d-1*u92
  u368=u39*u387
  u60=-1.6d1*u368
  u79=+u81*(u98+2.0d0*(u99+u60))
  u396=1.84246577801628391d4*u53
  u386=(1.1d1+u47)
  u87=+u62*u386
  u72=(u396+u87)
  u65=u72*u88
  u371=u363*u65
  u78=1.41047395886939072d-1*u371
  u5=u78*u364
  u369=u79*u364
  value=value+(c(21)*((u93+u364*(u5+u381))*u294+u364*(u89+u369)+u384))
  u90=1.125d1*u96
  u46=-4.23142187660817215d-1*u398
  u85=(u46+u5)*u294
  value=value+(c(22)*(u84*(u85+u369+u90)))
  u49=2.1d1*u63
  u97=u92*(u49+u54)
  u91=1.125d1*u97
  u75=+u81*(5.02490666731713794d3*u63+2.0d0*(u99-4.8d1*u368))
  u370=u75*u364
  value=value+(c(23)*(u319*(u85+u370+u91)))
  u96=5.58322963035237549d2*u53
  u85=8.05985119353937553d-3*u92*(5.0d0*(u96-u95)*u387-u68)
  u48=+u81*(u98+2.0d0*(2.65868077635827404d1*u41+u60))
  u98=1.41047395886939072d-1*u398
  u369=(5.02490666731713794d3*u53+u52*(6.0d0+pd2))
  u60=2.82094791773878144d-1*u92*(2.0d0*u369*u387+a3*u6)
  u6=-1.41047395886939072d-1*u371
  u371=u6*u364
  value=value+(c(24)*(u346*(u48+u364*(u371+u60)+(u371+u98)*u346)+u364*(&
 &u48+u98*u364)+u85))
  u85=3.75d0*u97
  u60=-1.41047395886939072d-1*u398
  value=value+(c(25)*(u42*((u60+u5)*u294+u370+u85)))
  u6=-4.02992559676968777d-3*u92*(2.79161481517618774d3*u63+2.0d0*(1.2d&
 &2*u73+u68))
  u98=7.5d0*u92*(u49-u41)
  u398=u38*u387
  u49=-2.82094791773878144d-1*u92
  u52=+u49*(u99+3.0d0*u398)
  value=value+(c(26)*(u294*(u98+u52*u364+(u5+u60)*u294)+u85*u364+u6))
  u5=u78*u346
  u4=(u46+u5)*u294
  u372=u79*u346
  value=value+(c(27)*(u84*(u4+u372+u90)))
  u373=u75*u346
  value=value+(c(28)*(u319*((u60+u5)*u294+u373+u85)))
  u83=u78*u294
  value=value+(c(29)*(u84*(u294*(u52+u83)+u85)))
  u399=1.125d1*u92
  u92=u399*(3.5d1*u63+u54)
  u54=+u49*(u99+5.0d0*u398)
  u49=(u294*(u54+u83)+u92)
  value=value+(c(30)*(u319*u49))
  value=value+(c(31)*((u93+u346*(u5+u381))*u294+u346*(u89+u372)+u384))
  value=value+(c(32)*(u42*(u4+u373+u91)))
  value=value+(c(33)*(u294*(u98+u52*u346+(u5+u60)*u294)+u85*u346+u6))
  value=value+(c(34)*(u42*u49))
  u370=+u44*(4.65269135862697957d3*u63+2.0d0*(2.0d2*u73+u68))
  u60=u399*(1.05d2*u63-4.0d0*u41)
  u92=-2.4d2*u368
  u93=+u81*(2.51245333365856897d4*u63+2.0d0*(u99+u92))
  value=value+(c(35)*(u294*(u60+u294*(u83+u93))+u370))
  if ( lmax .eq. 4 ) return
  u370=A8_(p,pd,erfpd,exppd2)
  u60=u370*u387
  u93=a3*u80
  u78=-2.0d0*u93
  u98=u363*u394
  u99=u98*(7.0d0*u60+u78)
  u75=-5.625d1*u99
  u68=-3.2d1*u7
  u63=(1.67496888910571265d3*u370+u68)
  u41=p**8
  u80=u363*u63*u41
  u52=2.11571093830408608d0*u80
  u44=a3*(-u57-u43)
  u43=u98*(2.0d0*(u44-1.6d1*u71)+1.67496888910571265d3*u60)
  u49=1.41047395886939072d0*u43
  u399=1.84246577801628391d4*u370
  u384=u7*u386
  u381=(u399-3.2d1*u384)
  u6=u363*u381*u41
  u394=-1.41047395886939072d0*u6
  u386=a3*(-u57*(u47+9.0d0)-u376)
  u398=1.84246577801628391d4*u60
  u4=u384*u387
  u384=-1.6d1*u4
  u79=-1.41047395886939072d-1*u98
  u81=+u79*(u398+2.0d0*(u386+u384))
  u90=pd2**2
  u54=4.0d0*u90
  u5=(u47+1.1d1)
  u85=u363*(2.39520551142116908d5*u370+u68*(1.3d1*u5+u54))*u41
  u83=1.41047395886939072d-1*u85
  u46=u83*u364
  u374=u81*u364
  value=value+(c(36)*(x*((u52+u364*(u46+u394))*u294+u364*(u49+u374)+u75&
 &)))
  u76=-1.125d1*u99
  u59=4.23142187660817215d-1*u80
  u80=8.46284375321634431d-1*u43
  u51=-8.46284375321634431d-1*u6
  u91=(u59+u364*(u46+u51))*u294
  value=value+(c(37)*(y*(u91+u364*(u80+u374)+u76)))
  u385=u98*(2.1d1*u60+u78)
  u69=-1.125d1*u385
  u40=u98*(2.0d0*(u44-4.8d1*u71)+5.02490666731713794d3*u60)
  u58=8.46284375321634431d-1*u40
  u70=5.52739733404885173d4*u60
  u8=-4.8d1*u4
  u94=+u79*(u70+2.0d0*(u386+u8))
  u375=u94*u364
  value=value+(c(38)*(z*(u91+u364*(u58+u375)+u69)))
  u91=4.0d0*(1.86107654345079183d2*u370-4.0d0*u7)*u387
  u7=8.46284375321634431d-1*u98
  u66=u7*(u91+2.65868077635827404d1*u93)
  u382=a3*u9
  u9=-4.23142187660817215d-1*u98*(u398+2.0d0*(u382+u384))
  u398=4.23142187660817215d-1*u6
  u384=(1.3d1+u47)
  u383=u384*u387
  u56=+u79*(2.17745955583742644d4*u60+2.0d0*(u382+u57*u383))
  u57=(u399+u77*(4.0d0*u5+u90))
  u399=2.82094791773878144d-1*u98
  u77=u399*(8.0d0*u57*u387+a3*u45)
  u45=-1.41047395886939072d-1*u85
  u85=1.41047395886939072d-1*u6
  u373=u45*u364
  u395=u364*(u373+u77)
  value=value+(c(39)*(x*(u346*(u9+u395+(u373+u398)*u346)+u364*(u56+u85*&
 &u364)+u66)))
  u50=4.23142187660817215d-1*u40
  u388=-4.23142187660817215d-1*u6
  value=value+(c(40)*(u67*((u388+u46)*u294+u375+u50)))
  u375=u63*u387
  u63=u98*(3.0d0*u375+u44)
  u55=8.46284375321634431d-1*u63
  u48=1.41047395886939072d-1*u40
  u380=u381*u387
  u381=a3*(-u64-u376)
  u43=-2.82094791773878144d-1*u98
  u64=+u43*(u381+3.0d0*u380)
  u376=u64*u364
  u377=u48*u364
  value=value+(c(41)*(x*(u294*(u55+u376+(u46+u388)*u294)+u377+u69)))
  u374=(u373+u85)
  u378=u398*u364
  value=value+(c(42)*(y*(u346*(u56+u395+u374*u346)+u364*(u9+u378)+u66))&
 &)
  u77=u7*(u91+8.86226925452758014d0*u93)
  u56=-1.41047395886939072d-1*u363*(2.17745955583742644d4*u370+u68*u384&
 &)*u41
  u66=a3**(-1)
  u93=-2.82094791773878144d-1*u66*u397
  u7=+u79*(u70+2.0d0*(u382+u8))
  u70=2.25675833419102515d0*u363*u57*u41
  u57=2.82094791773878144d-1*u66*u379
  value=value+(c(43)*(z*(u294*(u56+u364*(u373+u70)+u374*u294)+(u93+u57*&
 &u364)*u346+u364*(u7+u378)+u77)))
  u77=-3.75d0*u385
  u56=2.82094791773878144d-1*u63
  u93=-1.41047395886939072d-1*u6
  u85=(u46+u93)*u294
  value=value+(c(44)*(y*(u294*(u56+u376+u85)+u377+u77)))
  u70=-1.125d1*u98*(3.5d1*u60+u78)
  u57=u399*(5.0d0*u375+u44)
  u7=4.23142187660817215d-1*u98
  u45=u7*(2.0d0*(u44-8.0d1*u71)+8.37484444552856323d3*u60)
  u385=+u43*(u381+5.0d0*u380)
  value=value+(c(45)*(z*(u294*(u57+u385*u364+u85)+u45*u364+u70)))
  u85=u83*u346
  u9=(u59+u346*(u85+u51))*u294
  u379=u81*u346
  value=value+(c(46)*(x*(u9+u346*(u80+u379)+u76)))
  u380=u94*u346
  value=value+(c(47)*(u67*((u388+u85)*u294+u380+u50)))
  u63=(u85+u93)*u294
  u381=u64*u346
  u382=u48*u346
  value=value+(c(48)*(x*(u294*(u56+u381+u63)+u382+u77)))
  u374=u83*u294
  value=value+(c(49)*(u67*(u294*(u385+u374)+u45)))
  u97=u7*(4.0d0*(u44-1.2d2*u71)+2.51245333365856897d4*u60)
  u7=+u79*(2.76369866702442587d5*u60+2.0d0*(u386-2.4d2*u4))
  u68=(u294*(u97+u294*(u374+u7))+u70)
  value=value+(c(50)*(x*u68))
  value=value+(c(51)*(y*((u52+u346*(u85+u394))*u294+u346*(u49+u379)+u75&
 &)))
  value=value+(c(52)*(z*(u9+u346*(u58+u380)+u69)))
  value=value+(c(53)*(y*(u294*(u55+u381+(u85+u388)*u294)+u382+u69)))
  value=value+(c(54)*(z*(u294*(u57+u385*u346+u63)+u45*u346+u70)))
  value=value+(c(55)*(y*u68))
  u9=-5.625d1*u98*(4.9d1*u60+u78)
  u85=7.05236979434695359d-1*u98*(4.0d0*(u44-1.68d2*u71)+3.517434667121&
 &99656d4*u60)
  u93=+u79*(3.86917813383419621d5*u60+2.0d0*(u386-3.36d2*u4))
  value=value+(c(56)*(z*(u294*(u85+u294*(u374+u93))+u9)))
  if ( lmax .eq. 5 ) return
  u9=A9_(p,pd,erfpd,exppd2)
  u85=u9*u387
  u93=1.17247822237399885d4*u85
  u83=1.86107654345079183d2*u53
  u68=1.6d1*exppd2
  u64=a3*(u68+u83)
  u374=-9.0d0*u64
  u63=u363*u88
  u69=3.35827133064140647d-2*u63
  u386=1.12d2*u73
  u44=u69*(2.0d0*(u386+u374)+u93)
  u97=p**9
  u370=-3.54375d3*u363*u9*u97
  u70=a3*u53
  u7=-2.0d0*u70
  u77=u63*(9.0d0*u85+u7)
  u78=-1.18125d3*u77
  u388=-6.4d1*u39
  u59=(1.84246577801628391d4*u9+u388)
  u48=u363*u59*u97
  u53=6.34713281491225823d0*u48
  u51=1.84246577801628391d4*u85
  u380=a3*(-u62-u74)
  u382=-3.2d1*u368
  u45=u63*(2.0d0*(u380+u382)+u51)
  u91=2.11571093830408608d0*u45
  u385=(2.39520551142116908d5*u9+u388*u384)
  u74=u363*u385*u97
  u378=-2.11571093830408608d0*u74
  u377=a3*(-u62*u5-u396)
  u98=-1.41047395886939072d-1*u63
  u76=u39*u383
  u4=+u98*(2.39520551142116908d5*u85+2.0d0*(u377-3.2d1*u76))
  u398=(u47+1.3d1)
  u375=u363*(3.59280826713175363d6*u9+u388*(1.5d1*u398+u54))*u97
  u75=1.41047395886939072d-1*u375
  u373=u75*u364
  u383=u4*u364
  value=value+(c(57)*((u370+u364*(u364*(u378+u373)+u53))*u294+u364*(u78&
 &+u364*(u383+u91))+u44))
  u79=-3.9375d2*u77
  u77=2.11571093830408608d0*u48
  u99=1.41047395886939072d0*u45
  u45=-1.41047395886939072d0*u74
  u58=(u77+u364*(u373+u45))*u294
  value=value+(c(58)*(u84*(u58+u364*(u99+u383)+u79)))
  u81=2.7d1*u85
  u94=u63*(u81+u7)
  u80=-3.9375d2*u94
  u383=2.0d0*(u380-9.6d1*u368)
  u57=5.52739733404885173d4*u85
  u41=u63*(u383+u57)
  u372=1.41047395886939072d0*u41
  u5=+u98*(7.18561653426350725d5*u85+2.0d0*(u377-9.6d1*u76))
  u384=u5*u364
  value=value+(c(59)*(u319*(u58+u364*(u372+u384)+u80)))
  u43=(1.67496888910571265d3*u9-u392)*u387
  u58=-1.34330853225656259d-2*u63*(9.0d0*a3*(u83+u68)+2.8d1*u43)
  u83=1.86107654345079183d2*u70
  u68=(u83+u382)
  u382=4.23142187660817215d-1*u63
  u50=u382*(2.0d0*u68+u51)
  u392=-4.23142187660817215d-1*u48
  u51=u382*(4.0d0*u68+3.51743466712199656d4*u85)
  u68=a3*u38
  u394=-4.23142187660817215d-1*u63
  u399=u39*(2.7d1+4.0d0*pd2)
  u89=+u394*(4.97465760064396656d5*u85+4.0d0*(u68-1.6d1*u399*u387))
  u8=8.46284375321634431d-1*u74
  u6=(1.65821920021465552d5*u9+u388*(9.0d0+pd2))
  u55=-2.82094791773878144d-1*u63
  u49=+u55*(u68+u6*u387)
  u60=a3*u72
  u56=5.02993157398445508d6*u85
  u395=u39*(2.1d1*u398+u54)
  u397=1.41047395886939072d-1*u63
  u54=u397*(2.0d0*(u60-3.2d1*u395*u387)+u56)
  u72=-1.41047395886939072d-1*u375
  u375=1.41047395886939072d-1*u74
  u52=u72*u364
  u40=u364*(u54+u52)
  value=value+(c(60)*(u346*(u50+u364*(u40+u89)+(u364*(u8+u52)+u392)*u34&
 &6)+u364*(u51+u364*(u375*u364+u49))+u58))
  u71=-7.875d1*u94
  u381=4.23142187660817215d-1*u48
  u48=8.46284375321634431d-1*u41
  u376=-8.46284375321634431d-1*u74
  value=value+(c(61)*(u42*((u381+u364*(u373+u376))*u294+u364*(u48+u384)&
 &+u71)))
  u46=2.01496279838484388d-2*u63*(2.0d0*(u386-3.0d0*u64)+u93)
  u93=-1.575d2*u63*(u81-u70)
  u81=-1.575d2*u94
  u94=u59*u387
  u64=u63*(3.0d0*u94+u380)
  u59=1.69256875064326886d0*u64
  u384=1.41047395886939072d-1*u41
  u41=u385*u387
  u379=a3*(-u87-u396)
  u87=+u55*(u379+3.0d0*u41)
  u385=u87*u364
  u386=u384*u364
  value=value+(c(62)*(u294*(u93+u364*(u385+u59)+(u364*(u376+u373)+u381)&
 &*u294)+u364*(u81+u386)+u46))
  u396=(8.37484444552856323d3*u9+u62)*u387
  u38=2.53885312596490329d0*u63
  u67=u38*(u396+u83)
  u371=2.76369866702442587d5*u85
  u83=u39*(1.5d1+u47)
  u47=u83*u387
  u62=+u394*(u371+2.0d0*(u68-3.2d1*u47))
  u394=4.23142187660817215d-1*u74
  u82=(1.19760275571058454d6*u9+u388*(5.0d0*u398+u90))
  u398=2.82094791773878144d-1*u63*(2.0d0*u82*u387+u60)
  u60=(u52+u394)
  u387=u394*u364
  value=value+(c(63)*(u84*(u346*(u62+u364*(u52+u398)+u60*u346)+u364*(u6&
 &2+u387)+u67)))
  u67=u38*(u396+6.2035884781693061d1*u70)
  u398=-4.23142187660817215d-1*u363
  u62=+u398*(2.76369866702442587d5*u9-6.4d1*u83)*u97
  u396=u66*u86
  u38=-8.46284375321634431d-1*u396
  u83=+u98*(8.2910960010732776d5*u85+2.0d0*(u68-9.6d1*u47))
  u68=5.64189583547756287d-1*u363*u82*u97
  u86=u66*u65
  u47=2.82094791773878144d-1*u86
  u388=u47*u364
  value=value+(c(64)*(u319*(u294*(u62+u364*(u52+u68)+u60*u294)+(u38+u38&
 &8)*u346+u364*(u83+u387)+u67)))
  u60=8.46284375321634431d-1*u64
  u64=-4.23142187660817215d-1*u74
  u90=(u373+u64)*u294
  value=value+(c(65)*(u84*(u294*(u60+u385+u90)+u386+u71)))
  u387=u63*(4.5d1*u85+u7)
  u82=-2.3625d2*u387
  u65=u63*(5.0d0*u94+u380)
  u61=8.46284375321634431d-1*u65
  u386=u382*(2.0d0*(u380-1.6d2*u368)+9.21232889008141955d4*u85)
  u385=+u55*(u379+5.0d0*u41)
  u389=u385*u364
  u390=u386*u364
  value=value+(c(66)*(u319*(u294*(u61+u389+u90)+u390+u82)))
  u90=(u375+u52)
  u391=u8*u364
  value=value+(c(67)*(u346*(u51+u364*(u391+u89)+u346*(u90*u346+u40+u49)&
 &)+u364*(u50+u392*u364)+u58))
  u392=u72*u346
  value=value+(c(68)*(u42*(u294*(u62+u346*(u392+u68)+(u392+u394)*u294)+&
 &u346*(u83+u388+u394*u346)+u38*u364+u67)))
  u68=-2.68661706451312518d-2*u63*(3.0d0*a3*(u95-u96)+1.4d1*u43)
  u62=4.23142187660817215d-1*u363*(3.51743466712199656d4*u9-1.28d2*u39)&
 &*u97
  u394=-2.82094791773878144d-1*u363*u6*u97
  u47=-2.82094791773878144d-1*u396
  u83=2.82094791773878144d-1*u396
  u38=+u398*(4.97465760064396656d5*u9-6.4d1*u399)*u97
  u396=u393*(5.02993157398445508d6*u9-6.4d1*u395)*u97
  u6=1.12837916709551257d0*u66*u369*u88
  u395=-2.82094791773878144d-1*u86
  u398=+u98*(u57+u383)
  u393=u395*u364
  value=value+(c(69)*(u294*(u62+u364*(u391+u38)+u294*(u90*u294+u364*(u3&
 &96+u52)+u394))+u346*(u47+u364*(u393+u6)+(u393+u83)*u346)+u364*(u384+u&
 &398*u364)+u68))
  u68=-7.875d1*u387
  u62=2.82094791773878144d-1*u65
  u65=-1.41047395886939072d-1*u74
  value=value+(c(70)*(u42*(u294*(u62+u389+(u373+u65)*u294)+u390+u68)))
  u47=6.71654266128281294d-3*u63*(2.0d0*(5.6d2*u73+u374)+5.862391111869&
 &99426d4*u85)
  u83=-7.875d1*u63*(1.35d2*u85-4.0d0*u70)
  u38=u397*(2.0d0*(u380-4.8d2*u368)+u371)
  u389=u382*(4.0d0*(u380+u92)+u371)
  u6=+u98*(3.59280826713175363d6*u85+2.0d0*(u377-4.8d2*u76))
  value=value+(c(71)*(u294*(u83+u389*u364+u294*((u65+u373)*u294+u6*u364&
 &+u38))+u68*u364+u47))
  u364=u75*u346
  u373=(u77+u346*(u364+u45))*u294
  u394=u4*u346
  value=value+(c(72)*(u84*(u373+u346*(u99+u394)+u79)))
  u395=u5*u346
  value=value+(c(73)*(u319*((u381+u346*(u364+u376))*u294+u346*(u48+u395&
 &)+u71)))
  u48=(u364+u64)*u294
  u396=u87*u346
  u397=u384*u346
  value=value+(c(74)*(u84*(u294*(u60+u396+u48)+u397+u71)))
  u398=u385*u346
  u399=u386*u346
  value=value+(c(75)*(u319*(u294*(u62+u398+(u364+u65)*u294)+u399+u68)))
  u363=u75*u294
  value=value+(c(76)*(u84*(u294*(u389+u294*(u363+u6))+u68)))
  u84=-3.9375d2*u63*(6.3d1*u85+u7)
  u39=7.05236979434695359d-1*u63*(4.0d0*(u380-3.36d2*u368)+3.8691781338&
 &3419621d5*u85)
  u8=+u98*(u56+2.0d0*(u377-6.72d2*u76))
  u377=(u294*(u39+u294*(u363+u8))+u84)
  value=value+(c(77)*(u319*u377))
  value=value+(c(78)*((u370+u346*(u346*(u378+u364)+u53))*u294+u346*(u78&
 &+u346*(u394+u91))+u44))
  value=value+(c(79)*(u42*(u373+u346*(u372+u395)+u80)))
  value=value+(c(80)*(u294*(u93+u346*(u396+u59)+(u346*(u376+u364)+u381)&
 &*u294)+u346*(u81+u397)+u46))
  value=value+(c(81)*(u42*(u294*(u61+u398+u48)+u399+u82)))
  value=value+(c(82)*(u294*(u83+u389*u346+u294*((u65+u364)*u294+u6*u346&
 &+u38))+u68*u346+u47))
  value=value+(c(83)*(u42*u377))
  u377=u69*(2.0d0*(7.84d2*u73+u374)+8.20734755661799196d4*u85)
  u374=-2.3625d3*u63*(4.2d1*u85-u70)
  u70=4.23142187660817215d0*u63*(7.0d0*u94+u380)
  u380=+u55*(u379+1.4d1*u41)
  value=value+(c(84)*(u294*(u374+u294*(u294*(u380+u363)+u70))+u377))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_R_opt_2_6=.false.
end function

  
  
!> Compute D X R coulomb integral 
!!
recursive function D_X_R_opt_2(a12,a3,r3,ny3,nz3,c,lmax,value)
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12      !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3       !< exponent of the third solid harmonics 
  integer     , intent(in)  :: ny3      !< y exponent of the cubic harmonic
  integer     , intent(in)  :: nz3      !< z exponent of the cubic harmonic
  real(kind=8), intent(in)  :: r3(3)    !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)     !< hermit coefficients of the two fisrt obrital produt
  integer     , intent(in)  :: lmax     !< l of maximum non 0 coeff in hermit decomposition
  real(kind=8), intent(out) :: value    !< result 
  logical                   :: D_X_R_opt_2 !< true if the routine knows how to calculate the requested value
  
  ! init return value
  D_X_R_opt_2 = .true.
  ! select case on m component
  select case ((ny3+nz3)*(ny3+nz3+1)/2+nz3+1)
    case (1)
      D_X_R_opt_2 = D_X_R_opt_2_1(a12,a3,r3,c,lmax,value)
    case (2)
      D_X_R_opt_2 = D_X_R_opt_2_2(a12,a3,r3,c,lmax,value)
    case (3)
      D_X_R_opt_2 = D_X_R_opt_2_3(a12,a3,r3,c,lmax,value)
    case (4)
      D_X_R_opt_2 = D_X_R_opt_2_4(a12,a3,r3,c,lmax,value)
    case (5)
      D_X_R_opt_2 = D_X_R_opt_2_5(a12,a3,r3,c,lmax,value)
    case (6)
      D_X_R_opt_2 = D_X_R_opt_2_6(a12,a3,r3,c,lmax,value)
    case default
      D_X_R_opt_2 = .false.
  end select
end function
  
end module
