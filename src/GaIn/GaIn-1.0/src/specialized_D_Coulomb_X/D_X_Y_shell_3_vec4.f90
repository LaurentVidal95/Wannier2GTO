module mod_D_X_Y_shell_optv4_3
  
  implicit none
  
  private
  public :: D_X_Y_shell_optv4_3
  integer, parameter :: NVEC=4
  
contains
  
  

!> Compute D X Y/R coulomb integral for the right hand shell 
!!
recursive function D_X_Y_shell_optv4_3(a12,a3,r3,c,lmax,value)
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12          !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3(NVEC)     !< exponent of the third solid harmonics 
  real(kind=8), intent(in)  :: r3(NVEC,3)   !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)         !< coeffcient of the hermit decomposition of the product of orbital 1 and 2
  integer     , intent(in)  :: lmax         !< l of maximum non 0 coeff in hermit decomposition~%
  real(kind=8), intent(out) :: value(*)     !< result value
  
  ! return value
  logical :: D_X_Y_shell_optv4_3
  
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
  real(kind=8) :: u1306(NVEC)
  real(kind=8) :: u1307(NVEC)
  real(kind=8) :: u1308(NVEC)
  real(kind=8) :: u2160(NVEC)
  real(kind=8) :: u310(NVEC)
  real(kind=8) :: u522(NVEC)
  real(kind=8) :: u523(NVEC)
  real(kind=8) :: u524(NVEC)
  real(kind=8) :: u525(NVEC)
  real(kind=8) :: u526(NVEC)
  real(kind=8) :: u527(NVEC)
  real(kind=8) :: u528(NVEC)
  real(kind=8) :: u529(NVEC)
  real(kind=8) :: u53(NVEC)
  real(kind=8) :: u530(NVEC)
  real(kind=8) :: u531(NVEC)
  real(kind=8) :: u532(NVEC)
  real(kind=8) :: u533(NVEC)
  real(kind=8) :: u534(NVEC)
  real(kind=8) :: u535(NVEC)
  real(kind=8) :: u536(NVEC)
  real(kind=8) :: u537(NVEC)
  real(kind=8) :: u538(NVEC)
  real(kind=8) :: u539(NVEC)
  real(kind=8) :: u54(NVEC)
  real(kind=8) :: u540(NVEC)
  real(kind=8) :: u541(NVEC)
  real(kind=8) :: u542(NVEC)
  real(kind=8) :: u543(NVEC)
  real(kind=8) :: u544(NVEC)
  real(kind=8) :: u545(NVEC)
  real(kind=8) :: u546(NVEC)
  real(kind=8) :: u547(NVEC)
  real(kind=8) :: u548(NVEC)
  real(kind=8) :: u549(NVEC)
  real(kind=8) :: u55(NVEC)
  real(kind=8) :: u550(NVEC)
  real(kind=8) :: u551(NVEC)
  real(kind=8) :: u552(NVEC)
  real(kind=8) :: u553(NVEC)
  real(kind=8) :: u554(NVEC)
  real(kind=8) :: u555(NVEC)
  real(kind=8) :: u556(NVEC)
  real(kind=8) :: u557(NVEC)
  real(kind=8) :: u558(NVEC)
  real(kind=8) :: u559(NVEC)
  real(kind=8) :: u56(NVEC)
  real(kind=8) :: u560(NVEC)
  real(kind=8) :: u561(NVEC)
  real(kind=8) :: u562(NVEC)
  real(kind=8) :: u563(NVEC)
  real(kind=8) :: u564(NVEC)
  real(kind=8) :: u565(NVEC)
  real(kind=8) :: u566(NVEC)
  real(kind=8) :: u567(NVEC)
  real(kind=8) :: u568(NVEC)
  real(kind=8) :: u57(NVEC)
  real(kind=8) :: u570(NVEC)
  real(kind=8) :: u571(NVEC)
  real(kind=8) :: u572(NVEC)
  real(kind=8) :: u573(NVEC)
  real(kind=8) :: u574(NVEC)
  real(kind=8) :: u575(NVEC)
  real(kind=8) :: u576(NVEC)
  real(kind=8) :: u577(NVEC)
  real(kind=8) :: u578(NVEC)
  real(kind=8) :: u579(NVEC)
  real(kind=8) :: u58(NVEC)
  real(kind=8) :: u580(NVEC)
  real(kind=8) :: u581(NVEC)
  real(kind=8) :: u582(NVEC)
  real(kind=8) :: u583(NVEC)
  real(kind=8) :: u584(NVEC)
  real(kind=8) :: u585(NVEC)
  real(kind=8) :: u586(NVEC)
  real(kind=8) :: u587(NVEC)
  real(kind=8) :: u588(NVEC)
  real(kind=8) :: u589(NVEC)
  real(kind=8) :: u59(NVEC)
  real(kind=8) :: u590(NVEC)
  real(kind=8) :: u591(NVEC)
  real(kind=8) :: u592(NVEC)
  real(kind=8) :: u593(NVEC)
  real(kind=8) :: u594(NVEC)
  real(kind=8) :: u595(NVEC)
  real(kind=8) :: u596(NVEC)
  real(kind=8) :: u597(NVEC)
  real(kind=8) :: u598(NVEC)
  real(kind=8) :: u599(NVEC)
  real(kind=8) :: u6(NVEC)
  real(kind=8) :: u60(NVEC)
  real(kind=8) :: u600(NVEC)
  real(kind=8) :: u601(NVEC)
  real(kind=8) :: u602(NVEC)
  real(kind=8) :: u603(NVEC)
  real(kind=8) :: u604(NVEC)
  real(kind=8) :: u605(NVEC)
  real(kind=8) :: u606(NVEC)
  real(kind=8) :: u607(NVEC)
  real(kind=8) :: u608(NVEC)
  real(kind=8) :: u609(NVEC)
  real(kind=8) :: u61(NVEC)
  real(kind=8) :: u610(NVEC)
  real(kind=8) :: u611(NVEC)
  real(kind=8) :: u612(NVEC)
  real(kind=8) :: u613(NVEC)
  real(kind=8) :: u614(NVEC)
  real(kind=8) :: u615(NVEC)
  real(kind=8) :: u616(NVEC)
  real(kind=8) :: u617(NVEC)
  real(kind=8) :: u618(NVEC)
  real(kind=8) :: u619(NVEC)
  real(kind=8) :: u62(NVEC)
  real(kind=8) :: u620(NVEC)
  real(kind=8) :: u621(NVEC)
  real(kind=8) :: u622(NVEC)
  real(kind=8) :: u623(NVEC)
  real(kind=8) :: u624(NVEC)
  real(kind=8) :: u625(NVEC)
  real(kind=8) :: u626(NVEC)
  real(kind=8) :: u627(NVEC)
  real(kind=8) :: u628(NVEC)
  real(kind=8) :: u629(NVEC)
  real(kind=8) :: u63(NVEC)
  real(kind=8) :: u630(NVEC)
  real(kind=8) :: u631(NVEC)
  real(kind=8) :: u632(NVEC)
  real(kind=8) :: u633(NVEC)
  real(kind=8) :: u634(NVEC)
  real(kind=8) :: u635(NVEC)
  real(kind=8) :: u636(NVEC)
  real(kind=8) :: u637(NVEC)
  real(kind=8) :: u638(NVEC)
  real(kind=8) :: u639(NVEC)
  real(kind=8) :: u64(NVEC)
  real(kind=8) :: u640(NVEC)
  real(kind=8) :: u641(NVEC)
  real(kind=8) :: u642(NVEC)
  real(kind=8) :: u643(NVEC)
  real(kind=8) :: u644(NVEC)
  real(kind=8) :: u645(NVEC)
  real(kind=8) :: u646(NVEC)
  real(kind=8) :: u647(NVEC)
  real(kind=8) :: u648(NVEC)
  real(kind=8) :: u649(NVEC)
  real(kind=8) :: u65(NVEC)
  real(kind=8) :: u650(NVEC)
  real(kind=8) :: u651(NVEC)
  real(kind=8) :: u652(NVEC)
  real(kind=8) :: u653(NVEC)
  real(kind=8) :: u654(NVEC)
  real(kind=8) :: u655(NVEC)
  real(kind=8) :: u656(NVEC)
  real(kind=8) :: u657(NVEC)
  real(kind=8) :: u658(NVEC)
  real(kind=8) :: u659(NVEC)
  real(kind=8) :: u66(NVEC)
  real(kind=8) :: u660(NVEC)
  real(kind=8) :: u661(NVEC)
  real(kind=8) :: u662(NVEC)
  real(kind=8) :: u663(NVEC)
  real(kind=8) :: u664(NVEC)
  real(kind=8) :: u665(NVEC)
  real(kind=8) :: u666(NVEC)
  real(kind=8) :: u667(NVEC)
  real(kind=8) :: u668(NVEC)
  real(kind=8) :: u669(NVEC)
  real(kind=8) :: u67(NVEC)
  real(kind=8) :: u670(NVEC)
  real(kind=8) :: u671(NVEC)
  real(kind=8) :: u672(NVEC)
  real(kind=8) :: u673(NVEC)
  real(kind=8) :: u674(NVEC)
  real(kind=8) :: u675(NVEC)
  real(kind=8) :: u676(NVEC)
  real(kind=8) :: u677(NVEC)
  real(kind=8) :: u678(NVEC)
  real(kind=8) :: u679(NVEC)
  real(kind=8) :: u68(NVEC)
  real(kind=8) :: u680(NVEC)
  real(kind=8) :: u681(NVEC)
  real(kind=8) :: u682(NVEC)
  real(kind=8) :: u683(NVEC)
  real(kind=8) :: u684(NVEC)
  real(kind=8) :: u685(NVEC)
  real(kind=8) :: u686(NVEC)
  real(kind=8) :: u687(NVEC)
  real(kind=8) :: u688(NVEC)
  real(kind=8) :: u689(NVEC)
  real(kind=8) :: u69(NVEC)
  real(kind=8) :: u690(NVEC)
  real(kind=8) :: u691(NVEC)
  real(kind=8) :: u692(NVEC)
  real(kind=8) :: u693(NVEC)
  real(kind=8) :: u694(NVEC)
  real(kind=8) :: u695(NVEC)
  real(kind=8) :: u696(NVEC)
  real(kind=8) :: u697(NVEC)
  real(kind=8) :: u698(NVEC)
  real(kind=8) :: u699(NVEC)
  real(kind=8) :: u7(NVEC)
  real(kind=8) :: u70(NVEC)
  real(kind=8) :: u700(NVEC)
  real(kind=8) :: u701(NVEC)
  real(kind=8) :: u702(NVEC)
  real(kind=8) :: u703(NVEC)
  real(kind=8) :: u704(NVEC)
  real(kind=8) :: u705(NVEC)
  real(kind=8) :: u706(NVEC)
  real(kind=8) :: u707(NVEC)
  real(kind=8) :: u708(NVEC)
  real(kind=8) :: u709(NVEC)
  real(kind=8) :: u71(NVEC)
  real(kind=8) :: u710(NVEC)
  real(kind=8) :: u711(NVEC)
  real(kind=8) :: u712(NVEC)
  real(kind=8) :: u713(NVEC)
  real(kind=8) :: u714(NVEC)
  real(kind=8) :: u715(NVEC)
  real(kind=8) :: u716(NVEC)
  real(kind=8) :: u717(NVEC)
  real(kind=8) :: u718(NVEC)
  real(kind=8) :: u719(NVEC)
  real(kind=8) :: u72(NVEC)
  real(kind=8) :: u720(NVEC)
  real(kind=8) :: u721(NVEC)
  real(kind=8) :: u722(NVEC)
  real(kind=8) :: u723(NVEC)
  real(kind=8) :: u724(NVEC)
  real(kind=8) :: u725(NVEC)
  real(kind=8) :: u726(NVEC)
  real(kind=8) :: u727(NVEC)
  real(kind=8) :: u728(NVEC)
  real(kind=8) :: u729(NVEC)
  real(kind=8) :: u73(NVEC)
  real(kind=8) :: u730(NVEC)
  real(kind=8) :: u731(NVEC)
  real(kind=8) :: u732(NVEC)
  real(kind=8) :: u733(NVEC)
  real(kind=8) :: u734(NVEC)
  real(kind=8) :: u735(NVEC)
  real(kind=8) :: u736(NVEC)
  real(kind=8) :: u737(NVEC)
  real(kind=8) :: u738(NVEC)
  real(kind=8) :: u739(NVEC)
  real(kind=8) :: u74(NVEC)
  real(kind=8) :: u740(NVEC)
  real(kind=8) :: u741(NVEC)
  real(kind=8) :: u742(NVEC)
  real(kind=8) :: u743(NVEC)
  real(kind=8) :: u744(NVEC)
  real(kind=8) :: u745(NVEC)
  real(kind=8) :: u746(NVEC)
  real(kind=8) :: u747(NVEC)
  real(kind=8) :: u748(NVEC)
  real(kind=8) :: u749(NVEC)
  real(kind=8) :: u75(NVEC)
  real(kind=8) :: u750(NVEC)
  real(kind=8) :: u751(NVEC)
  real(kind=8) :: u752(NVEC)
  real(kind=8) :: u753(NVEC)
  real(kind=8) :: u754(NVEC)
  real(kind=8) :: u755(NVEC)
  real(kind=8) :: u756(NVEC)
  real(kind=8) :: u757(NVEC)
  real(kind=8) :: u758(NVEC)
  real(kind=8) :: u759(NVEC)
  real(kind=8) :: u76(NVEC)
  real(kind=8) :: u760(NVEC)
  real(kind=8) :: u761(NVEC)
  real(kind=8) :: u762(NVEC)
  real(kind=8) :: u763(NVEC)
  real(kind=8) :: u764(NVEC)
  real(kind=8) :: u765(NVEC)
  real(kind=8) :: u766(NVEC)
  real(kind=8) :: u767(NVEC)
  real(kind=8) :: u768(NVEC)
  real(kind=8) :: u769(NVEC)
  real(kind=8) :: u77(NVEC)
  real(kind=8) :: u770(NVEC)
  real(kind=8) :: u771(NVEC)
  real(kind=8) :: u772(NVEC)
  real(kind=8) :: u773(NVEC)
  real(kind=8) :: u774(NVEC)
  real(kind=8) :: u775(NVEC)
  real(kind=8) :: u776(NVEC)
  real(kind=8) :: u777(NVEC)
  real(kind=8) :: u778(NVEC)
  real(kind=8) :: u779(NVEC)
  real(kind=8) :: u78(NVEC)
  real(kind=8) :: u780(NVEC)
  real(kind=8) :: u781(NVEC)
  real(kind=8) :: u782(NVEC)
  real(kind=8) :: u783(NVEC)
  real(kind=8) :: u784(NVEC)
  real(kind=8) :: u785(NVEC)
  real(kind=8) :: u786(NVEC)
  real(kind=8) :: u787(NVEC)
  real(kind=8) :: u788(NVEC)
  real(kind=8) :: u789(NVEC)
  real(kind=8) :: u79(NVEC)
  real(kind=8) :: u790(NVEC)
  real(kind=8) :: u791(NVEC)
  real(kind=8) :: u792(NVEC)
  real(kind=8) :: u793(NVEC)
  real(kind=8) :: u794(NVEC)
  real(kind=8) :: u795(NVEC)
  real(kind=8) :: u796(NVEC)
  real(kind=8) :: u797(NVEC)
  real(kind=8) :: u798(NVEC)
  real(kind=8) :: u799(NVEC)
  real(kind=8) :: u8(NVEC)
  real(kind=8) :: u80(NVEC)
  real(kind=8) :: u800(NVEC)
  real(kind=8) :: u801(NVEC)
  real(kind=8) :: u802(NVEC)
  real(kind=8) :: u803(NVEC)
  real(kind=8) :: u804(NVEC)
  real(kind=8) :: u805(NVEC)
  real(kind=8) :: u806(NVEC)
  real(kind=8) :: u807(NVEC)
  real(kind=8) :: u808(NVEC)
  real(kind=8) :: u809(NVEC)
  real(kind=8) :: u81(NVEC)
  real(kind=8) :: u810(NVEC)
  real(kind=8) :: u811(NVEC)
  real(kind=8) :: u812(NVEC)
  real(kind=8) :: u813(NVEC)
  real(kind=8) :: u814(NVEC)
  real(kind=8) :: u815(NVEC)
  real(kind=8) :: u816(NVEC)
  real(kind=8) :: u817(NVEC)
  real(kind=8) :: u818(NVEC)
  real(kind=8) :: u819(NVEC)
  real(kind=8) :: u82(NVEC)
  real(kind=8) :: u820(NVEC)
  real(kind=8) :: u821(NVEC)
  real(kind=8) :: u822(NVEC)
  real(kind=8) :: u823(NVEC)
  real(kind=8) :: u824(NVEC)
  real(kind=8) :: u825(NVEC)
  real(kind=8) :: u826(NVEC)
  real(kind=8) :: u827(NVEC)
  real(kind=8) :: u828(NVEC)
  real(kind=8) :: u829(NVEC)
  real(kind=8) :: u83(NVEC)
  real(kind=8) :: u830(NVEC)
  real(kind=8) :: u831(NVEC)
  real(kind=8) :: u832(NVEC)
  real(kind=8) :: u833(NVEC)
  real(kind=8) :: u834(NVEC)
  real(kind=8) :: u835(NVEC)
  real(kind=8) :: u836(NVEC)
  real(kind=8) :: u837(NVEC)
  real(kind=8) :: u838(NVEC)
  real(kind=8) :: u839(NVEC)
  real(kind=8) :: u84(NVEC)
  real(kind=8) :: u840(NVEC)
  real(kind=8) :: u841(NVEC)
  real(kind=8) :: u842(NVEC)
  real(kind=8) :: u843(NVEC)
  real(kind=8) :: u844(NVEC)
  real(kind=8) :: u845(NVEC)
  real(kind=8) :: u846(NVEC)
  real(kind=8) :: u847(NVEC)
  real(kind=8) :: u848(NVEC)
  real(kind=8) :: u849(NVEC)
  real(kind=8) :: u85(NVEC)
  real(kind=8) :: u850(NVEC)
  real(kind=8) :: u851(NVEC)
  real(kind=8) :: u852(NVEC)
  real(kind=8) :: u853(NVEC)
  real(kind=8) :: u854(NVEC)
  real(kind=8) :: u855(NVEC)
  real(kind=8) :: u856(NVEC)
  real(kind=8) :: u857(NVEC)
  real(kind=8) :: u858(NVEC)
  real(kind=8) :: u859(NVEC)
  real(kind=8) :: u86(NVEC)
  real(kind=8) :: u860(NVEC)
  real(kind=8) :: u861(NVEC)
  real(kind=8) :: u862(NVEC)
  real(kind=8) :: u863(NVEC)
  real(kind=8) :: u864(NVEC)
  real(kind=8) :: u865(NVEC)
  real(kind=8) :: u866(NVEC)
  real(kind=8) :: u867(NVEC)
  real(kind=8) :: u868(NVEC)
  real(kind=8) :: u869(NVEC)
  real(kind=8) :: u87(NVEC)
  real(kind=8) :: u870(NVEC)
  real(kind=8) :: u871(NVEC)
  real(kind=8) :: u872(NVEC)
  real(kind=8) :: u873(NVEC)
  real(kind=8) :: u874(NVEC)
  real(kind=8) :: u875(NVEC)
  real(kind=8) :: u876(NVEC)
  real(kind=8) :: u877(NVEC)
  real(kind=8) :: u878(NVEC)
  real(kind=8) :: u879(NVEC)
  real(kind=8) :: u88(NVEC)
  real(kind=8) :: u880(NVEC)
  real(kind=8) :: u881(NVEC)
  real(kind=8) :: u882(NVEC)
  real(kind=8) :: u883(NVEC)
  real(kind=8) :: u884(NVEC)
  real(kind=8) :: u885(NVEC)
  real(kind=8) :: u886(NVEC)
  real(kind=8) :: u887(NVEC)
  real(kind=8) :: u888(NVEC)
  real(kind=8) :: u889(NVEC)
  real(kind=8) :: u89(NVEC)
  real(kind=8) :: u890(NVEC)
  real(kind=8) :: u891(NVEC)
  real(kind=8) :: u892(NVEC)
  real(kind=8) :: u893(NVEC)
  real(kind=8) :: u894(NVEC)
  real(kind=8) :: u895(NVEC)
  real(kind=8) :: u896(NVEC)
  real(kind=8) :: u897(NVEC)
  real(kind=8) :: u898(NVEC)
  real(kind=8) :: u899(NVEC)
  real(kind=8) :: u9(NVEC)
  real(kind=8) :: u90(NVEC)
  real(kind=8) :: u900(NVEC)
  real(kind=8) :: u901(NVEC)
  real(kind=8) :: u902(NVEC)
  real(kind=8) :: u903(NVEC)
  real(kind=8) :: u904(NVEC)
  real(kind=8) :: u905(NVEC)
  real(kind=8) :: u906(NVEC)
  real(kind=8) :: u907(NVEC)
  real(kind=8) :: u908(NVEC)
  real(kind=8) :: u909(NVEC)
  real(kind=8) :: u91(NVEC)
  real(kind=8) :: u910(NVEC)
  real(kind=8) :: u911(NVEC)
  real(kind=8) :: u912(NVEC)
  real(kind=8) :: u913(NVEC)
  real(kind=8) :: u914(NVEC)
  real(kind=8) :: u915(NVEC)
  real(kind=8) :: u916(NVEC)
  real(kind=8) :: u917(NVEC)
  real(kind=8) :: u918(NVEC)
  real(kind=8) :: u919(NVEC)
  real(kind=8) :: u92(NVEC)
  real(kind=8) :: u920(NVEC)
  real(kind=8) :: u921(NVEC)
  real(kind=8) :: u922(NVEC)
  real(kind=8) :: u923(NVEC)
  real(kind=8) :: u924(NVEC)
  real(kind=8) :: u925(NVEC)
  real(kind=8) :: u926(NVEC)
  real(kind=8) :: u927(NVEC)
  real(kind=8) :: u928(NVEC)
  real(kind=8) :: u929(NVEC)
  real(kind=8) :: u93(NVEC)
  real(kind=8) :: u930(NVEC)
  real(kind=8) :: u931(NVEC)
  real(kind=8) :: u932(NVEC)
  real(kind=8) :: u933(NVEC)
  real(kind=8) :: u934(NVEC)
  real(kind=8) :: u935(NVEC)
  real(kind=8) :: u936(NVEC)
  real(kind=8) :: u937(NVEC)
  real(kind=8) :: u938(NVEC)
  real(kind=8) :: u939(NVEC)
  real(kind=8) :: u94(NVEC)
  real(kind=8) :: u940(NVEC)
  real(kind=8) :: u941(NVEC)
  real(kind=8) :: u942(NVEC)
  real(kind=8) :: u943(NVEC)
  real(kind=8) :: u944(NVEC)
  real(kind=8) :: u945(NVEC)
  real(kind=8) :: u946(NVEC)
  real(kind=8) :: u948(NVEC)
  real(kind=8) :: u949(NVEC)
  real(kind=8) :: u95(NVEC)
  real(kind=8) :: u950(NVEC)
  real(kind=8) :: u951(NVEC)
  real(kind=8) :: u952(NVEC)
  real(kind=8) :: u953(NVEC)
  real(kind=8) :: u955(NVEC)
  real(kind=8) :: u96(NVEC)
  real(kind=8) :: u960(NVEC)
  real(kind=8) :: u961(NVEC)
  real(kind=8) :: u962(NVEC)
  real(kind=8) :: u964(NVEC)
  real(kind=8) :: u965(NVEC)
  real(kind=8) :: u966(NVEC)
  real(kind=8) :: u967(NVEC)
  real(kind=8) :: u968(NVEC)
  real(kind=8) :: u969(NVEC)
  real(kind=8) :: u97(NVEC)
  real(kind=8) :: u970(NVEC)
  real(kind=8) :: u972(NVEC)
  real(kind=8) :: u975(NVEC)
  real(kind=8) :: u976(NVEC)
  real(kind=8) :: u977(NVEC)
  real(kind=8) :: u978(NVEC)
  real(kind=8) :: u979(NVEC)
  real(kind=8) :: u98(NVEC)
  real(kind=8) :: u980(NVEC)
  real(kind=8) :: u981(NVEC)
  real(kind=8) :: u982(NVEC)
  real(kind=8) :: u984(NVEC)
  real(kind=8) :: u985(NVEC)
  real(kind=8) :: u986(NVEC)
  real(kind=8) :: u99(NVEC)
  real(kind=8) :: u990(NVEC)
  real(kind=8) :: u991(NVEC)
  real(kind=8) :: u992(NVEC)
  real(kind=8) :: u993(NVEC)
  real(kind=8) :: u994(NVEC)
  real(kind=8) :: u995(NVEC)
  real(kind=8) :: u996(NVEC)
  real(kind=8) :: u997(NVEC)
  real(kind=8) :: u998(NVEC)
  real(kind=8) :: u999(NVEC)
  real(kind=8) :: value_1_(NVEC)
  real(kind=8) :: value_2_(NVEC)
  real(kind=8) :: value_3_(NVEC)
  real(kind=8) :: value_4_(NVEC)
  real(kind=8) :: value_5_(NVEC)
  real(kind=8) :: value_6_(NVEC)
  real(kind=8) :: value_7_(NVEC)

  
  ! init values
  value_1_ = 0.0d0
  value_2_ = 0.0d0
  value_3_ = 0.0d0
  value_4_ = 0.0d0
  value_5_ = 0.0d0
  value_6_ = 0.0d0
  value_7_ = 0.0d0
  
  ! init computation flag
  D_X_Y_shell_optv4_3=.true.
  
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
  u539=a3**(-3)
  u785=u539*p**4*pd*A5_v(pd,exppd2,erfpd)
  u860=1.10633173111245658d0*u785
  u784=-3.31899519333736975d0*u785
  u573=+u784
  u1306=x**2
  u1307=y**2
  value_1_=value_1_+(c(1)*(y*(u860*u1307+u573*u1306)))
  u573=-5.41989645495103886d0*u785
  u310=x*y
  u545=u310*z
  value_2_=value_2_+(c(1)*(u573*u545))
  u573=-3.42784349598349302d0*u785
  u859=8.56960873995873256d-1*u785
  u1308=z**2
  u2160=(u573*u1308+u859*u1307+u859*u1306)
  value_3_=value_3_+(c(1)*(y*u2160))
  u786=-1.39941124721293272d0*u785
  u649=2.09911687081939908d0*u785
  value_4_=value_4_+(c(1)*(z*(u786*u1308+u649*u1307+u649*u1306)))
  value_5_=value_5_+(c(1)*(x*u2160))
  u859=2.70994822747551943d0*u785
  u573=-u859
  value_6_=value_6_+(c(1)*((u859*u1307+u573*u1306)*z))
  u573=-u784
  u859=-u860
  value_7_=value_7_+(c(1)*(x*(u573*u1307+u859*u1306)))
  if ( lmax .eq. 0 ) go to 100
  u573=A7_v(pd,exppd2,erfpd)
  u859=3.2986722862692829d2*u573
  u2160=2.83592616144882564d1*exppd2
  u785=(u859+u2160)
  u649=u539*p**5*pd2
  u860=u649*u785
  u784=2.01232187092526948d-2*u860
  u786=u649*u573
  u654=7.74432211778719608d0*u786
  u896=-2.32329663533615882d1*u786
  u916=+u896
  u787=u654*u1307
  u990=u916*u1306
  u655=u787+u990
  value_1_=value_1_+(c(2)*(u310*(u655+u784)))
  u784=1.64305392733656736d-2*u860
  u562=-3.7939275184657272d1*u786
  u730=y*z
  value_2_=value_2_+(c(2)*((u562*u1306+u784)*u730))
  u773=(u2160+u859)
  u859=u649*u773
  u649=-5.19579272886834657d-3*u859
  u965=-2.39949044718844512d1*u786
  u533=5.99872611797111279d0*u786
  u847=u533*u1306
  u881=u533*u1307
  u658=u965*u1308+u881+u847
  u811=(u310*(u658+u649))
  value_3_=value_3_+(c(2)*u811)
  u650=-1.27270409949904333d-2*u859
  u862=-9.79587873049052903d0*u786
  u976=1.46938180957357935d1*u786
  u756=x*z
  u548=u862*u1308
  u886=u976*u1307
  u589=u976*u1306
  u861=(u548+u886+u589+u650)
  value_4_=value_4_+(c(2)*(u756*u861))
  u809=1.03915854577366931d-2*u860
  u651=-2.59789636443417329d-3*u859
  u652=-7.79368909330251986d-3*u859
  u950=u652+u847
  value_5_=value_5_+(c(2)*((u809+u965*u1306)*u1308+(u651+u847)*u1307+u1&
 &306*(u950)))
  u847=1.8969637592328636d1*u786
  u630=-u847
  u550=u847*u1307
  u796=u630*u1306
  u634=u550+u796
  value_6_=value_6_+(c(2)*(x*(u634+u784)*z))
  u653=-1.00616093546263474d-2*u859
  u810=1.00616093546263474d-2*u860
  u979=-u896
  u896=-u654
  u786=u896*u1306
  value_7_=value_7_+(c(2)*((u653+u979*u1306)*u1307+u1306*(u810+u786)))
  value_1_=value_1_+(c(3)*(u1307*(u653+u990+u787)+u810*u1306))
  value_2_=value_2_+(c(3)*(x*(u562*u1307+u784)*z))
  value_3_=value_3_+(c(3)*((u809+u965*u1307)*u1308+u1307*(u950+u881)+u6&
 &51*u1306))
  value_4_=value_4_+(c(3)*(u730*u861))
  value_5_=value_5_+(c(3)*u811)
  u811=-1.64305392733656736d-2*u859
  value_6_=value_6_+(c(3)*(y*(u634+u811)*z))
  u811=-2.01232187092526948d-2*u859
  u634=u979*u1307+u786
  value_7_=value_7_+(c(3)*(u310*(u634+u811)))
  value_1_=value_1_+(c(4)*(y*(u655)*z))
  value_2_=value_2_+(c(4)*(u310*(u562*u1308+u784)))
  u811=2.07831709154733863d-2*u860
  u916=(u658+u811)
  value_3_=value_3_+(c(4)*(u730*u916))
  u562=1.27270409949904333d-2*u860
  u654=-6.36352049749521663d-3*u859
  value_4_=value_4_+(c(4)*(u1308*(u562+u589+u886+u548)+u654*u1307+u654*&
 &u1306))
  value_5_=value_5_+(c(4)*(u756*u916))
  u562=-8.21526963668283677d-3*u859
  u862=8.21526963668283677d-3*u860
  value_6_=value_6_+(c(4)*((u796+u550)*u1308+u562*u1307+u862*u1306))
  value_7_=value_7_+(c(4)*(x*(u634)*z))
  if ( lmax .eq. 1 ) go to 100
  u562=u539*p**6*pd
  u862=u562*u773
  u630=-2.01232187092526948d-2*u862
  u896=u562*u573
  u533=-7.74432211778719608d0*u896
  u976=+u533
  u979=1.16164831766807941d2*u896
  u860=-5.67185232289765129d1*exppd2
  u859=+u860*pd2
  u773=u562*(2.96880505764235461d3*u573+u859)
  u573=2.3477088494128144d-2*u773
  u847=-7.04312654823844319d-2*u773
  u916=u847*u1306
  value_1_=value_1_+(c(5)*(y*((u976+u573*u1306)*u1307+u1306*(u979+u916)&
 &+u630)))
  u658=1.13817825553971816d2*u896
  u655=-1.15013774913559715d-1*u773
  u634=u655*u1306
  value_2_=value_2_+(c(5)*(x*(u634+u658)*u730))
  u861=u562*u785
  u785=5.19579272886834657d-3*u861
  u562=2.39949044718844512d1*u896
  u950=-5.99872611797111279d0*u896
  u787=-2.9993630589855564d1*u896
  u990=+u787
  u881=-7.2741098204156852d-2*u773
  u548=1.8185274551039213d-2*u773
  u886=u548*u1306
  u589=u881*u1306
  u550=(u562+u589)*u1308+(u950+u886)*u1307
  u796=(y*(u550+u1306*(u990+u886)+u785))
  value_3_=value_3_+(c(5)*u796)
  u786=1.27270409949904333d-2*u861
  u547=9.79587873049052903d0*u896
  u966=-1.46938180957357935d1*u896
  u802=-7.34690904786789677d1*u896
  u815=-2.96964289883110109d-2*u773
  u553=4.45446434824665164d-2*u773
  u683=u553*u1306
  value_4_=value_4_+(c(5)*(z*((u547+u815*u1306)*u1308+(u966+u683)*u1307&
 &+u1306*(u802+u683)+u786)))
  u812=1.55873781866050397d-2*u861
  u568=7.19847134156533535d1*u896
  u857=-1.79961783539133384d1*u896
  u659=-4.19910828257977895d1*u896
  value_5_=value_5_+(c(5)*(x*((u568+u589)*u1308+(u857+u886)*u1307+u1306&
 &*(u659+u886)+u812)))
  u589=-1.64305392733656736d-2*u862
  u951=-1.8969637592328636d1*u896
  u977=+u951
  u955=9.48481879616431801d1*u896
  u574=5.75068874567798574d-2*u773
  u816=-5.75068874567798574d-2*u773
  u867=u816*u1306
  value_6_=value_6_+(c(5)*(((u977+u574*u1306)*u1307+u1306*(u955+u867)+u&
 &589)*z))
  u564=-6.96988990600847647d1*u896
  u935=+u564
  u813=5.42102548245103726d1*u896
  u575=7.04312654823844319d-2*u773
  u595=-2.3477088494128144d-2*u773
  u660=u595*u1306
  u561=u575*u1306
  value_7_=value_7_+(c(5)*(x*((u935+u561)*u1307+u1306*(u813+u660)+u630)&
 &))
  u935=2.32329663533615882d1*u896
  u991=u573*u1307
  u897=u916+u991
  u814=u935*u1306
  value_1_=value_1_+(c(6)*(x*(u1307*(u935+u897)+u814+u630)))
  u769=3.7939275184657272d1*u896
  u822=(u769+u634)
  u634=u769*u1306+u589
  value_2_=value_2_+(c(6)*((u822*u1307+u634)*z))
  u678=u881*u1307
  u858=(u562+u678)*u1308
  u82=u548*u1307
  u789=u886+u82
  u967=u950*u1306
  u729=(x*(u858+u1307*(u990+u789)+u967+u785))
  value_3_=value_3_+(c(6)*u729)
  u750=-5.87752723829431742d1*u896
  u586=u815*u1308
  u592=u553*u1307
  value_4_=value_4_+(c(6)*(u545*(u586+u592+u683+u750)))
  value_5_=value_5_+(c(6)*u796)
  u750=u574*u1307
  value_6_=value_6_+(c(6)*(u310*(u750+u867)*z))
  u796=2.01232187092526948d-2*u861
  u982=-u935
  u566=(u982+u561)*u1307
  value_7_=value_7_+(c(6)*(y*(u566+u1306*(u982+u660)+u796)))
  u561=4.64659327067231765d1*u896
  value_1_=value_1_+(c(7)*(u310*(u991+u916+u561)*z))
  value_2_=value_2_+(c(7)*(y*(u822*u1308+u634)))
  u822=3.59923567078266768d1*u896
  u634=u881*u1308
  u916=(u545*(u634+u82+u886+u822))
  value_3_=value_3_+(c(7)*u916)
  u991=u683+u592
  u683=u991+u586
  u82=u966*u1306+u786
  u586=(u1308*(u683)+u966*u1307+u82)
  value_4_=value_4_+(c(7)*(x*u586))
  u592=-2.07831709154733863d-2*u862
  u777=-u787
  value_5_=value_5_+(c(7)*(z*(u550+u1306*(u777+u886)+u592)))
  u787=-u951
  u951=u867+u750
  u867=u787*u1306
  u886=u977*u1307+u867
  value_6_=value_6_+(c(7)*(x*((u769+u951)*u1308+u886+u589)))
  value_7_=value_7_+(c(7)*((u566+u1306*(u935+u660))*z))
  u750=-u813
  u813=-u564
  value_1_=value_1_+(c(8)*(y*(u1307*(u750+u897)+u813*u1306+u796)))
  u630=u655*u1307
  value_2_=value_2_+(c(8)*(u310*(u630+u658)*z))
  value_3_=value_3_+(c(8)*(y*((u568+u678)*u1308+u1307*(u659+u789)+u857*&
 &u1306+u812)))
  value_4_=value_4_+(c(8)*(z*((u547+u815*u1307)*u1308+u1307*(u802+u991)&
 &+u82)))
  value_5_=value_5_+(c(8)*u729)
  u813=1.64305392733656736d-2*u861
  u857=-u955
  value_6_=value_6_+(c(8)*((u1307*(u857+u951)+u867+u813)*z))
  u955=-u979
  u857=-u533
  u861=u575*u1307
  u533=u660+u861
  u991=u857*u1306
  value_7_=value_7_+(c(8)*(x*(u1307*(u955+u533)+u991+u796)))
  value_1_=value_1_+(c(9)*((u1307*(u982+u897)+u814)*z))
  value_2_=value_2_+(c(9)*(x*((u769+u630)*u1308+u769*u1307+u589)))
  u979=u967+u592
  value_3_=value_3_+(c(9)*(z*(u858+u1307*(u777+u789)+u979)))
  value_4_=value_4_+(c(9)*(y*u586))
  value_5_=value_5_+(c(9)*u916)
  u822=-u769
  value_6_=value_6_+(c(9)*(y*((u822+u951)*u1308+u886+u813)))
  u822=-u561
  value_7_=value_7_+(c(9)*(u310*(u861+u660+u822)*z))
  value_1_=value_1_+(c(10)*(y*((u897)*u1308+u976*u1307+u814)))
  value_2_=value_2_+(c(10)*(u545*(u655*u1308+u658)))
  u822=1.19974522359422256d2*u896
  u562=(u1308*(u822+u789+u634)+u950*u1307+u979)
  value_3_=value_3_+(c(10)*(y*u562))
  u777=-2.54540819899808665d-2*u862
  u658=6.85711511134337032d1*u896
  u787=-4.40814542872073806d1*u896
  value_4_=value_4_+(c(10)*(z*(u1308*(u658+u683)+u787*u1307+u787*u1306+&
 &u777)))
  value_5_=value_5_+(c(10)*(x*u562))
  u822=-5.69089127769859081d1*u896
  u658=+u822
  u787=-u822
  value_6_=value_6_+(c(10)*(z*((u951)*u1308+u658*u1307+u787*u1306)))
  value_7_=value_7_+(c(10)*(x*((u533)*u1308+u982*u1307+u991)))
  if ( lmax .eq. 2 ) go to 100
  u658=A9_v(pd,exppd2,erfpd)
  u787=2.96880505764235461d3*u658
  u822=u539*p**7
  u777=u822*pd2
  u857=u777*(-u860+u787)
  u562=-9.39083539765125759d-2*u857
  u896=pd2*u658
  u533=u822*u896
  u979=-2.09096697180254294d2*u533
  u955=+u979
  u630=1.88187027462228865d3*u533
  u935=3.26568556340659007d4*u658
  u862=-1.13437046457953026d2*exppd2
  u950=+u862*pd2
  u769=u777*(u935+u950)
  u561=2.3477088494128144d-2*u769
  u976=-7.04312654823844319d-2*u769
  u861=u976*u1306
  u592=u561*u1306
  value_1_=value_1_+(c(11)*(u310*((u955+u592)*u1307+u1306*(u630+u861)+u&
 &562)))
  u562=-3.83379249711865716d-2*u857
  u586=2.04872085997149269d3*u533
  u916=-1.15013774913559715d-1*u769
  u683=u916*u1306
  value_2_=value_2_+(c(11)*((u1306*(u586+u683)+u562)*u730))
  u858=u777*(u787-u860)
  u787=2.42470327347189507d-2*u858
  u897=6.47862420740880182d2*u533
  u789=-1.61965605185220046d2*u533
  u951=+u789
  u991=-4.85896815555660137d2*u533
  u886=-7.2741098204156852d-2*u769
  u547=1.8185274551039213d-2*u769
  u82=u547*u1306
  u802=u886*u1306
  u660=(u897+u802)*u1308+(u951+u82)*u1307
  u867=(u310*(u660+u1306*(u991+u82)+u787))
  value_3_=value_3_+(c(11)*u867)
  u814=5.93928579766220219d-2*u858
  u678=2.64488725723244284d2*u533
  u967=-3.96733088584866426d2*u533
  u750=-1.19019926575459928d3*u533
  u634=-2.96964289883110109d-2*u769
  u550=4.45446434824665164d-2*u769
  u566=u550*u1306
  u659=u634*u1306
  value_4_=value_4_+(c(11)*(u756*((u678+u659)*u1308+(u967+u566)*u1307+u&
 &1306*(u750+u566)+u814)))
  u568=2.96880505764235461d3*u896
  u575=2.0d0*pd2
  u881=u822*(u2160*(u575+9.0d0)+u568)
  u847=-1.73193090962278219d-3*u881
  u655=-2.42470327347189507d-2*u857
  u815=6.06175818367973767d-3*u858
  u816=5.4555823653117639d-2*u858
  u595=1.29572484148176036d3*u533
  u977=-3.23931210370440091d2*u533
  u729=-u897
  u564=u82+u729
  value_5_=value_5_+(c(11)*((u655+u1306*(u802+u595))*u1308+(u815+u1306*&
 &(u82+u977))*u1307+u1306*(u816+u1306*(u564))+u847))
  u982=-7.66758499423731432d-2*u857
  u594=-5.12180214992873173d2*u533
  u992=+u594
  u528=1.53654064497861952d3*u533
  u556=5.75068874567798574d-2*u769
  u591=-5.75068874567798574d-2*u769
  u914=u591*u1306
  u68=u556*u1306
  value_6_=value_6_+(c(11)*(x*((u992+u68)*u1307+u1306*(u528+u914)+u982)&
 &*z))
  u982=(9.0d0+u575)
  u961=u822*(u568+u2160*u982)
  u846=2.23591318991696609d-3*u961
  u631=2.3477088494128144d-2*u858
  u644=-7.04312654823844319d-2*u857
  u713=-1.25458018308152576d3*u533
  u993=+u713
  u885=8.36386788721017177d2*u533
  u557=7.04312654823844319d-2*u769
  u525=-2.3477088494128144d-2*u769
  u99=u525*u1306
  u804=u557*u1306
  value_7_=value_7_+(c(11)*((u631+u1306*(u804+u993))*u1307+u1306*(u644+&
 &u1306*(u99+u885))+u846))
  u631=-7.82569616470938133d-3*u857
  u644=-6.96988990600847647d1*u533
  u993=+u644
  u805=-3.91284808235469066d-2*u857
  u764=-u979
  u979=u764*u1306
  value_1_=value_1_+(c(12)*(u1307*(u631+u1306*(u861+u885)+(u592+u993)*u&
 &1307)+u1306*(u805+u979)+u846))
  u805=1.02436042998574635d3*u533
  u846=3.41453476661915448d2*u533
  u576=(u805+u683)
  u719=u846*u1306
  u917=u719+u562
  value_2_=value_2_+(c(12)*(x*(u576*u1307+u917)*z))
  u662=1.48440252882117731d4*u896
  u896=-1.73193090962278219d-3*u822*(u2160*(1.0d1*pd2+3.0d0)+u662)
  u888=1.48440252882117731d4*u658
  u788=1.8185274551039213d-2*u777*(u888-u860)
  u994=-2.69942675308700076d2*u533
  u878=-7.2741098204156852d-2*u777*(4.45320758646353191d4*u658+u950)
  u593=9.0926372755196065d-2*u769
  u797=u593*u1306
  u768=(u1307*(u788+u1306*(u797+u878)+(u797+u994)*u1307)+u1306*(u788+u9&
 &94*u1306)+u896)
  value_3_=value_3_+(c(12)*u768)
  u797=1.9797619325540674d-2*u858
  u615=8.81629085744147612d1*u533
  u657=-1.32244362861622142d2*u533
  u968=+u657
  u577=-9.25710540031354993d2*u533
  value_4_=value_4_+(c(12)*(u730*((u615+u659)*u1308+(u968+u566)*u1307+u&
 &1306*(u577+u566)+u797)))
  value_5_=value_5_+(c(12)*u867)
  u867=-1.70726738330957724d2*u533
  u995=+u867
  u870=-u594
  value_6_=value_6_+(c(12)*(y*((u995+u68)*u1307+u1306*(u870+u914))*z))
  u594=3.13027846588375253d-2*u858
  u668=-6.27290091540762882d2*u533
  u941=+u668
  u913=-u644
  u644=(u941+u804)*u1307
  value_7_=value_7_+(c(12)*(u310*(u644+u1306*(u913+u99)+u594)))
  u594=-1.56513923294187627d-2*u857
  u796=1.04548348590127147d3*u533
  value_1_=value_1_+(c(13)*(y*((u993+u592)*u1307+u1306*(u796+u861)+u594&
 &)*z))
  value_2_=value_2_+(c(13)*(u310*(u576*u1308+u917)))
  u576=-1.21235163673594753d-2*u857
  u917=2.15954140246960061d2*u533
  u592=-5.39885350617400152d1*u533
  u952=+u592
  u665=-u789
  u789=(u730*((u917+u802)*u1308+(u952+u82)*u1307+u1306*(u665+u82)+u576)&
 &)
  value_3_=value_3_+(c(13)*u789)
  u548=5.0d0*pd2
  u884=7.0705783305502407d-4*u822*(u662-u860*(u548-9.0d0))
  u656=-4.94940483138516849d-3*u857
  u599=-u657
  u574=-2.26874092915906051d2*exppd2
  u657=-4.94940483138516849d-3*u777*(u888+u574)
  u818=-u678
  u553=exppd2*pd2
  u778=u777*(8.90641517292706383d3*u658-2.83592616144882564d1*u553)
  u658=1.78178573929866066d-1*u778
  u722=-4.45446434824665164d-2*u769
  u819=u722*u1306
  u666=u819+u658
  value_4_=value_4_+(c(13)*(u1308*(u656+u1306*(u566+u818)+(u659+u615)*u&
 &1308)+u1307*(u968+u1306*(u666)+(u819+u599)*u1307)+u657*u1306+u884))
  u659=-3.6370549102078426d-2*u857
  u819=-u592
  value_5_=value_5_+(c(13)*(u756*(u660+u1306*(u819+u82)+u659)))
  u592=exppd2*(1.0d0+pd2)
  u604=8.21526963668283677d-3*u822*(u568+5.67185232289765129d1*u592)
  u660=-1.91689624855932858d-2*u857
  u771=-u867
  u867=u777*(-u862+u888)
  u605=-1.91689624855932858d-2*u867
  u767=2.3002754982711943d-1*u778
  value_6_=value_6_+(c(13)*((u660+u1306*(u914+u805))*u1308+u1307*(u995+&
 &u1306*(u914+u767)+(u914+u771)*u1307)+u1306*(u605+u719)+u604))
  u604=4.87892293420593353d2*u533
  value_7_=value_7_+(c(13)*(x*(u644+u1306*(u604+u99)+u594)*z))
  u605=-3.13027846588375253d-2*u857
  u767=-u668
  u668=u561*u1307
  u644=u861+u668
  u719=u767*u1306
  value_1_=value_1_+(c(14)*(u310*(u1307*(u993+u644)+u719+u605)))
  u605=(u846+u683)
  u683=u805*u1306
  u708=u683+u562
  value_2_=value_2_+(c(14)*(y*(u605*u1307+u708)*z))
  u667=u886*u1307
  u780=(u897+u667)*u1308
  u927=u547*u1307
  u89=u82+u927
  u985=u951*u1306
  u596=(u310*(u780+u1307*(u991+u89)+u985+u787))
  value_3_=value_3_+(c(14)*u596)
  u928=u550*u1307
  u67=u566+u928
  u899=u968*u1306+u797
  u62=u634*u1307
  value_4_=value_4_+(c(14)*(u756*((u615+u62)*u1308+u1307*(u577+u67)+u89&
 &9)))
  value_5_=value_5_+(c(14)*u768)
  u577=u556*u1307
  u896=u914+u577
  u768=u771*u1306
  value_6_=value_6_+(c(14)*(x*(u1307*(u992+u896)+u768)*z))
  u878=-2.23591318991696609d-3*u881
  u881=3.91284808235469066d-2*u858
  u798=7.82569616470938133d-3*u858
  u560=-u885
  u885=u913*u1306
  value_7_=value_7_+(c(14)*(u1307*(u881+u1306*(u99+u560)+(u804+u955)*u1&
 &307)+u1306*(u798+u885)+u878))
  u881=u979+u594
  value_1_=value_1_+(c(15)*(x*(u1307*(u764+u644)+u881)*z))
  u593=4.0d0*pd2
  u817=-1.0953692848910449d-2*u822*(1.41796308072441282d1*exppd2*(u593-&
 &3.0d0)+u568)
  u996=-u846
  u969=-4.6005509965423886d-1*u778
  u661=1.15013774913559715d-1*u769
  u680=u661*u1306
  value_2_=value_2_+(c(15)*(u1307*(u846+u1306*(u680+u969)+(u680+u996)*u&
 &1307)+u1306*(u846+u996*u1306)+u817))
  u817=u952*u1306+u576
  u969=(u756*((u917+u667)*u1308+u1307*(u665+u89)+u817))
  value_3_=value_3_+(c(15)*u969)
  u680=u634*u1308
  u926=u67+u680
  value_4_=value_4_+(c(15)*(u310*(u1308*(u818+u926)+u968*u1307+u899)))
  value_5_=value_5_+(c(15)*u789)
  value_6_=value_6_+(c(15)*(u310*((u896)*u1308+u995*u1307+u768)))
  u789=1.56513923294187627d-2*u858
  value_7_=value_7_+(c(15)*(y*((u955+u804)*u1307+u1306*(u955+u99)+u789)&
 &*z))
  u899=4.18193394360508588d2*u533
  value_1_=value_1_+(c(16)*(u310*((u899+u644)*u1308+u993*u1307+u881)))
  value_2_=value_2_+(c(16)*(u730*(u605*u1308+u708)))
  u605=9.71793631111320274d2*u533
  u708=u89+u886*u1308
  u881=(u310*(u1308*(u605+u708)+u952*u1307+u817))
  value_3_=value_3_+(c(16)*u881)
  u817=9.89880966277033698d-3*u858
  u768=3.52651634297659045d2*u533
  u960=u967*u1306
  u587=(u1308*(u768+u926)+u967*u1307+u960+u817)
  value_4_=value_4_+(c(16)*(u756*u587))
  u926=8.65965454811391095d-4*u822*(u662-u860*(1.2d1+u548))
  u662=-4.24323072857581637d-2*u857
  u663=-6.06175818367973767d-3*u867
  u664=7.2741098204156852d-2*u778
  u939=-1.8185274551039213d-2*u769
  u565=u939*u1306
  u919=u565+u664
  value_5_=value_5_+(c(16)*(u1308*(u662+u1306*(u82+u605)+(u802+u917)*u1&
 &308)+u1307*(u952+u1306*(u919)+(u565+u819)*u1307)+u663*u1306+u926))
  u82=u870*u1306
  u802=u992*u1307+u82
  value_6_=value_6_+(c(16)*(u756*((u846+u896)*u1308+u802+u562)))
  u565=-7.82569616470938133d-3*u777*(-u862+u935)
  u907=2.81725061929537728d-1*u778
  u563=2.78795596240339059d2*u533
  value_7_=value_7_+(c(16)*((u631+u1306*(u99+u899))*u1308+u1307*(u955+u&
 &1306*(u861+u907)+(u861+u764)*u1307)+u1306*(u565+u563*u1306)+u798))
  u565=7.04312654823844319d-2*u858
  u907=-2.3477088494128144d-2*u857
  u549=-u713
  value_1_=value_1_+(c(17)*(u1307*(u565+u549*u1306+u1307*(u668+u861+u56&
 &0))+u907*u1306+u878))
  u549=u916*u1307
  value_2_=value_2_+(c(17)*(x*(u1307*(u586+u549)+u562)*z))
  value_3_=value_3_+(c(17)*((u655+u1307*(u667+u595))*u1308+u1307*(u816+&
 &u977*u1306+u1307*(u927+u564))+u815*u1306+u847))
  value_4_=value_4_+(c(17)*(u730*((u678+u62)*u1308+u1307*(u750+u67)+u96&
 &0+u814)))
  value_5_=value_5_+(c(17)*u596)
  u595=7.66758499423731432d-2*u858
  u729=-u528
  value_6_=value_6_+(c(17)*(y*(u1307*(u729+u896)+u82+u595)*z))
  u595=9.39083539765125759d-2*u858
  u729=-u630
  u560=u99+u557*u1307
  value_7_=value_7_+(c(17)*(u310*(u1307*(u729+u560)+u979+u595)))
  u595=-u604
  value_1_=value_1_+(c(18)*(y*(u1307*(u595+u644)+u719+u789)*z))
  value_2_=value_2_+(c(18)*(u310*((u805+u549)*u1308+u846*u1307+u562)))
  value_3_=value_3_+(c(18)*(u730*(u780+u1307*(u819+u89)+u985+u659)))
  value_4_=value_4_+(c(18)*(u1308*(u656+u1307*(u928+u818)+(u62+u615)*u1&
 &308)+u1307*(u657+u1306*(u722*u1307+u666))+u1306*(u968+u599*u1306)+u88&
 &4))
  value_5_=value_5_+(c(18)*u969)
  u656=exppd2*(pd2+1.0d0)
  u595=-8.21526963668283677d-3*u822*(5.67185232289765129d1*u656+u568)
  u818=1.91689624855932858d-2*u858
  u599=1.91689624855932858d-2*u777*(u888-u862)
  u729=-u805
  u568=-2.3002754982711943d-1*u778
  value_6_=value_6_+(c(18)*((u818+u1307*(u577+u729))*u1308+u1307*(u599+&
 &u1306*(u68+u568)+(u68+u996)*u1307)+u1306*(u771+u995*u1306)+u595))
  u771=-u796
  u599=u885+u789
  value_7_=value_7_+(c(18)*(x*(u1307*(u771+u560)+u599)*z))
  u771=7.82569616470938133d-3*u777*(u935-u862)
  u595=-u899
  u568=-u563
  u665=-2.81725061929537728d-1*u778
  value_1_=value_1_+(c(19)*((u798+u1307*(u668+u595))*u1308+u1307*(u771+&
 &u1306*(u804+u665)+(u804+u568)*u1307)+u1306*(u764+u955*u1306)+u631))
  value_2_=value_2_+(c(19)*(u756*((u846+u549)*u1308+u805*u1307+u562)))
  value_3_=value_3_+(c(19)*(u1308*(u662+u1307*(u927+u605)+(u667+u917)*u&
 &1308)+u1307*(u663+u1306*(u939*u1307+u919))+u1306*(u952+u819*u1306)+u9&
 &26))
  value_4_=value_4_+(c(19)*(u730*u587))
  value_5_=value_5_+(c(19)*u881)
  u665=3.83379249711865716d-2*u858
  value_6_=value_6_+(c(19)*(u730*((u996+u896)*u1308+u802+u665)))
  u665=u955*u1307
  value_7_=value_7_+(c(19)*(u310*((u595+u560)*u1308+u665+u599)))
  value_1_=value_1_+(c(20)*(u730*((u644)*u1308+u665+u719)))
  value_2_=value_2_+(c(20)*(u310*(u1308*(u586+u916*u1308)+u562)))
  u665=-9.69881309388758027d-2*u857
  u605=1.94358726222264055d3*u533
  u586=(u1308*(u605+u708)+u951*u1307+u985+u665)
  value_3_=value_3_+(c(20)*(u730*u586))
  u917=2.82823133222009628d-3*u961
  u870=-8.90892869649330328d-2*u857
  u913=1.05795490289297714d3*u533
  u819=1.48482144941555055d-2*u858
  u771=-7.93466177169732851d2*u533
  value_4_=value_4_+(c(20)*(u1308*(u870+u771*u1306+u771*u1307+u1308*(u6&
 &80+u928+u566+u913))+u819*u1307+u819*u1306+u917))
  value_5_=value_5_+(c(20)*(u756*u586))
  value_6_=value_6_+(c(20)*(u1308*(u683+u729*u1307+(u577+u914)*u1308)+u&
 &818*u1307+u660*u1306))
  value_7_=value_7_+(c(20)*(u756*((u560)*u1308+u941*u1307+u979)))
  if ( lmax .eq. 3 ) go to 100
  u917=A11_v(pd,exppd2,erfpd)
  u870=pd2*u917
  u913=3.26568556340659007d4*u870
  u941=(u913-u860*(1.1d1+u575))
  u729=u539*p**8*pd
  u771=u729*u941
  u764=8.53712308877387053d-3*u771
  u605=3.26568556340659007d4*u917
  u586=pd2*(u605-u862)
  u767=u729*u586
  u595=6.4028423165804029d-3*u767
  u599=pd2*(-u862+u605)
  u768=u729*u599
  u560=-2.49710850346635713d-1*u768
  u533=u729*u870
  u568=-4.60012733796559447d3*u533
  u888=+u568
  u935=3.22008913657591613d4*u533
  u777=4.24539123242856709d5*u917
  u557=+u574*pd2
  u847=pd2*(u777+u557)
  u549=u729*u847
  u822=2.3477088494128144d-2*u549
  u878=-7.04312654823844319d-2*u549
  u630=u822*u1306
  u897=(u630+u888)
  u678=u878*u1306
  value_1_=value_1_+(c(21)*(y*((u595+u1306*u897)*u1307+u1306*(u560+u130&
 &6*(u678+u935))+u764)))
  u560=-1.56836965791217793d-1*u768
  u993=3.75598824328106993d4*u533
  u713=-1.15013774913559715d-1*u549
  u805=u713*u1306
  value_2_=value_2_+(c(21)*(x*(u1306*(u993+u805)+u560)*u730))
  u846=(-u860*(u575+1.1d1)+u913)
  u796=u729*u846
  u604=-2.20427570315626824d-3*u796
  u899=-1.98384813284064142d-2*u768
  u563=4.95962033210160355d-3*u767
  u596=6.44750643173208461d-2*u767
  u969=1.4252973256299364d4*u533
  u881=-3.563243314074841d3*u533
  u587=-8.31423439950796233d3*u533
  u907=+u587
  u708=-7.2741098204156852d-2*u549
  u780=1.8185274551039213d-2*u549
  u644=u780*u1306
  u896=u708*u1306
  u89=(u899+u1306*(u896+u969))*u1308+(u563+u1306*(u644+u881))*u1307
  u750=(y*(u89+u1306*(u596+u1306*(u644+u907))+u604))
  value_3_=value_3_+(c(21)*u750)
  u67=-5.39935072514745653d-3*u796
  u802=-8.0990260877211848d-3*u768
  u564=1.21485391315817772d-2*u767
  u666=1.57931008710563104d-1*u767
  u919=5.81875196591137424d3*u533
  u977=-8.72812794886706136d3*u533
  u884=+u977
  u926=-2.03656318806898098d4*u533
  u861=-2.96964289883110109d-2*u549
  u659=4.45446434824665164d-2*u549
  u914=u659*u1306
  u99=u564+u1306*(u914+u884)
  u804=u861*u1306
  value_4_=value_4_+(c(21)*(z*((u802+u1306*(u804+u919))*u1308+(u99)*u13&
 &07+u1306*(u666+u1306*(u914+u926))+u67)))
  u68=-1.10213785157813412d-2*u796
  u979=-9.91924066420320709d-2*u768
  u667=2.47981016605080177d-2*u767
  u668=1.23990508302540089d-1*u767
  u719=2.37549554271656067d4*u533
  u927=-5.93873885679140167d3*u533
  u985=-1.0689729942224523d4*u533
  u928=+u985
  value_5_=value_5_+(c(21)*(x*((u979+u1306*(u896+u719))*u1308+(u667+u13&
 &06*(u644+u927))*u1307+u1306*(u668+u1306*(u644+u928))+u68)))
  u62=6.97053181294301302d-3*u771
  u577=1.56836965791217793d-2*u767
  u565=-2.03888055528583131d-1*u768
  u683=-1.12679647298432098d4*u533
  u885=+u683
  u680=2.62919177029674895d4*u533
  u576=5.75068874567798574d-2*u549
  u82=-5.75068874567798574d-2*u549
  u615=u82*u1306
  u961=u576*u1306
  value_6_=value_6_+(c(21)*(((u577+u1306*(u961+u885))*u1307+u1306*(u565&
 &+u1306*(u615+u680))+u62)*z))
  u62=1.42285384812897842d-2*u771
  u565=9.60426347487060436d-2*u767
  u886=-1.60071057914510073d-1*u768
  u976=-2.30006366898279724d4*u533
  u591=+u976
  u634=1.38003820138967834d4*u533
  u525=7.04312654823844319d-2*u549
  u916=-2.3477088494128144d-2*u549
  u722=u916*u1306
  u939=u525*u1306
  value_7_=value_7_+(c(21)*(x*((u565+u1306*(u939+u591))*u1307+u1306*(u8&
 &86+u1306*(u722+u634))+u62)))
  u62=-5.76255808492236261d-2*u768
  u565=-2.30006366898279723d3*u533
  u886=+u565
  u591=1.84005093518623779d4*u533
  u960=-u565
  u565=u960*u1306
  value_1_=value_1_+(c(22)*(x*(u1307*(u62+u1306*(u678+u591)+(u630+u886)&
 &*u1307)+u1306*(u62+u565)+u764)))
  u764=3.48526590647150651d-3*u771
  u566=-3.13673931582435586d-2*u768
  u528=-6.27347863164871172d-2*u768
  u970=2.25359294596864196d4*u533
  u686=3.75598824328106993d3*u533
  u594=(u566+u1306*(u805+u970))
  u825=u686*u1306
  u880=u1306*(u528+u825)+u764
  value_2_=value_2_+(c(22)*((u594*u1307+u880)*z))
  u828=1.63284278170329504d5*u870
  u556=2.0d1*pd2
  u706=-4.40855140631253648d-3*u729*(u2160*(u556+1.1d1)+u828)
  u584=1.63284278170329504d5*u917
  u720=u729*pd2
  u597=5.4555823653117639d-2*u720*(u584-u862)
  u908=-8.90810828518710251d3*u533
  u598=3.47173423247112248d-2*u720*(1.01080743629251597d5*u917-u862)
  u832=pd2*(u584+u859)
  u633=-2.90964392816627408d-1*u729*u832
  u813=9.0926372755196065d-2*u549
  u769=-2.96936942839570083d3*u533
  u909=+u769
  u526=u813*u1306
  u531=u1306*(u526+u633)
  u57=(x*(u1307*(u597+u531+(u526+u908)*u1307)+u1306*(u598+u909*u1306)+u&
 &706))
  value_3_=value_3_+(c(22)*u57)
  u669=7.28912347894906632d-2*u767
  u747=2.90937598295568712d3*u533
  u925=-4.36406397443353068d3*u533
  u883=+u925
  u929=-1.60015679062562792d4*u533
  value_4_=value_4_+(c(22)*(u545*((u747+u804)*u1308+(u883+u914)*u1307+u&
 &1306*(u929+u914)+u669)))
  value_5_=value_5_+(c(22)*u750)
  u750=-5.6339823649216049d3*u533
  u891=+u750
  u570=1.31459588514837448d4*u533
  value_6_=value_6_+(c(22)*(u310*((u891+u961)*u1307+u1306*(u570+u615)+u&
 &566)*z))
  u69=-2.84570769625795685d-3*u796
  u578=1.92085269497412087d-2*u767
  u901=-u634
  u972=-u568
  u568=(u578+u1306*(u939+u901))*u1307
  u78=u1306*(u722+u972)
  value_7_=value_7_+(c(22)*(y*(u568+u1306*(u578+u78)+u69)))
  u765=-7.68341077989648348d-2*u768
  u855=2.07005730208451751d4*u533
  value_1_=value_1_+(c(23)*(u310*((u886+u630)*u1307+u1306*(u855+u678)+u&
 &765)*z))
  value_2_=value_2_+(c(23)*(y*(u594*u1308+u880)))
  u765=7.126486628149682d3*u533
  u594=-1.7816216570374205d3*u533
  u880=+u594
  u718=-5.93873885679140167d2*u533
  u892=+u718
  u915=(u545*((u765+u896)*u1308+(u880+u644)*u1307+u1306*(u892+u644)+u89&
 &9))
  value_3_=value_3_+(c(23)*u915)
  u685=-u862*(u548-2.2d1)
  u675=1.34983768128686413d-3*u729*(u828+u685)
  u949=-u925
  u925=pd2*(5.44280927234431679d4*u917+u574)
  u938=-1.21485391315817772d-2*u729*u925
  u8=pd2*(u605-1.41796308072441282d1*u553)
  u558=u729*u8
  u845=7.12714295719464262d-1*u558
  u81=-4.45446434824665164d-2*u549
  u6=u81*u1306
  u63=u6+u845
  u562=u1306*(u63)
  value_4_=value_4_+(c(23)*(x*(u1308*(u99+(u804+u747)*u1308)+u1307*(u88&
 &3+u562+(u6+u949)*u1307)+u938*u1306+u675)))
  u99=3.30641355473440237d-3*u771
  u823=-3.47173423247112248d-2*u768
  u789=-2.37549554271656067d3*u533
  u910=+u789
  value_5_=value_5_+(c(23)*(z*(u89+u1306*(u823+u1306*(u644+u910))+u99))&
 &)
  u89=4.24539123242856709d5*u870
  u714=1.3d1*pd2
  u717=(u89-u862*(2.2d1+u714))
  u700=1.74263295323575325d-3*u729*u717
  u724=-7.84184828956088966d-2*u768
  u984=-u750
  u750=2.50369226527838572d5*u917
  u660=4.53748185831812103d2*exppd2
  u715=-1.56836965791217793d-2*u720*(u660+u750)
  u579=1.87799412164053497d4*u533
  u77=9.20110199308477719d-1*u558
  u7=u1306*(u615+u77)
  value_6_=value_6_+(c(23)*(x*((u724+u1306*(u615+u579))*u1308+u1307*(u8&
 &91+u7+(u615+u984)*u1307)+u1306*(u715+u825)+u700)))
  u700=1.42285384812897842d-3*u771
  u715=9.20025467593118894d3*u533
  value_7_=value_7_+(c(23)*((u568+u1306*(u62+u1306*(u722+u715))+u700)*z&
 &))
  u568=2.84570769625795685d-3*u771
  u825=2.13428077219346764d-3*u767
  u62=-7.66687889660932411d2*u533
  u911=+u62
  u948=-8.32369501155452378d-2*u768
  u997=6.13350311728745929d3*u533
  u868=6.9001910069483917d3*u533
  u632=(u630+u911)*u1307
  u852=u868*u1306
  value_1_=value_1_+(c(24)*(y*(u1307*(u825+u1306*(u678+u997)+u632)+u130&
 &6*(u948+u852)+u568)))
  u825=-9.41021794747306759d-2*u768
  u948=-u683
  u683=(u948+u805)
  u752=u948*u1306
  u79=u752+u825
  value_2_=value_2_+(c(24)*(u310*(u683*u1307+u79)*z))
  u71=(y*(u1307*(u598+u531+(u526+u909)*u1307)+u1306*(u597+u908*u1306)+u&
 &706))
  value_3_=value_3_+(c(24)*u71)
  u531=(1.13437046457953026d2*u656+u913)
  u526=u729*u531
  u638=-9.89880966277033698d-3*u526
  u831=pd2*(u777-u574)
  u779=u729*u831
  u920=2.69967536257372827d-3*u779
  u716=-9.69791994318562374d2*u533
  u684=+u716
  u788=3.64456173947453316d-2*u767
  u934=-1.45468799147784356d3*u533
  u665=+u934
  u80=3.64456173947453316d-2*u720*(1.12484724961782547d5*u917-u862)
  u581=-4.75142863812976175d-1*u558
  u770=2.96964289883110109d-2*u549
  u894=-1.74562558977341227d4*u533
  u893=+u894
  u996=u883*u1306
  u774=u770*u1306
  value_4_=value_4_+(c(24)*(z*(u1308*(u920+u1306*(u774+u581)+(u774+u684&
 &)*u1308)+u1307*(u788+u1306*(u914+u893)+(u914+u665)*u1307)+u1306*(u80+&
 &u996)+u638)))
  value_5_=value_5_+(c(24)*u57)
  u665=-1.87799412164053497d3*u533
  u920=+u665
  u684=-1.56836965791217793d-2*u768
  u788=-u665
  u665=x**4
  u80=(u961+u920)
  value_6_=value_6_+(c(24)*((u1307*(u577+u82*u665+u80*u1307)+u1306*(u68&
 &4+u788*u1306))*z))
  u581=8.32369501155452378d-2*u767
  u893=-u868
  u877=-2.13428077219346764d-3*u768
  u638=-u997
  u997=-u62
  u62=u997*u1306
  value_7_=value_7_+(c(24)*(x*(u1307*(u581+u1306*(u722+u638)+(u939+u893&
 &)*u1307)+u1306*(u877+u62)+u69)))
  u581=-6.4028423165804029d-3*u768
  u877=-3.20142115829020145d-2*u768
  value_1_=value_1_+(c(25)*((u1307*(u581+u1306*(u678+u715)+u632)+u1306*&
 &(u877+u565)+u700)*z))
  u69=7.08981540362206411d0*exppd2
  u638=u69*(1.6d1*pd2-1.1d1)
  u700=(u638+u913)
  u57=-2.78821272517720521d-2*u729*u700
  u632=1.04557977194145195d-2*u779
  u774=-1.84022039861695544d0*u558
  u770=1.15013774913559715d-1*u549
  u902=-u686
  u546=u770*u1306
  u775=u1306*(u546+u774)
  value_2_=value_2_+(c(25)*(x*(u1307*(u948+u775+(u546+u885)*u1307)+u130&
 &6*(u632+u902*u1306)+u57)))
  u70=-1.10213785157813412d-3*u729*(-u860*(3.4d1*pd2-1.1d1)+5.551665457&
 &79120312d5*u870)
  u670=6.61282710946880473d-3*u779
  u671=1.48788609963048106d-2*u767
  u613=5.11624071600365778d5*u917
  u672=1.48788609963048106d-2*u720*(u613-u862)
  u582=-1.16385757126650963d0*u558
  u837=7.2741098204156852d-2*u549
  u930=-u765
  u931=-7.72036051382882217d3*u533
  u673=u837*u1306
  u580=(z*(u1308*(u670+u1306*(u673+u582)+(u673+u910)*u1308)+u1307*(u671&
 &+u1306*(u644+u930)+(u644+u892)*u1307)+u1306*(u672+u931*u1306)+u70))
  value_3_=value_3_+(c(25)*u580)
  u673=4.49945893762288045d-4*u729*(7.51107679583515716d5*u870-u862*(2.&
 &3d1*pd2-2.2d1))
  u694=-4.0495130438605924d-3*u768
  u998=-u716
  u716=pd2*(-u574+u777)
  u777=u729*u716
  u674=-4.0495130438605924d-3*u777
  u72=-u934
  u934=-4.0495130438605924d-3*u720*(9.47048813387911121d5*u917-u660)
  u932=-u747
  value_4_=value_4_+(c(25)*(y*(u1308*(u694+u1306*(u914+u932)+(u804+u998&
 &)*u1308)+u1307*(u674+u562+(u6+u72)*u1307)+u1306*(u934+u747*u1306)+u67&
 &3)))
  value_5_=value_5_+(c(25)*u915)
  u915=-5.22789885970725977d-3*u777
  u804=-1.56836965791217793d-2*u777
  u933=7.51197648656213986d3*u533
  u682=(u615+u788)
  value_6_=value_6_+(c(25)*(y*((u684+u1306*(u615+u948))*u1308+u1307*(u9&
 &15+u7+u682*u1307)+u1306*(u804+u933*u1306)+u577)))
  u804=2.56113692663216116d-2*u767
  value_7_=value_7_+(c(25)*(u310*((u893+u939)*u1307+u1306*(u997+u722)+u&
 &804)*z))
  u804=(u913+1.13437046457953026d2*u592)
  u915=u729*u804
  u7=7.82569616470938133d-3*u915
  u77=-1.92085269497412087d-2*u768
  u999=-2.13428077219346764d-3*u777
  u635=1.85055515259706771d5*u917
  u820=u720*(-u574+u635)
  u536=-1.92085269497412087d-2*u820
  u676=3.75633415906050304d-1*u558
  u621=u972*u1306
  value_1_=value_1_+(c(26)*(y*((u77+u1306*(u678+u634))*u1308+u1307*(u99&
 &9+u1306*(u722+u676)+(u722+u997)*u1307)+u1306*(u536+u621)+u7)))
  value_2_=value_2_+(c(26)*(u545*(u683*u1308+u79)))
  u7=6.06175818367973767d-3*u915
  u999=-u789
  u536=-1.65320677736720118d-3*u777
  u676=-u718
  u718=-4.95962033210160355d-3*u777
  u789=-u985
  u985=2.90964392816627408d-1*u558
  u683=-1.8185274551039213d-2*u549
  u79=1.18774777135828033d3*u533
  u537=u683*u1306
  u84=u537+u985
  u904=u1306*(u84)
  u797=(y*(u1308*(u823+u1306*(u644+u789)+(u896+u999)*u1308)+u1307*(u536&
 &+u904+(u537+u676)*u1307)+u1306*(u718+u79*u1306)+u7))
  value_3_=value_3_+(c(26)*u797)
  u844=8.16421390851647518d5*u870
  u547=2.5d1*pd2
  u681=-u862*(u547-1.1d1)
  u790=u729*(u844+u681)
  u863=4.49945893762288045d-4*u790
  u816=8.16421390851647518d5*u917
  u612=-2.69967536257372827d-3*u720*(9.64214894892600719d2*exppd2+u816)
  u677=2.42447998579640593d3*u533
  u818=-9.07496371663624206d2*exppd2
  u827=pd2*(u816+u818)
  u975=-4.0495130438605924d-3*u729*u827
  u723=2.72140463617215839d5*u917
  u631=pd2*(u723+u950)
  u858=8.90892869649330328d-2*u729*u631
  u856=-7.42410724707775273d-2*u549
  value_4_=value_4_+(c(26)*(z*(u1308*(u612+u858*u1306+(u856*u1306+u677)&
 &*u1308)+u975*u1306+u863)))
  u975=1.65320677736720118d-3*u729*(u828-u862*(1.1d1+u548))
  u858=-9.42327863099304674d-2*u768
  u856=-u594
  u594=4.89852834510988511d5*u917
  u766=pd2*(u660+u594)
  u687=-1.65320677736720118d-3*u729*u766
  u936=-u587
  value_5_=value_5_+(c(26)*(x*(u1308*(u858+u1306*(u644+u936)+(u896+u765&
 &)*u1308)+u1307*(u880+u904+(u537+u856)*u1307)+u687*u1306+u975)))
  u587=1.91689624855932858d-2*u915
  u896=6.53137112681318014d3*u917
  u915=-4.18231908776580781d-1*u720*(u69+u896)
  u921=-4.70510897373653379d-2*u820
  u721=1.15013774913559715d-1*u720*(3.59225411974724908d5*u917+u950)
  value_6_=value_6_+(c(26)*(z*(u1308*(u915+u1306*(u805+u721)+u682*u1308&
 &)+u1306*(u921+u752)+u587)))
  u587=9.47048813387911121d5*u870
  u915=2.9d1*pd2
  u921=(u587-u862*(1.1d1+u915))
  u721=7.11426924064489211d-4*u729*u921
  u788=1.60018592606922914d6*u917
  u820=1.81499274332724841d3*exppd2
  u682=-2.13428077219346764d-3*u720*(u820+u788)
  u603=7.66687889660932412d3*u533
  u712=1.12690024771815091d0*u558
  u978=3.06675155864372965d3*u533
  u73=u1306*(u678+u712)
  value_7_=value_7_+(c(26)*(x*((u877+u1306*(u722+u603))*u1308+u1307*(u8&
 &93+u73+(u678+u868)*u1307)+u1306*(u682+u978*u1306)+u721)))
  u721=u822*u1307
  u682=u721+u678
  u877=u634*u1306
  u649=u77*u1306
  value_1_=value_1_+(c(27)*(x*(u1307*(u77+u877+u1307*(u682+u888))+u649+&
 &u568)))
  u568=(u805+u686)
  u805=u528+u970*u1306
  u662=u566*u1306+u764
  value_2_=value_2_+(c(27)*((u1307*(u805+u568*u1307)+u662)*z))
  u922=u708*u1307
  u679=(u899+u1307*(u922+u969))*u1308
  u953=u780*u1307
  u887=u953+u644
  u628=u881*u1306
  u600=u563*u1306
  u821=(x*(u679+u1307*(u596+u628+u1307*(u887+u907))+u600+u604))
  value_3_=value_3_+(c(27)*u821)
  u945=u659*u1307
  u762=u914+u945
  u783=u861*u1307
  value_4_=value_4_+(c(27)*(u545*((u747+u783)*u1308+u1307*(u929+u762)+u&
 &996+u669)))
  value_5_=value_5_+(c(27)*u71)
  u71=3.13673931582435586d-2*u767
  u633=-u570
  u570=u576*u1307
  u813=u615+u570
  u706=u984*u1306
  value_6_=value_6_+(c(27)*(u310*(u1307*(u633+u813)+u706+u71)*z))
  u71=-8.53712308877387053d-3*u796
  u633=5.76255808492236261d-2*u767
  u602=-u591
  u591=(u939+u886)
  u840=u591*u1307
  value_7_=value_7_+(c(27)*(y*(u1307*(u633+u1306*(u722+u602)+u840)+u130&
 &6*(u633+u565)+u71)))
  u602=-2.56113692663216116d-2*u768
  u559=u678+u721
  value_1_=value_1_+(c(28)*(u310*(u1307*(u911+u559)+u852+u602)*z))
  u602=u632+u775
  u775=u1306*(u948+u885*u1306)+u57
  u791=(u546+u902)
  value_2_=value_2_+(c(28)*(y*(u1307*(u602+u791*u1307)+u775)))
  u546=u644+u953
  u908=u880*u1306
  u644=(u545*((u765+u922)*u1308+u1307*(u892+u546)+u908+u899))
  value_3_=value_3_+(c(28)*u644)
  value_4_=value_4_+(c(28)*(x*(u1308*(u694+u1307*(u945+u932)+(u783+u998&
 &)*u1308)+u1307*(u934+u562+(u6+u747)*u1307)+u1306*(u674+u72*u1306)+u67&
 &3)))
  value_5_=value_5_+(c(28)*u580)
  u72=1.56836965791217793d-2*u779
  u673=-u933
  u933=5.22789885970725977d-3*u779
  u580=-9.20110199308477719d-1*u558
  u672=u1306*(u961+u580)
  u837=u1307*(u570+u885)
  value_6_=value_6_+(c(28)*(x*((u577+u837)*u1308+u1307*(u72+u672+(u961+&
 &u673)*u1307)+u1306*(u933+u920*u1306)+u684)))
  u72=-1.42285384812897842d-3*u796
  u673=3.20142115829020145d-2*u767
  u933=-u715
  value_7_=value_7_+(c(28)*((u1307*(u673+u1306*(u722+u933)+u840)+u1306*&
 &(u595+u62)+u72)*z))
  u840=(u681+u844)
  u681=-7.11426924064489211d-4*u729
  u62=+u681*u840
  u582=1.14298994719230652d6*u917
  u694=pd2*(u582-u574)
  u674=6.4028423165804029d-3*u729*u694
  u934=-u603
  u952=6.4028423165804029d-3*u779
  u905=-u712
  u619=u1306*(u939+u905)
  u670=u1307*(u721+u888)
  value_1_=value_1_+(c(29)*(x*((u595+u670)*u1308+u1307*(u674+u619+(u939&
 &+u934)*u1307)+u1306*(u952+u886*u1306)+u62)))
  value_2_=value_2_+(c(29)*(z*(u1308*(u602+u791*u1308)+u775)))
  u62=(x*(u1308*(u823+u1307*(u953+u789)+(u922+u999)*u1308)+u1307*(u718+&
 &u904+(u537+u79)*u1307)+u1306*(u536+u676*u1306)+u7))
  value_3_=value_3_+(c(29)*u62)
  u674=3.23961043508847392d-2*u767
  u952=u861*u1308
  u602=u883*u1307
  value_4_=value_4_+(c(29)*(u545*(u1308*(u998+u762+u952)+u602+u996+u674&
 &)))
  value_5_=value_5_+(c(29)*u797)
  u797=u891*u1307
  value_6_=value_6_+(c(29)*(u545*((u813)*u1308+u797+u706)))
  u775=7.11426924064489211d-4*u790
  u762=-6.4028423165804029d-3*u777
  u813=pd2*(-u574+u582)
  u582=-6.4028423165804029d-3*u729*u813
  u791=(u678+u960)
  value_7_=value_7_+(c(29)*(y*((u581+u78)*u1308+u1307*(u762+u73+u791*u1&
 &307)+u1306*(u582+u603*u1306)+u775)))
  u582=-3.84170538994824174d-2*u768
  u762=u886*u1307
  value_1_=value_1_+(c(30)*(u545*((u972+u559)*u1308+u762+u852+u582)))
  value_2_=value_2_+(c(30)*(y*(u1308*(u805+u568*u1308)+u662)))
  u582=-6.94346846494224497d-2*u768
  u775=2.01917121130907657d4*u533
  u603=u708*u1308
  u712=u880*u1307
  u78=(u545*(u1308*(u775+u546+u603)+u712+u908+u582))
  value_3_=value_3_+(c(30)*u78)
  u73=-8.99891787524576089d-4*u796
  u931=-2.42970782631635544d-2*u768
  u546=-u977
  u977=u945+u914
  u568=u952+u977
  u805=u884*u1306
  u662=u564*u1306
  u678=(u1308*(u931+u805+u884*u1307+u1308*(u568+u546))+u564*u1307+u662+&
 &u73)
  value_4_=value_4_+(c(30)*(x*u678))
  u914=5.5106892578906706d-4*u729*(u844-u862*(8.8d1+u547))
  u908=-9.91924066420320709d-3*u720*(7.37340801976694667d2*exppd2+u723)
  u852=-u769
  u769=1.16631627264521074d5*u917
  u996=-3.47173423247112248d-2*u720*(-u574+u769)
  u706=pd2*(u816+u950)
  u952=3.6370549102078426d-2*u729*u706
  u773=-9.0926372755196065d-2*u549
  value_5_=value_5_+(c(30)*(z*(u1308*(u908+u952*u1306+(u773*u1306+u852)&
 &*u1308)+u996*u1306+u914)))
  u996=u570+u615
  u952=u752+u885*u1307
  u773=u684*u1306
  u70=u577*u1307+u773
  value_6_=value_6_+(c(30)*(x*(u1308*(u528+u952+(u996+u686)*u1308)+u70+&
 &u764)))
  u774=7.6199329812820435d4*u917
  u770=pd2*(-u860+u774)
  u559=-3.84170538994824174d-2*u729*u770
  u790=-1.92085269497412087d-2*u777
  u615=2.93911700706593106d5*u917
  u640=pd2*(u615+u950)
  u824=1.40862530964768864d-1*u729*u640
  u601=-9.39083539765125759d-2*u549
  value_7_=value_7_+(c(30)*(z*(u1308*(u559+u1306*(u601*u1306+u824)+u791&
 &*u1308)+u1306*(u790+u715*u1306)+u578)))
  u790=-1.42285384812897842d-2*u796
  u601=1.60071057914510073d-1*u767
  u824=-9.60426347487060436d-2*u768
  u559=-u976
  value_1_=value_1_+(c(31)*(y*(u1307*(u601+u559*u1306+u1307*(u682+u901)&
 &)+u824*u1306+u790)))
  u790=u713*u1307
  value_2_=value_2_+(c(31)*(u310*(u1307*(u993+u790)+u560)*z))
  value_3_=value_3_+(c(31)*(y*((u979+u1307*(u922+u719))*u1308+u1307*(u6&
 &68+u927*u1306+u1307*(u887+u928))+u667*u1306+u68)))
  value_4_=value_4_+(c(31)*(z*((u802+u1307*(u783+u919))*u1308+u1307*(u6&
 &66+u805+u1307*(u977+u926))+u662+u67)))
  value_5_=value_5_+(c(31)*u821)
  u976=-6.97053181294301302d-3*u796
  u824=2.03888055528583131d-1*u767
  u821=-u680
  value_6_=value_6_+(c(31)*((u1307*(u824+u752+u1307*(u996+u821))+u773+u&
 &976)*z))
  u976=2.49710850346635713d-1*u767
  u824=-u935
  u680=u525*u1307
  u821=u680+u722
  u919=u581*u1306
  value_7_=value_7_+(c(31)*(x*(u1307*(u976+u621+u1307*(u821+u824))+u919&
 &+u71)))
  value_1_=value_1_+(c(32)*((u1307*(u633+u877+u1307*(u682+u933))+u649+u&
 &72)*z))
  value_2_=value_2_+(c(32)*(x*((u566+u1307*(u790+u970))*u1308+u1307*(u5&
 &28+u686*u1307)+u764)))
  value_3_=value_3_+(c(32)*(z*(u679+u1307*(u823+u628+u1307*(u887+u910))&
 &+u600+u99)))
  u679=u1306*(u883+u949*u1306)
  value_4_=value_4_+(c(32)*(y*(u1308*(u564+u1307*(u945+u884)+(u783+u747&
 &)*u1308)+u1307*(u938+u1306*(u81*u1307+u63))+u679+u675)))
  value_5_=value_5_+(c(32)*u644)
  u824=(-u862*(u714+2.2d1)+u89)
  u89=-1.74263295323575325d-3*u729*u824
  u675=7.84184828956088966d-2*u767
  u747=1.56836965791217793d-2*u720*(u750+u660)
  u935=-u579
  u714=u1306*(u984+u891*u1306)
  value_6_=value_6_+(c(32)*(y*((u675+u1307*(u570+u935))*u1308+u1307*(u7&
 &47+u672+(u961+u902)*u1307)+u714+u89)))
  u89=7.68341077989648348d-2*u767
  u747=-u855
  u644=u722+u680
  value_7_=value_7_+(c(32)*(u310*(u1307*(u747+u644)+u565+u89)*z))
  u89=(-u862*(u915+1.1d1)+u587)
  u915=+u681*u89
  u747=2.13428077219346764d-3*u720*(u788+u820)
  u788=-u978
  u855=u1306*(u868+u893*u1306)
  value_1_=value_1_+(c(33)*(y*((u673+u1307*(u721+u934))*u1308+u1307*(u7&
 &47+u619+(u939+u788)*u1307)+u855+u915)))
  value_2_=value_2_+(c(33)*(u545*((u948+u790)*u1308+u948*u1307+u825)))
  u747=u1306*(u880+u856*u1306)
  value_3_=value_3_+(c(33)*(y*(u1308*(u858+u1307*(u953+u936)+(u922+u765&
 &)*u1308)+u1307*(u687+u1306*(u683*u1307+u84))+u747+u975)))
  value_4_=value_4_+(c(33)*(z*(u1308*(u612+u562+u1307*(u945+u998)+(u783&
 &+u6+u677)*u1308)+u1307*(u674+u602)+u679+u863)))
  value_5_=value_5_+(c(33)*u62)
  u677=-1.91689624855932858d-2*u526
  u788=4.18231908776580781d-1*u720*(u896+u69)
  u676=9.41021794747306759d-2*u767
  value_6_=value_6_+(c(33)*(z*(u1308*(u788+u672+u837+u80*u1308)+u1307*(&
 &u676+u797)+u714+u677)))
  u677=-7.82569616470938133d-3*u526
  u788=pd2*(u635-u574)
  u676=1.92085269497412087d-2*u729*u788
  u79=2.13428077219346764d-3*u779
  u949=-3.75633415906050304d-1*u558
  value_7_=value_7_+(c(33)*(x*((u578+u1307*(u680+u901))*u1308+u1307*(u6&
 &76+u1306*(u630+u949)+u897*u1307)+u1306*(u79+u911*u1306)+u677)))
  u79=3.84170538994824174d-2*u720*(u774-u860)
  u676=3.84170538994824174d-2*u767
  value_1_=value_1_+(c(34)*(z*(u1308*(u79+u619+u670+u591*u1308)+u1307*(&
 &u676+u762)+u855+u77)))
  value_2_=value_2_+(c(34)*(x*(u1308*(u528+u970*u1307+(u790+u686)*u1308&
 &)+u566*u1307+u764)))
  value_3_=value_3_+(c(34)*(z*(u1308*(u908+u904+u1307*(u953+u775)+(u922&
 &+u537+u852)*u1308)+u1307*(u582+u712)+u747+u914)))
  value_4_=value_4_+(c(34)*(y*u678))
  value_5_=value_5_+(c(34)*u78)
  u852=-3.48526590647150651d-3*u796
  u764=6.27347863164871172d-2*u767
  value_6_=value_6_+(c(34)*(y*(u1308*(u764+u952+(u996+u902)*u1308)+u70+&
 &u852)))
  value_7_=value_7_+(c(34)*(u545*((u888+u644)*u1308+u893*u1307+u565+u67&
 &6)))
  value_1_=value_1_+(c(35)*(y*(u1308*(u877+u888*u1307+(u682)*u1308)+u59&
 &5*u1307+u649)))
  value_2_=value_2_+(c(35)*(u545*(u1308*(u993+u713*u1308)+u560)))
  u852=8.81710281262507296d-3*u771
  u764=-2.57900257269283384d-1*u768
  u546=3.32569375980318493d4*u533
  u79=(u1308*(u764+u628+u881*u1307+u1308*(u603+u887+u546))+u563*u1307+u&
 &600+u852)
  value_3_=value_3_+(c(35)*(y*u79))
  u949=1.79978357504915218d-2*u771
  u856=-2.0247565219302962d-1*u768
  u775=-u894
  u677=6.0742695657908886d-2*u767
  u936=-1.45468799147784356d4*u533
  value_4_=value_4_+(c(35)*(z*(u1308*(u856+u936*u1306+u936*u1307+u1308*&
 &(u568+u775))+u677*u1307+u677*u1306+u949)))
  value_5_=value_5_+(c(35)*(x*u79))
  value_6_=value_6_+(c(35)*(z*(u1308*(u579*u1306+u935*u1307+(u996)*u130&
 &8)+u675*u1307+u724*u1306)))
  value_7_=value_7_+(c(35)*(x*(u1308*(u621+u901*u1307+(u821)*u1308)+u57&
 &8*u1307+u919)))
  if ( lmax .eq. 4 ) go to 100
  u949=u539*p**9
  u856=u949*u586
  u775=5.76255808492236261d-1*u856
  u546=u949*u870
  u852=1.15003183449139862d4*u546
  u79=-2.18506048553365737d5*u546
  u533=+u79
  u764=u949*u847
  u847=-2.3477088494128144d-1*u764
  u586=1.40862530964768864d0*u764
  u870=6.36808684864285064d6*u917
  u884=u949*pd2
  u677=u553*(1.5d1+u575)
  u774=u884*(u870-2.26874092915906051d2*u677)
  u893=2.3477088494128144d-2*u774
  u729=-7.04312654823844319d-2*u774
  u896=u729*u1306
  u928=u893*u1306
  value_1_=value_1_+(c(36)*(u310*((u852+u1306*(u928+u847))*u1307+u1306*&
 &(u533+u1306*(u896+u586))+u775)))
  u775=1.56836965791217793d-1*u856
  u533=-1.69019470947648147d5*u546
  u634=1.72520662370339572d0*u764
  u595=-1.15013774913559715d-1*u774
  u747=u595*u1306
  value_2_=value_2_+(c(36)*((u1306*(u533+u1306*(u747+u634))+u775)*u730)&
 &)
  u855=u949*u599
  u599=-1.48788609963048107d-1*u855
  u975=-3.563243314074841d4*u546
  u765=8.90810828518710251d3*u546
  u714=5.64180191395183158d4*u546
  u579=7.2741098204156852d-1*u764
  u858=-1.8185274551039213d-1*u764
  u767=-3.6370549102078426d-1*u764
  u914=-7.2741098204156852d-2*u774
  u894=1.8185274551039213d-2*u774
  u915=u894*u1306
  u978=u1306*(u915+u858)
  u587=u914*u1306
  u62=(u975+u1306*(u587+u579))*u1308+(u765+u978)*u1307
  u821=(u310*(u62+u1306*(u714+u1306*(u915+u767))+u599))
  value_3_=value_3_+(c(36)*u821)
  u678=-3.64456173947453316d-1*u855
  u673=-1.45468799147784356d4*u546
  u644=2.18203198721676534d4*u546
  u78=1.38195359190395138d5*u546
  u679=2.96964289883110109d-1*u764
  u568=-4.45446434824665164d-1*u764
  u562=-8.90892869649330328d-1*u764
  u984=-2.96964289883110109d-2*u774
  u672=4.45446434824665164d-2*u774
  u891=u672*u1306
  u837=u984*u1306
  value_4_=value_4_+(c(36)*(u756*((u673+u1306*(u837+u679))*u1308+(u644+&
 &u1306*(u891+u568))*u1307+u1306*(u78+u1306*(u891+u562))+u678)))
  u670=u949*u941
  u941=1.10213785157813412d-2*u670
  u790=9.91924066420320709d-2*u856
  u680=-2.47981016605080177d-2*u855
  u681=-4.71163931549652337d-1*u855
  u99=-1.0689729942224523d5*u546
  u934=2.67243248555613075d4*u546
  u880=9.79891911370581276d4*u546
  u682=1.09111647306235278d0*u764
  u887=-2.72779118265588195d-1*u764
  u996=-4.54631863775980325d-1*u764
  u563=u996+u915
  value_5_=value_5_+(c(36)*((u790+u1306*(u1306*(u682+u587)+u99))*u1308+&
 &(u680+u1306*(u1306*(u887+u915)+u934))*u1307+u1306*(u681+u1306*(u1306*&
 &(u563)+u880))+u941))
  u999=4.70510897373653379d-1*u856
  u70=2.81699118246080245d4*u546
  u919=-1.78409441555850822d5*u546
  u976=+u919
  u768=-5.75068874567798575d-1*u764
  u920=1.15013774913559715d0*u764
  u80=5.75068874567798574d-2*u774
  u883=-5.75068874567798574d-2*u774
  u935=u883*u1306
  u63=u80*u1306
  value_6_=value_6_+(c(36)*(x*((u70+u1306*(u63+u768))*u1307+u1306*(u976&
 &+u1306*(u935+u920))+u999)*z))
  u999=u949*u846
  u976=-1.42285384812897842d-2*u999
  u846=-9.60426347487060436d-2*u855
  u84=6.08270020075138276d-1*u856
  u722=1.03502865104225876d5*u546
  u939=-1.26503501794053848d5*u546
  u961=+u939
  u922=-1.05646898223576648d0*u764
  u630=+u922
  u565=5.869272123532036d-1*u764
  u6=7.04312654823844319d-2*u774
  u998=-2.3477088494128144d-2*u774
  u970=u998*u1306
  u945=u6*u1306
  value_7_=value_7_+(c(36)*((u846+u1306*(u1306*(u630+u945)+u722))*u1307&
 &+u1306*(u84+u1306*(u1306*(u565+u970)+u961))+u976))
  u976=-8.53712308877387053d-3*u999
  u846=5.76255808492236261d-2*u856
  u84=2.30006366898279723d3*u546
  u961=2.49710850346635713d-1*u856
  u630=-7.59021010764323088d4*u546
  u565=+u630
  u783=-1.40862530964768864d-1*u764
  u752=-3.22008913657591613d4*u546
  u721=+u752
  u570=9.15606451270997615d-1*u764
  u559=7.04312654823844319d-2*u764
  u621=u559*u1306
  value_1_=value_1_+(c(37)*(u1307*(u846+u1306*(u1306*(u570+u896)+u565)+&
 &(u1306*(u783+u928)+u84)*u1307)+u1306*(u961+u1306*(u621+u721))+u976))
  u976=-5.6339823649216049d4*u546
  u846=-3.75598824328106993d4*u546
  u961=+u846
  u565=1.15013774913559715d-1*u764
  u721=(u976+u1306*(u747+u920))
  u570=u565*u1306
  u877=u1306*(u961+u570)+u775
  value_2_=value_2_+(c(37)*(x*(u721*u1307+u877)*z))
  u649=4.40855140631253648d-3*u949*(u828+u2160*(1.1d1+u556))
  u628=-5.4555823653117639d-2*u884
  u600=+u628*(-u862+u584)
  u601=-4.95962033210160355d-3*u884*(1.92842978978520144d3*exppd2+3.428&
 &96984157691957d6*u917)
  u708=5.4555823653117639d-2*u884
  u901=-1.81499274332724841d3*u553
  u602=u708*(4.4086755105988966d6*u917+u901)
  u797=-5.4555823653117639d-1*u764
  u762=3.6370549102078426d-2*u884
  u578=-4.53748185831812103d2*u553
  u603=u762*(1.46955850353296553d6*u917+u578)
  u712=-1.8185274551039213d-2*u884
  u805=8.0d0*pd2
  u662=+u712*(4.45766079404999545d7*u917+u557*(1.05d2+u805))
  u773=9.0926372755196065d-2*u774
  u719=-9.0926372755196065d-2*u764
  u878=u773*u1306
  u82=u1306*(u662+u878)
  u916=(u1307*(u600+u1306*(u82+u602)+(u1306*(u797+u878)+u765)*u1307)+u1&
 &306*(u601+u1306*(u719*u1306+u603))+u649)
  value_3_=value_3_+(c(37)*u916)
  u683=-7.28912347894906632d-2*u855
  u713=-2.90937598295568712d3*u546
  u861=4.36406397443353068d3*u546
  u81=7.41890875653700216d4*u546
  u580=1.78178573929866066d-1*u764
  u905=-2.67267860894799098d-1*u764
  u526=-7.12714295719464262d-1*u764
  u931=u861+u1306*(u891+u905)
  u659=u1306*(u891+u526)
  value_4_=value_4_+(c(37)*(u730*((u713+u1306*(u837+u580))*u1308+(u931)&
 &*u1307+u1306*(u81+u659)+u683)))
  value_5_=value_5_+(c(37)*u821)
  u821=3.13673931582435586d-2*u856
  u564=5.6339823649216049d3*u546
  u525=-5.07058412842944442d4*u546
  u822=+u525
  u899=-3.45041324740679145d-1*u764
  u566=6.90082649481358289d-1*u764
  value_6_=value_6_+(c(37)*(y*((u564+u1306*(u63+u899))*u1307+u1306*(u82&
 &2+u1306*(u935+u566))+u821)*z))
  u822=-6.4028423165804029d-2*u855
  u528=3.45009550347419585d4*u546
  u985=-u852
  u676=-7.0431265482384432d-1*u764
  u684=2.81725061929537728d-1*u764
  u77=(u528+u1306*(u945+u676))*u1307
  value_7_=value_7_+(c(37)*(u310*(u77+u1306*(u985+u1306*(u970+u684))+u8&
 &22)))
  u822=7.68341077989648348d-2*u856
  u845=-8.97024830903290922d4*u546
  u791=+u845
  u977=9.86037716753382047d-1*u764
  u73=(u928+u783)
  value_1_=value_1_+(c(38)*(y*((u84+u1306*u73)*u1307+u1306*(u791+u1306*&
 &(u896+u977))+u822)*z))
  value_2_=value_2_+(c(38)*(u310*(u721*u1308+u877)))
  u791=1.98384813284064142d-2*u856
  u977=-7.126486628149682d3*u546
  u721=1.7816216570374205d3*u546
  u877=-5.34486497111226151d3*u546
  u581=4.36446589224941112d-1*u764
  u536=-1.09111647306235278d-1*u764
  u718=u1306*(u915+u536)
  u868=(u730*((u977+u1306*(u587+u581))*u1308+(u721+u718)*u1307+u1306*(u&
 &877+u718)+u791))
  value_3_=value_3_+(c(38)*u868)
  u718=-1.34983768128686413d-3*u949*(u685+u828)
  u685=-1.21485391315817772d-2*u855
  u7=-u861
  u823=7.28912347894906632d-2*u949*u925
  u925=3.05484478210347148d4*u546
  u771=u884*(8.81735102119779319d5*u917+u578)
  u908=-1.33633930447399549d-1*u771
  u635=2.67267860894799098d-1*u764
  u576=-4.45446434824665164d-2*u884*(u584+u557)
  u952=-4.8999107830713168d-1*u764
  u675=u553*(2.1d1+u575)
  u671=u884*(8.9153215880999909d6*u917-2.26874092915906051d2*u675)
  u897=4.45446434824665164d-2*u671
  u886=-4.45446434824665164d-2*u774
  u619=u886*u1306
  u904=u897+u619
  u796=u1306*(u904)
  value_4_=value_4_+(c(38)*(u1308*(u685+u1306*(u1306*(u952+u891)+u925)+&
 &(u1306*(u580+u837)+u713)*u1308)+u1307*(u861+u1306*(u796+u908)+(u1306*&
 &(u635+u619)+u7)*u1307)+u1306*(u823+u576*u1306)+u718))
  u993=-2.96936942839570083d3*u546
  value_5_=value_5_+(c(38)*(u756*(u62+u1306*(u993+u978)+u790)))
  u62=-1.74263295323575325d-3*u949*u824
  u824=7.84184828956088966d-2*u856
  u978=-u564
  u715=4.19873858152275866d4*u917
  u591=4.3914350421540982d-1*u884*(u715-u860)
  u720=-8.45097354738240736d4*u546
  u537=+u720
  u953=-1.72520662370339572d-1*u771
  u560=3.45041324740679145d-1*u764
  u825=u884*(1.07767623592417472d6*u917+u557)
  u582=-5.75068874567798574d-2*u825
  u686=8.62603311851697862d-1*u764
  u780=5.75068874567798574d-2*u671
  value_6_=value_6_+(c(38)*((u824+u1306*(u1306*(u686+u935)+u537))*u1308&
 &+u1307*(u564+u1306*(u1306*(u780+u935)+u953)+(u1306*(u560+u935)+u978)*&
 &u1307)+u1306*(u591+u1306*(u570+u582))+u62))
  u780=1.28056846331608058d-1*u856
  u591=-5.75015917245699309d4*u546
  u953=+u591
  u582=4.22587592894306592d-1*u764
  value_7_=value_7_+(c(38)*(x*(u77+u1306*(u953+u1306*(u970+u582))+u780)&
 &*z))
  u780=1.92085269497412087d-1*u856
  u62=-7.04312654823844319d-2*u764
  u77=4.6954176988256288d-1*u764
  u570=2.11293796447153296d-1*u764
  u927=(u928+u62)
  u863=u927*u1307
  u802=u570*u1306
  value_1_=value_1_+(c(39)*(u310*(u1307*(u985+u1306*(u896+u77)+u863)+u1&
 &306*(u953+u802)+u780)))
  u780=9.41021794747306759d-2*u856
  u77=-1.12679647298432098d4*u546
  u953=+u77
  u948=-6.76077883790592589d4*u546
  u979=+u948
  u938=(u953+u1306*(u747+u566))
  u724=u560*u1306
  u612=u1306*(u979+u724)+u780
  value_2_=value_2_+(c(39)*(y*(u938*u1307+u612)*z))
  u687=-8.92731659778288639d-2*u884*(-u862+u723)
  u789=2.44926417255494255d6*u917
  u972=-9.07496371663624206d2*u553
  u604=u708*(u789+u972)
  u60=1.06134780810714177d7*u917
  u951=u553*(2.5d1+u575)
  u806=u884*(u60-2.26874092915906051d2*u951)
  u74=-7.2741098204156852d-2*u806
  u554=u887*u1306
  u995=(u310*(u1307*(u604+u1306*(u878+u74)+(u878+u887)*u1307)+u1306*(u6&
 &04+u554)+u687))
  value_3_=value_3_+(c(39)*u995)
  u688=-1.45782469578981326d-1*u949*u770
  u770=u884*(u594+u557)
  u689=8.90892869649330328d-2*u770
  u555=-8.90892869649330328d-2*u764
  u697=4.80047037187688375d4*u546
  u624=-1.33633930447399549d-1*u764
  u690=4.45446434824665164d-2*u884*(1.40424479226483373d6*u917+u578)
  u926=u553*(1.0d1+pd2)
  u527=u884*(2.12269561621428355d6*u917-1.13437046457953026d2*u926)
  u833=-1.18785715953244044d-1*u527
  u834=2.96964289883110109d-2*u774
  u838=u659+(u891+u624)*u1307
  u659=u834*u1306
  u873=u624*u1306
  value_4_=value_4_+(c(39)*(u756*(u1308*(u689+u1306*(u659+u833)+(u659+u&
 &555)*u1308)+u1307*(u697+u838)+u1306*(u690+u873)+u688)))
  value_5_=value_5_+(c(39)*u916)
  u916=-1.72520662370339572d-1*u764
  u532=-1.31459588514837448d4*u546
  u889=+u532
  u691=2.3002754982711943d-1*u764
  u636=5.75068874567798574d-2*u764
  u753=(u63+u916)
  value_6_=value_6_+(c(39)*(x*(u1307*(u564+u1306*(u935+u691)+u753*u1307&
 &)+u1306*(u889+u636*u1306)+u821)*z))
  u889=2.84570769625795685d-3*u670
  u779=-8.32369501155452378d-2*u855
  u573=6.9001910069483917d3*u546
  u692=-1.92085269497412087d-2*u855
  u898=4.83013370486387419d4*u546
  u550=-4.22587592894306592d-1*u764
  u618=-4.60012733796559447d3*u546
  u693=+u618
  u637=2.3477088494128144d-2*u764
  u76=u637*u1306
  value_7_=value_7_+(c(39)*(u1307*(u779+u1306*(u1306*(u62+u970)+u898)+(&
 &u1306*(u550+u945)+u573)*u1307)+u1306*(u692+u1306*(u76+u693))+u889))
  u889=-2.07005730208451751d4*u546
  u779=+u889
  u693=5.63450123859075456d-1*u764
  value_1_=value_1_+(c(40)*(x*(u1307*(u779+u1306*(u896+u693)+u863)+u130&
 &6*(u779+u621)+u822)*z))
  u779=u949*(u913+u638)
  u913=2.78821272517720521d-2*u779
  u558=-u77
  u77=6.85793968315383915d5*u917
  u863=u884*(u77+u862)
  u638=-3.13673931582435586d-2*u863
  u639=3.45041324740679145d-1*u771
  u960=-6.90082649481358289d-1*u764
  u782=u949*u640
  u640=2.3002754982711943d-1*u782
  u58=-1.15013774913559715d-1*u671
  u895=1.15013774913559715d-1*u774
  u990=-1.15013774913559715d-1*u764
  u843=u895*u1306
  u991=u1306*(u58+u843)
  value_2_=value_2_+(c(40)*(u1307*(u953+u1306*(u991+u639)+(u1306*(u960+&
 &u843)+u558)*u1307)+u1306*(u638+u1306*(u990*u1306+u640))+u913))
  u857=-3.54490770181103205d0*exppd2
  u864=u884*(u605+u857)
  u605=-6.34831402509005254d-1*u864
  u606=2.18223294612470556d-1*u770
  u710=-2.18223294612470556d-1*u764
  u875=1.95978382274116255d4*u546
  u674=-5.4555823653117639d-2*u764
  u607=u708*(2.10092437912490628d6*u917+u972)
  u708=-2.90964392816627408d-1*u527
  u632=7.2741098204156852d-2*u774
  u942=-2.90964392816627408d-1*u764
  u709=-2.36408569163509769d-1*u764
  u622=u1306*(u915+u942)
  u836=u632*u1306
  u751=(u756*(u1308*(u606+u1306*(u836+u708)+(u836+u710)*u1308)+u1307*(u&
 &875+u622+(u915+u674)*u1307)+u1306*(u607+u709*u1306)+u605))
  value_3_=value_3_+(c(40)*u751)
  u705=3.64456173947453316d-2*u884
  u826=u705*(u615+u574)
  u608=8.90892869649330328d-2*u764
  u888=-1.33633930447399549d-1*u770
  u609=1.33633930447399549d-1*u764
  u658=4.02767886153479442d5*u917
  u620=-1.33633930447399549d-1*u884
  u68=+u620*(u658+u557)
  u95=1.78178573929866066d-1*u527
  u803=u1306*(u619+u95)
  value_4_=value_4_+(c(40)*(u310*(u1308*(u931+(u837+u608)*u1308)+u1307*&
 &(u888+u803+(u619+u609)*u1307)+u1306*(u68+u608*u1306)+u826)))
  value_5_=value_5_+(c(40)*u868)
  u868=u949*u694
  u694=1.56836965791217793d-2*u868
  u931=-u70
  u616=-1.72520662370339572d-1*u770
  u583=1.72520662370339572d-1*u764
  u792=u884*(5.98709019957874846d5*u917+u557)
  u571=-1.72520662370339572d-1*u792
  u610=5.75068874567798575d-1*u764
  u655=2.3002754982711943d-1*u527
  u912=(u935+u583)
  value_6_=value_6_+(c(40)*(u310*((u931+u1306*(u935+u610))*u1308+u1307*&
 &(u616+u1306*(u935+u655)+u912*u1307)+u1306*(u571+u691*u1306)+u694)))
  u694=-2.56113692663216116d-2*u855
  u571=1.40862530964768864d-1*u764
  u821=u1306*(u970+u571)
  value_7_=value_7_+(c(40)*(y*((u573+u1306*(u945+u550))*u1307+u1306*(u5&
 &73+u821)+u694)*z))
  u865=u884*(u594+u660)
  u594=1.92085269497412087d-2*u865
  u986=-u528
  u997=-7.04312654823844319d-2*u770
  u781=u884*(u816+u557)
  u793=-7.04312654823844319d-2*u781
  u611=7.0431265482384432d-1*u764
  u923=9.39083539765125759d-2*u527
  u727=u571*u1306
  value_1_=value_1_+(c(41)*(u310*((u986+u1306*(u896+u611))*u1308+u1307*&
 &(u997+u1306*(u970+u923)+(u970+u559)*u1307)+u1306*(u793+u727)+u594)))
  value_2_=value_2_+(c(41)*(u730*(u938*u1308+u612)))
  u594=4.95962033210160355d-3*u884
  u793=u594*(1.20830365846043833d6*u917+u660)
  u938=-3.38508114837109895d4*u546
  u612=2.18223294612470556d-1*u764
  u65=-5.4555823653117639d-2*u770
  u567=5.4555823653117639d-2*u764
  u955=+u628*(u613+u557)
  u613=2.54593843714548982d-1*u764
  u911=7.2741098204156852d-2*u527
  u668=-1.8185274551039213d-2*u774
  u614=3.6370549102078426d-2*u764
  u56=u668*u1306
  u643=u1306*(u56+u911)
  u842=(u310*(u1308*(u938+u1306*(u915+u613)+(u587+u612)*u1308)+u1307*(u&
 &65+u643+(u56+u567)*u1307)+u1306*(u955+u614*u1306)+u793))
  value_3_=value_3_+(c(41)*u842)
  u869=u949*u827
  u827=1.21485391315817772d-2*u869
  u799=u949*u631
  u631=-2.67267860894799098d-1*u799
  u695=2.22723217412332582d-1*u764
  u798=+u620*(u723+u557)
  u543=3.0d0*pd2
  u691=u553*(2.5d1+u543)
  u937=2.96964289883110109d-2*u884*(u60-2.26874092915906051d2*u691)
  u59=-7.42410724707775273d-2*u774
  u867=u59*u1306
  value_4_=value_4_+(c(41)*(u756*(u1308*(u631+u937*u1306+(u867+u695)*u1&
 &308)+u798*u1306+u827)))
  u728=-1.65320677736720118d-3*u949*(-u862*(u548+1.1d1)+u828)
  u828=9.42327863099304674d-2*u856
  u819=-u721
  u829=9.91924066420320709d-3*u865
  u551=-5.87935146822348766d4*u546
  u549=-5.4555823653117639d-2*u771
  u641=1.09111647306235278d-1*u764
  u664=-1.8185274551039213d-2*u770
  u696=1.63667470959352917d-1*u764
  u879=1.8185274551039213d-2*u671
  u524=u879+u56
  u741=u1306*(u524)
  value_5_=value_5_+(c(41)*(u1308*(u828+u1306*(u1306*(u696+u915)+u551)+&
 &(u1306*(u581+u587)+u977)*u1308)+u1307*(u721+u1306*(u741+u549)+(u1306*&
 &(u641+u56)+u819)*u1307)+u1306*(u829+u664*u1306)+u728))
  u617=4.70510897373653379d-2*u865
  u830=u949*u832
  u832=-6.90082649481358289d-1*u830
  u962=-1.72520662370339572d-1*u781
  u865=u553*(1.5d1+pd2)
  u754=u884*(u870-2.26874092915906051d2*u865)
  u870=1.15013774913559715d-1*u754
  value_6_=value_6_+(c(41)*(u756*(u1308*(u832+u1306*(u747+u870)+u912*u1&
 &308)+u1306*(u962+u724)+u617)))
  u617=-7.11426924064489211d-4*u949*u89
  u832=3.20142115829020145d-2*u856
  u962=-u573
  u870=u949*u831
  u89=3.84170538994824174d-2*u870
  u777=-2.11293796447153296d-1*u771
  u912=7.51107679583515716d5*u917
  u831=u884*(u912+u557)
  u585=-7.04312654823844319d-2*u831
  u589=3.5215632741192216d-1*u764
  u98=7.04312654823844319d-2*u671
  u83=9.39083539765125759d-2*u764
  value_7_=value_7_+(c(41)*((u832+u1306*(u1306*(u589+u970)+u986))*u1308&
 &+u1307*(u573+u1306*(u1306*(u98+u896)+u777)+(u1306*(u582+u896)+u962)*u&
 &1307)+u1306*(u89+u1306*(u83*u1306+u585))+u617))
  u617=-2.84570769625795685d-3*u999
  u832=1.92085269497412087d-2*u856
  u89=-u618
  u777=-2.3477088494128144d-2*u764
  u585=8.32369501155452378d-2*u856
  u589=-u898
  u98=(u559+u896)
  u83=u582*u1306
  u898=u962*u1306
  value_1_=value_1_+(c(42)*(u1307*(u832+u1306*(u83+u589)+u1307*((u777+u&
 &928)*u1307+u1306*u98+u89))+u1306*(u585+u898)+u617))
  u617=(u747+u560)
  u89=u566*u1306
  u585=u979+u89
  u589=u953*u1306
  u618=u589+u780
  value_2_=value_2_+(c(42)*(x*(u1307*(u585+u617*u1307)+u618)*z))
  u725=u765*u1306
  u760=(u1307*(u601+u1306*(u797*u1306+u602)+u1307*((u719+u878)*u1307+u8&
 &2+u603))+u1306*(u600+u725)+u649)
  value_3_=value_3_+(c(42)*u760)
  u878=u834*u1307
  value_4_=value_4_+(c(42)*(u730*(u1308*(u689+u1307*(u878+u833)+(u878+u&
 &555)*u1308)+u1307*(u690+u838)+u1306*(u697+u873)+u688)))
  value_5_=value_5_+(c(42)*u995)
  u697=-3.13673931582435586d-2*u855
  u833=-u532
  u74=-5.75068874567798574d-2*u764
  u995=-2.3002754982711943d-1*u764
  u838=u583*u1306
  value_6_=value_6_+(c(42)*(y*(u1307*(u833+u1306*(u935+u995)+(u63+u74)*&
 &u1307)+u1306*(u978+u838)+u697)*z))
  u833=-1.92085269497412087d-1*u855
  u878=-u591
  u591=-2.11293796447153296d-1*u764
  u82=-4.6954176988256288d-1*u764
  u532=(u945+u591)
  u555=u532*u1307
  value_7_=value_7_+(c(42)*(u310*(u1307*(u878+u1306*(u970+u82)+u555)+u1&
 &306*(u852+u621)+u833)))
  u833=2.56113692663216116d-2*u856
  u82=7.66687889660932411d2*u546
  u701=-2.99008276967763641d4*u546
  u866=+u701
  u642=1.87816707953025152d-1*u764
  value_1_=value_1_+(c(43)*(y*(u1307*(u82+u1306*(u896+u642)+(u928+u777)&
 &*u1307)+u1306*(u866+u802)+u833)*z))
  u866=u884*(u584+u860)
  u642=-1.88204358949461352d-1*u866
  u584=3.45041324740679145d-1*u770
  u711=-4.6005509965423886d-1*u527
  u75=u584+u1306*(u843+u711)
  u604=u899*u1306
  u814=u1306*(u584+u604)+u642
  u968=(u843+u899)
  value_2_=value_2_+(c(43)*(u310*(u1307*(u75+u968*u1307)+u814)))
  u872=u632*u1307
  u874=u674*u1306
  u707=(u730*(u1308*(u606+u1307*(u872+u708)+(u872+u710)*u1308)+u1307*(u&
 &607+u622+(u915+u709)*u1307)+u1306*(u875+u874)+u605))
  value_3_=value_3_+(c(43)*u707)
  u622=u949*u700
  u700=-1.43982686003932174d-2*u622
  u872=-8.0990260877211848d-3*u863
  u817=5.93928579766220219d-2*u782
  u835=-2.96964289883110109d-2*u764
  u903=3.8875325221061687d-1*u864
  u794=u884*(u77+u557)
  u77=-4.45446434824665164d-2*u794
  u698=4.45446434824665164d-2*u764
  u661=2.42970782631635544d-2*u884*(u658+u860)
  u658=8.90892869649330328d-2*u771
  u726=-2.96964289883110109d-2*u671
  u645=-2.67267860894799098d-1*u771
  u61=u553*(2.7d1+u575)
  u940=u884*(1.14625563275571312d7*u917-2.26874092915906051d2*u61)
  u702=4.45446434824665164d-2*u940
  u800=u884*(6.20480257047252114d5*u917+u557)
  u871=-4.45446434824665164d-2*u800
  u699=-1.78178573929866066d-1*u764
  u801=-8.90892869649330328d-2*u774
  value_4_=value_4_+(c(43)*(u1308*(u872+u1306*(u699*u1306+u658)+u1308*(&
 &(u835+u659)*u1308+u1306*(u726+u659)+u817))+u1307*(u903+u1306*(u1306*(&
 &u702+u619)+u645)+u1307*((u698+u619)*u1307+u1306*(u702+u801*u1306)+u77&
 &))+u1306*(u661+u1306*(u698*u1306+u871))+u700))
  value_5_=value_5_+(c(43)*u751)
  u645=3.13673931582435586d-2*u866
  u871=-5.75068874567798574d-2*u770
  u699=-3.13673931582435586d-2*u866
  u726=5.75068874567798574d-2*u770
  u903=(u636+u935)
  value_6_=value_6_+(c(43)*(u1307*(u645+u665*(u63+u883)+u1307*(u903*u13&
 &07+u63+u871))+u1306*(u699+u1306*(u74*u1306+u726))))
  u726=-u701
  u645=-u82
  u871=-1.87816707953025152d-1*u764
  value_7_=value_7_+(c(43)*(x*(u1307*(u726+u1306*(u970+u871)+u555)+u130&
 &6*(u645+u76)+u694)*z))
  u726=1.13828307850318274d-2*u779
  u645=9.79705669021977021d4*u917
  u871=u884*(u2160+u645)
  u699=-2.56113692663216116d-2*u871
  u801=u884*(u615+u557)
  u700=-2.3477088494128144d-2*u801
  u872=u884*(1.41513041080952236d5*u917-u2160)
  u817=-7.68341077989648348d-2*u872
  u835=1.40862530964768864d-1*u771
  u77=1.27361736972857013d6*u917
  u661=u884*(u77+u557*(3.0d0+u575))
  u658=2.3477088494128144d-2*u661
  u834=u884*(5.55166545779120312d5*u917+u557)
  u615=7.04312654823844319d-2*u834
  u82=u884*(8.06624334161427748d6*u917+u557*(1.9d1+u575))
  u701=-7.04312654823844319d-2*u82
  u702=4.6954176988256288d-2*u774
  u751=u62*u1306
  value_1_=value_1_+(c(44)*(u1307*(u699+u1306*(u1306*(u701+u945)+u835)+&
 &u1307*((u637+u970)*u1307+u1306*(u658+u702*u1306)+u700))+u1306*(u817+u&
 &1306*(u751+u615))+u726))
  value_2_=value_2_+(c(44)*(u756*(u1308*(u75+u968*u1308)+u814)))
  u726=8.81710281262507296d-3*u779
  u699=-1.98384813284064142d-2*u863
  u700=1.45482196408313704d-1*u782
  u817=-7.2741098204156852d-2*u764
  u835=1.58707850627251313d-1*u864
  u658=-1.8185274551039213d-2*u794
  u615=1.8185274551039213d-2*u764
  u701=-1.98384813284064142d-2*u871
  u702=2.18223294612470556d-1*u771
  u75=-7.2741098204156852d-2*u671
  u871=-1.09111647306235278d-1*u771
  u555=1.8185274551039213d-2*u940
  u968=-1.8185274551039213d-2*u801
  u659=-4.36446589224941112d-1*u764
  u76=-3.6370549102078426d-2*u774
  u779=(u1308*(u699+u1306*(u659*u1306+u702)+u1308*((u817+u836)*u1308+u1&
 &306*(u75+u836)+u700))+u1307*(u835+u1306*(u1306*(u555+u56)+u871)+u1307&
 &*((u615+u56)*u1307+u1306*(u555+u76*u1306)+u658))+u1306*(u701+u1306*(u&
 &615*u1306+u968))+u726)
  value_3_=value_3_+(c(44)*u779)
  u836=4.0495130438605924d-3*u884*(2.97177386269999696d6*u917+u818)
  u814=-2.67267860894799098d-1*u884*(2.64883384587423417d5*u917+u950)
  u703=1.6333035943571056d-1*u764
  u876=1.16375039318227485d4*u546
  u704=2.96964289883110109d-2*u764
  u590=u886*u1307
  value_4_=value_4_+(c(44)*(u730*(u1308*(u814+u1306*(u891+u704)+u1307*(&
 &u590+u95)+(u590+u837+u703)*u1308)+u1307*(u888+u609*u1307)+u1306*(u876&
 &+u873)+u836)))
  value_5_=value_5_+(c(44)*u842)
  u842=1.69019470947648147d4*u546
  u837=u949*u8
  u8=-2.76033059792543316d0*u837
  u873=+u8
  u786=-3.38038941895296294d4*u546
  u839=+u786
  u795=u883*u1307
  value_6_=value_6_+(c(44)*(u730*(u1308*(u873+u1306*(u935+u560)+u1307*(&
 &u795+u655)+(u795+u583)*u1308)+u1307*(u616+u583*u1307)+u1306*(u839+u83&
 &8)+u842)))
  u873=u884*(3.10240128523626057d6*u917-u660)
  u839=6.4028423165804029d-3*u873
  u655=-2.11293796447153296d-1*u770
  u838=u884*(5.26138229659950623d5*u917+u557)
  u795=-2.11293796447153296d-1*u838
  u616=2.3477088494128144d-1*u764
  u698=2.81725061929537728d-1*u527
  u623=(u896+u570)
  u850=u616*u1306
  value_7_=value_7_+(c(44)*(u310*((u985+u1306*(u970+u616))*u1308+u1307*&
 &(u655+u1306*(u896+u698)+u623*u1307)+u1306*(u795+u850)+u839)))
  u839=1.92085269497412087d-2*u870
  u795=-1.40862530964768864d-1*u782
  u980=-4.14011460416903502d4*u546
  u694=+u980
  u841=u998*u1307
  value_1_=value_1_+(c(45)*(u730*(u1308*(u795+u1306*(u896+u582)+u1307*(&
 &u841+u923)+(u841+u559)*u1308)+u1307*(u997+u559*u1307)+u1306*(u694+u80&
 &2)+u839)))
  value_2_=value_2_+(c(45)*(u310*(u1308*(u585+u617*u1308)+u618)))
  u795=3.47173423247112248d-2*u884*(2.19267459257299619d5*u917-u574)
  u694=-1.09111647306235278d-1*u884*(4.46310360332233976d5*u917+u950)
  u617=1.27296921857274491d-1*u764
  u585=-2.4942703198523887d4*u546
  u618=6.18299334735333242d-1*u764
  u687=u668*u1307
  u802=(u730*(u1308*(u694+u1306*(u915+u618)+u1307*(u687+u911)+(u687+u58&
 &7+u617)*u1308)+u1307*(u65+u567*u1307)+u1306*(u585+u874)+u795))
  value_3_=value_3_+(c(45)*u802)
  u587=-4.49945893762288045d-4*u949*u840
  u840=u705*(u723-u862)
  u583=-4.45446434824665164d-2*u781
  u705=7.42410724707775273d-2*u764
  u874=4.0495130438605924d-3*u869
  u841=+u620*(u816+u578)
  u816=4.45446434824665164d-2*u806
  value_4_=value_4_+(c(45)*(u1308*(u840+u841*u1306+u1308*((u705+u867)*u&
 &1308+u816*u1306+u583))+u874*u1306+u587))
  u874=u884*(u769-u574)
  u841=1.04152026974133675d-1*u874
  u816=-1.09111647306235278d-1*u949*u706
  u706=2.72779118265588195d-1*u764
  u769=-5.4555823653117639d-2*u781
  u869=u553*(2.5d1+pd2)
  u620=u762*(u60-2.26874092915906051d2*u869)
  u60=-9.0926372755196065d-2*u774
  u762=u60*u1306
  value_5_=value_5_+(c(45)*(u756*(u1308*(u816+u620*u1306+(u762+u706)*u1&
 &308)+u769*u1306+u841)))
  u806=-1.91689624855932858d-2*u949*u531
  u997=1.76347020423955864d5*u917
  u867=7.84184828956088966d-2*u884*(u997-u862)
  u923=-5.75068874567798574d-2*u794
  u531=4.70510897373653379d-2*u949*u788
  u788=-1.72520662370339572d-1*u884*(u77+u578)
  u77=5.75068874567798574d-2*u940
  value_6_=value_6_+(c(45)*(u1308*(u867+u1306*(u89+u788)+u1308*(u903*u1&
 &308+u1306*(u77+u747)+u923))+u1306*(u531+u589)+u806))
  u806=1.92085269497412087d-2*u868
  u867=-4.22587592894306592d-1*u799
  u923=-2.11293796447153296d-1*u792
  u531=u884*(1.48588693134999848d7*u917+u557*(3.5d1+u543))
  u788=4.6954176988256288d-2*u531
  u77=-9.39083539765125759d-2*u774
  u543=u77*u1306
  value_7_=value_7_+(c(45)*(u756*(u1308*(u867+u1306*(u543+u788)+u623*u1&
 &308)+u1306*(u923+u684*u1306)+u806)))
  u806=6.4028423165804029d-2*u856
  u867=-2.81725061929537728d-1*u764
  u923=u893*u1307
  u788=u923+u896
  u903=u611*u1306
  u623=u986*u1306
  value_1_=value_1_+(c(46)*(u310*(u1307*(u852+u903+u1307*(u788+u867))+u&
 &623+u806)))
  u806=(u747+u565)
  u940=u961+u920*u1306
  u747=u976*u1306+u775
  value_2_=value_2_+(c(46)*(y*(u1307*(u940+u806*u1307)+u747)*z))
  u89=u914*u1307
  u589=(u975+u1307*(u89+u579))*u1308
  u723=u894*u1307
  u637=u723+u915
  u915=u858*u1306
  u966=(u310*(u589+u1307*(u714+u915+u1307*(u637+u767))+u725+u599))
  value_3_=value_3_+(c(46)*u966)
  u572=u672*u1307
  u627=u572+u891
  u743=u984*u1307
  u776=u905*u1306
  u785=u861*u1306
  value_4_=value_4_+(c(46)*(u756*((u713+u1307*(u743+u580))*u1308+u1307*&
 &(u81+u776+u1307*(u627+u526))+u785+u683)))
  value_5_=value_5_+(c(46)*u760)
  u868=-u525
  u81=u80*u1307
  u525=u81+u935
  u760=u978*u1306
  value_6_=value_6_+(c(46)*(x*(u1307*(u868+u724+u1307*(u525+u960))+u760&
 &+u697)*z))
  u868=8.53712308877387053d-3*u670
  u649=-2.49710850346635713d-1*u855
  u600=-u752
  u866=-5.76255808492236261d-2*u855
  u601=-u630
  u684=-9.15606451270997615d-1*u764
  u719=-u84
  u636=(u62+u945)
  u752=u719*u1306
  value_7_=value_7_+(c(46)*(u1307*(u649+u1306*(u727+u601)+u1307*(u636*u&
 &1307+u1306*(u684+u970)+u600))+u1306*(u866+u752)+u868))
  value_1_=value_1_+(c(47)*(x*(u1307*(u962+u83+u1307*(u788+u783))+u898+&
 &u833)*z))
  u868=u638+u1306*(u960*u1306+u639)
  u649=u991+u640
  u600=u1306*(u953+u558*u1306)+u913
  u866=(u990+u843)
  value_2_=value_2_+(c(47)*(u1307*(u868+u1307*(u866*u1307+u649))+u600))
  u601=u536*u1306
  u833=u721*u1306
  u684=(u756*((u977+u1307*(u89+u581))*u1308+u1307*(u877+u601+u1307*(u63&
 &7+u536))+u833+u791))
  value_3_=value_3_+(c(47)*u684)
  u991=u1306*(u888+u609*u1306)
  value_4_=value_4_+(c(47)*(u310*(u1308*(u861+u1307*(u572+u905)+(u743+u&
 &608)*u1308)+u1307*(u68+u803+(u619+u608)*u1307)+u991+u826)))
  value_5_=value_5_+(c(47)*u707)
  u875=u949*u813
  u708=-1.56836965791217793d-2*u875
  u813=1.72520662370339572d-1*u792
  u707=1.72520662370339572d-1*u770
  u710=-2.3002754982711943d-1*u527
  u662=u1306*(u63+u710)
  u942=u1306*(u707+u916*u1306)
  value_6_=value_6_+(c(47)*(u310*((u70+u1307*(u81+u768))*u1308+u1307*(u&
 &813+u662+(u63+u995)*u1307)+u942+u708)))
  u708=-7.68341077989648348d-2*u855
  u813=-u889
  u889=-5.63450123859075456d-1*u764
  value_7_=value_7_+(c(47)*(y*(u1307*(u813+u1306*(u970+u889)+(u945+u62)&
 &*u1307)+u1306*(u813+u621)+u708)*z))
  u889=-6.4028423165804029d-3*u873
  u813=2.11293796447153296d-1*u838
  u709=2.11293796447153296d-1*u770
  u526=-2.81725061929537728d-1*u527
  u843=u1306*(u945+u526)
  u797=u1306*(u709+u591*u1306)
  value_1_=value_1_+(c(48)*(u310*((u852+u1307*(u923+u847))*u1308+u1307*&
 &(u813+u843+(u945+u847)*u1307)+u797+u889)))
  u889=u899*u1307
  u838=u895*u1307
  value_2_=value_2_+(c(48)*(u730*(u1308*(u584+u1307*(u838+u711)+(u838+u&
 &899)*u1308)+u1307*(u584+u889)+u642)))
  u711=u1306*(u65+u567*u1306)
  u813=(u310*(u1308*(u938+u1307*(u723+u613)+(u89+u612)*u1308)+u1307*(u9&
 &55+u643+(u56+u614)*u1307)+u711+u793))
  value_3_=value_3_+(c(48)*u813)
  value_4_=value_4_+(c(48)*(u756*(u1308*(u814+u803+u1307*(u572+u704)+(u&
 &743+u619+u703)*u1308)+u1307*(u876+u624*u1307)+u991+u836)))
  value_5_=value_5_+(c(48)*u779)
  u726=-u842
  u876=-u8
  u555=-u786
  value_6_=value_6_+(c(48)*(u756*(u1308*(u876+u662+u1307*(u81+u899)+u75&
 &3*u1308)+u1307*(u555+u916*u1307)+u942+u726)))
  u707=-1.13828307850318274d-2*u622
  u876=7.68341077989648348d-2*u872
  u555=-7.04312654823844319d-2*u834
  u842=2.56113692663216116d-2*u884*(u645+u2160)
  u645=-1.40862530964768864d-1*u771
  u726=7.04312654823844319d-2*u82
  u872=2.3477088494128144d-2*u801
  u710=-2.3477088494128144d-2*u661
  u609=-4.6954176988256288d-2*u774
  value_7_=value_7_+(c(48)*(u1307*(u876+u1306*(u1306*(u710+u928)+u645)+&
 &u1307*(u98*u1307+u1306*(u726+u609*u1306)+u555))+u1306*(u842+u1306*(u7&
 &77*u1306+u872))+u707))
  u609=1.01236252465604292d6*u917
  u876=u884*(u609+u574)
  u661=-1.92085269497412087d-2*u876
  u842=u884*(u750+u950)
  u710=4.22587592894306592d-1*u842
  u872=1.38003820138967834d4*u546
  value_1_=value_1_+(c(49)*(u756*(u1308*(u710+u843+u1307*(u923+u783)+u5&
 &32*u1308)+u1307*(u872+u62*u1307)+u797+u661)))
  value_2_=value_2_+(c(49)*(u1308*(u868+u1308*(u866*u1308+u649))+u600))
  u661=(u756*(u1308*(u694+u643+u1307*(u723+u618)+(u89+u56+u617)*u1308)+&
 &u1307*(u585+u674*u1307)+u711+u795))
  value_3_=value_3_+(c(49)*u661)
  u710=-3.23961043508847392d-2*u855
  u703=8.72812794886706136d3*u546
  u82=u984*u1308
  u645=u82+u627
  u726=u905*u1307
  u555=u861*u1307
  value_4_=value_4_+(c(49)*(u310*(u1308*(u703+u776+u726+u1308*(u645+u58&
 &0))+u555+u785+u710)))
  value_5_=value_5_+(c(49)*u802)
  u776=u564*u1307
  value_6_=value_6_+(c(49)*(u310*(u1308*(u724+u889+(u525)*u1308)+u776+u&
 &760)))
  u802=1.92085269497412087d-2*u876
  u662=-4.22587592894306592d-1*u842
  u866=-u872
  u724=u729*u1307
  value_7_=value_7_+(c(49)*(u730*(u1308*(u662+u821+u1307*(u724+u698)+(u&
 &724+u570)*u1308)+u1307*(u655+u570*u1307)+u1306*(u866+u621)+u802)))
  u802=3.84170538994824174d-2*u856
  u662=-2.76007640277935668d4*u546
  u655=+u662
  u724=u783*u1307
  u866=u84*u1307
  value_1_=value_1_+(c(50)*(u310*(u1308*(u655+u83+u724+(u788+u571)*u130&
 &8)+u866+u898+u802)))
  value_2_=value_2_+(c(50)*(u730*(u1308*(u940+u806*u1308)+u747)))
  u802=6.94346846494224497d-2*u856
  u655=-8.55178395377961841d4*u546
  u643=9.82004825756117503d-1*u764
  u702=u914*u1308+u637
  u940=u536*u1307
  u83=u721*u1307
  u872=(u310*(u1308*(u655+u601+u940+u1308*(u702+u643))+u83+u833+u802))
  value_3_=value_3_+(c(50)*u872)
  u843=4.0495130438605924d-2*u856
  u797=-4.36406397443353068d4*u546
  u711=4.45446434824665164d-1*u764
  u898=u568*u1306
  u786=u644*u1306
  u704=(u1308*(u797+u898+u568*u1307+u1308*(u645+u711))+u644*u1307+u786+&
 &u843)
  value_4_=value_4_+(c(50)*(u756*u704))
  u645=-5.5106892578906706d-4*u949*(-u862*(u547+8.8d1)+u844)
  u844=u594*(u789+6.01216346227151036d3*exppd2)
  u594=+u712*(u789+u557)
  u712=9.0926372755196065d-2*u764
  u608=3.47173423247112248d-2*u874
  u868=+u628*(u789+u578)
  u874=u553*(7.5d1+u575)
  u760=1.8185274551039213d-2*u884*(3.18404342432142532d7*u917-2.2687409&
 &2915906051d2*u874)
  value_5_=value_5_+(c(50)*(u1308*(u844+u868*u1306+u1308*((u712+u762)*u&
 &1308+u760*u1306+u594))+u608*u1306+u645))
  u762=u610*u1306
  u868=u762+u768*u1307
  u760=u931*u1306
  u608=u70*u1307+u760
  value_6_=value_6_+(c(50)*(u756*(u1308*(u961+u868+(u525+u565)*u1308)+u&
 &608+u775)))
  u960=1.92085269497412087d-2*u884*(u912-u862)
  u707=-7.04312654823844319d-2*u800
  u747=-2.11293796447153296d-1*u884*(u609+u578)
  u584=7.04312654823844319d-2*u884*(9.76439983458570431d6*u917+u557*(2.&
 &3d1+u575))
  u609=-9.20025467593118894d3*u546
  value_7_=value_7_+(c(50)*(u1308*(u960+u1306*(u693*u1306+u747)+u1308*(&
 &u98*u1308+u1306*(u584+u543)+u707))+u1306*(u839+u609*u1306)+u692))
  u960=1.42285384812897842d-2*u670
  u609=-6.08270020075138276d-1*u855
  u747=-u939
  u584=-5.869272123532036d-1*u764
  u707=9.60426347487060436d-2*u856
  u939=-u722
  u722=-u922
  value_1_=value_1_+(c(51)*(u1307*(u609+u939*u1306+u1307*(u1307*(u584+u&
 &896+u923)+u722*u1306+u747))+u707*u1306+u960))
  u609=u595*u1307
  value_2_=value_2_+(c(51)*(x*(u1307*(u533+u1307*(u609+u634))+u775)*z))
  value_3_=value_3_+(c(51)*((u790+u1307*(u1307*(u682+u89)+u99))*u1308+u&
 &1307*(u681+u934*u1306+u1307*(u1307*(u563+u723)+u554+u880))+u680*u1306&
 &+u941))
  value_4_=value_4_+(c(51)*(u730*((u673+u1307*(u743+u679))*u1308+u1307*&
 &(u78+u898+u1307*(u627+u562))+u786+u678)))
  value_5_=value_5_+(c(51)*u966)
  u747=-4.70510897373653379d-1*u855
  u786=-u919
  u714=-u920
  value_6_=value_6_+(c(51)*(y*(u1307*(u786+u762+u1307*(u525+u714))+u760&
 &+u747)*z))
  u714=-5.76255808492236261d-1*u855
  u786=-u79
  u681=-u586
  u747=u6*u1307
  u922=u747+u970
  u627=u985*u1306
  value_7_=value_7_+(c(51)*(u310*(u1307*(u786+u850+u1307*(u922+u681))+u&
 &627+u714)))
  u786=-1.28056846331608058d-1*u855
  value_1_=value_1_+(c(52)*(y*(u1307*(u878+u903+u1307*(u788+u550))+u623&
 &+u786)*z))
  value_2_=value_2_+(c(52)*(u310*((u976+u1307*(u609+u920))*u1308+u1307*&
 &(u961+u565*u1307)+u775)))
  value_3_=value_3_+(c(52)*(u730*(u589+u1307*(u993+u915+u1307*(u637+u85&
 &8))+u725+u790)))
  u786=u1306*(u635*u1306+u908)
  u966=u1306*(u861+u7*u1306)
  value_4_=value_4_+(c(52)*(u1308*(u685+u1307*(u1307*(u952+u572)+u925)+&
 &(u1307*(u580+u743)+u713)*u1308)+u1307*(u823+u786+u1307*(u1306*(u904+u&
 &590)+u576))+u966+u718))
  value_5_=value_5_+(c(52)*u684)
  u637=1.74263295323575325d-3*u949*u717
  u713=-7.84184828956088966d-2*u855
  u877=-4.3914350421540982d-1*u884*(-u860+u715)
  u681=-u720
  u586=5.75068874567798574d-2*u825
  u934=-8.62603311851697862d-1*u764
  u714=1.72520662370339572d-1*u771
  u78=-5.75068874567798574d-2*u671
  u975=u1306*(u604+u714)
  u939=u1306*(u78+u63)
  u722=u1306*(u978+u564*u1306)
  value_6_=value_6_+(c(52)*((u713+u1307*(u1307*(u934+u81)+u681))*u1308+&
 &u1307*(u877+u975+u1307*((u990+u63)*u1307+u939+u586))+u722+u637))
  u586=-u845
  u637=-9.86037716753382047d-1*u764
  value_7_=value_7_+(c(52)*(x*(u1307*(u586+u727+u1307*(u922+u637))+u752&
 &+u708)*z))
  u586=7.11426924064489211d-4*u949*u921
  u637=-3.20142115829020145d-2*u855
  u877=u949*u716
  u878=-3.84170538994824174d-2*u877
  u717=7.04312654823844319d-2*u831
  u921=-3.5215632741192216d-1*u764
  u716=-9.39083539765125759d-2*u764
  u715=2.11293796447153296d-1*u771
  u79=-7.04312654823844319d-2*u671
  u682=u1306*(u550*u1306+u715)
  u960=u1306*(u79+u945)
  u584=u1306*(u962+u573*u1306)
  value_1_=value_1_+(c(53)*((u637+u1307*(u1307*(u921+u923)+u528))*u1308&
 &+u1307*(u878+u682+u1307*((u716+u945)*u1307+u960+u717))+u584+u586))
  value_2_=value_2_+(c(53)*(u756*((u953+u1307*(u609+u566))*u1308+u1307*&
 &(u979+u560*u1307)+u780)))
  u586=u1306*(u641*u1306+u549)
  u637=u1306*(u721+u819*u1306)
  value_3_=value_3_+(c(53)*(u1308*(u828+u1307*(u1307*(u696+u723)+u551)+&
 &(u1307*(u581+u89)+u977)*u1308)+u1307*(u829+u586+u1307*(u1306*(u524+u6&
 &87)+u664))+u637+u728))
  value_4_=value_4_+(c(53)*(u730*(u1308*(u631+u937*u1307+(u59*u1307+u69&
 &5)*u1308)+u798*u1307+u827)))
  value_5_=value_5_+(c(53)*u813)
  u878=u949*u766
  u717=-4.70510897373653379d-2*u878
  u921=6.90082649481358289d-1*u830
  u716=1.72520662370339572d-1*u781
  u720=-1.15013774913559715d-1*u754
  value_6_=value_6_+(c(53)*(u730*(u1308*(u921+u1307*(u838+u720)+(u81+u9&
 &16)*u1308)+u1307*(u716+u889)+u717)))
  u717=-1.92085269497412087d-2*u878
  u921=7.04312654823844319d-2*u781
  u716=7.04312654823844319d-2*u770
  u720=-9.39083539765125759d-2*u527
  u937=u1306*(u928+u720)
  u581=u1306*(u716+u751)
  value_7_=value_7_+(c(53)*(u310*((u528+u1307*(u747+u676))*u1308+u1307*&
 &(u921+u937+u73*u1307)+u581+u717)))
  u717=-1.92085269497412087d-2*u875
  u921=4.22587592894306592d-1*u799
  u766=2.11293796447153296d-1*u792
  u718=-4.6954176988256288d-2*u531
  u845=9.39083539765125759d-2*u774
  value_1_=value_1_+(c(54)*(u730*(u1308*(u921+u1307*(u845*u1307+u718)+(&
 &u747+u591)*u1308)+u1307*(u766+u867*u1307)+u717)))
  value_2_=value_2_+(c(54)*(u310*(u1308*(u979+u566*u1307+(u609+u560)*u1&
 &308)+u953*u1307+u780)))
  value_3_=value_3_+(c(54)*(u730*(u1308*(u816+u620*u1307+(u60*u1307+u70&
 &6)*u1308)+u769*u1307+u841)))
  value_4_=value_4_+(c(54)*(u1308*(u840+u786+u1307*(u726+u703)+u1308*((&
 &u705+u619+u743)*u1308+u1307*(u580+u572)+u796+u583))+u1307*(u710+u555)&
 &+u966+u587))
  value_5_=value_5_+(c(54)*u661)
  u717=1.91689624855932858d-2*u949*u804
  u921=-7.84184828956088966d-2*u884*(-u862+u997)
  u766=5.75068874567798574d-2*u794
  u718=-9.41021794747306759d-2*u855
  u845=-u948
  value_6_=value_6_+(c(54)*(u1308*(u921+u975+u1307*(u889+u845)+u1308*((&
 &u74+u63)*u1308+u1307*(u899+u81)+u939+u766))+u1307*(u718+u776)+u722+u7&
 &17))
  u845=-1.92085269497412087d-2*u877
  u766=1.40862530964768864d-1*u782
  u717=-u980
  value_7_=value_7_+(c(54)*(u756*(u1308*(u766+u937+u1307*(u747+u550)+u9&
 &27*u1308)+u1307*(u717+u591*u1307)+u581+u845)))
  u845=-1.92085269497412087d-2*u884*(-u862+u912)
  u766=7.04312654823844319d-2*u800
  u717=-3.84170538994824174d-2*u855
  u718=-u662
  value_1_=value_1_+(c(55)*(u1308*(u845+u682+u1307*(u724+u718)+u1308*(u&
 &636*u1308+u1307*(u783+u923)+u960+u766))+u1307*(u717+u866)+u584+u832))
  value_2_=value_2_+(c(55)*(u756*(u1308*(u961+u920*u1307+(u609+u565)*u1&
 &308)+u976*u1307+u775)))
  value_3_=value_3_+(c(55)*(u1308*(u844+u586+u1307*(u940+u655)+u1308*((&
 &u712+u56+u89)*u1308+u1307*(u643+u723)+u741+u594))+u1307*(u802+u83)+u6&
 &37+u645))
  value_4_=value_4_+(c(55)*(u730*u704))
  value_5_=value_5_+(c(55)*u872)
  u845=-1.56836965791217793d-1*u855
  u766=-u846
  value_6_=value_6_+(c(55)*(u730*(u1308*(u766+u868+(u525+u990)*u1308)+u&
 &608+u845)))
  value_7_=value_7_+(c(55)*(u310*(u1308*(u718+u727+u550*u1307+(u922+u78&
 &3)*u1308)+u573*u1307+u752+u717)))
  value_1_=value_1_+(c(56)*(u730*(u1308*(u903+u847*u1307+(u788)*u1308)+&
 &u852*u1307+u623)))
  value_2_=value_2_+(c(56)*(u310*(u1308*(u533+u1308*(u595*u1308+u634))+&
 &u775)))
  u845=5.95154439852192426d-1*u856
  u766=-2.25672076558073263d5*u546
  u718=1.45482196408313704d0*u764
  u819=(u1308*(u766+u915+u858*u1307+u1308*(u702+u718))+u765*u1307+u725+&
 &u845)
  value_3_=value_3_+(c(56)*(u730*u819))
  u921=-1.79978357504915218d-2*u999
  u567=7.69407478333512556d-1*u856
  u644=-1.60015679062562792d5*u546
  u960=7.42410724707775274d-1*u764
  u719=-6.0742695657908886d-2*u855
  u962=6.54609596165029602d4*u546
  u662=-6.68169652236997746d-1*u764
  value_4_=value_4_+(c(56)*(u1308*(u567+u962*u1306+u962*u1307+u1308*(u1&
 &308*(u960+u891+u572+u82)+u662*u1307+u662*u1306+u644))+u719*u1307+u719&
 &*u1306+u921))
  value_5_=value_5_+(c(56)*(u756*u819))
  value_6_=value_6_+(c(56)*(u1308*(u537*u1306+u681*u1307+u1308*((u935+u&
 &81)*u1308+u934*u1307+u686*u1306))+u713*u1307+u824*u1306))
  value_7_=value_7_+(c(56)*(u756*(u1308*(u850+u676*u1307+(u922)*u1308)+&
 &u528*u1307+u627)))
  if ( lmax .eq. 5 ) go to 100
  u686=A13_v(pd,exppd2,erfpd)
  u567=pd2*u686
  u960=4.24539123242856709d5*u567
  u819=p**10*pd
  u662=u539*u819
  u962=u662*(-u862*(u575+1.3d1)+u960)
  u80=-4.4327369884018174d-2*u962
  u644=4.24539123242856709d5*u686
  u824=u662*pd2
  u766=u824*(-u574+u644)
  u637=-2.70889482624555508d-2*u766
  u765=u824*(u644-u574)
  u644=2.03167111968416631d0*u765
  u581=u662*u567
  u934=4.48512415451645461d5*u581
  u661=-5.83066140087139099d6*u581
  u921=+u661
  u635=6.36808684864285064d6*u686
  u546=u824*(u635+u578)
  u681=-3.5215632741192216d-1*u546
  u719=1.90164416802437966d0*u546
  u798=u824*(1.08257476426928461d8*u686+u578*(1.7d1+u575))
  u609=2.3477088494128144d-2*u798
  u797=-7.04312654823844319d-2*u798
  u794=u797*u1306
  u769=u609*u1306
  value_1_=value_1_+(c(57)*(y*((u637+u1306*(u1306*(u681+u769)+u934))*u1&
 &307+u1306*(u644+u1306*(u1306*(u719+u794)+u921))+u80)))
  u644=9.28957412763366928d-1*u765
  u921=-5.12692395207866046d6*u581
  u997=2.41528927318475401d0*u546
  u6=-1.15013774913559715d-1*u798
  u549=u6*u1306
  value_2_=value_2_+(c(57)*(x*(u1306*(u921+u1306*(u549+u997))+u644)*u73&
 &0))
  u868=u662*(u960-u862*(1.3d1+u575))
  u687=1.1445277689465239d-2*u868
  u585=8.39320363894117524d-2*u765
  u608=-2.09830090973529381d-2*u766
  u884=-5.24575227433823452d-1*u766
  u894=-1.38966489248918799d6*u581
  u949=3.47416223122296998d5*u581
  u844=1.50547030019662032d6*u581
  u694=1.09111647306235278d0*u546
  u564=-2.72779118265588195d-1*u546
  u682=-4.91002412878058751d-1*u546
  u770=-7.2741098204156852d-2*u798
  u872=1.8185274551039213d-2*u798
  u676=u872*u1306
  u705=u770*u1306
  u571=(u585+u1306*(u1306*(u694+u705)+u894))*u1308+(u608+u1306*(u1306*(&
 &u564+u676)+u949))*u1307
  u528=(y*(u571+u1306*(u884+u1306*(u1306*(u682+u676)+u844))+u687))
  value_3_=value_3_+(c(57)*u528)
  u695=2.80350903036502551d-2*u868
  u720=3.42651103711280895d-2*u765
  u866=-5.13976655566921343d-2*u766
  u573=-1.28494163891730336d0*u766
  u937=-5.67328316676358988d5*u581
  u980=8.50992475014538483d5*u581
  u594=3.68763405839633343d6*u581
  u636=4.45446434824665164d-1*u546
  u975=-6.68169652236997746d-1*u546
  u584=-1.20270537402659594d0*u546
  u7=-2.96964289883110109d-2*u798
  u838=4.45446434824665164d-2*u798
  u618=u838*u1306
  u867=u7*u1306
  value_4_=value_4_+(c(57)*(z*((u720+u1306*(u1306*(u636+u867)+u937))*u1&
 &308+(u866+u1306*(u1306*(u975+u618)+u980))*u1307+u1306*(u573+u1306*(u1&
 &306*(u584+u618)+u594))+u695)))
  u920=8.01169438262566727d-2*u868
  u721=5.87524254725882266d-1*u765
  u813=-1.46881063681470567d-1*u766
  u786=-1.3219295731332351d0*u766
  u938=-3.24255141580810531d6*u581
  u764=8.10637853952026328d5*u581
  u877=2.43191356185607899d6*u581
  u966=1.52756306228729389d0*u546
  u768=-3.81890765571823473d-1*u546
  u707=-6.00114060184294029d-1*u546
  u976=u1306*(u768+u676)
  value_5_=value_5_+(c(57)*(x*((u721+u1306*(u1306*(u966+u705)+u938))*u1&
 &308+(u813+u1306*(u976+u764))*u1307+u1306*(u786+u1306*(u1306*(u707+u67&
 &6)+u877))+u920)))
  u931=-3.61931459518194907d-2*u962
  u789=-6.63541009116690663d-2*u766
  u589=1.65885252279172666d0*u765
  u855=1.09862656115971296d6*u581
  u856=-4.76071509835875614d6*u581
  u613=+u856
  u703=-8.62603311851697862d-1*u546
  u704=1.55268596133305615d0*u546
  u702=5.75068874567798574d-2*u798
  u655=-5.75068874567798574d-2*u798
  u888=u655*u1306
  u617=u702*u1306
  value_6_=value_6_+(c(57)*(((u789+u1306*(u1306*(u703+u617)+u855))*u130&
 &7+u1306*(u589+u1306*(u1306*(u704+u888)+u613))+u931)*z))
  u931=-1.03430529729375739d-1*u962
  u589=-5.68867913511566566d-1*u766
  u613=1.7066037405346997d0*u765
  u60=3.13958690816151823d6*u581
  u939=-u60
  u566=-1.47905657513007307d0*u546
  u722=+u566
  u955=7.74743920306228751d-1*u546
  u814=7.04312654823844319d-2*u798
  u880=-2.3477088494128144d-2*u798
  u978=u880*u1306
  u73=u814*u1306
  value_7_=value_7_+(c(57)*(x*((u589+u1306*(u1306*(u722+u73)+u60))*u130&
 &7+u1306*(u613+u1306*(u1306*(u955+u978)+u939))+u931)))
  u931=4.06334223936833262d-1*u765
  u613=1.4950413848388182d5*u581
  u722=5.14690016986655464d-1*u765
  u955=-2.54157035422599095d6*u581
  u586=+u955
  u599=-2.3477088494128144d-1*u546
  u634=-5.98016553935527281d5*u581
  u922=+u634
  u775=1.33819404416530421d0*u546
  u605=7.04312654823844319d-2*u546
  u878=u605*u1306
  value_1_=value_1_+(c(58)*(x*(u1307*(u931+u1306*(u1306*(u775+u794)+u58&
 &6)+(u1306*(u599+u769)+u613)*u1307)+u1306*(u722+u1306*(u878+u922))+u80&
 &)))
  u931=-1.20643819839398302d-2*u962
  u586=1.32708201823338133d-1*u765
  u551=3.98124605470014398d-1*u765
  u895=-2.19725312231942591d6*u581
  u672=-7.32417707439808637d5*u581
  u879=+u672
  u698=1.72520662370339572d0*u546
  u991=1.15013774913559715d-1*u546
  u684=(u586+u1306*(u1306*(u698+u549)+u895))
  u595=u991*u1306
  u604=u1306*(u551+u1306*(u595+u879))+u931
  value_2_=value_2_+(c(58)*((u684*u1307+u604)*z))
  u747=1.06134780810714177d7*u567
  u993=u662*(u747-u862*(3.9d1+5.0d1*pd2))
  u616=3.81509256315507965d-3*u993
  u627=2.12269561621428355d6*u686
  u970=-2.72779118265588195d-1*u824*(-u574+u627)
  u80=5.79027038537161663d5*u581
  u631=5.73127816377856558d7*u686
  u886=-6.99433636578431269d-3*u824*(7.0330968803930876d3*exppd2+u631)
  u98=-3.62998548665449682d3*u553
  u911=9.0926372755196065d-2*u824*(7.0048955335071357d7*u686+u98)
  u559=-9.0926372755196065d-1*u546
  u558=1.45482196408313704d-1*u824*(u635+u557)
  u713=-1.8185274551039213d-2*u824
  u985=+u713*(8.59691724566784836d8*u686+u578*(1.35d2+u805))
  u563=9.0926372755196065d-2*u798
  u692=-9.0926372755196065d-2*u546
  u524=u563*u1306
  u925=u1306*(u985+u524)
  u649=(x*(u1307*(u970+u1306*(u925+u911)+(u1306*(u559+u524)+u80)*u1307)&
 &+u1306*(u886+u1306*(u692*u1306+u558))+u616))
  value_3_=value_3_+(c(58)*u649)
  u796=-4.11181324453537075d-1*u766
  u896=-1.89109438892119663d5*u581
  u799=2.83664158338179494d5*u581
  u792=2.17475854725937612d6*u581
  u779=2.96964289883110109d-1*u546
  u871=-4.45446434824665164d-1*u546
  u837=-9.79982156614263361d-1*u546
  value_4_=value_4_+(c(58)*(u545*((u896+u1306*(u867+u779))*u1308+(u799+&
 &u1306*(u618+u871))*u1307+u1306*(u792+u1306*(u618+u837))+u796)))
  value_5_=value_5_+(c(58)*u528)
  u528=2.65416403646676265d-1*u765
  u875=3.66208853719904319d5*u581
  u601=-1.83104426859952159d6*u581
  u935=+u601
  u996=-5.75068874567798575d-1*u546
  u767=1.03512397422203744d0*u546
  value_6_=value_6_+(c(58)*(u310*((u875+u1306*(u617+u996))*u1307+u1306*&
 &(u935+u1306*(u888+u767))+u528)*z))
  u528=4.92526332044646377d-3*u868
  u935=-8.12668447873666523d-2*u766
  u587=2.70889482624555508d-2*u765
  u611=1.34553724635493638d6*u581
  u945=-7.47520692419409102d5*u581
  u912=+u945
  u63=-1.05646898223576648d0*u546
  u693=+u63
  u619=4.46064681388434736d-1*u546
  u56=(u935+u1306*(u1306*(u693+u73)+u611))*u1307
  value_7_=value_7_+(c(58)*(y*(u56+u1306*(u587+u1306*(u1306*(u619+u978)&
 &+u912))+u528)))
  u619=4.87601068724199913d-1*u765
  u927=-2.84057863119375459d6*u581
  u897=+u927
  u833=1.40862530964768864d0*u546
  value_1_=value_1_+(c(59)*(u310*((u613+u1306*(u769+u599))*u1307+u1306*&
 &(u897+u1306*(u794+u833))+u619)*z))
  value_2_=value_2_+(c(59)*(y*(u684*u1308+u604)))
  u619=4.19660181947058762d-2*u765
  u897=-4.6322163082972933d5*u581
  u684=1.15805407707432333d5*u581
  u604=7.2741098204156852d-1*u546
  u791=-1.8185274551039213d-1*u546
  u952=-2.18223294612470556d-1*u546
  u845=(u545*((u897+u1306*(u705+u604))*u1308+(u684+u1306*(u676+u791))*u&
 &1307+u1306*(u684+u1306*(u676+u952))+u619))
  value_3_=value_3_+(c(59)*u845)
  u89=2.12269561621428355d6*u567
  u81=-4.67251505060837585d-3*u662*(-u574*(u548-3.9d1)+u89)
  u861=-1.54192996670076403d-1*u766
  u923=-u799
  u723=3.42651103711280895d-2*u824*(u627-u820)
  u820=1.22921135279877781d6*u581
  u572=u824*(1.40097910670142714d7*u686+u972)
  u639=-2.22723217412332582d-1*u572
  u645=u824*(u627+u578)
  u658=-4.45446434824665164d-2*u645
  u968=-7.57258939201930779d-1*u546
  u548=u824*(1.71938344913356967d8*u686-4.53748185831812103d2*u61)
  u61=4.45446434824665164d-2*u548
  u788=-4.45446434824665164d-2*u798
  u795=u788*u1306
  u714=u61+u795
  u715=u1306*(u714)
  value_4_=value_4_+(c(59)*(x*(u1308*(u861+u1306*(u1306*(u968+u618)+u82&
 &0)+(u1306*(u779+u867)+u896)*u1308)+u1307*(u799+u1306*(u715+u639)+(u13&
 &06*(u636+u795)+u923)*u1307)+u1306*(u723+u658*u1306)+u81)))
  u82=-7.6301851263101593d-3*u962
  u724=1.0491504548676469d-1*u765
  u659=-3.09149667367666621d-1*u546
  value_5_=value_5_+(c(59)*(z*(u571+u1306*(u724+u1306*(u1306*(u659+u676&
 &)+u949))+u82)))
  u777=7.21716509512856406d6*u567
  u571=1.7d1*pd2
  u821=-6.03219099196991512d-3*u662*(-u574*(u571+3.9d1)+u777)
  u725=4.64478706381683464d-1*u765
  u889=-u875
  u641=1.1578339724805183d5*u686
  u717=3.89277392015125189d0*u824
  u780=u717*(u641+u2160)
  u782=-2.56346197603933023d6*u581
  u940=+u782
  u726=-2.87534437283899287d-1*u572
  u891=5.75068874567798575d-1*u546
  u606=1.91042605459285519d7*u686
  u620=u824*(u606+u578)
  u831=-5.75068874567798574d-2*u620
  u793=1.20764463659237701d0*u546
  u83=5.75068874567798574d-2*u548
  value_6_=value_6_+(c(59)*(x*((u725+u1306*(u1306*(u793+u888)+u940))*u1&
 &308+u1307*(u875+u1306*(u1306*(u83+u888)+u726)+(u1306*(u891+u888)+u889&
 &)*u1307)+u1306*(u780+u1306*(u595+u831))+u821)))
  u83=-9.85052664089292755d-3*u962
  u780=-1.64454552332270002d6*u581
  u821=+u780
  u726=5.869272123532036d-1*u546
  value_7_=value_7_+(c(59)*((u56+u1306*(u722+u1306*(u1306*(u726+u978)+u&
 &821))+u83)*z))
  u83=-1.47757899613393913d-2*u962
  u821=2.99008276967763641d4*u581
  u726=5.68867913511566566d-1*u765
  u831=-1.40862530964768864d-1*u546
  u56=-1.19603310787105456d6*u581
  u595=+u56
  u898=8.2169809729448504d-1*u546
  u600=2.11293796447153296d-1*u546
  u632=(u831+u769)
  u79=u1306*u632
  u640=(u79+u821)*u1307
  u74=u600*u1306
  value_1_=value_1_+(c(60)*(y*(u1307*(u587+u1306*(u1306*(u898+u794)+u91&
 &2)+u640)+u1306*(u726+u1306*(u74+u595))+u83)))
  u595=-1.46483541487961727d6*u581
  u898=+u595
  u610=1.15013774913559715d0*u546
  u873=3.45041324740679145d-1*u546
  u623=(u879+u1306*(u549+u610))
  u706=u873*u1306
  u915=u1306*(u898+u706)+u551
  value_2_=value_2_+(c(60)*(u310*(u623*u1307+u915)*z))
  u752=2.28905553789304779d-3*u993
  u993=3.18404342432142532d7*u686
  u984=-4.19660181947058762d-3*u824*(5.21810413706583918d3*exppd2+u993)
  u999=-2.6438591462664702d-1*u824*(-u574+2.52701859073128994d6*u686)
  u858=5.4555823653117639d-2*u824*(7.42943465674999241d7*u686+u98)
  u908=-5.4555823653117639d-1*u546
  u852=1.06134780810714177d7*u686
  u554=u824*(u852+u578)
  u674=2.18223294612470556d-1*u554
  u751=+u713*(7.9601085608035633d8*u686+u578*(1.25d2+u805))
  u805=u1306*(u751+u524)
  u841=u564*u1306
  u762=(y*(u1307*(u984+u1306*(u805+u858)+(u1306*(u908+u524)+u684)*u1307&
 &)+u1306*(u999+u1306*(u841+u674))+u752))
  value_3_=value_3_+(c(60)*u762)
  u781=2.97177386269999696d6*u567
  u543=3.73801204048670068d-3*u662*(u781-u860*(3.9d1+2.8d1*pd2))
  u776=u824*(u660+u635)
  u785=-6.85302207422561791d-3*u776
  u760=3.78218877784239326d4*u581
  u990=-1.13074864224722696d-1*u766
  u827=5.67328316676358989d4*u581
  u630=-1.02795331113384269d-2*u824*(4.31060776540221498d3*exppd2+3.014&
 &22777502428264d7*u686)
  u722=1.48588693134999848d7*u686
  u555=u824*(u722+u972)
  u883=8.90892869649330328d-2*u555
  u750=-1.78178573929866066d-1*u546
  u526=1.53178645502616927d6*u581
  u77=-2.67267860894799098d-1*u546
  u773=8.90892869649330328d-2*u824*(1.23116345740428446d7*u686+u578)
  u893=1.59202171216071266d8*u686
  u718=u824*(u893-4.53748185831812103d2*u951)
  u951=-2.96964289883110109d-2*u718
  u928=2.96964289883110109d-2*u798
  u839=-9.35437513131796845d-1*u546
  u582=-1.33633930447399549d-1*u546
  u95=u928*u1306
  u806=u1306*(u951+u95)
  u628=u1306*(u839+u618)
  value_4_=value_4_+(c(60)*(z*(u1308*(u785+u1306*(u806+u883)+(u1306*(u7&
 &50+u95)+u760)*u1308)+u1307*(u990+u1306*(u628+u526)+(u1306*(u77+u618)+&
 &u827)*u1307)+u1306*(u630+u1306*(u582*u1306+u773))+u543)))
  value_5_=value_5_+(c(60)*u649)
  u649=-2.41287639678796605d-3*u962
  u607=-1.32708201823338133d-2*u766
  u948=7.32417707439808637d4*u581
  u678=1.19437381641004319d-1*u765
  u525=-2.19725312231942591d5*u581
  u531=+u525
  u612=-3.45041324740679145d-1*u546
  u550=-2.92967082975923455d5*u581
  u941=+u550
  u998=5.17561987111018717d-1*u546
  u670=5.75068874567798574d-2*u546
  u919=u1306*(u612+u617)
  value_6_=value_6_+(c(60)*((u1307*(u607+u1306*(u1306*(u998+u888)+u531)&
 &+(u919+u948)*u1307)+u1306*(u678+u1306*(u670*u1306+u941))+u649)*z))
  u649=-2.97978430887011058d-1*u766
  u607=-u945
  u678=-7.0431265482384432d-1*u546
  u531=-1.19603310787105456d5*u581
  u941=+u531
  u998=2.3477088494128144d-2*u546
  u945=(u605+u978)
  u741=u998*u1306
  value_7_=value_7_+(c(60)*(x*(u1307*(u649+u1306*(u1306*u945+u607)+(u13&
 &06*(u678+u73)+u934)*u1307)+u1306*(u587+u1306*(u741+u941))+u528)))
  u649=-5.91031598453575652d-3*u962
  u528=4.87601068724199913d-2*u765
  u59=2.11293796447153296d-1*u765
  u643=-9.86727313993620014d5*u581
  u671=+u643
  u802=-4.18611587754869097d5*u581
  u840=+u802
  u846=9.15606451270997615d-1*u546
  value_1_=value_1_+(c(61)*((u1307*(u528+u1306*(u1306*(u846+u794)+u671)&
 &+u640)+u1306*(u59+u1306*(u878+u840))+u649)*z))
  u649=-u860*(u556-1.3d1)
  u528=u662*(u89+u649)
  u59=2.41287639678796605d-2*u528
  u556=-u672
  u840=1.14625563275571312d7*u686
  u846=u824*(-u574+u840)
  u671=-4.42360672744460442d-2*u846
  u672=5.75068874567798575d-1*u572
  u640=-u610
  u70=9.20110199308477719d-1*u824*(1.27361736972857013d6*u686+u859)
  u832=-1.15013774913559715d-1*u548
  u664=1.15013774913559715d-1*u798
  u576=-1.15013774913559715d-1*u546
  u816=u664*u1306
  u800=u1306*(u832+u816)
  value_2_=value_2_+(c(61)*(x*(u1307*(u879+u1306*(u800+u672)+(u1306*(u6&
 &40+u816)+u556)*u1307)+u1306*(u671+u1306*(u576*u1306+u70))+u59)))
  u712=4.88331848083850195d-2*u662*(u960-u857*(6.4d1*pd2-1.3d1))
  u857=-1.67864072778823505d-2*u776
  u995=9.26443261659458661d4*u581
  u562=-4.61626200141764638d-2*u766
  u583=2.31610815414864665d4*u581
  u68=-3.77694163752352886d-2*u824*(-u574+1.51418953956618893d7*u686)
  u602=2.18223294612470556d-1*u555
  u817=-4.36446589224941112d-1*u546
  u614=6.25349201620134597d5*u581
  u685=-1.09111647306235278d-1*u546
  u986=1.82551822994428385d7*u686
  u673=1.09111647306235278d-1*u824*(u986+u972)
  u8=-7.2741098204156852d-2*u718
  u668=7.2741098204156852d-2*u798
  u570=-2.36408569163509769d-1*u546
  u753=u668*u1306
  u743=u1306*(u8+u753)
  u850=(z*(u1308*(u857+u1306*(u743+u602)+(u1306*(u817+u753)+u995)*u1308&
 &)+u1307*(u562+u1306*(u976+u614)+(u1306*(u685+u676)+u583)*u1307)+u1306&
 &*(u68+u1306*(u570*u1306+u673))+u712))
  value_3_=value_3_+(c(61)*u850)
  u84=-2.80350903036502551d-3*u662*(-u574*(9.0d0*pd2-1.3d1)+3.820852109&
 &18571038d6*u567)
  u811=-1.02795331113384269d-2*u766
  u942=-u760
  u774=u824*(u635+u660)
  u621=1.02795331113384269d-2*u774
  u913=-u827
  u727=4.60523083387961523d0*u824*(6.06484461775509585d4*u686-u69)
  u859=3.97129821673451292d5*u581
  u892=1.78178573929866066d-1*u546
  u680=-1.33633930447399549d-1*u555
  u992=2.67267860894799098d-1*u546
  u527=-1.33633930447399549d-1*u824
  u532=+u527*(6.65111293080475511d6*u686+u578)
  u533=-4.8999107830713168d-1*u546
  u669=4.45446434824665164d-2*u718
  u909=8.90892869649330328d-2*u546
  u657=u1306*(u669+u795)
  value_4_=value_4_+(c(61)*(y*(u1308*(u811+u1306*(u1306*(u533+u618)+u85&
 &9)+(u1306*(u892+u867)+u942)*u1308)+u1307*(u621+u1306*(u657+u680)+(u13&
 &06*(u992+u795)+u913)*u1307)+u1306*(u727+u1306*(u909*u1306+u532))+u84)&
 &))
  value_5_=value_5_+(c(61)*u845)
  u580=1.48588693134999848d7*u567
  u845=3.5d1*pd2
  u754=u662*(-u574*(u845+1.3d1)+u580)
  u716=-1.20643819839398302d-3*u754
  u622=6.63541009116690663d-2*u765
  u728=1.32708201823338133d-2*u774
  u903=-u948
  u899=7.96249210940028796d-2*u774
  u914=-u855
  u536=-1.72520662370339572d-1*u555
  u99=-1.72520662370339572d-1*u554
  u910=8.62603311851697862d-1*u546
  u62=5.75068874567798574d-2*u718
  u666=2.3002754982711943d-1*u546
  u537=u1306*(u873+u888)
  u863=(u537+u903)
  value_6_=value_6_+(c(61)*(y*((u622+u1306*(u1306*(u910+u888)+u914))*u1&
 &308+u1307*(u728+u1306*(u1306*(u62+u888)+u536)+u863*u1307)+u1306*(u899&
 &+u1306*(u666*u1306+u99))+u716)))
  u716=-5.41778965249111015d-2*u766
  u899=-u613
  u969=2.81725061929537728d-1*u546
  value_7_=value_7_+(c(61)*(u310*((u934+u1306*(u73+u678))*u1307+u1306*(&
 &u899+u1306*(u978+u969))+u716)*z))
  u716=6.36808684864285064d6*u567
  u711=1.5d1*pd2
  u65=u662*(-u574*(u711+2.6d1)+u716)
  u829=-1.47757899613393913d-3*u65
  u590=8.12668447873666523d-2*u765
  u729=5.41778965249111015d-3*u774
  u904=-u821
  u902=3.03242230887754792d5*u686
  u803=u824*(u902-u860)
  u843=9.10188661618506505d-1*u803
  u905=-u611
  u932=-7.04312654823844319d-2*u555
  u697=1.40862530964768864d-1*u546
  u591=u824*(u722+u578)
  u961=-7.04312654823844319d-2*u591
  u633=-u63
  u63=2.3477088494128144d-2*u718
  u778=u1306*(u697+u978)
  u979=u697*u1306
  value_1_=value_1_+(c(62)*(y*((u590+u1306*(u1306*(u633+u794)+u905))*u1&
 &308+u1307*(u729+u1306*(u1306*(u63+u978)+u932)+(u778+u904)*u1307)+u130&
 &6*(u843+u1306*(u979+u961))+u829)))
  value_2_=value_2_+(c(62)*(u545*(u623*u1308+u915)))
  u961=1.57079475599856982d7*u567
  u843=3.7d1*pd2
  u829=u662*(-u574*(u843+2.6d1)+u961)
  u994=-3.81509256315507965d-4*u829
  u623=7.97354345699411647d-2*u765
  u915=-u995
  u568=4.19660181947058762d-3*u774
  u887=-u583
  u624=5.03592218336470514d-2*u824*(2.97177386269999696d6*u686-u862)
  u916=-7.64315690869053396d5*u581
  u72=4.36446589224941112d-1*u546
  u708=-5.4555823653117639d-2*u555
  u826=1.09111647306235278d-1*u546
  u756=-5.4555823653117639d-2*u824
  u842=+u756*(8.06624334161427748d6*u686+u578)
  u933=1.63667470959352917d-1*u546
  u907=1.8185274551039213d-2*u718
  u565=-1.8185274551039213d-2*u798
  u71=3.6370549102078426d-2*u546
  u683=u565*u1306
  u710=u1306*(u907+u683)
  u783=(y*(u1308*(u623+u1306*(u1306*(u933+u676)+u916)+(u1306*(u72+u705)&
 &+u915)*u1308)+u1307*(u568+u1306*(u710+u708)+(u1306*(u826+u683)+u887)*&
 &u1307)+u1306*(u624+u1306*(u71*u1306+u842))+u994))
  value_3_=value_3_+(c(62)*u783)
  u812=-9.34503010121675169d-4*u662*(-u574*(u547-5.2d1)+u747)
  u560=9.64861643733765248d5*u686
  u847=u824*(u560-u862)
  u730=7.5383242816481797d-2*u847
  u917=-9.45547194460598314d4*u581
  u804=u824*(1.51621115443877396d6*u686-u660)
  u731=1.43913463558737976d-1*u804
  u918=3.53782602702380591d6*u686
  u732=u824*(u918+u557)
  u787=-5.34535721789598197d-1*u732
  u876=+u527*(u918+u578)
  u85=2.96964289883110109d-2*u824*(u893-4.53748185831812103d2*u691)
  u691=-7.42410724707775273d-2*u798
  u835=u691*u1306
  u744=u876*u1306
  value_4_=value_4_+(c(62)*(z*(u1308*(u730+u1306*(u85*u1306+u787)+(u130&
 &6*(u636+u835)+u917)*u1308)+u1306*(u731+u744)+u812)))
  u86=-1.90754628157753983d-3*u65
  u733=3.56711154654999948d-1*u765
  u906=-u684
  u734=1.39886727315686254d-2*u774
  u943=-1.04224866936689099d6*u581
  u57=-9.0926372755196065d-2*u572
  u930=1.8185274551039213d-1*u546
  u679=-1.8185274551039213d-2*u546
  u885=5.4555823653117639d-2*u546
  u64=1.8185274551039213d-2*u548
  u78=u64+u683
  u953=u1306*(u78)
  value_5_=value_5_+(c(62)*(x*(u1308*(u733+u1306*(u1306*(u885+u676)+u94&
 &3)+(u1306*(u604+u705)+u897)*u1308)+u1307*(u684+u1306*(u953+u57)+(u130&
 &6*(u930+u683)+u906)*u1307)+u1306*(u734+u679*u1306)+u86)))
  u87=-3.61931459518194907d-3*u65
  u65=3.71582965105346771d-1*u803
  u977=2.22949779063208063d0*u803
  u882=-3.45041324740679145d-1*u554
  u803=-1.72520662370339572d-1*u591
  u836=u824*(u993+u950*(2.0d1+pd2))
  u588=4.6005509965423886d-1*u836
  value_6_=value_6_+(c(62)*(z*(u1308*(u65+u1306*(u1306*(u588+u549)+u882&
 &)+u863*u1308)+u1306*(u977+u1306*(u706+u803))+u87)))
  u65=-2.46263166022323189d-3*u829
  u977=1.89622637837188855d-1*u765
  u882=-u934
  u863=7.21716509512856406d6*u686
  u588=5.41778965249111015d-2*u824*(u863-u818)
  u803=-1.04652896938717274d6*u581
  u561=+u803
  u829=-3.5215632741192216d-1*u572
  u881=7.0431265482384432d-1*u546
  u735=u824*(1.3160712820528558d7*u686+u578)
  u900=-7.04312654823844319d-2*u735
  u603=4.93018858376691024d-1*u546
  u964=7.04312654823844319d-2*u548
  u529=9.39083539765125759d-2*u546
  value_7_=value_7_+(c(62)*(x*((u977+u1306*(u1306*(u603+u978)+u561))*u1&
 &308+u1307*(u934+u1306*(u1306*(u964+u794)+u829)+(u1306*(u881+u794)+u88&
 &2)*u1307)+u1306*(u588+u1306*(u529*u1306+u900))+u65)))
  u65=2.43800534362099957d-1*u765
  u977=-7.04312654823844319d-2*u546
  u588=1.35444741312277754d-1*u765
  u561=3.5215632741192216d-1*u546
  u829=-8.97024830903290922d4*u581
  u900=+u829
  u603=4.22587592894306592d-1*u546
  u964=(u769+u977)
  u529=u603*u1306
  u755=u900*u1306
  value_1_=value_1_+(c(63)*(x*(u1307*(u65+u1306*((u561+u794)*u1307+u529&
 &+u905)+u964*y**4)+u1306*(u588+u755)+u83)))
  u65=-7.23862919036389814d-3*u962
  u588=1.59249842188005759d-1*u765
  u772=-1.46483541487961727d5*u581
  u890=+u772
  u924=-1.75780249785554073d6*u581
  u870=6.90082649481358289d-1*u546
  u638=(u1306*(u870+u549)+u890)
  u946=u870*u1306
  u822=u588+u1306*(u946+u924)
  u815=u890*u1306
  u93=u1306*(u588+u815)+u65
  value_2_=value_2_+(c(63)*((u1307*(u822+u638*u1307)+u93)*z))
  u83=u684*u1306
  u759=(x*(u1307*(u999+u1306*(u908*u1306+u858)+u1307*((u564+u524)*u1307&
 &+u805+u674))+u1306*(u984+u83)+u752))
  value_3_=value_3_+(c(63)*u759)
  u805=u824*(-u818+u918)
  u823=-6.16771986680305612d-2*u805
  u736=u824*(u993+u972)
  u579=4.45446434824665164d-2*u736
  u981=-2.22723217412332582d-1*u546
  u92=u824*(u893-4.53748185831812103d2*u869)
  u893=-5.93928579766220219d-2*u92
  u88=7.42410724707775273d-2*u798
  u737=u88*u1306
  value_4_=value_4_+(c(63)*(u310*(u1307*(u579+u1306*(u737+u893)+(u737+u&
 &981)*u1307)+u1306*(u579+u981*u1306)+u823)*z))
  value_5_=value_5_+(c(63)*u762)
  u823=-1.72520662370339572d-1*u546
  u893=1.72520662370339572d-1*u546
  value_6_=value_6_+(c(63)*(u310*(u1307*(u875+u655*u665+(u617+u823)*u13&
 &07)+u1306*(u889+u893*u1306))*z))
  u762=1.47757899613393913d-2*u868
  u737=-1.35444741312277754d-1*u766
  u981=-u829
  u829=-2.43800534362099957d-1*u766
  u579=-4.22587592894306592d-1*u546
  u648=u1306*(u579+u73)
  u598=(u648+u981)*u1307
  value_7_=value_7_+(c(63)*(y*(u1307*(u737+u1306*(u1306*(u681+u978)+u61&
 &1)+u598)+u1306*(u829+u605*u665)+u762)))
  u737=1.62533689574733305d-1*u765
  u829=4.6954176988256288d-1*u546
  value_1_=value_1_+(c(64)*(u310*(u1307*(u899+u1306*(u794+u829)+u964*u1&
 &307)+u1306*(u912+u74)+u737)*z))
  u829=1.44772583807277963d-2*u528
  u74=-2.65416403646676265d-2*u776
  u964=-u772
  u772=u824*(u852+u574)
  u967=-7.96249210940028796d-2*u772
  u688=3.45041324740679145d-1*u555
  u790=-6.90082649481358289d-1*u546
  u552=u824*(u627+u950)
  u689=1.38016529896271658d0*u552
  u9=-1.15013774913559715d-1*u718
  u757=u1306*(u9+u816)
  u825=u74+u1306*(u757+u688)
  u535=u612*u1306
  u809=u1306*(u967+u1306*(u535+u689))+u829
  u53=(u1306*(u790+u816)+u964)
  value_2_=value_2_+(c(64)*(y*(u1307*(u825+u53*u1307)+u809)))
  u784=-2.51796109168235257d-2*u772
  u646=u824*(u993+u901)
  u596=5.4555823653117639d-2*u646
  u538=(u310*(u1307*(u596+u1306*(u524+u8)+(u524+u564)*u1307)+u1306*(u59&
 &6+u841)+u784)*z)
  value_3_=value_3_+(c(64)*u538)
  u771=u662*(u649+u89)
  u89=-1.1214036121460102d-2*u771
  u649=-2.05590662226768537d-2*u772
  u650=3.56357147859732131d-1*u552
  u522=-8.90892869649330328d-2*u546
  u848=u824*(4.95295643783332827d6*u686+u860)
  u738=1.23354397336061122d-1*u848
  u577=-1.33633930447399549d-1*u554
  u690=1.33633930447399549d-1*u546
  u530=4.54863346331632188d6*u686
  u739=4.79711545195793254d-2*u824*(u530-u862)
  u740=u824*(u631+u98)
  u98=-8.90892869649330328d-2*u740
  u950=u824*(2.22883039702499772d8*u686+u578*(3.5d1+u575))
  u66=4.45446434824665164d-2*u950
  u808=-4.45446434824665164d-2*u554
  u544=u824*(2.10146866005214071d8*u686+u578*(3.3d1+u575))
  u90=4.45446434824665164d-2*u544
  u854=-8.90892869649330328d-2*u798
  u651=4.45446434824665164d-2*u546
  u97=u854*u1306
  value_4_=value_4_+(c(64)*(x*(u1308*(u649+u1306*(u750*u1306+u883)+u130&
 &8*((u522+u95)*u1308+u806+u650))+u1307*(u738+u1306*(u1306*(u90+u795)+u&
 &98)+u1307*((u690+u795)*u1307+u1306*(u66+u97)+u577))+u1306*(u739+u1306&
 &*(u651*u1306+u808))+u89)))
  value_5_=value_5_+(c(64)*u850)
  u665=4.82575279357593209d-3*u528
  u806=u824*(7.07565205404761182d5*u686+u862)
  u850=7.96249210940028796d-2*u806
  u849=u824*(u606+u862)
  u58=-8.84721345488920883d-3*u849
  u94=9.55213027296427596d7*u686
  u542=u824*(u94-4.53748185831812103d2*u677)
  u597=5.75068874567798574d-2*u542
  u625=u824*(8.9153215880999909d6*u686+u578)
  u96=5.75068874567798574d-2*u625
  u853=u824*(1.33729823821499863d8*u686-4.53748185831812103d2*u675)
  u944=-5.75068874567798574d-2*u853
  u801=-5.75068874567798574d-2*u546
  u667=(u893+u888)
  value_6_=value_6_+(c(64)*(x*(u1307*(u850+u1306*(u1306*(u944+u617)+u66&
 &6)+u1307*(u667*u1307+u597*u1306+u823))+u1306*(u58+u1306*(u801*u1306+u&
 &96))+u665)))
  u665=1.97010532817858551d-3*u868
  u850=-7.04312654823844319d-2*u766
  u58=-1.62533689574733305d-2*u766
  u597=6.27917381632303645d5*u581
  u96=-5.98016553935527281d4*u581
  u944=+u96
  value_7_=value_7_+(c(64)*((u1307*(u850+u1306*(u1306*(u977+u978)+u597)&
 &+u598)+u1306*(u58+u1306*(u741+u944))+u665)*z))
  u665=1.18206319690715131d-2*u528
  u850=u824*(u2160+u627)
  u58=-1.30026951659786644d-1*u850
  u598=-7.04312654823844319d-2*u645
  u741=u824*(2.75950430107856861d7*u686+u901)
  u666=1.40862530964768864d-1*u741
  u901=-u578*(u575-5.0d0)
  u626=-2.3477088494128144d-2*u824*(u901+u993)
  u742=u824*(9.76439983458570431d6*u686+u578)
  u653=7.04312654823844319d-2*u742
  u540=-7.04312654823844319d-2*u718
  u91=4.6954176988256288d-2*u798
  u965=u1306*(u540+u73)
  u88=u977*u1306
  u523=u91*u1306
  value_1_=value_1_+(c(65)*(x*(u1307*(u58+u1306*(u965+u666)+u1307*(u945&
 &*u1307+u1306*(u626+u523)+u598))+u1306*(u58+u1306*(u88+u653))+u665)))
  value_2_=value_2_+(c(65)*(z*(u1308*(u825+u53*u1308)+u809)))
  u665=3.05207405052406372d-3*u528
  u598=-5.03592218336470514d-2*u772
  u666=8.72893178449882224d-1*u552
  u626=5.03592218336470514d-2*u848
  u825=-5.4555823653117639d-2*u554
  u945=u539*exppd2*u819
  u819=u945*pd2
  u848=-4.12576766819373212d0*u819
  u58=-3.6370549102078426d-2*u740
  u809=1.8185274551039213d-2*u950
  u53=1.8185274551039213d-2*u544
  u740=-3.6370549102078426d-2*u798
  u67=1.8185274551039213d-2*u546
  u663=u740*u1306
  u539=(x*(u1308*(u598+u1306*(u817*u1306+u602)+u1308*((u952+u753)*u1308&
 &+u743+u666))+u1307*(u626+u1306*(u1306*(u53+u683)+u58)+u1307*((u885+u6&
 &83)*u1307+u1306*(u809+u663)+u825))+u1306*(u848+u1306*(u67*u1306+u679)&
 &)+u665))
  value_3_=value_3_+(c(65)*u539)
  u743=7.1956731779368988d-2*u804
  u699=-1.48482144941555055d-1*u546
  u936=u824*(u993-2.26874092915906051d2*u926)
  u675=1.78178573929866066d-1*u936
  u615=u795+u675
  value_4_=value_4_+(c(65)*(u545*(u1308*(u917+u1306*(u618+u699)+(u867+u&
 &909)*u1308)+u1307*(u582+u1306*(u615)+(u795+u690)*u1307)+u744+u743)))
  value_5_=value_5_+(c(65)*u783)
  u929=3.98124605470014398d-2*u774
  u541=2.3002754982711943d-1*u936
  value_6_=value_6_+(c(65)*(u545*((u889+u1306*(u888+u891))*u1308+u1307*&
 &(u823+u1306*(u888+u541)+(u888+u893)*u1307)+u1306*(u99+u706)+u929)))
  u929=4.03312167080713874d7*u567
  u541=-u574*(9.5d1*pd2-2.6d1)
  u783=-4.92526332044646377d-4*u662*(u541+u929)
  u744=1.62533689574733305d-2*u774
  u851=u824*(u627+u2160)
  u627=2.60053903319573287d-1*u851
  u834=-2.11293796447153296d-1*u555
  u745=u824*(9.19834767026189537d6*u686+u578)
  u534=-2.11293796447153296d-1*u745
  u54=7.04312654823844319d-2*u718
  u677=2.3477088494128144d-1*u546
  u629=u1306*(u603+u794)
  u696=(u629+u900)
  u869=u677*u1306
  value_7_=value_7_+(c(65)*(y*((u587+u1306*(u1306*(u561+u978)+u882))*u1&
 &308+u1307*(u744+u1306*(u1306*(u54+u794)+u834)+u696*u1307)+u1306*(u627&
 &+u1306*(u869+u534))+u783)))
  u783=u824*(u918-u818)
  u627=4.87601068724199913d-2*u783
  u918=-7.04312654823844319d-2*u620
  u926=9.39083539765125759d-2*u936
  u76=u969*u1306
  value_1_=value_1_+(c(66)*(u545*((u882+u1306*(u794+u881))*u1308+u1307*&
 &(u977+u1306*(u978+u926)+(u978+u605)*u1307)+u1306*(u918+u76)+u627)))
  value_2_=value_2_+(c(66)*(y*(u1308*(u822+u638*u1308)+u93)))
  u627=3.77694163752352886d-2*u783
  u918=-u764
  u926=2.18223294612470556d-1*u546
  u638=-5.4555823653117639d-2*u546
  u822=5.4555823653117639d-1*u546
  u652=7.2741098204156852d-2*u936
  u654=u683+u652
  u93=u825*u1306
  u807=(u545*(u1308*(u918+u1306*(u676+u822)+(u705+u926)*u1308)+u1307*(u&
 &638+u1306*(u654)+(u683+u885)*u1307)+u93+u627))
  value_3_=value_3_+(c(66)*u807)
  u647=3.08385993340152806d-2*u772
  u810=2.22723217412332582d-1*u546
  u746=2.39855772597896627d-2*u804
  u700=+u527*(u852+u972)
  value_4_=value_4_+(c(66)*(x*(u1308*(u647+u700*u1306+u1308*((u810+u835&
 &)*u1308+u669*u1306+u577))+u746*u1306+u812)))
  u55=-1.1445277689465239d-3*u662*(-u574*(u547+9.1d1)+u747)
  u747=8.39320363894117523d-3*u824*(u852+4.19717071894426195d3*exppd2)
  u748=7.55388327504705771d-2*u783
  u701=-2.18223294612470556d-1*u824*(u852+u557)
  u547=3.6370549102078426d-2*u92
  u830=-9.0926372755196065d-2*u798
  u804=u830*u1306
  value_5_=value_5_+(c(66)*(z*(u1308*(u747+u1306*(u547*u1306+u701)+(u13&
 &06*(u822+u804)+u906)*u1308)+u1306*(u748+u93)+u55)))
  u557=1.65098547927777609d6*u686
  u852=u824*(u557-u574)
  u92=3.58312144923012958d-1*u852
  u93=3.98124605470014398d-2*u783
  u828=2.3349651778357119d7*u686
  u749=u824*(u828+u972)
  u758=-1.72520662370339572d-1*u749
  u763=5.75068874567798574d-2*u950
  value_6_=value_6_+(c(66)*(x*(u1308*(u92+u1306*(u946+u758)+u1308*(u667&
 &*u1308+u1306*(u763+u549)+u99))+u1306*(u93+u815)+u87)))
  u92=-1.47757899613393913d-3*u754
  u93=1.19191372354804423d-1*u847
  u758=9.75202137448399827d-2*u774
  u847=-1.69035037157722637d0*u552
  u763=+u847
  u667=-2.11293796447153296d-1*u554
  u946=u824*(u94-4.53748185831812103d2*u865)
  u865=1.40862530964768864d-1*u946
  u815=-9.39083539765125759d-2*u798
  u709=u815*u1306
  value_7_=value_7_+(c(66)*(z*(u1308*(u93+u1306*(u1306*(u865+u709)+u763&
 &)+u696*u1308)+u1306*(u758+u1306*(u76+u667))+u92)))
  u93=-4.92526332044646377d-3*u962
  u865=-u531
  u763=-2.3477088494128144d-2*u546
  u667=2.97978430887011058d-1*u765
  u758=(u763+u769)*u1307
  u696=u882*u1306
  u87=u881*u1306
  value_1_=value_1_+(c(67)*(y*(u1307*(u637+u1306*(u87+u912)+u1307*(u758&
 &+u1306*(u977+u794)+u865))+u1306*(u667+u696)+u93)))
  u667=(u549+u873)
  u531=u610*u1306
  u99=u898+u531
  u761=u879*u1306
  u75=u761+u551
  value_2_=value_2_+(c(67)*(u310*(u1307*(u99+u667*u1307)+u75)*z))
  u864=(y*(u1307*(u886+u1306*(u559*u1306+u911)+u1307*((u692+u524)*u1307&
 &+u925+u558))+u1306*(u970+u80*u1306)+u616))
  value_3_=value_3_+(c(67)*u864)
  u925=u928*u1307
  u524=u1307*(u951+u925)
  u912=u77*u1306
  u642=u827*u1306
  value_4_=value_4_+(c(67)*(z*(u1308*(u785+u1307*(u524+u883)+(u1307*(u7&
 &50+u925)+u760)*u1308)+u1307*(u630+u1306*(u912+u526)+u1307*((u582+u618&
 &)*u1307+u628+u773))+u1306*(u990+u642)+u543)))
  value_5_=value_5_+(c(67)*u759)
  u526=2.41287639678796605d-3*u868
  u628=-1.19437381641004319d-1*u766
  u858=-u550
  u984=1.32708201823338133d-2*u765
  u999=-u525
  u525=-5.17561987111018717d-1*u546
  u674=(u801+u617)
  u543=u903*u1306
  value_6_=value_6_+(c(67)*((u1307*(u628+u1306*(u706+u999)+u1307*(u674*&
 &u1307+u1306*(u525+u888)+u858))+u1306*(u984+u543)+u526)*z))
  u525=-u56
  u550=-2.11293796447153296d-1*u546
  u526=-8.2169809729448504d-1*u546
  u984=(u550+u73)
  u999=u904*u1306
  value_7_=value_7_+(c(67)*(x*(u1307*(u589+u1306*(u979+u607)+u1307*(u98&
 &4*u1307+u1306*(u526+u978)+u525))+u1306*(u637+u999)+u762)))
  u525=-1.97010532817858551d-3*u962
  u526=1.62533689574733305d-2*u765
  u858=-u96
  u56=7.04312654823844319d-2*u765
  u96=-u597
  u589=(u605+u794)
  value_1_=value_1_+(c(68)*((u1307*(u526+u1306*(u529+u96)+u1307*(u758+u&
 &1306*u589+u858))+u1306*(u56+u755)+u525)*z))
  u525=u967+u1306*(u790*u1306+u688)
  u526=u757+u689
  u56=u1306*(u74+u964*u1306)+u829
  u96=(u612+u816)
  value_2_=value_2_+(c(68)*(x*(u1307*(u525+u1307*(u96*u1307+u526))+u56)&
 &))
  u762=u668*u1307
  u628=u1307*(u8+u762)
  u757=u685*u1306
  u758=u583*u1306
  u597=(z*(u1308*(u857+u1307*(u628+u602)+(u1307*(u817+u762)+u995)*u1308&
 &)+u1307*(u68+u1306*(u757+u614)+u1307*((u570+u676)*u1307+u976+u673))+u&
 &1306*(u562+u758)+u712))
  value_3_=value_3_+(c(68)*u597)
  u759=u690*u1306
  value_4_=value_4_+(c(68)*(y*(u1308*(u649+u1307*(u750*u1307+u883)+u130&
 &8*((u522+u925)*u1308+u524+u650))+u1307*(u739+u1306*(u1306*(u66+u795)+&
 &u98)+u1307*((u651+u795)*u1307+u1306*(u90+u97)+u808))+u1306*(u738+u130&
 &6*(u759+u577))+u89)))
  value_5_=value_5_+(c(68)*u538)
  u784=-4.82575279357593209d-3*u771
  u97=8.84721345488920883d-3*u849
  u524=-5.75068874567798574d-2*u625
  u760=-7.96249210940028796d-2*u806
  u925=-2.3002754982711943d-1*u546
  u854=5.75068874567798574d-2*u853
  u853=-5.75068874567798574d-2*u542
  u542=(u670+u888)
  u849=u823*u1306
  value_6_=value_6_+(c(68)*(y*(u1307*(u97+u1306*(u1306*(u853+u617)+u925&
 &)+u1307*(u542*u1307+u854*u1306+u524))+u1306*(u760+u1306*(u849+u893))+&
 &u784)))
  u784=-1.62533689574733305d-1*u766
  u524=-4.6954176988256288d-1*u546
  value_7_=value_7_+(c(68)*(u310*(u1307*(u607+u1306*(u978+u524)+(u73+u5&
 &50)*u1307)+u1306*(u613+u878)+u784)*z))
  u524=3.94021065635717102d-3*u528
  u760=-5.32634315641177742d0*u819
  u878=+u760
  u853=u824*(2.35855068468253727d5*u686-u69)
  u69=-1.56032341991743972d0*u853
  u97=+u69
  u854=4.6954176988256288d-2*u646
  u752=u824*(u631+u578*u982)
  u982=2.3477088494128144d-2*u752
  u750=u824*(7.783217259452373d6*u686+u578)
  u522=2.11293796447153296d-1*u750
  u98=6.0d0*pd2
  u808=u824*(4.13925645161785292d8*u686+u578*(6.5d1+u98))
  u751=-2.3477088494128144d-2*u808
  u538=u550*u1306
  value_1_=value_1_+(c(69)*(y*(u1307*(u878+u1306*(u1306*(u751+u73)+u854&
 &)+u1307*((u998+u978)*u1307+u1306*(u982+u523)+u763))+u1306*(u97+u1306*&
 &(u538+u522))+u524)))
  u91=-4.77749526564017278d-1*u806
  u854=-4.6005509965423886d-1*u936
  value_2_=value_2_+(c(69)*(u310*(u1307*(u873+u1306*(u816+u854)+(u816+u&
 &612)*u1307)+u1306*(u873+u535)+u91)*z))
  u806=u885*u1306
  u854=(y*(u1308*(u598+u1307*(u817*u1307+u602)+u1308*((u952+u762)*u1308&
 &+u628+u666))+u1307*(u848+u1306*(u1306*(u809+u683)+u58)+u1307*((u67+u6&
 &83)*u1307+u1306*(u53+u663)+u679))+u1306*(u626+u1306*(u806+u825))+u665&
 &))
  value_3_=value_3_+(c(69)*u854)
  u91=-4.98401605398226757d-3*u771
  u751=6.73735038956467824d0*u819
  u524=2.96964289883110109d-2*u546
  u878=-2.96964289883110109d-2*u546
  u982=-5.34535721789598197d-1*u552
  u762=3.08385993340152806d-2*u824*(4.48124630089682082d6*u686+u574)
  u523=1.48482144941555055d-2*u824
  u596=u523*(u635+u972)
  u97=-5.93928579766220219d-2*u824*(u631+u578*(9.0d0+pd2))
  u635=-2.67267860894799098d-1*u732
  u631=u523*(3.50244776675356785d8*u686+u578*(5.5d1+u98))
  u523=-1.48482144941555055d-2*u798
  u732=u788*u1307
  u628=u1307*(u669+u732)
  value_4_=value_4_+(c(69)*(z*(u1308*(u751+u1306*(u1306*(u631+u795)+u59&
 &6)+u1307*(u628+u680)+u1308*((u878+u95)*u1308+u1307*(u992+u732)+u1306*&
 &(u97+u523*u1306)+u524))+u1307*(u647+u1307*(u690*u1307+u982))+u1306*(u&
 &762+u1306*(u759+u635))+u91)))
  value_5_=value_5_+(c(69)*u539)
  u524=3.98124605470014398d-2*u772
  u523=-6.90082649481358289d-1*u552
  u91=-3.98124605470014398d-2*u772
  u928=1.72520662370339572d-1*u555
  u751=6.90082649481358289d-1*u552
  u762=-5.75068874567798574d-2*u718
  u596=u1306*(u762+u617)
  u97=u655*u1307
  u635=u1307*(u62+u97)
  value_6_=value_6_+(c(69)*(z*(u1308*(u1306*(u596+u928)+u1307*(u635+u53&
 &6)+(u1307*(u873+u97)+u919)*u1308)+u1307*(u524+u1307*(u893*u1307+u523)&
 &)+u1306*(u91+u1306*(u849+u751)))))
  u919=-3.94021065635717102d-3*u771
  u523=-u69
  u69=-2.11293796447153296d-1*u750
  u751=-u760
  u524=-4.6954176988256288d-2*u646
  u853=2.3477088494128144d-2*u808
  u91=-2.3477088494128144d-2*u752
  u646=-4.6954176988256288d-2*u798
  u819=(u600+u794)
  u539=u646*u1306
  value_7_=value_7_+(c(69)*(x*(u1307*(u523+u1306*(u1306*(u91+u769)+u524&
 &)+u1307*(u819*u1307+u1306*(u853+u539)+u69))+u1306*(u751+u1306*(u763*u&
 &1306+u998))+u919)))
  u853=5.91031598453575652d-3*u528
  u523=-1.08355793049822203d-2*u776
  u751=1.62533689574733305d-2*u772
  u524=-2.81725061929537728d-1*u552
  u91=-4.87601068724199913d-2*u772
  u878=2.11293796447153296d-1*u555
  u752=8.45175185788613183d-1*u552
  u808=u880*u1307
  u919=u1307*(u63+u808)
  value_1_=value_1_+(c(70)*(z*(u1308*(u523+u1306*(u965+u878)+u1307*(u91&
 &9+u932)+(u1307*(u697+u808)+u648+u858)*u1308)+u1307*(u751+u1307*(u605*&
 &u1307+u524))+u1306*(u91+u1306*(u538+u752))+u853)))
  value_2_=value_2_+(c(70)*(x*(u1308*(u525+u1308*(u96*u1308+u526))+u56)&
 &))
  u853=1.0682259176834223d-2*u528
  u523=-2.79773454631372508d-3*u824*(-u818+u94)
  u751=u824*(u606+u972)
  u524=3.6370549102078426d-2*u751
  u91=-7.2741098204156852d-2*u546
  u752=1.25898054584117628d-2*u772
  u982=-2.18223294612470556d-1*u552
  u525=-1.25898054584117628d-2*u824*(-u574+2.61799125999761637d7*u686)
  u606=1.8185274551039213d-2*u824
  u526=u606*(2.35619213399785474d8*u686-1.54274383182816115d4*u553)
  u998=-3.6370549102078426d-2*u824*(3.24772429280785383d8*u686+u578*(5.&
 &1d1+u593))
  u648=2.47546060091623927d1*u945*pd2**2
  u750=+u713*(-u578*(u575-1.5d1)+u94)
  u94=5.4555823653117639d-2*u798
  u593=u565*u1307
  u96=u1307*(u907+u593)
  u553=(z*(u1308*(u523+u1306*(u1306*(u750+u683)+u526)+u1307*(u96+u708)+&
 &u1308*((u91+u753)*u1308+u1307*(u826+u593)+u1306*(u998+u94*u1306)+u524&
 &))+u1307*(u752+u1307*(u885*u1307+u982))+u1306*(u525+u1306*(u806+u648)&
 &)+u853))
  value_3_=value_3_+(c(70)*u553)
  u95=-4.0495130438605924d-3*u662*(-u574*(7.0d0*pd2-4.0d0)+u781)
  u753=1.02795331113384269d-2*u824*(3.353859073618568d7*u686-u574)
  u858=+u527*(9.48137375242379984d6*u686+u578)
  u527=1.6333035943571056d-1*u546
  u528=-2.74120882969024716d-2*u766
  u781=1.13465663335271798d5*u581
  value_4_=value_4_+(c(70)*(y*(u1308*(u753+u1306*(u912+u781)+u1307*(u99&
 &2*u1307+u680)+u1308*((u527+u867+u732)*u1308+u628+u1306*(u892+u618)+u8&
 &58))+u1307*(u621+u913*u1307)+u1306*(u528+u642)+u95)))
  value_5_=value_5_+(c(70)*u807)
  u807=-3.98124605470014398d-2*u766
  u628=3.98124605470014398d-2*u824*(u840-u574)
  u945=-1.72520662370339572d-1*u625
  u912=7.96249210940028796d-2*u765
  u760=-8.78901248927770365d5*u581
  u56=+u760
  u867=u873*u1307
  value_6_=value_6_+(c(70)*(y*(u1308*(u628+u1306*(u706+u56)+u1307*(u867&
 &+u536)+u1308*((u893+u97)*u1308+u635+u537+u945))+u1307*(u728+u903*u130&
 &7)+u1306*(u912+u543)+u807)))
  u807=u824*(u722+u818)
  u728=1.62533689574733305d-2*u807
  u56=2.81725061929537728d-1*u936
  value_7_=value_7_+(c(70)*(u545*((u899+u1306*(u978+u677))*u1308+u1307*&
 &(u550+u1306*(u794+u56)+(u794+u600)*u1307)+u1306*(u69+u76)+u728)))
  u728=u662*(2.26874092915906051d2*u656+u960)
  u56=-1.92085269497412087d-2*u728
  u656=2.6321425641057116d6*u686
  u642=8.12668447873666523d-2*u824
  u628=u642*(u656-u574)
  u903=-7.04312654823844319d-2*u742
  u62=9.75202137448399827d-2*u765
  u742=-1.07642979708394911d6*u581
  u969=+u742
  value_1_=value_1_+(c(71)*(y*(u1308*(u628+u1306*(u529+u969)+u1307*(u69&
 &7*u1307+u932)+u1308*((u605+u808)*u1308+u919+u629+u903))+u1307*(u729+u&
 &904*u1307)+u1306*(u62+u755)+u56)))
  value_2_=value_2_+(c(71)*(u545*(u1308*(u99+u667*u1308)+u75)))
  u56=-3.81509256315507965d-4*u662*(-u574*(4.7d1*pd2+9.1d1)+1.995333879&
 &24142653d7*u567)
  u628=1.25898054584117628d-2*u824*(u986+3.85685957957040288d3*exppd2)
  u62=+u756*(1.65570258064714117d7*u686+u578)
  u969=1.27296921857274491d-1*u546
  u629=5.87524254725882266d-2*u765
  u919=-1.11173191399135039d6*u581
  u635=9.82004825756117503d-1*u546
  u63=(y*(u1308*(u628+u1306*(u757+u919)+u1307*(u826*u1307+u708)+u1308*(&
 &(u969+u705+u593)*u1308+u96+u1306*(u635+u676)+u62))+u1307*(u568+u887*u&
 &1307)+u1306*(u629+u758)+u56))
  value_3_=value_3_+(c(71)*u63)
  u96=-1.55750501686945862d-3*u754
  u754=1.54192996670076403d-1*u852
  u912=-4.45446434824665164d-2*u591
  u529=7.42410724707775273d-2*u546
  u755=1.71325551855640448d-2*u807
  u75=-2.22723217412332582d-1*u555
  value_4_=value_4_+(c(71)*(z*(u1308*(u754+u75*u1306+u1308*((u529+u835)&
 &*u1308+u66*u1306+u912))+u755*u1306+u96)))
  u755=8.81286382088823399d-2*u824*(u530+1.58811865041134236d3*exppd2)
  u75=+u756*(u993+u578)
  u530=2.72779118265588195d-1*u546
  u756=1.25898054584117628d-2*u783
  u76=-5.4555823653117639d-2*u736
  u97=u606*(4.77606513648213798d8*u686-4.53748185831812103d2*u874)
  value_5_=value_5_+(c(71)*(x*(u1308*(u755+u76*u1306+u1308*((u530+u804)&
 &*u1308+u97*u1306+u75))+u756*u1306+u55)))
  u706=-7.84184828956088966d-2*u728
  u69=2.00139872385918163d6*u686
  u729=1.54826235460561155d-1*u824*(u69-u574)
  u757=u824*(u840+u578)
  u736=-5.75068874567798574d-2*u757
  u840=1.99062302735007199d-1*u824*(2.68874778053809249d6*u686+u660)
  u758=-2.87534437283899287d-1*u751
  u99=5.75068874567798574d-2*u544
  value_6_=value_6_+(c(71)*(z*(u1308*(u729+u1306*(u531+u758)+u1308*(u54&
 &2*u1308+u1306*(u99+u549)+u736))+u1306*(u840+u761)+u706)))
  u670=5.15511792509183147d6*u686
  u729=1.13773582702313313d-1*u824*(u670-u574)
  u706=3.62998548665449682d3*exppd2
  u840=5.41778965249111015d-3*u824*(u828+u706)
  u758=u824*(1.76891301351190296d7*u686+u972)
  u66=-2.11293796447153296d-1*u758
  u705=u824*(5.41287382134642304d8*u686+u578*(8.5d1+u98))
  u736=2.3477088494128144d-2*u705
  u783=5.63450123859075456d-1*u546
  value_7_=value_7_+(c(71)*(x*(u1308*(u729+u1306*(u783*u1306+u66)+u1308&
 &*(u819*u1308+u1306*(u736+u709)+u534))+u1306*(u840+u941*u1306)+u92)))
  u729=-4.46064681388434736d-1*u546
  u840=u609*u1307
  u728=u794+u840
  u736=u905*u1306
  u783=u633*u1306
  u66=u590*u1306
  value_1_=value_1_+(c(72)*(x*(u1307*(u637+u736+u1307*(u1307*(u729+u728&
 &)+u783+u607))+u66+u93)))
  u729=(u991+u549)
  u819=u551+u895*u1306
  u99=u698*u1306+u879
  u722=u586*u1306+u931
  value_2_=value_2_+(c(72)*((u1307*(u819+u1307*(u729*u1307+u99))+u722)*&
 &z))
  u92=u770*u1307
  u607=(u585+u1307*(u1307*(u694+u92)+u894))*u1308
  u941=u872*u1307
  u93=u676+u941
  u852=u949*u1306
  u986=u608*u1306
  u549=(x*(u607+u1307*(u884+u852+u1307*(u1307*(u682+u93)+u841+u844))+u9&
 &86+u687))
  value_3_=value_3_+(c(72)*u549)
  u751=u838*u1307
  u98=u751+u618
  u874=u7*u1307
  u835=u871*u1306
  u544=u799*u1306
  value_4_=value_4_+(c(72)*(u545*((u896+u1307*(u874+u779))*u1308+u1307*&
 &(u792+u835+u1307*(u98+u837))+u544+u796)))
  value_5_=value_5_+(c(72)*u864)
  u542=-2.65416403646676265d-1*u766
  u531=-u601
  u601=-u767
  u558=u702*u1307
  u761=u558+u888
  u543=u891*u1306
  u808=u889*u1306
  value_6_=value_6_+(c(72)*(u310*(u1307*(u531+u543+u1307*(u761+u601))+u&
 &808+u542)*z))
  u601=4.4327369884018174d-2*u868
  u531=-5.14690016986655464d-1*u766
  u767=-u634
  u542=-4.06334223936833262d-1*u766
  u563=-u955
  u945=-u775
  u775=(u977+u73)
  u634=u775*u1307
  u955=u899*u1306
  value_7_=value_7_+(c(72)*(y*(u1307*(u531+u1306*(u869+u563)+u1307*(u63&
 &4+u1306*(u945+u978)+u767))+u1306*(u542+u955)+u601)))
  u542=5.41778965249111015d-2*u765
  u563=-2.81725061929537728d-1*u546
  u945=u840+u794
  value_1_=value_1_+(c(73)*(u310*(u1307*(u613+u87+u1307*(u945+u563))+u6&
 &96+u542)*z))
  u542=u671+u1306*(u640*u1306+u672)
  u763=u800+u70
  u794=u1306*(u879+u556*u1306)+u59
  u800=(u576+u816)
  value_2_=value_2_+(c(73)*(y*(u1307*(u542+u1307*(u800*u1307+u763))+u79&
 &4)))
  u537=u941+u676
  u559=u791*u1306
  u676=(u545*((u897+u1307*(u92+u604))*u1308+u1307*(u684+u559+u1307*(u53&
 &7+u952))+u83+u619))
  value_3_=value_3_+(c(73)*u676)
  u616=u1306*(u992*u1306+u680)
  u816=u1306*(u621+u913*u1306)
  value_4_=value_4_+(c(73)*(x*(u1308*(u811+u1307*(u1307*(u533+u751)+u85&
 &9)+(u1307*(u892+u874)+u942)*u1308)+u1307*(u727+u616+u1307*((u909+u795&
 &)*u1307+u657+u532))+u816+u84)))
  value_5_=value_5_+(c(73)*u597)
  u597=u662*(u580-u574*(1.3d1+u845))
  u570=1.20643819839398302d-3*u597
  u536=-7.96249210940028796d-2*u776
  u614=1.72520662370339572d-1*u554
  u532=-1.32708201823338133d-2*u776
  u562=u1306*(u535+u928)
  u535=u1306*(u532+u948*u1306)
  value_6_=value_6_+(c(73)*(x*((u789+u1307*(u1307*(u703+u558)+u855))*u1&
 &308+u1307*(u536+u562+u1307*((u925+u617)*u1307+u596+u614))+u535+u570))&
 &)
  u570=5.91031598453575652d-3*u868
  u536=-2.11293796447153296d-1*u766
  u533=-u802
  u534=-4.87601068724199913d-2*u766
  u802=-u643
  u925=-9.15606451270997615d-1*u546
  value_7_=value_7_+(c(73)*((u1307*(u536+u1306*(u979+u802)+u1307*(u634+&
 &u1306*(u925+u978)+u533))+u1306*(u534+u999)+u570)*z))
  u634=4.92526332044646377d-4*u662*(u929+u541)
  u536=-2.60053903319573287d-1*u850
  u533=2.11293796447153296d-1*u745
  u534=-1.62533689574733305d-2*u776
  u643=u1306*(u579*u1306+u878)
  u925=u1306*(u534+u981*u1306)
  value_1_=value_1_+(c(74)*(x*((u637+u1307*(u1307*(u681+u840)+u934))*u1&
 &308+u1307*(u536+u643+u1307*((u599+u73)*u1307+u965+u533))+u925+u634)))
  u850=u664*u1307
  u536=u1307*(u9+u850)
  u929=u612*u1307
  value_2_=value_2_+(c(74)*(z*(u1308*(u74+u1307*(u536+u688)+(u1307*(u79&
 &0+u850)+u964)*u1308)+u1307*(u967+u1307*(u929+u689))+u829)))
  u634=u1306*(u826*u1306+u708)
  u802=u1306*(u568+u887*u1306)
  u745=(x*(u1308*(u623+u1307*(u1307*(u933+u941)+u916)+(u1307*(u72+u92)+&
 &u915)*u1308)+u1307*(u624+u634+u1307*((u71+u683)*u1307+u710+u842))+u80&
 &2+u994))
  value_3_=value_3_+(c(74)*u745)
  value_4_=value_4_+(c(74)*(u545*(u1308*(u917+u1307*(u751+u699)+(u874+u&
 &909)*u1308)+u1307*(u876+u1306*(u732+u615))+u1306*(u582+u759)+u743)))
  value_5_=value_5_+(c(74)*u854)
  u854=-3.98124605470014398d-2*u776
  u759=-2.3002754982711943d-1*u936
  value_6_=value_6_+(c(74)*(u545*((u875+u1307*(u558+u996))*u1308+u1307*&
 &(u614+u1306*(u617+u759)+(u617+u612)*u1307)+u1306*(u893+u849)+u854)))
  u854=-1.18206319690715131d-2*u771
  u759=1.30026951659786644d-1*u851
  u851=-1.40862530964768864d-1*u741
  u541=7.04312654823844319d-2*u645
  u645=2.3477088494128144d-2*u824*(u993+u901)
  value_7_=value_7_+(c(74)*(y*(u1307*(u759+u1306*(u1306*(u645+u769)+u85&
 &1)+u1307*(u589*u1307+u1306*(u54+u539)+u903))+u1306*(u759+u1306*(u88+u&
 &541))+u854)))
  u854=-1.62533689574733305d-2*u807
  u807=-2.81725061929537728d-1*u936
  value_1_=value_1_+(c(75)*(u545*((u613+u1307*(u840+u599))*u1308+u1307*&
 &(u522+u1306*(u73+u807)+(u73+u563)*u1307)+u1306*(u600+u538)+u854)))
  u522=u790*u1307
  u807=u964*u1307
  value_2_=value_2_+(c(75)*(y*(u1308*(u967+u1307*(u522+u688)+u1308*((u6&
 &12+u850)*u1308+u536+u689))+u1307*(u74+u807)+u829)))
  u688=(u545*(u1308*(u918+u1307*(u941+u822)+(u92+u926)*u1308)+u1307*(u8&
 &25+u1306*(u593+u654))+u1306*(u638+u806)+u627))
  value_3_=value_3_+(c(75)*u688)
  u654=u795+u874
  value_4_=value_4_+(c(75)*(x*(u1308*(u753+u616+u1307*(u77*u1307+u781)+&
 &u1308*((u527+u654)*u1308+u1307*(u892+u751)+u657+u858))+u1307*(u528+u8&
 &27*u1307)+u816+u95)))
  value_5_=value_5_+(c(75)*u553)
  u538=3.98124605470014398d-2*u765
  u648=-3.98124605470014398d-2*u846
  u541=1.72520662370339572d-1*u625
  u854=-7.96249210940028796d-2*u766
  u539=-u760
  u625=u1307*(u612+u558)
  value_6_=value_6_+(c(75)*(x*(u1308*(u648+u562+u1307*(u929+u539)+u1308&
 &*((u823+u617)*u1308+u625+u596+u541))+u1307*(u854+u948*u1307)+u535+u53&
 &8)))
  u538=-5.91031598453575652d-3*u771
  u648=1.08355793049822203d-2*u774
  u539=4.87601068724199913d-2*u772
  u854=-8.45175185788613183d-1*u552
  u535=-1.62533689574733305d-2*u772
  u772=7.04312654823844319d-2*u555
  u760=2.81725061929537728d-1*u552
  u555=-2.3477088494128144d-2*u718
  u718=u1306*(u555+u769)
  u762=u797*u1307
  u541=u1307*(u54+u762)
  value_7_=value_7_+(c(75)*(z*(u1308*(u648+u1306*(u718+u772)+u1307*(u54&
 &1+u834)+(u1307*(u603+u762)+u79+u944)*u1308)+u1307*(u539+u1307*(u600*u&
 &1307+u854))+u1306*(u535+u1306*(u88+u760))+u538)))
  u538=1.3160712820528558d7*u567
  u648=-u574*(3.1d1*pd2-1.3d1)
  u539=1.47757899613393913d-3*u662*(u538+u648)
  u854=u824*(3.26895124896999666d7*u686+u574)
  u944=-1.62533689574733305d-2*u854
  u760=u824*(8.63229550593808642d6*u686+u578)
  u535=2.11293796447153296d-1*u760
  u536=-3.25067379149466609d-2*u766
  u79=3.58809932361316369d5*u581
  u567=u1307*(u831+u840)
  value_1_=value_1_+(c(76)*(x*(u1308*(u944+u643+u1307*(u831*u1307+u79)+&
 &u1308*(u984*u1308+u567+u965+u535))+u1307*(u536+u821*u1307)+u925+u539)&
 &))
  value_2_=value_2_+(c(76)*(z*(u1308*(u542+u1308*(u800*u1308+u763))+u79&
 &4)))
  u59=u683+u92
  u683=(x*(u1308*(u628+u634+u1307*(u685*u1307+u919)+u1308*((u969+u59)*u&
 &1308+u1307*(u635+u941)+u710+u62))+u1307*(u629+u583*u1307)+u802+u56))
  value_3_=value_3_+(c(76)*u683)
  u535=-6.85302207422561791d-2*u766
  u536=3.56357147859732131d-1*u546
  u710=u871*u1307
  u542=u7*u1308
  u802=u799*u1307
  value_4_=value_4_+(c(76)*(u545*(u1308*(u896+u835+u710+u1308*(u542+u98&
 &+u536))+u802+u544+u535)))
  value_5_=value_5_+(c(76)*u63)
  u835=u996*u1307
  u63=u875*u1307
  value_6_=value_6_+(c(76)*(u545*(u1308*(u543+u835+(u761)*u1308)+u63+u8&
 &08)))
  u761=-1.47757899613393913d-3*u662*(u648+u538)
  u98=1.62533689574733305d-2*u854
  u648=-2.11293796447153296d-1*u760
  u925=3.25067379149466609d-2*u765
  u763=-u79
  value_7_=value_7_+(c(76)*(y*(u1308*(u98+u1306*(u979+u763)+u1307*(u603&
 &*u1307+u834)+u1308*((u600+u762)*u1308+u541+u778+u648))+u1307*(u744+u9&
 &00*u1307)+u1306*(u925+u999)+u761)))
  u600=u599*u1307
  u54=u613*u1307
  value_1_=value_1_+(c(77)*(u545*(u1308*(u922+u87+u600+(u945+u697)*u130&
 &8)+u54+u696+u737)))
  value_2_=value_2_+(c(77)*(y*(u1308*(u819+u1308*(u729*u1308+u99))+u722&
 &)))
  u648=4.61626200141764638d-1*u765
  u925=-2.77932978497837598d6*u581
  u697=1.41845141498105861d0*u546
  u729=u791*u1307
  u737=u770*u1308
  u87=u684*u1307
  u722=(u545*(u1308*(u925+u559+u729+u1308*(u737+u537+u697))+u87+u83+u64&
 &8))
  value_3_=value_3_+(c(77)*u722)
  u98=-3.11501003373891723d-3*u962
  u761=3.42651103711280895d-1*u765
  u945=-1.5128755111369573d6*u581
  u537=6.53321437742842241d-1*u546
  u763=u618+u751
  u618=u763+u542
  u543=u980*u1306
  u808=u975*u1306
  u559=u866*u1306
  u83=(u1308*(u761+u543+u980*u1307+u1308*(u1308*(u537+u618)+u975*u1307+&
 &u808+u945))+u866*u1307+u559+u98)
  value_4_=value_4_+(c(77)*(x*u83))
  u99=-1.90754628157753983d-3*u662*(-u574*(u845+1.56d2)+u580)
  u542=4.45766079404999545d7*u686
  u762=6.99433636578431269d-3*u824*(u542+1.79230533403565781d4*exppd2)
  u819=+u713*(u542+u578)
  u538=9.0926372755196065d-2*u546
  u845=2.30813100070882319d-1*u824*(1.35080630122727135d6*u686+u660)
  u541=-9.0926372755196065d-2*u824*(u542+u972)
  u660=u606*(6.68649119107499317d8*u686+u578*(1.05d2+u575))
  value_5_=value_5_+(c(77)*(z*(u1308*(u762+u541*u1306+u1308*((u538+u804&
 &)*u1308+u660*u1306+u819))+u845*u1306+u99)))
  u804=u888+u558
  u541=u914*u1306
  u606=u541+u855*u1307
  u845=u910*u1306
  u660=u703*u1307+u845
  u542=u622*u1306
  u713=u789*u1307+u542
  value_6_=value_6_+(c(77)*(x*(u1308*(u551+u606+u1308*((u991+u804)*u130&
 &8+u660+u879))+u713+u931)))
  u888=4.10387819134761486d6*u686
  u580=u642*(u888-u574)
  u540=-7.04312654823844319d-2*u554
  u642=8.12668447873666523d-2*u774
  u539=-3.5215632741192216d-1*u824*(1.57079475599856982d7*u686+u972)
  u972=7.04312654823844319d-2*u824*(1.84674518610642669d8*u686+u578*(2.&
 &9d1+u575))
  u544=9.39083539765125759d-1*u546
  value_7_=value_7_+(c(77)*(z*(u1308*(u580+u1306*(u544*u1306+u539)+u130&
 &8*(u589*u1308+u1306*(u972+u709)+u540))+u1306*(u642+u922*u1306)+u935))&
 &)
  u580=1.03430529729375739d-1*u868
  u539=-1.7066037405346997d0*u766
  u544=-7.74743920306228751d-1*u546
  u540=-u566
  value_1_=value_1_+(c(78)*(y*(u1307*(u539+u939*u1306+u1307*(u1307*(u54&
 &4+u728)+u540*u1306+u60))+u726*u1306+u580)))
  u60=u6*u1307
  value_2_=value_2_+(c(78)*(u310*(u1307*(u921+u1307*(u60+u997))+u644)*z&
 &))
  value_3_=value_3_+(c(78)*(y*((u721+u1307*(u1307*(u966+u92)+u938))*u13&
 &08+u1307*(u786+u764*u1306+u1307*(u1307*(u707+u93)+u768*u1306+u877))+u&
 &813*u1306+u920)))
  value_4_=value_4_+(c(78)*(z*((u720+u1307*(u1307*(u636+u874)+u937))*u1&
 &308+u1307*(u573+u543+u1307*(u1307*(u584+u763)+u808+u594))+u559+u695))&
 &)
  value_5_=value_5_+(c(78)*u549)
  u539=3.61931459518194907d-2*u868
  u544=-1.65885252279172666d0*u766
  u540=-u856
  u856=-u704
  value_6_=value_6_+(c(78)*((u1307*(u544+u541+u1307*(u1307*(u856+u804)+&
 &u845+u540))+u542+u539)*z))
  u539=-2.03167111968416631d0*u766
  u544=-u661
  u542=-u719
  u661=u814*u1307
  u719=u978+u661
  u845=u561*u1306
  u856=u587*u1306
  value_7_=value_7_+(c(78)*(x*(u1307*(u539+u696+u1307*(u1307*(u542+u719&
 &)+u845+u544))+u856+u601)))
  u539=9.85052664089292755d-3*u868
  u544=-u780
  u542=-5.869272123532036d-1*u546
  value_1_=value_1_+(c(79)*((u1307*(u531+u736+u1307*(u1307*(u542+u728)+&
 &u783+u544))+u66+u539)*z))
  value_2_=value_2_+(c(79)*(x*((u586+u1307*(u1307*(u698+u60)+u895))*u13&
 &08+u1307*(u551+u1307*(u991*u1307+u879))+u931)))
  value_3_=value_3_+(c(79)*(z*(u607+u1307*(u724+u852+u1307*(u1307*(u659&
 &+u93)+u841+u949))+u986+u82)))
  u585=u1306*(u636*u1306+u639)
  u607=u1306*(u799+u923*u1306)
  value_4_=value_4_+(c(79)*(y*(u1308*(u861+u1307*(u1307*(u968+u751)+u82&
 &0)+(u1307*(u779+u874)+u896)*u1308)+u1307*(u723+u585+u1307*(u1306*(u71&
 &4+u732)+u658))+u607+u81)))
  value_5_=value_5_+(c(79)*u676)
  u542=6.03219099196991512d-3*u662*(u777-u574*(3.9d1+u571))
  u539=-4.64478706381683464d-1*u766
  u544=-u717*(u2160+u641)
  u641=-u782
  u2160=5.75068874567798574d-2*u620
  u619=-u793
  u540=2.87534437283899287d-1*u572
  u571=-5.75068874567798574d-2*u548
  u782=u1306*(u996*u1306+u540)
  u717=u1306*(u571+u617)
  u777=u1306*(u889+u875*u1306)
  value_6_=value_6_+(c(79)*(y*((u539+u1307*(u1307*(u619+u558)+u641))*u1&
 &308+u1307*(u544+u782+u1307*((u576+u617)*u1307+u717+u2160))+u777+u542)&
 &))
  u542=-4.87601068724199913d-1*u766
  u544=-u927
  u2160=-u833
  u833=u661+u978
  value_7_=value_7_+(c(79)*(u310*(u1307*(u544+u869+u1307*(u833+u2160))+&
 &u955+u542)*z))
  u542=2.46263166022323189d-3*u662*(u961-u574*(2.6d1+u843))
  u544=-1.89622637837188855d-1*u766
  u2160=-5.41778965249111015d-2*u824*(-u818+u863)
  u310=-u803
  u543=7.04312654823844319d-2*u735
  u808=-4.93018858376691024d-1*u546
  u763=-9.39083539765125759d-2*u546
  u541=3.5215632741192216d-1*u572
  u803=-7.04312654823844319d-2*u548
  u548=u1306*(u678*u1306+u541)
  u572=u1306*(u803+u73)
  u735=u1306*(u882+u934*u1306)
  value_1_=value_1_+(c(80)*(y*((u544+u1307*(u1307*(u808+u840)+u310))*u1&
 &308+u1307*(u2160+u548+u1307*((u763+u73)*u1307+u572+u543))+u735+u542))&
 &)
  value_2_=value_2_+(c(80)*(u545*((u879+u1307*(u60+u610))*u1308+u1307*(&
 &u898+u867)+u551)))
  u867=u1306*(u930*u1306+u57)
  u73=u1306*(u684+u906*u1306)
  value_3_=value_3_+(c(80)*(y*(u1308*(u733+u1307*(u1307*(u885+u941)+u94&
 &3)+(u1307*(u604+u92)+u897)*u1308)+u1307*(u734+u867+u1307*(u1306*(u78+&
 &u593)+u679))+u73+u86)))
  u2160=u691*u1307
  value_4_=value_4_+(c(80)*(z*(u1308*(u730+u1307*(u85*u1307+u787)+(u130&
 &7*(u636+u2160)+u917)*u1308)+u1307*(u731+u876*u1307)+u812)))
  value_5_=value_5_+(c(80)*u745)
  u310=u662*(u716-u574*(2.6d1+u711))
  u543=3.61931459518194907d-3*u310
  u808=u824*(-u860+u902)
  u860=-3.71582965105346771d-1*u808
  u902=-2.22949779063208063d0*u808
  u544=3.45041324740679145d-1*u554
  u542=1.72520662370339572d-1*u591
  u763=-4.6005509965423886d-1*u836
  value_6_=value_6_+(c(80)*(z*(u1308*(u860+u1307*(u1307*(u763+u850)+u54&
 &4)+(u625+u948)*u1308)+u1307*(u902+u1307*(u929+u542))+u543)))
  u860=1.47757899613393913d-3*u310
  u763=-9.10188661618506505d-1*u808
  u544=7.04312654823844319d-2*u591
  u542=-5.41778965249111015d-3*u776
  u591=u1306*(u831*u1306+u772)
  u808=u1306*(u542+u821*u1306)
  value_7_=value_7_+(c(80)*(x*((u935+u1307*(u1307*(u693+u661)+u611))*u1&
 &308+u1307*(u763+u591+u1307*(u632*u1307+u718+u544))+u808+u860)))
  u632=1.47757899613393913d-3*u597
  u597=-1.19191372354804423d-1*u824*(-u862+u560)
  u544=-9.75202137448399827d-2*u776
  u776=-u847
  u763=2.11293796447153296d-1*u554
  u310=-1.40862530964768864d-1*u946
  u946=9.39083539765125759d-2*u798
  u798=u1307*(u579+u661)
  u552=u946*u1307
  value_1_=value_1_+(c(81)*(z*(u1308*(u597+u1307*(u1307*(u310+u552)+u77&
 &6)+(u798+u981)*u1308)+u1307*(u544+u1307*(u563*u1307+u763))+u632)))
  value_2_=value_2_+(c(81)*(x*(u1308*(u588+u1307*(u870*u1307+u924)+(u13&
 &07*(u870+u60)+u890)*u1308)+u1307*(u588+u890*u1307)+u65)))
  u588=u830*u1307
  value_3_=value_3_+(c(81)*(z*(u1308*(u747+u1307*(u547*u1307+u701)+(u13&
 &07*(u822+u588)+u906)*u1308)+u1307*(u748+u825*u1307)+u55)))
  value_4_=value_4_+(c(81)*(y*(u1308*(u647+u700*u1307+u1308*((u810+u216&
 &0)*u1308+u669*u1307+u577))+u746*u1307+u812)))
  value_5_=value_5_+(c(81)*u688)
  u544=-3.58312144923012958d-1*u824*(-u574+u557)
  u557=-3.98124605470014398d-2*u805
  u763=1.72520662370339572d-1*u749
  u310=-5.75068874567798574d-2*u950
  value_6_=value_6_+(c(81)*(y*(u1308*(u544+u1307*(u522+u763)+u1308*((u8&
 &23+u558)*u1308+u1307*(u310+u850)+u614))+u1307*(u557+u807)+u543)))
  u310=-4.87601068724199913d-2*u805
  u544=7.04312654823844319d-2*u620
  u763=-9.39083539765125759d-2*u936
  value_7_=value_7_+(c(81)*(u545*((u934+u1307*(u661+u678))*u1308+u1307*&
 &(u544+u1306*(u769+u763)+(u769+u563)*u1307)+u1306*(u605+u88)+u310)))
  u310=-1.13773582702313313d-1*u824*(-u574+u670)
  u544=-5.41778965249111015d-3*u824*(u706+u828)
  u706=2.11293796447153296d-1*u758
  u763=-2.3477088494128144d-2*u705
  u543=-5.63450123859075456d-1*u546
  value_1_=value_1_+(c(82)*(y*(u1308*(u310+u1307*(u543*u1307+u706)+u130&
 &8*((u550+u661)*u1308+u1307*(u763+u552)+u533))+u1307*(u544+u865*u1307)&
 &+u632)))
  value_2_=value_2_+(c(82)*(u545*(u1308*(u898+u610*u1307+(u60+u873)*u13&
 &08)+u879*u1307+u551)))
  value_3_=value_3_+(c(82)*(y*(u1308*(u755+u76*u1307+u1308*((u530+u588)&
 &*u1308+u97*u1307+u75))+u756*u1307+u55)))
  value_4_=value_4_+(c(82)*(z*(u1308*(u754+u585+u1307*(u710+u896)+u1308&
 &*((u529+u654)*u1308+u1307*(u536+u751)+u715+u912))+u1307*(u535+u802)+u&
 &607+u96)))
  value_5_=value_5_+(c(82)*u683)
  u310=u662*(u960+2.26874092915906051d2*u592)
  u544=7.84184828956088966d-2*u310
  u592=-1.54826235460561155d-1*u824*(-u574+u69)
  u763=5.75068874567798574d-2*u757
  u543=-3.98124605470014398d-1*u766
  u946=-u595
  value_6_=value_6_+(c(82)*(z*(u1308*(u592+u782+u1307*(u835+u946)+u1308&
 &*(u674*u1308+u625+u717+u763))+u1307*(u543+u63)+u777+u544)))
  u544=1.92085269497412087d-2*u310
  u310=-8.12668447873666523d-2*u824
  u763=+u310*(-u574+u656)
  u946=-9.75202137448399827d-2*u766
  u656=-u742
  value_7_=value_7_+(c(82)*(x*(u1308*(u763+u591+u1307*(u579*u1307+u656)&
 &+u1308*((u977+u769)*u1308+u798+u718+u653))+u1307*(u946+u981*u1307)+u8&
 &08+u544)))
  u763=+u310*(-u574+u888)
  u946=7.04312654823844319d-2*u554
  value_1_=value_1_+(c(83)*(z*(u1308*(u763+u548+u1307*(u600+u767)+u1308&
 &*(u775*u1308+u567+u572+u946))+u1307*(u784+u54)+u735+u590)))
  value_2_=value_2_+(c(83)*(x*(u1308*(u551+u895*u1307+u1308*((u991+u60)&
 &*u1308+u698*u1307+u879))+u586*u1307+u931)))
  value_3_=value_3_+(c(83)*(z*(u1308*(u762+u867+u1307*(u729+u925)+u1308&
 &*((u538+u59)*u1308+u1307*(u697+u941)+u953+u819))+u1307*(u648+u87)+u73&
 &+u99)))
  value_4_=value_4_+(c(83)*(y*u83))
  value_5_=value_5_+(c(83)*u722)
  u763=1.20643819839398302d-2*u868
  value_6_=value_6_+(c(83)*(y*(u1308*(u543+u606+u1308*((u576+u804)*u130&
 &8+u660+u556))+u713+u763)))
  u763=u934*u1307
  value_7_=value_7_+(c(83)*(u545*(u1308*(u767+u869+u678*u1307+(u833+u83&
 &1)*u1308)+u763+u955+u784)))
  value_1_=value_1_+(c(84)*(y*(u1308*(u736+u763+u1308*((u728)*u1308+u68&
 &1*u1307+u783))+u637*u1307+u66)))
  value_2_=value_2_+(c(84)*(u545*(u1308*(u921+u1308*(u6*u1308+u997))+u6&
 &44)))
  u644=-4.57811107578609558d-2*u962
  u763=2.09830090973529381d0*u765
  u946=-6.02188120078648129d6*u581
  u544=1.96400965151223501d0*u546
  u6=(u1308*(u763+u852+u949*u1307+u1308*(u1308*(u544+u93+u737)+u564*u13&
 &07+u841+u946))+u608*u1307+u986+u644)
  value_3_=value_3_+(c(84)*(y*u6))
  u737=-1.30830421417034524d-1*u962
  u962=2.15870195338106964d0*u765
  u765=-3.97129821673451292d6*u581
  u841=9.79982156614263361d-1*u546
  u545=-3.5978365889684494d-1*u766
  u766=1.98564910836725646d6*u581
  value_4_=value_4_+(c(84)*(z*(u1308*(u962+u766*u1306+u766*u1307+u1308*&
 &(u1308*(u841+u618)+u839*u1307+u839*u1306+u765))+u545*u1307+u545*u1306&
 &+u737)))
  value_5_=value_5_+(c(84)*(x*u6))
  value_6_=value_6_+(c(84)*(z*(u1308*(u940*u1306+u641*u1307+u1308*((u80&
 &4)*u1308+u619*u1307+u793*u1306))+u539*u1307+u725*u1306)))
  value_7_=value_7_+(c(84)*(x*(u1308*(u696+u611*u1307+u1308*((u719)*u13&
 &08+u693*u1307+u845))+u935*u1307+u856)))
  if ( lmax .eq. 6 ) go to 100

  ! set computation flag
  D_X_Y_shell_optv4_3 = .false.

  ! setting results
100  continue
  value(1:1+7*(NVEC-1):7) = value_1_
  value(2:2+7*(NVEC-1):7) = value_2_
  value(3:3+7*(NVEC-1):7) = value_3_
  value(4:4+7*(NVEC-1):7) = value_4_
  value(5:5+7*(NVEC-1):7) = value_5_
  value(6:6+7*(NVEC-1):7) = value_6_
  value(7:7+7*(NVEC-1):7) = value_7_

  contains

    include "Av_functions.h"

end function

  
end module
