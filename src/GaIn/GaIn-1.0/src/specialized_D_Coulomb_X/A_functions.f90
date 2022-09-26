
!> A functions defs.
!!
!! Author: I. Duchemin Octobre 2017
!!
module mod_A_functions

  implicit none

contains

!> A1(p,d) function defined as A0 = erf(d*p)/d
!!
recursive function A1_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A1_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A1_ =(1.5367449609015076d1*p*(1.42197314817798802d21*pd**8+1.21014227936723349d22*pd**6+3.33621490093070486d23*pd**4 &
+1.01717215836593919d24*pd**2+8.44762600642816689d24))/(1.15048620849450273d26+5.22024563447052628d25*pd**2 &
+1.04395636569385441d25*pd**4+1.16367164849543718d24*pd**6+7.35832271449147856d22*pd**8+2.20162444818591312d21*pd**10 &
+6.50278745612876125d18*pd**12)
  else
    A1_ =(p*erfpd)/pd
  end if
  
end function
  
!> A2(p,d) function defined as exp(d**2*p**2)*p**2*d**2*A2(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,2) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A2_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A2_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A2_ =(5.64189583547756287d-1*(1.82582646302925289d25*pd**7-1.15871850579653057d26*pd**5+1.10772110463889964d27*pd**3 &
-3.61493255600866731d27*pd))/((-2.71119941700650048d27)-7.95928821724725561d26*pd**2+1.65101226746749749d25*pd**4 &
+4.35337583276996732d25*pd**6+9.17318797690122307d24*pd**8+9.23279993542417062d23*pd**10+4.169055683850628d22*pd**12)
  else
    A2_ =-(5.64189583547756287d-1*(2.0d0*exppd2*pd-1.77245385090551603d0*erfpd))/pd**2
  end if
  
end function
  
!> A3(p,d) function defined as exp(d**2*p**2)*p**3*d**3*A3(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,3) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A3_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A3_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A3_ =-(5.64189583547756287d-1*(2.65274169003999195d30*pd**8-1.47840191711236527d32*pd**6+1.05510042934435257d33*pd**4 &
-9.75886020377926998d33*pd**2))/(1.82978628820861312d34+1.10915887536122898d34*pd**2+3.11702406542754856d33*pd**4 &
+5.26676727838812304d32*pd**6+5.739453775663442d31*pd**8+3.91277123132952498d30*pd**10+1.34704377000441565d29*pd**12)
  else
    A3_ =-(1.88063194515918762d-1*(2.0d0*exppd2*pd*(3.0d0+2.0d0*pd**2)-5.31736155271654808d0*erfpd))/pd**3
  end if
  
end function
  
!> A4(p,d) function defined as exp(d**2*p**2)*p**4*d**4*A4(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,4) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A4_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A4_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A4_ =-(5.64189583547756287d-1*(2.65274169003999195d30*pd**7-1.47840191711236527d32*pd**5+1.05510042934435257d33*pd**3 &
-9.75886020377926998d33*pd))/(1.82978628820861312d34+1.10915887536122898d34*pd**2+3.11702406542754856d33*pd**4 &
+5.26676727838812304d32*pd**6+5.739453775663442d31*pd**8+3.91277123132952498d30*pd**10+1.34704377000441565d29*pd**12)
  else
    A4_ =-(1.88063194515918762d-1*(2.0d0*exppd2*pd*(3.0d0+2.0d0*pd**2)-5.31736155271654808d0*erfpd))/pd**4
  end if
  
end function
  
!> A5(p,d) function defined as exp(d**2*p**2)*p**5*d**5*A5(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,5) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A5_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A5_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A5_ =-(5.64189583547756287d-1*(8.6971861193065171d31*pd**8-2.39972302190683709d33*pd**6+2.38326793983765042d34*pd**4 &
-1.38752674047482269d35*pd**2))/(9.10564423436602394d35+5.51814815232178275d35*pd**2+1.55212439246665151d35*pd**4 &
+2.6289912083424309d34*pd**6+2.87835632231077669d33*pd**8+1.97824276620716992d32*pd**10+6.90423175818625841d30*pd**12)
  else
    A5_ =-(3.76126389031837525d-2*(2.0d0*exppd2*pd*(1.5d1+1.0d1*pd**2+4.0d0*pd**4)-2.65868077635827404d1*erfpd))/pd**5
  end if
  
end function
  
!> A6(p,d) function defined as exp(d**2*p**2)*p**6*d**6*A6(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,6) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A6_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A6_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A6_ =-(5.64189583547756287d-1*(8.6971861193065171d31*pd**7-2.39972302190683709d33*pd**5+2.38326793983765042d34*pd**3 &
-1.38752674047482269d35*pd))/(9.10564423436602394d35+5.51814815232178275d35*pd**2+1.55212439246665151d35*pd**4 &
+2.6289912083424309d34*pd**6+2.87835632231077669d33*pd**8+1.97824276620716992d32*pd**10+6.90423175818625841d30*pd**12)
  else
    A6_ =-(3.76126389031837525d-2*(2.0d0*exppd2*pd*(1.5d1+1.0d1*pd**2+4.0d0*pd**4)-2.65868077635827404d1*erfpd))/pd**6
  end if
  
end function
  
!> A7(p,d) function defined as exp(d**2*p**2)*p**7*d**7*A7(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,7) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A7_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A7_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A7_ =-(5.64189583547756287d-1*(1.10211899535810335d31*pd**8-2.85391251683435624d32*pd**6+2.95001048565421505d33*pd**4 &
-1.40840046022945179d34*pd**2))/(4.15918260911509982d35+2.53179261773123038d35*pd**2+7.16029234246026915d34*pd**4 &
+1.2211591420851263d34*pd**6+1.34893338952549771d33*pd**8+9.38217298225845933d31*pd**10+3.3290719623770952d30*pd**12)
  else
    A7_ =-(5.37323412902625035d-3*(2.0d0*exppd2*pd*(1.05d2+7.0d1*pd**2+2.8d1*pd**4+8.0d0*pd**6)-1.86107654345079183d2*erfpd) &
)/pd**7
  end if
  
end function
  
!> A8(p,d) function defined as exp(d**2*p**2)*p**8*d**8*A8(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,8) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A8_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A8_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A8_ =-(5.64189583547756287d-1*(1.10211899535810335d31*pd**7-2.85391251683435624d32*pd**5+2.95001048565421505d33*pd**3 &
-1.40840046022945179d34*pd))/(4.15918260911509982d35+2.53179261773123038d35*pd**2+7.16029234246026915d34*pd**4 &
+1.2211591420851263d34*pd**6+1.34893338952549771d33*pd**8+9.38217298225845933d31*pd**10+3.3290719623770952d30*pd**12)
  else
    A8_ =-(5.37323412902625035d-3*(2.0d0*exppd2*pd*(1.05d2+7.0d1*pd**2+2.8d1*pd**4+8.0d0*pd**6)-1.86107654345079183d2*erfpd) &
)/pd**8
  end if
  
end function
  
!> A9(p,d) function defined as exp(d**2*p**2)*p**9*d**9*A9(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,9) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A9_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A9_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A9_ =-(5.64189583547756287d-1*(4.48394660940764824d32*pd**8-1.1356560762359877d34*pd**6+1.15043717987469144d35*pd**4 &
-4.9009320808152677d35*pd**2))/(7.96018577813667308d37+4.86698017401989612d37*pd**2+1.38393459734206249d37*pd**4 &
+2.37630682514844971d36*pd**6+2.647781585593613d35*pd**8+1.8625233880266517d34*pd**10+6.70916526463472235d32*pd**12)
  else
    A9_ =-(5.97026014336250039d-4*(2.0d0*exppd2*pd*(9.45d2+6.3d2*pd**2+2.52d2*pd**4+7.2d1*pd**6+1.6d1*pd**8) &
-1.67496888910571265d3*erfpd))/pd**9
  end if
  
end function
  
!> A10(p,d) function defined as exp(d**2*p**2)*p**10*d**10*A10(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,10) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A10_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A10_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A10_ =-(5.64189583547756287d-1*(4.48394660940764824d32*pd**7-1.1356560762359877d34*pd**5+1.15043717987469144d35*pd**3 &
-4.9009320808152677d35*pd))/(7.96018577813667308d37+4.86698017401989612d37*pd**2+1.38393459734206249d37*pd**4 &
+2.37630682514844971d36*pd**6+2.647781585593613d35*pd**8+1.8625233880266517d34*pd**10+6.70916526463472235d32*pd**12)
  else
    A10_ =-(5.97026014336250039d-4*(2.0d0*exppd2*pd*(9.45d2+6.3d2*pd**2+2.52d2*pd**4+7.2d1*pd**6+1.6d1*pd**8) &
-1.67496888910571265d3*erfpd))/pd**10
  end if
  
end function
  
!> A11(p,d) function defined as exp(d**2*p**2)*p**11*d**11*A11(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,11) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A11_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A11_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A11_ =-(5.64189583547756287d-1*(8.95585831463431248d34*pd**8-2.22573232738978417d36*pd**6+2.1904617771719888d37*pd**4 &
-8.67698884653933058d37*pd**2))/(9.16066318575853467d40+5.6266785193907037d40*pd**2+1.60882815621427956d40*pd**4 &
+2.78122606949093125d39*pd**6+3.12512465884678099d38*pd**8+2.22174723449297551d37*pd**10+8.1129217075017753d35*pd**12)
  else
    A11_ =-(5.42750922123863672d-5*(2.0d0*exppd2*pd*(1.0395d4+6.93d3*pd**2+2.772d3*pd**4+7.92d2*pd**6+1.76d2*pd**8+3.2d1*pd**10) &
-1.84246577801628391d4*erfpd))/pd**11
  end if
  
end function
  
!> A12(p,d) function defined as exp(d**2*p**2)*p**12*d**12*A12(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,12) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A12_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A12_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A12_ =-(5.64189583547756287d-1*(8.95585831463431248d34*pd**7-2.22573232738978417d36*pd**5+2.1904617771719888d37*pd**3 &
-8.67698884653933058d37*pd))/(9.16066318575853467d40+5.6266785193907037d40*pd**2+1.60882815621427956d40*pd**4 &
+2.78122606949093125d39*pd**6+3.12512465884678099d38*pd**8+2.22174723449297551d37*pd**10+8.1129217075017753d35*pd**12)
  else
    A12_ =-(5.42750922123863672d-5*(2.0d0*exppd2*pd*(1.0395d4+6.93d3*pd**2+2.772d3*pd**4+7.92d2*pd**6+1.76d2*pd**8+3.2d1*pd**10) &
-1.84246577801628391d4*erfpd))/pd**12
  end if
  
end function
  
!> A13(p,d) function defined as exp(d**2*p**2)*p**13*d**13*A13(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,13) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A13_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A13_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A13_ =-(5.64189583547756287d-1*(3.53717112978255505d39*pd**8-8.63925003222460506d40*pd**6+8.28121297589040641d41*pd**4 &
-3.12133761781311205d42*pd**2))/(2.47149585341704041d46+1.52501969269398265d46*pd**2+4.38421298201829905d45*pd**4 &
+7.6285805546202612d44*pd**6+8.63971613770652692d43*pd**8+6.20199405191995975d42*pd**10+2.29218461181251397d41*pd**12)
  else
    A13_ =-(4.17500709326048979d-6*(2.0d0*exppd2*pd*(1.35135d5+9.009d4*pd**2+3.6036d4*pd**4+1.0296d4*pd**6+2.288d3*pd**8 &
+4.16d2*pd**10+6.4d1*pd**12)-2.39520551142116908d5*erfpd))/pd**13
  end if
  
end function
  
!> A14(p,d) function defined as exp(d**2*p**2)*p**14*d**14*A14(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,14) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A14_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A14_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A14_ =-(5.64189583547756287d-1*(3.53717112978255505d39*pd**7-8.63925003222460506d40*pd**5+8.28121297589040641d41*pd**3 &
-3.12133761781311205d42*pd))/(2.47149585341704041d46+1.52501969269398265d46*pd**2+4.38421298201829905d45*pd**4 &
+7.6285805546202612d44*pd**6+8.63971613770652692d43*pd**8+6.20199405191995975d42*pd**10+2.29218461181251397d41*pd**12)
  else
    A14_ =-(4.17500709326048979d-6*(2.0d0*exppd2*pd*(1.35135d5+9.009d4*pd**2+3.6036d4*pd**4+1.0296d4*pd**6+2.288d3*pd**8 &
+4.16d2*pd**10+6.4d1*pd**12)-2.39520551142116908d5*erfpd))/pd**14
  end if
  
end function
  
!> A15(p,d) function defined as exp(d**2*p**2)*p**15*d**15*A15(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,15) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A15_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A15_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A15_ =-(5.64189583547756287d-1*(3.43825302446881269d39*pd**8-8.26997767393479215d40*pd**6+7.75449475587599174d41*pd**4 &
-2.82016827010081629d42*pd**2))/(1.89807376935388323d47+1.17637139287100748d47*pd**2+3.39934768478446547d46*pd**4 &
+5.95084504759202579d45*pd**6+6.78831838800430899d44*pd**8+4.9153323621316744d43*pd**10+1.83587153986720443d42*pd**12)
  else
    A15_ =-(2.78333806217365986d-7*(2.0d0*exppd2*pd*(2.027025d6+1.35135d6*pd**2+5.4054d5*pd**4+1.5444d5*pd**6+3.432d4*pd**8 &
+6.24d3*pd**10+9.6d2*pd**12+1.28d2*pd**14)-3.59280826713175363d6*erfpd))/pd**15
  end if
  
end function
  
!> A16(p,d) function defined as exp(d**2*p**2)*p**16*d**16*A16(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,16) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A16_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A16_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A16_ =-(5.64189583547756287d-1*(3.43825302446881269d39*pd**7-8.26997767393479215d40*pd**5+7.75449475587599174d41*pd**3 &
-2.82016827010081629d42*pd))/(1.89807376935388323d47+1.17637139287100748d47*pd**2+3.39934768478446547d46*pd**4 &
+5.95084504759202579d45*pd**6+6.78831838800430899d44*pd**8+4.9153323621316744d43*pd**10+1.83587153986720443d42*pd**12)
  else
    A16_ =-(2.78333806217365986d-7*(2.0d0*exppd2*pd*(2.027025d6+1.35135d6*pd**2+5.4054d5*pd**4+1.5444d5*pd**6+3.432d4*pd**8 &
+6.24d3*pd**10+9.6d2*pd**12+1.28d2*pd**14)-3.59280826713175363d6*erfpd))/pd**16
  end if
  
end function
  
!> A17(p,d) function defined as exp(d**2*p**2)*p**17*d**17*A17(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,17) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A17_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A17_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A17_ =-(5.64189583547756287d-1*(4.65320843743694271d42*pd**8-1.10455297187199966d44*pd**6+1.01714323266554697d45*pd**4 &
-3.60153853675944384d45*pd**2))/(2.30276561987242594d51+1.43310963983723162d51*pd**2+4.16104020747690795d50*pd**4 &
+7.32469348242190987d49*pd**6+8.40984295654744502d48*pd**8+6.13626951219897579d47*pd**10+2.31293573296737808d46*pd**12)
  else
    A17_ =-(1.63725768363156462d-8*(2.0d0*exppd2*pd*(3.4459425d7+2.297295d7*pd**2+9.18918d6*pd**4+2.62548d6*pd**6+5.8344d5*pd**8 &
+1.0608d5*pd**10+1.632d4*pd**12+2.176d3*pd**14+2.56d2*pd**16)-6.10777405412398116d7*erfpd))/pd**17
  end if
  
end function
  
!> A18(p,d) function defined as exp(d**2*p**2)*p**18*d**18*A18(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,18) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A18_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A18_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A18_ =-(5.64189583547756287d-1*(4.65320843743694271d42*pd**7-1.10455297187199966d44*pd**5+1.01714323266554697d45*pd**3 &
-3.60153853675944384d45*pd))/(2.30276561987242594d51+1.43310963983723162d51*pd**2+4.16104020747690795d50*pd**4 &
+7.32469348242190987d49*pd**6+8.40984295654744502d48*pd**8+6.13626951219897579d47*pd**10+2.31293573296737808d46*pd**12)
  else
    A18_ =-(1.63725768363156462d-8*(2.0d0*exppd2*pd*(3.4459425d7+2.297295d7*pd**2+9.18918d6*pd**4+2.62548d6*pd**6+5.8344d5*pd**8 &
+1.0608d5*pd**10+1.632d4*pd**12+2.176d3*pd**14+2.56d2*pd**16)-6.10777405412398116d7*erfpd))/pd**18
  end if
  
end function
  
!> A19(p,d) function defined as exp(d**2*p**2)*p**19*d**19*A19(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,19) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A19_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A19_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A19_ =-(5.64189583547756287d-1*(2.41773942854358771d45*pd**8-5.67470081412272937d46*pd**6+5.14852245007932701d47*pd**4 &
-1.78574969244499597d48*pd**2))/(1.19886851223813383d55+7.48971448615457741d54*pd**2+2.18415998723475595d54*pd**4 &
+3.86410192273039844d53*pd**6+4.46234578740070557d52*pd**8+3.27802313800833343d51*pd**10+1.24542734840601065d50*pd**12)
  else
    A19_ =-(8.61714570332402432d-10*(2.0d0*exppd2*pd*(6.54729075d8+4.3648605d8*pd**2+1.7459442d8*pd**4+4.988412d7*pd**6 &
+1.108536d7*pd**8+2.01552d6*pd**10+3.1008d5*pd**12+4.1344d4*pd**14+4.864d3*pd**16+5.12d2*pd**18)-1.16047707028355642d9*erfpd) &
)/pd**19
  end if
  
end function
  
!> A20(p,d) function defined as exp(d**2*p**2)*p**20*d**20*A20(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,20) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A20_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A20_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A20_ =-(5.64189583547756287d-1*(2.41773942854358771d45*pd**7-5.67470081412272937d46*pd**5+5.14852245007932701d47*pd**3 &
-1.78574969244499597d48*pd))/(1.19886851223813383d55+7.48971448615457741d54*pd**2+2.18415998723475595d54*pd**4 &
+3.86410192273039844d53*pd**6+4.46234578740070557d52*pd**8+3.27802313800833343d51*pd**10+1.24542734840601065d50*pd**12)
  else
    A20_ =-(8.61714570332402432d-10*(2.0d0*exppd2*pd*(6.54729075d8+4.3648605d8*pd**2+1.7459442d8*pd**4+4.988412d7*pd**6 &
+1.108536d7*pd**8+2.01552d6*pd**10+3.1008d5*pd**12+4.1344d4*pd**14+4.864d3*pd**16+5.12d2*pd**18)-1.16047707028355642d9*erfpd) &
)/pd**20
  end if
  
end function
  
!> A21(p,d) function defined as exp(d**2*p**2)*p**21*d**21*A21(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,21) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A21_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A21_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A21_ =-(5.64189583547756287d-1*(1.238359249657546d44*pd**8-2.87865498234057073d45*pd**6+2.57975635844756864d46*pd**4 &
-8.80264922212406924d46*pd**2))/(6.79613826872221131d54+4.26073079671927606d54*pd**2+1.24746913721343112d54*pd**4 &
+2.21695185444514265d53*pd**6+2.57344988337424055d52*pd**8+1.90173754335146702d51*pd**10+7.27543808750904725d49*pd**12)
  else
    A21_ =-(4.10340271586858301d-11*(2.0d0*exppd2*pd*(1.3749310575d10+9.16620705d9*pd**2+3.66648282d9*pd**4+1.04756652d9*pd**6 &
+2.3279256d8*pd**8+4.232592d7*pd**10+6.51168d6*pd**12+8.68224d5*pd**14+1.02144d5*pd**16+1.0752d4*pd**18+1.024d3*pd**20) &
-2.43700184759546848d10*erfpd))/pd**21
  end if
  
end function
  
!> A22(p,d) function defined as exp(d**2*p**2)*p**22*d**22*A22(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,22) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A22_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A22_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A22_ =-(5.64189583547756287d-1*(1.238359249657546d44*pd**7-2.87865498234057073d45*pd**5+2.57975635844756864d46*pd**3 &
-8.80264922212406924d46*pd))/(6.79613826872221131d54+4.26073079671927606d54*pd**2+1.24746913721343112d54*pd**4 &
+2.21695185444514265d53*pd**6+2.57344988337424055d52*pd**8+1.90173754335146702d51*pd**10+7.27543808750904725d49*pd**12)
  else
    A22_ =-(4.10340271586858301d-11*(2.0d0*exppd2*pd*(1.3749310575d10+9.16620705d9*pd**2+3.66648282d9*pd**4+1.04756652d9*pd**6 &
+2.3279256d8*pd**8+4.232592d7*pd**10+6.51168d6*pd**12+8.68224d5*pd**14+1.02144d5*pd**16+1.0752d4*pd**18+1.024d3*pd**20) &
-2.43700184759546848d10*erfpd))/pd**22
  end if
  
end function
  
!> A23(p,d) function defined as exp(d**2*p**2)*p**23*d**23*A23(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,23) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A23_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A23_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A23_ =-(5.64189583547756287d-1*(6.38413382959839544d46*pd**8-1.471793100201934d48*pd**6+1.3054373321495368d49*pd**4 &
-4.39586087812068545d49*pd**2))/(4.24231353310788303d58+2.66822965766413233d58*pd**2+7.84037936160313396d57*pd**4 &
+1.39904357352637024d57*pd**6+1.63152918201379481d56*pd**8+1.21203766681588607d55*pd**10+4.66499108910114236d53*pd**12)
  else
    A23_ =-(1.78408813733416653d-12*(2.0d0*exppd2*pd*(3.16234143225d11+2.1082276215d11*pd**2+8.432910486d10*pd**4 &
+2.409402996d10*pd**6+5.35422888d9*pd**8+9.7349616d8*pd**10+1.4976864d8*pd**12+1.9969152d7*pd**14+2.349312d6*pd**16 &
+2.47296d5*pd**18+2.3552d4*pd**20+2.048d3*pd**22)-5.60510424946957751d11*erfpd))/pd**23
  end if
  
end function
  
!> A24(p,d) function defined as exp(d**2*p**2)*p**24*d**24*A24(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,24) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A24_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A24_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A24_ =-(5.64189583547756287d-1*(6.38413382959839544d46*pd**7-1.471793100201934d48*pd**5+1.3054373321495368d49*pd**3 &
-4.39586087812068545d49*pd))/(4.24231353310788303d58+2.66822965766413233d58*pd**2+7.84037936160313396d57*pd**4 &
+1.39904357352637024d57*pd**6+1.63152918201379481d56*pd**8+1.21203766681588607d55*pd**10+4.66499108910114236d53*pd**12)
  else
    A24_ =-(1.78408813733416653d-12*(2.0d0*exppd2*pd*(3.16234143225d11+2.1082276215d11*pd**2+8.432910486d10*pd**4 &
+2.409402996d10*pd**6+5.35422888d9*pd**8+9.7349616d8*pd**10+1.4976864d8*pd**12+1.9969152d7*pd**14+2.349312d6*pd**16 &
+2.47296d5*pd**18+2.3552d4*pd**20+2.048d3*pd**22)-5.60510424946957751d11*erfpd))/pd**24
  end if
  
end function
  
!> A25(p,d) function defined as exp(d**2*p**2)*p**25*d**25*A25(p,d) + taylor(exp(d**2*p**2)*erf(d*p),d,0,25) = exp(d**2*p**2)*erf(d*p)
!!
recursive function A25_(p,pd,erfpd,exppd2)
  
  implicit none
  
  real(kind=8), intent(in)                :: p     
  real(kind=8), intent(in)                :: pd    
  real(kind=8), intent(in)                :: erfpd 
  real(kind=8), intent(in)                :: exppd2
  real(kind=8)                            :: A25_ !< result 
  
  ! test pd
  if ( abs(pd) < 1.2d0 ) then
    ! use order 20 padde approximant
    A25_ =-(5.64189583547756287d-1*(2.29589445564618788d46*pd**8-5.25519682398259314d47*pd**6+4.62070507741081858d48*pd**4 &
-1.53908703275995222d49*pd**2))/(2.00519111131986647d59+1.26489603882904571d59*pd**2+3.72900486070231213d58*pd**4 &
+6.67853458340932311d57*pd**6+7.82056479629040562d56*pd**8+5.83697908899783841d55*pd**10+2.25856065283650253d54*pd**12)
  else
    A25_ =-(7.13635254933666611d-14*(2.0d0*exppd2*pd*(7.905853580625d12+5.27056905375d12*pd**2+2.1082276215d12*pd**4 &
+6.02350749d11*pd**6+1.33855722d11*pd**8+2.4337404d10*pd**10+3.744216d9*pd**12+4.992288d8*pd**14+5.87328d7*pd**16+6.1824d6*pd**18 &
+5.888d5*pd**20+5.12d4*pd**22+4.096d3*pd**24)-1.40127606236739438d13*erfpd))/pd**25
  end if
  
end function

end module

