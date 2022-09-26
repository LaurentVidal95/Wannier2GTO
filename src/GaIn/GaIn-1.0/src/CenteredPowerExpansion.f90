
!> This module defines utility routines for (x-xi)^n power expansion:
!!
!! These routines are generated automatically through maxima scripting.
!!
!! Author: I. Duchemin July 2015
!!
module mod_CenteredPowerExpansion

  implicit none

contains

!!> expand the product (x-x1)^n1 * (x-x2)^n2 as sum_k=0^(n1+n2) ci(k+1)*(x-x3)^k
!!  
!!  uses the binomial expression of (x-xi)^ni = sum_k=0^ni  ni!/(k!(ni-k)!) * x^(ni-k) * (-xi)^k
!!  
recursive subroutine expand_centered_product(x1,n1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n1
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  select case(n1)
    case(0)
      call expand_centered_product_0(x1,x2,n2,x3,ci)
    case(1)
      call expand_centered_product_1(x1,x2,n2,x3,ci)
    case(2)
      call expand_centered_product_2(x1,x2,n2,x3,ci)
    case(3)
      call expand_centered_product_3(x1,x2,n2,x3,ci)
    case(4)
      call expand_centered_product_4(x1,x2,n2,x3,ci)
    case(5)
      call expand_centered_product_5(x1,x2,n2,x3,ci)
    case(6)
      call expand_centered_product_6(x1,x2,n2,x3,ci)
    case(7)
      call expand_centered_product_7(x1,x2,n2,x3,ci)
    case(8)
      call expand_centered_product_8(x1,x2,n2,x3,ci)
    case(9)
      call expand_centered_product_9(x1,x2,n2,x3,ci)
    case(10)
      call expand_centered_product_10(x1,x2,n2,x3,ci)
    case(11)
      call expand_centered_product_11(x1,x2,n2,x3,ci)
    case(12)
      call expand_centered_product_12(x1,x2,n2,x3,ci)
    case(13)
      call expand_centered_product_13(x1,x2,n2,x3,ci)
    case(14)
      call expand_centered_product_14(x1,x2,n2,x3,ci)
    case(15)
      call expand_centered_product_15(x1,x2,n2,x3,ci)
    case(16)
      call expand_centered_product_16(x1,x2,n2,x3,ci)
    case(17)
      call expand_centered_product_17(x1,x2,n2,x3,ci)
    case(18)
      call expand_centered_product_18(x1,x2,n2,x3,ci)
    case(19)
      call expand_centered_product_19(x1,x2,n2,x3,ci)
    case(20)
      call expand_centered_product_20(x1,x2,n2,x3,ci)
    case(21)
      call expand_centered_product_21(x1,x2,n2,x3,ci)
    case(22)
      call expand_centered_product_22(x1,x2,n2,x3,ci)
    case(23)
      call expand_centered_product_23(x1,x2,n2,x3,ci)
    case(24)
      call expand_centered_product_24(x1,x2,n2,x3,ci)
    case default
      print*,'Error: expand_centered_product for n1>24, here n1=',n1
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^0 * (x-x2)^n2 as sum_k=0^(0+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_0(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=1.0d0
    case(1)
      ci(1)=dx2
      ci(2)=1.0d0
    case(2)
      ci(1)=dx2**2
      ci(2)=2.0d0*dx2
      ci(3)=1.0d0
    case(3)
      ci(1)=dx2**3
      ci(2)=3.0d0*dx2**2
      ci(3)=3.0d0*dx2
      ci(4)=1.0d0
    case(4)
      ci(1)=dx2**4
      ci(2)=4.0d0*dx2**3
      ci(3)=6.0d0*dx2**2
      ci(4)=4.0d0*dx2
      ci(5)=1.0d0
    case(5)
      ci(1)=dx2**5
      ci(2)=5.0d0*dx2**4
      ci(3)=10.0d0*dx2**3
      ci(4)=10.0d0*dx2**2
      ci(5)=5.0d0*dx2
      ci(6)=1.0d0
    case(6)
      ci(1)=dx2**6
      ci(2)=6.0d0*dx2**5
      ci(3)=15.0d0*dx2**4
      ci(4)=20.0d0*dx2**3
      ci(5)=15.0d0*dx2**2
      ci(6)=6.0d0*dx2
      ci(7)=1.0d0
    case(7)
      ci(1)=dx2**7
      ci(2)=7.0d0*dx2**6
      ci(3)=21.0d0*dx2**5
      ci(4)=35.0d0*dx2**4
      ci(5)=35.0d0*dx2**3
      ci(6)=21.0d0*dx2**2
      ci(7)=7.0d0*dx2
      ci(8)=1.0d0
    case(8)
      ci(1)=dx2**8
      ci(2)=8.0d0*dx2**7
      ci(3)=28.0d0*dx2**6
      ci(4)=56.0d0*dx2**5
      ci(5)=70.0d0*dx2**4
      ci(6)=56.0d0*dx2**3
      ci(7)=28.0d0*dx2**2
      ci(8)=8.0d0*dx2
      ci(9)=1.0d0
    case(9)
      ci(1)=dx2**9
      ci(2)=9.0d0*dx2**8
      ci(3)=36.0d0*dx2**7
      ci(4)=84.0d0*dx2**6
      ci(5)=126.0d0*dx2**5
      ci(6)=126.0d0*dx2**4
      ci(7)=84.0d0*dx2**3
      ci(8)=36.0d0*dx2**2
      ci(9)=9.0d0*dx2
      ci(10)=1.0d0
    case(10)
      ci(1)=dx2**10
      ci(2)=10.0d0*dx2**9
      ci(3)=45.0d0*dx2**8
      ci(4)=120.0d0*dx2**7
      ci(5)=210.0d0*dx2**6
      ci(6)=252.0d0*dx2**5
      ci(7)=210.0d0*dx2**4
      ci(8)=120.0d0*dx2**3
      ci(9)=45.0d0*dx2**2
      ci(10)=10.0d0*dx2
      ci(11)=1.0d0
    case(11)
      ci(1)=dx2**11
      ci(2)=11.0d0*dx2**10
      ci(3)=55.0d0*dx2**9
      ci(4)=165.0d0*dx2**8
      ci(5)=330.0d0*dx2**7
      ci(6)=462.0d0*dx2**6
      ci(7)=462.0d0*dx2**5
      ci(8)=330.0d0*dx2**4
      ci(9)=165.0d0*dx2**3
      ci(10)=55.0d0*dx2**2
      ci(11)=11.0d0*dx2
      ci(12)=1.0d0
    case(12)
      ci(1)=dx2**12
      ci(2)=12.0d0*dx2**11
      ci(3)=66.0d0*dx2**10
      ci(4)=220.0d0*dx2**9
      ci(5)=495.0d0*dx2**8
      ci(6)=792.0d0*dx2**7
      ci(7)=924.0d0*dx2**6
      ci(8)=792.0d0*dx2**5
      ci(9)=495.0d0*dx2**4
      ci(10)=220.0d0*dx2**3
      ci(11)=66.0d0*dx2**2
      ci(12)=12.0d0*dx2
      ci(13)=1.0d0
    case(13)
      ci(1)=dx2**13
      ci(2)=13.0d0*dx2**12
      ci(3)=78.0d0*dx2**11
      ci(4)=286.0d0*dx2**10
      ci(5)=715.0d0*dx2**9
      ci(6)=1287.0d0*dx2**8
      ci(7)=1716.0d0*dx2**7
      ci(8)=1716.0d0*dx2**6
      ci(9)=1287.0d0*dx2**5
      ci(10)=715.0d0*dx2**4
      ci(11)=286.0d0*dx2**3
      ci(12)=78.0d0*dx2**2
      ci(13)=13.0d0*dx2
      ci(14)=1.0d0
    case(14)
      ci(1)=dx2**14
      ci(2)=14.0d0*dx2**13
      ci(3)=91.0d0*dx2**12
      ci(4)=364.0d0*dx2**11
      ci(5)=1001.0d0*dx2**10
      ci(6)=2002.0d0*dx2**9
      ci(7)=3003.0d0*dx2**8
      ci(8)=3432.0d0*dx2**7
      ci(9)=3003.0d0*dx2**6
      ci(10)=2002.0d0*dx2**5
      ci(11)=1001.0d0*dx2**4
      ci(12)=364.0d0*dx2**3
      ci(13)=91.0d0*dx2**2
      ci(14)=14.0d0*dx2
      ci(15)=1.0d0
    case(15)
      ci(1)=dx2**15
      ci(2)=15.0d0*dx2**14
      ci(3)=105.0d0*dx2**13
      ci(4)=455.0d0*dx2**12
      ci(5)=1365.0d0*dx2**11
      ci(6)=3003.0d0*dx2**10
      ci(7)=5005.0d0*dx2**9
      ci(8)=6435.0d0*dx2**8
      ci(9)=6435.0d0*dx2**7
      ci(10)=5005.0d0*dx2**6
      ci(11)=3003.0d0*dx2**5
      ci(12)=1365.0d0*dx2**4
      ci(13)=455.0d0*dx2**3
      ci(14)=105.0d0*dx2**2
      ci(15)=15.0d0*dx2
      ci(16)=1.0d0
    case(16)
      ci(1)=dx2**16
      ci(2)=16.0d0*dx2**15
      ci(3)=120.0d0*dx2**14
      ci(4)=560.0d0*dx2**13
      ci(5)=1820.0d0*dx2**12
      ci(6)=4368.0d0*dx2**11
      ci(7)=8008.0d0*dx2**10
      ci(8)=11440.0d0*dx2**9
      ci(9)=12870.0d0*dx2**8
      ci(10)=11440.0d0*dx2**7
      ci(11)=8008.0d0*dx2**6
      ci(12)=4368.0d0*dx2**5
      ci(13)=1820.0d0*dx2**4
      ci(14)=560.0d0*dx2**3
      ci(15)=120.0d0*dx2**2
      ci(16)=16.0d0*dx2
      ci(17)=1.0d0
    case(17)
      ci(1)=dx2**17
      ci(2)=17.0d0*dx2**16
      ci(3)=136.0d0*dx2**15
      ci(4)=680.0d0*dx2**14
      ci(5)=2380.0d0*dx2**13
      ci(6)=6188.0d0*dx2**12
      ci(7)=12376.0d0*dx2**11
      ci(8)=19448.0d0*dx2**10
      ci(9)=24310.0d0*dx2**9
      ci(10)=24310.0d0*dx2**8
      ci(11)=19448.0d0*dx2**7
      ci(12)=12376.0d0*dx2**6
      ci(13)=6188.0d0*dx2**5
      ci(14)=2380.0d0*dx2**4
      ci(15)=680.0d0*dx2**3
      ci(16)=136.0d0*dx2**2
      ci(17)=17.0d0*dx2
      ci(18)=1.0d0
    case(18)
      ci(1)=dx2**18
      ci(2)=18.0d0*dx2**17
      ci(3)=153.0d0*dx2**16
      ci(4)=816.0d0*dx2**15
      ci(5)=3060.0d0*dx2**14
      ci(6)=8568.0d0*dx2**13
      ci(7)=18564.0d0*dx2**12
      ci(8)=31824.0d0*dx2**11
      ci(9)=43758.0d0*dx2**10
      ci(10)=48620.0d0*dx2**9
      ci(11)=43758.0d0*dx2**8
      ci(12)=31824.0d0*dx2**7
      ci(13)=18564.0d0*dx2**6
      ci(14)=8568.0d0*dx2**5
      ci(15)=3060.0d0*dx2**4
      ci(16)=816.0d0*dx2**3
      ci(17)=153.0d0*dx2**2
      ci(18)=18.0d0*dx2
      ci(19)=1.0d0
    case(19)
      ci(1)=dx2**19
      ci(2)=19.0d0*dx2**18
      ci(3)=171.0d0*dx2**17
      ci(4)=969.0d0*dx2**16
      ci(5)=3876.0d0*dx2**15
      ci(6)=11628.0d0*dx2**14
      ci(7)=27132.0d0*dx2**13
      ci(8)=50388.0d0*dx2**12
      ci(9)=75582.0d0*dx2**11
      ci(10)=92378.0d0*dx2**10
      ci(11)=92378.0d0*dx2**9
      ci(12)=75582.0d0*dx2**8
      ci(13)=50388.0d0*dx2**7
      ci(14)=27132.0d0*dx2**6
      ci(15)=11628.0d0*dx2**5
      ci(16)=3876.0d0*dx2**4
      ci(17)=969.0d0*dx2**3
      ci(18)=171.0d0*dx2**2
      ci(19)=19.0d0*dx2
      ci(20)=1.0d0
    case(20)
      ci(1)=dx2**20
      ci(2)=20.0d0*dx2**19
      ci(3)=190.0d0*dx2**18
      ci(4)=1140.0d0*dx2**17
      ci(5)=4845.0d0*dx2**16
      ci(6)=15504.0d0*dx2**15
      ci(7)=38760.0d0*dx2**14
      ci(8)=77520.0d0*dx2**13
      ci(9)=125970.0d0*dx2**12
      ci(10)=167960.0d0*dx2**11
      ci(11)=184756.0d0*dx2**10
      ci(12)=167960.0d0*dx2**9
      ci(13)=125970.0d0*dx2**8
      ci(14)=77520.0d0*dx2**7
      ci(15)=38760.0d0*dx2**6
      ci(16)=15504.0d0*dx2**5
      ci(17)=4845.0d0*dx2**4
      ci(18)=1140.0d0*dx2**3
      ci(19)=190.0d0*dx2**2
      ci(20)=20.0d0*dx2
      ci(21)=1.0d0
    case(21)
      ci(1)=dx2**21
      ci(2)=21.0d0*dx2**20
      ci(3)=210.0d0*dx2**19
      ci(4)=1330.0d0*dx2**18
      ci(5)=5985.0d0*dx2**17
      ci(6)=20349.0d0*dx2**16
      ci(7)=54264.0d0*dx2**15
      ci(8)=116280.0d0*dx2**14
      ci(9)=203490.0d0*dx2**13
      ci(10)=293930.0d0*dx2**12
      ci(11)=352716.0d0*dx2**11
      ci(12)=352716.0d0*dx2**10
      ci(13)=293930.0d0*dx2**9
      ci(14)=203490.0d0*dx2**8
      ci(15)=116280.0d0*dx2**7
      ci(16)=54264.0d0*dx2**6
      ci(17)=20349.0d0*dx2**5
      ci(18)=5985.0d0*dx2**4
      ci(19)=1330.0d0*dx2**3
      ci(20)=210.0d0*dx2**2
      ci(21)=21.0d0*dx2
      ci(22)=1.0d0
    case(22)
      ci(1)=dx2**22
      ci(2)=22.0d0*dx2**21
      ci(3)=231.0d0*dx2**20
      ci(4)=1540.0d0*dx2**19
      ci(5)=7315.0d0*dx2**18
      ci(6)=26334.0d0*dx2**17
      ci(7)=74613.0d0*dx2**16
      ci(8)=170544.0d0*dx2**15
      ci(9)=319770.0d0*dx2**14
      ci(10)=497420.0d0*dx2**13
      ci(11)=646646.0d0*dx2**12
      ci(12)=705432.0d0*dx2**11
      ci(13)=646646.0d0*dx2**10
      ci(14)=497420.0d0*dx2**9
      ci(15)=319770.0d0*dx2**8
      ci(16)=170544.0d0*dx2**7
      ci(17)=74613.0d0*dx2**6
      ci(18)=26334.0d0*dx2**5
      ci(19)=7315.0d0*dx2**4
      ci(20)=1540.0d0*dx2**3
      ci(21)=231.0d0*dx2**2
      ci(22)=22.0d0*dx2
      ci(23)=1.0d0
    case(23)
      ci(1)=dx2**23
      ci(2)=23.0d0*dx2**22
      ci(3)=253.0d0*dx2**21
      ci(4)=1771.0d0*dx2**20
      ci(5)=8855.0d0*dx2**19
      ci(6)=33649.0d0*dx2**18
      ci(7)=100947.0d0*dx2**17
      ci(8)=245157.0d0*dx2**16
      ci(9)=490314.0d0*dx2**15
      ci(10)=817190.0d0*dx2**14
      ci(11)=1144066.0d0*dx2**13
      ci(12)=1352078.0d0*dx2**12
      ci(13)=1352078.0d0*dx2**11
      ci(14)=1144066.0d0*dx2**10
      ci(15)=817190.0d0*dx2**9
      ci(16)=490314.0d0*dx2**8
      ci(17)=245157.0d0*dx2**7
      ci(18)=100947.0d0*dx2**6
      ci(19)=33649.0d0*dx2**5
      ci(20)=8855.0d0*dx2**4
      ci(21)=1771.0d0*dx2**3
      ci(22)=253.0d0*dx2**2
      ci(23)=23.0d0*dx2
      ci(24)=1.0d0
    case(24)
      ci(1)=dx2**24
      ci(2)=24.0d0*dx2**23
      ci(3)=276.0d0*dx2**22
      ci(4)=2024.0d0*dx2**21
      ci(5)=10626.0d0*dx2**20
      ci(6)=42504.0d0*dx2**19
      ci(7)=134596.0d0*dx2**18
      ci(8)=346104.0d0*dx2**17
      ci(9)=735471.0d0*dx2**16
      ci(10)=1307504.0d0*dx2**15
      ci(11)=1961256.0d0*dx2**14
      ci(12)=2496144.0d0*dx2**13
      ci(13)=2704156.0d0*dx2**12
      ci(14)=2496144.0d0*dx2**11
      ci(15)=1961256.0d0*dx2**10
      ci(16)=1307504.0d0*dx2**9
      ci(17)=735471.0d0*dx2**8
      ci(18)=346104.0d0*dx2**7
      ci(19)=134596.0d0*dx2**6
      ci(20)=42504.0d0*dx2**5
      ci(21)=10626.0d0*dx2**4
      ci(22)=2024.0d0*dx2**3
      ci(23)=276.0d0*dx2**2
      ci(24)=24.0d0*dx2
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>24, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^1 * (x-x2)^n2 as sum_k=0^(1+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_1(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1
      ci(2)=1.0d0
    case(1)
      ci(1)=dx1*dx2
      ci(2)=dx2+dx1
      ci(3)=1.0d0
    case(2)
      ci(1)=dx1*dx2**2
      ci(2)=dx2*(dx2+2.0d0*dx1)
      ci(3)=2.0d0*dx2+dx1
      ci(4)=1.0d0
    case(3)
      ci(1)=dx1*dx2**3
      ci(2)=dx2**2*(dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx2*(dx2+dx1)
      ci(4)=3.0d0*dx2+dx1
      ci(5)=1.0d0
    case(4)
      ci(1)=dx1*dx2**4
      ci(2)=dx2**3*(dx2+4.0d0*dx1)
      ci(3)=2.0d0*dx2**2*(2.0d0*dx2+3.0d0*dx1)
      ci(4)=2.0d0*dx2*(3.0d0*dx2+2.0d0*dx1)
      ci(5)=4.0d0*dx2+dx1
      ci(6)=1.0d0
    case(5)
      ci(1)=dx1*dx2**5
      ci(2)=dx2**4*(dx2+5.0d0*dx1)
      ci(3)=5.0d0*dx2**3*(dx2+2.0d0*dx1)
      ci(4)=10.0d0*dx2**2*(dx2+dx1)
      ci(5)=5.0d0*dx2*(2.0d0*dx2+dx1)
      ci(6)=5.0d0*dx2+dx1
      ci(7)=1.0d0
    case(6)
      ci(1)=dx1*dx2**6
      ci(2)=dx2**5*(dx2+6.0d0*dx1)
      ci(3)=3.0d0*dx2**4*(2.0d0*dx2+5.0d0*dx1)
      ci(4)=5.0d0*dx2**3*(3.0d0*dx2+4.0d0*dx1)
      ci(5)=5.0d0*dx2**2*(4.0d0*dx2+3.0d0*dx1)
      ci(6)=3.0d0*dx2*(5.0d0*dx2+2.0d0*dx1)
      ci(7)=6.0d0*dx2+dx1
      ci(8)=1.0d0
    case(7)
      ci(1)=dx1*dx2**7
      ci(2)=dx2**6*(dx2+7.0d0*dx1)
      ci(3)=7.0d0*dx2**5*(dx2+3.0d0*dx1)
      ci(4)=7.0d0*dx2**4*(3.0d0*dx2+5.0d0*dx1)
      ci(5)=35.0d0*dx2**3*(dx2+dx1)
      ci(6)=7.0d0*dx2**2*(5.0d0*dx2+3.0d0*dx1)
      ci(7)=7.0d0*dx2*(3.0d0*dx2+dx1)
      ci(8)=7.0d0*dx2+dx1
      ci(9)=1.0d0
    case(8)
      ci(1)=dx1*dx2**8
      ci(2)=dx2**7*(dx2+8.0d0*dx1)
      ci(3)=4.0d0*dx2**6*(2.0d0*dx2+7.0d0*dx1)
      ci(4)=28.0d0*dx2**5*(dx2+2.0d0*dx1)
      ci(5)=14.0d0*dx2**4*(4.0d0*dx2+5.0d0*dx1)
      ci(6)=14.0d0*dx2**3*(5.0d0*dx2+4.0d0*dx1)
      ci(7)=28.0d0*dx2**2*(2.0d0*dx2+dx1)
      ci(8)=4.0d0*dx2*(7.0d0*dx2+2.0d0*dx1)
      ci(9)=8.0d0*dx2+dx1
      ci(10)=1.0d0
    case(9)
      ci(1)=dx1*dx2**9
      ci(2)=dx2**8*(dx2+9.0d0*dx1)
      ci(3)=9.0d0*dx2**7*(dx2+4.0d0*dx1)
      ci(4)=12.0d0*dx2**6*(3.0d0*dx2+7.0d0*dx1)
      ci(5)=42.0d0*dx2**5*(2.0d0*dx2+3.0d0*dx1)
      ci(6)=126.0d0*dx2**4*(dx2+dx1)
      ci(7)=42.0d0*dx2**3*(3.0d0*dx2+2.0d0*dx1)
      ci(8)=12.0d0*dx2**2*(7.0d0*dx2+3.0d0*dx1)
      ci(9)=9.0d0*dx2*(4.0d0*dx2+dx1)
      ci(10)=9.0d0*dx2+dx1
      ci(11)=1.0d0
    case(10)
      ci(1)=dx1*dx2**10
      ci(2)=dx2**9*(dx2+10.0d0*dx1)
      ci(3)=5.0d0*dx2**8*(2.0d0*dx2+9.0d0*dx1)
      ci(4)=15.0d0*dx2**7*(3.0d0*dx2+8.0d0*dx1)
      ci(5)=30.0d0*dx2**6*(4.0d0*dx2+7.0d0*dx1)
      ci(6)=42.0d0*dx2**5*(5.0d0*dx2+6.0d0*dx1)
      ci(7)=42.0d0*dx2**4*(6.0d0*dx2+5.0d0*dx1)
      ci(8)=30.0d0*dx2**3*(7.0d0*dx2+4.0d0*dx1)
      ci(9)=15.0d0*dx2**2*(8.0d0*dx2+3.0d0*dx1)
      ci(10)=5.0d0*dx2*(9.0d0*dx2+2.0d0*dx1)
      ci(11)=10.0d0*dx2+dx1
      ci(12)=1.0d0
    case(11)
      ci(1)=dx1*dx2**11
      ci(2)=dx2**10*(dx2+11.0d0*dx1)
      ci(3)=11.0d0*dx2**9*(dx2+5.0d0*dx1)
      ci(4)=55.0d0*dx2**8*(dx2+3.0d0*dx1)
      ci(5)=165.0d0*dx2**7*(dx2+2.0d0*dx1)
      ci(6)=66.0d0*dx2**6*(5.0d0*dx2+7.0d0*dx1)
      ci(7)=462.0d0*dx2**5*(dx2+dx1)
      ci(8)=66.0d0*dx2**4*(7.0d0*dx2+5.0d0*dx1)
      ci(9)=165.0d0*dx2**3*(2.0d0*dx2+dx1)
      ci(10)=55.0d0*dx2**2*(3.0d0*dx2+dx1)
      ci(11)=11.0d0*dx2*(5.0d0*dx2+dx1)
      ci(12)=11.0d0*dx2+dx1
      ci(13)=1.0d0
    case(12)
      ci(1)=dx1*dx2**12
      ci(2)=dx2**11*(dx2+12.0d0*dx1)
      ci(3)=6.0d0*dx2**10*(2.0d0*dx2+11.0d0*dx1)
      ci(4)=22.0d0*dx2**9*(3.0d0*dx2+10.0d0*dx1)
      ci(5)=55.0d0*dx2**8*(4.0d0*dx2+9.0d0*dx1)
      ci(6)=99.0d0*dx2**7*(5.0d0*dx2+8.0d0*dx1)
      ci(7)=132.0d0*dx2**6*(6.0d0*dx2+7.0d0*dx1)
      ci(8)=132.0d0*dx2**5*(7.0d0*dx2+6.0d0*dx1)
      ci(9)=99.0d0*dx2**4*(8.0d0*dx2+5.0d0*dx1)
      ci(10)=55.0d0*dx2**3*(9.0d0*dx2+4.0d0*dx1)
      ci(11)=22.0d0*dx2**2*(10.0d0*dx2+3.0d0*dx1)
      ci(12)=6.0d0*dx2*(11.0d0*dx2+2.0d0*dx1)
      ci(13)=12.0d0*dx2+dx1
      ci(14)=1.0d0
    case(13)
      ci(1)=dx1*dx2**13
      ci(2)=dx2**12*(dx2+13.0d0*dx1)
      ci(3)=13.0d0*dx2**11*(dx2+6.0d0*dx1)
      ci(4)=26.0d0*dx2**10*(3.0d0*dx2+11.0d0*dx1)
      ci(5)=143.0d0*dx2**9*(2.0d0*dx2+5.0d0*dx1)
      ci(6)=143.0d0*dx2**8*(5.0d0*dx2+9.0d0*dx1)
      ci(7)=429.0d0*dx2**7*(3.0d0*dx2+4.0d0*dx1)
      ci(8)=1716.0d0*dx2**6*(dx2+dx1)
      ci(9)=429.0d0*dx2**5*(4.0d0*dx2+3.0d0*dx1)
      ci(10)=143.0d0*dx2**4*(9.0d0*dx2+5.0d0*dx1)
      ci(11)=143.0d0*dx2**3*(5.0d0*dx2+2.0d0*dx1)
      ci(12)=26.0d0*dx2**2*(11.0d0*dx2+3.0d0*dx1)
      ci(13)=13.0d0*dx2*(6.0d0*dx2+dx1)
      ci(14)=13.0d0*dx2+dx1
      ci(15)=1.0d0
    case(14)
      ci(1)=dx1*dx2**14
      ci(2)=dx2**13*(dx2+14.0d0*dx1)
      ci(3)=7.0d0*dx2**12*(2.0d0*dx2+13.0d0*dx1)
      ci(4)=91.0d0*dx2**11*(dx2+4.0d0*dx1)
      ci(5)=91.0d0*dx2**10*(4.0d0*dx2+11.0d0*dx1)
      ci(6)=1001.0d0*dx2**9*(dx2+2.0d0*dx1)
      ci(7)=1001.0d0*dx2**8*(2.0d0*dx2+3.0d0*dx1)
      ci(8)=429.0d0*dx2**7*(7.0d0*dx2+8.0d0*dx1)
      ci(9)=429.0d0*dx2**6*(8.0d0*dx2+7.0d0*dx1)
      ci(10)=1001.0d0*dx2**5*(3.0d0*dx2+2.0d0*dx1)
      ci(11)=1001.0d0*dx2**4*(2.0d0*dx2+dx1)
      ci(12)=91.0d0*dx2**3*(11.0d0*dx2+4.0d0*dx1)
      ci(13)=91.0d0*dx2**2*(4.0d0*dx2+dx1)
      ci(14)=7.0d0*dx2*(13.0d0*dx2+2.0d0*dx1)
      ci(15)=14.0d0*dx2+dx1
      ci(16)=1.0d0
    case(15)
      ci(1)=dx1*dx2**15
      ci(2)=dx2**14*(dx2+15.0d0*dx1)
      ci(3)=15.0d0*dx2**13*(dx2+7.0d0*dx1)
      ci(4)=35.0d0*dx2**12*(3.0d0*dx2+13.0d0*dx1)
      ci(5)=455.0d0*dx2**11*(dx2+3.0d0*dx1)
      ci(6)=273.0d0*dx2**10*(5.0d0*dx2+11.0d0*dx1)
      ci(7)=1001.0d0*dx2**9*(3.0d0*dx2+5.0d0*dx1)
      ci(8)=715.0d0*dx2**8*(7.0d0*dx2+9.0d0*dx1)
      ci(9)=6435.0d0*dx2**7*(dx2+dx1)
      ci(10)=715.0d0*dx2**6*(9.0d0*dx2+7.0d0*dx1)
      ci(11)=1001.0d0*dx2**5*(5.0d0*dx2+3.0d0*dx1)
      ci(12)=273.0d0*dx2**4*(11.0d0*dx2+5.0d0*dx1)
      ci(13)=455.0d0*dx2**3*(3.0d0*dx2+dx1)
      ci(14)=35.0d0*dx2**2*(13.0d0*dx2+3.0d0*dx1)
      ci(15)=15.0d0*dx2*(7.0d0*dx2+dx1)
      ci(16)=15.0d0*dx2+dx1
      ci(17)=1.0d0
    case(16)
      ci(1)=dx1*dx2**16
      ci(2)=dx2**15*(dx2+16.0d0*dx1)
      ci(3)=8.0d0*dx2**14*(2.0d0*dx2+15.0d0*dx1)
      ci(4)=40.0d0*dx2**13*(3.0d0*dx2+14.0d0*dx1)
      ci(5)=140.0d0*dx2**12*(4.0d0*dx2+13.0d0*dx1)
      ci(6)=364.0d0*dx2**11*(5.0d0*dx2+12.0d0*dx1)
      ci(7)=728.0d0*dx2**10*(6.0d0*dx2+11.0d0*dx1)
      ci(8)=1144.0d0*dx2**9*(7.0d0*dx2+10.0d0*dx1)
      ci(9)=1430.0d0*dx2**8*(8.0d0*dx2+9.0d0*dx1)
      ci(10)=1430.0d0*dx2**7*(9.0d0*dx2+8.0d0*dx1)
      ci(11)=1144.0d0*dx2**6*(10.0d0*dx2+7.0d0*dx1)
      ci(12)=728.0d0*dx2**5*(11.0d0*dx2+6.0d0*dx1)
      ci(13)=364.0d0*dx2**4*(12.0d0*dx2+5.0d0*dx1)
      ci(14)=140.0d0*dx2**3*(13.0d0*dx2+4.0d0*dx1)
      ci(15)=40.0d0*dx2**2*(14.0d0*dx2+3.0d0*dx1)
      ci(16)=8.0d0*dx2*(15.0d0*dx2+2.0d0*dx1)
      ci(17)=16.0d0*dx2+dx1
      ci(18)=1.0d0
    case(17)
      ci(1)=dx1*dx2**17
      ci(2)=dx2**16*(dx2+17.0d0*dx1)
      ci(3)=17.0d0*dx2**15*(dx2+8.0d0*dx1)
      ci(4)=136.0d0*dx2**14*(dx2+5.0d0*dx1)
      ci(5)=340.0d0*dx2**13*(2.0d0*dx2+7.0d0*dx1)
      ci(6)=476.0d0*dx2**12*(5.0d0*dx2+13.0d0*dx1)
      ci(7)=6188.0d0*dx2**11*(dx2+2.0d0*dx1)
      ci(8)=1768.0d0*dx2**10*(7.0d0*dx2+11.0d0*dx1)
      ci(9)=4862.0d0*dx2**9*(4.0d0*dx2+5.0d0*dx1)
      ci(10)=24310.0d0*dx2**8*(dx2+dx1)
      ci(11)=4862.0d0*dx2**7*(5.0d0*dx2+4.0d0*dx1)
      ci(12)=1768.0d0*dx2**6*(11.0d0*dx2+7.0d0*dx1)
      ci(13)=6188.0d0*dx2**5*(2.0d0*dx2+dx1)
      ci(14)=476.0d0*dx2**4*(13.0d0*dx2+5.0d0*dx1)
      ci(15)=340.0d0*dx2**3*(7.0d0*dx2+2.0d0*dx1)
      ci(16)=136.0d0*dx2**2*(5.0d0*dx2+dx1)
      ci(17)=17.0d0*dx2*(8.0d0*dx2+dx1)
      ci(18)=17.0d0*dx2+dx1
      ci(19)=1.0d0
    case(18)
      ci(1)=dx1*dx2**18
      ci(2)=dx2**17*(dx2+18.0d0*dx1)
      ci(3)=9.0d0*dx2**16*(2.0d0*dx2+17.0d0*dx1)
      ci(4)=51.0d0*dx2**15*(3.0d0*dx2+16.0d0*dx1)
      ci(5)=204.0d0*dx2**14*(4.0d0*dx2+15.0d0*dx1)
      ci(6)=612.0d0*dx2**13*(5.0d0*dx2+14.0d0*dx1)
      ci(7)=1428.0d0*dx2**12*(6.0d0*dx2+13.0d0*dx1)
      ci(8)=2652.0d0*dx2**11*(7.0d0*dx2+12.0d0*dx1)
      ci(9)=3978.0d0*dx2**10*(8.0d0*dx2+11.0d0*dx1)
      ci(10)=4862.0d0*dx2**9*(9.0d0*dx2+10.0d0*dx1)
      ci(11)=4862.0d0*dx2**8*(10.0d0*dx2+9.0d0*dx1)
      ci(12)=3978.0d0*dx2**7*(11.0d0*dx2+8.0d0*dx1)
      ci(13)=2652.0d0*dx2**6*(12.0d0*dx2+7.0d0*dx1)
      ci(14)=1428.0d0*dx2**5*(13.0d0*dx2+6.0d0*dx1)
      ci(15)=612.0d0*dx2**4*(14.0d0*dx2+5.0d0*dx1)
      ci(16)=204.0d0*dx2**3*(15.0d0*dx2+4.0d0*dx1)
      ci(17)=51.0d0*dx2**2*(16.0d0*dx2+3.0d0*dx1)
      ci(18)=9.0d0*dx2*(17.0d0*dx2+2.0d0*dx1)
      ci(19)=18.0d0*dx2+dx1
      ci(20)=1.0d0
    case(19)
      ci(1)=dx1*dx2**19
      ci(2)=dx2**18*(dx2+19.0d0*dx1)
      ci(3)=19.0d0*dx2**17*(dx2+9.0d0*dx1)
      ci(4)=57.0d0*dx2**16*(3.0d0*dx2+17.0d0*dx1)
      ci(5)=969.0d0*dx2**15*(dx2+4.0d0*dx1)
      ci(6)=3876.0d0*dx2**14*(dx2+3.0d0*dx1)
      ci(7)=3876.0d0*dx2**13*(3.0d0*dx2+7.0d0*dx1)
      ci(8)=3876.0d0*dx2**12*(7.0d0*dx2+13.0d0*dx1)
      ci(9)=25194.0d0*dx2**11*(2.0d0*dx2+3.0d0*dx1)
      ci(10)=8398.0d0*dx2**10*(9.0d0*dx2+11.0d0*dx1)
      ci(11)=92378.0d0*dx2**9*(dx2+dx1)
      ci(12)=8398.0d0*dx2**8*(11.0d0*dx2+9.0d0*dx1)
      ci(13)=25194.0d0*dx2**7*(3.0d0*dx2+2.0d0*dx1)
      ci(14)=3876.0d0*dx2**6*(13.0d0*dx2+7.0d0*dx1)
      ci(15)=3876.0d0*dx2**5*(7.0d0*dx2+3.0d0*dx1)
      ci(16)=3876.0d0*dx2**4*(3.0d0*dx2+dx1)
      ci(17)=969.0d0*dx2**3*(4.0d0*dx2+dx1)
      ci(18)=57.0d0*dx2**2*(17.0d0*dx2+3.0d0*dx1)
      ci(19)=19.0d0*dx2*(9.0d0*dx2+dx1)
      ci(20)=19.0d0*dx2+dx1
      ci(21)=1.0d0
    case(20)
      ci(1)=dx1*dx2**20
      ci(2)=dx2**19*(dx2+20.0d0*dx1)
      ci(3)=10.0d0*dx2**18*(2.0d0*dx2+19.0d0*dx1)
      ci(4)=190.0d0*dx2**17*(dx2+6.0d0*dx1)
      ci(5)=285.0d0*dx2**16*(4.0d0*dx2+17.0d0*dx1)
      ci(6)=969.0d0*dx2**15*(5.0d0*dx2+16.0d0*dx1)
      ci(7)=7752.0d0*dx2**14*(2.0d0*dx2+5.0d0*dx1)
      ci(8)=38760.0d0*dx2**13*(dx2+2.0d0*dx1)
      ci(9)=9690.0d0*dx2**12*(8.0d0*dx2+13.0d0*dx1)
      ci(10)=41990.0d0*dx2**11*(3.0d0*dx2+4.0d0*dx1)
      ci(11)=16796.0d0*dx2**10*(10.0d0*dx2+11.0d0*dx1)
      ci(12)=16796.0d0*dx2**9*(11.0d0*dx2+10.0d0*dx1)
      ci(13)=41990.0d0*dx2**8*(4.0d0*dx2+3.0d0*dx1)
      ci(14)=9690.0d0*dx2**7*(13.0d0*dx2+8.0d0*dx1)
      ci(15)=38760.0d0*dx2**6*(2.0d0*dx2+dx1)
      ci(16)=7752.0d0*dx2**5*(5.0d0*dx2+2.0d0*dx1)
      ci(17)=969.0d0*dx2**4*(16.0d0*dx2+5.0d0*dx1)
      ci(18)=285.0d0*dx2**3*(17.0d0*dx2+4.0d0*dx1)
      ci(19)=190.0d0*dx2**2*(6.0d0*dx2+dx1)
      ci(20)=10.0d0*dx2*(19.0d0*dx2+2.0d0*dx1)
      ci(21)=20.0d0*dx2+dx1
      ci(22)=1.0d0
    case(21)
      ci(1)=dx1*dx2**21
      ci(2)=dx2**20*(dx2+21.0d0*dx1)
      ci(3)=21.0d0*dx2**19*(dx2+10.0d0*dx1)
      ci(4)=70.0d0*dx2**18*(3.0d0*dx2+19.0d0*dx1)
      ci(5)=665.0d0*dx2**17*(2.0d0*dx2+9.0d0*dx1)
      ci(6)=1197.0d0*dx2**16*(5.0d0*dx2+17.0d0*dx1)
      ci(7)=6783.0d0*dx2**15*(3.0d0*dx2+8.0d0*dx1)
      ci(8)=7752.0d0*dx2**14*(7.0d0*dx2+15.0d0*dx1)
      ci(9)=29070.0d0*dx2**13*(4.0d0*dx2+7.0d0*dx1)
      ci(10)=22610.0d0*dx2**12*(9.0d0*dx2+13.0d0*dx1)
      ci(11)=58786.0d0*dx2**11*(5.0d0*dx2+6.0d0*dx1)
      ci(12)=352716.0d0*dx2**10*(dx2+dx1)
      ci(13)=58786.0d0*dx2**9*(6.0d0*dx2+5.0d0*dx1)
      ci(14)=22610.0d0*dx2**8*(13.0d0*dx2+9.0d0*dx1)
      ci(15)=29070.0d0*dx2**7*(7.0d0*dx2+4.0d0*dx1)
      ci(16)=7752.0d0*dx2**6*(15.0d0*dx2+7.0d0*dx1)
      ci(17)=6783.0d0*dx2**5*(8.0d0*dx2+3.0d0*dx1)
      ci(18)=1197.0d0*dx2**4*(17.0d0*dx2+5.0d0*dx1)
      ci(19)=665.0d0*dx2**3*(9.0d0*dx2+2.0d0*dx1)
      ci(20)=70.0d0*dx2**2*(19.0d0*dx2+3.0d0*dx1)
      ci(21)=21.0d0*dx2*(10.0d0*dx2+dx1)
      ci(22)=21.0d0*dx2+dx1
      ci(23)=1.0d0
    case(22)
      ci(1)=dx1*dx2**22
      ci(2)=dx2**21*(dx2+22.0d0*dx1)
      ci(3)=11.0d0*dx2**20*(2.0d0*dx2+21.0d0*dx1)
      ci(4)=77.0d0*dx2**19*(3.0d0*dx2+20.0d0*dx1)
      ci(5)=385.0d0*dx2**18*(4.0d0*dx2+19.0d0*dx1)
      ci(6)=1463.0d0*dx2**17*(5.0d0*dx2+18.0d0*dx1)
      ci(7)=4389.0d0*dx2**16*(6.0d0*dx2+17.0d0*dx1)
      ci(8)=10659.0d0*dx2**15*(7.0d0*dx2+16.0d0*dx1)
      ci(9)=21318.0d0*dx2**14*(8.0d0*dx2+15.0d0*dx1)
      ci(10)=35530.0d0*dx2**13*(9.0d0*dx2+14.0d0*dx1)
      ci(11)=49742.0d0*dx2**12*(10.0d0*dx2+13.0d0*dx1)
      ci(12)=58786.0d0*dx2**11*(11.0d0*dx2+12.0d0*dx1)
      ci(13)=58786.0d0*dx2**10*(12.0d0*dx2+11.0d0*dx1)
      ci(14)=49742.0d0*dx2**9*(13.0d0*dx2+10.0d0*dx1)
      ci(15)=35530.0d0*dx2**8*(14.0d0*dx2+9.0d0*dx1)
      ci(16)=21318.0d0*dx2**7*(15.0d0*dx2+8.0d0*dx1)
      ci(17)=10659.0d0*dx2**6*(16.0d0*dx2+7.0d0*dx1)
      ci(18)=4389.0d0*dx2**5*(17.0d0*dx2+6.0d0*dx1)
      ci(19)=1463.0d0*dx2**4*(18.0d0*dx2+5.0d0*dx1)
      ci(20)=385.0d0*dx2**3*(19.0d0*dx2+4.0d0*dx1)
      ci(21)=77.0d0*dx2**2*(20.0d0*dx2+3.0d0*dx1)
      ci(22)=11.0d0*dx2*(21.0d0*dx2+2.0d0*dx1)
      ci(23)=22.0d0*dx2+dx1
      ci(24)=1.0d0
    case(23)
      ci(1)=dx1*dx2**23
      ci(2)=dx2**22*(dx2+23.0d0*dx1)
      ci(3)=23.0d0*dx2**21*(dx2+11.0d0*dx1)
      ci(4)=253.0d0*dx2**20*(dx2+7.0d0*dx1)
      ci(5)=1771.0d0*dx2**19*(dx2+5.0d0*dx1)
      ci(6)=1771.0d0*dx2**18*(5.0d0*dx2+19.0d0*dx1)
      ci(7)=33649.0d0*dx2**17*(dx2+3.0d0*dx1)
      ci(8)=14421.0d0*dx2**16*(7.0d0*dx2+17.0d0*dx1)
      ci(9)=245157.0d0*dx2**15*(dx2+2.0d0*dx1)
      ci(10)=163438.0d0*dx2**14*(3.0d0*dx2+5.0d0*dx1)
      ci(11)=163438.0d0*dx2**13*(5.0d0*dx2+7.0d0*dx1)
      ci(12)=104006.0d0*dx2**12*(11.0d0*dx2+13.0d0*dx1)
      ci(13)=1352078.0d0*dx2**11*(dx2+dx1)
      ci(14)=104006.0d0*dx2**10*(13.0d0*dx2+11.0d0*dx1)
      ci(15)=163438.0d0*dx2**9*(7.0d0*dx2+5.0d0*dx1)
      ci(16)=163438.0d0*dx2**8*(5.0d0*dx2+3.0d0*dx1)
      ci(17)=245157.0d0*dx2**7*(2.0d0*dx2+dx1)
      ci(18)=14421.0d0*dx2**6*(17.0d0*dx2+7.0d0*dx1)
      ci(19)=33649.0d0*dx2**5*(3.0d0*dx2+dx1)
      ci(20)=1771.0d0*dx2**4*(19.0d0*dx2+5.0d0*dx1)
      ci(21)=1771.0d0*dx2**3*(5.0d0*dx2+dx1)
      ci(22)=253.0d0*dx2**2*(7.0d0*dx2+dx1)
      ci(23)=23.0d0*dx2*(11.0d0*dx2+dx1)
      ci(24)=23.0d0*dx2+dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>23, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^2 * (x-x2)^n2 as sum_k=0^(2+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_2(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**2
      ci(2)=2.0d0*dx1
      ci(3)=1.0d0
    case(1)
      ci(1)=dx1**2*dx2
      ci(2)=dx1*(2.0d0*dx2+dx1)
      ci(3)=dx2+2.0d0*dx1
      ci(4)=1.0d0
    case(2)
      ci(1)=dx1**2*dx2**2
      ci(2)=2.0d0*dx1*dx2*(dx2+dx1)
      ci(3)=dx2**2+dx1*(4.0d0*dx2+dx1)
      ci(4)=2.0d0*(dx2+dx1)
      ci(5)=1.0d0
    case(3)
      ci(1)=dx1**2*dx2**3
      ci(2)=dx1*dx2**2*(2.0d0*dx2+3.0d0*dx1)
      ci(3)=dx2*(dx2**2+3.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=3.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)
      ci(5)=3.0d0*dx2+2.0d0*dx1
      ci(6)=1.0d0
    case(4)
      ci(1)=dx1**2*dx2**4
      ci(2)=2.0d0*dx1*dx2**3*(dx2+2.0d0*dx1)
      ci(3)=dx2**2*(dx2**2+2.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(4)=4.0d0*dx2*(dx2**2+dx1*(3.0d0*dx2+dx1))
      ci(5)=6.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)
      ci(6)=2.0d0*(2.0d0*dx2+dx1)
      ci(7)=1.0d0
    case(5)
      ci(1)=dx1**2*dx2**5
      ci(2)=dx1*dx2**4*(2.0d0*dx2+5.0d0*dx1)
      ci(3)=dx2**3*(dx2**2+10.0d0*dx1*(dx2+dx1))
      ci(4)=5.0d0*dx2**2*(dx2**2+2.0d0*dx1*(2.0d0*dx2+dx1))
      ci(5)=5.0d0*dx2*(2.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))
      ci(6)=10.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)
      ci(7)=5.0d0*dx2+2.0d0*dx1
      ci(8)=1.0d0
    case(6)
      ci(1)=dx1**2*dx2**6
      ci(2)=2.0d0*dx1*dx2**5*(dx2+3.0d0*dx1)
      ci(3)=dx2**4*(dx2**2+3.0d0*dx1*(4.0d0*dx2+5.0d0*dx1))
      ci(4)=2.0d0*dx2**3*(3.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(5)=5.0d0*dx2**2*(3.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(6)=2.0d0*dx2*(10.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+dx1))
      ci(7)=15.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)
      ci(8)=2.0d0*(3.0d0*dx2+dx1)
      ci(9)=1.0d0
    case(7)
      ci(1)=dx1**2*dx2**7
      ci(2)=dx1*dx2**6*(2.0d0*dx2+7.0d0*dx1)
      ci(3)=dx2**5*(dx2**2+7.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(4)=7.0d0*dx2**4*(dx2**2+dx1*(6.0d0*dx2+5.0d0*dx1))
      ci(5)=7.0d0*dx2**3*(3.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))
      ci(6)=7.0d0*dx2**2*(5.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1))
      ci(7)=7.0d0*dx2*(5.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(8)=21.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)
      ci(9)=7.0d0*dx2+2.0d0*dx1
      ci(10)=1.0d0
    case(8)
      ci(1)=dx1**2*dx2**8
      ci(2)=2.0d0*dx1*dx2**7*(dx2+4.0d0*dx1)
      ci(3)=dx2**6*(dx2**2+4.0d0*dx1*(4.0d0*dx2+7.0d0*dx1))
      ci(4)=8.0d0*dx2**5*(dx2**2+7.0d0*dx1*(dx2+dx1))
      ci(5)=14.0d0*dx2**4*(2.0d0*dx2**2+dx1*(8.0d0*dx2+5.0d0*dx1))
      ci(6)=28.0d0*dx2**3*(2.0d0*dx2**2+dx1*(5.0d0*dx2+2.0d0*dx1))
      ci(7)=14.0d0*dx2**2*(5.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+dx1))
      ci(8)=8.0d0*dx2*(7.0d0*dx2**2+dx1*(7.0d0*dx2+dx1))
      ci(9)=28.0d0*dx2**2+dx1*(16.0d0*dx2+dx1)
      ci(10)=2.0d0*(4.0d0*dx2+dx1)
      ci(11)=1.0d0
    case(9)
      ci(1)=dx1**2*dx2**9
      ci(2)=dx1*dx2**8*(2.0d0*dx2+9.0d0*dx1)
      ci(3)=dx2**7*(dx2**2+18.0d0*dx1*(dx2+2.0d0*dx1))
      ci(4)=3.0d0*dx2**6*(3.0d0*dx2**2+4.0d0*dx1*(6.0d0*dx2+7.0d0*dx1))
      ci(5)=6.0d0*dx2**5*(6.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(6)=42.0d0*dx2**4*(2.0d0*dx2**2+3.0d0*dx1*(2.0d0*dx2+dx1))
      ci(7)=42.0d0*dx2**3*(3.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))
      ci(8)=6.0d0*dx2**2*(21.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+3.0d0*dx1))
      ci(9)=3.0d0*dx2*(28.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+dx1))
      ci(10)=36.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)
      ci(11)=9.0d0*dx2+2.0d0*dx1
      ci(12)=1.0d0
    case(10)
      ci(1)=dx1**2*dx2**10
      ci(2)=2.0d0*dx1*dx2**9*(dx2+5.0d0*dx1)
      ci(3)=dx2**8*(dx2**2+5.0d0*dx1*(4.0d0*dx2+9.0d0*dx1))
      ci(4)=10.0d0*dx2**7*(dx2**2+3.0d0*dx1*(3.0d0*dx2+4.0d0*dx1))
      ci(5)=15.0d0*dx2**6*(3.0d0*dx2**2+2.0d0*dx1*(8.0d0*dx2+7.0d0*dx1))
      ci(6)=12.0d0*dx2**5*(10.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))
      ci(7)=42.0d0*dx2**4*(5.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(8)=12.0d0*dx2**3*(21.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))
      ci(9)=15.0d0*dx2**2*(14.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))
      ci(10)=10.0d0*dx2*(12.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))
      ci(11)=45.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)
      ci(12)=2.0d0*(5.0d0*dx2+dx1)
      ci(13)=1.0d0
    case(11)
      ci(1)=dx1**2*dx2**11
      ci(2)=dx1*dx2**10*(2.0d0*dx2+11.0d0*dx1)
      ci(3)=dx2**9*(dx2**2+11.0d0*dx1*(2.0d0*dx2+5.0d0*dx1))
      ci(4)=11.0d0*dx2**8*(dx2**2+5.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(5)=55.0d0*dx2**7*(dx2**2+6.0d0*dx1*(dx2+dx1))
      ci(6)=33.0d0*dx2**6*(5.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+7.0d0*dx1))
      ci(7)=66.0d0*dx2**5*(5.0d0*dx2**2+7.0d0*dx1*(2.0d0*dx2+dx1))
      ci(8)=66.0d0*dx2**4*(7.0d0*dx2**2+dx1*(14.0d0*dx2+5.0d0*dx1))
      ci(9)=33.0d0*dx2**3*(14.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))
      ci(10)=55.0d0*dx2**2*(6.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(11)=11.0d0*dx2*(15.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))
      ci(12)=55.0d0*dx2**2+dx1*(22.0d0*dx2+dx1)
      ci(13)=11.0d0*dx2+2.0d0*dx1
      ci(14)=1.0d0
    case(12)
      ci(1)=dx1**2*dx2**12
      ci(2)=2.0d0*dx1*dx2**11*(dx2+6.0d0*dx1)
      ci(3)=dx2**10*(dx2**2+6.0d0*dx1*(4.0d0*dx2+11.0d0*dx1))
      ci(4)=4.0d0*dx2**9*(3.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+5.0d0*dx1))
      ci(5)=11.0d0*dx2**8*(6.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+9.0d0*dx1))
      ci(6)=22.0d0*dx2**7*(10.0d0*dx2**2+9.0d0*dx1*(5.0d0*dx2+4.0d0*dx1))
      ci(7)=33.0d0*dx2**6*(15.0d0*dx2**2+4.0d0*dx1*(12.0d0*dx2+7.0d0*dx1))
      ci(8)=264.0d0*dx2**5*(3.0d0*dx2**2+dx1*(7.0d0*dx2+3.0d0*dx1))
      ci(9)=33.0d0*dx2**4*(28.0d0*dx2**2+3.0d0*dx1*(16.0d0*dx2+5.0d0*dx1))
      ci(10)=22.0d0*dx2**3*(36.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+2.0d0*dx1))
      ci(11)=11.0d0*dx2**2*(45.0d0*dx2**2+2.0d0*dx1*(20.0d0*dx2+3.0d0*dx1))
      ci(12)=4.0d0*dx2*(55.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+dx1))
      ci(13)=66.0d0*dx2**2+dx1*(24.0d0*dx2+dx1)
      ci(14)=2.0d0*(6.0d0*dx2+dx1)
      ci(15)=1.0d0
    case(13)
      ci(1)=dx1**2*dx2**13
      ci(2)=dx1*dx2**12*(2.0d0*dx2+13.0d0*dx1)
      ci(3)=dx2**11*(dx2**2+26.0d0*dx1*(dx2+3.0d0*dx1))
      ci(4)=13.0d0*dx2**10*(dx2**2+2.0d0*dx1*(6.0d0*dx2+11.0d0*dx1))
      ci(5)=13.0d0*dx2**9*(6.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+5.0d0*dx1))
      ci(6)=143.0d0*dx2**8*(2.0d0*dx2**2+dx1*(10.0d0*dx2+9.0d0*dx1))
      ci(7)=143.0d0*dx2**7*(5.0d0*dx2**2+6.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(8)=429.0d0*dx2**6*(3.0d0*dx2**2+4.0d0*dx1*(2.0d0*dx2+dx1))
      ci(9)=429.0d0*dx2**5*(4.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(10)=143.0d0*dx2**4*(12.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1))
      ci(11)=143.0d0*dx2**3*(9.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1))
      ci(12)=13.0d0*dx2**2*(55.0d0*dx2**2+2.0d0*dx1*(22.0d0*dx2+3.0d0*dx1))
      ci(13)=13.0d0*dx2*(22.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))
      ci(14)=78.0d0*dx2**2+dx1*(26.0d0*dx2+dx1)
      ci(15)=13.0d0*dx2+2.0d0*dx1
      ci(16)=1.0d0
    case(14)
      ci(1)=dx1**2*dx2**14
      ci(2)=2.0d0*dx1*dx2**13*(dx2+7.0d0*dx1)
      ci(3)=dx2**12*(dx2**2+7.0d0*dx1*(4.0d0*dx2+13.0d0*dx1))
      ci(4)=14.0d0*dx2**11*(dx2**2+13.0d0*dx1*(dx2+2.0d0*dx1))
      ci(5)=91.0d0*dx2**10*(dx2**2+dx1*(8.0d0*dx2+11.0d0*dx1))
      ci(6)=182.0d0*dx2**9*(2.0d0*dx2**2+11.0d0*dx1*(dx2+dx1))
      ci(7)=1001.0d0*dx2**8*(dx2**2+dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(8)=286.0d0*dx2**7*(7.0d0*dx2**2+3.0d0*dx1*(7.0d0*dx2+4.0d0*dx1))
      ci(9)=429.0d0*dx2**6*(7.0d0*dx2**2+dx1*(16.0d0*dx2+7.0d0*dx1))
      ci(10)=286.0d0*dx2**5*(12.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+dx1))
      ci(11)=1001.0d0*dx2**4*(3.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))
      ci(12)=182.0d0*dx2**3*(11.0d0*dx2**2+dx1*(11.0d0*dx2+2.0d0*dx1))
      ci(13)=91.0d0*dx2**2*(11.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(14)=14.0d0*dx2*(26.0d0*dx2**2+dx1*(13.0d0*dx2+dx1))
      ci(15)=91.0d0*dx2**2+dx1*(28.0d0*dx2+dx1)
      ci(16)=2.0d0*(7.0d0*dx2+dx1)
      ci(17)=1.0d0
    case(15)
      ci(1)=dx1**2*dx2**15
      ci(2)=dx1*dx2**14*(2.0d0*dx2+15.0d0*dx1)
      ci(3)=dx2**13*(dx2**2+15.0d0*dx1*(2.0d0*dx2+7.0d0*dx1))
      ci(4)=5.0d0*dx2**12*(3.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+13.0d0*dx1))
      ci(5)=35.0d0*dx2**11*(3.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(6)=91.0d0*dx2**10*(5.0d0*dx2**2+3.0d0*dx1*(10.0d0*dx2+11.0d0*dx1))
      ci(7)=91.0d0*dx2**9*(15.0d0*dx2**2+11.0d0*dx1*(6.0d0*dx2+5.0d0*dx1))
      ci(8)=143.0d0*dx2**8*(21.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+9.0d0*dx1))
      ci(9)=715.0d0*dx2**7*(7.0d0*dx2**2+9.0d0*dx1*(2.0d0*dx2+dx1))
      ci(10)=715.0d0*dx2**6*(9.0d0*dx2**2+dx1*(18.0d0*dx2+7.0d0*dx1))
      ci(11)=143.0d0*dx2**5*(45.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+3.0d0*dx1))
      ci(12)=91.0d0*dx2**4*(55.0d0*dx2**2+3.0d0*dx1*(22.0d0*dx2+5.0d0*dx1))
      ci(13)=91.0d0*dx2**3*(33.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+dx1))
      ci(14)=35.0d0*dx2**2*(39.0d0*dx2**2+dx1*(26.0d0*dx2+3.0d0*dx1))
      ci(15)=5.0d0*dx2*(91.0d0*dx2**2+3.0d0*dx1*(14.0d0*dx2+dx1))
      ci(16)=105.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)
      ci(17)=15.0d0*dx2+2.0d0*dx1
      ci(18)=1.0d0
    case(16)
      ci(1)=dx1**2*dx2**16
      ci(2)=2.0d0*dx1*dx2**15*(dx2+8.0d0*dx1)
      ci(3)=dx2**14*(dx2**2+8.0d0*dx1*(4.0d0*dx2+15.0d0*dx1))
      ci(4)=16.0d0*dx2**13*(dx2**2+5.0d0*dx1*(3.0d0*dx2+7.0d0*dx1))
      ci(5)=20.0d0*dx2**12*(6.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+13.0d0*dx1))
      ci(6)=56.0d0*dx2**11*(10.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+6.0d0*dx1))
      ci(7)=364.0d0*dx2**10*(5.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+11.0d0*dx1))
      ci(8)=208.0d0*dx2**9*(21.0d0*dx2**2+11.0d0*dx1*(7.0d0*dx2+5.0d0*dx1))
      ci(9)=286.0d0*dx2**8*(28.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+9.0d0*dx1))
      ci(10)=2860.0d0*dx2**7*(4.0d0*dx2**2+dx1*(9.0d0*dx2+4.0d0*dx1))
      ci(11)=286.0d0*dx2**6*(45.0d0*dx2**2+4.0d0*dx1*(20.0d0*dx2+7.0d0*dx1))
      ci(12)=208.0d0*dx2**5*(55.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+3.0d0*dx1))
      ci(13)=364.0d0*dx2**4*(22.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1))
      ci(14)=56.0d0*dx2**3*(78.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+2.0d0*dx1))
      ci(15)=20.0d0*dx2**2*(91.0d0*dx2**2+2.0d0*dx1*(28.0d0*dx2+3.0d0*dx1))
      ci(16)=16.0d0*dx2*(35.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))
      ci(17)=120.0d0*dx2**2+dx1*(32.0d0*dx2+dx1)
      ci(18)=2.0d0*(8.0d0*dx2+dx1)
      ci(19)=1.0d0
    case(17)
      ci(1)=dx1**2*dx2**17
      ci(2)=dx1*dx2**16*(2.0d0*dx2+17.0d0*dx1)
      ci(3)=dx2**15*(dx2**2+34.0d0*dx1*(dx2+4.0d0*dx1))
      ci(4)=17.0d0*dx2**14*(dx2**2+8.0d0*dx1*(2.0d0*dx2+5.0d0*dx1))
      ci(5)=68.0d0*dx2**13*(2.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+7.0d0*dx1))
      ci(6)=68.0d0*dx2**12*(10.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+13.0d0*dx1))
      ci(7)=476.0d0*dx2**11*(5.0d0*dx2**2+26.0d0*dx1*(dx2+dx1))
      ci(8)=884.0d0*dx2**10*(7.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+11.0d0*dx1))
      ci(9)=442.0d0*dx2**9*(28.0d0*dx2**2+11.0d0*dx1*(8.0d0*dx2+5.0d0*dx1))
      ci(10)=4862.0d0*dx2**8*(4.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))
      ci(11)=4862.0d0*dx2**7*(5.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))
      ci(12)=442.0d0*dx2**6*(55.0d0*dx2**2+4.0d0*dx1*(22.0d0*dx2+7.0d0*dx1))
      ci(13)=884.0d0*dx2**5*(22.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+dx1))
      ci(14)=476.0d0*dx2**4*(26.0d0*dx2**2+dx1*(26.0d0*dx2+5.0d0*dx1))
      ci(15)=68.0d0*dx2**3*(91.0d0*dx2**2+10.0d0*dx1*(7.0d0*dx2+dx1))
      ci(16)=68.0d0*dx2**2*(35.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+dx1))
      ci(17)=17.0d0*dx2*(40.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))
      ci(18)=136.0d0*dx2**2+dx1*(34.0d0*dx2+dx1)
      ci(19)=17.0d0*dx2+2.0d0*dx1
      ci(20)=1.0d0
    case(18)
      ci(1)=dx1**2*dx2**18
      ci(2)=2.0d0*dx1*dx2**17*(dx2+9.0d0*dx1)
      ci(3)=dx2**16*(dx2**2+9.0d0*dx1*(4.0d0*dx2+17.0d0*dx1))
      ci(4)=6.0d0*dx2**15*(3.0d0*dx2**2+17.0d0*dx1*(3.0d0*dx2+8.0d0*dx1))
      ci(5)=51.0d0*dx2**14*(3.0d0*dx2**2+4.0d0*dx1*(8.0d0*dx2+15.0d0*dx1))
      ci(6)=408.0d0*dx2**13*(2.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+7.0d0*dx1))
      ci(7)=204.0d0*dx2**12*(15.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+13.0d0*dx1))
      ci(8)=408.0d0*dx2**11*(21.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+6.0d0*dx1))
      ci(9)=1326.0d0*dx2**10*(14.0d0*dx2**2+3.0d0*dx1*(16.0d0*dx2+11.0d0*dx1))
      ci(10)=884.0d0*dx2**9*(36.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+5.0d0*dx1))
      ci(11)=4862.0d0*dx2**8*(9.0d0*dx2**2+dx1*(20.0d0*dx2+9.0d0*dx1))
      ci(12)=884.0d0*dx2**7*(55.0d0*dx2**2+9.0d0*dx1*(11.0d0*dx2+4.0d0*dx1))
      ci(13)=1326.0d0*dx2**6*(33.0d0*dx2**2+2.0d0*dx1*(24.0d0*dx2+7.0d0*dx1))
      ci(14)=408.0d0*dx2**5*(78.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+3.0d0*dx1))
      ci(15)=204.0d0*dx2**4*(91.0d0*dx2**2+3.0d0*dx1*(28.0d0*dx2+5.0d0*dx1))
      ci(16)=408.0d0*dx2**3*(21.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1))
      ci(17)=51.0d0*dx2**2*(60.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1))
      ci(18)=6.0d0*dx2*(136.0d0*dx2**2+3.0d0*dx1*(17.0d0*dx2+dx1))
      ci(19)=153.0d0*dx2**2+dx1*(36.0d0*dx2+dx1)
      ci(20)=2.0d0*(9.0d0*dx2+dx1)
      ci(21)=1.0d0
    case(19)
      ci(1)=dx1**2*dx2**19
      ci(2)=dx1*dx2**18*(2.0d0*dx2+19.0d0*dx1)
      ci(3)=dx2**17*(dx2**2+19.0d0*dx1*(2.0d0*dx2+9.0d0*dx1))
      ci(4)=19.0d0*dx2**16*(dx2**2+3.0d0*dx1*(6.0d0*dx2+17.0d0*dx1))
      ci(5)=57.0d0*dx2**15*(3.0d0*dx2**2+34.0d0*dx1*(dx2+2.0d0*dx1))
      ci(6)=969.0d0*dx2**14*(dx2**2+4.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(7)=3876.0d0*dx2**13*(dx2**2+dx1*(6.0d0*dx2+7.0d0*dx1))
      ci(8)=3876.0d0*dx2**12*(3.0d0*dx2**2+dx1*(14.0d0*dx2+13.0d0*dx1))
      ci(9)=1938.0d0*dx2**11*(14.0d0*dx2**2+13.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(10)=8398.0d0*dx2**10*(6.0d0*dx2**2+dx1*(18.0d0*dx2+11.0d0*dx1))
      ci(11)=8398.0d0*dx2**9*(9.0d0*dx2**2+11.0d0*dx1*(2.0d0*dx2+dx1))
      ci(12)=8398.0d0*dx2**8*(11.0d0*dx2**2+dx1*(22.0d0*dx2+9.0d0*dx1))
      ci(13)=8398.0d0*dx2**7*(11.0d0*dx2**2+6.0d0*dx1*(3.0d0*dx2+dx1))
      ci(14)=1938.0d0*dx2**6*(39.0d0*dx2**2+2.0d0*dx1*(26.0d0*dx2+7.0d0*dx1))
      ci(15)=3876.0d0*dx2**5*(13.0d0*dx2**2+dx1*(14.0d0*dx2+3.0d0*dx1))
      ci(16)=3876.0d0*dx2**4*(7.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(17)=969.0d0*dx2**3*(12.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(18)=57.0d0*dx2**2*(68.0d0*dx2**2+dx1*(34.0d0*dx2+3.0d0*dx1))
      ci(19)=19.0d0*dx2*(51.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))
      ci(20)=171.0d0*dx2**2+dx1*(38.0d0*dx2+dx1)
      ci(21)=19.0d0*dx2+2.0d0*dx1
      ci(22)=1.0d0
    case(20)
      ci(1)=dx1**2*dx2**20
      ci(2)=2.0d0*dx1*dx2**19*(dx2+10.0d0*dx1)
      ci(3)=dx2**18*(dx2**2+10.0d0*dx1*(4.0d0*dx2+19.0d0*dx1))
      ci(4)=20.0d0*dx2**17*(dx2**2+19.0d0*dx1*(dx2+3.0d0*dx1))
      ci(5)=95.0d0*dx2**16*(2.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+17.0d0*dx1))
      ci(6)=114.0d0*dx2**15*(10.0d0*dx2**2+17.0d0*dx1*(5.0d0*dx2+8.0d0*dx1))
      ci(7)=969.0d0*dx2**14*(5.0d0*dx2**2+8.0d0*dx1*(4.0d0*dx2+5.0d0*dx1))
      ci(8)=15504.0d0*dx2**13*(dx2**2+5.0d0*dx1*(dx2+dx1))
      ci(9)=9690.0d0*dx2**12*(4.0d0*dx2**2+dx1*(16.0d0*dx2+13.0d0*dx1))
      ci(10)=6460.0d0*dx2**11*(12.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(11)=8398.0d0*dx2**10*(15.0d0*dx2**2+2.0d0*dx1*(20.0d0*dx2+11.0d0*dx1))
      ci(12)=33592.0d0*dx2**9*(5.0d0*dx2**2+dx1*(11.0d0*dx2+5.0d0*dx1))
      ci(13)=8398.0d0*dx2**8*(22.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(14)=6460.0d0*dx2**7*(26.0d0*dx2**2+3.0d0*dx1*(13.0d0*dx2+4.0d0*dx1))
      ci(15)=9690.0d0*dx2**6*(13.0d0*dx2**2+4.0d0*dx1*(4.0d0*dx2+dx1))
      ci(16)=15504.0d0*dx2**5*(5.0d0*dx2**2+dx1*(5.0d0*dx2+dx1))
      ci(17)=969.0d0*dx2**4*(40.0d0*dx2**2+dx1*(32.0d0*dx2+5.0d0*dx1))
      ci(18)=114.0d0*dx2**3*(136.0d0*dx2**2+5.0d0*dx1*(17.0d0*dx2+2.0d0*dx1))
      ci(19)=95.0d0*dx2**2*(51.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+dx1))
      ci(20)=20.0d0*dx2*(57.0d0*dx2**2+dx1*(19.0d0*dx2+dx1))
      ci(21)=190.0d0*dx2**2+dx1*(40.0d0*dx2+dx1)
      ci(22)=2.0d0*(10.0d0*dx2+dx1)
      ci(23)=1.0d0
    case(21)
      ci(1)=dx1**2*dx2**21
      ci(2)=dx1*dx2**20*(2.0d0*dx2+21.0d0*dx1)
      ci(3)=dx2**19*(dx2**2+42.0d0*dx1*(dx2+5.0d0*dx1))
      ci(4)=7.0d0*dx2**18*(3.0d0*dx2**2+10.0d0*dx1*(6.0d0*dx2+19.0d0*dx1))
      ci(5)=35.0d0*dx2**17*(6.0d0*dx2**2+19.0d0*dx1*(4.0d0*dx2+9.0d0*dx1))
      ci(6)=133.0d0*dx2**16*(10.0d0*dx2**2+9.0d0*dx1*(10.0d0*dx2+17.0d0*dx1))
      ci(7)=399.0d0*dx2**15*(15.0d0*dx2**2+34.0d0*dx1*(3.0d0*dx2+4.0d0*dx1))
      ci(8)=969.0d0*dx2**14*(21.0d0*dx2**2+8.0d0*dx1*(14.0d0*dx2+15.0d0*dx1))
      ci(9)=1938.0d0*dx2**13*(28.0d0*dx2**2+15.0d0*dx1*(8.0d0*dx2+7.0d0*dx1))
      ci(10)=3230.0d0*dx2**12*(36.0d0*dx2**2+7.0d0*dx1*(18.0d0*dx2+13.0d0*dx1))
      ci(11)=4522.0d0*dx2**11*(45.0d0*dx2**2+26.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))
      ci(12)=58786.0d0*dx2**10*(5.0d0*dx2**2+6.0d0*dx1*(2.0d0*dx2+dx1))
      ci(13)=58786.0d0*dx2**9*(6.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(14)=4522.0d0*dx2**8*(78.0d0*dx2**2+5.0d0*dx1*(26.0d0*dx2+9.0d0*dx1))
      ci(15)=3230.0d0*dx2**7*(91.0d0*dx2**2+18.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))
      ci(16)=1938.0d0*dx2**6*(105.0d0*dx2**2+4.0d0*dx1*(30.0d0*dx2+7.0d0*dx1))
      ci(17)=969.0d0*dx2**5*(120.0d0*dx2**2+7.0d0*dx1*(16.0d0*dx2+3.0d0*dx1))
      ci(18)=399.0d0*dx2**4*(136.0d0*dx2**2+3.0d0*dx1*(34.0d0*dx2+5.0d0*dx1))
      ci(19)=133.0d0*dx2**3*(153.0d0*dx2**2+10.0d0*dx1*(9.0d0*dx2+dx1))
      ci(20)=35.0d0*dx2**2*(171.0d0*dx2**2+2.0d0*dx1*(38.0d0*dx2+3.0d0*dx1))
      ci(21)=7.0d0*dx2*(190.0d0*dx2**2+3.0d0*dx1*(20.0d0*dx2+dx1))
      ci(22)=210.0d0*dx2**2+dx1*(42.0d0*dx2+dx1)
      ci(23)=21.0d0*dx2+2.0d0*dx1
      ci(24)=1.0d0
    case(22)
      ci(1)=dx1**2*dx2**22
      ci(2)=2.0d0*dx1*dx2**21*(dx2+11.0d0*dx1)
      ci(3)=dx2**20*(dx2**2+11.0d0*dx1*(4.0d0*dx2+21.0d0*dx1))
      ci(4)=22.0d0*dx2**19*(dx2**2+7.0d0*dx1*(3.0d0*dx2+10.0d0*dx1))
      ci(5)=77.0d0*dx2**18*(3.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+19.0d0*dx1))
      ci(6)=154.0d0*dx2**17*(10.0d0*dx2**2+19.0d0*dx1*(5.0d0*dx2+9.0d0*dx1))
      ci(7)=1463.0d0*dx2**16*(5.0d0*dx2**2+3.0d0*dx1*(12.0d0*dx2+17.0d0*dx1))
      ci(8)=1254.0d0*dx2**15*(21.0d0*dx2**2+17.0d0*dx1*(7.0d0*dx2+8.0d0*dx1))
      ci(9)=10659.0d0*dx2**14*(7.0d0*dx2**2+2.0d0*dx1*(16.0d0*dx2+15.0d0*dx1))
      ci(10)=14212.0d0*dx2**13*(12.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+7.0d0*dx1))
      ci(11)=7106.0d0*dx2**12*(45.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+13.0d0*dx1))
      ci(12)=9044.0d0*dx2**11*(55.0d0*dx2**2+13.0d0*dx1*(11.0d0*dx2+6.0d0*dx1))
      ci(13)=58786.0d0*dx2**10*(11.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1))
      ci(14)=9044.0d0*dx2**9*(78.0d0*dx2**2+11.0d0*dx1*(13.0d0*dx2+5.0d0*dx1))
      ci(15)=7106.0d0*dx2**8*(91.0d0*dx2**2+5.0d0*dx1*(28.0d0*dx2+9.0d0*dx1))
      ci(16)=14212.0d0*dx2**7*(35.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+4.0d0*dx1))
      ci(17)=10659.0d0*dx2**6*(30.0d0*dx2**2+dx1*(32.0d0*dx2+7.0d0*dx1))
      ci(18)=1254.0d0*dx2**5*(136.0d0*dx2**2+7.0d0*dx1*(17.0d0*dx2+3.0d0*dx1))
      ci(19)=1463.0d0*dx2**4*(51.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))
      ci(20)=154.0d0*dx2**3*(171.0d0*dx2**2+5.0d0*dx1*(19.0d0*dx2+2.0d0*dx1))
      ci(21)=77.0d0*dx2**2*(95.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))
      ci(22)=22.0d0*dx2*(70.0d0*dx2**2+dx1*(21.0d0*dx2+dx1))
      ci(23)=231.0d0*dx2**2+dx1*(44.0d0*dx2+dx1)
      ci(24)=2.0d0*(11.0d0*dx2+dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>22, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^3 * (x-x2)^n2 as sum_k=0^(3+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_3(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**3
      ci(2)=3.0d0*dx1**2
      ci(3)=3.0d0*dx1
      ci(4)=1.0d0
    case(1)
      ci(1)=dx1**3*dx2
      ci(2)=dx1**2*(3.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1*(dx2+dx1)
      ci(4)=dx2+3.0d0*dx1
      ci(5)=1.0d0
    case(2)
      ci(1)=dx1**3*dx2**2
      ci(2)=dx1**2*dx2*(3.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1*(3.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(4)=dx2**2+3.0d0*dx1*(2.0d0*dx2+dx1)
      ci(5)=2.0d0*dx2+3.0d0*dx1
      ci(6)=1.0d0
    case(3)
      ci(1)=dx1**3*dx2**3
      ci(2)=3.0d0*dx1**2*dx2**2*(dx2+dx1)
      ci(3)=3.0d0*dx1*dx2*(dx2**2+dx1*(3.0d0*dx2+dx1))
      ci(4)=dx2**3+dx1*(9.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))
      ci(5)=3.0d0*(dx2**2+dx1*(3.0d0*dx2+dx1))
      ci(6)=3.0d0*(dx2+dx1)
      ci(7)=1.0d0
    case(4)
      ci(1)=dx1**3*dx2**4
      ci(2)=dx1**2*dx2**3*(3.0d0*dx2+4.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**2*(dx2**2+2.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=dx2*(dx2**3+2.0d0*dx1*(6.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(5)=4.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))
      ci(6)=3.0d0*(2.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))
      ci(7)=4.0d0*dx2+3.0d0*dx1
      ci(8)=1.0d0
    case(5)
      ci(1)=dx1**3*dx2**5
      ci(2)=dx1**2*dx2**4*(3.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1*dx2**3*(3.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(4)=dx2**2*(dx2**3+5.0d0*dx1*(3.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx2*(dx2**3+dx1*(6.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))
      ci(6)=10.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))
      ci(7)=10.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+dx1)
      ci(8)=5.0d0*dx2+3.0d0*dx1
      ci(9)=1.0d0
    case(6)
      ci(1)=dx1**3*dx2**6
      ci(2)=3.0d0*dx1**2*dx2**5*(dx2+2.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**4*(dx2**2+dx1*(6.0d0*dx2+5.0d0*dx1))
      ci(4)=dx2**3*(dx2**3+dx1*(18.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(5)=3.0d0*dx2**2*(2.0d0*dx2**3+5.0d0*dx1*(3.0d0*dx2**2+dx1*(4.0d0*dx2+dx1)))
      ci(6)=3.0d0*dx2*(5.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(7)=20.0d0*dx2**3+dx1*(45.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))
      ci(8)=3.0d0*(5.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(9)=3.0d0*(2.0d0*dx2+dx1)
      ci(10)=1.0d0
    case(7)
      ci(1)=dx1**3*dx2**7
      ci(2)=dx1**2*dx2**6*(3.0d0*dx2+7.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**5*(dx2**2+7.0d0*dx1*(dx2+dx1))
      ci(4)=dx2**4*(dx2**3+7.0d0*dx1*(3.0d0*dx2**2+dx1*(9.0d0*dx2+5.0d0*dx1)))
      ci(5)=7.0d0*dx2**3*(dx2**3+dx1*(9.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(6)=21.0d0*dx2**2*(dx2**3+dx1*(5.0d0*dx2**2+dx1*(5.0d0*dx2+dx1)))
      ci(7)=7.0d0*dx2*(5.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(8)=35.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(21.0d0*dx2+dx1))
      ci(9)=3.0d0*(7.0d0*dx2**2+dx1*(7.0d0*dx2+dx1))
      ci(10)=7.0d0*dx2+3.0d0*dx1
      ci(11)=1.0d0
    case(8)
      ci(1)=dx1**3*dx2**8
      ci(2)=dx1**2*dx2**7*(3.0d0*dx2+8.0d0*dx1)
      ci(3)=dx1*dx2**6*(3.0d0*dx2**2+4.0d0*dx1*(6.0d0*dx2+7.0d0*dx1))
      ci(4)=dx2**5*(dx2**3+4.0d0*dx1*(6.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(5)=2.0d0*dx2**4*(4.0d0*dx2**3+7.0d0*dx1*(6.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1)))
      ci(6)=14.0d0*dx2**3*(2.0d0*dx2**3+dx1*(12.0d0*dx2**2+dx1*(15.0d0*dx2+4.0d0*dx1)))
      ci(7)=14.0d0*dx2**2*(4.0d0*dx2**3+dx1*(15.0d0*dx2**2+2.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(8)=2.0d0*dx2*(35.0d0*dx2**3+2.0d0*dx1*(42.0d0*dx2**2+dx1*(21.0d0*dx2+2.0d0*dx1)))
      ci(9)=56.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))
      ci(10)=28.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+dx1)
      ci(11)=8.0d0*dx2+3.0d0*dx1
      ci(12)=1.0d0
    case(9)
      ci(1)=dx1**3*dx2**9
      ci(2)=3.0d0*dx1**2*dx2**8*(dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**7*(dx2**2+3.0d0*dx1*(3.0d0*dx2+4.0d0*dx1))
      ci(4)=dx2**6*(dx2**3+3.0d0*dx1*(9.0d0*dx2**2+4.0d0*dx1*(9.0d0*dx2+7.0d0*dx1)))
      ci(5)=9.0d0*dx2**5*(dx2**3+2.0d0*dx1*(6.0d0*dx2**2+7.0d0*dx1*(2.0d0*dx2+dx1)))
      ci(6)=18.0d0*dx2**4*(2.0d0*dx2**3+7.0d0*dx1*(2.0d0*dx2**2+dx1*(3.0d0*dx2+dx1)))
      ci(7)=42.0d0*dx2**3*(2.0d0*dx2**3+dx1*(9.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(8)=18.0d0*dx2**2*(7.0d0*dx2**3+dx1*(21.0d0*dx2**2+2.0d0*dx1*(7.0d0*dx2+dx1)))
      ci(9)=9.0d0*dx2*(14.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(10)=84.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1*(27.0d0*dx2+dx1))
      ci(11)=3.0d0*(12.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))
      ci(12)=3.0d0*(3.0d0*dx2+dx1)
      ci(13)=1.0d0
    case(10)
      ci(1)=dx1**3*dx2**10
      ci(2)=dx1**2*dx2**9*(3.0d0*dx2+10.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**8*(dx2**2+5.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(4)=dx2**7*(dx2**3+15.0d0*dx1*(2.0d0*dx2**2+dx1*(9.0d0*dx2+8.0d0*dx1)))
      ci(5)=5.0d0*dx2**6*(2.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)))
      ci(6)=9.0d0*dx2**5*(5.0d0*dx2**3+2.0d0*dx1*(20.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)))
      ci(7)=6.0d0*dx2**4*(20.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1)))
      ci(8)=6.0d0*dx2**3*(35.0d0*dx2**3+dx1*(126.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+4.0d0*dx1)))
      ci(9)=9.0d0*dx2**2*(28.0d0*dx2**3+5.0d0*dx1*(14.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))
      ci(10)=5.0d0*dx2*(42.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(27.0d0*dx2+2.0d0*dx1)))
      ci(11)=120.0d0*dx2**3+dx1*(135.0d0*dx2**2+dx1*(30.0d0*dx2+dx1))
      ci(12)=3.0d0*(15.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))
      ci(13)=10.0d0*dx2+3.0d0*dx1
      ci(14)=1.0d0
    case(11)
      ci(1)=dx1**3*dx2**11
      ci(2)=dx1**2*dx2**10*(3.0d0*dx2+11.0d0*dx1)
      ci(3)=dx1*dx2**9*(3.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+5.0d0*dx1))
      ci(4)=dx2**8*(dx2**3+33.0d0*dx1*(dx2**2+5.0d0*dx1*(dx2+dx1)))
      ci(5)=11.0d0*dx2**7*(dx2**3+15.0d0*dx1*(dx2**2+dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(6)=11.0d0*dx2**6*(5.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+7.0d0*dx1)))
      ci(7)=33.0d0*dx2**5*(5.0d0*dx2**3+2.0d0*dx1*(15.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(8)=66.0d0*dx2**4*(5.0d0*dx2**3+dx1*(21.0d0*dx2**2+dx1*(21.0d0*dx2+5.0d0*dx1)))
      ci(9)=33.0d0*dx2**3*(14.0d0*dx2**3+dx1*(42.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(10)=11.0d0*dx2**2*(42.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(11)=11.0d0*dx2*(30.0d0*dx2**3+dx1*(45.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))
      ci(12)=165.0d0*dx2**3+dx1*(165.0d0*dx2**2+dx1*(33.0d0*dx2+dx1))
      ci(13)=55.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+dx1)
      ci(14)=11.0d0*dx2+3.0d0*dx1
      ci(15)=1.0d0
    case(12)
      ci(1)=dx1**3*dx2**12
      ci(2)=3.0d0*dx1**2*dx2**11*(dx2+4.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**10*(dx2**2+2.0d0*dx1*(6.0d0*dx2+11.0d0*dx1))
      ci(4)=dx2**9*(dx2**3+2.0d0*dx1*(18.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+10.0d0*dx1)))
      ci(5)=3.0d0*dx2**8*(4.0d0*dx2**3+11.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+3.0d0*dx1)))
      ci(6)=33.0d0*dx2**7*(2.0d0*dx2**3+dx1*(20.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+8.0d0*dx1)))
      ci(7)=11.0d0*dx2**6*(20.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+4.0d0*dx1*(18.0d0*dx2+7.0d0*dx1)))
      ci(8)=99.0d0*dx2**5*(5.0d0*dx2**3+4.0d0*dx1*(6.0d0*dx2**2+dx1*(7.0d0*dx2+2.0d0*dx1)))
      ci(9)=99.0d0*dx2**4*(8.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1)))
      ci(10)=11.0d0*dx2**3*(84.0d0*dx2**3+dx1*(216.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+4.0d0*dx1)))
      ci(11)=33.0d0*dx2**2*(24.0d0*dx2**3+dx1*(45.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+dx1)))
      ci(12)=3.0d0*dx2*(165.0d0*dx2**3+2.0d0*dx1*(110.0d0*dx2**2+dx1*(33.0d0*dx2+2.0d0*dx1)))
      ci(13)=220.0d0*dx2**3+dx1*(198.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))
      ci(14)=3.0d0*(22.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))
      ci(15)=3.0d0*(4.0d0*dx2+dx1)
      ci(16)=1.0d0
    case(13)
      ci(1)=dx1**3*dx2**13
      ci(2)=dx1**2*dx2**12*(3.0d0*dx2+13.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**11*(dx2**2+13.0d0*dx1*(dx2+2.0d0*dx1))
      ci(4)=dx2**10*(dx2**3+13.0d0*dx1*(3.0d0*dx2**2+2.0d0*dx1*(9.0d0*dx2+11.0d0*dx1)))
      ci(5)=13.0d0*dx2**9*(dx2**3+dx1*(18.0d0*dx2**2+11.0d0*dx1*(6.0d0*dx2+5.0d0*dx1)))
      ci(6)=39.0d0*dx2**8*(2.0d0*dx2**3+11.0d0*dx1*(2.0d0*dx2**2+dx1*(5.0d0*dx2+3.0d0*dx1)))
      ci(7)=143.0d0*dx2**7*(2.0d0*dx2**3+3.0d0*dx1*(5.0d0*dx2**2+dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(8)=143.0d0*dx2**6*(5.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+4.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(9)=1287.0d0*dx2**5*(dx2**3+dx1*(4.0d0*dx2**2+dx1*(4.0d0*dx2+dx1)))
      ci(10)=143.0d0*dx2**4*(12.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(27.0d0*dx2+5.0d0*dx1)))
      ci(11)=143.0d0*dx2**3*(12.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(12)=39.0d0*dx2**2*(33.0d0*dx2**3+dx1*(55.0d0*dx2**2+2.0d0*dx1*(11.0d0*dx2+dx1)))
      ci(13)=13.0d0*dx2*(55.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))
      ci(14)=286.0d0*dx2**3+dx1*(234.0d0*dx2**2+dx1*(39.0d0*dx2+dx1))
      ci(15)=3.0d0*(26.0d0*dx2**2+dx1*(13.0d0*dx2+dx1))
      ci(16)=13.0d0*dx2+3.0d0*dx1
      ci(17)=1.0d0
    case(14)
      ci(1)=dx1**3*dx2**14
      ci(2)=dx1**2*dx2**13*(3.0d0*dx2+14.0d0*dx1)
      ci(3)=dx1*dx2**12*(3.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+13.0d0*dx1))
      ci(4)=dx2**11*(dx2**3+7.0d0*dx1*(6.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+4.0d0*dx1)))
      ci(5)=7.0d0*dx2**10*(2.0d0*dx2**3+13.0d0*dx1*(3.0d0*dx2**2+dx1*(12.0d0*dx2+11.0d0*dx1)))
      ci(6)=91.0d0*dx2**9*(dx2**3+dx1*(12.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(7)=91.0d0*dx2**8*(4.0d0*dx2**3+33.0d0*dx1*(dx2**2+dx1*(2.0d0*dx2+dx1)))
      ci(8)=143.0d0*dx2**7*(7.0d0*dx2**3+3.0d0*dx1*(14.0d0*dx2**2+dx1*(21.0d0*dx2+8.0d0*dx1)))
      ci(9)=143.0d0*dx2**6*(14.0d0*dx2**3+3.0d0*dx1*(21.0d0*dx2**2+dx1*(24.0d0*dx2+7.0d0*dx1)))
      ci(10)=143.0d0*dx2**5*(21.0d0*dx2**3+dx1*(72.0d0*dx2**2+7.0d0*dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(11)=143.0d0*dx2**4*(24.0d0*dx2**3+7.0d0*dx1*(9.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))
      ci(12)=91.0d0*dx2**3*(33.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(33.0d0*dx2+4.0d0*dx1)))
      ci(13)=91.0d0*dx2**2*(22.0d0*dx2**3+dx1*(33.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(14)=7.0d0*dx2*(143.0d0*dx2**3+dx1*(156.0d0*dx2**2+dx1*(39.0d0*dx2+2.0d0*dx1)))
      ci(15)=364.0d0*dx2**3+dx1*(273.0d0*dx2**2+dx1*(42.0d0*dx2+dx1))
      ci(16)=91.0d0*dx2**2+3.0d0*dx1*(14.0d0*dx2+dx1)
      ci(17)=14.0d0*dx2+3.0d0*dx1
      ci(18)=1.0d0
    case(15)
      ci(1)=dx1**3*dx2**15
      ci(2)=3.0d0*dx1**2*dx2**14*(dx2+5.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**13*(dx2**2+5.0d0*dx1*(3.0d0*dx2+7.0d0*dx1))
      ci(4)=dx2**12*(dx2**3+5.0d0*dx1*(9.0d0*dx2**2+7.0d0*dx1*(9.0d0*dx2+13.0d0*dx1)))
      ci(5)=15.0d0*dx2**11*(dx2**3+7.0d0*dx1*(3.0d0*dx2**2+13.0d0*dx1*(dx2+dx1)))
      ci(6)=21.0d0*dx2**10*(5.0d0*dx2**3+13.0d0*dx1*(5.0d0*dx2**2+dx1*(15.0d0*dx2+11.0d0*dx1)))
      ci(7)=91.0d0*dx2**9*(5.0d0*dx2**3+dx1*(45.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+5.0d0*dx1)))
      ci(8)=39.0d0*dx2**8*(35.0d0*dx2**3+11.0d0*dx1*(21.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+3.0d0*dx1)))
      ci(9)=429.0d0*dx2**7*(7.0d0*dx2**3+5.0d0*dx1*(7.0d0*dx2**2+3.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(10)=715.0d0*dx2**6*(7.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(27.0d0*dx2+7.0d0*dx1)))
      ci(11)=429.0d0*dx2**5*(15.0d0*dx2**3+dx1*(45.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+dx1)))
      ci(12)=39.0d0*dx2**4*(165.0d0*dx2**3+7.0d0*dx1*(55.0d0*dx2**2+dx1*(33.0d0*dx2+5.0d0*dx1)))
      ci(13)=91.0d0*dx2**3*(55.0d0*dx2**3+dx1*(99.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+dx1)))
      ci(14)=21.0d0*dx2**2*(143.0d0*dx2**3+5.0d0*dx1*(39.0d0*dx2**2+dx1*(13.0d0*dx2+dx1)))
      ci(15)=15.0d0*dx2*(91.0d0*dx2**3+dx1*(91.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)))
      ci(16)=455.0d0*dx2**3+dx1*(315.0d0*dx2**2+dx1*(45.0d0*dx2+dx1))
      ci(17)=3.0d0*(35.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))
      ci(18)=3.0d0*(5.0d0*dx2+dx1)
      ci(19)=1.0d0
    case(16)
      ci(1)=dx1**3*dx2**16
      ci(2)=dx1**2*dx2**15*(3.0d0*dx2+16.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**14*(dx2**2+8.0d0*dx1*(2.0d0*dx2+5.0d0*dx1))
      ci(4)=dx2**13*(dx2**3+8.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+14.0d0*dx1)))
      ci(5)=4.0d0*dx2**12*(4.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+13.0d0*dx1)))
      ci(6)=12.0d0*dx2**11*(10.0d0*dx2**3+7.0d0*dx1*(20.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+4.0d0*dx1)))
      ci(7)=28.0d0*dx2**10*(20.0d0*dx2**3+13.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(18.0d0*dx2+11.0d0*dx1)))
      ci(8)=52.0d0*dx2**9*(35.0d0*dx2**3+2.0d0*dx1*(126.0d0*dx2**2+11.0d0*dx1*(21.0d0*dx2+10.0d0*dx1)))
      ci(9)=78.0d0*dx2**8*(56.0d0*dx2**3+11.0d0*dx1*(28.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+3.0d0*dx1)))
      ci(10)=286.0d0*dx2**7*(28.0d0*dx2**3+5.0d0*dx1*(24.0d0*dx2**2+dx1*(27.0d0*dx2+8.0d0*dx1)))
      ci(11)=286.0d0*dx2**6*(40.0d0*dx2**3+dx1*(135.0d0*dx2**2+4.0d0*dx1*(30.0d0*dx2+7.0d0*dx1)))
      ci(12)=78.0d0*dx2**5*(165.0d0*dx2**3+4.0d0*dx1*(110.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+2.0d0*dx1)))
      ci(13)=52.0d0*dx2**4*(220.0d0*dx2**3+7.0d0*dx1*(66.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1)))
      ci(14)=28.0d0*dx2**3*(286.0d0*dx2**3+dx1*(468.0d0*dx2**2+5.0d0*dx1*(39.0d0*dx2+4.0d0*dx1)))
      ci(15)=12.0d0*dx2**2*(364.0d0*dx2**3+5.0d0*dx1*(91.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+dx1)))
      ci(16)=4.0d0*dx2*(455.0d0*dx2**3+2.0d0*dx1*(210.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))
      ci(17)=560.0d0*dx2**3+dx1*(360.0d0*dx2**2+dx1*(48.0d0*dx2+dx1))
      ci(18)=3.0d0*(40.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))
      ci(19)=16.0d0*dx2+3.0d0*dx1
      ci(20)=1.0d0
    case(17)
      ci(1)=dx1**3*dx2**17
      ci(2)=dx1**2*dx2**16*(3.0d0*dx2+17.0d0*dx1)
      ci(3)=dx1*dx2**15*(3.0d0*dx2**2+17.0d0*dx1*(3.0d0*dx2+8.0d0*dx1))
      ci(4)=dx2**14*(dx2**3+17.0d0*dx1*(3.0d0*dx2**2+8.0d0*dx1*(3.0d0*dx2+5.0d0*dx1)))
      ci(5)=17.0d0*dx2**13*(dx2**3+4.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+7.0d0*dx1)))
      ci(6)=68.0d0*dx2**12*(2.0d0*dx2**3+dx1*(30.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+13.0d0*dx1)))
      ci(7)=68.0d0*dx2**11*(10.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(8)=68.0d0*dx2**10*(35.0d0*dx2**3+13.0d0*dx1*(21.0d0*dx2**2+2.0d0*dx1*(21.0d0*dx2+11.0d0*dx1)))
      ci(9)=442.0d0*dx2**9*(14.0d0*dx2**3+dx1*(84.0d0*dx2**2+11.0d0*dx1*(12.0d0*dx2+5.0d0*dx1)))
      ci(10)=442.0d0*dx2**8*(28.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(11)=4862.0d0*dx2**7*(4.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(15.0d0*dx2+4.0d0*dx1)))
      ci(12)=442.0d0*dx2**6*(55.0d0*dx2**3+dx1*(165.0d0*dx2**2+4.0d0*dx1*(33.0d0*dx2+7.0d0*dx1)))
      ci(13)=442.0d0*dx2**5*(55.0d0*dx2**3+2.0d0*dx1*(66.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(14)=68.0d0*dx2**4*(286.0d0*dx2**3+7.0d0*dx1*(78.0d0*dx2**2+dx1*(39.0d0*dx2+5.0d0*dx1)))
      ci(15)=68.0d0*dx2**3*(182.0d0*dx2**3+dx1*(273.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+2.0d0*dx1)))
      ci(16)=68.0d0*dx2**2*(91.0d0*dx2**3+dx1*(105.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+dx1)))
      ci(17)=17.0d0*dx2*(140.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(24.0d0*dx2+dx1)))
      ci(18)=680.0d0*dx2**3+dx1*(408.0d0*dx2**2+dx1*(51.0d0*dx2+dx1))
      ci(19)=136.0d0*dx2**2+3.0d0*dx1*(17.0d0*dx2+dx1)
      ci(20)=17.0d0*dx2+3.0d0*dx1
      ci(21)=1.0d0
    case(18)
      ci(1)=dx1**3*dx2**18
      ci(2)=3.0d0*dx1**2*dx2**17*(dx2+6.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**16*(dx2**2+3.0d0*dx1*(6.0d0*dx2+17.0d0*dx1))
      ci(4)=dx2**15*(dx2**3+3.0d0*dx1*(18.0d0*dx2**2+17.0d0*dx1*(9.0d0*dx2+16.0d0*dx1)))
      ci(5)=9.0d0*dx2**14*(2.0d0*dx2**3+17.0d0*dx1*(3.0d0*dx2**2+4.0d0*dx1*(4.0d0*dx2+5.0d0*dx1)))
      ci(6)=153.0d0*dx2**13*(dx2**3+4.0d0*dx1*(4.0d0*dx2**2+dx1*(15.0d0*dx2+14.0d0*dx1)))
      ci(7)=204.0d0*dx2**12*(4.0d0*dx2**3+dx1*(45.0d0*dx2**2+7.0d0*dx1*(18.0d0*dx2+13.0d0*dx1)))
      ci(8)=612.0d0*dx2**11*(5.0d0*dx2**3+dx1*(42.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+4.0d0*dx1)))
      ci(9)=306.0d0*dx2**10*(28.0d0*dx2**3+13.0d0*dx1*(14.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1)))
      ci(10)=442.0d0*dx2**9*(42.0d0*dx2**3+dx1*(216.0d0*dx2**2+11.0d0*dx1*(27.0d0*dx2+10.0d0*dx1)))
      ci(11)=1326.0d0*dx2**8*(24.0d0*dx2**3+11.0d0*dx1*(9.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1)))
      ci(12)=1326.0d0*dx2**7*(33.0d0*dx2**3+dx1*(110.0d0*dx2**2+3.0d0*dx1*(33.0d0*dx2+8.0d0*dx1)))
      ci(13)=442.0d0*dx2**6*(110.0d0*dx2**3+3.0d0*dx1*(99.0d0*dx2**2+2.0d0*dx1*(36.0d0*dx2+7.0d0*dx1)))
      ci(14)=306.0d0*dx2**5*(143.0d0*dx2**3+2.0d0*dx1*(156.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+2.0d0*dx1)))
      ci(15)=612.0d0*dx2**4*(52.0d0*dx2**3+dx1*(91.0d0*dx2**2+dx1*(42.0d0*dx2+5.0d0*dx1)))
      ci(16)=204.0d0*dx2**3*(91.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(45.0d0*dx2+4.0d0*dx1)))
      ci(17)=153.0d0*dx2**2*(56.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(16.0d0*dx2+dx1)))
      ci(18)=9.0d0*dx2*(340.0d0*dx2**3+dx1*(272.0d0*dx2**2+dx1*(51.0d0*dx2+2.0d0*dx1)))
      ci(19)=816.0d0*dx2**3+dx1*(459.0d0*dx2**2+dx1*(54.0d0*dx2+dx1))
      ci(20)=3.0d0*(51.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))
      ci(21)=3.0d0*(6.0d0*dx2+dx1)
      ci(22)=1.0d0
    case(19)
      ci(1)=dx1**3*dx2**19
      ci(2)=dx1**2*dx2**18*(3.0d0*dx2+19.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**17*(dx2**2+19.0d0*dx1*(dx2+3.0d0*dx1))
      ci(4)=dx2**16*(dx2**3+57.0d0*dx1*(dx2**2+dx1*(9.0d0*dx2+17.0d0*dx1)))
      ci(5)=19.0d0*dx2**15*(dx2**3+3.0d0*dx1*(9.0d0*dx2**2+17.0d0*dx1*(3.0d0*dx2+4.0d0*dx1)))
      ci(6)=171.0d0*dx2**14*(dx2**3+17.0d0*dx1*(dx2**2+4.0d0*dx1*(dx2+dx1)))
      ci(7)=969.0d0*dx2**13*(dx2**3+4.0d0*dx1*(3.0d0*dx2**2+dx1*(9.0d0*dx2+7.0d0*dx1)))
      ci(8)=3876.0d0*dx2**12*(dx2**3+dx1*(9.0d0*dx2**2+dx1*(21.0d0*dx2+13.0d0*dx1)))
      ci(9)=5814.0d0*dx2**11*(2.0d0*dx2**3+dx1*(14.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+dx1)))
      ci(10)=646.0d0*dx2**10*(42.0d0*dx2**3+13.0d0*dx1*(18.0d0*dx2**2+dx1*(27.0d0*dx2+11.0d0*dx1)))
      ci(11)=8398.0d0*dx2**9*(6.0d0*dx2**3+dx1*(27.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(12)=25194.0d0*dx2**8*(3.0d0*dx2**3+dx1*(11.0d0*dx2**2+dx1*(11.0d0*dx2+3.0d0*dx1)))
      ci(13)=8398.0d0*dx2**7*(11.0d0*dx2**3+3.0d0*dx1*(11.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(14)=646.0d0*dx2**6*(143.0d0*dx2**3+3.0d0*dx1*(117.0d0*dx2**2+2.0d0*dx1*(39.0d0*dx2+7.0d0*dx1)))
      ci(15)=5814.0d0*dx2**5*(13.0d0*dx2**3+2.0d0*dx1*(13.0d0*dx2**2+dx1*(7.0d0*dx2+dx1)))
      ci(16)=3876.0d0*dx2**4*(13.0d0*dx2**3+dx1*(21.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(17)=969.0d0*dx2**3*(28.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(18)=171.0d0*dx2**2*(68.0d0*dx2**3+dx1*(68.0d0*dx2**2+dx1*(17.0d0*dx2+dx1)))
      ci(19)=19.0d0*dx2*(204.0d0*dx2**3+dx1*(153.0d0*dx2**2+dx1*(27.0d0*dx2+dx1)))
      ci(20)=969.0d0*dx2**3+dx1*(513.0d0*dx2**2+dx1*(57.0d0*dx2+dx1))
      ci(21)=3.0d0*(57.0d0*dx2**2+dx1*(19.0d0*dx2+dx1))
      ci(22)=19.0d0*dx2+3.0d0*dx1
      ci(23)=1.0d0
    case(20)
      ci(1)=dx1**3*dx2**20
      ci(2)=dx1**2*dx2**19*(3.0d0*dx2+20.0d0*dx1)
      ci(3)=dx1*dx2**18*(3.0d0*dx2**2+10.0d0*dx1*(6.0d0*dx2+19.0d0*dx1))
      ci(4)=dx2**17*(dx2**3+30.0d0*dx1*(2.0d0*dx2**2+19.0d0*dx1*(dx2+2.0d0*dx1)))
      ci(5)=5.0d0*dx2**16*(4.0d0*dx2**3+57.0d0*dx1*(2.0d0*dx2**2+dx1*(12.0d0*dx2+17.0d0*dx1)))
      ci(6)=19.0d0*dx2**15*(10.0d0*dx2**3+3.0d0*dx1*(60.0d0*dx2**2+17.0d0*dx1*(15.0d0*dx2+16.0d0*dx1)))
      ci(7)=57.0d0*dx2**14*(20.0d0*dx2**3+17.0d0*dx1*(15.0d0*dx2**2+8.0d0*dx1*(6.0d0*dx2+5.0d0*dx1)))
      ci(8)=969.0d0*dx2**13*(5.0d0*dx2**3+8.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(9)=1938.0d0*dx2**12*(8.0d0*dx2**3+5.0d0*dx1*(12.0d0*dx2**2+dx1*(24.0d0*dx2+13.0d0*dx1)))
      ci(10)=3230.0d0*dx2**11*(12.0d0*dx2**3+dx1*(72.0d0*dx2**2+13.0d0*dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(11)=646.0d0*dx2**10*(120.0d0*dx2**3+13.0d0*dx1*(45.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2+11.0d0*dx1)))
      ci(12)=8398.0d0*dx2**9*(15.0d0*dx2**3+2.0d0*dx1*(30.0d0*dx2**2+dx1*(33.0d0*dx2+10.0d0*dx1)))
      ci(13)=8398.0d0*dx2**8*(20.0d0*dx2**3+3.0d0*dx1*(22.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1)))
      ci(14)=646.0d0*dx2**7*(286.0d0*dx2**3+15.0d0*dx1*(52.0d0*dx2**2+dx1*(39.0d0*dx2+8.0d0*dx1)))
      ci(15)=3230.0d0*dx2**6*(52.0d0*dx2**3+3.0d0*dx1*(39.0d0*dx2**2+4.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(16)=1938.0d0*dx2**5*(65.0d0*dx2**3+4.0d0*dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(17)=969.0d0*dx2**4*(80.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(48.0d0*dx2+5.0d0*dx1)))
      ci(18)=57.0d0*dx2**3*(680.0d0*dx2**3+dx1*(816.0d0*dx2**2+5.0d0*dx1*(51.0d0*dx2+4.0d0*dx1)))
      ci(19)=19.0d0*dx2**2*(816.0d0*dx2**3+5.0d0*dx1*(153.0d0*dx2**2+2.0d0*dx1*(18.0d0*dx2+dx1)))
      ci(20)=5.0d0*dx2*(969.0d0*dx2**3+2.0d0*dx1*(342.0d0*dx2**2+dx1*(57.0d0*dx2+2.0d0*dx1)))
      ci(21)=1140.0d0*dx2**3+dx1*(570.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))
      ci(22)=190.0d0*dx2**2+3.0d0*dx1*(20.0d0*dx2+dx1)
      ci(23)=20.0d0*dx2+3.0d0*dx1
      ci(24)=1.0d0
    case(21)
      ci(1)=dx1**3*dx2**21
      ci(2)=3.0d0*dx1**2*dx2**20*(dx2+7.0d0*dx1)
      ci(3)=3.0d0*dx1*dx2**19*(dx2**2+7.0d0*dx1*(3.0d0*dx2+10.0d0*dx1))
      ci(4)=dx2**18*(dx2**3+7.0d0*dx1*(9.0d0*dx2**2+10.0d0*dx1*(9.0d0*dx2+19.0d0*dx1)))
      ci(5)=21.0d0*dx2**17*(dx2**3+5.0d0*dx1*(6.0d0*dx2**2+19.0d0*dx1*(2.0d0*dx2+3.0d0*dx1)))
      ci(6)=21.0d0*dx2**16*(10.0d0*dx2**3+19.0d0*dx1*(10.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+17.0d0*dx1)))
      ci(7)=133.0d0*dx2**15*(10.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+17.0d0*dx1*(9.0d0*dx2+8.0d0*dx1)))
      ci(8)=171.0d0*dx2**14*(35.0d0*dx2**3+17.0d0*dx1*(21.0d0*dx2**2+8.0d0*dx1*(7.0d0*dx2+5.0d0*dx1)))
      ci(9)=2907.0d0*dx2**13*(7.0d0*dx2**3+2.0d0*dx1*(28.0d0*dx2**2+5.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)))
      ci(10)=646.0d0*dx2**12*(84.0d0*dx2**3+5.0d0*dx1*(108.0d0*dx2**2+7.0d0*dx1*(27.0d0*dx2+13.0d0*dx1)))
      ci(11)=1938.0d0*dx2**11*(60.0d0*dx2**3+7.0d0*dx1*(45.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)))
      ci(12)=13566.0d0*dx2**10*(15.0d0*dx2**3+13.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(13)=58786.0d0*dx2**9*(5.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1)))
      ci(14)=13566.0d0*dx2**8*(26.0d0*dx2**3+dx1*(78.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+3.0d0*dx1)))
      ci(15)=1938.0d0*dx2**7*(182.0d0*dx2**3+5.0d0*dx1*(91.0d0*dx2**2+3.0d0*dx1*(21.0d0*dx2+4.0d0*dx1)))
      ci(16)=646.0d0*dx2**6*(455.0d0*dx2**3+3.0d0*dx1*(315.0d0*dx2**2+4.0d0*dx1*(45.0d0*dx2+7.0d0*dx1)))
      ci(17)=2907.0d0*dx2**5*(70.0d0*dx2**3+dx1*(120.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+dx1)))
      ci(18)=171.0d0*dx2**4*(680.0d0*dx2**3+7.0d0*dx1*(136.0d0*dx2**2+dx1*(51.0d0*dx2+5.0d0*dx1)))
      ci(19)=133.0d0*dx2**3*(408.0d0*dx2**3+dx1*(459.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+2.0d0*dx1)))
      ci(20)=21.0d0*dx2**2*(969.0d0*dx2**3+5.0d0*dx1*(171.0d0*dx2**2+2.0d0*dx1*(19.0d0*dx2+dx1)))
      ci(21)=21.0d0*dx2*(285.0d0*dx2**3+dx1*(190.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))
      ci(22)=1330.0d0*dx2**3+dx1*(630.0d0*dx2**2+dx1*(63.0d0*dx2+dx1))
      ci(23)=3.0d0*(70.0d0*dx2**2+dx1*(21.0d0*dx2+dx1))
      ci(24)=3.0d0*(7.0d0*dx2+dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>21, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^4 * (x-x2)^n2 as sum_k=0^(4+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_4(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**4
      ci(2)=4.0d0*dx1**3
      ci(3)=6.0d0*dx1**2
      ci(4)=4.0d0*dx1
      ci(5)=1.0d0
    case(1)
      ci(1)=dx1**4*dx2
      ci(2)=dx1**3*(4.0d0*dx2+dx1)
      ci(3)=2.0d0*dx1**2*(3.0d0*dx2+2.0d0*dx1)
      ci(4)=2.0d0*dx1*(2.0d0*dx2+3.0d0*dx1)
      ci(5)=dx2+4.0d0*dx1
      ci(6)=1.0d0
    case(2)
      ci(1)=dx1**4*dx2**2
      ci(2)=2.0d0*dx1**3*dx2*(2.0d0*dx2+dx1)
      ci(3)=dx1**2*(6.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(4)=4.0d0*dx1*(dx2**2+dx1*(3.0d0*dx2+dx1))
      ci(5)=dx2**2+2.0d0*dx1*(4.0d0*dx2+3.0d0*dx1)
      ci(6)=2.0d0*(dx2+2.0d0*dx1)
      ci(7)=1.0d0
    case(3)
      ci(1)=dx1**4*dx2**3
      ci(2)=dx1**3*dx2**2*(4.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**2*dx2*(2.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))
      ci(4)=dx1*(4.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(5)=dx2**3+2.0d0*dx1*(6.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1))
      ci(6)=3.0d0*(dx2**2+2.0d0*dx1*(2.0d0*dx2+dx1))
      ci(7)=3.0d0*dx2+4.0d0*dx1
      ci(8)=1.0d0
    case(4)
      ci(1)=dx1**4*dx2**4
      ci(2)=4.0d0*dx1**3*dx2**3*(dx2+dx1)
      ci(3)=2.0d0*dx1**2*dx2**2*(3.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2*(dx2**3+dx1*(6.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))
      ci(5)=dx2**4+dx1*(16.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(16.0d0*dx2+dx1)))
      ci(6)=4.0d0*(dx2**3+dx1*(6.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))
      ci(7)=2.0d0*(3.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(8)=4.0d0*(dx2+dx1)
      ci(9)=1.0d0
    case(5)
      ci(1)=dx1**4*dx2**5
      ci(2)=dx1**3*dx2**4*(4.0d0*dx2+5.0d0*dx1)
      ci(3)=2.0d0*dx1**2*dx2**3*(3.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=2.0d0*dx1*dx2**2*(2.0d0*dx2**3+5.0d0*dx1*(3.0d0*dx2**2+dx1*(4.0d0*dx2+dx1)))
      ci(5)=dx2*(dx2**4+5.0d0*dx1*(4.0d0*dx2**3+dx1*(12.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(6)=5.0d0*dx2**4+dx1*(40.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)))
      ci(7)=2.0d0*(5.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(8)=2.0d0*(5.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1))
      ci(9)=5.0d0*dx2+4.0d0*dx1
      ci(10)=1.0d0
    case(6)
      ci(1)=dx1**4*dx2**6
      ci(2)=2.0d0*dx1**3*dx2**5*(2.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**2*dx2**4*(2.0d0*dx2**2+dx1*(8.0d0*dx2+5.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2**3*(dx2**3+dx1*(9.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(5)=dx2**2*(dx2**4+dx1*(24.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(6)=6.0d0*dx2*(dx2**4+dx1*(10.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(7)=15.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(24.0d0*dx2+dx1)))
      ci(8)=4.0d0*(5.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(9)=3.0d0*(5.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+dx1))
      ci(10)=2.0d0*(3.0d0*dx2+2.0d0*dx1)
      ci(11)=1.0d0
    case(7)
      ci(1)=dx1**4*dx2**7
      ci(2)=dx1**3*dx2**6*(4.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**2*dx2**5*(6.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(4)=dx1*dx2**4*(4.0d0*dx2**3+7.0d0*dx1*(6.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1)))
      ci(5)=dx2**3*(dx2**4+7.0d0*dx1*(4.0d0*dx2**3+dx1*(18.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx2**2*(dx2**4+dx1*(12.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1))))
      ci(7)=7.0d0*dx2*(3.0d0*dx2**4+dx1*(20.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(8)=35.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(28.0d0*dx2+dx1)))
      ci(9)=35.0d0*dx2**3+2.0d0*dx1*(42.0d0*dx2**2+dx1*(21.0d0*dx2+2.0d0*dx1))
      ci(10)=21.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+3.0d0*dx1)
      ci(11)=7.0d0*dx2+4.0d0*dx1
      ci(12)=1.0d0
    case(8)
      ci(1)=dx1**4*dx2**8
      ci(2)=4.0d0*dx1**3*dx2**7*(dx2+2.0d0*dx1)
      ci(3)=2.0d0*dx1**2*dx2**6*(3.0d0*dx2**2+2.0d0*dx1*(8.0d0*dx2+7.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2**5*(dx2**3+2.0d0*dx1*(6.0d0*dx2**2+7.0d0*dx1*(2.0d0*dx2+dx1)))
      ci(5)=dx2**4*(dx2**4+2.0d0*dx1*(16.0d0*dx2**3+7.0d0*dx1*(12.0d0*dx2**2+dx1*(16.0d0*dx2+5.0d0*dx1))))
      ci(6)=8.0d0*dx2**3*(dx2**4+7.0d0*dx1*(2.0d0*dx2**3+dx1*(6.0d0*dx2**2+dx1*(5.0d0*dx2+dx1))))
      ci(7)=28.0d0*dx2**2*(dx2**4+dx1*(8.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(8)=8.0d0*dx2*(7.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(14.0d0*dx2+dx1))))
      ci(9)=70.0d0*dx2**4+dx1*(224.0d0*dx2**3+dx1*(168.0d0*dx2**2+dx1*(32.0d0*dx2+dx1)))
      ci(10)=4.0d0*(14.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(11)=2.0d0*(14.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))
      ci(12)=4.0d0*(2.0d0*dx2+dx1)
      ci(13)=1.0d0
    case(9)
      ci(1)=dx1**4*dx2**9
      ci(2)=dx1**3*dx2**8*(4.0d0*dx2+9.0d0*dx1)
      ci(3)=6.0d0*dx1**2*dx2**7*(dx2**2+6.0d0*dx1*(dx2+dx1))
      ci(4)=2.0d0*dx1*dx2**6*(2.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)))
      ci(5)=dx2**5*(dx2**4+6.0d0*dx1*(6.0d0*dx2**3+dx1*(36.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(6)=9.0d0*dx2**4*(dx2**4+2.0d0*dx1*(8.0d0*dx2**3+7.0d0*dx1*(4.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))))
      ci(7)=12.0d0*dx2**3*(3.0d0*dx2**4+7.0d0*dx1*(4.0d0*dx2**3+dx1*(9.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))
      ci(8)=12.0d0*dx2**2*(7.0d0*dx2**4+dx1*(42.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1))))
      ci(9)=9.0d0*dx2*(14.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(56.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(10)=126.0d0*dx2**4+dx1*(336.0d0*dx2**3+dx1*(216.0d0*dx2**2+dx1*(36.0d0*dx2+dx1)))
      ci(11)=2.0d0*(42.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(27.0d0*dx2+2.0d0*dx1)))
      ci(12)=6.0d0*(6.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(13)=9.0d0*dx2+4.0d0*dx1
      ci(14)=1.0d0
    case(10)
      ci(1)=dx1**4*dx2**10
      ci(2)=2.0d0*dx1**3*dx2**9*(2.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**2*dx2**8*(6.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+9.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2**7*(dx2**3+15.0d0*dx1*(dx2**2+dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx2**6*(dx2**4+10.0d0*dx1*(4.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+dx1*(16.0d0*dx2+7.0d0*dx1))))
      ci(6)=2.0d0*dx2**5*(5.0d0*dx2**4+6.0d0*dx1*(15.0d0*dx2**3+dx1*(60.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+3.0d0*dx1))))
      ci(7)=3.0d0*dx2**4*(15.0d0*dx2**4+2.0d0*dx1*(80.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1))))
      ci(8)=24.0d0*dx2**3*(5.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1*(63.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+dx1))))
      ci(9)=3.0d0*dx2**2*(70.0d0*dx2**4+dx1*(336.0d0*dx2**3+5.0d0*dx1*(84.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1))))
      ci(10)=2.0d0*dx2*(126.0d0*dx2**4+5.0d0*dx1*(84.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))))
      ci(11)=210.0d0*dx2**4+dx1*(480.0d0*dx2**3+dx1*(270.0d0*dx2**2+dx1*(40.0d0*dx2+dx1)))
      ci(12)=4.0d0*(30.0d0*dx2**3+dx1*(45.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))
      ci(13)=45.0d0*dx2**2+2.0d0*dx1*(20.0d0*dx2+3.0d0*dx1)
      ci(14)=2.0d0*(5.0d0*dx2+2.0d0*dx1)
      ci(15)=1.0d0
    case(11)
      ci(1)=dx1**4*dx2**11
      ci(2)=dx1**3*dx2**10*(4.0d0*dx2+11.0d0*dx1)
      ci(3)=dx1**2*dx2**9*(6.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+5.0d0*dx1))
      ci(4)=dx1*dx2**8*(4.0d0*dx2**3+11.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+3.0d0*dx1)))
      ci(5)=dx2**7*(dx2**4+22.0d0*dx1*(2.0d0*dx2**3+15.0d0*dx1*(dx2**2+dx1*(2.0d0*dx2+dx1))))
      ci(6)=11.0d0*dx2**6*(dx2**4+2.0d0*dx1*(10.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+dx1*(20.0d0*dx2+7.0d0*dx1))))
      ci(7)=11.0d0*dx2**5*(5.0d0*dx2**4+6.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(8)=33.0d0*dx2**4*(5.0d0*dx2**4+2.0d0*dx1*(20.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(28.0d0*dx2+5.0d0*dx1))))
      ci(9)=33.0d0*dx2**3*(10.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(84.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(10)=11.0d0*dx2**2*(42.0d0*dx2**4+dx1*(168.0d0*dx2**3+5.0d0*dx1*(36.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(11)=11.0d0*dx2*(42.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(12)=330.0d0*dx2**4+dx1*(660.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1*(44.0d0*dx2+dx1)))
      ci(13)=165.0d0*dx2**3+2.0d0*dx1*(110.0d0*dx2**2+dx1*(33.0d0*dx2+2.0d0*dx1))
      ci(14)=55.0d0*dx2**2+2.0d0*dx1*(22.0d0*dx2+3.0d0*dx1)
      ci(15)=11.0d0*dx2+4.0d0*dx1
      ci(16)=1.0d0
    case(12)
      ci(1)=dx1**4*dx2**12
      ci(2)=4.0d0*dx1**3*dx2**11*(dx2+3.0d0*dx1)
      ci(3)=6.0d0*dx1**2*dx2**10*(dx2**2+dx1*(8.0d0*dx2+11.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2**9*(dx2**3+dx1*(18.0d0*dx2**2+11.0d0*dx1*(6.0d0*dx2+5.0d0*dx1)))
      ci(5)=dx2**8*(dx2**4+dx1*(48.0d0*dx2**3+11.0d0*dx1*(36.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+9.0d0*dx1))))
      ci(6)=12.0d0*dx2**7*(dx2**4+11.0d0*dx1*(2.0d0*dx2**3+dx1*(10.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))))
      ci(7)=22.0d0*dx2**6*(3.0d0*dx2**4+dx1*(40.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+2.0d0*dx1*(24.0d0*dx2+7.0d0*dx1))))
      ci(8)=44.0d0*dx2**5*(5.0d0*dx2**4+3.0d0*dx1*(15.0d0*dx2**3+2.0d0*dx1*(18.0d0*dx2**2+dx1*(14.0d0*dx2+3.0d0*dx1))))
      ci(9)=99.0d0*dx2**4*(5.0d0*dx2**4+dx1*(32.0d0*dx2**3+dx1*(56.0d0*dx2**2+dx1*(32.0d0*dx2+5.0d0*dx1))))
      ci(10)=44.0d0*dx2**3*(18.0d0*dx2**4+dx1*(84.0d0*dx2**3+dx1*(108.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+dx1))))
      ci(11)=22.0d0*dx2**2*(42.0d0*dx2**4+dx1*(144.0d0*dx2**3+dx1*(135.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))))
      ci(12)=12.0d0*dx2*(66.0d0*dx2**4+dx1*(165.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))))
      ci(13)=495.0d0*dx2**4+dx1*(880.0d0*dx2**3+dx1*(396.0d0*dx2**2+dx1*(48.0d0*dx2+dx1)))
      ci(14)=4.0d0*(55.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))
      ci(15)=6.0d0*(11.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(16)=4.0d0*(3.0d0*dx2+dx1)
      ci(17)=1.0d0
    case(13)
      ci(1)=dx1**4*dx2**13
      ci(2)=dx1**3*dx2**12*(4.0d0*dx2+13.0d0*dx1)
      ci(3)=2.0d0*dx1**2*dx2**11*(3.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(4)=2.0d0*dx1*dx2**10*(2.0d0*dx2**3+13.0d0*dx1*(3.0d0*dx2**2+dx1*(12.0d0*dx2+11.0d0*dx1)))
      ci(5)=dx2**9*(dx2**4+13.0d0*dx1*(4.0d0*dx2**3+dx1*(36.0d0*dx2**2+11.0d0*dx1*(8.0d0*dx2+5.0d0*dx1))))
      ci(6)=13.0d0*dx2**8*(dx2**4+dx1*(24.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+dx1*(20.0d0*dx2+9.0d0*dx1))))
      ci(7)=26.0d0*dx2**7*(3.0d0*dx2**4+11.0d0*dx1*(4.0d0*dx2**3+3.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))))
      ci(8)=286.0d0*dx2**6*(dx2**4+dx1*(10.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(9)=143.0d0*dx2**5*(5.0d0*dx2**4+3.0d0*dx1*(12.0d0*dx2**3+dx1*(24.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(10)=143.0d0*dx2**4*(9.0d0*dx2**4+dx1*(48.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))))
      ci(11)=286.0d0*dx2**3*(6.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(12)=26.0d0*dx2**2*(66.0d0*dx2**4+dx1*(198.0d0*dx2**3+dx1*(165.0d0*dx2**2+dx1*(44.0d0*dx2+3.0d0*dx1))))
      ci(13)=13.0d0*dx2*(99.0d0*dx2**4+dx1*(220.0d0*dx2**3+dx1*(132.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(14)=715.0d0*dx2**4+dx1*(1144.0d0*dx2**3+dx1*(468.0d0*dx2**2+dx1*(52.0d0*dx2+dx1)))
      ci(15)=2.0d0*(143.0d0*dx2**3+dx1*(156.0d0*dx2**2+dx1*(39.0d0*dx2+2.0d0*dx1)))
      ci(16)=2.0d0*(39.0d0*dx2**2+dx1*(26.0d0*dx2+3.0d0*dx1))
      ci(17)=13.0d0*dx2+4.0d0*dx1
      ci(18)=1.0d0
    case(14)
      ci(1)=dx1**4*dx2**14
      ci(2)=2.0d0*dx1**3*dx2**13*(2.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**2*dx2**12*(6.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+13.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2**11*(dx2**3+7.0d0*dx1*(3.0d0*dx2**2+13.0d0*dx1*(dx2+dx1)))
      ci(5)=dx2**10*(dx2**4+7.0d0*dx1*(8.0d0*dx2**3+13.0d0*dx1*(6.0d0*dx2**2+dx1*(16.0d0*dx2+11.0d0*dx1))))
      ci(6)=14.0d0*dx2**9*(dx2**4+13.0d0*dx1*(2.0d0*dx2**3+dx1*(12.0d0*dx2**2+11.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(7)=91.0d0*dx2**8*(dx2**4+dx1*(16.0d0*dx2**3+11.0d0*dx1*(6.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(8)=52.0d0*dx2**7*(7.0d0*dx2**4+11.0d0*dx1*(7.0d0*dx2**3+3.0d0*dx1*(7.0d0*dx2**2+dx1*(7.0d0*dx2+2.0d0*dx1))))
      ci(9)=143.0d0*dx2**6*(7.0d0*dx2**4+dx1*(56.0d0*dx2**3+3.0d0*dx1*(42.0d0*dx2**2+dx1*(32.0d0*dx2+7.0d0*dx1))))
      ci(10)=286.0d0*dx2**5*(7.0d0*dx2**4+dx1*(42.0d0*dx2**3+dx1*(72.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+dx1))))
      ci(11)=143.0d0*dx2**4*(21.0d0*dx2**4+dx1*(96.0d0*dx2**3+7.0d0*dx1*(18.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(12)=52.0d0*dx2**3*(66.0d0*dx2**4+7.0d0*dx1*(33.0d0*dx2**3+dx1*(33.0d0*dx2**2+dx1*(11.0d0*dx2+dx1))))
      ci(13)=91.0d0*dx2**2*(33.0d0*dx2**4+dx1*(88.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(14)=14.0d0*dx2*(143.0d0*dx2**4+dx1*(286.0d0*dx2**3+dx1*(156.0d0*dx2**2+dx1*(26.0d0*dx2+dx1))))
      ci(15)=1001.0d0*dx2**4+dx1*(1456.0d0*dx2**3+dx1*(546.0d0*dx2**2+dx1*(56.0d0*dx2+dx1)))
      ci(16)=4.0d0*(91.0d0*dx2**3+dx1*(91.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)))
      ci(17)=91.0d0*dx2**2+2.0d0*dx1*(28.0d0*dx2+3.0d0*dx1)
      ci(18)=2.0d0*(7.0d0*dx2+2.0d0*dx1)
      ci(19)=1.0d0
    case(15)
      ci(1)=dx1**4*dx2**15
      ci(2)=dx1**3*dx2**14*(4.0d0*dx2+15.0d0*dx1)
      ci(3)=3.0d0*dx1**2*dx2**13*(2.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+7.0d0*dx1))
      ci(4)=dx1*dx2**12*(4.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+13.0d0*dx1)))
      ci(5)=dx2**11*(dx2**4+5.0d0*dx1*(12.0d0*dx2**3+7.0d0*dx1*(18.0d0*dx2**2+13.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))))
      ci(6)=3.0d0*dx2**10*(5.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+13.0d0*dx1*(10.0d0*dx2**2+dx1*(20.0d0*dx2+11.0d0*dx1))))
      ci(7)=7.0d0*dx2**9*(15.0d0*dx2**4+13.0d0*dx1*(20.0d0*dx2**3+dx1*(90.0d0*dx2**2+11.0d0*dx1*(12.0d0*dx2+5.0d0*dx1))))
      ci(8)=13.0d0*dx2**8*(35.0d0*dx2**4+dx1*(420.0d0*dx2**3+11.0d0*dx1*(126.0d0*dx2**2+5.0d0*dx1*(28.0d0*dx2+9.0d0*dx1))))
      ci(9)=39.0d0*dx2**7*(35.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+5.0d0*dx1*(14.0d0*dx2**2+3.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(10)=143.0d0*dx2**6*(21.0d0*dx2**4+5.0d0*dx1*(28.0d0*dx2**3+dx1*(54.0d0*dx2**2+dx1*(36.0d0*dx2+7.0d0*dx1))))
      ci(11)=143.0d0*dx2**5*(35.0d0*dx2**4+dx1*(180.0d0*dx2**3+dx1*(270.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+3.0d0*dx1))))
      ci(12)=39.0d0*dx2**4*(165.0d0*dx2**4+dx1*(660.0d0*dx2**3+7.0d0*dx1*(110.0d0*dx2**2+dx1*(44.0d0*dx2+5.0d0*dx1))))
      ci(13)=13.0d0*dx2**3*(495.0d0*dx2**4+7.0d0*dx1*(220.0d0*dx2**3+dx1*(198.0d0*dx2**2+5.0d0*dx1*(12.0d0*dx2+dx1))))
      ci(14)=7.0d0*dx2**2*(715.0d0*dx2**4+dx1*(1716.0d0*dx2**3+5.0d0*dx1*(234.0d0*dx2**2+dx1*(52.0d0*dx2+3.0d0*dx1))))
      ci(15)=3.0d0*dx2*(1001.0d0*dx2**4+5.0d0*dx1*(364.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(28.0d0*dx2+dx1))))
      ci(16)=1365.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1*(630.0d0*dx2**2+dx1*(60.0d0*dx2+dx1)))
      ci(17)=455.0d0*dx2**3+2.0d0*dx1*(210.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1))
      ci(18)=3.0d0*(35.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+dx1))
      ci(19)=15.0d0*dx2+4.0d0*dx1
      ci(20)=1.0d0
    case(16)
      ci(1)=dx1**4*dx2**16
      ci(2)=4.0d0*dx1**3*dx2**15*(dx2+4.0d0*dx1)
      ci(3)=2.0d0*dx1**2*dx2**14*(3.0d0*dx2**2+4.0d0*dx1*(8.0d0*dx2+15.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2**13*(dx2**3+4.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+7.0d0*dx1)))
      ci(5)=dx2**12*(dx2**4+4.0d0*dx1*(16.0d0*dx2**3+5.0d0*dx1*(36.0d0*dx2**2+7.0d0*dx1*(16.0d0*dx2+13.0d0*dx1))))
      ci(6)=16.0d0*dx2**11*(dx2**4+dx1*(30.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))))
      ci(7)=8.0d0*dx2**10*(15.0d0*dx2**4+7.0d0*dx1*(40.0d0*dx2**3+13.0d0*dx1*(15.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1))))
      ci(8)=16.0d0*dx2**9*(35.0d0*dx2**4+13.0d0*dx1*(35.0d0*dx2**3+dx1*(126.0d0*dx2**2+11.0d0*dx1*(14.0d0*dx2+5.0d0*dx1))))
      ci(9)=26.0d0*dx2**8*(70.0d0*dx2**4+dx1*(672.0d0*dx2**3+11.0d0*dx1*(168.0d0*dx2**2+5.0d0*dx1*(32.0d0*dx2+9.0d0*dx1))))
      ci(10)=104.0d0*dx2**7*(42.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+5.0d0*dx1*(12.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1))))
      ci(11)=572.0d0*dx2**6*(14.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(135.0d0*dx2**2+2.0d0*dx1*(40.0d0*dx2+7.0d0*dx1))))
      ci(12)=104.0d0*dx2**5*(110.0d0*dx2**4+dx1*(495.0d0*dx2**3+2.0d0*dx1*(330.0d0*dx2**2+7.0d0*dx1*(22.0d0*dx2+3.0d0*dx1))))
      ci(13)=26.0d0*dx2**4*(495.0d0*dx2**4+2.0d0*dx1*(880.0d0*dx2**3+7.0d0*dx1*(132.0d0*dx2**2+dx1*(48.0d0*dx2+5.0d0*dx1))))
      ci(14)=16.0d0*dx2**3*(715.0d0*dx2**4+7.0d0*dx1*(286.0d0*dx2**3+dx1*(234.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+dx1))))
      ci(15)=8.0d0*dx2**2*(1001.0d0*dx2**4+dx1*(2184.0d0*dx2**3+5.0d0*dx1*(273.0d0*dx2**2+dx1*(56.0d0*dx2+3.0d0*dx1))))
      ci(16)=16.0d0*dx2*(273.0d0*dx2**4+dx1*(455.0d0*dx2**3+dx1*(210.0d0*dx2**2+dx1*(30.0d0*dx2+dx1))))
      ci(17)=1820.0d0*dx2**4+dx1*(2240.0d0*dx2**3+dx1*(720.0d0*dx2**2+dx1*(64.0d0*dx2+dx1)))
      ci(18)=4.0d0*(140.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(24.0d0*dx2+dx1)))
      ci(19)=2.0d0*(60.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1))
      ci(20)=4.0d0*(4.0d0*dx2+dx1)
      ci(21)=1.0d0
    case(17)
      ci(1)=dx1**4*dx2**17
      ci(2)=dx1**3*dx2**16*(4.0d0*dx2+17.0d0*dx1)
      ci(3)=2.0d0*dx1**2*dx2**15*(3.0d0*dx2**2+34.0d0*dx1*(dx2+2.0d0*dx1))
      ci(4)=2.0d0*dx1*dx2**14*(2.0d0*dx2**3+17.0d0*dx1*(3.0d0*dx2**2+4.0d0*dx1*(4.0d0*dx2+5.0d0*dx1)))
      ci(5)=dx2**13*(dx2**4+68.0d0*dx1*(dx2**3+dx1*(12.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+7.0d0*dx1))))
      ci(6)=17.0d0*dx2**12*(dx2**4+4.0d0*dx1*(8.0d0*dx2**3+dx1*(60.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+13.0d0*dx1))))
      ci(7)=136.0d0*dx2**11*(dx2**4+dx1*(20.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(8)=136.0d0*dx2**10*(5.0d0*dx2**4+dx1*(70.0d0*dx2**3+13.0d0*dx1*(21.0d0*dx2**2+dx1*(28.0d0*dx2+11.0d0*dx1))))
      ci(9)=34.0d0*dx2**9*(70.0d0*dx2**4+13.0d0*dx1*(56.0d0*dx2**3+dx1*(168.0d0*dx2**2+11.0d0*dx1*(16.0d0*dx2+5.0d0*dx1))))
      ci(10)=442.0d0*dx2**8*(14.0d0*dx2**4+dx1*(112.0d0*dx2**3+11.0d0*dx1*(24.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(11)=884.0d0*dx2**7*(14.0d0*dx2**4+11.0d0*dx1*(8.0d0*dx2**3+dx1*(15.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1))))
      ci(12)=884.0d0*dx2**6*(22.0d0*dx2**4+dx1*(110.0d0*dx2**3+dx1*(165.0d0*dx2**2+2.0d0*dx1*(44.0d0*dx2+7.0d0*dx1))))
      ci(13)=442.0d0*dx2**5*(55.0d0*dx2**4+2.0d0*dx1*(110.0d0*dx2**3+dx1*(132.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(14)=34.0d0*dx2**4*(715.0d0*dx2**4+2.0d0*dx1*(1144.0d0*dx2**3+7.0d0*dx1*(156.0d0*dx2**2+dx1*(52.0d0*dx2+5.0d0*dx1))))
      ci(15)=136.0d0*dx2**3*(143.0d0*dx2**4+dx1*(364.0d0*dx2**3+dx1*(273.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+dx1))))
      ci(16)=136.0d0*dx2**2*(91.0d0*dx2**4+dx1*(182.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(17)=17.0d0*dx2*(364.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))
      ci(18)=2380.0d0*dx2**4+dx1*(2720.0d0*dx2**3+dx1*(816.0d0*dx2**2+dx1*(68.0d0*dx2+dx1)))
      ci(19)=2.0d0*(340.0d0*dx2**3+dx1*(272.0d0*dx2**2+dx1*(51.0d0*dx2+2.0d0*dx1)))
      ci(20)=2.0d0*(68.0d0*dx2**2+dx1*(34.0d0*dx2+3.0d0*dx1))
      ci(21)=17.0d0*dx2+4.0d0*dx1
      ci(22)=1.0d0
    case(18)
      ci(1)=dx1**4*dx2**18
      ci(2)=2.0d0*dx1**3*dx2**17*(2.0d0*dx2+9.0d0*dx1)
      ci(3)=3.0d0*dx1**2*dx2**16*(2.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+17.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2**15*(dx2**3+3.0d0*dx1*(9.0d0*dx2**2+17.0d0*dx1*(3.0d0*dx2+4.0d0*dx1)))
      ci(5)=dx2**14*(dx2**4+6.0d0*dx1*(12.0d0*dx2**3+17.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(16.0d0*dx2+15.0d0*dx1))))
      ci(6)=18.0d0*dx2**13*(dx2**4+34.0d0*dx1*(dx2**3+2.0d0*dx1*(4.0d0*dx2**2+dx1*(10.0d0*dx2+7.0d0*dx1))))
      ci(7)=51.0d0*dx2**12*(3.0d0*dx2**4+4.0d0*dx1*(16.0d0*dx2**3+dx1*(90.0d0*dx2**2+7.0d0*dx1*(24.0d0*dx2+13.0d0*dx1))))
      ci(8)=816.0d0*dx2**11*(dx2**4+dx1*(15.0d0*dx2**3+dx1*(63.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+3.0d0*dx1))))
      ci(9)=306.0d0*dx2**10*(10.0d0*dx2**4+dx1*(112.0d0*dx2**3+13.0d0*dx1*(28.0d0*dx2**2+dx1*(32.0d0*dx2+11.0d0*dx1))))
      ci(10)=68.0d0*dx2**9*(126.0d0*dx2**4+13.0d0*dx1*(84.0d0*dx2**3+dx1*(216.0d0*dx2**2+11.0d0*dx1*(18.0d0*dx2+5.0d0*dx1))))
      ci(11)=442.0d0*dx2**8*(42.0d0*dx2**4+dx1*(288.0d0*dx2**3+11.0d0*dx1*(54.0d0*dx2**2+dx1*(40.0d0*dx2+9.0d0*dx1))))
      ci(12)=5304.0d0*dx2**7*(6.0d0*dx2**4+dx1*(33.0d0*dx2**3+dx1*(55.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+2.0d0*dx1))))
      ci(13)=442.0d0*dx2**6*(99.0d0*dx2**4+2.0d0*dx1*(220.0d0*dx2**3+3.0d0*dx1*(99.0d0*dx2**2+dx1*(48.0d0*dx2+7.0d0*dx1))))
      ci(14)=68.0d0*dx2**5*(715.0d0*dx2**4+6.0d0*dx1*(429.0d0*dx2**3+dx1*(468.0d0*dx2**2+7.0d0*dx1*(26.0d0*dx2+3.0d0*dx1))))
      ci(15)=306.0d0*dx2**4*(143.0d0*dx2**4+2.0d0*dx1*(208.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(56.0d0*dx2+5.0d0*dx1))))
      ci(16)=816.0d0*dx2**3*(39.0d0*dx2**4+dx1*(91.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))))
      ci(17)=51.0d0*dx2**2*(364.0d0*dx2**4+dx1*(672.0d0*dx2**3+dx1*(360.0d0*dx2**2+dx1*(64.0d0*dx2+3.0d0*dx1))))
      ci(18)=18.0d0*dx2*(476.0d0*dx2**4+dx1*(680.0d0*dx2**3+dx1*(272.0d0*dx2**2+dx1*(34.0d0*dx2+dx1))))
      ci(19)=3060.0d0*dx2**4+dx1*(3264.0d0*dx2**3+dx1*(918.0d0*dx2**2+dx1*(72.0d0*dx2+dx1)))
      ci(20)=4.0d0*(204.0d0*dx2**3+dx1*(153.0d0*dx2**2+dx1*(27.0d0*dx2+dx1)))
      ci(21)=3.0d0*(51.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+dx1))
      ci(22)=2.0d0*(9.0d0*dx2+2.0d0*dx1)
      ci(23)=1.0d0
    case(19)
      ci(1)=dx1**4*dx2**19
      ci(2)=dx1**3*dx2**18*(4.0d0*dx2+19.0d0*dx1)
      ci(3)=dx1**2*dx2**17*(6.0d0*dx2**2+19.0d0*dx1*(4.0d0*dx2+9.0d0*dx1))
      ci(4)=dx1*dx2**16*(4.0d0*dx2**3+57.0d0*dx1*(2.0d0*dx2**2+dx1*(12.0d0*dx2+17.0d0*dx1)))
      ci(5)=dx2**15*(dx2**4+38.0d0*dx1*(2.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+34.0d0*dx1*(dx2+dx1))))
      ci(6)=19.0d0*dx2**14*(dx2**4+6.0d0*dx1*(6.0d0*dx2**3+17.0d0*dx1*(3.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))))
      ci(7)=57.0d0*dx2**13*(3.0d0*dx2**4+68.0d0*dx1*(dx2**3+dx1*(6.0d0*dx2**2+dx1*(12.0d0*dx2+7.0d0*dx1))))
      ci(8)=969.0d0*dx2**12*(dx2**4+4.0d0*dx1*(4.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(28.0d0*dx2+13.0d0*dx1))))
      ci(9)=1938.0d0*dx2**11*(2.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(84.0d0*dx2**2+13.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(10)=646.0d0*dx2**10*(18.0d0*dx2**4+dx1*(168.0d0*dx2**3+13.0d0*dx1*(36.0d0*dx2**2+dx1*(36.0d0*dx2+11.0d0*dx1))))
      ci(11)=646.0d0*dx2**9*(42.0d0*dx2**4+13.0d0*dx1*(24.0d0*dx2**3+dx1*(54.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(12)=8398.0d0*dx2**8*(6.0d0*dx2**4+dx1*(36.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(44.0d0*dx2+9.0d0*dx1))))
      ci(13)=8398.0d0*dx2**7*(9.0d0*dx2**4+2.0d0*dx1*(22.0d0*dx2**3+3.0d0*dx1*(11.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))
      ci(14)=646.0d0*dx2**6*(143.0d0*dx2**4+2.0d0*dx1*(286.0d0*dx2**3+3.0d0*dx1*(117.0d0*dx2**2+dx1*(52.0d0*dx2+7.0d0*dx1))))
      ci(15)=646.0d0*dx2**5*(143.0d0*dx2**4+6.0d0*dx1*(78.0d0*dx2**3+dx1*(78.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1))))
      ci(16)=1938.0d0*dx2**4*(39.0d0*dx2**4+2.0d0*dx1*(52.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(17)=969.0d0*dx2**3*(52.0d0*dx2**4+dx1*(112.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(18)=57.0d0*dx2**2*(476.0d0*dx2**4+dx1*(816.0d0*dx2**3+dx1*(408.0d0*dx2**2+dx1*(68.0d0*dx2+3.0d0*dx1))))
      ci(19)=19.0d0*dx2*(612.0d0*dx2**4+dx1*(816.0d0*dx2**3+dx1*(306.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))
      ci(20)=3876.0d0*dx2**4+dx1*(3876.0d0*dx2**3+dx1*(1026.0d0*dx2**2+dx1*(76.0d0*dx2+dx1)))
      ci(21)=969.0d0*dx2**3+2.0d0*dx1*(342.0d0*dx2**2+dx1*(57.0d0*dx2+2.0d0*dx1))
      ci(22)=171.0d0*dx2**2+2.0d0*dx1*(38.0d0*dx2+3.0d0*dx1)
      ci(23)=19.0d0*dx2+4.0d0*dx1
      ci(24)=1.0d0
    case(20)
      ci(1)=dx1**4*dx2**20
      ci(2)=4.0d0*dx1**3*dx2**19*(dx2+5.0d0*dx1)
      ci(3)=2.0d0*dx1**2*dx2**18*(3.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+19.0d0*dx1))
      ci(4)=4.0d0*dx1*dx2**17*(dx2**3+5.0d0*dx1*(6.0d0*dx2**2+19.0d0*dx1*(2.0d0*dx2+3.0d0*dx1)))
      ci(5)=dx2**16*(dx2**4+5.0d0*dx1*(16.0d0*dx2**3+57.0d0*dx1*(4.0d0*dx2**2+dx1*(16.0d0*dx2+17.0d0*dx1))))
      ci(6)=4.0d0*dx2**15*(5.0d0*dx2**4+19.0d0*dx1*(10.0d0*dx2**3+3.0d0*dx1*(30.0d0*dx2**2+17.0d0*dx1*(5.0d0*dx2+4.0d0*dx1))))
      ci(7)=38.0d0*dx2**14*(5.0d0*dx2**4+3.0d0*dx1*(40.0d0*dx2**3+17.0d0*dx1*(15.0d0*dx2**2+4.0d0*dx1*(8.0d0*dx2+5.0d0*dx1))))
      ci(8)=228.0d0*dx2**13*(5.0d0*dx2**4+17.0d0*dx1*(5.0d0*dx2**3+4.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(9)=969.0d0*dx2**12*(5.0d0*dx2**4+2.0d0*dx1*(32.0d0*dx2**3+5.0d0*dx1*(24.0d0*dx2**2+dx1*(32.0d0*dx2+13.0d0*dx1))))
      ci(10)=2584.0d0*dx2**11*(6.0d0*dx2**4+5.0d0*dx1*(12.0d0*dx2**3+dx1*(36.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+dx1))))
      ci(11)=1292.0d0*dx2**10*(30.0d0*dx2**4+dx1*(240.0d0*dx2**3+13.0d0*dx1*(45.0d0*dx2**2+dx1*(40.0d0*dx2+11.0d0*dx1))))
      ci(12)=2584.0d0*dx2**9*(30.0d0*dx2**4+13.0d0*dx1*(15.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(22.0d0*dx2+5.0d0*dx1))))
      ci(13)=8398.0d0*dx2**8*(15.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(132.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(14)=2584.0d0*dx2**7*(65.0d0*dx2**4+dx1*(286.0d0*dx2**3+15.0d0*dx1*(26.0d0*dx2**2+dx1*(13.0d0*dx2+2.0d0*dx1))))
      ci(15)=1292.0d0*dx2**6*(143.0d0*dx2**4+5.0d0*dx1*(104.0d0*dx2**3+3.0d0*dx1*(39.0d0*dx2**2+2.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(16)=2584.0d0*dx2**5*(65.0d0*dx2**4+3.0d0*dx1*(65.0d0*dx2**3+2.0d0*dx1*(30.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(17)=969.0d0*dx2**4*(130.0d0*dx2**4+dx1*(320.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(64.0d0*dx2+5.0d0*dx1))))
      ci(18)=228.0d0*dx2**3*(340.0d0*dx2**4+dx1*(680.0d0*dx2**3+dx1*(408.0d0*dx2**2+5.0d0*dx1*(17.0d0*dx2+dx1))))
      ci(19)=38.0d0*dx2**2*(1020.0d0*dx2**4+dx1*(1632.0d0*dx2**3+5.0d0*dx1*(153.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(20)=4.0d0*dx2*(3876.0d0*dx2**4+5.0d0*dx1*(969.0d0*dx2**3+dx1*(342.0d0*dx2**2+dx1*(38.0d0*dx2+dx1))))
      ci(21)=4845.0d0*dx2**4+dx1*(4560.0d0*dx2**3+dx1*(1140.0d0*dx2**2+dx1*(80.0d0*dx2+dx1)))
      ci(22)=4.0d0*(285.0d0*dx2**3+dx1*(190.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))
      ci(23)=2.0d0*(95.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))
      ci(24)=4.0d0*(5.0d0*dx2+dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>20, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^5 * (x-x2)^n2 as sum_k=0^(5+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_5(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**5
      ci(2)=5.0d0*dx1**4
      ci(3)=10.0d0*dx1**3
      ci(4)=10.0d0*dx1**2
      ci(5)=5.0d0*dx1
      ci(6)=1.0d0
    case(1)
      ci(1)=dx1**5*dx2
      ci(2)=dx1**4*(5.0d0*dx2+dx1)
      ci(3)=5.0d0*dx1**3*(2.0d0*dx2+dx1)
      ci(4)=10.0d0*dx1**2*(dx2+dx1)
      ci(5)=5.0d0*dx1*(dx2+2.0d0*dx1)
      ci(6)=dx2+5.0d0*dx1
      ci(7)=1.0d0
    case(2)
      ci(1)=dx1**5*dx2**2
      ci(2)=dx1**4*dx2*(5.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**3*(10.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))
      ci(4)=5.0d0*dx1**2*(2.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))
      ci(5)=5.0d0*dx1*(dx2**2+2.0d0*dx1*(2.0d0*dx2+dx1))
      ci(6)=dx2**2+10.0d0*dx1*(dx2+dx1)
      ci(7)=2.0d0*dx2+5.0d0*dx1
      ci(8)=1.0d0
    case(3)
      ci(1)=dx1**5*dx2**3
      ci(2)=dx1**4*dx2**2*(5.0d0*dx2+3.0d0*dx1)
      ci(3)=dx1**3*dx2*(10.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+dx1))
      ci(4)=dx1**2*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1*(dx2**3+dx1*(6.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))
      ci(6)=dx2**3+5.0d0*dx1*(3.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))
      ci(7)=3.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)
      ci(8)=3.0d0*dx2+5.0d0*dx1
      ci(9)=1.0d0
    case(4)
      ci(1)=dx1**5*dx2**4
      ci(2)=dx1**4*dx2**3*(5.0d0*dx2+4.0d0*dx1)
      ci(3)=2.0d0*dx1**3*dx2**2*(5.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1))
      ci(4)=2.0d0*dx1**2*dx2*(5.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1*(5.0d0*dx2**4+dx1*(40.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(6)=dx2**4+5.0d0*dx1*(4.0d0*dx2**3+dx1*(12.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))
      ci(7)=2.0d0*(2.0d0*dx2**3+5.0d0*dx1*(3.0d0*dx2**2+dx1*(4.0d0*dx2+dx1)))
      ci(8)=2.0d0*(3.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))
      ci(9)=4.0d0*dx2+5.0d0*dx1
      ci(10)=1.0d0
    case(5)
      ci(1)=dx1**5*dx2**5
      ci(2)=5.0d0*dx1**4*dx2**4*(dx2+dx1)
      ci(3)=5.0d0*dx1**3*dx2**3*(2.0d0*dx2**2+dx1*(5.0d0*dx2+2.0d0*dx1))
      ci(4)=10.0d0*dx1**2*dx2**2*(dx2**3+dx1*(5.0d0*dx2**2+dx1*(5.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1*dx2*(dx2**4+dx1*(10.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(6)=dx2**5+dx1*(25.0d0*dx2**4+dx1*(100.0d0*dx2**3+dx1*(100.0d0*dx2**2+dx1*(25.0d0*dx2+dx1))))
      ci(7)=5.0d0*(dx2**4+dx1*(10.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(8)=10.0d0*(dx2**3+dx1*(5.0d0*dx2**2+dx1*(5.0d0*dx2+dx1)))
      ci(9)=5.0d0*(2.0d0*dx2**2+dx1*(5.0d0*dx2+2.0d0*dx1))
      ci(10)=5.0d0*(dx2+dx1)
      ci(11)=1.0d0
    case(6)
      ci(1)=dx1**5*dx2**6
      ci(2)=dx1**4*dx2**5*(5.0d0*dx2+6.0d0*dx1)
      ci(3)=5.0d0*dx1**3*dx2**4*(2.0d0*dx2**2+3.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=5.0d0*dx1**2*dx2**3*(2.0d0*dx2**3+dx1*(12.0d0*dx2**2+dx1*(15.0d0*dx2+4.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**2*(dx2**4+dx1*(12.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1))))
      ci(6)=dx2*(dx2**5+dx1*(30.0d0*dx2**4+dx1*(150.0d0*dx2**3+dx1*(200.0d0*dx2**2+3.0d0*dx1*(25.0d0*dx2+2.0d0*dx1)))))
      ci(7)=6.0d0*dx2**5+dx1*(75.0d0*dx2**4+dx1*(200.0d0*dx2**3+dx1*(150.0d0*dx2**2+dx1*(30.0d0*dx2+dx1))))
      ci(8)=5.0d0*(3.0d0*dx2**4+dx1*(20.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(9)=5.0d0*(4.0d0*dx2**3+dx1*(15.0d0*dx2**2+2.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(10)=5.0d0*(3.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))
      ci(11)=6.0d0*dx2+5.0d0*dx1
      ci(12)=1.0d0
    case(7)
      ci(1)=dx1**5*dx2**7
      ci(2)=dx1**4*dx2**6*(5.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**3*dx2**5*(10.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))
      ci(4)=5.0d0*dx1**2*dx2**4*(2.0d0*dx2**3+7.0d0*dx1*(2.0d0*dx2**2+dx1*(3.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1*dx2**3*(dx2**4+7.0d0*dx1*(2.0d0*dx2**3+dx1*(6.0d0*dx2**2+dx1*(5.0d0*dx2+dx1))))
      ci(6)=dx2**2*(dx2**5+7.0d0*dx1*(5.0d0*dx2**4+dx1*(30.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(25.0d0*dx2+3.0d0*dx1)))))
      ci(7)=7.0d0*dx2*(dx2**5+dx1*(15.0d0*dx2**4+dx1*(50.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(8)=21.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(350.0d0*dx2**3+dx1*(210.0d0*dx2**2+dx1*(35.0d0*dx2+dx1))))
      ci(9)=5.0d0*(7.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(14.0d0*dx2+dx1))))
      ci(10)=5.0d0*(7.0d0*dx2**3+dx1*(21.0d0*dx2**2+2.0d0*dx1*(7.0d0*dx2+dx1)))
      ci(11)=21.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+2.0d0*dx1)
      ci(12)=7.0d0*dx2+5.0d0*dx1
      ci(13)=1.0d0
    case(8)
      ci(1)=dx1**5*dx2**8
      ci(2)=dx1**4*dx2**7*(5.0d0*dx2+8.0d0*dx1)
      ci(3)=2.0d0*dx1**3*dx2**6*(5.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+7.0d0*dx1))
      ci(4)=2.0d0*dx1**2*dx2**5*(5.0d0*dx2**3+2.0d0*dx1*(20.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**4*(dx2**4+2.0d0*dx1*(8.0d0*dx2**3+7.0d0*dx1*(4.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))))
      ci(6)=dx2**3*(dx2**5+2.0d0*dx1*(20.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+dx1*(40.0d0*dx2**2+dx1*(25.0d0*dx2+4.0d0*dx1)))))
      ci(7)=4.0d0*dx2**2*(2.0d0*dx2**5+7.0d0*dx1*(5.0d0*dx2**4+dx1*(20.0d0*dx2**3+dx1*(25.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(8)=4.0d0*dx2*(7.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(175.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1)))))
      ci(9)=56.0d0*dx2**5+dx1*(350.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(280.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))
      ci(10)=5.0d0*(14.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(56.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(11)=2.0d0*(28.0d0*dx2**3+5.0d0*dx1*(14.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))
      ci(12)=2.0d0*(14.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))
      ci(13)=8.0d0*dx2+5.0d0*dx1
      ci(14)=1.0d0
    case(9)
      ci(1)=dx1**5*dx2**9
      ci(2)=dx1**4*dx2**8*(5.0d0*dx2+9.0d0*dx1)
      ci(3)=dx1**3*dx2**7*(10.0d0*dx2**2+9.0d0*dx1*(5.0d0*dx2+4.0d0*dx1))
      ci(4)=2.0d0*dx1**2*dx2**6*(5.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+7.0d0*dx1)))
      ci(5)=dx1*dx2**5*(5.0d0*dx2**4+6.0d0*dx1*(15.0d0*dx2**3+dx1*(60.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+3.0d0*dx1))))
      ci(6)=dx2**4*(dx2**5+3.0d0*dx1*(15.0d0*dx2**4+2.0d0*dx1*(60.0d0*dx2**3+7.0d0*dx1*(20.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+dx1)))) &
)
      ci(7)=3.0d0*dx2**3*(3.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2 &
+2.0d0*dx1)))))
      ci(8)=12.0d0*dx2**2*(3.0d0*dx2**5+dx1*(35.0d0*dx2**4+dx1*(105.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1)))))
      ci(9)=3.0d0*dx2*(28.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(280.0d0*dx2**2+3.0d0*dx1*(20.0d0*dx2+dx1)))))
      ci(10)=126.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1*(840.0d0*dx2**3+dx1*(360.0d0*dx2**2+dx1*(45.0d0*dx2+dx1))))
      ci(11)=126.0d0*dx2**4+5.0d0*dx1*(84.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))
      ci(12)=2.0d0*(42.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(13)=36.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+2.0d0*dx1)
      ci(14)=9.0d0*dx2+5.0d0*dx1
      ci(15)=1.0d0
    case(10)
      ci(1)=dx1**5*dx2**10
      ci(2)=5.0d0*dx1**4*dx2**9*(dx2+2.0d0*dx1)
      ci(3)=5.0d0*dx1**3*dx2**8*(2.0d0*dx2**2+dx1*(10.0d0*dx2+9.0d0*dx1))
      ci(4)=5.0d0*dx1**2*dx2**7*(2.0d0*dx2**3+dx1*(20.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+8.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**6*(dx2**4+2.0d0*dx1*(10.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+dx1*(20.0d0*dx2+7.0d0*dx1))))
      ci(6)=dx2**5*(dx2**5+2.0d0*dx1*(25.0d0*dx2**4+3.0d0*dx1*(75.0d0*dx2**3+dx1*(200.0d0*dx2**2+7.0d0*dx1*(25.0d0*dx2+6.0d0*dx1)) &
)))
      ci(7)=5.0d0*dx2**4*(2.0d0*dx2**5+3.0d0*dx1*(15.0d0*dx2**4+2.0d0*dx1*(40.0d0*dx2**3+7.0d0*dx1*(10.0d0*dx2**2+dx1*(6.0d0*dx2 &
+dx1)))))
      ci(8)=15.0d0*dx2**3*(3.0d0*dx2**5+2.0d0*dx1*(20.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(35.0d0*dx2+4.0d0*dx1)) &
)))
      ci(9)=15.0d0*dx2**2*(8.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1)))))
      ci(10)=5.0d0*dx2*(42.0d0*dx2**5+dx1*(252.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))))
      ci(11)=252.0d0*dx2**5+dx1*(1050.0d0*dx2**4+dx1*(1200.0d0*dx2**3+dx1*(450.0d0*dx2**2+dx1*(50.0d0*dx2+dx1))))
      ci(12)=5.0d0*(42.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(13)=5.0d0*(24.0d0*dx2**3+dx1*(45.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+dx1)))
      ci(14)=5.0d0*(9.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1))
      ci(15)=5.0d0*(2.0d0*dx2+dx1)
      ci(16)=1.0d0
    case(11)
      ci(1)=dx1**5*dx2**11
      ci(2)=dx1**4*dx2**10*(5.0d0*dx2+11.0d0*dx1)
      ci(3)=5.0d0*dx1**3*dx2**9*(2.0d0*dx2**2+11.0d0*dx1*(dx2+dx1))
      ci(4)=5.0d0*dx1**2*dx2**8*(2.0d0*dx2**3+11.0d0*dx1*(2.0d0*dx2**2+dx1*(5.0d0*dx2+3.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**7*(dx2**4+11.0d0*dx1*(2.0d0*dx2**3+dx1*(10.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))))
      ci(6)=dx2**6*(dx2**5+11.0d0*dx1*(5.0d0*dx2**4+2.0d0*dx1*(25.0d0*dx2**3+3.0d0*dx1*(25.0d0*dx2**2+dx1*(25.0d0*dx2+7.0d0*dx1))) &
))
      ci(7)=11.0d0*dx2**5*(dx2**5+dx1*(25.0d0*dx2**4+6.0d0*dx1*(25.0d0*dx2**3+dx1*(50.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+dx1)))))
      ci(8)=55.0d0*dx2**4*(dx2**5+3.0d0*dx1*(5.0d0*dx2**4+2.0d0*dx1*(10.0d0*dx2**3+dx1*(14.0d0*dx2**2+dx1*(7.0d0*dx2+dx1)))))
      ci(9)=165.0d0*dx2**3*(dx2**5+dx1*(10.0d0*dx2**4+dx1*(28.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(10)=55.0d0*dx2**2*(6.0d0*dx2**5+dx1*(42.0d0*dx2**4+dx1*(84.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(11)=11.0d0*dx2*(42.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(300.0d0*dx2**3+dx1*(150.0d0*dx2**2+dx1*(25.0d0*dx2+dx1)))))
      ci(12)=462.0d0*dx2**5+dx1*(1650.0d0*dx2**4+dx1*(1650.0d0*dx2**3+dx1*(550.0d0*dx2**2+dx1*(55.0d0*dx2+dx1))))
      ci(13)=5.0d0*(66.0d0*dx2**4+dx1*(165.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))))
      ci(14)=5.0d0*(33.0d0*dx2**3+dx1*(55.0d0*dx2**2+2.0d0*dx1*(11.0d0*dx2+dx1)))
      ci(15)=5.0d0*(11.0d0*dx2**2+dx1*(11.0d0*dx2+2.0d0*dx1))
      ci(16)=11.0d0*dx2+5.0d0*dx1
      ci(17)=1.0d0
    case(12)
      ci(1)=dx1**5*dx2**12
      ci(2)=dx1**4*dx2**11*(5.0d0*dx2+12.0d0*dx1)
      ci(3)=2.0d0*dx1**3*dx2**10*(5.0d0*dx2**2+3.0d0*dx1*(10.0d0*dx2+11.0d0*dx1))
      ci(4)=10.0d0*dx1**2*dx2**9*(dx2**3+dx1*(12.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**8*(dx2**4+dx1*(24.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+dx1*(20.0d0*dx2+9.0d0*dx1))))
      ci(6)=dx2**7*(dx2**5+dx1*(60.0d0*dx2**4+11.0d0*dx1*(60.0d0*dx2**3+dx1*(200.0d0*dx2**2+9.0d0*dx1*(25.0d0*dx2+8.0d0*dx1)))))
      ci(7)=2.0d0*dx2**6*(6.0d0*dx2**5+11.0d0*dx1*(15.0d0*dx2**4+dx1*(100.0d0*dx2**3+3.0d0*dx1*(75.0d0*dx2**2+2.0d0*dx1* &
(30.0d0*dx2+7.0d0*dx1)))))
      ci(8)=22.0d0*dx2**5*(3.0d0*dx2**5+dx1*(50.0d0*dx2**4+3.0d0*dx1*(75.0d0*dx2**3+2.0d0*dx1*(60.0d0*dx2**2+dx1*(35.0d0*dx2 &
+6.0d0*dx1)))))
      ci(9)=55.0d0*dx2**4*(4.0d0*dx2**5+3.0d0*dx1*(15.0d0*dx2**4+dx1*(48.0d0*dx2**3+dx1*(56.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+dx1))) &
))
      ci(10)=55.0d0*dx2**3*(9.0d0*dx2**5+dx1*(72.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(144.0d0*dx2**2+dx1*(45.0d0*dx2+4.0d0*dx1)))))
      ci(11)=22.0d0*dx2**2*(36.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(360.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1*(50.0d0*dx2+3.0d0*dx1))) &
))
      ci(12)=2.0d0*dx2*(462.0d0*dx2**5+dx1*(1980.0d0*dx2**4+dx1*(2475.0d0*dx2**3+dx1*(1100.0d0*dx2**2+3.0d0*dx1*(55.0d0*dx2 &
+2.0d0*dx1)))))
      ci(13)=792.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(2200.0d0*dx2**3+dx1*(660.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))
      ci(14)=5.0d0*(99.0d0*dx2**4+dx1*(220.0d0*dx2**3+dx1*(132.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(15)=10.0d0*(22.0d0*dx2**3+dx1*(33.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(16)=2.0d0*(33.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+dx1))
      ci(17)=12.0d0*dx2+5.0d0*dx1
      ci(18)=1.0d0
    case(13)
      ci(1)=dx1**5*dx2**13
      ci(2)=dx1**4*dx2**12*(5.0d0*dx2+13.0d0*dx1)
      ci(3)=dx1**3*dx2**11*(10.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+6.0d0*dx1))
      ci(4)=2.0d0*dx1**2*dx2**10*(5.0d0*dx2**3+13.0d0*dx1*(5.0d0*dx2**2+dx1*(15.0d0*dx2+11.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**9*(dx2**4+13.0d0*dx1*(2.0d0*dx2**3+dx1*(12.0d0*dx2**2+11.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(6)=dx2**8*(dx2**5+13.0d0*dx1*(5.0d0*dx2**4+dx1*(60.0d0*dx2**3+11.0d0*dx1*(20.0d0*dx2**2+dx1*(25.0d0*dx2+9.0d0*dx1)))))
      ci(7)=13.0d0*dx2**7*(dx2**5+dx1*(30.0d0*dx2**4+11.0d0*dx1*(20.0d0*dx2**3+dx1*(50.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+4.0d0*dx1) &
))))
      ci(8)=26.0d0*dx2**6*(3.0d0*dx2**5+11.0d0*dx1*(5.0d0*dx2**4+dx1*(25.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2 &
+dx1)))))
      ci(9)=143.0d0*dx2**5*(2.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(30.0d0*dx2**3+dx1*(40.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1) &
))))
      ci(10)=715.0d0*dx2**4*(dx2**5+dx1*(9.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(24.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))))
      ci(11)=143.0d0*dx2**3*(9.0d0*dx2**5+dx1*(60.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(25.0d0*dx2+2.0d0*dx1)))))
      ci(12)=26.0d0*dx2**2*(66.0d0*dx2**5+dx1*(330.0d0*dx2**4+dx1*(495.0d0*dx2**3+dx1*(275.0d0*dx2**2+dx1*(55.0d0*dx2+3.0d0*dx1))) &
))
      ci(13)=13.0d0*dx2*(132.0d0*dx2**5+dx1*(495.0d0*dx2**4+dx1*(550.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))
      ci(14)=1287.0d0*dx2**5+dx1*(3575.0d0*dx2**4+dx1*(2860.0d0*dx2**3+dx1*(780.0d0*dx2**2+dx1*(65.0d0*dx2+dx1))))
      ci(15)=5.0d0*(143.0d0*dx2**4+dx1*(286.0d0*dx2**3+dx1*(156.0d0*dx2**2+dx1*(26.0d0*dx2+dx1))))
      ci(16)=2.0d0*(143.0d0*dx2**3+5.0d0*dx1*(39.0d0*dx2**2+dx1*(13.0d0*dx2+dx1)))
      ci(17)=78.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+2.0d0*dx1)
      ci(18)=13.0d0*dx2+5.0d0*dx1
      ci(19)=1.0d0
    case(14)
      ci(1)=dx1**5*dx2**14
      ci(2)=dx1**4*dx2**13*(5.0d0*dx2+14.0d0*dx1)
      ci(3)=dx1**3*dx2**12*(10.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+13.0d0*dx1))
      ci(4)=dx1**2*dx2**11*(10.0d0*dx2**3+7.0d0*dx1*(20.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+4.0d0*dx1)))
      ci(5)=dx1*dx2**10*(5.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+13.0d0*dx1*(10.0d0*dx2**2+dx1*(20.0d0*dx2+11.0d0*dx1))))
      ci(6)=dx2**9*(dx2**5+7.0d0*dx1*(10.0d0*dx2**4+13.0d0*dx1*(10.0d0*dx2**3+dx1*(40.0d0*dx2**2+11.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)) &
)))
      ci(7)=7.0d0*dx2**8*(2.0d0*dx2**5+13.0d0*dx1*(5.0d0*dx2**4+dx1*(40.0d0*dx2**3+11.0d0*dx1*(10.0d0*dx2**2+dx1*(10.0d0*dx2 &
+3.0d0*dx1)))))
      ci(8)=13.0d0*dx2**7*(7.0d0*dx2**5+dx1*(140.0d0*dx2**4+11.0d0*dx1*(70.0d0*dx2**3+dx1*(140.0d0*dx2**2+3.0d0*dx1*(35.0d0*dx2 &
+8.0d0*dx1)))))
      ci(9)=13.0d0*dx2**6*(28.0d0*dx2**5+11.0d0*dx1*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+3.0d0*dx1*(70.0d0*dx2**2+dx1*(40.0d0*dx2 &
+7.0d0*dx1)))))
      ci(10)=143.0d0*dx2**5*(7.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(210.0d0*dx2**3+dx1*(240.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2 &
+2.0d0*dx1)))))
      ci(11)=143.0d0*dx2**4*(14.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(240.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))) &
))
      ci(12)=13.0d0*dx2**3*(231.0d0*dx2**5+dx1*(1320.0d0*dx2**4+7.0d0*dx1*(330.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(55.0d0*dx2 &
+4.0d0*dx1)))))
      ci(13)=13.0d0*dx2**2*(264.0d0*dx2**5+7.0d0*dx1*(165.0d0*dx2**4+dx1*(220.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)) &
)))
      ci(14)=7.0d0*dx2*(429.0d0*dx2**5+dx1*(1430.0d0*dx2**4+dx1*(1430.0d0*dx2**3+dx1*(520.0d0*dx2**2+dx1*(65.0d0*dx2+2.0d0*dx1)))) &
)
      ci(15)=2002.0d0*dx2**5+dx1*(5005.0d0*dx2**4+dx1*(3640.0d0*dx2**3+dx1*(910.0d0*dx2**2+dx1*(70.0d0*dx2+dx1))))
      ci(16)=1001.0d0*dx2**4+5.0d0*dx1*(364.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(28.0d0*dx2+dx1)))
      ci(17)=364.0d0*dx2**3+5.0d0*dx1*(91.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+dx1))
      ci(18)=91.0d0*dx2**2+10.0d0*dx1*(7.0d0*dx2+dx1)
      ci(19)=14.0d0*dx2+5.0d0*dx1
      ci(20)=1.0d0
    case(15)
      ci(1)=dx1**5*dx2**15
      ci(2)=5.0d0*dx1**4*dx2**14*(dx2+3.0d0*dx1)
      ci(3)=5.0d0*dx1**3*dx2**13*(2.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+7.0d0*dx1))
      ci(4)=5.0d0*dx1**2*dx2**12*(2.0d0*dx2**3+dx1*(30.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+13.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**11*(dx2**4+dx1*(30.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))))
      ci(6)=dx2**10*(dx2**5+dx1*(75.0d0*dx2**4+7.0d0*dx1*(150.0d0*dx2**3+13.0d0*dx1*(50.0d0*dx2**2+3.0d0*dx1*(25.0d0*dx2 &
+11.0d0*dx1)))))
      ci(7)=5.0d0*dx2**9*(3.0d0*dx2**5+7.0d0*dx1*(15.0d0*dx2**4+13.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2 &
+dx1)))))
      ci(8)=5.0d0*dx2**8*(21.0d0*dx2**5+13.0d0*dx1*(35.0d0*dx2**4+dx1*(210.0d0*dx2**3+11.0d0*dx1*(42.0d0*dx2**2+dx1*(35.0d0*dx2 &
+9.0d0*dx1)))))
      ci(9)=65.0d0*dx2**7*(7.0d0*dx2**5+dx1*(105.0d0*dx2**4+11.0d0*dx1*(42.0d0*dx2**3+dx1*(70.0d0*dx2**2+9.0d0*dx1*(5.0d0*dx2+dx1) &
))))
      ci(10)=65.0d0*dx2**6*(21.0d0*dx2**5+11.0d0*dx1*(21.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(45.0d0*dx2 &
+7.0d0*dx1)))))
      ci(11)=143.0d0*dx2**5*(21.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(450.0d0*dx2**3+dx1*(450.0d0*dx2**2+7.0d0*dx1*(25.0d0*dx2 &
+3.0d0*dx1)))))
      ci(12)=65.0d0*dx2**4*(77.0d0*dx2**5+dx1*(495.0d0*dx2**4+dx1*(990.0d0*dx2**3+7.0d0*dx1*(110.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2 &
+dx1)))))
      ci(13)=65.0d0*dx2**3*(99.0d0*dx2**5+dx1*(495.0d0*dx2**4+7.0d0*dx1*(110.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))) &
)
      ci(14)=5.0d0*dx2**2*(1287.0d0*dx2**5+7.0d0*dx1*(715.0d0*dx2**4+dx1*(858.0d0*dx2**3+dx1*(390.0d0*dx2**2+dx1*(65.0d0*dx2 &
+3.0d0*dx1)))))
      ci(15)=5.0d0*dx2*(1001.0d0*dx2**5+dx1*(3003.0d0*dx2**4+dx1*(2730.0d0*dx2**3+dx1*(910.0d0*dx2**2+3.0d0*dx1*(35.0d0*dx2+dx1))) &
))
      ci(16)=3003.0d0*dx2**5+dx1*(6825.0d0*dx2**4+dx1*(4550.0d0*dx2**3+dx1*(1050.0d0*dx2**2+dx1*(75.0d0*dx2+dx1))))
      ci(17)=5.0d0*(273.0d0*dx2**4+dx1*(455.0d0*dx2**3+dx1*(210.0d0*dx2**2+dx1*(30.0d0*dx2+dx1))))
      ci(18)=5.0d0*(91.0d0*dx2**3+dx1*(105.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+dx1)))
      ci(19)=5.0d0*(21.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1))
      ci(20)=5.0d0*(3.0d0*dx2+dx1)
      ci(21)=1.0d0
    case(16)
      ci(1)=dx1**5*dx2**16
      ci(2)=dx1**4*dx2**15*(5.0d0*dx2+16.0d0*dx1)
      ci(3)=10.0d0*dx1**3*dx2**14*(dx2**2+4.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(4)=10.0d0*dx1**2*dx2**13*(dx2**3+4.0d0*dx1*(4.0d0*dx2**2+dx1*(15.0d0*dx2+14.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**12*(dx2**4+4.0d0*dx1*(8.0d0*dx2**3+dx1*(60.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+13.0d0*dx1))))
      ci(6)=dx2**11*(dx2**5+4.0d0*dx1*(20.0d0*dx2**4+dx1*(300.0d0*dx2**3+7.0d0*dx1*(200.0d0*dx2**2+13.0d0*dx1*(25.0d0*dx2 &
+12.0d0*dx1)))))
      ci(7)=8.0d0*dx2**10*(2.0d0*dx2**5+dx1*(75.0d0*dx2**4+7.0d0*dx1*(100.0d0*dx2**3+13.0d0*dx1*(25.0d0*dx2**2+dx1*(30.0d0*dx2 &
+11.0d0*dx1)))))
      ci(8)=40.0d0*dx2**9*(3.0d0*dx2**5+dx1*(70.0d0*dx2**4+13.0d0*dx1*(35.0d0*dx2**3+dx1*(84.0d0*dx2**2+11.0d0*dx1*(7.0d0*dx2 &
+2.0d0*dx1)))))
      ci(9)=10.0d0*dx2**8*(56.0d0*dx2**5+13.0d0*dx1*(70.0d0*dx2**4+dx1*(336.0d0*dx2**3+11.0d0*dx1*(56.0d0*dx2**2+dx1*(40.0d0*dx2 &
+9.0d0*dx1)))))
      ci(10)=130.0d0*dx2**7*(14.0d0*dx2**5+dx1*(168.0d0*dx2**4+11.0d0*dx1*(56.0d0*dx2**3+dx1*(80.0d0*dx2**2+dx1*(45.0d0*dx2 &
+8.0d0*dx1)))))
      ci(11)=52.0d0*dx2**6*(84.0d0*dx2**5+11.0d0*dx1*(70.0d0*dx2**4+dx1*(200.0d0*dx2**3+dx1*(225.0d0*dx2**2+2.0d0*dx1*(50.0d0*dx2 &
+7.0d0*dx1)))))
      ci(12)=52.0d0*dx2**5*(154.0d0*dx2**5+dx1*(1100.0d0*dx2**4+dx1*(2475.0d0*dx2**3+2.0d0*dx1*(1100.0d0*dx2**2+7.0d0*dx1* &
(55.0d0*dx2+6.0d0*dx1)))))
      ci(13)=130.0d0*dx2**4*(88.0d0*dx2**5+dx1*(495.0d0*dx2**4+2.0d0*dx1*(440.0d0*dx2**3+7.0d0*dx1*(44.0d0*dx2**2+dx1*(12.0d0*dx2 &
+dx1)))))
      ci(14)=10.0d0*dx2**3*(1287.0d0*dx2**5+2.0d0*dx1*(2860.0d0*dx2**4+7.0d0*dx1*(572.0d0*dx2**3+dx1*(312.0d0*dx2**2+dx1* &
(65.0d0*dx2+4.0d0*dx1)))))
      ci(15)=40.0d0*dx2**2*(286.0d0*dx2**5+dx1*(1001.0d0*dx2**4+dx1*(1092.0d0*dx2**3+dx1*(455.0d0*dx2**2+dx1*(70.0d0*dx2+3.0d0*dx1 &
)))))
      ci(16)=8.0d0*dx2*(1001.0d0*dx2**5+dx1*(2730.0d0*dx2**4+dx1*(2275.0d0*dx2**3+dx1*(700.0d0*dx2**2+dx1*(75.0d0*dx2+2.0d0*dx1))) &
))
      ci(17)=4368.0d0*dx2**5+dx1*(9100.0d0*dx2**4+dx1*(5600.0d0*dx2**3+dx1*(1200.0d0*dx2**2+dx1*(80.0d0*dx2+dx1))))
      ci(18)=5.0d0*(364.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))
      ci(19)=10.0d0*(56.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(16.0d0*dx2+dx1)))
      ci(20)=10.0d0*(12.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(21)=16.0d0*dx2+5.0d0*dx1
      ci(22)=1.0d0
    case(17)
      ci(1)=dx1**5*dx2**17
      ci(2)=dx1**4*dx2**16*(5.0d0*dx2+17.0d0*dx1)
      ci(3)=dx1**3*dx2**15*(10.0d0*dx2**2+17.0d0*dx1*(5.0d0*dx2+8.0d0*dx1))
      ci(4)=10.0d0*dx1**2*dx2**14*(dx2**3+17.0d0*dx1*(dx2**2+4.0d0*dx1*(dx2+dx1)))
      ci(5)=5.0d0*dx1*dx2**13*(dx2**4+34.0d0*dx1*(dx2**3+2.0d0*dx1*(4.0d0*dx2**2+dx1*(10.0d0*dx2+7.0d0*dx1))))
      ci(6)=dx2**12*(dx2**5+17.0d0*dx1*(5.0d0*dx2**4+4.0d0*dx1*(20.0d0*dx2**3+dx1*(100.0d0*dx2**2+7.0d0*dx1*(25.0d0*dx2+13.0d0*dx1 &
)))))
      ci(7)=17.0d0*dx2**11*(dx2**5+4.0d0*dx1*(10.0d0*dx2**4+dx1*(100.0d0*dx2**3+7.0d0*dx1*(50.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2 &
+2.0d0*dx1)))))
      ci(8)=136.0d0*dx2**10*(dx2**5+dx1*(25.0d0*dx2**4+dx1*(175.0d0*dx2**3+13.0d0*dx1*(35.0d0*dx2**2+dx1*(35.0d0*dx2+11.0d0*dx1))) &
))
      ci(9)=170.0d0*dx2**9*(4.0d0*dx2**5+dx1*(70.0d0*dx2**4+13.0d0*dx1*(28.0d0*dx2**3+dx1*(56.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+dx1 &
)))))
      ci(10)=170.0d0*dx2**8*(14.0d0*dx2**5+13.0d0*dx1*(14.0d0*dx2**4+dx1*(56.0d0*dx2**3+11.0d0*dx1*(8.0d0*dx2**2+dx1*(5.0d0*dx2 &
+dx1)))))
      ci(11)=442.0d0*dx2**7*(14.0d0*dx2**5+dx1*(140.0d0*dx2**4+11.0d0*dx1*(40.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(25.0d0*dx2 &
+4.0d0*dx1)))))
      ci(12)=884.0d0*dx2**6*(14.0d0*dx2**5+dx1*(110.0d0*dx2**4+dx1*(275.0d0*dx2**3+dx1*(275.0d0*dx2**2+2.0d0*dx1*(55.0d0*dx2 &
+7.0d0*dx1)))))
      ci(13)=442.0d0*dx2**5*(44.0d0*dx2**5+dx1*(275.0d0*dx2**4+2.0d0*dx1*(275.0d0*dx2**3+dx1*(220.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2 &
+dx1)))))
      ci(14)=170.0d0*dx2**4*(143.0d0*dx2**5+dx1*(715.0d0*dx2**4+2.0d0*dx1*(572.0d0*dx2**3+7.0d0*dx1*(52.0d0*dx2**2+dx1*(13.0d0*dx2 &
+dx1)))))
      ci(15)=170.0d0*dx2**3*(143.0d0*dx2**5+2.0d0*dx1*(286.0d0*dx2**4+dx1*(364.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(35.0d0*dx2 &
+2.0d0*dx1)))))
      ci(16)=136.0d0*dx2**2*(143.0d0*dx2**5+dx1*(455.0d0*dx2**4+dx1*(455.0d0*dx2**3+dx1*(175.0d0*dx2**2+dx1*(25.0d0*dx2+dx1)))))
      ci(17)=17.0d0*dx2*(728.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1400.0d0*dx2**3+dx1*(400.0d0*dx2**2+dx1*(40.0d0*dx2+dx1)))))
      ci(18)=6188.0d0*dx2**5+dx1*(11900.0d0*dx2**4+dx1*(6800.0d0*dx2**3+dx1*(1360.0d0*dx2**2+dx1*(85.0d0*dx2+dx1))))
      ci(19)=5.0d0*(476.0d0*dx2**4+dx1*(680.0d0*dx2**3+dx1*(272.0d0*dx2**2+dx1*(34.0d0*dx2+dx1))))
      ci(20)=10.0d0*(68.0d0*dx2**3+dx1*(68.0d0*dx2**2+dx1*(17.0d0*dx2+dx1)))
      ci(21)=136.0d0*dx2**2+5.0d0*dx1*(17.0d0*dx2+2.0d0*dx1)
      ci(22)=17.0d0*dx2+5.0d0*dx1
      ci(23)=1.0d0
    case(18)
      ci(1)=dx1**5*dx2**18
      ci(2)=dx1**4*dx2**17*(5.0d0*dx2+18.0d0*dx1)
      ci(3)=dx1**3*dx2**16*(10.0d0*dx2**2+9.0d0*dx1*(10.0d0*dx2+17.0d0*dx1))
      ci(4)=dx1**2*dx2**15*(10.0d0*dx2**3+3.0d0*dx1*(60.0d0*dx2**2+17.0d0*dx1*(15.0d0*dx2+16.0d0*dx1)))
      ci(5)=5.0d0*dx1*dx2**14*(dx2**4+6.0d0*dx1*(6.0d0*dx2**3+17.0d0*dx1*(3.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))))
      ci(6)=dx2**13*(dx2**5+6.0d0*dx1*(15.0d0*dx2**4+17.0d0*dx1*(15.0d0*dx2**3+2.0d0*dx1*(40.0d0*dx2**2+3.0d0*dx1*(25.0d0*dx2 &
+14.0d0*dx1)))))
      ci(7)=3.0d0*dx2**12*(6.0d0*dx2**5+17.0d0*dx1*(15.0d0*dx2**4+4.0d0*dx1*(40.0d0*dx2**3+dx1*(150.0d0*dx2**2+7.0d0*dx1* &
(30.0d0*dx2+13.0d0*dx1)))))
      ci(8)=51.0d0*dx2**11*(3.0d0*dx2**5+4.0d0*dx1*(20.0d0*dx2**4+dx1*(150.0d0*dx2**3+dx1*(420.0d0*dx2**2+13.0d0*dx1*(35.0d0*dx2 &
+12.0d0*dx1)))))
      ci(9)=102.0d0*dx2**10*(8.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(840.0d0*dx2**3+13.0d0*dx1*(140.0d0*dx2**2+3.0d0*dx1*(40.0d0*dx2 &
+11.0d0*dx1)))))
      ci(10)=170.0d0*dx2**9*(18.0d0*dx2**5+dx1*(252.0d0*dx2**4+13.0d0*dx1*(84.0d0*dx2**3+dx1*(144.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2 &
+2.0d0*dx1)))))
      ci(11)=34.0d0*dx2**8*(252.0d0*dx2**5+13.0d0*dx1*(210.0d0*dx2**4+dx1*(720.0d0*dx2**3+11.0d0*dx1*(90.0d0*dx2**2+dx1* &
(50.0d0*dx2+9.0d0*dx1)))))
      ci(12)=442.0d0*dx2**7*(42.0d0*dx2**5+dx1*(360.0d0*dx2**4+dx1*(990.0d0*dx2**3+dx1*(1100.0d0*dx2**2+9.0d0*dx1*(55.0d0*dx2 &
+8.0d0*dx1)))))
      ci(13)=442.0d0*dx2**6*(72.0d0*dx2**5+dx1*(495.0d0*dx2**4+2.0d0*dx1*(550.0d0*dx2**3+3.0d0*dx1*(165.0d0*dx2**2+dx1*(60.0d0*dx2 &
+7.0d0*dx1)))))
      ci(14)=34.0d0*dx2**5*(1287.0d0*dx2**5+2.0d0*dx1*(3575.0d0*dx2**4+3.0d0*dx1*(2145.0d0*dx2**3+dx1*(1560.0d0*dx2**2+7.0d0*dx1* &
(65.0d0*dx2+6.0d0*dx1)))))
      ci(15)=170.0d0*dx2**4*(286.0d0*dx2**5+3.0d0*dx1*(429.0d0*dx2**4+2.0d0*dx1*(312.0d0*dx2**3+dx1*(182.0d0*dx2**2+3.0d0*dx1* &
(14.0d0*dx2+dx1)))))
      ci(16)=102.0d0*dx2**3*(429.0d0*dx2**5+2.0d0*dx1*(780.0d0*dx2**4+dx1*(910.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(75.0d0*dx2 &
+4.0d0*dx1)))))
      ci(17)=51.0d0*dx2**2*(624.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1680.0d0*dx2**3+dx1*(600.0d0*dx2**2+dx1*(80.0d0*dx2+3.0d0*dx1 &
)))))
      ci(18)=3.0d0*dx2*(6188.0d0*dx2**5+dx1*(14280.0d0*dx2**4+dx1*(10200.0d0*dx2**3+dx1*(2720.0d0*dx2**2+3.0d0*dx1*(85.0d0*dx2 &
+2.0d0*dx1)))))
      ci(19)=8568.0d0*dx2**5+dx1*(15300.0d0*dx2**4+dx1*(8160.0d0*dx2**3+dx1*(1530.0d0*dx2**2+dx1*(90.0d0*dx2+dx1))))
      ci(20)=5.0d0*(612.0d0*dx2**4+dx1*(816.0d0*dx2**3+dx1*(306.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))
      ci(21)=816.0d0*dx2**3+5.0d0*dx1*(153.0d0*dx2**2+2.0d0*dx1*(18.0d0*dx2+dx1))
      ci(22)=153.0d0*dx2**2+10.0d0*dx1*(9.0d0*dx2+dx1)
      ci(23)=18.0d0*dx2+5.0d0*dx1
      ci(24)=1.0d0
    case(19)
      ci(1)=dx1**5*dx2**19
      ci(2)=dx1**4*dx2**18*(5.0d0*dx2+19.0d0*dx1)
      ci(3)=dx1**3*dx2**17*(10.0d0*dx2**2+19.0d0*dx1*(5.0d0*dx2+9.0d0*dx1))
      ci(4)=dx1**2*dx2**16*(10.0d0*dx2**3+19.0d0*dx1*(10.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+17.0d0*dx1)))
      ci(5)=dx1*dx2**15*(5.0d0*dx2**4+19.0d0*dx1*(10.0d0*dx2**3+3.0d0*dx1*(30.0d0*dx2**2+17.0d0*dx1*(5.0d0*dx2+4.0d0*dx1))))
      ci(6)=dx2**14*(dx2**5+19.0d0*dx1*(5.0d0*dx2**4+6.0d0*dx1*(15.0d0*dx2**3+17.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=19.0d0*dx2**13*(dx2**5+3.0d0*dx1*(15.0d0*dx2**4+34.0d0*dx1*(5.0d0*dx2**3+2.0d0*dx1*(10.0d0*dx2**2+dx1*(15.0d0*dx2 &
+7.0d0*dx1)))))
      ci(8)=57.0d0*dx2**12*(3.0d0*dx2**5+17.0d0*dx1*(5.0d0*dx2**4+4.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(35.0d0*dx2 &
+13.0d0*dx1)))))
      ci(9)=969.0d0*dx2**11*(dx2**5+2.0d0*dx1*(10.0d0*dx2**4+dx1*(60.0d0*dx2**3+dx1*(140.0d0*dx2**2+13.0d0*dx1*(10.0d0*dx2 &
+3.0d0*dx1)))))
      ci(10)=646.0d0*dx2**10*(6.0d0*dx2**5+dx1*(90.0d0*dx2**4+dx1*(420.0d0*dx2**3+13.0d0*dx1*(60.0d0*dx2**2+dx1*(45.0d0*dx2 &
+11.0d0*dx1)))))
      ci(11)=646.0d0*dx2**9*(18.0d0*dx2**5+dx1*(210.0d0*dx2**4+13.0d0*dx1*(60.0d0*dx2**3+dx1*(90.0d0*dx2**2+11.0d0*dx1*(5.0d0*dx2 &
+dx1)))))
      ci(12)=646.0d0*dx2**8*(42.0d0*dx2**5+13.0d0*dx1*(30.0d0*dx2**4+dx1*(90.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(55.0d0*dx2 &
+9.0d0*dx1)))))
      ci(13)=8398.0d0*dx2**7*(6.0d0*dx2**5+dx1*(45.0d0*dx2**4+dx1*(110.0d0*dx2**3+dx1*(110.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2 &
+2.0d0*dx1)))))
      ci(14)=646.0d0*dx2**6*(117.0d0*dx2**5+dx1*(715.0d0*dx2**4+2.0d0*dx1*(715.0d0*dx2**3+3.0d0*dx1*(195.0d0*dx2**2+dx1* &
(65.0d0*dx2+7.0d0*dx1)))))
      ci(15)=646.0d0*dx2**5*(143.0d0*dx2**5+dx1*(715.0d0*dx2**4+6.0d0*dx1*(195.0d0*dx2**3+dx1*(130.0d0*dx2**2+dx1*(35.0d0*dx2 &
+3.0d0*dx1)))))
      ci(16)=646.0d0*dx2**4*(143.0d0*dx2**5+3.0d0*dx1*(195.0d0*dx2**4+2.0d0*dx1*(130.0d0*dx2**3+dx1*(70.0d0*dx2**2+dx1*(15.0d0*dx2 &
+dx1)))))
      ci(17)=969.0d0*dx2**3*(78.0d0*dx2**5+dx1*(260.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)))))
      ci(18)=57.0d0*dx2**2*(884.0d0*dx2**5+dx1*(2380.0d0*dx2**4+dx1*(2040.0d0*dx2**3+dx1*(680.0d0*dx2**2+dx1*(85.0d0*dx2+3.0d0*dx1 &
)))))
      ci(19)=19.0d0*dx2*(1428.0d0*dx2**5+dx1*(3060.0d0*dx2**4+dx1*(2040.0d0*dx2**3+dx1*(510.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))))
      ci(20)=11628.0d0*dx2**5+dx1*(19380.0d0*dx2**4+dx1*(9690.0d0*dx2**3+dx1*(1710.0d0*dx2**2+dx1*(95.0d0*dx2+dx1))))
      ci(21)=3876.0d0*dx2**4+5.0d0*dx1*(969.0d0*dx2**3+dx1*(342.0d0*dx2**2+dx1*(38.0d0*dx2+dx1)))
      ci(22)=969.0d0*dx2**3+5.0d0*dx1*(171.0d0*dx2**2+2.0d0*dx1*(19.0d0*dx2+dx1))
      ci(23)=171.0d0*dx2**2+5.0d0*dx1*(19.0d0*dx2+2.0d0*dx1)
      ci(24)=19.0d0*dx2+5.0d0*dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>19, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^6 * (x-x2)^n2 as sum_k=0^(6+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_6(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**6
      ci(2)=6.0d0*dx1**5
      ci(3)=15.0d0*dx1**4
      ci(4)=20.0d0*dx1**3
      ci(5)=15.0d0*dx1**2
      ci(6)=6.0d0*dx1
      ci(7)=1.0d0
    case(1)
      ci(1)=dx1**6*dx2
      ci(2)=dx1**5*(6.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**4*(5.0d0*dx2+2.0d0*dx1)
      ci(4)=5.0d0*dx1**3*(4.0d0*dx2+3.0d0*dx1)
      ci(5)=5.0d0*dx1**2*(3.0d0*dx2+4.0d0*dx1)
      ci(6)=3.0d0*dx1*(2.0d0*dx2+5.0d0*dx1)
      ci(7)=dx2+6.0d0*dx1
      ci(8)=1.0d0
    case(2)
      ci(1)=dx1**6*dx2**2
      ci(2)=2.0d0*dx1**5*dx2*(3.0d0*dx2+dx1)
      ci(3)=dx1**4*(15.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))
      ci(4)=2.0d0*dx1**3*(10.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+dx1))
      ci(5)=5.0d0*dx1**2*(3.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(6)=2.0d0*dx1*(3.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(7)=dx2**2+3.0d0*dx1*(4.0d0*dx2+5.0d0*dx1)
      ci(8)=2.0d0*(dx2+3.0d0*dx1)
      ci(9)=1.0d0
    case(3)
      ci(1)=dx1**6*dx2**3
      ci(2)=3.0d0*dx1**5*dx2**2*(2.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**4*dx2*(5.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(4)=dx1**3*(20.0d0*dx2**3+dx1*(45.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))
      ci(5)=3.0d0*dx1**2*(5.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(6)=3.0d0*dx1*(2.0d0*dx2**3+5.0d0*dx1*(3.0d0*dx2**2+dx1*(4.0d0*dx2+dx1)))
      ci(7)=dx2**3+dx1*(18.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+4.0d0*dx1))
      ci(8)=3.0d0*(dx2**2+dx1*(6.0d0*dx2+5.0d0*dx1))
      ci(9)=3.0d0*(dx2+2.0d0*dx1)
      ci(10)=1.0d0
    case(4)
      ci(1)=dx1**6*dx2**4
      ci(2)=2.0d0*dx1**5*dx2**3*(3.0d0*dx2+2.0d0*dx1)
      ci(3)=3.0d0*dx1**4*dx2**2*(5.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+dx1))
      ci(4)=4.0d0*dx1**3*dx2*(5.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(5)=dx1**2*(15.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(6)=6.0d0*dx1*(dx2**4+dx1*(10.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(7)=dx2**4+dx1*(24.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1)))
      ci(8)=4.0d0*(dx2**3+dx1*(9.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(9)=3.0d0*(2.0d0*dx2**2+dx1*(8.0d0*dx2+5.0d0*dx1))
      ci(10)=2.0d0*(2.0d0*dx2+3.0d0*dx1)
      ci(11)=1.0d0
    case(5)
      ci(1)=dx1**6*dx2**5
      ci(2)=dx1**5*dx2**4*(6.0d0*dx2+5.0d0*dx1)
      ci(3)=5.0d0*dx1**4*dx2**3*(3.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))
      ci(4)=5.0d0*dx1**3*dx2**2*(4.0d0*dx2**3+dx1*(15.0d0*dx2**2+2.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**2*dx2*(3.0d0*dx2**4+dx1*(20.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(6)=dx1*(6.0d0*dx2**5+dx1*(75.0d0*dx2**4+dx1*(200.0d0*dx2**3+dx1*(150.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))
      ci(7)=dx2**5+dx1*(30.0d0*dx2**4+dx1*(150.0d0*dx2**3+dx1*(200.0d0*dx2**2+3.0d0*dx1*(25.0d0*dx2+2.0d0*dx1))))
      ci(8)=5.0d0*(dx2**4+dx1*(12.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1))))
      ci(9)=5.0d0*(2.0d0*dx2**3+dx1*(12.0d0*dx2**2+dx1*(15.0d0*dx2+4.0d0*dx1)))
      ci(10)=5.0d0*(2.0d0*dx2**2+3.0d0*dx1*(2.0d0*dx2+dx1))
      ci(11)=5.0d0*dx2+6.0d0*dx1
      ci(12)=1.0d0
    case(6)
      ci(1)=dx1**6*dx2**6
      ci(2)=6.0d0*dx1**5*dx2**5*(dx2+dx1)
      ci(3)=3.0d0*dx1**4*dx2**4*(5.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(4)=10.0d0*dx1**3*dx2**3*(2.0d0*dx2**3+dx1*(9.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(5)=15.0d0*dx1**2*dx2**2*(dx2**4+dx1*(8.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(6)=6.0d0*dx1*dx2*(dx2**5+dx1*(15.0d0*dx2**4+dx1*(50.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(7)=dx2**6+dx1*(36.0d0*dx2**5+dx1*(225.0d0*dx2**4+dx1*(400.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1*(36.0d0*dx2+dx1)))))
      ci(8)=6.0d0*(dx2**5+dx1*(15.0d0*dx2**4+dx1*(50.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(9)=15.0d0*(dx2**4+dx1*(8.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(10)=10.0d0*(2.0d0*dx2**3+dx1*(9.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(11)=3.0d0*(5.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(12)=6.0d0*(dx2+dx1)
      ci(13)=1.0d0
    case(7)
      ci(1)=dx1**6*dx2**7
      ci(2)=dx1**5*dx2**6*(6.0d0*dx2+7.0d0*dx1)
      ci(3)=3.0d0*dx1**4*dx2**5*(5.0d0*dx2**2+7.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=dx1**3*dx2**4*(20.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1)))
      ci(5)=5.0d0*dx1**2*dx2**3*(3.0d0*dx2**4+7.0d0*dx1*(4.0d0*dx2**3+dx1*(9.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))
      ci(6)=3.0d0*dx1*dx2**2*(2.0d0*dx2**5+7.0d0*dx1*(5.0d0*dx2**4+dx1*(20.0d0*dx2**3+dx1*(25.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(7)=dx2*(dx2**6+7.0d0*dx1*(6.0d0*dx2**5+dx1*(45.0d0*dx2**4+dx1*(100.0d0*dx2**3+dx1*(75.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))))) &
)
      ci(8)=7.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(525.0d0*dx2**4+dx1*(700.0d0*dx2**3+dx1*(315.0d0*dx2**2+dx1*(42.0d0*dx2+dx1)))))
      ci(9)=3.0d0*(7.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(175.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1)))))
      ci(10)=5.0d0*(7.0d0*dx2**4+dx1*(42.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1))))
      ci(11)=35.0d0*dx2**3+dx1*(126.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+4.0d0*dx1))
      ci(12)=3.0d0*(7.0d0*dx2**2+dx1*(14.0d0*dx2+5.0d0*dx1))
      ci(13)=7.0d0*dx2+6.0d0*dx1
      ci(14)=1.0d0
    case(8)
      ci(1)=dx1**6*dx2**8
      ci(2)=2.0d0*dx1**5*dx2**7*(3.0d0*dx2+4.0d0*dx1)
      ci(3)=dx1**4*dx2**6*(15.0d0*dx2**2+4.0d0*dx1*(12.0d0*dx2+7.0d0*dx1))
      ci(4)=4.0d0*dx1**3*dx2**5*(5.0d0*dx2**3+2.0d0*dx1*(15.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(5)=dx1**2*dx2**4*(15.0d0*dx2**4+2.0d0*dx1*(80.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1))))
      ci(6)=2.0d0*dx1*dx2**3*(3.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=dx2**2*(dx2**6+2.0d0*dx1*(24.0d0*dx2**5+7.0d0*dx1*(30.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(75.0d0*dx2**2+2.0d0*dx1* &
(12.0d0*dx2+dx1))))))
      ci(8)=8.0d0*dx2*(dx2**6+dx1*(21.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(175.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)) &
))))
      ci(9)=28.0d0*dx2**6+dx1*(336.0d0*dx2**5+dx1*(1050.0d0*dx2**4+dx1*(1120.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(48.0d0*dx2+dx1))) &
))
      ci(10)=2.0d0*(28.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(280.0d0*dx2**2+3.0d0*dx1*(20.0d0*dx2+dx1)))))
      ci(11)=70.0d0*dx2**4+dx1*(336.0d0*dx2**3+5.0d0*dx1*(84.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1)))
      ci(12)=4.0d0*(14.0d0*dx2**3+dx1*(42.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(13)=28.0d0*dx2**2+3.0d0*dx1*(16.0d0*dx2+5.0d0*dx1)
      ci(14)=2.0d0*(4.0d0*dx2+3.0d0*dx1)
      ci(15)=1.0d0
    case(9)
      ci(1)=dx1**6*dx2**9
      ci(2)=3.0d0*dx1**5*dx2**8*(2.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**4*dx2**7*(5.0d0*dx2**2+6.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(4)=dx1**3*dx2**6*(20.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+4.0d0*dx1*(18.0d0*dx2+7.0d0*dx1)))
      ci(5)=3.0d0*dx1**2*dx2**5*(5.0d0*dx2**4+6.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(6)=3.0d0*dx1*dx2**4*(2.0d0*dx2**5+3.0d0*dx1*(15.0d0*dx2**4+2.0d0*dx1*(40.0d0*dx2**3+7.0d0*dx1*(10.0d0*dx2**2+dx1* &
(6.0d0*dx2+dx1)))))
      ci(7)=dx2**3*(dx2**6+6.0d0*dx1*(9.0d0*dx2**5+dx1*(90.0d0*dx2**4+7.0d0*dx1*(40.0d0*dx2**3+dx1*(45.0d0*dx2**2+2.0d0*dx1* &
(9.0d0*dx2+dx1))))))
      ci(8)=9.0d0*dx2**2*(dx2**6+2.0d0*dx1*(12.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(105.0d0*dx2**2+2.0d0*dx1* &
(14.0d0*dx2+dx1))))))
      ci(9)=9.0d0*dx2*(4.0d0*dx2**6+dx1*(56.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(24.0d0*dx2 &
+dx1))))))
      ci(10)=84.0d0*dx2**6+dx1*(756.0d0*dx2**5+dx1*(1890.0d0*dx2**4+dx1*(1680.0d0*dx2**3+dx1*(540.0d0*dx2**2+dx1*(54.0d0*dx2+dx1)) &
)))
      ci(11)=3.0d0*(42.0d0*dx2**5+dx1*(252.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))))
      ci(12)=3.0d0*(42.0d0*dx2**4+dx1*(168.0d0*dx2**3+5.0d0*dx1*(36.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(13)=84.0d0*dx2**3+dx1*(216.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+4.0d0*dx1))
      ci(14)=3.0d0*(12.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1))
      ci(15)=3.0d0*(3.0d0*dx2+2.0d0*dx1)
      ci(16)=1.0d0
    case(10)
      ci(1)=dx1**6*dx2**10
      ci(2)=2.0d0*dx1**5*dx2**9*(3.0d0*dx2+5.0d0*dx1)
      ci(3)=15.0d0*dx1**4*dx2**8*(dx2**2+dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(4)=10.0d0*dx1**3*dx2**7*(2.0d0*dx2**3+3.0d0*dx1*(5.0d0*dx2**2+dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(5)=5.0d0*dx1**2*dx2**6*(3.0d0*dx2**4+dx1*(40.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+2.0d0*dx1*(24.0d0*dx2+7.0d0*dx1))))
      ci(6)=6.0d0*dx1*dx2**5*(dx2**5+dx1*(25.0d0*dx2**4+6.0d0*dx1*(25.0d0*dx2**3+dx1*(50.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+dx1)))))
      ci(7)=dx2**4*(dx2**6+3.0d0*dx1*(20.0d0*dx2**5+dx1*(225.0d0*dx2**4+2.0d0*dx1*(400.0d0*dx2**3+7.0d0*dx1*(75.0d0*dx2**2+dx1* &
(36.0d0*dx2+5.0d0*dx1))))))
      ci(8)=10.0d0*dx2**3*(dx2**6+3.0d0*dx1*(9.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1* &
(21.0d0*dx2+2.0d0*dx1))))))
      ci(9)=45.0d0*dx2**2*(dx2**6+dx1*(16.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(112.0d0*dx2**3+dx1*(70.0d0*dx2**2+dx1*(16.0d0*dx2+dx1 &
))))))
      ci(10)=10.0d0*dx2*(12.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(378.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1* &
(27.0d0*dx2+dx1))))))
      ci(11)=210.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1*(3150.0d0*dx2**4+dx1*(2400.0d0*dx2**3+dx1*(675.0d0*dx2**2+dx1*(60.0d0*dx2+dx1 &
)))))
      ci(12)=6.0d0*(42.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(300.0d0*dx2**3+dx1*(150.0d0*dx2**2+dx1*(25.0d0*dx2+dx1)))))
      ci(13)=5.0d0*(42.0d0*dx2**4+dx1*(144.0d0*dx2**3+dx1*(135.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))))
      ci(14)=10.0d0*(12.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(15)=15.0d0*(3.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))
      ci(16)=2.0d0*(5.0d0*dx2+3.0d0*dx1)
      ci(17)=1.0d0
    case(11)
      ci(1)=dx1**6*dx2**11
      ci(2)=dx1**5*dx2**10*(6.0d0*dx2+11.0d0*dx1)
      ci(3)=dx1**4*dx2**9*(15.0d0*dx2**2+11.0d0*dx1*(6.0d0*dx2+5.0d0*dx1))
      ci(4)=5.0d0*dx1**3*dx2**8*(4.0d0*dx2**3+33.0d0*dx1*(dx2**2+dx1*(2.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**2*dx2**7*(3.0d0*dx2**4+11.0d0*dx1*(4.0d0*dx2**3+3.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))))
      ci(6)=dx1*dx2**6*(6.0d0*dx2**5+11.0d0*dx1*(15.0d0*dx2**4+dx1*(100.0d0*dx2**3+3.0d0*dx1*(75.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2 &
+7.0d0*dx1)))))
      ci(7)=dx2**5*(dx2**6+33.0d0*dx1*(2.0d0*dx2**5+dx1*(25.0d0*dx2**4+2.0d0*dx1*(50.0d0*dx2**3+dx1*(75.0d0*dx2**2+7.0d0*dx1* &
(6.0d0*dx2+dx1))))))
      ci(8)=11.0d0*dx2**4*(dx2**6+3.0d0*dx1*(10.0d0*dx2**5+dx1*(75.0d0*dx2**4+2.0d0*dx1*(100.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1* &
(42.0d0*dx2+5.0d0*dx1))))))
      ci(9)=55.0d0*dx2**3*(dx2**6+3.0d0*dx1*(6.0d0*dx2**5+dx1*(30.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(12.0d0*dx2 &
+dx1))))))
      ci(10)=55.0d0*dx2**2*(3.0d0*dx2**6+dx1*(36.0d0*dx2**5+dx1*(126.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1* &
(18.0d0*dx2+dx1))))))
      ci(11)=11.0d0*dx2*(30.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1*(600.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1* &
(30.0d0*dx2+dx1))))))
      ci(12)=462.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(4950.0d0*dx2**4+dx1*(3300.0d0*dx2**3+dx1*(825.0d0*dx2**2+dx1*(66.0d0*dx2+dx1 &
)))))
      ci(13)=462.0d0*dx2**5+dx1*(1980.0d0*dx2**4+dx1*(2475.0d0*dx2**3+dx1*(1100.0d0*dx2**2+3.0d0*dx1*(55.0d0*dx2+2.0d0*dx1))))
      ci(14)=5.0d0*(66.0d0*dx2**4+dx1*(198.0d0*dx2**3+dx1*(165.0d0*dx2**2+dx1*(44.0d0*dx2+3.0d0*dx1))))
      ci(15)=5.0d0*(33.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(33.0d0*dx2+4.0d0*dx1)))
      ci(16)=55.0d0*dx2**2+3.0d0*dx1*(22.0d0*dx2+5.0d0*dx1)
      ci(17)=11.0d0*dx2+6.0d0*dx1
      ci(18)=1.0d0
    case(12)
      ci(1)=dx1**6*dx2**12
      ci(2)=6.0d0*dx1**5*dx2**11*(dx2+2.0d0*dx1)
      ci(3)=3.0d0*dx1**4*dx2**10*(5.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+11.0d0*dx1))
      ci(4)=4.0d0*dx1**3*dx2**9*(5.0d0*dx2**3+dx1*(45.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+5.0d0*dx1)))
      ci(5)=15.0d0*dx1**2*dx2**8*(dx2**4+dx1*(16.0d0*dx2**3+11.0d0*dx1*(6.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(6)=6.0d0*dx1*dx2**7*(dx2**5+dx1*(30.0d0*dx2**4+11.0d0*dx1*(20.0d0*dx2**3+dx1*(50.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2 &
+4.0d0*dx1)))))
      ci(7)=dx2**6*(dx2**6+dx1*(72.0d0*dx2**5+11.0d0*dx1*(90.0d0*dx2**4+dx1*(400.0d0*dx2**3+3.0d0*dx1*(225.0d0*dx2**2+4.0d0*dx1* &
(36.0d0*dx2+7.0d0*dx1))))))
      ci(8)=12.0d0*dx2**5*(dx2**6+11.0d0*dx1*(3.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(25.0d0*dx2**3+2.0d0*dx1*(15.0d0*dx2**2 &
+dx1*(7.0d0*dx2+dx1))))))
      ci(9)=33.0d0*dx2**4*(2.0d0*dx2**6+dx1*(40.0d0*dx2**5+3.0d0*dx1*(75.0d0*dx2**4+dx1*(160.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1* &
(48.0d0*dx2+5.0d0*dx1))))))
      ci(10)=110.0d0*dx2**3*(2.0d0*dx2**6+dx1*(27.0d0*dx2**5+dx1*(108.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1* &
(27.0d0*dx2+2.0d0*dx1))))))
      ci(11)=33.0d0*dx2**2*(15.0d0*dx2**6+dx1*(144.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1*(480.0d0*dx2**3+dx1*(225.0d0*dx2**2 &
+2.0d0*dx1*(20.0d0*dx2+dx1))))))
      ci(12)=12.0d0*dx2*(66.0d0*dx2**6+dx1*(462.0d0*dx2**5+dx1*(990.0d0*dx2**4+dx1*(825.0d0*dx2**3+dx1*(275.0d0*dx2**2+dx1* &
(33.0d0*dx2+dx1))))))
      ci(13)=924.0d0*dx2**6+dx1*(4752.0d0*dx2**5+dx1*(7425.0d0*dx2**4+dx1*(4400.0d0*dx2**3+dx1*(990.0d0*dx2**2+dx1*(72.0d0*dx2+dx1 &
)))))
      ci(14)=6.0d0*(132.0d0*dx2**5+dx1*(495.0d0*dx2**4+dx1*(550.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))
      ci(15)=15.0d0*(33.0d0*dx2**4+dx1*(88.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(16)=4.0d0*(55.0d0*dx2**3+dx1*(99.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+dx1)))
      ci(17)=3.0d0*(22.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1))
      ci(18)=6.0d0*(2.0d0*dx2+dx1)
      ci(19)=1.0d0
    case(13)
      ci(1)=dx1**6*dx2**13
      ci(2)=dx1**5*dx2**12*(6.0d0*dx2+13.0d0*dx1)
      ci(3)=3.0d0*dx1**4*dx2**11*(5.0d0*dx2**2+26.0d0*dx1*(dx2+dx1))
      ci(4)=dx1**3*dx2**10*(20.0d0*dx2**3+13.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(18.0d0*dx2+11.0d0*dx1)))
      ci(5)=dx1**2*dx2**9*(15.0d0*dx2**4+13.0d0*dx1*(20.0d0*dx2**3+dx1*(90.0d0*dx2**2+11.0d0*dx1*(12.0d0*dx2+5.0d0*dx1))))
      ci(6)=3.0d0*dx1*dx2**8*(2.0d0*dx2**5+13.0d0*dx1*(5.0d0*dx2**4+dx1*(40.0d0*dx2**3+11.0d0*dx1*(10.0d0*dx2**2+dx1*(10.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=dx2**7*(dx2**6+13.0d0*dx1*(6.0d0*dx2**5+dx1*(90.0d0*dx2**4+11.0d0*dx1*(40.0d0*dx2**3+3.0d0*dx1*(25.0d0*dx2**2 &
+2.0d0*dx1*(9.0d0*dx2+2.0d0*dx1))))))
      ci(8)=13.0d0*dx2**6*(dx2**6+dx1*(36.0d0*dx2**5+11.0d0*dx1*(30.0d0*dx2**4+dx1*(100.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2 &
+4.0d0*dx1*(6.0d0*dx2+dx1))))))
      ci(9)=39.0d0*dx2**5*(2.0d0*dx2**6+11.0d0*dx1*(4.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(20.0d0*dx2**3+dx1*(20.0d0*dx2**2 &
+dx1*(8.0d0*dx2+dx1))))))
      ci(10)=143.0d0*dx2**4*(2.0d0*dx2**6+dx1*(30.0d0*dx2**5+dx1*(135.0d0*dx2**4+dx1*(240.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1* &
(54.0d0*dx2+5.0d0*dx1))))))
      ci(11)=143.0d0*dx2**3*(5.0d0*dx2**6+dx1*(54.0d0*dx2**5+dx1*(180.0d0*dx2**4+dx1*(240.0d0*dx2**3+dx1*(135.0d0*dx2**2 &
+2.0d0*dx1*(15.0d0*dx2+dx1))))))
      ci(12)=39.0d0*dx2**2*(33.0d0*dx2**6+dx1*(264.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(660.0d0*dx2**3+dx1*(275.0d0*dx2**2 &
+2.0d0*dx1*(22.0d0*dx2+dx1))))))
      ci(13)=13.0d0*dx2*(132.0d0*dx2**6+dx1*(792.0d0*dx2**5+dx1*(1485.0d0*dx2**4+dx1*(1100.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1* &
(36.0d0*dx2+dx1))))))
      ci(14)=1716.0d0*dx2**6+dx1*(7722.0d0*dx2**5+dx1*(10725.0d0*dx2**4+dx1*(5720.0d0*dx2**3+dx1*(1170.0d0*dx2**2+dx1*(78.0d0*dx2 &
+dx1)))))
      ci(15)=3.0d0*(429.0d0*dx2**5+dx1*(1430.0d0*dx2**4+dx1*(1430.0d0*dx2**3+dx1*(520.0d0*dx2**2+dx1*(65.0d0*dx2+2.0d0*dx1)))))
      ci(16)=715.0d0*dx2**4+dx1*(1716.0d0*dx2**3+5.0d0*dx1*(234.0d0*dx2**2+dx1*(52.0d0*dx2+3.0d0*dx1)))
      ci(17)=286.0d0*dx2**3+dx1*(468.0d0*dx2**2+5.0d0*dx1*(39.0d0*dx2+4.0d0*dx1))
      ci(18)=3.0d0*(26.0d0*dx2**2+dx1*(26.0d0*dx2+5.0d0*dx1))
      ci(19)=13.0d0*dx2+6.0d0*dx1
      ci(20)=1.0d0
    case(14)
      ci(1)=dx1**6*dx2**14
      ci(2)=2.0d0*dx1**5*dx2**13*(3.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**4*dx2**12*(15.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+13.0d0*dx1))
      ci(4)=2.0d0*dx1**3*dx2**11*(10.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**2*dx2**10*(15.0d0*dx2**4+7.0d0*dx1*(40.0d0*dx2**3+13.0d0*dx1*(15.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1))))
      ci(6)=2.0d0*dx1*dx2**9*(3.0d0*dx2**5+7.0d0*dx1*(15.0d0*dx2**4+13.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+11.0d0*dx1* &
(3.0d0*dx2+dx1)))))
      ci(7)=dx2**8*(dx2**6+7.0d0*dx1*(12.0d0*dx2**5+13.0d0*dx1*(15.0d0*dx2**4+dx1*(80.0d0*dx2**3+33.0d0*dx1*(5.0d0*dx2**2+dx1* &
(4.0d0*dx2+dx1))))))
      ci(8)=2.0d0*dx2**7*(7.0d0*dx2**6+13.0d0*dx1*(21.0d0*dx2**5+dx1*(210.0d0*dx2**4+11.0d0*dx1*(70.0d0*dx2**3+3.0d0*dx1* &
(35.0d0*dx2**2+dx1*(21.0d0*dx2+4.0d0*dx1))))))
      ci(9)=13.0d0*dx2**6*(7.0d0*dx2**6+dx1*(168.0d0*dx2**5+11.0d0*dx1*(105.0d0*dx2**4+dx1*(280.0d0*dx2**3+3.0d0*dx1* &
(105.0d0*dx2**2+dx1*(48.0d0*dx2+7.0d0*dx1))))))
      ci(10)=26.0d0*dx2**5*(14.0d0*dx2**6+11.0d0*dx1*(21.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(210.0d0*dx2**3+dx1*(180.0d0*dx2**2 &
+7.0d0*dx1*(9.0d0*dx2+dx1))))))
      ci(11)=143.0d0*dx2**4*(7.0d0*dx2**6+dx1*(84.0d0*dx2**5+dx1*(315.0d0*dx2**4+dx1*(480.0d0*dx2**3+7.0d0*dx1*(45.0d0*dx2**2+dx1* &
(12.0d0*dx2+dx1))))))
      ci(12)=26.0d0*dx2**3*(77.0d0*dx2**6+dx1*(693.0d0*dx2**5+dx1*(1980.0d0*dx2**4+7.0d0*dx1*(330.0d0*dx2**3+dx1*(165.0d0*dx2**2 &
+dx1*(33.0d0*dx2+2.0d0*dx1))))))
      ci(13)=13.0d0*dx2**2*(231.0d0*dx2**6+dx1*(1584.0d0*dx2**5+7.0d0*dx1*(495.0d0*dx2**4+dx1*(440.0d0*dx2**3+dx1*(165.0d0*dx2**2 &
+dx1*(24.0d0*dx2+dx1))))))
      ci(14)=2.0d0*dx2*(1716.0d0*dx2**6+7.0d0*dx1*(1287.0d0*dx2**5+dx1*(2145.0d0*dx2**4+dx1*(1430.0d0*dx2**3+dx1*(390.0d0*dx2**2 &
+dx1*(39.0d0*dx2+dx1))))))
      ci(15)=3003.0d0*dx2**6+dx1*(12012.0d0*dx2**5+dx1*(15015.0d0*dx2**4+dx1*(7280.0d0*dx2**3+dx1*(1365.0d0*dx2**2+dx1*(84.0d0*dx2 &
+dx1)))))
      ci(16)=2.0d0*(1001.0d0*dx2**5+dx1*(3003.0d0*dx2**4+dx1*(2730.0d0*dx2**3+dx1*(910.0d0*dx2**2+3.0d0*dx1*(35.0d0*dx2+dx1)))))
      ci(17)=1001.0d0*dx2**4+dx1*(2184.0d0*dx2**3+5.0d0*dx1*(273.0d0*dx2**2+dx1*(56.0d0*dx2+3.0d0*dx1)))
      ci(18)=2.0d0*(182.0d0*dx2**3+dx1*(273.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+2.0d0*dx1)))
      ci(19)=91.0d0*dx2**2+3.0d0*dx1*(28.0d0*dx2+5.0d0*dx1)
      ci(20)=2.0d0*(7.0d0*dx2+3.0d0*dx1)
      ci(21)=1.0d0
    case(15)
      ci(1)=dx1**6*dx2**15
      ci(2)=3.0d0*dx1**5*dx2**14*(2.0d0*dx2+5.0d0*dx1)
      ci(3)=15.0d0*dx1**4*dx2**13*(dx2**2+dx1*(6.0d0*dx2+7.0d0*dx1))
      ci(4)=5.0d0*dx1**3*dx2**12*(4.0d0*dx2**3+dx1*(45.0d0*dx2**2+7.0d0*dx1*(18.0d0*dx2+13.0d0*dx1)))
      ci(5)=15.0d0*dx1**2*dx2**11*(dx2**4+dx1*(20.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(6)=3.0d0*dx1*dx2**10*(2.0d0*dx2**5+dx1*(75.0d0*dx2**4+7.0d0*dx1*(100.0d0*dx2**3+13.0d0*dx1*(25.0d0*dx2**2+dx1*(30.0d0*dx2 &
+11.0d0*dx1)))))
      ci(7)=dx2**9*(dx2**6+dx1*(90.0d0*dx2**5+7.0d0*dx1*(225.0d0*dx2**4+13.0d0*dx1*(100.0d0*dx2**3+dx1*(225.0d0*dx2**2+11.0d0*dx1* &
(18.0d0*dx2+5.0d0*dx1))))))
      ci(8)=15.0d0*dx2**8*(dx2**6+dx1*(42.0d0*dx2**5+13.0d0*dx1*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+11.0d0*dx1*(21.0d0*dx2**2+dx1* &
(14.0d0*dx2+3.0d0*dx1))))))
      ci(9)=15.0d0*dx2**7*(7.0d0*dx2**6+13.0d0*dx1*(14.0d0*dx2**5+dx1*(105.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+dx1*(35.0d0*dx2**2 &
+3.0d0*dx1*(6.0d0*dx2+dx1))))))
      ci(10)=65.0d0*dx2**6*(7.0d0*dx2**6+dx1*(126.0d0*dx2**5+11.0d0*dx1*(63.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(135.0d0*dx2**2 &
+dx1*(54.0d0*dx2+7.0d0*dx1))))))
      ci(11)=39.0d0*dx2**5*(35.0d0*dx2**6+11.0d0*dx1*(42.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(300.0d0*dx2**3+dx1*(225.0d0*dx2**2 &
+7.0d0*dx1*(10.0d0*dx2+dx1))))))
      ci(12)=39.0d0*dx2**4*(77.0d0*dx2**6+dx1*(770.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(3300.0d0*dx2**3+7.0d0*dx1*(275.0d0*dx2**2 &
+dx1*(66.0d0*dx2+5.0d0*dx1))))))
      ci(13)=65.0d0*dx2**3*(77.0d0*dx2**6+dx1*(594.0d0*dx2**5+dx1*(1485.0d0*dx2**4+7.0d0*dx1*(220.0d0*dx2**3+dx1*(99.0d0*dx2**2 &
+dx1*(18.0d0*dx2+dx1))))))
      ci(14)=15.0d0*dx2**2*(429.0d0*dx2**6+dx1*(2574.0d0*dx2**5+7.0d0*dx1*(715.0d0*dx2**4+dx1*(572.0d0*dx2**3+dx1*(195.0d0*dx2**2 &
+dx1*(26.0d0*dx2+dx1))))))
      ci(15)=15.0d0*dx2*(429.0d0*dx2**6+dx1*(2002.0d0*dx2**5+dx1*(3003.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1*(455.0d0*dx2**2+dx1* &
(42.0d0*dx2+dx1))))))
      ci(16)=5005.0d0*dx2**6+dx1*(18018.0d0*dx2**5+dx1*(20475.0d0*dx2**4+dx1*(9100.0d0*dx2**3+dx1*(1575.0d0*dx2**2+dx1*(90.0d0*dx2 &
+dx1)))))
      ci(17)=3.0d0*(1001.0d0*dx2**5+dx1*(2730.0d0*dx2**4+dx1*(2275.0d0*dx2**3+dx1*(700.0d0*dx2**2+dx1*(75.0d0*dx2+2.0d0*dx1)))))
      ci(18)=15.0d0*(91.0d0*dx2**4+dx1*(182.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(19)=5.0d0*(91.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(45.0d0*dx2+4.0d0*dx1)))
      ci(20)=15.0d0*(7.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(21)=3.0d0*(5.0d0*dx2+2.0d0*dx1)
      ci(22)=1.0d0
    case(16)
      ci(1)=dx1**6*dx2**16
      ci(2)=2.0d0*dx1**5*dx2**15*(3.0d0*dx2+8.0d0*dx1)
      ci(3)=3.0d0*dx1**4*dx2**14*(5.0d0*dx2**2+8.0d0*dx1*(4.0d0*dx2+5.0d0*dx1))
      ci(4)=20.0d0*dx1**3*dx2**13*(dx2**3+4.0d0*dx1*(3.0d0*dx2**2+dx1*(9.0d0*dx2+7.0d0*dx1)))
      ci(5)=5.0d0*dx1**2*dx2**12*(3.0d0*dx2**4+4.0d0*dx1*(16.0d0*dx2**3+dx1*(90.0d0*dx2**2+7.0d0*dx1*(24.0d0*dx2+13.0d0*dx1))))
      ci(6)=6.0d0*dx1*dx2**11*(dx2**5+4.0d0*dx1*(10.0d0*dx2**4+dx1*(100.0d0*dx2**3+7.0d0*dx1*(50.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=dx2**10*(dx2**6+4.0d0*dx1*(24.0d0*dx2**5+dx1*(450.0d0*dx2**4+7.0d0*dx1*(400.0d0*dx2**3+13.0d0*dx1*(75.0d0*dx2**2 &
+2.0d0*dx1*(36.0d0*dx2+11.0d0*dx1))))))
      ci(8)=16.0d0*dx2**9*(dx2**6+dx1*(45.0d0*dx2**5+dx1*(525.0d0*dx2**4+13.0d0*dx1*(175.0d0*dx2**3+dx1*(315.0d0*dx2**2 &
+11.0d0*dx1*(21.0d0*dx2+5.0d0*dx1))))))
      ci(9)=30.0d0*dx2**8*(4.0d0*dx2**6+dx1*(112.0d0*dx2**5+13.0d0*dx1*(70.0d0*dx2**4+dx1*(224.0d0*dx2**3+11.0d0*dx1* &
(28.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))))))
      ci(10)=20.0d0*dx2**7*(28.0d0*dx2**6+13.0d0*dx1*(42.0d0*dx2**5+dx1*(252.0d0*dx2**4+11.0d0*dx1*(56.0d0*dx2**3+dx1* &
(60.0d0*dx2**2+dx1*(27.0d0*dx2+4.0d0*dx1))))))
      ci(11)=26.0d0*dx2**6*(70.0d0*dx2**6+dx1*(1008.0d0*dx2**5+11.0d0*dx1*(420.0d0*dx2**4+dx1*(800.0d0*dx2**3+dx1*(675.0d0*dx2**2 &
+4.0d0*dx1*(60.0d0*dx2+7.0d0*dx1))))))
      ci(12)=312.0d0*dx2**5*(14.0d0*dx2**6+dx1*(154.0d0*dx2**5+dx1*(550.0d0*dx2**4+dx1*(825.0d0*dx2**3+2.0d0*dx1*(275.0d0*dx2**2 &
+7.0d0*dx1*(11.0d0*dx2+dx1))))))
      ci(13)=26.0d0*dx2**4*(308.0d0*dx2**6+dx1*(2640.0d0*dx2**5+dx1*(7425.0d0*dx2**4+2.0d0*dx1*(4400.0d0*dx2**3+7.0d0*dx1* &
(330.0d0*dx2**2+dx1*(72.0d0*dx2+5.0d0*dx1))))))
      ci(14)=20.0d0*dx2**3*(572.0d0*dx2**6+dx1*(3861.0d0*dx2**5+2.0d0*dx1*(4290.0d0*dx2**4+7.0d0*dx1*(572.0d0*dx2**3+dx1* &
(234.0d0*dx2**2+dx1*(39.0d0*dx2+2.0d0*dx1))))))
      ci(15)=30.0d0*dx2**2*(429.0d0*dx2**6+2.0d0*dx1*(1144.0d0*dx2**5+dx1*(2002.0d0*dx2**4+dx1*(1456.0d0*dx2**3+dx1* &
(455.0d0*dx2**2+2.0d0*dx1*(28.0d0*dx2+dx1))))))
      ci(16)=16.0d0*dx2*(715.0d0*dx2**6+dx1*(3003.0d0*dx2**5+dx1*(4095.0d0*dx2**4+dx1*(2275.0d0*dx2**3+dx1*(525.0d0*dx2**2+dx1* &
(45.0d0*dx2+dx1))))))
      ci(17)=8008.0d0*dx2**6+dx1*(26208.0d0*dx2**5+dx1*(27300.0d0*dx2**4+dx1*(11200.0d0*dx2**3+dx1*(1800.0d0*dx2**2+dx1* &
(96.0d0*dx2+dx1)))))
      ci(18)=6.0d0*(728.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1400.0d0*dx2**3+dx1*(400.0d0*dx2**2+dx1*(40.0d0*dx2+dx1)))))
      ci(19)=5.0d0*(364.0d0*dx2**4+dx1*(672.0d0*dx2**3+dx1*(360.0d0*dx2**2+dx1*(64.0d0*dx2+3.0d0*dx1))))
      ci(20)=20.0d0*(28.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(21)=3.0d0*(40.0d0*dx2**2+dx1*(32.0d0*dx2+5.0d0*dx1))
      ci(22)=2.0d0*(8.0d0*dx2+3.0d0*dx1)
      ci(23)=1.0d0
    case(17)
      ci(1)=dx1**6*dx2**17
      ci(2)=dx1**5*dx2**16*(6.0d0*dx2+17.0d0*dx1)
      ci(3)=dx1**4*dx2**15*(15.0d0*dx2**2+34.0d0*dx1*(3.0d0*dx2+4.0d0*dx1))
      ci(4)=dx1**3*dx2**14*(20.0d0*dx2**3+17.0d0*dx1*(15.0d0*dx2**2+8.0d0*dx1*(6.0d0*dx2+5.0d0*dx1)))
      ci(5)=5.0d0*dx1**2*dx2**13*(3.0d0*dx2**4+68.0d0*dx1*(dx2**3+dx1*(6.0d0*dx2**2+dx1*(12.0d0*dx2+7.0d0*dx1))))
      ci(6)=dx1*dx2**12*(6.0d0*dx2**5+17.0d0*dx1*(15.0d0*dx2**4+4.0d0*dx1*(40.0d0*dx2**3+dx1*(150.0d0*dx2**2+7.0d0*dx1*(30.0d0*dx2 &
+13.0d0*dx1)))))
      ci(7)=dx2**11*(dx2**6+34.0d0*dx1*(3.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+dx1*(200.0d0*dx2**3+7.0d0*dx1*(75.0d0*dx2**2 &
+26.0d0*dx1*(3.0d0*dx2+dx1))))))
      ci(8)=17.0d0*dx2**10*(dx2**6+4.0d0*dx1*(12.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(700.0d0*dx2**3+13.0d0*dx1*(105.0d0*dx2**2 &
+2.0d0*dx1*(42.0d0*dx2+11.0d0*dx1))))))
      ci(9)=34.0d0*dx2**9*(4.0d0*dx2**6+dx1*(120.0d0*dx2**5+dx1*(1050.0d0*dx2**4+13.0d0*dx1*(280.0d0*dx2**3+dx1*(420.0d0*dx2**2 &
+11.0d0*dx1*(24.0d0*dx2+5.0d0*dx1))))))
      ci(10)=170.0d0*dx2**8*(4.0d0*dx2**6+dx1*(84.0d0*dx2**5+13.0d0*dx1*(42.0d0*dx2**4+dx1*(112.0d0*dx2**3+11.0d0*dx1* &
(12.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))))
      ci(11)=34.0d0*dx2**7*(70.0d0*dx2**6+13.0d0*dx1*(84.0d0*dx2**5+dx1*(420.0d0*dx2**4+11.0d0*dx1*(80.0d0*dx2**3+dx1* &
(75.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+2.0d0*dx1))))))
      ci(12)=442.0d0*dx2**6*(14.0d0*dx2**6+dx1*(168.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(1100.0d0*dx2**3+dx1*(825.0d0*dx2**2 &
+4.0d0*dx1*(66.0d0*dx2+7.0d0*dx1))))))
      ci(13)=442.0d0*dx2**5*(28.0d0*dx2**6+dx1*(264.0d0*dx2**5+dx1*(825.0d0*dx2**4+2.0d0*dx1*(550.0d0*dx2**3+dx1*(330.0d0*dx2**2 &
+7.0d0*dx1*(12.0d0*dx2+dx1))))))
      ci(14)=34.0d0*dx2**4*(572.0d0*dx2**6+dx1*(4290.0d0*dx2**5+dx1*(10725.0d0*dx2**4+2.0d0*dx1*(5720.0d0*dx2**3+7.0d0*dx1* &
(390.0d0*dx2**2+dx1*(78.0d0*dx2+5.0d0*dx1))))))
      ci(15)=170.0d0*dx2**3*(143.0d0*dx2**6+2.0d0*dx1*(429.0d0*dx2**5+dx1*(858.0d0*dx2**4+dx1*(728.0d0*dx2**3+dx1*(273.0d0*dx2**2 &
+2.0d0*dx1*(21.0d0*dx2+dx1))))))
      ci(16)=34.0d0*dx2**2*(715.0d0*dx2**6+2.0d0*dx1*(1716.0d0*dx2**5+dx1*(2730.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1* &
(525.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2+dx1))))))
      ci(17)=17.0d0*dx2*(1144.0d0*dx2**6+dx1*(4368.0d0*dx2**5+dx1*(5460.0d0*dx2**4+dx1*(2800.0d0*dx2**3+dx1*(600.0d0*dx2**2+dx1* &
(48.0d0*dx2+dx1))))))
      ci(18)=12376.0d0*dx2**6+dx1*(37128.0d0*dx2**5+dx1*(35700.0d0*dx2**4+dx1*(13600.0d0*dx2**3+dx1*(2040.0d0*dx2**2+dx1* &
(102.0d0*dx2+dx1)))))
      ci(19)=6188.0d0*dx2**5+dx1*(14280.0d0*dx2**4+dx1*(10200.0d0*dx2**3+dx1*(2720.0d0*dx2**2+3.0d0*dx1*(85.0d0*dx2+2.0d0*dx1))))
      ci(20)=5.0d0*(476.0d0*dx2**4+dx1*(816.0d0*dx2**3+dx1*(408.0d0*dx2**2+dx1*(68.0d0*dx2+3.0d0*dx1))))
      ci(21)=680.0d0*dx2**3+dx1*(816.0d0*dx2**2+5.0d0*dx1*(51.0d0*dx2+4.0d0*dx1))
      ci(22)=136.0d0*dx2**2+3.0d0*dx1*(34.0d0*dx2+5.0d0*dx1)
      ci(23)=17.0d0*dx2+6.0d0*dx1
      ci(24)=1.0d0
    case(18)
      ci(1)=dx1**6*dx2**18
      ci(2)=6.0d0*dx1**5*dx2**17*(dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**4*dx2**16*(5.0d0*dx2**2+3.0d0*dx1*(12.0d0*dx2+17.0d0*dx1))
      ci(4)=2.0d0*dx1**3*dx2**15*(10.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+17.0d0*dx1*(9.0d0*dx2+8.0d0*dx1)))
      ci(5)=3.0d0*dx1**2*dx2**14*(5.0d0*dx2**4+3.0d0*dx1*(40.0d0*dx2**3+17.0d0*dx1*(15.0d0*dx2**2+4.0d0*dx1*(8.0d0*dx2+5.0d0*dx1)) &
))
      ci(6)=6.0d0*dx1*dx2**13*(dx2**5+3.0d0*dx1*(15.0d0*dx2**4+34.0d0*dx1*(5.0d0*dx2**3+2.0d0*dx1*(10.0d0*dx2**2+dx1*(15.0d0*dx2 &
+7.0d0*dx1)))))
      ci(7)=dx2**12*(dx2**6+3.0d0*dx1*(36.0d0*dx2**5+17.0d0*dx1*(45.0d0*dx2**4+4.0d0*dx1*(80.0d0*dx2**3+dx1*(225.0d0*dx2**2 &
+7.0d0*dx1*(36.0d0*dx2+13.0d0*dx1))))))
      ci(8)=18.0d0*dx2**11*(dx2**6+17.0d0*dx1*(3.0d0*dx2**5+4.0d0*dx1*(10.0d0*dx2**4+dx1*(50.0d0*dx2**3+dx1*(105.0d0*dx2**2 &
+13.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))))))
      ci(9)=153.0d0*dx2**10*(dx2**6+2.0d0*dx1*(16.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(560.0d0*dx2**3+13.0d0*dx1*(70.0d0*dx2**2 &
+dx1*(48.0d0*dx2+11.0d0*dx1))))))
      ci(10)=68.0d0*dx2**9*(12.0d0*dx2**6+dx1*(270.0d0*dx2**5+dx1*(1890.0d0*dx2**4+13.0d0*dx1*(420.0d0*dx2**3+dx1*(540.0d0*dx2**2 &
+11.0d0*dx1*(27.0d0*dx2+5.0d0*dx1))))))
      ci(11)=102.0d0*dx2**8*(30.0d0*dx2**6+dx1*(504.0d0*dx2**5+13.0d0*dx1*(210.0d0*dx2**4+dx1*(480.0d0*dx2**3+11.0d0*dx1* &
(45.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1))))))
      ci(12)=204.0d0*dx2**7*(42.0d0*dx2**6+13.0d0*dx1*(42.0d0*dx2**5+dx1*(180.0d0*dx2**4+dx1*(330.0d0*dx2**3+dx1*(275.0d0*dx2**2 &
+3.0d0*dx1*(33.0d0*dx2+4.0d0*dx1))))))
      ci(13)=442.0d0*dx2**6*(42.0d0*dx2**6+dx1*(432.0d0*dx2**5+dx1*(1485.0d0*dx2**4+dx1*(2200.0d0*dx2**3+3.0d0*dx1*(495.0d0*dx2**2 &
+2.0d0*dx1*(72.0d0*dx2+7.0d0*dx1))))))
      ci(14)=204.0d0*dx2**5*(156.0d0*dx2**6+dx1*(1287.0d0*dx2**5+dx1*(3575.0d0*dx2**4+6.0d0*dx1*(715.0d0*dx2**3+dx1* &
(390.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+dx1))))))
      ci(15)=102.0d0*dx2**4*(429.0d0*dx2**6+dx1*(2860.0d0*dx2**5+3.0d0*dx1*(2145.0d0*dx2**4+2.0d0*dx1*(1040.0d0*dx2**3+dx1* &
(455.0d0*dx2**2+dx1*(84.0d0*dx2+5.0d0*dx1))))))
      ci(16)=68.0d0*dx2**3*(715.0d0*dx2**6+3.0d0*dx1*(1287.0d0*dx2**5+2.0d0*dx1*(1170.0d0*dx2**4+dx1*(910.0d0*dx2**3+dx1* &
(315.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1))))))
      ci(17)=153.0d0*dx2**2*(286.0d0*dx2**6+dx1*(1248.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1120.0d0*dx2**3+dx1*(300.0d0*dx2**2 &
+dx1*(32.0d0*dx2+dx1))))))
      ci(18)=18.0d0*dx2*(1768.0d0*dx2**6+dx1*(6188.0d0*dx2**5+dx1*(7140.0d0*dx2**4+dx1*(3400.0d0*dx2**3+dx1*(680.0d0*dx2**2+dx1* &
(51.0d0*dx2+dx1))))))
      ci(19)=18564.0d0*dx2**6+dx1*(51408.0d0*dx2**5+dx1*(45900.0d0*dx2**4+dx1*(16320.0d0*dx2**3+dx1*(2295.0d0*dx2**2+dx1* &
(108.0d0*dx2+dx1)))))
      ci(20)=6.0d0*(1428.0d0*dx2**5+dx1*(3060.0d0*dx2**4+dx1*(2040.0d0*dx2**3+dx1*(510.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))))
      ci(21)=3.0d0*(1020.0d0*dx2**4+dx1*(1632.0d0*dx2**3+5.0d0*dx1*(153.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(22)=2.0d0*(408.0d0*dx2**3+dx1*(459.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+2.0d0*dx1)))
      ci(23)=3.0d0*(51.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))
      ci(24)=6.0d0*(3.0d0*dx2+dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>18, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^7 * (x-x2)^n2 as sum_k=0^(7+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_7(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**7
      ci(2)=7.0d0*dx1**6
      ci(3)=21.0d0*dx1**5
      ci(4)=35.0d0*dx1**4
      ci(5)=35.0d0*dx1**3
      ci(6)=21.0d0*dx1**2
      ci(7)=7.0d0*dx1
      ci(8)=1.0d0
    case(1)
      ci(1)=dx1**7*dx2
      ci(2)=dx1**6*(7.0d0*dx2+dx1)
      ci(3)=7.0d0*dx1**5*(3.0d0*dx2+dx1)
      ci(4)=7.0d0*dx1**4*(5.0d0*dx2+3.0d0*dx1)
      ci(5)=35.0d0*dx1**3*(dx2+dx1)
      ci(6)=7.0d0*dx1**2*(3.0d0*dx2+5.0d0*dx1)
      ci(7)=7.0d0*dx1*(dx2+3.0d0*dx1)
      ci(8)=dx2+7.0d0*dx1
      ci(9)=1.0d0
    case(2)
      ci(1)=dx1**7*dx2**2
      ci(2)=dx1**6*dx2*(7.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**5*(21.0d0*dx2**2+dx1*(14.0d0*dx2+dx1))
      ci(4)=7.0d0*dx1**4*(5.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(5)=7.0d0*dx1**3*(5.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1))
      ci(6)=7.0d0*dx1**2*(3.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))
      ci(7)=7.0d0*dx1*(dx2**2+dx1*(6.0d0*dx2+5.0d0*dx1))
      ci(8)=dx2**2+7.0d0*dx1*(2.0d0*dx2+3.0d0*dx1)
      ci(9)=2.0d0*dx2+7.0d0*dx1
      ci(10)=1.0d0
    case(3)
      ci(1)=dx1**7*dx2**3
      ci(2)=dx1**6*dx2**2*(7.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**5*dx2*(7.0d0*dx2**2+dx1*(7.0d0*dx2+dx1))
      ci(4)=dx1**4*(35.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)))
      ci(5)=7.0d0*dx1**3*(5.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(6)=21.0d0*dx1**2*(dx2**3+dx1*(5.0d0*dx2**2+dx1*(5.0d0*dx2+dx1)))
      ci(7)=7.0d0*dx1*(dx2**3+dx1*(9.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(8)=dx2**3+7.0d0*dx1*(3.0d0*dx2**2+dx1*(9.0d0*dx2+5.0d0*dx1))
      ci(9)=3.0d0*(dx2**2+7.0d0*dx1*(dx2+dx1))
      ci(10)=3.0d0*dx2+7.0d0*dx1
      ci(11)=1.0d0
    case(4)
      ci(1)=dx1**7*dx2**4
      ci(2)=dx1**6*dx2**3*(7.0d0*dx2+4.0d0*dx1)
      ci(3)=dx1**5*dx2**2*(21.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+3.0d0*dx1))
      ci(4)=dx1**4*dx2*(35.0d0*dx2**3+2.0d0*dx1*(42.0d0*dx2**2+dx1*(21.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**3*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(28.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**2*(3.0d0*dx2**4+dx1*(20.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(7)=7.0d0*dx1*(dx2**4+dx1*(12.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1))))
      ci(8)=dx2**4+7.0d0*dx1*(4.0d0*dx2**3+dx1*(18.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1)))
      ci(9)=4.0d0*dx2**3+7.0d0*dx1*(6.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(10)=6.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+3.0d0*dx1)
      ci(11)=4.0d0*dx2+7.0d0*dx1
      ci(12)=1.0d0
    case(5)
      ci(1)=dx1**7*dx2**5
      ci(2)=dx1**6*dx2**4*(7.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**5*dx2**3*(21.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))
      ci(4)=5.0d0*dx1**4*dx2**2*(7.0d0*dx2**3+dx1*(21.0d0*dx2**2+2.0d0*dx1*(7.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**3*dx2*(7.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(14.0d0*dx2+dx1))))
      ci(6)=dx1**2*(21.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(350.0d0*dx2**3+dx1*(210.0d0*dx2**2+dx1*(35.0d0*dx2+dx1)))))
      ci(7)=7.0d0*dx1*(dx2**5+dx1*(15.0d0*dx2**4+dx1*(50.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(8)=dx2**5+7.0d0*dx1*(5.0d0*dx2**4+dx1*(30.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(25.0d0*dx2+3.0d0*dx1))))
      ci(9)=5.0d0*(dx2**4+7.0d0*dx1*(2.0d0*dx2**3+dx1*(6.0d0*dx2**2+dx1*(5.0d0*dx2+dx1))))
      ci(10)=5.0d0*(2.0d0*dx2**3+7.0d0*dx1*(2.0d0*dx2**2+dx1*(3.0d0*dx2+dx1)))
      ci(11)=10.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+3.0d0*dx1)
      ci(12)=5.0d0*dx2+7.0d0*dx1
      ci(13)=1.0d0
    case(6)
      ci(1)=dx1**7*dx2**6
      ci(2)=dx1**6*dx2**5*(7.0d0*dx2+6.0d0*dx1)
      ci(3)=3.0d0*dx1**5*dx2**4*(7.0d0*dx2**2+dx1*(14.0d0*dx2+5.0d0*dx1))
      ci(4)=dx1**4*dx2**3*(35.0d0*dx2**3+dx1*(126.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+4.0d0*dx1)))
      ci(5)=5.0d0*dx1**3*dx2**2*(7.0d0*dx2**4+dx1*(42.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1))))
      ci(6)=3.0d0*dx1**2*dx2*(7.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(175.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1))) &
))
      ci(7)=dx1*(7.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(525.0d0*dx2**4+dx1*(700.0d0*dx2**3+dx1*(315.0d0*dx2**2+dx1*(42.0d0*dx2+dx1) &
)))))
      ci(8)=dx2**6+7.0d0*dx1*(6.0d0*dx2**5+dx1*(45.0d0*dx2**4+dx1*(100.0d0*dx2**3+dx1*(75.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))))
      ci(9)=3.0d0*(2.0d0*dx2**5+7.0d0*dx1*(5.0d0*dx2**4+dx1*(20.0d0*dx2**3+dx1*(25.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(10)=5.0d0*(3.0d0*dx2**4+7.0d0*dx1*(4.0d0*dx2**3+dx1*(9.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))
      ci(11)=20.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1))
      ci(12)=3.0d0*(5.0d0*dx2**2+7.0d0*dx1*(2.0d0*dx2+dx1))
      ci(13)=6.0d0*dx2+7.0d0*dx1
      ci(14)=1.0d0
    case(7)
      ci(1)=dx1**7*dx2**7
      ci(2)=7.0d0*dx1**6*dx2**6*(dx2+dx1)
      ci(3)=7.0d0*dx1**5*dx2**5*(3.0d0*dx2**2+dx1*(7.0d0*dx2+3.0d0*dx1))
      ci(4)=7.0d0*dx1**4*dx2**4*(5.0d0*dx2**3+dx1*(21.0d0*dx2**2+dx1*(21.0d0*dx2+5.0d0*dx1)))
      ci(5)=7.0d0*dx1**3*dx2**3*(5.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1*(63.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**2*dx2**2*(3.0d0*dx2**5+dx1*(35.0d0*dx2**4+dx1*(105.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1 &
)))))
      ci(7)=7.0d0*dx1*dx2*(dx2**6+dx1*(21.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(175.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(21.0d0*dx2 &
+dx1))))))
      ci(8)=dx2**7+dx1*(49.0d0*dx2**6+dx1*(441.0d0*dx2**5+dx1*(1225.0d0*dx2**4+dx1*(1225.0d0*dx2**3+dx1*(441.0d0*dx2**2+dx1* &
(49.0d0*dx2+dx1))))))
      ci(9)=7.0d0*(dx2**6+dx1*(21.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(175.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(21.0d0*dx2+dx1))))))
      ci(10)=7.0d0*(3.0d0*dx2**5+dx1*(35.0d0*dx2**4+dx1*(105.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1)))))
      ci(11)=7.0d0*(5.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1*(63.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+dx1))))
      ci(12)=7.0d0*(5.0d0*dx2**3+dx1*(21.0d0*dx2**2+dx1*(21.0d0*dx2+5.0d0*dx1)))
      ci(13)=7.0d0*(3.0d0*dx2**2+dx1*(7.0d0*dx2+3.0d0*dx1))
      ci(14)=7.0d0*(dx2+dx1)
      ci(15)=1.0d0
    case(8)
      ci(1)=dx1**7*dx2**8
      ci(2)=dx1**6*dx2**7*(7.0d0*dx2+8.0d0*dx1)
      ci(3)=7.0d0*dx1**5*dx2**6*(3.0d0*dx2**2+4.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=7.0d0*dx1**4*dx2**5*(5.0d0*dx2**3+4.0d0*dx1*(6.0d0*dx2**2+dx1*(7.0d0*dx2+2.0d0*dx1)))
      ci(5)=7.0d0*dx1**3*dx2**4*(5.0d0*dx2**4+2.0d0*dx1*(20.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(28.0d0*dx2+5.0d0*dx1))))
      ci(6)=7.0d0*dx1**2*dx2**3*(3.0d0*dx2**5+2.0d0*dx1*(20.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(35.0d0*dx2 &
+4.0d0*dx1)))))
      ci(7)=7.0d0*dx1*dx2**2*(dx2**6+2.0d0*dx1*(12.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(105.0d0*dx2**2 &
+2.0d0*dx1*(14.0d0*dx2+dx1))))))
      ci(8)=dx2*(dx2**7+2.0d0*dx1*(28.0d0*dx2**6+dx1*(294.0d0*dx2**5+dx1*(980.0d0*dx2**4+dx1*(1225.0d0*dx2**3+2.0d0*dx1* &
(294.0d0*dx2**2+dx1*(49.0d0*dx2+2.0d0*dx1)))))))
      ci(9)=8.0d0*dx2**7+dx1*(196.0d0*dx2**6+dx1*(1176.0d0*dx2**5+dx1*(2450.0d0*dx2**4+dx1*(1960.0d0*dx2**3+dx1*(588.0d0*dx2**2 &
+dx1*(56.0d0*dx2+dx1))))))
      ci(10)=7.0d0*(4.0d0*dx2**6+dx1*(56.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(24.0d0*dx2 &
+dx1))))))
      ci(11)=7.0d0*(8.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1)))))
      ci(12)=7.0d0*(10.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(84.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(13)=7.0d0*(8.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1)))
      ci(14)=7.0d0*(4.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(15)=8.0d0*dx2+7.0d0*dx1
      ci(16)=1.0d0
    case(9)
      ci(1)=dx1**7*dx2**9
      ci(2)=dx1**6*dx2**8*(7.0d0*dx2+9.0d0*dx1)
      ci(3)=3.0d0*dx1**5*dx2**7*(7.0d0*dx2**2+3.0d0*dx1*(7.0d0*dx2+4.0d0*dx1))
      ci(4)=7.0d0*dx1**4*dx2**6*(5.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+4.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(5)=7.0d0*dx1**3*dx2**5*(5.0d0*dx2**4+3.0d0*dx1*(15.0d0*dx2**3+2.0d0*dx1*(18.0d0*dx2**2+dx1*(14.0d0*dx2+3.0d0*dx1))))
      ci(6)=21.0d0*dx1**2*dx2**4*(dx2**5+3.0d0*dx1*(5.0d0*dx2**4+2.0d0*dx1*(10.0d0*dx2**3+dx1*(14.0d0*dx2**2+dx1*(7.0d0*dx2+dx1))) &
))
      ci(7)=7.0d0*dx1*dx2**3*(dx2**6+3.0d0*dx1*(9.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1* &
(21.0d0*dx2+2.0d0*dx1))))))
      ci(8)=dx2**2*(dx2**7+3.0d0*dx1*(21.0d0*dx2**6+2.0d0*dx1*(126.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(735.0d0*dx2**3+dx1* &
(441.0d0*dx2**2+2.0d0*dx1*(49.0d0*dx2+3.0d0*dx1)))))))
      ci(9)=9.0d0*dx2*(dx2**7+dx1*(28.0d0*dx2**6+dx1*(196.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(490.0d0*dx2**3+dx1*(196.0d0*dx2**2 &
+dx1*(28.0d0*dx2+dx1)))))))
      ci(10)=36.0d0*dx2**7+dx1*(588.0d0*dx2**6+dx1*(2646.0d0*dx2**5+dx1*(4410.0d0*dx2**4+dx1*(2940.0d0*dx2**3+dx1*(756.0d0*dx2**2 &
+dx1*(63.0d0*dx2+dx1))))))
      ci(11)=7.0d0*(12.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(378.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1*(27.0d0*dx2 &
+dx1))))))
      ci(12)=21.0d0*(6.0d0*dx2**5+dx1*(42.0d0*dx2**4+dx1*(84.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(13)=7.0d0*(18.0d0*dx2**4+dx1*(84.0d0*dx2**3+dx1*(108.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+dx1))))
      ci(14)=7.0d0*(12.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(27.0d0*dx2+5.0d0*dx1)))
      ci(15)=3.0d0*(12.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+dx1))
      ci(16)=9.0d0*dx2+7.0d0*dx1
      ci(17)=1.0d0
    case(10)
      ci(1)=dx1**7*dx2**10
      ci(2)=dx1**6*dx2**9*(7.0d0*dx2+10.0d0*dx1)
      ci(3)=dx1**5*dx2**8*(21.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+9.0d0*dx1))
      ci(4)=5.0d0*dx1**4*dx2**7*(7.0d0*dx2**3+3.0d0*dx1*(14.0d0*dx2**2+dx1*(21.0d0*dx2+8.0d0*dx1)))
      ci(5)=35.0d0*dx1**3*dx2**6*(dx2**4+dx1*(10.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**2*dx2**5*(3.0d0*dx2**5+dx1*(50.0d0*dx2**4+3.0d0*dx1*(75.0d0*dx2**3+2.0d0*dx1*(60.0d0*dx2**2+dx1*(35.0d0*dx2 &
+6.0d0*dx1)))))
      ci(7)=7.0d0*dx1*dx2**4*(dx2**6+3.0d0*dx1*(10.0d0*dx2**5+dx1*(75.0d0*dx2**4+2.0d0*dx1*(100.0d0*dx2**3+dx1*(105.0d0*dx2**2 &
+dx1*(42.0d0*dx2+5.0d0*dx1))))))
      ci(8)=dx2**3*(dx2**7+dx1*(70.0d0*dx2**6+3.0d0*dx1*(315.0d0*dx2**5+2.0d0*dx1*(700.0d0*dx2**4+dx1*(1225.0d0*dx2**3+dx1* &
(882.0d0*dx2**2+5.0d0*dx1*(49.0d0*dx2+4.0d0*dx1)))))))
      ci(9)=5.0d0*dx2**2*(2.0d0*dx2**7+3.0d0*dx1*(21.0d0*dx2**6+dx1*(168.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(588.0d0*dx2**3+dx1* &
(294.0d0*dx2**2+dx1*(56.0d0*dx2+3.0d0*dx1)))))))
      ci(10)=5.0d0*dx2*(9.0d0*dx2**7+dx1*(168.0d0*dx2**6+dx1*(882.0d0*dx2**5+dx1*(1764.0d0*dx2**4+dx1*(1470.0d0*dx2**3+dx1* &
(504.0d0*dx2**2+dx1*(63.0d0*dx2+2.0d0*dx1)))))))
      ci(11)=120.0d0*dx2**7+dx1*(1470.0d0*dx2**6+dx1*(5292.0d0*dx2**5+dx1*(7350.0d0*dx2**4+dx1*(4200.0d0*dx2**3+dx1* &
(945.0d0*dx2**2+dx1*(70.0d0*dx2+dx1))))))
      ci(12)=7.0d0*(30.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1*(600.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1*(30.0d0*dx2 &
+dx1))))))
      ci(13)=7.0d0*(36.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(360.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1*(50.0d0*dx2+3.0d0*dx1)))))
      ci(14)=35.0d0*(6.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(15)=5.0d0*(24.0d0*dx2**3+7.0d0*dx1*(9.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))
      ci(16)=45.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+3.0d0*dx1)
      ci(17)=10.0d0*dx2+7.0d0*dx1
      ci(18)=1.0d0
    case(11)
      ci(1)=dx1**7*dx2**11
      ci(2)=dx1**6*dx2**10*(7.0d0*dx2+11.0d0*dx1)
      ci(3)=dx1**5*dx2**9*(21.0d0*dx2**2+11.0d0*dx1*(7.0d0*dx2+5.0d0*dx1))
      ci(4)=dx1**4*dx2**8*(35.0d0*dx2**3+11.0d0*dx1*(21.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+3.0d0*dx1)))
      ci(5)=5.0d0*dx1**3*dx2**7*(7.0d0*dx2**4+11.0d0*dx1*(7.0d0*dx2**3+3.0d0*dx1*(7.0d0*dx2**2+dx1*(7.0d0*dx2+2.0d0*dx1))))
      ci(6)=7.0d0*dx1**2*dx2**6*(3.0d0*dx2**5+11.0d0*dx1*(5.0d0*dx2**4+dx1*(25.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1* &
(5.0d0*dx2+dx1)))))
      ci(7)=7.0d0*dx1*dx2**5*(dx2**6+11.0d0*dx1*(3.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(25.0d0*dx2**3+2.0d0*dx1*(15.0d0*dx2**2 &
+dx1*(7.0d0*dx2+dx1))))))
      ci(8)=dx2**4*(dx2**7+11.0d0*dx1*(7.0d0*dx2**6+3.0d0*dx1*(35.0d0*dx2**5+dx1*(175.0d0*dx2**4+2.0d0*dx1*(175.0d0*dx2**3+dx1* &
(147.0d0*dx2**2+dx1*(49.0d0*dx2+5.0d0*dx1)))))))
      ci(9)=11.0d0*dx2**3*(dx2**7+dx1*(35.0d0*dx2**6+3.0d0*dx1*(105.0d0*dx2**5+dx1*(350.0d0*dx2**4+dx1*(490.0d0*dx2**3+dx1* &
(294.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+dx1)))))))
      ci(10)=55.0d0*dx2**2*(dx2**7+dx1*(21.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(294.0d0*dx2**4+dx1*(294.0d0*dx2**3+dx1* &
(126.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)))))))
      ci(11)=11.0d0*dx2*(15.0d0*dx2**7+dx1*(210.0d0*dx2**6+dx1*(882.0d0*dx2**5+dx1*(1470.0d0*dx2**4+dx1*(1050.0d0*dx2**3+dx1* &
(315.0d0*dx2**2+dx1*(35.0d0*dx2+dx1)))))))
      ci(12)=330.0d0*dx2**7+dx1*(3234.0d0*dx2**6+dx1*(9702.0d0*dx2**5+dx1*(11550.0d0*dx2**4+dx1*(5775.0d0*dx2**3+dx1* &
(1155.0d0*dx2**2+dx1*(77.0d0*dx2+dx1))))))
      ci(13)=7.0d0*(66.0d0*dx2**6+dx1*(462.0d0*dx2**5+dx1*(990.0d0*dx2**4+dx1*(825.0d0*dx2**3+dx1*(275.0d0*dx2**2+dx1*(33.0d0*dx2 &
+dx1))))))
      ci(14)=7.0d0*(66.0d0*dx2**5+dx1*(330.0d0*dx2**4+dx1*(495.0d0*dx2**3+dx1*(275.0d0*dx2**2+dx1*(55.0d0*dx2+3.0d0*dx1)))))
      ci(15)=5.0d0*(66.0d0*dx2**4+7.0d0*dx1*(33.0d0*dx2**3+dx1*(33.0d0*dx2**2+dx1*(11.0d0*dx2+dx1))))
      ci(16)=165.0d0*dx2**3+7.0d0*dx1*(55.0d0*dx2**2+dx1*(33.0d0*dx2+5.0d0*dx1))
      ci(17)=55.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+3.0d0*dx1)
      ci(18)=11.0d0*dx2+7.0d0*dx1
      ci(19)=1.0d0
    case(12)
      ci(1)=dx1**7*dx2**12
      ci(2)=dx1**6*dx2**11*(7.0d0*dx2+12.0d0*dx1)
      ci(3)=3.0d0*dx1**5*dx2**10*(7.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+11.0d0*dx1))
      ci(4)=dx1**4*dx2**9*(35.0d0*dx2**3+2.0d0*dx1*(126.0d0*dx2**2+11.0d0*dx1*(21.0d0*dx2+10.0d0*dx1)))
      ci(5)=dx1**3*dx2**8*(35.0d0*dx2**4+dx1*(420.0d0*dx2**3+11.0d0*dx1*(126.0d0*dx2**2+5.0d0*dx1*(28.0d0*dx2+9.0d0*dx1))))
      ci(6)=3.0d0*dx1**2*dx2**7*(7.0d0*dx2**5+dx1*(140.0d0*dx2**4+11.0d0*dx1*(70.0d0*dx2**3+dx1*(140.0d0*dx2**2+3.0d0*dx1* &
(35.0d0*dx2+8.0d0*dx1)))))
      ci(7)=7.0d0*dx1*dx2**6*(dx2**6+dx1*(36.0d0*dx2**5+11.0d0*dx1*(30.0d0*dx2**4+dx1*(100.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2 &
+4.0d0*dx1*(6.0d0*dx2+dx1))))))
      ci(8)=dx2**5*(dx2**7+dx1*(84.0d0*dx2**6+11.0d0*dx1*(126.0d0*dx2**5+dx1*(700.0d0*dx2**4+3.0d0*dx1*(525.0d0*dx2**3+4.0d0*dx1* &
(126.0d0*dx2**2+dx1*(49.0d0*dx2+6.0d0*dx1)))))))
      ci(9)=3.0d0*dx2**4*(4.0d0*dx2**7+11.0d0*dx1*(14.0d0*dx2**6+dx1*(140.0d0*dx2**5+3.0d0*dx1*(175.0d0*dx2**4+dx1*(280.0d0*dx2**3 &
+dx1*(196.0d0*dx2**2+dx1*(56.0d0*dx2+5.0d0*dx1)))))))
      ci(10)=11.0d0*dx2**3*(6.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(945.0d0*dx2**5+dx1*(2520.0d0*dx2**4+dx1*(2940.0d0*dx2**3+dx1* &
(1512.0d0*dx2**2+5.0d0*dx1*(63.0d0*dx2+4.0d0*dx1)))))))
      ci(11)=11.0d0*dx2**2*(20.0d0*dx2**7+dx1*(315.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1*(2520.0d0*dx2**3+dx1* &
(945.0d0*dx2**2+2.0d0*dx1*(70.0d0*dx2+3.0d0*dx1)))))))
      ci(12)=3.0d0*dx2*(165.0d0*dx2**7+dx1*(1848.0d0*dx2**6+dx1*(6468.0d0*dx2**5+dx1*(9240.0d0*dx2**4+dx1*(5775.0d0*dx2**3 &
+2.0d0*dx1*(770.0d0*dx2**2+dx1*(77.0d0*dx2+2.0d0*dx1)))))))
      ci(13)=792.0d0*dx2**7+dx1*(6468.0d0*dx2**6+dx1*(16632.0d0*dx2**5+dx1*(17325.0d0*dx2**4+dx1*(7700.0d0*dx2**3+dx1* &
(1386.0d0*dx2**2+dx1*(84.0d0*dx2+dx1))))))
      ci(14)=7.0d0*(132.0d0*dx2**6+dx1*(792.0d0*dx2**5+dx1*(1485.0d0*dx2**4+dx1*(1100.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1* &
(36.0d0*dx2+dx1))))))
      ci(15)=3.0d0*(264.0d0*dx2**5+7.0d0*dx1*(165.0d0*dx2**4+dx1*(220.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)))))
      ci(16)=495.0d0*dx2**4+7.0d0*dx1*(220.0d0*dx2**3+dx1*(198.0d0*dx2**2+5.0d0*dx1*(12.0d0*dx2+dx1)))
      ci(17)=220.0d0*dx2**3+7.0d0*dx1*(66.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))
      ci(18)=3.0d0*(22.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+dx1))
      ci(19)=12.0d0*dx2+7.0d0*dx1
      ci(20)=1.0d0
    case(13)
      ci(1)=dx1**7*dx2**13
      ci(2)=dx1**6*dx2**12*(7.0d0*dx2+13.0d0*dx1)
      ci(3)=dx1**5*dx2**11*(21.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+6.0d0*dx1))
      ci(4)=dx1**4*dx2**10*(35.0d0*dx2**3+13.0d0*dx1*(21.0d0*dx2**2+2.0d0*dx1*(21.0d0*dx2+11.0d0*dx1)))
      ci(5)=dx1**3*dx2**9*(35.0d0*dx2**4+13.0d0*dx1*(35.0d0*dx2**3+dx1*(126.0d0*dx2**2+11.0d0*dx1*(14.0d0*dx2+5.0d0*dx1))))
      ci(6)=dx1**2*dx2**8*(21.0d0*dx2**5+13.0d0*dx1*(35.0d0*dx2**4+dx1*(210.0d0*dx2**3+11.0d0*dx1*(42.0d0*dx2**2+dx1*(35.0d0*dx2 &
+9.0d0*dx1)))))
      ci(7)=dx1*dx2**7*(7.0d0*dx2**6+13.0d0*dx1*(21.0d0*dx2**5+dx1*(210.0d0*dx2**4+11.0d0*dx1*(70.0d0*dx2**3+3.0d0*dx1* &
(35.0d0*dx2**2+dx1*(21.0d0*dx2+4.0d0*dx1))))))
      ci(8)=dx2**6*(dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(126.0d0*dx2**5+11.0d0*dx1*(70.0d0*dx2**4+dx1*(175.0d0*dx2**3+3.0d0*dx1* &
(63.0d0*dx2**2+4.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(9)=13.0d0*dx2**5*(dx2**7+dx1*(42.0d0*dx2**6+11.0d0*dx1*(42.0d0*dx2**5+dx1*(175.0d0*dx2**4+3.0d0*dx1*(105.0d0*dx2**3+dx1* &
(84.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1)))))))
      ci(10)=13.0d0*dx2**4*(6.0d0*dx2**7+11.0d0*dx1*(14.0d0*dx2**6+dx1*(105.0d0*dx2**5+dx1*(315.0d0*dx2**4+dx1*(420.0d0*dx2**3 &
+dx1*(252.0d0*dx2**2+dx1*(63.0d0*dx2+5.0d0*dx1)))))))
      ci(11)=143.0d0*dx2**3*(2.0d0*dx2**7+dx1*(35.0d0*dx2**6+dx1*(189.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1* &
(189.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1)))))))
      ci(12)=13.0d0*dx2**2*(55.0d0*dx2**7+dx1*(693.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(4620.0d0*dx2**4+dx1*(3465.0d0*dx2**3+dx1* &
(1155.0d0*dx2**2+2.0d0*dx1*(77.0d0*dx2+3.0d0*dx1)))))))
      ci(13)=13.0d0*dx2*(99.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(3465.0d0*dx2**4+dx1*(1925.0d0*dx2**3+dx1* &
(462.0d0*dx2**2+dx1*(42.0d0*dx2+dx1)))))))
      ci(14)=1716.0d0*dx2**7+dx1*(12012.0d0*dx2**6+dx1*(27027.0d0*dx2**5+dx1*(25025.0d0*dx2**4+dx1*(10010.0d0*dx2**3+dx1* &
(1638.0d0*dx2**2+dx1*(91.0d0*dx2+dx1))))))
      ci(15)=1716.0d0*dx2**6+7.0d0*dx1*(1287.0d0*dx2**5+dx1*(2145.0d0*dx2**4+dx1*(1430.0d0*dx2**3+dx1*(390.0d0*dx2**2+dx1* &
(39.0d0*dx2+dx1)))))
      ci(16)=1287.0d0*dx2**5+7.0d0*dx1*(715.0d0*dx2**4+dx1*(858.0d0*dx2**3+dx1*(390.0d0*dx2**2+dx1*(65.0d0*dx2+3.0d0*dx1))))
      ci(17)=715.0d0*dx2**4+7.0d0*dx1*(286.0d0*dx2**3+dx1*(234.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+dx1)))
      ci(18)=286.0d0*dx2**3+7.0d0*dx1*(78.0d0*dx2**2+dx1*(39.0d0*dx2+5.0d0*dx1))
      ci(19)=78.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+3.0d0*dx1)
      ci(20)=13.0d0*dx2+7.0d0*dx1
      ci(21)=1.0d0
    case(14)
      ci(1)=dx1**7*dx2**14
      ci(2)=7.0d0*dx1**6*dx2**13*(dx2+2.0d0*dx1)
      ci(3)=7.0d0*dx1**5*dx2**12*(3.0d0*dx2**2+dx1*(14.0d0*dx2+13.0d0*dx1))
      ci(4)=7.0d0*dx1**4*dx2**11*(5.0d0*dx2**3+dx1*(42.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+4.0d0*dx1)))
      ci(5)=7.0d0*dx1**3*dx2**10*(5.0d0*dx2**4+dx1*(70.0d0*dx2**3+13.0d0*dx1*(21.0d0*dx2**2+dx1*(28.0d0*dx2+11.0d0*dx1))))
      ci(6)=7.0d0*dx1**2*dx2**9*(3.0d0*dx2**5+dx1*(70.0d0*dx2**4+13.0d0*dx1*(35.0d0*dx2**3+dx1*(84.0d0*dx2**2+11.0d0*dx1* &
(7.0d0*dx2+2.0d0*dx1)))))
      ci(7)=7.0d0*dx1*dx2**8*(dx2**6+dx1*(42.0d0*dx2**5+13.0d0*dx1*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+11.0d0*dx1*(21.0d0*dx2**2 &
+dx1*(14.0d0*dx2+3.0d0*dx1))))))
      ci(8)=dx2**7*(dx2**7+dx1*(98.0d0*dx2**6+13.0d0*dx1*(147.0d0*dx2**5+dx1*(980.0d0*dx2**4+11.0d0*dx1*(245.0d0*dx2**3+3.0d0*dx1* &
(98.0d0*dx2**2+dx1*(49.0d0*dx2+8.0d0*dx1)))))))
      ci(9)=7.0d0*dx2**6*(2.0d0*dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(84.0d0*dx2**5+11.0d0*dx1*(35.0d0*dx2**4+dx1*(70.0d0*dx2**3 &
+3.0d0*dx1*(21.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))))))
      ci(10)=91.0d0*dx2**5*(dx2**7+dx1*(28.0d0*dx2**6+11.0d0*dx1*(21.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(105.0d0*dx2**3+dx1* &
(72.0d0*dx2**2+dx1*(21.0d0*dx2+2.0d0*dx1)))))))
      ci(11)=91.0d0*dx2**4*(4.0d0*dx2**7+11.0d0*dx1*(7.0d0*dx2**6+dx1*(42.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1* &
(63.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)))))))
      ci(12)=91.0d0*dx2**3*(11.0d0*dx2**7+dx1*(154.0d0*dx2**6+dx1*(693.0d0*dx2**5+dx1*(1320.0d0*dx2**4+dx1*(1155.0d0*dx2**3+dx1* &
(462.0d0*dx2**2+dx1*(77.0d0*dx2+4.0d0*dx1)))))))
      ci(13)=91.0d0*dx2**2*(22.0d0*dx2**7+dx1*(231.0d0*dx2**6+dx1*(792.0d0*dx2**5+dx1*(1155.0d0*dx2**4+dx1*(770.0d0*dx2**3+dx1* &
(231.0d0*dx2**2+dx1*(28.0d0*dx2+dx1)))))))
      ci(14)=7.0d0*dx2*(429.0d0*dx2**7+dx1*(3432.0d0*dx2**6+dx1*(9009.0d0*dx2**5+dx1*(10010.0d0*dx2**4+dx1*(5005.0d0*dx2**3+dx1* &
(1092.0d0*dx2**2+dx1*(91.0d0*dx2+2.0d0*dx1)))))))
      ci(15)=3432.0d0*dx2**7+dx1*(21021.0d0*dx2**6+dx1*(42042.0d0*dx2**5+dx1*(35035.0d0*dx2**4+dx1*(12740.0d0*dx2**3+dx1* &
(1911.0d0*dx2**2+dx1*(98.0d0*dx2+dx1))))))
      ci(16)=7.0d0*(429.0d0*dx2**6+dx1*(2002.0d0*dx2**5+dx1*(3003.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1*(455.0d0*dx2**2+dx1* &
(42.0d0*dx2+dx1))))))
      ci(17)=7.0d0*(286.0d0*dx2**5+dx1*(1001.0d0*dx2**4+dx1*(1092.0d0*dx2**3+dx1*(455.0d0*dx2**2+dx1*(70.0d0*dx2+3.0d0*dx1)))))
      ci(18)=7.0d0*(143.0d0*dx2**4+dx1*(364.0d0*dx2**3+dx1*(273.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+dx1))))
      ci(19)=7.0d0*(52.0d0*dx2**3+dx1*(91.0d0*dx2**2+dx1*(42.0d0*dx2+5.0d0*dx1)))
      ci(20)=7.0d0*(13.0d0*dx2**2+dx1*(14.0d0*dx2+3.0d0*dx1))
      ci(21)=7.0d0*(2.0d0*dx2+dx1)
      ci(22)=1.0d0
    case(15)
      ci(1)=dx1**7*dx2**15
      ci(2)=dx1**6*dx2**14*(7.0d0*dx2+15.0d0*dx1)
      ci(3)=21.0d0*dx1**5*dx2**13*(dx2**2+5.0d0*dx1*(dx2+dx1))
      ci(4)=35.0d0*dx1**4*dx2**12*(dx2**3+dx1*(9.0d0*dx2**2+dx1*(21.0d0*dx2+13.0d0*dx1)))
      ci(5)=35.0d0*dx1**3*dx2**11*(dx2**4+dx1*(15.0d0*dx2**3+dx1*(63.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+3.0d0*dx1))))
      ci(6)=21.0d0*dx1**2*dx2**10*(dx2**5+dx1*(25.0d0*dx2**4+dx1*(175.0d0*dx2**3+13.0d0*dx1*(35.0d0*dx2**2+dx1*(35.0d0*dx2 &
+11.0d0*dx1)))))
      ci(7)=7.0d0*dx1*dx2**9*(dx2**6+dx1*(45.0d0*dx2**5+dx1*(525.0d0*dx2**4+13.0d0*dx1*(175.0d0*dx2**3+dx1*(315.0d0*dx2**2 &
+11.0d0*dx1*(21.0d0*dx2+5.0d0*dx1))))))
      ci(8)=dx2**8*(dx2**7+dx1*(105.0d0*dx2**6+dx1*(2205.0d0*dx2**5+13.0d0*dx1*(1225.0d0*dx2**4+dx1*(3675.0d0*dx2**3+11.0d0*dx1* &
(441.0d0*dx2**2+5.0d0*dx1*(49.0d0*dx2+9.0d0*dx1)))))))
      ci(9)=15.0d0*dx2**7*(dx2**7+dx1*(49.0d0*dx2**6+13.0d0*dx1*(49.0d0*dx2**5+dx1*(245.0d0*dx2**4+11.0d0*dx1*(49.0d0*dx2**3+dx1* &
(49.0d0*dx2**2+3.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(10)=35.0d0*dx2**6*(3.0d0*dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(63.0d0*dx2**5+11.0d0*dx1*(21.0d0*dx2**4+dx1*(35.0d0*dx2**3 &
+dx1*(27.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))))))
      ci(11)=91.0d0*dx2**5*(5.0d0*dx2**7+dx1*(105.0d0*dx2**6+11.0d0*dx1*(63.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(225.0d0*dx2**3 &
+dx1*(135.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1)))))))
      ci(12)=273.0d0*dx2**4*(5.0d0*dx2**7+dx1*(77.0d0*dx2**6+dx1*(385.0d0*dx2**5+dx1*(825.0d0*dx2**4+dx1*(825.0d0*dx2**3+dx1* &
(385.0d0*dx2**2+dx1*(77.0d0*dx2+5.0d0*dx1)))))))
      ci(13)=91.0d0*dx2**3*(33.0d0*dx2**7+dx1*(385.0d0*dx2**6+dx1*(1485.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(1925.0d0*dx2**3+dx1* &
(693.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+dx1)))))))
      ci(14)=35.0d0*dx2**2*(143.0d0*dx2**7+dx1*(1287.0d0*dx2**6+dx1*(3861.0d0*dx2**5+dx1*(5005.0d0*dx2**4+dx1*(3003.0d0*dx2**3 &
+dx1*(819.0d0*dx2**2+dx1*(91.0d0*dx2+3.0d0*dx1)))))))
      ci(15)=15.0d0*dx2*(429.0d0*dx2**7+dx1*(3003.0d0*dx2**6+dx1*(7007.0d0*dx2**5+dx1*(7007.0d0*dx2**4+dx1*(3185.0d0*dx2**3+dx1* &
(637.0d0*dx2**2+dx1*(49.0d0*dx2+dx1)))))))
      ci(16)=6435.0d0*dx2**7+dx1*(35035.0d0*dx2**6+dx1*(63063.0d0*dx2**5+dx1*(47775.0d0*dx2**4+dx1*(15925.0d0*dx2**3+dx1* &
(2205.0d0*dx2**2+dx1*(105.0d0*dx2+dx1))))))
      ci(17)=7.0d0*(715.0d0*dx2**6+dx1*(3003.0d0*dx2**5+dx1*(4095.0d0*dx2**4+dx1*(2275.0d0*dx2**3+dx1*(525.0d0*dx2**2+dx1* &
(45.0d0*dx2+dx1))))))
      ci(18)=21.0d0*(143.0d0*dx2**5+dx1*(455.0d0*dx2**4+dx1*(455.0d0*dx2**3+dx1*(175.0d0*dx2**2+dx1*(25.0d0*dx2+dx1)))))
      ci(19)=35.0d0*(39.0d0*dx2**4+dx1*(91.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))))
      ci(20)=35.0d0*(13.0d0*dx2**3+dx1*(21.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(21)=21.0d0*(5.0d0*dx2**2+dx1*(5.0d0*dx2+dx1))
      ci(22)=15.0d0*dx2+7.0d0*dx1
      ci(23)=1.0d0
    case(16)
      ci(1)=dx1**7*dx2**16
      ci(2)=dx1**6*dx2**15*(7.0d0*dx2+16.0d0*dx1)
      ci(3)=dx1**5*dx2**14*(21.0d0*dx2**2+8.0d0*dx1*(14.0d0*dx2+15.0d0*dx1))
      ci(4)=7.0d0*dx1**4*dx2**13*(5.0d0*dx2**3+8.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(5)=35.0d0*dx1**3*dx2**12*(dx2**4+4.0d0*dx1*(4.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(28.0d0*dx2+13.0d0*dx1))))
      ci(6)=7.0d0*dx1**2*dx2**11*(3.0d0*dx2**5+4.0d0*dx1*(20.0d0*dx2**4+dx1*(150.0d0*dx2**3+dx1*(420.0d0*dx2**2+13.0d0*dx1* &
(35.0d0*dx2+12.0d0*dx1)))))
      ci(7)=7.0d0*dx1*dx2**10*(dx2**6+4.0d0*dx1*(12.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(700.0d0*dx2**3+13.0d0*dx1*(105.0d0*dx2**2 &
+2.0d0*dx1*(42.0d0*dx2+11.0d0*dx1))))))
      ci(8)=dx2**9*(dx2**7+4.0d0*dx1*(28.0d0*dx2**6+dx1*(630.0d0*dx2**5+dx1*(4900.0d0*dx2**4+13.0d0*dx1*(1225.0d0*dx2**3 &
+2.0d0*dx1*(882.0d0*dx2**2+11.0d0*dx1*(49.0d0*dx2+10.0d0*dx1)))))))
      ci(9)=2.0d0*dx2**8*(8.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(5880.0d0*dx2**5+13.0d0*dx1*(2450.0d0*dx2**4+dx1*(5880.0d0*dx2**3 &
+11.0d0*dx1*(588.0d0*dx2**2+5.0d0*dx1*(56.0d0*dx2+9.0d0*dx1)))))))
      ci(10)=10.0d0*dx2**7*(12.0d0*dx2**7+dx1*(392.0d0*dx2**6+13.0d0*dx1*(294.0d0*dx2**5+dx1*(1176.0d0*dx2**4+11.0d0*dx1* &
(196.0d0*dx2**3+dx1*(168.0d0*dx2**2+dx1*(63.0d0*dx2+8.0d0*dx1)))))))
      ci(11)=14.0d0*dx2**6*(40.0d0*dx2**7+13.0d0*dx1*(70.0d0*dx2**6+dx1*(504.0d0*dx2**5+11.0d0*dx1*(140.0d0*dx2**4+dx1* &
(200.0d0*dx2**3+dx1*(135.0d0*dx2**2+4.0d0*dx1*(10.0d0*dx2+dx1)))))))
      ci(12)=182.0d0*dx2**5*(10.0d0*dx2**7+dx1*(168.0d0*dx2**6+dx1*(924.0d0*dx2**5+dx1*(2200.0d0*dx2**4+dx1*(2475.0d0*dx2**3 &
+4.0d0*dx1*(330.0d0*dx2**2+dx1*(77.0d0*dx2+6.0d0*dx1)))))))
      ci(13)=182.0d0*dx2**4*(24.0d0*dx2**7+dx1*(308.0d0*dx2**6+dx1*(1320.0d0*dx2**5+dx1*(2475.0d0*dx2**4+2.0d0*dx1* &
(1100.0d0*dx2**3+dx1*(462.0d0*dx2**2+dx1*(84.0d0*dx2+5.0d0*dx1)))))))
      ci(14)=14.0d0*dx2**3*(572.0d0*dx2**7+dx1*(5720.0d0*dx2**6+dx1*(19305.0d0*dx2**5+2.0d0*dx1*(14300.0d0*dx2**4+dx1* &
(10010.0d0*dx2**3+dx1*(3276.0d0*dx2**2+5.0d0*dx1*(91.0d0*dx2+4.0d0*dx1)))))))
      ci(15)=10.0d0*dx2**2*(1144.0d0*dx2**7+dx1*(9009.0d0*dx2**6+2.0d0*dx1*(12012.0d0*dx2**5+dx1*(14014.0d0*dx2**4+dx1* &
(7644.0d0*dx2**3+dx1*(1911.0d0*dx2**2+2.0d0*dx1*(98.0d0*dx2+3.0d0*dx1)))))))
      ci(16)=2.0d0*dx2*(6435.0d0*dx2**7+2.0d0*dx1*(20020.0d0*dx2**6+dx1*(42042.0d0*dx2**5+dx1*(38220.0d0*dx2**4+dx1* &
(15925.0d0*dx2**3+2.0d0*dx1*(1470.0d0*dx2**2+dx1*(105.0d0*dx2+2.0d0*dx1)))))))
      ci(17)=11440.0d0*dx2**7+dx1*(56056.0d0*dx2**6+dx1*(91728.0d0*dx2**5+dx1*(63700.0d0*dx2**4+dx1*(19600.0d0*dx2**3+dx1* &
(2520.0d0*dx2**2+dx1*(112.0d0*dx2+dx1))))))
      ci(18)=7.0d0*(1144.0d0*dx2**6+dx1*(4368.0d0*dx2**5+dx1*(5460.0d0*dx2**4+dx1*(2800.0d0*dx2**3+dx1*(600.0d0*dx2**2+dx1* &
(48.0d0*dx2+dx1))))))
      ci(19)=7.0d0*(624.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1680.0d0*dx2**3+dx1*(600.0d0*dx2**2+dx1*(80.0d0*dx2+3.0d0*dx1)))))
      ci(20)=35.0d0*(52.0d0*dx2**4+dx1*(112.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(21)=7.0d0*(80.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(48.0d0*dx2+5.0d0*dx1)))
      ci(22)=120.0d0*dx2**2+7.0d0*dx1*(16.0d0*dx2+3.0d0*dx1)
      ci(23)=16.0d0*dx2+7.0d0*dx1
      ci(24)=1.0d0
    case(17)
      ci(1)=dx1**7*dx2**17
      ci(2)=dx1**6*dx2**16*(7.0d0*dx2+17.0d0*dx1)
      ci(3)=dx1**5*dx2**15*(21.0d0*dx2**2+17.0d0*dx1*(7.0d0*dx2+8.0d0*dx1))
      ci(4)=dx1**4*dx2**14*(35.0d0*dx2**3+17.0d0*dx1*(21.0d0*dx2**2+8.0d0*dx1*(7.0d0*dx2+5.0d0*dx1)))
      ci(5)=7.0d0*dx1**3*dx2**13*(5.0d0*dx2**4+17.0d0*dx1*(5.0d0*dx2**3+4.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**2*dx2**12*(3.0d0*dx2**5+17.0d0*dx1*(5.0d0*dx2**4+4.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1* &
(35.0d0*dx2+13.0d0*dx1)))))
      ci(7)=7.0d0*dx1*dx2**11*(dx2**6+17.0d0*dx1*(3.0d0*dx2**5+4.0d0*dx1*(10.0d0*dx2**4+dx1*(50.0d0*dx2**3+dx1*(105.0d0*dx2**2 &
+13.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))))))
      ci(8)=dx2**10*(dx2**7+17.0d0*dx1*(7.0d0*dx2**6+4.0d0*dx1*(42.0d0*dx2**5+dx1*(350.0d0*dx2**4+dx1*(1225.0d0*dx2**3+13.0d0*dx1* &
(147.0d0*dx2**2+2.0d0*dx1*(49.0d0*dx2+11.0d0*dx1)))))))
      ci(9)=17.0d0*dx2**9*(dx2**7+2.0d0*dx1*(28.0d0*dx2**6+dx1*(420.0d0*dx2**5+dx1*(2450.0d0*dx2**4+13.0d0*dx1*(490.0d0*dx2**3 &
+dx1*(588.0d0*dx2**2+11.0d0*dx1*(28.0d0*dx2+5.0d0*dx1)))))))
      ci(10)=34.0d0*dx2**8*(4.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(1470.0d0*dx2**5+13.0d0*dx1*(490.0d0*dx2**4+dx1*(980.0d0*dx2**3 &
+11.0d0*dx1*(84.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(11)=34.0d0*dx2**7*(20.0d0*dx2**7+dx1*(490.0d0*dx2**6+13.0d0*dx1*(294.0d0*dx2**5+dx1*(980.0d0*dx2**4+11.0d0*dx1* &
(140.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(35.0d0*dx2+4.0d0*dx1)))))))
      ci(12)=238.0d0*dx2**6*(10.0d0*dx2**7+13.0d0*dx1*(14.0d0*dx2**6+dx1*(84.0d0*dx2**5+dx1*(220.0d0*dx2**4+dx1*(275.0d0*dx2**3 &
+dx1*(165.0d0*dx2**2+4.0d0*dx1*(11.0d0*dx2+dx1)))))))
      ci(13)=3094.0d0*dx2**5*(2.0d0*dx2**7+dx1*(28.0d0*dx2**6+dx1*(132.0d0*dx2**5+dx1*(275.0d0*dx2**4+dx1*(275.0d0*dx2**3 &
+2.0d0*dx1*(66.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)))))))
      ci(14)=238.0d0*dx2**4*(52.0d0*dx2**7+dx1*(572.0d0*dx2**6+dx1*(2145.0d0*dx2**5+dx1*(3575.0d0*dx2**4+2.0d0*dx1* &
(1430.0d0*dx2**3+dx1*(546.0d0*dx2**2+dx1*(91.0d0*dx2+5.0d0*dx1)))))))
      ci(15)=34.0d0*dx2**3*(572.0d0*dx2**7+dx1*(5005.0d0*dx2**6+dx1*(15015.0d0*dx2**5+2.0d0*dx1*(10010.0d0*dx2**4+dx1* &
(6370.0d0*dx2**3+dx1*(1911.0d0*dx2**2+5.0d0*dx1*(49.0d0*dx2+2.0d0*dx1)))))))
      ci(16)=34.0d0*dx2**2*(715.0d0*dx2**7+dx1*(5005.0d0*dx2**6+2.0d0*dx1*(6006.0d0*dx2**5+dx1*(6370.0d0*dx2**4+dx1* &
(3185.0d0*dx2**3+dx1*(735.0d0*dx2**2+2.0d0*dx1*(35.0d0*dx2+dx1)))))))
      ci(17)=17.0d0*dx2*(1430.0d0*dx2**7+dx1*(8008.0d0*dx2**6+dx1*(15288.0d0*dx2**5+dx1*(12740.0d0*dx2**4+dx1*(4900.0d0*dx2**3 &
+dx1*(840.0d0*dx2**2+dx1*(56.0d0*dx2+dx1)))))))
      ci(18)=19448.0d0*dx2**7+dx1*(86632.0d0*dx2**6+dx1*(129948.0d0*dx2**5+dx1*(83300.0d0*dx2**4+dx1*(23800.0d0*dx2**3+dx1* &
(2856.0d0*dx2**2+dx1*(119.0d0*dx2+dx1))))))
      ci(19)=7.0d0*(1768.0d0*dx2**6+dx1*(6188.0d0*dx2**5+dx1*(7140.0d0*dx2**4+dx1*(3400.0d0*dx2**3+dx1*(680.0d0*dx2**2+dx1* &
(51.0d0*dx2+dx1))))))
      ci(20)=7.0d0*(884.0d0*dx2**5+dx1*(2380.0d0*dx2**4+dx1*(2040.0d0*dx2**3+dx1*(680.0d0*dx2**2+dx1*(85.0d0*dx2+3.0d0*dx1)))))
      ci(21)=7.0d0*(340.0d0*dx2**4+dx1*(680.0d0*dx2**3+dx1*(408.0d0*dx2**2+5.0d0*dx1*(17.0d0*dx2+dx1))))
      ci(22)=680.0d0*dx2**3+7.0d0*dx1*(136.0d0*dx2**2+dx1*(51.0d0*dx2+5.0d0*dx1))
      ci(23)=136.0d0*dx2**2+7.0d0*dx1*(17.0d0*dx2+3.0d0*dx1)
      ci(24)=17.0d0*dx2+7.0d0*dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>17, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^8 * (x-x2)^n2 as sum_k=0^(8+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_8(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**8
      ci(2)=8.0d0*dx1**7
      ci(3)=28.0d0*dx1**6
      ci(4)=56.0d0*dx1**5
      ci(5)=70.0d0*dx1**4
      ci(6)=56.0d0*dx1**3
      ci(7)=28.0d0*dx1**2
      ci(8)=8.0d0*dx1
      ci(9)=1.0d0
    case(1)
      ci(1)=dx1**8*dx2
      ci(2)=dx1**7*(8.0d0*dx2+dx1)
      ci(3)=4.0d0*dx1**6*(7.0d0*dx2+2.0d0*dx1)
      ci(4)=28.0d0*dx1**5*(2.0d0*dx2+dx1)
      ci(5)=14.0d0*dx1**4*(5.0d0*dx2+4.0d0*dx1)
      ci(6)=14.0d0*dx1**3*(4.0d0*dx2+5.0d0*dx1)
      ci(7)=28.0d0*dx1**2*(dx2+2.0d0*dx1)
      ci(8)=4.0d0*dx1*(2.0d0*dx2+7.0d0*dx1)
      ci(9)=dx2+8.0d0*dx1
      ci(10)=1.0d0
    case(2)
      ci(1)=dx1**8*dx2**2
      ci(2)=2.0d0*dx1**7*dx2*(4.0d0*dx2+dx1)
      ci(3)=dx1**6*(28.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))
      ci(4)=8.0d0*dx1**5*(7.0d0*dx2**2+dx1*(7.0d0*dx2+dx1))
      ci(5)=14.0d0*dx1**4*(5.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+dx1))
      ci(6)=28.0d0*dx1**3*(2.0d0*dx2**2+dx1*(5.0d0*dx2+2.0d0*dx1))
      ci(7)=14.0d0*dx1**2*(2.0d0*dx2**2+dx1*(8.0d0*dx2+5.0d0*dx1))
      ci(8)=8.0d0*dx1*(dx2**2+7.0d0*dx1*(dx2+dx1))
      ci(9)=dx2**2+4.0d0*dx1*(4.0d0*dx2+7.0d0*dx1)
      ci(10)=2.0d0*(dx2+4.0d0*dx1)
      ci(11)=1.0d0
    case(3)
      ci(1)=dx1**8*dx2**3
      ci(2)=dx1**7*dx2**2*(8.0d0*dx2+3.0d0*dx1)
      ci(3)=dx1**6*dx2*(28.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+dx1))
      ci(4)=dx1**5*(56.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(24.0d0*dx2+dx1)))
      ci(5)=2.0d0*dx1**4*(35.0d0*dx2**3+2.0d0*dx1*(42.0d0*dx2**2+dx1*(21.0d0*dx2+2.0d0*dx1)))
      ci(6)=14.0d0*dx1**3*(4.0d0*dx2**3+dx1*(15.0d0*dx2**2+2.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(7)=14.0d0*dx1**2*(2.0d0*dx2**3+dx1*(12.0d0*dx2**2+dx1*(15.0d0*dx2+4.0d0*dx1)))
      ci(8)=2.0d0*dx1*(4.0d0*dx2**3+7.0d0*dx1*(6.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1)))
      ci(9)=dx2**3+4.0d0*dx1*(6.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(10)=3.0d0*dx2**2+4.0d0*dx1*(6.0d0*dx2+7.0d0*dx1)
      ci(11)=3.0d0*dx2+8.0d0*dx1
      ci(12)=1.0d0
    case(4)
      ci(1)=dx1**8*dx2**4
      ci(2)=4.0d0*dx1**7*dx2**3*(2.0d0*dx2+dx1)
      ci(3)=2.0d0*dx1**6*dx2**2*(14.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))
      ci(4)=4.0d0*dx1**5*dx2*(14.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(5)=dx1**4*(70.0d0*dx2**4+dx1*(224.0d0*dx2**3+dx1*(168.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))
      ci(6)=8.0d0*dx1**3*(7.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(14.0d0*dx2+dx1))))
      ci(7)=28.0d0*dx1**2*(dx2**4+dx1*(8.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(8)=8.0d0*dx1*(dx2**4+7.0d0*dx1*(2.0d0*dx2**3+dx1*(6.0d0*dx2**2+dx1*(5.0d0*dx2+dx1))))
      ci(9)=dx2**4+2.0d0*dx1*(16.0d0*dx2**3+7.0d0*dx1*(12.0d0*dx2**2+dx1*(16.0d0*dx2+5.0d0*dx1)))
      ci(10)=4.0d0*(dx2**3+2.0d0*dx1*(6.0d0*dx2**2+7.0d0*dx1*(2.0d0*dx2+dx1)))
      ci(11)=2.0d0*(3.0d0*dx2**2+2.0d0*dx1*(8.0d0*dx2+7.0d0*dx1))
      ci(12)=4.0d0*(dx2+2.0d0*dx1)
      ci(13)=1.0d0
    case(5)
      ci(1)=dx1**8*dx2**5
      ci(2)=dx1**7*dx2**4*(8.0d0*dx2+5.0d0*dx1)
      ci(3)=2.0d0*dx1**6*dx2**3*(14.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))
      ci(4)=2.0d0*dx1**5*dx2**2*(28.0d0*dx2**3+5.0d0*dx1*(14.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**4*dx2*(14.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(56.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(6)=dx1**3*(56.0d0*dx2**5+dx1*(350.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(280.0d0*dx2**2+dx1*(40.0d0*dx2+dx1)))))
      ci(7)=4.0d0*dx1**2*(7.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(175.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1)))))
      ci(8)=4.0d0*dx1*(2.0d0*dx2**5+7.0d0*dx1*(5.0d0*dx2**4+dx1*(20.0d0*dx2**3+dx1*(25.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(9)=dx2**5+2.0d0*dx1*(20.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+dx1*(40.0d0*dx2**2+dx1*(25.0d0*dx2+4.0d0*dx1))))
      ci(10)=5.0d0*(dx2**4+2.0d0*dx1*(8.0d0*dx2**3+7.0d0*dx1*(4.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))))
      ci(11)=2.0d0*(5.0d0*dx2**3+2.0d0*dx1*(20.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)))
      ci(12)=2.0d0*(5.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+7.0d0*dx1))
      ci(13)=5.0d0*dx2+8.0d0*dx1
      ci(14)=1.0d0
    case(6)
      ci(1)=dx1**8*dx2**6
      ci(2)=2.0d0*dx1**7*dx2**5*(4.0d0*dx2+3.0d0*dx1)
      ci(3)=dx1**6*dx2**4*(28.0d0*dx2**2+3.0d0*dx1*(16.0d0*dx2+5.0d0*dx1))
      ci(4)=4.0d0*dx1**5*dx2**3*(14.0d0*dx2**3+dx1*(42.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(5)=dx1**4*dx2**2*(70.0d0*dx2**4+dx1*(336.0d0*dx2**3+5.0d0*dx1*(84.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1))))
      ci(6)=2.0d0*dx1**3*dx2*(28.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(280.0d0*dx2**2+3.0d0*dx1*(20.0d0*dx2+dx1) &
))))
      ci(7)=dx1**2*(28.0d0*dx2**6+dx1*(336.0d0*dx2**5+dx1*(1050.0d0*dx2**4+dx1*(1120.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1* &
(48.0d0*dx2+dx1))))))
      ci(8)=8.0d0*dx1*(dx2**6+dx1*(21.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(175.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)) &
))))
      ci(9)=dx2**6+2.0d0*dx1*(24.0d0*dx2**5+7.0d0*dx1*(30.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(75.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2 &
+dx1)))))
      ci(10)=2.0d0*(3.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1))) &
))
      ci(11)=15.0d0*dx2**4+2.0d0*dx1*(80.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1)))
      ci(12)=4.0d0*(5.0d0*dx2**3+2.0d0*dx1*(15.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(13)=15.0d0*dx2**2+4.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)
      ci(14)=2.0d0*(3.0d0*dx2+4.0d0*dx1)
      ci(15)=1.0d0
    case(7)
      ci(1)=dx1**8*dx2**7
      ci(2)=dx1**7*dx2**6*(8.0d0*dx2+7.0d0*dx1)
      ci(3)=7.0d0*dx1**6*dx2**5*(4.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(4)=7.0d0*dx1**5*dx2**4*(8.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1)))
      ci(5)=7.0d0*dx1**4*dx2**3*(10.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(84.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**3*dx2**2*(8.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1 &
)))))
      ci(7)=7.0d0*dx1**2*dx2*(4.0d0*dx2**6+dx1*(56.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1* &
(24.0d0*dx2+dx1))))))
      ci(8)=dx1*(8.0d0*dx2**7+dx1*(196.0d0*dx2**6+dx1*(1176.0d0*dx2**5+dx1*(2450.0d0*dx2**4+dx1*(1960.0d0*dx2**3+dx1* &
(588.0d0*dx2**2+dx1*(56.0d0*dx2+dx1)))))))
      ci(9)=dx2**7+2.0d0*dx1*(28.0d0*dx2**6+dx1*(294.0d0*dx2**5+dx1*(980.0d0*dx2**4+dx1*(1225.0d0*dx2**3+2.0d0*dx1*(294.0d0*dx2**2 &
+dx1*(49.0d0*dx2+2.0d0*dx1))))))
      ci(10)=7.0d0*(dx2**6+2.0d0*dx1*(12.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(105.0d0*dx2**2+2.0d0*dx1* &
(14.0d0*dx2+dx1))))))
      ci(11)=7.0d0*(3.0d0*dx2**5+2.0d0*dx1*(20.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(35.0d0*dx2+4.0d0*dx1)))))
      ci(12)=7.0d0*(5.0d0*dx2**4+2.0d0*dx1*(20.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(28.0d0*dx2+5.0d0*dx1))))
      ci(13)=7.0d0*(5.0d0*dx2**3+4.0d0*dx1*(6.0d0*dx2**2+dx1*(7.0d0*dx2+2.0d0*dx1)))
      ci(14)=7.0d0*(3.0d0*dx2**2+4.0d0*dx1*(2.0d0*dx2+dx1))
      ci(15)=7.0d0*dx2+8.0d0*dx1
      ci(16)=1.0d0
    case(8)
      ci(1)=dx1**8*dx2**8
      ci(2)=8.0d0*dx1**7*dx2**7*(dx2+dx1)
      ci(3)=4.0d0*dx1**6*dx2**6*(7.0d0*dx2**2+dx1*(16.0d0*dx2+7.0d0*dx1))
      ci(4)=56.0d0*dx1**5*dx2**5*(dx2**3+dx1*(4.0d0*dx2**2+dx1*(4.0d0*dx2+dx1)))
      ci(5)=14.0d0*dx1**4*dx2**4*(5.0d0*dx2**4+dx1*(32.0d0*dx2**3+dx1*(56.0d0*dx2**2+dx1*(32.0d0*dx2+5.0d0*dx1))))
      ci(6)=56.0d0*dx1**3*dx2**3*(dx2**5+dx1*(10.0d0*dx2**4+dx1*(28.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(7)=28.0d0*dx1**2*dx2**2*(dx2**6+dx1*(16.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(112.0d0*dx2**3+dx1*(70.0d0*dx2**2+dx1* &
(16.0d0*dx2+dx1))))))
      ci(8)=8.0d0*dx1*dx2*(dx2**7+dx1*(28.0d0*dx2**6+dx1*(196.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(490.0d0*dx2**3+dx1* &
(196.0d0*dx2**2+dx1*(28.0d0*dx2+dx1)))))))
      ci(9)=dx2**8+dx1*(64.0d0*dx2**7+dx1*(784.0d0*dx2**6+dx1*(3136.0d0*dx2**5+dx1*(4900.0d0*dx2**4+dx1*(3136.0d0*dx2**3+dx1* &
(784.0d0*dx2**2+dx1*(64.0d0*dx2+dx1)))))))
      ci(10)=8.0d0*(dx2**7+dx1*(28.0d0*dx2**6+dx1*(196.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(490.0d0*dx2**3+dx1*(196.0d0*dx2**2+dx1* &
(28.0d0*dx2+dx1)))))))
      ci(11)=28.0d0*(dx2**6+dx1*(16.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(112.0d0*dx2**3+dx1*(70.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))))
      ci(12)=56.0d0*(dx2**5+dx1*(10.0d0*dx2**4+dx1*(28.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(13)=14.0d0*(5.0d0*dx2**4+dx1*(32.0d0*dx2**3+dx1*(56.0d0*dx2**2+dx1*(32.0d0*dx2+5.0d0*dx1))))
      ci(14)=56.0d0*(dx2**3+dx1*(4.0d0*dx2**2+dx1*(4.0d0*dx2+dx1)))
      ci(15)=4.0d0*(7.0d0*dx2**2+dx1*(16.0d0*dx2+7.0d0*dx1))
      ci(16)=8.0d0*(dx2+dx1)
      ci(17)=1.0d0
    case(9)
      ci(1)=dx1**8*dx2**9
      ci(2)=dx1**7*dx2**8*(8.0d0*dx2+9.0d0*dx1)
      ci(3)=4.0d0*dx1**6*dx2**7*(7.0d0*dx2**2+9.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=4.0d0*dx1**5*dx2**6*(14.0d0*dx2**3+3.0d0*dx1*(21.0d0*dx2**2+dx1*(24.0d0*dx2+7.0d0*dx1)))
      ci(5)=14.0d0*dx1**4*dx2**5*(5.0d0*dx2**4+3.0d0*dx1*(12.0d0*dx2**3+dx1*(24.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(6)=14.0d0*dx1**3*dx2**4*(4.0d0*dx2**5+3.0d0*dx1*(15.0d0*dx2**4+dx1*(48.0d0*dx2**3+dx1*(56.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2 &
+dx1)))))
      ci(7)=28.0d0*dx1**2*dx2**3*(dx2**6+3.0d0*dx1*(6.0d0*dx2**5+dx1*(30.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1* &
(12.0d0*dx2+dx1))))))
      ci(8)=4.0d0*dx1*dx2**2*(2.0d0*dx2**7+3.0d0*dx1*(21.0d0*dx2**6+dx1*(168.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(588.0d0*dx2**3 &
+dx1*(294.0d0*dx2**2+dx1*(56.0d0*dx2+3.0d0*dx1)))))))
      ci(9)=dx2*(dx2**8+3.0d0*dx1*(24.0d0*dx2**7+dx1*(336.0d0*dx2**6+dx1*(1568.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(2352.0d0*dx2**3+dx1*(784.0d0*dx2**2+3.0d0*dx1*(32.0d0*dx2+dx1))))))))
      ci(10)=9.0d0*dx2**8+dx1*(288.0d0*dx2**7+dx1*(2352.0d0*dx2**6+dx1*(7056.0d0*dx2**5+dx1*(8820.0d0*dx2**4+dx1*(4704.0d0*dx2**3 &
+dx1*(1008.0d0*dx2**2+dx1*(72.0d0*dx2+dx1)))))))
      ci(11)=4.0d0*(9.0d0*dx2**7+dx1*(168.0d0*dx2**6+dx1*(882.0d0*dx2**5+dx1*(1764.0d0*dx2**4+dx1*(1470.0d0*dx2**3+dx1* &
(504.0d0*dx2**2+dx1*(63.0d0*dx2+2.0d0*dx1)))))))
      ci(12)=28.0d0*(3.0d0*dx2**6+dx1*(36.0d0*dx2**5+dx1*(126.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(18.0d0*dx2 &
+dx1))))))
      ci(13)=14.0d0*(9.0d0*dx2**5+dx1*(72.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(144.0d0*dx2**2+dx1*(45.0d0*dx2+4.0d0*dx1)))))
      ci(14)=14.0d0*(9.0d0*dx2**4+dx1*(48.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))))
      ci(15)=4.0d0*(21.0d0*dx2**3+dx1*(72.0d0*dx2**2+7.0d0*dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(16)=4.0d0*(9.0d0*dx2**2+dx1*(18.0d0*dx2+7.0d0*dx1))
      ci(17)=9.0d0*dx2+8.0d0*dx1
      ci(18)=1.0d0
    case(10)
      ci(1)=dx1**8*dx2**10
      ci(2)=2.0d0*dx1**7*dx2**9*(4.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**6*dx2**8*(28.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+9.0d0*dx1))
      ci(4)=8.0d0*dx1**5*dx2**7*(7.0d0*dx2**3+5.0d0*dx1*(7.0d0*dx2**2+3.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(5)=10.0d0*dx1**4*dx2**6*(7.0d0*dx2**4+dx1*(56.0d0*dx2**3+3.0d0*dx1*(42.0d0*dx2**2+dx1*(32.0d0*dx2+7.0d0*dx1))))
      ci(6)=28.0d0*dx1**3*dx2**5*(2.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(30.0d0*dx2**3+dx1*(40.0d0*dx2**2+dx1*(20.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=14.0d0*dx1**2*dx2**4*(2.0d0*dx2**6+dx1*(40.0d0*dx2**5+3.0d0*dx1*(75.0d0*dx2**4+dx1*(160.0d0*dx2**3+dx1*(140.0d0*dx2**2 &
+dx1*(48.0d0*dx2+5.0d0*dx1))))))
      ci(8)=8.0d0*dx1*dx2**3*(dx2**7+dx1*(35.0d0*dx2**6+3.0d0*dx1*(105.0d0*dx2**5+dx1*(350.0d0*dx2**4+dx1*(490.0d0*dx2**3+dx1* &
(294.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+dx1)))))))
      ci(9)=dx2**2*(dx2**8+dx1*(80.0d0*dx2**7+3.0d0*dx1*(420.0d0*dx2**6+dx1*(2240.0d0*dx2**5+dx1*(4900.0d0*dx2**4+dx1* &
(4704.0d0*dx2**3+5.0d0*dx1*(392.0d0*dx2**2+dx1*(64.0d0*dx2+3.0d0*dx1))))))))
      ci(10)=10.0d0*dx2*(dx2**8+dx1*(36.0d0*dx2**7+dx1*(336.0d0*dx2**6+dx1*(1176.0d0*dx2**5+dx1*(1764.0d0*dx2**4+dx1* &
(1176.0d0*dx2**3+dx1*(336.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))))))
      ci(11)=45.0d0*dx2**8+dx1*(960.0d0*dx2**7+dx1*(5880.0d0*dx2**6+dx1*(14112.0d0*dx2**5+dx1*(14700.0d0*dx2**4+dx1* &
(6720.0d0*dx2**3+dx1*(1260.0d0*dx2**2+dx1*(80.0d0*dx2+dx1)))))))
      ci(12)=8.0d0*(15.0d0*dx2**7+dx1*(210.0d0*dx2**6+dx1*(882.0d0*dx2**5+dx1*(1470.0d0*dx2**4+dx1*(1050.0d0*dx2**3+dx1* &
(315.0d0*dx2**2+dx1*(35.0d0*dx2+dx1)))))))
      ci(13)=14.0d0*(15.0d0*dx2**6+dx1*(144.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1*(480.0d0*dx2**3+dx1*(225.0d0*dx2**2+2.0d0*dx1* &
(20.0d0*dx2+dx1))))))
      ci(14)=28.0d0*(9.0d0*dx2**5+dx1*(60.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(25.0d0*dx2+2.0d0*dx1)))))
      ci(15)=10.0d0*(21.0d0*dx2**4+dx1*(96.0d0*dx2**3+7.0d0*dx1*(18.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(16)=8.0d0*(15.0d0*dx2**3+dx1*(45.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+dx1)))
      ci(17)=45.0d0*dx2**2+4.0d0*dx1*(20.0d0*dx2+7.0d0*dx1)
      ci(18)=2.0d0*(5.0d0*dx2+4.0d0*dx1)
      ci(19)=1.0d0
    case(11)
      ci(1)=dx1**8*dx2**11
      ci(2)=dx1**7*dx2**10*(8.0d0*dx2+11.0d0*dx1)
      ci(3)=dx1**6*dx2**9*(28.0d0*dx2**2+11.0d0*dx1*(8.0d0*dx2+5.0d0*dx1))
      ci(4)=dx1**5*dx2**8*(56.0d0*dx2**3+11.0d0*dx1*(28.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+3.0d0*dx1)))
      ci(5)=2.0d0*dx1**4*dx2**7*(35.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+5.0d0*dx1*(14.0d0*dx2**2+3.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(6)=2.0d0*dx1**3*dx2**6*(28.0d0*dx2**5+11.0d0*dx1*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+3.0d0*dx1*(70.0d0*dx2**2+dx1* &
(40.0d0*dx2+7.0d0*dx1)))))
      ci(7)=14.0d0*dx1**2*dx2**5*(2.0d0*dx2**6+11.0d0*dx1*(4.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(20.0d0*dx2**3+dx1* &
(20.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))))
      ci(8)=2.0d0*dx1*dx2**4*(4.0d0*dx2**7+11.0d0*dx1*(14.0d0*dx2**6+dx1*(140.0d0*dx2**5+3.0d0*dx1*(175.0d0*dx2**4+dx1* &
(280.0d0*dx2**3+dx1*(196.0d0*dx2**2+dx1*(56.0d0*dx2+5.0d0*dx1)))))))
      ci(9)=dx2**3*(dx2**8+11.0d0*dx1*(8.0d0*dx2**7+dx1*(140.0d0*dx2**6+3.0d0*dx1*(280.0d0*dx2**5+dx1*(700.0d0*dx2**4+dx1* &
(784.0d0*dx2**3+dx1*(392.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+dx1))))))))
      ci(10)=11.0d0*dx2**2*(dx2**8+dx1*(40.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(1680.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(2352.0d0*dx2**3+5.0d0*dx1*(168.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))))))
      ci(11)=11.0d0*dx2*(5.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(840.0d0*dx2**6+dx1*(2352.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(1680.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))))))
      ci(12)=165.0d0*dx2**8+dx1*(2640.0d0*dx2**7+dx1*(12936.0d0*dx2**6+dx1*(25872.0d0*dx2**5+dx1*(23100.0d0*dx2**4+dx1* &
(9240.0d0*dx2**3+dx1*(1540.0d0*dx2**2+dx1*(88.0d0*dx2+dx1)))))))
      ci(13)=2.0d0*(165.0d0*dx2**7+dx1*(1848.0d0*dx2**6+dx1*(6468.0d0*dx2**5+dx1*(9240.0d0*dx2**4+dx1*(5775.0d0*dx2**3+2.0d0*dx1* &
(770.0d0*dx2**2+dx1*(77.0d0*dx2+2.0d0*dx1)))))))
      ci(14)=14.0d0*(33.0d0*dx2**6+dx1*(264.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(660.0d0*dx2**3+dx1*(275.0d0*dx2**2+2.0d0*dx1* &
(22.0d0*dx2+dx1))))))
      ci(15)=2.0d0*(231.0d0*dx2**5+dx1*(1320.0d0*dx2**4+7.0d0*dx1*(330.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(55.0d0*dx2+4.0d0*dx1))) &
))
      ci(16)=2.0d0*(165.0d0*dx2**4+dx1*(660.0d0*dx2**3+7.0d0*dx1*(110.0d0*dx2**2+dx1*(44.0d0*dx2+5.0d0*dx1))))
      ci(17)=165.0d0*dx2**3+4.0d0*dx1*(110.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+2.0d0*dx1))
      ci(18)=55.0d0*dx2**2+4.0d0*dx1*(22.0d0*dx2+7.0d0*dx1)
      ci(19)=11.0d0*dx2+8.0d0*dx1
      ci(20)=1.0d0
    case(12)
      ci(1)=dx1**8*dx2**12
      ci(2)=4.0d0*dx1**7*dx2**11*(2.0d0*dx2+3.0d0*dx1)
      ci(3)=2.0d0*dx1**6*dx2**10*(14.0d0*dx2**2+3.0d0*dx1*(16.0d0*dx2+11.0d0*dx1))
      ci(4)=4.0d0*dx1**5*dx2**9*(14.0d0*dx2**3+dx1*(84.0d0*dx2**2+11.0d0*dx1*(12.0d0*dx2+5.0d0*dx1)))
      ci(5)=dx1**4*dx2**8*(70.0d0*dx2**4+dx1*(672.0d0*dx2**3+11.0d0*dx1*(168.0d0*dx2**2+5.0d0*dx1*(32.0d0*dx2+9.0d0*dx1))))
      ci(6)=8.0d0*dx1**3*dx2**7*(7.0d0*dx2**5+dx1*(105.0d0*dx2**4+11.0d0*dx1*(42.0d0*dx2**3+dx1*(70.0d0*dx2**2+9.0d0*dx1* &
(5.0d0*dx2+dx1)))))
      ci(7)=4.0d0*dx1**2*dx2**6*(7.0d0*dx2**6+dx1*(168.0d0*dx2**5+11.0d0*dx1*(105.0d0*dx2**4+dx1*(280.0d0*dx2**3+3.0d0*dx1* &
(105.0d0*dx2**2+dx1*(48.0d0*dx2+7.0d0*dx1))))))
      ci(8)=8.0d0*dx1*dx2**5*(dx2**7+dx1*(42.0d0*dx2**6+11.0d0*dx1*(42.0d0*dx2**5+dx1*(175.0d0*dx2**4+3.0d0*dx1*(105.0d0*dx2**3 &
+dx1*(84.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1)))))))
      ci(9)=dx2**4*(dx2**8+dx1*(96.0d0*dx2**7+11.0d0*dx1*(168.0d0*dx2**6+dx1*(1120.0d0*dx2**5+3.0d0*dx1*(1050.0d0*dx2**4+dx1* &
(1344.0d0*dx2**3+dx1*(784.0d0*dx2**2+3.0d0*dx1*(64.0d0*dx2+5.0d0*dx1))))))))
      ci(10)=4.0d0*dx2**3*(3.0d0*dx2**8+11.0d0*dx1*(12.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(630.0d0*dx2**5+dx1*(1260.0d0*dx2**4 &
+dx1*(1176.0d0*dx2**3+dx1*(504.0d0*dx2**2+5.0d0*dx1*(18.0d0*dx2+dx1))))))))
      ci(11)=22.0d0*dx2**2*(3.0d0*dx2**8+dx1*(80.0d0*dx2**7+dx1*(630.0d0*dx2**6+dx1*(2016.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(2016.0d0*dx2**3+dx1*(630.0d0*dx2**2+dx1*(80.0d0*dx2+3.0d0*dx1))))))))
      ci(12)=4.0d0*dx2*(55.0d0*dx2**8+dx1*(990.0d0*dx2**7+dx1*(5544.0d0*dx2**6+dx1*(12936.0d0*dx2**5+dx1*(13860.0d0*dx2**4+dx1* &
(6930.0d0*dx2**3+dx1*(1540.0d0*dx2**2+3.0d0*dx1*(44.0d0*dx2+dx1))))))))
      ci(13)=495.0d0*dx2**8+dx1*(6336.0d0*dx2**7+dx1*(25872.0d0*dx2**6+dx1*(44352.0d0*dx2**5+dx1*(34650.0d0*dx2**4+dx1* &
(12320.0d0*dx2**3+dx1*(1848.0d0*dx2**2+dx1*(96.0d0*dx2+dx1)))))))
      ci(14)=8.0d0*(99.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(3465.0d0*dx2**4+dx1*(1925.0d0*dx2**3+dx1* &
(462.0d0*dx2**2+dx1*(42.0d0*dx2+dx1)))))))
      ci(15)=4.0d0*(231.0d0*dx2**6+dx1*(1584.0d0*dx2**5+7.0d0*dx1*(495.0d0*dx2**4+dx1*(440.0d0*dx2**3+dx1*(165.0d0*dx2**2+dx1* &
(24.0d0*dx2+dx1))))))
      ci(16)=8.0d0*(99.0d0*dx2**5+dx1*(495.0d0*dx2**4+7.0d0*dx1*(110.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(17)=495.0d0*dx2**4+2.0d0*dx1*(880.0d0*dx2**3+7.0d0*dx1*(132.0d0*dx2**2+dx1*(48.0d0*dx2+5.0d0*dx1)))
      ci(18)=4.0d0*(55.0d0*dx2**3+2.0d0*dx1*(66.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(19)=2.0d0*(33.0d0*dx2**2+2.0d0*dx1*(24.0d0*dx2+7.0d0*dx1))
      ci(20)=4.0d0*(3.0d0*dx2+2.0d0*dx1)
      ci(21)=1.0d0
    case(13)
      ci(1)=dx1**8*dx2**13
      ci(2)=dx1**7*dx2**12*(8.0d0*dx2+13.0d0*dx1)
      ci(3)=2.0d0*dx1**6*dx2**11*(14.0d0*dx2**2+13.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(4)=2.0d0*dx1**5*dx2**10*(28.0d0*dx2**3+13.0d0*dx1*(14.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1)))
      ci(5)=dx1**4*dx2**9*(70.0d0*dx2**4+13.0d0*dx1*(56.0d0*dx2**3+dx1*(168.0d0*dx2**2+11.0d0*dx1*(16.0d0*dx2+5.0d0*dx1))))
      ci(6)=dx1**3*dx2**8*(56.0d0*dx2**5+13.0d0*dx1*(70.0d0*dx2**4+dx1*(336.0d0*dx2**3+11.0d0*dx1*(56.0d0*dx2**2+dx1*(40.0d0*dx2 &
+9.0d0*dx1)))))
      ci(7)=4.0d0*dx1**2*dx2**7*(7.0d0*dx2**6+13.0d0*dx1*(14.0d0*dx2**5+dx1*(105.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+dx1* &
(35.0d0*dx2**2+3.0d0*dx1*(6.0d0*dx2+dx1))))))
      ci(8)=4.0d0*dx1*dx2**6*(2.0d0*dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(84.0d0*dx2**5+11.0d0*dx1*(35.0d0*dx2**4+dx1* &
(70.0d0*dx2**3+3.0d0*dx1*(21.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))))))
      ci(9)=dx2**5*(dx2**8+13.0d0*dx1*(8.0d0*dx2**7+dx1*(168.0d0*dx2**6+11.0d0*dx1*(112.0d0*dx2**5+dx1*(350.0d0*dx2**4+3.0d0*dx1* &
(168.0d0*dx2**3+dx1*(112.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1))))))))
      ci(10)=13.0d0*dx2**4*(dx2**8+dx1*(48.0d0*dx2**7+11.0d0*dx1*(56.0d0*dx2**6+dx1*(280.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1* &
(672.0d0*dx2**3+dx1*(336.0d0*dx2**2+dx1*(72.0d0*dx2+5.0d0*dx1))))))))
      ci(11)=26.0d0*dx2**3*(3.0d0*dx2**8+11.0d0*dx1*(8.0d0*dx2**7+dx1*(70.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1* &
(336.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))))))
      ci(12)=26.0d0*dx2**2*(11.0d0*dx2**8+dx1*(220.0d0*dx2**7+dx1*(1386.0d0*dx2**6+dx1*(3696.0d0*dx2**5+dx1*(4620.0d0*dx2**4+dx1* &
(2772.0d0*dx2**3+dx1*(770.0d0*dx2**2+dx1*(88.0d0*dx2+3.0d0*dx1))))))))
      ci(13)=13.0d0*dx2*(55.0d0*dx2**8+dx1*(792.0d0*dx2**7+dx1*(3696.0d0*dx2**6+dx1*(7392.0d0*dx2**5+dx1*(6930.0d0*dx2**4+dx1* &
(3080.0d0*dx2**3+dx1*(616.0d0*dx2**2+dx1*(48.0d0*dx2+dx1))))))))
      ci(14)=1287.0d0*dx2**8+dx1*(13728.0d0*dx2**7+dx1*(48048.0d0*dx2**6+dx1*(72072.0d0*dx2**5+dx1*(50050.0d0*dx2**4+dx1* &
(16016.0d0*dx2**3+dx1*(2184.0d0*dx2**2+dx1*(104.0d0*dx2+dx1)))))))
      ci(15)=4.0d0*(429.0d0*dx2**7+dx1*(3432.0d0*dx2**6+dx1*(9009.0d0*dx2**5+dx1*(10010.0d0*dx2**4+dx1*(5005.0d0*dx2**3+dx1* &
(1092.0d0*dx2**2+dx1*(91.0d0*dx2+2.0d0*dx1)))))))
      ci(16)=4.0d0*(429.0d0*dx2**6+dx1*(2574.0d0*dx2**5+7.0d0*dx1*(715.0d0*dx2**4+dx1*(572.0d0*dx2**3+dx1*(195.0d0*dx2**2+dx1* &
(26.0d0*dx2+dx1))))))
      ci(17)=1287.0d0*dx2**5+2.0d0*dx1*(2860.0d0*dx2**4+7.0d0*dx1*(572.0d0*dx2**3+dx1*(312.0d0*dx2**2+dx1*(65.0d0*dx2+4.0d0*dx1))) &
)
      ci(18)=715.0d0*dx2**4+2.0d0*dx1*(1144.0d0*dx2**3+7.0d0*dx1*(156.0d0*dx2**2+dx1*(52.0d0*dx2+5.0d0*dx1)))
      ci(19)=2.0d0*(143.0d0*dx2**3+2.0d0*dx1*(156.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+2.0d0*dx1)))
      ci(20)=2.0d0*(39.0d0*dx2**2+2.0d0*dx1*(26.0d0*dx2+7.0d0*dx1))
      ci(21)=13.0d0*dx2+8.0d0*dx1
      ci(22)=1.0d0
    case(14)
      ci(1)=dx1**8*dx2**14
      ci(2)=2.0d0*dx1**7*dx2**13*(4.0d0*dx2+7.0d0*dx1)
      ci(3)=7.0d0*dx1**6*dx2**12*(4.0d0*dx2**2+dx1*(16.0d0*dx2+13.0d0*dx1))
      ci(4)=28.0d0*dx1**5*dx2**11*(2.0d0*dx2**3+dx1*(14.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+dx1)))
      ci(5)=7.0d0*dx1**4*dx2**10*(10.0d0*dx2**4+dx1*(112.0d0*dx2**3+13.0d0*dx1*(28.0d0*dx2**2+dx1*(32.0d0*dx2+11.0d0*dx1))))
      ci(6)=14.0d0*dx1**3*dx2**9*(4.0d0*dx2**5+dx1*(70.0d0*dx2**4+13.0d0*dx1*(28.0d0*dx2**3+dx1*(56.0d0*dx2**2+11.0d0*dx1* &
(4.0d0*dx2+dx1)))))
      ci(7)=7.0d0*dx1**2*dx2**8*(4.0d0*dx2**6+dx1*(112.0d0*dx2**5+13.0d0*dx1*(70.0d0*dx2**4+dx1*(224.0d0*dx2**3+11.0d0*dx1* &
(28.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))))))
      ci(8)=8.0d0*dx1*dx2**7*(dx2**7+dx1*(49.0d0*dx2**6+13.0d0*dx1*(49.0d0*dx2**5+dx1*(245.0d0*dx2**4+11.0d0*dx1*(49.0d0*dx2**3 &
+dx1*(49.0d0*dx2**2+3.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(9)=dx2**6*(dx2**8+dx1*(112.0d0*dx2**7+13.0d0*dx1*(196.0d0*dx2**6+dx1*(1568.0d0*dx2**5+11.0d0*dx1*(490.0d0*dx2**4+dx1* &
(784.0d0*dx2**3+3.0d0*dx1*(196.0d0*dx2**2+dx1*(64.0d0*dx2+7.0d0*dx1))))))))
      ci(10)=14.0d0*dx2**5*(dx2**8+13.0d0*dx1*(4.0d0*dx2**7+dx1*(56.0d0*dx2**6+11.0d0*dx1*(28.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1* &
(84.0d0*dx2**3+dx1*(48.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))))))
      ci(11)=91.0d0*dx2**4*(dx2**8+dx1*(32.0d0*dx2**7+11.0d0*dx1*(28.0d0*dx2**6+dx1*(112.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1* &
(192.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))))))
      ci(12)=364.0d0*dx2**3*(dx2**8+dx1*(22.0d0*dx2**7+dx1*(154.0d0*dx2**6+dx1*(462.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1* &
(462.0d0*dx2**3+dx1*(154.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))))))))
      ci(13)=91.0d0*dx2**2*(11.0d0*dx2**8+dx1*(176.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2112.0d0*dx2**5+dx1*(2310.0d0*dx2**4+dx1* &
(1232.0d0*dx2**3+dx1*(308.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))))))
      ci(14)=14.0d0*dx2*(143.0d0*dx2**8+dx1*(1716.0d0*dx2**7+dx1*(6864.0d0*dx2**6+dx1*(12012.0d0*dx2**5+dx1*(10010.0d0*dx2**4+dx1* &
(4004.0d0*dx2**3+dx1*(728.0d0*dx2**2+dx1*(52.0d0*dx2+dx1))))))))
      ci(15)=3003.0d0*dx2**8+dx1*(27456.0d0*dx2**7+dx1*(84084.0d0*dx2**6+dx1*(112112.0d0*dx2**5+dx1*(70070.0d0*dx2**4+dx1* &
(20384.0d0*dx2**3+dx1*(2548.0d0*dx2**2+dx1*(112.0d0*dx2+dx1)))))))
      ci(16)=8.0d0*(429.0d0*dx2**7+dx1*(3003.0d0*dx2**6+dx1*(7007.0d0*dx2**5+dx1*(7007.0d0*dx2**4+dx1*(3185.0d0*dx2**3+dx1* &
(637.0d0*dx2**2+dx1*(49.0d0*dx2+dx1)))))))
      ci(17)=7.0d0*(429.0d0*dx2**6+2.0d0*dx1*(1144.0d0*dx2**5+dx1*(2002.0d0*dx2**4+dx1*(1456.0d0*dx2**3+dx1*(455.0d0*dx2**2 &
+2.0d0*dx1*(28.0d0*dx2+dx1))))))
      ci(18)=14.0d0*(143.0d0*dx2**5+2.0d0*dx1*(286.0d0*dx2**4+dx1*(364.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1))) &
))
      ci(19)=7.0d0*(143.0d0*dx2**4+2.0d0*dx1*(208.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(56.0d0*dx2+5.0d0*dx1))))
      ci(20)=28.0d0*(13.0d0*dx2**3+2.0d0*dx1*(13.0d0*dx2**2+dx1*(7.0d0*dx2+dx1)))
      ci(21)=7.0d0*(13.0d0*dx2**2+4.0d0*dx1*(4.0d0*dx2+dx1))
      ci(22)=2.0d0*(7.0d0*dx2+4.0d0*dx1)
      ci(23)=1.0d0
    case(15)
      ci(1)=dx1**8*dx2**15
      ci(2)=dx1**7*dx2**14*(8.0d0*dx2+15.0d0*dx1)
      ci(3)=dx1**6*dx2**13*(28.0d0*dx2**2+15.0d0*dx1*(8.0d0*dx2+7.0d0*dx1))
      ci(4)=7.0d0*dx1**5*dx2**12*(8.0d0*dx2**3+5.0d0*dx1*(12.0d0*dx2**2+dx1*(24.0d0*dx2+13.0d0*dx1)))
      ci(5)=35.0d0*dx1**4*dx2**11*(2.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(84.0d0*dx2**2+13.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(6)=7.0d0*dx1**3*dx2**10*(8.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(840.0d0*dx2**3+13.0d0*dx1*(140.0d0*dx2**2+3.0d0*dx1* &
(40.0d0*dx2+11.0d0*dx1)))))
      ci(7)=7.0d0*dx1**2*dx2**9*(4.0d0*dx2**6+dx1*(120.0d0*dx2**5+dx1*(1050.0d0*dx2**4+13.0d0*dx1*(280.0d0*dx2**3+dx1* &
(420.0d0*dx2**2+11.0d0*dx1*(24.0d0*dx2+5.0d0*dx1))))))
      ci(8)=dx1*dx2**8*(8.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(5880.0d0*dx2**5+13.0d0*dx1*(2450.0d0*dx2**4+dx1*(5880.0d0*dx2**3 &
+11.0d0*dx1*(588.0d0*dx2**2+5.0d0*dx1*(56.0d0*dx2+9.0d0*dx1)))))))
      ci(9)=dx2**7*(dx2**8+dx1*(120.0d0*dx2**7+dx1*(2940.0d0*dx2**6+13.0d0*dx1*(1960.0d0*dx2**5+dx1*(7350.0d0*dx2**4+11.0d0*dx1* &
(1176.0d0*dx2**3+5.0d0*dx1*(196.0d0*dx2**2+9.0d0*dx1*(8.0d0*dx2+dx1))))))))
      ci(10)=5.0d0*dx2**6*(3.0d0*dx2**8+dx1*(168.0d0*dx2**7+13.0d0*dx1*(196.0d0*dx2**6+dx1*(1176.0d0*dx2**5+11.0d0*dx1* &
(294.0d0*dx2**4+dx1*(392.0d0*dx2**3+dx1*(252.0d0*dx2**2+dx1*(72.0d0*dx2+7.0d0*dx1))))))))
      ci(11)=7.0d0*dx2**5*(15.0d0*dx2**8+13.0d0*dx1*(40.0d0*dx2**7+dx1*(420.0d0*dx2**6+11.0d0*dx1*(168.0d0*dx2**5+dx1* &
(350.0d0*dx2**4+dx1*(360.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))))))))
      ci(12)=91.0d0*dx2**4*(5.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(3080.0d0*dx2**5+dx1*(4950.0d0*dx2**4+dx1* &
(3960.0d0*dx2**3+dx1*(1540.0d0*dx2**2+3.0d0*dx1*(88.0d0*dx2+5.0d0*dx1))))))))
      ci(13)=91.0d0*dx2**3*(15.0d0*dx2**8+dx1*(264.0d0*dx2**7+dx1*(1540.0d0*dx2**6+dx1*(3960.0d0*dx2**5+dx1*(4950.0d0*dx2**4+dx1* &
(3080.0d0*dx2**3+dx1*(924.0d0*dx2**2+5.0d0*dx1*(24.0d0*dx2+dx1))))))))
      ci(14)=7.0d0*dx2**2*(429.0d0*dx2**8+dx1*(5720.0d0*dx2**7+dx1*(25740.0d0*dx2**6+dx1*(51480.0d0*dx2**5+dx1*(50050.0d0*dx2**4 &
+dx1*(24024.0d0*dx2**3+5.0d0*dx1*(1092.0d0*dx2**2+dx1*(104.0d0*dx2+3.0d0*dx1))))))))
      ci(15)=5.0d0*dx2*(1001.0d0*dx2**8+dx1*(10296.0d0*dx2**7+dx1*(36036.0d0*dx2**6+dx1*(56056.0d0*dx2**5+dx1*(42042.0d0*dx2**4 &
+dx1*(15288.0d0*dx2**3+dx1*(2548.0d0*dx2**2+3.0d0*dx1*(56.0d0*dx2+dx1))))))))
      ci(16)=6435.0d0*dx2**8+dx1*(51480.0d0*dx2**7+dx1*(140140.0d0*dx2**6+dx1*(168168.0d0*dx2**5+dx1*(95550.0d0*dx2**4+dx1* &
(25480.0d0*dx2**3+dx1*(2940.0d0*dx2**2+dx1*(120.0d0*dx2+dx1)))))))
      ci(17)=6435.0d0*dx2**7+2.0d0*dx1*(20020.0d0*dx2**6+dx1*(42042.0d0*dx2**5+dx1*(38220.0d0*dx2**4+dx1*(15925.0d0*dx2**3 &
+2.0d0*dx1*(1470.0d0*dx2**2+dx1*(105.0d0*dx2+2.0d0*dx1))))))
      ci(18)=7.0d0*(715.0d0*dx2**6+2.0d0*dx1*(1716.0d0*dx2**5+dx1*(2730.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1*(525.0d0*dx2**2 &
+2.0d0*dx1*(30.0d0*dx2+dx1))))))
      ci(19)=7.0d0*(429.0d0*dx2**5+2.0d0*dx1*(780.0d0*dx2**4+dx1*(910.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(75.0d0*dx2+4.0d0*dx1)))) &
)
      ci(20)=35.0d0*(39.0d0*dx2**4+2.0d0*dx1*(52.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(21)=7.0d0*(65.0d0*dx2**3+4.0d0*dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(22)=105.0d0*dx2**2+4.0d0*dx1*(30.0d0*dx2+7.0d0*dx1)
      ci(23)=15.0d0*dx2+8.0d0*dx1
      ci(24)=1.0d0
    case(16)
      ci(1)=dx1**8*dx2**16
      ci(2)=8.0d0*dx1**7*dx2**15*(dx2+2.0d0*dx1)
      ci(3)=4.0d0*dx1**6*dx2**14*(7.0d0*dx2**2+2.0d0*dx1*(16.0d0*dx2+15.0d0*dx1))
      ci(4)=8.0d0*dx1**5*dx2**13*(7.0d0*dx2**3+2.0d0*dx1*(28.0d0*dx2**2+5.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)))
      ci(5)=14.0d0*dx1**4*dx2**12*(5.0d0*dx2**4+2.0d0*dx1*(32.0d0*dx2**3+5.0d0*dx1*(24.0d0*dx2**2+dx1*(32.0d0*dx2+13.0d0*dx1))))
      ci(6)=56.0d0*dx1**3*dx2**11*(dx2**5+2.0d0*dx1*(10.0d0*dx2**4+dx1*(60.0d0*dx2**3+dx1*(140.0d0*dx2**2+13.0d0*dx1*(10.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=28.0d0*dx1**2*dx2**10*(dx2**6+2.0d0*dx1*(16.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(560.0d0*dx2**3+13.0d0*dx1* &
(70.0d0*dx2**2+dx1*(48.0d0*dx2+11.0d0*dx1))))))
      ci(8)=8.0d0*dx1*dx2**9*(dx2**7+2.0d0*dx1*(28.0d0*dx2**6+dx1*(420.0d0*dx2**5+dx1*(2450.0d0*dx2**4+13.0d0*dx1*(490.0d0*dx2**3 &
+dx1*(588.0d0*dx2**2+11.0d0*dx1*(28.0d0*dx2+5.0d0*dx1)))))))
      ci(9)=dx2**8*(dx2**8+2.0d0*dx1*(64.0d0*dx2**7+dx1*(1680.0d0*dx2**6+dx1*(15680.0d0*dx2**5+13.0d0*dx1*(4900.0d0*dx2**4+dx1* &
(9408.0d0*dx2**3+11.0d0*dx1*(784.0d0*dx2**2+5.0d0*dx1*(64.0d0*dx2+9.0d0*dx1))))))))
      ci(10)=16.0d0*dx2**7*(dx2**8+dx1*(60.0d0*dx2**7+dx1*(980.0d0*dx2**6+13.0d0*dx1*(490.0d0*dx2**5+dx1*(1470.0d0*dx2**4 &
+11.0d0*dx1*(196.0d0*dx2**3+5.0d0*dx1*(28.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))))))))
      ci(11)=8.0d0*dx2**6*(15.0d0*dx2**8+dx1*(560.0d0*dx2**7+13.0d0*dx1*(490.0d0*dx2**6+dx1*(2352.0d0*dx2**5+11.0d0*dx1* &
(490.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(315.0d0*dx2**2+dx1*(80.0d0*dx2+7.0d0*dx1))))))))
      ci(12)=112.0d0*dx2**5*(5.0d0*dx2**8+13.0d0*dx1*(10.0d0*dx2**7+dx1*(84.0d0*dx2**6+dx1*(308.0d0*dx2**5+dx1*(550.0d0*dx2**4 &
+dx1*(495.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(44.0d0*dx2+3.0d0*dx1))))))))
      ci(13)=364.0d0*dx2**4*(5.0d0*dx2**8+dx1*(96.0d0*dx2**7+dx1*(616.0d0*dx2**6+dx1*(1760.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1* &
(1760.0d0*dx2**3+dx1*(616.0d0*dx2**2+dx1*(96.0d0*dx2+5.0d0*dx1))))))))
      ci(14)=112.0d0*dx2**3*(39.0d0*dx2**8+dx1*(572.0d0*dx2**7+dx1*(2860.0d0*dx2**6+dx1*(6435.0d0*dx2**5+dx1*(7150.0d0*dx2**4+dx1* &
(4004.0d0*dx2**3+dx1*(1092.0d0*dx2**2+5.0d0*dx1*(26.0d0*dx2+dx1))))))))
      ci(15)=8.0d0*dx2**2*(1001.0d0*dx2**8+dx1*(11440.0d0*dx2**7+dx1*(45045.0d0*dx2**6+dx1*(80080.0d0*dx2**5+dx1*(70070.0d0*dx2**4 &
+dx1*(30576.0d0*dx2**3+5.0d0*dx1*(1274.0d0*dx2**2+dx1*(112.0d0*dx2+3.0d0*dx1))))))))
      ci(16)=16.0d0*dx2*(715.0d0*dx2**8+dx1*(6435.0d0*dx2**7+dx1*(20020.0d0*dx2**6+dx1*(28028.0d0*dx2**5+dx1*(19110.0d0*dx2**4 &
+dx1*(6370.0d0*dx2**3+dx1*(980.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))))))
      ci(17)=12870.0d0*dx2**8+dx1*(91520.0d0*dx2**7+dx1*(224224.0d0*dx2**6+dx1*(244608.0d0*dx2**5+dx1*(127400.0d0*dx2**4+dx1* &
(31360.0d0*dx2**3+dx1*(3360.0d0*dx2**2+dx1*(128.0d0*dx2+dx1)))))))
      ci(18)=8.0d0*(1430.0d0*dx2**7+dx1*(8008.0d0*dx2**6+dx1*(15288.0d0*dx2**5+dx1*(12740.0d0*dx2**4+dx1*(4900.0d0*dx2**3+dx1* &
(840.0d0*dx2**2+dx1*(56.0d0*dx2+dx1)))))))
      ci(19)=28.0d0*(286.0d0*dx2**6+dx1*(1248.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1120.0d0*dx2**3+dx1*(300.0d0*dx2**2+dx1* &
(32.0d0*dx2+dx1))))))
      ci(20)=56.0d0*(78.0d0*dx2**5+dx1*(260.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)))))
      ci(21)=14.0d0*(130.0d0*dx2**4+dx1*(320.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(64.0d0*dx2+5.0d0*dx1))))
      ci(22)=8.0d0*(70.0d0*dx2**3+dx1*(120.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+dx1)))
      ci(23)=4.0d0*(30.0d0*dx2**2+dx1*(32.0d0*dx2+7.0d0*dx1))
      ci(24)=8.0d0*(2.0d0*dx2+dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>16, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^9 * (x-x2)^n2 as sum_k=0^(9+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_9(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**9
      ci(2)=9.0d0*dx1**8
      ci(3)=36.0d0*dx1**7
      ci(4)=84.0d0*dx1**6
      ci(5)=126.0d0*dx1**5
      ci(6)=126.0d0*dx1**4
      ci(7)=84.0d0*dx1**3
      ci(8)=36.0d0*dx1**2
      ci(9)=9.0d0*dx1
      ci(10)=1.0d0
    case(1)
      ci(1)=dx1**9*dx2
      ci(2)=dx1**8*(9.0d0*dx2+dx1)
      ci(3)=9.0d0*dx1**7*(4.0d0*dx2+dx1)
      ci(4)=12.0d0*dx1**6*(7.0d0*dx2+3.0d0*dx1)
      ci(5)=42.0d0*dx1**5*(3.0d0*dx2+2.0d0*dx1)
      ci(6)=126.0d0*dx1**4*(dx2+dx1)
      ci(7)=42.0d0*dx1**3*(2.0d0*dx2+3.0d0*dx1)
      ci(8)=12.0d0*dx1**2*(3.0d0*dx2+7.0d0*dx1)
      ci(9)=9.0d0*dx1*(dx2+4.0d0*dx1)
      ci(10)=dx2+9.0d0*dx1
      ci(11)=1.0d0
    case(2)
      ci(1)=dx1**9*dx2**2
      ci(2)=dx1**8*dx2*(9.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**7*(36.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))
      ci(4)=3.0d0*dx1**6*(28.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+dx1))
      ci(5)=6.0d0*dx1**5*(21.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+3.0d0*dx1))
      ci(6)=42.0d0*dx1**4*(3.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))
      ci(7)=42.0d0*dx1**3*(2.0d0*dx2**2+3.0d0*dx1*(2.0d0*dx2+dx1))
      ci(8)=6.0d0*dx1**2*(6.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(9)=3.0d0*dx1*(3.0d0*dx2**2+4.0d0*dx1*(6.0d0*dx2+7.0d0*dx1))
      ci(10)=dx2**2+18.0d0*dx1*(dx2+2.0d0*dx1)
      ci(11)=2.0d0*dx2+9.0d0*dx1
      ci(12)=1.0d0
    case(3)
      ci(1)=dx1**9*dx2**3
      ci(2)=3.0d0*dx1**8*dx2**2*(3.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**7*dx2*(12.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))
      ci(4)=dx1**6*(84.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1*(27.0d0*dx2+dx1)))
      ci(5)=9.0d0*dx1**5*(14.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(6)=18.0d0*dx1**4*(7.0d0*dx2**3+dx1*(21.0d0*dx2**2+2.0d0*dx1*(7.0d0*dx2+dx1)))
      ci(7)=42.0d0*dx1**3*(2.0d0*dx2**3+dx1*(9.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(8)=18.0d0*dx1**2*(2.0d0*dx2**3+7.0d0*dx1*(2.0d0*dx2**2+dx1*(3.0d0*dx2+dx1)))
      ci(9)=9.0d0*dx1*(dx2**3+2.0d0*dx1*(6.0d0*dx2**2+7.0d0*dx1*(2.0d0*dx2+dx1)))
      ci(10)=dx2**3+3.0d0*dx1*(9.0d0*dx2**2+4.0d0*dx1*(9.0d0*dx2+7.0d0*dx1))
      ci(11)=3.0d0*(dx2**2+3.0d0*dx1*(3.0d0*dx2+4.0d0*dx1))
      ci(12)=3.0d0*(dx2+3.0d0*dx1)
      ci(13)=1.0d0
    case(4)
      ci(1)=dx1**9*dx2**4
      ci(2)=dx1**8*dx2**3*(9.0d0*dx2+4.0d0*dx1)
      ci(3)=6.0d0*dx1**7*dx2**2*(6.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(4)=2.0d0*dx1**6*dx2*(42.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(27.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**5*(126.0d0*dx2**4+dx1*(336.0d0*dx2**3+dx1*(216.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))
      ci(6)=9.0d0*dx1**4*(14.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(56.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(7)=12.0d0*dx1**3*(7.0d0*dx2**4+dx1*(42.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1))))
      ci(8)=12.0d0*dx1**2*(3.0d0*dx2**4+7.0d0*dx1*(4.0d0*dx2**3+dx1*(9.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))
      ci(9)=9.0d0*dx1*(dx2**4+2.0d0*dx1*(8.0d0*dx2**3+7.0d0*dx1*(4.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))))
      ci(10)=dx2**4+6.0d0*dx1*(6.0d0*dx2**3+dx1*(36.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+3.0d0*dx1)))
      ci(11)=2.0d0*(2.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)))
      ci(12)=6.0d0*(dx2**2+6.0d0*dx1*(dx2+dx1))
      ci(13)=4.0d0*dx2+9.0d0*dx1
      ci(14)=1.0d0
    case(5)
      ci(1)=dx1**9*dx2**5
      ci(2)=dx1**8*dx2**4*(9.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**7*dx2**3*(36.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+2.0d0*dx1))
      ci(4)=2.0d0*dx1**6*dx2**2*(42.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(5)=dx1**5*dx2*(126.0d0*dx2**4+5.0d0*dx1*(84.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))))
      ci(6)=dx1**4*(126.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1*(840.0d0*dx2**3+dx1*(360.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))))
      ci(7)=3.0d0*dx1**3*(28.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(280.0d0*dx2**2+3.0d0*dx1*(20.0d0*dx2+dx1)))))
      ci(8)=12.0d0*dx1**2*(3.0d0*dx2**5+dx1*(35.0d0*dx2**4+dx1*(105.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1)))))
      ci(9)=3.0d0*dx1*(3.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1 &
)))))
      ci(10)=dx2**5+3.0d0*dx1*(15.0d0*dx2**4+2.0d0*dx1*(60.0d0*dx2**3+7.0d0*dx1*(20.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+dx1))))
      ci(11)=5.0d0*dx2**4+6.0d0*dx1*(15.0d0*dx2**3+dx1*(60.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+3.0d0*dx1)))
      ci(12)=2.0d0*(5.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+7.0d0*dx1)))
      ci(13)=10.0d0*dx2**2+9.0d0*dx1*(5.0d0*dx2+4.0d0*dx1)
      ci(14)=5.0d0*dx2+9.0d0*dx1
      ci(15)=1.0d0
    case(6)
      ci(1)=dx1**9*dx2**6
      ci(2)=3.0d0*dx1**8*dx2**5*(3.0d0*dx2+2.0d0*dx1)
      ci(3)=3.0d0*dx1**7*dx2**4*(12.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1))
      ci(4)=dx1**6*dx2**3*(84.0d0*dx2**3+dx1*(216.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+4.0d0*dx1)))
      ci(5)=3.0d0*dx1**5*dx2**2*(42.0d0*dx2**4+dx1*(168.0d0*dx2**3+5.0d0*dx1*(36.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(6)=3.0d0*dx1**4*dx2*(42.0d0*dx2**5+dx1*(252.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1) &
))))
      ci(7)=dx1**3*(84.0d0*dx2**6+dx1*(756.0d0*dx2**5+dx1*(1890.0d0*dx2**4+dx1*(1680.0d0*dx2**3+dx1*(540.0d0*dx2**2+dx1* &
(54.0d0*dx2+dx1))))))
      ci(8)=9.0d0*dx1**2*(4.0d0*dx2**6+dx1*(56.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1* &
(24.0d0*dx2+dx1))))))
      ci(9)=9.0d0*dx1*(dx2**6+2.0d0*dx1*(12.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(105.0d0*dx2**2+2.0d0*dx1* &
(14.0d0*dx2+dx1))))))
      ci(10)=dx2**6+6.0d0*dx1*(9.0d0*dx2**5+dx1*(90.0d0*dx2**4+7.0d0*dx1*(40.0d0*dx2**3+dx1*(45.0d0*dx2**2+2.0d0*dx1*(9.0d0*dx2 &
+dx1)))))
      ci(11)=3.0d0*(2.0d0*dx2**5+3.0d0*dx1*(15.0d0*dx2**4+2.0d0*dx1*(40.0d0*dx2**3+7.0d0*dx1*(10.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))) &
)
      ci(12)=3.0d0*(5.0d0*dx2**4+6.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(13)=20.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+4.0d0*dx1*(18.0d0*dx2+7.0d0*dx1))
      ci(14)=3.0d0*(5.0d0*dx2**2+6.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(15)=3.0d0*(2.0d0*dx2+3.0d0*dx1)
      ci(16)=1.0d0
    case(7)
      ci(1)=dx1**9*dx2**7
      ci(2)=dx1**8*dx2**6*(9.0d0*dx2+7.0d0*dx1)
      ci(3)=3.0d0*dx1**7*dx2**5*(12.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+dx1))
      ci(4)=7.0d0*dx1**6*dx2**4*(12.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(27.0d0*dx2+5.0d0*dx1)))
      ci(5)=7.0d0*dx1**5*dx2**3*(18.0d0*dx2**4+dx1*(84.0d0*dx2**3+dx1*(108.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+dx1))))
      ci(6)=21.0d0*dx1**4*dx2**2*(6.0d0*dx2**5+dx1*(42.0d0*dx2**4+dx1*(84.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(7)=7.0d0*dx1**3*dx2*(12.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(378.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1* &
(27.0d0*dx2+dx1))))))
      ci(8)=dx1**2*(36.0d0*dx2**7+dx1*(588.0d0*dx2**6+dx1*(2646.0d0*dx2**5+dx1*(4410.0d0*dx2**4+dx1*(2940.0d0*dx2**3+dx1* &
(756.0d0*dx2**2+dx1*(63.0d0*dx2+dx1)))))))
      ci(9)=9.0d0*dx1*(dx2**7+dx1*(28.0d0*dx2**6+dx1*(196.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(490.0d0*dx2**3+dx1*(196.0d0*dx2**2 &
+dx1*(28.0d0*dx2+dx1)))))))
      ci(10)=dx2**7+3.0d0*dx1*(21.0d0*dx2**6+2.0d0*dx1*(126.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(735.0d0*dx2**3+dx1*(441.0d0*dx2**2 &
+2.0d0*dx1*(49.0d0*dx2+3.0d0*dx1))))))
      ci(11)=7.0d0*(dx2**6+3.0d0*dx1*(9.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(21.0d0*dx2 &
+2.0d0*dx1))))))
      ci(12)=21.0d0*(dx2**5+3.0d0*dx1*(5.0d0*dx2**4+2.0d0*dx1*(10.0d0*dx2**3+dx1*(14.0d0*dx2**2+dx1*(7.0d0*dx2+dx1)))))
      ci(13)=7.0d0*(5.0d0*dx2**4+3.0d0*dx1*(15.0d0*dx2**3+2.0d0*dx1*(18.0d0*dx2**2+dx1*(14.0d0*dx2+3.0d0*dx1))))
      ci(14)=7.0d0*(5.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+4.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(15)=3.0d0*(7.0d0*dx2**2+3.0d0*dx1*(7.0d0*dx2+4.0d0*dx1))
      ci(16)=7.0d0*dx2+9.0d0*dx1
      ci(17)=1.0d0
    case(8)
      ci(1)=dx1**9*dx2**8
      ci(2)=dx1**8*dx2**7*(9.0d0*dx2+8.0d0*dx1)
      ci(3)=4.0d0*dx1**7*dx2**6*(9.0d0*dx2**2+dx1*(18.0d0*dx2+7.0d0*dx1))
      ci(4)=4.0d0*dx1**6*dx2**5*(21.0d0*dx2**3+dx1*(72.0d0*dx2**2+7.0d0*dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(5)=14.0d0*dx1**5*dx2**4*(9.0d0*dx2**4+dx1*(48.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))))
      ci(6)=14.0d0*dx1**4*dx2**3*(9.0d0*dx2**5+dx1*(72.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(144.0d0*dx2**2+dx1*(45.0d0*dx2 &
+4.0d0*dx1)))))
      ci(7)=28.0d0*dx1**3*dx2**2*(3.0d0*dx2**6+dx1*(36.0d0*dx2**5+dx1*(126.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1* &
(18.0d0*dx2+dx1))))))
      ci(8)=4.0d0*dx1**2*dx2*(9.0d0*dx2**7+dx1*(168.0d0*dx2**6+dx1*(882.0d0*dx2**5+dx1*(1764.0d0*dx2**4+dx1*(1470.0d0*dx2**3+dx1* &
(504.0d0*dx2**2+dx1*(63.0d0*dx2+2.0d0*dx1)))))))
      ci(9)=dx1*(9.0d0*dx2**8+dx1*(288.0d0*dx2**7+dx1*(2352.0d0*dx2**6+dx1*(7056.0d0*dx2**5+dx1*(8820.0d0*dx2**4+dx1* &
(4704.0d0*dx2**3+dx1*(1008.0d0*dx2**2+dx1*(72.0d0*dx2+dx1))))))))
      ci(10)=dx2**8+3.0d0*dx1*(24.0d0*dx2**7+dx1*(336.0d0*dx2**6+dx1*(1568.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1*(2352.0d0*dx2**3 &
+dx1*(784.0d0*dx2**2+3.0d0*dx1*(32.0d0*dx2+dx1)))))))
      ci(11)=4.0d0*(2.0d0*dx2**7+3.0d0*dx1*(21.0d0*dx2**6+dx1*(168.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(588.0d0*dx2**3+dx1* &
(294.0d0*dx2**2+dx1*(56.0d0*dx2+3.0d0*dx1)))))))
      ci(12)=28.0d0*(dx2**6+3.0d0*dx1*(6.0d0*dx2**5+dx1*(30.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)) &
))))
      ci(13)=14.0d0*(4.0d0*dx2**5+3.0d0*dx1*(15.0d0*dx2**4+dx1*(48.0d0*dx2**3+dx1*(56.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+dx1)))))
      ci(14)=14.0d0*(5.0d0*dx2**4+3.0d0*dx1*(12.0d0*dx2**3+dx1*(24.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(15)=4.0d0*(14.0d0*dx2**3+3.0d0*dx1*(21.0d0*dx2**2+dx1*(24.0d0*dx2+7.0d0*dx1)))
      ci(16)=4.0d0*(7.0d0*dx2**2+9.0d0*dx1*(2.0d0*dx2+dx1))
      ci(17)=8.0d0*dx2+9.0d0*dx1
      ci(18)=1.0d0
    case(9)
      ci(1)=dx1**9*dx2**9
      ci(2)=9.0d0*dx1**8*dx2**8*(dx2+dx1)
      ci(3)=9.0d0*dx1**7*dx2**7*(4.0d0*dx2**2+dx1*(9.0d0*dx2+4.0d0*dx1))
      ci(4)=12.0d0*dx1**6*dx2**6*(7.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(27.0d0*dx2+7.0d0*dx1)))
      ci(5)=18.0d0*dx1**5*dx2**5*(7.0d0*dx2**4+dx1*(42.0d0*dx2**3+dx1*(72.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+dx1))))
      ci(6)=126.0d0*dx1**4*dx2**4*(dx2**5+dx1*(9.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(24.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))))
      ci(7)=42.0d0*dx1**3*dx2**3*(2.0d0*dx2**6+dx1*(27.0d0*dx2**5+dx1*(108.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1* &
(27.0d0*dx2+2.0d0*dx1))))))
      ci(8)=36.0d0*dx1**2*dx2**2*(dx2**7+dx1*(21.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(294.0d0*dx2**4+dx1*(294.0d0*dx2**3+dx1* &
(126.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)))))))
      ci(9)=9.0d0*dx1*dx2*(dx2**8+dx1*(36.0d0*dx2**7+dx1*(336.0d0*dx2**6+dx1*(1176.0d0*dx2**5+dx1*(1764.0d0*dx2**4+dx1* &
(1176.0d0*dx2**3+dx1*(336.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))))))
      ci(10)=dx2**9+dx1*(81.0d0*dx2**8+dx1*(1296.0d0*dx2**7+dx1*(7056.0d0*dx2**6+dx1*(15876.0d0*dx2**5+dx1*(15876.0d0*dx2**4+dx1* &
(7056.0d0*dx2**3+dx1*(1296.0d0*dx2**2+dx1*(81.0d0*dx2+dx1))))))))
      ci(11)=9.0d0*(dx2**8+dx1*(36.0d0*dx2**7+dx1*(336.0d0*dx2**6+dx1*(1176.0d0*dx2**5+dx1*(1764.0d0*dx2**4+dx1*(1176.0d0*dx2**3 &
+dx1*(336.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))))))
      ci(12)=36.0d0*(dx2**7+dx1*(21.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(294.0d0*dx2**4+dx1*(294.0d0*dx2**3+dx1*(126.0d0*dx2**2 &
+dx1*(21.0d0*dx2+dx1)))))))
      ci(13)=42.0d0*(2.0d0*dx2**6+dx1*(27.0d0*dx2**5+dx1*(108.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1*(27.0d0*dx2 &
+2.0d0*dx1))))))
      ci(14)=126.0d0*(dx2**5+dx1*(9.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(24.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))))
      ci(15)=18.0d0*(7.0d0*dx2**4+dx1*(42.0d0*dx2**3+dx1*(72.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+dx1))))
      ci(16)=12.0d0*(7.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(27.0d0*dx2+7.0d0*dx1)))
      ci(17)=9.0d0*(4.0d0*dx2**2+dx1*(9.0d0*dx2+4.0d0*dx1))
      ci(18)=9.0d0*(dx2+dx1)
      ci(19)=1.0d0
    case(10)
      ci(1)=dx1**9*dx2**10
      ci(2)=dx1**8*dx2**9*(9.0d0*dx2+10.0d0*dx1)
      ci(3)=9.0d0*dx1**7*dx2**8*(4.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=3.0d0*dx1**6*dx2**7*(28.0d0*dx2**3+5.0d0*dx1*(24.0d0*dx2**2+dx1*(27.0d0*dx2+8.0d0*dx1)))
      ci(5)=6.0d0*dx1**5*dx2**6*(21.0d0*dx2**4+5.0d0*dx1*(28.0d0*dx2**3+dx1*(54.0d0*dx2**2+dx1*(36.0d0*dx2+7.0d0*dx1))))
      ci(6)=18.0d0*dx1**4*dx2**5*(7.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(210.0d0*dx2**3+dx1*(240.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=42.0d0*dx1**3*dx2**4*(2.0d0*dx2**6+dx1*(30.0d0*dx2**5+dx1*(135.0d0*dx2**4+dx1*(240.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1* &
(54.0d0*dx2+5.0d0*dx1))))))
      ci(8)=6.0d0*dx1**2*dx2**3*(6.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(945.0d0*dx2**5+dx1*(2520.0d0*dx2**4+dx1*(2940.0d0*dx2**3 &
+dx1*(1512.0d0*dx2**2+5.0d0*dx1*(63.0d0*dx2+4.0d0*dx1)))))))
      ci(9)=9.0d0*dx1*dx2**2*(dx2**8+dx1*(40.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(1680.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(2352.0d0*dx2**3+5.0d0*dx1*(168.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))))))
      ci(10)=dx2*(dx2**9+dx1*(90.0d0*dx2**8+dx1*(1620.0d0*dx2**7+dx1*(10080.0d0*dx2**6+dx1*(26460.0d0*dx2**5+dx1*(31752.0d0*dx2**4 &
+5.0d0*dx1*(3528.0d0*dx2**3+dx1*(864.0d0*dx2**2+dx1*(81.0d0*dx2+2.0d0*dx1)))))))))
      ci(11)=10.0d0*dx2**9+dx1*(405.0d0*dx2**8+dx1*(4320.0d0*dx2**7+dx1*(17640.0d0*dx2**6+dx1*(31752.0d0*dx2**5+dx1* &
(26460.0d0*dx2**4+dx1*(10080.0d0*dx2**3+dx1*(1620.0d0*dx2**2+dx1*(90.0d0*dx2+dx1))))))))
      ci(12)=9.0d0*(5.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(840.0d0*dx2**6+dx1*(2352.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(1680.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))))))
      ci(13)=6.0d0*(20.0d0*dx2**7+dx1*(315.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1*(2520.0d0*dx2**3+dx1* &
(945.0d0*dx2**2+2.0d0*dx1*(70.0d0*dx2+3.0d0*dx1)))))))
      ci(14)=42.0d0*(5.0d0*dx2**6+dx1*(54.0d0*dx2**5+dx1*(180.0d0*dx2**4+dx1*(240.0d0*dx2**3+dx1*(135.0d0*dx2**2+2.0d0*dx1* &
(15.0d0*dx2+dx1))))))
      ci(15)=18.0d0*(14.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(240.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(16)=6.0d0*(35.0d0*dx2**4+dx1*(180.0d0*dx2**3+dx1*(270.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+3.0d0*dx1))))
      ci(17)=3.0d0*(40.0d0*dx2**3+dx1*(135.0d0*dx2**2+4.0d0*dx1*(30.0d0*dx2+7.0d0*dx1)))
      ci(18)=9.0d0*(5.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))
      ci(19)=10.0d0*dx2+9.0d0*dx1
      ci(20)=1.0d0
    case(11)
      ci(1)=dx1**9*dx2**11
      ci(2)=dx1**8*dx2**10*(9.0d0*dx2+11.0d0*dx1)
      ci(3)=dx1**7*dx2**9*(36.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+5.0d0*dx1))
      ci(4)=3.0d0*dx1**6*dx2**8*(28.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(5)=3.0d0*dx1**5*dx2**7*(42.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+5.0d0*dx1*(12.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1))))
      ci(6)=6.0d0*dx1**4*dx2**6*(21.0d0*dx2**5+11.0d0*dx1*(21.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(45.0d0*dx2 &
+7.0d0*dx1)))))
      ci(7)=6.0d0*dx1**3*dx2**5*(14.0d0*dx2**6+11.0d0*dx1*(21.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(210.0d0*dx2**3+dx1* &
(180.0d0*dx2**2+7.0d0*dx1*(9.0d0*dx2+dx1))))))
      ci(8)=6.0d0*dx1**2*dx2**4*(6.0d0*dx2**7+11.0d0*dx1*(14.0d0*dx2**6+dx1*(105.0d0*dx2**5+dx1*(315.0d0*dx2**4+dx1* &
(420.0d0*dx2**3+dx1*(252.0d0*dx2**2+dx1*(63.0d0*dx2+5.0d0*dx1)))))))
      ci(9)=3.0d0*dx1*dx2**3*(3.0d0*dx2**8+11.0d0*dx1*(12.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(630.0d0*dx2**5+dx1*(1260.0d0*dx2**4 &
+dx1*(1176.0d0*dx2**3+dx1*(504.0d0*dx2**2+5.0d0*dx1*(18.0d0*dx2+dx1))))))))
      ci(10)=dx2**2*(dx2**9+11.0d0*dx1*(9.0d0*dx2**8+dx1*(180.0d0*dx2**7+dx1*(1260.0d0*dx2**6+dx1*(3780.0d0*dx2**5+dx1* &
(5292.0d0*dx2**4+dx1*(3528.0d0*dx2**3+5.0d0*dx1*(216.0d0*dx2**2+dx1*(27.0d0*dx2+dx1)))))))))
      ci(11)=11.0d0*dx2*(dx2**9+dx1*(45.0d0*dx2**8+dx1*(540.0d0*dx2**7+dx1*(2520.0d0*dx2**6+dx1*(5292.0d0*dx2**5+dx1* &
(5292.0d0*dx2**4+dx1*(2520.0d0*dx2**3+dx1*(540.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))))))))
      ci(12)=55.0d0*dx2**9+dx1*(1485.0d0*dx2**8+dx1*(11880.0d0*dx2**7+dx1*(38808.0d0*dx2**6+dx1*(58212.0d0*dx2**5+dx1* &
(41580.0d0*dx2**4+dx1*(13860.0d0*dx2**3+dx1*(1980.0d0*dx2**2+dx1*(99.0d0*dx2+dx1))))))))
      ci(13)=3.0d0*(55.0d0*dx2**8+dx1*(990.0d0*dx2**7+dx1*(5544.0d0*dx2**6+dx1*(12936.0d0*dx2**5+dx1*(13860.0d0*dx2**4+dx1* &
(6930.0d0*dx2**3+dx1*(1540.0d0*dx2**2+3.0d0*dx1*(44.0d0*dx2+dx1))))))))
      ci(14)=6.0d0*(55.0d0*dx2**7+dx1*(693.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(4620.0d0*dx2**4+dx1*(3465.0d0*dx2**3+dx1* &
(1155.0d0*dx2**2+2.0d0*dx1*(77.0d0*dx2+3.0d0*dx1)))))))
      ci(15)=6.0d0*(77.0d0*dx2**6+dx1*(693.0d0*dx2**5+dx1*(1980.0d0*dx2**4+7.0d0*dx1*(330.0d0*dx2**3+dx1*(165.0d0*dx2**2+dx1* &
(33.0d0*dx2+2.0d0*dx1))))))
      ci(16)=6.0d0*(77.0d0*dx2**5+dx1*(495.0d0*dx2**4+dx1*(990.0d0*dx2**3+7.0d0*dx1*(110.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+dx1)))))
      ci(17)=3.0d0*(110.0d0*dx2**4+dx1*(495.0d0*dx2**3+2.0d0*dx1*(330.0d0*dx2**2+7.0d0*dx1*(22.0d0*dx2+3.0d0*dx1))))
      ci(18)=3.0d0*(55.0d0*dx2**3+dx1*(165.0d0*dx2**2+4.0d0*dx1*(33.0d0*dx2+7.0d0*dx1)))
      ci(19)=55.0d0*dx2**2+9.0d0*dx1*(11.0d0*dx2+4.0d0*dx1)
      ci(20)=11.0d0*dx2+9.0d0*dx1
      ci(21)=1.0d0
    case(12)
      ci(1)=dx1**9*dx2**12
      ci(2)=3.0d0*dx1**8*dx2**11*(3.0d0*dx2+4.0d0*dx1)
      ci(3)=6.0d0*dx1**7*dx2**10*(6.0d0*dx2**2+dx1*(18.0d0*dx2+11.0d0*dx1))
      ci(4)=2.0d0*dx1**6*dx2**9*(42.0d0*dx2**3+dx1*(216.0d0*dx2**2+11.0d0*dx1*(27.0d0*dx2+10.0d0*dx1)))
      ci(5)=9.0d0*dx1**5*dx2**8*(14.0d0*dx2**4+dx1*(112.0d0*dx2**3+11.0d0*dx1*(24.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(6)=9.0d0*dx1**4*dx2**7*(14.0d0*dx2**5+dx1*(168.0d0*dx2**4+11.0d0*dx1*(56.0d0*dx2**3+dx1*(80.0d0*dx2**2+dx1*(45.0d0*dx2 &
+8.0d0*dx1)))))
      ci(7)=12.0d0*dx1**3*dx2**6*(7.0d0*dx2**6+dx1*(126.0d0*dx2**5+11.0d0*dx1*(63.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1* &
(135.0d0*dx2**2+dx1*(54.0d0*dx2+7.0d0*dx1))))))
      ci(8)=36.0d0*dx1**2*dx2**5*(dx2**7+dx1*(28.0d0*dx2**6+11.0d0*dx1*(21.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(105.0d0*dx2**3+dx1* &
(72.0d0*dx2**2+dx1*(21.0d0*dx2+2.0d0*dx1)))))))
      ci(9)=9.0d0*dx1*dx2**4*(dx2**8+dx1*(48.0d0*dx2**7+11.0d0*dx1*(56.0d0*dx2**6+dx1*(280.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1* &
(672.0d0*dx2**3+dx1*(336.0d0*dx2**2+dx1*(72.0d0*dx2+5.0d0*dx1))))))))
      ci(10)=dx2**3*(dx2**9+dx1*(108.0d0*dx2**8+11.0d0*dx1*(216.0d0*dx2**7+dx1*(1680.0d0*dx2**6+dx1*(5670.0d0*dx2**5+dx1* &
(9072.0d0*dx2**4+dx1*(7056.0d0*dx2**3+dx1*(2592.0d0*dx2**2+5.0d0*dx1*(81.0d0*dx2+4.0d0*dx1)))))))))
      ci(11)=6.0d0*dx2**2*(2.0d0*dx2**9+11.0d0*dx1*(9.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(630.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1* &
(1764.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(270.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))))))
      ci(12)=6.0d0*dx2*(11.0d0*dx2**9+dx1*(330.0d0*dx2**8+dx1*(2970.0d0*dx2**7+dx1*(11088.0d0*dx2**6+dx1*(19404.0d0*dx2**5+dx1* &
(16632.0d0*dx2**4+dx1*(6930.0d0*dx2**3+dx1*(1320.0d0*dx2**2+dx1*(99.0d0*dx2+2.0d0*dx1)))))))))
      ci(13)=220.0d0*dx2**9+dx1*(4455.0d0*dx2**8+dx1*(28512.0d0*dx2**7+dx1*(77616.0d0*dx2**6+dx1*(99792.0d0*dx2**5+dx1* &
(62370.0d0*dx2**4+dx1*(18480.0d0*dx2**3+dx1*(2376.0d0*dx2**2+dx1*(108.0d0*dx2+dx1))))))))
      ci(14)=9.0d0*(55.0d0*dx2**8+dx1*(792.0d0*dx2**7+dx1*(3696.0d0*dx2**6+dx1*(7392.0d0*dx2**5+dx1*(6930.0d0*dx2**4+dx1* &
(3080.0d0*dx2**3+dx1*(616.0d0*dx2**2+dx1*(48.0d0*dx2+dx1))))))))
      ci(15)=36.0d0*(22.0d0*dx2**7+dx1*(231.0d0*dx2**6+dx1*(792.0d0*dx2**5+dx1*(1155.0d0*dx2**4+dx1*(770.0d0*dx2**3+dx1* &
(231.0d0*dx2**2+dx1*(28.0d0*dx2+dx1)))))))
      ci(16)=12.0d0*(77.0d0*dx2**6+dx1*(594.0d0*dx2**5+dx1*(1485.0d0*dx2**4+7.0d0*dx1*(220.0d0*dx2**3+dx1*(99.0d0*dx2**2+dx1* &
(18.0d0*dx2+dx1))))))
      ci(17)=9.0d0*(88.0d0*dx2**5+dx1*(495.0d0*dx2**4+2.0d0*dx1*(440.0d0*dx2**3+7.0d0*dx1*(44.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))))
      ci(18)=9.0d0*(55.0d0*dx2**4+2.0d0*dx1*(110.0d0*dx2**3+dx1*(132.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(19)=2.0d0*(110.0d0*dx2**3+3.0d0*dx1*(99.0d0*dx2**2+2.0d0*dx1*(36.0d0*dx2+7.0d0*dx1)))
      ci(20)=6.0d0*(11.0d0*dx2**2+6.0d0*dx1*(3.0d0*dx2+dx1))
      ci(21)=3.0d0*(4.0d0*dx2+3.0d0*dx1)
      ci(22)=1.0d0
    case(13)
      ci(1)=dx1**9*dx2**13
      ci(2)=dx1**8*dx2**12*(9.0d0*dx2+13.0d0*dx1)
      ci(3)=3.0d0*dx1**7*dx2**11*(12.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(4)=2.0d0*dx1**6*dx2**10*(42.0d0*dx2**3+13.0d0*dx1*(18.0d0*dx2**2+dx1*(27.0d0*dx2+11.0d0*dx1)))
      ci(5)=dx1**5*dx2**9*(126.0d0*dx2**4+13.0d0*dx1*(84.0d0*dx2**3+dx1*(216.0d0*dx2**2+11.0d0*dx1*(18.0d0*dx2+5.0d0*dx1))))
      ci(6)=9.0d0*dx1**4*dx2**8*(14.0d0*dx2**5+13.0d0*dx1*(14.0d0*dx2**4+dx1*(56.0d0*dx2**3+11.0d0*dx1*(8.0d0*dx2**2+dx1* &
(5.0d0*dx2+dx1)))))
      ci(7)=3.0d0*dx1**3*dx2**7*(28.0d0*dx2**6+13.0d0*dx1*(42.0d0*dx2**5+dx1*(252.0d0*dx2**4+11.0d0*dx1*(56.0d0*dx2**3+dx1* &
(60.0d0*dx2**2+dx1*(27.0d0*dx2+4.0d0*dx1))))))
      ci(8)=12.0d0*dx1**2*dx2**6*(3.0d0*dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(63.0d0*dx2**5+11.0d0*dx1*(21.0d0*dx2**4+dx1* &
(35.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))))))
      ci(9)=9.0d0*dx1*dx2**5*(dx2**8+13.0d0*dx1*(4.0d0*dx2**7+dx1*(56.0d0*dx2**6+11.0d0*dx1*(28.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1* &
(84.0d0*dx2**3+dx1*(48.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))))))
      ci(10)=dx2**4*(dx2**9+13.0d0*dx1*(9.0d0*dx2**8+dx1*(216.0d0*dx2**7+11.0d0*dx1*(168.0d0*dx2**6+dx1*(630.0d0*dx2**5+dx1* &
(1134.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(432.0d0*dx2**2+dx1*(81.0d0*dx2+5.0d0*dx1)))))))))
      ci(11)=13.0d0*dx2**3*(dx2**9+dx1*(54.0d0*dx2**8+11.0d0*dx1*(72.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(1134.0d0*dx2**5+dx1* &
(1512.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(324.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))))))))
      ci(12)=78.0d0*dx2**2*(dx2**9+dx1*(33.0d0*dx2**8+dx1*(330.0d0*dx2**7+dx1*(1386.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1* &
(2772.0d0*dx2**4+dx1*(1386.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1*(33.0d0*dx2+dx1)))))))))
      ci(13)=13.0d0*dx2*(22.0d0*dx2**9+dx1*(495.0d0*dx2**8+dx1*(3564.0d0*dx2**7+dx1*(11088.0d0*dx2**6+dx1*(16632.0d0*dx2**5+dx1* &
(12474.0d0*dx2**4+dx1*(4620.0d0*dx2**3+dx1*(792.0d0*dx2**2+dx1*(54.0d0*dx2+dx1)))))))))
      ci(14)=715.0d0*dx2**9+dx1*(11583.0d0*dx2**8+dx1*(61776.0d0*dx2**7+dx1*(144144.0d0*dx2**6+dx1*(162162.0d0*dx2**5+dx1* &
(90090.0d0*dx2**4+dx1*(24024.0d0*dx2**3+dx1*(2808.0d0*dx2**2+dx1*(117.0d0*dx2+dx1))))))))
      ci(15)=9.0d0*(143.0d0*dx2**8+dx1*(1716.0d0*dx2**7+dx1*(6864.0d0*dx2**6+dx1*(12012.0d0*dx2**5+dx1*(10010.0d0*dx2**4+dx1* &
(4004.0d0*dx2**3+dx1*(728.0d0*dx2**2+dx1*(52.0d0*dx2+dx1))))))))
      ci(16)=12.0d0*(143.0d0*dx2**7+dx1*(1287.0d0*dx2**6+dx1*(3861.0d0*dx2**5+dx1*(5005.0d0*dx2**4+dx1*(3003.0d0*dx2**3+dx1* &
(819.0d0*dx2**2+dx1*(91.0d0*dx2+3.0d0*dx1)))))))
      ci(17)=3.0d0*(572.0d0*dx2**6+dx1*(3861.0d0*dx2**5+2.0d0*dx1*(4290.0d0*dx2**4+7.0d0*dx1*(572.0d0*dx2**3+dx1*(234.0d0*dx2**2 &
+dx1*(39.0d0*dx2+2.0d0*dx1))))))
      ci(18)=9.0d0*(143.0d0*dx2**5+dx1*(715.0d0*dx2**4+2.0d0*dx1*(572.0d0*dx2**3+7.0d0*dx1*(52.0d0*dx2**2+dx1*(13.0d0*dx2+dx1)))))
      ci(19)=715.0d0*dx2**4+6.0d0*dx1*(429.0d0*dx2**3+dx1*(468.0d0*dx2**2+7.0d0*dx1*(26.0d0*dx2+3.0d0*dx1)))
      ci(20)=2.0d0*(143.0d0*dx2**3+3.0d0*dx1*(117.0d0*dx2**2+2.0d0*dx1*(39.0d0*dx2+7.0d0*dx1)))
      ci(21)=3.0d0*(26.0d0*dx2**2+3.0d0*dx1*(13.0d0*dx2+4.0d0*dx1))
      ci(22)=13.0d0*dx2+9.0d0*dx1
      ci(23)=1.0d0
    case(14)
      ci(1)=dx1**9*dx2**14
      ci(2)=dx1**8*dx2**13*(9.0d0*dx2+14.0d0*dx1)
      ci(3)=dx1**7*dx2**12*(36.0d0*dx2**2+7.0d0*dx1*(18.0d0*dx2+13.0d0*dx1))
      ci(4)=7.0d0*dx1**6*dx2**11*(12.0d0*dx2**3+dx1*(72.0d0*dx2**2+13.0d0*dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(5)=7.0d0*dx1**5*dx2**10*(18.0d0*dx2**4+dx1*(168.0d0*dx2**3+13.0d0*dx1*(36.0d0*dx2**2+dx1*(36.0d0*dx2+11.0d0*dx1))))
      ci(6)=7.0d0*dx1**4*dx2**9*(18.0d0*dx2**5+dx1*(252.0d0*dx2**4+13.0d0*dx1*(84.0d0*dx2**3+dx1*(144.0d0*dx2**2+11.0d0*dx1* &
(9.0d0*dx2+2.0d0*dx1)))))
      ci(7)=21.0d0*dx1**3*dx2**8*(4.0d0*dx2**6+dx1*(84.0d0*dx2**5+13.0d0*dx1*(42.0d0*dx2**4+dx1*(112.0d0*dx2**3+11.0d0*dx1* &
(12.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))))
      ci(8)=3.0d0*dx1**2*dx2**7*(12.0d0*dx2**7+dx1*(392.0d0*dx2**6+13.0d0*dx1*(294.0d0*dx2**5+dx1*(1176.0d0*dx2**4+11.0d0*dx1* &
(196.0d0*dx2**3+dx1*(168.0d0*dx2**2+dx1*(63.0d0*dx2+8.0d0*dx1)))))))
      ci(9)=3.0d0*dx1*dx2**6*(3.0d0*dx2**8+dx1*(168.0d0*dx2**7+13.0d0*dx1*(196.0d0*dx2**6+dx1*(1176.0d0*dx2**5+11.0d0*dx1* &
(294.0d0*dx2**4+dx1*(392.0d0*dx2**3+dx1*(252.0d0*dx2**2+dx1*(72.0d0*dx2+7.0d0*dx1))))))))
      ci(10)=dx2**5*(dx2**9+dx1*(126.0d0*dx2**8+13.0d0*dx1*(252.0d0*dx2**7+dx1*(2352.0d0*dx2**6+11.0d0*dx1*(882.0d0*dx2**5+dx1* &
(1764.0d0*dx2**4+dx1*(1764.0d0*dx2**3+dx1*(864.0d0*dx2**2+7.0d0*dx1*(27.0d0*dx2+2.0d0*dx1)))))))))
      ci(11)=7.0d0*dx2**4*(2.0d0*dx2**9+13.0d0*dx1*(9.0d0*dx2**8+dx1*(144.0d0*dx2**7+11.0d0*dx1*(84.0d0*dx2**6+dx1*(252.0d0*dx2**5 &
+dx1*(378.0d0*dx2**4+dx1*(288.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))))))))
      ci(12)=91.0d0*dx2**3*(dx2**9+dx1*(36.0d0*dx2**8+dx1*(396.0d0*dx2**7+dx1*(1848.0d0*dx2**6+dx1*(4158.0d0*dx2**5+dx1* &
(4752.0d0*dx2**4+dx1*(2772.0d0*dx2**3+dx1*(792.0d0*dx2**2+dx1*(99.0d0*dx2+4.0d0*dx1)))))))))
      ci(13)=91.0d0*dx2**2*(4.0d0*dx2**9+dx1*(99.0d0*dx2**8+dx1*(792.0d0*dx2**7+dx1*(2772.0d0*dx2**6+dx1*(4752.0d0*dx2**5+dx1* &
(4158.0d0*dx2**4+dx1*(1848.0d0*dx2**3+dx1*(396.0d0*dx2**2+dx1*(36.0d0*dx2+dx1)))))))))
      ci(14)=7.0d0*dx2*(143.0d0*dx2**9+dx1*(2574.0d0*dx2**8+dx1*(15444.0d0*dx2**7+dx1*(41184.0d0*dx2**6+dx1*(54054.0d0*dx2**5+dx1* &
(36036.0d0*dx2**4+dx1*(12012.0d0*dx2**3+dx1*(1872.0d0*dx2**2+dx1*(117.0d0*dx2+2.0d0*dx1)))))))))
      ci(15)=2002.0d0*dx2**9+dx1*(27027.0d0*dx2**8+dx1*(123552.0d0*dx2**7+dx1*(252252.0d0*dx2**6+dx1*(252252.0d0*dx2**5+dx1* &
(126126.0d0*dx2**4+dx1*(30576.0d0*dx2**3+dx1*(3276.0d0*dx2**2+dx1*(126.0d0*dx2+dx1))))))))
      ci(16)=3.0d0*(1001.0d0*dx2**8+dx1*(10296.0d0*dx2**7+dx1*(36036.0d0*dx2**6+dx1*(56056.0d0*dx2**5+dx1*(42042.0d0*dx2**4+dx1* &
(15288.0d0*dx2**3+dx1*(2548.0d0*dx2**2+3.0d0*dx1*(56.0d0*dx2+dx1))))))))
      ci(17)=3.0d0*(1144.0d0*dx2**7+dx1*(9009.0d0*dx2**6+2.0d0*dx1*(12012.0d0*dx2**5+dx1*(14014.0d0*dx2**4+dx1*(7644.0d0*dx2**3 &
+dx1*(1911.0d0*dx2**2+2.0d0*dx1*(98.0d0*dx2+3.0d0*dx1)))))))
      ci(18)=21.0d0*(143.0d0*dx2**6+2.0d0*dx1*(429.0d0*dx2**5+dx1*(858.0d0*dx2**4+dx1*(728.0d0*dx2**3+dx1*(273.0d0*dx2**2 &
+2.0d0*dx1*(21.0d0*dx2+dx1))))))
      ci(19)=7.0d0*(286.0d0*dx2**5+3.0d0*dx1*(429.0d0*dx2**4+2.0d0*dx1*(312.0d0*dx2**3+dx1*(182.0d0*dx2**2+3.0d0*dx1*(14.0d0*dx2 &
+dx1)))))
      ci(20)=7.0d0*(143.0d0*dx2**4+6.0d0*dx1*(78.0d0*dx2**3+dx1*(78.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1))))
      ci(21)=7.0d0*(52.0d0*dx2**3+3.0d0*dx1*(39.0d0*dx2**2+4.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(22)=91.0d0*dx2**2+18.0d0*dx1*(7.0d0*dx2+2.0d0*dx1)
      ci(23)=14.0d0*dx2+9.0d0*dx1
      ci(24)=1.0d0
    case(15)
      ci(1)=dx1**9*dx2**15
      ci(2)=3.0d0*dx1**8*dx2**14*(3.0d0*dx2+5.0d0*dx1)
      ci(3)=3.0d0*dx1**7*dx2**13*(12.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+7.0d0*dx1))
      ci(4)=dx1**6*dx2**12*(84.0d0*dx2**3+5.0d0*dx1*(108.0d0*dx2**2+7.0d0*dx1*(27.0d0*dx2+13.0d0*dx1)))
      ci(5)=21.0d0*dx1**5*dx2**11*(6.0d0*dx2**4+5.0d0*dx1*(12.0d0*dx2**3+dx1*(36.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+dx1))))
      ci(6)=21.0d0*dx1**4*dx2**10*(6.0d0*dx2**5+dx1*(90.0d0*dx2**4+dx1*(420.0d0*dx2**3+13.0d0*dx1*(60.0d0*dx2**2+dx1*(45.0d0*dx2 &
+11.0d0*dx1)))))
      ci(7)=7.0d0*dx1**3*dx2**9*(12.0d0*dx2**6+dx1*(270.0d0*dx2**5+dx1*(1890.0d0*dx2**4+13.0d0*dx1*(420.0d0*dx2**3+dx1* &
(540.0d0*dx2**2+11.0d0*dx1*(27.0d0*dx2+5.0d0*dx1))))))
      ci(8)=9.0d0*dx1**2*dx2**8*(4.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(1470.0d0*dx2**5+13.0d0*dx1*(490.0d0*dx2**4+dx1* &
(980.0d0*dx2**3+11.0d0*dx1*(84.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(9)=9.0d0*dx1*dx2**7*(dx2**8+dx1*(60.0d0*dx2**7+dx1*(980.0d0*dx2**6+13.0d0*dx1*(490.0d0*dx2**5+dx1*(1470.0d0*dx2**4 &
+11.0d0*dx1*(196.0d0*dx2**3+5.0d0*dx1*(28.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))))))))
      ci(10)=dx2**6*(dx2**9+dx1*(135.0d0*dx2**8+dx1*(3780.0d0*dx2**7+13.0d0*dx1*(2940.0d0*dx2**6+dx1*(13230.0d0*dx2**5+11.0d0*dx1* &
(2646.0d0*dx2**4+5.0d0*dx1*(588.0d0*dx2**3+dx1*(324.0d0*dx2**2+dx1*(81.0d0*dx2+7.0d0*dx1)))))))))
      ci(11)=3.0d0*dx2**5*(5.0d0*dx2**9+dx1*(315.0d0*dx2**8+13.0d0*dx1*(420.0d0*dx2**7+dx1*(2940.0d0*dx2**6+11.0d0*dx1* &
(882.0d0*dx2**5+dx1*(1470.0d0*dx2**4+dx1*(1260.0d0*dx2**3+dx1*(540.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+dx1)))))))))
      ci(12)=21.0d0*dx2**4*(5.0d0*dx2**9+13.0d0*dx1*(15.0d0*dx2**8+dx1*(180.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2310.0d0*dx2**5 &
+dx1*(2970.0d0*dx2**4+dx1*(1980.0d0*dx2**3+dx1*(660.0d0*dx2**2+dx1*(99.0d0*dx2+5.0d0*dx1)))))))))
      ci(13)=91.0d0*dx2**3*(5.0d0*dx2**9+dx1*(135.0d0*dx2**8+dx1*(1188.0d0*dx2**7+dx1*(4620.0d0*dx2**6+dx1*(8910.0d0*dx2**5+dx1* &
(8910.0d0*dx2**4+dx1*(4620.0d0*dx2**3+dx1*(1188.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+dx1)))))))))
      ci(14)=21.0d0*dx2**2*(65.0d0*dx2**9+dx1*(1287.0d0*dx2**8+dx1*(8580.0d0*dx2**7+dx1*(25740.0d0*dx2**6+dx1*(38610.0d0*dx2**5 &
+dx1*(30030.0d0*dx2**4+dx1*(12012.0d0*dx2**3+5.0d0*dx1*(468.0d0*dx2**2+dx1*(39.0d0*dx2+dx1)))))))))
      ci(15)=3.0d0*dx2*(1001.0d0*dx2**9+dx1*(15015.0d0*dx2**8+dx1*(77220.0d0*dx2**7+dx1*(180180.0d0*dx2**6+dx1*(210210.0d0*dx2**5 &
+dx1*(126126.0d0*dx2**4+5.0d0*dx1*(7644.0d0*dx2**3+dx1*(1092.0d0*dx2**2+dx1*(63.0d0*dx2+dx1)))))))))
      ci(16)=5005.0d0*dx2**9+dx1*(57915.0d0*dx2**8+dx1*(231660.0d0*dx2**7+dx1*(420420.0d0*dx2**6+dx1*(378378.0d0*dx2**5+dx1* &
(171990.0d0*dx2**4+dx1*(38220.0d0*dx2**3+dx1*(3780.0d0*dx2**2+dx1*(135.0d0*dx2+dx1))))))))
      ci(17)=9.0d0*(715.0d0*dx2**8+dx1*(6435.0d0*dx2**7+dx1*(20020.0d0*dx2**6+dx1*(28028.0d0*dx2**5+dx1*(19110.0d0*dx2**4+dx1* &
(6370.0d0*dx2**3+dx1*(980.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))))))
      ci(18)=9.0d0*(715.0d0*dx2**7+dx1*(5005.0d0*dx2**6+2.0d0*dx1*(6006.0d0*dx2**5+dx1*(6370.0d0*dx2**4+dx1*(3185.0d0*dx2**3+dx1* &
(735.0d0*dx2**2+2.0d0*dx1*(35.0d0*dx2+dx1)))))))
      ci(19)=7.0d0*(715.0d0*dx2**6+3.0d0*dx1*(1287.0d0*dx2**5+2.0d0*dx1*(1170.0d0*dx2**4+dx1*(910.0d0*dx2**3+dx1*(315.0d0*dx2**2 &
+dx1*(45.0d0*dx2+2.0d0*dx1))))))
      ci(20)=21.0d0*(143.0d0*dx2**5+3.0d0*dx1*(195.0d0*dx2**4+2.0d0*dx1*(130.0d0*dx2**3+dx1*(70.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))) &
)
      ci(21)=21.0d0*(65.0d0*dx2**4+3.0d0*dx1*(65.0d0*dx2**3+2.0d0*dx1*(30.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(22)=455.0d0*dx2**3+3.0d0*dx1*(315.0d0*dx2**2+4.0d0*dx1*(45.0d0*dx2+7.0d0*dx1))
      ci(23)=3.0d0*(35.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+4.0d0*dx1))
      ci(24)=3.0d0*(5.0d0*dx2+3.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>15, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^10 * (x-x2)^n2 as sum_k=0^(10+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_10(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**10
      ci(2)=10.0d0*dx1**9
      ci(3)=45.0d0*dx1**8
      ci(4)=120.0d0*dx1**7
      ci(5)=210.0d0*dx1**6
      ci(6)=252.0d0*dx1**5
      ci(7)=210.0d0*dx1**4
      ci(8)=120.0d0*dx1**3
      ci(9)=45.0d0*dx1**2
      ci(10)=10.0d0*dx1
      ci(11)=1.0d0
    case(1)
      ci(1)=dx1**10*dx2
      ci(2)=dx1**9*(10.0d0*dx2+dx1)
      ci(3)=5.0d0*dx1**8*(9.0d0*dx2+2.0d0*dx1)
      ci(4)=15.0d0*dx1**7*(8.0d0*dx2+3.0d0*dx1)
      ci(5)=30.0d0*dx1**6*(7.0d0*dx2+4.0d0*dx1)
      ci(6)=42.0d0*dx1**5*(6.0d0*dx2+5.0d0*dx1)
      ci(7)=42.0d0*dx1**4*(5.0d0*dx2+6.0d0*dx1)
      ci(8)=30.0d0*dx1**3*(4.0d0*dx2+7.0d0*dx1)
      ci(9)=15.0d0*dx1**2*(3.0d0*dx2+8.0d0*dx1)
      ci(10)=5.0d0*dx1*(2.0d0*dx2+9.0d0*dx1)
      ci(11)=dx2+10.0d0*dx1
      ci(12)=1.0d0
    case(2)
      ci(1)=dx1**10*dx2**2
      ci(2)=2.0d0*dx1**9*dx2*(5.0d0*dx2+dx1)
      ci(3)=dx1**8*(45.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))
      ci(4)=10.0d0*dx1**7*(12.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))
      ci(5)=15.0d0*dx1**6*(14.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))
      ci(6)=12.0d0*dx1**5*(21.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))
      ci(7)=42.0d0*dx1**4*(5.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(8)=12.0d0*dx1**3*(10.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))
      ci(9)=15.0d0*dx1**2*(3.0d0*dx2**2+2.0d0*dx1*(8.0d0*dx2+7.0d0*dx1))
      ci(10)=10.0d0*dx1*(dx2**2+3.0d0*dx1*(3.0d0*dx2+4.0d0*dx1))
      ci(11)=dx2**2+5.0d0*dx1*(4.0d0*dx2+9.0d0*dx1)
      ci(12)=2.0d0*(dx2+5.0d0*dx1)
      ci(13)=1.0d0
    case(3)
      ci(1)=dx1**10*dx2**3
      ci(2)=dx1**9*dx2**2*(10.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**8*dx2*(15.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))
      ci(4)=dx1**7*(120.0d0*dx2**3+dx1*(135.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**6*(42.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(27.0d0*dx2+2.0d0*dx1)))
      ci(6)=9.0d0*dx1**5*(28.0d0*dx2**3+5.0d0*dx1*(14.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))
      ci(7)=6.0d0*dx1**4*(35.0d0*dx2**3+dx1*(126.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+4.0d0*dx1)))
      ci(8)=6.0d0*dx1**3*(20.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1)))
      ci(9)=9.0d0*dx1**2*(5.0d0*dx2**3+2.0d0*dx1*(20.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)))
      ci(10)=5.0d0*dx1*(2.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)))
      ci(11)=dx2**3+15.0d0*dx1*(2.0d0*dx2**2+dx1*(9.0d0*dx2+8.0d0*dx1))
      ci(12)=3.0d0*(dx2**2+5.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(13)=3.0d0*dx2+10.0d0*dx1
      ci(14)=1.0d0
    case(4)
      ci(1)=dx1**10*dx2**4
      ci(2)=2.0d0*dx1**9*dx2**3*(5.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**8*dx2**2*(45.0d0*dx2**2+2.0d0*dx1*(20.0d0*dx2+3.0d0*dx1))
      ci(4)=4.0d0*dx1**7*dx2*(30.0d0*dx2**3+dx1*(45.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))
      ci(5)=dx1**6*(210.0d0*dx2**4+dx1*(480.0d0*dx2**3+dx1*(270.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))
      ci(6)=2.0d0*dx1**5*(126.0d0*dx2**4+5.0d0*dx1*(84.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))))
      ci(7)=3.0d0*dx1**4*(70.0d0*dx2**4+dx1*(336.0d0*dx2**3+5.0d0*dx1*(84.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1))))
      ci(8)=24.0d0*dx1**3*(5.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1*(63.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+dx1))))
      ci(9)=3.0d0*dx1**2*(15.0d0*dx2**4+2.0d0*dx1*(80.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1))))
      ci(10)=2.0d0*dx1*(5.0d0*dx2**4+6.0d0*dx1*(15.0d0*dx2**3+dx1*(60.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+3.0d0*dx1))))
      ci(11)=dx2**4+10.0d0*dx1*(4.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+dx1*(16.0d0*dx2+7.0d0*dx1)))
      ci(12)=4.0d0*(dx2**3+15.0d0*dx1*(dx2**2+dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(13)=6.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+9.0d0*dx1)
      ci(14)=2.0d0*(2.0d0*dx2+5.0d0*dx1)
      ci(15)=1.0d0
    case(5)
      ci(1)=dx1**10*dx2**5
      ci(2)=5.0d0*dx1**9*dx2**4*(2.0d0*dx2+dx1)
      ci(3)=5.0d0*dx1**8*dx2**3*(9.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1))
      ci(4)=5.0d0*dx1**7*dx2**2*(24.0d0*dx2**3+dx1*(45.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**6*dx2*(42.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(6)=dx1**5*(252.0d0*dx2**5+dx1*(1050.0d0*dx2**4+dx1*(1200.0d0*dx2**3+dx1*(450.0d0*dx2**2+dx1*(50.0d0*dx2+dx1)))))
      ci(7)=5.0d0*dx1**4*(42.0d0*dx2**5+dx1*(252.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))))
      ci(8)=15.0d0*dx1**3*(8.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1)))))
      ci(9)=15.0d0*dx1**2*(3.0d0*dx2**5+2.0d0*dx1*(20.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(35.0d0*dx2+4.0d0*dx1)) &
)))
      ci(10)=5.0d0*dx1*(2.0d0*dx2**5+3.0d0*dx1*(15.0d0*dx2**4+2.0d0*dx1*(40.0d0*dx2**3+7.0d0*dx1*(10.0d0*dx2**2+dx1*(6.0d0*dx2+dx1 &
)))))
      ci(11)=dx2**5+2.0d0*dx1*(25.0d0*dx2**4+3.0d0*dx1*(75.0d0*dx2**3+dx1*(200.0d0*dx2**2+7.0d0*dx1*(25.0d0*dx2+6.0d0*dx1))))
      ci(12)=5.0d0*(dx2**4+2.0d0*dx1*(10.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+dx1*(20.0d0*dx2+7.0d0*dx1))))
      ci(13)=5.0d0*(2.0d0*dx2**3+dx1*(20.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+8.0d0*dx1)))
      ci(14)=5.0d0*(2.0d0*dx2**2+dx1*(10.0d0*dx2+9.0d0*dx1))
      ci(15)=5.0d0*(dx2+2.0d0*dx1)
      ci(16)=1.0d0
    case(6)
      ci(1)=dx1**10*dx2**6
      ci(2)=2.0d0*dx1**9*dx2**5*(5.0d0*dx2+3.0d0*dx1)
      ci(3)=15.0d0*dx1**8*dx2**4*(3.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))
      ci(4)=10.0d0*dx1**7*dx2**3*(12.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(5)=5.0d0*dx1**6*dx2**2*(42.0d0*dx2**4+dx1*(144.0d0*dx2**3+dx1*(135.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))))
      ci(6)=6.0d0*dx1**5*dx2*(42.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(300.0d0*dx2**3+dx1*(150.0d0*dx2**2+dx1*(25.0d0*dx2+dx1)))))
      ci(7)=dx1**4*(210.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1*(3150.0d0*dx2**4+dx1*(2400.0d0*dx2**3+dx1*(675.0d0*dx2**2+dx1* &
(60.0d0*dx2+dx1))))))
      ci(8)=10.0d0*dx1**3*(12.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(378.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1* &
(27.0d0*dx2+dx1))))))
      ci(9)=45.0d0*dx1**2*(dx2**6+dx1*(16.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(112.0d0*dx2**3+dx1*(70.0d0*dx2**2+dx1*(16.0d0*dx2+dx1 &
))))))
      ci(10)=10.0d0*dx1*(dx2**6+3.0d0*dx1*(9.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1* &
(21.0d0*dx2+2.0d0*dx1))))))
      ci(11)=dx2**6+3.0d0*dx1*(20.0d0*dx2**5+dx1*(225.0d0*dx2**4+2.0d0*dx1*(400.0d0*dx2**3+7.0d0*dx1*(75.0d0*dx2**2+dx1* &
(36.0d0*dx2+5.0d0*dx1)))))
      ci(12)=6.0d0*(dx2**5+dx1*(25.0d0*dx2**4+6.0d0*dx1*(25.0d0*dx2**3+dx1*(50.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+dx1)))))
      ci(13)=5.0d0*(3.0d0*dx2**4+dx1*(40.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+2.0d0*dx1*(24.0d0*dx2+7.0d0*dx1))))
      ci(14)=10.0d0*(2.0d0*dx2**3+3.0d0*dx1*(5.0d0*dx2**2+dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(15)=15.0d0*(dx2**2+dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(16)=2.0d0*(3.0d0*dx2+5.0d0*dx1)
      ci(17)=1.0d0
    case(7)
      ci(1)=dx1**10*dx2**7
      ci(2)=dx1**9*dx2**6*(10.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**8*dx2**5*(45.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+3.0d0*dx1))
      ci(4)=5.0d0*dx1**7*dx2**4*(24.0d0*dx2**3+7.0d0*dx1*(9.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))
      ci(5)=35.0d0*dx1**6*dx2**3*(6.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**5*dx2**2*(36.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(360.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1*(50.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=7.0d0*dx1**4*dx2*(30.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1*(600.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1* &
(30.0d0*dx2+dx1))))))
      ci(8)=dx1**3*(120.0d0*dx2**7+dx1*(1470.0d0*dx2**6+dx1*(5292.0d0*dx2**5+dx1*(7350.0d0*dx2**4+dx1*(4200.0d0*dx2**3+dx1* &
(945.0d0*dx2**2+dx1*(70.0d0*dx2+dx1)))))))
      ci(9)=5.0d0*dx1**2*(9.0d0*dx2**7+dx1*(168.0d0*dx2**6+dx1*(882.0d0*dx2**5+dx1*(1764.0d0*dx2**4+dx1*(1470.0d0*dx2**3+dx1* &
(504.0d0*dx2**2+dx1*(63.0d0*dx2+2.0d0*dx1)))))))
      ci(10)=5.0d0*dx1*(2.0d0*dx2**7+3.0d0*dx1*(21.0d0*dx2**6+dx1*(168.0d0*dx2**5+dx1*(490.0d0*dx2**4+dx1*(588.0d0*dx2**3+dx1* &
(294.0d0*dx2**2+dx1*(56.0d0*dx2+3.0d0*dx1)))))))
      ci(11)=dx2**7+dx1*(70.0d0*dx2**6+3.0d0*dx1*(315.0d0*dx2**5+2.0d0*dx1*(700.0d0*dx2**4+dx1*(1225.0d0*dx2**3+dx1* &
(882.0d0*dx2**2+5.0d0*dx1*(49.0d0*dx2+4.0d0*dx1))))))
      ci(12)=7.0d0*(dx2**6+3.0d0*dx1*(10.0d0*dx2**5+dx1*(75.0d0*dx2**4+2.0d0*dx1*(100.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1* &
(42.0d0*dx2+5.0d0*dx1))))))
      ci(13)=7.0d0*(3.0d0*dx2**5+dx1*(50.0d0*dx2**4+3.0d0*dx1*(75.0d0*dx2**3+2.0d0*dx1*(60.0d0*dx2**2+dx1*(35.0d0*dx2+6.0d0*dx1))) &
))
      ci(14)=35.0d0*(dx2**4+dx1*(10.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(15)=5.0d0*(7.0d0*dx2**3+3.0d0*dx1*(14.0d0*dx2**2+dx1*(21.0d0*dx2+8.0d0*dx1)))
      ci(16)=21.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+9.0d0*dx1)
      ci(17)=7.0d0*dx2+10.0d0*dx1
      ci(18)=1.0d0
    case(8)
      ci(1)=dx1**10*dx2**8
      ci(2)=2.0d0*dx1**9*dx2**7*(5.0d0*dx2+4.0d0*dx1)
      ci(3)=dx1**8*dx2**6*(45.0d0*dx2**2+4.0d0*dx1*(20.0d0*dx2+7.0d0*dx1))
      ci(4)=8.0d0*dx1**7*dx2**5*(15.0d0*dx2**3+dx1*(45.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+dx1)))
      ci(5)=10.0d0*dx1**6*dx2**4*(21.0d0*dx2**4+dx1*(96.0d0*dx2**3+7.0d0*dx1*(18.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(6)=28.0d0*dx1**5*dx2**3*(9.0d0*dx2**5+dx1*(60.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(25.0d0*dx2+2.0d0*dx1 &
)))))
      ci(7)=14.0d0*dx1**4*dx2**2*(15.0d0*dx2**6+dx1*(144.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1*(480.0d0*dx2**3+dx1*(225.0d0*dx2**2 &
+2.0d0*dx1*(20.0d0*dx2+dx1))))))
      ci(8)=8.0d0*dx1**3*dx2*(15.0d0*dx2**7+dx1*(210.0d0*dx2**6+dx1*(882.0d0*dx2**5+dx1*(1470.0d0*dx2**4+dx1*(1050.0d0*dx2**3+dx1* &
(315.0d0*dx2**2+dx1*(35.0d0*dx2+dx1)))))))
      ci(9)=dx1**2*(45.0d0*dx2**8+dx1*(960.0d0*dx2**7+dx1*(5880.0d0*dx2**6+dx1*(14112.0d0*dx2**5+dx1*(14700.0d0*dx2**4+dx1* &
(6720.0d0*dx2**3+dx1*(1260.0d0*dx2**2+dx1*(80.0d0*dx2+dx1))))))))
      ci(10)=10.0d0*dx1*(dx2**8+dx1*(36.0d0*dx2**7+dx1*(336.0d0*dx2**6+dx1*(1176.0d0*dx2**5+dx1*(1764.0d0*dx2**4+dx1* &
(1176.0d0*dx2**3+dx1*(336.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))))))
      ci(11)=dx2**8+dx1*(80.0d0*dx2**7+3.0d0*dx1*(420.0d0*dx2**6+dx1*(2240.0d0*dx2**5+dx1*(4900.0d0*dx2**4+dx1*(4704.0d0*dx2**3 &
+5.0d0*dx1*(392.0d0*dx2**2+dx1*(64.0d0*dx2+3.0d0*dx1)))))))
      ci(12)=8.0d0*(dx2**7+dx1*(35.0d0*dx2**6+3.0d0*dx1*(105.0d0*dx2**5+dx1*(350.0d0*dx2**4+dx1*(490.0d0*dx2**3+dx1* &
(294.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+dx1)))))))
      ci(13)=14.0d0*(2.0d0*dx2**6+dx1*(40.0d0*dx2**5+3.0d0*dx1*(75.0d0*dx2**4+dx1*(160.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1* &
(48.0d0*dx2+5.0d0*dx1))))))
      ci(14)=28.0d0*(2.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(30.0d0*dx2**3+dx1*(40.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1)))))
      ci(15)=10.0d0*(7.0d0*dx2**4+dx1*(56.0d0*dx2**3+3.0d0*dx1*(42.0d0*dx2**2+dx1*(32.0d0*dx2+7.0d0*dx1))))
      ci(16)=8.0d0*(7.0d0*dx2**3+5.0d0*dx1*(7.0d0*dx2**2+3.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(17)=28.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+9.0d0*dx1)
      ci(18)=2.0d0*(4.0d0*dx2+5.0d0*dx1)
      ci(19)=1.0d0
    case(9)
      ci(1)=dx1**10*dx2**9
      ci(2)=dx1**9*dx2**8*(10.0d0*dx2+9.0d0*dx1)
      ci(3)=9.0d0*dx1**8*dx2**7*(5.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))
      ci(4)=3.0d0*dx1**7*dx2**6*(40.0d0*dx2**3+dx1*(135.0d0*dx2**2+4.0d0*dx1*(30.0d0*dx2+7.0d0*dx1)))
      ci(5)=6.0d0*dx1**6*dx2**5*(35.0d0*dx2**4+dx1*(180.0d0*dx2**3+dx1*(270.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+3.0d0*dx1))))
      ci(6)=18.0d0*dx1**5*dx2**4*(14.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(240.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+dx1*(10.0d0*dx2 &
+dx1)))))
      ci(7)=42.0d0*dx1**4*dx2**3*(5.0d0*dx2**6+dx1*(54.0d0*dx2**5+dx1*(180.0d0*dx2**4+dx1*(240.0d0*dx2**3+dx1*(135.0d0*dx2**2 &
+2.0d0*dx1*(15.0d0*dx2+dx1))))))
      ci(8)=6.0d0*dx1**3*dx2**2*(20.0d0*dx2**7+dx1*(315.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1*(2520.0d0*dx2**3 &
+dx1*(945.0d0*dx2**2+2.0d0*dx1*(70.0d0*dx2+3.0d0*dx1)))))))
      ci(9)=9.0d0*dx1**2*dx2*(5.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(840.0d0*dx2**6+dx1*(2352.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(1680.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))))))
      ci(10)=dx1*(10.0d0*dx2**9+dx1*(405.0d0*dx2**8+dx1*(4320.0d0*dx2**7+dx1*(17640.0d0*dx2**6+dx1*(31752.0d0*dx2**5+dx1* &
(26460.0d0*dx2**4+dx1*(10080.0d0*dx2**3+dx1*(1620.0d0*dx2**2+dx1*(90.0d0*dx2+dx1)))))))))
      ci(11)=dx2**9+dx1*(90.0d0*dx2**8+dx1*(1620.0d0*dx2**7+dx1*(10080.0d0*dx2**6+dx1*(26460.0d0*dx2**5+dx1*(31752.0d0*dx2**4 &
+5.0d0*dx1*(3528.0d0*dx2**3+dx1*(864.0d0*dx2**2+dx1*(81.0d0*dx2+2.0d0*dx1))))))))
      ci(12)=9.0d0*(dx2**8+dx1*(40.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(1680.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1*(2352.0d0*dx2**3 &
+5.0d0*dx1*(168.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))))))
      ci(13)=6.0d0*(6.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(945.0d0*dx2**5+dx1*(2520.0d0*dx2**4+dx1*(2940.0d0*dx2**3+dx1* &
(1512.0d0*dx2**2+5.0d0*dx1*(63.0d0*dx2+4.0d0*dx1)))))))
      ci(14)=42.0d0*(2.0d0*dx2**6+dx1*(30.0d0*dx2**5+dx1*(135.0d0*dx2**4+dx1*(240.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1*(54.0d0*dx2 &
+5.0d0*dx1))))))
      ci(15)=18.0d0*(7.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(210.0d0*dx2**3+dx1*(240.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+2.0d0*dx1)))))
      ci(16)=6.0d0*(21.0d0*dx2**4+5.0d0*dx1*(28.0d0*dx2**3+dx1*(54.0d0*dx2**2+dx1*(36.0d0*dx2+7.0d0*dx1))))
      ci(17)=3.0d0*(28.0d0*dx2**3+5.0d0*dx1*(24.0d0*dx2**2+dx1*(27.0d0*dx2+8.0d0*dx1)))
      ci(18)=9.0d0*(4.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))
      ci(19)=9.0d0*dx2+10.0d0*dx1
      ci(20)=1.0d0
    case(10)
      ci(1)=dx1**10*dx2**10
      ci(2)=10.0d0*dx1**9*dx2**9*(dx2+dx1)
      ci(3)=5.0d0*dx1**8*dx2**8*(9.0d0*dx2**2+dx1*(20.0d0*dx2+9.0d0*dx1))
      ci(4)=30.0d0*dx1**7*dx2**7*(4.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(15.0d0*dx2+4.0d0*dx1)))
      ci(5)=15.0d0*dx1**6*dx2**6*(14.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(135.0d0*dx2**2+2.0d0*dx1*(40.0d0*dx2+7.0d0*dx1))))
      ci(6)=12.0d0*dx1**5*dx2**5*(21.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(450.0d0*dx2**3+dx1*(450.0d0*dx2**2+7.0d0*dx1*(25.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=30.0d0*dx1**4*dx2**4*(7.0d0*dx2**6+dx1*(84.0d0*dx2**5+dx1*(315.0d0*dx2**4+dx1*(480.0d0*dx2**3+7.0d0*dx1*(45.0d0*dx2**2 &
+dx1*(12.0d0*dx2+dx1))))))
      ci(8)=60.0d0*dx1**3*dx2**3*(2.0d0*dx2**7+dx1*(35.0d0*dx2**6+dx1*(189.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1* &
(189.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1)))))))
      ci(9)=15.0d0*dx1**2*dx2**2*(3.0d0*dx2**8+dx1*(80.0d0*dx2**7+dx1*(630.0d0*dx2**6+dx1*(2016.0d0*dx2**5+dx1*(2940.0d0*dx2**4 &
+dx1*(2016.0d0*dx2**3+dx1*(630.0d0*dx2**2+dx1*(80.0d0*dx2+3.0d0*dx1))))))))
      ci(10)=10.0d0*dx1*dx2*(dx2**9+dx1*(45.0d0*dx2**8+dx1*(540.0d0*dx2**7+dx1*(2520.0d0*dx2**6+dx1*(5292.0d0*dx2**5+dx1* &
(5292.0d0*dx2**4+dx1*(2520.0d0*dx2**3+dx1*(540.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))))))))
      ci(11)=dx2**10+dx1*(100.0d0*dx2**9+dx1*(2025.0d0*dx2**8+dx1*(14400.0d0*dx2**7+dx1*(44100.0d0*dx2**6+dx1*(63504.0d0*dx2**5 &
+dx1*(44100.0d0*dx2**4+dx1*(14400.0d0*dx2**3+dx1*(2025.0d0*dx2**2+dx1*(100.0d0*dx2+dx1)))))))))
      ci(12)=10.0d0*(dx2**9+dx1*(45.0d0*dx2**8+dx1*(540.0d0*dx2**7+dx1*(2520.0d0*dx2**6+dx1*(5292.0d0*dx2**5+dx1*(5292.0d0*dx2**4 &
+dx1*(2520.0d0*dx2**3+dx1*(540.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))))))))
      ci(13)=15.0d0*(3.0d0*dx2**8+dx1*(80.0d0*dx2**7+dx1*(630.0d0*dx2**6+dx1*(2016.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(2016.0d0*dx2**3+dx1*(630.0d0*dx2**2+dx1*(80.0d0*dx2+3.0d0*dx1))))))))
      ci(14)=60.0d0*(2.0d0*dx2**7+dx1*(35.0d0*dx2**6+dx1*(189.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1* &
(189.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1)))))))
      ci(15)=30.0d0*(7.0d0*dx2**6+dx1*(84.0d0*dx2**5+dx1*(315.0d0*dx2**4+dx1*(480.0d0*dx2**3+7.0d0*dx1*(45.0d0*dx2**2+dx1* &
(12.0d0*dx2+dx1))))))
      ci(16)=12.0d0*(21.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(450.0d0*dx2**3+dx1*(450.0d0*dx2**2+7.0d0*dx1*(25.0d0*dx2+3.0d0*dx1)))) &
)
      ci(17)=15.0d0*(14.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(135.0d0*dx2**2+2.0d0*dx1*(40.0d0*dx2+7.0d0*dx1))))
      ci(18)=30.0d0*(4.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(15.0d0*dx2+4.0d0*dx1)))
      ci(19)=5.0d0*(9.0d0*dx2**2+dx1*(20.0d0*dx2+9.0d0*dx1))
      ci(20)=10.0d0*(dx2+dx1)
      ci(21)=1.0d0
    case(11)
      ci(1)=dx1**10*dx2**11
      ci(2)=dx1**9*dx2**10*(10.0d0*dx2+11.0d0*dx1)
      ci(3)=5.0d0*dx1**8*dx2**9*(9.0d0*dx2**2+11.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=5.0d0*dx1**7*dx2**8*(24.0d0*dx2**3+11.0d0*dx1*(9.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1)))
      ci(5)=15.0d0*dx1**6*dx2**7*(14.0d0*dx2**4+11.0d0*dx1*(8.0d0*dx2**3+dx1*(15.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1))))
      ci(6)=3.0d0*dx1**5*dx2**6*(84.0d0*dx2**5+11.0d0*dx1*(70.0d0*dx2**4+dx1*(200.0d0*dx2**3+dx1*(225.0d0*dx2**2+2.0d0*dx1* &
(50.0d0*dx2+7.0d0*dx1)))))
      ci(7)=6.0d0*dx1**4*dx2**5*(35.0d0*dx2**6+11.0d0*dx1*(42.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(300.0d0*dx2**3+dx1* &
(225.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+dx1))))))
      ci(8)=30.0d0*dx1**3*dx2**4*(4.0d0*dx2**7+11.0d0*dx1*(7.0d0*dx2**6+dx1*(42.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(120.0d0*dx2**3 &
+dx1*(63.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)))))))
      ci(9)=15.0d0*dx1**2*dx2**3*(3.0d0*dx2**8+11.0d0*dx1*(8.0d0*dx2**7+dx1*(70.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1*(420.0d0*dx2**4 &
+dx1*(336.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))))))
      ci(10)=5.0d0*dx1*dx2**2*(2.0d0*dx2**9+11.0d0*dx1*(9.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(630.0d0*dx2**6+dx1*(1512.0d0*dx2**5 &
+dx1*(1764.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(270.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))))))
      ci(11)=dx2*(dx2**10+11.0d0*dx1*(10.0d0*dx2**9+dx1*(225.0d0*dx2**8+dx1*(1800.0d0*dx2**7+dx1*(6300.0d0*dx2**6+dx1* &
(10584.0d0*dx2**5+dx1*(8820.0d0*dx2**4+dx1*(3600.0d0*dx2**3+dx1*(675.0d0*dx2**2+dx1*(50.0d0*dx2+dx1))))))))))
      ci(12)=11.0d0*dx2**10+dx1*(550.0d0*dx2**9+dx1*(7425.0d0*dx2**8+dx1*(39600.0d0*dx2**7+dx1*(97020.0d0*dx2**6+dx1* &
(116424.0d0*dx2**5+dx1*(69300.0d0*dx2**4+dx1*(19800.0d0*dx2**3+dx1*(2475.0d0*dx2**2+dx1*(110.0d0*dx2+dx1)))))))))
      ci(13)=5.0d0*(11.0d0*dx2**9+dx1*(330.0d0*dx2**8+dx1*(2970.0d0*dx2**7+dx1*(11088.0d0*dx2**6+dx1*(19404.0d0*dx2**5+dx1* &
(16632.0d0*dx2**4+dx1*(6930.0d0*dx2**3+dx1*(1320.0d0*dx2**2+dx1*(99.0d0*dx2+2.0d0*dx1)))))))))
      ci(14)=15.0d0*(11.0d0*dx2**8+dx1*(220.0d0*dx2**7+dx1*(1386.0d0*dx2**6+dx1*(3696.0d0*dx2**5+dx1*(4620.0d0*dx2**4+dx1* &
(2772.0d0*dx2**3+dx1*(770.0d0*dx2**2+dx1*(88.0d0*dx2+3.0d0*dx1))))))))
      ci(15)=30.0d0*(11.0d0*dx2**7+dx1*(154.0d0*dx2**6+dx1*(693.0d0*dx2**5+dx1*(1320.0d0*dx2**4+dx1*(1155.0d0*dx2**3+dx1* &
(462.0d0*dx2**2+dx1*(77.0d0*dx2+4.0d0*dx1)))))))
      ci(16)=6.0d0*(77.0d0*dx2**6+dx1*(770.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(3300.0d0*dx2**3+7.0d0*dx1*(275.0d0*dx2**2+dx1* &
(66.0d0*dx2+5.0d0*dx1))))))
      ci(17)=3.0d0*(154.0d0*dx2**5+dx1*(1100.0d0*dx2**4+dx1*(2475.0d0*dx2**3+2.0d0*dx1*(1100.0d0*dx2**2+7.0d0*dx1*(55.0d0*dx2 &
+6.0d0*dx1)))))
      ci(18)=15.0d0*(22.0d0*dx2**4+dx1*(110.0d0*dx2**3+dx1*(165.0d0*dx2**2+2.0d0*dx1*(44.0d0*dx2+7.0d0*dx1))))
      ci(19)=5.0d0*(33.0d0*dx2**3+dx1*(110.0d0*dx2**2+3.0d0*dx1*(33.0d0*dx2+8.0d0*dx1)))
      ci(20)=5.0d0*(11.0d0*dx2**2+dx1*(22.0d0*dx2+9.0d0*dx1))
      ci(21)=11.0d0*dx2+10.0d0*dx1
      ci(22)=1.0d0
    case(12)
      ci(1)=dx1**10*dx2**12
      ci(2)=2.0d0*dx1**9*dx2**11*(5.0d0*dx2+6.0d0*dx1)
      ci(3)=3.0d0*dx1**8*dx2**10*(15.0d0*dx2**2+2.0d0*dx1*(20.0d0*dx2+11.0d0*dx1))
      ci(4)=20.0d0*dx1**7*dx2**9*(6.0d0*dx2**3+dx1*(27.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**6*dx2**8*(42.0d0*dx2**4+dx1*(288.0d0*dx2**3+11.0d0*dx1*(54.0d0*dx2**2+dx1*(40.0d0*dx2+9.0d0*dx1))))
      ci(6)=18.0d0*dx1**5*dx2**7*(14.0d0*dx2**5+dx1*(140.0d0*dx2**4+11.0d0*dx1*(40.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(25.0d0*dx2 &
+4.0d0*dx1)))))
      ci(7)=3.0d0*dx1**4*dx2**6*(70.0d0*dx2**6+dx1*(1008.0d0*dx2**5+11.0d0*dx1*(420.0d0*dx2**4+dx1*(800.0d0*dx2**3+dx1* &
(675.0d0*dx2**2+4.0d0*dx1*(60.0d0*dx2+7.0d0*dx1))))))
      ci(8)=24.0d0*dx1**3*dx2**5*(5.0d0*dx2**7+dx1*(105.0d0*dx2**6+11.0d0*dx1*(63.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1* &
(225.0d0*dx2**3+dx1*(135.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1)))))))
      ci(9)=45.0d0*dx1**2*dx2**4*(dx2**8+dx1*(32.0d0*dx2**7+11.0d0*dx1*(28.0d0*dx2**6+dx1*(112.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1* &
(192.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))))))
      ci(10)=10.0d0*dx1*dx2**3*(dx2**9+dx1*(54.0d0*dx2**8+11.0d0*dx1*(72.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(1134.0d0*dx2**5+dx1* &
(1512.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(324.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))))))))
      ci(11)=dx2**2*(dx2**10+dx1*(120.0d0*dx2**9+11.0d0*dx1*(270.0d0*dx2**8+dx1*(2400.0d0*dx2**7+dx1*(9450.0d0*dx2**6+dx1* &
(18144.0d0*dx2**5+dx1*(17640.0d0*dx2**4+dx1*(8640.0d0*dx2**3+dx1*(2025.0d0*dx2**2+2.0d0*dx1*(100.0d0*dx2+3.0d0*dx1))))))))))
      ci(12)=12.0d0*dx2*(dx2**10+dx1*(55.0d0*dx2**9+dx1*(825.0d0*dx2**8+dx1*(4950.0d0*dx2**7+dx1*(13860.0d0*dx2**6+dx1* &
(19404.0d0*dx2**5+dx1*(13860.0d0*dx2**4+dx1*(4950.0d0*dx2**3+dx1*(825.0d0*dx2**2+dx1*(55.0d0*dx2+dx1))))))))))
      ci(13)=66.0d0*dx2**10+dx1*(2200.0d0*dx2**9+dx1*(22275.0d0*dx2**8+dx1*(95040.0d0*dx2**7+dx1*(194040.0d0*dx2**6+dx1* &
(199584.0d0*dx2**5+dx1*(103950.0d0*dx2**4+dx1*(26400.0d0*dx2**3+dx1*(2970.0d0*dx2**2+dx1*(120.0d0*dx2+dx1)))))))))
      ci(14)=10.0d0*(22.0d0*dx2**9+dx1*(495.0d0*dx2**8+dx1*(3564.0d0*dx2**7+dx1*(11088.0d0*dx2**6+dx1*(16632.0d0*dx2**5+dx1* &
(12474.0d0*dx2**4+dx1*(4620.0d0*dx2**3+dx1*(792.0d0*dx2**2+dx1*(54.0d0*dx2+dx1)))))))))
      ci(15)=45.0d0*(11.0d0*dx2**8+dx1*(176.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2112.0d0*dx2**5+dx1*(2310.0d0*dx2**4+dx1* &
(1232.0d0*dx2**3+dx1*(308.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))))))
      ci(16)=24.0d0*(33.0d0*dx2**7+dx1*(385.0d0*dx2**6+dx1*(1485.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(1925.0d0*dx2**3+dx1* &
(693.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+dx1)))))))
      ci(17)=3.0d0*(308.0d0*dx2**6+dx1*(2640.0d0*dx2**5+dx1*(7425.0d0*dx2**4+2.0d0*dx1*(4400.0d0*dx2**3+7.0d0*dx1*(330.0d0*dx2**2 &
+dx1*(72.0d0*dx2+5.0d0*dx1))))))
      ci(18)=18.0d0*(44.0d0*dx2**5+dx1*(275.0d0*dx2**4+2.0d0*dx1*(275.0d0*dx2**3+dx1*(220.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+dx1)))) &
)
      ci(19)=5.0d0*(99.0d0*dx2**4+2.0d0*dx1*(220.0d0*dx2**3+3.0d0*dx1*(99.0d0*dx2**2+dx1*(48.0d0*dx2+7.0d0*dx1))))
      ci(20)=20.0d0*(11.0d0*dx2**3+3.0d0*dx1*(11.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(21)=3.0d0*(22.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(22)=2.0d0*(6.0d0*dx2+5.0d0*dx1)
      ci(23)=1.0d0
    case(13)
      ci(1)=dx1**10*dx2**13
      ci(2)=dx1**9*dx2**12*(10.0d0*dx2+13.0d0*dx1)
      ci(3)=dx1**8*dx2**11*(45.0d0*dx2**2+26.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))
      ci(4)=dx1**7*dx2**10*(120.0d0*dx2**3+13.0d0*dx1*(45.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2+11.0d0*dx1)))
      ci(5)=5.0d0*dx1**6*dx2**9*(42.0d0*dx2**4+13.0d0*dx1*(24.0d0*dx2**3+dx1*(54.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(6)=dx1**5*dx2**8*(252.0d0*dx2**5+13.0d0*dx1*(210.0d0*dx2**4+dx1*(720.0d0*dx2**3+11.0d0*dx1*(90.0d0*dx2**2+dx1*(50.0d0*dx2 &
+9.0d0*dx1)))))
      ci(7)=3.0d0*dx1**4*dx2**7*(70.0d0*dx2**6+13.0d0*dx1*(84.0d0*dx2**5+dx1*(420.0d0*dx2**4+11.0d0*dx1*(80.0d0*dx2**3+dx1* &
(75.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+2.0d0*dx1))))))
      ci(8)=3.0d0*dx1**3*dx2**6*(40.0d0*dx2**7+13.0d0*dx1*(70.0d0*dx2**6+dx1*(504.0d0*dx2**5+11.0d0*dx1*(140.0d0*dx2**4+dx1* &
(200.0d0*dx2**3+dx1*(135.0d0*dx2**2+4.0d0*dx1*(10.0d0*dx2+dx1)))))))
      ci(9)=3.0d0*dx1**2*dx2**5*(15.0d0*dx2**8+13.0d0*dx1*(40.0d0*dx2**7+dx1*(420.0d0*dx2**6+11.0d0*dx1*(168.0d0*dx2**5+dx1* &
(350.0d0*dx2**4+dx1*(360.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))))))))
      ci(10)=5.0d0*dx1*dx2**4*(2.0d0*dx2**9+13.0d0*dx1*(9.0d0*dx2**8+dx1*(144.0d0*dx2**7+11.0d0*dx1*(84.0d0*dx2**6+dx1* &
(252.0d0*dx2**5+dx1*(378.0d0*dx2**4+dx1*(288.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))))))))
      ci(11)=dx2**3*(dx2**10+13.0d0*dx1*(10.0d0*dx2**9+dx1*(270.0d0*dx2**8+11.0d0*dx1*(240.0d0*dx2**7+dx1*(1050.0d0*dx2**6+dx1* &
(2268.0d0*dx2**5+dx1*(2520.0d0*dx2**4+dx1*(1440.0d0*dx2**3+dx1*(405.0d0*dx2**2+2.0d0*dx1*(25.0d0*dx2+dx1))))))))))
      ci(12)=13.0d0*dx2**2*(dx2**10+dx1*(60.0d0*dx2**9+dx1*(990.0d0*dx2**8+dx1*(6600.0d0*dx2**7+dx1*(20790.0d0*dx2**6+dx1* &
(33264.0d0*dx2**5+dx1*(27720.0d0*dx2**4+dx1*(11880.0d0*dx2**3+dx1*(2475.0d0*dx2**2+2.0d0*dx1*(110.0d0*dx2+3.0d0*dx1))))))))))
      ci(13)=13.0d0*dx2*(6.0d0*dx2**10+dx1*(220.0d0*dx2**9+dx1*(2475.0d0*dx2**8+dx1*(11880.0d0*dx2**7+dx1*(27720.0d0*dx2**6+dx1* &
(33264.0d0*dx2**5+dx1*(20790.0d0*dx2**4+dx1*(6600.0d0*dx2**3+dx1*(990.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))))))))
      ci(14)=286.0d0*dx2**10+dx1*(7150.0d0*dx2**9+dx1*(57915.0d0*dx2**8+dx1*(205920.0d0*dx2**7+dx1*(360360.0d0*dx2**6+dx1* &
(324324.0d0*dx2**5+dx1*(150150.0d0*dx2**4+dx1*(34320.0d0*dx2**3+dx1*(3510.0d0*dx2**2+dx1*(130.0d0*dx2+dx1)))))))))
      ci(15)=5.0d0*(143.0d0*dx2**9+dx1*(2574.0d0*dx2**8+dx1*(15444.0d0*dx2**7+dx1*(41184.0d0*dx2**6+dx1*(54054.0d0*dx2**5+dx1* &
(36036.0d0*dx2**4+dx1*(12012.0d0*dx2**3+dx1*(1872.0d0*dx2**2+dx1*(117.0d0*dx2+2.0d0*dx1)))))))))
      ci(16)=3.0d0*(429.0d0*dx2**8+dx1*(5720.0d0*dx2**7+dx1*(25740.0d0*dx2**6+dx1*(51480.0d0*dx2**5+dx1*(50050.0d0*dx2**4+dx1* &
(24024.0d0*dx2**3+5.0d0*dx1*(1092.0d0*dx2**2+dx1*(104.0d0*dx2+3.0d0*dx1))))))))
      ci(17)=3.0d0*(572.0d0*dx2**7+dx1*(5720.0d0*dx2**6+dx1*(19305.0d0*dx2**5+2.0d0*dx1*(14300.0d0*dx2**4+dx1*(10010.0d0*dx2**3 &
+dx1*(3276.0d0*dx2**2+5.0d0*dx1*(91.0d0*dx2+4.0d0*dx1)))))))
      ci(18)=3.0d0*(572.0d0*dx2**6+dx1*(4290.0d0*dx2**5+dx1*(10725.0d0*dx2**4+2.0d0*dx1*(5720.0d0*dx2**3+7.0d0*dx1*(390.0d0*dx2**2 &
+dx1*(78.0d0*dx2+5.0d0*dx1))))))
      ci(19)=1287.0d0*dx2**5+2.0d0*dx1*(3575.0d0*dx2**4+3.0d0*dx1*(2145.0d0*dx2**3+dx1*(1560.0d0*dx2**2+7.0d0*dx1*(65.0d0*dx2 &
+6.0d0*dx1))))
      ci(20)=5.0d0*(143.0d0*dx2**4+2.0d0*dx1*(286.0d0*dx2**3+3.0d0*dx1*(117.0d0*dx2**2+dx1*(52.0d0*dx2+7.0d0*dx1))))
      ci(21)=286.0d0*dx2**3+15.0d0*dx1*(52.0d0*dx2**2+dx1*(39.0d0*dx2+8.0d0*dx1))
      ci(22)=78.0d0*dx2**2+5.0d0*dx1*(26.0d0*dx2+9.0d0*dx1)
      ci(23)=13.0d0*dx2+10.0d0*dx1
      ci(24)=1.0d0
    case(14)
      ci(1)=dx1**10*dx2**14
      ci(2)=2.0d0*dx1**9*dx2**13*(5.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**8*dx2**12*(45.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+13.0d0*dx1))
      ci(4)=2.0d0*dx1**7*dx2**11*(60.0d0*dx2**3+7.0d0*dx1*(45.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)))
      ci(5)=7.0d0*dx1**6*dx2**10*(30.0d0*dx2**4+dx1*(240.0d0*dx2**3+13.0d0*dx1*(45.0d0*dx2**2+dx1*(40.0d0*dx2+11.0d0*dx1))))
      ci(6)=14.0d0*dx1**5*dx2**9*(18.0d0*dx2**5+dx1*(210.0d0*dx2**4+13.0d0*dx1*(60.0d0*dx2**3+dx1*(90.0d0*dx2**2+11.0d0*dx1* &
(5.0d0*dx2+dx1)))))
      ci(7)=7.0d0*dx1**4*dx2**8*(30.0d0*dx2**6+dx1*(504.0d0*dx2**5+13.0d0*dx1*(210.0d0*dx2**4+dx1*(480.0d0*dx2**3+11.0d0*dx1* &
(45.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1))))))
      ci(8)=6.0d0*dx1**3*dx2**7*(20.0d0*dx2**7+dx1*(490.0d0*dx2**6+13.0d0*dx1*(294.0d0*dx2**5+dx1*(980.0d0*dx2**4+11.0d0*dx1* &
(140.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(35.0d0*dx2+4.0d0*dx1)))))))
      ci(9)=3.0d0*dx1**2*dx2**6*(15.0d0*dx2**8+dx1*(560.0d0*dx2**7+13.0d0*dx1*(490.0d0*dx2**6+dx1*(2352.0d0*dx2**5+11.0d0*dx1* &
(490.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(315.0d0*dx2**2+dx1*(80.0d0*dx2+7.0d0*dx1))))))))
      ci(10)=2.0d0*dx1*dx2**5*(5.0d0*dx2**9+dx1*(315.0d0*dx2**8+13.0d0*dx1*(420.0d0*dx2**7+dx1*(2940.0d0*dx2**6+11.0d0*dx1* &
(882.0d0*dx2**5+dx1*(1470.0d0*dx2**4+dx1*(1260.0d0*dx2**3+dx1*(540.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+dx1)))))))))
      ci(11)=dx2**4*(dx2**10+dx1*(140.0d0*dx2**9+13.0d0*dx1*(315.0d0*dx2**8+dx1*(3360.0d0*dx2**7+11.0d0*dx1*(1470.0d0*dx2**6+dx1* &
(3528.0d0*dx2**5+dx1*(4410.0d0*dx2**4+dx1*(2880.0d0*dx2**3+7.0d0*dx1*(135.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))))))))
      ci(12)=14.0d0*dx2**3*(dx2**10+13.0d0*dx1*(5.0d0*dx2**9+dx1*(90.0d0*dx2**8+dx1*(660.0d0*dx2**7+dx1*(2310.0d0*dx2**6+dx1* &
(4158.0d0*dx2**5+dx1*(3960.0d0*dx2**4+dx1*(1980.0d0*dx2**3+dx1*(495.0d0*dx2**2+dx1*(55.0d0*dx2+2.0d0*dx1))))))))))
      ci(13)=91.0d0*dx2**2*(dx2**10+dx1*(40.0d0*dx2**9+dx1*(495.0d0*dx2**8+dx1*(2640.0d0*dx2**7+dx1*(6930.0d0*dx2**6+dx1* &
(9504.0d0*dx2**5+dx1*(6930.0d0*dx2**4+dx1*(2640.0d0*dx2**3+dx1*(495.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))))))))
      ci(14)=14.0d0*dx2*(26.0d0*dx2**10+dx1*(715.0d0*dx2**9+dx1*(6435.0d0*dx2**8+dx1*(25740.0d0*dx2**7+dx1*(51480.0d0*dx2**6+dx1* &
(54054.0d0*dx2**5+dx1*(30030.0d0*dx2**4+dx1*(8580.0d0*dx2**3+dx1*(1170.0d0*dx2**2+dx1*(65.0d0*dx2+dx1))))))))))
      ci(15)=1001.0d0*dx2**10+dx1*(20020.0d0*dx2**9+dx1*(135135.0d0*dx2**8+dx1*(411840.0d0*dx2**7+dx1*(630630.0d0*dx2**6+dx1* &
(504504.0d0*dx2**5+dx1*(210210.0d0*dx2**4+dx1*(43680.0d0*dx2**3+dx1*(4095.0d0*dx2**2+dx1*(140.0d0*dx2+dx1)))))))))
      ci(16)=2.0d0*(1001.0d0*dx2**9+dx1*(15015.0d0*dx2**8+dx1*(77220.0d0*dx2**7+dx1*(180180.0d0*dx2**6+dx1*(210210.0d0*dx2**5+dx1* &
(126126.0d0*dx2**4+5.0d0*dx1*(7644.0d0*dx2**3+dx1*(1092.0d0*dx2**2+dx1*(63.0d0*dx2+dx1)))))))))
      ci(17)=3.0d0*(1001.0d0*dx2**8+dx1*(11440.0d0*dx2**7+dx1*(45045.0d0*dx2**6+dx1*(80080.0d0*dx2**5+dx1*(70070.0d0*dx2**4+dx1* &
(30576.0d0*dx2**3+5.0d0*dx1*(1274.0d0*dx2**2+dx1*(112.0d0*dx2+3.0d0*dx1))))))))
      ci(18)=6.0d0*(572.0d0*dx2**7+dx1*(5005.0d0*dx2**6+dx1*(15015.0d0*dx2**5+2.0d0*dx1*(10010.0d0*dx2**4+dx1*(6370.0d0*dx2**3 &
+dx1*(1911.0d0*dx2**2+5.0d0*dx1*(49.0d0*dx2+2.0d0*dx1)))))))
      ci(19)=7.0d0*(429.0d0*dx2**6+dx1*(2860.0d0*dx2**5+3.0d0*dx1*(2145.0d0*dx2**4+2.0d0*dx1*(1040.0d0*dx2**3+dx1*(455.0d0*dx2**2 &
+dx1*(84.0d0*dx2+5.0d0*dx1))))))
      ci(20)=14.0d0*(143.0d0*dx2**5+dx1*(715.0d0*dx2**4+6.0d0*dx1*(195.0d0*dx2**3+dx1*(130.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1))) &
))
      ci(21)=7.0d0*(143.0d0*dx2**4+5.0d0*dx1*(104.0d0*dx2**3+3.0d0*dx1*(39.0d0*dx2**2+2.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(22)=2.0d0*(182.0d0*dx2**3+5.0d0*dx1*(91.0d0*dx2**2+3.0d0*dx1*(21.0d0*dx2+4.0d0*dx1)))
      ci(23)=91.0d0*dx2**2+5.0d0*dx1*(28.0d0*dx2+9.0d0*dx1)
      ci(24)=2.0d0*(7.0d0*dx2+5.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>14, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^11 * (x-x2)^n2 as sum_k=0^(11+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_11(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**11
      ci(2)=11.0d0*dx1**10
      ci(3)=55.0d0*dx1**9
      ci(4)=165.0d0*dx1**8
      ci(5)=330.0d0*dx1**7
      ci(6)=462.0d0*dx1**6
      ci(7)=462.0d0*dx1**5
      ci(8)=330.0d0*dx1**4
      ci(9)=165.0d0*dx1**3
      ci(10)=55.0d0*dx1**2
      ci(11)=11.0d0*dx1
      ci(12)=1.0d0
    case(1)
      ci(1)=dx1**11*dx2
      ci(2)=dx1**10*(11.0d0*dx2+dx1)
      ci(3)=11.0d0*dx1**9*(5.0d0*dx2+dx1)
      ci(4)=55.0d0*dx1**8*(3.0d0*dx2+dx1)
      ci(5)=165.0d0*dx1**7*(2.0d0*dx2+dx1)
      ci(6)=66.0d0*dx1**6*(7.0d0*dx2+5.0d0*dx1)
      ci(7)=462.0d0*dx1**5*(dx2+dx1)
      ci(8)=66.0d0*dx1**4*(5.0d0*dx2+7.0d0*dx1)
      ci(9)=165.0d0*dx1**3*(dx2+2.0d0*dx1)
      ci(10)=55.0d0*dx1**2*(dx2+3.0d0*dx1)
      ci(11)=11.0d0*dx1*(dx2+5.0d0*dx1)
      ci(12)=dx2+11.0d0*dx1
      ci(13)=1.0d0
    case(2)
      ci(1)=dx1**11*dx2**2
      ci(2)=dx1**10*dx2*(11.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**9*(55.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))
      ci(4)=11.0d0*dx1**8*(15.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))
      ci(5)=55.0d0*dx1**7*(6.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(6)=33.0d0*dx1**6*(14.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))
      ci(7)=66.0d0*dx1**5*(7.0d0*dx2**2+dx1*(14.0d0*dx2+5.0d0*dx1))
      ci(8)=66.0d0*dx1**4*(5.0d0*dx2**2+7.0d0*dx1*(2.0d0*dx2+dx1))
      ci(9)=33.0d0*dx1**3*(5.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+7.0d0*dx1))
      ci(10)=55.0d0*dx1**2*(dx2**2+6.0d0*dx1*(dx2+dx1))
      ci(11)=11.0d0*dx1*(dx2**2+5.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(12)=dx2**2+11.0d0*dx1*(2.0d0*dx2+5.0d0*dx1)
      ci(13)=2.0d0*dx2+11.0d0*dx1
      ci(14)=1.0d0
    case(3)
      ci(1)=dx1**11*dx2**3
      ci(2)=dx1**10*dx2**2*(11.0d0*dx2+3.0d0*dx1)
      ci(3)=dx1**9*dx2*(55.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+dx1))
      ci(4)=dx1**8*(165.0d0*dx2**3+dx1*(165.0d0*dx2**2+dx1*(33.0d0*dx2+dx1)))
      ci(5)=11.0d0*dx1**7*(30.0d0*dx2**3+dx1*(45.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))
      ci(6)=11.0d0*dx1**6*(42.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(7)=33.0d0*dx1**5*(14.0d0*dx2**3+dx1*(42.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(8)=66.0d0*dx1**4*(5.0d0*dx2**3+dx1*(21.0d0*dx2**2+dx1*(21.0d0*dx2+5.0d0*dx1)))
      ci(9)=33.0d0*dx1**3*(5.0d0*dx2**3+2.0d0*dx1*(15.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(10)=11.0d0*dx1**2*(5.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+7.0d0*dx1)))
      ci(11)=11.0d0*dx1*(dx2**3+15.0d0*dx1*(dx2**2+dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(12)=dx2**3+33.0d0*dx1*(dx2**2+5.0d0*dx1*(dx2+dx1))
      ci(13)=3.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+5.0d0*dx1)
      ci(14)=3.0d0*dx2+11.0d0*dx1
      ci(15)=1.0d0
    case(4)
      ci(1)=dx1**11*dx2**4
      ci(2)=dx1**10*dx2**3*(11.0d0*dx2+4.0d0*dx1)
      ci(3)=dx1**9*dx2**2*(55.0d0*dx2**2+2.0d0*dx1*(22.0d0*dx2+3.0d0*dx1))
      ci(4)=dx1**8*dx2*(165.0d0*dx2**3+2.0d0*dx1*(110.0d0*dx2**2+dx1*(33.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**7*(330.0d0*dx2**4+dx1*(660.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1*(44.0d0*dx2+dx1))))
      ci(6)=11.0d0*dx1**6*(42.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(7)=11.0d0*dx1**5*(42.0d0*dx2**4+dx1*(168.0d0*dx2**3+5.0d0*dx1*(36.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(8)=33.0d0*dx1**4*(10.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(84.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(9)=33.0d0*dx1**3*(5.0d0*dx2**4+2.0d0*dx1*(20.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(28.0d0*dx2+5.0d0*dx1))))
      ci(10)=11.0d0*dx1**2*(5.0d0*dx2**4+6.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(11)=11.0d0*dx1*(dx2**4+2.0d0*dx1*(10.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+dx1*(20.0d0*dx2+7.0d0*dx1))))
      ci(12)=dx2**4+22.0d0*dx1*(2.0d0*dx2**3+15.0d0*dx1*(dx2**2+dx1*(2.0d0*dx2+dx1)))
      ci(13)=4.0d0*dx2**3+11.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(14)=6.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+5.0d0*dx1)
      ci(15)=4.0d0*dx2+11.0d0*dx1
      ci(16)=1.0d0
    case(5)
      ci(1)=dx1**11*dx2**5
      ci(2)=dx1**10*dx2**4*(11.0d0*dx2+5.0d0*dx1)
      ci(3)=5.0d0*dx1**9*dx2**3*(11.0d0*dx2**2+dx1*(11.0d0*dx2+2.0d0*dx1))
      ci(4)=5.0d0*dx1**8*dx2**2*(33.0d0*dx2**3+dx1*(55.0d0*dx2**2+2.0d0*dx1*(11.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**7*dx2*(66.0d0*dx2**4+dx1*(165.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))))
      ci(6)=dx1**6*(462.0d0*dx2**5+dx1*(1650.0d0*dx2**4+dx1*(1650.0d0*dx2**3+dx1*(550.0d0*dx2**2+dx1*(55.0d0*dx2+dx1)))))
      ci(7)=11.0d0*dx1**5*(42.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(300.0d0*dx2**3+dx1*(150.0d0*dx2**2+dx1*(25.0d0*dx2+dx1)))))
      ci(8)=55.0d0*dx1**4*(6.0d0*dx2**5+dx1*(42.0d0*dx2**4+dx1*(84.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(9)=165.0d0*dx1**3*(dx2**5+dx1*(10.0d0*dx2**4+dx1*(28.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(10.0d0*dx2+dx1)))))
      ci(10)=55.0d0*dx1**2*(dx2**5+3.0d0*dx1*(5.0d0*dx2**4+2.0d0*dx1*(10.0d0*dx2**3+dx1*(14.0d0*dx2**2+dx1*(7.0d0*dx2+dx1)))))
      ci(11)=11.0d0*dx1*(dx2**5+dx1*(25.0d0*dx2**4+6.0d0*dx1*(25.0d0*dx2**3+dx1*(50.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+dx1)))))
      ci(12)=dx2**5+11.0d0*dx1*(5.0d0*dx2**4+2.0d0*dx1*(25.0d0*dx2**3+3.0d0*dx1*(25.0d0*dx2**2+dx1*(25.0d0*dx2+7.0d0*dx1))))
      ci(13)=5.0d0*(dx2**4+11.0d0*dx1*(2.0d0*dx2**3+dx1*(10.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))))
      ci(14)=5.0d0*(2.0d0*dx2**3+11.0d0*dx1*(2.0d0*dx2**2+dx1*(5.0d0*dx2+3.0d0*dx1)))
      ci(15)=5.0d0*(2.0d0*dx2**2+11.0d0*dx1*(dx2+dx1))
      ci(16)=5.0d0*dx2+11.0d0*dx1
      ci(17)=1.0d0
    case(6)
      ci(1)=dx1**11*dx2**6
      ci(2)=dx1**10*dx2**5*(11.0d0*dx2+6.0d0*dx1)
      ci(3)=dx1**9*dx2**4*(55.0d0*dx2**2+3.0d0*dx1*(22.0d0*dx2+5.0d0*dx1))
      ci(4)=5.0d0*dx1**8*dx2**3*(33.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(33.0d0*dx2+4.0d0*dx1)))
      ci(5)=5.0d0*dx1**7*dx2**2*(66.0d0*dx2**4+dx1*(198.0d0*dx2**3+dx1*(165.0d0*dx2**2+dx1*(44.0d0*dx2+3.0d0*dx1))))
      ci(6)=dx1**6*dx2*(462.0d0*dx2**5+dx1*(1980.0d0*dx2**4+dx1*(2475.0d0*dx2**3+dx1*(1100.0d0*dx2**2+3.0d0*dx1*(55.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=dx1**5*(462.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(4950.0d0*dx2**4+dx1*(3300.0d0*dx2**3+dx1*(825.0d0*dx2**2+dx1* &
(66.0d0*dx2+dx1))))))
      ci(8)=11.0d0*dx1**4*(30.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1*(600.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1* &
(30.0d0*dx2+dx1))))))
      ci(9)=55.0d0*dx1**3*(3.0d0*dx2**6+dx1*(36.0d0*dx2**5+dx1*(126.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1* &
(18.0d0*dx2+dx1))))))
      ci(10)=55.0d0*dx1**2*(dx2**6+3.0d0*dx1*(6.0d0*dx2**5+dx1*(30.0d0*dx2**4+dx1*(56.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1* &
(12.0d0*dx2+dx1))))))
      ci(11)=11.0d0*dx1*(dx2**6+3.0d0*dx1*(10.0d0*dx2**5+dx1*(75.0d0*dx2**4+2.0d0*dx1*(100.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1* &
(42.0d0*dx2+5.0d0*dx1))))))
      ci(12)=dx2**6+33.0d0*dx1*(2.0d0*dx2**5+dx1*(25.0d0*dx2**4+2.0d0*dx1*(50.0d0*dx2**3+dx1*(75.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2 &
+dx1)))))
      ci(13)=6.0d0*dx2**5+11.0d0*dx1*(15.0d0*dx2**4+dx1*(100.0d0*dx2**3+3.0d0*dx1*(75.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2+7.0d0*dx1)) &
))
      ci(14)=5.0d0*(3.0d0*dx2**4+11.0d0*dx1*(4.0d0*dx2**3+3.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))))
      ci(15)=5.0d0*(4.0d0*dx2**3+33.0d0*dx1*(dx2**2+dx1*(2.0d0*dx2+dx1)))
      ci(16)=15.0d0*dx2**2+11.0d0*dx1*(6.0d0*dx2+5.0d0*dx1)
      ci(17)=6.0d0*dx2+11.0d0*dx1
      ci(18)=1.0d0
    case(7)
      ci(1)=dx1**11*dx2**7
      ci(2)=dx1**10*dx2**6*(11.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**9*dx2**5*(55.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+3.0d0*dx1))
      ci(4)=dx1**8*dx2**4*(165.0d0*dx2**3+7.0d0*dx1*(55.0d0*dx2**2+dx1*(33.0d0*dx2+5.0d0*dx1)))
      ci(5)=5.0d0*dx1**7*dx2**3*(66.0d0*dx2**4+7.0d0*dx1*(33.0d0*dx2**3+dx1*(33.0d0*dx2**2+dx1*(11.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**6*dx2**2*(66.0d0*dx2**5+dx1*(330.0d0*dx2**4+dx1*(495.0d0*dx2**3+dx1*(275.0d0*dx2**2+dx1*(55.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=7.0d0*dx1**5*dx2*(66.0d0*dx2**6+dx1*(462.0d0*dx2**5+dx1*(990.0d0*dx2**4+dx1*(825.0d0*dx2**3+dx1*(275.0d0*dx2**2+dx1* &
(33.0d0*dx2+dx1))))))
      ci(8)=dx1**4*(330.0d0*dx2**7+dx1*(3234.0d0*dx2**6+dx1*(9702.0d0*dx2**5+dx1*(11550.0d0*dx2**4+dx1*(5775.0d0*dx2**3+dx1* &
(1155.0d0*dx2**2+dx1*(77.0d0*dx2+dx1)))))))
      ci(9)=11.0d0*dx1**3*(15.0d0*dx2**7+dx1*(210.0d0*dx2**6+dx1*(882.0d0*dx2**5+dx1*(1470.0d0*dx2**4+dx1*(1050.0d0*dx2**3+dx1* &
(315.0d0*dx2**2+dx1*(35.0d0*dx2+dx1)))))))
      ci(10)=55.0d0*dx1**2*(dx2**7+dx1*(21.0d0*dx2**6+dx1*(126.0d0*dx2**5+dx1*(294.0d0*dx2**4+dx1*(294.0d0*dx2**3+dx1* &
(126.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)))))))
      ci(11)=11.0d0*dx1*(dx2**7+dx1*(35.0d0*dx2**6+3.0d0*dx1*(105.0d0*dx2**5+dx1*(350.0d0*dx2**4+dx1*(490.0d0*dx2**3+dx1* &
(294.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+dx1)))))))
      ci(12)=dx2**7+11.0d0*dx1*(7.0d0*dx2**6+3.0d0*dx1*(35.0d0*dx2**5+dx1*(175.0d0*dx2**4+2.0d0*dx1*(175.0d0*dx2**3+dx1* &
(147.0d0*dx2**2+dx1*(49.0d0*dx2+5.0d0*dx1))))))
      ci(13)=7.0d0*(dx2**6+11.0d0*dx1*(3.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(25.0d0*dx2**3+2.0d0*dx1*(15.0d0*dx2**2+dx1* &
(7.0d0*dx2+dx1))))))
      ci(14)=7.0d0*(3.0d0*dx2**5+11.0d0*dx1*(5.0d0*dx2**4+dx1*(25.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1)))) &
)
      ci(15)=5.0d0*(7.0d0*dx2**4+11.0d0*dx1*(7.0d0*dx2**3+3.0d0*dx1*(7.0d0*dx2**2+dx1*(7.0d0*dx2+2.0d0*dx1))))
      ci(16)=35.0d0*dx2**3+11.0d0*dx1*(21.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+3.0d0*dx1))
      ci(17)=21.0d0*dx2**2+11.0d0*dx1*(7.0d0*dx2+5.0d0*dx1)
      ci(18)=7.0d0*dx2+11.0d0*dx1
      ci(19)=1.0d0
    case(8)
      ci(1)=dx1**11*dx2**8
      ci(2)=dx1**10*dx2**7*(11.0d0*dx2+8.0d0*dx1)
      ci(3)=dx1**9*dx2**6*(55.0d0*dx2**2+4.0d0*dx1*(22.0d0*dx2+7.0d0*dx1))
      ci(4)=dx1**8*dx2**5*(165.0d0*dx2**3+4.0d0*dx1*(110.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+2.0d0*dx1)))
      ci(5)=2.0d0*dx1**7*dx2**4*(165.0d0*dx2**4+dx1*(660.0d0*dx2**3+7.0d0*dx1*(110.0d0*dx2**2+dx1*(44.0d0*dx2+5.0d0*dx1))))
      ci(6)=2.0d0*dx1**6*dx2**3*(231.0d0*dx2**5+dx1*(1320.0d0*dx2**4+7.0d0*dx1*(330.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(55.0d0*dx2 &
+4.0d0*dx1)))))
      ci(7)=14.0d0*dx1**5*dx2**2*(33.0d0*dx2**6+dx1*(264.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(660.0d0*dx2**3+dx1*(275.0d0*dx2**2 &
+2.0d0*dx1*(22.0d0*dx2+dx1))))))
      ci(8)=2.0d0*dx1**4*dx2*(165.0d0*dx2**7+dx1*(1848.0d0*dx2**6+dx1*(6468.0d0*dx2**5+dx1*(9240.0d0*dx2**4+dx1*(5775.0d0*dx2**3 &
+2.0d0*dx1*(770.0d0*dx2**2+dx1*(77.0d0*dx2+2.0d0*dx1)))))))
      ci(9)=dx1**3*(165.0d0*dx2**8+dx1*(2640.0d0*dx2**7+dx1*(12936.0d0*dx2**6+dx1*(25872.0d0*dx2**5+dx1*(23100.0d0*dx2**4+dx1* &
(9240.0d0*dx2**3+dx1*(1540.0d0*dx2**2+dx1*(88.0d0*dx2+dx1))))))))
      ci(10)=11.0d0*dx1**2*(5.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(840.0d0*dx2**6+dx1*(2352.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(1680.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))))))
      ci(11)=11.0d0*dx1*(dx2**8+dx1*(40.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(1680.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(2352.0d0*dx2**3+5.0d0*dx1*(168.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))))))
      ci(12)=dx2**8+11.0d0*dx1*(8.0d0*dx2**7+dx1*(140.0d0*dx2**6+3.0d0*dx1*(280.0d0*dx2**5+dx1*(700.0d0*dx2**4+dx1*(784.0d0*dx2**3 &
+dx1*(392.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+dx1)))))))
      ci(13)=2.0d0*(4.0d0*dx2**7+11.0d0*dx1*(14.0d0*dx2**6+dx1*(140.0d0*dx2**5+3.0d0*dx1*(175.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1* &
(196.0d0*dx2**2+dx1*(56.0d0*dx2+5.0d0*dx1)))))))
      ci(14)=14.0d0*(2.0d0*dx2**6+11.0d0*dx1*(4.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(20.0d0*dx2**3+dx1*(20.0d0*dx2**2+dx1* &
(8.0d0*dx2+dx1))))))
      ci(15)=2.0d0*(28.0d0*dx2**5+11.0d0*dx1*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+3.0d0*dx1*(70.0d0*dx2**2+dx1*(40.0d0*dx2+7.0d0*dx1 &
)))))
      ci(16)=2.0d0*(35.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+5.0d0*dx1*(14.0d0*dx2**2+3.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(17)=56.0d0*dx2**3+11.0d0*dx1*(28.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(18)=28.0d0*dx2**2+11.0d0*dx1*(8.0d0*dx2+5.0d0*dx1)
      ci(19)=8.0d0*dx2+11.0d0*dx1
      ci(20)=1.0d0
    case(9)
      ci(1)=dx1**11*dx2**9
      ci(2)=dx1**10*dx2**8*(11.0d0*dx2+9.0d0*dx1)
      ci(3)=dx1**9*dx2**7*(55.0d0*dx2**2+9.0d0*dx1*(11.0d0*dx2+4.0d0*dx1))
      ci(4)=3.0d0*dx1**8*dx2**6*(55.0d0*dx2**3+dx1*(165.0d0*dx2**2+4.0d0*dx1*(33.0d0*dx2+7.0d0*dx1)))
      ci(5)=3.0d0*dx1**7*dx2**5*(110.0d0*dx2**4+dx1*(495.0d0*dx2**3+2.0d0*dx1*(330.0d0*dx2**2+7.0d0*dx1*(22.0d0*dx2+3.0d0*dx1))))
      ci(6)=6.0d0*dx1**6*dx2**4*(77.0d0*dx2**5+dx1*(495.0d0*dx2**4+dx1*(990.0d0*dx2**3+7.0d0*dx1*(110.0d0*dx2**2+3.0d0*dx1* &
(11.0d0*dx2+dx1)))))
      ci(7)=6.0d0*dx1**5*dx2**3*(77.0d0*dx2**6+dx1*(693.0d0*dx2**5+dx1*(1980.0d0*dx2**4+7.0d0*dx1*(330.0d0*dx2**3+dx1* &
(165.0d0*dx2**2+dx1*(33.0d0*dx2+2.0d0*dx1))))))
      ci(8)=6.0d0*dx1**4*dx2**2*(55.0d0*dx2**7+dx1*(693.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(4620.0d0*dx2**4+dx1*(3465.0d0*dx2**3 &
+dx1*(1155.0d0*dx2**2+2.0d0*dx1*(77.0d0*dx2+3.0d0*dx1)))))))
      ci(9)=3.0d0*dx1**3*dx2*(55.0d0*dx2**8+dx1*(990.0d0*dx2**7+dx1*(5544.0d0*dx2**6+dx1*(12936.0d0*dx2**5+dx1*(13860.0d0*dx2**4 &
+dx1*(6930.0d0*dx2**3+dx1*(1540.0d0*dx2**2+3.0d0*dx1*(44.0d0*dx2+dx1))))))))
      ci(10)=dx1**2*(55.0d0*dx2**9+dx1*(1485.0d0*dx2**8+dx1*(11880.0d0*dx2**7+dx1*(38808.0d0*dx2**6+dx1*(58212.0d0*dx2**5+dx1* &
(41580.0d0*dx2**4+dx1*(13860.0d0*dx2**3+dx1*(1980.0d0*dx2**2+dx1*(99.0d0*dx2+dx1)))))))))
      ci(11)=11.0d0*dx1*(dx2**9+dx1*(45.0d0*dx2**8+dx1*(540.0d0*dx2**7+dx1*(2520.0d0*dx2**6+dx1*(5292.0d0*dx2**5+dx1* &
(5292.0d0*dx2**4+dx1*(2520.0d0*dx2**3+dx1*(540.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))))))))
      ci(12)=dx2**9+11.0d0*dx1*(9.0d0*dx2**8+dx1*(180.0d0*dx2**7+dx1*(1260.0d0*dx2**6+dx1*(3780.0d0*dx2**5+dx1*(5292.0d0*dx2**4 &
+dx1*(3528.0d0*dx2**3+5.0d0*dx1*(216.0d0*dx2**2+dx1*(27.0d0*dx2+dx1))))))))
      ci(13)=3.0d0*(3.0d0*dx2**8+11.0d0*dx1*(12.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(630.0d0*dx2**5+dx1*(1260.0d0*dx2**4+dx1* &
(1176.0d0*dx2**3+dx1*(504.0d0*dx2**2+5.0d0*dx1*(18.0d0*dx2+dx1))))))))
      ci(14)=6.0d0*(6.0d0*dx2**7+11.0d0*dx1*(14.0d0*dx2**6+dx1*(105.0d0*dx2**5+dx1*(315.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1* &
(252.0d0*dx2**2+dx1*(63.0d0*dx2+5.0d0*dx1)))))))
      ci(15)=6.0d0*(14.0d0*dx2**6+11.0d0*dx1*(21.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(210.0d0*dx2**3+dx1*(180.0d0*dx2**2+7.0d0*dx1* &
(9.0d0*dx2+dx1))))))
      ci(16)=6.0d0*(21.0d0*dx2**5+11.0d0*dx1*(21.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(45.0d0*dx2+7.0d0*dx1)))))
      ci(17)=3.0d0*(42.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+5.0d0*dx1*(12.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1))))
      ci(18)=3.0d0*(28.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(19)=36.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+5.0d0*dx1)
      ci(20)=9.0d0*dx2+11.0d0*dx1
      ci(21)=1.0d0
    case(10)
      ci(1)=dx1**11*dx2**10
      ci(2)=dx1**10*dx2**9*(11.0d0*dx2+10.0d0*dx1)
      ci(3)=5.0d0*dx1**9*dx2**8*(11.0d0*dx2**2+dx1*(22.0d0*dx2+9.0d0*dx1))
      ci(4)=5.0d0*dx1**8*dx2**7*(33.0d0*dx2**3+dx1*(110.0d0*dx2**2+3.0d0*dx1*(33.0d0*dx2+8.0d0*dx1)))
      ci(5)=15.0d0*dx1**7*dx2**6*(22.0d0*dx2**4+dx1*(110.0d0*dx2**3+dx1*(165.0d0*dx2**2+2.0d0*dx1*(44.0d0*dx2+7.0d0*dx1))))
      ci(6)=3.0d0*dx1**6*dx2**5*(154.0d0*dx2**5+dx1*(1100.0d0*dx2**4+dx1*(2475.0d0*dx2**3+2.0d0*dx1*(1100.0d0*dx2**2+7.0d0*dx1* &
(55.0d0*dx2+6.0d0*dx1)))))
      ci(7)=6.0d0*dx1**5*dx2**4*(77.0d0*dx2**6+dx1*(770.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(3300.0d0*dx2**3+7.0d0*dx1* &
(275.0d0*dx2**2+dx1*(66.0d0*dx2+5.0d0*dx1))))))
      ci(8)=30.0d0*dx1**4*dx2**3*(11.0d0*dx2**7+dx1*(154.0d0*dx2**6+dx1*(693.0d0*dx2**5+dx1*(1320.0d0*dx2**4+dx1*(1155.0d0*dx2**3 &
+dx1*(462.0d0*dx2**2+dx1*(77.0d0*dx2+4.0d0*dx1)))))))
      ci(9)=15.0d0*dx1**3*dx2**2*(11.0d0*dx2**8+dx1*(220.0d0*dx2**7+dx1*(1386.0d0*dx2**6+dx1*(3696.0d0*dx2**5+dx1*(4620.0d0*dx2**4 &
+dx1*(2772.0d0*dx2**3+dx1*(770.0d0*dx2**2+dx1*(88.0d0*dx2+3.0d0*dx1))))))))
      ci(10)=5.0d0*dx1**2*dx2*(11.0d0*dx2**9+dx1*(330.0d0*dx2**8+dx1*(2970.0d0*dx2**7+dx1*(11088.0d0*dx2**6+dx1*(19404.0d0*dx2**5 &
+dx1*(16632.0d0*dx2**4+dx1*(6930.0d0*dx2**3+dx1*(1320.0d0*dx2**2+dx1*(99.0d0*dx2+2.0d0*dx1)))))))))
      ci(11)=dx1*(11.0d0*dx2**10+dx1*(550.0d0*dx2**9+dx1*(7425.0d0*dx2**8+dx1*(39600.0d0*dx2**7+dx1*(97020.0d0*dx2**6+dx1* &
(116424.0d0*dx2**5+dx1*(69300.0d0*dx2**4+dx1*(19800.0d0*dx2**3+dx1*(2475.0d0*dx2**2+dx1*(110.0d0*dx2+dx1))))))))))
      ci(12)=dx2**10+11.0d0*dx1*(10.0d0*dx2**9+dx1*(225.0d0*dx2**8+dx1*(1800.0d0*dx2**7+dx1*(6300.0d0*dx2**6+dx1*(10584.0d0*dx2**5 &
+dx1*(8820.0d0*dx2**4+dx1*(3600.0d0*dx2**3+dx1*(675.0d0*dx2**2+dx1*(50.0d0*dx2+dx1)))))))))
      ci(13)=5.0d0*(2.0d0*dx2**9+11.0d0*dx1*(9.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(630.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1* &
(1764.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(270.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))))))
      ci(14)=15.0d0*(3.0d0*dx2**8+11.0d0*dx1*(8.0d0*dx2**7+dx1*(70.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1* &
(336.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))))))
      ci(15)=30.0d0*(4.0d0*dx2**7+11.0d0*dx1*(7.0d0*dx2**6+dx1*(42.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1* &
(63.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)))))))
      ci(16)=6.0d0*(35.0d0*dx2**6+11.0d0*dx1*(42.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(300.0d0*dx2**3+dx1*(225.0d0*dx2**2+7.0d0*dx1* &
(10.0d0*dx2+dx1))))))
      ci(17)=3.0d0*(84.0d0*dx2**5+11.0d0*dx1*(70.0d0*dx2**4+dx1*(200.0d0*dx2**3+dx1*(225.0d0*dx2**2+2.0d0*dx1*(50.0d0*dx2 &
+7.0d0*dx1)))))
      ci(18)=15.0d0*(14.0d0*dx2**4+11.0d0*dx1*(8.0d0*dx2**3+dx1*(15.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1))))
      ci(19)=5.0d0*(24.0d0*dx2**3+11.0d0*dx1*(9.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1)))
      ci(20)=5.0d0*(9.0d0*dx2**2+11.0d0*dx1*(2.0d0*dx2+dx1))
      ci(21)=10.0d0*dx2+11.0d0*dx1
      ci(22)=1.0d0
    case(11)
      ci(1)=dx1**11*dx2**11
      ci(2)=11.0d0*dx1**10*dx2**10*(dx2+dx1)
      ci(3)=11.0d0*dx1**9*dx2**9*(5.0d0*dx2**2+dx1*(11.0d0*dx2+5.0d0*dx1))
      ci(4)=55.0d0*dx1**8*dx2**8*(3.0d0*dx2**3+dx1*(11.0d0*dx2**2+dx1*(11.0d0*dx2+3.0d0*dx1)))
      ci(5)=55.0d0*dx1**7*dx2**7*(6.0d0*dx2**4+dx1*(33.0d0*dx2**3+dx1*(55.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+2.0d0*dx1))))
      ci(6)=33.0d0*dx1**6*dx2**6*(14.0d0*dx2**5+dx1*(110.0d0*dx2**4+dx1*(275.0d0*dx2**3+dx1*(275.0d0*dx2**2+2.0d0*dx1*(55.0d0*dx2 &
+7.0d0*dx1)))))
      ci(7)=33.0d0*dx1**5*dx2**5*(14.0d0*dx2**6+dx1*(154.0d0*dx2**5+dx1*(550.0d0*dx2**4+dx1*(825.0d0*dx2**3+2.0d0*dx1* &
(275.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+dx1))))))
      ci(8)=66.0d0*dx1**4*dx2**4*(5.0d0*dx2**7+dx1*(77.0d0*dx2**6+dx1*(385.0d0*dx2**5+dx1*(825.0d0*dx2**4+dx1*(825.0d0*dx2**3+dx1* &
(385.0d0*dx2**2+dx1*(77.0d0*dx2+5.0d0*dx1)))))))
      ci(9)=165.0d0*dx1**3*dx2**3*(dx2**8+dx1*(22.0d0*dx2**7+dx1*(154.0d0*dx2**6+dx1*(462.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1* &
(462.0d0*dx2**3+dx1*(154.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))))))))
      ci(10)=55.0d0*dx1**2*dx2**2*(dx2**9+dx1*(33.0d0*dx2**8+dx1*(330.0d0*dx2**7+dx1*(1386.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1* &
(2772.0d0*dx2**4+dx1*(1386.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1*(33.0d0*dx2+dx1)))))))))
      ci(11)=11.0d0*dx1*dx2*(dx2**10+dx1*(55.0d0*dx2**9+dx1*(825.0d0*dx2**8+dx1*(4950.0d0*dx2**7+dx1*(13860.0d0*dx2**6+dx1* &
(19404.0d0*dx2**5+dx1*(13860.0d0*dx2**4+dx1*(4950.0d0*dx2**3+dx1*(825.0d0*dx2**2+dx1*(55.0d0*dx2+dx1))))))))))
      ci(12)=dx2**11+dx1*(121.0d0*dx2**10+dx1*(3025.0d0*dx2**9+dx1*(27225.0d0*dx2**8+dx1*(108900.0d0*dx2**7+dx1*(213444.0d0*dx2**6 &
+dx1*(213444.0d0*dx2**5+dx1*(108900.0d0*dx2**4+dx1*(27225.0d0*dx2**3+dx1*(3025.0d0*dx2**2+dx1*(121.0d0*dx2+dx1))))))))))
      ci(13)=11.0d0*(dx2**10+dx1*(55.0d0*dx2**9+dx1*(825.0d0*dx2**8+dx1*(4950.0d0*dx2**7+dx1*(13860.0d0*dx2**6+dx1* &
(19404.0d0*dx2**5+dx1*(13860.0d0*dx2**4+dx1*(4950.0d0*dx2**3+dx1*(825.0d0*dx2**2+dx1*(55.0d0*dx2+dx1))))))))))
      ci(14)=55.0d0*(dx2**9+dx1*(33.0d0*dx2**8+dx1*(330.0d0*dx2**7+dx1*(1386.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(2772.0d0*dx2**4 &
+dx1*(1386.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1*(33.0d0*dx2+dx1)))))))))
      ci(15)=165.0d0*(dx2**8+dx1*(22.0d0*dx2**7+dx1*(154.0d0*dx2**6+dx1*(462.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(462.0d0*dx2**3 &
+dx1*(154.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))))))))
      ci(16)=66.0d0*(5.0d0*dx2**7+dx1*(77.0d0*dx2**6+dx1*(385.0d0*dx2**5+dx1*(825.0d0*dx2**4+dx1*(825.0d0*dx2**3+dx1* &
(385.0d0*dx2**2+dx1*(77.0d0*dx2+5.0d0*dx1)))))))
      ci(17)=33.0d0*(14.0d0*dx2**6+dx1*(154.0d0*dx2**5+dx1*(550.0d0*dx2**4+dx1*(825.0d0*dx2**3+2.0d0*dx1*(275.0d0*dx2**2 &
+7.0d0*dx1*(11.0d0*dx2+dx1))))))
      ci(18)=33.0d0*(14.0d0*dx2**5+dx1*(110.0d0*dx2**4+dx1*(275.0d0*dx2**3+dx1*(275.0d0*dx2**2+2.0d0*dx1*(55.0d0*dx2+7.0d0*dx1)))) &
)
      ci(19)=55.0d0*(6.0d0*dx2**4+dx1*(33.0d0*dx2**3+dx1*(55.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+2.0d0*dx1))))
      ci(20)=55.0d0*(3.0d0*dx2**3+dx1*(11.0d0*dx2**2+dx1*(11.0d0*dx2+3.0d0*dx1)))
      ci(21)=11.0d0*(5.0d0*dx2**2+dx1*(11.0d0*dx2+5.0d0*dx1))
      ci(22)=11.0d0*(dx2+dx1)
      ci(23)=1.0d0
    case(12)
      ci(1)=dx1**11*dx2**12
      ci(2)=dx1**10*dx2**11*(11.0d0*dx2+12.0d0*dx1)
      ci(3)=11.0d0*dx1**9*dx2**10*(5.0d0*dx2**2+6.0d0*dx1*(2.0d0*dx2+dx1))
      ci(4)=11.0d0*dx1**8*dx2**9*(15.0d0*dx2**3+2.0d0*dx1*(30.0d0*dx2**2+dx1*(33.0d0*dx2+10.0d0*dx1)))
      ci(5)=55.0d0*dx1**7*dx2**8*(6.0d0*dx2**4+dx1*(36.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(44.0d0*dx2+9.0d0*dx1))))
      ci(6)=11.0d0*dx1**6*dx2**7*(42.0d0*dx2**5+dx1*(360.0d0*dx2**4+dx1*(990.0d0*dx2**3+dx1*(1100.0d0*dx2**2+9.0d0*dx1*(55.0d0*dx2 &
+8.0d0*dx1)))))
      ci(7)=33.0d0*dx1**5*dx2**6*(14.0d0*dx2**6+dx1*(168.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(1100.0d0*dx2**3+dx1*(825.0d0*dx2**2 &
+4.0d0*dx1*(66.0d0*dx2+7.0d0*dx1))))))
      ci(8)=33.0d0*dx1**4*dx2**5*(10.0d0*dx2**7+dx1*(168.0d0*dx2**6+dx1*(924.0d0*dx2**5+dx1*(2200.0d0*dx2**4+dx1*(2475.0d0*dx2**3 &
+4.0d0*dx1*(330.0d0*dx2**2+dx1*(77.0d0*dx2+6.0d0*dx1)))))))
      ci(9)=33.0d0*dx1**3*dx2**4*(5.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(3080.0d0*dx2**5+dx1*(4950.0d0*dx2**4 &
+dx1*(3960.0d0*dx2**3+dx1*(1540.0d0*dx2**2+3.0d0*dx1*(88.0d0*dx2+5.0d0*dx1))))))))
      ci(10)=55.0d0*dx1**2*dx2**3*(dx2**9+dx1*(36.0d0*dx2**8+dx1*(396.0d0*dx2**7+dx1*(1848.0d0*dx2**6+dx1*(4158.0d0*dx2**5+dx1* &
(4752.0d0*dx2**4+dx1*(2772.0d0*dx2**3+dx1*(792.0d0*dx2**2+dx1*(99.0d0*dx2+4.0d0*dx1)))))))))
      ci(11)=11.0d0*dx1*dx2**2*(dx2**10+dx1*(60.0d0*dx2**9+dx1*(990.0d0*dx2**8+dx1*(6600.0d0*dx2**7+dx1*(20790.0d0*dx2**6+dx1* &
(33264.0d0*dx2**5+dx1*(27720.0d0*dx2**4+dx1*(11880.0d0*dx2**3+dx1*(2475.0d0*dx2**2+2.0d0*dx1*(110.0d0*dx2+3.0d0*dx1))))))))))
      ci(12)=dx2*(dx2**11+dx1*(132.0d0*dx2**10+dx1*(3630.0d0*dx2**9+dx1*(36300.0d0*dx2**8+dx1*(163350.0d0*dx2**7+dx1* &
(365904.0d0*dx2**6+dx1*(426888.0d0*dx2**5+dx1*(261360.0d0*dx2**4+dx1*(81675.0d0*dx2**3+2.0d0*dx1*(6050.0d0*dx2**2+3.0d0*dx1* &
(121.0d0*dx2+2.0d0*dx1)))))))))))
      ci(13)=12.0d0*dx2**11+dx1*(726.0d0*dx2**10+dx1*(12100.0d0*dx2**9+dx1*(81675.0d0*dx2**8+dx1*(261360.0d0*dx2**7+dx1* &
(426888.0d0*dx2**6+dx1*(365904.0d0*dx2**5+dx1*(163350.0d0*dx2**4+dx1*(36300.0d0*dx2**3+dx1*(3630.0d0*dx2**2+dx1*(132.0d0*dx2+dx1)) &
))))))))
      ci(14)=11.0d0*(6.0d0*dx2**10+dx1*(220.0d0*dx2**9+dx1*(2475.0d0*dx2**8+dx1*(11880.0d0*dx2**7+dx1*(27720.0d0*dx2**6+dx1* &
(33264.0d0*dx2**5+dx1*(20790.0d0*dx2**4+dx1*(6600.0d0*dx2**3+dx1*(990.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))))))))
      ci(15)=55.0d0*(4.0d0*dx2**9+dx1*(99.0d0*dx2**8+dx1*(792.0d0*dx2**7+dx1*(2772.0d0*dx2**6+dx1*(4752.0d0*dx2**5+dx1* &
(4158.0d0*dx2**4+dx1*(1848.0d0*dx2**3+dx1*(396.0d0*dx2**2+dx1*(36.0d0*dx2+dx1)))))))))
      ci(16)=33.0d0*(15.0d0*dx2**8+dx1*(264.0d0*dx2**7+dx1*(1540.0d0*dx2**6+dx1*(3960.0d0*dx2**5+dx1*(4950.0d0*dx2**4+dx1* &
(3080.0d0*dx2**3+dx1*(924.0d0*dx2**2+5.0d0*dx1*(24.0d0*dx2+dx1))))))))
      ci(17)=33.0d0*(24.0d0*dx2**7+dx1*(308.0d0*dx2**6+dx1*(1320.0d0*dx2**5+dx1*(2475.0d0*dx2**4+2.0d0*dx1*(1100.0d0*dx2**3+dx1* &
(462.0d0*dx2**2+dx1*(84.0d0*dx2+5.0d0*dx1)))))))
      ci(18)=33.0d0*(28.0d0*dx2**6+dx1*(264.0d0*dx2**5+dx1*(825.0d0*dx2**4+2.0d0*dx1*(550.0d0*dx2**3+dx1*(330.0d0*dx2**2 &
+7.0d0*dx1*(12.0d0*dx2+dx1))))))
      ci(19)=11.0d0*(72.0d0*dx2**5+dx1*(495.0d0*dx2**4+2.0d0*dx1*(550.0d0*dx2**3+3.0d0*dx1*(165.0d0*dx2**2+dx1*(60.0d0*dx2 &
+7.0d0*dx1)))))
      ci(20)=55.0d0*(9.0d0*dx2**4+2.0d0*dx1*(22.0d0*dx2**3+3.0d0*dx1*(11.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))
      ci(21)=11.0d0*(20.0d0*dx2**3+3.0d0*dx1*(22.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1)))
      ci(22)=11.0d0*(6.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(23)=12.0d0*dx2+11.0d0*dx1
      ci(24)=1.0d0
    case(13)
      ci(1)=dx1**11*dx2**13
      ci(2)=dx1**10*dx2**12*(11.0d0*dx2+13.0d0*dx1)
      ci(3)=dx1**9*dx2**11*(55.0d0*dx2**2+13.0d0*dx1*(11.0d0*dx2+6.0d0*dx1))
      ci(4)=11.0d0*dx1**8*dx2**10*(15.0d0*dx2**3+13.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(5)=11.0d0*dx1**7*dx2**9*(30.0d0*dx2**4+13.0d0*dx1*(15.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(22.0d0*dx2+5.0d0*dx1))))
      ci(6)=11.0d0*dx1**6*dx2**8*(42.0d0*dx2**5+13.0d0*dx1*(30.0d0*dx2**4+dx1*(90.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(55.0d0*dx2 &
+9.0d0*dx1)))))
      ci(7)=11.0d0*dx1**5*dx2**7*(42.0d0*dx2**6+13.0d0*dx1*(42.0d0*dx2**5+dx1*(180.0d0*dx2**4+dx1*(330.0d0*dx2**3+dx1* &
(275.0d0*dx2**2+3.0d0*dx1*(33.0d0*dx2+4.0d0*dx1))))))
      ci(8)=33.0d0*dx1**4*dx2**6*(10.0d0*dx2**7+13.0d0*dx1*(14.0d0*dx2**6+dx1*(84.0d0*dx2**5+dx1*(220.0d0*dx2**4+dx1* &
(275.0d0*dx2**3+dx1*(165.0d0*dx2**2+4.0d0*dx1*(11.0d0*dx2+dx1)))))))
      ci(9)=33.0d0*dx1**3*dx2**5*(5.0d0*dx2**8+13.0d0*dx1*(10.0d0*dx2**7+dx1*(84.0d0*dx2**6+dx1*(308.0d0*dx2**5+dx1* &
(550.0d0*dx2**4+dx1*(495.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(44.0d0*dx2+3.0d0*dx1))))))))
      ci(10)=11.0d0*dx1**2*dx2**4*(5.0d0*dx2**9+13.0d0*dx1*(15.0d0*dx2**8+dx1*(180.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1* &
(2310.0d0*dx2**5+dx1*(2970.0d0*dx2**4+dx1*(1980.0d0*dx2**3+dx1*(660.0d0*dx2**2+dx1*(99.0d0*dx2+5.0d0*dx1)))))))))
      ci(11)=11.0d0*dx1*dx2**3*(dx2**10+13.0d0*dx1*(5.0d0*dx2**9+dx1*(90.0d0*dx2**8+dx1*(660.0d0*dx2**7+dx1*(2310.0d0*dx2**6+dx1* &
(4158.0d0*dx2**5+dx1*(3960.0d0*dx2**4+dx1*(1980.0d0*dx2**3+dx1*(495.0d0*dx2**2+dx1*(55.0d0*dx2+2.0d0*dx1))))))))))
      ci(12)=dx2**2*(dx2**11+13.0d0*dx1*(11.0d0*dx2**10+dx1*(330.0d0*dx2**9+dx1*(3630.0d0*dx2**8+dx1*(18150.0d0*dx2**7+dx1* &
(45738.0d0*dx2**6+dx1*(60984.0d0*dx2**5+dx1*(43560.0d0*dx2**4+dx1*(16335.0d0*dx2**3+dx1*(3025.0d0*dx2**2+2.0d0*dx1*(121.0d0*dx2 &
+3.0d0*dx1)))))))))))
      ci(13)=13.0d0*dx2*(dx2**11+dx1*(66.0d0*dx2**10+dx1*(1210.0d0*dx2**9+dx1*(9075.0d0*dx2**8+dx1*(32670.0d0*dx2**7+dx1* &
(60984.0d0*dx2**6+dx1*(60984.0d0*dx2**5+dx1*(32670.0d0*dx2**4+dx1*(9075.0d0*dx2**3+dx1*(1210.0d0*dx2**2+dx1*(66.0d0*dx2+dx1))))))) &
))))
      ci(14)=78.0d0*dx2**11+dx1*(3146.0d0*dx2**10+dx1*(39325.0d0*dx2**9+dx1*(212355.0d0*dx2**8+dx1*(566280.0d0*dx2**7+dx1* &
(792792.0d0*dx2**6+dx1*(594594.0d0*dx2**5+dx1*(235950.0d0*dx2**4+dx1*(47190.0d0*dx2**3+dx1*(4290.0d0*dx2**2+dx1*(143.0d0*dx2+dx1)) &
))))))))
      ci(15)=11.0d0*(26.0d0*dx2**10+dx1*(715.0d0*dx2**9+dx1*(6435.0d0*dx2**8+dx1*(25740.0d0*dx2**7+dx1*(51480.0d0*dx2**6+dx1* &
(54054.0d0*dx2**5+dx1*(30030.0d0*dx2**4+dx1*(8580.0d0*dx2**3+dx1*(1170.0d0*dx2**2+dx1*(65.0d0*dx2+dx1))))))))))
      ci(16)=11.0d0*(65.0d0*dx2**9+dx1*(1287.0d0*dx2**8+dx1*(8580.0d0*dx2**7+dx1*(25740.0d0*dx2**6+dx1*(38610.0d0*dx2**5+dx1* &
(30030.0d0*dx2**4+dx1*(12012.0d0*dx2**3+5.0d0*dx1*(468.0d0*dx2**2+dx1*(39.0d0*dx2+dx1)))))))))
      ci(17)=33.0d0*(39.0d0*dx2**8+dx1*(572.0d0*dx2**7+dx1*(2860.0d0*dx2**6+dx1*(6435.0d0*dx2**5+dx1*(7150.0d0*dx2**4+dx1* &
(4004.0d0*dx2**3+dx1*(1092.0d0*dx2**2+5.0d0*dx1*(26.0d0*dx2+dx1))))))))
      ci(18)=33.0d0*(52.0d0*dx2**7+dx1*(572.0d0*dx2**6+dx1*(2145.0d0*dx2**5+dx1*(3575.0d0*dx2**4+2.0d0*dx1*(1430.0d0*dx2**3+dx1* &
(546.0d0*dx2**2+dx1*(91.0d0*dx2+5.0d0*dx1)))))))
      ci(19)=11.0d0*(156.0d0*dx2**6+dx1*(1287.0d0*dx2**5+dx1*(3575.0d0*dx2**4+6.0d0*dx1*(715.0d0*dx2**3+dx1*(390.0d0*dx2**2 &
+7.0d0*dx1*(13.0d0*dx2+dx1))))))
      ci(20)=11.0d0*(117.0d0*dx2**5+dx1*(715.0d0*dx2**4+2.0d0*dx1*(715.0d0*dx2**3+3.0d0*dx1*(195.0d0*dx2**2+dx1*(65.0d0*dx2 &
+7.0d0*dx1)))))
      ci(21)=11.0d0*(65.0d0*dx2**4+dx1*(286.0d0*dx2**3+15.0d0*dx1*(26.0d0*dx2**2+dx1*(13.0d0*dx2+2.0d0*dx1))))
      ci(22)=11.0d0*(26.0d0*dx2**3+dx1*(78.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+3.0d0*dx1)))
      ci(23)=78.0d0*dx2**2+11.0d0*dx1*(13.0d0*dx2+5.0d0*dx1)
      ci(24)=13.0d0*dx2+11.0d0*dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>13, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^12 * (x-x2)^n2 as sum_k=0^(12+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_12(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**12
      ci(2)=12.0d0*dx1**11
      ci(3)=66.0d0*dx1**10
      ci(4)=220.0d0*dx1**9
      ci(5)=495.0d0*dx1**8
      ci(6)=792.0d0*dx1**7
      ci(7)=924.0d0*dx1**6
      ci(8)=792.0d0*dx1**5
      ci(9)=495.0d0*dx1**4
      ci(10)=220.0d0*dx1**3
      ci(11)=66.0d0*dx1**2
      ci(12)=12.0d0*dx1
      ci(13)=1.0d0
    case(1)
      ci(1)=dx1**12*dx2
      ci(2)=dx1**11*(12.0d0*dx2+dx1)
      ci(3)=6.0d0*dx1**10*(11.0d0*dx2+2.0d0*dx1)
      ci(4)=22.0d0*dx1**9*(10.0d0*dx2+3.0d0*dx1)
      ci(5)=55.0d0*dx1**8*(9.0d0*dx2+4.0d0*dx1)
      ci(6)=99.0d0*dx1**7*(8.0d0*dx2+5.0d0*dx1)
      ci(7)=132.0d0*dx1**6*(7.0d0*dx2+6.0d0*dx1)
      ci(8)=132.0d0*dx1**5*(6.0d0*dx2+7.0d0*dx1)
      ci(9)=99.0d0*dx1**4*(5.0d0*dx2+8.0d0*dx1)
      ci(10)=55.0d0*dx1**3*(4.0d0*dx2+9.0d0*dx1)
      ci(11)=22.0d0*dx1**2*(3.0d0*dx2+10.0d0*dx1)
      ci(12)=6.0d0*dx1*(2.0d0*dx2+11.0d0*dx1)
      ci(13)=dx2+12.0d0*dx1
      ci(14)=1.0d0
    case(2)
      ci(1)=dx1**12*dx2**2
      ci(2)=2.0d0*dx1**11*dx2*(6.0d0*dx2+dx1)
      ci(3)=dx1**10*(66.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))
      ci(4)=4.0d0*dx1**9*(55.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+dx1))
      ci(5)=11.0d0*dx1**8*(45.0d0*dx2**2+2.0d0*dx1*(20.0d0*dx2+3.0d0*dx1))
      ci(6)=22.0d0*dx1**7*(36.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+2.0d0*dx1))
      ci(7)=33.0d0*dx1**6*(28.0d0*dx2**2+3.0d0*dx1*(16.0d0*dx2+5.0d0*dx1))
      ci(8)=264.0d0*dx1**5*(3.0d0*dx2**2+dx1*(7.0d0*dx2+3.0d0*dx1))
      ci(9)=33.0d0*dx1**4*(15.0d0*dx2**2+4.0d0*dx1*(12.0d0*dx2+7.0d0*dx1))
      ci(10)=22.0d0*dx1**3*(10.0d0*dx2**2+9.0d0*dx1*(5.0d0*dx2+4.0d0*dx1))
      ci(11)=11.0d0*dx1**2*(6.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+9.0d0*dx1))
      ci(12)=4.0d0*dx1*(3.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+5.0d0*dx1))
      ci(13)=dx2**2+6.0d0*dx1*(4.0d0*dx2+11.0d0*dx1)
      ci(14)=2.0d0*(dx2+6.0d0*dx1)
      ci(15)=1.0d0
    case(3)
      ci(1)=dx1**12*dx2**3
      ci(2)=3.0d0*dx1**11*dx2**2*(4.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**10*dx2*(22.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))
      ci(4)=dx1**9*(220.0d0*dx2**3+dx1*(198.0d0*dx2**2+dx1*(36.0d0*dx2+dx1)))
      ci(5)=3.0d0*dx1**8*(165.0d0*dx2**3+2.0d0*dx1*(110.0d0*dx2**2+dx1*(33.0d0*dx2+2.0d0*dx1)))
      ci(6)=33.0d0*dx1**7*(24.0d0*dx2**3+dx1*(45.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+dx1)))
      ci(7)=11.0d0*dx1**6*(84.0d0*dx2**3+dx1*(216.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+4.0d0*dx1)))
      ci(8)=99.0d0*dx1**5*(8.0d0*dx2**3+dx1*(28.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1)))
      ci(9)=99.0d0*dx1**4*(5.0d0*dx2**3+4.0d0*dx1*(6.0d0*dx2**2+dx1*(7.0d0*dx2+2.0d0*dx1)))
      ci(10)=11.0d0*dx1**3*(20.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+4.0d0*dx1*(18.0d0*dx2+7.0d0*dx1)))
      ci(11)=33.0d0*dx1**2*(2.0d0*dx2**3+dx1*(20.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+8.0d0*dx1)))
      ci(12)=3.0d0*dx1*(4.0d0*dx2**3+11.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+3.0d0*dx1)))
      ci(13)=dx2**3+2.0d0*dx1*(18.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+10.0d0*dx1))
      ci(14)=3.0d0*(dx2**2+2.0d0*dx1*(6.0d0*dx2+11.0d0*dx1))
      ci(15)=3.0d0*(dx2+4.0d0*dx1)
      ci(16)=1.0d0
    case(4)
      ci(1)=dx1**12*dx2**4
      ci(2)=4.0d0*dx1**11*dx2**3*(3.0d0*dx2+dx1)
      ci(3)=6.0d0*dx1**10*dx2**2*(11.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(4)=4.0d0*dx1**9*dx2*(55.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))
      ci(5)=dx1**8*(495.0d0*dx2**4+dx1*(880.0d0*dx2**3+dx1*(396.0d0*dx2**2+dx1*(48.0d0*dx2+dx1))))
      ci(6)=12.0d0*dx1**7*(66.0d0*dx2**4+dx1*(165.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))))
      ci(7)=22.0d0*dx1**6*(42.0d0*dx2**4+dx1*(144.0d0*dx2**3+dx1*(135.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))))
      ci(8)=44.0d0*dx1**5*(18.0d0*dx2**4+dx1*(84.0d0*dx2**3+dx1*(108.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+dx1))))
      ci(9)=99.0d0*dx1**4*(5.0d0*dx2**4+dx1*(32.0d0*dx2**3+dx1*(56.0d0*dx2**2+dx1*(32.0d0*dx2+5.0d0*dx1))))
      ci(10)=44.0d0*dx1**3*(5.0d0*dx2**4+3.0d0*dx1*(15.0d0*dx2**3+2.0d0*dx1*(18.0d0*dx2**2+dx1*(14.0d0*dx2+3.0d0*dx1))))
      ci(11)=22.0d0*dx1**2*(3.0d0*dx2**4+dx1*(40.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+2.0d0*dx1*(24.0d0*dx2+7.0d0*dx1))))
      ci(12)=12.0d0*dx1*(dx2**4+11.0d0*dx1*(2.0d0*dx2**3+dx1*(10.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))))
      ci(13)=dx2**4+dx1*(48.0d0*dx2**3+11.0d0*dx1*(36.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+9.0d0*dx1)))
      ci(14)=4.0d0*(dx2**3+dx1*(18.0d0*dx2**2+11.0d0*dx1*(6.0d0*dx2+5.0d0*dx1)))
      ci(15)=6.0d0*(dx2**2+dx1*(8.0d0*dx2+11.0d0*dx1))
      ci(16)=4.0d0*(dx2+3.0d0*dx1)
      ci(17)=1.0d0
    case(5)
      ci(1)=dx1**12*dx2**5
      ci(2)=dx1**11*dx2**4*(12.0d0*dx2+5.0d0*dx1)
      ci(3)=2.0d0*dx1**10*dx2**3*(33.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+dx1))
      ci(4)=10.0d0*dx1**9*dx2**2*(22.0d0*dx2**3+dx1*(33.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**8*dx2*(99.0d0*dx2**4+dx1*(220.0d0*dx2**3+dx1*(132.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(6)=dx1**7*(792.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(2200.0d0*dx2**3+dx1*(660.0d0*dx2**2+dx1*(60.0d0*dx2+dx1)))))
      ci(7)=2.0d0*dx1**6*(462.0d0*dx2**5+dx1*(1980.0d0*dx2**4+dx1*(2475.0d0*dx2**3+dx1*(1100.0d0*dx2**2+3.0d0*dx1*(55.0d0*dx2 &
+2.0d0*dx1)))))
      ci(8)=22.0d0*dx1**5*(36.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1*(360.0d0*dx2**3+dx1*(225.0d0*dx2**2+dx1*(50.0d0*dx2+3.0d0*dx1)))) &
)
      ci(9)=55.0d0*dx1**4*(9.0d0*dx2**5+dx1*(72.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(144.0d0*dx2**2+dx1*(45.0d0*dx2+4.0d0*dx1)))))
      ci(10)=55.0d0*dx1**3*(4.0d0*dx2**5+3.0d0*dx1*(15.0d0*dx2**4+dx1*(48.0d0*dx2**3+dx1*(56.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+dx1)) &
)))
      ci(11)=22.0d0*dx1**2*(3.0d0*dx2**5+dx1*(50.0d0*dx2**4+3.0d0*dx1*(75.0d0*dx2**3+2.0d0*dx1*(60.0d0*dx2**2+dx1*(35.0d0*dx2 &
+6.0d0*dx1)))))
      ci(12)=2.0d0*dx1*(6.0d0*dx2**5+11.0d0*dx1*(15.0d0*dx2**4+dx1*(100.0d0*dx2**3+3.0d0*dx1*(75.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2 &
+7.0d0*dx1)))))
      ci(13)=dx2**5+dx1*(60.0d0*dx2**4+11.0d0*dx1*(60.0d0*dx2**3+dx1*(200.0d0*dx2**2+9.0d0*dx1*(25.0d0*dx2+8.0d0*dx1))))
      ci(14)=5.0d0*(dx2**4+dx1*(24.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+dx1*(20.0d0*dx2+9.0d0*dx1))))
      ci(15)=10.0d0*(dx2**3+dx1*(12.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(16)=2.0d0*(5.0d0*dx2**2+3.0d0*dx1*(10.0d0*dx2+11.0d0*dx1))
      ci(17)=5.0d0*dx2+12.0d0*dx1
      ci(18)=1.0d0
    case(6)
      ci(1)=dx1**12*dx2**6
      ci(2)=6.0d0*dx1**11*dx2**5*(2.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**10*dx2**4*(22.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1))
      ci(4)=4.0d0*dx1**9*dx2**3*(55.0d0*dx2**3+dx1*(99.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+dx1)))
      ci(5)=15.0d0*dx1**8*dx2**2*(33.0d0*dx2**4+dx1*(88.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(6)=6.0d0*dx1**7*dx2*(132.0d0*dx2**5+dx1*(495.0d0*dx2**4+dx1*(550.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))
      ci(7)=dx1**6*(924.0d0*dx2**6+dx1*(4752.0d0*dx2**5+dx1*(7425.0d0*dx2**4+dx1*(4400.0d0*dx2**3+dx1*(990.0d0*dx2**2+dx1* &
(72.0d0*dx2+dx1))))))
      ci(8)=12.0d0*dx1**5*(66.0d0*dx2**6+dx1*(462.0d0*dx2**5+dx1*(990.0d0*dx2**4+dx1*(825.0d0*dx2**3+dx1*(275.0d0*dx2**2+dx1* &
(33.0d0*dx2+dx1))))))
      ci(9)=33.0d0*dx1**4*(15.0d0*dx2**6+dx1*(144.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1*(480.0d0*dx2**3+dx1*(225.0d0*dx2**2 &
+2.0d0*dx1*(20.0d0*dx2+dx1))))))
      ci(10)=110.0d0*dx1**3*(2.0d0*dx2**6+dx1*(27.0d0*dx2**5+dx1*(108.0d0*dx2**4+dx1*(168.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1* &
(27.0d0*dx2+2.0d0*dx1))))))
      ci(11)=33.0d0*dx1**2*(2.0d0*dx2**6+dx1*(40.0d0*dx2**5+3.0d0*dx1*(75.0d0*dx2**4+dx1*(160.0d0*dx2**3+dx1*(140.0d0*dx2**2+dx1* &
(48.0d0*dx2+5.0d0*dx1))))))
      ci(12)=12.0d0*dx1*(dx2**6+11.0d0*dx1*(3.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(25.0d0*dx2**3+2.0d0*dx1*(15.0d0*dx2**2+dx1* &
(7.0d0*dx2+dx1))))))
      ci(13)=dx2**6+dx1*(72.0d0*dx2**5+11.0d0*dx1*(90.0d0*dx2**4+dx1*(400.0d0*dx2**3+3.0d0*dx1*(225.0d0*dx2**2+4.0d0*dx1* &
(36.0d0*dx2+7.0d0*dx1)))))
      ci(14)=6.0d0*(dx2**5+dx1*(30.0d0*dx2**4+11.0d0*dx1*(20.0d0*dx2**3+dx1*(50.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+4.0d0*dx1)))))
      ci(15)=15.0d0*(dx2**4+dx1*(16.0d0*dx2**3+11.0d0*dx1*(6.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(16)=4.0d0*(5.0d0*dx2**3+dx1*(45.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+5.0d0*dx1)))
      ci(17)=3.0d0*(5.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+11.0d0*dx1))
      ci(18)=6.0d0*(dx2+2.0d0*dx1)
      ci(19)=1.0d0
    case(7)
      ci(1)=dx1**12*dx2**7
      ci(2)=dx1**11*dx2**6*(12.0d0*dx2+7.0d0*dx1)
      ci(3)=3.0d0*dx1**10*dx2**5*(22.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+dx1))
      ci(4)=dx1**9*dx2**4*(220.0d0*dx2**3+7.0d0*dx1*(66.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1)))
      ci(5)=dx1**8*dx2**3*(495.0d0*dx2**4+7.0d0*dx1*(220.0d0*dx2**3+dx1*(198.0d0*dx2**2+5.0d0*dx1*(12.0d0*dx2+dx1))))
      ci(6)=3.0d0*dx1**7*dx2**2*(264.0d0*dx2**5+7.0d0*dx1*(165.0d0*dx2**4+dx1*(220.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(20.0d0*dx2 &
+dx1)))))
      ci(7)=7.0d0*dx1**6*dx2*(132.0d0*dx2**6+dx1*(792.0d0*dx2**5+dx1*(1485.0d0*dx2**4+dx1*(1100.0d0*dx2**3+dx1*(330.0d0*dx2**2 &
+dx1*(36.0d0*dx2+dx1))))))
      ci(8)=dx1**5*(792.0d0*dx2**7+dx1*(6468.0d0*dx2**6+dx1*(16632.0d0*dx2**5+dx1*(17325.0d0*dx2**4+dx1*(7700.0d0*dx2**3+dx1* &
(1386.0d0*dx2**2+dx1*(84.0d0*dx2+dx1)))))))
      ci(9)=3.0d0*dx1**4*(165.0d0*dx2**7+dx1*(1848.0d0*dx2**6+dx1*(6468.0d0*dx2**5+dx1*(9240.0d0*dx2**4+dx1*(5775.0d0*dx2**3 &
+2.0d0*dx1*(770.0d0*dx2**2+dx1*(77.0d0*dx2+2.0d0*dx1)))))))
      ci(10)=11.0d0*dx1**3*(20.0d0*dx2**7+dx1*(315.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1*(2520.0d0*dx2**3+dx1* &
(945.0d0*dx2**2+2.0d0*dx1*(70.0d0*dx2+3.0d0*dx1)))))))
      ci(11)=11.0d0*dx1**2*(6.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(945.0d0*dx2**5+dx1*(2520.0d0*dx2**4+dx1*(2940.0d0*dx2**3+dx1* &
(1512.0d0*dx2**2+5.0d0*dx1*(63.0d0*dx2+4.0d0*dx1)))))))
      ci(12)=3.0d0*dx1*(4.0d0*dx2**7+11.0d0*dx1*(14.0d0*dx2**6+dx1*(140.0d0*dx2**5+3.0d0*dx1*(175.0d0*dx2**4+dx1*(280.0d0*dx2**3 &
+dx1*(196.0d0*dx2**2+dx1*(56.0d0*dx2+5.0d0*dx1)))))))
      ci(13)=dx2**7+dx1*(84.0d0*dx2**6+11.0d0*dx1*(126.0d0*dx2**5+dx1*(700.0d0*dx2**4+3.0d0*dx1*(525.0d0*dx2**3+4.0d0*dx1* &
(126.0d0*dx2**2+dx1*(49.0d0*dx2+6.0d0*dx1))))))
      ci(14)=7.0d0*(dx2**6+dx1*(36.0d0*dx2**5+11.0d0*dx1*(30.0d0*dx2**4+dx1*(100.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+4.0d0*dx1* &
(6.0d0*dx2+dx1))))))
      ci(15)=3.0d0*(7.0d0*dx2**5+dx1*(140.0d0*dx2**4+11.0d0*dx1*(70.0d0*dx2**3+dx1*(140.0d0*dx2**2+3.0d0*dx1*(35.0d0*dx2+8.0d0*dx1 &
)))))
      ci(16)=35.0d0*dx2**4+dx1*(420.0d0*dx2**3+11.0d0*dx1*(126.0d0*dx2**2+5.0d0*dx1*(28.0d0*dx2+9.0d0*dx1)))
      ci(17)=35.0d0*dx2**3+2.0d0*dx1*(126.0d0*dx2**2+11.0d0*dx1*(21.0d0*dx2+10.0d0*dx1))
      ci(18)=3.0d0*(7.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+11.0d0*dx1))
      ci(19)=7.0d0*dx2+12.0d0*dx1
      ci(20)=1.0d0
    case(8)
      ci(1)=dx1**12*dx2**8
      ci(2)=4.0d0*dx1**11*dx2**7*(3.0d0*dx2+2.0d0*dx1)
      ci(3)=2.0d0*dx1**10*dx2**6*(33.0d0*dx2**2+2.0d0*dx1*(24.0d0*dx2+7.0d0*dx1))
      ci(4)=4.0d0*dx1**9*dx2**5*(55.0d0*dx2**3+2.0d0*dx1*(66.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(5)=dx1**8*dx2**4*(495.0d0*dx2**4+2.0d0*dx1*(880.0d0*dx2**3+7.0d0*dx1*(132.0d0*dx2**2+dx1*(48.0d0*dx2+5.0d0*dx1))))
      ci(6)=8.0d0*dx1**7*dx2**3*(99.0d0*dx2**5+dx1*(495.0d0*dx2**4+7.0d0*dx1*(110.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(15.0d0*dx2 &
+dx1)))))
      ci(7)=4.0d0*dx1**6*dx2**2*(231.0d0*dx2**6+dx1*(1584.0d0*dx2**5+7.0d0*dx1*(495.0d0*dx2**4+dx1*(440.0d0*dx2**3+dx1* &
(165.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))))
      ci(8)=8.0d0*dx1**5*dx2*(99.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(3465.0d0*dx2**4+dx1*(1925.0d0*dx2**3 &
+dx1*(462.0d0*dx2**2+dx1*(42.0d0*dx2+dx1)))))))
      ci(9)=dx1**4*(495.0d0*dx2**8+dx1*(6336.0d0*dx2**7+dx1*(25872.0d0*dx2**6+dx1*(44352.0d0*dx2**5+dx1*(34650.0d0*dx2**4+dx1* &
(12320.0d0*dx2**3+dx1*(1848.0d0*dx2**2+dx1*(96.0d0*dx2+dx1))))))))
      ci(10)=4.0d0*dx1**3*(55.0d0*dx2**8+dx1*(990.0d0*dx2**7+dx1*(5544.0d0*dx2**6+dx1*(12936.0d0*dx2**5+dx1*(13860.0d0*dx2**4+dx1* &
(6930.0d0*dx2**3+dx1*(1540.0d0*dx2**2+3.0d0*dx1*(44.0d0*dx2+dx1))))))))
      ci(11)=22.0d0*dx1**2*(3.0d0*dx2**8+dx1*(80.0d0*dx2**7+dx1*(630.0d0*dx2**6+dx1*(2016.0d0*dx2**5+dx1*(2940.0d0*dx2**4+dx1* &
(2016.0d0*dx2**3+dx1*(630.0d0*dx2**2+dx1*(80.0d0*dx2+3.0d0*dx1))))))))
      ci(12)=4.0d0*dx1*(3.0d0*dx2**8+11.0d0*dx1*(12.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(630.0d0*dx2**5+dx1*(1260.0d0*dx2**4+dx1* &
(1176.0d0*dx2**3+dx1*(504.0d0*dx2**2+5.0d0*dx1*(18.0d0*dx2+dx1))))))))
      ci(13)=dx2**8+dx1*(96.0d0*dx2**7+11.0d0*dx1*(168.0d0*dx2**6+dx1*(1120.0d0*dx2**5+3.0d0*dx1*(1050.0d0*dx2**4+dx1* &
(1344.0d0*dx2**3+dx1*(784.0d0*dx2**2+3.0d0*dx1*(64.0d0*dx2+5.0d0*dx1)))))))
      ci(14)=8.0d0*(dx2**7+dx1*(42.0d0*dx2**6+11.0d0*dx1*(42.0d0*dx2**5+dx1*(175.0d0*dx2**4+3.0d0*dx1*(105.0d0*dx2**3+dx1* &
(84.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1)))))))
      ci(15)=4.0d0*(7.0d0*dx2**6+dx1*(168.0d0*dx2**5+11.0d0*dx1*(105.0d0*dx2**4+dx1*(280.0d0*dx2**3+3.0d0*dx1*(105.0d0*dx2**2+dx1* &
(48.0d0*dx2+7.0d0*dx1))))))
      ci(16)=8.0d0*(7.0d0*dx2**5+dx1*(105.0d0*dx2**4+11.0d0*dx1*(42.0d0*dx2**3+dx1*(70.0d0*dx2**2+9.0d0*dx1*(5.0d0*dx2+dx1)))))
      ci(17)=70.0d0*dx2**4+dx1*(672.0d0*dx2**3+11.0d0*dx1*(168.0d0*dx2**2+5.0d0*dx1*(32.0d0*dx2+9.0d0*dx1)))
      ci(18)=4.0d0*(14.0d0*dx2**3+dx1*(84.0d0*dx2**2+11.0d0*dx1*(12.0d0*dx2+5.0d0*dx1)))
      ci(19)=2.0d0*(14.0d0*dx2**2+3.0d0*dx1*(16.0d0*dx2+11.0d0*dx1))
      ci(20)=4.0d0*(2.0d0*dx2+3.0d0*dx1)
      ci(21)=1.0d0
    case(9)
      ci(1)=dx1**12*dx2**9
      ci(2)=3.0d0*dx1**11*dx2**8*(4.0d0*dx2+3.0d0*dx1)
      ci(3)=6.0d0*dx1**10*dx2**7*(11.0d0*dx2**2+6.0d0*dx1*(3.0d0*dx2+dx1))
      ci(4)=2.0d0*dx1**9*dx2**6*(110.0d0*dx2**3+3.0d0*dx1*(99.0d0*dx2**2+2.0d0*dx1*(36.0d0*dx2+7.0d0*dx1)))
      ci(5)=9.0d0*dx1**8*dx2**5*(55.0d0*dx2**4+2.0d0*dx1*(110.0d0*dx2**3+dx1*(132.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(6)=9.0d0*dx1**7*dx2**4*(88.0d0*dx2**5+dx1*(495.0d0*dx2**4+2.0d0*dx1*(440.0d0*dx2**3+7.0d0*dx1*(44.0d0*dx2**2+dx1* &
(12.0d0*dx2+dx1)))))
      ci(7)=12.0d0*dx1**6*dx2**3*(77.0d0*dx2**6+dx1*(594.0d0*dx2**5+dx1*(1485.0d0*dx2**4+7.0d0*dx1*(220.0d0*dx2**3+dx1* &
(99.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))))))
      ci(8)=36.0d0*dx1**5*dx2**2*(22.0d0*dx2**7+dx1*(231.0d0*dx2**6+dx1*(792.0d0*dx2**5+dx1*(1155.0d0*dx2**4+dx1*(770.0d0*dx2**3 &
+dx1*(231.0d0*dx2**2+dx1*(28.0d0*dx2+dx1)))))))
      ci(9)=9.0d0*dx1**4*dx2*(55.0d0*dx2**8+dx1*(792.0d0*dx2**7+dx1*(3696.0d0*dx2**6+dx1*(7392.0d0*dx2**5+dx1*(6930.0d0*dx2**4 &
+dx1*(3080.0d0*dx2**3+dx1*(616.0d0*dx2**2+dx1*(48.0d0*dx2+dx1))))))))
      ci(10)=dx1**3*(220.0d0*dx2**9+dx1*(4455.0d0*dx2**8+dx1*(28512.0d0*dx2**7+dx1*(77616.0d0*dx2**6+dx1*(99792.0d0*dx2**5+dx1* &
(62370.0d0*dx2**4+dx1*(18480.0d0*dx2**3+dx1*(2376.0d0*dx2**2+dx1*(108.0d0*dx2+dx1)))))))))
      ci(11)=6.0d0*dx1**2*(11.0d0*dx2**9+dx1*(330.0d0*dx2**8+dx1*(2970.0d0*dx2**7+dx1*(11088.0d0*dx2**6+dx1*(19404.0d0*dx2**5+dx1* &
(16632.0d0*dx2**4+dx1*(6930.0d0*dx2**3+dx1*(1320.0d0*dx2**2+dx1*(99.0d0*dx2+2.0d0*dx1)))))))))
      ci(12)=6.0d0*dx1*(2.0d0*dx2**9+11.0d0*dx1*(9.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(630.0d0*dx2**6+dx1*(1512.0d0*dx2**5+dx1* &
(1764.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(270.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))))))
      ci(13)=dx2**9+dx1*(108.0d0*dx2**8+11.0d0*dx1*(216.0d0*dx2**7+dx1*(1680.0d0*dx2**6+dx1*(5670.0d0*dx2**5+dx1*(9072.0d0*dx2**4 &
+dx1*(7056.0d0*dx2**3+dx1*(2592.0d0*dx2**2+5.0d0*dx1*(81.0d0*dx2+4.0d0*dx1))))))))
      ci(14)=9.0d0*(dx2**8+dx1*(48.0d0*dx2**7+11.0d0*dx1*(56.0d0*dx2**6+dx1*(280.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1* &
(672.0d0*dx2**3+dx1*(336.0d0*dx2**2+dx1*(72.0d0*dx2+5.0d0*dx1))))))))
      ci(15)=36.0d0*(dx2**7+dx1*(28.0d0*dx2**6+11.0d0*dx1*(21.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(105.0d0*dx2**3+dx1*(72.0d0*dx2**2 &
+dx1*(21.0d0*dx2+2.0d0*dx1)))))))
      ci(16)=12.0d0*(7.0d0*dx2**6+dx1*(126.0d0*dx2**5+11.0d0*dx1*(63.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(135.0d0*dx2**2+dx1* &
(54.0d0*dx2+7.0d0*dx1))))))
      ci(17)=9.0d0*(14.0d0*dx2**5+dx1*(168.0d0*dx2**4+11.0d0*dx1*(56.0d0*dx2**3+dx1*(80.0d0*dx2**2+dx1*(45.0d0*dx2+8.0d0*dx1)))))
      ci(18)=9.0d0*(14.0d0*dx2**4+dx1*(112.0d0*dx2**3+11.0d0*dx1*(24.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(19)=2.0d0*(42.0d0*dx2**3+dx1*(216.0d0*dx2**2+11.0d0*dx1*(27.0d0*dx2+10.0d0*dx1)))
      ci(20)=6.0d0*(6.0d0*dx2**2+dx1*(18.0d0*dx2+11.0d0*dx1))
      ci(21)=3.0d0*(3.0d0*dx2+4.0d0*dx1)
      ci(22)=1.0d0
    case(10)
      ci(1)=dx1**12*dx2**10
      ci(2)=2.0d0*dx1**11*dx2**9*(6.0d0*dx2+5.0d0*dx1)
      ci(3)=3.0d0*dx1**10*dx2**8*(22.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(4)=20.0d0*dx1**9*dx2**7*(11.0d0*dx2**3+3.0d0*dx1*(11.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(5)=5.0d0*dx1**8*dx2**6*(99.0d0*dx2**4+2.0d0*dx1*(220.0d0*dx2**3+3.0d0*dx1*(99.0d0*dx2**2+dx1*(48.0d0*dx2+7.0d0*dx1))))
      ci(6)=18.0d0*dx1**7*dx2**5*(44.0d0*dx2**5+dx1*(275.0d0*dx2**4+2.0d0*dx1*(275.0d0*dx2**3+dx1*(220.0d0*dx2**2+7.0d0*dx1* &
(10.0d0*dx2+dx1)))))
      ci(7)=3.0d0*dx1**6*dx2**4*(308.0d0*dx2**6+dx1*(2640.0d0*dx2**5+dx1*(7425.0d0*dx2**4+2.0d0*dx1*(4400.0d0*dx2**3+7.0d0*dx1* &
(330.0d0*dx2**2+dx1*(72.0d0*dx2+5.0d0*dx1))))))
      ci(8)=24.0d0*dx1**5*dx2**3*(33.0d0*dx2**7+dx1*(385.0d0*dx2**6+dx1*(1485.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(1925.0d0*dx2**3 &
+dx1*(693.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+dx1)))))))
      ci(9)=45.0d0*dx1**4*dx2**2*(11.0d0*dx2**8+dx1*(176.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2112.0d0*dx2**5+dx1*(2310.0d0*dx2**4 &
+dx1*(1232.0d0*dx2**3+dx1*(308.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))))))
      ci(10)=10.0d0*dx1**3*dx2*(22.0d0*dx2**9+dx1*(495.0d0*dx2**8+dx1*(3564.0d0*dx2**7+dx1*(11088.0d0*dx2**6+dx1*(16632.0d0*dx2**5 &
+dx1*(12474.0d0*dx2**4+dx1*(4620.0d0*dx2**3+dx1*(792.0d0*dx2**2+dx1*(54.0d0*dx2+dx1)))))))))
      ci(11)=dx1**2*(66.0d0*dx2**10+dx1*(2200.0d0*dx2**9+dx1*(22275.0d0*dx2**8+dx1*(95040.0d0*dx2**7+dx1*(194040.0d0*dx2**6+dx1* &
(199584.0d0*dx2**5+dx1*(103950.0d0*dx2**4+dx1*(26400.0d0*dx2**3+dx1*(2970.0d0*dx2**2+dx1*(120.0d0*dx2+dx1))))))))))
      ci(12)=12.0d0*dx1*(dx2**10+dx1*(55.0d0*dx2**9+dx1*(825.0d0*dx2**8+dx1*(4950.0d0*dx2**7+dx1*(13860.0d0*dx2**6+dx1* &
(19404.0d0*dx2**5+dx1*(13860.0d0*dx2**4+dx1*(4950.0d0*dx2**3+dx1*(825.0d0*dx2**2+dx1*(55.0d0*dx2+dx1))))))))))
      ci(13)=dx2**10+dx1*(120.0d0*dx2**9+11.0d0*dx1*(270.0d0*dx2**8+dx1*(2400.0d0*dx2**7+dx1*(9450.0d0*dx2**6+dx1* &
(18144.0d0*dx2**5+dx1*(17640.0d0*dx2**4+dx1*(8640.0d0*dx2**3+dx1*(2025.0d0*dx2**2+2.0d0*dx1*(100.0d0*dx2+3.0d0*dx1)))))))))
      ci(14)=10.0d0*(dx2**9+dx1*(54.0d0*dx2**8+11.0d0*dx1*(72.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(1134.0d0*dx2**5+dx1* &
(1512.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(324.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))))))))
      ci(15)=45.0d0*(dx2**8+dx1*(32.0d0*dx2**7+11.0d0*dx1*(28.0d0*dx2**6+dx1*(112.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1* &
(192.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))))))
      ci(16)=24.0d0*(5.0d0*dx2**7+dx1*(105.0d0*dx2**6+11.0d0*dx1*(63.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(225.0d0*dx2**3+dx1* &
(135.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1)))))))
      ci(17)=3.0d0*(70.0d0*dx2**6+dx1*(1008.0d0*dx2**5+11.0d0*dx1*(420.0d0*dx2**4+dx1*(800.0d0*dx2**3+dx1*(675.0d0*dx2**2 &
+4.0d0*dx1*(60.0d0*dx2+7.0d0*dx1))))))
      ci(18)=18.0d0*(14.0d0*dx2**5+dx1*(140.0d0*dx2**4+11.0d0*dx1*(40.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(25.0d0*dx2+4.0d0*dx1)))))
      ci(19)=5.0d0*(42.0d0*dx2**4+dx1*(288.0d0*dx2**3+11.0d0*dx1*(54.0d0*dx2**2+dx1*(40.0d0*dx2+9.0d0*dx1))))
      ci(20)=20.0d0*(6.0d0*dx2**3+dx1*(27.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(21)=3.0d0*(15.0d0*dx2**2+2.0d0*dx1*(20.0d0*dx2+11.0d0*dx1))
      ci(22)=2.0d0*(5.0d0*dx2+6.0d0*dx1)
      ci(23)=1.0d0
    case(11)
      ci(1)=dx1**12*dx2**11
      ci(2)=dx1**11*dx2**10*(12.0d0*dx2+11.0d0*dx1)
      ci(3)=11.0d0*dx1**10*dx2**9*(6.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(4)=11.0d0*dx1**9*dx2**8*(20.0d0*dx2**3+3.0d0*dx1*(22.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1)))
      ci(5)=55.0d0*dx1**8*dx2**7*(9.0d0*dx2**4+2.0d0*dx1*(22.0d0*dx2**3+3.0d0*dx1*(11.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))
      ci(6)=11.0d0*dx1**7*dx2**6*(72.0d0*dx2**5+dx1*(495.0d0*dx2**4+2.0d0*dx1*(550.0d0*dx2**3+3.0d0*dx1*(165.0d0*dx2**2+dx1* &
(60.0d0*dx2+7.0d0*dx1)))))
      ci(7)=33.0d0*dx1**6*dx2**5*(28.0d0*dx2**6+dx1*(264.0d0*dx2**5+dx1*(825.0d0*dx2**4+2.0d0*dx1*(550.0d0*dx2**3+dx1* &
(330.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+dx1))))))
      ci(8)=33.0d0*dx1**5*dx2**4*(24.0d0*dx2**7+dx1*(308.0d0*dx2**6+dx1*(1320.0d0*dx2**5+dx1*(2475.0d0*dx2**4+2.0d0*dx1* &
(1100.0d0*dx2**3+dx1*(462.0d0*dx2**2+dx1*(84.0d0*dx2+5.0d0*dx1)))))))
      ci(9)=33.0d0*dx1**4*dx2**3*(15.0d0*dx2**8+dx1*(264.0d0*dx2**7+dx1*(1540.0d0*dx2**6+dx1*(3960.0d0*dx2**5+dx1*(4950.0d0*dx2**4 &
+dx1*(3080.0d0*dx2**3+dx1*(924.0d0*dx2**2+5.0d0*dx1*(24.0d0*dx2+dx1))))))))
      ci(10)=55.0d0*dx1**3*dx2**2*(4.0d0*dx2**9+dx1*(99.0d0*dx2**8+dx1*(792.0d0*dx2**7+dx1*(2772.0d0*dx2**6+dx1*(4752.0d0*dx2**5 &
+dx1*(4158.0d0*dx2**4+dx1*(1848.0d0*dx2**3+dx1*(396.0d0*dx2**2+dx1*(36.0d0*dx2+dx1)))))))))
      ci(11)=11.0d0*dx1**2*dx2*(6.0d0*dx2**10+dx1*(220.0d0*dx2**9+dx1*(2475.0d0*dx2**8+dx1*(11880.0d0*dx2**7+dx1*(27720.0d0*dx2**6 &
+dx1*(33264.0d0*dx2**5+dx1*(20790.0d0*dx2**4+dx1*(6600.0d0*dx2**3+dx1*(990.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))))))))
      ci(12)=dx1*(12.0d0*dx2**11+dx1*(726.0d0*dx2**10+dx1*(12100.0d0*dx2**9+dx1*(81675.0d0*dx2**8+dx1*(261360.0d0*dx2**7+dx1* &
(426888.0d0*dx2**6+dx1*(365904.0d0*dx2**5+dx1*(163350.0d0*dx2**4+dx1*(36300.0d0*dx2**3+dx1*(3630.0d0*dx2**2+dx1*(132.0d0*dx2+dx1)) &
)))))))))
      ci(13)=dx2**11+dx1*(132.0d0*dx2**10+dx1*(3630.0d0*dx2**9+dx1*(36300.0d0*dx2**8+dx1*(163350.0d0*dx2**7+dx1*(365904.0d0*dx2**6 &
+dx1*(426888.0d0*dx2**5+dx1*(261360.0d0*dx2**4+dx1*(81675.0d0*dx2**3+2.0d0*dx1*(6050.0d0*dx2**2+3.0d0*dx1*(121.0d0*dx2+2.0d0*dx1)) &
))))))))
      ci(14)=11.0d0*(dx2**10+dx1*(60.0d0*dx2**9+dx1*(990.0d0*dx2**8+dx1*(6600.0d0*dx2**7+dx1*(20790.0d0*dx2**6+dx1* &
(33264.0d0*dx2**5+dx1*(27720.0d0*dx2**4+dx1*(11880.0d0*dx2**3+dx1*(2475.0d0*dx2**2+2.0d0*dx1*(110.0d0*dx2+3.0d0*dx1))))))))))
      ci(15)=55.0d0*(dx2**9+dx1*(36.0d0*dx2**8+dx1*(396.0d0*dx2**7+dx1*(1848.0d0*dx2**6+dx1*(4158.0d0*dx2**5+dx1*(4752.0d0*dx2**4 &
+dx1*(2772.0d0*dx2**3+dx1*(792.0d0*dx2**2+dx1*(99.0d0*dx2+4.0d0*dx1)))))))))
      ci(16)=33.0d0*(5.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(3080.0d0*dx2**5+dx1*(4950.0d0*dx2**4+dx1* &
(3960.0d0*dx2**3+dx1*(1540.0d0*dx2**2+3.0d0*dx1*(88.0d0*dx2+5.0d0*dx1))))))))
      ci(17)=33.0d0*(10.0d0*dx2**7+dx1*(168.0d0*dx2**6+dx1*(924.0d0*dx2**5+dx1*(2200.0d0*dx2**4+dx1*(2475.0d0*dx2**3+4.0d0*dx1* &
(330.0d0*dx2**2+dx1*(77.0d0*dx2+6.0d0*dx1)))))))
      ci(18)=33.0d0*(14.0d0*dx2**6+dx1*(168.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(1100.0d0*dx2**3+dx1*(825.0d0*dx2**2+4.0d0*dx1* &
(66.0d0*dx2+7.0d0*dx1))))))
      ci(19)=11.0d0*(42.0d0*dx2**5+dx1*(360.0d0*dx2**4+dx1*(990.0d0*dx2**3+dx1*(1100.0d0*dx2**2+9.0d0*dx1*(55.0d0*dx2+8.0d0*dx1))) &
))
      ci(20)=55.0d0*(6.0d0*dx2**4+dx1*(36.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(44.0d0*dx2+9.0d0*dx1))))
      ci(21)=11.0d0*(15.0d0*dx2**3+2.0d0*dx1*(30.0d0*dx2**2+dx1*(33.0d0*dx2+10.0d0*dx1)))
      ci(22)=11.0d0*(5.0d0*dx2**2+6.0d0*dx1*(2.0d0*dx2+dx1))
      ci(23)=11.0d0*dx2+12.0d0*dx1
      ci(24)=1.0d0
    case(12)
      ci(1)=dx1**12*dx2**12
      ci(2)=12.0d0*dx1**11*dx2**11*(dx2+dx1)
      ci(3)=6.0d0*dx1**10*dx2**10*(11.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1))
      ci(4)=44.0d0*dx1**9*dx2**9*(5.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1)))
      ci(5)=33.0d0*dx1**8*dx2**8*(15.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(132.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(6)=132.0d0*dx1**7*dx2**7*(6.0d0*dx2**5+dx1*(45.0d0*dx2**4+dx1*(110.0d0*dx2**3+dx1*(110.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=22.0d0*dx1**6*dx2**6*(42.0d0*dx2**6+dx1*(432.0d0*dx2**5+dx1*(1485.0d0*dx2**4+dx1*(2200.0d0*dx2**3+3.0d0*dx1* &
(495.0d0*dx2**2+2.0d0*dx1*(72.0d0*dx2+7.0d0*dx1))))))
      ci(8)=396.0d0*dx1**5*dx2**5*(2.0d0*dx2**7+dx1*(28.0d0*dx2**6+dx1*(132.0d0*dx2**5+dx1*(275.0d0*dx2**4+dx1*(275.0d0*dx2**3 &
+2.0d0*dx1*(66.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)))))))
      ci(9)=99.0d0*dx1**4*dx2**4*(5.0d0*dx2**8+dx1*(96.0d0*dx2**7+dx1*(616.0d0*dx2**6+dx1*(1760.0d0*dx2**5+dx1*(2475.0d0*dx2**4 &
+dx1*(1760.0d0*dx2**3+dx1*(616.0d0*dx2**2+dx1*(96.0d0*dx2+5.0d0*dx1))))))))
      ci(10)=44.0d0*dx1**3*dx2**3*(5.0d0*dx2**9+dx1*(135.0d0*dx2**8+dx1*(1188.0d0*dx2**7+dx1*(4620.0d0*dx2**6+dx1*(8910.0d0*dx2**5 &
+dx1*(8910.0d0*dx2**4+dx1*(4620.0d0*dx2**3+dx1*(1188.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+dx1)))))))))
      ci(11)=66.0d0*dx1**2*dx2**2*(dx2**10+dx1*(40.0d0*dx2**9+dx1*(495.0d0*dx2**8+dx1*(2640.0d0*dx2**7+dx1*(6930.0d0*dx2**6+dx1* &
(9504.0d0*dx2**5+dx1*(6930.0d0*dx2**4+dx1*(2640.0d0*dx2**3+dx1*(495.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))))))))
      ci(12)=12.0d0*dx1*dx2*(dx2**11+dx1*(66.0d0*dx2**10+dx1*(1210.0d0*dx2**9+dx1*(9075.0d0*dx2**8+dx1*(32670.0d0*dx2**7+dx1* &
(60984.0d0*dx2**6+dx1*(60984.0d0*dx2**5+dx1*(32670.0d0*dx2**4+dx1*(9075.0d0*dx2**3+dx1*(1210.0d0*dx2**2+dx1*(66.0d0*dx2+dx1))))))) &
))))
      ci(13)=dx2**12+dx1*(144.0d0*dx2**11+dx1*(4356.0d0*dx2**10+dx1*(48400.0d0*dx2**9+dx1*(245025.0d0*dx2**8+dx1* &
(627264.0d0*dx2**7+dx1*(853776.0d0*dx2**6+dx1*(627264.0d0*dx2**5+dx1*(245025.0d0*dx2**4+dx1*(48400.0d0*dx2**3+dx1*(4356.0d0*dx2**2 &
+dx1*(144.0d0*dx2+dx1)))))))))))
      ci(14)=12.0d0*(dx2**11+dx1*(66.0d0*dx2**10+dx1*(1210.0d0*dx2**9+dx1*(9075.0d0*dx2**8+dx1*(32670.0d0*dx2**7+dx1* &
(60984.0d0*dx2**6+dx1*(60984.0d0*dx2**5+dx1*(32670.0d0*dx2**4+dx1*(9075.0d0*dx2**3+dx1*(1210.0d0*dx2**2+dx1*(66.0d0*dx2+dx1))))))) &
))))
      ci(15)=66.0d0*(dx2**10+dx1*(40.0d0*dx2**9+dx1*(495.0d0*dx2**8+dx1*(2640.0d0*dx2**7+dx1*(6930.0d0*dx2**6+dx1*(9504.0d0*dx2**5 &
+dx1*(6930.0d0*dx2**4+dx1*(2640.0d0*dx2**3+dx1*(495.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))))))))
      ci(16)=44.0d0*(5.0d0*dx2**9+dx1*(135.0d0*dx2**8+dx1*(1188.0d0*dx2**7+dx1*(4620.0d0*dx2**6+dx1*(8910.0d0*dx2**5+dx1* &
(8910.0d0*dx2**4+dx1*(4620.0d0*dx2**3+dx1*(1188.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+dx1)))))))))
      ci(17)=99.0d0*(5.0d0*dx2**8+dx1*(96.0d0*dx2**7+dx1*(616.0d0*dx2**6+dx1*(1760.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1* &
(1760.0d0*dx2**3+dx1*(616.0d0*dx2**2+dx1*(96.0d0*dx2+5.0d0*dx1))))))))
      ci(18)=396.0d0*(2.0d0*dx2**7+dx1*(28.0d0*dx2**6+dx1*(132.0d0*dx2**5+dx1*(275.0d0*dx2**4+dx1*(275.0d0*dx2**3+2.0d0*dx1* &
(66.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)))))))
      ci(19)=22.0d0*(42.0d0*dx2**6+dx1*(432.0d0*dx2**5+dx1*(1485.0d0*dx2**4+dx1*(2200.0d0*dx2**3+3.0d0*dx1*(495.0d0*dx2**2 &
+2.0d0*dx1*(72.0d0*dx2+7.0d0*dx1))))))
      ci(20)=132.0d0*(6.0d0*dx2**5+dx1*(45.0d0*dx2**4+dx1*(110.0d0*dx2**3+dx1*(110.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+2.0d0*dx1)))))
      ci(21)=33.0d0*(15.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(132.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(22)=44.0d0*(5.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1)))
      ci(23)=6.0d0*(11.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1))
      ci(24)=12.0d0*(dx2+dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>12, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^13 * (x-x2)^n2 as sum_k=0^(13+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_13(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**13
      ci(2)=13.0d0*dx1**12
      ci(3)=78.0d0*dx1**11
      ci(4)=286.0d0*dx1**10
      ci(5)=715.0d0*dx1**9
      ci(6)=1287.0d0*dx1**8
      ci(7)=1716.0d0*dx1**7
      ci(8)=1716.0d0*dx1**6
      ci(9)=1287.0d0*dx1**5
      ci(10)=715.0d0*dx1**4
      ci(11)=286.0d0*dx1**3
      ci(12)=78.0d0*dx1**2
      ci(13)=13.0d0*dx1
      ci(14)=1.0d0
    case(1)
      ci(1)=dx1**13*dx2
      ci(2)=dx1**12*(13.0d0*dx2+dx1)
      ci(3)=13.0d0*dx1**11*(6.0d0*dx2+dx1)
      ci(4)=26.0d0*dx1**10*(11.0d0*dx2+3.0d0*dx1)
      ci(5)=143.0d0*dx1**9*(5.0d0*dx2+2.0d0*dx1)
      ci(6)=143.0d0*dx1**8*(9.0d0*dx2+5.0d0*dx1)
      ci(7)=429.0d0*dx1**7*(4.0d0*dx2+3.0d0*dx1)
      ci(8)=1716.0d0*dx1**6*(dx2+dx1)
      ci(9)=429.0d0*dx1**5*(3.0d0*dx2+4.0d0*dx1)
      ci(10)=143.0d0*dx1**4*(5.0d0*dx2+9.0d0*dx1)
      ci(11)=143.0d0*dx1**3*(2.0d0*dx2+5.0d0*dx1)
      ci(12)=26.0d0*dx1**2*(3.0d0*dx2+11.0d0*dx1)
      ci(13)=13.0d0*dx1*(dx2+6.0d0*dx1)
      ci(14)=dx2+13.0d0*dx1
      ci(15)=1.0d0
    case(2)
      ci(1)=dx1**13*dx2**2
      ci(2)=dx1**12*dx2*(13.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**11*(78.0d0*dx2**2+dx1*(26.0d0*dx2+dx1))
      ci(4)=13.0d0*dx1**10*(22.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))
      ci(5)=13.0d0*dx1**9*(55.0d0*dx2**2+2.0d0*dx1*(22.0d0*dx2+3.0d0*dx1))
      ci(6)=143.0d0*dx1**8*(9.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1))
      ci(7)=143.0d0*dx1**7*(12.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1))
      ci(8)=429.0d0*dx1**6*(4.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(9)=429.0d0*dx1**5*(3.0d0*dx2**2+4.0d0*dx1*(2.0d0*dx2+dx1))
      ci(10)=143.0d0*dx1**4*(5.0d0*dx2**2+6.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(11)=143.0d0*dx1**3*(2.0d0*dx2**2+dx1*(10.0d0*dx2+9.0d0*dx1))
      ci(12)=13.0d0*dx1**2*(6.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+5.0d0*dx1))
      ci(13)=13.0d0*dx1*(dx2**2+2.0d0*dx1*(6.0d0*dx2+11.0d0*dx1))
      ci(14)=dx2**2+26.0d0*dx1*(dx2+3.0d0*dx1)
      ci(15)=2.0d0*dx2+13.0d0*dx1
      ci(16)=1.0d0
    case(3)
      ci(1)=dx1**13*dx2**3
      ci(2)=dx1**12*dx2**2*(13.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**11*dx2*(26.0d0*dx2**2+dx1*(13.0d0*dx2+dx1))
      ci(4)=dx1**10*(286.0d0*dx2**3+dx1*(234.0d0*dx2**2+dx1*(39.0d0*dx2+dx1)))
      ci(5)=13.0d0*dx1**9*(55.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))
      ci(6)=39.0d0*dx1**8*(33.0d0*dx2**3+dx1*(55.0d0*dx2**2+2.0d0*dx1*(11.0d0*dx2+dx1)))
      ci(7)=143.0d0*dx1**7*(12.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(8)=143.0d0*dx1**6*(12.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(27.0d0*dx2+5.0d0*dx1)))
      ci(9)=1287.0d0*dx1**5*(dx2**3+dx1*(4.0d0*dx2**2+dx1*(4.0d0*dx2+dx1)))
      ci(10)=143.0d0*dx1**4*(5.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+4.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(11)=143.0d0*dx1**3*(2.0d0*dx2**3+3.0d0*dx1*(5.0d0*dx2**2+dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(12)=39.0d0*dx1**2*(2.0d0*dx2**3+11.0d0*dx1*(2.0d0*dx2**2+dx1*(5.0d0*dx2+3.0d0*dx1)))
      ci(13)=13.0d0*dx1*(dx2**3+dx1*(18.0d0*dx2**2+11.0d0*dx1*(6.0d0*dx2+5.0d0*dx1)))
      ci(14)=dx2**3+13.0d0*dx1*(3.0d0*dx2**2+2.0d0*dx1*(9.0d0*dx2+11.0d0*dx1))
      ci(15)=3.0d0*(dx2**2+13.0d0*dx1*(dx2+2.0d0*dx1))
      ci(16)=3.0d0*dx2+13.0d0*dx1
      ci(17)=1.0d0
    case(4)
      ci(1)=dx1**13*dx2**4
      ci(2)=dx1**12*dx2**3*(13.0d0*dx2+4.0d0*dx1)
      ci(3)=2.0d0*dx1**11*dx2**2*(39.0d0*dx2**2+dx1*(26.0d0*dx2+3.0d0*dx1))
      ci(4)=2.0d0*dx1**10*dx2*(143.0d0*dx2**3+dx1*(156.0d0*dx2**2+dx1*(39.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**9*(715.0d0*dx2**4+dx1*(1144.0d0*dx2**3+dx1*(468.0d0*dx2**2+dx1*(52.0d0*dx2+dx1))))
      ci(6)=13.0d0*dx1**8*(99.0d0*dx2**4+dx1*(220.0d0*dx2**3+dx1*(132.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(7)=26.0d0*dx1**7*(66.0d0*dx2**4+dx1*(198.0d0*dx2**3+dx1*(165.0d0*dx2**2+dx1*(44.0d0*dx2+3.0d0*dx1))))
      ci(8)=286.0d0*dx1**6*(6.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(9)=143.0d0*dx1**5*(9.0d0*dx2**4+dx1*(48.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))))
      ci(10)=143.0d0*dx1**4*(5.0d0*dx2**4+3.0d0*dx1*(12.0d0*dx2**3+dx1*(24.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(11)=286.0d0*dx1**3*(dx2**4+dx1*(10.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(12)=26.0d0*dx1**2*(3.0d0*dx2**4+11.0d0*dx1*(4.0d0*dx2**3+3.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1))))
      ci(13)=13.0d0*dx1*(dx2**4+dx1*(24.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+dx1*(20.0d0*dx2+9.0d0*dx1))))
      ci(14)=dx2**4+13.0d0*dx1*(4.0d0*dx2**3+dx1*(36.0d0*dx2**2+11.0d0*dx1*(8.0d0*dx2+5.0d0*dx1)))
      ci(15)=2.0d0*(2.0d0*dx2**3+13.0d0*dx1*(3.0d0*dx2**2+dx1*(12.0d0*dx2+11.0d0*dx1)))
      ci(16)=2.0d0*(3.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(17)=4.0d0*dx2+13.0d0*dx1
      ci(18)=1.0d0
    case(5)
      ci(1)=dx1**13*dx2**5
      ci(2)=dx1**12*dx2**4*(13.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**11*dx2**3*(78.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+2.0d0*dx1))
      ci(4)=2.0d0*dx1**10*dx2**2*(143.0d0*dx2**3+5.0d0*dx1*(39.0d0*dx2**2+dx1*(13.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**9*dx2*(143.0d0*dx2**4+dx1*(286.0d0*dx2**3+dx1*(156.0d0*dx2**2+dx1*(26.0d0*dx2+dx1))))
      ci(6)=dx1**8*(1287.0d0*dx2**5+dx1*(3575.0d0*dx2**4+dx1*(2860.0d0*dx2**3+dx1*(780.0d0*dx2**2+dx1*(65.0d0*dx2+dx1)))))
      ci(7)=13.0d0*dx1**7*(132.0d0*dx2**5+dx1*(495.0d0*dx2**4+dx1*(550.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))))
      ci(8)=26.0d0*dx1**6*(66.0d0*dx2**5+dx1*(330.0d0*dx2**4+dx1*(495.0d0*dx2**3+dx1*(275.0d0*dx2**2+dx1*(55.0d0*dx2+3.0d0*dx1)))) &
)
      ci(9)=143.0d0*dx1**5*(9.0d0*dx2**5+dx1*(60.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(25.0d0*dx2+2.0d0*dx1)))))
      ci(10)=715.0d0*dx1**4*(dx2**5+dx1*(9.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(24.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))))
      ci(11)=143.0d0*dx1**3*(2.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(30.0d0*dx2**3+dx1*(40.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1 &
)))))
      ci(12)=26.0d0*dx1**2*(3.0d0*dx2**5+11.0d0*dx1*(5.0d0*dx2**4+dx1*(25.0d0*dx2**3+3.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2 &
+dx1)))))
      ci(13)=13.0d0*dx1*(dx2**5+dx1*(30.0d0*dx2**4+11.0d0*dx1*(20.0d0*dx2**3+dx1*(50.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+4.0d0*dx1))) &
))
      ci(14)=dx2**5+13.0d0*dx1*(5.0d0*dx2**4+dx1*(60.0d0*dx2**3+11.0d0*dx1*(20.0d0*dx2**2+dx1*(25.0d0*dx2+9.0d0*dx1))))
      ci(15)=5.0d0*(dx2**4+13.0d0*dx1*(2.0d0*dx2**3+dx1*(12.0d0*dx2**2+11.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(16)=2.0d0*(5.0d0*dx2**3+13.0d0*dx1*(5.0d0*dx2**2+dx1*(15.0d0*dx2+11.0d0*dx1)))
      ci(17)=10.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+6.0d0*dx1)
      ci(18)=5.0d0*dx2+13.0d0*dx1
      ci(19)=1.0d0
    case(6)
      ci(1)=dx1**13*dx2**6
      ci(2)=dx1**12*dx2**5*(13.0d0*dx2+6.0d0*dx1)
      ci(3)=3.0d0*dx1**11*dx2**4*(26.0d0*dx2**2+dx1*(26.0d0*dx2+5.0d0*dx1))
      ci(4)=dx1**10*dx2**3*(286.0d0*dx2**3+dx1*(468.0d0*dx2**2+5.0d0*dx1*(39.0d0*dx2+4.0d0*dx1)))
      ci(5)=dx1**9*dx2**2*(715.0d0*dx2**4+dx1*(1716.0d0*dx2**3+5.0d0*dx1*(234.0d0*dx2**2+dx1*(52.0d0*dx2+3.0d0*dx1))))
      ci(6)=3.0d0*dx1**8*dx2*(429.0d0*dx2**5+dx1*(1430.0d0*dx2**4+dx1*(1430.0d0*dx2**3+dx1*(520.0d0*dx2**2+dx1*(65.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=dx1**7*(1716.0d0*dx2**6+dx1*(7722.0d0*dx2**5+dx1*(10725.0d0*dx2**4+dx1*(5720.0d0*dx2**3+dx1*(1170.0d0*dx2**2+dx1* &
(78.0d0*dx2+dx1))))))
      ci(8)=13.0d0*dx1**6*(132.0d0*dx2**6+dx1*(792.0d0*dx2**5+dx1*(1485.0d0*dx2**4+dx1*(1100.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1* &
(36.0d0*dx2+dx1))))))
      ci(9)=39.0d0*dx1**5*(33.0d0*dx2**6+dx1*(264.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(660.0d0*dx2**3+dx1*(275.0d0*dx2**2 &
+2.0d0*dx1*(22.0d0*dx2+dx1))))))
      ci(10)=143.0d0*dx1**4*(5.0d0*dx2**6+dx1*(54.0d0*dx2**5+dx1*(180.0d0*dx2**4+dx1*(240.0d0*dx2**3+dx1*(135.0d0*dx2**2 &
+2.0d0*dx1*(15.0d0*dx2+dx1))))))
      ci(11)=143.0d0*dx1**3*(2.0d0*dx2**6+dx1*(30.0d0*dx2**5+dx1*(135.0d0*dx2**4+dx1*(240.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1* &
(54.0d0*dx2+5.0d0*dx1))))))
      ci(12)=39.0d0*dx1**2*(2.0d0*dx2**6+11.0d0*dx1*(4.0d0*dx2**5+dx1*(25.0d0*dx2**4+3.0d0*dx1*(20.0d0*dx2**3+dx1*(20.0d0*dx2**2 &
+dx1*(8.0d0*dx2+dx1))))))
      ci(13)=13.0d0*dx1*(dx2**6+dx1*(36.0d0*dx2**5+11.0d0*dx1*(30.0d0*dx2**4+dx1*(100.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2 &
+4.0d0*dx1*(6.0d0*dx2+dx1))))))
      ci(14)=dx2**6+13.0d0*dx1*(6.0d0*dx2**5+dx1*(90.0d0*dx2**4+11.0d0*dx1*(40.0d0*dx2**3+3.0d0*dx1*(25.0d0*dx2**2+2.0d0*dx1* &
(9.0d0*dx2+2.0d0*dx1)))))
      ci(15)=3.0d0*(2.0d0*dx2**5+13.0d0*dx1*(5.0d0*dx2**4+dx1*(40.0d0*dx2**3+11.0d0*dx1*(10.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1)) &
)))
      ci(16)=15.0d0*dx2**4+13.0d0*dx1*(20.0d0*dx2**3+dx1*(90.0d0*dx2**2+11.0d0*dx1*(12.0d0*dx2+5.0d0*dx1)))
      ci(17)=20.0d0*dx2**3+13.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(18.0d0*dx2+11.0d0*dx1))
      ci(18)=3.0d0*(5.0d0*dx2**2+26.0d0*dx1*(dx2+dx1))
      ci(19)=6.0d0*dx2+13.0d0*dx1
      ci(20)=1.0d0
    case(7)
      ci(1)=dx1**13*dx2**7
      ci(2)=dx1**12*dx2**6*(13.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**11*dx2**5*(78.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+3.0d0*dx1))
      ci(4)=dx1**10*dx2**4*(286.0d0*dx2**3+7.0d0*dx1*(78.0d0*dx2**2+dx1*(39.0d0*dx2+5.0d0*dx1)))
      ci(5)=dx1**9*dx2**3*(715.0d0*dx2**4+7.0d0*dx1*(286.0d0*dx2**3+dx1*(234.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+dx1))))
      ci(6)=dx1**8*dx2**2*(1287.0d0*dx2**5+7.0d0*dx1*(715.0d0*dx2**4+dx1*(858.0d0*dx2**3+dx1*(390.0d0*dx2**2+dx1*(65.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=dx1**7*dx2*(1716.0d0*dx2**6+7.0d0*dx1*(1287.0d0*dx2**5+dx1*(2145.0d0*dx2**4+dx1*(1430.0d0*dx2**3+dx1*(390.0d0*dx2**2 &
+dx1*(39.0d0*dx2+dx1))))))
      ci(8)=dx1**6*(1716.0d0*dx2**7+dx1*(12012.0d0*dx2**6+dx1*(27027.0d0*dx2**5+dx1*(25025.0d0*dx2**4+dx1*(10010.0d0*dx2**3+dx1* &
(1638.0d0*dx2**2+dx1*(91.0d0*dx2+dx1)))))))
      ci(9)=13.0d0*dx1**5*(99.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(3465.0d0*dx2**4+dx1*(1925.0d0*dx2**3+dx1* &
(462.0d0*dx2**2+dx1*(42.0d0*dx2+dx1)))))))
      ci(10)=13.0d0*dx1**4*(55.0d0*dx2**7+dx1*(693.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1*(4620.0d0*dx2**4+dx1*(3465.0d0*dx2**3+dx1* &
(1155.0d0*dx2**2+2.0d0*dx1*(77.0d0*dx2+3.0d0*dx1)))))))
      ci(11)=143.0d0*dx1**3*(2.0d0*dx2**7+dx1*(35.0d0*dx2**6+dx1*(189.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1*(420.0d0*dx2**3+dx1* &
(189.0d0*dx2**2+dx1*(35.0d0*dx2+2.0d0*dx1)))))))
      ci(12)=13.0d0*dx1**2*(6.0d0*dx2**7+11.0d0*dx1*(14.0d0*dx2**6+dx1*(105.0d0*dx2**5+dx1*(315.0d0*dx2**4+dx1*(420.0d0*dx2**3 &
+dx1*(252.0d0*dx2**2+dx1*(63.0d0*dx2+5.0d0*dx1)))))))
      ci(13)=13.0d0*dx1*(dx2**7+dx1*(42.0d0*dx2**6+11.0d0*dx1*(42.0d0*dx2**5+dx1*(175.0d0*dx2**4+3.0d0*dx1*(105.0d0*dx2**3+dx1* &
(84.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1)))))))
      ci(14)=dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(126.0d0*dx2**5+11.0d0*dx1*(70.0d0*dx2**4+dx1*(175.0d0*dx2**3+3.0d0*dx1* &
(63.0d0*dx2**2+4.0d0*dx1*(7.0d0*dx2+dx1))))))
      ci(15)=7.0d0*dx2**6+13.0d0*dx1*(21.0d0*dx2**5+dx1*(210.0d0*dx2**4+11.0d0*dx1*(70.0d0*dx2**3+3.0d0*dx1*(35.0d0*dx2**2+dx1* &
(21.0d0*dx2+4.0d0*dx1)))))
      ci(16)=21.0d0*dx2**5+13.0d0*dx1*(35.0d0*dx2**4+dx1*(210.0d0*dx2**3+11.0d0*dx1*(42.0d0*dx2**2+dx1*(35.0d0*dx2+9.0d0*dx1))))
      ci(17)=35.0d0*dx2**4+13.0d0*dx1*(35.0d0*dx2**3+dx1*(126.0d0*dx2**2+11.0d0*dx1*(14.0d0*dx2+5.0d0*dx1)))
      ci(18)=35.0d0*dx2**3+13.0d0*dx1*(21.0d0*dx2**2+2.0d0*dx1*(21.0d0*dx2+11.0d0*dx1))
      ci(19)=21.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+6.0d0*dx1)
      ci(20)=7.0d0*dx2+13.0d0*dx1
      ci(21)=1.0d0
    case(8)
      ci(1)=dx1**13*dx2**8
      ci(2)=dx1**12*dx2**7*(13.0d0*dx2+8.0d0*dx1)
      ci(3)=2.0d0*dx1**11*dx2**6*(39.0d0*dx2**2+2.0d0*dx1*(26.0d0*dx2+7.0d0*dx1))
      ci(4)=2.0d0*dx1**10*dx2**5*(143.0d0*dx2**3+2.0d0*dx1*(156.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**9*dx2**4*(715.0d0*dx2**4+2.0d0*dx1*(1144.0d0*dx2**3+7.0d0*dx1*(156.0d0*dx2**2+dx1*(52.0d0*dx2+5.0d0*dx1))))
      ci(6)=dx1**8*dx2**3*(1287.0d0*dx2**5+2.0d0*dx1*(2860.0d0*dx2**4+7.0d0*dx1*(572.0d0*dx2**3+dx1*(312.0d0*dx2**2+dx1* &
(65.0d0*dx2+4.0d0*dx1)))))
      ci(7)=4.0d0*dx1**7*dx2**2*(429.0d0*dx2**6+dx1*(2574.0d0*dx2**5+7.0d0*dx1*(715.0d0*dx2**4+dx1*(572.0d0*dx2**3+dx1* &
(195.0d0*dx2**2+dx1*(26.0d0*dx2+dx1))))))
      ci(8)=4.0d0*dx1**6*dx2*(429.0d0*dx2**7+dx1*(3432.0d0*dx2**6+dx1*(9009.0d0*dx2**5+dx1*(10010.0d0*dx2**4+dx1*(5005.0d0*dx2**3 &
+dx1*(1092.0d0*dx2**2+dx1*(91.0d0*dx2+2.0d0*dx1)))))))
      ci(9)=dx1**5*(1287.0d0*dx2**8+dx1*(13728.0d0*dx2**7+dx1*(48048.0d0*dx2**6+dx1*(72072.0d0*dx2**5+dx1*(50050.0d0*dx2**4+dx1* &
(16016.0d0*dx2**3+dx1*(2184.0d0*dx2**2+dx1*(104.0d0*dx2+dx1))))))))
      ci(10)=13.0d0*dx1**4*(55.0d0*dx2**8+dx1*(792.0d0*dx2**7+dx1*(3696.0d0*dx2**6+dx1*(7392.0d0*dx2**5+dx1*(6930.0d0*dx2**4+dx1* &
(3080.0d0*dx2**3+dx1*(616.0d0*dx2**2+dx1*(48.0d0*dx2+dx1))))))))
      ci(11)=26.0d0*dx1**3*(11.0d0*dx2**8+dx1*(220.0d0*dx2**7+dx1*(1386.0d0*dx2**6+dx1*(3696.0d0*dx2**5+dx1*(4620.0d0*dx2**4+dx1* &
(2772.0d0*dx2**3+dx1*(770.0d0*dx2**2+dx1*(88.0d0*dx2+3.0d0*dx1))))))))
      ci(12)=26.0d0*dx1**2*(3.0d0*dx2**8+11.0d0*dx1*(8.0d0*dx2**7+dx1*(70.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1*(420.0d0*dx2**4+dx1* &
(336.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))))))
      ci(13)=13.0d0*dx1*(dx2**8+dx1*(48.0d0*dx2**7+11.0d0*dx1*(56.0d0*dx2**6+dx1*(280.0d0*dx2**5+dx1*(630.0d0*dx2**4+dx1* &
(672.0d0*dx2**3+dx1*(336.0d0*dx2**2+dx1*(72.0d0*dx2+5.0d0*dx1))))))))
      ci(14)=dx2**8+13.0d0*dx1*(8.0d0*dx2**7+dx1*(168.0d0*dx2**6+11.0d0*dx1*(112.0d0*dx2**5+dx1*(350.0d0*dx2**4+3.0d0*dx1* &
(168.0d0*dx2**3+dx1*(112.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1)))))))
      ci(15)=4.0d0*(2.0d0*dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(84.0d0*dx2**5+11.0d0*dx1*(35.0d0*dx2**4+dx1*(70.0d0*dx2**3 &
+3.0d0*dx1*(21.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))))))
      ci(16)=4.0d0*(7.0d0*dx2**6+13.0d0*dx1*(14.0d0*dx2**5+dx1*(105.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+dx1*(35.0d0*dx2**2 &
+3.0d0*dx1*(6.0d0*dx2+dx1))))))
      ci(17)=56.0d0*dx2**5+13.0d0*dx1*(70.0d0*dx2**4+dx1*(336.0d0*dx2**3+11.0d0*dx1*(56.0d0*dx2**2+dx1*(40.0d0*dx2+9.0d0*dx1))))
      ci(18)=70.0d0*dx2**4+13.0d0*dx1*(56.0d0*dx2**3+dx1*(168.0d0*dx2**2+11.0d0*dx1*(16.0d0*dx2+5.0d0*dx1)))
      ci(19)=2.0d0*(28.0d0*dx2**3+13.0d0*dx1*(14.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1)))
      ci(20)=2.0d0*(14.0d0*dx2**2+13.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(21)=8.0d0*dx2+13.0d0*dx1
      ci(22)=1.0d0
    case(9)
      ci(1)=dx1**13*dx2**9
      ci(2)=dx1**12*dx2**8*(13.0d0*dx2+9.0d0*dx1)
      ci(3)=3.0d0*dx1**11*dx2**7*(26.0d0*dx2**2+3.0d0*dx1*(13.0d0*dx2+4.0d0*dx1))
      ci(4)=2.0d0*dx1**10*dx2**6*(143.0d0*dx2**3+3.0d0*dx1*(117.0d0*dx2**2+2.0d0*dx1*(39.0d0*dx2+7.0d0*dx1)))
      ci(5)=dx1**9*dx2**5*(715.0d0*dx2**4+6.0d0*dx1*(429.0d0*dx2**3+dx1*(468.0d0*dx2**2+7.0d0*dx1*(26.0d0*dx2+3.0d0*dx1))))
      ci(6)=9.0d0*dx1**8*dx2**4*(143.0d0*dx2**5+dx1*(715.0d0*dx2**4+2.0d0*dx1*(572.0d0*dx2**3+7.0d0*dx1*(52.0d0*dx2**2+dx1* &
(13.0d0*dx2+dx1)))))
      ci(7)=3.0d0*dx1**7*dx2**3*(572.0d0*dx2**6+dx1*(3861.0d0*dx2**5+2.0d0*dx1*(4290.0d0*dx2**4+7.0d0*dx1*(572.0d0*dx2**3+dx1* &
(234.0d0*dx2**2+dx1*(39.0d0*dx2+2.0d0*dx1))))))
      ci(8)=12.0d0*dx1**6*dx2**2*(143.0d0*dx2**7+dx1*(1287.0d0*dx2**6+dx1*(3861.0d0*dx2**5+dx1*(5005.0d0*dx2**4+dx1* &
(3003.0d0*dx2**3+dx1*(819.0d0*dx2**2+dx1*(91.0d0*dx2+3.0d0*dx1)))))))
      ci(9)=9.0d0*dx1**5*dx2*(143.0d0*dx2**8+dx1*(1716.0d0*dx2**7+dx1*(6864.0d0*dx2**6+dx1*(12012.0d0*dx2**5+dx1*(10010.0d0*dx2**4 &
+dx1*(4004.0d0*dx2**3+dx1*(728.0d0*dx2**2+dx1*(52.0d0*dx2+dx1))))))))
      ci(10)=dx1**4*(715.0d0*dx2**9+dx1*(11583.0d0*dx2**8+dx1*(61776.0d0*dx2**7+dx1*(144144.0d0*dx2**6+dx1*(162162.0d0*dx2**5+dx1* &
(90090.0d0*dx2**4+dx1*(24024.0d0*dx2**3+dx1*(2808.0d0*dx2**2+dx1*(117.0d0*dx2+dx1)))))))))
      ci(11)=13.0d0*dx1**3*(22.0d0*dx2**9+dx1*(495.0d0*dx2**8+dx1*(3564.0d0*dx2**7+dx1*(11088.0d0*dx2**6+dx1*(16632.0d0*dx2**5 &
+dx1*(12474.0d0*dx2**4+dx1*(4620.0d0*dx2**3+dx1*(792.0d0*dx2**2+dx1*(54.0d0*dx2+dx1)))))))))
      ci(12)=78.0d0*dx1**2*(dx2**9+dx1*(33.0d0*dx2**8+dx1*(330.0d0*dx2**7+dx1*(1386.0d0*dx2**6+dx1*(2772.0d0*dx2**5+dx1* &
(2772.0d0*dx2**4+dx1*(1386.0d0*dx2**3+dx1*(330.0d0*dx2**2+dx1*(33.0d0*dx2+dx1)))))))))
      ci(13)=13.0d0*dx1*(dx2**9+dx1*(54.0d0*dx2**8+11.0d0*dx1*(72.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(1134.0d0*dx2**5+dx1* &
(1512.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(324.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))))))))
      ci(14)=dx2**9+13.0d0*dx1*(9.0d0*dx2**8+dx1*(216.0d0*dx2**7+11.0d0*dx1*(168.0d0*dx2**6+dx1*(630.0d0*dx2**5+dx1* &
(1134.0d0*dx2**4+dx1*(1008.0d0*dx2**3+dx1*(432.0d0*dx2**2+dx1*(81.0d0*dx2+5.0d0*dx1))))))))
      ci(15)=9.0d0*(dx2**8+13.0d0*dx1*(4.0d0*dx2**7+dx1*(56.0d0*dx2**6+11.0d0*dx1*(28.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1* &
(84.0d0*dx2**3+dx1*(48.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))))))
      ci(16)=12.0d0*(3.0d0*dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(63.0d0*dx2**5+11.0d0*dx1*(21.0d0*dx2**4+dx1*(35.0d0*dx2**3+dx1* &
(27.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))))))
      ci(17)=3.0d0*(28.0d0*dx2**6+13.0d0*dx1*(42.0d0*dx2**5+dx1*(252.0d0*dx2**4+11.0d0*dx1*(56.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1* &
(27.0d0*dx2+4.0d0*dx1))))))
      ci(18)=9.0d0*(14.0d0*dx2**5+13.0d0*dx1*(14.0d0*dx2**4+dx1*(56.0d0*dx2**3+11.0d0*dx1*(8.0d0*dx2**2+dx1*(5.0d0*dx2+dx1)))))
      ci(19)=126.0d0*dx2**4+13.0d0*dx1*(84.0d0*dx2**3+dx1*(216.0d0*dx2**2+11.0d0*dx1*(18.0d0*dx2+5.0d0*dx1)))
      ci(20)=2.0d0*(42.0d0*dx2**3+13.0d0*dx1*(18.0d0*dx2**2+dx1*(27.0d0*dx2+11.0d0*dx1)))
      ci(21)=3.0d0*(12.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(22)=9.0d0*dx2+13.0d0*dx1
      ci(23)=1.0d0
    case(10)
      ci(1)=dx1**13*dx2**10
      ci(2)=dx1**12*dx2**9*(13.0d0*dx2+10.0d0*dx1)
      ci(3)=dx1**11*dx2**8*(78.0d0*dx2**2+5.0d0*dx1*(26.0d0*dx2+9.0d0*dx1))
      ci(4)=dx1**10*dx2**7*(286.0d0*dx2**3+15.0d0*dx1*(52.0d0*dx2**2+dx1*(39.0d0*dx2+8.0d0*dx1)))
      ci(5)=5.0d0*dx1**9*dx2**6*(143.0d0*dx2**4+2.0d0*dx1*(286.0d0*dx2**3+3.0d0*dx1*(117.0d0*dx2**2+dx1*(52.0d0*dx2+7.0d0*dx1))))
      ci(6)=dx1**8*dx2**5*(1287.0d0*dx2**5+2.0d0*dx1*(3575.0d0*dx2**4+3.0d0*dx1*(2145.0d0*dx2**3+dx1*(1560.0d0*dx2**2+7.0d0*dx1* &
(65.0d0*dx2+6.0d0*dx1)))))
      ci(7)=3.0d0*dx1**7*dx2**4*(572.0d0*dx2**6+dx1*(4290.0d0*dx2**5+dx1*(10725.0d0*dx2**4+2.0d0*dx1*(5720.0d0*dx2**3+7.0d0*dx1* &
(390.0d0*dx2**2+dx1*(78.0d0*dx2+5.0d0*dx1))))))
      ci(8)=3.0d0*dx1**6*dx2**3*(572.0d0*dx2**7+dx1*(5720.0d0*dx2**6+dx1*(19305.0d0*dx2**5+2.0d0*dx1*(14300.0d0*dx2**4+dx1* &
(10010.0d0*dx2**3+dx1*(3276.0d0*dx2**2+5.0d0*dx1*(91.0d0*dx2+4.0d0*dx1)))))))
      ci(9)=3.0d0*dx1**5*dx2**2*(429.0d0*dx2**8+dx1*(5720.0d0*dx2**7+dx1*(25740.0d0*dx2**6+dx1*(51480.0d0*dx2**5+dx1* &
(50050.0d0*dx2**4+dx1*(24024.0d0*dx2**3+5.0d0*dx1*(1092.0d0*dx2**2+dx1*(104.0d0*dx2+3.0d0*dx1))))))))
      ci(10)=5.0d0*dx1**4*dx2*(143.0d0*dx2**9+dx1*(2574.0d0*dx2**8+dx1*(15444.0d0*dx2**7+dx1*(41184.0d0*dx2**6+dx1* &
(54054.0d0*dx2**5+dx1*(36036.0d0*dx2**4+dx1*(12012.0d0*dx2**3+dx1*(1872.0d0*dx2**2+dx1*(117.0d0*dx2+2.0d0*dx1)))))))))
      ci(11)=dx1**3*(286.0d0*dx2**10+dx1*(7150.0d0*dx2**9+dx1*(57915.0d0*dx2**8+dx1*(205920.0d0*dx2**7+dx1*(360360.0d0*dx2**6+dx1* &
(324324.0d0*dx2**5+dx1*(150150.0d0*dx2**4+dx1*(34320.0d0*dx2**3+dx1*(3510.0d0*dx2**2+dx1*(130.0d0*dx2+dx1))))))))))
      ci(12)=13.0d0*dx1**2*(6.0d0*dx2**10+dx1*(220.0d0*dx2**9+dx1*(2475.0d0*dx2**8+dx1*(11880.0d0*dx2**7+dx1*(27720.0d0*dx2**6 &
+dx1*(33264.0d0*dx2**5+dx1*(20790.0d0*dx2**4+dx1*(6600.0d0*dx2**3+dx1*(990.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))))))))
      ci(13)=13.0d0*dx1*(dx2**10+dx1*(60.0d0*dx2**9+dx1*(990.0d0*dx2**8+dx1*(6600.0d0*dx2**7+dx1*(20790.0d0*dx2**6+dx1* &
(33264.0d0*dx2**5+dx1*(27720.0d0*dx2**4+dx1*(11880.0d0*dx2**3+dx1*(2475.0d0*dx2**2+2.0d0*dx1*(110.0d0*dx2+3.0d0*dx1))))))))))
      ci(14)=dx2**10+13.0d0*dx1*(10.0d0*dx2**9+dx1*(270.0d0*dx2**8+11.0d0*dx1*(240.0d0*dx2**7+dx1*(1050.0d0*dx2**6+dx1* &
(2268.0d0*dx2**5+dx1*(2520.0d0*dx2**4+dx1*(1440.0d0*dx2**3+dx1*(405.0d0*dx2**2+2.0d0*dx1*(25.0d0*dx2+dx1)))))))))
      ci(15)=5.0d0*(2.0d0*dx2**9+13.0d0*dx1*(9.0d0*dx2**8+dx1*(144.0d0*dx2**7+11.0d0*dx1*(84.0d0*dx2**6+dx1*(252.0d0*dx2**5+dx1* &
(378.0d0*dx2**4+dx1*(288.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))))))))
      ci(16)=3.0d0*(15.0d0*dx2**8+13.0d0*dx1*(40.0d0*dx2**7+dx1*(420.0d0*dx2**6+11.0d0*dx1*(168.0d0*dx2**5+dx1*(350.0d0*dx2**4 &
+dx1*(360.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))))))))
      ci(17)=3.0d0*(40.0d0*dx2**7+13.0d0*dx1*(70.0d0*dx2**6+dx1*(504.0d0*dx2**5+11.0d0*dx1*(140.0d0*dx2**4+dx1*(200.0d0*dx2**3 &
+dx1*(135.0d0*dx2**2+4.0d0*dx1*(10.0d0*dx2+dx1)))))))
      ci(18)=3.0d0*(70.0d0*dx2**6+13.0d0*dx1*(84.0d0*dx2**5+dx1*(420.0d0*dx2**4+11.0d0*dx1*(80.0d0*dx2**3+dx1*(75.0d0*dx2**2 &
+2.0d0*dx1*(15.0d0*dx2+2.0d0*dx1))))))
      ci(19)=252.0d0*dx2**5+13.0d0*dx1*(210.0d0*dx2**4+dx1*(720.0d0*dx2**3+11.0d0*dx1*(90.0d0*dx2**2+dx1*(50.0d0*dx2+9.0d0*dx1))))
      ci(20)=5.0d0*(42.0d0*dx2**4+13.0d0*dx1*(24.0d0*dx2**3+dx1*(54.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(21)=120.0d0*dx2**3+13.0d0*dx1*(45.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2+11.0d0*dx1))
      ci(22)=45.0d0*dx2**2+26.0d0*dx1*(5.0d0*dx2+3.0d0*dx1)
      ci(23)=10.0d0*dx2+13.0d0*dx1
      ci(24)=1.0d0
    case(11)
      ci(1)=dx1**13*dx2**11
      ci(2)=dx1**12*dx2**10*(13.0d0*dx2+11.0d0*dx1)
      ci(3)=dx1**11*dx2**9*(78.0d0*dx2**2+11.0d0*dx1*(13.0d0*dx2+5.0d0*dx1))
      ci(4)=11.0d0*dx1**10*dx2**8*(26.0d0*dx2**3+dx1*(78.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+3.0d0*dx1)))
      ci(5)=11.0d0*dx1**9*dx2**7*(65.0d0*dx2**4+dx1*(286.0d0*dx2**3+15.0d0*dx1*(26.0d0*dx2**2+dx1*(13.0d0*dx2+2.0d0*dx1))))
      ci(6)=11.0d0*dx1**8*dx2**6*(117.0d0*dx2**5+dx1*(715.0d0*dx2**4+2.0d0*dx1*(715.0d0*dx2**3+3.0d0*dx1*(195.0d0*dx2**2+dx1* &
(65.0d0*dx2+7.0d0*dx1)))))
      ci(7)=11.0d0*dx1**7*dx2**5*(156.0d0*dx2**6+dx1*(1287.0d0*dx2**5+dx1*(3575.0d0*dx2**4+6.0d0*dx1*(715.0d0*dx2**3+dx1* &
(390.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+dx1))))))
      ci(8)=33.0d0*dx1**6*dx2**4*(52.0d0*dx2**7+dx1*(572.0d0*dx2**6+dx1*(2145.0d0*dx2**5+dx1*(3575.0d0*dx2**4+2.0d0*dx1* &
(1430.0d0*dx2**3+dx1*(546.0d0*dx2**2+dx1*(91.0d0*dx2+5.0d0*dx1)))))))
      ci(9)=33.0d0*dx1**5*dx2**3*(39.0d0*dx2**8+dx1*(572.0d0*dx2**7+dx1*(2860.0d0*dx2**6+dx1*(6435.0d0*dx2**5+dx1*(7150.0d0*dx2**4 &
+dx1*(4004.0d0*dx2**3+dx1*(1092.0d0*dx2**2+5.0d0*dx1*(26.0d0*dx2+dx1))))))))
      ci(10)=11.0d0*dx1**4*dx2**2*(65.0d0*dx2**9+dx1*(1287.0d0*dx2**8+dx1*(8580.0d0*dx2**7+dx1*(25740.0d0*dx2**6+dx1* &
(38610.0d0*dx2**5+dx1*(30030.0d0*dx2**4+dx1*(12012.0d0*dx2**3+5.0d0*dx1*(468.0d0*dx2**2+dx1*(39.0d0*dx2+dx1)))))))))
      ci(11)=11.0d0*dx1**3*dx2*(26.0d0*dx2**10+dx1*(715.0d0*dx2**9+dx1*(6435.0d0*dx2**8+dx1*(25740.0d0*dx2**7+dx1* &
(51480.0d0*dx2**6+dx1*(54054.0d0*dx2**5+dx1*(30030.0d0*dx2**4+dx1*(8580.0d0*dx2**3+dx1*(1170.0d0*dx2**2+dx1*(65.0d0*dx2+dx1))))))) &
)))
      ci(12)=dx1**2*(78.0d0*dx2**11+dx1*(3146.0d0*dx2**10+dx1*(39325.0d0*dx2**9+dx1*(212355.0d0*dx2**8+dx1*(566280.0d0*dx2**7+dx1* &
(792792.0d0*dx2**6+dx1*(594594.0d0*dx2**5+dx1*(235950.0d0*dx2**4+dx1*(47190.0d0*dx2**3+dx1*(4290.0d0*dx2**2+dx1*(143.0d0*dx2+dx1)) &
)))))))))
      ci(13)=13.0d0*dx1*(dx2**11+dx1*(66.0d0*dx2**10+dx1*(1210.0d0*dx2**9+dx1*(9075.0d0*dx2**8+dx1*(32670.0d0*dx2**7+dx1* &
(60984.0d0*dx2**6+dx1*(60984.0d0*dx2**5+dx1*(32670.0d0*dx2**4+dx1*(9075.0d0*dx2**3+dx1*(1210.0d0*dx2**2+dx1*(66.0d0*dx2+dx1))))))) &
))))
      ci(14)=dx2**11+13.0d0*dx1*(11.0d0*dx2**10+dx1*(330.0d0*dx2**9+dx1*(3630.0d0*dx2**8+dx1*(18150.0d0*dx2**7+dx1* &
(45738.0d0*dx2**6+dx1*(60984.0d0*dx2**5+dx1*(43560.0d0*dx2**4+dx1*(16335.0d0*dx2**3+dx1*(3025.0d0*dx2**2+2.0d0*dx1*(121.0d0*dx2 &
+3.0d0*dx1))))))))))
      ci(15)=11.0d0*(dx2**10+13.0d0*dx1*(5.0d0*dx2**9+dx1*(90.0d0*dx2**8+dx1*(660.0d0*dx2**7+dx1*(2310.0d0*dx2**6+dx1* &
(4158.0d0*dx2**5+dx1*(3960.0d0*dx2**4+dx1*(1980.0d0*dx2**3+dx1*(495.0d0*dx2**2+dx1*(55.0d0*dx2+2.0d0*dx1))))))))))
      ci(16)=11.0d0*(5.0d0*dx2**9+13.0d0*dx1*(15.0d0*dx2**8+dx1*(180.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2310.0d0*dx2**5+dx1* &
(2970.0d0*dx2**4+dx1*(1980.0d0*dx2**3+dx1*(660.0d0*dx2**2+dx1*(99.0d0*dx2+5.0d0*dx1)))))))))
      ci(17)=33.0d0*(5.0d0*dx2**8+13.0d0*dx1*(10.0d0*dx2**7+dx1*(84.0d0*dx2**6+dx1*(308.0d0*dx2**5+dx1*(550.0d0*dx2**4+dx1* &
(495.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(44.0d0*dx2+3.0d0*dx1))))))))
      ci(18)=33.0d0*(10.0d0*dx2**7+13.0d0*dx1*(14.0d0*dx2**6+dx1*(84.0d0*dx2**5+dx1*(220.0d0*dx2**4+dx1*(275.0d0*dx2**3+dx1* &
(165.0d0*dx2**2+4.0d0*dx1*(11.0d0*dx2+dx1)))))))
      ci(19)=11.0d0*(42.0d0*dx2**6+13.0d0*dx1*(42.0d0*dx2**5+dx1*(180.0d0*dx2**4+dx1*(330.0d0*dx2**3+dx1*(275.0d0*dx2**2 &
+3.0d0*dx1*(33.0d0*dx2+4.0d0*dx1))))))
      ci(20)=11.0d0*(42.0d0*dx2**5+13.0d0*dx1*(30.0d0*dx2**4+dx1*(90.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(55.0d0*dx2+9.0d0*dx1)))))
      ci(21)=11.0d0*(30.0d0*dx2**4+13.0d0*dx1*(15.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(22.0d0*dx2+5.0d0*dx1))))
      ci(22)=11.0d0*(15.0d0*dx2**3+13.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(23)=55.0d0*dx2**2+13.0d0*dx1*(11.0d0*dx2+6.0d0*dx1)
      ci(24)=11.0d0*dx2+13.0d0*dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>11, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^14 * (x-x2)^n2 as sum_k=0^(14+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_14(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**14
      ci(2)=14.0d0*dx1**13
      ci(3)=91.0d0*dx1**12
      ci(4)=364.0d0*dx1**11
      ci(5)=1001.0d0*dx1**10
      ci(6)=2002.0d0*dx1**9
      ci(7)=3003.0d0*dx1**8
      ci(8)=3432.0d0*dx1**7
      ci(9)=3003.0d0*dx1**6
      ci(10)=2002.0d0*dx1**5
      ci(11)=1001.0d0*dx1**4
      ci(12)=364.0d0*dx1**3
      ci(13)=91.0d0*dx1**2
      ci(14)=14.0d0*dx1
      ci(15)=1.0d0
    case(1)
      ci(1)=dx1**14*dx2
      ci(2)=dx1**13*(14.0d0*dx2+dx1)
      ci(3)=7.0d0*dx1**12*(13.0d0*dx2+2.0d0*dx1)
      ci(4)=91.0d0*dx1**11*(4.0d0*dx2+dx1)
      ci(5)=91.0d0*dx1**10*(11.0d0*dx2+4.0d0*dx1)
      ci(6)=1001.0d0*dx1**9*(2.0d0*dx2+dx1)
      ci(7)=1001.0d0*dx1**8*(3.0d0*dx2+2.0d0*dx1)
      ci(8)=429.0d0*dx1**7*(8.0d0*dx2+7.0d0*dx1)
      ci(9)=429.0d0*dx1**6*(7.0d0*dx2+8.0d0*dx1)
      ci(10)=1001.0d0*dx1**5*(2.0d0*dx2+3.0d0*dx1)
      ci(11)=1001.0d0*dx1**4*(dx2+2.0d0*dx1)
      ci(12)=91.0d0*dx1**3*(4.0d0*dx2+11.0d0*dx1)
      ci(13)=91.0d0*dx1**2*(dx2+4.0d0*dx1)
      ci(14)=7.0d0*dx1*(2.0d0*dx2+13.0d0*dx1)
      ci(15)=dx2+14.0d0*dx1
      ci(16)=1.0d0
    case(2)
      ci(1)=dx1**14*dx2**2
      ci(2)=2.0d0*dx1**13*dx2*(7.0d0*dx2+dx1)
      ci(3)=dx1**12*(91.0d0*dx2**2+dx1*(28.0d0*dx2+dx1))
      ci(4)=14.0d0*dx1**11*(26.0d0*dx2**2+dx1*(13.0d0*dx2+dx1))
      ci(5)=91.0d0*dx1**10*(11.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(6)=182.0d0*dx1**9*(11.0d0*dx2**2+dx1*(11.0d0*dx2+2.0d0*dx1))
      ci(7)=1001.0d0*dx1**8*(3.0d0*dx2**2+dx1*(4.0d0*dx2+dx1))
      ci(8)=286.0d0*dx1**7*(12.0d0*dx2**2+7.0d0*dx1*(3.0d0*dx2+dx1))
      ci(9)=429.0d0*dx1**6*(7.0d0*dx2**2+dx1*(16.0d0*dx2+7.0d0*dx1))
      ci(10)=286.0d0*dx1**5*(7.0d0*dx2**2+3.0d0*dx1*(7.0d0*dx2+4.0d0*dx1))
      ci(11)=1001.0d0*dx1**4*(dx2**2+dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(12)=182.0d0*dx1**3*(2.0d0*dx2**2+11.0d0*dx1*(dx2+dx1))
      ci(13)=91.0d0*dx1**2*(dx2**2+dx1*(8.0d0*dx2+11.0d0*dx1))
      ci(14)=14.0d0*dx1*(dx2**2+13.0d0*dx1*(dx2+2.0d0*dx1))
      ci(15)=dx2**2+7.0d0*dx1*(4.0d0*dx2+13.0d0*dx1)
      ci(16)=2.0d0*(dx2+7.0d0*dx1)
      ci(17)=1.0d0
    case(3)
      ci(1)=dx1**14*dx2**3
      ci(2)=dx1**13*dx2**2*(14.0d0*dx2+3.0d0*dx1)
      ci(3)=dx1**12*dx2*(91.0d0*dx2**2+3.0d0*dx1*(14.0d0*dx2+dx1))
      ci(4)=dx1**11*(364.0d0*dx2**3+dx1*(273.0d0*dx2**2+dx1*(42.0d0*dx2+dx1)))
      ci(5)=7.0d0*dx1**10*(143.0d0*dx2**3+dx1*(156.0d0*dx2**2+dx1*(39.0d0*dx2+2.0d0*dx1)))
      ci(6)=91.0d0*dx1**9*(22.0d0*dx2**3+dx1*(33.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(7)=91.0d0*dx1**8*(33.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(33.0d0*dx2+4.0d0*dx1)))
      ci(8)=143.0d0*dx1**7*(24.0d0*dx2**3+7.0d0*dx1*(9.0d0*dx2**2+dx1*(6.0d0*dx2+dx1)))
      ci(9)=143.0d0*dx1**6*(21.0d0*dx2**3+dx1*(72.0d0*dx2**2+7.0d0*dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(10)=143.0d0*dx1**5*(14.0d0*dx2**3+3.0d0*dx1*(21.0d0*dx2**2+dx1*(24.0d0*dx2+7.0d0*dx1)))
      ci(11)=143.0d0*dx1**4*(7.0d0*dx2**3+3.0d0*dx1*(14.0d0*dx2**2+dx1*(21.0d0*dx2+8.0d0*dx1)))
      ci(12)=91.0d0*dx1**3*(4.0d0*dx2**3+33.0d0*dx1*(dx2**2+dx1*(2.0d0*dx2+dx1)))
      ci(13)=91.0d0*dx1**2*(dx2**3+dx1*(12.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(14)=7.0d0*dx1*(2.0d0*dx2**3+13.0d0*dx1*(3.0d0*dx2**2+dx1*(12.0d0*dx2+11.0d0*dx1)))
      ci(15)=dx2**3+7.0d0*dx1*(6.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+4.0d0*dx1))
      ci(16)=3.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+13.0d0*dx1)
      ci(17)=3.0d0*dx2+14.0d0*dx1
      ci(18)=1.0d0
    case(4)
      ci(1)=dx1**14*dx2**4
      ci(2)=2.0d0*dx1**13*dx2**3*(7.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**12*dx2**2*(91.0d0*dx2**2+2.0d0*dx1*(28.0d0*dx2+3.0d0*dx1))
      ci(4)=4.0d0*dx1**11*dx2*(91.0d0*dx2**3+dx1*(91.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)))
      ci(5)=dx1**10*(1001.0d0*dx2**4+dx1*(1456.0d0*dx2**3+dx1*(546.0d0*dx2**2+dx1*(56.0d0*dx2+dx1))))
      ci(6)=14.0d0*dx1**9*(143.0d0*dx2**4+dx1*(286.0d0*dx2**3+dx1*(156.0d0*dx2**2+dx1*(26.0d0*dx2+dx1))))
      ci(7)=91.0d0*dx1**8*(33.0d0*dx2**4+dx1*(88.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(8)=52.0d0*dx1**7*(66.0d0*dx2**4+7.0d0*dx1*(33.0d0*dx2**3+dx1*(33.0d0*dx2**2+dx1*(11.0d0*dx2+dx1))))
      ci(9)=143.0d0*dx1**6*(21.0d0*dx2**4+dx1*(96.0d0*dx2**3+7.0d0*dx1*(18.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))))
      ci(10)=286.0d0*dx1**5*(7.0d0*dx2**4+dx1*(42.0d0*dx2**3+dx1*(72.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+dx1))))
      ci(11)=143.0d0*dx1**4*(7.0d0*dx2**4+dx1*(56.0d0*dx2**3+3.0d0*dx1*(42.0d0*dx2**2+dx1*(32.0d0*dx2+7.0d0*dx1))))
      ci(12)=52.0d0*dx1**3*(7.0d0*dx2**4+11.0d0*dx1*(7.0d0*dx2**3+3.0d0*dx1*(7.0d0*dx2**2+dx1*(7.0d0*dx2+2.0d0*dx1))))
      ci(13)=91.0d0*dx1**2*(dx2**4+dx1*(16.0d0*dx2**3+11.0d0*dx1*(6.0d0*dx2**2+dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(14)=14.0d0*dx1*(dx2**4+13.0d0*dx1*(2.0d0*dx2**3+dx1*(12.0d0*dx2**2+11.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(15)=dx2**4+7.0d0*dx1*(8.0d0*dx2**3+13.0d0*dx1*(6.0d0*dx2**2+dx1*(16.0d0*dx2+11.0d0*dx1)))
      ci(16)=4.0d0*(dx2**3+7.0d0*dx1*(3.0d0*dx2**2+13.0d0*dx1*(dx2+dx1)))
      ci(17)=6.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+13.0d0*dx1)
      ci(18)=2.0d0*(2.0d0*dx2+7.0d0*dx1)
      ci(19)=1.0d0
    case(5)
      ci(1)=dx1**14*dx2**5
      ci(2)=dx1**13*dx2**4*(14.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**12*dx2**3*(91.0d0*dx2**2+10.0d0*dx1*(7.0d0*dx2+dx1))
      ci(4)=dx1**11*dx2**2*(364.0d0*dx2**3+5.0d0*dx1*(91.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+dx1)))
      ci(5)=dx1**10*dx2*(1001.0d0*dx2**4+5.0d0*dx1*(364.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(28.0d0*dx2+dx1))))
      ci(6)=dx1**9*(2002.0d0*dx2**5+dx1*(5005.0d0*dx2**4+dx1*(3640.0d0*dx2**3+dx1*(910.0d0*dx2**2+dx1*(70.0d0*dx2+dx1)))))
      ci(7)=7.0d0*dx1**8*(429.0d0*dx2**5+dx1*(1430.0d0*dx2**4+dx1*(1430.0d0*dx2**3+dx1*(520.0d0*dx2**2+dx1*(65.0d0*dx2+2.0d0*dx1)) &
)))
      ci(8)=13.0d0*dx1**7*(264.0d0*dx2**5+7.0d0*dx1*(165.0d0*dx2**4+dx1*(220.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))) &
))
      ci(9)=13.0d0*dx1**6*(231.0d0*dx2**5+dx1*(1320.0d0*dx2**4+7.0d0*dx1*(330.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(55.0d0*dx2 &
+4.0d0*dx1)))))
      ci(10)=143.0d0*dx1**5*(14.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(240.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))) &
))
      ci(11)=143.0d0*dx1**4*(7.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(210.0d0*dx2**3+dx1*(240.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2 &
+2.0d0*dx1)))))
      ci(12)=13.0d0*dx1**3*(28.0d0*dx2**5+11.0d0*dx1*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+3.0d0*dx1*(70.0d0*dx2**2+dx1*(40.0d0*dx2 &
+7.0d0*dx1)))))
      ci(13)=13.0d0*dx1**2*(7.0d0*dx2**5+dx1*(140.0d0*dx2**4+11.0d0*dx1*(70.0d0*dx2**3+dx1*(140.0d0*dx2**2+3.0d0*dx1*(35.0d0*dx2 &
+8.0d0*dx1)))))
      ci(14)=7.0d0*dx1*(2.0d0*dx2**5+13.0d0*dx1*(5.0d0*dx2**4+dx1*(40.0d0*dx2**3+11.0d0*dx1*(10.0d0*dx2**2+dx1*(10.0d0*dx2 &
+3.0d0*dx1)))))
      ci(15)=dx2**5+7.0d0*dx1*(10.0d0*dx2**4+13.0d0*dx1*(10.0d0*dx2**3+dx1*(40.0d0*dx2**2+11.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))))
      ci(16)=5.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+13.0d0*dx1*(10.0d0*dx2**2+dx1*(20.0d0*dx2+11.0d0*dx1)))
      ci(17)=10.0d0*dx2**3+7.0d0*dx1*(20.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+4.0d0*dx1))
      ci(18)=10.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+13.0d0*dx1)
      ci(19)=5.0d0*dx2+14.0d0*dx1
      ci(20)=1.0d0
    case(6)
      ci(1)=dx1**14*dx2**6
      ci(2)=2.0d0*dx1**13*dx2**5*(7.0d0*dx2+3.0d0*dx1)
      ci(3)=dx1**12*dx2**4*(91.0d0*dx2**2+3.0d0*dx1*(28.0d0*dx2+5.0d0*dx1))
      ci(4)=2.0d0*dx1**11*dx2**3*(182.0d0*dx2**3+dx1*(273.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**10*dx2**2*(1001.0d0*dx2**4+dx1*(2184.0d0*dx2**3+5.0d0*dx1*(273.0d0*dx2**2+dx1*(56.0d0*dx2+3.0d0*dx1))))
      ci(6)=2.0d0*dx1**9*dx2*(1001.0d0*dx2**5+dx1*(3003.0d0*dx2**4+dx1*(2730.0d0*dx2**3+dx1*(910.0d0*dx2**2+3.0d0*dx1*(35.0d0*dx2 &
+dx1)))))
      ci(7)=dx1**8*(3003.0d0*dx2**6+dx1*(12012.0d0*dx2**5+dx1*(15015.0d0*dx2**4+dx1*(7280.0d0*dx2**3+dx1*(1365.0d0*dx2**2+dx1* &
(84.0d0*dx2+dx1))))))
      ci(8)=2.0d0*dx1**7*(1716.0d0*dx2**6+7.0d0*dx1*(1287.0d0*dx2**5+dx1*(2145.0d0*dx2**4+dx1*(1430.0d0*dx2**3+dx1*(390.0d0*dx2**2 &
+dx1*(39.0d0*dx2+dx1))))))
      ci(9)=13.0d0*dx1**6*(231.0d0*dx2**6+dx1*(1584.0d0*dx2**5+7.0d0*dx1*(495.0d0*dx2**4+dx1*(440.0d0*dx2**3+dx1*(165.0d0*dx2**2 &
+dx1*(24.0d0*dx2+dx1))))))
      ci(10)=26.0d0*dx1**5*(77.0d0*dx2**6+dx1*(693.0d0*dx2**5+dx1*(1980.0d0*dx2**4+7.0d0*dx1*(330.0d0*dx2**3+dx1*(165.0d0*dx2**2 &
+dx1*(33.0d0*dx2+2.0d0*dx1))))))
      ci(11)=143.0d0*dx1**4*(7.0d0*dx2**6+dx1*(84.0d0*dx2**5+dx1*(315.0d0*dx2**4+dx1*(480.0d0*dx2**3+7.0d0*dx1*(45.0d0*dx2**2+dx1* &
(12.0d0*dx2+dx1))))))
      ci(12)=26.0d0*dx1**3*(14.0d0*dx2**6+11.0d0*dx1*(21.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(210.0d0*dx2**3+dx1*(180.0d0*dx2**2 &
+7.0d0*dx1*(9.0d0*dx2+dx1))))))
      ci(13)=13.0d0*dx1**2*(7.0d0*dx2**6+dx1*(168.0d0*dx2**5+11.0d0*dx1*(105.0d0*dx2**4+dx1*(280.0d0*dx2**3+3.0d0*dx1* &
(105.0d0*dx2**2+dx1*(48.0d0*dx2+7.0d0*dx1))))))
      ci(14)=2.0d0*dx1*(7.0d0*dx2**6+13.0d0*dx1*(21.0d0*dx2**5+dx1*(210.0d0*dx2**4+11.0d0*dx1*(70.0d0*dx2**3+3.0d0*dx1* &
(35.0d0*dx2**2+dx1*(21.0d0*dx2+4.0d0*dx1))))))
      ci(15)=dx2**6+7.0d0*dx1*(12.0d0*dx2**5+13.0d0*dx1*(15.0d0*dx2**4+dx1*(80.0d0*dx2**3+33.0d0*dx1*(5.0d0*dx2**2+dx1*(4.0d0*dx2 &
+dx1)))))
      ci(16)=2.0d0*(3.0d0*dx2**5+7.0d0*dx1*(15.0d0*dx2**4+13.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+dx1)) &
)))
      ci(17)=15.0d0*dx2**4+7.0d0*dx1*(40.0d0*dx2**3+13.0d0*dx1*(15.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1)))
      ci(18)=2.0d0*(10.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(19)=15.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+13.0d0*dx1)
      ci(20)=2.0d0*(3.0d0*dx2+7.0d0*dx1)
      ci(21)=1.0d0
    case(7)
      ci(1)=dx1**14*dx2**7
      ci(2)=7.0d0*dx1**13*dx2**6*(2.0d0*dx2+dx1)
      ci(3)=7.0d0*dx1**12*dx2**5*(13.0d0*dx2**2+dx1*(14.0d0*dx2+3.0d0*dx1))
      ci(4)=7.0d0*dx1**11*dx2**4*(52.0d0*dx2**3+dx1*(91.0d0*dx2**2+dx1*(42.0d0*dx2+5.0d0*dx1)))
      ci(5)=7.0d0*dx1**10*dx2**3*(143.0d0*dx2**4+dx1*(364.0d0*dx2**3+dx1*(273.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**9*dx2**2*(286.0d0*dx2**5+dx1*(1001.0d0*dx2**4+dx1*(1092.0d0*dx2**3+dx1*(455.0d0*dx2**2+dx1*(70.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=7.0d0*dx1**8*dx2*(429.0d0*dx2**6+dx1*(2002.0d0*dx2**5+dx1*(3003.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1*(455.0d0*dx2**2 &
+dx1*(42.0d0*dx2+dx1))))))
      ci(8)=dx1**7*(3432.0d0*dx2**7+dx1*(21021.0d0*dx2**6+dx1*(42042.0d0*dx2**5+dx1*(35035.0d0*dx2**4+dx1*(12740.0d0*dx2**3+dx1* &
(1911.0d0*dx2**2+dx1*(98.0d0*dx2+dx1)))))))
      ci(9)=7.0d0*dx1**6*(429.0d0*dx2**7+dx1*(3432.0d0*dx2**6+dx1*(9009.0d0*dx2**5+dx1*(10010.0d0*dx2**4+dx1*(5005.0d0*dx2**3+dx1* &
(1092.0d0*dx2**2+dx1*(91.0d0*dx2+2.0d0*dx1)))))))
      ci(10)=91.0d0*dx1**5*(22.0d0*dx2**7+dx1*(231.0d0*dx2**6+dx1*(792.0d0*dx2**5+dx1*(1155.0d0*dx2**4+dx1*(770.0d0*dx2**3+dx1* &
(231.0d0*dx2**2+dx1*(28.0d0*dx2+dx1)))))))
      ci(11)=91.0d0*dx1**4*(11.0d0*dx2**7+dx1*(154.0d0*dx2**6+dx1*(693.0d0*dx2**5+dx1*(1320.0d0*dx2**4+dx1*(1155.0d0*dx2**3+dx1* &
(462.0d0*dx2**2+dx1*(77.0d0*dx2+4.0d0*dx1)))))))
      ci(12)=91.0d0*dx1**3*(4.0d0*dx2**7+11.0d0*dx1*(7.0d0*dx2**6+dx1*(42.0d0*dx2**5+dx1*(105.0d0*dx2**4+dx1*(120.0d0*dx2**3+dx1* &
(63.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)))))))
      ci(13)=91.0d0*dx1**2*(dx2**7+dx1*(28.0d0*dx2**6+11.0d0*dx1*(21.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1*(105.0d0*dx2**3+dx1* &
(72.0d0*dx2**2+dx1*(21.0d0*dx2+2.0d0*dx1)))))))
      ci(14)=7.0d0*dx1*(2.0d0*dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(84.0d0*dx2**5+11.0d0*dx1*(35.0d0*dx2**4+dx1*(70.0d0*dx2**3 &
+3.0d0*dx1*(21.0d0*dx2**2+dx1*(8.0d0*dx2+dx1)))))))
      ci(15)=dx2**7+dx1*(98.0d0*dx2**6+13.0d0*dx1*(147.0d0*dx2**5+dx1*(980.0d0*dx2**4+11.0d0*dx1*(245.0d0*dx2**3+3.0d0*dx1* &
(98.0d0*dx2**2+dx1*(49.0d0*dx2+8.0d0*dx1))))))
      ci(16)=7.0d0*(dx2**6+dx1*(42.0d0*dx2**5+13.0d0*dx1*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+11.0d0*dx1*(21.0d0*dx2**2+dx1* &
(14.0d0*dx2+3.0d0*dx1))))))
      ci(17)=7.0d0*(3.0d0*dx2**5+dx1*(70.0d0*dx2**4+13.0d0*dx1*(35.0d0*dx2**3+dx1*(84.0d0*dx2**2+11.0d0*dx1*(7.0d0*dx2+2.0d0*dx1)) &
)))
      ci(18)=7.0d0*(5.0d0*dx2**4+dx1*(70.0d0*dx2**3+13.0d0*dx1*(21.0d0*dx2**2+dx1*(28.0d0*dx2+11.0d0*dx1))))
      ci(19)=7.0d0*(5.0d0*dx2**3+dx1*(42.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+4.0d0*dx1)))
      ci(20)=7.0d0*(3.0d0*dx2**2+dx1*(14.0d0*dx2+13.0d0*dx1))
      ci(21)=7.0d0*(dx2+2.0d0*dx1)
      ci(22)=1.0d0
    case(8)
      ci(1)=dx1**14*dx2**8
      ci(2)=2.0d0*dx1**13*dx2**7*(7.0d0*dx2+4.0d0*dx1)
      ci(3)=7.0d0*dx1**12*dx2**6*(13.0d0*dx2**2+4.0d0*dx1*(4.0d0*dx2+dx1))
      ci(4)=28.0d0*dx1**11*dx2**5*(13.0d0*dx2**3+2.0d0*dx1*(13.0d0*dx2**2+dx1*(7.0d0*dx2+dx1)))
      ci(5)=7.0d0*dx1**10*dx2**4*(143.0d0*dx2**4+2.0d0*dx1*(208.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(56.0d0*dx2+5.0d0*dx1))))
      ci(6)=14.0d0*dx1**9*dx2**3*(143.0d0*dx2**5+2.0d0*dx1*(286.0d0*dx2**4+dx1*(364.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(35.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=7.0d0*dx1**8*dx2**2*(429.0d0*dx2**6+2.0d0*dx1*(1144.0d0*dx2**5+dx1*(2002.0d0*dx2**4+dx1*(1456.0d0*dx2**3+dx1* &
(455.0d0*dx2**2+2.0d0*dx1*(28.0d0*dx2+dx1))))))
      ci(8)=8.0d0*dx1**7*dx2*(429.0d0*dx2**7+dx1*(3003.0d0*dx2**6+dx1*(7007.0d0*dx2**5+dx1*(7007.0d0*dx2**4+dx1*(3185.0d0*dx2**3 &
+dx1*(637.0d0*dx2**2+dx1*(49.0d0*dx2+dx1)))))))
      ci(9)=dx1**6*(3003.0d0*dx2**8+dx1*(27456.0d0*dx2**7+dx1*(84084.0d0*dx2**6+dx1*(112112.0d0*dx2**5+dx1*(70070.0d0*dx2**4+dx1* &
(20384.0d0*dx2**3+dx1*(2548.0d0*dx2**2+dx1*(112.0d0*dx2+dx1))))))))
      ci(10)=14.0d0*dx1**5*(143.0d0*dx2**8+dx1*(1716.0d0*dx2**7+dx1*(6864.0d0*dx2**6+dx1*(12012.0d0*dx2**5+dx1*(10010.0d0*dx2**4 &
+dx1*(4004.0d0*dx2**3+dx1*(728.0d0*dx2**2+dx1*(52.0d0*dx2+dx1))))))))
      ci(11)=91.0d0*dx1**4*(11.0d0*dx2**8+dx1*(176.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2112.0d0*dx2**5+dx1*(2310.0d0*dx2**4+dx1* &
(1232.0d0*dx2**3+dx1*(308.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))))))
      ci(12)=364.0d0*dx1**3*(dx2**8+dx1*(22.0d0*dx2**7+dx1*(154.0d0*dx2**6+dx1*(462.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1* &
(462.0d0*dx2**3+dx1*(154.0d0*dx2**2+dx1*(22.0d0*dx2+dx1))))))))
      ci(13)=91.0d0*dx1**2*(dx2**8+dx1*(32.0d0*dx2**7+11.0d0*dx1*(28.0d0*dx2**6+dx1*(112.0d0*dx2**5+dx1*(210.0d0*dx2**4+dx1* &
(192.0d0*dx2**3+dx1*(84.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))))))
      ci(14)=14.0d0*dx1*(dx2**8+13.0d0*dx1*(4.0d0*dx2**7+dx1*(56.0d0*dx2**6+11.0d0*dx1*(28.0d0*dx2**5+dx1*(70.0d0*dx2**4+dx1* &
(84.0d0*dx2**3+dx1*(48.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))))))
      ci(15)=dx2**8+dx1*(112.0d0*dx2**7+13.0d0*dx1*(196.0d0*dx2**6+dx1*(1568.0d0*dx2**5+11.0d0*dx1*(490.0d0*dx2**4+dx1* &
(784.0d0*dx2**3+3.0d0*dx1*(196.0d0*dx2**2+dx1*(64.0d0*dx2+7.0d0*dx1)))))))
      ci(16)=8.0d0*(dx2**7+dx1*(49.0d0*dx2**6+13.0d0*dx1*(49.0d0*dx2**5+dx1*(245.0d0*dx2**4+11.0d0*dx1*(49.0d0*dx2**3+dx1* &
(49.0d0*dx2**2+3.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(17)=7.0d0*(4.0d0*dx2**6+dx1*(112.0d0*dx2**5+13.0d0*dx1*(70.0d0*dx2**4+dx1*(224.0d0*dx2**3+11.0d0*dx1*(28.0d0*dx2**2+dx1* &
(16.0d0*dx2+3.0d0*dx1))))))
      ci(18)=14.0d0*(4.0d0*dx2**5+dx1*(70.0d0*dx2**4+13.0d0*dx1*(28.0d0*dx2**3+dx1*(56.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+dx1)))))
      ci(19)=7.0d0*(10.0d0*dx2**4+dx1*(112.0d0*dx2**3+13.0d0*dx1*(28.0d0*dx2**2+dx1*(32.0d0*dx2+11.0d0*dx1))))
      ci(20)=28.0d0*(2.0d0*dx2**3+dx1*(14.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+dx1)))
      ci(21)=7.0d0*(4.0d0*dx2**2+dx1*(16.0d0*dx2+13.0d0*dx1))
      ci(22)=2.0d0*(4.0d0*dx2+7.0d0*dx1)
      ci(23)=1.0d0
    case(9)
      ci(1)=dx1**14*dx2**9
      ci(2)=dx1**13*dx2**8*(14.0d0*dx2+9.0d0*dx1)
      ci(3)=dx1**12*dx2**7*(91.0d0*dx2**2+18.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))
      ci(4)=7.0d0*dx1**11*dx2**6*(52.0d0*dx2**3+3.0d0*dx1*(39.0d0*dx2**2+4.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(5)=7.0d0*dx1**10*dx2**5*(143.0d0*dx2**4+6.0d0*dx1*(78.0d0*dx2**3+dx1*(78.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1))))
      ci(6)=7.0d0*dx1**9*dx2**4*(286.0d0*dx2**5+3.0d0*dx1*(429.0d0*dx2**4+2.0d0*dx1*(312.0d0*dx2**3+dx1*(182.0d0*dx2**2+3.0d0*dx1* &
(14.0d0*dx2+dx1)))))
      ci(7)=21.0d0*dx1**8*dx2**3*(143.0d0*dx2**6+2.0d0*dx1*(429.0d0*dx2**5+dx1*(858.0d0*dx2**4+dx1*(728.0d0*dx2**3+dx1* &
(273.0d0*dx2**2+2.0d0*dx1*(21.0d0*dx2+dx1))))))
      ci(8)=3.0d0*dx1**7*dx2**2*(1144.0d0*dx2**7+dx1*(9009.0d0*dx2**6+2.0d0*dx1*(12012.0d0*dx2**5+dx1*(14014.0d0*dx2**4+dx1* &
(7644.0d0*dx2**3+dx1*(1911.0d0*dx2**2+2.0d0*dx1*(98.0d0*dx2+3.0d0*dx1)))))))
      ci(9)=3.0d0*dx1**6*dx2*(1001.0d0*dx2**8+dx1*(10296.0d0*dx2**7+dx1*(36036.0d0*dx2**6+dx1*(56056.0d0*dx2**5+dx1* &
(42042.0d0*dx2**4+dx1*(15288.0d0*dx2**3+dx1*(2548.0d0*dx2**2+3.0d0*dx1*(56.0d0*dx2+dx1))))))))
      ci(10)=dx1**5*(2002.0d0*dx2**9+dx1*(27027.0d0*dx2**8+dx1*(123552.0d0*dx2**7+dx1*(252252.0d0*dx2**6+dx1*(252252.0d0*dx2**5 &
+dx1*(126126.0d0*dx2**4+dx1*(30576.0d0*dx2**3+dx1*(3276.0d0*dx2**2+dx1*(126.0d0*dx2+dx1)))))))))
      ci(11)=7.0d0*dx1**4*(143.0d0*dx2**9+dx1*(2574.0d0*dx2**8+dx1*(15444.0d0*dx2**7+dx1*(41184.0d0*dx2**6+dx1*(54054.0d0*dx2**5 &
+dx1*(36036.0d0*dx2**4+dx1*(12012.0d0*dx2**3+dx1*(1872.0d0*dx2**2+dx1*(117.0d0*dx2+2.0d0*dx1)))))))))
      ci(12)=91.0d0*dx1**3*(4.0d0*dx2**9+dx1*(99.0d0*dx2**8+dx1*(792.0d0*dx2**7+dx1*(2772.0d0*dx2**6+dx1*(4752.0d0*dx2**5+dx1* &
(4158.0d0*dx2**4+dx1*(1848.0d0*dx2**3+dx1*(396.0d0*dx2**2+dx1*(36.0d0*dx2+dx1)))))))))
      ci(13)=91.0d0*dx1**2*(dx2**9+dx1*(36.0d0*dx2**8+dx1*(396.0d0*dx2**7+dx1*(1848.0d0*dx2**6+dx1*(4158.0d0*dx2**5+dx1* &
(4752.0d0*dx2**4+dx1*(2772.0d0*dx2**3+dx1*(792.0d0*dx2**2+dx1*(99.0d0*dx2+4.0d0*dx1)))))))))
      ci(14)=7.0d0*dx1*(2.0d0*dx2**9+13.0d0*dx1*(9.0d0*dx2**8+dx1*(144.0d0*dx2**7+11.0d0*dx1*(84.0d0*dx2**6+dx1*(252.0d0*dx2**5 &
+dx1*(378.0d0*dx2**4+dx1*(288.0d0*dx2**3+dx1*(108.0d0*dx2**2+dx1*(18.0d0*dx2+dx1)))))))))
      ci(15)=dx2**9+dx1*(126.0d0*dx2**8+13.0d0*dx1*(252.0d0*dx2**7+dx1*(2352.0d0*dx2**6+11.0d0*dx1*(882.0d0*dx2**5+dx1* &
(1764.0d0*dx2**4+dx1*(1764.0d0*dx2**3+dx1*(864.0d0*dx2**2+7.0d0*dx1*(27.0d0*dx2+2.0d0*dx1))))))))
      ci(16)=3.0d0*(3.0d0*dx2**8+dx1*(168.0d0*dx2**7+13.0d0*dx1*(196.0d0*dx2**6+dx1*(1176.0d0*dx2**5+11.0d0*dx1*(294.0d0*dx2**4 &
+dx1*(392.0d0*dx2**3+dx1*(252.0d0*dx2**2+dx1*(72.0d0*dx2+7.0d0*dx1))))))))
      ci(17)=3.0d0*(12.0d0*dx2**7+dx1*(392.0d0*dx2**6+13.0d0*dx1*(294.0d0*dx2**5+dx1*(1176.0d0*dx2**4+11.0d0*dx1*(196.0d0*dx2**3 &
+dx1*(168.0d0*dx2**2+dx1*(63.0d0*dx2+8.0d0*dx1)))))))
      ci(18)=21.0d0*(4.0d0*dx2**6+dx1*(84.0d0*dx2**5+13.0d0*dx1*(42.0d0*dx2**4+dx1*(112.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+dx1* &
(6.0d0*dx2+dx1))))))
      ci(19)=7.0d0*(18.0d0*dx2**5+dx1*(252.0d0*dx2**4+13.0d0*dx1*(84.0d0*dx2**3+dx1*(144.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2 &
+2.0d0*dx1)))))
      ci(20)=7.0d0*(18.0d0*dx2**4+dx1*(168.0d0*dx2**3+13.0d0*dx1*(36.0d0*dx2**2+dx1*(36.0d0*dx2+11.0d0*dx1))))
      ci(21)=7.0d0*(12.0d0*dx2**3+dx1*(72.0d0*dx2**2+13.0d0*dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(22)=36.0d0*dx2**2+7.0d0*dx1*(18.0d0*dx2+13.0d0*dx1)
      ci(23)=9.0d0*dx2+14.0d0*dx1
      ci(24)=1.0d0
    case(10)
      ci(1)=dx1**14*dx2**10
      ci(2)=2.0d0*dx1**13*dx2**9*(7.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**12*dx2**8*(91.0d0*dx2**2+5.0d0*dx1*(28.0d0*dx2+9.0d0*dx1))
      ci(4)=2.0d0*dx1**11*dx2**7*(182.0d0*dx2**3+5.0d0*dx1*(91.0d0*dx2**2+3.0d0*dx1*(21.0d0*dx2+4.0d0*dx1)))
      ci(5)=7.0d0*dx1**10*dx2**6*(143.0d0*dx2**4+5.0d0*dx1*(104.0d0*dx2**3+3.0d0*dx1*(39.0d0*dx2**2+2.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(6)=14.0d0*dx1**9*dx2**5*(143.0d0*dx2**5+dx1*(715.0d0*dx2**4+6.0d0*dx1*(195.0d0*dx2**3+dx1*(130.0d0*dx2**2+dx1*(35.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=7.0d0*dx1**8*dx2**4*(429.0d0*dx2**6+dx1*(2860.0d0*dx2**5+3.0d0*dx1*(2145.0d0*dx2**4+2.0d0*dx1*(1040.0d0*dx2**3+dx1* &
(455.0d0*dx2**2+dx1*(84.0d0*dx2+5.0d0*dx1))))))
      ci(8)=6.0d0*dx1**7*dx2**3*(572.0d0*dx2**7+dx1*(5005.0d0*dx2**6+dx1*(15015.0d0*dx2**5+2.0d0*dx1*(10010.0d0*dx2**4+dx1* &
(6370.0d0*dx2**3+dx1*(1911.0d0*dx2**2+5.0d0*dx1*(49.0d0*dx2+2.0d0*dx1)))))))
      ci(9)=3.0d0*dx1**6*dx2**2*(1001.0d0*dx2**8+dx1*(11440.0d0*dx2**7+dx1*(45045.0d0*dx2**6+dx1*(80080.0d0*dx2**5+dx1* &
(70070.0d0*dx2**4+dx1*(30576.0d0*dx2**3+5.0d0*dx1*(1274.0d0*dx2**2+dx1*(112.0d0*dx2+3.0d0*dx1))))))))
      ci(10)=2.0d0*dx1**5*dx2*(1001.0d0*dx2**9+dx1*(15015.0d0*dx2**8+dx1*(77220.0d0*dx2**7+dx1*(180180.0d0*dx2**6+dx1* &
(210210.0d0*dx2**5+dx1*(126126.0d0*dx2**4+5.0d0*dx1*(7644.0d0*dx2**3+dx1*(1092.0d0*dx2**2+dx1*(63.0d0*dx2+dx1)))))))))
      ci(11)=dx1**4*(1001.0d0*dx2**10+dx1*(20020.0d0*dx2**9+dx1*(135135.0d0*dx2**8+dx1*(411840.0d0*dx2**7+dx1*(630630.0d0*dx2**6 &
+dx1*(504504.0d0*dx2**5+dx1*(210210.0d0*dx2**4+dx1*(43680.0d0*dx2**3+dx1*(4095.0d0*dx2**2+dx1*(140.0d0*dx2+dx1))))))))))
      ci(12)=14.0d0*dx1**3*(26.0d0*dx2**10+dx1*(715.0d0*dx2**9+dx1*(6435.0d0*dx2**8+dx1*(25740.0d0*dx2**7+dx1*(51480.0d0*dx2**6 &
+dx1*(54054.0d0*dx2**5+dx1*(30030.0d0*dx2**4+dx1*(8580.0d0*dx2**3+dx1*(1170.0d0*dx2**2+dx1*(65.0d0*dx2+dx1))))))))))
      ci(13)=91.0d0*dx1**2*(dx2**10+dx1*(40.0d0*dx2**9+dx1*(495.0d0*dx2**8+dx1*(2640.0d0*dx2**7+dx1*(6930.0d0*dx2**6+dx1* &
(9504.0d0*dx2**5+dx1*(6930.0d0*dx2**4+dx1*(2640.0d0*dx2**3+dx1*(495.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))))))))))
      ci(14)=14.0d0*dx1*(dx2**10+13.0d0*dx1*(5.0d0*dx2**9+dx1*(90.0d0*dx2**8+dx1*(660.0d0*dx2**7+dx1*(2310.0d0*dx2**6+dx1* &
(4158.0d0*dx2**5+dx1*(3960.0d0*dx2**4+dx1*(1980.0d0*dx2**3+dx1*(495.0d0*dx2**2+dx1*(55.0d0*dx2+2.0d0*dx1))))))))))
      ci(15)=dx2**10+dx1*(140.0d0*dx2**9+13.0d0*dx1*(315.0d0*dx2**8+dx1*(3360.0d0*dx2**7+11.0d0*dx1*(1470.0d0*dx2**6+dx1* &
(3528.0d0*dx2**5+dx1*(4410.0d0*dx2**4+dx1*(2880.0d0*dx2**3+7.0d0*dx1*(135.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)))))))))
      ci(16)=2.0d0*(5.0d0*dx2**9+dx1*(315.0d0*dx2**8+13.0d0*dx1*(420.0d0*dx2**7+dx1*(2940.0d0*dx2**6+11.0d0*dx1*(882.0d0*dx2**5 &
+dx1*(1470.0d0*dx2**4+dx1*(1260.0d0*dx2**3+dx1*(540.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+dx1)))))))))
      ci(17)=3.0d0*(15.0d0*dx2**8+dx1*(560.0d0*dx2**7+13.0d0*dx1*(490.0d0*dx2**6+dx1*(2352.0d0*dx2**5+11.0d0*dx1*(490.0d0*dx2**4 &
+dx1*(560.0d0*dx2**3+dx1*(315.0d0*dx2**2+dx1*(80.0d0*dx2+7.0d0*dx1))))))))
      ci(18)=6.0d0*(20.0d0*dx2**7+dx1*(490.0d0*dx2**6+13.0d0*dx1*(294.0d0*dx2**5+dx1*(980.0d0*dx2**4+11.0d0*dx1*(140.0d0*dx2**3 &
+dx1*(105.0d0*dx2**2+dx1*(35.0d0*dx2+4.0d0*dx1)))))))
      ci(19)=7.0d0*(30.0d0*dx2**6+dx1*(504.0d0*dx2**5+13.0d0*dx1*(210.0d0*dx2**4+dx1*(480.0d0*dx2**3+11.0d0*dx1*(45.0d0*dx2**2 &
+dx1*(20.0d0*dx2+3.0d0*dx1))))))
      ci(20)=14.0d0*(18.0d0*dx2**5+dx1*(210.0d0*dx2**4+13.0d0*dx1*(60.0d0*dx2**3+dx1*(90.0d0*dx2**2+11.0d0*dx1*(5.0d0*dx2+dx1)))))
      ci(21)=7.0d0*(30.0d0*dx2**4+dx1*(240.0d0*dx2**3+13.0d0*dx1*(45.0d0*dx2**2+dx1*(40.0d0*dx2+11.0d0*dx1))))
      ci(22)=2.0d0*(60.0d0*dx2**3+7.0d0*dx1*(45.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)))
      ci(23)=45.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+13.0d0*dx1)
      ci(24)=2.0d0*(5.0d0*dx2+7.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>10, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^15 * (x-x2)^n2 as sum_k=0^(15+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_15(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**15
      ci(2)=15.0d0*dx1**14
      ci(3)=105.0d0*dx1**13
      ci(4)=455.0d0*dx1**12
      ci(5)=1365.0d0*dx1**11
      ci(6)=3003.0d0*dx1**10
      ci(7)=5005.0d0*dx1**9
      ci(8)=6435.0d0*dx1**8
      ci(9)=6435.0d0*dx1**7
      ci(10)=5005.0d0*dx1**6
      ci(11)=3003.0d0*dx1**5
      ci(12)=1365.0d0*dx1**4
      ci(13)=455.0d0*dx1**3
      ci(14)=105.0d0*dx1**2
      ci(15)=15.0d0*dx1
      ci(16)=1.0d0
    case(1)
      ci(1)=dx1**15*dx2
      ci(2)=dx1**14*(15.0d0*dx2+dx1)
      ci(3)=15.0d0*dx1**13*(7.0d0*dx2+dx1)
      ci(4)=35.0d0*dx1**12*(13.0d0*dx2+3.0d0*dx1)
      ci(5)=455.0d0*dx1**11*(3.0d0*dx2+dx1)
      ci(6)=273.0d0*dx1**10*(11.0d0*dx2+5.0d0*dx1)
      ci(7)=1001.0d0*dx1**9*(5.0d0*dx2+3.0d0*dx1)
      ci(8)=715.0d0*dx1**8*(9.0d0*dx2+7.0d0*dx1)
      ci(9)=6435.0d0*dx1**7*(dx2+dx1)
      ci(10)=715.0d0*dx1**6*(7.0d0*dx2+9.0d0*dx1)
      ci(11)=1001.0d0*dx1**5*(3.0d0*dx2+5.0d0*dx1)
      ci(12)=273.0d0*dx1**4*(5.0d0*dx2+11.0d0*dx1)
      ci(13)=455.0d0*dx1**3*(dx2+3.0d0*dx1)
      ci(14)=35.0d0*dx1**2*(3.0d0*dx2+13.0d0*dx1)
      ci(15)=15.0d0*dx1*(dx2+7.0d0*dx1)
      ci(16)=dx2+15.0d0*dx1
      ci(17)=1.0d0
    case(2)
      ci(1)=dx1**15*dx2**2
      ci(2)=dx1**14*dx2*(15.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**13*(105.0d0*dx2**2+dx1*(30.0d0*dx2+dx1))
      ci(4)=5.0d0*dx1**12*(91.0d0*dx2**2+3.0d0*dx1*(14.0d0*dx2+dx1))
      ci(5)=35.0d0*dx1**11*(39.0d0*dx2**2+dx1*(26.0d0*dx2+3.0d0*dx1))
      ci(6)=91.0d0*dx1**10*(33.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+dx1))
      ci(7)=91.0d0*dx1**9*(55.0d0*dx2**2+3.0d0*dx1*(22.0d0*dx2+5.0d0*dx1))
      ci(8)=143.0d0*dx1**8*(45.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+3.0d0*dx1))
      ci(9)=715.0d0*dx1**7*(9.0d0*dx2**2+dx1*(18.0d0*dx2+7.0d0*dx1))
      ci(10)=715.0d0*dx1**6*(7.0d0*dx2**2+9.0d0*dx1*(2.0d0*dx2+dx1))
      ci(11)=143.0d0*dx1**5*(21.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+9.0d0*dx1))
      ci(12)=91.0d0*dx1**4*(15.0d0*dx2**2+11.0d0*dx1*(6.0d0*dx2+5.0d0*dx1))
      ci(13)=91.0d0*dx1**3*(5.0d0*dx2**2+3.0d0*dx1*(10.0d0*dx2+11.0d0*dx1))
      ci(14)=35.0d0*dx1**2*(3.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(15)=5.0d0*dx1*(3.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+13.0d0*dx1))
      ci(16)=dx2**2+15.0d0*dx1*(2.0d0*dx2+7.0d0*dx1)
      ci(17)=2.0d0*dx2+15.0d0*dx1
      ci(18)=1.0d0
    case(3)
      ci(1)=dx1**15*dx2**3
      ci(2)=3.0d0*dx1**14*dx2**2*(5.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**13*dx2*(35.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))
      ci(4)=dx1**12*(455.0d0*dx2**3+dx1*(315.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))
      ci(5)=15.0d0*dx1**11*(91.0d0*dx2**3+dx1*(91.0d0*dx2**2+dx1*(21.0d0*dx2+dx1)))
      ci(6)=21.0d0*dx1**10*(143.0d0*dx2**3+5.0d0*dx1*(39.0d0*dx2**2+dx1*(13.0d0*dx2+dx1)))
      ci(7)=91.0d0*dx1**9*(55.0d0*dx2**3+dx1*(99.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+dx1)))
      ci(8)=39.0d0*dx1**8*(165.0d0*dx2**3+7.0d0*dx1*(55.0d0*dx2**2+dx1*(33.0d0*dx2+5.0d0*dx1)))
      ci(9)=429.0d0*dx1**7*(15.0d0*dx2**3+dx1*(45.0d0*dx2**2+7.0d0*dx1*(5.0d0*dx2+dx1)))
      ci(10)=715.0d0*dx1**6*(7.0d0*dx2**3+dx1*(27.0d0*dx2**2+dx1*(27.0d0*dx2+7.0d0*dx1)))
      ci(11)=429.0d0*dx1**5*(7.0d0*dx2**3+5.0d0*dx1*(7.0d0*dx2**2+3.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(12)=39.0d0*dx1**4*(35.0d0*dx2**3+11.0d0*dx1*(21.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+3.0d0*dx1)))
      ci(13)=91.0d0*dx1**3*(5.0d0*dx2**3+dx1*(45.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+5.0d0*dx1)))
      ci(14)=21.0d0*dx1**2*(5.0d0*dx2**3+13.0d0*dx1*(5.0d0*dx2**2+dx1*(15.0d0*dx2+11.0d0*dx1)))
      ci(15)=15.0d0*dx1*(dx2**3+7.0d0*dx1*(3.0d0*dx2**2+13.0d0*dx1*(dx2+dx1)))
      ci(16)=dx2**3+5.0d0*dx1*(9.0d0*dx2**2+7.0d0*dx1*(9.0d0*dx2+13.0d0*dx1))
      ci(17)=3.0d0*(dx2**2+5.0d0*dx1*(3.0d0*dx2+7.0d0*dx1))
      ci(18)=3.0d0*(dx2+5.0d0*dx1)
      ci(19)=1.0d0
    case(4)
      ci(1)=dx1**15*dx2**4
      ci(2)=dx1**14*dx2**3*(15.0d0*dx2+4.0d0*dx1)
      ci(3)=3.0d0*dx1**13*dx2**2*(35.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+dx1))
      ci(4)=dx1**12*dx2*(455.0d0*dx2**3+2.0d0*dx1*(210.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**11*(1365.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1*(630.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))
      ci(6)=3.0d0*dx1**10*(1001.0d0*dx2**4+5.0d0*dx1*(364.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(28.0d0*dx2+dx1))))
      ci(7)=7.0d0*dx1**9*(715.0d0*dx2**4+dx1*(1716.0d0*dx2**3+5.0d0*dx1*(234.0d0*dx2**2+dx1*(52.0d0*dx2+3.0d0*dx1))))
      ci(8)=13.0d0*dx1**8*(495.0d0*dx2**4+7.0d0*dx1*(220.0d0*dx2**3+dx1*(198.0d0*dx2**2+5.0d0*dx1*(12.0d0*dx2+dx1))))
      ci(9)=39.0d0*dx1**7*(165.0d0*dx2**4+dx1*(660.0d0*dx2**3+7.0d0*dx1*(110.0d0*dx2**2+dx1*(44.0d0*dx2+5.0d0*dx1))))
      ci(10)=143.0d0*dx1**6*(35.0d0*dx2**4+dx1*(180.0d0*dx2**3+dx1*(270.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+3.0d0*dx1))))
      ci(11)=143.0d0*dx1**5*(21.0d0*dx2**4+5.0d0*dx1*(28.0d0*dx2**3+dx1*(54.0d0*dx2**2+dx1*(36.0d0*dx2+7.0d0*dx1))))
      ci(12)=39.0d0*dx1**4*(35.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+5.0d0*dx1*(14.0d0*dx2**2+3.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(13)=13.0d0*dx1**3*(35.0d0*dx2**4+dx1*(420.0d0*dx2**3+11.0d0*dx1*(126.0d0*dx2**2+5.0d0*dx1*(28.0d0*dx2+9.0d0*dx1))))
      ci(14)=7.0d0*dx1**2*(15.0d0*dx2**4+13.0d0*dx1*(20.0d0*dx2**3+dx1*(90.0d0*dx2**2+11.0d0*dx1*(12.0d0*dx2+5.0d0*dx1))))
      ci(15)=3.0d0*dx1*(5.0d0*dx2**4+7.0d0*dx1*(20.0d0*dx2**3+13.0d0*dx1*(10.0d0*dx2**2+dx1*(20.0d0*dx2+11.0d0*dx1))))
      ci(16)=dx2**4+5.0d0*dx1*(12.0d0*dx2**3+7.0d0*dx1*(18.0d0*dx2**2+13.0d0*dx1*(4.0d0*dx2+3.0d0*dx1)))
      ci(17)=4.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+13.0d0*dx1))
      ci(18)=3.0d0*(2.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+7.0d0*dx1))
      ci(19)=4.0d0*dx2+15.0d0*dx1
      ci(20)=1.0d0
    case(5)
      ci(1)=dx1**15*dx2**5
      ci(2)=5.0d0*dx1**14*dx2**4*(3.0d0*dx2+dx1)
      ci(3)=5.0d0*dx1**13*dx2**3*(21.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1))
      ci(4)=5.0d0*dx1**12*dx2**2*(91.0d0*dx2**3+dx1*(105.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**11*dx2*(273.0d0*dx2**4+dx1*(455.0d0*dx2**3+dx1*(210.0d0*dx2**2+dx1*(30.0d0*dx2+dx1))))
      ci(6)=dx1**10*(3003.0d0*dx2**5+dx1*(6825.0d0*dx2**4+dx1*(4550.0d0*dx2**3+dx1*(1050.0d0*dx2**2+dx1*(75.0d0*dx2+dx1)))))
      ci(7)=5.0d0*dx1**9*(1001.0d0*dx2**5+dx1*(3003.0d0*dx2**4+dx1*(2730.0d0*dx2**3+dx1*(910.0d0*dx2**2+3.0d0*dx1*(35.0d0*dx2+dx1) &
))))
      ci(8)=5.0d0*dx1**8*(1287.0d0*dx2**5+7.0d0*dx1*(715.0d0*dx2**4+dx1*(858.0d0*dx2**3+dx1*(390.0d0*dx2**2+dx1*(65.0d0*dx2 &
+3.0d0*dx1)))))
      ci(9)=65.0d0*dx1**7*(99.0d0*dx2**5+dx1*(495.0d0*dx2**4+7.0d0*dx1*(110.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(15.0d0*dx2+dx1)))))
      ci(10)=65.0d0*dx1**6*(77.0d0*dx2**5+dx1*(495.0d0*dx2**4+dx1*(990.0d0*dx2**3+7.0d0*dx1*(110.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2 &
+dx1)))))
      ci(11)=143.0d0*dx1**5*(21.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(450.0d0*dx2**3+dx1*(450.0d0*dx2**2+7.0d0*dx1*(25.0d0*dx2 &
+3.0d0*dx1)))))
      ci(12)=65.0d0*dx1**4*(21.0d0*dx2**5+11.0d0*dx1*(21.0d0*dx2**4+dx1*(70.0d0*dx2**3+dx1*(90.0d0*dx2**2+dx1*(45.0d0*dx2 &
+7.0d0*dx1)))))
      ci(13)=65.0d0*dx1**3*(7.0d0*dx2**5+dx1*(105.0d0*dx2**4+11.0d0*dx1*(42.0d0*dx2**3+dx1*(70.0d0*dx2**2+9.0d0*dx1*(5.0d0*dx2+dx1 &
)))))
      ci(14)=5.0d0*dx1**2*(21.0d0*dx2**5+13.0d0*dx1*(35.0d0*dx2**4+dx1*(210.0d0*dx2**3+11.0d0*dx1*(42.0d0*dx2**2+dx1*(35.0d0*dx2 &
+9.0d0*dx1)))))
      ci(15)=5.0d0*dx1*(3.0d0*dx2**5+7.0d0*dx1*(15.0d0*dx2**4+13.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2 &
+dx1)))))
      ci(16)=dx2**5+dx1*(75.0d0*dx2**4+7.0d0*dx1*(150.0d0*dx2**3+13.0d0*dx1*(50.0d0*dx2**2+3.0d0*dx1*(25.0d0*dx2+11.0d0*dx1))))
      ci(17)=5.0d0*(dx2**4+dx1*(30.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))))
      ci(18)=5.0d0*(2.0d0*dx2**3+dx1*(30.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+13.0d0*dx1)))
      ci(19)=5.0d0*(2.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+7.0d0*dx1))
      ci(20)=5.0d0*(dx2+3.0d0*dx1)
      ci(21)=1.0d0
    case(6)
      ci(1)=dx1**15*dx2**6
      ci(2)=3.0d0*dx1**14*dx2**5*(5.0d0*dx2+2.0d0*dx1)
      ci(3)=15.0d0*dx1**13*dx2**4*(7.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(4)=5.0d0*dx1**12*dx2**3*(91.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(45.0d0*dx2+4.0d0*dx1)))
      ci(5)=15.0d0*dx1**11*dx2**2*(91.0d0*dx2**4+dx1*(182.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(6)=3.0d0*dx1**10*dx2*(1001.0d0*dx2**5+dx1*(2730.0d0*dx2**4+dx1*(2275.0d0*dx2**3+dx1*(700.0d0*dx2**2+dx1*(75.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=dx1**9*(5005.0d0*dx2**6+dx1*(18018.0d0*dx2**5+dx1*(20475.0d0*dx2**4+dx1*(9100.0d0*dx2**3+dx1*(1575.0d0*dx2**2+dx1* &
(90.0d0*dx2+dx1))))))
      ci(8)=15.0d0*dx1**8*(429.0d0*dx2**6+dx1*(2002.0d0*dx2**5+dx1*(3003.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1*(455.0d0*dx2**2+dx1* &
(42.0d0*dx2+dx1))))))
      ci(9)=15.0d0*dx1**7*(429.0d0*dx2**6+dx1*(2574.0d0*dx2**5+7.0d0*dx1*(715.0d0*dx2**4+dx1*(572.0d0*dx2**3+dx1*(195.0d0*dx2**2 &
+dx1*(26.0d0*dx2+dx1))))))
      ci(10)=65.0d0*dx1**6*(77.0d0*dx2**6+dx1*(594.0d0*dx2**5+dx1*(1485.0d0*dx2**4+7.0d0*dx1*(220.0d0*dx2**3+dx1*(99.0d0*dx2**2 &
+dx1*(18.0d0*dx2+dx1))))))
      ci(11)=39.0d0*dx1**5*(77.0d0*dx2**6+dx1*(770.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(3300.0d0*dx2**3+7.0d0*dx1*(275.0d0*dx2**2 &
+dx1*(66.0d0*dx2+5.0d0*dx1))))))
      ci(12)=39.0d0*dx1**4*(35.0d0*dx2**6+11.0d0*dx1*(42.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(300.0d0*dx2**3+dx1*(225.0d0*dx2**2 &
+7.0d0*dx1*(10.0d0*dx2+dx1))))))
      ci(13)=65.0d0*dx1**3*(7.0d0*dx2**6+dx1*(126.0d0*dx2**5+11.0d0*dx1*(63.0d0*dx2**4+dx1*(140.0d0*dx2**3+dx1*(135.0d0*dx2**2 &
+dx1*(54.0d0*dx2+7.0d0*dx1))))))
      ci(14)=15.0d0*dx1**2*(7.0d0*dx2**6+13.0d0*dx1*(14.0d0*dx2**5+dx1*(105.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+dx1* &
(35.0d0*dx2**2+3.0d0*dx1*(6.0d0*dx2+dx1))))))
      ci(15)=15.0d0*dx1*(dx2**6+dx1*(42.0d0*dx2**5+13.0d0*dx1*(35.0d0*dx2**4+dx1*(140.0d0*dx2**3+11.0d0*dx1*(21.0d0*dx2**2+dx1* &
(14.0d0*dx2+3.0d0*dx1))))))
      ci(16)=dx2**6+dx1*(90.0d0*dx2**5+7.0d0*dx1*(225.0d0*dx2**4+13.0d0*dx1*(100.0d0*dx2**3+dx1*(225.0d0*dx2**2+11.0d0*dx1* &
(18.0d0*dx2+5.0d0*dx1)))))
      ci(17)=3.0d0*(2.0d0*dx2**5+dx1*(75.0d0*dx2**4+7.0d0*dx1*(100.0d0*dx2**3+13.0d0*dx1*(25.0d0*dx2**2+dx1*(30.0d0*dx2+11.0d0*dx1 &
)))))
      ci(18)=15.0d0*(dx2**4+dx1*(20.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(19)=5.0d0*(4.0d0*dx2**3+dx1*(45.0d0*dx2**2+7.0d0*dx1*(18.0d0*dx2+13.0d0*dx1)))
      ci(20)=15.0d0*(dx2**2+dx1*(6.0d0*dx2+7.0d0*dx1))
      ci(21)=3.0d0*(2.0d0*dx2+5.0d0*dx1)
      ci(22)=1.0d0
    case(7)
      ci(1)=dx1**15*dx2**7
      ci(2)=dx1**14*dx2**6*(15.0d0*dx2+7.0d0*dx1)
      ci(3)=21.0d0*dx1**13*dx2**5*(5.0d0*dx2**2+dx1*(5.0d0*dx2+dx1))
      ci(4)=35.0d0*dx1**12*dx2**4*(13.0d0*dx2**3+dx1*(21.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(5)=35.0d0*dx1**11*dx2**3*(39.0d0*dx2**4+dx1*(91.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))))
      ci(6)=21.0d0*dx1**10*dx2**2*(143.0d0*dx2**5+dx1*(455.0d0*dx2**4+dx1*(455.0d0*dx2**3+dx1*(175.0d0*dx2**2+dx1*(25.0d0*dx2+dx1) &
))))
      ci(7)=7.0d0*dx1**9*dx2*(715.0d0*dx2**6+dx1*(3003.0d0*dx2**5+dx1*(4095.0d0*dx2**4+dx1*(2275.0d0*dx2**3+dx1*(525.0d0*dx2**2 &
+dx1*(45.0d0*dx2+dx1))))))
      ci(8)=dx1**8*(6435.0d0*dx2**7+dx1*(35035.0d0*dx2**6+dx1*(63063.0d0*dx2**5+dx1*(47775.0d0*dx2**4+dx1*(15925.0d0*dx2**3+dx1* &
(2205.0d0*dx2**2+dx1*(105.0d0*dx2+dx1)))))))
      ci(9)=15.0d0*dx1**7*(429.0d0*dx2**7+dx1*(3003.0d0*dx2**6+dx1*(7007.0d0*dx2**5+dx1*(7007.0d0*dx2**4+dx1*(3185.0d0*dx2**3+dx1* &
(637.0d0*dx2**2+dx1*(49.0d0*dx2+dx1)))))))
      ci(10)=35.0d0*dx1**6*(143.0d0*dx2**7+dx1*(1287.0d0*dx2**6+dx1*(3861.0d0*dx2**5+dx1*(5005.0d0*dx2**4+dx1*(3003.0d0*dx2**3 &
+dx1*(819.0d0*dx2**2+dx1*(91.0d0*dx2+3.0d0*dx1)))))))
      ci(11)=91.0d0*dx1**5*(33.0d0*dx2**7+dx1*(385.0d0*dx2**6+dx1*(1485.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1*(1925.0d0*dx2**3+dx1* &
(693.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+dx1)))))))
      ci(12)=273.0d0*dx1**4*(5.0d0*dx2**7+dx1*(77.0d0*dx2**6+dx1*(385.0d0*dx2**5+dx1*(825.0d0*dx2**4+dx1*(825.0d0*dx2**3+dx1* &
(385.0d0*dx2**2+dx1*(77.0d0*dx2+5.0d0*dx1)))))))
      ci(13)=91.0d0*dx1**3*(5.0d0*dx2**7+dx1*(105.0d0*dx2**6+11.0d0*dx1*(63.0d0*dx2**5+dx1*(175.0d0*dx2**4+dx1*(225.0d0*dx2**3 &
+dx1*(135.0d0*dx2**2+dx1*(35.0d0*dx2+3.0d0*dx1)))))))
      ci(14)=35.0d0*dx1**2*(3.0d0*dx2**7+13.0d0*dx1*(7.0d0*dx2**6+dx1*(63.0d0*dx2**5+11.0d0*dx1*(21.0d0*dx2**4+dx1*(35.0d0*dx2**3 &
+dx1*(27.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))))))
      ci(15)=15.0d0*dx1*(dx2**7+dx1*(49.0d0*dx2**6+13.0d0*dx1*(49.0d0*dx2**5+dx1*(245.0d0*dx2**4+11.0d0*dx1*(49.0d0*dx2**3+dx1* &
(49.0d0*dx2**2+3.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(16)=dx2**7+dx1*(105.0d0*dx2**6+dx1*(2205.0d0*dx2**5+13.0d0*dx1*(1225.0d0*dx2**4+dx1*(3675.0d0*dx2**3+11.0d0*dx1* &
(441.0d0*dx2**2+5.0d0*dx1*(49.0d0*dx2+9.0d0*dx1))))))
      ci(17)=7.0d0*(dx2**6+dx1*(45.0d0*dx2**5+dx1*(525.0d0*dx2**4+13.0d0*dx1*(175.0d0*dx2**3+dx1*(315.0d0*dx2**2+11.0d0*dx1* &
(21.0d0*dx2+5.0d0*dx1))))))
      ci(18)=21.0d0*(dx2**5+dx1*(25.0d0*dx2**4+dx1*(175.0d0*dx2**3+13.0d0*dx1*(35.0d0*dx2**2+dx1*(35.0d0*dx2+11.0d0*dx1)))))
      ci(19)=35.0d0*(dx2**4+dx1*(15.0d0*dx2**3+dx1*(63.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+3.0d0*dx1))))
      ci(20)=35.0d0*(dx2**3+dx1*(9.0d0*dx2**2+dx1*(21.0d0*dx2+13.0d0*dx1)))
      ci(21)=21.0d0*(dx2**2+5.0d0*dx1*(dx2+dx1))
      ci(22)=7.0d0*dx2+15.0d0*dx1
      ci(23)=1.0d0
    case(8)
      ci(1)=dx1**15*dx2**8
      ci(2)=dx1**14*dx2**7*(15.0d0*dx2+8.0d0*dx1)
      ci(3)=dx1**13*dx2**6*(105.0d0*dx2**2+4.0d0*dx1*(30.0d0*dx2+7.0d0*dx1))
      ci(4)=7.0d0*dx1**12*dx2**5*(65.0d0*dx2**3+4.0d0*dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(5)=35.0d0*dx1**11*dx2**4*(39.0d0*dx2**4+2.0d0*dx1*(52.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**10*dx2**3*(429.0d0*dx2**5+2.0d0*dx1*(780.0d0*dx2**4+dx1*(910.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(75.0d0*dx2 &
+4.0d0*dx1)))))
      ci(7)=7.0d0*dx1**9*dx2**2*(715.0d0*dx2**6+2.0d0*dx1*(1716.0d0*dx2**5+dx1*(2730.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1* &
(525.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2+dx1))))))
      ci(8)=dx1**8*dx2*(6435.0d0*dx2**7+2.0d0*dx1*(20020.0d0*dx2**6+dx1*(42042.0d0*dx2**5+dx1*(38220.0d0*dx2**4+dx1* &
(15925.0d0*dx2**3+2.0d0*dx1*(1470.0d0*dx2**2+dx1*(105.0d0*dx2+2.0d0*dx1)))))))
      ci(9)=dx1**7*(6435.0d0*dx2**8+dx1*(51480.0d0*dx2**7+dx1*(140140.0d0*dx2**6+dx1*(168168.0d0*dx2**5+dx1*(95550.0d0*dx2**4+dx1* &
(25480.0d0*dx2**3+dx1*(2940.0d0*dx2**2+dx1*(120.0d0*dx2+dx1))))))))
      ci(10)=5.0d0*dx1**6*(1001.0d0*dx2**8+dx1*(10296.0d0*dx2**7+dx1*(36036.0d0*dx2**6+dx1*(56056.0d0*dx2**5+dx1*(42042.0d0*dx2**4 &
+dx1*(15288.0d0*dx2**3+dx1*(2548.0d0*dx2**2+3.0d0*dx1*(56.0d0*dx2+dx1))))))))
      ci(11)=7.0d0*dx1**5*(429.0d0*dx2**8+dx1*(5720.0d0*dx2**7+dx1*(25740.0d0*dx2**6+dx1*(51480.0d0*dx2**5+dx1*(50050.0d0*dx2**4 &
+dx1*(24024.0d0*dx2**3+5.0d0*dx1*(1092.0d0*dx2**2+dx1*(104.0d0*dx2+3.0d0*dx1))))))))
      ci(12)=91.0d0*dx1**4*(15.0d0*dx2**8+dx1*(264.0d0*dx2**7+dx1*(1540.0d0*dx2**6+dx1*(3960.0d0*dx2**5+dx1*(4950.0d0*dx2**4+dx1* &
(3080.0d0*dx2**3+dx1*(924.0d0*dx2**2+5.0d0*dx1*(24.0d0*dx2+dx1))))))))
      ci(13)=91.0d0*dx1**3*(5.0d0*dx2**8+dx1*(120.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(3080.0d0*dx2**5+dx1*(4950.0d0*dx2**4+dx1* &
(3960.0d0*dx2**3+dx1*(1540.0d0*dx2**2+3.0d0*dx1*(88.0d0*dx2+5.0d0*dx1))))))))
      ci(14)=7.0d0*dx1**2*(15.0d0*dx2**8+13.0d0*dx1*(40.0d0*dx2**7+dx1*(420.0d0*dx2**6+11.0d0*dx1*(168.0d0*dx2**5+dx1* &
(350.0d0*dx2**4+dx1*(360.0d0*dx2**3+dx1*(180.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))))))))
      ci(15)=5.0d0*dx1*(3.0d0*dx2**8+dx1*(168.0d0*dx2**7+13.0d0*dx1*(196.0d0*dx2**6+dx1*(1176.0d0*dx2**5+11.0d0*dx1* &
(294.0d0*dx2**4+dx1*(392.0d0*dx2**3+dx1*(252.0d0*dx2**2+dx1*(72.0d0*dx2+7.0d0*dx1))))))))
      ci(16)=dx2**8+dx1*(120.0d0*dx2**7+dx1*(2940.0d0*dx2**6+13.0d0*dx1*(1960.0d0*dx2**5+dx1*(7350.0d0*dx2**4+11.0d0*dx1* &
(1176.0d0*dx2**3+5.0d0*dx1*(196.0d0*dx2**2+9.0d0*dx1*(8.0d0*dx2+dx1)))))))
      ci(17)=8.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(5880.0d0*dx2**5+13.0d0*dx1*(2450.0d0*dx2**4+dx1*(5880.0d0*dx2**3+11.0d0*dx1* &
(588.0d0*dx2**2+5.0d0*dx1*(56.0d0*dx2+9.0d0*dx1))))))
      ci(18)=7.0d0*(4.0d0*dx2**6+dx1*(120.0d0*dx2**5+dx1*(1050.0d0*dx2**4+13.0d0*dx1*(280.0d0*dx2**3+dx1*(420.0d0*dx2**2 &
+11.0d0*dx1*(24.0d0*dx2+5.0d0*dx1))))))
      ci(19)=7.0d0*(8.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(840.0d0*dx2**3+13.0d0*dx1*(140.0d0*dx2**2+3.0d0*dx1*(40.0d0*dx2 &
+11.0d0*dx1)))))
      ci(20)=35.0d0*(2.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(84.0d0*dx2**2+13.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(21)=7.0d0*(8.0d0*dx2**3+5.0d0*dx1*(12.0d0*dx2**2+dx1*(24.0d0*dx2+13.0d0*dx1)))
      ci(22)=28.0d0*dx2**2+15.0d0*dx1*(8.0d0*dx2+7.0d0*dx1)
      ci(23)=8.0d0*dx2+15.0d0*dx1
      ci(24)=1.0d0
    case(9)
      ci(1)=dx1**15*dx2**9
      ci(2)=3.0d0*dx1**14*dx2**8*(5.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**13*dx2**7*(35.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+4.0d0*dx1))
      ci(4)=dx1**12*dx2**6*(455.0d0*dx2**3+3.0d0*dx1*(315.0d0*dx2**2+4.0d0*dx1*(45.0d0*dx2+7.0d0*dx1)))
      ci(5)=21.0d0*dx1**11*dx2**5*(65.0d0*dx2**4+3.0d0*dx1*(65.0d0*dx2**3+2.0d0*dx1*(30.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(6)=21.0d0*dx1**10*dx2**4*(143.0d0*dx2**5+3.0d0*dx1*(195.0d0*dx2**4+2.0d0*dx1*(130.0d0*dx2**3+dx1*(70.0d0*dx2**2+dx1* &
(15.0d0*dx2+dx1)))))
      ci(7)=7.0d0*dx1**9*dx2**3*(715.0d0*dx2**6+3.0d0*dx1*(1287.0d0*dx2**5+2.0d0*dx1*(1170.0d0*dx2**4+dx1*(910.0d0*dx2**3+dx1* &
(315.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1))))))
      ci(8)=9.0d0*dx1**8*dx2**2*(715.0d0*dx2**7+dx1*(5005.0d0*dx2**6+2.0d0*dx1*(6006.0d0*dx2**5+dx1*(6370.0d0*dx2**4+dx1* &
(3185.0d0*dx2**3+dx1*(735.0d0*dx2**2+2.0d0*dx1*(35.0d0*dx2+dx1)))))))
      ci(9)=9.0d0*dx1**7*dx2*(715.0d0*dx2**8+dx1*(6435.0d0*dx2**7+dx1*(20020.0d0*dx2**6+dx1*(28028.0d0*dx2**5+dx1* &
(19110.0d0*dx2**4+dx1*(6370.0d0*dx2**3+dx1*(980.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))))))
      ci(10)=dx1**6*(5005.0d0*dx2**9+dx1*(57915.0d0*dx2**8+dx1*(231660.0d0*dx2**7+dx1*(420420.0d0*dx2**6+dx1*(378378.0d0*dx2**5 &
+dx1*(171990.0d0*dx2**4+dx1*(38220.0d0*dx2**3+dx1*(3780.0d0*dx2**2+dx1*(135.0d0*dx2+dx1)))))))))
      ci(11)=3.0d0*dx1**5*(1001.0d0*dx2**9+dx1*(15015.0d0*dx2**8+dx1*(77220.0d0*dx2**7+dx1*(180180.0d0*dx2**6+dx1* &
(210210.0d0*dx2**5+dx1*(126126.0d0*dx2**4+5.0d0*dx1*(7644.0d0*dx2**3+dx1*(1092.0d0*dx2**2+dx1*(63.0d0*dx2+dx1)))))))))
      ci(12)=21.0d0*dx1**4*(65.0d0*dx2**9+dx1*(1287.0d0*dx2**8+dx1*(8580.0d0*dx2**7+dx1*(25740.0d0*dx2**6+dx1*(38610.0d0*dx2**5 &
+dx1*(30030.0d0*dx2**4+dx1*(12012.0d0*dx2**3+5.0d0*dx1*(468.0d0*dx2**2+dx1*(39.0d0*dx2+dx1)))))))))
      ci(13)=91.0d0*dx1**3*(5.0d0*dx2**9+dx1*(135.0d0*dx2**8+dx1*(1188.0d0*dx2**7+dx1*(4620.0d0*dx2**6+dx1*(8910.0d0*dx2**5+dx1* &
(8910.0d0*dx2**4+dx1*(4620.0d0*dx2**3+dx1*(1188.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+dx1)))))))))
      ci(14)=21.0d0*dx1**2*(5.0d0*dx2**9+13.0d0*dx1*(15.0d0*dx2**8+dx1*(180.0d0*dx2**7+dx1*(924.0d0*dx2**6+dx1*(2310.0d0*dx2**5 &
+dx1*(2970.0d0*dx2**4+dx1*(1980.0d0*dx2**3+dx1*(660.0d0*dx2**2+dx1*(99.0d0*dx2+5.0d0*dx1)))))))))
      ci(15)=3.0d0*dx1*(5.0d0*dx2**9+dx1*(315.0d0*dx2**8+13.0d0*dx1*(420.0d0*dx2**7+dx1*(2940.0d0*dx2**6+11.0d0*dx1* &
(882.0d0*dx2**5+dx1*(1470.0d0*dx2**4+dx1*(1260.0d0*dx2**3+dx1*(540.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+dx1)))))))))
      ci(16)=dx2**9+dx1*(135.0d0*dx2**8+dx1*(3780.0d0*dx2**7+13.0d0*dx1*(2940.0d0*dx2**6+dx1*(13230.0d0*dx2**5+11.0d0*dx1* &
(2646.0d0*dx2**4+5.0d0*dx1*(588.0d0*dx2**3+dx1*(324.0d0*dx2**2+dx1*(81.0d0*dx2+7.0d0*dx1))))))))
      ci(17)=9.0d0*(dx2**8+dx1*(60.0d0*dx2**7+dx1*(980.0d0*dx2**6+13.0d0*dx1*(490.0d0*dx2**5+dx1*(1470.0d0*dx2**4+11.0d0*dx1* &
(196.0d0*dx2**3+5.0d0*dx1*(28.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))))))))
      ci(18)=9.0d0*(4.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(1470.0d0*dx2**5+13.0d0*dx1*(490.0d0*dx2**4+dx1*(980.0d0*dx2**3 &
+11.0d0*dx1*(84.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(19)=7.0d0*(12.0d0*dx2**6+dx1*(270.0d0*dx2**5+dx1*(1890.0d0*dx2**4+13.0d0*dx1*(420.0d0*dx2**3+dx1*(540.0d0*dx2**2 &
+11.0d0*dx1*(27.0d0*dx2+5.0d0*dx1))))))
      ci(20)=21.0d0*(6.0d0*dx2**5+dx1*(90.0d0*dx2**4+dx1*(420.0d0*dx2**3+13.0d0*dx1*(60.0d0*dx2**2+dx1*(45.0d0*dx2+11.0d0*dx1)))))
      ci(21)=21.0d0*(6.0d0*dx2**4+5.0d0*dx1*(12.0d0*dx2**3+dx1*(36.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+dx1))))
      ci(22)=84.0d0*dx2**3+5.0d0*dx1*(108.0d0*dx2**2+7.0d0*dx1*(27.0d0*dx2+13.0d0*dx1))
      ci(23)=3.0d0*(12.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+7.0d0*dx1))
      ci(24)=3.0d0*(3.0d0*dx2+5.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>9, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^16 * (x-x2)^n2 as sum_k=0^(16+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_16(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**16
      ci(2)=16.0d0*dx1**15
      ci(3)=120.0d0*dx1**14
      ci(4)=560.0d0*dx1**13
      ci(5)=1820.0d0*dx1**12
      ci(6)=4368.0d0*dx1**11
      ci(7)=8008.0d0*dx1**10
      ci(8)=11440.0d0*dx1**9
      ci(9)=12870.0d0*dx1**8
      ci(10)=11440.0d0*dx1**7
      ci(11)=8008.0d0*dx1**6
      ci(12)=4368.0d0*dx1**5
      ci(13)=1820.0d0*dx1**4
      ci(14)=560.0d0*dx1**3
      ci(15)=120.0d0*dx1**2
      ci(16)=16.0d0*dx1
      ci(17)=1.0d0
    case(1)
      ci(1)=dx1**16*dx2
      ci(2)=dx1**15*(16.0d0*dx2+dx1)
      ci(3)=8.0d0*dx1**14*(15.0d0*dx2+2.0d0*dx1)
      ci(4)=40.0d0*dx1**13*(14.0d0*dx2+3.0d0*dx1)
      ci(5)=140.0d0*dx1**12*(13.0d0*dx2+4.0d0*dx1)
      ci(6)=364.0d0*dx1**11*(12.0d0*dx2+5.0d0*dx1)
      ci(7)=728.0d0*dx1**10*(11.0d0*dx2+6.0d0*dx1)
      ci(8)=1144.0d0*dx1**9*(10.0d0*dx2+7.0d0*dx1)
      ci(9)=1430.0d0*dx1**8*(9.0d0*dx2+8.0d0*dx1)
      ci(10)=1430.0d0*dx1**7*(8.0d0*dx2+9.0d0*dx1)
      ci(11)=1144.0d0*dx1**6*(7.0d0*dx2+10.0d0*dx1)
      ci(12)=728.0d0*dx1**5*(6.0d0*dx2+11.0d0*dx1)
      ci(13)=364.0d0*dx1**4*(5.0d0*dx2+12.0d0*dx1)
      ci(14)=140.0d0*dx1**3*(4.0d0*dx2+13.0d0*dx1)
      ci(15)=40.0d0*dx1**2*(3.0d0*dx2+14.0d0*dx1)
      ci(16)=8.0d0*dx1*(2.0d0*dx2+15.0d0*dx1)
      ci(17)=dx2+16.0d0*dx1
      ci(18)=1.0d0
    case(2)
      ci(1)=dx1**16*dx2**2
      ci(2)=2.0d0*dx1**15*dx2*(8.0d0*dx2+dx1)
      ci(3)=dx1**14*(120.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))
      ci(4)=16.0d0*dx1**13*(35.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))
      ci(5)=20.0d0*dx1**12*(91.0d0*dx2**2+2.0d0*dx1*(28.0d0*dx2+3.0d0*dx1))
      ci(6)=56.0d0*dx1**11*(78.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+2.0d0*dx1))
      ci(7)=364.0d0*dx1**10*(22.0d0*dx2**2+dx1*(24.0d0*dx2+5.0d0*dx1))
      ci(8)=208.0d0*dx1**9*(55.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+3.0d0*dx1))
      ci(9)=286.0d0*dx1**8*(45.0d0*dx2**2+4.0d0*dx1*(20.0d0*dx2+7.0d0*dx1))
      ci(10)=2860.0d0*dx1**7*(4.0d0*dx2**2+dx1*(9.0d0*dx2+4.0d0*dx1))
      ci(11)=286.0d0*dx1**6*(28.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+9.0d0*dx1))
      ci(12)=208.0d0*dx1**5*(21.0d0*dx2**2+11.0d0*dx1*(7.0d0*dx2+5.0d0*dx1))
      ci(13)=364.0d0*dx1**4*(5.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+11.0d0*dx1))
      ci(14)=56.0d0*dx1**3*(10.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+6.0d0*dx1))
      ci(15)=20.0d0*dx1**2*(6.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+13.0d0*dx1))
      ci(16)=16.0d0*dx1*(dx2**2+5.0d0*dx1*(3.0d0*dx2+7.0d0*dx1))
      ci(17)=dx2**2+8.0d0*dx1*(4.0d0*dx2+15.0d0*dx1)
      ci(18)=2.0d0*(dx2+8.0d0*dx1)
      ci(19)=1.0d0
    case(3)
      ci(1)=dx1**16*dx2**3
      ci(2)=dx1**15*dx2**2*(16.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**14*dx2*(40.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))
      ci(4)=dx1**13*(560.0d0*dx2**3+dx1*(360.0d0*dx2**2+dx1*(48.0d0*dx2+dx1)))
      ci(5)=4.0d0*dx1**12*(455.0d0*dx2**3+2.0d0*dx1*(210.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1)))
      ci(6)=12.0d0*dx1**11*(364.0d0*dx2**3+5.0d0*dx1*(91.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+dx1)))
      ci(7)=28.0d0*dx1**10*(286.0d0*dx2**3+dx1*(468.0d0*dx2**2+5.0d0*dx1*(39.0d0*dx2+4.0d0*dx1)))
      ci(8)=52.0d0*dx1**9*(220.0d0*dx2**3+7.0d0*dx1*(66.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1)))
      ci(9)=78.0d0*dx1**8*(165.0d0*dx2**3+4.0d0*dx1*(110.0d0*dx2**2+7.0d0*dx1*(11.0d0*dx2+2.0d0*dx1)))
      ci(10)=286.0d0*dx1**7*(40.0d0*dx2**3+dx1*(135.0d0*dx2**2+4.0d0*dx1*(30.0d0*dx2+7.0d0*dx1)))
      ci(11)=286.0d0*dx1**6*(28.0d0*dx2**3+5.0d0*dx1*(24.0d0*dx2**2+dx1*(27.0d0*dx2+8.0d0*dx1)))
      ci(12)=78.0d0*dx1**5*(56.0d0*dx2**3+11.0d0*dx1*(28.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+3.0d0*dx1)))
      ci(13)=52.0d0*dx1**4*(35.0d0*dx2**3+2.0d0*dx1*(126.0d0*dx2**2+11.0d0*dx1*(21.0d0*dx2+10.0d0*dx1)))
      ci(14)=28.0d0*dx1**3*(20.0d0*dx2**3+13.0d0*dx1*(15.0d0*dx2**2+2.0d0*dx1*(18.0d0*dx2+11.0d0*dx1)))
      ci(15)=12.0d0*dx1**2*(10.0d0*dx2**3+7.0d0*dx1*(20.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+4.0d0*dx1)))
      ci(16)=4.0d0*dx1*(4.0d0*dx2**3+5.0d0*dx1*(18.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+13.0d0*dx1)))
      ci(17)=dx2**3+8.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+14.0d0*dx1))
      ci(18)=3.0d0*(dx2**2+8.0d0*dx1*(2.0d0*dx2+5.0d0*dx1))
      ci(19)=3.0d0*dx2+16.0d0*dx1
      ci(20)=1.0d0
    case(4)
      ci(1)=dx1**16*dx2**4
      ci(2)=4.0d0*dx1**15*dx2**3*(4.0d0*dx2+dx1)
      ci(3)=2.0d0*dx1**14*dx2**2*(60.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1))
      ci(4)=4.0d0*dx1**13*dx2*(140.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(24.0d0*dx2+dx1)))
      ci(5)=dx1**12*(1820.0d0*dx2**4+dx1*(2240.0d0*dx2**3+dx1*(720.0d0*dx2**2+dx1*(64.0d0*dx2+dx1))))
      ci(6)=16.0d0*dx1**11*(273.0d0*dx2**4+dx1*(455.0d0*dx2**3+dx1*(210.0d0*dx2**2+dx1*(30.0d0*dx2+dx1))))
      ci(7)=8.0d0*dx1**10*(1001.0d0*dx2**4+dx1*(2184.0d0*dx2**3+5.0d0*dx1*(273.0d0*dx2**2+dx1*(56.0d0*dx2+3.0d0*dx1))))
      ci(8)=16.0d0*dx1**9*(715.0d0*dx2**4+7.0d0*dx1*(286.0d0*dx2**3+dx1*(234.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+dx1))))
      ci(9)=26.0d0*dx1**8*(495.0d0*dx2**4+2.0d0*dx1*(880.0d0*dx2**3+7.0d0*dx1*(132.0d0*dx2**2+dx1*(48.0d0*dx2+5.0d0*dx1))))
      ci(10)=104.0d0*dx1**7*(110.0d0*dx2**4+dx1*(495.0d0*dx2**3+2.0d0*dx1*(330.0d0*dx2**2+7.0d0*dx1*(22.0d0*dx2+3.0d0*dx1))))
      ci(11)=572.0d0*dx1**6*(14.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(135.0d0*dx2**2+2.0d0*dx1*(40.0d0*dx2+7.0d0*dx1))))
      ci(12)=104.0d0*dx1**5*(42.0d0*dx2**4+11.0d0*dx1*(28.0d0*dx2**3+5.0d0*dx1*(12.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1))))
      ci(13)=26.0d0*dx1**4*(70.0d0*dx2**4+dx1*(672.0d0*dx2**3+11.0d0*dx1*(168.0d0*dx2**2+5.0d0*dx1*(32.0d0*dx2+9.0d0*dx1))))
      ci(14)=16.0d0*dx1**3*(35.0d0*dx2**4+13.0d0*dx1*(35.0d0*dx2**3+dx1*(126.0d0*dx2**2+11.0d0*dx1*(14.0d0*dx2+5.0d0*dx1))))
      ci(15)=8.0d0*dx1**2*(15.0d0*dx2**4+7.0d0*dx1*(40.0d0*dx2**3+13.0d0*dx1*(15.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1))))
      ci(16)=16.0d0*dx1*(dx2**4+dx1*(30.0d0*dx2**3+7.0d0*dx1*(30.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))))
      ci(17)=dx2**4+4.0d0*dx1*(16.0d0*dx2**3+5.0d0*dx1*(36.0d0*dx2**2+7.0d0*dx1*(16.0d0*dx2+13.0d0*dx1)))
      ci(18)=4.0d0*(dx2**3+4.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+7.0d0*dx1)))
      ci(19)=2.0d0*(3.0d0*dx2**2+4.0d0*dx1*(8.0d0*dx2+15.0d0*dx1))
      ci(20)=4.0d0*(dx2+4.0d0*dx1)
      ci(21)=1.0d0
    case(5)
      ci(1)=dx1**16*dx2**5
      ci(2)=dx1**15*dx2**4*(16.0d0*dx2+5.0d0*dx1)
      ci(3)=10.0d0*dx1**14*dx2**3*(12.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(4)=10.0d0*dx1**13*dx2**2*(56.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(16.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**12*dx2*(364.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))
      ci(6)=dx1**11*(4368.0d0*dx2**5+dx1*(9100.0d0*dx2**4+dx1*(5600.0d0*dx2**3+dx1*(1200.0d0*dx2**2+dx1*(80.0d0*dx2+dx1)))))
      ci(7)=8.0d0*dx1**10*(1001.0d0*dx2**5+dx1*(2730.0d0*dx2**4+dx1*(2275.0d0*dx2**3+dx1*(700.0d0*dx2**2+dx1*(75.0d0*dx2+2.0d0*dx1 &
)))))
      ci(8)=40.0d0*dx1**9*(286.0d0*dx2**5+dx1*(1001.0d0*dx2**4+dx1*(1092.0d0*dx2**3+dx1*(455.0d0*dx2**2+dx1*(70.0d0*dx2+3.0d0*dx1) &
))))
      ci(9)=10.0d0*dx1**8*(1287.0d0*dx2**5+2.0d0*dx1*(2860.0d0*dx2**4+7.0d0*dx1*(572.0d0*dx2**3+dx1*(312.0d0*dx2**2+dx1* &
(65.0d0*dx2+4.0d0*dx1)))))
      ci(10)=130.0d0*dx1**7*(88.0d0*dx2**5+dx1*(495.0d0*dx2**4+2.0d0*dx1*(440.0d0*dx2**3+7.0d0*dx1*(44.0d0*dx2**2+dx1*(12.0d0*dx2 &
+dx1)))))
      ci(11)=52.0d0*dx1**6*(154.0d0*dx2**5+dx1*(1100.0d0*dx2**4+dx1*(2475.0d0*dx2**3+2.0d0*dx1*(1100.0d0*dx2**2+7.0d0*dx1* &
(55.0d0*dx2+6.0d0*dx1)))))
      ci(12)=52.0d0*dx1**5*(84.0d0*dx2**5+11.0d0*dx1*(70.0d0*dx2**4+dx1*(200.0d0*dx2**3+dx1*(225.0d0*dx2**2+2.0d0*dx1*(50.0d0*dx2 &
+7.0d0*dx1)))))
      ci(13)=130.0d0*dx1**4*(14.0d0*dx2**5+dx1*(168.0d0*dx2**4+11.0d0*dx1*(56.0d0*dx2**3+dx1*(80.0d0*dx2**2+dx1*(45.0d0*dx2 &
+8.0d0*dx1)))))
      ci(14)=10.0d0*dx1**3*(56.0d0*dx2**5+13.0d0*dx1*(70.0d0*dx2**4+dx1*(336.0d0*dx2**3+11.0d0*dx1*(56.0d0*dx2**2+dx1*(40.0d0*dx2 &
+9.0d0*dx1)))))
      ci(15)=40.0d0*dx1**2*(3.0d0*dx2**5+dx1*(70.0d0*dx2**4+13.0d0*dx1*(35.0d0*dx2**3+dx1*(84.0d0*dx2**2+11.0d0*dx1*(7.0d0*dx2 &
+2.0d0*dx1)))))
      ci(16)=8.0d0*dx1*(2.0d0*dx2**5+dx1*(75.0d0*dx2**4+7.0d0*dx1*(100.0d0*dx2**3+13.0d0*dx1*(25.0d0*dx2**2+dx1*(30.0d0*dx2 &
+11.0d0*dx1)))))
      ci(17)=dx2**5+4.0d0*dx1*(20.0d0*dx2**4+dx1*(300.0d0*dx2**3+7.0d0*dx1*(200.0d0*dx2**2+13.0d0*dx1*(25.0d0*dx2+12.0d0*dx1))))
      ci(18)=5.0d0*(dx2**4+4.0d0*dx1*(8.0d0*dx2**3+dx1*(60.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+13.0d0*dx1))))
      ci(19)=10.0d0*(dx2**3+4.0d0*dx1*(4.0d0*dx2**2+dx1*(15.0d0*dx2+14.0d0*dx1)))
      ci(20)=10.0d0*(dx2**2+4.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(21)=5.0d0*dx2+16.0d0*dx1
      ci(22)=1.0d0
    case(6)
      ci(1)=dx1**16*dx2**6
      ci(2)=2.0d0*dx1**15*dx2**5*(8.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**14*dx2**4*(40.0d0*dx2**2+dx1*(32.0d0*dx2+5.0d0*dx1))
      ci(4)=20.0d0*dx1**13*dx2**3*(28.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**12*dx2**2*(364.0d0*dx2**4+dx1*(672.0d0*dx2**3+dx1*(360.0d0*dx2**2+dx1*(64.0d0*dx2+3.0d0*dx1))))
      ci(6)=6.0d0*dx1**11*dx2*(728.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1400.0d0*dx2**3+dx1*(400.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))) &
))
      ci(7)=dx1**10*(8008.0d0*dx2**6+dx1*(26208.0d0*dx2**5+dx1*(27300.0d0*dx2**4+dx1*(11200.0d0*dx2**3+dx1*(1800.0d0*dx2**2+dx1* &
(96.0d0*dx2+dx1))))))
      ci(8)=16.0d0*dx1**9*(715.0d0*dx2**6+dx1*(3003.0d0*dx2**5+dx1*(4095.0d0*dx2**4+dx1*(2275.0d0*dx2**3+dx1*(525.0d0*dx2**2+dx1* &
(45.0d0*dx2+dx1))))))
      ci(9)=30.0d0*dx1**8*(429.0d0*dx2**6+2.0d0*dx1*(1144.0d0*dx2**5+dx1*(2002.0d0*dx2**4+dx1*(1456.0d0*dx2**3+dx1*(455.0d0*dx2**2 &
+2.0d0*dx1*(28.0d0*dx2+dx1))))))
      ci(10)=20.0d0*dx1**7*(572.0d0*dx2**6+dx1*(3861.0d0*dx2**5+2.0d0*dx1*(4290.0d0*dx2**4+7.0d0*dx1*(572.0d0*dx2**3+dx1* &
(234.0d0*dx2**2+dx1*(39.0d0*dx2+2.0d0*dx1))))))
      ci(11)=26.0d0*dx1**6*(308.0d0*dx2**6+dx1*(2640.0d0*dx2**5+dx1*(7425.0d0*dx2**4+2.0d0*dx1*(4400.0d0*dx2**3+7.0d0*dx1* &
(330.0d0*dx2**2+dx1*(72.0d0*dx2+5.0d0*dx1))))))
      ci(12)=312.0d0*dx1**5*(14.0d0*dx2**6+dx1*(154.0d0*dx2**5+dx1*(550.0d0*dx2**4+dx1*(825.0d0*dx2**3+2.0d0*dx1*(275.0d0*dx2**2 &
+7.0d0*dx1*(11.0d0*dx2+dx1))))))
      ci(13)=26.0d0*dx1**4*(70.0d0*dx2**6+dx1*(1008.0d0*dx2**5+11.0d0*dx1*(420.0d0*dx2**4+dx1*(800.0d0*dx2**3+dx1*(675.0d0*dx2**2 &
+4.0d0*dx1*(60.0d0*dx2+7.0d0*dx1))))))
      ci(14)=20.0d0*dx1**3*(28.0d0*dx2**6+13.0d0*dx1*(42.0d0*dx2**5+dx1*(252.0d0*dx2**4+11.0d0*dx1*(56.0d0*dx2**3+dx1* &
(60.0d0*dx2**2+dx1*(27.0d0*dx2+4.0d0*dx1))))))
      ci(15)=30.0d0*dx1**2*(4.0d0*dx2**6+dx1*(112.0d0*dx2**5+13.0d0*dx1*(70.0d0*dx2**4+dx1*(224.0d0*dx2**3+11.0d0*dx1* &
(28.0d0*dx2**2+dx1*(16.0d0*dx2+3.0d0*dx1))))))
      ci(16)=16.0d0*dx1*(dx2**6+dx1*(45.0d0*dx2**5+dx1*(525.0d0*dx2**4+13.0d0*dx1*(175.0d0*dx2**3+dx1*(315.0d0*dx2**2+11.0d0*dx1* &
(21.0d0*dx2+5.0d0*dx1))))))
      ci(17)=dx2**6+4.0d0*dx1*(24.0d0*dx2**5+dx1*(450.0d0*dx2**4+7.0d0*dx1*(400.0d0*dx2**3+13.0d0*dx1*(75.0d0*dx2**2+2.0d0*dx1* &
(36.0d0*dx2+11.0d0*dx1)))))
      ci(18)=6.0d0*(dx2**5+4.0d0*dx1*(10.0d0*dx2**4+dx1*(100.0d0*dx2**3+7.0d0*dx1*(50.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)) &
)))
      ci(19)=5.0d0*(3.0d0*dx2**4+4.0d0*dx1*(16.0d0*dx2**3+dx1*(90.0d0*dx2**2+7.0d0*dx1*(24.0d0*dx2+13.0d0*dx1))))
      ci(20)=20.0d0*(dx2**3+4.0d0*dx1*(3.0d0*dx2**2+dx1*(9.0d0*dx2+7.0d0*dx1)))
      ci(21)=3.0d0*(5.0d0*dx2**2+8.0d0*dx1*(4.0d0*dx2+5.0d0*dx1))
      ci(22)=2.0d0*(3.0d0*dx2+8.0d0*dx1)
      ci(23)=1.0d0
    case(7)
      ci(1)=dx1**16*dx2**7
      ci(2)=dx1**15*dx2**6*(16.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**14*dx2**5*(120.0d0*dx2**2+7.0d0*dx1*(16.0d0*dx2+3.0d0*dx1))
      ci(4)=7.0d0*dx1**13*dx2**4*(80.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(48.0d0*dx2+5.0d0*dx1)))
      ci(5)=35.0d0*dx1**12*dx2**3*(52.0d0*dx2**4+dx1*(112.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**11*dx2**2*(624.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1680.0d0*dx2**3+dx1*(600.0d0*dx2**2+dx1*(80.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=7.0d0*dx1**10*dx2*(1144.0d0*dx2**6+dx1*(4368.0d0*dx2**5+dx1*(5460.0d0*dx2**4+dx1*(2800.0d0*dx2**3+dx1*(600.0d0*dx2**2 &
+dx1*(48.0d0*dx2+dx1))))))
      ci(8)=dx1**9*(11440.0d0*dx2**7+dx1*(56056.0d0*dx2**6+dx1*(91728.0d0*dx2**5+dx1*(63700.0d0*dx2**4+dx1*(19600.0d0*dx2**3+dx1* &
(2520.0d0*dx2**2+dx1*(112.0d0*dx2+dx1)))))))
      ci(9)=2.0d0*dx1**8*(6435.0d0*dx2**7+2.0d0*dx1*(20020.0d0*dx2**6+dx1*(42042.0d0*dx2**5+dx1*(38220.0d0*dx2**4+dx1* &
(15925.0d0*dx2**3+2.0d0*dx1*(1470.0d0*dx2**2+dx1*(105.0d0*dx2+2.0d0*dx1)))))))
      ci(10)=10.0d0*dx1**7*(1144.0d0*dx2**7+dx1*(9009.0d0*dx2**6+2.0d0*dx1*(12012.0d0*dx2**5+dx1*(14014.0d0*dx2**4+dx1* &
(7644.0d0*dx2**3+dx1*(1911.0d0*dx2**2+2.0d0*dx1*(98.0d0*dx2+3.0d0*dx1)))))))
      ci(11)=14.0d0*dx1**6*(572.0d0*dx2**7+dx1*(5720.0d0*dx2**6+dx1*(19305.0d0*dx2**5+2.0d0*dx1*(14300.0d0*dx2**4+dx1* &
(10010.0d0*dx2**3+dx1*(3276.0d0*dx2**2+5.0d0*dx1*(91.0d0*dx2+4.0d0*dx1)))))))
      ci(12)=182.0d0*dx1**5*(24.0d0*dx2**7+dx1*(308.0d0*dx2**6+dx1*(1320.0d0*dx2**5+dx1*(2475.0d0*dx2**4+2.0d0*dx1* &
(1100.0d0*dx2**3+dx1*(462.0d0*dx2**2+dx1*(84.0d0*dx2+5.0d0*dx1)))))))
      ci(13)=182.0d0*dx1**4*(10.0d0*dx2**7+dx1*(168.0d0*dx2**6+dx1*(924.0d0*dx2**5+dx1*(2200.0d0*dx2**4+dx1*(2475.0d0*dx2**3 &
+4.0d0*dx1*(330.0d0*dx2**2+dx1*(77.0d0*dx2+6.0d0*dx1)))))))
      ci(14)=14.0d0*dx1**3*(40.0d0*dx2**7+13.0d0*dx1*(70.0d0*dx2**6+dx1*(504.0d0*dx2**5+11.0d0*dx1*(140.0d0*dx2**4+dx1* &
(200.0d0*dx2**3+dx1*(135.0d0*dx2**2+4.0d0*dx1*(10.0d0*dx2+dx1)))))))
      ci(15)=10.0d0*dx1**2*(12.0d0*dx2**7+dx1*(392.0d0*dx2**6+13.0d0*dx1*(294.0d0*dx2**5+dx1*(1176.0d0*dx2**4+11.0d0*dx1* &
(196.0d0*dx2**3+dx1*(168.0d0*dx2**2+dx1*(63.0d0*dx2+8.0d0*dx1)))))))
      ci(16)=2.0d0*dx1*(8.0d0*dx2**7+dx1*(420.0d0*dx2**6+dx1*(5880.0d0*dx2**5+13.0d0*dx1*(2450.0d0*dx2**4+dx1*(5880.0d0*dx2**3 &
+11.0d0*dx1*(588.0d0*dx2**2+5.0d0*dx1*(56.0d0*dx2+9.0d0*dx1)))))))
      ci(17)=dx2**7+4.0d0*dx1*(28.0d0*dx2**6+dx1*(630.0d0*dx2**5+dx1*(4900.0d0*dx2**4+13.0d0*dx1*(1225.0d0*dx2**3+2.0d0*dx1* &
(882.0d0*dx2**2+11.0d0*dx1*(49.0d0*dx2+10.0d0*dx1))))))
      ci(18)=7.0d0*(dx2**6+4.0d0*dx1*(12.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(700.0d0*dx2**3+13.0d0*dx1*(105.0d0*dx2**2+2.0d0*dx1* &
(42.0d0*dx2+11.0d0*dx1))))))
      ci(19)=7.0d0*(3.0d0*dx2**5+4.0d0*dx1*(20.0d0*dx2**4+dx1*(150.0d0*dx2**3+dx1*(420.0d0*dx2**2+13.0d0*dx1*(35.0d0*dx2 &
+12.0d0*dx1)))))
      ci(20)=35.0d0*(dx2**4+4.0d0*dx1*(4.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(28.0d0*dx2+13.0d0*dx1))))
      ci(21)=7.0d0*(5.0d0*dx2**3+8.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(22)=21.0d0*dx2**2+8.0d0*dx1*(14.0d0*dx2+15.0d0*dx1)
      ci(23)=7.0d0*dx2+16.0d0*dx1
      ci(24)=1.0d0
    case(8)
      ci(1)=dx1**16*dx2**8
      ci(2)=8.0d0*dx1**15*dx2**7*(2.0d0*dx2+dx1)
      ci(3)=4.0d0*dx1**14*dx2**6*(30.0d0*dx2**2+dx1*(32.0d0*dx2+7.0d0*dx1))
      ci(4)=8.0d0*dx1**13*dx2**5*(70.0d0*dx2**3+dx1*(120.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+dx1)))
      ci(5)=14.0d0*dx1**12*dx2**4*(130.0d0*dx2**4+dx1*(320.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(64.0d0*dx2+5.0d0*dx1))))
      ci(6)=56.0d0*dx1**11*dx2**3*(78.0d0*dx2**5+dx1*(260.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)) &
)))
      ci(7)=28.0d0*dx1**10*dx2**2*(286.0d0*dx2**6+dx1*(1248.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1120.0d0*dx2**3+dx1* &
(300.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))))
      ci(8)=8.0d0*dx1**9*dx2*(1430.0d0*dx2**7+dx1*(8008.0d0*dx2**6+dx1*(15288.0d0*dx2**5+dx1*(12740.0d0*dx2**4+dx1* &
(4900.0d0*dx2**3+dx1*(840.0d0*dx2**2+dx1*(56.0d0*dx2+dx1)))))))
      ci(9)=dx1**8*(12870.0d0*dx2**8+dx1*(91520.0d0*dx2**7+dx1*(224224.0d0*dx2**6+dx1*(244608.0d0*dx2**5+dx1*(127400.0d0*dx2**4 &
+dx1*(31360.0d0*dx2**3+dx1*(3360.0d0*dx2**2+dx1*(128.0d0*dx2+dx1))))))))
      ci(10)=16.0d0*dx1**7*(715.0d0*dx2**8+dx1*(6435.0d0*dx2**7+dx1*(20020.0d0*dx2**6+dx1*(28028.0d0*dx2**5+dx1*(19110.0d0*dx2**4 &
+dx1*(6370.0d0*dx2**3+dx1*(980.0d0*dx2**2+dx1*(60.0d0*dx2+dx1))))))))
      ci(11)=8.0d0*dx1**6*(1001.0d0*dx2**8+dx1*(11440.0d0*dx2**7+dx1*(45045.0d0*dx2**6+dx1*(80080.0d0*dx2**5+dx1*(70070.0d0*dx2**4 &
+dx1*(30576.0d0*dx2**3+5.0d0*dx1*(1274.0d0*dx2**2+dx1*(112.0d0*dx2+3.0d0*dx1))))))))
      ci(12)=112.0d0*dx1**5*(39.0d0*dx2**8+dx1*(572.0d0*dx2**7+dx1*(2860.0d0*dx2**6+dx1*(6435.0d0*dx2**5+dx1*(7150.0d0*dx2**4+dx1* &
(4004.0d0*dx2**3+dx1*(1092.0d0*dx2**2+5.0d0*dx1*(26.0d0*dx2+dx1))))))))
      ci(13)=364.0d0*dx1**4*(5.0d0*dx2**8+dx1*(96.0d0*dx2**7+dx1*(616.0d0*dx2**6+dx1*(1760.0d0*dx2**5+dx1*(2475.0d0*dx2**4+dx1* &
(1760.0d0*dx2**3+dx1*(616.0d0*dx2**2+dx1*(96.0d0*dx2+5.0d0*dx1))))))))
      ci(14)=112.0d0*dx1**3*(5.0d0*dx2**8+13.0d0*dx1*(10.0d0*dx2**7+dx1*(84.0d0*dx2**6+dx1*(308.0d0*dx2**5+dx1*(550.0d0*dx2**4 &
+dx1*(495.0d0*dx2**3+dx1*(220.0d0*dx2**2+dx1*(44.0d0*dx2+3.0d0*dx1))))))))
      ci(15)=8.0d0*dx1**2*(15.0d0*dx2**8+dx1*(560.0d0*dx2**7+13.0d0*dx1*(490.0d0*dx2**6+dx1*(2352.0d0*dx2**5+11.0d0*dx1* &
(490.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(315.0d0*dx2**2+dx1*(80.0d0*dx2+7.0d0*dx1))))))))
      ci(16)=16.0d0*dx1*(dx2**8+dx1*(60.0d0*dx2**7+dx1*(980.0d0*dx2**6+13.0d0*dx1*(490.0d0*dx2**5+dx1*(1470.0d0*dx2**4+11.0d0*dx1* &
(196.0d0*dx2**3+5.0d0*dx1*(28.0d0*dx2**2+dx1*(9.0d0*dx2+dx1))))))))
      ci(17)=dx2**8+2.0d0*dx1*(64.0d0*dx2**7+dx1*(1680.0d0*dx2**6+dx1*(15680.0d0*dx2**5+13.0d0*dx1*(4900.0d0*dx2**4+dx1* &
(9408.0d0*dx2**3+11.0d0*dx1*(784.0d0*dx2**2+5.0d0*dx1*(64.0d0*dx2+9.0d0*dx1)))))))
      ci(18)=8.0d0*(dx2**7+2.0d0*dx1*(28.0d0*dx2**6+dx1*(420.0d0*dx2**5+dx1*(2450.0d0*dx2**4+13.0d0*dx1*(490.0d0*dx2**3+dx1* &
(588.0d0*dx2**2+11.0d0*dx1*(28.0d0*dx2+5.0d0*dx1)))))))
      ci(19)=28.0d0*(dx2**6+2.0d0*dx1*(16.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(560.0d0*dx2**3+13.0d0*dx1*(70.0d0*dx2**2+dx1* &
(48.0d0*dx2+11.0d0*dx1))))))
      ci(20)=56.0d0*(dx2**5+2.0d0*dx1*(10.0d0*dx2**4+dx1*(60.0d0*dx2**3+dx1*(140.0d0*dx2**2+13.0d0*dx1*(10.0d0*dx2+3.0d0*dx1)))))
      ci(21)=14.0d0*(5.0d0*dx2**4+2.0d0*dx1*(32.0d0*dx2**3+5.0d0*dx1*(24.0d0*dx2**2+dx1*(32.0d0*dx2+13.0d0*dx1))))
      ci(22)=8.0d0*(7.0d0*dx2**3+2.0d0*dx1*(28.0d0*dx2**2+5.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)))
      ci(23)=4.0d0*(7.0d0*dx2**2+2.0d0*dx1*(16.0d0*dx2+15.0d0*dx1))
      ci(24)=8.0d0*(dx2+2.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>8, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^17 * (x-x2)^n2 as sum_k=0^(17+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_17(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**17
      ci(2)=17.0d0*dx1**16
      ci(3)=136.0d0*dx1**15
      ci(4)=680.0d0*dx1**14
      ci(5)=2380.0d0*dx1**13
      ci(6)=6188.0d0*dx1**12
      ci(7)=12376.0d0*dx1**11
      ci(8)=19448.0d0*dx1**10
      ci(9)=24310.0d0*dx1**9
      ci(10)=24310.0d0*dx1**8
      ci(11)=19448.0d0*dx1**7
      ci(12)=12376.0d0*dx1**6
      ci(13)=6188.0d0*dx1**5
      ci(14)=2380.0d0*dx1**4
      ci(15)=680.0d0*dx1**3
      ci(16)=136.0d0*dx1**2
      ci(17)=17.0d0*dx1
      ci(18)=1.0d0
    case(1)
      ci(1)=dx1**17*dx2
      ci(2)=dx1**16*(17.0d0*dx2+dx1)
      ci(3)=17.0d0*dx1**15*(8.0d0*dx2+dx1)
      ci(4)=136.0d0*dx1**14*(5.0d0*dx2+dx1)
      ci(5)=340.0d0*dx1**13*(7.0d0*dx2+2.0d0*dx1)
      ci(6)=476.0d0*dx1**12*(13.0d0*dx2+5.0d0*dx1)
      ci(7)=6188.0d0*dx1**11*(2.0d0*dx2+dx1)
      ci(8)=1768.0d0*dx1**10*(11.0d0*dx2+7.0d0*dx1)
      ci(9)=4862.0d0*dx1**9*(5.0d0*dx2+4.0d0*dx1)
      ci(10)=24310.0d0*dx1**8*(dx2+dx1)
      ci(11)=4862.0d0*dx1**7*(4.0d0*dx2+5.0d0*dx1)
      ci(12)=1768.0d0*dx1**6*(7.0d0*dx2+11.0d0*dx1)
      ci(13)=6188.0d0*dx1**5*(dx2+2.0d0*dx1)
      ci(14)=476.0d0*dx1**4*(5.0d0*dx2+13.0d0*dx1)
      ci(15)=340.0d0*dx1**3*(2.0d0*dx2+7.0d0*dx1)
      ci(16)=136.0d0*dx1**2*(dx2+5.0d0*dx1)
      ci(17)=17.0d0*dx1*(dx2+8.0d0*dx1)
      ci(18)=dx2+17.0d0*dx1
      ci(19)=1.0d0
    case(2)
      ci(1)=dx1**17*dx2**2
      ci(2)=dx1**16*dx2*(17.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**15*(136.0d0*dx2**2+dx1*(34.0d0*dx2+dx1))
      ci(4)=17.0d0*dx1**14*(40.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))
      ci(5)=68.0d0*dx1**13*(35.0d0*dx2**2+2.0d0*dx1*(10.0d0*dx2+dx1))
      ci(6)=68.0d0*dx1**12*(91.0d0*dx2**2+10.0d0*dx1*(7.0d0*dx2+dx1))
      ci(7)=476.0d0*dx1**11*(26.0d0*dx2**2+dx1*(26.0d0*dx2+5.0d0*dx1))
      ci(8)=884.0d0*dx1**10*(22.0d0*dx2**2+7.0d0*dx1*(4.0d0*dx2+dx1))
      ci(9)=442.0d0*dx1**9*(55.0d0*dx2**2+4.0d0*dx1*(22.0d0*dx2+7.0d0*dx1))
      ci(10)=4862.0d0*dx1**8*(5.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+2.0d0*dx1))
      ci(11)=4862.0d0*dx1**7*(4.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))
      ci(12)=442.0d0*dx1**6*(28.0d0*dx2**2+11.0d0*dx1*(8.0d0*dx2+5.0d0*dx1))
      ci(13)=884.0d0*dx1**5*(7.0d0*dx2**2+2.0d0*dx1*(14.0d0*dx2+11.0d0*dx1))
      ci(14)=476.0d0*dx1**4*(5.0d0*dx2**2+26.0d0*dx1*(dx2+dx1))
      ci(15)=68.0d0*dx1**3*(10.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2+13.0d0*dx1))
      ci(16)=68.0d0*dx1**2*(2.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+7.0d0*dx1))
      ci(17)=17.0d0*dx1*(dx2**2+8.0d0*dx1*(2.0d0*dx2+5.0d0*dx1))
      ci(18)=dx2**2+34.0d0*dx1*(dx2+4.0d0*dx1)
      ci(19)=2.0d0*dx2+17.0d0*dx1
      ci(20)=1.0d0
    case(3)
      ci(1)=dx1**17*dx2**3
      ci(2)=dx1**16*dx2**2*(17.0d0*dx2+3.0d0*dx1)
      ci(3)=dx1**15*dx2*(136.0d0*dx2**2+3.0d0*dx1*(17.0d0*dx2+dx1))
      ci(4)=dx1**14*(680.0d0*dx2**3+dx1*(408.0d0*dx2**2+dx1*(51.0d0*dx2+dx1)))
      ci(5)=17.0d0*dx1**13*(140.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(24.0d0*dx2+dx1)))
      ci(6)=68.0d0*dx1**12*(91.0d0*dx2**3+dx1*(105.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+dx1)))
      ci(7)=68.0d0*dx1**11*(182.0d0*dx2**3+dx1*(273.0d0*dx2**2+5.0d0*dx1*(21.0d0*dx2+2.0d0*dx1)))
      ci(8)=68.0d0*dx1**10*(286.0d0*dx2**3+7.0d0*dx1*(78.0d0*dx2**2+dx1*(39.0d0*dx2+5.0d0*dx1)))
      ci(9)=442.0d0*dx1**9*(55.0d0*dx2**3+2.0d0*dx1*(66.0d0*dx2**2+7.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(10)=442.0d0*dx1**8*(55.0d0*dx2**3+dx1*(165.0d0*dx2**2+4.0d0*dx1*(33.0d0*dx2+7.0d0*dx1)))
      ci(11)=4862.0d0*dx1**7*(4.0d0*dx2**3+dx1*(15.0d0*dx2**2+dx1*(15.0d0*dx2+4.0d0*dx1)))
      ci(12)=442.0d0*dx1**6*(28.0d0*dx2**3+11.0d0*dx1*(12.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(13)=442.0d0*dx1**5*(14.0d0*dx2**3+dx1*(84.0d0*dx2**2+11.0d0*dx1*(12.0d0*dx2+5.0d0*dx1)))
      ci(14)=68.0d0*dx1**4*(35.0d0*dx2**3+13.0d0*dx1*(21.0d0*dx2**2+2.0d0*dx1*(21.0d0*dx2+11.0d0*dx1)))
      ci(15)=68.0d0*dx1**3*(10.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(16)=68.0d0*dx1**2*(2.0d0*dx2**3+dx1*(30.0d0*dx2**2+7.0d0*dx1*(15.0d0*dx2+13.0d0*dx1)))
      ci(17)=17.0d0*dx1*(dx2**3+4.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(6.0d0*dx2+7.0d0*dx1)))
      ci(18)=dx2**3+17.0d0*dx1*(3.0d0*dx2**2+8.0d0*dx1*(3.0d0*dx2+5.0d0*dx1))
      ci(19)=3.0d0*dx2**2+17.0d0*dx1*(3.0d0*dx2+8.0d0*dx1)
      ci(20)=3.0d0*dx2+17.0d0*dx1
      ci(21)=1.0d0
    case(4)
      ci(1)=dx1**17*dx2**4
      ci(2)=dx1**16*dx2**3*(17.0d0*dx2+4.0d0*dx1)
      ci(3)=2.0d0*dx1**15*dx2**2*(68.0d0*dx2**2+dx1*(34.0d0*dx2+3.0d0*dx1))
      ci(4)=2.0d0*dx1**14*dx2*(340.0d0*dx2**3+dx1*(272.0d0*dx2**2+dx1*(51.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**13*(2380.0d0*dx2**4+dx1*(2720.0d0*dx2**3+dx1*(816.0d0*dx2**2+dx1*(68.0d0*dx2+dx1))))
      ci(6)=17.0d0*dx1**12*(364.0d0*dx2**4+dx1*(560.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(32.0d0*dx2+dx1))))
      ci(7)=136.0d0*dx1**11*(91.0d0*dx2**4+dx1*(182.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(20.0d0*dx2+dx1))))
      ci(8)=136.0d0*dx1**10*(143.0d0*dx2**4+dx1*(364.0d0*dx2**3+dx1*(273.0d0*dx2**2+5.0d0*dx1*(14.0d0*dx2+dx1))))
      ci(9)=34.0d0*dx1**9*(715.0d0*dx2**4+2.0d0*dx1*(1144.0d0*dx2**3+7.0d0*dx1*(156.0d0*dx2**2+dx1*(52.0d0*dx2+5.0d0*dx1))))
      ci(10)=442.0d0*dx1**8*(55.0d0*dx2**4+2.0d0*dx1*(110.0d0*dx2**3+dx1*(132.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(11)=884.0d0*dx1**7*(22.0d0*dx2**4+dx1*(110.0d0*dx2**3+dx1*(165.0d0*dx2**2+2.0d0*dx1*(44.0d0*dx2+7.0d0*dx1))))
      ci(12)=884.0d0*dx1**6*(14.0d0*dx2**4+11.0d0*dx1*(8.0d0*dx2**3+dx1*(15.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+dx1))))
      ci(13)=442.0d0*dx1**5*(14.0d0*dx2**4+dx1*(112.0d0*dx2**3+11.0d0*dx1*(24.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(14)=34.0d0*dx1**4*(70.0d0*dx2**4+13.0d0*dx1*(56.0d0*dx2**3+dx1*(168.0d0*dx2**2+11.0d0*dx1*(16.0d0*dx2+5.0d0*dx1))))
      ci(15)=136.0d0*dx1**3*(5.0d0*dx2**4+dx1*(70.0d0*dx2**3+13.0d0*dx1*(21.0d0*dx2**2+dx1*(28.0d0*dx2+11.0d0*dx1))))
      ci(16)=136.0d0*dx1**2*(dx2**4+dx1*(20.0d0*dx2**3+7.0d0*dx1*(15.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(17)=17.0d0*dx1*(dx2**4+4.0d0*dx1*(8.0d0*dx2**3+dx1*(60.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+13.0d0*dx1))))
      ci(18)=dx2**4+68.0d0*dx1*(dx2**3+dx1*(12.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+7.0d0*dx1)))
      ci(19)=2.0d0*(2.0d0*dx2**3+17.0d0*dx1*(3.0d0*dx2**2+4.0d0*dx1*(4.0d0*dx2+5.0d0*dx1)))
      ci(20)=2.0d0*(3.0d0*dx2**2+34.0d0*dx1*(dx2+2.0d0*dx1))
      ci(21)=4.0d0*dx2+17.0d0*dx1
      ci(22)=1.0d0
    case(5)
      ci(1)=dx1**17*dx2**5
      ci(2)=dx1**16*dx2**4*(17.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**15*dx2**3*(136.0d0*dx2**2+5.0d0*dx1*(17.0d0*dx2+2.0d0*dx1))
      ci(4)=10.0d0*dx1**14*dx2**2*(68.0d0*dx2**3+dx1*(68.0d0*dx2**2+dx1*(17.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**13*dx2*(476.0d0*dx2**4+dx1*(680.0d0*dx2**3+dx1*(272.0d0*dx2**2+dx1*(34.0d0*dx2+dx1))))
      ci(6)=dx1**12*(6188.0d0*dx2**5+dx1*(11900.0d0*dx2**4+dx1*(6800.0d0*dx2**3+dx1*(1360.0d0*dx2**2+dx1*(85.0d0*dx2+dx1)))))
      ci(7)=17.0d0*dx1**11*(728.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1400.0d0*dx2**3+dx1*(400.0d0*dx2**2+dx1*(40.0d0*dx2+dx1)))))
      ci(8)=136.0d0*dx1**10*(143.0d0*dx2**5+dx1*(455.0d0*dx2**4+dx1*(455.0d0*dx2**3+dx1*(175.0d0*dx2**2+dx1*(25.0d0*dx2+dx1)))))
      ci(9)=170.0d0*dx1**9*(143.0d0*dx2**5+2.0d0*dx1*(286.0d0*dx2**4+dx1*(364.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(35.0d0*dx2 &
+2.0d0*dx1)))))
      ci(10)=170.0d0*dx1**8*(143.0d0*dx2**5+dx1*(715.0d0*dx2**4+2.0d0*dx1*(572.0d0*dx2**3+7.0d0*dx1*(52.0d0*dx2**2+dx1*(13.0d0*dx2 &
+dx1)))))
      ci(11)=442.0d0*dx1**7*(44.0d0*dx2**5+dx1*(275.0d0*dx2**4+2.0d0*dx1*(275.0d0*dx2**3+dx1*(220.0d0*dx2**2+7.0d0*dx1*(10.0d0*dx2 &
+dx1)))))
      ci(12)=884.0d0*dx1**6*(14.0d0*dx2**5+dx1*(110.0d0*dx2**4+dx1*(275.0d0*dx2**3+dx1*(275.0d0*dx2**2+2.0d0*dx1*(55.0d0*dx2 &
+7.0d0*dx1)))))
      ci(13)=442.0d0*dx1**5*(14.0d0*dx2**5+dx1*(140.0d0*dx2**4+11.0d0*dx1*(40.0d0*dx2**3+dx1*(50.0d0*dx2**2+dx1*(25.0d0*dx2 &
+4.0d0*dx1)))))
      ci(14)=170.0d0*dx1**4*(14.0d0*dx2**5+13.0d0*dx1*(14.0d0*dx2**4+dx1*(56.0d0*dx2**3+11.0d0*dx1*(8.0d0*dx2**2+dx1*(5.0d0*dx2 &
+dx1)))))
      ci(15)=170.0d0*dx1**3*(4.0d0*dx2**5+dx1*(70.0d0*dx2**4+13.0d0*dx1*(28.0d0*dx2**3+dx1*(56.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2 &
+dx1)))))
      ci(16)=136.0d0*dx1**2*(dx2**5+dx1*(25.0d0*dx2**4+dx1*(175.0d0*dx2**3+13.0d0*dx1*(35.0d0*dx2**2+dx1*(35.0d0*dx2+11.0d0*dx1))) &
))
      ci(17)=17.0d0*dx1*(dx2**5+4.0d0*dx1*(10.0d0*dx2**4+dx1*(100.0d0*dx2**3+7.0d0*dx1*(50.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2 &
+2.0d0*dx1)))))
      ci(18)=dx2**5+17.0d0*dx1*(5.0d0*dx2**4+4.0d0*dx1*(20.0d0*dx2**3+dx1*(100.0d0*dx2**2+7.0d0*dx1*(25.0d0*dx2+13.0d0*dx1))))
      ci(19)=5.0d0*(dx2**4+34.0d0*dx1*(dx2**3+2.0d0*dx1*(4.0d0*dx2**2+dx1*(10.0d0*dx2+7.0d0*dx1))))
      ci(20)=10.0d0*(dx2**3+17.0d0*dx1*(dx2**2+4.0d0*dx1*(dx2+dx1)))
      ci(21)=10.0d0*dx2**2+17.0d0*dx1*(5.0d0*dx2+8.0d0*dx1)
      ci(22)=5.0d0*dx2+17.0d0*dx1
      ci(23)=1.0d0
    case(6)
      ci(1)=dx1**17*dx2**6
      ci(2)=dx1**16*dx2**5*(17.0d0*dx2+6.0d0*dx1)
      ci(3)=dx1**15*dx2**4*(136.0d0*dx2**2+3.0d0*dx1*(34.0d0*dx2+5.0d0*dx1))
      ci(4)=dx1**14*dx2**3*(680.0d0*dx2**3+dx1*(816.0d0*dx2**2+5.0d0*dx1*(51.0d0*dx2+4.0d0*dx1)))
      ci(5)=5.0d0*dx1**13*dx2**2*(476.0d0*dx2**4+dx1*(816.0d0*dx2**3+dx1*(408.0d0*dx2**2+dx1*(68.0d0*dx2+3.0d0*dx1))))
      ci(6)=dx1**12*dx2*(6188.0d0*dx2**5+dx1*(14280.0d0*dx2**4+dx1*(10200.0d0*dx2**3+dx1*(2720.0d0*dx2**2+3.0d0*dx1*(85.0d0*dx2 &
+2.0d0*dx1)))))
      ci(7)=dx1**11*(12376.0d0*dx2**6+dx1*(37128.0d0*dx2**5+dx1*(35700.0d0*dx2**4+dx1*(13600.0d0*dx2**3+dx1*(2040.0d0*dx2**2+dx1* &
(102.0d0*dx2+dx1))))))
      ci(8)=17.0d0*dx1**10*(1144.0d0*dx2**6+dx1*(4368.0d0*dx2**5+dx1*(5460.0d0*dx2**4+dx1*(2800.0d0*dx2**3+dx1*(600.0d0*dx2**2 &
+dx1*(48.0d0*dx2+dx1))))))
      ci(9)=34.0d0*dx1**9*(715.0d0*dx2**6+2.0d0*dx1*(1716.0d0*dx2**5+dx1*(2730.0d0*dx2**4+dx1*(1820.0d0*dx2**3+dx1*(525.0d0*dx2**2 &
+2.0d0*dx1*(30.0d0*dx2+dx1))))))
      ci(10)=170.0d0*dx1**8*(143.0d0*dx2**6+2.0d0*dx1*(429.0d0*dx2**5+dx1*(858.0d0*dx2**4+dx1*(728.0d0*dx2**3+dx1*(273.0d0*dx2**2 &
+2.0d0*dx1*(21.0d0*dx2+dx1))))))
      ci(11)=34.0d0*dx1**7*(572.0d0*dx2**6+dx1*(4290.0d0*dx2**5+dx1*(10725.0d0*dx2**4+2.0d0*dx1*(5720.0d0*dx2**3+7.0d0*dx1* &
(390.0d0*dx2**2+dx1*(78.0d0*dx2+5.0d0*dx1))))))
      ci(12)=442.0d0*dx1**6*(28.0d0*dx2**6+dx1*(264.0d0*dx2**5+dx1*(825.0d0*dx2**4+2.0d0*dx1*(550.0d0*dx2**3+dx1*(330.0d0*dx2**2 &
+7.0d0*dx1*(12.0d0*dx2+dx1))))))
      ci(13)=442.0d0*dx1**5*(14.0d0*dx2**6+dx1*(168.0d0*dx2**5+dx1*(660.0d0*dx2**4+dx1*(1100.0d0*dx2**3+dx1*(825.0d0*dx2**2 &
+4.0d0*dx1*(66.0d0*dx2+7.0d0*dx1))))))
      ci(14)=34.0d0*dx1**4*(70.0d0*dx2**6+13.0d0*dx1*(84.0d0*dx2**5+dx1*(420.0d0*dx2**4+11.0d0*dx1*(80.0d0*dx2**3+dx1* &
(75.0d0*dx2**2+2.0d0*dx1*(15.0d0*dx2+2.0d0*dx1))))))
      ci(15)=170.0d0*dx1**3*(4.0d0*dx2**6+dx1*(84.0d0*dx2**5+13.0d0*dx1*(42.0d0*dx2**4+dx1*(112.0d0*dx2**3+11.0d0*dx1* &
(12.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))))
      ci(16)=34.0d0*dx1**2*(4.0d0*dx2**6+dx1*(120.0d0*dx2**5+dx1*(1050.0d0*dx2**4+13.0d0*dx1*(280.0d0*dx2**3+dx1*(420.0d0*dx2**2 &
+11.0d0*dx1*(24.0d0*dx2+5.0d0*dx1))))))
      ci(17)=17.0d0*dx1*(dx2**6+4.0d0*dx1*(12.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(700.0d0*dx2**3+13.0d0*dx1*(105.0d0*dx2**2 &
+2.0d0*dx1*(42.0d0*dx2+11.0d0*dx1))))))
      ci(18)=dx2**6+34.0d0*dx1*(3.0d0*dx2**5+2.0d0*dx1*(30.0d0*dx2**4+dx1*(200.0d0*dx2**3+7.0d0*dx1*(75.0d0*dx2**2+26.0d0*dx1* &
(3.0d0*dx2+dx1)))))
      ci(19)=6.0d0*dx2**5+17.0d0*dx1*(15.0d0*dx2**4+4.0d0*dx1*(40.0d0*dx2**3+dx1*(150.0d0*dx2**2+7.0d0*dx1*(30.0d0*dx2+13.0d0*dx1) &
)))
      ci(20)=5.0d0*(3.0d0*dx2**4+68.0d0*dx1*(dx2**3+dx1*(6.0d0*dx2**2+dx1*(12.0d0*dx2+7.0d0*dx1))))
      ci(21)=20.0d0*dx2**3+17.0d0*dx1*(15.0d0*dx2**2+8.0d0*dx1*(6.0d0*dx2+5.0d0*dx1))
      ci(22)=15.0d0*dx2**2+34.0d0*dx1*(3.0d0*dx2+4.0d0*dx1)
      ci(23)=6.0d0*dx2+17.0d0*dx1
      ci(24)=1.0d0
    case(7)
      ci(1)=dx1**17*dx2**7
      ci(2)=dx1**16*dx2**6*(17.0d0*dx2+7.0d0*dx1)
      ci(3)=dx1**15*dx2**5*(136.0d0*dx2**2+7.0d0*dx1*(17.0d0*dx2+3.0d0*dx1))
      ci(4)=dx1**14*dx2**4*(680.0d0*dx2**3+7.0d0*dx1*(136.0d0*dx2**2+dx1*(51.0d0*dx2+5.0d0*dx1)))
      ci(5)=7.0d0*dx1**13*dx2**3*(340.0d0*dx2**4+dx1*(680.0d0*dx2**3+dx1*(408.0d0*dx2**2+5.0d0*dx1*(17.0d0*dx2+dx1))))
      ci(6)=7.0d0*dx1**12*dx2**2*(884.0d0*dx2**5+dx1*(2380.0d0*dx2**4+dx1*(2040.0d0*dx2**3+dx1*(680.0d0*dx2**2+dx1*(85.0d0*dx2 &
+3.0d0*dx1)))))
      ci(7)=7.0d0*dx1**11*dx2*(1768.0d0*dx2**6+dx1*(6188.0d0*dx2**5+dx1*(7140.0d0*dx2**4+dx1*(3400.0d0*dx2**3+dx1*(680.0d0*dx2**2 &
+dx1*(51.0d0*dx2+dx1))))))
      ci(8)=dx1**10*(19448.0d0*dx2**7+dx1*(86632.0d0*dx2**6+dx1*(129948.0d0*dx2**5+dx1*(83300.0d0*dx2**4+dx1*(23800.0d0*dx2**3 &
+dx1*(2856.0d0*dx2**2+dx1*(119.0d0*dx2+dx1)))))))
      ci(9)=17.0d0*dx1**9*(1430.0d0*dx2**7+dx1*(8008.0d0*dx2**6+dx1*(15288.0d0*dx2**5+dx1*(12740.0d0*dx2**4+dx1*(4900.0d0*dx2**3 &
+dx1*(840.0d0*dx2**2+dx1*(56.0d0*dx2+dx1)))))))
      ci(10)=34.0d0*dx1**8*(715.0d0*dx2**7+dx1*(5005.0d0*dx2**6+2.0d0*dx1*(6006.0d0*dx2**5+dx1*(6370.0d0*dx2**4+dx1* &
(3185.0d0*dx2**3+dx1*(735.0d0*dx2**2+2.0d0*dx1*(35.0d0*dx2+dx1)))))))
      ci(11)=34.0d0*dx1**7*(572.0d0*dx2**7+dx1*(5005.0d0*dx2**6+dx1*(15015.0d0*dx2**5+2.0d0*dx1*(10010.0d0*dx2**4+dx1* &
(6370.0d0*dx2**3+dx1*(1911.0d0*dx2**2+5.0d0*dx1*(49.0d0*dx2+2.0d0*dx1)))))))
      ci(12)=238.0d0*dx1**6*(52.0d0*dx2**7+dx1*(572.0d0*dx2**6+dx1*(2145.0d0*dx2**5+dx1*(3575.0d0*dx2**4+2.0d0*dx1* &
(1430.0d0*dx2**3+dx1*(546.0d0*dx2**2+dx1*(91.0d0*dx2+5.0d0*dx1)))))))
      ci(13)=3094.0d0*dx1**5*(2.0d0*dx2**7+dx1*(28.0d0*dx2**6+dx1*(132.0d0*dx2**5+dx1*(275.0d0*dx2**4+dx1*(275.0d0*dx2**3 &
+2.0d0*dx1*(66.0d0*dx2**2+dx1*(14.0d0*dx2+dx1)))))))
      ci(14)=238.0d0*dx1**4*(10.0d0*dx2**7+13.0d0*dx1*(14.0d0*dx2**6+dx1*(84.0d0*dx2**5+dx1*(220.0d0*dx2**4+dx1*(275.0d0*dx2**3 &
+dx1*(165.0d0*dx2**2+4.0d0*dx1*(11.0d0*dx2+dx1)))))))
      ci(15)=34.0d0*dx1**3*(20.0d0*dx2**7+dx1*(490.0d0*dx2**6+13.0d0*dx1*(294.0d0*dx2**5+dx1*(980.0d0*dx2**4+11.0d0*dx1* &
(140.0d0*dx2**3+dx1*(105.0d0*dx2**2+dx1*(35.0d0*dx2+4.0d0*dx1)))))))
      ci(16)=34.0d0*dx1**2*(4.0d0*dx2**7+dx1*(140.0d0*dx2**6+dx1*(1470.0d0*dx2**5+13.0d0*dx1*(490.0d0*dx2**4+dx1*(980.0d0*dx2**3 &
+11.0d0*dx1*(84.0d0*dx2**2+5.0d0*dx1*(7.0d0*dx2+dx1)))))))
      ci(17)=17.0d0*dx1*(dx2**7+2.0d0*dx1*(28.0d0*dx2**6+dx1*(420.0d0*dx2**5+dx1*(2450.0d0*dx2**4+13.0d0*dx1*(490.0d0*dx2**3+dx1* &
(588.0d0*dx2**2+11.0d0*dx1*(28.0d0*dx2+5.0d0*dx1)))))))
      ci(18)=dx2**7+17.0d0*dx1*(7.0d0*dx2**6+4.0d0*dx1*(42.0d0*dx2**5+dx1*(350.0d0*dx2**4+dx1*(1225.0d0*dx2**3+13.0d0*dx1* &
(147.0d0*dx2**2+2.0d0*dx1*(49.0d0*dx2+11.0d0*dx1))))))
      ci(19)=7.0d0*(dx2**6+17.0d0*dx1*(3.0d0*dx2**5+4.0d0*dx1*(10.0d0*dx2**4+dx1*(50.0d0*dx2**3+dx1*(105.0d0*dx2**2+13.0d0*dx1* &
(7.0d0*dx2+2.0d0*dx1))))))
      ci(20)=7.0d0*(3.0d0*dx2**5+17.0d0*dx1*(5.0d0*dx2**4+4.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(35.0d0*dx2+13.0d0*dx1)) &
)))
      ci(21)=7.0d0*(5.0d0*dx2**4+17.0d0*dx1*(5.0d0*dx2**3+4.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(22)=35.0d0*dx2**3+17.0d0*dx1*(21.0d0*dx2**2+8.0d0*dx1*(7.0d0*dx2+5.0d0*dx1))
      ci(23)=21.0d0*dx2**2+17.0d0*dx1*(7.0d0*dx2+8.0d0*dx1)
      ci(24)=7.0d0*dx2+17.0d0*dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>7, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^18 * (x-x2)^n2 as sum_k=0^(18+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_18(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**18
      ci(2)=18.0d0*dx1**17
      ci(3)=153.0d0*dx1**16
      ci(4)=816.0d0*dx1**15
      ci(5)=3060.0d0*dx1**14
      ci(6)=8568.0d0*dx1**13
      ci(7)=18564.0d0*dx1**12
      ci(8)=31824.0d0*dx1**11
      ci(9)=43758.0d0*dx1**10
      ci(10)=48620.0d0*dx1**9
      ci(11)=43758.0d0*dx1**8
      ci(12)=31824.0d0*dx1**7
      ci(13)=18564.0d0*dx1**6
      ci(14)=8568.0d0*dx1**5
      ci(15)=3060.0d0*dx1**4
      ci(16)=816.0d0*dx1**3
      ci(17)=153.0d0*dx1**2
      ci(18)=18.0d0*dx1
      ci(19)=1.0d0
    case(1)
      ci(1)=dx1**18*dx2
      ci(2)=dx1**17*(18.0d0*dx2+dx1)
      ci(3)=9.0d0*dx1**16*(17.0d0*dx2+2.0d0*dx1)
      ci(4)=51.0d0*dx1**15*(16.0d0*dx2+3.0d0*dx1)
      ci(5)=204.0d0*dx1**14*(15.0d0*dx2+4.0d0*dx1)
      ci(6)=612.0d0*dx1**13*(14.0d0*dx2+5.0d0*dx1)
      ci(7)=1428.0d0*dx1**12*(13.0d0*dx2+6.0d0*dx1)
      ci(8)=2652.0d0*dx1**11*(12.0d0*dx2+7.0d0*dx1)
      ci(9)=3978.0d0*dx1**10*(11.0d0*dx2+8.0d0*dx1)
      ci(10)=4862.0d0*dx1**9*(10.0d0*dx2+9.0d0*dx1)
      ci(11)=4862.0d0*dx1**8*(9.0d0*dx2+10.0d0*dx1)
      ci(12)=3978.0d0*dx1**7*(8.0d0*dx2+11.0d0*dx1)
      ci(13)=2652.0d0*dx1**6*(7.0d0*dx2+12.0d0*dx1)
      ci(14)=1428.0d0*dx1**5*(6.0d0*dx2+13.0d0*dx1)
      ci(15)=612.0d0*dx1**4*(5.0d0*dx2+14.0d0*dx1)
      ci(16)=204.0d0*dx1**3*(4.0d0*dx2+15.0d0*dx1)
      ci(17)=51.0d0*dx1**2*(3.0d0*dx2+16.0d0*dx1)
      ci(18)=9.0d0*dx1*(2.0d0*dx2+17.0d0*dx1)
      ci(19)=dx2+18.0d0*dx1
      ci(20)=1.0d0
    case(2)
      ci(1)=dx1**18*dx2**2
      ci(2)=2.0d0*dx1**17*dx2*(9.0d0*dx2+dx1)
      ci(3)=dx1**16*(153.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))
      ci(4)=6.0d0*dx1**15*(136.0d0*dx2**2+3.0d0*dx1*(17.0d0*dx2+dx1))
      ci(5)=51.0d0*dx1**14*(60.0d0*dx2**2+dx1*(32.0d0*dx2+3.0d0*dx1))
      ci(6)=408.0d0*dx1**13*(21.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1))
      ci(7)=204.0d0*dx1**12*(91.0d0*dx2**2+3.0d0*dx1*(28.0d0*dx2+5.0d0*dx1))
      ci(8)=408.0d0*dx1**11*(78.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+3.0d0*dx1))
      ci(9)=1326.0d0*dx1**10*(33.0d0*dx2**2+2.0d0*dx1*(24.0d0*dx2+7.0d0*dx1))
      ci(10)=884.0d0*dx1**9*(55.0d0*dx2**2+9.0d0*dx1*(11.0d0*dx2+4.0d0*dx1))
      ci(11)=4862.0d0*dx1**8*(9.0d0*dx2**2+dx1*(20.0d0*dx2+9.0d0*dx1))
      ci(12)=884.0d0*dx1**7*(36.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2+5.0d0*dx1))
      ci(13)=1326.0d0*dx1**6*(14.0d0*dx2**2+3.0d0*dx1*(16.0d0*dx2+11.0d0*dx1))
      ci(14)=408.0d0*dx1**5*(21.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+6.0d0*dx1))
      ci(15)=204.0d0*dx1**4*(15.0d0*dx2**2+7.0d0*dx1*(12.0d0*dx2+13.0d0*dx1))
      ci(16)=408.0d0*dx1**3*(2.0d0*dx2**2+3.0d0*dx1*(5.0d0*dx2+7.0d0*dx1))
      ci(17)=51.0d0*dx1**2*(3.0d0*dx2**2+4.0d0*dx1*(8.0d0*dx2+15.0d0*dx1))
      ci(18)=6.0d0*dx1*(3.0d0*dx2**2+17.0d0*dx1*(3.0d0*dx2+8.0d0*dx1))
      ci(19)=dx2**2+9.0d0*dx1*(4.0d0*dx2+17.0d0*dx1)
      ci(20)=2.0d0*(dx2+9.0d0*dx1)
      ci(21)=1.0d0
    case(3)
      ci(1)=dx1**18*dx2**3
      ci(2)=3.0d0*dx1**17*dx2**2*(6.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**16*dx2*(51.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))
      ci(4)=dx1**15*(816.0d0*dx2**3+dx1*(459.0d0*dx2**2+dx1*(54.0d0*dx2+dx1)))
      ci(5)=9.0d0*dx1**14*(340.0d0*dx2**3+dx1*(272.0d0*dx2**2+dx1*(51.0d0*dx2+2.0d0*dx1)))
      ci(6)=153.0d0*dx1**13*(56.0d0*dx2**3+dx1*(60.0d0*dx2**2+dx1*(16.0d0*dx2+dx1)))
      ci(7)=204.0d0*dx1**12*(91.0d0*dx2**3+dx1*(126.0d0*dx2**2+dx1*(45.0d0*dx2+4.0d0*dx1)))
      ci(8)=612.0d0*dx1**11*(52.0d0*dx2**3+dx1*(91.0d0*dx2**2+dx1*(42.0d0*dx2+5.0d0*dx1)))
      ci(9)=306.0d0*dx1**10*(143.0d0*dx2**3+2.0d0*dx1*(156.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+2.0d0*dx1)))
      ci(10)=442.0d0*dx1**9*(110.0d0*dx2**3+3.0d0*dx1*(99.0d0*dx2**2+2.0d0*dx1*(36.0d0*dx2+7.0d0*dx1)))
      ci(11)=1326.0d0*dx1**8*(33.0d0*dx2**3+dx1*(110.0d0*dx2**2+3.0d0*dx1*(33.0d0*dx2+8.0d0*dx1)))
      ci(12)=1326.0d0*dx1**7*(24.0d0*dx2**3+11.0d0*dx1*(9.0d0*dx2**2+dx1*(10.0d0*dx2+3.0d0*dx1)))
      ci(13)=442.0d0*dx1**6*(42.0d0*dx2**3+dx1*(216.0d0*dx2**2+11.0d0*dx1*(27.0d0*dx2+10.0d0*dx1)))
      ci(14)=306.0d0*dx1**5*(28.0d0*dx2**3+13.0d0*dx1*(14.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1)))
      ci(15)=612.0d0*dx1**4*(5.0d0*dx2**3+dx1*(42.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+4.0d0*dx1)))
      ci(16)=204.0d0*dx1**3*(4.0d0*dx2**3+dx1*(45.0d0*dx2**2+7.0d0*dx1*(18.0d0*dx2+13.0d0*dx1)))
      ci(17)=153.0d0*dx1**2*(dx2**3+4.0d0*dx1*(4.0d0*dx2**2+dx1*(15.0d0*dx2+14.0d0*dx1)))
      ci(18)=9.0d0*dx1*(2.0d0*dx2**3+17.0d0*dx1*(3.0d0*dx2**2+4.0d0*dx1*(4.0d0*dx2+5.0d0*dx1)))
      ci(19)=dx2**3+3.0d0*dx1*(18.0d0*dx2**2+17.0d0*dx1*(9.0d0*dx2+16.0d0*dx1))
      ci(20)=3.0d0*(dx2**2+3.0d0*dx1*(6.0d0*dx2+17.0d0*dx1))
      ci(21)=3.0d0*(dx2+6.0d0*dx1)
      ci(22)=1.0d0
    case(4)
      ci(1)=dx1**18*dx2**4
      ci(2)=2.0d0*dx1**17*dx2**3*(9.0d0*dx2+2.0d0*dx1)
      ci(3)=3.0d0*dx1**16*dx2**2*(51.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+dx1))
      ci(4)=4.0d0*dx1**15*dx2*(204.0d0*dx2**3+dx1*(153.0d0*dx2**2+dx1*(27.0d0*dx2+dx1)))
      ci(5)=dx1**14*(3060.0d0*dx2**4+dx1*(3264.0d0*dx2**3+dx1*(918.0d0*dx2**2+dx1*(72.0d0*dx2+dx1))))
      ci(6)=18.0d0*dx1**13*(476.0d0*dx2**4+dx1*(680.0d0*dx2**3+dx1*(272.0d0*dx2**2+dx1*(34.0d0*dx2+dx1))))
      ci(7)=51.0d0*dx1**12*(364.0d0*dx2**4+dx1*(672.0d0*dx2**3+dx1*(360.0d0*dx2**2+dx1*(64.0d0*dx2+3.0d0*dx1))))
      ci(8)=816.0d0*dx1**11*(39.0d0*dx2**4+dx1*(91.0d0*dx2**3+dx1*(63.0d0*dx2**2+dx1*(15.0d0*dx2+dx1))))
      ci(9)=306.0d0*dx1**10*(143.0d0*dx2**4+2.0d0*dx1*(208.0d0*dx2**3+dx1*(182.0d0*dx2**2+dx1*(56.0d0*dx2+5.0d0*dx1))))
      ci(10)=68.0d0*dx1**9*(715.0d0*dx2**4+6.0d0*dx1*(429.0d0*dx2**3+dx1*(468.0d0*dx2**2+7.0d0*dx1*(26.0d0*dx2+3.0d0*dx1))))
      ci(11)=442.0d0*dx1**8*(99.0d0*dx2**4+2.0d0*dx1*(220.0d0*dx2**3+3.0d0*dx1*(99.0d0*dx2**2+dx1*(48.0d0*dx2+7.0d0*dx1))))
      ci(12)=5304.0d0*dx1**7*(6.0d0*dx2**4+dx1*(33.0d0*dx2**3+dx1*(55.0d0*dx2**2+3.0d0*dx1*(11.0d0*dx2+2.0d0*dx1))))
      ci(13)=442.0d0*dx1**6*(42.0d0*dx2**4+dx1*(288.0d0*dx2**3+11.0d0*dx1*(54.0d0*dx2**2+dx1*(40.0d0*dx2+9.0d0*dx1))))
      ci(14)=68.0d0*dx1**5*(126.0d0*dx2**4+13.0d0*dx1*(84.0d0*dx2**3+dx1*(216.0d0*dx2**2+11.0d0*dx1*(18.0d0*dx2+5.0d0*dx1))))
      ci(15)=306.0d0*dx1**4*(10.0d0*dx2**4+dx1*(112.0d0*dx2**3+13.0d0*dx1*(28.0d0*dx2**2+dx1*(32.0d0*dx2+11.0d0*dx1))))
      ci(16)=816.0d0*dx1**3*(dx2**4+dx1*(15.0d0*dx2**3+dx1*(63.0d0*dx2**2+13.0d0*dx1*(7.0d0*dx2+3.0d0*dx1))))
      ci(17)=51.0d0*dx1**2*(3.0d0*dx2**4+4.0d0*dx1*(16.0d0*dx2**3+dx1*(90.0d0*dx2**2+7.0d0*dx1*(24.0d0*dx2+13.0d0*dx1))))
      ci(18)=18.0d0*dx1*(dx2**4+34.0d0*dx1*(dx2**3+2.0d0*dx1*(4.0d0*dx2**2+dx1*(10.0d0*dx2+7.0d0*dx1))))
      ci(19)=dx2**4+6.0d0*dx1*(12.0d0*dx2**3+17.0d0*dx1*(9.0d0*dx2**2+2.0d0*dx1*(16.0d0*dx2+15.0d0*dx1)))
      ci(20)=4.0d0*(dx2**3+3.0d0*dx1*(9.0d0*dx2**2+17.0d0*dx1*(3.0d0*dx2+4.0d0*dx1)))
      ci(21)=3.0d0*(2.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+17.0d0*dx1))
      ci(22)=2.0d0*(2.0d0*dx2+9.0d0*dx1)
      ci(23)=1.0d0
    case(5)
      ci(1)=dx1**18*dx2**5
      ci(2)=dx1**17*dx2**4*(18.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**16*dx2**3*(153.0d0*dx2**2+10.0d0*dx1*(9.0d0*dx2+dx1))
      ci(4)=dx1**15*dx2**2*(816.0d0*dx2**3+5.0d0*dx1*(153.0d0*dx2**2+2.0d0*dx1*(18.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**14*dx2*(612.0d0*dx2**4+dx1*(816.0d0*dx2**3+dx1*(306.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))
      ci(6)=dx1**13*(8568.0d0*dx2**5+dx1*(15300.0d0*dx2**4+dx1*(8160.0d0*dx2**3+dx1*(1530.0d0*dx2**2+dx1*(90.0d0*dx2+dx1)))))
      ci(7)=3.0d0*dx1**12*(6188.0d0*dx2**5+dx1*(14280.0d0*dx2**4+dx1*(10200.0d0*dx2**3+dx1*(2720.0d0*dx2**2+3.0d0*dx1*(85.0d0*dx2 &
+2.0d0*dx1)))))
      ci(8)=51.0d0*dx1**11*(624.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1680.0d0*dx2**3+dx1*(600.0d0*dx2**2+dx1*(80.0d0*dx2+3.0d0*dx1 &
)))))
      ci(9)=102.0d0*dx1**10*(429.0d0*dx2**5+2.0d0*dx1*(780.0d0*dx2**4+dx1*(910.0d0*dx2**3+dx1*(420.0d0*dx2**2+dx1*(75.0d0*dx2 &
+4.0d0*dx1)))))
      ci(10)=170.0d0*dx1**9*(286.0d0*dx2**5+3.0d0*dx1*(429.0d0*dx2**4+2.0d0*dx1*(312.0d0*dx2**3+dx1*(182.0d0*dx2**2+3.0d0*dx1* &
(14.0d0*dx2+dx1)))))
      ci(11)=34.0d0*dx1**8*(1287.0d0*dx2**5+2.0d0*dx1*(3575.0d0*dx2**4+3.0d0*dx1*(2145.0d0*dx2**3+dx1*(1560.0d0*dx2**2+7.0d0*dx1* &
(65.0d0*dx2+6.0d0*dx1)))))
      ci(12)=442.0d0*dx1**7*(72.0d0*dx2**5+dx1*(495.0d0*dx2**4+2.0d0*dx1*(550.0d0*dx2**3+3.0d0*dx1*(165.0d0*dx2**2+dx1*(60.0d0*dx2 &
+7.0d0*dx1)))))
      ci(13)=442.0d0*dx1**6*(42.0d0*dx2**5+dx1*(360.0d0*dx2**4+dx1*(990.0d0*dx2**3+dx1*(1100.0d0*dx2**2+9.0d0*dx1*(55.0d0*dx2 &
+8.0d0*dx1)))))
      ci(14)=34.0d0*dx1**5*(252.0d0*dx2**5+13.0d0*dx1*(210.0d0*dx2**4+dx1*(720.0d0*dx2**3+11.0d0*dx1*(90.0d0*dx2**2+dx1* &
(50.0d0*dx2+9.0d0*dx1)))))
      ci(15)=170.0d0*dx1**4*(18.0d0*dx2**5+dx1*(252.0d0*dx2**4+13.0d0*dx1*(84.0d0*dx2**3+dx1*(144.0d0*dx2**2+11.0d0*dx1*(9.0d0*dx2 &
+2.0d0*dx1)))))
      ci(16)=102.0d0*dx1**3*(8.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(840.0d0*dx2**3+13.0d0*dx1*(140.0d0*dx2**2+3.0d0*dx1*(40.0d0*dx2 &
+11.0d0*dx1)))))
      ci(17)=51.0d0*dx1**2*(3.0d0*dx2**5+4.0d0*dx1*(20.0d0*dx2**4+dx1*(150.0d0*dx2**3+dx1*(420.0d0*dx2**2+13.0d0*dx1*(35.0d0*dx2 &
+12.0d0*dx1)))))
      ci(18)=3.0d0*dx1*(6.0d0*dx2**5+17.0d0*dx1*(15.0d0*dx2**4+4.0d0*dx1*(40.0d0*dx2**3+dx1*(150.0d0*dx2**2+7.0d0*dx1*(30.0d0*dx2 &
+13.0d0*dx1)))))
      ci(19)=dx2**5+6.0d0*dx1*(15.0d0*dx2**4+17.0d0*dx1*(15.0d0*dx2**3+2.0d0*dx1*(40.0d0*dx2**2+3.0d0*dx1*(25.0d0*dx2+14.0d0*dx1)) &
))
      ci(20)=5.0d0*(dx2**4+6.0d0*dx1*(6.0d0*dx2**3+17.0d0*dx1*(3.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))))
      ci(21)=10.0d0*dx2**3+3.0d0*dx1*(60.0d0*dx2**2+17.0d0*dx1*(15.0d0*dx2+16.0d0*dx1))
      ci(22)=10.0d0*dx2**2+9.0d0*dx1*(10.0d0*dx2+17.0d0*dx1)
      ci(23)=5.0d0*dx2+18.0d0*dx1
      ci(24)=1.0d0
    case(6)
      ci(1)=dx1**18*dx2**6
      ci(2)=6.0d0*dx1**17*dx2**5*(3.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**16*dx2**4*(51.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))
      ci(4)=2.0d0*dx1**15*dx2**3*(408.0d0*dx2**3+dx1*(459.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+2.0d0*dx1)))
      ci(5)=3.0d0*dx1**14*dx2**2*(1020.0d0*dx2**4+dx1*(1632.0d0*dx2**3+5.0d0*dx1*(153.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(6)=6.0d0*dx1**13*dx2*(1428.0d0*dx2**5+dx1*(3060.0d0*dx2**4+dx1*(2040.0d0*dx2**3+dx1*(510.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)) &
)))
      ci(7)=dx1**12*(18564.0d0*dx2**6+dx1*(51408.0d0*dx2**5+dx1*(45900.0d0*dx2**4+dx1*(16320.0d0*dx2**3+dx1*(2295.0d0*dx2**2+dx1* &
(108.0d0*dx2+dx1))))))
      ci(8)=18.0d0*dx1**11*(1768.0d0*dx2**6+dx1*(6188.0d0*dx2**5+dx1*(7140.0d0*dx2**4+dx1*(3400.0d0*dx2**3+dx1*(680.0d0*dx2**2 &
+dx1*(51.0d0*dx2+dx1))))))
      ci(9)=153.0d0*dx1**10*(286.0d0*dx2**6+dx1*(1248.0d0*dx2**5+dx1*(1820.0d0*dx2**4+dx1*(1120.0d0*dx2**3+dx1*(300.0d0*dx2**2 &
+dx1*(32.0d0*dx2+dx1))))))
      ci(10)=68.0d0*dx1**9*(715.0d0*dx2**6+3.0d0*dx1*(1287.0d0*dx2**5+2.0d0*dx1*(1170.0d0*dx2**4+dx1*(910.0d0*dx2**3+dx1* &
(315.0d0*dx2**2+dx1*(45.0d0*dx2+2.0d0*dx1))))))
      ci(11)=102.0d0*dx1**8*(429.0d0*dx2**6+dx1*(2860.0d0*dx2**5+3.0d0*dx1*(2145.0d0*dx2**4+2.0d0*dx1*(1040.0d0*dx2**3+dx1* &
(455.0d0*dx2**2+dx1*(84.0d0*dx2+5.0d0*dx1))))))
      ci(12)=204.0d0*dx1**7*(156.0d0*dx2**6+dx1*(1287.0d0*dx2**5+dx1*(3575.0d0*dx2**4+6.0d0*dx1*(715.0d0*dx2**3+dx1* &
(390.0d0*dx2**2+7.0d0*dx1*(13.0d0*dx2+dx1))))))
      ci(13)=442.0d0*dx1**6*(42.0d0*dx2**6+dx1*(432.0d0*dx2**5+dx1*(1485.0d0*dx2**4+dx1*(2200.0d0*dx2**3+3.0d0*dx1*(495.0d0*dx2**2 &
+2.0d0*dx1*(72.0d0*dx2+7.0d0*dx1))))))
      ci(14)=204.0d0*dx1**5*(42.0d0*dx2**6+13.0d0*dx1*(42.0d0*dx2**5+dx1*(180.0d0*dx2**4+dx1*(330.0d0*dx2**3+dx1*(275.0d0*dx2**2 &
+3.0d0*dx1*(33.0d0*dx2+4.0d0*dx1))))))
      ci(15)=102.0d0*dx1**4*(30.0d0*dx2**6+dx1*(504.0d0*dx2**5+13.0d0*dx1*(210.0d0*dx2**4+dx1*(480.0d0*dx2**3+11.0d0*dx1* &
(45.0d0*dx2**2+dx1*(20.0d0*dx2+3.0d0*dx1))))))
      ci(16)=68.0d0*dx1**3*(12.0d0*dx2**6+dx1*(270.0d0*dx2**5+dx1*(1890.0d0*dx2**4+13.0d0*dx1*(420.0d0*dx2**3+dx1*(540.0d0*dx2**2 &
+11.0d0*dx1*(27.0d0*dx2+5.0d0*dx1))))))
      ci(17)=153.0d0*dx1**2*(dx2**6+2.0d0*dx1*(16.0d0*dx2**5+dx1*(150.0d0*dx2**4+dx1*(560.0d0*dx2**3+13.0d0*dx1*(70.0d0*dx2**2 &
+dx1*(48.0d0*dx2+11.0d0*dx1))))))
      ci(18)=18.0d0*dx1*(dx2**6+17.0d0*dx1*(3.0d0*dx2**5+4.0d0*dx1*(10.0d0*dx2**4+dx1*(50.0d0*dx2**3+dx1*(105.0d0*dx2**2 &
+13.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))))))
      ci(19)=dx2**6+3.0d0*dx1*(36.0d0*dx2**5+17.0d0*dx1*(45.0d0*dx2**4+4.0d0*dx1*(80.0d0*dx2**3+dx1*(225.0d0*dx2**2+7.0d0*dx1* &
(36.0d0*dx2+13.0d0*dx1)))))
      ci(20)=6.0d0*(dx2**5+3.0d0*dx1*(15.0d0*dx2**4+34.0d0*dx1*(5.0d0*dx2**3+2.0d0*dx1*(10.0d0*dx2**2+dx1*(15.0d0*dx2+7.0d0*dx1))) &
))
      ci(21)=3.0d0*(5.0d0*dx2**4+3.0d0*dx1*(40.0d0*dx2**3+17.0d0*dx1*(15.0d0*dx2**2+4.0d0*dx1*(8.0d0*dx2+5.0d0*dx1))))
      ci(22)=2.0d0*(10.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+17.0d0*dx1*(9.0d0*dx2+8.0d0*dx1)))
      ci(23)=3.0d0*(5.0d0*dx2**2+3.0d0*dx1*(12.0d0*dx2+17.0d0*dx1))
      ci(24)=6.0d0*(dx2+3.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>6, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^19 * (x-x2)^n2 as sum_k=0^(19+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_19(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**19
      ci(2)=19.0d0*dx1**18
      ci(3)=171.0d0*dx1**17
      ci(4)=969.0d0*dx1**16
      ci(5)=3876.0d0*dx1**15
      ci(6)=11628.0d0*dx1**14
      ci(7)=27132.0d0*dx1**13
      ci(8)=50388.0d0*dx1**12
      ci(9)=75582.0d0*dx1**11
      ci(10)=92378.0d0*dx1**10
      ci(11)=92378.0d0*dx1**9
      ci(12)=75582.0d0*dx1**8
      ci(13)=50388.0d0*dx1**7
      ci(14)=27132.0d0*dx1**6
      ci(15)=11628.0d0*dx1**5
      ci(16)=3876.0d0*dx1**4
      ci(17)=969.0d0*dx1**3
      ci(18)=171.0d0*dx1**2
      ci(19)=19.0d0*dx1
      ci(20)=1.0d0
    case(1)
      ci(1)=dx1**19*dx2
      ci(2)=dx1**18*(19.0d0*dx2+dx1)
      ci(3)=19.0d0*dx1**17*(9.0d0*dx2+dx1)
      ci(4)=57.0d0*dx1**16*(17.0d0*dx2+3.0d0*dx1)
      ci(5)=969.0d0*dx1**15*(4.0d0*dx2+dx1)
      ci(6)=3876.0d0*dx1**14*(3.0d0*dx2+dx1)
      ci(7)=3876.0d0*dx1**13*(7.0d0*dx2+3.0d0*dx1)
      ci(8)=3876.0d0*dx1**12*(13.0d0*dx2+7.0d0*dx1)
      ci(9)=25194.0d0*dx1**11*(3.0d0*dx2+2.0d0*dx1)
      ci(10)=8398.0d0*dx1**10*(11.0d0*dx2+9.0d0*dx1)
      ci(11)=92378.0d0*dx1**9*(dx2+dx1)
      ci(12)=8398.0d0*dx1**8*(9.0d0*dx2+11.0d0*dx1)
      ci(13)=25194.0d0*dx1**7*(2.0d0*dx2+3.0d0*dx1)
      ci(14)=3876.0d0*dx1**6*(7.0d0*dx2+13.0d0*dx1)
      ci(15)=3876.0d0*dx1**5*(3.0d0*dx2+7.0d0*dx1)
      ci(16)=3876.0d0*dx1**4*(dx2+3.0d0*dx1)
      ci(17)=969.0d0*dx1**3*(dx2+4.0d0*dx1)
      ci(18)=57.0d0*dx1**2*(3.0d0*dx2+17.0d0*dx1)
      ci(19)=19.0d0*dx1*(dx2+9.0d0*dx1)
      ci(20)=dx2+19.0d0*dx1
      ci(21)=1.0d0
    case(2)
      ci(1)=dx1**19*dx2**2
      ci(2)=dx1**18*dx2*(19.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**17*(171.0d0*dx2**2+dx1*(38.0d0*dx2+dx1))
      ci(4)=19.0d0*dx1**16*(51.0d0*dx2**2+dx1*(18.0d0*dx2+dx1))
      ci(5)=57.0d0*dx1**15*(68.0d0*dx2**2+dx1*(34.0d0*dx2+3.0d0*dx1))
      ci(6)=969.0d0*dx1**14*(12.0d0*dx2**2+dx1*(8.0d0*dx2+dx1))
      ci(7)=3876.0d0*dx1**13*(7.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))
      ci(8)=3876.0d0*dx1**12*(13.0d0*dx2**2+dx1*(14.0d0*dx2+3.0d0*dx1))
      ci(9)=1938.0d0*dx1**11*(39.0d0*dx2**2+2.0d0*dx1*(26.0d0*dx2+7.0d0*dx1))
      ci(10)=8398.0d0*dx1**10*(11.0d0*dx2**2+6.0d0*dx1*(3.0d0*dx2+dx1))
      ci(11)=8398.0d0*dx1**9*(11.0d0*dx2**2+dx1*(22.0d0*dx2+9.0d0*dx1))
      ci(12)=8398.0d0*dx1**8*(9.0d0*dx2**2+11.0d0*dx1*(2.0d0*dx2+dx1))
      ci(13)=8398.0d0*dx1**7*(6.0d0*dx2**2+dx1*(18.0d0*dx2+11.0d0*dx1))
      ci(14)=1938.0d0*dx1**6*(14.0d0*dx2**2+13.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))
      ci(15)=3876.0d0*dx1**5*(3.0d0*dx2**2+dx1*(14.0d0*dx2+13.0d0*dx1))
      ci(16)=3876.0d0*dx1**4*(dx2**2+dx1*(6.0d0*dx2+7.0d0*dx1))
      ci(17)=969.0d0*dx1**3*(dx2**2+4.0d0*dx1*(2.0d0*dx2+3.0d0*dx1))
      ci(18)=57.0d0*dx1**2*(3.0d0*dx2**2+34.0d0*dx1*(dx2+2.0d0*dx1))
      ci(19)=19.0d0*dx1*(dx2**2+3.0d0*dx1*(6.0d0*dx2+17.0d0*dx1))
      ci(20)=dx2**2+19.0d0*dx1*(2.0d0*dx2+9.0d0*dx1)
      ci(21)=2.0d0*dx2+19.0d0*dx1
      ci(22)=1.0d0
    case(3)
      ci(1)=dx1**19*dx2**3
      ci(2)=dx1**18*dx2**2*(19.0d0*dx2+3.0d0*dx1)
      ci(3)=3.0d0*dx1**17*dx2*(57.0d0*dx2**2+dx1*(19.0d0*dx2+dx1))
      ci(4)=dx1**16*(969.0d0*dx2**3+dx1*(513.0d0*dx2**2+dx1*(57.0d0*dx2+dx1)))
      ci(5)=19.0d0*dx1**15*(204.0d0*dx2**3+dx1*(153.0d0*dx2**2+dx1*(27.0d0*dx2+dx1)))
      ci(6)=171.0d0*dx1**14*(68.0d0*dx2**3+dx1*(68.0d0*dx2**2+dx1*(17.0d0*dx2+dx1)))
      ci(7)=969.0d0*dx1**13*(28.0d0*dx2**3+dx1*(36.0d0*dx2**2+dx1*(12.0d0*dx2+dx1)))
      ci(8)=3876.0d0*dx1**12*(13.0d0*dx2**3+dx1*(21.0d0*dx2**2+dx1*(9.0d0*dx2+dx1)))
      ci(9)=5814.0d0*dx1**11*(13.0d0*dx2**3+2.0d0*dx1*(13.0d0*dx2**2+dx1*(7.0d0*dx2+dx1)))
      ci(10)=646.0d0*dx1**10*(143.0d0*dx2**3+3.0d0*dx1*(117.0d0*dx2**2+2.0d0*dx1*(39.0d0*dx2+7.0d0*dx1)))
      ci(11)=8398.0d0*dx1**9*(11.0d0*dx2**3+3.0d0*dx1*(11.0d0*dx2**2+dx1*(9.0d0*dx2+2.0d0*dx1)))
      ci(12)=25194.0d0*dx1**8*(3.0d0*dx2**3+dx1*(11.0d0*dx2**2+dx1*(11.0d0*dx2+3.0d0*dx1)))
      ci(13)=8398.0d0*dx1**7*(6.0d0*dx2**3+dx1*(27.0d0*dx2**2+11.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(14)=646.0d0*dx1**6*(42.0d0*dx2**3+13.0d0*dx1*(18.0d0*dx2**2+dx1*(27.0d0*dx2+11.0d0*dx1)))
      ci(15)=5814.0d0*dx1**5*(2.0d0*dx2**3+dx1*(14.0d0*dx2**2+13.0d0*dx1*(2.0d0*dx2+dx1)))
      ci(16)=3876.0d0*dx1**4*(dx2**3+dx1*(9.0d0*dx2**2+dx1*(21.0d0*dx2+13.0d0*dx1)))
      ci(17)=969.0d0*dx1**3*(dx2**3+4.0d0*dx1*(3.0d0*dx2**2+dx1*(9.0d0*dx2+7.0d0*dx1)))
      ci(18)=171.0d0*dx1**2*(dx2**3+17.0d0*dx1*(dx2**2+4.0d0*dx1*(dx2+dx1)))
      ci(19)=19.0d0*dx1*(dx2**3+3.0d0*dx1*(9.0d0*dx2**2+17.0d0*dx1*(3.0d0*dx2+4.0d0*dx1)))
      ci(20)=dx2**3+57.0d0*dx1*(dx2**2+dx1*(9.0d0*dx2+17.0d0*dx1))
      ci(21)=3.0d0*(dx2**2+19.0d0*dx1*(dx2+3.0d0*dx1))
      ci(22)=3.0d0*dx2+19.0d0*dx1
      ci(23)=1.0d0
    case(4)
      ci(1)=dx1**19*dx2**4
      ci(2)=dx1**18*dx2**3*(19.0d0*dx2+4.0d0*dx1)
      ci(3)=dx1**17*dx2**2*(171.0d0*dx2**2+2.0d0*dx1*(38.0d0*dx2+3.0d0*dx1))
      ci(4)=dx1**16*dx2*(969.0d0*dx2**3+2.0d0*dx1*(342.0d0*dx2**2+dx1*(57.0d0*dx2+2.0d0*dx1)))
      ci(5)=dx1**15*(3876.0d0*dx2**4+dx1*(3876.0d0*dx2**3+dx1*(1026.0d0*dx2**2+dx1*(76.0d0*dx2+dx1))))
      ci(6)=19.0d0*dx1**14*(612.0d0*dx2**4+dx1*(816.0d0*dx2**3+dx1*(306.0d0*dx2**2+dx1*(36.0d0*dx2+dx1))))
      ci(7)=57.0d0*dx1**13*(476.0d0*dx2**4+dx1*(816.0d0*dx2**3+dx1*(408.0d0*dx2**2+dx1*(68.0d0*dx2+3.0d0*dx1))))
      ci(8)=969.0d0*dx1**12*(52.0d0*dx2**4+dx1*(112.0d0*dx2**3+dx1*(72.0d0*dx2**2+dx1*(16.0d0*dx2+dx1))))
      ci(9)=1938.0d0*dx1**11*(39.0d0*dx2**4+2.0d0*dx1*(52.0d0*dx2**3+dx1*(42.0d0*dx2**2+dx1*(12.0d0*dx2+dx1))))
      ci(10)=646.0d0*dx1**10*(143.0d0*dx2**4+6.0d0*dx1*(78.0d0*dx2**3+dx1*(78.0d0*dx2**2+dx1*(28.0d0*dx2+3.0d0*dx1))))
      ci(11)=646.0d0*dx1**9*(143.0d0*dx2**4+2.0d0*dx1*(286.0d0*dx2**3+3.0d0*dx1*(117.0d0*dx2**2+dx1*(52.0d0*dx2+7.0d0*dx1))))
      ci(12)=8398.0d0*dx1**8*(9.0d0*dx2**4+2.0d0*dx1*(22.0d0*dx2**3+3.0d0*dx1*(11.0d0*dx2**2+dx1*(6.0d0*dx2+dx1))))
      ci(13)=8398.0d0*dx1**7*(6.0d0*dx2**4+dx1*(36.0d0*dx2**3+dx1*(66.0d0*dx2**2+dx1*(44.0d0*dx2+9.0d0*dx1))))
      ci(14)=646.0d0*dx1**6*(42.0d0*dx2**4+13.0d0*dx1*(24.0d0*dx2**3+dx1*(54.0d0*dx2**2+11.0d0*dx1*(4.0d0*dx2+dx1))))
      ci(15)=646.0d0*dx1**5*(18.0d0*dx2**4+dx1*(168.0d0*dx2**3+13.0d0*dx1*(36.0d0*dx2**2+dx1*(36.0d0*dx2+11.0d0*dx1))))
      ci(16)=1938.0d0*dx1**4*(2.0d0*dx2**4+dx1*(24.0d0*dx2**3+dx1*(84.0d0*dx2**2+13.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))))
      ci(17)=969.0d0*dx1**3*(dx2**4+4.0d0*dx1*(4.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(28.0d0*dx2+13.0d0*dx1))))
      ci(18)=57.0d0*dx1**2*(3.0d0*dx2**4+68.0d0*dx1*(dx2**3+dx1*(6.0d0*dx2**2+dx1*(12.0d0*dx2+7.0d0*dx1))))
      ci(19)=19.0d0*dx1*(dx2**4+6.0d0*dx1*(6.0d0*dx2**3+17.0d0*dx1*(3.0d0*dx2**2+2.0d0*dx1*(4.0d0*dx2+3.0d0*dx1))))
      ci(20)=dx2**4+38.0d0*dx1*(2.0d0*dx2**3+3.0d0*dx1*(9.0d0*dx2**2+34.0d0*dx1*(dx2+dx1)))
      ci(21)=4.0d0*dx2**3+57.0d0*dx1*(2.0d0*dx2**2+dx1*(12.0d0*dx2+17.0d0*dx1))
      ci(22)=6.0d0*dx2**2+19.0d0*dx1*(4.0d0*dx2+9.0d0*dx1)
      ci(23)=4.0d0*dx2+19.0d0*dx1
      ci(24)=1.0d0
    case(5)
      ci(1)=dx1**19*dx2**5
      ci(2)=dx1**18*dx2**4*(19.0d0*dx2+5.0d0*dx1)
      ci(3)=dx1**17*dx2**3*(171.0d0*dx2**2+5.0d0*dx1*(19.0d0*dx2+2.0d0*dx1))
      ci(4)=dx1**16*dx2**2*(969.0d0*dx2**3+5.0d0*dx1*(171.0d0*dx2**2+2.0d0*dx1*(19.0d0*dx2+dx1)))
      ci(5)=dx1**15*dx2*(3876.0d0*dx2**4+5.0d0*dx1*(969.0d0*dx2**3+dx1*(342.0d0*dx2**2+dx1*(38.0d0*dx2+dx1))))
      ci(6)=dx1**14*(11628.0d0*dx2**5+dx1*(19380.0d0*dx2**4+dx1*(9690.0d0*dx2**3+dx1*(1710.0d0*dx2**2+dx1*(95.0d0*dx2+dx1)))))
      ci(7)=19.0d0*dx1**13*(1428.0d0*dx2**5+dx1*(3060.0d0*dx2**4+dx1*(2040.0d0*dx2**3+dx1*(510.0d0*dx2**2+dx1*(45.0d0*dx2+dx1)))))
      ci(8)=57.0d0*dx1**12*(884.0d0*dx2**5+dx1*(2380.0d0*dx2**4+dx1*(2040.0d0*dx2**3+dx1*(680.0d0*dx2**2+dx1*(85.0d0*dx2+3.0d0*dx1 &
)))))
      ci(9)=969.0d0*dx1**11*(78.0d0*dx2**5+dx1*(260.0d0*dx2**4+dx1*(280.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(20.0d0*dx2+dx1)))))
      ci(10)=646.0d0*dx1**10*(143.0d0*dx2**5+3.0d0*dx1*(195.0d0*dx2**4+2.0d0*dx1*(130.0d0*dx2**3+dx1*(70.0d0*dx2**2+dx1* &
(15.0d0*dx2+dx1)))))
      ci(11)=646.0d0*dx1**9*(143.0d0*dx2**5+dx1*(715.0d0*dx2**4+6.0d0*dx1*(195.0d0*dx2**3+dx1*(130.0d0*dx2**2+dx1*(35.0d0*dx2 &
+3.0d0*dx1)))))
      ci(12)=646.0d0*dx1**8*(117.0d0*dx2**5+dx1*(715.0d0*dx2**4+2.0d0*dx1*(715.0d0*dx2**3+3.0d0*dx1*(195.0d0*dx2**2+dx1* &
(65.0d0*dx2+7.0d0*dx1)))))
      ci(13)=8398.0d0*dx1**7*(6.0d0*dx2**5+dx1*(45.0d0*dx2**4+dx1*(110.0d0*dx2**3+dx1*(110.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2 &
+2.0d0*dx1)))))
      ci(14)=646.0d0*dx1**6*(42.0d0*dx2**5+13.0d0*dx1*(30.0d0*dx2**4+dx1*(90.0d0*dx2**3+dx1*(110.0d0*dx2**2+dx1*(55.0d0*dx2 &
+9.0d0*dx1)))))
      ci(15)=646.0d0*dx1**5*(18.0d0*dx2**5+dx1*(210.0d0*dx2**4+13.0d0*dx1*(60.0d0*dx2**3+dx1*(90.0d0*dx2**2+11.0d0*dx1*(5.0d0*dx2 &
+dx1)))))
      ci(16)=646.0d0*dx1**4*(6.0d0*dx2**5+dx1*(90.0d0*dx2**4+dx1*(420.0d0*dx2**3+13.0d0*dx1*(60.0d0*dx2**2+dx1*(45.0d0*dx2 &
+11.0d0*dx1)))))
      ci(17)=969.0d0*dx1**3*(dx2**5+2.0d0*dx1*(10.0d0*dx2**4+dx1*(60.0d0*dx2**3+dx1*(140.0d0*dx2**2+13.0d0*dx1*(10.0d0*dx2 &
+3.0d0*dx1)))))
      ci(18)=57.0d0*dx1**2*(3.0d0*dx2**5+17.0d0*dx1*(5.0d0*dx2**4+4.0d0*dx1*(10.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(35.0d0*dx2 &
+13.0d0*dx1)))))
      ci(19)=19.0d0*dx1*(dx2**5+3.0d0*dx1*(15.0d0*dx2**4+34.0d0*dx1*(5.0d0*dx2**3+2.0d0*dx1*(10.0d0*dx2**2+dx1*(15.0d0*dx2 &
+7.0d0*dx1)))))
      ci(20)=dx2**5+19.0d0*dx1*(5.0d0*dx2**4+6.0d0*dx1*(15.0d0*dx2**3+17.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))))
      ci(21)=5.0d0*dx2**4+19.0d0*dx1*(10.0d0*dx2**3+3.0d0*dx1*(30.0d0*dx2**2+17.0d0*dx1*(5.0d0*dx2+4.0d0*dx1)))
      ci(22)=10.0d0*dx2**3+19.0d0*dx1*(10.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+17.0d0*dx1))
      ci(23)=10.0d0*dx2**2+19.0d0*dx1*(5.0d0*dx2+9.0d0*dx1)
      ci(24)=5.0d0*dx2+19.0d0*dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>5, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^20 * (x-x2)^n2 as sum_k=0^(20+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_20(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**20
      ci(2)=20.0d0*dx1**19
      ci(3)=190.0d0*dx1**18
      ci(4)=1140.0d0*dx1**17
      ci(5)=4845.0d0*dx1**16
      ci(6)=15504.0d0*dx1**15
      ci(7)=38760.0d0*dx1**14
      ci(8)=77520.0d0*dx1**13
      ci(9)=125970.0d0*dx1**12
      ci(10)=167960.0d0*dx1**11
      ci(11)=184756.0d0*dx1**10
      ci(12)=167960.0d0*dx1**9
      ci(13)=125970.0d0*dx1**8
      ci(14)=77520.0d0*dx1**7
      ci(15)=38760.0d0*dx1**6
      ci(16)=15504.0d0*dx1**5
      ci(17)=4845.0d0*dx1**4
      ci(18)=1140.0d0*dx1**3
      ci(19)=190.0d0*dx1**2
      ci(20)=20.0d0*dx1
      ci(21)=1.0d0
    case(1)
      ci(1)=dx1**20*dx2
      ci(2)=dx1**19*(20.0d0*dx2+dx1)
      ci(3)=10.0d0*dx1**18*(19.0d0*dx2+2.0d0*dx1)
      ci(4)=190.0d0*dx1**17*(6.0d0*dx2+dx1)
      ci(5)=285.0d0*dx1**16*(17.0d0*dx2+4.0d0*dx1)
      ci(6)=969.0d0*dx1**15*(16.0d0*dx2+5.0d0*dx1)
      ci(7)=7752.0d0*dx1**14*(5.0d0*dx2+2.0d0*dx1)
      ci(8)=38760.0d0*dx1**13*(2.0d0*dx2+dx1)
      ci(9)=9690.0d0*dx1**12*(13.0d0*dx2+8.0d0*dx1)
      ci(10)=41990.0d0*dx1**11*(4.0d0*dx2+3.0d0*dx1)
      ci(11)=16796.0d0*dx1**10*(11.0d0*dx2+10.0d0*dx1)
      ci(12)=16796.0d0*dx1**9*(10.0d0*dx2+11.0d0*dx1)
      ci(13)=41990.0d0*dx1**8*(3.0d0*dx2+4.0d0*dx1)
      ci(14)=9690.0d0*dx1**7*(8.0d0*dx2+13.0d0*dx1)
      ci(15)=38760.0d0*dx1**6*(dx2+2.0d0*dx1)
      ci(16)=7752.0d0*dx1**5*(2.0d0*dx2+5.0d0*dx1)
      ci(17)=969.0d0*dx1**4*(5.0d0*dx2+16.0d0*dx1)
      ci(18)=285.0d0*dx1**3*(4.0d0*dx2+17.0d0*dx1)
      ci(19)=190.0d0*dx1**2*(dx2+6.0d0*dx1)
      ci(20)=10.0d0*dx1*(2.0d0*dx2+19.0d0*dx1)
      ci(21)=dx2+20.0d0*dx1
      ci(22)=1.0d0
    case(2)
      ci(1)=dx1**20*dx2**2
      ci(2)=2.0d0*dx1**19*dx2*(10.0d0*dx2+dx1)
      ci(3)=dx1**18*(190.0d0*dx2**2+dx1*(40.0d0*dx2+dx1))
      ci(4)=20.0d0*dx1**17*(57.0d0*dx2**2+dx1*(19.0d0*dx2+dx1))
      ci(5)=95.0d0*dx1**16*(51.0d0*dx2**2+2.0d0*dx1*(12.0d0*dx2+dx1))
      ci(6)=114.0d0*dx1**15*(136.0d0*dx2**2+5.0d0*dx1*(17.0d0*dx2+2.0d0*dx1))
      ci(7)=969.0d0*dx1**14*(40.0d0*dx2**2+dx1*(32.0d0*dx2+5.0d0*dx1))
      ci(8)=15504.0d0*dx1**13*(5.0d0*dx2**2+dx1*(5.0d0*dx2+dx1))
      ci(9)=9690.0d0*dx1**12*(13.0d0*dx2**2+4.0d0*dx1*(4.0d0*dx2+dx1))
      ci(10)=6460.0d0*dx1**11*(26.0d0*dx2**2+3.0d0*dx1*(13.0d0*dx2+4.0d0*dx1))
      ci(11)=8398.0d0*dx1**10*(22.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+3.0d0*dx1))
      ci(12)=33592.0d0*dx1**9*(5.0d0*dx2**2+dx1*(11.0d0*dx2+5.0d0*dx1))
      ci(13)=8398.0d0*dx1**8*(15.0d0*dx2**2+2.0d0*dx1*(20.0d0*dx2+11.0d0*dx1))
      ci(14)=6460.0d0*dx1**7*(12.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+2.0d0*dx1))
      ci(15)=9690.0d0*dx1**6*(4.0d0*dx2**2+dx1*(16.0d0*dx2+13.0d0*dx1))
      ci(16)=15504.0d0*dx1**5*(dx2**2+5.0d0*dx1*(dx2+dx1))
      ci(17)=969.0d0*dx1**4*(5.0d0*dx2**2+8.0d0*dx1*(4.0d0*dx2+5.0d0*dx1))
      ci(18)=114.0d0*dx1**3*(10.0d0*dx2**2+17.0d0*dx1*(5.0d0*dx2+8.0d0*dx1))
      ci(19)=95.0d0*dx1**2*(2.0d0*dx2**2+3.0d0*dx1*(8.0d0*dx2+17.0d0*dx1))
      ci(20)=20.0d0*dx1*(dx2**2+19.0d0*dx1*(dx2+3.0d0*dx1))
      ci(21)=dx2**2+10.0d0*dx1*(4.0d0*dx2+19.0d0*dx1)
      ci(22)=2.0d0*(dx2+10.0d0*dx1)
      ci(23)=1.0d0
    case(3)
      ci(1)=dx1**20*dx2**3
      ci(2)=dx1**19*dx2**2*(20.0d0*dx2+3.0d0*dx1)
      ci(3)=dx1**18*dx2*(190.0d0*dx2**2+3.0d0*dx1*(20.0d0*dx2+dx1))
      ci(4)=dx1**17*(1140.0d0*dx2**3+dx1*(570.0d0*dx2**2+dx1*(60.0d0*dx2+dx1)))
      ci(5)=5.0d0*dx1**16*(969.0d0*dx2**3+2.0d0*dx1*(342.0d0*dx2**2+dx1*(57.0d0*dx2+2.0d0*dx1)))
      ci(6)=19.0d0*dx1**15*(816.0d0*dx2**3+5.0d0*dx1*(153.0d0*dx2**2+2.0d0*dx1*(18.0d0*dx2+dx1)))
      ci(7)=57.0d0*dx1**14*(680.0d0*dx2**3+dx1*(816.0d0*dx2**2+5.0d0*dx1*(51.0d0*dx2+4.0d0*dx1)))
      ci(8)=969.0d0*dx1**13*(80.0d0*dx2**3+dx1*(120.0d0*dx2**2+dx1*(48.0d0*dx2+5.0d0*dx1)))
      ci(9)=1938.0d0*dx1**12*(65.0d0*dx2**3+4.0d0*dx1*(30.0d0*dx2**2+dx1*(15.0d0*dx2+2.0d0*dx1)))
      ci(10)=3230.0d0*dx1**11*(52.0d0*dx2**3+3.0d0*dx1*(39.0d0*dx2**2+4.0d0*dx1*(6.0d0*dx2+dx1)))
      ci(11)=646.0d0*dx1**10*(286.0d0*dx2**3+15.0d0*dx1*(52.0d0*dx2**2+dx1*(39.0d0*dx2+8.0d0*dx1)))
      ci(12)=8398.0d0*dx1**9*(20.0d0*dx2**3+3.0d0*dx1*(22.0d0*dx2**2+5.0d0*dx1*(4.0d0*dx2+dx1)))
      ci(13)=8398.0d0*dx1**8*(15.0d0*dx2**3+2.0d0*dx1*(30.0d0*dx2**2+dx1*(33.0d0*dx2+10.0d0*dx1)))
      ci(14)=646.0d0*dx1**7*(120.0d0*dx2**3+13.0d0*dx1*(45.0d0*dx2**2+2.0d0*dx1*(30.0d0*dx2+11.0d0*dx1)))
      ci(15)=3230.0d0*dx1**6*(12.0d0*dx2**3+dx1*(72.0d0*dx2**2+13.0d0*dx1*(9.0d0*dx2+4.0d0*dx1)))
      ci(16)=1938.0d0*dx1**5*(8.0d0*dx2**3+5.0d0*dx1*(12.0d0*dx2**2+dx1*(24.0d0*dx2+13.0d0*dx1)))
      ci(17)=969.0d0*dx1**4*(5.0d0*dx2**3+8.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(3.0d0*dx2+2.0d0*dx1)))
      ci(18)=57.0d0*dx1**3*(20.0d0*dx2**3+17.0d0*dx1*(15.0d0*dx2**2+8.0d0*dx1*(6.0d0*dx2+5.0d0*dx1)))
      ci(19)=19.0d0*dx1**2*(10.0d0*dx2**3+3.0d0*dx1*(60.0d0*dx2**2+17.0d0*dx1*(15.0d0*dx2+16.0d0*dx1)))
      ci(20)=5.0d0*dx1*(4.0d0*dx2**3+57.0d0*dx1*(2.0d0*dx2**2+dx1*(12.0d0*dx2+17.0d0*dx1)))
      ci(21)=dx2**3+30.0d0*dx1*(2.0d0*dx2**2+19.0d0*dx1*(dx2+2.0d0*dx1))
      ci(22)=3.0d0*dx2**2+10.0d0*dx1*(6.0d0*dx2+19.0d0*dx1)
      ci(23)=3.0d0*dx2+20.0d0*dx1
      ci(24)=1.0d0
    case(4)
      ci(1)=dx1**20*dx2**4
      ci(2)=4.0d0*dx1**19*dx2**3*(5.0d0*dx2+dx1)
      ci(3)=2.0d0*dx1**18*dx2**2*(95.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))
      ci(4)=4.0d0*dx1**17*dx2*(285.0d0*dx2**3+dx1*(190.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))
      ci(5)=dx1**16*(4845.0d0*dx2**4+dx1*(4560.0d0*dx2**3+dx1*(1140.0d0*dx2**2+dx1*(80.0d0*dx2+dx1))))
      ci(6)=4.0d0*dx1**15*(3876.0d0*dx2**4+5.0d0*dx1*(969.0d0*dx2**3+dx1*(342.0d0*dx2**2+dx1*(38.0d0*dx2+dx1))))
      ci(7)=38.0d0*dx1**14*(1020.0d0*dx2**4+dx1*(1632.0d0*dx2**3+5.0d0*dx1*(153.0d0*dx2**2+dx1*(24.0d0*dx2+dx1))))
      ci(8)=228.0d0*dx1**13*(340.0d0*dx2**4+dx1*(680.0d0*dx2**3+dx1*(408.0d0*dx2**2+5.0d0*dx1*(17.0d0*dx2+dx1))))
      ci(9)=969.0d0*dx1**12*(130.0d0*dx2**4+dx1*(320.0d0*dx2**3+dx1*(240.0d0*dx2**2+dx1*(64.0d0*dx2+5.0d0*dx1))))
      ci(10)=2584.0d0*dx1**11*(65.0d0*dx2**4+3.0d0*dx1*(65.0d0*dx2**3+2.0d0*dx1*(30.0d0*dx2**2+dx1*(10.0d0*dx2+dx1))))
      ci(11)=1292.0d0*dx1**10*(143.0d0*dx2**4+5.0d0*dx1*(104.0d0*dx2**3+3.0d0*dx1*(39.0d0*dx2**2+2.0d0*dx1*(8.0d0*dx2+dx1))))
      ci(12)=2584.0d0*dx1**9*(65.0d0*dx2**4+dx1*(286.0d0*dx2**3+15.0d0*dx1*(26.0d0*dx2**2+dx1*(13.0d0*dx2+2.0d0*dx1))))
      ci(13)=8398.0d0*dx1**8*(15.0d0*dx2**4+dx1*(80.0d0*dx2**3+dx1*(132.0d0*dx2**2+5.0d0*dx1*(16.0d0*dx2+3.0d0*dx1))))
      ci(14)=2584.0d0*dx1**7*(30.0d0*dx2**4+13.0d0*dx1*(15.0d0*dx2**3+dx1*(30.0d0*dx2**2+dx1*(22.0d0*dx2+5.0d0*dx1))))
      ci(15)=1292.0d0*dx1**6*(30.0d0*dx2**4+dx1*(240.0d0*dx2**3+13.0d0*dx1*(45.0d0*dx2**2+dx1*(40.0d0*dx2+11.0d0*dx1))))
      ci(16)=2584.0d0*dx1**5*(6.0d0*dx2**4+5.0d0*dx1*(12.0d0*dx2**3+dx1*(36.0d0*dx2**2+13.0d0*dx1*(3.0d0*dx2+dx1))))
      ci(17)=969.0d0*dx1**4*(5.0d0*dx2**4+2.0d0*dx1*(32.0d0*dx2**3+5.0d0*dx1*(24.0d0*dx2**2+dx1*(32.0d0*dx2+13.0d0*dx1))))
      ci(18)=228.0d0*dx1**3*(5.0d0*dx2**4+17.0d0*dx1*(5.0d0*dx2**3+4.0d0*dx1*(6.0d0*dx2**2+5.0d0*dx1*(2.0d0*dx2+dx1))))
      ci(19)=38.0d0*dx1**2*(5.0d0*dx2**4+3.0d0*dx1*(40.0d0*dx2**3+17.0d0*dx1*(15.0d0*dx2**2+4.0d0*dx1*(8.0d0*dx2+5.0d0*dx1))))
      ci(20)=4.0d0*dx1*(5.0d0*dx2**4+19.0d0*dx1*(10.0d0*dx2**3+3.0d0*dx1*(30.0d0*dx2**2+17.0d0*dx1*(5.0d0*dx2+4.0d0*dx1))))
      ci(21)=dx2**4+5.0d0*dx1*(16.0d0*dx2**3+57.0d0*dx1*(4.0d0*dx2**2+dx1*(16.0d0*dx2+17.0d0*dx1)))
      ci(22)=4.0d0*(dx2**3+5.0d0*dx1*(6.0d0*dx2**2+19.0d0*dx1*(2.0d0*dx2+3.0d0*dx1)))
      ci(23)=2.0d0*(3.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+19.0d0*dx1))
      ci(24)=4.0d0*(dx2+5.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>4, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^21 * (x-x2)^n2 as sum_k=0^(21+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_21(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**21
      ci(2)=21.0d0*dx1**20
      ci(3)=210.0d0*dx1**19
      ci(4)=1330.0d0*dx1**18
      ci(5)=5985.0d0*dx1**17
      ci(6)=20349.0d0*dx1**16
      ci(7)=54264.0d0*dx1**15
      ci(8)=116280.0d0*dx1**14
      ci(9)=203490.0d0*dx1**13
      ci(10)=293930.0d0*dx1**12
      ci(11)=352716.0d0*dx1**11
      ci(12)=352716.0d0*dx1**10
      ci(13)=293930.0d0*dx1**9
      ci(14)=203490.0d0*dx1**8
      ci(15)=116280.0d0*dx1**7
      ci(16)=54264.0d0*dx1**6
      ci(17)=20349.0d0*dx1**5
      ci(18)=5985.0d0*dx1**4
      ci(19)=1330.0d0*dx1**3
      ci(20)=210.0d0*dx1**2
      ci(21)=21.0d0*dx1
      ci(22)=1.0d0
    case(1)
      ci(1)=dx1**21*dx2
      ci(2)=dx1**20*(21.0d0*dx2+dx1)
      ci(3)=21.0d0*dx1**19*(10.0d0*dx2+dx1)
      ci(4)=70.0d0*dx1**18*(19.0d0*dx2+3.0d0*dx1)
      ci(5)=665.0d0*dx1**17*(9.0d0*dx2+2.0d0*dx1)
      ci(6)=1197.0d0*dx1**16*(17.0d0*dx2+5.0d0*dx1)
      ci(7)=6783.0d0*dx1**15*(8.0d0*dx2+3.0d0*dx1)
      ci(8)=7752.0d0*dx1**14*(15.0d0*dx2+7.0d0*dx1)
      ci(9)=29070.0d0*dx1**13*(7.0d0*dx2+4.0d0*dx1)
      ci(10)=22610.0d0*dx1**12*(13.0d0*dx2+9.0d0*dx1)
      ci(11)=58786.0d0*dx1**11*(6.0d0*dx2+5.0d0*dx1)
      ci(12)=352716.0d0*dx1**10*(dx2+dx1)
      ci(13)=58786.0d0*dx1**9*(5.0d0*dx2+6.0d0*dx1)
      ci(14)=22610.0d0*dx1**8*(9.0d0*dx2+13.0d0*dx1)
      ci(15)=29070.0d0*dx1**7*(4.0d0*dx2+7.0d0*dx1)
      ci(16)=7752.0d0*dx1**6*(7.0d0*dx2+15.0d0*dx1)
      ci(17)=6783.0d0*dx1**5*(3.0d0*dx2+8.0d0*dx1)
      ci(18)=1197.0d0*dx1**4*(5.0d0*dx2+17.0d0*dx1)
      ci(19)=665.0d0*dx1**3*(2.0d0*dx2+9.0d0*dx1)
      ci(20)=70.0d0*dx1**2*(3.0d0*dx2+19.0d0*dx1)
      ci(21)=21.0d0*dx1*(dx2+10.0d0*dx1)
      ci(22)=dx2+21.0d0*dx1
      ci(23)=1.0d0
    case(2)
      ci(1)=dx1**21*dx2**2
      ci(2)=dx1**20*dx2*(21.0d0*dx2+2.0d0*dx1)
      ci(3)=dx1**19*(210.0d0*dx2**2+dx1*(42.0d0*dx2+dx1))
      ci(4)=7.0d0*dx1**18*(190.0d0*dx2**2+3.0d0*dx1*(20.0d0*dx2+dx1))
      ci(5)=35.0d0*dx1**17*(171.0d0*dx2**2+2.0d0*dx1*(38.0d0*dx2+3.0d0*dx1))
      ci(6)=133.0d0*dx1**16*(153.0d0*dx2**2+10.0d0*dx1*(9.0d0*dx2+dx1))
      ci(7)=399.0d0*dx1**15*(136.0d0*dx2**2+3.0d0*dx1*(34.0d0*dx2+5.0d0*dx1))
      ci(8)=969.0d0*dx1**14*(120.0d0*dx2**2+7.0d0*dx1*(16.0d0*dx2+3.0d0*dx1))
      ci(9)=1938.0d0*dx1**13*(105.0d0*dx2**2+4.0d0*dx1*(30.0d0*dx2+7.0d0*dx1))
      ci(10)=3230.0d0*dx1**12*(91.0d0*dx2**2+18.0d0*dx1*(7.0d0*dx2+2.0d0*dx1))
      ci(11)=4522.0d0*dx1**11*(78.0d0*dx2**2+5.0d0*dx1*(26.0d0*dx2+9.0d0*dx1))
      ci(12)=58786.0d0*dx1**10*(6.0d0*dx2**2+dx1*(12.0d0*dx2+5.0d0*dx1))
      ci(13)=58786.0d0*dx1**9*(5.0d0*dx2**2+6.0d0*dx1*(2.0d0*dx2+dx1))
      ci(14)=4522.0d0*dx1**8*(45.0d0*dx2**2+26.0d0*dx1*(5.0d0*dx2+3.0d0*dx1))
      ci(15)=3230.0d0*dx1**7*(36.0d0*dx2**2+7.0d0*dx1*(18.0d0*dx2+13.0d0*dx1))
      ci(16)=1938.0d0*dx1**6*(28.0d0*dx2**2+15.0d0*dx1*(8.0d0*dx2+7.0d0*dx1))
      ci(17)=969.0d0*dx1**5*(21.0d0*dx2**2+8.0d0*dx1*(14.0d0*dx2+15.0d0*dx1))
      ci(18)=399.0d0*dx1**4*(15.0d0*dx2**2+34.0d0*dx1*(3.0d0*dx2+4.0d0*dx1))
      ci(19)=133.0d0*dx1**3*(10.0d0*dx2**2+9.0d0*dx1*(10.0d0*dx2+17.0d0*dx1))
      ci(20)=35.0d0*dx1**2*(6.0d0*dx2**2+19.0d0*dx1*(4.0d0*dx2+9.0d0*dx1))
      ci(21)=7.0d0*dx1*(3.0d0*dx2**2+10.0d0*dx1*(6.0d0*dx2+19.0d0*dx1))
      ci(22)=dx2**2+42.0d0*dx1*(dx2+5.0d0*dx1)
      ci(23)=2.0d0*dx2+21.0d0*dx1
      ci(24)=1.0d0
    case(3)
      ci(1)=dx1**21*dx2**3
      ci(2)=3.0d0*dx1**20*dx2**2*(7.0d0*dx2+dx1)
      ci(3)=3.0d0*dx1**19*dx2*(70.0d0*dx2**2+dx1*(21.0d0*dx2+dx1))
      ci(4)=dx1**18*(1330.0d0*dx2**3+dx1*(630.0d0*dx2**2+dx1*(63.0d0*dx2+dx1)))
      ci(5)=21.0d0*dx1**17*(285.0d0*dx2**3+dx1*(190.0d0*dx2**2+dx1*(30.0d0*dx2+dx1)))
      ci(6)=21.0d0*dx1**16*(969.0d0*dx2**3+5.0d0*dx1*(171.0d0*dx2**2+2.0d0*dx1*(19.0d0*dx2+dx1)))
      ci(7)=133.0d0*dx1**15*(408.0d0*dx2**3+dx1*(459.0d0*dx2**2+5.0d0*dx1*(27.0d0*dx2+2.0d0*dx1)))
      ci(8)=171.0d0*dx1**14*(680.0d0*dx2**3+7.0d0*dx1*(136.0d0*dx2**2+dx1*(51.0d0*dx2+5.0d0*dx1)))
      ci(9)=2907.0d0*dx1**13*(70.0d0*dx2**3+dx1*(120.0d0*dx2**2+7.0d0*dx1*(8.0d0*dx2+dx1)))
      ci(10)=646.0d0*dx1**12*(455.0d0*dx2**3+3.0d0*dx1*(315.0d0*dx2**2+4.0d0*dx1*(45.0d0*dx2+7.0d0*dx1)))
      ci(11)=1938.0d0*dx1**11*(182.0d0*dx2**3+5.0d0*dx1*(91.0d0*dx2**2+3.0d0*dx1*(21.0d0*dx2+4.0d0*dx1)))
      ci(12)=13566.0d0*dx1**10*(26.0d0*dx2**3+dx1*(78.0d0*dx2**2+5.0d0*dx1*(13.0d0*dx2+3.0d0*dx1)))
      ci(13)=58786.0d0*dx1**9*(5.0d0*dx2**3+dx1*(18.0d0*dx2**2+dx1*(18.0d0*dx2+5.0d0*dx1)))
      ci(14)=13566.0d0*dx1**8*(15.0d0*dx2**3+13.0d0*dx1*(5.0d0*dx2**2+2.0d0*dx1*(3.0d0*dx2+dx1)))
      ci(15)=1938.0d0*dx1**7*(60.0d0*dx2**3+7.0d0*dx1*(45.0d0*dx2**2+13.0d0*dx1*(5.0d0*dx2+2.0d0*dx1)))
      ci(16)=646.0d0*dx1**6*(84.0d0*dx2**3+5.0d0*dx1*(108.0d0*dx2**2+7.0d0*dx1*(27.0d0*dx2+13.0d0*dx1)))
      ci(17)=2907.0d0*dx1**5*(7.0d0*dx2**3+2.0d0*dx1*(28.0d0*dx2**2+5.0d0*dx1*(12.0d0*dx2+7.0d0*dx1)))
      ci(18)=171.0d0*dx1**4*(35.0d0*dx2**3+17.0d0*dx1*(21.0d0*dx2**2+8.0d0*dx1*(7.0d0*dx2+5.0d0*dx1)))
      ci(19)=133.0d0*dx1**3*(10.0d0*dx2**3+3.0d0*dx1*(45.0d0*dx2**2+17.0d0*dx1*(9.0d0*dx2+8.0d0*dx1)))
      ci(20)=21.0d0*dx1**2*(10.0d0*dx2**3+19.0d0*dx1*(10.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+17.0d0*dx1)))
      ci(21)=21.0d0*dx1*(dx2**3+5.0d0*dx1*(6.0d0*dx2**2+19.0d0*dx1*(2.0d0*dx2+3.0d0*dx1)))
      ci(22)=dx2**3+7.0d0*dx1*(9.0d0*dx2**2+10.0d0*dx1*(9.0d0*dx2+19.0d0*dx1))
      ci(23)=3.0d0*(dx2**2+7.0d0*dx1*(3.0d0*dx2+10.0d0*dx1))
      ci(24)=3.0d0*(dx2+7.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>3, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^22 * (x-x2)^n2 as sum_k=0^(22+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_22(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**22
      ci(2)=22.0d0*dx1**21
      ci(3)=231.0d0*dx1**20
      ci(4)=1540.0d0*dx1**19
      ci(5)=7315.0d0*dx1**18
      ci(6)=26334.0d0*dx1**17
      ci(7)=74613.0d0*dx1**16
      ci(8)=170544.0d0*dx1**15
      ci(9)=319770.0d0*dx1**14
      ci(10)=497420.0d0*dx1**13
      ci(11)=646646.0d0*dx1**12
      ci(12)=705432.0d0*dx1**11
      ci(13)=646646.0d0*dx1**10
      ci(14)=497420.0d0*dx1**9
      ci(15)=319770.0d0*dx1**8
      ci(16)=170544.0d0*dx1**7
      ci(17)=74613.0d0*dx1**6
      ci(18)=26334.0d0*dx1**5
      ci(19)=7315.0d0*dx1**4
      ci(20)=1540.0d0*dx1**3
      ci(21)=231.0d0*dx1**2
      ci(22)=22.0d0*dx1
      ci(23)=1.0d0
    case(1)
      ci(1)=dx1**22*dx2
      ci(2)=dx1**21*(22.0d0*dx2+dx1)
      ci(3)=11.0d0*dx1**20*(21.0d0*dx2+2.0d0*dx1)
      ci(4)=77.0d0*dx1**19*(20.0d0*dx2+3.0d0*dx1)
      ci(5)=385.0d0*dx1**18*(19.0d0*dx2+4.0d0*dx1)
      ci(6)=1463.0d0*dx1**17*(18.0d0*dx2+5.0d0*dx1)
      ci(7)=4389.0d0*dx1**16*(17.0d0*dx2+6.0d0*dx1)
      ci(8)=10659.0d0*dx1**15*(16.0d0*dx2+7.0d0*dx1)
      ci(9)=21318.0d0*dx1**14*(15.0d0*dx2+8.0d0*dx1)
      ci(10)=35530.0d0*dx1**13*(14.0d0*dx2+9.0d0*dx1)
      ci(11)=49742.0d0*dx1**12*(13.0d0*dx2+10.0d0*dx1)
      ci(12)=58786.0d0*dx1**11*(12.0d0*dx2+11.0d0*dx1)
      ci(13)=58786.0d0*dx1**10*(11.0d0*dx2+12.0d0*dx1)
      ci(14)=49742.0d0*dx1**9*(10.0d0*dx2+13.0d0*dx1)
      ci(15)=35530.0d0*dx1**8*(9.0d0*dx2+14.0d0*dx1)
      ci(16)=21318.0d0*dx1**7*(8.0d0*dx2+15.0d0*dx1)
      ci(17)=10659.0d0*dx1**6*(7.0d0*dx2+16.0d0*dx1)
      ci(18)=4389.0d0*dx1**5*(6.0d0*dx2+17.0d0*dx1)
      ci(19)=1463.0d0*dx1**4*(5.0d0*dx2+18.0d0*dx1)
      ci(20)=385.0d0*dx1**3*(4.0d0*dx2+19.0d0*dx1)
      ci(21)=77.0d0*dx1**2*(3.0d0*dx2+20.0d0*dx1)
      ci(22)=11.0d0*dx1*(2.0d0*dx2+21.0d0*dx1)
      ci(23)=dx2+22.0d0*dx1
      ci(24)=1.0d0
    case(2)
      ci(1)=dx1**22*dx2**2
      ci(2)=2.0d0*dx1**21*dx2*(11.0d0*dx2+dx1)
      ci(3)=dx1**20*(231.0d0*dx2**2+dx1*(44.0d0*dx2+dx1))
      ci(4)=22.0d0*dx1**19*(70.0d0*dx2**2+dx1*(21.0d0*dx2+dx1))
      ci(5)=77.0d0*dx1**18*(95.0d0*dx2**2+dx1*(40.0d0*dx2+3.0d0*dx1))
      ci(6)=154.0d0*dx1**17*(171.0d0*dx2**2+5.0d0*dx1*(19.0d0*dx2+2.0d0*dx1))
      ci(7)=1463.0d0*dx1**16*(51.0d0*dx2**2+dx1*(36.0d0*dx2+5.0d0*dx1))
      ci(8)=1254.0d0*dx1**15*(136.0d0*dx2**2+7.0d0*dx1*(17.0d0*dx2+3.0d0*dx1))
      ci(9)=10659.0d0*dx1**14*(30.0d0*dx2**2+dx1*(32.0d0*dx2+7.0d0*dx1))
      ci(10)=14212.0d0*dx1**13*(35.0d0*dx2**2+3.0d0*dx1*(15.0d0*dx2+4.0d0*dx1))
      ci(11)=7106.0d0*dx1**12*(91.0d0*dx2**2+5.0d0*dx1*(28.0d0*dx2+9.0d0*dx1))
      ci(12)=9044.0d0*dx1**11*(78.0d0*dx2**2+11.0d0*dx1*(13.0d0*dx2+5.0d0*dx1))
      ci(13)=58786.0d0*dx1**10*(11.0d0*dx2**2+dx1*(24.0d0*dx2+11.0d0*dx1))
      ci(14)=9044.0d0*dx1**9*(55.0d0*dx2**2+13.0d0*dx1*(11.0d0*dx2+6.0d0*dx1))
      ci(15)=7106.0d0*dx1**8*(45.0d0*dx2**2+7.0d0*dx1*(20.0d0*dx2+13.0d0*dx1))
      ci(16)=14212.0d0*dx1**7*(12.0d0*dx2**2+5.0d0*dx1*(9.0d0*dx2+7.0d0*dx1))
      ci(17)=10659.0d0*dx1**6*(7.0d0*dx2**2+2.0d0*dx1*(16.0d0*dx2+15.0d0*dx1))
      ci(18)=1254.0d0*dx1**5*(21.0d0*dx2**2+17.0d0*dx1*(7.0d0*dx2+8.0d0*dx1))
      ci(19)=1463.0d0*dx1**4*(5.0d0*dx2**2+3.0d0*dx1*(12.0d0*dx2+17.0d0*dx1))
      ci(20)=154.0d0*dx1**3*(10.0d0*dx2**2+19.0d0*dx1*(5.0d0*dx2+9.0d0*dx1))
      ci(21)=77.0d0*dx1**2*(3.0d0*dx2**2+5.0d0*dx1*(8.0d0*dx2+19.0d0*dx1))
      ci(22)=22.0d0*dx1*(dx2**2+7.0d0*dx1*(3.0d0*dx2+10.0d0*dx1))
      ci(23)=dx2**2+11.0d0*dx1*(4.0d0*dx2+21.0d0*dx1)
      ci(24)=2.0d0*(dx2+11.0d0*dx1)
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>2, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^23 * (x-x2)^n2 as sum_k=0^(23+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_23(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**23
      ci(2)=23.0d0*dx1**22
      ci(3)=253.0d0*dx1**21
      ci(4)=1771.0d0*dx1**20
      ci(5)=8855.0d0*dx1**19
      ci(6)=33649.0d0*dx1**18
      ci(7)=100947.0d0*dx1**17
      ci(8)=245157.0d0*dx1**16
      ci(9)=490314.0d0*dx1**15
      ci(10)=817190.0d0*dx1**14
      ci(11)=1144066.0d0*dx1**13
      ci(12)=1352078.0d0*dx1**12
      ci(13)=1352078.0d0*dx1**11
      ci(14)=1144066.0d0*dx1**10
      ci(15)=817190.0d0*dx1**9
      ci(16)=490314.0d0*dx1**8
      ci(17)=245157.0d0*dx1**7
      ci(18)=100947.0d0*dx1**6
      ci(19)=33649.0d0*dx1**5
      ci(20)=8855.0d0*dx1**4
      ci(21)=1771.0d0*dx1**3
      ci(22)=253.0d0*dx1**2
      ci(23)=23.0d0*dx1
      ci(24)=1.0d0
    case(1)
      ci(1)=dx1**23*dx2
      ci(2)=dx1**22*(23.0d0*dx2+dx1)
      ci(3)=23.0d0*dx1**21*(11.0d0*dx2+dx1)
      ci(4)=253.0d0*dx1**20*(7.0d0*dx2+dx1)
      ci(5)=1771.0d0*dx1**19*(5.0d0*dx2+dx1)
      ci(6)=1771.0d0*dx1**18*(19.0d0*dx2+5.0d0*dx1)
      ci(7)=33649.0d0*dx1**17*(3.0d0*dx2+dx1)
      ci(8)=14421.0d0*dx1**16*(17.0d0*dx2+7.0d0*dx1)
      ci(9)=245157.0d0*dx1**15*(2.0d0*dx2+dx1)
      ci(10)=163438.0d0*dx1**14*(5.0d0*dx2+3.0d0*dx1)
      ci(11)=163438.0d0*dx1**13*(7.0d0*dx2+5.0d0*dx1)
      ci(12)=104006.0d0*dx1**12*(13.0d0*dx2+11.0d0*dx1)
      ci(13)=1352078.0d0*dx1**11*(dx2+dx1)
      ci(14)=104006.0d0*dx1**10*(11.0d0*dx2+13.0d0*dx1)
      ci(15)=163438.0d0*dx1**9*(5.0d0*dx2+7.0d0*dx1)
      ci(16)=163438.0d0*dx1**8*(3.0d0*dx2+5.0d0*dx1)
      ci(17)=245157.0d0*dx1**7*(dx2+2.0d0*dx1)
      ci(18)=14421.0d0*dx1**6*(7.0d0*dx2+17.0d0*dx1)
      ci(19)=33649.0d0*dx1**5*(dx2+3.0d0*dx1)
      ci(20)=1771.0d0*dx1**4*(5.0d0*dx2+19.0d0*dx1)
      ci(21)=1771.0d0*dx1**3*(dx2+5.0d0*dx1)
      ci(22)=253.0d0*dx1**2*(dx2+7.0d0*dx1)
      ci(23)=23.0d0*dx1*(dx2+11.0d0*dx1)
      ci(24)=dx2+23.0d0*dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>1, here n2=',n2
      stop
  end select
 
end subroutine

!!> expand the product (x-x1)^24 * (x-x2)^n2 as sum_k=0^(24+n2) ci(k+1)*(x-x3)^k
!!  
recursive subroutine expand_centered_product_24(x1,x2,n2,x3,ci)
 
  implicit none
 
  ! input variables
  integer     , intent(in):: n2
  real(kind=8), intent(in):: x1
  real(kind=8), intent(in):: x2
  real(kind=8), intent(in):: x3
  real(kind=8), intent(inout), dimension(*):: ci
 
  ! local variables
  real(kind=8):: dx1
  real(kind=8):: dx2
 
  ! compute displacements
  dx1=x3-x1
  dx2=x3-x2
 
  select case(n2)
    case(0)
      ci(1)=dx1**24
      ci(2)=24.0d0*dx1**23
      ci(3)=276.0d0*dx1**22
      ci(4)=2024.0d0*dx1**21
      ci(5)=10626.0d0*dx1**20
      ci(6)=42504.0d0*dx1**19
      ci(7)=134596.0d0*dx1**18
      ci(8)=346104.0d0*dx1**17
      ci(9)=735471.0d0*dx1**16
      ci(10)=1307504.0d0*dx1**15
      ci(11)=1961256.0d0*dx1**14
      ci(12)=2496144.0d0*dx1**13
      ci(13)=2704156.0d0*dx1**12
      ci(14)=2496144.0d0*dx1**11
      ci(15)=1961256.0d0*dx1**10
      ci(16)=1307504.0d0*dx1**9
      ci(17)=735471.0d0*dx1**8
      ci(18)=346104.0d0*dx1**7
      ci(19)=134596.0d0*dx1**6
      ci(20)=42504.0d0*dx1**5
      ci(21)=10626.0d0*dx1**4
      ci(22)=2024.0d0*dx1**3
      ci(23)=276.0d0*dx1**2
      ci(24)=24.0d0*dx1
      ci(25)=1.0d0
    case default
      print*,'Error: expand_centered_product for n2>0, here n2=',n2
      stop
  end select
 
end subroutine


end module
