
!> @file R2Moment.f90
!!
!! Defines utility routines for computing overlap
!! integrals between gaussian basis elements, weighted by r^2 moment powers .
!!
!! Author: I. Duchemin July 2015
!!


!> Two centers weighted overlap integral for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) |r-R_2|^{n_{R_2}}
!!  \f$
!!  
!!  nx1+ny1+nz1+nx2+ny2+nz2+nr2 must stay <= 14
!!  
recursive function C_R2Moment_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,nr2)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: nr2      !< |r-R_2|^nr2 power for second cubic Harmonic rdial weight (nr2 should be either 0,2,4 or 6 at most)
  
  ! return value
  real(kind=8)  :: C_R2Moment_C
  
  ! local variable
  integer       :: l4
  integer       :: l5
  integer       :: ir1
  integer       :: ir2
  integer       :: ir3
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: c4(680)
  real(kind=8)  :: c5(680)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: alpha3
  real(kind=8)  :: alpha4
  real(kind=8)  :: alpha5
  
  ! check that exponent is even
  if ( mod(nr2,2).ne.0 ) then
    print *,"Error: inconsistent exponent nr2, detected in - C_R2Moment_C -"
    print *,"       nr2 should be an even number"
    print *,"       actual value is :",nr2
    stop -1
  end if
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 1st member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2=0.0d0
  c2(ir2)=1.0d0
  
  ! forge |r-R_2| coeff in the R basis
  c3=0.0d0
  select case (nr2)
    case (0)
      ir3 =ir_index(0,0,0)
      c3(ir3)=1.0d0
    case (2)
      ir3 =ir_index(2,0,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,2,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,0,2)
      c3(ir3)=1.0d0
    case (4)
      ir3 =ir_index(4,0,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,4,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,0,4)
      c3(ir3)=1.0d0
      ir3 =ir_index(2,2,0)
      c3(ir3)=2.0d0
      ir3 =ir_index(2,0,2)
      c3(ir3)=2.0d0
      ir3 =ir_index(0,2,2)
      c3(ir3)=2.0d0
    case (6)
      ir3 =ir_index(6,0,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,6,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,0,6)
      c3(ir3)=1.0d0
      ir3 =ir_index(4,2,0)
      c3(ir3)=3.0d0
      ir3 =ir_index(4,0,2)
      c3(ir3)=3.0d0
      ir3 =ir_index(2,4,0)
      c3(ir3)=3.0d0
      ir3 =ir_index(0,4,2)
      c3(ir3)=3.0d0
      ir3 =ir_index(2,0,4)
      c3(ir3)=3.0d0
      ir3 =ir_index(0,2,4)
      c3(ir3)=3.0d0
      ir3 =ir_index(2,2,2)
      c3(ir3)=6.0d0
    case default
      print *,"Error: invalid exponent nr2, detected in - C_R2Moment_C -"
      print *,"       nr2 should be <= 6"
      print *,"       actual value is :",nr2
      stop -1
  end select
  
  ! get product of |r-R_2|^nr2 times the second cubic harmonic
  r3    =r2
  alpha3=0.0d0
  call RxR_to_R(alpha2,r2,c2,nx2+ny2+nz2,alpha3,r3,c3,nr2,alpha4,r4,c4,l4)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha4,r4,c4,l4,alpha5,r5,c5,l5)
  
  ! compute 1 norm of decomposition
  C_R2Moment_C =R_1_norm(alpha5,c5,l5)
  
end function
  
  
!> Two centers weighted overlap integral for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) |r-R_2|^{n_{R_2}}
!!  \f$
!!  
!!  nx1+ny1+nz1+nx2+ny2+nz2+nr2 must stay <= 14
!!  
recursive function Y_R2Moment_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2,nr2)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: nr2      !< |r-R_2|^nr2 power for second cubic Harmonic rdial weight (nr2 should be either 0,2,4 or 6 at most)
  
  ! return value
  real(kind=8)  :: Y_R2Moment_Y
  
  ! local variable
  integer       :: l4
  integer       :: l5
  integer       :: ir1
  integer       :: ir2
  integer       :: ir3
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: c4(680)
  real(kind=8)  :: c5(680)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: alpha3
  real(kind=8)  :: alpha4
  real(kind=8)  :: alpha5
  
  ! check that exponent is even
  if ( mod(nr2,2).ne.0 ) then
    print *,"Error: inconsistent exponent nr2, detected in - C_R2Moment_C -"
    print *,"       nr2 should be an even number"
    print *,"       actual value is :",nr2
    stop -1
  end if
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
    ! forge |r-R_2| coeff in the R basis
  c3=0.0d0
  select case (nr2)
    case (0)
      ir3 =ir_index(0,0,0)
      c3(ir3)=1.0d0
    case (2)
      ir3 =ir_index(2,0,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,2,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,0,2)
      c3(ir3)=1.0d0
    case (4)
      ir3 =ir_index(4,0,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,4,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,0,4)
      c3(ir3)=1.0d0
      ir3 =ir_index(2,2,0)
      c3(ir3)=2.0d0
      ir3 =ir_index(2,0,2)
      c3(ir3)=2.0d0
      ir3 =ir_index(0,2,2)
      c3(ir3)=2.0d0
    case (6)
      ir3 =ir_index(6,0,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,6,0)
      c3(ir3)=1.0d0
      ir3 =ir_index(0,0,6)
      c3(ir3)=1.0d0
      ir3 =ir_index(4,2,0)
      c3(ir3)=3.0d0
      ir3 =ir_index(4,0,2)
      c3(ir3)=3.0d0
      ir3 =ir_index(2,4,0)
      c3(ir3)=3.0d0
      ir3 =ir_index(0,4,2)
      c3(ir3)=3.0d0
      ir3 =ir_index(2,0,4)
      c3(ir3)=3.0d0
      ir3 =ir_index(0,2,4)
      c3(ir3)=3.0d0
      ir3 =ir_index(2,2,2)
      c3(ir3)=6.0d0
    case default
      print *,"Error: invalid exponent nr2, detected in - C_R2Moment_C -"
      print *,"       nr2 should be <= 6"
      print *,"       actual value is :",nr2
      stop -1
  end select
  
  ! get product of |r-R_2|^nr2 times the second cubic harmonic
  r3    =r2
  alpha3=0.0d0
  call RxR_to_R(alpha2,r2,c2,l2,alpha3,r3,c3,nr2,alpha4,r4,c4,l4)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha4,r4,c4,l4,alpha5,r5,c5,l5)
  
  ! compute 1 norm of decomposition
  Y_R2Moment_Y =R_1_norm(alpha5,c5,l5)
  
end function

