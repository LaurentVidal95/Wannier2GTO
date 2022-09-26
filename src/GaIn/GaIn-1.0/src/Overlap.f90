
!> @file Overlap.f90
!!
!! Defines utility routines for computing overlap
!! integrals between gaussian basis elements.
!!
!! Author: I. Duchemin July 2015
!!


!> Two centers overlap integral for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2)
!!  \f$
!!  
!!  nx+ny+nz must stay <= 6 for each spherical harmonics
!!  
recursive function C_Overlap_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2)
  
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
  
  ! return value
  real(kind=8)  :: C_Overlap_C
  
  ! local variable
  integer       :: l3
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: a3
  
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
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,a3,r3,c3,l3)
  
  ! compute 1 norm of decomposition
  C_Overlap_C =R_1_norm(a3,c3,l3)
  
end function
  
!> Two centers overlap integral for Cubic/Solid Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2)
!!  \f$
!!  
!!  nx+ny+nz must stay <= 6 for each spherical harmonics
!!  
recursive function C_Overlap_Y(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2)
  
  use mod_R_1_norm
  use mod_R_from_Y
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
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: C_Overlap_Y
  
  ! local variable
  integer       :: l3
  integer       :: ir1
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: a3
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! set 1st member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,l2,a3,r3,c3,l3)
  
  ! compute 1 norm of decomposition
  C_Overlap_Y =R_1_norm(a3,c3,l3)
  
end function


!> Two centers overlap integral for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2)
!!  \f$
!!  
!!  nx+ny+nz must stay <= 6 for each spherical harmonics
!!  
recursive function Y_Overlap_C(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: Y_Overlap_C
  
  ! local variable
  integer       :: l3
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: a3
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! set 2nd member coeff in the R basis
  c2=0.0d0
  c2(ir2)=1.0d0
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,nx2+ny2+nz2,a3,r3,c3,l3)
  
  ! compute 1 norm of decomposition
  Y_Overlap_C =R_1_norm(a3,c3,l3)
  
end function


!> Two centers overlap integral for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2)
!!  \f$
!!  
!!  l must stay <= 6 for each spherical harmonics
!!  
recursive function Y_Overlap_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2)
  
  use mod_R_1_norm
  use mod_R_from_Y
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
  
  ! return value
  real(kind=8)  :: Y_Overlap_Y
  
  ! local variable
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: a3
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,a3,r3,c3,l3)
  
  ! compute 1 norm of decomposition
  Y_Overlap_Y =R_1_norm(a3,c3,l3)
  
end function

!> Three centers overlap integral for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) Y_{xyz}^{(3)}(r-R_3)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function CC_Overlap_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,nx3,ny3,nz3)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  real(kind=8) :: alpha3   !< exponent for third cubic Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  
  ! return value
  real(kind=8)  :: CC_Overlap_C
  
  ! local variable
  integer       :: l4
  integer       :: l5
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: a4
  real(kind=8)  :: a5
  
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
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,a4,r4,c3,l4)
  
  ! get 3rd member index in the R basis
  ir1 =ir_index(nx3,ny3,nz3)
  
  ! set 3rd member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! get product of the three cubic harmonics
  call RxR_to_R(a4,r4,c3,l4,alpha3,r3,c1,nx3+ny3+nz3,a5,r5,c2,l5)
  
  ! compute 1 norm of decomposition
  CC_Overlap_C =R_1_norm(a5,c2,l5)
  
end function


!> Three centers overlap integral for Cubic/Solid Spherical Harmonics (C/Y)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) Y_{lm}^{(3)}(r-R_3)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function CC_Overlap_Y(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,l3,m3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  real(kind=8) :: alpha3   !< exponent for third cubic Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  
  ! return value
  real(kind=8)  :: CC_Overlap_Y
  
  ! local variable
  integer       :: l4
  integer       :: l5
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: a4
  real(kind=8)  :: a5
  
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
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,a4,r4,c3,l4)
  
  ! decomposition of the third solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c1)
  
  ! get product of the three cubic harmonics
  call RxR_to_R(a4,r4,c3,l4,alpha3,r3,c1,l3,a5,r5,c2,l5)
  
  ! compute 1 norm of decomposition
  CC_Overlap_Y =R_1_norm(a5,c2,l5)
  
end function


!> Three centers overlap integral for Cubic/Solid Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) Y_{xyz}^{(3)}(r-R_3)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function CY_Overlap_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2,alpha3,r3,nx3,ny3,nz3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  real(kind=8) :: alpha3   !< exponent for third cubic Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  
  ! return value
  real(kind=8)  :: CY_Overlap_C
  
  ! local variable
  integer       :: l4
  integer       :: l5
  integer       :: ir1
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: a4
  real(kind=8)  :: a5
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! set 1st member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,l2,a4,r4,c3,l4)
  
  ! get 3rd member index in the R basis
  ir1 =ir_index(nx3,ny3,nz3)
  
  ! set 3rd member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! get product of the three cubic harmonics
  call RxR_to_R(a4,r4,c3,l4,alpha3,r3,c1,nx3+ny3+nz3,a5,r5,c2,l5)
  
  ! compute 1 norm of decomposition
  CY_Overlap_C =R_1_norm(a5,c2,l5)
  
end function


!> Three centers overlap integral for Cubic/Solid Spherical Harmonics (C/Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) Y_{xyz}^{(3)}(r-R_3)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function YC_Overlap_C(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,nx3,ny3,nz3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  real(kind=8) :: alpha3   !< exponent for third cubic Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  
  ! return value
  real(kind=8)  :: YC_Overlap_C
  
  ! local variable
  integer       :: l4
  integer       :: l5
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: a4
  real(kind=8)  :: a5
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! set 2nd member coeff in the R basis
  c2=0.0d0
  c2(ir2)=1.0d0
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,nx2+ny2+nz2,a4,r4,c3,l4)
  
  ! get 3rd member index in the R basis
  ir1 =ir_index(nx3,ny3,nz3)
  
  ! set 3rd member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! get product of the three cubic harmonics
  call RxR_to_R(a4,r4,c3,l4,alpha3,r3,c1,nx3+ny3+nz3,a5,r5,c2,l5)
  
  ! compute 1 norm of decomposition
  YC_Overlap_C =R_1_norm(a5,c2,l5)
  
end function

!> Three centers overlap integral for Cubic/Solid Spherical Harmonics (C/Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) Y_{xyz}^{(3)}(r-R_3)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function YY_Overlap_C(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,nx3,ny3,nz3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  
  ! return value
  real(kind=8)  :: YY_Overlap_C
  
  ! local variable
  integer       :: ir1
  integer       :: l4
  integer       :: l5
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: a4
  real(kind=8)  :: a5
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,a4,r4,c3,l4)
  
  ! get 3rd member index in the R basis
  ir1 =ir_index(nx3,ny3,nz3)
  
  ! set 3rd member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha3,r3,c1,nx3+ny3+nz3,a4,r4,c3,l4,a5,r5,c2,l5)
  
  ! compute 1 norm of decomposition
  YY_Overlap_C =R_1_norm(a5,c2,l5)
  
end function

!> Three centers overlap integral for Cubic/Solid Spherical Harmonics (C/Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) Y_{lm}^{(3)}(r-R_3)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function YC_Overlap_Y(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,l3,m3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: YC_Overlap_Y
  
  ! local variable
  integer       :: ir2
  integer       :: l4
  integer       :: l5
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: a4
  real(kind=8)  :: a5
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! set 2nd member coeff in the R basis
  c2=0.0d0
  c2(ir2)=1.0d0
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,nx2+ny2+nz2,a4,r4,c3,l4)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c1)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha3,r3,c1,l3,a4,r4,c3,l4,a5,r5,c2,l5)
  
  ! compute 1 norm of decomposition
  YC_Overlap_Y =R_1_norm(a5,c2,l5)
  
end function


!> Three centers overlap integral for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) Y_{lm}^{(3)}(r-R_3)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function CY_Overlap_Y(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2,alpha3,r3,l3,m3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  
  ! return value
  real(kind=8)  :: CY_Overlap_Y
  
  ! local variable
  integer       :: ir1
  integer       :: l4
  integer       :: l5
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: a4
  real(kind=8)  :: a5
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! set 1st member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,l2,a4,r4,c3,l4)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c1)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha3,r3,c1,l3,a4,r4,c3,l4,a5,r5,c2,l5)
  
  ! compute 1 norm of decomposition
  CY_Overlap_Y =R_1_norm(a5,c2,l5)
  
end function

!> Three centers overlap integral for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) Y_{lm}^{(3)}(r-R_3)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function YY_Overlap_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,l3,m3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  
  ! return value
  real(kind=8)  :: YY_Overlap_Y
  
  ! local variable
  integer       :: l4
  integer       :: l5
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r4(3)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: a4
  real(kind=8)  :: a5
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,a4,r4,c3,l4)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c1)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha3,r3,c1,l3,a4,r4,c3,l4,a5,r5,c2,l5)
  
  ! compute 1 norm of decomposition
  YY_Overlap_Y =R_1_norm(a5,c2,l5)
  
end function




!> Four centers overlap integral for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) Y_{lm}^{(3)}(r-R_3) Y_{lm}^{(4)}(r-R_4)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!
recursive function YY_Overlap_YY(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,l3,m3,alpha4,r4,l4,m4)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: r4(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  real(kind=8) :: alpha4   !< exponent for third spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: l4       !< angular momentum for third spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  integer      :: m4       !< orbital momentum for third spherical Harmonic
  
  ! return value
  real(kind=8)  :: YY_Overlap_YY
  
  ! local variable
  integer       :: l5
  integer       :: l6
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: r6(3)
  real(kind=8)  :: a5
  real(kind=8)  :: a6
  
    ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,a5,r5,c3,l5)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c1)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha3,r3,c1,l3,a5,r5,c3,l5,a6,r6,c2,l6)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l4,m4,c1)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha4,r4,c1,l4,a6,r6,c2,l6,a5,r5,c3,l5)
  
  ! compute 1 norm of decomposition
  YY_Overlap_YY =R_1_norm(a5,c3,l5)
  
end function



!> Four centers overlap integral for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) Y_{xyz}^{(3)}(r-R_3) Y_{xyz}^{(4)}(r-R_4)
!!  \f$
!!  
!!  l_tot must stay <= 14
!!  
recursive function CC_Overlap_CC(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,nx3,ny3,nz3,alpha4,r4,nx4,ny4,nz4)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: r4(3)    !< center for fourth cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  real(kind=8) :: alpha3   !< exponent for third cubic Harmonic
  real(kind=8) :: alpha4   !< exponent for fourth cubic Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  integer      :: nx4      !< x power for fourth cubic Harmonic
  integer      :: ny4      !< y power for fourth cubic Harmonic
  integer      :: nz4      !< z power for fourth cubic Harmonic
  
  ! return value
  real(kind=8)  :: CC_Overlap_CC
  
  ! local variable
  integer       :: l5
  integer       :: l6
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(680)
  real(kind=8)  :: c2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r5(3)
  real(kind=8)  :: r6(3)
  real(kind=8)  :: a5
  real(kind=8)  :: a6

  
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
  
    ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,a5,r5,c3,l5)
  
  ! get 3rd member index in the R basis
  ir1 =ir_index(nx3,ny3,nz3)
  
  ! set 3rd member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha3,r3,c1,nx3+ny3+nz3,a5,r5,c3,l5,a6,r6,c2,l6)
  
  ! get 3rd member index in the R basis
  ir1 =ir_index(nx4,ny4,nz4)
  
  ! set 3rd member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha4,r4,c1,nx4+ny4+nz4,a6,r6,c2,l6,a5,r5,c3,l5)
  
  ! compute 1 norm of decomposition
  CC_Overlap_CC =R_1_norm(a5,c3,l5)
  
end function



!> Provide an estimate of the upper bound for two center overlap integral as
!!
!!  \f$
!!      \int dr \, |r-R_1|^l_1 e^{-a1|r-R_1|^2}  |r-R_2|^l_1 e^{-a2|r-R_2|^2}
!!  \f$
!!  
!!  when l1 ot l2 is odd, the bound is estimated through either l1+1 or l2+1
!!  l1 and l2 must stay < 6.
!!  
recursive function Overlap_upper_bound(alpha1,r1,l1,alpha2,r2,l2) result(value)
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: value
  
  ! local parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  
  ! local variables
  integer      :: lmin
  integer      :: lmax
  real(kind=8) :: a1
  real(kind=8) :: a2
  real(kind=8) :: d2
  
  ! sort min and max l channel
  if ( l1.le.l2 ) then
    lmin =l1
    lmax =l2
    a1   =alpha1
    a2   =alpha2
  else
    lmin =l2
    lmax =l1
    a1   =alpha2
    a2   =alpha1
  end if
  
  ! approx to even channel
  if ( MOD(lmin, 2).ne.0 ) then
    lmin=lmin+1
  end if
  if ( MOD(lmax, 2).ne.0 ) then
    lmax=lmax+1
  end if
  
  ! compute d^2
  d2 =sum((r1-r2)**2)
  
  ! compute value
  select case(lmax)
    case(0)
      select case(lmin)
        case(0)
          value=(PI**1*sqrt(PI)*exp(-(a1*a2*d2)/(a1+a2)))/((a1+a2)**1*sqrt((a1+a2)))
        case default
          print *,'Error: detected in-Overlap_upper_bound-'
          print *,'       passed lmin=',lmin
          print *,'       passed lmax=',lmax
          print *,'       lmin and lmax should be >=0 and <=6'
          stop 1
      end select
    case(2)
      select case(lmin)
        case(0)
          value=(PI**1*sqrt(PI)*(a1*(2.0d0*a1*d2+3.0d0)+3.0d0*a2)*exp(-(a1*a2*d2)/(a1+a2)))/(2.0d0*(a1+a2)**3*sqrt((a1+a2)))
        case(2)
          value=(PI**1*sqrt(PI)*(4.0d0*a1**2*a2**2*d2**2+(6.0d0*a2**3-2.0d0*a1*a2**2-2.0d0*a1**2*a2+6.0d0*a1**3)*d2+15.0d0*a2**2 &
+30.0d0*a1*a2+15.0d0*a1**2)*exp(-(a1*a2*d2)/(a1+a2)))/(4.0d0*(a1+a2)**5*sqrt((a1+a2)))
        case default
          print *,'Error: detected in-Overlap_upper_bound-'
          print *,'       passed lmin=',lmin
          print *,'       passed lmax=',lmax
          print *,'       lmin and lmax should be >=0 and <=6'
          stop 1
      end select
    case(4)
      select case(lmin)
        case(0)
          value=(PI**1*sqrt(PI)*(4.0d0*a1**4*d2**2+(20.0d0*a1**2*a2+20.0d0*a1**3)*d2+15.0d0*a2**2+30.0d0*a1*a2+15.0d0*a1**2)*exp(- &
(a1*a2*d2)/(a1+a2)))/(4.0d0*(a1+a2)**5*sqrt((a1+a2)))
        case(2)
          value=(PI**1*sqrt(PI)*(8.0d0*a1**4*a2**2*d2**3+(40.0d0*a1**2*a2**3+8.0d0*a1**3*a2**2-20.0d0*a1**4*a2+12.0d0*a1**5)*d2**2 &
+(30.0d0*a2**4-20.0d0*a1*a2**3-30.0d0*a1**2*a2**2+120.0d0*a1**3*a2+100.0d0*a1**4)*d2+105.0d0*a2**3+315.0d0*a1*a2**2 &
+315.0d0*a1**2*a2+105.0d0*a1**3)*exp(-(a1*a2*d2)/(a1+a2)))/(8.0d0*(a1+a2)**7*sqrt((a1+a2)))
        case(4)
          value=(PI**1*sqrt(PI)*(16.0d0*a1**4*a2**4*d2**4+(80.0d0*a1**2*a2**5-48.0d0*a1**3*a2**4-48.0d0*a1**4*a2**3 &
+80.0d0*a1**5*a2**2)*d2**3+(60.0d0*a2**6-200.0d0*a1*a2**5+172.0d0*a1**2*a2**4+864.0d0*a1**3*a2**3+172.0d0*a1**4*a2**2 &
-200.0d0*a1**5*a2+60.0d0*a1**6)*d2**2+(700.0d0*a2**5+980.0d0*a1*a2**4-560.0d0*a1**2*a2**3-560.0d0*a1**3*a2**2+980.0d0*a1**4*a2 &
+700.0d0*a1**5)*d2+945.0d0*a2**4+3780.0d0*a1*a2**3+5670.0d0*a1**2*a2**2+3780.0d0*a1**3*a2+945.0d0*a1**4)*exp(-(a1*a2*d2)/(a1+a2)) &
)/(16.0d0*(a1+a2)**9*sqrt((a1+a2)))
        case default
          print *,'Error: detected in-Overlap_upper_bound-'
          print *,'       passed lmin=',lmin
          print *,'       passed lmax=',lmax
          print *,'       lmin and lmax should be >=0 and <=6'
          stop 1
      end select
    case(6)
      select case(lmin)
        case(0)
          value=(PI**1*sqrt(PI)*(8.0d0*a1**6*d2**3+(84.0d0*a1**4*a2+84.0d0*a1**5)*d2**2+(210.0d0*a1**2*a2**2+420.0d0*a1**3*a2 &
+210.0d0*a1**4)*d2+105.0d0*a2**3+315.0d0*a1*a2**2+315.0d0*a1**2*a2+105.0d0*a1**3)*exp(-(a1*a2*d2)/(a1+a2)))/(8.0d0*(a1+a2)**7*sqrt &
((a1+a2)))
        case(2)
          value=(PI**1*sqrt(PI)*(16.0d0*a1**6*a2**2*d2**4+(168.0d0*a1**4*a2**3+72.0d0*a1**5*a2**2-72.0d0*a1**6*a2+24.0d0*a1**7 &
)*d2**3+(420.0d0*a1**2*a2**4+168.0d0*a1**3*a2**3-504.0d0*a1**4*a2**2+168.0d0*a1**5*a2+420.0d0*a1**6)*d2**2+(210.0d0*a2**5 &
-210.0d0*a1*a2**4-420.0d0*a1**2*a2**3+2100.0d0*a1**3*a2**2+3570.0d0*a1**4*a2+1470.0d0*a1**5)*d2+945.0d0*a2**4+3780.0d0*a1*a2**3 &
+5670.0d0*a1**2*a2**2+3780.0d0*a1**3*a2+945.0d0*a1**4)*exp(-(a1*a2*d2)/(a1+a2)))/(16.0d0*(a1+a2)**9*sqrt((a1+a2)))
        case(4)
          value=(PI**1*sqrt(PI)*(32.0d0*a1**6*a2**4*d2**5+(336.0d0*a1**4*a2**5-48.0d0*a1**5*a2**4-224.0d0*a1**6*a2**3 &
+160.0d0*a1**7*a2**2)*d2**4+(840.0d0*a1**2*a2**6-1008.0d0*a1**3*a2**5-1224.0d0*a1**4*a2**4+2976.0d0*a1**5*a2**3 &
+1512.0d0*a1**6*a2**2-720.0d0*a1**7*a2+120.0d0*a1**8)*d2**3+(420.0d0*a2**7-2100.0d0*a1*a2**6+2772.0d0*a1**2*a2**5 &
+15708.0d0*a1**3*a2**4+6132.0d0*a1**4*a2**3-7812.0d0*a1**5*a2**2-588.0d0*a1**6*a2+2940.0d0*a1**7)*d2**2+(6300.0d0*a2**6 &
+10080.0d0*a1*a2**5-9450.0d0*a1**2*a2**4-12600.0d0*a1**3*a2**3+25200.0d0*a1**4*a2**2+37800.0d0*a1**5*a2+13230.0d0*a1**6)*d2 &
+10395.0d0*a2**5+51975.0d0*a1*a2**4+103950.0d0*a1**2*a2**3+103950.0d0*a1**3*a2**2+51975.0d0*a1**4*a2+10395.0d0*a1**5)*exp(- &
(a1*a2*d2)/(a1+a2)))/(32.0d0*(a1+a2)**11*sqrt((a1+a2)))
        case(6)
          value=(PI**1*sqrt(PI)*(64.0d0*a1**6*a2**6*d2**6+(672.0d0*a1**4*a2**7-480.0d0*a1**5*a2**6-480.0d0*a1**6*a2**5 &
+672.0d0*a1**7*a2**4)*d2**5+(1680.0d0*a1**2*a2**8-4704.0d0*a1**3*a2**7+384.0d0*a1**4*a2**6+13536.0d0*a1**5*a2**5 &
+384.0d0*a1**6*a2**4-4704.0d0*a1**7*a2**3+1680.0d0*a1**8*a2**2)*d2**4+(840.0d0*a2**9-7560.0d0*a1*a2**8+24192.0d0*a1**2*a2**7 &
+46080.0d0*a1**3*a2**6-43200.0d0*a1**4*a2**5-43200.0d0*a1**5*a2**4+46080.0d0*a1**6*a2**3+24192.0d0*a1**7*a2**2-7560.0d0*a1**8*a2 &
+840.0d0*a1**9)*d2**3+(26460.0d0*a2**8-21168.0d0*a1*a2**7-115668.0d0*a1**2*a2**6+151200.0d0*a1**3*a2**5+438480.0d0*a1**4*a2**4 &
+151200.0d0*a1**5*a2**3-115668.0d0*a1**6*a2**2-21168.0d0*a1**7*a2+26460.0d0*a1**8)*d2**2+(145530.0d0*a2**7+478170.0d0*a1*a2**6 &
+353430.0d0*a1**2*a2**5-311850.0d0*a1**3*a2**4-311850.0d0*a1**4*a2**3+353430.0d0*a1**5*a2**2+478170.0d0*a1**6*a2+145530.0d0*a1**7 &
)*d2+135135.0d0*a2**6+810810.0d0*a1*a2**5+2027025.0d0*a1**2*a2**4+2702700.0d0*a1**3*a2**3+2027025.0d0*a1**4*a2**2 &
+810810.0d0*a1**5*a2+135135.0d0*a1**6)*exp(-(a1*a2*d2)/(a1+a2)))/(64.0d0*(a1+a2)**13*sqrt((a1+a2)))
        case default
          print *,'Error: detected in-Overlap_upper_bound-'
          print *,'       passed lmin=',lmin
          print *,'       passed lmax=',lmax
          print *,'       lmin and lmax should be >=0 and <=6'
          stop 1
      end select
    case default
      print *,'Error: detected in-Overlap_upper_bound-'
      print *,'       passed lmin=',lmin
      print *,'       passed lmax=',lmax
      print *,'       lmin and lmax should be >=0 and <=6'
      stop 1
  end select
  
end function
