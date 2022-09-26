
!> @file Laplacian.f90
!!
!! Defines utility routines for computing laplacian
!! integrals between gaussian basis elements.
!!
!! Author: I. Duchemin July 2015
!!

!> Two centers Laplacian integrals for Cubic Spherical Harmonics (C)
!!  
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) \left(\nabla^2 \, Y_{xyz}^{(2)}(r-R_2)\right)
!!  \f$
!!  
!!  nx+ny+nz must stay <= 6 for each spherical harmonics
!!  
recursive function C_Laplacian_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2)
  
  use mod_R_1_norm
  use mod_R_Laplacian
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: C_Laplacian_C
  
  ! local variable
  integer       :: l3
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: t2(455)
  real(kind=8)  :: c3(680)
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
  
  ! apply laplacian to the right side
  call R_Laplacian(alpha2,c2,nx2+ny2+nz2,t2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,t2,nx2+ny2+nz2+2,a3,r3,c3,l3)
  
  ! compute 1 norm of decomposition
  C_Laplacian_C =R_1_norm(a3,c3,l3)
  
end function

!> Two centers Laplacian integrals for Cubic/Solid Spherical Harmonics (C/Y)
!!  
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) \left(\nabla^2 \, Y_{lm}^{(2)}(r-R_2)\right)
!!  \f$
!!  
!!  nx+ny+nz must stay <= 6 for each spherical harmonics
!!  
recursive function C_Laplacian_Y(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_R_Laplacian
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: C_Laplacian_Y
  
  ! local variable
  integer       :: l3
  integer       :: ir1
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: t2(455)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: a3
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! set 1st member coeff in the R basis
  c1=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! apply laplacian to the right side
  call R_Laplacian(alpha2,c2,l2,t2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,t2,l2+2,a3,r3,c3,l3)
  
  ! compute 1 norm of decomposition
  C_Laplacian_Y =R_1_norm(a3,c3,l3)
  
end function


!> Two centers Laplacian integrals for Cubic/Solid Spherical Harmonics (C/Y)
!!  
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) \left(\nabla^2 \, Y_{xyz}^{(2)}(r-R_2)\right)
!!  \f$
!!  
!!  nx+ny+nz must stay <= 6 for each spherical harmonics
!!  
recursive function Y_Laplacian_C(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_R_Laplacian
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: Y_Laplacian_C
  
  ! local variable
  integer       :: l3
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: t2(455)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: a3
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 2nd member coeff in the R basis
  c2=0.0d0
  c2(ir2)=1.0d0
  
  ! apply laplacian to the right side
  call R_Laplacian(alpha2,c2,nx2+ny2+nz2,t2)
  
  ! get product of the two cubic harmonics
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,t2,nx2+ny2+nz2+2,a3,r3,c3,l3)
  
  ! compute 1 norm of decomposition
  Y_Laplacian_C =R_1_norm(a3,c3,l3)
  
end function

!> Two centers Laplacian integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) \left(\nabla^2 \, Y_{lm}^{(2)}(r-R_2)\right)
!!  \f$
!!  
!!  l must stay <= 6 for each spherical harmonics
!!  
recursive function Y_Laplacian_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2)
  
  use mod_R_1_norm
  use mod_R_Laplacian
  use mod_R_from_Y
  use mod_CubicHarmonicsProduct
  use mod_rlYlm_laplacian
  
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
  real(kind=8)  :: Y_Laplacian_Y
  
  ! local variable
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: t2(680)
  real(kind=8)  :: c3(680)
  real(kind=8)  :: r3(3)
  real(kind=8)  :: a3
  
  ! general routine for high order momenta
  if ( l1.gt.4 .or. l2.gt.4 ) then
    
    ! decomposition of the first solid harmonic into cubic harmonic
    call R_from_Y(l1,m1,c1)
    
    ! decomposition of the second solid harmonic into cubic harmonic
    call R_from_Y(l2,m2,c2)
    
    ! apply laplacian to the right side
    call R_Laplacian(alpha2,c2,l2,t2)
    
    ! get product of the two cubic harmonics
    call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,t2,l2+2,a3,r3,c3,l3)
    
    ! compute 1 norm of decomposition
    Y_Laplacian_Y = R_1_norm(a3,c3,l3)
    
  ! specialized (faster) routine for low momentum
  else
    
    Y_Laplacian_Y = rlYlm_laplacian(alpha1,r1,l1,m1,alpha2,r2,l2,m2)
    
  end if
  
end function
