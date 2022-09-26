
!> This module defines utility routines to expand RxR product on R basis
!!
!! These routines are generated automatically through maxima scripting.
!!
!! Author: I. Duchemin July 2015
!!
module mod_CubicHarmonicsProduct

  implicit none

  integer, dimension(680) :: nx_from_ir =(/0,1,0,0,2,1,1,0,0,0,3,2,2,1,1,1,0,0,0,0,4,3,3,2,2,2,1,1,1,1,0,0,0,0,0,5,4,4,3,3,3,2,2,2 &
,2,1,1,1,1,1,0,0,0,0,0,0,6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1 &
,1,1,1,0,0,0,0,0,0,0,0,8,7,7,6,6,6,5,5,5,5,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,9,8,8,7,7,7,6,6,6 &
,6,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,10,9,9,8,8,8,7,7,7,7,6,6,6,6,6,5,5,5 &
,5,5,5,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,11,10,10,9,9,9,8,8,8,8,7,7,7,7,7 &
,6,6,6,6,6,6,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,12 &
,11,11,10,10,10,9,9,9,9,8,8,8,8,8,7,7,7,7,7,7,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2 &
,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,13,12,12,11,11,11,10,10,10,10,9,9,9,9,9,8,8,8,8,8,8,7,7,7,7,7,7,7,6,6,6,6 &
,6,6,6,6,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0 &
,0,0,0,0,0,0,0,0,14,13,13,12,12,12,11,11,11,11,10,10,10,10,10,9,9,9,9,9,9,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,5,5,5,5 &
,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0 &
,0,0,0,0,0,0/)
  integer, dimension(680) :: ny_from_ir =(/0,0,1,0,0,1,0,2,1,0,0,1,0,2,1,0,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,0,1,0,2,1,0,3,2,1 &
,0,4,3,2,1,0,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3 &
,2,1,0,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1 &
,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2 &
,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,10,9,8,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4 &
,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,10,9,8,7,6,5,4,3,2,1,0,11,10,9,8,7,6,5,4,3,2,1,0,0,1 &
,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,10,9,8,7,6,5,4,3,2,1,0 &
,11,10,9,8,7,6,5,4,3,2,1,0,12,11,10,9,8,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8 &
,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,10,9,8,7,6,5,4,3,2,1,0,11,10,9,8,7,6,5,4,3,2,1,0,12,11,10,9,8,7,6,5,4,3,2,1,0,13,12,11,10,9,8 &
,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,10 &
,9,8,7,6,5,4,3,2,1,0,11,10,9,8,7,6,5,4,3,2,1,0,12,11,10,9,8,7,6,5,4,3,2,1,0,13,12,11,10,9,8,7,6,5,4,3,2,1,0,14,13,12,11,10,9,8,7,6 &
,5,4,3,2,1,0/)
  integer, dimension(680) :: nz_from_ir =(/0,0,0,1,0,0,1,0,1,2,0,0,1,0,1,2,0,1,2,3,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,0,1,0,1,2,0,1,2 &
,3,0,1,2,3,4,0,1,2,3,4,5,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3 &
,4,5,6,0,1,2,3,4,5,6,7,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,0,1,0,1,2,0,1,2 &
,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3 &
,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,10,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1 &
,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9,10,11,0,0 &
,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,10 &
,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,12,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0 &
,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,12,0,1,2,3,4,5,6,7 &
,8,9,10,11,12,13,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,1 &
,2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,12,0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3,4,5,6,7,8,9,10 &
,11,12,13,14/)
  integer, dimension(0:14) :: irmax_from_lmax =(/1,4,10,20,35,56,84,120,165,220,286,364,455,560,680/)

contains

!> expansion of (x-x1)^i1*(y-y1)^j1*(z-z1)^k1*exp(-a1*|r-r1|^2)*(x-x2)^i2*(y-y2)^j2*(z-z2)^k2*exp(-a2*|r-r2|^2)
!! as: prefac * ( sum_k coeffs(k) * (x-xc)^ik*(y-yc)^jk*(z-zc)^kk*exp(-alpha*|r-rc|^2) )
!!
recursive subroutine RxR_to_R(a1,r1,c1,l1max,a2,r2,c2,l2max,ac,rc,cc,l3max)
  
  use mod_CenteredPowerExpansion
  
  implicit none
  
  ! input parameters
  integer     , intent(in)  :: l1max     !< maximum l reached for the first orbital components
  integer     , intent(in)  :: l2max     !< maximum l reached for the second orbital components
  integer     , intent(out) :: l3max     !< maximum l reached for the combined orbital components
  real(kind=8), intent(in)  :: c1(*)     !< coeff of the first orbital in the r basis
  real(kind=8), intent(in)  :: c2(*)     !< coeff of the second orbital in the r basis
  real(kind=8), intent(in)  :: a1        !< gaussian exponent for first orbital
  real(kind=8), intent(in)  :: a2        !< gaussian exponent for second orbital
  real(kind=8), intent(in)  :: r1(3)     !< first orbital center
  real(kind=8), intent(in)  :: r2(3)     !< second orbital center
  real(kind=8), intent(out) :: ac        !< combined exponant          : a1+a2
  real(kind=8), intent(out) :: rc(3)     !< weighted center of mass    : (a1*r1+a2*r5)/(a1+a2) 
  real(kind=8), intent(out) :: cc(*)     !< coefficients of decomposition (untouched beyond irmax_from_lmax(l3max)
  
  ! local variable
  integer      :: l
  integer      :: i
  integer      :: j
  integer      :: k
  integer      :: p
  integer      :: ir1
  integer      :: ir2
  integer      :: nx1
  integer      :: nx2
  integer      :: ny1
  integer      :: ny2
  integer      :: nz1
  integer      :: nz2
  real(kind=8) :: d12
  real(kind=8) :: cx(15)       !< coefficient of expansion in terms of (x-xc)
  real(kind=8) :: cy(15)       !< coefficient of expansion in terms of (y-yc)
  real(kind=8) :: cz(15)       !< coefficient of expansion in terms of (z-zc)
  real(kind=8) :: prefac       !< prefactor of decomposition : exp(a1*a2*d12/(a1+a2))
  
  real(kind=8) :: c1_c2
  real(kind=8) :: c1_c2_cx
  
  real(kind=8) :: cc_ref(680) 
  
  
  ! check maxl
  if ( l1max+l2max>14 ) then
    print *,'Error: RxR_to_R not implemented for l1+l2>14, here l1+l2=',l1max+l2max
    stop
  end if
  
  ! compute squared distance between centers
  d12=sum((r1-r2)**2) 
  
  ! compute combined exponant
  ac=a1+a2
  
  ! compute center of mass
  if ( d12>0.0d0 ) then
    rc=(a1*r1+a2*r2)/ac
  else
    rc=r1
  end if
  
  ! compute prefactor of decomposition
  prefac=exp(-a1*a2/ac*d12)
  
  ! compute coefficients of expansion
  cc(1:455)=0.0d0
  cc(1:irmax_from_lmax(l1max+l2max))=0.0d0
  do ir1=1,irmax_from_lmax(l1max)
    if ( c1(ir1).ne.0.0d0 ) then 
      do ir2=1,irmax_from_lmax(l2max)
        if ( c2(ir2).ne.0.0d0 ) then 
          ! get corresponding x, y and z powers 
          nx1 =nx_from_ir(ir1)
          nx2 =nx_from_ir(ir2)
          ny1 =ny_from_ir(ir1)
          ny2 =ny_from_ir(ir2)
          nz1 =nz_from_ir(ir1)
          nz2 =nz_from_ir(ir2)
          ! find expansion of polynomials in terms of (r-rc)^k 
          cx=0.0d0
          cy=0.0d0
          cz=0.0d0
          call expand_centered_product(r1(1),nx1,r2(1),nx2,rc(1),cx)
          call expand_centered_product(r1(2),ny1,r2(2),ny2,rc(2),cy)
          call expand_centered_product(r1(3),nz1,r2(3),nz2,rc(3),cz)
          ! reconstruct full polynom coeffs 
          p=0
          c1_c2 = c1(ir1)*c2(ir2)
          do l=0,l1max+l2max
            do i=l,0,-1
              c1_c2_cx = c1_c2 * cx(i+1)
              do j=l-i,0,-1
                k=l-i-j
                p=p+1
                cc(p)=cc(p)+c1_c2_cx*cy(j+1)*cz(k+1)
              end do
            end do
          end do
          
        end if
      end do
    end if
  end do
  
  ! set prefactor
  cc(1:irmax_from_lmax(l1max+l2max))=cc(1:irmax_from_lmax(l1max+l2max))*prefac
  
  ! set l3max
  l3max=l1max+l2max
  
end subroutine
  
end module
