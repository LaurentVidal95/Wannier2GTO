
module mod_rlYlm_laplacian

  implicit none

  private
  public :: rlYlm_laplacian

contains

!> compute the laplacian matrix element between two real solid hamronics:
!!  
!!  /
!!  | |r-r2|**l1 * Yl2m2(r-r2) * exp(-a*|r-r2|**2) * ( Nabla**2 *|r-r1|**l1 * Yl1m1(r-r1) * exp(-a*|r-r1|**2) ) *  dr
!!  /
!!  
function rlYlm_laplacian(a1,r1,l1,m1,a2,r2,l2,m2)
  
  implicit none
  
  real(kind=8), intent(in)               :: a1
  real(kind=8), intent(in)               :: a2
  integer     , intent(in)               :: l1
  integer     , intent(in)               :: m1
  integer     , intent(in)               :: l2
  integer     , intent(in)               :: m2
  real(kind=8), intent(in) ,dimension(3) :: r1
  real(kind=8), intent(in) ,dimension(3) :: r2
  
  real(kind=8)                           :: rlYlm_laplacian
  
  real(kind=8)                :: E,d2
  real(kind=8)                :: dx,dy,dz
  
  dx=r1(1)-r2(1)
  dy=r1(2)-r2(2)
  dz=r1(3)-r2(3)
  d2=dx**2+dy**2+dz**2
  
  E=exp(-a1*a2/(a1+a2)*d2)
  
  ! selection on l1
  select case (l1)
    case (0)
      ! selection on m1: l1=0
      select case (m1)
        case (0)
          ! selection on l2: l1=0, m1=0
          select case (l2)
            case (0)
              ! selection on m2: l1=0, m1=0, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =8.86226925452758014d-1*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-4)*(1.0d0*a2*(2.0d0*a1*d2-3.0d0)-3.0d0*a1)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=0, m1=0, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =1.53499006191973273d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-5)*(1.0d0*a2*(2.0d0*a1*d2-5.0d0)-5.0d0*a1)*dy
                case (0)
                  rlYlm_laplacian =1.53499006191973273d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-5)*(1.0d0*a2*(2.0d0*a1*d2-5.0d0)-5.0d0*a1)*dz
                case (1)
                  rlYlm_laplacian =1.53499006191973273d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-5)*(1.0d0*a2*(2.0d0*a1*d2-5.0d0)-5.0d0*a1)*dx
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=0, m1=0, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =3.43234212323913373d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dx*dy
                case (-1)
                  rlYlm_laplacian =3.43234212323913373d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dy*dz
                case (0)
                  rlYlm_laplacian =-9.90831824401502754d-1*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2*(2.0d0*a1*d2-7.0d0) &
-7.0d0*a1)*(3.0d0*(dy**2+dx**2)-2.0d0*d2)
                case (1)
                  rlYlm_laplacian =3.43234212323913373d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dx*dz
                case (2)
                  rlYlm_laplacian = &
-1.71617106161956686d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2* &
(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*(dy+dx)*(dy-1.0d0*dx)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=0, m1=0, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-1.85367660741129178d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2* &
(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dy*(dy**2-3.0d0*dx**2)
                case (-2)
                  rlYlm_laplacian =9.08112367258215863d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*dy*dz
                case (-1)
                  rlYlm_laplacian = &
-1.43585172595163941d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2* &
(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dy*(5.0d0*(dy**2+dx**2)-4.0d0*d2)
                case (0)
                  rlYlm_laplacian = &
-1.17236802495868785d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2* &
(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*(5.0d0*(dy**2+dx**2)-2.0d0*d2)*dz
                case (1)
                  rlYlm_laplacian = &
-1.43585172595163941d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2* &
(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*(5.0d0*(dy**2+dx**2)-4.0d0*d2)
                case (2)
                  rlYlm_laplacian = &
-4.54056183629107931d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2* &
(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*(dy+dx)*(dy-1.0d0*dx)*dz
                case (3)
                  rlYlm_laplacian = &
-1.85367660741129178d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2* &
(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*(3.0d0*dy**2-1.0d0*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=0, m1=0, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-7.86448379536438834d0*E*a1**5*a2*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)
                case (-3)
                  rlYlm_laplacian = &
-5.56102982223387534d0*E*a1**5*a2*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dy*(dy**2-3.0d0*dx**2)*dz
                case (-2)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1**5*a2*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*dy*(7.0d0*(dy**2+dx**2)-6.0d0*d2)
                case (-1)
                  rlYlm_laplacian = &
-2.10187170614922326d0*E*a1**5*a2*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dy*(7.0d0*(dy**2+dx**2)-4.0d0*d2)*dz
                case (0)
                  rlYlm_laplacian =3.32335097044784255d-1*E*a1**5*a2*sqrt &
                  (a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1) &
-1.1d1*a1)*(3.5d1*dy**4+7.0d1*dx**2*dy**2-4.0d1*d2*dy**2+3.5d1*dx**4 &
-4.0d1*d2*dx**2+8.0d0*d2**2)
                case (1)
                  rlYlm_laplacian = &
-2.10187170614922326d0*E*a1**5*a2*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*(7.0d0*(dy**2+dx**2)-4.0d0*d2)*dz
                case (2)
                  rlYlm_laplacian =1.48624773660225413d0*E*a1**5*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*(dy+dx)* &
(dy-1.0d0*dx)*(7.0d0*(dy**2+dx**2)-6.0d0*d2)
                case (3)
                  rlYlm_laplacian = &
-5.56102982223387534d0*E*a1**5*a2*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*(3.0d0*dy**2-1.0d0*dx**2)*dz
                case (4)
                  rlYlm_laplacian =1.96612094884109708d0*E*a1**5*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*(dy**2 &
-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case default
          print *,'Error: rlYlm_overlap not implemented for l1=',l1 &
,'m1=',m1,'l2=',l2,'m2=',m2
          stop
      end select
    case (1)
      ! selection on m1: l1=1
      select case (m1)
        case (-1)
          ! selection on l2: l1=1, m1=-1
          select case (l2)
            case (0)
              ! selection on m2: l1=1, m1=-1, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-1.53499006191973273d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-5)*(1.0d0*a2* &
(2.0d0*a1*d2-5.0d0)-5.0d0*a1)*dy
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=1, m1=-1, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =-1.32934038817913702d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-6)*(4.0d0*a1**2*a2**2*d2*dy**2-1.4d1*a1*a2**2*dy**2 &
-1.4d1*a1**2*a2*dy**2-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+5.0d0*a2**2 &
+1.0d1*a1*a2+5.0d0*a1**2)
                case (0)
                  rlYlm_laplacian = &
-2.65868077635827404d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-6)* &
(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dy*dz
                case (1)
                  rlYlm_laplacian = &
-2.65868077635827404d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-6)* &
(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dx*dy
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=1, m1=-1, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-7)*dx* &
(4.0d0*a1**2*a2**2*d2*dy**2-1.8d1*a1*a2**2*dy**2-1.8d1*a1**2*a2*dy**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2 &
+7.0d0*a1**2)
                case (-1)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-7)* &
(4.0d0*a1**2*a2**2*d2*dy**2-1.8d1*a1*a2**2*dy**2-1.8d1*a1**2*a2*dy**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2 &
+7.0d0*a1**2)*dz
                case (0)
                  rlYlm_laplacian =1.71617106161956686d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-7)*dy*(6.0d0*a1**2*a2**2*d2*dy**2 &
-2.7d1*a1*a2**2*dy**2-2.7d1*a1**2*a2*dy**2+6.0d0*a1**2*a2**2*d2*dx**2 &
-2.7d1*a1*a2**2*dx**2-2.7d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.6d1*a1*a2**2*d2+1.6d1*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2 &
+7.0d0*a1**2)
                case (1)
                  rlYlm_laplacian = &
-5.94499094640901652d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)* &
(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*dy*dz
                case (2)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-7)*dy*(2.0d0*a1**2*a2**2*d2*dy**2 &
-9.0d0*a1*a2**2*dy**2-9.0d0*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+9.0d0*a1*a2**2*dx**2+9.0d0*a1**2*a2*dx**2-2.0d0*a1*a2**2*d2 &
-2.0d0*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2+7.0d0*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=1, m1=-1, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =1.60533103241913232d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(4.0d0*a1**2*a2**2*d2*dy**4-2.2d1*a1*a2**2*dy**4 &
-2.2d1*a1**2*a2*dy**4-1.2d1*a1**2*a2**2*d2*dx**2*dy**2 &
+6.6d1*a1*a2**2*dx**2*dy**2+6.6d1*a1**2*a2*dx**2*dy**2 &
-6.0d0*a1*a2**2*d2*dy**2-6.0d0*a1**2*a2*d2*dy**2+2.7d1*a2**2*dy**2 &
+5.4d1*a1*a2*dy**2+2.7d1*a1**2*dy**2+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-2.7d1*a2**2*dx**2-5.4d1*a1*a2*dx**2 &
-2.7d1*a1**2*dx**2)
                case (-2)
                  rlYlm_laplacian = &
-7.86448379536438834d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-8)*dx* &
(4.0d0*a1**2*a2**2*d2*dy**2-2.2d1*a1*a2**2*dy**2-2.2d1*a1**2*a2*dy**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)*dz
                case (-1)
                  rlYlm_laplacian =1.24348407074185167d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(2.0d1*a1**2*a2**2*d2*dy**4-1.1d2*a1*a2**2*dy**4 &
-1.1d2*a1**2*a2*dy**4+2.0d1*a1**2*a2**2*d2*dx**2*dy**2 &
-1.1d2*a1*a2**2*dx**2*dy**2-1.1d2*a1**2*a2*dx**2*dy**2 &
-1.6d1*a1**2*a2**2*d2**2*dy**2+7.4d1*a1*a2**2*d2*dy**2 &
+7.4d1*a1**2*a2*d2*dy**2+6.3d1*a2**2*dy**2+1.26d2*a1*a2*dy**2 &
+6.3d1*a1**2*dy**2-1.0d1*a1*a2**2*d2*dx**2-1.0d1*a1**2*a2*d2*dx**2 &
+4.5d1*a2**2*dx**2+9.0d1*a1*a2*dx**2+4.5d1*a1**2*dx**2 &
+8.0d0*a1*a2**2*d2**2+8.0d0*a1**2*a2*d2**2-3.6d1*a2**2*d2 &
-7.2d1*a1*a2*d2-3.6d1*a1**2*d2)
                case (0)
                  rlYlm_laplacian =2.03060098439762498d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.6d1*a1*a2**2*d2+1.6d1*a1**2*a2*d2+2.7d1*a2**2+5.4d1*a1*a2 &
+2.7d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =2.48696814148370333d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+4.2d1*a1*a2**2*d2+4.2d1*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)
                case (2)
                  rlYlm_laplacian =7.86448379536438834d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dy*(2.0d0*a1**2*a2**2*d2*dy**2 &
-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2-2.0d0*a1*a2**2*d2 &
-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2+9.0d0*a1**2)*dz
                case (3)
                  rlYlm_laplacian =3.21066206483826464d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(6.0d0*a1**2*a2**2*d2*dy**2 &
-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2-6.0d0*a1*a2**2*d2 &
-6.0d0*a1**2*a2*d2+2.7d1*a2**2+5.4d1*a1*a2+2.7d1*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=1, m1=-1, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =6.81084275443661897d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(4.0d0*a1**2*a2**2*d2*dy**4 &
-2.6d1*a1*a2**2*dy**4-2.6d1*a1**2*a2*dy**4 &
-4.0d0*a1**2*a2**2*d2*dx**2*dy**2+2.6d1*a1*a2**2*dx**2*dy**2 &
+2.6d1*a1**2*a2*dx**2*dy**2-6.0d0*a1*a2**2*d2*dy**2 &
-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2+6.6d1*a1*a2*dy**2 &
+3.3d1*a1**2*dy**2+2.0d0*a1*a2**2*d2*dx**2+2.0d0*a1**2*a2*d2*dx**2 &
-1.1d1*a2**2*dx**2-2.2d1*a1*a2*dx**2-1.1d1*a1**2*dx**2)
                case (-3)
                  rlYlm_laplacian =4.81599309725739696d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(4.0d0*a1**2*a2**2*d2*dy**4-2.6d1*a1*a2**2*dy**4 &
-2.6d1*a1**2*a2*dy**4-1.2d1*a1**2*a2**2*d2*dx**2*dy**2 &
+7.8d1*a1*a2**2*dx**2*dy**2+7.8d1*a1**2*a2*dx**2*dy**2 &
-6.0d0*a1*a2**2*d2*dy**2-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2 &
+6.6d1*a1*a2*dy**2+3.3d1*a1**2*dy**2+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2 &
-3.3d1*a1**2*dx**2)*dz
                case (-2)
                  rlYlm_laplacian =2.5742565924293503d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(2.8d1*a1**2*a2**2*d2*dy**4 &
-1.82d2*a1*a2**2*dy**4-1.82d2*a1**2*a2*dy**4 &
+2.8d1*a1**2*a2**2*d2*dx**2*dy**2-1.82d2*a1*a2**2*dx**2*dy**2 &
-1.82d2*a1**2*a2*dx**2*dy**2-2.4d1*a1**2*a2**2*d2**2*dy**2 &
+1.38d2*a1*a2**2*d2*dy**2+1.38d2*a1**2*a2*d2*dy**2+9.9d1*a2**2*dy**2 &
+1.98d2*a1*a2*dy**2+9.9d1*a1**2*dy**2-1.4d1*a1*a2**2*d2*dx**2 &
-1.4d1*a1**2*a2*d2*dx**2+7.7d1*a2**2*dx**2+1.54d2*a1*a2*dx**2 &
+7.7d1*a1**2*dx**2+1.2d1*a1*a2**2*d2**2+1.2d1*a1**2*a2*d2**2 &
-6.6d1*a2**2*d2-1.32d2*a1*a2*d2-6.6d1*a1**2*d2)
                case (-1)
                  rlYlm_laplacian =1.82027429302096805d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(2.8d1*a1**2*a2**2*d2*dy**4 &
-1.82d2*a1*a2**2*dy**4-1.82d2*a1**2*a2*dy**4 &
+2.8d1*a1**2*a2**2*d2*dx**2*dy**2-1.82d2*a1*a2**2*dx**2*dy**2 &
-1.82d2*a1**2*a2*dx**2*dy**2-1.6d1*a1**2*a2**2*d2**2*dy**2 &
+7.8d1*a1*a2**2*d2*dy**2+7.8d1*a1**2*a2*d2*dy**2+1.43d2*a2**2*dy**2 &
+2.86d2*a1*a2*dy**2+1.43d2*a1**2*dy**2-1.4d1*a1*a2**2*d2*dx**2 &
-1.4d1*a1**2*a2*d2*dx**2+7.7d1*a2**2*dx**2+1.54d2*a1*a2*dx**2 &
+7.7d1*a1**2*dx**2+8.0d0*a1*a2**2*d2**2+8.0d0*a1**2*a2*d2**2 &
-4.4d1*a2**2*d2-8.8d1*a1*a2*d2-4.4d1*a1**2*d2)*dz
                case (0)
                  rlYlm_laplacian =-5.75621273219899775d-1*E*a1**4*a2*sqrt &
                  (a2+a1)*(a2+a1)**(-9)*dy*(7.0d1*a1**2*a2**2*d2*dy**4 &
-4.55d2*a1*a2**2*dy**4-4.55d2*a1**2*a2*dy**4 &
+1.4d2*a1**2*a2**2*d2*dx**2*dy**2-9.1d2*a1*a2**2*dx**2*dy**2 &
-9.1d2*a1**2*a2*dx**2*dy**2-8.0d1*a1**2*a2**2*d2**2*dy**2 &
+4.6d2*a1*a2**2*d2*dy**2+4.6d2*a1**2*a2*d2*dy**2+3.3d2*a2**2*dy**2 &
+6.6d2*a1*a2*dy**2+3.3d2*a1**2*dy**2+7.0d1*a1**2*a2**2*d2*dx**4 &
-4.55d2*a1*a2**2*dx**4-4.55d2*a1**2*a2*dx**4 &
-8.0d1*a1**2*a2**2*d2**2*dx**2+4.6d2*a1*a2**2*d2*dx**2 &
+4.6d2*a1**2*a2*d2*dx**2+3.3d2*a2**2*dx**2+6.6d2*a1*a2*dx**2 &
+3.3d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3-5.6d1*a1*a2**2*d2**2 &
-5.6d1*a1**2*a2*d2**2-2.64d2*a2**2*d2-5.28d2*a1*a2*d2-2.64d2*a1**2*d2)
                case (1)
                  rlYlm_laplacian =3.6405485860419361d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(1.4d1*a1**2*a2**2*d2*dy**2 &
-9.1d1*a1*a2**2*dy**2-9.1d1*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-9.1d1*a1*a2**2*dx**2-9.1d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+4.6d1*a1*a2**2*d2+4.6d1*a1**2*a2*d2+3.3d1*a2**2+6.6d1*a1*a2 &
+3.3d1*a1**2)*dz
                case (2)
                  rlYlm_laplacian = &
-2.5742565924293503d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(1.4d1*a1**2*a2**2*d2*dy**4-9.1d1*a1*a2**2*dy**4-9.1d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2**2*dy**2+6.2d1*a1*a2**2*d2*dy**2 &
+6.2d1*a1**2*a2*d2*dy**2+8.8d1*a2**2*dy**2+1.76d2*a1*a2*dy**2 &
+8.8d1*a1**2*dy**2-1.4d1*a1**2*a2**2*d2*dx**4+9.1d1*a1*a2**2*dx**4 &
+9.1d1*a1**2*a2*dx**4+1.2d1*a1**2*a2**2*d2**2*dx**2 &
-9.0d1*a1*a2**2*d2*dx**2-9.0d1*a1**2*a2*d2*dx**2+6.6d1*a2**2*dx**2 &
+1.32d2*a1*a2*dx**2+6.6d1*a1**2*dx**2+1.2d1*a1*a2**2*d2**2 &
+1.2d1*a1**2*a2*d2**2-6.6d1*a2**2*d2-1.32d2*a1*a2*d2-6.6d1*a1**2*d2)
                case (3)
                  rlYlm_laplacian =9.63198619451479393d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(6.0d0*a1**2*a2**2*d2*dy**2 &
-3.9d1*a1*a2**2*dy**2-3.9d1*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+1.3d1*a1*a2**2*dx**2+1.3d1*a1**2*a2*dx**2-6.0d0*a1*a2**2*d2 &
-6.0d0*a1**2*a2*d2+3.3d1*a2**2+6.6d1*a1*a2+3.3d1*a1**2)*dz
                case (4)
                  rlYlm_laplacian = &
-3.40542137721830949d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(2.0d0*a1**2*a2**2*d2*dy**4-1.3d1*a1*a2**2*dy**4-1.3d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+7.8d1*a1*a2**2*dx**2*dy**2 &
+7.8d1*a1**2*a2*dx**2*dy**2-4.0d0*a1*a2**2*d2*dy**2 &
-4.0d0*a1**2*a2*d2*dy**2+2.2d1*a2**2*dy**2+4.4d1*a1*a2*dy**2 &
+2.2d1*a1**2*dy**2+2.0d0*a1**2*a2**2*d2*dx**4-1.3d1*a1*a2**2*dx**4 &
-1.3d1*a1**2*a2*dx**4+1.2d1*a1*a2**2*d2*dx**2+1.2d1*a1**2*a2*d2*dx**2 &
-6.6d1*a2**2*dx**2-1.32d2*a1*a2*dx**2-6.6d1*a1**2*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (0)
          ! selection on l2: l1=1, m1=0
          select case (l2)
            case (0)
              ! selection on m2: l1=1, m1=0, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-1.53499006191973273d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-5)*(1.0d0*a2* &
(2.0d0*a1*d2-5.0d0)-5.0d0*a1)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=1, m1=0, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-2.65868077635827404d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-6)* &
(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dy*dz
                case (0)
                  rlYlm_laplacian =1.32934038817913702d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-6)*(4.0d0*a1**2*a2**2*d2*dy**2-1.4d1*a1*a2**2*dy**2 &
-1.4d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2-1.4d1*a1*a2**2*dx**2 &
-1.4d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+1.6d1*a1*a2**2*d2 &
+1.6d1*a1**2*a2*d2-5.0d0*a2**2-1.0d1*a1*a2-5.0d0*a1**2)
                case (1)
                  rlYlm_laplacian = &
-2.65868077635827404d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-6)* &
(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dx*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=1, m1=0, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-5.94499094640901652d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)* &
(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*dy*dz
                case (-1)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-7)*dy*(4.0d0*a1**2*a2**2*d2*dy**2 &
-1.8d1*a1*a2**2*dy**2-1.8d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2 &
-1.8d1*a1*a2**2*dx**2-1.8d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+2.0d1*a1*a2**2*d2+2.0d1*a1**2*a2*d2-7.0d0*a2**2-1.4d1*a1*a2 &
-7.0d0*a1**2)
                case (0)
                  rlYlm_laplacian =1.71617106161956686d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-7)*(6.0d0*a1**2*a2**2*d2*dy**2-2.7d1*a1*a2**2*dy**2 &
-2.7d1*a1**2*a2*dy**2+6.0d0*a1**2*a2**2*d2*dx**2-2.7d1*a1*a2**2*dx**2 &
-2.7d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+2.2d1*a1*a2**2*d2 &
+2.2d1*a1**2*a2*d2-1.4d1*a2**2-2.8d1*a1*a2-1.4d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-7)*dx*(4.0d0*a1**2*a2**2*d2*dy**2 &
-1.8d1*a1*a2**2*dy**2-1.8d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2 &
-1.8d1*a1*a2**2*dx**2-1.8d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+2.0d1*a1*a2**2*d2+2.0d1*a1**2*a2*d2-7.0d0*a2**2-1.4d1*a1*a2 &
-7.0d0*a1**2)
                case (2)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*(dy+dx)* &
(dy-1.0d0*dx)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=1, m1=0, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =3.21066206483826464d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dy* &
(dy**2-3.0d0*dx**2)*dz
                case (-2)
                  rlYlm_laplacian =7.86448379536438834d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(4.0d0*a1**2*a2**2*d2*dy**2 &
-2.2d1*a1*a2**2*dy**2-2.2d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2 &
-2.2d1*a1*a2**2*dx**2-2.2d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+2.4d1*a1*a2**2*d2+2.4d1*a1**2*a2*d2-9.0d0*a2**2-1.8d1*a1*a2 &
-9.0d0*a1**2)
                case (-1)
                  rlYlm_laplacian =2.48696814148370333d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+5.2d1*a1*a2**2*d2+5.2d1*a1**2*a2*d2-3.6d1*a2**2-7.2d1*a1*a2 &
-3.6d1*a1**2)*dz
                case (0)
                  rlYlm_laplacian = &
-1.01530049219881249d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-8)* &
(2.0d1*a1**2*a2**2*d2*dy**4-1.1d2*a1*a2**2*dy**4-1.1d2*a1**2*a2*dy**4 &
+4.0d1*a1**2*a2**2*d2*dx**2*dy**2-2.2d2*a1*a2**2*dx**2*dy**2 &
-2.2d2*a1**2*a2*dx**2*dy**2-2.8d1*a1**2*a2**2*d2**2*dy**2 &
+1.72d2*a1*a2**2*d2*dy**2+1.72d2*a1**2*a2*d2*dy**2-8.1d1*a2**2*dy**2 &
-1.62d2*a1*a2*dy**2-8.1d1*a1**2*dy**2+2.0d1*a1**2*a2**2*d2*dx**4 &
-1.1d2*a1*a2**2*dx**4-1.1d2*a1**2*a2*dx**4 &
-2.8d1*a1**2*a2**2*d2**2*dx**2+1.72d2*a1*a2**2*d2*dx**2 &
+1.72d2*a1**2*a2*d2*dx**2-8.1d1*a2**2*dx**2-1.62d2*a1*a2*dx**2 &
-8.1d1*a1**2*dx**2+8.0d0*a1**2*a2**2*d2**3-5.6d1*a1*a2**2*d2**2 &
-5.6d1*a1**2*a2*d2**2+5.4d1*a2**2*d2+1.08d2*a1*a2*d2+5.4d1*a1**2*d2)
                case (1)
                  rlYlm_laplacian =2.48696814148370333d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+5.2d1*a1*a2**2*d2+5.2d1*a1**2*a2*d2-3.6d1*a2**2-7.2d1*a1*a2 &
-3.6d1*a1**2)*dz
                case (2)
                  rlYlm_laplacian = &
-3.93224189768219417d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-8)*(dy+dx)* &
(dy-1.0d0*dx)*(4.0d0*a1**2*a2**2*d2*dy**2-2.2d1*a1*a2**2*dy**2 &
-2.2d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2-2.2d1*a1*a2**2*dx**2 &
-2.2d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+2.4d1*a1*a2**2*d2 &
+2.4d1*a1**2*a2*d2-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)
                case (3)
                  rlYlm_laplacian =3.21066206483826464d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=1, m1=0, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =1.36216855088732379d1*E*a1**5*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*dz
                case (-3)
                  rlYlm_laplacian = &
-4.81599309725739696d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dy*(dy**2 &
-3.0d0*dx**2)*(4.0d0*a1**2*a2**2*d2*dy**2-2.6d1*a1*a2**2*dy**2 &
-2.6d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2-2.6d1*a1*a2**2*dx**2 &
-2.6d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+2.8d1*a1*a2**2*d2 &
+2.8d1*a1**2*a2*d2-1.1d1*a2**2-2.2d1*a1*a2-1.1d1*a1**2)
                case (-2)
                  rlYlm_laplacian =5.14851318485870059d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(1.4d1*a1**2*a2**2*d2*dy**2 &
-9.1d1*a1*a2**2*dy**2-9.1d1*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-9.1d1*a1*a2**2*dx**2-9.1d1*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2 &
+9.0d1*a1*a2**2*d2+9.0d1*a1**2*a2*d2-6.6d1*a2**2-1.32d2*a1*a2 &
-6.6d1*a1**2)*dz
                case (-1)
                  rlYlm_laplacian = &
-1.82027429302096805d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(2.8d1*a1**2*a2**2*d2*dy**4-1.82d2*a1*a2**2*dy**4 &
-1.82d2*a1**2*a2*dy**4+5.6d1*a1**2*a2**2*d2*dx**2*dy**2 &
-3.64d2*a1*a2**2*dx**2*dy**2-3.64d2*a1**2*a2*dx**2*dy**2 &
-4.4d1*a1**2*a2**2*d2**2*dy**2+3.16d2*a1*a2**2*d2*dy**2 &
+3.16d2*a1**2*a2*d2*dy**2-1.65d2*a2**2*dy**2-3.3d2*a1*a2*dy**2 &
-1.65d2*a1**2*dy**2+2.8d1*a1**2*a2**2*d2*dx**4-1.82d2*a1*a2**2*dx**4 &
-1.82d2*a1**2*a2*dx**4-4.4d1*a1**2*a2**2*d2**2*dx**2 &
+3.16d2*a1*a2**2*d2*dx**2+3.16d2*a1**2*a2*d2*dx**2-1.65d2*a2**2*dx**2 &
-3.3d2*a1*a2*dx**2-1.65d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3 &
-1.28d2*a1*a2**2*d2**2-1.28d2*a1**2*a2*d2**2+1.32d2*a2**2*d2 &
+2.64d2*a1*a2*d2+1.32d2*a1**2*d2)
                case (0)
                  rlYlm_laplacian =-5.75621273219899775d-1*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(7.0d1*a1**2*a2**2*d2*dy**4 &
-4.55d2*a1*a2**2*dy**4-4.55d2*a1**2*a2*dy**4 &
+1.4d2*a1**2*a2**2*d2*dx**2*dy**2-9.1d2*a1*a2**2*dx**2*dy**2 &
-9.1d2*a1**2*a2*dx**2*dy**2-8.0d1*a1**2*a2**2*d2**2*dy**2 &
+6.0d2*a1*a2**2*d2*dy**2+6.0d2*a1**2*a2*d2*dy**2-4.4d2*a2**2*dy**2 &
-8.8d2*a1*a2*dy**2-4.4d2*a1**2*dy**2+7.0d1*a1**2*a2**2*d2*dx**4 &
-4.55d2*a1*a2**2*dx**4-4.55d2*a1**2*a2*dx**4 &
-8.0d1*a1**2*a2**2*d2**2*dx**2+6.0d2*a1*a2**2*d2*dx**2 &
+6.0d2*a1**2*a2*d2*dx**2-4.4d2*a2**2*dx**2-8.8d2*a1*a2*dx**2 &
-4.4d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3-1.36d2*a1*a2**2*d2**2 &
-1.36d2*a1**2*a2*d2**2+1.76d2*a2**2*d2+3.52d2*a1*a2*d2+1.76d2*a1**2*d2 &
)*dz
                case (1)
                  rlYlm_laplacian = &
-1.82027429302096805d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(2.8d1*a1**2*a2**2*d2*dy**4-1.82d2*a1*a2**2*dy**4 &
-1.82d2*a1**2*a2*dy**4+5.6d1*a1**2*a2**2*d2*dx**2*dy**2 &
-3.64d2*a1*a2**2*dx**2*dy**2-3.64d2*a1**2*a2*dx**2*dy**2 &
-4.4d1*a1**2*a2**2*d2**2*dy**2+3.16d2*a1*a2**2*d2*dy**2 &
+3.16d2*a1**2*a2*d2*dy**2-1.65d2*a2**2*dy**2-3.3d2*a1*a2*dy**2 &
-1.65d2*a1**2*dy**2+2.8d1*a1**2*a2**2*d2*dx**4-1.82d2*a1*a2**2*dx**4 &
-1.82d2*a1**2*a2*dx**4-4.4d1*a1**2*a2**2*d2**2*dx**2 &
+3.16d2*a1*a2**2*d2*dx**2+3.16d2*a1**2*a2*d2*dx**2-1.65d2*a2**2*dx**2 &
-3.3d2*a1*a2*dx**2-1.65d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3 &
-1.28d2*a1*a2**2*d2**2-1.28d2*a1**2*a2*d2**2+1.32d2*a2**2*d2 &
+2.64d2*a1*a2*d2+1.32d2*a1**2*d2)
                case (2)
                  rlYlm_laplacian = &
-2.5742565924293503d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*(dy+dx)*(dy &
-1.0d0*dx)*(1.4d1*a1**2*a2**2*d2*dy**2-9.1d1*a1*a2**2*dy**2 &
-9.1d1*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2-9.1d1*a1*a2**2*dx**2 &
-9.1d1*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2+9.0d1*a1*a2**2*d2 &
+9.0d1*a1**2*a2*d2-6.6d1*a2**2-1.32d2*a1*a2-6.6d1*a1**2)*dz
                case (3)
                  rlYlm_laplacian = &
-4.81599309725739696d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*(4.0d0*a1**2*a2**2*d2*dy**2 &
-2.6d1*a1*a2**2*dy**2-2.6d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2 &
-2.6d1*a1*a2**2*dx**2-2.6d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+2.8d1*a1*a2**2*d2+2.8d1*a1**2*a2*d2-1.1d1*a2**2-2.2d1*a1*a2 &
-1.1d1*a1**2)
                case (4)
                  rlYlm_laplacian = &
-3.40542137721830949d0*E*a1**5*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)* &
(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*(dy**2-2.0d0*dx*dy-1.0d0*dx**2 &
)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (1)
          ! selection on l2: l1=1, m1=1
          select case (l2)
            case (0)
              ! selection on m2: l1=1, m1=1, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-1.53499006191973273d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-5)*(1.0d0*a2* &
(2.0d0*a1*d2-5.0d0)-5.0d0*a1)*dx
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=1, m1=1, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-2.65868077635827404d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-6)* &
(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dx*dy
                case (0)
                  rlYlm_laplacian = &
-2.65868077635827404d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-6)* &
(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dx*dz
                case (1)
                  rlYlm_laplacian =-1.32934038817913702d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-6)*(4.0d0*a1**2*a2**2*d2*dx**2-1.4d1*a1*a2**2*dx**2 &
-1.4d1*a1**2*a2*dx**2-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+5.0d0*a2**2 &
+1.0d1*a1*a2+5.0d0*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=1, m1=1, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-7)* &
(4.0d0*a1**2*a2**2*d2*dx**2-1.8d1*a1*a2**2*dx**2-1.8d1*a1**2*a2*dx**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2 &
+7.0d0*a1**2)*dy
                case (-1)
                  rlYlm_laplacian = &
-5.94499094640901652d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)* &
(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*dy*dz
                case (0)
                  rlYlm_laplacian =1.71617106161956686d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-7)*dx*(6.0d0*a1**2*a2**2*d2*dy**2 &
-2.7d1*a1*a2**2*dy**2-2.7d1*a1**2*a2*dy**2+6.0d0*a1**2*a2**2*d2*dx**2 &
-2.7d1*a1*a2**2*dx**2-2.7d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.6d1*a1*a2**2*d2+1.6d1*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2 &
+7.0d0*a1**2)
                case (1)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-7)* &
(4.0d0*a1**2*a2**2*d2*dx**2-1.8d1*a1*a2**2*dx**2-1.8d1*a1**2*a2*dx**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2 &
+7.0d0*a1**2)*dz
                case (2)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-7)*dx*(2.0d0*a1**2*a2**2*d2*dy**2 &
-9.0d0*a1*a2**2*dy**2-9.0d0*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+9.0d0*a1*a2**2*dx**2+9.0d0*a1**2*a2*dx**2+2.0d0*a1*a2**2*d2 &
+2.0d0*a1**2*a2*d2-7.0d0*a2**2-1.4d1*a1*a2-7.0d0*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=1, m1=1, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =3.21066206483826464d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(2.0d0*a1**2*a2**2*d2*dy**2 &
-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2-6.0d0*a1**2*a2**2*d2*dx**2 &
+3.3d1*a1*a2**2*dx**2+3.3d1*a1**2*a2*dx**2+6.0d0*a1*a2**2*d2 &
+6.0d0*a1**2*a2*d2-2.7d1*a2**2-5.4d1*a1*a2-2.7d1*a1**2)
                case (-2)
                  rlYlm_laplacian = &
-7.86448379536438834d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-8)* &
(4.0d0*a1**2*a2**2*d2*dx**2-2.2d1*a1*a2**2*dx**2-2.2d1*a1**2*a2*dx**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)*dy*dz
                case (-1)
                  rlYlm_laplacian =2.48696814148370333d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+4.2d1*a1*a2**2*d2+4.2d1*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)
                case (0)
                  rlYlm_laplacian =2.03060098439762498d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.6d1*a1*a2**2*d2+1.6d1*a1**2*a2*d2+2.7d1*a2**2+5.4d1*a1*a2 &
+2.7d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =1.24348407074185167d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(2.0d1*a1**2*a2**2*d2*dx**2*dy**2 &
-1.1d2*a1*a2**2*dx**2*dy**2-1.1d2*a1**2*a2*dx**2*dy**2 &
-1.0d1*a1*a2**2*d2*dy**2-1.0d1*a1**2*a2*d2*dy**2+4.5d1*a2**2*dy**2 &
+9.0d1*a1*a2*dy**2+4.5d1*a1**2*dy**2+2.0d1*a1**2*a2**2*d2*dx**4 &
-1.1d2*a1*a2**2*dx**4-1.1d2*a1**2*a2*dx**4 &
-1.6d1*a1**2*a2**2*d2**2*dx**2+7.4d1*a1*a2**2*d2*dx**2 &
+7.4d1*a1**2*a2*d2*dx**2+6.3d1*a2**2*dx**2+1.26d2*a1*a2*dx**2 &
+6.3d1*a1**2*dx**2+8.0d0*a1*a2**2*d2**2+8.0d0*a1**2*a2*d2**2 &
-3.6d1*a2**2*d2-7.2d1*a1*a2*d2-3.6d1*a1**2*d2)
                case (2)
                  rlYlm_laplacian =7.86448379536438834d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*(2.0d0*a1**2*a2**2*d2*dy**2 &
-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2+2.0d0*a1*a2**2*d2 &
+2.0d0*a1**2*a2*d2-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)*dz
                case (3)
                  rlYlm_laplacian =1.60533103241913232d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.2d1*a1**2*a2**2*d2*dx**2*dy**2 &
-6.6d1*a1*a2**2*dx**2*dy**2-6.6d1*a1**2*a2*dx**2*dy**2 &
-6.0d0*a1*a2**2*d2*dy**2-6.0d0*a1**2*a2*d2*dy**2+2.7d1*a2**2*dy**2 &
+5.4d1*a1*a2*dy**2+2.7d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4 &
+2.2d1*a1*a2**2*dx**4+2.2d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-2.7d1*a2**2*dx**2-5.4d1*a1*a2*dx**2 &
-2.7d1*a1**2*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=1, m1=1, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =6.81084275443661897d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(4.0d0*a1**2*a2**2*d2*dx**2*dy**2 &
-2.6d1*a1*a2**2*dx**2*dy**2-2.6d1*a1**2*a2*dx**2*dy**2 &
-2.0d0*a1*a2**2*d2*dy**2-2.0d0*a1**2*a2*d2*dy**2+1.1d1*a2**2*dy**2 &
+2.2d1*a1*a2*dy**2+1.1d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4 &
+2.6d1*a1*a2**2*dx**4+2.6d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2 &
-3.3d1*a1**2*dx**2)
                case (-3)
                  rlYlm_laplacian =9.63198619451479393d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(2.0d0*a1**2*a2**2*d2*dy**2 &
-1.3d1*a1*a2**2*dy**2-1.3d1*a1**2*a2*dy**2-6.0d0*a1**2*a2**2*d2*dx**2 &
+3.9d1*a1*a2**2*dx**2+3.9d1*a1**2*a2*dx**2+6.0d0*a1*a2**2*d2 &
+6.0d0*a1**2*a2*d2-3.3d1*a2**2-6.6d1*a1*a2-3.3d1*a1**2)*dz
                case (-2)
                  rlYlm_laplacian =2.5742565924293503d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(2.8d1*a1**2*a2**2*d2*dx**2*dy**2 &
-1.82d2*a1*a2**2*dx**2*dy**2-1.82d2*a1**2*a2*dx**2*dy**2 &
-1.4d1*a1*a2**2*d2*dy**2-1.4d1*a1**2*a2*d2*dy**2+7.7d1*a2**2*dy**2 &
+1.54d2*a1*a2*dy**2+7.7d1*a1**2*dy**2+2.8d1*a1**2*a2**2*d2*dx**4 &
-1.82d2*a1*a2**2*dx**4-1.82d2*a1**2*a2*dx**4 &
-2.4d1*a1**2*a2**2*d2**2*dx**2+1.38d2*a1*a2**2*d2*dx**2 &
+1.38d2*a1**2*a2*d2*dx**2+9.9d1*a2**2*dx**2+1.98d2*a1*a2*dx**2 &
+9.9d1*a1**2*dx**2+1.2d1*a1*a2**2*d2**2+1.2d1*a1**2*a2*d2**2 &
-6.6d1*a2**2*d2-1.32d2*a1*a2*d2-6.6d1*a1**2*d2)
                case (-1)
                  rlYlm_laplacian =3.6405485860419361d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(1.4d1*a1**2*a2**2*d2*dy**2 &
-9.1d1*a1*a2**2*dy**2-9.1d1*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-9.1d1*a1*a2**2*dx**2-9.1d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+4.6d1*a1*a2**2*d2+4.6d1*a1**2*a2*d2+3.3d1*a2**2+6.6d1*a1*a2 &
+3.3d1*a1**2)*dz
                case (0)
                  rlYlm_laplacian =-5.75621273219899775d-1*E*a1**4*a2*sqrt &
                  (a2+a1)*(a2+a1)**(-9)*dx*(7.0d1*a1**2*a2**2*d2*dy**4 &
-4.55d2*a1*a2**2*dy**4-4.55d2*a1**2*a2*dy**4 &
+1.4d2*a1**2*a2**2*d2*dx**2*dy**2-9.1d2*a1*a2**2*dx**2*dy**2 &
-9.1d2*a1**2*a2*dx**2*dy**2-8.0d1*a1**2*a2**2*d2**2*dy**2 &
+4.6d2*a1*a2**2*d2*dy**2+4.6d2*a1**2*a2*d2*dy**2+3.3d2*a2**2*dy**2 &
+6.6d2*a1*a2*dy**2+3.3d2*a1**2*dy**2+7.0d1*a1**2*a2**2*d2*dx**4 &
-4.55d2*a1*a2**2*dx**4-4.55d2*a1**2*a2*dx**4 &
-8.0d1*a1**2*a2**2*d2**2*dx**2+4.6d2*a1*a2**2*d2*dx**2 &
+4.6d2*a1**2*a2*d2*dx**2+3.3d2*a2**2*dx**2+6.6d2*a1*a2*dx**2 &
+3.3d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3-5.6d1*a1*a2**2*d2**2 &
-5.6d1*a1**2*a2*d2**2-2.64d2*a2**2*d2-5.28d2*a1*a2*d2-2.64d2*a1**2*d2)
                case (1)
                  rlYlm_laplacian =1.82027429302096805d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(2.8d1*a1**2*a2**2*d2*dx**2*dy**2 &
-1.82d2*a1*a2**2*dx**2*dy**2-1.82d2*a1**2*a2*dx**2*dy**2 &
-1.4d1*a1*a2**2*d2*dy**2-1.4d1*a1**2*a2*d2*dy**2+7.7d1*a2**2*dy**2 &
+1.54d2*a1*a2*dy**2+7.7d1*a1**2*dy**2+2.8d1*a1**2*a2**2*d2*dx**4 &
-1.82d2*a1*a2**2*dx**4-1.82d2*a1**2*a2*dx**4 &
-1.6d1*a1**2*a2**2*d2**2*dx**2+7.8d1*a1*a2**2*d2*dx**2 &
+7.8d1*a1**2*a2*d2*dx**2+1.43d2*a2**2*dx**2+2.86d2*a1*a2*dx**2 &
+1.43d2*a1**2*dx**2+8.0d0*a1*a2**2*d2**2+8.0d0*a1**2*a2*d2**2 &
-4.4d1*a2**2*d2-8.8d1*a1*a2*d2-4.4d1*a1**2*d2)*dz
                case (2)
                  rlYlm_laplacian = &
-2.5742565924293503d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(1.4d1*a1**2*a2**2*d2*dy**4-9.1d1*a1*a2**2*dy**4-9.1d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2**2*dy**2+9.0d1*a1*a2**2*d2*dy**2 &
+9.0d1*a1**2*a2*d2*dy**2-6.6d1*a2**2*dy**2-1.32d2*a1*a2*dy**2 &
-6.6d1*a1**2*dy**2-1.4d1*a1**2*a2**2*d2*dx**4+9.1d1*a1*a2**2*dx**4 &
+9.1d1*a1**2*a2*dx**4+1.2d1*a1**2*a2**2*d2**2*dx**2 &
-6.2d1*a1*a2**2*d2*dx**2-6.2d1*a1**2*a2*d2*dx**2-8.8d1*a2**2*dx**2 &
-1.76d2*a1*a2*dx**2-8.8d1*a1**2*dx**2-1.2d1*a1*a2**2*d2**2 &
-1.2d1*a1**2*a2*d2**2+6.6d1*a2**2*d2+1.32d2*a1*a2*d2+6.6d1*a1**2*d2)
                case (3)
                  rlYlm_laplacian =4.81599309725739696d0*E*a1**4*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(1.2d1*a1**2*a2**2*d2*dx**2*dy**2 &
-7.8d1*a1*a2**2*dx**2*dy**2-7.8d1*a1**2*a2*dx**2*dy**2 &
-6.0d0*a1*a2**2*d2*dy**2-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2 &
+6.6d1*a1*a2*dy**2+3.3d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4 &
+2.6d1*a1*a2**2*dx**4+2.6d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2 &
-3.3d1*a1**2*dx**2)*dz
                case (4)
                  rlYlm_laplacian = &
-3.40542137721830949d0*E*a1**4*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(2.0d0*a1**2*a2**2*d2*dy**4-1.3d1*a1*a2**2*dy**4-1.3d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+7.8d1*a1*a2**2*dx**2*dy**2 &
+7.8d1*a1**2*a2*dx**2*dy**2+1.2d1*a1*a2**2*d2*dy**2 &
+1.2d1*a1**2*a2*d2*dy**2-6.6d1*a2**2*dy**2-1.32d2*a1*a2*dy**2 &
-6.6d1*a1**2*dy**2+2.0d0*a1**2*a2**2*d2*dx**4-1.3d1*a1*a2**2*dx**4 &
-1.3d1*a1**2*a2*dx**4-4.0d0*a1*a2**2*d2*dx**2-4.0d0*a1**2*a2*d2*dx**2 &
+2.2d1*a2**2*dx**2+4.4d1*a1*a2*dx**2+2.2d1*a1**2*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case default
          print *,'Error: rlYlm_overlap not implemented for l1=',l1 &
,'m1=',m1,'l2=',l2,'m2=',m2
          stop
      end select
    case (2)
      ! selection on m1: l1=2
      select case (m1)
        case (-2)
          ! selection on l2: l1=2, m1=-2
          select case (l2)
            case (0)
              ! selection on m2: l1=2, m1=-2, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =3.43234212323913373d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dx*dy
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=2, m1=-2, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-7)*dx*(4.0d0*a1**2*a2**2*d2*dy**2 &
-1.8d1*a1*a2**2*dy**2-1.8d1*a1**2*a2*dy**2-2.0d0*a1*a2**2*d2 &
-2.0d0*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2+7.0d0*a1**2)
                case (0)
                  rlYlm_laplacian =5.94499094640901652d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*dy*dz
                case (1)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-7)*(4.0d0*a1**2*a2**2*d2*dx**2-1.8d1*a1*a2**2*dx**2 &
-1.8d1*a1**2*a2*dx**2-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+7.0d0*a2**2 &
+1.4d1*a1*a2+7.0d0*a1**2)*dy
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=2, m1=-2, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =3.32335097044784255d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
-4.4d1*a1**2*a2**3*dx**2*dy**2-4.4d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+1.8d1*a1*a2**3*dy**2+3.6d1*a1**2*a2**2*dy**2+1.8d1*a1**3*a2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+1.8d1*a1*a2**3*dx**2+3.6d1*a1**2*a2**2*dx**2+1.8d1*a1**3*a2*dx**2 &
+2.0d0*a1*a2**3*d2+4.0d0*a1**2*a2**2*d2+2.0d0*a1**3*a2*d2-7.0d0*a2**3 &
-2.1d1*a1*a2**2-2.1d1*a1**2*a2-7.0d0*a1**3)
                case (-1)
                  rlYlm_laplacian =6.6467019408956851d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*(4.0d0*a1**2*a2**2*d2*dy**2 &
-2.2d1*a1*a2**2*dy**2-2.2d1*a1**2*a2*dy**2-2.0d0*a1*a2**2*d2 &
-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2+9.0d0*a1**2)*dz
                case (0)
                  rlYlm_laplacian = &
-3.83747515479933183d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dx*dy* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.3d1*a1*a2**2*dx**2-3.3d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.8d1*a1*a2**2*d2+1.8d1*a1**2*a2*d2 &
+1.8d1*a2**2+3.6d1*a1*a2+1.8d1*a1**2)
                case (1)
                  rlYlm_laplacian =6.6467019408956851d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(4.0d0*a1**2*a2**2*d2*dx**2-2.2d1*a1*a2**2*dx**2 &
-2.2d1*a1**2*a2*dx**2-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+9.0d0*a2**2 &
+1.8d1*a1*a2+9.0d0*a1**2)*dy*dz
                case (2)
                  rlYlm_laplacian = &
-6.6467019408956851d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-8)* &
(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=2, m1=-2, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-3.58962931487909853d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(4.0d0*a1**3*a2**3*d2*dy**4-2.6d1*a1**2*a2**3*dy**4 &
-2.6d1*a1**3*a2**2*dy**4-1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
+7.8d1*a1**2*a2**3*dx**2*dy**2+7.8d1*a1**3*a2**2*dx**2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dy**2+6.0d0*a1**3*a2**2*d2*dy**2 &
-3.3d1*a1*a2**3*dy**2-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dx**2+6.0d0*a1**3*a2**2*d2*dx**2 &
-3.3d1*a1*a2**3*dx**2-6.6d1*a1**2*a2**2*dx**2-3.3d1*a1**3*a2*dx**2 &
-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2-6.0d0*a1**3*a2*d2+2.7d1*a2**3 &
+8.1d1*a1*a2**2+8.1d1*a1**2*a2+2.7d1*a1**3)
                case (-2)
                  rlYlm_laplacian =8.79276018719015889d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
-5.2d1*a1**2*a2**3*dx**2*dy**2-5.2d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+2.2d1*a1*a2**3*dy**2+4.4d1*a1**2*a2**2*dy**2+2.2d1*a1**3*a2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+2.2d1*a1*a2**3*dx**2+4.4d1*a1**2*a2**2*dx**2+2.2d1*a1**3*a2*dx**2 &
+2.0d0*a1*a2**3*d2+4.0d0*a1**2*a2**2*d2+2.0d0*a1**3*a2*d2-9.0d0*a2**3 &
-2.7d1*a1*a2**2-2.7d1*a1**2*a2-9.0d0*a1**3)*dz
                case (-1)
                  rlYlm_laplacian = &
-2.78051491111693767d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.3d2*a1**2*a2**3*dy**4 &
-1.3d2*a1**3*a2**2*dy**4+2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.3d2*a1**2*a2**3*dx**2*dy**2-1.3d2*a1**3*a2**2*dx**2*dy**2 &
-1.6d1*a1**3*a2**3*d2**2*dy**2+8.6d1*a1**2*a2**3*d2*dy**2 &
+8.6d1*a1**3*a2**2*d2*dy**2+9.9d1*a1*a2**3*dy**2 &
+1.98d2*a1**2*a2**2*dy**2+9.9d1*a1**3*a2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dx**2-1.0d1*a1**3*a2**2*d2*dx**2 &
+5.5d1*a1*a2**3*dx**2+1.1d2*a1**2*a2**2*dx**2+5.5d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.2d1*a1*a2**3*d2 &
-8.4d1*a1**2*a2**2*d2-4.2d1*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)
                case (0)
                  rlYlm_laplacian = &
-4.54056183629107931d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(1.0d1*a1**2*a2**2*d2*dy**2-6.5d1*a1*a2**2*dy**2-6.5d1*a1**2*a2*dy**2 &
+1.0d1*a1**2*a2**2*d2*dx**2-6.5d1*a1*a2**2*dx**2-6.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.4d1*a1*a2**2*d2+1.4d1*a1**2*a2*d2 &
+6.6d1*a2**2+1.32d2*a1*a2+6.6d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian = &
-2.78051491111693767d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(2.0d1*a1**3*a2**3*d2*dx**2*dy**2-1.3d2*a1**2*a2**3*dx**2*dy**2 &
-1.3d2*a1**3*a2**2*dx**2*dy**2-1.0d1*a1**2*a2**3*d2*dy**2 &
-1.0d1*a1**3*a2**2*d2*dy**2+5.5d1*a1*a2**3*dy**2 &
+1.1d2*a1**2*a2**2*dy**2+5.5d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.3d2*a1**2*a2**3*dx**4 &
-1.3d2*a1**3*a2**2*dx**4-1.6d1*a1**3*a2**3*d2**2*dx**2 &
+8.6d1*a1**2*a2**3*d2*dx**2+8.6d1*a1**3*a2**2*d2*dx**2 &
+9.9d1*a1*a2**3*dx**2+1.98d2*a1**2*a2**2*dx**2+9.9d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.2d1*a1*a2**3*d2 &
-8.4d1*a1**2*a2**2*d2-4.2d1*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)
                case (2)
                  rlYlm_laplacian = &
-1.75855203743803178d1*E*a1**4*a2**3*sqrt(a2+a1)*(a2+a1)**(-9)* &
(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)*dz
                case (3)
                  rlYlm_laplacian = &
-3.58962931487909853d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(1.2d1*a1**3*a2**3*d2*dx**2*dy**2-7.8d1*a1**2*a2**3*dx**2*dy**2 &
-7.8d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.3d1*a1*a2**3*dy**2 &
+6.6d1*a1**2*a2**2*dy**2+3.3d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+2.6d1*a1**2*a2**3*dx**4 &
+2.6d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.3d1*a1*a2**3*dx**2 &
+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-2.7d1*a2**3-8.1d1*a1*a2**2 &
-8.1d1*a1**2*a2-2.7d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=2, m1=-2, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-7.6147536914910937d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*(dy+dx)* &
(dy-1.0d0*dx)*(8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
-6.0d1*a1**2*a2**3*dx**2*dy**2-6.0d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+2.6d1*a1*a2**3*dy**2+5.2d1*a1**2*a2**2*dy**2+2.6d1*a1**3*a2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+2.6d1*a1*a2**3*dx**2+5.2d1*a1**2*a2**2*dx**2+2.6d1*a1**3*a2*dx**2 &
+6.0d0*a1*a2**3*d2+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3 &
-9.9d1*a1*a2**2-9.9d1*a1**2*a2-3.3d1*a1**3)
                case (-3)
                  rlYlm_laplacian = &
-1.07688879446372956d1*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(4.0d0*a1**3*a2**3*d2*dy**4-3.0d1*a1**2*a2**3*dy**4 &
-3.0d1*a1**3*a2**2*dy**4-1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**3*dx**2*dy**2+9.0d1*a1**3*a2**2*dx**2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dy**2+6.0d0*a1**3*a2**2*d2*dy**2 &
-3.9d1*a1*a2**3*dy**2-7.8d1*a1**2*a2**2*dy**2-3.9d1*a1**3*a2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dx**2+6.0d0*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2-6.0d0*a1**3*a2*d2+3.3d1*a2**3 &
+9.9d1*a1*a2**2+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (-2)
                  rlYlm_laplacian = &
-2.87810636609949888d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)* &
(5.6d1*a1**3*a2**3*d2*dx**2*dy**4-4.2d2*a1**2*a2**3*dx**2*dy**4 &
-4.2d2*a1**3*a2**2*dx**2*dy**4-2.8d1*a1**2*a2**3*d2*dy**4 &
-2.8d1*a1**3*a2**2*d2*dy**4+1.82d2*a1*a2**3*dy**4 &
+3.64d2*a1**2*a2**2*dy**4+1.82d2*a1**3*a2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**4*dy**2-4.2d2*a1**2*a2**3*dx**4*dy**2 &
-4.2d2*a1**3*a2**2*dx**4*dy**2-4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+2.88d2*a1**2*a2**3*d2*dx**2*dy**2+2.88d2*a1**3*a2**2*d2*dx**2*dy**2 &
+4.68d2*a1*a2**3*dx**2*dy**2+9.36d2*a1**2*a2**2*dx**2*dy**2 &
+4.68d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.38d2*a1*a2**3*d2*dy**2 &
-2.76d2*a1**2*a2**2*d2*dy**2-1.38d2*a1**3*a2*d2*dy**2 &
-9.9d1*a2**3*dy**2-2.97d2*a1*a2**2*dy**2-2.97d2*a1**2*a2*dy**2 &
-9.9d1*a1**3*dy**2-2.8d1*a1**2*a2**3*d2*dx**4 &
-2.8d1*a1**3*a2**2*d2*dx**4+1.82d2*a1*a2**3*dx**4 &
+3.64d2*a1**2*a2**2*dx**4+1.82d2*a1**3*a2*dx**4 &
+2.4d1*a1**2*a2**3*d2**2*dx**2+2.4d1*a1**3*a2**2*d2**2*dx**2 &
-1.38d2*a1*a2**3*d2*dx**2-2.76d2*a1**2*a2**2*d2*dx**2 &
-1.38d2*a1**3*a2*d2*dx**2-9.9d1*a2**3*dx**2-2.97d2*a1*a2**2*dx**2 &
-2.97d2*a1**2*a2*dx**2-9.9d1*a1**3*dx**2-1.2d1*a1*a2**3*d2**2 &
-2.4d1*a1**2*a2**2*d2**2-1.2d1*a1**3*a2*d2**2+6.6d1*a2**3*d2 &
+1.98d2*a1*a2**2*d2+1.98d2*a1**2*a2*d2+6.6d1*a1**3*d2)
                case (-1)
                  rlYlm_laplacian = &
-4.07025705689025559d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(2.8d1*a1**3*a2**3*d2*dy**4-2.1d2*a1**2*a2**3*dy**4 &
-2.1d2*a1**3*a2**2*dy**4+2.8d1*a1**3*a2**3*d2*dx**2*dy**2 &
-2.1d2*a1**2*a2**3*dx**2*dy**2-2.1d2*a1**3*a2**2*dx**2*dy**2 &
-1.6d1*a1**3*a2**3*d2**2*dy**2+8.2d1*a1**2*a2**3*d2*dy**2 &
+8.2d1*a1**3*a2**2*d2*dy**2+2.47d2*a1*a2**3*dy**2 &
+4.94d2*a1**2*a2**2*dy**2+2.47d2*a1**3*a2*dy**2 &
-1.4d1*a1**2*a2**3*d2*dx**2-1.4d1*a1**3*a2**2*d2*dx**2 &
+9.1d1*a1*a2**3*dx**2+1.82d2*a1**2*a2**2*dx**2+9.1d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.6d1*a1*a2**3*d2 &
-9.2d1*a1**2*a2**2*d2-4.6d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian =1.28712829621467515d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.25d2*a1**2*a2**3*dy**4-5.25d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.05d3*a1**2*a2**3*dx**2*dy**2 &
-1.05d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+4.8d2*a1**2*a2**3*d2*dy**2+4.8d2*a1**3*a2**2*d2*dy**2 &
+7.8d2*a1*a2**3*dy**2+1.56d3*a1**2*a2**2*dy**2+7.8d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.25d2*a1**2*a2**3*dx**4 &
-5.25d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+4.8d2*a1**2*a2**3*d2*dx**2+4.8d2*a1**3*a2**2*d2*dx**2 &
+7.8d2*a1*a2**3*dx**2+1.56d3*a1**2*a2**2*dx**2+7.8d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-2.4d1*a1**2*a2**3*d2**2 &
-2.4d1*a1**3*a2**2*d2**2-6.12d2*a1*a2**3*d2-1.224d3*a1**2*a2**2*d2 &
-6.12d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2-1.98d2*a1**2*a2 &
-6.6d1*a1**3)
                case (1)
                  rlYlm_laplacian = &
-4.07025705689025559d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(2.8d1*a1**3*a2**3*d2*dx**2*dy**2-2.1d2*a1**2*a2**3*dx**2*dy**2 &
-2.1d2*a1**3*a2**2*dx**2*dy**2-1.4d1*a1**2*a2**3*d2*dy**2 &
-1.4d1*a1**3*a2**2*d2*dy**2+9.1d1*a1*a2**3*dy**2 &
+1.82d2*a1**2*a2**2*dy**2+9.1d1*a1**3*a2*dy**2 &
+2.8d1*a1**3*a2**3*d2*dx**4-2.1d2*a1**2*a2**3*dx**4 &
-2.1d2*a1**3*a2**2*dx**4-1.6d1*a1**3*a2**3*d2**2*dx**2 &
+8.2d1*a1**2*a2**3*d2*dx**2+8.2d1*a1**3*a2**2*d2*dx**2 &
+2.47d2*a1*a2**3*dx**2+4.94d2*a1**2*a2**2*dx**2+2.47d2*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.6d1*a1*a2**3*d2 &
-9.2d1*a1**2*a2**2*d2-4.6d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian =5.75621273219899775d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**2*a2**2*d2*dy**2-1.05d2*a1*a2**2*dy**2 &
-1.05d2*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-1.05d2*a1*a2**2*dx**2-1.05d2*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2 &
+8.6d1*a1*a2**2*d2+8.6d1*a1**2*a2*d2+2.6d1*a2**2+5.2d1*a1*a2 &
+2.6d1*a1**2)
                case (3)
                  rlYlm_laplacian = &
-1.07688879446372956d1*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(1.2d1*a1**3*a2**3*d2*dx**2*dy**2-9.0d1*a1**2*a2**3*dx**2*dy**2 &
-9.0d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+3.0d1*a1**2*a2**3*dx**4 &
+3.0d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (4)
                  rlYlm_laplacian =7.6147536914910937d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.5d1*a1**2*a2**3*dy**4-1.5d1*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2*dx**2*dy**2+9.0d1*a1**2*a2**3*dx**2*dy**2 &
+9.0d1*a1**3*a2**2*dx**2*dy**2+8.0d0*a1**2*a2**3*d2*dy**2 &
+8.0d0*a1**3*a2**2*d2*dy**2-5.2d1*a1*a2**3*dy**2 &
-1.04d2*a1**2*a2**2*dy**2-5.2d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4+8.0d0*a1**2*a2**3*d2*dx**2 &
+8.0d0*a1**3*a2**2*d2*dx**2-5.2d1*a1*a2**3*dx**2 &
-1.04d2*a1**2*a2**2*dx**2-5.2d1*a1**3*a2*dx**2-1.2d1*a1*a2**3*d2 &
-2.4d1*a1**2*a2**2*d2-1.2d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (-1)
          ! selection on l2: l1=2, m1=-1
          select case (l2)
            case (0)
              ! selection on m2: l1=2, m1=-1, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =3.43234212323913373d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dy*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=2, m1=-1, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-7)*(4.0d0*a1**2*a2**2*d2*dy**2-1.8d1*a1*a2**2*dy**2 &
-1.8d1*a1**2*a2*dy**2-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+7.0d0*a2**2 &
+1.4d1*a1*a2+7.0d0*a1**2)*dz
                case (0)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)*dy* &
(4.0d0*a1**2*a2**2*d2*dy**2-1.8d1*a1*a2**2*dy**2-1.8d1*a1**2*a2*dy**2 &
+4.0d0*a1**2*a2**2*d2*dx**2-1.8d1*a1*a2**2*dx**2-1.8d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.0d1*a1*a2**2*d2+2.0d1*a1**2*a2*d2 &
-7.0d0*a2**2-1.4d1*a1*a2-7.0d0*a1**2)
                case (1)
                  rlYlm_laplacian =5.94499094640901652d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*dy*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=2, m1=-1, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =6.6467019408956851d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*(4.0d0*a1**2*a2**2*d2*dy**2 &
-2.2d1*a1*a2**2*dy**2-2.2d1*a1**2*a2*dy**2-2.0d0*a1*a2**2*d2 &
-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2+9.0d0*a1**2)*dz
                case (-1)
                  rlYlm_laplacian =-3.32335097044784255d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(8.0d0*a1**3*a2**3*d2*dy**4 &
-4.4d1*a1**2*a2**3*dy**4-4.4d1*a1**3*a2**2*dy**4 &
+8.0d0*a1**3*a2**3*d2*dx**2*dy**2-4.4d1*a1**2*a2**3*dx**2*dy**2 &
-4.4d1*a1**3*a2**2*dx**2*dy**2-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+4.4d1*a1**2*a2**3*d2*dy**2+4.4d1*a1**3*a2**2*d2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+1.8d1*a1*a2**3*dx**2+3.6d1*a1**2*a2**2*dx**2+1.8d1*a1**3*a2*dx**2 &
+4.0d0*a1**2*a2**3*d2**2+4.0d0*a1**3*a2**2*d2**2-2.0d1*a1*a2**3*d2 &
-4.0d1*a1**2*a2**2*d2-2.0d1*a1**3*a2*d2+7.0d0*a2**3+2.1d1*a1*a2**2 &
+2.1d1*a1**2*a2+7.0d0*a1**3)
                case (0)
                  rlYlm_laplacian = &
-3.83747515479933183d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dy* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.3d1*a1*a2**2*dx**2-3.3d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.4d1*a1*a2**2*d2+2.4d1*a1**2*a2*d2 &
-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)*dz
                case (1)
                  rlYlm_laplacian = &
-6.6467019408956851d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dx*dy* &
(4.0d0*a1**2*a2**2*d2*dy**2-2.2d1*a1*a2**2*dy**2-2.2d1*a1**2*a2*dy**2 &
+4.0d0*a1**2*a2**2*d2*dx**2-2.2d1*a1*a2**2*dx**2-2.2d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.4d1*a1*a2**2*d2+2.4d1*a1**2*a2*d2 &
-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)
                case (2)
                  rlYlm_laplacian = &
-6.6467019408956851d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dy* &
(2.0d0*a1**2*a2**2*d2*dy**2-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2 &
-2.0d0*a1**2*a2**2*d2*dx**2+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=2, m1=-1, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-3.58962931487909853d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)* &
(4.0d0*a1**2*a2**2*d2*dy**4-2.6d1*a1*a2**2*dy**4-2.6d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+7.8d1*a1*a2**2*dx**2*dy**2 &
+7.8d1*a1**2*a2*dx**2*dy**2-6.0d0*a1*a2**2*d2*dy**2 &
-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2+6.6d1*a1*a2*dy**2 &
+3.3d1*a1**2*dy**2+6.0d0*a1*a2**2*d2*dx**2+6.0d0*a1**2*a2*d2*dx**2 &
-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2-3.3d1*a1**2*dx**2)*dz
                case (-2)
                  rlYlm_laplacian = &
-8.79276018719015889d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(8.0d0*a1**3*a2**3*d2*dy**4-5.2d1*a1**2*a2**3*dy**4 &
-5.2d1*a1**3*a2**2*dy**4+8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
-5.2d1*a1**2*a2**3*dx**2*dy**2-5.2d1*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+5.2d1*a1**2*a2**3*d2*dy**2 &
+5.2d1*a1**3*a2**2*d2*dy**2-4.0d0*a1**2*a2**3*d2*dx**2 &
-4.0d0*a1**3*a2**2*d2*dx**2+2.2d1*a1*a2**3*dx**2 &
+4.4d1*a1**2*a2**2*dx**2+2.2d1*a1**3*a2*dx**2+4.0d0*a1**2*a2**3*d2**2 &
+4.0d0*a1**3*a2**2*d2**2-2.4d1*a1*a2**3*d2-4.8d1*a1**2*a2**2*d2 &
-2.4d1*a1**3*a2*d2+9.0d0*a2**3+2.7d1*a1*a2**2+2.7d1*a1**2*a2 &
+9.0d0*a1**3)
                case (-1)
                  rlYlm_laplacian = &
-2.78051491111693767d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-9)* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.3d2*a1**2*a2**3*dy**4 &
-1.3d2*a1**3*a2**2*dy**4+2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.3d2*a1**2*a2**3*dx**2*dy**2-1.3d2*a1**3*a2**2*dx**2*dy**2 &
-1.6d1*a1**3*a2**3*d2**2*dy**2+1.06d2*a1**2*a2**3*d2*dy**2 &
+1.06d2*a1**3*a2**2*d2*dy**2-1.1d1*a1*a2**3*dy**2 &
-2.2d1*a1**2*a2**2*dy**2-1.1d1*a1**3*a2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dx**2-1.0d1*a1**3*a2**2*d2*dx**2 &
+5.5d1*a1*a2**3*dx**2+1.1d2*a1**2*a2**2*dx**2+5.5d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.2d1*a1*a2**3*d2 &
-1.04d2*a1**2*a2**2*d2-5.2d1*a1**3*a2*d2+3.6d1*a2**3+1.08d2*a1*a2**2 &
+1.08d2*a1**2*a2+3.6d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian =2.27028091814553966d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(2.0d1*a1**3*a2**3*d2*dy**4 &
-1.3d2*a1**2*a2**3*dy**4-1.3d2*a1**3*a2**2*dy**4 &
+4.0d1*a1**3*a2**3*d2*dx**2*dy**2-2.6d2*a1**2*a2**3*dx**2*dy**2 &
-2.6d2*a1**3*a2**2*dx**2*dy**2-2.8d1*a1**3*a2**3*d2**2*dy**2 &
+1.88d2*a1**2*a2**3*d2*dy**2+1.88d2*a1**3*a2**2*d2*dy**2 &
-3.3d1*a1*a2**3*dy**2-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.3d2*a1**2*a2**3*dx**4 &
-1.3d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+1.88d2*a1**2*a2**3*d2*dx**2+1.88d2*a1**3*a2**2*d2*dx**2 &
-3.3d1*a1*a2**3*dx**2-6.6d1*a1**2*a2**2*dx**2-3.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-5.2d1*a1**2*a2**3*d2**2 &
-5.2d1*a1**3*a2**2*d2**2-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2 &
-6.0d0*a1**3*a2*d2+2.7d1*a2**3+8.1d1*a1*a2**2+8.1d1*a1**2*a2 &
+2.7d1*a1**3)
                case (1)
                  rlYlm_laplacian = &
-5.56102982223387534d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(1.0d1*a1**2*a2**2*d2*dy**2-6.5d1*a1*a2**2*dy**2-6.5d1*a1**2*a2*dy**2 &
+1.0d1*a1**2*a2**2*d2*dx**2-6.5d1*a1*a2**2*dx**2-6.5d1*a1**2*a2*dx**2 &
-8.0d0*a1**2*a2**2*d2**2+5.8d1*a1*a2**2*d2+5.8d1*a1**2*a2*d2 &
-3.3d1*a2**2-6.6d1*a1*a2-3.3d1*a1**2)*dz
                case (2)
                  rlYlm_laplacian =8.79276018719015889d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(4.0d0*a1**3*a2**3*d2*dy**4 &
-2.6d1*a1**2*a2**3*dy**4-2.6d1*a1**3*a2**2*dy**4 &
-4.0d0*a1**3*a2**3*d2**2*dy**2+2.4d1*a1**2*a2**3*d2*dy**2 &
+2.4d1*a1**3*a2**2*d2*dy**2+1.1d1*a1*a2**3*dy**2 &
+2.2d1*a1**2*a2**2*dy**2+1.1d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+2.6d1*a1**2*a2**3*dx**4 &
+2.6d1*a1**3*a2**2*dx**4+4.0d0*a1**3*a2**3*d2**2*dx**2 &
-3.2d1*a1**2*a2**3*d2*dx**2-3.2d1*a1**3*a2**2*d2*dx**2 &
+3.3d1*a1*a2**3*dx**2+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2 &
+4.0d0*a1**2*a2**3*d2**2+4.0d0*a1**3*a2**2*d2**2-2.4d1*a1*a2**3*d2 &
-4.8d1*a1**2*a2**2*d2-2.4d1*a1**3*a2*d2+9.0d0*a2**3+2.7d1*a1*a2**2 &
+2.7d1*a1**2*a2+9.0d0*a1**3)
                case (3)
                  rlYlm_laplacian = &
-7.17925862975819707d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.9d1*a1*a2**2*dy**2-3.9d1*a1**2*a2*dy**2 &
-2.0d0*a1**2*a2**2*d2*dx**2+1.3d1*a1*a2**2*dx**2+1.3d1*a1**2*a2*dx**2 &
-6.0d0*a1*a2**2*d2-6.0d0*a1**2*a2*d2+3.3d1*a2**2+6.6d1*a1*a2 &
+3.3d1*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=2, m1=-1, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-1.52295073829821874d1*E*a1**4*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(4.0d0*a1**2*a2**2*d2*dy**4-3.0d1*a1*a2**2*dy**4-3.0d1*a1**2*a2*dy**4 &
-4.0d0*a1**2*a2**2*d2*dx**2*dy**2+3.0d1*a1*a2**2*dx**2*dy**2 &
+3.0d1*a1**2*a2*dx**2*dy**2-6.0d0*a1*a2**2*d2*dy**2 &
-6.0d0*a1**2*a2*d2*dy**2+3.9d1*a2**2*dy**2+7.8d1*a1*a2*dy**2 &
+3.9d1*a1**2*dy**2+2.0d0*a1*a2**2*d2*dx**2+2.0d0*a1**2*a2*d2*dx**2 &
-1.3d1*a2**2*dx**2-2.6d1*a1*a2*dx**2-1.3d1*a1**2*dx**2)*dz
                case (-3)
                  rlYlm_laplacian =5.3844439723186478d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(8.0d0*a1**3*a2**3*d2*dy**6 &
-6.0d1*a1**2*a2**3*dy**6-6.0d1*a1**3*a2**2*dy**6 &
-1.6d1*a1**3*a2**3*d2*dx**2*dy**4+1.2d2*a1**2*a2**3*dx**2*dy**4 &
+1.2d2*a1**3*a2**2*dx**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**4 &
+5.2d1*a1**2*a2**3*d2*dy**4+5.2d1*a1**3*a2**2*d2*dy**4 &
+5.2d1*a1*a2**3*dy**4+1.04d2*a1**2*a2**2*dy**4+5.2d1*a1**3*a2*dy**4 &
-2.4d1*a1**3*a2**3*d2*dx**4*dy**2+1.8d2*a1**2*a2**3*dx**4*dy**2 &
+1.8d2*a1**3*a2**2*dx**4*dy**2+2.4d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-1.92d2*a1**2*a2**3*d2*dx**2*dy**2-1.92d2*a1**3*a2**2*d2*dx**2*dy**2 &
+7.8d1*a1*a2**3*dx**2*dy**2+1.56d2*a1**2*a2**2*dx**2*dy**2 &
+7.8d1*a1**3*a2*dx**2*dy**2+1.2d1*a1**2*a2**3*d2**2*dy**2 &
+1.2d1*a1**3*a2**2*d2**2*dy**2-8.4d1*a1*a2**3*d2*dy**2 &
-1.68d2*a1**2*a2**2*d2*dy**2-8.4d1*a1**3*a2*d2*dy**2+3.3d1*a2**3*dy**2 &
+9.9d1*a1*a2**2*dy**2+9.9d1*a1**2*a2*dy**2+3.3d1*a1**3*dy**2 &
+1.2d1*a1**2*a2**3*d2*dx**4+1.2d1*a1**3*a2**2*d2*dx**4 &
-7.8d1*a1*a2**3*dx**4-1.56d2*a1**2*a2**2*dx**4-7.8d1*a1**3*a2*dx**4 &
-1.2d1*a1**2*a2**3*d2**2*dx**2-1.2d1*a1**3*a2**2*d2**2*dx**2 &
+8.4d1*a1*a2**3*d2*dx**2+1.68d2*a1**2*a2**2*d2*dx**2 &
+8.4d1*a1**3*a2*d2*dx**2-3.3d1*a2**3*dx**2-9.9d1*a1*a2**2*dx**2 &
-9.9d1*a1**2*a2*dx**2-3.3d1*a1**3*dx**2)
                case (-2)
                  rlYlm_laplacian = &
-5.75621273219899775d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(2.8d1*a1**3*a2**3*d2*dy**4-2.1d2*a1**2*a2**3*dy**4 &
-2.1d2*a1**3*a2**2*dy**4+2.8d1*a1**3*a2**3*d2*dx**2*dy**2 &
-2.1d2*a1**2*a2**3*dx**2*dy**2-2.1d2*a1**3*a2**2*dx**2*dy**2 &
-2.4d1*a1**3*a2**3*d2**2*dy**2+1.86d2*a1**2*a2**3*d2*dy**2 &
+1.86d2*a1**3*a2**2*d2*dy**2-3.9d1*a1*a2**3*dy**2 &
-7.8d1*a1**2*a2**2*dy**2-3.9d1*a1**3*a2*dy**2 &
-1.4d1*a1**2*a2**3*d2*dx**2-1.4d1*a1**3*a2**2*d2*dx**2 &
+9.1d1*a1*a2**3*dx**2+1.82d2*a1**2*a2**2*dx**2+9.1d1*a1**3*a2*dx**2 &
+1.2d1*a1**2*a2**3*d2**2+1.2d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2 &
-1.8d2*a1**2*a2**2*d2-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =2.03512852844512779d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(5.6d1*a1**3*a2**3*d2*dy**6 &
-4.2d2*a1**2*a2**3*dy**6-4.2d2*a1**3*a2**2*dy**6 &
+1.12d2*a1**3*a2**3*d2*dx**2*dy**4-8.4d2*a1**2*a2**3*dx**2*dy**4 &
-8.4d2*a1**3*a2**2*dx**2*dy**4-8.8d1*a1**3*a2**3*d2**2*dy**4 &
+6.68d2*a1**2*a2**3*d2*dy**4+6.68d2*a1**3*a2**2*d2*dy**4 &
-5.2d1*a1*a2**3*dy**4-1.04d2*a1**2*a2**2*dy**4-5.2d1*a1**3*a2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**4*dy**2-4.2d2*a1**2*a2**3*dx**4*dy**2 &
-4.2d2*a1**3*a2**2*dx**4*dy**2-8.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+6.4d2*a1**2*a2**3*d2*dx**2*dy**2+6.4d2*a1**3*a2**2*d2*dx**2*dy**2 &
+1.3d2*a1*a2**3*dx**2*dy**2+2.6d2*a1**2*a2**2*dx**2*dy**2 &
+1.3d2*a1**3*a2*dx**2*dy**2+3.2d1*a1**3*a2**3*d2**3*dy**2 &
-2.2d2*a1**2*a2**3*d2**2*dy**2-2.2d2*a1**3*a2**2*d2**2*dy**2 &
-1.72d2*a1*a2**3*d2*dy**2-3.44d2*a1**2*a2**2*d2*dy**2 &
-1.72d2*a1**3*a2*d2*dy**2+2.31d2*a2**3*dy**2+6.93d2*a1*a2**2*dy**2 &
+6.93d2*a1**2*a2*dy**2+2.31d2*a1**3*dy**2-2.8d1*a1**2*a2**3*d2*dx**4 &
-2.8d1*a1**3*a2**2*d2*dx**4+1.82d2*a1*a2**3*dx**4 &
+3.64d2*a1**2*a2**2*dx**4+1.82d2*a1**3*a2*dx**4 &
+4.4d1*a1**2*a2**3*d2**2*dx**2+4.4d1*a1**3*a2**2*d2**2*dx**2 &
-3.16d2*a1*a2**3*d2*dx**2-6.32d2*a1**2*a2**2*d2*dx**2 &
-3.16d2*a1**3*a2*d2*dx**2+1.65d2*a2**3*dx**2+4.95d2*a1*a2**2*dx**2 &
+4.95d2*a1**2*a2*dx**2+1.65d2*a1**3*dx**2-1.6d1*a1**2*a2**3*d2**3 &
-1.6d1*a1**3*a2**2*d2**3+1.28d2*a1*a2**3*d2**2 &
+2.56d2*a1**2*a2**2*d2**2+1.28d2*a1**3*a2*d2**2-1.32d2*a2**3*d2 &
-3.96d2*a1*a2**2*d2-3.96d2*a1**2*a2*d2-1.32d2*a1**3*d2)
                case (0)
                  rlYlm_laplacian =1.28712829621467515d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.25d2*a1**2*a2**3*dy**4-5.25d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.05d3*a1**2*a2**3*dx**2*dy**2 &
-1.05d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.2d2*a1**2*a2**3*d2*dy**2+6.2d2*a1**3*a2**2*d2*dy**2 &
-1.3d2*a1*a2**3*dy**2-2.6d2*a1**2*a2**2*dy**2-1.3d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.25d2*a1**2*a2**3*dx**4 &
-5.25d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.2d2*a1**2*a2**3*d2*dx**2+6.2d2*a1**3*a2**2*d2*dx**2 &
-1.3d2*a1*a2**3*dx**2-2.6d2*a1**2*a2**2*dx**2-1.3d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.04d2*a1**2*a2**3*d2**2 &
-1.04d2*a1**3*a2**2*d2**2-1.52d2*a1*a2**3*d2-3.04d2*a1**2*a2**2*d2 &
-1.52d2*a1**3*a2*d2+2.64d2*a2**3+7.92d2*a1*a2**2+7.92d2*a1**2*a2 &
+2.64d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian =4.07025705689025559d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(2.8d1*a1**3*a2**3*d2*dy**4 &
-2.1d2*a1**2*a2**3*dy**4-2.1d2*a1**3*a2**2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**2*dy**2-4.2d2*a1**2*a2**3*dx**2*dy**2 &
-4.2d2*a1**3*a2**2*dx**2*dy**2-4.4d1*a1**3*a2**3*d2**2*dy**2 &
+3.48d2*a1**2*a2**3*d2*dy**2+3.48d2*a1**3*a2**2*d2*dy**2 &
-1.17d2*a1*a2**3*dy**2-2.34d2*a1**2*a2**2*dy**2-1.17d2*a1**3*a2*dy**2 &
+2.8d1*a1**3*a2**3*d2*dx**4-2.1d2*a1**2*a2**3*dx**4 &
-2.1d2*a1**3*a2**2*dx**4-4.4d1*a1**3*a2**3*d2**2*dx**2 &
+3.48d2*a1**2*a2**3*d2*dx**2+3.48d2*a1**3*a2**2*d2*dx**2 &
-1.17d2*a1*a2**3*dx**2-2.34d2*a1**2*a2**2*dx**2-1.17d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.32d2*a1**2*a2**3*d2**2 &
-1.32d2*a1**3*a2**2*d2**2+7.2d1*a1*a2**3*d2+1.44d2*a1**2*a2**2*d2 &
+7.2d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2+9.9d1*a1**2*a2 &
+3.3d1*a1**3)
                case (2)
                  rlYlm_laplacian =5.75621273219899775d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.05d2*a1**2*a2**3*dy**4-1.05d2*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2**2*dy**2+8.6d1*a1**2*a2**3*d2*dy**2 &
+8.6d1*a1**3*a2**2*d2*dy**2+2.6d1*a1*a2**3*dy**2 &
+5.2d1*a1**2*a2**2*dy**2+2.6d1*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.05d2*a1**2*a2**3*dx**4 &
+1.05d2*a1**3*a2**2*dx**4+1.2d1*a1**3*a2**3*d2**2*dx**2 &
-1.14d2*a1**2*a2**3*d2*dx**2-1.14d2*a1**3*a2**2*d2*dx**2 &
+1.56d2*a1*a2**3*dx**2+3.12d2*a1**2*a2**2*dx**2+1.56d2*a1**3*a2*dx**2 &
+1.2d1*a1**2*a2**3*d2**2+1.2d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2 &
-1.8d2*a1**2*a2**2*d2-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)*dz
                case (3)
                  rlYlm_laplacian =1.07688879446372956d1*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(1.2d1*a1**3*a2**3*d2*dy**4 &
-9.0d1*a1**2*a2**3*dy**4-9.0d1*a1**3*a2**2*dy**4 &
+8.0d0*a1**3*a2**3*d2*dx**2*dy**2-6.0d1*a1**2*a2**3*dx**2*dy**2 &
-6.0d1*a1**3*a2**2*dx**2*dy**2-1.2d1*a1**3*a2**3*d2**2*dy**2 &
+8.4d1*a1**2*a2**3*d2*dy**2+8.4d1*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+3.0d1*a1**2*a2**3*dx**4 &
+3.0d1*a1**3*a2**2*dx**4+4.0d0*a1**3*a2**3*d2**2*dx**2 &
-4.4d1*a1**2*a2**3*d2*dx**2-4.4d1*a1**3*a2**2*d2*dx**2 &
+9.1d1*a1*a2**3*dx**2+1.82d2*a1**2*a2**2*dx**2+9.1d1*a1**3*a2*dx**2 &
+1.2d1*a1**2*a2**3*d2**2+1.2d1*a1**3*a2**2*d2**2-8.4d1*a1*a2**3*d2 &
-1.68d2*a1**2*a2**2*d2-8.4d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)
                case (4)
                  rlYlm_laplacian =7.6147536914910937d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(2.0d0*a1**2*a2**2*d2*dy**4 &
-1.5d1*a1*a2**2*dy**4-1.5d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+9.0d1*a1*a2**2*dx**2*dy**2 &
+9.0d1*a1**2*a2*dx**2*dy**2-4.0d0*a1*a2**2*d2*dy**2 &
-4.0d0*a1**2*a2*d2*dy**2+2.6d1*a2**2*dy**2+5.2d1*a1*a2*dy**2 &
+2.6d1*a1**2*dy**2+2.0d0*a1**2*a2**2*d2*dx**4-1.5d1*a1*a2**2*dx**4 &
-1.5d1*a1**2*a2*dx**4+1.2d1*a1*a2**2*d2*dx**2+1.2d1*a1**2*a2*d2*dx**2 &
-7.8d1*a2**2*dx**2-1.56d2*a1*a2*dx**2-7.8d1*a1**2*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (0)
          ! selection on l2: l1=2, m1=0
          select case (l2)
            case (0)
              ! selection on m2: l1=2, m1=0, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =-9.90831824401502754d-1*E*a1*a2**3*sqrt &
                  (a2+a1)*(a2+a1)**(-6)*(1.0d0*a2*(2.0d0*a1*d2-7.0d0) &
-7.0d0*a1)*(3.0d0*(dy**2+dx**2)-2.0d0*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=2, m1=0, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-1.71617106161956686d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)*dy* &
(6.0d0*a1**2*a2**2*d2*dy**2-2.7d1*a1*a2**2*dy**2-2.7d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-2.7d1*a1*a2**2*dx**2-2.7d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.6d1*a1*a2**2*d2+1.6d1*a1**2*a2*d2 &
+7.0d0*a2**2+1.4d1*a1*a2+7.0d0*a1**2)
                case (0)
                  rlYlm_laplacian = &
-1.71617106161956686d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)* &
(6.0d0*a1**2*a2**2*d2*dy**2-2.7d1*a1*a2**2*dy**2-2.7d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-2.7d1*a1*a2**2*dx**2-2.7d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.2d1*a1*a2**2*d2+2.2d1*a1**2*a2*d2 &
-1.4d1*a2**2-2.8d1*a1*a2-1.4d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian = &
-1.71617106161956686d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)*dx* &
(6.0d0*a1**2*a2**2*d2*dy**2-2.7d1*a1*a2**2*dy**2-2.7d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-2.7d1*a1*a2**2*dx**2-2.7d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.6d1*a1*a2**2*d2+1.6d1*a1**2*a2*d2 &
+7.0d0*a2**2+1.4d1*a1*a2+7.0d0*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=2, m1=0, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-3.83747515479933183d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dx*dy* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.3d1*a1*a2**2*dx**2-3.3d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.8d1*a1*a2**2*d2+1.8d1*a1**2*a2*d2 &
+1.8d1*a2**2+3.6d1*a1*a2+1.8d1*a1**2)
                case (-1)
                  rlYlm_laplacian = &
-3.83747515479933183d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dy* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.3d1*a1*a2**2*dx**2-3.3d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.4d1*a1*a2**2*d2+2.4d1*a1**2*a2*d2 &
-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)*dz
                case (0)
                  rlYlm_laplacian =1.10778365681594752d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.8d1*a1**3*a2**3*d2*dy**4 &
-9.9d1*a1**2*a2**3*dy**4-9.9d1*a1**3*a2**2*dy**4 &
+3.6d1*a1**3*a2**3*d2*dx**2*dy**2-1.98d2*a1**2*a2**3*dx**2*dy**2 &
-1.98d2*a1**3*a2**2*dx**2*dy**2-2.4d1*a1**3*a2**3*d2**2*dy**2 &
+1.44d2*a1**2*a2**3*d2*dy**2+1.44d2*a1**3*a2**2*d2*dy**2 &
-5.4d1*a1*a2**3*dy**2-1.08d2*a1**2*a2**2*dy**2-5.4d1*a1**3*a2*dy**2 &
+1.8d1*a1**3*a2**3*d2*dx**4-9.9d1*a1**2*a2**3*dx**4 &
-9.9d1*a1**3*a2**2*dx**4-2.4d1*a1**3*a2**3*d2**2*dx**2 &
+1.44d2*a1**2*a2**3*d2*dx**2+1.44d2*a1**3*a2**2*d2*dx**2 &
-5.4d1*a1*a2**3*dx**2-1.08d2*a1**2*a2**2*dx**2-5.4d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-6.0d1*a1**2*a2**3*d2**2 &
-6.0d1*a1**3*a2**2*d2**2+7.8d1*a1*a2**3*d2+1.56d2*a1**2*a2**2*d2 &
+7.8d1*a1**3*a2*d2-2.1d1*a2**3-6.3d1*a1*a2**2-6.3d1*a1**2*a2 &
-2.1d1*a1**3)
                case (1)
                  rlYlm_laplacian = &
-3.83747515479933183d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dx* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.3d1*a1*a2**2*dx**2-3.3d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.4d1*a1*a2**2*d2+2.4d1*a1**2*a2*d2 &
-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)*dz
                case (2)
                  rlYlm_laplacian =1.91873757739966592d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(dy+dx)*(dy-1.0d0*dx)* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.3d1*a1*a2**2*dx**2-3.3d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.8d1*a1*a2**2*d2+1.8d1*a1**2*a2*d2 &
+1.8d1*a2**2+3.6d1*a1*a2+1.8d1*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=2, m1=0, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =2.07247345123641944d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(dy**2-3.0d0*dx**2)* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.9d1*a1*a2**2*dy**2-3.9d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.9d1*a1*a2**2*dx**2-3.9d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.0d1*a1*a2**2*d2+2.0d1*a1**2*a2*d2 &
+3.3d1*a2**2+6.6d1*a1*a2+3.3d1*a1**2)
                case (-2)
                  rlYlm_laplacian = &
-1.01530049219881249d1*E*a1**4*a2**3*sqrt(a2+a1)*(a2+a1)**(-9)* &
(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*dx*dy*(3.0d0*(dy**2+dx**2) &
-2.0d0*d2)*dz
                case (-1)
                  rlYlm_laplacian =1.60533103241913232d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(3.0d1*a1**3*a2**3*d2*dy**4 &
-1.95d2*a1**2*a2**3*dy**4-1.95d2*a1**3*a2**2*dy**4 &
+6.0d1*a1**3*a2**3*d2*dx**2*dy**2-3.9d2*a1**2*a2**3*dx**2*dy**2 &
-3.9d2*a1**3*a2**2*dx**2*dy**2-4.4d1*a1**3*a2**3*d2**2*dy**2 &
+3.04d2*a1**2*a2**3*d2*dy**2+3.04d2*a1**3*a2**2*d2*dy**2 &
-9.9d1*a1*a2**3*dy**2-1.98d2*a1**2*a2**2*dy**2-9.9d1*a1**3*a2*dy**2 &
+3.0d1*a1**3*a2**3*d2*dx**4-1.95d2*a1**2*a2**3*dx**4 &
-1.95d2*a1**3*a2**2*dx**4-4.4d1*a1**3*a2**3*d2**2*dx**2 &
+3.04d2*a1**2*a2**3*d2*dx**2+3.04d2*a1**3*a2**2*d2*dx**2 &
-9.9d1*a1*a2**3*dx**2-1.98d2*a1**2*a2**2*dx**2-9.9d1*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.28d2*a1**2*a2**3*d2**2 &
-1.28d2*a1**3*a2**2*d2**2+1.44d2*a1*a2**3*d2+2.88d2*a1**2*a2**2*d2 &
+1.44d2*a1**3*a2*d2-5.4d1*a2**3-1.62d2*a1*a2**2-1.62d2*a1**2*a2 &
-5.4d1*a1**3)
                case (0)
                  rlYlm_laplacian =1.31074729922739806d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(3.0d1*a1**3*a2**3*d2*dy**4 &
-1.95d2*a1**2*a2**3*dy**4-1.95d2*a1**3*a2**2*dy**4 &
+6.0d1*a1**3*a2**3*d2*dx**2*dy**2-3.9d2*a1**2*a2**3*dx**2*dy**2 &
-3.9d2*a1**3*a2**2*dx**2*dy**2-3.2d1*a1**3*a2**3*d2**2*dy**2 &
+2.32d2*a1**2*a2**3*d2*dy**2+2.32d2*a1**3*a2**2*d2*dy**2 &
-1.32d2*a1*a2**3*dy**2-2.64d2*a1**2*a2**2*dy**2-1.32d2*a1**3*a2*dy**2 &
+3.0d1*a1**3*a2**3*d2*dx**4-1.95d2*a1**2*a2**3*dx**4 &
-1.95d2*a1**3*a2**2*dx**4-3.2d1*a1**3*a2**3*d2**2*dx**2 &
+2.32d2*a1**2*a2**3*d2*dx**2+2.32d2*a1**3*a2**2*d2*dx**2 &
-1.32d2*a1*a2**3*dx**2-2.64d2*a1**2*a2**2*dx**2-1.32d2*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-7.6d1*a1**2*a2**3*d2**2 &
-7.6d1*a1**3*a2**2*d2**2+1.5d2*a1*a2**3*d2+3.0d2*a1**2*a2**2*d2 &
+1.5d2*a1**3*a2*d2-8.1d1*a2**3-2.43d2*a1*a2**2-2.43d2*a1**2*a2 &
-8.1d1*a1**3)*dz
                case (1)
                  rlYlm_laplacian =1.60533103241913232d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(3.0d1*a1**3*a2**3*d2*dy**4 &
-1.95d2*a1**2*a2**3*dy**4-1.95d2*a1**3*a2**2*dy**4 &
+6.0d1*a1**3*a2**3*d2*dx**2*dy**2-3.9d2*a1**2*a2**3*dx**2*dy**2 &
-3.9d2*a1**3*a2**2*dx**2*dy**2-4.4d1*a1**3*a2**3*d2**2*dy**2 &
+3.04d2*a1**2*a2**3*d2*dy**2+3.04d2*a1**3*a2**2*d2*dy**2 &
-9.9d1*a1*a2**3*dy**2-1.98d2*a1**2*a2**2*dy**2-9.9d1*a1**3*a2*dy**2 &
+3.0d1*a1**3*a2**3*d2*dx**4-1.95d2*a1**2*a2**3*dx**4 &
-1.95d2*a1**3*a2**2*dx**4-4.4d1*a1**3*a2**3*d2**2*dx**2 &
+3.04d2*a1**2*a2**3*d2*dx**2+3.04d2*a1**3*a2**2*d2*dx**2 &
-9.9d1*a1*a2**3*dx**2-1.98d2*a1**2*a2**2*dx**2-9.9d1*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.28d2*a1**2*a2**3*d2**2 &
-1.28d2*a1**3*a2**2*d2**2+1.44d2*a1*a2**3*d2+2.88d2*a1**2*a2**2*d2 &
+1.44d2*a1**3*a2*d2-5.4d1*a2**3-1.62d2*a1*a2**2-1.62d2*a1**2*a2 &
-5.4d1*a1**3)
                case (2)
                  rlYlm_laplacian =5.07650246099406246d0*E*a1**4*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-9)*(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*(dy+dx)* &
(dy-1.0d0*dx)*(3.0d0*(dy**2+dx**2)-2.0d0*d2)*dz
                case (3)
                  rlYlm_laplacian =2.07247345123641944d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(3.0d0*dy**2-1.0d0*dx**2)* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.9d1*a1*a2**2*dy**2-3.9d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.9d1*a1*a2**2*dx**2-3.9d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.0d1*a1*a2**2*d2+2.0d1*a1**2*a2*d2 &
+3.3d1*a2**2+6.6d1*a1*a2+3.3d1*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=2, m1=0, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =8.79276018719015889d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(6.0d0*a1**2*a2**2*d2*dy**2-4.5d1*a1*a2**2*dy**2-4.5d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-4.5d1*a1*a2**2*dx**2-4.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.2d1*a1*a2**2*d2+2.2d1*a1**2*a2*d2 &
+5.2d1*a2**2+1.04d2*a1*a2+5.2d1*a1**2)
                case (-3)
                  rlYlm_laplacian =6.21742035370925833d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(dy**2-3.0d0*dx**2)* &
(6.0d0*a1**2*a2**2*d2*dy**2-4.5d1*a1*a2**2*dy**2-4.5d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-4.5d1*a1*a2**2*dx**2-4.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.8d1*a1*a2**2*d2+2.8d1*a1**2*a2*d2 &
+1.3d1*a2**2+2.6d1*a1*a2+1.3d1*a1**2)*dz
                case (-2)
                  rlYlm_laplacian =3.32335097044784255d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(4.2d1*a1**3*a2**3*d2*dy**4 &
-3.15d2*a1**2*a2**3*dy**4-3.15d2*a1**3*a2**2*dy**4 &
+8.4d1*a1**3*a2**3*d2*dx**2*dy**2-6.3d2*a1**2*a2**3*dx**2*dy**2 &
-6.3d2*a1**3*a2**2*dx**2*dy**2-6.4d1*a1**3*a2**3*d2**2*dy**2 &
+4.96d2*a1**2*a2**3*d2*dy**2+4.96d2*a1**3*a2**2*d2*dy**2 &
-1.04d2*a1*a2**3*dy**2-2.08d2*a1**2*a2**2*dy**2-1.04d2*a1**3*a2*dy**2 &
+4.2d1*a1**3*a2**3*d2*dx**4-3.15d2*a1**2*a2**3*dx**4 &
-3.15d2*a1**3*a2**2*dx**4-6.4d1*a1**3*a2**3*d2**2*dx**2 &
+4.96d2*a1**2*a2**3*d2*dx**2+4.96d2*a1**3*a2**2*d2*dx**2 &
-1.04d2*a1*a2**3*dx**2-2.08d2*a1**2*a2**2*dx**2-1.04d2*a1**3*a2*dx**2 &
+2.4d1*a1**3*a2**3*d2**3-2.04d2*a1**2*a2**3*d2**2 &
-2.04d2*a1**3*a2**2*d2**2+1.74d2*a1*a2**3*d2+3.48d2*a1**2*a2**2*d2 &
+1.74d2*a1**3*a2*d2-9.9d1*a2**3-2.97d2*a1*a2**2-2.97d2*a1**2*a2 &
-9.9d1*a1**3)
                case (-1)
                  rlYlm_laplacian =2.34996400746656297d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(4.2d1*a1**3*a2**3*d2*dy**4 &
-3.15d2*a1**2*a2**3*dy**4-3.15d2*a1**3*a2**2*dy**4 &
+8.4d1*a1**3*a2**3*d2*dx**2*dy**2-6.3d2*a1**2*a2**3*dx**2*dy**2 &
-6.3d2*a1**3*a2**2*dx**2*dy**2-5.2d1*a1**3*a2**3*d2**2*dy**2 &
+4.24d2*a1**2*a2**3*d2*dy**2+4.24d2*a1**3*a2**2*d2*dy**2 &
-2.21d2*a1*a2**3*dy**2-4.42d2*a1**2*a2**2*dy**2-2.21d2*a1**3*a2*dy**2 &
+4.2d1*a1**3*a2**3*d2*dx**4-3.15d2*a1**2*a2**3*dx**4 &
-3.15d2*a1**3*a2**2*dx**4-5.2d1*a1**3*a2**3*d2**2*dx**2 &
+4.24d2*a1**2*a2**3*d2*dx**2+4.24d2*a1**3*a2**2*d2*dx**2 &
-2.21d2*a1*a2**3*dx**2-4.42d2*a1**2*a2**2*dx**2-2.21d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.6d2*a1**2*a2**3*d2**2 &
-1.6d2*a1**3*a2**2*d2**2+2.96d2*a1*a2**3*d2+5.92d2*a1**2*a2**2*d2 &
+2.96d2*a1**3*a2*d2-1.98d2*a2**3-5.94d2*a1*a2**2-5.94d2*a1**2*a2 &
-1.98d2*a1**3)*dz
                case (0)
                  rlYlm_laplacian =-3.71561934150563533d-1*E*a1**3*a2*sqrt &
                  (a2+a1)*(a2+a1)**(-10)*(2.1d2*a1**3*a2**3*d2*dy**6 &
-1.575d3*a1**2*a2**3*dy**6-1.575d3*a1**3*a2**2*dy**6 &
+6.3d2*a1**3*a2**3*d2*dx**2*dy**4-4.725d3*a1**2*a2**3*dx**2*dy**4 &
-4.725d3*a1**3*a2**2*dx**2*dy**4-3.8d2*a1**3*a2**3*d2**2*dy**4 &
+3.05d3*a1**2*a2**3*d2*dy**4+3.05d3*a1**3*a2**2*d2*dy**4 &
-1.3d3*a1*a2**3*dy**4-2.6d3*a1**2*a2**2*dy**4-1.3d3*a1**3*a2*dy**4 &
+6.3d2*a1**3*a2**3*d2*dx**4*dy**2-4.725d3*a1**2*a2**3*dx**4*dy**2 &
-4.725d3*a1**3*a2**2*dx**4*dy**2-7.6d2*a1**3*a2**3*d2**2*dx**2*dy**2 &
+6.1d3*a1**2*a2**3*d2*dx**2*dy**2+6.1d3*a1**3*a2**2*d2*dx**2*dy**2 &
-2.6d3*a1*a2**3*dx**2*dy**2-5.2d3*a1**2*a2**2*dx**2*dy**2 &
-2.6d3*a1**3*a2*dx**2*dy**2+2.08d2*a1**3*a2**3*d2**3*dy**2 &
-1.912d3*a1**2*a2**3*d2**2*dy**2-1.912d3*a1**3*a2**2*d2**2*dy**2 &
+2.504d3*a1*a2**3*d2*dy**2+5.008d3*a1**2*a2**2*d2*dy**2 &
+2.504d3*a1**3*a2*d2*dy**2-1.188d3*a2**3*dy**2-3.564d3*a1*a2**2*dy**2 &
-3.564d3*a1**2*a2*dy**2-1.188d3*a1**3*dy**2+2.1d2*a1**3*a2**3*d2*dx**6 &
-1.575d3*a1**2*a2**3*dx**6-1.575d3*a1**3*a2**2*dx**6 &
-3.8d2*a1**3*a2**3*d2**2*dx**4+3.05d3*a1**2*a2**3*d2*dx**4 &
+3.05d3*a1**3*a2**2*d2*dx**4-1.3d3*a1*a2**3*dx**4 &
-2.6d3*a1**2*a2**2*dx**4-1.3d3*a1**3*a2*dx**4 &
+2.08d2*a1**3*a2**3*d2**3*dx**2-1.912d3*a1**2*a2**3*d2**2*dx**2 &
-1.912d3*a1**3*a2**2*d2**2*dx**2+2.504d3*a1*a2**3*d2*dx**2 &
+5.008d3*a1**2*a2**2*d2*dx**2+2.504d3*a1**3*a2*d2*dx**2 &
-1.188d3*a2**3*dx**2-3.564d3*a1*a2**2*dx**2-3.564d3*a1**2*a2*dx**2 &
-1.188d3*a1**3*dx**2-3.2d1*a1**3*a2**3*d2**4+3.68d2*a1**2*a2**3*d2**3 &
+3.68d2*a1**3*a2**2*d2**3-9.76d2*a1*a2**3*d2**2 &
-1.952d3*a1**2*a2**2*d2**2-9.76d2*a1**3*a2*d2**2+7.92d2*a2**3*d2 &
+2.376d3*a1*a2**2*d2+2.376d3*a1**2*a2*d2+7.92d2*a1**3*d2)
                case (1)
                  rlYlm_laplacian =2.34996400746656297d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(4.2d1*a1**3*a2**3*d2*dy**4 &
-3.15d2*a1**2*a2**3*dy**4-3.15d2*a1**3*a2**2*dy**4 &
+8.4d1*a1**3*a2**3*d2*dx**2*dy**2-6.3d2*a1**2*a2**3*dx**2*dy**2 &
-6.3d2*a1**3*a2**2*dx**2*dy**2-5.2d1*a1**3*a2**3*d2**2*dy**2 &
+4.24d2*a1**2*a2**3*d2*dy**2+4.24d2*a1**3*a2**2*d2*dy**2 &
-2.21d2*a1*a2**3*dy**2-4.42d2*a1**2*a2**2*dy**2-2.21d2*a1**3*a2*dy**2 &
+4.2d1*a1**3*a2**3*d2*dx**4-3.15d2*a1**2*a2**3*dx**4 &
-3.15d2*a1**3*a2**2*dx**4-5.2d1*a1**3*a2**3*d2**2*dx**2 &
+4.24d2*a1**2*a2**3*d2*dx**2+4.24d2*a1**3*a2**2*d2*dx**2 &
-2.21d2*a1*a2**3*dx**2-4.42d2*a1**2*a2**2*dx**2-2.21d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.6d2*a1**2*a2**3*d2**2 &
-1.6d2*a1**3*a2**2*d2**2+2.96d2*a1*a2**3*d2+5.92d2*a1**2*a2**2*d2 &
+2.96d2*a1**3*a2*d2-1.98d2*a2**3-5.94d2*a1*a2**2-5.94d2*a1**2*a2 &
-1.98d2*a1**3)*dz
                case (2)
                  rlYlm_laplacian = &
-1.66167548522392128d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*(dy+dx)* &
(dy-1.0d0*dx)*(4.2d1*a1**3*a2**3*d2*dy**4-3.15d2*a1**2*a2**3*dy**4 &
-3.15d2*a1**3*a2**2*dy**4+8.4d1*a1**3*a2**3*d2*dx**2*dy**2 &
-6.3d2*a1**2*a2**3*dx**2*dy**2-6.3d2*a1**3*a2**2*dx**2*dy**2 &
-6.4d1*a1**3*a2**3*d2**2*dy**2+4.96d2*a1**2*a2**3*d2*dy**2 &
+4.96d2*a1**3*a2**2*d2*dy**2-1.04d2*a1*a2**3*dy**2 &
-2.08d2*a1**2*a2**2*dy**2-1.04d2*a1**3*a2*dy**2 &
+4.2d1*a1**3*a2**3*d2*dx**4-3.15d2*a1**2*a2**3*dx**4 &
-3.15d2*a1**3*a2**2*dx**4-6.4d1*a1**3*a2**3*d2**2*dx**2 &
+4.96d2*a1**2*a2**3*d2*dx**2+4.96d2*a1**3*a2**2*d2*dx**2 &
-1.04d2*a1*a2**3*dx**2-2.08d2*a1**2*a2**2*dx**2-1.04d2*a1**3*a2*dx**2 &
+2.4d1*a1**3*a2**3*d2**3-2.04d2*a1**2*a2**3*d2**2 &
-2.04d2*a1**3*a2**2*d2**2+1.74d2*a1*a2**3*d2+3.48d2*a1**2*a2**2*d2 &
+1.74d2*a1**3*a2*d2-9.9d1*a2**3-2.97d2*a1*a2**2-2.97d2*a1**2*a2 &
-9.9d1*a1**3)
                case (3)
                  rlYlm_laplacian =6.21742035370925833d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(3.0d0*dy**2-1.0d0*dx**2)* &
(6.0d0*a1**2*a2**2*d2*dy**2-4.5d1*a1*a2**2*dy**2-4.5d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-4.5d1*a1*a2**2*dx**2-4.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.8d1*a1*a2**2*d2+2.8d1*a1**2*a2*d2 &
+1.3d1*a2**2+2.6d1*a1*a2+1.3d1*a1**2)*dz
                case (4)
                  rlYlm_laplacian = &
-2.19819004679753972d0*E*a1**4*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*(dy**2 &
-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)* &
(6.0d0*a1**2*a2**2*d2*dy**2-4.5d1*a1*a2**2*dy**2-4.5d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-4.5d1*a1*a2**2*dx**2-4.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.2d1*a1*a2**2*d2+2.2d1*a1**2*a2*d2 &
+5.2d1*a2**2+1.04d2*a1*a2+5.2d1*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (1)
          ! selection on l2: l1=2, m1=1
          select case (l2)
            case (0)
              ! selection on m2: l1=2, m1=1, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =3.43234212323913373d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2*(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*dx*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=2, m1=1, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =5.94499094640901652d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*dy*dz
                case (0)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)*dx* &
(4.0d0*a1**2*a2**2*d2*dy**2-1.8d1*a1*a2**2*dy**2-1.8d1*a1**2*a2*dy**2 &
+4.0d0*a1**2*a2**2*d2*dx**2-1.8d1*a1*a2**2*dx**2-1.8d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.0d1*a1*a2**2*d2+2.0d1*a1**2*a2*d2 &
-7.0d0*a2**2-1.4d1*a1*a2-7.0d0*a1**2)
                case (1)
                  rlYlm_laplacian =2.97249547320450826d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-7)*(4.0d0*a1**2*a2**2*d2*dx**2-1.8d1*a1*a2**2*dx**2 &
-1.8d1*a1**2*a2*dx**2-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+7.0d0*a2**2 &
+1.4d1*a1*a2+7.0d0*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=2, m1=1, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =6.6467019408956851d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(4.0d0*a1**2*a2**2*d2*dx**2-2.2d1*a1*a2**2*dx**2 &
-2.2d1*a1**2*a2*dx**2-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+9.0d0*a2**2 &
+1.8d1*a1*a2+9.0d0*a1**2)*dy*dz
                case (-1)
                  rlYlm_laplacian = &
-6.6467019408956851d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dx*dy* &
(4.0d0*a1**2*a2**2*d2*dy**2-2.2d1*a1*a2**2*dy**2-2.2d1*a1**2*a2*dy**2 &
+4.0d0*a1**2*a2**2*d2*dx**2-2.2d1*a1*a2**2*dx**2-2.2d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.4d1*a1*a2**2*d2+2.4d1*a1**2*a2*d2 &
-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)
                case (0)
                  rlYlm_laplacian = &
-3.83747515479933183d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dx* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.3d1*a1*a2**2*dx**2-3.3d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.4d1*a1*a2**2*d2+2.4d1*a1**2*a2*d2 &
-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)*dz
                case (1)
                  rlYlm_laplacian =-3.32335097044784255d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
-4.4d1*a1**2*a2**3*dx**2*dy**2-4.4d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+1.8d1*a1*a2**3*dy**2+3.6d1*a1**2*a2**2*dy**2+1.8d1*a1**3*a2*dy**2 &
+8.0d0*a1**3*a2**3*d2*dx**4-4.4d1*a1**2*a2**3*dx**4 &
-4.4d1*a1**3*a2**2*dx**4-8.0d0*a1**3*a2**3*d2**2*dx**2 &
+4.4d1*a1**2*a2**3*d2*dx**2+4.4d1*a1**3*a2**2*d2*dx**2 &
+4.0d0*a1**2*a2**3*d2**2+4.0d0*a1**3*a2**2*d2**2-2.0d1*a1*a2**3*d2 &
-4.0d1*a1**2*a2**2*d2-2.0d1*a1**3*a2*d2+7.0d0*a2**3+2.1d1*a1*a2**2 &
+2.1d1*a1**2*a2+7.0d0*a1**3)
                case (2)
                  rlYlm_laplacian = &
-6.6467019408956851d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dx* &
(2.0d0*a1**2*a2**2*d2*dy**2-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2 &
-2.0d0*a1**2*a2**2*d2*dx**2+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2 &
+2.0d0*a1*a2**2*d2+2.0d0*a1**2*a2*d2-9.0d0*a2**2-1.8d1*a1*a2 &
-9.0d0*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=2, m1=1, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-7.17925862975819707d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(2.0d0*a1**2*a2**2*d2*dy**2-1.3d1*a1*a2**2*dy**2-1.3d1*a1**2*a2*dy**2 &
-6.0d0*a1**2*a2**2*d2*dx**2+3.9d1*a1*a2**2*dx**2+3.9d1*a1**2*a2*dx**2 &
+6.0d0*a1*a2**2*d2+6.0d0*a1**2*a2*d2-3.3d1*a2**2-6.6d1*a1*a2 &
-3.3d1*a1**2)*dz
                case (-2)
                  rlYlm_laplacian = &
-8.79276018719015889d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(8.0d0*a1**3*a2**3*d2*dx**2*dy**2-5.2d1*a1**2*a2**3*dx**2*dy**2 &
-5.2d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**2*a2**3*d2*dy**2 &
-4.0d0*a1**3*a2**2*d2*dy**2+2.2d1*a1*a2**3*dy**2 &
+4.4d1*a1**2*a2**2*dy**2+2.2d1*a1**3*a2*dy**2 &
+8.0d0*a1**3*a2**3*d2*dx**4-5.2d1*a1**2*a2**3*dx**4 &
-5.2d1*a1**3*a2**2*dx**4-8.0d0*a1**3*a2**3*d2**2*dx**2 &
+5.2d1*a1**2*a2**3*d2*dx**2+5.2d1*a1**3*a2**2*d2*dx**2 &
+4.0d0*a1**2*a2**3*d2**2+4.0d0*a1**3*a2**2*d2**2-2.4d1*a1*a2**3*d2 &
-4.8d1*a1**2*a2**2*d2-2.4d1*a1**3*a2*d2+9.0d0*a2**3+2.7d1*a1*a2**2 &
+2.7d1*a1**2*a2+9.0d0*a1**3)
                case (-1)
                  rlYlm_laplacian = &
-5.56102982223387534d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(1.0d1*a1**2*a2**2*d2*dy**2-6.5d1*a1*a2**2*dy**2-6.5d1*a1**2*a2*dy**2 &
+1.0d1*a1**2*a2**2*d2*dx**2-6.5d1*a1*a2**2*dx**2-6.5d1*a1**2*a2*dx**2 &
-8.0d0*a1**2*a2**2*d2**2+5.8d1*a1*a2**2*d2+5.8d1*a1**2*a2*d2 &
-3.3d1*a2**2-6.6d1*a1*a2-3.3d1*a1**2)*dz
                case (0)
                  rlYlm_laplacian =2.27028091814553966d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(2.0d1*a1**3*a2**3*d2*dy**4 &
-1.3d2*a1**2*a2**3*dy**4-1.3d2*a1**3*a2**2*dy**4 &
+4.0d1*a1**3*a2**3*d2*dx**2*dy**2-2.6d2*a1**2*a2**3*dx**2*dy**2 &
-2.6d2*a1**3*a2**2*dx**2*dy**2-2.8d1*a1**3*a2**3*d2**2*dy**2 &
+1.88d2*a1**2*a2**3*d2*dy**2+1.88d2*a1**3*a2**2*d2*dy**2 &
-3.3d1*a1*a2**3*dy**2-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.3d2*a1**2*a2**3*dx**4 &
-1.3d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+1.88d2*a1**2*a2**3*d2*dx**2+1.88d2*a1**3*a2**2*d2*dx**2 &
-3.3d1*a1*a2**3*dx**2-6.6d1*a1**2*a2**2*dx**2-3.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-5.2d1*a1**2*a2**3*d2**2 &
-5.2d1*a1**3*a2**2*d2**2-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2 &
-6.0d0*a1**3*a2*d2+2.7d1*a2**3+8.1d1*a1*a2**2+8.1d1*a1**2*a2 &
+2.7d1*a1**3)
                case (1)
                  rlYlm_laplacian = &
-2.78051491111693767d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-9)* &
(2.0d1*a1**3*a2**3*d2*dx**2*dy**2-1.3d2*a1**2*a2**3*dx**2*dy**2 &
-1.3d2*a1**3*a2**2*dx**2*dy**2-1.0d1*a1**2*a2**3*d2*dy**2 &
-1.0d1*a1**3*a2**2*d2*dy**2+5.5d1*a1*a2**3*dy**2 &
+1.1d2*a1**2*a2**2*dy**2+5.5d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.3d2*a1**2*a2**3*dx**4 &
-1.3d2*a1**3*a2**2*dx**4-1.6d1*a1**3*a2**3*d2**2*dx**2 &
+1.06d2*a1**2*a2**3*d2*dx**2+1.06d2*a1**3*a2**2*d2*dx**2 &
-1.1d1*a1*a2**3*dx**2-2.2d1*a1**2*a2**2*dx**2-1.1d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.2d1*a1*a2**3*d2 &
-1.04d2*a1**2*a2**2*d2-5.2d1*a1**3*a2*d2+3.6d1*a2**3+1.08d2*a1*a2**2 &
+1.08d2*a1**2*a2+3.6d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian =8.79276018719015889d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(4.0d0*a1**3*a2**3*d2*dy**4 &
-2.6d1*a1**2*a2**3*dy**4-2.6d1*a1**3*a2**2*dy**4 &
-4.0d0*a1**3*a2**3*d2**2*dy**2+3.2d1*a1**2*a2**3*d2*dy**2 &
+3.2d1*a1**3*a2**2*d2*dy**2-3.3d1*a1*a2**3*dy**2 &
-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+2.6d1*a1**2*a2**3*dx**4 &
+2.6d1*a1**3*a2**2*dx**4+4.0d0*a1**3*a2**3*d2**2*dx**2 &
-2.4d1*a1**2*a2**3*d2*dx**2-2.4d1*a1**3*a2**2*d2*dx**2 &
-1.1d1*a1*a2**3*dx**2-2.2d1*a1**2*a2**2*dx**2-1.1d1*a1**3*a2*dx**2 &
-4.0d0*a1**2*a2**3*d2**2-4.0d0*a1**3*a2**2*d2**2+2.4d1*a1*a2**3*d2 &
+4.8d1*a1**2*a2**2*d2+2.4d1*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)
                case (3)
                  rlYlm_laplacian = &
-3.58962931487909853d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)* &
(1.2d1*a1**2*a2**2*d2*dx**2*dy**2-7.8d1*a1*a2**2*dx**2*dy**2 &
-7.8d1*a1**2*a2*dx**2*dy**2-6.0d0*a1*a2**2*d2*dy**2 &
-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2+6.6d1*a1*a2*dy**2 &
+3.3d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4+2.6d1*a1*a2**2*dx**4 &
+2.6d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2+6.0d0*a1**2*a2*d2*dx**2 &
-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2-3.3d1*a1**2*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=2, m1=1, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-1.52295073829821874d1*E*a1**4*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(4.0d0*a1**2*a2**2*d2*dx**2*dy**2-3.0d1*a1*a2**2*dx**2*dy**2 &
-3.0d1*a1**2*a2*dx**2*dy**2-2.0d0*a1*a2**2*d2*dy**2 &
-2.0d0*a1**2*a2*d2*dy**2+1.3d1*a2**2*dy**2+2.6d1*a1*a2*dy**2 &
+1.3d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4+3.0d1*a1*a2**2*dx**4 &
+3.0d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2+6.0d0*a1**2*a2*d2*dx**2 &
-3.9d1*a2**2*dx**2-7.8d1*a1*a2*dx**2-3.9d1*a1**2*dx**2)*dz
                case (-3)
                  rlYlm_laplacian =1.07688879446372956d1*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(4.0d0*a1**3*a2**3*d2*dy**4 &
-3.0d1*a1**2*a2**3*dy**4-3.0d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2*dx**2*dy**2+6.0d1*a1**2*a2**3*dx**2*dy**2 &
+6.0d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**3*a2**3*d2**2*dy**2 &
+4.4d1*a1**2*a2**3*d2*dy**2+4.4d1*a1**3*a2**2*d2*dy**2 &
-9.1d1*a1*a2**3*dy**2-1.82d2*a1**2*a2**2*dy**2-9.1d1*a1**3*a2*dy**2 &
-1.2d1*a1**3*a2**3*d2*dx**4+9.0d1*a1**2*a2**3*dx**4 &
+9.0d1*a1**3*a2**2*dx**4+1.2d1*a1**3*a2**3*d2**2*dx**2 &
-8.4d1*a1**2*a2**3*d2*dx**2-8.4d1*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
-1.2d1*a1**2*a2**3*d2**2-1.2d1*a1**3*a2**2*d2**2+8.4d1*a1*a2**3*d2 &
+1.68d2*a1**2*a2**2*d2+8.4d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)
                case (-2)
                  rlYlm_laplacian = &
-5.75621273219899775d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(2.8d1*a1**3*a2**3*d2*dx**2*dy**2-2.1d2*a1**2*a2**3*dx**2*dy**2 &
-2.1d2*a1**3*a2**2*dx**2*dy**2-1.4d1*a1**2*a2**3*d2*dy**2 &
-1.4d1*a1**3*a2**2*d2*dy**2+9.1d1*a1*a2**3*dy**2 &
+1.82d2*a1**2*a2**2*dy**2+9.1d1*a1**3*a2*dy**2 &
+2.8d1*a1**3*a2**3*d2*dx**4-2.1d2*a1**2*a2**3*dx**4 &
-2.1d2*a1**3*a2**2*dx**4-2.4d1*a1**3*a2**3*d2**2*dx**2 &
+1.86d2*a1**2*a2**3*d2*dx**2+1.86d2*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
+1.2d1*a1**2*a2**3*d2**2+1.2d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2 &
-1.8d2*a1**2*a2**2*d2-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =4.07025705689025559d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(2.8d1*a1**3*a2**3*d2*dy**4 &
-2.1d2*a1**2*a2**3*dy**4-2.1d2*a1**3*a2**2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**2*dy**2-4.2d2*a1**2*a2**3*dx**2*dy**2 &
-4.2d2*a1**3*a2**2*dx**2*dy**2-4.4d1*a1**3*a2**3*d2**2*dy**2 &
+3.48d2*a1**2*a2**3*d2*dy**2+3.48d2*a1**3*a2**2*d2*dy**2 &
-1.17d2*a1*a2**3*dy**2-2.34d2*a1**2*a2**2*dy**2-1.17d2*a1**3*a2*dy**2 &
+2.8d1*a1**3*a2**3*d2*dx**4-2.1d2*a1**2*a2**3*dx**4 &
-2.1d2*a1**3*a2**2*dx**4-4.4d1*a1**3*a2**3*d2**2*dx**2 &
+3.48d2*a1**2*a2**3*d2*dx**2+3.48d2*a1**3*a2**2*d2*dx**2 &
-1.17d2*a1*a2**3*dx**2-2.34d2*a1**2*a2**2*dx**2-1.17d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.32d2*a1**2*a2**3*d2**2 &
-1.32d2*a1**3*a2**2*d2**2+7.2d1*a1*a2**3*d2+1.44d2*a1**2*a2**2*d2 &
+7.2d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2+9.9d1*a1**2*a2 &
+3.3d1*a1**3)
                case (0)
                  rlYlm_laplacian =1.28712829621467515d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.25d2*a1**2*a2**3*dy**4-5.25d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.05d3*a1**2*a2**3*dx**2*dy**2 &
-1.05d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.2d2*a1**2*a2**3*d2*dy**2+6.2d2*a1**3*a2**2*d2*dy**2 &
-1.3d2*a1*a2**3*dy**2-2.6d2*a1**2*a2**2*dy**2-1.3d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.25d2*a1**2*a2**3*dx**4 &
-5.25d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.2d2*a1**2*a2**3*d2*dx**2+6.2d2*a1**3*a2**2*d2*dx**2 &
-1.3d2*a1*a2**3*dx**2-2.6d2*a1**2*a2**2*dx**2-1.3d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.04d2*a1**2*a2**3*d2**2 &
-1.04d2*a1**3*a2**2*d2**2-1.52d2*a1*a2**3*d2-3.04d2*a1**2*a2**2*d2 &
-1.52d2*a1**3*a2*d2+2.64d2*a2**3+7.92d2*a1*a2**2+7.92d2*a1**2*a2 &
+2.64d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian =2.03512852844512779d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(5.6d1*a1**3*a2**3*d2*dx**2*dy**4 &
-4.2d2*a1**2*a2**3*dx**2*dy**4-4.2d2*a1**3*a2**2*dx**2*dy**4 &
-2.8d1*a1**2*a2**3*d2*dy**4-2.8d1*a1**3*a2**2*d2*dy**4 &
+1.82d2*a1*a2**3*dy**4+3.64d2*a1**2*a2**2*dy**4+1.82d2*a1**3*a2*dy**4 &
+1.12d2*a1**3*a2**3*d2*dx**4*dy**2-8.4d2*a1**2*a2**3*dx**4*dy**2 &
-8.4d2*a1**3*a2**2*dx**4*dy**2-8.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+6.4d2*a1**2*a2**3*d2*dx**2*dy**2+6.4d2*a1**3*a2**2*d2*dx**2*dy**2 &
+1.3d2*a1*a2**3*dx**2*dy**2+2.6d2*a1**2*a2**2*dx**2*dy**2 &
+1.3d2*a1**3*a2*dx**2*dy**2+4.4d1*a1**2*a2**3*d2**2*dy**2 &
+4.4d1*a1**3*a2**2*d2**2*dy**2-3.16d2*a1*a2**3*d2*dy**2 &
-6.32d2*a1**2*a2**2*d2*dy**2-3.16d2*a1**3*a2*d2*dy**2 &
+1.65d2*a2**3*dy**2+4.95d2*a1*a2**2*dy**2+4.95d2*a1**2*a2*dy**2 &
+1.65d2*a1**3*dy**2+5.6d1*a1**3*a2**3*d2*dx**6-4.2d2*a1**2*a2**3*dx**6 &
-4.2d2*a1**3*a2**2*dx**6-8.8d1*a1**3*a2**3*d2**2*dx**4 &
+6.68d2*a1**2*a2**3*d2*dx**4+6.68d2*a1**3*a2**2*d2*dx**4 &
-5.2d1*a1*a2**3*dx**4-1.04d2*a1**2*a2**2*dx**4-5.2d1*a1**3*a2*dx**4 &
+3.2d1*a1**3*a2**3*d2**3*dx**2-2.2d2*a1**2*a2**3*d2**2*dx**2 &
-2.2d2*a1**3*a2**2*d2**2*dx**2-1.72d2*a1*a2**3*d2*dx**2 &
-3.44d2*a1**2*a2**2*d2*dx**2-1.72d2*a1**3*a2*d2*dx**2 &
+2.31d2*a2**3*dx**2+6.93d2*a1*a2**2*dx**2+6.93d2*a1**2*a2*dx**2 &
+2.31d2*a1**3*dx**2-1.6d1*a1**2*a2**3*d2**3-1.6d1*a1**3*a2**2*d2**3 &
+1.28d2*a1*a2**3*d2**2+2.56d2*a1**2*a2**2*d2**2+1.28d2*a1**3*a2*d2**2 &
-1.32d2*a2**3*d2-3.96d2*a1*a2**2*d2-3.96d2*a1**2*a2*d2-1.32d2*a1**3*d2 &
)
                case (2)
                  rlYlm_laplacian =5.75621273219899775d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.05d2*a1**2*a2**3*dy**4-1.05d2*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2**2*dy**2+1.14d2*a1**2*a2**3*d2*dy**2 &
+1.14d2*a1**3*a2**2*d2*dy**2-1.56d2*a1*a2**3*dy**2 &
-3.12d2*a1**2*a2**2*dy**2-1.56d2*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.05d2*a1**2*a2**3*dx**4 &
+1.05d2*a1**3*a2**2*dx**4+1.2d1*a1**3*a2**3*d2**2*dx**2 &
-8.6d1*a1**2*a2**3*d2*dx**2-8.6d1*a1**3*a2**2*d2*dx**2 &
-2.6d1*a1*a2**3*dx**2-5.2d1*a1**2*a2**2*dx**2-2.6d1*a1**3*a2*dx**2 &
-1.2d1*a1**2*a2**3*d2**2-1.2d1*a1**3*a2**2*d2**2+9.0d1*a1*a2**3*d2 &
+1.8d2*a1**2*a2**2*d2+9.0d1*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2 &
-1.98d2*a1**2*a2-6.6d1*a1**3)*dz
                case (3)
                  rlYlm_laplacian =5.3844439723186478d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(2.4d1*a1**3*a2**3*d2*dx**2*dy**4 &
-1.8d2*a1**2*a2**3*dx**2*dy**4-1.8d2*a1**3*a2**2*dx**2*dy**4 &
-1.2d1*a1**2*a2**3*d2*dy**4-1.2d1*a1**3*a2**2*d2*dy**4 &
+7.8d1*a1*a2**3*dy**4+1.56d2*a1**2*a2**2*dy**4+7.8d1*a1**3*a2*dy**4 &
+1.6d1*a1**3*a2**3*d2*dx**4*dy**2-1.2d2*a1**2*a2**3*dx**4*dy**2 &
-1.2d2*a1**3*a2**2*dx**4*dy**2-2.4d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+1.92d2*a1**2*a2**3*d2*dx**2*dy**2+1.92d2*a1**3*a2**2*d2*dx**2*dy**2 &
-7.8d1*a1*a2**3*dx**2*dy**2-1.56d2*a1**2*a2**2*dx**2*dy**2 &
-7.8d1*a1**3*a2*dx**2*dy**2+1.2d1*a1**2*a2**3*d2**2*dy**2 &
+1.2d1*a1**3*a2**2*d2**2*dy**2-8.4d1*a1*a2**3*d2*dy**2 &
-1.68d2*a1**2*a2**2*d2*dy**2-8.4d1*a1**3*a2*d2*dy**2+3.3d1*a2**3*dy**2 &
+9.9d1*a1*a2**2*dy**2+9.9d1*a1**2*a2*dy**2+3.3d1*a1**3*dy**2 &
-8.0d0*a1**3*a2**3*d2*dx**6+6.0d1*a1**2*a2**3*dx**6 &
+6.0d1*a1**3*a2**2*dx**6+8.0d0*a1**3*a2**3*d2**2*dx**4 &
-5.2d1*a1**2*a2**3*d2*dx**4-5.2d1*a1**3*a2**2*d2*dx**4 &
-5.2d1*a1*a2**3*dx**4-1.04d2*a1**2*a2**2*dx**4-5.2d1*a1**3*a2*dx**4 &
-1.2d1*a1**2*a2**3*d2**2*dx**2-1.2d1*a1**3*a2**2*d2**2*dx**2 &
+8.4d1*a1*a2**3*d2*dx**2+1.68d2*a1**2*a2**2*d2*dx**2 &
+8.4d1*a1**3*a2*d2*dx**2-3.3d1*a2**3*dx**2-9.9d1*a1*a2**2*dx**2 &
-9.9d1*a1**2*a2*dx**2-3.3d1*a1**3*dx**2)
                case (4)
                  rlYlm_laplacian =7.6147536914910937d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(2.0d0*a1**2*a2**2*d2*dy**4 &
-1.5d1*a1*a2**2*dy**4-1.5d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+9.0d1*a1*a2**2*dx**2*dy**2 &
+9.0d1*a1**2*a2*dx**2*dy**2+1.2d1*a1*a2**2*d2*dy**2 &
+1.2d1*a1**2*a2*d2*dy**2-7.8d1*a2**2*dy**2-1.56d2*a1*a2*dy**2 &
-7.8d1*a1**2*dy**2+2.0d0*a1**2*a2**2*d2*dx**4-1.5d1*a1*a2**2*dx**4 &
-1.5d1*a1**2*a2*dx**4-4.0d0*a1*a2**2*d2*dx**2-4.0d0*a1**2*a2*d2*dx**2 &
+2.6d1*a2**2*dx**2+5.2d1*a1*a2*dx**2+2.6d1*a1**2*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (2)
          ! selection on l2: l1=2, m1=2
          select case (l2)
            case (0)
              ! selection on m2: l1=2, m1=2, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-1.71617106161956686d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-6)*(1.0d0*a2* &
(2.0d0*a1*d2-7.0d0)-7.0d0*a1)*(dy+dx)*(dy-1.0d0*dx)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=2, m1=2, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)*dy* &
(2.0d0*a1**2*a2**2*d2*dy**2-9.0d0*a1*a2**2*dy**2-9.0d0*a1**2*a2*dy**2 &
-2.0d0*a1**2*a2**2*d2*dx**2+9.0d0*a1*a2**2*dx**2+9.0d0*a1**2*a2*dx**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+7.0d0*a2**2+1.4d1*a1*a2 &
+7.0d0*a1**2)
                case (0)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-7)* &
(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*(dy+dx)*(dy-1.0d0*dx)*dz
                case (1)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-7)*dx* &
(2.0d0*a1**2*a2**2*d2*dy**2-9.0d0*a1*a2**2*dy**2-9.0d0*a1**2*a2*dy**2 &
-2.0d0*a1**2*a2**2*d2*dx**2+9.0d0*a1*a2**2*dx**2+9.0d0*a1**2*a2*dx**2 &
+2.0d0*a1*a2**2*d2+2.0d0*a1**2*a2*d2-7.0d0*a2**2-1.4d1*a1*a2 &
-7.0d0*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=2, m1=2, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-6.6467019408956851d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-8)* &
(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)
                case (-1)
                  rlYlm_laplacian = &
-6.6467019408956851d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dy* &
(2.0d0*a1**2*a2**2*d2*dy**2-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2 &
-2.0d0*a1**2*a2**2*d2*dx**2+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)*dz
                case (0)
                  rlYlm_laplacian =1.91873757739966592d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(dy+dx)*(dy-1.0d0*dx)* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-3.3d1*a1*a2**2*dx**2-3.3d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.8d1*a1*a2**2*d2+1.8d1*a1**2*a2*d2 &
+1.8d1*a2**2+3.6d1*a1*a2+1.8d1*a1**2)
                case (1)
                  rlYlm_laplacian = &
-6.6467019408956851d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-8)*dx* &
(2.0d0*a1**2*a2**2*d2*dy**2-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2 &
-2.0d0*a1**2*a2**2*d2*dx**2+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2 &
+2.0d0*a1*a2**2*d2+2.0d0*a1**2*a2*d2-9.0d0*a2**2-1.8d1*a1*a2 &
-9.0d0*a1**2)*dz
                case (2)
                  rlYlm_laplacian =3.32335097044784255d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-8)*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.1d1*a1**2*a2**3*dy**4-1.1d1*a1**3*a2**2*dy**4 &
-4.0d0*a1**3*a2**3*d2*dx**2*dy**2+2.2d1*a1**2*a2**3*dx**2*dy**2 &
+2.2d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**2*a2**3*d2*dy**2 &
-4.0d0*a1**3*a2**2*d2*dy**2+1.8d1*a1*a2**3*dy**2 &
+3.6d1*a1**2*a2**2*dy**2+1.8d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.1d1*a1**2*a2**3*dx**4 &
-1.1d1*a1**3*a2**2*dx**4-4.0d0*a1**2*a2**3*d2*dx**2 &
-4.0d0*a1**3*a2**2*d2*dx**2+1.8d1*a1*a2**3*dx**2 &
+3.6d1*a1**2*a2**2*dx**2+1.8d1*a1**3*a2*dx**2+2.0d0*a1*a2**3*d2 &
+4.0d0*a1**2*a2**2*d2+2.0d0*a1**3*a2*d2-7.0d0*a2**3-2.1d1*a1*a2**2 &
-2.1d1*a1**2*a2-7.0d0*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=2, m1=2, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =3.58962931487909853d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.3d1*a1**2*a2**3*dy**4-1.3d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2*dx**2*dy**2+5.2d1*a1**2*a2**3*dx**2*dy**2 &
+5.2d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.3d1*a1*a2**3*dy**2 &
+6.6d1*a1**2*a2**2*dy**2+3.3d1*a1**3*a2*dy**2 &
+6.0d0*a1**3*a2**3*d2*dx**4-3.9d1*a1**2*a2**3*dx**4 &
-3.9d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.3d1*a1*a2**3*dx**2 &
+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-2.7d1*a2**3-8.1d1*a1*a2**2 &
-8.1d1*a1**2*a2-2.7d1*a1**3)
                case (-2)
                  rlYlm_laplacian = &
-1.75855203743803178d1*E*a1**4*a2**3*sqrt(a2+a1)*(a2+a1)**(-9)* &
(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)*dz
                case (-1)
                  rlYlm_laplacian =2.78051491111693767d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(1.0d1*a1**3*a2**3*d2*dy**4 &
-6.5d1*a1**2*a2**3*dy**4-6.5d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+3.8d1*a1**2*a2**3*d2*dy**2 &
+3.8d1*a1**3*a2**2*d2*dy**2+7.7d1*a1*a2**3*dy**2 &
+1.54d2*a1**2*a2**2*dy**2+7.7d1*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+6.5d1*a1**2*a2**3*dx**4 &
+6.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-5.8d1*a1**2*a2**3*d2*dx**2-5.8d1*a1**3*a2**2*d2*dx**2 &
+3.3d1*a1*a2**3*dx**2+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.2d1*a1*a2**3*d2 &
-8.4d1*a1**2*a2**2*d2-4.2d1*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)
                case (0)
                  rlYlm_laplacian =2.27028091814553966d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(dy+dx)*(dy-1.0d0*dx)* &
(1.0d1*a1**2*a2**2*d2*dy**2-6.5d1*a1*a2**2*dy**2-6.5d1*a1**2*a2*dy**2 &
+1.0d1*a1**2*a2**2*d2*dx**2-6.5d1*a1*a2**2*dx**2-6.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.4d1*a1*a2**2*d2+1.4d1*a1**2*a2*d2 &
+6.6d1*a2**2+1.32d2*a1*a2+6.6d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =2.78051491111693767d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(1.0d1*a1**3*a2**3*d2*dy**4 &
-6.5d1*a1**2*a2**3*dy**4-6.5d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+5.8d1*a1**2*a2**3*d2*dy**2 &
+5.8d1*a1**3*a2**2*d2*dy**2-3.3d1*a1*a2**3*dy**2 &
-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+6.5d1*a1**2*a2**3*dx**4 &
+6.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-3.8d1*a1**2*a2**3*d2*dx**2-3.8d1*a1**3*a2**2*d2*dx**2 &
-7.7d1*a1*a2**3*dx**2-1.54d2*a1**2*a2**2*dx**2-7.7d1*a1**3*a2*dx**2 &
-8.0d0*a1**2*a2**3*d2**2-8.0d0*a1**3*a2**2*d2**2+4.2d1*a1*a2**3*d2 &
+8.4d1*a1**2*a2**2*d2+4.2d1*a1**3*a2*d2+9.0d0*a2**3+2.7d1*a1*a2**2 &
+2.7d1*a1**2*a2+9.0d0*a1**3)
                case (2)
                  rlYlm_laplacian =8.79276018719015889d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.3d1*a1**2*a2**3*dy**4-1.3d1*a1**3*a2**2*dy**4 &
-4.0d0*a1**3*a2**3*d2*dx**2*dy**2+2.6d1*a1**2*a2**3*dx**2*dy**2 &
+2.6d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**2*a2**3*d2*dy**2 &
-4.0d0*a1**3*a2**2*d2*dy**2+2.2d1*a1*a2**3*dy**2 &
+4.4d1*a1**2*a2**2*dy**2+2.2d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.3d1*a1**2*a2**3*dx**4 &
-1.3d1*a1**3*a2**2*dx**4-4.0d0*a1**2*a2**3*d2*dx**2 &
-4.0d0*a1**3*a2**2*d2*dx**2+2.2d1*a1*a2**3*dx**2 &
+4.4d1*a1**2*a2**2*dx**2+2.2d1*a1**3*a2*dx**2+2.0d0*a1*a2**3*d2 &
+4.0d0*a1**2*a2**2*d2+2.0d0*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)*dz
                case (3)
                  rlYlm_laplacian =3.58962931487909853d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(6.0d0*a1**3*a2**3*d2*dy**4 &
-3.9d1*a1**2*a2**3*dy**4-3.9d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2*dx**2*dy**2+5.2d1*a1**2*a2**3*dx**2*dy**2 &
+5.2d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.3d1*a1*a2**3*dy**2 &
+6.6d1*a1**2*a2**2*dy**2+3.3d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.3d1*a1**2*a2**3*dx**4 &
-1.3d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.3d1*a1*a2**3*dx**2 &
+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-2.7d1*a2**3-8.1d1*a1*a2**2 &
-8.1d1*a1**2*a2-2.7d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=2, m1=2, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =1.52295073829821874d1*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.5d1*a1**2*a2**3*dy**4-1.5d1*a1**3*a2**2*dy**4 &
-4.0d0*a1**3*a2**3*d2*dx**2*dy**2+3.0d1*a1**2*a2**3*dx**2*dy**2 &
+3.0d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**2*a2**3*d2*dy**2 &
-4.0d0*a1**3*a2**2*d2*dy**2+2.6d1*a1*a2**3*dy**2 &
+5.2d1*a1**2*a2**2*dy**2+2.6d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4-4.0d0*a1**2*a2**3*d2*dx**2 &
-4.0d0*a1**3*a2**2*d2*dx**2+2.6d1*a1*a2**3*dx**2 &
+5.2d1*a1**2*a2**2*dx**2+2.6d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)
                case (-3)
                  rlYlm_laplacian =1.07688879446372956d1*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.5d1*a1**2*a2**3*dy**4-1.5d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2*dx**2*dy**2+6.0d1*a1**2*a2**3*dx**2*dy**2 &
+6.0d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+6.0d0*a1**3*a2**3*d2*dx**4-4.5d1*a1**2*a2**3*dx**4 &
-4.5d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (-2)
                  rlYlm_laplacian =5.75621273219899775d0*E*a1**4*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**2*a2**2*d2*dy**2-1.05d2*a1*a2**2*dy**2 &
-1.05d2*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-1.05d2*a1*a2**2*dx**2-1.05d2*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2 &
+8.6d1*a1*a2**2*d2+8.6d1*a1**2*a2*d2+2.6d1*a2**2+5.2d1*a1*a2 &
+2.6d1*a1**2)
                case (-1)
                  rlYlm_laplacian =4.07025705689025559d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.05d2*a1**2*a2**3*dy**4-1.05d2*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+3.4d1*a1**2*a2**3*d2*dy**2 &
+3.4d1*a1**3*a2**2*d2*dy**2+1.69d2*a1*a2**3*dy**2 &
+3.38d2*a1**2*a2**2*dy**2+1.69d2*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.05d2*a1**2*a2**3*dx**4 &
+1.05d2*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-6.2d1*a1**2*a2**3*d2*dx**2-6.2d1*a1**3*a2**2*d2*dx**2 &
+1.3d1*a1*a2**3*dx**2+2.6d1*a1**2*a2**2*dx**2+1.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.6d1*a1*a2**3*d2 &
-9.2d1*a1**2*a2**2*d2-4.6d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian =-6.43564148107337574d-1*E*a1**3*a2*sqrt &
                  (a2+a1)*(a2+a1)**(-10)*(dy+dx)*(dy-1.0d0*dx)* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.25d2*a1**2*a2**3*dy**4 &
-5.25d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.05d3*a1**2*a2**3*dx**2*dy**2-1.05d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+4.8d2*a1**2*a2**3*d2*dy**2 &
+4.8d2*a1**3*a2**2*d2*dy**2+7.8d2*a1*a2**3*dy**2 &
+1.56d3*a1**2*a2**2*dy**2+7.8d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.25d2*a1**2*a2**3*dx**4 &
-5.25d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+4.8d2*a1**2*a2**3*d2*dx**2+4.8d2*a1**3*a2**2*d2*dx**2 &
+7.8d2*a1*a2**3*dx**2+1.56d3*a1**2*a2**2*dx**2+7.8d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-2.4d1*a1**2*a2**3*d2**2 &
-2.4d1*a1**3*a2**2*d2**2-6.12d2*a1*a2**3*d2-1.224d3*a1**2*a2**2*d2 &
-6.12d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2-1.98d2*a1**2*a2 &
-6.6d1*a1**3)
                case (1)
                  rlYlm_laplacian =4.07025705689025559d0*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.05d2*a1**2*a2**3*dy**4-1.05d2*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+6.2d1*a1**2*a2**3*d2*dy**2 &
+6.2d1*a1**3*a2**2*d2*dy**2-1.3d1*a1*a2**3*dy**2 &
-2.6d1*a1**2*a2**2*dy**2-1.3d1*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.05d2*a1**2*a2**3*dx**4 &
+1.05d2*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-3.4d1*a1**2*a2**3*d2*dx**2-3.4d1*a1**3*a2**2*d2*dx**2 &
-1.69d2*a1*a2**3*dx**2-3.38d2*a1**2*a2**2*dx**2-1.69d2*a1**3*a2*dx**2 &
-8.0d0*a1**2*a2**3*d2**2-8.0d0*a1**3*a2**2*d2**2+4.6d1*a1*a2**3*d2 &
+9.2d1*a1**2*a2**2*d2+4.6d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian = &
-2.87810636609949888d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)* &
(1.4d1*a1**3*a2**3*d2*dy**6-1.05d2*a1**2*a2**3*dy**6 &
-1.05d2*a1**3*a2**2*dy**6-1.4d1*a1**3*a2**3*d2*dx**2*dy**4 &
+1.05d2*a1**2*a2**3*dx**2*dy**4+1.05d2*a1**3*a2**2*dx**2*dy**4 &
-1.2d1*a1**3*a2**3*d2**2*dy**4+5.8d1*a1**2*a2**3*d2*dy**4 &
+5.8d1*a1**3*a2**2*d2*dy**4+2.08d2*a1*a2**3*dy**4 &
+4.16d2*a1**2*a2**2*dy**4+2.08d2*a1**3*a2*dy**4 &
-1.4d1*a1**3*a2**3*d2*dx**4*dy**2+1.05d2*a1**2*a2**3*dx**4*dy**2 &
+1.05d2*a1**3*a2**2*dx**4*dy**2+2.4d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-2.28d2*a1**2*a2**3*d2*dx**2*dy**2-2.28d2*a1**3*a2**2*d2*dx**2*dy**2 &
+3.12d2*a1*a2**3*dx**2*dy**2+6.24d2*a1**2*a2**2*dx**2*dy**2 &
+3.12d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.38d2*a1*a2**3*d2*dy**2 &
-2.76d2*a1**2*a2**2*d2*dy**2-1.38d2*a1**3*a2*d2*dy**2 &
-9.9d1*a2**3*dy**2-2.97d2*a1*a2**2*dy**2-2.97d2*a1**2*a2*dy**2 &
-9.9d1*a1**3*dy**2+1.4d1*a1**3*a2**3*d2*dx**6-1.05d2*a1**2*a2**3*dx**6 &
-1.05d2*a1**3*a2**2*dx**6-1.2d1*a1**3*a2**3*d2**2*dx**4 &
+5.8d1*a1**2*a2**3*d2*dx**4+5.8d1*a1**3*a2**2*d2*dx**4 &
+2.08d2*a1*a2**3*dx**4+4.16d2*a1**2*a2**2*dx**4+2.08d2*a1**3*a2*dx**4 &
+2.4d1*a1**2*a2**3*d2**2*dx**2+2.4d1*a1**3*a2**2*d2**2*dx**2 &
-1.38d2*a1*a2**3*d2*dx**2-2.76d2*a1**2*a2**2*d2*dx**2 &
-1.38d2*a1**3*a2*d2*dx**2-9.9d1*a2**3*dx**2-2.97d2*a1*a2**2*dx**2 &
-2.97d2*a1**2*a2*dx**2-9.9d1*a1**3*dx**2-1.2d1*a1*a2**3*d2**2 &
-2.4d1*a1**2*a2**2*d2**2-1.2d1*a1**3*a2*d2**2+6.6d1*a2**3*d2 &
+1.98d2*a1*a2**2*d2+1.98d2*a1**2*a2*d2+6.6d1*a1**3*d2)
                case (3)
                  rlYlm_laplacian =1.07688879446372956d1*E*a1**3*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(6.0d0*a1**3*a2**3*d2*dy**4 &
-4.5d1*a1**2*a2**3*dy**4-4.5d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2*dx**2*dy**2+6.0d1*a1**2*a2**3*dx**2*dy**2 &
+6.0d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (4)
                  rlYlm_laplacian = &
-3.80737684574554685d0*E*a1**3*a2*sqrt(a2+a1)*(a2+a1)**(-10)*(dy+dx)* &
(dy-1.0d0*dx)*(2.0d0*a1**3*a2**3*d2*dy**4-1.5d1*a1**2*a2**3*dy**4 &
-1.5d1*a1**3*a2**2*dy**4-1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**3*dx**2*dy**2+9.0d1*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**2*a2**3*d2*dy**2-8.0d0*a1**3*a2**2*d2*dy**2 &
+5.2d1*a1*a2**3*dy**2+1.04d2*a1**2*a2**2*dy**2+5.2d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4-8.0d0*a1**2*a2**3*d2*dx**2 &
-8.0d0*a1**3*a2**2*d2*dx**2+5.2d1*a1*a2**3*dx**2 &
+1.04d2*a1**2*a2**2*dx**2+5.2d1*a1**3*a2*dx**2+1.2d1*a1*a2**3*d2 &
+2.4d1*a1**2*a2**2*d2+1.2d1*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2 &
-1.98d2*a1**2*a2-6.6d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case default
          print *,'Error: rlYlm_overlap not implemented for l1=',l1 &
,'m1=',m1,'l2=',l2,'m2=',m2
          stop
      end select
    case (3)
      ! selection on m1: l1=3
      select case (m1)
        case (-3)
          ! selection on l2: l1=3, m1=-3
          select case (l2)
            case (0)
              ! selection on m2: l1=3, m1=-3, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =1.85367660741129178d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dy* &
(dy**2-3.0d0*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=3, m1=-3, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =1.60533103241913232d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*(4.0d0*a1**2*a2**2*d2*dy**4-2.2d1*a1*a2**2*dy**4 &
-2.2d1*a1**2*a2*dy**4-1.2d1*a1**2*a2**2*d2*dx**2*dy**2 &
+6.6d1*a1*a2**2*dx**2*dy**2+6.6d1*a1**2*a2*dx**2*dy**2 &
-6.0d0*a1*a2**2*d2*dy**2-6.0d0*a1**2*a2*d2*dy**2+2.7d1*a2**2*dy**2 &
+5.4d1*a1*a2*dy**2+2.7d1*a1**2*dy**2+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-2.7d1*a2**2*dx**2-5.4d1*a1*a2*dx**2 &
-2.7d1*a1**2*dx**2)
                case (0)
                  rlYlm_laplacian =3.21066206483826464d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dy* &
(dy**2-3.0d0*dx**2)*dz
                case (1)
                  rlYlm_laplacian =3.21066206483826464d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(2.0d0*a1**2*a2**2*d2*dy**2 &
-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2-6.0d0*a1**2*a2**2*d2*dx**2 &
+3.3d1*a1*a2**2*dx**2+3.3d1*a1**2*a2*dx**2+6.0d0*a1*a2**2*d2 &
+6.0d0*a1**2*a2*d2-2.7d1*a2**2-5.4d1*a1*a2-2.7d1*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=3, m1=-3, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =3.58962931487909853d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(4.0d0*a1**3*a2**3*d2*dy**4 &
-2.6d1*a1**2*a2**3*dy**4-2.6d1*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2*dx**2*dy**2+7.8d1*a1**2*a2**3*dx**2*dy**2 &
+7.8d1*a1**3*a2**2*dx**2*dy**2+6.0d0*a1**2*a2**3*d2*dy**2 &
+6.0d0*a1**3*a2**2*d2*dy**2-3.3d1*a1*a2**3*dy**2 &
-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dx**2+6.0d0*a1**3*a2**2*d2*dx**2 &
-3.3d1*a1*a2**3*dx**2-6.6d1*a1**2*a2**2*dx**2-3.3d1*a1**3*a2*dx**2 &
-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2-6.0d0*a1**3*a2*d2+2.7d1*a2**3 &
+8.1d1*a1*a2**2+8.1d1*a1**2*a2+2.7d1*a1**3)
                case (-1)
                  rlYlm_laplacian =3.58962931487909853d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-9)*(4.0d0*a1**2*a2**2*d2*dy**4-2.6d1*a1*a2**2*dy**4 &
-2.6d1*a1**2*a2*dy**4-1.2d1*a1**2*a2**2*d2*dx**2*dy**2 &
+7.8d1*a1*a2**2*dx**2*dy**2+7.8d1*a1**2*a2*dx**2*dy**2 &
-6.0d0*a1*a2**2*d2*dy**2-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2 &
+6.6d1*a1*a2*dy**2+3.3d1*a1**2*dy**2+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2 &
-3.3d1*a1**2*dx**2)*dz
                case (0)
                  rlYlm_laplacian = &
-2.07247345123641944d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(dy**2-3.0d0*dx**2)*(6.0d0*a1**2*a2**2*d2*dy**2-3.9d1*a1*a2**2*dy**2 &
-3.9d1*a1**2*a2*dy**2+6.0d0*a1**2*a2**2*d2*dx**2-3.9d1*a1*a2**2*dx**2 &
-3.9d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+2.0d1*a1*a2**2*d2 &
+2.0d1*a1**2*a2*d2+3.3d1*a2**2+6.6d1*a1*a2+3.3d1*a1**2)
                case (1)
                  rlYlm_laplacian =7.17925862975819707d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(2.0d0*a1**2*a2**2*d2*dy**2 &
-1.3d1*a1*a2**2*dy**2-1.3d1*a1**2*a2*dy**2-6.0d0*a1**2*a2**2*d2*dx**2 &
+3.9d1*a1*a2**2*dx**2+3.9d1*a1**2*a2*dx**2+6.0d0*a1*a2**2*d2 &
+6.0d0*a1**2*a2*d2-3.3d1*a2**2-6.6d1*a1*a2-3.3d1*a1**2)*dz
                case (2)
                  rlYlm_laplacian = &
-3.58962931487909853d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(2.0d0*a1**3*a2**3*d2*dy**4-1.3d1*a1**2*a2**3*dy**4 &
-1.3d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
+5.2d1*a1**2*a2**3*dx**2*dy**2+5.2d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.3d1*a1*a2**3*dy**2+6.6d1*a1**2*a2**2*dy**2+3.3d1*a1**3*a2*dy**2 &
+6.0d0*a1**3*a2**3*d2*dx**4-3.9d1*a1**2*a2**3*dx**4 &
-3.9d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.3d1*a1*a2**3*dx**2 &
+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-2.7d1*a2**3-8.1d1*a1*a2**2 &
-8.1d1*a1**2*a2-2.7d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=3, m1=-3, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =-1.93862139942790816d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(4.0d0*a1**4*a2**4*d2*dy**6 &
-3.0d1*a1**3*a2**4*dy**6-3.0d1*a1**4*a2**3*dy**6 &
-2.4d1*a1**4*a2**4*d2*dx**2*dy**4+1.8d2*a1**3*a2**4*dx**2*dy**4 &
+1.8d2*a1**4*a2**3*dx**2*dy**4-1.8d1*a1**3*a2**4*d2*dy**4 &
-1.8d1*a1**4*a2**3*d2*dy**4+1.17d2*a1**2*a2**4*dy**4 &
+2.34d2*a1**3*a2**3*dy**4+1.17d2*a1**4*a2**2*dy**4 &
+3.6d1*a1**4*a2**4*d2*dx**4*dy**2-2.7d2*a1**3*a2**4*dx**4*dy**2 &
-2.7d2*a1**4*a2**3*dx**4*dy**2-3.6d1*a1**3*a2**4*d2*dx**2*dy**2 &
-3.6d1*a1**4*a2**3*d2*dx**2*dy**2+2.34d2*a1**2*a2**4*dx**2*dy**2 &
+4.68d2*a1**3*a2**3*dx**2*dy**2+2.34d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**2*a2**4*d2*dy**2+7.2d1*a1**3*a2**3*d2*dy**2 &
+3.6d1*a1**4*a2**2*d2*dy**2-1.98d2*a1*a2**4*dy**2 &
-5.94d2*a1**2*a2**3*dy**2-5.94d2*a1**3*a2**2*dy**2 &
-1.98d2*a1**4*a2*dy**2-1.8d1*a1**3*a2**4*d2*dx**4 &
-1.8d1*a1**4*a2**3*d2*dx**4+1.17d2*a1**2*a2**4*dx**4 &
+2.34d2*a1**3*a2**3*dx**4+1.17d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**2*a2**4*d2*dx**2+7.2d1*a1**3*a2**3*d2*dx**2 &
+3.6d1*a1**4*a2**2*d2*dx**2-1.98d2*a1*a2**4*dx**2 &
-5.94d2*a1**2*a2**3*dx**2-5.94d2*a1**3*a2**2*dx**2 &
-1.98d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+5.4d1*a2**4+2.16d2*a1*a2**3 &
+3.24d2*a1**2*a2**2+2.16d2*a1**3*a2+5.4d1*a1**4)
                case (-2)
                  rlYlm_laplacian =9.49726646607726303d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(4.0d0*a1**3*a2**3*d2*dy**4 &
-3.0d1*a1**2*a2**3*dy**4-3.0d1*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2*dx**2*dy**2+9.0d1*a1**2*a2**3*dx**2*dy**2 &
+9.0d1*a1**3*a2**2*dx**2*dy**2+6.0d0*a1**2*a2**3*d2*dy**2 &
+6.0d0*a1**3*a2**2*d2*dy**2-3.9d1*a1*a2**3*dy**2 &
-7.8d1*a1**2*a2**2*dy**2-3.9d1*a1**3*a2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dx**2+6.0d0*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2-6.0d0*a1**3*a2*d2+3.3d1*a2**3 &
+9.9d1*a1*a2**2+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian = &
-1.50164967891712101d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)* &
(2.0d1*a1**3*a2**3*d2*dy**6-1.5d2*a1**2*a2**3*dy**6 &
-1.5d2*a1**3*a2**2*dy**6-4.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+3.0d2*a1**2*a2**3*dx**2*dy**4+3.0d2*a1**3*a2**2*dx**2*dy**4 &
-1.6d1*a1**3*a2**3*d2**2*dy**4+7.8d1*a1**2*a2**3*d2*dy**4 &
+7.8d1*a1**3*a2**2*d2*dy**4+2.73d2*a1*a2**3*dy**4 &
+5.46d2*a1**2*a2**2*dy**4+2.73d2*a1**3*a2*dy**4 &
-6.0d1*a1**3*a2**3*d2*dx**4*dy**2+4.5d2*a1**2*a2**3*dx**4*dy**2 &
+4.5d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-3.24d2*a1**2*a2**3*d2*dx**2*dy**2-3.24d2*a1**3*a2**2*d2*dx**2*dy**2 &
-2.34d2*a1*a2**3*dx**2*dy**2-4.68d2*a1**2*a2**2*dx**2*dy**2 &
-2.34d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.44d2*a1*a2**3*d2*dy**2 &
-2.88d2*a1**2*a2**2*d2*dy**2-1.44d2*a1**3*a2*d2*dy**2 &
-6.6d1*a2**3*dy**2-1.98d2*a1*a2**2*dy**2-1.98d2*a1**2*a2*dy**2 &
-6.6d1*a1**3*dy**2+3.0d1*a1**2*a2**3*d2*dx**4 &
+3.0d1*a1**3*a2**2*d2*dx**4-1.95d2*a1*a2**3*dx**4 &
-3.9d2*a1**2*a2**2*dx**4-1.95d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.44d2*a1*a2**3*d2*dx**2+2.88d2*a1**2*a2**2*d2*dx**2 &
+1.44d2*a1**3*a2*d2*dx**2+6.6d1*a2**3*dx**2+1.98d2*a1*a2**2*dx**2 &
+1.98d2*a1**2*a2*dx**2+6.6d1*a1**3*dx**2)
                case (0)
                  rlYlm_laplacian = &
-2.45218365717409381d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(dy**2-3.0d0*dx**2)*(1.0d1*a1**2*a2**2*d2*dy**2-7.5d1*a1*a2**2*dy**2 &
-7.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2-7.5d1*a1*a2**2*dx**2 &
-7.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+1.2d1*a1*a2**2*d2 &
+1.2d1*a1**2*a2*d2+1.17d2*a2**2+2.34d2*a1*a2+1.17d2*a1**2)*dz
                case (1)
                  rlYlm_laplacian = &
-3.00329935783424201d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(1.0d1*a1**3*a2**3*d2*dy**4-7.5d1*a1**2*a2**3*dy**4 &
-7.5d1*a1**3*a2**2*dy**4-2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
+1.5d2*a1**2*a2**3*dx**2*dy**2+1.5d2*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+8.4d1*a1**2*a2**3*d2*dy**2 &
+8.4d1*a1**3*a2**2*d2*dy**2-1.56d2*a1*a2**3*dy**2 &
-3.12d2*a1**2*a2**2*dy**2-1.56d2*a1**3*a2*dy**2 &
-3.0d1*a1**3*a2**3*d2*dx**4+2.25d2*a1**2*a2**3*dx**4 &
+2.25d2*a1**3*a2**2*dx**4+2.4d1*a1**3*a2**3*d2**2*dx**2 &
-1.32d2*a1**2*a2**3*d2*dx**2-1.32d2*a1**3*a2**2*d2*dx**2 &
-3.12d2*a1*a2**3*dx**2-6.24d2*a1**2*a2**2*dx**2-3.12d2*a1**3*a2*dx**2 &
-2.4d1*a1**2*a2**3*d2**2-2.4d1*a1**3*a2**2*d2**2+1.44d2*a1*a2**3*d2 &
+2.88d2*a1**2*a2**2*d2+1.44d2*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)
                case (2)
                  rlYlm_laplacian = &
-9.49726646607726303d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(2.0d0*a1**3*a2**3*d2*dy**4-1.5d1*a1**2*a2**3*dy**4 &
-1.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
+6.0d1*a1**2*a2**3*dx**2*dy**2+6.0d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+6.0d0*a1**3*a2**3*d2*dx**4-4.5d1*a1**2*a2**3*dx**4 &
-4.5d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (3)
                  rlYlm_laplacian = &
-3.87724279885581631d0*E*a1**4*a2**4*sqrt(a2+a1)*(a2+a1)**(-10)* &
(1.0d0*a2*(2.0d0*a1*d2-1.5d1)-1.5d1*a1)*dx*dy*(dy**2-3.0d0*dx**2)* &
(3.0d0*dy**2-1.0d0*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=3, m1=-3, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-8.2248740261329704d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(4.0d0*a1**4*a2**4*d2*dy**6-3.4d1*a1**3*a2**4*dy**6 &
-3.4d1*a1**4*a2**3*dy**6-1.6d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.36d2*a1**3*a2**4*dx**2*dy**4+1.36d2*a1**4*a2**3*dx**2*dy**4 &
-6.0d0*a1**3*a2**4*d2*dy**4-6.0d0*a1**4*a2**3*d2*dy**4 &
+4.5d1*a1**2*a2**4*dy**4+9.0d1*a1**3*a2**3*dy**4 &
+4.5d1*a1**4*a2**2*dy**4+1.2d1*a1**4*a2**4*d2*dx**4*dy**2 &
-1.02d2*a1**3*a2**4*dx**4*dy**2-1.02d2*a1**4*a2**3*dx**4*dy**2 &
-1.2d1*a1**3*a2**4*d2*dx**2*dy**2-1.2d1*a1**4*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**4*dx**2*dy**2+1.8d2*a1**3*a2**3*dx**2*dy**2 &
+9.0d1*a1**4*a2**2*dx**2*dy**2+1.8d1*a1**2*a2**4*d2*dy**2 &
+3.6d1*a1**3*a2**3*d2*dy**2+1.8d1*a1**4*a2**2*d2*dy**2 &
-1.17d2*a1*a2**4*dy**2-3.51d2*a1**2*a2**3*dy**2 &
-3.51d2*a1**3*a2**2*dy**2-1.17d2*a1**4*a2*dy**2 &
-6.0d0*a1**3*a2**4*d2*dx**4-6.0d0*a1**4*a2**3*d2*dx**4 &
+4.5d1*a1**2*a2**4*dx**4+9.0d1*a1**3*a2**3*dx**4 &
+4.5d1*a1**4*a2**2*dx**4+1.8d1*a1**2*a2**4*d2*dx**2 &
+3.6d1*a1**3*a2**3*d2*dx**2+1.8d1*a1**4*a2**2*d2*dx**2 &
-1.17d2*a1*a2**4*dx**2-3.51d2*a1**2*a2**3*dx**2 &
-3.51d2*a1**3*a2**2*dx**2-1.17d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2 &
-3.6d1*a1**2*a2**3*d2-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2 &
+6.6d1*a2**4+2.64d2*a1*a2**3+3.96d2*a1**2*a2**2+2.64d2*a1**3*a2 &
+6.6d1*a1**4)
                case (-3)
                  rlYlm_laplacian = &
-5.81586419828372447d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(4.0d0*a1**4*a2**4*d2*dy**6-3.4d1*a1**3*a2**4*dy**6 &
-3.4d1*a1**4*a2**3*dy**6-2.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
+2.04d2*a1**3*a2**4*dx**2*dy**4+2.04d2*a1**4*a2**3*dx**2*dy**4 &
-1.8d1*a1**3*a2**4*d2*dy**4-1.8d1*a1**4*a2**3*d2*dy**4 &
+1.35d2*a1**2*a2**4*dy**4+2.7d2*a1**3*a2**3*dy**4 &
+1.35d2*a1**4*a2**2*dy**4+3.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
-3.06d2*a1**3*a2**4*dx**4*dy**2-3.06d2*a1**4*a2**3*dx**4*dy**2 &
-3.6d1*a1**3*a2**4*d2*dx**2*dy**2-3.6d1*a1**4*a2**3*d2*dx**2*dy**2 &
+2.7d2*a1**2*a2**4*dx**2*dy**2+5.4d2*a1**3*a2**3*dx**2*dy**2 &
+2.7d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**2*a2**4*d2*dy**2 &
+7.2d1*a1**3*a2**3*d2*dy**2+3.6d1*a1**4*a2**2*d2*dy**2 &
-2.34d2*a1*a2**4*dy**2-7.02d2*a1**2*a2**3*dy**2 &
-7.02d2*a1**3*a2**2*dy**2-2.34d2*a1**4*a2*dy**2 &
-1.8d1*a1**3*a2**4*d2*dx**4-1.8d1*a1**4*a2**3*d2*dx**4 &
+1.35d2*a1**2*a2**4*dx**4+2.7d2*a1**3*a2**3*dx**4 &
+1.35d2*a1**4*a2**2*dx**4+3.6d1*a1**2*a2**4*d2*dx**2 &
+7.2d1*a1**3*a2**3*d2*dx**2+3.6d1*a1**4*a2**2*d2*dx**2 &
-2.34d2*a1*a2**4*dx**2-7.02d2*a1**2*a2**3*dx**2 &
-7.02d2*a1**3*a2**2*dx**2-2.34d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2 &
-3.6d1*a1**2*a2**3*d2-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2 &
+6.6d1*a2**4+2.64d2*a1*a2**3+3.96d2*a1**2*a2**2+2.64d2*a1**3*a2 &
+6.6d1*a1**4)*dz
                case (-2)
                  rlYlm_laplacian = &
-3.10871017685462917d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(2.8d1*a1**4*a2**4*d2*dy**6-2.38d2*a1**3*a2**4*dy**6 &
-2.38d2*a1**4*a2**3*dy**6-5.6d1*a1**4*a2**4*d2*dx**2*dy**4 &
+4.76d2*a1**3*a2**4*dx**2*dy**4+4.76d2*a1**4*a2**3*dx**2*dy**4 &
-2.4d1*a1**4*a2**4*d2**2*dy**4+2.34d2*a1**3*a2**4*d2*dy**4 &
+2.34d2*a1**4*a2**3*d2*dy**4-2.25d2*a1**2*a2**4*dy**4 &
-4.5d2*a1**3*a2**3*dy**4-2.25d2*a1**4*a2**2*dy**4 &
-8.4d1*a1**4*a2**4*d2*dx**4*dy**2+7.14d2*a1**3*a2**4*dx**4*dy**2 &
+7.14d2*a1**4*a2**3*dx**4*dy**2+7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-4.92d2*a1**3*a2**4*d2*dx**2*dy**2-4.92d2*a1**4*a2**3*d2*dx**2*dy**2 &
-9.0d2*a1**2*a2**4*dx**2*dy**2-1.8d3*a1**3*a2**3*dx**2*dy**2 &
-9.0d2*a1**4*a2**2*dx**2*dy**2-3.6d1*a1**3*a2**4*d2**2*dy**2 &
-3.6d1*a1**4*a2**3*d2**2*dy**2+2.16d2*a1**2*a2**4*d2*dy**2 &
+4.32d2*a1**3*a2**3*d2*dy**2+2.16d2*a1**4*a2**2*d2*dy**2 &
+3.51d2*a1*a2**4*dy**2+1.053d3*a1**2*a2**3*dy**2 &
+1.053d3*a1**3*a2**2*dy**2+3.51d2*a1**4*a2*dy**2 &
+4.2d1*a1**3*a2**4*d2*dx**4+4.2d1*a1**4*a2**3*d2*dx**4 &
-3.15d2*a1**2*a2**4*dx**4-6.3d2*a1**3*a2**3*dx**4 &
-3.15d2*a1**4*a2**2*dx**4-3.6d1*a1**3*a2**4*d2**2*dx**2 &
-3.6d1*a1**4*a2**3*d2**2*dx**2+2.16d2*a1**2*a2**4*d2*dx**2 &
+4.32d2*a1**3*a2**3*d2*dx**2+2.16d2*a1**4*a2**2*d2*dx**2 &
+3.51d2*a1*a2**4*dx**2+1.053d3*a1**2*a2**3*dx**2 &
+1.053d3*a1**3*a2**2*dx**2+3.51d2*a1**4*a2*dx**2 &
+3.6d1*a1**2*a2**4*d2**2+7.2d1*a1**3*a2**3*d2**2 &
+3.6d1*a1**4*a2**2*d2**2-2.28d2*a1*a2**4*d2-6.84d2*a1**2*a2**3*d2 &
-6.84d2*a1**3*a2**2*d2-2.28d2*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (-1)
                  rlYlm_laplacian = &
-2.19819004679753972d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(2.8d1*a1**3*a2**3*d2*dy**6-2.38d2*a1**2*a2**3*dy**6 &
-2.38d2*a1**3*a2**2*dy**6-5.6d1*a1**3*a2**3*d2*dx**2*dy**4 &
+4.76d2*a1**2*a2**3*dx**2*dy**4+4.76d2*a1**3*a2**2*dx**2*dy**4 &
-1.6d1*a1**3*a2**3*d2**2*dy**4+5.8d1*a1**2*a2**3*d2*dy**4 &
+5.8d1*a1**3*a2**2*d2*dy**4+5.85d2*a1*a2**3*dy**4 &
+1.17d3*a1**2*a2**2*dy**4+5.85d2*a1**3*a2*dy**4 &
-8.4d1*a1**3*a2**3*d2*dx**4*dy**2+7.14d2*a1**2*a2**3*dx**4*dy**2 &
+7.14d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-3.0d2*a1**2*a2**3*d2*dx**2*dy**2-3.0d2*a1**3*a2**2*d2*dx**2*dy**2 &
-8.1d2*a1*a2**3*dx**2*dy**2-1.62d3*a1**2*a2**2*dx**2*dy**2 &
-8.1d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.44d2*a1*a2**3*d2*dy**2 &
-2.88d2*a1**2*a2**2*d2*dy**2-1.44d2*a1**3*a2*d2*dy**2 &
-2.34d2*a2**3*dy**2-7.02d2*a1*a2**2*dy**2-7.02d2*a1**2*a2*dy**2 &
-2.34d2*a1**3*dy**2+4.2d1*a1**2*a2**3*d2*dx**4 &
+4.2d1*a1**3*a2**2*d2*dx**4-3.15d2*a1*a2**3*dx**4 &
-6.3d2*a1**2*a2**2*dx**4-3.15d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.44d2*a1*a2**3*d2*dx**2+2.88d2*a1**2*a2**2*d2*dx**2 &
+1.44d2*a1**3*a2*d2*dx**2+2.34d2*a2**3*dx**2+7.02d2*a1*a2**2*dx**2 &
+7.02d2*a1**2*a2*dx**2+2.34d2*a1**3*dx**2)*dz
                case (0)
                  rlYlm_laplacian =6.95128727779234418d-1*E*a1**3*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*dy*(dy**2-3.0d0*dx**2)* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+5.0d2*a1**2*a2**3*d2*dy**2 &
+5.0d2*a1**3*a2**2*d2*dy**2+1.35d3*a1*a2**3*dy**2 &
+2.7d3*a1**2*a2**2*dy**2+1.35d3*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+5.0d2*a1**2*a2**3*d2*dx**2+5.0d2*a1**3*a2**2*d2*dx**2 &
+1.35d3*a1*a2**3*dx**2+2.7d3*a1**2*a2**2*dx**2+1.35d3*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3+8.0d0*a1**2*a2**3*d2**2 &
+8.0d0*a1**3*a2**2*d2**2-1.044d3*a1*a2**3*d2-2.088d3*a1**2*a2**2*d2 &
-1.044d3*a1**3*a2*d2-2.34d2*a2**3-7.02d2*a1*a2**2-7.02d2*a1**2*a2 &
-2.34d2*a1**3)
                case (1)
                  rlYlm_laplacian = &
-4.39638009359507944d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(1.4d1*a1**3*a2**3*d2*dy**4-1.19d2*a1**2*a2**3*dy**4 &
-1.19d2*a1**3*a2**2*dy**4-2.8d1*a1**3*a2**3*d2*dx**2*dy**2 &
+2.38d2*a1**2*a2**3*dx**2*dy**2+2.38d2*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+9.2d1*a1**2*a2**3*d2*dy**2 &
+9.2d1*a1**3*a2**2*d2*dy**2-1.8d2*a1*a2**3*dy**2 &
-3.6d2*a1**2*a2**2*dy**2-1.8d2*a1**3*a2*dy**2 &
-4.2d1*a1**3*a2**3*d2*dx**4+3.57d2*a1**2*a2**3*dx**4 &
+3.57d2*a1**3*a2**2*dx**4+2.4d1*a1**3*a2**3*d2**2*dx**2 &
-1.08d2*a1**2*a2**3*d2*dx**2-1.08d2*a1**3*a2**2*d2*dx**2 &
-7.2d2*a1*a2**3*dx**2-1.44d3*a1**2*a2**2*dx**2-7.2d2*a1**3*a2*dx**2 &
-2.4d1*a1**2*a2**3*d2**2-2.4d1*a1**3*a2**2*d2**2+1.44d2*a1*a2**3*d2 &
+2.88d2*a1**2*a2**2*d2+1.44d2*a1**3*a2*d2+2.34d2*a2**3+7.02d2*a1*a2**2 &
+7.02d2*a1**2*a2+2.34d2*a1**3)*dz
                case (2)
                  rlYlm_laplacian =3.10871017685462917d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(1.4d1*a1**4*a2**4*d2*dy**6 &
-1.19d2*a1**3*a2**4*dy**6-1.19d2*a1**4*a2**3*dy**6 &
-4.2d1*a1**4*a2**4*d2*dx**2*dy**4+3.57d2*a1**3*a2**4*dx**2*dy**4 &
+3.57d2*a1**4*a2**3*dx**2*dy**4-1.2d1*a1**4*a2**4*d2**2*dy**4 &
+5.4d1*a1**3*a2**4*d2*dy**4+5.4d1*a1**4*a2**3*d2*dy**4 &
+3.6d2*a1**2*a2**4*dy**4+7.2d2*a1**3*a2**3*dy**4 &
+3.6d2*a1**4*a2**2*dy**4-1.4d1*a1**4*a2**4*d2*dx**4*dy**2 &
+1.19d2*a1**3*a2**4*dx**4*dy**2+1.19d2*a1**4*a2**3*dx**4*dy**2 &
+4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2-4.68d2*a1**3*a2**4*d2*dx**2*dy**2 &
-4.68d2*a1**4*a2**3*d2*dx**2*dy**2+4.5d2*a1**2*a2**4*dx**2*dy**2 &
+9.0d2*a1**3*a2**3*dx**2*dy**2+4.5d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**3*a2**4*d2**2*dy**2+3.6d1*a1**4*a2**3*d2**2*dy**2 &
-2.16d2*a1**2*a2**4*d2*dy**2-4.32d2*a1**3*a2**3*d2*dy**2 &
-2.16d2*a1**4*a2**2*d2*dy**2-3.51d2*a1*a2**4*dy**2 &
-1.053d3*a1**2*a2**3*dy**2-1.053d3*a1**3*a2**2*dy**2 &
-3.51d2*a1**4*a2*dy**2+4.2d1*a1**4*a2**4*d2*dx**6 &
-3.57d2*a1**3*a2**4*dx**6-3.57d2*a1**4*a2**3*dx**6 &
-3.6d1*a1**4*a2**4*d2**2*dx**4+2.46d2*a1**3*a2**4*d2*dx**4 &
+2.46d2*a1**4*a2**3*d2*dx**4+4.5d2*a1**2*a2**4*dx**4 &
+9.0d2*a1**3*a2**3*dx**4+4.5d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**3*a2**4*d2**2*dx**2+3.6d1*a1**4*a2**3*d2**2*dx**2 &
-2.16d2*a1**2*a2**4*d2*dx**2-4.32d2*a1**3*a2**3*d2*dx**2 &
-2.16d2*a1**4*a2**2*d2*dx**2-3.51d2*a1*a2**4*dx**2 &
-1.053d3*a1**2*a2**3*dx**2-1.053d3*a1**3*a2**2*dx**2 &
-3.51d2*a1**4*a2*dx**2-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.28d2*a1*a2**4*d2+6.84d2*a1**2*a2**3*d2 &
+6.84d2*a1**3*a2**2*d2+2.28d2*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case (3)
                  rlYlm_laplacian = &
-1.16317283965674489d1*E*a1**5*a2**4*sqrt(a2+a1)*(a2+a1)**(-11)* &
(1.0d0*a2*(2.0d0*a1*d2-1.7d1)-1.7d1*a1)*dx*dy*(dy**2-3.0d0*dx**2)* &
(3.0d0*dy**2-1.0d0*dx**2)*dz
                case (4)
                  rlYlm_laplacian =4.1124370130664852d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(2.0d0*a1**4*a2**4*d2*dy**6 &
-1.7d1*a1**3*a2**4*dy**6-1.7d1*a1**4*a2**3*dy**6 &
-1.8d1*a1**4*a2**4*d2*dx**2*dy**4+1.53d2*a1**3*a2**4*dx**2*dy**4 &
+1.53d2*a1**4*a2**3*dx**2*dy**4-1.2d1*a1**3*a2**4*d2*dy**4 &
-1.2d1*a1**4*a2**3*d2*dy**4+9.0d1*a1**2*a2**4*dy**4 &
+1.8d2*a1**3*a2**3*dy**4+9.0d1*a1**4*a2**2*dy**4 &
+3.8d1*a1**4*a2**4*d2*dx**4*dy**2-3.23d2*a1**3*a2**4*dx**4*dy**2 &
-3.23d2*a1**4*a2**3*dx**4*dy**2-2.4d1*a1**3*a2**4*d2*dx**2*dy**2 &
-2.4d1*a1**4*a2**3*d2*dx**2*dy**2+1.8d2*a1**2*a2**4*dx**2*dy**2 &
+3.6d2*a1**3*a2**3*dx**2*dy**2+1.8d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**2*a2**4*d2*dy**2+7.2d1*a1**3*a2**3*d2*dy**2 &
+3.6d1*a1**4*a2**2*d2*dy**2-2.34d2*a1*a2**4*dy**2 &
-7.02d2*a1**2*a2**3*dy**2-7.02d2*a1**3*a2**2*dy**2 &
-2.34d2*a1**4*a2*dy**2-6.0d0*a1**4*a2**4*d2*dx**6 &
+5.1d1*a1**3*a2**4*dx**6+5.1d1*a1**4*a2**3*dx**6 &
-1.2d1*a1**3*a2**4*d2*dx**4-1.2d1*a1**4*a2**3*d2*dx**4 &
+9.0d1*a1**2*a2**4*dx**4+1.8d2*a1**3*a2**3*dx**4 &
+9.0d1*a1**4*a2**2*dx**4+3.6d1*a1**2*a2**4*d2*dx**2 &
+7.2d1*a1**3*a2**3*d2*dx**2+3.6d1*a1**4*a2**2*d2*dx**2 &
-2.34d2*a1*a2**4*dx**2-7.02d2*a1**2*a2**3*dx**2 &
-7.02d2*a1**3*a2**2*dx**2-2.34d2*a1**4*a2*dx**2-2.4d1*a1*a2**4*d2 &
-7.2d1*a1**2*a2**3*d2-7.2d1*a1**3*a2**2*d2-2.4d1*a1**4*a2*d2 &
+1.32d2*a2**4+5.28d2*a1*a2**3+7.92d2*a1**2*a2**2+5.28d2*a1**3*a2 &
+1.32d2*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (-2)
          ! selection on l2: l1=3, m1=-2
          select case (l2)
            case (0)
              ! selection on m2: l1=3, m1=-2, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-9.08112367258215863d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2* &
(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx*dy*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=3, m1=-2, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-7.86448379536438834d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-8)*dx* &
(4.0d0*a1**2*a2**2*d2*dy**2-2.2d1*a1*a2**2*dy**2-2.2d1*a1**2*a2*dy**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)*dz
                case (0)
                  rlYlm_laplacian =7.86448379536438834d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(4.0d0*a1**2*a2**2*d2*dy**2 &
-2.2d1*a1*a2**2*dy**2-2.2d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2 &
-2.2d1*a1*a2**2*dx**2-2.2d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+2.4d1*a1*a2**2*d2+2.4d1*a1**2*a2*d2-9.0d0*a2**2-1.8d1*a1*a2 &
-9.0d0*a1**2)
                case (1)
                  rlYlm_laplacian = &
-7.86448379536438834d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-8)* &
(4.0d0*a1**2*a2**2*d2*dx**2-2.2d1*a1*a2**2*dx**2-2.2d1*a1**2*a2*dx**2 &
-2.0d0*a1*a2**2*d2-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)*dy*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=3, m1=-2, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-8.79276018719015889d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)* &
(8.0d0*a1**3*a2**3*d2*dx**2*dy**2-5.2d1*a1**2*a2**3*dx**2*dy**2 &
-5.2d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**2*a2**3*d2*dy**2 &
-4.0d0*a1**3*a2**2*d2*dy**2+2.2d1*a1*a2**3*dy**2 &
+4.4d1*a1**2*a2**2*dy**2+2.2d1*a1**3*a2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+2.2d1*a1*a2**3*dx**2+4.4d1*a1**2*a2**2*dx**2+2.2d1*a1**3*a2*dx**2 &
+2.0d0*a1*a2**3*d2+4.0d0*a1**2*a2**2*d2+2.0d0*a1**3*a2*d2-9.0d0*a2**3 &
-2.7d1*a1*a2**2-2.7d1*a1**2*a2-9.0d0*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =8.79276018719015889d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(8.0d0*a1**3*a2**3*d2*dy**4 &
-5.2d1*a1**2*a2**3*dy**4-5.2d1*a1**3*a2**2*dy**4 &
+8.0d0*a1**3*a2**3*d2*dx**2*dy**2-5.2d1*a1**2*a2**3*dx**2*dy**2 &
-5.2d1*a1**3*a2**2*dx**2*dy**2-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+5.2d1*a1**2*a2**3*d2*dy**2+5.2d1*a1**3*a2**2*d2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+2.2d1*a1*a2**3*dx**2+4.4d1*a1**2*a2**2*dx**2+2.2d1*a1**3*a2*dx**2 &
+4.0d0*a1**2*a2**3*d2**2+4.0d0*a1**3*a2**2*d2**2-2.4d1*a1*a2**3*d2 &
-4.8d1*a1**2*a2**2*d2-2.4d1*a1**3*a2*d2+9.0d0*a2**3+2.7d1*a1*a2**2 &
+2.7d1*a1**2*a2+9.0d0*a1**3)
                case (0)
                  rlYlm_laplacian =1.01530049219881249d1*E*a1**3*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*dx*dy* &
(3.0d0*(dy**2+dx**2)-2.0d0*d2)*dz
                case (1)
                  rlYlm_laplacian =8.79276018719015889d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
-5.2d1*a1**2*a2**3*dx**2*dy**2-5.2d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+2.2d1*a1*a2**3*dy**2+4.4d1*a1**2*a2**2*dy**2+2.2d1*a1**3*a2*dy**2 &
+8.0d0*a1**3*a2**3*d2*dx**4-5.2d1*a1**2*a2**3*dx**4 &
-5.2d1*a1**3*a2**2*dx**4-8.0d0*a1**3*a2**3*d2**2*dx**2 &
+5.2d1*a1**2*a2**3*d2*dx**2+5.2d1*a1**3*a2**2*d2*dx**2 &
+4.0d0*a1**2*a2**3*d2**2+4.0d0*a1**3*a2**2*d2**2-2.4d1*a1*a2**3*d2 &
-4.8d1*a1**2*a2**2*d2-2.4d1*a1**3*a2*d2+9.0d0*a2**3+2.7d1*a1*a2**2 &
+2.7d1*a1**2*a2+9.0d0*a1**3)
                case (2)
                  rlYlm_laplacian =1.75855203743803178d1*E*a1**3*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=3, m1=-2, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =9.49726646607726303d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(4.0d0*a1**3*a2**3*d2*dy**4 &
-3.0d1*a1**2*a2**3*dy**4-3.0d1*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2*dx**2*dy**2+9.0d1*a1**2*a2**3*dx**2*dy**2 &
+9.0d1*a1**3*a2**2*dx**2*dy**2+6.0d0*a1**2*a2**3*d2*dy**2 &
+6.0d0*a1**3*a2**2*d2*dy**2-3.9d1*a1*a2**3*dy**2 &
-7.8d1*a1**2*a2**2*dy**2-3.9d1*a1**3*a2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dx**2+6.0d0*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2-6.0d0*a1**3*a2*d2+3.3d1*a2**3 &
+9.9d1*a1*a2**2+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (-2)
                  rlYlm_laplacian =1.16317283965674489d1*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(1.6d1*a1**4*a2**4*d2*dx**2*dy**4 &
-1.2d2*a1**3*a2**4*dx**2*dy**4-1.2d2*a1**4*a2**3*dx**2*dy**4 &
-8.0d0*a1**3*a2**4*d2*dy**4-8.0d0*a1**4*a2**3*d2*dy**4 &
+5.2d1*a1**2*a2**4*dy**4+1.04d2*a1**3*a2**3*dy**4 &
+5.2d1*a1**4*a2**2*dy**4+1.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
-1.2d2*a1**3*a2**4*dx**4*dy**2-1.2d2*a1**4*a2**3*dx**4*dy**2 &
-1.6d1*a1**4*a2**4*d2**2*dx**2*dy**2+1.12d2*a1**3*a2**4*d2*dx**2*dy**2 &
+1.12d2*a1**4*a2**3*d2*dx**2*dy**2+5.2d1*a1**2*a2**4*dx**2*dy**2 &
+1.04d2*a1**3*a2**3*dx**2*dy**2+5.2d1*a1**4*a2**2*dx**2*dy**2 &
+8.0d0*a1**3*a2**4*d2**2*dy**2+8.0d0*a1**4*a2**3*d2**2*dy**2 &
-5.2d1*a1**2*a2**4*d2*dy**2-1.04d2*a1**3*a2**3*d2*dy**2 &
-5.2d1*a1**4*a2**2*d2*dy**2-8.0d0*a1**3*a2**4*d2*dx**4 &
-8.0d0*a1**4*a2**3*d2*dx**4+5.2d1*a1**2*a2**4*dx**4 &
+1.04d2*a1**3*a2**3*dx**4+5.2d1*a1**4*a2**2*dx**4 &
+8.0d0*a1**3*a2**4*d2**2*dx**2+8.0d0*a1**4*a2**3*d2**2*dx**2 &
-5.2d1*a1**2*a2**4*d2*dx**2-1.04d2*a1**3*a2**3*d2*dx**2 &
-5.2d1*a1**4*a2**2*d2*dx**2-4.0d0*a1**2*a2**4*d2**2 &
-8.0d0*a1**3*a2**3*d2**2-4.0d0*a1**4*a2**2*d2**2+2.4d1*a1*a2**4*d2 &
+7.2d1*a1**2*a2**3*d2+7.2d1*a1**3*a2**2*d2+2.4d1*a1**4*a2*d2 &
-9.0d0*a2**4-3.6d1*a1*a2**3-5.4d1*a1**2*a2**2-3.6d1*a1**3*a2 &
-9.0d0*a1**4)
                case (-1)
                  rlYlm_laplacian =7.35655097152228142d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(2.0d1*a1**3*a2**3*d2*dy**4 &
-1.5d2*a1**2*a2**3*dy**4-1.5d2*a1**3*a2**2*dy**4 &
+2.0d1*a1**3*a2**3*d2*dx**2*dy**2-1.5d2*a1**2*a2**3*dx**2*dy**2 &
-1.5d2*a1**3*a2**2*dx**2*dy**2-1.6d1*a1**3*a2**3*d2**2*dy**2 &
+1.18d2*a1**2*a2**3*d2*dy**2+1.18d2*a1**3*a2**2*d2*dy**2 &
+1.3d1*a1*a2**3*dy**2+2.6d1*a1**2*a2**2*dy**2+1.3d1*a1**3*a2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dx**2-1.0d1*a1**3*a2**2*d2*dx**2 &
+6.5d1*a1*a2**3*dx**2+1.3d2*a1**2*a2**2*dx**2+6.5d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.8d1*a1*a2**3*d2 &
-1.16d2*a1**2*a2**2*d2-5.8d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian = &
-6.00659871566848402d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.5d2*a1**2*a2**3*dy**4 &
-1.5d2*a1**3*a2**2*dy**4+4.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.0d2*a1**2*a2**3*dx**2*dy**2-3.0d2*a1**3*a2**2*dx**2*dy**2 &
-2.8d1*a1**3*a2**3*d2**2*dy**2+2.04d2*a1**2*a2**3*d2*dy**2 &
+2.04d2*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.5d2*a1**2*a2**3*dx**4 &
-1.5d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+2.04d2*a1**2*a2**3*d2*dx**2+2.04d2*a1**3*a2**2*d2*dx**2 &
+3.9d1*a1*a2**3*dx**2+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-4.8d1*a1**2*a2**3*d2**2 &
-4.8d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2-1.8d2*a1**2*a2**2*d2 &
-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2+1.98d2*a1**2*a2 &
+6.6d1*a1**3)
                case (1)
                  rlYlm_laplacian =7.35655097152228142d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.5d2*a1**2*a2**3*dx**2*dy**2-1.5d2*a1**3*a2**2*dx**2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dy**2-1.0d1*a1**3*a2**2*d2*dy**2 &
+6.5d1*a1*a2**3*dy**2+1.3d2*a1**2*a2**2*dy**2+6.5d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.5d2*a1**2*a2**3*dx**4 &
-1.5d2*a1**3*a2**2*dx**4-1.6d1*a1**3*a2**3*d2**2*dx**2 &
+1.18d2*a1**2*a2**3*d2*dx**2+1.18d2*a1**3*a2**2*d2*dx**2 &
+1.3d1*a1*a2**3*dx**2+2.6d1*a1**2*a2**2*dx**2+1.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.8d1*a1*a2**3*d2 &
-1.16d2*a1**2*a2**2*d2-5.8d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian = &
-2.32634567931348979d1*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(4.0d0*a1**2*a2**2*d2*dy**2-3.0d1*a1*a2**2*dy**2 &
-3.0d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2-3.0d1*a1*a2**2*dx**2 &
-3.0d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+3.2d1*a1*a2**2*d2 &
+3.2d1*a1**2*a2*d2-1.3d1*a2**2-2.6d1*a1*a2-1.3d1*a1**2)
                case (3)
                  rlYlm_laplacian =9.49726646607726303d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
-9.0d1*a1**2*a2**3*dx**2*dy**2-9.0d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+3.0d1*a1**2*a2**3*dx**4 &
+3.0d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=3, m1=-2, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =2.01467445626964921d1*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(dy+dx)*(dy-1.0d0*dx)* &
(8.0d0*a1**3*a2**3*d2*dx**2*dy**2-6.8d1*a1**2*a2**3*dx**2*dy**2 &
-6.8d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**2*a2**3*d2*dy**2 &
-4.0d0*a1**3*a2**2*d2*dy**2+3.0d1*a1*a2**3*dy**2 &
+6.0d1*a1**2*a2**2*dy**2+3.0d1*a1**3*a2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+3.0d1*a1*a2**3*dx**2+6.0d1*a1**2*a2**2*dx**2+3.0d1*a1**3*a2*dx**2 &
+6.0d0*a1*a2**3*d2+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.9d1*a2**3 &
-1.17d2*a1*a2**2-1.17d2*a1**2*a2-3.9d1*a1**3)*dz
                case (-3)
                  rlYlm_laplacian = &
-1.42458996991158946d1*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(8.0d0*a1**4*a2**4*d2*dy**6-6.8d1*a1**3*a2**4*dy**6 &
-6.8d1*a1**4*a2**3*dy**6-1.6d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.36d2*a1**3*a2**4*dx**2*dy**4+1.36d2*a1**4*a2**3*dx**2*dy**4 &
-8.0d0*a1**4*a2**4*d2**2*dy**4+8.4d1*a1**3*a2**4*d2*dy**4 &
+8.4d1*a1**4*a2**3*d2*dy**4-1.2d2*a1**2*a2**4*dy**4 &
-2.4d2*a1**3*a2**3*dy**4-1.2d2*a1**4*a2**2*dy**4 &
-2.4d1*a1**4*a2**4*d2*dx**4*dy**2+2.04d2*a1**3*a2**4*dx**4*dy**2 &
+2.04d2*a1**4*a2**3*dx**4*dy**2+2.4d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-1.92d2*a1**3*a2**4*d2*dx**2*dy**2-1.92d2*a1**4*a2**3*d2*dx**2*dy**2 &
-9.0d1*a1**2*a2**4*dx**2*dy**2-1.8d2*a1**3*a2**3*dx**2*dy**2 &
-9.0d1*a1**4*a2**2*dx**2*dy**2-1.2d1*a1**3*a2**4*d2**2*dy**2 &
-1.2d1*a1**4*a2**3*d2**2*dy**2+8.4d1*a1**2*a2**4*d2*dy**2 &
+1.68d2*a1**3*a2**3*d2*dy**2+8.4d1*a1**4*a2**2*d2*dy**2 &
+3.9d1*a1*a2**4*dy**2+1.17d2*a1**2*a2**3*dy**2 &
+1.17d2*a1**3*a2**2*dy**2+3.9d1*a1**4*a2*dy**2 &
+1.2d1*a1**3*a2**4*d2*dx**4+1.2d1*a1**4*a2**3*d2*dx**4 &
-9.0d1*a1**2*a2**4*dx**4-1.8d2*a1**3*a2**3*dx**4 &
-9.0d1*a1**4*a2**2*dx**4-1.2d1*a1**3*a2**4*d2**2*dx**2 &
-1.2d1*a1**4*a2**3*d2**2*dx**2+8.4d1*a1**2*a2**4*d2*dx**2 &
+1.68d2*a1**3*a2**3*d2*dx**2+8.4d1*a1**4*a2**2*d2*dx**2 &
+3.9d1*a1*a2**4*dx**2+1.17d2*a1**2*a2**3*dx**2 &
+1.17d2*a1**3*a2**2*dx**2+3.9d1*a1**4*a2*dx**2+1.2d1*a1**2*a2**4*d2**2 &
+2.4d1*a1**3*a2**3*d2**2+1.2d1*a1**4*a2**2*d2**2-8.4d1*a1*a2**4*d2 &
-2.52d2*a1**2*a2**3*d2-2.52d2*a1**3*a2**2*d2-8.4d1*a1**4*a2*d2 &
+3.3d1*a2**4+1.32d2*a1*a2**3+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2 &
+3.3d1*a1**4)
                case (-2)
                  rlYlm_laplacian =7.6147536914910937d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(5.6d1*a1**4*a2**4*d2*dx**2*dy**4 &
-4.76d2*a1**3*a2**4*dx**2*dy**4-4.76d2*a1**4*a2**3*dx**2*dy**4 &
-2.8d1*a1**3*a2**4*d2*dy**4-2.8d1*a1**4*a2**3*d2*dy**4 &
+2.1d2*a1**2*a2**4*dy**4+4.2d2*a1**3*a2**3*dy**4 &
+2.1d2*a1**4*a2**2*dy**4+5.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
-4.76d2*a1**3*a2**4*dx**4*dy**2-4.76d2*a1**4*a2**3*dx**4*dy**2 &
-4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2+3.84d2*a1**3*a2**4*d2*dx**2*dy**2 &
+3.84d2*a1**4*a2**3*d2*dx**2*dy**2+1.8d2*a1**2*a2**4*dx**2*dy**2 &
+3.6d2*a1**3*a2**3*dx**2*dy**2+1.8d2*a1**4*a2**2*dx**2*dy**2 &
+2.4d1*a1**3*a2**4*d2**2*dy**2+2.4d1*a1**4*a2**3*d2**2*dy**2 &
-1.86d2*a1**2*a2**4*d2*dy**2-3.72d2*a1**3*a2**3*d2*dy**2 &
-1.86d2*a1**4*a2**2*d2*dy**2+3.9d1*a1*a2**4*dy**2 &
+1.17d2*a1**2*a2**3*dy**2+1.17d2*a1**3*a2**2*dy**2 &
+3.9d1*a1**4*a2*dy**2-2.8d1*a1**3*a2**4*d2*dx**4 &
-2.8d1*a1**4*a2**3*d2*dx**4+2.1d2*a1**2*a2**4*dx**4 &
+4.2d2*a1**3*a2**3*dx**4+2.1d2*a1**4*a2**2*dx**4 &
+2.4d1*a1**3*a2**4*d2**2*dx**2+2.4d1*a1**4*a2**3*d2**2*dx**2 &
-1.86d2*a1**2*a2**4*d2*dx**2-3.72d2*a1**3*a2**3*d2*dx**2 &
-1.86d2*a1**4*a2**2*d2*dx**2+3.9d1*a1*a2**4*dx**2 &
+1.17d2*a1**2*a2**3*dx**2+1.17d2*a1**3*a2**2*dx**2 &
+3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2-2.4d1*a1**3*a2**3*d2**2 &
-1.2d1*a1**4*a2**2*d2**2+9.0d1*a1*a2**4*d2+2.7d2*a1**2*a2**3*d2 &
+2.7d2*a1**3*a2**2*d2+9.0d1*a1**4*a2*d2-6.6d1*a2**4-2.64d2*a1*a2**3 &
-3.96d2*a1**2*a2**2-2.64d2*a1**3*a2-6.6d1*a1**4)*dz
                case (-1)
                  rlYlm_laplacian = &
-5.3844439723186478d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(5.6d1*a1**4*a2**4*d2*dy**6-4.76d2*a1**3*a2**4*dy**6 &
-4.76d2*a1**4*a2**3*dy**6+1.12d2*a1**4*a2**4*d2*dx**2*dy**4 &
-9.52d2*a1**3*a2**4*dx**2*dy**4-9.52d2*a1**4*a2**3*dx**2*dy**4 &
-8.8d1*a1**4*a2**4*d2**2*dy**4+7.32d2*a1**3*a2**4*d2*dy**4 &
+7.32d2*a1**4*a2**3*d2*dy**4+1.2d2*a1**2*a2**4*dy**4 &
+2.4d2*a1**3*a2**3*dy**4+1.2d2*a1**4*a2**2*dy**4 &
+5.6d1*a1**4*a2**4*d2*dx**4*dy**2-4.76d2*a1**3*a2**4*dx**4*dy**2 &
-4.76d2*a1**4*a2**3*dx**4*dy**2-8.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+7.04d2*a1**3*a2**4*d2*dx**2*dy**2+7.04d2*a1**4*a2**3*d2*dx**2*dy**2 &
+3.3d2*a1**2*a2**4*dx**2*dy**2+6.6d2*a1**3*a2**3*dx**2*dy**2 &
+3.3d2*a1**4*a2**2*dx**2*dy**2+3.2d1*a1**4*a2**4*d2**3*dy**2 &
-2.28d2*a1**3*a2**4*d2**2*dy**2-2.28d2*a1**4*a2**3*d2**2*dy**2 &
-3.72d2*a1**2*a2**4*d2*dy**2-7.44d2*a1**3*a2**3*d2*dy**2 &
-3.72d2*a1**4*a2**2*d2*dy**2+2.73d2*a1*a2**4*dy**2 &
+8.19d2*a1**2*a2**3*dy**2+8.19d2*a1**3*a2**2*dy**2 &
+2.73d2*a1**4*a2*dy**2-2.8d1*a1**3*a2**4*d2*dx**4 &
-2.8d1*a1**4*a2**3*d2*dx**4+2.1d2*a1**2*a2**4*dx**4 &
+4.2d2*a1**3*a2**3*dx**4+2.1d2*a1**4*a2**2*dx**4 &
+4.4d1*a1**3*a2**4*d2**2*dx**2+4.4d1*a1**4*a2**3*d2**2*dx**2 &
-3.48d2*a1**2*a2**4*d2*dx**2-6.96d2*a1**3*a2**3*d2*dx**2 &
-3.48d2*a1**4*a2**2*d2*dx**2+1.17d2*a1*a2**4*dx**2 &
+3.51d2*a1**2*a2**3*dx**2+3.51d2*a1**3*a2**2*dx**2 &
+1.17d2*a1**4*a2*dx**2-1.6d1*a1**3*a2**4*d2**3-1.6d1*a1**4*a2**3*d2**3 &
+1.32d2*a1**2*a2**4*d2**2+2.64d2*a1**3*a2**3*d2**2 &
+1.32d2*a1**4*a2**2*d2**2-7.2d1*a1*a2**4*d2-2.16d2*a1**2*a2**3*d2 &
-2.16d2*a1**3*a2**2*d2-7.2d1*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (0)
                  rlYlm_laplacian = &
-3.40542137721830949d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+6.4d2*a1**2*a2**3*d2*dy**2 &
+6.4d2*a1**3*a2**2*d2*dy**2+3.0d2*a1*a2**3*dy**2 &
+6.0d2*a1**2*a2**2*dy**2+3.0d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.4d2*a1**2*a2**3*d2*dx**2+6.4d2*a1**3*a2**2*d2*dx**2 &
+3.0d2*a1*a2**3*dx**2+6.0d2*a1**2*a2**2*dx**2+3.0d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-7.2d1*a1**2*a2**3*d2**2 &
-7.2d1*a1**3*a2**2*d2**2-5.64d2*a1*a2**3*d2-1.128d3*a1**2*a2**2*d2 &
-5.64d2*a1**3*a2*d2+5.46d2*a2**3+1.638d3*a1*a2**2+1.638d3*a1**2*a2 &
+5.46d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian = &
-5.3844439723186478d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(5.6d1*a1**4*a2**4*d2*dx**2*dy**4-4.76d2*a1**3*a2**4*dx**2*dy**4 &
-4.76d2*a1**4*a2**3*dx**2*dy**4-2.8d1*a1**3*a2**4*d2*dy**4 &
-2.8d1*a1**4*a2**3*d2*dy**4+2.1d2*a1**2*a2**4*dy**4 &
+4.2d2*a1**3*a2**3*dy**4+2.1d2*a1**4*a2**2*dy**4 &
+1.12d2*a1**4*a2**4*d2*dx**4*dy**2-9.52d2*a1**3*a2**4*dx**4*dy**2 &
-9.52d2*a1**4*a2**3*dx**4*dy**2-8.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+7.04d2*a1**3*a2**4*d2*dx**2*dy**2+7.04d2*a1**4*a2**3*d2*dx**2*dy**2 &
+3.3d2*a1**2*a2**4*dx**2*dy**2+6.6d2*a1**3*a2**3*dx**2*dy**2 &
+3.3d2*a1**4*a2**2*dx**2*dy**2+4.4d1*a1**3*a2**4*d2**2*dy**2 &
+4.4d1*a1**4*a2**3*d2**2*dy**2-3.48d2*a1**2*a2**4*d2*dy**2 &
-6.96d2*a1**3*a2**3*d2*dy**2-3.48d2*a1**4*a2**2*d2*dy**2 &
+1.17d2*a1*a2**4*dy**2+3.51d2*a1**2*a2**3*dy**2 &
+3.51d2*a1**3*a2**2*dy**2+1.17d2*a1**4*a2*dy**2 &
+5.6d1*a1**4*a2**4*d2*dx**6-4.76d2*a1**3*a2**4*dx**6 &
-4.76d2*a1**4*a2**3*dx**6-8.8d1*a1**4*a2**4*d2**2*dx**4 &
+7.32d2*a1**3*a2**4*d2*dx**4+7.32d2*a1**4*a2**3*d2*dx**4 &
+1.2d2*a1**2*a2**4*dx**4+2.4d2*a1**3*a2**3*dx**4 &
+1.2d2*a1**4*a2**2*dx**4+3.2d1*a1**4*a2**4*d2**3*dx**2 &
-2.28d2*a1**3*a2**4*d2**2*dx**2-2.28d2*a1**4*a2**3*d2**2*dx**2 &
-3.72d2*a1**2*a2**4*d2*dx**2-7.44d2*a1**3*a2**3*d2*dx**2 &
-3.72d2*a1**4*a2**2*d2*dx**2+2.73d2*a1*a2**4*dx**2 &
+8.19d2*a1**2*a2**3*dx**2+8.19d2*a1**3*a2**2*dx**2 &
+2.73d2*a1**4*a2*dx**2-1.6d1*a1**3*a2**4*d2**3-1.6d1*a1**4*a2**3*d2**3 &
+1.32d2*a1**2*a2**4*d2**2+2.64d2*a1**3*a2**3*d2**2 &
+1.32d2*a1**4*a2**2*d2**2-7.2d1*a1*a2**4*d2-2.16d2*a1**2*a2**3*d2 &
-2.16d2*a1**3*a2**2*d2-7.2d1*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (2)
                  rlYlm_laplacian = &
-1.52295073829821874d1*E*a1**4*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(1.4d1*a1**2*a2**2*d2*dy**2 &
-1.19d2*a1*a2**2*dy**2-1.19d2*a1**2*a2*dy**2 &
+1.4d1*a1**2*a2**2*d2*dx**2-1.19d2*a1*a2**2*dx**2 &
-1.19d2*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2+1.1d2*a1*a2**2*d2 &
+1.1d2*a1**2*a2*d2-6.0d1*a2**2-1.2d2*a1*a2-6.0d1*a1**2)*dz
                case (3)
                  rlYlm_laplacian = &
-1.42458996991158946d1*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(2.4d1*a1**4*a2**4*d2*dx**2*dy**4-2.04d2*a1**3*a2**4*dx**2*dy**4 &
-2.04d2*a1**4*a2**3*dx**2*dy**4-1.2d1*a1**3*a2**4*d2*dy**4 &
-1.2d1*a1**4*a2**3*d2*dy**4+9.0d1*a1**2*a2**4*dy**4 &
+1.8d2*a1**3*a2**3*dy**4+9.0d1*a1**4*a2**2*dy**4 &
+1.6d1*a1**4*a2**4*d2*dx**4*dy**2-1.36d2*a1**3*a2**4*dx**4*dy**2 &
-1.36d2*a1**4*a2**3*dx**4*dy**2-2.4d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.92d2*a1**3*a2**4*d2*dx**2*dy**2+1.92d2*a1**4*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**4*dx**2*dy**2+1.8d2*a1**3*a2**3*dx**2*dy**2 &
+9.0d1*a1**4*a2**2*dx**2*dy**2+1.2d1*a1**3*a2**4*d2**2*dy**2 &
+1.2d1*a1**4*a2**3*d2**2*dy**2-8.4d1*a1**2*a2**4*d2*dy**2 &
-1.68d2*a1**3*a2**3*d2*dy**2-8.4d1*a1**4*a2**2*d2*dy**2 &
-3.9d1*a1*a2**4*dy**2-1.17d2*a1**2*a2**3*dy**2 &
-1.17d2*a1**3*a2**2*dy**2-3.9d1*a1**4*a2*dy**2 &
-8.0d0*a1**4*a2**4*d2*dx**6+6.8d1*a1**3*a2**4*dx**6 &
+6.8d1*a1**4*a2**3*dx**6+8.0d0*a1**4*a2**4*d2**2*dx**4 &
-8.4d1*a1**3*a2**4*d2*dx**4-8.4d1*a1**4*a2**3*d2*dx**4 &
+1.2d2*a1**2*a2**4*dx**4+2.4d2*a1**3*a2**3*dx**4 &
+1.2d2*a1**4*a2**2*dx**4+1.2d1*a1**3*a2**4*d2**2*dx**2 &
+1.2d1*a1**4*a2**3*d2**2*dx**2-8.4d1*a1**2*a2**4*d2*dx**2 &
-1.68d2*a1**3*a2**3*d2*dx**2-8.4d1*a1**4*a2**2*d2*dx**2 &
-3.9d1*a1*a2**4*dx**2-1.17d2*a1**2*a2**3*dx**2 &
-1.17d2*a1**3*a2**2*dx**2-3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2 &
-2.4d1*a1**3*a2**3*d2**2-1.2d1*a1**4*a2**2*d2**2+8.4d1*a1*a2**4*d2 &
+2.52d2*a1**2*a2**3*d2+2.52d2*a1**3*a2**2*d2+8.4d1*a1**4*a2*d2 &
-3.3d1*a2**4-1.32d2*a1*a2**3-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2 &
-3.3d1*a1**4)
                case (4)
                  rlYlm_laplacian = &
-2.01467445626964921d1*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(2.0d0*a1**3*a2**3*d2*dy**4-1.7d1*a1**2*a2**3*dy**4 &
-1.7d1*a1**3*a2**2*dy**4-1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
+1.02d2*a1**2*a2**3*dx**2*dy**2+1.02d2*a1**3*a2**2*dx**2*dy**2 &
+8.0d0*a1**2*a2**3*d2*dy**2+8.0d0*a1**3*a2**2*d2*dy**2 &
-6.0d1*a1*a2**3*dy**2-1.2d2*a1**2*a2**2*dy**2-6.0d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.7d1*a1**2*a2**3*dx**4 &
-1.7d1*a1**3*a2**2*dx**4+8.0d0*a1**2*a2**3*d2*dx**2 &
+8.0d0*a1**3*a2**2*d2*dx**2-6.0d1*a1*a2**3*dx**2 &
-1.2d2*a1**2*a2**2*dx**2-6.0d1*a1**3*a2*dx**2-1.2d1*a1*a2**3*d2 &
-2.4d1*a1**2*a2**2*d2-1.2d1*a1**3*a2*d2+7.8d1*a2**3+2.34d2*a1*a2**2 &
+2.34d2*a1**2*a2+7.8d1*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (-1)
          ! selection on l2: l1=3, m1=-1
          select case (l2)
            case (0)
              ! selection on m2: l1=3, m1=-1, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =1.43585172595163941d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dy* &
(5.0d0*(dy**2+dx**2)-4.0d0*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=3, m1=-1, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =1.24348407074185167d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*(2.0d1*a1**2*a2**2*d2*dy**4-1.1d2*a1*a2**2*dy**4 &
-1.1d2*a1**2*a2*dy**4+2.0d1*a1**2*a2**2*d2*dx**2*dy**2 &
-1.1d2*a1*a2**2*dx**2*dy**2-1.1d2*a1**2*a2*dx**2*dy**2 &
-1.6d1*a1**2*a2**2*d2**2*dy**2+7.4d1*a1*a2**2*d2*dy**2 &
+7.4d1*a1**2*a2*d2*dy**2+6.3d1*a2**2*dy**2+1.26d2*a1*a2*dy**2 &
+6.3d1*a1**2*dy**2-1.0d1*a1*a2**2*d2*dx**2-1.0d1*a1**2*a2*d2*dx**2 &
+4.5d1*a2**2*dx**2+9.0d1*a1*a2*dx**2+4.5d1*a1**2*dx**2 &
+8.0d0*a1*a2**2*d2**2+8.0d0*a1**2*a2*d2**2-3.6d1*a2**2*d2 &
-7.2d1*a1*a2*d2-3.6d1*a1**2*d2)
                case (0)
                  rlYlm_laplacian =2.48696814148370333d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+5.2d1*a1*a2**2*d2+5.2d1*a1**2*a2*d2-3.6d1*a2**2-7.2d1*a1*a2 &
-3.6d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =2.48696814148370333d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+4.2d1*a1*a2**2*d2+4.2d1*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=3, m1=-1, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =2.78051491111693767d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(2.0d1*a1**3*a2**3*d2*dy**4 &
-1.3d2*a1**2*a2**3*dy**4-1.3d2*a1**3*a2**2*dy**4 &
+2.0d1*a1**3*a2**3*d2*dx**2*dy**2-1.3d2*a1**2*a2**3*dx**2*dy**2 &
-1.3d2*a1**3*a2**2*dx**2*dy**2-1.6d1*a1**3*a2**3*d2**2*dy**2 &
+8.6d1*a1**2*a2**3*d2*dy**2+8.6d1*a1**3*a2**2*d2*dy**2 &
+9.9d1*a1*a2**3*dy**2+1.98d2*a1**2*a2**2*dy**2+9.9d1*a1**3*a2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dx**2-1.0d1*a1**3*a2**2*d2*dx**2 &
+5.5d1*a1*a2**3*dx**2+1.1d2*a1**2*a2**2*dx**2+5.5d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.2d1*a1*a2**3*d2 &
-8.4d1*a1**2*a2**2*d2-4.2d1*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)
                case (-1)
                  rlYlm_laplacian =2.78051491111693767d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(2.0d1*a1**3*a2**3*d2*dy**4 &
-1.3d2*a1**2*a2**3*dy**4-1.3d2*a1**3*a2**2*dy**4 &
+2.0d1*a1**3*a2**3*d2*dx**2*dy**2-1.3d2*a1**2*a2**3*dx**2*dy**2 &
-1.3d2*a1**3*a2**2*dx**2*dy**2-1.6d1*a1**3*a2**3*d2**2*dy**2 &
+1.06d2*a1**2*a2**3*d2*dy**2+1.06d2*a1**3*a2**2*d2*dy**2 &
-1.1d1*a1*a2**3*dy**2-2.2d1*a1**2*a2**2*dy**2-1.1d1*a1**3*a2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dx**2-1.0d1*a1**3*a2**2*d2*dx**2 &
+5.5d1*a1*a2**3*dx**2+1.1d2*a1**2*a2**2*dx**2+5.5d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.2d1*a1*a2**3*d2 &
-1.04d2*a1**2*a2**2*d2-5.2d1*a1**3*a2*d2+3.6d1*a2**3+1.08d2*a1*a2**2 &
+1.08d2*a1**2*a2+3.6d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian = &
-1.60533103241913232d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(3.0d1*a1**3*a2**3*d2*dy**4-1.95d2*a1**2*a2**3*dy**4 &
-1.95d2*a1**3*a2**2*dy**4+6.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.9d2*a1**2*a2**3*dx**2*dy**2-3.9d2*a1**3*a2**2*dx**2*dy**2 &
-4.4d1*a1**3*a2**3*d2**2*dy**2+3.04d2*a1**2*a2**3*d2*dy**2 &
+3.04d2*a1**3*a2**2*d2*dy**2-9.9d1*a1*a2**3*dy**2 &
-1.98d2*a1**2*a2**2*dy**2-9.9d1*a1**3*a2*dy**2 &
+3.0d1*a1**3*a2**3*d2*dx**4-1.95d2*a1**2*a2**3*dx**4 &
-1.95d2*a1**3*a2**2*dx**4-4.4d1*a1**3*a2**3*d2**2*dx**2 &
+3.04d2*a1**2*a2**3*d2*dx**2+3.04d2*a1**3*a2**2*d2*dx**2 &
-9.9d1*a1*a2**3*dx**2-1.98d2*a1**2*a2**2*dx**2-9.9d1*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.28d2*a1**2*a2**3*d2**2 &
-1.28d2*a1**3*a2**2*d2**2+1.44d2*a1*a2**3*d2+2.88d2*a1**2*a2**2*d2 &
+1.44d2*a1**3*a2*d2-5.4d1*a2**3-1.62d2*a1*a2**2-1.62d2*a1**2*a2 &
-5.4d1*a1**3)
                case (1)
                  rlYlm_laplacian =5.56102982223387534d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-6.5d1*a1*a2**2*dy**2-6.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-6.5d1*a1*a2**2*dx**2-6.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+5.8d1*a1*a2**2*d2+5.8d1*a1**2*a2*d2-3.3d1*a2**2-6.6d1*a1*a2 &
-3.3d1*a1**2)*dz
                case (2)
                  rlYlm_laplacian = &
-2.78051491111693767d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(1.0d1*a1**3*a2**3*d2*dy**4-6.5d1*a1**2*a2**3*dy**4 &
-6.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+3.8d1*a1**2*a2**3*d2*dy**2+3.8d1*a1**3*a2**2*d2*dy**2 &
+7.7d1*a1*a2**3*dy**2+1.54d2*a1**2*a2**2*dy**2+7.7d1*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+6.5d1*a1**2*a2**3*dx**4 &
+6.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-5.8d1*a1**2*a2**3*d2*dx**2-5.8d1*a1**3*a2**2*d2*dx**2 &
+3.3d1*a1*a2**3*dx**2+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.2d1*a1*a2**3*d2 &
-8.4d1*a1**2*a2**2*d2-4.2d1*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=3, m1=-1, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-1.50164967891712101d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)* &
(2.0d1*a1**3*a2**3*d2*dy**6-1.5d2*a1**2*a2**3*dy**6 &
-1.5d2*a1**3*a2**2*dy**6-4.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+3.0d2*a1**2*a2**3*dx**2*dy**4+3.0d2*a1**3*a2**2*dx**2*dy**4 &
-1.6d1*a1**3*a2**3*d2**2*dy**4+7.8d1*a1**2*a2**3*d2*dy**4 &
+7.8d1*a1**3*a2**2*d2*dy**4+2.73d2*a1*a2**3*dy**4 &
+5.46d2*a1**2*a2**2*dy**4+2.73d2*a1**3*a2*dy**4 &
-6.0d1*a1**3*a2**3*d2*dx**4*dy**2+4.5d2*a1**2*a2**3*dx**4*dy**2 &
+4.5d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-3.24d2*a1**2*a2**3*d2*dx**2*dy**2-3.24d2*a1**3*a2**2*d2*dx**2*dy**2 &
-2.34d2*a1*a2**3*dx**2*dy**2-4.68d2*a1**2*a2**2*dx**2*dy**2 &
-2.34d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.44d2*a1*a2**3*d2*dy**2 &
-2.88d2*a1**2*a2**2*d2*dy**2-1.44d2*a1**3*a2*d2*dy**2 &
-6.6d1*a2**3*dy**2-1.98d2*a1*a2**2*dy**2-1.98d2*a1**2*a2*dy**2 &
-6.6d1*a1**3*dy**2+3.0d1*a1**2*a2**3*d2*dx**4 &
+3.0d1*a1**3*a2**2*d2*dx**4-1.95d2*a1*a2**3*dx**4 &
-3.9d2*a1**2*a2**2*dx**4-1.95d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.44d2*a1*a2**3*d2*dx**2+2.88d2*a1**2*a2**2*d2*dx**2 &
+1.44d2*a1**3*a2*d2*dx**2+6.6d1*a2**3*dx**2+1.98d2*a1*a2**2*dx**2 &
+1.98d2*a1**2*a2*dx**2+6.6d1*a1**3*dx**2)
                case (-2)
                  rlYlm_laplacian =7.35655097152228142d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(2.0d1*a1**3*a2**3*d2*dy**4 &
-1.5d2*a1**2*a2**3*dy**4-1.5d2*a1**3*a2**2*dy**4 &
+2.0d1*a1**3*a2**3*d2*dx**2*dy**2-1.5d2*a1**2*a2**3*dx**2*dy**2 &
-1.5d2*a1**3*a2**2*dx**2*dy**2-1.6d1*a1**3*a2**3*d2**2*dy**2 &
+1.18d2*a1**2*a2**3*d2*dy**2+1.18d2*a1**3*a2**2*d2*dy**2 &
+1.3d1*a1*a2**3*dy**2+2.6d1*a1**2*a2**2*dy**2+1.3d1*a1**3*a2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dx**2-1.0d1*a1**3*a2**2*d2*dx**2 &
+6.5d1*a1*a2**3*dx**2+1.3d2*a1**2*a2**2*dx**2+6.5d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.8d1*a1*a2**3*d2 &
-1.16d2*a1**2*a2**2*d2-5.8d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =-1.16317283965674489d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(1.0d2*a1**4*a2**4*d2*dy**6 &
-7.5d2*a1**3*a2**4*dy**6-7.5d2*a1**4*a2**3*dy**6 &
+2.0d2*a1**4*a2**4*d2*dx**2*dy**4-1.5d3*a1**3*a2**4*dx**2*dy**4 &
-1.5d3*a1**4*a2**3*dx**2*dy**4-1.6d2*a1**4*a2**4*d2**2*dy**4 &
+1.23d3*a1**3*a2**4*d2*dy**4+1.23d3*a1**4*a2**3*d2*dy**4 &
-1.95d2*a1**2*a2**4*dy**4-3.9d2*a1**3*a2**3*dy**4 &
-1.95d2*a1**4*a2**2*dy**4+1.0d2*a1**4*a2**4*d2*dx**4*dy**2 &
-7.5d2*a1**3*a2**4*dx**4*dy**2-7.5d2*a1**4*a2**3*dx**4*dy**2 &
-1.6d2*a1**4*a2**4*d2**2*dx**2*dy**2+1.18d3*a1**3*a2**4*d2*dx**2*dy**2 &
+1.18d3*a1**4*a2**3*d2*dx**2*dy**2+1.3d2*a1**2*a2**4*dx**2*dy**2 &
+2.6d2*a1**3*a2**3*dx**2*dy**2+1.3d2*a1**4*a2**2*dx**2*dy**2 &
+6.4d1*a1**4*a2**4*d2**3*dy**2-4.96d2*a1**3*a2**4*d2**2*dy**2 &
-4.96d2*a1**4*a2**3*d2**2*dy**2+9.2d1*a1**2*a2**4*d2*dy**2 &
+1.84d2*a1**3*a2**3*d2*dy**2+9.2d1*a1**4*a2**2*d2*dy**2 &
+6.6d1*a1*a2**4*dy**2+1.98d2*a1**2*a2**3*dy**2 &
+1.98d2*a1**3*a2**2*dy**2+6.6d1*a1**4*a2*dy**2 &
-5.0d1*a1**3*a2**4*d2*dx**4-5.0d1*a1**4*a2**3*d2*dx**4 &
+3.25d2*a1**2*a2**4*dx**4+6.5d2*a1**3*a2**3*dx**4 &
+3.25d2*a1**4*a2**2*dx**4+8.0d1*a1**3*a2**4*d2**2*dx**2 &
+8.0d1*a1**4*a2**3*d2**2*dx**2-5.8d2*a1**2*a2**4*d2*dx**2 &
-1.16d3*a1**3*a2**3*d2*dx**2-5.8d2*a1**4*a2**2*d2*dx**2 &
+3.3d2*a1*a2**4*dx**2+9.9d2*a1**2*a2**3*dx**2+9.9d2*a1**3*a2**2*dx**2 &
+3.3d2*a1**4*a2*dx**2-3.2d1*a1**3*a2**4*d2**3-3.2d1*a1**4*a2**3*d2**3 &
+2.72d2*a1**2*a2**4*d2**2+5.44d2*a1**3*a2**3*d2**2 &
+2.72d2*a1**4*a2**2*d2**2-3.72d2*a1*a2**4*d2-1.116d3*a1**2*a2**3*d2 &
-1.116d3*a1**3*a2**2*d2-3.72d2*a1**4*a2*d2+9.0d1*a2**4+3.6d2*a1*a2**3 &
+5.4d2*a1**2*a2**2+3.6d2*a1**3*a2+9.0d1*a1**4)
                case (0)
                  rlYlm_laplacian = &
-1.89945329321545261d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(5.0d1*a1**3*a2**3*d2*dy**4-3.75d2*a1**2*a2**3*dy**4 &
-3.75d2*a1**3*a2**2*dy**4+1.0d2*a1**3*a2**3*d2*dx**2*dy**2 &
-7.5d2*a1**2*a2**3*dx**2*dy**2-7.5d2*a1**3*a2**2*dx**2*dy**2 &
-6.0d1*a1**3*a2**3*d2**2*dy**2+4.8d2*a1**2*a2**3*d2*dy**2 &
+4.8d2*a1**3*a2**2*d2*dy**2-1.95d2*a1*a2**3*dy**2 &
-3.9d2*a1**2*a2**2*dy**2-1.95d2*a1**3*a2*dy**2 &
+5.0d1*a1**3*a2**3*d2*dx**4-3.75d2*a1**2*a2**3*dx**4 &
-3.75d2*a1**3*a2**2*dx**4-6.0d1*a1**3*a2**3*d2**2*dx**2 &
+4.8d2*a1**2*a2**3*d2*dx**2+4.8d2*a1**3*a2**2*d2*dx**2 &
-1.95d2*a1*a2**3*dx**2-3.9d2*a1**2*a2**2*dx**2-1.95d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.44d2*a1**2*a2**3*d2**2 &
-1.44d2*a1**3*a2**2*d2**2+1.68d2*a1*a2**3*d2+3.36d2*a1**2*a2**2*d2 &
+1.68d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2-1.98d2*a1**2*a2 &
-6.6d1*a1**3)*dz
                case (1)
                  rlYlm_laplacian = &
-2.32634567931348979d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(5.0d1*a1**3*a2**3*d2*dy**4-3.75d2*a1**2*a2**3*dy**4 &
-3.75d2*a1**3*a2**2*dy**4+1.0d2*a1**3*a2**3*d2*dx**2*dy**2 &
-7.5d2*a1**2*a2**3*dx**2*dy**2-7.5d2*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+6.4d2*a1**2*a2**3*d2*dy**2 &
+6.4d2*a1**3*a2**2*d2*dy**2-2.6d2*a1*a2**3*dy**2 &
-5.2d2*a1**2*a2**2*dy**2-2.6d2*a1**3*a2*dy**2 &
+5.0d1*a1**3*a2**3*d2*dx**4-3.75d2*a1**2*a2**3*dx**4 &
-3.75d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.4d2*a1**2*a2**3*d2*dx**2+6.4d2*a1**3*a2**2*d2*dx**2 &
-2.6d2*a1*a2**3*dx**2-5.2d2*a1**2*a2**2*dx**2-2.6d2*a1**3*a2*dx**2 &
+3.2d1*a1**3*a2**3*d2**3-2.88d2*a1**2*a2**3*d2**2 &
-2.88d2*a1**3*a2**2*d2**2+3.36d2*a1*a2**3*d2+6.72d2*a1**2*a2**2*d2 &
+3.36d2*a1**3*a2*d2-1.32d2*a2**3-3.96d2*a1*a2**2-3.96d2*a1**2*a2 &
-1.32d2*a1**3)
                case (2)
                  rlYlm_laplacian = &
-7.35655097152228142d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(1.0d1*a1**3*a2**3*d2*dy**4-7.5d1*a1**2*a2**3*dy**4 &
-7.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+5.4d1*a1**2*a2**3*d2*dy**2+5.4d1*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+7.5d1*a1**2*a2**3*dx**4 &
+7.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-7.4d1*a1**2*a2**3*d2*dx**2-7.4d1*a1**3*a2**2*d2*dx**2 &
+9.1d1*a1*a2**3*dx**2+1.82d2*a1**2*a2**2*dx**2+9.1d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.8d1*a1*a2**3*d2 &
-1.16d2*a1**2*a2**2*d2-5.8d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (3)
                  rlYlm_laplacian = &
-3.00329935783424201d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(3.0d1*a1**3*a2**3*d2*dy**4-2.25d2*a1**2*a2**3*dy**4 &
-2.25d2*a1**3*a2**2*dy**4+2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.5d2*a1**2*a2**3*dx**2*dy**2-1.5d2*a1**3*a2**2*dx**2*dy**2 &
-2.4d1*a1**3*a2**3*d2**2*dy**2+1.32d2*a1**2*a2**3*d2*dy**2 &
+1.32d2*a1**3*a2**2*d2*dy**2+3.12d2*a1*a2**3*dy**2 &
+6.24d2*a1**2*a2**2*dy**2+3.12d2*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+7.5d1*a1**2*a2**3*dx**4 &
+7.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-8.4d1*a1**2*a2**3*d2*dx**2-8.4d1*a1**3*a2**2*d2*dx**2 &
+1.56d2*a1*a2**3*dx**2+3.12d2*a1**2*a2**2*dx**2+1.56d2*a1**3*a2*dx**2 &
+2.4d1*a1**2*a2**3*d2**2+2.4d1*a1**3*a2**2*d2**2-1.44d2*a1*a2**3*d2 &
-2.88d2*a1**2*a2**2*d2-1.44d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2 &
-1.98d2*a1**2*a2-6.6d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=3, m1=-1, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-6.37096002557338818d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(2.0d1*a1**3*a2**3*d2*dy**6-1.7d2*a1**2*a2**3*dy**6 &
-1.7d2*a1**3*a2**2*dy**6-1.6d1*a1**3*a2**3*d2**2*dy**4 &
+9.0d1*a1**2*a2**3*d2*dy**4+9.0d1*a1**3*a2**2*d2*dy**4 &
+3.45d2*a1*a2**3*dy**4+6.9d2*a1**2*a2**2*dy**4+3.45d2*a1**3*a2*dy**4 &
-2.0d1*a1**3*a2**3*d2*dx**4*dy**2+1.7d2*a1**2*a2**3*dx**4*dy**2 &
+1.7d2*a1**3*a2**2*dx**4*dy**2+1.6d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-1.4d2*a1**2*a2**3*d2*dx**2*dy**2-1.4d2*a1**3*a2**2*d2*dx**2*dy**2 &
+3.0d1*a1*a2**3*dx**2*dy**2+6.0d1*a1**2*a2**2*dx**2*dy**2 &
+3.0d1*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.62d2*a1*a2**3*d2*dy**2 &
-3.24d2*a1**2*a2**2*d2*dy**2-1.62d2*a1**3*a2*d2*dy**2 &
-1.17d2*a2**3*dy**2-3.51d2*a1*a2**2*dy**2-3.51d2*a1**2*a2*dy**2 &
-1.17d2*a1**3*dy**2+1.0d1*a1**2*a2**3*d2*dx**4 &
+1.0d1*a1**3*a2**2*d2*dx**4-7.5d1*a1*a2**3*dx**4 &
-1.5d2*a1**2*a2**2*dx**4-7.5d1*a1**3*a2*dx**4 &
-8.0d0*a1**2*a2**3*d2**2*dx**2-8.0d0*a1**3*a2**2*d2**2*dx**2 &
+5.4d1*a1*a2**3*d2*dx**2+1.08d2*a1**2*a2**2*d2*dx**2 &
+5.4d1*a1**3*a2*d2*dx**2+3.9d1*a2**3*dx**2+1.17d2*a1*a2**2*dx**2 &
+1.17d2*a1**2*a2*dx**2+3.9d1*a1**3*dx**2)
                case (-3)
                  rlYlm_laplacian = &
-4.50494903675136302d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(2.0d1*a1**3*a2**3*d2*dy**6-1.7d2*a1**2*a2**3*dy**6 &
-1.7d2*a1**3*a2**2*dy**6-4.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+3.4d2*a1**2*a2**3*dx**2*dy**4+3.4d2*a1**3*a2**2*dx**2*dy**4 &
-1.6d1*a1**3*a2**3*d2**2*dy**4+1.1d2*a1**2*a2**3*d2*dy**4 &
+1.1d2*a1**3*a2**2*d2*dy**4+1.95d2*a1*a2**3*dy**4 &
+3.9d2*a1**2*a2**2*dy**4+1.95d2*a1**3*a2*dy**4 &
-6.0d1*a1**3*a2**3*d2*dx**4*dy**2+5.1d2*a1**2*a2**3*dx**4*dy**2 &
+5.1d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-4.2d2*a1**2*a2**3*d2*dx**2*dy**2-4.2d2*a1**3*a2**2*d2*dx**2*dy**2 &
+9.0d1*a1*a2**3*dx**2*dy**2+1.8d2*a1**2*a2**2*dx**2*dy**2 &
+9.0d1*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.92d2*a1*a2**3*d2*dy**2 &
-3.84d2*a1**2*a2**2*d2*dy**2-1.92d2*a1**3*a2*d2*dy**2 &
+7.8d1*a2**3*dy**2+2.34d2*a1*a2**2*dy**2+2.34d2*a1**2*a2*dy**2 &
+7.8d1*a1**3*dy**2+3.0d1*a1**2*a2**3*d2*dx**4 &
+3.0d1*a1**3*a2**2*d2*dx**4-2.25d2*a1*a2**3*dx**4 &
-4.5d2*a1**2*a2**2*dx**4-2.25d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.92d2*a1*a2**3*d2*dx**2+3.84d2*a1**2*a2**2*d2*dx**2 &
+1.92d2*a1**3*a2*d2*dx**2-7.8d1*a2**3*dx**2-2.34d2*a1*a2**2*dx**2 &
-2.34d2*a1**2*a2*dx**2-7.8d1*a1**3*dx**2)*dz
                case (-2)
                  rlYlm_laplacian = &
-2.40799654862869848d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(1.4d2*a1**4*a2**4*d2*dy**6-1.19d3*a1**3*a2**4*dy**6 &
-1.19d3*a1**4*a2**3*dy**6+2.8d2*a1**4*a2**4*d2*dx**2*dy**4 &
-2.38d3*a1**3*a2**4*dx**2*dy**4-2.38d3*a1**4*a2**3*dx**2*dy**4 &
-2.32d2*a1**4*a2**4*d2**2*dy**4+2.01d3*a1**3*a2**4*d2*dy**4 &
+2.01d3*a1**4*a2**3*d2*dy**4-2.85d2*a1**2*a2**4*dy**4 &
-5.7d2*a1**3*a2**3*dy**4-2.85d2*a1**4*a2**2*dy**4 &
+1.4d2*a1**4*a2**4*d2*dx**4*dy**2-1.19d3*a1**3*a2**4*dx**4*dy**2 &
-1.19d3*a1**4*a2**3*dx**4*dy**2-2.32d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.94d3*a1**3*a2**4*d2*dx**2*dy**2+1.94d3*a1**4*a2**3*d2*dx**2*dy**2 &
+2.4d2*a1**2*a2**4*dx**2*dy**2+4.8d2*a1**3*a2**3*dx**2*dy**2 &
+2.4d2*a1**4*a2**2*dx**2*dy**2+9.6d1*a1**4*a2**4*d2**3*dy**2 &
-8.28d2*a1**3*a2**4*d2**2*dy**2-8.28d2*a1**4*a2**3*d2**2*dy**2 &
+8.4d1*a1**2*a2**4*d2*dy**2+1.68d2*a1**3*a2**3*d2*dy**2 &
+8.4d1*a1**4*a2**2*d2*dy**2+3.9d1*a1*a2**4*dy**2 &
+1.17d2*a1**2*a2**3*dy**2+1.17d2*a1**3*a2**2*dy**2 &
+3.9d1*a1**4*a2*dy**2-7.0d1*a1**3*a2**4*d2*dx**4 &
-7.0d1*a1**4*a2**3*d2*dx**4+5.25d2*a1**2*a2**4*dx**4 &
+1.05d3*a1**3*a2**3*dx**4+5.25d2*a1**4*a2**2*dx**4 &
+1.16d2*a1**3*a2**4*d2**2*dx**2+1.16d2*a1**4*a2**3*d2**2*dx**2 &
-9.48d2*a1**2*a2**4*d2*dx**2-1.896d3*a1**3*a2**3*d2*dx**2 &
-9.48d2*a1**4*a2**2*d2*dx**2+5.07d2*a1*a2**4*dx**2 &
+1.521d3*a1**2*a2**3*dx**2+1.521d3*a1**3*a2**2*dx**2 &
+5.07d2*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+4.44d2*a1**2*a2**4*d2**2+8.88d2*a1**3*a2**3*d2**2 &
+4.44d2*a1**4*a2**2*d2**2-5.76d2*a1*a2**4*d2-1.728d3*a1**2*a2**3*d2 &
-1.728d3*a1**3*a2**2*d2-5.76d2*a1**4*a2*d2+1.65d2*a2**4+6.6d2*a1*a2**3 &
+9.9d2*a1**2*a2**2+6.6d2*a1**3*a2+1.65d2*a1**4)
                case (-1)
                  rlYlm_laplacian = &
-1.70271068860915474d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(1.4d2*a1**4*a2**4*d2*dy**6-1.19d3*a1**3*a2**4*dy**6 &
-1.19d3*a1**4*a2**3*dy**6+2.8d2*a1**4*a2**4*d2*dx**2*dy**4 &
-2.38d3*a1**3*a2**4*dx**2*dy**4-2.38d3*a1**4*a2**3*dx**2*dy**4 &
-1.92d2*a1**4*a2**4*d2**2*dy**4+1.69d3*a1**3*a2**4*d2*dy**4 &
+1.69d3*a1**4*a2**3*d2*dy**4-4.35d2*a1**2*a2**4*dy**4 &
-8.7d2*a1**3*a2**3*dy**4-4.35d2*a1**4*a2**2*dy**4 &
+1.4d2*a1**4*a2**4*d2*dx**4*dy**2-1.19d3*a1**3*a2**4*dx**4*dy**2 &
-1.19d3*a1**4*a2**3*dx**4*dy**2-1.92d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.62d3*a1**3*a2**4*d2*dx**2*dy**2+1.62d3*a1**4*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**4*dx**2*dy**2+1.8d2*a1**3*a2**3*dx**2*dy**2 &
+9.0d1*a1**4*a2**2*dx**2*dy**2+6.4d1*a1**4*a2**4*d2**3*dy**2 &
-5.76d2*a1**3*a2**4*d2**2*dy**2-5.76d2*a1**4*a2**3*d2**2*dy**2 &
+2.28d2*a1**2*a2**4*d2*dy**2+4.56d2*a1**3*a2**3*d2*dy**2 &
+2.28d2*a1**4*a2**2*d2*dy**2+7.8d1*a1*a2**4*dy**2 &
+2.34d2*a1**2*a2**3*dy**2+2.34d2*a1**3*a2**2*dy**2 &
+7.8d1*a1**4*a2*dy**2-7.0d1*a1**3*a2**4*d2*dx**4 &
-7.0d1*a1**4*a2**3*d2*dx**4+5.25d2*a1**2*a2**4*dx**4 &
+1.05d3*a1**3*a2**3*dx**4+5.25d2*a1**4*a2**2*dx**4 &
+9.6d1*a1**3*a2**4*d2**2*dx**2+9.6d1*a1**4*a2**3*d2**2*dx**2 &
-8.28d2*a1**2*a2**4*d2*dx**2-1.656d3*a1**3*a2**3*d2*dx**2 &
-8.28d2*a1**4*a2**2*d2*dx**2+7.02d2*a1*a2**4*dx**2 &
+2.106d3*a1**2*a2**3*dx**2+2.106d3*a1**3*a2**2*dx**2 &
+7.02d2*a1**4*a2*dx**2-3.2d1*a1**3*a2**4*d2**3-3.2d1*a1**4*a2**3*d2**3 &
+3.36d2*a1**2*a2**4*d2**2+6.72d2*a1**3*a2**3*d2**2 &
+3.36d2*a1**4*a2**2*d2**2-6.84d2*a1*a2**4*d2-2.052d3*a1**2*a2**3*d2 &
-2.052d3*a1**3*a2**2*d2-6.84d2*a1**4*a2*d2+3.3d2*a2**4+1.32d3*a1*a2**3 &
+1.98d3*a1**2*a2**2+1.32d3*a1**3*a2+3.3d2*a1**4)*dz
                case (0)
                  rlYlm_laplacian =5.3844439723186478d-1*E*a1**2*a2*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*dy* &
(3.5d2*a1**4*a2**4*d2*dy**6-2.975d3*a1**3*a2**4*dy**6 &
-2.975d3*a1**4*a2**3*dy**6+1.05d3*a1**4*a2**4*d2*dx**2*dy**4 &
-8.925d3*a1**3*a2**4*dx**2*dy**4-8.925d3*a1**4*a2**3*dx**2*dy**4 &
-6.8d2*a1**4*a2**4*d2**2*dy**4+6.0d3*a1**3*a2**4*d2*dy**4 &
+6.0d3*a1**4*a2**3*d2*dy**4-1.65d3*a1**2*a2**4*dy**4 &
-3.3d3*a1**3*a2**3*dy**4-1.65d3*a1**4*a2**2*dy**4 &
+1.05d3*a1**4*a2**4*d2*dx**4*dy**2-8.925d3*a1**3*a2**4*dx**4*dy**2 &
-8.925d3*a1**4*a2**3*dx**4*dy**2-1.36d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.2d4*a1**3*a2**4*d2*dx**2*dy**2+1.2d4*a1**4*a2**3*d2*dx**2*dy**2 &
-3.3d3*a1**2*a2**4*dx**2*dy**2-6.6d3*a1**3*a2**3*dx**2*dy**2 &
-3.3d3*a1**4*a2**2*dx**2*dy**2+4.0d2*a1**4*a2**4*d2**3*dy**2 &
-3.72d3*a1**3*a2**4*d2**2*dy**2-3.72d3*a1**4*a2**3*d2**2*dy**2 &
+2.46d3*a1**2*a2**4*d2*dy**2+4.92d3*a1**3*a2**3*d2*dy**2 &
+2.46d3*a1**4*a2**2*d2*dy**2-3.9d2*a1*a2**4*dy**2 &
-1.17d3*a1**2*a2**3*dy**2-1.17d3*a1**3*a2**2*dy**2 &
-3.9d2*a1**4*a2*dy**2+3.5d2*a1**4*a2**4*d2*dx**6 &
-2.975d3*a1**3*a2**4*dx**6-2.975d3*a1**4*a2**3*dx**6 &
-6.8d2*a1**4*a2**4*d2**2*dx**4+6.0d3*a1**3*a2**4*d2*dx**4 &
+6.0d3*a1**4*a2**3*d2*dx**4-1.65d3*a1**2*a2**4*dx**4 &
-3.3d3*a1**3*a2**3*dx**4-1.65d3*a1**4*a2**2*dx**4 &
+4.0d2*a1**4*a2**4*d2**3*dx**2-3.72d3*a1**3*a2**4*d2**2*dx**2 &
-3.72d3*a1**4*a2**3*d2**2*dx**2+2.46d3*a1**2*a2**4*d2*dx**2 &
+4.92d3*a1**3*a2**3*d2*dx**2+2.46d3*a1**4*a2**2*d2*dx**2 &
-3.9d2*a1*a2**4*dx**2-1.17d3*a1**2*a2**3*dx**2 &
-1.17d3*a1**3*a2**2*dx**2-3.9d2*a1**4*a2*dx**2-6.4d1*a1**4*a2**4*d2**4 &
+6.08d2*a1**3*a2**4*d2**3+6.08d2*a1**4*a2**3*d2**3 &
-3.84d2*a1**2*a2**4*d2**2-7.68d2*a1**3*a2**3*d2**2 &
-3.84d2*a1**4*a2**2*d2**2-7.44d2*a1*a2**4*d2-2.232d3*a1**2*a2**3*d2 &
-2.232d3*a1**3*a2**2*d2-7.44d2*a1**4*a2*d2+6.6d2*a2**4+2.64d3*a1*a2**3 &
+3.96d3*a1**2*a2**2+2.64d3*a1**3*a2+6.6d2*a1**4)
                case (1)
                  rlYlm_laplacian = &
-3.40542137721830949d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-9.6d1*a1**3*a2**3*d2**2*dy**2+8.8d2*a1**2*a2**3*d2*dy**2 &
+8.8d2*a1**3*a2**2*d2*dy**2-4.8d2*a1*a2**3*dy**2 &
-9.6d2*a1**2*a2**2*dy**2-4.8d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-9.6d1*a1**3*a2**3*d2**2*dx**2 &
+8.8d2*a1**2*a2**3*d2*dx**2+8.8d2*a1**3*a2**2*d2*dx**2 &
-4.8d2*a1*a2**3*dx**2-9.6d2*a1**2*a2**2*dx**2-4.8d2*a1**3*a2*dx**2 &
+3.2d1*a1**3*a2**3*d2**3-3.36d2*a1**2*a2**3*d2**2 &
-3.36d2*a1**3*a2**2*d2**2+5.28d2*a1*a2**3*d2+1.056d3*a1**2*a2**2*d2 &
+5.28d2*a1**3*a2*d2-3.12d2*a2**3-9.36d2*a1*a2**2-9.36d2*a1**2*a2 &
-3.12d2*a1**3)*dz
                case (2)
                  rlYlm_laplacian =2.40799654862869848d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(7.0d1*a1**4*a2**4*d2*dy**6 &
-5.95d2*a1**3*a2**4*dy**6-5.95d2*a1**4*a2**3*dy**6 &
+7.0d1*a1**4*a2**4*d2*dx**2*dy**4-5.95d2*a1**3*a2**4*dx**2*dy**4 &
-5.95d2*a1**4*a2**3*dx**2*dy**4-1.16d2*a1**4*a2**4*d2**2*dy**4 &
+9.7d2*a1**3*a2**4*d2*dy**4+9.7d2*a1**4*a2**3*d2*dy**4 &
+1.2d2*a1**2*a2**4*dy**4+2.4d2*a1**3*a2**3*dy**4 &
+1.2d2*a1**4*a2**2*dy**4-7.0d1*a1**4*a2**4*d2*dx**4*dy**2 &
+5.95d2*a1**3*a2**4*dx**4*dy**2+5.95d2*a1**4*a2**3*dx**4*dy**2 &
-1.4d2*a1**3*a2**4*d2*dx**2*dy**2-1.4d2*a1**4*a2**3*d2*dx**2*dy**2 &
+1.05d3*a1**2*a2**4*dx**2*dy**2+2.1d3*a1**3*a2**3*dx**2*dy**2 &
+1.05d3*a1**4*a2**2*dx**2*dy**2+4.8d1*a1**4*a2**4*d2**3*dy**2 &
-3.56d2*a1**3*a2**4*d2**2*dy**2-3.56d2*a1**4*a2**3*d2**2*dy**2 &
-4.32d2*a1**2*a2**4*d2*dy**2-8.64d2*a1**3*a2**3*d2*dy**2 &
-4.32d2*a1**4*a2**2*d2*dy**2+2.73d2*a1*a2**4*dy**2 &
+8.19d2*a1**2*a2**3*dy**2+8.19d2*a1**3*a2**2*dy**2 &
+2.73d2*a1**4*a2*dy**2-7.0d1*a1**4*a2**4*d2*dx**6 &
+5.95d2*a1**3*a2**4*dx**6+5.95d2*a1**4*a2**3*dx**6 &
+1.16d2*a1**4*a2**4*d2**2*dx**4-1.11d3*a1**3*a2**4*d2*dx**4 &
-1.11d3*a1**4*a2**3*d2*dx**4+9.3d2*a1**2*a2**4*dx**4 &
+1.86d3*a1**3*a2**3*dx**4+9.3d2*a1**4*a2**2*dx**4 &
-4.8d1*a1**4*a2**4*d2**3*dx**2+5.88d2*a1**3*a2**4*d2**2*dx**2 &
+5.88d2*a1**4*a2**3*d2**2*dx**2-1.464d3*a1**2*a2**4*d2*dx**2 &
-2.928d3*a1**3*a2**3*d2*dx**2-1.464d3*a1**4*a2**2*d2*dx**2 &
+7.41d2*a1*a2**4*dx**2+2.223d3*a1**2*a2**3*dx**2 &
+2.223d3*a1**3*a2**2*dx**2+7.41d2*a1**4*a2*dx**2 &
-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+4.44d2*a1**2*a2**4*d2**2+8.88d2*a1**3*a2**3*d2**2 &
+4.44d2*a1**4*a2**2*d2**2-5.76d2*a1*a2**4*d2-1.728d3*a1**2*a2**3*d2 &
-1.728d3*a1**3*a2**2*d2-5.76d2*a1**4*a2*d2+1.65d2*a2**4+6.6d2*a1*a2**3 &
+9.9d2*a1**2*a2**2+6.6d2*a1**3*a2+1.65d2*a1**4)
                case (3)
                  rlYlm_laplacian = &
-9.00989807350272604d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(3.0d1*a1**3*a2**3*d2*dy**4-2.55d2*a1**2*a2**3*dy**4 &
-2.55d2*a1**3*a2**2*dy**4+2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.7d2*a1**2*a2**3*dx**2*dy**2-1.7d2*a1**3*a2**2*dx**2*dy**2 &
-2.4d1*a1**3*a2**3*d2**2*dy**2+1.8d2*a1**2*a2**3*d2*dy**2 &
+1.8d2*a1**3*a2**2*d2*dy**2+1.8d2*a1*a2**3*dy**2 &
+3.6d2*a1**2*a2**2*dy**2+1.8d2*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+8.5d1*a1**2*a2**3*dx**4 &
+8.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-1.0d2*a1**2*a2**3*d2*dx**2-1.0d2*a1**3*a2**2*d2*dx**2 &
+2.4d2*a1*a2**3*dx**2+4.8d2*a1**2*a2**2*dx**2+2.4d2*a1**3*a2*dx**2 &
+2.4d1*a1**2*a2**3*d2**2+2.4d1*a1**3*a2**2*d2**2-1.92d2*a1*a2**3*d2 &
-3.84d2*a1**2*a2**2*d2-1.92d2*a1**3*a2*d2+7.8d1*a2**3+2.34d2*a1*a2**2 &
+2.34d2*a1**2*a2+7.8d1*a1**3)*dz
                case (4)
                  rlYlm_laplacian =3.18548001278669409d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(1.0d1*a1**3*a2**3*d2*dy**6 &
-8.5d1*a1**2*a2**3*dy**6-8.5d1*a1**3*a2**2*dy**6 &
-5.0d1*a1**3*a2**3*d2*dx**2*dy**4+4.25d2*a1**2*a2**3*dx**2*dy**4 &
+4.25d2*a1**3*a2**2*dx**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**4 &
+4.0d1*a1**2*a2**3*d2*dy**4+4.0d1*a1**3*a2**2*d2*dy**4 &
+2.1d2*a1*a2**3*dy**4+4.2d2*a1**2*a2**2*dy**4+2.1d2*a1**3*a2*dy**4 &
-5.0d1*a1**3*a2**3*d2*dx**4*dy**2+4.25d2*a1**2*a2**3*dx**4*dy**2 &
+4.25d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-3.2d2*a1**2*a2**3*d2*dx**2*dy**2-3.2d2*a1**3*a2**2*d2*dx**2*dy**2 &
-6.6d2*a1*a2**3*dx**2*dy**2-1.32d3*a1**2*a2**2*dx**2*dy**2 &
-6.6d2*a1**3*a2*dx**2*dy**2+1.6d1*a1**2*a2**3*d2**2*dy**2 &
+1.6d1*a1**3*a2**2*d2**2*dy**2-1.08d2*a1*a2**3*d2*dy**2 &
-2.16d2*a1**2*a2**2*d2*dy**2-1.08d2*a1**3*a2*d2*dy**2 &
-7.8d1*a2**3*dy**2-2.34d2*a1*a2**2*dy**2-2.34d2*a1**2*a2*dy**2 &
-7.8d1*a1**3*dy**2+1.0d1*a1**3*a2**3*d2*dx**6-8.5d1*a1**2*a2**3*dx**6 &
-8.5d1*a1**3*a2**2*dx**6-8.0d0*a1**3*a2**3*d2**2*dx**4 &
+1.2d2*a1**2*a2**3*d2*dx**4+1.2d2*a1**3*a2**2*d2*dx**4 &
-3.9d2*a1*a2**3*dx**4-7.8d2*a1**2*a2**2*dx**4-3.9d2*a1**3*a2*dx**4 &
-4.8d1*a1**2*a2**3*d2**2*dx**2-4.8d1*a1**3*a2**2*d2**2*dx**2 &
+3.24d2*a1*a2**3*d2*dx**2+6.48d2*a1**2*a2**2*d2*dx**2 &
+3.24d2*a1**3*a2*d2*dx**2+2.34d2*a2**3*dx**2+7.02d2*a1*a2**2*dx**2 &
+7.02d2*a1**2*a2*dx**2+2.34d2*a1**3*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (0)
          ! selection on l2: l1=3, m1=0
          select case (l2)
            case (0)
              ! selection on m2: l1=3, m1=0, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =1.17236802495868785d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*(5.0d0* &
(dy**2+dx**2)-2.0d0*d2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=3, m1=0, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =2.03060098439762498d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.6d1*a1*a2**2*d2+1.6d1*a1**2*a2*d2+2.7d1*a2**2+5.4d1*a1*a2 &
+2.7d1*a1**2)*dz
                case (0)
                  rlYlm_laplacian = &
-1.01530049219881249d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-8)* &
(2.0d1*a1**2*a2**2*d2*dy**4-1.1d2*a1*a2**2*dy**4-1.1d2*a1**2*a2*dy**4 &
+4.0d1*a1**2*a2**2*d2*dx**2*dy**2-2.2d2*a1*a2**2*dx**2*dy**2 &
-2.2d2*a1**2*a2*dx**2*dy**2-2.8d1*a1**2*a2**2*d2**2*dy**2 &
+1.72d2*a1*a2**2*d2*dy**2+1.72d2*a1**2*a2*d2*dy**2-8.1d1*a2**2*dy**2 &
-1.62d2*a1*a2*dy**2-8.1d1*a1**2*dy**2+2.0d1*a1**2*a2**2*d2*dx**4 &
-1.1d2*a1*a2**2*dx**4-1.1d2*a1**2*a2*dx**4 &
-2.8d1*a1**2*a2**2*d2**2*dx**2+1.72d2*a1*a2**2*d2*dx**2 &
+1.72d2*a1**2*a2*d2*dx**2-8.1d1*a2**2*dx**2-1.62d2*a1*a2*dx**2 &
-8.1d1*a1**2*dx**2+8.0d0*a1**2*a2**2*d2**3-5.6d1*a1*a2**2*d2**2 &
-5.6d1*a1**2*a2*d2**2+5.4d1*a2**2*d2+1.08d2*a1*a2*d2+5.4d1*a1**2*d2)
                case (1)
                  rlYlm_laplacian =2.03060098439762498d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.6d1*a1*a2**2*d2+1.6d1*a1**2*a2*d2+2.7d1*a2**2+5.4d1*a1*a2 &
+2.7d1*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=3, m1=0, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =4.54056183629107931d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-6.5d1*a1*a2**2*dy**2-6.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-6.5d1*a1*a2**2*dx**2-6.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.4d1*a1*a2**2*d2+1.4d1*a1**2*a2*d2+6.6d1*a2**2+1.32d2*a1*a2 &
+6.6d1*a1**2)*dz
                case (-1)
                  rlYlm_laplacian = &
-2.27028091814553966d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.3d2*a1**2*a2**3*dy**4 &
-1.3d2*a1**3*a2**2*dy**4+4.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-2.6d2*a1**2*a2**3*dx**2*dy**2-2.6d2*a1**3*a2**2*dx**2*dy**2 &
-2.8d1*a1**3*a2**3*d2**2*dy**2+1.88d2*a1**2*a2**3*d2*dy**2 &
+1.88d2*a1**3*a2**2*d2*dy**2-3.3d1*a1*a2**3*dy**2 &
-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.3d2*a1**2*a2**3*dx**4 &
-1.3d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+1.88d2*a1**2*a2**3*d2*dx**2+1.88d2*a1**3*a2**2*d2*dx**2 &
-3.3d1*a1*a2**3*dx**2-6.6d1*a1**2*a2**2*dx**2-3.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-5.2d1*a1**2*a2**3*d2**2 &
-5.2d1*a1**3*a2**2*d2**2-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2 &
-6.0d0*a1**3*a2*d2+2.7d1*a2**3+8.1d1*a1*a2**2+8.1d1*a1**2*a2 &
+2.7d1*a1**3)
                case (0)
                  rlYlm_laplacian = &
-1.31074729922739806d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)* &
(3.0d1*a1**3*a2**3*d2*dy**4-1.95d2*a1**2*a2**3*dy**4 &
-1.95d2*a1**3*a2**2*dy**4+6.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.9d2*a1**2*a2**3*dx**2*dy**2-3.9d2*a1**3*a2**2*dx**2*dy**2 &
-3.2d1*a1**3*a2**3*d2**2*dy**2+2.32d2*a1**2*a2**3*d2*dy**2 &
+2.32d2*a1**3*a2**2*d2*dy**2-1.32d2*a1*a2**3*dy**2 &
-2.64d2*a1**2*a2**2*dy**2-1.32d2*a1**3*a2*dy**2 &
+3.0d1*a1**3*a2**3*d2*dx**4-1.95d2*a1**2*a2**3*dx**4 &
-1.95d2*a1**3*a2**2*dx**4-3.2d1*a1**3*a2**3*d2**2*dx**2 &
+2.32d2*a1**2*a2**3*d2*dx**2+2.32d2*a1**3*a2**2*d2*dx**2 &
-1.32d2*a1*a2**3*dx**2-2.64d2*a1**2*a2**2*dx**2-1.32d2*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-7.6d1*a1**2*a2**3*d2**2 &
-7.6d1*a1**3*a2**2*d2**2+1.5d2*a1*a2**3*d2+3.0d2*a1**2*a2**2*d2 &
+1.5d2*a1**3*a2*d2-8.1d1*a2**3-2.43d2*a1*a2**2-2.43d2*a1**2*a2 &
-8.1d1*a1**3)*dz
                case (1)
                  rlYlm_laplacian = &
-2.27028091814553966d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.3d2*a1**2*a2**3*dy**4 &
-1.3d2*a1**3*a2**2*dy**4+4.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-2.6d2*a1**2*a2**3*dx**2*dy**2-2.6d2*a1**3*a2**2*dx**2*dy**2 &
-2.8d1*a1**3*a2**3*d2**2*dy**2+1.88d2*a1**2*a2**3*d2*dy**2 &
+1.88d2*a1**3*a2**2*d2*dy**2-3.3d1*a1*a2**3*dy**2 &
-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.3d2*a1**2*a2**3*dx**4 &
-1.3d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+1.88d2*a1**2*a2**3*d2*dx**2+1.88d2*a1**3*a2**2*d2*dx**2 &
-3.3d1*a1*a2**3*dx**2-6.6d1*a1**2*a2**2*dx**2-3.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-5.2d1*a1**2*a2**3*d2**2 &
-5.2d1*a1**3*a2**2*d2**2-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2 &
-6.0d0*a1**3*a2*d2+2.7d1*a2**3+8.1d1*a1*a2**2+8.1d1*a1**2*a2 &
+2.7d1*a1**3)
                case (2)
                  rlYlm_laplacian = &
-2.27028091814553966d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-9)*(dy+dx &
)*(dy-1.0d0*dx)*(1.0d1*a1**2*a2**2*d2*dy**2-6.5d1*a1*a2**2*dy**2 &
-6.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2-6.5d1*a1*a2**2*dx**2 &
-6.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+1.4d1*a1*a2**2*d2 &
+1.4d1*a1**2*a2*d2+6.6d1*a2**2+1.32d2*a1*a2+6.6d1*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=3, m1=0, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-2.45218365717409381d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(dy**2-3.0d0*dx**2)*(1.0d1*a1**2*a2**2*d2*dy**2-7.5d1*a1*a2**2*dy**2 &
-7.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2-7.5d1*a1*a2**2*dx**2 &
-7.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+1.2d1*a1*a2**2*d2 &
+1.2d1*a1**2*a2*d2+1.17d2*a2**2+2.34d2*a1*a2+1.17d2*a1**2)*dz
                case (-2)
                  rlYlm_laplacian = &
-6.00659871566848402d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.5d2*a1**2*a2**3*dy**4 &
-1.5d2*a1**3*a2**2*dy**4+4.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.0d2*a1**2*a2**3*dx**2*dy**2-3.0d2*a1**3*a2**2*dx**2*dy**2 &
-2.8d1*a1**3*a2**3*d2**2*dy**2+2.04d2*a1**2*a2**3*d2*dy**2 &
+2.04d2*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.5d2*a1**2*a2**3*dx**4 &
-1.5d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+2.04d2*a1**2*a2**3*d2*dx**2+2.04d2*a1**3*a2**2*d2*dx**2 &
+3.9d1*a1*a2**3*dx**2+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-4.8d1*a1**2*a2**3*d2**2 &
-4.8d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2-1.8d2*a1**2*a2**2*d2 &
-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2+1.98d2*a1**2*a2 &
+6.6d1*a1**3)
                case (-1)
                  rlYlm_laplacian = &
-1.89945329321545261d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(5.0d1*a1**3*a2**3*d2*dy**4-3.75d2*a1**2*a2**3*dy**4 &
-3.75d2*a1**3*a2**2*dy**4+1.0d2*a1**3*a2**3*d2*dx**2*dy**2 &
-7.5d2*a1**2*a2**3*dx**2*dy**2-7.5d2*a1**3*a2**2*dx**2*dy**2 &
-6.0d1*a1**3*a2**3*d2**2*dy**2+4.8d2*a1**2*a2**3*d2*dy**2 &
+4.8d2*a1**3*a2**2*d2*dy**2-1.95d2*a1*a2**3*dy**2 &
-3.9d2*a1**2*a2**2*dy**2-1.95d2*a1**3*a2*dy**2 &
+5.0d1*a1**3*a2**3*d2*dx**4-3.75d2*a1**2*a2**3*dx**4 &
-3.75d2*a1**3*a2**2*dx**4-6.0d1*a1**3*a2**3*d2**2*dx**2 &
+4.8d2*a1**2*a2**3*d2*dx**2+4.8d2*a1**3*a2**2*d2*dx**2 &
-1.95d2*a1*a2**3*dx**2-3.9d2*a1**2*a2**2*dx**2-1.95d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.44d2*a1**2*a2**3*d2**2 &
-1.44d2*a1**3*a2**2*d2**2+1.68d2*a1*a2**3*d2+3.36d2*a1**2*a2**2*d2 &
+1.68d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2-1.98d2*a1**2*a2 &
-6.6d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian =7.75448559771163262d-1*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(1.0d2*a1**4*a2**4*d2*dy**6 &
-7.5d2*a1**3*a2**4*dy**6-7.5d2*a1**4*a2**3*dy**6 &
+3.0d2*a1**4*a2**4*d2*dx**2*dy**4-2.25d3*a1**3*a2**4*dx**2*dy**4 &
-2.25d3*a1**4*a2**3*dx**2*dy**4-1.8d2*a1**4*a2**4*d2**2*dy**4 &
+1.44d3*a1**3*a2**4*d2*dy**4+1.44d3*a1**4*a2**3*d2*dy**4 &
-5.85d2*a1**2*a2**4*dy**4-1.17d3*a1**3*a2**3*dy**4 &
-5.85d2*a1**4*a2**2*dy**4+3.0d2*a1**4*a2**4*d2*dx**4*dy**2 &
-2.25d3*a1**3*a2**4*dx**4*dy**2-2.25d3*a1**4*a2**3*dx**4*dy**2 &
-3.6d2*a1**4*a2**4*d2**2*dx**2*dy**2+2.88d3*a1**3*a2**4*d2*dx**2*dy**2 &
+2.88d3*a1**4*a2**3*d2*dx**2*dy**2-1.17d3*a1**2*a2**4*dx**2*dy**2 &
-2.34d3*a1**3*a2**3*dx**2*dy**2-1.17d3*a1**4*a2**2*dx**2*dy**2 &
+9.6d1*a1**4*a2**4*d2**3*dy**2-8.64d2*a1**3*a2**4*d2**2*dy**2 &
-8.64d2*a1**4*a2**3*d2**2*dy**2+1.008d3*a1**2*a2**4*d2*dy**2 &
+2.016d3*a1**3*a2**3*d2*dy**2+1.008d3*a1**4*a2**2*d2*dy**2 &
-3.96d2*a1*a2**4*dy**2-1.188d3*a1**2*a2**3*dy**2 &
-1.188d3*a1**3*a2**2*dy**2-3.96d2*a1**4*a2*dy**2 &
+1.0d2*a1**4*a2**4*d2*dx**6-7.5d2*a1**3*a2**4*dx**6 &
-7.5d2*a1**4*a2**3*dx**6-1.8d2*a1**4*a2**4*d2**2*dx**4 &
+1.44d3*a1**3*a2**4*d2*dx**4+1.44d3*a1**4*a2**3*d2*dx**4 &
-5.85d2*a1**2*a2**4*dx**4-1.17d3*a1**3*a2**3*dx**4 &
-5.85d2*a1**4*a2**2*dx**4+9.6d1*a1**4*a2**4*d2**3*dx**2 &
-8.64d2*a1**3*a2**4*d2**2*dx**2-8.64d2*a1**4*a2**3*d2**2*dx**2 &
+1.008d3*a1**2*a2**4*d2*dx**2+2.016d3*a1**3*a2**3*d2*dx**2 &
+1.008d3*a1**4*a2**2*d2*dx**2-3.96d2*a1*a2**4*dx**2 &
-1.188d3*a1**2*a2**3*dx**2-1.188d3*a1**3*a2**2*dx**2 &
-3.96d2*a1**4*a2*dx**2-1.6d1*a1**4*a2**4*d2**4 &
+1.92d2*a1**3*a2**4*d2**3+1.92d2*a1**4*a2**3*d2**3 &
-5.76d2*a1**2*a2**4*d2**2-1.152d3*a1**3*a2**3*d2**2 &
-5.76d2*a1**4*a2**2*d2**2+6.24d2*a1*a2**4*d2+1.872d3*a1**2*a2**3*d2 &
+1.872d3*a1**3*a2**2*d2+6.24d2*a1**4*a2*d2-1.35d2*a2**4-5.4d2*a1*a2**3 &
-8.1d2*a1**2*a2**2-5.4d2*a1**3*a2-1.35d2*a1**4)
                case (1)
                  rlYlm_laplacian = &
-1.89945329321545261d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(5.0d1*a1**3*a2**3*d2*dy**4-3.75d2*a1**2*a2**3*dy**4 &
-3.75d2*a1**3*a2**2*dy**4+1.0d2*a1**3*a2**3*d2*dx**2*dy**2 &
-7.5d2*a1**2*a2**3*dx**2*dy**2-7.5d2*a1**3*a2**2*dx**2*dy**2 &
-6.0d1*a1**3*a2**3*d2**2*dy**2+4.8d2*a1**2*a2**3*d2*dy**2 &
+4.8d2*a1**3*a2**2*d2*dy**2-1.95d2*a1*a2**3*dy**2 &
-3.9d2*a1**2*a2**2*dy**2-1.95d2*a1**3*a2*dy**2 &
+5.0d1*a1**3*a2**3*d2*dx**4-3.75d2*a1**2*a2**3*dx**4 &
-3.75d2*a1**3*a2**2*dx**4-6.0d1*a1**3*a2**3*d2**2*dx**2 &
+4.8d2*a1**2*a2**3*d2*dx**2+4.8d2*a1**3*a2**2*d2*dx**2 &
-1.95d2*a1*a2**3*dx**2-3.9d2*a1**2*a2**2*dx**2-1.95d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.44d2*a1**2*a2**3*d2**2 &
-1.44d2*a1**3*a2**2*d2**2+1.68d2*a1*a2**3*d2+3.36d2*a1**2*a2**2*d2 &
+1.68d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2-1.98d2*a1**2*a2 &
-6.6d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian =3.00329935783424201d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(dy+dx)*(dy-1.0d0*dx)* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.5d2*a1**2*a2**3*dy**4 &
-1.5d2*a1**3*a2**2*dy**4+4.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.0d2*a1**2*a2**3*dx**2*dy**2-3.0d2*a1**3*a2**2*dx**2*dy**2 &
-2.8d1*a1**3*a2**3*d2**2*dy**2+2.04d2*a1**2*a2**3*d2*dy**2 &
+2.04d2*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.5d2*a1**2*a2**3*dx**4 &
-1.5d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+2.04d2*a1**2*a2**3*d2*dx**2+2.04d2*a1**3*a2**2*d2*dx**2 &
+3.9d1*a1*a2**3*dx**2+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-4.8d1*a1**2*a2**3*d2**2 &
-4.8d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2-1.8d2*a1**2*a2**2*d2 &
-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2+1.98d2*a1**2*a2 &
+6.6d1*a1**3)
                case (3)
                  rlYlm_laplacian = &
-2.45218365717409381d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*(1.0d1*a1**2*a2**2*d2*dy**2 &
-7.5d1*a1*a2**2*dy**2-7.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-7.5d1*a1*a2**2*dx**2-7.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.2d1*a1*a2**2*d2+1.2d1*a1**2*a2*d2+1.17d2*a2**2+2.34d2*a1*a2 &
+1.17d2*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=3, m1=0, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-1.04037341562157789d1*E*a1**4*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(1.0d1*a1**2*a2**2*d2*dy**2-8.5d1*a1*a2**2*dy**2 &
-8.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2-8.5d1*a1*a2**2*dx**2 &
-8.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+1.0d1*a1*a2**2*d2 &
+1.0d1*a1**2*a2*d2+1.8d2*a2**2+3.6d2*a1*a2+1.8d2*a1**2)*dz
                case (-3)
                  rlYlm_laplacian =3.67827548576114071d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(dy**2-3.0d0*dx**2)* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.7d2*a1**2*a2**3*dy**4 &
-1.7d2*a1**3*a2**2*dy**4+4.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.4d2*a1**2*a2**3*dx**2*dy**2-3.4d2*a1**3*a2**2*dx**2*dy**2 &
-2.8d1*a1**3*a2**3*d2**2*dy**2+2.2d2*a1**2*a2**3*d2*dy**2 &
+2.2d2*a1**3*a2**2*d2*dy**2+1.35d2*a1*a2**3*dy**2 &
+2.7d2*a1**2*a2**2*dy**2+1.35d2*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.7d2*a1**2*a2**3*dx**4 &
-1.7d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+2.2d2*a1**2*a2**3*d2*dx**2+2.2d2*a1**3*a2**2*d2*dx**2 &
+1.35d2*a1*a2**3*dx**2+2.7d2*a1**2*a2**2*dx**2+1.35d2*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-4.4d1*a1**2*a2**3*d2**2 &
-4.4d1*a1**3*a2**2*d2**2-1.98d2*a1*a2**3*d2-3.96d2*a1**2*a2**2*d2 &
-1.98d2*a1**3*a2*d2+1.17d2*a2**3+3.51d2*a1*a2**2+3.51d2*a1**2*a2 &
+1.17d2*a1**3)
                case (-2)
                  rlYlm_laplacian = &
-3.93224189768219417d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.8d1*a1**3*a2**3*d2**2*dy**2+7.6d2*a1**2*a2**3*d2*dy**2 &
+7.6d2*a1**3*a2**2*d2*dy**2-9.0d1*a1*a2**3*dy**2 &
-1.8d2*a1**2*a2**2*dy**2-9.0d1*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.8d1*a1**3*a2**3*d2**2*dx**2 &
+7.6d2*a1**2*a2**3*d2*dx**2+7.6d2*a1**3*a2**2*d2*dx**2 &
-9.0d1*a1*a2**3*dx**2-1.8d2*a1**2*a2**2*dx**2-9.0d1*a1**3*a2*dx**2 &
+2.4d1*a1**3*a2**3*d2**3-2.04d2*a1**2*a2**3*d2**2 &
-2.04d2*a1**3*a2**2*d2**2-1.8d1*a1*a2**3*d2-3.6d1*a1**2*a2**2*d2 &
-1.8d1*a1**3*a2*d2+1.17d2*a2**3+3.51d2*a1*a2**2+3.51d2*a1**2*a2 &
+1.17d2*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =1.39025745555846884d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(1.4d2*a1**4*a2**4*d2*dy**6 &
-1.19d3*a1**3*a2**4*dy**6-1.19d3*a1**4*a2**3*dy**6 &
+4.2d2*a1**4*a2**4*d2*dx**2*dy**4-3.57d3*a1**3*a2**4*dx**2*dy**4 &
-3.57d3*a1**4*a2**3*dx**2*dy**4-2.76d2*a1**4*a2**4*d2**2*dy**4 &
+2.46d3*a1**3*a2**4*d2*dy**4+2.46d3*a1**4*a2**3*d2*dy**4 &
-8.55d2*a1**2*a2**4*dy**4-1.71d3*a1**3*a2**3*dy**4 &
-8.55d2*a1**4*a2**2*dy**4+4.2d2*a1**4*a2**4*d2*dx**4*dy**2 &
-3.57d3*a1**3*a2**4*dx**4*dy**2-3.57d3*a1**4*a2**3*dx**4*dy**2 &
-5.52d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+4.92d3*a1**3*a2**4*d2*dx**2*dy**2+4.92d3*a1**4*a2**3*d2*dx**2*dy**2 &
-1.71d3*a1**2*a2**4*dx**2*dy**2-3.42d3*a1**3*a2**3*dx**2*dy**2 &
-1.71d3*a1**4*a2**2*dx**2*dy**2+1.68d2*a1**4*a2**4*d2**3*dy**2 &
-1.62d3*a1**3*a2**4*d2**2*dy**2-1.62d3*a1**4*a2**3*d2**2*dy**2 &
+1.53d3*a1**2*a2**4*d2*dy**2+3.06d3*a1**3*a2**3*d2*dy**2 &
+1.53d3*a1**4*a2**2*d2*dy**2-5.85d2*a1*a2**4*dy**2 &
-1.755d3*a1**2*a2**3*dy**2-1.755d3*a1**3*a2**2*dy**2 &
-5.85d2*a1**4*a2*dy**2+1.4d2*a1**4*a2**4*d2*dx**6 &
-1.19d3*a1**3*a2**4*dx**6-1.19d3*a1**4*a2**3*dx**6 &
-2.76d2*a1**4*a2**4*d2**2*dx**4+2.46d3*a1**3*a2**4*d2*dx**4 &
+2.46d3*a1**4*a2**3*d2*dx**4-8.55d2*a1**2*a2**4*dx**4 &
-1.71d3*a1**3*a2**3*dx**4-8.55d2*a1**4*a2**2*dx**4 &
+1.68d2*a1**4*a2**4*d2**3*dx**2-1.62d3*a1**3*a2**4*d2**2*dx**2 &
-1.62d3*a1**4*a2**3*d2**2*dx**2+1.53d3*a1**2*a2**4*d2*dx**2 &
+3.06d3*a1**3*a2**3*d2*dx**2+1.53d3*a1**4*a2**2*d2*dx**2 &
-5.85d2*a1*a2**4*dx**2-1.755d3*a1**2*a2**3*dx**2 &
-1.755d3*a1**3*a2**2*dx**2-5.85d2*a1**4*a2*dx**2 &
-3.2d1*a1**4*a2**4*d2**4+3.68d2*a1**3*a2**4*d2**3 &
+3.68d2*a1**4*a2**3*d2**3-8.64d2*a1**2*a2**4*d2**2 &
-1.728d3*a1**3*a2**3*d2**2-8.64d2*a1**4*a2**2*d2**2+9.96d2*a1*a2**4*d2 &
+2.988d3*a1**2*a2**3*d2+2.988d3*a1**3*a2**2*d2+9.96d2*a1**4*a2*d2 &
-3.3d2*a2**4-1.32d3*a1*a2**3-1.98d3*a1**2*a2**2-1.32d3*a1**3*a2 &
-3.3d2*a1**4)
                case (0)
                  rlYlm_laplacian =4.39638009359507944d-1*E*a1**2*a2*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*(3.5d2*a1**4*a2**4*d2*dy**6 &
-2.975d3*a1**3*a2**4*dy**6-2.975d3*a1**4*a2**3*dy**6 &
+1.05d3*a1**4*a2**4*d2*dx**2*dy**4-8.925d3*a1**3*a2**4*dx**2*dy**4 &
-8.925d3*a1**4*a2**3*dx**2*dy**4-5.4d2*a1**4*a2**4*d2**2*dy**4 &
+4.95d3*a1**3*a2**4*d2*dy**4+4.95d3*a1**4*a2**3*d2*dy**4 &
-2.7d3*a1**2*a2**4*dy**4-5.4d3*a1**3*a2**3*dy**4 &
-2.7d3*a1**4*a2**2*dy**4+1.05d3*a1**4*a2**4*d2*dx**4*dy**2 &
-8.925d3*a1**3*a2**4*dx**4*dy**2-8.925d3*a1**4*a2**3*dx**4*dy**2 &
-1.08d3*a1**4*a2**4*d2**2*dx**2*dy**2+9.9d3*a1**3*a2**4*d2*dx**2*dy**2 &
+9.9d3*a1**4*a2**3*d2*dx**2*dy**2-5.4d3*a1**2*a2**4*dx**2*dy**2 &
-1.08d4*a1**3*a2**3*dx**2*dy**2-5.4d3*a1**4*a2**2*dx**2*dy**2 &
+2.4d2*a1**4*a2**4*d2**3*dy**2-2.52d3*a1**3*a2**4*d2**2*dy**2 &
-2.52d3*a1**4*a2**3*d2**2*dy**2+3.96d3*a1**2*a2**4*d2*dy**2 &
+7.92d3*a1**3*a2**3*d2*dy**2+3.96d3*a1**4*a2**2*d2*dy**2 &
-2.34d3*a1*a2**4*dy**2-7.02d3*a1**2*a2**3*dy**2 &
-7.02d3*a1**3*a2**2*dy**2-2.34d3*a1**4*a2*dy**2 &
+3.5d2*a1**4*a2**4*d2*dx**6-2.975d3*a1**3*a2**4*dx**6 &
-2.975d3*a1**4*a2**3*dx**6-5.4d2*a1**4*a2**4*d2**2*dx**4 &
+4.95d3*a1**3*a2**4*d2*dx**4+4.95d3*a1**4*a2**3*d2*dx**4 &
-2.7d3*a1**2*a2**4*dx**4-5.4d3*a1**3*a2**3*dx**4 &
-2.7d3*a1**4*a2**2*dx**4+2.4d2*a1**4*a2**4*d2**3*dx**2 &
-2.52d3*a1**3*a2**4*d2**2*dx**2-2.52d3*a1**4*a2**3*d2**2*dx**2 &
+3.96d3*a1**2*a2**4*d2*dx**2+7.92d3*a1**3*a2**3*d2*dx**2 &
+3.96d3*a1**4*a2**2*d2*dx**2-2.34d3*a1*a2**4*dx**2 &
-7.02d3*a1**2*a2**3*dx**2-7.02d3*a1**3*a2**2*dx**2 &
-2.34d3*a1**4*a2*dx**2-3.2d1*a1**4*a2**4*d2**4 &
+4.64d2*a1**3*a2**4*d2**3+4.64d2*a1**4*a2**3*d2**3 &
-1.872d3*a1**2*a2**4*d2**2-3.744d3*a1**3*a2**3*d2**2 &
-1.872d3*a1**4*a2**2*d2**2+3.048d3*a1*a2**4*d2+9.144d3*a1**2*a2**3*d2 &
+9.144d3*a1**3*a2**2*d2+3.048d3*a1**4*a2*d2-1.32d3*a2**4 &
-5.28d3*a1*a2**3-7.92d3*a1**2*a2**2-5.28d3*a1**3*a2-1.32d3*a1**4)*dz
                case (1)
                  rlYlm_laplacian =1.39025745555846884d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(1.4d2*a1**4*a2**4*d2*dy**6 &
-1.19d3*a1**3*a2**4*dy**6-1.19d3*a1**4*a2**3*dy**6 &
+4.2d2*a1**4*a2**4*d2*dx**2*dy**4-3.57d3*a1**3*a2**4*dx**2*dy**4 &
-3.57d3*a1**4*a2**3*dx**2*dy**4-2.76d2*a1**4*a2**4*d2**2*dy**4 &
+2.46d3*a1**3*a2**4*d2*dy**4+2.46d3*a1**4*a2**3*d2*dy**4 &
-8.55d2*a1**2*a2**4*dy**4-1.71d3*a1**3*a2**3*dy**4 &
-8.55d2*a1**4*a2**2*dy**4+4.2d2*a1**4*a2**4*d2*dx**4*dy**2 &
-3.57d3*a1**3*a2**4*dx**4*dy**2-3.57d3*a1**4*a2**3*dx**4*dy**2 &
-5.52d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+4.92d3*a1**3*a2**4*d2*dx**2*dy**2+4.92d3*a1**4*a2**3*d2*dx**2*dy**2 &
-1.71d3*a1**2*a2**4*dx**2*dy**2-3.42d3*a1**3*a2**3*dx**2*dy**2 &
-1.71d3*a1**4*a2**2*dx**2*dy**2+1.68d2*a1**4*a2**4*d2**3*dy**2 &
-1.62d3*a1**3*a2**4*d2**2*dy**2-1.62d3*a1**4*a2**3*d2**2*dy**2 &
+1.53d3*a1**2*a2**4*d2*dy**2+3.06d3*a1**3*a2**3*d2*dy**2 &
+1.53d3*a1**4*a2**2*d2*dy**2-5.85d2*a1*a2**4*dy**2 &
-1.755d3*a1**2*a2**3*dy**2-1.755d3*a1**3*a2**2*dy**2 &
-5.85d2*a1**4*a2*dy**2+1.4d2*a1**4*a2**4*d2*dx**6 &
-1.19d3*a1**3*a2**4*dx**6-1.19d3*a1**4*a2**3*dx**6 &
-2.76d2*a1**4*a2**4*d2**2*dx**4+2.46d3*a1**3*a2**4*d2*dx**4 &
+2.46d3*a1**4*a2**3*d2*dx**4-8.55d2*a1**2*a2**4*dx**4 &
-1.71d3*a1**3*a2**3*dx**4-8.55d2*a1**4*a2**2*dx**4 &
+1.68d2*a1**4*a2**4*d2**3*dx**2-1.62d3*a1**3*a2**4*d2**2*dx**2 &
-1.62d3*a1**4*a2**3*d2**2*dx**2+1.53d3*a1**2*a2**4*d2*dx**2 &
+3.06d3*a1**3*a2**3*d2*dx**2+1.53d3*a1**4*a2**2*d2*dx**2 &
-5.85d2*a1*a2**4*dx**2-1.755d3*a1**2*a2**3*dx**2 &
-1.755d3*a1**3*a2**2*dx**2-5.85d2*a1**4*a2*dx**2 &
-3.2d1*a1**4*a2**4*d2**4+3.68d2*a1**3*a2**4*d2**3 &
+3.68d2*a1**4*a2**3*d2**3-8.64d2*a1**2*a2**4*d2**2 &
-1.728d3*a1**3*a2**3*d2**2-8.64d2*a1**4*a2**2*d2**2+9.96d2*a1*a2**4*d2 &
+2.988d3*a1**2*a2**3*d2+2.988d3*a1**3*a2**2*d2+9.96d2*a1**4*a2*d2 &
-3.3d2*a2**4-1.32d3*a1*a2**3-1.98d3*a1**2*a2**2-1.32d3*a1**3*a2 &
-3.3d2*a1**4)
                case (2)
                  rlYlm_laplacian =1.96612094884109708d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(dy+dx)*(dy-1.0d0*dx)* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.8d1*a1**3*a2**3*d2**2*dy**2+7.6d2*a1**2*a2**3*d2*dy**2 &
+7.6d2*a1**3*a2**2*d2*dy**2-9.0d1*a1*a2**3*dy**2 &
-1.8d2*a1**2*a2**2*dy**2-9.0d1*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.8d1*a1**3*a2**3*d2**2*dx**2 &
+7.6d2*a1**2*a2**3*d2*dx**2+7.6d2*a1**3*a2**2*d2*dx**2 &
-9.0d1*a1*a2**3*dx**2-1.8d2*a1**2*a2**2*dx**2-9.0d1*a1**3*a2*dx**2 &
+2.4d1*a1**3*a2**3*d2**3-2.04d2*a1**2*a2**3*d2**2 &
-2.04d2*a1**3*a2**2*d2**2-1.8d1*a1*a2**3*d2-3.6d1*a1**2*a2**2*d2 &
-1.8d1*a1**3*a2*d2+1.17d2*a2**3+3.51d2*a1*a2**2+3.51d2*a1**2*a2 &
+1.17d2*a1**3)*dz
                case (3)
                  rlYlm_laplacian =3.67827548576114071d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(3.0d0*dy**2-1.0d0*dx**2)* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.7d2*a1**2*a2**3*dy**4 &
-1.7d2*a1**3*a2**2*dy**4+4.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.4d2*a1**2*a2**3*dx**2*dy**2-3.4d2*a1**3*a2**2*dx**2*dy**2 &
-2.8d1*a1**3*a2**3*d2**2*dy**2+2.2d2*a1**2*a2**3*d2*dy**2 &
+2.2d2*a1**3*a2**2*d2*dy**2+1.35d2*a1*a2**3*dy**2 &
+2.7d2*a1**2*a2**2*dy**2+1.35d2*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.7d2*a1**2*a2**3*dx**4 &
-1.7d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+2.2d2*a1**2*a2**3*d2*dx**2+2.2d2*a1**3*a2**2*d2*dx**2 &
+1.35d2*a1*a2**3*dx**2+2.7d2*a1**2*a2**2*dx**2+1.35d2*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-4.4d1*a1**2*a2**3*d2**2 &
-4.4d1*a1**3*a2**2*d2**2-1.98d2*a1*a2**3*d2-3.96d2*a1**2*a2**2*d2 &
-1.98d2*a1**3*a2*d2+1.17d2*a2**3+3.51d2*a1*a2**2+3.51d2*a1**2*a2 &
+1.17d2*a1**3)
                case (4)
                  rlYlm_laplacian =2.60093353905394473d0*E*a1**4*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*(dy**2-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2 &
+2.0d0*dx*dy-1.0d0*dx**2)*(1.0d1*a1**2*a2**2*d2*dy**2 &
-8.5d1*a1*a2**2*dy**2-8.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-8.5d1*a1*a2**2*dx**2-8.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.0d1*a1*a2**2*d2+1.0d1*a1**2*a2*d2+1.8d2*a2**2+3.6d2*a1*a2 &
+1.8d2*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (1)
          ! selection on l2: l1=3, m1=1
          select case (l2)
            case (0)
              ! selection on m2: l1=3, m1=1, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =1.43585172595163941d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx* &
(5.0d0*(dy**2+dx**2)-4.0d0*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=3, m1=1, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =2.48696814148370333d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+4.2d1*a1*a2**2*d2+4.2d1*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2 &
+9.0d0*a1**2)
                case (0)
                  rlYlm_laplacian =2.48696814148370333d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*(1.0d1*a1**2*a2**2*d2*dy**2 &
-5.5d1*a1*a2**2*dy**2-5.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-5.5d1*a1*a2**2*dx**2-5.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+5.2d1*a1*a2**2*d2+5.2d1*a1**2*a2*d2-3.6d1*a2**2-7.2d1*a1*a2 &
-3.6d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =1.24348407074185167d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*(2.0d1*a1**2*a2**2*d2*dx**2*dy**2 &
-1.1d2*a1*a2**2*dx**2*dy**2-1.1d2*a1**2*a2*dx**2*dy**2 &
-1.0d1*a1*a2**2*d2*dy**2-1.0d1*a1**2*a2*d2*dy**2+4.5d1*a2**2*dy**2 &
+9.0d1*a1*a2*dy**2+4.5d1*a1**2*dy**2+2.0d1*a1**2*a2**2*d2*dx**4 &
-1.1d2*a1*a2**2*dx**4-1.1d2*a1**2*a2*dx**4 &
-1.6d1*a1**2*a2**2*d2**2*dx**2+7.4d1*a1*a2**2*d2*dx**2 &
+7.4d1*a1**2*a2*d2*dx**2+6.3d1*a2**2*dx**2+1.26d2*a1*a2*dx**2 &
+6.3d1*a1**2*dx**2+8.0d0*a1*a2**2*d2**2+8.0d0*a1**2*a2*d2**2 &
-3.6d1*a2**2*d2-7.2d1*a1*a2*d2-3.6d1*a1**2*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=3, m1=1, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =2.78051491111693767d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.3d2*a1**2*a2**3*dx**2*dy**2-1.3d2*a1**3*a2**2*dx**2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dy**2-1.0d1*a1**3*a2**2*d2*dy**2 &
+5.5d1*a1*a2**3*dy**2+1.1d2*a1**2*a2**2*dy**2+5.5d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.3d2*a1**2*a2**3*dx**4 &
-1.3d2*a1**3*a2**2*dx**4-1.6d1*a1**3*a2**3*d2**2*dx**2 &
+8.6d1*a1**2*a2**3*d2*dx**2+8.6d1*a1**3*a2**2*d2*dx**2 &
+9.9d1*a1*a2**3*dx**2+1.98d2*a1**2*a2**2*dx**2+9.9d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.2d1*a1*a2**3*d2 &
-8.4d1*a1**2*a2**2*d2-4.2d1*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)
                case (-1)
                  rlYlm_laplacian =5.56102982223387534d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(1.0d1*a1**2*a2**2*d2*dy**2 &
-6.5d1*a1*a2**2*dy**2-6.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-6.5d1*a1*a2**2*dx**2-6.5d1*a1**2*a2*dx**2-8.0d0*a1**2*a2**2*d2**2 &
+5.8d1*a1*a2**2*d2+5.8d1*a1**2*a2*d2-3.3d1*a2**2-6.6d1*a1*a2 &
-3.3d1*a1**2)*dz
                case (0)
                  rlYlm_laplacian = &
-1.60533103241913232d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(3.0d1*a1**3*a2**3*d2*dy**4-1.95d2*a1**2*a2**3*dy**4 &
-1.95d2*a1**3*a2**2*dy**4+6.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.9d2*a1**2*a2**3*dx**2*dy**2-3.9d2*a1**3*a2**2*dx**2*dy**2 &
-4.4d1*a1**3*a2**3*d2**2*dy**2+3.04d2*a1**2*a2**3*d2*dy**2 &
+3.04d2*a1**3*a2**2*d2*dy**2-9.9d1*a1*a2**3*dy**2 &
-1.98d2*a1**2*a2**2*dy**2-9.9d1*a1**3*a2*dy**2 &
+3.0d1*a1**3*a2**3*d2*dx**4-1.95d2*a1**2*a2**3*dx**4 &
-1.95d2*a1**3*a2**2*dx**4-4.4d1*a1**3*a2**3*d2**2*dx**2 &
+3.04d2*a1**2*a2**3*d2*dx**2+3.04d2*a1**3*a2**2*d2*dx**2 &
-9.9d1*a1*a2**3*dx**2-1.98d2*a1**2*a2**2*dx**2-9.9d1*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.28d2*a1**2*a2**3*d2**2 &
-1.28d2*a1**3*a2**2*d2**2+1.44d2*a1*a2**3*d2+2.88d2*a1**2*a2**2*d2 &
+1.44d2*a1**3*a2*d2-5.4d1*a2**3-1.62d2*a1*a2**2-1.62d2*a1**2*a2 &
-5.4d1*a1**3)
                case (1)
                  rlYlm_laplacian =2.78051491111693767d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*(2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.3d2*a1**2*a2**3*dx**2*dy**2-1.3d2*a1**3*a2**2*dx**2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dy**2-1.0d1*a1**3*a2**2*d2*dy**2 &
+5.5d1*a1*a2**3*dy**2+1.1d2*a1**2*a2**2*dy**2+5.5d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.3d2*a1**2*a2**3*dx**4 &
-1.3d2*a1**3*a2**2*dx**4-1.6d1*a1**3*a2**3*d2**2*dx**2 &
+1.06d2*a1**2*a2**3*d2*dx**2+1.06d2*a1**3*a2**2*d2*dx**2 &
-1.1d1*a1*a2**3*dx**2-2.2d1*a1**2*a2**2*dx**2-1.1d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.2d1*a1*a2**3*d2 &
-1.04d2*a1**2*a2**2*d2-5.2d1*a1**3*a2*d2+3.6d1*a2**3+1.08d2*a1*a2**2 &
+1.08d2*a1**2*a2+3.6d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian = &
-2.78051491111693767d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(1.0d1*a1**3*a2**3*d2*dy**4-6.5d1*a1**2*a2**3*dy**4 &
-6.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+5.8d1*a1**2*a2**3*d2*dy**2+5.8d1*a1**3*a2**2*d2*dy**2 &
-3.3d1*a1*a2**3*dy**2-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+6.5d1*a1**2*a2**3*dx**4 &
+6.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-3.8d1*a1**2*a2**3*d2*dx**2-3.8d1*a1**3*a2**2*d2*dx**2 &
-7.7d1*a1*a2**3*dx**2-1.54d2*a1**2*a2**2*dx**2-7.7d1*a1**3*a2*dx**2 &
-8.0d0*a1**2*a2**3*d2**2-8.0d0*a1**3*a2**2*d2**2+4.2d1*a1*a2**3*d2 &
+8.4d1*a1**2*a2**2*d2+4.2d1*a1**3*a2*d2+9.0d0*a2**3+2.7d1*a1*a2**2 &
+2.7d1*a1**2*a2+9.0d0*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=3, m1=1, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-3.00329935783424201d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(1.0d1*a1**3*a2**3*d2*dy**4-7.5d1*a1**2*a2**3*dy**4 &
-7.5d1*a1**3*a2**2*dy**4-2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
+1.5d2*a1**2*a2**3*dx**2*dy**2+1.5d2*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+8.4d1*a1**2*a2**3*d2*dy**2 &
+8.4d1*a1**3*a2**2*d2*dy**2-1.56d2*a1*a2**3*dy**2 &
-3.12d2*a1**2*a2**2*dy**2-1.56d2*a1**3*a2*dy**2 &
-3.0d1*a1**3*a2**3*d2*dx**4+2.25d2*a1**2*a2**3*dx**4 &
+2.25d2*a1**3*a2**2*dx**4+2.4d1*a1**3*a2**3*d2**2*dx**2 &
-1.32d2*a1**2*a2**3*d2*dx**2-1.32d2*a1**3*a2**2*d2*dx**2 &
-3.12d2*a1*a2**3*dx**2-6.24d2*a1**2*a2**2*dx**2-3.12d2*a1**3*a2*dx**2 &
-2.4d1*a1**2*a2**3*d2**2-2.4d1*a1**3*a2**2*d2**2+1.44d2*a1*a2**3*d2 &
+2.88d2*a1**2*a2**2*d2+1.44d2*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)
                case (-2)
                  rlYlm_laplacian =7.35655097152228142d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.5d2*a1**2*a2**3*dx**2*dy**2-1.5d2*a1**3*a2**2*dx**2*dy**2 &
-1.0d1*a1**2*a2**3*d2*dy**2-1.0d1*a1**3*a2**2*d2*dy**2 &
+6.5d1*a1*a2**3*dy**2+1.3d2*a1**2*a2**2*dy**2+6.5d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.5d2*a1**2*a2**3*dx**4 &
-1.5d2*a1**3*a2**2*dx**4-1.6d1*a1**3*a2**3*d2**2*dx**2 &
+1.18d2*a1**2*a2**3*d2*dx**2+1.18d2*a1**3*a2**2*d2*dx**2 &
+1.3d1*a1*a2**3*dx**2+2.6d1*a1**2*a2**2*dx**2+1.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.8d1*a1*a2**3*d2 &
-1.16d2*a1**2*a2**2*d2-5.8d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian = &
-2.32634567931348979d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(5.0d1*a1**3*a2**3*d2*dy**4-3.75d2*a1**2*a2**3*dy**4 &
-3.75d2*a1**3*a2**2*dy**4+1.0d2*a1**3*a2**3*d2*dx**2*dy**2 &
-7.5d2*a1**2*a2**3*dx**2*dy**2-7.5d2*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+6.4d2*a1**2*a2**3*d2*dy**2 &
+6.4d2*a1**3*a2**2*d2*dy**2-2.6d2*a1*a2**3*dy**2 &
-5.2d2*a1**2*a2**2*dy**2-2.6d2*a1**3*a2*dy**2 &
+5.0d1*a1**3*a2**3*d2*dx**4-3.75d2*a1**2*a2**3*dx**4 &
-3.75d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.4d2*a1**2*a2**3*d2*dx**2+6.4d2*a1**3*a2**2*d2*dx**2 &
-2.6d2*a1*a2**3*dx**2-5.2d2*a1**2*a2**2*dx**2-2.6d2*a1**3*a2*dx**2 &
+3.2d1*a1**3*a2**3*d2**3-2.88d2*a1**2*a2**3*d2**2 &
-2.88d2*a1**3*a2**2*d2**2+3.36d2*a1*a2**3*d2+6.72d2*a1**2*a2**2*d2 &
+3.36d2*a1**3*a2*d2-1.32d2*a2**3-3.96d2*a1*a2**2-3.96d2*a1**2*a2 &
-1.32d2*a1**3)
                case (0)
                  rlYlm_laplacian = &
-1.89945329321545261d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(5.0d1*a1**3*a2**3*d2*dy**4-3.75d2*a1**2*a2**3*dy**4 &
-3.75d2*a1**3*a2**2*dy**4+1.0d2*a1**3*a2**3*d2*dx**2*dy**2 &
-7.5d2*a1**2*a2**3*dx**2*dy**2-7.5d2*a1**3*a2**2*dx**2*dy**2 &
-6.0d1*a1**3*a2**3*d2**2*dy**2+4.8d2*a1**2*a2**3*d2*dy**2 &
+4.8d2*a1**3*a2**2*d2*dy**2-1.95d2*a1*a2**3*dy**2 &
-3.9d2*a1**2*a2**2*dy**2-1.95d2*a1**3*a2*dy**2 &
+5.0d1*a1**3*a2**3*d2*dx**4-3.75d2*a1**2*a2**3*dx**4 &
-3.75d2*a1**3*a2**2*dx**4-6.0d1*a1**3*a2**3*d2**2*dx**2 &
+4.8d2*a1**2*a2**3*d2*dx**2+4.8d2*a1**3*a2**2*d2*dx**2 &
-1.95d2*a1*a2**3*dx**2-3.9d2*a1**2*a2**2*dx**2-1.95d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.44d2*a1**2*a2**3*d2**2 &
-1.44d2*a1**3*a2**2*d2**2+1.68d2*a1*a2**3*d2+3.36d2*a1**2*a2**2*d2 &
+1.68d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2-1.98d2*a1**2*a2 &
-6.6d1*a1**3)*dz
                case (1)
                  rlYlm_laplacian =-1.16317283965674489d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(1.0d2*a1**4*a2**4*d2*dx**2*dy**4 &
-7.5d2*a1**3*a2**4*dx**2*dy**4-7.5d2*a1**4*a2**3*dx**2*dy**4 &
-5.0d1*a1**3*a2**4*d2*dy**4-5.0d1*a1**4*a2**3*d2*dy**4 &
+3.25d2*a1**2*a2**4*dy**4+6.5d2*a1**3*a2**3*dy**4 &
+3.25d2*a1**4*a2**2*dy**4+2.0d2*a1**4*a2**4*d2*dx**4*dy**2 &
-1.5d3*a1**3*a2**4*dx**4*dy**2-1.5d3*a1**4*a2**3*dx**4*dy**2 &
-1.6d2*a1**4*a2**4*d2**2*dx**2*dy**2+1.18d3*a1**3*a2**4*d2*dx**2*dy**2 &
+1.18d3*a1**4*a2**3*d2*dx**2*dy**2+1.3d2*a1**2*a2**4*dx**2*dy**2 &
+2.6d2*a1**3*a2**3*dx**2*dy**2+1.3d2*a1**4*a2**2*dx**2*dy**2 &
+8.0d1*a1**3*a2**4*d2**2*dy**2+8.0d1*a1**4*a2**3*d2**2*dy**2 &
-5.8d2*a1**2*a2**4*d2*dy**2-1.16d3*a1**3*a2**3*d2*dy**2 &
-5.8d2*a1**4*a2**2*d2*dy**2+3.3d2*a1*a2**4*dy**2 &
+9.9d2*a1**2*a2**3*dy**2+9.9d2*a1**3*a2**2*dy**2+3.3d2*a1**4*a2*dy**2 &
+1.0d2*a1**4*a2**4*d2*dx**6-7.5d2*a1**3*a2**4*dx**6 &
-7.5d2*a1**4*a2**3*dx**6-1.6d2*a1**4*a2**4*d2**2*dx**4 &
+1.23d3*a1**3*a2**4*d2*dx**4+1.23d3*a1**4*a2**3*d2*dx**4 &
-1.95d2*a1**2*a2**4*dx**4-3.9d2*a1**3*a2**3*dx**4 &
-1.95d2*a1**4*a2**2*dx**4+6.4d1*a1**4*a2**4*d2**3*dx**2 &
-4.96d2*a1**3*a2**4*d2**2*dx**2-4.96d2*a1**4*a2**3*d2**2*dx**2 &
+9.2d1*a1**2*a2**4*d2*dx**2+1.84d2*a1**3*a2**3*d2*dx**2 &
+9.2d1*a1**4*a2**2*d2*dx**2+6.6d1*a1*a2**4*dx**2 &
+1.98d2*a1**2*a2**3*dx**2+1.98d2*a1**3*a2**2*dx**2 &
+6.6d1*a1**4*a2*dx**2-3.2d1*a1**3*a2**4*d2**3-3.2d1*a1**4*a2**3*d2**3 &
+2.72d2*a1**2*a2**4*d2**2+5.44d2*a1**3*a2**3*d2**2 &
+2.72d2*a1**4*a2**2*d2**2-3.72d2*a1*a2**4*d2-1.116d3*a1**2*a2**3*d2 &
-1.116d3*a1**3*a2**2*d2-3.72d2*a1**4*a2*d2+9.0d1*a2**4+3.6d2*a1*a2**3 &
+5.4d2*a1**2*a2**2+3.6d2*a1**3*a2+9.0d1*a1**4)
                case (2)
                  rlYlm_laplacian = &
-7.35655097152228142d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(1.0d1*a1**3*a2**3*d2*dy**4-7.5d1*a1**2*a2**3*dy**4 &
-7.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+7.4d1*a1**2*a2**3*d2*dy**2+7.4d1*a1**3*a2**2*d2*dy**2 &
-9.1d1*a1*a2**3*dy**2-1.82d2*a1**2*a2**2*dy**2-9.1d1*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+7.5d1*a1**2*a2**3*dx**4 &
+7.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-5.4d1*a1**2*a2**3*d2*dx**2-5.4d1*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
-8.0d0*a1**2*a2**3*d2**2-8.0d0*a1**3*a2**2*d2**2+5.8d1*a1*a2**3*d2 &
+1.16d2*a1**2*a2**2*d2+5.8d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (3)
                  rlYlm_laplacian = &
-1.50164967891712101d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)* &
(6.0d1*a1**3*a2**3*d2*dx**2*dy**4-4.5d2*a1**2*a2**3*dx**2*dy**4 &
-4.5d2*a1**3*a2**2*dx**2*dy**4-3.0d1*a1**2*a2**3*d2*dy**4 &
-3.0d1*a1**3*a2**2*d2*dy**4+1.95d2*a1*a2**3*dy**4 &
+3.9d2*a1**2*a2**2*dy**4+1.95d2*a1**3*a2*dy**4 &
+4.0d1*a1**3*a2**3*d2*dx**4*dy**2-3.0d2*a1**2*a2**3*dx**4*dy**2 &
-3.0d2*a1**3*a2**2*dx**4*dy**2-4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+3.24d2*a1**2*a2**3*d2*dx**2*dy**2+3.24d2*a1**3*a2**2*d2*dx**2*dy**2 &
+2.34d2*a1*a2**3*dx**2*dy**2+4.68d2*a1**2*a2**2*dx**2*dy**2 &
+2.34d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.44d2*a1*a2**3*d2*dy**2 &
-2.88d2*a1**2*a2**2*d2*dy**2-1.44d2*a1**3*a2*d2*dy**2 &
-6.6d1*a2**3*dy**2-1.98d2*a1*a2**2*dy**2-1.98d2*a1**2*a2*dy**2 &
-6.6d1*a1**3*dy**2-2.0d1*a1**3*a2**3*d2*dx**6+1.5d2*a1**2*a2**3*dx**6 &
+1.5d2*a1**3*a2**2*dx**6+1.6d1*a1**3*a2**3*d2**2*dx**4 &
-7.8d1*a1**2*a2**3*d2*dx**4-7.8d1*a1**3*a2**2*d2*dx**4 &
-2.73d2*a1*a2**3*dx**4-5.46d2*a1**2*a2**2*dx**4-2.73d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.44d2*a1*a2**3*d2*dx**2+2.88d2*a1**2*a2**2*d2*dx**2 &
+1.44d2*a1**3*a2*d2*dx**2+6.6d1*a2**3*dx**2+1.98d2*a1*a2**2*dx**2 &
+1.98d2*a1**2*a2*dx**2+6.6d1*a1**3*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=3, m1=1, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-6.37096002557338818d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(2.0d1*a1**3*a2**3*d2*dx**2*dy**4-1.7d2*a1**2*a2**3*dx**2*dy**4 &
-1.7d2*a1**3*a2**2*dx**2*dy**4-1.0d1*a1**2*a2**3*d2*dy**4 &
-1.0d1*a1**3*a2**2*d2*dy**4+7.5d1*a1*a2**3*dy**4 &
+1.5d2*a1**2*a2**2*dy**4+7.5d1*a1**3*a2*dy**4 &
-1.6d1*a1**3*a2**3*d2**2*dx**2*dy**2+1.4d2*a1**2*a2**3*d2*dx**2*dy**2 &
+1.4d2*a1**3*a2**2*d2*dx**2*dy**2-3.0d1*a1*a2**3*dx**2*dy**2 &
-6.0d1*a1**2*a2**2*dx**2*dy**2-3.0d1*a1**3*a2*dx**2*dy**2 &
+8.0d0*a1**2*a2**3*d2**2*dy**2+8.0d0*a1**3*a2**2*d2**2*dy**2 &
-5.4d1*a1*a2**3*d2*dy**2-1.08d2*a1**2*a2**2*d2*dy**2 &
-5.4d1*a1**3*a2*d2*dy**2-3.9d1*a2**3*dy**2-1.17d2*a1*a2**2*dy**2 &
-1.17d2*a1**2*a2*dy**2-3.9d1*a1**3*dy**2-2.0d1*a1**3*a2**3*d2*dx**6 &
+1.7d2*a1**2*a2**3*dx**6+1.7d2*a1**3*a2**2*dx**6 &
+1.6d1*a1**3*a2**3*d2**2*dx**4-9.0d1*a1**2*a2**3*d2*dx**4 &
-9.0d1*a1**3*a2**2*d2*dx**4-3.45d2*a1*a2**3*dx**4 &
-6.9d2*a1**2*a2**2*dx**4-3.45d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.62d2*a1*a2**3*d2*dx**2+3.24d2*a1**2*a2**2*d2*dx**2 &
+1.62d2*a1**3*a2*d2*dx**2+1.17d2*a2**3*dx**2+3.51d2*a1*a2**2*dx**2 &
+3.51d2*a1**2*a2*dx**2+1.17d2*a1**3*dx**2)
                case (-3)
                  rlYlm_laplacian = &
-9.00989807350272604d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(1.0d1*a1**3*a2**3*d2*dy**4-8.5d1*a1**2*a2**3*dy**4 &
-8.5d1*a1**3*a2**2*dy**4-2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
+1.7d2*a1**2*a2**3*dx**2*dy**2+1.7d2*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+1.0d2*a1**2*a2**3*d2*dy**2 &
+1.0d2*a1**3*a2**2*d2*dy**2-2.4d2*a1*a2**3*dy**2 &
-4.8d2*a1**2*a2**2*dy**2-2.4d2*a1**3*a2*dy**2 &
-3.0d1*a1**3*a2**3*d2*dx**4+2.55d2*a1**2*a2**3*dx**4 &
+2.55d2*a1**3*a2**2*dx**4+2.4d1*a1**3*a2**3*d2**2*dx**2 &
-1.8d2*a1**2*a2**3*d2*dx**2-1.8d2*a1**3*a2**2*d2*dx**2 &
-1.8d2*a1*a2**3*dx**2-3.6d2*a1**2*a2**2*dx**2-1.8d2*a1**3*a2*dx**2 &
-2.4d1*a1**2*a2**3*d2**2-2.4d1*a1**3*a2**2*d2**2+1.92d2*a1*a2**3*d2 &
+3.84d2*a1**2*a2**2*d2+1.92d2*a1**3*a2*d2-7.8d1*a2**3-2.34d2*a1*a2**2 &
-2.34d2*a1**2*a2-7.8d1*a1**3)*dz
                case (-2)
                  rlYlm_laplacian = &
-2.40799654862869848d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(1.4d2*a1**4*a2**4*d2*dx**2*dy**4-1.19d3*a1**3*a2**4*dx**2*dy**4 &
-1.19d3*a1**4*a2**3*dx**2*dy**4-7.0d1*a1**3*a2**4*d2*dy**4 &
-7.0d1*a1**4*a2**3*d2*dy**4+5.25d2*a1**2*a2**4*dy**4 &
+1.05d3*a1**3*a2**3*dy**4+5.25d2*a1**4*a2**2*dy**4 &
+2.8d2*a1**4*a2**4*d2*dx**4*dy**2-2.38d3*a1**3*a2**4*dx**4*dy**2 &
-2.38d3*a1**4*a2**3*dx**4*dy**2-2.32d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.94d3*a1**3*a2**4*d2*dx**2*dy**2+1.94d3*a1**4*a2**3*d2*dx**2*dy**2 &
+2.4d2*a1**2*a2**4*dx**2*dy**2+4.8d2*a1**3*a2**3*dx**2*dy**2 &
+2.4d2*a1**4*a2**2*dx**2*dy**2+1.16d2*a1**3*a2**4*d2**2*dy**2 &
+1.16d2*a1**4*a2**3*d2**2*dy**2-9.48d2*a1**2*a2**4*d2*dy**2 &
-1.896d3*a1**3*a2**3*d2*dy**2-9.48d2*a1**4*a2**2*d2*dy**2 &
+5.07d2*a1*a2**4*dy**2+1.521d3*a1**2*a2**3*dy**2 &
+1.521d3*a1**3*a2**2*dy**2+5.07d2*a1**4*a2*dy**2 &
+1.4d2*a1**4*a2**4*d2*dx**6-1.19d3*a1**3*a2**4*dx**6 &
-1.19d3*a1**4*a2**3*dx**6-2.32d2*a1**4*a2**4*d2**2*dx**4 &
+2.01d3*a1**3*a2**4*d2*dx**4+2.01d3*a1**4*a2**3*d2*dx**4 &
-2.85d2*a1**2*a2**4*dx**4-5.7d2*a1**3*a2**3*dx**4 &
-2.85d2*a1**4*a2**2*dx**4+9.6d1*a1**4*a2**4*d2**3*dx**2 &
-8.28d2*a1**3*a2**4*d2**2*dx**2-8.28d2*a1**4*a2**3*d2**2*dx**2 &
+8.4d1*a1**2*a2**4*d2*dx**2+1.68d2*a1**3*a2**3*d2*dx**2 &
+8.4d1*a1**4*a2**2*d2*dx**2+3.9d1*a1*a2**4*dx**2 &
+1.17d2*a1**2*a2**3*dx**2+1.17d2*a1**3*a2**2*dx**2 &
+3.9d1*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+4.44d2*a1**2*a2**4*d2**2+8.88d2*a1**3*a2**3*d2**2 &
+4.44d2*a1**4*a2**2*d2**2-5.76d2*a1*a2**4*d2-1.728d3*a1**2*a2**3*d2 &
-1.728d3*a1**3*a2**2*d2-5.76d2*a1**4*a2*d2+1.65d2*a2**4+6.6d2*a1*a2**3 &
+9.9d2*a1**2*a2**2+6.6d2*a1**3*a2+1.65d2*a1**4)
                case (-1)
                  rlYlm_laplacian = &
-3.40542137721830949d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-9.6d1*a1**3*a2**3*d2**2*dy**2+8.8d2*a1**2*a2**3*d2*dy**2 &
+8.8d2*a1**3*a2**2*d2*dy**2-4.8d2*a1*a2**3*dy**2 &
-9.6d2*a1**2*a2**2*dy**2-4.8d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-9.6d1*a1**3*a2**3*d2**2*dx**2 &
+8.8d2*a1**2*a2**3*d2*dx**2+8.8d2*a1**3*a2**2*d2*dx**2 &
-4.8d2*a1*a2**3*dx**2-9.6d2*a1**2*a2**2*dx**2-4.8d2*a1**3*a2*dx**2 &
+3.2d1*a1**3*a2**3*d2**3-3.36d2*a1**2*a2**3*d2**2 &
-3.36d2*a1**3*a2**2*d2**2+5.28d2*a1*a2**3*d2+1.056d3*a1**2*a2**2*d2 &
+5.28d2*a1**3*a2*d2-3.12d2*a2**3-9.36d2*a1*a2**2-9.36d2*a1**2*a2 &
-3.12d2*a1**3)*dz
                case (0)
                  rlYlm_laplacian =5.3844439723186478d-1*E*a1**2*a2*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*dx* &
(3.5d2*a1**4*a2**4*d2*dy**6-2.975d3*a1**3*a2**4*dy**6 &
-2.975d3*a1**4*a2**3*dy**6+1.05d3*a1**4*a2**4*d2*dx**2*dy**4 &
-8.925d3*a1**3*a2**4*dx**2*dy**4-8.925d3*a1**4*a2**3*dx**2*dy**4 &
-6.8d2*a1**4*a2**4*d2**2*dy**4+6.0d3*a1**3*a2**4*d2*dy**4 &
+6.0d3*a1**4*a2**3*d2*dy**4-1.65d3*a1**2*a2**4*dy**4 &
-3.3d3*a1**3*a2**3*dy**4-1.65d3*a1**4*a2**2*dy**4 &
+1.05d3*a1**4*a2**4*d2*dx**4*dy**2-8.925d3*a1**3*a2**4*dx**4*dy**2 &
-8.925d3*a1**4*a2**3*dx**4*dy**2-1.36d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.2d4*a1**3*a2**4*d2*dx**2*dy**2+1.2d4*a1**4*a2**3*d2*dx**2*dy**2 &
-3.3d3*a1**2*a2**4*dx**2*dy**2-6.6d3*a1**3*a2**3*dx**2*dy**2 &
-3.3d3*a1**4*a2**2*dx**2*dy**2+4.0d2*a1**4*a2**4*d2**3*dy**2 &
-3.72d3*a1**3*a2**4*d2**2*dy**2-3.72d3*a1**4*a2**3*d2**2*dy**2 &
+2.46d3*a1**2*a2**4*d2*dy**2+4.92d3*a1**3*a2**3*d2*dy**2 &
+2.46d3*a1**4*a2**2*d2*dy**2-3.9d2*a1*a2**4*dy**2 &
-1.17d3*a1**2*a2**3*dy**2-1.17d3*a1**3*a2**2*dy**2 &
-3.9d2*a1**4*a2*dy**2+3.5d2*a1**4*a2**4*d2*dx**6 &
-2.975d3*a1**3*a2**4*dx**6-2.975d3*a1**4*a2**3*dx**6 &
-6.8d2*a1**4*a2**4*d2**2*dx**4+6.0d3*a1**3*a2**4*d2*dx**4 &
+6.0d3*a1**4*a2**3*d2*dx**4-1.65d3*a1**2*a2**4*dx**4 &
-3.3d3*a1**3*a2**3*dx**4-1.65d3*a1**4*a2**2*dx**4 &
+4.0d2*a1**4*a2**4*d2**3*dx**2-3.72d3*a1**3*a2**4*d2**2*dx**2 &
-3.72d3*a1**4*a2**3*d2**2*dx**2+2.46d3*a1**2*a2**4*d2*dx**2 &
+4.92d3*a1**3*a2**3*d2*dx**2+2.46d3*a1**4*a2**2*d2*dx**2 &
-3.9d2*a1*a2**4*dx**2-1.17d3*a1**2*a2**3*dx**2 &
-1.17d3*a1**3*a2**2*dx**2-3.9d2*a1**4*a2*dx**2-6.4d1*a1**4*a2**4*d2**4 &
+6.08d2*a1**3*a2**4*d2**3+6.08d2*a1**4*a2**3*d2**3 &
-3.84d2*a1**2*a2**4*d2**2-7.68d2*a1**3*a2**3*d2**2 &
-3.84d2*a1**4*a2**2*d2**2-7.44d2*a1*a2**4*d2-2.232d3*a1**2*a2**3*d2 &
-2.232d3*a1**3*a2**2*d2-7.44d2*a1**4*a2*d2+6.6d2*a2**4+2.64d3*a1*a2**3 &
+3.96d3*a1**2*a2**2+2.64d3*a1**3*a2+6.6d2*a1**4)
                case (1)
                  rlYlm_laplacian = &
-1.70271068860915474d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(1.4d2*a1**4*a2**4*d2*dx**2*dy**4-1.19d3*a1**3*a2**4*dx**2*dy**4 &
-1.19d3*a1**4*a2**3*dx**2*dy**4-7.0d1*a1**3*a2**4*d2*dy**4 &
-7.0d1*a1**4*a2**3*d2*dy**4+5.25d2*a1**2*a2**4*dy**4 &
+1.05d3*a1**3*a2**3*dy**4+5.25d2*a1**4*a2**2*dy**4 &
+2.8d2*a1**4*a2**4*d2*dx**4*dy**2-2.38d3*a1**3*a2**4*dx**4*dy**2 &
-2.38d3*a1**4*a2**3*dx**4*dy**2-1.92d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.62d3*a1**3*a2**4*d2*dx**2*dy**2+1.62d3*a1**4*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**4*dx**2*dy**2+1.8d2*a1**3*a2**3*dx**2*dy**2 &
+9.0d1*a1**4*a2**2*dx**2*dy**2+9.6d1*a1**3*a2**4*d2**2*dy**2 &
+9.6d1*a1**4*a2**3*d2**2*dy**2-8.28d2*a1**2*a2**4*d2*dy**2 &
-1.656d3*a1**3*a2**3*d2*dy**2-8.28d2*a1**4*a2**2*d2*dy**2 &
+7.02d2*a1*a2**4*dy**2+2.106d3*a1**2*a2**3*dy**2 &
+2.106d3*a1**3*a2**2*dy**2+7.02d2*a1**4*a2*dy**2 &
+1.4d2*a1**4*a2**4*d2*dx**6-1.19d3*a1**3*a2**4*dx**6 &
-1.19d3*a1**4*a2**3*dx**6-1.92d2*a1**4*a2**4*d2**2*dx**4 &
+1.69d3*a1**3*a2**4*d2*dx**4+1.69d3*a1**4*a2**3*d2*dx**4 &
-4.35d2*a1**2*a2**4*dx**4-8.7d2*a1**3*a2**3*dx**4 &
-4.35d2*a1**4*a2**2*dx**4+6.4d1*a1**4*a2**4*d2**3*dx**2 &
-5.76d2*a1**3*a2**4*d2**2*dx**2-5.76d2*a1**4*a2**3*d2**2*dx**2 &
+2.28d2*a1**2*a2**4*d2*dx**2+4.56d2*a1**3*a2**3*d2*dx**2 &
+2.28d2*a1**4*a2**2*d2*dx**2+7.8d1*a1*a2**4*dx**2 &
+2.34d2*a1**2*a2**3*dx**2+2.34d2*a1**3*a2**2*dx**2 &
+7.8d1*a1**4*a2*dx**2-3.2d1*a1**3*a2**4*d2**3-3.2d1*a1**4*a2**3*d2**3 &
+3.36d2*a1**2*a2**4*d2**2+6.72d2*a1**3*a2**3*d2**2 &
+3.36d2*a1**4*a2**2*d2**2-6.84d2*a1*a2**4*d2-2.052d3*a1**2*a2**3*d2 &
-2.052d3*a1**3*a2**2*d2-6.84d2*a1**4*a2*d2+3.3d2*a2**4+1.32d3*a1*a2**3 &
+1.98d3*a1**2*a2**2+1.32d3*a1**3*a2+3.3d2*a1**4)*dz
                case (2)
                  rlYlm_laplacian =2.40799654862869848d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(7.0d1*a1**4*a2**4*d2*dy**6 &
-5.95d2*a1**3*a2**4*dy**6-5.95d2*a1**4*a2**3*dy**6 &
+7.0d1*a1**4*a2**4*d2*dx**2*dy**4-5.95d2*a1**3*a2**4*dx**2*dy**4 &
-5.95d2*a1**4*a2**3*dx**2*dy**4-1.16d2*a1**4*a2**4*d2**2*dy**4 &
+1.11d3*a1**3*a2**4*d2*dy**4+1.11d3*a1**4*a2**3*d2*dy**4 &
-9.3d2*a1**2*a2**4*dy**4-1.86d3*a1**3*a2**3*dy**4 &
-9.3d2*a1**4*a2**2*dy**4-7.0d1*a1**4*a2**4*d2*dx**4*dy**2 &
+5.95d2*a1**3*a2**4*dx**4*dy**2+5.95d2*a1**4*a2**3*dx**4*dy**2 &
+1.4d2*a1**3*a2**4*d2*dx**2*dy**2+1.4d2*a1**4*a2**3*d2*dx**2*dy**2 &
-1.05d3*a1**2*a2**4*dx**2*dy**2-2.1d3*a1**3*a2**3*dx**2*dy**2 &
-1.05d3*a1**4*a2**2*dx**2*dy**2+4.8d1*a1**4*a2**4*d2**3*dy**2 &
-5.88d2*a1**3*a2**4*d2**2*dy**2-5.88d2*a1**4*a2**3*d2**2*dy**2 &
+1.464d3*a1**2*a2**4*d2*dy**2+2.928d3*a1**3*a2**3*d2*dy**2 &
+1.464d3*a1**4*a2**2*d2*dy**2-7.41d2*a1*a2**4*dy**2 &
-2.223d3*a1**2*a2**3*dy**2-2.223d3*a1**3*a2**2*dy**2 &
-7.41d2*a1**4*a2*dy**2-7.0d1*a1**4*a2**4*d2*dx**6 &
+5.95d2*a1**3*a2**4*dx**6+5.95d2*a1**4*a2**3*dx**6 &
+1.16d2*a1**4*a2**4*d2**2*dx**4-9.7d2*a1**3*a2**4*d2*dx**4 &
-9.7d2*a1**4*a2**3*d2*dx**4-1.2d2*a1**2*a2**4*dx**4 &
-2.4d2*a1**3*a2**3*dx**4-1.2d2*a1**4*a2**2*dx**4 &
-4.8d1*a1**4*a2**4*d2**3*dx**2+3.56d2*a1**3*a2**4*d2**2*dx**2 &
+3.56d2*a1**4*a2**3*d2**2*dx**2+4.32d2*a1**2*a2**4*d2*dx**2 &
+8.64d2*a1**3*a2**3*d2*dx**2+4.32d2*a1**4*a2**2*d2*dx**2 &
-2.73d2*a1*a2**4*dx**2-8.19d2*a1**2*a2**3*dx**2 &
-8.19d2*a1**3*a2**2*dx**2-2.73d2*a1**4*a2*dx**2 &
+4.8d1*a1**3*a2**4*d2**3+4.8d1*a1**4*a2**3*d2**3 &
-4.44d2*a1**2*a2**4*d2**2-8.88d2*a1**3*a2**3*d2**2 &
-4.44d2*a1**4*a2**2*d2**2+5.76d2*a1*a2**4*d2+1.728d3*a1**2*a2**3*d2 &
+1.728d3*a1**3*a2**2*d2+5.76d2*a1**4*a2*d2-1.65d2*a2**4-6.6d2*a1*a2**3 &
-9.9d2*a1**2*a2**2-6.6d2*a1**3*a2-1.65d2*a1**4)
                case (3)
                  rlYlm_laplacian = &
-4.50494903675136302d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(6.0d1*a1**3*a2**3*d2*dx**2*dy**4-5.1d2*a1**2*a2**3*dx**2*dy**4 &
-5.1d2*a1**3*a2**2*dx**2*dy**4-3.0d1*a1**2*a2**3*d2*dy**4 &
-3.0d1*a1**3*a2**2*d2*dy**4+2.25d2*a1*a2**3*dy**4 &
+4.5d2*a1**2*a2**2*dy**4+2.25d2*a1**3*a2*dy**4 &
+4.0d1*a1**3*a2**3*d2*dx**4*dy**2-3.4d2*a1**2*a2**3*dx**4*dy**2 &
-3.4d2*a1**3*a2**2*dx**4*dy**2-4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+4.2d2*a1**2*a2**3*d2*dx**2*dy**2+4.2d2*a1**3*a2**2*d2*dx**2*dy**2 &
-9.0d1*a1*a2**3*dx**2*dy**2-1.8d2*a1**2*a2**2*dx**2*dy**2 &
-9.0d1*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.92d2*a1*a2**3*d2*dy**2 &
-3.84d2*a1**2*a2**2*d2*dy**2-1.92d2*a1**3*a2*d2*dy**2 &
+7.8d1*a2**3*dy**2+2.34d2*a1*a2**2*dy**2+2.34d2*a1**2*a2*dy**2 &
+7.8d1*a1**3*dy**2-2.0d1*a1**3*a2**3*d2*dx**6+1.7d2*a1**2*a2**3*dx**6 &
+1.7d2*a1**3*a2**2*dx**6+1.6d1*a1**3*a2**3*d2**2*dx**4 &
-1.1d2*a1**2*a2**3*d2*dx**4-1.1d2*a1**3*a2**2*d2*dx**4 &
-1.95d2*a1*a2**3*dx**4-3.9d2*a1**2*a2**2*dx**4-1.95d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.92d2*a1*a2**3*d2*dx**2+3.84d2*a1**2*a2**2*d2*dx**2 &
+1.92d2*a1**3*a2*d2*dx**2-7.8d1*a2**3*dx**2-2.34d2*a1*a2**2*dx**2 &
-2.34d2*a1**2*a2*dx**2-7.8d1*a1**3*dx**2)*dz
                case (4)
                  rlYlm_laplacian =3.18548001278669409d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(1.0d1*a1**3*a2**3*d2*dy**6 &
-8.5d1*a1**2*a2**3*dy**6-8.5d1*a1**3*a2**2*dy**6 &
-5.0d1*a1**3*a2**3*d2*dx**2*dy**4+4.25d2*a1**2*a2**3*dx**2*dy**4 &
+4.25d2*a1**3*a2**2*dx**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**4 &
+1.2d2*a1**2*a2**3*d2*dy**4+1.2d2*a1**3*a2**2*d2*dy**4 &
-3.9d2*a1*a2**3*dy**4-7.8d2*a1**2*a2**2*dy**4-3.9d2*a1**3*a2*dy**4 &
-5.0d1*a1**3*a2**3*d2*dx**4*dy**2+4.25d2*a1**2*a2**3*dx**4*dy**2 &
+4.25d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-3.2d2*a1**2*a2**3*d2*dx**2*dy**2-3.2d2*a1**3*a2**2*d2*dx**2*dy**2 &
-6.6d2*a1*a2**3*dx**2*dy**2-1.32d3*a1**2*a2**2*dx**2*dy**2 &
-6.6d2*a1**3*a2*dx**2*dy**2-4.8d1*a1**2*a2**3*d2**2*dy**2 &
-4.8d1*a1**3*a2**2*d2**2*dy**2+3.24d2*a1*a2**3*d2*dy**2 &
+6.48d2*a1**2*a2**2*d2*dy**2+3.24d2*a1**3*a2*d2*dy**2 &
+2.34d2*a2**3*dy**2+7.02d2*a1*a2**2*dy**2+7.02d2*a1**2*a2*dy**2 &
+2.34d2*a1**3*dy**2+1.0d1*a1**3*a2**3*d2*dx**6-8.5d1*a1**2*a2**3*dx**6 &
-8.5d1*a1**3*a2**2*dx**6-8.0d0*a1**3*a2**3*d2**2*dx**4 &
+4.0d1*a1**2*a2**3*d2*dx**4+4.0d1*a1**3*a2**2*d2*dx**4 &
+2.1d2*a1*a2**3*dx**4+4.2d2*a1**2*a2**2*dx**4+2.1d2*a1**3*a2*dx**4 &
+1.6d1*a1**2*a2**3*d2**2*dx**2+1.6d1*a1**3*a2**2*d2**2*dx**2 &
-1.08d2*a1*a2**3*d2*dx**2-2.16d2*a1**2*a2**2*d2*dx**2 &
-1.08d2*a1**3*a2*d2*dx**2-7.8d1*a2**3*dx**2-2.34d2*a1*a2**2*dx**2 &
-2.34d2*a1**2*a2*dx**2-7.8d1*a1**3*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (2)
          ! selection on l2: l1=3, m1=2
          select case (l2)
            case (0)
              ! selection on m2: l1=3, m1=2, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =4.54056183629107931d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*(dy+dx)* &
(dy-1.0d0*dx)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=3, m1=2, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =7.86448379536438834d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dy*(2.0d0*a1**2*a2**2*d2*dy**2 &
-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2-2.0d0*a1*a2**2*d2 &
-2.0d0*a1**2*a2*d2+9.0d0*a2**2+1.8d1*a1*a2+9.0d0*a1**2)*dz
                case (0)
                  rlYlm_laplacian = &
-3.93224189768219417d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-8)*(dy+dx)* &
(dy-1.0d0*dx)*(4.0d0*a1**2*a2**2*d2*dy**2-2.2d1*a1*a2**2*dy**2 &
-2.2d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2-2.2d1*a1*a2**2*dx**2 &
-2.2d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+2.4d1*a1*a2**2*d2 &
+2.4d1*a1**2*a2*d2-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)
                case (1)
                  rlYlm_laplacian =7.86448379536438834d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*(2.0d0*a1**2*a2**2*d2*dy**2 &
-1.1d1*a1*a2**2*dy**2-1.1d1*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2+2.0d0*a1*a2**2*d2 &
+2.0d0*a1**2*a2*d2-9.0d0*a2**2-1.8d1*a1*a2-9.0d0*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=3, m1=2, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =1.75855203743803178d1*E*a1**3*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*dz
                case (-1)
                  rlYlm_laplacian = &
-8.79276018719015889d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(4.0d0*a1**3*a2**3*d2*dy**4-2.6d1*a1**2*a2**3*dy**4 &
-2.6d1*a1**3*a2**2*dy**4-4.0d0*a1**3*a2**3*d2**2*dy**2 &
+2.4d1*a1**2*a2**3*d2*dy**2+2.4d1*a1**3*a2**2*d2*dy**2 &
+1.1d1*a1*a2**3*dy**2+2.2d1*a1**2*a2**2*dy**2+1.1d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+2.6d1*a1**2*a2**3*dx**4 &
+2.6d1*a1**3*a2**2*dx**4+4.0d0*a1**3*a2**3*d2**2*dx**2 &
-3.2d1*a1**2*a2**3*d2*dx**2-3.2d1*a1**3*a2**2*d2*dx**2 &
+3.3d1*a1*a2**3*dx**2+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2 &
+4.0d0*a1**2*a2**3*d2**2+4.0d0*a1**3*a2**2*d2**2-2.4d1*a1*a2**3*d2 &
-4.8d1*a1**2*a2**2*d2-2.4d1*a1**3*a2*d2+9.0d0*a2**3+2.7d1*a1*a2**2 &
+2.7d1*a1**2*a2+9.0d0*a1**3)
                case (0)
                  rlYlm_laplacian = &
-5.07650246099406246d0*E*a1**3*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)* &
(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*(dy+dx)*(dy-1.0d0*dx)*(3.0d0* &
(dy**2+dx**2)-2.0d0*d2)*dz
                case (1)
                  rlYlm_laplacian = &
-8.79276018719015889d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(4.0d0*a1**3*a2**3*d2*dy**4-2.6d1*a1**2*a2**3*dy**4 &
-2.6d1*a1**3*a2**2*dy**4-4.0d0*a1**3*a2**3*d2**2*dy**2 &
+3.2d1*a1**2*a2**3*d2*dy**2+3.2d1*a1**3*a2**2*d2*dy**2 &
-3.3d1*a1*a2**3*dy**2-6.6d1*a1**2*a2**2*dy**2-3.3d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+2.6d1*a1**2*a2**3*dx**4 &
+2.6d1*a1**3*a2**2*dx**4+4.0d0*a1**3*a2**3*d2**2*dx**2 &
-2.4d1*a1**2*a2**3*d2*dx**2-2.4d1*a1**3*a2**2*d2*dx**2 &
-1.1d1*a1*a2**3*dx**2-2.2d1*a1**2*a2**2*dx**2-1.1d1*a1**3*a2*dx**2 &
-4.0d0*a1**2*a2**3*d2**2-4.0d0*a1**3*a2**2*d2**2+2.4d1*a1*a2**3*d2 &
+4.8d1*a1**2*a2**2*d2+2.4d1*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)
                case (2)
                  rlYlm_laplacian = &
-8.79276018719015889d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)* &
(2.0d0*a1**3*a2**3*d2*dy**4-1.3d1*a1**2*a2**3*dy**4 &
-1.3d1*a1**3*a2**2*dy**4-4.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
+2.6d1*a1**2*a2**3*dx**2*dy**2+2.6d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+2.2d1*a1*a2**3*dy**2+4.4d1*a1**2*a2**2*dy**2+2.2d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.3d1*a1**2*a2**3*dx**4 &
-1.3d1*a1**3*a2**2*dx**4-4.0d0*a1**2*a2**3*d2*dx**2 &
-4.0d0*a1**3*a2**2*d2*dx**2+2.2d1*a1*a2**3*dx**2 &
+4.4d1*a1**2*a2**2*dx**2+2.2d1*a1**3*a2*dx**2+2.0d0*a1*a2**3*d2 &
+4.0d0*a1**2*a2**2*d2+2.0d0*a1**3*a2*d2-9.0d0*a2**3-2.7d1*a1*a2**2 &
-2.7d1*a1**2*a2-9.0d0*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=3, m1=2, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-9.49726646607726303d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(2.0d0*a1**3*a2**3*d2*dy**4-1.5d1*a1**2*a2**3*dy**4 &
-1.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
+6.0d1*a1**2*a2**3*dx**2*dy**2+6.0d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+6.0d0*a1**3*a2**3*d2*dx**4-4.5d1*a1**2*a2**3*dx**4 &
-4.5d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (-2)
                  rlYlm_laplacian = &
-2.32634567931348979d1*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(4.0d0*a1**2*a2**2*d2*dy**2-3.0d1*a1*a2**2*dy**2 &
-3.0d1*a1**2*a2*dy**2+4.0d0*a1**2*a2**2*d2*dx**2-3.0d1*a1*a2**2*dx**2 &
-3.0d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2+3.2d1*a1*a2**2*d2 &
+3.2d1*a1**2*a2*d2-1.3d1*a2**2-2.6d1*a1*a2-1.3d1*a1**2)
                case (-1)
                  rlYlm_laplacian = &
-7.35655097152228142d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(1.0d1*a1**3*a2**3*d2*dy**4-7.5d1*a1**2*a2**3*dy**4 &
-7.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+5.4d1*a1**2*a2**3*d2*dy**2+5.4d1*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+7.5d1*a1**2*a2**3*dx**4 &
+7.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-7.4d1*a1**2*a2**3*d2*dx**2-7.4d1*a1**3*a2**2*d2*dx**2 &
+9.1d1*a1*a2**3*dx**2+1.82d2*a1**2*a2**2*dx**2+9.1d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-5.8d1*a1*a2**3*d2 &
-1.16d2*a1**2*a2**2*d2-5.8d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian =3.00329935783424201d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(dy+dx)*(dy-1.0d0*dx)* &
(2.0d1*a1**3*a2**3*d2*dy**4-1.5d2*a1**2*a2**3*dy**4 &
-1.5d2*a1**3*a2**2*dy**4+4.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-3.0d2*a1**2*a2**3*dx**2*dy**2-3.0d2*a1**3*a2**2*dx**2*dy**2 &
-2.8d1*a1**3*a2**3*d2**2*dy**2+2.04d2*a1**2*a2**3*d2*dy**2 &
+2.04d2*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.5d2*a1**2*a2**3*dx**4 &
-1.5d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+2.04d2*a1**2*a2**3*d2*dx**2+2.04d2*a1**3*a2**2*d2*dx**2 &
+3.9d1*a1*a2**3*dx**2+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-4.8d1*a1**2*a2**3*d2**2 &
-4.8d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2-1.8d2*a1**2*a2**2*d2 &
-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2+1.98d2*a1**2*a2 &
+6.6d1*a1**3)
                case (1)
                  rlYlm_laplacian = &
-7.35655097152228142d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(1.0d1*a1**3*a2**3*d2*dy**4-7.5d1*a1**2*a2**3*dy**4 &
-7.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+7.4d1*a1**2*a2**3*d2*dy**2+7.4d1*a1**3*a2**2*d2*dy**2 &
-9.1d1*a1*a2**3*dy**2-1.82d2*a1**2*a2**2*dy**2-9.1d1*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+7.5d1*a1**2*a2**3*dx**4 &
+7.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-5.4d1*a1**2*a2**3*d2*dx**2-5.4d1*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
-8.0d0*a1**2*a2**3*d2**2-8.0d0*a1**3*a2**2*d2**2+5.8d1*a1*a2**3*d2 &
+1.16d2*a1**2*a2**2*d2+5.8d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian =1.16317283965674489d1*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(4.0d0*a1**4*a2**4*d2*dy**6 &
-3.0d1*a1**3*a2**4*dy**6-3.0d1*a1**4*a2**3*dy**6 &
-4.0d0*a1**4*a2**4*d2*dx**2*dy**4+3.0d1*a1**3*a2**4*dx**2*dy**4 &
+3.0d1*a1**4*a2**3*dx**2*dy**4-4.0d0*a1**4*a2**4*d2**2*dy**4 &
+2.4d1*a1**3*a2**4*d2*dy**4+2.4d1*a1**4*a2**3*d2*dy**4 &
+3.9d1*a1**2*a2**4*dy**4+7.8d1*a1**3*a2**3*dy**4 &
+3.9d1*a1**4*a2**2*dy**4-4.0d0*a1**4*a2**4*d2*dx**4*dy**2 &
+3.0d1*a1**3*a2**4*dx**4*dy**2+3.0d1*a1**4*a2**3*dx**4*dy**2 &
+8.0d0*a1**4*a2**4*d2**2*dx**2*dy**2-8.0d1*a1**3*a2**4*d2*dx**2*dy**2 &
-8.0d1*a1**4*a2**3*d2*dx**2*dy**2+1.3d2*a1**2*a2**4*dx**2*dy**2 &
+2.6d2*a1**3*a2**3*dx**2*dy**2+1.3d2*a1**4*a2**2*dx**2*dy**2 &
+8.0d0*a1**3*a2**4*d2**2*dy**2+8.0d0*a1**4*a2**3*d2**2*dy**2 &
-5.2d1*a1**2*a2**4*d2*dy**2-1.04d2*a1**3*a2**3*d2*dy**2 &
-5.2d1*a1**4*a2**2*d2*dy**2+4.0d0*a1**4*a2**4*d2*dx**6 &
-3.0d1*a1**3*a2**4*dx**6-3.0d1*a1**4*a2**3*dx**6 &
-4.0d0*a1**4*a2**4*d2**2*dx**4+2.4d1*a1**3*a2**4*d2*dx**4 &
+2.4d1*a1**4*a2**3*d2*dx**4+3.9d1*a1**2*a2**4*dx**4 &
+7.8d1*a1**3*a2**3*dx**4+3.9d1*a1**4*a2**2*dx**4 &
+8.0d0*a1**3*a2**4*d2**2*dx**2+8.0d0*a1**4*a2**3*d2**2*dx**2 &
-5.2d1*a1**2*a2**4*d2*dx**2-1.04d2*a1**3*a2**3*d2*dx**2 &
-5.2d1*a1**4*a2**2*d2*dx**2-4.0d0*a1**2*a2**4*d2**2 &
-8.0d0*a1**3*a2**3*d2**2-4.0d0*a1**4*a2**2*d2**2+2.4d1*a1*a2**4*d2 &
+7.2d1*a1**2*a2**3*d2+7.2d1*a1**3*a2**2*d2+2.4d1*a1**4*a2*d2 &
-9.0d0*a2**4-3.6d1*a1*a2**3-5.4d1*a1**2*a2**2-3.6d1*a1**3*a2 &
-9.0d0*a1**4)
                case (3)
                  rlYlm_laplacian = &
-9.49726646607726303d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(6.0d0*a1**3*a2**3*d2*dy**4-4.5d1*a1**2*a2**3*dy**4 &
-4.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
+6.0d1*a1**2*a2**3*dx**2*dy**2+6.0d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=3, m1=2, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-4.02934891253929843d1*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(2.0d0*a1**3*a2**3*d2*dy**4-1.7d1*a1**2*a2**3*dy**4 &
-1.7d1*a1**3*a2**2*dy**4-4.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
+3.4d1*a1**2*a2**3*dx**2*dy**2+3.4d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+3.0d1*a1*a2**3*dy**2+6.0d1*a1**2*a2**2*dy**2+3.0d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.7d1*a1**2*a2**3*dx**4 &
-1.7d1*a1**3*a2**2*dx**4-4.0d0*a1**2*a2**3*d2*dx**2 &
-4.0d0*a1**3*a2**2*d2*dx**2+3.0d1*a1*a2**3*dx**2 &
+6.0d1*a1**2*a2**2*dx**2+3.0d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.9d1*a2**3-1.17d2*a1*a2**2 &
-1.17d2*a1**2*a2-3.9d1*a1**3)*dz
                case (-3)
                  rlYlm_laplacian =1.42458996991158946d1*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(4.0d0*a1**4*a2**4*d2*dy**6 &
-3.4d1*a1**3*a2**4*dy**6-3.4d1*a1**4*a2**3*dy**6 &
-1.2d1*a1**4*a2**4*d2*dx**2*dy**4+1.02d2*a1**3*a2**4*dx**2*dy**4 &
+1.02d2*a1**4*a2**3*dx**2*dy**4-4.0d0*a1**4*a2**4*d2**2*dy**4 &
+2.4d1*a1**3*a2**4*d2*dy**4+2.4d1*a1**4*a2**3*d2*dy**4 &
+7.5d1*a1**2*a2**4*dy**4+1.5d2*a1**3*a2**3*dy**4 &
+7.5d1*a1**4*a2**2*dy**4-4.0d0*a1**4*a2**4*d2*dx**4*dy**2 &
+3.4d1*a1**3*a2**4*dx**4*dy**2+3.4d1*a1**4*a2**3*dx**4*dy**2 &
+1.6d1*a1**4*a2**4*d2**2*dx**2*dy**2-1.68d2*a1**3*a2**4*d2*dx**2*dy**2 &
-1.68d2*a1**4*a2**3*d2*dx**2*dy**2+2.4d2*a1**2*a2**4*dx**2*dy**2 &
+4.8d2*a1**3*a2**3*dx**2*dy**2+2.4d2*a1**4*a2**2*dx**2*dy**2 &
+1.2d1*a1**3*a2**4*d2**2*dy**2+1.2d1*a1**4*a2**3*d2**2*dy**2 &
-8.4d1*a1**2*a2**4*d2*dy**2-1.68d2*a1**3*a2**3*d2*dy**2 &
-8.4d1*a1**4*a2**2*d2*dy**2-3.9d1*a1*a2**4*dy**2 &
-1.17d2*a1**2*a2**3*dy**2-1.17d2*a1**3*a2**2*dy**2 &
-3.9d1*a1**4*a2*dy**2+1.2d1*a1**4*a2**4*d2*dx**6 &
-1.02d2*a1**3*a2**4*dx**6-1.02d2*a1**4*a2**3*dx**6 &
-1.2d1*a1**4*a2**4*d2**2*dx**4+9.6d1*a1**3*a2**4*d2*dx**4 &
+9.6d1*a1**4*a2**3*d2*dx**4+4.5d1*a1**2*a2**4*dx**4 &
+9.0d1*a1**3*a2**3*dx**4+4.5d1*a1**4*a2**2*dx**4 &
+1.2d1*a1**3*a2**4*d2**2*dx**2+1.2d1*a1**4*a2**3*d2**2*dx**2 &
-8.4d1*a1**2*a2**4*d2*dx**2-1.68d2*a1**3*a2**3*d2*dx**2 &
-8.4d1*a1**4*a2**2*d2*dx**2-3.9d1*a1*a2**4*dx**2 &
-1.17d2*a1**2*a2**3*dx**2-1.17d2*a1**3*a2**2*dx**2 &
-3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2-2.4d1*a1**3*a2**3*d2**2 &
-1.2d1*a1**4*a2**2*d2**2+8.4d1*a1*a2**4*d2+2.52d2*a1**2*a2**3*d2 &
+2.52d2*a1**3*a2**2*d2+8.4d1*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (-2)
                  rlYlm_laplacian = &
-1.52295073829821874d1*E*a1**4*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(1.4d1*a1**2*a2**2*d2*dy**2 &
-1.19d2*a1*a2**2*dy**2-1.19d2*a1**2*a2*dy**2 &
+1.4d1*a1**2*a2**2*d2*dx**2-1.19d2*a1*a2**2*dx**2 &
-1.19d2*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2+1.1d2*a1*a2**2*d2 &
+1.1d2*a1**2*a2*d2-6.0d1*a2**2-1.2d2*a1*a2-6.0d1*a1**2)*dz
                case (-1)
                  rlYlm_laplacian =5.3844439723186478d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(2.8d1*a1**4*a2**4*d2*dy**6 &
-2.38d2*a1**3*a2**4*dy**6-2.38d2*a1**4*a2**3*dy**6 &
+2.8d1*a1**4*a2**4*d2*dx**2*dy**4-2.38d2*a1**3*a2**4*dx**2*dy**4 &
-2.38d2*a1**4*a2**3*dx**2*dy**4-4.4d1*a1**4*a2**4*d2**2*dy**4 &
+3.52d2*a1**3*a2**4*d2*dy**4+3.52d2*a1**4*a2**3*d2*dy**4 &
+1.65d2*a1**2*a2**4*dy**4+3.3d2*a1**3*a2**3*dy**4 &
+1.65d2*a1**4*a2**2*dy**4-2.8d1*a1**4*a2**4*d2*dx**4*dy**2 &
+2.38d2*a1**3*a2**4*dx**4*dy**2+2.38d2*a1**4*a2**3*dx**4*dy**2 &
-5.6d1*a1**3*a2**4*d2*dx**2*dy**2-5.6d1*a1**4*a2**3*d2*dx**2*dy**2 &
+4.2d2*a1**2*a2**4*dx**2*dy**2+8.4d2*a1**3*a2**3*dx**2*dy**2 &
+4.2d2*a1**4*a2**2*dx**2*dy**2+1.6d1*a1**4*a2**4*d2**3*dy**2 &
-9.2d1*a1**3*a2**4*d2**2*dy**2-9.2d1*a1**4*a2**3*d2**2*dy**2 &
-3.6d2*a1**2*a2**4*d2*dy**2-7.2d2*a1**3*a2**3*d2*dy**2 &
-3.6d2*a1**4*a2**2*d2*dy**2+1.95d2*a1*a2**4*dy**2 &
+5.85d2*a1**2*a2**3*dy**2+5.85d2*a1**3*a2**2*dy**2 &
+1.95d2*a1**4*a2*dy**2-2.8d1*a1**4*a2**4*d2*dx**6 &
+2.38d2*a1**3*a2**4*dx**6+2.38d2*a1**4*a2**3*dx**6 &
+4.4d1*a1**4*a2**4*d2**2*dx**4-4.08d2*a1**3*a2**4*d2*dx**4 &
-4.08d2*a1**4*a2**3*d2*dx**4+2.55d2*a1**2*a2**4*dx**4 &
+5.1d2*a1**3*a2**3*dx**4+2.55d2*a1**4*a2**2*dx**4 &
-1.6d1*a1**4*a2**4*d2**3*dx**2+1.8d2*a1**3*a2**4*d2**2*dx**2 &
+1.8d2*a1**4*a2**3*d2**2*dx**2-3.36d2*a1**2*a2**4*d2*dx**2 &
-6.72d2*a1**3*a2**3*d2*dx**2-3.36d2*a1**4*a2**2*d2*dx**2 &
+3.9d1*a1*a2**4*dx**2+1.17d2*a1**2*a2**3*dx**2 &
+1.17d2*a1**3*a2**2*dx**2+3.9d1*a1**4*a2*dx**2-1.6d1*a1**3*a2**4*d2**3 &
-1.6d1*a1**4*a2**3*d2**3+1.32d2*a1**2*a2**4*d2**2 &
+2.64d2*a1**3*a2**3*d2**2+1.32d2*a1**4*a2**2*d2**2-7.2d1*a1*a2**4*d2 &
-2.16d2*a1**2*a2**3*d2-2.16d2*a1**3*a2**2*d2-7.2d1*a1**4*a2*d2 &
-3.3d1*a2**4-1.32d2*a1*a2**3-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2 &
-3.3d1*a1**4)
                case (0)
                  rlYlm_laplacian =1.70271068860915474d0*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(dy+dx)*(dy-1.0d0*dx)* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+6.4d2*a1**2*a2**3*d2*dy**2 &
+6.4d2*a1**3*a2**2*d2*dy**2+3.0d2*a1*a2**3*dy**2 &
+6.0d2*a1**2*a2**2*dy**2+3.0d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.4d2*a1**2*a2**3*d2*dx**2+6.4d2*a1**3*a2**2*d2*dx**2 &
+3.0d2*a1*a2**3*dx**2+6.0d2*a1**2*a2**2*dx**2+3.0d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-7.2d1*a1**2*a2**3*d2**2 &
-7.2d1*a1**3*a2**2*d2**2-5.64d2*a1*a2**3*d2-1.128d3*a1**2*a2**2*d2 &
-5.64d2*a1**3*a2*d2+5.46d2*a2**3+1.638d3*a1*a2**2+1.638d3*a1**2*a2 &
+5.46d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian =5.3844439723186478d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(2.8d1*a1**4*a2**4*d2*dy**6 &
-2.38d2*a1**3*a2**4*dy**6-2.38d2*a1**4*a2**3*dy**6 &
+2.8d1*a1**4*a2**4*d2*dx**2*dy**4-2.38d2*a1**3*a2**4*dx**2*dy**4 &
-2.38d2*a1**4*a2**3*dx**2*dy**4-4.4d1*a1**4*a2**4*d2**2*dy**4 &
+4.08d2*a1**3*a2**4*d2*dy**4+4.08d2*a1**4*a2**3*d2*dy**4 &
-2.55d2*a1**2*a2**4*dy**4-5.1d2*a1**3*a2**3*dy**4 &
-2.55d2*a1**4*a2**2*dy**4-2.8d1*a1**4*a2**4*d2*dx**4*dy**2 &
+2.38d2*a1**3*a2**4*dx**4*dy**2+2.38d2*a1**4*a2**3*dx**4*dy**2 &
+5.6d1*a1**3*a2**4*d2*dx**2*dy**2+5.6d1*a1**4*a2**3*d2*dx**2*dy**2 &
-4.2d2*a1**2*a2**4*dx**2*dy**2-8.4d2*a1**3*a2**3*dx**2*dy**2 &
-4.2d2*a1**4*a2**2*dx**2*dy**2+1.6d1*a1**4*a2**4*d2**3*dy**2 &
-1.8d2*a1**3*a2**4*d2**2*dy**2-1.8d2*a1**4*a2**3*d2**2*dy**2 &
+3.36d2*a1**2*a2**4*d2*dy**2+6.72d2*a1**3*a2**3*d2*dy**2 &
+3.36d2*a1**4*a2**2*d2*dy**2-3.9d1*a1*a2**4*dy**2 &
-1.17d2*a1**2*a2**3*dy**2-1.17d2*a1**3*a2**2*dy**2 &
-3.9d1*a1**4*a2*dy**2-2.8d1*a1**4*a2**4*d2*dx**6 &
+2.38d2*a1**3*a2**4*dx**6+2.38d2*a1**4*a2**3*dx**6 &
+4.4d1*a1**4*a2**4*d2**2*dx**4-3.52d2*a1**3*a2**4*d2*dx**4 &
-3.52d2*a1**4*a2**3*d2*dx**4-1.65d2*a1**2*a2**4*dx**4 &
-3.3d2*a1**3*a2**3*dx**4-1.65d2*a1**4*a2**2*dx**4 &
-1.6d1*a1**4*a2**4*d2**3*dx**2+9.2d1*a1**3*a2**4*d2**2*dx**2 &
+9.2d1*a1**4*a2**3*d2**2*dx**2+3.6d2*a1**2*a2**4*d2*dx**2 &
+7.2d2*a1**3*a2**3*d2*dx**2+3.6d2*a1**4*a2**2*d2*dx**2 &
-1.95d2*a1*a2**4*dx**2-5.85d2*a1**2*a2**3*dx**2 &
-5.85d2*a1**3*a2**2*dx**2-1.95d2*a1**4*a2*dx**2 &
+1.6d1*a1**3*a2**4*d2**3+1.6d1*a1**4*a2**3*d2**3 &
-1.32d2*a1**2*a2**4*d2**2-2.64d2*a1**3*a2**3*d2**2 &
-1.32d2*a1**4*a2**2*d2**2+7.2d1*a1*a2**4*d2+2.16d2*a1**2*a2**3*d2 &
+2.16d2*a1**3*a2**2*d2+7.2d1*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case (2)
                  rlYlm_laplacian =7.6147536914910937d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(1.4d1*a1**4*a2**4*d2*dy**6 &
-1.19d2*a1**3*a2**4*dy**6-1.19d2*a1**4*a2**3*dy**6 &
-1.4d1*a1**4*a2**4*d2*dx**2*dy**4+1.19d2*a1**3*a2**4*dx**2*dy**4 &
+1.19d2*a1**4*a2**3*dx**2*dy**4-1.2d1*a1**4*a2**4*d2**2*dy**4 &
+8.2d1*a1**3*a2**4*d2*dy**4+8.2d1*a1**4*a2**3*d2*dy**4 &
+1.5d2*a1**2*a2**4*dy**4+3.0d2*a1**3*a2**3*dy**4 &
+1.5d2*a1**4*a2**2*dy**4-1.4d1*a1**4*a2**4*d2*dx**4*dy**2 &
+1.19d2*a1**3*a2**4*dx**4*dy**2+1.19d2*a1**4*a2**3*dx**4*dy**2 &
+2.4d1*a1**4*a2**4*d2**2*dx**2*dy**2-2.76d2*a1**3*a2**4*d2*dx**2*dy**2 &
-2.76d2*a1**4*a2**3*d2*dx**2*dy**2+5.4d2*a1**2*a2**4*dx**2*dy**2 &
+1.08d3*a1**3*a2**3*dx**2*dy**2+5.4d2*a1**4*a2**2*dx**2*dy**2 &
+2.4d1*a1**3*a2**4*d2**2*dy**2+2.4d1*a1**4*a2**3*d2**2*dy**2 &
-1.86d2*a1**2*a2**4*d2*dy**2-3.72d2*a1**3*a2**3*d2*dy**2 &
-1.86d2*a1**4*a2**2*d2*dy**2+3.9d1*a1*a2**4*dy**2 &
+1.17d2*a1**2*a2**3*dy**2+1.17d2*a1**3*a2**2*dy**2 &
+3.9d1*a1**4*a2*dy**2+1.4d1*a1**4*a2**4*d2*dx**6 &
-1.19d2*a1**3*a2**4*dx**6-1.19d2*a1**4*a2**3*dx**6 &
-1.2d1*a1**4*a2**4*d2**2*dx**4+8.2d1*a1**3*a2**4*d2*dx**4 &
+8.2d1*a1**4*a2**3*d2*dx**4+1.5d2*a1**2*a2**4*dx**4 &
+3.0d2*a1**3*a2**3*dx**4+1.5d2*a1**4*a2**2*dx**4 &
+2.4d1*a1**3*a2**4*d2**2*dx**2+2.4d1*a1**4*a2**3*d2**2*dx**2 &
-1.86d2*a1**2*a2**4*d2*dx**2-3.72d2*a1**3*a2**3*d2*dx**2 &
-1.86d2*a1**4*a2**2*d2*dx**2+3.9d1*a1*a2**4*dx**2 &
+1.17d2*a1**2*a2**3*dx**2+1.17d2*a1**3*a2**2*dx**2 &
+3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2-2.4d1*a1**3*a2**3*d2**2 &
-1.2d1*a1**4*a2**2*d2**2+9.0d1*a1*a2**4*d2+2.7d2*a1**2*a2**3*d2 &
+2.7d2*a1**3*a2**2*d2+9.0d1*a1**4*a2*d2-6.6d1*a2**4-2.64d2*a1*a2**3 &
-3.96d2*a1**2*a2**2-2.64d2*a1**3*a2-6.6d1*a1**4)*dz
                case (3)
                  rlYlm_laplacian =1.42458996991158946d1*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(1.2d1*a1**4*a2**4*d2*dy**6 &
-1.02d2*a1**3*a2**4*dy**6-1.02d2*a1**4*a2**3*dy**6 &
-4.0d0*a1**4*a2**4*d2*dx**2*dy**4+3.4d1*a1**3*a2**4*dx**2*dy**4 &
+3.4d1*a1**4*a2**3*dx**2*dy**4-1.2d1*a1**4*a2**4*d2**2*dy**4 &
+9.6d1*a1**3*a2**4*d2*dy**4+9.6d1*a1**4*a2**3*d2*dy**4 &
+4.5d1*a1**2*a2**4*dy**4+9.0d1*a1**3*a2**3*dy**4 &
+4.5d1*a1**4*a2**2*dy**4-1.2d1*a1**4*a2**4*d2*dx**4*dy**2 &
+1.02d2*a1**3*a2**4*dx**4*dy**2+1.02d2*a1**4*a2**3*dx**4*dy**2 &
+1.6d1*a1**4*a2**4*d2**2*dx**2*dy**2-1.68d2*a1**3*a2**4*d2*dx**2*dy**2 &
-1.68d2*a1**4*a2**3*d2*dx**2*dy**2+2.4d2*a1**2*a2**4*dx**2*dy**2 &
+4.8d2*a1**3*a2**3*dx**2*dy**2+2.4d2*a1**4*a2**2*dx**2*dy**2 &
+1.2d1*a1**3*a2**4*d2**2*dy**2+1.2d1*a1**4*a2**3*d2**2*dy**2 &
-8.4d1*a1**2*a2**4*d2*dy**2-1.68d2*a1**3*a2**3*d2*dy**2 &
-8.4d1*a1**4*a2**2*d2*dy**2-3.9d1*a1*a2**4*dy**2 &
-1.17d2*a1**2*a2**3*dy**2-1.17d2*a1**3*a2**2*dy**2 &
-3.9d1*a1**4*a2*dy**2+4.0d0*a1**4*a2**4*d2*dx**6 &
-3.4d1*a1**3*a2**4*dx**6-3.4d1*a1**4*a2**3*dx**6 &
-4.0d0*a1**4*a2**4*d2**2*dx**4+2.4d1*a1**3*a2**4*d2*dx**4 &
+2.4d1*a1**4*a2**3*d2*dx**4+7.5d1*a1**2*a2**4*dx**4 &
+1.5d2*a1**3*a2**3*dx**4+7.5d1*a1**4*a2**2*dx**4 &
+1.2d1*a1**3*a2**4*d2**2*dx**2+1.2d1*a1**4*a2**3*d2**2*dx**2 &
-8.4d1*a1**2*a2**4*d2*dx**2-1.68d2*a1**3*a2**3*d2*dx**2 &
-8.4d1*a1**4*a2**2*d2*dx**2-3.9d1*a1*a2**4*dx**2 &
-1.17d2*a1**2*a2**3*dx**2-1.17d2*a1**3*a2**2*dx**2 &
-3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2-2.4d1*a1**3*a2**3*d2**2 &
-1.2d1*a1**4*a2**2*d2**2+8.4d1*a1*a2**4*d2+2.52d2*a1**2*a2**3*d2 &
+2.52d2*a1**3*a2**2*d2+8.4d1*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (4)
                  rlYlm_laplacian =1.00733722813482461d1*E*a1**3*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(dy+dx)*(dy-1.0d0*dx)* &
(2.0d0*a1**3*a2**3*d2*dy**4-1.7d1*a1**2*a2**3*dy**4 &
-1.7d1*a1**3*a2**2*dy**4-1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
+1.02d2*a1**2*a2**3*dx**2*dy**2+1.02d2*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**2*a2**3*d2*dy**2-8.0d0*a1**3*a2**2*d2*dy**2 &
+6.0d1*a1*a2**3*dy**2+1.2d2*a1**2*a2**2*dy**2+6.0d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.7d1*a1**2*a2**3*dx**4 &
-1.7d1*a1**3*a2**2*dx**4-8.0d0*a1**2*a2**3*d2*dx**2 &
-8.0d0*a1**3*a2**2*d2*dx**2+6.0d1*a1*a2**3*dx**2 &
+1.2d2*a1**2*a2**2*dx**2+6.0d1*a1**3*a2*dx**2+1.2d1*a1*a2**3*d2 &
+2.4d1*a1**2*a2**2*d2+1.2d1*a1**3*a2*d2-7.8d1*a2**3-2.34d2*a1*a2**2 &
-2.34d2*a1**2*a2-7.8d1*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (3)
          ! selection on l2: l1=3, m1=3
          select case (l2)
            case (0)
              ! selection on m2: l1=3, m1=3, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =1.85367660741129178d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-7)*(1.0d0*a2*(2.0d0*a1*d2-9.0d0)-9.0d0*a1)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=3, m1=3, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =3.21066206483826464d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*dx*dy*(6.0d0*a1**2*a2**2*d2*dy**2 &
-3.3d1*a1*a2**2*dy**2-3.3d1*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+1.1d1*a1*a2**2*dx**2+1.1d1*a1**2*a2*dx**2-6.0d0*a1*a2**2*d2 &
-6.0d0*a1**2*a2*d2+2.7d1*a2**2+5.4d1*a1*a2+2.7d1*a1**2)
                case (0)
                  rlYlm_laplacian =3.21066206483826464d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*dz
                case (1)
                  rlYlm_laplacian =1.60533103241913232d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.2d1*a1**2*a2**2*d2*dx**2*dy**2 &
-6.6d1*a1*a2**2*dx**2*dy**2-6.6d1*a1**2*a2*dx**2*dy**2 &
-6.0d0*a1*a2**2*d2*dy**2-6.0d0*a1**2*a2*d2*dy**2+2.7d1*a2**2*dy**2 &
+5.4d1*a1*a2*dy**2+2.7d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4 &
+2.2d1*a1*a2**2*dx**4+2.2d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-2.7d1*a2**2*dx**2-5.4d1*a1*a2*dx**2 &
-2.7d1*a1**2*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=3, m1=3, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =3.58962931487909853d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
-7.8d1*a1**2*a2**3*dx**2*dy**2-7.8d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.3d1*a1*a2**3*dy**2+6.6d1*a1**2*a2**2*dy**2+3.3d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+2.6d1*a1**2*a2**3*dx**4 &
+2.6d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.3d1*a1*a2**3*dx**2 &
+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-2.7d1*a2**3-8.1d1*a1*a2**2 &
-8.1d1*a1**2*a2-2.7d1*a1**3)
                case (-1)
                  rlYlm_laplacian =7.17925862975819707d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*dy*(6.0d0*a1**2*a2**2*d2*dy**2 &
-3.9d1*a1*a2**2*dy**2-3.9d1*a1**2*a2*dy**2-2.0d0*a1**2*a2**2*d2*dx**2 &
+1.3d1*a1*a2**2*dx**2+1.3d1*a1**2*a2*dx**2-6.0d0*a1*a2**2*d2 &
-6.0d0*a1**2*a2*d2+3.3d1*a2**2+6.6d1*a1*a2+3.3d1*a1**2)*dz
                case (0)
                  rlYlm_laplacian = &
-2.07247345123641944d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*(6.0d0*a1**2*a2**2*d2*dy**2 &
-3.9d1*a1*a2**2*dy**2-3.9d1*a1**2*a2*dy**2+6.0d0*a1**2*a2**2*d2*dx**2 &
-3.9d1*a1*a2**2*dx**2-3.9d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+2.0d1*a1*a2**2*d2+2.0d1*a1**2*a2*d2+3.3d1*a2**2+6.6d1*a1*a2 &
+3.3d1*a1**2)
                case (1)
                  rlYlm_laplacian =3.58962931487909853d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-9)*(1.2d1*a1**2*a2**2*d2*dx**2*dy**2 &
-7.8d1*a1*a2**2*dx**2*dy**2-7.8d1*a1**2*a2*dx**2*dy**2 &
-6.0d0*a1*a2**2*d2*dy**2-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2 &
+6.6d1*a1*a2*dy**2+3.3d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4 &
+2.6d1*a1*a2**2*dx**4+2.6d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2 &
+6.0d0*a1**2*a2*d2*dx**2-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2 &
-3.3d1*a1**2*dx**2)*dz
                case (2)
                  rlYlm_laplacian = &
-3.58962931487909853d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(6.0d0*a1**3*a2**3*d2*dy**4-3.9d1*a1**2*a2**3*dy**4 &
-3.9d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
+5.2d1*a1**2*a2**3*dx**2*dy**2+5.2d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.3d1*a1*a2**3*dy**2+6.6d1*a1**2*a2**2*dy**2+3.3d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.3d1*a1**2*a2**3*dx**4 &
-1.3d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.3d1*a1*a2**3*dx**2 &
+6.6d1*a1**2*a2**2*dx**2+3.3d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-2.7d1*a2**3-8.1d1*a1*a2**2 &
-8.1d1*a1**2*a2-2.7d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=3, m1=3, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-3.87724279885581631d0*E*a1**4*a2**4*sqrt(a2+a1)*(a2+a1)**(-10)* &
(1.0d0*a2*(2.0d0*a1*d2-1.5d1)-1.5d1*a1)*dx*dy*(dy**2-3.0d0*dx**2)* &
(3.0d0*dy**2-1.0d0*dx**2)
                case (-2)
                  rlYlm_laplacian =9.49726646607726303d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
-9.0d1*a1**2*a2**3*dx**2*dy**2-9.0d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+3.0d1*a1**2*a2**3*dx**4 &
+3.0d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian = &
-3.00329935783424201d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx*dy* &
(3.0d1*a1**3*a2**3*d2*dy**4-2.25d2*a1**2*a2**3*dy**4 &
-2.25d2*a1**3*a2**2*dy**4+2.0d1*a1**3*a2**3*d2*dx**2*dy**2 &
-1.5d2*a1**2*a2**3*dx**2*dy**2-1.5d2*a1**3*a2**2*dx**2*dy**2 &
-2.4d1*a1**3*a2**3*d2**2*dy**2+1.32d2*a1**2*a2**3*d2*dy**2 &
+1.32d2*a1**3*a2**2*d2*dy**2+3.12d2*a1*a2**3*dy**2 &
+6.24d2*a1**2*a2**2*dy**2+3.12d2*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+7.5d1*a1**2*a2**3*dx**4 &
+7.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-8.4d1*a1**2*a2**3*d2*dx**2-8.4d1*a1**3*a2**2*d2*dx**2 &
+1.56d2*a1*a2**3*dx**2+3.12d2*a1**2*a2**2*dx**2+1.56d2*a1**3*a2*dx**2 &
+2.4d1*a1**2*a2**3*d2**2+2.4d1*a1**3*a2**2*d2**2-1.44d2*a1*a2**3*d2 &
-2.88d2*a1**2*a2**2*d2-1.44d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2 &
-1.98d2*a1**2*a2-6.6d1*a1**3)
                case (0)
                  rlYlm_laplacian = &
-2.45218365717409381d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*(1.0d1*a1**2*a2**2*d2*dy**2 &
-7.5d1*a1*a2**2*dy**2-7.5d1*a1**2*a2*dy**2+1.0d1*a1**2*a2**2*d2*dx**2 &
-7.5d1*a1*a2**2*dx**2-7.5d1*a1**2*a2*dx**2-4.0d0*a1**2*a2**2*d2**2 &
+1.2d1*a1*a2**2*d2+1.2d1*a1**2*a2*d2+1.17d2*a2**2+2.34d2*a1*a2 &
+1.17d2*a1**2)*dz
                case (1)
                  rlYlm_laplacian = &
-1.50164967891712101d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)* &
(6.0d1*a1**3*a2**3*d2*dx**2*dy**4-4.5d2*a1**2*a2**3*dx**2*dy**4 &
-4.5d2*a1**3*a2**2*dx**2*dy**4-3.0d1*a1**2*a2**3*d2*dy**4 &
-3.0d1*a1**3*a2**2*d2*dy**4+1.95d2*a1*a2**3*dy**4 &
+3.9d2*a1**2*a2**2*dy**4+1.95d2*a1**3*a2*dy**4 &
+4.0d1*a1**3*a2**3*d2*dx**4*dy**2-3.0d2*a1**2*a2**3*dx**4*dy**2 &
-3.0d2*a1**3*a2**2*dx**4*dy**2-4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+3.24d2*a1**2*a2**3*d2*dx**2*dy**2+3.24d2*a1**3*a2**2*d2*dx**2*dy**2 &
+2.34d2*a1*a2**3*dx**2*dy**2+4.68d2*a1**2*a2**2*dx**2*dy**2 &
+2.34d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.44d2*a1*a2**3*d2*dy**2 &
-2.88d2*a1**2*a2**2*d2*dy**2-1.44d2*a1**3*a2*d2*dy**2 &
-6.6d1*a2**3*dy**2-1.98d2*a1*a2**2*dy**2-1.98d2*a1**2*a2*dy**2 &
-6.6d1*a1**3*dy**2-2.0d1*a1**3*a2**3*d2*dx**6+1.5d2*a1**2*a2**3*dx**6 &
+1.5d2*a1**3*a2**2*dx**6+1.6d1*a1**3*a2**3*d2**2*dx**4 &
-7.8d1*a1**2*a2**3*d2*dx**4-7.8d1*a1**3*a2**2*d2*dx**4 &
-2.73d2*a1*a2**3*dx**4-5.46d2*a1**2*a2**2*dx**4-2.73d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.44d2*a1*a2**3*d2*dx**2+2.88d2*a1**2*a2**2*d2*dx**2 &
+1.44d2*a1**3*a2*d2*dx**2+6.6d1*a2**3*dx**2+1.98d2*a1*a2**2*dx**2 &
+1.98d2*a1**2*a2*dx**2+6.6d1*a1**3*dx**2)
                case (2)
                  rlYlm_laplacian = &
-9.49726646607726303d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(6.0d0*a1**3*a2**3*d2*dy**4-4.5d1*a1**2*a2**3*dy**4 &
-4.5d1*a1**3*a2**2*dy**4-8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
+6.0d1*a1**2*a2**3*dx**2*dy**2+6.0d1*a1**3*a2**2*dx**2*dy**2 &
-6.0d0*a1**2*a2**3*d2*dy**2-6.0d0*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (3)
                  rlYlm_laplacian =-1.93862139942790816d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-10)*(3.6d1*a1**4*a2**4*d2*dx**2*dy**4 &
-2.7d2*a1**3*a2**4*dx**2*dy**4-2.7d2*a1**4*a2**3*dx**2*dy**4 &
-1.8d1*a1**3*a2**4*d2*dy**4-1.8d1*a1**4*a2**3*d2*dy**4 &
+1.17d2*a1**2*a2**4*dy**4+2.34d2*a1**3*a2**3*dy**4 &
+1.17d2*a1**4*a2**2*dy**4-2.4d1*a1**4*a2**4*d2*dx**4*dy**2 &
+1.8d2*a1**3*a2**4*dx**4*dy**2+1.8d2*a1**4*a2**3*dx**4*dy**2 &
-3.6d1*a1**3*a2**4*d2*dx**2*dy**2-3.6d1*a1**4*a2**3*d2*dx**2*dy**2 &
+2.34d2*a1**2*a2**4*dx**2*dy**2+4.68d2*a1**3*a2**3*dx**2*dy**2 &
+2.34d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**2*a2**4*d2*dy**2 &
+7.2d1*a1**3*a2**3*d2*dy**2+3.6d1*a1**4*a2**2*d2*dy**2 &
-1.98d2*a1*a2**4*dy**2-5.94d2*a1**2*a2**3*dy**2 &
-5.94d2*a1**3*a2**2*dy**2-1.98d2*a1**4*a2*dy**2 &
+4.0d0*a1**4*a2**4*d2*dx**6-3.0d1*a1**3*a2**4*dx**6 &
-3.0d1*a1**4*a2**3*dx**6-1.8d1*a1**3*a2**4*d2*dx**4 &
-1.8d1*a1**4*a2**3*d2*dx**4+1.17d2*a1**2*a2**4*dx**4 &
+2.34d2*a1**3*a2**3*dx**4+1.17d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**2*a2**4*d2*dx**2+7.2d1*a1**3*a2**3*d2*dx**2 &
+3.6d1*a1**4*a2**2*d2*dx**2-1.98d2*a1*a2**4*dx**2 &
-5.94d2*a1**2*a2**3*dx**2-5.94d2*a1**3*a2**2*dx**2 &
-1.98d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+5.4d1*a2**4+2.16d2*a1*a2**3 &
+3.24d2*a1**2*a2**2+2.16d2*a1**3*a2+5.4d1*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=3, m1=3, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-8.2248740261329704d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(1.2d1*a1**4*a2**4*d2*dx**2*dy**4-1.02d2*a1**3*a2**4*dx**2*dy**4 &
-1.02d2*a1**4*a2**3*dx**2*dy**4-6.0d0*a1**3*a2**4*d2*dy**4 &
-6.0d0*a1**4*a2**3*d2*dy**4+4.5d1*a1**2*a2**4*dy**4 &
+9.0d1*a1**3*a2**3*dy**4+4.5d1*a1**4*a2**2*dy**4 &
-1.6d1*a1**4*a2**4*d2*dx**4*dy**2+1.36d2*a1**3*a2**4*dx**4*dy**2 &
+1.36d2*a1**4*a2**3*dx**4*dy**2-1.2d1*a1**3*a2**4*d2*dx**2*dy**2 &
-1.2d1*a1**4*a2**3*d2*dx**2*dy**2+9.0d1*a1**2*a2**4*dx**2*dy**2 &
+1.8d2*a1**3*a2**3*dx**2*dy**2+9.0d1*a1**4*a2**2*dx**2*dy**2 &
+1.8d1*a1**2*a2**4*d2*dy**2+3.6d1*a1**3*a2**3*d2*dy**2 &
+1.8d1*a1**4*a2**2*d2*dy**2-1.17d2*a1*a2**4*dy**2 &
-3.51d2*a1**2*a2**3*dy**2-3.51d2*a1**3*a2**2*dy**2 &
-1.17d2*a1**4*a2*dy**2+4.0d0*a1**4*a2**4*d2*dx**6 &
-3.4d1*a1**3*a2**4*dx**6-3.4d1*a1**4*a2**3*dx**6 &
-6.0d0*a1**3*a2**4*d2*dx**4-6.0d0*a1**4*a2**3*d2*dx**4 &
+4.5d1*a1**2*a2**4*dx**4+9.0d1*a1**3*a2**3*dx**4 &
+4.5d1*a1**4*a2**2*dx**4+1.8d1*a1**2*a2**4*d2*dx**2 &
+3.6d1*a1**3*a2**3*d2*dx**2+1.8d1*a1**4*a2**2*d2*dx**2 &
-1.17d2*a1*a2**4*dx**2-3.51d2*a1**2*a2**3*dx**2 &
-3.51d2*a1**3*a2**2*dx**2-1.17d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2 &
-3.6d1*a1**2*a2**3*d2-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2 &
+6.6d1*a2**4+2.64d2*a1*a2**3+3.96d2*a1**2*a2**2+2.64d2*a1**3*a2 &
+6.6d1*a1**4)
                case (-3)
                  rlYlm_laplacian = &
-1.16317283965674489d1*E*a1**5*a2**4*sqrt(a2+a1)*(a2+a1)**(-11)* &
(1.0d0*a2*(2.0d0*a1*d2-1.7d1)-1.7d1*a1)*dx*dy*(dy**2-3.0d0*dx**2)* &
(3.0d0*dy**2-1.0d0*dx**2)*dz
                case (-2)
                  rlYlm_laplacian = &
-3.10871017685462917d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(8.4d1*a1**4*a2**4*d2*dx**2*dy**4-7.14d2*a1**3*a2**4*dx**2*dy**4 &
-7.14d2*a1**4*a2**3*dx**2*dy**4-4.2d1*a1**3*a2**4*d2*dy**4 &
-4.2d1*a1**4*a2**3*d2*dy**4+3.15d2*a1**2*a2**4*dy**4 &
+6.3d2*a1**3*a2**3*dy**4+3.15d2*a1**4*a2**2*dy**4 &
+5.6d1*a1**4*a2**4*d2*dx**4*dy**2-4.76d2*a1**3*a2**4*dx**4*dy**2 &
-4.76d2*a1**4*a2**3*dx**4*dy**2-7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+4.92d2*a1**3*a2**4*d2*dx**2*dy**2+4.92d2*a1**4*a2**3*d2*dx**2*dy**2 &
+9.0d2*a1**2*a2**4*dx**2*dy**2+1.8d3*a1**3*a2**3*dx**2*dy**2 &
+9.0d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**3*a2**4*d2**2*dy**2 &
+3.6d1*a1**4*a2**3*d2**2*dy**2-2.16d2*a1**2*a2**4*d2*dy**2 &
-4.32d2*a1**3*a2**3*d2*dy**2-2.16d2*a1**4*a2**2*d2*dy**2 &
-3.51d2*a1*a2**4*dy**2-1.053d3*a1**2*a2**3*dy**2 &
-1.053d3*a1**3*a2**2*dy**2-3.51d2*a1**4*a2*dy**2 &
-2.8d1*a1**4*a2**4*d2*dx**6+2.38d2*a1**3*a2**4*dx**6 &
+2.38d2*a1**4*a2**3*dx**6+2.4d1*a1**4*a2**4*d2**2*dx**4 &
-2.34d2*a1**3*a2**4*d2*dx**4-2.34d2*a1**4*a2**3*d2*dx**4 &
+2.25d2*a1**2*a2**4*dx**4+4.5d2*a1**3*a2**3*dx**4 &
+2.25d2*a1**4*a2**2*dx**4+3.6d1*a1**3*a2**4*d2**2*dx**2 &
+3.6d1*a1**4*a2**3*d2**2*dx**2-2.16d2*a1**2*a2**4*d2*dx**2 &
-4.32d2*a1**3*a2**3*d2*dx**2-2.16d2*a1**4*a2**2*d2*dx**2 &
-3.51d2*a1*a2**4*dx**2-1.053d3*a1**2*a2**3*dx**2 &
-1.053d3*a1**3*a2**2*dx**2-3.51d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.28d2*a1*a2**4*d2+6.84d2*a1**2*a2**3*d2 &
+6.84d2*a1**3*a2**2*d2+2.28d2*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case (-1)
                  rlYlm_laplacian = &
-4.39638009359507944d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx*dy* &
(4.2d1*a1**3*a2**3*d2*dy**4-3.57d2*a1**2*a2**3*dy**4 &
-3.57d2*a1**3*a2**2*dy**4+2.8d1*a1**3*a2**3*d2*dx**2*dy**2 &
-2.38d2*a1**2*a2**3*dx**2*dy**2-2.38d2*a1**3*a2**2*dx**2*dy**2 &
-2.4d1*a1**3*a2**3*d2**2*dy**2+1.08d2*a1**2*a2**3*d2*dy**2 &
+1.08d2*a1**3*a2**2*d2*dy**2+7.2d2*a1*a2**3*dy**2 &
+1.44d3*a1**2*a2**2*dy**2+7.2d2*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.19d2*a1**2*a2**3*dx**4 &
+1.19d2*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-9.2d1*a1**2*a2**3*d2*dx**2-9.2d1*a1**3*a2**2*d2*dx**2 &
+1.8d2*a1*a2**3*dx**2+3.6d2*a1**2*a2**2*dx**2+1.8d2*a1**3*a2*dx**2 &
+2.4d1*a1**2*a2**3*d2**2+2.4d1*a1**3*a2**2*d2**2-1.44d2*a1*a2**3*d2 &
-2.88d2*a1**2*a2**2*d2-1.44d2*a1**3*a2*d2-2.34d2*a2**3-7.02d2*a1*a2**2 &
-7.02d2*a1**2*a2-2.34d2*a1**3)*dz
                case (0)
                  rlYlm_laplacian =6.95128727779234418d-1*E*a1**3*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*dx*(3.0d0*dy**2 &
-1.0d0*dx**2)*(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+5.0d2*a1**2*a2**3*d2*dy**2 &
+5.0d2*a1**3*a2**2*d2*dy**2+1.35d3*a1*a2**3*dy**2 &
+2.7d3*a1**2*a2**2*dy**2+1.35d3*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+5.0d2*a1**2*a2**3*d2*dx**2+5.0d2*a1**3*a2**2*d2*dx**2 &
+1.35d3*a1*a2**3*dx**2+2.7d3*a1**2*a2**2*dx**2+1.35d3*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3+8.0d0*a1**2*a2**3*d2**2 &
+8.0d0*a1**3*a2**2*d2**2-1.044d3*a1*a2**3*d2-2.088d3*a1**2*a2**2*d2 &
-1.044d3*a1**3*a2*d2-2.34d2*a2**3-7.02d2*a1*a2**2-7.02d2*a1**2*a2 &
-2.34d2*a1**3)
                case (1)
                  rlYlm_laplacian = &
-2.19819004679753972d0*E*a1**3*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(8.4d1*a1**3*a2**3*d2*dx**2*dy**4-7.14d2*a1**2*a2**3*dx**2*dy**4 &
-7.14d2*a1**3*a2**2*dx**2*dy**4-4.2d1*a1**2*a2**3*d2*dy**4 &
-4.2d1*a1**3*a2**2*d2*dy**4+3.15d2*a1*a2**3*dy**4 &
+6.3d2*a1**2*a2**2*dy**4+3.15d2*a1**3*a2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**4*dy**2-4.76d2*a1**2*a2**3*dx**4*dy**2 &
-4.76d2*a1**3*a2**2*dx**4*dy**2-4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+3.0d2*a1**2*a2**3*d2*dx**2*dy**2+3.0d2*a1**3*a2**2*d2*dx**2*dy**2 &
+8.1d2*a1*a2**3*dx**2*dy**2+1.62d3*a1**2*a2**2*dx**2*dy**2 &
+8.1d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.44d2*a1*a2**3*d2*dy**2 &
-2.88d2*a1**2*a2**2*d2*dy**2-1.44d2*a1**3*a2*d2*dy**2 &
-2.34d2*a2**3*dy**2-7.02d2*a1*a2**2*dy**2-7.02d2*a1**2*a2*dy**2 &
-2.34d2*a1**3*dy**2-2.8d1*a1**3*a2**3*d2*dx**6 &
+2.38d2*a1**2*a2**3*dx**6+2.38d2*a1**3*a2**2*dx**6 &
+1.6d1*a1**3*a2**3*d2**2*dx**4-5.8d1*a1**2*a2**3*d2*dx**4 &
-5.8d1*a1**3*a2**2*d2*dx**4-5.85d2*a1*a2**3*dx**4 &
-1.17d3*a1**2*a2**2*dx**4-5.85d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.44d2*a1*a2**3*d2*dx**2+2.88d2*a1**2*a2**2*d2*dx**2 &
+1.44d2*a1**3*a2*d2*dx**2+2.34d2*a2**3*dx**2+7.02d2*a1*a2**2*dx**2 &
+7.02d2*a1**2*a2*dx**2+2.34d2*a1**3*dx**2)*dz
                case (2)
                  rlYlm_laplacian =3.10871017685462917d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(4.2d1*a1**4*a2**4*d2*dy**6 &
-3.57d2*a1**3*a2**4*dy**6-3.57d2*a1**4*a2**3*dy**6 &
-1.4d1*a1**4*a2**4*d2*dx**2*dy**4+1.19d2*a1**3*a2**4*dx**2*dy**4 &
+1.19d2*a1**4*a2**3*dx**2*dy**4-3.6d1*a1**4*a2**4*d2**2*dy**4 &
+2.46d2*a1**3*a2**4*d2*dy**4+2.46d2*a1**4*a2**3*d2*dy**4 &
+4.5d2*a1**2*a2**4*dy**4+9.0d2*a1**3*a2**3*dy**4 &
+4.5d2*a1**4*a2**2*dy**4-4.2d1*a1**4*a2**4*d2*dx**4*dy**2 &
+3.57d2*a1**3*a2**4*dx**4*dy**2+3.57d2*a1**4*a2**3*dx**4*dy**2 &
+4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2-4.68d2*a1**3*a2**4*d2*dx**2*dy**2 &
-4.68d2*a1**4*a2**3*d2*dx**2*dy**2+4.5d2*a1**2*a2**4*dx**2*dy**2 &
+9.0d2*a1**3*a2**3*dx**2*dy**2+4.5d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**3*a2**4*d2**2*dy**2+3.6d1*a1**4*a2**3*d2**2*dy**2 &
-2.16d2*a1**2*a2**4*d2*dy**2-4.32d2*a1**3*a2**3*d2*dy**2 &
-2.16d2*a1**4*a2**2*d2*dy**2-3.51d2*a1*a2**4*dy**2 &
-1.053d3*a1**2*a2**3*dy**2-1.053d3*a1**3*a2**2*dy**2 &
-3.51d2*a1**4*a2*dy**2+1.4d1*a1**4*a2**4*d2*dx**6 &
-1.19d2*a1**3*a2**4*dx**6-1.19d2*a1**4*a2**3*dx**6 &
-1.2d1*a1**4*a2**4*d2**2*dx**4+5.4d1*a1**3*a2**4*d2*dx**4 &
+5.4d1*a1**4*a2**3*d2*dx**4+3.6d2*a1**2*a2**4*dx**4 &
+7.2d2*a1**3*a2**3*dx**4+3.6d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**3*a2**4*d2**2*dx**2+3.6d1*a1**4*a2**3*d2**2*dx**2 &
-2.16d2*a1**2*a2**4*d2*dx**2-4.32d2*a1**3*a2**3*d2*dx**2 &
-2.16d2*a1**4*a2**2*d2*dx**2-3.51d2*a1*a2**4*dx**2 &
-1.053d3*a1**2*a2**3*dx**2-1.053d3*a1**3*a2**2*dx**2 &
-3.51d2*a1**4*a2*dx**2-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.28d2*a1*a2**4*d2+6.84d2*a1**2*a2**3*d2 &
+6.84d2*a1**3*a2**2*d2+2.28d2*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case (3)
                  rlYlm_laplacian = &
-5.81586419828372447d0*E*a1**2*a2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(3.6d1*a1**4*a2**4*d2*dx**2*dy**4-3.06d2*a1**3*a2**4*dx**2*dy**4 &
-3.06d2*a1**4*a2**3*dx**2*dy**4-1.8d1*a1**3*a2**4*d2*dy**4 &
-1.8d1*a1**4*a2**3*d2*dy**4+1.35d2*a1**2*a2**4*dy**4 &
+2.7d2*a1**3*a2**3*dy**4+1.35d2*a1**4*a2**2*dy**4 &
-2.4d1*a1**4*a2**4*d2*dx**4*dy**2+2.04d2*a1**3*a2**4*dx**4*dy**2 &
+2.04d2*a1**4*a2**3*dx**4*dy**2-3.6d1*a1**3*a2**4*d2*dx**2*dy**2 &
-3.6d1*a1**4*a2**3*d2*dx**2*dy**2+2.7d2*a1**2*a2**4*dx**2*dy**2 &
+5.4d2*a1**3*a2**3*dx**2*dy**2+2.7d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**2*a2**4*d2*dy**2+7.2d1*a1**3*a2**3*d2*dy**2 &
+3.6d1*a1**4*a2**2*d2*dy**2-2.34d2*a1*a2**4*dy**2 &
-7.02d2*a1**2*a2**3*dy**2-7.02d2*a1**3*a2**2*dy**2 &
-2.34d2*a1**4*a2*dy**2+4.0d0*a1**4*a2**4*d2*dx**6 &
-3.4d1*a1**3*a2**4*dx**6-3.4d1*a1**4*a2**3*dx**6 &
-1.8d1*a1**3*a2**4*d2*dx**4-1.8d1*a1**4*a2**3*d2*dx**4 &
+1.35d2*a1**2*a2**4*dx**4+2.7d2*a1**3*a2**3*dx**4 &
+1.35d2*a1**4*a2**2*dx**4+3.6d1*a1**2*a2**4*d2*dx**2 &
+7.2d1*a1**3*a2**3*d2*dx**2+3.6d1*a1**4*a2**2*d2*dx**2 &
-2.34d2*a1*a2**4*dx**2-7.02d2*a1**2*a2**3*dx**2 &
-7.02d2*a1**3*a2**2*dx**2-2.34d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2 &
-3.6d1*a1**2*a2**3*d2-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2 &
+6.6d1*a2**4+2.64d2*a1*a2**3+3.96d2*a1**2*a2**2+2.64d2*a1**3*a2 &
+6.6d1*a1**4)*dz
                case (4)
                  rlYlm_laplacian =4.1124370130664852d0*E*a1**2*a2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(6.0d0*a1**4*a2**4*d2*dy**6 &
-5.1d1*a1**3*a2**4*dy**6-5.1d1*a1**4*a2**3*dy**6 &
-3.8d1*a1**4*a2**4*d2*dx**2*dy**4+3.23d2*a1**3*a2**4*dx**2*dy**4 &
+3.23d2*a1**4*a2**3*dx**2*dy**4+1.2d1*a1**3*a2**4*d2*dy**4 &
+1.2d1*a1**4*a2**3*d2*dy**4-9.0d1*a1**2*a2**4*dy**4 &
-1.8d2*a1**3*a2**3*dy**4-9.0d1*a1**4*a2**2*dy**4 &
+1.8d1*a1**4*a2**4*d2*dx**4*dy**2-1.53d2*a1**3*a2**4*dx**4*dy**2 &
-1.53d2*a1**4*a2**3*dx**4*dy**2+2.4d1*a1**3*a2**4*d2*dx**2*dy**2 &
+2.4d1*a1**4*a2**3*d2*dx**2*dy**2-1.8d2*a1**2*a2**4*dx**2*dy**2 &
-3.6d2*a1**3*a2**3*dx**2*dy**2-1.8d2*a1**4*a2**2*dx**2*dy**2 &
-3.6d1*a1**2*a2**4*d2*dy**2-7.2d1*a1**3*a2**3*d2*dy**2 &
-3.6d1*a1**4*a2**2*d2*dy**2+2.34d2*a1*a2**4*dy**2 &
+7.02d2*a1**2*a2**3*dy**2+7.02d2*a1**3*a2**2*dy**2 &
+2.34d2*a1**4*a2*dy**2-2.0d0*a1**4*a2**4*d2*dx**6 &
+1.7d1*a1**3*a2**4*dx**6+1.7d1*a1**4*a2**3*dx**6 &
+1.2d1*a1**3*a2**4*d2*dx**4+1.2d1*a1**4*a2**3*d2*dx**4 &
-9.0d1*a1**2*a2**4*dx**4-1.8d2*a1**3*a2**3*dx**4 &
-9.0d1*a1**4*a2**2*dx**4-3.6d1*a1**2*a2**4*d2*dx**2 &
-7.2d1*a1**3*a2**3*d2*dx**2-3.6d1*a1**4*a2**2*d2*dx**2 &
+2.34d2*a1*a2**4*dx**2+7.02d2*a1**2*a2**3*dx**2 &
+7.02d2*a1**3*a2**2*dx**2+2.34d2*a1**4*a2*dx**2+2.4d1*a1*a2**4*d2 &
+7.2d1*a1**2*a2**3*d2+7.2d1*a1**3*a2**2*d2+2.4d1*a1**4*a2*d2 &
-1.32d2*a2**4-5.28d2*a1*a2**3-7.92d2*a1**2*a2**2-5.28d2*a1**3*a2 &
-1.32d2*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case default
          print *,'Error: rlYlm_overlap not implemented for l1=',l1 &
,'m1=',m1,'l2=',l2,'m2=',m2
          stop
      end select
    case (4)
      ! selection on m1: l1=4
      select case (m1)
        case (-4)
          ! selection on l2: l1=4, m1=-4
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=-4, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-7.86448379536438834d0*E*a1*a2**5*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=-4, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-6.81084275443661897d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(4.0d0*a1**2*a2**2*d2*dy**4-2.6d1*a1*a2**2*dy**4-2.6d1*a1**2*a2*dy**4 &
-4.0d0*a1**2*a2**2*d2*dx**2*dy**2+2.6d1*a1*a2**2*dx**2*dy**2 &
+2.6d1*a1**2*a2*dx**2*dy**2-6.0d0*a1*a2**2*d2*dy**2 &
-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2+6.6d1*a1*a2*dy**2 &
+3.3d1*a1**2*dy**2+2.0d0*a1*a2**2*d2*dx**2+2.0d0*a1**2*a2*d2*dx**2 &
-1.1d1*a2**2*dx**2-2.2d1*a1*a2*dx**2-1.1d1*a1**2*dx**2)
                case (0)
                  rlYlm_laplacian = &
-1.36216855088732379d1*E*a1**2*a2**5*sqrt(a2+a1)*(a2+a1)**(-9)* &
(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)*dz
                case (1)
                  rlYlm_laplacian = &
-6.81084275443661897d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(4.0d0*a1**2*a2**2*d2*dx**2*dy**2-2.6d1*a1*a2**2*dx**2*dy**2 &
-2.6d1*a1**2*a2*dx**2*dy**2-2.0d0*a1*a2**2*d2*dy**2 &
-2.0d0*a1**2*a2*d2*dy**2+1.1d1*a2**2*dy**2+2.2d1*a1*a2*dy**2 &
+1.1d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4+2.6d1*a1*a2**2*dx**4 &
+2.6d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2+6.0d0*a1**2*a2*d2*dx**2 &
-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2-3.3d1*a1**2*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=-4, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-7.6147536914910937d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*(dy+dx)* &
(dy-1.0d0*dx)*(8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
-6.0d1*a1**2*a2**3*dx**2*dy**2-6.0d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+2.6d1*a1*a2**3*dy**2+5.2d1*a1**2*a2**2*dy**2+2.6d1*a1**3*a2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+2.6d1*a1*a2**3*dx**2+5.2d1*a1**2*a2**2*dx**2+2.6d1*a1**3*a2*dx**2 &
+6.0d0*a1*a2**3*d2+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3 &
-9.9d1*a1*a2**2-9.9d1*a1**2*a2-3.3d1*a1**3)
                case (-1)
                  rlYlm_laplacian = &
-1.52295073829821874d1*E*a1**2*a2**4*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(4.0d0*a1**2*a2**2*d2*dy**4-3.0d1*a1*a2**2*dy**4-3.0d1*a1**2*a2*dy**4 &
-4.0d0*a1**2*a2**2*d2*dx**2*dy**2+3.0d1*a1*a2**2*dx**2*dy**2 &
+3.0d1*a1**2*a2*dx**2*dy**2-6.0d0*a1*a2**2*d2*dy**2 &
-6.0d0*a1**2*a2*d2*dy**2+3.9d1*a2**2*dy**2+7.8d1*a1*a2*dy**2 &
+3.9d1*a1**2*dy**2+2.0d0*a1*a2**2*d2*dx**2+2.0d0*a1**2*a2*d2*dx**2 &
-1.3d1*a2**2*dx**2-2.6d1*a1*a2*dx**2-1.3d1*a1**2*dx**2)*dz
                case (0)
                  rlYlm_laplacian =8.79276018719015889d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(6.0d0*a1**2*a2**2*d2*dy**2-4.5d1*a1*a2**2*dy**2-4.5d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-4.5d1*a1*a2**2*dx**2-4.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.2d1*a1*a2**2*d2+2.2d1*a1**2*a2*d2 &
+5.2d1*a2**2+1.04d2*a1*a2+5.2d1*a1**2)
                case (1)
                  rlYlm_laplacian = &
-1.52295073829821874d1*E*a1**2*a2**4*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(4.0d0*a1**2*a2**2*d2*dx**2*dy**2-3.0d1*a1*a2**2*dx**2*dy**2 &
-3.0d1*a1**2*a2*dx**2*dy**2-2.0d0*a1*a2**2*d2*dy**2 &
-2.0d0*a1**2*a2*d2*dy**2+1.3d1*a2**2*dy**2+2.6d1*a1*a2*dy**2 &
+1.3d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4+3.0d1*a1*a2**2*dx**4 &
+3.0d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2+6.0d0*a1**2*a2*d2*dx**2 &
-3.9d1*a2**2*dx**2-7.8d1*a1*a2*dx**2-3.9d1*a1**2*dx**2)*dz
                case (2)
                  rlYlm_laplacian =1.52295073829821874d1*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.5d1*a1**2*a2**3*dy**4-1.5d1*a1**3*a2**2*dy**4 &
-4.0d0*a1**3*a2**3*d2*dx**2*dy**2+3.0d1*a1**2*a2**3*dx**2*dy**2 &
+3.0d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**2*a2**3*d2*dy**2 &
-4.0d0*a1**3*a2**2*d2*dy**2+2.6d1*a1*a2**3*dy**2 &
+5.2d1*a1**2*a2**2*dy**2+2.6d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4-4.0d0*a1**2*a2**3*d2*dx**2 &
-4.0d0*a1**3*a2**2*d2*dx**2+2.6d1*a1*a2**3*dx**2 &
+5.2d1*a1**2*a2**2*dx**2+2.6d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=-4, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =8.2248740261329704d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(4.0d0*a1**4*a2**4*d2*dy**6 &
-3.4d1*a1**3*a2**4*dy**6-3.4d1*a1**4*a2**3*dy**6 &
-1.6d1*a1**4*a2**4*d2*dx**2*dy**4+1.36d2*a1**3*a2**4*dx**2*dy**4 &
+1.36d2*a1**4*a2**3*dx**2*dy**4-6.0d0*a1**3*a2**4*d2*dy**4 &
-6.0d0*a1**4*a2**3*d2*dy**4+4.5d1*a1**2*a2**4*dy**4 &
+9.0d1*a1**3*a2**3*dy**4+4.5d1*a1**4*a2**2*dy**4 &
+1.2d1*a1**4*a2**4*d2*dx**4*dy**2-1.02d2*a1**3*a2**4*dx**4*dy**2 &
-1.02d2*a1**4*a2**3*dx**4*dy**2-1.2d1*a1**3*a2**4*d2*dx**2*dy**2 &
-1.2d1*a1**4*a2**3*d2*dx**2*dy**2+9.0d1*a1**2*a2**4*dx**2*dy**2 &
+1.8d2*a1**3*a2**3*dx**2*dy**2+9.0d1*a1**4*a2**2*dx**2*dy**2 &
+1.8d1*a1**2*a2**4*d2*dy**2+3.6d1*a1**3*a2**3*d2*dy**2 &
+1.8d1*a1**4*a2**2*d2*dy**2-1.17d2*a1*a2**4*dy**2 &
-3.51d2*a1**2*a2**3*dy**2-3.51d2*a1**3*a2**2*dy**2 &
-1.17d2*a1**4*a2*dy**2-6.0d0*a1**3*a2**4*d2*dx**4 &
-6.0d0*a1**4*a2**3*d2*dx**4+4.5d1*a1**2*a2**4*dx**4 &
+9.0d1*a1**3*a2**3*dx**4+4.5d1*a1**4*a2**2*dx**4 &
+1.8d1*a1**2*a2**4*d2*dx**2+3.6d1*a1**3*a2**3*d2*dx**2 &
+1.8d1*a1**4*a2**2*d2*dx**2-1.17d2*a1*a2**4*dx**2 &
-3.51d2*a1**2*a2**3*dx**2-3.51d2*a1**3*a2**2*dx**2 &
-1.17d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+6.6d1*a2**4+2.64d2*a1*a2**3 &
+3.96d2*a1**2*a2**2+2.64d2*a1**3*a2+6.6d1*a1**4)
                case (-2)
                  rlYlm_laplacian = &
-2.01467445626964921d1*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*(dy+dx &
)*(dy-1.0d0*dx)*(8.0d0*a1**3*a2**3*d2*dx**2*dy**2 &
-6.8d1*a1**2*a2**3*dx**2*dy**2-6.8d1*a1**3*a2**2*dx**2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dy**2-4.0d0*a1**3*a2**2*d2*dy**2 &
+3.0d1*a1*a2**3*dy**2+6.0d1*a1**2*a2**2*dy**2+3.0d1*a1**3*a2*dy**2 &
-4.0d0*a1**2*a2**3*d2*dx**2-4.0d0*a1**3*a2**2*d2*dx**2 &
+3.0d1*a1*a2**3*dx**2+6.0d1*a1**2*a2**2*dx**2+3.0d1*a1**3*a2*dx**2 &
+6.0d0*a1*a2**3*d2+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.9d1*a2**3 &
-1.17d2*a1*a2**2-1.17d2*a1**2*a2-3.9d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =6.37096002557338818d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(2.0d1*a1**3*a2**3*d2*dy**6 &
-1.7d2*a1**2*a2**3*dy**6-1.7d2*a1**3*a2**2*dy**6 &
-1.6d1*a1**3*a2**3*d2**2*dy**4+9.0d1*a1**2*a2**3*d2*dy**4 &
+9.0d1*a1**3*a2**2*d2*dy**4+3.45d2*a1*a2**3*dy**4 &
+6.9d2*a1**2*a2**2*dy**4+3.45d2*a1**3*a2*dy**4 &
-2.0d1*a1**3*a2**3*d2*dx**4*dy**2+1.7d2*a1**2*a2**3*dx**4*dy**2 &
+1.7d2*a1**3*a2**2*dx**4*dy**2+1.6d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-1.4d2*a1**2*a2**3*d2*dx**2*dy**2-1.4d2*a1**3*a2**2*d2*dx**2*dy**2 &
+3.0d1*a1*a2**3*dx**2*dy**2+6.0d1*a1**2*a2**2*dx**2*dy**2 &
+3.0d1*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.62d2*a1*a2**3*d2*dy**2 &
-3.24d2*a1**2*a2**2*d2*dy**2-1.62d2*a1**3*a2*d2*dy**2 &
-1.17d2*a2**3*dy**2-3.51d2*a1*a2**2*dy**2-3.51d2*a1**2*a2*dy**2 &
-1.17d2*a1**3*dy**2+1.0d1*a1**2*a2**3*d2*dx**4 &
+1.0d1*a1**3*a2**2*d2*dx**4-7.5d1*a1*a2**3*dx**4 &
-1.5d2*a1**2*a2**2*dx**4-7.5d1*a1**3*a2*dx**4 &
-8.0d0*a1**2*a2**3*d2**2*dx**2-8.0d0*a1**3*a2**2*d2**2*dx**2 &
+5.4d1*a1*a2**3*d2*dx**2+1.08d2*a1**2*a2**2*d2*dx**2 &
+5.4d1*a1**3*a2*d2*dx**2+3.9d1*a2**3*dx**2+1.17d2*a1*a2**2*dx**2 &
+1.17d2*a1**2*a2*dx**2+3.9d1*a1**3*dx**2)
                case (0)
                  rlYlm_laplacian =1.04037341562157789d1*E*a1**3*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(1.0d1*a1**2*a2**2*d2*dy**2-8.5d1*a1*a2**2*dy**2-8.5d1*a1**2*a2*dy**2 &
+1.0d1*a1**2*a2**2*d2*dx**2-8.5d1*a1*a2**2*dx**2-8.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.0d1*a1*a2**2*d2+1.0d1*a1**2*a2*d2 &
+1.8d2*a2**2+3.6d2*a1*a2+1.8d2*a1**2)*dz
                case (1)
                  rlYlm_laplacian =6.37096002557338818d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(2.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
-1.7d2*a1**2*a2**3*dx**2*dy**4-1.7d2*a1**3*a2**2*dx**2*dy**4 &
-1.0d1*a1**2*a2**3*d2*dy**4-1.0d1*a1**3*a2**2*d2*dy**4 &
+7.5d1*a1*a2**3*dy**4+1.5d2*a1**2*a2**2*dy**4+7.5d1*a1**3*a2*dy**4 &
-1.6d1*a1**3*a2**3*d2**2*dx**2*dy**2+1.4d2*a1**2*a2**3*d2*dx**2*dy**2 &
+1.4d2*a1**3*a2**2*d2*dx**2*dy**2-3.0d1*a1*a2**3*dx**2*dy**2 &
-6.0d1*a1**2*a2**2*dx**2*dy**2-3.0d1*a1**3*a2*dx**2*dy**2 &
+8.0d0*a1**2*a2**3*d2**2*dy**2+8.0d0*a1**3*a2**2*d2**2*dy**2 &
-5.4d1*a1*a2**3*d2*dy**2-1.08d2*a1**2*a2**2*d2*dy**2 &
-5.4d1*a1**3*a2*d2*dy**2-3.9d1*a2**3*dy**2-1.17d2*a1*a2**2*dy**2 &
-1.17d2*a1**2*a2*dy**2-3.9d1*a1**3*dy**2-2.0d1*a1**3*a2**3*d2*dx**6 &
+1.7d2*a1**2*a2**3*dx**6+1.7d2*a1**3*a2**2*dx**6 &
+1.6d1*a1**3*a2**3*d2**2*dx**4-9.0d1*a1**2*a2**3*d2*dx**4 &
-9.0d1*a1**3*a2**2*d2*dx**4-3.45d2*a1*a2**3*dx**4 &
-6.9d2*a1**2*a2**2*dx**4-3.45d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.62d2*a1*a2**3*d2*dx**2+3.24d2*a1**2*a2**2*d2*dx**2 &
+1.62d2*a1**3*a2*d2*dx**2+1.17d2*a2**3*dx**2+3.51d2*a1*a2**2*dx**2 &
+3.51d2*a1**2*a2*dx**2+1.17d2*a1**3*dx**2)
                case (2)
                  rlYlm_laplacian =4.02934891253929843d1*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.7d1*a1**2*a2**3*dy**4-1.7d1*a1**3*a2**2*dy**4 &
-4.0d0*a1**3*a2**3*d2*dx**2*dy**2+3.4d1*a1**2*a2**3*dx**2*dy**2 &
+3.4d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**2*a2**3*d2*dy**2 &
-4.0d0*a1**3*a2**2*d2*dy**2+3.0d1*a1*a2**3*dy**2 &
+6.0d1*a1**2*a2**2*dy**2+3.0d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.7d1*a1**2*a2**3*dx**4 &
-1.7d1*a1**3*a2**2*dx**4-4.0d0*a1**2*a2**3*d2*dx**2 &
-4.0d0*a1**3*a2**2*d2*dx**2+3.0d1*a1*a2**3*dx**2 &
+6.0d1*a1**2*a2**2*dx**2+3.0d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.9d1*a2**3-1.17d2*a1*a2**2 &
-1.17d2*a1**2*a2-3.9d1*a1**3)*dz
                case (3)
                  rlYlm_laplacian =8.2248740261329704d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(1.2d1*a1**4*a2**4*d2*dx**2*dy**4 &
-1.02d2*a1**3*a2**4*dx**2*dy**4-1.02d2*a1**4*a2**3*dx**2*dy**4 &
-6.0d0*a1**3*a2**4*d2*dy**4-6.0d0*a1**4*a2**3*d2*dy**4 &
+4.5d1*a1**2*a2**4*dy**4+9.0d1*a1**3*a2**3*dy**4 &
+4.5d1*a1**4*a2**2*dy**4-1.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
+1.36d2*a1**3*a2**4*dx**4*dy**2+1.36d2*a1**4*a2**3*dx**4*dy**2 &
-1.2d1*a1**3*a2**4*d2*dx**2*dy**2-1.2d1*a1**4*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**4*dx**2*dy**2+1.8d2*a1**3*a2**3*dx**2*dy**2 &
+9.0d1*a1**4*a2**2*dx**2*dy**2+1.8d1*a1**2*a2**4*d2*dy**2 &
+3.6d1*a1**3*a2**3*d2*dy**2+1.8d1*a1**4*a2**2*d2*dy**2 &
-1.17d2*a1*a2**4*dy**2-3.51d2*a1**2*a2**3*dy**2 &
-3.51d2*a1**3*a2**2*dy**2-1.17d2*a1**4*a2*dy**2 &
+4.0d0*a1**4*a2**4*d2*dx**6-3.4d1*a1**3*a2**4*dx**6 &
-3.4d1*a1**4*a2**3*dx**6-6.0d0*a1**3*a2**4*d2*dx**4 &
-6.0d0*a1**4*a2**3*d2*dx**4+4.5d1*a1**2*a2**4*dx**4 &
+9.0d1*a1**3*a2**3*dx**4+4.5d1*a1**4*a2**2*dx**4 &
+1.8d1*a1**2*a2**4*d2*dx**2+3.6d1*a1**3*a2**3*d2*dx**2 &
+1.8d1*a1**4*a2**2*d2*dx**2-1.17d2*a1*a2**4*dx**2 &
-3.51d2*a1**2*a2**3*dx**2-3.51d2*a1**3*a2**2*dx**2 &
-1.17d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+6.6d1*a2**4+2.64d2*a1*a2**3 &
+3.96d2*a1**2*a2**2+2.64d2*a1**3*a2+6.6d1*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=-4, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =1.74475925948511734d1*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(8.0d0*a1**5*a2**5*d2*dx**2*dy**6 &
-7.6d1*a1**4*a2**5*dx**2*dy**6-7.6d1*a1**5*a2**4*dx**2*dy**6 &
-4.0d0*a1**4*a2**5*d2*dy**6-4.0d0*a1**5*a2**4*d2*dy**6 &
+3.4d1*a1**3*a2**5*dy**6+6.8d1*a1**4*a2**4*dy**6 &
+3.4d1*a1**5*a2**3*dy**6-1.6d1*a1**5*a2**5*d2*dx**4*dy**4 &
+1.52d2*a1**4*a2**5*dx**4*dy**4+1.52d2*a1**5*a2**4*dx**4*dy**4 &
-1.2d1*a1**4*a2**5*d2*dx**2*dy**4-1.2d1*a1**5*a2**4*d2*dx**2*dy**4 &
+1.02d2*a1**3*a2**5*dx**2*dy**4+2.04d2*a1**4*a2**4*dx**2*dy**4 &
+1.02d2*a1**5*a2**3*dx**2*dy**4+1.8d1*a1**3*a2**5*d2*dy**4 &
+3.6d1*a1**4*a2**4*d2*dy**4+1.8d1*a1**5*a2**3*d2*dy**4 &
-1.35d2*a1**2*a2**5*dy**4-4.05d2*a1**3*a2**4*dy**4 &
-4.05d2*a1**4*a2**3*dy**4-1.35d2*a1**5*a2**2*dy**4 &
+8.0d0*a1**5*a2**5*d2*dx**6*dy**2-7.6d1*a1**4*a2**5*dx**6*dy**2 &
-7.6d1*a1**5*a2**4*dx**6*dy**2-1.2d1*a1**4*a2**5*d2*dx**4*dy**2 &
-1.2d1*a1**5*a2**4*d2*dx**4*dy**2+1.02d2*a1**3*a2**5*dx**4*dy**2 &
+2.04d2*a1**4*a2**4*dx**4*dy**2+1.02d2*a1**5*a2**3*dx**4*dy**2 &
+3.6d1*a1**3*a2**5*d2*dx**2*dy**2+7.2d1*a1**4*a2**4*d2*dx**2*dy**2 &
+3.6d1*a1**5*a2**3*d2*dx**2*dy**2-2.7d2*a1**2*a2**5*dx**2*dy**2 &
-8.1d2*a1**3*a2**4*dx**2*dy**2-8.1d2*a1**4*a2**3*dx**2*dy**2 &
-2.7d2*a1**5*a2**2*dx**2*dy**2-2.4d1*a1**2*a2**5*d2*dy**2 &
-7.2d1*a1**3*a2**4*d2*dy**2-7.2d1*a1**4*a2**3*d2*dy**2 &
-2.4d1*a1**5*a2**2*d2*dy**2+1.56d2*a1*a2**5*dy**2 &
+6.24d2*a1**2*a2**4*dy**2+9.36d2*a1**3*a2**3*dy**2 &
+6.24d2*a1**4*a2**2*dy**2+1.56d2*a1**5*a2*dy**2 &
-4.0d0*a1**4*a2**5*d2*dx**6-4.0d0*a1**5*a2**4*d2*dx**6 &
+3.4d1*a1**3*a2**5*dx**6+6.8d1*a1**4*a2**4*dx**6 &
+3.4d1*a1**5*a2**3*dx**6+1.8d1*a1**3*a2**5*d2*dx**4 &
+3.6d1*a1**4*a2**4*d2*dx**4+1.8d1*a1**5*a2**3*d2*dx**4 &
-1.35d2*a1**2*a2**5*dx**4-4.05d2*a1**3*a2**4*dx**4 &
-4.05d2*a1**4*a2**3*dx**4-1.35d2*a1**5*a2**2*dx**4 &
-2.4d1*a1**2*a2**5*d2*dx**2-7.2d1*a1**3*a2**4*d2*dx**2 &
-7.2d1*a1**4*a2**3*d2*dx**2-2.4d1*a1**5*a2**2*d2*dx**2 &
+1.56d2*a1*a2**5*dx**2+6.24d2*a1**2*a2**4*dx**2 &
+9.36d2*a1**3*a2**3*dx**2+6.24d2*a1**4*a2**2*dx**2 &
+1.56d2*a1**5*a2*dx**2+6.0d0*a1*a2**5*d2+2.4d1*a1**2*a2**4*d2 &
+3.6d1*a1**3*a2**3*d2+2.4d1*a1**4*a2**2*d2+6.0d0*a1**5*a2*d2 &
-3.3d1*a2**5-1.65d2*a1*a2**4-3.3d2*a1**2*a2**3-3.3d2*a1**3*a2**2 &
-1.65d2*a1**4*a2-3.3d1*a1**5)
                case (-3)
                  rlYlm_laplacian =2.46746220783989112d1*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dx*(4.0d0*a1**4*a2**4*d2*dy**6 &
-3.8d1*a1**3*a2**4*dy**6-3.8d1*a1**4*a2**3*dy**6 &
-1.6d1*a1**4*a2**4*d2*dx**2*dy**4+1.52d2*a1**3*a2**4*dx**2*dy**4 &
+1.52d2*a1**4*a2**3*dx**2*dy**4-6.0d0*a1**3*a2**4*d2*dy**4 &
-6.0d0*a1**4*a2**3*d2*dy**4+5.1d1*a1**2*a2**4*dy**4 &
+1.02d2*a1**3*a2**3*dy**4+5.1d1*a1**4*a2**2*dy**4 &
+1.2d1*a1**4*a2**4*d2*dx**4*dy**2-1.14d2*a1**3*a2**4*dx**4*dy**2 &
-1.14d2*a1**4*a2**3*dx**4*dy**2-1.2d1*a1**3*a2**4*d2*dx**2*dy**2 &
-1.2d1*a1**4*a2**3*d2*dx**2*dy**2+1.02d2*a1**2*a2**4*dx**2*dy**2 &
+2.04d2*a1**3*a2**3*dx**2*dy**2+1.02d2*a1**4*a2**2*dx**2*dy**2 &
+1.8d1*a1**2*a2**4*d2*dy**2+3.6d1*a1**3*a2**3*d2*dy**2 &
+1.8d1*a1**4*a2**2*d2*dy**2-1.35d2*a1*a2**4*dy**2 &
-4.05d2*a1**2*a2**3*dy**2-4.05d2*a1**3*a2**2*dy**2 &
-1.35d2*a1**4*a2*dy**2-6.0d0*a1**3*a2**4*d2*dx**4 &
-6.0d0*a1**4*a2**3*d2*dx**4+5.1d1*a1**2*a2**4*dx**4 &
+1.02d2*a1**3*a2**3*dx**4+5.1d1*a1**4*a2**2*dx**4 &
+1.8d1*a1**2*a2**4*d2*dx**2+3.6d1*a1**3*a2**3*d2*dx**2 &
+1.8d1*a1**4*a2**2*d2*dx**2-1.35d2*a1*a2**4*dx**2 &
-4.05d2*a1**2*a2**3*dx**2-4.05d2*a1**3*a2**2*dx**2 &
-1.35d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+7.8d1*a2**4+3.12d2*a1*a2**3 &
+4.68d2*a1**2*a2**2+3.12d2*a1**3*a2+7.8d1*a1**4)*dz
                case (-2)
                  rlYlm_laplacian =6.59457014039261917d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(dy+dx)*(dy-1.0d0*dx)* &
(5.6d1*a1**4*a2**4*d2*dx**2*dy**4-5.32d2*a1**3*a2**4*dx**2*dy**4 &
-5.32d2*a1**4*a2**3*dx**2*dy**4-2.8d1*a1**3*a2**4*d2*dy**4 &
-2.8d1*a1**4*a2**3*d2*dy**4+2.38d2*a1**2*a2**4*dy**4 &
+4.76d2*a1**3*a2**3*dy**4+2.38d2*a1**4*a2**2*dy**4 &
+5.6d1*a1**4*a2**4*d2*dx**4*dy**2-5.32d2*a1**3*a2**4*dx**4*dy**2 &
-5.32d2*a1**4*a2**3*dx**4*dy**2-4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+3.68d2*a1**3*a2**4*d2*dx**2*dy**2+3.68d2*a1**4*a2**3*d2*dx**2*dy**2 &
+7.48d2*a1**2*a2**4*dx**2*dy**2+1.496d3*a1**3*a2**3*dx**2*dy**2 &
+7.48d2*a1**4*a2**2*dx**2*dy**2+2.4d1*a1**3*a2**4*d2**2*dy**2 &
+2.4d1*a1**4*a2**3*d2**2*dy**2-1.5d2*a1**2*a2**4*d2*dy**2 &
-3.0d2*a1**3*a2**3*d2*dy**2-1.5d2*a1**4*a2**2*d2*dy**2 &
-4.05d2*a1*a2**4*dy**2-1.215d3*a1**2*a2**3*dy**2 &
-1.215d3*a1**3*a2**2*dy**2-4.05d2*a1**4*a2*dy**2 &
-2.8d1*a1**3*a2**4*d2*dx**4-2.8d1*a1**4*a2**3*d2*dx**4 &
+2.38d2*a1**2*a2**4*dx**4+4.76d2*a1**3*a2**3*dx**4 &
+2.38d2*a1**4*a2**2*dx**4+2.4d1*a1**3*a2**4*d2**2*dx**2 &
+2.4d1*a1**4*a2**3*d2**2*dx**2-1.5d2*a1**2*a2**4*d2*dx**2 &
-3.0d2*a1**3*a2**3*d2*dx**2-1.5d2*a1**4*a2**2*d2*dx**2 &
-4.05d2*a1*a2**4*dx**2-1.215d3*a1**2*a2**3*dx**2 &
-1.215d3*a1**3*a2**2*dx**2-4.05d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.58d2*a1*a2**4*d2+7.74d2*a1**2*a2**3*d2 &
+7.74d2*a1**3*a2**2*d2+2.58d2*a1**4*a2*d2+7.8d1*a2**4+3.12d2*a1*a2**3 &
+4.68d2*a1**2*a2**2+3.12d2*a1**3*a2+7.8d1*a1**4)
                case (-1)
                  rlYlm_laplacian =9.3261305305638875d0*E*a1**3*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-12)*dx*(2.8d1*a1**3*a2**3*d2*dy**6 &
-2.66d2*a1**2*a2**3*dy**6-2.66d2*a1**3*a2**2*dy**6 &
-1.6d1*a1**3*a2**3*d2**2*dy**4+6.2d1*a1**2*a2**3*d2*dy**4 &
+6.2d1*a1**3*a2**2*d2*dy**4+7.65d2*a1*a2**3*dy**4 &
+1.53d3*a1**2*a2**2*dy**4+7.65d2*a1**3*a2*dy**4 &
-2.8d1*a1**3*a2**3*d2*dx**4*dy**2+2.66d2*a1**2*a2**3*dx**4*dy**2 &
+2.66d2*a1**3*a2**2*dx**4*dy**2+1.6d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-1.32d2*a1**2*a2**3*d2*dx**2*dy**2-1.32d2*a1**3*a2**2*d2*dx**2*dy**2 &
-1.7d2*a1*a2**3*dx**2*dy**2-3.4d2*a1**2*a2**2*dx**2*dy**2 &
-1.7d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.5d2*a1*a2**3*d2*dy**2 &
-3.0d2*a1**2*a2**2*d2*dy**2-1.5d2*a1**3*a2*d2*dy**2-4.05d2*a2**3*dy**2 &
-1.215d3*a1*a2**2*dy**2-1.215d3*a1**2*a2*dy**2-4.05d2*a1**3*dy**2 &
+1.4d1*a1**2*a2**3*d2*dx**4+1.4d1*a1**3*a2**2*d2*dx**4 &
-1.19d2*a1*a2**3*dx**4-2.38d2*a1**2*a2**2*dx**4-1.19d2*a1**3*a2*dx**4 &
-8.0d0*a1**2*a2**3*d2**2*dx**2-8.0d0*a1**3*a2**2*d2**2*dx**2 &
+5.0d1*a1*a2**3*d2*dx**2+1.0d2*a1**2*a2**2*d2*dx**2 &
+5.0d1*a1**3*a2*d2*dx**2+1.35d2*a2**3*dx**2+4.05d2*a1*a2**2*dx**2 &
+4.05d2*a1**2*a2*dx**2+1.35d2*a1**3*dx**2)*dz
                case (0)
                  rlYlm_laplacian = &
-2.94918142326164563d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(7.0d1*a1**3*a2**3*d2*dy**4 &
-6.65d2*a1**2*a2**3*dy**4-6.65d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.33d3*a1**2*a2**3*dx**2*dy**2 &
-1.33d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+5.2d2*a1**2*a2**3*d2*dy**2+5.2d2*a1**3*a2**2*d2*dy**2 &
+2.04d3*a1*a2**3*dy**2+4.08d3*a1**2*a2**2*dy**2+2.04d3*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-6.65d2*a1**2*a2**3*dx**4 &
-6.65d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+5.2d2*a1**2*a2**3*d2*dx**2+5.2d2*a1**3*a2**2*d2*dx**2 &
+2.04d3*a1*a2**3*dx**2+4.08d3*a1**2*a2**2*dx**2+2.04d3*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3+4.0d1*a1**2*a2**3*d2**2 &
+4.0d1*a1**3*a2**2*d2**2-1.56d3*a1*a2**3*d2-3.12d3*a1**2*a2**2*d2 &
-1.56d3*a1**3*a2*d2-5.4d2*a2**3-1.62d3*a1*a2**2-1.62d3*a1**2*a2 &
-5.4d2*a1**3)
                case (1)
                  rlYlm_laplacian =9.3261305305638875d0*E*a1**3*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-12)*dy*(2.8d1*a1**3*a2**3*d2*dx**2*dy**4 &
-2.66d2*a1**2*a2**3*dx**2*dy**4-2.66d2*a1**3*a2**2*dx**2*dy**4 &
-1.4d1*a1**2*a2**3*d2*dy**4-1.4d1*a1**3*a2**2*d2*dy**4 &
+1.19d2*a1*a2**3*dy**4+2.38d2*a1**2*a2**2*dy**4+1.19d2*a1**3*a2*dy**4 &
-1.6d1*a1**3*a2**3*d2**2*dx**2*dy**2+1.32d2*a1**2*a2**3*d2*dx**2*dy**2 &
+1.32d2*a1**3*a2**2*d2*dx**2*dy**2+1.7d2*a1*a2**3*dx**2*dy**2 &
+3.4d2*a1**2*a2**2*dx**2*dy**2+1.7d2*a1**3*a2*dx**2*dy**2 &
+8.0d0*a1**2*a2**3*d2**2*dy**2+8.0d0*a1**3*a2**2*d2**2*dy**2 &
-5.0d1*a1*a2**3*d2*dy**2-1.0d2*a1**2*a2**2*d2*dy**2 &
-5.0d1*a1**3*a2*d2*dy**2-1.35d2*a2**3*dy**2-4.05d2*a1*a2**2*dy**2 &
-4.05d2*a1**2*a2*dy**2-1.35d2*a1**3*dy**2-2.8d1*a1**3*a2**3*d2*dx**6 &
+2.66d2*a1**2*a2**3*dx**6+2.66d2*a1**3*a2**2*dx**6 &
+1.6d1*a1**3*a2**3*d2**2*dx**4-6.2d1*a1**2*a2**3*d2*dx**4 &
-6.2d1*a1**3*a2**2*d2*dx**4-7.65d2*a1*a2**3*dx**4 &
-1.53d3*a1**2*a2**2*dx**4-7.65d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.5d2*a1*a2**3*d2*dx**2+3.0d2*a1**2*a2**2*d2*dx**2 &
+1.5d2*a1**3*a2*d2*dx**2+4.05d2*a2**3*dx**2+1.215d3*a1*a2**2*dx**2 &
+1.215d3*a1**2*a2*dx**2+4.05d2*a1**3*dx**2)*dz
                case (2)
                  rlYlm_laplacian = &
-1.31891402807852383d1*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.33d2*a1**3*a2**4*dy**6 &
-1.33d2*a1**4*a2**3*dy**6-1.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.33d2*a1**3*a2**4*dx**2*dy**4+1.33d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+7.8d1*a1**3*a2**4*d2*dy**4 &
+7.8d1*a1**4*a2**3*d2*dy**4+3.06d2*a1**2*a2**4*dy**4 &
+6.12d2*a1**3*a2**3*dy**4+3.06d2*a1**4*a2**2*dy**4 &
-1.4d1*a1**4*a2**4*d2*dx**4*dy**2+1.33d2*a1**3*a2**4*dx**4*dy**2 &
+1.33d2*a1**4*a2**3*dx**4*dy**2+2.4d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-2.68d2*a1**3*a2**4*d2*dx**2*dy**2-2.68d2*a1**4*a2**3*d2*dx**2*dy**2 &
+3.4d2*a1**2*a2**4*dx**2*dy**2+6.8d2*a1**3*a2**3*dx**2*dy**2 &
+3.4d2*a1**4*a2**2*dx**2*dy**2+2.4d1*a1**3*a2**4*d2**2*dy**2 &
+2.4d1*a1**4*a2**3*d2**2*dy**2-1.5d2*a1**2*a2**4*d2*dy**2 &
-3.0d2*a1**3*a2**3*d2*dy**2-1.5d2*a1**4*a2**2*d2*dy**2 &
-4.05d2*a1*a2**4*dy**2-1.215d3*a1**2*a2**3*dy**2 &
-1.215d3*a1**3*a2**2*dy**2-4.05d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.33d2*a1**3*a2**4*dx**6 &
-1.33d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+7.8d1*a1**3*a2**4*d2*dx**4+7.8d1*a1**4*a2**3*d2*dx**4 &
+3.06d2*a1**2*a2**4*dx**4+6.12d2*a1**3*a2**3*dx**4 &
+3.06d2*a1**4*a2**2*dx**4+2.4d1*a1**3*a2**4*d2**2*dx**2 &
+2.4d1*a1**4*a2**3*d2**2*dx**2-1.5d2*a1**2*a2**4*d2*dx**2 &
-3.0d2*a1**3*a2**3*d2*dx**2-1.5d2*a1**4*a2**2*d2*dx**2 &
-4.05d2*a1*a2**4*dx**2-1.215d3*a1**2*a2**3*dx**2 &
-1.215d3*a1**3*a2**2*dx**2-4.05d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.58d2*a1*a2**4*d2+7.74d2*a1**2*a2**3*d2 &
+7.74d2*a1**3*a2**2*d2+2.58d2*a1**4*a2*d2+7.8d1*a2**4+3.12d2*a1*a2**3 &
+4.68d2*a1**2*a2**2+3.12d2*a1**3*a2+7.8d1*a1**4)
                case (3)
                  rlYlm_laplacian =2.46746220783989112d1*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dy*(1.2d1*a1**4*a2**4*d2*dx**2*dy**4 &
-1.14d2*a1**3*a2**4*dx**2*dy**4-1.14d2*a1**4*a2**3*dx**2*dy**4 &
-6.0d0*a1**3*a2**4*d2*dy**4-6.0d0*a1**4*a2**3*d2*dy**4 &
+5.1d1*a1**2*a2**4*dy**4+1.02d2*a1**3*a2**3*dy**4 &
+5.1d1*a1**4*a2**2*dy**4-1.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
+1.52d2*a1**3*a2**4*dx**4*dy**2+1.52d2*a1**4*a2**3*dx**4*dy**2 &
-1.2d1*a1**3*a2**4*d2*dx**2*dy**2-1.2d1*a1**4*a2**3*d2*dx**2*dy**2 &
+1.02d2*a1**2*a2**4*dx**2*dy**2+2.04d2*a1**3*a2**3*dx**2*dy**2 &
+1.02d2*a1**4*a2**2*dx**2*dy**2+1.8d1*a1**2*a2**4*d2*dy**2 &
+3.6d1*a1**3*a2**3*d2*dy**2+1.8d1*a1**4*a2**2*d2*dy**2 &
-1.35d2*a1*a2**4*dy**2-4.05d2*a1**2*a2**3*dy**2 &
-4.05d2*a1**3*a2**2*dy**2-1.35d2*a1**4*a2*dy**2 &
+4.0d0*a1**4*a2**4*d2*dx**6-3.8d1*a1**3*a2**4*dx**6 &
-3.8d1*a1**4*a2**3*dx**6-6.0d0*a1**3*a2**4*d2*dx**4 &
-6.0d0*a1**4*a2**3*d2*dx**4+5.1d1*a1**2*a2**4*dx**4 &
+1.02d2*a1**3*a2**3*dx**4+5.1d1*a1**4*a2**2*dx**4 &
+1.8d1*a1**2*a2**4*d2*dx**2+3.6d1*a1**3*a2**3*d2*dx**2 &
+1.8d1*a1**4*a2**2*d2*dx**2-1.35d2*a1*a2**4*dx**2 &
-4.05d2*a1**2*a2**3*dx**2-4.05d2*a1**3*a2**2*dx**2 &
-1.35d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+7.8d1*a2**4+3.12d2*a1*a2**3 &
+4.68d2*a1**2*a2**2+3.12d2*a1**3*a2+7.8d1*a1**4)*dz
                case (4)
                  rlYlm_laplacian = &
-1.74475925948511734d1*E*a1**5*a2**5*sqrt(a2+a1)*(a2+a1)**(-12)* &
(1.0d0*a2*(2.0d0*a1*d2-1.9d1)-1.9d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(dy**2-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (-3)
          ! selection on l2: l1=4, m1=-3
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=-3, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-5.56102982223387534d0*E*a1*a2**5*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dy*(dy**2-3.0d0*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=-3, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-4.81599309725739696d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)* &
(4.0d0*a1**2*a2**2*d2*dy**4-2.6d1*a1*a2**2*dy**4-2.6d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+7.8d1*a1*a2**2*dx**2*dy**2 &
+7.8d1*a1**2*a2*dx**2*dy**2-6.0d0*a1*a2**2*d2*dy**2 &
-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2+6.6d1*a1*a2*dy**2 &
+3.3d1*a1**2*dy**2+6.0d0*a1*a2**2*d2*dx**2+6.0d0*a1**2*a2*d2*dx**2 &
-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2-3.3d1*a1**2*dx**2)*dz
                case (0)
                  rlYlm_laplacian =4.81599309725739696d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(dy**2-3.0d0*dx**2)* &
(4.0d0*a1**2*a2**2*d2*dy**2-2.6d1*a1*a2**2*dy**2-2.6d1*a1**2*a2*dy**2 &
+4.0d0*a1**2*a2**2*d2*dx**2-2.6d1*a1*a2**2*dx**2-2.6d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.8d1*a1*a2**2*d2+2.8d1*a1**2*a2*d2 &
-1.1d1*a2**2-2.2d1*a1*a2-1.1d1*a1**2)
                case (1)
                  rlYlm_laplacian = &
-9.63198619451479393d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(2.0d0*a1**2*a2**2*d2*dy**2-1.3d1*a1*a2**2*dy**2-1.3d1*a1**2*a2*dy**2 &
-6.0d0*a1**2*a2**2*d2*dx**2+3.9d1*a1*a2**2*dx**2+3.9d1*a1**2*a2*dx**2 &
+6.0d0*a1*a2**2*d2+6.0d0*a1**2*a2*d2-3.3d1*a2**2-6.6d1*a1*a2 &
-3.3d1*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=-3, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-1.07688879446372956d1*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(4.0d0*a1**3*a2**3*d2*dy**4-3.0d1*a1**2*a2**3*dy**4 &
-3.0d1*a1**3*a2**2*dy**4-1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**3*dx**2*dy**2+9.0d1*a1**3*a2**2*dx**2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dy**2+6.0d0*a1**3*a2**2*d2*dy**2 &
-3.9d1*a1*a2**3*dy**2-7.8d1*a1**2*a2**2*dy**2-3.9d1*a1**3*a2*dy**2 &
+6.0d0*a1**2*a2**3*d2*dx**2+6.0d0*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
-6.0d0*a1*a2**3*d2-1.2d1*a1**2*a2**2*d2-6.0d0*a1**3*a2*d2+3.3d1*a2**3 &
+9.9d1*a1*a2**2+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =5.3844439723186478d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*(8.0d0*a1**3*a2**3*d2*dy**6 &
-6.0d1*a1**2*a2**3*dy**6-6.0d1*a1**3*a2**2*dy**6 &
-1.6d1*a1**3*a2**3*d2*dx**2*dy**4+1.2d2*a1**2*a2**3*dx**2*dy**4 &
+1.2d2*a1**3*a2**2*dx**2*dy**4-8.0d0*a1**3*a2**3*d2**2*dy**4 &
+5.2d1*a1**2*a2**3*d2*dy**4+5.2d1*a1**3*a2**2*d2*dy**4 &
+5.2d1*a1*a2**3*dy**4+1.04d2*a1**2*a2**2*dy**4+5.2d1*a1**3*a2*dy**4 &
-2.4d1*a1**3*a2**3*d2*dx**4*dy**2+1.8d2*a1**2*a2**3*dx**4*dy**2 &
+1.8d2*a1**3*a2**2*dx**4*dy**2+2.4d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-1.92d2*a1**2*a2**3*d2*dx**2*dy**2-1.92d2*a1**3*a2**2*d2*dx**2*dy**2 &
+7.8d1*a1*a2**3*dx**2*dy**2+1.56d2*a1**2*a2**2*dx**2*dy**2 &
+7.8d1*a1**3*a2*dx**2*dy**2+1.2d1*a1**2*a2**3*d2**2*dy**2 &
+1.2d1*a1**3*a2**2*d2**2*dy**2-8.4d1*a1*a2**3*d2*dy**2 &
-1.68d2*a1**2*a2**2*d2*dy**2-8.4d1*a1**3*a2*d2*dy**2+3.3d1*a2**3*dy**2 &
+9.9d1*a1*a2**2*dy**2+9.9d1*a1**2*a2*dy**2+3.3d1*a1**3*dy**2 &
+1.2d1*a1**2*a2**3*d2*dx**4+1.2d1*a1**3*a2**2*d2*dx**4 &
-7.8d1*a1*a2**3*dx**4-1.56d2*a1**2*a2**2*dx**4-7.8d1*a1**3*a2*dx**4 &
-1.2d1*a1**2*a2**3*d2**2*dx**2-1.2d1*a1**3*a2**2*d2**2*dx**2 &
+8.4d1*a1*a2**3*d2*dx**2+1.68d2*a1**2*a2**2*d2*dx**2 &
+8.4d1*a1**3*a2*d2*dx**2-3.3d1*a2**3*dx**2-9.9d1*a1*a2**2*dx**2 &
-9.9d1*a1**2*a2*dx**2-3.3d1*a1**3*dx**2)
                case (0)
                  rlYlm_laplacian =6.21742035370925833d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(dy**2-3.0d0*dx**2)* &
(6.0d0*a1**2*a2**2*d2*dy**2-4.5d1*a1*a2**2*dy**2-4.5d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-4.5d1*a1*a2**2*dx**2-4.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.8d1*a1*a2**2*d2+2.8d1*a1**2*a2*d2 &
+1.3d1*a2**2+2.6d1*a1*a2+1.3d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =1.07688879446372956d1*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(4.0d0*a1**3*a2**3*d2*dy**4 &
-3.0d1*a1**2*a2**3*dy**4-3.0d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2*dx**2*dy**2+6.0d1*a1**2*a2**3*dx**2*dy**2 &
+6.0d1*a1**3*a2**2*dx**2*dy**2-4.0d0*a1**3*a2**3*d2**2*dy**2 &
+4.4d1*a1**2*a2**3*d2*dy**2+4.4d1*a1**3*a2**2*d2*dy**2 &
-9.1d1*a1*a2**3*dy**2-1.82d2*a1**2*a2**2*dy**2-9.1d1*a1**3*a2*dy**2 &
-1.2d1*a1**3*a2**3*d2*dx**4+9.0d1*a1**2*a2**3*dx**4 &
+9.0d1*a1**3*a2**2*dx**4+1.2d1*a1**3*a2**3*d2**2*dx**2 &
-8.4d1*a1**2*a2**3*d2*dx**2-8.4d1*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
-1.2d1*a1**2*a2**3*d2**2-1.2d1*a1**3*a2**2*d2**2+8.4d1*a1*a2**3*d2 &
+1.68d2*a1**2*a2**2*d2+8.4d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)
                case (2)
                  rlYlm_laplacian =1.07688879446372956d1*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.5d1*a1**2*a2**3*dy**4-1.5d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2*dx**2*dy**2+6.0d1*a1**2*a2**3*dx**2*dy**2 &
+6.0d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+6.0d0*a1**3*a2**3*d2*dx**4-4.5d1*a1**2*a2**3*dx**4 &
-4.5d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=-3, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =5.81586419828372447d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(4.0d0*a1**4*a2**4*d2*dy**6 &
-3.4d1*a1**3*a2**4*dy**6-3.4d1*a1**4*a2**3*dy**6 &
-2.4d1*a1**4*a2**4*d2*dx**2*dy**4+2.04d2*a1**3*a2**4*dx**2*dy**4 &
+2.04d2*a1**4*a2**3*dx**2*dy**4-1.8d1*a1**3*a2**4*d2*dy**4 &
-1.8d1*a1**4*a2**3*d2*dy**4+1.35d2*a1**2*a2**4*dy**4 &
+2.7d2*a1**3*a2**3*dy**4+1.35d2*a1**4*a2**2*dy**4 &
+3.6d1*a1**4*a2**4*d2*dx**4*dy**2-3.06d2*a1**3*a2**4*dx**4*dy**2 &
-3.06d2*a1**4*a2**3*dx**4*dy**2-3.6d1*a1**3*a2**4*d2*dx**2*dy**2 &
-3.6d1*a1**4*a2**3*d2*dx**2*dy**2+2.7d2*a1**2*a2**4*dx**2*dy**2 &
+5.4d2*a1**3*a2**3*dx**2*dy**2+2.7d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**2*a2**4*d2*dy**2+7.2d1*a1**3*a2**3*d2*dy**2 &
+3.6d1*a1**4*a2**2*d2*dy**2-2.34d2*a1*a2**4*dy**2 &
-7.02d2*a1**2*a2**3*dy**2-7.02d2*a1**3*a2**2*dy**2 &
-2.34d2*a1**4*a2*dy**2-1.8d1*a1**3*a2**4*d2*dx**4 &
-1.8d1*a1**4*a2**3*d2*dx**4+1.35d2*a1**2*a2**4*dx**4 &
+2.7d2*a1**3*a2**3*dx**4+1.35d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**2*a2**4*d2*dx**2+7.2d1*a1**3*a2**3*d2*dx**2 &
+3.6d1*a1**4*a2**2*d2*dx**2-2.34d2*a1*a2**4*dx**2 &
-7.02d2*a1**2*a2**3*dx**2-7.02d2*a1**3*a2**2*dx**2 &
-2.34d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+6.6d1*a2**4+2.64d2*a1*a2**3 &
+3.96d2*a1**2*a2**2+2.64d2*a1**3*a2+6.6d1*a1**4)*dz
                case (-2)
                  rlYlm_laplacian =1.42458996991158946d1*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(8.0d0*a1**4*a2**4*d2*dy**6 &
-6.8d1*a1**3*a2**4*dy**6-6.8d1*a1**4*a2**3*dy**6 &
-1.6d1*a1**4*a2**4*d2*dx**2*dy**4+1.36d2*a1**3*a2**4*dx**2*dy**4 &
+1.36d2*a1**4*a2**3*dx**2*dy**4-8.0d0*a1**4*a2**4*d2**2*dy**4 &
+8.4d1*a1**3*a2**4*d2*dy**4+8.4d1*a1**4*a2**3*d2*dy**4 &
-1.2d2*a1**2*a2**4*dy**4-2.4d2*a1**3*a2**3*dy**4 &
-1.2d2*a1**4*a2**2*dy**4-2.4d1*a1**4*a2**4*d2*dx**4*dy**2 &
+2.04d2*a1**3*a2**4*dx**4*dy**2+2.04d2*a1**4*a2**3*dx**4*dy**2 &
+2.4d1*a1**4*a2**4*d2**2*dx**2*dy**2-1.92d2*a1**3*a2**4*d2*dx**2*dy**2 &
-1.92d2*a1**4*a2**3*d2*dx**2*dy**2-9.0d1*a1**2*a2**4*dx**2*dy**2 &
-1.8d2*a1**3*a2**3*dx**2*dy**2-9.0d1*a1**4*a2**2*dx**2*dy**2 &
-1.2d1*a1**3*a2**4*d2**2*dy**2-1.2d1*a1**4*a2**3*d2**2*dy**2 &
+8.4d1*a1**2*a2**4*d2*dy**2+1.68d2*a1**3*a2**3*d2*dy**2 &
+8.4d1*a1**4*a2**2*d2*dy**2+3.9d1*a1*a2**4*dy**2 &
+1.17d2*a1**2*a2**3*dy**2+1.17d2*a1**3*a2**2*dy**2 &
+3.9d1*a1**4*a2*dy**2+1.2d1*a1**3*a2**4*d2*dx**4 &
+1.2d1*a1**4*a2**3*d2*dx**4-9.0d1*a1**2*a2**4*dx**4 &
-1.8d2*a1**3*a2**3*dx**4-9.0d1*a1**4*a2**2*dx**4 &
-1.2d1*a1**3*a2**4*d2**2*dx**2-1.2d1*a1**4*a2**3*d2**2*dx**2 &
+8.4d1*a1**2*a2**4*d2*dx**2+1.68d2*a1**3*a2**3*d2*dx**2 &
+8.4d1*a1**4*a2**2*d2*dx**2+3.9d1*a1*a2**4*dx**2 &
+1.17d2*a1**2*a2**3*dx**2+1.17d2*a1**3*a2**2*dx**2 &
+3.9d1*a1**4*a2*dx**2+1.2d1*a1**2*a2**4*d2**2+2.4d1*a1**3*a2**3*d2**2 &
+1.2d1*a1**4*a2**2*d2**2-8.4d1*a1*a2**4*d2-2.52d2*a1**2*a2**3*d2 &
-2.52d2*a1**3*a2**2*d2-8.4d1*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case (-1)
                  rlYlm_laplacian =4.50494903675136302d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*(2.0d1*a1**3*a2**3*d2*dy**6 &
-1.7d2*a1**2*a2**3*dy**6-1.7d2*a1**3*a2**2*dy**6 &
-4.0d1*a1**3*a2**3*d2*dx**2*dy**4+3.4d2*a1**2*a2**3*dx**2*dy**4 &
+3.4d2*a1**3*a2**2*dx**2*dy**4-1.6d1*a1**3*a2**3*d2**2*dy**4 &
+1.1d2*a1**2*a2**3*d2*dy**4+1.1d2*a1**3*a2**2*d2*dy**4 &
+1.95d2*a1*a2**3*dy**4+3.9d2*a1**2*a2**2*dy**4+1.95d2*a1**3*a2*dy**4 &
-6.0d1*a1**3*a2**3*d2*dx**4*dy**2+5.1d2*a1**2*a2**3*dx**4*dy**2 &
+5.1d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-4.2d2*a1**2*a2**3*d2*dx**2*dy**2-4.2d2*a1**3*a2**2*d2*dx**2*dy**2 &
+9.0d1*a1*a2**3*dx**2*dy**2+1.8d2*a1**2*a2**2*dx**2*dy**2 &
+9.0d1*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.92d2*a1*a2**3*d2*dy**2 &
-3.84d2*a1**2*a2**2*d2*dy**2-1.92d2*a1**3*a2*d2*dy**2 &
+7.8d1*a2**3*dy**2+2.34d2*a1*a2**2*dy**2+2.34d2*a1**2*a2*dy**2 &
+7.8d1*a1**3*dy**2+3.0d1*a1**2*a2**3*d2*dx**4 &
+3.0d1*a1**3*a2**2*d2*dx**4-2.25d2*a1*a2**3*dx**4 &
-4.5d2*a1**2*a2**2*dx**4-2.25d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.92d2*a1*a2**3*d2*dx**2+3.84d2*a1**2*a2**2*d2*dx**2 &
+1.92d2*a1**3*a2*d2*dx**2-7.8d1*a2**3*dx**2-2.34d2*a1*a2**2*dx**2 &
-2.34d2*a1**2*a2*dx**2-7.8d1*a1**3*dx**2)*dz
                case (0)
                  rlYlm_laplacian = &
-3.67827548576114071d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(dy**2-3.0d0*dx**2)*(2.0d1*a1**3*a2**3*d2*dy**4 &
-1.7d2*a1**2*a2**3*dy**4-1.7d2*a1**3*a2**2*dy**4 &
+4.0d1*a1**3*a2**3*d2*dx**2*dy**2-3.4d2*a1**2*a2**3*dx**2*dy**2 &
-3.4d2*a1**3*a2**2*dx**2*dy**2-2.8d1*a1**3*a2**3*d2**2*dy**2 &
+2.2d2*a1**2*a2**3*d2*dy**2+2.2d2*a1**3*a2**2*d2*dy**2 &
+1.35d2*a1*a2**3*dy**2+2.7d2*a1**2*a2**2*dy**2+1.35d2*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.7d2*a1**2*a2**3*dx**4 &
-1.7d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+2.2d2*a1**2*a2**3*d2*dx**2+2.2d2*a1**3*a2**2*d2*dx**2 &
+1.35d2*a1*a2**3*dx**2+2.7d2*a1**2*a2**2*dx**2+1.35d2*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-4.4d1*a1**2*a2**3*d2**2 &
-4.4d1*a1**3*a2**2*d2**2-1.98d2*a1*a2**3*d2-3.96d2*a1**2*a2**2*d2 &
-1.98d2*a1**3*a2*d2+1.17d2*a2**3+3.51d2*a1*a2**2+3.51d2*a1**2*a2 &
+1.17d2*a1**3)
                case (1)
                  rlYlm_laplacian =9.00989807350272604d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(1.0d1*a1**3*a2**3*d2*dy**4 &
-8.5d1*a1**2*a2**3*dy**4-8.5d1*a1**3*a2**2*dy**4 &
-2.0d1*a1**3*a2**3*d2*dx**2*dy**2+1.7d2*a1**2*a2**3*dx**2*dy**2 &
+1.7d2*a1**3*a2**2*dx**2*dy**2-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+1.0d2*a1**2*a2**3*d2*dy**2+1.0d2*a1**3*a2**2*d2*dy**2 &
-2.4d2*a1*a2**3*dy**2-4.8d2*a1**2*a2**2*dy**2-2.4d2*a1**3*a2*dy**2 &
-3.0d1*a1**3*a2**3*d2*dx**4+2.55d2*a1**2*a2**3*dx**4 &
+2.55d2*a1**3*a2**2*dx**4+2.4d1*a1**3*a2**3*d2**2*dx**2 &
-1.8d2*a1**2*a2**3*d2*dx**2-1.8d2*a1**3*a2**2*d2*dx**2 &
-1.8d2*a1*a2**3*dx**2-3.6d2*a1**2*a2**2*dx**2-1.8d2*a1**3*a2*dx**2 &
-2.4d1*a1**2*a2**3*d2**2-2.4d1*a1**3*a2**2*d2**2+1.92d2*a1*a2**3*d2 &
+3.84d2*a1**2*a2**2*d2+1.92d2*a1**3*a2*d2-7.8d1*a2**3-2.34d2*a1*a2**2 &
-2.34d2*a1**2*a2-7.8d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian = &
-1.42458996991158946d1*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(4.0d0*a1**4*a2**4*d2*dy**6-3.4d1*a1**3*a2**4*dy**6 &
-3.4d1*a1**4*a2**3*dy**6-1.2d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.02d2*a1**3*a2**4*dx**2*dy**4+1.02d2*a1**4*a2**3*dx**2*dy**4 &
-4.0d0*a1**4*a2**4*d2**2*dy**4+2.4d1*a1**3*a2**4*d2*dy**4 &
+2.4d1*a1**4*a2**3*d2*dy**4+7.5d1*a1**2*a2**4*dy**4 &
+1.5d2*a1**3*a2**3*dy**4+7.5d1*a1**4*a2**2*dy**4 &
-4.0d0*a1**4*a2**4*d2*dx**4*dy**2+3.4d1*a1**3*a2**4*dx**4*dy**2 &
+3.4d1*a1**4*a2**3*dx**4*dy**2+1.6d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-1.68d2*a1**3*a2**4*d2*dx**2*dy**2-1.68d2*a1**4*a2**3*d2*dx**2*dy**2 &
+2.4d2*a1**2*a2**4*dx**2*dy**2+4.8d2*a1**3*a2**3*dx**2*dy**2 &
+2.4d2*a1**4*a2**2*dx**2*dy**2+1.2d1*a1**3*a2**4*d2**2*dy**2 &
+1.2d1*a1**4*a2**3*d2**2*dy**2-8.4d1*a1**2*a2**4*d2*dy**2 &
-1.68d2*a1**3*a2**3*d2*dy**2-8.4d1*a1**4*a2**2*d2*dy**2 &
-3.9d1*a1*a2**4*dy**2-1.17d2*a1**2*a2**3*dy**2 &
-1.17d2*a1**3*a2**2*dy**2-3.9d1*a1**4*a2*dy**2 &
+1.2d1*a1**4*a2**4*d2*dx**6-1.02d2*a1**3*a2**4*dx**6 &
-1.02d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+9.6d1*a1**3*a2**4*d2*dx**4+9.6d1*a1**4*a2**3*d2*dx**4 &
+4.5d1*a1**2*a2**4*dx**4+9.0d1*a1**3*a2**3*dx**4 &
+4.5d1*a1**4*a2**2*dx**4+1.2d1*a1**3*a2**4*d2**2*dx**2 &
+1.2d1*a1**4*a2**3*d2**2*dx**2-8.4d1*a1**2*a2**4*d2*dx**2 &
-1.68d2*a1**3*a2**3*d2*dx**2-8.4d1*a1**4*a2**2*d2*dx**2 &
-3.9d1*a1*a2**4*dx**2-1.17d2*a1**2*a2**3*dx**2 &
-1.17d2*a1**3*a2**2*dx**2-3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2 &
-2.4d1*a1**3*a2**3*d2**2-1.2d1*a1**4*a2**2*d2**2+8.4d1*a1*a2**4*d2 &
+2.52d2*a1**2*a2**3*d2+2.52d2*a1**3*a2**2*d2+8.4d1*a1**4*a2*d2 &
-3.3d1*a2**4-1.32d2*a1*a2**3-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2 &
-3.3d1*a1**4)
                case (3)
                  rlYlm_laplacian =1.16317283965674489d1*E*a1**4*a2**5*sqrt &
(a2+a1)*(a2+a1)**(-11)*(1.0d0*a2*(2.0d0*a1*d2-1.7d1)-1.7d1*a1)*dx*dy* &
(dy**2-3.0d0*dx**2)*(3.0d0*dy**2-1.0d0*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=-3, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =2.46746220783989112d1*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dx*(4.0d0*a1**4*a2**4*d2*dy**6 &
-3.8d1*a1**3*a2**4*dy**6-3.8d1*a1**4*a2**3*dy**6 &
-1.6d1*a1**4*a2**4*d2*dx**2*dy**4+1.52d2*a1**3*a2**4*dx**2*dy**4 &
+1.52d2*a1**4*a2**3*dx**2*dy**4-6.0d0*a1**3*a2**4*d2*dy**4 &
-6.0d0*a1**4*a2**3*d2*dy**4+5.1d1*a1**2*a2**4*dy**4 &
+1.02d2*a1**3*a2**3*dy**4+5.1d1*a1**4*a2**2*dy**4 &
+1.2d1*a1**4*a2**4*d2*dx**4*dy**2-1.14d2*a1**3*a2**4*dx**4*dy**2 &
-1.14d2*a1**4*a2**3*dx**4*dy**2-1.2d1*a1**3*a2**4*d2*dx**2*dy**2 &
-1.2d1*a1**4*a2**3*d2*dx**2*dy**2+1.02d2*a1**2*a2**4*dx**2*dy**2 &
+2.04d2*a1**3*a2**3*dx**2*dy**2+1.02d2*a1**4*a2**2*dx**2*dy**2 &
+1.8d1*a1**2*a2**4*d2*dy**2+3.6d1*a1**3*a2**3*d2*dy**2 &
+1.8d1*a1**4*a2**2*d2*dy**2-1.35d2*a1*a2**4*dy**2 &
-4.05d2*a1**2*a2**3*dy**2-4.05d2*a1**3*a2**2*dy**2 &
-1.35d2*a1**4*a2*dy**2-6.0d0*a1**3*a2**4*d2*dx**4 &
-6.0d0*a1**4*a2**3*d2*dx**4+5.1d1*a1**2*a2**4*dx**4 &
+1.02d2*a1**3*a2**3*dx**4+5.1d1*a1**4*a2**2*dx**4 &
+1.8d1*a1**2*a2**4*d2*dx**2+3.6d1*a1**3*a2**3*d2*dx**2 &
+1.8d1*a1**4*a2**2*d2*dx**2-1.35d2*a1*a2**4*dx**2 &
-4.05d2*a1**2*a2**3*dx**2-4.05d2*a1**3*a2**2*dx**2 &
-1.35d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+7.8d1*a2**4+3.12d2*a1*a2**3 &
+4.68d2*a1**2*a2**2+3.12d2*a1**3*a2+7.8d1*a1**4)*dz
                case (-3)
                  rlYlm_laplacian =-8.7237962974255867d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(8.0d0*a1**5*a2**5*d2*dy**8 &
-7.6d1*a1**4*a2**5*dy**8-7.6d1*a1**5*a2**4*dy**8 &
-4.0d1*a1**5*a2**5*d2*dx**2*dy**6+3.8d2*a1**4*a2**5*dx**2*dy**6 &
+3.8d2*a1**5*a2**4*dx**2*dy**6-8.0d0*a1**5*a2**5*d2**2*dy**6 &
+4.4d1*a1**4*a2**5*d2*dy**6+4.4d1*a1**5*a2**4*d2*dy**6 &
+2.72d2*a1**3*a2**5*dy**6+5.44d2*a1**4*a2**4*dy**6 &
+2.72d2*a1**5*a2**3*dy**6+2.4d1*a1**5*a2**5*d2*dx**4*dy**4 &
-2.28d2*a1**4*a2**5*dx**4*dy**4-2.28d2*a1**5*a2**4*dx**4*dy**4 &
+4.8d1*a1**5*a2**5*d2**2*dx**2*dy**4-5.88d2*a1**4*a2**5*d2*dx**2*dy**4 &
-5.88d2*a1**5*a2**4*d2*dx**2*dy**4+1.122d3*a1**3*a2**5*dx**2*dy**4 &
+2.244d3*a1**4*a2**4*dx**2*dy**4+1.122d3*a1**5*a2**3*dx**2*dy**4 &
+3.6d1*a1**4*a2**5*d2**2*dy**4+3.6d1*a1**5*a2**4*d2**2*dy**4 &
-2.52d2*a1**3*a2**5*d2*dy**4-5.04d2*a1**4*a2**4*d2*dy**4 &
-2.52d2*a1**5*a2**3*d2*dy**4-4.05d2*a1**2*a2**5*dy**4 &
-1.215d3*a1**3*a2**4*dy**4-1.215d3*a1**4*a2**3*dy**4 &
-4.05d2*a1**5*a2**2*dy**4+7.2d1*a1**5*a2**5*d2*dx**6*dy**2 &
-6.84d2*a1**4*a2**5*dx**6*dy**2-6.84d2*a1**5*a2**4*dx**6*dy**2 &
-7.2d1*a1**5*a2**5*d2**2*dx**4*dy**2+6.12d2*a1**4*a2**5*d2*dx**4*dy**2 &
+6.12d2*a1**5*a2**4*d2*dx**4*dy**2+6.12d2*a1**3*a2**5*dx**4*dy**2 &
+1.224d3*a1**4*a2**4*dx**4*dy**2+6.12d2*a1**5*a2**3*dx**4*dy**2 &
+7.2d1*a1**4*a2**5*d2**2*dx**2*dy**2 &
+7.2d1*a1**5*a2**4*d2**2*dx**2*dy**2-5.04d2*a1**3*a2**5*d2*dx**2*dy**2 &
-1.008d3*a1**4*a2**4*d2*dx**2*dy**2-5.04d2*a1**5*a2**3*d2*dx**2*dy**2 &
-8.1d2*a1**2*a2**5*dx**2*dy**2-2.43d3*a1**3*a2**4*dx**2*dy**2 &
-2.43d3*a1**4*a2**3*dx**2*dy**2-8.1d2*a1**5*a2**2*dx**2*dy**2 &
-7.2d1*a1**3*a2**5*d2**2*dy**2-1.44d2*a1**4*a2**4*d2**2*dy**2 &
-7.2d1*a1**5*a2**3*d2**2*dy**2+5.52d2*a1**2*a2**5*d2*dy**2 &
+1.656d3*a1**3*a2**4*d2*dy**2+1.656d3*a1**4*a2**3*d2*dy**2 &
+5.52d2*a1**5*a2**2*d2*dy**2-7.8d1*a1*a2**5*dy**2 &
-3.12d2*a1**2*a2**4*dy**2-4.68d2*a1**3*a2**3*dy**2 &
-3.12d2*a1**4*a2**2*dy**2-7.8d1*a1**5*a2*dy**2 &
-3.6d1*a1**4*a2**5*d2*dx**6-3.6d1*a1**5*a2**4*d2*dx**6 &
+3.06d2*a1**3*a2**5*dx**6+6.12d2*a1**4*a2**4*dx**6 &
+3.06d2*a1**5*a2**3*dx**6+3.6d1*a1**4*a2**5*d2**2*dx**4 &
+3.6d1*a1**5*a2**4*d2**2*dx**4-2.52d2*a1**3*a2**5*d2*dx**4 &
-5.04d2*a1**4*a2**4*d2*dx**4-2.52d2*a1**5*a2**3*d2*dx**4 &
-4.05d2*a1**2*a2**5*dx**4-1.215d3*a1**3*a2**4*dx**4 &
-1.215d3*a1**4*a2**3*dx**4-4.05d2*a1**5*a2**2*dx**4 &
-7.2d1*a1**3*a2**5*d2**2*dx**2-1.44d2*a1**4*a2**4*d2**2*dx**2 &
-7.2d1*a1**5*a2**3*d2**2*dx**2+5.52d2*a1**2*a2**5*d2*dx**2 &
+1.656d3*a1**3*a2**4*d2*dx**2+1.656d3*a1**4*a2**3*d2*dx**2 &
+5.52d2*a1**5*a2**2*d2*dx**2-7.8d1*a1*a2**5*dx**2 &
-3.12d2*a1**2*a2**4*dx**2-4.68d2*a1**3*a2**3*dx**2 &
-3.12d2*a1**4*a2**2*dx**2-7.8d1*a1**5*a2*dx**2+2.4d1*a1**2*a2**5*d2**2 &
+7.2d1*a1**3*a2**4*d2**2+7.2d1*a1**4*a2**3*d2**2 &
+2.4d1*a1**5*a2**2*d2**2-1.68d2*a1*a2**5*d2-6.72d2*a1**2*a2**4*d2 &
-1.008d3*a1**3*a2**3*d2-6.72d2*a1**4*a2**2*d2-1.68d2*a1**5*a2*d2 &
+6.6d1*a2**5+3.3d2*a1*a2**4+6.6d2*a1**2*a2**3+6.6d2*a1**3*a2**2 &
+3.3d2*a1**4*a2+6.6d1*a1**5)
                case (-2)
                  rlYlm_laplacian =9.3261305305638875d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dx*(2.8d1*a1**4*a2**4*d2*dy**6 &
-2.66d2*a1**3*a2**4*dy**6-2.66d2*a1**4*a2**3*dy**6 &
-5.6d1*a1**4*a2**4*d2*dx**2*dy**4+5.32d2*a1**3*a2**4*dx**2*dy**4 &
+5.32d2*a1**4*a2**3*dx**2*dy**4-2.4d1*a1**4*a2**4*d2**2*dy**4 &
+2.82d2*a1**3*a2**4*d2*dy**4+2.82d2*a1**4*a2**3*d2*dy**4 &
-4.59d2*a1**2*a2**4*dy**4-9.18d2*a1**3*a2**3*dy**4 &
-4.59d2*a1**4*a2**2*dy**4-8.4d1*a1**4*a2**4*d2*dx**4*dy**2 &
+7.98d2*a1**3*a2**4*dx**4*dy**2+7.98d2*a1**4*a2**3*dx**4*dy**2 &
+7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2-6.36d2*a1**3*a2**4*d2*dx**2*dy**2 &
-6.36d2*a1**4*a2**3*d2*dx**2*dy**2-4.08d2*a1**2*a2**4*dx**2*dy**2 &
-8.16d2*a1**3*a2**3*dx**2*dy**2-4.08d2*a1**4*a2**2*dx**2*dy**2 &
-3.6d1*a1**3*a2**4*d2**2*dy**2-3.6d1*a1**4*a2**3*d2**2*dy**2 &
+2.88d2*a1**2*a2**4*d2*dy**2+5.76d2*a1**3*a2**3*d2*dy**2 &
+2.88d2*a1**4*a2**2*d2*dy**2+1.35d2*a1*a2**4*dy**2 &
+4.05d2*a1**2*a2**3*dy**2+4.05d2*a1**3*a2**2*dy**2 &
+1.35d2*a1**4*a2*dy**2+4.2d1*a1**3*a2**4*d2*dx**4 &
+4.2d1*a1**4*a2**3*d2*dx**4-3.57d2*a1**2*a2**4*dx**4 &
-7.14d2*a1**3*a2**3*dx**4-3.57d2*a1**4*a2**2*dx**4 &
-3.6d1*a1**3*a2**4*d2**2*dx**2-3.6d1*a1**4*a2**3*d2**2*dx**2 &
+2.88d2*a1**2*a2**4*d2*dx**2+5.76d2*a1**3*a2**3*d2*dx**2 &
+2.88d2*a1**4*a2**2*d2*dx**2+1.35d2*a1*a2**4*dx**2 &
+4.05d2*a1**2*a2**3*dx**2+4.05d2*a1**3*a2**2*dx**2 &
+1.35d2*a1**4*a2*dx**2+3.6d1*a1**2*a2**4*d2**2+7.2d1*a1**3*a2**3*d2**2 &
+3.6d1*a1**4*a2**2*d2**2-3.0d2*a1*a2**4*d2-9.0d2*a1**2*a2**3*d2 &
-9.0d2*a1**3*a2**2*d2-3.0d2*a1**4*a2*d2+1.95d2*a2**4+7.8d2*a1*a2**3 &
+1.17d3*a1**2*a2**2+7.8d2*a1**3*a2+1.95d2*a1**4)*dz
                case (-1)
                  rlYlm_laplacian = &
-3.29728507019630958d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)* &
(5.6d1*a1**4*a2**4*d2*dy**8-5.32d2*a1**3*a2**4*dy**8 &
-5.32d2*a1**4*a2**3*dy**8-5.6d1*a1**4*a2**4*d2*dx**2*dy**6 &
+5.32d2*a1**3*a2**4*dx**2*dy**6+5.32d2*a1**4*a2**3*dx**2*dy**6 &
-8.8d1*a1**4*a2**4*d2**2*dy**6+7.4d2*a1**3*a2**4*d2*dy**6 &
+7.4d2*a1**4*a2**3*d2*dy**6+8.16d2*a1**2*a2**4*dy**6 &
+1.632d3*a1**3*a2**3*dy**6+8.16d2*a1**4*a2**2*dy**6 &
-2.8d2*a1**4*a2**4*d2*dx**4*dy**4+2.66d3*a1**3*a2**4*dx**4*dy**4 &
+2.66d3*a1**4*a2**3*dx**4*dy**4+1.76d2*a1**4*a2**4*d2**2*dx**2*dy**4 &
-1.732d3*a1**3*a2**4*d2*dx**2*dy**4-1.732d3*a1**4*a2**3*d2*dx**2*dy**4 &
+5.1d2*a1**2*a2**4*dx**2*dy**4+1.02d3*a1**3*a2**3*dx**2*dy**4 &
+5.1d2*a1**4*a2**2*dx**2*dy**4+3.2d1*a1**4*a2**4*d2**3*dy**4 &
-1.48d2*a1**3*a2**4*d2**2*dy**4-1.48d2*a1**4*a2**3*d2**2*dy**4 &
-1.38d3*a1**2*a2**4*d2*dy**4-2.76d3*a1**3*a2**3*d2*dy**4 &
-1.38d3*a1**4*a2**2*d2*dy**4+4.05d2*a1*a2**4*dy**4 &
+1.215d3*a1**2*a2**3*dy**4+1.215d3*a1**3*a2**2*dy**4 &
+4.05d2*a1**4*a2*dy**4-1.68d2*a1**4*a2**4*d2*dx**6*dy**2 &
+1.596d3*a1**3*a2**4*dx**6*dy**2+1.596d3*a1**4*a2**3*dx**6*dy**2 &
+2.64d2*a1**4*a2**4*d2**2*dx**4*dy**2 &
-2.388d3*a1**3*a2**4*d2*dx**4*dy**2-2.388d3*a1**4*a2**3*d2*dx**4*dy**2 &
-1.02d3*a1**2*a2**4*dx**4*dy**2-2.04d3*a1**3*a2**3*dx**4*dy**2 &
-1.02d3*a1**4*a2**2*dx**4*dy**2-9.6d1*a1**4*a2**4*d2**3*dx**2*dy**2 &
+8.4d2*a1**3*a2**4*d2**2*dx**2*dy**2 &
+8.4d2*a1**4*a2**3*d2**2*dx**2*dy**2+7.2d2*a1**2*a2**4*d2*dx**2*dy**2 &
+1.44d3*a1**3*a2**3*d2*dx**2*dy**2+7.2d2*a1**4*a2**2*d2*dx**2*dy**2 &
-8.1d2*a1*a2**4*dx**2*dy**2-2.43d3*a1**2*a2**3*dx**2*dy**2 &
-2.43d3*a1**3*a2**2*dx**2*dy**2-8.1d2*a1**4*a2*dx**2*dy**2 &
-4.8d1*a1**3*a2**4*d2**3*dy**2-4.8d1*a1**4*a2**3*d2**3*dy**2 &
+4.08d2*a1**2*a2**4*d2**2*dy**2+8.16d2*a1**3*a2**3*d2**2*dy**2 &
+4.08d2*a1**4*a2**2*d2**2*dy**2+3.6d1*a1*a2**4*d2*dy**2 &
+1.08d2*a1**2*a2**3*d2*dy**2+1.08d2*a1**3*a2**2*d2*dy**2 &
+3.6d1*a1**4*a2*d2*dy**2-2.34d2*a2**4*dy**2-9.36d2*a1*a2**3*dy**2 &
-1.404d3*a1**2*a2**2*dy**2-9.36d2*a1**3*a2*dy**2-2.34d2*a1**4*dy**2 &
+8.4d1*a1**3*a2**4*d2*dx**6+8.4d1*a1**4*a2**3*d2*dx**6 &
-7.14d2*a1**2*a2**4*dx**6-1.428d3*a1**3*a2**3*dx**6 &
-7.14d2*a1**4*a2**2*dx**6-1.32d2*a1**3*a2**4*d2**2*dx**4 &
-1.32d2*a1**4*a2**3*d2**2*dx**4+1.14d3*a1**2*a2**4*d2*dx**4 &
+2.28d3*a1**3*a2**3*d2*dx**4+1.14d3*a1**4*a2**2*d2*dx**4 &
-1.35d2*a1*a2**4*dx**4-4.05d2*a1**2*a2**3*dx**4 &
-4.05d2*a1**3*a2**2*dx**4-1.35d2*a1**4*a2*dx**4 &
+4.8d1*a1**3*a2**4*d2**3*dx**2+4.8d1*a1**4*a2**3*d2**3*dx**2 &
-4.08d2*a1**2*a2**4*d2**2*dx**2-8.16d2*a1**3*a2**3*d2**2*dx**2 &
-4.08d2*a1**4*a2**2*d2**2*dx**2-3.6d1*a1*a2**4*d2*dx**2 &
-1.08d2*a1**2*a2**3*d2*dx**2-1.08d2*a1**3*a2**2*d2*dx**2 &
-3.6d1*a1**4*a2*d2*dx**2+2.34d2*a2**4*dx**2+9.36d2*a1*a2**3*dx**2 &
+1.404d3*a1**2*a2**2*dx**2+9.36d2*a1**3*a2*dx**2+2.34d2*a1**4*dx**2)
                case (0)
                  rlYlm_laplacian = &
-2.08538618333770325d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(dy**2-3.0d0*dx**2)*(7.0d1*a1**3*a2**3*d2*dy**4 &
-6.65d2*a1**2*a2**3*dy**4-6.65d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.33d3*a1**2*a2**3*dx**2*dy**2 &
-1.33d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.6d2*a1**2*a2**3*d2*dy**2+6.6d2*a1**3*a2**2*d2*dy**2 &
+8.5d2*a1*a2**3*dy**2+1.7d3*a1**2*a2**2*dy**2+8.5d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-6.65d2*a1**2*a2**3*dx**4 &
-6.65d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.6d2*a1**2*a2**3*d2*dx**2+6.6d2*a1**3*a2**2*d2*dx**2 &
+8.5d2*a1*a2**3*dx**2+1.7d3*a1**2*a2**2*dx**2+8.5d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-4.0d1*a1**2*a2**3*d2**2 &
-4.0d1*a1**3*a2**2*d2**2-1.06d3*a1*a2**3*d2-2.12d3*a1**2*a2**2*d2 &
-1.06d3*a1**3*a2*d2+8.1d2*a2**3+2.43d3*a1*a2**2+2.43d3*a1**2*a2 &
+8.1d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian = &
-6.59457014039261917d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(2.8d1*a1**4*a2**4*d2*dy**6-2.66d2*a1**3*a2**4*dy**6 &
-2.66d2*a1**4*a2**3*dy**6-2.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
+2.66d2*a1**3*a2**4*dx**2*dy**4+2.66d2*a1**4*a2**3*dx**2*dy**4 &
-4.4d1*a1**4*a2**4*d2**2*dy**4+4.96d2*a1**3*a2**4*d2*dy**4 &
+4.96d2*a1**4*a2**3*d2*dy**4-6.63d2*a1**2*a2**4*dy**4 &
-1.326d3*a1**3*a2**3*dy**4-6.63d2*a1**4*a2**2*dy**4 &
-1.4d2*a1**4*a2**4*d2*dx**4*dy**2+1.33d3*a1**3*a2**4*dx**4*dy**2 &
+1.33d3*a1**4*a2**3*dx**4*dy**2+8.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-6.56d2*a1**3*a2**4*d2*dx**2*dy**2-6.56d2*a1**4*a2**3*d2*dx**2*dy**2 &
-1.53d3*a1**2*a2**4*dx**2*dy**2-3.06d3*a1**3*a2**3*dx**2*dy**2 &
-1.53d3*a1**4*a2**2*dx**2*dy**2+1.6d1*a1**4*a2**4*d2**3*dy**2 &
-2.72d2*a1**3*a2**4*d2**2*dy**2-2.72d2*a1**4*a2**3*d2**2*dy**2 &
+1.02d3*a1**2*a2**4*d2*dy**2+2.04d3*a1**3*a2**3*d2*dy**2 &
+1.02d3*a1**4*a2**2*d2*dy**2-8.4d1*a1**4*a2**4*d2*dx**6 &
+7.98d2*a1**3*a2**4*dx**6+7.98d2*a1**4*a2**3*dx**6 &
+1.32d2*a1**4*a2**4*d2**2*dx**4-1.152d3*a1**3*a2**4*d2*dx**4 &
-1.152d3*a1**4*a2**3*d2*dx**4-8.67d2*a1**2*a2**4*dx**4 &
-1.734d3*a1**3*a2**3*dx**4-8.67d2*a1**4*a2**2*dx**4 &
-4.8d1*a1**4*a2**4*d2**3*dx**2+2.88d2*a1**3*a2**4*d2**2*dx**2 &
+2.88d2*a1**4*a2**3*d2**2*dx**2+1.5d3*a1**2*a2**4*d2*dx**2 &
+3.0d3*a1**3*a2**3*d2*dx**2+1.5d3*a1**4*a2**2*d2*dx**2 &
-5.4d2*a1*a2**4*dx**2-1.62d3*a1**2*a2**3*dx**2 &
-1.62d3*a1**3*a2**2*dx**2-5.4d2*a1**4*a2*dx**2+4.8d1*a1**3*a2**4*d2**3 &
+4.8d1*a1**4*a2**3*d2**3-4.08d2*a1**2*a2**4*d2**2 &
-8.16d2*a1**3*a2**3*d2**2-4.08d2*a1**4*a2**2*d2**2-3.6d1*a1*a2**4*d2 &
-1.08d2*a1**2*a2**3*d2-1.08d2*a1**3*a2**2*d2-3.6d1*a1**4*a2*d2 &
+2.34d2*a2**4+9.36d2*a1*a2**3+1.404d3*a1**2*a2**2+9.36d2*a1**3*a2 &
+2.34d2*a1**4)
                case (2)
                  rlYlm_laplacian = &
-9.3261305305638875d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.33d2*a1**3*a2**4*dy**6 &
-1.33d2*a1**4*a2**3*dy**6-4.2d1*a1**4*a2**4*d2*dx**2*dy**4 &
+3.99d2*a1**3*a2**4*dx**2*dy**4+3.99d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+7.8d1*a1**3*a2**4*d2*dy**4 &
+7.8d1*a1**4*a2**3*d2*dy**4+3.06d2*a1**2*a2**4*dy**4 &
+6.12d2*a1**3*a2**3*dy**4+3.06d2*a1**4*a2**2*dy**4 &
-1.4d1*a1**4*a2**4*d2*dx**4*dy**2+1.33d2*a1**3*a2**4*dx**4*dy**2 &
+1.33d2*a1**4*a2**3*dx**4*dy**2+4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-5.64d2*a1**3*a2**4*d2*dx**2*dy**2-5.64d2*a1**4*a2**3*d2*dx**2*dy**2 &
+9.18d2*a1**2*a2**4*dx**2*dy**2+1.836d3*a1**3*a2**3*dx**2*dy**2 &
+9.18d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**3*a2**4*d2**2*dy**2 &
+3.6d1*a1**4*a2**3*d2**2*dy**2-2.88d2*a1**2*a2**4*d2*dy**2 &
-5.76d2*a1**3*a2**3*d2*dy**2-2.88d2*a1**4*a2**2*d2*dy**2 &
-1.35d2*a1*a2**4*dy**2-4.05d2*a1**2*a2**3*dy**2 &
-4.05d2*a1**3*a2**2*dy**2-1.35d2*a1**4*a2*dy**2 &
+4.2d1*a1**4*a2**4*d2*dx**6-3.99d2*a1**3*a2**4*dx**6 &
-3.99d2*a1**4*a2**3*dx**6-3.6d1*a1**4*a2**4*d2**2*dx**4 &
+3.18d2*a1**3*a2**4*d2*dx**4+3.18d2*a1**4*a2**3*d2*dx**4 &
+2.04d2*a1**2*a2**4*dx**4+4.08d2*a1**3*a2**3*dx**4 &
+2.04d2*a1**4*a2**2*dx**4+3.6d1*a1**3*a2**4*d2**2*dx**2 &
+3.6d1*a1**4*a2**3*d2**2*dx**2-2.88d2*a1**2*a2**4*d2*dx**2 &
-5.76d2*a1**3*a2**3*d2*dx**2-2.88d2*a1**4*a2**2*d2*dx**2 &
-1.35d2*a1*a2**4*dx**2-4.05d2*a1**2*a2**3*dx**2 &
-4.05d2*a1**3*a2**2*dx**2-1.35d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+3.0d2*a1*a2**4*d2+9.0d2*a1**2*a2**3*d2 &
+9.0d2*a1**3*a2**2*d2+3.0d2*a1**4*a2*d2-1.95d2*a2**4-7.8d2*a1*a2**3 &
-1.17d3*a1**2*a2**2-7.8d2*a1**3*a2-1.95d2*a1**4)*dz
                case (3)
                  rlYlm_laplacian = &
-1.74475925948511734d1*E*a1**4*a2**4*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(dy**2-3.0d0*dx**2)*(3.0d0*dy**2-1.0d0*dx**2)* &
(4.0d0*a1**2*a2**2*d2*dy**2-3.8d1*a1*a2**2*dy**2-3.8d1*a1**2*a2*dy**2 &
+4.0d0*a1**2*a2**2*d2*dx**2-3.8d1*a1*a2**2*dx**2-3.8d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+4.0d1*a1*a2**2*d2+4.0d1*a1**2*a2*d2 &
-1.7d1*a2**2-3.4d1*a1*a2-1.7d1*a1**2)
                case (4)
                  rlYlm_laplacian = &
-1.23373110391994556d1*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(2.0d0*a1**4*a2**4*d2*dy**6-1.9d1*a1**3*a2**4*dy**6 &
-1.9d1*a1**4*a2**3*dy**6-1.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.71d2*a1**3*a2**4*dx**2*dy**4+1.71d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**3*a2**4*d2*dy**4-1.2d1*a1**4*a2**3*d2*dy**4 &
+1.02d2*a1**2*a2**4*dy**4+2.04d2*a1**3*a2**3*dy**4 &
+1.02d2*a1**4*a2**2*dy**4+3.8d1*a1**4*a2**4*d2*dx**4*dy**2 &
-3.61d2*a1**3*a2**4*dx**4*dy**2-3.61d2*a1**4*a2**3*dx**4*dy**2 &
-2.4d1*a1**3*a2**4*d2*dx**2*dy**2-2.4d1*a1**4*a2**3*d2*dx**2*dy**2 &
+2.04d2*a1**2*a2**4*dx**2*dy**2+4.08d2*a1**3*a2**3*dx**2*dy**2 &
+2.04d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**2*a2**4*d2*dy**2 &
+7.2d1*a1**3*a2**3*d2*dy**2+3.6d1*a1**4*a2**2*d2*dy**2 &
-2.7d2*a1*a2**4*dy**2-8.1d2*a1**2*a2**3*dy**2-8.1d2*a1**3*a2**2*dy**2 &
-2.7d2*a1**4*a2*dy**2-6.0d0*a1**4*a2**4*d2*dx**6 &
+5.7d1*a1**3*a2**4*dx**6+5.7d1*a1**4*a2**3*dx**6 &
-1.2d1*a1**3*a2**4*d2*dx**4-1.2d1*a1**4*a2**3*d2*dx**4 &
+1.02d2*a1**2*a2**4*dx**4+2.04d2*a1**3*a2**3*dx**4 &
+1.02d2*a1**4*a2**2*dx**4+3.6d1*a1**2*a2**4*d2*dx**2 &
+7.2d1*a1**3*a2**3*d2*dx**2+3.6d1*a1**4*a2**2*d2*dx**2 &
-2.7d2*a1*a2**4*dx**2-8.1d2*a1**2*a2**3*dx**2-8.1d2*a1**3*a2**2*dx**2 &
-2.7d2*a1**4*a2*dx**2-2.4d1*a1*a2**4*d2-7.2d1*a1**2*a2**3*d2 &
-7.2d1*a1**3*a2**2*d2-2.4d1*a1**4*a2*d2+1.56d2*a2**4+6.24d2*a1*a2**3 &
+9.36d2*a1**2*a2**2+6.24d2*a1**3*a2+1.56d2*a1**4)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (-2)
          ! selection on l2: l1=4, m1=-2
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=-2, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-2.97249547320450826d0*E*a1*a2**5*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*dy*(7.0d0*(dy**2+dx**2)-6.0d0*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=-2, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-2.5742565924293503d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dx* &
(2.8d1*a1**2*a2**2*d2*dy**4-1.82d2*a1*a2**2*dy**4 &
-1.82d2*a1**2*a2*dy**4+2.8d1*a1**2*a2**2*d2*dx**2*dy**2 &
-1.82d2*a1*a2**2*dx**2*dy**2-1.82d2*a1**2*a2*dx**2*dy**2 &
-2.4d1*a1**2*a2**2*d2**2*dy**2+1.38d2*a1*a2**2*d2*dy**2 &
+1.38d2*a1**2*a2*d2*dy**2+9.9d1*a2**2*dy**2+1.98d2*a1*a2*dy**2 &
+9.9d1*a1**2*dy**2-1.4d1*a1*a2**2*d2*dx**2-1.4d1*a1**2*a2*d2*dx**2 &
+7.7d1*a2**2*dx**2+1.54d2*a1*a2*dx**2+7.7d1*a1**2*dx**2 &
+1.2d1*a1*a2**2*d2**2+1.2d1*a1**2*a2*d2**2-6.6d1*a2**2*d2 &
-1.32d2*a1*a2*d2-6.6d1*a1**2*d2)
                case (0)
                  rlYlm_laplacian = &
-5.14851318485870059d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(1.4d1*a1**2*a2**2*d2*dy**2-9.1d1*a1*a2**2*dy**2-9.1d1*a1**2*a2*dy**2 &
+1.4d1*a1**2*a2**2*d2*dx**2-9.1d1*a1*a2**2*dx**2-9.1d1*a1**2*a2*dx**2 &
-1.2d1*a1**2*a2**2*d2**2+9.0d1*a1*a2**2*d2+9.0d1*a1**2*a2*d2 &
-6.6d1*a2**2-1.32d2*a1*a2-6.6d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian = &
-2.5742565924293503d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dy* &
(2.8d1*a1**2*a2**2*d2*dx**2*dy**2-1.82d2*a1*a2**2*dx**2*dy**2 &
-1.82d2*a1**2*a2*dx**2*dy**2-1.4d1*a1*a2**2*d2*dy**2 &
-1.4d1*a1**2*a2*d2*dy**2+7.7d1*a2**2*dy**2+1.54d2*a1*a2*dy**2 &
+7.7d1*a1**2*dy**2+2.8d1*a1**2*a2**2*d2*dx**4-1.82d2*a1*a2**2*dx**4 &
-1.82d2*a1**2*a2*dx**4-2.4d1*a1**2*a2**2*d2**2*dx**2 &
+1.38d2*a1*a2**2*d2*dx**2+1.38d2*a1**2*a2*d2*dx**2+9.9d1*a2**2*dx**2 &
+1.98d2*a1*a2*dx**2+9.9d1*a1**2*dx**2+1.2d1*a1*a2**2*d2**2 &
+1.2d1*a1**2*a2*d2**2-6.6d1*a2**2*d2-1.32d2*a1*a2*d2-6.6d1*a1**2*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=-2, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-2.87810636609949888d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)* &
(5.6d1*a1**3*a2**3*d2*dx**2*dy**4-4.2d2*a1**2*a2**3*dx**2*dy**4 &
-4.2d2*a1**3*a2**2*dx**2*dy**4-2.8d1*a1**2*a2**3*d2*dy**4 &
-2.8d1*a1**3*a2**2*d2*dy**4+1.82d2*a1*a2**3*dy**4 &
+3.64d2*a1**2*a2**2*dy**4+1.82d2*a1**3*a2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**4*dy**2-4.2d2*a1**2*a2**3*dx**4*dy**2 &
-4.2d2*a1**3*a2**2*dx**4*dy**2-4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+2.88d2*a1**2*a2**3*d2*dx**2*dy**2+2.88d2*a1**3*a2**2*d2*dx**2*dy**2 &
+4.68d2*a1*a2**3*dx**2*dy**2+9.36d2*a1**2*a2**2*dx**2*dy**2 &
+4.68d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.38d2*a1*a2**3*d2*dy**2 &
-2.76d2*a1**2*a2**2*d2*dy**2-1.38d2*a1**3*a2*d2*dy**2 &
-9.9d1*a2**3*dy**2-2.97d2*a1*a2**2*dy**2-2.97d2*a1**2*a2*dy**2 &
-9.9d1*a1**3*dy**2-2.8d1*a1**2*a2**3*d2*dx**4 &
-2.8d1*a1**3*a2**2*d2*dx**4+1.82d2*a1*a2**3*dx**4 &
+3.64d2*a1**2*a2**2*dx**4+1.82d2*a1**3*a2*dx**4 &
+2.4d1*a1**2*a2**3*d2**2*dx**2+2.4d1*a1**3*a2**2*d2**2*dx**2 &
-1.38d2*a1*a2**3*d2*dx**2-2.76d2*a1**2*a2**2*d2*dx**2 &
-1.38d2*a1**3*a2*d2*dx**2-9.9d1*a2**3*dx**2-2.97d2*a1*a2**2*dx**2 &
-2.97d2*a1**2*a2*dx**2-9.9d1*a1**3*dx**2-1.2d1*a1*a2**3*d2**2 &
-2.4d1*a1**2*a2**2*d2**2-1.2d1*a1**3*a2*d2**2+6.6d1*a2**3*d2 &
+1.98d2*a1*a2**2*d2+1.98d2*a1**2*a2*d2+6.6d1*a1**3*d2)
                case (-1)
                  rlYlm_laplacian = &
-5.75621273219899775d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(2.8d1*a1**3*a2**3*d2*dy**4-2.1d2*a1**2*a2**3*dy**4 &
-2.1d2*a1**3*a2**2*dy**4+2.8d1*a1**3*a2**3*d2*dx**2*dy**2 &
-2.1d2*a1**2*a2**3*dx**2*dy**2-2.1d2*a1**3*a2**2*dx**2*dy**2 &
-2.4d1*a1**3*a2**3*d2**2*dy**2+1.86d2*a1**2*a2**3*d2*dy**2 &
+1.86d2*a1**3*a2**2*d2*dy**2-3.9d1*a1*a2**3*dy**2 &
-7.8d1*a1**2*a2**2*dy**2-3.9d1*a1**3*a2*dy**2 &
-1.4d1*a1**2*a2**3*d2*dx**2-1.4d1*a1**3*a2**2*d2*dx**2 &
+9.1d1*a1*a2**3*dx**2+1.82d2*a1**2*a2**2*dx**2+9.1d1*a1**3*a2*dx**2 &
+1.2d1*a1**2*a2**3*d2**2+1.2d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2 &
-1.8d2*a1**2*a2**2*d2-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian =3.32335097044784255d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(4.2d1*a1**3*a2**3*d2*dy**4 &
-3.15d2*a1**2*a2**3*dy**4-3.15d2*a1**3*a2**2*dy**4 &
+8.4d1*a1**3*a2**3*d2*dx**2*dy**2-6.3d2*a1**2*a2**3*dx**2*dy**2 &
-6.3d2*a1**3*a2**2*dx**2*dy**2-6.4d1*a1**3*a2**3*d2**2*dy**2 &
+4.96d2*a1**2*a2**3*d2*dy**2+4.96d2*a1**3*a2**2*d2*dy**2 &
-1.04d2*a1*a2**3*dy**2-2.08d2*a1**2*a2**2*dy**2-1.04d2*a1**3*a2*dy**2 &
+4.2d1*a1**3*a2**3*d2*dx**4-3.15d2*a1**2*a2**3*dx**4 &
-3.15d2*a1**3*a2**2*dx**4-6.4d1*a1**3*a2**3*d2**2*dx**2 &
+4.96d2*a1**2*a2**3*d2*dx**2+4.96d2*a1**3*a2**2*d2*dx**2 &
-1.04d2*a1*a2**3*dx**2-2.08d2*a1**2*a2**2*dx**2-1.04d2*a1**3*a2*dx**2 &
+2.4d1*a1**3*a2**3*d2**3-2.04d2*a1**2*a2**3*d2**2 &
-2.04d2*a1**3*a2**2*d2**2+1.74d2*a1*a2**3*d2+3.48d2*a1**2*a2**2*d2 &
+1.74d2*a1**3*a2*d2-9.9d1*a2**3-2.97d2*a1*a2**2-2.97d2*a1**2*a2 &
-9.9d1*a1**3)
                case (1)
                  rlYlm_laplacian = &
-5.75621273219899775d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(2.8d1*a1**3*a2**3*d2*dx**2*dy**2-2.1d2*a1**2*a2**3*dx**2*dy**2 &
-2.1d2*a1**3*a2**2*dx**2*dy**2-1.4d1*a1**2*a2**3*d2*dy**2 &
-1.4d1*a1**3*a2**2*d2*dy**2+9.1d1*a1*a2**3*dy**2 &
+1.82d2*a1**2*a2**2*dy**2+9.1d1*a1**3*a2*dy**2 &
+2.8d1*a1**3*a2**3*d2*dx**4-2.1d2*a1**2*a2**3*dx**4 &
-2.1d2*a1**3*a2**2*dx**4-2.4d1*a1**3*a2**3*d2**2*dx**2 &
+1.86d2*a1**2*a2**3*d2*dx**2+1.86d2*a1**3*a2**2*d2*dx**2 &
-3.9d1*a1*a2**3*dx**2-7.8d1*a1**2*a2**2*dx**2-3.9d1*a1**3*a2*dx**2 &
+1.2d1*a1**2*a2**3*d2**2+1.2d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2 &
-1.8d2*a1**2*a2**2*d2-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian =5.75621273219899775d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**2*a2**2*d2*dy**2-1.05d2*a1*a2**2*dy**2 &
-1.05d2*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-1.05d2*a1*a2**2*dx**2-1.05d2*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2 &
+8.6d1*a1*a2**2*d2+8.6d1*a1**2*a2*d2+2.6d1*a2**2+5.2d1*a1*a2 &
+2.6d1*a1**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=-2, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =3.10871017685462917d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(2.8d1*a1**4*a2**4*d2*dy**6 &
-2.38d2*a1**3*a2**4*dy**6-2.38d2*a1**4*a2**3*dy**6 &
-5.6d1*a1**4*a2**4*d2*dx**2*dy**4+4.76d2*a1**3*a2**4*dx**2*dy**4 &
+4.76d2*a1**4*a2**3*dx**2*dy**4-2.4d1*a1**4*a2**4*d2**2*dy**4 &
+2.34d2*a1**3*a2**4*d2*dy**4+2.34d2*a1**4*a2**3*d2*dy**4 &
-2.25d2*a1**2*a2**4*dy**4-4.5d2*a1**3*a2**3*dy**4 &
-2.25d2*a1**4*a2**2*dy**4-8.4d1*a1**4*a2**4*d2*dx**4*dy**2 &
+7.14d2*a1**3*a2**4*dx**4*dy**2+7.14d2*a1**4*a2**3*dx**4*dy**2 &
+7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2-4.92d2*a1**3*a2**4*d2*dx**2*dy**2 &
-4.92d2*a1**4*a2**3*d2*dx**2*dy**2-9.0d2*a1**2*a2**4*dx**2*dy**2 &
-1.8d3*a1**3*a2**3*dx**2*dy**2-9.0d2*a1**4*a2**2*dx**2*dy**2 &
-3.6d1*a1**3*a2**4*d2**2*dy**2-3.6d1*a1**4*a2**3*d2**2*dy**2 &
+2.16d2*a1**2*a2**4*d2*dy**2+4.32d2*a1**3*a2**3*d2*dy**2 &
+2.16d2*a1**4*a2**2*d2*dy**2+3.51d2*a1*a2**4*dy**2 &
+1.053d3*a1**2*a2**3*dy**2+1.053d3*a1**3*a2**2*dy**2 &
+3.51d2*a1**4*a2*dy**2+4.2d1*a1**3*a2**4*d2*dx**4 &
+4.2d1*a1**4*a2**3*d2*dx**4-3.15d2*a1**2*a2**4*dx**4 &
-6.3d2*a1**3*a2**3*dx**4-3.15d2*a1**4*a2**2*dx**4 &
-3.6d1*a1**3*a2**4*d2**2*dx**2-3.6d1*a1**4*a2**3*d2**2*dx**2 &
+2.16d2*a1**2*a2**4*d2*dx**2+4.32d2*a1**3*a2**3*d2*dx**2 &
+2.16d2*a1**4*a2**2*d2*dx**2+3.51d2*a1*a2**4*dx**2 &
+1.053d3*a1**2*a2**3*dx**2+1.053d3*a1**3*a2**2*dx**2 &
+3.51d2*a1**4*a2*dx**2+3.6d1*a1**2*a2**4*d2**2+7.2d1*a1**3*a2**3*d2**2 &
+3.6d1*a1**4*a2**2*d2**2-2.28d2*a1*a2**4*d2-6.84d2*a1**2*a2**3*d2 &
-6.84d2*a1**3*a2**2*d2-2.28d2*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (-2)
                  rlYlm_laplacian = &
-7.6147536914910937d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(5.6d1*a1**4*a2**4*d2*dx**2*dy**4-4.76d2*a1**3*a2**4*dx**2*dy**4 &
-4.76d2*a1**4*a2**3*dx**2*dy**4-2.8d1*a1**3*a2**4*d2*dy**4 &
-2.8d1*a1**4*a2**3*d2*dy**4+2.1d2*a1**2*a2**4*dy**4 &
+4.2d2*a1**3*a2**3*dy**4+2.1d2*a1**4*a2**2*dy**4 &
+5.6d1*a1**4*a2**4*d2*dx**4*dy**2-4.76d2*a1**3*a2**4*dx**4*dy**2 &
-4.76d2*a1**4*a2**3*dx**4*dy**2-4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+3.84d2*a1**3*a2**4*d2*dx**2*dy**2+3.84d2*a1**4*a2**3*d2*dx**2*dy**2 &
+1.8d2*a1**2*a2**4*dx**2*dy**2+3.6d2*a1**3*a2**3*dx**2*dy**2 &
+1.8d2*a1**4*a2**2*dx**2*dy**2+2.4d1*a1**3*a2**4*d2**2*dy**2 &
+2.4d1*a1**4*a2**3*d2**2*dy**2-1.86d2*a1**2*a2**4*d2*dy**2 &
-3.72d2*a1**3*a2**3*d2*dy**2-1.86d2*a1**4*a2**2*d2*dy**2 &
+3.9d1*a1*a2**4*dy**2+1.17d2*a1**2*a2**3*dy**2 &
+1.17d2*a1**3*a2**2*dy**2+3.9d1*a1**4*a2*dy**2 &
-2.8d1*a1**3*a2**4*d2*dx**4-2.8d1*a1**4*a2**3*d2*dx**4 &
+2.1d2*a1**2*a2**4*dx**4+4.2d2*a1**3*a2**3*dx**4 &
+2.1d2*a1**4*a2**2*dx**4+2.4d1*a1**3*a2**4*d2**2*dx**2 &
+2.4d1*a1**4*a2**3*d2**2*dx**2-1.86d2*a1**2*a2**4*d2*dx**2 &
-3.72d2*a1**3*a2**3*d2*dx**2-1.86d2*a1**4*a2**2*d2*dx**2 &
+3.9d1*a1*a2**4*dx**2+1.17d2*a1**2*a2**3*dx**2 &
+1.17d2*a1**3*a2**2*dx**2+3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2 &
-2.4d1*a1**3*a2**3*d2**2-1.2d1*a1**4*a2**2*d2**2+9.0d1*a1*a2**4*d2 &
+2.7d2*a1**2*a2**3*d2+2.7d2*a1**3*a2**2*d2+9.0d1*a1**4*a2*d2 &
-6.6d1*a2**4-2.64d2*a1*a2**3-3.96d2*a1**2*a2**2-2.64d2*a1**3*a2 &
-6.6d1*a1**4)*dz
                case (-1)
                  rlYlm_laplacian =2.40799654862869848d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(1.4d2*a1**4*a2**4*d2*dy**6 &
-1.19d3*a1**3*a2**4*dy**6-1.19d3*a1**4*a2**3*dy**6 &
+2.8d2*a1**4*a2**4*d2*dx**2*dy**4-2.38d3*a1**3*a2**4*dx**2*dy**4 &
-2.38d3*a1**4*a2**3*dx**2*dy**4-2.32d2*a1**4*a2**4*d2**2*dy**4 &
+2.01d3*a1**3*a2**4*d2*dy**4+2.01d3*a1**4*a2**3*d2*dy**4 &
-2.85d2*a1**2*a2**4*dy**4-5.7d2*a1**3*a2**3*dy**4 &
-2.85d2*a1**4*a2**2*dy**4+1.4d2*a1**4*a2**4*d2*dx**4*dy**2 &
-1.19d3*a1**3*a2**4*dx**4*dy**2-1.19d3*a1**4*a2**3*dx**4*dy**2 &
-2.32d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.94d3*a1**3*a2**4*d2*dx**2*dy**2+1.94d3*a1**4*a2**3*d2*dx**2*dy**2 &
+2.4d2*a1**2*a2**4*dx**2*dy**2+4.8d2*a1**3*a2**3*dx**2*dy**2 &
+2.4d2*a1**4*a2**2*dx**2*dy**2+9.6d1*a1**4*a2**4*d2**3*dy**2 &
-8.28d2*a1**3*a2**4*d2**2*dy**2-8.28d2*a1**4*a2**3*d2**2*dy**2 &
+8.4d1*a1**2*a2**4*d2*dy**2+1.68d2*a1**3*a2**3*d2*dy**2 &
+8.4d1*a1**4*a2**2*d2*dy**2+3.9d1*a1*a2**4*dy**2 &
+1.17d2*a1**2*a2**3*dy**2+1.17d2*a1**3*a2**2*dy**2 &
+3.9d1*a1**4*a2*dy**2-7.0d1*a1**3*a2**4*d2*dx**4 &
-7.0d1*a1**4*a2**3*d2*dx**4+5.25d2*a1**2*a2**4*dx**4 &
+1.05d3*a1**3*a2**3*dx**4+5.25d2*a1**4*a2**2*dx**4 &
+1.16d2*a1**3*a2**4*d2**2*dx**2+1.16d2*a1**4*a2**3*d2**2*dx**2 &
-9.48d2*a1**2*a2**4*d2*dx**2-1.896d3*a1**3*a2**3*d2*dx**2 &
-9.48d2*a1**4*a2**2*d2*dx**2+5.07d2*a1*a2**4*dx**2 &
+1.521d3*a1**2*a2**3*dx**2+1.521d3*a1**3*a2**2*dx**2 &
+5.07d2*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+4.44d2*a1**2*a2**4*d2**2+8.88d2*a1**3*a2**3*d2**2 &
+4.44d2*a1**4*a2**2*d2**2-5.76d2*a1*a2**4*d2-1.728d3*a1**2*a2**3*d2 &
-1.728d3*a1**3*a2**2*d2-5.76d2*a1**4*a2*d2+1.65d2*a2**4+6.6d2*a1*a2**3 &
+9.9d2*a1**2*a2**2+6.6d2*a1**3*a2+1.65d2*a1**4)
                case (0)
                  rlYlm_laplacian =3.93224189768219417d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.95d2*a1**2*a2**3*dy**4-5.95d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.19d3*a1**2*a2**3*dx**2*dy**2 &
-1.19d3*a1**3*a2**2*dx**2*dy**2-8.8d1*a1**3*a2**3*d2**2*dy**2 &
+7.6d2*a1**2*a2**3*d2*dy**2+7.6d2*a1**3*a2**2*d2*dy**2 &
-9.0d1*a1*a2**3*dy**2-1.8d2*a1**2*a2**2*dy**2-9.0d1*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.8d1*a1**3*a2**3*d2**2*dx**2 &
+7.6d2*a1**2*a2**3*d2*dx**2+7.6d2*a1**3*a2**2*d2*dx**2 &
-9.0d1*a1*a2**3*dx**2-1.8d2*a1**2*a2**2*dx**2-9.0d1*a1**3*a2*dx**2 &
+2.4d1*a1**3*a2**3*d2**3-2.04d2*a1**2*a2**3*d2**2 &
-2.04d2*a1**3*a2**2*d2**2-1.8d1*a1*a2**3*d2-3.6d1*a1**2*a2**2*d2 &
-1.8d1*a1**3*a2*d2+1.17d2*a2**3+3.51d2*a1*a2**2+3.51d2*a1**2*a2 &
+1.17d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian =2.40799654862869848d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(1.4d2*a1**4*a2**4*d2*dx**2*dy**4 &
-1.19d3*a1**3*a2**4*dx**2*dy**4-1.19d3*a1**4*a2**3*dx**2*dy**4 &
-7.0d1*a1**3*a2**4*d2*dy**4-7.0d1*a1**4*a2**3*d2*dy**4 &
+5.25d2*a1**2*a2**4*dy**4+1.05d3*a1**3*a2**3*dy**4 &
+5.25d2*a1**4*a2**2*dy**4+2.8d2*a1**4*a2**4*d2*dx**4*dy**2 &
-2.38d3*a1**3*a2**4*dx**4*dy**2-2.38d3*a1**4*a2**3*dx**4*dy**2 &
-2.32d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.94d3*a1**3*a2**4*d2*dx**2*dy**2+1.94d3*a1**4*a2**3*d2*dx**2*dy**2 &
+2.4d2*a1**2*a2**4*dx**2*dy**2+4.8d2*a1**3*a2**3*dx**2*dy**2 &
+2.4d2*a1**4*a2**2*dx**2*dy**2+1.16d2*a1**3*a2**4*d2**2*dy**2 &
+1.16d2*a1**4*a2**3*d2**2*dy**2-9.48d2*a1**2*a2**4*d2*dy**2 &
-1.896d3*a1**3*a2**3*d2*dy**2-9.48d2*a1**4*a2**2*d2*dy**2 &
+5.07d2*a1*a2**4*dy**2+1.521d3*a1**2*a2**3*dy**2 &
+1.521d3*a1**3*a2**2*dy**2+5.07d2*a1**4*a2*dy**2 &
+1.4d2*a1**4*a2**4*d2*dx**6-1.19d3*a1**3*a2**4*dx**6 &
-1.19d3*a1**4*a2**3*dx**6-2.32d2*a1**4*a2**4*d2**2*dx**4 &
+2.01d3*a1**3*a2**4*d2*dx**4+2.01d3*a1**4*a2**3*d2*dx**4 &
-2.85d2*a1**2*a2**4*dx**4-5.7d2*a1**3*a2**3*dx**4 &
-2.85d2*a1**4*a2**2*dx**4+9.6d1*a1**4*a2**4*d2**3*dx**2 &
-8.28d2*a1**3*a2**4*d2**2*dx**2-8.28d2*a1**4*a2**3*d2**2*dx**2 &
+8.4d1*a1**2*a2**4*d2*dx**2+1.68d2*a1**3*a2**3*d2*dx**2 &
+8.4d1*a1**4*a2**2*d2*dx**2+3.9d1*a1*a2**4*dx**2 &
+1.17d2*a1**2*a2**3*dx**2+1.17d2*a1**3*a2**2*dx**2 &
+3.9d1*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+4.44d2*a1**2*a2**4*d2**2+8.88d2*a1**3*a2**3*d2**2 &
+4.44d2*a1**4*a2**2*d2**2-5.76d2*a1*a2**4*d2-1.728d3*a1**2*a2**3*d2 &
-1.728d3*a1**3*a2**2*d2-5.76d2*a1**4*a2*d2+1.65d2*a2**4+6.6d2*a1*a2**3 &
+9.9d2*a1**2*a2**2+6.6d2*a1**3*a2+1.65d2*a1**4)
                case (2)
                  rlYlm_laplacian =1.52295073829821874d1*E*a1**3*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**2*a2**2*d2*dy**2-1.19d2*a1*a2**2*dy**2 &
-1.19d2*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-1.19d2*a1*a2**2*dx**2-1.19d2*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2 &
+1.1d2*a1*a2**2*d2+1.1d2*a1**2*a2*d2-6.0d1*a2**2-1.2d2*a1*a2 &
-6.0d1*a1**2)*dz
                case (3)
                  rlYlm_laplacian =3.10871017685462917d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(8.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
-7.14d2*a1**3*a2**4*dx**2*dy**4-7.14d2*a1**4*a2**3*dx**2*dy**4 &
-4.2d1*a1**3*a2**4*d2*dy**4-4.2d1*a1**4*a2**3*d2*dy**4 &
+3.15d2*a1**2*a2**4*dy**4+6.3d2*a1**3*a2**3*dy**4 &
+3.15d2*a1**4*a2**2*dy**4+5.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
-4.76d2*a1**3*a2**4*dx**4*dy**2-4.76d2*a1**4*a2**3*dx**4*dy**2 &
-7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2+4.92d2*a1**3*a2**4*d2*dx**2*dy**2 &
+4.92d2*a1**4*a2**3*d2*dx**2*dy**2+9.0d2*a1**2*a2**4*dx**2*dy**2 &
+1.8d3*a1**3*a2**3*dx**2*dy**2+9.0d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**3*a2**4*d2**2*dy**2+3.6d1*a1**4*a2**3*d2**2*dy**2 &
-2.16d2*a1**2*a2**4*d2*dy**2-4.32d2*a1**3*a2**3*d2*dy**2 &
-2.16d2*a1**4*a2**2*d2*dy**2-3.51d2*a1*a2**4*dy**2 &
-1.053d3*a1**2*a2**3*dy**2-1.053d3*a1**3*a2**2*dy**2 &
-3.51d2*a1**4*a2*dy**2-2.8d1*a1**4*a2**4*d2*dx**6 &
+2.38d2*a1**3*a2**4*dx**6+2.38d2*a1**4*a2**3*dx**6 &
+2.4d1*a1**4*a2**4*d2**2*dx**4-2.34d2*a1**3*a2**4*d2*dx**4 &
-2.34d2*a1**4*a2**3*d2*dx**4+2.25d2*a1**2*a2**4*dx**4 &
+4.5d2*a1**3*a2**3*dx**4+2.25d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**3*a2**4*d2**2*dx**2+3.6d1*a1**4*a2**3*d2**2*dx**2 &
-2.16d2*a1**2*a2**4*d2*dx**2-4.32d2*a1**3*a2**3*d2*dx**2 &
-2.16d2*a1**4*a2**2*d2*dx**2-3.51d2*a1*a2**4*dx**2 &
-1.053d3*a1**2*a2**3*dx**2-1.053d3*a1**3*a2**2*dx**2 &
-3.51d2*a1**4*a2*dx**2-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.28d2*a1*a2**4*d2+6.84d2*a1**2*a2**3*d2 &
+6.84d2*a1**3*a2**2*d2+2.28d2*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=-2, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =6.59457014039261917d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(dy+dx)*(dy-1.0d0*dx)* &
(5.6d1*a1**4*a2**4*d2*dx**2*dy**4-5.32d2*a1**3*a2**4*dx**2*dy**4 &
-5.32d2*a1**4*a2**3*dx**2*dy**4-2.8d1*a1**3*a2**4*d2*dy**4 &
-2.8d1*a1**4*a2**3*d2*dy**4+2.38d2*a1**2*a2**4*dy**4 &
+4.76d2*a1**3*a2**3*dy**4+2.38d2*a1**4*a2**2*dy**4 &
+5.6d1*a1**4*a2**4*d2*dx**4*dy**2-5.32d2*a1**3*a2**4*dx**4*dy**2 &
-5.32d2*a1**4*a2**3*dx**4*dy**2-4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+3.68d2*a1**3*a2**4*d2*dx**2*dy**2+3.68d2*a1**4*a2**3*d2*dx**2*dy**2 &
+7.48d2*a1**2*a2**4*dx**2*dy**2+1.496d3*a1**3*a2**3*dx**2*dy**2 &
+7.48d2*a1**4*a2**2*dx**2*dy**2+2.4d1*a1**3*a2**4*d2**2*dy**2 &
+2.4d1*a1**4*a2**3*d2**2*dy**2-1.5d2*a1**2*a2**4*d2*dy**2 &
-3.0d2*a1**3*a2**3*d2*dy**2-1.5d2*a1**4*a2**2*d2*dy**2 &
-4.05d2*a1*a2**4*dy**2-1.215d3*a1**2*a2**3*dy**2 &
-1.215d3*a1**3*a2**2*dy**2-4.05d2*a1**4*a2*dy**2 &
-2.8d1*a1**3*a2**4*d2*dx**4-2.8d1*a1**4*a2**3*d2*dx**4 &
+2.38d2*a1**2*a2**4*dx**4+4.76d2*a1**3*a2**3*dx**4 &
+2.38d2*a1**4*a2**2*dx**4+2.4d1*a1**3*a2**4*d2**2*dx**2 &
+2.4d1*a1**4*a2**3*d2**2*dx**2-1.5d2*a1**2*a2**4*d2*dx**2 &
-3.0d2*a1**3*a2**3*d2*dx**2-1.5d2*a1**4*a2**2*d2*dx**2 &
-4.05d2*a1*a2**4*dx**2-1.215d3*a1**2*a2**3*dx**2 &
-1.215d3*a1**3*a2**2*dx**2-4.05d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.58d2*a1*a2**4*d2+7.74d2*a1**2*a2**3*d2 &
+7.74d2*a1**3*a2**2*d2+2.58d2*a1**4*a2*d2+7.8d1*a2**4+3.12d2*a1*a2**3 &
+4.68d2*a1**2*a2**2+3.12d2*a1**3*a2+7.8d1*a1**4)
                case (-3)
                  rlYlm_laplacian =9.3261305305638875d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dx*(2.8d1*a1**4*a2**4*d2*dy**6 &
-2.66d2*a1**3*a2**4*dy**6-2.66d2*a1**4*a2**3*dy**6 &
-5.6d1*a1**4*a2**4*d2*dx**2*dy**4+5.32d2*a1**3*a2**4*dx**2*dy**4 &
+5.32d2*a1**4*a2**3*dx**2*dy**4-2.4d1*a1**4*a2**4*d2**2*dy**4 &
+2.82d2*a1**3*a2**4*d2*dy**4+2.82d2*a1**4*a2**3*d2*dy**4 &
-4.59d2*a1**2*a2**4*dy**4-9.18d2*a1**3*a2**3*dy**4 &
-4.59d2*a1**4*a2**2*dy**4-8.4d1*a1**4*a2**4*d2*dx**4*dy**2 &
+7.98d2*a1**3*a2**4*dx**4*dy**2+7.98d2*a1**4*a2**3*dx**4*dy**2 &
+7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2-6.36d2*a1**3*a2**4*d2*dx**2*dy**2 &
-6.36d2*a1**4*a2**3*d2*dx**2*dy**2-4.08d2*a1**2*a2**4*dx**2*dy**2 &
-8.16d2*a1**3*a2**3*dx**2*dy**2-4.08d2*a1**4*a2**2*dx**2*dy**2 &
-3.6d1*a1**3*a2**4*d2**2*dy**2-3.6d1*a1**4*a2**3*d2**2*dy**2 &
+2.88d2*a1**2*a2**4*d2*dy**2+5.76d2*a1**3*a2**3*d2*dy**2 &
+2.88d2*a1**4*a2**2*d2*dy**2+1.35d2*a1*a2**4*dy**2 &
+4.05d2*a1**2*a2**3*dy**2+4.05d2*a1**3*a2**2*dy**2 &
+1.35d2*a1**4*a2*dy**2+4.2d1*a1**3*a2**4*d2*dx**4 &
+4.2d1*a1**4*a2**3*d2*dx**4-3.57d2*a1**2*a2**4*dx**4 &
-7.14d2*a1**3*a2**3*dx**4-3.57d2*a1**4*a2**2*dx**4 &
-3.6d1*a1**3*a2**4*d2**2*dx**2-3.6d1*a1**4*a2**3*d2**2*dx**2 &
+2.88d2*a1**2*a2**4*d2*dx**2+5.76d2*a1**3*a2**3*d2*dx**2 &
+2.88d2*a1**4*a2**2*d2*dx**2+1.35d2*a1*a2**4*dx**2 &
+4.05d2*a1**2*a2**3*dx**2+4.05d2*a1**3*a2**2*dx**2 &
+1.35d2*a1**4*a2*dx**2+3.6d1*a1**2*a2**4*d2**2+7.2d1*a1**3*a2**3*d2**2 &
+3.6d1*a1**4*a2**2*d2**2-3.0d2*a1*a2**4*d2-9.0d2*a1**2*a2**3*d2 &
-9.0d2*a1**3*a2**2*d2-3.0d2*a1**4*a2*d2+1.95d2*a2**4+7.8d2*a1*a2**3 &
+1.17d3*a1**2*a2**2+7.8d2*a1**3*a2+1.95d2*a1**4)*dz
                case (-2)
                  rlYlm_laplacian =2.49251322783588191d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(3.92d2*a1**5*a2**5*d2*dx**2*dy**6 &
-3.724d3*a1**4*a2**5*dx**2*dy**6-3.724d3*a1**5*a2**4*dx**2*dy**6 &
-1.96d2*a1**4*a2**5*d2*dy**6-1.96d2*a1**5*a2**4*d2*dy**6 &
+1.666d3*a1**3*a2**5*dy**6+3.332d3*a1**4*a2**4*dy**6 &
+1.666d3*a1**5*a2**3*dy**6+7.84d2*a1**5*a2**5*d2*dx**4*dy**4 &
-7.448d3*a1**4*a2**5*dx**4*dy**4-7.448d3*a1**5*a2**4*dx**4*dy**4 &
-6.72d2*a1**5*a2**5*d2**2*dx**2*dy**4 &
+6.132d3*a1**4*a2**5*d2*dx**2*dy**4+6.132d3*a1**5*a2**4*d2*dx**2*dy**4 &
+2.142d3*a1**3*a2**5*dx**2*dy**4+4.284d3*a1**4*a2**4*dx**2*dy**4 &
+2.142d3*a1**5*a2**3*dx**2*dy**4+3.36d2*a1**4*a2**5*d2**2*dy**4 &
+3.36d2*a1**5*a2**4*d2**2*dy**4-2.982d3*a1**3*a2**5*d2*dy**4 &
-5.964d3*a1**4*a2**4*d2*dy**4-2.982d3*a1**5*a2**3*d2*dy**4 &
+9.45d2*a1**2*a2**5*dy**4+2.835d3*a1**3*a2**4*dy**4 &
+2.835d3*a1**4*a2**3*dy**4+9.45d2*a1**5*a2**2*dy**4 &
+3.92d2*a1**5*a2**5*d2*dx**6*dy**2-3.724d3*a1**4*a2**5*dx**6*dy**2 &
-3.724d3*a1**5*a2**4*dx**6*dy**2-6.72d2*a1**5*a2**5*d2**2*dx**4*dy**2 &
+6.132d3*a1**4*a2**5*d2*dx**4*dy**2+6.132d3*a1**5*a2**4*d2*dx**4*dy**2 &
+2.142d3*a1**3*a2**5*dx**4*dy**2+4.284d3*a1**4*a2**4*dx**4*dy**2 &
+2.142d3*a1**5*a2**3*dx**4*dy**2+2.88d2*a1**5*a2**5*d2**3*dx**2*dy**2 &
-2.448d3*a1**4*a2**5*d2**2*dx**2*dy**2 &
-2.448d3*a1**5*a2**4*d2**2*dx**2*dy**2 &
-2.484d3*a1**3*a2**5*d2*dx**2*dy**2-4.968d3*a1**4*a2**4*d2*dx**2*dy**2 &
-2.484d3*a1**5*a2**3*d2*dx**2*dy**2+2.7d2*a1**2*a2**5*dx**2*dy**2 &
+8.1d2*a1**3*a2**4*dx**2*dy**2+8.1d2*a1**4*a2**3*dx**2*dy**2 &
+2.7d2*a1**5*a2**2*dx**2*dy**2-1.44d2*a1**4*a2**5*d2**3*dy**2 &
-1.44d2*a1**5*a2**4*d2**3*dy**2+1.296d3*a1**3*a2**5*d2**2*dy**2 &
+2.592d3*a1**4*a2**4*d2**2*dy**2+1.296d3*a1**5*a2**3*d2**2*dy**2 &
-4.92d2*a1**2*a2**5*d2*dy**2-1.476d3*a1**3*a2**4*d2*dy**2 &
-1.476d3*a1**4*a2**3*d2*dy**2-4.92d2*a1**5*a2**2*d2*dy**2 &
-3.12d2*a1*a2**5*dy**2-1.248d3*a1**2*a2**4*dy**2 &
-1.872d3*a1**3*a2**3*dy**2-1.248d3*a1**4*a2**2*dy**2 &
-3.12d2*a1**5*a2*dy**2-1.96d2*a1**4*a2**5*d2*dx**6 &
-1.96d2*a1**5*a2**4*d2*dx**6+1.666d3*a1**3*a2**5*dx**6 &
+3.332d3*a1**4*a2**4*dx**6+1.666d3*a1**5*a2**3*dx**6 &
+3.36d2*a1**4*a2**5*d2**2*dx**4+3.36d2*a1**5*a2**4*d2**2*dx**4 &
-2.982d3*a1**3*a2**5*d2*dx**4-5.964d3*a1**4*a2**4*d2*dx**4 &
-2.982d3*a1**5*a2**3*d2*dx**4+9.45d2*a1**2*a2**5*dx**4 &
+2.835d3*a1**3*a2**4*dx**4+2.835d3*a1**4*a2**3*dx**4 &
+9.45d2*a1**5*a2**2*dx**4-1.44d2*a1**4*a2**5*d2**3*dx**2 &
-1.44d2*a1**5*a2**4*d2**3*dx**2+1.296d3*a1**3*a2**5*d2**2*dx**2 &
+2.592d3*a1**4*a2**4*d2**2*dx**2+1.296d3*a1**5*a2**3*d2**2*dx**2 &
-4.92d2*a1**2*a2**5*d2*dx**2-1.476d3*a1**3*a2**4*d2*dx**2 &
-1.476d3*a1**4*a2**3*d2*dx**2-4.92d2*a1**5*a2**2*d2*dx**2 &
-3.12d2*a1*a2**5*dx**2-1.248d3*a1**2*a2**4*dx**2 &
-1.872d3*a1**3*a2**3*dx**2-1.248d3*a1**4*a2**2*dx**2 &
-3.12d2*a1**5*a2*dx**2+7.2d1*a1**3*a2**5*d2**3 &
+1.44d2*a1**4*a2**4*d2**3+7.2d1*a1**5*a2**3*d2**3 &
-6.84d2*a1**2*a2**5*d2**2-2.052d3*a1**3*a2**4*d2**2 &
-2.052d3*a1**4*a2**3*d2**2-6.84d2*a1**5*a2**2*d2**2+9.78d2*a1*a2**5*d2 &
+3.912d3*a1**2*a2**4*d2+5.868d3*a1**3*a2**3*d2+3.912d3*a1**4*a2**2*d2 &
+9.78d2*a1**5*a2*d2-2.31d2*a2**5-1.155d3*a1*a2**4-2.31d3*a1**2*a2**3 &
-2.31d3*a1**3*a2**2-1.155d3*a1**4*a2-2.31d2*a1**5)
                case (-1)
                  rlYlm_laplacian =3.52494601119984446d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dx*(1.96d2*a1**4*a2**4*d2*dy**6 &
-1.862d3*a1**3*a2**4*dy**6-1.862d3*a1**4*a2**3*dy**6 &
+3.92d2*a1**4*a2**4*d2*dx**2*dy**4-3.724d3*a1**3*a2**4*dx**2*dy**4 &
-3.724d3*a1**4*a2**3*dx**2*dy**4-2.8d2*a1**4*a2**4*d2**2*dy**4 &
+2.702d3*a1**3*a2**4*d2*dy**4+2.702d3*a1**4*a2**3*d2*dy**4 &
-3.57d2*a1**2*a2**4*dy**4-7.14d2*a1**3*a2**3*dy**4 &
-3.57d2*a1**4*a2**2*dy**4+1.96d2*a1**4*a2**4*d2*dx**4*dy**2 &
-1.862d3*a1**3*a2**4*dx**4*dy**2-1.862d3*a1**4*a2**3*dx**4*dy**2 &
-2.8d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+2.604d3*a1**3*a2**4*d2*dx**2*dy**2+2.604d3*a1**4*a2**3*d2*dx**2*dy**2 &
+4.76d2*a1**2*a2**4*dx**2*dy**2+9.52d2*a1**3*a2**3*dx**2*dy**2 &
+4.76d2*a1**4*a2**2*dx**2*dy**2+9.6d1*a1**4*a2**4*d2**3*dy**2 &
-9.0d2*a1**3*a2**4*d2**2*dy**2-9.0d2*a1**4*a2**3*d2**2*dy**2 &
-1.56d2*a1**2*a2**4*d2*dy**2-3.12d2*a1**3*a2**3*d2*dy**2 &
-1.56d2*a1**4*a2**2*d2*dy**2+4.05d2*a1*a2**4*dy**2 &
+1.215d3*a1**2*a2**3*dy**2+1.215d3*a1**3*a2**2*dy**2 &
+4.05d2*a1**4*a2*dy**2-9.8d1*a1**3*a2**4*d2*dx**4 &
-9.8d1*a1**4*a2**3*d2*dx**4+8.33d2*a1**2*a2**4*dx**4 &
+1.666d3*a1**3*a2**3*dx**4+8.33d2*a1**4*a2**2*dx**4 &
+1.4d2*a1**3*a2**4*d2**2*dx**2+1.4d2*a1**4*a2**3*d2**2*dx**2 &
-1.316d3*a1**2*a2**4*d2*dx**2-2.632d3*a1**3*a2**3*d2*dx**2 &
-1.316d3*a1**4*a2**2*d2*dx**2+9.45d2*a1*a2**4*dx**2 &
+2.835d3*a1**2*a2**3*dx**2+2.835d3*a1**3*a2**2*dx**2 &
+9.45d2*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+5.16d2*a1**2*a2**4*d2**2+1.032d3*a1**3*a2**3*d2**2 &
+5.16d2*a1**4*a2**2*d2**2-8.64d2*a1*a2**4*d2-2.592d3*a1**2*a2**3*d2 &
-2.592d3*a1**3*a2**2*d2-8.64d2*a1**4*a2*d2+3.51d2*a2**4 &
+1.404d3*a1*a2**3+2.106d3*a1**2*a2**2+1.404d3*a1**3*a2+3.51d2*a1**4 &
)*dz
                case (0)
                  rlYlm_laplacian = &
-1.1146858024516906d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(4.9d2*a1**4*a2**4*d2*dy**6-4.655d3*a1**3*a2**4*dy**6 &
-4.655d3*a1**4*a2**3*dy**6+1.47d3*a1**4*a2**4*d2*dx**2*dy**4 &
-1.3965d4*a1**3*a2**4*dx**2*dy**4-1.3965d4*a1**4*a2**3*dx**2*dy**4 &
-9.8d2*a1**4*a2**4*d2**2*dy**4+9.31d3*a1**3*a2**4*d2*dy**4 &
+9.31d3*a1**4*a2**3*d2*dy**4+1.47d3*a1**4*a2**4*d2*dx**4*dy**2 &
-1.3965d4*a1**3*a2**4*dx**4*dy**2-1.3965d4*a1**4*a2**3*dx**4*dy**2 &
-1.96d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.862d4*a1**3*a2**4*d2*dx**2*dy**2+1.862d4*a1**4*a2**3*d2*dx**2*dy**2 &
+5.92d2*a1**4*a2**4*d2**3*dy**2-5.48d3*a1**3*a2**4*d2**2*dy**2 &
-5.48d3*a1**4*a2**3*d2**2*dy**2-1.62d3*a1**2*a2**4*d2*dy**2 &
-3.24d3*a1**3*a2**3*d2*dy**2-1.62d3*a1**4*a2**2*d2*dy**2 &
+2.97d3*a1*a2**4*dy**2+8.91d3*a1**2*a2**3*dy**2 &
+8.91d3*a1**3*a2**2*dy**2+2.97d3*a1**4*a2*dy**2 &
+4.9d2*a1**4*a2**4*d2*dx**6-4.655d3*a1**3*a2**4*dx**6 &
-4.655d3*a1**4*a2**3*dx**6-9.8d2*a1**4*a2**4*d2**2*dx**4 &
+9.31d3*a1**3*a2**4*d2*dx**4+9.31d3*a1**4*a2**3*d2*dx**4 &
+5.92d2*a1**4*a2**4*d2**3*dx**2-5.48d3*a1**3*a2**4*d2**2*dx**2 &
-5.48d3*a1**4*a2**3*d2**2*dx**2-1.62d3*a1**2*a2**4*d2*dx**2 &
-3.24d3*a1**3*a2**3*d2*dx**2-1.62d3*a1**4*a2**2*d2*dx**2 &
+2.97d3*a1*a2**4*dx**2+8.91d3*a1**2*a2**3*dx**2 &
+8.91d3*a1**3*a2**2*dx**2+2.97d3*a1**4*a2*dx**2 &
-9.6d1*a1**4*a2**4*d2**4+7.2d2*a1**3*a2**4*d2**3 &
+7.2d2*a1**4*a2**3*d2**3+2.28d3*a1**2*a2**4*d2**2 &
+4.56d3*a1**3*a2**3*d2**2+2.28d3*a1**4*a2**2*d2**2-5.22d3*a1*a2**4*d2 &
-1.566d4*a1**2*a2**3*d2-1.566d4*a1**3*a2**2*d2-5.22d3*a1**4*a2*d2 &
+2.34d3*a2**4+9.36d3*a1*a2**3+1.404d4*a1**2*a2**2+9.36d3*a1**3*a2 &
+2.34d3*a1**4)
                case (1)
                  rlYlm_laplacian =3.52494601119984446d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dy*(1.96d2*a1**4*a2**4*d2*dx**2*dy**4 &
-1.862d3*a1**3*a2**4*dx**2*dy**4-1.862d3*a1**4*a2**3*dx**2*dy**4 &
-9.8d1*a1**3*a2**4*d2*dy**4-9.8d1*a1**4*a2**3*d2*dy**4 &
+8.33d2*a1**2*a2**4*dy**4+1.666d3*a1**3*a2**3*dy**4 &
+8.33d2*a1**4*a2**2*dy**4+3.92d2*a1**4*a2**4*d2*dx**4*dy**2 &
-3.724d3*a1**3*a2**4*dx**4*dy**2-3.724d3*a1**4*a2**3*dx**4*dy**2 &
-2.8d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+2.604d3*a1**3*a2**4*d2*dx**2*dy**2+2.604d3*a1**4*a2**3*d2*dx**2*dy**2 &
+4.76d2*a1**2*a2**4*dx**2*dy**2+9.52d2*a1**3*a2**3*dx**2*dy**2 &
+4.76d2*a1**4*a2**2*dx**2*dy**2+1.4d2*a1**3*a2**4*d2**2*dy**2 &
+1.4d2*a1**4*a2**3*d2**2*dy**2-1.316d3*a1**2*a2**4*d2*dy**2 &
-2.632d3*a1**3*a2**3*d2*dy**2-1.316d3*a1**4*a2**2*d2*dy**2 &
+9.45d2*a1*a2**4*dy**2+2.835d3*a1**2*a2**3*dy**2 &
+2.835d3*a1**3*a2**2*dy**2+9.45d2*a1**4*a2*dy**2 &
+1.96d2*a1**4*a2**4*d2*dx**6-1.862d3*a1**3*a2**4*dx**6 &
-1.862d3*a1**4*a2**3*dx**6-2.8d2*a1**4*a2**4*d2**2*dx**4 &
+2.702d3*a1**3*a2**4*d2*dx**4+2.702d3*a1**4*a2**3*d2*dx**4 &
-3.57d2*a1**2*a2**4*dx**4-7.14d2*a1**3*a2**3*dx**4 &
-3.57d2*a1**4*a2**2*dx**4+9.6d1*a1**4*a2**4*d2**3*dx**2 &
-9.0d2*a1**3*a2**4*d2**2*dx**2-9.0d2*a1**4*a2**3*d2**2*dx**2 &
-1.56d2*a1**2*a2**4*d2*dx**2-3.12d2*a1**3*a2**3*d2*dx**2 &
-1.56d2*a1**4*a2**2*d2*dx**2+4.05d2*a1*a2**4*dx**2 &
+1.215d3*a1**2*a2**3*dx**2+1.215d3*a1**3*a2**2*dx**2 &
+4.05d2*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+5.16d2*a1**2*a2**4*d2**2+1.032d3*a1**3*a2**3*d2**2 &
+5.16d2*a1**4*a2**2*d2**2-8.64d2*a1*a2**4*d2-2.592d3*a1**2*a2**3*d2 &
-2.592d3*a1**3*a2**2*d2-8.64d2*a1**4*a2*d2+3.51d2*a2**4 &
+1.404d3*a1*a2**3+2.106d3*a1**2*a2**2+1.404d3*a1**3*a2+3.51d2*a1**4 &
)*dz
                case (2)
                  rlYlm_laplacian = &
-4.98502645567176383d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(9.8d1*a1**3*a2**3*d2*dy**4 &
-9.31d2*a1**2*a2**3*dy**4-9.31d2*a1**3*a2**2*dy**4 &
+1.96d2*a1**3*a2**3*d2*dx**2*dy**2-1.862d3*a1**2*a2**3*dx**2*dy**2 &
-1.862d3*a1**3*a2**2*dx**2*dy**2-1.68d2*a1**3*a2**3*d2**2*dy**2 &
+1.68d3*a1**2*a2**3*d2*dy**2+1.68d3*a1**3*a2**2*d2*dy**2 &
-7.14d2*a1*a2**3*dy**2-1.428d3*a1**2*a2**2*dy**2-7.14d2*a1**3*a2*dy**2 &
+9.8d1*a1**3*a2**3*d2*dx**4-9.31d2*a1**2*a2**3*dx**4 &
-9.31d2*a1**3*a2**2*dx**4-1.68d2*a1**3*a2**3*d2**2*dx**2 &
+1.68d3*a1**2*a2**3*d2*dx**2+1.68d3*a1**3*a2**2*d2*dx**2 &
-7.14d2*a1*a2**3*dx**2-1.428d3*a1**2*a2**2*dx**2-7.14d2*a1**3*a2*dx**2 &
+7.2d1*a1**3*a2**3*d2**3-7.8d2*a1**2*a2**3*d2**2 &
-7.8d2*a1**3*a2**2*d2**2+8.7d2*a1*a2**3*d2+1.74d3*a1**2*a2**2*d2 &
+8.7d2*a1**3*a2*d2-4.05d2*a2**3-1.215d3*a1*a2**2-1.215d3*a1**2*a2 &
-4.05d2*a1**3)
                case (3)
                  rlYlm_laplacian =9.3261305305638875d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dy*(8.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
-7.98d2*a1**3*a2**4*dx**2*dy**4-7.98d2*a1**4*a2**3*dx**2*dy**4 &
-4.2d1*a1**3*a2**4*d2*dy**4-4.2d1*a1**4*a2**3*d2*dy**4 &
+3.57d2*a1**2*a2**4*dy**4+7.14d2*a1**3*a2**3*dy**4 &
+3.57d2*a1**4*a2**2*dy**4+5.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
-5.32d2*a1**3*a2**4*dx**4*dy**2-5.32d2*a1**4*a2**3*dx**4*dy**2 &
-7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2+6.36d2*a1**3*a2**4*d2*dx**2*dy**2 &
+6.36d2*a1**4*a2**3*d2*dx**2*dy**2+4.08d2*a1**2*a2**4*dx**2*dy**2 &
+8.16d2*a1**3*a2**3*dx**2*dy**2+4.08d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**3*a2**4*d2**2*dy**2+3.6d1*a1**4*a2**3*d2**2*dy**2 &
-2.88d2*a1**2*a2**4*d2*dy**2-5.76d2*a1**3*a2**3*d2*dy**2 &
-2.88d2*a1**4*a2**2*d2*dy**2-1.35d2*a1*a2**4*dy**2 &
-4.05d2*a1**2*a2**3*dy**2-4.05d2*a1**3*a2**2*dy**2 &
-1.35d2*a1**4*a2*dy**2-2.8d1*a1**4*a2**4*d2*dx**6 &
+2.66d2*a1**3*a2**4*dx**6+2.66d2*a1**4*a2**3*dx**6 &
+2.4d1*a1**4*a2**4*d2**2*dx**4-2.82d2*a1**3*a2**4*d2*dx**4 &
-2.82d2*a1**4*a2**3*d2*dx**4+4.59d2*a1**2*a2**4*dx**4 &
+9.18d2*a1**3*a2**3*dx**4+4.59d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**3*a2**4*d2**2*dx**2+3.6d1*a1**4*a2**3*d2**2*dx**2 &
-2.88d2*a1**2*a2**4*d2*dx**2-5.76d2*a1**3*a2**3*d2*dx**2 &
-2.88d2*a1**4*a2**2*d2*dx**2-1.35d2*a1*a2**4*dx**2 &
-4.05d2*a1**2*a2**3*dx**2-4.05d2*a1**3*a2**2*dx**2 &
-1.35d2*a1**4*a2*dx**2-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+3.0d2*a1*a2**4*d2+9.0d2*a1**2*a2**3*d2 &
+9.0d2*a1**3*a2**2*d2+3.0d2*a1**4*a2*d2-1.95d2*a2**4-7.8d2*a1*a2**3 &
-1.17d3*a1**2*a2**2-7.8d2*a1**3*a2-1.95d2*a1**4)*dz
                case (4)
                  rlYlm_laplacian = &
-6.59457014039261917d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.33d2*a1**3*a2**4*dy**6 &
-1.33d2*a1**4*a2**3*dy**6-7.0d1*a1**4*a2**4*d2*dx**2*dy**4 &
+6.65d2*a1**3*a2**4*dx**2*dy**4+6.65d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+1.62d2*a1**3*a2**4*d2*dy**4 &
+1.62d2*a1**4*a2**3*d2*dy**4-4.08d2*a1**2*a2**4*dy**4 &
-8.16d2*a1**3*a2**3*dy**4-4.08d2*a1**4*a2**2*dy**4 &
-7.0d1*a1**4*a2**4*d2*dx**4*dy**2+6.65d2*a1**3*a2**4*dx**4*dy**2 &
+6.65d2*a1**4*a2**3*dx**4*dy**2+7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-5.24d2*a1**3*a2**4*d2*dx**2*dy**2-5.24d2*a1**4*a2**3*d2*dx**2*dy**2 &
-1.36d3*a1**2*a2**4*dx**2*dy**2-2.72d3*a1**3*a2**3*dx**2*dy**2 &
-1.36d3*a1**4*a2**2*dx**2*dy**2-4.8d1*a1**3*a2**4*d2**2*dy**2 &
-4.8d1*a1**4*a2**3*d2**2*dy**2+3.0d2*a1**2*a2**4*d2*dy**2 &
+6.0d2*a1**3*a2**3*d2*dy**2+3.0d2*a1**4*a2**2*d2*dy**2 &
+8.1d2*a1*a2**4*dy**2+2.43d3*a1**2*a2**3*dy**2 &
+2.43d3*a1**3*a2**2*dy**2+8.1d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.33d2*a1**3*a2**4*dx**6 &
-1.33d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+1.62d2*a1**3*a2**4*d2*dx**4+1.62d2*a1**4*a2**3*d2*dx**4 &
-4.08d2*a1**2*a2**4*dx**4-8.16d2*a1**3*a2**3*dx**4 &
-4.08d2*a1**4*a2**2*dx**4-4.8d1*a1**3*a2**4*d2**2*dx**2 &
-4.8d1*a1**4*a2**3*d2**2*dx**2+3.0d2*a1**2*a2**4*d2*dx**2 &
+6.0d2*a1**3*a2**3*d2*dx**2+3.0d2*a1**4*a2**2*d2*dx**2 &
+8.1d2*a1*a2**4*dx**2+2.43d3*a1**2*a2**3*dx**2 &
+2.43d3*a1**3*a2**2*dx**2+8.1d2*a1**4*a2*dx**2+7.2d1*a1**2*a2**4*d2**2 &
+1.44d2*a1**3*a2**3*d2**2+7.2d1*a1**4*a2**2*d2**2-5.16d2*a1*a2**4*d2 &
-1.548d3*a1**2*a2**3*d2-1.548d3*a1**3*a2**2*d2-5.16d2*a1**4*a2*d2 &
-1.56d2*a2**4-6.24d2*a1*a2**3-9.36d2*a1**2*a2**2-6.24d2*a1**3*a2 &
-1.56d2*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (-1)
          ! selection on l2: l1=4, m1=-1
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=-1, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-2.10187170614922326d0*E*a1*a2**5*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dy*(7.0d0*(dy**2+dx**2)-4.0d0*d2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=-1, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-1.82027429302096805d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)* &
(2.8d1*a1**2*a2**2*d2*dy**4-1.82d2*a1*a2**2*dy**4 &
-1.82d2*a1**2*a2*dy**4+2.8d1*a1**2*a2**2*d2*dx**2*dy**2 &
-1.82d2*a1*a2**2*dx**2*dy**2-1.82d2*a1**2*a2*dx**2*dy**2 &
-1.6d1*a1**2*a2**2*d2**2*dy**2+7.8d1*a1*a2**2*d2*dy**2 &
+7.8d1*a1**2*a2*d2*dy**2+1.43d2*a2**2*dy**2+2.86d2*a1*a2*dy**2 &
+1.43d2*a1**2*dy**2-1.4d1*a1*a2**2*d2*dx**2-1.4d1*a1**2*a2*d2*dx**2 &
+7.7d1*a2**2*dx**2+1.54d2*a1*a2*dx**2+7.7d1*a1**2*dx**2 &
+8.0d0*a1*a2**2*d2**2+8.0d0*a1**2*a2*d2**2-4.4d1*a2**2*d2 &
-8.8d1*a1*a2*d2-4.4d1*a1**2*d2)*dz
                case (0)
                  rlYlm_laplacian =1.82027429302096805d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(2.8d1*a1**2*a2**2*d2*dy**4 &
-1.82d2*a1*a2**2*dy**4-1.82d2*a1**2*a2*dy**4 &
+5.6d1*a1**2*a2**2*d2*dx**2*dy**2-3.64d2*a1*a2**2*dx**2*dy**2 &
-3.64d2*a1**2*a2*dx**2*dy**2-4.4d1*a1**2*a2**2*d2**2*dy**2 &
+3.16d2*a1*a2**2*d2*dy**2+3.16d2*a1**2*a2*d2*dy**2-1.65d2*a2**2*dy**2 &
-3.3d2*a1*a2*dy**2-1.65d2*a1**2*dy**2+2.8d1*a1**2*a2**2*d2*dx**4 &
-1.82d2*a1*a2**2*dx**4-1.82d2*a1**2*a2*dx**4 &
-4.4d1*a1**2*a2**2*d2**2*dx**2+3.16d2*a1*a2**2*d2*dx**2 &
+3.16d2*a1**2*a2*d2*dx**2-1.65d2*a2**2*dx**2-3.3d2*a1*a2*dx**2 &
-1.65d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3-1.28d2*a1*a2**2*d2**2 &
-1.28d2*a1**2*a2*d2**2+1.32d2*a2**2*d2+2.64d2*a1*a2*d2+1.32d2*a1**2*d2 &
)
                case (1)
                  rlYlm_laplacian = &
-3.6405485860419361d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(1.4d1*a1**2*a2**2*d2*dy**2-9.1d1*a1*a2**2*dy**2-9.1d1*a1**2*a2*dy**2 &
+1.4d1*a1**2*a2**2*d2*dx**2-9.1d1*a1*a2**2*dx**2-9.1d1*a1**2*a2*dx**2 &
-8.0d0*a1**2*a2**2*d2**2+4.6d1*a1*a2**2*d2+4.6d1*a1**2*a2*d2 &
+3.3d1*a2**2+6.6d1*a1*a2+3.3d1*a1**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=-1, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-4.07025705689025559d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dx* &
(2.8d1*a1**3*a2**3*d2*dy**4-2.1d2*a1**2*a2**3*dy**4 &
-2.1d2*a1**3*a2**2*dy**4+2.8d1*a1**3*a2**3*d2*dx**2*dy**2 &
-2.1d2*a1**2*a2**3*dx**2*dy**2-2.1d2*a1**3*a2**2*dx**2*dy**2 &
-1.6d1*a1**3*a2**3*d2**2*dy**2+8.2d1*a1**2*a2**3*d2*dy**2 &
+8.2d1*a1**3*a2**2*d2*dy**2+2.47d2*a1*a2**3*dy**2 &
+4.94d2*a1**2*a2**2*dy**2+2.47d2*a1**3*a2*dy**2 &
-1.4d1*a1**2*a2**3*d2*dx**2-1.4d1*a1**3*a2**2*d2*dx**2 &
+9.1d1*a1*a2**3*dx**2+1.82d2*a1**2*a2**2*dx**2+9.1d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.6d1*a1*a2**3*d2 &
-9.2d1*a1**2*a2**2*d2-4.6d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =2.03512852844512779d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*(5.6d1*a1**3*a2**3*d2*dy**6 &
-4.2d2*a1**2*a2**3*dy**6-4.2d2*a1**3*a2**2*dy**6 &
+1.12d2*a1**3*a2**3*d2*dx**2*dy**4-8.4d2*a1**2*a2**3*dx**2*dy**4 &
-8.4d2*a1**3*a2**2*dx**2*dy**4-8.8d1*a1**3*a2**3*d2**2*dy**4 &
+6.68d2*a1**2*a2**3*d2*dy**4+6.68d2*a1**3*a2**2*d2*dy**4 &
-5.2d1*a1*a2**3*dy**4-1.04d2*a1**2*a2**2*dy**4-5.2d1*a1**3*a2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**4*dy**2-4.2d2*a1**2*a2**3*dx**4*dy**2 &
-4.2d2*a1**3*a2**2*dx**4*dy**2-8.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+6.4d2*a1**2*a2**3*d2*dx**2*dy**2+6.4d2*a1**3*a2**2*d2*dx**2*dy**2 &
+1.3d2*a1*a2**3*dx**2*dy**2+2.6d2*a1**2*a2**2*dx**2*dy**2 &
+1.3d2*a1**3*a2*dx**2*dy**2+3.2d1*a1**3*a2**3*d2**3*dy**2 &
-2.2d2*a1**2*a2**3*d2**2*dy**2-2.2d2*a1**3*a2**2*d2**2*dy**2 &
-1.72d2*a1*a2**3*d2*dy**2-3.44d2*a1**2*a2**2*d2*dy**2 &
-1.72d2*a1**3*a2*d2*dy**2+2.31d2*a2**3*dy**2+6.93d2*a1*a2**2*dy**2 &
+6.93d2*a1**2*a2*dy**2+2.31d2*a1**3*dy**2-2.8d1*a1**2*a2**3*d2*dx**4 &
-2.8d1*a1**3*a2**2*d2*dx**4+1.82d2*a1*a2**3*dx**4 &
+3.64d2*a1**2*a2**2*dx**4+1.82d2*a1**3*a2*dx**4 &
+4.4d1*a1**2*a2**3*d2**2*dx**2+4.4d1*a1**3*a2**2*d2**2*dx**2 &
-3.16d2*a1*a2**3*d2*dx**2-6.32d2*a1**2*a2**2*d2*dx**2 &
-3.16d2*a1**3*a2*d2*dx**2+1.65d2*a2**3*dx**2+4.95d2*a1*a2**2*dx**2 &
+4.95d2*a1**2*a2*dx**2+1.65d2*a1**3*dx**2-1.6d1*a1**2*a2**3*d2**3 &
-1.6d1*a1**3*a2**2*d2**3+1.28d2*a1*a2**3*d2**2 &
+2.56d2*a1**2*a2**2*d2**2+1.28d2*a1**3*a2*d2**2-1.32d2*a2**3*d2 &
-3.96d2*a1*a2**2*d2-3.96d2*a1**2*a2*d2-1.32d2*a1**3*d2)
                case (0)
                  rlYlm_laplacian =2.34996400746656297d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(4.2d1*a1**3*a2**3*d2*dy**4 &
-3.15d2*a1**2*a2**3*dy**4-3.15d2*a1**3*a2**2*dy**4 &
+8.4d1*a1**3*a2**3*d2*dx**2*dy**2-6.3d2*a1**2*a2**3*dx**2*dy**2 &
-6.3d2*a1**3*a2**2*dx**2*dy**2-5.2d1*a1**3*a2**3*d2**2*dy**2 &
+4.24d2*a1**2*a2**3*d2*dy**2+4.24d2*a1**3*a2**2*d2*dy**2 &
-2.21d2*a1*a2**3*dy**2-4.42d2*a1**2*a2**2*dy**2-2.21d2*a1**3*a2*dy**2 &
+4.2d1*a1**3*a2**3*d2*dx**4-3.15d2*a1**2*a2**3*dx**4 &
-3.15d2*a1**3*a2**2*dx**4-5.2d1*a1**3*a2**3*d2**2*dx**2 &
+4.24d2*a1**2*a2**3*d2*dx**2+4.24d2*a1**3*a2**2*d2*dx**2 &
-2.21d2*a1*a2**3*dx**2-4.42d2*a1**2*a2**2*dx**2-2.21d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.6d2*a1**2*a2**3*d2**2 &
-1.6d2*a1**3*a2**2*d2**2+2.96d2*a1*a2**3*d2+5.92d2*a1**2*a2**2*d2 &
+2.96d2*a1**3*a2*d2-1.98d2*a2**3-5.94d2*a1*a2**2-5.94d2*a1**2*a2 &
-1.98d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian =4.07025705689025559d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(2.8d1*a1**3*a2**3*d2*dy**4 &
-2.1d2*a1**2*a2**3*dy**4-2.1d2*a1**3*a2**2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**2*dy**2-4.2d2*a1**2*a2**3*dx**2*dy**2 &
-4.2d2*a1**3*a2**2*dx**2*dy**2-4.4d1*a1**3*a2**3*d2**2*dy**2 &
+3.48d2*a1**2*a2**3*d2*dy**2+3.48d2*a1**3*a2**2*d2*dy**2 &
-1.17d2*a1*a2**3*dy**2-2.34d2*a1**2*a2**2*dy**2-1.17d2*a1**3*a2*dy**2 &
+2.8d1*a1**3*a2**3*d2*dx**4-2.1d2*a1**2*a2**3*dx**4 &
-2.1d2*a1**3*a2**2*dx**4-4.4d1*a1**3*a2**3*d2**2*dx**2 &
+3.48d2*a1**2*a2**3*d2*dx**2+3.48d2*a1**3*a2**2*d2*dx**2 &
-1.17d2*a1*a2**3*dx**2-2.34d2*a1**2*a2**2*dx**2-1.17d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.32d2*a1**2*a2**3*d2**2 &
-1.32d2*a1**3*a2**2*d2**2+7.2d1*a1*a2**3*d2+1.44d2*a1**2*a2**2*d2 &
+7.2d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2+9.9d1*a1**2*a2 &
+3.3d1*a1**3)
                case (2)
                  rlYlm_laplacian =4.07025705689025559d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.05d2*a1**2*a2**3*dy**4-1.05d2*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+3.4d1*a1**2*a2**3*d2*dy**2 &
+3.4d1*a1**3*a2**2*d2*dy**2+1.69d2*a1*a2**3*dy**2 &
+3.38d2*a1**2*a2**2*dy**2+1.69d2*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.05d2*a1**2*a2**3*dx**4 &
+1.05d2*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-6.2d1*a1**2*a2**3*d2*dx**2-6.2d1*a1**3*a2**2*d2*dx**2 &
+1.3d1*a1*a2**3*dx**2+2.6d1*a1**2*a2**2*dx**2+1.3d1*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.6d1*a1*a2**3*d2 &
-9.2d1*a1**2*a2**2*d2-4.6d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=-1, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =2.19819004679753972d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*(2.8d1*a1**3*a2**3*d2*dy**6 &
-2.38d2*a1**2*a2**3*dy**6-2.38d2*a1**3*a2**2*dy**6 &
-5.6d1*a1**3*a2**3*d2*dx**2*dy**4+4.76d2*a1**2*a2**3*dx**2*dy**4 &
+4.76d2*a1**3*a2**2*dx**2*dy**4-1.6d1*a1**3*a2**3*d2**2*dy**4 &
+5.8d1*a1**2*a2**3*d2*dy**4+5.8d1*a1**3*a2**2*d2*dy**4 &
+5.85d2*a1*a2**3*dy**4+1.17d3*a1**2*a2**2*dy**4+5.85d2*a1**3*a2*dy**4 &
-8.4d1*a1**3*a2**3*d2*dx**4*dy**2+7.14d2*a1**2*a2**3*dx**4*dy**2 &
+7.14d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-3.0d2*a1**2*a2**3*d2*dx**2*dy**2-3.0d2*a1**3*a2**2*d2*dx**2*dy**2 &
-8.1d2*a1*a2**3*dx**2*dy**2-1.62d3*a1**2*a2**2*dx**2*dy**2 &
-8.1d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.44d2*a1*a2**3*d2*dy**2 &
-2.88d2*a1**2*a2**2*d2*dy**2-1.44d2*a1**3*a2*d2*dy**2 &
-2.34d2*a2**3*dy**2-7.02d2*a1*a2**2*dy**2-7.02d2*a1**2*a2*dy**2 &
-2.34d2*a1**3*dy**2+4.2d1*a1**2*a2**3*d2*dx**4 &
+4.2d1*a1**3*a2**2*d2*dx**4-3.15d2*a1*a2**3*dx**4 &
-6.3d2*a1**2*a2**2*dx**4-3.15d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.44d2*a1*a2**3*d2*dx**2+2.88d2*a1**2*a2**2*d2*dx**2 &
+1.44d2*a1**3*a2*d2*dx**2+2.34d2*a2**3*dx**2+7.02d2*a1*a2**2*dx**2 &
+7.02d2*a1**2*a2*dx**2+2.34d2*a1**3*dx**2)*dz
                case (-2)
                  rlYlm_laplacian =5.3844439723186478d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*(5.6d1*a1**4*a2**4*d2*dy**6 &
-4.76d2*a1**3*a2**4*dy**6-4.76d2*a1**4*a2**3*dy**6 &
+1.12d2*a1**4*a2**4*d2*dx**2*dy**4-9.52d2*a1**3*a2**4*dx**2*dy**4 &
-9.52d2*a1**4*a2**3*dx**2*dy**4-8.8d1*a1**4*a2**4*d2**2*dy**4 &
+7.32d2*a1**3*a2**4*d2*dy**4+7.32d2*a1**4*a2**3*d2*dy**4 &
+1.2d2*a1**2*a2**4*dy**4+2.4d2*a1**3*a2**3*dy**4 &
+1.2d2*a1**4*a2**2*dy**4+5.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
-4.76d2*a1**3*a2**4*dx**4*dy**2-4.76d2*a1**4*a2**3*dx**4*dy**2 &
-8.8d1*a1**4*a2**4*d2**2*dx**2*dy**2+7.04d2*a1**3*a2**4*d2*dx**2*dy**2 &
+7.04d2*a1**4*a2**3*d2*dx**2*dy**2+3.3d2*a1**2*a2**4*dx**2*dy**2 &
+6.6d2*a1**3*a2**3*dx**2*dy**2+3.3d2*a1**4*a2**2*dx**2*dy**2 &
+3.2d1*a1**4*a2**4*d2**3*dy**2-2.28d2*a1**3*a2**4*d2**2*dy**2 &
-2.28d2*a1**4*a2**3*d2**2*dy**2-3.72d2*a1**2*a2**4*d2*dy**2 &
-7.44d2*a1**3*a2**3*d2*dy**2-3.72d2*a1**4*a2**2*d2*dy**2 &
+2.73d2*a1*a2**4*dy**2+8.19d2*a1**2*a2**3*dy**2 &
+8.19d2*a1**3*a2**2*dy**2+2.73d2*a1**4*a2*dy**2 &
-2.8d1*a1**3*a2**4*d2*dx**4-2.8d1*a1**4*a2**3*d2*dx**4 &
+2.1d2*a1**2*a2**4*dx**4+4.2d2*a1**3*a2**3*dx**4 &
+2.1d2*a1**4*a2**2*dx**4+4.4d1*a1**3*a2**4*d2**2*dx**2 &
+4.4d1*a1**4*a2**3*d2**2*dx**2-3.48d2*a1**2*a2**4*d2*dx**2 &
-6.96d2*a1**3*a2**3*d2*dx**2-3.48d2*a1**4*a2**2*d2*dx**2 &
+1.17d2*a1*a2**4*dx**2+3.51d2*a1**2*a2**3*dx**2 &
+3.51d2*a1**3*a2**2*dx**2+1.17d2*a1**4*a2*dx**2 &
-1.6d1*a1**3*a2**4*d2**3-1.6d1*a1**4*a2**3*d2**3 &
+1.32d2*a1**2*a2**4*d2**2+2.64d2*a1**3*a2**3*d2**2 &
+1.32d2*a1**4*a2**2*d2**2-7.2d1*a1*a2**4*d2-2.16d2*a1**2*a2**3*d2 &
-2.16d2*a1**3*a2**2*d2-7.2d1*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (-1)
                  rlYlm_laplacian =1.70271068860915474d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(1.4d2*a1**4*a2**4*d2*dy**6 &
-1.19d3*a1**3*a2**4*dy**6-1.19d3*a1**4*a2**3*dy**6 &
+2.8d2*a1**4*a2**4*d2*dx**2*dy**4-2.38d3*a1**3*a2**4*dx**2*dy**4 &
-2.38d3*a1**4*a2**3*dx**2*dy**4-1.92d2*a1**4*a2**4*d2**2*dy**4 &
+1.69d3*a1**3*a2**4*d2*dy**4+1.69d3*a1**4*a2**3*d2*dy**4 &
-4.35d2*a1**2*a2**4*dy**4-8.7d2*a1**3*a2**3*dy**4 &
-4.35d2*a1**4*a2**2*dy**4+1.4d2*a1**4*a2**4*d2*dx**4*dy**2 &
-1.19d3*a1**3*a2**4*dx**4*dy**2-1.19d3*a1**4*a2**3*dx**4*dy**2 &
-1.92d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.62d3*a1**3*a2**4*d2*dx**2*dy**2+1.62d3*a1**4*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**4*dx**2*dy**2+1.8d2*a1**3*a2**3*dx**2*dy**2 &
+9.0d1*a1**4*a2**2*dx**2*dy**2+6.4d1*a1**4*a2**4*d2**3*dy**2 &
-5.76d2*a1**3*a2**4*d2**2*dy**2-5.76d2*a1**4*a2**3*d2**2*dy**2 &
+2.28d2*a1**2*a2**4*d2*dy**2+4.56d2*a1**3*a2**3*d2*dy**2 &
+2.28d2*a1**4*a2**2*d2*dy**2+7.8d1*a1*a2**4*dy**2 &
+2.34d2*a1**2*a2**3*dy**2+2.34d2*a1**3*a2**2*dy**2 &
+7.8d1*a1**4*a2*dy**2-7.0d1*a1**3*a2**4*d2*dx**4 &
-7.0d1*a1**4*a2**3*d2*dx**4+5.25d2*a1**2*a2**4*dx**4 &
+1.05d3*a1**3*a2**3*dx**4+5.25d2*a1**4*a2**2*dx**4 &
+9.6d1*a1**3*a2**4*d2**2*dx**2+9.6d1*a1**4*a2**3*d2**2*dx**2 &
-8.28d2*a1**2*a2**4*d2*dx**2-1.656d3*a1**3*a2**3*d2*dx**2 &
-8.28d2*a1**4*a2**2*d2*dx**2+7.02d2*a1*a2**4*dx**2 &
+2.106d3*a1**2*a2**3*dx**2+2.106d3*a1**3*a2**2*dx**2 &
+7.02d2*a1**4*a2*dx**2-3.2d1*a1**3*a2**4*d2**3-3.2d1*a1**4*a2**3*d2**3 &
+3.36d2*a1**2*a2**4*d2**2+6.72d2*a1**3*a2**3*d2**2 &
+3.36d2*a1**4*a2**2*d2**2-6.84d2*a1*a2**4*d2-2.052d3*a1**2*a2**3*d2 &
-2.052d3*a1**3*a2**2*d2-6.84d2*a1**4*a2*d2+3.3d2*a2**4+1.32d3*a1*a2**3 &
+1.98d3*a1**2*a2**2+1.32d3*a1**3*a2+3.3d2*a1**4)*dz
                case (0)
                  rlYlm_laplacian = &
-1.39025745555846884d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(1.4d2*a1**4*a2**4*d2*dy**6-1.19d3*a1**3*a2**4*dy**6 &
-1.19d3*a1**4*a2**3*dy**6+4.2d2*a1**4*a2**4*d2*dx**2*dy**4 &
-3.57d3*a1**3*a2**4*dx**2*dy**4-3.57d3*a1**4*a2**3*dx**2*dy**4 &
-2.76d2*a1**4*a2**4*d2**2*dy**4+2.46d3*a1**3*a2**4*d2*dy**4 &
+2.46d3*a1**4*a2**3*d2*dy**4-8.55d2*a1**2*a2**4*dy**4 &
-1.71d3*a1**3*a2**3*dy**4-8.55d2*a1**4*a2**2*dy**4 &
+4.2d2*a1**4*a2**4*d2*dx**4*dy**2-3.57d3*a1**3*a2**4*dx**4*dy**2 &
-3.57d3*a1**4*a2**3*dx**4*dy**2-5.52d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+4.92d3*a1**3*a2**4*d2*dx**2*dy**2+4.92d3*a1**4*a2**3*d2*dx**2*dy**2 &
-1.71d3*a1**2*a2**4*dx**2*dy**2-3.42d3*a1**3*a2**3*dx**2*dy**2 &
-1.71d3*a1**4*a2**2*dx**2*dy**2+1.68d2*a1**4*a2**4*d2**3*dy**2 &
-1.62d3*a1**3*a2**4*d2**2*dy**2-1.62d3*a1**4*a2**3*d2**2*dy**2 &
+1.53d3*a1**2*a2**4*d2*dy**2+3.06d3*a1**3*a2**3*d2*dy**2 &
+1.53d3*a1**4*a2**2*d2*dy**2-5.85d2*a1*a2**4*dy**2 &
-1.755d3*a1**2*a2**3*dy**2-1.755d3*a1**3*a2**2*dy**2 &
-5.85d2*a1**4*a2*dy**2+1.4d2*a1**4*a2**4*d2*dx**6 &
-1.19d3*a1**3*a2**4*dx**6-1.19d3*a1**4*a2**3*dx**6 &
-2.76d2*a1**4*a2**4*d2**2*dx**4+2.46d3*a1**3*a2**4*d2*dx**4 &
+2.46d3*a1**4*a2**3*d2*dx**4-8.55d2*a1**2*a2**4*dx**4 &
-1.71d3*a1**3*a2**3*dx**4-8.55d2*a1**4*a2**2*dx**4 &
+1.68d2*a1**4*a2**4*d2**3*dx**2-1.62d3*a1**3*a2**4*d2**2*dx**2 &
-1.62d3*a1**4*a2**3*d2**2*dx**2+1.53d3*a1**2*a2**4*d2*dx**2 &
+3.06d3*a1**3*a2**3*d2*dx**2+1.53d3*a1**4*a2**2*d2*dx**2 &
-5.85d2*a1*a2**4*dx**2-1.755d3*a1**2*a2**3*dx**2 &
-1.755d3*a1**3*a2**2*dx**2-5.85d2*a1**4*a2*dx**2 &
-3.2d1*a1**4*a2**4*d2**4+3.68d2*a1**3*a2**4*d2**3 &
+3.68d2*a1**4*a2**3*d2**3-8.64d2*a1**2*a2**4*d2**2 &
-1.728d3*a1**3*a2**3*d2**2-8.64d2*a1**4*a2**2*d2**2+9.96d2*a1*a2**4*d2 &
+2.988d3*a1**2*a2**3*d2+2.988d3*a1**3*a2**2*d2+9.96d2*a1**4*a2*d2 &
-3.3d2*a2**4-1.32d3*a1*a2**3-1.98d3*a1**2*a2**2-1.32d3*a1**3*a2 &
-3.3d2*a1**4)
                case (1)
                  rlYlm_laplacian =3.40542137721830949d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.95d2*a1**2*a2**3*dy**4-5.95d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.19d3*a1**2*a2**3*dx**2*dy**2 &
-1.19d3*a1**3*a2**2*dx**2*dy**2-9.6d1*a1**3*a2**3*d2**2*dy**2 &
+8.8d2*a1**2*a2**3*d2*dy**2+8.8d2*a1**3*a2**2*d2*dy**2 &
-4.8d2*a1*a2**3*dy**2-9.6d2*a1**2*a2**2*dy**2-4.8d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-9.6d1*a1**3*a2**3*d2**2*dx**2 &
+8.8d2*a1**2*a2**3*d2*dx**2+8.8d2*a1**3*a2**2*d2*dx**2 &
-4.8d2*a1*a2**3*dx**2-9.6d2*a1**2*a2**2*dx**2-4.8d2*a1**3*a2*dx**2 &
+3.2d1*a1**3*a2**3*d2**3-3.36d2*a1**2*a2**3*d2**2 &
-3.36d2*a1**3*a2**2*d2**2+5.28d2*a1*a2**3*d2+1.056d3*a1**2*a2**2*d2 &
+5.28d2*a1**3*a2*d2-3.12d2*a2**3-9.36d2*a1*a2**2-9.36d2*a1**2*a2 &
-3.12d2*a1**3)*dz
                case (2)
                  rlYlm_laplacian = &
-5.3844439723186478d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(2.8d1*a1**4*a2**4*d2*dy**6-2.38d2*a1**3*a2**4*dy**6 &
-2.38d2*a1**4*a2**3*dy**6+2.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
-2.38d2*a1**3*a2**4*dx**2*dy**4-2.38d2*a1**4*a2**3*dx**2*dy**4 &
-4.4d1*a1**4*a2**4*d2**2*dy**4+3.52d2*a1**3*a2**4*d2*dy**4 &
+3.52d2*a1**4*a2**3*d2*dy**4+1.65d2*a1**2*a2**4*dy**4 &
+3.3d2*a1**3*a2**3*dy**4+1.65d2*a1**4*a2**2*dy**4 &
-2.8d1*a1**4*a2**4*d2*dx**4*dy**2+2.38d2*a1**3*a2**4*dx**4*dy**2 &
+2.38d2*a1**4*a2**3*dx**4*dy**2-5.6d1*a1**3*a2**4*d2*dx**2*dy**2 &
-5.6d1*a1**4*a2**3*d2*dx**2*dy**2+4.2d2*a1**2*a2**4*dx**2*dy**2 &
+8.4d2*a1**3*a2**3*dx**2*dy**2+4.2d2*a1**4*a2**2*dx**2*dy**2 &
+1.6d1*a1**4*a2**4*d2**3*dy**2-9.2d1*a1**3*a2**4*d2**2*dy**2 &
-9.2d1*a1**4*a2**3*d2**2*dy**2-3.6d2*a1**2*a2**4*d2*dy**2 &
-7.2d2*a1**3*a2**3*d2*dy**2-3.6d2*a1**4*a2**2*d2*dy**2 &
+1.95d2*a1*a2**4*dy**2+5.85d2*a1**2*a2**3*dy**2 &
+5.85d2*a1**3*a2**2*dy**2+1.95d2*a1**4*a2*dy**2 &
-2.8d1*a1**4*a2**4*d2*dx**6+2.38d2*a1**3*a2**4*dx**6 &
+2.38d2*a1**4*a2**3*dx**6+4.4d1*a1**4*a2**4*d2**2*dx**4 &
-4.08d2*a1**3*a2**4*d2*dx**4-4.08d2*a1**4*a2**3*d2*dx**4 &
+2.55d2*a1**2*a2**4*dx**4+5.1d2*a1**3*a2**3*dx**4 &
+2.55d2*a1**4*a2**2*dx**4-1.6d1*a1**4*a2**4*d2**3*dx**2 &
+1.8d2*a1**3*a2**4*d2**2*dx**2+1.8d2*a1**4*a2**3*d2**2*dx**2 &
-3.36d2*a1**2*a2**4*d2*dx**2-6.72d2*a1**3*a2**3*d2*dx**2 &
-3.36d2*a1**4*a2**2*d2*dx**2+3.9d1*a1*a2**4*dx**2 &
+1.17d2*a1**2*a2**3*dx**2+1.17d2*a1**3*a2**2*dx**2 &
+3.9d1*a1**4*a2*dx**2-1.6d1*a1**3*a2**4*d2**3-1.6d1*a1**4*a2**3*d2**3 &
+1.32d2*a1**2*a2**4*d2**2+2.64d2*a1**3*a2**3*d2**2 &
+1.32d2*a1**4*a2**2*d2**2-7.2d1*a1*a2**4*d2-2.16d2*a1**2*a2**3*d2 &
-2.16d2*a1**3*a2**2*d2-7.2d1*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (3)
                  rlYlm_laplacian =4.39638009359507944d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(4.2d1*a1**3*a2**3*d2*dy**4 &
-3.57d2*a1**2*a2**3*dy**4-3.57d2*a1**3*a2**2*dy**4 &
+2.8d1*a1**3*a2**3*d2*dx**2*dy**2-2.38d2*a1**2*a2**3*dx**2*dy**2 &
-2.38d2*a1**3*a2**2*dx**2*dy**2-2.4d1*a1**3*a2**3*d2**2*dy**2 &
+1.08d2*a1**2*a2**3*d2*dy**2+1.08d2*a1**3*a2**2*d2*dy**2 &
+7.2d2*a1*a2**3*dy**2+1.44d3*a1**2*a2**2*dy**2+7.2d2*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.19d2*a1**2*a2**3*dx**4 &
+1.19d2*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-9.2d1*a1**2*a2**3*d2*dx**2-9.2d1*a1**3*a2**2*d2*dx**2 &
+1.8d2*a1*a2**3*dx**2+3.6d2*a1**2*a2**2*dx**2+1.8d2*a1**3*a2*dx**2 &
+2.4d1*a1**2*a2**3*d2**2+2.4d1*a1**3*a2**2*d2**2-1.44d2*a1*a2**3*d2 &
-2.88d2*a1**2*a2**2*d2-1.44d2*a1**3*a2*d2-2.34d2*a2**3-7.02d2*a1*a2**2 &
-7.02d2*a1**2*a2-2.34d2*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=-1, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =9.3261305305638875d0*E*a1**3*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-12)*dx*(2.8d1*a1**3*a2**3*d2*dy**6 &
-2.66d2*a1**2*a2**3*dy**6-2.66d2*a1**3*a2**2*dy**6 &
-1.6d1*a1**3*a2**3*d2**2*dy**4+6.2d1*a1**2*a2**3*d2*dy**4 &
+6.2d1*a1**3*a2**2*d2*dy**4+7.65d2*a1*a2**3*dy**4 &
+1.53d3*a1**2*a2**2*dy**4+7.65d2*a1**3*a2*dy**4 &
-2.8d1*a1**3*a2**3*d2*dx**4*dy**2+2.66d2*a1**2*a2**3*dx**4*dy**2 &
+2.66d2*a1**3*a2**2*dx**4*dy**2+1.6d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-1.32d2*a1**2*a2**3*d2*dx**2*dy**2-1.32d2*a1**3*a2**2*d2*dx**2*dy**2 &
-1.7d2*a1*a2**3*dx**2*dy**2-3.4d2*a1**2*a2**2*dx**2*dy**2 &
-1.7d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.5d2*a1*a2**3*d2*dy**2 &
-3.0d2*a1**2*a2**2*d2*dy**2-1.5d2*a1**3*a2*d2*dy**2-4.05d2*a2**3*dy**2 &
-1.215d3*a1*a2**2*dy**2-1.215d3*a1**2*a2*dy**2-4.05d2*a1**3*dy**2 &
+1.4d1*a1**2*a2**3*d2*dx**4+1.4d1*a1**3*a2**2*d2*dx**4 &
-1.19d2*a1*a2**3*dx**4-2.38d2*a1**2*a2**2*dx**4-1.19d2*a1**3*a2*dx**4 &
-8.0d0*a1**2*a2**3*d2**2*dx**2-8.0d0*a1**3*a2**2*d2**2*dx**2 &
+5.0d1*a1*a2**3*d2*dx**2+1.0d2*a1**2*a2**2*d2*dx**2 &
+5.0d1*a1**3*a2*d2*dx**2+1.35d2*a2**3*dx**2+4.05d2*a1*a2**2*dx**2 &
+4.05d2*a1**2*a2*dx**2+1.35d2*a1**3*dx**2)*dz
                case (-3)
                  rlYlm_laplacian = &
-3.29728507019630958d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)* &
(5.6d1*a1**4*a2**4*d2*dy**8-5.32d2*a1**3*a2**4*dy**8 &
-5.32d2*a1**4*a2**3*dy**8-5.6d1*a1**4*a2**4*d2*dx**2*dy**6 &
+5.32d2*a1**3*a2**4*dx**2*dy**6+5.32d2*a1**4*a2**3*dx**2*dy**6 &
-8.8d1*a1**4*a2**4*d2**2*dy**6+7.4d2*a1**3*a2**4*d2*dy**6 &
+7.4d2*a1**4*a2**3*d2*dy**6+8.16d2*a1**2*a2**4*dy**6 &
+1.632d3*a1**3*a2**3*dy**6+8.16d2*a1**4*a2**2*dy**6 &
-2.8d2*a1**4*a2**4*d2*dx**4*dy**4+2.66d3*a1**3*a2**4*dx**4*dy**4 &
+2.66d3*a1**4*a2**3*dx**4*dy**4+1.76d2*a1**4*a2**4*d2**2*dx**2*dy**4 &
-1.732d3*a1**3*a2**4*d2*dx**2*dy**4-1.732d3*a1**4*a2**3*d2*dx**2*dy**4 &
+5.1d2*a1**2*a2**4*dx**2*dy**4+1.02d3*a1**3*a2**3*dx**2*dy**4 &
+5.1d2*a1**4*a2**2*dx**2*dy**4+3.2d1*a1**4*a2**4*d2**3*dy**4 &
-1.48d2*a1**3*a2**4*d2**2*dy**4-1.48d2*a1**4*a2**3*d2**2*dy**4 &
-1.38d3*a1**2*a2**4*d2*dy**4-2.76d3*a1**3*a2**3*d2*dy**4 &
-1.38d3*a1**4*a2**2*d2*dy**4+4.05d2*a1*a2**4*dy**4 &
+1.215d3*a1**2*a2**3*dy**4+1.215d3*a1**3*a2**2*dy**4 &
+4.05d2*a1**4*a2*dy**4-1.68d2*a1**4*a2**4*d2*dx**6*dy**2 &
+1.596d3*a1**3*a2**4*dx**6*dy**2+1.596d3*a1**4*a2**3*dx**6*dy**2 &
+2.64d2*a1**4*a2**4*d2**2*dx**4*dy**2 &
-2.388d3*a1**3*a2**4*d2*dx**4*dy**2-2.388d3*a1**4*a2**3*d2*dx**4*dy**2 &
-1.02d3*a1**2*a2**4*dx**4*dy**2-2.04d3*a1**3*a2**3*dx**4*dy**2 &
-1.02d3*a1**4*a2**2*dx**4*dy**2-9.6d1*a1**4*a2**4*d2**3*dx**2*dy**2 &
+8.4d2*a1**3*a2**4*d2**2*dx**2*dy**2 &
+8.4d2*a1**4*a2**3*d2**2*dx**2*dy**2+7.2d2*a1**2*a2**4*d2*dx**2*dy**2 &
+1.44d3*a1**3*a2**3*d2*dx**2*dy**2+7.2d2*a1**4*a2**2*d2*dx**2*dy**2 &
-8.1d2*a1*a2**4*dx**2*dy**2-2.43d3*a1**2*a2**3*dx**2*dy**2 &
-2.43d3*a1**3*a2**2*dx**2*dy**2-8.1d2*a1**4*a2*dx**2*dy**2 &
-4.8d1*a1**3*a2**4*d2**3*dy**2-4.8d1*a1**4*a2**3*d2**3*dy**2 &
+4.08d2*a1**2*a2**4*d2**2*dy**2+8.16d2*a1**3*a2**3*d2**2*dy**2 &
+4.08d2*a1**4*a2**2*d2**2*dy**2+3.6d1*a1*a2**4*d2*dy**2 &
+1.08d2*a1**2*a2**3*d2*dy**2+1.08d2*a1**3*a2**2*d2*dy**2 &
+3.6d1*a1**4*a2*d2*dy**2-2.34d2*a2**4*dy**2-9.36d2*a1*a2**3*dy**2 &
-1.404d3*a1**2*a2**2*dy**2-9.36d2*a1**3*a2*dy**2-2.34d2*a1**4*dy**2 &
+8.4d1*a1**3*a2**4*d2*dx**6+8.4d1*a1**4*a2**3*d2*dx**6 &
-7.14d2*a1**2*a2**4*dx**6-1.428d3*a1**3*a2**3*dx**6 &
-7.14d2*a1**4*a2**2*dx**6-1.32d2*a1**3*a2**4*d2**2*dx**4 &
-1.32d2*a1**4*a2**3*d2**2*dx**4+1.14d3*a1**2*a2**4*d2*dx**4 &
+2.28d3*a1**3*a2**3*d2*dx**4+1.14d3*a1**4*a2**2*d2*dx**4 &
-1.35d2*a1*a2**4*dx**4-4.05d2*a1**2*a2**3*dx**4 &
-4.05d2*a1**3*a2**2*dx**4-1.35d2*a1**4*a2*dx**4 &
+4.8d1*a1**3*a2**4*d2**3*dx**2+4.8d1*a1**4*a2**3*d2**3*dx**2 &
-4.08d2*a1**2*a2**4*d2**2*dx**2-8.16d2*a1**3*a2**3*d2**2*dx**2 &
-4.08d2*a1**4*a2**2*d2**2*dx**2-3.6d1*a1*a2**4*d2*dx**2 &
-1.08d2*a1**2*a2**3*d2*dx**2-1.08d2*a1**3*a2**2*d2*dx**2 &
-3.6d1*a1**4*a2*d2*dx**2+2.34d2*a2**4*dx**2+9.36d2*a1*a2**3*dx**2 &
+1.404d3*a1**2*a2**2*dx**2+9.36d2*a1**3*a2*dx**2+2.34d2*a1**4*dx**2)
                case (-2)
                  rlYlm_laplacian =3.52494601119984446d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dx*(1.96d2*a1**4*a2**4*d2*dy**6 &
-1.862d3*a1**3*a2**4*dy**6-1.862d3*a1**4*a2**3*dy**6 &
+3.92d2*a1**4*a2**4*d2*dx**2*dy**4-3.724d3*a1**3*a2**4*dx**2*dy**4 &
-3.724d3*a1**4*a2**3*dx**2*dy**4-2.8d2*a1**4*a2**4*d2**2*dy**4 &
+2.702d3*a1**3*a2**4*d2*dy**4+2.702d3*a1**4*a2**3*d2*dy**4 &
-3.57d2*a1**2*a2**4*dy**4-7.14d2*a1**3*a2**3*dy**4 &
-3.57d2*a1**4*a2**2*dy**4+1.96d2*a1**4*a2**4*d2*dx**4*dy**2 &
-1.862d3*a1**3*a2**4*dx**4*dy**2-1.862d3*a1**4*a2**3*dx**4*dy**2 &
-2.8d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+2.604d3*a1**3*a2**4*d2*dx**2*dy**2+2.604d3*a1**4*a2**3*d2*dx**2*dy**2 &
+4.76d2*a1**2*a2**4*dx**2*dy**2+9.52d2*a1**3*a2**3*dx**2*dy**2 &
+4.76d2*a1**4*a2**2*dx**2*dy**2+9.6d1*a1**4*a2**4*d2**3*dy**2 &
-9.0d2*a1**3*a2**4*d2**2*dy**2-9.0d2*a1**4*a2**3*d2**2*dy**2 &
-1.56d2*a1**2*a2**4*d2*dy**2-3.12d2*a1**3*a2**3*d2*dy**2 &
-1.56d2*a1**4*a2**2*d2*dy**2+4.05d2*a1*a2**4*dy**2 &
+1.215d3*a1**2*a2**3*dy**2+1.215d3*a1**3*a2**2*dy**2 &
+4.05d2*a1**4*a2*dy**2-9.8d1*a1**3*a2**4*d2*dx**4 &
-9.8d1*a1**4*a2**3*d2*dx**4+8.33d2*a1**2*a2**4*dx**4 &
+1.666d3*a1**3*a2**3*dx**4+8.33d2*a1**4*a2**2*dx**4 &
+1.4d2*a1**3*a2**4*d2**2*dx**2+1.4d2*a1**4*a2**3*d2**2*dx**2 &
-1.316d3*a1**2*a2**4*d2*dx**2-2.632d3*a1**3*a2**3*d2*dx**2 &
-1.316d3*a1**4*a2**2*d2*dx**2+9.45d2*a1*a2**4*dx**2 &
+2.835d3*a1**2*a2**3*dx**2+2.835d3*a1**3*a2**2*dx**2 &
+9.45d2*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+5.16d2*a1**2*a2**4*d2**2+1.032d3*a1**3*a2**3*d2**2 &
+5.16d2*a1**4*a2**2*d2**2-8.64d2*a1*a2**4*d2-2.592d3*a1**2*a2**3*d2 &
-2.592d3*a1**3*a2**2*d2-8.64d2*a1**4*a2*d2+3.51d2*a2**4 &
+1.404d3*a1*a2**3+2.106d3*a1**2*a2**2+1.404d3*a1**3*a2+3.51d2*a1**4 &
)*dz
                case (-1)
                  rlYlm_laplacian =-1.24625661391794096d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(3.92d2*a1**5*a2**5*d2*dy**8 &
-3.724d3*a1**4*a2**5*dy**8-3.724d3*a1**5*a2**4*dy**8 &
+1.176d3*a1**5*a2**5*d2*dx**2*dy**6-1.1172d4*a1**4*a2**5*dx**2*dy**6 &
-1.1172d4*a1**5*a2**4*dx**2*dy**6-8.4d2*a1**5*a2**5*d2**2*dy**6 &
+8.204d3*a1**4*a2**5*d2*dy**6+8.204d3*a1**5*a2**4*d2*dy**6 &
-1.904d3*a1**3*a2**5*dy**6-3.808d3*a1**4*a2**4*dy**6 &
-1.904d3*a1**5*a2**3*dy**6+1.176d3*a1**5*a2**5*d2*dx**4*dy**4 &
-1.1172d4*a1**4*a2**5*dx**4*dy**4-1.1172d4*a1**5*a2**4*dx**4*dy**4 &
-1.68d3*a1**5*a2**5*d2**2*dx**2*dy**4 &
+1.6212d4*a1**4*a2**5*d2*dx**2*dy**4 &
+1.6212d4*a1**5*a2**4*d2*dx**2*dy**4-2.142d3*a1**3*a2**5*dx**2*dy**4 &
-4.284d3*a1**4*a2**4*dx**2*dy**4-2.142d3*a1**5*a2**3*dx**2*dy**4 &
+5.76d2*a1**5*a2**5*d2**3*dy**4-5.82d3*a1**4*a2**5*d2**2*dy**4 &
-5.82d3*a1**5*a2**4*d2**2*dy**4+3.012d3*a1**3*a2**5*d2*dy**4 &
+6.024d3*a1**4*a2**4*d2*dy**4+3.012d3*a1**5*a2**3*d2*dy**4 &
-4.05d2*a1**2*a2**5*dy**4-1.215d3*a1**3*a2**4*dy**4 &
-1.215d3*a1**4*a2**3*dy**4-4.05d2*a1**5*a2**2*dy**4 &
+3.92d2*a1**5*a2**5*d2*dx**6*dy**2-3.724d3*a1**4*a2**5*dx**6*dy**2 &
-3.724d3*a1**5*a2**4*dx**6*dy**2-8.4d2*a1**5*a2**5*d2**2*dx**4*dy**2 &
+7.812d3*a1**4*a2**5*d2*dx**4*dy**2+7.812d3*a1**5*a2**4*d2*dx**4*dy**2 &
+1.428d3*a1**3*a2**5*dx**4*dy**2+2.856d3*a1**4*a2**4*dx**4*dy**2 &
+1.428d3*a1**5*a2**3*dx**4*dy**2+5.76d2*a1**5*a2**5*d2**3*dx**2*dy**2 &
-5.4d3*a1**4*a2**5*d2**2*dx**2*dy**2 &
-5.4d3*a1**5*a2**4*d2**2*dx**2*dy**2-9.36d2*a1**3*a2**5*d2*dx**2*dy**2 &
-1.872d3*a1**4*a2**4*d2*dx**2*dy**2-9.36d2*a1**5*a2**3*d2*dx**2*dy**2 &
+2.43d3*a1**2*a2**5*dx**2*dy**2+7.29d3*a1**3*a2**4*dx**2*dy**2 &
+7.29d3*a1**4*a2**3*dx**2*dy**2+2.43d3*a1**5*a2**2*dx**2*dy**2 &
-1.28d2*a1**5*a2**5*d2**4*dy**2+1.312d3*a1**4*a2**5*d2**3*dy**2 &
+1.312d3*a1**5*a2**4*d2**3*dy**2-7.44d2*a1**3*a2**5*d2**2*dy**2 &
-1.488d3*a1**4*a2**4*d2**2*dy**2-7.44d2*a1**5*a2**3*d2**2*dy**2 &
-6.24d2*a1**2*a2**5*d2*dy**2-1.872d3*a1**3*a2**4*d2*dy**2 &
-1.872d3*a1**4*a2**3*d2*dy**2-6.24d2*a1**5*a2**2*d2*dy**2 &
+5.46d2*a1*a2**5*dy**2+2.184d3*a1**2*a2**4*dy**2 &
+3.276d3*a1**3*a2**3*dy**2+2.184d3*a1**4*a2**2*dy**2 &
+5.46d2*a1**5*a2*dy**2-1.96d2*a1**4*a2**5*d2*dx**6 &
-1.96d2*a1**5*a2**4*d2*dx**6+1.666d3*a1**3*a2**5*dx**6 &
+3.332d3*a1**4*a2**4*dx**6+1.666d3*a1**5*a2**3*dx**6 &
+4.2d2*a1**4*a2**5*d2**2*dx**4+4.2d2*a1**5*a2**4*d2**2*dx**4 &
-3.948d3*a1**3*a2**5*d2*dx**4-7.896d3*a1**4*a2**4*d2*dx**4 &
-3.948d3*a1**5*a2**3*d2*dx**4+2.835d3*a1**2*a2**5*dx**4 &
+8.505d3*a1**3*a2**4*dx**4+8.505d3*a1**4*a2**3*dx**4 &
+2.835d3*a1**5*a2**2*dx**4-2.88d2*a1**4*a2**5*d2**3*dx**2 &
-2.88d2*a1**5*a2**4*d2**3*dx**2+3.096d3*a1**3*a2**5*d2**2*dx**2 &
+6.192d3*a1**4*a2**4*d2**2*dx**2+3.096d3*a1**5*a2**3*d2**2*dx**2 &
-5.184d3*a1**2*a2**5*d2*dx**2-1.5552d4*a1**3*a2**4*d2*dx**2 &
-1.5552d4*a1**4*a2**3*d2*dx**2-5.184d3*a1**5*a2**2*d2*dx**2 &
+2.106d3*a1*a2**5*dx**2+8.424d3*a1**2*a2**4*dx**2 &
+1.2636d4*a1**3*a2**3*dx**2+8.424d3*a1**4*a2**2*dx**2 &
+2.106d3*a1**5*a2*dx**2+6.4d1*a1**4*a2**5*d2**4 &
+6.4d1*a1**5*a2**4*d2**4-8.32d2*a1**3*a2**5*d2**3 &
-1.664d3*a1**4*a2**4*d2**3-8.32d2*a1**5*a2**3*d2**3 &
+2.52d3*a1**2*a2**5*d2**2+7.56d3*a1**3*a2**4*d2**2 &
+7.56d3*a1**4*a2**3*d2**2+2.52d3*a1**5*a2**2*d2**2-2.424d3*a1*a2**5*d2 &
-9.696d3*a1**2*a2**4*d2-1.4544d4*a1**3*a2**3*d2-9.696d3*a1**4*a2**2*d2 &
-2.424d3*a1**5*a2*d2+4.62d2*a2**5+2.31d3*a1*a2**4+4.62d3*a1**2*a2**3 &
+4.62d3*a1**3*a2**2+2.31d3*a1**4*a2+4.62d2*a1**5)
                case (0)
                  rlYlm_laplacian =-7.88201889805958723d-1*E*a1**2*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-12)*dy* &
(4.9d2*a1**4*a2**4*d2*dy**6-4.655d3*a1**3*a2**4*dy**6 &
-4.655d3*a1**4*a2**3*dy**6+1.47d3*a1**4*a2**4*d2*dx**2*dy**4 &
-1.3965d4*a1**3*a2**4*dx**2*dy**4-1.3965d4*a1**4*a2**3*dx**2*dy**4 &
-8.4d2*a1**4*a2**4*d2**2*dy**4+8.4d3*a1**3*a2**4*d2*dy**4 &
+8.4d3*a1**4*a2**3*d2*dy**4-3.57d3*a1**2*a2**4*dy**4 &
-7.14d3*a1**3*a2**3*dy**4-3.57d3*a1**4*a2**2*dy**4 &
+1.47d3*a1**4*a2**4*d2*dx**4*dy**2-1.3965d4*a1**3*a2**4*dx**4*dy**2 &
-1.3965d4*a1**4*a2**3*dx**4*dy**2-1.68d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.68d4*a1**3*a2**4*d2*dx**2*dy**2+1.68d4*a1**4*a2**3*d2*dx**2*dy**2 &
-7.14d3*a1**2*a2**4*dx**2*dy**2-1.428d4*a1**3*a2**3*dx**2*dy**2 &
-7.14d3*a1**4*a2**2*dx**2*dy**2+4.32d2*a1**4*a2**4*d2**3*dy**2 &
-4.68d3*a1**3*a2**4*d2**2*dy**2-4.68d3*a1**4*a2**3*d2**2*dy**2 &
+5.22d3*a1**2*a2**4*d2*dy**2+1.044d4*a1**3*a2**3*d2*dy**2 &
+5.22d3*a1**4*a2**2*d2*dy**2-2.43d3*a1*a2**4*dy**2 &
-7.29d3*a1**2*a2**3*dy**2-7.29d3*a1**3*a2**2*dy**2 &
-2.43d3*a1**4*a2*dy**2+4.9d2*a1**4*a2**4*d2*dx**6 &
-4.655d3*a1**3*a2**4*dx**6-4.655d3*a1**4*a2**3*dx**6 &
-8.4d2*a1**4*a2**4*d2**2*dx**4+8.4d3*a1**3*a2**4*d2*dx**4 &
+8.4d3*a1**4*a2**3*d2*dx**4-3.57d3*a1**2*a2**4*dx**4 &
-7.14d3*a1**3*a2**3*dx**4-3.57d3*a1**4*a2**2*dx**4 &
+4.32d2*a1**4*a2**4*d2**3*dx**2-4.68d3*a1**3*a2**4*d2**2*dx**2 &
-4.68d3*a1**4*a2**3*d2**2*dx**2+5.22d3*a1**2*a2**4*d2*dx**2 &
+1.044d4*a1**3*a2**3*d2*dx**2+5.22d3*a1**4*a2**2*d2*dx**2 &
-2.43d3*a1*a2**4*dx**2-7.29d3*a1**2*a2**3*dx**2 &
-7.29d3*a1**3*a2**2*dx**2-2.43d3*a1**4*a2*dx**2 &
-6.4d1*a1**4*a2**4*d2**4+8.0d2*a1**3*a2**4*d2**3 &
+8.0d2*a1**4*a2**3*d2**3-1.92d3*a1**2*a2**4*d2**2 &
-3.84d3*a1**3*a2**3*d2**2-1.92d3*a1**4*a2**2*d2**2+2.28d3*a1*a2**4*d2 &
+6.84d3*a1**2*a2**3*d2+6.84d3*a1**3*a2**2*d2+2.28d3*a1**4*a2*d2 &
-7.8d2*a2**4-3.12d3*a1*a2**3-4.68d3*a1**2*a2**2-3.12d3*a1**3*a2 &
-7.8d2*a1**4)*dz
                case (1)
                  rlYlm_laplacian = &
-2.49251322783588191d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(1.96d2*a1**4*a2**4*d2*dy**6-1.862d3*a1**3*a2**4*dy**6 &
-1.862d3*a1**4*a2**3*dy**6+5.88d2*a1**4*a2**4*d2*dx**2*dy**4 &
-5.586d3*a1**3*a2**4*dx**2*dy**4-5.586d3*a1**4*a2**3*dx**2*dy**4 &
-4.2d2*a1**4*a2**4*d2**2*dy**4+4.2d3*a1**3*a2**4*d2*dy**4 &
+4.2d3*a1**4*a2**3*d2*dy**4-1.785d3*a1**2*a2**4*dy**4 &
-3.57d3*a1**3*a2**3*dy**4-1.785d3*a1**4*a2**2*dy**4 &
+5.88d2*a1**4*a2**4*d2*dx**4*dy**2-5.586d3*a1**3*a2**4*dx**4*dy**2 &
-5.586d3*a1**4*a2**3*dx**4*dy**2-8.4d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+8.4d3*a1**3*a2**4*d2*dx**2*dy**2+8.4d3*a1**4*a2**3*d2*dx**2*dy**2 &
-3.57d3*a1**2*a2**4*dx**2*dy**2-7.14d3*a1**3*a2**3*dx**2*dy**2 &
-3.57d3*a1**4*a2**2*dx**2*dy**2+2.88d2*a1**4*a2**4*d2**3*dy**2 &
-3.12d3*a1**3*a2**4*d2**2*dy**2-3.12d3*a1**4*a2**3*d2**2*dy**2 &
+3.48d3*a1**2*a2**4*d2*dy**2+6.96d3*a1**3*a2**3*d2*dy**2 &
+3.48d3*a1**4*a2**2*d2*dy**2-1.62d3*a1*a2**4*dy**2 &
-4.86d3*a1**2*a2**3*dy**2-4.86d3*a1**3*a2**2*dy**2 &
-1.62d3*a1**4*a2*dy**2+1.96d2*a1**4*a2**4*d2*dx**6 &
-1.862d3*a1**3*a2**4*dx**6-1.862d3*a1**4*a2**3*dx**6 &
-4.2d2*a1**4*a2**4*d2**2*dx**4+4.2d3*a1**3*a2**4*d2*dx**4 &
+4.2d3*a1**4*a2**3*d2*dx**4-1.785d3*a1**2*a2**4*dx**4 &
-3.57d3*a1**3*a2**3*dx**4-1.785d3*a1**4*a2**2*dx**4 &
+2.88d2*a1**4*a2**4*d2**3*dx**2-3.12d3*a1**3*a2**4*d2**2*dx**2 &
-3.12d3*a1**4*a2**3*d2**2*dx**2+3.48d3*a1**2*a2**4*d2*dx**2 &
+6.96d3*a1**3*a2**3*d2*dx**2+3.48d3*a1**4*a2**2*d2*dx**2 &
-1.62d3*a1*a2**4*dx**2-4.86d3*a1**2*a2**3*dx**2 &
-4.86d3*a1**3*a2**2*dx**2-1.62d3*a1**4*a2*dx**2 &
-6.4d1*a1**4*a2**4*d2**4+8.0d2*a1**3*a2**4*d2**3 &
+8.0d2*a1**4*a2**3*d2**3-1.92d3*a1**2*a2**4*d2**2 &
-3.84d3*a1**3*a2**3*d2**2-1.92d3*a1**4*a2**2*d2**2+2.28d3*a1*a2**4*d2 &
+6.84d3*a1**2*a2**3*d2+6.84d3*a1**3*a2**2*d2+2.28d3*a1**4*a2*d2 &
-7.8d2*a2**4-3.12d3*a1*a2**3-4.68d3*a1**2*a2**2-3.12d3*a1**3*a2 &
-7.8d2*a1**4)
                case (2)
                  rlYlm_laplacian = &
-3.52494601119984446d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(9.8d1*a1**4*a2**4*d2*dy**6-9.31d2*a1**3*a2**4*dy**6 &
-9.31d2*a1**4*a2**3*dy**6+9.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
-9.31d2*a1**3*a2**4*dx**2*dy**4-9.31d2*a1**4*a2**3*dx**2*dy**4 &
-1.4d2*a1**4*a2**4*d2**2*dy**4+1.302d3*a1**3*a2**4*d2*dy**4 &
+1.302d3*a1**4*a2**3*d2*dy**4+2.38d2*a1**2*a2**4*dy**4 &
+4.76d2*a1**3*a2**3*dy**4+2.38d2*a1**4*a2**2*dy**4 &
-9.8d1*a1**4*a2**4*d2*dx**4*dy**2+9.31d2*a1**3*a2**4*dx**4*dy**2 &
+9.31d2*a1**4*a2**3*dx**4*dy**2-1.96d2*a1**3*a2**4*d2*dx**2*dy**2 &
-1.96d2*a1**4*a2**3*d2*dx**2*dy**2+1.666d3*a1**2*a2**4*dx**2*dy**2 &
+3.332d3*a1**3*a2**3*dx**2*dy**2+1.666d3*a1**4*a2**2*dx**2*dy**2 &
+4.8d1*a1**4*a2**4*d2**3*dy**2-3.8d2*a1**3*a2**4*d2**2*dy**2 &
-3.8d2*a1**4*a2**3*d2**2*dy**2-7.36d2*a1**2*a2**4*d2*dy**2 &
-1.472d3*a1**3*a2**3*d2*dy**2-7.36d2*a1**4*a2**2*d2*dy**2 &
+6.75d2*a1*a2**4*dy**2+2.025d3*a1**2*a2**3*dy**2 &
+2.025d3*a1**3*a2**2*dy**2+6.75d2*a1**4*a2*dy**2 &
-9.8d1*a1**4*a2**4*d2*dx**6+9.31d2*a1**3*a2**4*dx**6 &
+9.31d2*a1**4*a2**3*dx**6+1.4d2*a1**4*a2**4*d2**2*dx**4 &
-1.498d3*a1**3*a2**4*d2*dx**4-1.498d3*a1**4*a2**3*d2*dx**4 &
+1.428d3*a1**2*a2**4*dx**4+2.856d3*a1**3*a2**3*dx**4 &
+1.428d3*a1**4*a2**2*dx**4-4.8d1*a1**4*a2**4*d2**3*dx**2 &
+6.6d2*a1**3*a2**4*d2**2*dx**2+6.6d2*a1**4*a2**3*d2**2*dx**2 &
-1.896d3*a1**2*a2**4*d2*dx**2-3.792d3*a1**3*a2**3*d2*dx**2 &
-1.896d3*a1**4*a2**2*d2*dx**2+1.215d3*a1*a2**4*dx**2 &
+3.645d3*a1**2*a2**3*dx**2+3.645d3*a1**3*a2**2*dx**2 &
+1.215d3*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3 &
-4.8d1*a1**4*a2**3*d2**3+5.16d2*a1**2*a2**4*d2**2 &
+1.032d3*a1**3*a2**3*d2**2+5.16d2*a1**4*a2**2*d2**2-8.64d2*a1*a2**4*d2 &
-2.592d3*a1**2*a2**3*d2-2.592d3*a1**3*a2**2*d2-8.64d2*a1**4*a2*d2 &
+3.51d2*a2**4+1.404d3*a1*a2**3+2.106d3*a1**2*a2**2+1.404d3*a1**3*a2 &
+3.51d2*a1**4)*dz
                case (3)
                  rlYlm_laplacian = &
-6.59457014039261917d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(8.4d1*a1**4*a2**4*d2*dy**6-7.98d2*a1**3*a2**4*dy**6 &
-7.98d2*a1**4*a2**3*dy**6+1.4d2*a1**4*a2**4*d2*dx**2*dy**4 &
-1.33d3*a1**3*a2**4*dx**2*dy**4-1.33d3*a1**4*a2**3*dx**2*dy**4 &
-1.32d2*a1**4*a2**4*d2**2*dy**4+1.152d3*a1**3*a2**4*d2*dy**4 &
+1.152d3*a1**4*a2**3*d2*dy**4+8.67d2*a1**2*a2**4*dy**4 &
+1.734d3*a1**3*a2**3*dy**4+8.67d2*a1**4*a2**2*dy**4 &
+2.8d1*a1**4*a2**4*d2*dx**4*dy**2-2.66d2*a1**3*a2**4*dx**4*dy**2 &
-2.66d2*a1**4*a2**3*dx**4*dy**2-8.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+6.56d2*a1**3*a2**4*d2*dx**2*dy**2+6.56d2*a1**4*a2**3*d2*dx**2*dy**2 &
+1.53d3*a1**2*a2**4*dx**2*dy**2+3.06d3*a1**3*a2**3*dx**2*dy**2 &
+1.53d3*a1**4*a2**2*dx**2*dy**2+4.8d1*a1**4*a2**4*d2**3*dy**2 &
-2.88d2*a1**3*a2**4*d2**2*dy**2-2.88d2*a1**4*a2**3*d2**2*dy**2 &
-1.5d3*a1**2*a2**4*d2*dy**2-3.0d3*a1**3*a2**3*d2*dy**2 &
-1.5d3*a1**4*a2**2*d2*dy**2+5.4d2*a1*a2**4*dy**2 &
+1.62d3*a1**2*a2**3*dy**2+1.62d3*a1**3*a2**2*dy**2 &
+5.4d2*a1**4*a2*dy**2-2.8d1*a1**4*a2**4*d2*dx**6 &
+2.66d2*a1**3*a2**4*dx**6+2.66d2*a1**4*a2**3*dx**6 &
+4.4d1*a1**4*a2**4*d2**2*dx**4-4.96d2*a1**3*a2**4*d2*dx**4 &
-4.96d2*a1**4*a2**3*d2*dx**4+6.63d2*a1**2*a2**4*dx**4 &
+1.326d3*a1**3*a2**3*dx**4+6.63d2*a1**4*a2**2*dx**4 &
-1.6d1*a1**4*a2**4*d2**3*dx**2+2.72d2*a1**3*a2**4*d2**2*dx**2 &
+2.72d2*a1**4*a2**3*d2**2*dx**2-1.02d3*a1**2*a2**4*d2*dx**2 &
-2.04d3*a1**3*a2**3*d2*dx**2-1.02d3*a1**4*a2**2*d2*dx**2 &
-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+4.08d2*a1**2*a2**4*d2**2+8.16d2*a1**3*a2**3*d2**2 &
+4.08d2*a1**4*a2**2*d2**2+3.6d1*a1*a2**4*d2+1.08d2*a1**2*a2**3*d2 &
+1.08d2*a1**3*a2**2*d2+3.6d1*a1**4*a2*d2-2.34d2*a2**4-9.36d2*a1*a2**3 &
-1.404d3*a1**2*a2**2-9.36d2*a1**3*a2-2.34d2*a1**4)
                case (4)
                  rlYlm_laplacian = &
-4.66306526528194375d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(1.4d1*a1**3*a2**3*d2*dy**6-1.33d2*a1**2*a2**3*dy**6 &
-1.33d2*a1**3*a2**2*dy**6-7.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+6.65d2*a1**2*a2**3*dx**2*dy**4+6.65d2*a1**3*a2**2*dx**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**4+2.4d1*a1**2*a2**3*d2*dy**4 &
+2.4d1*a1**3*a2**2*d2*dy**4+4.42d2*a1*a2**3*dy**4 &
+8.84d2*a1**2*a2**2*dy**4+4.42d2*a1**3*a2*dy**4 &
-7.0d1*a1**3*a2**3*d2*dx**4*dy**2+6.65d2*a1**2*a2**3*dx**4*dy**2 &
+6.65d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-2.56d2*a1**2*a2**3*d2*dx**2*dy**2-2.56d2*a1**3*a2**2*d2*dx**2*dy**2 &
-1.7d3*a1*a2**3*dx**2*dy**2-3.4d3*a1**2*a2**2*dx**2*dy**2 &
-1.7d3*a1**3*a2*dx**2*dy**2+1.6d1*a1**2*a2**3*d2**2*dy**2 &
+1.6d1*a1**3*a2**2*d2**2*dy**2-1.0d2*a1*a2**3*d2*dy**2 &
-2.0d2*a1**2*a2**2*d2*dy**2-1.0d2*a1**3*a2*d2*dy**2-2.7d2*a2**3*dy**2 &
-8.1d2*a1*a2**2*dy**2-8.1d2*a1**2*a2*dy**2-2.7d2*a1**3*dy**2 &
+1.4d1*a1**3*a2**3*d2*dx**6-1.33d2*a1**2*a2**3*dx**6 &
-1.33d2*a1**3*a2**2*dx**6-8.0d0*a1**3*a2**3*d2**2*dx**4 &
+1.36d2*a1**2*a2**3*d2*dx**4+1.36d2*a1**3*a2**2*d2*dx**4 &
-5.1d2*a1*a2**3*dx**4-1.02d3*a1**2*a2**2*dx**4-5.1d2*a1**3*a2*dx**4 &
-4.8d1*a1**2*a2**3*d2**2*dx**2-4.8d1*a1**3*a2**2*d2**2*dx**2 &
+3.0d2*a1*a2**3*d2*dx**2+6.0d2*a1**2*a2**2*d2*dx**2 &
+3.0d2*a1**3*a2*d2*dx**2+8.1d2*a2**3*dx**2+2.43d3*a1*a2**2*dx**2 &
+2.43d3*a1**2*a2*dx**2+8.1d2*a1**3*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (0)
          ! selection on l2: l1=4, m1=0
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=0, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =3.32335097044784255d-1*E*a1*a2**5*sqrt &
                  (a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1) &
-1.1d1*a1)*(3.5d1*dy**4+7.0d1*dx**2*dy**2-4.0d1*d2*dy**2+3.5d1*dx**4 &
-4.0d1*d2*dx**2+8.0d0*d2**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=0, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =5.75621273219899775d-1*E*a1*a2**4*sqrt &
                  (a2+a1)*(a2+a1)**(-9)*dy*(7.0d1*a1**2*a2**2*d2*dy**4 &
-4.55d2*a1*a2**2*dy**4-4.55d2*a1**2*a2*dy**4 &
+1.4d2*a1**2*a2**2*d2*dx**2*dy**2-9.1d2*a1*a2**2*dx**2*dy**2 &
-9.1d2*a1**2*a2*dx**2*dy**2-8.0d1*a1**2*a2**2*d2**2*dy**2 &
+4.6d2*a1*a2**2*d2*dy**2+4.6d2*a1**2*a2*d2*dy**2+3.3d2*a2**2*dy**2 &
+6.6d2*a1*a2*dy**2+3.3d2*a1**2*dy**2+7.0d1*a1**2*a2**2*d2*dx**4 &
-4.55d2*a1*a2**2*dx**4-4.55d2*a1**2*a2*dx**4 &
-8.0d1*a1**2*a2**2*d2**2*dx**2+4.6d2*a1*a2**2*d2*dx**2 &
+4.6d2*a1**2*a2*d2*dx**2+3.3d2*a2**2*dx**2+6.6d2*a1*a2*dx**2 &
+3.3d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3-5.6d1*a1*a2**2*d2**2 &
-5.6d1*a1**2*a2*d2**2-2.64d2*a2**2*d2-5.28d2*a1*a2*d2-2.64d2*a1**2*d2)
                case (0)
                  rlYlm_laplacian =5.75621273219899775d-1*E*a1*a2**4*sqrt &
                  (a2+a1)*(a2+a1)**(-9)*(7.0d1*a1**2*a2**2*d2*dy**4 &
-4.55d2*a1*a2**2*dy**4-4.55d2*a1**2*a2*dy**4 &
+1.4d2*a1**2*a2**2*d2*dx**2*dy**2-9.1d2*a1*a2**2*dx**2*dy**2 &
-9.1d2*a1**2*a2*dx**2*dy**2-8.0d1*a1**2*a2**2*d2**2*dy**2 &
+6.0d2*a1*a2**2*d2*dy**2+6.0d2*a1**2*a2*d2*dy**2-4.4d2*a2**2*dy**2 &
-8.8d2*a1*a2*dy**2-4.4d2*a1**2*dy**2+7.0d1*a1**2*a2**2*d2*dx**4 &
-4.55d2*a1*a2**2*dx**4-4.55d2*a1**2*a2*dx**4 &
-8.0d1*a1**2*a2**2*d2**2*dx**2+6.0d2*a1*a2**2*d2*dx**2 &
+6.0d2*a1**2*a2*d2*dx**2-4.4d2*a2**2*dx**2-8.8d2*a1*a2*dx**2 &
-4.4d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3-1.36d2*a1*a2**2*d2**2 &
-1.36d2*a1**2*a2*d2**2+1.76d2*a2**2*d2+3.52d2*a1*a2*d2+1.76d2*a1**2*d2 &
)*dz
                case (1)
                  rlYlm_laplacian =5.75621273219899775d-1*E*a1*a2**4*sqrt &
                  (a2+a1)*(a2+a1)**(-9)*dx*(7.0d1*a1**2*a2**2*d2*dy**4 &
-4.55d2*a1*a2**2*dy**4-4.55d2*a1**2*a2*dy**4 &
+1.4d2*a1**2*a2**2*d2*dx**2*dy**2-9.1d2*a1*a2**2*dx**2*dy**2 &
-9.1d2*a1**2*a2*dx**2*dy**2-8.0d1*a1**2*a2**2*d2**2*dy**2 &
+4.6d2*a1*a2**2*d2*dy**2+4.6d2*a1**2*a2*d2*dy**2+3.3d2*a2**2*dy**2 &
+6.6d2*a1*a2*dy**2+3.3d2*a1**2*dy**2+7.0d1*a1**2*a2**2*d2*dx**4 &
-4.55d2*a1*a2**2*dx**4-4.55d2*a1**2*a2*dx**4 &
-8.0d1*a1**2*a2**2*d2**2*dx**2+4.6d2*a1*a2**2*d2*dx**2 &
+4.6d2*a1**2*a2*d2*dx**2+3.3d2*a2**2*dx**2+6.6d2*a1*a2*dx**2 &
+3.3d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3-5.6d1*a1*a2**2*d2**2 &
-5.6d1*a1**2*a2*d2**2-2.64d2*a2**2*d2-5.28d2*a1*a2*d2-2.64d2*a1**2*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=0, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =1.28712829621467515d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.25d2*a1**2*a2**3*dy**4-5.25d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.05d3*a1**2*a2**3*dx**2*dy**2 &
-1.05d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+4.8d2*a1**2*a2**3*d2*dy**2+4.8d2*a1**3*a2**2*d2*dy**2 &
+7.8d2*a1*a2**3*dy**2+1.56d3*a1**2*a2**2*dy**2+7.8d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.25d2*a1**2*a2**3*dx**4 &
-5.25d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+4.8d2*a1**2*a2**3*d2*dx**2+4.8d2*a1**3*a2**2*d2*dx**2 &
+7.8d2*a1*a2**3*dx**2+1.56d3*a1**2*a2**2*dx**2+7.8d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-2.4d1*a1**2*a2**3*d2**2 &
-2.4d1*a1**3*a2**2*d2**2-6.12d2*a1*a2**3*d2-1.224d3*a1**2*a2**2*d2 &
-6.12d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2-1.98d2*a1**2*a2 &
-6.6d1*a1**3)
                case (-1)
                  rlYlm_laplacian =1.28712829621467515d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.25d2*a1**2*a2**3*dy**4-5.25d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.05d3*a1**2*a2**3*dx**2*dy**2 &
-1.05d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.2d2*a1**2*a2**3*d2*dy**2+6.2d2*a1**3*a2**2*d2*dy**2 &
-1.3d2*a1*a2**3*dy**2-2.6d2*a1**2*a2**2*dy**2-1.3d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.25d2*a1**2*a2**3*dx**4 &
-5.25d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.2d2*a1**2*a2**3*d2*dx**2+6.2d2*a1**3*a2**2*d2*dx**2 &
-1.3d2*a1*a2**3*dx**2-2.6d2*a1**2*a2**2*dx**2-1.3d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.04d2*a1**2*a2**3*d2**2 &
-1.04d2*a1**3*a2**2*d2**2-1.52d2*a1*a2**3*d2-3.04d2*a1**2*a2**2*d2 &
-1.52d2*a1**3*a2*d2+2.64d2*a2**3+7.92d2*a1*a2**2+7.92d2*a1**2*a2 &
+2.64d2*a1**3)*dz
                case (0)
                  rlYlm_laplacian =-3.71561934150563533d-1*E*a1*a2**3*sqrt &
                  (a2+a1)*(a2+a1)**(-10)*(2.1d2*a1**3*a2**3*d2*dy**6 &
-1.575d3*a1**2*a2**3*dy**6-1.575d3*a1**3*a2**2*dy**6 &
+6.3d2*a1**3*a2**3*d2*dx**2*dy**4-4.725d3*a1**2*a2**3*dx**2*dy**4 &
-4.725d3*a1**3*a2**2*dx**2*dy**4-3.8d2*a1**3*a2**3*d2**2*dy**4 &
+3.05d3*a1**2*a2**3*d2*dy**4+3.05d3*a1**3*a2**2*d2*dy**4 &
-1.3d3*a1*a2**3*dy**4-2.6d3*a1**2*a2**2*dy**4-1.3d3*a1**3*a2*dy**4 &
+6.3d2*a1**3*a2**3*d2*dx**4*dy**2-4.725d3*a1**2*a2**3*dx**4*dy**2 &
-4.725d3*a1**3*a2**2*dx**4*dy**2-7.6d2*a1**3*a2**3*d2**2*dx**2*dy**2 &
+6.1d3*a1**2*a2**3*d2*dx**2*dy**2+6.1d3*a1**3*a2**2*d2*dx**2*dy**2 &
-2.6d3*a1*a2**3*dx**2*dy**2-5.2d3*a1**2*a2**2*dx**2*dy**2 &
-2.6d3*a1**3*a2*dx**2*dy**2+2.08d2*a1**3*a2**3*d2**3*dy**2 &
-1.912d3*a1**2*a2**3*d2**2*dy**2-1.912d3*a1**3*a2**2*d2**2*dy**2 &
+2.504d3*a1*a2**3*d2*dy**2+5.008d3*a1**2*a2**2*d2*dy**2 &
+2.504d3*a1**3*a2*d2*dy**2-1.188d3*a2**3*dy**2-3.564d3*a1*a2**2*dy**2 &
-3.564d3*a1**2*a2*dy**2-1.188d3*a1**3*dy**2+2.1d2*a1**3*a2**3*d2*dx**6 &
-1.575d3*a1**2*a2**3*dx**6-1.575d3*a1**3*a2**2*dx**6 &
-3.8d2*a1**3*a2**3*d2**2*dx**4+3.05d3*a1**2*a2**3*d2*dx**4 &
+3.05d3*a1**3*a2**2*d2*dx**4-1.3d3*a1*a2**3*dx**4 &
-2.6d3*a1**2*a2**2*dx**4-1.3d3*a1**3*a2*dx**4 &
+2.08d2*a1**3*a2**3*d2**3*dx**2-1.912d3*a1**2*a2**3*d2**2*dx**2 &
-1.912d3*a1**3*a2**2*d2**2*dx**2+2.504d3*a1*a2**3*d2*dx**2 &
+5.008d3*a1**2*a2**2*d2*dx**2+2.504d3*a1**3*a2*d2*dx**2 &
-1.188d3*a2**3*dx**2-3.564d3*a1*a2**2*dx**2-3.564d3*a1**2*a2*dx**2 &
-1.188d3*a1**3*dx**2-3.2d1*a1**3*a2**3*d2**4+3.68d2*a1**2*a2**3*d2**3 &
+3.68d2*a1**3*a2**2*d2**3-9.76d2*a1*a2**3*d2**2 &
-1.952d3*a1**2*a2**2*d2**2-9.76d2*a1**3*a2*d2**2+7.92d2*a2**3*d2 &
+2.376d3*a1*a2**2*d2+2.376d3*a1**2*a2*d2+7.92d2*a1**3*d2)
                case (1)
                  rlYlm_laplacian =1.28712829621467515d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.25d2*a1**2*a2**3*dy**4-5.25d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.05d3*a1**2*a2**3*dx**2*dy**2 &
-1.05d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.2d2*a1**2*a2**3*d2*dy**2+6.2d2*a1**3*a2**2*d2*dy**2 &
-1.3d2*a1*a2**3*dy**2-2.6d2*a1**2*a2**2*dy**2-1.3d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.25d2*a1**2*a2**3*dx**4 &
-5.25d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.2d2*a1**2*a2**3*d2*dx**2+6.2d2*a1**3*a2**2*d2*dx**2 &
-1.3d2*a1*a2**3*dx**2-2.6d2*a1**2*a2**2*dx**2-1.3d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.04d2*a1**2*a2**3*d2**2 &
-1.04d2*a1**3*a2**2*d2**2-1.52d2*a1*a2**3*d2-3.04d2*a1**2*a2**2*d2 &
-1.52d2*a1**3*a2*d2+2.64d2*a2**3+7.92d2*a1*a2**2+7.92d2*a1**2*a2 &
+2.64d2*a1**3)*dz
                case (2)
                  rlYlm_laplacian =-6.43564148107337574d-1*E*a1*a2**3*sqrt &
                  (a2+a1)*(a2+a1)**(-10)*(dy+dx)*(dy-1.0d0*dx)* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.25d2*a1**2*a2**3*dy**4 &
-5.25d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.05d3*a1**2*a2**3*dx**2*dy**2-1.05d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+4.8d2*a1**2*a2**3*d2*dy**2 &
+4.8d2*a1**3*a2**2*d2*dy**2+7.8d2*a1*a2**3*dy**2 &
+1.56d3*a1**2*a2**2*dy**2+7.8d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.25d2*a1**2*a2**3*dx**4 &
-5.25d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+4.8d2*a1**2*a2**3*d2*dx**2+4.8d2*a1**3*a2**2*d2*dx**2 &
+7.8d2*a1*a2**3*dx**2+1.56d3*a1**2*a2**2*dx**2+7.8d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-2.4d1*a1**2*a2**3*d2**2 &
-2.4d1*a1**3*a2**2*d2**2-6.12d2*a1*a2**3*d2-1.224d3*a1**2*a2**2*d2 &
-6.12d2*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2-1.98d2*a1**2*a2 &
-6.6d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=0, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =-6.95128727779234418d-1*E*a1**2*a2**3*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*dy*(dy**2-3.0d0*dx**2)* &
(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+5.0d2*a1**2*a2**3*d2*dy**2 &
+5.0d2*a1**3*a2**2*d2*dy**2+1.35d3*a1*a2**3*dy**2 &
+2.7d3*a1**2*a2**2*dy**2+1.35d3*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+5.0d2*a1**2*a2**3*d2*dx**2+5.0d2*a1**3*a2**2*d2*dx**2 &
+1.35d3*a1*a2**3*dx**2+2.7d3*a1**2*a2**2*dx**2+1.35d3*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3+8.0d0*a1**2*a2**3*d2**2 &
+8.0d0*a1**3*a2**2*d2**2-1.044d3*a1*a2**3*d2-2.088d3*a1**2*a2**2*d2 &
-1.044d3*a1**3*a2*d2-2.34d2*a2**3-7.02d2*a1*a2**2-7.02d2*a1**2*a2 &
-2.34d2*a1**3)
                case (-2)
                  rlYlm_laplacian =3.40542137721830949d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.95d2*a1**2*a2**3*dy**4-5.95d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.19d3*a1**2*a2**3*dx**2*dy**2 &
-1.19d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.4d2*a1**2*a2**3*d2*dy**2+6.4d2*a1**3*a2**2*d2*dy**2 &
+3.0d2*a1*a2**3*dy**2+6.0d2*a1**2*a2**2*dy**2+3.0d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.4d2*a1**2*a2**3*d2*dx**2+6.4d2*a1**3*a2**2*d2*dx**2 &
+3.0d2*a1*a2**3*dx**2+6.0d2*a1**2*a2**2*dx**2+3.0d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-7.2d1*a1**2*a2**3*d2**2 &
-7.2d1*a1**3*a2**2*d2**2-5.64d2*a1*a2**3*d2-1.128d3*a1**2*a2**2*d2 &
-5.64d2*a1**3*a2*d2+5.46d2*a2**3+1.638d3*a1*a2**2+1.638d3*a1**2*a2 &
+5.46d2*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =-5.3844439723186478d-1*E*a1*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*dy* &
(3.5d2*a1**4*a2**4*d2*dy**6-2.975d3*a1**3*a2**4*dy**6 &
-2.975d3*a1**4*a2**3*dy**6+1.05d3*a1**4*a2**4*d2*dx**2*dy**4 &
-8.925d3*a1**3*a2**4*dx**2*dy**4-8.925d3*a1**4*a2**3*dx**2*dy**4 &
-6.8d2*a1**4*a2**4*d2**2*dy**4+6.0d3*a1**3*a2**4*d2*dy**4 &
+6.0d3*a1**4*a2**3*d2*dy**4-1.65d3*a1**2*a2**4*dy**4 &
-3.3d3*a1**3*a2**3*dy**4-1.65d3*a1**4*a2**2*dy**4 &
+1.05d3*a1**4*a2**4*d2*dx**4*dy**2-8.925d3*a1**3*a2**4*dx**4*dy**2 &
-8.925d3*a1**4*a2**3*dx**4*dy**2-1.36d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.2d4*a1**3*a2**4*d2*dx**2*dy**2+1.2d4*a1**4*a2**3*d2*dx**2*dy**2 &
-3.3d3*a1**2*a2**4*dx**2*dy**2-6.6d3*a1**3*a2**3*dx**2*dy**2 &
-3.3d3*a1**4*a2**2*dx**2*dy**2+4.0d2*a1**4*a2**4*d2**3*dy**2 &
-3.72d3*a1**3*a2**4*d2**2*dy**2-3.72d3*a1**4*a2**3*d2**2*dy**2 &
+2.46d3*a1**2*a2**4*d2*dy**2+4.92d3*a1**3*a2**3*d2*dy**2 &
+2.46d3*a1**4*a2**2*d2*dy**2-3.9d2*a1*a2**4*dy**2 &
-1.17d3*a1**2*a2**3*dy**2-1.17d3*a1**3*a2**2*dy**2 &
-3.9d2*a1**4*a2*dy**2+3.5d2*a1**4*a2**4*d2*dx**6 &
-2.975d3*a1**3*a2**4*dx**6-2.975d3*a1**4*a2**3*dx**6 &
-6.8d2*a1**4*a2**4*d2**2*dx**4+6.0d3*a1**3*a2**4*d2*dx**4 &
+6.0d3*a1**4*a2**3*d2*dx**4-1.65d3*a1**2*a2**4*dx**4 &
-3.3d3*a1**3*a2**3*dx**4-1.65d3*a1**4*a2**2*dx**4 &
+4.0d2*a1**4*a2**4*d2**3*dx**2-3.72d3*a1**3*a2**4*d2**2*dx**2 &
-3.72d3*a1**4*a2**3*d2**2*dx**2+2.46d3*a1**2*a2**4*d2*dx**2 &
+4.92d3*a1**3*a2**3*d2*dx**2+2.46d3*a1**4*a2**2*d2*dx**2 &
-3.9d2*a1*a2**4*dx**2-1.17d3*a1**2*a2**3*dx**2 &
-1.17d3*a1**3*a2**2*dx**2-3.9d2*a1**4*a2*dx**2-6.4d1*a1**4*a2**4*d2**4 &
+6.08d2*a1**3*a2**4*d2**3+6.08d2*a1**4*a2**3*d2**3 &
-3.84d2*a1**2*a2**4*d2**2-7.68d2*a1**3*a2**3*d2**2 &
-3.84d2*a1**4*a2**2*d2**2-7.44d2*a1*a2**4*d2-2.232d3*a1**2*a2**3*d2 &
-2.232d3*a1**3*a2**2*d2-7.44d2*a1**4*a2*d2+6.6d2*a2**4+2.64d3*a1*a2**3 &
+3.96d3*a1**2*a2**2+2.64d3*a1**3*a2+6.6d2*a1**4)
                case (0)
                  rlYlm_laplacian =-4.39638009359507944d-1*E*a1*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*(3.5d2*a1**4*a2**4*d2*dy**6 &
-2.975d3*a1**3*a2**4*dy**6-2.975d3*a1**4*a2**3*dy**6 &
+1.05d3*a1**4*a2**4*d2*dx**2*dy**4-8.925d3*a1**3*a2**4*dx**2*dy**4 &
-8.925d3*a1**4*a2**3*dx**2*dy**4-5.4d2*a1**4*a2**4*d2**2*dy**4 &
+4.95d3*a1**3*a2**4*d2*dy**4+4.95d3*a1**4*a2**3*d2*dy**4 &
-2.7d3*a1**2*a2**4*dy**4-5.4d3*a1**3*a2**3*dy**4 &
-2.7d3*a1**4*a2**2*dy**4+1.05d3*a1**4*a2**4*d2*dx**4*dy**2 &
-8.925d3*a1**3*a2**4*dx**4*dy**2-8.925d3*a1**4*a2**3*dx**4*dy**2 &
-1.08d3*a1**4*a2**4*d2**2*dx**2*dy**2+9.9d3*a1**3*a2**4*d2*dx**2*dy**2 &
+9.9d3*a1**4*a2**3*d2*dx**2*dy**2-5.4d3*a1**2*a2**4*dx**2*dy**2 &
-1.08d4*a1**3*a2**3*dx**2*dy**2-5.4d3*a1**4*a2**2*dx**2*dy**2 &
+2.4d2*a1**4*a2**4*d2**3*dy**2-2.52d3*a1**3*a2**4*d2**2*dy**2 &
-2.52d3*a1**4*a2**3*d2**2*dy**2+3.96d3*a1**2*a2**4*d2*dy**2 &
+7.92d3*a1**3*a2**3*d2*dy**2+3.96d3*a1**4*a2**2*d2*dy**2 &
-2.34d3*a1*a2**4*dy**2-7.02d3*a1**2*a2**3*dy**2 &
-7.02d3*a1**3*a2**2*dy**2-2.34d3*a1**4*a2*dy**2 &
+3.5d2*a1**4*a2**4*d2*dx**6-2.975d3*a1**3*a2**4*dx**6 &
-2.975d3*a1**4*a2**3*dx**6-5.4d2*a1**4*a2**4*d2**2*dx**4 &
+4.95d3*a1**3*a2**4*d2*dx**4+4.95d3*a1**4*a2**3*d2*dx**4 &
-2.7d3*a1**2*a2**4*dx**4-5.4d3*a1**3*a2**3*dx**4 &
-2.7d3*a1**4*a2**2*dx**4+2.4d2*a1**4*a2**4*d2**3*dx**2 &
-2.52d3*a1**3*a2**4*d2**2*dx**2-2.52d3*a1**4*a2**3*d2**2*dx**2 &
+3.96d3*a1**2*a2**4*d2*dx**2+7.92d3*a1**3*a2**3*d2*dx**2 &
+3.96d3*a1**4*a2**2*d2*dx**2-2.34d3*a1*a2**4*dx**2 &
-7.02d3*a1**2*a2**3*dx**2-7.02d3*a1**3*a2**2*dx**2 &
-2.34d3*a1**4*a2*dx**2-3.2d1*a1**4*a2**4*d2**4 &
+4.64d2*a1**3*a2**4*d2**3+4.64d2*a1**4*a2**3*d2**3 &
-1.872d3*a1**2*a2**4*d2**2-3.744d3*a1**3*a2**3*d2**2 &
-1.872d3*a1**4*a2**2*d2**2+3.048d3*a1*a2**4*d2+9.144d3*a1**2*a2**3*d2 &
+9.144d3*a1**3*a2**2*d2+3.048d3*a1**4*a2*d2-1.32d3*a2**4 &
-5.28d3*a1*a2**3-7.92d3*a1**2*a2**2-5.28d3*a1**3*a2-1.32d3*a1**4)*dz
                case (1)
                  rlYlm_laplacian =-5.3844439723186478d-1*E*a1*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*dx* &
(3.5d2*a1**4*a2**4*d2*dy**6-2.975d3*a1**3*a2**4*dy**6 &
-2.975d3*a1**4*a2**3*dy**6+1.05d3*a1**4*a2**4*d2*dx**2*dy**4 &
-8.925d3*a1**3*a2**4*dx**2*dy**4-8.925d3*a1**4*a2**3*dx**2*dy**4 &
-6.8d2*a1**4*a2**4*d2**2*dy**4+6.0d3*a1**3*a2**4*d2*dy**4 &
+6.0d3*a1**4*a2**3*d2*dy**4-1.65d3*a1**2*a2**4*dy**4 &
-3.3d3*a1**3*a2**3*dy**4-1.65d3*a1**4*a2**2*dy**4 &
+1.05d3*a1**4*a2**4*d2*dx**4*dy**2-8.925d3*a1**3*a2**4*dx**4*dy**2 &
-8.925d3*a1**4*a2**3*dx**4*dy**2-1.36d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.2d4*a1**3*a2**4*d2*dx**2*dy**2+1.2d4*a1**4*a2**3*d2*dx**2*dy**2 &
-3.3d3*a1**2*a2**4*dx**2*dy**2-6.6d3*a1**3*a2**3*dx**2*dy**2 &
-3.3d3*a1**4*a2**2*dx**2*dy**2+4.0d2*a1**4*a2**4*d2**3*dy**2 &
-3.72d3*a1**3*a2**4*d2**2*dy**2-3.72d3*a1**4*a2**3*d2**2*dy**2 &
+2.46d3*a1**2*a2**4*d2*dy**2+4.92d3*a1**3*a2**3*d2*dy**2 &
+2.46d3*a1**4*a2**2*d2*dy**2-3.9d2*a1*a2**4*dy**2 &
-1.17d3*a1**2*a2**3*dy**2-1.17d3*a1**3*a2**2*dy**2 &
-3.9d2*a1**4*a2*dy**2+3.5d2*a1**4*a2**4*d2*dx**6 &
-2.975d3*a1**3*a2**4*dx**6-2.975d3*a1**4*a2**3*dx**6 &
-6.8d2*a1**4*a2**4*d2**2*dx**4+6.0d3*a1**3*a2**4*d2*dx**4 &
+6.0d3*a1**4*a2**3*d2*dx**4-1.65d3*a1**2*a2**4*dx**4 &
-3.3d3*a1**3*a2**3*dx**4-1.65d3*a1**4*a2**2*dx**4 &
+4.0d2*a1**4*a2**4*d2**3*dx**2-3.72d3*a1**3*a2**4*d2**2*dx**2 &
-3.72d3*a1**4*a2**3*d2**2*dx**2+2.46d3*a1**2*a2**4*d2*dx**2 &
+4.92d3*a1**3*a2**3*d2*dx**2+2.46d3*a1**4*a2**2*d2*dx**2 &
-3.9d2*a1*a2**4*dx**2-1.17d3*a1**2*a2**3*dx**2 &
-1.17d3*a1**3*a2**2*dx**2-3.9d2*a1**4*a2*dx**2-6.4d1*a1**4*a2**4*d2**4 &
+6.08d2*a1**3*a2**4*d2**3+6.08d2*a1**4*a2**3*d2**3 &
-3.84d2*a1**2*a2**4*d2**2-7.68d2*a1**3*a2**3*d2**2 &
-3.84d2*a1**4*a2**2*d2**2-7.44d2*a1*a2**4*d2-2.232d3*a1**2*a2**3*d2 &
-2.232d3*a1**3*a2**2*d2-7.44d2*a1**4*a2*d2+6.6d2*a2**4+2.64d3*a1*a2**3 &
+3.96d3*a1**2*a2**2+2.64d3*a1**3*a2+6.6d2*a1**4)
                case (2)
                  rlYlm_laplacian = &
-1.70271068860915474d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*(dy+dx &
)*(dy-1.0d0*dx)*(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+6.4d2*a1**2*a2**3*d2*dy**2 &
+6.4d2*a1**3*a2**2*d2*dy**2+3.0d2*a1*a2**3*dy**2 &
+6.0d2*a1**2*a2**2*dy**2+3.0d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.4d2*a1**2*a2**3*d2*dx**2+6.4d2*a1**3*a2**2*d2*dx**2 &
+3.0d2*a1*a2**3*dx**2+6.0d2*a1**2*a2**2*dx**2+3.0d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-7.2d1*a1**2*a2**3*d2**2 &
-7.2d1*a1**3*a2**2*d2**2-5.64d2*a1*a2**3*d2-1.128d3*a1**2*a2**2*d2 &
-5.64d2*a1**3*a2*d2+5.46d2*a2**3+1.638d3*a1*a2**2+1.638d3*a1**2*a2 &
+5.46d2*a1**3)*dz
                case (3)
                  rlYlm_laplacian =-6.95128727779234418d-1*E*a1**2*a2**3*sqrt &
                  (a2+a1)*(a2+a1)**(-11)*dx*(3.0d0*dy**2 &
-1.0d0*dx**2)*(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+5.0d2*a1**2*a2**3*d2*dy**2 &
+5.0d2*a1**3*a2**2*d2*dy**2+1.35d3*a1*a2**3*dy**2 &
+2.7d3*a1**2*a2**2*dy**2+1.35d3*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+5.0d2*a1**2*a2**3*d2*dx**2+5.0d2*a1**3*a2**2*d2*dx**2 &
+1.35d3*a1*a2**3*dx**2+2.7d3*a1**2*a2**2*dx**2+1.35d3*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3+8.0d0*a1**2*a2**3*d2**2 &
+8.0d0*a1**3*a2**2*d2**2-1.044d3*a1*a2**3*d2-2.088d3*a1**2*a2**2*d2 &
-1.044d3*a1**3*a2*d2-2.34d2*a2**3-7.02d2*a1*a2**2-7.02d2*a1**2*a2 &
-2.34d2*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=0, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-2.94918142326164563d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(7.0d1*a1**3*a2**3*d2*dy**4 &
-6.65d2*a1**2*a2**3*dy**4-6.65d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.33d3*a1**2*a2**3*dx**2*dy**2 &
-1.33d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+5.2d2*a1**2*a2**3*d2*dy**2+5.2d2*a1**3*a2**2*d2*dy**2 &
+2.04d3*a1*a2**3*dy**2+4.08d3*a1**2*a2**2*dy**2+2.04d3*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-6.65d2*a1**2*a2**3*dx**4 &
-6.65d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+5.2d2*a1**2*a2**3*d2*dx**2+5.2d2*a1**3*a2**2*d2*dx**2 &
+2.04d3*a1*a2**3*dx**2+4.08d3*a1**2*a2**2*dx**2+2.04d3*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3+4.0d1*a1**2*a2**3*d2**2 &
+4.0d1*a1**3*a2**2*d2**2-1.56d3*a1*a2**3*d2-3.12d3*a1**2*a2**2*d2 &
-1.56d3*a1**3*a2*d2-5.4d2*a2**3-1.62d3*a1*a2**2-1.62d3*a1**2*a2 &
-5.4d2*a1**3)
                case (-3)
                  rlYlm_laplacian = &
-2.08538618333770325d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(dy**2-3.0d0*dx**2)*(7.0d1*a1**3*a2**3*d2*dy**4 &
-6.65d2*a1**2*a2**3*dy**4-6.65d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.33d3*a1**2*a2**3*dx**2*dy**2 &
-1.33d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.6d2*a1**2*a2**3*d2*dy**2+6.6d2*a1**3*a2**2*d2*dy**2 &
+8.5d2*a1*a2**3*dy**2+1.7d3*a1**2*a2**2*dy**2+8.5d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-6.65d2*a1**2*a2**3*dx**4 &
-6.65d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.6d2*a1**2*a2**3*d2*dx**2+6.6d2*a1**3*a2**2*d2*dx**2 &
+8.5d2*a1*a2**3*dx**2+1.7d3*a1**2*a2**2*dx**2+8.5d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-4.0d1*a1**2*a2**3*d2**2 &
-4.0d1*a1**3*a2**2*d2**2-1.06d3*a1*a2**3*d2-2.12d3*a1**2*a2**2*d2 &
-1.06d3*a1**3*a2*d2+8.1d2*a2**3+2.43d3*a1*a2**2+2.43d3*a1**2*a2 &
+8.1d2*a1**3)*dz
                case (-2)
                  rlYlm_laplacian = &
-1.1146858024516906d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(4.9d2*a1**4*a2**4*d2*dy**6-4.655d3*a1**3*a2**4*dy**6 &
-4.655d3*a1**4*a2**3*dy**6+1.47d3*a1**4*a2**4*d2*dx**2*dy**4 &
-1.3965d4*a1**3*a2**4*dx**2*dy**4-1.3965d4*a1**4*a2**3*dx**2*dy**4 &
-9.8d2*a1**4*a2**4*d2**2*dy**4+9.31d3*a1**3*a2**4*d2*dy**4 &
+9.31d3*a1**4*a2**3*d2*dy**4+1.47d3*a1**4*a2**4*d2*dx**4*dy**2 &
-1.3965d4*a1**3*a2**4*dx**4*dy**2-1.3965d4*a1**4*a2**3*dx**4*dy**2 &
-1.96d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.862d4*a1**3*a2**4*d2*dx**2*dy**2+1.862d4*a1**4*a2**3*d2*dx**2*dy**2 &
+5.92d2*a1**4*a2**4*d2**3*dy**2-5.48d3*a1**3*a2**4*d2**2*dy**2 &
-5.48d3*a1**4*a2**3*d2**2*dy**2-1.62d3*a1**2*a2**4*d2*dy**2 &
-3.24d3*a1**3*a2**3*d2*dy**2-1.62d3*a1**4*a2**2*d2*dy**2 &
+2.97d3*a1*a2**4*dy**2+8.91d3*a1**2*a2**3*dy**2 &
+8.91d3*a1**3*a2**2*dy**2+2.97d3*a1**4*a2*dy**2 &
+4.9d2*a1**4*a2**4*d2*dx**6-4.655d3*a1**3*a2**4*dx**6 &
-4.655d3*a1**4*a2**3*dx**6-9.8d2*a1**4*a2**4*d2**2*dx**4 &
+9.31d3*a1**3*a2**4*d2*dx**4+9.31d3*a1**4*a2**3*d2*dx**4 &
+5.92d2*a1**4*a2**4*d2**3*dx**2-5.48d3*a1**3*a2**4*d2**2*dx**2 &
-5.48d3*a1**4*a2**3*d2**2*dx**2-1.62d3*a1**2*a2**4*d2*dx**2 &
-3.24d3*a1**3*a2**3*d2*dx**2-1.62d3*a1**4*a2**2*d2*dx**2 &
+2.97d3*a1*a2**4*dx**2+8.91d3*a1**2*a2**3*dx**2 &
+8.91d3*a1**3*a2**2*dx**2+2.97d3*a1**4*a2*dx**2 &
-9.6d1*a1**4*a2**4*d2**4+7.2d2*a1**3*a2**4*d2**3 &
+7.2d2*a1**4*a2**3*d2**3+2.28d3*a1**2*a2**4*d2**2 &
+4.56d3*a1**3*a2**3*d2**2+2.28d3*a1**4*a2**2*d2**2-5.22d3*a1*a2**4*d2 &
-1.566d4*a1**2*a2**3*d2-1.566d4*a1**3*a2**2*d2-5.22d3*a1**4*a2*d2 &
+2.34d3*a2**4+9.36d3*a1*a2**3+1.404d4*a1**2*a2**2+9.36d3*a1**3*a2 &
+2.34d3*a1**4)
                case (-1)
                  rlYlm_laplacian =-7.88201889805958723d-1*E*a1**2*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-12)*dy* &
(4.9d2*a1**4*a2**4*d2*dy**6-4.655d3*a1**3*a2**4*dy**6 &
-4.655d3*a1**4*a2**3*dy**6+1.47d3*a1**4*a2**4*d2*dx**2*dy**4 &
-1.3965d4*a1**3*a2**4*dx**2*dy**4-1.3965d4*a1**4*a2**3*dx**2*dy**4 &
-8.4d2*a1**4*a2**4*d2**2*dy**4+8.4d3*a1**3*a2**4*d2*dy**4 &
+8.4d3*a1**4*a2**3*d2*dy**4-3.57d3*a1**2*a2**4*dy**4 &
-7.14d3*a1**3*a2**3*dy**4-3.57d3*a1**4*a2**2*dy**4 &
+1.47d3*a1**4*a2**4*d2*dx**4*dy**2-1.3965d4*a1**3*a2**4*dx**4*dy**2 &
-1.3965d4*a1**4*a2**3*dx**4*dy**2-1.68d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.68d4*a1**3*a2**4*d2*dx**2*dy**2+1.68d4*a1**4*a2**3*d2*dx**2*dy**2 &
-7.14d3*a1**2*a2**4*dx**2*dy**2-1.428d4*a1**3*a2**3*dx**2*dy**2 &
-7.14d3*a1**4*a2**2*dx**2*dy**2+4.32d2*a1**4*a2**4*d2**3*dy**2 &
-4.68d3*a1**3*a2**4*d2**2*dy**2-4.68d3*a1**4*a2**3*d2**2*dy**2 &
+5.22d3*a1**2*a2**4*d2*dy**2+1.044d4*a1**3*a2**3*d2*dy**2 &
+5.22d3*a1**4*a2**2*d2*dy**2-2.43d3*a1*a2**4*dy**2 &
-7.29d3*a1**2*a2**3*dy**2-7.29d3*a1**3*a2**2*dy**2 &
-2.43d3*a1**4*a2*dy**2+4.9d2*a1**4*a2**4*d2*dx**6 &
-4.655d3*a1**3*a2**4*dx**6-4.655d3*a1**4*a2**3*dx**6 &
-8.4d2*a1**4*a2**4*d2**2*dx**4+8.4d3*a1**3*a2**4*d2*dx**4 &
+8.4d3*a1**4*a2**3*d2*dx**4-3.57d3*a1**2*a2**4*dx**4 &
-7.14d3*a1**3*a2**3*dx**4-3.57d3*a1**4*a2**2*dx**4 &
+4.32d2*a1**4*a2**4*d2**3*dx**2-4.68d3*a1**3*a2**4*d2**2*dx**2 &
-4.68d3*a1**4*a2**3*d2**2*dx**2+5.22d3*a1**2*a2**4*d2*dx**2 &
+1.044d4*a1**3*a2**3*d2*dx**2+5.22d3*a1**4*a2**2*d2*dx**2 &
-2.43d3*a1*a2**4*dx**2-7.29d3*a1**2*a2**3*dx**2 &
-7.29d3*a1**3*a2**2*dx**2-2.43d3*a1**4*a2*dx**2 &
-6.4d1*a1**4*a2**4*d2**4+8.0d2*a1**3*a2**4*d2**3 &
+8.0d2*a1**4*a2**3*d2**3-1.92d3*a1**2*a2**4*d2**2 &
-3.84d3*a1**3*a2**3*d2**2-1.92d3*a1**4*a2**2*d2**2+2.28d3*a1*a2**4*d2 &
+6.84d3*a1**2*a2**3*d2+6.84d3*a1**3*a2**2*d2+2.28d3*a1**4*a2*d2 &
-7.8d2*a2**4-3.12d3*a1*a2**3-4.68d3*a1**2*a2**2-3.12d3*a1**3*a2 &
-7.8d2*a1**4)*dz
                case (0)
                  rlYlm_laplacian =1.24625661391794096d-1*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(2.45d3*a1**5*a2**5*d2*dy**8 &
-2.3275d4*a1**4*a2**5*dy**8-2.3275d4*a1**5*a2**4*dy**8 &
+9.8d3*a1**5*a2**5*d2*dx**2*dy**6-9.31d4*a1**4*a2**5*dx**2*dy**6 &
-9.31d4*a1**5*a2**4*dx**2*dy**6-5.6d3*a1**5*a2**5*d2**2*dy**6 &
+5.6d4*a1**4*a2**5*d2*dy**6+5.6d4*a1**5*a2**4*d2*dy**6 &
-2.38d4*a1**3*a2**5*dy**6-4.76d4*a1**4*a2**4*dy**6 &
-2.38d4*a1**5*a2**3*dy**6+1.47d4*a1**5*a2**5*d2*dx**4*dy**4 &
-1.3965d5*a1**4*a2**5*dx**4*dy**4-1.3965d5*a1**5*a2**4*dx**4*dy**4 &
-1.68d4*a1**5*a2**5*d2**2*dx**2*dy**4 &
+1.68d5*a1**4*a2**5*d2*dx**2*dy**4+1.68d5*a1**5*a2**4*d2*dx**2*dy**4 &
-7.14d4*a1**3*a2**5*dx**2*dy**4-1.428d5*a1**4*a2**4*dx**2*dy**4 &
-7.14d4*a1**5*a2**3*dx**2*dy**4+4.32d3*a1**5*a2**5*d2**3*dy**4 &
-4.68d4*a1**4*a2**5*d2**2*dy**4-4.68d4*a1**5*a2**4*d2**2*dy**4 &
+5.22d4*a1**3*a2**5*d2*dy**4+1.044d5*a1**4*a2**4*d2*dy**4 &
+5.22d4*a1**5*a2**3*d2*dy**4-2.43d4*a1**2*a2**5*dy**4 &
-7.29d4*a1**3*a2**4*dy**4-7.29d4*a1**4*a2**3*dy**4 &
-2.43d4*a1**5*a2**2*dy**4+9.8d3*a1**5*a2**5*d2*dx**6*dy**2 &
-9.31d4*a1**4*a2**5*dx**6*dy**2-9.31d4*a1**5*a2**4*dx**6*dy**2 &
-1.68d4*a1**5*a2**5*d2**2*dx**4*dy**2 &
+1.68d5*a1**4*a2**5*d2*dx**4*dy**2+1.68d5*a1**5*a2**4*d2*dx**4*dy**2 &
-7.14d4*a1**3*a2**5*dx**4*dy**2-1.428d5*a1**4*a2**4*dx**4*dy**2 &
-7.14d4*a1**5*a2**3*dx**4*dy**2+8.64d3*a1**5*a2**5*d2**3*dx**2*dy**2 &
-9.36d4*a1**4*a2**5*d2**2*dx**2*dy**2 &
-9.36d4*a1**5*a2**4*d2**2*dx**2*dy**2 &
+1.044d5*a1**3*a2**5*d2*dx**2*dy**2+2.088d5*a1**4*a2**4*d2*dx**2*dy**2 &
+1.044d5*a1**5*a2**3*d2*dx**2*dy**2-4.86d4*a1**2*a2**5*dx**2*dy**2 &
-1.458d5*a1**3*a2**4*dx**2*dy**2-1.458d5*a1**4*a2**3*dx**2*dy**2 &
-4.86d4*a1**5*a2**2*dx**2*dy**2-1.28d3*a1**5*a2**5*d2**4*dy**2 &
+1.6d4*a1**4*a2**5*d2**3*dy**2+1.6d4*a1**5*a2**4*d2**3*dy**2 &
-3.84d4*a1**3*a2**5*d2**2*dy**2-7.68d4*a1**4*a2**4*d2**2*dy**2 &
-3.84d4*a1**5*a2**3*d2**2*dy**2+4.56d4*a1**2*a2**5*d2*dy**2 &
+1.368d5*a1**3*a2**4*d2*dy**2+1.368d5*a1**4*a2**3*d2*dy**2 &
+4.56d4*a1**5*a2**2*d2*dy**2-1.56d4*a1*a2**5*dy**2 &
-6.24d4*a1**2*a2**4*dy**2-9.36d4*a1**3*a2**3*dy**2 &
-6.24d4*a1**4*a2**2*dy**2-1.56d4*a1**5*a2*dy**2 &
+2.45d3*a1**5*a2**5*d2*dx**8-2.3275d4*a1**4*a2**5*dx**8 &
-2.3275d4*a1**5*a2**4*dx**8-5.6d3*a1**5*a2**5*d2**2*dx**6 &
+5.6d4*a1**4*a2**5*d2*dx**6+5.6d4*a1**5*a2**4*d2*dx**6 &
-2.38d4*a1**3*a2**5*dx**6-4.76d4*a1**4*a2**4*dx**6 &
-2.38d4*a1**5*a2**3*dx**6+4.32d3*a1**5*a2**5*d2**3*dx**4 &
-4.68d4*a1**4*a2**5*d2**2*dx**4-4.68d4*a1**5*a2**4*d2**2*dx**4 &
+5.22d4*a1**3*a2**5*d2*dx**4+1.044d5*a1**4*a2**4*d2*dx**4 &
+5.22d4*a1**5*a2**3*d2*dx**4-2.43d4*a1**2*a2**5*dx**4 &
-7.29d4*a1**3*a2**4*dx**4-7.29d4*a1**4*a2**3*dx**4 &
-2.43d4*a1**5*a2**2*dx**4-1.28d3*a1**5*a2**5*d2**4*dx**2 &
+1.6d4*a1**4*a2**5*d2**3*dx**2+1.6d4*a1**5*a2**4*d2**3*dx**2 &
-3.84d4*a1**3*a2**5*d2**2*dx**2-7.68d4*a1**4*a2**4*d2**2*dx**2 &
-3.84d4*a1**5*a2**3*d2**2*dx**2+4.56d4*a1**2*a2**5*d2*dx**2 &
+1.368d5*a1**3*a2**4*d2*dx**2+1.368d5*a1**4*a2**3*d2*dx**2 &
+4.56d4*a1**5*a2**2*d2*dx**2-1.56d4*a1*a2**5*dx**2 &
-6.24d4*a1**2*a2**4*dx**2-9.36d4*a1**3*a2**3*dx**2 &
-6.24d4*a1**4*a2**2*dx**2-1.56d4*a1**5*a2*dx**2 &
+1.28d2*a1**5*a2**5*d2**5-2.24d3*a1**4*a2**5*d2**4 &
-2.24d3*a1**5*a2**4*d2**4+1.216d4*a1**3*a2**5*d2**3 &
+2.432d4*a1**4*a2**4*d2**3+1.216d4*a1**5*a2**3*d2**3 &
-2.976d4*a1**2*a2**5*d2**2-8.928d4*a1**3*a2**4*d2**2 &
-8.928d4*a1**4*a2**3*d2**2-2.976d4*a1**5*a2**2*d2**2 &
+2.58d4*a1*a2**5*d2+1.032d5*a1**2*a2**4*d2+1.548d5*a1**3*a2**3*d2 &
+1.032d5*a1**4*a2**2*d2+2.58d4*a1**5*a2*d2-4.62d3*a2**5 &
-2.31d4*a1*a2**4-4.62d4*a1**2*a2**3-4.62d4*a1**3*a2**2-2.31d4*a1**4*a2 &
-4.62d3*a1**5)
                case (1)
                  rlYlm_laplacian =-7.88201889805958723d-1*E*a1**2*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-12)*dx* &
(4.9d2*a1**4*a2**4*d2*dy**6-4.655d3*a1**3*a2**4*dy**6 &
-4.655d3*a1**4*a2**3*dy**6+1.47d3*a1**4*a2**4*d2*dx**2*dy**4 &
-1.3965d4*a1**3*a2**4*dx**2*dy**4-1.3965d4*a1**4*a2**3*dx**2*dy**4 &
-8.4d2*a1**4*a2**4*d2**2*dy**4+8.4d3*a1**3*a2**4*d2*dy**4 &
+8.4d3*a1**4*a2**3*d2*dy**4-3.57d3*a1**2*a2**4*dy**4 &
-7.14d3*a1**3*a2**3*dy**4-3.57d3*a1**4*a2**2*dy**4 &
+1.47d3*a1**4*a2**4*d2*dx**4*dy**2-1.3965d4*a1**3*a2**4*dx**4*dy**2 &
-1.3965d4*a1**4*a2**3*dx**4*dy**2-1.68d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.68d4*a1**3*a2**4*d2*dx**2*dy**2+1.68d4*a1**4*a2**3*d2*dx**2*dy**2 &
-7.14d3*a1**2*a2**4*dx**2*dy**2-1.428d4*a1**3*a2**3*dx**2*dy**2 &
-7.14d3*a1**4*a2**2*dx**2*dy**2+4.32d2*a1**4*a2**4*d2**3*dy**2 &
-4.68d3*a1**3*a2**4*d2**2*dy**2-4.68d3*a1**4*a2**3*d2**2*dy**2 &
+5.22d3*a1**2*a2**4*d2*dy**2+1.044d4*a1**3*a2**3*d2*dy**2 &
+5.22d3*a1**4*a2**2*d2*dy**2-2.43d3*a1*a2**4*dy**2 &
-7.29d3*a1**2*a2**3*dy**2-7.29d3*a1**3*a2**2*dy**2 &
-2.43d3*a1**4*a2*dy**2+4.9d2*a1**4*a2**4*d2*dx**6 &
-4.655d3*a1**3*a2**4*dx**6-4.655d3*a1**4*a2**3*dx**6 &
-8.4d2*a1**4*a2**4*d2**2*dx**4+8.4d3*a1**3*a2**4*d2*dx**4 &
+8.4d3*a1**4*a2**3*d2*dx**4-3.57d3*a1**2*a2**4*dx**4 &
-7.14d3*a1**3*a2**3*dx**4-3.57d3*a1**4*a2**2*dx**4 &
+4.32d2*a1**4*a2**4*d2**3*dx**2-4.68d3*a1**3*a2**4*d2**2*dx**2 &
-4.68d3*a1**4*a2**3*d2**2*dx**2+5.22d3*a1**2*a2**4*d2*dx**2 &
+1.044d4*a1**3*a2**3*d2*dx**2+5.22d3*a1**4*a2**2*d2*dx**2 &
-2.43d3*a1*a2**4*dx**2-7.29d3*a1**2*a2**3*dx**2 &
-7.29d3*a1**3*a2**2*dx**2-2.43d3*a1**4*a2*dx**2 &
-6.4d1*a1**4*a2**4*d2**4+8.0d2*a1**3*a2**4*d2**3 &
+8.0d2*a1**4*a2**3*d2**3-1.92d3*a1**2*a2**4*d2**2 &
-3.84d3*a1**3*a2**3*d2**2-1.92d3*a1**4*a2**2*d2**2+2.28d3*a1*a2**4*d2 &
+6.84d3*a1**2*a2**3*d2+6.84d3*a1**3*a2**2*d2+2.28d3*a1**4*a2*d2 &
-7.8d2*a2**4-3.12d3*a1*a2**3-4.68d3*a1**2*a2**2-3.12d3*a1**3*a2 &
-7.8d2*a1**4)*dz
                case (2)
                  rlYlm_laplacian =5.57342901225845299d-1*E*a1**2*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-12)*(dy+dx)*(dy-1.0d0*dx)* &
(4.9d2*a1**4*a2**4*d2*dy**6-4.655d3*a1**3*a2**4*dy**6 &
-4.655d3*a1**4*a2**3*dy**6+1.47d3*a1**4*a2**4*d2*dx**2*dy**4 &
-1.3965d4*a1**3*a2**4*dx**2*dy**4-1.3965d4*a1**4*a2**3*dx**2*dy**4 &
-9.8d2*a1**4*a2**4*d2**2*dy**4+9.31d3*a1**3*a2**4*d2*dy**4 &
+9.31d3*a1**4*a2**3*d2*dy**4+1.47d3*a1**4*a2**4*d2*dx**4*dy**2 &
-1.3965d4*a1**3*a2**4*dx**4*dy**2-1.3965d4*a1**4*a2**3*dx**4*dy**2 &
-1.96d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.862d4*a1**3*a2**4*d2*dx**2*dy**2+1.862d4*a1**4*a2**3*d2*dx**2*dy**2 &
+5.92d2*a1**4*a2**4*d2**3*dy**2-5.48d3*a1**3*a2**4*d2**2*dy**2 &
-5.48d3*a1**4*a2**3*d2**2*dy**2-1.62d3*a1**2*a2**4*d2*dy**2 &
-3.24d3*a1**3*a2**3*d2*dy**2-1.62d3*a1**4*a2**2*d2*dy**2 &
+2.97d3*a1*a2**4*dy**2+8.91d3*a1**2*a2**3*dy**2 &
+8.91d3*a1**3*a2**2*dy**2+2.97d3*a1**4*a2*dy**2 &
+4.9d2*a1**4*a2**4*d2*dx**6-4.655d3*a1**3*a2**4*dx**6 &
-4.655d3*a1**4*a2**3*dx**6-9.8d2*a1**4*a2**4*d2**2*dx**4 &
+9.31d3*a1**3*a2**4*d2*dx**4+9.31d3*a1**4*a2**3*d2*dx**4 &
+5.92d2*a1**4*a2**4*d2**3*dx**2-5.48d3*a1**3*a2**4*d2**2*dx**2 &
-5.48d3*a1**4*a2**3*d2**2*dx**2-1.62d3*a1**2*a2**4*d2*dx**2 &
-3.24d3*a1**3*a2**3*d2*dx**2-1.62d3*a1**4*a2**2*d2*dx**2 &
+2.97d3*a1*a2**4*dx**2+8.91d3*a1**2*a2**3*dx**2 &
+8.91d3*a1**3*a2**2*dx**2+2.97d3*a1**4*a2*dx**2 &
-9.6d1*a1**4*a2**4*d2**4+7.2d2*a1**3*a2**4*d2**3 &
+7.2d2*a1**4*a2**3*d2**3+2.28d3*a1**2*a2**4*d2**2 &
+4.56d3*a1**3*a2**3*d2**2+2.28d3*a1**4*a2**2*d2**2-5.22d3*a1*a2**4*d2 &
-1.566d4*a1**2*a2**3*d2-1.566d4*a1**3*a2**2*d2-5.22d3*a1**4*a2*d2 &
+2.34d3*a2**4+9.36d3*a1*a2**3+1.404d4*a1**2*a2**2+9.36d3*a1**3*a2 &
+2.34d3*a1**4)
                case (3)
                  rlYlm_laplacian = &
-2.08538618333770325d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*(7.0d1*a1**3*a2**3*d2*dy**4 &
-6.65d2*a1**2*a2**3*dy**4-6.65d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.33d3*a1**2*a2**3*dx**2*dy**2 &
-1.33d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.6d2*a1**2*a2**3*d2*dy**2+6.6d2*a1**3*a2**2*d2*dy**2 &
+8.5d2*a1*a2**3*dy**2+1.7d3*a1**2*a2**2*dy**2+8.5d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-6.65d2*a1**2*a2**3*dx**4 &
-6.65d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.6d2*a1**2*a2**3*d2*dx**2+6.6d2*a1**3*a2**2*d2*dx**2 &
+8.5d2*a1*a2**3*dx**2+1.7d3*a1**2*a2**2*dx**2+8.5d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-4.0d1*a1**2*a2**3*d2**2 &
-4.0d1*a1**3*a2**2*d2**2-1.06d3*a1*a2**3*d2-2.12d3*a1**2*a2**2*d2 &
-1.06d3*a1**3*a2*d2+8.1d2*a2**3+2.43d3*a1*a2**2+2.43d3*a1**2*a2 &
+8.1d2*a1**3)*dz
                case (4)
                  rlYlm_laplacian =7.37295355815411407d-1*E*a1**3*a2**3*sqrt &
                  (a2+a1)*(a2+a1)**(-12)*(dy**2-2.0d0*dx*dy &
-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)* &
(7.0d1*a1**3*a2**3*d2*dy**4-6.65d2*a1**2*a2**3*dy**4 &
-6.65d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.33d3*a1**2*a2**3*dx**2*dy**2-1.33d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+5.2d2*a1**2*a2**3*d2*dy**2 &
+5.2d2*a1**3*a2**2*d2*dy**2+2.04d3*a1*a2**3*dy**2 &
+4.08d3*a1**2*a2**2*dy**2+2.04d3*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-6.65d2*a1**2*a2**3*dx**4 &
-6.65d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+5.2d2*a1**2*a2**3*d2*dx**2+5.2d2*a1**3*a2**2*d2*dx**2 &
+2.04d3*a1*a2**3*dx**2+4.08d3*a1**2*a2**2*dx**2+2.04d3*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3+4.0d1*a1**2*a2**3*d2**2 &
+4.0d1*a1**3*a2**2*d2**2-1.56d3*a1*a2**3*d2-3.12d3*a1**2*a2**2*d2 &
-1.56d3*a1**3*a2*d2-5.4d2*a2**3-1.62d3*a1*a2**2-1.62d3*a1**2*a2 &
-5.4d2*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (1)
          ! selection on l2: l1=4, m1=1
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=1, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-2.10187170614922326d0*E*a1*a2**5*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*(7.0d0*(dy**2+dx**2)-4.0d0*d2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=1, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-3.6405485860419361d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(1.4d1*a1**2*a2**2*d2*dy**2-9.1d1*a1*a2**2*dy**2-9.1d1*a1**2*a2*dy**2 &
+1.4d1*a1**2*a2**2*d2*dx**2-9.1d1*a1*a2**2*dx**2-9.1d1*a1**2*a2*dx**2 &
-8.0d0*a1**2*a2**2*d2**2+4.6d1*a1*a2**2*d2+4.6d1*a1**2*a2*d2 &
+3.3d1*a2**2+6.6d1*a1*a2+3.3d1*a1**2)*dz
                case (0)
                  rlYlm_laplacian =1.82027429302096805d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(2.8d1*a1**2*a2**2*d2*dy**4 &
-1.82d2*a1*a2**2*dy**4-1.82d2*a1**2*a2*dy**4 &
+5.6d1*a1**2*a2**2*d2*dx**2*dy**2-3.64d2*a1*a2**2*dx**2*dy**2 &
-3.64d2*a1**2*a2*dx**2*dy**2-4.4d1*a1**2*a2**2*d2**2*dy**2 &
+3.16d2*a1*a2**2*d2*dy**2+3.16d2*a1**2*a2*d2*dy**2-1.65d2*a2**2*dy**2 &
-3.3d2*a1*a2*dy**2-1.65d2*a1**2*dy**2+2.8d1*a1**2*a2**2*d2*dx**4 &
-1.82d2*a1*a2**2*dx**4-1.82d2*a1**2*a2*dx**4 &
-4.4d1*a1**2*a2**2*d2**2*dx**2+3.16d2*a1*a2**2*d2*dx**2 &
+3.16d2*a1**2*a2*d2*dx**2-1.65d2*a2**2*dx**2-3.3d2*a1*a2*dx**2 &
-1.65d2*a1**2*dx**2+1.6d1*a1**2*a2**2*d2**3-1.28d2*a1*a2**2*d2**2 &
-1.28d2*a1**2*a2*d2**2+1.32d2*a2**2*d2+2.64d2*a1*a2*d2+1.32d2*a1**2*d2 &
)
                case (1)
                  rlYlm_laplacian = &
-1.82027429302096805d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)* &
(2.8d1*a1**2*a2**2*d2*dx**2*dy**2-1.82d2*a1*a2**2*dx**2*dy**2 &
-1.82d2*a1**2*a2*dx**2*dy**2-1.4d1*a1*a2**2*d2*dy**2 &
-1.4d1*a1**2*a2*d2*dy**2+7.7d1*a2**2*dy**2+1.54d2*a1*a2*dy**2 &
+7.7d1*a1**2*dy**2+2.8d1*a1**2*a2**2*d2*dx**4-1.82d2*a1*a2**2*dx**4 &
-1.82d2*a1**2*a2*dx**4-1.6d1*a1**2*a2**2*d2**2*dx**2 &
+7.8d1*a1*a2**2*d2*dx**2+7.8d1*a1**2*a2*d2*dx**2+1.43d2*a2**2*dx**2 &
+2.86d2*a1*a2*dx**2+1.43d2*a1**2*dx**2+8.0d0*a1*a2**2*d2**2 &
+8.0d0*a1**2*a2*d2**2-4.4d1*a2**2*d2-8.8d1*a1*a2*d2-4.4d1*a1**2*d2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=1, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-4.07025705689025559d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(2.8d1*a1**3*a2**3*d2*dx**2*dy**2-2.1d2*a1**2*a2**3*dx**2*dy**2 &
-2.1d2*a1**3*a2**2*dx**2*dy**2-1.4d1*a1**2*a2**3*d2*dy**2 &
-1.4d1*a1**3*a2**2*d2*dy**2+9.1d1*a1*a2**3*dy**2 &
+1.82d2*a1**2*a2**2*dy**2+9.1d1*a1**3*a2*dy**2 &
+2.8d1*a1**3*a2**3*d2*dx**4-2.1d2*a1**2*a2**3*dx**4 &
-2.1d2*a1**3*a2**2*dx**4-1.6d1*a1**3*a2**3*d2**2*dx**2 &
+8.2d1*a1**2*a2**3*d2*dx**2+8.2d1*a1**3*a2**2*d2*dx**2 &
+2.47d2*a1*a2**3*dx**2+4.94d2*a1**2*a2**2*dx**2+2.47d2*a1**3*a2*dx**2 &
+8.0d0*a1**2*a2**3*d2**2+8.0d0*a1**3*a2**2*d2**2-4.6d1*a1*a2**3*d2 &
-9.2d1*a1**2*a2**2*d2-4.6d1*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =4.07025705689025559d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(2.8d1*a1**3*a2**3*d2*dy**4 &
-2.1d2*a1**2*a2**3*dy**4-2.1d2*a1**3*a2**2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**2*dy**2-4.2d2*a1**2*a2**3*dx**2*dy**2 &
-4.2d2*a1**3*a2**2*dx**2*dy**2-4.4d1*a1**3*a2**3*d2**2*dy**2 &
+3.48d2*a1**2*a2**3*d2*dy**2+3.48d2*a1**3*a2**2*d2*dy**2 &
-1.17d2*a1*a2**3*dy**2-2.34d2*a1**2*a2**2*dy**2-1.17d2*a1**3*a2*dy**2 &
+2.8d1*a1**3*a2**3*d2*dx**4-2.1d2*a1**2*a2**3*dx**4 &
-2.1d2*a1**3*a2**2*dx**4-4.4d1*a1**3*a2**3*d2**2*dx**2 &
+3.48d2*a1**2*a2**3*d2*dx**2+3.48d2*a1**3*a2**2*d2*dx**2 &
-1.17d2*a1*a2**3*dx**2-2.34d2*a1**2*a2**2*dx**2-1.17d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.32d2*a1**2*a2**3*d2**2 &
-1.32d2*a1**3*a2**2*d2**2+7.2d1*a1*a2**3*d2+1.44d2*a1**2*a2**2*d2 &
+7.2d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2+9.9d1*a1**2*a2 &
+3.3d1*a1**3)
                case (0)
                  rlYlm_laplacian =2.34996400746656297d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(4.2d1*a1**3*a2**3*d2*dy**4 &
-3.15d2*a1**2*a2**3*dy**4-3.15d2*a1**3*a2**2*dy**4 &
+8.4d1*a1**3*a2**3*d2*dx**2*dy**2-6.3d2*a1**2*a2**3*dx**2*dy**2 &
-6.3d2*a1**3*a2**2*dx**2*dy**2-5.2d1*a1**3*a2**3*d2**2*dy**2 &
+4.24d2*a1**2*a2**3*d2*dy**2+4.24d2*a1**3*a2**2*d2*dy**2 &
-2.21d2*a1*a2**3*dy**2-4.42d2*a1**2*a2**2*dy**2-2.21d2*a1**3*a2*dy**2 &
+4.2d1*a1**3*a2**3*d2*dx**4-3.15d2*a1**2*a2**3*dx**4 &
-3.15d2*a1**3*a2**2*dx**4-5.2d1*a1**3*a2**3*d2**2*dx**2 &
+4.24d2*a1**2*a2**3*d2*dx**2+4.24d2*a1**3*a2**2*d2*dx**2 &
-2.21d2*a1*a2**3*dx**2-4.42d2*a1**2*a2**2*dx**2-2.21d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-1.6d2*a1**2*a2**3*d2**2 &
-1.6d2*a1**3*a2**2*d2**2+2.96d2*a1*a2**3*d2+5.92d2*a1**2*a2**2*d2 &
+2.96d2*a1**3*a2*d2-1.98d2*a2**3-5.94d2*a1*a2**2-5.94d2*a1**2*a2 &
-1.98d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian =2.03512852844512779d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*(5.6d1*a1**3*a2**3*d2*dx**2*dy**4 &
-4.2d2*a1**2*a2**3*dx**2*dy**4-4.2d2*a1**3*a2**2*dx**2*dy**4 &
-2.8d1*a1**2*a2**3*d2*dy**4-2.8d1*a1**3*a2**2*d2*dy**4 &
+1.82d2*a1*a2**3*dy**4+3.64d2*a1**2*a2**2*dy**4+1.82d2*a1**3*a2*dy**4 &
+1.12d2*a1**3*a2**3*d2*dx**4*dy**2-8.4d2*a1**2*a2**3*dx**4*dy**2 &
-8.4d2*a1**3*a2**2*dx**4*dy**2-8.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+6.4d2*a1**2*a2**3*d2*dx**2*dy**2+6.4d2*a1**3*a2**2*d2*dx**2*dy**2 &
+1.3d2*a1*a2**3*dx**2*dy**2+2.6d2*a1**2*a2**2*dx**2*dy**2 &
+1.3d2*a1**3*a2*dx**2*dy**2+4.4d1*a1**2*a2**3*d2**2*dy**2 &
+4.4d1*a1**3*a2**2*d2**2*dy**2-3.16d2*a1*a2**3*d2*dy**2 &
-6.32d2*a1**2*a2**2*d2*dy**2-3.16d2*a1**3*a2*d2*dy**2 &
+1.65d2*a2**3*dy**2+4.95d2*a1*a2**2*dy**2+4.95d2*a1**2*a2*dy**2 &
+1.65d2*a1**3*dy**2+5.6d1*a1**3*a2**3*d2*dx**6-4.2d2*a1**2*a2**3*dx**6 &
-4.2d2*a1**3*a2**2*dx**6-8.8d1*a1**3*a2**3*d2**2*dx**4 &
+6.68d2*a1**2*a2**3*d2*dx**4+6.68d2*a1**3*a2**2*d2*dx**4 &
-5.2d1*a1*a2**3*dx**4-1.04d2*a1**2*a2**2*dx**4-5.2d1*a1**3*a2*dx**4 &
+3.2d1*a1**3*a2**3*d2**3*dx**2-2.2d2*a1**2*a2**3*d2**2*dx**2 &
-2.2d2*a1**3*a2**2*d2**2*dx**2-1.72d2*a1*a2**3*d2*dx**2 &
-3.44d2*a1**2*a2**2*d2*dx**2-1.72d2*a1**3*a2*d2*dx**2 &
+2.31d2*a2**3*dx**2+6.93d2*a1*a2**2*dx**2+6.93d2*a1**2*a2*dx**2 &
+2.31d2*a1**3*dx**2-1.6d1*a1**2*a2**3*d2**3-1.6d1*a1**3*a2**2*d2**3 &
+1.28d2*a1*a2**3*d2**2+2.56d2*a1**2*a2**2*d2**2+1.28d2*a1**3*a2*d2**2 &
-1.32d2*a2**3*d2-3.96d2*a1*a2**2*d2-3.96d2*a1**2*a2*d2-1.32d2*a1**3*d2 &
)
                case (2)
                  rlYlm_laplacian =4.07025705689025559d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.05d2*a1**2*a2**3*dy**4-1.05d2*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**2+6.2d1*a1**2*a2**3*d2*dy**2 &
+6.2d1*a1**3*a2**2*d2*dy**2-1.3d1*a1*a2**3*dy**2 &
-2.6d1*a1**2*a2**2*dy**2-1.3d1*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.05d2*a1**2*a2**3*dx**4 &
+1.05d2*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-3.4d1*a1**2*a2**3*d2*dx**2-3.4d1*a1**3*a2**2*d2*dx**2 &
-1.69d2*a1*a2**3*dx**2-3.38d2*a1**2*a2**2*dx**2-1.69d2*a1**3*a2*dx**2 &
-8.0d0*a1**2*a2**3*d2**2-8.0d0*a1**3*a2**2*d2**2+4.6d1*a1*a2**3*d2 &
+9.2d1*a1**2*a2**2*d2+4.6d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=1, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =4.39638009359507944d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.19d2*a1**2*a2**3*dy**4-1.19d2*a1**3*a2**2*dy**4 &
-2.8d1*a1**3*a2**3*d2*dx**2*dy**2+2.38d2*a1**2*a2**3*dx**2*dy**2 &
+2.38d2*a1**3*a2**2*dx**2*dy**2-8.0d0*a1**3*a2**3*d2**2*dy**2 &
+9.2d1*a1**2*a2**3*d2*dy**2+9.2d1*a1**3*a2**2*d2*dy**2 &
-1.8d2*a1*a2**3*dy**2-3.6d2*a1**2*a2**2*dy**2-1.8d2*a1**3*a2*dy**2 &
-4.2d1*a1**3*a2**3*d2*dx**4+3.57d2*a1**2*a2**3*dx**4 &
+3.57d2*a1**3*a2**2*dx**4+2.4d1*a1**3*a2**3*d2**2*dx**2 &
-1.08d2*a1**2*a2**3*d2*dx**2-1.08d2*a1**3*a2**2*d2*dx**2 &
-7.2d2*a1*a2**3*dx**2-1.44d3*a1**2*a2**2*dx**2-7.2d2*a1**3*a2*dx**2 &
-2.4d1*a1**2*a2**3*d2**2-2.4d1*a1**3*a2**2*d2**2+1.44d2*a1*a2**3*d2 &
+2.88d2*a1**2*a2**2*d2+1.44d2*a1**3*a2*d2+2.34d2*a2**3+7.02d2*a1*a2**2 &
+7.02d2*a1**2*a2+2.34d2*a1**3)*dz
                case (-2)
                  rlYlm_laplacian =5.3844439723186478d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(5.6d1*a1**4*a2**4*d2*dx**2*dy**4 &
-4.76d2*a1**3*a2**4*dx**2*dy**4-4.76d2*a1**4*a2**3*dx**2*dy**4 &
-2.8d1*a1**3*a2**4*d2*dy**4-2.8d1*a1**4*a2**3*d2*dy**4 &
+2.1d2*a1**2*a2**4*dy**4+4.2d2*a1**3*a2**3*dy**4 &
+2.1d2*a1**4*a2**2*dy**4+1.12d2*a1**4*a2**4*d2*dx**4*dy**2 &
-9.52d2*a1**3*a2**4*dx**4*dy**2-9.52d2*a1**4*a2**3*dx**4*dy**2 &
-8.8d1*a1**4*a2**4*d2**2*dx**2*dy**2+7.04d2*a1**3*a2**4*d2*dx**2*dy**2 &
+7.04d2*a1**4*a2**3*d2*dx**2*dy**2+3.3d2*a1**2*a2**4*dx**2*dy**2 &
+6.6d2*a1**3*a2**3*dx**2*dy**2+3.3d2*a1**4*a2**2*dx**2*dy**2 &
+4.4d1*a1**3*a2**4*d2**2*dy**2+4.4d1*a1**4*a2**3*d2**2*dy**2 &
-3.48d2*a1**2*a2**4*d2*dy**2-6.96d2*a1**3*a2**3*d2*dy**2 &
-3.48d2*a1**4*a2**2*d2*dy**2+1.17d2*a1*a2**4*dy**2 &
+3.51d2*a1**2*a2**3*dy**2+3.51d2*a1**3*a2**2*dy**2 &
+1.17d2*a1**4*a2*dy**2+5.6d1*a1**4*a2**4*d2*dx**6 &
-4.76d2*a1**3*a2**4*dx**6-4.76d2*a1**4*a2**3*dx**6 &
-8.8d1*a1**4*a2**4*d2**2*dx**4+7.32d2*a1**3*a2**4*d2*dx**4 &
+7.32d2*a1**4*a2**3*d2*dx**4+1.2d2*a1**2*a2**4*dx**4 &
+2.4d2*a1**3*a2**3*dx**4+1.2d2*a1**4*a2**2*dx**4 &
+3.2d1*a1**4*a2**4*d2**3*dx**2-2.28d2*a1**3*a2**4*d2**2*dx**2 &
-2.28d2*a1**4*a2**3*d2**2*dx**2-3.72d2*a1**2*a2**4*d2*dx**2 &
-7.44d2*a1**3*a2**3*d2*dx**2-3.72d2*a1**4*a2**2*d2*dx**2 &
+2.73d2*a1*a2**4*dx**2+8.19d2*a1**2*a2**3*dx**2 &
+8.19d2*a1**3*a2**2*dx**2+2.73d2*a1**4*a2*dx**2 &
-1.6d1*a1**3*a2**4*d2**3-1.6d1*a1**4*a2**3*d2**3 &
+1.32d2*a1**2*a2**4*d2**2+2.64d2*a1**3*a2**3*d2**2 &
+1.32d2*a1**4*a2**2*d2**2-7.2d1*a1*a2**4*d2-2.16d2*a1**2*a2**3*d2 &
-2.16d2*a1**3*a2**2*d2-7.2d1*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (-1)
                  rlYlm_laplacian =3.40542137721830949d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(7.0d1*a1**3*a2**3*d2*dy**4 &
-5.95d2*a1**2*a2**3*dy**4-5.95d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.19d3*a1**2*a2**3*dx**2*dy**2 &
-1.19d3*a1**3*a2**2*dx**2*dy**2-9.6d1*a1**3*a2**3*d2**2*dy**2 &
+8.8d2*a1**2*a2**3*d2*dy**2+8.8d2*a1**3*a2**2*d2*dy**2 &
-4.8d2*a1*a2**3*dy**2-9.6d2*a1**2*a2**2*dy**2-4.8d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-9.6d1*a1**3*a2**3*d2**2*dx**2 &
+8.8d2*a1**2*a2**3*d2*dx**2+8.8d2*a1**3*a2**2*d2*dx**2 &
-4.8d2*a1*a2**3*dx**2-9.6d2*a1**2*a2**2*dx**2-4.8d2*a1**3*a2*dx**2 &
+3.2d1*a1**3*a2**3*d2**3-3.36d2*a1**2*a2**3*d2**2 &
-3.36d2*a1**3*a2**2*d2**2+5.28d2*a1*a2**3*d2+1.056d3*a1**2*a2**2*d2 &
+5.28d2*a1**3*a2*d2-3.12d2*a2**3-9.36d2*a1*a2**2-9.36d2*a1**2*a2 &
-3.12d2*a1**3)*dz
                case (0)
                  rlYlm_laplacian = &
-1.39025745555846884d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(1.4d2*a1**4*a2**4*d2*dy**6-1.19d3*a1**3*a2**4*dy**6 &
-1.19d3*a1**4*a2**3*dy**6+4.2d2*a1**4*a2**4*d2*dx**2*dy**4 &
-3.57d3*a1**3*a2**4*dx**2*dy**4-3.57d3*a1**4*a2**3*dx**2*dy**4 &
-2.76d2*a1**4*a2**4*d2**2*dy**4+2.46d3*a1**3*a2**4*d2*dy**4 &
+2.46d3*a1**4*a2**3*d2*dy**4-8.55d2*a1**2*a2**4*dy**4 &
-1.71d3*a1**3*a2**3*dy**4-8.55d2*a1**4*a2**2*dy**4 &
+4.2d2*a1**4*a2**4*d2*dx**4*dy**2-3.57d3*a1**3*a2**4*dx**4*dy**2 &
-3.57d3*a1**4*a2**3*dx**4*dy**2-5.52d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+4.92d3*a1**3*a2**4*d2*dx**2*dy**2+4.92d3*a1**4*a2**3*d2*dx**2*dy**2 &
-1.71d3*a1**2*a2**4*dx**2*dy**2-3.42d3*a1**3*a2**3*dx**2*dy**2 &
-1.71d3*a1**4*a2**2*dx**2*dy**2+1.68d2*a1**4*a2**4*d2**3*dy**2 &
-1.62d3*a1**3*a2**4*d2**2*dy**2-1.62d3*a1**4*a2**3*d2**2*dy**2 &
+1.53d3*a1**2*a2**4*d2*dy**2+3.06d3*a1**3*a2**3*d2*dy**2 &
+1.53d3*a1**4*a2**2*d2*dy**2-5.85d2*a1*a2**4*dy**2 &
-1.755d3*a1**2*a2**3*dy**2-1.755d3*a1**3*a2**2*dy**2 &
-5.85d2*a1**4*a2*dy**2+1.4d2*a1**4*a2**4*d2*dx**6 &
-1.19d3*a1**3*a2**4*dx**6-1.19d3*a1**4*a2**3*dx**6 &
-2.76d2*a1**4*a2**4*d2**2*dx**4+2.46d3*a1**3*a2**4*d2*dx**4 &
+2.46d3*a1**4*a2**3*d2*dx**4-8.55d2*a1**2*a2**4*dx**4 &
-1.71d3*a1**3*a2**3*dx**4-8.55d2*a1**4*a2**2*dx**4 &
+1.68d2*a1**4*a2**4*d2**3*dx**2-1.62d3*a1**3*a2**4*d2**2*dx**2 &
-1.62d3*a1**4*a2**3*d2**2*dx**2+1.53d3*a1**2*a2**4*d2*dx**2 &
+3.06d3*a1**3*a2**3*d2*dx**2+1.53d3*a1**4*a2**2*d2*dx**2 &
-5.85d2*a1*a2**4*dx**2-1.755d3*a1**2*a2**3*dx**2 &
-1.755d3*a1**3*a2**2*dx**2-5.85d2*a1**4*a2*dx**2 &
-3.2d1*a1**4*a2**4*d2**4+3.68d2*a1**3*a2**4*d2**3 &
+3.68d2*a1**4*a2**3*d2**3-8.64d2*a1**2*a2**4*d2**2 &
-1.728d3*a1**3*a2**3*d2**2-8.64d2*a1**4*a2**2*d2**2+9.96d2*a1*a2**4*d2 &
+2.988d3*a1**2*a2**3*d2+2.988d3*a1**3*a2**2*d2+9.96d2*a1**4*a2*d2 &
-3.3d2*a2**4-1.32d3*a1*a2**3-1.98d3*a1**2*a2**2-1.32d3*a1**3*a2 &
-3.3d2*a1**4)
                case (1)
                  rlYlm_laplacian =1.70271068860915474d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(1.4d2*a1**4*a2**4*d2*dx**2*dy**4 &
-1.19d3*a1**3*a2**4*dx**2*dy**4-1.19d3*a1**4*a2**3*dx**2*dy**4 &
-7.0d1*a1**3*a2**4*d2*dy**4-7.0d1*a1**4*a2**3*d2*dy**4 &
+5.25d2*a1**2*a2**4*dy**4+1.05d3*a1**3*a2**3*dy**4 &
+5.25d2*a1**4*a2**2*dy**4+2.8d2*a1**4*a2**4*d2*dx**4*dy**2 &
-2.38d3*a1**3*a2**4*dx**4*dy**2-2.38d3*a1**4*a2**3*dx**4*dy**2 &
-1.92d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.62d3*a1**3*a2**4*d2*dx**2*dy**2+1.62d3*a1**4*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**4*dx**2*dy**2+1.8d2*a1**3*a2**3*dx**2*dy**2 &
+9.0d1*a1**4*a2**2*dx**2*dy**2+9.6d1*a1**3*a2**4*d2**2*dy**2 &
+9.6d1*a1**4*a2**3*d2**2*dy**2-8.28d2*a1**2*a2**4*d2*dy**2 &
-1.656d3*a1**3*a2**3*d2*dy**2-8.28d2*a1**4*a2**2*d2*dy**2 &
+7.02d2*a1*a2**4*dy**2+2.106d3*a1**2*a2**3*dy**2 &
+2.106d3*a1**3*a2**2*dy**2+7.02d2*a1**4*a2*dy**2 &
+1.4d2*a1**4*a2**4*d2*dx**6-1.19d3*a1**3*a2**4*dx**6 &
-1.19d3*a1**4*a2**3*dx**6-1.92d2*a1**4*a2**4*d2**2*dx**4 &
+1.69d3*a1**3*a2**4*d2*dx**4+1.69d3*a1**4*a2**3*d2*dx**4 &
-4.35d2*a1**2*a2**4*dx**4-8.7d2*a1**3*a2**3*dx**4 &
-4.35d2*a1**4*a2**2*dx**4+6.4d1*a1**4*a2**4*d2**3*dx**2 &
-5.76d2*a1**3*a2**4*d2**2*dx**2-5.76d2*a1**4*a2**3*d2**2*dx**2 &
+2.28d2*a1**2*a2**4*d2*dx**2+4.56d2*a1**3*a2**3*d2*dx**2 &
+2.28d2*a1**4*a2**2*d2*dx**2+7.8d1*a1*a2**4*dx**2 &
+2.34d2*a1**2*a2**3*dx**2+2.34d2*a1**3*a2**2*dx**2 &
+7.8d1*a1**4*a2*dx**2-3.2d1*a1**3*a2**4*d2**3-3.2d1*a1**4*a2**3*d2**3 &
+3.36d2*a1**2*a2**4*d2**2+6.72d2*a1**3*a2**3*d2**2 &
+3.36d2*a1**4*a2**2*d2**2-6.84d2*a1*a2**4*d2-2.052d3*a1**2*a2**3*d2 &
-2.052d3*a1**3*a2**2*d2-6.84d2*a1**4*a2*d2+3.3d2*a2**4+1.32d3*a1*a2**3 &
+1.98d3*a1**2*a2**2+1.32d3*a1**3*a2+3.3d2*a1**4)*dz
                case (2)
                  rlYlm_laplacian = &
-5.3844439723186478d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(2.8d1*a1**4*a2**4*d2*dy**6-2.38d2*a1**3*a2**4*dy**6 &
-2.38d2*a1**4*a2**3*dy**6+2.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
-2.38d2*a1**3*a2**4*dx**2*dy**4-2.38d2*a1**4*a2**3*dx**2*dy**4 &
-4.4d1*a1**4*a2**4*d2**2*dy**4+4.08d2*a1**3*a2**4*d2*dy**4 &
+4.08d2*a1**4*a2**3*d2*dy**4-2.55d2*a1**2*a2**4*dy**4 &
-5.1d2*a1**3*a2**3*dy**4-2.55d2*a1**4*a2**2*dy**4 &
-2.8d1*a1**4*a2**4*d2*dx**4*dy**2+2.38d2*a1**3*a2**4*dx**4*dy**2 &
+2.38d2*a1**4*a2**3*dx**4*dy**2+5.6d1*a1**3*a2**4*d2*dx**2*dy**2 &
+5.6d1*a1**4*a2**3*d2*dx**2*dy**2-4.2d2*a1**2*a2**4*dx**2*dy**2 &
-8.4d2*a1**3*a2**3*dx**2*dy**2-4.2d2*a1**4*a2**2*dx**2*dy**2 &
+1.6d1*a1**4*a2**4*d2**3*dy**2-1.8d2*a1**3*a2**4*d2**2*dy**2 &
-1.8d2*a1**4*a2**3*d2**2*dy**2+3.36d2*a1**2*a2**4*d2*dy**2 &
+6.72d2*a1**3*a2**3*d2*dy**2+3.36d2*a1**4*a2**2*d2*dy**2 &
-3.9d1*a1*a2**4*dy**2-1.17d2*a1**2*a2**3*dy**2 &
-1.17d2*a1**3*a2**2*dy**2-3.9d1*a1**4*a2*dy**2 &
-2.8d1*a1**4*a2**4*d2*dx**6+2.38d2*a1**3*a2**4*dx**6 &
+2.38d2*a1**4*a2**3*dx**6+4.4d1*a1**4*a2**4*d2**2*dx**4 &
-3.52d2*a1**3*a2**4*d2*dx**4-3.52d2*a1**4*a2**3*d2*dx**4 &
-1.65d2*a1**2*a2**4*dx**4-3.3d2*a1**3*a2**3*dx**4 &
-1.65d2*a1**4*a2**2*dx**4-1.6d1*a1**4*a2**4*d2**3*dx**2 &
+9.2d1*a1**3*a2**4*d2**2*dx**2+9.2d1*a1**4*a2**3*d2**2*dx**2 &
+3.6d2*a1**2*a2**4*d2*dx**2+7.2d2*a1**3*a2**3*d2*dx**2 &
+3.6d2*a1**4*a2**2*d2*dx**2-1.95d2*a1*a2**4*dx**2 &
-5.85d2*a1**2*a2**3*dx**2-5.85d2*a1**3*a2**2*dx**2 &
-1.95d2*a1**4*a2*dx**2+1.6d1*a1**3*a2**4*d2**3+1.6d1*a1**4*a2**3*d2**3 &
-1.32d2*a1**2*a2**4*d2**2-2.64d2*a1**3*a2**3*d2**2 &
-1.32d2*a1**4*a2**2*d2**2+7.2d1*a1*a2**4*d2+2.16d2*a1**2*a2**3*d2 &
+2.16d2*a1**3*a2**2*d2+7.2d1*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case (3)
                  rlYlm_laplacian =2.19819004679753972d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*(8.4d1*a1**3*a2**3*d2*dx**2*dy**4 &
-7.14d2*a1**2*a2**3*dx**2*dy**4-7.14d2*a1**3*a2**2*dx**2*dy**4 &
-4.2d1*a1**2*a2**3*d2*dy**4-4.2d1*a1**3*a2**2*d2*dy**4 &
+3.15d2*a1*a2**3*dy**4+6.3d2*a1**2*a2**2*dy**4+3.15d2*a1**3*a2*dy**4 &
+5.6d1*a1**3*a2**3*d2*dx**4*dy**2-4.76d2*a1**2*a2**3*dx**4*dy**2 &
-4.76d2*a1**3*a2**2*dx**4*dy**2-4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+3.0d2*a1**2*a2**3*d2*dx**2*dy**2+3.0d2*a1**3*a2**2*d2*dx**2*dy**2 &
+8.1d2*a1*a2**3*dx**2*dy**2+1.62d3*a1**2*a2**2*dx**2*dy**2 &
+8.1d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.44d2*a1*a2**3*d2*dy**2 &
-2.88d2*a1**2*a2**2*d2*dy**2-1.44d2*a1**3*a2*d2*dy**2 &
-2.34d2*a2**3*dy**2-7.02d2*a1*a2**2*dy**2-7.02d2*a1**2*a2*dy**2 &
-2.34d2*a1**3*dy**2-2.8d1*a1**3*a2**3*d2*dx**6 &
+2.38d2*a1**2*a2**3*dx**6+2.38d2*a1**3*a2**2*dx**6 &
+1.6d1*a1**3*a2**3*d2**2*dx**4-5.8d1*a1**2*a2**3*d2*dx**4 &
-5.8d1*a1**3*a2**2*d2*dx**4-5.85d2*a1*a2**3*dx**4 &
-1.17d3*a1**2*a2**2*dx**4-5.85d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.44d2*a1*a2**3*d2*dx**2+2.88d2*a1**2*a2**2*d2*dx**2 &
+1.44d2*a1**3*a2*d2*dx**2+2.34d2*a2**3*dx**2+7.02d2*a1*a2**2*dx**2 &
+7.02d2*a1**2*a2*dx**2+2.34d2*a1**3*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=1, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =9.3261305305638875d0*E*a1**3*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-12)*dy*(2.8d1*a1**3*a2**3*d2*dx**2*dy**4 &
-2.66d2*a1**2*a2**3*dx**2*dy**4-2.66d2*a1**3*a2**2*dx**2*dy**4 &
-1.4d1*a1**2*a2**3*d2*dy**4-1.4d1*a1**3*a2**2*d2*dy**4 &
+1.19d2*a1*a2**3*dy**4+2.38d2*a1**2*a2**2*dy**4+1.19d2*a1**3*a2*dy**4 &
-1.6d1*a1**3*a2**3*d2**2*dx**2*dy**2+1.32d2*a1**2*a2**3*d2*dx**2*dy**2 &
+1.32d2*a1**3*a2**2*d2*dx**2*dy**2+1.7d2*a1*a2**3*dx**2*dy**2 &
+3.4d2*a1**2*a2**2*dx**2*dy**2+1.7d2*a1**3*a2*dx**2*dy**2 &
+8.0d0*a1**2*a2**3*d2**2*dy**2+8.0d0*a1**3*a2**2*d2**2*dy**2 &
-5.0d1*a1*a2**3*d2*dy**2-1.0d2*a1**2*a2**2*d2*dy**2 &
-5.0d1*a1**3*a2*d2*dy**2-1.35d2*a2**3*dy**2-4.05d2*a1*a2**2*dy**2 &
-4.05d2*a1**2*a2*dy**2-1.35d2*a1**3*dy**2-2.8d1*a1**3*a2**3*d2*dx**6 &
+2.66d2*a1**2*a2**3*dx**6+2.66d2*a1**3*a2**2*dx**6 &
+1.6d1*a1**3*a2**3*d2**2*dx**4-6.2d1*a1**2*a2**3*d2*dx**4 &
-6.2d1*a1**3*a2**2*d2*dx**4-7.65d2*a1*a2**3*dx**4 &
-1.53d3*a1**2*a2**2*dx**4-7.65d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.5d2*a1*a2**3*d2*dx**2+3.0d2*a1**2*a2**2*d2*dx**2 &
+1.5d2*a1**3*a2*d2*dx**2+4.05d2*a2**3*dx**2+1.215d3*a1*a2**2*dx**2 &
+1.215d3*a1**2*a2*dx**2+4.05d2*a1**3*dx**2)*dz
                case (-3)
                  rlYlm_laplacian = &
-6.59457014039261917d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(2.8d1*a1**4*a2**4*d2*dy**6-2.66d2*a1**3*a2**4*dy**6 &
-2.66d2*a1**4*a2**3*dy**6-2.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
+2.66d2*a1**3*a2**4*dx**2*dy**4+2.66d2*a1**4*a2**3*dx**2*dy**4 &
-4.4d1*a1**4*a2**4*d2**2*dy**4+4.96d2*a1**3*a2**4*d2*dy**4 &
+4.96d2*a1**4*a2**3*d2*dy**4-6.63d2*a1**2*a2**4*dy**4 &
-1.326d3*a1**3*a2**3*dy**4-6.63d2*a1**4*a2**2*dy**4 &
-1.4d2*a1**4*a2**4*d2*dx**4*dy**2+1.33d3*a1**3*a2**4*dx**4*dy**2 &
+1.33d3*a1**4*a2**3*dx**4*dy**2+8.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-6.56d2*a1**3*a2**4*d2*dx**2*dy**2-6.56d2*a1**4*a2**3*d2*dx**2*dy**2 &
-1.53d3*a1**2*a2**4*dx**2*dy**2-3.06d3*a1**3*a2**3*dx**2*dy**2 &
-1.53d3*a1**4*a2**2*dx**2*dy**2+1.6d1*a1**4*a2**4*d2**3*dy**2 &
-2.72d2*a1**3*a2**4*d2**2*dy**2-2.72d2*a1**4*a2**3*d2**2*dy**2 &
+1.02d3*a1**2*a2**4*d2*dy**2+2.04d3*a1**3*a2**3*d2*dy**2 &
+1.02d3*a1**4*a2**2*d2*dy**2-8.4d1*a1**4*a2**4*d2*dx**6 &
+7.98d2*a1**3*a2**4*dx**6+7.98d2*a1**4*a2**3*dx**6 &
+1.32d2*a1**4*a2**4*d2**2*dx**4-1.152d3*a1**3*a2**4*d2*dx**4 &
-1.152d3*a1**4*a2**3*d2*dx**4-8.67d2*a1**2*a2**4*dx**4 &
-1.734d3*a1**3*a2**3*dx**4-8.67d2*a1**4*a2**2*dx**4 &
-4.8d1*a1**4*a2**4*d2**3*dx**2+2.88d2*a1**3*a2**4*d2**2*dx**2 &
+2.88d2*a1**4*a2**3*d2**2*dx**2+1.5d3*a1**2*a2**4*d2*dx**2 &
+3.0d3*a1**3*a2**3*d2*dx**2+1.5d3*a1**4*a2**2*d2*dx**2 &
-5.4d2*a1*a2**4*dx**2-1.62d3*a1**2*a2**3*dx**2 &
-1.62d3*a1**3*a2**2*dx**2-5.4d2*a1**4*a2*dx**2+4.8d1*a1**3*a2**4*d2**3 &
+4.8d1*a1**4*a2**3*d2**3-4.08d2*a1**2*a2**4*d2**2 &
-8.16d2*a1**3*a2**3*d2**2-4.08d2*a1**4*a2**2*d2**2-3.6d1*a1*a2**4*d2 &
-1.08d2*a1**2*a2**3*d2-1.08d2*a1**3*a2**2*d2-3.6d1*a1**4*a2*d2 &
+2.34d2*a2**4+9.36d2*a1*a2**3+1.404d3*a1**2*a2**2+9.36d2*a1**3*a2 &
+2.34d2*a1**4)
                case (-2)
                  rlYlm_laplacian =3.52494601119984446d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dy*(1.96d2*a1**4*a2**4*d2*dx**2*dy**4 &
-1.862d3*a1**3*a2**4*dx**2*dy**4-1.862d3*a1**4*a2**3*dx**2*dy**4 &
-9.8d1*a1**3*a2**4*d2*dy**4-9.8d1*a1**4*a2**3*d2*dy**4 &
+8.33d2*a1**2*a2**4*dy**4+1.666d3*a1**3*a2**3*dy**4 &
+8.33d2*a1**4*a2**2*dy**4+3.92d2*a1**4*a2**4*d2*dx**4*dy**2 &
-3.724d3*a1**3*a2**4*dx**4*dy**2-3.724d3*a1**4*a2**3*dx**4*dy**2 &
-2.8d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+2.604d3*a1**3*a2**4*d2*dx**2*dy**2+2.604d3*a1**4*a2**3*d2*dx**2*dy**2 &
+4.76d2*a1**2*a2**4*dx**2*dy**2+9.52d2*a1**3*a2**3*dx**2*dy**2 &
+4.76d2*a1**4*a2**2*dx**2*dy**2+1.4d2*a1**3*a2**4*d2**2*dy**2 &
+1.4d2*a1**4*a2**3*d2**2*dy**2-1.316d3*a1**2*a2**4*d2*dy**2 &
-2.632d3*a1**3*a2**3*d2*dy**2-1.316d3*a1**4*a2**2*d2*dy**2 &
+9.45d2*a1*a2**4*dy**2+2.835d3*a1**2*a2**3*dy**2 &
+2.835d3*a1**3*a2**2*dy**2+9.45d2*a1**4*a2*dy**2 &
+1.96d2*a1**4*a2**4*d2*dx**6-1.862d3*a1**3*a2**4*dx**6 &
-1.862d3*a1**4*a2**3*dx**6-2.8d2*a1**4*a2**4*d2**2*dx**4 &
+2.702d3*a1**3*a2**4*d2*dx**4+2.702d3*a1**4*a2**3*d2*dx**4 &
-3.57d2*a1**2*a2**4*dx**4-7.14d2*a1**3*a2**3*dx**4 &
-3.57d2*a1**4*a2**2*dx**4+9.6d1*a1**4*a2**4*d2**3*dx**2 &
-9.0d2*a1**3*a2**4*d2**2*dx**2-9.0d2*a1**4*a2**3*d2**2*dx**2 &
-1.56d2*a1**2*a2**4*d2*dx**2-3.12d2*a1**3*a2**3*d2*dx**2 &
-1.56d2*a1**4*a2**2*d2*dx**2+4.05d2*a1*a2**4*dx**2 &
+1.215d3*a1**2*a2**3*dx**2+1.215d3*a1**3*a2**2*dx**2 &
+4.05d2*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+5.16d2*a1**2*a2**4*d2**2+1.032d3*a1**3*a2**3*d2**2 &
+5.16d2*a1**4*a2**2*d2**2-8.64d2*a1*a2**4*d2-2.592d3*a1**2*a2**3*d2 &
-2.592d3*a1**3*a2**2*d2-8.64d2*a1**4*a2*d2+3.51d2*a2**4 &
+1.404d3*a1*a2**3+2.106d3*a1**2*a2**2+1.404d3*a1**3*a2+3.51d2*a1**4 &
)*dz
                case (-1)
                  rlYlm_laplacian = &
-2.49251322783588191d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(1.96d2*a1**4*a2**4*d2*dy**6-1.862d3*a1**3*a2**4*dy**6 &
-1.862d3*a1**4*a2**3*dy**6+5.88d2*a1**4*a2**4*d2*dx**2*dy**4 &
-5.586d3*a1**3*a2**4*dx**2*dy**4-5.586d3*a1**4*a2**3*dx**2*dy**4 &
-4.2d2*a1**4*a2**4*d2**2*dy**4+4.2d3*a1**3*a2**4*d2*dy**4 &
+4.2d3*a1**4*a2**3*d2*dy**4-1.785d3*a1**2*a2**4*dy**4 &
-3.57d3*a1**3*a2**3*dy**4-1.785d3*a1**4*a2**2*dy**4 &
+5.88d2*a1**4*a2**4*d2*dx**4*dy**2-5.586d3*a1**3*a2**4*dx**4*dy**2 &
-5.586d3*a1**4*a2**3*dx**4*dy**2-8.4d2*a1**4*a2**4*d2**2*dx**2*dy**2 &
+8.4d3*a1**3*a2**4*d2*dx**2*dy**2+8.4d3*a1**4*a2**3*d2*dx**2*dy**2 &
-3.57d3*a1**2*a2**4*dx**2*dy**2-7.14d3*a1**3*a2**3*dx**2*dy**2 &
-3.57d3*a1**4*a2**2*dx**2*dy**2+2.88d2*a1**4*a2**4*d2**3*dy**2 &
-3.12d3*a1**3*a2**4*d2**2*dy**2-3.12d3*a1**4*a2**3*d2**2*dy**2 &
+3.48d3*a1**2*a2**4*d2*dy**2+6.96d3*a1**3*a2**3*d2*dy**2 &
+3.48d3*a1**4*a2**2*d2*dy**2-1.62d3*a1*a2**4*dy**2 &
-4.86d3*a1**2*a2**3*dy**2-4.86d3*a1**3*a2**2*dy**2 &
-1.62d3*a1**4*a2*dy**2+1.96d2*a1**4*a2**4*d2*dx**6 &
-1.862d3*a1**3*a2**4*dx**6-1.862d3*a1**4*a2**3*dx**6 &
-4.2d2*a1**4*a2**4*d2**2*dx**4+4.2d3*a1**3*a2**4*d2*dx**4 &
+4.2d3*a1**4*a2**3*d2*dx**4-1.785d3*a1**2*a2**4*dx**4 &
-3.57d3*a1**3*a2**3*dx**4-1.785d3*a1**4*a2**2*dx**4 &
+2.88d2*a1**4*a2**4*d2**3*dx**2-3.12d3*a1**3*a2**4*d2**2*dx**2 &
-3.12d3*a1**4*a2**3*d2**2*dx**2+3.48d3*a1**2*a2**4*d2*dx**2 &
+6.96d3*a1**3*a2**3*d2*dx**2+3.48d3*a1**4*a2**2*d2*dx**2 &
-1.62d3*a1*a2**4*dx**2-4.86d3*a1**2*a2**3*dx**2 &
-4.86d3*a1**3*a2**2*dx**2-1.62d3*a1**4*a2*dx**2 &
-6.4d1*a1**4*a2**4*d2**4+8.0d2*a1**3*a2**4*d2**3 &
+8.0d2*a1**4*a2**3*d2**3-1.92d3*a1**2*a2**4*d2**2 &
-3.84d3*a1**3*a2**3*d2**2-1.92d3*a1**4*a2**2*d2**2+2.28d3*a1*a2**4*d2 &
+6.84d3*a1**2*a2**3*d2+6.84d3*a1**3*a2**2*d2+2.28d3*a1**4*a2*d2 &
-7.8d2*a2**4-3.12d3*a1*a2**3-4.68d3*a1**2*a2**2-3.12d3*a1**3*a2 &
-7.8d2*a1**4)
                case (0)
                  rlYlm_laplacian =-7.88201889805958723d-1*E*a1**2*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-12)*dx* &
(4.9d2*a1**4*a2**4*d2*dy**6-4.655d3*a1**3*a2**4*dy**6 &
-4.655d3*a1**4*a2**3*dy**6+1.47d3*a1**4*a2**4*d2*dx**2*dy**4 &
-1.3965d4*a1**3*a2**4*dx**2*dy**4-1.3965d4*a1**4*a2**3*dx**2*dy**4 &
-8.4d2*a1**4*a2**4*d2**2*dy**4+8.4d3*a1**3*a2**4*d2*dy**4 &
+8.4d3*a1**4*a2**3*d2*dy**4-3.57d3*a1**2*a2**4*dy**4 &
-7.14d3*a1**3*a2**3*dy**4-3.57d3*a1**4*a2**2*dy**4 &
+1.47d3*a1**4*a2**4*d2*dx**4*dy**2-1.3965d4*a1**3*a2**4*dx**4*dy**2 &
-1.3965d4*a1**4*a2**3*dx**4*dy**2-1.68d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.68d4*a1**3*a2**4*d2*dx**2*dy**2+1.68d4*a1**4*a2**3*d2*dx**2*dy**2 &
-7.14d3*a1**2*a2**4*dx**2*dy**2-1.428d4*a1**3*a2**3*dx**2*dy**2 &
-7.14d3*a1**4*a2**2*dx**2*dy**2+4.32d2*a1**4*a2**4*d2**3*dy**2 &
-4.68d3*a1**3*a2**4*d2**2*dy**2-4.68d3*a1**4*a2**3*d2**2*dy**2 &
+5.22d3*a1**2*a2**4*d2*dy**2+1.044d4*a1**3*a2**3*d2*dy**2 &
+5.22d3*a1**4*a2**2*d2*dy**2-2.43d3*a1*a2**4*dy**2 &
-7.29d3*a1**2*a2**3*dy**2-7.29d3*a1**3*a2**2*dy**2 &
-2.43d3*a1**4*a2*dy**2+4.9d2*a1**4*a2**4*d2*dx**6 &
-4.655d3*a1**3*a2**4*dx**6-4.655d3*a1**4*a2**3*dx**6 &
-8.4d2*a1**4*a2**4*d2**2*dx**4+8.4d3*a1**3*a2**4*d2*dx**4 &
+8.4d3*a1**4*a2**3*d2*dx**4-3.57d3*a1**2*a2**4*dx**4 &
-7.14d3*a1**3*a2**3*dx**4-3.57d3*a1**4*a2**2*dx**4 &
+4.32d2*a1**4*a2**4*d2**3*dx**2-4.68d3*a1**3*a2**4*d2**2*dx**2 &
-4.68d3*a1**4*a2**3*d2**2*dx**2+5.22d3*a1**2*a2**4*d2*dx**2 &
+1.044d4*a1**3*a2**3*d2*dx**2+5.22d3*a1**4*a2**2*d2*dx**2 &
-2.43d3*a1*a2**4*dx**2-7.29d3*a1**2*a2**3*dx**2 &
-7.29d3*a1**3*a2**2*dx**2-2.43d3*a1**4*a2*dx**2 &
-6.4d1*a1**4*a2**4*d2**4+8.0d2*a1**3*a2**4*d2**3 &
+8.0d2*a1**4*a2**3*d2**3-1.92d3*a1**2*a2**4*d2**2 &
-3.84d3*a1**3*a2**3*d2**2-1.92d3*a1**4*a2**2*d2**2+2.28d3*a1*a2**4*d2 &
+6.84d3*a1**2*a2**3*d2+6.84d3*a1**3*a2**2*d2+2.28d3*a1**4*a2*d2 &
-7.8d2*a2**4-3.12d3*a1*a2**3-4.68d3*a1**2*a2**2-3.12d3*a1**3*a2 &
-7.8d2*a1**4)*dz
                case (1)
                  rlYlm_laplacian =-1.24625661391794096d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(3.92d2*a1**5*a2**5*d2*dx**2*dy**6 &
-3.724d3*a1**4*a2**5*dx**2*dy**6-3.724d3*a1**5*a2**4*dx**2*dy**6 &
-1.96d2*a1**4*a2**5*d2*dy**6-1.96d2*a1**5*a2**4*d2*dy**6 &
+1.666d3*a1**3*a2**5*dy**6+3.332d3*a1**4*a2**4*dy**6 &
+1.666d3*a1**5*a2**3*dy**6+1.176d3*a1**5*a2**5*d2*dx**4*dy**4 &
-1.1172d4*a1**4*a2**5*dx**4*dy**4-1.1172d4*a1**5*a2**4*dx**4*dy**4 &
-8.4d2*a1**5*a2**5*d2**2*dx**2*dy**4 &
+7.812d3*a1**4*a2**5*d2*dx**2*dy**4+7.812d3*a1**5*a2**4*d2*dx**2*dy**4 &
+1.428d3*a1**3*a2**5*dx**2*dy**4+2.856d3*a1**4*a2**4*dx**2*dy**4 &
+1.428d3*a1**5*a2**3*dx**2*dy**4+4.2d2*a1**4*a2**5*d2**2*dy**4 &
+4.2d2*a1**5*a2**4*d2**2*dy**4-3.948d3*a1**3*a2**5*d2*dy**4 &
-7.896d3*a1**4*a2**4*d2*dy**4-3.948d3*a1**5*a2**3*d2*dy**4 &
+2.835d3*a1**2*a2**5*dy**4+8.505d3*a1**3*a2**4*dy**4 &
+8.505d3*a1**4*a2**3*dy**4+2.835d3*a1**5*a2**2*dy**4 &
+1.176d3*a1**5*a2**5*d2*dx**6*dy**2-1.1172d4*a1**4*a2**5*dx**6*dy**2 &
-1.1172d4*a1**5*a2**4*dx**6*dy**2-1.68d3*a1**5*a2**5*d2**2*dx**4*dy**2 &
+1.6212d4*a1**4*a2**5*d2*dx**4*dy**2 &
+1.6212d4*a1**5*a2**4*d2*dx**4*dy**2-2.142d3*a1**3*a2**5*dx**4*dy**2 &
-4.284d3*a1**4*a2**4*dx**4*dy**2-2.142d3*a1**5*a2**3*dx**4*dy**2 &
+5.76d2*a1**5*a2**5*d2**3*dx**2*dy**2 &
-5.4d3*a1**4*a2**5*d2**2*dx**2*dy**2 &
-5.4d3*a1**5*a2**4*d2**2*dx**2*dy**2-9.36d2*a1**3*a2**5*d2*dx**2*dy**2 &
-1.872d3*a1**4*a2**4*d2*dx**2*dy**2-9.36d2*a1**5*a2**3*d2*dx**2*dy**2 &
+2.43d3*a1**2*a2**5*dx**2*dy**2+7.29d3*a1**3*a2**4*dx**2*dy**2 &
+7.29d3*a1**4*a2**3*dx**2*dy**2+2.43d3*a1**5*a2**2*dx**2*dy**2 &
-2.88d2*a1**4*a2**5*d2**3*dy**2-2.88d2*a1**5*a2**4*d2**3*dy**2 &
+3.096d3*a1**3*a2**5*d2**2*dy**2+6.192d3*a1**4*a2**4*d2**2*dy**2 &
+3.096d3*a1**5*a2**3*d2**2*dy**2-5.184d3*a1**2*a2**5*d2*dy**2 &
-1.5552d4*a1**3*a2**4*d2*dy**2-1.5552d4*a1**4*a2**3*d2*dy**2 &
-5.184d3*a1**5*a2**2*d2*dy**2+2.106d3*a1*a2**5*dy**2 &
+8.424d3*a1**2*a2**4*dy**2+1.2636d4*a1**3*a2**3*dy**2 &
+8.424d3*a1**4*a2**2*dy**2+2.106d3*a1**5*a2*dy**2 &
+3.92d2*a1**5*a2**5*d2*dx**8-3.724d3*a1**4*a2**5*dx**8 &
-3.724d3*a1**5*a2**4*dx**8-8.4d2*a1**5*a2**5*d2**2*dx**6 &
+8.204d3*a1**4*a2**5*d2*dx**6+8.204d3*a1**5*a2**4*d2*dx**6 &
-1.904d3*a1**3*a2**5*dx**6-3.808d3*a1**4*a2**4*dx**6 &
-1.904d3*a1**5*a2**3*dx**6+5.76d2*a1**5*a2**5*d2**3*dx**4 &
-5.82d3*a1**4*a2**5*d2**2*dx**4-5.82d3*a1**5*a2**4*d2**2*dx**4 &
+3.012d3*a1**3*a2**5*d2*dx**4+6.024d3*a1**4*a2**4*d2*dx**4 &
+3.012d3*a1**5*a2**3*d2*dx**4-4.05d2*a1**2*a2**5*dx**4 &
-1.215d3*a1**3*a2**4*dx**4-1.215d3*a1**4*a2**3*dx**4 &
-4.05d2*a1**5*a2**2*dx**4-1.28d2*a1**5*a2**5*d2**4*dx**2 &
+1.312d3*a1**4*a2**5*d2**3*dx**2+1.312d3*a1**5*a2**4*d2**3*dx**2 &
-7.44d2*a1**3*a2**5*d2**2*dx**2-1.488d3*a1**4*a2**4*d2**2*dx**2 &
-7.44d2*a1**5*a2**3*d2**2*dx**2-6.24d2*a1**2*a2**5*d2*dx**2 &
-1.872d3*a1**3*a2**4*d2*dx**2-1.872d3*a1**4*a2**3*d2*dx**2 &
-6.24d2*a1**5*a2**2*d2*dx**2+5.46d2*a1*a2**5*dx**2 &
+2.184d3*a1**2*a2**4*dx**2+3.276d3*a1**3*a2**3*dx**2 &
+2.184d3*a1**4*a2**2*dx**2+5.46d2*a1**5*a2*dx**2 &
+6.4d1*a1**4*a2**5*d2**4+6.4d1*a1**5*a2**4*d2**4 &
-8.32d2*a1**3*a2**5*d2**3-1.664d3*a1**4*a2**4*d2**3 &
-8.32d2*a1**5*a2**3*d2**3+2.52d3*a1**2*a2**5*d2**2 &
+7.56d3*a1**3*a2**4*d2**2+7.56d3*a1**4*a2**3*d2**2 &
+2.52d3*a1**5*a2**2*d2**2-2.424d3*a1*a2**5*d2-9.696d3*a1**2*a2**4*d2 &
-1.4544d4*a1**3*a2**3*d2-9.696d3*a1**4*a2**2*d2-2.424d3*a1**5*a2*d2 &
+4.62d2*a2**5+2.31d3*a1*a2**4+4.62d3*a1**2*a2**3+4.62d3*a1**3*a2**2 &
+2.31d3*a1**4*a2+4.62d2*a1**5)
                case (2)
                  rlYlm_laplacian = &
-3.52494601119984446d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(9.8d1*a1**4*a2**4*d2*dy**6-9.31d2*a1**3*a2**4*dy**6 &
-9.31d2*a1**4*a2**3*dy**6+9.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
-9.31d2*a1**3*a2**4*dx**2*dy**4-9.31d2*a1**4*a2**3*dx**2*dy**4 &
-1.4d2*a1**4*a2**4*d2**2*dy**4+1.498d3*a1**3*a2**4*d2*dy**4 &
+1.498d3*a1**4*a2**3*d2*dy**4-1.428d3*a1**2*a2**4*dy**4 &
-2.856d3*a1**3*a2**3*dy**4-1.428d3*a1**4*a2**2*dy**4 &
-9.8d1*a1**4*a2**4*d2*dx**4*dy**2+9.31d2*a1**3*a2**4*dx**4*dy**2 &
+9.31d2*a1**4*a2**3*dx**4*dy**2+1.96d2*a1**3*a2**4*d2*dx**2*dy**2 &
+1.96d2*a1**4*a2**3*d2*dx**2*dy**2-1.666d3*a1**2*a2**4*dx**2*dy**2 &
-3.332d3*a1**3*a2**3*dx**2*dy**2-1.666d3*a1**4*a2**2*dx**2*dy**2 &
+4.8d1*a1**4*a2**4*d2**3*dy**2-6.6d2*a1**3*a2**4*d2**2*dy**2 &
-6.6d2*a1**4*a2**3*d2**2*dy**2+1.896d3*a1**2*a2**4*d2*dy**2 &
+3.792d3*a1**3*a2**3*d2*dy**2+1.896d3*a1**4*a2**2*d2*dy**2 &
-1.215d3*a1*a2**4*dy**2-3.645d3*a1**2*a2**3*dy**2 &
-3.645d3*a1**3*a2**2*dy**2-1.215d3*a1**4*a2*dy**2 &
-9.8d1*a1**4*a2**4*d2*dx**6+9.31d2*a1**3*a2**4*dx**6 &
+9.31d2*a1**4*a2**3*dx**6+1.4d2*a1**4*a2**4*d2**2*dx**4 &
-1.302d3*a1**3*a2**4*d2*dx**4-1.302d3*a1**4*a2**3*d2*dx**4 &
-2.38d2*a1**2*a2**4*dx**4-4.76d2*a1**3*a2**3*dx**4 &
-2.38d2*a1**4*a2**2*dx**4-4.8d1*a1**4*a2**4*d2**3*dx**2 &
+3.8d2*a1**3*a2**4*d2**2*dx**2+3.8d2*a1**4*a2**3*d2**2*dx**2 &
+7.36d2*a1**2*a2**4*d2*dx**2+1.472d3*a1**3*a2**3*d2*dx**2 &
+7.36d2*a1**4*a2**2*d2*dx**2-6.75d2*a1*a2**4*dx**2 &
-2.025d3*a1**2*a2**3*dx**2-2.025d3*a1**3*a2**2*dx**2 &
-6.75d2*a1**4*a2*dx**2+4.8d1*a1**3*a2**4*d2**3+4.8d1*a1**4*a2**3*d2**3 &
-5.16d2*a1**2*a2**4*d2**2-1.032d3*a1**3*a2**3*d2**2 &
-5.16d2*a1**4*a2**2*d2**2+8.64d2*a1*a2**4*d2+2.592d3*a1**2*a2**3*d2 &
+2.592d3*a1**3*a2**2*d2+8.64d2*a1**4*a2*d2-3.51d2*a2**4 &
-1.404d3*a1*a2**3-2.106d3*a1**2*a2**2-1.404d3*a1**3*a2-3.51d2*a1**4 &
)*dz
                case (3)
                  rlYlm_laplacian = &
-3.29728507019630958d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)* &
(1.68d2*a1**4*a2**4*d2*dx**2*dy**6-1.596d3*a1**3*a2**4*dx**2*dy**6 &
-1.596d3*a1**4*a2**3*dx**2*dy**6-8.4d1*a1**3*a2**4*d2*dy**6 &
-8.4d1*a1**4*a2**3*d2*dy**6+7.14d2*a1**2*a2**4*dy**6 &
+1.428d3*a1**3*a2**3*dy**6+7.14d2*a1**4*a2**2*dy**6 &
+2.8d2*a1**4*a2**4*d2*dx**4*dy**4-2.66d3*a1**3*a2**4*dx**4*dy**4 &
-2.66d3*a1**4*a2**3*dx**4*dy**4-2.64d2*a1**4*a2**4*d2**2*dx**2*dy**4 &
+2.388d3*a1**3*a2**4*d2*dx**2*dy**4+2.388d3*a1**4*a2**3*d2*dx**2*dy**4 &
+1.02d3*a1**2*a2**4*dx**2*dy**4+2.04d3*a1**3*a2**3*dx**2*dy**4 &
+1.02d3*a1**4*a2**2*dx**2*dy**4+1.32d2*a1**3*a2**4*d2**2*dy**4 &
+1.32d2*a1**4*a2**3*d2**2*dy**4-1.14d3*a1**2*a2**4*d2*dy**4 &
-2.28d3*a1**3*a2**3*d2*dy**4-1.14d3*a1**4*a2**2*d2*dy**4 &
+1.35d2*a1*a2**4*dy**4+4.05d2*a1**2*a2**3*dy**4 &
+4.05d2*a1**3*a2**2*dy**4+1.35d2*a1**4*a2*dy**4 &
+5.6d1*a1**4*a2**4*d2*dx**6*dy**2-5.32d2*a1**3*a2**4*dx**6*dy**2 &
-5.32d2*a1**4*a2**3*dx**6*dy**2-1.76d2*a1**4*a2**4*d2**2*dx**4*dy**2 &
+1.732d3*a1**3*a2**4*d2*dx**4*dy**2+1.732d3*a1**4*a2**3*d2*dx**4*dy**2 &
-5.1d2*a1**2*a2**4*dx**4*dy**2-1.02d3*a1**3*a2**3*dx**4*dy**2 &
-5.1d2*a1**4*a2**2*dx**4*dy**2+9.6d1*a1**4*a2**4*d2**3*dx**2*dy**2 &
-8.4d2*a1**3*a2**4*d2**2*dx**2*dy**2 &
-8.4d2*a1**4*a2**3*d2**2*dx**2*dy**2-7.2d2*a1**2*a2**4*d2*dx**2*dy**2 &
-1.44d3*a1**3*a2**3*d2*dx**2*dy**2-7.2d2*a1**4*a2**2*d2*dx**2*dy**2 &
+8.1d2*a1*a2**4*dx**2*dy**2+2.43d3*a1**2*a2**3*dx**2*dy**2 &
+2.43d3*a1**3*a2**2*dx**2*dy**2+8.1d2*a1**4*a2*dx**2*dy**2 &
-4.8d1*a1**3*a2**4*d2**3*dy**2-4.8d1*a1**4*a2**3*d2**3*dy**2 &
+4.08d2*a1**2*a2**4*d2**2*dy**2+8.16d2*a1**3*a2**3*d2**2*dy**2 &
+4.08d2*a1**4*a2**2*d2**2*dy**2+3.6d1*a1*a2**4*d2*dy**2 &
+1.08d2*a1**2*a2**3*d2*dy**2+1.08d2*a1**3*a2**2*d2*dy**2 &
+3.6d1*a1**4*a2*d2*dy**2-2.34d2*a2**4*dy**2-9.36d2*a1*a2**3*dy**2 &
-1.404d3*a1**2*a2**2*dy**2-9.36d2*a1**3*a2*dy**2-2.34d2*a1**4*dy**2 &
-5.6d1*a1**4*a2**4*d2*dx**8+5.32d2*a1**3*a2**4*dx**8 &
+5.32d2*a1**4*a2**3*dx**8+8.8d1*a1**4*a2**4*d2**2*dx**6 &
-7.4d2*a1**3*a2**4*d2*dx**6-7.4d2*a1**4*a2**3*d2*dx**6 &
-8.16d2*a1**2*a2**4*dx**6-1.632d3*a1**3*a2**3*dx**6 &
-8.16d2*a1**4*a2**2*dx**6-3.2d1*a1**4*a2**4*d2**3*dx**4 &
+1.48d2*a1**3*a2**4*d2**2*dx**4+1.48d2*a1**4*a2**3*d2**2*dx**4 &
+1.38d3*a1**2*a2**4*d2*dx**4+2.76d3*a1**3*a2**3*d2*dx**4 &
+1.38d3*a1**4*a2**2*d2*dx**4-4.05d2*a1*a2**4*dx**4 &
-1.215d3*a1**2*a2**3*dx**4-1.215d3*a1**3*a2**2*dx**4 &
-4.05d2*a1**4*a2*dx**4+4.8d1*a1**3*a2**4*d2**3*dx**2 &
+4.8d1*a1**4*a2**3*d2**3*dx**2-4.08d2*a1**2*a2**4*d2**2*dx**2 &
-8.16d2*a1**3*a2**3*d2**2*dx**2-4.08d2*a1**4*a2**2*d2**2*dx**2 &
-3.6d1*a1*a2**4*d2*dx**2-1.08d2*a1**2*a2**3*d2*dx**2 &
-1.08d2*a1**3*a2**2*d2*dx**2-3.6d1*a1**4*a2*d2*dx**2 &
+2.34d2*a2**4*dx**2+9.36d2*a1*a2**3*dx**2+1.404d3*a1**2*a2**2*dx**2 &
+9.36d2*a1**3*a2*dx**2+2.34d2*a1**4*dx**2)
                case (4)
                  rlYlm_laplacian = &
-4.66306526528194375d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(1.4d1*a1**3*a2**3*d2*dy**6-1.33d2*a1**2*a2**3*dy**6 &
-1.33d2*a1**3*a2**2*dy**6-7.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+6.65d2*a1**2*a2**3*dx**2*dy**4+6.65d2*a1**3*a2**2*dx**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**4+1.36d2*a1**2*a2**3*d2*dy**4 &
+1.36d2*a1**3*a2**2*d2*dy**4-5.1d2*a1*a2**3*dy**4 &
-1.02d3*a1**2*a2**2*dy**4-5.1d2*a1**3*a2*dy**4 &
-7.0d1*a1**3*a2**3*d2*dx**4*dy**2+6.65d2*a1**2*a2**3*dx**4*dy**2 &
+6.65d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-2.56d2*a1**2*a2**3*d2*dx**2*dy**2-2.56d2*a1**3*a2**2*d2*dx**2*dy**2 &
-1.7d3*a1*a2**3*dx**2*dy**2-3.4d3*a1**2*a2**2*dx**2*dy**2 &
-1.7d3*a1**3*a2*dx**2*dy**2-4.8d1*a1**2*a2**3*d2**2*dy**2 &
-4.8d1*a1**3*a2**2*d2**2*dy**2+3.0d2*a1*a2**3*d2*dy**2 &
+6.0d2*a1**2*a2**2*d2*dy**2+3.0d2*a1**3*a2*d2*dy**2+8.1d2*a2**3*dy**2 &
+2.43d3*a1*a2**2*dy**2+2.43d3*a1**2*a2*dy**2+8.1d2*a1**3*dy**2 &
+1.4d1*a1**3*a2**3*d2*dx**6-1.33d2*a1**2*a2**3*dx**6 &
-1.33d2*a1**3*a2**2*dx**6-8.0d0*a1**3*a2**3*d2**2*dx**4 &
+2.4d1*a1**2*a2**3*d2*dx**4+2.4d1*a1**3*a2**2*d2*dx**4 &
+4.42d2*a1*a2**3*dx**4+8.84d2*a1**2*a2**2*dx**4+4.42d2*a1**3*a2*dx**4 &
+1.6d1*a1**2*a2**3*d2**2*dx**2+1.6d1*a1**3*a2**2*d2**2*dx**2 &
-1.0d2*a1*a2**3*d2*dx**2-2.0d2*a1**2*a2**2*d2*dx**2 &
-1.0d2*a1**3*a2*d2*dx**2-2.7d2*a2**3*dx**2-8.1d2*a1*a2**2*dx**2 &
-8.1d2*a1**2*a2*dx**2-2.7d2*a1**3*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (2)
          ! selection on l2: l1=4, m1=2
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=2, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =1.48624773660225413d0*E*a1*a2**5*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*(dy+dx)* &
(dy-1.0d0*dx)*(7.0d0*(dy**2+dx**2)-6.0d0*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=2, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =2.5742565924293503d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(1.4d1*a1**2*a2**2*d2*dy**4 &
-9.1d1*a1*a2**2*dy**4-9.1d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2**2*dy**2+6.2d1*a1*a2**2*d2*dy**2 &
+6.2d1*a1**2*a2*d2*dy**2+8.8d1*a2**2*dy**2+1.76d2*a1*a2*dy**2 &
+8.8d1*a1**2*dy**2-1.4d1*a1**2*a2**2*d2*dx**4+9.1d1*a1*a2**2*dx**4 &
+9.1d1*a1**2*a2*dx**4+1.2d1*a1**2*a2**2*d2**2*dx**2 &
-9.0d1*a1*a2**2*d2*dx**2-9.0d1*a1**2*a2*d2*dx**2+6.6d1*a2**2*dx**2 &
+1.32d2*a1*a2*dx**2+6.6d1*a1**2*dx**2+1.2d1*a1*a2**2*d2**2 &
+1.2d1*a1**2*a2*d2**2-6.6d1*a2**2*d2-1.32d2*a1*a2*d2-6.6d1*a1**2*d2)
                case (0)
                  rlYlm_laplacian =2.5742565924293503d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**2*a2**2*d2*dy**2-9.1d1*a1*a2**2*dy**2-9.1d1*a1**2*a2*dy**2 &
+1.4d1*a1**2*a2**2*d2*dx**2-9.1d1*a1*a2**2*dx**2-9.1d1*a1**2*a2*dx**2 &
-1.2d1*a1**2*a2**2*d2**2+9.0d1*a1*a2**2*d2+9.0d1*a1**2*a2*d2 &
-6.6d1*a2**2-1.32d2*a1*a2-6.6d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =2.5742565924293503d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(1.4d1*a1**2*a2**2*d2*dy**4 &
-9.1d1*a1*a2**2*dy**4-9.1d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2**2*dy**2+9.0d1*a1*a2**2*d2*dy**2 &
+9.0d1*a1**2*a2*d2*dy**2-6.6d1*a2**2*dy**2-1.32d2*a1*a2*dy**2 &
-6.6d1*a1**2*dy**2-1.4d1*a1**2*a2**2*d2*dx**4+9.1d1*a1*a2**2*dx**4 &
+9.1d1*a1**2*a2*dx**4+1.2d1*a1**2*a2**2*d2**2*dx**2 &
-6.2d1*a1*a2**2*d2*dx**2-6.2d1*a1**2*a2*d2*dx**2-8.8d1*a2**2*dx**2 &
-1.76d2*a1*a2*dx**2-8.8d1*a1**2*dx**2-1.2d1*a1*a2**2*d2**2 &
-1.2d1*a1**2*a2*d2**2+6.6d1*a2**2*d2+1.32d2*a1*a2*d2+6.6d1*a1**2*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=2, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =5.75621273219899775d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**2*a2**2*d2*dy**2-1.05d2*a1*a2**2*dy**2 &
-1.05d2*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-1.05d2*a1*a2**2*dx**2-1.05d2*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2 &
+8.6d1*a1*a2**2*d2+8.6d1*a1**2*a2*d2+2.6d1*a2**2+5.2d1*a1*a2 &
+2.6d1*a1**2)
                case (-1)
                  rlYlm_laplacian =5.75621273219899775d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.05d2*a1**2*a2**3*dy**4-1.05d2*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2**2*dy**2+8.6d1*a1**2*a2**3*d2*dy**2 &
+8.6d1*a1**3*a2**2*d2*dy**2+2.6d1*a1*a2**3*dy**2 &
+5.2d1*a1**2*a2**2*dy**2+2.6d1*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.05d2*a1**2*a2**3*dx**4 &
+1.05d2*a1**3*a2**2*dx**4+1.2d1*a1**3*a2**3*d2**2*dx**2 &
-1.14d2*a1**2*a2**3*d2*dx**2-1.14d2*a1**3*a2**2*d2*dx**2 &
+1.56d2*a1*a2**3*dx**2+3.12d2*a1**2*a2**2*dx**2+1.56d2*a1**3*a2*dx**2 &
+1.2d1*a1**2*a2**3*d2**2+1.2d1*a1**3*a2**2*d2**2-9.0d1*a1*a2**3*d2 &
-1.8d2*a1**2*a2**2*d2-9.0d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian = &
-1.66167548522392128d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*(dy+dx)* &
(dy-1.0d0*dx)*(4.2d1*a1**3*a2**3*d2*dy**4-3.15d2*a1**2*a2**3*dy**4 &
-3.15d2*a1**3*a2**2*dy**4+8.4d1*a1**3*a2**3*d2*dx**2*dy**2 &
-6.3d2*a1**2*a2**3*dx**2*dy**2-6.3d2*a1**3*a2**2*dx**2*dy**2 &
-6.4d1*a1**3*a2**3*d2**2*dy**2+4.96d2*a1**2*a2**3*d2*dy**2 &
+4.96d2*a1**3*a2**2*d2*dy**2-1.04d2*a1*a2**3*dy**2 &
-2.08d2*a1**2*a2**2*dy**2-1.04d2*a1**3*a2*dy**2 &
+4.2d1*a1**3*a2**3*d2*dx**4-3.15d2*a1**2*a2**3*dx**4 &
-3.15d2*a1**3*a2**2*dx**4-6.4d1*a1**3*a2**3*d2**2*dx**2 &
+4.96d2*a1**2*a2**3*d2*dx**2+4.96d2*a1**3*a2**2*d2*dx**2 &
-1.04d2*a1*a2**3*dx**2-2.08d2*a1**2*a2**2*dx**2-1.04d2*a1**3*a2*dx**2 &
+2.4d1*a1**3*a2**3*d2**3-2.04d2*a1**2*a2**3*d2**2 &
-2.04d2*a1**3*a2**2*d2**2+1.74d2*a1*a2**3*d2+3.48d2*a1**2*a2**2*d2 &
+1.74d2*a1**3*a2*d2-9.9d1*a2**3-2.97d2*a1*a2**2-2.97d2*a1**2*a2 &
-9.9d1*a1**3)
                case (1)
                  rlYlm_laplacian =5.75621273219899775d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(1.4d1*a1**3*a2**3*d2*dy**4 &
-1.05d2*a1**2*a2**3*dy**4-1.05d2*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2**2*dy**2+1.14d2*a1**2*a2**3*d2*dy**2 &
+1.14d2*a1**3*a2**2*d2*dy**2-1.56d2*a1*a2**3*dy**2 &
-3.12d2*a1**2*a2**2*dy**2-1.56d2*a1**3*a2*dy**2 &
-1.4d1*a1**3*a2**3*d2*dx**4+1.05d2*a1**2*a2**3*dx**4 &
+1.05d2*a1**3*a2**2*dx**4+1.2d1*a1**3*a2**3*d2**2*dx**2 &
-8.6d1*a1**2*a2**3*d2*dx**2-8.6d1*a1**3*a2**2*d2*dx**2 &
-2.6d1*a1*a2**3*dx**2-5.2d1*a1**2*a2**2*dx**2-2.6d1*a1**3*a2*dx**2 &
-1.2d1*a1**2*a2**3*d2**2-1.2d1*a1**3*a2**2*d2**2+9.0d1*a1*a2**3*d2 &
+1.8d2*a1**2*a2**2*d2+9.0d1*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2 &
-1.98d2*a1**2*a2-6.6d1*a1**3)*dz
                case (2)
                  rlYlm_laplacian = &
-2.87810636609949888d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)* &
(1.4d1*a1**3*a2**3*d2*dy**6-1.05d2*a1**2*a2**3*dy**6 &
-1.05d2*a1**3*a2**2*dy**6-1.4d1*a1**3*a2**3*d2*dx**2*dy**4 &
+1.05d2*a1**2*a2**3*dx**2*dy**4+1.05d2*a1**3*a2**2*dx**2*dy**4 &
-1.2d1*a1**3*a2**3*d2**2*dy**4+5.8d1*a1**2*a2**3*d2*dy**4 &
+5.8d1*a1**3*a2**2*d2*dy**4+2.08d2*a1*a2**3*dy**4 &
+4.16d2*a1**2*a2**2*dy**4+2.08d2*a1**3*a2*dy**4 &
-1.4d1*a1**3*a2**3*d2*dx**4*dy**2+1.05d2*a1**2*a2**3*dx**4*dy**2 &
+1.05d2*a1**3*a2**2*dx**4*dy**2+2.4d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-2.28d2*a1**2*a2**3*d2*dx**2*dy**2-2.28d2*a1**3*a2**2*d2*dx**2*dy**2 &
+3.12d2*a1*a2**3*dx**2*dy**2+6.24d2*a1**2*a2**2*dx**2*dy**2 &
+3.12d2*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.38d2*a1*a2**3*d2*dy**2 &
-2.76d2*a1**2*a2**2*d2*dy**2-1.38d2*a1**3*a2*d2*dy**2 &
-9.9d1*a2**3*dy**2-2.97d2*a1*a2**2*dy**2-2.97d2*a1**2*a2*dy**2 &
-9.9d1*a1**3*dy**2+1.4d1*a1**3*a2**3*d2*dx**6-1.05d2*a1**2*a2**3*dx**6 &
-1.05d2*a1**3*a2**2*dx**6-1.2d1*a1**3*a2**3*d2**2*dx**4 &
+5.8d1*a1**2*a2**3*d2*dx**4+5.8d1*a1**3*a2**2*d2*dx**4 &
+2.08d2*a1*a2**3*dx**4+4.16d2*a1**2*a2**2*dx**4+2.08d2*a1**3*a2*dx**4 &
+2.4d1*a1**2*a2**3*d2**2*dx**2+2.4d1*a1**3*a2**2*d2**2*dx**2 &
-1.38d2*a1*a2**3*d2*dx**2-2.76d2*a1**2*a2**2*d2*dx**2 &
-1.38d2*a1**3*a2*d2*dx**2-9.9d1*a2**3*dx**2-2.97d2*a1*a2**2*dx**2 &
-2.97d2*a1**2*a2*dx**2-9.9d1*a1**3*dx**2-1.2d1*a1*a2**3*d2**2 &
-2.4d1*a1**2*a2**2*d2**2-1.2d1*a1**3*a2*d2**2+6.6d1*a2**3*d2 &
+1.98d2*a1*a2**2*d2+1.98d2*a1**2*a2*d2+6.6d1*a1**3*d2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=2, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-3.10871017685462917d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.19d2*a1**3*a2**4*dy**6 &
-1.19d2*a1**4*a2**3*dy**6-4.2d1*a1**4*a2**4*d2*dx**2*dy**4 &
+3.57d2*a1**3*a2**4*dx**2*dy**4+3.57d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+5.4d1*a1**3*a2**4*d2*dy**4 &
+5.4d1*a1**4*a2**3*d2*dy**4+3.6d2*a1**2*a2**4*dy**4 &
+7.2d2*a1**3*a2**3*dy**4+3.6d2*a1**4*a2**2*dy**4 &
-1.4d1*a1**4*a2**4*d2*dx**4*dy**2+1.19d2*a1**3*a2**4*dx**4*dy**2 &
+1.19d2*a1**4*a2**3*dx**4*dy**2+4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-4.68d2*a1**3*a2**4*d2*dx**2*dy**2-4.68d2*a1**4*a2**3*d2*dx**2*dy**2 &
+4.5d2*a1**2*a2**4*dx**2*dy**2+9.0d2*a1**3*a2**3*dx**2*dy**2 &
+4.5d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**3*a2**4*d2**2*dy**2 &
+3.6d1*a1**4*a2**3*d2**2*dy**2-2.16d2*a1**2*a2**4*d2*dy**2 &
-4.32d2*a1**3*a2**3*d2*dy**2-2.16d2*a1**4*a2**2*d2*dy**2 &
-3.51d2*a1*a2**4*dy**2-1.053d3*a1**2*a2**3*dy**2 &
-1.053d3*a1**3*a2**2*dy**2-3.51d2*a1**4*a2*dy**2 &
+4.2d1*a1**4*a2**4*d2*dx**6-3.57d2*a1**3*a2**4*dx**6 &
-3.57d2*a1**4*a2**3*dx**6-3.6d1*a1**4*a2**4*d2**2*dx**4 &
+2.46d2*a1**3*a2**4*d2*dx**4+2.46d2*a1**4*a2**3*d2*dx**4 &
+4.5d2*a1**2*a2**4*dx**4+9.0d2*a1**3*a2**3*dx**4 &
+4.5d2*a1**4*a2**2*dx**4+3.6d1*a1**3*a2**4*d2**2*dx**2 &
+3.6d1*a1**4*a2**3*d2**2*dx**2-2.16d2*a1**2*a2**4*d2*dx**2 &
-4.32d2*a1**3*a2**3*d2*dx**2-2.16d2*a1**4*a2**2*d2*dx**2 &
-3.51d2*a1*a2**4*dx**2-1.053d3*a1**2*a2**3*dx**2 &
-1.053d3*a1**3*a2**2*dx**2-3.51d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.28d2*a1*a2**4*d2+6.84d2*a1**2*a2**3*d2 &
+6.84d2*a1**3*a2**2*d2+2.28d2*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case (-2)
                  rlYlm_laplacian =1.52295073829821874d1*E*a1**3*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**2*a2**2*d2*dy**2-1.19d2*a1*a2**2*dy**2 &
-1.19d2*a1**2*a2*dy**2+1.4d1*a1**2*a2**2*d2*dx**2 &
-1.19d2*a1*a2**2*dx**2-1.19d2*a1**2*a2*dx**2-1.2d1*a1**2*a2**2*d2**2 &
+1.1d2*a1*a2**2*d2+1.1d2*a1**2*a2*d2-6.0d1*a2**2-1.2d2*a1*a2 &
-6.0d1*a1**2)*dz
                case (-1)
                  rlYlm_laplacian = &
-2.40799654862869848d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(7.0d1*a1**4*a2**4*d2*dy**6-5.95d2*a1**3*a2**4*dy**6 &
-5.95d2*a1**4*a2**3*dy**6+7.0d1*a1**4*a2**4*d2*dx**2*dy**4 &
-5.95d2*a1**3*a2**4*dx**2*dy**4-5.95d2*a1**4*a2**3*dx**2*dy**4 &
-1.16d2*a1**4*a2**4*d2**2*dy**4+9.7d2*a1**3*a2**4*d2*dy**4 &
+9.7d2*a1**4*a2**3*d2*dy**4+1.2d2*a1**2*a2**4*dy**4 &
+2.4d2*a1**3*a2**3*dy**4+1.2d2*a1**4*a2**2*dy**4 &
-7.0d1*a1**4*a2**4*d2*dx**4*dy**2+5.95d2*a1**3*a2**4*dx**4*dy**2 &
+5.95d2*a1**4*a2**3*dx**4*dy**2-1.4d2*a1**3*a2**4*d2*dx**2*dy**2 &
-1.4d2*a1**4*a2**3*d2*dx**2*dy**2+1.05d3*a1**2*a2**4*dx**2*dy**2 &
+2.1d3*a1**3*a2**3*dx**2*dy**2+1.05d3*a1**4*a2**2*dx**2*dy**2 &
+4.8d1*a1**4*a2**4*d2**3*dy**2-3.56d2*a1**3*a2**4*d2**2*dy**2 &
-3.56d2*a1**4*a2**3*d2**2*dy**2-4.32d2*a1**2*a2**4*d2*dy**2 &
-8.64d2*a1**3*a2**3*d2*dy**2-4.32d2*a1**4*a2**2*d2*dy**2 &
+2.73d2*a1*a2**4*dy**2+8.19d2*a1**2*a2**3*dy**2 &
+8.19d2*a1**3*a2**2*dy**2+2.73d2*a1**4*a2*dy**2 &
-7.0d1*a1**4*a2**4*d2*dx**6+5.95d2*a1**3*a2**4*dx**6 &
+5.95d2*a1**4*a2**3*dx**6+1.16d2*a1**4*a2**4*d2**2*dx**4 &
-1.11d3*a1**3*a2**4*d2*dx**4-1.11d3*a1**4*a2**3*d2*dx**4 &
+9.3d2*a1**2*a2**4*dx**4+1.86d3*a1**3*a2**3*dx**4 &
+9.3d2*a1**4*a2**2*dx**4-4.8d1*a1**4*a2**4*d2**3*dx**2 &
+5.88d2*a1**3*a2**4*d2**2*dx**2+5.88d2*a1**4*a2**3*d2**2*dx**2 &
-1.464d3*a1**2*a2**4*d2*dx**2-2.928d3*a1**3*a2**3*d2*dx**2 &
-1.464d3*a1**4*a2**2*d2*dx**2+7.41d2*a1*a2**4*dx**2 &
+2.223d3*a1**2*a2**3*dx**2+2.223d3*a1**3*a2**2*dx**2 &
+7.41d2*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+4.44d2*a1**2*a2**4*d2**2+8.88d2*a1**3*a2**3*d2**2 &
+4.44d2*a1**4*a2**2*d2**2-5.76d2*a1*a2**4*d2-1.728d3*a1**2*a2**3*d2 &
-1.728d3*a1**3*a2**2*d2-5.76d2*a1**4*a2*d2+1.65d2*a2**4+6.6d2*a1*a2**3 &
+9.9d2*a1**2*a2**2+6.6d2*a1**3*a2+1.65d2*a1**4)
                case (0)
                  rlYlm_laplacian = &
-1.96612094884109708d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*(dy+dx &
)*(dy-1.0d0*dx)*(7.0d1*a1**3*a2**3*d2*dy**4-5.95d2*a1**2*a2**3*dy**4 &
-5.95d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.19d3*a1**2*a2**3*dx**2*dy**2-1.19d3*a1**3*a2**2*dx**2*dy**2 &
-8.8d1*a1**3*a2**3*d2**2*dy**2+7.6d2*a1**2*a2**3*d2*dy**2 &
+7.6d2*a1**3*a2**2*d2*dy**2-9.0d1*a1*a2**3*dy**2 &
-1.8d2*a1**2*a2**2*dy**2-9.0d1*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-5.95d2*a1**2*a2**3*dx**4 &
-5.95d2*a1**3*a2**2*dx**4-8.8d1*a1**3*a2**3*d2**2*dx**2 &
+7.6d2*a1**2*a2**3*d2*dx**2+7.6d2*a1**3*a2**2*d2*dx**2 &
-9.0d1*a1*a2**3*dx**2-1.8d2*a1**2*a2**2*dx**2-9.0d1*a1**3*a2*dx**2 &
+2.4d1*a1**3*a2**3*d2**3-2.04d2*a1**2*a2**3*d2**2 &
-2.04d2*a1**3*a2**2*d2**2-1.8d1*a1*a2**3*d2-3.6d1*a1**2*a2**2*d2 &
-1.8d1*a1**3*a2*d2+1.17d2*a2**3+3.51d2*a1*a2**2+3.51d2*a1**2*a2 &
+1.17d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian = &
-2.40799654862869848d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(7.0d1*a1**4*a2**4*d2*dy**6-5.95d2*a1**3*a2**4*dy**6 &
-5.95d2*a1**4*a2**3*dy**6+7.0d1*a1**4*a2**4*d2*dx**2*dy**4 &
-5.95d2*a1**3*a2**4*dx**2*dy**4-5.95d2*a1**4*a2**3*dx**2*dy**4 &
-1.16d2*a1**4*a2**4*d2**2*dy**4+1.11d3*a1**3*a2**4*d2*dy**4 &
+1.11d3*a1**4*a2**3*d2*dy**4-9.3d2*a1**2*a2**4*dy**4 &
-1.86d3*a1**3*a2**3*dy**4-9.3d2*a1**4*a2**2*dy**4 &
-7.0d1*a1**4*a2**4*d2*dx**4*dy**2+5.95d2*a1**3*a2**4*dx**4*dy**2 &
+5.95d2*a1**4*a2**3*dx**4*dy**2+1.4d2*a1**3*a2**4*d2*dx**2*dy**2 &
+1.4d2*a1**4*a2**3*d2*dx**2*dy**2-1.05d3*a1**2*a2**4*dx**2*dy**2 &
-2.1d3*a1**3*a2**3*dx**2*dy**2-1.05d3*a1**4*a2**2*dx**2*dy**2 &
+4.8d1*a1**4*a2**4*d2**3*dy**2-5.88d2*a1**3*a2**4*d2**2*dy**2 &
-5.88d2*a1**4*a2**3*d2**2*dy**2+1.464d3*a1**2*a2**4*d2*dy**2 &
+2.928d3*a1**3*a2**3*d2*dy**2+1.464d3*a1**4*a2**2*d2*dy**2 &
-7.41d2*a1*a2**4*dy**2-2.223d3*a1**2*a2**3*dy**2 &
-2.223d3*a1**3*a2**2*dy**2-7.41d2*a1**4*a2*dy**2 &
-7.0d1*a1**4*a2**4*d2*dx**6+5.95d2*a1**3*a2**4*dx**6 &
+5.95d2*a1**4*a2**3*dx**6+1.16d2*a1**4*a2**4*d2**2*dx**4 &
-9.7d2*a1**3*a2**4*d2*dx**4-9.7d2*a1**4*a2**3*d2*dx**4 &
-1.2d2*a1**2*a2**4*dx**4-2.4d2*a1**3*a2**3*dx**4 &
-1.2d2*a1**4*a2**2*dx**4-4.8d1*a1**4*a2**4*d2**3*dx**2 &
+3.56d2*a1**3*a2**4*d2**2*dx**2+3.56d2*a1**4*a2**3*d2**2*dx**2 &
+4.32d2*a1**2*a2**4*d2*dx**2+8.64d2*a1**3*a2**3*d2*dx**2 &
+4.32d2*a1**4*a2**2*d2*dx**2-2.73d2*a1*a2**4*dx**2 &
-8.19d2*a1**2*a2**3*dx**2-8.19d2*a1**3*a2**2*dx**2 &
-2.73d2*a1**4*a2*dx**2+4.8d1*a1**3*a2**4*d2**3+4.8d1*a1**4*a2**3*d2**3 &
-4.44d2*a1**2*a2**4*d2**2-8.88d2*a1**3*a2**3*d2**2 &
-4.44d2*a1**4*a2**2*d2**2+5.76d2*a1*a2**4*d2+1.728d3*a1**2*a2**3*d2 &
+1.728d3*a1**3*a2**2*d2+5.76d2*a1**4*a2*d2-1.65d2*a2**4-6.6d2*a1*a2**3 &
-9.9d2*a1**2*a2**2-6.6d2*a1**3*a2-1.65d2*a1**4)
                case (2)
                  rlYlm_laplacian = &
-7.6147536914910937d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.19d2*a1**3*a2**4*dy**6 &
-1.19d2*a1**4*a2**3*dy**6-1.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.19d2*a1**3*a2**4*dx**2*dy**4+1.19d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+8.2d1*a1**3*a2**4*d2*dy**4 &
+8.2d1*a1**4*a2**3*d2*dy**4+1.5d2*a1**2*a2**4*dy**4 &
+3.0d2*a1**3*a2**3*dy**4+1.5d2*a1**4*a2**2*dy**4 &
-1.4d1*a1**4*a2**4*d2*dx**4*dy**2+1.19d2*a1**3*a2**4*dx**4*dy**2 &
+1.19d2*a1**4*a2**3*dx**4*dy**2+2.4d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-2.76d2*a1**3*a2**4*d2*dx**2*dy**2-2.76d2*a1**4*a2**3*d2*dx**2*dy**2 &
+5.4d2*a1**2*a2**4*dx**2*dy**2+1.08d3*a1**3*a2**3*dx**2*dy**2 &
+5.4d2*a1**4*a2**2*dx**2*dy**2+2.4d1*a1**3*a2**4*d2**2*dy**2 &
+2.4d1*a1**4*a2**3*d2**2*dy**2-1.86d2*a1**2*a2**4*d2*dy**2 &
-3.72d2*a1**3*a2**3*d2*dy**2-1.86d2*a1**4*a2**2*d2*dy**2 &
+3.9d1*a1*a2**4*dy**2+1.17d2*a1**2*a2**3*dy**2 &
+1.17d2*a1**3*a2**2*dy**2+3.9d1*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.19d2*a1**3*a2**4*dx**6 &
-1.19d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+8.2d1*a1**3*a2**4*d2*dx**4+8.2d1*a1**4*a2**3*d2*dx**4 &
+1.5d2*a1**2*a2**4*dx**4+3.0d2*a1**3*a2**3*dx**4 &
+1.5d2*a1**4*a2**2*dx**4+2.4d1*a1**3*a2**4*d2**2*dx**2 &
+2.4d1*a1**4*a2**3*d2**2*dx**2-1.86d2*a1**2*a2**4*d2*dx**2 &
-3.72d2*a1**3*a2**3*d2*dx**2-1.86d2*a1**4*a2**2*d2*dx**2 &
+3.9d1*a1*a2**4*dx**2+1.17d2*a1**2*a2**3*dx**2 &
+1.17d2*a1**3*a2**2*dx**2+3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2 &
-2.4d1*a1**3*a2**3*d2**2-1.2d1*a1**4*a2**2*d2**2+9.0d1*a1*a2**4*d2 &
+2.7d2*a1**2*a2**3*d2+2.7d2*a1**3*a2**2*d2+9.0d1*a1**4*a2*d2 &
-6.6d1*a2**4-2.64d2*a1*a2**3-3.96d2*a1**2*a2**2-2.64d2*a1**3*a2 &
-6.6d1*a1**4)*dz
                case (3)
                  rlYlm_laplacian = &
-3.10871017685462917d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(4.2d1*a1**4*a2**4*d2*dy**6-3.57d2*a1**3*a2**4*dy**6 &
-3.57d2*a1**4*a2**3*dy**6-1.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.19d2*a1**3*a2**4*dx**2*dy**4+1.19d2*a1**4*a2**3*dx**2*dy**4 &
-3.6d1*a1**4*a2**4*d2**2*dy**4+2.46d2*a1**3*a2**4*d2*dy**4 &
+2.46d2*a1**4*a2**3*d2*dy**4+4.5d2*a1**2*a2**4*dy**4 &
+9.0d2*a1**3*a2**3*dy**4+4.5d2*a1**4*a2**2*dy**4 &
-4.2d1*a1**4*a2**4*d2*dx**4*dy**2+3.57d2*a1**3*a2**4*dx**4*dy**2 &
+3.57d2*a1**4*a2**3*dx**4*dy**2+4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-4.68d2*a1**3*a2**4*d2*dx**2*dy**2-4.68d2*a1**4*a2**3*d2*dx**2*dy**2 &
+4.5d2*a1**2*a2**4*dx**2*dy**2+9.0d2*a1**3*a2**3*dx**2*dy**2 &
+4.5d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**3*a2**4*d2**2*dy**2 &
+3.6d1*a1**4*a2**3*d2**2*dy**2-2.16d2*a1**2*a2**4*d2*dy**2 &
-4.32d2*a1**3*a2**3*d2*dy**2-2.16d2*a1**4*a2**2*d2*dy**2 &
-3.51d2*a1*a2**4*dy**2-1.053d3*a1**2*a2**3*dy**2 &
-1.053d3*a1**3*a2**2*dy**2-3.51d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.19d2*a1**3*a2**4*dx**6 &
-1.19d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+5.4d1*a1**3*a2**4*d2*dx**4+5.4d1*a1**4*a2**3*d2*dx**4 &
+3.6d2*a1**2*a2**4*dx**4+7.2d2*a1**3*a2**3*dx**4 &
+3.6d2*a1**4*a2**2*dx**4+3.6d1*a1**3*a2**4*d2**2*dx**2 &
+3.6d1*a1**4*a2**3*d2**2*dx**2-2.16d2*a1**2*a2**4*d2*dx**2 &
-4.32d2*a1**3*a2**3*d2*dx**2-2.16d2*a1**4*a2**2*d2*dx**2 &
-3.51d2*a1*a2**4*dx**2-1.053d3*a1**2*a2**3*dx**2 &
-1.053d3*a1**3*a2**2*dx**2-3.51d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.28d2*a1*a2**4*d2+6.84d2*a1**2*a2**3*d2 &
+6.84d2*a1**3*a2**2*d2+2.28d2*a1**4*a2*d2+3.3d1*a2**4+1.32d2*a1*a2**3 &
+1.98d2*a1**2*a2**2+1.32d2*a1**3*a2+3.3d1*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=2, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-1.31891402807852383d1*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.33d2*a1**3*a2**4*dy**6 &
-1.33d2*a1**4*a2**3*dy**6-1.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.33d2*a1**3*a2**4*dx**2*dy**4+1.33d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+7.8d1*a1**3*a2**4*d2*dy**4 &
+7.8d1*a1**4*a2**3*d2*dy**4+3.06d2*a1**2*a2**4*dy**4 &
+6.12d2*a1**3*a2**3*dy**4+3.06d2*a1**4*a2**2*dy**4 &
-1.4d1*a1**4*a2**4*d2*dx**4*dy**2+1.33d2*a1**3*a2**4*dx**4*dy**2 &
+1.33d2*a1**4*a2**3*dx**4*dy**2+2.4d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-2.68d2*a1**3*a2**4*d2*dx**2*dy**2-2.68d2*a1**4*a2**3*d2*dx**2*dy**2 &
+3.4d2*a1**2*a2**4*dx**2*dy**2+6.8d2*a1**3*a2**3*dx**2*dy**2 &
+3.4d2*a1**4*a2**2*dx**2*dy**2+2.4d1*a1**3*a2**4*d2**2*dy**2 &
+2.4d1*a1**4*a2**3*d2**2*dy**2-1.5d2*a1**2*a2**4*d2*dy**2 &
-3.0d2*a1**3*a2**3*d2*dy**2-1.5d2*a1**4*a2**2*d2*dy**2 &
-4.05d2*a1*a2**4*dy**2-1.215d3*a1**2*a2**3*dy**2 &
-1.215d3*a1**3*a2**2*dy**2-4.05d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.33d2*a1**3*a2**4*dx**6 &
-1.33d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+7.8d1*a1**3*a2**4*d2*dx**4+7.8d1*a1**4*a2**3*d2*dx**4 &
+3.06d2*a1**2*a2**4*dx**4+6.12d2*a1**3*a2**3*dx**4 &
+3.06d2*a1**4*a2**2*dx**4+2.4d1*a1**3*a2**4*d2**2*dx**2 &
+2.4d1*a1**4*a2**3*d2**2*dx**2-1.5d2*a1**2*a2**4*d2*dx**2 &
-3.0d2*a1**3*a2**3*d2*dx**2-1.5d2*a1**4*a2**2*d2*dx**2 &
-4.05d2*a1*a2**4*dx**2-1.215d3*a1**2*a2**3*dx**2 &
-1.215d3*a1**3*a2**2*dx**2-4.05d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+2.58d2*a1*a2**4*d2+7.74d2*a1**2*a2**3*d2 &
+7.74d2*a1**3*a2**2*d2+2.58d2*a1**4*a2*d2+7.8d1*a2**4+3.12d2*a1*a2**3 &
+4.68d2*a1**2*a2**2+3.12d2*a1**3*a2+7.8d1*a1**4)
                case (-3)
                  rlYlm_laplacian = &
-9.3261305305638875d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.33d2*a1**3*a2**4*dy**6 &
-1.33d2*a1**4*a2**3*dy**6-4.2d1*a1**4*a2**4*d2*dx**2*dy**4 &
+3.99d2*a1**3*a2**4*dx**2*dy**4+3.99d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+7.8d1*a1**3*a2**4*d2*dy**4 &
+7.8d1*a1**4*a2**3*d2*dy**4+3.06d2*a1**2*a2**4*dy**4 &
+6.12d2*a1**3*a2**3*dy**4+3.06d2*a1**4*a2**2*dy**4 &
-1.4d1*a1**4*a2**4*d2*dx**4*dy**2+1.33d2*a1**3*a2**4*dx**4*dy**2 &
+1.33d2*a1**4*a2**3*dx**4*dy**2+4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-5.64d2*a1**3*a2**4*d2*dx**2*dy**2-5.64d2*a1**4*a2**3*d2*dx**2*dy**2 &
+9.18d2*a1**2*a2**4*dx**2*dy**2+1.836d3*a1**3*a2**3*dx**2*dy**2 &
+9.18d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**3*a2**4*d2**2*dy**2 &
+3.6d1*a1**4*a2**3*d2**2*dy**2-2.88d2*a1**2*a2**4*d2*dy**2 &
-5.76d2*a1**3*a2**3*d2*dy**2-2.88d2*a1**4*a2**2*d2*dy**2 &
-1.35d2*a1*a2**4*dy**2-4.05d2*a1**2*a2**3*dy**2 &
-4.05d2*a1**3*a2**2*dy**2-1.35d2*a1**4*a2*dy**2 &
+4.2d1*a1**4*a2**4*d2*dx**6-3.99d2*a1**3*a2**4*dx**6 &
-3.99d2*a1**4*a2**3*dx**6-3.6d1*a1**4*a2**4*d2**2*dx**4 &
+3.18d2*a1**3*a2**4*d2*dx**4+3.18d2*a1**4*a2**3*d2*dx**4 &
+2.04d2*a1**2*a2**4*dx**4+4.08d2*a1**3*a2**3*dx**4 &
+2.04d2*a1**4*a2**2*dx**4+3.6d1*a1**3*a2**4*d2**2*dx**2 &
+3.6d1*a1**4*a2**3*d2**2*dx**2-2.88d2*a1**2*a2**4*d2*dx**2 &
-5.76d2*a1**3*a2**3*d2*dx**2-2.88d2*a1**4*a2**2*d2*dx**2 &
-1.35d2*a1*a2**4*dx**2-4.05d2*a1**2*a2**3*dx**2 &
-4.05d2*a1**3*a2**2*dx**2-1.35d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+3.0d2*a1*a2**4*d2+9.0d2*a1**2*a2**3*d2 &
+9.0d2*a1**3*a2**2*d2+3.0d2*a1**4*a2*d2-1.95d2*a2**4-7.8d2*a1*a2**3 &
-1.17d3*a1**2*a2**2-7.8d2*a1**3*a2-1.95d2*a1**4)*dz
                case (-2)
                  rlYlm_laplacian = &
-4.98502645567176383d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(dy+dx)*(dy-1.0d0*dx)*(9.8d1*a1**3*a2**3*d2*dy**4 &
-9.31d2*a1**2*a2**3*dy**4-9.31d2*a1**3*a2**2*dy**4 &
+1.96d2*a1**3*a2**3*d2*dx**2*dy**2-1.862d3*a1**2*a2**3*dx**2*dy**2 &
-1.862d3*a1**3*a2**2*dx**2*dy**2-1.68d2*a1**3*a2**3*d2**2*dy**2 &
+1.68d3*a1**2*a2**3*d2*dy**2+1.68d3*a1**3*a2**2*d2*dy**2 &
-7.14d2*a1*a2**3*dy**2-1.428d3*a1**2*a2**2*dy**2-7.14d2*a1**3*a2*dy**2 &
+9.8d1*a1**3*a2**3*d2*dx**4-9.31d2*a1**2*a2**3*dx**4 &
-9.31d2*a1**3*a2**2*dx**4-1.68d2*a1**3*a2**3*d2**2*dx**2 &
+1.68d3*a1**2*a2**3*d2*dx**2+1.68d3*a1**3*a2**2*d2*dx**2 &
-7.14d2*a1*a2**3*dx**2-1.428d3*a1**2*a2**2*dx**2-7.14d2*a1**3*a2*dx**2 &
+7.2d1*a1**3*a2**3*d2**3-7.8d2*a1**2*a2**3*d2**2 &
-7.8d2*a1**3*a2**2*d2**2+8.7d2*a1*a2**3*d2+1.74d3*a1**2*a2**2*d2 &
+8.7d2*a1**3*a2*d2-4.05d2*a2**3-1.215d3*a1*a2**2-1.215d3*a1**2*a2 &
-4.05d2*a1**3)
                case (-1)
                  rlYlm_laplacian = &
-3.52494601119984446d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(9.8d1*a1**4*a2**4*d2*dy**6-9.31d2*a1**3*a2**4*dy**6 &
-9.31d2*a1**4*a2**3*dy**6+9.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
-9.31d2*a1**3*a2**4*dx**2*dy**4-9.31d2*a1**4*a2**3*dx**2*dy**4 &
-1.4d2*a1**4*a2**4*d2**2*dy**4+1.302d3*a1**3*a2**4*d2*dy**4 &
+1.302d3*a1**4*a2**3*d2*dy**4+2.38d2*a1**2*a2**4*dy**4 &
+4.76d2*a1**3*a2**3*dy**4+2.38d2*a1**4*a2**2*dy**4 &
-9.8d1*a1**4*a2**4*d2*dx**4*dy**2+9.31d2*a1**3*a2**4*dx**4*dy**2 &
+9.31d2*a1**4*a2**3*dx**4*dy**2-1.96d2*a1**3*a2**4*d2*dx**2*dy**2 &
-1.96d2*a1**4*a2**3*d2*dx**2*dy**2+1.666d3*a1**2*a2**4*dx**2*dy**2 &
+3.332d3*a1**3*a2**3*dx**2*dy**2+1.666d3*a1**4*a2**2*dx**2*dy**2 &
+4.8d1*a1**4*a2**4*d2**3*dy**2-3.8d2*a1**3*a2**4*d2**2*dy**2 &
-3.8d2*a1**4*a2**3*d2**2*dy**2-7.36d2*a1**2*a2**4*d2*dy**2 &
-1.472d3*a1**3*a2**3*d2*dy**2-7.36d2*a1**4*a2**2*d2*dy**2 &
+6.75d2*a1*a2**4*dy**2+2.025d3*a1**2*a2**3*dy**2 &
+2.025d3*a1**3*a2**2*dy**2+6.75d2*a1**4*a2*dy**2 &
-9.8d1*a1**4*a2**4*d2*dx**6+9.31d2*a1**3*a2**4*dx**6 &
+9.31d2*a1**4*a2**3*dx**6+1.4d2*a1**4*a2**4*d2**2*dx**4 &
-1.498d3*a1**3*a2**4*d2*dx**4-1.498d3*a1**4*a2**3*d2*dx**4 &
+1.428d3*a1**2*a2**4*dx**4+2.856d3*a1**3*a2**3*dx**4 &
+1.428d3*a1**4*a2**2*dx**4-4.8d1*a1**4*a2**4*d2**3*dx**2 &
+6.6d2*a1**3*a2**4*d2**2*dx**2+6.6d2*a1**4*a2**3*d2**2*dx**2 &
-1.896d3*a1**2*a2**4*d2*dx**2-3.792d3*a1**3*a2**3*d2*dx**2 &
-1.896d3*a1**4*a2**2*d2*dx**2+1.215d3*a1*a2**4*dx**2 &
+3.645d3*a1**2*a2**3*dx**2+3.645d3*a1**3*a2**2*dx**2 &
+1.215d3*a1**4*a2*dx**2-4.8d1*a1**3*a2**4*d2**3 &
-4.8d1*a1**4*a2**3*d2**3+5.16d2*a1**2*a2**4*d2**2 &
+1.032d3*a1**3*a2**3*d2**2+5.16d2*a1**4*a2**2*d2**2-8.64d2*a1*a2**4*d2 &
-2.592d3*a1**2*a2**3*d2-2.592d3*a1**3*a2**2*d2-8.64d2*a1**4*a2*d2 &
+3.51d2*a2**4+1.404d3*a1*a2**3+2.106d3*a1**2*a2**2+1.404d3*a1**3*a2 &
+3.51d2*a1**4)*dz
                case (0)
                  rlYlm_laplacian =5.57342901225845299d-1*E*a1**2*a2**2*sqrt &
                  (a2+a1)*(a2+a1)**(-12)*(dy+dx)*(dy-1.0d0*dx)* &
(4.9d2*a1**4*a2**4*d2*dy**6-4.655d3*a1**3*a2**4*dy**6 &
-4.655d3*a1**4*a2**3*dy**6+1.47d3*a1**4*a2**4*d2*dx**2*dy**4 &
-1.3965d4*a1**3*a2**4*dx**2*dy**4-1.3965d4*a1**4*a2**3*dx**2*dy**4 &
-9.8d2*a1**4*a2**4*d2**2*dy**4+9.31d3*a1**3*a2**4*d2*dy**4 &
+9.31d3*a1**4*a2**3*d2*dy**4+1.47d3*a1**4*a2**4*d2*dx**4*dy**2 &
-1.3965d4*a1**3*a2**4*dx**4*dy**2-1.3965d4*a1**4*a2**3*dx**4*dy**2 &
-1.96d3*a1**4*a2**4*d2**2*dx**2*dy**2 &
+1.862d4*a1**3*a2**4*d2*dx**2*dy**2+1.862d4*a1**4*a2**3*d2*dx**2*dy**2 &
+5.92d2*a1**4*a2**4*d2**3*dy**2-5.48d3*a1**3*a2**4*d2**2*dy**2 &
-5.48d3*a1**4*a2**3*d2**2*dy**2-1.62d3*a1**2*a2**4*d2*dy**2 &
-3.24d3*a1**3*a2**3*d2*dy**2-1.62d3*a1**4*a2**2*d2*dy**2 &
+2.97d3*a1*a2**4*dy**2+8.91d3*a1**2*a2**3*dy**2 &
+8.91d3*a1**3*a2**2*dy**2+2.97d3*a1**4*a2*dy**2 &
+4.9d2*a1**4*a2**4*d2*dx**6-4.655d3*a1**3*a2**4*dx**6 &
-4.655d3*a1**4*a2**3*dx**6-9.8d2*a1**4*a2**4*d2**2*dx**4 &
+9.31d3*a1**3*a2**4*d2*dx**4+9.31d3*a1**4*a2**3*d2*dx**4 &
+5.92d2*a1**4*a2**4*d2**3*dx**2-5.48d3*a1**3*a2**4*d2**2*dx**2 &
-5.48d3*a1**4*a2**3*d2**2*dx**2-1.62d3*a1**2*a2**4*d2*dx**2 &
-3.24d3*a1**3*a2**3*d2*dx**2-1.62d3*a1**4*a2**2*d2*dx**2 &
+2.97d3*a1*a2**4*dx**2+8.91d3*a1**2*a2**3*dx**2 &
+8.91d3*a1**3*a2**2*dx**2+2.97d3*a1**4*a2*dx**2 &
-9.6d1*a1**4*a2**4*d2**4+7.2d2*a1**3*a2**4*d2**3 &
+7.2d2*a1**4*a2**3*d2**3+2.28d3*a1**2*a2**4*d2**2 &
+4.56d3*a1**3*a2**3*d2**2+2.28d3*a1**4*a2**2*d2**2-5.22d3*a1*a2**4*d2 &
-1.566d4*a1**2*a2**3*d2-1.566d4*a1**3*a2**2*d2-5.22d3*a1**4*a2*d2 &
+2.34d3*a2**4+9.36d3*a1*a2**3+1.404d4*a1**2*a2**2+9.36d3*a1**3*a2 &
+2.34d3*a1**4)
                case (1)
                  rlYlm_laplacian = &
-3.52494601119984446d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(9.8d1*a1**4*a2**4*d2*dy**6-9.31d2*a1**3*a2**4*dy**6 &
-9.31d2*a1**4*a2**3*dy**6+9.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
-9.31d2*a1**3*a2**4*dx**2*dy**4-9.31d2*a1**4*a2**3*dx**2*dy**4 &
-1.4d2*a1**4*a2**4*d2**2*dy**4+1.498d3*a1**3*a2**4*d2*dy**4 &
+1.498d3*a1**4*a2**3*d2*dy**4-1.428d3*a1**2*a2**4*dy**4 &
-2.856d3*a1**3*a2**3*dy**4-1.428d3*a1**4*a2**2*dy**4 &
-9.8d1*a1**4*a2**4*d2*dx**4*dy**2+9.31d2*a1**3*a2**4*dx**4*dy**2 &
+9.31d2*a1**4*a2**3*dx**4*dy**2+1.96d2*a1**3*a2**4*d2*dx**2*dy**2 &
+1.96d2*a1**4*a2**3*d2*dx**2*dy**2-1.666d3*a1**2*a2**4*dx**2*dy**2 &
-3.332d3*a1**3*a2**3*dx**2*dy**2-1.666d3*a1**4*a2**2*dx**2*dy**2 &
+4.8d1*a1**4*a2**4*d2**3*dy**2-6.6d2*a1**3*a2**4*d2**2*dy**2 &
-6.6d2*a1**4*a2**3*d2**2*dy**2+1.896d3*a1**2*a2**4*d2*dy**2 &
+3.792d3*a1**3*a2**3*d2*dy**2+1.896d3*a1**4*a2**2*d2*dy**2 &
-1.215d3*a1*a2**4*dy**2-3.645d3*a1**2*a2**3*dy**2 &
-3.645d3*a1**3*a2**2*dy**2-1.215d3*a1**4*a2*dy**2 &
-9.8d1*a1**4*a2**4*d2*dx**6+9.31d2*a1**3*a2**4*dx**6 &
+9.31d2*a1**4*a2**3*dx**6+1.4d2*a1**4*a2**4*d2**2*dx**4 &
-1.302d3*a1**3*a2**4*d2*dx**4-1.302d3*a1**4*a2**3*d2*dx**4 &
-2.38d2*a1**2*a2**4*dx**4-4.76d2*a1**3*a2**3*dx**4 &
-2.38d2*a1**4*a2**2*dx**4-4.8d1*a1**4*a2**4*d2**3*dx**2 &
+3.8d2*a1**3*a2**4*d2**2*dx**2+3.8d2*a1**4*a2**3*d2**2*dx**2 &
+7.36d2*a1**2*a2**4*d2*dx**2+1.472d3*a1**3*a2**3*d2*dx**2 &
+7.36d2*a1**4*a2**2*d2*dx**2-6.75d2*a1*a2**4*dx**2 &
-2.025d3*a1**2*a2**3*dx**2-2.025d3*a1**3*a2**2*dx**2 &
-6.75d2*a1**4*a2*dx**2+4.8d1*a1**3*a2**4*d2**3+4.8d1*a1**4*a2**3*d2**3 &
-5.16d2*a1**2*a2**4*d2**2-1.032d3*a1**3*a2**3*d2**2 &
-5.16d2*a1**4*a2**2*d2**2+8.64d2*a1*a2**4*d2+2.592d3*a1**2*a2**3*d2 &
+2.592d3*a1**3*a2**2*d2+8.64d2*a1**4*a2*d2-3.51d2*a2**4 &
-1.404d3*a1*a2**3-2.106d3*a1**2*a2**2-1.404d3*a1**3*a2-3.51d2*a1**4 &
)*dz
                case (2)
                  rlYlm_laplacian =2.49251322783588191d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(9.8d1*a1**5*a2**5*d2*dy**8 &
-9.31d2*a1**4*a2**5*dy**8-9.31d2*a1**5*a2**4*dy**8 &
-1.68d2*a1**5*a2**5*d2**2*dy**6+1.484d3*a1**4*a2**5*d2*dy**6 &
+1.484d3*a1**5*a2**4*d2*dy**6+9.52d2*a1**3*a2**5*dy**6 &
+1.904d3*a1**4*a2**4*dy**6+9.52d2*a1**5*a2**3*dy**6 &
-1.96d2*a1**5*a2**5*d2*dx**4*dy**4+1.862d3*a1**4*a2**5*dx**4*dy**4 &
+1.862d3*a1**5*a2**4*dx**4*dy**4+1.68d2*a1**5*a2**5*d2**2*dx**2*dy**4 &
-2.268d3*a1**4*a2**5*d2*dx**2*dy**4-2.268d3*a1**5*a2**4*d2*dx**2*dy**4 &
+5.712d3*a1**3*a2**5*dx**2*dy**4+1.1424d4*a1**4*a2**4*dx**2*dy**4 &
+5.712d3*a1**5*a2**3*dx**2*dy**4+7.2d1*a1**5*a2**5*d2**3*dy**4 &
-4.44d2*a1**4*a2**5*d2**2*dy**4-4.44d2*a1**5*a2**4*d2**2*dy**4 &
-2.112d3*a1**3*a2**5*d2*dy**4-4.224d3*a1**4*a2**4*d2*dy**4 &
-2.112d3*a1**5*a2**3*d2*dy**4+5.4d2*a1**2*a2**5*dy**4 &
+1.62d3*a1**3*a2**4*dy**4+1.62d3*a1**4*a2**3*dy**4 &
+5.4d2*a1**5*a2**2*dy**4+1.68d2*a1**5*a2**5*d2**2*dx**4*dy**2 &
-2.268d3*a1**4*a2**5*d2*dx**4*dy**2-2.268d3*a1**5*a2**4*d2*dx**4*dy**2 &
+5.712d3*a1**3*a2**5*dx**4*dy**2+1.1424d4*a1**4*a2**4*dx**4*dy**2 &
+5.712d3*a1**5*a2**3*dx**4*dy**2-1.44d2*a1**5*a2**5*d2**3*dx**2*dy**2 &
+2.232d3*a1**4*a2**5*d2**2*dx**2*dy**2 &
+2.232d3*a1**5*a2**4*d2**2*dx**2*dy**2 &
-7.704d3*a1**3*a2**5*d2*dx**2*dy**2 &
-1.5408d4*a1**4*a2**4*d2*dx**2*dy**2 &
-7.704d3*a1**5*a2**3*d2*dx**2*dy**2+2.7d3*a1**2*a2**5*dx**2*dy**2 &
+8.1d3*a1**3*a2**4*dx**2*dy**2+8.1d3*a1**4*a2**3*dx**2*dy**2 &
+2.7d3*a1**5*a2**2*dx**2*dy**2-1.44d2*a1**4*a2**5*d2**3*dy**2 &
-1.44d2*a1**5*a2**4*d2**3*dy**2+1.296d3*a1**3*a2**5*d2**2*dy**2 &
+2.592d3*a1**4*a2**4*d2**2*dy**2+1.296d3*a1**5*a2**3*d2**2*dy**2 &
-4.92d2*a1**2*a2**5*d2*dy**2-1.476d3*a1**3*a2**4*d2*dy**2 &
-1.476d3*a1**4*a2**3*d2*dy**2-4.92d2*a1**5*a2**2*d2*dy**2 &
-3.12d2*a1*a2**5*dy**2-1.248d3*a1**2*a2**4*dy**2 &
-1.872d3*a1**3*a2**3*dy**2-1.248d3*a1**4*a2**2*dy**2 &
-3.12d2*a1**5*a2*dy**2+9.8d1*a1**5*a2**5*d2*dx**8 &
-9.31d2*a1**4*a2**5*dx**8-9.31d2*a1**5*a2**4*dx**8 &
-1.68d2*a1**5*a2**5*d2**2*dx**6+1.484d3*a1**4*a2**5*d2*dx**6 &
+1.484d3*a1**5*a2**4*d2*dx**6+9.52d2*a1**3*a2**5*dx**6 &
+1.904d3*a1**4*a2**4*dx**6+9.52d2*a1**5*a2**3*dx**6 &
+7.2d1*a1**5*a2**5*d2**3*dx**4-4.44d2*a1**4*a2**5*d2**2*dx**4 &
-4.44d2*a1**5*a2**4*d2**2*dx**4-2.112d3*a1**3*a2**5*d2*dx**4 &
-4.224d3*a1**4*a2**4*d2*dx**4-2.112d3*a1**5*a2**3*d2*dx**4 &
+5.4d2*a1**2*a2**5*dx**4+1.62d3*a1**3*a2**4*dx**4 &
+1.62d3*a1**4*a2**3*dx**4+5.4d2*a1**5*a2**2*dx**4 &
-1.44d2*a1**4*a2**5*d2**3*dx**2-1.44d2*a1**5*a2**4*d2**3*dx**2 &
+1.296d3*a1**3*a2**5*d2**2*dx**2+2.592d3*a1**4*a2**4*d2**2*dx**2 &
+1.296d3*a1**5*a2**3*d2**2*dx**2-4.92d2*a1**2*a2**5*d2*dx**2 &
-1.476d3*a1**3*a2**4*d2*dx**2-1.476d3*a1**4*a2**3*d2*dx**2 &
-4.92d2*a1**5*a2**2*d2*dx**2-3.12d2*a1*a2**5*dx**2 &
-1.248d3*a1**2*a2**4*dx**2-1.872d3*a1**3*a2**3*dx**2 &
-1.248d3*a1**4*a2**2*dx**2-3.12d2*a1**5*a2*dx**2 &
+7.2d1*a1**3*a2**5*d2**3+1.44d2*a1**4*a2**4*d2**3 &
+7.2d1*a1**5*a2**3*d2**3-6.84d2*a1**2*a2**5*d2**2 &
-2.052d3*a1**3*a2**4*d2**2-2.052d3*a1**4*a2**3*d2**2 &
-6.84d2*a1**5*a2**2*d2**2+9.78d2*a1*a2**5*d2+3.912d3*a1**2*a2**4*d2 &
+5.868d3*a1**3*a2**3*d2+3.912d3*a1**4*a2**2*d2+9.78d2*a1**5*a2*d2 &
-2.31d2*a2**5-1.155d3*a1*a2**4-2.31d3*a1**2*a2**3-2.31d3*a1**3*a2**2 &
-1.155d3*a1**4*a2-2.31d2*a1**5)
                case (3)
                  rlYlm_laplacian = &
-9.3261305305638875d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(4.2d1*a1**4*a2**4*d2*dy**6-3.99d2*a1**3*a2**4*dy**6 &
-3.99d2*a1**4*a2**3*dy**6-1.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.33d2*a1**3*a2**4*dx**2*dy**4+1.33d2*a1**4*a2**3*dx**2*dy**4 &
-3.6d1*a1**4*a2**4*d2**2*dy**4+3.18d2*a1**3*a2**4*d2*dy**4 &
+3.18d2*a1**4*a2**3*d2*dy**4+2.04d2*a1**2*a2**4*dy**4 &
+4.08d2*a1**3*a2**3*dy**4+2.04d2*a1**4*a2**2*dy**4 &
-4.2d1*a1**4*a2**4*d2*dx**4*dy**2+3.99d2*a1**3*a2**4*dx**4*dy**2 &
+3.99d2*a1**4*a2**3*dx**4*dy**2+4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-5.64d2*a1**3*a2**4*d2*dx**2*dy**2-5.64d2*a1**4*a2**3*d2*dx**2*dy**2 &
+9.18d2*a1**2*a2**4*dx**2*dy**2+1.836d3*a1**3*a2**3*dx**2*dy**2 &
+9.18d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**3*a2**4*d2**2*dy**2 &
+3.6d1*a1**4*a2**3*d2**2*dy**2-2.88d2*a1**2*a2**4*d2*dy**2 &
-5.76d2*a1**3*a2**3*d2*dy**2-2.88d2*a1**4*a2**2*d2*dy**2 &
-1.35d2*a1*a2**4*dy**2-4.05d2*a1**2*a2**3*dy**2 &
-4.05d2*a1**3*a2**2*dy**2-1.35d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.33d2*a1**3*a2**4*dx**6 &
-1.33d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+7.8d1*a1**3*a2**4*d2*dx**4+7.8d1*a1**4*a2**3*d2*dx**4 &
+3.06d2*a1**2*a2**4*dx**4+6.12d2*a1**3*a2**3*dx**4 &
+3.06d2*a1**4*a2**2*dx**4+3.6d1*a1**3*a2**4*d2**2*dx**2 &
+3.6d1*a1**4*a2**3*d2**2*dx**2-2.88d2*a1**2*a2**4*d2*dx**2 &
-5.76d2*a1**3*a2**3*d2*dx**2-2.88d2*a1**4*a2**2*d2*dx**2 &
-1.35d2*a1*a2**4*dx**2-4.05d2*a1**2*a2**3*dx**2 &
-4.05d2*a1**3*a2**2*dx**2-1.35d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+3.0d2*a1*a2**4*d2+9.0d2*a1**2*a2**3*d2 &
+9.0d2*a1**3*a2**2*d2+3.0d2*a1**4*a2*d2-1.95d2*a2**4-7.8d2*a1*a2**3 &
-1.17d3*a1**2*a2**2-7.8d2*a1**3*a2-1.95d2*a1**4)*dz
                case (4)
                  rlYlm_laplacian =3.29728507019630958d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.33d2*a1**3*a2**4*dy**6 &
-1.33d2*a1**4*a2**3*dy**6-7.0d1*a1**4*a2**4*d2*dx**2*dy**4 &
+6.65d2*a1**3*a2**4*dx**2*dy**4+6.65d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+5.0d1*a1**3*a2**4*d2*dy**4 &
+5.0d1*a1**4*a2**3*d2*dy**4+5.44d2*a1**2*a2**4*dy**4 &
+1.088d3*a1**3*a2**3*dy**4+5.44d2*a1**4*a2**2*dy**4 &
-7.0d1*a1**4*a2**4*d2*dx**4*dy**2+6.65d2*a1**3*a2**4*dx**4*dy**2 &
+6.65d2*a1**4*a2**3*dx**4*dy**2+7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-7.48d2*a1**3*a2**4*d2*dx**2*dy**2-7.48d2*a1**4*a2**3*d2*dx**2*dy**2 &
+5.44d2*a1**2*a2**4*dx**2*dy**2+1.088d3*a1**3*a2**3*dx**2*dy**2 &
+5.44d2*a1**4*a2**2*dx**2*dy**2+4.8d1*a1**3*a2**4*d2**2*dy**2 &
+4.8d1*a1**4*a2**3*d2**2*dy**2-3.0d2*a1**2*a2**4*d2*dy**2 &
-6.0d2*a1**3*a2**3*d2*dy**2-3.0d2*a1**4*a2**2*d2*dy**2 &
-8.1d2*a1*a2**4*dy**2-2.43d3*a1**2*a2**3*dy**2 &
-2.43d3*a1**3*a2**2*dy**2-8.1d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.33d2*a1**3*a2**4*dx**6 &
-1.33d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+5.0d1*a1**3*a2**4*d2*dx**4+5.0d1*a1**4*a2**3*d2*dx**4 &
+5.44d2*a1**2*a2**4*dx**4+1.088d3*a1**3*a2**3*dx**4 &
+5.44d2*a1**4*a2**2*dx**4+4.8d1*a1**3*a2**4*d2**2*dx**2 &
+4.8d1*a1**4*a2**3*d2**2*dx**2-3.0d2*a1**2*a2**4*d2*dx**2 &
-6.0d2*a1**3*a2**3*d2*dx**2-3.0d2*a1**4*a2**2*d2*dx**2 &
-8.1d2*a1*a2**4*dx**2-2.43d3*a1**2*a2**3*dx**2 &
-2.43d3*a1**3*a2**2*dx**2-8.1d2*a1**4*a2*dx**2-7.2d1*a1**2*a2**4*d2**2 &
-1.44d2*a1**3*a2**3*d2**2-7.2d1*a1**4*a2**2*d2**2+5.16d2*a1*a2**4*d2 &
+1.548d3*a1**2*a2**3*d2+1.548d3*a1**3*a2**2*d2+5.16d2*a1**4*a2*d2 &
+1.56d2*a2**4+6.24d2*a1*a2**3+9.36d2*a1**2*a2**2+6.24d2*a1**3*a2 &
+1.56d2*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (3)
          ! selection on l2: l1=4, m1=3
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=3, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian = &
-5.56102982223387534d0*E*a1*a2**5*sqrt(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2* &
(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*dx*(3.0d0*dy**2-1.0d0*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=3, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian = &
-9.63198619451479393d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)*dx*dy* &
(6.0d0*a1**2*a2**2*d2*dy**2-3.9d1*a1*a2**2*dy**2-3.9d1*a1**2*a2*dy**2 &
-2.0d0*a1**2*a2**2*d2*dx**2+1.3d1*a1*a2**2*dx**2+1.3d1*a1**2*a2*dx**2 &
-6.0d0*a1*a2**2*d2-6.0d0*a1**2*a2*d2+3.3d1*a2**2+6.6d1*a1*a2 &
+3.3d1*a1**2)*dz
                case (0)
                  rlYlm_laplacian =4.81599309725739696d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(3.0d0*dy**2-1.0d0*dx**2)* &
(4.0d0*a1**2*a2**2*d2*dy**2-2.6d1*a1*a2**2*dy**2-2.6d1*a1**2*a2*dy**2 &
+4.0d0*a1**2*a2**2*d2*dx**2-2.6d1*a1*a2**2*dx**2-2.6d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.8d1*a1*a2**2*d2+2.8d1*a1**2*a2*d2 &
-1.1d1*a2**2-2.2d1*a1*a2-1.1d1*a1**2)
                case (1)
                  rlYlm_laplacian = &
-4.81599309725739696d0*E*a1*a2**4*sqrt(a2+a1)*(a2+a1)**(-9)* &
(1.2d1*a1**2*a2**2*d2*dx**2*dy**2-7.8d1*a1*a2**2*dx**2*dy**2 &
-7.8d1*a1**2*a2*dx**2*dy**2-6.0d0*a1*a2**2*d2*dy**2 &
-6.0d0*a1**2*a2*d2*dy**2+3.3d1*a2**2*dy**2+6.6d1*a1*a2*dy**2 &
+3.3d1*a1**2*dy**2-4.0d0*a1**2*a2**2*d2*dx**4+2.6d1*a1*a2**2*dx**4 &
+2.6d1*a1**2*a2*dx**4+6.0d0*a1*a2**2*d2*dx**2+6.0d0*a1**2*a2*d2*dx**2 &
-3.3d1*a2**2*dx**2-6.6d1*a1*a2*dx**2-3.3d1*a1**2*dx**2)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=3, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian = &
-1.07688879446372956d1*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*dy* &
(1.2d1*a1**3*a2**3*d2*dx**2*dy**2-9.0d1*a1**2*a2**3*dx**2*dy**2 &
-9.0d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+3.0d1*a1**2*a2**3*dx**4 &
+3.0d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian =1.07688879446372956d1*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(1.2d1*a1**3*a2**3*d2*dy**4 &
-9.0d1*a1**2*a2**3*dy**4-9.0d1*a1**3*a2**2*dy**4 &
+8.0d0*a1**3*a2**3*d2*dx**2*dy**2-6.0d1*a1**2*a2**3*dx**2*dy**2 &
-6.0d1*a1**3*a2**2*dx**2*dy**2-1.2d1*a1**3*a2**3*d2**2*dy**2 &
+8.4d1*a1**2*a2**3*d2*dy**2+8.4d1*a1**3*a2**2*d2*dy**2 &
+3.9d1*a1*a2**3*dy**2+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
-4.0d0*a1**3*a2**3*d2*dx**4+3.0d1*a1**2*a2**3*dx**4 &
+3.0d1*a1**3*a2**2*dx**4+4.0d0*a1**3*a2**3*d2**2*dx**2 &
-4.4d1*a1**2*a2**3*d2*dx**2-4.4d1*a1**3*a2**2*d2*dx**2 &
+9.1d1*a1*a2**3*dx**2+1.82d2*a1**2*a2**2*dx**2+9.1d1*a1**3*a2*dx**2 &
+1.2d1*a1**2*a2**3*d2**2+1.2d1*a1**3*a2**2*d2**2-8.4d1*a1*a2**3*d2 &
-1.68d2*a1**2*a2**2*d2-8.4d1*a1**3*a2*d2+3.3d1*a2**3+9.9d1*a1*a2**2 &
+9.9d1*a1**2*a2+3.3d1*a1**3)
                case (0)
                  rlYlm_laplacian =6.21742035370925833d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(3.0d0*dy**2-1.0d0*dx**2)* &
(6.0d0*a1**2*a2**2*d2*dy**2-4.5d1*a1*a2**2*dy**2-4.5d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-4.5d1*a1*a2**2*dx**2-4.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.8d1*a1*a2**2*d2+2.8d1*a1**2*a2*d2 &
+1.3d1*a2**2+2.6d1*a1*a2+1.3d1*a1**2)*dz
                case (1)
                  rlYlm_laplacian =5.3844439723186478d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*(2.4d1*a1**3*a2**3*d2*dx**2*dy**4 &
-1.8d2*a1**2*a2**3*dx**2*dy**4-1.8d2*a1**3*a2**2*dx**2*dy**4 &
-1.2d1*a1**2*a2**3*d2*dy**4-1.2d1*a1**3*a2**2*d2*dy**4 &
+7.8d1*a1*a2**3*dy**4+1.56d2*a1**2*a2**2*dy**4+7.8d1*a1**3*a2*dy**4 &
+1.6d1*a1**3*a2**3*d2*dx**4*dy**2-1.2d2*a1**2*a2**3*dx**4*dy**2 &
-1.2d2*a1**3*a2**2*dx**4*dy**2-2.4d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+1.92d2*a1**2*a2**3*d2*dx**2*dy**2+1.92d2*a1**3*a2**2*d2*dx**2*dy**2 &
-7.8d1*a1*a2**3*dx**2*dy**2-1.56d2*a1**2*a2**2*dx**2*dy**2 &
-7.8d1*a1**3*a2*dx**2*dy**2+1.2d1*a1**2*a2**3*d2**2*dy**2 &
+1.2d1*a1**3*a2**2*d2**2*dy**2-8.4d1*a1*a2**3*d2*dy**2 &
-1.68d2*a1**2*a2**2*d2*dy**2-8.4d1*a1**3*a2*d2*dy**2+3.3d1*a2**3*dy**2 &
+9.9d1*a1*a2**2*dy**2+9.9d1*a1**2*a2*dy**2+3.3d1*a1**3*dy**2 &
-8.0d0*a1**3*a2**3*d2*dx**6+6.0d1*a1**2*a2**3*dx**6 &
+6.0d1*a1**3*a2**2*dx**6+8.0d0*a1**3*a2**3*d2**2*dx**4 &
-5.2d1*a1**2*a2**3*d2*dx**4-5.2d1*a1**3*a2**2*d2*dx**4 &
-5.2d1*a1*a2**3*dx**4-1.04d2*a1**2*a2**2*dx**4-5.2d1*a1**3*a2*dx**4 &
-1.2d1*a1**2*a2**3*d2**2*dx**2-1.2d1*a1**3*a2**2*d2**2*dx**2 &
+8.4d1*a1*a2**3*d2*dx**2+1.68d2*a1**2*a2**2*d2*dx**2 &
+8.4d1*a1**3*a2*d2*dx**2-3.3d1*a2**3*dx**2-9.9d1*a1*a2**2*dx**2 &
-9.9d1*a1**2*a2*dx**2-3.3d1*a1**3*dx**2)
                case (2)
                  rlYlm_laplacian =1.07688879446372956d1*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(6.0d0*a1**3*a2**3*d2*dy**4 &
-4.5d1*a1**2*a2**3*dy**4-4.5d1*a1**3*a2**2*dy**4 &
-8.0d0*a1**3*a2**3*d2*dx**2*dy**2+6.0d1*a1**2*a2**3*dx**2*dy**2 &
+6.0d1*a1**3*a2**2*dx**2*dy**2-6.0d0*a1**2*a2**3*d2*dy**2 &
-6.0d0*a1**3*a2**2*d2*dy**2+3.9d1*a1*a2**3*dy**2 &
+7.8d1*a1**2*a2**2*dy**2+3.9d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4-6.0d0*a1**2*a2**3*d2*dx**2 &
-6.0d0*a1**3*a2**2*d2*dx**2+3.9d1*a1*a2**3*dx**2 &
+7.8d1*a1**2*a2**2*dx**2+3.9d1*a1**3*a2*dx**2+6.0d0*a1*a2**3*d2 &
+1.2d1*a1**2*a2**2*d2+6.0d0*a1**3*a2*d2-3.3d1*a2**3-9.9d1*a1*a2**2 &
-9.9d1*a1**2*a2-3.3d1*a1**3)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=3, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian =1.16317283965674489d1*E*a1**4*a2**5*sqrt &
(a2+a1)*(a2+a1)**(-11)*(1.0d0*a2*(2.0d0*a1*d2-1.7d1)-1.7d1*a1)*dx*dy* &
(dy**2-3.0d0*dx**2)*(3.0d0*dy**2-1.0d0*dx**2)*dz
                case (-2)
                  rlYlm_laplacian =1.42458996991158946d1*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*dy*(2.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
-2.04d2*a1**3*a2**4*dx**2*dy**4-2.04d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**3*a2**4*d2*dy**4-1.2d1*a1**4*a2**3*d2*dy**4 &
+9.0d1*a1**2*a2**4*dy**4+1.8d2*a1**3*a2**3*dy**4 &
+9.0d1*a1**4*a2**2*dy**4+1.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
-1.36d2*a1**3*a2**4*dx**4*dy**2-1.36d2*a1**4*a2**3*dx**4*dy**2 &
-2.4d1*a1**4*a2**4*d2**2*dx**2*dy**2+1.92d2*a1**3*a2**4*d2*dx**2*dy**2 &
+1.92d2*a1**4*a2**3*d2*dx**2*dy**2+9.0d1*a1**2*a2**4*dx**2*dy**2 &
+1.8d2*a1**3*a2**3*dx**2*dy**2+9.0d1*a1**4*a2**2*dx**2*dy**2 &
+1.2d1*a1**3*a2**4*d2**2*dy**2+1.2d1*a1**4*a2**3*d2**2*dy**2 &
-8.4d1*a1**2*a2**4*d2*dy**2-1.68d2*a1**3*a2**3*d2*dy**2 &
-8.4d1*a1**4*a2**2*d2*dy**2-3.9d1*a1*a2**4*dy**2 &
-1.17d2*a1**2*a2**3*dy**2-1.17d2*a1**3*a2**2*dy**2 &
-3.9d1*a1**4*a2*dy**2-8.0d0*a1**4*a2**4*d2*dx**6 &
+6.8d1*a1**3*a2**4*dx**6+6.8d1*a1**4*a2**3*dx**6 &
+8.0d0*a1**4*a2**4*d2**2*dx**4-8.4d1*a1**3*a2**4*d2*dx**4 &
-8.4d1*a1**4*a2**3*d2*dx**4+1.2d2*a1**2*a2**4*dx**4 &
+2.4d2*a1**3*a2**3*dx**4+1.2d2*a1**4*a2**2*dx**4 &
+1.2d1*a1**3*a2**4*d2**2*dx**2+1.2d1*a1**4*a2**3*d2**2*dx**2 &
-8.4d1*a1**2*a2**4*d2*dx**2-1.68d2*a1**3*a2**3*d2*dx**2 &
-8.4d1*a1**4*a2**2*d2*dx**2-3.9d1*a1*a2**4*dx**2 &
-1.17d2*a1**2*a2**3*dx**2-1.17d2*a1**3*a2**2*dx**2 &
-3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2-2.4d1*a1**3*a2**3*d2**2 &
-1.2d1*a1**4*a2**2*d2**2+8.4d1*a1*a2**4*d2+2.52d2*a1**2*a2**3*d2 &
+2.52d2*a1**3*a2**2*d2+8.4d1*a1**4*a2*d2-3.3d1*a2**4-1.32d2*a1*a2**3 &
-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2-3.3d1*a1**4)
                case (-1)
                  rlYlm_laplacian =9.00989807350272604d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(3.0d1*a1**3*a2**3*d2*dy**4 &
-2.55d2*a1**2*a2**3*dy**4-2.55d2*a1**3*a2**2*dy**4 &
+2.0d1*a1**3*a2**3*d2*dx**2*dy**2-1.7d2*a1**2*a2**3*dx**2*dy**2 &
-1.7d2*a1**3*a2**2*dx**2*dy**2-2.4d1*a1**3*a2**3*d2**2*dy**2 &
+1.8d2*a1**2*a2**3*d2*dy**2+1.8d2*a1**3*a2**2*d2*dy**2 &
+1.8d2*a1*a2**3*dy**2+3.6d2*a1**2*a2**2*dy**2+1.8d2*a1**3*a2*dy**2 &
-1.0d1*a1**3*a2**3*d2*dx**4+8.5d1*a1**2*a2**3*dx**4 &
+8.5d1*a1**3*a2**2*dx**4+8.0d0*a1**3*a2**3*d2**2*dx**2 &
-1.0d2*a1**2*a2**3*d2*dx**2-1.0d2*a1**3*a2**2*d2*dx**2 &
+2.4d2*a1*a2**3*dx**2+4.8d2*a1**2*a2**2*dx**2+2.4d2*a1**3*a2*dx**2 &
+2.4d1*a1**2*a2**3*d2**2+2.4d1*a1**3*a2**2*d2**2-1.92d2*a1*a2**3*d2 &
-3.84d2*a1**2*a2**2*d2-1.92d2*a1**3*a2*d2+7.8d1*a2**3+2.34d2*a1*a2**2 &
+2.34d2*a1**2*a2+7.8d1*a1**3)*dz
                case (0)
                  rlYlm_laplacian = &
-3.67827548576114071d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*(2.0d1*a1**3*a2**3*d2*dy**4 &
-1.7d2*a1**2*a2**3*dy**4-1.7d2*a1**3*a2**2*dy**4 &
+4.0d1*a1**3*a2**3*d2*dx**2*dy**2-3.4d2*a1**2*a2**3*dx**2*dy**2 &
-3.4d2*a1**3*a2**2*dx**2*dy**2-2.8d1*a1**3*a2**3*d2**2*dy**2 &
+2.2d2*a1**2*a2**3*d2*dy**2+2.2d2*a1**3*a2**2*d2*dy**2 &
+1.35d2*a1*a2**3*dy**2+2.7d2*a1**2*a2**2*dy**2+1.35d2*a1**3*a2*dy**2 &
+2.0d1*a1**3*a2**3*d2*dx**4-1.7d2*a1**2*a2**3*dx**4 &
-1.7d2*a1**3*a2**2*dx**4-2.8d1*a1**3*a2**3*d2**2*dx**2 &
+2.2d2*a1**2*a2**3*d2*dx**2+2.2d2*a1**3*a2**2*d2*dx**2 &
+1.35d2*a1*a2**3*dx**2+2.7d2*a1**2*a2**2*dx**2+1.35d2*a1**3*a2*dx**2 &
+8.0d0*a1**3*a2**3*d2**3-4.4d1*a1**2*a2**3*d2**2 &
-4.4d1*a1**3*a2**2*d2**2-1.98d2*a1*a2**3*d2-3.96d2*a1**2*a2**2*d2 &
-1.98d2*a1**3*a2*d2+1.17d2*a2**3+3.51d2*a1*a2**2+3.51d2*a1**2*a2 &
+1.17d2*a1**3)
                case (1)
                  rlYlm_laplacian =4.50494903675136302d0*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*(6.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
-5.1d2*a1**2*a2**3*dx**2*dy**4-5.1d2*a1**3*a2**2*dx**2*dy**4 &
-3.0d1*a1**2*a2**3*d2*dy**4-3.0d1*a1**3*a2**2*d2*dy**4 &
+2.25d2*a1*a2**3*dy**4+4.5d2*a1**2*a2**2*dy**4+2.25d2*a1**3*a2*dy**4 &
+4.0d1*a1**3*a2**3*d2*dx**4*dy**2-3.4d2*a1**2*a2**3*dx**4*dy**2 &
-3.4d2*a1**3*a2**2*dx**4*dy**2-4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
+4.2d2*a1**2*a2**3*d2*dx**2*dy**2+4.2d2*a1**3*a2**2*d2*dx**2*dy**2 &
-9.0d1*a1*a2**3*dx**2*dy**2-1.8d2*a1**2*a2**2*dx**2*dy**2 &
-9.0d1*a1**3*a2*dx**2*dy**2+2.4d1*a1**2*a2**3*d2**2*dy**2 &
+2.4d1*a1**3*a2**2*d2**2*dy**2-1.92d2*a1*a2**3*d2*dy**2 &
-3.84d2*a1**2*a2**2*d2*dy**2-1.92d2*a1**3*a2*d2*dy**2 &
+7.8d1*a2**3*dy**2+2.34d2*a1*a2**2*dy**2+2.34d2*a1**2*a2*dy**2 &
+7.8d1*a1**3*dy**2-2.0d1*a1**3*a2**3*d2*dx**6+1.7d2*a1**2*a2**3*dx**6 &
+1.7d2*a1**3*a2**2*dx**6+1.6d1*a1**3*a2**3*d2**2*dx**4 &
-1.1d2*a1**2*a2**3*d2*dx**4-1.1d2*a1**3*a2**2*d2*dx**4 &
-1.95d2*a1*a2**3*dx**4-3.9d2*a1**2*a2**2*dx**4-1.95d2*a1**3*a2*dx**4 &
-2.4d1*a1**2*a2**3*d2**2*dx**2-2.4d1*a1**3*a2**2*d2**2*dx**2 &
+1.92d2*a1*a2**3*d2*dx**2+3.84d2*a1**2*a2**2*d2*dx**2 &
+1.92d2*a1**3*a2*d2*dx**2-7.8d1*a2**3*dx**2-2.34d2*a1*a2**2*dx**2 &
-2.34d2*a1**2*a2*dx**2-7.8d1*a1**3*dx**2)*dz
                case (2)
                  rlYlm_laplacian = &
-1.42458996991158946d1*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(1.2d1*a1**4*a2**4*d2*dy**6-1.02d2*a1**3*a2**4*dy**6 &
-1.02d2*a1**4*a2**3*dy**6-4.0d0*a1**4*a2**4*d2*dx**2*dy**4 &
+3.4d1*a1**3*a2**4*dx**2*dy**4+3.4d1*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+9.6d1*a1**3*a2**4*d2*dy**4 &
+9.6d1*a1**4*a2**3*d2*dy**4+4.5d1*a1**2*a2**4*dy**4 &
+9.0d1*a1**3*a2**3*dy**4+4.5d1*a1**4*a2**2*dy**4 &
-1.2d1*a1**4*a2**4*d2*dx**4*dy**2+1.02d2*a1**3*a2**4*dx**4*dy**2 &
+1.02d2*a1**4*a2**3*dx**4*dy**2+1.6d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-1.68d2*a1**3*a2**4*d2*dx**2*dy**2-1.68d2*a1**4*a2**3*d2*dx**2*dy**2 &
+2.4d2*a1**2*a2**4*dx**2*dy**2+4.8d2*a1**3*a2**3*dx**2*dy**2 &
+2.4d2*a1**4*a2**2*dx**2*dy**2+1.2d1*a1**3*a2**4*d2**2*dy**2 &
+1.2d1*a1**4*a2**3*d2**2*dy**2-8.4d1*a1**2*a2**4*d2*dy**2 &
-1.68d2*a1**3*a2**3*d2*dy**2-8.4d1*a1**4*a2**2*d2*dy**2 &
-3.9d1*a1*a2**4*dy**2-1.17d2*a1**2*a2**3*dy**2 &
-1.17d2*a1**3*a2**2*dy**2-3.9d1*a1**4*a2*dy**2 &
+4.0d0*a1**4*a2**4*d2*dx**6-3.4d1*a1**3*a2**4*dx**6 &
-3.4d1*a1**4*a2**3*dx**6-4.0d0*a1**4*a2**4*d2**2*dx**4 &
+2.4d1*a1**3*a2**4*d2*dx**4+2.4d1*a1**4*a2**3*d2*dx**4 &
+7.5d1*a1**2*a2**4*dx**4+1.5d2*a1**3*a2**3*dx**4 &
+7.5d1*a1**4*a2**2*dx**4+1.2d1*a1**3*a2**4*d2**2*dx**2 &
+1.2d1*a1**4*a2**3*d2**2*dx**2-8.4d1*a1**2*a2**4*d2*dx**2 &
-1.68d2*a1**3*a2**3*d2*dx**2-8.4d1*a1**4*a2**2*d2*dx**2 &
-3.9d1*a1*a2**4*dx**2-1.17d2*a1**2*a2**3*dx**2 &
-1.17d2*a1**3*a2**2*dx**2-3.9d1*a1**4*a2*dx**2-1.2d1*a1**2*a2**4*d2**2 &
-2.4d1*a1**3*a2**3*d2**2-1.2d1*a1**4*a2**2*d2**2+8.4d1*a1*a2**4*d2 &
+2.52d2*a1**2*a2**3*d2+2.52d2*a1**3*a2**2*d2+8.4d1*a1**4*a2*d2 &
-3.3d1*a2**4-1.32d2*a1*a2**3-1.98d2*a1**2*a2**2-1.32d2*a1**3*a2 &
-3.3d1*a1**4)
                case (3)
                  rlYlm_laplacian =5.81586419828372447d0*E*a1*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-11)*(3.6d1*a1**4*a2**4*d2*dx**2*dy**4 &
-3.06d2*a1**3*a2**4*dx**2*dy**4-3.06d2*a1**4*a2**3*dx**2*dy**4 &
-1.8d1*a1**3*a2**4*d2*dy**4-1.8d1*a1**4*a2**3*d2*dy**4 &
+1.35d2*a1**2*a2**4*dy**4+2.7d2*a1**3*a2**3*dy**4 &
+1.35d2*a1**4*a2**2*dy**4-2.4d1*a1**4*a2**4*d2*dx**4*dy**2 &
+2.04d2*a1**3*a2**4*dx**4*dy**2+2.04d2*a1**4*a2**3*dx**4*dy**2 &
-3.6d1*a1**3*a2**4*d2*dx**2*dy**2-3.6d1*a1**4*a2**3*d2*dx**2*dy**2 &
+2.7d2*a1**2*a2**4*dx**2*dy**2+5.4d2*a1**3*a2**3*dx**2*dy**2 &
+2.7d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**2*a2**4*d2*dy**2 &
+7.2d1*a1**3*a2**3*d2*dy**2+3.6d1*a1**4*a2**2*d2*dy**2 &
-2.34d2*a1*a2**4*dy**2-7.02d2*a1**2*a2**3*dy**2 &
-7.02d2*a1**3*a2**2*dy**2-2.34d2*a1**4*a2*dy**2 &
+4.0d0*a1**4*a2**4*d2*dx**6-3.4d1*a1**3*a2**4*dx**6 &
-3.4d1*a1**4*a2**3*dx**6-1.8d1*a1**3*a2**4*d2*dx**4 &
-1.8d1*a1**4*a2**3*d2*dx**4+1.35d2*a1**2*a2**4*dx**4 &
+2.7d2*a1**3*a2**3*dx**4+1.35d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**2*a2**4*d2*dx**2+7.2d1*a1**3*a2**3*d2*dx**2 &
+3.6d1*a1**4*a2**2*d2*dx**2-2.34d2*a1*a2**4*dx**2 &
-7.02d2*a1**2*a2**3*dx**2-7.02d2*a1**3*a2**2*dx**2 &
-2.34d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+6.6d1*a2**4+2.64d2*a1*a2**3 &
+3.96d2*a1**2*a2**2+2.64d2*a1**3*a2+6.6d1*a1**4)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=3, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian =2.46746220783989112d1*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dy*(1.2d1*a1**4*a2**4*d2*dx**2*dy**4 &
-1.14d2*a1**3*a2**4*dx**2*dy**4-1.14d2*a1**4*a2**3*dx**2*dy**4 &
-6.0d0*a1**3*a2**4*d2*dy**4-6.0d0*a1**4*a2**3*d2*dy**4 &
+5.1d1*a1**2*a2**4*dy**4+1.02d2*a1**3*a2**3*dy**4 &
+5.1d1*a1**4*a2**2*dy**4-1.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
+1.52d2*a1**3*a2**4*dx**4*dy**2+1.52d2*a1**4*a2**3*dx**4*dy**2 &
-1.2d1*a1**3*a2**4*d2*dx**2*dy**2-1.2d1*a1**4*a2**3*d2*dx**2*dy**2 &
+1.02d2*a1**2*a2**4*dx**2*dy**2+2.04d2*a1**3*a2**3*dx**2*dy**2 &
+1.02d2*a1**4*a2**2*dx**2*dy**2+1.8d1*a1**2*a2**4*d2*dy**2 &
+3.6d1*a1**3*a2**3*d2*dy**2+1.8d1*a1**4*a2**2*d2*dy**2 &
-1.35d2*a1*a2**4*dy**2-4.05d2*a1**2*a2**3*dy**2 &
-4.05d2*a1**3*a2**2*dy**2-1.35d2*a1**4*a2*dy**2 &
+4.0d0*a1**4*a2**4*d2*dx**6-3.8d1*a1**3*a2**4*dx**6 &
-3.8d1*a1**4*a2**3*dx**6-6.0d0*a1**3*a2**4*d2*dx**4 &
-6.0d0*a1**4*a2**3*d2*dx**4+5.1d1*a1**2*a2**4*dx**4 &
+1.02d2*a1**3*a2**3*dx**4+5.1d1*a1**4*a2**2*dx**4 &
+1.8d1*a1**2*a2**4*d2*dx**2+3.6d1*a1**3*a2**3*d2*dx**2 &
+1.8d1*a1**4*a2**2*d2*dx**2-1.35d2*a1*a2**4*dx**2 &
-4.05d2*a1**2*a2**3*dx**2-4.05d2*a1**3*a2**2*dx**2 &
-1.35d2*a1**4*a2*dx**2-1.2d1*a1*a2**4*d2-3.6d1*a1**2*a2**3*d2 &
-3.6d1*a1**3*a2**2*d2-1.2d1*a1**4*a2*d2+7.8d1*a2**4+3.12d2*a1*a2**3 &
+4.68d2*a1**2*a2**2+3.12d2*a1**3*a2+7.8d1*a1**4)*dz
                case (-3)
                  rlYlm_laplacian = &
-1.74475925948511734d1*E*a1**4*a2**4*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(dy**2-3.0d0*dx**2)*(3.0d0*dy**2-1.0d0*dx**2)* &
(4.0d0*a1**2*a2**2*d2*dy**2-3.8d1*a1*a2**2*dy**2-3.8d1*a1**2*a2*dy**2 &
+4.0d0*a1**2*a2**2*d2*dx**2-3.8d1*a1*a2**2*dx**2-3.8d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+4.0d1*a1*a2**2*d2+4.0d1*a1**2*a2*d2 &
-1.7d1*a2**2-3.4d1*a1*a2-1.7d1*a1**2)
                case (-2)
                  rlYlm_laplacian =9.3261305305638875d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*dy*(8.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
-7.98d2*a1**3*a2**4*dx**2*dy**4-7.98d2*a1**4*a2**3*dx**2*dy**4 &
-4.2d1*a1**3*a2**4*d2*dy**4-4.2d1*a1**4*a2**3*d2*dy**4 &
+3.57d2*a1**2*a2**4*dy**4+7.14d2*a1**3*a2**3*dy**4 &
+3.57d2*a1**4*a2**2*dy**4+5.6d1*a1**4*a2**4*d2*dx**4*dy**2 &
-5.32d2*a1**3*a2**4*dx**4*dy**2-5.32d2*a1**4*a2**3*dx**4*dy**2 &
-7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2+6.36d2*a1**3*a2**4*d2*dx**2*dy**2 &
+6.36d2*a1**4*a2**3*d2*dx**2*dy**2+4.08d2*a1**2*a2**4*dx**2*dy**2 &
+8.16d2*a1**3*a2**3*dx**2*dy**2+4.08d2*a1**4*a2**2*dx**2*dy**2 &
+3.6d1*a1**3*a2**4*d2**2*dy**2+3.6d1*a1**4*a2**3*d2**2*dy**2 &
-2.88d2*a1**2*a2**4*d2*dy**2-5.76d2*a1**3*a2**3*d2*dy**2 &
-2.88d2*a1**4*a2**2*d2*dy**2-1.35d2*a1*a2**4*dy**2 &
-4.05d2*a1**2*a2**3*dy**2-4.05d2*a1**3*a2**2*dy**2 &
-1.35d2*a1**4*a2*dy**2-2.8d1*a1**4*a2**4*d2*dx**6 &
+2.66d2*a1**3*a2**4*dx**6+2.66d2*a1**4*a2**3*dx**6 &
+2.4d1*a1**4*a2**4*d2**2*dx**4-2.82d2*a1**3*a2**4*d2*dx**4 &
-2.82d2*a1**4*a2**3*d2*dx**4+4.59d2*a1**2*a2**4*dx**4 &
+9.18d2*a1**3*a2**3*dx**4+4.59d2*a1**4*a2**2*dx**4 &
+3.6d1*a1**3*a2**4*d2**2*dx**2+3.6d1*a1**4*a2**3*d2**2*dx**2 &
-2.88d2*a1**2*a2**4*d2*dx**2-5.76d2*a1**3*a2**3*d2*dx**2 &
-2.88d2*a1**4*a2**2*d2*dx**2-1.35d2*a1*a2**4*dx**2 &
-4.05d2*a1**2*a2**3*dx**2-4.05d2*a1**3*a2**2*dx**2 &
-1.35d2*a1**4*a2*dx**2-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+3.0d2*a1*a2**4*d2+9.0d2*a1**2*a2**3*d2 &
+9.0d2*a1**3*a2**2*d2+3.0d2*a1**4*a2*d2-1.95d2*a2**4-7.8d2*a1*a2**3 &
-1.17d3*a1**2*a2**2-7.8d2*a1**3*a2-1.95d2*a1**4)*dz
                case (-1)
                  rlYlm_laplacian = &
-6.59457014039261917d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(8.4d1*a1**4*a2**4*d2*dy**6-7.98d2*a1**3*a2**4*dy**6 &
-7.98d2*a1**4*a2**3*dy**6+1.4d2*a1**4*a2**4*d2*dx**2*dy**4 &
-1.33d3*a1**3*a2**4*dx**2*dy**4-1.33d3*a1**4*a2**3*dx**2*dy**4 &
-1.32d2*a1**4*a2**4*d2**2*dy**4+1.152d3*a1**3*a2**4*d2*dy**4 &
+1.152d3*a1**4*a2**3*d2*dy**4+8.67d2*a1**2*a2**4*dy**4 &
+1.734d3*a1**3*a2**3*dy**4+8.67d2*a1**4*a2**2*dy**4 &
+2.8d1*a1**4*a2**4*d2*dx**4*dy**2-2.66d2*a1**3*a2**4*dx**4*dy**2 &
-2.66d2*a1**4*a2**3*dx**4*dy**2-8.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
+6.56d2*a1**3*a2**4*d2*dx**2*dy**2+6.56d2*a1**4*a2**3*d2*dx**2*dy**2 &
+1.53d3*a1**2*a2**4*dx**2*dy**2+3.06d3*a1**3*a2**3*dx**2*dy**2 &
+1.53d3*a1**4*a2**2*dx**2*dy**2+4.8d1*a1**4*a2**4*d2**3*dy**2 &
-2.88d2*a1**3*a2**4*d2**2*dy**2-2.88d2*a1**4*a2**3*d2**2*dy**2 &
-1.5d3*a1**2*a2**4*d2*dy**2-3.0d3*a1**3*a2**3*d2*dy**2 &
-1.5d3*a1**4*a2**2*d2*dy**2+5.4d2*a1*a2**4*dy**2 &
+1.62d3*a1**2*a2**3*dy**2+1.62d3*a1**3*a2**2*dy**2 &
+5.4d2*a1**4*a2*dy**2-2.8d1*a1**4*a2**4*d2*dx**6 &
+2.66d2*a1**3*a2**4*dx**6+2.66d2*a1**4*a2**3*dx**6 &
+4.4d1*a1**4*a2**4*d2**2*dx**4-4.96d2*a1**3*a2**4*d2*dx**4 &
-4.96d2*a1**4*a2**3*d2*dx**4+6.63d2*a1**2*a2**4*dx**4 &
+1.326d3*a1**3*a2**3*dx**4+6.63d2*a1**4*a2**2*dx**4 &
-1.6d1*a1**4*a2**4*d2**3*dx**2+2.72d2*a1**3*a2**4*d2**2*dx**2 &
+2.72d2*a1**4*a2**3*d2**2*dx**2-1.02d3*a1**2*a2**4*d2*dx**2 &
-2.04d3*a1**3*a2**3*d2*dx**2-1.02d3*a1**4*a2**2*d2*dx**2 &
-4.8d1*a1**3*a2**4*d2**3-4.8d1*a1**4*a2**3*d2**3 &
+4.08d2*a1**2*a2**4*d2**2+8.16d2*a1**3*a2**3*d2**2 &
+4.08d2*a1**4*a2**2*d2**2+3.6d1*a1*a2**4*d2+1.08d2*a1**2*a2**3*d2 &
+1.08d2*a1**3*a2**2*d2+3.6d1*a1**4*a2*d2-2.34d2*a2**4-9.36d2*a1*a2**3 &
-1.404d3*a1**2*a2**2-9.36d2*a1**3*a2-2.34d2*a1**4)
                case (0)
                  rlYlm_laplacian = &
-2.08538618333770325d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(3.0d0*dy**2-1.0d0*dx**2)*(7.0d1*a1**3*a2**3*d2*dy**4 &
-6.65d2*a1**2*a2**3*dy**4-6.65d2*a1**3*a2**2*dy**4 &
+1.4d2*a1**3*a2**3*d2*dx**2*dy**2-1.33d3*a1**2*a2**3*dx**2*dy**2 &
-1.33d3*a1**3*a2**2*dx**2*dy**2-8.0d1*a1**3*a2**3*d2**2*dy**2 &
+6.6d2*a1**2*a2**3*d2*dy**2+6.6d2*a1**3*a2**2*d2*dy**2 &
+8.5d2*a1*a2**3*dy**2+1.7d3*a1**2*a2**2*dy**2+8.5d2*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-6.65d2*a1**2*a2**3*dx**4 &
-6.65d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+6.6d2*a1**2*a2**3*d2*dx**2+6.6d2*a1**3*a2**2*d2*dx**2 &
+8.5d2*a1*a2**3*dx**2+1.7d3*a1**2*a2**2*dx**2+8.5d2*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3-4.0d1*a1**2*a2**3*d2**2 &
-4.0d1*a1**3*a2**2*d2**2-1.06d3*a1*a2**3*d2-2.12d3*a1**2*a2**2*d2 &
-1.06d3*a1**3*a2*d2+8.1d2*a2**3+2.43d3*a1*a2**2+2.43d3*a1**2*a2 &
+8.1d2*a1**3)*dz
                case (1)
                  rlYlm_laplacian = &
-3.29728507019630958d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)* &
(1.68d2*a1**4*a2**4*d2*dx**2*dy**6-1.596d3*a1**3*a2**4*dx**2*dy**6 &
-1.596d3*a1**4*a2**3*dx**2*dy**6-8.4d1*a1**3*a2**4*d2*dy**6 &
-8.4d1*a1**4*a2**3*d2*dy**6+7.14d2*a1**2*a2**4*dy**6 &
+1.428d3*a1**3*a2**3*dy**6+7.14d2*a1**4*a2**2*dy**6 &
+2.8d2*a1**4*a2**4*d2*dx**4*dy**4-2.66d3*a1**3*a2**4*dx**4*dy**4 &
-2.66d3*a1**4*a2**3*dx**4*dy**4-2.64d2*a1**4*a2**4*d2**2*dx**2*dy**4 &
+2.388d3*a1**3*a2**4*d2*dx**2*dy**4+2.388d3*a1**4*a2**3*d2*dx**2*dy**4 &
+1.02d3*a1**2*a2**4*dx**2*dy**4+2.04d3*a1**3*a2**3*dx**2*dy**4 &
+1.02d3*a1**4*a2**2*dx**2*dy**4+1.32d2*a1**3*a2**4*d2**2*dy**4 &
+1.32d2*a1**4*a2**3*d2**2*dy**4-1.14d3*a1**2*a2**4*d2*dy**4 &
-2.28d3*a1**3*a2**3*d2*dy**4-1.14d3*a1**4*a2**2*d2*dy**4 &
+1.35d2*a1*a2**4*dy**4+4.05d2*a1**2*a2**3*dy**4 &
+4.05d2*a1**3*a2**2*dy**4+1.35d2*a1**4*a2*dy**4 &
+5.6d1*a1**4*a2**4*d2*dx**6*dy**2-5.32d2*a1**3*a2**4*dx**6*dy**2 &
-5.32d2*a1**4*a2**3*dx**6*dy**2-1.76d2*a1**4*a2**4*d2**2*dx**4*dy**2 &
+1.732d3*a1**3*a2**4*d2*dx**4*dy**2+1.732d3*a1**4*a2**3*d2*dx**4*dy**2 &
-5.1d2*a1**2*a2**4*dx**4*dy**2-1.02d3*a1**3*a2**3*dx**4*dy**2 &
-5.1d2*a1**4*a2**2*dx**4*dy**2+9.6d1*a1**4*a2**4*d2**3*dx**2*dy**2 &
-8.4d2*a1**3*a2**4*d2**2*dx**2*dy**2 &
-8.4d2*a1**4*a2**3*d2**2*dx**2*dy**2-7.2d2*a1**2*a2**4*d2*dx**2*dy**2 &
-1.44d3*a1**3*a2**3*d2*dx**2*dy**2-7.2d2*a1**4*a2**2*d2*dx**2*dy**2 &
+8.1d2*a1*a2**4*dx**2*dy**2+2.43d3*a1**2*a2**3*dx**2*dy**2 &
+2.43d3*a1**3*a2**2*dx**2*dy**2+8.1d2*a1**4*a2*dx**2*dy**2 &
-4.8d1*a1**3*a2**4*d2**3*dy**2-4.8d1*a1**4*a2**3*d2**3*dy**2 &
+4.08d2*a1**2*a2**4*d2**2*dy**2+8.16d2*a1**3*a2**3*d2**2*dy**2 &
+4.08d2*a1**4*a2**2*d2**2*dy**2+3.6d1*a1*a2**4*d2*dy**2 &
+1.08d2*a1**2*a2**3*d2*dy**2+1.08d2*a1**3*a2**2*d2*dy**2 &
+3.6d1*a1**4*a2*d2*dy**2-2.34d2*a2**4*dy**2-9.36d2*a1*a2**3*dy**2 &
-1.404d3*a1**2*a2**2*dy**2-9.36d2*a1**3*a2*dy**2-2.34d2*a1**4*dy**2 &
-5.6d1*a1**4*a2**4*d2*dx**8+5.32d2*a1**3*a2**4*dx**8 &
+5.32d2*a1**4*a2**3*dx**8+8.8d1*a1**4*a2**4*d2**2*dx**6 &
-7.4d2*a1**3*a2**4*d2*dx**6-7.4d2*a1**4*a2**3*d2*dx**6 &
-8.16d2*a1**2*a2**4*dx**6-1.632d3*a1**3*a2**3*dx**6 &
-8.16d2*a1**4*a2**2*dx**6-3.2d1*a1**4*a2**4*d2**3*dx**4 &
+1.48d2*a1**3*a2**4*d2**2*dx**4+1.48d2*a1**4*a2**3*d2**2*dx**4 &
+1.38d3*a1**2*a2**4*d2*dx**4+2.76d3*a1**3*a2**3*d2*dx**4 &
+1.38d3*a1**4*a2**2*d2*dx**4-4.05d2*a1*a2**4*dx**4 &
-1.215d3*a1**2*a2**3*dx**4-1.215d3*a1**3*a2**2*dx**4 &
-4.05d2*a1**4*a2*dx**4+4.8d1*a1**3*a2**4*d2**3*dx**2 &
+4.8d1*a1**4*a2**3*d2**3*dx**2-4.08d2*a1**2*a2**4*d2**2*dx**2 &
-8.16d2*a1**3*a2**3*d2**2*dx**2-4.08d2*a1**4*a2**2*d2**2*dx**2 &
-3.6d1*a1*a2**4*d2*dx**2-1.08d2*a1**2*a2**3*d2*dx**2 &
-1.08d2*a1**3*a2**2*d2*dx**2-3.6d1*a1**4*a2*d2*dx**2 &
+2.34d2*a2**4*dx**2+9.36d2*a1*a2**3*dx**2+1.404d3*a1**2*a2**2*dx**2 &
+9.36d2*a1**3*a2*dx**2+2.34d2*a1**4*dx**2)
                case (2)
                  rlYlm_laplacian = &
-9.3261305305638875d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(4.2d1*a1**4*a2**4*d2*dy**6-3.99d2*a1**3*a2**4*dy**6 &
-3.99d2*a1**4*a2**3*dy**6-1.4d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.33d2*a1**3*a2**4*dx**2*dy**4+1.33d2*a1**4*a2**3*dx**2*dy**4 &
-3.6d1*a1**4*a2**4*d2**2*dy**4+3.18d2*a1**3*a2**4*d2*dy**4 &
+3.18d2*a1**4*a2**3*d2*dy**4+2.04d2*a1**2*a2**4*dy**4 &
+4.08d2*a1**3*a2**3*dy**4+2.04d2*a1**4*a2**2*dy**4 &
-4.2d1*a1**4*a2**4*d2*dx**4*dy**2+3.99d2*a1**3*a2**4*dx**4*dy**2 &
+3.99d2*a1**4*a2**3*dx**4*dy**2+4.8d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-5.64d2*a1**3*a2**4*d2*dx**2*dy**2-5.64d2*a1**4*a2**3*d2*dx**2*dy**2 &
+9.18d2*a1**2*a2**4*dx**2*dy**2+1.836d3*a1**3*a2**3*dx**2*dy**2 &
+9.18d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**3*a2**4*d2**2*dy**2 &
+3.6d1*a1**4*a2**3*d2**2*dy**2-2.88d2*a1**2*a2**4*d2*dy**2 &
-5.76d2*a1**3*a2**3*d2*dy**2-2.88d2*a1**4*a2**2*d2*dy**2 &
-1.35d2*a1*a2**4*dy**2-4.05d2*a1**2*a2**3*dy**2 &
-4.05d2*a1**3*a2**2*dy**2-1.35d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.33d2*a1**3*a2**4*dx**6 &
-1.33d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+7.8d1*a1**3*a2**4*d2*dx**4+7.8d1*a1**4*a2**3*d2*dx**4 &
+3.06d2*a1**2*a2**4*dx**4+6.12d2*a1**3*a2**3*dx**4 &
+3.06d2*a1**4*a2**2*dx**4+3.6d1*a1**3*a2**4*d2**2*dx**2 &
+3.6d1*a1**4*a2**3*d2**2*dx**2-2.88d2*a1**2*a2**4*d2*dx**2 &
-5.76d2*a1**3*a2**3*d2*dx**2-2.88d2*a1**4*a2**2*d2*dx**2 &
-1.35d2*a1*a2**4*dx**2-4.05d2*a1**2*a2**3*dx**2 &
-4.05d2*a1**3*a2**2*dx**2-1.35d2*a1**4*a2*dx**2 &
-3.6d1*a1**2*a2**4*d2**2-7.2d1*a1**3*a2**3*d2**2 &
-3.6d1*a1**4*a2**2*d2**2+3.0d2*a1*a2**4*d2+9.0d2*a1**2*a2**3*d2 &
+9.0d2*a1**3*a2**2*d2+3.0d2*a1**4*a2*d2-1.95d2*a2**4-7.8d2*a1*a2**3 &
-1.17d3*a1**2*a2**2-7.8d2*a1**3*a2-1.95d2*a1**4)*dz
                case (3)
                  rlYlm_laplacian =-8.7237962974255867d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(7.2d1*a1**5*a2**5*d2*dx**2*dy**6 &
-6.84d2*a1**4*a2**5*dx**2*dy**6-6.84d2*a1**5*a2**4*dx**2*dy**6 &
-3.6d1*a1**4*a2**5*d2*dy**6-3.6d1*a1**5*a2**4*d2*dy**6 &
+3.06d2*a1**3*a2**5*dy**6+6.12d2*a1**4*a2**4*dy**6 &
+3.06d2*a1**5*a2**3*dy**6+2.4d1*a1**5*a2**5*d2*dx**4*dy**4 &
-2.28d2*a1**4*a2**5*dx**4*dy**4-2.28d2*a1**5*a2**4*dx**4*dy**4 &
-7.2d1*a1**5*a2**5*d2**2*dx**2*dy**4+6.12d2*a1**4*a2**5*d2*dx**2*dy**4 &
+6.12d2*a1**5*a2**4*d2*dx**2*dy**4+6.12d2*a1**3*a2**5*dx**2*dy**4 &
+1.224d3*a1**4*a2**4*dx**2*dy**4+6.12d2*a1**5*a2**3*dx**2*dy**4 &
+3.6d1*a1**4*a2**5*d2**2*dy**4+3.6d1*a1**5*a2**4*d2**2*dy**4 &
-2.52d2*a1**3*a2**5*d2*dy**4-5.04d2*a1**4*a2**4*d2*dy**4 &
-2.52d2*a1**5*a2**3*d2*dy**4-4.05d2*a1**2*a2**5*dy**4 &
-1.215d3*a1**3*a2**4*dy**4-1.215d3*a1**4*a2**3*dy**4 &
-4.05d2*a1**5*a2**2*dy**4-4.0d1*a1**5*a2**5*d2*dx**6*dy**2 &
+3.8d2*a1**4*a2**5*dx**6*dy**2+3.8d2*a1**5*a2**4*dx**6*dy**2 &
+4.8d1*a1**5*a2**5*d2**2*dx**4*dy**2-5.88d2*a1**4*a2**5*d2*dx**4*dy**2 &
-5.88d2*a1**5*a2**4*d2*dx**4*dy**2+1.122d3*a1**3*a2**5*dx**4*dy**2 &
+2.244d3*a1**4*a2**4*dx**4*dy**2+1.122d3*a1**5*a2**3*dx**4*dy**2 &
+7.2d1*a1**4*a2**5*d2**2*dx**2*dy**2 &
+7.2d1*a1**5*a2**4*d2**2*dx**2*dy**2-5.04d2*a1**3*a2**5*d2*dx**2*dy**2 &
-1.008d3*a1**4*a2**4*d2*dx**2*dy**2-5.04d2*a1**5*a2**3*d2*dx**2*dy**2 &
-8.1d2*a1**2*a2**5*dx**2*dy**2-2.43d3*a1**3*a2**4*dx**2*dy**2 &
-2.43d3*a1**4*a2**3*dx**2*dy**2-8.1d2*a1**5*a2**2*dx**2*dy**2 &
-7.2d1*a1**3*a2**5*d2**2*dy**2-1.44d2*a1**4*a2**4*d2**2*dy**2 &
-7.2d1*a1**5*a2**3*d2**2*dy**2+5.52d2*a1**2*a2**5*d2*dy**2 &
+1.656d3*a1**3*a2**4*d2*dy**2+1.656d3*a1**4*a2**3*d2*dy**2 &
+5.52d2*a1**5*a2**2*d2*dy**2-7.8d1*a1*a2**5*dy**2 &
-3.12d2*a1**2*a2**4*dy**2-4.68d2*a1**3*a2**3*dy**2 &
-3.12d2*a1**4*a2**2*dy**2-7.8d1*a1**5*a2*dy**2 &
+8.0d0*a1**5*a2**5*d2*dx**8-7.6d1*a1**4*a2**5*dx**8 &
-7.6d1*a1**5*a2**4*dx**8-8.0d0*a1**5*a2**5*d2**2*dx**6 &
+4.4d1*a1**4*a2**5*d2*dx**6+4.4d1*a1**5*a2**4*d2*dx**6 &
+2.72d2*a1**3*a2**5*dx**6+5.44d2*a1**4*a2**4*dx**6 &
+2.72d2*a1**5*a2**3*dx**6+3.6d1*a1**4*a2**5*d2**2*dx**4 &
+3.6d1*a1**5*a2**4*d2**2*dx**4-2.52d2*a1**3*a2**5*d2*dx**4 &
-5.04d2*a1**4*a2**4*d2*dx**4-2.52d2*a1**5*a2**3*d2*dx**4 &
-4.05d2*a1**2*a2**5*dx**4-1.215d3*a1**3*a2**4*dx**4 &
-1.215d3*a1**4*a2**3*dx**4-4.05d2*a1**5*a2**2*dx**4 &
-7.2d1*a1**3*a2**5*d2**2*dx**2-1.44d2*a1**4*a2**4*d2**2*dx**2 &
-7.2d1*a1**5*a2**3*d2**2*dx**2+5.52d2*a1**2*a2**5*d2*dx**2 &
+1.656d3*a1**3*a2**4*d2*dx**2+1.656d3*a1**4*a2**3*d2*dx**2 &
+5.52d2*a1**5*a2**2*d2*dx**2-7.8d1*a1*a2**5*dx**2 &
-3.12d2*a1**2*a2**4*dx**2-4.68d2*a1**3*a2**3*dx**2 &
-3.12d2*a1**4*a2**2*dx**2-7.8d1*a1**5*a2*dx**2+2.4d1*a1**2*a2**5*d2**2 &
+7.2d1*a1**3*a2**4*d2**2+7.2d1*a1**4*a2**3*d2**2 &
+2.4d1*a1**5*a2**2*d2**2-1.68d2*a1*a2**5*d2-6.72d2*a1**2*a2**4*d2 &
-1.008d3*a1**3*a2**3*d2-6.72d2*a1**4*a2**2*d2-1.68d2*a1**5*a2*d2 &
+6.6d1*a2**5+3.3d2*a1*a2**4+6.6d2*a1**2*a2**3+6.6d2*a1**3*a2**2 &
+3.3d2*a1**4*a2+6.6d1*a1**5)
                case (4)
                  rlYlm_laplacian = &
-1.23373110391994556d1*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(6.0d0*a1**4*a2**4*d2*dy**6-5.7d1*a1**3*a2**4*dy**6 &
-5.7d1*a1**4*a2**3*dy**6-3.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
+3.61d2*a1**3*a2**4*dx**2*dy**4+3.61d2*a1**4*a2**3*dx**2*dy**4 &
+1.2d1*a1**3*a2**4*d2*dy**4+1.2d1*a1**4*a2**3*d2*dy**4 &
-1.02d2*a1**2*a2**4*dy**4-2.04d2*a1**3*a2**3*dy**4 &
-1.02d2*a1**4*a2**2*dy**4+1.8d1*a1**4*a2**4*d2*dx**4*dy**2 &
-1.71d2*a1**3*a2**4*dx**4*dy**2-1.71d2*a1**4*a2**3*dx**4*dy**2 &
+2.4d1*a1**3*a2**4*d2*dx**2*dy**2+2.4d1*a1**4*a2**3*d2*dx**2*dy**2 &
-2.04d2*a1**2*a2**4*dx**2*dy**2-4.08d2*a1**3*a2**3*dx**2*dy**2 &
-2.04d2*a1**4*a2**2*dx**2*dy**2-3.6d1*a1**2*a2**4*d2*dy**2 &
-7.2d1*a1**3*a2**3*d2*dy**2-3.6d1*a1**4*a2**2*d2*dy**2 &
+2.7d2*a1*a2**4*dy**2+8.1d2*a1**2*a2**3*dy**2+8.1d2*a1**3*a2**2*dy**2 &
+2.7d2*a1**4*a2*dy**2-2.0d0*a1**4*a2**4*d2*dx**6 &
+1.9d1*a1**3*a2**4*dx**6+1.9d1*a1**4*a2**3*dx**6 &
+1.2d1*a1**3*a2**4*d2*dx**4+1.2d1*a1**4*a2**3*d2*dx**4 &
-1.02d2*a1**2*a2**4*dx**4-2.04d2*a1**3*a2**3*dx**4 &
-1.02d2*a1**4*a2**2*dx**4-3.6d1*a1**2*a2**4*d2*dx**2 &
-7.2d1*a1**3*a2**3*d2*dx**2-3.6d1*a1**4*a2**2*d2*dx**2 &
+2.7d2*a1*a2**4*dx**2+8.1d2*a1**2*a2**3*dx**2+8.1d2*a1**3*a2**2*dx**2 &
+2.7d2*a1**4*a2*dx**2+2.4d1*a1*a2**4*d2+7.2d1*a1**2*a2**3*d2 &
+7.2d1*a1**3*a2**2*d2+2.4d1*a1**4*a2*d2-1.56d2*a2**4-6.24d2*a1*a2**3 &
-9.36d2*a1**2*a2**2-6.24d2*a1**3*a2-1.56d2*a1**4)*dz
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case (4)
          ! selection on l2: l1=4, m1=4
          select case (l2)
            case (0)
              ! selection on m2: l1=4, m1=4, l2=0
              select case (m2)
                case (0)
                  rlYlm_laplacian =1.96612094884109708d0*E*a1*a2**5*sqrt &
(a2+a1)*(a2+a1)**(-8)*(1.0d0*a2*(2.0d0*a1*d2-1.1d1)-1.1d1*a1)*(dy**2 &
-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (1)
              ! selection on m2: l1=4, m1=4, l2=1
              select case (m2)
                case (-1)
                  rlYlm_laplacian =3.40542137721830949d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*dy*(2.0d0*a1**2*a2**2*d2*dy**4 &
-1.3d1*a1*a2**2*dy**4-1.3d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+7.8d1*a1*a2**2*dx**2*dy**2 &
+7.8d1*a1**2*a2*dx**2*dy**2-4.0d0*a1*a2**2*d2*dy**2 &
-4.0d0*a1**2*a2*d2*dy**2+2.2d1*a2**2*dy**2+4.4d1*a1*a2*dy**2 &
+2.2d1*a1**2*dy**2+2.0d0*a1**2*a2**2*d2*dx**4-1.3d1*a1*a2**2*dx**4 &
-1.3d1*a1**2*a2*dx**4+1.2d1*a1*a2**2*d2*dx**2+1.2d1*a1**2*a2*d2*dx**2 &
-6.6d1*a2**2*dx**2-1.32d2*a1*a2*dx**2-6.6d1*a1**2*dx**2)
                case (0)
                  rlYlm_laplacian =3.40542137721830949d0*E*a1**2*a2**5*sqrt &
(a2+a1)*(a2+a1)**(-9)*(1.0d0*a2*(2.0d0*a1*d2-1.3d1)-1.3d1*a1)*(dy**2 &
-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)*dz
                case (1)
                  rlYlm_laplacian =3.40542137721830949d0*E*a1*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-9)*dx*(2.0d0*a1**2*a2**2*d2*dy**4 &
-1.3d1*a1*a2**2*dy**4-1.3d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+7.8d1*a1*a2**2*dx**2*dy**2 &
+7.8d1*a1**2*a2*dx**2*dy**2+1.2d1*a1*a2**2*d2*dy**2 &
+1.2d1*a1**2*a2*d2*dy**2-6.6d1*a2**2*dy**2-1.32d2*a1*a2*dy**2 &
-6.6d1*a1**2*dy**2+2.0d0*a1**2*a2**2*d2*dx**4-1.3d1*a1*a2**2*dx**4 &
-1.3d1*a1**2*a2*dx**4-4.0d0*a1*a2**2*d2*dx**2-4.0d0*a1**2*a2*d2*dx**2 &
+2.2d1*a2**2*dx**2+4.4d1*a1*a2*dx**2+2.2d1*a1**2*dx**2)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (2)
              ! selection on m2: l1=4, m1=4, l2=2
              select case (m2)
                case (-2)
                  rlYlm_laplacian =7.6147536914910937d0*E*a1*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.5d1*a1**2*a2**3*dy**4-1.5d1*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2*dx**2*dy**2+9.0d1*a1**2*a2**3*dx**2*dy**2 &
+9.0d1*a1**3*a2**2*dx**2*dy**2+8.0d0*a1**2*a2**3*d2*dy**2 &
+8.0d0*a1**3*a2**2*d2*dy**2-5.2d1*a1*a2**3*dy**2 &
-1.04d2*a1**2*a2**2*dy**2-5.2d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4+8.0d0*a1**2*a2**3*d2*dx**2 &
+8.0d0*a1**3*a2**2*d2*dx**2-5.2d1*a1*a2**3*dx**2 &
-1.04d2*a1**2*a2**2*dx**2-5.2d1*a1**3*a2*dx**2-1.2d1*a1*a2**3*d2 &
-2.4d1*a1**2*a2**2*d2-1.2d1*a1**3*a2*d2+6.6d1*a2**3+1.98d2*a1*a2**2 &
+1.98d2*a1**2*a2+6.6d1*a1**3)
                case (-1)
                  rlYlm_laplacian =7.6147536914910937d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-10)*dy*(2.0d0*a1**2*a2**2*d2*dy**4 &
-1.5d1*a1*a2**2*dy**4-1.5d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+9.0d1*a1*a2**2*dx**2*dy**2 &
+9.0d1*a1**2*a2*dx**2*dy**2-4.0d0*a1*a2**2*d2*dy**2 &
-4.0d0*a1**2*a2*d2*dy**2+2.6d1*a2**2*dy**2+5.2d1*a1*a2*dy**2 &
+2.6d1*a1**2*dy**2+2.0d0*a1**2*a2**2*d2*dx**4-1.5d1*a1*a2**2*dx**4 &
-1.5d1*a1**2*a2*dx**4+1.2d1*a1*a2**2*d2*dx**2+1.2d1*a1**2*a2*d2*dx**2 &
-7.8d1*a2**2*dx**2-1.56d2*a1*a2*dx**2-7.8d1*a1**2*dx**2)*dz
                case (0)
                  rlYlm_laplacian = &
-2.19819004679753972d0*E*a1**2*a2**4*sqrt(a2+a1)*(a2+a1)**(-10)*(dy**2 &
-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)* &
(6.0d0*a1**2*a2**2*d2*dy**2-4.5d1*a1*a2**2*dy**2-4.5d1*a1**2*a2*dy**2 &
+6.0d0*a1**2*a2**2*d2*dx**2-4.5d1*a1*a2**2*dx**2-4.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+2.2d1*a1*a2**2*d2+2.2d1*a1**2*a2*d2 &
+5.2d1*a2**2+1.04d2*a1*a2+5.2d1*a1**2)
                case (1)
                  rlYlm_laplacian =7.6147536914910937d0*E*a1**2*a2**4*sqrt &
(a2+a1)*(a2+a1)**(-10)*dx*(2.0d0*a1**2*a2**2*d2*dy**4 &
-1.5d1*a1*a2**2*dy**4-1.5d1*a1**2*a2*dy**4 &
-1.2d1*a1**2*a2**2*d2*dx**2*dy**2+9.0d1*a1*a2**2*dx**2*dy**2 &
+9.0d1*a1**2*a2*dx**2*dy**2+1.2d1*a1*a2**2*d2*dy**2 &
+1.2d1*a1**2*a2*d2*dy**2-7.8d1*a2**2*dy**2-1.56d2*a1*a2*dy**2 &
-7.8d1*a1**2*dy**2+2.0d0*a1**2*a2**2*d2*dx**4-1.5d1*a1*a2**2*dx**4 &
-1.5d1*a1**2*a2*dx**4-4.0d0*a1*a2**2*d2*dx**2-4.0d0*a1**2*a2*d2*dx**2 &
+2.6d1*a2**2*dx**2+5.2d1*a1*a2*dx**2+2.6d1*a1**2*dx**2)*dz
                case (2)
                  rlYlm_laplacian = &
-3.80737684574554685d0*E*a1*a2**3*sqrt(a2+a1)*(a2+a1)**(-10)*(dy+dx)* &
(dy-1.0d0*dx)*(2.0d0*a1**3*a2**3*d2*dy**4-1.5d1*a1**2*a2**3*dy**4 &
-1.5d1*a1**3*a2**2*dy**4-1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
+9.0d1*a1**2*a2**3*dx**2*dy**2+9.0d1*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**2*a2**3*d2*dy**2-8.0d0*a1**3*a2**2*d2*dy**2 &
+5.2d1*a1*a2**3*dy**2+1.04d2*a1**2*a2**2*dy**2+5.2d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.5d1*a1**2*a2**3*dx**4 &
-1.5d1*a1**3*a2**2*dx**4-8.0d0*a1**2*a2**3*d2*dx**2 &
-8.0d0*a1**3*a2**2*d2*dx**2+5.2d1*a1*a2**3*dx**2 &
+1.04d2*a1**2*a2**2*dx**2+5.2d1*a1**3*a2*dx**2+1.2d1*a1*a2**3*d2 &
+2.4d1*a1**2*a2**2*d2+1.2d1*a1**3*a2*d2-6.6d1*a2**3-1.98d2*a1*a2**2 &
-1.98d2*a1**2*a2-6.6d1*a1**3)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (3)
              ! selection on m2: l1=4, m1=4, l2=3
              select case (m2)
                case (-3)
                  rlYlm_laplacian = &
-4.1124370130664852d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(2.0d0*a1**4*a2**4*d2*dy**6-1.7d1*a1**3*a2**4*dy**6 &
-1.7d1*a1**4*a2**3*dy**6-1.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.53d2*a1**3*a2**4*dx**2*dy**4+1.53d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**3*a2**4*d2*dy**4-1.2d1*a1**4*a2**3*d2*dy**4 &
+9.0d1*a1**2*a2**4*dy**4+1.8d2*a1**3*a2**3*dy**4 &
+9.0d1*a1**4*a2**2*dy**4+3.8d1*a1**4*a2**4*d2*dx**4*dy**2 &
-3.23d2*a1**3*a2**4*dx**4*dy**2-3.23d2*a1**4*a2**3*dx**4*dy**2 &
-2.4d1*a1**3*a2**4*d2*dx**2*dy**2-2.4d1*a1**4*a2**3*d2*dx**2*dy**2 &
+1.8d2*a1**2*a2**4*dx**2*dy**2+3.6d2*a1**3*a2**3*dx**2*dy**2 &
+1.8d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**2*a2**4*d2*dy**2 &
+7.2d1*a1**3*a2**3*d2*dy**2+3.6d1*a1**4*a2**2*d2*dy**2 &
-2.34d2*a1*a2**4*dy**2-7.02d2*a1**2*a2**3*dy**2 &
-7.02d2*a1**3*a2**2*dy**2-2.34d2*a1**4*a2*dy**2 &
-6.0d0*a1**4*a2**4*d2*dx**6+5.1d1*a1**3*a2**4*dx**6 &
+5.1d1*a1**4*a2**3*dx**6-1.2d1*a1**3*a2**4*d2*dx**4 &
-1.2d1*a1**4*a2**3*d2*dx**4+9.0d1*a1**2*a2**4*dx**4 &
+1.8d2*a1**3*a2**3*dx**4+9.0d1*a1**4*a2**2*dx**4 &
+3.6d1*a1**2*a2**4*d2*dx**2+7.2d1*a1**3*a2**3*d2*dx**2 &
+3.6d1*a1**4*a2**2*d2*dx**2-2.34d2*a1*a2**4*dx**2 &
-7.02d2*a1**2*a2**3*dx**2-7.02d2*a1**3*a2**2*dx**2 &
-2.34d2*a1**4*a2*dx**2-2.4d1*a1*a2**4*d2-7.2d1*a1**2*a2**3*d2 &
-7.2d1*a1**3*a2**2*d2-2.4d1*a1**4*a2*d2+1.32d2*a2**4+5.28d2*a1*a2**3 &
+7.92d2*a1**2*a2**2+5.28d2*a1**3*a2+1.32d2*a1**4)
                case (-2)
                  rlYlm_laplacian =2.01467445626964921d1*E*a1**2*a2**3*sqrt &
(a2+a1)*(a2+a1)**(-11)*dx*dy*(2.0d0*a1**3*a2**3*d2*dy**4 &
-1.7d1*a1**2*a2**3*dy**4-1.7d1*a1**3*a2**2*dy**4 &
-1.2d1*a1**3*a2**3*d2*dx**2*dy**2+1.02d2*a1**2*a2**3*dx**2*dy**2 &
+1.02d2*a1**3*a2**2*dx**2*dy**2+8.0d0*a1**2*a2**3*d2*dy**2 &
+8.0d0*a1**3*a2**2*d2*dy**2-6.0d1*a1*a2**3*dy**2 &
-1.2d2*a1**2*a2**2*dy**2-6.0d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.7d1*a1**2*a2**3*dx**4 &
-1.7d1*a1**3*a2**2*dx**4+8.0d0*a1**2*a2**3*d2*dx**2 &
+8.0d0*a1**3*a2**2*d2*dx**2-6.0d1*a1*a2**3*dx**2 &
-1.2d2*a1**2*a2**2*dx**2-6.0d1*a1**3*a2*dx**2-1.2d1*a1*a2**3*d2 &
-2.4d1*a1**2*a2**2*d2-1.2d1*a1**3*a2*d2+7.8d1*a2**3+2.34d2*a1*a2**2 &
+2.34d2*a1**2*a2+7.8d1*a1**3)*dz
                case (-1)
                  rlYlm_laplacian = &
-3.18548001278669409d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*dy* &
(1.0d1*a1**3*a2**3*d2*dy**6-8.5d1*a1**2*a2**3*dy**6 &
-8.5d1*a1**3*a2**2*dy**6-5.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+4.25d2*a1**2*a2**3*dx**2*dy**4+4.25d2*a1**3*a2**2*dx**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**4+4.0d1*a1**2*a2**3*d2*dy**4 &
+4.0d1*a1**3*a2**2*d2*dy**4+2.1d2*a1*a2**3*dy**4 &
+4.2d2*a1**2*a2**2*dy**4+2.1d2*a1**3*a2*dy**4 &
-5.0d1*a1**3*a2**3*d2*dx**4*dy**2+4.25d2*a1**2*a2**3*dx**4*dy**2 &
+4.25d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-3.2d2*a1**2*a2**3*d2*dx**2*dy**2-3.2d2*a1**3*a2**2*d2*dx**2*dy**2 &
-6.6d2*a1*a2**3*dx**2*dy**2-1.32d3*a1**2*a2**2*dx**2*dy**2 &
-6.6d2*a1**3*a2*dx**2*dy**2+1.6d1*a1**2*a2**3*d2**2*dy**2 &
+1.6d1*a1**3*a2**2*d2**2*dy**2-1.08d2*a1*a2**3*d2*dy**2 &
-2.16d2*a1**2*a2**2*d2*dy**2-1.08d2*a1**3*a2*d2*dy**2 &
-7.8d1*a2**3*dy**2-2.34d2*a1*a2**2*dy**2-2.34d2*a1**2*a2*dy**2 &
-7.8d1*a1**3*dy**2+1.0d1*a1**3*a2**3*d2*dx**6-8.5d1*a1**2*a2**3*dx**6 &
-8.5d1*a1**3*a2**2*dx**6-8.0d0*a1**3*a2**3*d2**2*dx**4 &
+1.2d2*a1**2*a2**3*d2*dx**4+1.2d2*a1**3*a2**2*d2*dx**4 &
-3.9d2*a1*a2**3*dx**4-7.8d2*a1**2*a2**2*dx**4-3.9d2*a1**3*a2*dx**4 &
-4.8d1*a1**2*a2**3*d2**2*dx**2-4.8d1*a1**3*a2**2*d2**2*dx**2 &
+3.24d2*a1*a2**3*d2*dx**2+6.48d2*a1**2*a2**2*d2*dx**2 &
+3.24d2*a1**3*a2*d2*dx**2+2.34d2*a2**3*dx**2+7.02d2*a1*a2**2*dx**2 &
+7.02d2*a1**2*a2*dx**2+2.34d2*a1**3*dx**2)
                case (0)
                  rlYlm_laplacian = &
-2.60093353905394473d0*E*a1**3*a2**4*sqrt(a2+a1)*(a2+a1)**(-11)*(dy**2 &
-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)* &
(1.0d1*a1**2*a2**2*d2*dy**2-8.5d1*a1*a2**2*dy**2-8.5d1*a1**2*a2*dy**2 &
+1.0d1*a1**2*a2**2*d2*dx**2-8.5d1*a1*a2**2*dx**2-8.5d1*a1**2*a2*dx**2 &
-4.0d0*a1**2*a2**2*d2**2+1.0d1*a1*a2**2*d2+1.0d1*a1**2*a2*d2 &
+1.8d2*a2**2+3.6d2*a1*a2+1.8d2*a1**2)*dz
                case (1)
                  rlYlm_laplacian = &
-3.18548001278669409d0*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(1.0d1*a1**3*a2**3*d2*dy**6-8.5d1*a1**2*a2**3*dy**6 &
-8.5d1*a1**3*a2**2*dy**6-5.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+4.25d2*a1**2*a2**3*dx**2*dy**4+4.25d2*a1**3*a2**2*dx**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**4+1.2d2*a1**2*a2**3*d2*dy**4 &
+1.2d2*a1**3*a2**2*d2*dy**4-3.9d2*a1*a2**3*dy**4 &
-7.8d2*a1**2*a2**2*dy**4-3.9d2*a1**3*a2*dy**4 &
-5.0d1*a1**3*a2**3*d2*dx**4*dy**2+4.25d2*a1**2*a2**3*dx**4*dy**2 &
+4.25d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-3.2d2*a1**2*a2**3*d2*dx**2*dy**2-3.2d2*a1**3*a2**2*d2*dx**2*dy**2 &
-6.6d2*a1*a2**3*dx**2*dy**2-1.32d3*a1**2*a2**2*dx**2*dy**2 &
-6.6d2*a1**3*a2*dx**2*dy**2-4.8d1*a1**2*a2**3*d2**2*dy**2 &
-4.8d1*a1**3*a2**2*d2**2*dy**2+3.24d2*a1*a2**3*d2*dy**2 &
+6.48d2*a1**2*a2**2*d2*dy**2+3.24d2*a1**3*a2*d2*dy**2 &
+2.34d2*a2**3*dy**2+7.02d2*a1*a2**2*dy**2+7.02d2*a1**2*a2*dy**2 &
+2.34d2*a1**3*dy**2+1.0d1*a1**3*a2**3*d2*dx**6-8.5d1*a1**2*a2**3*dx**6 &
-8.5d1*a1**3*a2**2*dx**6-8.0d0*a1**3*a2**3*d2**2*dx**4 &
+4.0d1*a1**2*a2**3*d2*dx**4+4.0d1*a1**3*a2**2*d2*dx**4 &
+2.1d2*a1*a2**3*dx**4+4.2d2*a1**2*a2**2*dx**4+2.1d2*a1**3*a2*dx**4 &
+1.6d1*a1**2*a2**3*d2**2*dx**2+1.6d1*a1**3*a2**2*d2**2*dx**2 &
-1.08d2*a1*a2**3*d2*dx**2-2.16d2*a1**2*a2**2*d2*dx**2 &
-1.08d2*a1**3*a2*d2*dx**2-7.8d1*a2**3*dx**2-2.34d2*a1*a2**2*dx**2 &
-2.34d2*a1**2*a2*dx**2-7.8d1*a1**3*dx**2)
                case (2)
                  rlYlm_laplacian = &
-1.00733722813482461d1*E*a1**2*a2**3*sqrt(a2+a1)*(a2+a1)**(-11)*(dy+dx &
)*(dy-1.0d0*dx)*(2.0d0*a1**3*a2**3*d2*dy**4-1.7d1*a1**2*a2**3*dy**4 &
-1.7d1*a1**3*a2**2*dy**4-1.2d1*a1**3*a2**3*d2*dx**2*dy**2 &
+1.02d2*a1**2*a2**3*dx**2*dy**2+1.02d2*a1**3*a2**2*dx**2*dy**2 &
-8.0d0*a1**2*a2**3*d2*dy**2-8.0d0*a1**3*a2**2*d2*dy**2 &
+6.0d1*a1*a2**3*dy**2+1.2d2*a1**2*a2**2*dy**2+6.0d1*a1**3*a2*dy**2 &
+2.0d0*a1**3*a2**3*d2*dx**4-1.7d1*a1**2*a2**3*dx**4 &
-1.7d1*a1**3*a2**2*dx**4-8.0d0*a1**2*a2**3*d2*dx**2 &
-8.0d0*a1**3*a2**2*d2*dx**2+6.0d1*a1*a2**3*dx**2 &
+1.2d2*a1**2*a2**2*dx**2+6.0d1*a1**3*a2*dx**2+1.2d1*a1*a2**3*d2 &
+2.4d1*a1**2*a2**2*d2+1.2d1*a1**3*a2*d2-7.8d1*a2**3-2.34d2*a1*a2**2 &
-2.34d2*a1**2*a2-7.8d1*a1**3)*dz
                case (3)
                  rlYlm_laplacian = &
-4.1124370130664852d0*E*a1*a2**2*sqrt(a2+a1)*(a2+a1)**(-11)*dx* &
(6.0d0*a1**4*a2**4*d2*dy**6-5.1d1*a1**3*a2**4*dy**6 &
-5.1d1*a1**4*a2**3*dy**6-3.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
+3.23d2*a1**3*a2**4*dx**2*dy**4+3.23d2*a1**4*a2**3*dx**2*dy**4 &
+1.2d1*a1**3*a2**4*d2*dy**4+1.2d1*a1**4*a2**3*d2*dy**4 &
-9.0d1*a1**2*a2**4*dy**4-1.8d2*a1**3*a2**3*dy**4 &
-9.0d1*a1**4*a2**2*dy**4+1.8d1*a1**4*a2**4*d2*dx**4*dy**2 &
-1.53d2*a1**3*a2**4*dx**4*dy**2-1.53d2*a1**4*a2**3*dx**4*dy**2 &
+2.4d1*a1**3*a2**4*d2*dx**2*dy**2+2.4d1*a1**4*a2**3*d2*dx**2*dy**2 &
-1.8d2*a1**2*a2**4*dx**2*dy**2-3.6d2*a1**3*a2**3*dx**2*dy**2 &
-1.8d2*a1**4*a2**2*dx**2*dy**2-3.6d1*a1**2*a2**4*d2*dy**2 &
-7.2d1*a1**3*a2**3*d2*dy**2-3.6d1*a1**4*a2**2*d2*dy**2 &
+2.34d2*a1*a2**4*dy**2+7.02d2*a1**2*a2**3*dy**2 &
+7.02d2*a1**3*a2**2*dy**2+2.34d2*a1**4*a2*dy**2 &
-2.0d0*a1**4*a2**4*d2*dx**6+1.7d1*a1**3*a2**4*dx**6 &
+1.7d1*a1**4*a2**3*dx**6+1.2d1*a1**3*a2**4*d2*dx**4 &
+1.2d1*a1**4*a2**3*d2*dx**4-9.0d1*a1**2*a2**4*dx**4 &
-1.8d2*a1**3*a2**3*dx**4-9.0d1*a1**4*a2**2*dx**4 &
-3.6d1*a1**2*a2**4*d2*dx**2-7.2d1*a1**3*a2**3*d2*dx**2 &
-3.6d1*a1**4*a2**2*d2*dx**2+2.34d2*a1*a2**4*dx**2 &
+7.02d2*a1**2*a2**3*dx**2+7.02d2*a1**3*a2**2*dx**2 &
+2.34d2*a1**4*a2*dx**2+2.4d1*a1*a2**4*d2+7.2d1*a1**2*a2**3*d2 &
+7.2d1*a1**3*a2**2*d2+2.4d1*a1**4*a2*d2-1.32d2*a2**4-5.28d2*a1*a2**3 &
-7.92d2*a1**2*a2**2-5.28d2*a1**3*a2-1.32d2*a1**4)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case (4)
              ! selection on m2: l1=4, m1=4, l2=4
              select case (m2)
                case (-4)
                  rlYlm_laplacian = &
-1.74475925948511734d1*E*a1**5*a2**5*sqrt(a2+a1)*(a2+a1)**(-12)* &
(1.0d0*a2*(2.0d0*a1*d2-1.9d1)-1.9d1*a1)*dx*dy*(dy+dx)*(dy-1.0d0*dx)* &
(dy**2-2.0d0*dx*dy-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)
                case (-3)
                  rlYlm_laplacian = &
-1.23373110391994556d1*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(2.0d0*a1**4*a2**4*d2*dy**6-1.9d1*a1**3*a2**4*dy**6 &
-1.9d1*a1**4*a2**3*dy**6-1.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
+1.71d2*a1**3*a2**4*dx**2*dy**4+1.71d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**3*a2**4*d2*dy**4-1.2d1*a1**4*a2**3*d2*dy**4 &
+1.02d2*a1**2*a2**4*dy**4+2.04d2*a1**3*a2**3*dy**4 &
+1.02d2*a1**4*a2**2*dy**4+3.8d1*a1**4*a2**4*d2*dx**4*dy**2 &
-3.61d2*a1**3*a2**4*dx**4*dy**2-3.61d2*a1**4*a2**3*dx**4*dy**2 &
-2.4d1*a1**3*a2**4*d2*dx**2*dy**2-2.4d1*a1**4*a2**3*d2*dx**2*dy**2 &
+2.04d2*a1**2*a2**4*dx**2*dy**2+4.08d2*a1**3*a2**3*dx**2*dy**2 &
+2.04d2*a1**4*a2**2*dx**2*dy**2+3.6d1*a1**2*a2**4*d2*dy**2 &
+7.2d1*a1**3*a2**3*d2*dy**2+3.6d1*a1**4*a2**2*d2*dy**2 &
-2.7d2*a1*a2**4*dy**2-8.1d2*a1**2*a2**3*dy**2-8.1d2*a1**3*a2**2*dy**2 &
-2.7d2*a1**4*a2*dy**2-6.0d0*a1**4*a2**4*d2*dx**6 &
+5.7d1*a1**3*a2**4*dx**6+5.7d1*a1**4*a2**3*dx**6 &
-1.2d1*a1**3*a2**4*d2*dx**4-1.2d1*a1**4*a2**3*d2*dx**4 &
+1.02d2*a1**2*a2**4*dx**4+2.04d2*a1**3*a2**3*dx**4 &
+1.02d2*a1**4*a2**2*dx**4+3.6d1*a1**2*a2**4*d2*dx**2 &
+7.2d1*a1**3*a2**3*d2*dx**2+3.6d1*a1**4*a2**2*d2*dx**2 &
-2.7d2*a1*a2**4*dx**2-8.1d2*a1**2*a2**3*dx**2-8.1d2*a1**3*a2**2*dx**2 &
-2.7d2*a1**4*a2*dx**2-2.4d1*a1*a2**4*d2-7.2d1*a1**2*a2**3*d2 &
-7.2d1*a1**3*a2**2*d2-2.4d1*a1**4*a2*d2+1.56d2*a2**4+6.24d2*a1*a2**3 &
+9.36d2*a1**2*a2**2+6.24d2*a1**3*a2+1.56d2*a1**4)*dz
                case (-2)
                  rlYlm_laplacian = &
-6.59457014039261917d0*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx*dy* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.33d2*a1**3*a2**4*dy**6 &
-1.33d2*a1**4*a2**3*dy**6-7.0d1*a1**4*a2**4*d2*dx**2*dy**4 &
+6.65d2*a1**3*a2**4*dx**2*dy**4+6.65d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+1.62d2*a1**3*a2**4*d2*dy**4 &
+1.62d2*a1**4*a2**3*d2*dy**4-4.08d2*a1**2*a2**4*dy**4 &
-8.16d2*a1**3*a2**3*dy**4-4.08d2*a1**4*a2**2*dy**4 &
-7.0d1*a1**4*a2**4*d2*dx**4*dy**2+6.65d2*a1**3*a2**4*dx**4*dy**2 &
+6.65d2*a1**4*a2**3*dx**4*dy**2+7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-5.24d2*a1**3*a2**4*d2*dx**2*dy**2-5.24d2*a1**4*a2**3*d2*dx**2*dy**2 &
-1.36d3*a1**2*a2**4*dx**2*dy**2-2.72d3*a1**3*a2**3*dx**2*dy**2 &
-1.36d3*a1**4*a2**2*dx**2*dy**2-4.8d1*a1**3*a2**4*d2**2*dy**2 &
-4.8d1*a1**4*a2**3*d2**2*dy**2+3.0d2*a1**2*a2**4*d2*dy**2 &
+6.0d2*a1**3*a2**3*d2*dy**2+3.0d2*a1**4*a2**2*d2*dy**2 &
+8.1d2*a1*a2**4*dy**2+2.43d3*a1**2*a2**3*dy**2 &
+2.43d3*a1**3*a2**2*dy**2+8.1d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.33d2*a1**3*a2**4*dx**6 &
-1.33d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+1.62d2*a1**3*a2**4*d2*dx**4+1.62d2*a1**4*a2**3*d2*dx**4 &
-4.08d2*a1**2*a2**4*dx**4-8.16d2*a1**3*a2**3*dx**4 &
-4.08d2*a1**4*a2**2*dx**4-4.8d1*a1**3*a2**4*d2**2*dx**2 &
-4.8d1*a1**4*a2**3*d2**2*dx**2+3.0d2*a1**2*a2**4*d2*dx**2 &
+6.0d2*a1**3*a2**3*d2*dx**2+3.0d2*a1**4*a2**2*d2*dx**2 &
+8.1d2*a1*a2**4*dx**2+2.43d3*a1**2*a2**3*dx**2 &
+2.43d3*a1**3*a2**2*dx**2+8.1d2*a1**4*a2*dx**2+7.2d1*a1**2*a2**4*d2**2 &
+1.44d2*a1**3*a2**3*d2**2+7.2d1*a1**4*a2**2*d2**2-5.16d2*a1*a2**4*d2 &
-1.548d3*a1**2*a2**3*d2-1.548d3*a1**3*a2**2*d2-5.16d2*a1**4*a2*d2 &
-1.56d2*a2**4-6.24d2*a1*a2**3-9.36d2*a1**2*a2**2-6.24d2*a1**3*a2 &
-1.56d2*a1**4)
                case (-1)
                  rlYlm_laplacian = &
-4.66306526528194375d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dy* &
(1.4d1*a1**3*a2**3*d2*dy**6-1.33d2*a1**2*a2**3*dy**6 &
-1.33d2*a1**3*a2**2*dy**6-7.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+6.65d2*a1**2*a2**3*dx**2*dy**4+6.65d2*a1**3*a2**2*dx**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**4+2.4d1*a1**2*a2**3*d2*dy**4 &
+2.4d1*a1**3*a2**2*d2*dy**4+4.42d2*a1*a2**3*dy**4 &
+8.84d2*a1**2*a2**2*dy**4+4.42d2*a1**3*a2*dy**4 &
-7.0d1*a1**3*a2**3*d2*dx**4*dy**2+6.65d2*a1**2*a2**3*dx**4*dy**2 &
+6.65d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-2.56d2*a1**2*a2**3*d2*dx**2*dy**2-2.56d2*a1**3*a2**2*d2*dx**2*dy**2 &
-1.7d3*a1*a2**3*dx**2*dy**2-3.4d3*a1**2*a2**2*dx**2*dy**2 &
-1.7d3*a1**3*a2*dx**2*dy**2+1.6d1*a1**2*a2**3*d2**2*dy**2 &
+1.6d1*a1**3*a2**2*d2**2*dy**2-1.0d2*a1*a2**3*d2*dy**2 &
-2.0d2*a1**2*a2**2*d2*dy**2-1.0d2*a1**3*a2*d2*dy**2-2.7d2*a2**3*dy**2 &
-8.1d2*a1*a2**2*dy**2-8.1d2*a1**2*a2*dy**2-2.7d2*a1**3*dy**2 &
+1.4d1*a1**3*a2**3*d2*dx**6-1.33d2*a1**2*a2**3*dx**6 &
-1.33d2*a1**3*a2**2*dx**6-8.0d0*a1**3*a2**3*d2**2*dx**4 &
+1.36d2*a1**2*a2**3*d2*dx**4+1.36d2*a1**3*a2**2*d2*dx**4 &
-5.1d2*a1*a2**3*dx**4-1.02d3*a1**2*a2**2*dx**4-5.1d2*a1**3*a2*dx**4 &
-4.8d1*a1**2*a2**3*d2**2*dx**2-4.8d1*a1**3*a2**2*d2**2*dx**2 &
+3.0d2*a1*a2**3*d2*dx**2+6.0d2*a1**2*a2**2*d2*dx**2 &
+3.0d2*a1**3*a2*d2*dx**2+8.1d2*a2**3*dx**2+2.43d3*a1*a2**2*dx**2 &
+2.43d3*a1**2*a2*dx**2+8.1d2*a1**3*dx**2)*dz
                case (0)
                  rlYlm_laplacian =7.37295355815411407d-1*E*a1**3*a2**3*sqrt &
                  (a2+a1)*(a2+a1)**(-12)*(dy**2-2.0d0*dx*dy &
-1.0d0*dx**2)*(dy**2+2.0d0*dx*dy-1.0d0*dx**2)* &
(7.0d1*a1**3*a2**3*d2*dy**4-6.65d2*a1**2*a2**3*dy**4 &
-6.65d2*a1**3*a2**2*dy**4+1.4d2*a1**3*a2**3*d2*dx**2*dy**2 &
-1.33d3*a1**2*a2**3*dx**2*dy**2-1.33d3*a1**3*a2**2*dx**2*dy**2 &
-8.0d1*a1**3*a2**3*d2**2*dy**2+5.2d2*a1**2*a2**3*d2*dy**2 &
+5.2d2*a1**3*a2**2*d2*dy**2+2.04d3*a1*a2**3*dy**2 &
+4.08d3*a1**2*a2**2*dy**2+2.04d3*a1**3*a2*dy**2 &
+7.0d1*a1**3*a2**3*d2*dx**4-6.65d2*a1**2*a2**3*dx**4 &
-6.65d2*a1**3*a2**2*dx**4-8.0d1*a1**3*a2**3*d2**2*dx**2 &
+5.2d2*a1**2*a2**3*d2*dx**2+5.2d2*a1**3*a2**2*d2*dx**2 &
+2.04d3*a1*a2**3*dx**2+4.08d3*a1**2*a2**2*dx**2+2.04d3*a1**3*a2*dx**2 &
+1.6d1*a1**3*a2**3*d2**3+4.0d1*a1**2*a2**3*d2**2 &
+4.0d1*a1**3*a2**2*d2**2-1.56d3*a1*a2**3*d2-3.12d3*a1**2*a2**2*d2 &
-1.56d3*a1**3*a2*d2-5.4d2*a2**3-1.62d3*a1*a2**2-1.62d3*a1**2*a2 &
-5.4d2*a1**3)
                case (1)
                  rlYlm_laplacian = &
-4.66306526528194375d0*E*a1**3*a2**3*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(1.4d1*a1**3*a2**3*d2*dy**6-1.33d2*a1**2*a2**3*dy**6 &
-1.33d2*a1**3*a2**2*dy**6-7.0d1*a1**3*a2**3*d2*dx**2*dy**4 &
+6.65d2*a1**2*a2**3*dx**2*dy**4+6.65d2*a1**3*a2**2*dx**2*dy**4 &
-8.0d0*a1**3*a2**3*d2**2*dy**4+1.36d2*a1**2*a2**3*d2*dy**4 &
+1.36d2*a1**3*a2**2*d2*dy**4-5.1d2*a1*a2**3*dy**4 &
-1.02d3*a1**2*a2**2*dy**4-5.1d2*a1**3*a2*dy**4 &
-7.0d1*a1**3*a2**3*d2*dx**4*dy**2+6.65d2*a1**2*a2**3*dx**4*dy**2 &
+6.65d2*a1**3*a2**2*dx**4*dy**2+4.8d1*a1**3*a2**3*d2**2*dx**2*dy**2 &
-2.56d2*a1**2*a2**3*d2*dx**2*dy**2-2.56d2*a1**3*a2**2*d2*dx**2*dy**2 &
-1.7d3*a1*a2**3*dx**2*dy**2-3.4d3*a1**2*a2**2*dx**2*dy**2 &
-1.7d3*a1**3*a2*dx**2*dy**2-4.8d1*a1**2*a2**3*d2**2*dy**2 &
-4.8d1*a1**3*a2**2*d2**2*dy**2+3.0d2*a1*a2**3*d2*dy**2 &
+6.0d2*a1**2*a2**2*d2*dy**2+3.0d2*a1**3*a2*d2*dy**2+8.1d2*a2**3*dy**2 &
+2.43d3*a1*a2**2*dy**2+2.43d3*a1**2*a2*dy**2+8.1d2*a1**3*dy**2 &
+1.4d1*a1**3*a2**3*d2*dx**6-1.33d2*a1**2*a2**3*dx**6 &
-1.33d2*a1**3*a2**2*dx**6-8.0d0*a1**3*a2**3*d2**2*dx**4 &
+2.4d1*a1**2*a2**3*d2*dx**4+2.4d1*a1**3*a2**2*d2*dx**4 &
+4.42d2*a1*a2**3*dx**4+8.84d2*a1**2*a2**2*dx**4+4.42d2*a1**3*a2*dx**4 &
+1.6d1*a1**2*a2**3*d2**2*dx**2+1.6d1*a1**3*a2**2*d2**2*dx**2 &
-1.0d2*a1*a2**3*d2*dx**2-2.0d2*a1**2*a2**2*d2*dx**2 &
-1.0d2*a1**3*a2*d2*dx**2-2.7d2*a2**3*dx**2-8.1d2*a1*a2**2*dx**2 &
-8.1d2*a1**2*a2*dx**2-2.7d2*a1**3*dx**2)*dz
                case (2)
                  rlYlm_laplacian =3.29728507019630958d0*E*a1**2*a2**2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(dy+dx)*(dy-1.0d0*dx)* &
(1.4d1*a1**4*a2**4*d2*dy**6-1.33d2*a1**3*a2**4*dy**6 &
-1.33d2*a1**4*a2**3*dy**6-7.0d1*a1**4*a2**4*d2*dx**2*dy**4 &
+6.65d2*a1**3*a2**4*dx**2*dy**4+6.65d2*a1**4*a2**3*dx**2*dy**4 &
-1.2d1*a1**4*a2**4*d2**2*dy**4+5.0d1*a1**3*a2**4*d2*dy**4 &
+5.0d1*a1**4*a2**3*d2*dy**4+5.44d2*a1**2*a2**4*dy**4 &
+1.088d3*a1**3*a2**3*dy**4+5.44d2*a1**4*a2**2*dy**4 &
-7.0d1*a1**4*a2**4*d2*dx**4*dy**2+6.65d2*a1**3*a2**4*dx**4*dy**2 &
+6.65d2*a1**4*a2**3*dx**4*dy**2+7.2d1*a1**4*a2**4*d2**2*dx**2*dy**2 &
-7.48d2*a1**3*a2**4*d2*dx**2*dy**2-7.48d2*a1**4*a2**3*d2*dx**2*dy**2 &
+5.44d2*a1**2*a2**4*dx**2*dy**2+1.088d3*a1**3*a2**3*dx**2*dy**2 &
+5.44d2*a1**4*a2**2*dx**2*dy**2+4.8d1*a1**3*a2**4*d2**2*dy**2 &
+4.8d1*a1**4*a2**3*d2**2*dy**2-3.0d2*a1**2*a2**4*d2*dy**2 &
-6.0d2*a1**3*a2**3*d2*dy**2-3.0d2*a1**4*a2**2*d2*dy**2 &
-8.1d2*a1*a2**4*dy**2-2.43d3*a1**2*a2**3*dy**2 &
-2.43d3*a1**3*a2**2*dy**2-8.1d2*a1**4*a2*dy**2 &
+1.4d1*a1**4*a2**4*d2*dx**6-1.33d2*a1**3*a2**4*dx**6 &
-1.33d2*a1**4*a2**3*dx**6-1.2d1*a1**4*a2**4*d2**2*dx**4 &
+5.0d1*a1**3*a2**4*d2*dx**4+5.0d1*a1**4*a2**3*d2*dx**4 &
+5.44d2*a1**2*a2**4*dx**4+1.088d3*a1**3*a2**3*dx**4 &
+5.44d2*a1**4*a2**2*dx**4+4.8d1*a1**3*a2**4*d2**2*dx**2 &
+4.8d1*a1**4*a2**3*d2**2*dx**2-3.0d2*a1**2*a2**4*d2*dx**2 &
-6.0d2*a1**3*a2**3*d2*dx**2-3.0d2*a1**4*a2**2*d2*dx**2 &
-8.1d2*a1*a2**4*dx**2-2.43d3*a1**2*a2**3*dx**2 &
-2.43d3*a1**3*a2**2*dx**2-8.1d2*a1**4*a2*dx**2-7.2d1*a1**2*a2**4*d2**2 &
-1.44d2*a1**3*a2**3*d2**2-7.2d1*a1**4*a2**2*d2**2+5.16d2*a1*a2**4*d2 &
+1.548d3*a1**2*a2**3*d2+1.548d3*a1**3*a2**2*d2+5.16d2*a1**4*a2*d2 &
+1.56d2*a2**4+6.24d2*a1*a2**3+9.36d2*a1**2*a2**2+6.24d2*a1**3*a2 &
+1.56d2*a1**4)
                case (3)
                  rlYlm_laplacian = &
-1.23373110391994556d1*E*a1**2*a2**2*sqrt(a2+a1)*(a2+a1)**(-12)*dx* &
(6.0d0*a1**4*a2**4*d2*dy**6-5.7d1*a1**3*a2**4*dy**6 &
-5.7d1*a1**4*a2**3*dy**6-3.8d1*a1**4*a2**4*d2*dx**2*dy**4 &
+3.61d2*a1**3*a2**4*dx**2*dy**4+3.61d2*a1**4*a2**3*dx**2*dy**4 &
+1.2d1*a1**3*a2**4*d2*dy**4+1.2d1*a1**4*a2**3*d2*dy**4 &
-1.02d2*a1**2*a2**4*dy**4-2.04d2*a1**3*a2**3*dy**4 &
-1.02d2*a1**4*a2**2*dy**4+1.8d1*a1**4*a2**4*d2*dx**4*dy**2 &
-1.71d2*a1**3*a2**4*dx**4*dy**2-1.71d2*a1**4*a2**3*dx**4*dy**2 &
+2.4d1*a1**3*a2**4*d2*dx**2*dy**2+2.4d1*a1**4*a2**3*d2*dx**2*dy**2 &
-2.04d2*a1**2*a2**4*dx**2*dy**2-4.08d2*a1**3*a2**3*dx**2*dy**2 &
-2.04d2*a1**4*a2**2*dx**2*dy**2-3.6d1*a1**2*a2**4*d2*dy**2 &
-7.2d1*a1**3*a2**3*d2*dy**2-3.6d1*a1**4*a2**2*d2*dy**2 &
+2.7d2*a1*a2**4*dy**2+8.1d2*a1**2*a2**3*dy**2+8.1d2*a1**3*a2**2*dy**2 &
+2.7d2*a1**4*a2*dy**2-2.0d0*a1**4*a2**4*d2*dx**6 &
+1.9d1*a1**3*a2**4*dx**6+1.9d1*a1**4*a2**3*dx**6 &
+1.2d1*a1**3*a2**4*d2*dx**4+1.2d1*a1**4*a2**3*d2*dx**4 &
-1.02d2*a1**2*a2**4*dx**4-2.04d2*a1**3*a2**3*dx**4 &
-1.02d2*a1**4*a2**2*dx**4-3.6d1*a1**2*a2**4*d2*dx**2 &
-7.2d1*a1**3*a2**3*d2*dx**2-3.6d1*a1**4*a2**2*d2*dx**2 &
+2.7d2*a1*a2**4*dx**2+8.1d2*a1**2*a2**3*dx**2+8.1d2*a1**3*a2**2*dx**2 &
+2.7d2*a1**4*a2*dx**2+2.4d1*a1*a2**4*d2+7.2d1*a1**2*a2**3*d2 &
+7.2d1*a1**3*a2**2*d2+2.4d1*a1**4*a2*d2-1.56d2*a2**4-6.24d2*a1*a2**3 &
-9.36d2*a1**2*a2**2-6.24d2*a1**3*a2-1.56d2*a1**4)*dz
                case (4)
                  rlYlm_laplacian =4.36189814871279335d0*E*a1*a2*sqrt &
(a2+a1)*(a2+a1)**(-12)*(2.0d0*a1**5*a2**5*d2*dy**8 &
-1.9d1*a1**4*a2**5*dy**8-1.9d1*a1**5*a2**4*dy**8 &
-2.4d1*a1**5*a2**5*d2*dx**2*dy**6+2.28d2*a1**4*a2**5*dx**2*dy**6 &
+2.28d2*a1**5*a2**4*dx**2*dy**6-1.6d1*a1**4*a2**5*d2*dy**6 &
-1.6d1*a1**5*a2**4*d2*dy**6+1.36d2*a1**3*a2**5*dy**6 &
+2.72d2*a1**4*a2**4*dy**6+1.36d2*a1**5*a2**3*dy**6 &
+7.6d1*a1**5*a2**5*d2*dx**4*dy**4-7.22d2*a1**4*a2**5*dx**4*dy**4 &
-7.22d2*a1**5*a2**4*dx**4*dy**4-4.8d1*a1**4*a2**5*d2*dx**2*dy**4 &
-4.8d1*a1**5*a2**4*d2*dx**2*dy**4+4.08d2*a1**3*a2**5*dx**2*dy**4 &
+8.16d2*a1**4*a2**4*dx**2*dy**4+4.08d2*a1**5*a2**3*dx**2*dy**4 &
+7.2d1*a1**3*a2**5*d2*dy**4+1.44d2*a1**4*a2**4*d2*dy**4 &
+7.2d1*a1**5*a2**3*d2*dy**4-5.4d2*a1**2*a2**5*dy**4 &
-1.62d3*a1**3*a2**4*dy**4-1.62d3*a1**4*a2**3*dy**4 &
-5.4d2*a1**5*a2**2*dy**4-2.4d1*a1**5*a2**5*d2*dx**6*dy**2 &
+2.28d2*a1**4*a2**5*dx**6*dy**2+2.28d2*a1**5*a2**4*dx**6*dy**2 &
-4.8d1*a1**4*a2**5*d2*dx**4*dy**2-4.8d1*a1**5*a2**4*d2*dx**4*dy**2 &
+4.08d2*a1**3*a2**5*dx**4*dy**2+8.16d2*a1**4*a2**4*dx**4*dy**2 &
+4.08d2*a1**5*a2**3*dx**4*dy**2+1.44d2*a1**3*a2**5*d2*dx**2*dy**2 &
+2.88d2*a1**4*a2**4*d2*dx**2*dy**2+1.44d2*a1**5*a2**3*d2*dx**2*dy**2 &
-1.08d3*a1**2*a2**5*dx**2*dy**2-3.24d3*a1**3*a2**4*dx**2*dy**2 &
-3.24d3*a1**4*a2**3*dx**2*dy**2-1.08d3*a1**5*a2**2*dx**2*dy**2 &
-9.6d1*a1**2*a2**5*d2*dy**2-2.88d2*a1**3*a2**4*d2*dy**2 &
-2.88d2*a1**4*a2**3*d2*dy**2-9.6d1*a1**5*a2**2*d2*dy**2 &
+6.24d2*a1*a2**5*dy**2+2.496d3*a1**2*a2**4*dy**2 &
+3.744d3*a1**3*a2**3*dy**2+2.496d3*a1**4*a2**2*dy**2 &
+6.24d2*a1**5*a2*dy**2+2.0d0*a1**5*a2**5*d2*dx**8 &
-1.9d1*a1**4*a2**5*dx**8-1.9d1*a1**5*a2**4*dx**8 &
-1.6d1*a1**4*a2**5*d2*dx**6-1.6d1*a1**5*a2**4*d2*dx**6 &
+1.36d2*a1**3*a2**5*dx**6+2.72d2*a1**4*a2**4*dx**6 &
+1.36d2*a1**5*a2**3*dx**6+7.2d1*a1**3*a2**5*d2*dx**4 &
+1.44d2*a1**4*a2**4*d2*dx**4+7.2d1*a1**5*a2**3*d2*dx**4 &
-5.4d2*a1**2*a2**5*dx**4-1.62d3*a1**3*a2**4*dx**4 &
-1.62d3*a1**4*a2**3*dx**4-5.4d2*a1**5*a2**2*dx**4 &
-9.6d1*a1**2*a2**5*d2*dx**2-2.88d2*a1**3*a2**4*d2*dx**2 &
-2.88d2*a1**4*a2**3*d2*dx**2-9.6d1*a1**5*a2**2*d2*dx**2 &
+6.24d2*a1*a2**5*dx**2+2.496d3*a1**2*a2**4*dx**2 &
+3.744d3*a1**3*a2**3*dx**2+2.496d3*a1**4*a2**2*dx**2 &
+6.24d2*a1**5*a2*dx**2+2.4d1*a1*a2**5*d2+9.6d1*a1**2*a2**4*d2 &
+1.44d2*a1**3*a2**3*d2+9.6d1*a1**4*a2**2*d2+2.4d1*a1**5*a2*d2 &
-1.32d2*a2**5-6.6d2*a1*a2**4-1.32d3*a1**2*a2**3-1.32d3*a1**3*a2**2 &
-6.6d2*a1**4*a2-1.32d2*a1**5)
                case default
                  print * &
,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=',m1,'l2=',l2 &
,'m2=',m2
                  stop
              end select
            case default
              print *,'Error: rlYlm_overlap not implemented for l1=' &
,l1,'m1=',m1,'l2=',l2,'m2=',m2
              stop
          end select
        case default
          print *,'Error: rlYlm_overlap not implemented for l1=',l1 &
,'m1=',m1,'l2=',l2,'m2=',m2
          stop
      end select
    case default
      print *,'Error: rlYlm_overlap not implemented for l1=',l1,'m1=' &
,m1,'l2=',l2,'m2=',m2
      stop
  end select
  
end function

end module mod_rlYlm_laplacian
