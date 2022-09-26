!> @file Electrostatics.f90
!!
!! Defines utility routines for computing coulomb potential and field
!! issued from gaussian basis elements.
!!
!! Author: I. Duchemin July 2015
!!

!> compute the potential generated at position \f$r_{test}\f$ by a Cubic Harmonic Orbital centered in r.
!!
!!  \f$
!!      \int dr' \, \frac{1}{|r_{test}-r'|} Y_{xyz}(r'-r)
!!  \f$
!! 
recursive function potential_from_C(r_test,alpha,r,nx,ny,nz) result(value)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: nx         !< x power for first cubic Harmonic
  integer      :: ny         !< y power for first cubic Harmonic
  integer      :: nz         !< z power for first cubic Harmonic
  
  ! return value
  real(kind=8) :: value
  
  ! parameters
  integer , parameter     :: l_test       =0
  integer , parameter     :: m_test       =0
  real(kind=8), parameter :: alpha_test   =1.0d16
  real(kind=8), parameter :: charge_norm  =2.0d0/(4.0d0*atan(1.0d0))*(alpha_test*sqrt(alpha_test))
  integer, parameter      :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! local variable
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the S solid harmonic into cubic harmonic
  call R_from_Y(0,0,c1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx,ny,nz)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx+ny+nz))=0.0d0
  c2(ir2)=1.0d0
  
  ! set value
  value = charge_norm * R_X_R(alpha_test,r_test,c1,0,alpha,r,c2,nx+ny+nz)
  
end function

!> compute the potential generated at position r_test by a Solid Harmonic Orbital centered in r.
!!
!!  \f$
!!      \int dr' \, \frac{1}{|r_{test}-r'|} Y_{lm}(r'-r)
!!  \f$
!!
recursive function potential_from_Y(r_test,alpha,r,l,m) result(value)
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: l          !< l momentum of source the charge distribution
  integer      :: m          !< m momentum of source the charge distribution
  
  ! return value
  real(kind=8) :: value
  
  ! parameters
  integer , parameter     :: l_test       =0
  integer , parameter     :: m_test       =0
  real(kind=8), parameter :: alpha_test   =1.0d16
  real(kind=8), parameter :: charge_norm  =2.0d0/(4.0d0*atan(1.0d0))*(alpha_test*sqrt(alpha_test))
  
  ! coulomb interaction function
  real(kind=8) :: Y_Coulomb_Y
  
  ! set value
  value = charge_norm * Y_Coulomb_Y(alpha_test, r_test, l_test, m_test, alpha, r, l, m)
  
end function
  
!> compute the three component of the electric field generated at position r_test by a Solid Harmonic Orbital centered in r
!!
!!  \f$
!!     \nabla \int dr' \, \frac{1}{|r_{test}-r'|} Y_{lm}(r'-r)
!!  \f$
!!
recursive subroutine field_from_Y(r_test,alpha,r,l,m,E_test)
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: l          !< l momentum of source the charge distribution
  integer      :: m          !< m momentum of source the charge distribution
  
  ! return value
  real(kind=8) :: E_test(3)  !< electric field at test position
  
  ! parameters
  real(kind=8), parameter :: alpha_test   =1.0d16
  real(kind=8), parameter :: charge_norm  =-4.0d0/(4.0d0*atan(1.0d0))/sqrt(3.0d0)*(alpha_test**2*sqrt(alpha_test))
  
  ! coulomb interaction function
  real(kind=8) :: Y_Coulomb_Y
  
  ! set value
  E_test(1) = charge_norm * Y_Coulomb_Y(alpha_test, r_test, 1, 1,alpha,r,l,m)
  E_test(2) = charge_norm * Y_Coulomb_Y(alpha_test, r_test, 1,-1,alpha,r,l,m)
  E_test(3) = charge_norm * Y_Coulomb_Y(alpha_test, r_test, 1, 0,alpha,r,l,m)
  
end subroutine  
  
!> compute the three component of the electric field generated at position r_test by a Cubic Harmonic Orbital centered in r
!!
!!  \f$
!!     \nabla \int dr' \, \frac{1}{|r_{test}-r'|} Y_{xyz}(r'-r)
!!  \f$
!!
recursive subroutine field_from_C(r_test,alpha,r,nx,ny,nz,E_test)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: nx         !< x power for first cubic Harmonic
  integer      :: ny         !< y power for first cubic Harmonic
  integer      :: nz         !< z power for first cubic Harmonic
  
  ! return value
  real(kind=8) :: E_test(3)  !< electric field at test position
  
  ! parameters
  real(kind=8), parameter :: alpha_test   =1.0d16
  real(kind=8), parameter :: charge_norm  =-4.0d0/(4.0d0*atan(1.0d0))/sqrt(3.0d0)*(alpha_test**2*sqrt(alpha_test))
  integer, parameter      :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! local variable
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx,ny,nz)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx+ny+nz))=0.0d0
  c2(ir2)=1.0d0
  
  ! decomposition of the Px solid harmonic into cubic harmonic
  call R_from_Y(1,1,c1)
  
  ! set value
  E_test(1) = charge_norm * R_X_R(alpha_test,r_test,c1,1,alpha,r,c2,nx+ny+nz)
  
  ! decomposition of the Py solid harmonic into cubic harmonic
  call R_from_Y(1,-1,c1)
  
  ! set value
  E_test(2) = charge_norm * R_X_R(alpha_test,r_test,c1,1,alpha,r,c2,nx+ny+nz)
  
  ! decomposition of the Pz solid harmonic into cubic harmonic
  call R_from_Y(1, 0,c1)
  
  ! set value
  E_test(3) = charge_norm * R_X_R(alpha_test,r_test,c1,1,alpha,r,c2,nx+ny+nz)
  
end subroutine




!!> returns field and potential coming form a cubic harmonic charge distribution 
!!!
recursive subroutine electrostatics_from_C(r_test,alpha,r,nx,ny,nz,V_test,E_test)
  
  use mod_D_X_V_and_E
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: nx         !< x power for first cubic Harmonic
  integer      :: ny         !< y power for first cubic Harmonic
  integer      :: nz         !< z power for first cubic Harmonic
  
  ! return value
  real(kind=8) :: V_test     !< electric potential at test position
  real(kind=8) :: E_test(3)  !< electric field at test position
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer, parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! temp variables
  real(kind=8) :: norm 
  real(kind=8) :: V(455) 
  real(kind=8) :: Ex(455) 
  real(kind=8) :: Ey(455) 
  real(kind=8) :: Ez(455) 
  real(kind=8) :: coeffs(455) ! holds the D decomp of source orbital 
  
  ! compute contribution from this index
  coeffs=0.0d0
  call R_to_D(alpha,ir_index(nx,ny,nz),coeffs)
  
  ! call elctrostatic routine
  if ( .not.D_X_V_and_E(alpha, nx+ny+nz, r_test-r, V, Ex, Ey, Ez) ) then
    print *,"problem"
    stop  1
  end if
  
  ! get state norm 
  norm = sqrt(PI/alpha)*(PI/alpha)
  
  ! compute result
  V_test    = sum( coeffs(1:imax(nx+ny+nz)) * V(1:imax(nx+ny+nz))  ) * norm
  E_test(1) =-sum( coeffs(1:imax(nx+ny+nz)) * Ex(1:imax(nx+ny+nz)) ) * norm
  E_test(2) =-sum( coeffs(1:imax(nx+ny+nz)) * Ey(1:imax(nx+ny+nz)) ) * norm
  E_test(3) =-sum( coeffs(1:imax(nx+ny+nz)) * Ez(1:imax(nx+ny+nz)) ) * norm
  
end subroutine
  
  
!!> returns field and potential coming form a solid harmonic charge distribution 
!!!
recursive subroutine electrostatics_from_Y(r_test,alpha,r,l,m,V_test,E_test)
  
  use mod_D_X_V_and_E
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: l          !< l momentum of solid Harmonic
  integer      :: m          !< m momentum of solid Harmonic
  
  ! return value
  real(kind=8) :: V_test     !< electric potential at test position
  real(kind=8) :: E_test(3)  !< electric field at test position
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer, parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! temp variables
  integer      :: i
  real(kind=8) :: norm 
  real(kind=8) :: V(455) 
  real(kind=8) :: Ex(455) 
  real(kind=8) :: Ey(455) 
  real(kind=8) :: Ez(455) 
  real(kind=8) :: rD(3)       ! center for the derivatives basis
  real(kind=8) :: alphaD      ! exponent for the derivatives basis
  integer      :: iD(455)     ! corresponding derivative index array
  real(kind=8) :: cD(455)     ! corresponding derivative coefficient array
  integer      :: nD          ! number of corresponding derivatives
  
  ! compute contribution from this index
  call Y_to_D(alpha,r,l,m,alphaD,rD,cD,iD,nD)
  
  ! call elctrostatic routine
  if ( .not.D_X_V_and_E(alphaD, l, r_test-rD, V, Ex, Ey, Ez) ) then
    print *,"problem"
    stop  1
  end if
  
  ! get state norm 
  norm = sqrt(PI/alphaD)*(PI/alphaD)
  
  ! compute result
  V_test=0.0d0
  E_test=0.0d0
  do i=1,nD
    V_test    = V_test    + cD(i) * V( iD(i) )
    E_test(1) = E_test(1) + cD(i) * Ex( iD(i) )
    E_test(2) = E_test(2) + cD(i) * Ey( iD(i) )
    E_test(3) = E_test(3) + cD(i) * Ez( iD(i) )
  end do
  V_test = V_test * norm
  E_test =-E_test * norm
  
end subroutine


!!> returns total coulomb energy form interaction of hermit harmonics
!!! charge distribution and weighted source charges and dipoles 
!!!
recursive subroutine cumulative_electrostatics_on_D(r_source,  & !< source charge and dipole location
                                                   q_source,  & !< source charge value
                                                   dx_source, & !< source x dipole value
                                                   dy_source, & !< source y dipole value
                                                   dz_source, & !< source z dipole value
                                                   n_source,  & !< number of source charges and dipoles
                                                   alpha,     & !< exponent of cubic harmonic test charge distribution
                                                   r,         & !< center of cubic harmonic test charge distribution
                                                   l,         & !< max l for hermit harmonic test charge distribution
                                                   V_test,    & !< cubic harmonic test charge distribution cumulated interaction with source charge
                                                   Ex_test,   & !< cubic harmonic test charge distribution cumulated interaction with source x dipoles
                                                   Ey_test,   & !< cubic harmonic test charge distribution cumulated interaction with source y dipoles
                                                   Ez_test)     !< cubic harmonic test charge distribution cumulated interaction with source z dipoles
  
  use mod_D_X_V_and_E
  use mod_D_X_V_and_E_vec1
  use mod_D_X_V_and_E_vec4
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! input variables
  integer      :: n_source              !< number of source charges and dipoles
  real(kind=8) :: r_source(3,n_source)  !< position of source charges and dipoles
  real(kind=8) :: q_source(n_source)    !< source charge value
  real(kind=8) :: dx_source(n_source)   !< source x dipole value
  real(kind=8) :: dy_source(n_source)   !< source y dipole value
  real(kind=8) :: dz_source(n_source)   !< source z dipole value
  real(kind=8) :: alpha                 !< exponent of solid harmonic test charge distribution
  real(kind=8) :: r(3)                  !< center of solid harmonic test charge distribution
  integer      :: l                     !< max l for hermit harmonic test charge distribution
  
  ! return value
  real(kind=8) :: V_test(*)   !< cubic harmonic test charge distribution cumulated interaction with source charge
  real(kind=8) :: Ex_test(*)  !< cubic harmonic test charge distribution cumulated interaction with source x dipoles
  real(kind=8) :: Ey_test(*)  !< cubic harmonic test charge distribution cumulated interaction with source y dipoles
  real(kind=8) :: Ez_test(*)  !< cubic harmonic test charge distribution cumulated interaction with source z dipoles
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer, parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! temp variables
  integer      :: i
  real(kind=8) :: rc(3)
  real(kind=8) :: rc_v1(1,3)
  real(kind=8) :: rc_v4(4,3)
  real(kind=8) :: wv_v1(1)
  real(kind=8) :: wx_v1(1)
  real(kind=8) :: wy_v1(1)
  real(kind=8) :: wz_v1(1)
  real(kind=8) :: wv_v4(4)
  real(kind=8) :: wx_v4(4)
  real(kind=8) :: wy_v4(4)
  real(kind=8) :: wz_v4(4)
  real(kind=8) :: norm 
  real(kind=8) :: V(455) 
  real(kind=8) :: Ex(455) 
  real(kind=8) :: Ey(455) 
  real(kind=8) :: Ez(455) 
  
  ! get state norm 
  norm = sqrt(PI/alpha)*(PI/alpha)
  
  ! init result
   V_test(1:imax(l))=0.0d0
  Ex_test(1:imax(l))=0.0d0
  Ey_test(1:imax(l))=0.0d0
  Ez_test(1:imax(l))=0.0d0
  
  ! loop on sources
  ! 
  ! uses vectorized routine if possible
  if ( l.le.8 ) then
    ! initial loop
    do i = 1,modulo(n_source,4)
      ! set rc
      rc_v1(1,1) = r_source(1,i)-r(1)
      rc_v1(1,2) = r_source(2,i)-r(2)
      rc_v1(1,3) = r_source(3,i)-r(3)
      ! set w
      wv_v1(1) =  q_source(i)
      wx_v1(1) = dx_source(i)
      wy_v1(1) = dy_source(i)
      wz_v1(1) = dz_source(i)
      ! call non-vectorized routine
      if ( .not.D_X_V_and_E_vec1(alpha, l, rc_v1, wv_v1, wx_v1, wy_v1, wz_v1, V, Ex, Ey, Ez) ) then
        print *,'problem with D_X_V_and_E_vec1 in - cumulative_electrostatics_on_D -'
        stop 1
      end if
      ! cumulate result
       V_test(1:imax(l)) =  V_test(1:imax(l)) +  V(1:imax(l)) * norm
      Ex_test(1:imax(l)) = Ex_test(1:imax(l)) - Ex(1:imax(l)) * norm
      Ey_test(1:imax(l)) = Ey_test(1:imax(l)) - Ey(1:imax(l)) * norm
      Ez_test(1:imax(l)) = Ez_test(1:imax(l)) - Ez(1:imax(l)) * norm
    end do
    ! vectorized loop
    do i = modulo(n_source,4)+1,n_source-1,4
      ! set rc
      rc_v4(:,1) = r_source(1,i:i+3)-r(1)
      rc_v4(:,2) = r_source(2,i:i+3)-r(2)
      rc_v4(:,3) = r_source(3,i:i+3)-r(3)
      ! set w
      wv_v4(1:4) =  q_source(i:i+3)
      wx_v4(1:4) = dx_source(i:i+3)
      wy_v4(1:4) = dy_source(i:i+3)
      wz_v4(1:4) = dz_source(i:i+3)
      ! call vectorized routine
      if ( .not.D_X_V_and_E_vec4(alpha, l, rc_v4, wv_v4, wx_v4, wy_v4, wz_v4, V, Ex, Ey, Ez) ) then
        print *,'problem with D_X_V_and_E_vec4 in - cumulative_electrostatics_on_D -'
        stop 1
      end if
      ! cumulate result
       V_test(1:imax(l)) =  V_test(1:imax(l)) +  V(1:imax(l)) * norm
      Ex_test(1:imax(l)) = Ex_test(1:imax(l)) - Ex(1:imax(l)) * norm
      Ey_test(1:imax(l)) = Ey_test(1:imax(l)) - Ey(1:imax(l)) * norm
      Ez_test(1:imax(l)) = Ez_test(1:imax(l)) - Ez(1:imax(l)) * norm
    end do
  ! 
  ! otherwise, use general purpose routine
  ! 
  else
    
    ! initial loop
    do i = 1,n_source
      ! set rc
      rc(1) = r_source(1,i)-r(1)
      rc(2) = r_source(2,i)-r(2)
      rc(3) = r_source(3,i)-r(3)
      ! call non-vectorized routine
      if ( .not.D_X_V_and_E(alpha, l, rc, V, Ex, Ey, Ez) ) then
        print *,'problem with D_X_V_and_E in - cumulative_electrostatics_on_D -'
        stop 1
      end if
      ! cumulate result
       V_test(1:imax(l)) =  V_test(1:imax(l)) +  V(1:imax(l)) *  q_source(i) * norm
      Ex_test(1:imax(l)) = Ex_test(1:imax(l)) - Ex(1:imax(l)) * dx_source(i) * norm
      Ey_test(1:imax(l)) = Ey_test(1:imax(l)) - Ey(1:imax(l)) * dy_source(i) * norm
      Ez_test(1:imax(l)) = Ez_test(1:imax(l)) - Ez(1:imax(l)) * dz_source(i) * norm
    end do
    
  end if
  
end subroutine







!!> returns total coulomb energy form interaction of a cubic harmonic
!!!  charge distribution and weighted source charges and dipoles 
!!!
recursive subroutine cumulative_electrostatics_on_C(r_source,  & !< source charge and dipole location
                                                   q_source,  & !< source charge value
                                                   dx_source, & !< source x dipole value
                                                   dy_source, & !< source y dipole value
                                                   dz_source, & !< source z dipole value
                                                   n_source,  & !< number of source charges and dipoles
                                                   alpha,     & !< exponent of cubic harmonic test charge distribution
                                                   r,         & !< center of cubic harmonic test charge distribution
                                                   nx,        & !< x power for cubic harmonic test charge distribution
                                                   ny,        & !< y power for cubic harmonic test charge distribution
                                                   nz,        & !< z power for cubic harmonic test charge distribution
                                                   V_test,    & !< cubic harmonic test charge distribution cumulated interaction with source charge
                                                   Ex_test,   & !< solid harmonic test charge distribution cumulated interaction with x source dipoles
                                                   Ey_test,   & !< solid harmonic test charge distribution cumulated interaction with y source dipoles
                                                   Ez_test)     !< solid harmonic test charge distribution cumulated interaction with z source dipoles
  
  use mod_D_X_V_and_E
  use mod_D_X_V_and_E_vec1
  use mod_D_X_V_and_E_vec4
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! input variables
  integer      :: n_source              !< number of source charges and dipoles
  real(kind=8) :: r_source(3,n_source)  !< position of source charges and dipoles
  real(kind=8) :: q_source(n_source)    !< source charge value
  real(kind=8) :: dx_source(n_source)   !< source x dipole value
  real(kind=8) :: dy_source(n_source)   !< source y dipole value
  real(kind=8) :: dz_source(n_source)   !< source z dipole value
  real(kind=8) :: alpha                 !< exponent of solid harmonic test charge distribution
  real(kind=8) :: r(3)                  !< center of solid harmonic test charge distribution
  integer      :: nx                    !< x power for cubic harmonic test charge distribution
  integer      :: ny                    !< y power for cubic harmonic test charge distribution
  integer      :: nz                    !< z power for cubic harmonic test charge distribution
  
  ! return value
  real(kind=8) :: V_test     !< cubic harmonic test charge distribution cumulated interaction with source charge
  real(kind=8) :: Ex_test    !< solid harmonic test charge distribution cumulated interaction with x source dipoles
  real(kind=8) :: Ey_test    !< solid harmonic test charge distribution cumulated interaction with y source dipoles
  real(kind=8) :: Ez_test    !< solid harmonic test charge distribution cumulated interaction with z source dipoles
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer, parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! temp variables
  integer      :: i
  real(kind=8) :: rc(3)
  real(kind=8) :: rc_v1(1,3)
  real(kind=8) :: rc_v4(4,3)
  real(kind=8) :: wv_v1(1)
  real(kind=8) :: wx_v1(1)
  real(kind=8) :: wy_v1(1)
  real(kind=8) :: wz_v1(1)
  real(kind=8) :: wv_v4(4)
  real(kind=8) :: wx_v4(4)
  real(kind=8) :: wy_v4(4)
  real(kind=8) :: wz_v4(4)
  real(kind=8) :: norm 
  real(kind=8) :: V(455) 
  real(kind=8) :: Ex(455) 
  real(kind=8) :: Ey(455) 
  real(kind=8) :: Ez(455) 
  real(kind=8) :: coeffs(455) ! holds the D decomp of source orbital
  
  ! compute contribution from this index
  coeffs=0.0d0
  call R_to_D(alpha,ir_index(nx,ny,nz),coeffs)
  
  ! get state norm 
  norm = sqrt(PI/alpha)*(PI/alpha)
  
  ! call D routine
  call cumulative_electrostatics_on_D(r_source,  &
                                     q_source,  & 
                                     dx_source, & 
                                     dy_source, & 
                                     dz_source, & 
                                     n_source,  &
                                     alpha,     & 
                                     r,         & 
                                     nx+ny+nz,  & 
                                     V,         & 
                                     Ex,        & 
                                     Ey,        & 
                                     Ez)
  
  ! set results 
  V_test    = sum( coeffs(1:imax(nx+ny+nz)) * V(1:imax(nx+ny+nz))  ) ! * norm  now taken into account in cumulative_electrostatics_on_D
  Ex_test = sum( coeffs(1:imax(nx+ny+nz)) * Ex(1:imax(nx+ny+nz)) )   ! * norm 
  Ey_test = sum( coeffs(1:imax(nx+ny+nz)) * Ey(1:imax(nx+ny+nz)) )   ! * norm
  Ez_test = sum( coeffs(1:imax(nx+ny+nz)) * Ez(1:imax(nx+ny+nz)) )   ! * norm
  
!!
!! DIRECT COMPUTATION WITHOUT INVOKING D ROUTINE, SHOULD NOT BE MUCH FASTER
!!
!!  ! init result
!!  V_test=0.0d0
!!  E_test=0.0d0
!!  
!!  ! loop on sources
!!  ! 
!!  ! uses vectorized routine if possible
!!  if ( nx+ny+nz.le.8 ) then
!!    ! initial loop
!!    do i = 1,modulo(n_source,4)
!!      ! set rc
!!      rc_v1(1,1) = r_source(1,i)-r(1)
!!      rc_v1(1,2) = r_source(2,i)-r(2)
!!      rc_v1(1,3) = r_source(3,i)-r(3)
!!      ! set w
!!      wv_v1(1) =  q_source(i)
!!      wx_v1(1) = dx_source(i)
!!      wy_v1(1) = dy_source(i)
!!      wz_v1(1) = dz_source(i)
!!      ! call non-vectorized routine
!!      if ( .not.D_X_V_and_E_vec1(alpha, nx+ny+nz, rc_v1, wv_v1, wx_v1, wy_v1, wz_v1, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E_vec1 in - cumulative_electrostatics_on_C -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      V_test    = V_test   +sum( coeffs(1:imax(nx+ny+nz)) * V(1:imax(nx+ny+nz))  ) * norm
!!      Ex_test = Ex_test-sum( coeffs(1:imax(nx+ny+nz)) * Ex(1:imax(nx+ny+nz)) ) * norm
!!      Ey_test = Ey_test-sum( coeffs(1:imax(nx+ny+nz)) * Ey(1:imax(nx+ny+nz)) ) * norm
!!      Ez_test = Ez_test-sum( coeffs(1:imax(nx+ny+nz)) * Ez(1:imax(nx+ny+nz)) ) * norm
!!    end do
!!    ! vectorized loop
!!    do i = modulo(n_source,4)+1,n_source-1,4
!!      ! set rc
!!      rc_v4(:,1) = r_source(1,i:i+3)-r(1)
!!      rc_v4(:,2) = r_source(2,i:i+3)-r(2)
!!      rc_v4(:,3) = r_source(3,i:i+3)-r(3)
!!      ! set w
!!      wv_v4(1:4) =  q_source(i:i+3)
!!      wx_v4(1:4) = dx_source(i:i+3)
!!      wy_v4(1:4) = dy_source(i:i+3)
!!      wz_v4(1:4) = dz_source(i:i+3)
!!      ! call vectorized routine
!!      if ( .not.D_X_V_and_E_vec4(alpha, nx+ny+nz, rc_v4, wv_v4, wx_v4, wy_v4, wz_v4, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E_vec4 in - cumulative_electrostatics_on_C -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      V_test    = V_test   +sum( coeffs(1:imax(nx+ny+nz)) * V(1:imax(nx+ny+nz))  ) * norm
!!      Ex_test = Ex_test-sum( coeffs(1:imax(nx+ny+nz)) * Ex(1:imax(nx+ny+nz)) ) * norm
!!      Ey_test = Ey_test-sum( coeffs(1:imax(nx+ny+nz)) * Ey(1:imax(nx+ny+nz)) ) * norm
!!      Ez_test = Ez_test-sum( coeffs(1:imax(nx+ny+nz)) * Ez(1:imax(nx+ny+nz)) ) * norm
!!    end do
!!  ! 
!!  ! otherwise, use general purpose routine
!!  ! 
!!  else
!!    
!!    ! initial loop
!!    do i = 1,n_source
!!      ! set rc
!!      rc(1) = r_source(1,i)-r(1)
!!      rc(2) = r_source(2,i)-r(2)
!!      rc(3) = r_source(3,i)-r(3)
!!      ! call non-vectorized routine
!!      if ( .not.D_X_V_and_E(alpha, nx+ny+nz, rc, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E in - cumulative_electrostatics_on_C -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      V_test    = V_test   +sum( coeffs(1:imax(nx+ny+nz)) * V(1:imax(nx+ny+nz))  ) * norm *  q_source(i)
!!      Ex_test = Ex_test-sum( coeffs(1:imax(nx+ny+nz)) * Ex(1:imax(nx+ny+nz)) ) * norm * dx_source(i)
!!      Ey_test = Ey_test-sum( coeffs(1:imax(nx+ny+nz)) * Ey(1:imax(nx+ny+nz)) ) * norm * dy_source(i)
!!      Ez_test = Ez_test-sum( coeffs(1:imax(nx+ny+nz)) * Ez(1:imax(nx+ny+nz)) ) * norm * dz_source(i)
!!    end do
!!    
!!  end if
  
end subroutine





!!> returns total coulomb energy form interaction of a cubic harmonic
!!!  charge distribution and weighted source charges and dipoles 
!!!
recursive subroutine cumulative_electrostatics_on_Y(r_source,  & !< source charge and dipole location
                                                   q_source,  & !< source charge value
                                                   dx_source, & !< source x dipole value
                                                   dy_source, & !< source y dipole value
                                                   dz_source, & !< source z dipole value
                                                   n_source,  & !< number of source charges and dipoles
                                                   alpha,     & !< exponent of solid harmonic test charge distribution
                                                   r,         & !< center of solid harmonic test charge distribution
                                                   l,         & !< l momentum of solid harmonic test charge distribution
                                                   m,         & !< m momentum of solid harmonic test charge distribution
                                                   V_test,    & !< solid harmonic test charge distribution cumulated interaction with source charge
                                                   Ex_test,   & !< solid harmonic test charge distribution cumulated interaction with x source dipoles
                                                   Ey_test,   & !< solid harmonic test charge distribution cumulated interaction with y source dipoles
                                                   Ez_test)     !< solid harmonic test charge distribution cumulated interaction with z source dipoles
  
  use mod_D_X_V_and_E
  use mod_D_X_V_and_E_vec1
  use mod_D_X_V_and_E_vec4
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! input variables
  integer      :: n_source              !< number of source charges and dipoles
  real(kind=8) :: r_source(3,n_source)  !< position of source charges and dipoles
  real(kind=8) :: q_source(n_source)    !< source charge value
  real(kind=8) :: dx_source(n_source)   !< source x dipole value
  real(kind=8) :: dy_source(n_source)   !< source y dipole value
  real(kind=8) :: dz_source(n_source)   !< source z dipole value
  real(kind=8) :: alpha                 !< exponent of solid harmonic test charge distribution
  real(kind=8) :: r(3)                  !< center of solid harmonic test charge distribution
  integer      :: l                     !< l momentum of solid harmonic test charge distribution
  integer      :: m                     !< m momentum of solid harmonic test charge distribution
  
  ! return value
  real(kind=8) :: V_test     !< solid harmonic test charge distribution cumulated interaction with source charge
  real(kind=8) :: Ex_test    !< solid harmonic test charge distribution cumulated interaction with x source dipoles
  real(kind=8) :: Ey_test    !< solid harmonic test charge distribution cumulated interaction with y source dipoles
  real(kind=8) :: Ez_test    !< solid harmonic test charge distribution cumulated interaction with z source dipoles
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer, parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! temp variables
  integer      :: i
  integer      :: j
  real(kind=8) :: rc(3)
  real(kind=8) :: rc_v1(1,3)
  real(kind=8) :: rc_v4(4,3)
  real(kind=8) :: wv_v1(1)
  real(kind=8) :: wx_v1(1)
  real(kind=8) :: wy_v1(1)
  real(kind=8) :: wz_v1(1)
  real(kind=8) :: wv_v4(4)
  real(kind=8) :: wx_v4(4)
  real(kind=8) :: wy_v4(4)
  real(kind=8) :: wz_v4(4)
  real(kind=8) :: norm 
  real(kind=8) :: V(455) 
  real(kind=8) :: Ex(455) 
  real(kind=8) :: Ey(455) 
  real(kind=8) :: Ez(455) 
  real(kind=8) :: coeffs(455) ! holds the D decomp of source orbital
  real(kind=8) :: rD(3)       ! center for the derivatives basis
  real(kind=8) :: alphaD      ! exponent for the derivatives basis
  integer      :: iD(455)     ! corresponding derivative index array
  real(kind=8) :: cD(455)     ! corresponding derivative coefficient array
  integer      :: nD          ! number of corresponding derivatives
  
  ! compute contribution from this index
  call Y_to_D(alpha,r,l,m,alphaD,rD,cD,iD,nD)
  
  ! get state norm 
  norm = sqrt(PI/alphaD)*(PI/alphaD)
  
  ! call D routine
  call cumulative_electrostatics_on_D(r_source,  &
                                     q_source,  & 
                                     dx_source, & 
                                     dy_source, & 
                                     dz_source, & 
                                     n_source,  &
                                     alphaD,    & 
                                     rD,        & 
                                     l,         & 
                                     V,         & 
                                     Ex,        & 
                                     Ey,        & 
                                     Ez)
  
  ! set results 
  V_test=0.0d0
  Ex_test=0.0d0
  Ey_test=0.0d0
  Ez_test=0.0d0
  do j=1,nD
    V_test    = V_test    + cD(j) *  V( iD(j) ) ! * norm  now taken into account in cumulative_electrostatics_on_D
    Ex_test = Ex_test + cD(j) * Ex( iD(j) )     ! * norm
    Ey_test = Ey_test + cD(j) * Ey( iD(j) )     ! * norm
    Ez_test = Ez_test + cD(j) * Ez( iD(j) )     ! * norm
  end do
  
  
!!
!! DIRECT COMPUTATION WITHOUT INVOKING D ROUTINE, SHOULD NOT BE MUCH FASTER
!!
!!  ! init result
!!  V_test=0.0d0
!!  E_test=0.0d0
!!  
!!  ! loop on sources
!!  ! 
!!  ! uses vectorized routine if possible
!!  if ( l.le.8 ) then
!!    ! initial loop
!!    do i = 1,modulo(n_source,4)
!!      ! set rc
!!      rc_v1(1,1) = r_source(1,i)-rD(1)
!!      rc_v1(1,2) = r_source(2,i)-rD(2)
!!      rc_v1(1,3) = r_source(3,i)-rD(3)
!!      ! set w
!!      wv_v1(1) =  q_source(i)
!!      wx_v1(1) = dx_source(i)
!!      wy_v1(1) = dy_source(i)
!!      wz_v1(1) = dz_source(i)
!!      ! call non-vectorized routine
!!      if ( .not.D_X_V_and_E_vec1(alphaD, l, rc_v1, wv_v1, wx_v1, wy_v1, wz_v1, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E_vec1 in - cumulative_electrostatics_on_Y -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      do j=1,nD
!!        V_test    = V_test    + cD(j) *  V( iD(j) ) * norm
!!        Ex_test = Ex_test - cD(j) * Ex( iD(j) ) * norm
!!        Ey_test = Ey_test - cD(j) * Ey( iD(j) ) * norm
!!        Ez_test = Ez_test - cD(j) * Ez( iD(j) ) * norm
!!      end do
!!    end do
!!    ! vectorized loop
!!    do i = modulo(n_source,4)+1,n_source-1,4
!!      ! set rc
!!      rc_v4(:,1) = r_source(1,i:i+3)-rD(1)
!!      rc_v4(:,2) = r_source(2,i:i+3)-rD(2)
!!      rc_v4(:,3) = r_source(3,i:i+3)-rD(3)
!!      ! set w
!!      wv_v4(1:4) =  q_source(i:i+3)
!!      wx_v4(1:4) = dx_source(i:i+3)
!!      wy_v4(1:4) = dy_source(i:i+3)
!!      wz_v4(1:4) = dz_source(i:i+3)
!!      ! call vectorized routine
!!      if ( .not.D_X_V_and_E_vec4(alphaD, l, rc_v4, wv_v4, wx_v4, wy_v4, wz_v4, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E_vec4 in - cumulative_electrostatics_on_Y -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      do j=1,nD
!!        V_test    = V_test    + cD(j) *  V( iD(j) ) * norm
!!        Ex_test = Ex_test - cD(j) * Ex( iD(j) ) * norm
!!        Ey_test = Ey_test - cD(j) * Ey( iD(j) ) * norm
!!        Ez_test = Ez_test - cD(j) * Ez( iD(j) ) * norm
!!      end do
!!    end do
!!  ! 
!!  ! otherwise, use general purpose routine
!!  ! 
!!  else
!!    
!!    ! initial loop
!!    do i = 1,n_source
!!      ! set rc
!!      rc(1) = r_source(1,i)-rD(1)
!!      rc(2) = r_source(2,i)-rD(2)
!!      rc(3) = r_source(3,i)-rD(3)
!!      ! call non-vectorized routine
!!      if ( .not.D_X_V_and_E(alphaD, l, rc, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E in - cumulative_electrostatics_on_Y -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      do j=1,nD
!!        V_test    = V_test    + cD(j) *  V( iD(j) ) * norm *  q_source(i)
!!        Ex_test = Ex_test - cD(j) * Ex( iD(j) ) * norm * dx_source(i)
!!        Ey_test = Ey_test - cD(j) * Ey( iD(j) ) * norm * dy_source(i)
!!        Ez_test = Ez_test - cD(j) * Ez( iD(j) ) * norm * dz_source(i)
!!      end do
!!    end do
!!    
!!  end if
  
end subroutine






!!> returns total coulomb energy form interaction of a cubic harmonic
!!!  charge distribution and weighted source charges and dipoles 
!!!
recursive subroutine cumulative_electrostatics_on_YY(r_source,  & !< source charge and dipole location
                                                    q_source,  & !< source charge value
                                                    dx_source, & !< source x dipole value
                                                    dy_source, & !< source y dipole value
                                                    dz_source, & !< source z dipole value
                                                    n_source,  & !< number of source charges and dipoles
                                                    alpha1,    & !< exponent of first solid harmonic test charge distribution
                                                    r1,        & !< center of first solid harmonic test charge distribution
                                                    l1,        & !< l momentum of first solid harmonic test charge distribution
                                                    m1,        & !< m momentum of first solid harmonic test charge distribution
                                                    alpha2,    & !< exponent of second solid harmonic test charge distribution
                                                    r2,        & !< center of second solid harmonic test charge distribution
                                                    l2,        & !< l momentum of second solid harmonic test charge distribution
                                                    m2,        & !< m momentum of second solid harmonic test charge distribution
                                                    V_test,    & !< solid harmonic test charge distribution cumulated interaction with source charge
                                                    Ex_test,   & !< solid harmonic test charge distribution cumulated interaction with x source dipoles
                                                    Ey_test,   & !< solid harmonic test charge distribution cumulated interaction with y source dipoles
                                                    Ez_test)     !< solid harmonic test charge distribution cumulated interaction with z source dipoles
  
  use mod_D_X_V_and_E
  use mod_D_X_V_and_E_vec1
  use mod_D_X_V_and_E_vec4
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! input variables
  integer      :: n_source              !< number of source charges and dipoles
  real(kind=8) :: r_source(3,n_source)  !< position of source charges and dipoles
  real(kind=8) :: q_source(n_source)    !< source charge value
  real(kind=8) :: dx_source(n_source)   !< source x dipole value
  real(kind=8) :: dy_source(n_source)   !< source y dipole value
  real(kind=8) :: dz_source(n_source)   !< source z dipole value
  real(kind=8) :: alpha1                !< exponent of solid harmonic test charge distribution
  real(kind=8) :: r1(3)                 !< center of solid harmonic test charge distribution
  integer      :: l1                    !< l momentum of solid harmonic test charge distribution
  integer      :: m1                    !< m momentum of solid harmonic test charge distribution
  real(kind=8) :: alpha2                !< exponent of solid harmonic test charge distribution
  real(kind=8) :: r2(3)                 !< center of solid harmonic test charge distribution
  integer      :: l2                    !< l momentum of solid harmonic test charge distribution
  integer      :: m2                    !< m momentum of solid harmonic test charge distribution
  
  ! return value
  real(kind=8) :: V_test     !< solid harmonic test charge distribution cumulated interaction with source charge
  real(kind=8) :: Ex_test    !< solid harmonic test charge distribution cumulated interaction with x source dipoles
  real(kind=8) :: Ey_test    !< solid harmonic test charge distribution cumulated interaction with y source dipoles
  real(kind=8) :: Ez_test    !< solid harmonic test charge distribution cumulated interaction with z source dipoles
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer, parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! temp variables
  integer      :: i
  integer      :: j
  real(kind=8) :: rc(3)
  real(kind=8) :: rc_v1(1,3)
  real(kind=8) :: rc_v4(4,3)
  real(kind=8) :: wv_v1(1)
  real(kind=8) :: wx_v1(1)
  real(kind=8) :: wy_v1(1)
  real(kind=8) :: wz_v1(1)
  real(kind=8) :: wv_v4(4)
  real(kind=8) :: wx_v4(4)
  real(kind=8) :: wy_v4(4)
  real(kind=8) :: wz_v4(4)
  real(kind=8) :: norm 
  real(kind=8) :: V(455) 
  real(kind=8) :: Ex(455) 
  real(kind=8) :: Ey(455) 
  real(kind=8) :: Ez(455) 
  real(kind=8) :: coeffs(455) ! holds the D decomp of source orbital
  real(kind=8) :: rD(3)       ! center for the derivatives basis
  real(kind=8) :: alphaD      ! exponent for the derivatives basis
  integer      :: iD(455)     ! corresponding derivative index array
  real(kind=8) :: cD(455)     ! corresponding derivative coefficient array
  integer      :: nD          ! number of corresponding derivatives
  
  ! compute contribution from this index
  call YY_to_D(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alphaD,rD,cD,iD,nD)
  
  ! get state norm 
  norm = sqrt(PI/alphaD)*(PI/alphaD)
  
  ! call D routine
  call cumulative_electrostatics_on_D(r_source,  &
                                     q_source,  & 
                                     dx_source, & 
                                     dy_source, & 
                                     dz_source, & 
                                     n_source,  &
                                     alphaD,    & 
                                     rD,        & 
                                     l1+l2,     & 
                                     V,         & 
                                     Ex,        & 
                                     Ey,        & 
                                     Ez)
  
  ! set results 
  V_test=0.0d0
  Ex_test=0.0d0
  Ey_test=0.0d0
  Ez_test=0.0d0
  do j=1,nD
    V_test    = V_test    + cD(j) *  V( iD(j) ) ! * norm  now taken into account in cumulative_electrostatics_on_D
    Ex_test = Ex_test + cD(j) * Ex( iD(j) )     ! * norm
    Ey_test = Ey_test + cD(j) * Ey( iD(j) )     ! * norm
    Ez_test = Ez_test + cD(j) * Ez( iD(j) )     ! * norm
  end do
  
  
!!
!! DIRECT COMPUTATION WITHOUT INVOKING D ROUTINE, SHOULD NOT BE MUCH FASTER
!!
!!  ! init result
!!  V_test=0.0d0
!!  E_test=0.0d0
!!  
!!  ! loop on sources
!!  ! 
!!  ! uses vectorized routine if possible
!!  if ( l.le.8 ) then
!!    ! initial loop
!!    do i = 1,modulo(n_source,4)
!!      ! set rc
!!      rc_v1(1,1) = r_source(1,i)-rD(1)
!!      rc_v1(1,2) = r_source(2,i)-rD(2)
!!      rc_v1(1,3) = r_source(3,i)-rD(3)
!!      ! set w
!!      wv_v1(1) =  q_source(i)
!!      wx_v1(1) = dx_source(i)
!!      wy_v1(1) = dy_source(i)
!!      wz_v1(1) = dz_source(i)
!!      ! call non-vectorized routine
!!      if ( .not.D_X_V_and_E_vec1(alphaD, l, rc_v1, wv_v1, wx_v1, wy_v1, wz_v1, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E_vec1 in - cumulative_electrostatics_on_Y -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      do j=1,nD
!!        V_test    = V_test    + cD(j) *  V( iD(j) ) * norm
!!        Ex_test = Ex_test - cD(j) * Ex( iD(j) ) * norm
!!        Ey_test = Ey_test - cD(j) * Ey( iD(j) ) * norm
!!        Ez_test = Ez_test - cD(j) * Ez( iD(j) ) * norm
!!      end do
!!    end do
!!    ! vectorized loop
!!    do i = modulo(n_source,4)+1,n_source-1,4
!!      ! set rc
!!      rc_v4(:,1) = r_source(1,i:i+3)-rD(1)
!!      rc_v4(:,2) = r_source(2,i:i+3)-rD(2)
!!      rc_v4(:,3) = r_source(3,i:i+3)-rD(3)
!!      ! set w
!!      wv_v4(1:4) =  q_source(i:i+3)
!!      wx_v4(1:4) = dx_source(i:i+3)
!!      wy_v4(1:4) = dy_source(i:i+3)
!!      wz_v4(1:4) = dz_source(i:i+3)
!!      ! call vectorized routine
!!      if ( .not.D_X_V_and_E_vec4(alphaD, l, rc_v4, wv_v4, wx_v4, wy_v4, wz_v4, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E_vec4 in - cumulative_electrostatics_on_Y -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      do j=1,nD
!!        V_test    = V_test    + cD(j) *  V( iD(j) ) * norm
!!        Ex_test = Ex_test - cD(j) * Ex( iD(j) ) * norm
!!        Ey_test = Ey_test - cD(j) * Ey( iD(j) ) * norm
!!        Ez_test = Ez_test - cD(j) * Ez( iD(j) ) * norm
!!      end do
!!    end do
!!  ! 
!!  ! otherwise, use general purpose routine
!!  ! 
!!  else
!!    
!!    ! initial loop
!!    do i = 1,n_source
!!      ! set rc
!!      rc(1) = r_source(1,i)-rD(1)
!!      rc(2) = r_source(2,i)-rD(2)
!!      rc(3) = r_source(3,i)-rD(3)
!!      ! call non-vectorized routine
!!      if ( .not.D_X_V_and_E(alphaD, l, rc, V, Ex, Ey, Ez) ) then
!!        print *,'problem with D_X_V_and_E in - cumulative_electrostatics_on_Y -'
!!        stop 1
!!      end if
!!      ! cumulate result
!!      do j=1,nD
!!        V_test    = V_test    + cD(j) *  V( iD(j) ) * norm *  q_source(i)
!!        Ex_test = Ex_test - cD(j) * Ex( iD(j) ) * norm * dx_source(i)
!!        Ey_test = Ey_test - cD(j) * Ey( iD(j) ) * norm * dy_source(i)
!!        Ez_test = Ez_test - cD(j) * Ez( iD(j) ) * norm * dz_source(i)
!!      end do
!!    end do
!!    
!!  end if
  
end subroutine



