
!> @file Coulomb.f90
!!
!! Defines utility routines for computing coulomb 
!! integrals between gaussian basis elements.
!!
!! Author: I. Duchemin July 2015
!!
recursive function R_X_R(alpha1,r1,c1,l1max,alpha2,r2,c2,l2max)
  
  use mod_CoulombUtils
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)   !< center for first cubic Harmonic
  real(kind=8) :: r2(3)   !< center for second cubic Harmonic
  real(kind=8) :: alpha1  !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2  !< exponent for second cubic Harmonic
  integer      :: l1max   !< max order for first cubic Harmonic
  integer      :: l2max   !< max order for second cubic Harmonic
  real(kind=8) :: c1(455) !< coefficients of first cubic Harmonic in the R basis
  real(kind=8) :: c2(455) !< coefficients of second cubic Harmonic in the R basis
  
  ! return value
  real(kind=8)  :: R_X_R
  
  ! local variables
  real(kind=8) :: coeffs1(455)
  real(kind=8) :: coeffs2(455)
  real(kind=8) :: coeffs_tmp(455)
  real(kind=8) :: x,y,z,q
  real(kind=8) :: p,d,erfpd,exppd
  integer      :: il,ir
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer     , parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! form coeffs for derivatives on the left side
  coeffs1=0.0d0
  do il=1,imax(l1max)
    if ( c1(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha1,il,coeffs_tmp)
      ! add contrib
      coeffs1(1:il)=coeffs1(1:il)+c1(il)*coeffs_tmp(1:il)
    end if
  end do
  
  ! form coeffs for derivatives on the right side
  coeffs2=0.0d0
  do ir=1,imax(l2max)
    if ( c2(ir).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha2,ir,coeffs_tmp)
      ! add contrib
      coeffs2(1:ir)=coeffs2(1:ir)+c2(ir)*coeffs_tmp(1:ir)
    end if
  end do
  
  ! compute combined exponant
  q=alpha1*alpha2/(alpha1+alpha2)
  
  ! compute center distance
  x=r2(1)-r1(1)
  y=r2(2)-r1(2)
  z=r2(3)-r1(3)
  
  ! compute once and for all p=sqrt(q), d, erf(p*d) and exp(-p^2*d^2)
  p=sqrt(q)
  d=sqrt(x**2+y**2+z**2)
  erfpd=erf(p*d)
  exppd=exp(-p**2*d**2)
  
  ! loop on the components
  R_X_R=0.0d0
  do il=1,imax(l1max)
    if ( coeffs1(il).ne.0.0d0 ) then
      do ir=1,imax(l2max)
        if ( coeffs2(ir).ne.0.0d0 ) then
          
!          print *,il,ir,coeffs1(il),coeffs2(ir),D_X_D(il,ir,x,y,z,q,p,d,erfpd,exppd)
          
          R_X_R =R_X_R + coeffs1(il)*coeffs2(ir)*D_X_D(il,ir,x,y,z,q,p,d,erfpd,exppd)
          
!          print *,R_X_R*sqrt(PI/alpha1)*(PI/alpha1)*sqrt(PI/alpha2)*(PI/alpha2)
        end if
      end do
    end if
  end do
  
  ! norm and return value
  R_X_R =R_X_R*sqrt(PI/alpha1)*(PI/alpha1)*sqrt(PI/alpha2)*(PI/alpha2)
  
end function
  
!!> attenuated coulomb integral i.e. using erf(omega*r)/r instead of 1/r kernel
!!
recursive function attenuated_R_X_R(alpha1,r1,c1,l1max,alpha2,r2,c2,l2max,omega)
  
  use mod_CoulombUtils
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)   !< center for first cubic Harmonic
  real(kind=8) :: r2(3)   !< center for second cubic Harmonic
  real(kind=8) :: alpha1  !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2  !< exponent for second cubic Harmonic
  real(kind=8) :: omega   !< exponent for second cubic Harmonic
  integer      :: l1max   !< max order for first cubic Harmonic
  integer      :: l2max   !< max order for second cubic Harmonic
  real(kind=8) :: c1(455) !< coefficients of first cubic Harmonic in the R basis
  real(kind=8) :: c2(455) !< coefficients of second cubic Harmonic in the R basis
  
  ! return value
  real(kind=8)  :: attenuated_R_X_R
  
  ! local variables
  real(kind=8) :: coeffs1(455)
  real(kind=8) :: coeffs2(455)
  real(kind=8) :: coeffs_tmp(455)
  real(kind=8) :: x,y,z,q
  real(kind=8) :: p,d,erfpd,exppd
  integer      :: il,ir
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer     , parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! form coeffs for derivatives on the left side
  coeffs1=0.0d0
  do il=1,imax(l1max)
    if ( c1(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha1,il,coeffs_tmp)
      ! add contrib
      coeffs1(1:il)=coeffs1(1:il)+c1(il)*coeffs_tmp(1:il)
    end if
  end do
  
  ! form coeffs for derivatives on the right side
  coeffs2=0.0d0
  do ir=1,imax(l2max)
    if ( c2(ir).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha2,ir,coeffs_tmp)
      ! add contrib
      coeffs2(1:ir)=coeffs2(1:ir)+c2(ir)*coeffs_tmp(1:ir)
    end if
  end do
  
  ! compute combined exponent
  q=(alpha1*alpha2*omega**2)/((alpha1+alpha2)*omega**2+alpha1*alpha2)
  
  ! compute center distance
  x=r2(1)-r1(1)
  y=r2(2)-r1(2)
  z=r2(3)-r1(3)
  
  ! compute once and for all p=sqrt(q), d, erf(p*d) and exp(-p^2*d^2)
  p=sqrt(q)
  d=sqrt(x**2+y**2+z**2)
  erfpd=erf(p*d)
  exppd=exp(-p**2*d**2)
  
  ! loop on the components
  attenuated_R_X_R=0.0d0
  do il=1,imax(l1max)
    if ( coeffs1(il).ne.0.0d0 ) then
      do ir=1,imax(l2max)
        if ( coeffs2(ir).ne.0.0d0 ) then
          attenuated_R_X_R =attenuated_R_X_R + coeffs1(il)*coeffs2(ir)*D_X_D(il,ir,x,y,z,q,p,d,erfpd,exppd)
        end if
      end do
    end if
  end do
  
  ! norm and return value
  attenuated_R_X_R =attenuated_R_X_R*sqrt(PI/alpha1)*(PI/alpha1)*sqrt(PI/alpha2)*(PI/alpha2)
  
end function



!> @file Coulomb.f90
!!
!! Computing coulomb integrals between gaussian basis elements, defined as 
!! derivatives of (s||s) kernel.
!!
!! Author: I. Duchemin July 2015
!!
recursive function D_Coulomb_D(alphaD1,rD1,cD1,iD1,nD1,alphaD2,rD2,cD2,iD2,nD2)
   
  use mod_CoulombUtilsExpert
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: rD1(3)    !< center for the derivatives basis
  real(kind=8) :: alphaD1   !< exponent for the derivatives basis
  integer      :: iD1(*)    !< corresponding derivative index array
  real(kind=8) :: cD1(*)    !< corresponding derivative coefficient array
  integer      :: nD1       !< number of corresponding derivatives
  real(kind=8) :: rD2(3)    !< center for the derivatives basis
  real(kind=8) :: alphaD2   !< exponent for the derivatives basis
  integer      :: iD2(*)    !< corresponding derivative index array
  real(kind=8) :: cD2(*)    !< corresponding derivative coefficient array
  integer      :: nD2       !< number of corresponding derivatives
  
  ! return value
  real(kind=8)  :: D_Coulomb_D
  
  ! local variables
  real(kind=8) :: x,y,z,q
  real(kind=8) :: p,d,erfpd,exppd
  integer      :: il,ir
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  
  ! compute combined exponant
  q=alphaD1*alphaD2/(alphaD1+alphaD2)
  
  ! compute center distance
  x=rD2(1)-rD1(1)
  y=rD2(2)-rD1(2)
  z=rD2(3)-rD1(3)
  
  ! compute once and for all p=sqrt(q), d, erf(p*d) and exp(-p^2*d^2)
  p=sqrt(q)
  d=sqrt(x**2+y**2+z**2)
  erfpd=erf(p*d)
  exppd=exp(-p**2*d**2)
  
  ! loop on the components
  D_Coulomb_D = D_X_D_expert_util(cD1,iD1,nD1,cD2,iD2,nD2,x,y,z,q,p,d,erfpd,exppd)
  
  ! norm and return value
  D_Coulomb_D = D_Coulomb_D * sqrt(PI*PI/(alphaD1*alphaD2)) * PI*PI/(alphaD1*alphaD2)
  
end function


!> @file Coulomb.f90
!!
!! Computing coulomb integrals between gaussian basis elements, defined as 
!! derivatives of (s||s) kernel.
!!
!! Author: I. Duchemin July 2015
!!
recursive function D_Coulomb_D_unpacked(alphaD1,rD1,cD1,lD1,alphaD2,rD2,cD2,lD2)
   
  use mod_CoulombUtils
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: rD1(3)    !< center for the derivatives basis
  real(kind=8) :: alphaD1   !< exponent for the derivatives basis
  real(kind=8) :: cD1(*)    !< corresponding derivative coefficient array
  integer      :: lD1       !< max l momentum
  real(kind=8) :: rD2(3)    !< center for the derivatives basis
  real(kind=8) :: alphaD2   !< exponent for the derivatives basis
  real(kind=8) :: cD2(*)    !< corresponding derivative coefficient array
  integer      :: lD2       !< max l momentum
  
  ! return value
  real(kind=8)  :: D_Coulomb_D_unpacked
  
  ! local variables
  real(kind=8) :: x,y,z,q
  real(kind=8) :: p,d,erfpd,exppd
  integer      :: il,ir
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer     , parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! compute combined exponant
  q=alphaD1*alphaD2/(alphaD1+alphaD2)
  
  ! compute center distance
  x=rD2(1)-rD1(1)
  y=rD2(2)-rD1(2)
  z=rD2(3)-rD1(3)
  
  ! compute once and for all p=sqrt(q), d, erf(p*d) and exp(-p^2*d^2)
  p=sqrt(q)
  d=sqrt(x**2+y**2+z**2)
  erfpd=erf(p*d)
  exppd=exp(-p**2*d**2)
  
  ! loop on the components
  D_Coulomb_D_unpacked = 0.0d0
  do il=1,imax(lD1)
    if ( cD1(il).ne.0.0d0 ) then
      do ir=1,imax(lD2)
        if ( cD2(ir).ne.0.0d0 ) then
          D_Coulomb_D_unpacked = D_Coulomb_D_unpacked + cD1(il)*cD2(ir)*D_X_D(il,ir,x,y,z,q,p,d,erfpd,exppd)
        end if
      end do
    end if
  end do
  
  ! norm and return value
  D_Coulomb_D_unpacked = D_Coulomb_D_unpacked * sqrt(PI*PI/(alphaD1*alphaD2)) * PI*PI/(alphaD1*alphaD2)
  
end function


!> @file Coulomb.f90
!!
!! Computing coulomb integrals between gaussian basis elements, defined as 
!! derivatives of (s||s) kernel.
!!
!! Author: I. Duchemin December 2017
!!
recursive function D_Coulomb_Y(alphaD,rD,cD,lD,alpha,r,l,m,value)
  
  use mod_D_X_Y_opt_0
  use mod_D_X_Y_opt_1
  use mod_D_X_Y_opt_2
  use mod_D_X_Y_opt_3
  use mod_D_X_Y_opt_4
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha   !< exponent for third spherical Harmonic
  real(kind=8) :: value   !< value of the three center integral upon return 
  integer      :: lD      !< max angular momentum for hermit decomposition
  integer      :: l       !< l momenta of solid Harmonic
  integer      :: m       !< m momenta of solid Harmonic
  real(kind=8) :: rD(3)  
  real(kind=8) :: alphaD 
  real(kind=8) :: cD(*)
  
  ! return value
  logical      :: D_Coulomb_Y
  
  ! local variables
  real(kind=8) :: rc(3)  
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  
  ! center
  rc = r - rD
  
  ! fast return if possible
  if ( lD>6 ) then
    D_Coulomb_Y = .false.
    return
  end if
  
  ! call derivative routine
  select case (l)
    case (0)
      D_Coulomb_Y = D_X_Y_opt_0(alphaD,alpha,rc,m,cD,lD,value)
    case (1)
      D_Coulomb_Y = D_X_Y_opt_1(alphaD,alpha,rc,m,cD,lD,value)
    case (2)
      D_Coulomb_Y = D_X_Y_opt_2(alphaD,alpha,rc,m,cD,lD,value)
    case (3)
      D_Coulomb_Y = D_X_Y_opt_3(alphaD,alpha,rc,m,cD,lD,value)
    case (4)
      D_Coulomb_Y = D_X_Y_opt_4(alphaD,alpha,rc,m,cD,lD,value)
    case default
      D_Coulomb_Y = .false.
      return 
  end select
  
  ! norm
  value =  value * sqrt(PI*PI/(alphaD*alpha)) * PI*PI/(alphaD*alpha)
  
end function



!> @file Coulomb.f90
!!
!! Computing coulomb integrals between gaussian basis elements, defined as 
!! derivatives of (s||s) kernel.
!!
!! Author: I. Duchemin December 2017
!!
recursive function D_Coulomb_C(alphaD,rD,cD,lD,alpha,r,nx,ny,nz,value)
  
  use mod_D_X_R_opt_0
  use mod_D_X_R_opt_1
  use mod_D_X_R_opt_2
  use mod_D_X_R_opt_3
  use mod_D_X_R_opt_4
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha   !< exponent for third spherical Harmonic
  real(kind=8) :: value   !< value of the three center integral upon return
  integer      :: lD      !< max angular momentum for hermit decomposition
  integer      :: nx      !< x exponent of cubic Harmonic
  integer      :: ny      !< y exponent of cubic Harmonic
  integer      :: nz      !< z exponent of cubic Harmonic
  real(kind=8) :: rD(3)  
  real(kind=8) :: alphaD 
  real(kind=8) :: cD(*)
  
  ! return value
  logical      :: D_Coulomb_C
  
  ! local variables
  real(kind=8) :: rc(3)  
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  
  ! center
  rc = r - rD
  
  ! fast return if possible
  if ( lD>6 ) then
    D_Coulomb_C = .false.
    return
  end if
  
  ! call derivative routine
  value=0.0d0
  select case (nx+ny+nz)
    case (0)
      D_Coulomb_C = D_X_R_opt_0(alphaD,alpha,rc,ny,nz,cD,lD,value)
    case (1)
      D_Coulomb_C = D_X_R_opt_1(alphaD,alpha,rc,ny,nz,cD,lD,value)
    case (2)
      D_Coulomb_C = D_X_R_opt_2(alphaD,alpha,rc,ny,nz,cD,lD,value)
    case (3)
      D_Coulomb_C = D_X_R_opt_3(alphaD,alpha,rc,ny,nz,cD,lD,value)
    case (4)
      D_Coulomb_C = D_X_R_opt_4(alphaD,alpha,rc,ny,nz,cD,lD,value)
    case default
      D_Coulomb_C = .false.
      return 
  end select
  
  ! norm
  value =  value * sqrt(PI*PI/(alphaD*alpha)) * PI*PI/(alphaD*alpha)
  
end function



!> @file Coulomb.f90
!!
!! Computing coulomb integrals between gaussian basis elements, defined as 
!! derivatives of (s||s) kernel.
!!
!! Author: I. Duchemin December 2017
!!
recursive function D_Coulomb_Y_shell(alphaD,rD,cD,lD,alpha,r,l,value)
  
  use mod_D_X_Y_shell_opt_0
  use mod_D_X_Y_shell_opt_1
  use mod_D_X_Y_shell_opt_2
  use mod_D_X_Y_shell_opt_3
  use mod_D_X_Y_shell_opt_4
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r(3)     !< center for third spherical Harmonic
  real(kind=8) :: alpha    !< exponent for third spherical Harmonic
  real(kind=8) :: value(*) !< value of the three center integral upon return 
  integer      :: lD       !< angular momentum for first spherical Harmonic
  integer      :: l        !< l momenta of cubic Harmonic shell
  real(kind=8) :: rD(3)
  real(kind=8) :: alphaD
  real(kind=8) :: cD(*)
  
  ! return value
  logical :: D_Coulomb_Y_shell
  
  ! local variables
  real(kind=8) :: rc(3)
  integer      :: is
  integer      :: m
  real(kind=8) :: cr(455)
  real(kind=8) :: D_Coulomb_D_unpacked
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  
  ! center rD
  rc = r - rD
  
  ! call derivative routine
  select case (l)
    case (0)
      D_Coulomb_Y_shell = D_X_Y_shell_opt_0(alphaD,alpha,rc,cD,lD,value)
    case (1)
      D_Coulomb_Y_shell = D_X_Y_shell_opt_1(alphaD,alpha,rc,cD,lD,value)
    case (2)
      D_Coulomb_Y_shell = D_X_Y_shell_opt_2(alphaD,alpha,rc,cD,lD,value)
    case (3)
      D_Coulomb_Y_shell = D_X_Y_shell_opt_3(alphaD,alpha,rc,cD,lD,value)
    case (4)
      D_Coulomb_Y_shell = D_X_Y_shell_opt_4(alphaD,alpha,rc,cD,lD,value)
    case default
      D_Coulomb_Y_shell = .false.
  end select
  
  ! norm
  value(1:(l+1)*(l+2)/2) =  value(1:(l+1)*(l+2)/2) * sqrt(PI*PI/(alphaD*alpha)) * PI*PI/(alphaD*alpha)
  
  ! fall back on standard procedure if not possible...
  if ( .not.D_Coulomb_Y_shell ) then
    
    D_Coulomb_Y_shell = .true.
    
    ! loop on shell components
    is=0
    do m=-l,l
      
      ! increment shell component index
      is=is+1
      
      ! compute contribution from this index
      call Y_to_D_unpacked(alpha,l,m,cr)
      
      ! compute integral in D rep
      value(is) = D_Coulomb_D_unpacked(alphaD,rD,cD,lD,alpha,r,cr,l)
      
    end do
    
  end if
  
end function





!> @file Coulomb.f90
!!
!! Computing coulomb integrals between gaussian basis elements, defined as 
!! derivatives of (s||s) kernel.
!!
!! Author: I. Duchemin December 2017
!!
recursive function D_Coulomb_C_shell(alphaD,rD,cD,lD,alpha,r,l,value)
  
  use mod_D_X_R_shell_opt_0
  use mod_D_X_R_shell_opt_1
  use mod_D_X_R_shell_opt_2
  use mod_D_X_R_shell_opt_3
  use mod_D_X_R_shell_opt_4
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r(3)     !< center for third spherical Harmonic
  real(kind=8) :: alpha    !< exponent for third spherical Harmonic
  real(kind=8) :: value(*) !< value of the three center integral upon return 
  integer      :: lD       !< angular momentum for first spherical Harmonic
  integer      :: l        !< l momenta of cubic Harmonic shell
  real(kind=8) :: rD(3)
  real(kind=8) :: alphaD
  real(kind=8) :: cD(*)
  
  ! return value
  logical      :: D_Coulomb_C_shell
  
  ! local variables
  real(kind=8)  :: rc(3)
  integer       :: ir
  integer       :: is
  integer       :: nx
  integer       :: ny
  integer       :: nz
  real(kind=8)  :: cr(455)
  real(kind=8)  :: D_Coulomb_D_unpacked
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! center rD
  rc = r - rD
  
  ! call derivative routine
  select case (l)
    case (0)
      D_Coulomb_C_shell = D_X_R_shell_opt_0(alphaD,alpha,rc,cD,lD,value)
    case (1)
      D_Coulomb_C_shell = D_X_R_shell_opt_1(alphaD,alpha,rc,cD,lD,value)
    case (2)
      D_Coulomb_C_shell = D_X_R_shell_opt_2(alphaD,alpha,rc,cD,lD,value)
    case (3)
      D_Coulomb_C_shell = D_X_R_shell_opt_3(alphaD,alpha,rc,cD,lD,value)
    case (4)
      D_Coulomb_C_shell = D_X_R_shell_opt_4(alphaD,alpha,rc,cD,lD,value)
    case default
      D_Coulomb_C_shell = .false.
  end select
  
  ! norm
  value(1:(l+1)*(l+2)/2) =  value(1:(l+1)*(l+2)/2) * sqrt(PI*PI/(alphaD*alpha)) * PI*PI/(alphaD*alpha)
  
  ! fall back on standard procedure if not possible...
  if ( .not.D_Coulomb_C_shell ) then
    
    ! set flag
    D_Coulomb_C_shell = .true.
    
    ! loop on shell components
    is=0
    do nx=l,0,-1
      do ny=l-nx,0,-1
        
        ! get nz exponent
        nz=l-nx-ny
        
        ! increment shell component index
        is=is+1
        
        ! compute contribution from this index
        call R_to_D(alpha,ir_index(nx,ny,nz),cr)
        
        ! compute integral in D rep
        value(is) = D_Coulomb_D_unpacked(alphaD,rD,cD,lD,alpha,r,cr,l)
        
      end do
    end do
    
  end if
  
end function


!> @file Coulomb.f90
!!
!! Computing coulomb integrals between gaussian basis elements, defined as 
!! derivatives of (s||s) kernel.
!!
!! Author: I. Duchemin December 2017
!!
recursive function D_Coulomb_C_shells(alphaD,rD,cD,lD,alpha,rx,ry,rz,l,nshell,value)
  
  use mod_D_X_R_shell_optv1_0
  use mod_D_X_R_shell_optv4_0
  use mod_D_X_R_shell_optv1_1
  use mod_D_X_R_shell_optv4_1
  use mod_D_X_R_shell_optv1_2
  use mod_D_X_R_shell_optv4_2
  use mod_D_X_R_shell_optv1_3
  use mod_D_X_R_shell_optv4_3
  use mod_D_X_R_shell_optv1_4
  use mod_D_X_R_shell_optv4_4
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: rx(*)    !< center for third spherical Harmonics
  real(kind=8) :: ry(*)    !< center for third spherical Harmonics
  real(kind=8) :: rz(*)    !< center for third spherical Harmonics
  real(kind=8) :: alpha(*) !< exponent for third spherical Harmonics
  real(kind=8) :: value(*) !< value of the three center integral upon return 
  integer      :: lD       !< angular momentum for first spherical Harmonic
  integer      :: l        !< l momenta of cubic Harmonic shells
  integer      :: nshell   !< number of cubic Harmonic shells to be treated
  real(kind=8) :: rD(3)
  real(kind=8) :: alphaD
  real(kind=8) :: cD(*)
  
  ! return value
  logical :: D_Coulomb_C_shells
  
  ! local variables
  integer      :: i
  integer      :: s
  real(kind=8) :: rc_v4(4,3)
  real(kind=8) :: rc_v1(1,3)
  logical      :: D_Coulomb_C_shell

  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  
  ! set size of the shell
  s = (l+1)*(l+2)/2
  
  ! call derivative routine
  select case (l)
    case (0)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv1_0(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv4_0(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
    case (1)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv1_1(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv4_1(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
    case (2)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv1_2(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv4_2(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
    case (3)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv1_3(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv4_3(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
    case (4)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv1_4(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_C_shells = D_X_R_shell_optv4_4(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_C_shells ) return
      end do
    case default
      D_Coulomb_C_shells = .false.
  end select
  
  ! norm
  do i = 1,nshell
    value((i-1)*s+1:i*s) =  value((i-1)*s+1:i*s) * sqrt(PI*PI/(alphaD*alpha(i))) * PI*PI/(alphaD*alpha(i))
  end do
  
  ! try falling back to standard procedure if somethng went wrong
  if ( .not.D_Coulomb_C_shells ) then
    ! compute shells one by one
    do i = 1,nshell
      ! set rc
      rc_v1(1,1) = rx(i)
      rc_v1(1,2) = ry(i)
      rc_v1(1,3) = rz(i)
      ! call routine
      D_Coulomb_C_shells = D_Coulomb_C_shell(alphaD,rD,cD,lD,alpha(i),rc_v1,l,value((i-1)*s+1:(i+0)*s))
      ! test that everything works so far
      if ( .not.D_Coulomb_C_shells ) return
    end do
  end if
  
end function


!> @file Coulomb.f90
!!
!! Computing coulomb integrals between gaussian basis elements, defined as 
!! derivatives of (s||s) kernel.
!!
!! Author: I. Duchemin December 2017
!!
recursive function D_Coulomb_Y_shells(alphaD,rD,cD,lD,alpha,rx,ry,rz,l,nshell,value)
  
  use mod_D_X_Y_shell_optv1_0
  use mod_D_X_Y_shell_optv4_0
  use mod_D_X_Y_shell_optv1_1
  use mod_D_X_Y_shell_optv4_1
  use mod_D_X_Y_shell_optv1_2
  use mod_D_X_Y_shell_optv4_2
  use mod_D_X_Y_shell_optv1_3
  use mod_D_X_Y_shell_optv4_3
  use mod_D_X_Y_shell_optv1_4
  use mod_D_X_Y_shell_optv4_4
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: rx(*)    !< center for third spherical Harmonics
  real(kind=8) :: ry(*)    !< center for third spherical Harmonics
  real(kind=8) :: rz(*)    !< center for third spherical Harmonics
  real(kind=8) :: alpha(*) !< exponent for third spherical Harmonics
  real(kind=8) :: value(*) !< value of the three center integral upon return 
  integer      :: lD       !< angular momentum for first spherical Harmonic
  integer      :: l        !< l momenta of cubic Harmonic shells
  integer      :: nshell   !< number of cubic Harmonic shells to be treated
  real(kind=8) :: rD(3)
  real(kind=8) :: alphaD
  real(kind=8) :: cD(*)
  
  ! return value
  logical :: D_Coulomb_Y_shells
  
  ! local variables
  integer      :: i
  integer      :: s
  real(kind=8) :: rc_v4(4,3)
  real(kind=8) :: rc_v1(1,3)
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  
  ! set size of the shell
  s = 2*l+1
  
  ! call derivative routine
  select case (l)
    case (0)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv1_0(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv4_0(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
    case (1)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv1_1(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv4_1(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
    case (2)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv1_2(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv4_2(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
    case (3)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv1_3(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv4_3(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
    case (4)
      ! initial loop
      do i = 1,modulo(nshell,4)
        ! set rc
        rc_v1(1,1) = rx(i)-rD(1)
        rc_v1(1,2) = ry(i)-rD(2)
        rc_v1(1,3) = rz(i)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv1_4(alphaD,alpha(i:i+0),rc_v1,cD,lD,value((i-1)*s+1:(i+0)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
      ! vectorized loop
      do i = modulo(nshell,4)+1,nshell-1,4
        ! set rc
        rc_v4(:,1) = rx(i:i+3)-rD(1)
        rc_v4(:,2) = ry(i:i+3)-rD(2)
        rc_v4(:,3) = rz(i:i+3)-rD(3)
        ! call routine
        D_Coulomb_Y_shells = D_X_Y_shell_optv4_4(alphaD,alpha(i:i+3),rc_v4,cD,lD,value((i-1)*s+1:(i+3)*s))
        if ( .not.D_Coulomb_Y_shells ) return
      end do
    case default
      D_Coulomb_Y_shells = .false.
      return
  end select
  
  ! norm
  do i = 1,nshell
    value((i-1)*s+1:i*s) =  value((i-1)*s+1:i*s) * sqrt(PI*PI/(alphaD*alpha(i))) * PI*PI/(alphaD*alpha(i))
  end do
  
end function







!> Move Cubic Spherical Harmonics to derivative representation
!!
recursive subroutine C_to_D(alpha,r,nx,ny,nz,alphaD,rD,cD,iD,nD)
  
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r(3)     !< center for the cubic Harmonic
  real(kind=8) :: alpha    !< exponent for the cubic Harmonic
  integer      :: nx       !< x power for the cubic Harmonic
  integer      :: ny       !< y power for the cubic Harmonic
  integer      :: nz       !< z power for the cubic Harmonic
  real(kind=8) :: rD(3)    !< center for the derivatives basis
  real(kind=8) :: alphaD   !< exponent for the derivatives basis
  integer      :: iD(*)    !< corresponding derivative index array
  real(kind=8) :: cD(*)    !< corresponding derivative coefficient array
  integer      :: nD       !< number of corresponding derivatives
  
  ! local variable
  integer       :: il
  integer       :: ir1
  real(kind=8)  :: c1(455)
  real(kind=8)  :: coeffs1(455)
  real(kind=8)  :: coeffs_tmp(455)
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! set center and exponant
  rD(:) =r(:)
  alphaD=alpha
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx,ny,nz)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx+ny+nz))=0.0d0
  c1(ir1)=1.0d0
  
  ! form coeffs for derivatives
  coeffs1=0.0d0
  do il=1,imax(nx+ny+nz)
    if ( c1(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha,il,coeffs_tmp)
      ! add contrib
      coeffs1(1:il)=coeffs1(1:il)+c1(il)*coeffs_tmp(1:il)
    end if
  end do
  
  ! compress results in result arrays
  nD=0
  do il=1,imax(nx+ny+nz)
    if ( coeffs1(il).ne.0.0d0 ) then
      nD=nD+1
      iD(nD)=il
      cD(nD)=coeffs1(il)
    end if
  end do
  
end subroutine
  
!> Move Solid Spherical Harmonics to derivative representation
!!
recursive subroutine Y_to_D_unpacked(alpha,l,m,cD)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: alpha    !< exponent for the spherical Harmonic
  integer      :: l        !< angular momentum for the spherical Harmonic
  integer      :: m        !< orbital momentum for the spherical Harmonic
  real(kind=8) :: cD(*)    !< corresponding derivative coefficient array
  
  ! local variable
  integer       :: il
  real(kind=8)  :: c1(455)
  real(kind=8)  :: coeffs_tmp(455)
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l,m,c1)
  
  ! form coeffs for derivatives
  cD(1:imax(l))=0.0d0
  do il=1,imax(l)
    if ( c1(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha,il,coeffs_tmp)
      ! add contrib
      cD(1:il)=cD(1:il)+c1(il)*coeffs_tmp(1:il)
    end if
  end do
  
end subroutine


!> Move Solid Spherical Harmonics to derivative representation
!!
recursive subroutine Y_to_D(alpha,r,l,m,alphaD,rD,cD,iD,nD)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r(3)     !< center for the spherical Harmonic
  real(kind=8) :: alpha    !< exponent for the spherical Harmonic
  integer      :: l        !< angular momentum for the spherical Harmonic
  integer      :: m        !< orbital momentum for the spherical Harmonic
  real(kind=8) :: rD(3)    !< center for the derivatives basis
  real(kind=8) :: alphaD   !< exponent for the derivatives basis
  integer      :: iD(*)    !< corresponding derivative index array
  real(kind=8) :: cD(*)    !< corresponding derivative coefficient array
  integer      :: nD       !< number of corresponding derivatives
  
  ! local variable
  integer       :: il
  integer       :: ir1
  real(kind=8)  :: c1(455)
  real(kind=8)  :: coeffs1(455)
  real(kind=8)  :: coeffs_tmp(455)
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! set center and exponant
  rD(:) =r(:)
  alphaD=alpha
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l,m,c1)
  
  ! form coeffs for derivatives
  coeffs1=0.0d0
  do il=1,imax(l)
    if ( c1(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha,il,coeffs_tmp)
      ! add contrib
      coeffs1(1:il)=coeffs1(1:il)+c1(il)*coeffs_tmp(1:il)
    end if
  end do
  
  ! compress results in result arrays
  nD=0
  do il=1,imax(l)
    if ( coeffs1(il).ne.0.0d0 ) then
      nD=nD+1
      iD(nD)=il
      cD(nD)=coeffs1(il)
    end if
  end do
  
end subroutine
  
  
!> Move Cubic Spherical Harmonics product to derivative representation
!!
recursive subroutine CC_to_D(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alphaD,rD,cD,iD,nD)
  
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)     !< center for the cubic Harmonic
  real(kind=8) :: alpha1    !< exponent for the cubic Harmonic
  integer      :: nx1       !< x power for the cubic Harmonic
  integer      :: ny1       !< y power for the cubic Harmonic
  integer      :: nz1       !< z power for the cubic Harmonic
  real(kind=8) :: r2(3)     !< center for the cubic Harmonic
  real(kind=8) :: alpha2    !< exponent for the cubic Harmonic
  integer      :: nx2       !< x power for the cubic Harmonic
  integer      :: ny2       !< y power for the cubic Harmonic
  integer      :: nz2       !< z power for the cubic Harmonic
  real(kind=8) :: rD(3)     !< center for the derivatives basis
  real(kind=8) :: alphaD    !< exponent for the derivatives basis
  integer      :: iD(*)     !< corresponding derivative index array
  real(kind=8) :: cD(*)     !< corresponding derivative coefficient array
  integer      :: nD        !< number of corresponding derivatives
  
  ! local variable
  integer       :: il
  integer       :: ir1
  integer       :: ir2
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: alpha3
  real(kind=8)  :: r3(3)
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,alpha3,r3,c3,l3)
  
  ! set center and exponant
  rD(:) =r3(:)
  alphaD=alpha3
  
  ! form coeffs for derivatives
  c1=0.0d0
  do il=1,imax(l3)
    if ( c3(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha3,il,c2)
      ! add contrib
      c1(1:il)=c1(1:il)+c3(il)*c2(1:il)
    end if
  end do
  
  ! compress results in result arrays
  nD=0
  do il=1,imax(l3)
    if ( c1(il).ne.0.0d0 ) then
      nD=nD+1
      iD(nD)=il
      cD(nD)=c1(il)
    end if
  end do
  
end subroutine
  

!> Move Cubic times Spherical  Spherical Harmonics product to derivative representation
!!
recursive subroutine YC_to_D(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2,alphaD,rD,cD,iD,nD)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)     !< center for the spherical Harmonic
  real(kind=8) :: alpha1    !< exponent for the spherical Harmonic
  integer      :: l1        !< angular momentum for the spherical Harmonic
  integer      :: m1        !< orbital momentum for the spherical Harmonic
  real(kind=8) :: r2(3)     !< center for the cubic Harmonic
  real(kind=8) :: alpha2    !< exponent for the cubic Harmonic
  integer      :: nx2       !< x power for the cubic Harmonic
  integer      :: ny2       !< y power for the cubic Harmonic
  integer      :: nz2       !< z power for the cubic Harmonic
  real(kind=8) :: rD(3)     !< center for the derivatives basis
  real(kind=8) :: alphaD    !< exponent for the derivatives basis
  integer      :: iD(*)     !< corresponding derivative index array
  real(kind=8) :: cD(*)     !< corresponding derivative coefficient array
  integer      :: nD        !< number of corresponding derivatives
  
  ! local variable
  integer       :: il
  integer       :: ir2
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: alpha3
  real(kind=8)  :: r3(3)
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,nx2+ny2+nz2,alpha3,r3,c3,l3)
  
  ! set center and exponant
  rD(:) =r3(:)
  alphaD=alpha3
  
  ! form coeffs for derivatives
  c1=0.0d0
  do il=1,imax(l3)
    if ( c3(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha3,il,c2)
      ! add contrib
      c1(1:il)=c1(1:il)+c3(il)*c2(1:il)
    end if
  end do
  
  ! compress results in result arrays
  nD=0
  do il=1,imax(l3)
    if ( c1(il).ne.0.0d0 ) then
      nD=nD+1
      iD(nD)=il
      cD(nD)=c1(il)
    end if
  end do
  
end subroutine
  

!> Move Cubic Spherical Harmonics product to derivative representation
!!
recursive subroutine CY_to_D(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2,alphaD,rD,cD,iD,nD)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)     !< center for the cubic Harmonic
  real(kind=8) :: alpha1    !< exponent for the cubic Harmonic
  integer      :: nx1       !< x power for the cubic Harmonic
  integer      :: ny1       !< y power for the cubic Harmonic
  integer      :: nz1       !< z power for the cubic Harmonic
  real(kind=8) :: r2(3)     !< center for the spherical Harmonic
  real(kind=8) :: alpha2    !< exponent for the spherical Harmonic
  integer      :: l2        !< angular momentum for the spherical Harmonic
  integer      :: m2        !< orbital momentum for the spherical Harmonic
  real(kind=8) :: rD(3)     !< center for the derivatives basis
  real(kind=8) :: alphaD    !< exponent for the derivatives basis
  integer      :: iD(*)     !< corresponding derivative index array
  real(kind=8) :: cD(*)     !< corresponding derivative coefficient array
  integer      :: nD        !< number of corresponding derivatives
  
  ! local variable
  integer       :: il
  integer       :: ir1
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: alpha3
  real(kind=8)  :: r3(3)
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,l2,alpha3,r3,c3,l3)
  
  ! set center and exponant
  rD(:) =r3(:)
  alphaD=alpha3
  
  ! form coeffs for derivatives
  c1=0.0d0
  do il=1,imax(l3)
    if ( c3(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha3,il,c2)
      ! add contrib
      c1(1:il)=c1(1:il)+c3(il)*c2(1:il)
    end if
  end do
  
  ! compress results in result arrays
  nD=0
  do il=1,imax(l3)
    if ( c1(il).ne.0.0d0 ) then
      nD=nD+1
      iD(nD)=il
      cD(nD)=c1(il)
    end if
  end do
  
end subroutine
  
  
!> Move Spherical times Spherical  Spherical Harmonics product to derivative representation
!!
recursive subroutine YY_to_D(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alphaD,rD,cD,iD,nD)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_R_to_D_Conversion
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)     !< center for the spherical Harmonic
  real(kind=8) :: alpha1    !< exponent for the spherical Harmonic
  integer      :: l1        !< angular momentum for the spherical Harmonic
  integer      :: m1        !< orbital momentum for the spherical Harmonic
  real(kind=8) :: r2(3)     !< center for the cubic Harmonic
  real(kind=8) :: alpha2    !< exponent for the cubic Harmonic
  integer      :: l2        !< angular momentum for the spherical Harmonic
  integer      :: m2        !< orbital momentum for the spherical Harmonic
  real(kind=8) :: rD(3)     !< center for the derivatives basis
  real(kind=8) :: alphaD    !< exponent for the derivatives basis
  integer      :: iD(*)     !< corresponding derivative index array
  real(kind=8) :: cD(*)     !< corresponding derivative coefficient array
  integer      :: nD        !< number of corresponding derivatives
  
  ! local variable
  integer       :: il
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: alpha3
  real(kind=8)  :: r3(3)
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,alpha3,r3,c3,l3)
  
  ! set center and exponant
  rD(:) =r3(:)
  alphaD=alpha3
  
  ! form coeffs for derivatives
  c1=0.0d0
  do il=1,imax(l3)
    if ( c3(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha3,il,c2)
      ! add contrib
      c1(1:il)=c1(1:il)+c3(il)*c2(1:il)
    end if
  end do
  
  ! compress results in result arrays
  nD=0
  do il=1,imax(l3)
    if ( c1(il).ne.0.0d0 ) then
      nD=nD+1
      iD(nD)=il
      cD(nD)=c1(il)
    end if
  end do
  
end subroutine
  
!> Two centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) \frac{1}{|r-r'|} Y_{xyz}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function C_Coulomb_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2)
  
  use mod_XZY_power_to_ir
  
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
  real(kind=8)  :: C_Coulomb_C
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! compute coulomb integral
  C_Coulomb_C =R_X_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2)
  
end function

!> Two centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) \frac{erf(omega|r-r'|)}{|r-r'|} Y_{xyz}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function attenuated_C_Coulomb_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,omega)
  
  use mod_XZY_power_to_ir
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  real(kind=8) :: omega    !< coefficient for attenuation
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: attenuated_C_Coulomb_C
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: attenuated_R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! compute coulomb integral
  attenuated_C_Coulomb_C =attenuated_R_X_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,omega)
  
end function
  
  
!> Two centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) \frac{1}{|r-r'|} Y_{lm}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function C_Coulomb_Y(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  
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
  real(kind=8)  :: C_Coulomb_Y
  
  ! local variable
  integer       :: ir1
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! compute coulomb integral
  C_Coulomb_Y =R_X_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,l2)
  
end function
  
!> Two centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) \frac{erf(omega|r-r'|)}{|r-r'|} Y_{lm}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function attenuated_C_Coulomb_Y(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2,omega)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  real(kind=8) :: omega    !< coefficient for attenuation
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: attenuated_C_Coulomb_Y
  
  ! local variable
  integer       :: ir1
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: attenuated_R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! compute coulomb integral
  attenuated_C_Coulomb_Y =attenuated_R_X_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,l2,omega)
  
end function

!> Two centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) \frac{1}{|r-r'|} Y_{xyz}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function Y_Coulomb_C(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  
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
  real(kind=8)  :: Y_Coulomb_C
  
  ! local variable
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! compute coulomb integral
  Y_Coulomb_C =R_X_R(alpha1,r1,c1,l1,alpha2,r2,c2,nx2+ny2+nz2)
  
end function
  

!> Two centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) \frac{erf(omega|r-r'|)}{|r-r'|} Y_{xyz}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function attenuated_Y_Coulomb_C(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2,omega)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  real(kind=8) :: omega    !< coefficient for attenuation
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: attenuated_Y_Coulomb_C
  
  ! local variable
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: attenuated_R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! compute coulomb integral
  attenuated_Y_Coulomb_C =attenuated_R_X_R(alpha1,r1,c1,l1,alpha2,r2,c2,nx2+ny2+nz2,omega)
  
end function

!> Two centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) \frac{1}{|r-r'|} Y_{lm}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function Y_Coulomb_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2)
  
  use mod_R_from_Y
  
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
  real(kind=8)  :: Y_Coulomb_Y
  
  ! local variable
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! compute coulomb integral
  Y_Coulomb_Y =R_X_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2)
  
end function
  

!> Two centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) \frac{erf(omega|r-r'|)}{|r-r'|} Y_{lm}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function attenuated_Y_Coulomb_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2,omega)
  
  use mod_R_from_Y
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: omega    !< coefficient for attenuation
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: attenuated_Y_Coulomb_Y
  
  ! local variable
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: attenuated_R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! compute coulomb integral
  attenuated_Y_Coulomb_Y =attenuated_R_X_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,omega)
  
end function

!> Two centers ionic integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-R_{ion}|} 
!!  \f$
!!
recursive function CC_Coulomb_Ion(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,rion)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: rion(3)  !< ion position
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: CC_Coulomb_Ion
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha3
  real(kind=8)  :: r3(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  real(kind=8)       :: alpha_ion=1.0e16
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,alpha3,r3,c3,l3)
  
  ! set ion coeff in the R basis
  c4(1)=1.0d0/R_1_norm(alpha_ion,(/1.0d0/),0)
  
  ! compute coulomb integral
  CC_Coulomb_Ion =R_X_R(alpha3,r3,c3,l3,alpha_ion,rion,c4,0)
  
end function

!> Two centers ionic integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-R_{ion}|} 
!!  \f$
!!
recursive function YY_Coulomb_Ion(alpha1,r1,l1,m1,alpha2,r2,l2,m2,rion)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: rion(3)  !< ion position
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: YY_Coulomb_Ion
  
  ! local variable
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha3
  real(kind=8)  :: r3(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  real(kind=8)  :: alpha_ion=1.0e16
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,alpha3,r3,c3,l3)
  
  ! set ion coeff in the R basis
  c4(1)=1.0d0/R_1_norm(alpha_ion,(/1.0d0/),0)
  
  ! compute coulomb integral
  YY_Coulomb_Ion =R_X_R(alpha3,r3,c3,l3,alpha_ion,rion,c4,0)
  
end function


!> Three centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{xyz}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function CC_Coulomb_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,nx3,ny3,nz3)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
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
  real(kind=8)  :: CC_Coulomb_C
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  integer       :: ir3
  integer       :: l4
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! get 3rd member index in the R basis
  ir3 =ir_index(nx3,ny3,nz3)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! set 3rd member coeff in the R basis
  c3(1:imax(nx3+ny3+nz3))=0.0d0
  c3(ir3)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,alpha4,r4,c4,l4)
  
  ! compute coulomb integral
  CC_Coulomb_C =R_X_R(alpha3,r3,c3,nx3+ny3+nz3,alpha4,r4,c4,l4)
  
end function


!> Three centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{lm}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function CC_Coulomb_Y(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,l3,m3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  
  ! return value
  real(kind=8)  :: CC_Coulomb_Y
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  integer       :: l4
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,alpha4,r4,c4,l4)
  
  ! compute coulomb integral
  CC_Coulomb_Y =R_X_R(alpha3,r3,c3,l3,alpha4,r4,c4,l4)
  
end function

!> Three centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{xyz}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function CY_Coulomb_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2,alpha3,r3,nx3,ny3,nz3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  
  ! return value
  real(kind=8)  :: CY_Coulomb_C
  
  ! local variable
  integer       :: ir1
  integer       :: ir3
  integer       :: l4
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 3rd member index in the R basis
  ir3 =ir_index(nx3,ny3,nz3)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! set 3rd member coeff in the R basis
  c3(1:imax(nx3+ny3+nz3))=0.0d0
  c3(ir3)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,l2,alpha4,r4,c4,l4)
  
  ! compute coulomb integral
  CY_Coulomb_C =R_X_R(alpha3,r3,c3,nx3+ny3+nz3,alpha4,r4,c4,l4)
  
end function


!> Three centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{xyz}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function YC_Coulomb_C(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,nx3,ny3,nz3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  
  ! return value
  real(kind=8)  :: YC_Coulomb_C
  
  ! local variable
  integer       :: ir2
  integer       :: ir3
  integer       :: l4
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! get 3rd member index in the R basis
  ir3 =ir_index(nx3,ny3,nz3)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! set 3rd member coeff in the R basis
  c3(1:imax(nx3+ny3+nz3))=0.0d0
  c3(ir3)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,nx2+ny2+nz2,alpha4,r4,c4,l4)
  
  ! compute coulomb integral
  YC_Coulomb_C =R_X_R(alpha3,r3,c3,nx3+ny3+nz3,alpha4,r4,c4,l4)
  
end function

!> Three centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{xyz}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function YY_Coulomb_C(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,nx3,ny3,nz3) result(value)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  use mod_R_to_D_Conversion
  
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
  real(kind=8)  :: value
  
  ! local variable
  integer       :: il
  integer       :: l3
  integer       :: l4
  integer       :: ir3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  logical       :: D_Coulomb_C
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,alpha4,r4,c4,l4)
  
  ! get l3
  l3 = nx3+ny3+nz3
  
  ! if we may use specialized routine
  if ( ( l4.le.6 ).and.( l3.le.4 ) ) then
    
    ! form coeffs in the derivative basis
    c1=0.0d0
    do il=1,imax(l4)
      if ( c4(il).ne.0.0d0 ) then
        ! compute contribution from this index
        call R_to_D(alpha4,il,c2)
        ! add contrib
        c1(1:il)=c1(1:il)+c4(il)*c2(1:il)
      end if
    end do
    
    ! if we can get the value through specialized routine
    if ( D_Coulomb_C(alpha4,r4,c1,l4,alpha3,r3,nx3,ny3,nz3,value) ) then
      return
    end if
    
  end if
  
  ! otherwise, get 3rd member index in the R basis
  ir3 =ir_index(nx3,ny3,nz3)
  
  ! set 3rd member coeff in the R basis
  c3(1:imax(nx3+ny3+nz3))=0.0d0
  c3(ir3)=1.0d0
  
  ! compute coulomb integral with general routine
  value =R_X_R(alpha3,r3,c3,nx3+ny3+nz3,alpha4,r4,c4,l4)
  
end function


!> Three centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{lm}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function YC_Coulomb_Y(alpha1,r1,l1,m1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,l3,m3)
  
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
  real(kind=8)  :: YC_Coulomb_Y
  
  ! local variable
  integer       :: l4
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,nx2+ny2+nz2,alpha4,r4,c4,l4)
  
  ! compute coulomb integral
  YC_Coulomb_Y =R_X_R(alpha3,r3,c3,l3,alpha4,r4,c4,l4)
  
end function


!> Three centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{lm}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function CY_Coulomb_Y(alpha1,r1,nx1,ny1,nz1,alpha2,r2,l2,m2,alpha3,r3,l3,m3)
  
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
  real(kind=8)  :: CY_Coulomb_Y
  
  ! local variable
  integer       :: l4
  integer       :: ir1
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,l2,alpha4,r4,c4,l4)
  
  ! compute coulomb integral
  CY_Coulomb_Y =R_X_R(alpha3,r3,c3,l3,alpha4,r4,c4,l4)
  
end function


!> Three centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{lm}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function YY_Coulomb_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,l3,m3) result(value)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_CubicHarmonicsProduct
  use mod_R_to_D_Conversion
  
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
  logical      :: D_Coulomb_Y
  
  ! return value
  real(kind=8)  :: value
  
  ! local variable
  integer       :: il
  integer       :: l4
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,alpha4,r4,c4,l4)
  
  ! if we may use specialized routine
  if ( ( l4.le.6 ).and.( l3.le.4 ) ) then
    
    ! form coeffs in the derivative basis
    c1=0.0d0
    do il=1,imax(l4)
      if ( c4(il).ne.0.0d0 ) then
        ! compute contribution from this index
        call R_to_D(alpha4,il,c2)
        ! add contrib
        c1(1:il)=c1(1:il)+c4(il)*c2(1:il)
      end if
    end do
    
    ! if we can get the value through specialized routine
    if ( D_Coulomb_Y(alpha4,r4,c1,l4,alpha3,r3,l3,m3,value) ) then
      return
    end if
    
  end if
  
  ! otherwise, get decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! compute coulomb integral trhough general routine
  value=R_X_R(alpha3,r3,c3,l3,alpha4,r4,c4,l4)
  
end function


!> Four centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{xyz}^{(3)}(r'-R_3) Y_{xyz}^{(4)}(r'-R_4)
!!  \f$
!!
recursive function CC_Coulomb_CC(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,nx3,ny3,nz3,alpha4,r4,nx4,ny4,nz4)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: r4(3)    !< center for fourth cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  real(kind=8) :: alpha4   !< exponent for fourth spherical Harmonic
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
  real(kind=8)  :: CC_Coulomb_CC
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  integer       :: ir3
  integer       :: ir4
  integer       :: l5
  integer       :: l6
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: c5(455)
  real(kind=8)  :: c6(455)
  real(kind=8)  :: alpha5
  real(kind=8)  :: alpha6
  real(kind=8)  :: r5(3)
  real(kind=8)  :: r6(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! get 3rd member index in the R basis
  ir3 =ir_index(nx3,ny3,nz3)
  
  ! get 4th member index in the R basis
  ir4 =ir_index(nx4,ny4,nz4)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! set 3rd member coeff in the R basis
  c3(1:imax(nx3+ny3+nz3))=0.0d0
  c3(ir3)=1.0d0
  
  ! set 4th member coeff in the R basis
  c4(1:imax(nx4+ny4+nz4))=0.0d0
  c4(ir4)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,alpha5,r5,c5,l5)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha3,r3,c3,nx3+ny3+nz3,alpha4,r4,c4,nx4+ny4+nz4,alpha6,r6,c6,l6)
  
  ! compute coulomb integral
  CC_Coulomb_CC =R_X_R(alpha5,r5,c5,l5,alpha6,r6,c6,l6)
  
end function

!> Four centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{lm}^{(3)}(r'-R_3) Y_{lm}^{(4)}(r'-R_4)
!!  \f$
!!
recursive function YY_Coulomb_YY(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,l3,m3,alpha4,r4,l4,m4)
  
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
  real(kind=8)  :: YY_Coulomb_YY
  
  ! local variable
  integer       :: l5
  integer       :: l6
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: c5(455)
  real(kind=8)  :: c6(455)
  real(kind=8)  :: alpha5
  real(kind=8)  :: alpha6
  real(kind=8)  :: r5(3)
  real(kind=8)  :: r6(3)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l4,m4,c4)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,alpha5,r5,c5,l5)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha3,r3,c3,l3,alpha4,r4,c4,l4,alpha6,r6,c6,l6)
  
  ! compute coulomb integral
  YY_Coulomb_YY =R_X_R(alpha5,r5,c5,l5,alpha6,r6,c6,l6)
  
end function


!> Two centers Modified Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!         /
!! compute | dr1 dr2 Ylm1(r1-R1) * Y00(r1-R2) * 1/|r1-r2| * Ylm2(r2-R2)
!!         /
!!
recursive function Y_ModCoulomb_Y(acut,alpha1,r1,l1,m1,alpha2,r2,l2,m2)
  
  use mod_CubicHarmonicsProduct
  use mod_R_from_Y
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: acut     !< exponent for cutoff spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: Y_ModCoulomb_Y
  
  ! local variable
  real(kind=8)  :: c1(455),c2(455),c3(455)
  real(kind=8)  :: ac
  integer       :: lmax
  real(kind=8)  :: rc(3)
  real(kind=8)  :: cc(455)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! decomposition of the cutoff solid harmonic into cubic harmonic
  c3   =0.0d0
  c3(1)=1.0d0
  
  ! get product of first solid harmonic times cutoff solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,acut,r2,c3,0,ac,rc,cc,lmax)
  
  ! compute coulomb integral
  Y_ModCoulomb_Y =R_X_R(ac,rc,cc,lmax,alpha2,r2,c2,l2)
  
end function


!> Three centers Electron Modified Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!         /
!! compute | dr1 dr2 Ylm1(r1-R1) * Ylm2(r1-R2) * Y00(r1-R3) * 1/|r1-r2| * Ylm3(r2-R3)
!!         /
!!
recursive function YY_ModCoulomb_Y(acut,alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,l3,m3)
  
  use mod_CubicHarmonicsProduct
  use mod_R_from_Y
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: acut     !< exponent for cutoff spherical Harmonic
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
  real(kind=8)  :: YY_ModCoulomb_Y
  
  ! local variable
  real(kind=8)  :: c1(455),c2(455),c3(455)
  real(kind=8)  :: ac,atmp
  integer       :: lmax,ltmp
  real(kind=8)  :: rc(3),rtmp(3)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,atmp,rtmp,c3,ltmp)
  
  ! get product with cutoff solid harmonic
  c2(1)=1.0d0
  call RxR_to_R(atmp,rtmp,c3,ltmp,acut,r3,c2, 0,ac,rc,c1,lmax)
  
  ! decomposition of the third solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! compute coulomb integral
  YY_ModCoulomb_Y =R_X_R(ac,rc,c1,lmax,alpha3,r3,c3,l3)
  
end function

