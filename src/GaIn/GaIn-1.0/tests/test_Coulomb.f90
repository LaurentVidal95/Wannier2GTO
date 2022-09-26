program test_Coulomb
  
  use mod_R_1_norm
  
  implicit none
  
  ! local variables
  integer, parameter      :: nmax = 3
  integer, parameter      :: ng   = 127
  real(kind=8), parameter :: a    = 1.0d16
  real(kind=8), parameter :: b    = 3.5d0
  real(kind=8), parameter :: d    = 2.0d0
  real(kind=8), parameter :: e    = 1.5d0
  real(kind=8), parameter :: lg   = 7.0d0
  real(kind=8), parameter :: x0   = -12.0d0/2
  real(kind=8), parameter :: y0   = -13.5d0/2
  real(kind=8), parameter :: z0   =  15.0d0/2
  integer, parameter      :: nx2  = 2
  integer, parameter      :: ny2  = 0
  integer, parameter      :: nz2  = 1
  integer, parameter      :: nx3  = 0
  integer, parameter      :: ny3  = 1
  integer, parameter      :: nz3  = 1
  integer, parameter      :: l1   = 1
  integer, parameter      :: m1   = 0
  integer, parameter      :: l2   = 2
  integer, parameter      :: m2   =-1
  integer, parameter      :: l3   = 3
  integer, parameter      :: m3   = 2
  integer, parameter      :: lmax = 6
  real(kind=8), parameter :: r1(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  real(kind=8), parameter :: r2(3) = (/ 0.5d0,-0.3d0, 0.2d0 /)
  real(kind=8), parameter :: r3(3) = (/ 0.1d0, 0.0d0, 0.5d0 /)
  integer      :: i
  integer      :: j
  integer      :: k
  integer      :: l
  integer      :: m
  integer      :: nx
  integer      :: ny
  integer      :: nz
  integer      :: nx_
  integer      :: ny_
  integer      :: nz_
  real(kind=8) :: c(455)
  real(kind=8) :: x
  real(kind=8) :: y
  real(kind=8) :: z
  real(kind=8) :: tmp
  real(kind=8) :: V
  real(kind=8) :: S1
  real(kind=8) :: S2
  real(kind=8) :: S3
  real(kind=8) :: S1_
  real(kind=8) :: S1__
  real(kind=8) :: C_Coulomb_C
  real(kind=8) :: attenuated_C_Coulomb_C
  real(kind=8) :: Y_Coulomb_Y
  real(kind=8) :: attenuated_Y_Coulomb_Y
  real(kind=8) :: Y_Coulomb_C
  real(kind=8) :: C_Coulomb_Y
  real(kind=8) :: CC_Coulomb_Ion
  real(kind=8) :: CC_Coulomb_C
  real(kind=8) :: CC_Coulomb_Y
  real(kind=8) :: CY_Coulomb_C
  real(kind=8) :: YC_Coulomb_C
  real(kind=8) :: CC_Coulomb_CC
  real(kind=8) :: YY_Coulomb_Ion
  real(kind=8) :: YY_Coulomb_Y
  real(kind=8) :: YY_Coulomb_C
  real(kind=8) :: YC_Coulomb_Y
  real(kind=8) :: CY_Coulomb_Y
  real(kind=8) :: YY_Coulomb_YY
  real(kind=8) :: Y_Value
  real(kind=8) :: D_Coulomb_D
  real(kind=8) :: norme
  real(kind=8) :: omega
  real(kind=8) :: dist
  real(kind=8) :: alphaD
  real(kind=8) :: cD(455)
  real(kind=8) :: iD(455)
  real(kind=8) :: rD(3)
  integer      :: nd
  real(kind=8) :: alphaD2
  real(kind=8) :: cD2(455)
  real(kind=8) :: iD2(455)
  real(kind=8) :: rD2(3)
  integer      :: nd2
  
  INTEGER :: count_rate
  INTEGER :: start_count,end_count
            
  ! compute 1 norm for delta function
  c=0.0d0
  c(1)=1.0d0
  norme=R_1_norm(a,c,0)
  
!  print *,'tic'
!  CALL SYSTEM_CLOCK(start_count, count_rate)
!  do i=1,1000000
!!    call  Y_to_D(b,r2,l2,m2,alphaD,rD,cD,iD,nD) 
!    call  Y_to_D(d,r3,l3,m3,alphaD2,rD2,cD2,iD2,nD2)
!    call  YY_to_D(b,r2,l2,m2,d,r3,l3,m3,alphaD,rD,cD,iD,nD)
!  end do
!  CALL SYSTEM_CLOCK(end_count, count_rate)
!  print *,'tac',(1.0d0*(end_count-start_count))/count_rate
!  CALL SYSTEM_CLOCK(start_count, count_rate)
!  do i=1,1000000
!!    S1 = D_Coulomb_D(alphaD,rD,cD,iD,nD,alphaD2,rD2,cD2,iD2,nD2)
!    S1 = D_Coulomb_D(alphaD,rD,cD,iD,nD,alphaD2,rD2,cD2,iD2,nD2)
!!    S1 = D_Coulomb_D(alphaD,rD,cD,iD,nD,alphaD,rD,cD,iD,nD)
!  end do
!  CALL SYSTEM_CLOCK(end_count, count_rate)
!  print *,'tac',(1.0d0*(end_count-start_count))/count_rate
!  CALL SYSTEM_CLOCK(start_count, count_rate)
!  do i=1,1000000
!!    S1 = Y_Coulomb_Y(b,r2,l2,m2,d,r3,l3,m3)
!    S1 = YY_Coulomb_Y(b,r2,l2,m2,d,r3,l3,m3,d,r3,l3,m3)
!!    S1 = YY_Coulomb_YY(b,r2,l2,m2,d,r3,l3,m3,b,r2,l2,m2,d,r3,l3,m3)
!  end do
!  CALL SYSTEM_CLOCK(end_count, count_rate)
!  print *,'toc',(1.0d0*(end_count-start_count))/count_rate
!  
!  stop 1
  
  ! check consistency between expert and standard drivers
  do nx=0,6
  do ny=0,6-nx
  do nz=0,6-nx-ny
    
    ! compute D basis decomposition
    call  C_to_D(b,r2,nx,ny,nz,alphaD,rD,cD,iD,nD)
    call  C_to_D(d,r3,nx2,ny2,nz2,alphaD2,rD2,cD2,iD2,nD2)
    
    ! check integral 
    S1 = D_Coulomb_D(alphaD,rD,cD,iD,nD,alphaD2,rD2,cD2,iD2,nD2)
    S2 = C_Coulomb_C(b,r2,nx,ny,nz,d,r3,nx2,ny2,nz2)
    
    write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  expert driver=',S1,'  C_Coulomb_C result=',S2
        
    if ( 2.0* abs( ( S1 - S2 ) / ( S1 + S2 ) ) > 1.0d-10 .or. S1.ne.S1 ) then
      print *,'Error: expert <CC> Two center Coulomb failed'
      print *,''
      print *,'test_Coulomb failed'
      print *,''
      stop 1
    end if
    
  end do
  end do
  end do
  
  print *,''
  write (*,'(A)') 'test_Coulomb passed up to l=6'
  print *,''
  
  ! check consistency between expert and standard drivers
  do l=0,6
  do m=-l,l
    
    ! compute D basis decomposition
    call  Y_to_D(b,r2,l,m,alphaD,rD,cD,iD,nD)
    call  C_to_D(d,r3,nx2,ny2,nz2,alphaD2,rD2,cD2,iD2,nD2)
    
    ! check integral 
    S1 = D_Coulomb_D(alphaD,rD,cD,iD,nD,alphaD2,rD2,cD2,iD2,nD2)
    S2 = Y_Coulomb_C(b,r2,l,m,d,r3,nx2,ny2,nz2)
    
    write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'l=',l,'  m=',m,'  expert driver=',S1,'   Y_Coulomb_C result=',S2
        
    if ( 2.0* abs( ( S1 - S2 ) / ( S1 + S2 ) ) > 1.0d-10 .or. S1.ne.S1 ) then
      print *,'Error: expert <CC> Two center Coulomb failed'
      print *,''
      print *,'test_Coulomb failed'
      print *,''
      stop 1
    end if
    
    call  YY_to_D(b,r2,l,m,d,r3,l2,m2,alphaD,rD,cD,iD,nD)
    S1 = D_Coulomb_D(alphaD,rD,cD,iD,nD,alphaD2,rD2,cD2,iD2,nD2)
    S2 = YY_Coulomb_C(b,r2,l,m,d,r3,l2,m2,d,r3,nx2,ny2,nz2)
    
    write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'l=',l,'  m=',m,'  expert driver=',S1,'  YY_Coulomb_C result=',S2
        
    if ( 2.0* abs( ( S1 - S2 ) / ( S1 + S2 ) ) > 1.0d-10 .or. S1.ne.S1 ) then
      print *,'Error: expert <CC> Two center Coulomb failed'
      print *,''
      print *,'test_Coulomb failed'
      print *,''
      stop 1
    end if
    
  end do
  end do
  
  print *,''
  write (*,'(A)') 'test_Coulomb passed up to l=6'
  print *,''
  
  ! check consistency between expert and standard drivers
  do nx=0,6
  do ny=0,6-nx
  do nz=0,6-nx-ny
    
    ! compute D basis decomposition
    call  C_to_D(b,r2,nx,ny,nz,alphaD2,rD2,cD2,iD2,nD2)
    call  CC_to_D(b,r2,nx,ny,nz,d,r3,nx2,ny2,nz2,alphaD,rD,cD,iD,nD)
    
    ! check integral 
    S1 = D_Coulomb_D(alphaD,rD,cD,iD,nD,alphaD,rD,cD,iD,nD)
    S2 = CC_Coulomb_CC(b,r2,nx,ny,nz,d,r3,nx2,ny2,nz2,b,r2,nx,ny,nz,d,r3,nx2,ny2,nz2)
    
    write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  expert driver=',S1,'  CC_Coulomb_CC result=',S2
      
    if ( 2.0* abs( ( S1 - S2 ) / ( S1 + S2 ) ) > 1.0d-10 .or. S1.ne.S1 ) then
      print *,'Error: expert <CC> Four center Coulomb failed'
      print *,''
      print *,'test_Coulomb failed'
      print *,''
      stop 1
    end if
    
  end do
  end do
  end do
  
  print *,''
  write (*,'(A)') 'test_Coulomb passed up to l=6'
  print *,''
  
  ! check consistency between expert and standard drivers
  do nx=0,6
  do ny=0,6-nx
  do nz=0,6-nx-ny
    
    ! compute D basis decomposition
    call  C_to_D(b,r2,nx,ny,nz,alphaD2,rD2,cD2,iD2,nD2)
    call  CY_to_D(b,r2,nx,ny,nz,d,r3,l2,m2,alphaD,rD,cD,iD,nD)
    S1 = D_Coulomb_D(alphaD,rD,cD,iD,nD,alphaD2,rD2,cD2,iD2,nD2)
    call  YC_to_D(d,r3,l2,m2,b,r2,nx,ny,nz,alphaD,rD,cD,iD,nD)
    S3 = D_Coulomb_D(alphaD2,rD2,cD2,iD2,nD2,alphaD,rD,cD,iD,nD)
    
    ! check integral 
    S2 = CY_Coulomb_C(b,r2,nx,ny,nz,d,r3,l2,m2,b,r2,nx,ny,nz)
    
    write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  expert driver=',S1,'  CY_Coulomb_C result=',S2
    
    if ( 2.0* abs( ( S1 - S2 ) / ( S1 + S2 ) ) > 1.0d-10 .or. S1.ne.S1 ) then
      print *,'Error: expert <CC> Four center Coulomb failed'
      print *,''
      print *,'test_Coulomb failed'
      print *,''
      stop 1
    end if
    
  end do
  end do
  end do
  
  print *,''
  write (*,'(A)') 'test_Coulomb passed up to l=6'
  print *,''
  
  !!
  !! check for C routines
  !!
  norme=R_1_norm(a,c,0)
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! explicit calculation
        S1=0.0d0
        
        do k=0,ng
          z=-lg+(2.0d0*lg*k)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do i=0,ng
              x=-lg+(2.0d0*lg*i)/ng
              V  = 1.0d0/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
              tmp= x**nx * y**ny * z**nz * exp( -b * (x**2+y**2+z**2) )
              S1 = S1 + tmp * V
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        
        ! get corresponding overlap
        S1_= C_Coulomb_C(a,(/x0,y0,z0/),0,0,0,b,r1,nx,ny,nz)/norme
        
        ! get corresponding symmetric overlap
        S1__= C_Coulomb_C(b,r1,nx,ny,nz,a,(/x0,y0,z0/),0,0,0)/norme
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S1,'  C_Coulomb_C result=',S1_,S1__
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-5 ) then
          print *,'Error: Coulomb differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1__)/S1) > 1.0d-5 ) then
          print *,'Error: Coulomb differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Coulomb passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''
  
  !!
  !! check for Y routine
  !!
  norme = 2.0d0/(4.0d0*atan(1.0d0))*(a*sqrt(a))
  do l=0,lmax
    do m=-l,l
      ! explicit calculation
      S1=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            V  = 1.0d0/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
            tmp=Y_Value(b,l,m,(/x,y,z/))
            S1 = S1 + tmp * V
          end do
        end do
      end do
      S1=S1*8.0d0*lg**3/ng**3
      
      ! get corresponding interaction
      S1_= Y_Coulomb_Y(a,(/x0,y0,z0/),0,0,b,r1,l,m)*norme
      
      ! get corresponding symmetric interaction
      S1__= Y_Coulomb_Y(b,r1,l,m,a,(/x0,y0,z0/),0,0)*norme
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S1,'  Y_Coulomb_Y result=',S1_,S1__
        
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-5 ) then
        print *,'Error: Coulomb differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1__)/S1) > 1.0d-5 ) then
        print *,'Error: Coulomb differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Coulomb passed up to lmax=',lmax
  print *,''
  
  !!
  !! check for C routines
  !!
  norme=R_1_norm(a,c,0)
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! explicit calculation
        S1=0.0d0
        omega=0.1D0
        do k=0,ng
          z=-lg+(2.0d0*lg*k)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do i=0,ng
              x=-lg+(2.0d0*lg*i)/ng
              dist = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
              V    = erf(omega*dist)/dist
              tmp  = x**nx * y**ny * z**nz * exp( -b * (x**2+y**2+z**2) )
              S1   = S1 + tmp * V
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        
        ! get corresponding overlap
        S1_= attenuated_C_Coulomb_C(a,(/x0,y0,z0/),0,0,0,b,r1,nx,ny,nz,omega)/norme
        
        ! get corresponding symmetric overlap
        S1__= attenuated_C_Coulomb_C(b,r1,nx,ny,nz,a,(/x0,y0,z0/),0,0,0,omega)/norme
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S1,'  attenuated_C_Coulomb_C result=',S1_,S1__
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-4 ) then
          print *,'Error: Coulomb differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1__)/S1) > 1.0d-4 ) then
          print *,'Error: Coulomb differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Coulomb passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''

  
  !!
  !! check for Y routine
  !!
  norme = 2.0d0/(4.0d0*atan(1.0d0))*(a*sqrt(a))
  do l=0,lmax
    do m=-l,l
      ! explicit calculation
      omega=0.1D0
      S1=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            dist = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
            V    = erf(omega*dist)/dist
            tmp=Y_Value(b,l,m,(/x,y,z/))
            S1 = S1 + tmp * V
          end do
        end do
      end do
      S1=S1*8.0d0*lg**3/ng**3
      
      ! get corresponding interaction
      S1_= attenuated_Y_Coulomb_Y(a,(/x0,y0,z0/),0,0,b,r1,l,m,omega)*norme
      
      ! get corresponding symmetric interaction
      S1__= attenuated_Y_Coulomb_Y(b,r1,l,m,a,(/x0,y0,z0/),0,0,omega)*norme
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S1,'  attenuated_Y_Coulomb_Y result=',S1_,S1__
        
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-4 ) then
        print *,'Error: Coulomb differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1__)/S1) > 1.0d-4 ) then
        print *,'Error: Coulomb differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Coulomb passed up to lmax=',lmax
  print *,''
  
  
  !!
  !! check for C routines
  !!
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! explicit calculation
        S1=0.0d0
        
        do k=0,ng
          z=-lg+(2.0d0*lg*k)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do i=0,ng
              x=-lg+(2.0d0*lg*i)/ng
              V  = 1.0d0/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
              tmp= x**nx * y**ny * z**nz * exp( -b * (x**2+y**2+z**2) )
              tmp=tmp * (x-r2(1))**nx2 * (y-r2(2))**ny2 * (z-r2(3))**nz2 * exp( -d * ((x-r2(1))**2+(y-r2(2))**2+(z-r2(3))**2) )
              S1 = S1 + tmp * V
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        
        ! get corresponding overlap
        S1_= CC_Coulomb_Ion(d,r2,nx2,ny2,nz2,b,r1,nx,ny,nz,(/x0,y0,z0/))
        
        ! get corresponding symmetric overlap
        S1__= CC_Coulomb_Ion(b,r1,nx,ny,nz,d,r2,nx2,ny2,nz2,(/x0,y0,z0/))
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S1,'  CC_Coulomb_Ion result=',S1_,S1__
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-5 ) then
          print *,'Error: Coulomb differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1__)/S1) > 1.0d-5 ) then
          print *,'Error: Coulomb differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Coulomb passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''
  
  !!
  !! check for Y routine
  !!
  do l=0,lmax
    do m=-l,l
      ! explicit calculation
      S1=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            V  = 1.0d0/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
            tmp=Y_Value(b,l,m,(/x,y,z/))
            tmp=tmp*Y_Value(d,l2,m2,(/x-r2(1),y-r2(2),z-r2(3)/))
            S1 = S1 + tmp * V
          end do
        end do
      end do
      S1=S1*8.0d0*lg**3/ng**3
      
      ! get corresponding interaction
      S1_= YY_Coulomb_Ion(d,r2,l2,m2,b,r1,l,m,(/x0,y0,z0/))
      
      ! get corresponding symmetric interaction
      S1__= YY_Coulomb_Ion(b,r1,l,m,d,r2,l2,m2,(/x0,y0,z0/))
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S1,'  YY_Coulomb_Ion result=',S1_,S1__
        
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-5 ) then
        print *,'Error: Coulomb differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1__)/S1) > 1.0d-5 ) then
        print *,'Error: Coulomb differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Coulomb passed up to lmax=',lmax
  print *,''
  
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! test three center case
        S1=C_Coulomb_C(b,r1,nx,ny,nz,d,r2,nx2,ny2,nz2)
        S1_=CC_Coulomb_C(b,r1,nx,ny,nz,0.0d0,(/x0,y0,z0/),0,0,0,d,r2,nx2,ny2,nz2)
        S1__=CC_Coulomb_C(0.0d0,(/x0,y0,z0/),0,0,0,b,r1,nx,ny,nz,d,r2,nx2,ny2,nz2)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  C_Coulomb_C=',S1,'  CC_Coulomb_C  result=',S1_,S1__
        
        if ( abs((S1-S1_)/S1) > 1.0d-10 .or. abs((S1-S1__)/S1) > 1.0d-10 ) then
          print *,'Error: Three center integrals not consistent with two center case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed'
          print *,''
          stop 1
        end if
        
        ! test four center case
        S1=C_Coulomb_C(b,r1,nx,ny,nz,d,r2,nx2,ny2,nz2)
        
        S1_=CC_Coulomb_CC(b,r1,nx,ny,nz, &
                          0.0d0,(/x0,y0,z0/),0,0,0, &
                          d,r2,nx2,ny2,nz2, &
                          0.0d0,(/x0,y0,z0/),0,0,0)
                          
        S1__=CC_Coulomb_CC(0.0d0,(/x0,y0,z0/),0,0,0,b,r1,nx,ny,nz,d,r2,nx2,ny2,nz2,0.0d0,(/x0,y0,z0/),0,0,0)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  C_Coulomb_C=',S1,'  CC_Coulomb_CC result=',S1_,S1__
        
        if ( abs((S1-S1_)/S1) > 1.0d-10 .or. abs((S1-S1__)/S1) > 1.0d-10 ) then
          print *,'Error: Four center integrals not consistent with two center case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed ici'
          print *,''
          stop 1
        end if
        
        
        S1=C_Coulomb_C(b,r1,nx,ny,nz,d,r2,nx2,ny2,nz2)
        S1_=CC_Coulomb_CC(b,r1,nx,ny,nz,0.0d0,(/x0,y0,z0/),0,0,0,0.0d0,(/x0,y0,z0/),0,0,0,d,r2,nx2,ny2,nz2)
        S1__=CC_Coulomb_CC(0.0d0,(/x0,y0,z0/),0,0,0,b,r1,nx,ny,nz,0.0d0,(/x0,y0,z0/),0,0,0,d,r2,nx2,ny2,nz2)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  C_Coulomb_C=',S1,'  CC_Coulomb_CC result=',S1_,S1__
        
        if ( abs((S1-S1_)/S1) > 1.0d-10 .or. abs((S1-S1__)/S1) > 1.0d-10 ) then
          print *,'Error: Four center integrals not consistent with two center case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Coulomb failed la'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Coulomb passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''

  !!
  !! check for Y routine
  !!
  do l=0,lmax
    do m=-l,l
      
      ! test three center case
      S1=Y_Coulomb_Y(b,r1,l,m,d,r2,l2,m2)
      S1_=YY_Coulomb_Y(b,r1,l,m,0.0d0,(/x0,y0,z0/),0,0,d,r2,l2,m2)/(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))
      S1__=YY_Coulomb_Y(0.0d0,(/x0,y0,z0/),0,0,b,r1,l,m,d,r2,l2,m2)/(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  Y_Coulomb_Y=',S1,'  YY_Coulomb_Y  result=',S1_,S1__
        
      if ( abs((S1-S1_)/S1) > 1.0d-8 .or. abs((S1-S1__)/S1) > 1.0d-8 ) then
        print *,'Error: Three center integrals not consistent with two center case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      
      ! test four center case
      S1=Y_Coulomb_Y(b,r1,l,m,d,r2,l2,m2)
      S1_=YY_Coulomb_YY(b,r1,l,m,0.0d0,(/x0,y0,z0/),0,0,d,r2,l2,m2,0.0d0,(/x0,y0,z0/),0,0)&
         /(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))/(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))
      S1__=YY_Coulomb_YY(b,r1,l,m,0.0d0,(/x0,y0,z0/),0,0,0.0d0,(/x0,y0,z0/),0,0,d,r2,l2,m2)&
         /(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))/(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  Y_Coulomb_Y=',S1,'  YY_Coulomb_YY result=',S1_,S1__
        
      if ( abs((S1-S1_)/S1) > 1.0d-10 .or. abs((S1-S1__)/S1) > 1.0d-10 ) then
        print *,'Error: Four center integrals not consistent with two center case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      S1_=YY_Coulomb_YY(0.0d0,(/x0,y0,z0/),0,0,b,r1,l,m,d,r2,l2,m2,0.0d0,(/x0,y0,z0/),0,0)&
         /(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))/(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))
      S1__=YY_Coulomb_YY(0.0d0,(/x0,y0,z0/),0,0,b,r1,l,m,0.0d0,(/x0,y0,z0/),0,0,d,r2,l2,m2)&
         /(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))/(1.0d0/2.0d0/sqrt(4.0d0*atan(1.0d0)))
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  Y_Coulomb_Y=',S1,'  YY_Coulomb_YY result=',S1_,S1__
        
      if ( abs((S1-S1_)/S1) > 1.0d-10 .or. abs((S1-S1__)/S1) > 1.0d-10 ) then
        print *,'Error: Four center integrals not consistent with two center case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Coulomb failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Coulomb passed up to lmax=',lmax
  print *,''
  
  ! check consistency for two and three centers routine
  S1=Y_Coulomb_C(b,r2,l2,m2,e,r3,nx3,ny3,nz3)
  S2=C_Coulomb_Y(e,r3,nx3,ny3,nz3,b,r2,l2,m2)
  ! test result
  if ( abs((S1-S2)/S1) > 1.0d-8 ) then
    print *,S1,S2
    print *,'Error: <YC> Coulomb failed'
    print *,''
    print *,'test_Overlap failed'
    print *,''
    stop 1
  end if
  
  ! check consistency for two and three centers routine
  S1=YC_Coulomb_Y(d,r1,l1,m1,e,r3,nx3,ny3,nz3,b,r2,l2,m2)
  S2=CY_Coulomb_Y(e,r3,nx3,ny3,nz3,d,r1,l1,m1,b,r2,l2,m2)
  ! test result
  if ( abs((S1-S2)/S1) > 1.0d-8 ) then
    print *,S1,S2,S3
    print *,'Error: <YYC> Three center Coulomb failed'
    print *,''
    print *,'test_Coulomb failed'
    print *,''
    stop 1
  end if
  
  ! check consistency for two and three centers routine
  S1=YC_Coulomb_C(d,r1,l1,m1,b,r2,nx2,ny2,nz2,e,r3,nx3,ny3,nz3)
  S2=CY_Coulomb_C(b,r2,nx2,ny2,nz2,d,r1,l1,m1,e,r3,nx3,ny3,nz3)
  ! test result
  if ( abs((S1-S2)/S1) > 1.0d-8 ) then
    print *,S1,S2,S3
    print *,'Error: <CCY> Three center Coulomb failed'
    print *,''
    print *,'test_Coulomb failed'
    print *,''
    stop 1
  end if
  
  print *,''
  write (*,'(A)') 'test_Coulomb passed'
  print *,''
  
end program

