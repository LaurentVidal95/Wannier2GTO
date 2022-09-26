program test_Overlap
  
  implicit none
  
  ! local variables
  integer, parameter      :: nmax = 2
  integer, parameter      :: ng   = 63
  real(kind=8), parameter :: a    = 2.0d0
  real(kind=8), parameter :: b    = 3.5d0
  real(kind=8), parameter :: c    = 1.5d0
  real(kind=8), parameter :: lg   = 5.0d0
  integer, parameter      :: nx1  = 0
  integer, parameter      :: ny1  = 0
  integer, parameter      :: nz1  = 0
  integer, parameter      :: nx2  = 1
  integer, parameter      :: ny2  = 0
  integer, parameter      :: nz2  = 2
  integer, parameter      :: nx3  = 0
  integer, parameter      :: ny3  = 1
  integer, parameter      :: nz3  = 1
  integer, parameter      :: nx4  = 2
  integer, parameter      :: ny4  = 2
  integer, parameter      :: nz4  = 1
  integer, parameter      :: l1   = 0
  integer, parameter      :: m1   = 0
  integer, parameter      :: l2   = 1
  integer, parameter      :: m2   = 0
  integer, parameter      :: l3   = 3
  integer, parameter      :: m3   = 2
  integer, parameter      :: l4   = 6
  integer, parameter      :: m4   =-1
  integer, parameter      :: lmax = 5
  real(kind=8), parameter :: r1(3) = (/ 0.2d0,-0.5d0, 0.3d0 /)
  real(kind=8), parameter :: r2(3) = (/ 0.2d0, 0.5d0, 0.1d0 /)
  real(kind=8), parameter :: r3(3) = (/ 0.1d0, 0.0d0, 0.5d0 /)
  real(kind=8), parameter :: r4(3) = (/-0.2d0,-0.1d0, 0.4d0 /)
  integer      :: i
  integer      :: j
  integer      :: k
  integer      :: l
  integer      :: ll1
  integer      :: ll2
  integer      :: m
  integer      :: nx
  integer      :: ny
  integer      :: nz
  real(kind=8) :: x
  real(kind=8) :: y
  real(kind=8) :: z
  real(kind=8) :: d1
  real(kind=8) :: d2
  real(kind=8) :: tmp
  real(kind=8) :: S1,S2,S3,S4
  real(kind=8) :: S1_,S2_,S3_,S4_
  real(kind=8) :: S1__,S2__,S3__,S4__
  real(kind=8) :: S1___,S2___,S3___,S4___
  real(kind=8) :: C_Overlap_C
  real(kind=8) :: Y_Overlap_Y
  real(kind=8) :: C_Overlap_Y
  real(kind=8) :: Y_Overlap_C
  real(kind=8) :: CC_Overlap_C
  real(kind=8) :: YY_Overlap_Y
  real(kind=8) :: YY_Overlap_C
  real(kind=8) :: YC_Overlap_Y
  real(kind=8) :: CY_Overlap_Y
  real(kind=8) :: CC_Overlap_Y
  real(kind=8) :: CY_Overlap_C
  real(kind=8) :: YC_Overlap_C
  real(kind=8) :: YY_Overlap_YY
  real(kind=8) :: CC_Overlap_CC
  
  real(kind=8) :: Y_Value
  real(kind=8) :: Overlap_upper_bound
  
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! explicit calculation
        S1=0.0d0
        S2=0.0d0
        S3=0.0d0
        S4=0.0d0
        do i=0,ng
          x=-lg+(2.0d0*lg*i)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do k=0,ng
              z=-lg+(2.0d0*lg*k)/ng
              S1 = S1 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r1(1))**nx1 * (y-r1(2))**ny1 * (z-r1(3))**nz1 &
                      * exp( -b * ( (x-r1(1))**2 + (y-r1(2))**2 + (z-r1(3))**2 ) )
              S2 = S2 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r2(1))**nx2 * (y-r2(2))**ny2 * (z-r2(3))**nz2 &
                      * exp( -b * ( (x-r2(1))**2 + (y-r2(2))**2 + (z-r2(3))**2 ) )
              S3 = S3 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r3(1))**nx3 * (y-r3(2))**ny3 * (z-r3(3))**nz3 &
                      * exp( -b * ( (x-r3(1))**2 + (y-r3(2))**2 + (z-r3(3))**2 ) )
              S4 = S4 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r4(1))**nx4 * (y-r4(2))**ny4 * (z-r4(3))**nz4 &
                      * exp( -b * ( (x-r4(1))**2 + (y-r4(2))**2 + (z-r4(3))**2 ) )
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        S2=S2*8.0d0*lg**3/ng**3
        S3=S3*8.0d0*lg**3/ng**3
        S4=S4*8.0d0*lg**3/ng**3
        
        ! get corresponding overlap
        S1_= C_Overlap_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r1,nx1,ny1,nz1)
        S2_= C_Overlap_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r2,nx2,ny2,nz2)
        S3_= C_Overlap_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r3,nx3,ny3,nz3)
        S4_= C_Overlap_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r4,nx4,ny4,nz4)
        
        ! get corresponding symmetric overlap
        S1__= C_Overlap_C(b,r1,nx1,ny1,nz1,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        S2__= C_Overlap_C(b,r2,nx2,ny2,nz2,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        S3__= C_Overlap_C(b,r3,nx3,ny3,nz3,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        S4__= C_Overlap_C(b,r4,nx4,ny4,nz4,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit overlap=',S1,'  C_Overlap_C result=',S1_
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit overlap=',S2,'  C_Overlap_C result=',S2_
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit overlap=',S3,'  C_Overlap_C result=',S3_
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit overlap=',S4,'  C_Overlap_C result=',S4_
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S3-S3_)/S3) > 1.0d-8 .and. ( abs(S3)>1.0d-10 .or. abs(S3_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_R_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S4-S4_)/S4) > 1.0d-8 .and. ( abs(S4)>1.0d-10 .or. abs(S4_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1__)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1__)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2__)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2__)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S3-S3__)/S3) > 1.0d-8 .and. ( abs(S3)>1.0d-10 .or. abs(S3__)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S4-S4__)/S4) > 1.0d-8 .and. ( abs(S4)>1.0d-10 .or. abs(S4__)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Overlap passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''
  
  !!
  !! check for Y routine
  !!
  do l=0,lmax
    do m=-l,l
      ! explicit calculation
      S1=0.0d0
      S2=0.0d0
      S3=0.0d0
      S4=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            
            tmp=Y_Value(a,l,m,(/x,y,z/))
            
            S1=S1+Y_Value(b,l1,m1,(/x-r1(1),y-r1(2),z-r1(3)/))*tmp
            S2=S2+Y_Value(b,l2,m2,(/x-r2(1),y-r2(2),z-r2(3)/))*tmp
            S3=S3+Y_Value(b,l3,m3,(/x-r3(1),y-r3(2),z-r3(3)/))*tmp
            S4=S4+Y_Value(b,l4,m4,(/x-r4(1),y-r4(2),z-r4(3)/))*tmp
            
          end do
        end do
      end do
      
      S1=S1*8.0d0*lg**3/ng**3
      S2=S2*8.0d0*lg**3/ng**3
      S3=S3*8.0d0*lg**3/ng**3
      S4=S4*8.0d0*lg**3/ng**3
      
      ! get corresponding Overlap
      S1_= Y_Overlap_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r1,l1,m1)
      S2_= Y_Overlap_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r2,l2,m2)
      S3_= Y_Overlap_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r3,l3,m3)
      S4_= Y_Overlap_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r4,l4,m4)
      
      ! get corresponding symmetric Overlap
      S1__= Y_Overlap_Y(b,r1,l1,m1,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      S2__= Y_Overlap_Y(b,r2,l2,m2,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      S3__= Y_Overlap_Y(b,r3,l3,m3,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      S4__= Y_Overlap_Y(b,r4,l4,m4,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit overlap=',S1,'  Y_Overlap_Y result=',S1_,S1__
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit overlap=',S2,'  Y_Overlap_Y result=',S2_,S2__
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit overlap=',S3,'  Y_Overlap_Y result=',S3_,S3__
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit overlap=',S4,'  Y_Overlap_Y result=',S4_,S4__
      
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S3-S3_)/S3) > 1.0d-8 .and. ( abs(S3)>1.0d-10 .or. abs(S3_)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S4-S4_)/S4) > 1.0d-8 .and. ( abs(S4)>1.0d-10 .or. abs(S4_)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1__)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2__)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S3-S3__)/S3) > 1.0d-8 .and. ( abs(S3)>1.0d-10 .or. abs(S3__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S4-S4__)/S4) > 1.0d-8 .and. ( abs(S4)>1.0d-10 .or. abs(S4__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Overlap passed up to lmax=',lmax
  print *,''
  
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! explicit calculation
        S1=0.0d0
        S2=0.0d0
        do i=0,ng
          x=-lg+(2.0d0*lg*i)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do k=0,ng
              z=-lg+(2.0d0*lg*k)/ng
              S1 = S1 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r1(1))**nx1 * (y-r1(2))**ny1 * (z-r1(3))**nz1 &
                      * exp( -b * ( (x-r1(1))**2 + (y-r1(2))**2 + (z-r1(3))**2 ) ) &
                      * (x-r2(1))**nx2 * (y-r2(2))**ny2 * (z-r2(3))**nz2 &
                      * exp( -c * ( (x-r2(1))**2 + (y-r2(2))**2 + (z-r2(3))**2 ) )
              S2 = S2 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r3(1))**nx3 * (y-r3(2))**ny3 * (z-r3(3))**nz3 &
                      * exp( -b * ( (x-r3(1))**2 + (y-r3(2))**2 + (z-r3(3))**2 ) ) &
                      * (x-r4(1))**nx4 * (y-r4(2))**ny4 * (z-r4(3))**nz4 &
                      * exp( -c * ( (x-r4(1))**2 + (y-r4(2))**2 + (z-r4(3))**2 ) )
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        S2=S2*8.0d0*lg**3/ng**3
        
        ! get corresponding overlap
        S1_= CC_Overlap_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r1,nx1,ny1,nz1,c,r2,nx2,ny2,nz2)
        S2_= CC_Overlap_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r3,nx3,ny3,nz3,c,r4,nx4,ny4,nz4)
        
        ! get corresponding symmetric overlap
        S1__= CC_Overlap_C(b,r1,nx1,ny1,nz1,c,r2,nx2,ny2,nz2,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        S2__= CC_Overlap_C(b,r3,nx3,ny3,nz3,c,r4,nx4,ny4,nz4,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        
        ! get corresponding symmetric overlap
        S1___= CC_Overlap_C(c,r2,nx2,ny2,nz2,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r1,nx1,ny1,nz1)
        S2___= CC_Overlap_C(c,r4,nx4,ny4,nz4,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r3,nx3,ny3,nz3)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit overlap=',S1,'  CC_Overlap_C result=',S1_,S1__,S1___
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit overlap=',S2,'  CC_Overlap_C result=',S2_,S2__,S2___
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1__)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2__)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1___)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2___)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Overlap passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''
  
  
  !!
  !! check for Y routine
  !!
  do l=0,lmax
    do m=-l,l
      ! explicit calculation
      S1=0.0d0
      S2=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            
            tmp=Y_Value(a,l,m,(/x,y,z/))
            S1=S1+Y_Value(b,l1,m1,(/x-r1(1),y-r1(2),z-r1(3)/))*Y_Value(c,l2,m2,(/x-r2(1),y-r2(2),z-r2(3)/))*tmp
            S2=S2+Y_Value(b,l3,m3,(/x-r3(1),y-r3(2),z-r3(3)/))*Y_Value(c,l4,m4,(/x-r4(1),y-r4(2),z-r4(3)/))*tmp
            
          end do
        end do
      end do
      
      S1=S1*8.0d0*lg**3/ng**3
      S2=S2*8.0d0*lg**3/ng**3
      
      ! get corresponding Overlap
      S1_= YY_Overlap_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r1,l1,m1,c,r2,l2,m2)
      S2_= YY_Overlap_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r3,l3,m3,c,r4,l4,m4)
      
      ! get corresponding symmetric Overlap
      S1__= YY_Overlap_Y(b,r1,l1,m1,c,r2,l2,m2,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      S2__= YY_Overlap_Y(b,r3,l3,m3,c,r4,l4,m4,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      
      ! get corresponding symmetric Overlap
      S1___= YY_Overlap_Y(c,r2,l2,m2,a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r1,l1,m1)
      S2___= YY_Overlap_Y(c,r4,l4,m4,a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r3,l3,m3)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit overlap=',S1,'  YY_Overlap_Y result=',S1_,S1__,S1___
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit overlap=',S2,'  YY_Overlap_Y result=',S2_,S2__,S2___
      
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1__)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2__)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1___)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2___)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Overlap passed up to lmax=',lmax
  print *,''
  
  
  !!
  !! check for Overlap upper bound routine
  !!
  do ll1=0,6,2
    do ll2=0,ll1,2
      ! explicit calculation
      S1=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            d1=sum((/x-r1(1),y-r1(2),z-r1(3)/)**2)
            d2=sum((/x-r2(1),y-r2(2),z-r2(3)/)**2)
            S1=S1+ d1**(ll1/2)*exp(-a*d1) * d2**(ll2/2)*exp(-b*d2)
          end do
        end do
      end do
      
      S1=S1*8.0d0*lg**3/ng**3
      
      ! get corresponding Overlap
      S1_= Overlap_upper_bound(a,r1,ll1,b,r2,ll2)
      
      ! get corresponding symmetric Overlap
      S1__= Overlap_upper_bound(b,r2,ll2,a,r1,ll1)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l1=',ll1,'  l2=',ll2,'  explicit bound=',S1,'  Overlap_upper_bound result=',S1_,S1__
      
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
        print *,'Error: Overlap_upper_bound differs from explicit calculation for case'
        print *,'       l1 =',ll1
        print *,'       l1 =',ll2
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1)') 'test_Overlap passed up to lmax=',ll1
  print *,''
  
  ! check consistency for two and three centers routine
  S1=Y_Overlap_C(b,r2,l2,m2,c,r3,nx3,ny3,nz3)
  S2=C_Overlap_Y(c,r3,nx3,ny3,nz3,b,r2,l2,m2)
  ! test result
  if ( abs((S1-S2)/S1) > 1.0d-8 ) then
    print *,S1,S2
    print *,'Error: <YC> Three center overlap failed'
    print *,''
    print *,'test_Overlap failed'
    print *,''
    stop 1
  end if
  
  ! check consistency for two and three centers routine
  S1=YY_Overlap_C(a,r1,l1,m1,b,r2,l2,m2,c,r3,nx3,ny3,nz3)
  S2=YC_Overlap_Y(a,r1,l1,m1,c,r3,nx3,ny3,nz3,b,r2,l2,m2)
  S3=CY_Overlap_Y(c,r3,nx3,ny3,nz3,b,r2,l2,m2,a,r1,l1,m1)
  ! test result
  if ( abs((S1-S2)/S1) > 1.0d-8 .or. abs((S1-S3)/S1) > 1.0d-8 ) then
    print *,S1,S2,S3
    print *,'Error: <YYC> Three center overlap failed'
    print *,''
    print *,'test_Overlap failed'
    print *,''
    stop 1
  end if
  
  ! check consistency for two and three centers routine
  S1=YC_Overlap_C(a,r1,l1,m1,b,r2,nx2,ny2,nz2,c,r3,nx3,ny3,nz3)
  S2=CY_Overlap_C(b,r2,nx2,ny2,nz2,a,r1,l1,m1,c,r3,nx3,ny3,nz3)
  S3=CC_Overlap_Y(b,r2,nx2,ny2,nz2,c,r3,nx3,ny3,nz3,a,r1,l1,m1)
  ! test result
  if ( abs((S1-S2)/S1) > 1.0d-8 .or. abs((S1-S3)/S1) > 1.0d-8 ) then
    print *,S1,S2,S3
    print *,'Error: <CCY> Three center overlap failed'
    print *,''
    print *,'test_Overlap failed'
    print *,''
    stop 1
  end if
  
  print *,''
  write (*,'(A)') 'test_Overlap passed'
  print *,''
  
  
  
  !!
  !! check for Y routine
  !!
  do l=0,4
    do m=-l,l
      ! explicit calculation
      S1=0.0d0
      S2=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            
            tmp=Y_Value(a,l,m,(/x,y,z/))
            S1=S1+Y_Value(b,l1,m1,(/x-r1(1),y-r1(2),z-r1(3)/)) &
                 *Y_Value(c,l2,m2,(/x-r2(1),y-r2(2),z-r2(3)/)) &
                 *Y_Value(b,l3,m3,(/x-r3(1),y-r3(2),z-r3(3)/)) &
                 *tmp
            S2=S2+Y_Value(c,l2,m2,(/x-r2(1),y-r2(2),z-r2(3)/)) &
                 *Y_Value(b,l3,m3,(/x-r3(1),y-r3(2),z-r3(3)/)) &
                 *Y_Value(c,l4,m4,(/x-r4(1),y-r4(2),z-r4(3)/)) &
                 *tmp
            
          end do
        end do
      end do
      
      S1=S1*8.0d0*lg**3/ng**3
      S2=S2*8.0d0*lg**3/ng**3
      
      ! get corresponding Overlap
      S1_= YY_Overlap_YY(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r1,l1,m1,c,r2,l2,m2,b,r3,l3,m3)
      S2_= YY_Overlap_YY(a,(/0.0d0,0.0d0,0.0d0/),l,m,c,r2,l2,m2,b,r3,l3,m3,c,r4,l4,m4)
      
      ! get corresponding symmetric Overlap
      S1__= YY_Overlap_YY(b,r1,l1,m1,c,r2,l2,m2,b,r3,l3,m3,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      S2__= YY_Overlap_YY(c,r2,l2,m2,b,r3,l3,m3,c,r4,l4,m4,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      
      ! get corresponding symmetric Overlap
      S1___= YY_Overlap_YY(c,r2,l2,m2,a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r1,l1,m1,b,r3,l3,m3)
      S2___= YY_Overlap_YY(c,r4,l4,m4,a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r3,l3,m3,c,r2,l2,m2)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit overlap=',S1,'  YY_Overlap_YY result=',S1_,S1__,S1___
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
        'l=',l,'  m=',m,'  explicit overlap=',S2,'  YY_Overlap_YY result=',S2_,S2__,S2___
      
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1__)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2__)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1___)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2___)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2__)>1.0d-10 ) ) then
        print *,'Error: Overlap differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Overlap failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Overlap passed up to lmax=',lmax
  print *,''
  
  
  
  
  
  do nx=0,4
    do ny=0,4-nx
      do nz=0,4-nx-ny
        
        ! explicit calculation
        S1=0.0d0
        S2=0.0d0
        do i=0,ng
          x=-lg+(2.0d0*lg*i)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do k=0,ng
              z=-lg+(2.0d0*lg*k)/ng
              S1 = S1 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r1(1))**nx1 * (y-r1(2))**ny1 * (z-r1(3))**nz1 &
                      * exp( -b * ( (x-r1(1))**2 + (y-r1(2))**2 + (z-r1(3))**2 ) ) &
                      * (x-r2(1))**nx2 * (y-r2(2))**ny2 * (z-r2(3))**nz2 &
                      * exp( -c * ( (x-r2(1))**2 + (y-r2(2))**2 + (z-r2(3))**2 ) ) &
                      * (x-r3(1))**nx3 * (y-r3(2))**ny3 * (z-r3(3))**nz3 &
                      * exp( -b * ( (x-r3(1))**2 + (y-r3(2))**2 + (z-r3(3))**2 ) ) 
              S2 = S2 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r2(1))**nx2 * (y-r2(2))**ny2 * (z-r2(3))**nz2 &
                      * exp( -c * ( (x-r2(1))**2 + (y-r2(2))**2 + (z-r2(3))**2 ) ) &
                      * (x-r3(1))**nx3 * (y-r3(2))**ny3 * (z-r3(3))**nz3 &
                      * exp( -b * ( (x-r3(1))**2 + (y-r3(2))**2 + (z-r3(3))**2 ) ) &
                      * (x-r4(1))**nx4 * (y-r4(2))**ny4 * (z-r4(3))**nz4 &
                      * exp( -c * ( (x-r4(1))**2 + (y-r4(2))**2 + (z-r4(3))**2 ) )
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        S2=S2*8.0d0*lg**3/ng**3
        
        ! get corresponding overlap
        S1_= CC_Overlap_CC(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r1,nx1,ny1,nz1,c,r2,nx2,ny2,nz2,b,r3,nx3,ny3,nz3)
        S2_= CC_Overlap_CC(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,c,r2,nx2,ny2,nz2,b,r3,nx3,ny3,nz3,c,r4,nx4,ny4,nz4)
        
        ! get corresponding symmetric overlap
        S1__= CC_Overlap_CC(b,r1,nx1,ny1,nz1,c,r2,nx2,ny2,nz2,b,r3,nx3,ny3,nz3,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        S2__= CC_Overlap_CC(c,r2,nx2,ny2,nz2,b,r3,nx3,ny3,nz3,c,r4,nx4,ny4,nz4,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        
        ! get corresponding symmetric overlap
        S1___= CC_Overlap_CC(c,r2,nx2,ny2,nz2,b,r3,nx3,ny3,nz3,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r1,nx1,ny1,nz1)
        S2___= CC_Overlap_CC(c,r4,nx4,ny4,nz4,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r3,nx3,ny3,nz3,c,r2,nx2,ny2,nz2)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit overlap=',S1,'  CC_Overlap_CC result=',S1_,S1__,S1___
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit overlap=',S2,'  CC_Overlap_CC result=',S2_,S2__,S2___
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1__)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2__)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1___)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2___)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: Overlap differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Overlap failed'
          print *,''
          stop 1
        end if
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Overlap passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''

  
  
  
end program 

