program test_Laplacian
  
  implicit none
  
  ! local variables
  integer, parameter      :: nmax = 2
  integer, parameter      :: ng   = 63
  real(kind=8), parameter :: a    = 2.0d0
  real(kind=8), parameter :: b    = 3.5d0
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
  integer, parameter      :: lmax = 6
  real(kind=8), parameter :: delta= 0.0001d0
  real(kind=8), parameter :: r1(3) = (/ 0.2d0,-0.5d0, 0.3d0 /)
  real(kind=8), parameter :: r2(3) = (/ 0.2d0, 0.5d0, 0.1d0 /)
  real(kind=8), parameter :: r3(3) = (/ 0.1d0, 0.0d0, 0.5d0 /)
  real(kind=8), parameter :: r4(3) = (/-0.2d0,-0.1d0, 0.4d0 /)
  integer      :: l
  integer      :: m
  integer      :: i
  integer      :: j
  integer      :: k
  integer      :: nx
  integer      :: ny
  integer      :: nz
  real(kind=8) :: x
  real(kind=8) :: y
  real(kind=8) :: z
  real(kind=8) :: tmp
  real(kind=8) :: tmp_x
  real(kind=8) :: tmp_y
  real(kind=8) :: tmp_z
  real(kind=8) :: S1,S2,S3,S4
  real(kind=8) :: S1_,S2_,S3_,S4_
  real(kind=8) :: S1__,S2__,S3__,S4__
  real(kind=8) :: C_Laplacian_C
  real(kind=8) :: Y_Laplacian_Y
  real(kind=8) :: Y_Laplacian_C
  real(kind=8) :: C_Laplacian_Y
  real(kind=8) :: Y_Value
  real(kind=8) :: V_000
  real(kind=8) :: V_p00
  real(kind=8) :: V_0p0
  real(kind=8) :: V_00p
  real(kind=8) :: V_m00
  real(kind=8) :: V_0m0
  real(kind=8) :: V_00m
  
  !!
  !! check for C routine
  !!
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
              
              if ( ny==0 ) then
                tmp_y = 2*a*x**nx*(2*a*y**2-1)*z**nz*exp(-a*(z**2+y**2+x**2))
              else if ( ny==1 ) then
                tmp_y = 2*a*y*(2*a*y**2-3)*x**nx*z**nz*exp(-a*(z**2+y**2+x**2))
              else
                tmp_y = x**nx*y**(ny-2)*(2*a*y**2*(2*a*y**2-(2*ny+1))+(ny-1)*ny)*z**nz*exp(-a*(z**2+y**2+x**2))
              end if 
              
              if ( nx==0 ) then
                tmp_x = 2*a*(2*a*x**2-1)*y**ny*z**nz*exp(-a*(z**2+y**2+x**2))
              else if ( nx==1 ) then
                tmp_x = 2*a*x*(2*a*x**2-3)*y**ny*z**nz*exp(-a*(z**2+y**2+x**2))
              else
                tmp_x = x**(nx-2)*(2*a*x**2*(2*a*x**2-(2*nx+1))+(nx-1)*nx)*y**ny*z**nz*exp(-a*(z**2+y**2+x**2))
              end if 
              
              if ( nz==0 ) then
                tmp_z = 2*a*x**nx*y**ny*(2*a*z**2-1)*exp(-a*(z**2+y**2+x**2))
              else if ( nz==1 ) then
                tmp_z = 2*a*z*(2*a*z**2-3)*x**nx*y**ny*exp(-a*(z**2+y**2+x**2))
              else
                tmp_z = x**nx*y**ny*z**(nz-2)*(2*a*z**2*(2*a*z**2-(2*nz+1))+(nz-1)*nz)*exp(-a*(z**2+y**2+x**2))
              end if 
              
              ! compute laplacian
              tmp = tmp_x + tmp_y + tmp_z
              
              S1 = S1 + tmp * (x-r1(1))**nx1 * (y-r1(2))**ny1 * (z-r1(3))**nz1 &
                      * exp( -b * ( (x-r1(1))**2 + (y-r1(2))**2 + (z-r1(3))**2 ) )
              S2 = S2 + tmp * (x-r2(1))**nx2 * (y-r2(2))**ny2 * (z-r2(3))**nz2 &
                      * exp( -b * ( (x-r2(1))**2 + (y-r2(2))**2 + (z-r2(3))**2 ) )
              S3 = S3 + tmp * (x-r3(1))**nx3 * (y-r3(2))**ny3 * (z-r3(3))**nz3 &
                      * exp( -b * ( (x-r3(1))**2 + (y-r3(2))**2 + (z-r3(3))**2 ) )
              S4 = S4 + tmp * (x-r4(1))**nx4 * (y-r4(2))**ny4 * (z-r4(3))**nz4 &
                      * exp( -b * ( (x-r4(1))**2 + (y-r4(2))**2 + (z-r4(3))**2 ) )
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        S2=S2*8.0d0*lg**3/ng**3
        S3=S3*8.0d0*lg**3/ng**3
        S4=S4*8.0d0*lg**3/ng**3
        
        ! get corresponding laplacian
        S1_= C_Laplacian_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r1,nx1,ny1,nz1)
        S2_= C_Laplacian_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r2,nx2,ny2,nz2)
        S3_= C_Laplacian_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r3,nx3,ny3,nz3)
        S4_= C_Laplacian_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r4,nx4,ny4,nz4)
        
        ! get corresponding symmetric laplacian
        S1__= C_Laplacian_C(b,r1,nx1,ny1,nz1,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        S2__= C_Laplacian_C(b,r2,nx2,ny2,nz2,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        S3__= C_Laplacian_C(b,r3,nx3,ny3,nz3,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        S4__= C_Laplacian_C(b,r4,nx4,ny4,nz4,a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit laplacian=',S1,'  C_Laplacian_C result=',S1_,S1__
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit laplacian=',S2,'  C_Laplacian_C result=',S2_,S2__
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit laplacian=',S3,'  C_Laplacian_C result=',S3_,S3__
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit laplacian=',S4,'  C_Laplacian_C result=',S4_,S4__
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: Laplacian differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Laplacian failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: Laplacian differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Laplacian failed'
          print *,''
          stop 1
        end if
        if ( abs((S3-S3_)/S3) > 1.0d-8 .and. ( abs(S3)>1.0d-10 .or. abs(S3_)>1.0d-10 ) ) then
          print *,'Error: Laplacian differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Laplacian failed'
          print *,''
          stop 1
        end if
        if ( abs((S4-S4_)/S4) > 1.0d-8 .and. ( abs(S4)>1.0d-10 .or. abs(S4_)>1.0d-10 ) ) then
          print *,'Error: Laplacian differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Laplacian failed'
          print *,''
          stop 1
        end if
        if ( abs((S1-S1__)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1__)>1.0d-10 ) ) then
          print *,'Error: Laplacian differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Laplacian failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2__)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2__)>1.0d-10 ) ) then
          print *,'Error: Laplacian differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Laplacian failed'
          print *,''
          stop 1
        end if
        if ( abs((S3-S3__)/S3) > 1.0d-8 .and. ( abs(S3)>1.0d-10 .or. abs(S3__)>1.0d-10 ) ) then
          print *,'Error: Laplacian differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Laplacian failed'
          print *,''
          stop 1
        end if
        if ( abs((S4-S4__)/S4) > 1.0d-8 .and. ( abs(S4)>1.0d-10 .or. abs(S4__)>1.0d-10 ) ) then
          print *,'Error: Laplacian differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Laplacian failed'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Laplacian passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
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
            
            V_000=Y_Value(a,l,m,(/x,y,z/))
            V_p00=Y_Value(a,l,m,(/x+delta,y,z/))
            V_m00=Y_Value(a,l,m,(/x-delta,y,z/))
            V_0p0=Y_Value(a,l,m,(/x,y+delta,z/))
            V_0m0=Y_Value(a,l,m,(/x,y-delta,z/))
            V_00p=Y_Value(a,l,m,(/x,y,z+delta/))
            V_00m=Y_Value(a,l,m,(/x,y,z-delta/))
            tmp=(V_p00+V_0p0+V_00p+V_m00+V_0m0+V_00m-6.0d0*V_000)/delta**2
            
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
      
      ! get corresponding laplacian
      S1_= Y_Laplacian_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r1,l1,m1)
      S2_= Y_Laplacian_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r2,l2,m2)
      S3_= Y_Laplacian_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r3,l3,m3)
      S4_= Y_Laplacian_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r4,l4,m4)
      
      ! get corresponding symmetric laplacian
      S1__= Y_Laplacian_Y(b,r1,l1,m1,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      S2__= Y_Laplacian_Y(b,r2,l2,m2,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      S3__= Y_Laplacian_Y(b,r3,l3,m3,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      S4__= Y_Laplacian_Y(b,r4,l4,m4,a,(/0.0d0,0.0d0,0.0d0/),l,m)
      
      write (*,'(A,I2,A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'l=',l1,'  m=',m1,'  explicit laplacian=',S1,'  Y_Laplacian_Y result=',S1_,S1__
      write (*,'(A,I2,A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'l=',l2,'  m=',m2,'  explicit laplacian=',S2,'  Y_Laplacian_Y result=',S2_,S2__
      write (*,'(A,I2,A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'l=',l3,'  m=',m3,'  explicit laplacian=',S3,'  Y_Laplacian_Y result=',S3_,S3__
      write (*,'(A,I2,A,I2,A,I2,A,I2,A,E22.16,A,E22.16,E24.16)')  & 
        'l=',l,'  m=',m,'l=',l4,'  m=',m4,'  explicit laplacian=',S4,'  Y_Laplacian_Y result=',S4_,S4__
      
      ! test result
      if ( abs((S1-S1_)/S1) > delta .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
        print *,'Error: Laplacian differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Laplacian failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2_)/S2) > delta .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
        print *,'Error: Laplacian differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Laplacian failed'
        print *,''
        stop 1
      end if
      if ( abs((S3-S3_)/S3) > delta .and. ( abs(S3)>1.0d-10 .or. abs(S3_)>1.0d-10 ) ) then
        print *,'Error: Laplacian differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Laplacian failed'
        print *,''
        stop 1
      end if
      if ( abs((S4-S4_)/S4) > delta .and. ( abs(S4)>1.0d-10 .or. abs(S4_)>1.0d-10 ) ) then
        print *,'Error: Laplacian differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Laplacian failed'
        print *,''
        stop 1
      end if
      if ( abs((S1-S1__)/S1) > delta .and. ( abs(S1)>1.0d-10 .or. abs(S1__)>1.0d-10 ) ) then
        print *,'Error: Laplacian differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Laplacian failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2__)/S2) > delta .and. ( abs(S2)>1.0d-10 .or. abs(S2__)>1.0d-10 ) ) then
        print *,'Error: Laplacian differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Laplacian failed'
        print *,''
        stop 1
      end if
      if ( abs((S3-S3__)/S3) > delta .and. ( abs(S3)>1.0d-10 .or. abs(S3__)>1.0d-10 ) ) then
        print *,'Error: Laplacian differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Laplacian failed'
        print *,''
        stop 1
      end if
      if ( abs((S4-S4__)/S4) > delta .and. ( abs(S4)>1.0d-10 .or. abs(S4__)>1.0d-10 ) ) then
        print *,'Error: Laplacian differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Laplacian failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Laplacian passed up to lmax=',lmax
  print *,''
  
  ! check consistency for two and three centers routine
  S1=Y_Laplacian_C(b,r2,l2,m2,a,r3,nx3,ny3,nz3)
  S2=C_Laplacian_Y(a,r3,nx3,ny3,nz3,b,r2,l2,m2)
  ! test result
  if ( abs((S1-S2)/S1) > 1.0d-8 ) then
    print *,S1,S2
    print *,'Error: <YC> Laplacian failed'
    print *,''
    print *,'test_Laplacian failed'
    print *,''
    stop 1
  end if
  
end program 

