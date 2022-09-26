program test_R2Moment
  
  implicit none
  
  ! local variables
  integer, parameter      :: nmax = 3
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
  integer, parameter      :: l4   = 2
  integer, parameter      :: m4   =-1
  integer, parameter      :: lmax = 4
  real(kind=8), parameter :: r1(3) = (/ 0.2d0,-0.5d0, 0.3d0 /)
  real(kind=8), parameter :: r2(3) = (/ 0.2d0, 0.5d0, 0.1d0 /)
  real(kind=8), parameter :: r3(3) = (/ 0.1d0, 0.0d0, 0.5d0 /)
  real(kind=8), parameter :: r4(3) = (/-0.2d0,-0.1d0, 0.4d0 /)
  integer      :: i
  integer      :: j
  integer      :: k
  integer      :: l
  integer      :: m
  integer      :: nx
  integer      :: ny
  integer      :: nz
  real(kind=8) :: x
  real(kind=8) :: y
  real(kind=8) :: z
  real(kind=8) :: tmp
  real(kind=8) :: S1,S2,S3,S4
  real(kind=8) :: S1_,S2_,S3_,S4_
  real(kind=8) :: C_R2Moment_C
  real(kind=8) :: Y_R2Moment_Y
  real(kind=8) :: Y_Value
  
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        if ( nx+ny+nz.le.3 ) then
        
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
              tmp= ( (x-r1(1))**2 + (y-r1(2))**2 + (z-r1(3))**2 )
              S1 = S1 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r1(1))**nx1 * (y-r1(2))**ny1 * (z-r1(3))**nz1 &
                      * exp( -b * tmp ) * tmp**0
              tmp= ( (x-r2(1))**2 + (y-r2(2))**2 + (z-r2(3))**2 )
              S2 = S2 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r2(1))**nx2 * (y-r2(2))**ny2 * (z-r2(3))**nz2 &
                      * exp( -b * tmp ) * tmp**1
              tmp= ( (x-r3(1))**2 + (y-r3(2))**2 + (z-r3(3))**2 )
              S3 = S3 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r3(1))**nx3 * (y-r3(2))**ny3 * (z-r3(3))**nz3 &
                      * exp( -b * tmp ) * tmp**2
              tmp= ( (x-r4(1))**2 + (y-r4(2))**2 + (z-r4(3))**2 )
              S4 = S4 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) ) &
                      * (x-r4(1))**nx4 * (y-r4(2))**ny4 * (z-r4(3))**nz4 &
                      * exp( -b * tmp ) * tmp**3
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        S2=S2*8.0d0*lg**3/ng**3
        S3=S3*8.0d0*lg**3/ng**3
        S4=S4*8.0d0*lg**3/ng**3
        
        ! get corresponding overlap
        S1_= C_R2Moment_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r1,nx1,ny1,nz1,0)
        S2_= C_R2Moment_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r2,nx2,ny2,nz2,2)
        S3_= C_R2Moment_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r3,nx3,ny3,nz3,4)
        S4_= C_R2Moment_C(a,(/0.0d0,0.0d0,0.0d0/),nx,ny,nz,b,r4,nx4,ny4,nz4,6)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit R2Moment=',S1,'  C_R2Moment_C result=',S1_
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit R2Moment=',S2,'  C_R2Moment_C result=',S2_
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit R2Moment=',S3,'  C_R2Moment_C result=',S3_
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit R2Moment=',S4,'  C_R2Moment_C result=',S4_
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
          print *,'Error: R2Moment differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_R2Moment failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
          print *,'Error: R2Moment differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_R2Moment failed'
          print *,''
          stop 1
        end if
        if ( abs((S3-S3_)/S3) > 1.0d-8 .and. ( abs(S3)>1.0d-10 .or. abs(S3_)>1.0d-10 ) ) then
          print *,'Error: R2Moment differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_R2Moment failed'
          print *,''
          stop 1
        end if
        if ( abs((S4-S4_)/S4) > 1.0d-8 .and. ( abs(S4)>1.0d-10 .or. abs(S4_)>1.0d-10 ) ) then
          print *,'Error: R2Moment differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_R2Moment failed'
          print *,''
          stop 1
        end if
        
        end if
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_R2Moment passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
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
            
            S1=S1+Y_Value(b,l1,m1,(/x-r1(1),y-r1(2),z-r1(3)/))*tmp &
                 *( (x-r1(1))**2 + (y-r1(2))**2 + (z-r1(3))**2 )**0
            S2=S2+Y_Value(b,l2,m2,(/x-r2(1),y-r2(2),z-r2(3)/))*tmp &
                 *( (x-r2(1))**2 + (y-r2(2))**2 + (z-r2(3))**2 )**1
            S3=S3+Y_Value(b,l3,m3,(/x-r3(1),y-r3(2),z-r3(3)/))*tmp &
                 *( (x-r3(1))**2 + (y-r3(2))**2 + (z-r3(3))**2 )**2
            S4=S4+Y_Value(b,l4,m4,(/x-r4(1),y-r4(2),z-r4(3)/))*tmp &
                 *( (x-r4(1))**2 + (y-r4(2))**2 + (z-r4(3))**2 )**3
            
          end do
        end do
      end do
      
      S1=S1*8.0d0*lg**3/ng**3
      S2=S2*8.0d0*lg**3/ng**3
      S3=S3*8.0d0*lg**3/ng**3
      S4=S4*8.0d0*lg**3/ng**3
      
      ! get corresponding Overlap
      S1_= Y_R2Moment_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r1,l1,m1,0)
      S2_= Y_R2Moment_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r2,l2,m2,2)
      S3_= Y_R2Moment_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r3,l3,m3,4)
      S4_= Y_R2Moment_Y(a,(/0.0d0,0.0d0,0.0d0/),l,m,b,r4,l4,m4,6)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit R2Moment=',S1,'  Y_R2Moment_Y result=',S1_
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit R2Moment=',S2,'  Y_R2Moment_Y result=',S2_
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit R2Moment=',S3,'  Y_R2Moment_Y result=',S3_
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit R2Moment=',S4,'  Y_R2Moment_Y result=',S4_
      
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-8 .and. ( abs(S1)>1.0d-10 .or. abs(S1_)>1.0d-10 ) ) then
        print *,'Error: R2Moment differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_R2Moment failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-S2_)/S2) > 1.0d-8 .and. ( abs(S2)>1.0d-10 .or. abs(S2_)>1.0d-10 ) ) then
        print *,'Error: R2Moment differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_R2Moment failed'
        print *,''
        stop 1
      end if
      if ( abs((S3-S3_)/S3) > 1.0d-8 .and. ( abs(S3)>1.0d-10 .or. abs(S3_)>1.0d-10 ) ) then
        print *,'Error: R2Moment differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_R2Moment failed'
        print *,''
        stop 1
      end if
      if ( abs((S4-S4_)/S4) > 1.0d-8 .and. ( abs(S4)>1.0d-10 .or. abs(S4_)>1.0d-10 ) ) then
        print *,'Error: R2Moment differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_R2Moment failed'
        print *,''
        stop 1
      end if
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_R2Moment passed up to lmax=',lmax
  print *,''
  
end program 

