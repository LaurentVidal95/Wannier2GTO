





program test_Electrostatics
  
  use mod_R_1_norm
  
  implicit none
  
  ! local variables
  integer, parameter      :: n_source = 18
  integer, parameter      :: nmax = 3
  integer, parameter      :: ng   = 127
  real(kind=8), parameter :: a    = 1.0d16
  real(kind=8), parameter :: b    = 3.5d0
  real(kind=8), parameter :: d    = 2.0d0
  real(kind=8), parameter :: lg   = 7.0d0
  real(kind=8), parameter :: x0   = -12.0d0
  real(kind=8), parameter :: y0   = -13.5d0
  real(kind=8), parameter :: z0   =  15.0d0
  integer, parameter      :: nx2  = 2
  integer, parameter      :: ny2  = 0
  integer, parameter      :: nz2  = 1
  integer, parameter      :: l2   = 2
  integer, parameter      :: m2   =-1
  integer, parameter      :: lmax = 6
  real(kind=8), parameter :: r1(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  real(kind=8), parameter :: r2(3) = (/ 0.5d0,-0.3d0, 0.2d0 /)
  integer      :: i
  integer      :: j
  integer      :: k
  integer      :: l
  integer      :: m
  integer      :: l1
  integer      :: m1
  integer      :: nx
  integer      :: ny
  integer      :: nz
  real(kind=8) :: Ex
  real(kind=8) :: Ey
  real(kind=8) :: Ez
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
  real(kind=8) :: E_(3)
  real(kind=8) :: E__(3)
  real(kind=8) :: potential_from_C
  real(kind=8) :: potential_from_Y
  real(kind=8) :: YY_Coulomb_Ion
  real(kind=8) :: Y_Value
  real(kind=8) :: r_source(3,n_source)
  real(kind=8) :: wv_source(n_source)
  real(kind=8) :: wx_source(n_source)
  real(kind=8) :: wy_source(n_source)
  real(kind=8) :: wz_source(n_source)
  
  !!
  !! check for C routines
  !!
  do nx=0,12
    do ny=0,12-nx
      do nz=0,12-nx-ny
        
        ! get potential
        S1_= potential_from_C((/x0,y0,z0/),b,r1,nx,ny,nz)
        
        ! get field
        call field_from_C((/x0,y0,z0/),b,r1,nx,ny,nz,E_)
        
        ! try electrostatic routine
        call electrostatics_from_C((/x0,y0,z0/),b,r1,nx,ny,nz,S1__,E__)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  potential_from_C=',S1_,'  potential_from_C=',S1__
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  field_from_C    =',E_(1),'  field_from_C    =',E__(1)
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  field_from_C    =',E_(2),'  field_from_C    =',E__(2)
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  field_from_C    =',E_(3),'  field_from_C    =',E__(3)
        
        ! test result
        if ( abs((S1_-S1__)/S1_) > 1.0d-8 ) then
          print *,'Error: potential_from_C differs from electrostatic_from_C'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        
        ! test result
        do i=1,3
          if ( abs((E_(i)-E__(i))/E_(i)) > 1.0d-8 ) then
            print *,'Error: field_from_C differs from electrostatic_from_C'
            print *,'       nx =',nx
            print *,'       ny =',ny
            print *,'       nz =',nz
            print *,''
            print *,'test_Electrostatics failed'
            print *,''
            stop 1
          end if
        end do
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to lmax=',12
  print *,''
  
  !!
  !! check for Y routine
  !!
  do l=0,6
    do m=-l,l
      
      ! get corresponding interaction
      S1_= potential_from_Y((/x0,y0,z0/),b,r1,l,m)
      
      ! get field
      call field_from_Y((/x0,y0,z0/),b,r1,l,m,E_)
        
      ! try electrostatic routine
      call electrostatics_from_Y((/x0,y0,z0/),b,r1,l,m,S1__,E__)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  potential_from_Y=',S1_,'  potential_from_Y=',S1__
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  field_from_Y    =',E_(1),'  field_from_Y    =',E__(1)
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  field_from_Y    =',E_(2),'  field_from_Y    =',E__(2)
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  field_from_Y    =',E_(3),'  field_from_Y    =',E__(3)
      
      ! test result
      if ( abs((S1_-S1__)/S1_) > 1.0d-8 ) then
        print *,'Error: potential_from_Y differs from electrostatic_from_Y'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      
      ! test result
      do i=1,3
        if ( abs((E_(i)-E__(i))/E_(i)) > 1.0d-8 ) then
          print *,'Error: field_from_C differs from electrostatic_from_C'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
      end do
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to lmax=',12
  print *,''
  
  !!
  !! check for cumulative C electrostatic routines
  !!
  call random_number(r_source)
  r_source = 10.0d0*r_source-5.0d0
  call random_number(wv_source)
  call random_number(wx_source)
  call random_number(wy_source)
  call random_number(wz_source)
  do nx=0,12
    do ny=0,12-nx
      do nz=0,12-nx-ny
        
        ! compute values explicitely 
        S1_=0.0d0
        E_=0.0d0
        do i=1,n_source
          ! get potential
          S1__= potential_from_C(r_source(:,i),b,r1,nx,ny,nz)
          S1_ = S1_ + wv_source(i) * S1__ 
          ! get field
          call field_from_C(r_source(:,i),b,r1,nx,ny,nz,E__)
          E_(1) = E_(1) + wx_source(i) * E__(1) 
          E_(2) = E_(2) + wy_source(i) * E__(2) 
          E_(3) = E_(3) + wz_source(i) * E__(3) 
        end do
        
        ! try cumulative routine
        call cumulative_electrostatics_on_C(r_source, &
                                            wv_source, &
                                            wx_source, &
                                            wy_source, &
                                            wz_source, &
                                            n_source, &
                                            b,r1,nx,ny,nz, &
                                            S1__,E__(1),E__(2),E__(3))
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  cumulated potential_from_C=',S1_,'  cumulated potential_from_C=',S1__
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  cumulated field_from_C    =',E_(1),'  cumulated field_from_C    =',E__(1)
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  cumulated field_from_C    =',E_(2),'  cumulated field_from_C    =',E__(2)
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  cumulated field_from_C    =',E_(3),'  cumulated field_from_C    =',E__(3)
        
        ! test result
        if ( abs((S1_-S1__)/S1_) > 1.0d-8 ) then
          print *,'Error: cumulated potential_from_C differs from cumulative_electrostatic_on_C'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        
        ! test result
        do i=1,3
          if ( abs((E_(i)-E__(i))/E_(i)) > 1.0d-8 ) then
            print *,'Error: cumulated field_from_C differs from cumulative_electrostatic_from_C'
            print *,'       nx =',nx
            print *,'       ny =',ny
            print *,'       nz =',nz
            print *,''
            print *,'test_Electrostatics failed'
            print *,''
            stop 1
          end if
        end do
        
      end do
    end do
  end do
  
  !!
  !! check for cumulative Y electrostatic routines
  !!
  do l=0,6
    do m=-l,l
      
      ! compute values explicitely 
      S1_=0.0d0
      E_=0.0d0
      do i=1,n_source
        ! get potential
        S1__= potential_from_Y(r_source(:,i),b,r1,l,m)
        S1_ = S1_ + wv_source(i) * S1__ 
        ! get field
        call field_from_Y(r_source(:,i),b,r1,l,m,E__)
        E_(1) = E_(1) + wx_source(i) * E__(1) 
        E_(2) = E_(2) + wy_source(i) * E__(2) 
        E_(3) = E_(3) + wz_source(i) * E__(3) 
      end do
      
      ! try cumulative routine
      call cumulative_electrostatics_on_Y(r_source, &
                                          wv_source, &
                                          wx_source, &
                                          wy_source, &
                                          wz_source, &
                                          n_source, &
                                          b,r1,l,m, &
                                          S1__,E__(1),E__(2),E__(3))
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  cumulated potential_from_Y=',S1_,'  cumulated potential_from_Y=',S1__
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  cumulated field_from_Y    =',E_(1),'  cumulated field_from_Y    =',E__(1)
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  cumulated field_from_Y    =',E_(2),'  cumulated field_from_Y    =',E__(2)
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  cumulated field_from_Y    =',E_(3),'  cumulated field_from_Y    =',E__(3)
      
      ! test result
      if ( abs((S1_-S1__)/S1_) > 1.0d-8 ) then
        print *,'Error: cumulated potential_from_Y differs from cumulative_electrostatic_on_Y'
        print *,'       nx =',nx
        print *,'       ny =',ny
        print *,'       nz =',nz
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      
      ! test result
      do i=1,3
        if ( abs((E_(i)-E__(i))/E_(i)) > 1.0d-8 ) then
          print *,'Error: cumulated field_from_Y differs from cumulative_electrostatic_from_Y'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
      end do
      
    end do
  end do
  
  
  
  !!
  !! check for cumulative YY electrostatic routines
  !!
  do l=0,6
    do m=-l,l
      do l1=0,6
        do m1=-l1,l1
          
          ! compute values explicitely 
          S1_=0.0d0
          E_=0.0d0
          do i=1,n_source
            ! get potential
            S1__= YY_Coulomb_Ion(d,r1,l1,m1,b,r2,l,m,r_source(:,i))
            S1_ = S1_ + wv_source(i) * S1__ 
          end do
          
          ! try cumulative routine
          call cumulative_electrostatics_on_YY(r_source, &
                                               wv_source, &
                                               wx_source, &
                                               wy_source, &
                                               wz_source, &
                                               n_source, &
                                               b,r2 ,l ,m , &
                                               d,r1,l1,m1, &
                                               S1__,E__(1),E__(2),E__(3))
          
          write (*,'(A,I2,A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
            'l1=',l1,'  m1=',m1,'l2=',l,'  m2=',m,'  cumulated potential_from_YY=',S1_,'  cumulated potential_from_YY=',S1__
          
          ! test result
          if ( abs((S1_-S1__)/S1_) > 1.0d-8 ) then
            print *,'Error: cumulated potential_from_Y differs from cumulative_electrostatic_on_Y'
            print *,'       l1 =',l1
            print *,'       m1 =',m1
            print *,'       l2 =',l
            print *,'       m2 =',m
            print *,''
            print *,'test_Electrostatics failed'
            print *,''
            stop 1
          end if
          
        end do
      end do
    end do
  end do
  
  
  
  
  
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
              S1 = S1 + tmp * V
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        
        ! get potential
        S1_= potential_from_C((/x0,y0,z0/),b,r1,nx,ny,nz)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S1,'  potential_from_C=',S1_        
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-6 ) then
          print *,'Error: potential_from_C differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
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
            S1 = S1 + tmp * V
          end do
        end do
      end do
      S1=S1*8.0d0*lg**3/ng**3
      
      ! get corresponding interaction
      S1_= potential_from_Y((/x0,y0,z0/),b,r1,l,m)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S1,'  potential_from_Y=',S1_
        
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-5 ) then
        print *,'Error: potential_from_Y differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to lmax=',lmax
  print *,''
  
  !!
  !! check for C routines
  !!
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! explicit calculation
        S1=0.0d0
        S2=0.0d0
        S3=0.0d0
        do k=0,ng
          z=-lg+(2.0d0*lg*k)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do i=0,ng
              x=-lg+(2.0d0*lg*i)/ng
              Ex  = -(x-x0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
              Ey  = -(y-y0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
              Ez  = -(z-z0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
              tmp= x**nx * y**ny * z**nz * exp( -b * (x**2+y**2+z**2) )
              S1 = S1 + tmp * Ex
              S2 = S2 + tmp * Ey
              S3 = S3 + tmp * Ez
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        S2=S2*8.0d0*lg**3/ng**3
        S3=S3*8.0d0*lg**3/ng**3
        
        ! get field
        call field_from_C((/x0,y0,z0/),b,r1,nx,ny,nz,E_)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S1,'  field_from_C(1)=',E_(1)
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S2,'  field_from_C(2)=',E_(2)
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S3,'  field_from_C(3)=',E_(3)
        
        ! test result
        if ( abs((S1-E_(1))/S1) > 1.0d-6 ) then
          print *,'Error: field_from_C differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-E_(2))/S2) > 1.0d-6 ) then
          print *,'Error: field_from_C differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        if ( abs((S3-E_(3))/S3) > 1.0d-6 ) then
          print *,'Error: field_from_C differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
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
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            Ex  = -(x-x0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
            Ey  = -(y-y0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
            Ez  = -(z-z0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
            tmp=Y_Value(b,l,m,(/x,y,z/))
            S1 = S1 + tmp * Ex
            S2 = S2 + tmp * Ey
            S3 = S3 + tmp * Ez
          end do
        end do
      end do
      S1=S1*8.0d0*lg**3/ng**3
      S2=S2*8.0d0*lg**3/ng**3
      S3=S3*8.0d0*lg**3/ng**3
      
      ! get field
      call field_from_Y((/x0,y0,z0/),b,r1,l,m,E_)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S1,'  field_from_Y(1)=',E_(1)
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S2,'  field_from_Y(2)=',E_(2)
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S3,'  field_from_Y(3)=',E_(3)
      
      ! test result
      if ( abs((S1-E_(1))/S1) > 1.0d-5 ) then
        print *,'Error: field_from_Y differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-E_(2))/S2) > 1.0d-5 ) then
        print *,'Error: field_from_Y differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      if ( abs((S3-E_(3))/S3) > 1.0d-5 ) then
        print *,'Error: field_from_Y differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to lmax=',lmax
  print *,''
  
end program

