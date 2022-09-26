program test_Overlap
  
  implicit none
  
  ! local variables
  integer     , parameter :: lmax = 6
  real(kind=8), parameter :: PI = 4.0d0*atan(1.0d0)
  real(kind=8) :: a
  real(kind=8) :: d2
  real(kind=8) :: r1(3)
  real(kind=8) :: r2(3)
  real(kind=8) :: r3(3)
  integer      :: i
  integer      :: j
  integer      :: l
  integer      :: m
  integer      :: nx
  integer      :: ny
  integer      :: nz
  integer      :: Rl_dim
  real(kind=8) :: theta
  real(kind=8) :: phi
  real(kind=8) :: psi
  real(kind=8) :: R(9)
  real(kind=8) :: Rl(28*28)
  real(kind=8) :: val_ref(28)
  real(kind=8) :: val_tst(28)
  real(kind=8) :: val_rot(28)
  real(kind=8) :: Y_Value
  
  do l=0,lmax
    
    ! choose a random a 
    call RANDOM_NUMBER(a)
    
    ! choose an arbitrary r1
    call RANDOM_NUMBER(r1)
    
    ! choose an arbitrary r2
    call RANDOM_NUMBER(r2)
    
    ! compute value at r2
    d2=sum((r2-r1)**2)
    
    ! choose an arbitrary rotation angle
    call RANDOM_NUMBER(theta)
    call RANDOM_NUMBER(phi)
    call RANDOM_NUMBER(psi)
    
    ! generate corresponding rotation matrix
    R(1+0*3) = cos(2*PI*phi)*cos(2*PI*theta)
    R(1+1*3) = cos(2*PI*psi)*sin(2*PI*theta)-sin(2*PI*phi)*sin(2*PI*psi)*cos(2*PI*theta)
    R(1+2*3) = sin(2*PI*psi)*sin(2*PI*theta)+sin(2*PI*phi)*cos(2*PI*psi)*cos(2*PI*theta)
    R(2+0*3) =-cos(2*PI*phi)*sin(2*PI*theta)
    R(2+1*3) = sin(2*PI*phi)*sin(2*PI*psi)*sin(2*PI*theta)+cos(2*PI*psi)*cos(2*PI*theta)
    R(2+2*3) = sin(2*PI*psi)*cos(2*PI*theta)-sin(2*PI*phi)*cos(2*PI*psi)*sin(2*PI*theta)
    R(3+0*3) =-sin(2*PI*phi)
    R(3+1*3) =-cos(2*PI*phi)*sin(2*PI*psi)
    R(3+2*3) = cos(2*PI*phi)*cos(2*PI*psi)
    
    ! compute r3 = rotated r2 point
    r3 = 0.0d0
    do i=1,3
      do j=1,3
        r3(i) = r3(i)+R(i+(j-1)*3)*(r2(j)-r1(j))
      end do
    end do
    r3=r3+r1
    
    ! get ref cubic values
    i=0
    val_ref=0.0d0
    do nx=l,0,-1
      do ny=l-nx,0,-1
        ! increment i
        i=i+1
        ! set nz
        nz=l-nx-ny
        ! get corresponding reference value
        val_ref(i) = exp(-d2*a) * (r2(1)-r1(1))**nx * (r2(2)-r1(2))**ny * (r2(3)-r1(3))**nz
      end do
    end do
    
    ! compute cubic values at r3
    i=0
    val_tst=0.0d0
    do nx=l,0,-1
      do ny=l-nx,0,-1
        ! increment i
        i=i+1
        ! set nz
        nz=l-nx-ny
        ! get corresponding reference value
        val_tst(i) = exp(-d2*a) * (r3(1)-r1(1))**nx * (r3(2)-r1(2))**ny * (r3(3)-r1(3))**nz
      end do
    end do
    
    ! compute orbital rotation matrix
    call R_Rotation_Matrix(l,R,Rl,Rl_dim)
    
    ! compute rotated orbital value
    val_rot=0.0d0
    do i=1,Rl_dim
      do j=1,Rl_dim
        val_rot(i) = val_rot(i)+Rl(i+(j-1)*Rl_dim)*val_tst(j)
      end do
    end do
    
    ! test result
    do i=1,Rl_dim
      if ( abs((val_ref(i)-val_rot(i))/(val_rot(i)+val_ref(i))) > 1.0d-7 ) then
        print *,'Error: Rotated cubic orbital differs from explicit calculation for case l=',l
        print *,''
        print *,val_ref(i)
        print *,''
        print *,val_rot(i)
        print *,''
        print *,'test_Transform failed'
        print *,''
        stop 1
      end if
    end do
    
    ! get ref solid values
    i=0
    val_ref=0.0d0
    do m=-l,l
      ! increment i
      i=i+1
      ! get corresponding reference value
      val_ref(i) = Y_Value(a,l,m,r2-r1)
    end do
    
    ! compute solid values at r3
    i=0
    val_tst=0.0d0
    do m=-l,l
      ! increment i
      i=i+1
      ! get corresponding reference value
      val_tst(i) = Y_Value(a,l,m,r3-r1)
    end do
    
    ! compute orbital rotation matrix
    call Y_Rotation_Matrix(l,R,Rl,Rl_dim)
    
    ! compute rotated orbital value
    val_rot=0.0d0
    do i=1,Rl_dim
      do j=1,Rl_dim
        val_rot(i) = val_rot(i)+Rl(i+(j-1)*Rl_dim)*val_tst(j)
      end do
    end do
    
    ! test result
    do i=1,Rl_dim
      if ( abs((val_ref(i)-val_rot(i))/(val_rot(i)+val_ref(i))) > 1.0d-7 ) then
        print *,'Error: Rotated solid orbital differs from explicit calculation for case l=',l
        print *,''
        print *,val_ref(i)
        print *,''
        print *,val_rot(i)
        print *,''
        print *,'test_Transform failed'
        print *,''
        stop 1
      end if
    end do
    
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Transform passed up to l=',lmax
  print *,''
  
end program 

