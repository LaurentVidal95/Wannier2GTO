program test_R_1_norm
  
  use mod_XZY_power_to_ir
  use mod_R_1_norm
  
  implicit none
  
  ! local variables
  integer, parameter      :: nmax = 8
  integer, parameter      :: ng   = 100
  real(kind=8), parameter :: a    = 2.0d0
  real(kind=8), parameter :: lg   = 5.0d0
  integer      :: i
  integer      :: j
  integer      :: k
  integer      :: ir
  integer      :: nx=0
  integer      :: ny=2
  integer      :: nz=4
  integer      :: lmax
  real(kind=8) :: x
  real(kind=8) :: y
  real(kind=8) :: z
  real(kind=8) :: norm_1
  real(kind=8) :: norm_1_
  real(kind=8) :: c(2925)
  INTEGER :: count_rate
  INTEGER :: start_count,end_count

!  print *,'tic'
!  CALL SYSTEM_CLOCK(start_count, count_rate)
!  ir = ir_index(nx,ny,nz)
!  ! get corresponding 1_norm
!  c      = 0.0d0
!  c(ir)  = 1.0d0
!  lmax   = nx+ny+nz
!  do i=1,100000000
!    norm_1_= R_1_norm(a,c,lmax)
!  end do
!  CALL SYSTEM_CLOCK(end_count, count_rate)
!  print *,'tac',(1.0d0*(end_count-start_count))/count_rate
  
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        if ( mod(nx,2)==0 .and. mod(ny,2)==0 .and. mod(nz,2)==0 ) then
          
          ! explicit calculation
          norm_1=0.0d0
          do i=0,ng
            x=-lg+(2.0d0*lg*i)/ng
            do j=0,ng
              y=-lg+(2.0d0*lg*j)/ng
              do k=0,ng
                z=-lg+(2.0d0*lg*k)/ng
                norm_1 = norm_1 + x**nx * y**ny * z**nz * exp( -a * (x**2+y**2+z**2) )
              end do
            end do
          end do
          norm_1=norm_1*8.0d0*lg**3/ng**3
          
          ! get corresponding ir index
          ir = ir_index(nx,ny,nz)
          
          ! get corresponding 1_norm
          c      = 0.0d0
          c(ir)  = 1.0d0
          lmax   = nx+ny+nz
          norm_1_= R_1_norm(a,c,lmax)
          
          write (*,'(A,I2,A,I2,A,I2,A,I4,A,E22.16,A,E22.16)')  & 
            'nx=',nx,'  ny=',ny,'  nz=',nz,'  ir=',ir,'  explicit 1 norm=',norm_1,'  R_1_norm result=',norm_1_
          
          ! test result
          if ( abs(norm_1-norm_1_)/norm_1 > 1.0d-8 ) then
            print *,'Error: 1 norm differs from explicit calculation for case'
            print *,'       nx =',nx
            print *,'       ny =',ny
            print *,'       nz =',nz
            print *,''
            print *,'test_R_1_norm failed'
            print *,''
            stop 1
          end if
          
        end if
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_R_1_norm passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''
  
end program 

