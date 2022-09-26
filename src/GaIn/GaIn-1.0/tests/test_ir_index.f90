program test_ir_index
  
  use mod_XZY_power_to_ir
  
  implicit none
  
  integer :: nx,ny,nz,l
  integer :: nx_,ny_,nz_
  integer :: ir
  integer, parameter :: lmax=24
  integer, parameter :: lmax_=12
  
  ir=0
  do l=0,lmax
    do nx=l,0,-1
      do ny=l-nx,0,-1
        ! increment ir
        ir=ir+1
        ! compute nz
        nz=l-nx-ny
        ! test ir calculation 
        if ( ir .ne. ir_index(nx,ny,nz) ) then
          print *,'Error: test_ir_index failed for the following entries'
          print *,'       nx=',nx
          print *,'       ny=',ny
          print *,'       nz=',nz
          print *,'       ir=',ir
          print *,'       ir_index returned value ',ir_index(nx,ny,nz)
          stop 1
        end if
      end do
    end do
  end do
  
  ir=0
  do l=0,lmax_
    do nx=l,0,-1
      do ny=l-nx,0,-1
        ! increment ir
        ir=ir+1
        ! compute nz
        nz=l-nx-ny
        ! test inverse ir calculation 
        call XYZ_power(ir,nx_,ny_,nz_)
        if ( nx .ne. nx_ .or. ny .ne. ny_ .or. nz .ne. nz_ ) then
          print *,'Error: test_ir_index failed for the following entries'
          print *,'       nx=',nx
          print *,'       ny=',ny
          print *,'       nz=',nz
          print *,'       ir=',ir
          print *,'       XYZ_power returned values ',nx_,ny_,nz_
          stop 1
        end if
      end do
    end do
  end do

  
  ! inform 
  print *,''
  write (*,'(A,I2)') 'test_ir_index passed up to order lmax=',lmax
  print *,''
  
end program
