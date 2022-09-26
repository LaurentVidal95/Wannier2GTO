
module mod_XZY_power_to_ir


implicit none

  integer, dimension(455) :: nx_from_ir =(/0,1,0,0,2,1,1,0,0,0,3,2,2,1,1,1,0,0,0,0,4,3,3,2,2,2,1,1,1,1,0,0,0,0,0,5,4,4,3,3,3,2,2,&
                                            2,2,1,1,1,1,1,0,0,0,0,0,0,6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,7,6 &
,&
                                            6,5,5,5,4,4,4,4,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,8,7,7,6,6,6,5,5,5 &
,&
                                            5,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,9,8,8,7,7,7,6 &
,&
                                            6,6,6,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0 &
,&
                                            0,0,0,0,0,10,9,9,8,8,8,7,7,7,7,6,6,6,6,6,5,5,5,5,5,5,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,2,&
                                            2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,11,10,10,9,9,9,8,8,8,8,7,7,&
                                            7,7,7,6,6,6,6,6,6,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2 &
,&
                                            1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,12,11,11,10,10,10,9,9,9,9,8,8,8,8,8,7,7 &
,&
                                            7,7,7,7,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2 &
,&
                                            2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
  integer, dimension(455) :: ny_from_ir =(/0,0,1,0,0,1,0,2,1,0,0,1,0,2,1,0,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,0,1,0,2,1,0,3,2,&
                                            1,0,4,3,2,1,0,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,0,1 &
,&
                                            0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1 &
,&
                                            0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3 &
,&
                                            2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,9,8,7,6,5 &
,&
                                            4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7 &
,&
                                            6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,10,9,8,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,&
                                            0,5,4,3,2,1,0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,10,&
                                            9,8,7,6,5,4,3,2,1,0,11,10,9,8,7,6,5,4,3,2,1,0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1 &
,&
                                            0,6,5,4,3,2,1,0,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,10,9,8,7,6,5,4,&
                                            3,2,1,0,11,10,9,8,7,6,5,4,3,2,1,0,12,11,10,9,8,7,6,5,4,3,2,1,0/)
  integer, dimension(455) :: nz_from_ir =(/0,0,0,1,0,0,1,0,1,2,0,0,1,0,1,2,0,1,2,3,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,0,1,0,1,2,0,1,&
                                            2,3,0,1,2,3,4,0,1,2,3,4,5,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,0 &
,&
                                            1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,0,1,0,1,2,0,1,2 &
,&
                                            3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,0,1,0,1,2,0 &
,&
                                            1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,1,2,3,4 &
,&
                                            5,6,7,8,9,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1 &
,&
                                            2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,10,0,0,1,0,1,2,0,1,2,3,0,1,2,3,&
                                            4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,1 &
,&
                                            2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9,10,11,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,&
                                            5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7 &
,&
                                            8,9,10,0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7,8,9,10,11,12/)

private
public :: ir_index
public :: XYZ_power

contains

!!> return (n+k-1)!/n!/(k-1)!
!!!
recursive function C(n,k)
  
  ! inputs argument
  integer, intent(in) :: n
  integer, intent(in) :: k
  
  ! ouput
  integer             :: C
  
  ! local variables
  integer             :: i
  
  C=1
  do i=n+1,n+k-1
    C=C*i
  end do
  do i=2,k-1
    C=C/i
  end do
  
end function

!!> return the offset in the ir basis of the l=nx+ny+nz components
!!!
recursive function offset_l(order)
  
  implicit none
  
  ! inputs argument
  integer, intent(in) :: order
  
  ! ouput
  integer             :: offset_l
  
  ! local variables
  integer             :: i
  
  ! use precomputed values for low order
  select case (order)
    case (0)
      offset_l=0
    case (1)
      offset_l=1
    case (2)
      offset_l=4
    case (3)
      offset_l=10
    case (4)
      offset_l=20
    case (5)
      offset_l=35
    case (6)
      offset_l=56
    case (7)
      offset_l=84
    case (8)
      offset_l=120
    case (9)
      offset_l=165
    case (10)
      offset_l=220
    case (11)
      offset_l=286
    case (12)
      offset_l=364
    case default
      ! general formula
      offset_l=0
      do i=0,order-1
        offset_l=offset_l+C(i,3)
      end do
  end select
  
end function


!!> return the local offset in the ir basis at l=nx+ny+nz du to nx
!!!
recursive function offset_x(order)
  
  implicit none
  
  ! inputs argument
  integer, intent(in) :: order
  
  ! ouput
  integer             :: offset_x
  
  ! local variables
  integer             :: i
  
  ! use precomputed values for low order
  select case (order)
    case (0)
      offset_x=0
    case (1)
      offset_x=1
    case (2)
      offset_x=3
    case (3)
      offset_x=6
    case (4)
      offset_x=10
    case (5)
      offset_x=15
    case (6)
      offset_x=21
    case (7)
      offset_x=28
    case (8)
      offset_x=36
    case (9)
      offset_x=45
    case (10)
      offset_x=55
    case (11)
      offset_x=66
    case (12)
      offset_x=78
    case default
      ! general formula
      offset_x=0
      do i=0,order-1
        offset_x=offset_x+C(i,2)
      end do
  end select
  
end function

!!> return index in the ir basis corresponding to component x^nx * y^ny * z^nz
!!!
!!! offset_l = (nx+ny+nz)*(nx+ny+nz+1)/2
!!! ir = offset_l + (ny+nz)*(ny+nz+1)/2+ny+nz + 1
!!!
recursive function ir_index(nx,ny,nz)
  
  implicit none
  
  ! input arguments
  integer, intent(in) :: nx
  integer, intent(in) :: ny
  integer, intent(in) :: nz
  
  ! return value
  integer             :: ir_index
  
  ! get offset from total l
  ir_index=offset_l(nx+ny+nz)
  
  ! add offset from local x
  ir_index=ir_index+offset_x(ny+nz)
  
  ! add offset from local y
  ir_index=ir_index+nz
  
  ! add 1 to get fortran indexing convention
  ir_index=ir_index+1
  
end function

!!> return powers of component x^nx * y^ny * z^nz given by index ir in the R basis
!!!
recursive subroutine XYZ_power(ir,nx,ny,nz)
  
  implicit none
  
  ! input arguments
  integer, intent(in)  :: ir
  integer, intent(out) :: nx
  integer, intent(out) :: ny
  integer, intent(out) :: nz
  
  ! test input
  if ( ir>455 ) then
    print *,'Error: XYZ_power not implemented for ir>455'
    stop
  end if
  
  ! set values
  nx=nx_from_ir(ir)
  ny=ny_from_ir(ir)
  nz=nz_from_ir(ir)
  
end subroutine

end module
