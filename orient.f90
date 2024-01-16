!==============================================================================
subroutine setorigin(geom)
  implicit none
  real*8, intent(inout) :: geom(3,11)
  real*8 :: dx(3)
  integer :: i
  dx=geom(:,5)
  do i=1,11
    geom(:,i)=geom(:,i)-dx
  end do
  return
end
!==============================================================================
subroutine orientation(geom,R,info)
  implicit none
  real*8, intent(inout) :: geom(3,11)
  real*8, intent(out) :: R(3,3)
  integer, intent(out) :: info

  call setorigin(geom)
  call bfxyz(geom(:,5),geom(:,6),geom(:,3),R,info)
  R=transpose(R)
  geom=matmul(R,geom)

  return
end
!==============================================================================
subroutine bfxyz(a1,a2,a3,xyz,info)
  implicit none
  real*8, intent(in) :: a1(3),a2(3),a3(3)
  real*8, intent(out) :: xyz(3,3)
  integer, intent(out) :: info
  real*8 :: dx(3),norm

  info=0
  !x
  dx=a2-a1
  norm=sqrt(dot_product(dx,dx))
  if(norm .lt. 1.d-6) then
    write(*,"('Small norm of vector in determining new x axis: ',e15.6)") norm
    info=-1
  end if
  xyz(:,1)=dx/norm

  !z
  call calc_cross_product(a2-a1,a3-a1,dx)
  norm=sqrt(dot_product(dx,dx))
  if(norm .lt. 1.d-6) then
    write(*,"('Small norm of vector in determining new z axis: ',e15.6)") norm
    print*, 'Warning: atoms may be colinear or nearly colinear.'
    info=-2
  end if
  xyz(:,3)=dx/norm

  !y
  call calc_cross_product(xyz(:,3),xyz(:,1),dx)
  norm=sqrt(dot_product(dx,dx))
  xyz(:,2)=dx/norm

  return
end
!==============================================================================
subroutine calc_cross_product(a,b,c)
  implicit none
  real*8, intent(in) :: a(3),b(3)
  real*8, intent(out) :: c(3)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=-a(1)*b(3)+a(3)*b(1)
  c(3)=a(1)*b(2)-a(2)*b(1)
  return
end
!==============================================================================
subroutine cart2dist(natoms,x,r)
  implicit none
  integer, intent(in) :: natoms
  real*8, intent(in) :: x(3,natoms)
  real*8, intent(out) :: r(natoms*(natoms-1)/2)
  real*8 :: dx(3)
  integer :: i,j,k
  k=0
  do i=1,natoms-1
    do j=i+1,natoms
      k=k+1
      dx=x(:,i)-x(:,j)
      r(k)=sqrt(dot_product(dx,dx))
    end do
  end do
  return
end
!==============================================================================
