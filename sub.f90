!================================================================================
subroutine idx(i,j,nstates,id)
  !get array index for the matrix element (i,j)
  !in the order of upper triangle matrix
  implicit none
  integer, intent(in) :: i,j,nstates
  integer, intent(out) :: id
  integer :: ii,jj,k

  ii=min(i,j)
  jj=max(i,j)

  id=0
  do k=1,ii-1
    id=id+nstates-k+1
  end do
  id=id+jj-ii+1

  return
end
!=================================================================================
subroutine FLUnit(fid)
      !get an available UNIT index for input/output
      implicit none
      integer, intent(out) :: fid
      integer :: i
      logical :: unitex,unitop
      fid=0
      do i=15,99999
            inquire(UNIT=i,EXIST=unitex,OPENED=unitop)
            if(unitex .and. .not. unitop) then
                  fid=i
                  exit
            end if
      end do
      if(fid.eq.0) stop "FLUnit: failed to find an available unit."
      return
end subroutine FLUnit
!=================================================================================
!---------------bubble sort the elements of an integer array-------------------------
SUBROUTINE bsort(nelements,A,ordering)
 implicit none
 integer,intent(IN)                         ::  nelements
 integer,dimension(nelements),intent(INOUT) ::  A
 integer,intent(IN)                         ::  ordering
 integer      :: i,j,temp

 do i=nelements-1,1,-1
  do j=1,i
    if(A(j)*ordering .gt. A(j+1)*ordering) then
      temp=A(j)
      A(j)=A(j+1)
      A(j+1)=temp
    end if
  end do
 end do

 return
END SUBROUTINE bsort
!=================================================================================
!M x N real matrix A with full row or column rank
!return the transpose of its generalized inverse
subroutine dGeneralizedInverseTranspose(A, M, N)
  integer,intent(in)::M,N
  real*8,dimension(M,N),intent(inout)::A
  real*8,dimension(M,M)::AAT
  integer :: info
  AAT=matmul(A,transpose(A))
  call dpotrf('L',M,AAT,M,info)
  if(info.ne.0) stop 'info.ne.0 dpotrf in dGeneralizedInverse!'
  call dpotri('L',M,AAT,M,info)
  if(info.ne.0) stop 'info.ne.0 dpotri in dGeneralizedInverse!'
  call syL2U(AAT,M)
  A=matmul(AAT,A)
  return
end subroutine dGeneralizedInverseTranspose
!=================================================================================
!N order matrix A, strictly upper triangle is blank,
!copy strictly lower triangle elements to strictly upper triangle
subroutine syL2U(A, N)
  integer,intent(in)::N
  real*8,dimension(N,N),intent(inout)::A
  integer::i,j
  forall(i=1:N-1,j=2:N,i<j)
    A(i,j)=A(j,i)
  end forall
  return
end subroutine syL2U
!=================================================================================
subroutine orthgh0(natoms,g,h,gort,hort)
  !Yarkony orthogonalization JCP 112, 2111 (2000)
  implicit none
  integer, intent(in) :: natoms
  real*8, intent(in) :: g(3,natoms),h(3,natoms)
  real*8, intent(out) :: gort(3,natoms),hort(3,natoms)
  real*8 :: g2,h2,gh,beta

  call scalar(natoms,g,h,gh)
  if(abs(gh).lt.1.d-10) then
    print*,"gh very small. No rotation is performed."
    return
  end if

  call scalar(natoms,g,g,g2)
  call scalar(natoms,h,h,h2)
  beta = 0.5d0*datan(2.d0*gh/(h2-g2))
  gort =  g*dcos(beta) - h*dsin(beta)
  hort =  g*dsin(beta) + h*dcos(beta)
  return

end subroutine
!===============================================================================
subroutine scalar(natoms,v1,v2,q)
  implicit none
  integer, intent(in) :: natoms
  real*8, intent(in) :: v1(3,natoms),v2(3,natoms)
  real*8, intent(out) :: q
  integer :: i
  q=0.d0
  do i=1,natoms
    q=q+dot_product(v1(1:3,i),v2(1:3,i))
  end do
  return
end subroutine
!===============================================================================
subroutine bdchk(t)
  implicit none
  real*8, intent(inout) :: t
  if(abs(t) .gt. 1.d0) then
    if(t.gt.0.d0) then
      t=1.d0
    else
      t=-1.d0
    end if
  end if
  return
end
!===============================================================================
