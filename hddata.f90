!=================================================================================
MODULE HdDATA
  IMPLICIT NONE
  INTEGER :: ncoord   !total number of internal coords
  INTEGER :: nstates  !number of electronic states
CONTAINS
!---------------------------------------------------------------------------------
SUBROUTINE EvaluateHd(dispgeom,nvibs,hmat,dhmat)
  use progdata, only: abpoint,natoms
  use DiabaticHamiltonian
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nvibs
  TYPE(abpoint),INTENT(IN) :: dispgeom
  DOUBLE PRECISION,DIMENSION(nstates,nstates),INTENT(OUT) :: hmat
  DOUBLE PRECISION,DIMENSION(nvibs,nstates,nstates),INTENT(OUT) :: dhmat
  real*8, allocatable :: q(:),t(:),dt(:,:)
  real*8, allocatable :: dhtmp(:,:,:),dtmp(:)
  integer :: i,j

  allocate(q(Hd_intdim),t(NHdExpansionBasis))
  allocate(dt(Hd_intdim,NHdExpansionBasis))
  allocate(dhtmp(ncoord,nstates,nstates),dtmp(3*natoms))

  q(1:Hd_intdim)=dispgeom%igeom(1:Hd_intdim)

  do i=1,NHdExpansionBasis
    t(i)=ExpansionBasis(q,i)
    dt(:,i)=ExpansionBasisGradient(q,i)
  end do

  dhtmp=0.d0
  do j=1,nstates
    do i=j,nstates
      hmat(i,j)=dot_product(t,Hd_HdEC(i,j)%Array)
      dhtmp(1:Hd_intdim,i,j)=matmul(dt,Hd_HdEC(i,j)%Array)
      !fill up
      if(i.ne.j) then
        hmat(j,i)=hmat(i,j)
        dhtmp(:,j,i)=dhtmp(:,i,j)
      end if
    end do
  end do

  do i=1,nstates
     do j=i,nstates
       call dgemv('T',ncoord,3*natoms,1.d0,dispgeom%bmat,ncoord,&
                      dhtmp(:,i,j),1,0.d0,dtmp,1)
       dhmat(1:nvibs,i,j)=dtmp(1:nvibs)
       if(j.ne.i) dhmat(1:nvibs,j,i)=dhmat(1:nvibs,i,j)
     end do
  end do

  return
END SUBROUTINE EvaluateHd
!---------------------------------------------------------------------------------
END MODULE HdDATA
!=================================================================================
