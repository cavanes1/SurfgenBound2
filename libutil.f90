!----------------------------------------------------------------------------------
! Partition states at a certain geometry into groups that are energetically close.
! Group states into minimal sets so that the energy difference between any two states
! from different groups is greater than de. Input energies are supposed to be sorted.
! Energy groups are only generated for the range state1 to state2
SUBROUTINE genEnerGroups(pt,de)
  use progdata, only: abpoint
  IMPLICIT NONE
  type(abpoint),INTENT(INOUT)     :: pt
  DOUBLE PRECISION,INTENT(IN)     :: de
  integer                         :: i,ngrp,cnt
  integer,dimension(pt%ub,2)      :: deg_table

  ngrp=0
  cnt=1
  do i=pt%lb+1,pt%ub
    if((pt%energy(i,i)-pt%energy(i-1,i-1))>de)then
      if(cnt>1)then
        ngrp=ngrp+1
        deg_table(ngrp,2)=cnt
        deg_table(ngrp,1)=i-cnt
      end if
      cnt=0
    end if
    cnt=cnt+1
  end do
  if(cnt>1)then
    ngrp=ngrp+1
    deg_table(ngrp,2)=cnt
    deg_table(ngrp,1)=pt%ub-cnt+1
  end if

  if(allocated(pt%deg_groups)) deallocate(pt%deg_groups)
  allocate(pt%deg_groups(ngrp,2))
  pt%ndeggrp=ngrp
  pt%deg_groups=deg_table(1:ngrp,:)

  return
END SUBROUTINE
!--------------------------------------------------------------------------
!Obtain the rotation angle beta that orthogonalize g and h
SUBROUTINE orthgh(n,g,h,beta)
  IMPLICIT NONE
  INTEGER,intent(IN)          :: n
  DOUBLE PRECISION,intent(IN) :: g(n),h(n)
  DOUBLE PRECISION,intent(OUT):: beta
  double precision :: denom,num
  num=dot_product(g,h)*2.d0
  denom=-dot_product(g,g)+dot_product(h,h)
  if(abs(denom)<1D-10)then
      PRINT *,"DENORM VERY SMALL.  NO ROTATION IS PERFORMED"
      beta=0D0
  else
      beta=0.25d0*atan(num/denom)
  end if
  return
END SUBROUTINE orthgh
!---------------------------------------------------------------------
! Chop() sets the value of the variable to 0 if it is smaller
! than a threshold
SUBROUTINE chop(dpvar,threshold)
  implicit none
  double precision, intent(INOUT) :: dpvar
  double precision, intent(in) :: threshold

  if(abs(dpvar)<abs(threshold)) dpvar=0.d0

  return
END SUBROUTINE chop
!---------------------------------------------------------------------
!---------------sort the elements of an array-------------------------
SUBROUTINE sort(nelements,data,ordering)
 implicit none
 integer,intent(IN)                         ::  nelements
 integer,dimension(nelements),intent(INOUT) ::  data
 integer,intent(IN)                         ::  ordering

 integer      :: i,j,temp

 do i=nelements-1,1,-1
  do j=i,nelements-1
    if(data(j)*ordering>data(j+1)*ordering)then
      temp=data(j)
      data(j)=data(j+1)
      data(j+1)=temp
    end if
  end do
 end do

 return
END SUBROUTINE sort
!-------------------------------------------------------------------
! SCHDMIT() Performs Schdmit orthorgonalization for rows of matrix A
!-------------------------------------------------------------------
! Each row in A will be orthogonalized to all row vectors above it.
! After removal, if the norm of residue vector is smaller than $(tol),
! it will be permuted to the bottom of A and considered 0. Otherwise,
! The residue vector will be normalized
!-------------------------------------------------------------------
! m,n       (input) INTEGER
!           Number of rows/columns
! A         (input/output) DOUBLE PRECISION,dimension(m,n)
!           On entry, rows of A contain vectors to be orthogonalized
!           On exit, first $(rank) rows of A are orthonormal and the
!           other rows are zero
! d         (input) INTEGER
!           Dimension of vector.  Only first $(d) columns will be 
!           used when calculating inner products and included in 
!           orthogonalization.  Other columns will simply be carried.  
! tol       (input) DOUBLE PRECISION
!           Tolerance for the norm of residue to be considered zero.
! rank      (output) INTEGER
!           Number of orthonormal vectors contained in A
! nstart    (input) INTEGER
!           The number of row to start orthogonalization.  The rows
!           above nstart will be supposed to be orthogonal and will
!           only be normalized.
SUBROUTINE Schdmit(m,n,A,d,tol,rank,nstart)
  IMPLICIT NONE
  INTEGER,intent(IN)                           :: m,n,d,nstart
  DOUBLE PRECISION,dimension(m,n),intent(INOUT):: A
  DOUBLE PRECISION,intent(IN)                  :: tol
  INTEGER,intent(OUT)                          :: rank
  integer           :: i
  double precision  :: p(m),nrm
  double precision,external :: DNRM2
  rank=m
  do i=1,nstart-1
    A(i,:)=A(i,:)/DNRM2(d,A(i,1:d),int(1))
  end do
  i=nstart
  do while(i<=rank)
    !Compute the inner product with all row vectors above
    if(i>1)then
      CALL DGEMV('N',i-1,d,dble(1),A,m,A(i,1),m,dble(0),p,int(1))
      CALL DGEMV('T',i-1,n,dble(-1),A,m,p,int(1),dble(1),A(i,1),m)
    end if!(i>1)
    nrm=DNRM2(d,A(i,1),m)
    if(nrm<tol)then
      CALL DSWAP(n,A(i,1),m,A(rank,1),m)
      rank=rank-1
    else!if(nrm<tol)
      A(i,:)=A(i,:)/nrm
      i=i+1
    end if!(nrm<tol)
  end do!while(i<m)

END SUBROUTINE Schdmit
!----------------------------------------------------------------------------
!DSYINV inverts a symmetric double precision matrix
SUBROUTINE DSYINV(UPLO,N,A,LDA,INVA,LDI,TOL,NNEG,NZERO,ENFPD)
  IMPLICIT NONE
  CHARACTER(1),intent(IN)     :: UPLO
  INTEGER,intent(IN)          :: N,LDA,LDI
  DOUBLE PRECISION,intent(IN) :: A(LDA,N),TOL
  DOUBLE PRECISION,intent(OUT):: INVA(LDI,N)
  LOGICAL,intent(IN)          :: ENFPD
  INTEGER,INTENT(OUT)         :: NNEG,NZERO
  double precision,dimension(:),allocatable     :: WORK
  integer,dimension(:),allocatable              :: IWORK
  integer          :: LWORK,LIWORK,ISUPPZ(2*N),INFO,itst(1),nev,j
  double precision :: dtst(1),eval(N),evec(N,N)

  CALL DSYEVR('V','A',UPLO,n,A,LDA,dble(0),dble(0),0,0,TOL,nev, &
         eval,evec,n,ISUPPZ,dtst,int(-1),itst,int(-1), INFO )
  if(INFO/=0)STOP"DSYINV: Workspace query failed."
  LWORK  = int(dtst(1))
  LIWORK = itst(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  CALL DSYEVR('V','A',UPLO,n,A,LDA,dble(0),dble(0),0,0,TOL,nev, &
             eval,evec,n,ISUPPZ,WORK,LWORK,IWORK,LIWORK, INFO )
  if(INFO/=0)STOP"DSYINV: eigenvalue decomposition failed."
  nneg=0
  nzero=n-nev
  invA=dble(0)
  do j=1,nev
    if(abs(eval(j))>1D-4)then
      if(eval(j)<0)nneg=nneg+1
      if(eval(j)<0.and. ENFPD)eval(j)=dble(1)
      CALL DSYR(UPLO,n,1/eval(j),evec(:,j),int(1),invA,LDI)
    else
      nzero=nzero+1
    end if
  end do
  deallocate(WORK)
  deallocate(IWORK)

END SUBROUTINE
!-----------------------------------------------------------------------------
! Rotate degenerate ab initio adiabatic data so that g and h vectors between any
! pairs of two degenerate states are orthogonalized.  The rotation is made so that
! total norm of off diagonal gradients (h vectors) is minimized.
! This subroutine is used as interface for external programs. 
!dh             : gradients in the current representation
!h              : hamiltonian in the current representation
!lstate,ustate  : range of degenerate states
!nvibs          : vibrational degrees of freedom
!
!maxiter: maximum number of jacobi rotations
!toler  : convergence tolerance for rotation angle beta
!hasGrad: specifies if gradient/couple data is available for a block
![method]
! This subroutine tries orthogonalize all g-h vectors by an algorithm similar to 
! the Jacobi's method for eigenvalues.  The Jacobi transformations here extremize 
! the norm of the corresponding coupling block rather than eliminating them.   
SUBROUTINE OrthogonalizeGH(h,dh,lstate,ustate,nvibs,maxiter,toler)
  USE hddata, only: nstates
  USE progdata, only: abpoint
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nvibs,nstates,nstates) :: dh
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nstates,nstates)       :: h
  INTEGER,INTENT(IN)                                    :: maxiter,lstate,ustate,nvibs
  DOUBLE PRECISION,INTENT(IN)                           :: toler

  integer           ::  i,j,iter
  integer           ::  mi,mj  !location of maximum rotation
  double precision, dimension(nvibs,nstates,nstates)    :: newgrad
  double precision, dimension(nstates,nstates)          :: newener 
  double precision, dimension(nstates,nstates)          :: beta   !required rotation 
  double precision           :: max_b,t, c,s

  beta=dble(0)  
  max_b=-1
  iter=0
  !build normalized g and h vectors from Hd and compute rotation angles the first time
  do i=lstate,ustate 
      do j=lstate,i-1
        beta(i,j) = getBeta(i,j)
        if(abs(beta(i,j))>max_b)then
          max_b=abs(beta(i,j))
          mi=i
          mj=j
        end if!(abs(beta(i,j))>max_b)
      end do!j=lstate,i-1
    end do!i=lstate,ustate
    do while(iter<maxiter.and.max_b>toler)
      iter=iter+1
      t=beta(mi,mj)
      c = cos(t)
      s = sin(t)
! Gnew = J^T.G.J.   Gnew_ij=Jki*Gkl*Jlj
      newgrad = dh 
      newener = h
      do i=1,nstates
        newgrad(:,i,mi)=dh(:,i,mi)*c-dh(:,i,mj)*s
        newgrad(:,i,mj)=dh(:,i,mj)*c+dh(:,i,mi)*s
        newener(i,mi)=h(i,mi)*c-h(i,mj)*s
        newener(i,mj)=h(i,mj)*c+h(i,mi)*s
      end do
      dh     = newgrad
      h      = newener
      do i=1,nstates
        dh(:,mi,i)=newgrad(:,mi,i)*c-newgrad(:,mj,i)*s
        dh(:,mj,i)=newgrad(:,mj,i)*c+newgrad(:,mi,i)*s 
        h(mi,i)=newener(mi,i)*c-newener(mj,i)*s
        h(mj,i)=newener(mj,i)*c+newener(mi,i)*s 
      end do
      !update rotation angles
      do i=mj+1,ustate
        beta(i,mj)=getBeta(i,mj)
      end do
      do j=lstate,mj-1
        beta(mi,j)=getBeta(mi,j)
      end do
      do j=mj+1,mi-1
        beta(mi,j)=getBeta(mi,j)
      end do
      max_b=-1
      do i=lstate,ustate
        do j=lstate,i-1
          if(abs(beta(i,j))>max_b)then
            max_b=abs(beta(i,j))
            mi=i
            mj=j
          end if!(abs(beta(i,j))>max_b)
        end do!j=lstate,i-1
      end do!i=lstate,ustate 
    end do!while(iter<maxiter.and.max_b>toler)do
CONTAINS
  FUNCTION getBeta(i,j) RESULT(beta)
    USE progdata, only : abpoint,pi
    IMPLICIT NONE
    INTEGER,INTENT(IN)         :: i,j
    DOUBLE PRECISION           :: beta

    double precision, dimension(nvibs)                :: g,h

    g=(dh(:,i,i)-dh(:,j,j))/2
    h=dh(:,i,j)
    CALL orthgh(nvibs,g,h,beta)
  END FUNCTION getBeta
END SUBROUTINE OrthogonalizeGH 
!------------------------------------------------------------------------------
! Rotate degenerate ab initio adiabatic data so that g and h vectors between any
! pairs of two degenerate states are orthogonalized.  The rotation is made so that
! total norm of off diagonal gradients (h vectors) is minimized.
!pt     : ab initio data of the point to perform transformation 
!maxiter: maximum number of jacobi rotations
!toler  : convergence tolerance for rotation angle beta
!hasGrad: specifies if gradient/couple data is available for a block
![method]
! This subroutine tries orthogonalize all g-h vectors by an algorithm similar to 
! the Jacobi's method for eigenvalues.  The Jacobi transformations here extremize 
! the norm of the corresponding coupling block rather than eliminating them.   
SUBROUTINE OrthGH_ab(pt,maxiter,toler,hasGrad,needall)
  USE hddata, only: nstates
  USE progdata, only: abpoint,printlvl
  IMPLICIT NONE
  TYPE(abpoint),INTENT(INOUT)                           :: pt
  INTEGER,INTENT(IN)                                    :: maxiter
  LOGICAL,DIMENSION(nstates,nstates),INTENT(in)         :: hasGrad
  DOUBLE PRECISION,INTENT(IN)                           :: toler
  LOGICAL,INTENT(in)                                    :: needall

  integer :: igrp,i,j,iter,ldeg,udeg
  integer :: mi,mj  !location of maximum rotation
  double precision, dimension(pt%nvibs,nstates,nstates) :: newgrad
  double precision, dimension(nstates,nstates)          :: newener 
  double precision, dimension(nstates,nstates)          :: beta   !required rotation 
  double precision :: max_b,t, c,s
  !allowedRot stores the infomation whether rotation between two specific states 
  !will not cause a gradient data to mix into a gradient that has not data
  LOGICAL,dimension(nstates,nstates) :: allowedRot

  beta=dble(0)  
  if(pt%ndeggrp>0.and.printlvl>1) &
    print *,"     Transforming ab initio data for point ",pt%id
  if(printlvl>2.and. pt%ndeggrp>0)then
    print *,"      Couplings and gradients before transformation" 
    do i=1,nstates
      do j=1,i
        print 1002,i,j
        if(hasGrad(i,j))then
            print 1003,pt%grads(:pt%nvibs,i,j)
        else
            print "(10X,A)","Data not available"
        end if
      end do
    end do!i=1,nstates
  end if!(printlvl>2)
grploop: do igrp=1,pt%ndeggrp
    !in case there are multiple groups of degeneracies, loops over all of them       
    max_b=-1.d0
    ldeg=pt%deg_groups(igrp,1)
    udeg=pt%deg_groups(igrp,2)+pt%deg_groups(igrp,1)-1
    iter=0
    !check if rotations among these states are allowed.
    if(all(hasGrad)) then
      allowedRot(ldeg:udeg,ldeg:udeg)=.true.
      do i=ldeg,udeg
        allowedRot(i,i)=.false.
      end do
    else
      allowedRot(ldeg:udeg,ldeg:udeg)=.false.
      if(printlvl>2)write (*,"(6X,A)",advance='no') "Allowed rotations:"
      do i=ldeg,udeg
          do j=ldeg,i-1
            !both gradients and coupling between them has to be present
            !for them to be rotatable
            if(hasGrad(i,i).and.hasGrad(i,j).and.hasGrad(j,j)) then
              !the available gradients list has to be identical for them
              if(all(hasGrad(i,:).eqv.hasGrad(j,:)))then
                  if(printlvl>2)write (*,"(' (',I1,',',I1,') ')",advance='no') I,J
                  allowedRot(i,j)=.true.
                  allowedRot(j,i)=.true.
              end if
            end if!hasGrad ii, ij, jj
          end do!j
      end do!i=ldeg,udeg
      if(.not.any(allowedRot(ldeg:udeg,ldeg:udeg)))then
        print *," none"
        print "(4X,A)","Missing data forbit any rotations.  Skipping transformation of degenerate group."
        cycle grploop
      else
        print *,""
      end if
      if(needall)then
        do i=ldeg,udeg
          do j=ldeg,i-1
            if(.not.allowedRot(i,j))then
              print *,"  Not all rotations are included.  Skipping block."
              cycle grploop
            end if
          end do
        end do
      end if!need all
    end if
    print *,"building g and h vectors"

    !build normalized g and h vectors from Hd
    !and compute rotation angles the first time
    do i=ldeg,udeg 
      do j=ldeg,i-1

        print *,"i,j=",i,j
        if(.not.allowedRot(i,j))then
          beta(i,j) = 0.d0
          cycle
        end if
        beta(i,j) = getBeta(i,j)
        if(abs(beta(i,j))>max_b)then
          max_b=abs(beta(i,j))
          mi=i
          mj=j
        end if!(abs(beta(i,j))>max_b)
      end do!j=ldeg,i-1
    end do!i=ldeg,udeg
    if(printlvl>2) print "(6X,A,E12.5,2(A,I1))","max(|beta|)=",max_b,&
                         " ,block:",mi,",",mj
    do while(iter<maxiter.and.max_b>toler)
      iter=iter+1
      t=beta(mi,mj)
      c = cos(t)
      s = sin(t)
! Gnew = J^T.G.J.   Gnew_ij=Jki*Gkl*Jlj
      newgrad = pt%grads(:pt%nvibs,:,:)
      newener = pt%energy
      do i=1,nstates
        newgrad(:,i,mi)=pt%grads(:pt%nvibs,i,mi)*c-pt%grads(:pt%nvibs,i,mj)*s
        newgrad(:,i,mj)=pt%grads(:pt%nvibs,i,mj)*c+pt%grads(:pt%nvibs,i,mi)*s
        newener(i,mi)=pt%energy(i,mi)*c-pt%energy(i,mj)*s
        newener(i,mj)=pt%energy(i,mj)*c+pt%energy(i,mi)*s
      end do
      pt%grads(:pt%nvibs,:,:) = newgrad
      pt%energy               = newener
      do i=1,nstates
        pt%grads(:pt%nvibs,mi,i)=newgrad(:,mi,i)*c-newgrad(:,mj,i)*s
        pt%grads(:pt%nvibs,mj,i)=newgrad(:,mj,i)*c+newgrad(:,mi,i)*s 
        pt%energy(mi,i)=newener(mi,i)*c-newener(mj,i)*s
        pt%energy(mj,i)=newener(mj,i)*c+newener(mi,i)*s 
      end do
      !update rotation angles
      do i=mj+1,mi-1
        if(allowedRot(i,mj)) beta(i,mj)=getBeta(i,mj)
      end do
      do i=mi+1,udeg
        if(allowedRot(i,mi)) beta(i,mi)=getBeta(i,mi)
        if(allowedRot(i,mj)) beta(i,mj)=getBeta(i,mj)
      end do
      do j=ldeg,mi-1
        if(allowedRot(mi,j)) beta(mi,j)=getBeta(mi,j)
      end do
      do j=ldeg,mj-1
        if(allowedRot(mj,j)) beta(mj,j)=getBeta(mj,j)
      end do
      max_b=-1.d0
      do i=ldeg,udeg
        do j=ldeg,i-1
          if(.not.allowedRot(i,j)) cycle
          if(abs(beta(i,j))>max_b) then
            max_b=abs(beta(i,j))
            mi=i
            mj=j
          end if!(abs(beta(i,j))>max_b)
        end do!j=ldeg,i-1
      end do!i=ldeg,udeg      
      if(printlvl>2)print *,"   iteration ",iter,", max(|beta|)=",max_b
    end do!while(iter<maxiter.and.max_b>toler)do
    if(printlvl>1.and.iter>0)then
        if(max_b<toler)then
          print 1001,iter
        else
          print 1000,iter,max_b
        end if!(max_b<toler)
        print *,"      Energies after transformation" 
        do i=1,nstates
          print "(20F12.7)",pt%energy(i,:)
        end do
        print *,"      Couplings and gradients after transformation" 
        do i=1,nstates
          do j=1,i
            write (*,1002) i,j
            if(hasGrad(i,j))then
                print 1003,pt%grads(:pt%nvibs,i,j)
            else
                print "(10X,A)","Data not available"
            end if
          end do
        end do!i=1,nstates
    end if!(printlvl>0)
  end do grploop!igrp=1,pt%ndeggrp
1000 format(7X,"no convergence after ",I4," iterations, max residue angle=",F8.2)
1001 format(7X,"convergence after ",I4," iterations")
1002 format(14X,'states(',I2,',',I2,'):')
1003 format(12X,3E16.7)
CONTAINS
  FUNCTION getBeta(i,j) RESULT(beta)
    USE progdata, only : abpoint,pi
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: i,j
    DOUBLE PRECISION :: beta
    double precision, dimension(pt%nvibs) :: g,h

    g=(pt%grads(:pt%nvibs,i,i)-pt%grads(:pt%nvibs,j,j))*0.5d0
    h=pt%grads(:pt%nvibs,i,j)
    if(printlvl>2) then
      print "(A,I2,A,I2)","Calculating rotation between states ",i," and ",j
      print *,"G vector:"
      print "(3F12.8)",g
      print *,"H vector:"
      print "(3F12.8)",h
    end if
    CALL orthgh(pt%nvibs,g,h,beta)
    if(printlvl>2) print *,"beta=",beta
  END FUNCTION getBeta
END SUBROUTINE OrthGH_Ab 
!------------------------------------------------------------------------------
! Rotate Hd generated eigenvectors so that the resulting g and h vectors between any
! pairs of two degenerate states are orthogonalized.  g and h vectors are identified 
! by their overlap with the ab initio g and h vectors
! pt     : ab initio data of the point to perform transformation, already orthongalized 
! dhmat  : gradient of diabatic Hamiltonian matrix
! ckl    : eigenvectors corresponding to rotation that orthogonalize g and h vectors 
! maxiter: maximum number of jacobi rotations
! toler  : convergence tolerance for rotation angle beta
![method]
! This subroutine tries orthogonalize all g-h vectors by an algorithm similar to 
! the Jacobi's method for eigenvalues.  The Jacobi transformations here extremize 
! the norm of the corresponding coupling block rather than eliminating them.   
SUBROUTINE OrthGH_Hd(pt,dhmat,ckl,maxiter,toler)
  USE hddata, only: nstates
  USE progdata, only: abpoint,printlvl
  IMPLICIT NONE
  TYPE(abpoint),INTENT(IN) :: pt
  INTEGER,INTENT(IN) :: maxiter
  DOUBLE PRECISION,INTENT(IN) :: toler
  DOUBLE PRECISION,dimension(pt%nvibs,nstates,nstates),INTENT(IN) :: dhmat
  DOUBLE PRECISION,dimension(nstates,nstates),INTENT(INOUT) :: ckl

  integer :: igrp,i,j,iter,ldeg,udeg
  integer :: mi,mj  !location of maximum rotation
  !desired rotation angles and gh overlaps
  double precision, dimension(nstates,nstates) :: beta, gh
  !Hd gradients after rotation
  double precision, dimension(pt%nvibs,nstates,nstates) :: gradnew
  double precision :: max_gh,t

  if(pt%ndeggrp .le. 0) return

  !initialize eigenvectors and rotated gradients
  do i=1,pt%nvibs
    gradnew(i,:,:)=matmul(matmul(transpose(ckl),dhmat(i,:,:)),ckl)
  end do
  beta=0.d0
  if(printlvl>2)print *,"     transforming Hd eigenvectors for point ",pt%id
  if(printlvl>2)then
    print *,"      gradients before transformation"
    do i=1,nstates
      do j=1,i
        print 1002,i,j
        print 1003,gradnew(:,i,j)
      end do
    end do!i=1,nstates
  end if!(printlvl>2)

! generate rotations
grploop:  do igrp=1,pt%ndeggrp  ! in case there are multiple groups of degeneracies, loops over all of them
    max_gh=-1.d0
    ldeg=pt%deg_groups(igrp,1)
    udeg=pt%deg_groups(igrp,2)+pt%deg_groups(igrp,1)-1
    if(printlvl>2) print "(3(A,I3))","  Degenerate group ",igrp,", states ",ldeg," to ",udeg
    iter=0

    !build normalized g and h vectors from Hd and compute rotation angles the first time
    do i=ldeg,udeg 
      do j=ldeg,i-1
        call getBetaHd(i,j,beta(i,j),gh(i,j))
        if(gh(i,j)>max_gh)then
          max_gh=gh(i,j)
          mi=i
          mj=j
        end if!(abs(beta(i,j))>max_gh)
      end do!j=ldeg,i-1
    end do!i=ldeg,udeg
    if(printlvl>2)print "(6x,A,E15.7,A,I3,A,I3,A)","max(g.h)=",max_gh,&
                        "(",mi,",",mj,")"
    do while(iter<maxiter.and.max_gh>toler)
      iter=iter+1
      t=beta(mi,mj)
      call GivensStep(pt%nvibs,nstates,gradnew,mi,mj,t,ckl)
      !update rotation angles
      do i=mj+1,mi-1 
        call getBetaHd(i,mj,beta(i,mj),gh(i,mj))
      end do
      do i=mi+1,udeg
        call getBetaHd(i,mj,beta(i,mj),gh(i,mj))
        call getBetaHd(i,mi,beta(i,mi),gh(i,mi))
      end do
      do j=ldeg,mi-1
        call getBetaHd(mi,j,beta(mi,j),gh(mi,j))
      end do
      do j=ldeg,mj-1
        call getBetaHd(mj,j,beta(mj,j),gh(mj,j))
      end do
      max_gh=-1.d0
      do i=ldeg,udeg
        do j=ldeg,i-1
          if(gh(i,j)>max_gh)then
            max_gh=gh(i,j)
            mi=i
            mj=j
          end if!(abs(beta(i,j))>max_gh)
        end do!j=ldeg,i-1
      end do!i=ldeg,udeg      
      if(printlvl>2)print "(4x,A,I4,A,E15.7,A,I3,A,I3,A)","iteration ",iter,", max(g.h)=",max_gh,&
                        "(",mi,",",mj,")"
    end do!while(iter<maxiter.and.max_gh>toler)do

    if(printlvl>2.and.iter>0) then
        if(max_gh<toler) then
          print 1001,iter
        else
          print 1000,iter,max_gh
        end if!(max_gh<toler)
        print *,"      gradients after transformation" 
        do i=1,nstates
          do j=1,i
            write (*,1002) i,j
            write (*,1003) gradnew(:pt%nvibs,i,j)
          end do
        end do!i=1,nstates
    end if!(printlvl>0)

  end do grploop!igrp=1,pt%ndeggrp

1000 format(7X,"no convergence after ",I4," iterations, max residue g.h=",F8.2)
1001 format(7X,"convergence after ",I4," iterations")
1002 format(14X,'states(',I2,',',I2,'):')
1003 format(12X,3E20.11)
CONTAINS
  SUBROUTINE getBetaHd(i,j,beta,gh) 
    USE progdata, only : abpoint,pi
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: i,j
    DOUBLE PRECISION,INTENT(OUT) :: beta,gh
    double precision, dimension(pt%nvibs) :: g,h

    g=(gradnew(:,i,i)-gradnew(:,j,j))*0.5d0
    h=gradnew(:,i,j)
    if(printlvl>2) then
      print "(A,I2,A,I2)","Calculating rotation between states ",i," and ",j
      print *,"G vector norm=",sqrt(dot_product(g,g))
      print "(3F12.8)",g
      print *,"H vector norm=",sqrt(dot_product(h,h))
      print "(3F12.8)",h
    end if
    CALL orthgh(pt%nvibs,g,h,beta)
    gh=abs(dot_product(g,h))
    if(printlvl>2) print *,"beta=",beta
  END SUBROUTINE getBetaHd
END SUBROUTINE OrthGH_Hd 
!------------------------------------------------------------------------------
SUBROUTINE GivensStep(nvibs,nst,grad,mi,mj,beta,ckl)
  IMPLICIT NONE
  integer,intent(IN) :: nvibs,nst,mi,mj
  double precision,intent(IN)  :: beta
  double precision, dimension(nvibs,nst,nst),intent(INOUT)::  grad
  double precision, dimension(nst,nst),intent(INOUT)::  ckl 
  integer :: i
  double precision :: jacobi(nst,nst)

  jacobi=dble(0)  !construct Givens rotation matrix
  do i=1,nst
     jacobi(i,i)=dble(1)
  end do
  jacobi(mi,mi)=cos(beta)
  jacobi(mj,mj)=jacobi(mi,mi)
  jacobi(mi,mj)=sin(beta)
  jacobi(mj,mi)=-jacobi(mi,mj)
  do i=1,nvibs
    grad(i,:,:)=matmul(matmul(transpose(jacobi),grad(i,:,:)),jacobi)
  end do
  ckl = matmul(ckl,jacobi) !rotate eigenvectors
  return
END SUBROUTINE GivensStep
!----------------------------------------------------------------------------------
!adjust signs to match vectors
SUBROUTINE ConformSgn(nvibs,nst,grd,ref,ckl)
  implicit none
  integer,intent(in) :: nvibs,nst
  double precision,intent(in)   :: ref(nvibs,nst,nst)
  double precision,intent(inout):: ckl(nst,nst)
  double precision,intent(inout):: grd(nvibs,nst,nst)

  integer :: i,j,k,signs(nst),signs_best(nst)
  double precision :: err,err_best,grdtest(nvibs,nst,nst)
 
  signs=1
  err_best=-1.d0
  do i=1,2**(nst-1)
    do j=1,nst
      do k=1,nst
        grdtest(:,j,k)=grd(:,j,k)*signs(j)*signs(k)
      end do
    end do
    err=sum((grdtest-ref)**2)
    if(err<err_best.or.err_best<0.d0)then
      err_best=err
      signs_best=signs
    end if
    j=1
    do while (signs(j)==-1)
      signs(j)=1
      j=j+1
    end do 
    signs(j)=-1
  end do

! use the best signs
  do j=1,nst
    do k=1,nst
      grdtest(:,j,k)=grd(:,j,k)*signs_best(j)*signs_best(k)
    end do
    ckl(:,j)=ckl(:,j)*signs_best(j)
  end do
  grd=grdtest
END SUBROUTINE ConformSgn
!--------------------------------------------------------------------------------
! Rotate the representation so that the gradients of Hd conforms with a
! reference set of values.
! st1 and st2 defines the degenerate range
SUBROUTINE ConformVec(nvibs,nst,st1,st2,grd,ref,ckl,allowedRot)
  IMPLICIT NONE
  integer,intent(in) :: nvibs,nst,st1,st2
  double precision,intent(in)   :: ref(nvibs,nst,nst)
  logical,intent(in)            :: allowedRot(nst,nst)
  double precision,intent(inout):: ckl(nst,nst)
  double precision,intent(inout):: grd(nvibs,nst,nst)
  
  integer :: i,j
  double precision :: err, err_new
  double precision,parameter :: terr=1.d-3
! Simply rotate along each direction until difference between reference and new
! gradients is at its minimum.
  call getGHErr(nvibs,nst,st1,st2,grd,ref,err_new)
  err=err_new*2.d0
  do while(err-err_new>terr)
    err=err_new
    do i=st1,st2-1
      do j=i+1,st2
        if(.not.allowedRot(i,j))cycle
        call findRotMin(nvibs,nst,st1,st2,grd,ref,ckl,i,j)
      end do
    end do
    call getGHErr(nvibs,nst,st1,st2,grd,ref,err_new)
  end do 
END SUBROUTINE ConformVec
!----------------------------------------------------------------------------------
! rotate states s1 and s2 until error is minimized
! st1 and st2 defines the range of degenerate states
SUBROUTINE findRotMin(nvibs,nst,st1,st2,grd,ref,ckl,s1,s2)
  USE progdata, only : pi,printlvl
  implicit none
  integer,intent(in) :: nvibs,nst,s1,s2,st1,st2
  double precision,intent(in)   :: ref(nvibs,nst,nst)
  double precision,intent(inout):: ckl(nst,nst)
  double precision,intent(inout):: grd(nvibs,nst,nst)
  
  double precision,dimension(nvibs,nst,nst) :: grdtest
  double precision,dimension(nst,nst) :: ckltest
  double precision :: a,b,x1,x2,f_x1,f_x2
  integer :: k
  !golden ratio
  real*8, parameter :: tau=(sqrt(5.d0)-1.d0)*0.5d0
  integer,parameter :: nsteps=18
  integer,parameter :: iter=50
  double precision,parameter :: eps=1.d-3
  double precision,parameter :: step=pi/dble(nsteps)

  do k=-nsteps,nsteps
    x1=dble(k)*step
    f_x1=f(x1)
    if(f_x1 .lt. f_x2 .or. k .eq. -nsteps) then
      f_x2=f_x1
      x2=x1
    end if
  end do

  a=x2+step
  b=x2-step

  x1=a+(1.d0-tau)*(b-a)
  x2=a+tau*(b-a)
  f_x1=f(x1)
  f_x2=f(x2)

  write(*,"(3f15.9)") a,b,b-a

  k=0
  do while(abs(b-a) .gt. eps .and. k .lt. iter)
    k=k+1
    if(f_x1 .lt. f_x2) then
      b=x2
      x2=x1
      f_x2=f_x1
      x1=a+(1.d0-tau)*(b-a)
      f_x1=f(x1)
      write(*,"(3f15.9)") a,b,b-a
    else
      a=x1
      x1=x2
      f_x1=f_x2
      x2=a+tau*(b-a)
      f_x2=f(x2)
      write(*,"(3f15.9)") a,b,b-a
    end if
  end do

  if(f_x2.lt.f_x1) then
    x1=x2
    f_x1=f_x2
  end if

  if(printlvl>3) print "(A,2I3,3(A,F12.7))","Rotating states ",s1,s2,&
                       ", error : ",0.d0," -> ",f_x1," , angle=",x1
  call GivensStep(nvibs,nst,grd,s1,s2,x1,ckl)

  return

  contains
  function f(x)
    implicit none
    real*8 :: f
    real*8, intent(in) :: x
    grdtest=grd
    ckltest=ckl
    call GivensStep(nvibs,nst,grdtest,s1,s2,x,ckltest)
    call getGHErr(nvibs,nst,st1,st2,grdtest,ref,f)
    return
  end 
END SUBROUTINE findRotMin
!--------------------------------------------------------------------------------
! Adjust g and h vectors in the group of states from st1 to st2
SUBROUTINE AdjustGH(nvibs,nst,st1,st2,grd,ref,ckl)
  USE progdata, only : pi
  implicit none
  integer,intent(in) :: nvibs,nst,st1,st2
  double precision,intent(in)   :: ref(nvibs,nst,nst)
  double precision,intent(inout):: ckl(nst,nst)
  double precision,intent(inout):: grd(nvibs,nst,nst)
  double precision,dimension(nvibs,nst,nst):: grdtest,grdbest
  double precision :: errbest,err,ckltest(nst,nst),cklbest(nst,nst),rand
  integer i,j,k,nrot,rotlist(nst*(nst-1)/2,2)
  integer rotind(nst*(nst-1)/2)
  double precision,parameter :: t=pi/4.d0

!call getGHErr(nvibs,nst,grd,ref,err)
!print "('Angle change:',F10.4,' , Error:',E14.7)",0d0,err
!DO I=1,63
!  call GivensStep(nvibs,nst,grd,1,2,1d-1,ckl)  
!  call getGHErr(nvibs,nst,grd,ref,err)
!  print "('Angle change:',F10.4,' , Error:',E14.7)",1d-1*i,err
!
!
!END DO
!
!STOP

  call getGHErr(nvibs,nst,st1,st2,grd,ref,err)
  if(err<1d-2)then
    print *,"Error very small.  Skipping adjustments for GH switchings."
    return
  end if
  print *,"Adjusting for GH vector swithcings. Initial error=",err
  nrot=(st2-st1+1)*(st2-st1)/2
  k=1
  do i=st1,st2-1
    do j=i+1,st2
      rotlist(k,1)=i
      rotlist(k,2)=j
      k=k+1
    end do
  end do
  rotind=0
  grdtest=grd
  ckltest=ckl
  errbest=-1d0
  do i=1,4**nrot
    call getGHErr(nvibs,nst,st1,st2,grdtest,ref,err)
    print "(A,I3,A,E14.7,A,5I3)","i=",i,", Error=",err,", Rotation string:",rotind
    if(errbest<0.or.err<errbest)then
      grdbest=grdtest
      errbest=err
      cklbest=ckltest
    end if
    call GivensStep(nvibs,nst,grdtest,rotlist(1,1),rotlist(1,2),t,ckltest)
    rotind(1)=rotind(1)+1
    do j=2,nrot
      if(rotind(j-1)==4)then
        rotind(j-1)=0
        rotind(j)=rotind(j)+1
        call GivensStep(nvibs,nst,grdtest,rotlist(j,1),rotlist(j,2),t,ckltest)
      end if
    end do
  end do
  k=0
  ! do random rotations if no good combination is found
!  if(errbest>1d-2)print *,"Error too large.  Performing random rotations."
!  do while(errbest>1d-2 .and. k<500)
!    k=k+1
!    grdtest=grdbest
!    ckltest=cklbest
!    do i=1,5
!      call random_number(rand)
!      j=int(rand*3)+1
!      call random_number(rand)
!      rand=rand*2*pi
!      call GivensStep(nvibs,nst,grdtest,rotlist(j,1),rotlist(j,2),rand,ckltest)
!      call getGHErr(nvibs,nst,grdtest,ref,err)
!      if(errbest<0.or.err<errbest)then
!        print "(2(A,I3),A,E14.7)","Found lower error. k=",k,", i=",i,", Error=",err
!        grdbest=grdtest
!        errbest=err
!        cklbest=ckltest
!      end if
!    end do
!  end do
  grd=grdbest
  ckl=cklbest
END SUBROUTINE AdjustGH
!-------------------------------------------------------------------------------
! calculate the g h vector error in the group of states from st1 to st2
SUBROUTINE getGHErr(nvibs,nst,st1,st2,grd,ref,err)
  implicit none
  integer,intent(in) :: nvibs,nst,st1,st2
  double precision,dimension(nvibs,nst,nst),intent(in) :: grd,ref
  double precision,intent(out)  :: err
  double precision :: vec(nvibs)
  integer :: i,j
  double precision :: e1,e2

  err = 0.d0
  do i=st1,st2
    vec=grd(:,i,i)-ref(:,i,i)
    err=err+dot_product(vec,vec)
    do j=i+1,st2
      vec=grd(:,i,j)-ref(:,i,j)
      e1=dot_product(vec,vec)
      vec=grd(:,i,j)+ref(:,i,j)
      e2=dot_product(vec,vec)
      err=err+min(e1,e2)
    end do
  end do
  err=sqrt(err)
  return
END SUBROUTINE getGHErr
!------------------------------------------------------------------------------
