!===============================================================================
subroutine calc_perf0(P)
  use hddata
  use makesurfdata
  implicit none
  real*8, intent(out) :: P
  real*8, allocatable :: hmat(:,:), dhmat(:,:,:)
  real*8, allocatable :: work(:),vtmp(:)
  !variables that specify equation
  integer :: pt, s1, s2, g, ieq
  integer :: lwork,info
  integer :: i,j,l

  lwork=999
  allocate(hmat(nstates,nstates),dhmat(nvibs,nstates,nstates),vtmp(nstates))
  allocate(work(lwork))

  !coef to Hij
  call Hij2coef(1)

  do i=1,npoints
    call EvaluateHd(dispgeoms(i),nvibs,hmat,dhmat)
    !eigenvector
    ckl(i,:,:)=hmat
    call dsyev('V','U',nstates,ckl(i,:,:),nstates,vtmp,work,lwork,info)
    if(info .ne. 0) stop 'dsyev failed in calc_perf!'
    !fitE
    fitE(i,:,:)=matmul(transpose(ckl(i,:,:)),matmul(hmat,ckl(i,:,:)))
    !fitG
    do j=1,nvibs
      fitG(i,j,:,:)=matmul(transpose(ckl(i,:,:)),matmul(dhmat(j,:,:),ckl(i,:,:)))
    end do
    !fixphase
    l=dispgeoms(i)%nvibs
    call fixphase(l,fitG(i,1:l,:,:),dispgeoms(i)%grads(1:l,:,:),ckl(i,:,:),&
                  phaseList,incGrad(i,:,:),gWeight(i,:,:))
  end do !i

  !calculate difference vectors
  do ieq = 1, neqs
    pt=eqmap(ieq,1)
    s1=eqmap(ieq,2)
    s2=eqmap(ieq,3)
    g =eqmap(ieq,4)
    if(g==0) then !is an energy
      if(s2>0) then !its adiabatic energy or adiabatic off diagonal terms(0)
        err(ieq) = fitE(pt,s1,s2) - dispgeoms(pt)%energy(s1,s2)
      else 
        stop 's2.le.0 in clac_perf0!'
      end if !(s2>0)
    else !(g==0)   is gradient
      err(ieq) = fitG(pt,g,s1,s2) - dispgeoms(pt)%grads(g,s1,s2)
    end if!(g==0)
  end do !ieq
  err=err*weight
  P = 0.5d0*(dot_product(err,err)+lambda*dot_product(coef,coef)) !regularization

  !calculate RMS Error
  call getError

  return
end
!===============================================================================
subroutine calc_jac0
  use progdata, only:natoms
  use hddata
  use makesurfdata
  use DiabaticHamiltonian
  implicit none
  real*8, allocatable :: dpdh(:,:,:)
  real*8, allocatable :: dpdhdx(:,:,:,:)   !w.r.t. local coordinates
  real*8, allocatable :: dpdhdr(:,:,:,:)   !w.r.t. internal coordinates
  real*8, allocatable :: dtmp(:)
  !W^IJ_i =<f^I|W_i|f^J>
  !D^IJ_i = W^IJ_i / (E^J-E^I)
  !V^IJ_i =<f^I|del W_i|f^J>
  real*8, allocatable :: WIJ(:,:,:), DIJ(:,:,:)
  real*8, allocatable :: VIJ(:,:,:,:)
  real*8, allocatable :: cc(:)
  real*8, allocatable :: q(:),t(:),dt(:,:)

  integer :: i,j,k,l,m,n,R,ico,igd,ipt,nt,id
  !variables that specify equation
  integer :: ieq,pt,s1,s2,g

  nt=nstates*(nstates+1)/2
  allocate(dpdh(npoints,ncons,nt))
  allocate(dpdhdx(npoints,ncons,nvibs,nt))
  allocate(dpdhdr(ncons,ncoord,nstates,nstates))
  allocate(dtmp(3*natoms))
  allocate(WIJ(npoints,nstates,nstates))
  allocate(VIJ(npoints,nvibs,nstates,nstates))
  allocate(DIJ(npoints,nstates,nstates))
  allocate(cc(npoints))

  allocate(q(Hd_intdim),t(NHdExpansionBasis))
  allocate(dt(Hd_intdim,NHdExpansionBasis))

  !coef to Hij
  call Hij2coef(1)

  dpdh=0.d0
  do i=1,npoints
    q=0.d0
    q(1:ncoord)=dispgeoms(i)%igeom(1:ncoord)
    do j=1,NHdExpansionBasis
      t(j)=ExpansionBasis(q,j)
      dt(:,j)=ExpansionBasisGradient(q,j)
    end do

    dpdhdr=0.d0
    if(fixref) then
      g=NHdExpansionBasis-nzf
    else
      g=NHdExpansionBasis
    end if
    k=0
    do j=1,NHdExpansionBasis
      !fix the coeffs of zeroth and first order if fixref
      if((Hd_EBNR(j)%order.eq.0 .or. Hd_EBNR(j)%order.eq.1) .and. fixref) cycle
      k=k+1
      l=0
      do m=1,nstates
        do n=m,nstates
          l=l+1
          dpdh(i,k+g*(l-1),l)=t(j)
          dpdhdr(k+g*(l-1),1:ncoord,m,n)=dt(1:ncoord,j)
        end do
      end do
      if(l.ne.nt) stop 'l.ne.nt in calc_jac0!'
    end do
    if(k.ne.g) stop 'k.ne.g in calc_jac0!'
    !fill up
    do m=1,nstates-1
      do n=m+1,nstates
        dpdhdr(:,:,n,m)=dpdhdr(:,:,m,n)
      end do
    end do

    do l=1,ncons
      id=0
      do j=1,nstates
        do k=j,nstates
          id=id+1
          call dgemv('T',ncoord,3*natoms,1.d0,dispgeoms(i)%bmat,&
                         ncoord,dpdhdr(l,:,j,k),1,0.d0,dtmp,1)
          dpdhdx(i,l,1:nvibs,id)=dtmp(1:nvibs)
        end do
      end do
      if(id.ne.nt) stop 'id.ne.nt in calc_perf0!'
    end do

  end do !i

  !WIJ,VIJ,DIJ
  do ico=1,ncons
    WIJ=0.d0
    VIJ=0.d0
    DIJ=0.d0
    !construct vectors W,V and D for the current coefficient
    do i=1,nstates
      do j=i,nstates
        do l=1,nstates
          do k=1,nstates
            call idx(k,l,nstates,id)
            !calculate the multiplication of ckl
            cc=ckl(:,k,i)*ckl(:,l,j)
            WIJ(:,i,j)=WIJ(:,i,j)+cc*dpdh(:,ico,id)
            !construct the Hamiltonian contribution of VIJ 
            do igd=1,nvibs
              VIJ(:,igd,i,j)=VIJ(:,igd,i,j)+cc*dpdhdx(:,ico,igd,id)
            end do!igd
          end do!K
        end do!L

        !construct DIJ from WIJ. non zero only on the off diagonal
        if(j.ne.i) then
          DIJ(:,i,j)=WIJ(:,i,j)/(fitE(:,j,j)-fitE(:,i,i))
          !fill the lower triangle of D and W
          WIJ(:,j,i)   = WIJ(:,i,j)
          DIJ(:,j,i)   = -DIJ(:,i,j)
          VIJ(:,:,j,i) = VIJ(:,:,i,j)
        end if !(J.ne.I)
      end do!J
    end do!I

    !update the wave function dependent part of Wij and Vij
    do i=1,nstates
      do j=i,nstates
        do k=1,nstates
          if(k.ne.i) WIJ(:,i,j)=WIJ(:,i,j)+DIJ(:,k,i)*fitE(:,k,j)
          if(k.ne.j) WIJ(:,i,j)=WIJ(:,i,j)+DIJ(:,k,j)*fitE(:,k,i)
        end do
        if(i.ne.j) WIJ(:,j,i)=WIJ(:,i,j)
        !calculate the second piece of VIJ with DIJ and hIJ(=fitG)
        do igd = 1,nvibs
          do k = 1,nstates
            if(k.ne.i) VIJ(:,igd,i,j) = VIJ(:,igd,i,j) + DIJ(:,k,i)*fitG(:,igd,k,j)
            if(k.ne.j) VIJ(:,igd,i,j) = VIJ(:,igd,i,j) + DIJ(:,k,j)*fitG(:,igd,k,i)
          end do ! K = 1,nstates
          !complete lower triangle
          if(i.ne.j)  VIJ(:,igd,j,i) = VIJ(:,igd,i,j)
        end do!igd=1,nvibs 
      end do !J=I,nstates
    end do !I=1,nstates

    !calculate derivative vector with respect to the current variable
    do ieq = 1, neqs
      pt=eqmap(ieq,1)
      s1=eqmap(ieq,2)
      s2=eqmap(ieq,3)
      g =eqmap(ieq,4)
      if(g==0) then !is an energy
        if(s2>0) then !its adiabatic energy or adiabatic off diagonal terms(0)
          Jac(ico,ieq) = WIJ(pt,s1,s2)
        else !(s2>0)
          stop 's2.le.0 in clac_jac0!'
        end if !(s2>0)
      else !(g==0) is gradient
        Jac(ico,ieq) = VIJ(pt,g,s1,s2)
      end if!(g==0)
      Jac(ico,ieq)=Jac(ico,ieq)*weight(ieq)
    end do !ieq

  end do!ico

  return
end
!===============================================================================
