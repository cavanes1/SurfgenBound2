!==============================================================================
subroutine makesurf
  use progdata, only:printlvl
  use hddata
  use makesurfdata
  use DiabaticHamiltonian
  implicit none
  character*32 :: FileName
  real*8 :: P,P0
  integer :: ep,i,j,k,info
  logical :: here

  call initmakesurf

  inquire(file='Hd.CheckPoint.old',exist=here)
  if(here) then
    FileName='Hd.CheckPoint.old'
    call ReadHdExpansionCoefficients(Hd_HdEC,FileName)
    call Hij2coef(0)
  else
    print*,'The reference point is ', enfDiab
    call setrefpt
    call Hij2coef(0)
    call WriteHdExpansionCoefficients(Hd_HdEC)
  end if

  mu=1.d0
  call calc_perf0(P0)
  do ep=1,epmax
    call copy_coef
    call calc_jac0
    call calc_jtj
    call calc_jte
100 call calc_jtjmui
    call rec_coef
    call update(info)
    if(info .ne. 0) then
      call rec_coef
      print*, 'Cholesky Decomposition Failed!'
      exit
    end if
    call calc_perf0(P)

    if(P .lt. P0) then
      P0=P
      mu=mu*0.5d0
    else
      mu=mu*2.d0
      if(mu .gt. 1.d10) then
        write(*,"('Large mu=',1x,e10.2,1x,'exit loop.' )") mu
        exit
      end if
      goto 100
    end if

    write(*,"('Epoch: ',i4,'/',i4,2x,'P=',es9.2,2x,'d[E]=',es9.2,2x,'<d[E]>=',&
               es9.2,2x,'d[g]=',es9.2,2x,'<d[g]>=',es9.2,2x,'d[hij]=',es9.2,2x,&
              'mu=',es7.0)") ep,epmax,P,rmsee,mue,rmseg,mueg,rmsec,mu
    call Hij2coef(1)
    call WriteHdExpansionCoefficients(Hd_HdEC)
  end do


  return
end
!==============================================================================
