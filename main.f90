!==============================================================================
program surfgen
  use progdata
  use makesurfdata
  implicit none

  call random_seed()
  call readinput
  call readdisps

  call makesurf

  stop
end program surfgen
!==============================================================================
SUBROUTINE readinput
  use progdata
  use hddata
  use makesurfdata
  use GeomTrans
  implicit none
  integer :: fid,i,j
  namelist /fitting/ npoints, enfDiab, epmax, w_energy, w_grad, w_fij, &
  gradcutoff, cpcutoff, deggrdbinding, deg_cap, lambda, eshift, energyT, &
  highEScale, nrmediff, ediffcutoff,fixref

  fixref=.true.
  call FLUnit(fid)
  open(fid,file='fit.in',delim='APOSTROPHE')
  read(unit=fid,nml=fitting)
  close(fid)

  deg_cap=deg_cap/au2cm
  gorder=deg_cap

  energyT=energyT/au2cm
  energyT_en=energyT
  highEScale_en=highEScale

  gradcutoff=gradcutoff/au2cm
  cpcutoff=cpcutoff/au2cm

  !progdata
  natoms=11
  printlvl=5
  nvibs=3*natoms-6

  !hddata
  nstates=2
  ncoord=DefineInternalCoordinate()
  print*,'System has ',ncoord, 'internal degrees of freedom.'
  if(ncoord.ne.nvibs) stop 'ncoord.ne.nvibs!'

  !print coordinate definition
  do i=1,ncoord
    do j=1,GeometryTransformation_IntCoordDef(i)%NMotions
      write(*,"(2i5,2x,A10,2x,e12.5,4i5)") i,j,&
                GeometryTransformation_IntCoordDef(i)%motion(j)%type,&
                GeometryTransformation_IntCoordDef(i)%motion(j)%coeff,&
                GeometryTransformation_IntCoordDef(i)%motion(j)%atom
    end do
  end do

  !readdisps
  enfptn   =  'energy.all'
  gmfptn   =  'geom.all'
  grdfptn  =  'cartgrd.drt1.state$.all'
  cpfptn   =  'cartgrd.nad.drt1.state$.drt1.state$.all'

  !makeLocalIntCoord
  intGradT = 1.d-8
  intGradS = intGradT

  !OrthGH_ab
  gcutoff = 1.d-8

  return
end SUBROUTINE readinput
!================================================================================
