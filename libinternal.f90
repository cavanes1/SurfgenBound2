!=============================================================================
subroutine buildWBmat(natom,ncoords,cgeom,igeom,Bmat)
  use GeomTrans
  implicit none
  integer, intent(in) :: natom, ncoords
  real*8, intent(in) :: cgeom(3*natom)
  real*8, intent(out) :: igeom(ncoords)
  real*8, intent(out) :: Bmat(ncoords,3*natom)

  call WilsonBMatrixAndInternalCoordinate(cgeom, Bmat, igeom, 3*natom, ncoords)

  return
end
!=============================================================================
