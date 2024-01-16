!Common geometry transformation applied in molecule computation:
!    Uniquify geometry
!    Cartesian <-> internal coodinate
!    Normal mode and vibrational frequency
!
!Reference: E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular viobrations: the theory of infrared and Raman vibrational spectra* (Dover, 1980)
module GeomTrans
    implicit none

!Derived type
    !Example: type(IntCoordDef),allocatable,dimension(:)::IntCoordDef
    !         IntCoordDef(iIntCoord) stands for the iIntCoord-th internal coordinate
    !         IntCoordDef(iIntCoord)%NMotions is the number of motions involved in this internal coordinate
    !         IntCoordDef(iIntCoord)%motion(iMotion) stands for its iMotion-th motion
    !         IntCoordDef(iIntCoord)%motion(iMotion)%type is its type
    !         IntCoordDef(iIntCoord)%motion(iMotion)%coeff is its normalized linear combination coefficient
    !         IntCoordDef(iIntCoord)%motion(iMotion)%atom(i) is its i-th involved atom
    type InvolvedMotion
        !Currently only support stretching, bending, torsion, OutOfPlane
        !stretching: the motion coordinate is bond length atom1_atom2
        !bending   : the motion coordinate is bond angle atom1_atom2_atom3, range [0,pi]
        !            derivative encounters singularity at pi
        !torsion   : the motion coordinate is dihedral angle atom1_atom2_atom3_atom4, range (-pi,pi]
        !            specifically, angle between plane 123 and plane 234
        !            dihedral angle has same sign to n_123 x n_234 . r_23
        !            where n_abc (the normal vector of plane abc) is the unit vector along r_ab x r_bc
        !            dihedral angle value encounters discontinuity at pi
        !OutOfPlane: the motion coordinate is out of plane angle atom1_atom2_atom3_atom4, range [-pi/2,pi/2]
        !            specifically, bond 12 out of plane 234
        character*10::type
        real*8::coeff
        integer,allocatable,dimension(:)::atom
    end type InvolvedMotion
    type IntCoordDef!short for INTernal COORDinate DEFinition
        integer::NMotions
        type(InvolvedMotion),allocatable,dimension(:)::motion
    end type IntCoordDef

!GeometryTransformation module only variable
    type(IntCoordDef),allocatable,dimension(:)::GeometryTransformation_IntCoordDef!short for INTernal COORDinate DEFinition

    !target internal coordinates in int->cart transformation
    real*8, allocatable :: trgt(:)
    integer :: ncart, nintc 

contains
!---------------------------------------------------------------------------------------
!cross_product(a,b) = a x b
function cross_product(a, b)
    real*8,dimension(3),intent(in)::a,b
    real*8,dimension(3)::cross_product
    cross_product(1)=a(2)*b(3)-a(3)*b(2)
    cross_product(2)=a(3)*b(1)-a(1)*b(3)
    cross_product(3)=a(1)*b(2)-a(2)*b(1)
end function cross_product
!---------------------------------------------------------------------------------------
!triple_product(a,b,c) = ( a x b ) . c
    real*8 function triple_product(a, b, c)
        real*8,dimension(3),intent(in)::a,b,c
        triple_product=c(1)*(a(2)*b(3)-a(3)*b(2))+c(2)*(a(3)*b(1)-a(1)*b(3))+c(3)*(a(1)*b(2)-a(2)*b(1))
    end function triple_product
!---------------------------------------------------------------------------------------
subroutine angle(a,b,c,phi,dphi)
  !return angle a-b-c and its derivative
  implicit none
  real*8, intent(in) :: a(3), b(3), c(3)
  real*8, intent(out) :: phi, dphi(9)
  real*8 :: ab(3), cb(3), ab1, cb1, t, t1, fac
  real*8 :: dab(3,9), dcb(3,9), dab1(9), dcb1(9), dt(9), dt1(9)
  integer :: i

  ab=a-b
  cb=c-b

  ab1=norm2(ab)
  cb1=norm2(cb)

  t=dot_product(ab,cb)
  t1=t/(ab1*cb1)

  call bdchk(t1)
  phi=dacos(t1)

  dab=0.d0
  dcb=0.d0
  dab(1,1)=1.d0
  dab(2,2)=1.d0
  dab(3,3)=1.d0
  dab(1,4)=-1.d0
  dab(2,5)=-1.d0
  dab(3,6)=-1.d0
  dcb(1,7)=1.d0
  dcb(2,8)=1.d0
  dcb(3,9)=1.d0
  dcb(1,4)=-1.d0
  dcb(2,5)=-1.d0
  dcb(3,6)=-1.d0

  do i=1,9
    dt(i)=dot_product(dab(:,i),cb)+dot_product(ab,dcb(:,i))
  end do

  dab1(1:3)=ab/ab1
  dab1(4:6)=-ab/ab1
  dcb1(4:6)=-cb/cb1
  dcb1(7:9)=cb/cb1

  dt1=-(dab1*cb1+ab1*dcb1)*t/(ab1*cb1)**2+dt/(ab1*cb1)

  fac=-1.d0/sqrt(abs(1.d0-t1**2))
  dphi=fac*dt1

  return
end
!---------------------------------------------------------------------------------------
subroutine norm_cross_product(a,b,c,V,dV)
  implicit none
  real*8, intent(in) :: a(3), b(3), c(3)
  real*8, intent(out) :: V(3), dV(3,9)
  real*8 :: ab(3), cb(3), V1, dV1(9)
  integer :: i

  ab=a-b
  cb=c-b

  V(1)=ab(2)*cb(3)-ab(3)*cb(2)
  V(2)=ab(3)*cb(1)-ab(1)*cb(3)
  V(3)=ab(1)*cb(2)-ab(2)*cb(1)

  V1=norm2(V)

  dV=0.d0
  dV(1,2)=cb(3)
  dV(1,3)=-cb(2)
  dV(1,5)=-cb(3)+ab(3)
  dV(1,6)=-ab(2)+cb(2)
  dV(1,8)=-ab(3)
  dV(1,9)=ab(2)

  dV(2,1)=-cb(3)
  dV(2,3)=cb(1)
  dV(2,4)=-ab(3)+cb(3)
  dV(2,6)=-cb(1)+ab(1)
  dV(2,7)=ab(3)
  dV(2,9)=-ab(1)

  dV(3,1)=cb(2)
  dV(3,2)=-cb(1)
  dV(3,4)=-cb(2)+ab(2)
  dV(3,5)=-ab(1)+cb(1)
  dV(3,7)=-ab(2)
  dV(3,8)=ab(1)

  do i=1,9
    dV1(i)=dot_product(V,dV(:,i))
  end do
  dV1=dV1/V1

  do i=1,9
    dV(:,i)=dV(:,i)/V1-V*dV1(i)/V1**2
  end do
  V=V/V1

  return
end
!---------------------------------------------------------------------------------------
!---------- Cartesian <-> Internal ----------
    !An interal coordinate is the linear combination of several translationally and rotationally invariant displacements
    !    but only displacements under same unit can be combined, i.e. you must treat length and angle separately
    !    unless appropriate metric tensor is applied
    !It is OK to define more than 3NAtoms-6 (or 3NAtoms-5 for linear molecule) internal coordinates,
    !    but only 3NAtoms-6 (or 3NAtoms-5 for linear molecule) partial derivatives are independent
    !Although the transformation from Cartesian coordinate to internal coordinate is not necessarily linear
    !    for infinitesimal displacement it is linear, corresponding to a matrix form: dq = B . dr
    !    where dq is internal coordinate differentiation, dr is Cartesian coordinate differentiation
    !    B is Jacobian(q,r) (historically called Wilson B matrix)
    !r is a 3NAtoms order vector with r[3*i-2:3*i] corresponding to the coordinate of i-th atom
    !Nomenclature:
    !    cartdim & intdim: Cartesian & internal space dimensionality
    !    cartgrad & intgrad: Cartesian & internal coordinate gradient (cartdim & intdim x NStates x NStates 3rd-order tensor)

    !Define internal coordinate, return the internal space dimensionality
    !Input:  format: internal coordinate definition format (Available: Columbus7, default)
    !        (optional) file: (default = 'intcfl' for Columbus7, 'IntCoordDef' for default) internal coordinate definition file name
    !Output: the internal space dimensionality
    !        also set the module-wide variable GeometryTransformation_IntCoordDef
    !        which will be refered by all routines in this section
    !See InvolvedMotion in 'Derived type' section for available types and ordering of atoms
    integer function DefineInternalCoordinate()
        integer::intdim
        if(allocated(GeometryTransformation_IntCoordDef)) deallocate(GeometryTransformation_IntCoordDef)
        call Columbus7()
        DefineInternalCoordinate=intdim
        contains
        !First line is always 'TEXAS'
        !New internal coordinate line starts with 'K'
        subroutine Columbus7()
            integer::NDef
            integer,allocatable,dimension(:)::NewLine
            character*10,allocatable,dimension(:)::MotionType
            character*24::chartemp; integer::i,j,k; real*8::dbletemp

            open(unit=99,file='intcfl',status='old')
                !The number of motion definition lines & internal coordinates
                    NDef=0; intdim=0; read(99,*)
                    do
                        read(99,'(A24)',iostat=i)chartemp
                        if(i/=0&!End of file or no definition
                        .or.(index(chartemp,'STRE')==0.and.index(chartemp,'BEND')==0&
                        .and.index(chartemp,'TORS')==0.and.index(chartemp,'OUT' )==0&
                        .and.index(chartemp,'LIN1')==0.and.index(chartemp,'LIN2')==0)) exit
                        NDef=NDef+1
                        if(scan(chartemp,'K')==1) intdim=intdim+1
                    end do; rewind 99
                !New internal coordinate lines & motions of line
                    allocate(NewLine(intdim+1)); NewLine(intdim+1)=NDef+1
                    allocate(MotionType(NDef))
                    k=1; read(99,*)
                    do i=1,NDef
                        read(99,'(A24)')chartemp
                        if(scan(chartemp,'K')==1) then; NewLine(k)=i; k=k+1; end if
                        if(index(chartemp,'STRE')>0) then; MotionType(i)='stretching'
                          else if(index(chartemp,'BEND')>0) then; MotionType(i)='bending'
                          else if(index(chartemp,'TORS')>0) then; MotionType(i)='torsion'
                          else if(index(chartemp,'OUT')>0) then; MotionType(i)='OutOfPlane'
                          else if(index(chartemp,'LIN1')>0) then; MotionType(i)='lin1'
                          else if(index(chartemp,'LIN2')>0) then; MotionType(i)='lin2'
                        end if
                    end do; rewind 99
                !Finally read internal coordinate definition. Linear combinations are normalized
                    allocate(GeometryTransformation_IntCoordDef(intdim))
                    k=1; read(99,*)
                    do i=1,intdim
                        GeometryTransformation_IntCoordDef(i)%NMotions=NewLine(i+1)-NewLine(i)
                        allocate(GeometryTransformation_IntCoordDef(i)%motion(GeometryTransformation_IntCoordDef(i)%NMotions))
                        if(GeometryTransformation_IntCoordDef(i)%NMotions==1) then
                            GeometryTransformation_IntCoordDef(i)%motion(1)%type=MotionType(k)
                            GeometryTransformation_IntCoordDef(i)%motion(1)%coeff=1d0
                            select case(MotionType(k))
                            case('stretching')
                                allocate(GeometryTransformation_IntCoordDef(i)%motion(1)%atom(2))
                                read(99,'(A28,I5,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom
                            case('bending')
                                allocate(GeometryTransformation_IntCoordDef(i)%motion(1)%atom(3))
                                read(99,'(A28,I6,1x,I9,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(1),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(3),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(2)
                            case('torsion')
                                allocate(GeometryTransformation_IntCoordDef(i)%motion(1)%atom(4))
                                read(99,'(A28,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom
                            case('OutOfPlane')
                                allocate(GeometryTransformation_IntCoordDef(i)%motion(1)%atom(4))
                                read(99,'(A28,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(1),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(3),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(4),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(2)
                            case('lin1')
                                allocate(GeometryTransformation_IntCoordDef(i)%motion(1)%atom(4))
                                read(99,'(A28,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(1),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(2),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(3),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(4)
                            case('lin2')
                                allocate(GeometryTransformation_IntCoordDef(i)%motion(1)%atom(4))
                                read(99,'(A28,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(1),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(2),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(3),&
                                GeometryTransformation_IntCoordDef(i)%motion(1)%atom(4)
                            case default; write(*,*)'Program abort: unsupported internal coordinate type '//trim(adjustl(MotionType(k))); stop
                            end select
                            k=k+1
                        else
                            dbletemp=0d0
                            do j=1,GeometryTransformation_IntCoordDef(i)%NMotions
                                GeometryTransformation_IntCoordDef(i)%motion(j)%type=MotionType(k)
                                select case(MotionType(k))
                                case('stretching')
                                    allocate(GeometryTransformation_IntCoordDef(i)%motion(j)%atom(2))
                                    read(99,'(A10,F10.7,8x,I5,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%coeff,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom
                                case('bending')
                                    allocate(GeometryTransformation_IntCoordDef(i)%motion(j)%atom(3))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%coeff,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(1),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(3),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(2)
                                case('torsion')
                                    allocate(GeometryTransformation_IntCoordDef(i)%motion(j)%atom(4))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%coeff,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom
                                case('OutOfPlane')
                                    allocate(GeometryTransformation_IntCoordDef(i)%motion(j)%atom(4))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%coeff,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(1),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(3),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(4),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(2)
                                case('lin1')
                                    allocate(GeometryTransformation_IntCoordDef(i)%motion(j)%atom(4))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%coeff,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(1),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(2),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(3),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(4)
                                case('lin2')
                                    allocate(GeometryTransformation_IntCoordDef(i)%motion(j)%atom(4))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%coeff,&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(1),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(2),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(3),&
                                    GeometryTransformation_IntCoordDef(i)%motion(j)%atom(4)
                                case default; write(*,*)'Program abort: unsupported internal coordinate type '//trim(adjustl(MotionType(k))); stop
                                end select
                                k=k+1
                                dbletemp=dbletemp+&
                                   GeometryTransformation_IntCoordDef(i)%motion(j)%coeff**2
                            end do
                            dbletemp=Sqrt(abs(dbletemp))
                            forall(j=1:GeometryTransformation_IntCoordDef(i)%NMotions)
                                GeometryTransformation_IntCoordDef(i)%motion(j)%coeff=&
                                GeometryTransformation_IntCoordDef(i)%motion(j)%coeff/dbletemp
                            end forall
                        end if
                    end do
            close(99)
            deallocate(NewLine); deallocate(MotionType)!Clean up
        end subroutine Columbus7
    end function DefineInternalCoordinate
!---------------------------------------------------------------------------------------
    !========== Cartesian -> Internal ==========
        !Convert r to q
        subroutine InternalCoordinate(r, q, cartdim, intdim)
            integer,intent(in)::cartdim,intdim
            real*8,dimension(cartdim),intent(in)::r
            real*8,dimension(intdim),intent(out)::q
            integer::iIntC,iMotion
            q=0d0
            do iIntC=1,intdim
                do iMotion=1,GeometryTransformation_IntCoordDef(iIntC)%NMotions
                    select case(GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%type)
                    case('stretching')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%coeff&
                            *stretching(r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('bending')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%coeff&
                            *bending(r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('torsion')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%coeff&
                            *torsion(r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('OutOfPlane')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%coeff&
                            *OutOfPlane(r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('lin1')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%coeff&
                            *lin1(r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('lin2')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%coeff&
                            *lin2(r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    end select
                end do
            end do
            contains
            !Transform from Cartesian coordinate r to a certain motion coordinate q, atom defines which atoms are involved
            !For stretching, q = bond length
            real*8 function stretching(r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(2),intent(in)::atom
                real*8,dimension(3)::r12
                r12=r(3*atom(2)-2:3*atom(2))-r(3*atom(1)-2:3*atom(1))
                stretching=Norm2(r12)
            end function stretching
            !For bending, q = bond angle
            real*8 function bending(r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(3),intent(in)::atom
                real*8,dimension(3)::runit21,runit23
                real*8 :: tt
                runit21=r(3*atom(1)-2:3*atom(1))-r(3*atom(2)-2:3*atom(2))
                    runit21=runit21/Norm2(runit21)
                runit23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                    runit23=runit23/Norm2(runit23)
                tt=dot_product(runit21,runit23)
                call bdchk(tt)
                bending=acos(tt)
            end function bending
            !For torsion, q = dihedral angle
            real*8 function torsion(r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8,dimension(3)::r12,r23,r34,n123,n234
                real*8 :: tt
                r12=r(3*atom(2)-2:3*atom(2))-r(3*atom(1)-2:3*atom(1))
                r23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                r34=r(3*atom(4)-2:3*atom(4))-r(3*atom(3)-2:3*atom(3))
                n123=cross_product(r12,r23); n123=n123/Norm2(n123)
                n234=cross_product(r23,r34); n234=n234/Norm2(n234)
                tt=dot_product(n123,n234)
                call bdchk(tt)
                torsion=acos(tt)
                if(triple_product(n123,n234,r23)<0d0) torsion=-torsion
            end function torsion
            !For out of plane, q = out of plane angle
            real*8 function OutOfPlane(r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8,dimension(3)::r21,r23,r24
                real*8 :: tt
                r21=r(3*atom(1)-2:3*atom(1))-r(3*atom(2)-2:3*atom(2))
                r23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                r24=r(3*atom(4)-2:3*atom(4))-r(3*atom(2)-2:3*atom(2))
                r23=cross_product(r23,r24)
                tt=dot_product(r23/norm2(r23),r21/norm2(r21))
                call bdchk(tt)
                OutOfPlane=asin(tt)
            end function OutOfPlane
            !colinear bending
            real*8 function lin1(r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8,dimension(3):: U,V,X
                real*8 :: phi1,phi2
                real*8 :: tt

                U=r(3*atom(1)-2:3*atom(1))-r(3*atom(3)-2:3*atom(3))
                U=U/norm2(U)

                V=r(3*atom(4)-2:3*atom(4))-r(3*atom(3)-2:3*atom(3))
                V=V/norm2(V)

                X=r(3*atom(2)-2:3*atom(2))-r(3*atom(3)-2:3*atom(3))
                X=X/norm2(X)

                tt=dot_product(U,V)
                call bdchk(tt)
                phi1=dacos(tt)

                tt=dot_product(X,V)
                call bdchk(tt)
                phi2=dacos(tt)
                lin1=dacos(-1.d0)-phi1-phi2

            end function lin1
            !perpendicular linear bending
            real*8 function lin2(r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8,dimension(3):: U,V,Z,X,W
                real*8 :: phi1,phi2
                real*8 :: tt

                U=r(3*atom(1)-2:3*atom(1))-r(3*atom(3)-2:3*atom(3))
                U=U/norm2(U)

                V=r(3*atom(4)-2:3*atom(4))-r(3*atom(3)-2:3*atom(3))
                V=V/norm2(V)

                Z=r(3*atom(2)-2:3*atom(2))-r(3*atom(3)-2:3*atom(3))
                Z=Z/norm2(Z)

                W=cross_product(V,U)
                W=W/norm2(W)

                X=cross_product(Z,V)
                X=X/norm2(X)

                tt=dot_product(U,W)
                call bdchk(tt)
                phi1=dacos(tt)

                tt=dot_product(Z,W)
                call bdchk(tt)
                phi2=dacos(tt)
                lin2=dacos(-1.d0)-phi1-phi2

            end function lin2
        end subroutine InternalCoordinate
!---------------------------------------------------------------------------------------
        !From r, generate B & q
        subroutine WilsonBMatrixAndInternalCoordinate(r, B, q, cartdim, intdim)
            integer,intent(in)::cartdim,intdim
            real*8,dimension(cartdim),intent(in)::r
            real*8,dimension(intdim,cartdim),intent(out)::B
            real*8,dimension(intdim),intent(out)::q
            integer::iIntC,iMotion; real*8::qMotion; real*8,dimension(cartdim)::BRowVector
            B=0d0; q=0d0
            do iIntC=1,intdim
                do iMotion=1,GeometryTransformation_IntCoordDef(iIntC)%NMotions
                    select case(GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%type)
                    case('stretching'); call bAndStretching(BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('bending')   ; call bAndBending   (BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('torsion')   ; call bAndTorsion   (BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('OutOfPlane'); call bAndOutOfPlane(BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('lin1'); call bAndlin1(BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    case('lin2'); call bAndlin2(BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%atom,cartdim)
                    end select
                    B(iIntC,:)=B(iIntC,:)+GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%coeff*BRowVector
                    q(iIntC)  =q(iIntC)  +GeometryTransformation_IntCoordDef(iIntC)%motion(iMotion)%coeff*qMotion
                end do
            end do
            contains
            !Generate the transformation vector b from dr to dq: b . dr = dq
            !Transform from Cartesian coordinate r to a certain motion coordinate q
            !Internal coordinate is the linear combination of several motions,
            !so b contributes (but not necessarily equals) to one row of Wilson B matrix
            ! d( i-th internal coordinate ) = ( i-th row vector of B ) . dr
            !For stretching, q = bond length
            subroutine bAndStretching(b, q, r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(2),intent(in)::atom
                real*8,dimension(3)::runit12
                b=0d0!Initialize
                runit12=r(3*atom(2)-2:3*atom(2))-r(3*atom(1)-2:3*atom(1))
                q=Norm2(runit12)
                runit12=runit12/q
                b(3*atom(1)-2:3*atom(1))=-runit12
                b(3*atom(2)-2:3*atom(2))=runit12
            end subroutine bAndStretching
            !For bending, q = bond angle
            subroutine bAndBending(b, q, r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(3),intent(in)::atom
                real*8::r21,r23,costheta,sintheta
                real*8,dimension(3)::runit21,runit23
                real*8 :: tt
                b=0d0!Initialize
                !Prepare
                runit21=r(3*atom(1)-2:3*atom(1))-r(3*atom(2)-2:3*atom(2))
                    r21=Norm2(runit21); runit21=runit21/r21
                runit23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                    r23=Norm2(runit23); runit23=runit23/r23
                costheta=dot_product(runit21,runit23); sintheta=dSqrt(abs(1.d0-costheta**2))
                !Output
                b(3*atom(1)-2:3*atom(1))=(costheta*runit21-runit23)/(sintheta*r21)
                b(3*atom(3)-2:3*atom(3))=(costheta*runit23-runit21)/(sintheta*r23)
                b(3*atom(2)-2:3*atom(2))=-b(3*atom(1)-2:3*atom(1))-b(3*atom(3)-2:3*atom(3))
                tt=costheta
                call bdchk(tt)
                q=acos(tt)
            end subroutine bAndBending
            !For torsion, q = dihedral angle
            subroutine bAndTorsion(b, q, r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8::r12,r23,r34,sin123,cos123,sin234,cos234
                real*8,dimension(3)::runit12,runit23,runit34,n123,n234
                real*8 :: tt
                b=0d0!Initialize
                !Prepare
                runit12=r(3*atom(2)-2:3*atom(2))-r(3*atom(1)-2:3*atom(1))
                r12=Norm2(runit12); runit12=runit12/r12
                runit23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                r23=Norm2(runit23); runit23=runit23/r23
                runit34=r(3*atom(4)-2:3*atom(4))-r(3*atom(3)-2:3*atom(3))
                r34=Norm2(runit34); runit34=runit34/r34
                cos123=-dot_product(runit12,runit23); sin123=dSqrt(abs(1d0-cos123*cos123))
                n123=cross_product(runit12,runit23)/sin123
                cos234=-dot_product(runit23,runit34); sin234=dSqrt(abs(1d0-cos234*cos234))
                n234=cross_product(runit23,runit34)/sin234
                !Output
                b(3*atom(1)-2:3*atom(1))=-n123/(r12*sin123)
                b(3*atom(2)-2:3*atom(2))=(r23-r12*cos123)/(r12*r23*sin123)*n123-cos234/(r23*sin234)*n234
                b(3*atom(3)-2:3*atom(3))=(r34*cos234-r23)/(r23*r34*sin234)*n234+cos123/(r23*sin123)*n123
                b(3*atom(4)-2:3*atom(4))= n234/(r34*sin234)
                tt=dot_product(n123,n234)
                call bdchk(tt)
                q=acos(tt)
                if(triple_product(n123,n234,runit23)<0d0) q=-q
            end subroutine bAndTorsion
            !For out of plane, q = out of plane angle
            subroutine bAndOutOfPlane(b, q, r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8::r21,r23,r24,sin324,cos324,sin324sq,sintheta,costheta,tantheta
                real*8,dimension(3)::runit21,runit23,runit24
                real*8 :: tt
                b=0d0!Initialize
                !Prepare
                runit21=r(3*atom(1)-2:3*atom(1))-r(3*atom(2)-2:3*atom(2))
                r21=Norm2(runit21); runit21=runit21/r21
                runit23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                r23=Norm2(runit23); runit23=runit23/r23
                runit24=r(3*atom(4)-2:3*atom(4))-r(3*atom(2)-2:3*atom(2))
                r24=Norm2(runit24); runit24=runit24/r24
                cos324=dot_product(runit23,runit24)
                sin324=dSqrt(abs(1d0-cos324*cos324)); sin324sq=sin324*sin324
                sintheta=triple_product(runit23,runit24,runit21)/sin324
                costheta=dSqrt(abs(1d0-sintheta*sintheta)); tantheta=sintheta/costheta
                !Output
                b(3*atom(1)-2:3*atom(1))=(cross_product(runit23,runit24)/costheta/sin324-tantheta*runit21)/r21
                b(3*atom(3)-2:3*atom(3))=(cross_product(runit24,runit21)/costheta/sin324-tantheta/sin324sq*(runit23-cos324*runit24))/r23
                b(3*atom(4)-2:3*atom(4))=(cross_product(runit21,runit23)/costheta/sin324-tantheta/sin324sq*(runit24-cos324*runit23))/r24
                b(3*atom(2)-2:3*atom(2))=-b(3*atom(1)-2:3*atom(1))-b(3*atom(3)-2:3*atom(3))-b(3*atom(4)-2:3*atom(4))
                tt=sintheta
                call bdchk(tt)
                q=asin(tt)
            end subroutine bAndOutOfPlane
            !colinear bending
            subroutine bAndlin1(b, q, r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8 :: phi1,dphi1(9),phi2,dphi2(9)
                call angle(r(3*atom(1)-2:3*atom(1)),r(3*atom(3)-2:3*atom(3)),r(3*atom(4)-2:3*atom(4)),phi1,dphi1)
                call angle(r(3*atom(2)-2:3*atom(2)),r(3*atom(3)-2:3*atom(3)),r(3*atom(4)-2:3*atom(4)),phi2,dphi2)
                q=dacos(-1.d0)-phi1-phi2
                b=0d0
                b(3*atom(1)-2:3*atom(1))=-dphi1(1:3)
                b(3*atom(2)-2:3*atom(2))=-dphi2(1:3)
                b(3*atom(3)-2:3*atom(3))=-dphi1(4:6)-dphi2(4:6)
                b(3*atom(4)-2:3*atom(4))=-dphi1(7:9)-dphi2(7:9)
                return
            end subroutine bAndlin1
            !perpendicular linear bending
            subroutine bAndlin2(b, q, r, atom, cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8,dimension(3):: U,Z,W
                real*8 :: U1,Z1
                real*8, allocatable :: dU(:,:), dZ(:,:), dW(:,:), dW1(:,:)
                real*8, allocatable :: dU1(:), dZ1(:)
                real*8 :: fac1, fac2, t1, t2
                integer :: i

                allocate(dU(3,cartdim),dZ(3,cartdim),dW(3,cartdim),dW1(3,9))
                allocate(dU1(cartdim),dZ1(cartdim))

                call norm_cross_product(r(3*atom(4)-2:3*atom(4)),r(3*atom(3)-2:3*atom(3)),&
                                        r(3*atom(1)-2:3*atom(1)),W,dW1)
                dW=0.d0
                dW(:,3*atom(4)-2:3*atom(4))=dW1(:,1:3)
                dW(:,3*atom(3)-2:3*atom(3))=dW1(:,4:6)
                dW(:,3*atom(1)-2:3*atom(1))=dW1(:,7:9)

                U=r(3*atom(1)-2:3*atom(1))-r(3*atom(3)-2:3*atom(3))
                U1=norm2(U)
                dU=0.d0
                dU(1,3*atom(1)-2)=1.d0
                dU(2,3*atom(1)-1)=1.d0
                dU(3,3*atom(1))=1.d0
                dU(1,3*atom(3)-2)=-1.d0
                dU(2,3*atom(3)-1)=-1.d0
                dU(3,3*atom(3))=-1.d0
                do i=1,cartdim
                  dU1(i)=dot_product(U,dU(:,i))
                end do
                dU1=dU1/U1

                do i=1,cartdim
                  dU(:,i)=dU(:,i)/U1-U*dU1(i)/U1**2
                end do
                U=U/U1

                Z=r(3*atom(2)-2:3*atom(2))-r(3*atom(3)-2:3*atom(3))
                Z1=norm2(Z)
                dZ=0.d0
                dZ(1,3*atom(2)-2)=1.d0
                dZ(2,3*atom(2)-1)=1.d0
                dZ(3,3*atom(2))=1.d0
                dZ(1,3*atom(3)-2)=-1.d0
                dZ(2,3*atom(3)-1)=-1.d0
                dZ(3,3*atom(3))=-1.d0
                do i=1,cartdim
                  dZ1(i)=dot_product(Z,dZ(:,i))
                end do
                dZ1=dZ1/Z1

                do i=1,cartdim
                  dZ(:,i)=dZ(:,i)/Z1-Z*dZ1(i)/Z1**2
                end do
                Z=Z/Z1

                t1=dot_product(U,W)
                call bdchk(t1)
                t2=dot_product(Z,W)
                call bdchk(t2)
                q=dacos(-1.d0)-dacos(t1)-dacos(t2)

                fac1=1.d0/sqrt(abs(1.d0-t1**2))
                fac2=1.d0/sqrt(abs(1.d0-t2**2))
                do i=1,cartdim
                  b(i)=fac1*(dot_product(dU(:,i),W)+dot_product(U,dW(:,i)))+&
                       fac2*(dot_product(dZ(:,i),W)+dot_product(Z,dW(:,i)))
                end do

                return
            end subroutine bAndlin2
        end subroutine WilsonBMatrixAndInternalCoordinate
!---------------------------------------------------------------------------------
!======= Cartesian <- Internal ==========
!Convert q to r, geometry is placed at standard orientation
subroutine CartesianCoordinate(q, r, intdim, cartdim)
    !Required argument
    integer,intent(in) :: intdim,cartdim
    real*8, dimension(intdim), intent(in) :: q
    real*8, dimension(cartdim), intent(inout) :: r
    real*8 :: Rot(3,3)
    real*8, allocatable :: x(:), fvec(:), fjac(:,:)
    real*8 :: tol
    integer :: info

    allocate(x(intdim),fvec(intdim),fjac(intdim,intdim))
    trgt(1:intdim)=q(1:intdim)
    call orientation(r,Rot,info)
    call xinout(cartdim,r,intdim,x,1)

    tol=1.d-16
    call hybrj1(hybrj1_f, intdim, x, fvec, fjac, intdim, tol, info)
    !if(info.ne.1) then
    !  print*,'Output INFO in CartesianCoordinate: ',info
    !  stop
    !end if

    call xinout(cartdim,r,intdim,x,0)

    return
end subroutine CartesianCoordinate
!-----------------------------------------------------------------------------------
subroutine xinout(cartdim,cart,intdim,xt,id)
  implicit none
  integer, intent(in) :: cartdim, intdim, id
  real*8, intent(inout) :: cart(cartdim), xt(intdim)
  integer :: i,k

  !id=0 xt->cart
  if(id.eq.0) then
    k=0
    do i=1,cartdim
      if(i.eq.9 .or. i.eq.13 .or. i.eq.14 .or. i.eq.15 .or. i.eq.17 .or. i.eq.18) then
        cart(i)=0.d0
      else
        k=k+1
        cart(i)=xt(k)
      end if
    end do
  end if

  !id=1 cart->xt
  if(id.eq.1) then
    k=0
    do i=1,cartdim
      if(i.eq.9 .or. i.eq.13 .or. i.eq.14 .or. i.eq.15 .or. i.eq.17 .or. i.eq.18) then
        cycle
      else
        k=k+1
        xt(k)=cart(i)
      end if
    end do
  end if

  return
end subroutine xinout
!---------------------------------------------------------------------------------
subroutine hybrj1_f(intdim, x, fvec, fjac, ldfjac, iflag)
  implicit none
  !function/jacobian subroutine for internal-to-cartesian transformation
  integer, intent(in) :: intdim, ldfjac, iflag
  real*8, intent(in) :: x(intdim)
  real*8, intent(out) :: fvec(intdim), fjac(ldfjac,intdim)
  real*8, allocatable :: cart(:),Bmat(:,:),xt(:),intc(:)
  integer :: i

  allocate(cart(ncart),Bmat(intdim,ncart),xt(intdim),intc(intdim))

  if(iflag.eq.0) then
    !do nothing
  else if(iflag.eq.1) then
    !fvec
    xt=x
    call xinout(ncart,cart,intdim,xt,0)
    call WilsonBMatrixAndInternalCoordinate(cart, Bmat, fvec, ncart, intdim)
    fvec(1:intdim)=fvec(1:intdim)-trgt(1:intdim)
  else if(iflag.eq.2) then
    !fjac
    xt=x
    call xinout(ncart,cart,intdim,xt,0)
    call WilsonBMatrixAndInternalCoordinate(cart, Bmat, intc, ncart, intdim)
    fjac=0.d0
    do i=1,intdim
      call xinout(ncart,Bmat(i,:),intdim,fjac(i,:),1)
    end do
  end if

  return
end
!---------------------------------------------------------------------------------
!Use Wilson GF method to obtain normal mode and vibrational frequency from
!Hessian in internal coordinate
!Input:     H    : internal coordinate Hessian
!           B    : Wilson B matrix
!          mass  : mass of each atom
!Output: freq  : vibrational angular frequencies (negative if imaginary)
!      intmode : internal coordinate normal modes contained in each column (Wilson L matrix)
!        Linv  : Wilson L^-1 matrix
!      cartmode: Cartesian coordinate normal modes contained in each column
!use Miyazawa method, see S. Califano Vibrational States and JCP 29 246(1958)

!---------------------------------------------------------------------------------
subroutine WilsonGFMethod(H, B, mass, freq, intmode, Linv, cartmode, intdim, NAtoms)
  integer,intent(in)::intdim,NAtoms
  real*8,dimension(intdim,intdim),intent(in)::H
  real*8,dimension(intdim,3*NAtoms),intent(in)::B
  real*8,dimension(NAtoms),intent(in)::mass
  real*8,dimension(intdim),intent(out)::freq
  real*8,dimension(intdim,intdim),intent(out)::intmode,Linv
  real*8,dimension(3*NAtoms,intdim),intent(out)::cartmode

  real*8,allocatable :: Btemp(:,:), T(:), Tinv(:,:)
  real*8,allocatable :: A(:,:), W(:,:), C(:,:)

  real*8, allocatable :: work(:)
  integer :: i, lwork, info

  allocate(Btemp(intdim,3*NAtoms),T(intdim),Tinv(intdim,intdim))
  allocate(A(intdim,intdim),W(intdim,intdim),C(intdim,intdim))
  lwork=5*intdim*intdim
  allocate(work(lwork))

  do i=1,NAtoms
    Btemp(:,3*i-2:3*i)=B(:,3*i-2:3*i)/mass(i)
  end do
  A=matmul(Btemp,transpose(B))
  call dsyev('V','U',intdim,A,intdim,T,work,lwork,info)
  if(minval(T).lt.0.d0) stop 'G matrix not positive definite in WilsonGFMethod!'

  T=abs(T)
  do i=1,intdim
    W(:,i)=A(:,i)*sqrt(T(i))
  end do

  C=matmul(transpose(W),matmul(H,W))
  call dsyev('V','U',intdim,C,intdim,freq,work,lwork,info)

  do i=1,intdim
    if(freq(i).lt.0.d0) then
      freq(i)=-sqrt(abs(freq(i)))
    else
      freq(i)=sqrt(abs(freq(i)))
    end if
  end do

  intmode=matmul(W,C)

  Tinv=0.d0
  do i=1,intdim
    Tinv(i,i)=1.d0/sqrt(T(i))
  end do
  Linv=matmul(transpose(C),matmul(Tinv,transpose(A)))

  !Convert internal coordinate normal mode to Cartesian coordinate normal mode
  do i=1,NAtoms
    Btemp(:,3*i-2:3*i)=B(:,3*i-2:3*i)/sqrt(abs(mass(i)))
  end do
  call dGeneralizedInverseTranspose(Btemp,intdim,3*NAtoms)
  cartmode=matmul(transpose(Btemp),intmode)

  return
end subroutine WilsonGFMethod
!--------------------------------------------------------------------------------------
end module GeomTrans
!=================================================================================
