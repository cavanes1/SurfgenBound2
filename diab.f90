!=================================================================================
module DiabaticHamiltonian
!Derived type
    !This defines the mapping rule between an actual basis and its serial number
    !Example: type(ExpansionBasisDefinition),allocatable,dimension(:)::EBNR
    !         EBNR(i).order is the polynomial order for i-th basis function
    !         This function is the product of q(EBNR(i).indice(1))*...*q(EBNR(i).indice(EBNR(i).order))
    !         Where q is the internal coordinate (for this program it's internal coordinate difference)
    type ExpansionBasisNumberingRule
        !In this program, expansion basis function is polynomial of internal coordinate
        integer::order!The order of this polynomial
        integer,allocatable,dimension(:)::indice!The product of which internal coordinates forms this polynomial
    end type ExpansionBasisNumberingRule

    type d2PArray
        real*8,allocatable,dimension(:) :: Array
    end type d2PArray

    type i2PArray
        integer,allocatable,dimension(:)::Array
    end type i2PArray

    !Global variable
    !Number of Hd expansion basis functions, number of expansion coefficients
    integer::NHdExpansionBasis,NHdExpansionCoefficients

    !DiabaticHamiltonian module only variable
    !Basic information of Hd:
    !    NState: number of states Hd describes
    !    intdim: dimension of internal space Hd takes into account
    !    EBNR: see derived type section above
    integer::Hd_NState,Hd_intdim
    type(ExpansionBasisNumberingRule),allocatable,dimension(:)::Hd_EBNR!short for Expansion Basis Numbering Rule
    type(d2PArray),allocatable,dimension(:,:) :: Hd_HdEC!Short for Hd Expansion Coefficient, use only lower triangle
contains

!-------------- Hd definition ---------------
!This version generates Hd in NadVibS format: cast the elements of Hd into
!simple polynomials of the internal coordinate difference from reference geometry
!This treatment is good only for bounded system, where the molecule is semi-rigid

subroutine InitializeExpansionBasisNumberingRule()!Hd_EBNR
    integer::i,j
    NHdExpansionBasis=0
    open(unit=99,file='basis.in',status='old')
        do!Count how many expansion basis functions there are
            read(99,*,iostat=i); if(i/=0) exit
            NHdExpansionBasis=NHdExpansionBasis+1
        end do
        rewind 99
        allocate(Hd_EBNR(NHdExpansionBasis))

        !Read the definition of expansion basis functions
        do i=1,NHdExpansionBasis
            read(99,'(I5)',advance='no') Hd_EBNR(i)%order
            allocate(Hd_EBNR(i)%indice(Hd_EBNR(i)%order))
            if(Hd_EBNR(i)%order>0) then
                do j=1,Hd_EBNR(i)%order-1
                    read(99,'(I5)',advance='no')Hd_EBNR(i)%indice(j)
                end do
                read(99,'(I5)')Hd_EBNR(i)%indice(Hd_EBNR(i)%order)
            else
                read(99,*)
            end if
        end do
    close(99)
    return
end subroutine InitializeExpansionBasisNumberingRule
!----------------------------------------------------------------------------------------
    !Load Hd expansion coefficient from file to HdEC
    !Optional: FileName: (default = 'Hd.CheckPoint') name of the input file
    subroutine ReadHdExpansionCoefficients(HdEC,FileName)
        type(d2PArray),dimension(Hd_NState,Hd_NState),intent(inout)::HdEC
        character*32,optional,intent(in)::FileName
        character*2::char2temp; character*28::char28temp
        integer::NState,NBasis,NOrder!The old Hd is not necessarily fitted under same condition
        integer::istate,jstate,i,j,order,location
        integer,allocatable,dimension(:)::indice
        real*8::dbletemp
        if(present(FileName)) then; open(unit=99,file=FileName,status='old')
        else; open(unit=99,file='Hd.CheckPoint',status='old'); end if
            read(99,'(A28,I2)')char28temp,NState!Get old Hd fitting condition
            read(99,*); read(99,*)dbletemp
            read(99,'(I5)')NOrder; allocate(indice(NOrder))
            NBasis=2!Get number of basis functions
            do
                read(99,'(A2)')char2temp; if(char2temp=='Hd') exit
                NBasis=NBasis+1
            end do
            NBasis=NBasis/2
            rewind 99!Number of basis functions has been gotten
            read(99,*)!Read old Hd
            do istate=1,NState
                do jstate=istate,NState
                    read(99,*)
                    do i=1,NBasis
                        read(99,*)dbletemp
                        read(99,'(I5)',advance='no')order
                        if(order>0) then
                            do j=1,order-1; read(99,'(I5)',advance='no')indice(j); end do; read(99,'(I5)')indice(order)
                        else
                            read(99,*)
                        end if
                        location=WhichExpansionBasis(order,indice(1:order))
                        if(location>0) then
                            HdEC(jstate,istate)%Array(location)=dbletemp
                        else
                            write(*,*)'Warning: an old basis function does not belong to current basis space:'
                            do j=1,order-1; write(*,'(I5)',advance='no')indice(j); end do; write(*,'(I5)')indice(order)
                            write(*,*)'Its coefficient = ',dbletemp
                        end if
                    end do
                end do
            end do

        !fill up
        do istate=1,NState-1
          do jstate=istate+1,NState
            HdEC(istate,jstate)%Array=HdEC(jstate,istate)%Array
          end do
        end do
        close(99)
        return
    end subroutine ReadHdExpansionCoefficients
!----------------------------------------------------------------------------------------
    integer function WhichExpansionBasis(order,indice)!Return the location of the specified basis in Hd_EBNR, 0 if not found
        integer,intent(in)::order
        integer,dimension(order),intent(in)::indice
        call bisect(1,NHdExpansionBasis)
        contains
        recursive subroutine bisect(low,up)
            integer,intent(in)::low,up
            integer::bisection,i
            if(up-low==1) then
                if(order==Hd_EBNR(low).order) then
                    do i=order,1,-1
                        if(indice(i)/=Hd_EBNR(low).indice(i)) exit
                    end do
                    if(i<1) then
                        WhichExpansionBasis=low
                        return
                    end if
                end if
                if(order==Hd_EBNR(up).order) then
                    do i=order,1,-1
                        if(indice(i)/=Hd_EBNR(up).indice(i)) exit
                    end do
                    if(i<1) then
                        WhichExpansionBasis=up
                        return
                    end if
                end if
                WhichExpansionBasis=0
            else
                bisection=(low+up)/2
                if(order>Hd_EBNR(bisection).order) then
                    call bisect(low,bisection)
                else if(order<Hd_EBNR(bisection).order) then
                    call bisect(bisection,up)
                else
                    do i=order,1,-1
                        if(indice(i)/=Hd_EBNR(bisection).indice(i)) exit
                    end do
                    if(i<1) then
                        WhichExpansionBasis=bisection
                    else
                        if(indice(i)>Hd_EBNR(bisection).indice(i)) then
                            call bisect(bisection,up)
                        else
                            call bisect(low,bisection)
                        end if
                    end if
                end if
            end if
        end subroutine bisect
    end function WhichExpansionBasis
!----------------------------------------------------------------------------------------
    subroutine OriginShift(shift)!Transform HdEC according to origin shift from q0 to q1: shift = q1 - q0
        !This is done by:
        !    1, select an Hd expansion basis (Hd_HdEC)
        !    2, under translation, the terms making up the multiplication become var + const, so
        !       we go through all combination of var & const and add the contribution to HdECtemp
        !    3, go to 1 until all Hd_HdEC are done
        !    4, copy HdECtemp to Hd_HdEC
        real*8,dimension(Hd_intdim),intent(in)::shift
        integer::location,nvar,nconst,ivar,iconst,n,i,j
        integer,dimension(Hd_EBNR(1).order)::usevar,indicevar,indiceconst
        real*8::coeff
        type(d2PArray),dimension(Hd_NState,Hd_Nstate)::HdECtemp
        do j=1,Hd_NState!Allocate work space
            do i=j,Hd_NState
                allocate(HdECtemp(i,j).Array(NHdExpansionBasis))
                HdECtemp(i,j).Array=0d0
            end do
        end do

        do n=1,NHdExpansionBasis!Main loop
            if(Hd_EBNR(n).order==0) then!Const term will not change under any transformation
                forall(i=1:Hd_NState,j=1:Hd_NState,i>=j)
                    HdECtemp(i,j).Array(n)=HdECtemp(i,j).Array(n)+Hd_HdEC(i,j).Array(n)
                end forall
            else!Go through all combination in a binary counter manner
                usevar(1:Hd_EBNR(n).order)=0!usevar(i)=1 means using const at i-th position in multiplication
                do while(usevar(Hd_EBNR(n).order)<2)!Done when the counter overflows
                    nvar=sum(usevar(1:Hd_EBNR(n).order))
                    nconst=Hd_EBNR(n).order-nvar
                    ivar=1
                    iconst=1
                    do i=1,Hd_EBNR(n).order
                        if(usevar(i)==1) then
                            indicevar(ivar)=Hd_EBNR(n).indice(i)
                            ivar=ivar+1
                        else
                            indiceconst(iconst)=Hd_EBNR(n).indice(i)
                            iconst=iconst+1
                        end if
                    end do
                    location=WhichExpansionBasis(nvar,indicevar(1:nvar))
                    if(location==0) then
                      print*,indicevar(1:nvar)
                      stop 'Program abort: basis space is not closed under origin shift'
                    end if
                    coeff=1d0
                    do i=1,nconst
                        coeff=coeff*shift(indiceconst(i))
                    end do
                    forall(i=1:Hd_NState,j=1:Hd_NState,i>=j)
                        HdECtemp(i,j).Array(location)=HdECtemp(i,j).Array(location)+coeff*Hd_HdEC(i,j).Array(n)
                    end forall
                    usevar(1)=usevar(1)+1!Add 1 to the binary counter
                    do i=1,Hd_EBNR(n).order-1
                        if(usevar(i)==2) then!Carry
                            usevar(i)=0
                            usevar(i+1)=usevar(i+1)+1
                        end if
                    end do
                end do
            end if
        end do

        Hd_HdEC=HdECtemp
        do j=1,Hd_NState!Clean up
            do i=j,Hd_NState
                deallocate(HdECtemp(i,j).Array)
            end do
        end do
        return
    end subroutine OriginShift
!----------------------------------------------------------------------------------------
!Initialize DiabaticHamiltonian module
!Required: NState & intdim
!Optional: NewHd: (default = false) if true, will return a blank Hd expansion coefficient
subroutine InitializeDiabaticHamiltonian(NState,intdim)
  integer,intent(in)::NState,intdim
  integer::istate,jstate,iorder,i,n
  Hd_NState=NState; Hd_intdim=intdim
  call InitializeExpansionBasisNumberingRule()!NHdExpansionBasis, EBNR
  NHdExpansionCoefficients=Hd_NState*(Hd_NState+1)/2*NHdExpansionBasis!Multiply the number of independent Hd elements
  allocate(Hd_HdEC(Hd_NState,Hd_NState))!HdEC
  do istate=1,Hd_NState
    do jstate=1,Hd_NState
      allocate(Hd_HdEC(jstate,istate)%Array(NHdExpansionBasis))
      Hd_HdEC(jstate,istate)%Array=0.d0
    end do
  end do
  return
end subroutine InitializeDiabaticHamiltonian
!----------------------------------------------------------------------------------------
    !Write Hd expansion coefficient and expansion basis specification to file
    !Optional: FileName: (default = 'Hd.CheckPoint') name of the output file
    subroutine WriteHdExpansionCoefficients(HdEC,FileName)
        type(d2PArray),dimension(Hd_NState,Hd_NState),intent(in)::HdEC
        character*32,optional,intent(in)::FileName
        integer::istate,jstate,i,j
        if(present(FileName)) then; open(unit=99,file=FileName,status='replace')
        else; open(unit=99,file='Hd.CheckPoint',status='replace'); end if
        write(99,'(A28,I2)')'Number of electronic states:',Hd_NState
            do istate=1,Hd_NState
                do jstate=istate,Hd_NState
                    write(99,'(A2,I2,I2)')'Hd',jstate,istate
                    do i=1,NHdExpansionBasis
                        write(99,*)HdEC(jstate,istate)%Array(i)
                        if(Hd_EBNR(i)%order>0) then
                            write(99,'(I5)',advance='no')Hd_EBNR(i)%order
                            do j=1,Hd_EBNR(i)%order-1; write(99,'(I5)',advance='no')Hd_EBNR(i)%indice(j); end do; write(99,'(I5)')Hd_EBNR(i)%indice(Hd_EBNR(i)%order)
                        else
                            write(99,'(I5)')Hd_EBNR(i)%order
                        end if
                    end do
                end do
            end do
        close(99)
    end subroutine WriteHdExpansionCoefficients
!----------------------------------------------------------------------------------------
!The value of n-th expansion basis function at some coordinate q
real*8 function ExpansionBasis(q,n)
    real*8,dimension(Hd_intdim),intent(in)::q
    integer,intent(in)::n
    integer::i
    ExpansionBasis=1d0
    do i=1,Hd_EBNR(n)%order
        ExpansionBasis=ExpansionBasis*q(Hd_EBNR(n)%indice(i))
    end do
end function ExpansionBasis
!------------------------------------------------------------------------------------
!The gradients of (n-th expansion basis function) at some coordinate q
function ExpansionBasisGradient(q,n)
    real*8,dimension(Hd_intdim)::ExpansionBasisGradient
    real*8,dimension(Hd_intdim),intent(in)::q
    integer,intent(in)::n
    integer::m,i,OrderCount
    do m=1,Hd_intdim
        OrderCount=0
        do i=1,Hd_EBNR(n)%order
            if(Hd_EBNR(n)%indice(i)==m) OrderCount=OrderCount+1
        end do
        if(OrderCount>0) then
            ExpansionBasisGradient(m)=dble(OrderCount)*q(m)**(OrderCount-1)
            do i=1,Hd_EBNR(n)%order
                if(Hd_EBNR(n)%indice(i)/=m) &
                ExpansionBasisGradient(m)=ExpansionBasisGradient(m)*q(Hd_EBNR(n)%indice(i))
            end do
        else
            ExpansionBasisGradient(m)=0d0
        end if
    end do
end function ExpansionBasisGradient
!-------------------------------------------------------------------------------
end module DiabaticHamiltonian
!=================================================================================
