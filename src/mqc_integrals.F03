      module MQC_Integrals
!
!     ***********************************************************************
!     ***********************************************************************
!     **                                                                   **
!     **               The Merced Quantum Chemistry Package                **
!     **                            (MQCPack)                              **
!     **                                                                   **
!     ** Written By:                                                       **
!     **    Lee M. Thompson, Xianghai Sheng, Andrew D. Mahler, Dave        **
!     **    Mullally, and Hrant P. Hratchian                               **
!     **                                                                   **
!     **                          Version 25.7.3                           **
!     **                          July 30, 2025                            **
!     **                                                                   **
!     **                                                                   **
!     ** Modules beloning to MQCPack:                                      **
!     **    1. MQC_General                                                 **
!     **    2. MQC_Binary                                                  **
!     **    3. MQC_DataStructures                                          **
!     **    4. MQC_Integrals                                               **
!     **    5. MQC_Algebra                                                 **
!     **    6. MQC_Files                                                   **
!     **    7. MQC_Molecule                                                **
!     **    8. MQC_EST                                                     **
!     **    9. MQC_Gaussian                                                **
!     **                                                                   **
!     ***********************************************************************
!     ***********************************************************************
!
!
!     This module contains objects and procedures for working with and
!     evaluating certain types of integrals over gaussian-type basis functions.
!
!
      use MQC_General
      use MQC_DataStructures
      use iso_fortran_env, only: int64, real64
      implicit none
      integer(kind=int64),parameter,private::iOut=6
      logical::MQC_Integrals_DEBUG=.false.
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     MQC_CGTF
!
!     The MQC_CGTF class contains a contracted Gaussian-type function.
      type MQC_CGTF
        logical,private::initialized=.false.
        integer(kind=int64),private::angularMomentum
        integer(kind=int64),private::nPrimMax,nPrim,nBasis=0
        integer(kind=int64),dimension(:,:),allocatable,public::lArrays
        real(kind=real64),dimension(3)::center
        real(kind=real64),dimension(:),allocatable::coefficient
        real(kind=real64),dimension(:),allocatable::alpha
        real(kind=real64),dimension(:,:),allocatable::primitiveNorm
      Contains
        procedure,public::init                  => MQC_CGTF_init
        procedure,public::print                 => MQC_CGTF_print
        procedure,public::shell2nBasis          => MQC_CGTF_nBasisFunctions
        procedure,public::primitiveSelfOverlap  => MQC_CGTF_primitiveSelfOverlap
      end type MQC_CGTF
!
!     MQC_BasisSet
!
!     The MQC_BasisSet class is an object of a list of MQC_CGTF members to build
!     up a full basis set.
      type MQC_basisSet
        logical::initialized=.false.
        integer(kind=int64)::nBasis,nShells,mxShells
        type(mqc_cgtf),dimension(:),allocatable::shells
      Contains
        procedure,public::init                  => MQC_basisSet_init
        procedure,public::print                 => MQC_basisSet_print
        procedure,public::addShell              => MQC_basisSet_addShell
      end type MQC_BasisSet
!
!
      CONTAINS
!
!--------------------------------------------------------------------------
!                                                                         |
!     OBJECT SET-UP, PRINTING, AND GENERAL MANIPULATION PROCEDURES        |
!                                                                         |
!--------------------------------------------------------------------------
!
!
!PROCEDURE MQC_CGTF_init
      subroutine MQC_CGTF_init(myCGTF,angularMomentum,center,coefficients,  &
        alphas)
!
!     This routine initializes a contracted Gaussian-type basis function object.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64),intent(in)::angularMomentum
      real(kind=real64),dimension(3),intent(in)::center
      real(kind=real64),dimension(:),intent(in)::coefficients,alphas
      real(kind=real64)::totalNorm
      real(kind=real64),dimension(:,:),allocatable::basisIntegrals
      integer(kind=int64)::i,j
!
!     Begin by making sure the object's allocatable arrays aren't already
!     allocated.
!
      if(Allocated(myCGTF%coefficient)) DeAllocate(myCGTF%coefficient)
      if(Allocated(myCGTF%alpha)) DeAllocate(myCGTF%alpha)
      if(Allocated(myCGTF%primitiveNorm)) DeAllocate(myCGTF%primitiveNorm)
!
!     Fill scalars and center array. Also, allocate the allocatable arrays in
!     the object. Also, set the initialization flage to TRUE.
!
      if(angularMomentum.lt.0)  &
        call mqc_error('MQC_CGTF_init: Negative angular momentum sent.')
      myCGTF%angularMomentum = angularMomentum
      myCGTF%nPrim           = Size(coefficients)
      myCGTF%center          = center
      Allocate(myCGTF%coefficient(myCGTF%nPrim),myCGTF%alpha(myCGTF%nPrim))
!
!     Load the primitive coefficients and alphas.
!
      if(Size(coefficients).ne.Size(alphas))  &
        call mqc_error('MQC_CGTF_init: Coefficients and exponents must be same size.')
      myCGTF%coefficient = coefficients
      myCGTF%alpha = alphas
!
!     Allocate primitiveNorm and lArrays arrays. Then fill in both arrays.
!
      myCGTF%nBasis = ((angularMomentum+1)*(angularMomentum+2))/2
      myCGTF%lArrays = MQC_CGTF_extendLArray(myCGTF)
      Allocate(myCGTF%primitiveNorm(myCGTF%nPrim,myCGTF%nBasis))
      do i = 1,myCGTF%nBasis
        do j = 1,myCGTF%nPrim
          myCGTF%primitiveNorm(j,i) = mqc_float(1)/  &
            sqrt(myCGTF%primitiveSelfOverlap(j,myCGTF%lArrays(:,i)))
        endDo
      endDo
      if(MQC_Integrals_DEBUG) call mqc_print(myCGTF%primitiveNorm,iOut,  &
        header='MQC_CGTF_init, primitiveNorm pre-norm:')
!
!     Set the initialization flag to TRUE.
!
      myCGTF%initialized = .true.
!
!     Finalize the normalization constant(s).
!
      call MQC_Overlap_CGFT(myCGTF,myCGTF,basisIntegrals)
      do i = 1,myCGTF%nBasis
        myCGTF%primitiveNorm(:,i) =  &
          myCGTF%primitiveNorm(:,i)/sqrt(basisIntegrals(i,i))
      endDo
      DeAllocate(basisIntegrals)
!
      return
      end subroutine MQC_CGTF_init

!
!PROCEDURE MQC_CGTG_print
      subroutine MQC_CGTF_print(myCGTF,iOut,header)
!
!     This routine prints the data in a MQC_CGTF object.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64),intent(in)::iOut
      Character(Len=*),Intent(In),optional::header
      integer(kind=int64)::i
!
 1000 format(1x,'l=',i4,' shell at (',f8.4,',',f8.4,',',f8.4,'). NPrim=',i4)
 2000 format(3x,'Primitive Exponents and Coefficients.')
 2100 format(6x,f15.8,4x,f15.8)
!
      if(PRESENT(header)) write(iOut,'(/,1x,A)') TRIM(header)
      write(iOut,1000) myCGTF%angularMomentum,myCGTF%center,myCGTF%nPrim
      write(iOut,2000)
      do i = 1,myCGTF%nPrim
        write(iOut,2100) myCGTF%alpha(i),myCGTF%coefficient(i)
      endDo
!
      return
      end subroutine MQC_CGTF_print

!
!PROCEDURE MQC_basisSet_init
      subroutine MQC_basisSet_init(myBasisSet,maxShells)
!
!     This routine initializes a basis set of contracted Gaussian-type basis
!     functions.
!
      implicit none
      class(MQC_basisSet)::myBasisSet
      integer(kind=int64),intent(in)::maxShells
!
!     Check if the init flag is TRUE. If it is, clear things out.
!
      if(myBasisSet%initialized) then
        myBasisSet%initialized=.false.
        myBasisSet%nbasis=0
        myBasisSet%nShells=0
        myBasisSet%mxShells=0
        if(Allocated(myBasisSet%shells)) DeAllocate(myBasisSet%shells)
      endIf
!
!     Set nBasis, nShells, and mxShells. Then, allocate the basis function
!     array.
!
      myBasisSet%nBasis = 0
      myBasisSet%nShells = 0
      myBasisSet%mxShells = maxShells
      Allocate(myBasisSet%shells(maxShells))
!
!     Set the initialization flag to TRUE.
!
      myBasisSet%initialized = .true.
!
      return
      end subroutine MQC_basisSet_init

!
!PROCEDURE MQC_basisSet_addShell
      subroutine MQC_basisSet_addShell(myBasisSet,angularMomentum,  &
        center,coefficients,alphas)
!
!     This routine adds a shell to the object myBasisSet.
!
!
!     H. P. Hratchian, 2025
!
      implicit none
      class(MQC_basisSet)::myBasisSet
      integer(kind=int64),intent(in)::angularMomentum
      real(kind=real64),dimension(3),intent(in)::center
      real(kind=real64),dimension(:)::coefficients,alphas
!
!     Check if the init flag is FALSE.
!
      if(.not.myBasisSet%initialized)  &
        call mqc_error('MQC_basisSet_addShell: basis set object not yet initialized.')
!
!     Ensure that adding this shell won't take us past the maximum number of
!     allowed shells. For now, in such cases we fail out.
!
      if(myBasisSet%nShells.ge.myBasisSet%mxShells)  &
        call mqc_error('MQC_basisSet_addShell: cannot add a new shell.')
!
!     Add the next shell to myBasisSet.
!
      myBasisSet%nShells = myBasisSet%nShells+1
      call MQC_CGTF_init(myBasisSet%shells(myBasisSet%nShells),  &
        angularMomentum,center,coefficients,alphas)
      myBasisSet%nBasis = myBasisSet%nBasis+  &
        myBasisSet%shells(myBasisSet%nShells)%shell2nBasis()
!
      return
      end subroutine MQC_basisSet_addShell

!
!PROCEDURE MQC_basisSet_print
      subroutine MQC_basisSet_print(myBasisSet,iOut)
!
!     This routine prints the data in a MQC_basisSet object.
!
      implicit none
      class(MQC_basisSet)::myBasisSet
      integer(kind=int64),intent(in)::iOut
      integer(kind=int64)::i
!
 1000 format(1x,'Shell number ',i5)
!
      do i = 1,myBasisSet%nShells
        write(iOut,1000) i
        call myBasisSet%shells(i)%print(iOut)
      endDo
!
      return
      end subroutine MQC_basisSet_print

!
!PROCEDURE MQC_Integrals_SetDEBUG
      subroutine MQC_Integrals_SetDEBUG(setDebugValue)
!
!     This subroutine is called to set a Module-Wide DEBUG flag to .TRUE. or
!     .FALSE.
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      logical,intent(IN)::setDebugValue
!
      MQC_Integrals_DEBUG = setDebugValue
!
      return
      end subroutine MQC_Integrals_SetDEBUG
!
!--------------------------------------------------------------------------
!                                                                         |
!     HIGH LEVEL PROCEDURES                                               |
!                                                                         |
!--------------------------------------------------------------------------
!
!PROCEDURE basisSetValuesList
      function basisSetValuesList(myBasisSet,xyz) result(valuesList)
!
!     This function computes the values of each basis function at a position
!     xyz.
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      class(MQC_basisSet),intent(in)::myBasisSet
      real(kind=real64),dimension(3),intent(in)::xyz
      real(kind=real64),dimension(:),allocatable::valuesList,  &
        valuesListPartial
      integer(kind=int64)::i,iStart,iEnd
!
!     Loop over shells in the basis set object and fill the values list in
!     blocks.
!
      Allocate(valuesList(myBasisSet%nBasis))
      iStart = 1
      do i = 1,myBasisSet%nShells
        call MQC_Value_CGFT(myBasisSet%shells(i),xyz,valuesListPartial)
        iEnd = iStart+Size(valuesListPartial)-1
        valuesList(iStart:iEnd) = valuesListPartial
        iStart = iEnd + 1
      endDo
!
      return
      end function basisSetValuesList

!
!PROCEDURE basisSetValuesList1
      subroutine basisSetValuesList1(myBasisSet,xyz,valuesList)
!
!     This routine computes the values of each basis function at a position
!     xyz. This is the same as Function basisSetValuesList, except that
!     valuesList comes here allocated and it avoids memory allocation costs.
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      class(MQC_basisSet),intent(in)::myBasisSet
      real(kind=real64),dimension(3),intent(in)::xyz
      real(kind=real64),dimension(:),allocatable::valuesList,  &
        valuesListPartial
      integer(kind=int64)::i,iStart,iEnd
!
!     Loop over shells in the basis set object and fill the values list in
!     blocks.
!
      if(.not.Allocated(valuesList)) Allocate(valuesList(myBasisSet%nBasis))
      iStart = 1
      do i = 1,myBasisSet%nShells
        call MQC_Value_CGFT(myBasisSet%shells(i),xyz,valuesListPartial)
        iEnd = iStart+Size(valuesListPartial)-1
        valuesList(iStart:iEnd) = valuesListPartial
        iStart = iEnd + 1
      endDo
!
      return
      end subroutine basisSetValuesList1

!
!PROCEDURE basisSetOverlapMatrix
      function basisSetOverlapMatrix(myBasisSet) result(overlapMatrix)
!
!     This function computes the overlap matrix based on basisSet in the bra and
!     ket.
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      class(MQC_basisSet)::myBasisSet
      real(kind=real64),dimension(:,:),allocatable::integralBlock,  &
        overlapMatrix
      integer(kind=int64)::i,j,iRowStart,iColumnStart,iRowEnd,iColumnEnd
!
!     Loop over shells in the basis set object and fill the overlap matrix in
!     blocks.
!
      Allocate(overlapMatrix(myBasisSet%nBasis,myBasisSet%nBasis))
      iRowStart = 1
      do i = 1,myBasisSet%nShells
        iColumnStart = 1
        do j = 1,myBasisSet%nShells
          call MQC_Overlap_CGFT(myBasisSet%shells(i),  &
            myBasisSet%shells(j),integralBlock)
          iRowEnd = iRowStart+Size(integralBlock,1)-1
          iColumnEnd = iColumnStart+Size(integralBlock,2)-1
          overlapMatrix(iRowStart:iRowEnd,iColumnStart:iColumnEnd) =  &
            integralBlock
          iColumnStart = iColumnEnd + 1
        endDo
        iRowStart = iRowEnd + 1
      endDo
!
      return
      end function basisSetOverlapMatrix

!
!PROCEDURE LinkedList_Return_CGTF_Value
      subroutine LinkedList_Return_CGTF_Value(LL_Cur,Return_Value)
!
!     This routine returns the value stored in a linked list node
!     when the value is of type MQC_CGTF.
!
!     H. P. Hratchian, 2025.
!
      implicit none
      type(MQC_LinkedList),pointer::LL_Cur
      type(MQC_CGTF),intent(out)::Return_Value
!
!     Check that the current node is associated and contains a CGTF object.
!
      if(Associated(LL_Cur)) then
        if(Allocated(LL_Cur%Val)) then
          select type (Val => LL_Cur%Val)
          type Is (MQC_CGTF)
            Return_Value = Val
          class default
            call mqc_error('Data types not compatible in LinkedList_Return_CGTF_Value')
          end Select
        endIf
      endIf
!
      return
      end subroutine LinkedList_Return_CGTF_Value

!
!PROCEDURE MQC_CGTG_nBasisFunctions
      function MQC_CGTF_nBasisFunctions(myCGTF)  result(nBasis)
!
!     This function returns the number of basis functions associated with the
!     shell described by object myCGTF.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64)::nBasis
!
      nBasis = myCGTF%nBasis
!
      return
      end function MQC_CGTF_nBasisFunctions

!
!PROCEDURE MQC_CGTF_extendLArray
      function MQC_CGTF_extendLArray(myCGTF,ordering)  result(lArrays)
!
!     This function builds the list of angular momentum arrays, which give
!     (lx,ly,lz) triples, for the contracted Gaussian type function myCGTF. The
!     returned array, lArrays, is allocated (3,nBasisFunctions).
!
!     There are multiple ways to do the ordering. The OPTIONAL dummy argument
!     <ordering> is used to chose an ordering algorithm. The available options
!     are:
!           lexicographic     This ordering is the standard used in libint,
!                             PySCF, and MolPro. As an example, Cartesian d
!                             functions are ordered as (0,0,2), (0,1,1),
!                             (0,2,0), (1,0,1), (1,1,0), and (2,0,0).
!           gaussian          This ordering is the standard used in Gaussian. As
!                             an example, Cartesian d functions are ordered as
!                             (2,0,0), (0,2,0), (0,0,2), (1,1,0), (1,0,1), and
!                             (0,1,1). This is the default.
!
!
!     - H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      character(len=64),optional::ordering
      integer(kind=int64),dimension(:,:),allocatable::lArrays
      integer(kind=int64)::nBasisFunctions,i,lx,ly,lz
      character(len=64)::myOrdering
!
!     Start by setting myOrdering.
!
      if(PRESENT(ordering)) then
        call String_Change_Case(ordering,'l',myOrdering)
      else
        myOrdering = 'gaussian'
      endIf
!
!     Allocate lArrays and then fill it.
!
      nBasisFunctions = myCGTF%shell2nBasis()
      Allocate(lArrays(3,nBasisFunctions))
      lArrays = -1
      select case(myOrdering)
      case('lexicographic')
        i = 1
        do lx = 0,myCGTF%angularMomentum
          do ly = 0,myCGTF%angularMomentum-lx
            lz = myCGTF%angularMomentum-lx-ly
            lArrays(1,i) = lx
            lArrays(2,i) = ly
            lArrays(3,i) = lz
            i = i+1
          endDo
        endDo
      case('gaussian')
        if(myCGTF%angularMomentum.eq.0) then
          lArrays(1,1) = 0
          lArrays(2,1) = 0
          lArrays(3,1) = 0
        elseIf(myCGTF%angularMomentum.eq.1) then
          lArrays(1,1) = 1
          lArrays(2,1) = 0
          lArrays(3,1) = 0
          lArrays(1,2) = 0
          lArrays(2,2) = 1
          lArrays(3,2) = 0
          lArrays(1,3) = 0
          lArrays(2,3) = 0
          lArrays(3,3) = 1
        elseIf(myCGTF%angularMomentum.eq.2) then
          lArrays(1,1) = 2
          lArrays(2,1) = 0
          lArrays(3,1) = 0
          lArrays(1,2) = 0
          lArrays(2,2) = 2
          lArrays(3,2) = 0
          lArrays(1,3) = 0
          lArrays(2,3) = 0
          lArrays(3,3) = 2
          lArrays(1,4) = 1
          lArrays(2,4) = 1
          lArrays(3,4) = 0
          lArrays(1,5) = 1
          lArrays(2,5) = 0
          lArrays(3,5) = 1
          lArrays(1,6) = 0
          lArrays(2,6) = 1
          lArrays(3,6) = 1
        elseIf(myCGTF%angularMomentum.eq.3) then
          lArrays(1,1) = 3
          lArrays(2,1) = 0
          lArrays(3,1) = 0
          lArrays(1,2) = 0
          lArrays(2,2) = 3
          lArrays(3,2) = 0
          lArrays(1,3) = 0
          lArrays(2,3) = 0
          lArrays(3,3) = 3
          lArrays(1,4) = 1
          lArrays(2,4) = 2
          lArrays(3,4) = 0
          lArrays(1,5) = 2
          lArrays(2,5) = 1
          lArrays(3,5) = 0
          lArrays(1,6) = 2
          lArrays(2,6) = 0
          lArrays(3,6) = 1
          lArrays(1,7) = 1
          lArrays(2,7) = 0
          lArrays(3,7) = 2
          lArrays(1,8) = 0
          lArrays(2,8) = 1
          lArrays(3,8) = 2
          lArrays(1,9) = 0
          lArrays(2,9) = 2
          lArrays(3,9) = 1
          lArrays(1,10) = 1
          lArrays(2,10) = 1
          lArrays(3,10) = 1
        elseIf(myCGTF%angularMomentum.gt.3) then
          call mqc_error('MQC_CGTF_extendLArray: Gaussian ordering NYI for L>3.')
        endIf
      case default
        call mqc_error('MQC_CGTF_extendLArray: Invalid ordering sent.')
      endSelect
!
      return
      end function MQC_CGTF_extendLArray

!
!PROCEDURE MQC_CGTG_primitiveSelfOverlap
      function MQC_CGTF_primitiveSelfOverlap(myCGTF,i,lVector)  result(selfOverlap)
!
!     This routine evaluates and returns the self-overlap of the i-th primitive
!     basis function in the shell of myCGTF with angular momentum vector
!     lVector.
!
!     -H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      integer(kind=int64),intent(in)::i
      integer(kind=int64),dimension(3),intent(in)::lVector
      real(kind=real64)::selfOverlap
      real(kind=real64)::fourAlpha,rootPiFourA3
!
!     Check that lVector is valid.
!
      if(Sum(lVector).ne.myCGTF%angularMomentum) call mqc_error('MQC_CGTF_selfOverlap: Invalid lVector sent.')
!
!     Evaluate the self-overlap.
!
      fourAlpha = 4*myCGTF%alpha(i)
      rootPiFourA3 = (SQRT(Pi/(2*myCGTF%alpha(i))))**3
      if(lVector(1).eq.0) then
        selfOverlap = mqc_float(1)
      else
        selfOverlap = mqc_doubleFactorial((2*lVector(1)-1))/(fourAlpha**lVector(1))
      endIf
      if(lVector(2).gt.0) selfOverlap = selfOverlap*(mqc_doubleFactorial((2*lVector(2)-1))/(fourAlpha**lVector(2)))
      if(lVector(3).gt.0) selfOverlap = selfOverlap*(mqc_doubleFactorial((2*lVector(3)-1))/(fourAlpha**lVector(3)))
      selfOverlap = selfOverlap*rootPiFourA3
!
      return
      end function MQC_CGTF_primitiveSelfOverlap

!
!--------------------------------------------------------------------------
!                                                                         |
!     INTERMEDIATE AND LOW LEVEL ROUTINES                                 |
!                                                                         |
!--------------------------------------------------------------------------
!
!PROCEDURE MQC_Overlap_Distribution_Constants
      subroutine MQC_Overlap_Distribution_Primitive_XYZ_Constants(Ga,Gb, &
        iPrimA,iPrimB,mu,p,xAB,xPA,xPB)
!
!     This routine is sent two primitive contracted Gaussian basis functions and
!     returns the Gaussian Overlap Distribution intermediate constants (mu, p,
!     xAB, xPA, and xPB) for the product of primitive number iPrimA from Ga and
!     iPrimB from Gb.
!
!
!     - H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::Ga,Gb
      integer(kind=int64),intent(in)::iPrimA,iPrimB
      real(kind=real64),intent(out)::mu,p
      real(kind=real64),dimension(3),intent(out)::xAB,xPA,xPB
!
!     Do the work...
!
      mu = ga%alpha(iPrimA)*gb%alpha(iPrimB)/  &
        (ga%alpha(iPrimA)+gb%alpha(iPrimB))
      p = ga%alpha(iPrimA)+gb%alpha(iPrimB)
      xAB = ga%center-gb%center
      xPA = -xAB*gb%alpha(iPrimB)/p
      xPB =  xAB*ga%alpha(iPrimA)/p
!
      return
      end subroutine MQC_Overlap_Distribution_Primitive_XYZ_Constants

!
!PROCEDURE MQC_Value_CGTF
      subroutine MQC_Value_CGFT(bfShell,xyz,basisValues)
!
!     This subroutine calculates the values of the basis functions associated
!     with contracted shell bfShell. Dummy argument xyz gives the position where
!     the function is to be evaluated.
!
!     - H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::bfShell
      real(kind=real64),dimension(3),intent(in)::xyz
      real(kind=real64),dimension(:),allocatable,intent(out)::basisValues
      integer(kind=int64)::i,j
      real(kind=real64)::radialValue
      real(kind=real64),dimension(:,:),allocatable::primCoeffs,  &
        primBasisValues
!
!     Calculate nBasisPairs and nPrimPairs. Then allocate some arrays.
!
      if(MQC_Integrals_DEBUG) then
        write(iOut,*)' In MQC_Value_CGFT, basis shell:'
        call MQC_CGTF_print(bfShell,iOut)
      endIf
      if(Allocated(basisValues)) DeAllocate(basisValues)
      Allocate(basisValues(bfShell%nBasis),  &
        primBasisValues(bfShell%nPrim,bfShell%nBasis),  &
        primCoeffs(bfShell%nPrim,bfShell%nBasis))
!
!     Work through the list of basis functions this shell expands to and build
!     up the function values at position xyz.
!
      do i = 1,bfShell%nBasis
        do j = 1,bfShell%nPrim
          primCoeffs(j,i) =  &
            bfShell%primitiveNorm(j,i)*  &
            bfShell%coefficient(j)
          radialValue = MQC_Value_Primitive_Radial(xyz,  &
            bfShell%center,bfShell%alpha(j))
          primBasisValues(j,i) = radialValue*  &
            MQC_Value_Primitive_Angular(xyz,bfShell%center,  &
            bfShell%alpha(j),bfShell%lArrays(:,i))
        endDo
      endDo
      do i = 1,bfShell%nBasis
        basisValues(i) =  &
          dot_product(primCoeffs(:,i),primBasisValues(:,i))
      endDo
!
      return
      end subroutine MQC_Value_CGFT

!
!PROCEDURE MQC_Value_Primitive_Radial
      function MQC_Value_Primitive_Radial(xyz,center,alpha)  &
        Result(primitiveValue)
!
!     This function evaluates the radial part of the function value of a
!     primitive Cartesian gaussian function at position xyz.
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      real(kind=real64),dimension(3),intent(in)::xyz,center
      real(kind=real64),intent(in)::alpha
      real(kind=real64)::primitiveValue
      real(kind=real64)::distanceSq
!
      distanceSq = dot_product(xyz-center,xyz-center)
      primitiveValue = exp(-alpha*distanceSq)
!
      return
      end function MQC_Value_Primitive_Radial

!
!PROCEDURE MQC_Value_Primitive_Angular
      function MQC_Value_Primitive_Angular(xyz,center,alpha,lArray)  &
        Result(primitiveValue)
!
!     This function evaluates the angular part of the function value of a
!     primitive Cartesian gaussian function at position xyz.
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      real(kind=real64),dimension(3),intent(in)::xyz,center
      real(kind=real64),intent(in)::alpha
      integer(kind=int64),dimension(3),intent(in)::lArray
      real(kind=real64)::primitiveValue
!
      primitiveValue = ((xyz(1)-center(1))**lArray(1))*  &
        ((xyz(2)-center(2))**lArray(2))*  &
        ((xyz(3)-center(3))**lArray(3))
!
      return
      end function MQC_Value_Primitive_Angular

!
!PROCEDURE MQC_Overlap_CGTF
      subroutine MQC_Overlap_CGFT(bra,ket,basisIntegrals)
!
!     This subroutine calculates an overlap integral between two contracted
!     Gaussian type functions. These are sent as MQC_CGTF objects. The output is
!     basisIntegrals, which is a matrix of the overlap integrals over the
!     expansion of basis functions that come from the contracted Gaussians in
!     bra and ket. At present, all overlap integrals are solved using the
!     Obara-Saika recursion relation.
!
!     - H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::bra,ket
      real(kind=real64),dimension(:,:),allocatable,intent(out)::basisIntegrals
      integer(kind=int64)::i,j,k,ixyz,iBasisBra,iBasisKet,nBasisBra,  &
        nBasisKet,nBasisPairs,nPrimPairs
      integer(kind=int64),dimension(:,:),allocatable::indexBasisPairs,  &
        indexPrimPairs
      real(kind=real64)::mu,p
      real(kind=real64),dimension(3)::xAB,xPA,xPB
      real(kind=real64),dimension(:,:),allocatable::tmpSij,primCoeffs,  &
        primBasisIntsTot
      real(kind=real64),dimension(:,:,:),allocatable::primBasisInts
      logical,dimension(:,:),allocatable::haveSij
!
 1000 format(1x,'Basis Pair: (',i2,' | ',i2,') = ',f8.5)
 1010 format(3x,'iPrim      IntX          IntY          IntZ        IntTot')
 1020 format(4x,i3,4x,f8.5,6x,f8.5,6x,f8.5,6x,f8.5)
 1021 format(4x,i3,4x,i3,4x,f8.5,6x,f8.5,6x,f8.5,6x,f8.5)
 1030 format(/,/)
!
!     Calculate nBasisPairs and nPrimPairs. Then allocate some arrays.
!
      if(MQC_Integrals_DEBUG) then
        write(iOut,*)' In MQC_Overlap_CGFT, bra and ket:'
        call MQC_CGTF_print(bra,iOut)
        call MQC_CGTF_print(ket,iOut)
      endIf
      nBasisBra = bra%nBasis
      nBasisKet = ket%nBasis
      nBasisPairs = nBasisBra*nBasisKet
      nPrimPairs  = bra%nPrim*ket%nPrim
      if(Allocated(basisIntegrals)) DeAllocate(basisIntegrals)
      Allocate(basisIntegrals(nBasisBra,nBasisKet))
      Allocate(indexBasisPairs(2,nBasisPairs),indexPrimPairs(2,nPrimPairs))
      Allocate(primBasisInts(nPrimPairs,3,nBasisPairs),  &
        primBasisIntsTot(nPrimPairs,nBasisPairs),  &
        primCoeffs(nPrimPairs,nBasisPairs))
!
!     Fill indexBasisPairs and indexPrimPairs.
!
      k = 0
      do j = 1,nBasisKet
        do i = 1,nBasisBra
          k = k+1
          indexBasisPairs(1,k) = i
          indexBasisPairs(2,k) = j
        endDo
      endDo
      k = 0
      do j = 1,ket%nPrim
        do i = 1,bra%nPrim
          k = k+1
          indexPrimPairs(1,k) = i
          indexPrimPairs(2,k) = j
        endDo
      endDo
!
!     Work through the basis pairs and prim pairs to build up all of the overlap
!     integrals. We do these integrals in x, y, and z directions separatly and
!     then contract.
!
      Allocate(tmpSij(bra%angularMomentum+1,ket%angularMomentum+1),  &
        haveSij(bra%angularMomentum+1,ket%angularMomentum+1))
      do i = 1,nBasisPairs
        do j = 1,nPrimPairs
          primCoeffs(j,i) =  &
            bra%primitiveNorm(indexPrimPairs(1,j),indexBasisPairs(1,i))*  &
            bra%coefficient(indexPrimPairs(1,j))*  &
            ket%primitiveNorm(indexPrimPairs(2,j),indexBasisPairs(2,i))*  &
            ket%coefficient(indexPrimPairs(2,j))
          call MQC_Overlap_Distribution_Primitive_XYZ_Constants(bra,ket, &
            indexPrimPairs(1,j),indexPrimPairs(2,j),mu,p,xAB,xPA,xPB)
          do ixyz = 1,3
            tmpSij  = mqc_float(-1)
            haveSij = .false.
            call MQC_Overlap_Primitive_XYZ_OS(  &
              bra%lArrays(ixyz,indexBasisPairs(1,i)),  &
              ket%lArrays(ixyz,indexBasisPairs(2,i)),  &
              mu,p,xAB(ixyz),xPA(ixyz),xPB(ixyz),tmpSij,haveSij)
            primBasisInts(j,ixyz,i) = tmpSij(  &
              bra%lArrays(ixyz,indexBasisPairs(1,i))+1,  &
              ket%lArrays(ixyz,indexBasisPairs(2,i))+1)
          endDo
          primBasisIntsTot(j,i) = primBasisInts(j,1,i)*  &
            primBasisInts(j,2,i)*primBasisInts(j,3,i)
        endDo
      endDo
      do i = 1,nBasisPairs
        basisIntegrals(indexBasisPairs(1,i),indexBasisPairs(2,i)) =  &
          dot_product(primCoeffs(:,i),primBasisIntsTot(:,i))
        if(MQC_Integrals_DEBUG) then
          write(iOut,1000) indexBasisPairs(1,i),indexBasisPairs(2,i),  &
            dot_product(primCoeffs(:,i),primBasisIntsTot(:,i))
          write(iOut,1010)
          do j = 1,nPrimPairs
            write(iOut,1020) j,primBasisInts(j,1,i),primBasisInts(j,2,i),  &
              primBasisInts(j,3,i),primBasisIntsTot(j,i)
          endDo
        endIf
      endDo
      if(MQC_Integrals_DEBUG) write(iOut,1030)
!
      return
      end subroutine MQC_Overlap_CGFT

!
!PROCEDURE MQC_Overlap_Primitive_XYZ_OS
      recursive subroutine MQC_Overlap_Primitive_XYZ_OS(ixyz,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
!
!     This recursive subroutine solves a 1-D Cartesian overlap between two
!     Cartesian Gaussian-type primitive functions using the Obara-Saika
!     recurrence relation.
!
!
!     -H. P. Hratchian, 2025.
!
      implicit none
      integer(kind=int64)::ixyz,jxyz
      real(kind=real64)::mu,p,XAB,XPA,XPB
      real(kind=real64),dimension(:,:)::Sij
      logical,dimension(:,:)::haveSij
!
!     Look at ixyz and jxyz to determine what we need to do here. There are 9
!     unique cases.
!
!     0,0
      if(ixyz.eq.0.and.jxyz.eq.0) then
        if(haveSij(1,1)) return
        Sij(1,1) = SQRT(pi/p)*exp(-mu*XAB**2)
        haveSij(1,1) = .true.
!
!     0,1
      elseIf(ixyz.eq.0.and.jxyz.eq.1) then
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(1,2) = XPB*Sij(1,1)
        haveSij(1,2) = .true.
!
!     0,>=2
      elseIf(ixyz.eq.0.and.jxyz.ge.2) then
        if(.not.haveSij(1,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,jxyz-1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-2,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(1,jxyz+1) = XPB*Sij(1,jxyz)+Sij(1,jxyz-1)*(jxyz-1)/(mqc_float(2)*p)
        haveSij(1,jxyz+1) = .true.
!
!     1,0
      elseIf(ixyz.eq.1.and.jxyz.eq.0) then
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,1) = XPA*Sij(1,1)
        haveSij(2,1) = .true.
!
!     >=2,0
      elseIf(ixyz.ge.2.and.jxyz.eq.0) then
        if(.not.haveSij(ixyz,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,1) = XPA*Sij(ixyz,1)+Sij(ixyz-1,1)*(ixyz-1)/(mqc_float(2)*p)
        haveSij(ixyz+1,1) = .true.
!
!     1,1  ==>  0,1
      elseIf(ixyz.eq.1.and.jxyz.eq.1) then
        if(.not.haveSij(1,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,2) = XPA*Sij(1,2)+Sij(1,1)/(mqc_float(2)*p)
        haveSij(2,2) = .true.
!
!     1,>=2
      elseIf(ixyz.eq.1.and.jxyz.ge.2) then
        if(.not.haveSij(1,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,jxyz+1) = XPA*Sij(1,jxyz+1)+jxyz*Sij(1,jxyz)/(mqc_float(2)*p)
        haveSij(2,jxyz+1) = .true.
!
!     >=2,1
      elseIf(ixyz.ge.2.and.jxyz.eq.1) then
        if(.not.haveSij(ixyz,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,2) = XPA*Sij(ixyz,2)+(ixyz*Sij(ixyz-1,2)+Sij(ixyz,1))/(mqc_float(2)*p)
        haveSij(ixyz+1,2) = .true.
!
!     >=2,>=2
      elseIf(ixyz.ge.2.and.jxyz.ge.2) then
        if(.not.haveSij(ixyz,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,jxyz+1) = XPA*Sij(ixyz,jxyz+1)+((ixyz-1)*Sij(ixyz-1,jxyz+1)+(jxyz)*Sij(ixyz,jxyz))/(mqc_float(2)*p)
        haveSij(ixyz+1,jxyz+1) = .true.
!
!     Uknown combination...fail out.
      else
        call mqc_error('MQC_Overlap_Primitive_XYZ_OS: Invalid i and j.')
      endIf
!
!     If DEBUG print is on, print out haveSij and Sij.
!
      if(MQC_Integrals_DEBUG) then
        write(iOut,*)
        call mqc_print(haveSij,iOut,header='haveSij')
        call mqc_print(Sij,iOut,header='Sij')
        write(iOut,*)
      endIf
!
      return
      end subroutine MQC_Overlap_Primitive_XYZ_OS
!
!
      end module MQC_Integrals
