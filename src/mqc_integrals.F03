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
!     **                          Version 26.9.2                           **
!     **                         September 3, 2026                         **
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
      integer(kind=int64),parameter,public::MQC_CGTF_CARTESIAN=0_int64
      integer(kind=int64),parameter,public::MQC_CGTF_REAL_PURE=1_int64
      private::MQC_CGTF_buildCartesianToRealPure
      private::MQC_CGTF_cartesianToRealPureCoefficient
      private::MQC_CGTF_paritySign
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
        integer(kind=int64),private::angularRepresentation=  &
          MQC_CGTF_CARTESIAN
        integer(kind=int64),private::nPrimMax,nPrim,nBasis=0,nCartesian=0
        integer(kind=int64),dimension(:,:),allocatable,public::lArrays
        real(kind=real64),dimension(3)::center
        real(kind=real64),dimension(:),allocatable::coefficient
        real(kind=real64),dimension(:),allocatable::alpha
        real(kind=real64),dimension(:,:),allocatable::primitiveNorm
        real(kind=real64),dimension(:,:),allocatable,private::  &
          cartesianToBasis
      Contains
        procedure,public::init                  => MQC_CGTF_init
        procedure,public::print                 => MQC_CGTF_print
        procedure,public::shell2nBasis          => MQC_CGTF_nBasisFunctions
        procedure,public::shell2nCartesian      =>  &
          MQC_CGTF_nCartesianFunctions
        procedure,public::getAngularRepresentation =>  &
          MQC_CGTF_getAngularRepresentation
        procedure,public::getCartesianToBasis   =>  &
          MQC_CGTF_getCartesianToBasis
        procedure,public::primitiveSelfOverlap  => MQC_CGTF_primitiveSelfOverlap
        procedure,public::primitiveValues       => MQC_Value_Primitive_CGTF
      end type MQC_CGTF
!
!     MQC_gtoBasisSet
!
!     The MQC_gtoBasisSet class is an object of a list of MQC_CGTF members to
!     build up a full Gaussian-type-orbital basis set.
      type MQC_gtoBasisSet
        logical::initialized=.false.
        integer(kind=int64)::nBasis,nCartesian,nShells,mxShells
        type(mqc_cgtf),dimension(:),allocatable::shells
      Contains
        procedure,public::init                  => MQC_gtoBasisSet_init
        procedure,public::print                 => MQC_gtoBasisSet_print
        procedure,public::addShell              => MQC_gtoBasisSet_addShell
      end type MQC_gtoBasisSet
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
        alphas,angularRepresentation)
!
!     This routine initializes a contracted Gaussian-type basis function object.
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64),intent(in)::angularMomentum
      integer(kind=int64),intent(in),optional::angularRepresentation
      real(kind=real64),dimension(3),intent(in)::center
      real(kind=real64),dimension(:),intent(in)::coefficients,alphas
      real(kind=real64)::totalNorm
      real(kind=real64),dimension(:,:),allocatable::basisIntegrals
      integer(kind=int64)::i,j,myAngularRepresentation
!
!     Begin by making sure the object's allocatable arrays aren't already
!     allocated.
!
      if(Allocated(myCGTF%coefficient)) DeAllocate(myCGTF%coefficient)
      if(Allocated(myCGTF%alpha)) DeAllocate(myCGTF%alpha)
      if(Allocated(myCGTF%primitiveNorm)) DeAllocate(myCGTF%primitiveNorm)
      if(Allocated(myCGTF%cartesianToBasis))  &
        DeAllocate(myCGTF%cartesianToBasis)
!
!     Fill scalars and center array. Also, allocate the allocatable arrays in
!     the object. Also, set the initialization flage to TRUE.
!
      if(angularMomentum.lt.0)  &
        call mqc_error('MQC_CGTF_init: Negative angular momentum sent.')
      myAngularRepresentation = MQC_CGTF_CARTESIAN
      if(Present(angularRepresentation))  &
        myAngularRepresentation = angularRepresentation
      if(myAngularRepresentation.ne.MQC_CGTF_CARTESIAN.and.  &
        myAngularRepresentation.ne.MQC_CGTF_REAL_PURE)  &
        call mqc_error('MQC_CGTF_init: Invalid angular representation sent.')
      myCGTF%angularMomentum = angularMomentum
      myCGTF%angularRepresentation = MQC_CGTF_CARTESIAN
      myCGTF%nPrim           = Size(coefficients)
      myCGTF%nCartesian      =  &
        ((angularMomentum+1)*(angularMomentum+2))/2
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
      myCGTF%nBasis = myCGTF%nCartesian
      myCGTF%lArrays = MQC_CGTF_extendLArray(myCGTF)
      Allocate(myCGTF%primitiveNorm(myCGTF%nPrim,myCGTF%nCartesian))
      do i = 1,myCGTF%nCartesian
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
!     Set the requested/exposed shell dimension and its transformation from
!     the complete Cartesian working shell. Real-pure S and P shells have the
!     same components as their Cartesian forms. Higher real-pure shells use
!     the Gaussian-order normalized Cartesian-to-real-pure transformation.
!
      myCGTF%angularRepresentation = myAngularRepresentation
      if(myCGTF%angularRepresentation.eq.MQC_CGTF_REAL_PURE)  &
        myCGTF%nBasis = 2*angularMomentum+1
      if(myCGTF%angularRepresentation.eq.MQC_CGTF_REAL_PURE.and.  &
        angularMomentum.gt.1) then
        myCGTF%cartesianToBasis =  &
          MQC_CGTF_buildCartesianToRealPure(myCGTF)
      else
        Allocate(myCGTF%cartesianToBasis(myCGTF%nCartesian,myCGTF%nBasis))
        myCGTF%cartesianToBasis = mqc_float(0)
        do i = 1,myCGTF%nBasis
          myCGTF%cartesianToBasis(i,i) = mqc_float(1)
        endDo
      endIf
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
!PROCEDURE MQC_gtoBasisSet_init
      subroutine MQC_gtoBasisSet_init(myBasisSet,maxShells)
!
!     This routine initializes a basis set of contracted Gaussian-type basis
!     functions.
!
      implicit none
      class(MQC_gtoBasisSet)::myBasisSet
      integer(kind=int64),intent(in)::maxShells
!
!     Check if the init flag is TRUE. If it is, clear things out.
!
      if(myBasisSet%initialized) then
        myBasisSet%initialized=.false.
        myBasisSet%nbasis=0
        myBasisSet%nCartesian=0
        myBasisSet%nShells=0
        myBasisSet%mxShells=0
        if(Allocated(myBasisSet%shells)) DeAllocate(myBasisSet%shells)
      endIf
!
!     Set nBasis, nShells, and mxShells. Then, allocate the basis function
!     array.
!
      myBasisSet%nBasis = 0
      myBasisSet%nCartesian = 0
      myBasisSet%nShells = 0
      myBasisSet%mxShells = maxShells
      Allocate(myBasisSet%shells(maxShells))
!
!     Set the initialization flag to TRUE.
!
      myBasisSet%initialized = .true.
!
      return
      end subroutine MQC_gtoBasisSet_init

!
!PROCEDURE MQC_gtoBasisSet_addShell
      subroutine MQC_gtoBasisSet_addShell(myBasisSet,angularMomentum,  &
        center,coefficients,alphas,angularRepresentation)
!
!     This routine adds a shell to the object myBasisSet.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      class(MQC_gtoBasisSet)::myBasisSet
      integer(kind=int64),intent(in)::angularMomentum
      integer(kind=int64),intent(in),optional::angularRepresentation
      real(kind=real64),dimension(3),intent(in)::center
      real(kind=real64),dimension(:)::coefficients,alphas
!
!     Check if the init flag is FALSE.
!
      if(.not.myBasisSet%initialized)  &
        call mqc_error('MQC_gtoBasisSet_addShell: basis set object not yet initialized.')
!
!     Ensure that adding this shell won't take us past the maximum number of
!     allowed shells. For now, in such cases we fail out.
!
      if(myBasisSet%nShells.ge.myBasisSet%mxShells)  &
        call mqc_error('MQC_gtoBasisSet_addShell: cannot add a new shell.')
!
!     Add the next shell to myBasisSet.
!
      myBasisSet%nShells = myBasisSet%nShells+1
      if(Present(angularRepresentation)) then
        call MQC_CGTF_init(myBasisSet%shells(myBasisSet%nShells),  &
          angularMomentum,center,coefficients,alphas,  &
          angularRepresentation)
      else
        call MQC_CGTF_init(myBasisSet%shells(myBasisSet%nShells),  &
          angularMomentum,center,coefficients,alphas)
      endIf
      myBasisSet%nBasis = myBasisSet%nBasis+  &
        myBasisSet%shells(myBasisSet%nShells)%shell2nBasis()
      myBasisSet%nCartesian = myBasisSet%nCartesian+  &
        myBasisSet%shells(myBasisSet%nShells)%shell2nCartesian()
!
      return
      end subroutine MQC_gtoBasisSet_addShell

!
!PROCEDURE MQC_gtoBasisSet_print
      subroutine MQC_gtoBasisSet_print(myBasisSet,iOut)
!
!     This routine prints the data in a MQC_gtoBasisSet object.
!
      implicit none
      class(MQC_gtoBasisSet)::myBasisSet
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
      end subroutine MQC_gtoBasisSet_print

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
      class(MQC_gtoBasisSet),intent(in)::myBasisSet
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
      class(MQC_gtoBasisSet),intent(in)::myBasisSet
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
      class(MQC_gtoBasisSet)::myBasisSet
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
!PROCEDURE MQC_CGTF_nCartesianFunctions
      function MQC_CGTF_nCartesianFunctions(myCGTF)  result(nCartesian)
!
!     This function returns the number of functions in the complete Cartesian
!     working shell associated with object myCGTF.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      integer(kind=int64)::nCartesian
!
      nCartesian = myCGTF%nCartesian
!
      return
      end function MQC_CGTF_nCartesianFunctions

!
!PROCEDURE MQC_CGTF_getAngularRepresentation
      function MQC_CGTF_getAngularRepresentation(myCGTF)  &
        result(angularRepresentation)
!
!     This function returns the angular representation identifier for the
!     requested/exposed functions in object myCGTF.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      integer(kind=int64)::angularRepresentation
!
      angularRepresentation = myCGTF%angularRepresentation
!
      return
      end function MQC_CGTF_getAngularRepresentation

!
!PROCEDURE MQC_CGTF_getCartesianToBasis
      function MQC_CGTF_getCartesianToBasis(myCGTF)  &
        result(cartesianToBasis)
!
!     This function returns the transformation from the complete Cartesian
!     working shell to the requested/exposed shell representation.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      real(kind=real64),dimension(:,:),allocatable::cartesianToBasis
!
      if(.not.myCGTF%initialized)  &
        call mqc_error('MQC_CGTF_getCartesianToBasis: Object not initialized.')
      cartesianToBasis = myCGTF%cartesianToBasis
!
      return
      end function MQC_CGTF_getCartesianToBasis

!
!PROCEDURE MQC_CGTF_buildCartesianToRealPure
      function MQC_CGTF_buildCartesianToRealPure(myCGTF)  &
        result(cartesianToRealPure)
!
!     This function constructs Gaussian's transformation from normalized
!     Cartesian Gaussian functions to normalized real-pure functions. The
!     Cartesian rows follow the lArrays order stored in myCGTF. Pure columns
!     follow Gaussian C2PGen: 0,1m,1p,2m,2p,... . Numerically, those paired
!     columns are generated by C2PCFF arguments +m and -m, respectively.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      real(kind=real64),dimension(:,:),allocatable::cartesianToRealPure
      real(kind=real64),dimension(:),allocatable::factorial
      integer(kind=int64)::i,l,m,nPure
!
      l = myCGTF%angularMomentum
      nPure = 2*l+1
      Allocate(cartesianToRealPure(myCGTF%nCartesian,nPure))
      Allocate(factorial(0:2*l))
      factorial(0) = mqc_float(1)
      do i = 1,2*l
        factorial(i) = mqc_float(i)*factorial(i-1)
      endDo
!
      do i = 1,myCGTF%nCartesian
        cartesianToRealPure(i,1) =  &
          MQC_CGTF_cartesianToRealPureCoefficient(l,0_int64,  &
            myCGTF%lArrays(1,i),myCGTF%lArrays(2,i),factorial)
        do m = 1,l
          cartesianToRealPure(i,2*m) =  &
            MQC_CGTF_cartesianToRealPureCoefficient(l,m,  &
              myCGTF%lArrays(1,i),myCGTF%lArrays(2,i),factorial)
          cartesianToRealPure(i,2*m+1) =  &
            MQC_CGTF_cartesianToRealPureCoefficient(l,-m,  &
              myCGTF%lArrays(1,i),myCGTF%lArrays(2,i),factorial)
        endDo
      endDo
!
      return
      end function MQC_CGTF_buildCartesianToRealPure

!
!PROCEDURE MQC_CGTF_cartesianToRealPureCoefficient
      function MQC_CGTF_cartesianToRealPureCoefficient(l,m,lx,ly,  &
        factorial) result(coefficient)
!
!     This function evaluates one normalized Cartesian-to-real-pure
!     coefficient using the closed-form expression in Gaussian C2PCFF.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::l,m,lx,ly
      real(kind=real64),dimension(0:),intent(in)::factorial
      real(kind=real64)::coefficient,a,b,c,d,e,g
      integer(kind=int64)::i,iLimit,j,k,lz,mAbsolute
!
      coefficient = mqc_float(0)
      lz = l-lx-ly
      mAbsolute = Abs(m)
      if(mAbsolute.gt.l.or.lz.lt.0) return
!
      a = Sqrt(factorial(2*lx)*factorial(2*ly)*factorial(2*lz)*  &
        factorial(l)/(factorial(2*l)*factorial(lx)*factorial(ly)*  &
        factorial(lz)))
      b = Sqrt(factorial(l-mAbsolute)/factorial(l+mAbsolute))/  &
        ((mqc_float(2)**l)*factorial(l))
      iLimit = (l-mAbsolute)/2
      j = (l-mAbsolute-lz)/2
      g = mqc_float(0)
!
      if(2*j.eq.l-mAbsolute-lz) then
        do i = 0,iLimit
          c = mqc_float(0)
          if(j.ge.0.and.j.le.i) then
            c = (factorial(l)/(factorial(i)*factorial(l-i)))*  &
              (factorial(2*l-2*i)*MQC_CGTF_paritySign(i)/  &
              factorial(l-mAbsolute-2*i))*  &
              (factorial(i)/(factorial(j)*factorial(i-j)))
            do k = 0,j
              d = mqc_float(0)
              if(lx-2*k.ge.0.and.lx-2*k.le.mAbsolute) then
                d = (factorial(j)/(factorial(k)*factorial(j-k)))*  &
                  (factorial(mAbsolute)/(factorial(lx-2*k)*  &
                  factorial(mAbsolute+2*k-lx)))
                e = mqc_float(0)
                if(m.eq.0.and.Mod(lx,2_int64).eq.0)  &
                  e = MQC_CGTF_paritySign(k-lx/2)
                if(m.gt.0.and.Mod(Abs(mAbsolute-lx),2_int64).eq.0)  &
                  e = Sqrt(mqc_float(2))*  &
                    MQC_CGTF_paritySign((2*k+mAbsolute-lx)/2)
                if(m.lt.0.and.Mod(Abs(mAbsolute-lx),2_int64).eq.1)  &
                  e = Sqrt(mqc_float(2))*  &
                    MQC_CGTF_paritySign((2*k+mAbsolute-lx)/2)
                g = g+c*d*e
              endIf
            endDo
          endIf
        endDo
      endIf
      coefficient = a*b*g
!
      return
      end function MQC_CGTF_cartesianToRealPureCoefficient

!
!PROCEDURE MQC_CGTF_paritySign
      function MQC_CGTF_paritySign(power) result(signValue)
!
!     This function evaluates (-1)**power for positive or negative integer
!     powers without relying on compiler handling of negative integer powers.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::power
      real(kind=real64)::signValue
!
      signValue = mqc_float(1)
      if(Mod(Abs(power),2_int64).eq.1) signValue = -mqc_float(1)
!
      return
      end function MQC_CGTF_paritySign

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
!     - H. P. Hratchian, 2025, 2026.
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
      nBasisFunctions = myCGTF%shell2nCartesian()
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
!
!         Gaussian uses its PRISM order for G and higher Cartesian shells.
!         This is the same order used by the lexicographic branch above.
!
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
!     - H. P. Hratchian, 2025, 2026.
!
      implicit none
      class(MQC_CGTF),intent(in)::bfShell
      real(kind=real64),dimension(3),intent(in)::xyz
      real(kind=real64),dimension(:),allocatable,intent(out)::basisValues
      integer(kind=int64)::i,j
      real(kind=real64)::radialValue
      real(kind=real64),dimension(:),allocatable::cartesianValues
      real(kind=real64),dimension(:,:),allocatable::primCoeffs,  &
        primCartesianValues
!
!     Calculate nBasisPairs and nPrimPairs. Then allocate some arrays.
!
      if(.not.bfShell%initialized)  &
        call mqc_error('MQC_Value_CGFT: Basis shell not initialized.')
      if(MQC_Integrals_DEBUG) then
        write(iOut,*)' In MQC_Value_CGFT, basis shell:'
        call MQC_CGTF_print(bfShell,iOut)
      endIf
      if(Allocated(basisValues)) DeAllocate(basisValues)
      Allocate(basisValues(bfShell%nBasis),  &
        cartesianValues(bfShell%nCartesian),  &
        primCartesianValues(bfShell%nPrim,bfShell%nCartesian),  &
        primCoeffs(bfShell%nPrim,bfShell%nCartesian))
!
!     Build the complete normalized contracted Cartesian shell at xyz.
!
      do i = 1,bfShell%nCartesian
        do j = 1,bfShell%nPrim
          primCoeffs(j,i) =  &
            bfShell%primitiveNorm(j,i)*  &
            bfShell%coefficient(j)
          radialValue = MQC_Value_Primitive_Radial(xyz,  &
            bfShell%center,bfShell%alpha(j))
          primCartesianValues(j,i) = radialValue*  &
            MQC_Value_Primitive_Angular(xyz,bfShell%center,  &
            bfShell%alpha(j),bfShell%lArrays(:,i))
        endDo
        cartesianValues(i) =  &
          dot_product(primCoeffs(:,i),primCartesianValues(:,i))
      endDo
!
!     Return either the Cartesian vector or the requested real-pure vector.
!     Cartesian shells carry an identity transformation, so the same path
!     preserves their historical values and ordering.
!
      basisValues = MatMul(Transpose(bfShell%cartesianToBasis),  &
        cartesianValues)
!
      return
      end subroutine MQC_Value_CGFT

!
!PROCEDURE MQC_Value_Primitive_CGTF
      subroutine MQC_Value_Primitive_CGTF(bfShell,iPrimitive,xyz,  &
        basisValues)
!
!     This subroutine evaluates one individually normalized primitive from
!     contracted shell bfShell at xyz. The returned vector uses the shell's
!     requested angular representation and excludes both the primitive
!     contraction coefficient and the contracted-shell renormalization.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      class(MQC_CGTF),intent(in)::bfShell
      integer(kind=int64),intent(in)::iPrimitive
      real(kind=real64),dimension(3),intent(in)::xyz
      real(kind=real64),dimension(:),allocatable,intent(out)::basisValues
      integer(kind=int64)::i
      real(kind=real64)::radialValue
      real(kind=real64),dimension(:),allocatable::cartesianValues
!
      if(.not.bfShell%initialized)  &
        call mqc_error(  &
          'MQC_Value_Primitive_CGTF: Basis shell not initialized.')
      if(iPrimitive.lt.1.or.iPrimitive.gt.bfShell%nPrim)  &
        call mqc_error(  &
          'MQC_Value_Primitive_CGTF: Invalid primitive index sent.')
      if(Allocated(basisValues)) DeAllocate(basisValues)
      Allocate(basisValues(bfShell%nBasis),  &
        cartesianValues(bfShell%nCartesian))
!
      radialValue = MQC_Value_Primitive_Radial(xyz,bfShell%center,  &
        bfShell%alpha(iPrimitive))
      do i = 1,bfShell%nCartesian
        cartesianValues(i) = radialValue*  &
          MQC_Value_Primitive_Angular(xyz,bfShell%center,  &
            bfShell%alpha(iPrimitive),bfShell%lArrays(:,i))/  &
          Sqrt(bfShell%primitiveSelfOverlap(iPrimitive,  &
            bfShell%lArrays(:,i)))
      endDo
      basisValues = MatMul(Transpose(bfShell%cartesianToBasis),  &
        cartesianValues)
!
      return
      end subroutine MQC_Value_Primitive_CGTF

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
!PROCEDURE MQC_Overlap_CGFT
      subroutine MQC_Overlap_CGFT(bra,ket,basisIntegrals)
!
!     This subroutine calculates an overlap integral between two contracted
!     Gaussian type functions. These are sent as MQC_CGTF objects. The output is
!     basisIntegrals, which is a matrix of the overlap integrals over the
!     expansion of basis functions that come from the contracted Gaussians in
!     bra and ket. At present, all overlap integrals are solved using the
!     Obara-Saika recursion relation.
!
!     - H. P. Hratchian, 2025, 2026.
!
      implicit none
      class(MQC_CGTF),intent(in)::bra,ket
      real(kind=real64),dimension(:,:),allocatable,intent(out)::basisIntegrals
      integer(kind=int64)::i,j,k,ixyz,nCartesianBra,nCartesianKet,  &
        nCartesianPairs,nPrimPairs
      integer(kind=int64),dimension(:,:),allocatable::indexBasisPairs,  &
        indexPrimPairs
      real(kind=real64)::mu,p
      real(kind=real64),dimension(3)::xAB,xPA,xPB
      real(kind=real64),dimension(:,:),allocatable::cartesianIntegrals,  &
        tmpSij,primCoeffs,primBasisIntsTot
      real(kind=real64),dimension(:,:,:),allocatable::primBasisInts
      logical,dimension(:,:),allocatable::haveSij
!
 1000 format(1x,'Basis Pair: (',i2,' | ',i2,') = ',f8.5)
 1010 format(3x,'iPrim      IntX          IntY          IntZ        IntTot')
 1020 format(4x,i3,4x,f8.5,6x,f8.5,6x,f8.5,6x,f8.5)
 1021 format(4x,i3,4x,i3,4x,f8.5,6x,f8.5,6x,f8.5,6x,f8.5)
 1030 format(/,/)
!
!     The analytical recurrence operates on complete Cartesian shells. Build
!     that working-basis block for every shell pair, then transform it to the
!     requested exposed representations before returning it.
!
!     Calculate the Cartesian basis-pair and primitive-pair counts. Then
!     allocate the Cartesian working arrays.
!
      if(MQC_Integrals_DEBUG) then
        write(iOut,*)' In MQC_Overlap_CGFT, bra and ket:'
        call MQC_CGTF_print(bra,iOut)
        call MQC_CGTF_print(ket,iOut)
      endIf
      nCartesianBra = bra%nCartesian
      nCartesianKet = ket%nCartesian
      nCartesianPairs = nCartesianBra*nCartesianKet
      nPrimPairs  = bra%nPrim*ket%nPrim
      Allocate(cartesianIntegrals(nCartesianBra,nCartesianKet))
      Allocate(indexBasisPairs(2,nCartesianPairs),  &
        indexPrimPairs(2,nPrimPairs))
      Allocate(primBasisInts(nPrimPairs,3,nCartesianPairs),  &
        primBasisIntsTot(nPrimPairs,nCartesianPairs),  &
        primCoeffs(nPrimPairs,nCartesianPairs))
!
!     Fill indexBasisPairs and indexPrimPairs.
!
      k = 0
      do j = 1,nCartesianKet
        do i = 1,nCartesianBra
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
      do i = 1,nCartesianPairs
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
      do i = 1,nCartesianPairs
        cartesianIntegrals(indexBasisPairs(1,i),indexBasisPairs(2,i)) =  &
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
!     Transform the complete Cartesian block into the exposed bra and ket
!     representations. During shell initialization the representation is
!     temporarily Cartesian and the transformation has not yet been built,
!     so retain a direct Cartesian path.
!
      if(bra%angularRepresentation.eq.MQC_CGTF_REAL_PURE) then
        if(ket%angularRepresentation.eq.MQC_CGTF_REAL_PURE) then
          basisIntegrals = MatMul(Transpose(bra%cartesianToBasis),  &
            MatMul(cartesianIntegrals,ket%cartesianToBasis))
        else
          basisIntegrals = MatMul(Transpose(bra%cartesianToBasis),  &
            cartesianIntegrals)
        endIf
      elseIf(ket%angularRepresentation.eq.MQC_CGTF_REAL_PURE) then
        basisIntegrals = MatMul(cartesianIntegrals,  &
          ket%cartesianToBasis)
      else
        call Move_Alloc(cartesianIntegrals,basisIntegrals)
      endIf
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
        Sij(ixyz+1,2) = XPA*Sij(ixyz,2)+  &
          ((ixyz-1)*Sij(ixyz-1,2)+Sij(ixyz,1))/(mqc_float(2)*p)
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
