      program unitTest08
!
!     Validate Cartesian and Gaussian-convention real-pure contracted shells,
!     including D/F/G transformations, primitive and contracted point values,
!     Gaussian PRISM Cartesian ordering, and mixed-shell basis assembly.
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env, only: int64,real64
      use MQC_Integrals
      implicit none
!
      type(MQC_CGTF)::cartDefault,cartExplicit,cartF,pureD,pureG
      type(MQC_gtoBasisSet)::mixedBasisSet
      integer(kind=int64),dimension(3,15)::expectedGOrder
      real(kind=real64),dimension(1)::coefficients,alphas
      real(kind=real64),dimension(3)::center,xyz
      real(kind=real64),dimension(:),allocatable::expectedMixedValues,  &
        valuesD,valuesDefault,valuesExplicit,valuesF,valuesG,valuesMixed
      real(kind=real64),dimension(:,:),allocatable::transformation
!
      center = [0.0_real64,0.0_real64,0.0_real64]
      xyz = [0.31_real64,-0.47_real64,0.83_real64]
      coefficients = [1.0_real64]
      alphas = [0.75_real64]
!
!     Omission of angularRepresentation must retain Cartesian behavior.
!
      call cartDefault%init(2_int64,center,coefficients,alphas)
      call assertInteger(cartDefault%shell2nBasis(),6_int64,  &
        'default Cartesian nBasis')
      call assertInteger(cartDefault%shell2nCartesian(),6_int64,  &
        'default Cartesian nCartesian')
      call assertInteger(cartDefault%getAngularRepresentation(),  &
        MQC_CGTF_CARTESIAN,'default Cartesian representation')
      transformation = cartDefault%getCartesianToBasis()
      call assertIdentity(transformation,'default Cartesian transformation')
!
      call cartExplicit%init(2_int64,center,coefficients,alphas,  &
        MQC_CGTF_CARTESIAN)
      call MQC_Value_CGFT(cartDefault,xyz,valuesDefault)
      call MQC_Value_CGFT(cartExplicit,xyz,valuesExplicit)
      call assertVectorClose(valuesDefault,valuesExplicit,  &
        'default and explicit Cartesian values')
!
!     Validate dimensionally equivalent S/P cases and real-pure D/F/G
!     dimensions, transformations, and values.
!
      call assertPureShell(0_int64,1_int64,'real-pure S')
      call assertPureShell(1_int64,3_int64,'real-pure P')
      call assertPureShell(2_int64,5_int64,'real-pure D')
      call assertPureShell(3_int64,7_int64,'real-pure F')
      call assertPureShell(4_int64,9_int64,'real-pure G')
      call assertMultiPrimitiveValues()
      call assertPrimitiveIndependentOfContraction()
!
!     Gaussian uses PRISM ordering for G and higher Cartesian working shells.
!
      call pureG%init(4_int64,center,coefficients,alphas,  &
        MQC_CGTF_REAL_PURE)
      expectedGOrder = Reshape([  &
        0_int64,0_int64,4_int64,  &
        0_int64,1_int64,3_int64,  &
        0_int64,2_int64,2_int64,  &
        0_int64,3_int64,1_int64,  &
        0_int64,4_int64,0_int64,  &
        1_int64,0_int64,3_int64,  &
        1_int64,1_int64,2_int64,  &
        1_int64,2_int64,1_int64,  &
        1_int64,3_int64,0_int64,  &
        2_int64,0_int64,2_int64,  &
        2_int64,1_int64,1_int64,  &
        2_int64,2_int64,0_int64,  &
        3_int64,0_int64,1_int64,  &
        3_int64,1_int64,0_int64,  &
        4_int64,0_int64,0_int64],[3,15])
      call assertIntegerMatrix(pureG%lArrays,expectedGOrder,  &
        'Gaussian G Cartesian order')
!
!     A mixed basis must expose pure dimensions while retaining complete
!     Cartesian working dimensions and concatenate values shell by shell.
!
      call mixedBasisSet%init(3_int64)
      call mixedBasisSet%addShell(2_int64,center,coefficients,alphas,  &
        MQC_CGTF_REAL_PURE)
      call mixedBasisSet%addShell(3_int64,center,coefficients,alphas,  &
        MQC_CGTF_CARTESIAN)
      call mixedBasisSet%addShell(4_int64,center,coefficients,alphas,  &
        MQC_CGTF_REAL_PURE)
      call assertInteger(mixedBasisSet%nBasis,24_int64,  &
        'D-pure/F-Cartesian/G-pure basis-set nBasis')
      call assertInteger(mixedBasisSet%nCartesian,31_int64,  &
        'D-pure/F-Cartesian/G-pure basis-set nCartesian')
!
      call pureD%init(2_int64,center,coefficients,alphas,  &
        MQC_CGTF_REAL_PURE)
      call cartF%init(3_int64,center,coefficients,alphas,  &
        MQC_CGTF_CARTESIAN)
      call MQC_Value_CGFT(pureD,xyz,valuesD)
      call MQC_Value_CGFT(cartF,xyz,valuesF)
      call MQC_Value_CGFT(pureG,xyz,valuesG)
      valuesMixed = basisSetValuesList(mixedBasisSet,xyz)
      expectedMixedValues = [valuesD,valuesF,valuesG]
      call assertVectorClose(valuesMixed,expectedMixedValues,  &
        'mixed basis-set point values')
!
      write(*,'(1x,A)') 'unitTest08: PASS'
!
      contains


!PROCEDURE assertPureShell
      subroutine assertPureShell(angularMomentum,nPure,label)
!
!     Validate one pure shell against an independently constructed normalized
!     Cartesian angular metric and the corresponding Cartesian shell values.
!
      implicit none
      integer(kind=int64),intent(in)::angularMomentum,nPure
      character(len=*),intent(in)::label
      type(MQC_CGTF)::cartesianShell,pureShell
      real(kind=real64),dimension(:,:),allocatable::cartesianMetric,  &
        pureOverlap,shellTransformation
      real(kind=real64),dimension(:),allocatable::cartesianValues,  &
        expectedValues,primitiveCartesianValues,primitivePureValues,  &
        pureValues
!
      call cartesianShell%init(angularMomentum,center,coefficients,alphas,  &
        MQC_CGTF_CARTESIAN)
      call pureShell%init(angularMomentum,center,coefficients,alphas,  &
        MQC_CGTF_REAL_PURE)
      call assertInteger(pureShell%shell2nBasis(),nPure,  &
        Trim(label)//' nBasis')
      call assertInteger(pureShell%shell2nCartesian(),  &
        ((angularMomentum+1)*(angularMomentum+2))/2,  &
        Trim(label)//' nCartesian')
      call assertInteger(pureShell%getAngularRepresentation(),  &
        MQC_CGTF_REAL_PURE,Trim(label)//' representation')
!
      shellTransformation = pureShell%getCartesianToBasis()
      if(Size(shellTransformation,1).ne.  &
        cartesianShell%shell2nCartesian().or.  &
        Size(shellTransformation,2).ne.nPure)  &
        call fail(Trim(label)//' transformation shape differs.')
!
      cartesianMetric = buildNormalizedCartesianMetric(  &
        cartesianShell%lArrays)
      pureOverlap = MatMul(Transpose(shellTransformation),  &
        MatMul(cartesianMetric,shellTransformation))
      call assertIdentity(pureOverlap,Trim(label)//' transformed overlap')
!
      call cartesianShell%primitiveValues(1_int64,xyz,  &
        primitiveCartesianValues)
      call pureShell%primitiveValues(1_int64,xyz,primitivePureValues)
      expectedValues = MatMul(Transpose(shellTransformation),  &
        primitiveCartesianValues)
      call assertVectorClose(primitivePureValues,expectedValues,  &
        Trim(label)//' primitive values')
!
      call MQC_Value_CGFT(cartesianShell,xyz,cartesianValues)
      call MQC_Value_CGFT(pureShell,xyz,pureValues)
      expectedValues = MatMul(Transpose(shellTransformation),  &
        cartesianValues)
      call assertVectorClose(pureValues,expectedValues,  &
        Trim(label)//' contracted values')
!
      return
      end subroutine assertPureShell


!PROCEDURE assertMultiPrimitiveValues
      subroutine assertMultiPrimitiveValues()
!
!     Validate transformation after a nontrivial three-primitive contraction.
!
      implicit none
      type(MQC_CGTF)::cartesianShell,pureShell
      real(kind=real64),dimension(3)::multiCoefficients,multiAlphas
      real(kind=real64),dimension(:),allocatable::cartesianValues,  &
        expectedValues,pureValues
      real(kind=real64),dimension(:,:),allocatable::shellTransformation
!
      multiCoefficients = [0.31_real64,-0.27_real64,0.83_real64]
      multiAlphas = [3.20_real64,0.91_real64,0.28_real64]
      call cartesianShell%init(3_int64,center,multiCoefficients,  &
        multiAlphas,MQC_CGTF_CARTESIAN)
      call pureShell%init(3_int64,center,multiCoefficients,multiAlphas,  &
        MQC_CGTF_REAL_PURE)
      shellTransformation = pureShell%getCartesianToBasis()
      call MQC_Value_CGFT(cartesianShell,xyz,cartesianValues)
      call MQC_Value_CGFT(pureShell,xyz,pureValues)
      expectedValues = MatMul(Transpose(shellTransformation),  &
        cartesianValues)
      call assertVectorClose(pureValues,expectedValues,  &
        'multi-primitive real-pure F contracted values')
!
      return
      end subroutine assertMultiPrimitiveValues


!PROCEDURE assertPrimitiveIndependentOfContraction
      subroutine assertPrimitiveIndependentOfContraction()
!
!     The value of one individually normalized primitive must not depend on
!     other primitives or coefficients in the contracted shell containing it.
!
      implicit none
      type(MQC_CGTF)::singleShell,multiShell
      real(kind=real64),dimension(1)::singleAlpha,singleCoefficient
      real(kind=real64),dimension(3)::multiAlphas,multiCoefficients
      real(kind=real64),dimension(:),allocatable::multiValues,singleValues
!
      singleAlpha = [0.91_real64]
      singleCoefficient = [1.0_real64]
      multiAlphas = [3.20_real64,0.91_real64,0.28_real64]
      multiCoefficients = [0.31_real64,-0.27_real64,0.83_real64]
      call singleShell%init(3_int64,center,singleCoefficient,singleAlpha,  &
        MQC_CGTF_REAL_PURE)
      call multiShell%init(3_int64,center,multiCoefficients,multiAlphas,  &
        MQC_CGTF_REAL_PURE)
      call singleShell%primitiveValues(1_int64,xyz,singleValues)
      call multiShell%primitiveValues(2_int64,xyz,multiValues)
      call assertVectorClose(multiValues,singleValues,  &
        'primitive value independent of contraction')
!
      return
      end subroutine assertPrimitiveIndependentOfContraction


!PROCEDURE buildNormalizedCartesianMetric
      function buildNormalizedCartesianMetric(lArrays) result(metric)
!
!     Build the exponent-independent angular metric between individually
!     normalized Cartesian components of one same-center shell.
!
      implicit none
      integer(kind=int64),dimension(:,:),intent(in)::lArrays
      real(kind=real64),dimension(:,:),allocatable::metric
      integer(kind=int64)::i,j,k,lTotal
      real(kind=real64)::normalizationI,normalizationJ
!
      Allocate(metric(Size(lArrays,2),Size(lArrays,2)))
      metric = 0.0_real64
      do i = 1,Size(lArrays,2)
        normalizationI = 1.0_real64
        do k = 1,3
          normalizationI = normalizationI*  &
            factorialReal(lArrays(k,i))/factorialReal(2*lArrays(k,i))
        endDo
        normalizationI = Sqrt(normalizationI)
        do j = 1,Size(lArrays,2)
          if(Any(Mod(lArrays(:,i)+lArrays(:,j),2_int64).ne.0)) cycle
          normalizationJ = 1.0_real64
          do k = 1,3
            normalizationJ = normalizationJ*  &
              factorialReal(lArrays(k,j))/  &
              factorialReal(2*lArrays(k,j))
          endDo
          normalizationJ = Sqrt(normalizationJ)
          metric(i,j) = normalizationI*normalizationJ
          do k = 1,3
            lTotal = lArrays(k,i)+lArrays(k,j)
            metric(i,j) = metric(i,j)*  &
              factorialReal(lTotal)/factorialReal(lTotal/2)
          endDo
        endDo
      endDo
!
      return
      end function buildNormalizedCartesianMetric


!PROCEDURE factorialReal
      function factorialReal(n) result(value)
      implicit none
      integer(kind=int64),intent(in)::n
      real(kind=real64)::value
      integer(kind=int64)::i
!
      value = 1.0_real64
      do i = 2,n
        value = value*Real(i,kind=real64)
      endDo
!
      return
      end function factorialReal


!PROCEDURE assertInteger
      subroutine assertInteger(actual,expected,label)
      implicit none
      integer(kind=int64),intent(in)::actual,expected
      character(len=*),intent(in)::label
!
      if(actual.ne.expected) call fail(Trim(label)//' differs.')
!
      return
      end subroutine assertInteger


!PROCEDURE assertIntegerMatrix
      subroutine assertIntegerMatrix(actual,expected,label)
      implicit none
      integer(kind=int64),dimension(:,:),intent(in)::actual,expected
      character(len=*),intent(in)::label
!
      if(Any(Shape(actual).ne.Shape(expected)))  &
        call fail(Trim(label)//' shape differs.')
      if(Any(actual.ne.expected)) call fail(Trim(label)//' differs.')
!
      return
      end subroutine assertIntegerMatrix


!PROCEDURE assertIdentity
      subroutine assertIdentity(matrix,label)
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::matrix
      character(len=*),intent(in)::label
      integer(kind=int64)::i,j
      real(kind=real64)::expected
!
      if(Size(matrix,1).ne.Size(matrix,2))  &
        call fail(Trim(label)//' is not square.')
      do i = 1,Size(matrix,1)
        do j = 1,Size(matrix,2)
          expected = 0.0_real64
          if(i.eq.j) expected = 1.0_real64
          if(Abs(matrix(i,j)-expected).gt.1.0e-11_real64)  &
            call fail(Trim(label)//' differs.')
        endDo
      endDo
!
      return
      end subroutine assertIdentity


!PROCEDURE assertVectorClose
      subroutine assertVectorClose(actual,expected,label)
      implicit none
      real(kind=real64),dimension(:),intent(in)::actual,expected
      character(len=*),intent(in)::label
!
      if(Size(actual).ne.Size(expected))  &
        call fail(Trim(label)//' size differs.')
      if(Any(Abs(actual-expected).gt.1.0e-11_real64))  &
        call fail(Trim(label)//' differs.')
!
      return
      end subroutine assertVectorClose


!PROCEDURE fail
      subroutine fail(message)
      implicit none
      character(len=*),intent(in)::message
!
      write(*,'(1x,A)') 'FAIL: '//Trim(message)
      error stop 1
!
      end subroutine fail


      end program unitTest08
