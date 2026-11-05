      program unitTest10
!
!     This program tests the reusable one-dimensional, spherical Lebedev, and
!     product-quadrature utilities used for PAD matrix-element integration.
!
!     Command line arguments:
!           No command line arguments are required for the normal test suite.
!           A positive argument selects an invalid-input case that must fail.
!
!     Usage:
!           ./unitTest10.exe
!
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env
      use mqc_general
      use mqc_integrals
      use mqc_integrals, only: LEBEDEV_EXACTNESS_DEGREES,  &
        LEBEDEV_POINT_COUNTS,N_LEBEDEV_RULES,lebedevRuleDegree,  &
        lebedevSelectPointCount
      implicit none
      integer(kind=int64),parameter::iOut=6_int64
      real(kind=real64),parameter::tolAngular=1.0e-10_real64
      real(kind=real64),parameter::tolPolynomial=2.0e-10_real64
      integer(kind=int64)::iRule,invalidInputCase
!
!     Format statements.
!
 1000 format(1x,'Starting unitTest10')
 2000 format(1x,'Lebedev order ',i4,3x,'angular weight sum = ',f20.10)
 2100 format(1x,a,3x,'nPoints = ',i8,3x,'volume = ',f20.10)
 2999 format(/,1x,'unitTest10 complete.')
!
!     Start the unit test program.
!
      invalidInputCase = 0_int64
      if(command_argument_count().gt.0)  &
        call mqc_get_command_argument_integer(1,invalidInputCase)
      if(invalidInputCase.gt.0) then
        call runInvalidInputCase(invalidInputCase)
        call mqc_error('unitTest10: invalid-input case was not rejected.')
      endIf
!
      write(iOut,1000)
      call mqc_version_print(iOut)
!
      call runOneDimensionalRuleCases()
      call runSemiInfiniteRadialRuleCase()
      do iRule = 1,N_LEBEDEV_RULES
        call runLebedevAngularCase(LEBEDEV_POINT_COUNTS(iRule),  &
          LEBEDEV_EXACTNESS_DEGREES(iRule))
      endDo
      call runLebedevSelectionCase()
      call runSphericalProductCase('trapezoid',QUAD_RADIAL_TRAPEZOID,  &
        801_int64,26_int64,1.0e-5_real64)
      call runSphericalProductCase('gauss-legendre',  &
        QUAD_RADIAL_GAUSS_LEGENDRE,6_int64,26_int64,1.0e-10_real64)
      call runSuppliedRuleProductCase()
!
      write(iOut,2999)
!
      contains

!PROCEDURE runOneDimensionalRuleCases
      subroutine runOneDimensionalRuleCases()
!
!     Check trapezoid constants and linear functions and verify that an
!     n-point Gauss-Legendre rule integrates monomials through degree 2*n-1.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),parameter::nGauss=5_int64,nTrap=17_int64
      integer(kind=int64)::power
      real(kind=real64)::actual,expected,lowerBound,stepSize,tolerance,  &
        upperBound
      real(kind=real64),dimension(nGauss)::gaussGrid,gaussWeights
      real(kind=real64),dimension(nTrap)::trapGrid,trapWeights
!
      lowerBound = -mqc_float(1)
      upperBound = mqc_float(3)
      stepSize = (upperBound-lowerBound)/mqc_float(nTrap-1)
      call setup_quadrature_trapezoid1d(nTrap,stepSize,lowerBound,  &
        trapGrid,trapWeights)
      call assertNear(SUM(trapWeights),upperBound-lowerBound,  &
        1.0e-13_real64,'trapezoid constant integral')
      expected = (upperBound*upperBound-lowerBound*lowerBound)/  &
        mqc_float(2)
      call assertNear(dot_product(trapWeights,trapGrid),expected,  &
        1.0e-13_real64,'trapezoid linear integral')
      call assertNear(trapWeights(1),stepSize/mqc_float(2),  &
        1.0e-14_real64,'trapezoid lower endpoint weight')
      call assertNear(trapWeights(nTrap),stepSize/mqc_float(2),  &
        1.0e-14_real64,'trapezoid upper endpoint weight')
!
      lowerBound = -mqc_float(7)/mqc_float(10)
      upperBound = mqc_float(13)/mqc_float(10)
      call setup_quadrature_gauss_legendre1d(nGauss,lowerBound,  &
        upperBound,gaussGrid,gaussWeights)
      do power = 0,2*nGauss-1
        actual = SUM(gaussWeights*gaussGrid**power)
        expected = (upperBound**(power+1)-lowerBound**(power+1))/  &
          mqc_float(power+1)
        tolerance = 2.0e-12_real64*max(mqc_float(1),abs(expected))
        call assertNear(actual,expected,tolerance,  &
          'Gauss-Legendre polynomial exactness')
      endDo
      write(iOut,'(1x,A)') 'One-dimensional quadrature checks passed.'
!
      return
      end subroutine runOneDimensionalRuleCases

!PROCEDURE runSemiInfiniteRadialRuleCase
      subroutine runSemiInfiniteRadialRuleCase()
!
!     Check the mapped Gauss-Legendre radial rule against exponential moments
!     on the semi-infinite interval. Exact integral r**n exp(-r) dr is n!.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),parameter::nRadial=96_int64
      integer(kind=int64)::power
      real(kind=real64)::actual,expected
      real(kind=real64),dimension(nRadial)::radialGrid,radialWeights
!
      call setup_quadrature_gauss_legendre_semi_infinite(nRadial,  &
        mqc_float(1),radialGrid,radialWeights)
      expected = mqc_float(1)
      do power = 0,4
        if(power.gt.0) expected = expected*mqc_float(power)
        actual = SUM(radialWeights*radialGrid**power*exp(-radialGrid))
        call assertNear(actual,expected,2.0e-12_real64,  &
          'semi-infinite radial exponential moment')
      endDo
      if(MINVAL(radialGrid).le.mqc_float(0).or.  &
        MINVAL(radialWeights).le.mqc_float(0))  &
        call mqc_error(  &
          'unitTest10: semi-infinite radial points and weights must be positive.')
      write(iOut,'(1x,A)')  &
        'Semi-infinite radial quadrature checks passed.'
!
      return
      end subroutine runSemiInfiniteRadialRuleCase

!PROCEDURE runLebedevAngularCase
      subroutine runLebedevAngularCase(lebedevPoints,exactnessDegree)
!
!     This routine checks one angular Lebedev grid's normalization, spherical
!     moments, and polynomial exactness through its claimed degree.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::lebedevPoints,exactnessDegree
      integer(kind=int64)::i
      real(kind=real64)::weightSum
      real(kind=real64),dimension(3)::firstMoment,secondMoment
      real(kind=real64),dimension(:),allocatable::angularWeights
      real(kind=real64),dimension(:,:),allocatable::angularGrid
!
      call setup_lebedev_angular_grid(lebedevPoints,angularGrid,  &
        angularWeights)
      if(Size(angularWeights).ne.lebedevPoints) then
        write(iOut,*)' Size(angularWeights) = ',Size(angularWeights)
        write(iOut,*)' lebedevPoints        = ',lebedevPoints
        call mqc_error('unitTest10: unexpected Lebedev grid size.')
      endIf
      weightSum = SUM(angularWeights)
      firstMoment = mqc_float(0)
      secondMoment = mqc_float(0)
      do i = 1,Size(angularWeights)
        call assertNear(sqrt(dot_product(angularGrid(:,i),angularGrid(:,i))),  &
          mqc_float(1),tolAngular,'Lebedev point norm')
        firstMoment = firstMoment+angularWeights(i)*angularGrid(:,i)
        secondMoment = secondMoment+angularWeights(i)*  &
          angularGrid(:,i)*angularGrid(:,i)
      endDo
      write(iOut,'(1x,''Lebedev points '',i4,3x,''degree = '',i3,3x,  &
       &''angular weight sum = '',f20.10)') lebedevPoints,exactnessDegree,  &
        weightSum
      call assertNear(weightSum,mqc_float(4)*Pi,tolAngular,  &
        'Lebedev angular weight sum')
      call assertVectorNear(firstMoment,[ mqc_float(0),mqc_float(0),  &
        mqc_float(0) ],tolAngular,'Lebedev angular first moment')
      call assertVectorNear(secondMoment,[ mqc_float(4)*Pi/mqc_float(3),  &
        mqc_float(4)*Pi/mqc_float(3),mqc_float(4)*Pi/mqc_float(3) ],  &
        tolAngular,'Lebedev angular second moment')
      if(lebedevRuleDegree(lebedevPoints).ne.exactnessDegree)  &
        call mqc_error('unitTest10: Lebedev degree lookup mismatch.')
      call checkLebedevPolynomialExactness(angularGrid,angularWeights,  &
        exactnessDegree)
!
      return
      end subroutine runLebedevAngularCase


!PROCEDURE runLebedevSelectionCase
      subroutine runLebedevSelectionCase()
!
!     This routine checks exact, upward, downward, and capped Lebedev point-
!     count selection.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
!
      if(lebedevSelectPointCount(26_int64).ne.26_int64)  &
        call mqc_error('unitTest10: exact Lebedev selection failed.')
      if(lebedevSelectPointCount(27_int64).ne.38_int64)  &
        call mqc_error('unitTest10: upward Lebedev selection failed.')
      if(lebedevSelectPointCount(27_int64,roundUp=.false.).ne.26_int64)  &
        call mqc_error('unitTest10: downward Lebedev selection failed.')
      if(lebedevSelectPointCount(3000_int64,  &
        capAtMaximum=.true.).ne.2030_int64)  &
        call mqc_error('unitTest10: capped Lebedev selection failed.')
      if(lebedevSelectPointCount(975_int64).ne.1202_int64) &
        call mqc_error('unitTest10: extended upward selection failed.')
      if(lebedevSelectPointCount(2000_int64,roundUp=.false.).ne.1730_int64) &
        call mqc_error('unitTest10: extended downward selection failed.')
      if(lebedevSelectPointCount(1_int64).ne.6_int64)  &
        call mqc_error('unitTest10: minimum Lebedev selection failed.')
!
      return
      end subroutine runLebedevSelectionCase


!PROCEDURE checkLebedevPolynomialExactness
      subroutine checkLebedevPolynomialExactness(angularGrid,angularWeights,  &
        exactnessDegree)
!
!     This routine integrates every Cartesian monomial through the claimed
!     Lebedev exactness degree and compares it with the analytical sphere
!     integral.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::exactnessDegree
      real(kind=real64),dimension(:,:),intent(in)::angularGrid
      real(kind=real64),dimension(:),intent(in)::angularWeights
      integer(kind=int64)::iX,iY,iZ
      real(kind=real64)::actual,error,expected,tolerance
!
      do iX = 0,exactnessDegree
        do iY = 0,exactnessDegree-iX
          do iZ = 0,exactnessDegree-iX-iY
            actual = SUM(angularWeights*angularGrid(1,:)**iX*  &
              angularGrid(2,:)**iY*angularGrid(3,:)**iZ)
            expected = sphericalMonomialIntegral(iX,iY,iZ)
            tolerance = tolPolynomial*MAX(mqc_float(1),abs(expected))
            error = abs(actual-expected)
            if(error.gt.tolerance) then
              write(iOut,*)' monomial powers = ',iX,iY,iZ
              write(iOut,*)' exactnessDegree = ',exactnessDegree
              write(iOut,*)' actual          = ',actual
              write(iOut,*)' expected        = ',expected
              write(iOut,*)' error           = ',error
              write(iOut,*)' tolerance       = ',tolerance
              call mqc_error('unitTest10: Lebedev polynomial exactness failed.')
            endIf
          endDo
        endDo
      endDo
!
      return
      end subroutine checkLebedevPolynomialExactness


!PROCEDURE sphericalMonomialIntegral
      function sphericalMonomialIntegral(iX,iY,iZ) result(integralValue)
!
!     This function returns the analytical unit-sphere integral of
!     x**iX*y**iY*z**iZ.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::iX,iY,iZ
      integer(kind=int64)::i,nHalf
      real(kind=real64)::integralValue
!
      if(mod(iX,2_int64).ne.0_int64.or.mod(iY,2_int64).ne.0_int64.or.  &
        mod(iZ,2_int64).ne.0_int64) then
        integralValue = mqc_float(0)
        return
      endIf
      integralValue = mqc_float(4)*Pi
      do i = 1,iX/2_int64
        integralValue = integralValue*mqc_float(2_int64*i-1_int64)
      endDo
      do i = 1,iY/2_int64
        integralValue = integralValue*mqc_float(2_int64*i-1_int64)
      endDo
      do i = 1,iZ/2_int64
        integralValue = integralValue*mqc_float(2_int64*i-1_int64)
      endDo
      nHalf = (iX+iY+iZ)/2_int64
      do i = 1,nHalf+1_int64
        integralValue = integralValue/mqc_float(2_int64*i-1_int64)
      endDo
!
      return
      end function sphericalMonomialIntegral


!PROCEDURE runSphericalProductCase
      subroutine runSphericalProductCase(label,radialRule,nRadial,  &
        lebedevOrder,tolerance)
!
!     This routine checks one-center spherical product quadrature volume and
!     central Cartesian moments.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::radialRule,nRadial,lebedevOrder
      real(kind=real64),intent(in)::tolerance
      real(kind=real64)::rMaxBohr
      real(kind=real64),dimension(3)::originBohr
      real(kind=real64),dimension(:),allocatable::quadWeights
      real(kind=real64),dimension(:,:),allocatable::quadGrid
!
      rMaxBohr = mqc_float(8)
      originBohr = [ mqc_float(2)/mqc_float(10),  &
        -mqc_float(3)/mqc_float(10),mqc_float(4)/mqc_float(10) ]
      call setup_quadrature_lebedev_spherical(originBohr,rMaxBohr,nRadial,  &
        lebedevOrder,radialRule,quadGrid,quadWeights)
      if(Size(quadWeights).ne.nRadial*lebedevOrder) then
        write(iOut,*)' Size(quadWeights) = ',Size(quadWeights)
        write(iOut,*)' expected          = ',nRadial*lebedevOrder
        call mqc_error('unitTest10: unexpected spherical grid size.')
      endIf
      call checkSphericalProductMoments(label,originBohr,rMaxBohr,quadGrid,  &
        quadWeights,tolerance)
!
      return
      end subroutine runSphericalProductCase


!PROCEDURE runSuppliedRuleProductCase
      subroutine runSuppliedRuleProductCase()
!
!     Build radial and angular rules independently and send them to the
!     low-level spherical-product constructor. This is the composition path
!     needed by future atom-centered molecular grids.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),parameter::lebedevPoints=26_int64,  &
        nRadial=6_int64
      real(kind=real64)::rMaxBohr
      real(kind=real64),dimension(3)::originBohr
      real(kind=real64),dimension(nRadial)::radialGrid,radialWeights
      real(kind=real64),dimension(:),allocatable::angularWeights,quadWeights
      real(kind=real64),dimension(:,:),allocatable::angularGrid,quadGrid
!
      rMaxBohr = mqc_float(8)
      originBohr = [ -mqc_float(1)/mqc_float(5),  &
        mqc_float(1)/mqc_float(10),mqc_float(3)/mqc_float(10) ]
      call setup_quadrature_gauss_legendre1d(nRadial,mqc_float(0),  &
        rMaxBohr,radialGrid,radialWeights)
      call setup_lebedev_angular_grid(lebedevPoints,angularGrid,  &
        angularWeights)
      call setup_quadrature_spherical_product(originBohr,radialGrid,  &
        radialWeights,angularGrid,angularWeights,quadGrid,quadWeights)
      if(Size(quadWeights).ne.nRadial*lebedevPoints)  &
        call mqc_error(  &
          'unitTest10: unexpected supplied-rule spherical grid size.')
      call checkSphericalProductMoments('supplied rules',originBohr,  &
        rMaxBohr,quadGrid,quadWeights,1.0e-10_real64)
!
      return
      end subroutine runSuppliedRuleProductCase


!PROCEDURE checkSphericalProductMoments
      subroutine checkSphericalProductMoments(label,originBohr,rMaxBohr,  &
        quadGrid,quadWeights,tolerance)
!
!     Check the volume and low Cartesian moments of a spherical volume grid.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::label
      real(kind=real64),dimension(3),intent(in)::originBohr
      real(kind=real64),intent(in)::rMaxBohr,tolerance
      real(kind=real64),dimension(:,:),intent(in)::quadGrid
      real(kind=real64),dimension(:),intent(in)::quadWeights
      integer(kind=int64)::i
      real(kind=real64)::expectedSecondMoment,expectedVolume,volume
      real(kind=real64),dimension(3)::centralPoint,firstMoment,secondMoment
!
      if(Size(quadGrid,1).ne.3.or.  &
        Size(quadGrid,2).ne.Size(quadWeights))  &
        call mqc_error('unitTest10: spherical point/weight shapes differ.')
      volume = SUM(quadWeights)
      expectedVolume = mqc_float(4)*Pi*rMaxBohr**3/mqc_float(3)
      expectedSecondMoment = mqc_float(4)*Pi*rMaxBohr**5/mqc_float(15)
      firstMoment = mqc_float(0)
      secondMoment = mqc_float(0)
      do i = 1,Size(quadWeights)
        centralPoint = quadGrid(:,i)-originBohr
        firstMoment = firstMoment+quadWeights(i)*centralPoint
        secondMoment = secondMoment+quadWeights(i)*centralPoint*centralPoint
      endDo
      write(iOut,'(1x,a,3x,''nPoints = '',i8,3x,''volume = '',f20.10)')  &
        label,Size(quadWeights),volume
      call assertNear(volume,expectedVolume,tolerance*expectedVolume,  &
        TRIM(label)//' spherical volume')
      call assertVectorNear(firstMoment,[ mqc_float(0),mqc_float(0),  &
        mqc_float(0) ],tolerance*expectedVolume,  &
        TRIM(label)//' central first moment')
      call assertVectorNear(secondMoment,  &
        [ expectedSecondMoment,expectedSecondMoment,  &
          expectedSecondMoment ],tolerance*expectedSecondMoment,  &
        TRIM(label)//' central second moment')
!
      return
      end subroutine checkSphericalProductMoments




!PROCEDURE runInvalidInputCase
      subroutine runInvalidInputCase(invalidInputCase)
!
!     Exercise one invalid quadrature input. Each case is invoked by the
!     makefile in a separate process and must terminate through mqc_error.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::invalidInputCase
      real(kind=real64),dimension(1)::angularWeights,grid1,weights1
      real(kind=real64),dimension(2)::radialGrid
      real(kind=real64),dimension(1)::radialWeights
      real(kind=real64),dimension(3)::origin
      real(kind=real64),dimension(3,1)::angularGrid
      real(kind=real64),dimension(:),allocatable::quadWeights
      real(kind=real64),dimension(:,:),allocatable::quadGrid
!
      origin = mqc_float(0)
      radialGrid = [ mqc_float(0),mqc_float(1) ]
      radialWeights = [ mqc_float(1) ]
      angularGrid(:,1) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      angularWeights = [ mqc_float(4)*Pi ]
!
      select case(invalidInputCase)
      case(1_int64)
        call setup_quadrature_trapezoid1d(1_int64,mqc_float(1),  &
          mqc_float(0),grid1,weights1)
      case(2_int64)
        call setup_quadrature_spherical_product(origin,radialGrid,  &
          radialWeights,angularGrid,angularWeights,quadGrid,quadWeights)
      case(3_int64)
        angularGrid(:,1) = [ mqc_float(2),mqc_float(0),mqc_float(0) ]
        call setup_quadrature_spherical_product(origin,radialGrid(1:1),  &
          radialWeights,angularGrid,angularWeights,quadGrid,quadWeights)
      case(4_int64)
        call setup_quadrature_gauss_legendre_semi_infinite(1_int64,  &
          mqc_float(0),grid1,weights1)
      case default
        call mqc_error('unitTest10: unknown invalid-input case.')
      end select
!
      return
      end subroutine runInvalidInputCase


!PROCEDURE assertNear
      subroutine assertNear(value,target,tolerance,label)
!
!     This routine checks whether two real values are equal within a specified
!     tolerance.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),intent(in)::value,target,tolerance
      character(len=*),intent(in)::label
!
      if(abs(value-target).gt.tolerance) then
        write(iOut,*)' label     = ',TRIM(label)
        write(iOut,*)' value     = ',value
        write(iOut,*)' target    = ',target
        write(iOut,*)' tolerance = ',tolerance
        call mqc_error('unitTest10: scalar comparison failed.')
      endIf
!
      return
      end subroutine assertNear


!PROCEDURE assertVectorNear
      subroutine assertVectorNear(value,target,tolerance,label)
!
!     This routine checks whether two real vectors are equal within a specified
!     component-wise tolerance.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::value,target
      real(kind=real64),intent(in)::tolerance
      character(len=*),intent(in)::label
!
      if(Size(value).ne.Size(target)) then
        write(iOut,*)' label       = ',TRIM(label)
        write(iOut,*)' Size(value) = ',Size(value)
        write(iOut,*)' Size(target)= ',Size(target)
        call mqc_error('unitTest10: vector shape comparison failed.')
      endIf
      if(maxval(abs(value-target)).gt.tolerance) then
        write(iOut,*)' label     = ',TRIM(label)
        call mqc_print(value,iOut,header='value')
        call mqc_print(target,iOut,header='target')
        write(iOut,*)' tolerance = ',tolerance
        call mqc_error('unitTest10: vector comparison failed.')
      endIf
!
      return
      end subroutine assertVectorNear

      end program unitTest10
