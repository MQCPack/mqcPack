      program unitTest11
!
!     This program validates the original-Becke molecular partition and the
!     unpruned atom-centered molecular-grid assembly.
!
!     Usage:
!           ./unitTest11.exe
!
!     A positive command-line argument selects an invalid-input case that
!     must terminate through mqc_error.
!
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env, only: int64,real64
      use mqc_general, only: angPBohr,mqc_distance,mqc_error,  &
        mqc_get_command_argument_integer,mqc_version_print
      use mqc_integrals, only: basisSetOverlapMatrix,  &
        basisSetValuesList1,MQC_gtoBasisSet
      use mqc_integrals, only: apply_becke1988_partition,  &
        becke1988_partition_weights,setup_quadrature_becke_molecular
      use mqc_integrals, only:  &
        becke1988_radii_from_atomic_numbers,  &
        setup_quadrature_becke_molecular_from_atomic_numbers
      use mqc_integrals, only: setup_lebedev_angular_grid,  &
        setup_quadrature_gauss_legendre1d
      implicit none
      integer(kind=int64),parameter::iOut=6_int64
      integer(kind=int64)::invalidInputCase
!
      invalidInputCase = 0_int64
      if(command_argument_count().gt.0)  &
        call mqc_get_command_argument_integer(1,invalidInputCase)
      if(invalidInputCase.gt.0) then
        call runInvalidInputCase(invalidInputCase)
        call mqc_error('unitTest11: invalid-input case was not rejected.')
      endIf
!
      write(iOut,'(1x,A)') 'Starting unitTest11'
      call mqc_version_print(iOut)
      call runSingleAtomCase()
      call runHomonuclearCase()
      call runHeteronuclearCase()
      call runInvarianceCases()
      call runMolecularGridAssemblyCase()
      call runGaussianOverlapCase()
      write(iOut,'(/,1x,A)') 'unitTest11 complete.'
!
      contains

!PROCEDURE runSingleAtomCase
      subroutine runSingleAtomCase()
!
!     Confirm that the one-atom molecular partition is exactly one at every
!     point, including the atomic center.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(3,1)::centerCoordsBohr
      real(kind=real64),dimension(1)::atomicRadiiBohr
      real(kind=real64),dimension(3,4)::pointsBohr
      real(kind=real64),dimension(:,:),allocatable::partitionWeights
!
      centerCoordsBohr(:,1) =  &
        [ 0.3_real64,-0.4_real64,0.2_real64 ]
      atomicRadiiBohr = [ 0.7_real64 ]
      pointsBohr(:,1) = centerCoordsBohr(:,1)
      pointsBohr(:,2) = [ 1.1_real64,0.2_real64,-0.7_real64 ]
      pointsBohr(:,3) = [ -2.0_real64,3.0_real64,0.5_real64 ]
      pointsBohr(:,4) = [ 20.0_real64,-15.0_real64,8.0_real64 ]
      call becke1988_partition_weights(centerCoordsBohr,atomicRadiiBohr,  &
        pointsBohr,partitionWeights)
      call assertMatrixNear(partitionWeights,  &
        reshape([ 1.0_real64,1.0_real64,1.0_real64,1.0_real64 ],  &
          [ 1,4 ]),0.0_real64,'one-atom partition')
      write(iOut,'(1x,A)') 'One-atom partition checks passed.'
!
      return
      end subroutine runSingleAtomCase

!PROCEDURE runHomonuclearCase
      subroutine runHomonuclearCase()
!
!     Check partition of unity, atomic-center limits, and equal bisector
!     weights for a homonuclear diatomic molecule.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),parameter::tolerance=2.0e-14_real64
      real(kind=real64),dimension(3,2)::centerCoordsBohr
      real(kind=real64),dimension(2)::atomicRadiiBohr
      real(kind=real64),dimension(3,5)::pointsBohr
      real(kind=real64),dimension(:,:),allocatable::partitionWeights
!
      centerCoordsBohr(:,1) = [ -1.0_real64,0.0_real64,0.0_real64 ]
      centerCoordsBohr(:,2) = [ 1.0_real64,0.0_real64,0.0_real64 ]
      atomicRadiiBohr = [ 0.7_real64,0.7_real64 ]
      pointsBohr(:,1) = centerCoordsBohr(:,1)
      pointsBohr(:,2) = [ 0.0_real64,0.0_real64,0.0_real64 ]
      pointsBohr(:,3) = centerCoordsBohr(:,2)
      pointsBohr(:,4) = [ 0.0_real64,0.4_real64,-0.3_real64 ]
      pointsBohr(:,5) = [ -0.35_real64,0.2_real64,0.1_real64 ]
      call becke1988_partition_weights(centerCoordsBohr,atomicRadiiBohr,  &
        pointsBohr,partitionWeights)
      call assertVectorNear(SUM(partitionWeights,dim=1),  &
        [ 1.0_real64,1.0_real64,1.0_real64,1.0_real64,1.0_real64 ],  &
        tolerance,'homonuclear partition of unity')
      call assertNear(partitionWeights(1,1),1.0_real64,tolerance,  &
        'first atomic-center weight')
      call assertNear(partitionWeights(2,1),0.0_real64,tolerance,  &
        'first atomic-center complement')
      call assertNear(partitionWeights(1,3),0.0_real64,tolerance,  &
        'second atomic-center complement')
      call assertNear(partitionWeights(2,3),1.0_real64,tolerance,  &
        'second atomic-center weight')
      call assertNear(partitionWeights(1,2),0.5_real64,tolerance,  &
        'bond midpoint weight')
      call assertNear(partitionWeights(1,4),0.5_real64,tolerance,  &
        'off-axis bisector weight')
      call assertPartitionBounds(partitionWeights,tolerance,  &
        'homonuclear partition bounds')
      write(iOut,'(1x,A)') 'Homonuclear partition checks passed.'
!
      return
      end subroutine runHomonuclearCase

!PROCEDURE runHeteronuclearCase
      subroutine runHeteronuclearCase()
!
!     Validate MQCPack radius lookup and a fixed original-Becke H--O
!     heteronuclear midpoint reference value.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),dimension(2)::atomicNumbers
      real(kind=real64),dimension(3,2)::centerCoordsBohr
      real(kind=real64),dimension(3,1)::pointsBohr
      real(kind=real64),dimension(:),allocatable::atomicRadiiBohr
      real(kind=real64),dimension(:,:),allocatable::partitionWeights
!
      atomicNumbers = [ 1_int64,8_int64 ]
      call becke1988_radii_from_atomic_numbers(atomicNumbers,  &
        atomicRadiiBohr)
      call assertNear(atomicRadiiBohr(1)*angPBohr,0.35_real64,  &
        1.0e-12_real64,'Becke hydrogen radius')
      call assertNear(atomicRadiiBohr(2)*angPBohr,0.60_real64,  &
        1.0e-12_real64,'Becke oxygen radius')
      centerCoordsBohr(:,1) = [ -1.0_real64,0.0_real64,0.0_real64 ]
      centerCoordsBohr(:,2) = [ 1.0_real64,0.0_real64,0.0_real64 ]
      pointsBohr(:,1) = [ 0.0_real64,0.0_real64,0.0_real64 ]
      call becke1988_partition_weights(centerCoordsBohr,atomicRadiiBohr,  &
        pointsBohr,partitionWeights)
      call assertNear(partitionWeights(1,1),  &
        0.11177633883081156_real64,2.0e-14_real64,  &
        'H--O midpoint hydrogen weight')
      call assertNear(partitionWeights(2,1),  &
        0.8882236611691885_real64,2.0e-14_real64,  &
        'H--O midpoint oxygen weight')
      write(iOut,'(1x,A)') 'Heteronuclear partition checks passed.'
!
      return
      end subroutine runHeteronuclearCase

!PROCEDURE runInvarianceCases
      subroutine runInvarianceCases()
!
!     Check atom-list permutation invariance and common molecular translation
!     and rotation invariance for a three-center partition.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),dimension(3)::permutation
      real(kind=real64),parameter::tolerance=5.0e-14_real64
      real(kind=real64),dimension(3,3)::centerCoordsBohr,  &
        permutedCentersBohr,rotatedCentersBohr,rotation,  &
        translatedCentersBohr
      real(kind=real64),dimension(3)::atomicRadiiBohr,  &
        permutedRadiiBohr,translationBohr
      real(kind=real64),dimension(3,5)::pointsBohr,rotatedPointsBohr,  &
        translatedPointsBohr
      real(kind=real64),dimension(:,:),allocatable::partitionWeights,  &
        permutedWeights,rotatedWeights,translatedWeights
      integer(kind=int64)::i
!
      centerCoordsBohr(:,1) = [ -0.8_real64,0.1_real64,0.0_real64 ]
      centerCoordsBohr(:,2) = [ 0.9_real64,-0.2_real64,0.3_real64 ]
      centerCoordsBohr(:,3) = [ 0.2_real64,1.1_real64,-0.4_real64 ]
      atomicRadiiBohr = [ 0.7_real64,0.6_real64,1.2_real64 ]
      pointsBohr(:,1) = [ 0.1_real64,0.2_real64,0.3_real64 ]
      pointsBohr(:,2) = [ -1.2_real64,0.5_real64,-0.1_real64 ]
      pointsBohr(:,3) = [ 0.7_real64,1.4_real64,0.8_real64 ]
      pointsBohr(:,4) = [ 2.1_real64,-1.3_real64,0.4_real64 ]
      pointsBohr(:,5) = [ -0.3_real64,0.9_real64,-1.1_real64 ]
      call becke1988_partition_weights(centerCoordsBohr,atomicRadiiBohr,  &
        pointsBohr,partitionWeights)
!     Independent direct-product evaluation of the published equations.
      call assertVectorNear(partitionWeights(:,1),  &
        [ 0.18211265927496753_real64,  &
          0.093195460072559147_real64,  &
          0.72469188065247325_real64 ],tolerance,  &
        'three-center reference weights')
!
      permutation = [ 3_int64,1_int64,2_int64 ]
      permutedCentersBohr = centerCoordsBohr(:,permutation)
      permutedRadiiBohr = atomicRadiiBohr(permutation)
      call becke1988_partition_weights(permutedCentersBohr,  &
        permutedRadiiBohr,pointsBohr,permutedWeights)
      do i = 1,3
        call assertVectorNear(permutedWeights(i,:),  &
          partitionWeights(permutation(i),:),tolerance,  &
          'atom-list permutation invariance')
      endDo
!
      translationBohr = [ 1.7_real64,-0.6_real64,2.3_real64 ]
      do i = 1,3
        translatedCentersBohr(i,:) =  &
          centerCoordsBohr(i,:)+translationBohr(i)
        translatedPointsBohr(i,:) = pointsBohr(i,:)+translationBohr(i)
      endDo
      call becke1988_partition_weights(translatedCentersBohr,  &
        atomicRadiiBohr,translatedPointsBohr,translatedWeights)
      call assertMatrixNear(translatedWeights,partitionWeights,tolerance,  &
        'translation invariance')
!
      rotation = reshape([ 0.0_real64,1.0_real64,0.0_real64,  &
        -1.0_real64,0.0_real64,0.0_real64,  &
        0.0_real64,0.0_real64,1.0_real64 ],[ 3,3 ])
      rotatedCentersBohr = matmul(rotation,centerCoordsBohr)
      rotatedPointsBohr = matmul(rotation,pointsBohr)
      call becke1988_partition_weights(rotatedCentersBohr,  &
        atomicRadiiBohr,rotatedPointsBohr,rotatedWeights)
      call assertMatrixNear(rotatedWeights,partitionWeights,tolerance,  &
        'rotation invariance')
      call assertVectorNear(SUM(partitionWeights,dim=1),  &
        [ 1.0_real64,1.0_real64,1.0_real64,1.0_real64,1.0_real64 ],  &
        tolerance,'three-center partition of unity')
      write(iOut,'(1x,A)') 'Partition invariance checks passed.'
!
      return
      end subroutine runInvarianceCases

!PROCEDURE runMolecularGridAssemblyCase
      subroutine runMolecularGridAssemblyCase()
!
!     Build a two-atom molecular grid through both the explicit-radii and
!     atomic-number interfaces and verify construction metadata and weights.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),parameter::nRadial=4_int64
      integer(kind=int64),dimension(2)::atomicNumbers
      integer(kind=int64),dimension(:),allocatable::generatingAtoms,  &
        generatingAtomsFromNumbers
      integer(kind=int64)::nAtomicPoints,nPoints
      real(kind=real64),dimension(3,2)::centerCoordsBohr
      real(kind=real64),dimension(nRadial)::radialGridBohr,  &
        radialWeightsBohr
      real(kind=real64),dimension(:),allocatable::angularWeights,  &
        atomicRadiiBohr,baseWeightsBohr3,baseWeightsFromNumbers,  &
        finalWeightsFromApply,partitionFromApply,partitionWeights,  &
        partitionWeightsFromNumbers,quadWeightsBohr3,  &
        quadWeightsFromNumbers
      real(kind=real64),dimension(:,:),allocatable::angularGrid,  &
        quadGridBohr,quadGridFromNumbers
!
      atomicNumbers = [ 1_int64,8_int64 ]
      centerCoordsBohr(:,1) = [ -0.8_real64,0.2_real64,-0.1_real64 ]
      centerCoordsBohr(:,2) = [ 1.0_real64,-0.3_real64,0.4_real64 ]
      call setup_quadrature_gauss_legendre1d(nRadial,0.0_real64,  &
        3.0_real64,radialGridBohr,radialWeightsBohr)
      call setup_lebedev_angular_grid(6_int64,angularGrid,angularWeights)
      call becke1988_radii_from_atomic_numbers(atomicNumbers,  &
        atomicRadiiBohr)
      call setup_quadrature_becke_molecular(centerCoordsBohr,  &
        atomicRadiiBohr,radialGridBohr,radialWeightsBohr,angularGrid,  &
        angularWeights,quadGridBohr,generatingAtoms,baseWeightsBohr3,  &
        partitionWeights,quadWeightsBohr3)
      call setup_quadrature_becke_molecular_from_atomic_numbers(  &
        atomicNumbers,centerCoordsBohr,radialGridBohr,radialWeightsBohr,  &
        angularGrid,angularWeights,quadGridFromNumbers,  &
        generatingAtomsFromNumbers,baseWeightsFromNumbers,  &
        partitionWeightsFromNumbers,quadWeightsFromNumbers)
!
      nAtomicPoints = nRadial*Size(angularWeights)
      nPoints = Size(atomicNumbers)*nAtomicPoints
      if(Size(quadWeightsBohr3).ne.nPoints)  &
        call mqc_error('unitTest11: unexpected molecular-grid size.')
      if(ANY(generatingAtoms(1:nAtomicPoints).ne.1_int64).or.  &
        ANY(generatingAtoms(nAtomicPoints+1:nPoints).ne.2_int64))  &
        call mqc_error('unitTest11: generating-atom blocks are incorrect.')
      call assertNear(mqc_distance(quadGridBohr(:,1),  &
        centerCoordsBohr(:,1)),radialGridBohr(1),1.0e-13_real64,  &
        'first atomic-grid origin')
      call assertNear(mqc_distance(quadGridBohr(:,nAtomicPoints+1),  &
        centerCoordsBohr(:,2)),radialGridBohr(1),1.0e-13_real64,  &
        'second atomic-grid origin')
      call assertVectorNear(quadWeightsBohr3,  &
        baseWeightsBohr3*partitionWeights,1.0e-14_real64,  &
        'assembled final weights')
      call assertPartitionBounds(  &
        reshape(partitionWeights,[ 1,nPoints ]),1.0e-14_real64,  &
        'assembled partition bounds')
      call assertMatrixNear(quadGridFromNumbers,quadGridBohr,  &
        1.0e-14_real64,'atomic-number adapter points')
      call assertIntegerVectorEqual(generatingAtomsFromNumbers,  &
        generatingAtoms,'atomic-number adapter generating atoms')
      call assertVectorNear(baseWeightsFromNumbers,baseWeightsBohr3,  &
        1.0e-14_real64,'atomic-number adapter base weights')
      call assertVectorNear(partitionWeightsFromNumbers,partitionWeights,  &
        1.0e-14_real64,'atomic-number adapter partition weights')
      call assertVectorNear(quadWeightsFromNumbers,quadWeightsBohr3,  &
        1.0e-14_real64,'atomic-number adapter final weights')
!
      call apply_becke1988_partition(centerCoordsBohr,atomicRadiiBohr,  &
        quadGridBohr,generatingAtoms,baseWeightsBohr3,partitionFromApply,  &
        finalWeightsFromApply)
      call assertVectorNear(partitionFromApply,partitionWeights,  &
        1.0e-14_real64,'partition application path')
      call assertVectorNear(finalWeightsFromApply,quadWeightsBohr3,  &
        1.0e-14_real64,'final-weight application path')
      write(iOut,'(1x,A)') 'Molecular-grid assembly checks passed.'
!
      return
      end subroutine runMolecularGridAssemblyCase

!PROCEDURE runGaussianOverlapCase
      subroutine runGaussianOverlapCase()
!
!     Compare the complete numerical overlap matrix for two normalized
!     atom-centered s Gaussians against MQCPack analytical integrals. This
!     simultaneously tests diagonal normalization and the off-center overlap.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),parameter::nRadial=80_int64
      type(MQC_gtoBasisSet)::basisSet
      integer(kind=int64),dimension(:),allocatable::generatingAtoms
      integer(kind=int64)::i,j,iPoint
      real(kind=real64),dimension(3,2)::centerCoordsBohr
      real(kind=real64),dimension(2)::atomicRadiiBohr
      real(kind=real64),dimension(1)::coefficients,exponents
      real(kind=real64),dimension(nRadial)::radialGridBohr,  &
        radialWeightsBohr
      real(kind=real64),dimension(:),allocatable::angularWeights,  &
        baseWeightsBohr3,basisValues,partitionWeights,  &
        quadWeightsBohr3
      real(kind=real64),dimension(:,:),allocatable::analyticalOverlap,  &
        angularGrid,numericalOverlap,quadGridBohr
!
      centerCoordsBohr(:,1) = [ 0.0_real64,0.0_real64,-0.7_real64 ]
      centerCoordsBohr(:,2) = [ 0.0_real64,0.0_real64,0.7_real64 ]
      atomicRadiiBohr = [ 0.7_real64,0.7_real64 ]
      coefficients = [ 1.0_real64 ]
      exponents = [ 1.0_real64 ]
      call basisSet%init(2_int64)
      call basisSet%addShell(0_int64,centerCoordsBohr(:,1),coefficients,  &
        exponents)
      call basisSet%addShell(0_int64,centerCoordsBohr(:,2),coefficients,  &
        exponents)
      analyticalOverlap = basisSetOverlapMatrix(basisSet)
!
      call setup_quadrature_gauss_legendre1d(nRadial,0.0_real64,  &
        8.0_real64,radialGridBohr,radialWeightsBohr)
      call setup_lebedev_angular_grid(302_int64,angularGrid,  &
        angularWeights)
      call setup_quadrature_becke_molecular(centerCoordsBohr,  &
        atomicRadiiBohr,radialGridBohr,radialWeightsBohr,angularGrid,  &
        angularWeights,quadGridBohr,generatingAtoms,baseWeightsBohr3,  &
        partitionWeights,quadWeightsBohr3)
      Allocate(numericalOverlap(basisSet%nBasis,basisSet%nBasis))
      numericalOverlap = 0.0_real64
      do iPoint = 1,Size(quadWeightsBohr3)
        call basisSetValuesList1(basisSet,quadGridBohr(:,iPoint),  &
          basisValues)
        do i = 1,basisSet%nBasis
          do j = 1,basisSet%nBasis
            numericalOverlap(i,j) = numericalOverlap(i,j)+  &
              quadWeightsBohr3(iPoint)*basisValues(i)*basisValues(j)
          endDo
        endDo
      endDo
      call assertMatrixNear(numericalOverlap,analyticalOverlap,  &
        2.0e-9_real64,'Gaussian overlap matrix')
      write(iOut,'(1x,A,ES12.4)')  &
        'Gaussian overlap maximum error = ',  &
        MAXVAL(abs(numericalOverlap-analyticalOverlap))
!
      return
      end subroutine runGaussianOverlapCase

!PROCEDURE runInvalidInputCase
      subroutine runInvalidInputCase(invalidInputCase)
!
!     Exercise one invalid partition or molecular-grid input. Each case is run
!     in a separate process by the makefile and must terminate through
!     mqc_error.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::invalidInputCase
      integer(kind=int64),dimension(1)::atomicNumbers,generatingAtoms
      real(kind=real64),dimension(2,1)::badCenters
      real(kind=real64),dimension(3,1)::oneCenter,pointsBohr
      real(kind=real64),dimension(3,2)::twoCenters
      real(kind=real64),dimension(1)::baseWeights,oneRadius
      real(kind=real64),dimension(2)::twoRadii
      real(kind=real64),dimension(:),allocatable::finalWeights,  &
        partitionVector,radii
      real(kind=real64),dimension(:,:),allocatable::partitionWeights
!
      badCenters = 0.0_real64
      oneCenter = 0.0_real64
      pointsBohr = 0.0_real64
      twoCenters = 0.0_real64
      twoCenters(:,2) = [ 1.0_real64,0.0_real64,0.0_real64 ]
      oneRadius = [ 1.0_real64 ]
      twoRadii = [ 1.0_real64,1.0_real64 ]
      generatingAtoms = [ 1_int64 ]
      baseWeights = [ 1.0_real64 ]
      atomicNumbers = [ 1_int64 ]
!
      select case(invalidInputCase)
      case(1_int64)
        call becke1988_partition_weights(badCenters,oneRadius,pointsBohr,  &
          partitionWeights)
      case(2_int64)
        call becke1988_partition_weights(twoCenters,oneRadius,pointsBohr,  &
          partitionWeights)
      case(3_int64)
        oneRadius = [ 0.0_real64 ]
        call becke1988_partition_weights(oneCenter,oneRadius,pointsBohr,  &
          partitionWeights)
      case(4_int64)
        twoCenters(:,2) = twoCenters(:,1)
        call becke1988_partition_weights(twoCenters,twoRadii,pointsBohr,  &
          partitionWeights)
      case(5_int64)
        generatingAtoms = [ 2_int64 ]
        call apply_becke1988_partition(oneCenter,oneRadius,pointsBohr,  &
          generatingAtoms,baseWeights,partitionVector,finalWeights)
      case(6_int64)
        atomicNumbers = [ 0_int64 ]
        call becke1988_radii_from_atomic_numbers(atomicNumbers,radii)
      case default
        call mqc_error('unitTest11: unknown invalid-input case.')
      end select
!
      return
      end subroutine runInvalidInputCase

!PROCEDURE assertPartitionBounds
      subroutine assertPartitionBounds(partitionWeights,tolerance,label)
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::partitionWeights
      real(kind=real64),intent(in)::tolerance
      character(len=*),intent(in)::label
!
      if(MINVAL(partitionWeights).lt.-tolerance.or.  &
        MAXVAL(partitionWeights).gt.1.0_real64+tolerance)  &
        call mqc_error('unitTest11: '//TRIM(label)//' failed.')
!
      return
      end subroutine assertPartitionBounds

!PROCEDURE assertNear
      subroutine assertNear(value,target,tolerance,label)
      implicit none
      real(kind=real64),intent(in)::value,target,tolerance
      character(len=*),intent(in)::label
!
      if(abs(value-target).gt.tolerance) then
        write(iOut,*)' label     = ',TRIM(label)
        write(iOut,*)' value     = ',value
        write(iOut,*)' target    = ',target
        write(iOut,*)' tolerance = ',tolerance
        call mqc_error('unitTest11: scalar comparison failed.')
      endIf
!
      return
      end subroutine assertNear

!PROCEDURE assertVectorNear
      subroutine assertVectorNear(value,target,tolerance,label)
      implicit none
      real(kind=real64),dimension(:),intent(in)::value,target
      real(kind=real64),intent(in)::tolerance
      character(len=*),intent(in)::label
!
      if(Size(value).ne.Size(target))  &
        call mqc_error('unitTest11: vector comparison shape mismatch.')
      if(MAXVAL(abs(value-target)).gt.tolerance) then
        write(iOut,*)' label     = ',TRIM(label)
        write(iOut,*)' max error = ',MAXVAL(abs(value-target))
        write(iOut,*)' tolerance = ',tolerance
        call mqc_error('unitTest11: vector comparison failed.')
      endIf
!
      return
      end subroutine assertVectorNear

!PROCEDURE assertMatrixNear
      subroutine assertMatrixNear(value,target,tolerance,label)
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::value,target
      real(kind=real64),intent(in)::tolerance
      character(len=*),intent(in)::label
!
      if(ANY(shape(value).ne.shape(target)))  &
        call mqc_error('unitTest11: matrix comparison shape mismatch.')
      if(MAXVAL(abs(value-target)).gt.tolerance) then
        write(iOut,*)' label     = ',TRIM(label)
        write(iOut,*)' max error = ',MAXVAL(abs(value-target))
        write(iOut,*)' tolerance = ',tolerance
        call mqc_error('unitTest11: matrix comparison failed.')
      endIf
!
      return
      end subroutine assertMatrixNear

!PROCEDURE assertIntegerVectorEqual
      subroutine assertIntegerVectorEqual(value,target,label)
      implicit none
      integer(kind=int64),dimension(:),intent(in)::value,target
      character(len=*),intent(in)::label
!
      if(Size(value).ne.Size(target).or.ANY(value.ne.target))  &
        call mqc_error('unitTest11: '//TRIM(label)//' failed.')
!
      return
      end subroutine assertIntegerVectorEqual

      end program unitTest11
