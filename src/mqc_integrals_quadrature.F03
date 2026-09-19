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
!     **                           Version 26.9.3                          **
!     **                         September 19, 2026                        **
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
!     Quadrature implementation included in MQC_Integrals after CONTAINS.
!     Owning grids, explicit rules, original-Becke partition, and Lebedev tables.
!     No standalone module or separate compilation unit.
!
!     H. P. Hratchian, 2026.
!
!PROCEDURE gridClear
      subroutine gridClear(self)
!     Dynamic dispatch through the base also clears molecular provenance.
      implicit none
      class(quadrature_grid),intent(inout)::self
      type(quadrature_grid)::empty
      type(molecular_quadrature_grid)::emptyMolecular
      select type(self)
      type is(quadrature_grid)
        self = empty
      type is(molecular_quadrature_grid)
        self = emptyMolecular
      class default
        call mqc_error('gridClear: unsupported extension.')
      end select
      end subroutine gridClear

!PROCEDURE gridIsInitialized
      logical function gridIsInitialized(self) result(ready)
      implicit none
      class(quadrature_grid),intent(in)::self
      ready = allocated(self%points).and.allocated(self%weights)
      end function gridIsInitialized

!PROCEDURE gridPointCount
      function gridPointCount(self) result(n)
      implicit none
      class(quadrature_grid),intent(in)::self
      integer(int64)::n
      n = 0
      if(allocated(self%weights)) n = size(self%weights,kind=int64)
      end function gridPointCount

!PROCEDURE requireGrid
      subroutine requireGrid(self)
      implicit none
      class(quadrature_grid),intent(in)::self
      if(.not.self%isInitialized()) call mqc_error('quadrature_grid: uninitialized grid.')
      end subroutine requireGrid

!PROCEDURE checkGridData
      subroutine checkGridData(points,weights,integrationDimension,coordinateUnits,weightUnits,rule,normalization)
!     Negative weights are valid (including some Lebedev rules).
      implicit none
      real(real64),intent(in)::points(:,:),weights(:)
      integer(int64),intent(in)::integrationDimension
      character(*),intent(in)::coordinateUnits,weightUnits,rule,normalization
      if(size(points,1).lt.1.or.size(weights).lt.1.or.size(points,2).ne.size(weights)) &
        call mqc_error('quadrature_grid: invalid point/weight shape.')
      if(integrationDimension.lt.1.or.integrationDimension.gt.size(points,1)) &
        call mqc_error('quadrature_grid: invalid integration dimension.')
      if(.not.all(ieee_is_finite(points)).or..not.all(ieee_is_finite(weights))) &
        call mqc_error('quadrature_grid: nonfinite point or weight.')
      if(len_trim(coordinateUnits).eq.0.or.len_trim(weightUnits).eq.0.or. &
        len_trim(rule).eq.0.or.len_trim(normalization).eq.0) &
        call mqc_error('quadrature_grid: metadata must be explicit.')
      end subroutine checkGridData

!PROCEDURE gridTakeArrays
      subroutine gridTakeArrays(grid,points,weights,integrationDimension,coordinateUnits,weightUnits,rule,normalization)
!     Transfer ownership without copying. Input allocatables become unallocated.
!     Exact base type prevents overwriting only the parent of a molecular grid.
      implicit none
      type(quadrature_grid),intent(out)::grid
      real(real64),allocatable,intent(inout)::points(:,:),weights(:)
      integer(int64),intent(in)::integrationDimension
      character(*),intent(in)::coordinateUnits,weightUnits,rule,normalization
      if(.not.allocated(points).or..not.allocated(weights)) &
        call mqc_error('gridTakeArrays: arrays must be allocated.')
      call checkGridData(points,weights,integrationDimension,coordinateUnits,weightUnits,rule,normalization)
      call move_alloc(points,grid%points)
      call move_alloc(weights,grid%weights)
      grid%integrationDimension = integrationDimension
      grid%coordinateUnits = trim(coordinateUnits)
      grid%weightUnits = trim(weightUnits)
      grid%rule = trim(rule)
      grid%normalization = trim(normalization)
      end subroutine gridTakeArrays

!PROCEDURE gridFromArrays
      subroutine gridFromArrays(grid,points,weights,integrationDimension,coordinateUnits,weightUnits,rule,normalization)
!     Explicit copying import for caller-owned arrays. No exactness is inferred.
      implicit none
      type(quadrature_grid),intent(out)::grid
      real(real64),intent(in)::points(:,:),weights(:)
      integer(int64),intent(in)::integrationDimension
      character(*),intent(in)::coordinateUnits,weightUnits,rule,normalization
      real(real64),allocatable::p(:,:),w(:)
      call checkGridData(points,weights,integrationDimension,coordinateUnits,weightUnits,rule,normalization)
      p = points
      w = weights
      call gridTakeArrays(grid,p,w,integrationDimension,coordinateUnits,weightUnits,rule,normalization)
      end subroutine gridFromArrays

!PROCEDURE gridGetMetadata
      subroutine gridGetMetadata(self,coordinateDimension,integrationDimension,coordinateUnits,weightUnits, &
        rule,normalization,exactnessDegree)
!     Allocatable strings returned here are independent copies.
      implicit none
      class(quadrature_grid),intent(in)::self
      integer(int64),optional,intent(out)::coordinateDimension,integrationDimension,exactnessDegree
      character(:),allocatable,optional,intent(out)::coordinateUnits,weightUnits,rule,normalization
      call requireGrid(self)
      if(present(coordinateDimension)) coordinateDimension = size(self%points,1,kind=int64)
      if(present(integrationDimension)) integrationDimension = self%integrationDimension
      if(present(exactnessDegree)) exactnessDegree = self%exactnessDegree
      if(present(coordinateUnits)) coordinateUnits = self%coordinateUnits
      if(present(weightUnits)) weightUnits = self%weightUnits
      if(present(rule)) rule = self%rule
      if(present(normalization)) normalization = self%normalization
      end subroutine gridGetMetadata

!PROCEDURE checkBlock
      subroutine checkBlock(self,first,n)
      implicit none
      class(quadrature_grid),intent(in)::self
      integer(int64),intent(in)::first,n
      call requireGrid(self)
      if(first.lt.1.or.n.lt.1.or.first.gt.self%pointCount()) &
        call mqc_error('quadrature_grid: invalid block range.')
      if(n.gt.self%pointCount()-first+1) call mqc_error('quadrature_grid: block exceeds grid.')
      end subroutine checkBlock

!PROCEDURE gridCopyBlock
      subroutine gridCopyBlock(self,first,points,weights)
!     Caller-owned buffers permit bounded reusable scratch; no aliases escape.
!     Number of requested points is the second extent of points.
      implicit none
      class(quadrature_grid),intent(in)::self
      integer(int64),intent(in)::first
      real(real64),intent(out)::points(:,:)
      real(real64),optional,intent(out)::weights(:)
      integer(int64)::n
      n = size(points,2,kind=int64)
      call checkBlock(self,first,n)
      if(size(points,1).ne.size(self%points,1)) call mqc_error('gridCopyBlock: coordinate shape mismatch.')
      points = self%points(:,first:first+n-1)
      if(present(weights)) then
        if(size(weights).ne.n) call mqc_error('gridCopyBlock: weight shape mismatch.')
        weights = self%weights(first:first+n-1)
      endIf
      end subroutine gridCopyBlock

!PROCEDURE gridIntegrate
      function gridIntegrate(self,values) result(integral)
!     Values include the integrand only; weights already include the measure.
      implicit none
      class(quadrature_grid),intent(in)::self
      real(real64),intent(in)::values(:)
      real(real64)::integral
      call requireGrid(self)
      if(size(values).ne.self%pointCount()) call mqc_error('gridIntegrate: value count mismatch.')
      if(.not.all(ieee_is_finite(values))) call mqc_error('gridIntegrate: nonfinite values.')
      integral = dot_product(self%weights,values)
      end function gridIntegrate

!PROCEDURE gridOverlapMatrix
      subroutine gridOverlapMatrix(self,values,overlap)
!     Real sampled functions are columns, points are rows. Preserve the
!     established dot_product summation and symmetry of quadratureTest.
      implicit none
      class(quadrature_grid),intent(in)::self
      real(real64),intent(in)::values(:,:)
      real(real64),allocatable,intent(out)::overlap(:,:)
      integer(int64)::i,j,n
      call requireGrid(self)
      if(size(values,1).ne.self%pointCount().or.size(values,2).lt.1) &
        call mqc_error('gridOverlapMatrix: value shape mismatch.')
      if(.not.all(ieee_is_finite(values))) call mqc_error('gridOverlapMatrix: nonfinite values.')
      n = size(values,2,kind=int64)
      allocate(overlap(n,n))
      do j = 1,n
        do i = 1,j
          overlap(i,j) = dot_product(self%weights,values(:,i)*values(:,j))
          overlap(j,i) = overlap(i,j)
        endDo
      endDo
      end subroutine gridOverlapMatrix

!PROCEDURE gridValidate
      subroutine gridValidate(self)
      implicit none
      class(quadrature_grid),intent(in)::self
      call requireGrid(self)
      call checkGridData(self%points,self%weights,self%integrationDimension, &
        self%coordinateUnits,self%weightUnits,self%rule,self%normalization)
      select type(self)
      class is(molecular_quadrature_grid)
        if(.not.allocated(self%centersBohr).or..not.allocated(self%generatingAtoms)) &
          call mqc_error('gridValidate: missing molecular provenance.')
        if(size(self%generatingAtoms).ne.self%pointCount()) &
          call mqc_error('gridValidate: molecular point count mismatch.')
      end select
      end subroutine gridValidate

!PROCEDURE gridTrapezoid1D
      subroutine gridTrapezoid1D(grid,n,step,origin,units)
      implicit none
      type(quadrature_grid),intent(out)::grid
      integer(int64),intent(in)::n
      real(real64),intent(in)::step,origin
      character(*),intent(in)::units
      real(real64),allocatable::p(:,:),w(:)
      if(n.lt.2) call mqc_error('gridTrapezoid1D: count must be at least two.')
      if(.not.ieee_is_finite(step).or..not.ieee_is_finite(origin)) &
        call mqc_error('gridTrapezoid1D: nonfinite parameter.')
      allocate(p(1,n),w(n))
      call setup_quadrature_trapezoid1d(n,step,origin,p(1,:),w)
      call gridTakeArrays(grid,p,w,1_int64,units,units,'trapezoid','line measure')
      grid%exactnessDegree = 1
      end subroutine gridTrapezoid1D

!PROCEDURE gridGaussLegendre1D
      subroutine gridGaussLegendre1D(grid,n,lower,upper,units)
      implicit none
      type(quadrature_grid),intent(out)::grid
      integer(int64),intent(in)::n
      real(real64),intent(in)::lower,upper
      character(*),intent(in)::units
      real(real64),allocatable::p(:,:),w(:)
      if(n.lt.1) call mqc_error('gridGaussLegendre1D: count must be positive.')
      if(.not.ieee_is_finite(lower).or..not.ieee_is_finite(upper)) &
        call mqc_error('gridGaussLegendre1D: nonfinite bound.')
      allocate(p(1,n),w(n))
      call setup_quadrature_gauss_legendre1d(n,lower,upper,p(1,:),w)
      call gridTakeArrays(grid,p,w,1_int64,units,units,'gauss-legendre','line measure')
      grid%exactnessDegree = 2*n-1
      end subroutine gridGaussLegendre1D

!PROCEDURE gridMappedRadial
      subroutine gridMappedRadial(grid,n,scale,units)
      implicit none
      type(quadrature_grid),intent(out)::grid
      integer(int64),intent(in)::n
      real(real64),intent(in)::scale
      character(*),intent(in)::units
      real(real64),allocatable::p(:,:),w(:)
      if(n.lt.1) call mqc_error('gridMappedRadial: count must be positive.')
      if(.not.ieee_is_finite(scale)) call mqc_error('gridMappedRadial: nonfinite scale.')
      allocate(p(1,n),w(n))
      call setup_quadrature_gauss_legendre_semi_infinite(n,scale,p(1,:),w)
      call gridTakeArrays(grid,p,w,1_int64,units,units,'mapped-gauss-legendre','dr; no radial Jacobian')
      end subroutine gridMappedRadial

!PROCEDURE gridCartesian
      subroutine gridCartesian(grid,n,step,origin,units)
      implicit none
      type(quadrature_grid),intent(out)::grid
      integer(int64),intent(in)::n
      real(real64),intent(in)::step,origin(3)
      character(*),intent(in)::units
      real(real64),allocatable::p(:,:),w(:)
      if(n.lt.2) call mqc_error('gridCartesian: count must be at least two.')
      if(n.gt.huge(n)/n) call mqc_error('gridCartesian: point count overflow.')
      if(n*n.gt.huge(n)/n) call mqc_error('gridCartesian: point count overflow.')
      if(.not.ieee_is_finite(step).or..not.all(ieee_is_finite(origin))) &
        call mqc_error('gridCartesian: nonfinite parameter.')
      allocate(p(3,n**3),w(n**3))
      call setup_quadrature_trapezoid3d(n,step,origin,p,w)
      call gridTakeArrays(grid,p,w,3_int64,units,trim(units)//'**3','cartesian-trapezoid','volume measure')
      grid%exactnessDegree = 1
      end subroutine gridCartesian

!PROCEDURE gridLebedev
      subroutine gridLebedev(grid,n)
      implicit none
      type(quadrature_grid),intent(out)::grid
      integer(int64),intent(in)::n
      real(real64),allocatable::p(:,:),w(:)
      call setup_lebedev_angular_grid(n,p,w)
      call gridTakeArrays(grid,p,w,2_int64,'dimensionless','steradian','lebedev-laikov','solid angle; sum=4*pi')
      grid%exactnessDegree = lebedevRuleDegree(n)
      end subroutine gridLebedev

!PROCEDURE gridFiniteSpherical
      subroutine gridFiniteSpherical(grid,originBohr,rMaxBohr,nRadial,nAngular,radialRule)
      implicit none
      type(quadrature_grid),intent(out)::grid
      integer(int64),intent(in)::nRadial,nAngular,radialRule
      real(real64),intent(in)::originBohr(3),rMaxBohr
      real(real64),allocatable::p(:,:),w(:)
      if(.not.all(ieee_is_finite(originBohr)).or..not.ieee_is_finite(rMaxBohr)) &
        call mqc_error('gridFiniteSpherical: nonfinite parameter.')
      call setup_quadrature_lebedev_spherical(originBohr,rMaxBohr,nRadial,nAngular,radialRule,p,w)
      call gridTakeArrays(grid,p,w,3_int64,'bohr','bohr**3','finite-lebedev-spherical','volume measure')
      end subroutine gridFiniteSpherical

!PROCEDURE gridSphericalProduct
      subroutine gridSphericalProduct(grid,origin,radialPoints,radialWeights,angularPoints,angularWeights,units)
      implicit none
      type(quadrature_grid),intent(out)::grid
      real(real64),intent(in)::origin(3),radialPoints(:),radialWeights(:),angularPoints(:,:),angularWeights(:)
      character(*),intent(in)::units
      real(real64),allocatable::p(:,:),w(:)
      call checkRules(radialPoints,radialWeights,angularPoints,angularWeights)
      if(.not.all(ieee_is_finite(origin))) call mqc_error('gridSphericalProduct: nonfinite origin.')
      call setup_quadrature_spherical_product(origin,radialPoints,radialWeights,angularPoints,angularWeights,p,w)
      call gridTakeArrays(grid,p,w,3_int64,units,trim(units)//'**3','supplied-spherical-product','volume measure')
      end subroutine gridSphericalProduct

!PROCEDURE checkRules
      subroutine checkRules(radialPoints,radialWeights,angularPoints,angularWeights)
      implicit none
      real(real64),intent(in)::radialPoints(:),radialWeights(:),angularPoints(:,:),angularWeights(:)
      integer(int64)::nr,na
      nr = size(radialPoints,kind=int64)
      na = size(angularWeights,kind=int64)
      if(nr.lt.1.or.na.lt.1) call mqc_error('quadrature_grid: empty supplied rule.')
      if(nr.gt.huge(nr)/na) call mqc_error('quadrature_grid: product count overflow.')
      if(.not.all(ieee_is_finite(radialPoints)).or..not.all(ieee_is_finite(radialWeights)).or. &
        .not.all(ieee_is_finite(angularPoints)).or..not.all(ieee_is_finite(angularWeights))) &
        call mqc_error('quadrature_grid: nonfinite supplied rule.')
      end subroutine checkRules

!PROCEDURE gridBeckeMolecular
      subroutine gridBeckeMolecular(grid,centersBohr,radiiBohr,radialPointsBohr,radialWeightsBohr, &
        angularPoints,angularWeights)
!     Retain the exact supplied atomic rules and radii for reproducibility.
!     Explicit radii keep this constructor independent of element lookup.
      implicit none
      type(molecular_quadrature_grid),intent(out)::grid
      real(real64),intent(in)::centersBohr(:,:),radiiBohr(:),radialPointsBohr(:),radialWeightsBohr(:)
      real(real64),intent(in)::angularPoints(:,:),angularWeights(:)
      integer(int64)::nAtomic
      call checkRules(radialPointsBohr,radialWeightsBohr,angularPoints,angularWeights)
      nAtomic = size(radialPointsBohr,kind=int64)*size(angularWeights,kind=int64)
      if(size(centersBohr,2,kind=int64).gt.huge(nAtomic)/nAtomic) &
        call mqc_error('gridBeckeMolecular: point count overflow.')
      if(.not.all(ieee_is_finite(centersBohr)).or..not.all(ieee_is_finite(radiiBohr))) &
        call mqc_error('gridBeckeMolecular: nonfinite molecular data.')
      call setup_quadrature_becke_molecular(centersBohr,radiiBohr,radialPointsBohr,radialWeightsBohr, &
        angularPoints,angularWeights,grid%points,grid%generatingAtoms,grid%baseWeights, &
        grid%partitionWeights,grid%weights)
      grid%centersBohr = centersBohr
      grid%radiiBohr = radiiBohr
      grid%radialPoints = radialPointsBohr
      grid%radialWeights = radialWeightsBohr
      grid%angularPoints = angularPoints
      grid%angularWeights = angularWeights
      grid%integrationDimension = 3
      grid%coordinateUnits = 'bohr'
      grid%weightUnits = 'bohr**3'
      grid%rule = 'becke-1988; unpruned; supplied atomic rules'
      grid%normalization = 'volume measure; base weight * partition weight'
      call grid%validate()
      end subroutine gridBeckeMolecular

!PROCEDURE molecularAtomCount
      function molecularAtomCount(self) result(n)
      implicit none
      class(molecular_quadrature_grid),intent(in)::self
      integer(int64)::n
      n = 0
      if(allocated(self%centersBohr)) n = size(self%centersBohr,2,kind=int64)
      end function molecularAtomCount

!PROCEDURE molecularCopyMolecule
      subroutine molecularCopyMolecule(self,centersBohr,radiiBohr)
      implicit none
      class(molecular_quadrature_grid),intent(in)::self
      real(real64),allocatable,intent(out)::centersBohr(:,:),radiiBohr(:)
      call requireGrid(self)
      centersBohr = self%centersBohr
      radiiBohr = self%radiiBohr
      end subroutine molecularCopyMolecule

!PROCEDURE molecularCopyAtomicRule
      subroutine molecularCopyAtomicRule(self,radialPointsBohr,radialWeightsBohr,angularPoints,angularWeights)
      implicit none
      class(molecular_quadrature_grid),intent(in)::self
      real(real64),allocatable,intent(out)::radialPointsBohr(:),radialWeightsBohr(:),angularPoints(:,:),angularWeights(:)
      call requireGrid(self)
      radialPointsBohr = self%radialPoints
      radialWeightsBohr = self%radialWeights
      angularPoints = self%angularPoints
      angularWeights = self%angularWeights
      end subroutine molecularCopyAtomicRule

!PROCEDURE molecularCopyPartitionBlock
      subroutine molecularCopyPartitionBlock(self,first,generatingAtoms,baseWeightsBohr3,partitionWeights)
      implicit none
      class(molecular_quadrature_grid),intent(in)::self
      integer(int64),intent(in)::first
      integer(int64),intent(out)::generatingAtoms(:)
      real(real64),intent(out)::baseWeightsBohr3(:),partitionWeights(:)
      integer(int64)::n
      n = size(generatingAtoms,kind=int64)
      call checkBlock(self,first,n)
      if(size(baseWeightsBohr3).ne.n.or.size(partitionWeights).ne.n) &
        call mqc_error('copyPartitionBlock: metadata shape mismatch.')
      generatingAtoms = self%generatingAtoms(first:first+n-1)
      baseWeightsBohr3 = self%baseWeights(first:first+n-1)
      partitionWeights = self%partitionWeights(first:first+n-1)
      end subroutine molecularCopyPartitionBlock

!PROCEDURE molecularPartitionRange
      subroutine molecularPartitionRange(self,minimum,maximum)
      implicit none
      class(molecular_quadrature_grid),intent(in)::self
      real(real64),intent(out)::minimum,maximum
      call requireGrid(self)
      minimum = minval(self%partitionWeights)
      maximum = maxval(self%partitionWeights)
      end subroutine molecularPartitionRange

!     Quadrature procedures.

!PROCEDURE setup_quadrature_trapezoid1d
      subroutine setup_quadrature_trapezoid1d(nPoints,stepSize,origin,  &
        quadGrid,quadWeights)
!
!     Generate a uniform one-dimensional trapezoid grid. The first point is
!     origin, the interval length is (nPoints-1)*stepSize, and the two endpoint
!     weights are one half of the interior weights.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::nPoints
      real(kind=real64),intent(in)::stepSize,origin
      real(kind=real64),dimension(nPoints),intent(out)::quadGrid,  &
        quadWeights
      integer(kind=int64)::i
!
      if(nPoints.lt.2)  &
        call mqc_error(  &
          'setup_quadrature_trapezoid1d: nPoints must be at least 2.')
      if(stepSize.le.mqc_float(0))  &
        call mqc_error(  &
          'setup_quadrature_trapezoid1d: stepSize must be positive.')
!
      quadGrid(1) = origin
      quadWeights(1) = stepSize/mqc_float(2)
      do i = 2,nPoints-1
        quadGrid(i) = quadGrid(i-1)+stepSize
        quadWeights(i) = stepSize
      endDo
      quadGrid(nPoints) = quadGrid(nPoints-1)+stepSize
      quadWeights(nPoints) = stepSize/mqc_float(2)
!
      return
      end subroutine setup_quadrature_trapezoid1d

!PROCEDURE setup_quadrature_trapezoid3d
      subroutine setup_quadrature_trapezoid3d(nPoints,stepSize,origin,  &
        quadGrid,quadWeights)
!
!     Generate the Cartesian product of three identical uniform trapezoid
!     grids. origin is the lower corner, coordinates use the caller's length
!     unit, and the product weights use the cube of that unit.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::nPoints
      real(kind=real64),intent(in)::stepSize
      real(kind=real64),dimension(3),intent(in)::origin
      real(kind=real64),dimension(3,nPoints**3),intent(out)::quadGrid
      real(kind=real64),dimension(nPoints**3),intent(out)::quadWeights
      integer(kind=int64)::i,j,k,idx
      real(kind=real64)::wx,wy,wz,x,y,z
!
      if(nPoints.lt.2)  &
        call mqc_error(  &
          'setup_quadrature_trapezoid3d: nPoints must be at least 2.')
      if(stepSize.le.mqc_float(0))  &
        call mqc_error(  &
          'setup_quadrature_trapezoid3d: stepSize must be positive.')
!
      idx = 0
      do k = 1,nPoints
        z = origin(3)+mqc_float(k-1)*stepSize
        wz = stepSize
        if(k.eq.1.or.k.eq.nPoints) wz = stepSize/mqc_float(2)
        do j = 1,nPoints
          y = origin(2)+mqc_float(j-1)*stepSize
          wy = stepSize
          if(j.eq.1.or.j.eq.nPoints) wy = stepSize/mqc_float(2)
          do i = 1,nPoints
            x = origin(1)+mqc_float(i-1)*stepSize
            wx = stepSize
            if(i.eq.1.or.i.eq.nPoints) wx = stepSize/mqc_float(2)
            idx = idx+1
            quadGrid(:,idx) = [ x,y,z ]
            quadWeights(idx) = wx*wy*wz
          endDo
        endDo
      endDo
!
      return
      end subroutine setup_quadrature_trapezoid3d

!PROCEDURE setup_quadrature_gauss_legendre1d
      subroutine setup_quadrature_gauss_legendre1d(nPoints,lowerBound,  &
        upperBound,quadGrid,quadWeights)
!
!     Generate an nPoints Gauss-Legendre grid on the finite interval
!     [lowerBound,upperBound].
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::nPoints
      real(kind=real64),intent(in)::lowerBound,upperBound
      real(kind=real64),dimension(nPoints),intent(out)::quadGrid,  &
        quadWeights
      integer(kind=int64)::i,j,m,nIterations
      real(kind=real64)::p1,p2,p3,pp,tol,xMid,xHalf,z,zOld
      logical::converged
!
      if(nPoints.lt.1)  &
        call mqc_error(  &
          'setup_quadrature_gauss_legendre1d: nPoints must be positive.')
      if(upperBound.le.lowerBound)  &
        call mqc_error(  &
          'setup_quadrature_gauss_legendre1d: invalid interval.')
!
      tol = 1.0e-14_real64
      m = (nPoints+1)/2
      xMid = (upperBound+lowerBound)/mqc_float(2)
      xHalf = (upperBound-lowerBound)/mqc_float(2)
      do i = 1,m
        z = cos(Pi*(mqc_float(i)-mqc_float(1)/mqc_float(4))/  &
          (mqc_float(nPoints)+mqc_float(1)/mqc_float(2)))
        converged = .false.
        do nIterations = 1,100
          p1 = mqc_float(1)
          p2 = mqc_float(0)
          do j = 1,nPoints
            p3 = p2
            p2 = p1
            p1 = ((mqc_float(2*j-1)*z*p2)-mqc_float(j-1)*p3)/  &
              mqc_float(j)
          endDo
          pp = mqc_float(nPoints)*(z*p1-p2)/(z*z-mqc_float(1))
          zOld = z
          z = zOld-p1/pp
          if(abs(z-zOld).le.tol) then
            converged = .true.
            exit
          endIf
        endDo
        if(.not.converged)  &
          call mqc_error(  &
            'setup_quadrature_gauss_legendre1d: root iteration failed.')
        quadGrid(i) = xMid-xHalf*z
        quadGrid(nPoints+1-i) = xMid+xHalf*z
        quadWeights(i) = mqc_float(2)*xHalf/  &
          ((mqc_float(1)-z*z)*pp*pp)
        quadWeights(nPoints+1-i) = quadWeights(i)
      endDo
!
      return
      end subroutine setup_quadrature_gauss_legendre1d

!PROCEDURE setup_quadrature_gauss_legendre_semi_infinite
      subroutine setup_quadrature_gauss_legendre_semi_infinite(nPoints,  &
        radialScale,quadGrid,quadWeights)
!
!     Map Gauss-Legendre quadrature from -1 < x < 1 onto 0 < r < infinity
!     with r=radialScale*(1+x)/(1-x). The returned weights integrate dr;
!     spherical-product construction applies the separate r**2 Jacobian.
!     radialScale and the returned grid use the caller's length unit.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::nPoints
      real(kind=real64),intent(in)::radialScale
      real(kind=real64),dimension(nPoints),intent(out)::quadGrid,  &
        quadWeights
      integer(kind=int64)::i
      real(kind=real64)::denominator
      real(kind=real64),dimension(nPoints)::mappedGrid,mappedWeights
!
      if(nPoints.lt.1)  &
        call mqc_error(  &
          'setup_quadrature_gauss_legendre_semi_infinite: nPoints must be positive.')
      if(radialScale.le.mqc_float(0))  &
        call mqc_error(  &
          'setup_quadrature_gauss_legendre_semi_infinite: radialScale must be positive.')
!
      call setup_quadrature_gauss_legendre1d(nPoints,-mqc_float(1),  &
        mqc_float(1),mappedGrid,mappedWeights)
      do i = 1,nPoints
        denominator = mqc_float(1)-mappedGrid(i)
        quadGrid(i) = radialScale*(mqc_float(1)+mappedGrid(i))/  &
          denominator
        quadWeights(i) = mappedWeights(i)*mqc_float(2)*radialScale/  &
          (denominator*denominator)
      endDo
!
      return
      end subroutine setup_quadrature_gauss_legendre_semi_infinite



!PROCEDURE setup_quadrature_spherical_product
      subroutine setup_quadrature_spherical_product(origin,radialGrid,  &
        radialWeights,angularGrid,angularWeights,quadGrid,quadWeights)
!
!     Build a one-center Cartesian point/weight list from caller-supplied
!     radial and unit-sphere angular rules. radialGrid and origin must use the
!     same length unit. The volume weights are radialWeight*r**2*angularWeight
!     and therefore use the cube of that length unit.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(3),intent(in)::origin
      real(kind=real64),dimension(:),intent(in)::radialGrid,radialWeights,  &
        angularWeights
      real(kind=real64),dimension(:,:),intent(in)::angularGrid
      real(kind=real64),dimension(:,:),allocatable,intent(out)::quadGrid
      real(kind=real64),dimension(:),allocatable,intent(out)::quadWeights
      integer(kind=int64)::iAngular,iRadial,idx,nAngular,nRadial
      real(kind=real64)::angularNorm
!
      nRadial = Size(radialGrid)
      nAngular = Size(angularWeights)
      if(nRadial.lt.1)  &
        call mqc_error(  &
          'setup_quadrature_spherical_product: no radial points supplied.')
      if(Size(radialWeights).ne.nRadial)  &
        call mqc_error(  &
          'setup_quadrature_spherical_product: radial grid/weight sizes differ.')
      if(MINVAL(radialGrid).lt.mqc_float(0))  &
        call mqc_error(  &
          'setup_quadrature_spherical_product: negative radius supplied.')
      if(Size(angularGrid,1).ne.3)  &
        call mqc_error(  &
          'setup_quadrature_spherical_product: angular grid must be 3 by n.')
      if(nAngular.lt.1)  &
        call mqc_error(  &
          'setup_quadrature_spherical_product: no angular points supplied.')
      if(Size(angularGrid,2).ne.nAngular)  &
        call mqc_error(  &
          'setup_quadrature_spherical_product: angular grid/weight sizes differ.')
      do iAngular = 1,nAngular
        angularNorm = sqrt(dot_product(angularGrid(:,iAngular),  &
          angularGrid(:,iAngular)))
        if(abs(angularNorm-mqc_float(1)).gt.1.0e-10_real64)  &
          call mqc_error(  &
            'setup_quadrature_spherical_product: angular direction is not unit length.')
      endDo
!
      Allocate(quadGrid(3,nRadial*nAngular),  &
        quadWeights(nRadial*nAngular))
      idx = 0
      do iRadial = 1,nRadial
        do iAngular = 1,nAngular
          idx = idx+1
          quadGrid(:,idx) = origin+radialGrid(iRadial)*  &
            angularGrid(:,iAngular)
          quadWeights(idx) = radialWeights(iRadial)*  &
            radialGrid(iRadial)*radialGrid(iRadial)*  &
            angularWeights(iAngular)
        endDo
      endDo
!
      return
      end subroutine setup_quadrature_spherical_product

!PROCEDURE setup_quadrature_lebedev_spherical
      subroutine setup_quadrature_lebedev_spherical(originBohr,rMaxBohr,  &
        nRadial,lebedevOrder,radialRule,quadGrid,quadWeights)
!
!     Build the existing finite one-center spherical product grid by preparing
!     a trapezoid or Gauss-Legendre radial rule and a selected Lebedev angular
!     rule, then delegating their Cartesian product to the general constructor.
!     Coordinates are in bohr and weights are in bohr**3.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(3),intent(in)::originBohr
      real(kind=real64),intent(in)::rMaxBohr
      integer(kind=int64),intent(in)::nRadial,lebedevOrder,radialRule
      real(kind=real64),dimension(:,:),allocatable,intent(out)::quadGrid
      real(kind=real64),dimension(:),allocatable,intent(out)::quadWeights
      real(kind=real64)::stepSize
      real(kind=real64),dimension(:),allocatable::angularWeights,  &
        radialGrid,radialWeights
      real(kind=real64),dimension(:,:),allocatable::angularGrid
!
      if(rMaxBohr.le.mqc_float(0))  &
        call mqc_error(  &
          'setup_quadrature_lebedev_spherical: rMaxBohr must be positive.')
      if(nRadial.lt.1)  &
        call mqc_error(  &
          'setup_quadrature_lebedev_spherical: nRadial must be positive.')
      if(radialRule.eq.QUAD_RADIAL_TRAPEZOID.and.nRadial.lt.2)  &
        call mqc_error(  &
          'setup_quadrature_lebedev_spherical: trapezoid radial grid needs nRadial >= 2.')
!
      Allocate(radialGrid(nRadial),radialWeights(nRadial))
      select case(radialRule)
      case(QUAD_RADIAL_TRAPEZOID)
        stepSize = rMaxBohr/mqc_float(nRadial-1)
        call setup_quadrature_trapezoid1d(nRadial,stepSize,  &
          mqc_float(0),radialGrid,radialWeights)
      case(QUAD_RADIAL_GAUSS_LEGENDRE)
        call setup_quadrature_gauss_legendre1d(nRadial,mqc_float(0),  &
          rMaxBohr,radialGrid,radialWeights)
      case default
        call mqc_error(  &
          'setup_quadrature_lebedev_spherical: unsupported radial rule.')
      end select
!
      call setup_lebedev_angular_grid(lebedevOrder,angularGrid,  &
        angularWeights)
      call setup_quadrature_spherical_product(originBohr,radialGrid,  &
        radialWeights,angularGrid,angularWeights,quadGrid,quadWeights)
!
      return
      end subroutine setup_quadrature_lebedev_spherical

!     Molecular Grid procedures.

!PROCEDURE prepareBecke1988PairData
      subroutine prepareBecke1988PairData(centerCoordsBohr,atomicRadiiBohr,  &
        pairDistancesBohr,sizeAdjustments)
!
!     Validate the molecular centers and original-Becke atomic radii, then
!     prepare the internuclear distances and antisymmetric heteronuclear size
!     adjustments. Coordinates and radii are in bohr.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::centerCoordsBohr
      real(kind=real64),dimension(:),intent(in)::atomicRadiiBohr
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        pairDistancesBohr,sizeAdjustments
      integer(kind=int64)::iAtom,jAtom,nAtoms
      real(kind=real64)::adjustment,radiusRatio
!
      if(Size(centerCoordsBohr,1).ne.3)  &
        call mqc_error(  &
          'prepareBecke1988PairData: centerCoordsBohr must be 3 by nAtoms.')
      nAtoms = Size(centerCoordsBohr,2)
      if(nAtoms.lt.1)  &
        call mqc_error(  &
          'prepareBecke1988PairData: at least one atom is required.')
      if(Size(atomicRadiiBohr).ne.nAtoms)  &
        call mqc_error(  &
          'prepareBecke1988PairData: center and radius counts differ.')
      if(MINVAL(atomicRadiiBohr).le.mqc_small)  &
        call mqc_error(  &
          'prepareBecke1988PairData: atomic radii must be positive.')
!
      pairDistancesBohr = mqc_distance_matrix(centerCoordsBohr)
      Allocate(sizeAdjustments(nAtoms,nAtoms))
      sizeAdjustments = mqc_float(0)
      do iAtom = 1,nAtoms-1
        do jAtom = iAtom+1,nAtoms
          if(pairDistancesBohr(iAtom,jAtom).le.mqc_small)  &
            call mqc_error(  &
              'prepareBecke1988PairData: atomic centers must be distinct.')
          radiusRatio = atomicRadiiBohr(iAtom)/atomicRadiiBohr(jAtom)
          adjustment = (mqc_float(1)/mqc_float(4))*  &
            (mqc_float(1)/radiusRatio-radiusRatio)
          adjustment = max(-mqc_float(1)/mqc_float(2),  &
            min(mqc_float(1)/mqc_float(2),adjustment))
          sizeAdjustments(iAtom,jAtom) = adjustment
          sizeAdjustments(jAtom,iAtom) = -adjustment
        endDo
      endDo
!
      return
      end subroutine prepareBecke1988PairData

!PROCEDURE becke1988SmoothStep
      function becke1988SmoothStep(value) result(smoothedValue)
!
!     Apply Becke's cubic switching polynomial three times. The input and
!     output are limited to their mathematical interval to suppress roundoff
!     excursions at atomic centers and very distant points.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),intent(in)::value
      real(kind=real64)::smoothedValue
      integer(kind=int64)::i
!
      smoothedValue = max(-mqc_float(1),min(mqc_float(1),value))
      do i = 1,3
        smoothedValue = (smoothedValue/  &
          mqc_float(2))*(mqc_float(3)-smoothedValue*smoothedValue)
      endDo
      smoothedValue = max(-mqc_float(1),  &
        min(mqc_float(1),smoothedValue))
!
      return
      end function becke1988SmoothStep

!PROCEDURE becke1988PartitionWeightsPrepared
      subroutine becke1988PartitionWeightsPrepared(centerCoordsBohr,  &
        pairDistancesBohr,sizeAdjustments,pointsBohr,partitionWeights)
!
!     Evaluate the normalized original-Becke partition for Cartesian points
!     using validated, precomputed atom-pair data. The cell products are
!     accumulated logarithmically so normalization remains stable for many
!     atoms. Coordinates and pair distances are in bohr.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::centerCoordsBohr,  &
        pairDistancesBohr,sizeAdjustments,pointsBohr
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        partitionWeights
      integer(kind=int64)::iAtom,iPoint,jAtom,nAtoms,nPoints
      real(kind=real64)::factorI,factorJ,maxLogWeight,muAdjusted,  &
        muPair,normalization,smoothedMu
      logical,dimension(:),allocatable::zeroCellWeight
      real(kind=real64),dimension(:),allocatable::cellWeights,  &
        logCellWeights,pointDistancesBohr
!
      if(Size(pointsBohr,1).ne.3)  &
        call mqc_error(  &
          'becke1988_partition_weights: pointsBohr must be 3 by nPoints.')
      nAtoms = Size(centerCoordsBohr,2)
      nPoints = Size(pointsBohr,2)
      if(nPoints.lt.1)  &
        call mqc_error(  &
          'becke1988_partition_weights: at least one point is required.')
!
      Allocate(partitionWeights(nAtoms,nPoints))
      if(nAtoms.eq.1) then
        partitionWeights = mqc_float(1)
        return
      endIf
      Allocate(cellWeights(nAtoms),logCellWeights(nAtoms),  &
        pointDistancesBohr(nAtoms),zeroCellWeight(nAtoms))
!
      do iPoint = 1,nPoints
        do iAtom = 1,nAtoms
          pointDistancesBohr(iAtom) = mqc_distance(pointsBohr(:,iPoint),  &
            centerCoordsBohr(:,iAtom))
        endDo
        logCellWeights = mqc_float(0)
        zeroCellWeight = .false.
        do iAtom = 1,nAtoms-1
          do jAtom = iAtom+1,nAtoms
            muPair = (pointDistancesBohr(iAtom)-  &
              pointDistancesBohr(jAtom))/  &
              pairDistancesBohr(iAtom,jAtom)
            muPair = max(-mqc_float(1),min(mqc_float(1),muPair))
            muAdjusted = muPair+sizeAdjustments(iAtom,jAtom)*  &
              (mqc_float(1)-muPair*muPair)
            smoothedMu = becke1988SmoothStep(muAdjusted)
            factorI = (mqc_float(1)-smoothedMu)/mqc_float(2)
            factorJ = mqc_float(1)-factorI
            factorI = max(mqc_float(0),min(mqc_float(1),factorI))
            factorJ = max(mqc_float(0),min(mqc_float(1),factorJ))
            if(factorI.le.mqc_float(0)) then
              zeroCellWeight(iAtom) = .true.
            else
              logCellWeights(iAtom) =  &
                logCellWeights(iAtom)+log(factorI)
            endIf
            if(factorJ.le.mqc_float(0)) then
              zeroCellWeight(jAtom) = .true.
            else
              logCellWeights(jAtom) =  &
                logCellWeights(jAtom)+log(factorJ)
            endIf
          endDo
        endDo
        if(.not.ANY(.not.zeroCellWeight))  &
          call mqc_error(  &
            'becke1988_partition_weights: all cell weights vanished.')
        maxLogWeight = MAXVAL(logCellWeights,mask=.not.zeroCellWeight)
        cellWeights = mqc_float(0)
        do iAtom = 1,nAtoms
          if(.not.zeroCellWeight(iAtom))  &
            cellWeights(iAtom) = exp(logCellWeights(iAtom)-maxLogWeight)
        endDo
        normalization = SUM(cellWeights)
        if(normalization.le.mqc_float(0))  &
          call mqc_error(  &
            'becke1988_partition_weights: invalid normalization.')
        partitionWeights(:,iPoint) = cellWeights/normalization
      endDo
!
      return
      end subroutine becke1988PartitionWeightsPrepared

!PROCEDURE becke1988_partition_weights
      subroutine becke1988_partition_weights(centerCoordsBohr,  &
        atomicRadiiBohr,pointsBohr,partitionWeights)
!
!     Return the normalized original-Becke molecular partition weights for a
!     caller-supplied Cartesian point list. Coordinates and radii are in bohr;
!     partitionWeights(iAtom,iPoint) is dimensionless.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::centerCoordsBohr,  &
        pointsBohr
      real(kind=real64),dimension(:),intent(in)::atomicRadiiBohr
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        partitionWeights
      real(kind=real64),dimension(:,:),allocatable::pairDistancesBohr,  &
        sizeAdjustments
!
      call prepareBecke1988PairData(centerCoordsBohr,atomicRadiiBohr,  &
        pairDistancesBohr,sizeAdjustments)
      call becke1988PartitionWeightsPrepared(centerCoordsBohr,  &
        pairDistancesBohr,sizeAdjustments,pointsBohr,partitionWeights)
!
      return
      end subroutine becke1988_partition_weights

!PROCEDURE apply_becke1988_partition
      subroutine apply_becke1988_partition(centerCoordsBohr,  &
        atomicRadiiBohr,pointsBohr,generatingAtoms,baseWeights,  &
        partitionWeights,finalWeights)
!
!     Apply the partition belonging to each point's generating atom to an
!     existing concatenated atomic grid. The base and final weights retain the
!     units supplied by the caller; partitionWeights is dimensionless.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::centerCoordsBohr,  &
        pointsBohr
      real(kind=real64),dimension(:),intent(in)::atomicRadiiBohr,  &
        baseWeights
      integer(kind=int64),dimension(:),intent(in)::generatingAtoms
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        partitionWeights,finalWeights
      integer(kind=int64)::iPoint,nAtoms,nPoints
      real(kind=real64),dimension(:,:),allocatable::allPartitionWeights
!
      nAtoms = Size(centerCoordsBohr,2)
      nPoints = Size(pointsBohr,2)
      if(nPoints.lt.1)  &
        call mqc_error(  &
          'apply_becke1988_partition: at least one point is required.')
      if(Size(generatingAtoms).ne.nPoints.or.  &
        Size(baseWeights).ne.nPoints)  &
        call mqc_error(  &
          'apply_becke1988_partition: point metadata sizes differ.')
      if(MINVAL(generatingAtoms).lt.1.or.  &
        MAXVAL(generatingAtoms).gt.nAtoms)  &
        call mqc_error(  &
          'apply_becke1988_partition: invalid generating atom index.')
!
      call becke1988_partition_weights(centerCoordsBohr,atomicRadiiBohr,  &
        pointsBohr,allPartitionWeights)
      Allocate(partitionWeights(nPoints),finalWeights(nPoints))
      do iPoint = 1,nPoints
        partitionWeights(iPoint) =  &
          allPartitionWeights(generatingAtoms(iPoint),iPoint)
      endDo
      finalWeights = baseWeights*partitionWeights
!
      return
      end subroutine apply_becke1988_partition

!PROCEDURE setup_quadrature_becke_molecular
      subroutine setup_quadrature_becke_molecular(centerCoordsBohr,  &
        atomicRadiiBohr,radialGridBohr,radialWeightsBohr,angularGrid,  &
        angularWeights,quadGridBohr,generatingAtoms,baseWeightsBohr3,  &
        partitionWeights,quadWeightsBohr3)
!
!     Build an unpruned atom-centered molecular grid. Each atomic block is
!     constructed with the existing spherical-product routine and retained in
!     the concatenated grid without geometric deduplication. The final volume
!     weight is the base atomic weight times the generating atom's original-
!     Becke partition weight.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::centerCoordsBohr,  &
        angularGrid
      real(kind=real64),dimension(:),intent(in)::atomicRadiiBohr,  &
        radialGridBohr,radialWeightsBohr,angularWeights
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        quadGridBohr
      integer(kind=int64),dimension(:),allocatable,intent(out)::  &
        generatingAtoms
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        baseWeightsBohr3,partitionWeights,quadWeightsBohr3
      integer(kind=int64)::iAtom,iEnd,iStart,nAtomicPoints,nAtoms,nPoints
      real(kind=real64),dimension(:),allocatable::atomicBaseWeights,  &
        atomicPartitionWeights
      real(kind=real64),dimension(:,:),allocatable::atomicGrid,  &
        blockPartitionWeights,pairDistancesBohr,sizeAdjustments
!
      call prepareBecke1988PairData(centerCoordsBohr,atomicRadiiBohr,  &
        pairDistancesBohr,sizeAdjustments)
      nAtoms = Size(centerCoordsBohr,2)
      nAtomicPoints = Size(radialGridBohr)*Size(angularWeights)
      if(nAtomicPoints.lt.1)  &
        call mqc_error(  &
          'setup_quadrature_becke_molecular: empty atomic quadrature.')
      nPoints = nAtoms*nAtomicPoints
      Allocate(quadGridBohr(3,nPoints),generatingAtoms(nPoints),  &
        baseWeightsBohr3(nPoints),partitionWeights(nPoints),  &
        quadWeightsBohr3(nPoints))
!
      do iAtom = 1,nAtoms
        call setup_quadrature_spherical_product(  &
          centerCoordsBohr(:,iAtom),radialGridBohr,radialWeightsBohr,  &
          angularGrid,angularWeights,atomicGrid,atomicBaseWeights)
        if(Size(atomicBaseWeights).ne.nAtomicPoints)  &
          call mqc_error(  &
            'setup_quadrature_becke_molecular: atomic grid size changed.')
        iStart = (iAtom-1)*nAtomicPoints+1
        iEnd = iStart+nAtomicPoints-1
        quadGridBohr(:,iStart:iEnd) = atomicGrid
        generatingAtoms(iStart:iEnd) = iAtom
        baseWeightsBohr3(iStart:iEnd) = atomicBaseWeights
        call becke1988PartitionWeightsPrepared(centerCoordsBohr,  &
          pairDistancesBohr,sizeAdjustments,atomicGrid,  &
          blockPartitionWeights)
        atomicPartitionWeights = blockPartitionWeights(iAtom,:)
        partitionWeights(iStart:iEnd) = atomicPartitionWeights
        quadWeightsBohr3(iStart:iEnd) =  &
          atomicBaseWeights*atomicPartitionWeights
      endDo
!
      return
      end subroutine setup_quadrature_becke_molecular

!     Molecular Grid Adapter procedures.

!PROCEDURE becke1988_radii_from_atomic_numbers
      subroutine becke1988_radii_from_atomic_numbers(atomicNumbers,  &
        atomicRadiiBohr)
!
!     Resolve the original-Becke atomic-size radii for true atomic numbers.
!     Ghost centers and elements without a tabulated Bragg-Slater radius are
!     rejected rather than being assigned an implicit physical radius.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),dimension(:),intent(in)::atomicNumbers
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        atomicRadiiBohr
      integer(kind=int64)::iAtom,nAtoms
!
      nAtoms = Size(atomicNumbers)
      if(nAtoms.lt.1)  &
        call mqc_error(  &
          'becke1988_radii_from_atomic_numbers: no atoms supplied.')
      Allocate(atomicRadiiBohr(nAtoms))
      do iAtom = 1,nAtoms
        if(.not.mqc_element_has_bragg_slater_radius(  &
          atomicNumbers(iAtom)))  &
          call mqc_error(  &
            'becke1988_radii_from_atomic_numbers: radius unavailable.')
        atomicRadiiBohr(iAtom) =  &
          mqc_element_becke_1988_radius(atomicNumbers(iAtom))
      endDo
!
      return
      end subroutine becke1988_radii_from_atomic_numbers

!PROCEDURE setup_quadrature_becke_molecular_from_atomic_numbers
      subroutine setup_quadrature_becke_molecular_from_atomic_numbers(  &
        atomicNumbers,centerCoordsBohr,radialGridBohr,radialWeightsBohr,  &
        angularGrid,angularWeights,quadGridBohr,generatingAtoms,  &
        baseWeightsBohr3,partitionWeights,quadWeightsBohr3)
!
!     Build the unpruned original-Becke molecular grid after resolving the
!     atomic-size radii from true atomic numbers through MQC_General.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),dimension(:),intent(in)::atomicNumbers
      real(kind=real64),dimension(:,:),intent(in)::centerCoordsBohr,  &
        angularGrid
      real(kind=real64),dimension(:),intent(in)::radialGridBohr,  &
        radialWeightsBohr,angularWeights
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        quadGridBohr
      integer(kind=int64),dimension(:),allocatable,intent(out)::  &
        generatingAtoms
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        baseWeightsBohr3,partitionWeights,quadWeightsBohr3
      real(kind=real64),dimension(:),allocatable::atomicRadiiBohr
!
      if(Size(atomicNumbers).ne.Size(centerCoordsBohr,2))  &
        call mqc_error(  &
          'setup_quadrature_becke_molecular_from_atomic_numbers: atom counts differ.')
      call becke1988_radii_from_atomic_numbers(atomicNumbers,  &
        atomicRadiiBohr)
      call setup_quadrature_becke_molecular(centerCoordsBohr,  &
        atomicRadiiBohr,radialGridBohr,radialWeightsBohr,angularGrid,  &
        angularWeights,quadGridBohr,generatingAtoms,baseWeightsBohr3,  &
        partitionWeights,quadWeightsBohr3)
!
      return
      end subroutine setup_quadrature_becke_molecular_from_atomic_numbers

!     Lebedev Grid procedures.
!
!     This module generates Lebedev-Laikov angular quadrature rules on
!     the unit sphere. The compact rule coefficients were transcribed
!     from the publicly distributed CCL Lebedev-Laikov Fortran tables:
!
!       https://www.ccl.net/cca/software/SOURCES/FORTRAN/
!         Lebedev-Laikov-Grids/Lebedev-Laikov.F
!
!     The tables originate with V. I. Lebedev and D. N. Laikov,
!     Doklady Mathematics 59, 477-481 (1999). Angular weights returned
!     here include the 4*pi solid-angle normalization.
!
!     H. P. Hratchian, 2026.
!

!PROCEDURE setup_lebedev_angular_grid
      subroutine setup_lebedev_angular_grid(lebedevPoints,angularGrid,  &
        angularWeights)
!
!     This routine builds one supported Lebedev-Laikov angular rule.
!     Points lie on the unit sphere and weights sum to 4*pi.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::lebedevPoints
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        angularGrid
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        angularWeights
      integer(kind=int64)::idx
!
      if(.not.lebedevRuleIsSupported(lebedevPoints))  &
        call mqc_error('Unsupported Lebedev angular point count.')
      Allocate(angularGrid(3,lebedevPoints),angularWeights(lebedevPoints))
      idx = 0_int64
      select case(lebedevPoints)
      case(6_int64)
        call buildLebedevRule0006(angularGrid,angularWeights,idx)
      case(14_int64)
        call buildLebedevRule0014(angularGrid,angularWeights,idx)
      case(26_int64)
        call buildLebedevRule0026(angularGrid,angularWeights,idx)
      case(38_int64)
        call buildLebedevRule0038(angularGrid,angularWeights,idx)
      case(50_int64)
        call buildLebedevRule0050(angularGrid,angularWeights,idx)
      case(74_int64)
        call buildLebedevRule0074(angularGrid,angularWeights,idx)
      case(86_int64)
        call buildLebedevRule0086(angularGrid,angularWeights,idx)
      case(110_int64)
        call buildLebedevRule0110(angularGrid,angularWeights,idx)
      case(146_int64)
        call buildLebedevRule0146(angularGrid,angularWeights,idx)
      case(170_int64)
        call buildLebedevRule0170(angularGrid,angularWeights,idx)
      case(194_int64)
        call buildLebedevRule0194(angularGrid,angularWeights,idx)
      case(230_int64)
        call buildLebedevRule0230(angularGrid,angularWeights,idx)
      case(266_int64)
        call buildLebedevRule0266(angularGrid,angularWeights,idx)
      case(302_int64)
        call buildLebedevRule0302(angularGrid,angularWeights,idx)
      case(350_int64)
        call buildLebedevRule0350(angularGrid,angularWeights,idx)
      case(434_int64)
        call buildLebedevRule0434(angularGrid,angularWeights,idx)
      case(590_int64)
        call buildLebedevRule0590(angularGrid,angularWeights,idx)
      case(770_int64)
        call buildLebedevRule0770(angularGrid,angularWeights,idx)
      case(974_int64)
        call buildLebedevRule0974(angularGrid,angularWeights,idx)
      case(1202_int64)
        call buildLebedevRule1202(angularGrid,angularWeights,idx)
      case(1454_int64)
        call buildLebedevRule1454(angularGrid,angularWeights,idx)
      case(1730_int64)
        call buildLebedevRule1730(angularGrid,angularWeights,idx)
      case(2030_int64)
        call buildLebedevRule2030(angularGrid,angularWeights,idx)
      end select
      if(idx.ne.lebedevPoints)  &
        call mqc_error('Lebedev rule generated an invalid point count.')
      angularWeights = mqc_float(4)*Pi*angularWeights
!
      return
      end subroutine setup_lebedev_angular_grid


!PROCEDURE lebedevRuleIsSupported
      function lebedevRuleIsSupported(lebedevPoints) result(isSupported)
!
!     This function reports whether a point count selects an available
!     Lebedev-Laikov angular rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::lebedevPoints
      logical::isSupported
!
      isSupported = ANY(LEBEDEV_POINT_COUNTS.eq.lebedevPoints)
!
      return
      end function lebedevRuleIsSupported


!PROCEDURE lebedevSelectPointCount
      function lebedevSelectPointCount(requestedPoints,roundUp,  &
        capAtMaximum) result(selectedPoints)
!
!     This function selects an available Lebedev angular point count. Exact
!     requests are returned unchanged. Otherwise, the default is to select the
!     next larger rule; roundUp=.false. selects the next smaller rule.
!     Requests above the implemented maximum are errors unless
!     capAtMaximum=.true., in which case the largest rule is returned.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::requestedPoints
      logical,optional,intent(in)::roundUp,capAtMaximum
      integer(kind=int64)::i,selectedPoints
      logical::capMaximum,selectUpward
!
      selectUpward = .true.
      if(Present(roundUp)) selectUpward = roundUp
      capMaximum = .false.
      if(Present(capAtMaximum)) capMaximum = capAtMaximum
!
      if(requestedPoints.lt.1_int64)  &
        call mqc_error('Lebedev angular point count must be positive.')
      if(requestedPoints.gt.LEBEDEV_POINT_COUNTS(N_LEBEDEV_RULES)) then
        if(capMaximum) then
          selectedPoints = LEBEDEV_POINT_COUNTS(N_LEBEDEV_RULES)
          return
        endIf
        call mqc_error('Requested Lebedev angular point count exceeds the implemented maximum.')
      endIf
!
      selectedPoints = -1_int64
      if(selectUpward) then
        do i = 1,N_LEBEDEV_RULES
          if(LEBEDEV_POINT_COUNTS(i).ge.requestedPoints) then
            selectedPoints = LEBEDEV_POINT_COUNTS(i)
            exit
          endIf
        endDo
      else
        do i = N_LEBEDEV_RULES,1,-1
          if(LEBEDEV_POINT_COUNTS(i).le.requestedPoints) then
            selectedPoints = LEBEDEV_POINT_COUNTS(i)
            exit
          endIf
        endDo
      endIf
      if(selectedPoints.lt.0_int64)  &
        call mqc_error('No supported Lebedev angular point count in the requested direction.')
!
      return
      end function lebedevSelectPointCount


!PROCEDURE lebedevRuleDegree
      function lebedevRuleDegree(lebedevPoints) result(exactnessDegree)
!
!     This function returns the polynomial exactness degree associated
!     with a supported Lebedev-Laikov point count.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::lebedevPoints
      integer(kind=int64)::exactnessDegree,i
!
      exactnessDegree = -1_int64
      do i = 1,N_LEBEDEV_RULES
        if(LEBEDEV_POINT_COUNTS(i).eq.lebedevPoints) then
          exactnessDegree = LEBEDEV_EXACTNESS_DEGREES(i)
          exit
        endIf
      endDo
      if(exactnessDegree.lt.0_int64)  &
        call mqc_error('Unsupported Lebedev angular point count.')
!
      return
      end function lebedevRuleDegree


!PROCEDURE buildLebedevRule0006
      subroutine buildLebedevRule0006(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 6-point, degree-3 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.1666666666666667e+0_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0006


!PROCEDURE buildLebedevRule0014
      subroutine buildLebedevRule0014(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 14-point, degree-5 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.6666666666666667e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.7500000000000000e-1_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0014


!PROCEDURE buildLebedevRule0026
      subroutine buildLebedevRule0026(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 26-point, degree-7 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.4761904761904762e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.3809523809523810e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.3214285714285714e-1_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0026


!PROCEDURE buildLebedevRule0038
      subroutine buildLebedevRule0038(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 38-point, degree-9 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.9523809523809524e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.3214285714285714e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.4597008433809831e+0_real64,0.0_real64,  &
        0.2857142857142857e-1_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0038


!PROCEDURE buildLebedevRule0050
      subroutine buildLebedevRule0050(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 50-point, degree-11 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.1269841269841270e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.2257495590828924e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.2109375000000000e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3015113445777636e+0_real64,0.0_real64,  &
        0.2017333553791887e-1_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0050


!PROCEDURE buildLebedevRule0074
      subroutine buildLebedevRule0074(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 74-point, degree-13 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.5130671797338464e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.1660406956574204e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        -0.2958603896103896e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4803844614152614e+0_real64,0.0_real64,  &
        0.2657620708215946e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3207726489807764e+0_real64,0.0_real64,  &
        0.1652217099371571e-1_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0074


!PROCEDURE buildLebedevRule0086
      subroutine buildLebedevRule0086(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 86-point, degree-15 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.1154401154401154e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.1194390908585628e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3696028464541502e+0_real64,0.0_real64,  &
        0.1111055571060340e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6943540066026664e+0_real64,0.0_real64,  &
        0.1187650129453714e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3742430390903412e+0_real64,0.0_real64,  &
        0.1181230374690448e-1_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0086


!PROCEDURE buildLebedevRule0110
      subroutine buildLebedevRule0110(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 110-point, degree-17 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.3828270494937162e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.9793737512487512e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1851156353447362e+0_real64,0.0_real64,  &
        0.8211737283191111e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6904210483822922e+0_real64,0.0_real64,  &
        0.9942814891178103e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3956894730559419e+0_real64,0.0_real64,  &
        0.9595471336070963e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.4783690288121502e+0_real64,0.0_real64,  &
        0.9694996361663028e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0110


!PROCEDURE buildLebedevRule0146
      subroutine buildLebedevRule0146(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 146-point, degree-19 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.5996313688621381e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.7372999718620756e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.7210515360144488e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6764410400114264e+0_real64,0.0_real64,  &
        0.7116355493117555e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4174961227965453e+0_real64,0.0_real64,  &
        0.6753829486314477e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1574676672039082e+0_real64,0.0_real64,  &
        0.7574394159054034e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1403553811713183e+0_real64,0.4493328323269557e+0_real64,  &
        0.6991087353303262e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0146


!PROCEDURE buildLebedevRule0170
      subroutine buildLebedevRule0170(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 170-point, degree-21 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.5544842902037365e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.6071332770670752e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.6383674773515093e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2551252621114134e+0_real64,0.0_real64,  &
        0.5183387587747790e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6743601460362766e+0_real64,0.0_real64,  &
        0.6317929009813725e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4318910696719410e+0_real64,0.0_real64,  &
        0.6201670006589077e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.2613931360335988e+0_real64,0.0_real64,  &
        0.5477143385137348e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4990453161796037e+0_real64,0.1446630744325115e+0_real64,  &
        0.5968383987681156e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0170


!PROCEDURE buildLebedevRule0194
      subroutine buildLebedevRule0194(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 194-point, degree-23 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.1782340447244611e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.5716905949977102e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.5573383178848738e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6712973442695226e+0_real64,0.0_real64,  &
        0.5608704082587997e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2892465627575439e+0_real64,0.0_real64,  &
        0.5158237711805383e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4446933178717437e+0_real64,0.0_real64,  &
        0.5518771467273614e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1299335447650067e+0_real64,0.0_real64,  &
        0.4106777028169394e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3457702197611283e+0_real64,0.0_real64,  &
        0.5051846064614808e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1590417105383530e+0_real64,0.8360360154824589e+0_real64,  &
        0.5530248916233094e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0194


!PROCEDURE buildLebedevRule0230
      subroutine buildLebedevRule0230(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 230-point, degree-25 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        -0.5522639919727325e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.4450274607445226e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4492044687397611e+0_real64,0.0_real64,  &
        0.4496841067921404e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2520419490210201e+0_real64,0.0_real64,  &
        0.5049153450478750e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6981906658447242e+0_real64,0.0_real64,  &
        0.3976408018051883e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6587405243460960e+0_real64,0.0_real64,  &
        0.4401400650381014e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4038544050097660e-1_real64,0.0_real64,  &
        0.1724544350544401e-1_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.5823842309715585e+0_real64,0.0_real64,  &
        0.4231083095357343e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3545877390518688e+0_real64,0.0_real64,  &
        0.5198069864064399e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2272181808998187e+0_real64,0.4864661535886647e+0_real64,  &
        0.4695720972568883e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0230


!PROCEDURE buildLebedevRule0266
      subroutine buildLebedevRule0266(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 266-point, degree-27 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        -0.1313769127326952e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        -0.2522728704859336e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.4186853881700583e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7039373391585475e+0_real64,0.0_real64,  &
        0.5315167977810885e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1012526248572414e+0_real64,0.0_real64,  &
        0.4047142377086219e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4647448726420539e+0_real64,0.0_real64,  &
        0.4112482394406990e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3277420654971629e+0_real64,0.0_real64,  &
        0.3595584899758782e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6620338663699974e+0_real64,0.0_real64,  &
        0.4256131351428158e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.8506508083520399e+0_real64,0.0_real64,  &
        0.4229582700647240e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3233484542692899e+0_real64,0.1153112011009701e+0_real64,  &
        0.4080914225780505e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2314790158712601e+0_real64,0.5244939240922365e+0_real64,  &
        0.4071467593830964e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0266


!PROCEDURE buildLebedevRule0302
      subroutine buildLebedevRule0302(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 302-point, degree-29 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.8545911725128148e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.3599119285025571e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3515640345570105e+0_real64,0.0_real64,  &
        0.3449788424305883e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6566329410219612e+0_real64,0.0_real64,  &
        0.3604822601419882e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4729054132581005e+0_real64,0.0_real64,  &
        0.3576729661743367e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.9618308522614784e-1_real64,0.0_real64,  &
        0.2352101413689164e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2219645236294178e+0_real64,0.0_real64,  &
        0.3108953122413675e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7011766416089545e+0_real64,0.0_real64,  &
        0.3650045807677255e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.2644152887060663e+0_real64,0.0_real64,  &
        0.2982344963171804e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.5718955891878961e+0_real64,0.0_real64,  &
        0.3600820932216460e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2510034751770465e+0_real64,0.8000727494073952e+0_real64,  &
        0.3571540554273387e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1233548532583327e+0_real64,0.4127724083168531e+0_real64,  &
        0.3392312205006170e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0302


!PROCEDURE buildLebedevRule0350
      subroutine buildLebedevRule0350(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 350-point, degree-31 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.3006796749453936e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.3050627745650771e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7068965463912316e+0_real64,0.0_real64,  &
        0.1621104600288991e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4794682625712025e+0_real64,0.0_real64,  &
        0.3005701484901752e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1927533154878019e+0_real64,0.0_real64,  &
        0.2990992529653774e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6930357961327123e+0_real64,0.0_real64,  &
        0.2982170644107595e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3608302115520091e+0_real64,0.0_real64,  &
        0.2721564237310992e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6498486161496169e+0_real64,0.0_real64,  &
        0.3033513795811141e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.1932945013230339e+0_real64,0.0_real64,  &
        0.3007949555218533e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3800494919899303e+0_real64,0.0_real64,  &
        0.2881964603055307e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2899558825499574e+0_real64,0.7934537856582316e+0_real64,  &
        0.2958357626535696e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.9684121455103957e-1_real64,0.8280801506686862e+0_real64,  &
        0.3036020026407088e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1833434647041659e+0_real64,0.9074658265305127e+0_real64,  &
        0.2832187403926303e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0350


!PROCEDURE buildLebedevRule0434
      subroutine buildLebedevRule0434(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 434-point, degree-35 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.5265897968224436e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.2548219972002607e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.2512317418927307e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6909346307509111e+0_real64,0.0_real64,  &
        0.2530403801186355e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1774836054609158e+0_real64,0.0_real64,  &
        0.2014279020918528e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4914342637784746e+0_real64,0.0_real64,  &
        0.2501725168402936e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6456664707424256e+0_real64,0.0_real64,  &
        0.2513267174597564e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2861289010307638e+0_real64,0.0_real64,  &
        0.2302694782227416e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7568084367178018e-1_real64,0.0_real64,  &
        0.1462495621594614e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3927259763368002e+0_real64,0.0_real64,  &
        0.2445373437312980e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.8818132877794288e+0_real64,0.0_real64,  &
        0.2417442375638981e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.9776428111182649e+0_real64,0.0_real64,  &
        0.1910951282179532e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2054823696403044e+0_real64,0.8689460322872412e+0_real64,  &
        0.2416930044324775e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5905157048925271e+0_real64,0.7999278543857286e+0_real64,  &
        0.2512236854563495e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5550152361076807e+0_real64,0.7717462626915901e+0_real64,  &
        0.2496644054553086e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.9371809858553722e+0_real64,0.3344363145343455e+0_real64,  &
        0.2236607760437849e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0434


!PROCEDURE buildLebedevRule0590
      subroutine buildLebedevRule0590(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 590-point, degree-41 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.3095121295306187e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.1852379698597489e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7040954938227469e+0_real64,0.0_real64,  &
        0.1871790639277744e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6807744066455243e+0_real64,0.0_real64,  &
        0.1858812585438317e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6372546939258752e+0_real64,0.0_real64,  &
        0.1852028828296213e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.5044419707800358e+0_real64,0.0_real64,  &
        0.1846715956151242e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4215761784010967e+0_real64,0.0_real64,  &
        0.1818471778162769e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3317920736472123e+0_real64,0.0_real64,  &
        0.1749564657281154e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2384736701421887e+0_real64,0.0_real64,  &
        0.1617210647254411e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1459036449157763e+0_real64,0.0_real64,  &
        0.1384737234851692e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6095034115507196e-1_real64,0.0_real64,  &
        0.9764331165051050e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.6116843442009876e+0_real64,0.0_real64,  &
        0.1857161196774078e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3964755348199858e+0_real64,0.0_real64,  &
        0.1705153996395864e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.1724782009907724e+0_real64,0.0_real64,  &
        0.1300321685886048e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5610263808622060e+0_real64,0.3518280927733519e+0_real64,  &
        0.1842866472905286e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4742392842551980e+0_real64,0.2634716655937950e+0_real64,  &
        0.1802658934377451e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5984126497885380e+0_real64,0.1816640840360209e+0_real64,  &
        0.1849830560443660e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3791035407695563e+0_real64,0.1720795225656878e+0_real64,  &
        0.1713904507106709e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2778673190586244e+0_real64,0.8213021581932511e-1_real64,  &
        0.1555213603396808e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5033564271075117e+0_real64,0.8999205842074875e-1_real64,  &
        0.1802239128008525e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0590


!PROCEDURE buildLebedevRule0770
      subroutine buildLebedevRule0770(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 770-point, degree-47 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.2192942088181184e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.1436433617319080e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.1421940344335877e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.5087204410502360e-1_real64,0.0_real64,  &
        0.6798123511050502e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1228198790178831e+0_real64,0.0_real64,  &
        0.9913184235294912e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2026890814408786e+0_real64,0.0_real64,  &
        0.1180207833238949e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2847745156464294e+0_real64,0.0_real64,  &
        0.1296599602080921e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3656719078978026e+0_real64,0.0_real64,  &
        0.1365871427428316e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4428264886713469e+0_real64,0.0_real64,  &
        0.1402988604775325e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.5140619627249735e+0_real64,0.0_real64,  &
        0.1418645563595609e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6306401219166803e+0_real64,0.0_real64,  &
        0.1421376741851662e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6716883332022612e+0_real64,0.0_real64,  &
        0.1423996475490962e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6979792685336881e+0_real64,0.0_real64,  &
        0.1431554042178567e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.1446865674195309e+0_real64,0.0_real64,  &
        0.9254401499865368e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3390263475411216e+0_real64,0.0_real64,  &
        0.1250239995053509e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.5335804651263506e+0_real64,0.0_real64,  &
        0.1394365843329230e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6944024393349413e-1_real64,0.2355187894242326e+0_real64,  &
        0.1127089094671749e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2269004109529460e+0_real64,0.4102182474045730e+0_real64,  &
        0.1345753760910670e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.8025574607775339e-1_real64,0.6214302417481605e+0_real64,  &
        0.1424957283316783e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1467999527896572e+0_real64,0.3245284345717394e+0_real64,  &
        0.1261523341237750e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1571507769824727e+0_real64,0.5224482189696630e+0_real64,  &
        0.1392547106052696e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2365702993157246e+0_real64,0.6017546634089558e+0_real64,  &
        0.1418761677877656e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.7714815866765732e-1_real64,0.4346575516141163e+0_real64,  &
        0.1338366684479554e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3062936666210730e+0_real64,0.4908826589037616e+0_real64,  &
        0.1393700862676131e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3822477379524787e+0_real64,0.5648768149099500e+0_real64,  &
        0.1415914757466932e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0770


!PROCEDURE buildLebedevRule0974
      subroutine buildLebedevRule0974(angularGrid,  &
        angularWeights,idx)
!
!     This routine appends the 974-point, degree-53 rule.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
!
      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.1438294190527431e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.1125772288287004e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4292963545341347e-1_real64,0.0_real64,  &
        0.4948029341949241e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1051426854086404e+0_real64,0.0_real64,  &
        0.7357990109125470e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1750024867623087e+0_real64,0.0_real64,  &
        0.8889132771304384e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2477653379650257e+0_real64,0.0_real64,  &
        0.9888347838921435e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3206567123955957e+0_real64,0.0_real64,  &
        0.1053299681709471e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3916520749849983e+0_real64,0.0_real64,  &
        0.1092778807014578e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4590825874187624e+0_real64,0.0_real64,  &
        0.1114389394063227e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.5214563888415861e+0_real64,0.0_real64,  &
        0.1123724788051555e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6253170244654199e+0_real64,0.0_real64,  &
        0.1125239325243814e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6637926744523170e+0_real64,0.0_real64,  &
        0.1126153271815905e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6910410398498301e+0_real64,0.0_real64,  &
        0.1130286931123841e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7052907007457760e+0_real64,0.0_real64,  &
        0.1134986534363955e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.1236686762657990e+0_real64,0.0_real64,  &
        0.6823367927109931e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.2940777114468387e+0_real64,0.0_real64,  &
        0.9454158160447096e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.4697753849207649e+0_real64,0.0_real64,  &
        0.1074429975385679e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.6334563241139567e+0_real64,0.0_real64,  &
        0.1129300086569132e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5974048614181342e-1_real64,0.2029128752777523e+0_real64,  &
        0.8436884500901954e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1375760408473636e+0_real64,0.4602621942484054e+0_real64,  &
        0.1075255720448885e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3391016526336286e+0_real64,0.5030673999662036e+0_real64,  &
        0.1108577236864462e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1271675191439820e+0_real64,0.2817606422442134e+0_real64,  &
        0.9566475323783357e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2693120740413512e+0_real64,0.4331561291720157e+0_real64,  &
        0.1080663250717391e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1419786452601918e+0_real64,0.6256167358580814e+0_real64,  &
        0.1126797131196295e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6709284600738255e-1_real64,0.3798395216859157e+0_real64,  &
        0.1022568715358061e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.7057738183256172e-1_real64,0.5517505421423520e+0_real64,  &
        0.1108960267713108e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2783888477882155e+0_real64,0.6029619156159187e+0_real64,  &
        0.1122790653435766e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1979578938917407e+0_real64,0.3589606329589096e+0_real64,  &
        0.1032401847117460e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2087307061103274e+0_real64,0.5348666438135476e+0_real64,  &
        0.1107249382283854e-2_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4055122137872836e+0_real64,0.5674997546074373e+0_real64,  &
        0.1121780048519972e-2_real64,angularGrid,angularWeights,idx)
!
      return
      end subroutine buildLebedevRule0974


!PROCEDURE buildLebedevRule1202
      subroutine buildLebedevRule1202(angularGrid,angularWeights,idx)
!
!     CCL Lebedev-Laikov orbit coefficients; H. P. Hratchian, 2026.
!
      implicit none
      real(real64),intent(inout)::angularGrid(:,:),angularWeights(:)
      integer(int64),intent(inout)::idx

      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.1105189233267572e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.9205232738090741e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.9133159786443561e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3712636449657089e-1_real64,0.0_real64,  &
        0.3690421898017899e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.9140060412262223e-1_real64,0.0_real64,  &
        0.5603990928680660e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1531077852469906e+0_real64,0.0_real64,  &
        0.6865297629282609e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2180928891660612e+0_real64,0.0_real64,  &
        0.7720338551145630e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2839874532200175e+0_real64,0.0_real64,  &
        0.8301545958894795e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3491177600963764e+0_real64,0.0_real64,  &
        0.8686692550179628e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4121431461444309e+0_real64,0.0_real64,  &
        0.8927076285846890e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4718993627149127e+0_real64,0.0_real64,  &
        0.9060820238568219e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.5273145452842337e+0_real64,0.0_real64,  &
        0.9119777254940867e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6209475332444019e+0_real64,0.0_real64,  &
        0.9128720138604181e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6569722711857291e+0_real64,0.0_real64,  &
        0.9130714935691735e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6841788309070143e+0_real64,0.0_real64,  &
        0.9152873784554116e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7012604330123631e+0_real64,0.0_real64,  &
        0.9187436274321654e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.1072382215478166e+0_real64,0.0_real64,  &
        0.5176977312965694e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.2582068959496968e+0_real64,0.0_real64,  &
        0.7331143682101417e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.4172752955306717e+0_real64,0.0_real64,  &
        0.8463232836379928e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.5700366911792503e+0_real64,0.0_real64,  &
        0.9031122694253992e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.9827986018263947e+0_real64,0.1771774022615325e+0_real64,  &
        0.6485778453163257e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.9624249230326228e+0_real64,0.2475716463426288e+0_real64,  &
        0.7435030910982369e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.9402007994128811e+0_real64,0.3354616289066489e+0_real64,  &
        0.7998527891839054e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.9320822040143202e+0_real64,0.3173615246611977e+0_real64,  &
        0.8101731497468018e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.9043674199393299e+0_real64,0.4090268427085357e+0_real64,  &
        0.8483389574594331e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.8912407560074747e+0_real64,0.3854291150669224e+0_real64,  &
        0.8556299257311812e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.8676435628462708e+0_real64,0.4932221184851285e+0_real64,  &
        0.8803208679738260e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.8581979986041619e+0_real64,0.4785320675922435e+0_real64,  &
        0.8811048182425720e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.8396753624049856e+0_real64,0.4507422593157064e+0_real64,  &
        0.8850282341265444e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.8165288564022188e+0_real64,0.5632123020762100e+0_real64,  &
        0.9021342299040653e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.8015469370783529e+0_real64,0.5434303569693900e+0_real64,  &
        0.9010091677105086e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.7773563069070351e+0_real64,0.5123518486419871e+0_real64,  &
        0.9022692938426915e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.7661621213900394e+0_real64,0.6394279634749102e+0_real64,  &
        0.9158016174693465e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.7553584143533510e+0_real64,0.6269805509024392e+0_real64,  &
        0.9131578003189435e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.7344305757559503e+0_real64,0.6031161693096310e+0_real64,  &
        0.9107813579482705e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.7043837184021765e+0_real64,0.5693702498468441e+0_real64,  &
        0.9105760258970126e-3_real64,angularGrid,angularWeights,idx)
      end subroutine buildLebedevRule1202


!PROCEDURE buildLebedevRule1454
      subroutine buildLebedevRule1454(angularGrid,angularWeights,idx)
!
!     CCL Lebedev-Laikov orbit coefficients; H. P. Hratchian, 2026.
!
      implicit none
      real(real64),intent(inout)::angularGrid(:,:),angularWeights(:)
      integer(int64),intent(inout)::idx

      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.7777160743261247e-4_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.7557646413004701e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3229290663413854e-1_real64,0.0_real64,  &
        0.2841633806090617e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.8036733271462222e-1_real64,0.0_real64,  &
        0.4374419127053555e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1354289960531653e+0_real64,0.0_real64,  &
        0.5417174740872172e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1938963861114426e+0_real64,0.0_real64,  &
        0.6148000891358593e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2537343715011275e+0_real64,0.0_real64,  &
        0.6664394485800705e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3135251434752570e+0_real64,0.0_real64,  &
        0.7025039356923220e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3721558339375338e+0_real64,0.0_real64,  &
        0.7268511789249627e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4286809575195696e+0_real64,0.0_real64,  &
        0.7422637534208629e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4822510128282994e+0_real64,0.0_real64,  &
        0.7509545035841214e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.5320679333566263e+0_real64,0.0_real64,  &
        0.7548535057718401e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6172998195394274e+0_real64,0.0_real64,  &
        0.7554088969774001e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6510679849127481e+0_real64,0.0_real64,  &
        0.7553147174442808e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6777315251687360e+0_real64,0.0_real64,  &
        0.7564767653292297e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6963109410648741e+0_real64,0.0_real64,  &
        0.7587991808518730e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7058935009831749e+0_real64,0.0_real64,  &
        0.7608261832033027e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.9955546194091857e+0_real64,0.0_real64,  &
        0.4021680447874916e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.9734115901794209e+0_real64,0.0_real64,  &
        0.5804871793945964e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.9275693732388626e+0_real64,0.0_real64,  &
        0.6792151955945159e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.8568022422795103e+0_real64,0.0_real64,  &
        0.7336741211286294e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.7623495553719372e+0_real64,0.0_real64,  &
        0.7581866300989608e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5707522908892223e+0_real64,0.4387028039889501e+0_real64,  &
        0.7538257859800743e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5196463388403083e+0_real64,0.3858908414762617e+0_real64,  &
        0.7483517247053123e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4646337531215351e+0_real64,0.3301937372343854e+0_real64,  &
        0.7371763661112059e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4063901697557691e+0_real64,0.2725423573563777e+0_real64,  &
        0.7183448895756934e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3456329466643087e+0_real64,0.2139510237495250e+0_real64,  &
        0.6895815529822191e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2831395121050332e+0_real64,0.1555922309786647e+0_real64,  &
        0.6480105801792886e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2197682022925330e+0_real64,0.9892878979686097e-1_real64,  &
        0.5897558896594636e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1564696098650355e+0_real64,0.4598642910675510e-1_real64,  &
        0.5095708849247346e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6027356673721295e+0_real64,0.3376625140173426e+0_real64,  &
        0.7536906428909755e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5496032320255096e+0_real64,0.2822301309727988e+0_real64,  &
        0.7472505965575118e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4921707755234567e+0_real64,0.2248632342592540e+0_real64,  &
        0.7343017132279698e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4309422998598483e+0_real64,0.1666224723456479e+0_real64,  &
        0.7130871582177445e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3664108182313672e+0_real64,0.1086964901822169e+0_real64,  &
        0.6817022032112776e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2990189057758436e+0_real64,0.5251989784120085e-1_real64,  &
        0.6380941145604121e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6268724013144998e+0_real64,0.2297523657550023e+0_real64,  &
        0.7550381377920310e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5707324144834607e+0_real64,0.1723080607093800e+0_real64,  &
        0.7478646640144802e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5096360901960365e+0_real64,0.1140238465390513e+0_real64,  &
        0.7335918720601220e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4438729938312456e+0_real64,0.5611522095882537e-1_real64,  &
        0.7110120527658118e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6419978471082389e+0_real64,0.1164174423140873e+0_real64,  &
        0.7571363978689501e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5817218061802611e+0_real64,0.5797589531445219e-1_real64,  &
        0.7489908329079234e-3_real64,angularGrid,angularWeights,idx)
      end subroutine buildLebedevRule1454


!PROCEDURE buildLebedevRule1730
      subroutine buildLebedevRule1730(angularGrid,angularWeights,idx)
!
!     CCL Lebedev-Laikov orbit coefficients; H. P. Hratchian, 2026.
!
      implicit none
      real(real64),intent(inout)::angularGrid(:,:),angularWeights(:)
      integer(int64),intent(inout)::idx

      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.6309049437420976e-4_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(2_int64,0.0_real64,0.0_real64,  &
        0.6398287705571748e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.6357185073530720e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2860923126194662e-1_real64,0.0_real64,  &
        0.2221207162188168e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7142556767711522e-1_real64,0.0_real64,  &
        0.3475784022286848e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1209199540995559e+0_real64,0.0_real64,  &
        0.4350742443589804e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1738673106594379e+0_real64,0.0_real64,  &
        0.4978569136522127e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2284645438467734e+0_real64,0.0_real64,  &
        0.5435036221998053e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2834807671701512e+0_real64,0.0_real64,  &
        0.5765913388219542e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3379680145467339e+0_real64,0.0_real64,  &
        0.6001200359226003e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3911355454819537e+0_real64,0.0_real64,  &
        0.6162178172717512e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4422860353001403e+0_real64,0.0_real64,  &
        0.6265218152438485e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4907781568726057e+0_real64,0.0_real64,  &
        0.6323987160974212e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.5360006153211468e+0_real64,0.0_real64,  &
        0.6350767851540569e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6142105973596603e+0_real64,0.0_real64,  &
        0.6354362775297107e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6459300387977504e+0_real64,0.0_real64,  &
        0.6352302462706235e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6718056125089225e+0_real64,0.0_real64,  &
        0.6358117881417972e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6910888533186254e+0_real64,0.0_real64,  &
        0.6373101590310117e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7030467416823252e+0_real64,0.0_real64,  &
        0.6390428961368665e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.8354951166354646e-1_real64,0.0_real64,  &
        0.3186913449946576e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.2050143009099486e+0_real64,0.0_real64,  &
        0.4678028558591711e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3370208290706637e+0_real64,0.0_real64,  &
        0.5538829697598626e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.4689051484233963e+0_real64,0.0_real64,  &
        0.6044475907190476e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.5939400424557334e+0_real64,0.0_real64,  &
        0.6313575103509012e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1394983311832261e+0_real64,0.4097581162050343e-1_real64,  &
        0.4078626431855630e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1967999180485014e+0_real64,0.8851987391293348e-1_real64,  &
        0.4759933057812725e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2546183732548967e+0_real64,0.1397680182969819e+0_real64,  &
        0.5268151186413440e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3121281074713875e+0_real64,0.1929452542226526e+0_real64,  &
        0.5643048560507316e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3685981078502492e+0_real64,0.2467898337061562e+0_real64,  &
        0.5914501076613073e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4233760321547856e+0_real64,0.3003104124785409e+0_real64,  &
        0.6104561257874195e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4758671236059246e+0_real64,0.3526684328175033e+0_real64,  &
        0.6230252860707806e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5255178579796463e+0_real64,0.4031134861145713e+0_real64,  &
        0.6305618761760796e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5718025633734589e+0_real64,0.4509426448342351e+0_real64,  &
        0.6343092767597889e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2686927772723415e+0_real64,0.4711322502423248e-1_real64,  &
        0.5176268945737826e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3306006819904809e+0_real64,0.9784487303942695e-1_real64,  &
        0.5564840313313692e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3904906850594983e+0_real64,0.1505395810025273e+0_real64,  &
        0.5856426671038980e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4479957951904390e+0_real64,0.2039728156296050e+0_real64,  &
        0.6066386925777091e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5027076848919780e+0_real64,0.2571529941121107e+0_real64,  &
        0.6208824962234458e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5542087392260217e+0_real64,0.3092191375815670e+0_real64,  &
        0.6296314297822907e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6020850887375187e+0_real64,0.3593807506130276e+0_real64,  &
        0.6340423756791859e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4019851409179594e+0_real64,0.5063389934378671e-1_real64,  &
        0.5829627677107342e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4635614567449800e+0_real64,0.1032422269160612e+0_real64,  &
        0.6048693376081110e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5215860931591575e+0_real64,0.1566322094006254e+0_real64,  &
        0.6202362317732461e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5758202499099271e+0_real64,0.2098082827491099e+0_real64,  &
        0.6299005328403779e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6259893683876795e+0_real64,0.2618824114553391e+0_real64,  &
        0.6347722390609353e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5313795124811891e+0_real64,0.5263245019338556e-1_real64,  &
        0.6203778981238834e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5893317955931995e+0_real64,0.1061059730982005e+0_real64,  &
        0.6308414671239979e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6426246321215801e+0_real64,0.1594171564034221e+0_real64,  &
        0.6362706466959498e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6511904367376113e+0_real64,0.5354789536565540e-1_real64,  &
        0.6375414170333233e-3_real64,angularGrid,angularWeights,idx)
      end subroutine buildLebedevRule1730


!PROCEDURE buildLebedevRule2030
      subroutine buildLebedevRule2030(angularGrid,angularWeights,idx)
!
!     CCL Lebedev-Laikov orbit coefficients; H. P. Hratchian, 2026.
!
      implicit none
      real(real64),intent(inout)::angularGrid(:,:),angularWeights(:)
      integer(int64),intent(inout)::idx

      call addLebedevOrbit(1_int64,0.0_real64,0.0_real64,  &
        0.4656031899197431e-4_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(3_int64,0.0_real64,0.0_real64,  &
        0.5421549195295507e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2540835336814348e-1_real64,0.0_real64,  &
        0.1778522133346553e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6399322800504915e-1_real64,0.0_real64,  &
        0.2811325405682796e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1088269469804125e+0_real64,0.0_real64,  &
        0.3548896312631459e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.1570670798818287e+0_real64,0.0_real64,  &
        0.4090310897173364e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2071163932282514e+0_real64,0.0_real64,  &
        0.4493286134169965e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.2578914044450844e+0_real64,0.0_real64,  &
        0.4793728447962723e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3085687558169623e+0_real64,0.0_real64,  &
        0.5015415319164265e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.3584719706267024e+0_real64,0.0_real64,  &
        0.5175127372677937e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4070135594428709e+0_real64,0.0_real64,  &
        0.5285522262081019e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4536618626222638e+0_real64,0.0_real64,  &
        0.5356832703713962e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.4979195686463577e+0_real64,0.0_real64,  &
        0.5397914736175170e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.5393075111126999e+0_real64,0.0_real64,  &
        0.5416899441599930e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6115617676843916e+0_real64,0.0_real64,  &
        0.5419308476889938e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6414308435160159e+0_real64,0.0_real64,  &
        0.5416936902030596e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6664099412721607e+0_real64,0.0_real64,  &
        0.5419544338703164e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6859161771214913e+0_real64,0.0_real64,  &
        0.5428983656630975e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.6993625593503890e+0_real64,0.0_real64,  &
        0.5442286500098193e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(4_int64,0.7062393387719380e+0_real64,0.0_real64,  &
        0.5452250345057301e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.7479028168349763e-1_real64,0.0_real64,  &
        0.2568002497728530e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.1848951153969366e+0_real64,0.0_real64,  &
        0.3827211700292145e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.3059529066581305e+0_real64,0.0_real64,  &
        0.4579491561917824e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.4285556101021362e+0_real64,0.0_real64,  &
        0.5042003969083574e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.5468758653496526e+0_real64,0.0_real64,  &
        0.5312708889976025e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(5_int64,0.6565821978343439e+0_real64,0.0_real64,  &
        0.5438401790747117e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1253901572367117e+0_real64,0.3681917226439641e-1_real64,  &
        0.3316041873197344e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.1775721510383941e+0_real64,0.7982487607213301e-1_real64,  &
        0.3899113567153771e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2305693358216114e+0_real64,0.1264640966592335e+0_real64,  &
        0.4343343327201309e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2836502845992063e+0_real64,0.1751585683418957e+0_real64,  &
        0.4679415262318919e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3361794746232590e+0_real64,0.2247995907632670e+0_real64,  &
        0.4930847981631031e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3875979172264824e+0_real64,0.2745299257422246e+0_real64,  &
        0.5115031867540091e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4374019316999074e+0_real64,0.3236373482441118e+0_real64,  &
        0.5245217148457367e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4851275843340022e+0_real64,0.3714967859436741e+0_real64,  &
        0.5332041499895321e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5303391803806868e+0_real64,0.4175353646321745e+0_real64,  &
        0.5384583126021542e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5726197380596287e+0_real64,0.4612084406355461e+0_real64,  &
        0.5411067210798852e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.2431520732564863e+0_real64,0.4258040133043952e-1_real64,  &
        0.4259797391468714e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3002096800895869e+0_real64,0.8869424306722721e-1_real64,  &
        0.4604931368460021e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3558554457457432e+0_real64,0.1368811706510655e+0_real64,  &
        0.4871814878255202e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4097782537048887e+0_real64,0.1860739985015033e+0_real64,  &
        0.5072242910074885e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4616337666067458e+0_real64,0.2354235077395853e+0_real64,  &
        0.5217069845235350e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5110707008417874e+0_real64,0.2842074921347011e+0_real64,  &
        0.5315785966280310e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5577415286163795e+0_real64,0.3317784414984102e+0_real64,  &
        0.5376833708758905e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6013060431366950e+0_real64,0.3775299002040700e+0_real64,  &
        0.5408032092069521e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.3661596767261781e+0_real64,0.4599367887164592e-1_real64,  &
        0.4842744917904866e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4237633153506581e+0_real64,0.9404893773654421e-1_real64,  &
        0.5048926076188130e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4786328454658452e+0_real64,0.1431377109091971e+0_real64,  &
        0.5202607980478373e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5305702076789774e+0_real64,0.1924186388843570e+0_real64,  &
        0.5309932388325743e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5793436224231788e+0_real64,0.2411590944775190e+0_real64,  &
        0.5377419770895208e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6247069017094747e+0_real64,0.2886871491583605e+0_real64,  &
        0.5411696331677717e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.4874315552535204e+0_real64,0.4804978774953206e-1_real64,  &
        0.5197996293282420e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5427337322059053e+0_real64,0.9716857199366665e-1_real64,  &
        0.5311120836622945e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.5943493747246700e+0_real64,0.1465205839795055e+0_real64,  &
        0.5384309319956951e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6421314033564943e+0_real64,0.1953579449803574e+0_real64,  &
        0.5421859504051886e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6020628374713980e+0_real64,0.4916375015738108e-1_real64,  &
        0.5390948355046314e-3_real64,angularGrid,angularWeights,idx)
      call addLebedevOrbit(6_int64,0.6529222529856881e+0_real64,0.9861621540127005e-1_real64,  &
        0.5433312705027845e-3_real64,angularGrid,angularWeights,idx)
      end subroutine buildLebedevRule2030


!PROCEDURE addLebedevOrbit
      subroutine addLebedevOrbit(orbitType,a,b,weight,angularGrid,  &
        angularWeights,idx)
!
!     This routine expands one representative point under octahedral
!     coordinate permutations and independent sign changes. Exact
!     comparisons remove duplicates created by zeros or repeated
!     coordinates from the same representative values.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::orbitType
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),intent(in)::a,b,weight
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
      integer(kind=int64)::iPerm,iPoint,iSignX,iSignY,iSignZ,nOrbit
      integer(kind=int64),dimension(3,6),parameter::permutations=  &
        RESHAPE([ 1_int64,2_int64,3_int64,1_int64,3_int64,2_int64,  &
          2_int64,1_int64,3_int64,2_int64,3_int64,1_int64,  &
          3_int64,1_int64,2_int64,3_int64,2_int64,1_int64 ],[ 3,6 ])
      logical::isNew
      real(kind=real64)::c
      real(kind=real64),dimension(3)::base,candidate
      real(kind=real64),dimension(3,48)::orbitPoints
!
      select case(orbitType)
      case(1_int64)
        base = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      case(2_int64)
        c = sqrt(mqc_float(1)/mqc_float(2))
        base = [ mqc_float(0),c,c ]
      case(3_int64)
        c = sqrt(mqc_float(1)/mqc_float(3))
        base = [ c,c,c ]
      case(4_int64)
        c = sqrt(mqc_float(1)-mqc_float(2)*a*a)
        base = [ a,a,c ]
      case(5_int64)
        c = sqrt(mqc_float(1)-a*a)
        base = [ a,c,mqc_float(0) ]
      case(6_int64)
        c = sqrt(mqc_float(1)-a*a-b*b)
        base = [ a,b,c ]
      case default
        call mqc_error('Invalid Lebedev octahedral orbit type.')
      end select
!
      nOrbit = 0_int64
      do iPerm = 1,6
        do iSignX = -1,1,2
          do iSignY = -1,1,2
            do iSignZ = -1,1,2
              candidate = base(permutations(:,iPerm))*  &
                [ mqc_float(iSignX),mqc_float(iSignY),  &
                  mqc_float(iSignZ) ]
              isNew = .true.
              do iPoint = 1,nOrbit
                if(ALL(candidate.eq.orbitPoints(:,iPoint))) then
                  isNew = .false.
                  exit
                endIf
              endDo
              if(isNew) then
                nOrbit = nOrbit+1_int64
                orbitPoints(:,nOrbit) = candidate
              endIf
            endDo
          endDo
        endDo
      endDo
!
      if(idx+nOrbit.gt.Size(angularWeights))  &
        call mqc_error('Lebedev orbit exceeds allocated angular grid.')
      do iPoint = 1,nOrbit
        idx = idx+1_int64
        angularGrid(:,idx) = orbitPoints(:,iPoint)
        angularWeights(idx) = weight
      endDo
!
      return
      end subroutine addLebedevOrbit
