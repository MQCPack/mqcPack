      program unitTest12
!
!     Object ownership, integration, metadata, procedural compatibility,
!     molecular provenance, and invalid-state contracts. H. P. Hratchian, 2026.
!
      use iso_fortran_env, only: int64,real64
      use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
      use mqc_general, only: mqc_error,mqc_version_print,Pi
      use mqc_integrals
      implicit none
      type(quadrature_grid)::grid,copy,radial,angular
      type(molecular_quadrature_grid)::molecule,moleculeCopy
      class(quadrature_grid),allocatable::polymorphic
      real(real64),allocatable::p(:,:),w(:),rp(:,:),rw(:),ap(:,:),aw(:),overlap(:,:),values(:,:)
      real(real64),allocatable::referencePoints(:,:),referenceWeights(:),centers(:,:),radii(:)
      real(real64),allocatable::base(:),partition(:),savedRadial(:),savedRW(:),savedAngular(:,:),savedAW(:)
      integer(int64),allocatable::atoms(:)
      integer(int64)::cd,id,degree,n,i,invalidCase
      real(real64)::x,before,minimum,maximum
      character(:),allocatable::units,weightUnits,rule,normalization
      character(32)::argument
      call mqc_version_print(6_int64)
      if(command_argument_count().gt.0) then
        call get_command_argument(1,argument)
        read(argument,*) invalidCase
        call invalidInput(invalidCase)
        stop
      endIf
      if(grid%isInitialized().or.grid%pointCount().ne.0) call mqc_error('unitTest12: initial state.')
      call grid%clear()
      call gridGaussLegendre1D(grid,4_int64,-1.0_real64,1.0_real64,'bohr')
      allocate(p(1,4),w(4))
      call grid%copyBlock(1_int64,p,w)
      do i = 0,7
        x = 0.0_real64
        if(mod(i,2_int64).eq.0) x = 2.0_real64/real(i+1,real64)
        call near(grid%integrate(p(1,:)**i),x,2.0e-13_real64,'Gauss-Legendre moment')
      endDo
      call grid%getMetadata(cd,id,units,weightUnits,rule,normalization,degree)
      if(cd.ne.1.or.id.ne.1.or.degree.ne.7.or.units.ne.'bohr') call mqc_error('unitTest12: metadata.')
      copy = grid
      call grid%clear()
      call near(copy%integrate(w*0+1),2.0_real64,2.0e-13_real64,'deep copy survives clear')
      call gridTakeArrays(grid,p,w,1_int64,'bohr','bohr','imported','line measure')
      if(allocated(p).or.allocated(w)) call mqc_error('unitTest12: ownership transfer.')
      call grid%validate()
      call gridTrapezoid1D(grid,3_int64,0.5_real64,0.0_real64,'bohr')
      call near(grid%integrate([0.0_real64,0.5_real64,1.0_real64]),0.5_real64,1.0e-14_real64,'trapezoid')
      call gridCartesian(grid,3_int64,0.5_real64,[0.0_real64,0.0_real64,0.0_real64],'bohr')
      allocate(p(3,27),w(27),referencePoints(3,27),referenceWeights(27))
      call grid%copyBlock(1_int64,p,w)
      call setup_quadrature_trapezoid3d(3_int64,0.5_real64,[0.0_real64,0.0_real64,0.0_real64], &
        referencePoints,referenceWeights)
      call near(maxval(abs(p-referencePoints)),0.0_real64,1.0e-14_real64,'Cartesian points')
      call near(maxval(abs(w-referenceWeights)),0.0_real64,1.0e-14_real64,'Cartesian weights')
      deallocate(p,w,referencePoints,referenceWeights)
      call gridLebedev(grid,74_int64)
      call grid%getMetadata(cd,id,units,weightUnits,rule,normalization,degree)
      if(cd.ne.3.or.id.ne.2.or.degree.ne.13) call mqc_error('unitTest12: sphere metadata.')
      allocate(p(3,74),w(74),values(74,2))
      call grid%copyBlock(1_int64,p,w)
      if(minval(w).ge.0) call mqc_error('unitTest12: signed Lebedev weights lost.')
      values(:,1) = 1.0_real64
      values(:,2) = p(1,:)
      call grid%overlapMatrix(values,overlap)
      call near(overlap(1,1),4*Pi,2.0e-12_real64,'sphere constant')
      call near(overlap(1,2),0.0_real64,2.0e-12_real64,'sphere odd moment')
      call near(overlap(2,2),4*Pi/3,2.0e-12_real64,'sphere second moment')
      before = grid%integrate(values(:,1))
      p = 99.0_real64
      w = 0.0_real64
      call near(grid%integrate(values(:,1)),before,1.0e-14_real64,'copy isolation')
      deallocate(p,w,values)
      call gridMappedRadial(radial,96_int64,1.0_real64,'bohr')
      allocate(rp(1,96),rw(96))
      call radial%copyBlock(1_int64,rp,rw)
      call near(radial%integrate(exp(-rp(1,:))),1.0_real64,2.0e-12_real64,'mapped radial')
      call gridLebedev(angular,26_int64)
      allocate(ap(3,26),aw(26))
      call angular%copyBlock(1_int64,ap,aw)
      call gridSphericalProduct(grid,[0.0_real64,0.0_real64,0.0_real64],rp(1,:),rw,ap,aw,'bohr')
      n = grid%pointCount()
      allocate(p(3,n),w(n),centers(3,1),radii(1),atoms(n),base(n),partition(n))
      call grid%copyBlock(1_int64,p,w)
      call near(grid%integrate(exp(-sum(p*p,dim=1))),Pi**1.5_real64,2.0e-10_real64,'Gaussian volume')
      centers = 0.0_real64
      radii = 1.0_real64
      call gridBeckeMolecular(molecule,centers,radii,rp(1,:),rw,ap,aw)
      call molecule%copyPartitionBlock(1_int64,atoms,base,partition)
      if(any(atoms.ne.1)) call mqc_error('unitTest12: generating atoms.')
      call near(maxval(abs(partition-1)),0.0_real64,1.0e-14_real64,'one atom partition')
      call near(maxval(abs(base*partition-w)),0.0_real64,1.0e-14_real64,'molecular weight factorization')
      call molecule%copyAtomicRule(savedRadial,savedRW,savedAngular,savedAW)
      call near(maxval(abs(savedRadial-rp(1,:))),0.0_real64,1.0e-14_real64,'radial provenance')
      call molecule%copyMolecule(centers,radii)
      if(molecule%atomCount().ne.1) call mqc_error('unitTest12: atom count.')
      moleculeCopy = molecule
      allocate(polymorphic,source=molecule)
      call molecule%clear()
      if(molecule%atomCount().ne.0.or.molecule%isInitialized()) call mqc_error('unitTest12: molecular clear.')
      call near(polymorphic%integrate(exp(-sum(p*p,dim=1))),Pi**1.5_real64,2.0e-10_real64,'polymorphic integral')
      call polymorphic%clear()
      if(polymorphic%pointCount().ne.0) call mqc_error('unitTest12: polymorphic clear.')
      call moleculeCopy%validate()
      deallocate(centers,radii,p,w,atoms,base,partition)
      allocate(centers(3,2),radii(2))
      centers = 0.0_real64
      centers(3,:) = [-1.0_real64,1.0_real64]
      radii = [1.0_real64,1.5_real64]
      call gridBeckeMolecular(molecule,centers,radii,rp(1,:),rw,ap,aw)
      n = molecule%pointCount()
      allocate(p(3,n),w(n),atoms(n),base(n),partition(n))
      call molecule%copyBlock(1_int64,p,w)
      call molecule%copyPartitionBlock(1_int64,atoms,base,partition)
      call near(maxval(abs(w-base*partition)/max(1.0_real64,abs(w))), &
        0.0_real64,2.0e-15_real64,'two atom factorization')
      if(any(atoms(1:n/2).ne.1).or.any(atoms(n/2+1:n).ne.2)) call mqc_error('unitTest12: atom order.')
      call molecule%partitionRange(minimum,maximum)
      if(minimum.lt.0.or.maximum.gt.1) call mqc_error('unitTest12: partition bounds.')
      call gridFiniteSpherical(grid,[0.0_real64,0.0_real64,0.0_real64],2.0_real64,8_int64,26_int64, &
        QUAD_RADIAL_GAUSS_LEGENDRE)
      call setup_quadrature_lebedev_spherical([0.0_real64,0.0_real64,0.0_real64],2.0_real64, &
        8_int64,26_int64,QUAD_RADIAL_GAUSS_LEGENDRE,referencePoints,referenceWeights)
      deallocate(p,w)
      allocate(p(3,grid%pointCount()),w(grid%pointCount()))
      call grid%copyBlock(1_int64,p,w)
      call near(maxval(abs(p-referencePoints)),0.0_real64,1.0e-14_real64,'finite product points')
      call near(maxval(abs(w-referenceWeights)),0.0_real64,1.0e-14_real64,'finite product weights')
      call molecularCovariance()
      write(6,'(a)') 'Quadrature object API checks passed.'
      contains
!PROCEDURE molecularCovariance
      subroutine molecularCovariance()
!
!     Rebuild a three-center object under translation and a general rotation.
!     Rotate the supplied angular directions with the molecule; this checks
!     coordinate covariance, not convergence with a fixed laboratory grid.
!     Compare the public geometry/weight metadata and sampled overlap matrix.
!     H. P. Hratchian, 2026.
!
      implicit none
      type(molecular_quadrature_grid)::original,transformed
      real(real64)::c(3,3),ct(3,3),radius(3),rotation(3,3),axis(3),shift(3), &
        skew(3,3),identity(3,3),angle
      real(real64),allocatable::r(:),wr(:),a(:,:),wa(:),points(:,:),weights(:), &
        moved(:,:),movedWeights(:),base0(:),base1(:),part0(:),part1(:), &
        samples(:,:),movedSamples(:,:),overlap0(:,:),overlap1(:,:),back(:,:)
      integer(int64),allocatable::atom0(:),atom1(:)
      integer(int64)::j,k,nPoints,mode
!
      c(:,1) = [-0.7_real64,0.2_real64,0.1_real64]
      c(:,2) = [1.1_real64,-0.3_real64,0.5_real64]
      c(:,3) = [0.2_real64,1.4_real64,-0.4_real64]
      radius = [0.7_real64,1.2_real64,1.5_real64]
      identity = 0.0_real64
      do j = 1,3
        identity(j,j) = 1.0_real64
      endDo
      axis = [1.0_real64,2.0_real64,3.0_real64]/sqrt(14.0_real64)
      skew = reshape([0.0_real64,axis(3),-axis(2),-axis(3),0.0_real64, &
        axis(1),axis(2),-axis(1),0.0_real64],[3,3])
      allocate(r(48),wr(48))
      call setup_quadrature_gauss_legendre_semi_infinite(48_int64,1.0_real64,r,wr)
      call setup_lebedev_angular_grid(110_int64,a,wa)
      call gridBeckeMolecular(original,c,radius,r,wr,a,wa)
      nPoints = original%pointCount()
      allocate(points(3,nPoints),moved(3,nPoints),weights(nPoints),movedWeights(nPoints), &
        base0(nPoints),base1(nPoints),part0(nPoints),part1(nPoints), &
        atom0(nPoints),atom1(nPoints),samples(nPoints,2),movedSamples(nPoints,2))
      call original%copyBlock(1_int64,points,weights)
      call original%copyPartitionBlock(1_int64,atom0,base0,part0)
      do j = 1,2
        samples(:,j) = exp(-real(j,real64)*sum((points-spread(c(:,j),2,nPoints))**2,dim=1))
      endDo
      call original%overlapMatrix(samples,overlap0)
      do mode = 1,2
        rotation = identity
        shift = [1.7_real64,-0.6_real64,2.3_real64]
        if(mode.eq.2) then
          angle = 0.37_real64
          rotation = cos(angle)*identity+sin(angle)*skew
          do j = 1,3
            do k = 1,3
              rotation(j,k) = rotation(j,k)+(1.0_real64-cos(angle))*axis(j)*axis(k)
            endDo
          endDo
          shift = 0.0_real64
        endIf
        ct = matmul(rotation,c)+spread(shift,2,3)
        call gridBeckeMolecular(transformed,ct,radius,r,wr,matmul(rotation,a),wa)
        call transformed%copyBlock(1_int64,moved,movedWeights)
        call transformed%copyPartitionBlock(1_int64,atom1,base1,part1)
        back = matmul(transpose(rotation),moved-spread(shift,2,nPoints))
        call near(maxval(abs(back-points)/max(1.0_real64,abs(points))), &
          0.0_real64,2.0e-11_real64,'molecular coordinate covariance')
!       Scale by the base volume weight: nearly zero partition factors at
!       distant points make relative errors in final weights ill-conditioned.
        call near(maxval(abs(movedWeights-weights)/max(1.0_real64,abs(base0))), &
          0.0_real64,2.0e-11_real64,'molecular weight covariance')
        call near(maxval(abs(base1-base0)/max(1.0_real64,abs(base0))), &
          0.0_real64,2.0e-14_real64,'molecular base weight covariance')
        call near(maxval(abs(part1-part0)),0.0_real64,2.0e-11_real64,'partition covariance')
        if(any(atom0.ne.atom1)) call mqc_error('unitTest12: covariance generating atoms.')
        do j = 1,2
          movedSamples(:,j) = exp(-real(j,real64)*sum((back-spread(c(:,j),2,nPoints))**2,dim=1))
        endDo
        call transformed%overlapMatrix(movedSamples,overlap1)
        call near(maxval(abs(overlap1-overlap0)),0.0_real64,2.0e-12_real64,'overlap covariance')
      endDo
      write(6,'(a)') 'Molecular object translation/rotation covariance checks passed.'
      end subroutine molecularCovariance
!PROCEDURE near
      subroutine near(actual,expected,tolerance,label)
      implicit none
      real(real64),intent(in)::actual,expected,tolerance
      character(*),intent(in)::label
      if(abs(actual-expected).gt.tolerance) then
        write(6,'(a,3es22.12)') trim(label)//': actual, expected, tolerance ',actual,expected,tolerance
        call mqc_error('unitTest12: '//label)
      endIf
      end subroutine near
!PROCEDURE invalidInput
      subroutine invalidInput(which)
      implicit none
      integer(int64),intent(in)::which
      real(real64)::buffer(3,1),v(1)
      type(quadrature_grid)::bad
      v = 1.0_real64
      buffer = 0.0_real64
      if(which.eq.1) then
        x = bad%integrate(v)
        return
      endIf
      call gridLebedev(bad,6_int64)
      select case(which)
      case(2)
        call bad%copyBlock(7_int64,buffer)
      case(3)
        x = bad%integrate(v)
      case(4)
        call gridFromArrays(bad,buffer,v,4_int64,'bohr','bohr**3','supplied','volume')
      case(5)
        buffer(1,1) = ieee_value(0.0_real64,ieee_quiet_nan)
        call gridFromArrays(bad,buffer,v,3_int64,'bohr','bohr**3','supplied','volume')
      case(6)
        call gridFromArrays(bad,buffer,v,3_int64,'','bohr**3','supplied','volume')
      case(7)
        call gridGaussLegendre1D(bad,0_int64,0.0_real64,1.0_real64,'bohr')
      case(8)
        call bad%copyBlock(1_int64,buffer(1:2,:))
      end select
      end subroutine invalidInput
      end program unitTest12
