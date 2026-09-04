      program unitTest01
!
!     Validate MQCPack general utilities. Every expectation is checked directly
!     so Automake receives a nonzero status if a result is inconsistent.
!
!     H. P. Hratchian, 2020, 2026.
!
      use iso_fortran_env, only: int64,real64
      use mqc_general
      implicit none
!
      integer(kind=int64)::major,minor,revision
      real(kind=real64),dimension(4)::point1,point2
      real(kind=real64),dimension(2,3)::points
      real(kind=real64),dimension(:,:),allocatable::distanceMatrix
      character(len=512)::versionString
      character(len=64)::expectedVersion
!
      call mqc_version(major,minor,revision,versionString)
      write(expectedVersion,'(I0,".",I0,".",I0)') major,minor,revision
!
      call assertTrue(major.ge.0_int64.and.major.le.99_int64,  &
        'The major version is not a two-digit calendar year.')
      call assertTrue(minor.ge.1_int64.and.minor.le.12_int64,  &
        'The minor version is not a calendar month.')
      call assertTrue(revision.ge.0_int64,  &
        'The revision counter is negative.')
      call assertTrue(TRIM(versionString).eq.TRIM(expectedVersion),  &
        'The full version string disagrees with its components.')
!
      call assertTrue(mqc_version_check(isMajor=major),  &
        'The exact major-version check failed.')
      call assertTrue(mqc_version_check(isMajor=major,isMinor=minor),  &
        'The exact major/minor-version check failed.')
      call assertTrue(mqc_version_check(isMajor=major,isMinor=minor,  &
        isRevision=revision),'The exact full-version check failed.')
      call assertTrue(.not.mqc_version_check(isMajor=major,isMinor=minor,  &
        isRevision=revision+1_int64),  &
        'An incorrect revision was accepted as an exact match.')
!
!     Older-than and newer-than comparisons are inclusive.
!
      call assertTrue(mqc_version_check(newerThanMajor=major),  &
        'The inclusive newer-than check rejected the current major version.')
      call assertTrue(mqc_version_check(olderThanMajor=major),  &
        'The inclusive older-than check rejected the current major version.')
      call assertTrue(mqc_version_check(newerThanMajor=major-1_int64),  &
        'The newer-than check rejected an older major version.')
      call assertTrue(mqc_version_check(olderThanMajor=major+1_int64),  &
        'The older-than check rejected a newer major version.')
      call assertTrue(.not.mqc_version_check(newerThanMajor=  &
        major+1_int64),'The newer-than check accepted a newer version.')
      call assertTrue(.not.mqc_version_check(olderThanMajor=  &
        major-1_int64),'The older-than check accepted an older version.')
!
!     Euclidean point distances support arbitrary positive dimension. Pairwise
!     distance matrices treat each column of the input matrix as one point.
!
      point1 = [0.0_real64,0.0_real64,0.0_real64,0.0_real64]
      point2 = [1.0_real64,2.0_real64,2.0_real64,4.0_real64]
      call assertNear(mqc_distance(point1,point2),5.0_real64,  &
        'The four-dimensional point distance is incorrect.')
!
      points(:,1) = [0.0_real64,0.0_real64]
      points(:,2) = [3.0_real64,4.0_real64]
      points(:,3) = [3.0_real64,0.0_real64]
      distanceMatrix = mqc_distance_matrix(points)
      call assertTrue(all(Shape(distanceMatrix).eq.[3,3]),  &
        'The pairwise distance matrix has the wrong shape.')
      call assertNear(distanceMatrix(1,1),0.0_real64,  &
        'A distance-matrix diagonal element is not zero.')
      call assertNear(distanceMatrix(1,2),5.0_real64,  &
        'The first pairwise distance is incorrect.')
      call assertNear(distanceMatrix(1,3),3.0_real64,  &
        'The second pairwise distance is incorrect.')
      call assertNear(distanceMatrix(2,3),4.0_real64,  &
        'The third pairwise distance is incorrect.')
      call assertNear(distanceMatrix(2,1),distanceMatrix(1,2),  &
        'The pairwise distance matrix is not symmetric.')
!
      write(*,'(1x,A)') 'unitTest01: PASS'
!
      contains


!PROCEDURE assertTrue
      subroutine assertTrue(condition,message)
      implicit none
      logical,intent(in)::condition
      character(len=*),intent(in)::message
!
      if(.not.condition) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(message)
        error stop 1
      endIf
!
      return
      end subroutine assertTrue


!PROCEDURE assertNear
      subroutine assertNear(actual,expected,message)
      implicit none
      real(kind=real64),intent(in)::actual,expected
      character(len=*),intent(in)::message
!
      if(Abs(actual-expected).gt.1.0e-12_real64) then
        write(*,'(1x,A,2(1x,ES24.16))') 'FAIL: '//Trim(message),  &
          actual,expected
        error stop 1
      endIf
!
      return
      end subroutine assertNear


      end program unitTest01
