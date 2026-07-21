      program unitTest01
!
!     Validate MQCPack version reporting and inclusive version comparisons.
!     Every expectation is checked directly so Automake receives a nonzero
!     status if the version routines become inconsistent.
!
!     H. P. Hratchian, 2020, 2026.
!
      use iso_fortran_env, only: int64
      use mqc_general
      implicit none
!
      integer(kind=int64)::major,minor,revision
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


      end program unitTest01
