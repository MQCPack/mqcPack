      Program unitTest01
!
!     unitTest01: This unit test program tests the version reporting and
!     checking routines in MQC_General.
!
!     -H. P. Hratchian, 2020.
!
!
!     USE Connections
!
      use mqc_general
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::major,minor,revision
      character(len=512)::versionString
!
!     Format Statements
!
 1000 Format(1x,'Major Version Number: ',I2,/,  &
        1x,'Minor Version Number: ',I2,2x,'Revision Number: ',I2)
 1100 Format(/,1x,'Full Version String: <',A,'>',/)
 2000 Format(1x,'Is 2020:     ',L1)
 2010 Format(1x,'Is 2020.5:   ',L1)
 2020 Format(1x,'Is 2020.5.1: ',L1)
 2100 Format(1x,'Is 2020:     ',L1)
 2110 Format(1x,'Is 2020.4:   ',L1)
 2120 Format(1x,'Is 2020.4.0: ',L1)
 3000 Format(1x,'Is NEWER than 2018: ',L1)
 3100 Format(1x,'Is OLDER than 2018: ',L1)
 9000 Format(/,1x,'END OF unitTest01')

!
!
      call mqc_version(major,minor,revision,versionString)
      write(*,1000) major,minor,revision
      write(*,1100) TRIM(versionString)
!
!
      write(*,2000) mqc_version_check(isMajor=20)
      write(*,2010) mqc_version_check(isMajor=20,isMinor=5)
      write(*,2020) mqc_version_check(isMajor=20,isMinor=5,isRevision=1)
!
!
      write(*,2100) mqc_version_check(isMajor=20)
      write(*,2110) mqc_version_check(isMajor=20,isMinor=4)
      write(*,2120) mqc_version_check(isMajor=20,isMinor=4,isRevision=0)
!
!
      write(*,3000) mqc_version_check(newerThanMajor=18)
      write(*,3100) mqc_version_check(olderThanMajor=18)
!
      write(*,9000)
      end program unitTest01
