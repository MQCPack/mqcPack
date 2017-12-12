! This program tests Reading everthing in a MartixFile
     program readWriteInfo
!
!     This program contains examples of the matrix file read and write routines
!
!     Variable Declarations...       
!
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
      Use MQC_EST
      use mqc_gaussian
      use mqc_gaussian_read_write

      implicit none

      character(len=:),allocatable::fileName
      class(mqc_link_list),allocatable::self
!
!     Command line arguments can be returned as strings of correctly allocated length using
!     mqc_get_command_argument
!
      call mqc_get_command_argument(1,FileName)
      allocate(self)
      call self%add(FileName)

      end program readWriteInfo

