! This program tests Reading everthing in a MartixFile
     program DataSummary
!
!     This program contains examples of the matrix file read and write routines
!
!     Variable Declarations...       
!
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
      Use MQC_EST
      use mqc_General
!hph      use mqc_FullWavefunction
      use mqc_matwrapper
      use iso_fortran_env
      implicit none

!hph+
!      character(len=1024)::program_name
!      character(len=1024)::no_program_name
!      character(len=1024)::FileName
!      integer(kind=int64) cmd_count
!      integer(kind=int64) IAllocated
!!
!      integer(kind=int64) iout, i, status
!      integer(kind=int64) Print_Options
!      integer(kind=int64) Storage_Options
!      class(mqc_link_list), allocatable:: mqc_input_file_list
!      class(mqc_READONLY_list), allocatable:: my_READONLY_list
!
! 1010 Format( A )
! 1020 Format( '   - ',A )
!hph-

!hph+
      Write(*,*)' Hrant - This example program is obsolescent.'
!      IOut = 6
!      allocate(mqc_input_file_list)
!      allocate(my_READONLY_list)
!      Print_Options = 1
!      Storage_Options = 0
!      Call MQC_Create_Wavefunctions( Print_Options, Storage_Options, &
!           my_READONLY_list, mqc_input_file_list, IOut)
!
!!
!! The rest of this is stubs for when we can read FCHK files.
!!
!      call  mqc_get_FormChk_Name_F2C(FileName, iout)
!      if ( trim(FileName).ne.'Done with FormChk' ) then
!         write(iout,1010) 'FormChk Files:'
!      endif
!      if ( trim(FileName).ne.'Done with FormChk' ) then
!         if ( .not.allocated(mqc_input_file_list) ) then
!            allocate(mqc_input_file_list)
!         endif
!      endif
!      Do While( trim(FileName).ne.'Done with FormChk' )
!         write(iout,1020) trim(FileName)
!         call  mqc_get_FormChk_Name_F2C(FileName, iout)
!      EndDo
!hph-

    end program DataSummary
