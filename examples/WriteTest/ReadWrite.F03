! This program tests Reading everthing in a MartixFile
     program ReadWrite
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
      use mqc_FullWavefunction
      use mqc_matwrapper
      use iso_fortran_env
      implicit none

      character(len=1024)::program_name
      character(len=1024)::no_program_name
      character(len=64)::FileName
!
      integer(kind=int64) iout, i, status
      integer(kind=int64) Print_Options
      class(MQC_ALL_Wavefunction),pointer::current_wavefunction
      class(mqc_link_list), allocatable:: mqc_input_file_list
      class(mqc_READONLY_list), allocatable:: my_READONLY_list

 1010 Format( A )
 1020 Format( '   - ',A )
 1030 Format( 2A )

      iout = 6
      allocate(mqc_input_file_list)
      allocate(my_READONLY_list)
      Print_Options=0

      Call MQC_Create_Wavefunctions( Print_Options, my_READONLY_list, &
           mqc_input_file_list, Iout)

      call mqc_link_list_start(mqc_input_file_list,Status,current_wavefunction)
      Do While( Status.eq.1 )
         write(FileName,1030) trim(current_wavefunction%FileName),".copy"
         Current_wavefunction%FileName=FileName
         call mqc_WRITE_all_wavefunction(current_wavefunction,iOut)
         call mqc_link_list_next(mqc_input_file_list,Status, &
              current_wavefunction)
      EndDo

    end program ReadWrite

