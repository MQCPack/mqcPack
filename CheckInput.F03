! This program tests Reading everthing in a MartixFile
     program CheckGauInput
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
      use iso_fortran_env

      implicit none

      character(len=1024)::program_name
      character(len=1024)::no_program_name
      character(len=1024)::FileName
      integer(kind=int64) cmd_count
!
      integer(kind=int64) iout, i, status

 1010 Format( A )
 1020 Format( '   - ',A )
      cmd_count = COMMAND_ARGUMENT_COUNT()
      iout = 6
      if ( cmd_count .eq. 0 ) then
         write( iout, 1010 ) "Error in number of argments to this program."
         write( iout, 1010 ) "Arguments should be:"
         write( iout, 1010 ) "     - Full path to a Gaussian file"
         write( iout, 1010 ) "        - This can either be FormChk File"
         write( iout, 1010 ) "        - This can either be a MatrixFile"
         write( iout, 1010 ) "        - This can be a Gaussian input file"
         write( iout, 1010 ) "     - Optional, only required if the first "
         write( iout, 1010 ) "          argument is a Gaussian input file."
         write( iout, 1010 ) "          The name of the binary for the "
         write( iout, 1010 ) "          program that generates the"
         write( iout, 1010 ) "          FormChk/MatrixFile, most likely"
         write( iout, 1010 ) "          'g16'. The environment should be"
         write( iout, 1010 ) "          setup to run the program."
         call mqc_error_i( "No arguments to this program", iout)
      endif
      if ( cmd_count.eq. 1 ) then
         no_program_name = 'Unknown'
         call get_command_argument(1, FileName)
         call mqc_create_File_Name_lists_F2C(FileName, no_program_name, iout)
      else
         call get_command_argument(cmd_count, program_name)
         call mqc_File_Is_Executable_F2C( program_name, Status, iout )
         if ( Status .eq. 1 ) then
            do i = 1, cmd_count-1
               call get_command_argument(i, FileName)
               call mqc_create_File_Name_lists_F2C(FileName,Program_name,iout)
            end do
         else
            no_program_name = 'Unknown'
            do i = 1, cmd_count
               call get_command_argument(i, FileName)
               call mqc_create_File_Name_lists_F2C(FileName,no_program_name,iout)
            end do
         endif
      endif
#ifdef MATRIXFILE
      call  mqc_get_MatrixFile_Name_F2C(FileName, iout)
      if ( trim(FileName).ne.'Done with MatrixFiles' ) then
         write(iout,1010) 'MatrixFiles:'
      endif
      Do While( trim(FileName).ne.'Done with MatrixFiles' )
         write(iout,1020) trim(FileName)
         call  mqc_get_MatrixFile_Name_F2C(FileName, iout)
      EndDo
#endif
      call  mqc_get_FormChk_Name_F2C(FileName, iout)
      if ( trim(FileName).ne.'Done with FormChk' ) then
         write(iout,1010) 'FormChk Files:'
      endif
      Do While( trim(FileName).ne.'Done with FormChk' )
         write(iout,1020) trim(FileName)
         call  mqc_get_FormChk_Name_F2C(FileName, iout)
      EndDo

    end program CheckGauInput
