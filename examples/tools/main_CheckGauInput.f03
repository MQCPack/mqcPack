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

      implicit none

      character(len=1024)::program_name
      character(len=1024)::FileName
      integer cmd_count
!
!
      integer iout, i

 1010 Format( A )
      cmd_count = COMMAND_ARGUMENT_COUNT()

      iout = 6

      if ( cmd_count .lt. 2 ) then
         write( iout, 1010 ) "Error in number of argments to this program."
         write( iout, 1010 ) "Arguments should be:"
         write( iout, 1010 ) "     - Full path to a Gaussian file"
         write( iout, 1010 ) "        - This can either be FormCheck File"
         write( iout, 1010 ) "        - This can either be a MatrixFile"
         write( iout, 1010 ) "        - This can be a Gaussian input file"
         write( iout, 1010 ) "     - The name of the binary for the program that"
         write( iout, 1010 ) "       generates the MatrixFile, most likely 'g16'."
         write( iout, 1010 ) "       The environment should be setup to run the program."
         if ( cmd_count .eq. 1 ) then
            call mqc_error_i( "No arguments to this program", iout)
         else
            call mqc_error_i(  "Only one argument to this program", iout)
         endif
      endif

      call get_command_argument(cmd_count, program_name)
      do i = 1, cmd_count-1
         call get_command_argument(i, FileName)
         call mqc_get_MatrixFile_Names_F2C(FileName, Program_name, iout)
      end do

    end program CheckGauInput
