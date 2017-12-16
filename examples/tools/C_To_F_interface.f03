!
! C to Fortran 2003 Interface
!
      subroutine mqc_error_i_c2f_0 ( in_message, iout ) bind ( C, name="mqc_error_i_c2f_0" )
! Interface from C to FORTRAN to call mqc_error_i
! with 0 of the optional arguments

        use iso_c_binding
        use mqc_general
        implicit none

        character (kind=c_char, len=1), dimension (2048), intent (in) :: in_message
        integer (kind=c_int), intent (in) :: iout
        character (len=2096) :: out_message
        integer :: i
        integer :: IOUT_F

        out_message = " "
        message_string: do i=1, 2096
           if ( in_message (i) == c_null_char ) then
              exit message_string
           else
              out_message(i:i) = in_message (i)
           end if
        end do message_string

        IOUT_F = iout
        call mqc_error_i( out_message, IOUT_F)
        return
      end subroutine mqc_error_i_c2f_0

      subroutine print_line_c2f ( in_message, iout ) bind ( C, name="print_line_c2f" )
! Interface from C to FORTRAN to call mqc_error_i
! with 0 of the optional arguments
! Insure that All I/O to unit 6 happens through Fortran

        use iso_c_binding
        implicit none

        character (kind=c_char, len=1), dimension (2048), intent (in) :: in_message
        integer (kind=c_int), intent (in) :: iout
        character (len=2096) :: out_message
        integer :: i
        integer :: IOUT_F

        out_message = " "
        message_string: do i=1, 2096
           if ( in_message (i) == c_null_char ) then
              exit message_string
           else
              out_message(i:i) = in_message (i)
           end if
        end do message_string
        IOUT_F = iout
        write( IOUT_F, 1010 ) out_message
 1010   format( A )

        return

      end subroutine print_line_c2f

      subroutine flush_c2f ( iout ) bind ( C, name="flush_c2f" )
        use iso_c_binding
        implicit none

        integer (kind=c_int), intent (in) :: iout
        integer :: IOUT_F

        flush(IOUT_F)
        return
      end subroutine flush_c2f

!
! Fortran 2003 to C Interface
!
    subroutine mqc_get_MatrixFile_Names_F2C(FileName_F, Program_F, iout_F)
      use iso_c_binding
      implicit none

      interface
         subroutine mqc_get_MatrixFile_Names(FileName_C, Program_C, iout_C) bind(C, name="mqc_get_MatrixFile_Names")
           import
           character(kind=c_char) :: FileName_C
           character(kind=c_char) :: Program_C
           integer(kind=c_int), value :: iout_C
         end subroutine mqc_get_MatrixFile_Names
      end interface

      character(len=*), intent(inout)::FileName_F
      character(len=*), intent(inout)::Program_F
      integer, intent(inout) :: iout_F
      character(len=len_trim(FileName_F)+1) :: FileName_CALL
      character(len=len_trim(Program_F)+1) :: Program_CALL
      integer(kind=c_int) :: iout_CALL

      FileName_CALL = trim(FileName_F) // c_null_char
      Program_CALL = trim(Program_F) // c_null_char
      iout_CALL = iout_F

      call mqc_get_MatrixFile_Names(FileName_CALL, &
           Program_CALL, iout_CALL)
    end subroutine mqc_get_MatrixFile_Names_F2C

