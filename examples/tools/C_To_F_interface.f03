
subroutine mqc_error_i_c2f_0 ( in_message ) bind ( C, name="mqc_error_i_c2f_0" )
! Interface from C to FORTRAN to call mqc_error_i
! with 0 of the optional arguments

   use iso_c_binding, only: C_CHAR, c_null_char
   use mqc_general
   implicit none

   character (kind=c_char, len=1), dimension (2048), intent (in) :: input_string
   character (len=2096) :: regular_string
   integer :: i

   regular_string = " "
   loop_string: do i=1, 2096
      if ( input_string (i) == c_null_char ) then
         exit loop_string
      else
         regular_string (i:i) = input_string (i)
      end if
   end do loop_string

mqc_error_i(
, 6)
return

end subroutine myfor
   write (*, *) ">", trim (regular_string), "<", len_trim (regular_string)
return

end subroutine myfortsub