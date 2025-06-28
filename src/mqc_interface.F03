!
!     ***********************************************************************
!     ***********************************************************************
!     **                                                                   **
!     **               The Merced Quantum Chemistry Package                **
!     **                            (MQCPack)                              **
!     **                                                                   **
!     ** Written By:                                                       **
!     **    Lee M. Thompson, Xianghai Sheng, Andrew D. Mahler, Dave        **
!     **    Mullally, and Hrant P. Hratchian                               **
!     **                                                                   **
!     **                          Version 25.6.3                           **
!     **                           June 18, 2025                           **
!     **                                                                   **
!     **                                                                   **
!     ** Modules beloning to MQCPack:                                      **
!     **    1. MQC_General                                                 **
!     **    2. MQC_Binary                                                  **
!     **    3. MQC_DataStructures                                          **
!     **    4. MQC_Integrals                                               **
!     **    5. MQC_Algebra                                                 **
!     **    6. MQC_Files                                                   **
!     **    7. MQC_Molecule                                                **
!     **    8. MQC_EST                                                     **
!     **    9. MQC_Gaussian                                                **
!     **                                                                   **
!     ***********************************************************************
!     ***********************************************************************
!                                                                           
!                                                                           
!                                                                           
!
! C to Fortran 2003 Interface
!
      subroutine mqc_error_i_c2f_0 ( in_message, iout ) bind ( C, name="mqc_error_i_c2f_0" )
! Interface from C to FORTRAN to call mqc_error_i
! with 0 of the optional arguments

        use iso_c_binding
        use mqc_general
        use iso_fortran_env, only: int32, int64, real64
        implicit none

        character (kind=c_char, len=1), dimension (2048), intent (inout) :: in_message
        integer (kind=c_int), intent (inout) :: iout
        character (len=2048) :: out_message
        integer( kind=int32 ) :: i
        integer( kind=int64 ) :: IOUT_F
 1010   Format( "IOut should be 0 or 6.  It is %d" )
        out_message = " "
        message_string: do i=1, 2048
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
        use iso_fortran_env, only: int32, int64, real64
        implicit none

        character (kind=c_char, len=1), dimension (2048), intent (in) :: in_message
        integer (kind=c_int), intent (in) :: iout
        character (len=2048) :: out_message
        integer( kind=int32 ) :: i
        integer( kind=int64 ) :: IOUT_F

        out_message = " "
        message_string: do i=1, 2048
           if ( in_message (i) == c_null_char ) then
              exit message_string
           else
              out_message(i:i) = in_message (i)
           end if
        end do message_string
        IOUT_F = iout
        return
      end subroutine print_line_c2f

      subroutine flush_c2f ( iout ) bind ( C, name="flush_c2f" )
        use iso_c_binding
        use iso_fortran_env, only: int32, int64, real64
        implicit none

        integer (kind=c_int), intent (in) :: iout
        integer( kind=int64 ) :: IOUT_F

        flush(IOUT_F)
        return
      end subroutine flush_c2f

!
! Fortran 2003 to C Interface
!
    subroutine mqc_create_File_Name_lists_F2C(FileName_F, Program_F, iout_F)
      use iso_c_binding
      use iso_fortran_env, only: int32, int64, real64
      implicit none

      interface
         subroutine mqc_create_File_Name_lists(FileName_C, Program_C, iout_C) bind(C, name="mqc_create_File_Name_lists")
           import
           character(kind=c_char) :: FileName_C
           character(kind=c_char) :: Program_C
           integer(kind=c_int), value :: iout_C
         end subroutine mqc_create_File_Name_lists
      end interface

      character(len=*), intent(inout)::FileName_F
      character(len=*), intent(inout)::Program_F
      integer( kind=int64 ), intent(inout) :: iout_F
      character(kind=c_char, len=len_trim(FileName_F)+1) :: FileName_CALL
      character(kind=c_char, len=2048) :: Program_CALL
      integer(kind=c_int) :: iout_CALL
      character(len=2048) :: Program_Local

      FileName_CALL = trim(FileName_F) // c_null_char
      Program_CALL = trim(Program_F) // c_null_char
      iout_CALL = iout_F

      call mqc_create_File_Name_lists(FileName_CALL, Program_CALL, iout_CALL)
           
    end subroutine mqc_create_File_Name_lists_F2C

#ifdef MATRIXFILE
    subroutine mqc_simplify_formula_F2C(Formula_F)
      use iso_c_binding
      use iso_fortran_env, only: int32, int64, real64
      implicit none

      interface
         subroutine mqc_simplify_formula(Formula_C) bind(C, name="mqc_Simplify_Formula")
           import
           character(kind=c_char) :: Formula_C
         end subroutine mqc_simplify_formula
      end interface

      character(len=*), intent(inout)::Formula_F
      character(kind=c_char, len=len_trim(Formula_F)+1) :: Formula_CALL
      integer i

      Formula_CALL = trim(Formula_F) // c_null_char

      call mqc_simplify_formula(Formula_CALL)
           
      Formula_F = " "
      DataType_string: do i=1, 2048
         if ( Formula_CALL(i:i) .eq. c_null_char ) then
            exit DataType_string
         else
            Formula_F(i:i) = Formula_CALL(i:i)
         end if
      end do DataType_string

    end subroutine mqc_simplify_formula_F2C

    subroutine mqc_Matfile_Guide_F2C(FileName_F, DataType_F, Char_Len_F, &
         rec_len_F, bufsize_F )

      use iso_c_binding
      use iso_fortran_env, only: int32, int64, real64
      implicit none

      interface
         subroutine mqc_Matfile_Guide(FileName_C, DataType_C, Char_Len_C, &
              rec_len_C, bufsize_C ) bind(C, name="mqc_Matfile_Guide")
           import
           character(kind=c_char) :: FileName_C
           character(kind=c_char) :: DataType_C
           integer(kind=c_int) :: Char_Len_C
           integer(kind=c_int) :: rec_len_C
           integer(kind=c_int) :: bufsize_C
         end subroutine mqc_Matfile_Guide
      end interface

      character(len=*), intent(in)::FileName_F
      character(len=*), intent(out)::DataType_F
      integer( kind=int64 ), intent(in) :: Char_Len_F
      integer( kind=int64 ), intent(out) :: rec_len_F
      integer( kind=int64 ), intent(out) :: bufsize_F

      character(kind=c_char, len=len_trim(FileName_F)+1) :: FileName_CALL
      character(kind=c_char, len=2048) :: DataType_CALL
      integer(kind=c_int) :: Char_Len_CALL
      integer(kind=c_int) :: rec_len_CALL
      integer(kind=c_int) :: bufsize_CALL
      integer(kind=c_int) :: i

      FileName_CALL = trim(FileName_F) // c_null_char
      Char_Len_CALL = Char_Len_F

      call mqc_Matfile_Guide(FileName_CALL, DataType_CALL, Char_Len_CALL, rec_len_CALL, bufsize_CALL )

      DataType_F = " "
      DataType_string: do i=1, 2048
         if ( DataType_CALL(i:i) .eq. c_null_char ) then
            exit DataType_string
         else
            DataType_F(i:i) = DataType_CALL(i:i)
         end if
      end do DataType_string
      rec_len_F = rec_len_CALL
      bufsize_F = bufsize_CALL

    end subroutine mqc_Matfile_Guide_F2C

    subroutine mqc_get_MatrixFile_Name_F2C(FileName_F, iout_F)
      use iso_c_binding
      use iso_fortran_env, only: int32, int64, real64
      implicit none

      interface
         subroutine mqc_get_MatrixFile_Name(FileName_C, iout_C) bind(C, name="mqc_get_MatrixFile_Name")
           import
           character(kind=c_char) :: FileName_C
           integer(kind=c_int), value :: iout_C
         end subroutine mqc_get_MatrixFile_Name
      end interface

      character(len=*), intent(inout)::FileName_F
      integer( kind=int64 ), intent(inout) :: iout_F
      character(kind=c_char, len=2048) :: FileName_CALL
      integer(kind=c_int) :: iout_CALL
      integer( kind=int32 ) :: i, j

      FileName_CALL = " "
      iout_CALL = iout_F

      call mqc_get_MatrixFile_Name(FileName_CALL, iout_CALL)

      FileName_F = " "
      FileName_string: do i=1, 2048
         if ( FileName_CALL(i:i) .eq. c_null_char ) then
            exit FileName_string
         else
            FileName_F(i:i) = FileName_CALL(i:i)
         end if
      end do FileName_string

    end subroutine mqc_get_MatrixFile_Name_F2C
#endif

    subroutine mqc_get_FormChk_Name_F2C(FileName_F, iout_F)
      use iso_c_binding
      use iso_fortran_env, only: int32, int64, real64
      implicit none

      interface
         subroutine mqc_get_FormChk_Name(FileName_C, iout_C) bind(C, name="mqc_get_FormChk_Name")
           import
           character(kind=c_char) :: FileName_C
           integer(kind=c_int), value :: iout_C
         end subroutine mqc_get_FormChk_Name
      end interface

      character(len=*), intent(inout)::FileName_F
      integer( kind=int64 ), intent(inout) :: iout_F
      character(kind=c_char, len=2048) :: FileName_CALL
      integer(kind=c_int) :: iout_CALL
      integer( kind=int32 ):: i, j

      FileName_CALL = " "
      iout_CALL = iout_F

      call mqc_get_FormChk_Name(FileName_CALL, iout_CALL)

      FileName_F = " "
      FileName_string: do i=1, 2048
         if ( FileName_CALL(i:i) .eq. c_null_char ) then
            exit FileName_string
         else
            FileName_F(i:i) = FileName_CALL(i:i)
         end if
      end do FileName_string

    end subroutine mqc_get_FormChk_Name_F2C

    subroutine mqc_File_Exists_F2C(FileName_F, Status_F)
      use iso_c_binding
      use iso_fortran_env, only: int32, int64, real64
      implicit none

      interface
         subroutine mqc_File_Exists(FileName_C, Status_C) bind(C, name="mqc_File_Exists")
           import
           character(kind=c_char) :: FileName_C
           integer(kind=c_int) :: Status_C
         end subroutine mqc_File_Exists
      end interface

      character(len=*), intent(inout)::FileName_F
      integer( kind=int64 ), intent(inout) :: Status_F
      character(kind=c_char, len=2048) :: FileName_CALL
      integer(kind=c_int) :: Status_CALL

      FileName_CALL = trim(FileName_F) // c_null_char
      call mqc_File_Exists(FileName_CALL, Status_CALL)
      Status_F = Status_CALL

    end subroutine mqc_File_Exists_F2C

    subroutine mqc_File_Is_Executable_F2C(FileName_F, Status_F, iout_F)
      use iso_c_binding
      use iso_fortran_env, only: int32, int64, real64
      implicit none

      interface
         subroutine mqc_File_Is_Executable(FileName_C, Status_C, iout_C) bind(C, name="mqc_File_Is_Executable")
           import
           character(kind=c_char) :: FileName_C
           integer(kind=c_int) :: Status_C
           integer(kind=c_int), value :: iout_C
         end subroutine mqc_File_Is_Executable
      end interface

      character(len=*), intent(inout)::FileName_F
      integer( kind=int64 ), intent(inout) :: Status_F
      integer( kind=int64 ), intent(inout) :: iout_F
      character(kind=c_char, len=2048) :: FileName_CALL
      integer(kind=c_int) :: Status_CALL
      integer(kind=c_int) :: iout_CALL

      FileName_CALL = trim(FileName_F) // c_null_char
      iout_CALL = iout_F
      call mqc_File_Is_Executable(FileName_CALL, Status_CALL, iout_CALL)
      Status_F = Status_CALL

      end subroutine mqc_File_Is_Executable_F2C

#ifdef MATRIXFILE
!
! This routine is here, because if violates typing.  Want to copy all 
! the bits of integerTmp to arrayTmp.  This is a way to get real and integer
! data from a single read of the MatrixFile
! arrayTmp is was real in the calling subroutine.
! integerTmp is was integer in the calling subroutine.
!
      subroutine mqc_copy_vec( integerTmp8, integerTmp4, arrayTmp, LR)
      use MQC_MatWrapper
      use iso_fortran_env, only: int32, int64, real64
      Implicit None
      integer(kind=int64),dimension(:),intent(inout)::integerTmp8
      integer(kind=int32),dimension(:),intent(in)::integerTmp4
      real(kind=real64),dimension(:),intent(out)::arrayTmp
      integer(kind=int64),intent(in):: LR
      integer(kind=int64):: i
      Character(Len=64)::Ignorable
      Integer(kind=int32)::IntPrecision
      Integer(kind=int32)::Call_Test

      Call_Test=0 
      call MQC_Mat_Precision(Ignorable, IntPrecision, Call_Test)

      write(6,1234) lr
1234  format( i10 )
        flush(6)
      write(6,1234) integerTmp8(1) 
        flush(6)
      return

      do i=1, LR
      write(6,1234) integerTmp8(i) 
        flush(6)
         arrayTmp(i) = integerTmp8(i) 
      end do

      return

      if ( IntPrecision .eq. 4 ) then
!
! These 2 arrays are the same, with different size declarations.
! Starting from the end allows us to expand the array without overwriting data
!
         do i = LR, 1, -1
            integerTmp8(i) = integerTmp4(i)
         end do
      endif
      return
      end subroutine mqc_copy_vec
#endif
