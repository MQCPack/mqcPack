!
!     **********************************************************************
!     **********************************************************************
!     **                                                                  **
!     **               The Merced Quantum Chemistry Package               **
!     **                            (MQCPack)                             **
!     **                       Development Version                        **
!     **                            Based On:                             **
!     **                     Development Version 0.1                      **
!     **                                                                  **
!     **                                                                  **
!     ** Written By:                                                      **
!     **    Lee M. Thompson, Xianghai Sheng, and Hrant P. Hratchian       **
!     **                                                                  **
!     **                                                                  **
!     **                      Version 1.0 Completed                       **
!     **                           May 1, 2017                            **
!     **                                                                  **
!     **                                                                  **
!     ** Modules beloning to MQCPack:                                     **
!     **    1. MQC_General                                                **
!     **    2. MQC_DataStructures                                         **
!     **    3. MQC_Algebra                                                **
!     **    4. MQC_Files                                                  **
!     **    5. MQC_Molecule                                               **
!     **    6. MQC_EST                                                    **
!     **    7. MQC_Gaussian                                               **
!     **                                                                  **
!     **********************************************************************
!     **********************************************************************
!                                                                         
!                                                                           
!                                                                           
!                                                                           
!     This module includes procedures that provide a variety of general purpose
!     utilities. The subroutines and functions provided by this module are 
!     grouped into the following sections:  
!
! C to Fortran 2003 Interface
!
      subroutine mqc_error_i_c2f_0 ( in_message, iout ) bind ( C, name="mqc_error_i_c2f_0" )
! Interface from C to FORTRAN to call mqc_error_i
! with 0 of the optional arguments

        use iso_c_binding
        use mqc_general
        use iso_fortran_env
        implicit none

        character (kind=c_char, len=1), dimension (2048), intent (inout) :: in_message
        integer (kind=c_int), intent (inout) :: iout
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
        call mqc_error_i( out_message, IOUT_F)
        return
      end subroutine mqc_error_i_c2f_0

      subroutine print_line_c2f ( in_message, iout ) bind ( C, name="print_line_c2f" )
! Interface from C to FORTRAN to call mqc_error_i
! with 0 of the optional arguments
! Insure that All I/O to unit 6 happens through Fortran

        use iso_c_binding
        use iso_fortran_env
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
        use iso_fortran_env
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
      use iso_fortran_env
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
    subroutine mqc_get_MatrixFile_Name_F2C(FileName_F, iout_F)
      use iso_c_binding
      use iso_fortran_env
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
      use iso_fortran_env
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
      use iso_fortran_env
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
      use iso_fortran_env
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

