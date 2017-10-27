      Module MQC_General
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
!           (1)  MQC suite control;
!           (2)  Printing;
!           (3)  Character conversion and manipulation;
!           (4)  Algebra; and
!           (5)  Other.
!
!
!
!
!     Define Types and Classes.
!
      implicit none
      real,public,parameter::angPBohr=0.52917706d0                                                          ! Angstom per Bohr
      real,public,parameter::kGPAMU=1.660538921d-27                                                         ! Kilograms per atomic mass unit
      real,public,parameter::planck=6.62606957d-34                                                          ! Planck constant
      real,public,parameter::avogadro=6.02214129d+23                                                        ! Avogadro's number
      real,public,parameter::jPCal=4.184d00                                                                 ! Joules per calorie
      real,public,parameter::jPHartree=4.35974434d-18                                                       ! Joules per Hartree
      real,public,parameter::sLight=2.99792458d10                                                           ! Speed of light
      real,public,parameter::boltzman=1.3806488d-23                                                         ! Boltzman Constant
      real,public,parameter::molVol=22.4139679d-3                                                           ! Molar volume of ideal gass in m**3 at 273.15 K 
      real,public,parameter::eMM=928.476430d-26                                                             ! Electron magnetic moment (J T**-1) (sign flipped).
      real,public,parameter::pRM=1.672621777d-27                                                            ! Proton rest mass (Kg)
      real,public,parameter::gFree=2.00231930436153d0                                                       ! Free-electron g-factor
      real,public,parameter::ESUPElec = 1.602176565d-19 * sLight / 10.0                                     ! Electrostatic units (ESU) per electron
      real,public,parameter::mPBohr = 0.52917706d0/1.0d10                                                   ! Meters per Bohr
      real,public,parameter::fineSC = 1.0d0/137.035999074d0                                                 ! Fine structure constant
      real,public,parameter::eMPKG = 4.35974434d-18*10000.0d0/(2.99792458d10*1.0d0/137.035999074d0)**2      ! Electron mass in Kg
      real,public,parameter::pi = 4.0d0 * ATan(1.0)                                                         ! Pi
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!
      interface mqc_print
        module procedure MQC_Print_Scalar_Integer
        module procedure MQC_Print_Scalar_Real
        module procedure MQC_Print_Vector_Array_Integer
        module procedure MQC_Print_Vector_Array_Real
        module procedure MQC_Print_Matrix_Array_Integer
        module procedure MQC_Print_Matrix_Array_Real
      end interface
!
      interface mqc_print_scalar
        module procedure MQC_Print_Scalar_Integer
        module procedure MQC_Print_Scalar_Real
      end interface
!
      interface mqc_print_vector
        module procedure MQC_Print_Vector_Array_Integer
        module procedure MQC_Print_Vector_Array_Real
      end interface
!
      interface mqc_print_matrix
        module procedure MQC_Print_Matrix_Array_Integer
        module procedure MQC_Print_Matrix_Array_Real
      end interface
!
      interface mqc_packedDiagonalMatrix2FullMatrix
        module procedure mqc_packedDiagonalMatrix2FullMatrix_integer
        module procedure mqc_packedDiagonalMatrix2FullMatrix_real
      end interface
!
!
!
!     Subroutines/Functions...
!
      CONTAINS
!
!
!----------------------------------------------------------------
!                                                               |
!     MQC Suite Control and System Interaction                  |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE mqc_error
      Subroutine mqc_error(Message,IOut)
!
!     This subroutine is used to kill a MQC job. The character string in
!     Message is printed. If possible, this routine will return an exit
!     value code of 999 to the operating system.
!
!     IOut is an OPTIONAL argument with the unit number corresponding to
!     the file where the error message should be written. If IOut is not
!     sent, unit number 6 is used.
!
!     -H. P. Hratchian, 2016
!
!
      implicit none
      character(LEN=*),intent(in),OPTIONAL::Message
      integer,intent(in),OPTIONAL::IOut
      integer::IJunk
!
 1000 Format(/,1x,'MQC ERROR: ',A,/)
!
      If(Present(Message)) then
        If(Present(IOut)) then
          Write(IOut,1000) TRIM(Message)
        else
          Write(6,1000) TRIM(Message)
        endIf
      endIf
      Stop 999
!
      Return
      End Subroutine MQC_Error


!
!PROCEDURE mqc_get_command_argument
      Subroutine mqc_get_command_argument(argNum,argument)
!
!     This subroutine is used to dynamically allocate command line arguments.
!     <argument> should be passed as deferred length allocatable character.  
!     
!     -L. M. Thompson, 2017.
!
!
      implicit none
      character(len=:),allocatable::argument
      integer,intent(in)::argNum
      integer::argLen
!
!     Do the work...
!
      call get_command_argument(argNum,length=argLen)
      allocate(character(len=argLen)::argument)
      call get_command_argument(argNum,argument)
!
      Return
      End Subroutine mqc_get_command_argument

!
!----------------------------------------------------------------
!                                                               |
!     Printing                                                  |
!                                                               |
!----------------------------------------------------------------
!
!
!PROCEDURE MQC_Print_Scalar_Integer
      Subroutine MQC_Print_Scalar_Integer(iOut,scalar,Header,Blank_At_Top, &
        Blank_At_Bottom,formatString)
!
!     This subroutine is used to print an integer scalar.
!
!     L. M. Thompson, 2016.
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Integer,Intent(In)::iOut
      Integer,Intent(In)::scalar
      Character(Len=*),Intent(In),optional::Header,formatString
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      character(len=256)::lineToPrint
!
      if(PRESENT(Blank_At_Top)) then
        if(Blank_At_Top) write(iOut,*)
      endIf
      if(PRESENT(header)) then
        lineToPrint = TRIM(header)
      else
        lineToPrint = ' '
      endIf
      if(PRESENT(formatString)) then
        lineToPrint = TRIM(lineToPrint)//' '//TRIM(integer2character(scalar,formatString))
      else
        lineToPrint = TRIM(lineToPrint)//' '//TRIM(integer2character(scalar))
      endIf
      lineToPrint = ADJUSTL(lineToPrint)
      write(iOut,'(A)') TRIM(lineToPrint)
      if(PRESENT(Blank_At_Bottom)) then
        if(Blank_At_Bottom) write(iOut,*)
      endIf
!
      return
      end subroutine MQC_Print_Scalar_Integer


!
!PROCEDURE MQC_Print_Scalar_Real
      Subroutine MQC_Print_Scalar_Real(iOut,scalar,Header,Blank_At_Top, &
        Blank_At_Bottom,formatString)
!
!     This subroutine is used to print an real scalar.
!
!     L. M. Thompson, 2016.
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Integer,Intent(In)::iOut
      Real,Intent(In)::scalar
      Character(Len=*),Intent(In),optional::Header,formatString
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      character(len=256)::lineToPrint
!
      if(PRESENT(Blank_At_Top)) then
        if(Blank_At_Top) write(iOut,*)
      endIf
      if(PRESENT(header)) then
        lineToPrint = TRIM(header)
      else
        lineToPrint = ' '
      endIf
      if(PRESENT(formatString)) then
        lineToPrint = TRIM(lineToPrint)//' '//TRIM(real2character(scalar,formatString))
      else
        lineToPrint = TRIM(lineToPrint)//' '//TRIM(real2character(scalar))
      endIf
      lineToPrint = ADJUSTL(lineToPrint)
      write(iOut,'(A)') TRIM(lineToPrint)
      if(PRESENT(Blank_At_Bottom)) then
        if(Blank_At_Bottom) write(iOut,*)
      endIf
!
      return
      end subroutine MQC_Print_Scalar_Real


!
!PROCEDURE MQC_Print_Matrix_Array_Integer
      Subroutine MQC_Print_Matrix_Array_Integer(IOut,Array,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print an integer Matrix array variable.
!
!     L. M. Thompson, 2016.
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Integer,Intent(In)::IOut
      Integer,dimension(:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer::I,J,NCols,NRows,IFirst,ILast
      Integer,Parameter::ColWidth=10 ! Xianghai: can we make it bigger?
!
!
 1000 Format(1x,A)
 1001 Format(5x,10(7x,I7))
 2001 Format(1x,I7,10I14)

      NRows = Size(Array,1)
      NCols = Size(Array,2)

      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,*)
      EndIf
      Write(IOut,1000) TRIM(Header)

      Do IFirst = 1,NCols,ColWidth
        ILast = Min(IFirst+ColWidth-1,NCols)
        Write(IOut,1001) (I,I=IFirst,ILast)
        Do I = 1,NRows
          Write(IOut,2001) I, (Array(I,J),J=IFirst,ILast)
        EndDo
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,*)
      EndIf
!
      Return
    End Subroutine MQC_Print_Matrix_Array_Integer


!PROCEDURE MQC_Print_Matrix_Array_Real
      Subroutine MQC_Print_Matrix_Array_Real(IOut,Array,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a real Matrix array variable.
!
!     L. M. Thompson, 2016.
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Integer,Intent(In)::IOut
      Real,dimension(:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer::I,J,NCols,NRows,IFirst,ILast
      Integer,Parameter::ColWidth=10 ! Xianghai: can we make it bigger?
!
 1000 Format(1x,A)
 1001 Format(5x,10(7x,I7))
 2002 Format(1x,I7,10F14.6)
!
      NRows = Size(Array,1)
      NCols = Size(Array,2)
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,*)
      EndIf
      Write(IOut,1000) TRIM(Header)

      Do IFirst = 1,NCols,ColWidth
        ILast = Min(IFirst+ColWidth-1,NCols)
        Write(IOut,1001) (I,I=IFirst,ILast)
        do I = 1,NRows
          Write(IOut,2002) I, (Array(I,J),J=IFirst,ILast)
        enddo
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,*)
      EndIf
!
      Return
      End Subroutine MQC_Print_Matrix_Array_Real


!PROCEDURE MQC_Print_Vector_Array_Integer
      subroutine MQC_Print_Vector_Array_Integer(IOut,Vector,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a Vector array type variable.
!
!     L. M. Thompson, 2016.
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      integer,Intent(In)::IOut
      integer,Dimension(:),Intent(In)::Vector
      character(Len=*),Intent(In),Optional::Header
      logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      integer::I,Length
!
 1000 Format(1x,A)
 1001 Format(1x,I7,2x,I14)
 1002 Format(1x,I7,2X,F14.6)
 1003 Format(1x,I7,2x,A)
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,*)
      EndIf
      if(Present(Header)) Write(IOut,1000) TRIM(Header)
      Length = Size(Vector)
      Do I = 1, Length
        Write(IOut,1001) I, Vector(I)
      EndDo
      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,*)
      EndIf
!
      Return
      End Subroutine MQC_Print_Vector_Array_Integer


!PROCEDURE MQC_Print_Vector_Array_Real
      Subroutine MQC_Print_Vector_Array_Real(IOut,Vector,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a Vector array type variable.
!
!     L. M. Thompson, 2016.
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Integer,Intent(In)::IOut
      Real,Dimension(:),Intent(In)::Vector
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer::I,Length
!
 1000 Format(1x,A)
 1002 Format(1x,I7,2X,F14.6)
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,*)
      EndIf
      Write(IOut,1000) TRIM(Header)
      Length = Size(Vector)
      Do I = 1, Length
        Write(IOut,1002) I, Vector(I)
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,*)
      EndIf
!
      Return
      End Subroutine MQC_Print_Vector_Array_Real


!
!----------------------------------------------------------------
!                                                               |
!     Character Conversion and Manipulation                     |
!                                                               |
!----------------------------------------------------------------
!
!
!PROCEDURE String_Change_Case
      Subroutine String_Change_Case(string,upperlower,stringOut)
!
!     This subroutine is used to change the case of a string, sent as in/out
!     dummy argument <string>. The input dummy argument <upperlower> is sent as
!     'U' or 'L' to indicate whether the routine should make <string> all upper-
!     case or lower-case. The output dummy argument <stringOut> is optional. If
!     <stringOut> is sent, then the modified version of <string> is returned in
!     <stringOut>. If <stringOut> is NOT sent, then <string> is replaced and
!     returned.
!     
!
!     -H. P. Hratchian, 2017
!
!
      implicit none
      character(len=*)::string
      character(len=1),intent(in)::upperlower
      character(len=*),optional::stringOut
!
      integer::i,charVal
!
!     Do the work...
!
      if(PRESENT(stringOut)) stringOut = ' '
      select case(upperlower)
      case('U','u')
        do i = 1,LEN(string)
          select case(string(i:i))
          case ('a':'z')
            charVal = ichar(string(i:i))
            if(PRESENT(stringOut)) then
              stringOut(i:i) = char(charVal-32)
            else
              string(i:i) = char(charVal-32)
            endIf
          case default
            if(PRESENT(stringOut)) stringOut(i:i) = string(i:i)
          endSelect
        endDo
      case('L','l')
        do i = 1,LEN(string)
          select case(string(i:i))
          case ('A':'Z')
            charVal = ichar(string(i:i))
            if(PRESENT(stringOut)) then
              stringOut(i:i) = char(charVal+32)
            else
              string(i:i) = char(charVal+32)
            endIf
          case default
            if(PRESENT(stringOut)) stringOut(i:i) = string(i:i)
          endSelect
        endDo
      case default
        call MQC_Error('Unknown upperLower in String_Chage_Case.')
      endSelect
!
      Return
      End Subroutine String_Change_Case


!PROCEDURE Integer2Character
      function integer2character(integerIn,formatString) result(integerString)
!
!     This function converts an integer variable and returns a character string.
!
!
!     -H. P. Hratchian, 2017.
!
!
      implicit none
      integer,intent(in)::integerIn
      character(*),intent(in),optional::formatString
      character(len=256)::myFormatString,integerString
!
      if(PRESENT(formatString)) then
        myFormatString = formatString
      else
        myFormatString = 'I10'
      endIf
      myFormatString = '('//TRIM(myFormatString)//')'
      write(integerString,myFormatString) integerIn
      integerString = ADJUSTL(integerString)
!
      return
      end function integer2character


!PROCEDURE real2character
      function real2character(realIn,formatString) result(realString)
!
!     This function converts a real variable and returns a character string. The
!     input dummy argument <formatString> is OPTIONAL and can be sent to specify
!     the desired formatting of the realIn. Note that <formatString> must
!     conform to standard fortran requirments.
!
!
!     -H. P. Hratchian, 2017.
!
!
      implicit none
      real,intent(in)::realIn
      character(*),intent(in),optional::formatString
      character(len=256)::myFormatString,realString
!
      if(PRESENT(formatString)) then
        myFormatString = formatString
      else
        myFormatString = 'f20.5'
      endIf
      myFormatString = '('//TRIM(myFormatString)//')'
      write(realString,myFormatString) realIn
      realString = ADJUSTL(realString)
!
      return
      end function real2character


!PROCEDURE Build_String_Add_Int
      Subroutine Build_String_Add_Int(IntIn,String,IntWidth,NDigits)
!
!     This subroutine is used to append an integer to a string.
!
!     -H. P. Hratchian, 2015
!
!
      Implicit None
      Integer,Intent(In)::IntIn
      Character(Len=*),Intent(InOut)::String
      Integer,Optional,Intent(In)::IntWidth,NDigits
      Integer::My_IntWidth
      Character(Len=512)::FormatString,TempChar

!     Set-up My_IntWidth and then fill FormatString.
!
      If(Present(NDigits).and..not.Present(IntWidth)) then
        My_IntWidth = NDigits
      else if(Present(IntWidth)) then
        My_IntWidth = IntWidth
      else
        My_IntWidth = 0
      endIf
      FormatString = '(I'
      If(My_IntWidth.gt.0) then
        Write(TempChar,'(I10)') My_IntWidth
        FormatString = TRIM(FormatString)//TRIM(ADJUSTL(TempChar))
      endIf
      If(Present(NDigits)) then
        Write(TempChar,'(I10)') NDigits
      endIf
      FormatString = TRIM(FormatString)//')'
      Write(TempChar,FormatString) IntIn
      String = TRIM(String)//TRIM(ADJUSTL(TempChar))
      Return
      End Subroutine Build_String_Add_Int


!
!----------------------------------------------------------------
!                                                               |
!     Algebra                                                   |
!                                                               |
!----------------------------------------------------------------
!
!
!PROCEDURE mqc_packedDiagonalMatrix2FullMatrix_integer
      Subroutine mqc_packedDiagonalMatrix2FullMatrix_integer(matrixDiagonal,  &
        matrixFull)
!
!     This subroutine accepts an input one-dimensional array of the diagonal
!     elements of a diagonal matrix (input dummy argument <matrixDiagonal>) and
!     outputs a full matrix (output dummy argument <matrixFull>). Note that
!     matrixFull is initialized in this routine. Furthermore, <matrixFull> must
!     be allocated prior to calling this routine.
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      integer,dimension(:),intent(in)::matrixDiagonal
      integer,dimension(:,:),intent(out)::matrixFull
      integer::nDim,i
!
      nDim = SIZE(matrixDiagonal)
      if(nDim.ne.SIZE(matrixFull,1).or.nDim.ne.SIZE(matrixFull,2))  &
        call mqc_error('mqc_packedDiagonalMatrix2FullMatrix_real: disagreement in input/output array size!')
      matrixFull = 0.0
      do i = 1,nDim
        matrixFull(i,i) = matrixDiagonal(i)
      endDo
!
      return
      end Subroutine mqc_packedDiagonalMatrix2FullMatrix_integer


!
!PROCEDURE mqc_packedDiagonalMatrix2FullMatrix_real
      Subroutine mqc_packedDiagonalMatrix2FullMatrix_real(matrixDiagonal,  &
        matrixFull)
!
!     This subroutine accepts an input one-dimensional array of the diagonal
!     elements of a diagonal matrix (input dummy argument <matrixDiagonal>) and
!     outputs a full matrix (output dummy argument <matrixFull>). Note that
!     matrixFull is initialized in this routine. Furthermore, <matrixFull> must
!     be allocated prior to calling this routine.
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      real,dimension(:),intent(in)::matrixDiagonal
      real,dimension(:,:),intent(out)::matrixFull
      integer::nDim,i
!
      nDim = SIZE(matrixDiagonal)
      if(nDim.ne.SIZE(matrixFull,1).or.nDim.ne.SIZE(matrixFull,2))  &
        call mqc_error('mqc_packedDiagonalMatrix2FullMatrix_real: disagreement in input/output array size!')
      matrixFull = 0.0
      do i = 1,nDim
        matrixFull(i,i) = matrixDiagonal(i)
      endDo
!
      return
      end Subroutine mqc_packedDiagonalMatrix2FullMatrix_real


!
!
!----------------------------------------------------------------
!                                                               |
!     Other                                                     |
!                                                               |
!----------------------------------------------------------------
!

!
!
      End Module MQC_General
