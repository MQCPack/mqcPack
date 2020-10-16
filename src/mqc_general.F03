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
      use iso_fortran_env
      implicit none
      real(kind=real64),public,parameter::angPBohr=0.52917706d0                                                          ! Angstom per Bohr
      real(kind=real64),public,parameter::kGPAMU=1.660538921d-27                                                         ! Kilograms per atomic mass unit
      real(kind=real64),public,parameter::planck=6.62606957d-34                                                          ! Planck constant
      real(kind=real64),public,parameter::avogadro=6.02214129d+23                                                        ! Avogadro's number
      real(kind=real64),public,parameter::jPCal=4.184d00                                                                 ! Joules per calorie
      real(kind=real64),public,parameter::jPHartree=4.35974434d-18                                                       ! Joules per Hartree
      real(kind=real64),public,parameter::sLight=2.99792458d10                                                           ! Speed of light
      real(kind=real64),public,parameter::boltzman=1.3806488d-23                                                         ! Boltzman Constant
      real(kind=real64),public,parameter::molVol=22.4139679d-3                                                           ! Molar volume of ideal gass in m**3 at 273.15 K 
      real(kind=real64),public,parameter::eMM=928.476430d-26                                                             ! Electron magnetic moment (J T**-1) (sign flipped).
      real(kind=real64),public,parameter::pRM=1.672621777d-27                                                            ! Proton rest mass (Kg)
      real(kind=real64),public,parameter::gFree=2.00231930436153d0                                                       ! Free-electron g-factor
      real(kind=real64),public,parameter::ESUPElec = 1.602176565d-19 * sLight / 10.0                                     ! Electrostatic units (ESU) per electron
      real(kind=real64),public,parameter::mPBohr = 0.52917706d0/1.0d10                                                   ! Meters per Bohr
      real(kind=real64),public,parameter::fineSC = 1.0d0/137.035999074d0                                                 ! Fine structure constant
      real(kind=real64),public,parameter::eMPKG = 4.35974434d-18*10000.0d0/(2.99792458d10*1.0d0/137.035999074d0)**2      ! Electron mass in Kg
      real(kind=real64),public,parameter::pi = 4.0d0 * ATan(1.0)                                                         ! Pi
      real(kind=real64),public,parameter::MQC_small = 1.0d-6
                     ! 10^-6, used as a threshold for near zero.
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
        module procedure MQC_Print_Vector_Array_Complex
        module procedure MQC_Print_Matrix_Array_Integer
        module procedure MQC_Print_Matrix_Array_Real
        module procedure MQC_Print_Matrix_Array_Complex
        module procedure MQC_Print_Rank4Tensor_Array_Real
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
        module procedure MQC_Print_Vector_Array_Complex
      end interface
!
      interface mqc_print_matrix
        module procedure MQC_Print_Matrix_Array_Integer
        module procedure MQC_Print_Matrix_Array_Real
        module procedure MQC_Print_Matrix_Array_Complex
      end interface
!
      interface mqc_print_rank4Tensor
        module procedure MQC_Print_Rank4Tensor_Array_Real
      end interface
!
      interface mqc_packedDiagonalMatrix2FullMatrix
        module procedure mqc_packedDiagonalMatrix2FullMatrix_integer
        module procedure mqc_packedDiagonalMatrix2FullMatrix_real
      end interface
!
      interface mqc_matrixSymm2Full
        module procedure mqc_matrixSymm2Full_integer
        module procedure mqc_matrixSymm2Full_real
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
!PROCEDURE mqc_version_check
      function mqc_version_check(isMajor,isMinor,isRevision,  &
        olderThanMajor,olderThanMinor,olderThanRevision,  &
        newerThanMajor,newerThanMinor,newerThanRevision) result(ok)
!
!     This function is used to test if the MQC package being used is new and/or
!     older than a specified version. The output argument is a logical, which is
!     returned TRUE or FALSE for the version limitations requested being met or
!     not. The dummy arguments <isMajor>, <isMinor>, <isRevision>,
!     <olderThanMajor>, <olderThanMinor>, <olderThanVevision>, <newerThanMajor>,
!     <newerThanMinor>, and <newerThanVevision> are all OPTIONAL and should be
!     sent with appropriate integers. The "olderThan" and "newerThan" arguments
!     are treated as being inclusive. Asking for a major revision with
!     oldThanMajor=2019 will be TRUE is the Major version is <= 2019, etc.
!
!
      implicit none
      integer,OPTIONAL,intent(in)::isMajor,isMinor,isRevision
      integer,OPTIONAL,intent(in)::olderThanMajor,olderThanMinor,olderThanRevision
      integer,OPTIONAL,intent(in)::newerThanMajor,newerThanMinor,newerThanRevision
      logical::ok
      integer::major,minor,revision
!
      ok = .true.
      call mqc_version(major,minor,revision)
!
      if(PRESENT(isMajor))           ok = ok.and.major.eq.isMajor
      if(PRESENT(isMinor))           ok = ok.and.minor.eq.isMinor
      if(PRESENT(isRevision))        ok = ok.and.revision.eq.isRevision
      if(PRESENT(olderThanMajor))    ok = ok.and.major.le.olderThanMajor
      if(PRESENT(olderThanMinor))    ok = ok.and.minor.le.olderThanMinor
      if(PRESENT(olderThanRevision)) ok = ok.and.revision.le.olderThanRevision
      if(PRESENT(newerThanMajor))    ok = ok.and.major.ge.newerThanMajor
      if(PRESENT(newerThanMinor))    ok = ok.and.minor.ge.newerThanMinor
      if(PRESENT(newerThanRevision)) ok = ok.and.revision.ge.newerThanRevision
!
      return
      end function mqc_version_check


!
!PROCEDURE mqc_version
      subroutine mqc_version(major,minor,revision,versionString)
!
!     This subroutine returns the current MQC version number. The full MQC
!     version number is comprised of three components: major.minor.revision. The
!     major version indicated is given by the last two digits of the versioning
!     year; the minor label is two digits indicating a versioning month; and the
!     last digit gives a revision number.
!
!
      implicit none
      integer,OPTIONAL,intent(out)::major,minor,revision
      character(len=*),OPTIONAL,intent(out)::versionString
!
      if(PRESENT(major))         major         = 20
      if(PRESENT(minor))         minor         = 5
      if(PRESENT(revision))      revision      = 1
      if(PRESENT(versionString)) versionString = '20.05.1'
!
      return
      end subroutine mqc_version


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
      integer(kind=int64),intent(in),OPTIONAL::IOut
      integer(kind=int64)::myIOut,IJunk
!
 1000 Format(/,1x,'MQC ERROR: ',A,/)
!
      myIOut = 6
      If(Present(IOut)) myIOut = IOut
      If(Present(Message)) Write(myIOut,1000) TRIM(Message)
      flush(myIOut)
      call mqc_abort()
!
      Return
      End Subroutine MQC_Error

      subroutine mqc_error_i( Message, IOut, Var_Name1, INT1, &
           Var_Name2, INT2, Var_Name3, INT3, Var_Name4, INT4, &
           Var_Name5, INT5, Var_Name6, INT6)
!
!     This subroutine is used to kill a MQC job. The character string in
!     Message is printed as are the names and values of variables related 
!     to the failure. If possible, this routine will abort, which should 
!     cause a stack trace.
!
!     The names and values of the integer variables are OPTIONAL arguments, 
!     so that one routine can be used for many types of cases
!
!     When modifing code use us this routine, if IOUT is not specified, use 
!     unit number 6.
!
!     -H. P. Hratchian, 2016
!
!
      implicit none
      character(LEN=*),intent(in)::Message
      integer(kind=int64),intent(in)::IOut
      character(LEN=*),intent(in),OPTIONAL::Var_Name1
      character(LEN=*),intent(in),OPTIONAL::Var_Name2
      character(LEN=*),intent(in),OPTIONAL::Var_Name3
      character(LEN=*),intent(in),OPTIONAL::Var_Name4
      character(LEN=*),intent(in),OPTIONAL::Var_Name5
      character(LEN=*),intent(in),OPTIONAL::Var_Name6
      integer(kind=int64),intent(in),OPTIONAL::INT1
      integer(kind=int64),intent(in),OPTIONAL::INT2
      integer(kind=int64),intent(in),OPTIONAL::INT3
      integer(kind=int64),intent(in),OPTIONAL::INT4
      integer(kind=int64),intent(in),OPTIONAL::INT5
      integer(kind=int64),intent(in),OPTIONAL::INT6
      integer(kind=int64)::IJunk
!
 1000 Format(4x,A,' = ',I15)
 1030 Format(1x,'MQC ERROR: ',A)
!
      Write(IOut,1030) TRIM(Message)
      If(Present(INT1)) then
         Write(IOut,1000) TRIM(Var_Name1), INT1
      endIf
      If(Present(INT2)) then
         Write(IOut,1000) TRIM(Var_Name2), INT2
      endIf
      If(Present(INT3)) then
         Write(IOut,1000) TRIM(Var_Name3), INT3
      endIf
      If(Present(INT4)) then
         Write(IOut,1000) TRIM(Var_Name4), INT4
      endIf
      If(Present(INT5)) then
         Write(IOut,1000) TRIM(Var_Name5), INT5
      endIf
      If(Present(INT6)) then
         Write(IOut,1000) TRIM(Var_Name6), INT6
      endIf
!
      flush(IOut)
      call mqc_abort()

      Return
      end Subroutine mqc_error_i
      !

      subroutine mqc_error_r( Message, IOut, Var_NAME1, FLT1, &
           Var_NAME2, FLT2, Var_NAME3, FLT3, Var_NAME4, FLT4, &
           Var_NAME5, FLT5, Var_NAME6, FLT6)
!
!     This subroutine is used to kill a MQC job. The character string in
!     Message is printed as are the names and values of variables related 
!     to the failure. If possible, this routine will abort, which should 
!     cause a stack trace.
!
!     The names and values of the integer variables are OPTIONAL arguments, 
!     so that one routine can be used for many types of cases
!
!     When modifing code use us this routine, if IOUT is not specified, use 
!     unit number 6.
!
!     -H. P. Hratchian, 2016
!
!
      implicit none
      character(LEN=*),intent(in)::Message
      integer(kind=int64),intent(in)::IOut
      character(LEN=*),intent(in),OPTIONAL::Var_NAME1
      character(LEN=*),intent(in),OPTIONAL::Var_NAME2
      character(LEN=*),intent(in),OPTIONAL::Var_NAME3
      character(LEN=*),intent(in),OPTIONAL::Var_NAME4
      character(LEN=*),intent(in),OPTIONAL::Var_NAME5
      character(LEN=*),intent(in),OPTIONAL::Var_NAME6
      real(kind=real64),intent(in),OPTIONAL::FLT1
      real(kind=real64),intent(in),OPTIONAL::FLT2
      real(kind=real64),intent(in),OPTIONAL::FLT3
      real(kind=real64),intent(in),OPTIONAL::FLT4
      real(kind=real64),intent(in),OPTIONAL::FLT5
      real(kind=real64),intent(in),OPTIONAL::FLT6
      integer(kind=int64)::IJunk
!
 1000 Format(4x,A,' = ',E20.14)
 1030 Format(1x,'MQC ERROR: ',A)
!
      Write(IOut,1030) TRIM(Message)

      If(Present(FLT1)) then
         Write(IOut,1000) TRIM(Var_NAME1), FLT1
      endIf
      If(Present(FLT2)) then
         Write(IOut,1000) TRIM(Var_NAME2), FLT2
      endIf
      If(Present(FLT3)) then
         Write(IOut,1000) TRIM(Var_NAME3), FLT3
      endIf
      If(Present(FLT4)) then
         Write(IOut,1000) TRIM(Var_NAME4), FLT4
      endIf
      If(Present(FLT5)) then
         Write(IOut,1000) TRIM(Var_NAME5), FLT5
      endIf
      If(Present(FLT6)) then
         Write(IOut,1000) TRIM(Var_NAME6), FLT6
      endIf
!
      flush(IOut)
      call mqc_abort()

      Return
      end Subroutine mqc_error_r
      !

      subroutine mqc_error_l( Message, IOut, Var_NAME1, LOG1, &
           Var_NAME2, LOG2, Var_NAME3, LOG3, Var_NAME4, LOG4, &
           Var_NAME5, LOG5, Var_NAME6, LOG6)
!
!     This subroutine is used to kill a MQC job. The character string in
!     Message is printed as are the names and values of variables related 
!     to the failure. If possible, this routine will abort, which should 
!     cause a stack trace.
!
!     The names and values of the integer variables are OPTIONAL arguments, 
!     so that one routine can be used for many types of cases
!
!     When modifing code use us this routine, if IOUT is not specified, use 
!     unit number 6.
!
!     -H. P. Hratchian, 2016
!
!
      implicit none
      character(LEN=*),intent(in)::Message
      integer(kind=int64),intent(in)::IOut
      character(LEN=*),intent(in),OPTIONAL::Var_NAME1
      character(LEN=*),intent(in),OPTIONAL::Var_NAME2
      character(LEN=*),intent(in),OPTIONAL::Var_NAME3
      character(LEN=*),intent(in),OPTIONAL::Var_NAME4
      character(LEN=*),intent(in),OPTIONAL::Var_NAME5
      character(LEN=*),intent(in),OPTIONAL::Var_NAME6
      logical,intent(in),OPTIONAL::LOG1
      logical,intent(in),OPTIONAL::LOG2
      logical,intent(in),OPTIONAL::LOG3
      logical,intent(in),OPTIONAL::LOG4
      logical,intent(in),OPTIONAL::LOG5
      logical,intent(in),OPTIONAL::LOG6
      integer(kind=int64)::IJunk
!
 1030 Format(1x,'MQC ERROR: ',A)
 1050 Format(4x,A,' = ','.TRUE.')
 1060 Format(4x,A,' = ','.FALSE.')
!
      Write(IOut,1030) TRIM(Message)
      If(Present(LOG1)) then
         If(LOG1) then
            Write(IOut,1050) TRIM(Var_NAME1)
         else 
            Write(IOut,1060) TRIM(Var_NAME1)
         endif
      endif
      If(Present(LOG2)) then
         If(LOG2) then
            Write(IOut,1050) TRIM(Var_NAME2)
         else 
            Write(IOut,1060) TRIM(Var_NAME2)
         endif
      endif
      If(Present(LOG3)) then
         If(LOG3) then
            Write(IOut,1050) TRIM(Var_NAME3)
         else 
            Write(IOut,1060) TRIM(Var_NAME3)
         endif
      endif
      If(Present(LOG4)) then
         If(LOG4) then
            Write(IOut,1050) TRIM(Var_NAME4)
         else 
            Write(IOut,1060) TRIM(Var_NAME4)
         endif
      endif

      If(Present(LOG5)) then
         If(LOG5) then
            Write(IOut,1050) TRIM(Var_NAME5)
         else 
            Write(IOut,1060) TRIM(Var_NAME5)
         endif
      endif

      If(Present(LOG6)) then
         If(LOG6) then
            Write(IOut,1050) TRIM(Var_NAME6)
         else 
            Write(IOut,1060) TRIM(Var_NAME6)
         endif
      endif
!
      flush(IOut)
      call mqc_abort()

      Return
      end Subroutine mqc_error_l
      !

      subroutine mqc_error_A( Message, IOut, CHAR1, Var_NAME1, &
           CHAR2, Var_NAME2, CHAR3, Var_NAME3, CHAR4, Var_NAME4, &
           CHAR5, Var_NAME5, CHAR6, Var_NAME6)
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
      character(LEN=*),intent(in)::Message
      integer(kind=int64),intent(in)::IOut
      character(LEN=*),intent(in),OPTIONAL::CHAR1
      character(LEN=*),intent(in),OPTIONAL::Var_NAME1
      character(LEN=*),intent(in),OPTIONAL::CHAR2
      character(LEN=*),intent(in),OPTIONAL::Var_NAME2
      character(LEN=*),intent(in),OPTIONAL::CHAR3
      character(LEN=*),intent(in),OPTIONAL::Var_NAME3
      character(LEN=*),intent(in),OPTIONAL::CHAR4
      character(LEN=*),intent(in),OPTIONAL::Var_NAME4
      character(LEN=*),intent(in),OPTIONAL::CHAR5
      character(LEN=*),intent(in),OPTIONAL::Var_NAME5
      character(LEN=*),intent(in),OPTIONAL::CHAR6
      character(LEN=*),intent(in),OPTIONAL::Var_NAME6
      integer(kind=int64)::IJunk
!
 1000 Format(4x,A,' = ',A)
 1030 Format(1x,'MQC ERROR: ',A)
!

      Write(IOut,1030) TRIM(Message)
      If(Present(Var_NAME1)) then
         Write(IOut,1000) TRIM(CHAR1), TRIM(Var_NAME1)
      endif
      If(Present(Var_NAME2)) then
         Write(IOut,1000) TRIM(CHAR2), TRIM(Var_NAME2)
      endif
      If(Present(Var_NAME3)) then
         Write(IOut,1000) TRIM(CHAR3), TRIM(Var_NAME3)
      endif
      If(Present(Var_NAME4)) then
         Write(IOut,1000) TRIM(CHAR4), TRIM(Var_NAME4)
      endif
      If(Present(Var_NAME5)) then
         Write(IOut,1000) TRIM(CHAR5), TRIM(Var_NAME5)
      endif
      If(Present(Var_NAME6)) then
         Write(IOut,1000) TRIM(CHAR6), TRIM(Var_NAME6)
      endif
!

      flush(IOut)
      call mqc_abort()

      Return
      end Subroutine mqc_error_A

      
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
      integer(kind=int64),intent(in)::argNum
      integer(kind=int64)::argLen
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
      Integer(kind=int64),Intent(In)::iOut
      Integer(kind=int64),Intent(In)::scalar
      Character(Len=*),Intent(In),optional::Header,formatString
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      character(len=256)::lineToPrint
!
 1020 Format(" ")
      if(PRESENT(Blank_At_Top)) then
        if(Blank_At_Top) write(iOut,1020)
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
        if(Blank_At_Bottom) write(iOut,1020)
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
      Integer(kind=int64),Intent(In)::iOut
      Real(kind=real64),Intent(In)::scalar
      Character(Len=*),Intent(In),optional::Header,formatString
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      character(len=256)::lineToPrint
!
 1020 Format(" ")
      if(PRESENT(Blank_At_Top)) then
        if(Blank_At_Top) write(iOut,1020)
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
        if(Blank_At_Bottom) write(iOut,1020)
      endIf
!
      return
      end subroutine MQC_Print_Scalar_Real
!
!
!PROCEDURE MQC_Print_Rank4Tensor_Array_Real
      Subroutine MQC_Print_Rank4Tensor_Array_Real(IOut,Array,Header,  &
        Blank_At_Top,Blank_At_Bottom)
!
!     This subroutine is used to print a real rank-4 tensor array.
!
!     H. P. Hratchian, 2020.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Real(kind=real64),dimension(:,:,:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::NI,NJ,NK,NL,i,j,k,l
!
 1000 Format(/,1x,A)
 1001 Format(5x,10(7x,I7))
 1020 Format(" ")
 2002 Format(1x,'I=',I4,' J=',I4,' K=',I4,' L=',I4,3x,10F14.6)
!
      NI = Size(Array,1)
      NJ = Size(Array,2)
      NK = Size(Array,3)
      NL = Size(Array,4)
!
      if(PRESENT(Blank_At_Top)) then
        if(Blank_At_Top) write(IOut,1020)
      endIf
      if(PRESENT(Header)) write(IOut,1000) TRIM(Header)
      do i = 1,NI
        do j = 1,NJ
          do k = 1,NK
            do l = 1,NL
              write(iOut,2002) i,j,k,l,Array(i,j,k,l)
            endDo
          endDo
        endDo
      endDo
      if(PRESENT(Blank_At_Bottom)) then
        if(Blank_At_Bottom) write(IOut,1020)
      endIf
!
      return
      end subroutine MQC_Print_Rank4Tensor_Array_Real

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
      Integer(kind=int64),Intent(In)::IOut
      Integer(kind=int64),dimension(:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64),Parameter::ColWidth=10
!
!
 1000 Format(/,1x,A)
 1001 Format(5x,10(7x,I7))
 1020 Format(" ")
 2001 Format(1x,I7,10I14)

      NRows = Size(Array,1)
      NCols = Size(Array,2)

      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
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
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
    End Subroutine MQC_Print_Matrix_Array_Integer
!
!
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
      Integer(kind=int64),Intent(In)::IOut
      Real(kind=real64),dimension(:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64),Parameter::ColWidth=10
!
 1000 Format(/,1x,A)
 1001 Format(5x,10(7x,I7))
 1020 Format(" ")
 2002 Format(1x,I7,10F14.6)
!
      NRows = Size(Array,1)
      NCols = Size(Array,2)
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
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
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Matrix_Array_Real

!
!PROCEDURE MQC_Print_Matrix_Array_Complex
      Subroutine MQC_Print_Matrix_Array_Complex(Array,IOut,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a complex Matrix array variable.
!
!     L. M. Thompson, 2019.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Complex(kind=real64),dimension(:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64),Parameter::ColWidth=10 
!
 1000 Format(/,1x,A)
 1001 Format(5x,10(7x,I7))
 1020 Format(" ")
 2003 Format(1x,I7,10(F12.5,F11.5,"i"))
!
      NRows = Size(Array,1)
      NCols = Size(Array,2)
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      Write(IOut,1000) TRIM(Header)

      Do IFirst = 1,NCols,ColWidth
        ILast = Min(IFirst+ColWidth-1,NCols)
        Write(IOut,1001) (I,I=IFirst,ILast)
        do I = 1,NRows
          Write(IOut,2003) I, (Array(I,J),J=IFirst,ILast)
        enddo
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Matrix_Array_Complex


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
      integer(kind=int64),Intent(In)::IOut
      integer(kind=int64),Dimension(:),Intent(In)::Vector
      character(Len=*),Intent(In),Optional::Header
      logical(kind=int64),Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      integer(kind=int64)::I,Length
!
 1000 Format(/,1x,A)
 1001 Format(1x,I7,2x,I14)
 1002 Format(1x,I7,2X,F14.6)
 1003 Format(1x,I7,2x,A)
 1020 Format(' ')
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      if(Present(Header)) Write(IOut,1000) TRIM(Header)
      Length = Size(Vector)
      Do I = 1, Length
        Write(IOut,1001) I, Vector(I)
      EndDo
      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
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
      Integer(kind=int64),Intent(In)::IOut
      Real(kind=real64),Dimension(:),Intent(In)::Vector
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,Length
!
 1000 Format(/,1x,A)
 1002 Format(1x,I7,2X,F14.6)
 1020 Format(" ")
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      Write(IOut,1000) TRIM(Header)
      Length = Size(Vector)
      Do I = 1, Length
        Write(IOut,1002) I, Vector(I)
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Vector_Array_Real


!PROCEDURE MQC_Print_Vector_Array_Complex
      Subroutine MQC_Print_Vector_Array_Complex(IOut,Vector,Header,Blank_At_Top, &
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
      Integer(kind=int64),Intent(In)::IOut
      Complex(kind=real64),Dimension(:),Intent(In)::Vector
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,Length
!
 1000 Format(/,1x,A)
 1002 Format(1x,I7,F12.5,F11.5,"i")
 1020 Format(" ")
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      Write(IOut,1000) TRIM(Header)
      Length = Size(Vector)
      Do I = 1, Length
        Write(IOut,1002) I, Vector(I)
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Vector_Array_Complex

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
      integer(kind=int64)::i,charVal
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
        call MQC_Error_A('Unknown upperLower in String_Chage_Case.', 6, &
             'upperlower', upperlower)
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
      integer(kind=int64),intent(in)::integerIn
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
      real(kind=real64),intent(in)::realIn
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
      Integer(kind=int64),Intent(In)::IntIn
      Character(Len=*),Intent(InOut)::String
      Integer(kind=int64),Optional,Intent(In)::IntWidth,NDigits
      Integer(kind=int64)::My_IntWidth
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
!PROCEDURE mqc_normalizeVector
      Subroutine mqc_normalizeVector(vector)
!
!     This subroutine accepts a vector and returns its normalized form.
!
!
!     H. P. Hratchian, 2019.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),intent(inOut)::vector
      real(kind=real64)::vectorMagnitude
!
      vectorMagnitude = SQRT(dot_product(vector,vector))
      if(vectorMagnitude.gt.MQC_small) then
        vector = vector/vectorMagnitude
      else
        vector = float(0)
      endIf
!
      return
      end subroutine mqc_normalizeVector
!
!
!PROCEDURE mqc_vectorTrimZero
      Subroutine mqc_vectorTrimZero(vectorIn,vectorOut,small)
!
!     This subroutine accepts a vector, <vectorIn>, and sets all entries that
!     have absolute values below a threshold to float(0).
!
!     If the optional dummy argument <vectorOut> is sent, then <vectorIn> is
!     left unchanges and the zero-d version of the vector is returned in
!     <vectorOut>. Otherwise, the zero-ing is done in place in <vectorIn>.
!
!     The optional dummy argument <small> may be sent with a user specified
!     threshold. Otherwise, the threshold is set to the MQC system value of
!     <MQC_small>.
!
!
!     H. P. Hratchian, 2020.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),intent(inOut)::vectorIn
      real(kind=real64),dimension(:),allocatable,optional,intent(Out)::vectorOut
      real(kind=real64),optional,,intent(in)::small
      integer(kind=int64)::i
      real(kind=real64)::mySmall,zero=0.
!
!
!     Initial set-up code...
!
      mySmall = MQC_small
      if(PRESENT(small)) mySmall = small
      if(PRESENT(vectorOut)) then
        if(Allocated(vectorOut)) then
          if(SIZE(vectorOut).ne.SIZE(vectorIn)) then
            DeAllocate(vectorOut)
            Allocate(vectorout(size(vectorin)))
          endIf
        else
          Allocate(vectorout(size(vectorin)))
        endIf
      endIf
!
!     Do trim work.
!
      if(PRESENT(vectorOut)) then
        vectorOut = vectorIn
        do i = 1,SIZE(vectorIn)
          if(abs(vectorOut(i)).lt.mySmall) vectorOut(i) = zero
        endDo
      else
        do i = 1,SIZE(vectorIn)
          if(abs(vectorIn(i)).lt.mySmall) vectorIn(i) = zero
        endDo
      endIf
!
      return
      end Subroutine mqc_vectorTrimZero
!
!
!PROCEDURE mqc_vectorPhase
      Subroutine mqc_vectorPhase(vector,largestPositive,trimZero)
!
!     This subroutine accepts a vector fixes the phase so that the largest
!     magniture value is either positive (largestPostive=.TRUE.) or negative
!     (largestPostive=.FALSE.). The vector phase is modified, as needed, in
!     place.
!
!     Dummy argument largestPostive is OPTIONAL and defaults to TRUE.
!
!     Dummy argument trimZero is OPTIONAL and defaults to TRUE.
!
!
!     H. P. Hratchian, 2019.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),intent(inOut)::vector
      logical,intent(in),optional::largestPositive,trimZero
      integer(kind=int64)::i,largestMagnitudePosition
      real(kind=real64)::largestMagnitudeValue
      logical::myLargestPositive,myTrimZero
!
!
!     Initial set-up code...
!
      myLargestPositive = .true.
      if(PRESENT(largestPositive)) myLargestPositive = largestPositive
      myTrimZero = .true.
      if(PRESENT(trimZero)) myTrimZero = trimZero
!
!     Do phase correction and trim work.
!
      if(myTrimZero.and.(abs(vector(1)).lt.MQC_small))  &
        vector(1) = float(0)
      largestMagnitudeValue = abs(vector(1)) + MQC_small
      largestMagnitudePosition = 1
      do i = 2,SIZE(vector)
        if(myTrimZero.and.(abs(vector(i)).lt.MQC_small))  &
          vector(i) = float(0)
        if(abs(vector(i)).gt.largestMagnitudeValue) then
          largestMagnitudeValue = abs(vector(i)) + MQC_small
          largestMagnitudePosition = i
        endIf
      endDo
      if(myLargestPositive) then
        if(vector(largestMagnitudePosition).lt.float(0)) vector = -vector
      else
        if(vector(largestMagnitudePosition).gt.float(0)) vector = -vector
      endIf
!
      return
      end subroutine mqc_vectorPhase


!PROCEDURE mqc_packedSymmetricMatrix2FullMatrix_integer
      Subroutine mqc_packedSymmetricMatrix2FullMatrix_integer(matrixSymmetric,  &
        matrixFull,upperLower)
!
!     This subroutine accepts an input one-dimensional array of a symmetric
!     matrix (input dummy argument <matrixSymmetric>) stored in packed form and
!     returns a full matrix. (output dummy argument <matrixFull>). Note that
!     <matrixFull> is initialized in this routine. Furthermore, <matrixFull>
!     must be allocated prior to calling this routine.
!
!     Optional input dummy argument <upperLower> is a character string with 'U'
!     or 'L' to denote whether <matrixSymmetric> is sent in an upper- or
!     lower-symmetric column-wise storage format. If not sent, the default value
!     for <upperLower> is 'U'.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),dimension(:),intent(in)::matrixSymmetric
      integer(kind=int64),dimension(:,:),intent(out)::matrixFull
      integer(kind=int64)::i,j,k,nDim,nSymm
      character(len=*),intent(in),optional::upperLower
      character(len=1)::myUpperLower
!
!
!     Error traps...
!
      if(SIZE(matrixFull,1).ne.SIZE(matrixFull,2))  &
        call mqc_error_i('mqc_packedSymmetricMatrix2FullMatrix_integer: matrixFull array is NOT square!',  &
        6,'nDim(1)',SIZE(matrixFull,1),'nDim(2)',SIZE(matrixFull,2))
      nDim = SIZE(matrixFull,1)
      nSymm = (nDim*(nDim+1))/2
      if(nSymm.ne.SIZE(matrixSymmetric,1))  &
        call mqc_error_i('mqc_packedSymmetricMatrix2FullMatrix_integer: matrixSymmetric has wrong dimension!',  &
        6,'nSymm',nSymm,'SIZE(matrixSymmetric)',SIZE(matrixSymmetric,1))
!
!     Set-up myUpperLower.
!
      if(PRESENT(upperLower)) then
        myUpperLower = TRIM(upperLower)
      else
        myUpperLower = 'U'
      endIf
!
!     Convert the symmetric matrix to full storage.
!
      matrixFull = 0
      select case(myUpperLower)
      case('U')
        k = 0
        do i = 1,nDim
          do j = 1,i
            k = k+1
            matrixFull(i,j) = matrixSymmetric(k)
            matrixFull(j,i) = matrixSymmetric(k)
          endDo
        endDo
      case('L')
        k = 0
        do i = 1,nDim
          do j = i,nDim
            k = k+1
            matrixFull(i,j) = matrixSymmetric(k)
            matrixFull(j,i) = matrixSymmetric(k)
          endDo
        endDo
      case default
        call mqc_error('mqc_packedSymmetricMatrix2FullMatrix_integer: Invalid upperLower.'//  &
          ' upperLower='//TRIM(upperLower)//'.')
      end select
!
      return
      end Subroutine mqc_packedSymmetricMatrix2FullMatrix_integer


!PROCEDURE mqc_packedSymmetricMatrix2FullMatrix_real
      Subroutine mqc_packedSymmetricMatrix2FullMatrix_real(matrixSymmetric,  &
        matrixFull,upperLower)
!
!     This subroutine accepts an input one-dimensional array of a symmetric
!     matrix (input dummy argument <matrixSymmetric>) stored in packed form and
!     returns a full matrix. (output dummy argument <matrixFull>). Note that
!     <matrixFull> is initialized in this routine. Furthermore, <matrixFull>
!     must be allocated prior to calling this routine.
!
!     Optional input dummy argument <upperLower> is a character string with 'U'
!     or 'L' to denote whether <matrixSymmetric> is sent in an upper- or
!     lower-symmetric column-wise storage format. If not sent, the default value
!     for <upperLower> is 'U'.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),intent(in)::matrixSymmetric
      real(kind=real64),dimension(:,:),intent(out)::matrixFull
      integer(kind=int64)::i,j,k,nDim,nSymm
      character(len=*),intent(in),optional::upperLower
      character(len=1)::myUpperLower
!
!
!     Error traps...
!
      if(SIZE(matrixFull,1).ne.SIZE(matrixFull,2))  &
        call mqc_error_i('mqc_packedSymmetricMatrix2FullMatrix_real: matrixFull array is NOT square!',  &
        6,'nDim(1)',SIZE(matrixFull,1),'nDim(2)',SIZE(matrixFull,2))
      nDim = SIZE(matrixFull,1)
      nSymm = (nDim*(nDim+1))/2
      if(nSymm.ne.SIZE(matrixSymmetric,1))  &
        call mqc_error_i('mqc_packedSymmetricMatrix2FullMatrix_real: matrixSymmetric has wrong dimension!',  &
        6,'nSymm',nSymm,'SIZE(matrixSymmetric)',SIZE(matrixSymmetric,1))
!
!     Set-up myUpperLower.
!
      if(PRESENT(upperLower)) then
        myUpperLower = TRIM(upperLower)
      else
        myUpperLower = 'U'
      endIf
!
!     Convert the symmetric matrix to full storage.
!
#ifdef PGI
      matrixFull = dfloat(0)
#else
      matrixFull = float(0)
#endif
      select case(myUpperLower)
      case('u','U')
        k = 0
        do i = 1,nDim
          do j = 1,i
            k = k+1
            matrixFull(i,j) = matrixSymmetric(k)
            matrixFull(j,i) = matrixSymmetric(k)
          endDo
        endDo
      case('l','L')
        k = 0
        do i = 1,nDim
          do j = i,nDim
            k = k+1
            matrixFull(i,j) = matrixSymmetric(k)
            matrixFull(j,i) = matrixSymmetric(k)
          endDo
        endDo
      case default
        call mqc_error('mqc_packedSymmetricMatrix2FullMatrix_real: Invalid upperLower.'//  &
          ' upperLower='//TRIM(upperLower)//'.')
      end select
!
      return
      end Subroutine mqc_packedSymmetricMatrix2FullMatrix_real


!PROCEDURE mqc_matrixSymm2Full_integer
      function mqc_matrixSymm2Full_integer(matrixSymmetric,upperLower) result(matrixFull)
!
!     This function accepts an input one-dimensional array of a symmetric matrix
!     (input dummy argument <matrixSymmetric>) stored in packed form and returns
!     a full matrix.
!
!     Optional input dummy argument <upperLower> is a character string with 'U'
!     or 'L' to denote whether <matrixSymmetric> is sent in an upper- or
!     lower-symmetric column-wise storage format. If not sent, the default value
!     for <upperLower> is 'U'.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),dimension(:),intent(in)::matrixSymmetric
      integer(kind=int64),dimension(:,:),allocatable::matrixFull
      integer(kind=int64)::iTmp,nDim,nSymm
      character(len=*),intent(in),optional::upperLower
      character(len=1)::myUpperLower
!
!     Set-up myUpperLower and alloacte matrixFull.
!
      if(PRESENT(upperLower)) then
        myUpperLower = TRIM(upperLower)
      else
        myUpperLower = 'U'
      endIf
      nSymm = SIZE(matrixSymmetric)
      iTmp = mqc_isqrt(1+8*nSymm) - 1
      nDim = iTmp/2
      allocate(matrixFull(nDim,nDim))
      call mqc_packedSymmetricMatrix2FullMatrix_integer(matrixSymmetric,  &
        matrixFull,myUpperLower)
!
      return
      end Function mqc_matrixSymm2Full_integer


!PROCEDURE mqc_matrixSymm2Full_real
      function mqc_matrixSymm2Full_real(matrixSymmetric,upperLower) result(matrixFull)
!
!     This function accepts an input one-dimensional array of a symmetric matrix
!     (input dummy argument <matrixSymmetric>) stored in packed form and returns
!     a full matrix.
!
!     Optional input dummy argument <upperLower> is a character string with 'U'
!     or 'L' to denote whether <matrixSymmetric> is sent in an upper- or
!     lower-symmetric column-wise storage format. If not sent, the default value
!     for <upperLower> is 'U'.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),intent(in)::matrixSymmetric
      real(kind=real64),dimension(:,:),allocatable::matrixFull
      integer(kind=int64)::iTmp,nDim,nSymm
      character(len=*),intent(in),optional::upperLower
      character(len=1)::myUpperLower
!
!     Set-up myUpperLower and alloacte matrixFull.
!
      if(PRESENT(upperLower)) then
        myUpperLower = TRIM(upperLower)
      else
        myUpperLower = 'U'
      endIf
      nSymm = SIZE(matrixSymmetric)
      iTmp = mqc_isqrt(1+8*nSymm) - 1
      nDim = iTmp/2
      allocate(matrixFull(nDim,nDim))
      call mqc_packedSymmetricMatrix2FullMatrix_real(matrixSymmetric,  &
        matrixFull,myUpperLower)
!
      return
      end Function mqc_matrixSymm2Full_real


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
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),dimension(:),intent(in)::matrixDiagonal
      integer(kind=int64),dimension(:,:),intent(out)::matrixFull
      integer(kind=int64)::nDim,i
!
      nDim = SIZE(matrixDiagonal)
      if(nDim.ne.SIZE(matrixFull,1).or.nDim.ne.SIZE(matrixFull,2))  &
        call mqc_error_i('mqc_packedDiagonalMatrix2FullMatrix_integer: disagreement in input/output array size!', 6, &
        'nDim', nDim, &
        'SIZE(matrixFull,1)', SIZE(matrixFull,1), &
        'SIZE(matrixFull,2)', SIZE(matrixFull,2))
      matrixFull = 0
      do i = 1,nDim
        matrixFull(i,i) = matrixDiagonal(i)
      endDo
!
      return
      end Subroutine mqc_packedDiagonalMatrix2FullMatrix_integer


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
      real(kind=real64),dimension(:),intent(in)::matrixDiagonal
      real(kind=real64),dimension(:,:),intent(out)::matrixFull
      integer(kind=int64)::nDim,i
!
      nDim = SIZE(matrixDiagonal)
      if(nDim.ne.SIZE(matrixFull,1).or.nDim.ne.SIZE(matrixFull,2))  &
        call mqc_error_i('mqc_packedDiagonalMatrix2FullMatrix_real: disagreement in input/output array size!', 6, &
        'nDim', nDim, &
        'SIZE(matrixFull,1)', SIZE(matrixFull,1), &
        'SIZE(matrixFull,2)', SIZE(matrixFull,2))
#ifdef PGI
      matrixFull = dfloat(0)
#else
      matrixFull = float(0)
#endif
      do i = 1,nDim
        matrixFull(i,i) = matrixDiagonal(i)
      endDo
!
      return
      end Subroutine mqc_packedDiagonalMatrix2FullMatrix_real
!
!
!PROCEDURE mqc_packedSymmetricSymmetricR4Tensor2Full_real
      subroutine mqc_packedSymmetricSymmetricR4Tensor2Full_real(r4TensorSymSym,  &
        r4TensorFull)
!
!     This subroutine accepts a rank-4 array packed in symmetric-symmetric form
!     and produces the corresponding fully un-packed rank-4 array. The tensor
!     must have the same size for all four dimensions.
!
!     A symmetric-symmetric rank-4 tensor has the following 8-fold equivalences:
!           (i,j,k,l) = (j,i,k,l)
!                     = (i,j,l,k)
!                     = (j,i,l,k)
!                     = (k,l,i,j)
!                     = (l,k,i,j)
!                     = (k,l,j,i)
!                     = (l,k,j,i)
!
!
!     H. P. Hratchian, 2020.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),intent(in)::r4TensorSymSym
      real(kind=real64),dimension(:,:,:,:),intent(out)::r4TensorFull
      integer(kind=int64)::i,j,k,l,ij,kl,ijkl,nDim,lenSymSym
      real(kind=real64)::temp
!
!     Error traps...
!
      if(SIZE(r4TensorFull,1).ne.SIZE(r4TensorFull,2).or.  &
        SIZE(r4TensorFull,1).ne.SIZE(r4TensorFull,3).or.   &
        SIZE(r4TensorFull,1).ne.SIZE(r4TensorFull,4))      &
        call mqc_error('mqc_packedSymmetricSymmetricR4Tensor2Full_real: Incorrect dimmensions!')
!
!     Do the work...
!
      nDim = SIZE(r4TensorFull,1)
      i = (nDim*(nDim+1))/2
      lenSymSym = (i*(i+1))/2
      do i = 1,nDim
        do j = 1,i
          IJ = (i*(i-1))/2+j
          do k = 1,nDim
            do l = 1,k
              KL = (k*(k-1))/2+l
              if(IJ.ge.KL) then
                IJKL = (IJ*(IJ-1))/2+KL
                temp = r4TensorSymSym(IJKL)
                r4TensorFull(i,j,k,l) = temp
                r4TensorFull(j,i,k,l) = temp
                r4TensorFull(i,j,l,k) = temp
                r4TensorFull(j,i,l,k) = temp
                r4TensorFull(k,l,i,j) = temp
                r4TensorFull(l,k,i,j) = temp
                r4TensorFull(k,l,j,i) = temp
                r4TensorFull(l,k,j,i) = temp
              endIf
            endDo
          endDo
        endDo
      endDo
!
      return
      end subroutine mqc_packedSymmetricSymmetricR4Tensor2Full_real
!
!
!PROCEDURE mqc_diagMatrixTimesFullMatrix
      subroutine mqc_diagMatrixTimesFullMatrix(diagMatrix,fullMatrix,outputMatrix)
!
!     This subroutine computes the product of a diagonal packed matrix, <diagMatrix>,
!     and a full matrix. The output is a full matrix with the same dimensions of
!     <fullMatrix>, which is return in dummy argument <outputMatrix>.
!     <outputMatrix> is an allocatable array.
!
!
!     H. P. Hratchian, 2020
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::diagMatrix
      real(kind=real64),dimension(:,:),intent(in)::fullMatrix
      real(kind=real64),dimension(:,:),allocatable,intent(out)::outputMatrix
!
      integer::i,M,N
!
!
!     Do some initial checking of the input arguments to ensure comformability
!     of <diagMatrix> and <fullMatrix>.
!
      M = SIZE(diagMatrix)
      N = SIZE(fullMatrix,2)
      if(SIZE(fullMatrix,1).ne.M) call mqc_error(  &
        'mqc_diagMatrixTimesFullMatrix: Diagonal and full matrices are NOT conformable.')
!
!     Allocate <outputMatrix>.
!
      if(Allocated(outputMatrix)) then
        if(SIZE(outputMatrix,1).ne.M.or.SIZE(outputMatrix,2).ne.N) then
          DeAllocate(outputMatrix)
          Allocate(outputMatrix(M,N))
        endIf
      else
        Allocate(outputMatrix(M,N))
      endIf
!
!     Do the matrix multiplication now...
!
      do i = 1,M
        outputMatrix(i,:) = diagMatrix(i)*fullMatrix(i,:)
      endDo
!
      return
      end subroutine mqc_diagMatrixTimesFullMatrix
!
!
!PROCEDURE mqc_matrixInverse_symmFull
      subroutine mqc_matrixInverse_symmFull(matrix,matrixInverse,  &
        eigenValues,eigenVectors)
!
!     This subroutine computes the inverse of input argument <matrix>, which is
!     a symmetric matrix sent in full/dense memory storage format. The output is
!     argument <matrixInverse>, which should be allocatable and will use
!     full/dense memory storage format.
!
!     The dummy arguments <eigenValues> and <eigenVectors> are OPTIONAL. If they
!     are sent, it is assumed they have already been formed. Otherwise, this
!     routine begins by forming the eigen-values and eigen-vectors of <matrix>
!     in scratch memory space.
!     
!
!     H. P. Hratchian, 2020
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::matrix
      real(kind=real64),dimension(:,:),allocatable,intent(out)::matrixInverse
      real(kind=real64),dimension(:),allocatable,OPTIONAL,intent(in)::eigenValues
      real(kind=real64),dimension(:,:),allocatable,OPTIONAL,intent(in)::eigenVectors
!
      integer::i,N
!
!
!     Do some initial checking of matrix to ensure it's square. Next, figure out
!     the dimension of <matrix>, <N>. Then, allocate space for <matrixInverse>
!     if needed.
!
      if(SIZE(matrix,1).ne.SIZE(matrix,2)) call mqc_error(  &
        'mqc_matrixInverse_symmFull: Rectangular matrix sent.')
      N = SIZE(matrix,1)
      if(Allocated(matrixInverse)) then
        if(SIZE(matrixInverse,1).ne.N.or.SIZE(matrixInverse,2).ne.N) then
          DeAllocate(matrixInverse)
          Allocate(matrixInverse(N,N))
        endIf
      else
        Allocate(matrixInverse(N,N))
      endIf
!!
!!    If <eigenValues> and <eigenVectors> have both been sent, form the
!!    generalized inverse of <matrix>. If not, diagonalize <matrix> and then form
!!    <matrixInverse>.
!!
!      If(PRESENT(eigenValues).and.PRESENT(eigenVectors)) then
!        do i = 1,N
!
!!
!!    Call LAPCK routine DSYEV. 
!!
!      call dsyev(JOBZ,'U',N,matrixCopy,N,eigenValues,WorkTmp,-1,INFO)    
!      if(INFO.ne.0) call mqc_error(  &
!        'mqc_eigensystem_symmMatrixFull: Failure in DSYEV getting WORK length.')
!      LWork = WorkTmp(1)
!      Allocate(Work(LWork))
!      call dsyev(JOBZ,'U',N,matrixCopy,N,eigenValues,Work,LWork,INFO)    
!      if(INFO.ne.0) call mqc_error(  &
!        'mqc_eigensystem_symmMatrixFull: Failure in DSYEV diagonalizing matrix.')
!      if(PRESENT(eigenVectors)) then
!        eigenVectors = matrixCopy
!        do i = 1,N
!          call mqc_vectorPhase(eigenVectors(:,i),.true.)
!        endDo
!      endIf
!
!      deAllocate(matrixCopy,Work)
!
      return
      end subroutine mqc_matrixInverse_symmFull
!




!
!
!----------------------------------------------------------------
!                                                               |
!     Other                                                     |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE mqc_isqrt
      function mqc_isqrt(iArg) result(iResult)
!
!     This function carries out integer sqrt.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),intent(in)::iArg
      integer(kind=int64)::iResult
      real(kind=real64)::tmp
!
#ifdef PGI
      tmp = dfloat(iArg)
#else
      tmp = float(iArg)
#endif
      iResult = int(sqrt(tmp))
!
      return
      end function mqc_isqrt
!
!
INCLUDE 'mqc_general_lapack.F03'
!
!
      End Module MQC_General
