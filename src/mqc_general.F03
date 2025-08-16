      Module MQC_General
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
!     **                          Version 25.8.0                           **
!     **                          August 1, 2025                           **
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
!      |          1         2         3         4         5         6     |
!      |123456789012345678901234567890123456789012345678901234567890123456|
!      |                                *                                 |
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
!           (4)  Algebra;
!           (5)  Other; and
!           (6)  BLAS and LAPACK wrappers.
!
!
!
!
!     Define Types and Classes.
!
      use iso_fortran_env, only: int64, real64
      implicit none
      real(kind=real64),public,parameter::angPBohr=0.52917706d0                                                          ! Angstom per Bohr
      real(kind=real64),public,parameter::kGPAMU=1.660538921d-27                                                         ! Kilograms per atomic mass unit
      real(kind=real64),public,parameter::planck=6.62606957d-34                                                          ! Planck constant (J.s)
      real(kind=real64),public,parameter::avogadro=6.02214129d+23                                                        ! Avogadro's number
      real(kind=real64),public,parameter::jPCal=4.184d00                                                                 ! Joules per calorie
      real(kind=real64),public,parameter::jPHartree=4.35974434d-18                                                       ! Joules per Hartree
      real(kind=real64),public,parameter::evPHartree=27.211399d0                                                         ! eV per Hartree
      real(kind=real64),public,parameter::cmM1PHartree=219474.63d0                                                       ! cm^-1 per Hartree
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
      real(kind=real64),public,parameter::MQC_small = 1.0d-6                                                             ! 10^-6, used as a threshold for near zero.
      complex(kind=real64),public,parameter::im = cmplx(0.0,1.0)                                                         ! Imaginary number.
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
        module procedure MQC_Print_Scalar_Complex
        module procedure MQC_Print_Vector_Array_Integer
        module procedure MQC_Print_Vector_Array_Real
        module procedure MQC_Print_Vector_Array_Complex
        module procedure MQC_Print_Matrix_Array_Integer
        module procedure MQC_Print_Matrix_Array_Real
        module procedure MQC_Print_Matrix_Array_Complex
        module procedure MQC_Print_Matrix_Array_Logical
        module procedure MQC_Print_R4Tensor_Array_Integer
        module procedure MQC_Print_R4Tensor_Array_Real
        module procedure MQC_Print_R4Tensor_Array_Complex
      end interface
!
      interface mqc_print_scalar
        module procedure MQC_Print_Scalar_Integer
        module procedure MQC_Print_Scalar_Real
        module procedure MQC_Print_Scalar_Complex
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
        module procedure MQC_Print_Matrix_Array_Logical
      end interface
!
      interface mqc_print_r4tensor
        module procedure MQC_Print_R4Tensor_Array_Integer
        module procedure MQC_Print_R4Tensor_Array_Real
        module procedure MQC_Print_R4Tensor_Array_Complex
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
      interface num2char
        module procedure integer2character
        module procedure real2character
        module procedure complex2character
      end interface
!
      interface sort
        module procedure mqc_bubbleSort_integer
        module procedure mqc_bubbleSort_real
      end interface
!
      interface matrixOrderedColumns
        module procedure mqc_matrixOrderedColumns_integer
        module procedure mqc_matrixOrderedColumns_real
      end interface
!
      interface flatten
        module procedure mqc_flattenR4Real
      end interface
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
      if(PRESENT(major))         major         = 25
      if(PRESENT(minor))         minor         = 8
      if(PRESENT(revision))      revision      = 0
      if(PRESENT(versionString)) versionString = '25.8.0'
!
      return
      end subroutine mqc_version
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
      integer::versionRequested,versionCurrent
      logical::ok
      integer::major,minor,revision
!
      ok = .true.
      call mqc_version(major,minor,revision)
!
!     If requested, do the test of match of requested and actual version.
!
      if(PRESENT(isMajor).or.PRESENT(isMinor).or.PRESENT(isRevision)) then
        versionRequested = 0
        versionCurrent   = 0
        if(PRESENT(isMajor)) then
          versionRequested = 10000*isMajor
          versionCurrent   = 10000*major
        endIf
        if(PRESENT(isMinor)) then
          versionRequested = versionRequested + 100*isMinor
          versionCurrent   = versionCurrent   + 100*minor
        endIf
        if(PRESENT(isRevision)) then
          versionRequested = versionRequested + isRevision
          versionCurrent   = versionCurrent   + revision
        endIf
        ok = (versionCurrent.eq.versionRequested)
      endIf
!
!     If requested, do an "older than" test.
!
      if(PRESENT(olderThanMajor).or.PRESENT(olderThanMinor).or.  &
        PRESENT(olderThanRevision)) then
        versionRequested = 0
        versionCurrent   = 0
        if(PRESENT(olderThanMajor)) then
          versionRequested = 10000*olderThanMajor
          versionCurrent   = 10000*major
        endIf
        if(PRESENT(olderThanMinor)) then
          versionRequested = versionRequested + 100*olderThanMinor
          versionCurrent   = versionCurrent   + 100*minor
        endIf
        if(PRESENT(olderThanRevision)) then
          versionRequested = versionRequested + olderThanRevision
          versionCurrent   = versionCurrent   + revision
        endIf
        ok = (versionCurrent.le.versionRequested)
      endIf
!
!     If requested, do a "newer than" test.
!
      if(PRESENT(newerThanMajor).or.PRESENT(newerThanMinor).or.  &
        PRESENT(newerThanRevision)) then
        versionRequested = 0
        versionCurrent   = 0
        if(PRESENT(newerThanMajor)) then
          versionRequested = 10000*newerThanMajor
          versionCurrent   = 10000*major
        endIf
        if(PRESENT(newerThanMinor)) then
          versionRequested = versionRequested + 100*newerThanMinor
          versionCurrent   = versionCurrent   + 100*minor
        endIf
        if(PRESENT(newerThanRevision)) then
          versionRequested = versionRequested + newerThanRevision
          versionCurrent   = versionCurrent   + revision
        endIf
        ok = ok.and.(versionCurrent.ge.versionRequested)
      endIf
!
      return
      end function mqc_version_check
!
!PROCEDURE mqc_version_print
      subroutine mqc_version_print(iUnit)
!
!     This subroutine prints the MQCPack version. Input dummy argument <iUnit>
!     is an integer giving the unit number of the file where the version
!     information should be printed.
!
!
      implicit none
      integer,intent(in)::iUnit
      integer::major,minor,revision
!
 1000 format(/,1x,'This program uses the Merced Quantum Chemistry Package (MQCPack).',/,  &
               1x,'MQCPack Version ',A,'.',A,'.',A,/)
!
!
!     Get the version information from routine mqc_version and then print it out
!     to file iUnit.
!
      call mqc_version(major,minor,revision)
      write(iUnit,1000) TRIM(num2char(major)),TRIM(num2char(minor)),TRIM(num2char(revision))
!
      return
      end subroutine mqc_version_print
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
      if(allocated(argument)) deallocate(argument)
      call get_command_argument(argNum,length=argLen)
      allocate(character(len=argLen)::argument)
      call get_command_argument(argNum,argument)
!
      Return
      End Subroutine mqc_get_command_argument

!
!PROCEDURE mqc_get_command_argument_integer
      Subroutine mqc_get_command_argument_integer(argNum,iArgument)
!
!     This subroutine is used to get an integer argument from the command line.
!     The output dummy argument iArgument is returned with the integer found as
!     command line arugment argNum. If this command line argument is NOT an
!     integer, the program will SegFault.
!     
!     -H. P. Hratchian, 2021.
!
!
      implicit none
      integer(kind=int64),intent(in)::argNum
      integer(kind=int64),intent(out)::iArgument
      character(len=:),allocatable::argument
!
!     Do the work...
!
      call mqc_get_command_argument(argNum,argument)
      read(argument,*) iArgument
!
      Return
      End Subroutine mqc_get_command_argument_integer

!
!PROCEDURE mqc_get_command_argument_real
      Subroutine mqc_get_command_argument_real(argNum,rArgument)
!
!     This subroutine is used to get a real argument from the command line. The
!     output dummy argument rArgument is returned with the real number found as
!     command line arugment argNum. If this command line argument is NOT an
!     real, the program will SegFault.
!     
!     -H. P. Hratchian, 2024.
!
!
      implicit none
      integer(kind=int64),intent(in)::argNum
      real(kind=real64),intent(out)::rArgument
      character(len=:),allocatable::argument
!
!     Do the work...
!
      call mqc_get_command_argument(argNum,argument)
      read(argument,*) rArgument
!
      Return
      End Subroutine mqc_get_command_argument_real


!
!----------------------------------------------------------------
!                                                               |
!     Printing                                                  |
!                                                               |
!----------------------------------------------------------------
!
!
!PROCEDURE MQC_Print_Scalar_Integer
      Subroutine MQC_Print_Scalar_Integer(scalar,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom,formatStr)
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
      Character(Len=*),Intent(In),optional::Header,formatStr
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
      if(PRESENT(formatStr)) then
        lineToPrint = TRIM(lineToPrint)//' '//TRIM(integer2character(scalar,formatStr))
      else
        lineToPrint = TRIM(lineToPrint)//' '//TRIM(integer2character(scalar))
      endIf
      lineToPrint = ADJUSTL(lineToPrint)
      write(iOut,'(1x,A)') TRIM(lineToPrint)
      if(PRESENT(Blank_At_Bottom)) then
        if(Blank_At_Bottom) write(iOut,1020)
      endIf
!
      return
      end subroutine MQC_Print_Scalar_Integer


!
!PROCEDURE MQC_Print_Scalar_Real
      Subroutine MQC_Print_Scalar_Real(scalar,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom,formatStr)
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
      Character(Len=*),Intent(In),optional::Header,formatStr
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
      if(PRESENT(formatStr)) then
        lineToPrint = TRIM(lineToPrint)//' '//TRIM(real2character(scalar,formatStr))
      else
        lineToPrint = TRIM(lineToPrint)//' '//TRIM(real2character(scalar))
      endIf
      lineToPrint = ADJUSTL(lineToPrint)
      write(iOut,'(1x,A)') TRIM(lineToPrint)
      if(PRESENT(Blank_At_Bottom)) then
        if(Blank_At_Bottom) write(iOut,1020)
      endIf
!
      return
      end subroutine MQC_Print_Scalar_Real
!
!
!PROCEDURE MQC_Print_Scalar_Complex
      Subroutine MQC_Print_Scalar_Complex(scalar,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom,formatStr)
!
!     This subroutine is used to print complex scalar.
!
!     L. M. Thompson, 2023.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::iOut
      Complex(kind=real64),Intent(In)::scalar
      Character(Len=*),Intent(In),optional::Header,formatStr
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      character(len=256)::lineToPrint,FormatLine
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
      if(PRESENT(formatStr)) then
        FormatLine = '(1X,A,1X,"=",1X,'//trim(FormatStr)//',SP,'//trim(FormatStr)//',"i")'
      else
        FormatLine = '(1x,A,1x,"=",1x,F14.6,SP,F14.6,"i")'
      endIf
      lineToPrint = ADJUSTL(lineToPrint)
      write(iOut,trim(FormatLine)) TRIM(lineToPrint), scalar
      if(PRESENT(Blank_At_Bottom)) then
        if(Blank_At_Bottom) write(iOut,1020)
      endIf
!
      return
      end subroutine MQC_Print_Scalar_Complex
!
!
!PROCEDURE MQC_Print_Rank4Tensor_Array_Real
      Subroutine MQC_Print_Rank4Tensor_Array_Real(Array,iOut,Header,  &
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
 1000 Format(1x,A)
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
      Subroutine MQC_Print_Matrix_Array_Integer(Array,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom,FormatStr)
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
      Character(Len=*),Intent(In),Optional::Header,FormatStr
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64),Parameter::ColWidth=10
      Character(Len=256)::FormatLine
!
!
 1000 Format(1x,A)
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
          If(present(FormatStr)) then
            FormatLine = '(1x,I7,2x,10'//trim(FormatStr)//')'
            Write(IOut,trim(FormatLine)) I, (Array(I,J),J=IFirst,ILast)
          Else
            Write(IOut,2001) I, (Array(I,J),J=IFirst,ILast)
          EndIf
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
      Subroutine MQC_Print_Matrix_Array_Real(Array,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom,FormatStr,PrintGT)
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
      Character(Len=*),Intent(In),Optional::Header,FormatStr
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      real(kind=real64),intent(in),optional::PrintGT
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64),Parameter::ColWidth=10
      Character(Len=256)::FormatLine
      logical::printLine
!
 1000 Format(1x,A)
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
          if(present(PrintGT)) then
            printLine = MAXVAL(Array(I,:)).gt.PrintGT
          else
            printLine = .true.
          endIf
          If(present(FormatStr)) then
            FormatLine = '(1x,I7,2x,10'//trim(FormatStr)//')'
            if(printLine) Write(IOut,trim(FormatLine)) I, (Array(I,J),J=IFirst,ILast)
          else
            if(printLine) Write(IOut,2002) I, (Array(I,J),J=IFirst,ILast)
          endIf
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
        Blank_At_Bottom,FormatStr)
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
      Character(Len=*),Intent(In),Optional::Header,FormatStr
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64),Parameter::ColWidth=10 
      Character(Len=256)::FormatLine
!
 1000 Format(1x,A)
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
          If(present(FormatStr)) then
            FormatLine = '(1x,I7,2x,10('//trim(FormatStr)//','//trim(FormatStr)//',"i"))'
            Write(IOut,trim(FormatLine)) I, (Array(I,J),J=IFirst,ILast)
          else
            Write(IOut,2003) I, (Array(I,J),J=IFirst,ILast)
          endIf
        enddo
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Matrix_Array_Complex

!
!PROCEDURE MQC_Print_Matrix_Array_Logical
      Subroutine MQC_Print_Matrix_Array_Logical(Array,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a logical Matrix array variable.
!
!     Adapted from MQC_Print_Matrix_Array_Integer.
!     H. P. Hratchian, 2025.
!

!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Logical,dimension(:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64),Parameter::ColWidth=10
      Character(Len=256)::FormatLine
!
!
 1000 Format(1x,A)
 1001 Format(5x,10(7x,I7))
 1020 Format(" ")
 2001 Format(1x,I7,10L14)
!
      NRows = Size(Array,1)
      NCols = Size(Array,2)
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      Write(IOut,1000) TRIM(Header)
!
      IFirst = 1
 10   Continue
      ILast = MIN(IFirst+9,NCols)
      Write(IOut,1001) (J,J=IFirst,ILast)
!
      Do I = 1,NRows
         Write(IOut,2001) I, (Array(I,J),J=IFirst,ILast)
      End Do
!
      IFirst = ILast+1
      If(IFirst <= NCols) Go To 10
!
      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Matrix_Array_Logical

!
!PROCEDURE MQC_Print_Vector_Array_Integer
      subroutine MQC_Print_Vector_Array_Integer(Vector,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom,FormatStr)
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
      character(Len=*),Intent(In),Optional::Header,FormatStr
      logical(kind=int64),Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      integer(kind=int64)::I,Length
      Character(Len=256)::FormatLine
!
 1000 Format(1x,A)
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
        If(present(FormatStr)) then
          FormatLine = '(1x,I7,2x,'//trim(FormatStr)//')'
          Write(IOut,trim(FormatLine)) I, Vector(I)
        Else
          Write(IOut,1001) I, Vector(I)
        EndIf
      EndDo
      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Vector_Array_Integer


!PROCEDURE MQC_Print_Vector_Array_Real
      Subroutine MQC_Print_Vector_Array_Real(Vector,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom,FormatStr,PrintGT)
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
      Character(Len=*),Intent(In),Optional::Header,FormatStr
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      real(kind=real64),intent(in),optional::PrintGT
      Integer(kind=int64)::I,Length
      Character(Len=256)::FormatLine
      logical::printLine
!
 1000 Format(1x,A)
 1002 Format(1x,I7,2X,F14.6)
 1020 Format(" ")
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      Write(IOut,1000) TRIM(Header)
      Length = Size(Vector)
      Do I = 1, Length
        if(present(PrintGT)) then
          printLine = Vector(I).gt.PrintGT
        else
          printLine = .true.
        endIf
        If(present(FormatStr)) then
          FormatLine = '(1x,I7,2x,'//trim(FormatStr)//')'
          if(printLine) Write(IOut,trim(FormatLine)) I, Vector(I)
        Else
          if(printLine) Write(IOut,1002) I, Vector(I)
        endIf
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Vector_Array_Real


!PROCEDURE MQC_Print_Vector_Array_Complex
      Subroutine MQC_Print_Vector_Array_Complex(vector,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom,FormatStr)
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
      Character(Len=*),Intent(In),Optional::Header,FormatStr
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,Length
      Character(Len=256)::FormatLine
!
 1000 Format(1x,A)
 1002 Format(1x,I7,F12.5,F11.5,"i")
 1020 Format(" ")
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      Write(IOut,1000) TRIM(Header)
      Length = Size(Vector)
      Do I = 1, Length
        If(present(FormatStr)) then
          FormatLine = '(1x,I7,2x,'//trim(FormatStr)//',SP,'//trim(FormatStr)//',"i")'
          Write(IOut,trim(FormatLine)) I, Vector(I)
        Else
          Write(IOut,1002) I, Vector(I)
        EndIf
      EndDo

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Vector_Array_Complex

!
!PROCEDURE MQC_Print_Vector_Array_Logical
      subroutine MQC_Print_Vector_Array_Logical(Vector,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a Vector array type variable.
!
!     H. P. Hratchian, 2025.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),Intent(In)::IOut
      logical,Dimension(:),Intent(In)::Vector
      character(Len=*),Intent(In),Optional::Header
      logical(kind=int64),Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      integer(kind=int64)::I,Length
      Character(Len=256)::FormatLine
!
 1000 Format(1x,A)
 1001 Format(1x,I7,2x,L14)
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
      End Subroutine MQC_Print_Vector_Array_Logical

!
!PROCEDURE MQC_Print_R4Tensor_Array_Integer
      Subroutine MQC_Print_R4Tensor_Array_Integer(Array,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print an integer 4-dimensional array variable.
!
!     L. M. Thompson, 2020.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Integer(kind=int64),dimension(:,:,:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,K,L,D1,D2,D3,D4
!
!
 1000 Format(1x,A)
 1020 Format(" ")
 1100 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',I10)

      D1 = size(Array,1)
      D2 = size(Array,2)
      D3 = size(Array,3)
      D4 = size(Array,4)

      If(Present(blank_at_top)) then
        If(blank_at_top) Write(IOut,1020)
      EndIf
      If(Present(Header)) Write(IOut,1000) Trim(Header)
      Do I = 1, D1
        Do J = 1, D2 
          Do K = 1, D3
            Do l = 1, D4
              Write(IOut,1100) I,J,K,L,Array(I,J,K,L)
            EndDo
          EndDo
        EndDo
      EndDo
      If(present(blank_at_bottom)) then
        If(blank_at_bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_R4Tensor_Array_Integer
!
!
!PROCEDURE MQC_Print_R4Tensor_Array_Real
      Subroutine MQC_Print_R4Tensor_Array_Real(Array,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a real 4-dimensional array variable.
!
!     L. M. Thompson, 2020.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Real(kind=real64),dimension(:,:,:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,K,L,D1,D2,D3,D4
!
!
 1000 Format(1x,A)
 1020 Format(" ")
 1200 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',F15.8)

      D1 = size(Array,1)
      D2 = size(Array,2)
      D3 = size(Array,3)
      D4 = size(Array,4)

      If(Present(blank_at_top)) then
        If(blank_at_top) Write(IOut,1020)
      EndIf
      If(Present(Header)) Write(IOut,1000) Trim(Header)
      Do I = 1, D1
        Do J = 1, D2 
          Do K = 1, D3
            Do l = 1, D4
              Write(IOut,1200) I,J,K,L,Array(I,J,K,L)
            EndDo
          EndDo
        EndDo
      EndDo
      If(present(blank_at_bottom)) then
        If(blank_at_bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_R4Tensor_Array_Real
!
!
!PROCEDURE MQC_Print_R4Tensor_Array_Complex
      Subroutine MQC_Print_R4Tensor_Array_Complex(Array,iOut,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a complex 4-dimensional array variable.
!
!     L. M. Thompson, 2020.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Complex(kind=real64),dimension(:,:,:,:),Intent(In)::Array
      Character(Len=*),Intent(In),Optional::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,K,L,D1,D2,D3,D4
!
!
 1000 Format(1x,A)
 1020 Format(" ")
 1300 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',F12.5,F11.5,"i")

      D1 = size(Array,1)
      D2 = size(Array,2)
      D3 = size(Array,3)
      D4 = size(Array,4)

      If(Present(blank_at_top)) then
        If(blank_at_top) Write(IOut,1020)
      EndIf
      If(Present(Header)) Write(IOut,1000) Trim(Header)
      Do I = 1, D1
        Do J = 1, D2 
          Do K = 1, D3
            Do l = 1, D4
              Write(IOut,1300) I,J,K,L,Array(I,J,K,L)
            EndDo
          EndDo
        EndDo
      EndDo
      If(present(blank_at_bottom)) then
        If(blank_at_bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_R4Tensor_Array_Complex
!
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

!PROCEDURE complex2character
      function complex2character(complexIn,formatString) result(complexString)
!
!     This function converts a complex variable and returns a character string. The
!     input dummy argument <formatString> is OPTIONAL and can be sent to specify
!     the desired formatting of the complexIn. Note that <formatString> must
!     conform to standard fortran requirements.
!
!     L. M. Thompson, 2020.
!
      implicit none
      complex(kind=real64),intent(in)::complexIn
      character(*),intent(in),optional::formatString
      character(len=256)::myFormatString,complexString
!
      if(PRESENT(formatString)) then
        myFormatString = formatString
      else
        myFormatString = 'f20.5'
      endIf
      myFormatString = '('//TRIM(myFormatString)//',SP,'//TRIM(myFormatString)//'"i")'
      write(complexString,myFormatString) complexIn
      complexString = ADJUSTL(complexString)
!
      return
      end function complex2character


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
!PROCEDURE mqc_float
      function mqc_float(valIn) result(valOut)
!
!     This function is a wrapper for the float function.
!
!
!     H. P. Hratchian, 2025.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64)::valIn
      real(kind=real64)::valOut
!
#ifdef PGI
      valOut = dfloat(valIn)
#else
      valOut = float(valIn)
#endif
!
      return
      end function mqc_float
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
#ifdef PGI
        vector = dfloat(0)
#else
        vector = float(0)
#endif
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
!     have absolute values below a threshold to float(0). Arguments <vectorIn>
!     and <vectorOut> CANNOT be the same.
!
!     If the optional dummy argument <vectorOut> is sent, then <vectorIn> is
!     left unchanged and the zero-d version of the vector is returned in
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
      real(kind=real64),optional,intent(in)::small
      integer(kind=int64)::i
      real(kind=real64)::mySmall,zero=0.0
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
        vector(1) = mqc_float(0)
      largestMagnitudeValue = abs(vector(1)) + MQC_small
      largestMagnitudePosition = 1
      do i = 2,SIZE(vector)
        if(myTrimZero.and.(abs(vector(i)).lt.MQC_small))  &
          vector(i) = mqc_float(0)
        if(abs(vector(i)).gt.largestMagnitudeValue) then
          largestMagnitudeValue = abs(vector(i)) + MQC_small
          largestMagnitudePosition = i
        endIf
      endDo
      if(myLargestPositive) then
        if(vector(largestMagnitudePosition).lt.mqc_float(0)) vector = -vector
      else
        if(vector(largestMagnitudePosition).gt.mqc_float(0)) vector = -vector
      endIf
      if(myTrimZero) then
        do i = 1,SIZE(vector)
          if(abs(vector(i)).lt.MQC_small) vector(i) = mqc_float(0)
        endDo
      endIf
!
      return
      end subroutine mqc_vectorPhase
!
!PROCEDURE mqc_crossProduct3D_real
      function mqc_crossProduct3D_real(vector1,vector2) result(vectorOut)
!
!     This function forms the 3D vector cross product between two vectors --
!     vector1 and vector2. The result is vector vectorOut.
!
!
!     H. P. Hratchian, 2025.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),intent(in)::vector1,vector2
      real(kind=real64),allocatable,dimension(:)::vectorOut
!
      integer(kind=int64)::n1,n2,i,j
!
!     Do the work...
!
      n1 = SIZE(vector1)
      n2 = SIZE(vector2)
      if((n1.ne.3).or.(n2.ne.3))  &
        call mqc_error('MQC_CrossProduct3D: Vectors are not length 3.')
      Allocate(vectorOut(n1))
      vectorOut(1) = vector1(2)*vector2(3)-vector1(3)*vector2(2)
      vectorOut(2) = vector1(3)*vector2(1)-vector1(1)*vector2(3)
      vectorOut(3) = vector1(1)*vector2(2)-vector1(2)*vector2(1)
!
      return
      end function mqc_crossProduct3D_real
!
!PROCEDURE mqc_outerProduct_real
      Function mqc_outerProduct_real(vector1,vector2,scalar) result(matrixOut)
!
!     This function forms the outer product of two vectors -- vector1 and
!     vector2. The result is a matrix that is returned in matrixOut, which is
!     assumed to be allocatable. Optionally, a scalar <scalar> can be sent to
!     multiply the outer product.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),intent(in)::vector1,vector2
      real(kind=real64),OPTIONAL,intent(in)::scalar
      real(kind=real64),allocatable,dimension(:,:)::matrixOut
!
      integer(kind=int64)::n1,n2,i,j
!
!     Do the work...
!
      n1 = SIZE(vector1)
      n2 = SIZE(vector2)
      Allocate(matrixOut(n1,n2))
      if(PRESENT(scalar)) then
        do i = 1,n1
          do j = 1,n2
            matrixOut(i,j) = scalar*vector1(i)*vector2(j)
          endDo
        endDo
      else
        do i = 1,n1
          do j = 1,n2
            matrixOut(i,j) = vector1(i)*vector2(j)
          endDo
        endDo
      endIf
!
      return
      end function mqc_outerProduct_real
!
!
!PROCEDURE mqc_matrixTrimZero
      Subroutine mqc_matrixTrimZero(matrixIn,matrixOut,small)
!
!     This subroutine accepts a vector, <matrixIn>, and sets all entries that
!     have absolute values below a threshold to float(0). Arguments <matrixIn>
!     and <matrixOut> CANNOT be the same.
!
!     If the optional dummy argument <matrixOut> is sent, then <matrixIn> is
!     left unchanged and the zero-d version of the vector is returned in
!     <matrixOut>. Otherwise, the zero-ing is done in place in <matrixIn>.
!
!     The optional dummy argument <small> may be sent with a user specified
!     threshold. Otherwise, the threshold is set to the MQC system value of
!     <MQC_small>.
!
!
!     H. P. Hratchian, 2025.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:,:),intent(inOut)::matrixIn
      real(kind=real64),dimension(:,:),allocatable,optional,intent(Out)::matrixOut
      real(kind=real64),optional,intent(in)::small
      integer(kind=int64)::i,j
      real(kind=real64)::mySmall,zero=0.0
!
!
!     Initial set-up code...
!
      mySmall = MQC_small
      if(PRESENT(small)) mySmall = small
      if(PRESENT(matrixOut)) then
        if(Allocated(matrixOut)) then
          if((SIZE(matrixOut,1).ne.SIZE(matrixIn,1)).or.  &
            (SIZE(matrixOut,2).ne.SIZE(matrixIn,2))) then
            DeAllocate(matrixOut)
            Allocate(matrixOut(size(matrixIn,1),size(matrixIn,2)))
          endIf
        else
          Allocate(matrixOut(size(matrixIn,1),size(matrixIn,2)))
        endIf
      endIf
!
!     Do trim work.
!
      if(PRESENT(matrixOut)) then
        matrixOut = matrixIn
        do j = 1,SIZE(matrixIn,2)
          do i = 1,SIZE(matrixIn,1)
            if(abs(matrixOut(i,j)).lt.mySmall) matrixOut(i,j) = zero
          endDo
        endDo
      else
        do j = 1,SIZE(matrixIn,2)
          do i = 1,SIZE(matrixIn,1)
            if(abs(matrixIn(i,j)).lt.mySmall) matrixIn(i,j) = zero
          endDo
        endDo
      endIf
!
      return
      end Subroutine mqc_matrixTrimZero

!
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

!
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
      matrixFull = mqc_float(0)
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

!
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

!
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

!
!PROCEDURE mqc_packedDiagonalMatrix2FullMatrix_real
      subroutine mqc_packedDiagonalMatrix2FullMatrix_real(matrixDiagonal,  &
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
      matrixFull = mqc_float(0)
      do i = 1,nDim
        matrixFull(i,i) = matrixDiagonal(i)
      endDo
!
      return
      end Subroutine mqc_packedDiagonalMatrix2FullMatrix_real

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
!PROCEDURE mqc_matrixTrace
      function mqc_matrixTrace(matrix) result(trace)
!
!     This function returns the trace of the matrix <matrix>. For now, this
!     function only works for full square matrices.
!
!
!     H. P. Hratchian, 2025
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::matrix
      real(kind=real64)::trace
      integer::i
!
!     Evaluate the trace of the matrix.
!
      if(Size(matrix,1).ne.Size(matrix,2))  &
        call mqc_error('mqc_matrixTrace: Non-square matrix sent.')
      trace = mqc_float(0)
      do i = 1,Size(matrix,1)
        trace = trace + matrix(i,i)
      endDo
!
      return
      end function mqc_matrixTrace

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
!PROCEDURE MQC_leastSquaresFit
      subroutine MQC_leastSquaresFit(x,y,slope,intercept,rSquared)
!
!     This subroutine evaluates a 1D least squares fit for the data set x and
!     y(x). The output of the routine is the slope, intercept, and R**2 value.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::x,y
      real(kind=real64),intent(out)::slope,intercept,rSquared
      integer(kind=int64)::n,i
      real(kind=real64)::xMean,yMean
      real(kind=real64)::SSxy,SSxx,SSyy,SSR
!
!     Set n, the number of data points, and then do error trapping.
!
      n = size(x)
      if(n.ne.size(y))  &
        call mqc_error('Error: x and y arrays must have the same size.')
      if(n.lt.2)  &
        call mqc_error('Error: At least two points required for least-squares fitting.')
!
!     Solve the linear least squares problem.
!
      xMean     = sum(x)/mqc_float(n)
      yMean     = sum(y)/mqc_float(n)
      SSxy      = sum((x-xMean)*(y-yMean))
      SSxx      = sum((x-xMean)**2)
      SSyy      = sum((y-yMean)**2)
      slope     = SSxy/SSxx
      intercept = yMean-slope*xMean
      SSR       = SSxx*slope**2
      rSquared  = SSR/SSyy
!
      return
      end subroutine MQC_leastSquaresFit

!
!
!----------------------------------------------------------------
!                                                               |
!     Other                                                     |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE mqc_bubbleSort_integer
      subroutine mqc_bubbleSort_integer(listIn,listOut,map,sortListIn,reverse)
!
!     This subroutine carries out a simple bubble sort algorithm to order the
!     values in the array <listIn>. The sorted list is returned in the optional
!     argument <listOut> if it is sent. The mapping of the unsorted list to the
!     sorted list is returned in the optional argument <map>. If neither
!     <listOut> or <map> is sent, then the sorted list will overwrite the input
!     data in listIn. If <listOut> and/or <map> are sent, it should already be
!     appropriately allocated. If <sortListIn> is sent and is TRUE, then listIn
!     is always returned with the sorted values. Optional argument <reverse> is
!     sent as TRUE if the list should be sorted from largest to smallest value.
!     The default value of <reverse> is FALSE.
!
!     H. P. Hratchian, 2023.
!
!
!     Variabile Declarations
      implicit none
      integer(kind=int64),dimension(:)::listIn
      integer(kind=int64),dimension(:),optional::listOut,map
      logical,optional::sortListIn,reverse
      integer(kind=int64)::i,j,nDim,nSwaps,valueTemp
      integer(kind=int64),dimension(:),allocatable::listTemp,listMapTemp
      logical::overwriteListIn,doReverse
!
!     Set up overwriteListIn and doReverse based on whether the related optional
!     dummy arguments have been sent by the calling program unit.
!
      overwriteListIn = .false.
      if(PRESENT(sortListIn)) overwriteListIn = sortListIn
      doReverse = .false.
      if(PRESENT(reverse)) doReverse = reverse
!
!     Allocate listTemp and copy listIn into it.
!
      nDim = Size(listIn)
      Allocate(listTemp(nDim),listMapTemp(nDim))
      listTemp = listIn
      if(doReverse) listTemp = -listTemp
      call mqc_seq(listMapTemp)
!
!     Carry out the bubble sort algorithm on listTemp.
!
      do i = 1,nDim
        nSwaps = 0
        do j = 1,nDim-i
          if(listTemp(j).gt.listTemp(j+1)) then
            nSwaps = nSwaps + 1
            valueTemp = listTemp(j)
            listTemp(j) = listTemp(j+1)
            listTemp(j+1) = valueTemp
            valueTemp = listMapTemp(j)
            listMapTemp(j) = listMapTemp(j+1)
            listMapTemp(j+1) = valueTemp
          endIf
        endDo
        if(nSwaps.eq.0) exit
      endDo
!
!     Put the sorted list into listOut or back into listIn.
!
      if(doReverse) listTemp = -listTemp
      if(PRESENT(listOut)) listOut = listTemp
      if(PRESENT(map)) map = listMapTemp
      if(overwriteListIn.or.(.not.(PRESENT(listOut).or.PRESENT(map))))  &
        listIn = listTemp
!
      return
      end subroutine mqc_bubbleSort_integer


!
!PROCEDURE mqc_bubbleSort_real
      subroutine mqc_bubbleSort_real(listIn,listOut,map,sortListIn,reverse)
!
!     This subroutine carries out a simple bubble sort algorithm to order the
!     values in the array <listIn>. The sorted list is returned in the optional
!     argument <listOut> if it is sent. The mapping of the unsorted list to the
!     sorted list is returned in the optional argument <map>. If neither
!     <listOut> or <map> is sent, then the sorted list will overwrite the input
!     data in listIn. If <listOut> and/or <map> are sent, it should already be
!     appropriately allocated. If <sortListIn> is sent and is TRUE, then listIn
!     is always returned with the sorted values. Optional argument <reverse> is
!     sent as TRUE if the list should be sorted from largest to smallest value.
!     The default value of <reverse> is FALSE.
!
!
!     H. P. Hratchian, 2023.
!
!
!     Variabile Declarations
      implicit none
      real(kind=real64),dimension(:)::listIn
      real(kind=real64),dimension(:),optional::listOut
      integer(kind=int64),dimension(:),optional::map
      logical,optional::sortListIn,reverse
      integer(kind=int64)::i,j,nDim,nSwaps
      integer(kind=int64),dimension(:),allocatable::listMapTemp
      real(kind=real64)::valueTemp
      real(kind=real64),dimension(:),allocatable::listTemp
      logical::overwriteListIn,doReverse
!
!     Set up overwriteListIn and doReverse based on whether the related optional
!     dummy arguments have been sent by the calling program unit.
!
      overwriteListIn = .false.
      if(PRESENT(sortListIn)) overwriteListIn = sortListIn
      doReverse = .false.
      if(PRESENT(reverse)) doReverse = reverse
!
!     Allocate listTemp and copy listIn into it.
!
      nDim = Size(listIn)
      Allocate(listTemp(nDim),listMapTemp(nDim))
      listTemp = listIn
      if(doReverse) listTemp = -listTemp
      call mqc_seq(listMapTemp)
!
!     Carry out the bubble sort algorithm on listTemp.
!
      do i = 1,nDim
        nSwaps = 0
        do j = 1,nDim-i
          if(listTemp(j).gt.listTemp(j+1)) then
            nSwaps = nSwaps + 1
            valueTemp = listTemp(j)
            listTemp(j) = listTemp(j+1)
            listTemp(j+1) = valueTemp
            valueTemp = listMapTemp(j)
            listMapTemp(j) = listMapTemp(j+1)
            listMapTemp(j+1) = valueTemp
          endIf
        endDo
        if(nSwaps.eq.0) exit
      endDo
!
!     Put the sorted list into listOut or back into listIn.
!
      if(doReverse) listTemp = -listTemp
      if(PRESENT(listOut)) listOut = listTemp
      if(PRESENT(map)) map = listMapTemp
      if(overwriteListIn.or.(.not.(PRESENT(listOut).or.PRESENT(map))))  &
        listIn = listTemp
!
      return
      end subroutine mqc_bubbleSort_real


!
!PROCEDURE mqc_matrixOrderedColumns_integer
      subroutine mqc_matrixOrderedColumns_integer(matrix,map)
!
!     This routine reordered the columns of the matrix <matrix> according to the
!     mapping array <map>.
!
!
!     H. P. Hratchian, 2023.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),dimension(:,:)::matrix
      integer(kind=int64),dimension(:)::map
      integer(kind=int64)::i
      integer(kind=int64),dimension(:,:),allocatable::tempMatrix
      
!
!     Begin by ensuring map is the same length as the number of columns in
!     matrix.
!
      if(SIZE(matrix,2).ne.SIZE(map))  &
        call mqc_error('mqc_matrixOrderedColumns: map is wrong length!')
!
!     Do the work.
!
      Allocate(tempMatrix(Size(matrix,1),Size(matrix,2)))
      do i = 1,Size(matrix,2)
        tempMatrix(:,i) = matrix(:,map(i))
      endDo
      matrix = tempMatrix
!
      return
      end subroutine mqc_matrixOrderedColumns_integer


!
!PROCEDURE mqc_matrixOrderedColumns_real
      subroutine mqc_matrixOrderedColumns_real(matrix,map)
!
!     This routine reordered the columns of the matrix <matrix> according to the
!     mapping array <map>.
!
!
!     H. P. Hratchian, 2023.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:,:)::matrix
      integer(kind=int64),dimension(:)::map
      integer(kind=int64)::i
      real(kind=real64),dimension(:,:),allocatable::tempMatrix
      
!
!     Begin by ensuring map is the same length as the number of columns in
!     matrix.
!
      if(SIZE(matrix,2).ne.SIZE(map))  &
        call mqc_error('mqc_matrixOrderedColumns: map is wrong length!')
!
!     Do the work.
!
      Allocate(tempMatrix(Size(matrix,1),Size(matrix,2)))
      do i = 1,Size(matrix,2)
        tempMatrix(:,i) = matrix(:,map(i))
      endDo
      matrix = tempMatrix
!
      return
      end subroutine mqc_matrixOrderedColumns_real


!
!PROCEDURE mqc_flatten
      function mqc_flattenR4Real(inArray) result(outArray)
!
!     This function flattens a rank-4 real array.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:,:,:,:)::inArray
      integer(kind=int64)::length
      real(kind=real64),dimension(:),allocatable::outArray
!
      length = SIZE(inArray)
      allocate(outArray(length))
      outArray = TRANSFER(inArray,inArray)
!
      return
      end function mqc_flattenR4Real


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
      tmp = mqc_float(iArg)
      iResult = int(sqrt(tmp))
!
      return
      end function mqc_isqrt


!
!PROCEDURE mqc_factorial
      function mqc_factorial(iArg) result(iResult)
!
!     This function returns iArg!.
!
!
!     H. P. Hratchian, 2024.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),intent(in)::iArg
      integer(kind=int64)::iResult
      integer(kind=int64)::i
!
      if(iArg.lt.0) call mqc_error('Factorial of negative number not allowed.')
      iResult = 1
      do i = iArg,1,-1
        iResult = iResult*i
      endDo
!
      return
      end function mqc_factorial


!
!PROCEDURE mqc_doubleFactorial
      function mqc_doubleFactorial(iArg) result(iResult)
!
!     This function returns iArg!!.
!
!
!     H. P. Hratchian, 2024.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),intent(in)::iArg
      integer(kind=int64)::iResult
      integer(kind=int64)::i
!
      if(iArg.eq.-1) then
        iResult = -1
      elseIf(iArg.eq.0) then
        iResult = 1
      elseIf(iArg.gt.0) then
        iResult = 1
        do i = iArg,1,-2
          iResult = iResult*i
        endDo
      else
        call mqc_error('Double factorial of negative number < -1 not allowed.')
      endIf
!
      return
      end function mqc_doubleFactorial
!
!PROCEDURE mqc_seq
      subroutine mqc_seq(list,start,step)
!
!     This subroutine fills <list> with a sequence of integers that start at
!     <start> and increment by <step>. The dummy arguments <start> and <step>
!     are optional. If <start> is not sent, it is set to 1. If <step> is not
!     sent, it is set to 1.
!
!
!     H. P. Hratchian, 2023.
!
!
!     Variable Declarations
      implicit none
      integer(kind=int64),dimension(:)::list
      integer(kind=int64),optional::start,step
      integer(kind=int64)::i,nDim,myStart,myStep
!
!     Do the work.
!
      nDim = SIZE(list)
      myStart = 1
      myStep  = 1
      if(PRESENT(start)) myStart = start
      if(PRESENT(step))  myStep  = step
      list(1) = myStart
      do i = 2,nDim
        list(i) = list(i-1) + myStep
      endDo
!
      return
      end subroutine mqc_seq
!
!PROCEDURE mqc_sphericalPolar2xyz
      subroutine mqc_sphericalPolar2xyz(r,theta,phi,x,y,z,inDegrees)
!
!     This routine converts spherical polar coordinates (sent as r, theta, phi)
!     to Cartesian coordiantes (returned in x, y, z). The logical argument
!     inDegrees is sent as .TRUE. or .FALSE. if theta and phi are sent in
!     degrees (TRUE) or radians (FALSE).
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),intent(in)::r,theta,phi
      real(kind=real64),intent(out)::x,y,z
      logical,intent(in)::inDegrees
      real(kind=real64)::deg2rad
!
      if(inDegrees) then
        deg2rad = pi/mqc_float(180)
        x = r*cos(deg2rad*theta)*sin(deg2rad*phi)
        y = r*sin(deg2rad*theta)*sin(deg2rad*phi)
        z = r*cos(deg2rad*phi)
      else
        x = r*cos(theta)*sin(phi)
        y = r*sin(theta)*sin(phi)
        z = r*cos(phi)
      endIf
!
      return
      end subroutine mqc_sphericalPolar2xyz
!
!
!----------------------------------------------------------------
!                                                               |
!     BLAS AND LAPACK WRAPPERS                                  |
!                                                               |
!----------------------------------------------------------------
!
!
!
INCLUDE 'mqc_general_lapack.F03'
!
!
      End Module MQC_General
