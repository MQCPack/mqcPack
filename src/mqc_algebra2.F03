      Module MQC_Algebra2
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
!     This module is a new (relative to mqc_algebra) module that provides
!     advanced numerical variable functionalities and easy-to-use algebra tools.
!
!     The functions and subroutines contained within this module are grouped
!     into 
!           1. MQC_Variable Management Procedures;
!           2. MQC_Variable Printing Procedures; and
!           3. MQC_Variable Value Setting & Basic Manipulation Procedures; and
!           4. MQC_Variable Mathematics Procedures;


!
      Use MQC_General
      use iso_fortran_env
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     Variable...
      type MQC_Variable
        integer(kind=int64),private::rank=-1
        character(len=64),private::dataType,storageFormat
        integer(kind=int64),dimension(10),private::dimensions=0
        real(kind=real64),dimension(:),allocatable,private::realArray
        integer(kind=int64),dimension(:),allocatable,private::integerArray
        logical,private::initialized=.false.
      Contains
        procedure,public::initialize    => MQC_Variable_initialize
        procedure,public::init          => MQC_Variable_initialize
        procedure,private::MQC_Variable_clear_mqc
        procedure,private::MQC_Variable_clear_integer
        procedure,private::MQC_Variable_clear_real
        generic,public   ::clear        => MQC_Variable_clear_mqc,  &
                                           MQC_Variable_clear_integer,  &
                                           MQC_Variable_clear_real
        procedure,private::MQC_Variable_put_MQC
        procedure,private::MQC_Variable_put_intrinsicInteger
        procedure,private::MQC_Variable_put_intrinsicReal
        generic,public   ::put          => MQC_Variable_put_MQC,  &
                                           MQC_Variable_put_intrinsicInteger,  &
                                           MQC_Variable_put_intrinsicReal
        procedure,public::getRank       => MQC_Variable_getRank
        procedure,public::getType       => MQC_Variable_getType
        procedure,public::print         => MQC_Print_mqcVariable
        procedure,public::isConformable => MQC_Variable_isConformable
      end type MQC_Variable
!
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!
!     Interface MQC Variables to the RANK intrinsic function.
      interface rank
        module procedure MQC_Variable_getRank
      end interface
!
!     Interface MQC Variables to the SIZE intrinsic function.
      interface size
        module procedure MQC_Variable_getSize
      end interface
!
!     Interface intrinsic real/integer to MQC variable.
      interface mqc
        module procedure MQC_Variable_from_intrinsicInteger
        module procedure MQC_Variable_from_intrinsicReal
      end interface
!
!     Interface MQC Varaiable --> Intrinsic conversion to INT and FLOAT.
      interface int
        module procedure MQC_Variable_INT_Scalar
      end interface
      interface float
        module procedure MQC_Variable_FLOAT_Scalar
      end interface
!
!     Interface MQC Varaiable full contraction for vectors --> dot_product.
      interface dot_product
        module procedure MQC_Variable_dot_product
      end interface
!
!     Interface routines that clear MQC variable objects.
      interface MQC_Variable_clear
        module procedure MQC_Variable_clear_mqc
        module procedure MQC_Variable_clear_integer
        module procedure MQC_Variable_clear_real
      end interface
!
!     Interface routines that fill MQC variable objects.
      interface mqc_variable_fillVal
        module procedure MQC_Variable_setVal_intrinsicScalar2mqcScalar_integer
        module procedure MQC_Variable_setVal_intrinsicScalar2mqcScalar_real
        module procedure MQC_Variable_setVal_intrinsicVector2mqcVector_integer
        module procedure MQC_Variable_setVal_intrinsicVector2mqcVector_real
        module procedure MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_integer
        module procedure MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_real
      end interface
!
!     Interface routines that fill specific elements of MQC variable objects.
      interface mqc_variable_put
        module procedure MQC_Variable_put_MQC
        module procedure MQC_Variable_put_intrinsicInteger
        module procedure MQC_Variable_put_intrinsicReal
      end interface
!
!     Interface print routines for MQC_Variable types.
      interface mqc_print
         module procedure MQC_Print_mqcVariable
      end interface
!
!
!----------------------------------------------------------------
!                                                               |
!     OPERATOR INTERFACES                                       |
!                                                               |
!----------------------------------------------------------------
!
!
!     Assignment and type conversion interface (=).
      interface assignment (=)
        module procedure MQC_Variable_setVal_intrinsicScalar2mqcScalar_integer
        module procedure MQC_Variable_setVal_intrinsicScalar2mqcScalar_real
        module procedure MQC_Variable_setVal_intrinsicVector2mqcVector_integer
        module procedure MQC_Variable_setVal_intrinsicVector2mqcVector_real
        module procedure MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_integer
        module procedure MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_real
        module procedure MQC_Variable_mqc2intrinsicIntegerScalar
        module procedure MQC_Variable_mqc2intrinsicRealScalar
      end interface
!
!     Whole-array addition (+).
      interface operator (+)
        module procedure MQC_Variable_Addition
      end interface
!
!     Whole-array addition (-).
      interface operator (-)
        module procedure MQC_Variable_Subtraction
      end interface
!
!     Whole-array multiplication (*).
      interface operator (*)
        module procedure MQC_Variable_Multiplication
      end interface
!
!     Whole-array division (/).
      interface operator (/)
        module procedure MQC_Variable_Division
      end interface
!
!     Conformable comparison operator (.conformable.)
      interface operator (.conformable.)
        module procedure MQC_Variable_isConformable
      end interface

!      Interface Operator (+)
!        Module Procedure MQC_ScalarAdd
!      End Interface
!      Interface Operator (-)
!        Module Procedure MQC_ScalarSubtract
!      End Interface
!      Interface Operator (*)
!        Module Procedure MQC_ScalarMultiply
!        Module Procedure MQC_ScalarVectorProduct
!        Module Procedure MQC_VectorScalarProduct
!        Module Procedure MQC_ScalarMatrixProduct
!        Module Procedure MQC_MatrixScalarProduct
!      End Interface
!      Interface Operator (/)
!        Module Procedure MQC_ScalarDivide
!      End Interface
!      Interface Operator (.ne.)
!        Module Procedure MQC_ScalarNE
!      End Interface
!      Interface Operator (.eq.)
!        Module Procedure MQC_ScalarEQ
!      End Interface
!      Interface Operator (.lt.)
!        Module Procedure MQC_ScalarLT
!      End Interface
!      Interface Operator (.gt.)
!        Module Procedure MQC_ScalarGT
!      End Interface
!      Interface Operator (.le.)
!        Module Procedure MQC_ScalarLE
!      End Interface
!      Interface Operator (.ge.)
!        Module Procedure MQC_ScalarGE
!      End Interface
!hph-
!
!
!
!----------------------------------------------------------------
!                                                               |
!     SUBROUTINES AND FUNCTIONS                                 |
!                                                               |
!----------------------------------------------------------------
!
!
      CONTAINS
!
!----------------------------------------------------------------
!                                                               |
!     MQC_Variable MANAGEMENT PROCEDURES                        |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE MQC_Variable_initialize
      subroutine MQC_Variable_initialize(mqcVariable,dataType,dimensions,  &
          rank,storageFormat)
!
!     This subroutine is used to initialize an MQC_Variable object.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(inOut)::mqcVariable
      character(len=*),intent(in)::dataType
      integer(kind=int64),dimension(:),intent(in),optional::dimensions
      integer(kind=int64),intent(in),optional::rank
      character(len=*),intent(in),optional::storageFormat
!
      character(len=64)::myDataType
      integer(kind=int64)::myRank,myLength
      character(len=64)::myStorageFormat
!
!
!     Set-up <myRank> and ensure the sent <rank> and <dimensions> are
!     compatable.
!
      call string_change_case(dataType,'u',myDataType)
      if(Present(rank)) then
        myRank = rank
      elseIf(Present(dimensions)) then
        myRank = Size(dimensions)
      else
        myRank = 0
      endIf
      if(Present(dimensions)) then
        if(Size(dimensions).ne.myRank)  &
          call MQC_Error_I('Logic error in MQC_Variable_initialize: Dimensions and Rank do not jive!', 6, &
          'Size(dimensions)', Size(dimensions), &
          'myRank', myRank )
      endIf
!
!     Set-up <myStorageFormat>.
!
      if(Present(storageFormat)) then
        myStorageFormat = storageFormat
      else
        myStorageFormat = 'FULL'
      endIf
      call string_change_case(myStorageFormat,'u')
      mqcVariable%storageFormat = TRIM(myStorageFormat)
!
!     Fill the dataType.
!
      select case(myDataType)
      case('REAL','R')
        mqcVariable%dataType = 'REAL'
      case('INTEGER','INT','I')
        mqcVariable%dataType = 'INTEGER'
      case default
        call MQC_Error_A('Initializing MQC_Variable with unknown data type.', 6, &
             'myDataType', myDataType)
      end select
!
!     Take care of rank and dimensions.
!
      if(mqcVariable%rank.ne.myRank) then
        mqcVariable%rank = myRank
        select case(myRank)
        case(0)
          mqcVariable%dimensions = 0
        case(1:10)
          if(.not.Present(dimensions))  &
            call MQC_Error_L('Illegal attempt to initialize MQC_Variable of rank > 0 without defined dimensions.', 6, &
            'Present(dimensions)', Present(dimensions) )
          mqcVariable%dimensions = 0
          mqcVariable%dimensions(1:myRank) = dimensions
        case default
          call MQC_Error_I('Illegal rank sent to MQC_Variable_initialize.', 6, &
               'myRank', myRank)
        end select
      elseIf(myRank.eq.0) then
        mqcVariable%dimensions = 0
      else
        mqcVariable%dimensions(1:myRank) = dimensions
      endIf
!
!     Now, deallocate/reallocate the data arrays if necessary.
!
      myLength = MQC_Variable_getLength(mqcVariable)
      select case(mqcVariable%dataType)
      case('REAL')
        if(Allocated(mqcVariable%integerArray)) DeAllocate(mqcVariable%integerArray)
        if(.not.Allocated(mqcVariable%realArray)) then
          Allocate(mqcVariable%realArray(myLength))
        elseIf(Size(mqcVariable%realArray).ne.myLength) then
          DeAllocate(mqcVariable%realArray)
          Allocate(mqcVariable%realArray(myLength))
        endIf
      case('INTEGER')
        if(Allocated(mqcVariable%realArray)) DeAllocate(mqcVariable%realArray)
        if(.not.Allocated(mqcVariable%integerArray)) then
          Allocate(mqcVariable%integerArray(myLength))
        elseIf(Size(mqcVariable%integerArray).ne.myLength) then
          DeAllocate(mqcVariable%integerArray)
          Allocate(mqcVariable%integerArray(myLength))
        endIf
      case default
        call MQC_Error_A('MQC_Variable_initialize is confused in allocation block.', 6, &
             'mqcVariable%dataType', mqcVariable%dataType)
      end select
!
!     Lastly, set the initialized flag.
!
      mqcVariable%initialized = .True.
!
      return
      end subroutine MQC_Variable_initialize


!
!PROCEDURE MQC_Variable_getArrayPosition
      function MQC_Variable_getArrayPosition(mqcVariable,arrayElement) result(k)
!
!     This function returns an integer giving the position in the object data
!     array for input MQC Variable object <mqcVariable> corresponding to the
!     array index given by input argument <arrayElement>.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64),dimension(:)::arrayElement
      integer(kind=int64)::k
      integer(kind=int64)::i,subLength
!
!
!     Do the work...
!
      select case(TRIM(mqcVariable%storageFormat))
      case('FULL')
        k = 0
        do i=mqcVariable%rank,2,-1
          subLength = SUM(mqcVariable%dimensions(1:i-1))
          k = k + (arrayElement(i)-1)*subLength
        endDo
        k = k + arrayElement(1)
      case default
        call mqc_error_a('MQC_Variable_getArrayPosition: Found an unknown storage format.', 6, &
             'TRIM(mqcVariable%storageFormat)', TRIM(mqcVariable%storageFormat))
      end select
!
      return
      end function MQC_Variable_getArrayPosition


!
!PROCEDURE MQC_Variable_getLength
      function MQC_Variable_getLength(mqcVariable) result(myLen)
!
!     This function returns an integer giving the total number of values that
!     should be stored for an MQC_Variable type matching the properties of input
!     dummy argument <mqcVariable>.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64)::myLen
!
      integer(kind=int64)::i
!
!
      select case(mqcVariable%rank)
      case(0)
        myLen = 1
      case(1:)
        myLen = 1
        do i = 1,mqcVariable%rank
          myLen = myLen*mqcVariable%dimensions(i)
        endDo
      case default
        call MQC_Error_I('Invalid rank found in MQC_Variable_getLength.', 6, &
             'mqcVariable%rank', mqcVariable%rank)
      end select
!
      return
      end function MQC_Variable_getLength


!
!PROCEDURE MQC_Variable_getRank
      function MQC_Variable_getRank(mqcVariable) result(myRank)
!
!     This function returns an integer giving the rank of the MQC_Variable input
!     dummy argument <mqcVariable>.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64)::myRank
!
!
      myRank = mqcVariable%rank
!
      return
      end function MQC_Variable_getRank


!
!PROCEDURE MQC_Variable_getSize
      function MQC_Variable_getSize(mqcVariable,iDimension) result(mySize)
!
!     This function operates on MQC_Variable opjects in an analogous manner to
!     the fortran intrinsic funciotn SIZE. When called, this function will
!     return an integer corresponding to the size of dimension iDimension. If
!     iDimension is NOT sent or is sent <= 0, this function returns the full
!     size of the variable (equivalent to MQC_Variable_getLen, without regard
!     for special storage formats.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64),intent(in),optional::iDimension
      integer(kind=int64)::mySize
!
      integer(kind=int64)::myDimension
!
!
      if(Present(iDimension)) then
        myDimension = iDimension
      else
        myDimension = 0
      endIf
      if(myDimension.gt.MQC_Variable_getRank(mqcVariable))  &
        call MQC_Error_I('MQC_Variable_getSize: Illegal iDimension sent.', 6,&
        'myDimension', myDimension, &
        'MQC_Variable_getRank(mqcVariable)', MQC_Variable_getRank(mqcVariable))
      select case(myDimension)
      case(1:)
        mySize = mqcVariable%dimensions(myDimension)
      case(:0)
        mySize = MQC_Variable_getLength(mqcVariable)
      case default
        call MQC_Error_I('Error in SIZE: Invalid dimension sent.', 6, &
             'myDimension', myDimension )
      end select
!
      return
      end function MQC_Variable_getSize


!
!PROCEDURE MQC_Variable_getType
      function MQC_Variable_getType(mqcVariable) result(myType)
!
!     This function returns a character string that indicates the data type of
!     an MQC_Variable object. Possible return values include:
!
!           'UNKNOWN'   ...   Indicates an unknown data type has been found.
!           'REAL'      ...   Indicates the type of mqcVariable is real.
!           'INTEGER'   ...   Indicates the type of mqcVariable is integer.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      character(len=64)::myType
!
!
!     Determine the right value for OUTPUT argument <myType>. Then, ensure the
!     resulting string is left-justified.
!
      select case(mqcVariable%dataType)
      case('REAL')
        myType = 'REAL'
      case('INTEGER')
        myType = 'INTEGER'
      case default
        myType = 'UNKNOWN'
      end select
      myType = ADJUSTL(myType)
!
      return
      end function MQC_Variable_getType


!
!PROCEDURE MQC_Variable_getTypeCode
      function MQC_Variable_getTypeCode(mqcVariable) result(myTypeCode)
!
!     This function returns an integer code that corresponds to a variable type.
!     This function is only used in this module to make type mixing/conversion
!     during operations consistent throughout the module.
!
!     NOTE: The type codes returned by this function are ordered by actual
!     underlying typing such that lower valued codes indicated higher-precedent
!     types. In other words, type-code 3 can always be treated as a type-code 2
!     value, but the reverse is not generally true.
!
!     The output values are:
!           0 ... UNKNOWN
!           2 ... REAL
!           3 ... INTEGER
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64)::myTypeCode
!
!
!     Determine the right value for OUTPUT argument <myTypeCode>.
!
      select case(mqcVariable%dataType)
      case('REAL')
        myTypeCode = 2
      case('INTEGER')
        myTypeCode = 3
      case default
        myTypeCode = 0
      end select
!
      return
      end function MQC_Variable_getTypeCode


!
!PROCEDURE MQC_Variable_isConformable
      function MQC_Variable_isConformable(mqcVariable1,mqcVariable2) result(conformable)
!
!     This function returns a logical flag indicating is input dummy arguments
!     <mqcVariable1> and <mqcVariable2> conform. This routine checks that rank
!     and dimension lengths of the two input dummy arguments match.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      logical::conformable
      integer(kind=int64)::i
!
!
!     Initialize <conformable> and then check ranks followed by dimension
!     lengths.
!
      conformable = mqcVariable1%getRank().eq.mqcVariable2%getRank()
      if(.not.conformable) return
      do i = 1,mqcVariable1%getRank()
        conformable = mqcVariable1%dimensions(i).eq.mqcVariable2%dimensions(i)
        if(.not.conformable) return
      endDo
!
      return
      end function MQC_Variable_isConformable


!
!PROCEDURE MQC_Variable_isInitialized
      function MQC_Variable_isInitialized(mqcVariable) result(myIsInitialized)
!
!     This function returns a logical value indicated whether input dummy
!     argument <mqcVariable> has been initialized (output is .true.) or not
!     (output is .false.).
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      logical::myIsInitialized
!
!
      myIsInitialized = mqcVariable%initialized
!
      return
      end function MQC_Variable_isInitialized


!
!----------------------------------------------------------------
!                                                               |
!     MQC_Variable PRINTING PROCEDURES                          |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE MQC_Print_mqcVariable
      Subroutine MQC_Print_mqcVariable(mqcVariable,iOut,header,  &
        blankAtTop,blankAtBottom,formatString)
!
!     This subroutine is used to print an MQC variable. Dummy arguments <iOut>,
!     <header>, <blankAtTop>, <blankAtBottom>, and <formatString> are OPTIONAL
!     input dummy arguments.
!
!     Dummy argument definitions (all dummy arguments are INPUT):
!           mqcVariable       This is the MQC variable object to be printed.
!           iOut              This is the output unit number. By default, this
!                             argument is taken to be unit 6.
!
!
!     L. M. Thompson, 2016.
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),optional,intent(In)::iOut
      character(len=*),intent(in),optional::header,formatString
      logical,intent(in),optional::blankAtTop,blankAtBottom
!
      integer(kind=int64)::myIOut
      logical::myBlankAtTop,myBlankAtBottom
!
!
!     Set-up local versions of optional dummy arguments.
!
      myIOut = 6
      myBlankAtTop = .false.
      myBlankAtBottom = .false.
      if(PRESENT(iOut)) myIOut = iOut
      if(PRESENT(blankAtTop)) myBlankAtTop = blankAtTop
      if(PRESENT(blankAtBottom)) myBlankAtBottom = blankAtBottom
!
!     Determine if we need to print a scalar or array, then call the underlying
!     MQC printing routine.
!
      select case(mqcVariable%getRank())
!
!     Print scalars...
      case(0)
        select case(mqcVariable%getType())
        case('REAL')
          if(PRESENT(header)) then
            call mqc_print_scalar(myIOut,mqcVariable%realArray(1),Header,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_scalar(myIOut,mqcVariable%realArray(1),  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case('INTEGER')
          if(PRESENT(header)) then
            call mqc_print_scalar(myIOut,mqcVariable%integerArray(1),Header,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_scalar(myIOut,mqcVariable%integerArray(1),  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case default
          call mqc_error_a('MQC variable print requested for unknown type.', 6, &
          'mqcVariable%getType()', mqcVariable%getType())
        end select
!
!     Print vectors...
      case(1)
        select case(mqcVariable%getType())
        case('REAL')
          if(PRESENT(header)) then
            call mqc_print_vector(myIOut,mqcVariable%realArray(:),Header,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_vector(myIOut,mqcVariable%realArray(:),  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case('INTEGER')
          if(PRESENT(header)) then
            call mqc_print_vector(myIOut,mqcVariable%integerArray(:),Header,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_vector(myIOut,mqcVariable%integerArray(:),  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case default
          call mqc_error_a('MQC variable print requested for unknown type.', 6, &
               'mqcVariable%getType()', mqcVariable%getType() )
        end select
!
!     Print matrices...
      case(2)
        select case(mqcVariable%getType())
        case('REAL')
          if(PRESENT(header)) then
            call mqc_print_matrix(myIOut,  &
              RESHAPE(mqcVariable%realArray(:),mqcVariable%dimensions(1:2)),  &
              Header,Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_matrix(myIOut,  &
              RESHAPE(mqcVariable%realArray,mqcVariable%dimensions(1:2)),  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case('INTEGER')
          if(PRESENT(header)) then
            call mqc_print_matrix(myIOut,  &
              RESHAPE(mqcVariable%integerArray,mqcVariable%dimensions(1:2)),  &
              Header,Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_matrix(myIOut,  &
              RESHAPE(mqcVariable%integerArray,mqcVariable%dimensions(1:2)),  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case default
          call mqc_error_a('MQC variable print requested for unknown type.', 6, &
               'mqcVariable%getType()', mqcVariable%getType())
        end select
      case default
        call mqc_error_i('MQC variable print requested for unknown rank.', 6, &
             'mqcVariable%getRank()', mqcVariable%getRank())
      end select
!
      return
      end subroutine MQC_Print_mqcVariable


!
!----------------------------------------------------------------
!                                                               |
!     MQC_Variable VALUE SETTING & BASIC MANIPULATION PROCEDURES|
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE MQC_Variable_put_MQC
      subroutine MQC_Variable_put_MQC(mqcVariable,valueMQC,arrayElement)
!
!     This subroutine is used to put a value into a specified element of an
!     MQC_Variable array. The value is sent as an MQC variable type in dummy
!     argument <valueMQC> and the array element is sent in dummy argument
!     <arrayElement>.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      class(MQC_Variable),intent(in)::valueMQC
      integer(kind=int64),dimension(:)::arrayElement
      integer(kind=int64)::k
!
!
!     Do the work...
!
      if(valueMQC%rank.ne.0) call mqc_error_i('MQC_Variable_Put: Sent MQC value must be a scalar.', 6, &
      'valueMQC%rank', valueMQC%rank)
      if(.not.mqcVariable%initialized) call mqc_error_l('Cannot use MQC_Variable_Put with uninitialized MQC variable.', 6, &
           'mqcVariable%initialized', mqcVariable%initialized )
      if(mqcVariable%rank.ne.SIZE(arrayElement))  &
        call mqc_error_i('MQC_Variable_Put: Rank of variable does NOT match number of array indices provided.', 6, &
        'mqcVariable%rank', mqcVariable%rank, &
        'SIZE(arrayElement)', SIZE(arrayElement))
      k = MQC_Variable_getArrayPosition(mqcVariable,arrayElement)
      if(k.gt.MQC_Variable_getLength(mqcVariable)) call mqc_error_i('MQC_Variable_Put: Invalid array index provided.', 6, &
           'k', k, &
           'MQC_Variable_getLength(mqcVariable)', MQC_Variable_getLength(mqcVariable))
      select case(TRIM(MQC_Variable_getType(mqcVariable)))
      case('REAL')
        mqcVariable%realArray(k) = INT(valueMQC)
      case('INTEGER')
        mqcVariable%integerArray(k) = FLOAT(valueMQC)
      case default
        call mqc_error_a('MQC_Variable_Put: Unknown MQC variable type.', 6, &
             'TRIM(MQC_Variable_getType(mqcVariable))', TRIM(MQC_Variable_getType(mqcVariable)) )
      end select
!
      return
      end subroutine MQC_Variable_put_MQC


!
!PROCEDURE MQC_Variable_put_intrinsicInteger
      subroutine MQC_Variable_put_intrinsicInteger(mqcVariable,valueInteger,arrayElement)
!
!     This subroutine is used to put a value into a specified element of an
!     MQC_Variable array. The value is sent as an intrinsic integer in dummy
!     argument <valueInteger> and the array element is sent in dummy argument
!     <arrayElement>.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),intent(in)::valueInteger
      integer(kind=int64),dimension(:)::arrayElement
      integer(kind=int64)::k
!
!
!     Do the work...
!
      if(.not.mqcVariable%initialized) call mqc_error_l('Cannot use MQC_Variable_Put with uninitialized MQC variable.', 6, &
           'mqcVariable%initialized', mqcVariable%initialized )
      if(mqcVariable%rank.ne.SIZE(arrayElement))  &
        call mqc_error_I('MQC_Variable_Put: Rank of variable does NOT match number of array indices provided.', 6, &
        'mqcVariable%rank', mqcVariable%rank, &
        'SIZE(arrayElement)', SIZE(arrayElement) )
      k = MQC_Variable_getArrayPosition(mqcVariable,arrayElement)
      if(k.gt.MQC_Variable_getLength(mqcVariable)) call mqc_error_i('MQC_Variable_Put: Invalid array index provided.', 6, &
           'k', k, &
           'MQC_Variable_getLength(mqcVariable)', MQC_Variable_getLength(mqcVariable))
      select case(TRIM(MQC_Variable_getType(mqcVariable)))
      case('REAL')
        mqcVariable%realArray(k) = valueInteger
      case('INTEGER')
        mqcVariable%integerArray(k) = valueInteger
      case default
        call mqc_error_A('MQC_Variable_Put: Unknown MQC variable type.', 6, &
             'TRIM(MQC_Variable_getType(mqcVariable))', TRIM(MQC_Variable_getType(mqcVariable)))
      end select
!
      return
      end subroutine MQC_Variable_put_intrinsicInteger


!
!PROCEDURE MQC_Variable_put_intrinsicReal
      subroutine MQC_Variable_put_intrinsicReal(mqcVariable,valueReal,arrayElement)
!
!     This subroutine is used to put a value into a specified element of an
!     MQC_Variable array. The value is sent as an intrinsic real in dummy
!     argument <valueReal> and the array element is sent in dummy argument
!     <arrayElement>.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      real(kind=real64),intent(in)::valueReal
      integer(kind=int64),dimension(:)::arrayElement
      integer(kind=int64)::k
!
!
!     Do the work...
!
      if(.not.mqcVariable%initialized) call mqc_error_L('Cannot use MQC_Variable_Put with uninitialized MQC variable.', 6, &
           'mqcVariable%initialized', mqcVariable%initialized )
      if(mqcVariable%rank.ne.SIZE(arrayElement))  &
        call mqc_error_i('MQC_Variable_Put: Rank of variable does NOT match number of array indices provided.', 6, &
        'mqcVariable%rank', mqcVariable%rank, &
        'SIZE(arrayElement)', SIZE(arrayElement) )
      k = MQC_Variable_getArrayPosition(mqcVariable,arrayElement)
      if(k.gt.MQC_Variable_getLength(mqcVariable)) call mqc_error_I('MQC_Variable_Put: Invalid array index provided.', 6, &
           'MQC_Variable_getLength(mqcVariable))', MQC_Variable_getLength(mqcVariable))
      select case(TRIM(MQC_Variable_getType(mqcVariable)))
      case('REAL')
        mqcVariable%realArray(k) = valueReal
      case('INTEGER')
        mqcVariable%integerArray(k) = valueReal
      case default
        call mqc_error_a('MQC_Variable_Put: Unknown MQC variable type.', 6, &
             'TRIM(MQC_Variable_getType(mqcVariable))', TRIM(MQC_Variable_getType(mqcVariable)))
      end select
!
      return
      end subroutine MQC_Variable_put_intrinsicReal


!
!PROCEDURE MQC_Variable_clear_general
      subroutine MQC_Variable_clear_general(mqcVariable,valueInteger,  &
        valueReal,dimensions)
!
!     This subroutine is used to clear and initialize the MQC variable
!     <mqcVariable>. Dummy arguments <valueInteger>, <valueReal>, and
!     <dimensions> are optional.
!
!     When this routine is called, <mqcVariable> does not already need to have
!     been initialized. All elements of mqcVariable are set equal to
!     <valueInteger> or <valueReal>, if sent. Otherwise, all elements of
!     <mqcVariable> are set equal to 0.0. If <dimensions> is sent, then
!     <mqcVariable> will be set (or re-set) to have those dimensions. Otherwise,
!     <mqcVariable> will be set as a scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),intent(in),optional::valueInteger
      real(kind=real64),intent(in),optional::valueReal
      integer(kind=int64),dimension(:),intent(in),optional::dimensions
!
      integer(kind=int64)::nInputs,myRank
      integer(kind=int64),dimension(:),allocatable::myDimensions
      logical::fillInteger,fillReal
      character(len=64)::myDataType
!
!
!     Begin by figuring out if we are working with integer or real numbers.
!
      nInputs = 0
      if(Present(valueInteger)) nInputs = nInputs+1
      if(Present(valueReal)) nInputs = nInputs+1
      if(nInputs.gt.1)  &
        call MQC_Error_I('MQC_Variable_clear Error: Multiple input arguments sent.', 6, &
        'nInputs', nInputs)
      if(nInputs.lt.1)  &
        call MQC_Error_I('MQC_Variable_clear Error: No input values sent.', 6, &
        'nInputs', nInputs)
      fillInteger = Present(valueInteger)
      fillReal = Present(valueReal)
!
!     Set myDataType, myRank, and myDimensions.
!
      if(fillInteger) then
        myDataType = 'INTEGER'
      elseIf(fillReal) then
        myDataType = 'REAL'
      else
        call MQC_Error_L('MQC_Variable_clear: Confused setting myDataType.', 6, &
             'fillInteger', fillInteger, &
             'fillReal', fillReal )
      endIf
      if(PRESENT(dimensions)) then
        myRank = SIZE(dimensions)
        allocate(myDimensions(myRank))
        myDimensions = dimensions
      else
        myRank = 0
      endIf
!
!     Initialize mqcVariable and then set the value appropriately.
!
      if(myRank.eq.0) then
        call MQC_Variable_initialize(mqcVariable,myDataType)
      else
        call MQC_Variable_initialize(mqcVariable,myDataType,myDimensions,  &
          myRank,'FULL')
      endIf
      if(fillInteger) then
        mqcVariable%integerArray(:) = valueInteger
      elseIf(fillReal) then
        mqcVariable%realArray(:) = valueReal
      else
        call MQC_Error_L('MQC_Variable_clear: Confused filling input value.', 6, &
             'fillInteger', fillInteger, &
             'fillReal', fillReal )
      endIf
!
      return
      end subroutine MQC_Variable_clear_general


!
!PROCEDURE MQC_Variable_clear_mqc
      subroutine MQC_Variable_clear_mqc(mqcVariable,valueMQC,dimensions)
!
!     This subroutine is a wrapper for Procedure MQC_Variable_clear_general,
!     which is used to clear and initialize the MQC variable <mqcVariable>.
!     Dummy argument <dimensions> is optional.
!
!     When this routine is called, <mqcVariable> does not already need to have
!     been initialized. All elements of mqcVariable are set equal to
!     <valueInteger>. If <dimensions> is sent, then <mqcVariable> will be set
!     (or re-set) to have those dimensions. Otherwise, <mqcVariable> will be set
!     as a scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      class(MQC_Variable),intent(in)::valueMQC
      integer(kind=int64),dimension(:),intent(in),optional::dimensions
!
!     Call Procedure MQC_Variable_clear_general.
!
      if(valueMQC%rank.ne.0) call mqc_error_i('MQC_Variable_clear: Sent MQC value must be a scalar.', 6, &
           'valueMQC%rank', valueMQC%rank )
      select case(TRIM(valueMQC%getType()))
      case('REAL')
        if(PRESENT(dimensions)) then
          call MQC_Variable_clear_general(mqcVariable,  &
            valueReal=valueMQC%realArray(1),dimensions=dimensions)
        else
          call MQC_Variable_clear_general(mqcVariable,  &
            valueReal=valueMQC%realArray(1))
        endIf
      case('INTEGER')
        if(PRESENT(dimensions)) then
          call MQC_Variable_clear_general(mqcVariable,  &
            valueInteger=valueMQC%integerArray(1),dimensions=dimensions)
        else
          call MQC_Variable_clear_general(mqcVariable,  &
            valueInteger=valueMQC%integerArray(1))
        endIf
      case default
        call mqc_error_A('MQC_Variable_clear: Unknown type of MQC variable sent.', 6, &
             'TRIM(valueMQC%getType())', TRIM(valueMQC%getType()))
      end select
!
      return
      end subroutine MQC_Variable_clear_mqc


!
!PROCEDURE MQC_Variable_clear_integer
      subroutine MQC_Variable_clear_integer(mqcVariable,valueInteger,dimensions)
!
!     This subroutine is a wrapper for Procedure MQC_Variable_clear_general,
!     which is used to clear and initialize the MQC variable <mqcVariable>.
!     Dummy arguments <valueInteger> and <dimensions> are optional.
!
!     When this routine is called, <mqcVariable> does not already need to have
!     been initialized. All elements of mqcVariable are set equal to
!     <valueInteger>. If <dimensions> is sent, then <mqcVariable> will be set
!     (or re-set) to have those dimensions. Otherwise, <mqcVariable> will be set
!     as a scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),intent(in)::valueInteger
      integer(kind=int64),dimension(:),intent(in),optional::dimensions
!
!     Call Procedure MQC_Variable_clear_general.
!
      if(PRESENT(dimensions)) then
        call MQC_Variable_clear_general(mqcVariable,valueInteger=valueInteger,  &
          dimensions=dimensions)
      else
        call MQC_Variable_clear_general(mqcVariable,valueInteger=valueInteger)
      endIf
!
      return
      end subroutine MQC_Variable_clear_integer


!
!PROCEDURE MQC_Variable_clear_real
      subroutine MQC_Variable_clear_real(mqcVariable,valueReal,dimensions)
!
!     This subroutine is a wrapper for Procedure MQC_Variable_clear_general,
!     which is used to clear and initialize the MQC variable <mqcVariable>.
!     Dummy arguments <valueReal> and <dimensions> are optional.
!
!     When this routine is called, <mqcVariable> does not already need to have
!     been initialized. All elements of mqcVariable are set equal to
!     <valueReal>. If <dimensions> is sent, then <mqcVariable> will be set
!     (or re-set) to have those dimensions. Otherwise, <mqcVariable> will be set
!     as a scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      real(kind=real64),intent(in)::valueReal
      integer(kind=int64),dimension(:),intent(in),optional::dimensions
!
!     Call Procedure MQC_Variable_clear_general.
!
      if(PRESENT(dimensions)) then
        call MQC_Variable_clear_general(mqcVariable,valueReal=valueReal,  &
          dimensions=dimensions)
      else
        call MQC_Variable_clear_general(mqcVariable,valueReal=valueReal)
      endIf
!
      return
      end subroutine MQC_Variable_clear_real


!
!PROCEDURE MQC_Variable_mqc2mqc
      subroutine MQC_Variable_mqc2mqc(mqcVariable,mqcIn)
!
!     This subroutine is used to set an MQC Variable equal to another MQC
!     Variable.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      class(MQC_Variable),intent(in)::mqcIn
!
!
!     Do the work...
!
      if(Allocated(mqcVariable%realArray)) DeAllocate(mqcVariable%realArray)
      if(Allocated(mqcVariable%integerArray)) DeAllocate(mqcVariable%integerArray)
      mqcVariable%rank = mqcIn%rank
      mqcVariable%dataType = mqcIn%dataType
      mqcVariable%storageFormat = mqcIn%storageFormat
      mqcVariable%dimensions = mqcIn%dimensions
      select case(TRIM(mqcIn%getType()))
      case('REAL')
        Allocate(mqcVariable%realArray(SIZE(mqcIn%realArray(:))))
        mqcVariable%realArray = mqcIn%realArray
      case('INTEGER')
        Allocate(mqcVariable%integerArray(SIZE(mqcIn%integerArray(:))))
        mqcVariable%integerArray = mqcIn%integerArray
      case default
        call mqc_error_A('MQC_Variable_mqc2mqc: Confused by data type.', 6, &
             'TRIM(mqcIn%getType())', TRIM(mqcIn%getType()) )
      end select
      mqcVariable%initialized = .true.
!
      return
      end subroutine MQC_Variable_mqc2mqc


!
!PROCEDURE MQC_Variable_setVal
      subroutine MQC_Variable_setVal(mqcVariable,scalarIntegerIn,  &
        scalarRealIn,arrayIntegerIn,arrayRealIn,dimensions,DEBUGPRINT)
!
!     This subroutine is used to set one or more values in the packed storage
!     array of MQC_Variable dummy argument <mqcVariable> to an scalar or array
!     given by one of <scalarRealIn>, <scalarIntegerIn>, <arrayRealIn>, or
!     <arrayIntegerIn>.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),intent(in),optional::scalarIntegerIn
      real(kind=real64),intent(in),optional::scalarRealIn
      integer(kind=int64),dimension(:),optional::arrayIntegerIn
      real(kind=real64),dimension(:),optional::arrayRealIn
      integer(kind=int64),dimension(:),optional::dimensions
      logical,optional::DEBUGPRINT
!
      integer(kind=int64)::nInputs,myRank
      integer(kind=int64),dimension(:),allocatable::myDimensions
      logical::DEBUG,fillInteger,fillReal,fillScalar,fillArray
      character(len=64)::myDataType
!
!
!     Begin by setting DEBUG and then figuring out if we are working with
!     integer/real and scalar/array. Make sure only one of the four
!     possibilities has been sent.
!
      if(Present(DEBUGPRINT)) then
        DEBUG = DEBUGPRINT
      else
        DEBUG = .False.
      endIf
      if(DEBUG) then
        write(*,*)
        write(*,*)' ---------------------------------------------------'
        write(*,*)' DEBUG PRINT FOR MQC_Variable_setVal'
        write(*,*)' ---------------------------------------------------'
        write(*,*)' INPUT ARGUMENTS:'
        write(*,*)'   PRESENT(scalarIntegerIn) = ',PRESENT(scalarIntegerIn)
        write(*,*)'   PRESENT(scalarRealIn)    = ',PRESENT(scalarRealIn)
        write(*,*)'   PRESENT(arrayIntegerIn)  = ',PRESENT(arrayIntegerIn)
        write(*,*)'   PRESENT(arrayRealIn)     = ',PRESENT(arrayRealIn)
        write(*,*)'   PRESENT(dimensions)      = ',PRESENT(dimensions)
        write(*,*)
      endIf
      nInputs = 0
      if(Present(scalarIntegerIn)) nInputs = nInputs+1
      if(Present(scalarRealIn)) nInputs = nInputs+1
      if(Present(arrayIntegerIn)) nInputs = nInputs+1
      if(Present(arrayRealIn)) nInputs = nInputs+1
      if(nInputs.gt.1)  &
        call MQC_Error_I('MQC_Variable_setVal Error: Multiple input values sent.', 6, &
        "nInputs", nInputs )
      if(nInputs.lt.1)  &
        call MQC_Error_I('MQC_Variable_setVal Error: No input values sent.', 6, &
        "nInputs", nInputs )
      fillInteger = Present(scalarIntegerIn).or.Present(arrayIntegerIn)
      fillReal = Present(scalarRealIn).or.Present(arrayRealIn)
      fillScalar = Present(scalarIntegerIn).or.Present(scalarRealIn)
      fillArray = Present(arrayIntegerIn).or.Present(arrayRealIn)
      if(DEBUG) then
        write(*,*)'   fillInteger=',fillInteger
        write(*,*)'   fillReal   =',fillReal
        write(*,*)'   fillScalar =',fillScalar
        write(*,*)'   fillArray  =',fillArray
      endIf
!
!     Set myDataType and myRank.
!
      if(fillInteger) then
        myDataType = 'INTEGER'
      elseIf(fillReal) then
        myDataType = 'REAL'
      else
        call MQC_Error_L('MQC_Variable_setVal: Confused setting myDataType.', 6, &
             'fillInteger', fillInteger, &
             'fillReal', fillReal )
      endIf
      if(DEBUG) write(*,*)'   myDataType = ',myDataType
      if(PRESENT(dimensions)) then
        myRank = SIZE(dimensions)
        if(myRank.eq.1.and.fillScalar) myRank = 0
        if(myRank.gt.0) then
          Allocate(myDimensions(myRank))
          myDimensions = dimensions
        endIf
      else
        if(fillScalar) then
          myRank = 0
        elseIf(fillArray) then
          myRank = 1
          Allocate(myDimensions(myRank))
          if(fillInteger) then
            myDimensions(1) = SIZE(arrayIntegerIn)
          elseIf(fillReal) then
            myDimensions(1) = SIZE(arrayRealIn)
          else
            call mqc_error_l('MQC_Variable_setVal: Confused setting myDimensions.', 6, &
              'fillScalar', fillScalar, &
              'fillArray', fillArray, &
              'fillInteger', fillInteger, &
              'fillReal', fillReal )
          endIf
        else
          call MQC_Error_L('MQC_Variable_setVal: Confused setting myRank.', 6, &
              'fillScalar', fillScalar, &
              'fillArray', fillArray )
        endIf
      endIf
      if(DEBUG) write(*,*)'   myRank = ',myRank
!
!     Initialize mqcVariable and then set the value appropriately.
!
      if(DEBUG) write(*,*)'   Calling Init...'
      if(myRank.eq.0) then
        call MQC_Variable_initialize(mqcVariable,myDataType)
      else
        call MQC_Variable_initialize(mqcVariable,myDataType,myDimensions,  &
          myRank,'FULL')
      endIf
      if(DEBUG) then
        write(*,*)'   Back from init.'
        write(*,*)'   rank = ',mqcVariable%rank
        write(*,*)'   dimensions = ',mqcVariable%dimensions
        write(*,*)'   ALLOCATED(realArray) ?    = ',ALLOCATED(mqcVariable%realArray)
        write(*,*)'   ALLOCATED(integerArray) ? = ',ALLOCATED(mqcVariable%integerArray)
        if(Allocated(mqcVariable%realArray)) write(*,*)'   SIZE(realArray) = ',SIZE(mqcVariable%realArray)
        if(Allocated(mqcVariable%integerArray)) write(*,*)'   SIZE(integerArray) = ',SIZE(mqcVariable%integerArray)
      endIf
      if(fillScalar) then
        if(fillInteger) then
          mqcVariable%integerArray(1) = scalarIntegerIn
        elseIf(fillReal) then
          mqcVariable%realArray(1) = scalarRealIn
        else
          call MQC_Error_L('MQC_Variable_setVal: Confused filling scalar value.', 6, &
               'fillScalar', fillScalar, &
               'fillInteger', fillInteger, &
               'fillReal', fillReal )
        endIf
      elseIf(fillArray) then
        if(fillInteger) then
          mqcVariable%integerArray(:) = arrayIntegerIn
        elseIf(fillReal) then
          mqcVariable%realArray(:) = arrayRealIn
        else
          call MQC_Error_L('MQC_Variable_setVal: Confused filling array values.', 6, &
               'fillArray', fillArray, &
               'fillInteger', fillInteger, &
               'fillReal', fillReal )
        endIf
      endIf
!
      if(DEBUG) then
        write(*,*)' ---------------------------------------------------'
        write(*,*)' ---------------------------------------------------'
        write(*,*)
        write(*,*)
      endIf
      return
      end subroutine MQC_Variable_setVal


!
!PROCEDURE MQC_Variable_setVal_intrinsicScalar2mqcScalar_integer
      subroutine MQC_Variable_setVal_intrinsicScalar2mqcScalar_integer(mqcVariable,intrinsicIn)
!
!     This subroutine is used to set an MQC Variable scalar equal to an
!     intrinsic integer fortran scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),intent(in)::intrinsicIn
!
!
!     Load the MQC variable.
!
      call MQC_Variable_setVal(mqcVariable,scalarIntegerIn=intrinsicIn)
!
      return
      end subroutine MQC_Variable_setVal_intrinsicScalar2mqcScalar_integer


!
!PROCEDURE MQC_Variable_setVal_intrinsicScalar2mqcScalar_real
      subroutine MQC_Variable_setVal_intrinsicScalar2mqcScalar_real(mqcVariable,intrinsicIn)
!
!     This subroutine is used to set an MQC Variable scalar equal to an
!     intrinsic real fortran scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      real(kind=real64),intent(in)::intrinsicIn
!
!
!     Load the MQC variable.
!
      call MQC_Variable_setVal(mqcVariable,scalarRealIn=intrinsicIn)
!
      return
      end subroutine MQC_Variable_setVal_intrinsicScalar2mqcScalar_real


!
!PROCEDURE MQC_Variable_setVal_intrinsicVector2mqcVector_integer
      subroutine MQC_Variable_setVal_intrinsicVector2mqcVector_integer(mqcVariable,intrinsicIn)
!
!     This subroutine is used to set an MQC Variable vector equal to an
!     intrinsic integer fortran vector.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),dimension(:),intent(in)::intrinsicIn
!
!
!     Load the MQC variable.
!
      call MQC_Variable_setVal(mqcVariable,arrayIntegerIn=intrinsicIn)
!
      return
      end subroutine MQC_Variable_setVal_intrinsicVector2mqcVector_integer


!
!PROCEDURE MQC_Variable_setVal_intrinsicVector2mqcVector_real
      subroutine MQC_Variable_setVal_intrinsicVector2mqcVector_real(mqcVariable,intrinsicIn)
!
!     This subroutine is used to set an MQC Variable vector equal to an
!     intrinsic real fortran vector.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      real(kind=real64),dimension(:),intent(in)::intrinsicIn
!
!
!     Begin by figuring out if we are working with integer/real and scalar/
!     array. Make sure only one of the four possibilities has been sent.
!
      call MQC_Variable_setVal(mqcVariable,arrayRealIn=intrinsicIn)
!
      return
      end subroutine MQC_Variable_setVal_intrinsicVector2mqcVector_real


!
!PROCEDURE MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_integer
      subroutine MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_integer(mqcVariable,intrinsicIn)
!
!     This subroutine is used to set an MQC Variable matrix (full storage form)
!     equal to an intrinsic integer fortran matrix.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),dimension(:,:),intent(in)::intrinsicIn
      integer(kind=int64)::myLength
      integer(kind=int64),dimension(2)::myDimensions
!
!
!     Load the MQC variable.
!
      myLength = SIZE(intrinsicIn)
      myDimensions(1) = SIZE(intrinsicIn,1)
      myDimensions(2) = SIZE(intrinsicIn,2)
      write(*,*)
      write(*,*)' hrant - new routine:'
      write(*,*)' myLength = ',myLength
      call mqc_print(6,myDimensions,'myDimensions:')

      write(*,*)
      write(*,*)' HRANT -- Calling SetVal routine....'
      write(*,*)
      call MQC_Variable_setVal(mqcVariable,arrayIntegerIn=RESHAPE(intrinsicIn,[myLength]),  &
        dimensions=myDimensions,DEBUGPRINT=.False.)
!
      return
      end subroutine MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_integer


!
!PROCEDURE MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_real
      subroutine MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_real(mqcVariable,intrinsicIn)
!
!     This subroutine is used to set an MQC Variable matrix (full storage form)
!     equal to an intrinsic integer fortran matrix.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      real(kind=real64),dimension(:,:),intent(in)::intrinsicIn
      integer(kind=int64)::myLength
      integer(kind=int64),dimension(2)::myDimensions
!
!
!     Load the MQC variable.
!
      myLength = SIZE(intrinsicIn)
      myDimensions(1) = SIZE(intrinsicIn,1)
      myDimensions(2) = SIZE(intrinsicIn,2)
      call MQC_Variable_setVal(mqcVariable,arrayRealIn=RESHAPE(intrinsicIn,[myLength]),  &
        dimensions=myDimensions,DEBUGPRINT=.False.)
!
      return
      end subroutine MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_real


!
!PROCEDURE MQC_Variable_mqc2intrinsicIntegerScalar
      subroutine MQC_Variable_mqc2intrinsicIntegerScalar(intrinsicOut,mqcVariable)
!
!     This subroutine is used to set an (output) intrinsic integer scalar equal
!     to an (input) MQC Variable scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64)::intrinsicOut
      class(MQC_Variable),intent(in)::mqcVariable
!
!
!     Do the work...
!
      if(Rank(mqcVariable).ne.0)  &
        call mqc_error_i('Attemp to convert MQC_Variable array to intrinsic SCALAR not allowed.', 6, &
        'Rank(mqcVariable)', Rank(mqcVariable) )
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        intrinsicOut = INT(mqcVariable%realArray(1))
      case(3)
        intrinsicOut = mqcVariable%integerArray(1)
      case default
        call mqc_error_I('MQC_Variable_mqc2intrinsicIntegerScalar: Unknown MQC_Variable type found.', 6, &
             'MQC_Variable_getTypeCode(mqcVariable)', MQC_Variable_getTypeCode(mqcVariable) )
      end select
!
      return
      end subroutine MQC_Variable_mqc2intrinsicIntegerScalar


!
!PROCEDURE MQC_Variable_mqc2intrinsicRealScalar
      subroutine MQC_Variable_mqc2intrinsicRealScalar(intrinsicOut,mqcVariable)
!
!     This subroutine is used to set an (output) intrinsic real scalar equal to
!     an (input) MQC Variable scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),intent(inOut)::intrinsicOut
      type(MQC_Variable),intent(in)::mqcVariable
!
!
!     Do the work...
!
      if(Rank(mqcVariable).ne.0)  &
        call mqc_error_I('Attemp to convert MQC_Variable array to intrinsic SCALAR not allowed.', 6, &
        'Rank(mqcVariable)', Rank(mqcVariable) )
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        intrinsicOut = mqcVariable%realArray(1)
      case(3)
        intrinsicOut = mqcVariable%integerArray(1)
      case default
        call mqc_error_I('MQC_Variable_mqc2intrinsicRealScalar: Unknown MQC_Variable type found.', 6, &
             'MQC_Variable_getTypeCode(mqcVariable)', MQC_Variable_getTypeCode(mqcVariable) )
      end select
!
      return
      end subroutine MQC_Variable_mqc2intrinsicRealScalar


!
!PROCEDURE MQC_Variable_from_intrinsicInteger
      function MQC_Variable_from_intrinsicInteger(intrinsicIn) result(mqcVariable)
!
!     This function is used to convert an intrinsic integer scalar to an MQC
!     Variable.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),intent(in)::intrinsicIn
      type(MQC_Variable)::mqcVariable
!
!
!     Load the MQC variable.
!
      call MQC_Variable_setVal(mqcVariable,scalarIntegerIn=intrinsicIn)
!
      return
      end function MQC_Variable_from_intrinsicInteger


!
!PROCEDURE MQC_Variable_from_intrinsicReal
      function MQC_Variable_from_intrinsicReal(intrinsicIn) result(mqcVariable)
!
!     This subroutine is used to set an MQC Variable scalar equal to an
!     intrinsic real fortran scalar.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),intent(in)::intrinsicIn
      type(MQC_Variable)::mqcVariable
!
!
!     Load the MQC variable.
!
      call MQC_Variable_setVal(mqcVariable,scalarRealIn=intrinsicIn)
!
      return
      end function MQC_Variable_from_intrinsicReal


!
!PROCEDURE MQC_Variable_INT_Scalar
      function MQC_Variable_INT_Scalar(mqcVariable) result(intrinsicOut)
!
!     This function provides INT as a function for MQC Variables.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64)::intrinsicOut
!
!
!     Do the work...
!
      call MQC_Variable_mqc2intrinsicIntegerScalar(intrinsicOut,mqcVariable)
!
      return
      end function MQC_Variable_INT_Scalar


!
!PROCEDURE MQC_Variable_FLOAT_Scalar
      function MQC_Variable_FLOAT_Scalar(mqcVariable) result(intrinsicOut)
!
!     This function provides FLOAT as a function for MQC Variables.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      real(kind=real64)::intrinsicOut
!
!
!     Do the work...
!
      call MQC_Variable_mqc2intrinsicRealScalar(intrinsicOut,mqcVariable)
!
      return
      end function MQC_Variable_FLOAT_Scalar









!hph+
!!
!!PROCEDURE MQC_Variable_getArrayIndex
!      function MQC_Variable_getArrayIndex(mqcVariable,fullStorageIndices)  &
!        result(linearIndex)
!!
!!     This function is used to convert an input full-storage set of array
!!     indices into the linear'ized index where the actual number of interest
!!     resides in mqcVariable%realArray or mqcVariable%integerArray.
!!
!!
!!     H. P. Hratchian, 2017.
!!
!!
!!     Variable Declarations.
!      implicit none
!      class(MQC_Variable)::mqcVariable
!      integer,dimension(:),intent(in)::fullStorageIndices
!!
!
!      
!      
!      
!      
!      select case(mqcVariable%rank)
!      case(0)
!        myLen = 1
!      case(1:)
!        myLen = 1
!        do i = 1,mqcVariable%rank
!          myLen = myLen*mqcVariable%dimensions(i)
!        endDo
!      case default
!        call MQC_Error('Invalid rank found in MQC_Variable_getLength.')
!      end select
!      
!!
!!
!!     Begin by figuring out if we are working with integer/real and scalar/
!!     array. Make sure only one of the four possibilities has been sent.
!!
!      nInputs = 0
!!
!      return
!      end subroutine MQC_Variable_getVal
!hph-



!
!----------------------------------------------------------------
!                                                               |
!     MQC_Variable MATHEMATICS                                  |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE MQC_Variable_Addition
      function MQC_Variable_Addition(mqcVariable1,mqcVariable2) result(mqcVariableOut)
!
!     This function carries out addition of two MQC_Variables. The output is an
!     MQC_Variable type. This function provides the same functionality as the
!     fortran standard provides for intrinsic scalar and whole-array addition.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode1,typeCode2
      integer(kind=int64),dimension(:),allocatable::sumVectorInteger
      real(kind=real64),dimension(:),allocatable::sumVectorReal
!
!
!     Start by ensuring the two MQC Variables are conformable.
!
      if(.not.MQC_Variable_isConformable(mqcVariable1,mqcVariable2))  &
        call mqc_error_l('MQC Addition only allowed for conformable values.', 6, &
        'MQC_Variable_isConformable(mqcVariable1,mqcVariable2)', &
        MQC_Variable_isConformable(mqcVariable1,mqcVariable2) )
!
!     Get the right combination of integer/real and build a temporary result
!     array. This will allow uses of this function to work for C = A + B and for
!     C = C + B cases.
!
!     To keep things simple in the code, we use an encoding scheme here where
!     the current situation is encoded as a two-digit integer. The ten's digit
!     corresponds to the type code of <mqcVariable1> and the one's digit
!     corresponds to <mqcVariable2>. The individual digit values follow the
!     convension in Procedure MQC_Variable_getTypeCode.
!
!     Here is some prelim work...
!
      typeCode1 = MQC_Variable_getTypeCode(mqcVariable1)
      typeCode2 = MQC_Variable_getTypeCode(mqcVariable2)
      if(typeCode1.eq.0)  &
        call mqc_error_I('MQC_Variable_Addition: Var1 is of UNKONWN type.', 6, &
        'typeCode1', typeCode1 )
      if(typeCode2.eq.0)  &
        call mqc_error_I('MQC_Variable_Addition: Var2 is of UNKONWN type.', 6, &
        'typeCode2', typeCode2 )
      select case(MIN(typeCode1,typeCode2))
      case(2)
        allocate(sumVectorReal(SIZE(mqcVariable1)))
      case(3)
        allocate(sumVectorInteger(SIZE(mqcVariable1)))
      case default
        call mqc_error_I('MQC_Variable_Addition: Result is of UNKNOWN type.', 6, &
             'typeCode1', typeCode1, &
             'typeCode2', typeCode2 )
      end select
!
!     Now, do the addition and set the output function value accordingly.
!
      select case(typeCode1*10 + typeCode2)
      case(22)
        sumVectorReal = mqcVariable1%realArray + mqcVariable2%realArray
        mqcVariableOut = sumVectorReal
      case(23)
        sumVectorReal = mqcVariable1%realArray + mqcVariable2%integerArray
        mqcVariableOut = sumVectorReal
      case(32)
        sumVectorReal = mqcVariable1%integerArray(:) + mqcVariable2%realArray(:)
        mqcVariableOut = sumVectorReal
      case(33)
        sumVectorInteger = mqcVariable1%integerArray + mqcVariable2%integerArray
        mqcVariableOut = sumVectorInteger
      case default
        call mqc_error_I('MQC_Variable_Addition: Combination of var1 and var 2 types is UNKNOWN.', 6, &
             'typeCode1*10', typeCode1*10, &
             'typeCode2', typeCode2)
      end select
!
      return
      end function MQC_Variable_Addition


!
!PROCEDURE MQC_Variable_Subtraction
      function MQC_Variable_Subtraction(mqcVariable1,mqcVariable2) result(mqcVariableOut)
!
!     This function carries out subtraction of two MQC_Variables. The output is
!     an MQC_Variable type. This function provides the same functionality as the
!     fortran standard provides for intrinsic scalar and whole-array
!     subtraction.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode1,typeCode2
      integer(kind=int64),dimension(:),allocatable::differenceVectorInteger
      real(kind=real64),dimension(:),allocatable::differenceVectorReal
!
!
!     Start by ensuring the two MQC Variables are conformable.
!
      if(.not.MQC_Variable_isConformable(mqcVariable1,mqcVariable2))  &
        call mqc_error_l('MQC Subtraction only allowed for conformable values.', 6, &
        'MQC_Variable_isConformable(mqcVariable1,mqcVariable2)', &
        MQC_Variable_isConformable(mqcVariable1,mqcVariable2) )
!
!     Get the right combination of integer/real and build a temporary result
!     array. This will allow uses of this function to work for C = A - B and for
!     C = C - B cases.
!
!     To keep things simple in the code, we use an encoding scheme here where
!     the current situation is encoded as a two-digit integer. The ten's digit
!     corresponds to the type code of <mqcVariable1> and the one's digit
!     corresponds to <mqcVariable2>. The individual digit values follow the
!     convension in Procedure MQC_Variable_getTypeCode.
!
!     Here is some prelim work...
!
      typeCode1 = MQC_Variable_getTypeCode(mqcVariable1)
      typeCode2 = MQC_Variable_getTypeCode(mqcVariable2)
      if(typeCode1.eq.0)  &
        call mqc_error_I('MQC_Variable_Subtraction: Var1 is of UNKONWN type.', 6, &
        'typeCode1', typeCode1 )
      if(typeCode2.eq.0)  &
        call mqc_error_I('MQC_Variable_Subtraction: Var2 is of UNKONWN type.', 6, &
        'typeCode2', typeCode2 )
      select case(MIN(typeCode1,typeCode2))
      case(2)
        allocate(differenceVectorReal(SIZE(mqcVariable1)))
      case(3)
        allocate(differenceVectorInteger(SIZE(mqcVariable1)))
      case default
        call mqc_error_I('MQC_Variable_Subtraction: Result is of UNKNOWN type.', 6, &
             'typeCode1', typeCode1, &
             'typeCode2', typeCode2 )
      end select
!
!     Now, do the addition and set the output function value accordingly.
!
      select case(typeCode1*10 + typeCode2)
      case(22)
        differenceVectorReal = mqcVariable1%realArray - mqcVariable2%realArray
        mqcVariableOut = differenceVectorReal
      case(23)
        differenceVectorReal = mqcVariable1%realArray - mqcVariable2%integerArray
        mqcVariableOut = differenceVectorReal
      case(32)
        differenceVectorReal = mqcVariable1%integerArray(:) - mqcVariable2%realArray(:)
        mqcVariableOut = differenceVectorReal
      case(33)
        differenceVectorInteger = mqcVariable1%integerArray - mqcVariable2%integerArray
        mqcVariableOut = differenceVectorInteger
      case default
        call mqc_error_I('MQC_Variable_Subtraction: Combination of var1 and var 2 types is UNKNOWN.', 6, &
             'typeCode1*10', typeCode1*10, &
             'typeCode2', typeCode2 )
      end select
!
      return
      end function MQC_Variable_Subtraction


!
!PROCEDURE MQC_Variable_Multiplication
      function MQC_Variable_Multiplication(mqcVariable1,mqcVariable2) result(mqcVariableOut)
!
!     This function carries out multiplication of two MQC_Variables. The output
!     is an MQC_Variable type. This function provides the same functionality as
!     the fortran standard provides for intrinsic scalar and whole-array
!     multiplication.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode1,typeCode2
      integer(kind=int64),dimension(:),allocatable::productVectorInteger
      real(kind=real64),dimension(:),allocatable::productVectorReal
!
!
!     Start by ensuring the two MQC Variables are conformable.
!
      if(.not.MQC_Variable_isConformable(mqcVariable1,mqcVariable2))  &
        call mqc_error_l('MQC Multiplication only allowed for conformable values.', 6, &
        'MQC_Variable_isConformable(mqcVariable1,mqcVariable2)', &
        MQC_Variable_isConformable(mqcVariable1,mqcVariable2) )
!
!     Get the right combination of integer/real and build a temporary result
!     array. This will allow uses of this function to work for C = A * B and for
!     C = C * B cases.
!
!     To keep things simple in the code, we use an encoding scheme here where
!     the current situation is encoded as a two-digit integer. The ten's digit
!     corresponds to the type code of <mqcVariable1> and the one's digit
!     corresponds to <mqcVariable2>. The individual digit values follow the
!     convension in Procedure MQC_Variable_getTypeCode.
!
!     Here is some prelim work...
!
      typeCode1 = MQC_Variable_getTypeCode(mqcVariable1)
      typeCode2 = MQC_Variable_getTypeCode(mqcVariable2)
      if(typeCode1.eq.0)  &
        call mqc_error_I('MQC_Variable_Multiplication: Var1 is of UNKONWN type.', 6, &
        'typeCode1', typeCode1)
      if(typeCode2.eq.0)  &
        call mqc_error_I('MQC_Variable_Multiplication: Var2 is of UNKONWN type.', 6, &
        'typeCode2' ,typeCode2);
      select case(MIN(typeCode1,typeCode2))
      case(2)
        allocate(productVectorReal(SIZE(mqcVariable1)))
      case(3)
        allocate(productVectorInteger(SIZE(mqcVariable1)))
      case default
        call mqc_error_i('MQC_Variable_Multiplication: Result is of UNKNOWN type.', 6, &
             'typeCode1', typeCode1, &
             'typeCode2', typeCode2 )
      end select
!
!     Now, do the addition and set the output function value accordingly.
!
      select case(typeCode1*10 + typeCode2)
      case(22)
        productVectorReal = mqcVariable1%realArray * mqcVariable2%realArray
        mqcVariableOut = productVectorReal
      case(23)
        productVectorReal = mqcVariable1%realArray * mqcVariable2%integerArray
        mqcVariableOut = productVectorReal
      case(32)
        productVectorReal = mqcVariable1%integerArray(:) * mqcVariable2%realArray(:)
        mqcVariableOut = productVectorReal
      case(33)
        productVectorInteger = mqcVariable1%integerArray * mqcVariable2%integerArray
        mqcVariableOut = productVectorInteger
      case default
        call mqc_error_I('MQC_Variable_Multiplication: Combination of var1 and var 2 types is UNKNOWN.', 6, &
             'typeCode1*10', typeCode1*10, &
             'typeCode2', typeCode2 )
      end select
!
      return
      end function MQC_Variable_Multiplication


!
!PROCEDURE MQC_Variable_Division
      function MQC_Variable_Division(mqcVariable1,mqcVariable2) result(mqcVariableOut)
!
!     This function carries out division of two MQC_Variables. The output is an
!     MQC_Variable type. This function provides the same functionality as the
!     fortran standard provides for intrinsic scalar and whole-array division.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode1,typeCode2
      integer(kind=int64),dimension(:),allocatable::quotientVectorInteger
      real(kind=real64),dimension(:),allocatable::quotientVectorReal
!
!
!     Start by ensuring the two MQC Variables are conformable.
!
      if(.not.MQC_Variable_isConformable(mqcVariable1,mqcVariable2))  &
        call mqc_error_l('MQC Division only allowed for conformable values.', 6, &
        'MQC_Variable_isConformable(mqcVariable1,mqcVariable2)', &
        MQC_Variable_isConformable(mqcVariable1,mqcVariable2) )
!
!     Get the right combination of integer/real and build a temporary result
!     array. This will allow uses of this function to work for C = A / B and for
!     C = C / B cases.
!
!     To keep things simple in the code, we use an encoding scheme here where
!     the current situation is encoded as a two-digit integer. The ten's digit
!     corresponds to the type code of <mqcVariable1> and the one's digit
!     corresponds to <mqcVariable2>. The individual digit values follow the
!     convension in Procedure MQC_Variable_getTypeCode.
!
!     Here is some prelim work...
!
      typeCode1 = MQC_Variable_getTypeCode(mqcVariable1)
      typeCode2 = MQC_Variable_getTypeCode(mqcVariable2)
      if(typeCode1.eq.0)  &
        call mqc_error_I('MQC_Variable_Division: Var1 is of UNKONWN type.', 6, &
        'typeCode1', typeCode1 )
      if(typeCode2.eq.0)  &
        call mqc_error_I('MQC_Variable_Division: Var2 is of UNKONWN type.', 6, &
        'typeCode2', typeCode2 )
      select case(MIN(typeCode1,typeCode2))
      case(2)
        allocate(quotientVectorReal(SIZE(mqcVariable1)))
      case(3)
        allocate(quotientVectorInteger(SIZE(mqcVariable1)))
      case default
        call mqc_error_I('MQC_Variable_Division: Result is of UNKNOWN type.', 6, &
             'typeCode1', typeCode1, &
             'typeCode2', typeCode2 )
      end select
!
!     Now, do the addition and set the output function value accordingly.
!
      select case(typeCode1*10 + typeCode2)
      case(22)
        quotientVectorReal = mqcVariable1%realArray / mqcVariable2%realArray
        mqcVariableOut = quotientVectorReal
      case(23)
        quotientVectorReal = mqcVariable1%realArray / mqcVariable2%integerArray
        mqcVariableOut = quotientVectorReal
      case(32)
        quotientVectorReal = mqcVariable1%integerArray(:) / mqcVariable2%realArray(:)
        mqcVariableOut = quotientVectorReal
      case(33)
        quotientVectorInteger = mqcVariable1%integerArray / mqcVariable2%integerArray
        mqcVariableOut = quotientVectorInteger
      case default
        call mqc_error_I('MQC_Variable_Division: Combination of var1 and var 2 types is UNKNOWN.', 6, &
             'typeCode1*10', typeCode1*10, &
             'typeCode2', typeCode2)
      end select
!
      return
      end function MQC_Variable_Division


!
!PROCEDURE MQC_Variable_Contraction_Full
      function MQC_Variable_Contraction_Full(mqcVariable1,mqcVariable2) result(mqcVariableOut)
!
!     This function carries out full contraction of two MQC variables. The two
!     input MQC variables must be conformable.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode1,typeCode2
      integer(kind=int64)::contractionInteger
      real(kind=real64)::contractionReal
!
!
!     Start by ensuring the two MQC Variables are conformable.
!
      if(.not.MQC_Variable_isConformable(mqcVariable1,mqcVariable2))  &
        call mqc_error_l('MQC Division only allowed for conformable values.', 6, &
        'MQC_Variable_isConformable(mqcVariable1,mqcVariable2)', &
        MQC_Variable_isConformable(mqcVariable1,mqcVariable2) )

!
!     Get the right combination of integer/real and form the result scalar.
!
!     To keep things simple in the code, we use an encoding scheme here where
!     the current situation is encoded as a two-digit integer. The ten's digit
!     corresponds to the type code of <mqcVariable1> and the one's digit
!     corresponds to <mqcVariable2>. The individual digit values follow the
!     convension in Procedure MQC_Variable_getTypeCode.
!
!     Here is some prelim work...
!
      typeCode1 = MQC_Variable_getTypeCode(mqcVariable1)
      typeCode2 = MQC_Variable_getTypeCode(mqcVariable2)
      if(typeCode1.eq.0)  &
        call mqc_error_I('MQC_Variable_Division: Var1 is of UNKONWN type.', 6, &
        'typeCode1', typeCode1 )

      if(typeCode2.eq.0)  &
        call mqc_error_I('MQC_Variable_Division: Var2 is of UNKONWN type.', 6, &
'typeCode2', typeCode2 )

!
!     Now, do the contraction and set the output function value accordingly.
!
      select case(typeCode1*10 + typeCode2)
      case(22)
        contractionReal = dot_product(mqcVariable1%realArray,mqcVariable2%realArray)
        mqcVariableOut = contractionReal
      case(23)
        contractionReal = dot_product(mqcVariable1%realArray,mqcVariable2%integerArray)
        mqcVariableOut = contractionReal
      case(32)
        contractionReal = dot_product(mqcVariable1%integerArray,mqcVariable2%realArray)
        mqcVariableOut = contractionReal
      case(33)
        contractionInteger = dot_product(mqcVariable1%integerArray,mqcVariable2%integerArray)
        mqcVariableOut = contractionInteger
      case default
        call mqc_error_I('MQC_Variable_Contraction_Full: Combination of var1 and var 2 types is UNKNOWN.', 6, &
             'typeCode1', typeCode1, &
             'typeCode2', typeCode2 )
      end select
!
      return
      end function MQC_Variable_Contraction_Full


!
!PROCEDURE MQC_Variable_dot_product
      function MQC_Variable_dot_product(mqcVariable1,mqcVariable2) result(mqcVariableOut)
!
!     This function solves the dot product of two vectors stored in the input
!     dummy MQC_Variables arguments <mqcVariable1> and <mqcVariable2>. If these
!     two arguments do not correspond to vectors, this function fails.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      type(MQC_Variable)::mqcVariableOut
!
!
!     Start by ensuring the two MQC Variables are vectors and are conformable.
!
      if(mqcVariable1%getRank().ne.1.or.mqcVariable2%getRank().ne.1)  &
        call mqc_error_i('dot_product must be between two vectors.',6, &
        'RANK(mqcVariable1)',mqcVariable1%getRank(),  &
        'RANK(mqcVariable2)',mqcVariable2%getRank())
      if(.not.MQC_Variable_isConformable(mqcVariable1,mqcVariable2))  &
        call mqc_error_l('MQC Division only allowed for conformable values.',6, &
        'MQC_Variable_isConformable(mqcVariable1,mqcVariable2)', &
        MQC_Variable_isConformable(mqcVariable1,mqcVariable2))
!
!     Calculate the dot product.
!
      mqcVariableOut = MQC_Variable_Contraction_Full(mqcVariable1,mqcVariable2)
!
      return
      end function MQC_Variable_dot_product


!
!
      end module MQC_Algebra2
