      Module MQC_Algebra2
!MD+
!#MQCPack: MQC_Algebra2
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
!     **                          Version 25.3.0                           **
!     **                          March 25, 2025                           **
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
!     This module is a new (relative to mqc_algebra) module that provides
!     advanced numerical variable functionalities and easy-to-use algebra tools.
!
!     The functions and subroutines contained within this module are grouped
!     into 
!           1. MQC_Variable Management Procedures;
!           2. MQC_Variable Printing Procedures; and
!           3. MQC_Variable Value Setting & Basic Manipulation Procedures; and
!           4. MQC_Variable Mathematics Procedures;
!MD-

!
      Use MQC_General
      use iso_fortran_env, only: int64, real64
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     Variable...
!
!     The MQC_Variable class is a flexible numerical structure. MQC_Variable
!     objects may represent scalars, vectors, matrices, or general tensors.
!     Additionally, MQC_Variable objects may represent such structures as
!     integer, real, or complex types. The goal is for such objects to be used
!     in general form, abstracting away storage format and the like.
!
!     Elements of the MQC_Variable class include:
!       - rank          The rank of the numerical array structure. This value
!                       defaults to -1 to indicate that the object has not yet
!                       been initialized.
!       - dataType      This character string encodes the data type of the
!                       object as "INTEGER", "REAL", or "COMPLEX".
      type MQC_Variable
        integer(kind=int64),private::rank=-1
        character(len=64),private::dataType,storageFormat
        integer(kind=int64),dimension(10),private::dimensions=0
        real(kind=real64),dimension(:),allocatable::realArray
        complex(kind=real64),dimension(:),allocatable::complexArray
        integer(kind=int64),dimension(:),allocatable,private::integerArray
        logical,private::initialized=.false.
      Contains
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
        procedure,public ::putMQC       => MQC_Variable_put_MQC
        generic,public   ::put          => MQC_Variable_put_MQC,  &
                                           MQC_Variable_put_intrinsicInteger,  &
                                           MQC_Variable_put_intrinsicReal
        procedure,public ::getVal        => MQC_Variable_get_MQC
        procedure,public ::getRank       => MQC_Variable_getRank
        procedure,public ::getType       => MQC_Variable_getType
        procedure,public ::print         => MQC_Print_mqcVariable
        procedure,public ::isConformable => MQC_Variable_isConformable
        procedure,public ::trimZero      => mqc_variable_trimZero
        procedure,public ::trace         => mqc_variable_matrixTrace
        procedure,public ::iPower        => MQC_Variable_integerPower
        procedure,public ::rPower        => MQC_Variable_realPower
        procedure,public ::column        => MQC_Variable_MatrixGetColumn
        procedure,public ::eigen         => MQC_Variable_Matrix_Eigensystem
        procedure,public ::diag          => MQC_Variable_DiagonalMatrix2FullMatrix
        procedure,public ::det           => MQC_Variable_Matrix_Determinant
        procedure,public ::subMatrix     => MQC_Variable_SubMatrix
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
!     Interface MQC Variables to the RESHAPE intrinsic function.
      interface reshape
        module procedure MQC_Variable_reshape
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
!     Interface MQC Varaiable full contraction.
      interface contraction
        module procedure MQC_Variable_Contraction_Full
      end interface
!
!     Interface MQC Varaiable full contraction for vectors --> dot_product.
      interface dot_product
        module procedure MQC_Variable_dot_product
      end interface
!
!     Interface MQC Varaiable vector-vector outer product --> outer_product.
      interface outer_product
        module procedure MQC_Variable_outer_product
      end interface
!
!     Interface MQC Varaiable full SUM.
      interface sum
        module procedure MQC_Variable_Array_Sum
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
        module procedure MQC_Variable_setVal_intrinsicRank42mqcRank4_real
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
!     Interface for transposing a matrix.
      interface transpose
        module procedure MQC_Variable_Matrix_Transpose
      end interface
!
!     Interface for MatMul.
      interface MatMul
        module procedure MQC_Variable_Matrix_Multiplication
      end interface
!
!     Interface for MQC_Variable_Power.
      interface MQC_Variable_Power
        module procedure MQC_Variable_IntegerPower
        module procedure MQC_Variable_RealPower
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
        module procedure MQC_Variable_setVal_intrinsicRank42mqcRank4_real
        module procedure MQC_Variable_mqc2intrinsicIntegerScalar
        module procedure MQC_Variable_mqc2intrinsicRealScalar
        module procedure MQC_Variable_mqc2intrinsicInteger1Array
        module procedure MQC_Variable_mqc2intrinsicReal1Array
        module procedure MQC_Variable_mqc2intrinsicInteger2Array
        module procedure MQC_Variable_mqc2intrinsicReal2Array
        module procedure MQC_Variable_mqc2intrinsicReal4Array
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
!----------------------------------------------------------------
!                                                               |
!     MODULE VARIABLES/CONSTANTS                                |
!                                                               |
!----------------------------------------------------------------
!
      Integer,Private::MQC_Algebra2_IOut=6
      Logical,Private::MQC_Algebra2_DEBUG=.false.
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
!     MODULE CONTROL PROCEDURES                                 |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE MQC_Algebra2_SetIOut
      subroutine MQC_Algebra2_SetIOut(setIOutValue)
!
!     This subroutine is called to set a Module-Wide default output unit number.
!     The default is unit 6.
!
!     H. P. Hratchian, 2020.
!
!
      implicit none
      integer(kind=int64),intent(IN)::setIOutValue
!
      MQC_Algebra2_IOut = setIOutValue
!
      return
      end subroutine MQC_Algebra2_SetIOut


!
!PROCEDURE MQC_Algebra2_SetDEBUG
      subroutine MQC_Algebra2_SetDEBUG(setDebugValue)
!
!     This subroutine is called to set a Module-Wide DEBUG flag to .TRUE. or
!     .FALSE.
!
!     H. P. Hratchian, 2020.
!
!
      implicit none
      logical,intent(IN)::setDebugValue
!
      MQC_Algebra2_DEBUG = setDebugValue
!
      return
      end subroutine MQC_Algebra2_SetDEBUG


!
!----------------------------------------------------------------
!                                                               |
!     MQC_Variable MANAGEMENT PROCEDURES                        |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE MQC_Variable_initialize
      subroutine MQC_Variable_initialize(mqcVariable,dataType,dimensions,  &
        rank,storageFormat,debug)
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
      logical,intent(in),optional::debug
!
      character(len=64)::myDataType
      integer(kind=int64)::myRank,myLength
      character(len=64)::myStorageFormat
      logical::myDebug
!
!
!     Set-up <myRank> and ensure the sent <rank> and <dimensions> are
!     compatable.
!
      if(Present(debug)) then
        myDebug = debug
      else
        myDebug = .False.
      endIf
      myDebug = myDebug.or.MQC_Algebra2_DEBUG
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
          call MQC_Error_I('Logic error in MQC_Variable_initialize:'// &
            ' Dimensions and Rank do not jive!', 6, &
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
        call MQC_Error_A('Initializing MQC_Variable with unknown'//  &
          ' data type.', 6,'myDataType', myDataType)
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
            call MQC_Error_L('Illegal attempt to initialize '//  &
              'MQC_Variable of rank > 0 without defined dimensions.',  &
              6,'Present(dimensions)', Present(dimensions) )
          mqcVariable%dimensions = 0
          mqcVariable%dimensions(1:myRank) = dimensions
        case default
          call MQC_Error_I('Illegal rank in MQC_Variable_initialize.',  &
            6,'myRank',myRank)
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

      if(myDebug) then
        write(*,*)
        write(*,*)' in Init: myLength = ',myLength
        write(*,*)
      endIf

      select case(mqcVariable%dataType)
      case('REAL')
        if(Allocated(mqcVariable%integerArray))  &
          DeAllocate(mqcVariable%integerArray)
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
        call MQC_Error_A('MQC_Variable_initialize is confused in '//  &
          'allocation block.',6,'mqcVariable%dataType',  &
          mqcVariable%dataType)
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
      function MQC_Variable_getArrayPosition(mqcVariable,arrayElement,DEBUG) result(k)
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
      logical,optional::DEBUG
      logical::DEBUG1
!
      DEBUG1 = .false.
      if(Present(DEBUG)) DEBUG1 = DEBUG
!
!
!     Do the work...
!
      select case(TRIM(mqcVariable%storageFormat))
      case('FULL')
        k = 0
        subLength = 1
        do i=1,mqcVariable%rank
          subLength = subLength*mqcVariable%dimensions(i)
        endDo
        do i=mqcVariable%rank,2,-1
          subLength = subLength/mqcVariable%dimensions(i)
          k = k + (arrayElement(i)-1)*subLength
          if(DEBUG1) write(*,*) i,subLength,k
        endDo
        k = k + arrayElement(1)
      case default
        call mqc_error_a('MQC_Variable_getArrayPosition: '//  &
          'Found an unknown storage format.',6,  &
          'TRIM(mqcVariable%storageFormat)',  &
          TRIM(mqcVariable%storageFormat))
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
        call MQC_Error_I('MQC_Variable_getSize: '//  &
        'Illegal iDimension sent.', 6,&
        'myDimension', myDimension, &
        'MQC_Variable_getRank(mqcVariable)',  &
        MQC_Variable_getRank(mqcVariable))
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
!PROCEDURE MQC_Variable_reshape
      function MQC_Variable_reshape(mqcVariable,newShape) result(mqcVariableNew)
!
!     This function reshapes an MQC variable array.
!
!
!     H. P. Hratchian, 2025.
!
!
!     Variable Declarations.
      implicit none
      type(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64),dimension(:),intent(in)::newShape
      type(MQC_Variable)::mqcVariableNew
!
!
      if(PRODUCT(newShape).ne.SIZE(mqcVariable))  &
        call mqc_error('RESHAPE: Inconsistend new shape and size.')
      mqcVariableNew = mqcVariable
      mqcVariableNew%dimensions=0
      mqcVariableNew%dimensions(1:SIZE(newShape))=newShape
      mqcVariableNew%rank=SIZE(newShape)
!
      return
      end function MQC_Variable_reshape


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
!PROCEDURE MQC_Variable_getShape
      subroutine MQC_Variable_getShape(mqcVariable,iShape)
!
!     This subroutine returns the allocatable (intrinsic) integer array
!     <iShape> with the shape of <mqcVariable>. Dummy variable <iShape>
!     must be allocated before calling this routine.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64),dimension(:)::iShape
      integer(kind=int64)::myRank,iShapeSize
!
!
!     Do the work...
!
      myRank = RANK(mqcVariable)
      iShapeSize = SIZE(iShape)
      if(iShapeSize.lt.myRank)  &
        call mqc_error('getShape: iShape is too small.')
      if(iShapeSize.gt.myRank) iShape(:iShapeSize) = 0
      select case(myRank)
      case(:0)
        iShape(1) = 0
      case(1:)
        iShape(1:myRank) = mqcVariable%dimensions(1:myRank)
      end select
!
      return
      end subroutine MQC_Variable_getShape


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
!     H. P. Hratchian, 2017, 2020.
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
            call mqc_print_scalar(mqcVariable%realArray(1),myIOut,Header,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_scalar(mqcVariable%realArray(1),myIOut,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case('INTEGER')
          if(PRESENT(header)) then
            call mqc_print_scalar(mqcVariable%integerArray(1),myIOut,Header,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_scalar(mqcVariable%integerArray(1),myIOut,  &
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
            call mqc_print_vector(mqcVariable%realArray(:),myIOut,Header,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_vector(mqcVariable%realArray(:),myIOut,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case('INTEGER')
          if(PRESENT(header)) then
            call mqc_print_vector(mqcVariable%integerArray(:),myIOut,Header,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_vector(mqcVariable%integerArray(:),myIOut,  &
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
            call mqc_print_matrix(  &
              RESHAPE(mqcVariable%realArray(:),mqcVariable%dimensions(1:2)),myIOut,  &
              Header,Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_matrix(  &
              RESHAPE(mqcVariable%realArray,mqcVariable%dimensions(1:2)),myIOut,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case('INTEGER')
          if(PRESENT(header)) then
            call mqc_print_matrix(  &
              RESHAPE(mqcVariable%integerArray,mqcVariable%dimensions(1:2)),myIOut,  &
              Header,Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_matrix(  &
              RESHAPE(mqcVariable%integerArray,mqcVariable%dimensions(1:2)),myIOut,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case default
          call mqc_error_a('MQC variable print requested for unknown type.', 6, &
               'mqcVariable%getType()', mqcVariable%getType())
        end select
!
!     Print rank-4 tensors...
      case(4)
        select case(mqcVariable%getType())
        case('REAL')
          if(PRESENT(header)) then
            call mqc_print_rank4Tensor_array_real(  &
              RESHAPE(mqcVariable%realArray(:),mqcVariable%dimensions(1:4)),myIOut,  &
              Header,Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          else
            call mqc_print_rank4Tensor_array_real(  &
              RESHAPE(mqcVariable%realArray,mqcVariable%dimensions(1:4)),myIOut,  &
              Blank_At_Top=myBlankAtTop,Blank_At_Bottom=myBlankAtBottom)
          endIf
        case default
          call mqc_error_a('MQC variable print requested for unknown type.', 6, &
               'mqcVariable%getType()', mqcVariable%getType())
        end select
!
!     Print requested for a general MQC variable with a rank not yet coded...
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
!     MQC_Variable VALUE SETTING, RETREVAL, & BASIC             |
!     MANIPULATION PROCEDURES                                   |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE MQC_Variable_get_MQC
      function MQC_Variable_get_MQC(mqcVariable,arrayElement) result(valueMQC)
!
!     This subroutine is used to get a value from a specified element of an
!     MQC_Variable array. The value is returned as an MQC variable type. The
!     array element is sent in dummy argument <arrayElement>.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      type(MQC_Variable)::valueMQC
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),dimension(:)::arrayElement
      integer(kind=int64)::k
!
!
!     Do the work...
!
      if(.not.mqcVariable%initialized) call mqc_error_l('Cannot use MQC_Variable_Get with uninitialized MQC variable.', 6, &
        'mqcVariable%initialized', mqcVariable%initialized )
      if(mqcVariable%rank.ne.SIZE(arrayElement))  &
        call mqc_error_i('MQC_Variable_Get: Rank of variable does NOT match number of array indices provided.', 6, &
        'mqcVariable%rank', mqcVariable%rank, &
        'SIZE(arrayElement)', SIZE(arrayElement))
      k = MQC_Variable_getArrayPosition(mqcVariable,arrayElement)
      if(k.gt.MQC_Variable_getLength(mqcVariable)) call mqc_error_i('MQC_Variable_Get: Invalid array index provided.', 6, &
        'k', k, &
        'MQC_Variable_getLength(mqcVariable)', MQC_Variable_getLength(mqcVariable))
      select case(TRIM(MQC_Variable_getType(mqcVariable)))
      case('REAL')
        valueMQC = mqcVariable%realArray(k)
      case('INTEGER')
        valueMQC = mqcVariable%integerArray(k)
      case default
        call mqc_error_a('MQC_Variable_Put: Unknown MQC variable type.', 6, &
          'TRIM(MQC_Variable_getType(mqcVariable))', TRIM(MQC_Variable_getType(mqcVariable)) )
      end select
!
      return
      end function MQC_Variable_get_MQC


!
!PROCEDURE MQC_Variable_MatrixGetColumn
      function MQC_Variable_MatrixGetColumn(mqcVariable,iColumn)  &
        result(columnVector)
!
!     This subroutine is used to cut out a single column from a matrix. The
!     output is a real vector.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),intent(in)::iColumn
      real(kind=real64),dimension(:),allocatable::columnVector
!
      integer(kind=int64)::nRows,nColumns
      real(kind=real64),dimension(:,:),allocatable::tmpMatrix
!
!
!     Ensure that <mqcVariable> is a matrix. Then allocate <columnVector> and
!     fill it with the requested column.
!
      if(.not.mqcVariable%initialized)  &
        call mqc_error('Cannot use MQC_Variable_MatrixGetColumn with uninitialized permuteColumns.')
      if(mqcVariable%rank.ne.2)  &
        call mqc_error('MQC_Variable_MatrixGetColumn: Rank of variable must be 2.')
      if(TRIM(mqcVariable%storageFormat).ne.'FULL')   &
        call mqc_error_a('MQC_Variable_MatrixGetColumn: Can only be used on FULL storage arrays.',  &
        6,'Current storage form:',TRIM(mqcVariable%storageFormat))
      nRows = SIZE(mqcVariable,1)
      nColumns = SIZE(mqcVariable,2)
      Allocate(columnVector(nRows))
      Allocate(tmpMatrix(nRows,nColumns))
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        tmpMatrix = RESHAPE(mqcVariable%realArray(:),mqcVariable%dimensions(1:2))
        columnVector = tmpMatrix(:,iColumn)
      case(3)
        tmpMatrix = RESHAPE(mqcVariable%integerArray(:),mqcVariable%dimensions(1:2))
        columnVector = tmpMatrix(:,iColumn)
      case default
        call mqc_error('GetColumn: Erroneous typing found.')
      end select
      Deallocate(tmpMatrix)
!
      return
      end function MQC_Variable_MatrixGetColumn


!
!PROCEDURE MQC_Variable_MatrixPermuteColumns
      subroutine MQC_Variable_MatrixPermuteColumns(mqcVariable,columnI,columnJ)
!
!     This subroutine is used to permute 2 columns in a matrix. If this routine
!     is called with an mqcVariable that isn't rank 2, it results in error.
!
!
!     H. P. Hratchian, 2021.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      integer(kind=int64),intent(in)::columnI,columnJ
      integer(kind=int64)::lenVector,iStart,jStart,iEnd,jEnd
      integer(kind=int64),dimension(:),allocatable::tmpArrayI
      real(kind=real64),dimension(:),allocatable::tmpArrayR
      integer(kind=int64)::k
!
!
!     Do the work...
!
      if(.not.mqcVariable%initialized) call mqc_error('Cannot use MQC_Variable_Get with uninitialized permuteColumns.')
      if(mqcVariable%rank.ne.2)  &
        call mqc_error('permuteColumns: Rank of variable must be 2.')
      if(TRIM(mqcVariable%storageFormat).ne.'FULL')   &
        call mqc_error_a('permuteColumns: Can only be used on FULL storage arrays.',  &
        6,'Current storage form:',TRIM(mqcVariable%storageFormat))
      lenVector = SIZE(mqcVariable,1)
      iStart = MQC_Variable_getArrayPosition(mqcVariable,[ 1, columnI ])
      jStart = MQC_Variable_getArrayPosition(mqcVariable,[ 1, columnJ ])
      iEnd = iStart + lenVector - 1
      jEnd = jStart + lenVector - 1
      select case(TRIM(MQC_Variable_getType(mqcVariable)))
      case('REAL')
        Allocate(tmpArrayR(lenVector))
        tmpArrayR = mqcVariable%realArray(iStart:iEnd)
        mqcVariable%realArray(iStart:iEnd) = mqcVariable%realArray(jStart:jEnd)
        mqcVariable%realArray(jStart:jEnd) = tmpArrayR
        DeAllocate(tmpArrayR)
      case('INTEGER')
        Allocate(tmpArrayI(lenVector))
        tmpArrayI = mqcVariable%integerArray(iStart:iEnd)
        mqcVariable%integerArray(iStart:iEnd) = mqcVariable%integerArray(jStart:jEnd)
        mqcVariable%integerArray(jStart:jEnd) = tmpArrayI
        DeAllocate(tmpArrayI)
      case default
        call mqc_error_a('permuteColumns: Unknown MQC variable type.', 6, &
          'TRIM(MQC_Variable_getType(mqcVariable))', TRIM(MQC_Variable_getType(mqcVariable)) )
      end select
!
      return
      end subroutine MQC_Variable_MatrixPermuteColumns


!hph+
!!
!!PROCEDURE MQC_Variable_MatrixOrderedColumns
!      subroutine MQC_Variable_MatrixOrderedColumns(mqcVariable,map)
!!
!!     This subroutine is used to order the columns of a matrix according to the
!!     array map. If this routine is called with an mqcVariable that isn't rank
!!     2, it results in error.
!!
!!
!!     H. P. Hratchian, 2023.
!!
!!
!!     Variable Declarations.
!      implicit none
!      class(MQC_Variable)::mqcVariable
!      integer(kind=int64),intent(in)::columnI,columnJ
!      integer(kind=int64)::lenVector,iStart,jStart,iEnd,jEnd
!      integer(kind=int64),dimension(:),allocatable::tmpArrayI
!      real(kind=real64),dimension(:),allocatable::tmpArrayR
!      integer(kind=int64)::k
!!
!!
!!     Do the work...
!!
!      if(.not.mqcVariable%initialized) call mqc_error('Cannot use MQC_Variable_Get with uninitialized permuteColumns.')
!      if(mqcVariable%rank.ne.2)  &
!        call mqc_error('permuteColumns: Rank of variable must be 2.')
!      if(TRIM(mqcVariable%storageFormat).ne.'FULL')   &
!        call mqc_error_a('permuteColumns: Can only be used on FULL storage arrays.',  &
!        6,'Current storage form:',TRIM(mqcVariable%storageFormat))
!      lenVector = SIZE(mqcVariable,1)
!      iStart = MQC_Variable_getArrayPosition(mqcVariable,[ 1, columnI ])
!      jStart = MQC_Variable_getArrayPosition(mqcVariable,[ 1, columnJ ])
!      iEnd = iStart + lenVector - 1
!      jEnd = jStart + lenVector - 1
!      select case(TRIM(MQC_Variable_getType(mqcVariable)))
!      case('REAL')
!        Allocate(tmpArrayR(lenVector))
!        tmpArrayR = mqcVariable%realArray(iStart:iEnd)
!        mqcVariable%realArray(iStart:iEnd) = mqcVariable%realArray(jStart:jEnd)
!        mqcVariable%realArray(jStart:jEnd) = tmpArrayR
!        DeAllocate(tmpArrayR)
!      case('INTEGER')
!        Allocate(tmpArrayI(lenVector))
!        tmpArrayI = mqcVariable%integerArray(iStart:iEnd)
!        mqcVariable%integerArray(iStart:iEnd) = mqcVariable%integerArray(jStart:jEnd)
!        mqcVariable%integerArray(jStart:jEnd) = tmpArrayI
!        DeAllocate(tmpArrayI)
!      case default
!        call mqc_error_a('permuteColumns: Unknown MQC variable type.', 6, &
!          'TRIM(MQC_Variable_getType(mqcVariable))', TRIM(MQC_Variable_getType(mqcVariable)) )
!      end select
!!
!      return
!      end subroutine MQC_Variable_MatrixPermuteColumns
!hph-


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
        mqcVariable%realArray(k) = Float(valueMQC)
      case('INTEGER')
        mqcVariable%integerArray(k) = Int(valueMQC)
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
        call mqc_error_A(  &
          'MQC_Variable_clear: Unknown type of MQC variable sent.', 6, &
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
!PROCEDURE MQC_Variable_mqc2mqc_subroutine
      subroutine MQC_Variable_mqc2mqc_subroutine(mqcVariable,mqcIn)
!
!     This subroutine is used to set an MQC Variable equal to another MQC
!     Variable. Note that the name of this procedure ends with
!     "_subroutine" to distinguish it from the closely related function
!     that calls it.
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
        call mqc_error_A('MQC_Variable_mqc2mqc: Confused by data type.',  &
           6,'TRIM(mqcIn%getType())', TRIM(mqcIn%getType()) )
      end select
      mqcVariable%initialized = .true.
!
      return
      end subroutine MQC_Variable_mqc2mqc_subroutine


!
!PROCEDURE MQC_Variable_mqc2mqc
      function MQC_Variable_mqc2mqc(mqcVariableIn) result(mqcVariableOut)
!
!     This function sets a new MQC_Variable equal to the one sent here.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariableIn
      type(MQC_Variable)::mqcVariableOut
!
!
!     Call Routine MQC_Variable_mqc2mqc_subroutine...
      call MQC_Variable_mqc2mqc_subroutine(mqcVariableIn,mqcVariableOut)
!
      return
      end function MQC_Variable_mqc2mqc


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
      integer(kind=int64)::i,nInputs,myRank
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
      DEBUG = DEBUG.or.MQC_Algebra2_DEBUG
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
        flush(MQC_Algebra2_IOut)
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
        flush(MQC_Algebra2_IOut)
      endIf
!
!     Set myDataType and myRank.
!
      if(fillInteger) then
        myDataType = 'INTEGER'
      elseIf(fillReal) then
        myDataType = 'REAL'
      else
        call MQC_Error_L(  &
          'MQC_Variable_setVal: Confused setting myDataType.', 6, &
          'fillInteger', fillInteger, &
          'fillReal', fillReal )
      endIf
      if(DEBUG) write(*,*)'   myDataType = ',TRIM(myDataType)
      if(PRESENT(dimensions)) then
        myRank = SIZE(dimensions)
        do i = SIZE(dimensions),1,-1
          if(dimensions(i).ne.0) exit
          myRank = myRank - 1
        endDo
        if(myRank.eq.1.and.fillScalar) myRank = 0
        if(myRank.gt.0) then
          Allocate(myDimensions(myRank))
          myDimensions = dimensions(1:myRank)
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
            call mqc_error_l(  &
              'MQC_Variable_setVal: Confused setting myDimensions.', 6, &
              'fillScalar', fillScalar, &
              'fillArray', fillArray, &
              'fillInteger', fillInteger, &
              'fillReal', fillReal )
          endIf
        else
          call MQC_Error_L(  &
            'MQC_Variable_setVal: Confused setting myRank.', 6, &
            'fillScalar', fillScalar, &
            'fillArray', fillArray )
        endIf
      endIf
!
!     Initialize mqcVariable and then set the value appropriately.
!
      if(DEBUG) write(*,*)'   Calling Init...'
      if(myRank.eq.0) then
        call MQC_Variable_initialize(mqcVariable,myDataType,debug=DEBUG)
      else
        call MQC_Variable_initialize(mqcVariable,myDataType,myDimensions,  &
          myRank,'FULL',debug=DEBUG)
      endIf
      if(DEBUG) then
        write(*,*)'   Back from init.'
        write(*,*)'   rank = ',mqcVariable%rank
        write(*,*)'   dimensions = ',mqcVariable%dimensions
        write(*,*)'   ALLOCATED(realArray) ?    = ',ALLOCATED(mqcVariable%realArray)
        write(*,*)'   ALLOCATED(integerArray) ? = ',ALLOCATED(mqcVariable%integerArray)
        if(Allocated(mqcVariable%realArray)) write(*,*)'   SIZE(realArray) = ',SIZE(mqcVariable%realArray)
        if(Allocated(mqcVariable%integerArray)) write(*,*)'   SIZE(integerArray) = ',SIZE(mqcVariable%integerArray)
        flush(MQC_Algebra2_IOut)
      endIf

      if(DEBUG) then
        write(*,*)
        write(*,*)' fillScalar = ',fillScalar
        write(*,*)' fillReal   = ',fillReal
        write(*,*)
      endIf

      if(fillScalar) then
        if(fillInteger) then
          mqcVariable%integerArray(1) = scalarIntegerIn
        elseIf(fillReal) then
          mqcVariable%realArray(1) = scalarRealIn
        else
          call MQC_Error_L(  &
            'MQC_Variable_setVal: Confused filling scalar value.', 6, &
            'fillScalar', fillScalar, &
            'fillInteger', fillInteger, &
            'fillReal', fillReal )
        endIf
      elseIf(fillArray) then
        if(fillInteger) then
          mqcVariable%integerArray(:) = arrayIntegerIn
        elseIf(fillReal) then

          if(DEBUG) then
            write(*,*)'   arrayRealIn(1) = ',arrayRealIn(1)
            write(*,*)'   length(arrayRealIn) = ',Size(arrayRealIn)
            write(*,*)'   length(mqcVariable%realArray(:)) = ',Size(mqcVariable%realArray(:))
            write(*,*)
          endIf
          mqcVariable%realArray(1:Size(arrayRealIn)) = arrayRealIn(1:Size(arrayRealIn))
          if(DEBUG) then
            write(*,*)' mqcVariable%realArray(1) = ',mqcVariable%realArray(1)
            write(*,*)
          endIf
        else
          call MQC_Error_L(  &
            'MQC_Variable_setVal: Confused filling array values.', 6, &
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
!
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
!     equal to an intrinsic real fortran matrix.
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
      logical::DEBUG=.false.
!
 1000 format(1x,'Converting an intrinsic real matrix to an mqcVar.')
 1100 format(3x,'myLength=',I8,'   dim 1=',I8,'   dim 2=',I8)
!
      DEBUG = DEBUG.or.MQC_Algebra2_DEBUG
!
!     Load the MQC variable.
!
      if(DEBUG) write(MQC_Algebra2_IOut,1000)
      myLength = SIZE(intrinsicIn)
      myDimensions(1) = SIZE(intrinsicIn,1)
      myDimensions(2) = SIZE(intrinsicIn,2)
      if(DEBUG) write(MQC_Algebra2_IOut,1100) myLength,  &
        myDimensions(1),myDimensions(2)
      call MQC_Variable_setVal(mqcVariable,arrayRealIn=RESHAPE(intrinsicIn,[myLength]),  &
        dimensions=myDimensions,DEBUGPRINT=DEBUG)
!
      return
      end subroutine MQC_Variable_setVal_intrinsicMatrix2mqcMatrix_real


!
!PROCEDURE MQC_Variable_setVal_intrinsicRank42mqcRank4_real
      subroutine MQC_Variable_setVal_intrinsicRank42mqcRank4_real(mqcVariable,intrinsicIn)
!
!     This subroutine is used to set an MQC Variable rank-4 tensor (full storage
!     form) equal to an intrinsic real fortran rank-4 tensor.
!
!
!     H. P. Hratchian, 2020.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable)::mqcVariable
      real(kind=real64),dimension(:,:,:,:),intent(in)::intrinsicIn
      integer(kind=int64)::myLength
      integer(kind=int64),dimension(4)::myDimensions
      logical::DEBUG=.false.
      integer(kind=int64),dimension(:),allocatable::tmpDimensions
      real(kind=real64),dimension(:),allocatable::tmpArray
!
 1000 format(1x,'Converting an intrinsic real rank-4 tensor to an mqcVar.')
 1100 format(3x,'myLength=',I8,'   dim 1=',I8,'   dim 2=',I8,  &
        '   dim 3=',I8,'   dim 4=',I8)
!
      DEBUG = DEBUG.or.MQC_Algebra2_DEBUG
!
!     Load the MQC variable.
!
      if(DEBUG) then
        write(MQC_Algebra2_IOut,1000)
        flush(MQC_Algebra2_IOut)
      endIf
      myLength = SIZE(intrinsicIn)
      myDimensions(1) = SIZE(intrinsicIn,1)
      myDimensions(2) = SIZE(intrinsicIn,2)
      myDimensions(3) = SIZE(intrinsicIn,3)
      myDimensions(4) = SIZE(intrinsicIn,4)
      if(DEBUG) write(MQC_Algebra2_IOut,1100) myLength,  &
        myDimensions(1),myDimensions(2),myDimensions(3),myDimensions(4)
      call MQC_Variable_setVal(mqcVariable,arrayRealIn=flatten(intrinsicIn),  &
        dimensions=myDimensions,DEBUGPRINT=DEBUG)
!
      return
      end subroutine MQC_Variable_setVal_intrinsicRank42mqcRank4_real


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
!PROCEDURE MQC_Variable_mqc2intrinsicInteger1Array
      subroutine MQC_Variable_mqc2intrinsicInteger1Array(intrinsicOut,mqcVariable)
!
!     This subroutine is used to set an (output) intrinsic integer vector equal
!     to an (input) MQC Variable vector.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),dimension(:),allocatable::intrinsicOut
      type(MQC_Variable),intent(in)::mqcVariable
!
!
!     Do the work...
!
      if(Rank(mqcVariable).ne.1)  &
        call mqc_error_I('Attemp to convert MQC_Variable array to '//  &
          'intrinsic SCALAR not allowed.',6,'Rank(mqcVariable)',  &
          Rank(mqcVariable))
      allocate(intrinsicOut(SIZE(mqcVariable,1)))
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        call mqc_error('MQC_Variable_mqc2intrinsicInteger1Array: '// &
          'REAL MQC_Variable type found when INTEGER required.',6)
      case(3)
        intrinsicOut = mqcVariable%integerArray(:)
      case default
        call mqc_error_I('MQC_Variable_mqc2intrinsicInteger1Array: '// &
          'Unknown MQC_Variable type found.',6, &
          'MQC_Variable_getTypeCode(mqcVariable)',  &
          MQC_Variable_getTypeCode(mqcVariable))
      end select
!
      return
      end subroutine MQC_Variable_mqc2intrinsicInteger1Array


!
!PROCEDURE MQC_Variable_mqc2intrinsicReal1Array
      subroutine MQC_Variable_mqc2intrinsicReal1Array(intrinsicOut,mqcVariable)
!
!     This subroutine is used to set an (output) intrinsic real vector equal to
!     an (input) MQC Variable vector.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:),allocatable::intrinsicOut
      type(MQC_Variable),intent(in)::mqcVariable
!
!
!     Do the work...
!
      if(Rank(mqcVariable).ne.1)  &
        call mqc_error_I('Attemp to convert MQC_Variable array to '//  &
          'intrinsic SCALAR not allowed.',6,'Rank(mqcVariable)',  &
          Rank(mqcVariable))
      if(allocated(intrinsicOut)) then
        if(SIZE(mqcVariable,1).ne.SIZE(intrinsicOut,1)) then
          deallocate(intrinsicOut)
          allocate(intrinsicOut(SIZE(mqcVariable,1)))
        endIf
      else
        allocate(intrinsicOut(SIZE(mqcVariable,1)))
      endIf
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        intrinsicOut = mqcVariable%realArray(:)
      case(3)
        intrinsicOut = mqcVariable%integerArray(:)
      case default
        call mqc_error_I('MQC_Variable_mqc2intrinsicReal1Array: '// &
          'Unknown MQC_Variable type found.',6, &
          'MQC_Variable_getTypeCode(mqcVariable)',  &
          MQC_Variable_getTypeCode(mqcVariable))
      end select
!
      return
      end subroutine MQC_Variable_mqc2intrinsicReal1Array


!
!PROCEDURE MQC_Variable_mqc2intrinsicInteger2Array
      subroutine MQC_Variable_mqc2intrinsicInteger2Array(intrinsicOut,mqcVariable)
!
!     This subroutine is used to set an (output) intrinsic integer matrix equal
!     to an (input) MQC Variable matrix.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      integer(kind=int64),dimension(:,:),allocatable,intent(out)::intrinsicOut
      type(MQC_Variable),intent(in)::mqcVariable
!
!
!     Do the work...
!
      if(Rank(mqcVariable).ne.2)  &
        call mqc_error_I('Attemp to convert MQC_Variable array to '//  &
          'intrinsic array with rank != 2 not allowed.',6,'Rank(mqcVariable)',  &
          Rank(mqcVariable))
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        call mqc_error('MQC_Variable_mqc2intrinsicInteger2Array: '// &
          'REAL MQC_Variable type found when INTEGER required.',6)
      case(3)
        allocate(intrinsicOut(SIZE(mqcVariable,1),SIZE(mqcVariable,2)))
        intrinsicOut = RESHAPE(mqcVariable%integerArray(:),  &
          mqcVariable%dimensions(1:2))
      case default
        call mqc_error_I('MQC_Variable_mqc2intrinsicInteger2Array: '// &
          'Unknown MQC_Variable type found.',6, &
          'MQC_Variable_getTypeCode(mqcVariable)',  &
          MQC_Variable_getTypeCode(mqcVariable))
      end select
!
      return
      end subroutine MQC_Variable_mqc2intrinsicInteger2Array


!
!PROCEDURE MQC_Variable_mqc2intrinsicReal2Array
      subroutine MQC_Variable_mqc2intrinsicReal2Array(intrinsicOut,mqcVariable)
!
!     This subroutine is used to set an (output) intrinsic real matrix equal to
!     an (input) MQC Variable matrix.
!
!
!     H. P. Hratchian, 2019.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:,:),allocatable,intent(out)::intrinsicOut
      type(MQC_Variable),intent(in)::mqcVariable
!
!
!     Do the work...
!
      if(Rank(mqcVariable).ne.2)  &
        call mqc_error_I('Attemp to convert MQC_Variable array to '//  &
          'intrinsic array with rank != 2 not allowed.',6,'Rank(mqcVariable)',  &
          Rank(mqcVariable))
      allocate(intrinsicOut(SIZE(mqcVariable,1),SIZE(mqcVariable,2)))
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        intrinsicOut = RESHAPE(mqcVariable%realArray(:),  &
          mqcVariable%dimensions(1:2))
      case(3)
        intrinsicOut = RESHAPE(mqcVariable%integerArray(:),  &
          mqcVariable%dimensions(1:2))
      case default
        call mqc_error_I('MQC_Variable_mqc2intrinsicReal2Array: '// &
          'Unknown MQC_Variable type found.',6, &
          'MQC_Variable_getTypeCode(mqcVariable)',  &
          MQC_Variable_getTypeCode(mqcVariable))
      end select
!
      return
      end subroutine MQC_Variable_mqc2intrinsicReal2Array


!
!PROCEDURE MQC_Variable_mqc2intrinsicReal4Array
      subroutine MQC_Variable_mqc2intrinsicReal4Array(intrinsicOut,mqcVariable)
!
!     This subroutine is used to set an (output) intrinsic real rank-4 array equal to
!     an (input) MQC Variable.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      real(kind=real64),dimension(:,:,:,:),allocatable,intent(out)::intrinsicOut
      type(MQC_Variable),intent(in)::mqcVariable
      real(kind=real64),dimension(:),allocatable::tmpArray
!
!
!     Do the work...
!
      if(Rank(mqcVariable).ne.4)  &
        call mqc_error_I('Attemp to convert MQC_Variable array to '//  &
          'intrinsic array with rank != 4 not allowed.',6,'Rank(mqcVariable)',  &
          Rank(mqcVariable))
      allocate(intrinsicOut(SIZE(mqcVariable,1),SIZE(mqcVariable,2),  &
        SIZE(mqcVariable,3),SIZE(mqcVariable,4)))
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        intrinsicOut = RESHAPE(mqcVariable%realArray(:),  &
          mqcVariable%dimensions(1:4))
      case(3)
        intrinsicOut = RESHAPE(mqcVariable%integerArray(:),  &
          mqcVariable%dimensions(1:4))
      case default
        call mqc_error_I('MQC_Variable_mqc2intrinsicReal4Array: '// &
          'Unknown MQC_Variable type found.',6, &
          'MQC_Variable_getTypeCode(mqcVariable)',  &
          MQC_Variable_getTypeCode(mqcVariable))
      end select
!
      return
      end subroutine MQC_Variable_mqc2intrinsicReal4Array


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


!
!PROCEDURE MQC_Variable_Matrix_Transpose
      function MQC_Variable_Matrix_Transpose(mqcVariable) result(mqcVariableOut)
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
      type(MQC_Variable)::mqcVariableOut
!
!
!     Do the work...
!
      if(RANK(mqcVariable).ne.2) call mqc_error_i('Transpose invoked on non-matrix. RANK=',  &
        6,'RANK=',RANK(mqcVariable))
      if(TRIM(mqcVariable%storageFormat).ne.'FULL')   &
        call mqc_error_a('Transpose invoked on non-FULL storage format matrix NYI.',6,'Storage form:',  &
        TRIM(mqcVariable%storageFormat))
      select case(MQC_Variable_getTypeCode(mqcVariable))
      case(2)
        mqcVariableOut = Transpose(RESHAPE(mqcVariable%realArray(:),mqcVariable%dimensions(1:2)))
      case(3)
        mqcVariableOut = Transpose(RESHAPE(mqcVariable%integerArray(:),mqcVariable%dimensions(1:2)))
      case default
        call mqc_error('Transpose: Erroneous typing found.')
      end select
!
      return
      end function MQC_Variable_Matrix_Transpose


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
      integer(kind=int64)::typeCode1,typeCode2,sumScalarInteger
      integer(kind=int64),dimension(:),allocatable::sumVectorInteger
      real(kind=real64),dimension(:),allocatable::sumVectorReal
      real(kind=real64)::sumScalarReal
!
!
!     Start by ensuring the two MQC Variables are conformable.
!
      if(.not.MQC_Variable_isConformable(mqcVariable1,mqcVariable2))  &
        call mqc_error_l('MQC Addition only allowed for conformable values.', 6, &
        'MQC_Variable_isConformable(mqcVariable1,mqcVariable2)', &
        MQC_Variable_isConformable(mqcVariable1,mqcVariable2) )
      if(mqcVariable1%storageFormat.ne.mqcVariable2%storageFormat)  &
        call mqc_error_a('MQC Addition only allowed for same storage types.', 6, &
        'mqcVariable1 storage type:',TRIM(mqcVariable1%storageFormat),  &
        'mqcVariable2 storage type:',TRIM(mqcVariable2%storageFormat))
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
      select case(mqcVariable1%rank)
      case(0)
        select case(typeCode1*10 + typeCode2)
        case(22)
          sumScalarReal = mqcVariable1%realArray(1) + mqcVariable2%realArray(1)
          call MQC_Variable_setVal(mqcVariableOut,scalarRealIn=sumScalarReal,  &
            dimensions=mqcVariable1%dimensions)
        case(23)
          sumScalarReal = mqcVariable1%realArray(1) + mqcVariable2%integerArray(1)
          call MQC_Variable_setVal(mqcVariableOut,scalarRealIn=sumScalarReal,  &
            dimensions=mqcVariable1%dimensions)
        case(32)
          sumScalarReal = mqcVariable1%integerArray(1) + mqcVariable2%realArray(1)
          call MQC_Variable_setVal(mqcVariableOut,scalarRealIn=sumScalarReal,  &
            dimensions=mqcVariable1%dimensions)
        case(33)
          sumScalarInteger = mqcVariable1%integerArray(1) + mqcVariable2%integerArray(1)
          call MQC_Variable_setVal(mqcVariableOut,scalarIntegerIn=sumScalarInteger,  &
            dimensions=mqcVariable1%dimensions)
        case default
          call mqc_error_I('MQC_Variable_Addition: Combination of var1 and var 2 types is UNKNOWN.', 6, &
              'typeCode1*10', typeCode1*10, &
              'typeCode2', typeCode2)
        end select
      case(1:)
        select case(typeCode1*10 + typeCode2)
        case(22)
          sumVectorReal = mqcVariable1%realArray + mqcVariable2%realArray
          call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=sumVectorReal,  &
            dimensions=mqcVariable1%dimensions)
        case(23)
          sumVectorReal = mqcVariable1%realArray + mqcVariable2%integerArray
          call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=sumVectorReal,  &
            dimensions=mqcVariable1%dimensions)
        case(32)
          sumVectorReal = mqcVariable1%integerArray(:) + mqcVariable2%realArray(:)
          call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=sumVectorReal,  &
            dimensions=mqcVariable1%dimensions)
        case(33)
          sumVectorInteger = mqcVariable1%integerArray + mqcVariable2%integerArray
          call MQC_Variable_setVal(mqcVariableOut,arrayIntegerIn=sumVectorInteger,  &
            dimensions=mqcVariable1%dimensions)
        case default
          call mqc_error_I('MQC_Variable_Addition: Combination of var1 and var 2 types is UNKNOWN.', 6, &
              'typeCode1*10', typeCode1*10, &
              'typeCode2', typeCode2)
        end select
      case default
        call mqc_error('MQC_Variable_Addition is confused!')
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
      integer(kind=int64)::typeCode1,typeCode2,differenceScalarInteger
      integer(kind=int64),dimension(:),allocatable::differenceVectorInteger
      real(kind=real64),dimension(:),allocatable::differenceVectorReal
      real(kind=real64)::differenceScalarReal
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
      select case(mqcVariable1%rank)
      case(0)
        select case(typeCode1*10 + typeCode2)
        case(22)
          differenceScalarReal = mqcVariable1%realArray(1) - mqcVariable2%realArray(1)
          call MQC_Variable_setVal(mqcVariableOut,scalarRealIn=differenceScalarReal,  &
            dimensions=mqcVariable1%dimensions)
        case(23)
          differenceScalarReal = mqcVariable1%realArray(1) - mqcVariable2%integerArray(1)
          call MQC_Variable_setVal(mqcVariableOut,scalarRealIn=differenceScalarReal,  &
            dimensions=mqcVariable1%dimensions)
        case(32)
          differenceScalarReal = mqcVariable1%integerArray(1) - mqcVariable2%realArray(1)
          call MQC_Variable_setVal(mqcVariableOut,scalarRealIn=differenceScalarReal,  &
            dimensions=mqcVariable1%dimensions)
        case(33)
          differenceScalarInteger = mqcVariable1%integerArray(1) - mqcVariable2%integerArray(1)
          call MQC_Variable_setVal(mqcVariableOut,scalarIntegerIn=differenceScalarInteger,  &
            dimensions=mqcVariable1%dimensions)
        case default
          call mqc_error_I('MQC_Variable_subtraction: Combination of var1 and var 2 types is UNKNOWN.', 6, &
              'typeCode1*10', typeCode1*10, &
              'typeCode2', typeCode2)
        end select
      case(1:)
        select case(typeCode1*10 + typeCode2)
        case(22)
          differenceVectorReal = mqcVariable1%realArray - mqcVariable2%realArray
          call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=differenceVectorReal,  &
            dimensions=mqcVariable1%dimensions)
        case(23)
          differenceVectorReal = mqcVariable1%realArray - mqcVariable2%integerArray
          call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=differenceVectorReal,  &
            dimensions=mqcVariable1%dimensions)
        case(32)
          differenceVectorReal = mqcVariable1%integerArray(:) - mqcVariable2%realArray(:)
          call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=differenceVectorReal,  &
            dimensions=mqcVariable1%dimensions)
        case(33)
          differenceVectorInteger = mqcVariable1%integerArray - mqcVariable2%integerArray
          call MQC_Variable_setVal(mqcVariableOut,arrayIntegerIn=differenceVectorInteger,  &
            dimensions=mqcVariable1%dimensions)
        case default
          call mqc_error_I('MQC_Variable_subtraction: Combination of var1 and var 2 types is UNKNOWN.', 6, &
              'typeCode1*10', typeCode1*10, &
              'typeCode2', typeCode2)
        end select
      case default
        call mqc_error('MQC_Variable_Subtraction is confused!')
      end select
!hph-

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
!     Now, do the mulitplication and set the output function value accordingly.
!
      select case(typeCode1*10 + typeCode2)
      case(22)
        productVectorReal = mqcVariable1%realArray * mqcVariable2%realArray
        call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=productVectorReal,  &
          dimensions=mqcVariable1%dimensions)
      case(23)
        productVectorReal = mqcVariable1%realArray * mqcVariable2%integerArray
        call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=productVectorReal,  &
          dimensions=mqcVariable1%dimensions)
      case(32)
        productVectorReal = mqcVariable1%integerArray(:) * mqcVariable2%realArray(:)
        call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=productVectorReal,  &
          dimensions=mqcVariable1%dimensions)
      case(33)
        productVectorInteger = mqcVariable1%integerArray * mqcVariable2%integerArray
        call MQC_Variable_setVal(mqcVariableOut,arrayIntegerIn=productVectorInteger,  &
          dimensions=mqcVariable1%dimensions)
      case default
        call mqc_error_I('MQC_Variable_Multiplication: Combination of var1 and var 2 types is UNKNOWN.', 6, &
             'typeCode1*10', typeCode1*10, &
             'typeCode2', typeCode2 )
      end select
!hph+
!!
!!     Ensure the other settings defining the MQC_Variable object info match the
!!     input arguments.
!!
!      mqcVariableOut%rank = mqcVariable1%rank
!      mqcVariableOut%dimensions = mqcVariable1%dimensions
!hph-

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
        call mqc_error_I(  &
          'MQC_Variable_Division: Result is of UNKNOWN type.', 6, &
          'typeCode1', typeCode1, &
          'typeCode2', typeCode2 )
      end select
!
!     Now, do the division and set the output function value accordingly.
!
      select case(typeCode1*10 + typeCode2)
      case(22)
        quotientVectorReal = mqcVariable1%realArray / mqcVariable2%realArray
        call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=quotientVectorReal,  &
          dimensions=mqcVariable1%dimensions)
      case(23)
        quotientVectorReal = mqcVariable1%realArray / mqcVariable2%integerArray
        call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=quotientVectorReal,  &
          dimensions=mqcVariable1%dimensions)
      case(32)
        quotientVectorReal = mqcVariable1%integerArray(:) / mqcVariable2%realArray(:)
        call MQC_Variable_setVal(mqcVariableOut,arrayRealIn=quotientVectorReal,  &
          dimensions=mqcVariable1%dimensions)
      case(33)
        quotientVectorInteger = mqcVariable1%integerArray / mqcVariable2%integerArray
        call MQC_Variable_setVal(mqcVariableOut,arrayIntegerIn=quotientVectorInteger,  &
          dimensions=mqcVariable1%dimensions)
      case default
        call mqc_error_I(  &
          'MQC_Variable_Division: Combination of var1 and var 2 '//  &
          'types is UNKNOWN.',6,'typeCode1*10',typeCode1*10, &
             'typeCode2', typeCode2)
      end select
!
      return
      end function MQC_Variable_Division


!
!PROCEDURE MQC_Variable_IntegerPower
      function MQC_Variable_IntegerPower(mqcVariable,intPower) result(mqcVariableOut)
!
!     This function raises all elements of the array <mqcVariable> to the power
!     <intPower>. Care is taken to leave zero values alone when <intPower> is
!     negative.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      integer(kind=int64),intent(in)::intPower
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode,nDim,i
!
!     Get the type-code that indicates whether <mqcVariable> is real or integer.
!
      typeCode = MQC_Variable_getTypeCode(mqcVariable)
      if(typeCode.eq.0)  &
        call mqc_error_I('MQC_Variable_IntegerPower: mqcVariable is of UNKONWN type.',6, &
        'typeCode',typeCode)
!
!     Based on the sign of <intPower>, we now do the work.
!
      mqcVariableOut = mqcVariable
      if(intPower.ge.0) then
        select case(typeCode)
        case(2)
          mqcVariableOut%realArray = mqcVariableOut%realArray**intPower
        case(3)
          mqcVariableOut%integerArray = mqcVariableOut%integerArray**intPower
        case default
          call mqc_error_I(  &
          'MQC_Variable_IntPower: mqcVariable is of UNKNOWN type.',6, &
          'typeCode',typeCode)
        end select
      else
        select case(typeCode)
        case(2)
          nDim = SIZE(mqcVariableOut%realArray)
          do i = 1,nDim
            if(ABS(mqcVariableOut%realArray(i)).gt.MQC_small)  &
              mqcVariableOut%realArray(i) = mqcVariableOut%realArray(i)**intPower
          endDo
        case(3)
          nDim = SIZE(mqcVariableOut%integerArray)
          do i = 1,nDim
            if(ABS(mqcVariableOut%integerArray(i)).gt.MQC_small)  &
              mqcVariableOut%integerArray(i) = mqcVariableOut%integerArray(i)**intPower
          endDo
        case default
          call mqc_error_I(  &
          'MQC_Variable_IntPower: mqcVariable is of UNKNOWN type.',6, &
          'typeCode',typeCode)
        end select
      endIf
!
      return
      end function MQC_Variable_IntegerPower


!
!PROCEDURE MQC_Variable_TrimZero
      Subroutine mqc_variable_trimZero(mqcVariable,small)
!
!     This subroutine accepts an MQC Variable (mqcVariable) and sets all entries that
!     have absolute values below a threshold to float(0).
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
      class(MQC_Variable),intent(inOut)::mqcVariable
      real(kind=real64),optional,intent(in)::small
      real(kind=real64)::mySmall
!
!     Do the requested work by calling mqc_vectorTrimZero from MQC_General.
!
      if(mqcVariable%getType().ne.'REAL')  &
        call mqc_error('mqc_variable_trimZero: Only real arrays are allowed.')
      if(PRESENT(small)) then
        mySmall = small
      else
        mySmall = MQC_small
      endIf
      call mqc_vectorTrimZero(mqcVariable%realArray,small=mySmall)
!
      return
      end subroutine mqc_variable_trimZero


!
!PROCEDURE MQC_Variable_RealPower
      function MQC_Variable_RealPower(mqcVariable,realPower) result(mqcVariableOut)
!
!     This function raises all elements of the array <mqcVariable> to the power
!     <realPower>. Care is taken to leave zero values alone when <realPower> is
!     negative.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      real(kind=real64),intent(in)::realPower
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode,nDim,i
!
!     Get the type-code that indicates whether <mqcVariable> is real or integer.
!
      typeCode = MQC_Variable_getTypeCode(mqcVariable)
      if(typeCode.eq.0)  &
        call mqc_error_I('MQC_Variable_realPower: mqcVariable is of UNKONWN type.',6, &
        'typeCode',typeCode)
!
!     Based on the sign of <realPower>, we now do the work.
!
      mqcVariableOut = mqcVariable
      if(realPower.ge.float(0)) then
        select case(typeCode)
        case(2)
          mqcVariableOut%realArray = mqcVariableOut%realArray**realPower
        case(3)
          mqcVariableOut%integerArray = mqcVariableOut%integerArray**realPower
        case default
          call mqc_error_I(  &
          'MQC_Variable_realPower: mqcVariableOut is of UNKNOWN type.',6, &
          'typeCode',typeCode)
        end select
      else
        select case(typeCode)
        case(2)
          nDim = SIZE(mqcVariableOut%realArray)
          do i = 1,nDim
            if(ABS(mqcVariableOut%realArray(i)).gt.MQC_small)  &
              mqcVariableOut%realArray(i) = mqcVariableOut%realArray(i)**realPower
          endDo
        case(3)
          nDim = SIZE(mqcVariableOut%integerArray)
          do i = 1,nDim
            if(ABS(mqcVariableOut%integerArray(i)).gt.MQC_small)  &
              mqcVariableOut%integerArray(i) = mqcVariableOut%integerArray(i)**realPower
          endDo
        case default
          call mqc_error_I(  &
          'MQC_Variable_realPower: mqcVariable is of UNKNOWN type.',6, &
          'typeCode',typeCode)
        end select
      endIf
!
      return
      end function MQC_Variable_RealPower


!
!PROCEDURE MQC_Variable_Array_Sum
      function MQC_Variable_Array_Sum(mqcVariableIn,iDimension)  &
        result(mqcVariableOut)
!
!     This function carries out a summation of all elements in an array or
!     possibly sums only across elements of one specified dimension.
!     Specifically, INPUT dummy variable <iDimension> is optional. If it is
!     NOT sent, then this function returns the sum of all elements in INPUT
!     dummy variable <mqcVariableIn>.
!
!     The output of this function is a new MQC_Variable with the correct
!     rank/dimensions.
!
!
!     H. P. Hratchian, 2018, 2019.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariableIn
      integer(kind=int64),optional::iDimension
      type(MQC_Variable)::mqcVariableOut
      real(kind=real64)::scalarReal
      real(kind=real64),dimension(:),allocatable::tempReal1
      integer(kind=int64)::itmp,myDimension,myRank,resultRank,myLen,  &
        resultLen,scalarInteger
      integer(kind=int64),dimension(:),allocatable::tempInteger1
      real(kind=real64),dimension(2,3)::junkMat

!
!
!     Ensure that iDimension is valid and determine the rank of the result.
!
      myRank = mqcVariableIn%getRank()
      if(PRESENT(iDimension)) then
        myDimension = iDimension
        resultRank = myRank-1
      else
        myDimension = 0
        resultRank = 0
      endIf
      if(myDimension.gt.myRank)  &
        call mqc_error_I('MQC_Variable_Array_Sum: myDimension > RANK.', &
        6,'myDimension=',myDimension,'RANK=',myRank)
!
!     Ensure we are working with full storage array forms. Other storage
!     options will be coded in the future.
!
      if(.not.TRIM(mqcVariableIn%storageFormat).eq.'FULL')  &
        call mqc_error('SUM of MQC Variables is only available for '//  &
          'FULL storage arrays.')
!
!     Take care of the special case where <mqcVariableIn> is a scalar.
!
      if(myRank.eq.0) then
        mqcVariableOut = mqcVariableIn
        return
      endIf
      mqcVariableOut = 0
      select case(resultRank)
      case(0)
        select case(MQC_Variable_getTypeCode(mqcVariableIn))
        case(2)
          scalarReal = SUM(mqcVariableIn%realArray(:))
          mqcVariableOut = scalarReal
        case(3)
          scalarInteger = SUM(mqcVariableIn%integerArray(:))
          mqcVariableOut = scalarInteger
        case default
          call mqc_error_I('MQC_Variable_Array_Sum: '// &
            'Unknown MQC_Variable type found.',6, &
            'MQC_Variable_getTypeCode(mqcVariable)',  &
            MQC_Variable_getTypeCode(mqcVariableIn))
        end select
      case(1)
        myLen = MQC_Variable_getLength(mqcVariableIn)
        if(myRank.ne.2)  &
          call mqc_error('MQC_Variable_Array_Sum: '// &
            'Confusion with the rank of input array.')
        select case(MQC_Variable_getTypeCode(mqcVariableIn))
        case(2) ! REALS
          select case(myDimension)
          case(1)
            ALLOCATE(tempReal1(Size(mqcVariableIn,2)))
            tempReal1 = SUM(  &
              RESHAPE(mqcVariableIn%realArray,  &
              [Size(mqcVariableIn,1),Size(mqcVariableIn,2)]),1)
            call MQC_Variable_setVal(mqcVariableOut,  &
              arrayRealIn=tempReal1,dimensions=[Size(tempReal1)])
            DEALLOCATE(tempReal1)
          case(2)
            ALLOCATE(tempReal1(Size(mqcVariableIn,1)))
            tempReal1 = SUM(  &
              RESHAPE(mqcVariableIn%realArray,  &
              [Size(mqcVariableIn,1),Size(mqcVariableIn,2)]),2)
            call MQC_Variable_setVal(mqcVariableOut,  &
              arrayRealIn=tempReal1,dimensions=[Size(tempReal1)])
            DEALLOCATE(tempReal1)
          case default
            call mqc_error('MQC_Variable_Array_Sum: '// &
              'Attempt to contract a matrix using dimension < 1 or > 2.')
          end select
        case(3) ! INTEGERS
          call mqc_error('MQC_Variable_Array_Sum: '// &
            'SUM with integers NYI.')
        case default
          call mqc_error('MQC_Variable_Array_Sum: '// &
            'Unknown array type (Real/Integer/etc.).')
        end select
      case default
          call mqc_error('MQC_Variable_Array_Sum: '// &
            'SUM with resultant of rank 2 or greater NYI.')
      end select
!
      return
      end function MQC_Variable_Array_Sum


!
!PROCEDURE MQC_Variable_SubMatrix
      function MQC_Variable_SubMatrix(mqcVariableIn,newRange1,newRange2)  &
        result(mqcVariableOut)
!
!     This function forms a sub-matrix from mqcVariableIn. The output matrix has
!     new ranges in the first and second dimension of <newRange1> and
!     <newRange2>, respectively.
!
!     Dummy argument iDimension is optional. If NOT sent, it defaults to 1.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      implicit none
      class(MQC_Variable),intent(in)::mqcVariableIn
      integer(kind=int64),dimension(2),optional::newRange1,newRange2
      type(MQC_Variable)::mqcVariableOut
!
      integer(kind=int64)::start1,end1,start2,end2
      integer(kind=int64),dimension(:,:),allocatable::intMatrixTmp
      real(kind=real64),dimension(:,:),allocatable::realMatrixTmp
!
!
!     Start by checking that we're working with a matrix. If not, die.
!
      if(RANK(mqcVariableIn).ne.2)  &
        call mqc_error_i('Error in MQC_Variable_Submatrix: Incorrect rank found.',  &
        6,'RANK(mqcVariableIn)=',RANK(mqcVariableIn))
!
!     Get the subarray starting and ending points figured out for the 2
!     dimensions.
!
      if(PRESENT(newRange1)) then
        start1 = newRange1(1)
        end1   = newRange1(2)
      else
        start1 = 1
        end1   = SIZE(mqcVariableIn,1)
      endIf
      if(PRESENT(newRange2)) then
        start2 = newRange2(1)
        end2   = newRange2(2)
      else
        start2 = 1
        end2   = SIZE(mqcVariableIn,2)
      endIf
!
!     Form the sub-matrix.
!
      select case(TRIM(MQC_Variable_getType(mqcVariableIn)))
      case('REAL')
        Allocate(realMatrixTmp(SIZE(mqcVariableIn,1),SIZE(mqcVariableIn,2)))
        realMatrixTmp = mqcVariableIn
        mqcVariableOut = realMatrixTmp(start1:end1,start2:end2)
        DeAllocate(realMatrixTmp)
      case('INTEGER')
        Allocate(intMatrixTmp(SIZE(mqcVariableIn,1),SIZE(mqcVariableIn,2)))
        intMatrixTmp = mqcVariableIn
        mqcVariableOut = intMatrixTmp(start1:end1,start2:end2)
        DeAllocate(intMatrixTmp)
      case default
        call mqc_error_a('MQC_Varaible_SubMatrix:: Unknown MQC variable type.',6, &
          'TRIM(MQC_Variable_getType(mqcVariableIn))',  &
          TRIM(MQC_Variable_getType(mqcVariableIn)))
      end select
!
      return
      end function MQC_Variable_SubMatrix

!
!PROCEDURE mqc_variable_matrixTrace
      function mqc_variable_matrixTrace(mqcVariable) result(trace)
!
!     This function returns the trace of matrix. For now, only real matrices are
!     allowed.
!
!
!     H. P. Hratchian, 2025
!
      implicit none
      class(MQC_Variable),intent(in)::mqcVariable
      real(kind=real64)::trace
!
!     Call the MQC_General routine and then return the trace.
!
      if(mqcVariable%getType().ne.'REAL')  &
        call mqc_error('mqc_variable_matrixTrace: Only real arrays are allowed.')
      trace = mqc_matrixTrace(RESHAPE(mqcVariable%realArray,  &
        [Size(mqcVariable,1),Size(mqcVariable,2)]))
!
      return
      end function mqc_variable_matrixTrace

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
!
!     Calculate the dot product.
!
      mqcVariableOut = MQC_Variable_Contraction_Full(mqcVariable1,mqcVariable2)
!
      return
      end function MQC_Variable_dot_product


!
!PROCEDURE MQC_Variable_outer_product
      function MQC_Variable_outer_product(mqcVariable1,mqcVariable2) result(mqcVariableOut)
!
!     This function solves the outer product of two vectors stored in the input
!     dummy MQC_Variables arguments <mqcVariable1> and <mqcVariable2>. If these
!     two arguments do not correspond to vectors, this function fails.
!
!
!     H. P. Hratchian, 2024.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode1,typeCode2,len1,len2
      integer(kind=int64),dimension(:),allocatable::vectorInteger1,  &
        vectorInteger2
      integer(kind=int64),dimension(:,:),allocatable::matrixInteger
      real(kind=real64),dimension(:),allocatable::vectorReal1,  &
        vectorReal2
      real(kind=real64),dimension(:,:),allocatable::matrixReal
!
!
!     Start by ensuring the two MQC Variables are vectors and are conformable.
!
      if(mqcVariable1%getRank().ne.1.or.mqcVariable2%getRank().ne.1)  &
        call mqc_error_i('outer_product must be between two vectors.',6, &
        'RANK(mqcVariable1)',mqcVariable1%getRank(),  &
        'RANK(mqcVariable2)',mqcVariable2%getRank())
!
!     Determine the type codes for the two vectors and ensure they're allowable.
!
      typeCode1 = MQC_Variable_getTypeCode(mqcVariable1)
      typeCode2 = MQC_Variable_getTypeCode(mqcVariable2)
      if(typeCode1.eq.0)  &
        call mqc_error_I('Outer product: Var1 is of UNKONWN type.', 6, &
          'typeCode1', typeCode1)
      if(typeCode2.eq.0)  &
        call mqc_error_I('Outer product: Var2 is of UNKONWN type.', 6, &
          'typeCode2', typeCode2)
!
!     Now, do the outer product and set the output function value accordingly.
!
      len1 = SIZE(mqcVariable1)
      len2 = SIZE(mqcVariable2)
      select case(typeCode1*10 + typeCode2)
      case(22)
        vectorReal1 = mqcVariable1
        vectorReal2 = mqcVariable2
        matrixReal = matmul(RESHAPE(vectorReal1,[ len1,1 ]),  &
          RESHAPE(vectorReal2,[ 1,len2 ]))
        mqcVariableOut = matrixReal
      case(23)
        vectorReal1 = mqcVariable1
        vectorInteger2 = mqcVariable2
        matrixReal = matmul(RESHAPE(vectorReal1,[ len1,1 ]),  &
          RESHAPE(vectorInteger2,[ 1,len2 ]))
        mqcVariableOut = matrixReal
      case(32)
        vectorInteger1 = mqcVariable1
        vectorReal2 = mqcVariable2
        matrixReal = matmul(RESHAPE(vectorInteger1,[ len1,1 ]),  &
          RESHAPE(vectorReal2,[ 1,len2 ]))
        mqcVariableOut = matrixReal
      case(33)
        vectorInteger1 = mqcVariable1
        vectorInteger2 = mqcVariable2
        matrixInteger = matmul(RESHAPE(vectorInteger1,[ len1,1 ]),  &
          RESHAPE(vectorInteger2,[ 1,len2 ]))
        mqcVariableOut = matrixInteger
      case default
        call mqc_error_I('Outer product: Combination of var1 and var 2 types is UNKNOWN.', 6, &
             'typeCode1', typeCode1, &
             'typeCode2', typeCode2 )
      end select
!
      return
      end function MQC_Variable_outer_product


!
!PROCEDURE MQC_Variable_MatrixVector
      function MQC_Variable_MatrixVector(mqcMatrix,mqcVector) result(mqcVariableOut)
!
!     This function carries out a matrix-vector product using two input
!     MQC_Variables. The output is an MQC_Variable type.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcMatrix,mqcVector
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode1,typeCode2
      integer(kind=int64),dimension(:),allocatable::productVectorInteger
      real(kind=real64),dimension(:),allocatable::productVectorReal
!
!
!     Start by ensuring we have a matrix and a vector. Then check that the two
!     MQC Variables are conformable.
!
      if(Rank(mqcMatrix).ne.2) call mqc_error('MatrixVector: Matrix is NOT rank=2.')
      if(Rank(mqcVector).ne.1) call mqc_error('MatrixVector: Vector is NOT rank=1.')
      if(Size(mqcMatrix,2).ne.Size(mqcVector))  &
        call mqc_error('MatrixVector: Matrix/Vector not conformable.')
!
!     Get the type codes figured out.
!
      typeCode1 = MQC_Variable_getTypeCode(mqcMatrix)
      typeCode2 = MQC_Variable_getTypeCode(mqcVector)
      if(typeCode1.eq.0)  &
        call mqc_error_I('MatrixVector: Matrix is of UNKONWN type.', 6, &
        'typeCode1', typeCode1)
      if(typeCode2.eq.0)  &
        call mqc_error_I('MatrixVector: Vector is of UNKONWN type.', 6, &
        'typeCode2' ,typeCode2);
      select case(MIN(typeCode1,typeCode2))
      case(2)
        allocate(productVectorReal(SIZE(mqcMatrix,1)))
      case(3)
        allocate(productVectorInteger(SIZE(mqcMatrix,1)))
      case default
        call mqc_error_i('MatrixVector: Result is of UNKNOWN type.', 6, &
             'typeCode1', typeCode1, &
             'typeCode2', typeCode2 )
      end select
!
!     Now, do the matrix-vector product and set the output function value
!     accordingly.
!
      select case(typeCode1*10 + typeCode2)
      case(22)
        productVectorReal = MatMul(  &
          RESHAPE(mqcMatrix%realArray(:),mqcMatrix%dimensions(1:2)),  &
          mqcVector%realArray(:))
        mqcVariableOut = productVectorReal
      case(23)
        productVectorReal = MatMul(  &
          RESHAPE(mqcMatrix%realArray(:),mqcMatrix%dimensions(1:2)),  &
          mqcVector%integerArray(:))
        mqcVariableOut = productVectorReal
      case(32)
        productVectorReal = MatMul(  &
          RESHAPE(mqcMatrix%integerArray(:),mqcMatrix%dimensions(1:2)),  &
          mqcVector%realArray(:))
        mqcVariableOut = productVectorReal
      case(33)
        productVectorInteger = MatMul(  &
          RESHAPE(mqcMatrix%integerArray(:),mqcMatrix%dimensions(1:2)),  &
          mqcVector%integerArray(:))
        mqcVariableOut = productVectorInteger
      case default
        call mqc_error_I('MQC_MatrixVector: Combination of var1 and var 2 types is UNKNOWN.', 6, &
             'typeCode1*10', typeCode1*10, &
             'typeCode2', typeCode2 )
      end select
!
      return
      end function MQC_Variable_MatrixVector


!
!PROCEDURE MQC_Variable_Matrix_Multiplication
      function MQC_Variable_Matrix_Multiplication(mqcVariable1,mqcVariable2) result(mqcVariableOut)
!
!     This function carries out matrix multiplication using two input
!     MQC_Variables. The output is an MQC_Variable type. This function provides
!     the same functionality as the fortran standard provides for intrinsic
!     MatMul.
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable1,mqcVariable2
      type(MQC_Variable)::mqcVariableOut
      integer(kind=int64)::typeCode1,typeCode2
      integer(kind=int64),dimension(:),allocatable::productVectorInteger
      integer(kind=int64),dimension(:,:),allocatable::productMatrixInteger
      real(kind=real64),dimension(:),allocatable::productVectorReal
      real(kind=real64),dimension(:,:),allocatable::productMatrixReal
!
!
!     Start by ensuring the two MQC Variables are matrices and have conformable
!     dimensions.
!
      if(Rank(mqcVariable1).ne.2)  &
        call mqc_error('MQC Matrix Multiplication: mqcVariable 1 is NOT a matrix.')
      if((Rank(mqcVariable2).ne.2).and.(Rank(mqcVariable2).ne.1))  &
        call mqc_error('MQC Matrix Multiplication: mqcVariable 2 is NOT a matrix or vector.')
      if(Size(mqcVariable1,2).ne.Size(mqcVariable2,1))  &
        call mqc_error('MQC Matrix Multiplication: Matrices do not have matching inner-dimensions.')
!
!     Get the right combination of integer/real and build a temporary result
!     array. This will allow uses of this function to work for C = MatMul(A,B),
!     C = MatMul(C,B), and C = MatMul(A,C).
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
        if(RANK(mqcVariable2).eq.2) then
          allocate(productMatrixReal(SIZE(mqcVariable1,1),SIZE(mqcVariable2,2)))
        else
          allocate(productVectorReal(SIZE(mqcVariable1,1)))
        endIf
      case(3)
        if(RANK(mqcVariable2).eq.2) then
          allocate(productMatrixInteger(SIZE(mqcVariable1,1),SIZE(mqcVariable2,2)))
        else
          allocate(productVectorInteger(SIZE(mqcVariable1,1)))
        endIf
      case default
        call mqc_error_i('MQC_Variable_Multiplication: Result is of UNKNOWN type.', 6, &
             'typeCode1', typeCode1, &
             'typeCode2', typeCode2 )
      end select
!
!     Now, do the multiplication and set the output function value accordingly.
!     Note that this is done in three blocks...matrix-matrix, matrix-vector, and
!     vector-matrix (NYI).
!
      if(RANK(mqcVariable1).eq.2.and.RANK(mqcVariable2).eq.2) then
        select case(typeCode1*10 + typeCode2)
        case(22)
          productMatrixReal = MatMul(  &
            RESHAPE(mqcVariable1%realArray(:),mqcVariable1%dimensions(1:2)),  &
            RESHAPE(mqcVariable2%realArray(:),mqcVariable2%dimensions(1:2)))
          mqcVariableOut = productMatrixReal
        case(23)
          productMatrixReal = MatMul(  &
            RESHAPE(mqcVariable1%realArray(:),mqcVariable1%dimensions(1:2)),  &
            RESHAPE(mqcVariable2%integerArray(:),mqcVariable2%dimensions(1:2)))
          mqcVariableOut = productMatrixReal
        case(32)
          productMatrixReal = MatMul(  &
            RESHAPE(mqcVariable1%integerArray(:),mqcVariable1%dimensions(1:2)),  &
            RESHAPE(mqcVariable2%realArray(:),mqcVariable2%dimensions(1:2)))
          mqcVariableOut = productMatrixReal
        case(33)
          productMatrixInteger = MatMul(  &
            RESHAPE(mqcVariable1%integerArray(:),mqcVariable1%dimensions(1:2)),  &
            RESHAPE(mqcVariable2%integerArray(:),mqcVariable2%dimensions(1:2)))
          mqcVariableOut = productMatrixInteger
        case default
          call mqc_error_I('MQC_Variable_Multiplication: Combination of var1 and var 2 types is UNKNOWN.', 6, &
              'typeCode1*10', typeCode1*10, &
              'typeCode2', typeCode2 )
        end select
      elseIf(RANK(mqcVariable1).eq.2.and.RANK(mqcVariable2).eq.1) then
        select case(typeCode1*10 + typeCode2)
        case(22)
          productVectorReal = MatMul(  &
            RESHAPE(mqcVariable1%realArray(:),mqcVariable1%dimensions(1:2)),  &
            RESHAPE(mqcVariable2%realArray(:),mqcVariable2%dimensions(1:1)))
          mqcVariableOut = productVectorReal
        case(23)
          productVectorReal = MatMul(  &
            RESHAPE(mqcVariable1%realArray(:),mqcVariable1%dimensions(1:2)),  &
            RESHAPE(mqcVariable2%integerArray(:),mqcVariable2%dimensions(1:1)))
          mqcVariableOut = productVectorReal
        case(32)
          productVectorReal = MatMul(  &
            RESHAPE(mqcVariable1%integerArray(:),mqcVariable1%dimensions(1:2)),  &
            RESHAPE(mqcVariable2%realArray(:),mqcVariable2%dimensions(1:1)))
          mqcVariableOut = productVectorReal
        case(33)
          productVectorInteger = MatMul(  &
            RESHAPE(mqcVariable1%integerArray(:),mqcVariable1%dimensions(1:2)),  &
            RESHAPE(mqcVariable2%integerArray(:),mqcVariable2%dimensions(1:1)))
          mqcVariableOut = productVectorInteger
        case default
          call mqc_error_I('MQC_Variable_Multiplication: Combination of var1 and var 2 types is UNKNOWN.', 6, &
              'typeCode1*10', typeCode1*10, &
              'typeCode2', typeCode2 )
        end select
      elseIf(RANK(mqcVariable1).eq.1.and.RANK(mqcVariable2).eq.2) then
        call mqc_error('MQC_Variable_Matrix_Multiplication: vector-matrix multiplication NYI.')
      else
        call mqc_error('MQC_Variable_Matrix_Multiplication: Unknow RANK combination.')
      endIf
!
      return
      end function MQC_Variable_Matrix_Multiplication


!
!PROCEDURE MQC_Variable_DiagonalMatrix2FullMatrix
      function MQC_Variable_DiagonalMatrix2FullMatrix(mqcVariable) result(mqcVariableOut)
!
!     This function converts a vector stored in an MQC_Variable argument
!     <mqcVariable> into a diagonal matrix using full storage.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable
      type(MQC_Variable)::mqcVariableOut
      real(kind=real64),dimension(:,:),allocatable::matrixTmp
!
!
!     Start by ensuring <mqcVariable> is a real vector.
!
      if(.not.RANK(mqcVariable).eq.1)  &
        call mqc_error_i(  &
        'MQC_Variable_DiagonalMatrix2FullMatrix: mqcVariable must be a vector.',  &
        6,'Rank',RANK(mqcVariable))
      if(MQC_Variable_getTypeCode(mqcVariable).ne.2)  &
        call mqc_error_i(  &
        'MQC_Variable_DiagonalMatrix2FullMatrix: mqcVariable must be a real vector.',  &
        6,'TypeCode',MQC_Variable_getTypeCode(mqcVariable))
!
!     Convert the vector into a matrix.
!
      Allocate(matrixTmp(SIZE(mqcVariable),SIZE(mqcVariable)))
      call mqc_packedDiagonalMatrix2FullMatrix_real(  &
        mqcVariable%realArray(:),matrixTmp)
      mqcVariableOut = matrixTmp
      Deallocate(matrixTmp)
!
      return
      end function MQC_Variable_DiagonalMatrix2FullMatrix


!
!PROCEDURE MQC_Variable_UnitMatrix
      function MQC_Variable_UnitMatrix(N,val) result(mqcVariableOut)
!
!     This function makes an (<N> x <N>) unit matrix.
!
!     If the optional dummy argument <val> is sent, its value is put on the
!     diagonal elements rather than one.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      integer(kind=int64),intent(in)::N
      real(kind=real64),optional,intent(in)::val
      type(MQC_Variable)::mqcVariableTmp,mqcVariableOut
      real(kind=real64),dimension(:),allocatable::vectorTmp
!
!     Do the work.
!
      Allocate(vectorTmp(N))
      if(PRESENT(val)) then
        vectorTmp = val
      else
        vectorTmp = 1
      endIf
      mqcVariableTmp = vectorTmp
      mqcVariableOut = mqcVariableTmp%diag()
      Deallocate(vectorTmp)
!
      return
      end function MQC_Variable_UnitMatrix


!
!PROCEDURE MQC_Variable_Matrix_Eigensystem
      subroutine MQC_Variable_Matrix_Eigensystem(mqcVariable,eigenValues,eigenVectors)
!
!     This routine carries out eigen decomposision of a matrix stored in input
!     argument <mqcVariable>.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable
      type(MQC_Variable),intent(out)::eigenValues,eigenVectors
      real(kind=real64),dimension(:),allocatable::eigenValuesTmp
      real(kind=real64),dimension(:,:),allocatable::eigenVectorsTmp
!
!
!     Start by ensuring <mqcVariable> is a real square matrix.
!
      if(.not.RANK(mqcVariable).eq.2)  &
        call mqc_error_i('Eigensystem: mqcVariable must be a matrix.',  &
        6,'Rank',RANK(mqcVariable))
      if(SIZE(mqcVariable,1).ne.SIZE(mqcVariable,2))  &
        call mqc_error_i('Eigensystem: mqcVariable must be a square matrix.',  &
        6,'M',SIZE(mqcVariable,1),'N',SIZE(mqcVariable,2))
      if(MQC_Variable_getTypeCode(mqcVariable).ne.2)  &
        call mqc_error_i('Eigensystem: mqcVariable must be a real matrix.',  &
        6,'TypeCode',MQC_Variable_getTypeCode(mqcVariable))
!
!     Call the LAPACK wrapper for eigen decomposition.
!
      call mqc_dsyev_eigensystem_symmFull(  &
        RESHAPE(mqcVariable%realArray(:),mqcVariable%dimensions(1:2)),  &
        eigenValuesTmp,eigenVectorsTmp)
      eigenValues  = eigenValuesTmp
      eigenVectors = eigenVectorsTmp
      Deallocate(eigenValuesTmp,eigenVectorsTmp)
!
      return
      end subroutine MQC_Variable_Matrix_Eigensystem


!
!PROCEDURE MQC_Variable_Matrix_Determinant
      function MQC_Variable_Matrix_Determinant(mqcVariable) result(det)
!
!     This routine determines the determinant of a matrix of a matrix stored in
!     input argument <mqcVariable>.
!
!
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
      class(MQC_Variable),intent(in)::mqcVariable
      type(MQC_Variable)::det
      type(MQC_Variable)::eigenValues,eigenVectors
      integer(kind=int64)::N,i
      real(kind=real64)::productReal,productImag,productReal1,  &
        productImag1
      real(kind=real64),dimension(:),allocatable::eigenValuesTmp,  &
        eigenValuesImagTmp
      real(kind=real64),dimension(:,:),allocatable::eigenVectorsLeftTmp,  &
        eigenVectorsRightTmp
!
!
!     Start by ensuring <mqcVariable> is a real square matrix.
!
      if(.not.RANK(mqcVariable).eq.2)  &
        call mqc_error_i('Determinant: mqcVariable must be a matrix.',  &
        6,'Rank',RANK(mqcVariable))
      if(SIZE(mqcVariable,1).ne.SIZE(mqcVariable,2))  &
        call mqc_error_i('Determinant: mqcVariable must be a square matrix.',  &
        6,'M',SIZE(mqcVariable,1),'N',SIZE(mqcVariable,2))
      if(MQC_Variable_getTypeCode(mqcVariable).ne.2)  &
        call mqc_error_i('Determinant: mqcVariable must be a real matrix.',  &
        6,'TypeCode',MQC_Variable_getTypeCode(mqcVariable))
!
!     Find the eigenvalues of the matrix.
!
      N = SIZE(mqcVariable,1)
      Allocate(eigenValuesTmp(N),eigenValuesImagTmp(N))
      Allocate(eigenVectorsLeftTmp(N,N),eigenVectorsRightTmp(N,N))
      call mqc_dgeev_eigensystem_Full(  &
        RESHAPE(mqcVariable%realArray(:),mqcVariable%dimensions(1:2)),  &
        eigenValuesTmp,eigenValuesImagTmp,eigenVectorsLeftTmp,  &
        eigenVectorsRightTmp)
!
!     Solve for the determinant now. We check if the final determinant is
!     complex. If it is, we die.
!
      productReal = eigenValuesTmp(1)
      productImag = eigenValuesImagTmp(1)
      do i = 2,N
        productReal1 = productReal*eigenValuesTmp(i)  &
          - productImag*eigenValuesImagTmp(i)
        productImag1 = productImag*eigenValuesTmp(i)  &
          + productReal*eigenValuesImagTmp(i)
        productReal = productReal1
        productImag = productImag1
      endDo
      if(ABS(productImag).gt.MQC_small)  &
        call mqc_error('MQC_Variable_Matrix_Determinant: Complex determinant. NYI.')
      det = productReal
      return
      end function MQC_Variable_Matrix_Determinant
!
!
      end module MQC_Algebra2
