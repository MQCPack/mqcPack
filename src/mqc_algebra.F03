      Module MQC_Algebra
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
      Use MQC_General
      use iso_fortran_env
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     Scalars...
      Type MQC_Scalar
        Real(kind=real64),Private,Allocatable::ScaR
        Integer(kind=int64),Private,Allocatable::ScaI
        Complex(Kind=real64),Private,Allocatable::ScaC
        Character(Len=64),Private::Data_Type ! Real,Integer,Complex
      Contains
        Procedure, Public::print => MQC_Print_Scalar_Algebra1
        Procedure, Public::rval => MQC_Scalar_Get_Intrinsic_Real
        Procedure, Public::ival => MQC_Scalar_Get_Intrinsic_Integer
        Procedure, Public::cval => MQC_Scalar_Get_Intrinsic_Complex
        Procedure, Public::absval => MQC_Scalar_Get_ABS_Value
      End Type MQC_Scalar
!
!     Vectors...
      Type MQC_Vector
        Integer(kind=int64)::Length=0
        Logical,Private::Column=.True.
        Character(Len=64)::Data_Type ! Real,Integer,Complex
        Real(kind=real64),Dimension(:),Allocatable::VecR
        Integer(kind=int64),Dimension(:),Allocatable::VecI
        Complex(Kind=real64),Dimension(:),Allocatable::VecC
      Contains
        Procedure, Public::print => MQC_Print_Vector_Algebra1
        Procedure, Public::initialize => MQC_Vector_Initialize
        Procedure, Public::size => MQC_Length_Vector
        Procedure, Public::init => MQC_Vector_Initialize
        Procedure, Public::norm => MQC_Vector_Norm
        Procedure, Public::transpose => MQC_Vector_Transpose
        Procedure, Public::dagger => MQC_Vector_Conjugate_Transpose
        Procedure, Public::at => MQC_Vector_Scalar_At
        Procedure, Public::vat => MQC_Vector_Vector_At
        Procedure, Public::put => MQC_Vector_Scalar_Put
        Procedure, Public::vput => MQC_Vector_Vector_Put
        Procedure, Public::push => MQC_Vector_Push
        Procedure, Public::unshift => MQC_Vector_Unshift
        Procedure, Public::pop => MQC_Vector_Pop
        Procedure, Public::shift => MQC_Vector_Shift
        Procedure, Public::maxval => MQC_Vector_MaxVal
        Procedure, Public::minval => MQC_Vector_MinVal
        Procedure, Public::sqrt => MQC_Vector_Sqrt
        Procedure, Public::power => MQC_Vector_Power
        Procedure, Public::diag => MQC_Matrix_DiagMatrix_Put_Vector
      End Type MQC_Vector
!
!     Matrices...
!     Lower triangular matrix is stored in a row-major manner so that it is stored one row at a time.
      Type MQC_Matrix
        Integer(kind=int64),Private::NCol=0,NRow=0
        Character(Len=64),Private::Data_Type ! Real,Integer,Complex
        Character(Len=64),Private::Storage ! StorFull,StorSymm,StorDiag,StorTriD
        Real(kind=real64),Dimension(:,:),Allocatable::MatR
        Integer(kind=int64),Dimension(:,:),Allocatable::MatI
        Complex(Kind=real64),Dimension(:,:),Allocatable::MatC
      Contains
        Procedure, Public::print => MQC_Print_Matrix_Algebra1
        Procedure, Public::initialize => MQC_Matrix_Initialize
        Procedure, Public::init => MQC_Matrix_Initialize
        Procedure, Public::identity => MQC_Matrix_Identity
        Procedure, Public::set => MQC_Matrix_Set
        Procedure, Public::norm => MQC_Matrix_Norm
        Procedure, Public::transpose => MQC_Matrix_Transpose
        Procedure, Public::dagger => MQC_Matrix_Conjugate_Transpose
        Procedure, Public::diag => MQC_Matrix_Diagonalize
        Procedure, Public::svd => MQC_Matrix_SVD
        Procedure, Public::eigensys => mqc_matrix_generalized_eigensystem
        Procedure, Public::inv => MQC_Matrix_Inverse
        Procedure, Public::det => mqc_matrix_determinant
        Procedure, Public::trace => mqc_matrix_trace
        Procedure, Public::rmsmax => MQC_Matrix_RMS_Max
        Procedure, Public::sqrt => MQC_Matrix_Sqrt
        Procedure, Public::at => MQC_Matrix_Scalar_At
        Procedure, Public::vat => MQC_Matrix_Vector_At
        Procedure, Public::mat => MQC_Matrix_Matrix_At
        Procedure, Public::put => MQC_Matrix_Scalar_Put
        Procedure, Public::vput => MQC_Matrix_Vector_Put
        Procedure, Public::mput => MQC_Matrix_Matrix_Put
      End Type MQC_Matrix

      Type MQC_R4Tensor
        Integer(kind=int64),Private::I,J,K,L
        Character(Len=64),Private::Data_Type ! Real,Integer,Complex
        Character(Len=64),Private::Storage ! StorFull
        Real(kind=real64),Allocatable,Private::RTen(:,:,:,:)
        Integer(kind=int64),Allocatable,Private::ITen(:,:,:,:)
        Complex(Kind=real64),Allocatable,Private::CTen(:,:,:,:)
      Contains
        Procedure, Public::print => MQC_Print_r4tensor_Algebra1
        Procedure, Public::at => MQC_R4Tensor_At
        Procedure, Public::put => MQC_R4Tensor_Put
        Procedure, Public::initialize => MQC_R4Tensor_Initialize
      End Type MQC_R4Tensor
!
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!
      Interface MQC_Print
        Module Procedure MQC_Print_Scalar_Algebra1
        Module Procedure MQC_Print_Vector_Algebra1
        Module Procedure MQC_Print_Matrix_Algebra1
        Module Procedure MQC_Print_r4tensor_Algebra1
      End Interface
!
      Interface Contraction
        Module Procedure MQC_Matrix_Matrix_Contraction
      End Interface

      interface conjg
        module procedure MQC_Scalar_Complex_Conjugate
      end interface

      interface mqc_have_real
        Module procedure mqc_vector_havereal
        module procedure mqc_matrix_havereal
      end interface

      interface mqc_have_int
        module procedure mqc_vector_haveinteger
        module procedure mqc_matrix_haveinteger
      end interface

      interface mqc_have_complex
        module procedure mqc_vector_havecomplex
        module procedure mqc_matrix_havecomplex
      end interface

      interface mqc_cast_real
        module procedure mqc_vector_cast_real
        module procedure mqc_matrix_cast_real
      end interface

      interface mqc_cast_complex
        module procedure mqc_vector_cast_complex
        module procedure mqc_matrix_cast_complex
      end interface
!
      Interface Matmul
        Module Procedure MQC_MatrixMatrixDotProduct
        Module Procedure MQC_MatrixVectorDotProduct
        Module Procedure MQC_VectorMatrixDotProduct
      End Interface
!
      Interface Transpose
        Module Procedure MQC_Vector_Transpose
        Module Procedure MQC_Matrix_Transpose
      End Interface
!
      Interface Dagger
        Module Procedure MQC_Vector_Conjugate_Transpose
        Module Procedure MQC_Matrix_Conjugate_Transpose
      End Interface
!      
      Interface Sqrt
        Module Procedure MQC_Scalar_Sqrt
      End Interface
!
      Interface MQC_Set_Array2Vector
        Module Procedure MQC_Set_Array2Vector_Integer
        Module Procedure MQC_Set_Array2Vector_Real
        Module Procedure MQC_Set_Array2Vector_Complex
      End Interface
!
      Interface MQC_Matrix_SymmMatrix_Put
        Module Procedure MQC_Matrix_SymmMatrix_Put_Integer,MQC_Matrix_SymmMatrix_Put_Real, MQC_Matrix_SymmMatrix_Put_Complex
      End Interface
!
      Interface MQC_Matrix_DiagMatrix_Put
        Module Procedure MQC_Matrix_DiagMatrix_Put_Integer,MQC_Matrix_DiagMatrix_Put_Real,MQC_Matrix_DiagMatrix_Put_Complex
        Module Procedure MQC_Matrix_DiagMatrix_Put_Vector
      End Interface
!
      Interface Matrix_Symm2Sq
        Module Procedure Matrix_Symm2Sq_Integer
        Module Procedure Matrix_Symm2Sq_Real
        Module Procedure Matrix_Symm2Sq_Complex
      End Interface
!
      Interface Dot_Product
        Module Procedure MQC_VectorVectorDotProduct
      End Interface
!
!----------------------------------------------------------------
!                                                               |
!     OPERATOR INTERFACES                                       |
!                                                               |
!----------------------------------------------------------------
!
!
!     Define Operators.
!
!     Scalar operators...
      Interface Assignment (=)
        Module Procedure MQC_Input_Integer_Scalar
        Module Procedure MQC_Input_Real_Scalar
        Module Procedure MQC_Input_Complex_Scalar
        Module Procedure MQC_Output_MQCScalar_Scalar
        Module Procedure MQC_Output_Integer_Scalar
        Module Procedure MQC_Output_Real_Scalar
        Module Procedure MQC_Output_Complex_Scalar
      End Interface
      Interface Operator (+)
        Module Procedure MQC_ScalarAdd
      End Interface
      Interface Operator (-)
        Module Procedure MQC_ScalarSubtract
      End Interface
      Interface Operator (*)
        Module Procedure MQC_ScalarMultiply
        Module Procedure MQC_ScalarVectorProduct
        Module Procedure MQC_VectorScalarProduct
        Module Procedure MQC_ScalarMatrixProduct
        Module Procedure MQC_MatrixScalarProduct
      End Interface
      Interface Operator (/)
        Module Procedure MQC_ScalarDivide
      End Interface
      Interface Operator (**)
        Module Procedure MQC_ScalarExponent
      End Interface
      Interface Operator (.ne.)
        Module Procedure MQC_ScalarNE
      End Interface
      Interface Operator (.eq.)
        Module Procedure MQC_ScalarEQ
      End Interface
      Interface Operator (.lt.)
        Module Procedure MQC_ScalarLT
      End Interface
      Interface Operator (.gt.)
        Module Procedure MQC_ScalarGT
      End Interface
      Interface Operator (.le.)
        Module Procedure MQC_ScalarLE
      End Interface
      Interface Operator (.ge.)
        Module Procedure MQC_ScalarGE
      End Interface
!
!     Vector operators...
      Interface Assignment (=)
        Module Procedure MQC_Set_Vector2Vector
        Module Procedure MQC_Set_Vector2IntegerArray, MQC_Set_Vector2RealArray, MQC_Set_Vector2ComplexArray
        Module Procedure MQC_Set_Array2Vector_Integer, MQC_Set_Array2Vector_Real, MQC_Set_Array2Vector_Complex 
      End Interface
      Interface Operator (.dot.)
        Module Procedure MQC_VectorVectorDotProduct
      End Interface
      Interface Operator (+)
        Module Procedure MQC_VectorVectorSum
        Module Procedure MQC_ScalarVectorSum
      End Interface
      Interface Operator (-)
        Module Procedure MQC_VectorVectorDifference
        Module Procedure MQC_ScalarVectorDifference
      End Interface
      Interface Operator (.ewp.)
        Module Procedure MQC_ElementVectorProduct
        Module Procedure MQC_ElementMatrixProduct
      End Interface
      Interface Operator (.ewd.)
        Module Procedure MQC_ElementMatrixDivide
      End Interface
      Interface Operator (.x.)
        Module Procedure MQC_CrossProduct
      End Interface
!
!     Matrix operators...
      Interface Assignment (=)
        Module Procedure MQC_Set_Matrix2Matrix
        Module Procedure MQC_Set_Matrix2IntegerArray
        Module Procedure MQC_Set_Matrix2RealArray
        Module Procedure MQC_Set_Matrix2ComplexArray
        Module Procedure MQC_Set_IntegerArray2Matrix
        Module Procedure MQC_Set_RealArray2Matrix
        Module Procedure MQC_Set_ComplexArray2Matrix
      End Interface
      Interface Operator (+)
        Module Procedure MQC_MatrixMatrixSum
      End Interface
      Interface Operator (-)
        Module Procedure MQC_MatrixMatrixSubtract
      End Interface
      Interface Operator (*)
        Module Procedure MQC_MatrixMatrixProduct
      End Interface
      Interface Operator (.dot.)
        Module Procedure MQC_VectorMatrixDotProduct
        Module Procedure MQC_MatrixVectorDotProduct
        Module Procedure MQC_MatrixMatrixDotProduct
      End Interface

!     Rank 4 tensor operators...
      Interface Assignment (=)
        Module Procedure MQC_Set_Array2tensor
      End Interface
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
!
!     PROCEDURE factorial
      function factorial(n)
!
!     This function computes the factorial of a number
!
!     Variable Declarations...
!
      Implicit None
      Integer(kind=int64),Intent(In)::N
      Integer(kind=int64)::I,Ans,Factorial
!
      ans = 1
      Do I = 1, N
        Ans = Ans * I
      EndDo
!
      Factorial = Ans
!
      End Function Factorial
!
!
!     PROCEDURE Bin_Coeff
      Function Bin_Coeff(N,K)
!
!     This function computes the binomial coefficient
!
!     Variable Declarations...
!
      Implicit None
      Integer(kind=int64),Intent(In)::N,K
      Integer(kind=int64)::I,J,Ans,Bin_Coeff
!
      If(K.lt.0.or.K.gt.N) then
        Ans = 0
      ElseIf(K.eq.0.or.K.eq.N) then
        Ans = 1
      Else
        J = Min(K, N-K)
        Ans = 1
        Do I = 0, J-1
          Ans = Ans * (N - I) / (I + 1)
        EndDo
      EndIf
!
      Bin_Coeff = Ans
!
      End Function Bin_Coeff
!
!----------------------------------------------------------------
!                                                               |
!     SCALAR PROCEDURES                                         |
!                                                               |
!----------------------------------------------------------------
!
!     PROCEDURE MQC_Allocate_Scalar
      Subroutine MQC_Allocate_Scalar(Scalar,Data_type)
!
!     This subroutine is used to allocate a scalar type variable
!     of the MQC_Scalar class.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Character(Len=*),Intent(In)::Data_Type
      Type(MQC_Scalar),Intent(InOut)::Scalar
!
      Call MQC_Deallocate_Scalar(Scalar)
      If(Data_Type.eq.'Real' .or. Data_Type.eq.'real') then
        Scalar%Data_type = 'Real'
        Allocate(Scalar%ScaR)
      ElseIf(Data_Type.eq.'Integer' .or. Data_Type.eq.'integer' .or. Data_Type.eq.'Int' .or. Data_Type.eq.'int') then
        Scalar%Data_type = 'Integer'
        Allocate(Scalar%ScaI)
      ElseIf(Data_Type.eq.'Complex' .or. Data_Type.eq.'complex') then
        Scalar%Data_type = 'Complex'
        Allocate(Scalar%ScaC)
      Else
        Call MQC_Error_A('ScalarIn type unspecified in MQC_Allocate_Scalar', 6, &
             'Data_Type', Data_Type )
      EndIf
!
      Return
      End Subroutine MQC_Allocate_Scalar
!
!
!     PROCEDURE MQC_Deallocate_Scalar
      Subroutine MQC_Deallocate_Scalar(Scalar)
!
!     This subroutine is used to deallocate a scalar type variable
!     of the MQC_Scalar class.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar),Intent(InOut)::Scalar
!
      If(Allocated(Scalar%ScaI)) Deallocate(Scalar%ScaI)
      If(Allocated(Scalar%ScaR)) Deallocate(Scalar%ScaR)
      If(Allocated(Scalar%ScaC)) Deallocate(Scalar%ScaC)
      If(.not.Allocated(Scalar%ScaR).and..not.Allocated(Scalar%ScaI).and..not.Allocated(Scalar%ScaC)) then
        Scalar%Data_Type = ''
      EndIf
!
      Return
      End Subroutine MQC_Deallocate_Scalar
!
!
!     PROCEDURE MQC_Scalar_IsAllocated
      Function MQC_Scalar_IsAllocated(Scalar)
!
!     This function is used to determine if an MQC_Scalar is Allocated.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar),Intent(InOut)::Scalar
      Logical::MQC_Scalar_IsAllocated
!
      MQC_Scalar_IsAllocated = .False.
      If(Allocated(Scalar%ScaI)) MQC_Scalar_IsAllocated = .True.
      If(Allocated(Scalar%ScaR)) MQC_Scalar_IsAllocated = .True.
      If(Allocated(Scalar%ScaC)) MQC_Scalar_IsAllocated = .True.
!
      Return
      End Function MQC_Scalar_IsAllocated
!
!
!     PROCEDURE MQC_Input_Integer_Scalar
      Subroutine MQC_Input_Integer_Scalar(ScalarOut,ScalarIn)
!
!     This subroutine is used to set an intrinsic scalar to
!     an MQC_Scalar.
!
!     L. M. Thompson, 2016.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar),Intent(InOut)::ScalarOut
      Integer(kind=int64),Intent(In)::ScalarIn
!
      Call MQC_Allocate_Scalar(ScalarOut,'Integer')
      ScalarOut%scai = ScalarIn
!
      Return
      End Subroutine MQC_Input_Integer_Scalar
!
!
!     PROCEDURE MQC_Input_Real_Scalar
      Subroutine MQC_Input_Real_Scalar(ScalarOut,ScalarIn)
!
!     This subroutine is used to set an intrinsic scalar to
!     an MQC_Scalar.
!
!     L. M. Thompson, 2016.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar),Intent(InOut)::ScalarOut
      Real(kind=real64),Intent(In)::ScalarIn
!
      Call MQC_Allocate_Scalar(ScalarOut,'Real')
      ScalarOut%scar = ScalarIn
!
      Return
      End Subroutine MQC_Input_Real_Scalar
!
!
!     PROCEDURE MQC_Input_Complex_Scalar
      Subroutine MQC_Input_Complex_Scalar(ScalarOut,ScalarIn)
!
!     This subroutine is used to set an intrinsic scalar to
!     an MQC_Scalar.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar),Intent(InOut)::ScalarOut
      Complex(Kind=real64),Intent(In)::ScalarIn
!
      Call MQC_Allocate_Scalar(ScalarOut,'Complex')
      ScalarOut%scac = ScalarIn
!
      Return
      End Subroutine MQC_Input_Complex_Scalar
!
!
!     PROCEDURE MQC_Output_MQCScalar_Scalar
      Subroutine MQC_Output_MQCScalar_Scalar(ScalarOut,ScalarIn)
!
!     This routine is used to output a scalar equal to an MQC scalar type.
!
!     L. M. Thompson, 2016.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar),Intent(InOut)::ScalarOut
      Type(MQC_Scalar),Intent(In)::ScalarIn
!
      Call MQC_Deallocate_Scalar(ScalarOut)
      If(ScalarIn%Data_type.eq.'Integer') then
        Allocate(ScalarOut%ScaI)
        ScalarOut%ScaI = ScalarIn%ScaI
        ScalarOut%Data_type = 'Integer'
      ElseIf(ScalarIn%Data_type.eq.'Real') then
        Allocate(ScalarOut%ScaR)
        ScalarOut%ScaR = ScalarIn%ScaR
        ScalarOut%Data_type = 'Real'
      ElseIf(ScalarIn%Data_type.eq.'Complex') then
        Allocate(ScalarOut%ScaC)
        ScalarOut%ScaC = ScalarIn%ScaC
        ScalarOut%Data_type = 'Complex'
      Else
        Call MQC_Error_A('ScalarIn type not assigned in MQC_Scalar2Scalar', 6, &
             'ScalarIn%Data_type', ScalarIn%Data_type )
      EndIf
!
      Return
      End Subroutine MQC_Output_MQCScalar_Scalar
!
!
!     PROCEDURE MQC_Output_Integer_Scalar
      Subroutine MQC_Output_Integer_Scalar(ScalarOut,ScalarIn)
!
!     This routine is used to output a scalar equal to an MQC scalar type.
!
!     L. M. Thompson, 2016.
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(InOut)::ScalarOut
      Type(MQC_Scalar),Intent(In)::ScalarIn
!
      If(ScalarIn%Data_type.eq.'Integer') then
        ScalarOut = ScalarIn%scai
      ElseIf(ScalarIn%Data_type.eq.'Real') then
        ScalarOut = Int(ScalarIn%scar)
      Else
        Call MQC_Error_A('ScalarIn type unspecified in MQC_Output_Scalar', 6, &
             'ScalarIn%Data_type', ScalarIn%Data_type )
      EndIf
!
      Return
      End Subroutine MQC_Output_Integer_Scalar
!
!
!     PROCEDURE MQC_Output_Real_Scalar
      Subroutine MQC_Output_Real_Scalar(ScalarOut,ScalarIn)
!
!     This routine is used to output a scalar equal to an MQC scalar type.
!
!     L. M. Thompson, 2016.
!
!     Variable Declarations.
      Implicit None
      Real(kind=real64),Intent(InOut)::ScalarOut
      Type(MQC_Scalar),Intent(In)::ScalarIn
!
      If(ScalarIn%Data_type.eq.'Real') then
        ScalarOut = ScalarIn%scar
      ElseIf(ScalarIn%Data_type.eq.'Integer') then
        ScalarOut = Dble(ScalarIn%scai)
      Else
        Call MQC_Error_A('ScalarIn type unspecified in MQC_Output_Scalar', 6, &
             'ScalarIn%Data_type', ScalarIn%Data_type )
      EndIf
!
      Return
      End Subroutine MQC_Output_Real_Scalar
!
!
!     PROCEDURE MQC_Output_Complex_Scalar
      Subroutine MQC_Output_Complex_Scalar(ScalarOut,ScalarIn)
!
!     This routine is used to output a scalar equal to an MQC scalar type.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      Implicit None
      Complex(Kind=real64),Intent(InOut)::ScalarOut
      Type(MQC_Scalar),Intent(In)::ScalarIn
!
      If(ScalarIn%Data_type.eq.'Real') then
        ScalarOut = cmplx(ScalarIn%scar,0.0)
      ElseIf(ScalarIn%Data_type.eq.'Integer') then
        ScalarOut = cmplx(ScalarIn%scai,0)
      ElseIf(ScalarIn%Data_type.eq.'Complex') then
        ScalarOut = ScalarIn%scac
      Else
        Call MQC_Error_A('ScalarIn type unspecified in MQC_Output_Scalar', 6, &
             'ScalarIn%Data_type', ScalarIn%Data_type )
      EndIf
!
      Return
      End Subroutine MQC_Output_Complex_Scalar
!
!
!     PROCEDURE MQC_Print_Scalar_Algebra1
      Subroutine MQC_Print_Scalar_Algebra1(Scalar,IOut,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a MQC_Scalar type variable.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Class(MQC_Scalar),Intent(In)::Scalar
      Character(Len=*),Intent(In)::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
!
 1001 Format(1x,A,1x,'=',1x,I14)
 1002 Format(1x,A,1x,'=',1x,F14.6)
 1003 Format(1x,A,1x,'=',1x,F14.6,SP,F14.6,"i")
 1020 Format( " " )
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      If(Scalar%Data_type.eq.'Integer') then
        Write(IOut,1001) TRIM(Header), Scalar%ScaI
      ElseIf(Scalar%Data_type.eq.'Real') then
        Write(IOut,1002) TRIM(Header), Scalar%ScaR
      ElseIf(Scalar%Data_type.eq.'Complex') then
        Write(IOut,1003) TRIM(Header), Scalar%ScaC
      Else
        Call MQC_Error_A('Scalar type unspecified in MQC_Print_Scalar_Algebra1', 6, &
             'Scalar%Data_type', Scalar%Data_type )
      EndIf
      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Scalar_Algebra1
!
!
!     PROCEDURE MQC_Scalar_Sqrt
      Function MQC_Scalar_Sqrt(Scalar)
!
!     This subroutine is used to obtain the square-root of
!     an MQC_Scalar type variable.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar)::MQC_Scalar_Sqrt
      Type(MQC_Scalar),Intent(In)::Scalar
!
      If(Scalar%Data_type.eq.'Real') then
        MQC_Scalar_Sqrt%ScaR = Sqrt(Scalar%ScaR)
        MQC_Scalar_Sqrt%Data_type = 'Real'
      ElseIf(Scalar%Data_Type.eq.'Integer') then
#ifdef PGI
        MQC_Scalar_Sqrt%ScaR = Sqrt(dFloat(Scalar%ScaI))
#else
        MQC_Scalar_Sqrt%ScaR = Sqrt(Float(Scalar%ScaI))
#endif
        MQC_Scalar_Sqrt%Data_type = 'Real'
      ElseIf(Scalar%Data_Type.eq.'Integer') then
        MQC_Scalar_Sqrt%ScaC = Sqrt(Scalar%ScaC)
        MQC_Scalar_Sqrt%Data_type = 'Complex'
      Else
        Call MQC_Error_A('Scalar type unspecified in MQC_Scalar_Sqrt', 6, &
             'Scalar%Data_Type', Scalar%Data_Type )
      EndIf
!
      Return
      End Function MQC_Scalar_Sqrt
!
!
!     PROCEDURE MQC_Scalar_HaveReal
      Function MQC_Scalar_HaveReal(Scalar)
!
!     This function returns TRUE or FALSE indicating whether Scalar has
!     a real value.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Scalar_HaveReal
      Type(MQC_Scalar),Intent(In)::Scalar
!
      MQC_Scalar_HaveReal = Allocated(Scalar%ScaR)
!
      Return
      End Function MQC_Scalar_HaveReal
!
!
!     PROCEDURE MQC_Scalar_HaveInteger
      Function MQC_Scalar_HaveInteger(Scalar)
!
!     This function returns TRUE or FALSE indicating whether Scalar has
!     a integer value.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Scalar_HaveInteger
      Type(MQC_Scalar),Intent(In)::Scalar
!
      MQC_Scalar_HaveInteger = Allocated(Scalar%ScaI)
!
      Return
      End Function MQC_Scalar_HaveInteger
!
!
!     PROCEDURE MQC_Scalar_HaveComplex
      Function MQC_Scalar_HaveComplex(Scalar)
!
!     This function returns TRUE or FALSE indicating whether Scalar has
!     a complex value.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Scalar_HaveComplex
      Type(MQC_Scalar),Intent(In)::Scalar
!
      MQC_Scalar_HaveComplex = Allocated(Scalar%ScaC)
!
      Return
      End Function MQC_Scalar_HaveComplex
!
!
!     PROCEDURE MQC_Scalar_Get_Intrinsic_Real
      Function MQC_Scalar_Get_Intrinsic_Real(Scalar) Result(Output)
!
!     This function returns the scalar value as a fortran real.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Scalar),Intent(In)::Scalar
      Real(kind=real64)::Output
!
      If(Scalar%Data_Type.eq.'Integer') then
        Output = Real(Scalar%ScaI)
      ElseIf(Scalar%Data_Type.eq.'Real') then
        Output = Scalar%ScaR
      ElseIf(Scalar%Data_Type.eq.'Complex') then
        Output = Real(Scalar%ScaC)
      Else
        call mqc_error_A('Data type unrecognised.', 6, &
             'Scalar%Data_Type', Scalar%Data_Type )
      EndIf
!
      Return
      End Function MQC_Scalar_Get_Intrinsic_Real
!
!
!     PROCEDURE MQC_Scalar_Get_Intrinsic_Integer
      Function MQC_Scalar_Get_Intrinsic_Integer(Scalar) Result(Output)
!
!     This function returns the scalar value as a fortran integer.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Scalar),Intent(In)::Scalar
      Integer(kind=int64)::Output
!
      If(Scalar%Data_Type.eq.'Integer') then
        Output = Scalar%ScaI
      ElseIf(Scalar%Data_Type.eq.'Real') then
        Output = int(Scalar%ScaR)
      ElseIf(Scalar%Data_Type.eq.'Complex') then
        Output = int(Scalar%ScaC)
      Else
        call mqc_error_A('Data type unrecognised.', 6, &
             'Scalar%Data_Type', Scalar%Data_Type )
      EndIf
!
      Return
      End Function MQC_Scalar_Get_Intrinsic_Integer
!
!
!     PROCEDURE MQC_Scalar_Get_Intrinsic_Complex
      Function MQC_Scalar_Get_Intrinsic_Complex(Scalar) Result(Output)
!
!     This function returns the scalar value as a fortran complex.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Scalar),Intent(In)::Scalar
      Complex(Kind=real64)::Output
!
      If(Scalar%Data_Type.eq.'Integer') then
        Output = cmplx(Scalar%ScaI,0)
      ElseIf(Scalar%Data_Type.eq.'Real') then
        Output = cmplx(Scalar%ScaR,0.0)
      ElseIf(Scalar%Data_Type.eq.'Complex') then
        Output = Scalar%ScaC
      Else
        call mqc_error_A('Data type unrecognised.', 6, &
             'Scalar%Data_Type', Scalar%Data_Type )
      EndIf
!
      Return
      End Function MQC_Scalar_Get_Intrinsic_Complex
!
!
!     PROCEDURE MQC_Scalar_Get_ABS_Value
      Function MQC_Scalar_Get_ABS_Value(Scalar) Result(Output)
!
!     This function returns the abs value of an MQC Scalar as
!     an intrinsic real
!
!     A. Mahler, 2018
!
!     Variable Declarations
      Implicit None
      Class(MQC_Scalar),Intent(In)::Scalar
      Real(kind=real64)::Output

      if(Scalar%Data_Type.eq.'Integer') then
        output = real(abs(Scalar%ScaI))
      elseif(Scalar%Data_Type.eq.'Real') then
        output = abs(Scalar%ScaR)
      elseif(Scalar%Data_Type.eq.'Complex') then
        output = abs(Scalar%ScaC)
      else
        call mqc_error_A('Data type unrecognized.', 6, &
          'Scalar%Data_Type', Scalar%Data_Type)
      endif
!
      return
      END FUNCTION MQC_Scalar_Get_ABS_Value

!
!
!     PROCEDURE MQC_ScalarAdd
      Function MQC_ScalarAdd(Scalar1,Scalar2)
!
!     This subroutine is used to add MQC_Scalar type variables.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar)::MQC_ScalarAdd
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
 1050 Format( 2A )
      Call MQC_Deallocate_Scalar(MQC_ScalarAdd)
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Real')
        MQC_ScalarAdd%ScaR = Scalar1%ScaR + Scalar2%ScaR
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Real')
#ifdef PGI
        MQC_ScalarAdd%ScaR = Scalar1%ScaR + dFloat(Scalar2%ScaI)
#else
        MQC_ScalarAdd%ScaR = Scalar1%ScaR + Float(Scalar2%ScaI)
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Real')
#ifdef PGI
        MQC_ScalarAdd%ScaR = dFloat(Scalar1%ScaI) + Scalar2%ScaR
#else
        MQC_ScalarAdd%ScaR = Float(Scalar1%ScaI) + Scalar2%ScaR
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Integer')
        MQC_ScalarAdd%ScaI = Scalar1%ScaI + Scalar2%ScaI
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Complex')
        MQC_ScalarAdd%ScaC = Scalar1%ScaC + Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Complex')
        MQC_ScalarAdd%ScaC = cmplx(Scalar1%ScaR,0.0) + Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Complex')
        MQC_ScalarAdd%ScaC = Scalar1%ScaC + cmplx(Scalar2%ScaR,0.0)
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Complex')
        MQC_ScalarAdd%ScaC = cmplx(Scalar1%ScaI,0) + Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarAdd,'Complex')
        MQC_ScalarAdd%ScaC = Scalar1%ScaC + cmplx(Scalar2%ScaI,0)
      Else
        write(*,1050)' Scalar type of #1 = ',TRIM(Scalar1%Data_type)
        write(*,1050)' Scalar type of #2 = ',TRIM(Scalar2%Data_type)
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarAdd', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarAdd
!
!
!     PROCEDURE MQC_ScalarSubtract
      Function MQC_ScalarSubtract(Scalar1,Scalar2)
!
!     This subroutine is used to subtract MQC_Scalar type variables.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar)::MQC_ScalarSubtract
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      Call MQC_Deallocate_Scalar(MQC_ScalarSubtract)
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Real')
        MQC_ScalarSubtract%ScaR = Scalar1%ScaR - Scalar2%ScaR
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Real')
#ifdef PGI
        MQC_ScalarSubtract%ScaR = Scalar1%ScaR - dFloat(Scalar2%ScaI)
#else
        MQC_ScalarSubtract%ScaR = Scalar1%ScaR - Float(Scalar2%ScaI)
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Real')
#ifdef PGI
        MQC_ScalarSubtract%ScaR = dFloat(Scalar1%ScaI) - Scalar2%ScaR
#else
        MQC_ScalarSubtract%ScaR = Float(Scalar1%ScaI) - Scalar2%ScaR
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Integer')
        MQC_ScalarSubtract%ScaI = Scalar1%ScaI - Scalar2%ScaI
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Complex')
        MQC_ScalarSubtract%ScaC = Scalar1%ScaC - Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Complex')
        MQC_ScalarSubtract%ScaC = cmplx(Scalar1%ScaR,0.0) - Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Complex')
        MQC_ScalarSubtract%ScaC = Scalar1%ScaC - cmplx(Scalar2%ScaR,0.0)
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Complex')
        MQC_ScalarSubtract%ScaC = cmplx(Scalar1%ScaI,0) - Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarSubtract,'Complex')
        MQC_ScalarSubtract%ScaC = Scalar1%ScaC - cmplx(Scalar2%ScaI,0)
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarSubtract', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarSubtract
!
!
!     PROCEDURE MQC_ScalarMultiply
      Function MQC_ScalarMultiply(Scalar1,Scalar2)
!
!     This subroutine is used to multiply MQC_Scalar type variables.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar)::MQC_ScalarMultiply
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      Call MQC_Deallocate_Scalar(MQC_ScalarMultiply)
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Real')
        MQC_ScalarMultiply%ScaR = Scalar1%ScaR * Scalar2%ScaR
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Real')
#ifdef PGI
        MQC_ScalarMultiply%ScaR = Scalar1%ScaR * dFloat(Scalar2%ScaI)
#else
        MQC_ScalarMultiply%ScaR = Scalar1%ScaR * Float(Scalar2%ScaI)
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Real')
#ifdef PGI
        MQC_ScalarMultiply%ScaR = dFloat(Scalar1%ScaI) * Scalar2%ScaR
#else
        MQC_ScalarMultiply%ScaR = Float(Scalar1%ScaI) * Scalar2%ScaR
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Integer')
        MQC_ScalarMultiply%ScaI = Scalar1%ScaI * Scalar2%ScaI
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Complex')
        MQC_ScalarMultiply%ScaC = Scalar1%ScaC * Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Complex')
        MQC_ScalarMultiply%ScaC = cmplx(Scalar1%ScaR,0.0) * Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Complex')
        MQC_ScalarMultiply%ScaC = Scalar1%ScaC * cmplx(Scalar2%ScaR,0.0)
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Complex')
        MQC_ScalarMultiply%ScaC = cmplx(Scalar1%ScaI,0) * Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarMultiply,'Complex')
        MQC_ScalarMultiply%ScaC = Scalar1%ScaC * cmplx(Scalar2%ScaI,0)
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarMultiply', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarMultiply
!
!
!     PROCEDURE MQC_ScalarDivide
      Function MQC_ScalarDivide(Scalar1,Scalar2)
!
!     This subroutine is used to divide MQC_Scalar type variables.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar)::MQC_ScalarDivide
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      Call MQC_Deallocate_Scalar(MQC_ScalarDivide)
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Real')
        MQC_ScalarDivide%ScaR = Scalar1%ScaR / Scalar2%ScaR
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Real')
#ifdef PGI
        MQC_ScalarDivide%ScaR = Scalar1%ScaR / dFloat(Scalar2%ScaI)
#else
        MQC_ScalarDivide%ScaR = Scalar1%ScaR / Float(Scalar2%ScaI)
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Real')
#ifdef PGI
        MQC_ScalarDivide%ScaR = dFloat(Scalar1%ScaI) / Scalar2%ScaR
#else
        MQC_ScalarDivide%ScaR = Float(Scalar1%ScaI) / Scalar2%ScaR
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Integer')
        MQC_ScalarDivide%ScaI = Scalar1%ScaI / Scalar2%ScaI
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Complex')
        MQC_ScalarDivide%ScaC = Scalar1%ScaC / Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Complex')
        MQC_ScalarDivide%ScaC = cmplx(Scalar1%ScaR,0.0) / Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Complex')
        MQC_ScalarDivide%ScaC = Scalar1%ScaC / cmplx(Scalar2%ScaR,0.0)
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Complex')
        MQC_ScalarDivide%ScaC = cmplx(Scalar1%ScaI,0) / Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarDivide,'Complex')
        MQC_ScalarDivide%ScaC = Scalar1%ScaC + cmplx(Scalar2%ScaI,0)
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarDivide', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarDivide
!
!
!     PROCEDURE MQC_ScalarExponent
      Function MQC_ScalarExponent(Scalar1,Scalar2)
!
!     This subroutine is used to raise MQC_Scalar to a power.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar)::MQC_ScalarExponent
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      Call MQC_Deallocate_Scalar(MQC_ScalarExponent)
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Real')
        MQC_ScalarExponent%ScaR = Scalar1%ScaR ** Scalar2%ScaR
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Real')
#ifdef PGI
        MQC_ScalarExponent%ScaR = Scalar1%ScaR ** dFloat(Scalar2%ScaI)
#else
        MQC_ScalarExponent%ScaR = Scalar1%ScaR ** Float(Scalar2%ScaI)
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Real')
#ifdef PGI
        MQC_ScalarExponent%ScaR = dFloat(Scalar1%ScaI) ** Scalar2%ScaR
#else
        MQC_ScalarExponent%ScaR = Float(Scalar1%ScaI) ** Scalar2%ScaR
#endif
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Integer')
        MQC_ScalarExponent%ScaI = Scalar1%ScaI ** Scalar2%ScaI
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Complex')
        MQC_ScalarExponent%ScaC = Scalar1%ScaC ** Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Complex')
        MQC_ScalarExponent%ScaC = cmplx(Scalar1%ScaR,0.0) ** Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Real') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Complex')
        MQC_ScalarExponent%ScaC = Scalar1%ScaC ** cmplx(Scalar2%ScaR,0.0)
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Complex') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Complex')
        MQC_ScalarExponent%ScaC = cmplx(Scalar1%ScaI,0) ** Scalar2%ScaC
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Integer') then
        Call MQC_Allocate_Scalar(MQC_ScalarExponent,'Complex')
        MQC_ScalarExponent%ScaC = Scalar1%ScaC ** cmplx(Scalar2%ScaI,0)
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarExponent', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarExponent
!
!
!     PROCEDURE MQC_ScalarNE
      Function MQC_ScalarNE(Scalar1,Scalar2)
!
!     This subroutine is used to return if two MQC_Scalar type variables
!     are not equal.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_ScalarNE
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scar.eq.Scalar2%scar) MQC_ScalarNE = .False.
        If(Scalar1%scar.ne.Scalar2%scar) MQC_ScalarNE = .True.
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scar.eq.Scalar2%scai) MQC_ScalarNE = .False.
        If(Scalar1%scar.ne.Scalar2%scai) MQC_ScalarNE = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scai.eq.Scalar2%scar) MQC_ScalarNE = .False.
        If(Scalar1%scai.ne.Scalar2%scar) MQC_ScalarNE = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scai.eq.Scalar2%scai) MQC_ScalarNE = .False.
        If(Scalar1%scai.ne.Scalar2%scai) MQC_ScalarNE = .True.
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Complex') then
        If(Scalar1%scac.eq.Scalar2%scac) MQC_ScalarNE = .False.
        If(Scalar1%scac.ne.Scalar2%scac) MQC_ScalarNE = .True.
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Complex') then
        If(Scalar1%scar.eq.Scalar2%scac) MQC_ScalarNE = .False.
        If(Scalar1%scar.ne.Scalar2%scac) MQC_ScalarNE = .True.
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scac.eq.Scalar2%scar) MQC_ScalarNE = .False.
        If(Scalar1%scac.ne.Scalar2%scar) MQC_ScalarNE = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Complex') then
        If(Scalar1%scai.eq.Scalar2%scac) MQC_ScalarNE = .False.
        If(Scalar1%scai.ne.Scalar2%scac) MQC_ScalarNE = .True.
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scac.eq.Scalar2%scai) MQC_ScalarNE = .False.
        If(Scalar1%scac.ne.Scalar2%scai) MQC_ScalarNE = .True.
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarNE', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarNE
!
!
!     PROCEDURE MQC_ScalarEQ
      Function MQC_ScalarEQ(Scalar1,Scalar2)
!
!     This subroutine is used to return if two MQC_Scalar type variables
!     are equal.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_ScalarEQ
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      MQC_ScalarEQ = .False.
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scar.eq.Scalar2%scar) MQC_ScalarEQ = .True.
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scar.eq.Scalar2%scai) MQC_ScalarEQ = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scai.eq.Scalar2%scar) MQC_ScalarEQ = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scai.eq.Scalar2%scai) MQC_ScalarEQ = .True.
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Complex') then
        If(Scalar1%scac.eq.Scalar2%scac) MQC_ScalarEQ = .True.
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Complex') then
        If(Scalar1%scar.eq.Scalar2%scac) MQC_ScalarEQ = .True.
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scac.eq.Scalar2%scar) MQC_ScalarEQ = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Complex') then
        If(Scalar1%scai.eq.Scalar2%scac) MQC_ScalarEQ = .True.
      ElseIf(Scalar1%Data_type.eq.'Complex'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scac.eq.Scalar2%scai) MQC_ScalarEQ = .True.
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarEQ', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarEQ
!
!
!     PROCEDURE MQC_ScalarLT
      Function MQC_ScalarLT(Scalar1,Scalar2)
!
!     This subroutine is used to return if the left  MQC_Scalar type variable
!     is less than the right MQC_Scalar variable.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_ScalarLT
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      MQC_ScalarLT = .False.
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scar.lt.Scalar2%scar) MQC_ScalarLT = .True.
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scar.lt.Scalar2%scai) MQC_ScalarLT = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scai.lt.Scalar2%scar) MQC_ScalarLT = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scai.lt.Scalar2%scai) MQC_ScalarLT = .True.
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarLT', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarLT
!
!
!     PROCEDURE MQC_ScalarGT
      Function MQC_ScalarGT(Scalar1,Scalar2)
!
!     This subroutine is used to return if the left  MQC_Scalar type variable
!     is greater than the right MQC_Scalar variable.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_ScalarGT
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      MQC_ScalarGT = .False.
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scar.gt.Scalar2%scar) MQC_ScalarGT = .True.
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scar.gt.Scalar2%scai) MQC_ScalarGT = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scai.gt.Scalar2%scar) MQC_ScalarGT = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scai.gt.Scalar2%scai) MQC_ScalarGT = .True.
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarGT', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarGT
!
!
!     PROCEDURE MQC_ScalarLE
      Function MQC_ScalarLE(Scalar1,Scalar2)
!
!     This subroutine is used to return if the left  MQC_Scalar type variable
!     is less than or equal to the right MQC_Scalar variable.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_ScalarLE
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      MQC_ScalarLE = .False.
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scar.le.Scalar2%scar) MQC_ScalarLE = .True.
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scar.le.Scalar2%scai) MQC_ScalarLE = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scai.le.Scalar2%scar) MQC_ScalarLE = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scai.le.Scalar2%scai) MQC_ScalarLE = .True.
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarLE', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarLE
!
!
!     PROCEDURE MQC_ScalarGE
      Function MQC_ScalarGE(Scalar1,Scalar2)
!
!     This subroutine is used to return if the left  MQC_Scalar type variable
!     is greater than or equal to the right MQC_Scalar variable.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_ScalarGE
      Type(MQC_Scalar),Intent(In)::Scalar1,Scalar2
!
      MQC_ScalarGE = .False.
      If(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scar.ge.Scalar2%scar) MQC_ScalarGE = .True.
      ElseIf(Scalar1%Data_type.eq.'Real'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scar.ge.Scalar2%scai) MQC_ScalarGE = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Real') then
        If(Scalar1%scai.ge.Scalar2%scar) MQC_ScalarGE = .True.
      ElseIf(Scalar1%Data_type.eq.'Integer'.and.Scalar2%Data_type.eq.'Integer') then
        If(Scalar1%scai.ge.Scalar2%scai) MQC_ScalarGE = .True.
      Else
        Call MQC_Error_A('Scalar types unspecified in MQC_ScalarGE', 6, &
             'Scalar1%Data_type', Scalar1%Data_type, &
             'Scalar2%Data_type', Scalar2%Data_type )
      EndIf
!
      Return
      End Function MQC_ScalarGE
!
!
!     PROCEDURE MQC_Scalar_Complex_Conjugate
      Function MQC_Scalar_Complex_Conjugate(ScalarIn) Result(ScalarOut)
!
!     This function returns the complex conjugate of an mqc scalar. 
!
!     L. M. Thompson, 2018.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Scalar),Intent(In)::ScalarIn
      Type(MQC_Scalar)::ScalarOut
      Character(Len=64)::datatype
!
      If(ScalarIn%Data_Type.eq.'Complex') then
        ScalarOut = conjg(ScalarIn%ScaC)
      Else
        ScalarOut = ScalarIn
      EndIf
!
      Return
      End Function MQC_Scalar_Complex_Conjugate
!
!----------------------------------------------------------------
!                                                               |
!     VECTOR PROCEDURES                                         |
!                                                               |
!----------------------------------------------------------------
!
!     PROCEDURE MQC_Allocate_Vector
      Subroutine MQC_Allocate_Vector(N,Vector,Data_Type)
!
!     This subroutine is used to allocate memory for a vector variable
!     belonging to the MQC_Vector class.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::N
      Character(Len=*),Intent(In)::Data_Type
      Type(MQC_Vector),Intent(InOut)::Vector
!
      Call MQC_Deallocate_Vector(Vector)
      If(Data_Type.eq.'Integer' .or. Data_Type.eq.'integer' .or. Data_Type.eq.'Int' .or. Data_Type.eq.'int') then
        Vector%Data_Type = 'Integer'
        If(Allocated(Vector%veci)) DeAllocate(Vector%veci)
        Allocate(Vector%veci(N))
      ElseIf(Data_Type.eq.'Real' .or. Data_Type.eq.'real') then
        Vector%Data_Type = 'Real'
        If(Allocated(Vector%vecr)) DeAllocate(Vector%vecr)
        Allocate(Vector%vecr(N))
      ElseIf(Data_Type.eq.'Complex' .or. Data_Type.eq.'complex') then
        Vector%Data_Type = 'Complex'
        If(Allocated(Vector%vecc)) DeAllocate(Vector%vecc)
        Allocate(Vector%vecc(N))
      EndIf
      Vector%Length = N
      Vector%Column = .True.
!
      Return
      End Subroutine MQC_Allocate_Vector
!
!
!     PROCEDURE MQC_Deallocate_Vector
      Subroutine MQC_Deallocate_Vector(Vector)
!
!     This subroutine is used to de-allocate memory for a vector variable
!     belonging to the MQC_Vector class.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector),Intent(InOut)::Vector
!
      If(Allocated(Vector%veci)) DeAllocate(Vector%veci)
      If(Allocated(Vector%vecr)) DeAllocate(Vector%vecr)
      If(Allocated(Vector%vecc)) DeAllocate(Vector%vecc)
      Vector%Data_Type = ''
      Vector%Length = 0
      Vector%Column = .True.
!
      Return
      End Subroutine MQC_Deallocate_Vector
!
!
!     PROCEDURE MQC_Length_Vector
      Function MQC_Length_Vector(Vector) Result(ValueOut)
!
!     This subroutine is used to return the length of an MQC vector. If the
!     vector is NOT allocated, the length is returned as 0.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Vector)::Vector
      Integer(kind=int64)::ValueOut
!
      ValueOut = Vector%Length
!
      Return
      End Function MQC_Length_Vector
!
!
!     PROCEDURE MQC_Vector_HaveReal
      Function MQC_Vector_HaveReal(Vector)
!
!     This function returns TRUE or FALSE indicating whether Vector has an
!     allocated real vector.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Vector_HaveReal
      Type(MQC_Vector),Intent(In)::Vector
!
      MQC_Vector_HaveReal = Allocated(Vector%vecr)
!
      Return
      End Function MQC_Vector_HaveReal
!
!
!     PROCEDURE MQC_Vector_HaveInteger
      Function MQC_Vector_HaveInteger(Vector)
!
!     This function returns TRUE or FALSE indicating whether Vector has an
!     allocated integer vector.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Vector_HaveInteger
      Type(MQC_Vector),Intent(In)::Vector
!
      MQC_Vector_HaveInteger = Allocated(Vector%veci)
!
      Return
      End Function MQC_Vector_HaveInteger
!
!
!     PROCEDURE MQC_Vector_HaveComplex
      Function MQC_Vector_HaveComplex(Vector)
!
!     This function returns TRUE or FALSE indicating whether Vector has an
!     allocated complex vector.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Vector_HaveComplex
      Type(MQC_Vector),Intent(In)::Vector
!
      MQC_Vector_HaveComplex = Allocated(Vector%vecc)
!
      Return
      End Function MQC_Vector_HaveComplex
!
!
!     PROCEDURE MQC_Vector_IsColumn
      Function MQC_Vector_IsColumn(Vector)
!
!     This function returns TRUE or FALSE indicating whether the vector is
!     a column or row respectively.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Vector_IsColumn
      Type(MQC_Vector),Intent(In)::Vector
!
      MQC_Vector_IsColumn = Vector%Column
!
      Return
      End Function MQC_Vector_IsColumn
!
!
!     PROCEDURE MQC_Vector_Copy_Int2Real
      Subroutine MQC_Vector_Copy_Int2Real(Vector)
!
!     This subroutine copies an integer MQC_Vector into it's real vector
!     space.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::Vector
!
      Call MQC_Allocate_Vector(MQC_Length_Vector(Vector),Vector,  &
        'Real')
      Vector%vecr = Vector%veci
      if(allocated(vector%veci)) deallocate(vector%veci)
!
      Return
      End Subroutine MQC_Vector_Copy_Int2Real
!
!
!     PROCEDURE MQC_Vector_Copy_Int2Complex
      Subroutine MQC_Vector_Copy_Int2Complex(Vector)
!
!     This subroutine copies an integer MQC_Vector into it's complex vector
!     space.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::Vector
!
      Call MQC_Allocate_Vector(MQC_Length_Vector(Vector),Vector,  &
        'Complex')
      Vector%vecc = cmplx(Vector%veci,0)
      if(allocated(vector%veci)) deallocate(vector%veci)
!
      Return
      End Subroutine MQC_Vector_Copy_Int2Complex
!
!
!     PROCEDURE MQC_Vector_Copy_Real2Int
      Subroutine MQC_Vector_Copy_Real2Int(Vector)
!
!     This subroutine copies a real MQC_Vector into it's integer vector
!     space.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::Vector
!
      Call MQC_Allocate_Vector(MQC_Length_Vector(Vector),Vector,  &
        'Integer')
      Vector%veci = Vector%vecr
      if(allocated(vector%vecr)) deallocate(vector%vecr)
!
      Return
      End Subroutine MQC_Vector_Copy_Real2Int
!
!
!     PROCEDURE MQC_Vector_Copy_Real2Complex
      Subroutine MQC_Vector_Copy_Real2Complex(Vector)
!
!     This subroutine copies a real MQC_Vector into it's complex vector
!     space.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::Vector
!
      Call MQC_Allocate_Vector(MQC_Length_Vector(Vector),Vector,  &
        'Complex')
      Vector%vecc = cmplx(Vector%vecr,0.0)
      if(allocated(vector%vecr)) deallocate(vector%vecr)
!
      Return
      End Subroutine MQC_Vector_Copy_Real2Complex
!
!
!     PROCEDURE MQC_Vector_Copy_Complex2Int
      Subroutine MQC_Vector_Copy_complex2Int(Vector)
!
!     This subroutine copies a complex MQC_Vector into it's integer vector
!     space.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::Vector
!
      Call MQC_Allocate_Vector(MQC_Length_Vector(Vector),Vector,  &
        'Integer')
      Vector%veci = real(Vector%vecc)
      if(allocated(vector%vecc)) deallocate(vector%vecc)
!
      Return
      End Subroutine MQC_Vector_Copy_Complex2Int
!
!
!     PROCEDURE MQC_Vector_Copy_Complex2Real
      Subroutine MQC_Vector_Copy_Complex2Real(Vector)
!
!     This subroutine copies a complex MQC_Vector into it's real vector
!     space.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::Vector
!
      Call MQC_Allocate_Vector(MQC_Length_Vector(Vector),Vector,  &
        'Real')
      Vector%vecr = real(Vector%vecc)
      if(allocated(vector%vecc)) deallocate(vector%vecc)
!
      Return
      End Subroutine MQC_Vector_Copy_Complex2Real
!
!
!     PROCEDURE MQC_Vector_Scalar_At
      Function MQC_Vector_Scalar_At(Vec,I) Result(Element)
!
!     A function that returns the element at (i) of a mqc_vector type as a mqc_scalar type.
!
!     X. Sheng 2017.02
!
      Implicit None
      Class(MQC_Vector)::Vec
      Type(MQC_Scalar)::Element
      Integer(kind=int64),Intent(In)::I
      Integer(kind=int64)::IndI

      IndI = I
      If(IndI.lt.0) IndI = Vec%Length + IndI + 1
      If (IndI.gt.Vec%Length.or.IndI.eq.0) Call MQC_Error_I('Index out of bounds in MQC_Vector_Scalar_At', 6, &
           'Vec%Length', Vec%Length, &
           'IndI', IndI )
      If (Vec%Data_Type.eq.'Integer') then
        Element = (Vec%VecI(IndI))
      ElseIf (Vec%Data_Type.eq.'Real') then
        Element = (Vec%VecR(IndI))
      ElseIf (Vec%Data_Type.eq.'Complex') then
        Element = (Vec%VecC(IndI))
      Else
        Call MQC_Error_A('Vector type not defined in MQC_Vector_Scalar_At', 6, &
             'Vec%Data_Type', Vec%Data_Type )
      EndIf

      End Function
!
!
!     PROCEDURE MQC_Vector_Vector_At
      Function MQC_Vector_Vector_At(Vec,I,J) Result(Vector)
!
!     A function that returns MQC_Vector Vector between I and J of an
!     MQC_Vector Vec. Negative I and J count from end of index.
!
!     L. M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Vector)::Vec
      Type(MQC_Vector)::Vector
      Integer(kind=int64),Intent(In)::I
      Integer(kind=int64),Optional,Intent(In)::J
      Integer(kind=int64)::Length,IndI,IndJ

      IndI = I
      If(Present(J)) then
        IndJ = J
      Else
        IndJ = Vec%Length
      EndIf
      If (IndI.lt.0) IndI = Vec%Length + IndI + 1
      If (IndJ.lt.0) IndJ = Vec%Length + IndJ + 1
      If (IndI.eq.0) then
        IndI = 1
        If(Present(J)) Call MQC_Error_L('Vector length badly defined in MQC_Vector_Vector_At', 6, &
             'Present(J)', Present(J) )
        IndJ = Vec%Length
      EndIf

      Length = IndJ-IndI+1
      If (Length.le.0.or.Length.gt.Vec%Length) Call MQC_Error_I('Vector length badly &
     &   defined in MQC_Vector_Vector_At', 6, &
     'Length', Length, &
     'Vec%Length', Vec%Length )
      If (IndI.le.0.or.IndI.gt.(Vec%Length-Length+1)) Call MQC_Error_I('Index I out of bounds &
     &   in MQC_Vector_Vector_At', 6, &
     'IndI', IndI, &
     'Vec%Length', Vec%Length, &
     'Length', Length )
      If (IndJ.lt.Length.or.IndJ.gt.Vec%Length) Call MQC_Error_I('Index J out of bounds &
     &   in MQC_Vector_Vector_At', 6, &
     'IndJ', IndJ, &
     'Length', Length, &
     'Vec%Length', Vec%Length )
      If (Vec%Data_Type.eq.'Integer') then
        Call MQC_Allocate_Vector(Length,Vector,'Integer')
        Vector%VecI = (Vec%VecI(IndI:IndJ))
      ElseIf (Vec%Data_Type.eq.'Real') then
        Call MQC_Allocate_Vector(Length,Vector,'Real')
        Vector%VecR = (Vec%VecR(IndI:IndJ))
      ElseIf (Vec%Data_Type.eq.'Complex') then
        Call MQC_Allocate_Vector(Length,Vector,'Complex')
        Vector%VecC = (Vec%VecC(IndI:IndJ))
      Else
        Call MQC_Error_A('Vector type not defined in MQC_Vector_Vector_At', 6, &
             'Vec%Data_Type', Vec%Data_Type )
      EndIf
      Vector%Column = Vec%Column

      End Function
!
!
!     PROCEDURE MQC_Set_Vector2IntegerArray
      Subroutine MQC_Set_Vector2IntegerArray(ArrayOut,VectorIn)
!
!     This routine is used to set a rank-1 array equal to an MQC vector.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Allocatable,Intent(InOut)::ArrayOut(:)
      Type(MQC_Vector),Intent(In)::VectorIn
      Integer(kind=int64)::LenVectorIn
!
      LenVectorIn = MQC_Length_Vector(VectorIn)
      If(.not.Allocated(ArrayOut)) then
        Allocate(ArrayOut(LenVectorIn))
      ElseIf(LenVectorIn.ne.Size(ArrayOut))  then
        Deallocate(ArrayOut)
        Allocate(ArrayOut(LenVectorIn))
      EndIf

      If(LenVectorIn.gt.0) then
        If(VectorIn%Data_Type.eq.'Integer') then
          ArrayOut = VectorIn%veci
        Elseif(VectorIn%Data_Type.eq.'Real') then
          ArrayOut = VectorIn%vecr
        Elseif(VectorIn%Data_Type.eq.'Complex') then
          ArrayOut = Real(VectorIn%vecc)
        Else
          Call MQC_Error_A('VectorIn type not defined in MQC_Vector2IntegerArray', 6, &
               'VectorIn%Data_Type', VectorIn%Data_Type )
        EndIf
      EndIf
!
      Return
      End Subroutine MQC_Set_Vector2IntegerArray
!
!
!     PROCEDURE MQC_Set_Vector2RealArray
      Subroutine MQC_Set_Vector2RealArray(ArrayOut,VectorIn)
!
!     This routine is used to set a rank-1 array equal to an MQC vector.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Real(kind=real64),Allocatable,Intent(InOut)::ArrayOut(:)
      Type(MQC_Vector),Intent(In)::VectorIn
      Integer(kind=int64)::LenVectorIn
!
      LenVectorIn = MQC_Length_Vector(VectorIn)
      If(.not.Allocated(ArrayOut)) then
        Allocate(ArrayOut(LenVectorIn))
      ElseIf(LenVectorIn.ne.Size(ArrayOut))  then
        Deallocate(ArrayOut)
        Allocate(ArrayOut(LenVectorIn))
      EndIf

      If(LenVectorIn.gt.0) then
        If(VectorIn%Data_Type.eq.'Integer') then
          ArrayOut = VectorIn%veci
        Elseif(VectorIn%Data_Type.eq.'Real') then
          ArrayOut = VectorIn%vecr
        Elseif(VectorIn%Data_Type.eq.'Complex') then
          ArrayOut = real(VectorIn%vecc)
        Else
          Call MQC_Error_A('VectorIn type not defined in MQC_Vector2RealArray', 6, &
               'VectorIn%Data_Type', VectorIn%Data_Type )
        EndIf
      EndIf
!
      Return
      End Subroutine MQC_Set_Vector2RealArray
!
!
!     PROCEDURE MQC_Set_Vector2ComplexArray
      Subroutine MQC_Set_Vector2ComplexArray(ArrayOut,VectorIn)
!
!     This routine is used to set a rank-1 array equal to an MQC vector.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Complex(Kind=real64),Allocatable,Intent(InOut)::ArrayOut(:)
      Type(MQC_Vector),Intent(In)::VectorIn
      Integer(kind=int64)::LenVectorIn
!
      LenVectorIn = MQC_Length_Vector(VectorIn)
      If(.not.Allocated(ArrayOut)) then
        Allocate(ArrayOut(LenVectorIn))
      ElseIf(LenVectorIn.ne.Size(ArrayOut))  then
        Deallocate(ArrayOut)
        Allocate(ArrayOut(LenVectorIn))
      EndIf

      If(LenVectorIn.gt.0) then
        If(VectorIn%Data_Type.eq.'Integer') then
          ArrayOut = VectorIn%veci
        Elseif(VectorIn%Data_Type.eq.'Real') then
          ArrayOut = VectorIn%vecr
        Elseif(VectorIn%Data_Type.eq.'Complex') then
          ArrayOut = VectorIn%vecc
        Else
          Call MQC_Error_A('VectorIn type not defined in MQC_Vector2ComplexArray', 6, &
               'VectorIn%Data_Type', VectorIn%Data_Type )
        EndIf
      EndIf
!
      Return
      End Subroutine MQC_Set_Vector2ComplexArray
!
!
!     PROCEDURE MQC_Set_Array2Vector_Integer
      Subroutine MQC_Set_Array2Vector_Integer(VectorOut,ArrayIn)
!
!     This routine is used to set an MQC vector equal to a rank-1 array.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector),Intent(InOut)::VectorOut
      Integer(kind=int64),Dimension(:),Intent(In)::ArrayIn
!
      Call MQC_Deallocate_Vector(VectorOut)
      Call MQC_Allocate_Vector(Size(ArrayIn),VectorOut,'Integer')
      VectorOut%veci = ArrayIn
!
      Return
      End Subroutine MQC_Set_Array2Vector_Integer
!
!
!     PROCEDURE MQC_Set_Array2Vector_Real
      Subroutine MQC_Set_Array2Vector_Real(VectorOut,ArrayIn)
!
!     This routine is used to set an MQC vector equal to a rank-1 array.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector),Intent(InOut)::VectorOut
      Real(kind=real64),Dimension(:),Intent(In)::ArrayIn
!
      Call MQC_Deallocate_Vector(VectorOut)
      Call MQC_Allocate_Vector(Size(ArrayIn),VectorOut,'Real')
      VectorOut%vecr = ArrayIn
!
      Return
      End Subroutine MQC_Set_Array2Vector_Real
!
!
!     PROCEDURE MQC_Set_Array2Vector_Complex
      Subroutine MQC_Set_Array2Vector_Complex(VectorOut,ArrayIn)
!
!     This routine is used to set an MQC vector equal to a rank-1 array.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector),Intent(InOut)::VectorOut
      Complex(Kind=real64),Dimension(:),Intent(In)::ArrayIn
!
      Call MQC_Deallocate_Vector(VectorOut)
      Call MQC_Allocate_Vector(Size(ArrayIn),VectorOut,'Complex')
      VectorOut%vecc = ArrayIn
!
      Return
      End Subroutine MQC_Set_Array2Vector_Complex
!
!
!     PROCEDURE MQC_Set_Vector2Vector
      Subroutine MQC_Set_Vector2Vector(VectorOut,VectorIn)
!
!     This routine is used to set an MQC vector equal to another MQC
!     vector.
!
!     H. P. Hratchian, 2016.
!     L. M. Thompson, 2016.
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Vector),Intent(InOut)::VectorOut
      Class(MQC_Vector),Intent(In)::VectorIn
!
      Call MQC_Deallocate_Vector(VectorOut)
      If(MQC_Vector_HaveReal(VectorIn)) then
        If(MQC_Length_Vector(VectorOut).eq.0 .or. &
          MQC_Length_Vector(VectorOut).ne.MQC_Length_Vector(VectorIn))  &
          Call MQC_Allocate_Vector(MQC_Length_Vector(VectorIn),VectorOut,'Real')
        VectorOut%vecr = VectorIn%vecr
      EndIf
      If(MQC_Vector_HaveInteger(VectorIn)) then
        If(MQC_Length_Vector(VectorOut).eq.0 .or. &
          MQC_Length_Vector(VectorOut).ne.MQC_Length_Vector(VectorIn)) &
          Call MQC_Allocate_Vector(MQC_Length_Vector(VectorIn),VectorOut,'Integer')
        VectorOut%veci = VectorIn%veci
      EndIf
      If(MQC_Vector_HaveComplex(VectorIn)) then
        If(MQC_Length_Vector(VectorOut).eq.0 .or. &
          MQC_Length_Vector(VectorOut).ne.MQC_Length_Vector(VectorIn)) &
          Call MQC_Allocate_Vector(MQC_Length_Vector(VectorIn),VectorOut,'Complex')
        VectorOut%vecc = VectorIn%vecc
      EndIf
!
      VectorOut%Column = VectorIn%Column
!
      Return
      End Subroutine MQC_Set_Vector2Vector
!
!
!     PROCEDURE MQC_VectorVectorSum
      Function MQC_VectorVectorSum(Vector1In,Vector2In)
!
!     This routine is used to add two MQC vectors and store them into
!     another MQC vector.
!
!     Lee M. Thompson 2016.
!
      Implicit None
      Type(MQC_Vector)::MQC_VectorVectorSum
      Type(MQC_Vector),Intent(In)::Vector1In, Vector2In
      Integer(kind=int64) :: Vector1Len, Vector2Len, VectorOutLen

      Vector1Len = MQC_Length_Vector(Vector1In)
      Vector2Len = MQC_Length_Vector(Vector2In)

      If (Vector1Len.ne.Vector2Len) then
        Call MQC_Error_I('Vector lengths unequal in MQC_VectorVectorSum', 6, &
             'Vector1Len', Vector1Len, &
             'Vector2Len', Vector2Len )
      EndIf

      If (Vector1In%column.neqv.Vector2In%column) then
        Call MQC_Error_L('Vector orientations are different in MQC_VectorVectorSum', 6, &
             'Vector1In%column', Vector1In%column, &
             'Vector2In%column', Vector2In%column )
      EndIf

      If(Vector1In%Data_type.eq.'Real') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Real')
          MQC_VectorVectorSum%vecr = Vector1In%vecr + Vector2In%vecr
          MQC_VectorVectorSum%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Real')
          MQC_VectorVectorSum%vecr = Vector1In%vecr + Vector2In%veci
          MQC_VectorVectorSum%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Complex')
          MQC_VectorVectorSum%vecc = Vector1In%vecr + Vector2In%vecc
          MQC_VectorVectorSum%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_VectorVectorSum', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      ElseIf(Vector1In%Data_type.eq.'Integer') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Real')
          MQC_VectorVectorSum%vecr = Vector1In%veci + Vector2In%vecr
          MQC_VectorVectorSum%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Integer')
          MQC_VectorVectorSum%veci = Vector1In%veci + Vector2In%veci
          MQC_VectorVectorSum%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Complex')
          MQC_VectorVectorSum%vecc = Vector1In%veci + Vector2In%vecc
          MQC_VectorVectorSum%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_VectorVectorSum', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      ElseIf(Vector1In%Data_type.eq.'Complex') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Complex')
          MQC_VectorVectorSum%vecc = Vector1In%vecc + Vector2In%vecr
          MQC_VectorVectorSum%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Complex')
          MQC_VectorVectorSum%vecc = Vector1In%vecc + Vector2In%veci
          MQC_VectorVectorSum%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorSum,'Complex')
          MQC_VectorVectorSum%vecc = Vector1In%vecc + Vector2In%vecc
          MQC_VectorVectorSum%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_VectorVectorSum', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      Else
        Call MQC_Error_A('Vector1In type unspecified in MQC_VectorVectorSum', 6, &
             'Vector1In%Data_type', Vector1In%Data_type )
      EndIf

      Return
      End Function MQC_VectorVectorSum
!
!
!     PROCEDURE MQC_VectorVectorDifference
      Function MQC_VectorVectorDifference(Vector1In,Vector2In)
!
!     This routine is used to subtract two MQC vectors and store them into
!     another MQC vector.
!
!     Lee M. Thompson 2016.
!
      Implicit None
      Type(MQC_Vector)::MQC_VectorVectorDifference
      Type(MQC_Vector),Intent(In)::Vector1In, Vector2In
      Integer(kind=int64) :: Vector1Len, Vector2Len, VectorOutLen

      Vector1Len = MQC_Length_Vector(Vector1In)
      Vector2Len = MQC_Length_Vector(Vector2In)

      If (Vector1Len.ne.Vector2Len) then
        Call MQC_Error_I('Vector lengths unequal in MQC_VectorVectorDifference', 6, &
             'Vector1Len', Vector1Len, &
             'Vector2Len', Vector2Len )
      EndIf

      If (Vector1In%column.neqv.Vector2In%column) then
        Call MQC_Error_L('Vector orientations are different in MQC_VectorVectorDifference', 6, &
             'Vector1In%column', Vector1In%column, &
             'Vector2In%column', Vector2In%column )
      EndIf

      If(Vector1In%Data_type.eq.'Real') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Real')
          MQC_VectorVectorDifference%vecr = Vector1In%vecr - Vector2In%vecr
          MQC_VectorVectorDifference%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Real')
          MQC_VectorVectorDifference%vecr = Vector1In%vecr - Vector2In%veci
          MQC_VectorVectorDifference%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Complex')
          MQC_VectorVectorDifference%vecc = Vector1In%vecr - Vector2In%vecc
          MQC_VectorVectorDifference%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_VectorVectorDifference', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      ElseIf(Vector1In%Data_type.eq.'Integer') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Real')
          MQC_VectorVectorDifference%vecr = Vector1In%veci - Vector2In%vecr
          MQC_VectorVectorDifference%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Integer')
          MQC_VectorVectorDifference%veci = Vector1In%veci - Vector2In%veci
          MQC_VectorVectorDifference%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Complex')
          MQC_VectorVectorDifference%vecc = Vector1In%veci - Vector2In%vecc
          MQC_VectorVectorDifference%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_VectorVectorDifference', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      ElseIf(Vector1In%Data_type.eq.'Complex') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Complex')
          MQC_VectorVectorDifference%vecc = Vector1In%vecc - Vector2In%vecr
          MQC_VectorVectorDifference%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Complex')
          MQC_VectorVectorDifference%vecc = Vector1In%vecc - Vector2In%veci
          MQC_VectorVectorDifference%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_VectorVectorDifference,'Complex')
          MQC_VectorVectorDifference%vecc = Vector1In%vecc - Vector2In%vecc
          MQC_VectorVectorDifference%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_VectorVectorDifference', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      Else
        Call MQC_Error_A('Vector1In type unspecified in MQC_VectorVectorDifference', 6, &
             'Vector1In%Data_type', Vector1In%Data_type )
      EndIf

      Return
      End Function MQC_VectorVectorDifference
!
!
!     PROCEDURE MQC_ScalarVectorSum
      Function MQC_ScalarVectorSum(ScalarIn,VectorIn)
!
!     This routine is used to add an MQC scalar to all elements in a MQC vector and
!     store it into another MQC vector.
!
!     Lee M. Thompson 2016.
!
      Implicit None
      Type(MQC_Vector)::MQC_ScalarVectorSum
      Type(MQC_Scalar),Intent(In)::ScalarIn
      Type(MQC_Vector),Intent(In)::VectorIn
      Integer(kind=int64) :: VectorLen

      VectorLen = MQC_Length_Vector(VectorIn)

      If(VectorIn%Data_type.eq.'Real') then
        If(ScalarIn%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Real')
          MQC_ScalarVectorSum%vecr = VectorIn%vecr + ScalarIn%scar
          MQC_ScalarVectorSum%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Real')
          MQC_ScalarVectorSum%vecr = VectorIn%vecr + ScalarIn%scai
          MQC_ScalarVectorSum%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Complex')
          MQC_ScalarVectorSum%vecc = VectorIn%vecr + ScalarIn%scac
          MQC_ScalarVectorSum%column = VectorIn%column
        Else
          Call MQC_Error_A('ScalarIn type unspecified in MQC_ScalarVectorSum', 6,&
               'ScalarIn%Data_type', ScalarIn%Data_type )
        EndIf
      ElseIf(VectorIn%Data_type.eq.'Integer') then
        If(ScalarIn%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Real')
          MQC_ScalarVectorSum%vecr = VectorIn%veci + ScalarIn%scar
          MQC_ScalarVectorSum%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Integer')
          MQC_ScalarVectorSum%veci = VectorIn%veci + ScalarIn%scai
          MQC_ScalarVectorSum%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Complex')
          MQC_ScalarVectorSum%vecc = VectorIn%veci + ScalarIn%scac
          MQC_ScalarVectorSum%column = VectorIn%column
        Else
          Call MQC_Error_A('ScalarIn type unspecified in MQC_ScalarVectorSum', 6, &
               'ScalarIn%Data_type', ScalarIn%Data_type )
        EndIf
      ElseIf(VectorIn%Data_type.eq.'Complex') then
        If(ScalarIn%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Complex')
          MQC_ScalarVectorSum%vecc = VectorIn%vecc + ScalarIn%scar
          MQC_ScalarVectorSum%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Complex')
          MQC_ScalarVectorSum%vecc = VectorIn%vecc + ScalarIn%scai
          MQC_ScalarVectorSum%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorSum,'Complex')
          MQC_ScalarVectorSum%vecc = VectorIn%vecc + ScalarIn%scac
          MQC_ScalarVectorSum%column = VectorIn%column
        Else
          Call MQC_Error_A('ScalarIn type unspecified in MQC_ScalarVectorSum', 6, &
               'ScalarIn%Data_type', ScalarIn%Data_type )
        EndIf
      Else
        Call MQC_Error_A('VectorIn type unspecified in MQC_ScalarVectorSum', 6, &
             'VectorIn%Data_type', VectorIn%Data_type )
      EndIf

      Return
      End Function MQC_ScalarVectorSum
!
!
!     PROCEDURE MQC_ScalarVectorDifference
      Function MQC_ScalarVectorDifference(ScalarIn,VectorIn)
!
!     This routine is used to subtract an MQC scalar to all elements in a MQC vector and
!     store it into another MQC vector.
!
!     Lee M. Thompson 2016.
!
      Implicit None
      Type(MQC_Vector)::MQC_ScalarVectorDifference
      Type(MQC_Scalar),Intent(In)::ScalarIn
      Type(MQC_Vector),Intent(In)::VectorIn
      Integer(kind=int64) :: VectorLen

      VectorLen = MQC_Length_Vector(VectorIn)

      If(VectorIn%Data_type.eq.'Real') then
        If(ScalarIn%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Real')
          MQC_ScalarVectorDifference%vecr = VectorIn%vecr - ScalarIn%scar
          MQC_ScalarVectorDifference%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Real')
          MQC_ScalarVectorDifference%vecr = VectorIn%vecr - ScalarIn%scai
          MQC_ScalarVectorDifference%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Complex')
          MQC_ScalarVectorDifference%vecc = VectorIn%vecr - ScalarIn%scac
          MQC_ScalarVectorDifference%column = VectorIn%column
        Else
          Call MQC_Error_A('ScalarIn type unspecified in MQC_ScalarVectorDifference', 6, &
               'ScalarIn%Data_type', ScalarIn%Data_type, &
               'VectorIn%Data_type', VectorIn%Data_type )
        EndIf
      ElseIf(VectorIn%Data_type.eq.'Integer') then
        If(ScalarIn%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Real')
          MQC_ScalarVectorDifference%vecr = VectorIn%veci - ScalarIn%scar
          MQC_ScalarVectorDifference%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Integer')
          MQC_ScalarVectorDifference%veci = VectorIn%veci - ScalarIn%scai
          MQC_ScalarVectorDifference%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Complex')
          MQC_ScalarVectorDifference%vecc = VectorIn%veci - ScalarIn%scac
          MQC_ScalarVectorDifference%column = VectorIn%column
        Else
          Call MQC_Error_A('ScalarIn type unspecified in MQC_ScalarVectorDifference', 6, &
               'ScalarIn%Data_type', ScalarIn%Data_type, &
               'VectorIn%Data_type', VectorIn%Data_type )
        EndIf
      ElseIf(VectorIn%Data_type.eq.'Complex') then
        If(ScalarIn%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Complex')
          MQC_ScalarVectorDifference%vecc = VectorIn%vecc - ScalarIn%scar
          MQC_ScalarVectorDifference%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Complex')
          MQC_ScalarVectorDifference%vecc = VectorIn%vecc - ScalarIn%scai
          MQC_ScalarVectorDifference%column = VectorIn%column
        ElseIf(ScalarIn%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(VectorLen,MQC_ScalarVectorDifference,'Complex')
          MQC_ScalarVectorDifference%vecc = VectorIn%vecc - ScalarIn%scac
          MQC_ScalarVectorDifference%column = VectorIn%column
        Else
          Call MQC_Error_A('ScalarIn type unspecified in MQC_ScalarVectorDifference', 6, &
               'ScalarIn%Data_type', ScalarIn%Data_type, &
               'VectorIn%Data_type', VectorIn%Data_type )
        EndIf
      Else
        Call MQC_Error_A('VectorIn type unspecified in MQC_ScalarVectorDifference', 6, &
             'VectorIn%Data_type', VectorIn%Data_type )
      EndIf

      Return
      End Function MQC_ScalarVectorDifference
!
!
!     PROCEDURE MQC_ElementVectorProduct
      Function MQC_ElementVectorProduct(Vector1In,Vector2In)
!
!     This routine is used to multiply two MQC vectors elementwise and store them
!     into another MQC vector.
!
!     Lee M. Thompson 2016.
!
      Implicit None
      Type(MQC_Vector)::MQC_ElementVectorProduct
      Type(MQC_Vector),Intent(In)::Vector1In, Vector2In
      Integer(kind=int64) :: Vector1Len, Vector2Len, VectorOutLen, I

      Vector1Len = MQC_Length_Vector(Vector1In)
      Vector2Len = MQC_Length_Vector(Vector2In)

      If (Vector1Len.ne.Vector2Len) then
        Call MQC_Error_I('Vector lengths unequal in MQC_ElementVectorProduct', 6, &
             'Vector1Len', Vector1Len, &
             'Vector2Len', Vector2Len )
      EndIf

      If (Vector1In%column.neqv.Vector2In%column) then
        Call MQC_Error_L('Vector orientations are different in ElementVectorProduct', 6, &
             'Vector1In%column', Vector1In%column, &
             'Vector2In%column', Vector2In%column )
      EndIf

      If(Vector1In%Data_type.eq.'Real') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Real')
          MQC_ElementVectorProduct%vecr = Vector1In%vecr * Vector2In%vecr
          MQC_ElementVectorProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Real')
          MQC_ElementVectorProduct%vecr = Vector1In%vecr * Vector2In%veci
          MQC_ElementVectorProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Complex')
          MQC_ElementVectorProduct%vecc = Vector1In%vecr * Vector2In%vecc
          MQC_ElementVectorProduct%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_ElementVectorProduct', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      ElseIf(Vector1In%Data_type.eq.'Integer') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Real')
          MQC_ElementVectorProduct%vecr = Vector1In%veci * Vector2In%vecr
          MQC_ElementVectorProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Integer')
          MQC_ElementVectorProduct%veci = Vector1In%veci * Vector2In%veci
          MQC_ElementVectorProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Complex')
          MQC_ElementVectorProduct%vecc = Vector1In%veci * Vector2In%vecc
          MQC_ElementVectorProduct%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_ElementVectorProduct', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      ElseIf(Vector1In%Data_type.eq.'Complex') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Complex')
          MQC_ElementVectorProduct%vecc = Vector1In%vecc * Vector2In%vecr
          MQC_ElementVectorProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Complex')
          MQC_ElementVectorProduct%vecc = Vector1In%vecc * Vector2In%veci
          MQC_ElementVectorProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_ElementVectorProduct,'Complex')
          MQC_ElementVectorProduct%vecc = Vector1In%vecc * Vector2In%vecc
          MQC_ElementVectorProduct%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_ElementVectorProduct', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      Else
        Call MQC_Error_A('Vector1In type unspecified in MQC_ElementVectorProduct', 6, &
             'Vector1In%Data_type', Vector1In%Data_type )
      EndIf

      Return
      End Function MQC_ElementVectorProduct
!
!
!     PROCEDURE MQC_Vector_Transpose
      Function MQC_Vector_Transpose(Vector)
!
!     This function returns the transpose of an MQC vector.
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::MQC_Vector_Transpose
      Class(MQC_Vector),Intent(In)::Vector
!
      MQC_Vector_Transpose = Vector
      MQC_Vector_Transpose%Column =  &
        .not.MQC_Vector_Transpose%Column
!
      Return
      End Function MQC_Vector_Transpose
!
!
!     PROCEDURE MQC_Vector_Conjugate_Transpose
      Function MQC_Vector_Conjugate_Transpose(Vector)
!
!     This function returns the conjugate transpose of an MQC vector.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::MQC_Vector_Conjugate_Transpose
      Class(MQC_Vector),Intent(In)::Vector
!
      MQC_Vector_Conjugate_Transpose = Vector
      If(MQC_Vector_HaveComplex(MQC_Vector_Conjugate_Transpose)) &
        MQC_Vector_Conjugate_Transpose%vecC = conjg(MQC_Vector_Conjugate_Transpose%vecC)
      MQC_Vector_Conjugate_Transpose%Column =  &
        .not.MQC_Vector_Conjugate_Transpose%Column
!
      Return
      End Function MQC_Vector_Conjugate_Transpose
!
!
!     PROCEDURE MQC_VectorVectorDotProduct
      Function MQC_VectorVectorDotProduct(Vector1,Vector2)
!
!     This function is used to compute the dot-product of two MQC vectors.
!
!     H. P. Hratchian, 2016.
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
!      Real,Intent(Out)::MQC_DotProduct
      Type(MQC_Scalar)::MQC_VectorVectorDotProduct
      Type(MQC_Vector),Intent(In)::Vector1,Vector2
!
      If (Vector1%Length.ne.Vector2%Length) then
        Call MQC_Error_I('Vector dimensions are different in MQC_VectorVectorDotProduct', 6, &
             'Vector1%Length', Vector1%Length, &
             'Vector2%Length', Vector2%Length )
      EndIf
      If (Vector1%Column.eqv.Vector2%Column) then
        Call MQC_Error_L('Vector orientations are the same in MQC_VectorVectorDotProduct', 6, &
             'Vector1%Column', Vector1%Column, &
             'Vector2%Column', Vector2%Column )
      EndIf
      if (vector1%column) then
        call mqc_error_L('Vector1 is column type and vector2 is not. &
     &    Use mqc_outer(Vector1,Vector2) for outer product', 6, &
     'vector1%column', vector1%column )
      endif
      If(MQC_Vector_HaveReal(Vector1)) then
        If(MQC_Vector_HaveReal(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%vecr,Vector2%vecr)
        ElseIf(MQC_Vector_HaveInteger(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%vecr,Vector2%veci)
        ElseIf(MQC_Vector_HaveComplex(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%vecr,Vector2%vecc)
        Else
          Call MQC_Error_L('Vector2In type unspecified in MQC_VectorVectorDotProduct', 6, &
               'MQC_Vector_HaveReal(Vector2)', MQC_Vector_HaveReal(Vector2), &
               'MQC_Vector_HaveInteger(Vector2)', MQC_Vector_HaveInteger(Vector2), &
               'MQC_Vector_HaveComplex(Vector2)', MQC_Vector_HaveComplex(Vector2) )
        EndIf
      ElseIf(MQC_Vector_HaveInteger(Vector1)) then
        If(MQC_Vector_HaveReal(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%veci,Vector2%vecr)
        ElseIf(MQC_Vector_HaveInteger(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%veci,Vector2%veci)
        ElseIf(MQC_Vector_HaveComplex(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%veci,Vector2%vecc)
        Else
          Call MQC_Error_L('Vector2In type unspecified in MQC_VectorVectorDotProduct', 6, &
               'MQC_Vector_HaveReal(Vector2)', MQC_Vector_HaveReal(Vector2), &
               'MQC_Vector_HaveInteger(Vector2)', MQC_Vector_HaveInteger(Vector2), &
               'MQC_Vector_HaveComplex(Vector2)', MQC_Vector_HaveComplex(Vector2) )
        EndIf
      ElseIf(MQC_Vector_HaveComplex(Vector1)) then
        If(MQC_Vector_HaveReal(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%vecc,Vector2%vecr)
        ElseIf(MQC_Vector_HaveInteger(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%vecc,Vector2%veci)
        ElseIf(MQC_Vector_HaveComplex(Vector2)) then
          MQC_VectorVectorDotProduct = dot_product(Vector1%vecc,Vector2%vecc)
        Else
          Call MQC_Error_L('Vector2In type unspecified in MQC_VectorVectorDotProduct', 6, &
               'MQC_Vector_HaveReal(Vector2)', MQC_Vector_HaveReal(Vector2), &
               'MQC_Vector_HaveInteger(Vector2)', MQC_Vector_HaveInteger(Vector2), &
               'MQC_Vector_HaveComplex(Vector2)', MQC_Vector_HaveComplex(Vector2) )
        EndIf
      Else
        Call MQC_Error_L('Vector1In type unspecified in MQC_VectorVectorDotProduct', 6, &
        'MQC_Vector_HaveComplex(Vector1)', MQC_Vector_HaveComplex(Vector1), &
        'MQC_Vector_HaveInteger(Vector1)', MQC_Vector_HaveInteger(Vector1), &
        'MQC_Vector_HaveReal(Vector1)', MQC_Vector_HaveReal(Vector1) )
      EndIf
!
      Return
    End Function MQC_VectorVectorDotProduct

!   Procedure mqc_outer
    Function mqc_outer(VA,VB) result(M)
!
!   This function calculates the outer product of mqc_vector VA and VB,
!   returns a mqc_matrix M.
!   X. Sheng 2017.02
!
      implicit none
      type(mqc_vector),intent(in)::VA,VB
      type(mqc_matrix)::M
      integer(kind=int64)::i,j
      type(mqc_vector)::VC

      If (VA%Column.eqv.VB%Column) then
        Call MQC_Error_L('Vector orientations are the same in MQC_Outer', 6, &
             'VA%Column', VA%Column, &
             'VB%Column', VB%Column )
      EndIf
      if (VB%column) then
        call mqc_error_l('Vector1 is row type and vector2 is column. &
     &    Use .dot. as a dot product operator for dot product. ', 6, &
     'VB%column', VB%column )
      endif
      if (mqc_have_int(VA) .and. mqc_have_int(VB)) then
        call mqc_allocate_matrix(VA%length,VB%length,M,'Integer','StorFull')
        call dger(VA%length,VB%length,1.0,VA%veci,1,VB%veci,1,M%mati,VA%length)
      elseif((mqc_have_int(VA).or.mqc_have_real(VA)).and.(mqc_have_int(VB).or.mqc_have_real(VB))) then
        call mqc_allocate_matrix(VA%length,VB%length,M,'Real','StorFull')
        if (mqc_have_int(VA)) then
          VC = mqc_cast_real(VA)
          call dger(VA%length,VB%length,1.0,VC%vecr,1,VB%vecr,1,M%matr,VA%length)
        elseif (mqc_have_int(VB)) then
          VC = mqc_cast_real(VB)
          call dger(VA%length,VB%length,1.0,VA%vecr,1,VC%vecr,1,M%matr,VA%length)
        else
          call dger(VA%length,VB%length,1.0,VA%vecr,1,VB%vecr,1,M%matr,VA%length)
        endif
      else
        call mqc_allocate_matrix(VA%length,VB%length,M,'Complex','StorFull')
        if (mqc_have_int(VA)) then
          VC = mqc_cast_complex(VA)
          call zgerc(VA%length,VB%length,1.0,VC%vecc,1,VB%vecc,1,M%matc,VA%length)
        elseif (mqc_have_int(VB)) then
          VC = mqc_cast_complex(VB)
          call dger(VA%length,VB%length,1.0,VA%vecc,1,VC%vecc,1,M%matc,VA%length)
        elseif (mqc_have_real(VA)) then
          VC = mqc_cast_complex(VA)
          call zgerc(VA%length,VB%length,1.0,VC%vecc,1,VB%vecc,1,M%matc,VA%length)
        elseif (mqc_have_real(VB)) then
          VC = mqc_cast_complex(VB)
          call dger(VA%length,VB%length,1.0,VA%vecc,1,VC%vecc,1,M%matc,VA%length)
        else
          call dger(VA%length,VB%length,1.0,VA%vecc,1,VB%vecc,1,M%matc,VA%length)
        endif
      endif

    return
    end function mqc_outer
!
!
!     PROCEDURE MQC_CrossProduct
      Function MQC_CrossProduct(Vector1In,Vector2In)
!
!     This function is used to compute the cross-product of two 3D MQC vectors.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Vector)::MQC_CrossProduct
      Type(MQC_Vector),Intent(In)::Vector1In,Vector2In
      Integer(kind=int64) :: Vector1Len, Vector2Len, VectorOutLen, I

      Vector1Len = MQC_Length_Vector(Vector1In)
      Vector2Len = MQC_Length_Vector(Vector2In)

      If (Vector1Len.ne.3.or.Vector2Len.ne.3) then
        Call MQC_Error_I('Vector not R3 in MQC_CrossProduct', 6, &
             'Vector1Len', Vector1Len, &
             'Vector2Len', Vector2Len )
      EndIf

      If (Vector1Len.ne.Vector2Len) then
        Call MQC_Error_I('Vector lengths unequal in MQC_CrossProduct', 6, &
             'Vector1Len', Vector1Len, &
             'Vector2Len', Vector2Len )
      EndIf

      If (Vector1In%column.neqv.Vector2In%column) then
        Call MQC_Error_L('Vector orientations are different in MQC_CrossProduct', 6, &
             'Vector1In%column', Vector1In%column, &
             'Vector2In%column', Vector2In%column )
      EndIf

      If(Vector1In%Data_type.eq.'Real') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Real')
          MQC_CrossProduct%vecr(1) = Vector1In%vecr(2)*Vector2In%vecr(3) - &
            Vector1In%vecr(3)*Vector2In%vecr(2)
          MQC_CrossProduct%vecr(2) = Vector1In%vecr(3)*Vector2In%vecr(1) - &
            Vector1In%vecr(1)*Vector2In%vecr(3)
          MQC_CrossProduct%vecr(3) = Vector1In%vecr(1)*Vector2In%vecr(2) - &
            Vector1In%vecr(2)*Vector2In%vecr(1)
          MQC_CrossProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Real')
          MQC_CrossProduct%vecr(1) = Vector1In%vecr(2)*Vector2In%veci(3) - &
            Vector1In%vecr(3)*Vector2In%veci(2)
          MQC_CrossProduct%vecr(2) = Vector1In%vecr(3)*Vector2In%veci(1) - &
            Vector1In%vecr(1)*Vector2In%veci(3)
          MQC_CrossProduct%vecr(3) = Vector1In%vecr(1)*Vector2In%veci(2) - &
            Vector1In%vecr(2)*Vector2In%veci(1)
          MQC_CrossProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Complex')
          MQC_CrossProduct%vecc(1) = Vector1In%vecr(2)*Vector2In%vecc(3) - &
            Vector1In%vecr(3)*Vector2In%vecc(2)
          MQC_CrossProduct%vecc(2) = Vector1In%vecr(3)*Vector2In%vecc(1) - &
            Vector1In%vecr(1)*Vector2In%vecc(3)
          MQC_CrossProduct%vecc(3) = Vector1In%vecr(1)*Vector2In%vecc(2) - &
            Vector1In%vecr(2)*Vector2In%vecc(1)
          MQC_CrossProduct%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_CrossProduct', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      ElseIf(Vector1In%Data_type.eq.'Integer') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Real')
          MQC_CrossProduct%vecr(1) = Vector1In%veci(2)*Vector2In%vecr(3) - &
            Vector1In%veci(3)*Vector2In%vecr(2)
          MQC_CrossProduct%vecr(2) = Vector1In%veci(3)*Vector2In%vecr(1) - &
            Vector1In%veci(1)*Vector2In%vecr(3)
          MQC_CrossProduct%vecr(3) = Vector1In%veci(1)*Vector2In%vecr(2) - &
            Vector1In%veci(2)*Vector2In%vecr(1)
          MQC_CrossProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Integer')
          MQC_CrossProduct%veci(1) = Vector1In%veci(2)*Vector2In%veci(3) - &
            Vector1In%veci(3)*Vector2In%veci(2)
          MQC_CrossProduct%veci(2) = Vector1In%veci(3)*Vector2In%veci(1) - &
            Vector1In%veci(1)*Vector2In%veci(3)
          MQC_CrossProduct%veci(3) = Vector1In%veci(1)*Vector2In%veci(2) - &
            Vector1In%veci(2)*Vector2In%veci(1)
          MQC_CrossProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Complex')
          MQC_CrossProduct%vecc(1) = Vector1In%veci(2)*Vector2In%vecc(3) - &
            Vector1In%veci(3)*Vector2In%vecc(2)
          MQC_CrossProduct%vecc(2) = Vector1In%veci(3)*Vector2In%vecc(1) - &
            Vector1In%veci(1)*Vector2In%vecc(3)
          MQC_CrossProduct%vecc(3) = Vector1In%veci(1)*Vector2In%vecc(2) - &
            Vector1In%veci(2)*Vector2In%vecc(1)
          MQC_CrossProduct%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_CrossProduct', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      ElseIf(Vector1In%Data_type.eq.'Complex') then
        If(Vector2In%Data_type.eq.'Real') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Complex')
          MQC_CrossProduct%vecc(1) = Vector1In%vecc(2)*Vector2In%vecr(3) - &
            Vector1In%vecc(3)*Vector2In%vecr(2)
          MQC_CrossProduct%vecc(2) = Vector1In%vecc(3)*Vector2In%vecr(1) - &
            Vector1In%vecc(1)*Vector2In%vecr(3)
          MQC_CrossProduct%vecc(3) = Vector1In%vecc(1)*Vector2In%vecr(2) - &
            Vector1In%vecc(2)*Vector2In%vecr(1)
          MQC_CrossProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Integer') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Complex')
          MQC_CrossProduct%vecc(1) = Vector1In%vecc(2)*Vector2In%veci(3) - &
            Vector1In%vecc(3)*Vector2In%veci(2)
          MQC_CrossProduct%vecc(2) = Vector1In%vecc(3)*Vector2In%veci(1) - &
            Vector1In%vecc(1)*Vector2In%veci(3)
          MQC_CrossProduct%vecc(3) = Vector1In%vecc(1)*Vector2In%veci(2) - &
            Vector1In%vecc(2)*Vector2In%veci(1)
          MQC_CrossProduct%column = Vector1In%column
        ElseIf(Vector2In%Data_type.eq.'Complex') then
          Call MQC_Allocate_Vector(Vector1Len,MQC_CrossProduct,'Complex')
          MQC_CrossProduct%vecc(1) = Vector1In%vecc(2)*Vector2In%vecc(3) - &
            Vector1In%vecc(3)*Vector2In%vecc(2)
          MQC_CrossProduct%vecc(2) = Vector1In%vecc(3)*Vector2In%vecc(1) - &
            Vector1In%vecc(1)*Vector2In%vecc(3)
          MQC_CrossProduct%vecc(3) = Vector1In%vecc(1)*Vector2In%vecc(2) - &
            Vector1In%vecc(2)*Vector2In%vecc(1)
          MQC_CrossProduct%column = Vector1In%column
        Else
          Call MQC_Error_A('Vector2In type unspecified in MQC_CrossProduct', 6, &
               'Vector2In%Data_type', Vector2In%Data_type )
        EndIf
      Else
        Call MQC_Error_A('Vector1In type unspecified in MQC_CrossProduct', 6, &
             'Vector1In%Data_type', Vector1In%Data_type )
      EndIf
      Return
      End Function MQC_CrossProduct
!
!
!     PROCEDURE MQC_Print_Vector_Algebra1
      Subroutine MQC_Print_Vector_Algebra1(Vector,IOut,Header,Verbose,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a MQC_Vector type variable.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Class(MQC_Vector),Intent(In)::Vector
      Character(Len=*),Intent(In)::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom,Verbose
      Integer(kind=int64)::I,Length
!
 1000 Format(1x,A)
 1001 Format(1x,I7,2x,I14)
 1002 Format(1x,I7,2X,F14.6)
 1003 Format(1x,I7,2X,F12.5,F11.5,"i")
 1020 Format( " " )
 1030 Format( A )
!
      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      Write(IOut,1000) TRIM(Header)

      Length = MQC_Length_Vector(Vector)
      if (Present(Verbose)) then
        if (verbose.and.vector%column) then
          write(IOut,1030) 'The vector being printed is a column vector'
        elseif (verbose) then
          write(IOut,1030) 'The vector being printed is a row vector'
        endif
      endIf
      If(Vector%Data_type.eq.'Integer') then
          Do I = 1,Length
            Write(IOut,1001) I, Vector%veci(I)
          EndDo
      ElseIf(Vector%Data_type.eq.'Real') then
        Do I = 1,Length
          Write(IOut,1002) I, Vector%vecr(I)
        EndDo
      ElseIf(Vector%Data_type.eq.'Complex') then
        Do I = 1,Length
          Write(IOut,1003) I, Vector%vecc(I)
        EndDo
      Else
        Call MQC_Error_A('Vector type unspecified in MQC_Print_Vector_Algebra1', 6, &
             'Vector%Data_type', Vector%Data_type )
      EndIf

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Vector_Algebra1
!
!
!     PROCEDURE MQC_vector_cast_real
      function MQC_vector_cast_real(VA) result(VB)
!     This function converts a mqc_vector to real. Only integer vector is allowed
!     because it would lose accuracy otherwise.
!
!     X. Sheng, 2017.02
!
      implicit none
      type(mqc_vector),intent(in)::VA
      type(mqc_vector)::VB

      call mqc_allocate_vector(VA%length,VB,'Real')
      VB%vecr = VA%veci

      end function MQC_vector_cast_real
!
!
!     PROCEDURE MQC_vector_cast_complex
      function MQC_vector_cast_complex(VA) result(VB)
!     This function converts a mqc_vector to complex. 
!
!     X. Sheng, 2017.02
!
      implicit none
      type(mqc_vector),intent(in)::VA
      type(mqc_vector)::VB

      call mqc_allocate_vector(VA%length,VB,'Complex')
      if(mqc_have_int(VA)) then
        VB%vecc = VA%veci
      elseif(mqc_have_real(VA)) then
        VB%vecc = VA%vecr
      endif

      end function MQC_vector_cast_complex
!
!
!     PROCEDURE MQC_Vector_Scalar_Put
      Subroutine MQC_Vector_Scalar_Put(Vector,Scalar,I)
!
!     This subroutine updates the value of MQC_Vector Vector with the value
!     in MQC_Scalar Scalar at element I.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      Implicit None
      Integer(kind=int64),Intent(In)::I
      Class(MQC_Vector),Intent(InOut)::Vector
      Type(MQC_Scalar),Intent(In)::Scalar
      Integer(kind=int64)::IndI

      IndI = I
      If (IndI.lt.0) IndI = Vector%Length + IndI + 1
      If (IndI.eq.0.or.IndI.gt.Vector%Length) Call MQC_Error_I('Index I badly specified in mqc_vector_scalar_put', 6, &
           'IndI', IndI, &
           'Vector%Length', Vector%Length )

      If (Vector%Data_Type.eq.'Integer') then
        If (Scalar%Data_Type.eq.'Integer') then
          Vector%VecI(IndI) = Scalar%ScaI
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Call MQC_Vector_Copy_Int2Real(Vector)
          Vector%VecR(IndI) = Scalar%ScaR
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Call MQC_Vector_Copy_Int2Complex(Vector)
          Vector%VecC(IndI) = Scalar%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Vector_Scalar_Put', 6, &
               'Scalar%Data_Type', Scalar%Data_Type )
        EndIf
      ElseIf (Vector%Data_Type.eq.'Real') then
        If (Scalar%Data_Type.eq.'Integer') then
          Vector%VecR(IndI) = Scalar%ScaI
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Vector%VecR(IndI) = Scalar%ScaR
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Call MQC_Vector_Copy_Real2Complex(Vector)
          Vector%VecC(IndI) = Scalar%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Vector_Scalar_Put', 6, &
               'Scalar%Data_Type', Scalar%Data_Type )
        EndIf
      ElseIf (Vector%Data_Type.eq.'Complex') then
        If (Scalar%Data_Type.eq.'Integer') then
          Vector%VecC(IndI) = Scalar%ScaI
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Vector%VecC(IndI) = Scalar%ScaR
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Vector%VecC(IndI) = Scalar%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Vector_Scalar_Put', 6, &
               'Scalar%Data_Type', Scalar%Data_Type )
        EndIf
      Else
        Call MQC_Error_A('Vector type not defined in MQC_Vector_Scalar_Put', 6, &
             'Vector%Data_Type', Vector%Data_Type )
      EndIf

      End Subroutine MQC_Vector_Scalar_Put
!
!
!     PROCEDURE MQC_Vector_Scalar_Increment
      Subroutine MQC_Vector_Scalar_Increment(Vector,Scalar,I)
!
!     This subroutine increments the value of MQC_Vector <Vector> by <scalar> at
!     element <I>.
!
!
!     Variable Declarations.
!
      Implicit None
      Integer(kind=int64),Intent(In)::I
      Class(MQC_Vector),Intent(InOut)::Vector
      Type(MQC_Scalar),Intent(In)::Scalar
      Integer(kind=int64)::IndI

      IndI = I
      If (IndI.lt.0) IndI = Vector%Length + IndI + 1
      If (IndI.eq.0.or.IndI.gt.Vector%Length) Call MQC_Error_I('Index I badly specified in mqc_vector_scalar_increment.', 6, &
           'IndI', IndI, &
           'Vector%Length', Vector%Length )
      If (Vector%Data_Type.eq.'Integer') then
        If (Scalar%Data_Type.eq.'Integer') then
          Vector%VecI(IndI) = Vector%VecI(IndI) + Scalar%ScaI
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Vector%VecI(IndI) = Vector%VecI(IndI) + int(Scalar%ScaR)
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Vector%VecI(IndI) = Vector%VecI(IndI) + Real(Scalar%ScaC)
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Vector_Scalar_Incrememt', 6, &
               'Scalar%Data_Type', Scalar%Data_Type )
        EndIf
      ElseIf (Vector%Data_Type.eq.'Real') then
        If (Scalar%Data_Type.eq.'Integer') then
#ifdef PGI
          Vector%VecR(IndI) = Vector%VecR(IndI) + dfloat(Scalar%ScaI)
#else
          Vector%VecR(IndI) = Vector%VecR(IndI) + float(Scalar%ScaI)
#endif
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Vector%VecR(IndI) = Vector%VecR(IndI) + Scalar%ScaR
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Vector%VecR(IndI) = Vector%VecR(IndI) + Real(Scalar%ScaC)
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Vector_Scalar_Increment', 6, &
               'Scalar%Data_Type', Scalar%Data_Type )
        EndIf
      ElseIf (Vector%Data_Type.eq.'Complex') then
        If (Scalar%Data_Type.eq.'Integer') then
          Vector%VecC(IndI) = Vector%VecC(IndI) + cmplx(Scalar%ScaI,0)
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Vector%VecC(IndI) = Vector%VecC(IndI) + cmplx(Scalar%ScaR,0.0)
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Vector%VecC(IndI) = Vector%VecC(IndI) + Scalar%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Vector_Scalar_Increment', 6, &
               'Scalar%Data_Type', Scalar%Data_Type )
        EndIf
      Else
        Call MQC_Error_A('Vector type not defined in MQC_Vector_Scalar_Increment', 6, &
'Vector%Data_Type', Vector%Data_Type )
      EndIf
!
      End Subroutine MQC_Vector_Scalar_Increment
!
!
!     PROCEDURE MQC_Vector_Vector_Put
      Subroutine MQC_Vector_Vector_Put(Vector,VectorIn,I)
!
!     This subroutine updates the value of MQC_Vector Vector with the value
!     in MQC_Vector VectorIn elements between I and J. Negative I and J count
!     from end of vector.
!
!     L. M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Vector),Intent(InOut)::Vector
      Type(MQC_Vector),Intent(In)::VectorIn
      Integer(kind=int64),Optional,Intent(In)::I
      Integer(kind=int64)::Length,IndI,IndJ

      If(Present(I)) then
        IndI = I
      Else
        IndI = 1
      EndIf
      If (IndI.lt.0) IndI = Vector%Length + IndI + 1
      If (IndI.eq.0) Call MQC_Error_I('Index I badly specified in mqc_vector_vector_put', 6, &
           'IndI', IndI )
      IndJ = IndI + MQC_Length_Vector(VectorIn) - 1

      Length = IndJ-IndI+1
      If (Length.le.0.or.Length.gt.Vector%Length.or.Length.ne.VectorIn%Length) Call MQC_Error_I('Vector length badly &
    &  defined in MQC_Vector_Vector_Put', 6, &
    'Length', Length, &
    'Vector%Length', Vector%Length, &
    'VectorIn%Length', VectorIn%Length )
      If (IndI.le.0.or.IndI.gt.(Vector%Length-Length+1)) Call MQC_Error_I('Index I out of bounds &
    &   in MQC_Vector_Vector_Put', 6, &
    'IndI', IndI, &
    'Vector%Length', Vector%Length, &
    'Length', Length )
      If (IndJ.lt.Length.or.IndJ.gt.Vector%Length) Call MQC_Error_I('Index J out of bounds &
    &   in MQC_Vector_Vector_Put', 6, &
    'IndJ', IndJ, &
    'Length', Length, &
    'Vector%Length', Vector%Length )
      If (Vector%Data_Type.eq.'Integer') then
        If (VectorIn%Data_Type.eq.'Integer') then
          Vector%VecI(IndI:IndJ) = VectorIn%VecI
        ElseIf (VectorIn%Data_Type.eq.'Real') then
          Call MQC_Vector_Copy_Int2Real(Vector)
          Vector%VecR(IndI:IndJ) = VectorIn%VecR
        ElseIf (VectorIn%Data_Type.eq.'Complex') then
          Call MQC_Vector_Copy_Int2Complex(Vector)
          Vector%VecC(IndI:IndJ) = VectorIn%VecC
        Else
          Call MQC_Error_A('VectorIn type not defined in MQC_Vector_Vector_Put', 6, &
               'VectorIn%Data_Type', VectorIn%Data_Type )
        EndIf
      ElseIf (Vector%Data_Type.eq.'Real') then
        If (VectorIn%Data_Type.eq.'Integer') then
          Vector%VecR(IndI:IndJ) = VectorIn%VecI
        ElseIf (VectorIn%Data_Type.eq.'Real') then
          Vector%VecR(IndI:IndJ) = VectorIn%VecR
        ElseIf (VectorIn%Data_Type.eq.'Complex') then
          Call MQC_Vector_Copy_Real2Complex(Vector)
          Vector%VecC(IndI:IndJ) = VectorIn%VecC
        Else
          Call MQC_Error_a('VectorIn type not defined in MQC_Vector_Vector_Put', 6, &
               'VectorIn%Data_Type', VectorIn%Data_Type )
        EndIf
      ElseIf (Vector%Data_Type.eq.'Complex') then
        If (VectorIn%Data_Type.eq.'Integer') then
          Vector%VecC(IndI:IndJ) = VectorIn%VecI
        ElseIf (VectorIn%Data_Type.eq.'Real') then
          Vector%VecC(IndI:IndJ) = VectorIn%VecR
        ElseIf (VectorIn%Data_Type.eq.'Complex') then
          Vector%VecC(IndI:IndJ) = VectorIn%VecC
        Else
          Call MQC_Error_A('VectorIn type not defined in MQC_Vector_Vector_Put', 6, &
               'VectorIn%Data_Type', VectorIn%Data_Type )
        EndIf
      Else
        Call MQC_Error_a('Vector type not defined in MQC_Vector_Vector_Put', 6, &
             'Vector%Data_Type', Vector%Data_Type )
      EndIf

      End Subroutine MQC_Vector_Vector_Put
!
!
!     PROCEDURE MQC_Vector_Initialize
      Subroutine MQC_Vector_Initialize(Vector,Length,Scalar)
!
!     This subroutine initializes an MQC_Vector Vector of dimension
!     Length with the value Scalar in all elements.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      Implicit None
      Real(kind=real64)::Zero=0.0d0
      Integer(kind=int64),Intent(In)::Length
      Class(*),Optional::Scalar
      Class(MQC_Vector),Intent(InOut)::Vector

      If(Present(Scalar)) then
        Select Type (Scalar)
        Type is (Integer)
          Call MQC_Allocate_Vector(Length,Vector,'Integer')
          Vector%VecI = Scalar
        Type is (Real)
          Call MQC_Allocate_Vector(Length,Vector,'Real')
          Vector%VecR = Scalar
        Type is (Complex(Kind=real64))
          Call MQC_Allocate_Vector(Length,Vector,'Complex')
          Vector%VecC = Scalar
        Type is (MQC_Scalar)
          If(MQC_Scalar_HaveInteger(Scalar)) then
            Call MQC_Allocate_Vector(Length,Vector,'Integer')
            Vector%VecI = Scalar%ival()
          ElseIf(MQC_Scalar_HaveReal(Scalar)) then
            Call MQC_Allocate_Vector(Length,Vector,'Real')
            Vector%VecR = Scalar%rval()
          ElseIf(MQC_Scalar_HaveComplex(Scalar)) then
            Call MQC_Allocate_Vector(Length,Vector,'Complex')
            Vector%VecC = Scalar%cval()
          EndIf
        Class Default
          Call MQC_Error_I('Scalar Type not defined in MQC_Vector_Initialize', 6 )
        End Select
      Else
        Call MQC_Allocate_Vector(Length,Vector,'Real')
        Vector%VecR = Zero
      EndIf

      End Subroutine MQC_Vector_Initialize
!
!
!     PROCEDURE MQC_ScalarVectorProduct
      Function MQC_ScalarVectorProduct(Scalar,Vector) Result(Vector_Res)
!     This function calculates the product of an mqc_scalar and an mqc_vector
!
!      X. Sheng 2017.02
!
        implicit none
        type(mqc_scalar),intent(in)::scalar
        type(mqc_vector),intent(in)::vector
        type(mqc_vector)::vector_res

        call mqc_deallocate_vector(vector_res)
        if (vector%Data_type.eq.'Integer' .and. scalar%Data_type.eq.'Integer') then
          call mqc_allocate_vector(vector%length,vector_res,'Integer')
          vector_res%veci = scalar%scai * vector%veci
        else if (vector%Data_type.eq.'Integer' .and. scalar%Data_type.eq.'Real') then
          call mqc_allocate_vector(vector%length,vector_res,'Real')
          vector_res%vecr = scalar%scar * vector%veci
        else if (vector%Data_type.eq.'Real' .and. scalar%Data_type.eq.'Integer') then
          call mqc_allocate_vector(vector%length,vector_res,'Real')
          vector_res%vecr = scalar%scai * vector%vecr
        else
          call mqc_allocate_vector(vector%length,vector_res,'Real')
          vector_res%vecr = scalar%scar * vector%vecr
        endif

      end function mqc_ScalarVectorProduct
!
!
!     PROCEDURE MQC_VectorScalarProduct
      function mqc_VectorScalarProduct(vector,scalar) result(vector_res)
!
!     This function calculates the product of an mqc_scalar and an mqc_matrix
!
!     X. Sheng 2017.02
!
        implicit none
        type(mqc_scalar),intent(in)::scalar
        type(mqc_vector),intent(in)::vector
        type(mqc_vector)::vector_res

        call mqc_deallocate_vector(vector_res)
        if (vector%Data_type.eq.'Integer') then
          if(scalar%Data_type.eq.'Integer') then
            call mqc_allocate_vector(vector%length,vector_res,'Integer')
            vector_res%veci = scalar%scai * vector%veci
          elseIf(scalar%Data_type.eq.'Real') then
            call mqc_allocate_vector(vector%length,vector_res,'Real')
            vector_res%vecr = scalar%scar * vector%veci
          elseIf(scalar%Data_type.eq.'Complex') then
            call mqc_allocate_vector(vector%length,vector_res,'Complex')
            vector_res%vecc = scalar%scac * vector%veci
          else
            call mqc_error_a('unrecognised data type in mqc_vectorscalarproduct', 6, &
                 'scalar%Data_type', scalar%Data_type )
          endIf
        elseIf(vector%Data_type.eq.'Real') then
          if(scalar%Data_type.eq.'Integer') then
            call mqc_allocate_vector(vector%length,vector_res,'Real')
            vector_res%vecr = scalar%scai * vector%vecr
          elseIf(scalar%Data_type.eq.'Real') then
            call mqc_allocate_vector(vector%length,vector_res,'Real')
            vector_res%vecr = scalar%scar * vector%vecr
          elseIf(scalar%Data_type.eq.'Complex') then
            call mqc_allocate_vector(vector%length,vector_res,'Complex')
            vector_res%vecc = scalar%scac * vector%vecr
          else
            call mqc_error_a('unrecognised data type in mqc_vectorscalarproduct', 6, &
                 'scalar%Data_type', scalar%Data_type )
          endIf
        elseIf(vector%Data_type.eq.'Complex') then
          if(scalar%Data_type.eq.'Integer') then
            call mqc_allocate_vector(vector%length,vector_res,'Complex')
            vector_res%vecc = scalar%scai * vector%vecc
          elseIf(scalar%Data_type.eq.'Real') then
            call mqc_allocate_vector(vector%length,vector_res,'Complex')
            vector_res%vecc = scalar%scar * vector%vecc
          elseIf(scalar%Data_type.eq.'Complex') then
            call mqc_allocate_vector(vector%length,vector_res,'Complex')
            vector_res%vecc = scalar%scac * vector%vecr
          else
            call mqc_error_A('unrecognised data type in mqc_vectorscalarproduct', 6, &
                 'scalar%Data_type', scalar%Data_type )
          endIf
        else
          call mqc_error_a('unrecognised data type in mqc_vectorscalarproduct', 6, &
               'vector%Data_type', vector%Data_type )
        endif

      end function mqc_VectorScalarProduct
!
!
!     PROCEDURE MQC_Vector_Norm
      Function MQC_Vector_Norm(vector,methodIn) result(norm)
!
!     This function returns the norm of an mqc vector.
!     Different norms can be computed with optional argument methodIn.
!     'M' - max(abs(A(i,j)))
!     '1' - one norm
!     'I' - infinty norm
!     'F' - Frobenius norm (default)
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Vector),intent(inOut)::vector
      character(len=1),optional,intent(in)::methodIn
      character(len=1)::method
      character(len=64)::TypeFlag
      real(kind=real64),dimension(:),allocatable::work
      real(kind=real64)::dlange,zlange
      real(kind=real64),dimension(:),allocatable::temp
      Complex(Kind=real64),dimension(:),allocatable::tempc
      Type(MQC_Scalar)::norm

      if(Present(methodIn)) method = methodIn
      if(.not.Present(methodIn)) method = 'F'
      if(method.eq.'I') allocate(work(max(1,MQC_Length_Vector(vector))))
      TypeFlag = vector%Data_Type
      If(vector%Data_Type.eq.'Integer'.or.vector%Data_Type.eq.'Real') then 
        If(vector%Data_Type.eq.'Integer') Call MQC_Vector_Copy_Int2Real(vector)
        Allocate(temp(MQC_Length_Vector(vector)))
        temp = vector%vecr
        norm = dlange(method,MQC_Length_Vector(vector),1,temp, &
          max(MQC_Length_Vector(vector),1),work)
        if(allocated(work)) deallocate(work)
        if(allocated(temp)) deallocate(temp)
        if(TypeFlag.eq.'Integer') Call MQC_Vector_Copy_Real2Int(vector)
      ElseIf(vector%Data_Type.eq.'Complex') then
        Allocate(tempc(MQC_Length_Vector(vector)))
        tempc = vector%vecc
        norm = zlange(method,MQC_Length_Vector(vector),1,tempc, &
          max(MQC_Length_Vector(vector),1),work)
        if(allocated(work)) deallocate(work)
        if(allocated(tempc)) deallocate(tempc)
      Else
        call mqc_error_a('unrecognised data type in mqc_vector_norm', 6, &
             'vector%Data_Type', vector%Data_Type )
      EndIf

      End Function MQC_Vector_Norm
!
!     PROCEDURE MQC_Vector_isAllocated
      Function MQC_Vector_isAllocated(Vector) Result(isAllocated)
!
!     This function returns if an mqc vector is allocated or not.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Vector),Intent(InOut)::Vector
      Logical::isAllocated
!
      isAllocated = .False.
      If(Allocated(Vector%VecR)) isAllocated = .True.
      If(Allocated(Vector%VecI)) isAllocated = .True.
      If(Allocated(Vector%VecC)) isAllocated = .True.
!
      Return
      End Function MQC_Vector_isAllocated
!
!
!     PROCEDURE MQC_Vector_Push
      Subroutine MQC_Vector_Push(Vector,Scalar) 
!
!     This function adds a value to the end of a vector.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      implicit none
      class(mqc_vector),intent(inOut)::vector
      type(mqc_scalar),intent(in)::scalar
      type(mqc_vector)::temp
!
      if(mqc_length_vector(vector).eq.0) then
        call vector%init(1,scalar)
      else
        temp = vector
        call mqc_deallocate_vector(vector)
        call vector%init(mqc_length_vector(temp)+1)
        call vector%vput(temp)
        call vector%put(scalar,-1)
        call mqc_deallocate_vector(temp)
      endIf
!
      return
      end subroutine mqc_vector_push
!
!
!     PROCEDURE MQC_Vector_Unshift
      Subroutine MQC_Vector_Unshift(Vector,Scalar) 
!
!     This function adds a value to the beginning of a vector.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      implicit none
      class(mqc_vector),intent(inOut)::vector
      type(mqc_scalar),intent(in)::scalar
      type(mqc_vector)::temp
!
      if(mqc_length_vector(vector).eq.0) then
        call vector%init(1,scalar)
      else
        temp = vector
        call mqc_deallocate_vector(vector)
        call vector%init(mqc_length_vector(temp)+1)
        call vector%vput(temp,2)
        call vector%put(scalar,1)
        call mqc_deallocate_vector(temp)
      endIf
!
      return
      end subroutine mqc_vector_unshift
!
!
!     PROCEDURE MQC_Vector_Pop
      Function MQC_Vector_Pop(Vector) result(Output)
!
!     This function removes a value from the end of a vector and returns it.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      implicit none
      class(mqc_vector),intent(inOut)::vector
      type(mqc_scalar)::output
      type(mqc_vector)::temp
!
      if(mqc_length_vector(vector).eq.0) then
        return
      else
        output = vector%at(-1)
        temp = vector%vat(1,-2)
        call mqc_deallocate_vector(vector)
        vector = temp
        call mqc_deallocate_vector(temp)
      endIf
!
      return
      end function mqc_vector_pop 
!
!
!     PROCEDURE MQC_Vector_Shift
      Function MQC_Vector_Shift(Vector) result(Output)
!
!     This function removes a value from the beginning of a vector and returns it.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      implicit none
      class(mqc_vector),intent(inOut)::vector
      type(mqc_scalar)::output
      type(mqc_vector)::temp
!
      if(mqc_length_vector(vector).eq.0) then
        return
      else
        output = vector%at(1)
        temp = vector%vat(2,-1)
        call mqc_deallocate_vector(vector)
        vector = temp
        call mqc_deallocate_vector(temp)
      endIf
!
      return
      end function mqc_vector_shift
!
!
!     PROCEDURE MQC_Vector_MaxVal
      Function MQC_Vector_MaxVal(Vector) result(Output)
!
!     This function returns the largest value in a vector.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      implicit none
      class(mqc_vector),intent(in)::vector
      type(mqc_scalar)::output
!
      if(mqc_length_vector(vector).eq.0) then
        return
      else
        if(mqc_have_real(vector)) then
          output = maxval(vector%vecr)
        elseIf(mqc_have_int(vector)) then
          output = maxval(vector%veci)
        endIf
      endIf
!
      return
      end function mqc_vector_maxval
!
!
!     PROCEDURE MQC_Vector_MinVal
      Function MQC_Vector_MinVal(Vector) result(Output)
!
!     This function returns the smallest value in a vector.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      implicit none
      class(mqc_vector),intent(in)::vector
      type(mqc_scalar)::output
!
      if(mqc_length_vector(vector).eq.0) then
        return
      else
        if(mqc_have_real(vector)) then
          output = minval(vector%vecr)
        elseIf(mqc_have_int(vector)) then
          output = minval(vector%veci)
        endIf
      endIf
!
      return
      end function mqc_vector_minval
!
!
!     PROCEDURE MQC_Vector_Sqrt
      subroutine mqc_vector_sqrt(A)
!
!     This subroutine returns the square root of INPUT vector A.
!
      Implicit None
      class(mqc_vector),intent(inOut)::A
      integer(kind=int64)::i
!
!     Do the work.
!
      do i=1,MQC_Length_Vector(A)
        call A%put(sqrt(A%at(i)),i)
      endDo
!
      end subroutine mqc_vector_sqrt
!
!
!     PROCEDURE MQC_Vector_Power
      subroutine mqc_vector_power(A,P)
!
!     This subroutine returns the elements of INPUT vector A raised
!     to the power P.
!
      Implicit None
      class(mqc_vector),intent(inOut)::A
      class(*)::P
      integer(kind=int64)::i
      character(len=64)::data_type
      type(mqc_scalar)::scalar
!
!     Do the work.
!
      data_type = A%data_type

      select type(P)
      type is (integer)
        scalar = P
      type is (real)
        scalar = P
      type is (complex)
        scalar = P
      type is (mqc_scalar)
        scalar = P
      class default
        call mqc_error_I('P type not defined in MQC_Vector_Power',6)
      end select
!
      do i=1,MQC_Length_Vector(A)
        call A%put(A%at(i)**scalar,i)
      endDo
!
      end subroutine mqc_vector_power
!
!
!----------------------------------------------------------------
!                                                               |
!     MATRIX PROCEDURES                                         |
!                                                               |
!----------------------------------------------------------------
!
!     PROCEDURE MQC_Matrix_Diagonalize
      Subroutine MQC_Matrix_Diagonalize(A,EVals,EVecs)
!     This subroutine uses lapack routines to diagonalize a symmetric (Hermitian) MQC matrix, 
!     and returns eigen vectors as a MQC matrix and eigen values as a MQC vector.
!     X. Sheng 2017.03
!     L. M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Matrix),Intent(In)::A
      Type(MQC_Matrix)::B
      Type(MQC_Matrix),Optional,Intent(InOut)::EVecs
      Type(MQC_Vector),Optional,Intent(InOut)::EVals
      Real(kind=real64),Dimension(:,:),allocatable::A_EVecs
      Complex(Kind=real64),Dimension(:,:),allocatable::AC_EVecs
      Real(kind=real64),Dimension(:),allocatable::A_EVals
!
      Integer(kind=int64)::NDim,i,j,k,IError=1
      Real(kind=real64),Dimension(:),Allocatable::A_Symm,Temp_Vector
      Complex(Kind=real64),Dimension(:),Allocatable::AC_Symm,TempC_Vector
      character(len=64)::typeFlag
!
      typeFlag = a%Data_Type
      if (A%data_type.eq.'Real'.or.A%data_type.eq.'Integer') then
        if(A%data_type.eq.'Integer') call mqc_matrix_copy_int2Real(A)
        if (A%storage .eq. 'StorFull') then
          if(.not.mqc_matrix_test_symmetric(A)) call mqc_error_l('Input matrix not symmetric in MQC_Matrix_Diagonalize', 6, &
               'mqc_matrix_test_symmetric(A)', mqc_matrix_test_symmetric(A) )
          if (A%NRow .eq. A%NCol) then
            Allocate(A_Symm((A%NRow*(A%NRow+1))/2),Temp_Vector(3*A%NRow),A_EVecs(A%NRow,A%NRow),A_EVals(A%NRow))
            k = 0
            Do i = 1,A%NRow
              Do j = 1,i
                k = k+1
                A_Symm(k) = A%matR(i,j)
              EndDo
            EndDo
            Call DSPEV('V','U',A%NRow,A_Symm,A_EVals,A_EVecs,A%NRow,Temp_Vector,IError)
            If(IError.ne.0) then
              call mqc_error_i('Diagonalization error: the lapack routine DSPEV failed', 6, &
                 'IError', IError )
            endif
            DeAllocate(A_Symm)
          else
            call mqc_error_I('Diagonalization error: the matrix being passed is not square', 6, &
                 'A%NRow', A%NRow, &
                 'A%NCol', A%NCol )
          endif
        else if (A%storage .eq. 'StorSymm') then
          Allocate(Temp_Vector(3*A%NRow),A_EVecs(A%NRow,A%NRow),A_EVals(A%NRow))
          call DSPEV('V','U',A%NRow,A%matR(:,1),A_EVals,A_EVecs,A%NRow,Temp_Vector,IError)
          If(IError.ne.0) then
            call mqc_error_I('Diagonalization error: the lapack routine DSPEV failed', 6, &
                 'IError', IError )
          endif
        endif
        if(typeFlag.eq.'Integer') Call MQC_Matrix_Copy_Real2Int(A)
      elseIf(A%data_type.eq.'Complex') then

        if(.not.mqc_matrix_test_symmetric(A,'hermitian')) call mqc_error_L('Input matrix not Hermitian in MQC_Matrix_Diagonalize', &
             6, &
        "mqc_matrix_test_symmetric(A,'hermitian')", mqc_matrix_test_symmetric(A,'hermitian') )
        if(A%storage.eq.'StorFull') then
          if (A%NRow .eq. A%NCol) then
            allocate(AC_Symm((A%NRow*(A%NRow+1))/2),TempC_Vector(Max(1,2*A%NRow-1)),Temp_Vector(Max(1,3*A%NRow-2)), &
              AC_EVecs(A%NRow,A%NRow),A_EVals(A%NRow))
            k = 0
            do i = 1,A%NRow
              do j = 1,i
                k = k+1
                AC_Symm(k) = A%matC(j,i)
              endDo
            endDo
            call ZHPEV('V','U',A%NRow,AC_Symm,A_EVals,AC_EVecs,A%NRow,TempC_Vector,Temp_Vector,IError)
            if(IError.ne.0) then
              call mqc_error_i('Diagonalization error: the lapack routine ZHPEV failed', 6, &
                 'IError', IError )
            endIf
            deallocate(AC_Symm)
          else
            call mqc_error_i('Diagonalization error: the matrix being passed is not square', 6, &
                 'A%NRow', A%NRow, &
                 'A%NCol', A%NCol )
          endif
        elseIf (A%storage .eq. 'StorSymm') then
          allocate(TempC_Vector(max(1,2*A%NRow-1)),Temp_Vector(max(1,3*A%NRow-2)),AC_EVecs(A%NRow,A%NRow),A_EVals(A%NRow))
          call ZHPEV('V','U',A%NRow,A%matC(:,1),A_EVals,AC_EVecs,A%NRow,TempC_Vector,Temp_Vector,IError)
          if(IError.ne.0) then
            call mqc_error_i('Diagonalization error: the lapack routine ZHPEV failed', 6, &
                 'IError', IError )
          endIf
        endIf
      else
        call mqc_error_A('data type not recognsed in mqc_matrix_diagonalize', 6, &
             'A%data_type', A%data_type )
      endIf

      if(present(EVecs)) then
        if(allocated(A_EVecs)) EVecs = A_EVecs
        if(allocated(AC_EVecs)) EVecs = AC_EVecs
      endIf
      if(present(EVals)) EVals = A_EVals
      if(allocated(A_EVecs)) deallocate(A_EVecs)
      if(allocated(AC_EVecs)) deallocate(AC_EVecs)
      if(allocated(A_EVals)) deallocate(A_EVals)
      if(allocated(Temp_Vector)) deallocate(Temp_Vector)
      if(allocated(TempC_Vector)) deallocate(TempC_Vector)
!
      Return
      End Subroutine MQC_Matrix_Diagonalize
!
!     PROCEDURE MQC_Matrix_Cast_Real
      Function MQC_Matrix_Cast_Real(MA) Result(MB)
!
!     This function converts a mqc_matrix to real. Only integer vector is allowed
!     because it would lose accuracy otherwise.
!
!     X. Sheng, 2017.02
!
      Implicit None
      Type(MQC_Matrix),Intent(In)::MA
      Type(MQC_Matrix)::MB

      Call MQC_Allocate_Matrix(MA%NRow,MA%NCol,MB,'Real',MA%Storage)
      if(mqc_have_int(MA)) then
        MB%MatR = MA%MatI
      elseif(mqc_have_complex(MA)) then
        MB%MatR = MA%MatC
      endif

      End Function MQC_Matrix_Cast_Real
!
!
!     PROCEDURE MQC_Matrix_Cast_Complex
      Function MQC_Matrix_Cast_Complex(MA) Result(MB)
!
!     This function converts a mqc_matrix to complex.
!
!     X. Sheng, 2017.02
!
      Implicit None
      Type(MQC_Matrix),Intent(In)::MA
      Type(MQC_Matrix)::MB

      Call MQC_Allocate_Matrix(MA%NRow,MA%NCol,MB,'Complex',MA%Storage)
      if(mqc_have_int(MA)) then
        MB%MatC = MA%MatI
      elseIf(mqc_have_real(MA)) then
        MB%MatC = MA%MatR
      endIf

      End Function MQC_Matrix_Cast_Complex
!
!
!     PROCEDURE MQC_Matrix_Scalar_At
      Function MQC_Matrix_Scalar_At(Mat,I,J) Result(Element)
!
!     A function that returns the element at (I,J) of a mqc_matrix type as a mqc_scalar type.
!
!     X. Sheng 2017.02
!
      Implicit None
      Class(MQC_Matrix),Intent(In)::Mat
      Type(MQC_Scalar)::Element
      Integer(kind=int64),Intent(In)::I,J
      Integer(kind=int64)::IndI,IndJ

      IndI = I
      IndJ = J
      If (IndI.lt.0) IndI = Mat%NRow + IndI + 1
      If (IndJ.lt.0) IndJ = Mat%NRow + IndJ + 1
      If (IndI.eq.0.or.IndI.gt.Mat%NRow) Call MQC_Error_I('Index I badly specified in mqc_matrix_scalar_at', 6, &
           'IndI', IndI, &
           'Mat%NRow', Mat%NRow )
      If (IndJ.eq.0.or.IndJ.gt.Mat%NRow) Call MQC_Error_I('Index J badly specified in mqc_matrix_scalar_at', 6, &
           'IndJ', IndJ, &
           'Mat%NRow', Mat%NRow )

      If (Mat%Storage.eq.'StorFull') then
        If (Mat%Data_Type.eq.'Integer') then
          Element = (Mat%MatI(IndI,IndJ))
        Elseif (Mat%Data_Type.eq.'Real') then
          Element = (Mat%MatR(IndI,IndJ))
        Elseif (Mat%Data_Type.eq.'Complex') then
          Element = (Mat%MatC(IndI,IndJ))
        Else
          Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Scalar_At', 6, &
                 'Mat%Data_Type', Mat%Data_Type )
        EndIf
      ElseIf (Mat%Storage.eq.'StorSymm') then
        If (IndI.ge.IndJ) then
          If (Mat%Data_Type.eq.'Integer') then
            Element = (Mat%MatI(IndI*(IndI-1)/2+IndJ,1))
          ElseIf (Mat%Data_Type.eq.'Real') then
            Element = (Mat%MatR(IndI*(IndI-1)/2+IndJ,1))
          ElseIf (Mat%Data_Type.eq.'Complex') then
            Element = (Mat%MatC(IndI*(IndI-1)/2+IndJ,1))
          Else
            Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Scalar_At', 6, &
                 'Mat%Data_Type', Mat%Data_Type )
          EndIf
        Else
          If (Mat%Data_Type.eq.'Integer') then
            Element = (Mat%MatI(IndJ*(IndJ-1)/2+IndI,1))
          ElseIf (Mat%Data_Type.eq.'Real') then
            Element = (Mat%MatR(IndJ*(IndJ-1)/2+IndI,1))
          ElseIf (Mat%Data_Type.eq.'Complex') then
            Element = (Mat%MatC(IndJ*(IndJ-1)/2+IndI,1))
          Else
            Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Scalar_At', 6, &
                 'Mat%Data_Type', Mat%Data_Type )
          EndIf
        EndIf
      ElseIf (Mat%Storage.eq.'StorDiag') then
        If (IndI.eq.IndJ) then
          If (Mat%Data_Type.eq.'Integer') then
            Element = (Mat%MatI(IndI,1))
          ElseIf (Mat%Data_Type.eq.'Real') then
            Element = (Mat%MatR(IndI,1))
          ElseIf (Mat%Data_Type.eq.'Complex') then
            Element = (Mat%MatR(IndI,1))
          Else
            Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Scalar_At', 6, &
                 'Mat%Data_Type', Mat%Data_Type )
          EndIf
        Else
          If (Mat%Data_Type.eq.'Integer') then
            Element = 0
          ElseIf (Mat%Data_Type.eq.'Real') then
            Element = 0.0
          ElseIf (Mat%Data_Type.eq.'Complex') then
            Element = cmplx(0.0,0.0,kind=real64)
          Else
            Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Scalar_At', 6, &
                 'Mat%Data_Type', Mat%Data_Type )
          EndIf
        EndIf
      Else
        Call MQC_Error_A('MQC_Matrix_Scalar_At: Only full, Symm, and Diagonal type matrices are supported', 6, &
             'Mat%Storage', Mat%Storage )
      EndIf

      End Function MQC_Matrix_Scalar_At
!
!
!     PROCEDURE MQC_Matrix_Vector_At
      Function MQC_Matrix_Vector_At(Mat,Rows,Cols) Result(Vector)
!
!     A function that returns the vector between (I,J) in the Kth column/row
!     (depending on whether optional argument ColumnIn is set to 'Row' or
!     'Column', column is the default) of an MQC_Matrix Mat as an MQC_Vector
!     Vector. If I, J or K is negative, it selects the (N-I+1)th index value.
!
!     L. M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Matrix),Intent(In)::Mat
      Type(MQC_Vector)::Vector
      Integer(kind=int64),Dimension(:),Intent(In)::Rows,Cols
      Integer(kind=int64)::Length,IndI,IndJ,IndK,Cnt,L,I,J,K
      Logical::Column

      If(Size(Rows).gt.1.and.Size(Cols).gt.1) then
        Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_vector_at', 6, &
             'Size(Rows)', Size(Rows), &
             'Size(Cols)', Size(Cols) )
      ElseIf(Size(Rows).eq.1.and.Size(Cols).eq.1) then
        If(Rows(1).eq.0.and.Cols(1).ne.0) then
          Column = .True.
          I = 1
          J = Mat%NRow
          K = Cols(1)
        EndIf
        If(Cols(1).eq.0.and.Rows(1).ne.0) then
          Column = .False.
          I = 1
          J = Mat%NCol
          K = Rows(1)
        EndIf
        If(Rows(1).eq.0.and.Cols(1).eq.0) Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_vector_at', 6, &
            'Rows(1)', Rows(1), &
             'Cols(1)', Cols(1) )
        If(Rows(1).ne.0.and.Cols(1).ne.0) Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_vector_at', 6, &
             'Rows(1)', Rows(1), &
             'Cols(1)', Cols(1) )
      ElseIf(Size(Rows).eq.2.and.Size(Cols).eq.1) then
        Column = .True.
        I = Rows(1)
        J = Rows(2)
        K = Cols(1)
      ElseIf(Size(Rows).eq.1.and.Size(Cols).eq.2) then
        Column = .False.
        I = Cols(1)
        J = Cols(2)
        K = Rows(1)
      Else
        Call MQC_Error_I('Unspecified boundaries in mqc_matrix_vector_at', 6, &
             'Size(Rows)', Size(Rows), &
             'Size(Cols)', Size(Cols) )
      EndIf

      IndK = K
      If (IndK.lt.0.and.Column) IndK = Mat%NCol + IndK + 1
      If (IndK.lt.0.and..not.Column) IndK = Mat%NRow + IndK + 1
      If (IndK.le.0.or.(IndK.gt.Mat%NCol.and.Column).or.(IndK.gt.Mat%NRow.and..not.Column)) &
        Call MQC_Error_I('Cannot select the Kth Row/Column in MQC_Matrix_Vector_At', 6, &
        'IndK', IndK, &
        'Mat%NCol', Mat%NCol, &
        'Mat%NRow', Mat%NRow )

      IndI = I
      IndJ = J
      If (IndI.lt.0.and.Column) IndI = Mat%NRow + IndI + 1
      If (IndJ.lt.0.and.Column) IndJ = Mat%NRow + IndJ + 1
      If (IndI.lt.0.and..not.Column) IndI = Mat%NCol + IndI + 1
      If (IndJ.lt.0.and..not.Column) IndJ = Mat%NCol + IndJ + 1

      Length = IndJ-IndI+1
      If(Column) then
        If (Length.le.0.or.Length.gt.Mat%NRow) Call MQC_Error_I('Vector length badly &
     &    defined in MQC_Matrix_Vector_At', 6, &
     'Length', Length, &
     'Mat%NRow', Mat%NRow )
        If (IndI.le.0.or.IndI.gt.(Mat%NRow-Length+1)) Call MQC_Error_I('Index I out of bounds &
    &    in MQC_Matrix_Vector_At', 6, &
    'IndI', IndI, &
    'Mat%NRow', Mat%NRow, &
    'Length', Length )
        If (IndJ.lt.Length.or.IndJ.gt.Mat%NRow) Call MQC_Error_I('Index J out of bounds &
    &     in MQC_Matrix_Vector_At', 6, &
    'IndJ', IndJ, &
    'Length', Length, &
    'Mat%NRow', Mat%NRow )
      ElseIf(.not.Column) then
        If (Length.le.0.or.Length.gt.Mat%NCol) Call MQC_Error_I('Vector length badly &
    &     defined in MQC_Matrix_Vector_At', 6, &
    'Length', Length, &
    'Mat%NCol', Mat%NCol )
        If (IndI.le.0.or.IndI.gt.(Mat%NCol-Length+1)) Call MQC_Error_I('Index I out of bounds &
    &     in MQC_Matrix_Vector_At', 6, &
    'IndI', IndI, &
    'Mat%NCol', Mat%NCol, &
    'Length', Length )
        If (IndJ.lt.Length.or.IndJ.gt.Mat%NCol) Call MQC_Error_I('Index J out of bounds &
    &     in MQC_Matrix_Vector_At', 6, &
    'IndJ', IndJ, &
    'Length', Length, &
    'Mat%NCol', Mat%NCol )
      EndIf

      If (Mat%Storage.eq.'StorFull') then
        If (Mat%Data_Type.eq.'Integer') then
          If (Column) then
            Call MQC_Allocate_Vector(Length,Vector,'Integer')
            Vector%VecI = (Mat%MatI(IndI:IndJ,IndK))
          ElseIf (.not.Column) then
            Call MQC_Allocate_Vector(Length,Vector,'Integer')
            Vector%VecI = (Mat%MatI(IndK,IndI:IndJ))
          EndIf
        Elseif (Mat%Data_Type.eq.'Real') then
          If (Column) then
            Call MQC_Allocate_Vector(Length,Vector,'Real')
            Vector%VecR = (Mat%MatR(IndI:IndJ,IndK))
          ElseIf (.not.Column) then
            Call MQC_Allocate_Vector(Length,Vector,'Real')
            Vector%VecR = (Mat%MatR(IndK,IndI:IndJ))
          EndIf
        Elseif (Mat%Data_Type.eq.'Complex') then
          If (Column) then
            Call MQC_Allocate_Vector(Length,Vector,'Complex')
            Vector%VecC = (Mat%MatC(IndI:IndJ,IndK))
          ElseIf (.not.Column) then
            Call MQC_Allocate_Vector(Length,Vector,'Complex')
            Vector%VecC = (Mat%MatC(IndK,IndI:IndJ))
          EndIf
        Else
          Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Vector_At', 6, &
               'Mat%Data_Type', Mat%Data_Type )
        EndIf
      ElseIf (Mat%Storage.eq.'StorSymm') then
        If (Mat%Data_Type.eq.'Integer') then
          Call MQC_Allocate_Vector(Length,Vector,'Integer')
        ElseIf (Mat%Data_Type.eq.'Real') then
          Call MQC_Allocate_Vector(Length,Vector,'Real')
        ElseIf (Mat%Data_Type.eq.'Complex') then
          Call MQC_Allocate_Vector(Length,Vector,'Complex')
        Else
          Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Vector_At', 6, &
               'Mat%Data_Type', Mat%Data_Type )
        EndIf
        L = 1
        Do Cnt = IndI, IndJ
          If (Cnt.gt.IndK) then
            If (Mat%Data_Type.eq.'Integer') Vector%VecI(L) = (Mat%MatI(Cnt*(Cnt-1)/2+IndK,1))
            If (Mat%Data_Type.eq.'Real') Vector%VecR(L) = (Mat%MatR(Cnt*(Cnt-1)/2+IndK,1))
            If (Mat%Data_Type.eq.'Complex') Vector%VecC(L) = (Mat%MatC(Cnt*(Cnt-1)/2+IndK,1))
          Else
            If (Mat%Data_Type.eq.'Integer') Vector%VecI(L) = (Mat%MatI(IndK*(IndK-1)/2+Cnt,1))
            If (Mat%Data_Type.eq.'Real') Vector%VecR(L) = (Mat%MatR(IndK*(IndK-1)/2+Cnt,1))
            If (Mat%Data_Type.eq.'Complex') Vector%VecC(L) = (Mat%MatC(IndK*(IndK-1)/2+Cnt,1))
          EndIf
          L = L + 1
        EndDo
      ElseIf (Mat%Storage.eq.'StorDiag') then
        If (Mat%Data_Type.eq.'Integer') then
          Call Vector%initialize(Length,0)
        ElseIf (Mat%Data_Type.eq.'Real') then
          Call Vector%initialize(Length,0.0)
        ElseIf (Mat%Data_Type.eq.'Complex') then
          Call Vector%initialize(Length,cmplx(0.0,0.0))
        Else
          Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Vector_At', 6, &
               'Mat%Data_Type', Mat%Data_Type )
        EndIf
        L = 1
        Do Cnt = IndI,IndJ
          If (Cnt.eq.IndK) then
            If (Mat%Data_Type.eq.'Integer') Vector%VecI(L) = Mat%MatI(IndK,1)
            If (Mat%Data_Type.eq.'Real') Vector%VecR(L) = Mat%MatR(IndK,1)
            If (Mat%Data_Type.eq.'Complex') Vector%VecC(L) = Mat%MatC(IndK,1)
          EndIf
          L = L + 1
        EndDo
      Else
        Call MQC_Error_A('MQC_Matrix_Vector_At: Only full, Symm, and Diagonal type matrices are supported', 6, &
             'Mat%Storage', Mat%Storage )
      EndIf

      Vector%Column = Column

      End Function MQC_Matrix_Vector_At
!
!
!     PROCEDURE MQC_Matrix_Vector_Put
!      Recursive Subroutine MQC_Matrix_Vector_Put(Mat,I,J,K,VectorIn,ColumnIn)
      Recursive Subroutine MQC_Matrix_Vector_Put(Mat,VectorIn,Rows,Cols)
!
!     This subroutine updates the value of MQC_Matrix Mat with the value
!     in MQC_Vector VectorIn elements between I and J in the Kth column/row
!     (depending on whether optional argument ColumnIn is set to 'Row' or
!     'Column', column is the default). If I, J or K are negative, it selects
!     the (N-I+1)th index value.
!
!     L. M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Matrix),Intent(InOut)::Mat
      Integer(kind=int64),Dimension(:),Intent(In)::Rows,Cols
      Type(MQC_Vector),Intent(In)::VectorIn
      Integer(kind=int64)::Length,IndI,IndJ,IndK,I,J,K
      Logical::Column

 1020 Format( " " )
      If(Size(Rows).gt.1.and.Size(Cols).gt.1) then
        Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_vector_put', 6, &
             'Size(Rows)', Size(Rows), &
             'Size(Cols)', Size(Cols) )
      ElseIf(Size(Rows).eq.1.and.Size(Cols).eq.1) then
        If(Rows(1).eq.0.and.Cols(1).ne.0) then
          Column = .True.
          I = 1
          J = Mat%NRow
          K = Cols(1)
        EndIf
        If(Cols(1).eq.0.and.Rows(1).ne.0) then
          Column = .False.
          I = 1
          J = Mat%NCol
          K = Rows(1)
        EndIf
        If(Rows(1).eq.0.and.Cols(1).eq.0) Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_vector_put', 6, &
             'Rows(1)', Rows(1), &
             'Cols(1)', Cols(1) )
        If(Rows(1).ne.0.and.Cols(1).ne.0) Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_vector_put', 6, &
             'Rows(1)', Rows(1), &
             'Cols(1)', Cols(1) )
      ElseIf(Size(Rows).eq.2.and.Size(Cols).eq.1) then
        Column = .True.
        I = Rows(1)
        J = Rows(2)
        K = Cols(1)
      ElseIf(Size(Rows).eq.1.and.Size(Cols).eq.2) then
        Column = .False.
        I = Cols(1)
        J = Cols(2)
        K = Rows(1)
      Else
        Call MQC_Error_I('Unspecified boundaries in mqc_matrix_vector_put', 6, &
             'Size(Rows)', Size(Rows), &
             'Size(Cols)', Size(Cols) )
      EndIf

      IndK = K
      If (IndK.lt.0.and.Column) IndK = Mat%NCol + IndK + 1
      If (IndK.lt.0.and..not.Column) IndK = Mat%NRow + IndK + 1
      If (IndK.le.0.or.(IndK.gt.Mat%NCol.and.Column).or.(IndK.gt.Mat%NRow.and..not.Column)) &
        Call MQC_Error_I('Cannot select the Kth Row/Column in MQC_Matrix_Vector_Put', 6, &
        'IndK', IndK, &
        'Mat%NCol', Mat%NCol, &
        'Mat%NRow', Mat%NRow )
      IndI = I
      IndJ = J
      If (IndI.lt.0.and.Column) IndI = Mat%NRow + IndI + 1
      If (IndJ.lt.0.and.Column) IndJ = Mat%NRow + IndJ + 1
      If (IndI.lt.0.and..not.Column) IndI = Mat%NCol + IndI + 1
      If (IndJ.lt.0.and..not.Column) IndJ = Mat%NCol + IndJ + 1

      Length = IndJ-IndI+1
      If(Column) then
        If (Length.le.0.or.Length.gt.Mat%NRow.or.Length.ne.VectorIn%Length) &
          Call MQC_Error_I('Vector length badly defined in MQC_Matrix_Vector_Put', 6, &
          'Length', Length, &
          'Mat%NRow', Mat%NRow, &
          'VectorIn%Length', VectorIn%Length )
        If (IndI.le.0.or.IndI.gt.(Mat%NRow-Length+1)) Call MQC_Error_I('Index I out of bounds &
    &    in MQC_Matrix_Vector_Put', 6, &
    'IndI', IndI, &
    'Mat%NRow', Mat%NRow, &
    'Length', Length )
        If (IndJ.lt.Length.or.IndJ.gt.Mat%NRow) Call MQC_Error_I('Index J out of bounds &
    &     in MQC_Matrix_Vector_Put', 6, &
    'IndJ', IndJ, &
    'Length', Length, &
    'Mat%NRow', Mat%NRow )
      ElseIf(.not.Column) then
        If (Length.le.0.or.Length.gt.Mat%NCol.or.Length.ne.VectorIn%Length) &
          Call MQC_Error_I('Vector length badly defined in MQC_Matrix_Vector_Put', 6, &
          'Length', Length, &
          'Mat%NCol', Mat%NCol, &
          'VectorIn%Length', VectorIn%Length )
        If (IndI.le.0.or.IndI.gt.(Mat%NCol-Length+1)) Call MQC_Error_I('Index I out of bounds &
    &     in MQC_Matrix_Vector_Put', 6, &
    'IndI', IndI, &
    'Mat%NCol', Mat%NCol, &
    'Length', Length )
        If (IndJ.lt.Length.or.IndJ.gt.Mat%NCol) Call MQC_Error_I('Index J out of bounds &
             &     in MQC_Matrix_Vector_Put', 6, &
             'IndJ', IndJ, &
             'Length', Length, &
             'Mat%NCol', Mat%NCol )
      EndIf

      If (Mat%Storage.eq.'StorFull') then
        If (Mat%Data_Type.eq.'Integer') then
          If (VectorIn%Data_Type.eq.'Integer') then
            If(Column) Mat%MatI(IndI:IndJ,IndK) = VectorIn%VecI
            If(.not.Column) Mat%MatI(IndK,IndI:IndJ) = VectorIn%VecI
          ElseIf (VectorIn%Data_Type.eq.'Real') then
            Call MQC_Matrix_Copy_Int2Real(Mat)
            If(Column) Mat%MatR(IndI:IndJ,IndK) = VectorIn%VecR
            If(.not.Column) Mat%MatR(IndK,IndI:IndJ) = VectorIn%VecR
          ElseIf (VectorIn%Data_Type.eq.'Complex') then
            Call MQC_Matrix_Copy_Int2Complex(Mat)
            If(Column) Mat%MatC(IndI:IndJ,IndK) = VectorIn%VecC
            If(.not.Column) Mat%MatC(IndK,IndI:IndJ) = VectorIn%VecC
          Else
            Call MQC_Error_A('VectorIn type not defined in MQC_Vector_Vector_Put', 6, &
                 'VectorIn%Data_Type', VectorIn%Data_Type )
          EndIf
        ElseIf (Mat%Data_Type.eq.'Real') then
          If (VectorIn%Data_Type.eq.'Integer') then
            If(Column) Mat%MatR(IndI:IndJ,IndK) = VectorIn%VecI
            If(.not.Column) Mat%MatR(IndK,IndI:IndJ) = VectorIn%VecI
          ElseIf (VectorIn%Data_Type.eq.'Real') then
            If(Column) Mat%MatR(IndI:IndJ,IndK) = VectorIn%VecR
            If(.not.Column) Mat%MatR(IndK,IndI:IndJ) = VectorIn%VecR
          ElseIf (VectorIn%Data_Type.eq.'Complex') then
            Call MQC_Matrix_Copy_Real2Complex(Mat)
            If(Column) Mat%MatC(IndI:IndJ,IndK) = VectorIn%VecC
            If(.not.Column) Mat%MatC(IndK,IndI:IndJ) = VectorIn%VecC
          Else
            Call MQC_Error_A('VectorIn type not defined in MQC_Vector_Vector_Put', 6, &
                 'VectorIn%Data_Type', VectorIn%Data_Type )
          EndIf
        ElseIf (Mat%Data_Type.eq.'Complex') then
          If (VectorIn%Data_Type.eq.'Integer') then
            If(Column) Mat%MatC(IndI:IndJ,IndK) = VectorIn%VecI
            If(.not.Column) Mat%MatC(IndK,IndI:IndJ) = VectorIn%VecI
          ElseIf (VectorIn%Data_Type.eq.'Real') then
            If(Column) Mat%MatC(IndI:IndJ,IndK) = VectorIn%VecR
            If(.not.Column) Mat%MatC(IndK,IndI:IndJ) = VectorIn%VecR
          ElseIf (VectorIn%Data_Type.eq.'Complex') then
            If(Column) Mat%MatC(IndI:IndJ,IndK) = VectorIn%VecC
            If(.not.Column) Mat%MatC(IndK,IndI:IndJ) = VectorIn%VecC
          Else
            Call MQC_Error_A('VectorIn type not defined in MQC_Vector_Vector_Put', 6, &
                 'VectorIn%Data_Type', VectorIn%Data_Type )
          EndIf
        Else
          Call MQC_Error_A('Mat type not defined in MQC_Vector_Vector_Put', 6, &
               'Mat%Data_Type', Mat%Data_Type )
        EndIf
      ElseIf (Mat%Storage.eq.'StorSymm') then
        Call MQC_Matrix_Symm2Full(Mat)
        If(Column) Call MQC_Matrix_Vector_Put(Mat,VectorIn,[IndI,IndJ],[IndK])
        If(.not.Column) Call MQC_Matrix_Vector_Put(Mat,VectorIn,[IndK],[IndI,IndJ])
        Write(*,1020)
        If(MQC_Matrix_Test_Symmetric(Mat)) Call MQC_Matrix_Full2Symm(Mat)
      ElseIf (Mat%Storage.eq.'StorDiag') then
        Call MQC_Matrix_Diag2Full(Mat)
        If(Column) Call MQC_Matrix_Vector_Put(Mat,VectorIn,[IndI,IndJ],[IndK])
        If(.not.Column) Call MQC_Matrix_Vector_Put(Mat,VectorIn,[IndK],[IndI,IndJ])
        If(MQC_Matrix_Test_Diagonal(Mat)) Call MQC_Matrix_Full2Diag(Mat)
      Else
        Call MQC_Error_A('Matrix type not recognosed in MQC_Matrix_Vector_Put', 6, &
             'Mat%Storage', Mat%Storage )
      EndIf

      End Subroutine MQC_Matrix_Vector_Put
!
!
!     PROCEDURE MQC_Matrix_Matrix_At
      Function MQC_Matrix_Matrix_At(Mat,Rows,Cols) Result(Matrix)
!
!     A function that returns the matrix between rows (I,J) and columns (K,L)
!     of an MQC_Matrix Mat as an MQC_Matrix Matrix. If I, J, K or L is negative,
!     it selects the (N-I+1)th index value.
!
!     L. M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Matrix),Intent(In)::Mat
      Type(MQC_Matrix)::Matrix
      Integer(kind=int64),Dimension(:),Intent(In)::Rows,Cols
      Integer(kind=int64)::LenCol,LenRow,IndI,IndJ,IndK,IndL
      Integer(kind=int64)::M,N,Cnt1,Cnt2,I,J,K,L

      If(Size(Rows).gt.2.and.Size(Cols).gt.2) then
        Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_matrix_at', 6, &
             'Size(Rows)', Size(Rows), &
             'Size(Cols)', Size(Cols) )
      ElseIf(Size(Rows).eq.2.and.Size(Cols).eq.2) then
          I = Rows(1)
          J = Rows(2)
          K = Cols(1)
          L = Cols(2)
      ElseIf(Size(Rows).eq.1.and.Size(Cols).eq.1) then
        If(Rows(1).eq.0.and.Cols(1).eq.0) then
          I = 1
          J = Mat%NRow
          K = 1
          L = Mat%NCol
        Else
          Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_matrix_at', 6, &
               'Rows(1)', Rows(1), &
               'Cols(1)', Cols(1), &
               'Size(Rows)', Size(Rows) , &
               'Size(Cols)', Size(Cols) )
        EndIf
      ElseIf(Size(Rows).eq.2.and.Size(Cols).eq.1) then
        If(Cols(1).eq.0) then
          I = Rows(1)
          J = Rows(2)
          K = 1
          L = Mat%NCol
        Else
          Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_matrix_at', 6, &
               'Cols(1)', Cols(1), &
               'Size(Rows)', Size(Rows) , &
               'Size(Cols)', Size(Cols) )
        EndIf
      ElseIf(Size(Rows).eq.1.and.Size(Cols).eq.2) then
        If(Rows(1).eq.0) then
          I = 1
          J = Mat%NRow
          K = Cols(1)
          L = Cols(2)
        Else
          Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_matrix_at', 6, &
               'Rows(1)', Rows(1), &
               'Size(Rows)', Size(Rows) , &
               'Size(Cols)', Size(Cols) )
        EndIf
      Else
        Call MQC_Error_I('Unspecified boundaries in mqc_matrix_matrix_at', 6, &
             'Size(Rows)', Size(Rows) , &
             'Size(Cols)', Size(Cols) )
      EndIf

      IndI = I
      IndJ = J
      If (IndI.lt.0) IndI = Mat%NRow + IndI + 1
      If (IndJ.lt.0) IndJ = Mat%NRow + IndJ + 1
      LenRow = IndJ-IndI+1
      If (LenRow.le.0.or.LenRow.gt.Mat%NRow) Call MQC_Error_I('Row length badly &
           &   defined in MQC_Matrix_Matrix_At', 6, &
           'LenRow', LenRow, &
           'Mat%NRow', Mat%NRow )
      If (IndI.le.0.or.IndI.gt.(Mat%NRow-LenRow+1)) Call MQC_Error_I('Index I out of bounds &
           &   in MQC_Matrix_Matrix_At', 6, &
           'IndI', IndI, &
           'Mat%NRow', Mat%NRow, &
           'LenRow', LenRow )
      If (IndJ.lt.LenRow.or.IndJ.gt.Mat%NRow) Call MQC_Error_I('Index J out of bounds &
           &   in MQC_Matrix_Matrix_At', 6, &
           'IndJ', IndJ, &
           'LenRow', LenRow, &
           'Mat%NRow', Mat%NRow )
      IndK = K
      IndL = L
      If (IndK.lt.0) IndK = Mat%NCol + IndK + 1
      If (IndL.lt.0) IndL = Mat%NCol + IndL + 1
      LenCol = IndL-IndK+1
      If (LenCol.le.0.or.LenCol.gt.Mat%NCol) Call MQC_Error_I('Column length badly &
           &  defined in MQC_Matrix_Matrix_At', 6, &
           'LenCol', LenCol, &
           'Mat%NCol', Mat%NCol )
      If (IndK.le.0.or.IndK.gt.(Mat%NCol-LenCol+1)) Call MQC_Error_I('Index K out of bounds &
           &   in MQC_Matrix_Matrix_At', 6, &
           'IndK', IndK, &
           'Mat%NCol', Mat%NCol, &
           'LenCol', LenCol )
      If (IndL.lt.LenCol.or.IndL.gt.Mat%NCol) Call MQC_Error_I('Index L out of bounds &
           &   in MQC_Matrix_Matrix_At', 6, &
           'IndL', IndL, &
           'LenCol', LenCol, &
           'Mat%NCol', Mat%NCol )
      If (Mat%Storage.eq.'StorFull') then
        If (Mat%Data_Type.eq.'Integer') then
          Call MQC_Allocate_Matrix(LenRow,LenCol,Matrix,'Integer','StorFull')
          Matrix%MatI = (Mat%MatI(IndI:IndJ,IndK:IndL))
        Elseif (Mat%Data_Type.eq.'Real') then
          Call MQC_Allocate_Matrix(LenRow,LenCol,Matrix,'Real','StorFull')
          Matrix%MatR = (Mat%MatR(IndI:IndJ,IndK:IndL))
        Elseif (Mat%Data_Type.eq.'Complex') then
          Call MQC_Allocate_Matrix(LenRow,LenCol,Matrix,'Complex','StorFull')
          Matrix%MatC = (Mat%MatC(IndI:IndJ,IndK:IndL))
        Else
          Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Vector_At', 6, &
               'Mat%Data_Type', Mat%Data_Type )
        EndIf
      ElseIf (Mat%Storage.eq.'StorSymm') then
        If (Mat%Data_Type.eq.'Integer') then
          Call MQC_Allocate_Matrix(LenRow,LenCol,Matrix,'Integer','StorFull')
        ElseIf (Mat%Data_Type.eq.'Real') then
          Call MQC_Allocate_Matrix(LenRow,LenCol,Matrix,'Real','StorFull')
        ElseIf (Mat%Data_Type.eq.'Complex') then
          Call MQC_Allocate_Matrix(LenRow,LenCol,Matrix,'Complex','StorFull')
        EndIf
        M = 1
        N = 1
        Do Cnt1 = IndI, IndJ
          Do Cnt2 = IndK, IndL
            If (Cnt1.gt.Cnt2) then
              If (Mat%Data_Type.eq.'Integer') Matrix%MatI(M,N) = (Mat%MatI(Cnt1*(Cnt1-1)/2+Cnt2,1))
              If (Mat%Data_Type.eq.'Real') Matrix%MatR(M,N) = (Mat%MatR(Cnt1*(Cnt1-1)/2+Cnt2,1))
              If (Mat%Data_Type.eq.'Complex') Matrix%MatC(M,N) = (Mat%MatC(Cnt1*(Cnt1-1)/2+Cnt2,1))
            Else
              If (Mat%Data_Type.eq.'Integer') Matrix%MatI(M,N) = (Mat%MatI(Cnt2*(Cnt2-1)/2+Cnt1,1))
              If (Mat%Data_Type.eq.'Real') Matrix%MatR(M,N) = (Mat%MatR(Cnt2*(Cnt2-1)/2+Cnt1,1))
              If (Mat%Data_Type.eq.'Complex') Matrix%MatC(M,N) = (Mat%MatC(Cnt2*(Cnt2-1)/2+Cnt1,1))
            EndIf
            If(N.eq.LenCol) then
              N = 1
            Else
              N = N + 1
            EndIf
          EndDo
          M = M + 1
        EndDo
        If(MQC_Matrix_Test_Symmetric(Matrix)) Call MQC_Matrix_Full2Symm(Matrix)
      ElseIf (Mat%Storage.eq.'StorDiag') then
        If (Mat%Data_Type.eq.'Integer') then
          Call Matrix%initialize(LenRow,LenCol,0)
        ElseIf (Mat%Data_Type.eq.'Real') then
          Call Matrix%initialize(LenRow,LenCol,0.0)
        ElseIf (Mat%Data_Type.eq.'Complex') then
          Call Matrix%initialize(LenRow,LenCol,cmplx(0.0,0.0))
        EndIf
        M = 1
        N = 1
        Do Cnt1 = IndI, IndJ
          Do Cnt2 = IndK, IndL
            If(Cnt1.eq.Cnt2) then
              If (Mat%Data_Type.eq.'Integer') Matrix%MatI(M,N) = Mat%MatI(Cnt1,1)
              If (Mat%data_Type.eq.'Real') Matrix%MatR(M,N) = Mat%MatR(Cnt1,1)
              If (Mat%data_Type.eq.'Complex') Matrix%MatR(M,N) = Mat%MatR(Cnt1,1)
            EndIf
            If(N.eq.LenCol) then
              N = 1
            Else
              N = N + 1
            EndIf
          EndDo
          M = M + 1
        EndDo
        If(MQC_Matrix_Test_Diagonal(Matrix)) Call MQC_Matrix_Full2Diag(Matrix)
      Else
        Call MQC_Error_A('Matrix type not recognised in MQC_Matrix_Matrix_At', 6, &
             'Mat%Storage', Mat%Storage )
      EndIf

      End Function MQC_Matrix_Matrix_At
!
!
!     PROCEDURE MQC_Matrix_DiagMatrix_Put_Vector 
      Subroutine MQC_Matrix_DiagMatrix_Put_Vector(diagVectorIn,mat)
!
!     This subroutine loads a (MQC) diagonal matrix <mat> from an input (mqc vector
!     type) diagonal vector <diagVectorIn>.
!
!     L. M. Thompson, 2018.
!
!
      implicit none
      class(MQC_Matrix),intent(inOut)::mat
      class(MQC_Vector),intent(in)::diagVectorIn
!
      integer(kind=int64)::n
      character(len=64)::data_type 
!
      n = diagVectorIn%length
      data_type = diagVectorIn%data_type
      call MQC_Allocate_Matrix(n,n,mat,data_type,'stordiag')
      select case (data_type)
      case ('Integer') 
        mat%matI(:,1) = diagVectorIn%vecI(:)
      case ('Real') 
        mat%matR(:,1) = diagVectorIn%vecR(:)
      case ('Complex') 
        mat%matC(:,1) = diagVectorIn%vecC(:)
      case default
        call mqc_error_A('diagVectorIn has unrecognised data type:',6,'data_type',data_type)
      end select
!
      return
      end subroutine MQC_Matrix_DiagMatrix_Put_Vector 
!
!
!     PROCEDURE MQC_Matrix_DiagMatrix_Put_Integer
      Subroutine MQC_Matrix_DiagMatrix_Put_Integer(mat,diagMatrixIn)
!
!     This subroutine loads a (MQC) diagonal matrix <mat> from an input (fortran
!     intrinsic type) integer diagonal matrix <diagMatrixIn>, which should be
!     passed as a rank-1 array with size n where n is the leading dimension of
!     the matrix.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      class(MQC_Matrix),intent(inOut)::mat
      integer(kind=int64),dimension(:),intent(in)::diagMatrixIn
!
      integer(kind=int64)::n
!
!
!     Begin by figuring out the leading dimension of the matrix. Then, ensure
!     mat is allocated correctly and then use mat%mput routine.
!
      n = SIZE(diagMatrixIn)
      call MQC_Allocate_Matrix(n,n,mat,'integer','stordiag')
      mat%matI(:,1) = diagMatrixIn(:)
!
      return
      end subroutine MQC_Matrix_DiagMatrix_Put_Integer
!
!
!     PROCEDURE MQC_Matrix_DiagMatrix_Put_Real
      Subroutine MQC_Matrix_DiagMatrix_Put_Real(mat,diagMatrixIn)
!
!     This subroutine loads a (MQC) diagonal matrix <mat> from an input (fortran
!     intrinsic type) integer diagonal matrix <diagMatrixIn>, which should be
!     passed as a rank-1 array with size n where n is the leading dimension of
!     the matrix.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      class(MQC_Matrix),intent(inOut)::mat
      real(kind=real64),dimension(:),intent(in)::diagMatrixIn
!
      integer(kind=int64)::n
!
!
!     Begin by figuring out the leading dimension of the matrix. Then, ensure
!     mat is allocated correctly and then use mat%mput routine.
!
      n = SIZE(diagMatrixIn)
      call MQC_Allocate_Matrix(n,n,mat,'real','stordiag')
      mat%matR(:,1) = diagMatrixIn(:)
!
      return
      end subroutine MQC_Matrix_DiagMatrix_Put_Real
!
!
!     PROCEDURE MQC_Matrix_DiagMatrix_Put_Complex
      Subroutine MQC_Matrix_DiagMatrix_Put_Complex(mat,diagMatrixIn)
!
!     This subroutine loads a (MQC) diagonal matrix <mat> from an input (fortran
!     intrinsic type) integer diagonal matrix <diagMatrixIn>, which should be
!     passed as a rank-1 array with size n where n is the leading dimension of
!     the matrix.
!
!     L. M. Thompson, 2017.
!
!
      implicit none
      class(MQC_Matrix),intent(inOut)::mat
      Complex(Kind=real64),dimension(:),intent(in)::diagMatrixIn
!
      integer(kind=int64)::n
!
!
!     Begin by figuring out the leading dimension of the matrix. Then, ensure
!     mat is allocated correctly and then use mat%mput routine.
!
      n = SIZE(diagMatrixIn)
      call MQC_Allocate_Matrix(n,n,mat,'Complex','stordiag')
      mat%matC(:,1) = diagMatrixIn(:)
!
      return
      end subroutine MQC_Matrix_DiagMatrix_Put_Complex
!
!
!     PROCEDURE MQC_Matrix_SymmMatrix_Put_Integer
      Subroutine MQC_Matrix_SymmMatrix_Put_Integer(mat,symmMatrixIn)
!
!     This subroutine loads a (MQC) symmetric matrix <mat> from an input
!     (fortran intrinsic type) integer packed-symmetric matrix <symmMatrixIn>,
!     which should be passed as a rank-1 array with size n(n+1)/2 where n is the
!     leading dimension of the matrix.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      class(MQC_Matrix),intent(inOut)::mat
      integer(kind=int64),dimension(:),intent(in)::symmMatrixIn
!
      integer(kind=int64)::n,nSymm,intTemp1
      real(kind=real64)::realTemp1,realTemp2
!
!
!     Begin by figuring out the leading dimension of the matrix. Then, ensure
!     mat is allocated correctly and then use mat%mput routine.
!
      nSymm = SIZE(symmMatrixIn)
      intTemp1 = 1+8*nSymm
#ifdef PGI
      realTemp1 = (sqrt(dfloat(intTemp1))-1)/dfloat(2)
#else
      realTemp1 = (sqrt(float(intTemp1))-1)/float(2)
#endif
      n = INT(realTemp1)
#ifdef PGI
      realTemp2 = ABS(realTemp1 - dfloat(n))
#else
      realTemp2 = ABS(realTemp1 - float(n))
#endif
      if(realTemp2.gt.1.d-5) call MQC_Error_R('Error in MQC_Matrix_SymmMatrix_Put', 6, &
           'realTemp2', realTemp2 )
      call MQC_Allocate_Matrix(n,n,mat,'integer','StorSymm')
      mat%matI(:,1) = symmMatrixIn(:)
!
      return
      end subroutine MQC_Matrix_SymmMatrix_Put_Integer
!
!
!     PROCEDURE MQC_Matrix_SymmMatrix_Put_Real
      Subroutine MQC_Matrix_SymmMatrix_Put_Real(mat,symmMatrixIn)
!
!     This subroutine loads a (MQC) symmetric matrix <mat> from an input
!     (fortran intrinsic type) real packed-symmetric matrix <symmMatrixIn>,
!     which should be passed as a rank-1 array with size n(n+1)/2 where n is the
!     leading dimension of the matrix.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      class(MQC_Matrix),intent(inOut)::mat
      real(kind=real64),dimension(:),intent(in)::symmMatrixIn
!
      integer(kind=int64)::n,nSymm,intTemp1
      real(kind=real64)::realTemp1,realTemp2
!
!
!     Begin by figuring out the leading dimension of the matrix. Then, ensure
!     mat is allocated correctly and then use mat%mput routine.
!
      nSymm = SIZE(symmMatrixIn)
      intTemp1 = 1+8*nSymm
#ifdef PGI
      realTemp1 = (sqrt(dfloat(intTemp1))-1)/dfloat(2)
#else
      realTemp1 = (sqrt(float(intTemp1))-1)/float(2)
#endif
      n = INT(realTemp1)
#ifdef PGI
      realTemp2 = ABS(realTemp1 - dfloat(n))
#else
      realTemp2 = ABS(realTemp1 - float(n))
#endif
      if(realTemp2.gt.1.d-5) call MQC_Error_R('Error in MQC_Matrix_SymmMatrix_Put', 6, &
           'realTemp2', realTemp2 )

      call MQC_Allocate_Matrix(n,n,mat,'real','symm')
      mat%matR(:,1) = symmMatrixIn(:)
!
      return
      end subroutine MQC_Matrix_SymmMatrix_Put_Real
!
!
!     PROCEDURE MQC_Matrix_SymmMatrix_Put_Complex
      Subroutine MQC_Matrix_SymmMatrix_Put_Complex(mat,symmMatrixIn)
!
!     This subroutine loads a (MQC) symmetric matrix <mat> from an input
!     (fortran intrinsic type) complex packed-symmetric matrix <symmMatrixIn>,
!     which should be passed as a rank-1 array with size n(n+1)/2 where n is the
!     leading dimension of the matrix.
!
!     H. P. Hratchian, 2017.
!     L. M. Thompson, 2017.
!
!
      implicit none
      class(MQC_Matrix),intent(inOut)::mat
      Complex(Kind=real64),dimension(:),intent(in)::symmMatrixIn
!
      integer(kind=int64)::n,nSymm,intTemp1
      real(kind=real64)::realTemp1,realTemp2
!
!
!     Begin by figuring out the leading dimension of the matrix. Then, ensure
!     mat is allocated correctly and then use mat%mput routine.
!
      nSymm = SIZE(symmMatrixIn)
      intTemp1 = 1+8*nSymm
#ifdef PGI
      realTemp1 = (sqrt(dfloat(intTemp1))-1)/dfloat(2)
#else
      realTemp1 = (sqrt(float(intTemp1))-1)/float(2)
#endif
      n = INT(realTemp1)
#ifdef PGI
      realTemp2 = ABS(realTemp1 - dfloat(n))
#else
      realTemp2 = ABS(realTemp1 - float(n))
#endif
      if(realTemp2.gt.1.d-5) call MQC_Error_R('Error in MQC_Matrix_SymmMatrix_Put', 6, &
           'realTemp2', realTemp2 )
      call MQC_Allocate_Matrix(n,n,mat,'Complex','symm')
      mat%matC(:,1) = symmMatrixIn(:)
!
      return
      end subroutine MQC_Matrix_SymmMatrix_Put_Complex
!
!
!     PROCEDURE MQC_Matrix_Matrix_Put
      Recursive Subroutine MQC_Matrix_Matrix_Put(Mat,MatrixIn,Rows,Cols)
!
!     This subroutine updates the value of MQC_Matrix Mat with the value
!     in MQC_Matrix MatrixIn elements between rows I and J and columns K
!     and L. If I, J, K or L are negative, it selects the (N-I+1)th index
!     value.
!
!     L. M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Matrix),Intent(InOut)::Mat
      Type(MQC_Matrix),Intent(In)::MatrixIn
      Integer(kind=int64),Dimension(:),Intent(In)::Rows,Cols
      Integer(kind=int64)::IndI,IndJ,IndK,IndL,LenCol,LenRow
      Integer(kind=int64)::M,N,Cnt1,Cnt2,I,J,K,L
!
!
!     Do the work...
!
 1020 Format( " " )
      If(Size(Rows).gt.2.and.Size(Cols).gt.2) then
        Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_matrix_put', 6, &
             'Size(Rows)', Size(Rows), &
             'Size(Cols)', Size(Cols) )
      ElseIf(Size(Rows).eq.2.and.Size(Cols).eq.2) then
          I = Rows(1)
          J = Rows(2)
          K = Cols(1)
          L = Cols(2)
      ElseIf(Size(Rows).eq.1.and.Size(Cols).eq.1) then
        If(Rows(1).eq.0.and.Cols(1).eq.0) then
          I = 1
          J = Mat%NRow
          K = 1
          L = Mat%NCol
        Else
          Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_matrix_put', 6, &
               'Rows(1)', Rows(1), &
               'Cols(1)', Cols(1) )
        EndIf
      ElseIf(Size(Rows).eq.2.and.Size(Cols).eq.1) then
        If(Cols(1).eq.0) then
          I = Rows(1)
          J = Rows(2)
          K = 1
          L = Mat%NCol
        Else
          Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_matrix_put', 6, &
               'Cols(1)', Cols(1) )
        EndIf
      ElseIf(Size(Rows).eq.1.and.Size(Cols).eq.2) then
        If(Rows(1).eq.0) then
          I = 1
          J = Mat%NRow
          K = Cols(1)
          L = Cols(2)
        Else
          Call MQC_Error_I('Vector bounds badly specified in mqc_matrix_matrix_put', 6, &
               'Rows(1)', Rows(1) )
        EndIf
      Else
        Call MQC_Error_I('Unspecified boundaries in mqc_matrix_matrix_put', 6, &
             'Size(Rows)', Size(Rows), &
             'Size(Cols)', Size(Cols) )
      EndIf
!
      IndI = I
      IndJ = J
      IndK = K
      IndL = L
      If (IndI.lt.0) IndI = Mat%NRow + IndI + 1
      If (IndJ.lt.0) IndJ = Mat%NRow + IndJ + 1
      If (IndK.lt.0) IndK = Mat%NCol + IndK + 1
      If (IndL.lt.0) IndL = Mat%NCol + IndL + 1
      LenRow = IndJ-IndI+1
      LenCol = IndL-IndK+1
      If (LenRow.le.0.or.LenRow.gt.Mat%NRow.or.LenRow.ne.MatrixIn%NRow) Call MQC_Error_I('Row length badly &
     &  defined in MQC_Matrix_Matrix_At', 6, &
     'LenRow', LenRow, &
     'Mat%NRow', Mat%NRow, &
     'MatrixIn%NRow', MatrixIn%NRow )
      If (IndI.le.0.or.IndI.gt.(Mat%NRow-LenRow+1)) Call MQC_Error_I('Index I out of bounds &
     &  in MQC_Matrix_Matrix_At', 6, &
     'IndI', IndI, &
     'Mat%NRow', Mat%NRow, &
     'LenRow', LenRow )
      If (IndJ.lt.LenRow.or.IndJ.gt.Mat%NRow) Call MQC_Error_I('Index J out of bounds &
     &  in MQC_Matrix_Matrix_At', 6, &
     'IndJ', IndJ, &
     'LenRow', LenRow, &
     'Mat%NRow', Mat%NRow )
      If (LenCol.le.0.or.LenCol.gt.Mat%NCol.or.LenCol.ne.MatrixIn%NCol) Call MQC_Error_I('Column length badly &
     &  defined in MQC_Matrix_Matrix_At', 6, &
     'LenCol', LenCol, &
     'Mat%NCol', Mat%NCol, &
     'MatrixIn%NCol', MatrixIn%NCol )
      If (IndK.le.0.or.IndK.gt.(Mat%NCol-LenCol+1)) Call MQC_Error_I('Index K out of bounds &
     &  in MQC_Matrix_Matrix_At', 6, &
     'IndK', IndK, &
     'Mat%NCol', Mat%NCol, &
     'LenCol', LenCol )
      If (IndL.lt.LenCol.or.IndL.gt.Mat%NCol) Call MQC_Error_I('Index L out of bounds &
     &  in MQC_Matrix_Matrix_At', 6, &
     'IndL', IndL, &
     'LenCol', LenCol, &
     'Mat%NCol', Mat%NCol )

      If (Mat%Storage.eq.'StorFull') then
        If (Mat%Data_Type.eq.'Integer') then
          If (MatrixIn%Data_Type.eq.'Integer') then
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatI(IndI:IndJ,IndK:IndL) = MatrixIn%MatI
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatI(IndI+N-1,IndK+M-1) = MatrixIn%MatI(L,1)
                  Mat%MatI(IndK+M-1,IndI+N-1) = MatrixIn%MatI(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatI(Cnt1,Cnt2) = (MatrixIn%MatI(M,1))
                  Else
                    Mat%MatI(Cnt1,Cnt2) = 0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage )
            EndIf
          ElseIf (MatrixIn%Data_Type.eq.'Real') then
            Call MQC_Matrix_Copy_Int2Real(Mat)
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatR(IndI:IndJ,IndK:IndL) = MatrixIn%MatR
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatI(IndI+N-1,IndK+M-1) = MatrixIn%MatR(L,1)
                  Mat%MatI(IndK+M-1,IndI+N-1) = MatrixIn%MatR(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatR(Cnt1,Cnt2) = (MatrixIn%MatR(M,1))
                  Else
                    Mat%MatR(Cnt1,Cnt2) = 0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage )
            EndIf
          ElseIf (MatrixIn%Data_Type.eq.'Complex') then
            Call MQC_Matrix_Copy_Int2Complex(Mat)
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatC(IndI:IndJ,IndK:IndL) = MatrixIn%MatC
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatI(IndI+N-1,IndK+M-1) = MatrixIn%MatC(L,1)
                  Mat%MatI(IndK+M-1,IndI+N-1) = MatrixIn%MatC(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatC(Cnt1,Cnt2) = (MatrixIn%MatC(M,1))
                  Else
                    Mat%MatC(Cnt1,Cnt2) = 0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage )
            EndIf
          Else
            Call MQC_Error_A('MatrixIn type not defined in MQC_Vector_Vector_Put', 6, &
                 'MatrixIn%Data_Type', MatrixIn%Data_Type )
          EndIf
        ElseIf (Mat%Data_Type.eq.'Real') then
          If (MatrixIn%Data_Type.eq.'Integer') then
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatR(IndI:IndJ,IndK:IndL) = MatrixIn%MatI
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatR(IndI+N-1,IndK+M-1) = MatrixIn%MatI(L,1)
                  Mat%MatR(IndK+M-1,IndI+N-1) = MatrixIn%MatI(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatR(Cnt1,Cnt2) = (MatrixIn%MatI(M,1))
                  Else
                    Mat%MatR(Cnt1,Cnt2) = 0.0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage )
            EndIf
          ElseIf (MatrixIn%Data_Type.eq.'Real') then
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatR(IndI:IndJ,IndK:IndL) = MatrixIn%MatR
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatR(IndI+N-1,IndK+M-1) = MatrixIn%MatR(L,1)
                  Mat%MatR(IndK+M-1,IndI+N-1) = MatrixIn%MatR(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatR(Cnt1,Cnt2) = (MatrixIn%MatR(M,1))
                  Else
                    Mat%MatR(Cnt1,Cnt2) = 0.0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage )
            EndIf
          ElseIf (MatrixIn%Data_Type.eq.'Complex') then
            Call MQC_Matrix_Copy_Real2Complex(Mat)
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatC(IndI:IndJ,IndK:IndL) = matrixIn%MatC
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatR(IndI+N-1,IndK+M-1) = MatrixIn%MatC(L,1)
                  Mat%MatR(IndK+M-1,IndI+N-1) = MatrixIn%MatC(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatC(Cnt1,Cnt2) = (MatrixIn%MatC(M,1))
                  Else
                    Mat%MatC(Cnt1,Cnt2) = 0.0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage )
            EndIf
          Else
            Call MQC_Error_A('MatrixIn type not defined in MQC_Vector_Vector_Put', 6, &
                 'MatrixIn%Data_Type', MatrixIn%Data_Type )
          EndIf
        ElseIf (Mat%Data_Type.eq.'Complex') then
          If (MatrixIn%Data_Type.eq.'Integer') then
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatC(IndI:IndJ,IndK:IndL) = MatrixIn%MatI
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatC(IndI+N-1,IndK+M-1) = MatrixIn%MatI(L,1)
                  Mat%MatC(IndK+M-1,IndI+N-1) = MatrixIn%MatI(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatC(Cnt1,Cnt2) = (MatrixIn%MatI(M,1))
                  Else
                    Mat%MatC(Cnt1,Cnt2) = 0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage )
            EndIf
          ElseIf (MatrixIn%Data_Type.eq.'Real') then
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatC(IndI:IndJ,IndK:IndL) = MatrixIn%MatR
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatC(IndI+N-1,IndK+M-1) = MatrixIn%MatR(L,1)
                  Mat%MatC(IndK+M-1,IndI+N-1) = MatrixIn%MatR(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatC(Cnt1,Cnt2) = (MatrixIn%MatR(M,1))
                  Else
                    Mat%MatC(Cnt1,Cnt2) = 0.0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage )
            EndIf
          ElseIf (MatrixIn%Data_Type.eq.'Complex') then
            If(MatrixIn%Storage.eq.'StorFull') then
              Mat%MatC(IndI:IndJ,IndK:IndL) = matrixIn%MatC
            ElseIf(MatrixIn%Storage.eq.'StorSymm') then
              L = 1
              Do M = 1, LenCol
                Do N = 1, M
                  Mat%MatC(IndI+N-1,IndK+M-1) = MatrixIn%MatC(L,1)
                  Mat%MatC(IndK+M-1,IndI+N-1) = MatrixIn%MatC(L,1)
                  L = L + 1
                endDo
              endDo
            ElseIf(MatrixIn%Storage.eq.'StorDiag') then
              M = 1
              N = 1
              Do Cnt1 = IndI, IndJ
                Do Cnt2 = IndK, IndL
                  If (M.eq.N) then
                    Mat%MatC(Cnt1,Cnt2) = (MatrixIn%MatC(M,1))
                  Else
                    Mat%MatC(Cnt1,Cnt2) = 0.0
                  EndIf
                  If(N.eq.LenCol) then
                    N = 1
                  Else
                    N = N + 1
                  EndIf
                EndDo
                M = M + 1
              EndDo
            Else
              Call MQC_Error_A('MatrixIn storage type not recognised in MQC_Matrix_Matrix_Put', 6, &
                   'MatrixIn%Storage', MatrixIn%Storage)
            EndIf
          Else
            Call MQC_Error_A('MatrixIn type not defined in MQC_Vector_Vector_Put', 6, &
                 'MatrixIn%Data_Type', MatrixIn%Data_Type )
          EndIf
        Else
          Call MQC_Error_A('Mat type not defined in MQC_Vector_Vector_Put', 6, &
               'Mat%Data_Type', Mat%Data_Type )
        EndIf
      ElseIf (Mat%Storage.eq.'StorSymm') then
        Call MQC_Matrix_Symm2Full(Mat)
        Call MQC_Matrix_Matrix_Put(Mat,MatrixIn,[IndI,IndJ],[IndK,IndL])
        Write(*,1020)
        If(MQC_Matrix_Test_Symmetric(Mat)) Call MQC_Matrix_Full2Symm(Mat)
      ElseIf (Mat%Storage.eq.'StorDiag') then
        Call MQC_Matrix_Diag2Full(Mat)
        Call MQC_Matrix_Matrix_Put(Mat,MatrixIn,[IndI,IndJ],[IndK,IndL])
        If(MQC_Matrix_Test_Diagonal(Mat)) Call MQC_Matrix_Full2Diag(Mat)
      Else
        Call MQC_Error_A('Mat storage type not recognosed in MQC_Matrix_Vector_Put', 6, &
             'Mat%Storage', Mat%Storage )
      EndIf

      End Subroutine MQC_Matrix_Matrix_Put
!
!
      Function symIndexHash(i,j,k,l) result(hash)
!     This function calculates the correct index for a symmetry packed matrix or rank 4 tensor.
!     If it is a matrix, it is assumed to be low triangular stored row-wise
!
!     X. Sheng 2017.03
!     L. M. Thompson, 2017

      implicit none
      integer(kind=int64),intent(in)::i,j
      integer(kind=int64),intent(in),optional::k,l
      integer(kind=int64)::ij,kl,ijkl,hash

      if(i.gt.j) then
        ij = i*(i-1)/2+j
      else
        ij = j*(j-1)/2+i
      endif
      if(present(k).and.present(l)) then
        if(k.gt.l) then
          kl = k*(k-1)/2+l 
        else
          kl = l*(l-1)/2+k
        endIf
        if(ij.gt.kl) then
          ijkl = ij*(ij-1)/2+kl 
        else
          ijkl = kl*(kl-1)/2+ij 
        endIf
        hash = ijkl
      else
        hash = ij
      endif

      end function symIndexHash
!
!
      Function MQC_ElementMatrixProduct(A,B) Result(C)
!
!     This functions calculates the element by element product of two MQC_matrix
!
!     X. Sheng, 2017.02
!
      implicit none
      type(mqc_matrix),intent(in)::A,B
      type(mqc_matrix)::C
      integer(kind=int64)::i,j

      if (A%NRow .ne. B%NRow .or. A%NCol .ne. B%NCol) then
        call mqc_error_I('The two matrices must have the same dimensions in MQC_elementMatrixProduct', 6, &
             'A%NRow', A%NRow, &
             'B%NRow', B%NRow, &
             'A%NCol', A%NCol, &
             'B%NCol', B%NCol )
      endif

      if (A%data_type .eq. 'Integer') then
        if (B%data_type .eq. 'Integer') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorFull')
              C%matI = A%matI * B%matI
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matI(i,j) = A%matI(i,j) * B%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              do i = 1, A%NRow
                C%matI(i,1) = A%matI(i,i) * B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matI(i,j) = B%matI(i,j) * A%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorSymm')
              C%matI(:,1) = A%matI(:,1) * B%matI(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              do i = 1, A%NRow
                  C%matI(i,1) = A%matI(symIndexHash(i,i),1) * B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              do i = 1, A%NRow
                C%matI(i,1) = B%matI(i,i) * A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              do i = 1, A%NRow
                  C%matI(i,1) = B%matI(symIndexHash(i,i),1) * A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              C%matI(:,1) = A%matI(:,1) * B%matI(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Real') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              C%matR = A%matI * B%matR
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = A%matI(i,j) * B%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = A%matI(i,i) * B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = B%matR(i,j) * A%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorSymm')
              C%matR(:,1) = A%matI(:,1) * B%matR(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = A%matI(symIndexHash(i,i),1) * B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = B%matR(i,i) * A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = B%matR(symIndexHash(i,i),1) * A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              C%matR(:,1) = A%matI(:,1) * B%matR(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Complex') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matI * B%matC
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matI(i,j) * B%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matI(i,i) * B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matC(i,j) * A%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matI(:,1) * B%matC(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matI(symIndexHash(i,i),1) * B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matC(i,i) * A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matC(symIndexHash(i,i),1) * A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matI(:,1) * B%matC(:,1)
            endif
          endif
        endif
      else if (A%data_type .eq. 'Real') then
        if (B%data_type .eq. 'Integer') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              C%matR = A%matR * B%matI
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = A%matR(i,j) * B%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = A%matR(i,i) * B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = B%matI(i,j) * A%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorSymm')
              C%matR(:,1) = A%matR(:,1) * B%matI(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = A%matR(symIndexHash(i,i),1) * B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = B%matI(i,i) * A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = B%matI(symIndexHash(i,i),1) * A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              C%matR(:,1) = A%matR(:,1) * B%matI(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Real') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              C%matR = A%matR * B%matR
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = A%matR(i,j) * B%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = A%matR(i,i) * B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = B%matR(i,j) * A%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorSymm')
              C%matR(:,1) = A%matR(:,1) * B%matR(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = A%matR(symIndexHash(i,i),1) * B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = B%matR(i,i) * A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = B%matR(symIndexHash(i,i),1) * A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              C%matR(:,1) = A%matR(:,1) * B%matR(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Complex') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matR * B%matC
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matR(i,j) * B%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matR(i,i) * B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matC(i,j) * A%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matR(:,1) * B%matC(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matR(symIndexHash(i,i),1) * B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matC(i,i) * A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matC(symIndexHash(i,i),1) * A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matR(:,1) * B%matC(:,1)
            endif
          endif
        endif
      else if (A%data_type .eq. 'Complex') then
        if (B%data_type .eq. 'Integer') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matC * B%matI
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matC(i,j) * B%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matC(i,i) * B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matI(i,j) * A%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matC(:,1) * B%matI(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matC(symIndexHash(i,i),1) * B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matI(i,i) * A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matI(symIndexHash(i,i),1) * A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matC(:,1) * B%matI(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Real') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matC * B%matR
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matC(i,j) * B%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matC(i,i) * B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matR(i,j) * A%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matC(:,1) * B%matR(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = A%matR(symIndexHash(i,i),1) * B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matR(i,i) * A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matR(symIndexHash(i,i),1) * A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matC(:,1) * B%matR(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Complex') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matC * B%matC
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matC(i,j) * B%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matC(i,i) * B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matC(i,j) * A%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matC(:,1) * B%matC(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matC(symIndexHash(i,i),1) * B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matC(i,i) * A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matC(symIndexHash(i,i),1) * A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matC(:,1) * B%matC(:,1)
            endif
          endif
        endif
      endif

      end function MQC_elementMatrixProduct

      Function MQC_ElementMatrixDivide(A,B) Result(C)
!
!     This functions calculates the element by element product of two MQC_matrix
!
!     X. Sheng, 2017.02
!
      implicit none
      type(mqc_matrix),intent(in)::A,B
      type(mqc_matrix)::C
      integer(kind=int64)::i,j

      if (A%NRow .ne. B%NRow .or. A%NCol .ne. B%NCol) then
        call mqc_error_I('The two matrices must have the same dimensions in MQC_elementMatrixProduct', 6, &
             'A%NRow', A%NRow, &
             'B%NRow', B%NRow, &
             'A%NCol', A%NCol, &
             'B%NCol', B%NCol )
      endif

      if (A%data_type .eq. 'Integer') then
        if (B%data_type .eq. 'Integer') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorFull')
              C%matI = A%matI / B%matI
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matI(i,j) = A%matI(i,j) / B%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              do i = 1, A%NRow
                C%matI(i,1) = A%matI(i,i) / B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matI(i,j) = B%matI(i,j) / A%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorSymm')
              C%matI(:,1) = A%matI(:,1) / B%matI(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              do i = 1, A%NRow
                  C%matI(i,1) = A%matI(symIndexHash(i,i),1) / B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              do i = 1, A%NRow
                C%matI(i,1) = B%matI(i,i) / A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              do i = 1, A%NRow
                  C%matI(i,1) = B%matI(symIndexHash(i,i),1) / A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Integer','StorDiag')
              C%matI(:,1) = A%matI(:,1) / B%matI(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Real') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              C%matR = A%matI / B%matR
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = A%matI(i,j) / B%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = A%matI(i,i) / B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = B%matR(i,j) / A%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorSymm')
              C%matR(:,1) = A%matI(:,1) / B%matR(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = A%matI(symIndexHash(i,i),1) / B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = B%matR(i,i) / A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = B%matR(symIndexHash(i,i),1) / A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              C%matR(:,1) = A%matI(:,1) / B%matR(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Complex') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matI / B%matC
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matI(i,j) / B%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matI(i,i) / B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matC(i,j) / A%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matI(:,1) / B%matC(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matI(symIndexHash(i,i),1) / B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matC(i,i) / A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matC(symIndexHash(i,i),1) / A%matI(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matI(:,1) / B%matC(:,1)
            endif
          endif
        endif
      else if (A%data_type .eq. 'Real') then
        if (B%data_type .eq. 'Integer') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              C%matR = A%matR / B%matI
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = A%matR(i,j) / B%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = A%matR(i,i) / B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = B%matI(i,j) / A%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorSymm')
              C%matR(:,1) = A%matR(:,1) / B%matI(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = A%matR(symIndexHash(i,i),1) / B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = B%matI(i,i) / A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = B%matI(symIndexHash(i,i),1) / A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              C%matR(:,1) = A%matR(:,1) / B%matI(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Real') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              C%matR = A%matR / B%matR
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = A%matR(i,j) / B%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = A%matR(i,i) / B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matR(i,j) = B%matR(i,j) / A%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorSymm')
              C%matR(:,1) = A%matR(:,1) / B%matR(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = A%matR(symIndexHash(i,i),1) / B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                C%matR(i,1) = B%matR(i,i) / A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              do i = 1, A%NRow
                  C%matR(i,1) = B%matR(symIndexHash(i,i),1) / A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Real','StorDiag')
              C%matR(:,1) = A%matR(:,1) / B%matR(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Complex') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matR / B%matC
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matR(i,j) / B%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matR(i,i) / B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matC(i,j) / A%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matR(:,1) / B%matC(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matR(symIndexHash(i,i),1) / B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matC(i,i) / A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matC(symIndexHash(i,i),1) / A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matR(:,1) / B%matC(:,1)
            endif
          endif
        endif
      else if (A%data_type .eq. 'Complex') then
        if (B%data_type .eq. 'Integer') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matC / B%matI
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matC(i,j) / B%matI(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matC(i,i) / B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matI(i,j) / A%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matC(:,1) / B%matI(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matC(symIndexHash(i,i),1) / B%matI(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matI(i,i) / A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matI(symIndexHash(i,i),1) / A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matC(:,1) / B%matI(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Real') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matC / B%matR
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matC(i,j) / B%matR(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matC(i,i) / B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matR(i,j) / A%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matC(:,1) / B%matR(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matC(symIndexHash(i,i),1) / B%matR(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matR(i,i) / A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matR(symIndexHash(i,i),1) / A%matC(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matC(:,1) / B%matR(:,1)
            endif
          endif
        else if (B%data_type .eq. 'Complex') then
          if (A%storage .eq. 'StorFull') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              C%matC = A%matC / B%matC
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = A%matC(i,j) / B%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = A%matC(i,i) / B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorSymm') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorFull')
              do i = 1, A%NRow
                do j = 1, A%NCol
                  C%matC(i,j) = B%matC(i,j) / A%matC(symIndexHash(i,j),1)
                enddo
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorSymm')
              C%matC(:,1) = A%matC(:,1) / B%matC(:,1)
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = A%matC(symIndexHash(i,i),1) / B%matC(i,1)
              enddo
            endif
          else if (A%storage .eq. 'StorDiag') then
            if (B%storage .eq. 'StorFull') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                C%matC(i,1) = B%matC(i,i) / A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorSymm') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              do i = 1, A%NRow
                  C%matC(i,1) = B%matC(symIndexHash(i,i),1) / A%matR(i,1)
              enddo
            else if (B%storage .eq. 'StorDiag') then
              call mqc_allocate_matrix(A%NRow,A%NCol,C,'Complex','StorDiag')
              C%matC(:,1) = A%matR(:,1) / B%matC(:,1)
            endif
          endif
        endif
      endif

    end function MQC_elementMatrixDivide
!
!
!     PROCEDURE MQC_Matrix_Test_Symmetric
      Function MQC_Matrix_Test_Symmetric(Matrix,Option) Result(Symmetric)
!
!     This function tests a matrix for symmetry related options and returns true 
!     if the matrix fufills the criteria and false if it does not. This contrasts
!     with the haveSymmetric subroutine which tests how a matrix is packed.
!
!     Options for this subroutine are:
!       'symmetric'      Test if matrix is symmetric
!       'antisymmetric'  Test if matrix is antisymmetric
!       'hermitian'      Test if matrix is hermitian
!       'antihermitian'  Test if matrix is antihermitian
!
!     Lee M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Matrix),Intent(In)::Matrix
      Character(len=*),Optional,Intent(In)::Option
      Logical::Symmetric
      Integer(kind=int64)::I,J,K
      Real(kind=real64)::Thresh=1e-10
      Character(len=64)::myOption

      If(Present(Option)) then
        call string_change_case(Option,'l',myOption)
      else
        myOption = 'symmetric'
      endIf

      Symmetric = .True.
      Select Case (myOption)
      Case ('symmetric')
        If(MQC_Matrix_Rows(Matrix).ne.MQC_Matrix_Columns(Matrix)) then
          Symmetric = .False.
          Return
        Else
          If(Matrix%Storage.eq.'StorSymm'.or.Matrix%Storage.eq.'StorDiag') then
            Return
          ElseIf(Matrix%Storage.eq.'StorFull') then
            Do I = 1, MQC_Matrix_Rows(Matrix)
              Do J = 1, I-1
                If(Matrix%Data_Type.eq.'Integer') then
                  If((Matrix%MatI(I,J) - Matrix%MatI(J,I)).ne.0) then
                    Symmetric = .False.
                    Return
                  EndIf
                elseIf(Matrix%Data_Type.eq.'Real') then
                  If(abs((Matrix%MatR(I,J) - Matrix%MatR(J,I))).gt.Thresh) then
                    Symmetric = .False.
                    Return
                  EndIf
                elseIf(Matrix%Data_Type.eq.'Complex') then
                  If((abs(Real(Matrix%MatC(I,J)) - Real(Matrix%MatC(J,I))).gt.Thresh) .or. &
                    (abs(Aimag(Matrix%MatC(I,J)) - Aimag(Matrix%MatC(J,I))).gt.Thresh)) then
                    Symmetric = .False.
                    Return
                  EndIf
                EndIf
              EndDo
            EndDo
          Else
            Call MQC_Error_A('Matrix type not identified in MQC_Matrix_Test_Symmetric', 6, &
                 'Matrix%Storage', Matrix%Storage )
          EndIf
        EndIf
      Case ('antisymmetric')
        If(MQC_Matrix_Rows(Matrix).ne.MQC_Matrix_Columns(Matrix)) then
          Symmetric = .False.
        Else
          If(Matrix%Storage.eq.'StorSymm') then
            Symmetric = .False.
          ElseIf(Matrix%Storage.eq.'StorDiag') then
            Symmetric = .True.
          ElseIf(Matrix%Storage.eq.'StorFull') then
            Do I = 1, MQC_Matrix_Rows(Matrix)
              Do J = 1, I-1
                If(Matrix%Data_Type.eq.'Integer') then
                  If((Matrix%MatI(I,J) + Matrix%MatI(J,I)).ne.0) then
                    Symmetric = .False.
                    Return
                  EndIf
                elseIf(Matrix%Data_Type.eq.'Real') then
                  If(abs(Matrix%MatR(I,J) + Matrix%MatR(J,I)).gt.Thresh) then
                    Symmetric = .False.
                    Return
                  EndIf
                elseIf(Matrix%Data_Type.eq.'Complex') then
                  If((abs(Real(Matrix%MatC(I,J)) + Real(Matrix%MatC(J,I))).gt.Thresh) .or. &
                    (abs(Aimag(Matrix%MatC(I,J)) + Aimag(Matrix%MatC(J,I))).gt.Thresh)) then
                    Symmetric = .False.
                    Return
                  EndIf
                EndIf
              EndDo
            EndDo
          Else
            Call MQC_Error_A('Matrix type not identified in MQC_Matrix_Test_Symmetric', 6, &
                 'Matrix%Storage', Matrix%Storage )
          EndIf
        EndIf
      Case ('hermitian')
        If(MQC_Matrix_Rows(Matrix).ne.MQC_Matrix_Columns(Matrix)) then
          Symmetric = .False.
        Else
          If(Matrix%Storage.eq.'StorSymm'.or.Matrix%Storage.eq.'StorDiag') then
            If(Matrix%Data_Type.eq.'Integer'.or.Matrix%Data_Type.eq.'Real') then
              Symmetric = .True.
            ElseIf(Matrix%Data_Type.eq.'Complex') then
              Do I = 1, MQC_Matrix_Rows(Matrix)
                Do J = 1, I-1
                  K = symIndexHash(I,J)
                  If(abs(Aimag(Matrix%MatC(K,1))).gt.Thresh) then
                    Symmetric = .False.
                    return
                  EndIf
                EndDo
              EndDo
            Else
              Call MQC_Error_A('Unrecognised data type in mqc_matrix_test_symmetric', 6, &
                   'Matrix%Data_Type', Matrix%Data_Type )
            EndIf
          ElseIf(Matrix%Storage.eq.'StorFull') then
            Do I = 1, MQC_Matrix_Rows(Matrix)
              Do J = 1, I-1
                If(Matrix%Data_Type.eq.'Integer') then
                  If((Matrix%MatI(I,J) - Matrix%MatI(J,I)).ne.0) then
                    Symmetric = .False.
                    Return
                  EndIf
                elseIf(Matrix%Data_Type.eq.'Real') then
                  If(abs(Matrix%MatR(I,J) - Matrix%MatR(J,I)).gt.Thresh) then
                    Symmetric = .False.
                    Return
                  EndIf
                elseIf(Matrix%Data_Type.eq.'Complex') then
                  If((abs(Real(Matrix%MatC(I,J)) - Real(Matrix%MatC(J,I))).gt.Thresh) .or. &
                    (abs(Aimag(Matrix%MatC(I,J)) + Aimag(Matrix%MatC(J,I))).gt.Thresh)) then
                    Symmetric = .False.
                    Return
                  EndIf
                EndIf
              EndDo
            EndDo
          Else
            Call MQC_Error_A('Matrix type not identified in MQC_Matrix_Test_Symmetric', 6, &
                 'Matrix%Storage', Matrix%Storage )
          EndIf
        EndIf
      Case ('antihermitian')
        If(MQC_Matrix_Rows(Matrix).ne.MQC_Matrix_Columns(Matrix)) then
          Symmetric = .False.
        Else
          If(Matrix%Storage.eq.'StorSymm'.or.Matrix%Storage.eq.'StorDiag') then
            Symmetric = .False.
          ElseIf(Matrix%Storage.eq.'StorFull') then
            Do I = 1, MQC_Matrix_Rows(Matrix)
              Do J = 1, I-1
                If(Matrix%Data_Type.eq.'Integer') then
                  If((Matrix%MatI(I,J) + Matrix%MatI(J,I)).ne.0) then
                    Symmetric = .False.
                    Return
                  EndIf
                elseIf(Matrix%Data_Type.eq.'Real') then
                  If((Matrix%MatR(I,J) + Matrix%MatR(J,I)).ne.0.0) then
                    Symmetric = .False.
                    Return
                  EndIf
                elseIf(Matrix%Data_Type.eq.'Complex') then
                  If((Real(Matrix%MatC(I,J)) + Real(Matrix%MatC(J,I))).ne.0.0 .or. &
                    (Aimag(Matrix%MatC(I,J)) - Aimag(Matrix%MatC(J,I))).ne.0.0) then
                    Symmetric = .False.
                    Return
                  EndIf
                EndIf
              EndDo
            EndDo
          Else
            Call MQC_Error_A('Matrix type not identified in MQC_Matrix_Test_Symmetric', 6, &
                 'Matrix%Storage', Matrix%Storage)

!            Call MQC_Error('Matrix type not identified in MQC_Matrix_Test_Symmetric')
          EndIf
        EndIf
      Case Default
        Call MQC_Error_A('Invalid option sent to mqc_matrix_test_symmetric', 6, &
             'myOption', myOption )
      End Select

      Return
      End Function MQC_Matrix_Test_Symmetric
!
!
!     PROCEDURE MQC_Matrix_Test_Diagonal
      Function MQC_Matrix_Test_Diagonal(Matrix) Result(Diagonal)
!
!     This function returns true if a matrix is diagonal and false if it is not.
!
!     Lee M. Thompson, 2017.
!
      Implicit None
      Class(MQC_Matrix),Intent(In)::Matrix
      Logical::Diagonal
      Integer(kind=int64)::I,J,II
      Real(kind=real64)::Thresh=1e-14

      Diagonal = .True.
      If(Matrix%Storage.eq.'StorDiag') then
        Return
      ElseIf(Matrix%Storage.eq.'StorSymm') then
        Do I = 1, MQC_Matrix_Rows(Matrix)
          II = (I*(I-1))/2
          Do J = 1, I-1
            If(Matrix%Data_Type.eq.'Integer') then
              If(Matrix%MatI(II+J,1).ne.0) then
                Diagonal = .False.
                Return
              EndIf
            ElseIf(Matrix%Data_Type.eq.'Real') then
              If(Matrix%MatR(II+J,1).ne.0.0) then
                Diagonal = .False.
                Return
              EndIf
            ElseIf(Matrix%Data_Type.eq.'Complex') then 
              If(Real(Matrix%MatC(II+J,1)).ne.0.0.or.Aimag(Matrix%MatC(II+J,1)).ne.0.0) then
                Diagonal = .False.
                Return
              EndIf
            EndIf
          EndDo
        EndDo
      ElseIf(Matrix%Storage.eq.'StorFull') then
        Do I = 1, MQC_Matrix_Rows(Matrix)
          Do J = 1 , MQC_Matrix_Columns(Matrix)
            If(I.eq.J) cycle
            If(Matrix%Data_Type.eq.'Integer') then
              If(Matrix%MatI(I,J).ne.0) then
                Diagonal = .False.
                Return
              EndIf
            ElseIf(Matrix%Data_Type.eq.'Real') then
              If(Matrix%MatR(I,J).ne.0.0) then
                Diagonal = .False.
                Return
              EndIf
            ElseIf(Matrix%Data_Type.eq.'Complex') then 
              If(Real(Matrix%MatC(I,J)).ne.0.0.or.Aimag(Matrix%MatC(I,J)).ne.0.0) then
                Diagonal = .False.
                Return
              EndIf
            EndIf
          EndDo
        EndDo
      Else
        Call MQC_Error_A('Matrix type not identified in MQC_Matrix_Test_Diagonal', 6, &
             'Matrix%Storage', Matrix%Storage )
      EndIf

      Return
      End Function MQC_Matrix_Test_Diagonal
!
!
!     PROCEDURE MQC_Allocate_Matrix
      Subroutine MQC_Allocate_Matrix(M,N,Matrix,Data_Type,Storage)
!
!     This subroutine is used to allocate memory for a matrix variable
!     belonging to the MQC_Matrix class.
!
!     H. P. Hratchian, 2016.
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::M,N
      Class(MQC_Matrix),Intent(InOut)::Matrix
      Character(Len=*),Intent(In)::Data_Type,Storage
      Integer(kind=int64)::MCur,NCur
!
      If (M.lt.0 .or. N.lt.0) then
        Call MQC_Error_I('Dimensions less than zero in MQC_Allocate_Matrix', 6, &
             'M', M, &
             'N', N )
      EndIf

      Call MQC_Deallocate_Matrix(Matrix)

      If(Storage.eq.'StorFull' .or. Storage.eq.'Full' .or. Storage.eq.'full' .or. Storage.eq.'storfull') then
        Matrix%Storage = 'StorFull'
        MCur = M
        NCur = N
      ElseIf(Storage.eq.'StorSymm' .or. Storage.eq.'Symm' .or. Storage.eq.'symm' .or. Storage.eq.'storsymm') then
         If(M.ne.N) Call MQC_Error_I('Attempting to assign non-square matrix to lower triangular in MQC_Allocate_Matrix', 6, &
              'M', M, &
              'N', N )
        Matrix%Storage = 'StorSymm'
        MCur = (M*(M+1))/2
        NCur = 1
      ElseIf(Storage.eq.'StorDiag' .or. Storage.eq.'Diag' .or. Storage.eq.'diag' .or. Storage.eq.'stordiag') then
        Matrix%Storage = 'StorDiag'
        MCur = Min(M,N)
        NCur = 1
      ElseIf(Storage.eq.'StorTriD' .or. Storage.eq.'TriD' .or. Storage.eq.'trid' .or. Storage.eq.'stortrid') then
        Matrix%Storage = 'StorTriD'
        If(M.ne.N) then
          MCur = 3*Min(M,N)-1
        Else
          MCur = 3*N-2
        EndIf
        NCur = 1
      Else
        Call MQC_Error_A('Storage type not recognised in MQC_Allocate_Matrix', 6, &
             'Storage', Storage )
      EndIf

      If(Data_Type.eq.'Integer' .or. Data_Type.eq.'integer' .or. Data_Type.eq.'Int' .or. Data_Type.eq.'int') then
        Matrix%Data_Type = 'Integer'
        If(Allocated(Matrix%MatI)) DeAllocate(Matrix%MatI)
        Allocate(Matrix%MatI(MCur,NCur))
      ElseIf(Data_Type.eq.'Real' .or. Data_Type.eq.'real') then
        Matrix%Data_Type = 'Real'
        If(Allocated(Matrix%MatR)) DeAllocate(Matrix%MatR)
        Allocate(Matrix%MatR(MCur,NCur))
      ElseIf(Data_Type.eq.'Complex' .or. Data_Type.eq.'complex') then
        Matrix%Data_Type = 'Complex'
        If(Allocated(Matrix%MatC)) DeAllocate(Matrix%MatC)
        Allocate(Matrix%MatC(MCur,NCur))
      EndIf

      Matrix%NRow = M
      Matrix%NCol = N
!
      Return
      End Subroutine MQC_Allocate_Matrix
!
!
!     PROCEDURE MQC_Deallocate_Matrix
      Subroutine MQC_Deallocate_Matrix(Matrix)
!
!     This subroutine is used to deallocate memory for a matrix variable
!     belonging to the MQC_Matrix class.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Matrix),Intent(InOut)::Matrix
!
      If(Allocated(Matrix%MatI)) DeAllocate(Matrix%MatI)
      If(Allocated(Matrix%MatR)) DeAllocate(Matrix%MatR)
      If(Allocated(Matrix%MatC)) DeAllocate(Matrix%MatC)
      If(.not.Allocated(Matrix%MatR).and..not.Allocated(Matrix%MatI).and..not.Allocated(Matrix%MatC)) then
        Matrix%NCol = 0
        Matrix%NRow = 0
        Matrix%Data_Type = ''
        Matrix%Storage = ''
      EndIf
!
      Return
      End Subroutine MQC_Deallocate_Matrix
!
!
!     PROCEDURE MQC_Matrix_isAllocated
      Function MQC_Matrix_isAllocated(Matrix) Result(isAllocated)
!
!     This function returns if an mqc matrix is allocated or not.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Matrix),Intent(InOut)::Matrix
      Logical::isAllocated
!
      isAllocated = .False.
      If(Allocated(Matrix%MatR)) isAllocated = .True.
      If(Allocated(Matrix%MatI)) isAllocated = .True.
      If(Allocated(Matrix%MatC)) isAllocated = .True.
!
      Return
      End Function MQC_Matrix_isAllocated
!
!
!     PROCEDURE MQC_Set_IntegerArray2Matrix
      Subroutine MQC_Set_IntegerArray2Matrix(MatrixOut,ArrayIn)
!
!     This routine is used to set an MQC matrix equal to an integer rank-2 array.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::MatrixOut
      Integer(kind=int64),Dimension(:,:),Intent(In)::ArrayIn
!
      Call MQC_Deallocate_Matrix(MatrixOut)
      Call MQC_Allocate_Matrix(Size(ArrayIn,1),Size(ArrayIn,2),MatrixOut,'Integer', &
            'StorFull')
      MatrixOut%MatI = ArrayIn
!
      Return
      End Subroutine MQC_Set_IntegerArray2Matrix
!
!
!     PROCEDURE MQC_Set_RealArray2Matrix
      Subroutine MQC_Set_RealArray2Matrix(MatrixOut,ArrayIn)
!
!     This routine is used to set an MQC matrix equal to a real rank-2 array.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::MatrixOut
      Real(kind=real64),Dimension(:,:),Intent(In)::ArrayIn
!
      Call MQC_Deallocate_Matrix(MatrixOut)
      Call MQC_Allocate_Matrix(Size(ArrayIn,1),Size(ArrayIn,2),MatrixOut,'Real', &
        'StorFull')
      MatrixOut%MatR = ArrayIn
!
      Return
      End Subroutine MQC_Set_RealArray2Matrix
!
!
!     PROCEDURE MQC_Set_ComplexArray2Matrix
      Subroutine MQC_Set_ComplexArray2Matrix(MatrixOut,ArrayIn)
!
!     This routine is used to set an MQC matrix equal to a complex rank-2 array.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::MatrixOut
      Complex(Kind=real64),Dimension(:,:),Intent(In)::ArrayIn
!
      Call MQC_Deallocate_Matrix(MatrixOut)
      Call MQC_Allocate_Matrix(Size(ArrayIn,1),Size(ArrayIn,2),MatrixOut,'Complex', &
        'StorFull')
      MatrixOut%MatC = ArrayIn
!
      Return
      End Subroutine MQC_Set_ComplexArray2Matrix
!
!     PROCEDURE MQC_Set_Matrix2IntegerArray
      Subroutine MQC_Set_Matrix2IntegerArray(ArrayOut,MatrixIn)
!
!     This routine is used to set a rank-2 allocatable array equal to an MQC matrix.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Dimension(:,:),Allocatable,Intent(InOut)::ArrayOut
      Type(MQC_Matrix),Intent(In)::MatrixIn
      Integer(kind=int64)::ColsMatIn,RowsMatIn,I,J
!
!     Need to have a way of determining output based on matrix packing which
!     depends on type of array in as well ( here assuming full).
!     Need to have a way of dealing with row vs column-wise arrays (here
!     assuming cols).
!
      RowsMatIn = MQC_Matrix_Rows(MatrixIn)
      If(MatrixIn%Storage.eq.'StorFull') then
        ColsMatIn = MQC_Matrix_Columns(MatrixIn)
      ElseIf(MatrixIn%Storage.eq.'StorTriD') then
        If(MQC_Matrix_Columns(MatrixIn).ne.RowsMatIn) then
          ColsMatIn = 3*Min(MQC_Matrix_Columns(MatrixIn),RowsMatIn)-1
        Else
          ColsMatIn = 3*RowsMatIn-2
        EndIf ! Xianghai: Is this the right dimensions for triD?
      Elseif (MatrixIn%Storage.eq.'StorSymm' .or. MatrixIn%Storage.eq.'StorDiag') then
        ColsMatIn = 1
      Else
        Call MQC_Error_A('MatrixIn storage type unkown in MQC_Matrix2Array', 6, &
             'MatrixIn%Storage', MatrixIn%Storage )
      EndIf
!
      If(.not.Allocated(ArrayOut)) then
        Allocate(ArrayOut(RowsMatIn,ColsMatIn))
      ElseIf(RowsMatIn.ne.Size(ArrayOut,1).or.ColsMatIn.ne.Size(ArrayOut,2)) then
        Deallocate(ArrayOut)
        Allocate(ArrayOut(RowsMatIn,ColsMatIn))
      EndIf

      If(MatrixIn%Data_Type.eq.'Integer') then
        ArrayOut = MatrixIn%MatI
      ElseIf(MatrixIn%Data_Type.eq.'Real') then
        ArrayOut = MatrixIn%MatR
      ElseIf(MatrixIn%Data_Type.eq.'Complex') then
        ArrayOut = MatrixIn%MatC
      Else
        Call MQC_Error_A('MatrixIn type unkown in MQC_Matrix2Array', 6, &
             'MatrixIn%Data_Type', MatrixIn%Data_Type )
      EndIf
!
      Return
      End Subroutine MQC_Set_Matrix2IntegerArray
!
!
!     PROCEDURE MQC_Set_Matrix2RealArray
      Subroutine MQC_Set_Matrix2RealArray(ArrayOut,MatrixIn)
!
!     This routine is used to set a rank-2 allocatable array equal to an MQC matrix.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Real(kind=real64),Dimension(:,:),Allocatable,Intent(InOut)::ArrayOut
      Type(MQC_Matrix),Intent(In)::MatrixIn
      Integer(kind=int64)::ColsMatIn,RowsMatIn,I,J
!
!     Need to have a way of determining output based on matrix packing which
!     depends on type of array in as well ( here assuming full).
!     Need to have a way of dealing with row vs column-wise arrays (here
!     assuming cols).
!
      RowsMatIn = MQC_Matrix_Rows(MatrixIn)
      If(MatrixIn%Storage.eq.'StorFull') then
        ColsMatIn = MQC_Matrix_Columns(MatrixIn)
      ElseIf(MatrixIn%Storage.eq.'StorTriD') then
        If(MQC_Matrix_Columns(MatrixIn).ne.RowsMatIn) then
          ColsMatIn = 3*Min(MQC_Matrix_Columns(MatrixIn),RowsMatIn)-1
        Else
          ColsMatIn = 3*RowsMatIn-2
        EndIf ! Xianghai: Is this the right dimensions for triD?
      Elseif (MatrixIn%Storage.eq.'StorSymm' .or. MatrixIn%Storage.eq.'StorDiag') then
        ColsMatIn = 1
      Else
        Call MQC_Error_A('MatrixIn storage type unkown in MQC_Matrix2Array', 6, &
             'MatrixIn%Storage', MatrixIn%Storage )
      EndIf
!
      If(.not.Allocated(ArrayOut)) then
        Allocate(ArrayOut(RowsMatIn,ColsMatIn))
      ElseIf(RowsMatIn.ne.Size(ArrayOut,1).or.ColsMatIn.ne.Size(ArrayOut,2)) then
        Deallocate(ArrayOut)
        Allocate(ArrayOut(RowsMatIn,ColsMatIn))
      EndIf

      If(RowsMatIn.gt.0.and.ColsMatIn.gt.0) then
        If(MatrixIn%Data_Type.eq.'Integer') then
          ArrayOut = MatrixIn%MatI
        ElseIf(MatrixIn%Data_Type.eq.'Real') then
          ArrayOut = MatrixIn%MatR
        ElseIf(MatrixIn%Data_Type.eq.'Complex') then
          ArrayOut = MatrixIn%MatC
        Else
          Call MQC_Error_A('MatrixIn type unkown in MQC_Matrix2Array', 6, &
               'MatrixIn%Data_Type', MatrixIn%Data_Type )
        EndIf
      EndIf
!
      Return
      End Subroutine MQC_Set_Matrix2RealArray
!
!
!     PROCEDURE MQC_Set_Matrix2ComplexArray
      Subroutine MQC_Set_Matrix2ComplexArray(ArrayOut,MatrixIn)
!
!     This routine is used to set a rank-2 allocatable array equal to an MQC matrix.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Complex(Kind=real64),Dimension(:,:),Allocatable,Intent(InOut)::ArrayOut
      Type(MQC_Matrix),Intent(In)::MatrixIn
      Integer(kind=int64)::ColsMatIn,RowsMatIn,I,J
!
!     Need to have a way of determining output based on matrix packing which
!     depends on type of array in as well ( here assuming full).
!     Need to have a way of dealing with row vs column-wise arrays (here
!     assuming cols).
!
      RowsMatIn = MQC_Matrix_Rows(MatrixIn)
      If(MatrixIn%Storage.eq.'StorFull') then
        ColsMatIn = MQC_Matrix_Columns(MatrixIn)
      ElseIf(MatrixIn%Storage.eq.'StorTriD') then
        If(MQC_Matrix_Columns(MatrixIn).ne.RowsMatIn) then
          ColsMatIn = 3*Min(MQC_Matrix_Columns(MatrixIn),RowsMatIn)-1
        Else
          ColsMatIn = 3*RowsMatIn-2
        EndIf ! Xianghai: Is this the right dimensions for triD?
      Elseif (MatrixIn%Storage.eq.'StorSymm' .or. MatrixIn%Storage.eq.'StorDiag') then
        ColsMatIn = 1
      Else
        Call MQC_Error_A('MatrixIn storage type unkown in MQC_Matrix2Array', 6, &
             'MatrixIn%Storage', MatrixIn%Storage )
      EndIf
!
      If(.not.Allocated(ArrayOut)) then
        Allocate(ArrayOut(RowsMatIn,ColsMatIn))
      ElseIf(RowsMatIn.ne.Size(ArrayOut,1).or.ColsMatIn.ne.Size(ArrayOut,2)) then
        Deallocate(ArrayOut)
        Allocate(ArrayOut(RowsMatIn,ColsMatIn))
      EndIf

      If(RowsMatIn.gt.0.and.ColsMatIn.gt.0) then
        If(MatrixIn%Data_Type.eq.'Integer') then
          ArrayOut = MatrixIn%MatI
        ElseIf(MatrixIn%Data_Type.eq.'Real') then
          ArrayOut = MatrixIn%MatR
        ElseIf(MatrixIn%Data_Type.eq.'Complex') then
          ArrayOut = MatrixIn%MatC
        Else
          Call MQC_Error_A('MatrixIn type unkown in MQC_Matrix2Array', 6, &
               'MatrixIn%Data_Type', MatrixIn%Data_Type )
        EndIf
      EndIf
!
      Return
      End Subroutine MQC_Set_Matrix2ComplexArray
!
!
!     PROCEDURE MQC_Set_Matrix2Matrix
      Subroutine MQC_Set_Matrix2Matrix(MatrixOut,MatrixIn)
!
!     This routine is used to set an MQC matrix equal to another MQC
!     matrix.
!
!     L. M. Thompson, 2016.
!
!     Variable Declarations.
      Implicit None
      Class(MQC_Matrix),Intent(InOut)::MatrixOut
      Class(MQC_Matrix),Intent(In)::MatrixIn
!
      Call MQC_Deallocate_Matrix(MatrixOut)
      If(MQC_Matrix_HaveReal(MatrixIn)) then
        Call MQC_Allocate_Matrix(MQC_Matrix_Rows(MatrixIn),MQC_Matrix_Columns(MatrixIn), &
          MatrixOut,'Real',MatrixIn%Storage)
        MatrixOut%MatR = MatrixIn%MatR
      EndIf
      If(MQC_Matrix_HaveInteger(MatrixIn)) then
        Call MQC_Allocate_Matrix(MQC_Matrix_Rows(MatrixIn),MQC_Matrix_Columns(MatrixIn), &
          MatrixOut,'Integer',MatrixIn%Storage)
        MatrixOut%MatI = MatrixIn%MatI
      EndIf
      If(MQC_Matrix_HaveComplex(MatrixIn)) then
        Call MQC_Allocate_Matrix(MQC_Matrix_Rows(MatrixIn),MQC_Matrix_Columns(MatrixIn), &
          MatrixOut,'Complex',MatrixIn%Storage)
        MatrixOut%MatC = MatrixIn%MatC
      EndIf
!
      !MatrixOut%Column = MatrixIn%Column
!
      Return
      End Subroutine MQC_Set_Matrix2Matrix
!
!
!     PROCEDURE MQC_Print_Matrix_Algebra1
      Subroutine MQC_Print_Matrix_Algebra1(Matrix,IOut,Header,Blank_At_Top, &
        Blank_At_Bottom)
!
!     This subroutine is used to print a MQC_Matrix type variable.
!     All matricies are converted to full before printing.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Class(MQC_Matrix),Intent(In)::Matrix
      Character(Len=*),Intent(In)::Header
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64)::I,J,II,ILim,NCols,NRows,IFirst,ILast
      Integer(kind=int64),Parameter::ColWidth=10 ! Xianghai: can we make it bigger?
!
 1000 Format(1x,A)
 1001 Format(5x,10(7x,I7))
 1002 Format(3x,'This is a diagonal matrix, so only diagonal elements are printed')
 1003 Format(10(17x,I7))
 1020 Format( " " )
 2001 Format(1x,I7,10I14)
 2002 Format(1x,I7,10F14.6)
 2003 Format(1x,I7,10(F12.5,F11.5,"i"))
!
!     Output will depend on matrix packing type (assumed full) or column/row
      NCols = MQC_Matrix_Columns(Matrix)
      NRows = MQC_Matrix_Rows(Matrix)

      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      Write(IOut,1000) TRIM(Header)

!     Printing for full storage
      If(Matrix%Storage.eq.'StorFull') then
        Do IFirst = 1,NCols,ColWidth
          ILast = Min(IFirst+ColWidth-1,NCols)
          If(Matrix%Data_Type.eq.'Complex') then
            Write(IOut,1003) (I,I=IFirst,ILast)
          Else
            Write(IOut,1001) (I,I=IFirst,ILast)
          EndIf
          Do I = 1,NRows
            If(Matrix%Data_Type.eq.'Integer') then
              Write(IOut,2001) I, (Matrix%MatI(I,J),J=IFirst,ILast)
              !If(Matrix%Column) Write(IOut,2001) I, (Matrix%MatI(I,J),J=IFirst,ILast)
              !If(.not.Matrix%Column) Write(IOut,2001) I, (Matrix%MatI(J,I),J=IFirst,ILast)
            ElseIf(Matrix%Data_Type.eq.'Real') then
              Write(IOut,2002) I, (Matrix%MatR(I,J),J=IFirst,ILast)
              !If(Matrix%Column) Write(IOut,2002) I, (Matrix%MatR(I,J),J=IFirst,ILast)
              !If(.not.Matrix%Column) Write(IOut,2002) I, (Matrix%MatR(J,I),J=IFirst,ILast)
            ElseIf(Matrix%Data_Type.eq.'Complex') then
              Write(IOut,2003) I, (Matrix%MatC(I,J),J=IFirst,ILast)
            Else
              Call MQC_Error_A('Matrix data type unspecified in MQC_Print_Matrix_Algebra1', 6, &
                 'Matrix%Data_Type', Matrix%Data_Type )
            EndIf
          EndDo
        EndDo

!     Printing for Symm storage
      ElseIf(Matrix%Storage.eq.'StorSymm') then
        Do IFirst = 1,NCols,ColWidth
          ILast = Min(IFirst+ColWidth-1,NCols)
          If(Matrix%Data_Type.eq.'Complex') then
            Write(IOut,1003) (I,I=IFirst,ILast)
          Else
            Write(IOut,1001) (I,I=IFirst,ILast)
          EndIf
          Do I = 1,NRows
            If(IFirst.gt.I) Cycle
            ILim = Min(I-IFirst+1,ColWidth)
            II = (I*(I-1))/2 + IFirst
!            Do J = 1,ILim
            If(Matrix%Data_Type.eq.'Integer') then
              Write(IOut,2001) I, (Matrix%MatI(II+J-1,1),J=1,ILim)
            ElseIf(Matrix%Data_Type.eq.'Real') then
              Write(IOut,2002) I, (Matrix%MatR(II+J-1,1),J=1,ILim)
            ElseIf(Matrix%Data_Type.eq.'Complex') then
              Write(IOut,2003) I, (Matrix%MatC(II+J-1,1),J=1,ILim)
            Else
              Call MQC_Error_A('Matrix data type unspecified in MQC_Print_Matrix_Algebra1', 6, &
                 'Matrix%Data_Type', Matrix%Data_Type )
            EndIf
          EndDo
        EndDo

!     Printing for diagonal storage
      ElseIf(Matrix%Storage.eq.'StorDiag') then
        Write(IOut,1002)
        Do I = 1,Min(NRows,NCols)
          If(Matrix%Data_Type.eq.'Integer') then
            Write(IOut,2001) I, Matrix%MatI(I,1)
          ElseIf(Matrix%Data_Type.eq.'Real') then
            Write(IOut,2002) I, Matrix%MatR(I,1)
          ElseIf(Matrix%Data_Type.eq.'Complex') then
            Write(IOut,2003) I, Matrix%MatC(I,1)
          Else
            Call MQC_Error_A('Matrix data type unspecified in MQC_Print_Matrix_Algebra1', 6, &
                 'Matrix%Data_Type', Matrix%Data_Type )
          EndIf
        EndDo

      Else
        Call MQC_Error_A('Matrix storage type unspecified in MQC_Print_Matrix_Algebra1', 6, &
             'Matrix%Storage', Matrix%Storage )
      EndIf

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
!
      Return
      End Subroutine MQC_Print_Matrix_Algebra1
!
!
!     PROCEDURE MQC_Matrix_Copy_Int2Real
      Subroutine MQC_Matrix_Copy_Int2Real(Matrix)
!
!     This subroutine copies an integer MQC_Matrix into it's real matrix
!     space.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::Matrix
!
      If(.not.MQC_Matrix_HaveInteger(Matrix)) Call MQC_Error_L('Incoming matrix not integer in MQC_Matrix_Copy_Int2Real', 6, &
           'MQC_Matrix_HaveInteger(Matrix)', MQC_Matrix_HaveInteger(Matrix) )
      Matrix%MatR = Matrix%MatI
      If(Allocated(Matrix%MatI)) Deallocate(Matrix%MatI)
      Matrix%Data_Type = 'Real'
!
      Return
      End Subroutine MQC_Matrix_Copy_Int2Real
!
!
!     PROCEDURE MQC_Matrix_Copy_Int2Complex
      Subroutine MQC_Matrix_Copy_Int2Complex(Matrix)
!
!     This subroutine copies an integer MQC_Matrix into it's complex matrix
!     space.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::Matrix
!
      If(.not.MQC_Matrix_HaveInteger(Matrix)) Call MQC_Error_L('Incoming matrix not integer in MQC_Matrix_Copy_Int2Complex', 6, &
           'MQC_Matrix_HaveInteger(Matrix)', MQC_Matrix_HaveInteger(Matrix) )
      Matrix%MatC = Matrix%MatI
      If(Allocated(Matrix%MatI)) Deallocate(Matrix%MatI)
      Matrix%Data_Type = 'Complex'
!
      Return
      End Subroutine MQC_Matrix_Copy_Int2Complex
!
!
!     PROCEDURE MQC_Matrix_Copy_Real2Int
      Subroutine MQC_Matrix_Copy_Real2Int(Matrix)
!
!     This subroutine copies a real MQC_Matrix into it's integer matrix
!     space.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::Matrix
!
      If(.not.MQC_Matrix_HaveReal(Matrix)) Call MQC_Error_L('Incoming matrix not real in MQC_Matrix_Copy_Real2Int', 6, &
           'MQC_Matrix_HaveReal(Matrix)', MQC_Matrix_HaveReal(Matrix) )
      Matrix%MatI = Matrix%MatR
      If(Allocated(Matrix%MatR)) Deallocate(Matrix%MatR)
      Matrix%Data_Type = 'Integer'
!
      Return
      End Subroutine MQC_Matrix_Copy_Real2Int
!
!
!     PROCEDURE MQC_Matrix_Copy_Real2Complex
      Subroutine MQC_Matrix_Copy_Real2Complex(Matrix)
!
!     This subroutine copies a real MQC_Matrix into it's complex matrix
!     space.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::Matrix
!
      If(.not.MQC_Matrix_HaveReal(Matrix)) Call MQC_Error_L('Incoming matrix not real in MQC_Matrix_Copy_Real2Complex', 6, &
           'MQC_Matrix_HaveReal(Matrix)', MQC_Matrix_HaveReal(Matrix) )
      Matrix%MatC = Matrix%MatR
      If(Allocated(Matrix%MatR)) Deallocate(Matrix%MatR)
      Matrix%Data_Type = 'Complex'
!
      Return
      End Subroutine MQC_Matrix_Copy_Real2Complex
!
!
!     PROCEDURE MQC_Matrix_Copy_Complex2Int
      Subroutine MQC_Matrix_Copy_Complex2Int(Matrix)
!
!     This subroutine copies a complex MQC_Matrix into it's integer matrix
!     space.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::Matrix
!
      If(.not.MQC_Matrix_HaveComplex(Matrix)) Call MQC_Error_L('Incoming matrix not complex in MQC_Matrix_Copy_Complex2Int', 6, &
           'MQC_Matrix_HaveComplex(Matrix)', MQC_Matrix_HaveComplex(Matrix) )
      Matrix%MatI = Real(Matrix%MatC)
      If(Allocated(Matrix%MatC)) Deallocate(Matrix%MatC)
      Matrix%Data_Type = 'Integer'
!
      Return
      End Subroutine MQC_Matrix_Copy_Complex2Int
!
!
!     PROCEDURE MQC_Matrix_Copy_Complex2Real
      Subroutine MQC_Matrix_Copy_Complex2Real(Matrix)
!
!     This subroutine copies a complex MQC_Matrix into it's real matrix
!     space.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::Matrix
!
      If(.not.MQC_Matrix_HaveComplex(Matrix)) Call MQC_Error_L('Incoming matrix not real in MQC_Matrix_Copy_Complex2Real', 6, &
           'MQC_Matrix_HaveComplex(Matrix)', MQC_Matrix_HaveComplex(Matrix) )
      Matrix%MatR = Real(Matrix%MatC)
      If(Allocated(Matrix%MatC)) Deallocate(Matrix%MatC)
      Matrix%Data_Type = 'Real'
!
      Return
      End Subroutine MQC_Matrix_Copy_Complex2Real
!
!
!     PROCEDURE MQC_Matrix_Rows
      Function MQC_Matrix_Rows(Matrix)
!
!     This subroutine is used to return the number of rows of an MQC matrix.
!     If the matrix is NOT allocated, the number of rows is returned as 0.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64)::MQC_Matrix_Rows
      Type(MQC_Matrix),Intent(In)::Matrix
!
      MQC_Matrix_Rows = Matrix%NRow
!
      Return
      End Function MQC_Matrix_Rows
!
!
!     PROCEDURE MQC_Matrix_Columns
      Function MQC_Matrix_Columns(Matrix)
!
!     This subroutine is used to return the number of columns of an MQC matrix.
!     If the matrix is NOT allocated, the number of columns is returned as 0.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Integer(kind=int64)::MQC_Matrix_Columns
      Type(MQC_Matrix),Intent(In)::Matrix
!
      MQC_Matrix_Columns = Matrix%NCol
!
      Return
      End Function MQC_Matrix_Columns
!
!
!     PROCEDURE MQC_Matrix_HaveReal
      Function MQC_Matrix_HaveReal(Matrix)
!
!     This function returns TRUE or FALSE indicating whether Matrix has an
!     allocated real matrix.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Matrix_HaveReal
      Type(MQC_Matrix),Intent(In)::Matrix
!
      MQC_Matrix_HaveReal = Allocated(Matrix%MatR)
!
      Return
      End Function MQC_Matrix_HaveReal
!
!
!     PROCEDURE MQC_Matrix_HaveInteger
      Function MQC_Matrix_HaveInteger(Matrix)
!
!     This function returns TRUE or FALSE indicating whether Matrix has an
!     allocated integer matrix.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Matrix_HaveInteger
      Type(MQC_Matrix),Intent(In)::Matrix
!
      MQC_Matrix_HaveInteger = Allocated(Matrix%MatI)
!
      Return
      End Function MQC_Matrix_HaveInteger
!
!
!     PROCEDURE MQC_Matrix_HaveComplex
      Function MQC_Matrix_HaveComplex(Matrix)
!
!     This function returns TRUE or FALSE indicating whether Matrix has an
!     allocated complex matrix.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Matrix_HaveComplex
      Type(MQC_Matrix),Intent(In)::Matrix
!
      MQC_Matrix_HaveComplex = Allocated(Matrix%MatC)
!
      Return
      End Function MQC_Matrix_HaveComplex
!
!
!     PROCEDURE MQC_Matrix_HaveFull   
      Function MQC_Matrix_HaveFull(Matrix)
!
!     This function returns TRUE or FALSE indicating whether Matrix is 
!     stored full.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Matrix_HaveFull   
      Type(MQC_Matrix),Intent(In)::Matrix
!
      If(Matrix%Storage.eq.'StorFull') then
        MQC_Matrix_HaveFull = .True. 
      Else
        MQC_Matrix_HaveFull = .False.
      EndIf
!
      Return
      End Function MQC_Matrix_HaveFull   
!
!
!     PROCEDURE MQC_Matrix_HaveSymmetric
      Function MQC_Matrix_HaveSymmetric(Matrix)
!
!     This function returns TRUE or FALSE indicating whether Matrix is 
!     stored symmetric.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Matrix_HaveSymmetric 
      Type(MQC_Matrix),Intent(In)::Matrix
!
      If(Matrix%Storage.eq.'StorSymm') then
        MQC_Matrix_HaveSymmetric = .True. 
      Else
        MQC_Matrix_HaveSymmetric = .False.
      EndIf
!
      Return
      End Function MQC_Matrix_HaveSymmetric   
!
!
!     PROCEDURE MQC_Matrix_HaveDiagonal 
      Function MQC_Matrix_HaveDiagonal(Matrix)
!
!     This function returns TRUE or FALSE indicating whether Matrix is 
!     stored diagonal.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Logical::MQC_Matrix_HaveDiagonal   
      Type(MQC_Matrix),Intent(In)::Matrix
!
      If(Matrix%Storage.eq.'StorDiag') then
        MQC_Matrix_HaveDiagonal = .True. 
      Else
        MQC_Matrix_HaveDiagonal = .False.
      EndIf
!
      Return
      End Function MQC_Matrix_HaveDiagonal   
!
!
!     PROCEDURE MQC_Matrix_Transpose
      Function MQC_Matrix_Transpose(Matrix)
!
!     This function returns the transpose of an MQC matrx.
!
!     L. M. Thompson, 2016.
!     X. Sheng, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::MQC_Matrix_Transpose
      Class(MQC_Matrix),Intent(In)::Matrix
!
      If(MQC_Matrix_HaveReal(Matrix)) then
        Call MQC_Allocate_Matrix(Matrix%NCol,Matrix%NRow,MQC_Matrix_Transpose, &
          'Real',Matrix%Storage)
        If (Matrix%Storage.eq.'StorFull') then
          MQC_Matrix_Transpose%MatR = Transpose(Matrix%MatR)
        ElseIf (Matrix%Storage.eq.'StorSymm' .or. Matrix%Storage.eq.'StorDiag') then
          MQC_Matrix_Transpose%MatR = Matrix%MatR
        Else
          Call MQC_Error_A('transposing a triD matrix is not supported yet', 6, &
               'Matrix%Storage', Matrix%Storage )
        EndIf
      ElseIf(MQC_Matrix_HaveInteger(Matrix)) then
        Call MQC_Allocate_Matrix(Matrix%NCol,Matrix%NRow,MQC_Matrix_Transpose, &
          'Integer',Matrix%Storage)
        If (Matrix%Storage.eq.'StorFull') then
          MQC_Matrix_Transpose%MatI = Transpose(Matrix%MatI)
        ElseIf (Matrix%Storage.eq.'StorSymm' .or. Matrix%Storage.eq.'StorDiag') then
          MQC_Matrix_Transpose%MatI = Matrix%MatI
        Else
          Call MQC_Error_A('transposing a triD matrix is not supported yet', 6, &
               'Matrix%Storage', Matrix%Storage )
        EndIf
      ElseIf(MQC_Matrix_HaveComplex(Matrix)) then
        Call MQC_Allocate_Matrix(Matrix%NCol,Matrix%NRow,MQC_Matrix_Transpose, &
          'Complex',Matrix%Storage)
        If (Matrix%Storage.eq.'StorFull') then
          MQC_Matrix_Transpose%MatC = Transpose(Matrix%MatC)
        ElseIf (Matrix%Storage.eq.'StorSymm' .or. Matrix%Storage.eq.'StorDiag') then
          MQC_Matrix_Transpose%MatC = Matrix%MatC
        Else
          Call MQC_Error_A('transposing a triD matrix is not supported yet', 6, &
          'Matrix%Storage', Matrix%Storage )
        EndIf
      EndIf
!
      Return
      End Function MQC_Matrix_Transpose
!
!
!     PROCEDURE MQC_Matrix_Conjugate_Transpose
      Function MQC_Matrix_Conjugate_Transpose(Matrix)
!
!     This function returns the conjugate transpose of an MQC matrx.
!
!     L. M. Thompson, 2016.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::MQC_Matrix_Conjugate_Transpose
      Class(MQC_Matrix),Intent(In)::Matrix
!
      If(MQC_Matrix_HaveReal(Matrix)) then
        Call MQC_Allocate_Matrix(Matrix%NCol,Matrix%NRow,MQC_Matrix_Conjugate_Transpose, &
          'Real',Matrix%Storage)
        If (Matrix%Storage.eq.'StorFull') then
          MQC_Matrix_Conjugate_Transpose%MatR = Transpose(Matrix%MatR)
        ElseIf (Matrix%Storage.eq.'StorSymm' .or. Matrix%Storage.eq.'StorDiag') then
          MQC_Matrix_Conjugate_Transpose%MatR = Matrix%MatR
        Else
          Call MQC_Error_A('Unrecognise matrix storage type in MQC_Matrix_Conjugate_Transpose', 6, &
               'Matrix%Storage', Matrix%Storage )
        EndIf
      ElseIf(MQC_Matrix_HaveInteger(Matrix)) then
        Call MQC_Allocate_Matrix(Matrix%NCol,Matrix%NRow,MQC_Matrix_Conjugate_Transpose, &
          'Integer',Matrix%Storage)
        If (Matrix%Storage.eq.'StorFull') then
          MQC_Matrix_Conjugate_Transpose%MatI = Transpose(Matrix%MatI)
        ElseIf (Matrix%Storage.eq.'StorSymm' .or. Matrix%Storage.eq.'StorDiag') then
          MQC_Matrix_Conjugate_Transpose%MatI = Matrix%MatI
        Else
          Call MQC_Error_A('Unrecognise matrix storage type in MQC_Matrix_Conjugate_Transpose', 6, &
               'Matrix%Storage', Matrix%Storage )
        EndIf
      ElseIf(MQC_Matrix_HaveComplex(Matrix)) then
        Call MQC_Allocate_Matrix(Matrix%NCol,Matrix%NRow,MQC_Matrix_Conjugate_Transpose, &
          'Complex',Matrix%Storage)
        If (Matrix%Storage.eq.'StorFull') then
          MQC_Matrix_Conjugate_Transpose%MatC = Transpose(conjg(Matrix%MatC))
        ElseIf (Matrix%Storage.eq.'StorSymm' .or. Matrix%Storage.eq.'StorDiag') then
          MQC_Matrix_Conjugate_Transpose%MatC = conjg(Matrix%MatC)
        Else
          Call MQC_Error_A('Unrecognise matrix storage type in MQC_Matrix_Conjugate_Transpose', 6, &
               'Matrix%Storage', Matrix%Storage )
        EndIf
      EndIf
!
      Return
      End Function MQC_Matrix_Conjugate_Transpose
!
!
!     PROCEDURE MQC_Matrix_Symmetrize
      Function MQC_Matrix_Symmetrize(Matrix)
!
!     This subroutine is used to symmetrize or antisymmetrize a full
!     packed MQC matrix.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix)::MQC_Matrix_Symmetrize
      Type(MQC_Matrix),Intent(In)::Matrix
      Integer(kind=int64)::I,J
      Real(kind=real64),Parameter::Pt5=0.5
!
      If(.not.Matrix%Storage.eq.'StorFull') Call MQC_Error_A('Matrix not full packed in MQC_Matrix_Symmetrize', 6, &
           'Matrix%Storage', Matrix%Storage )
      If(Matrix%NCol.ne.Matrix%NRow) Call MQC_Error_I('Only square matricies can be symmetrixed &
     &  in MQC_Matrix_Symmetrize', 6, &
     'Matrix%NCol', Matrix%NCol, &
     'Matrix%NRow', Matrix%NRow )
      Call MQC_Set_Matrix2Matrix(MQC_Matrix_Symmetrize,Matrix)
      If(Matrix%Data_Type.eq.'Real') then
        Do I = 1, Matrix%NCol
          Do J = 1, Matrix%NRow
            MQC_Matrix_Symmetrize%MatR(I,J) = Pt5*(Matrix%MatR(I,J)+Matrix%MatR(J,I))
          EndDo
        EndDo
      EndIf
      If(Matrix%Data_Type.eq.'Integer') then
        Call MQC_Matrix_Copy_Int2Real(MQC_Matrix_Symmetrize)
        Do I = 1, Matrix%NCol
          Do J = 1, Matrix%NRow
            MQC_Matrix_Symmetrize%MatR(I,J) = Pt5*(Matrix%MatI(I,J)+Matrix%MatI(J,I))
          EndDo
        EndDo
      EndIf
      If(Matrix%Data_Type.eq.'Complex') then
        Do I = 1, Matrix%NCol
          Do J = 1, Matrix%NRow
            MQC_Matrix_Symmetrize%MatC(I,J) = Pt5*(Matrix%MatC(I,J)+Matrix%MatC(J,I))
          EndDo
        EndDo
      EndIf
!
      Return
      End Function MQC_Matrix_Symmetrize
!
!
!     PROCEDURE MQC_Matrix_Full2Symm
      Subroutine MQC_Matrix_Full2Symm(Matrix)
!
!     This subroutine is used to convert a matrix from full to lower-
!     triangular packing.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::Matrix
      Real(kind=real64),Dimension(:,:),Allocatable::Temp
      Complex(Kind=real64),Dimension(:,:),Allocatable::TempC
      Integer(kind=int64)::I,J,II
!
      If(.not.Matrix%Storage.eq.'StorFull') Call MQC_Error_A('Input matrix must be full &
     &  packed in MQC_Matrix_Full2Symm', 6, &
     'Matrix%Storage', Matrix%Storage )
      If(.not.MQC_Matrix_Test_Symmetric(Matrix).and. &
        .not.MQC_Matrix_Test_Symmetric(Matrix,'antisymmetric').and. &
        .not.MQC_Matrix_Test_Symmetric(Matrix,'hermitian')) &
        Call MQC_Error_L('Input matrix must be symmetric/antisymmetric/hermitian in MQC_Matrix_Full2Symm', 6, &
        'MQC_Matrix_Test_Symmetric(Matrix)', MQC_Matrix_Test_Symmetric(Matrix), &
        "MQC_Matrix_Test_Symmetric(Matrix,'antisymmetric')", MQC_Matrix_Test_Symmetric(Matrix,'antisymmetric'), &
        "MQC_Matrix_Test_Symmetric(Matrix,'hermitian')", MQC_Matrix_Test_Symmetric(Matrix,'hermitian') )
      If(MQC_Matrix_HaveReal(Matrix)) then
        Allocate(Temp(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)))
        Temp = Matrix%MatR
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(Temp,1),Size(Temp,2),Matrix,'Real', &
          'StorSymm')
        Do I = 1, Size(Temp,1)
          II = (I*(I-1))/2
          Do J = 1, I
            Matrix%MatR(II+J,1) = Temp(I,J)
          EndDo
        EndDo
        Deallocate(Temp)
      ElseIf(MQC_Matrix_HaveInteger(Matrix)) then
        Allocate(Temp(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)))
        Temp = Matrix%MatI
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(Temp,1),Size(Temp,2),Matrix,'Integer', &
          'StorSymm')
        Do I = 1, Size(Temp,1)
          II = (I*(I-1))/2
          Do J = 1, I
            Matrix%MatI(II+J,1) = Temp(I,J)
          EndDo
        EndDo
        Deallocate(Temp)
      ElseIf(MQC_Matrix_HaveComplex(Matrix)) then
        Allocate(TempC(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)))
        TempC = Matrix%MatC
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(TempC,1),Size(TempC,2),Matrix,'Complex', &
          'StorSymm')
        Do I = 1, Size(TempC,1)
          II = (I*(I-1))/2
          Do J = 1, I
            Matrix%MatC(II+J,1) = TempC(I,J)
          EndDo
        EndDo
        Deallocate(TempC)
      EndIf
!
      Return
      End Subroutine MQC_Matrix_Full2Symm
!
!
!     PROCEDURE MQC_Matrix_Symm2Full
      Subroutine MQC_Matrix_Symm2Full(Matrix,Option)
!
!     This subroutine is used to convert a matrix from lower-triangular
!     to full packing.
!
!     Options for this subroutine are:
!       'symmetric'      Unpack as if matrix is symmetric
!       'antisymmetric'  Unpack as if matrix is antisymmetric
!       'hermitian'      Unpack as if matrix is hermitian
!       'antihermitian'  Unpack as if matrix is antihermitian
!
!     L. M. Thompson, 2016, 2018.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::Matrix
      Character(len=*),Optional,Intent(In)::Option
      Real(kind=real64),Dimension(:,:),Allocatable::Temp
      Complex(Kind=real64),Dimension(:,:),Allocatable::TempC
      Integer(kind=int64)::I,J,II,N,NTT
      Character(len=64)::myOption
!
      If(Present(Option)) then
        call string_change_case(Option,'l',myOption)
      else
        myOption = 'symmetric'
      endIf
      If(myOption.ne.'symmetric'.and.myOption.ne.'antisymmetric'.and.myOption.ne.'hermitian'.and. &
        myOption.ne.'antihermitian') Call MQC_Error_A('MQC_Matrix_Full2Symm has unrecognised &
     &  transformation option',6,'myOption',myOption)
!
      If(.not.Matrix%Storage.eq.'StorSymm') Call MQC_Error_A('Input matrix must be &
     &  lower-triangular packed in MQC_Matrix_Symm2Full', 6, &
     'Matrix%Storage', Matrix%Storage )
      N = MQC_Matrix_Columns(Matrix)
      NTT = (N*(N+1))/2
      If(Matrix%Data_Type.eq.'Real') then
        Allocate(Temp(NTT,1))
        Temp = Matrix%MatR
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(N,N,Matrix,'Real','StorFull')
        Do I = 1, N
          II = (I*(I-1))/2
          Do J = 1,I
            If(myOption.eq.'antisymmetric'.or.myOption.eq.'antihermitian') then
              Matrix%MatR(J,I) = -Temp(II+J,1)
            Else
              Matrix%MatR(J,I) = Temp(II+J,1)
            EndIf
            Matrix%MatR(I,J) = Temp(II+J,1)
          EndDo
        EndDo
        Deallocate(Temp)
      ElseIf(Matrix%Data_Type.eq.'Integer') then
        Allocate(Temp(NTT,1))
        Temp = Matrix%MatI
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(N,N,Matrix,'Integer','StorFull')
        Do I = 1, N
          II = (I*(I-1))/2
          Do J = 1,I
            If(myOption.eq.'antisymmetric'.or.myOption.eq.'antihermitian') then
              Matrix%MatI(J,I) = -Temp(II+J,1)
            Else
              Matrix%MatI(J,I) = Temp(II+J,1)
            EndIf
            Matrix%MatI(I,J) = Temp(II+J,1)
          EndDo
        EndDo
        Deallocate(Temp)
      ElseIf(Matrix%Data_Type.eq.'Complex') then
        Allocate(TempC(NTT,1))
        TempC = Matrix%MatC
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(N,N,Matrix,'Complex','StorFull')
        Do I = 1, N
          II = (I*(I-1))/2
          Do J = 1,I
            If(myOption.eq.'antisymmetric') then
              Matrix%MatC(J,I) = -TempC(II+J,1)
            ElseIf(myOption.eq.'antihermitian') then
              Matrix%MatC(J,I) = -conjg(TempC(II+J,1))
            ElseIf(myOption.eq.'hermitian') then
              Matrix%MatC(J,I) = conjg(TempC(II+J,1))
            Else
              Matrix%MatC(J,I) = TempC(II+J,1)
            EndIf
            Matrix%MatC(I,J) = TempC(II+J,1)
          EndDo
        EndDo
        Deallocate(TempC)
      EndIf
!
      Return
      End Subroutine MQC_Matrix_Symm2Full
!
!
!     PROCEDURE MQC_Matrix_Full2Diag
      Subroutine MQC_Matrix_Full2Diag(Matrix)
!
!     This subroutine is used to convert a matrix from full to
!     diagonal packing.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::Matrix
      Real(kind=real64),Dimension(:,:),Allocatable::Temp
      Complex(Kind=real64),Dimension(:,:),Allocatable::TempC
      Integer(kind=int64)::I
!
      If(.not.Matrix%Storage.eq.'StorFull') Call MQC_Error_A('Input matrix must be full &
     &  packed in MQC_Matrix_Full2Diag', 6, &
     'Matrix%Storage', Matrix%Storage )
      If(.not.MQC_Matrix_Test_Diagonal(Matrix)) Call MQC_Error_L('Input matrix must be &
     &  diagonal in MQC_Matrix_Full2Diag', 6, &
      'MQC_Matrix_Test_Diagonal(Matrix)', MQC_Matrix_Test_Diagonal(Matrix) )
      If(MQC_Matrix_HaveReal(Matrix)) then
        Allocate(Temp(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)))
        Temp = Matrix%MatR
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(Temp,1),Size(Temp,2),Matrix,'Real', &
          'StorDiag')
        Do I = 1, Min(Size(Temp,1),Size(Temp,2))
          Matrix%MatR(I,1) = Temp(I,I)
        EndDo
        Deallocate(Temp)
      ElseIf(MQC_Matrix_HaveInteger(Matrix)) then
        Allocate(Temp(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)))
        Temp = Matrix%MatI
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(Temp,1),Size(Temp,2),Matrix,'Integer', &
          'StorDiag')
        Do I = 1, Min(Size(Temp,1),Size(Temp,2))
          Matrix%MatI(I,1) = Temp(I,I)
        EndDo
        Deallocate(Temp)
      ElseIf(MQC_Matrix_HaveComplex(Matrix)) then
        Allocate(TempC(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)))
        TempC = Matrix%MatC
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(TempC,1),Size(TempC,2),Matrix,'Complex', &
          'StorDiag')
        Do I = 1, Min(Size(TempC,1),Size(TempC,2))
          Matrix%MatC(I,1) = TempC(I,I)
        EndDo
        Deallocate(TempC)
      Else
        Call MQC_Error_L('Matrix type not defined in MQC_Matrix_Full2Diag', 6, &
             'MQC_Matrix_HaveReal(Matrix)', MQC_Matrix_HaveReal(Matrix) )
      EndIf
!
      Return
      End Subroutine MQC_Matrix_Full2Diag
!
!
!     PROCEDURE MQC_Matrix_Diag2Full
      Subroutine MQC_Matrix_Diag2Full(Matrix)
!
!     This subroutine is used to convert a matrix from diagonal
!     to full packing.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::Matrix
      Real(kind=real64),Dimension(:,:),Allocatable::Temp
      Complex(Kind=real64),Dimension(:,:),Allocatable::TempC
      Integer(kind=int64)::I,N,Columns,Rows
!
      If(.not.Matrix%Storage.eq.'StorDiag') Call MQC_Error_A('Input matrix must be &
     &  diagonal packed in MQC_Matrix_Diag2Full', 6, &
     'Matrix%Storage', Matrix%Storage )
      N = Min(MQC_Matrix_Columns(Matrix),MQC_Matrix_Rows(Matrix))
      Columns = MQC_Matrix_Columns(Matrix)
      Rows = MQC_Matrix_Rows(Matrix)
      If(Matrix%Data_Type.eq.'Real') then
        Allocate(Temp(N,1))
        Temp = Matrix%MatR
        Call MQC_Deallocate_Matrix(Matrix)
        Call Matrix%initialize(Columns,Rows)
        Do I = 1, N
          Matrix%MatR(I,I) = Temp(I,1)
        EndDo
        Deallocate(Temp)
      ElseIf(Matrix%Data_Type.eq.'Integer') then
        Allocate(Temp(N,1))
        Temp = Matrix%MatI
        Call MQC_Deallocate_Matrix(Matrix)
        Call Matrix%initialize(Columns,Rows,0)
        Do I = 1, N
          Matrix%MatI(I,I) = Temp(I,1)
        EndDo
        Deallocate(Temp)
      ElseIf(Matrix%Data_Type.eq.'Complex') then
        Allocate(TempC(N,1))
        TempC = Matrix%MatC
        Call MQC_Deallocate_Matrix(Matrix)
        Call Matrix%initialize(Columns,Rows,cmplx(0,0))
        Do I = 1, N
          Matrix%MatC(I,I) = TempC(I,1)
        EndDo
        Deallocate(TempC)
      EndIf
!
      Return
      End Subroutine MQC_Matrix_Diag2Full
!
!
!     PROCEDURE MQC_Matrix_Symm2Diag
      Subroutine MQC_Matrix_Symm2Diag(Matrix)
!
!     This subroutine is used to convert a matrix from lower-triangular
!     to diagonal packing.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::Matrix
      Real(kind=real64),Dimension(:,:),Allocatable::Temp
      Complex(Kind=real64),Dimension(:,:),Allocatable::TempC
      Integer(kind=int64)::I,II,N,NTT
!
      If(.not.Matrix%Storage.eq.'StorSymm') Call MQC_Error_A('Input matrix must be &
     &  lower-triangular packed in MQC_Matrix_Symm2Diag', 6, &
     'Matrix%Storage', Matrix%Storage )
      If(.not.MQC_Matrix_Test_Diagonal(Matrix)) Call MQC_Error_L('Input matrix must be &
     &  diagonal in MQC_Matrix_Symm2Diag', 6, &
     'MQC_Matrix_Test_Diagonal(Matrix)', MQC_Matrix_Test_Diagonal(Matrix) )
      N = MQC_Matrix_Columns(Matrix)
      NTT = (N*(N+1))/2
      If(Matrix%Data_Type.eq.'Real') then
        Allocate(Temp(NTT,1))
        Temp = Matrix%MatR
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(N,N,Matrix,'Real','StorDiag')
        Do I = 1, N
          II = (I*(I-1))/2
          Matrix%MatR(I,1) = Temp(II+I,1)
        EndDo
        Deallocate(Temp)
      ElseIf(Matrix%Data_Type.eq.'Integer') then
        Allocate(Temp(NTT,1))
        Temp = Matrix%MatI
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(N,N,Matrix,'Integer','StorDiag')
        Do I = 1, N
          II = (I*(I-1))/2
          Matrix%MatI(I,1) = Temp(II+I,1)
        EndDo
        Deallocate(Temp)
      ElseIf(Matrix%Data_Type.eq.'Complex') then
        Allocate(TempC(NTT,1))
        TempC = Matrix%MatC
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(N,N,Matrix,'Complex','StorDiag')
        Do I = 1, N
          II = (I*(I-1))/2
          Matrix%MatC(I,1) = TempC(II+I,1)
        EndDo
        Deallocate(TempC)
      EndIf
!
      Return
      End Subroutine MQC_Matrix_Symm2Diag
!
!
!     PROCEDURE MQC_Matrix_Diag2Symm
      Subroutine MQC_Matrix_Diag2Symm(Matrix)
!
!     This subroutine is used to convert a matrix from diagonal to
!     lower-triangular packing.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(InOut)::Matrix
      Real(kind=real64),Dimension(:,:),Allocatable::Temp
      Complex(Kind=real64),Dimension(:,:),Allocatable::TempC
      Integer(kind=int64)::I,II
!
      If(.not.Matrix%Storage.eq.'StorDiag') Call MQC_Error_A('Input matrix must be &
     &  diagonal packed in MQC_Matrix_Diag2Symm', 6, &
     'Matrix%Storage', Matrix%Storage)
      If(MQC_Matrix_Columns(Matrix).ne.MQC_Matrix_Rows(Matrix)) &
        Call MQC_Error_I('Input matrix must be square in MQC_Matrix_Diag2Symm', 6, &
        'MQC_Matrix_Columns(Matrix)', MQC_Matrix_Columns(Matrix), &
        'MQC_Matrix_Rows(Matrix)', MQC_Matrix_Rows(Matrix) )
      If(Matrix%Data_Type.eq.'Integer') then
        Allocate(Temp(Min(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)),1))
        Temp = Matrix%MatI
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(Temp,1),Size(Temp,1),Matrix,'Integer', &
          'StorSymm')
        Matrix%MatI = 0
        Do I = 1, Size(Temp,1)
          II = (I*(I-1))/2
          Matrix%MatI(II+I,1) = Temp(I,1)
        EndDo
        Deallocate(Temp)
      ElseIf(Matrix%Data_Type.eq.'Real') then
        Allocate(Temp(Min(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)),1))
        Temp = Matrix%MatR
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(Temp,1),Size(Temp,1),Matrix,'Real', &
          'StorSymm')
        Matrix%MatR = 0.0d0
        Do I = 1, Size(Temp,1)
          II = (I*(I-1))/2
          Matrix%MatR(II+I,1) = Temp(I,1)
        EndDo
        Deallocate(Temp)
      ElseIf(Matrix%Data_Type.eq.'Complex') then
        Allocate(TempC(Min(MQC_Matrix_Rows(Matrix),MQC_Matrix_Columns(Matrix)),1))
        TempC = Matrix%MatC
        Call MQC_Deallocate_Matrix(Matrix)
        Call MQC_Allocate_Matrix(Size(TempC,1),Size(TempC,1),Matrix,'Complex', &
          'StorSymm')
        Matrix%MatC = Cmplx(0.0,0.0)
        Do I = 1, Size(TempC,1)
          II = (I*(I-1))/2
          Matrix%MatC(II+I,1) = TempC(I,1)
        EndDo
        Deallocate(TempC)
      EndIf
!
      Return
      End Subroutine MQC_Matrix_Diag2Symm
!
!
!     PROCEDURE MQC_Matrix_Symm2Full_Func
      Function MQC_Matrix_Symm2Full_Func(Matrix) Result(Temp)
!
!     This subroutine is used to convert a matrix from lower-triangular
!     to full packing as a function, for temperory use.
!     E.g., in MQC_MatrixMatrixDotProduct, a Symm stored matrix have to be converted
!     to a full symmetric matrix to work with lapack routines.
!
!     X. Sheng, 2017.
!
!     Variable Declarations.
      Implicit None
      Type(MQC_Matrix),Intent(In)::Matrix
      type(MQC_Matrix)::Temp
      Integer(kind=int64)::I,J,II,N
!
      If(.not.Matrix%Storage.eq.'StorSymm') Call MQC_Error_A('Input matrix must be &
     &  lower-triangular packed in MQC_Matrix_Symm2Full', 6, &
     'Matrix%Storage', Matrix%Storage )
      N = MQC_Matrix_Columns(Matrix)
      Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Temp,Matrix%Data_Type,'StorFull')
      If(Matrix%Data_Type.eq.'Real') then
        Do I = 1, N
          II = (I*(I-1))/2
          Do J = 1,I
            Temp%MatR(I,J) = Matrix%MatR(II+J,1)
            Temp%MatR(J,I) = Matrix%MatR(II+J,1)
          EndDo
        EndDo
      ElseIf(Matrix%Data_Type.eq.'Integer') then
        Do I = 1, N
          II = (I*(I-1))/2
          Do J = 1,I
            Temp%MatI(I,J) = Matrix%MatI(II+J,1)
            Temp%MatI(J,I) = Matrix%MatI(II+J,1)
          EndDo
        EndDo
      ElseIf(Matrix%Data_Type.eq.'Complex') then
        Do I = 1, N
          II = (I*(I-1))/2
          Do J = 1,I
            Temp%MatI(I,J) = Matrix%MatC(II+J,1)
            Temp%MatI(J,I) = Matrix%MatC(II+J,1)
          EndDo
        EndDo
      EndIf
!
      Return
      End Function MQC_Matrix_Symm2Full_Func
!
!
!     PROCEDURE Matrix_Symm2Sq_Integer
      Subroutine Matrix_Symm2Sq_Integer(N,I_Symm,I_Sq)
!
!     This subroutine takes a symmetric matrix (I_Symm) with leading
!     dimension N and expands it to a N x N square matrix (I_Sq).
!
      Implicit None
      Integer(kind=int64),Intent(In)::N
      Integer(kind=int64),Dimension(:),Intent(In)::I_Symm
      Integer(kind=int64),Dimension(N,N),Intent(Out)::I_Sq
!
      Integer(kind=int64)::i,j,k
!
!     Do the work...
!
      k = N*(N+1)/2
      Do j = N,1,-1
        Do i = j,1,-1
          I_Sq(i,j) = I_Symm(k)
          I_Sq(j,i) = I_Symm(k)
          k = k-1
        EndDo
      EndDo
!
      End Subroutine Matrix_Symm2Sq_Integer
!
!
!     PROCEDURE Matrix_Symm2Sq_Real
      Subroutine Matrix_Symm2Sq_Real(N,A_Symm,A_Sq)
!
!     This subroutine takes a symmetric matrix (A_Symm) with leading
!     dimension N and expands it to a N x N square matrix (A_Sq).
!
      Implicit None
      Integer(kind=int64),Intent(In)::N
      Real(kind=real64),Dimension(:),Intent(In)::A_Symm
      Real(kind=real64),Dimension(N,N),Intent(Out)::A_Sq
!
      Integer(kind=int64)::i,j,k
!
!     Do the work...
!
      k = N*(N+1)/2
      Do j = N,1,-1
        Do i = j,1,-1
          A_Sq(i,j) = A_Symm(k)
          A_Sq(j,i) = A_Symm(k)
          k = k-1
        EndDo
      EndDo
!
      End Subroutine Matrix_Symm2Sq_Real
!
!
!     PROCEDURE Matrix_Symm2Sq_Complex
      Subroutine Matrix_Symm2Sq_Complex(N,A_Symm,A_Sq)
!
!     This subroutine takes a symmetric matrix (A_Symm) with leading
!     dimension N and expands it to a N x N square matrix (A_Sq).
!
      Implicit None
      Integer(kind=int64),Intent(In)::N
      Complex(Kind=real64),Dimension(:),Intent(In)::A_Symm
      Complex(Kind=real64),Dimension(N,N),Intent(Out)::A_Sq
!
      Integer(kind=int64)::i,j,k
!
!     Do the work...
!
      k = N*(N+1)/2
      Do j = N,1,-1
        Do i = j,1,-1
          A_Sq(i,j) = A_Symm(k)
          A_Sq(j,i) = A_Symm(k)
          k = k-1
        EndDo
      EndDo
!
      End Subroutine Matrix_Symm2Sq_Complex
!
!
!     PROCEDURE MQC_Vector2DiagMatrix
      Function MQC_Vector2DiagMatrix(vector) Result(matrix)
!
!     This subroutine is used to convert a vector to a Diagonal matrix
!
!     X. Sheng, 2017.
!
      Implicit None
      Type(MQC_Matrix)::matrix
      type(MQC_Vector),Intent(In)::vector
!
      if (vector%Data_type.eq.'Real') then
        call mqc_allocate_matrix(vector%length,vector%length,matrix,'Real','StorDiag')
        matrix%matR(:,1) = vector%vecR
      elseif (vector%Data_type.eq.'Integer') then
        call mqc_allocate_matrix(vector%length,vector%length,matrix,'Integer','StorDiag')
        matrix%matI(:,1) = vector%vecI
      elseif (vector%Data_type.eq.'Complex') then
        call mqc_allocate_matrix(vector%length,vector%length,matrix,'Complex','StorDiag')
        matrix%matC(:,1) = vector%vecC
      endif
!
      End Function MQC_Vector2DiagMatrix
!
!
!     PROCEDURE MQC_MatrixMatrixSum
      Function MQC_MatrixMatrixSum(MA,MB) Result(MC)
!
!     This function returns the sum of two MQC_Matrix matrices.
!
!     H. P. Hratchian, 2017.
!     L. M. Thompson, 2018.
!
!
      implicit none
      type(MQC_Matrix),Intent(In)::MA,MB
      type(MQC_Matrix)::MC,Tmp1,Tmp2
!
!
      Tmp1 = MA
      Tmp2 = MB
      select case(Tmp1%storage)
      case('StorFull')
        select case(Tmp2%storage)
        case('StorSymm')
          call MQC_Matrix_Symm2Full(Tmp2)
        case('StorDiag')
          call MQC_Matrix_Diag2Full(Tmp2)
        case('StorFull')
          !-Empty Case-!
        case default
          call mqc_error_A('UNKOWN storage type in MQC_MatrixMatrixSum.', 6, &
               'Tmp2%storage', Tmp2%storage) 
        end select
      case('StorSymm')
        select case(Tmp2%storage)
        case('StorFull')
          Call MQC_Matrix_Symm2Full(Tmp1)
        case('StorDiag')
          Call MQC_Matrix_Diag2Symm(Tmp2)
        case('StorSymm')
          !-Empty Case-!
        case default
          call mqc_error_A('UNKOWN storage type in MQC_MatrixMatrixSum.', 6, &
               'Tmp2%storage', Tmp2%storage) 
        end select
      case('StorDiag')
        select case(Tmp2%storage)
        case('StorFull')
          Call MQC_Matrix_Diag2Full(Tmp1)
        case('StorSymm')
          Call MQC_Matrix_Diag2Symm(Tmp1)
        case('StorDiag')
          !-Empty Case-!
        case default
          call mqc_error_A('UNKOWN storage type in MQC_MatrixMatrixSum.', 6, &
               'Tmp2%storage', Tmp2%storage) 
        end select
      case default
        call mqc_error_A('UNKNOWN storage type in MQC_MatrixMatrixSum.', 6, &
             'Tmp1%storage', Tmp1%storage) 
      end select
!
      select case(Tmp1%data_type)
      case('Real')
        select case(Tmp2%data_type)
        case('Real')
          MC = Tmp2
          MC%MatR = Tmp1%MatR + Tmp2%MatR
        case('Integer')
          MC = Tmp1
          MC%MatR = Tmp1%MatR + Tmp2%MatI
        case('Complex')
          MC = Tmp2
          MC%MatC = Tmp1%MatR + Tmp2%MatC
        case default
          call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixSum.', 6, &
               'Tmp2%data_type', Tmp2%data_type) 
        end select
      case('Integer')
        select case(Tmp2%data_type)
        case('Real')
          MC = Tmp2
          MC%MatR = Tmp1%MatI + Tmp2%MatR
        case('Integer')
          MC = Tmp2
          MC%MatI = Tmp1%MatI + Tmp2%MatI
        case('Complex')
          MC = Tmp2
          MC%MatC = Tmp1%MatI + Tmp2%MatC
        case default
          call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixSum.', 6, &
               'Tmp2%data_type', Tmp2%data_type) 
        end select
      case('Complex')
        select case(Tmp2%data_type)
        case('Real')
          MC = Tmp1
          MC%MatC = Tmp1%MatC + Tmp2%MatR
        case('Integer')
          MC = Tmp1
          MC%MatC = Tmp1%MatC + Tmp2%MatI
        case('Complex')
          MC = Tmp2
          MC%MatC = Tmp1%MatC + Tmp2%MatC
        case default
          call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixSum.', 6, &
               'Tmp2%data_type', Tmp2%data_type) 
        end select
      case default
        call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixSum.', 6, &
             'Tmp1%data_type', Tmp1%data_type) 
      end select
!
      return
      end function MQC_MatrixMatrixSum
!
!
!     PROCEDURE MQC_MatrixMatrixSubtract
      Function MQC_MatrixMatrixSubtract(MA,MB) Result(MC)
!
!     This function returns the difference of two MQC_Matrix matrices.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      type(MQC_Matrix),Intent(In)::MA,MB
      type(MQC_Matrix)::MC,Tmp1,Tmp2

      Tmp1 = MA
      Tmp2 = MB
      select case(Tmp1%storage)
      case('StorFull')
        select case(Tmp2%storage)
        case('StorSymm')
          call MQC_Matrix_Symm2Full(Tmp2)
        case('StorDiag')
          call MQC_Matrix_Diag2Full(Tmp2)
        case('StorFull')
          !-Empty Case-!
        case default
          call mqc_error_A('UNKNOWN storage type in MQC_MatrixMatrixSubtract.', 6, &
               'Tmp2%storage', Tmp2%storage)
        end select
      case('StorSymm')
        select case(Tmp2%storage)
        case('StorFull')
          call MQC_Matrix_Symm2Full(Tmp1)
        case('StorDiag')
          call MQC_Matrix_Diag2Symm(Tmp2)
        case('StorSymm')
          !-Empty Case-!
        case default
          call mqc_error_A('UNKNOWN storage type in MQC_MatrixMatrixSubtract.', 6, &
               'Tmp2%storage', Tmp2%storage)
        end select
      case('StorDiag')
        select case(Tmp2%storage)
        case('StorFull')
          call MQC_Matrix_Diag2Full(Tmp1)
        case('StorSymm')
          call MQC_Matrix_Diag2Symm(Tmp1)
        case('StorDiag')
          !-Empty Case-!
        case default
          call mqc_error_A('UNKNOWN storage type in MQC_MatrixMatrixSubtract.', 6, &
               'Tmp2%storage', Tmp2%storage)
        end select
      case default
        call mqc_error_A('UNKNOWN storage type in MQC_MatrixMatrixSubtract.', 6, &
             'Tmp1%storage', Tmp1%storage)
      end select
!     
!     Do the work.
!
      select case(Tmp1%data_type)
      case('Real')
        select case(Tmp2%data_type)
        case('Real')
          MC = Tmp2
          MC%MatR = Tmp1%MatR - Tmp2%MatR
        case('Integer')
          MC = Tmp1
          MC%MatR = Tmp1%MatR - Tmp2%MatI
        case('Complex')
          MC = Tmp2
          MC%MatC = Tmp1%MatR - Tmp2%MatC
        case default
          call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixSubtract.', 6, &
               'Tmp2%data_type', Tmp2%data_type)
        end select
      case('Integer')
        select case(Tmp2%data_type)
        case('Real')
          MC = Tmp2
          MC%MatR = Tmp1%MatI - Tmp2%MatR
        case('Integer')
          MC = Tmp2
          MC%MatI = Tmp1%MatI - Tmp2%MatI
        case('Complex')
          MC = Tmp2
          MC%MatC = Tmp1%MatI - Tmp2%MatC
        case default
          call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixSubtract.', 6, &
               'Tmp2%data_type', Tmp2%data_type)
        end select
      case('Complex')
        select case(Tmp2%data_type)
        case('Real')
          MC = Tmp1
          MC%MatC = Tmp1%MatC - Tmp2%MatR
        case('Integer')
          MC = Tmp1
          MC%MatC = Tmp1%MatC - Tmp2%MatI
        case('Complex')
          MC = Tmp2
          MC%MatC = Tmp1%MatC - Tmp2%MatC
        case default
          call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixSubtract.', 6, &
               'Tmp2%data_type', Tmp2%data_type)
        end select
      case default
        call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixSubtract.', 6, &
             'Tmp1%data_type', Tmp1%data_type)
      end select
!
      return
      end function MQC_MatrixMatrixSubtract
!
!
!     PROCEDURE MQC_MatrixMatrixProduct
      Function MQC_MatrixMatrixProduct(MA,MB) Result(MC)
!
!     This function returns the (element-by-element) product of two MQC_Matrix
!     matrices.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      type(MQC_Matrix),Intent(In)::MA,MB
      type(MQC_Matrix)::MC
      type(MQC_Matrix)::MAreal,MBreal,Tmp1,Tmp2
      character(len=64)::matrixType
!
!     For now, addition of two MQC matrices is only allowed for same-type and
!     same-storage pairs. Start by ensuring this limitation is honored.
!
      if(TRIM(MA%data_type) /= TRIM(MB%data_type))  &
        call mqc_error_a('MQC_MatrixMatrixProduct: Data types do NOT match.', 6, &
        'TRIM(MA%data_type)', TRIM(MA%data_type), &
        'TRIM(MB%data_type)', TRIM(MB%data_type) )
      if(TRIM(MA%storage) /= TRIM(MB%storage))  &
        call mqc_error_a('MQC_MatrixMatrixProduct: Storage forms do NOT match.', 6, &
        'TRIM(MA%storage)', TRIM(MA%storage), &
        'TRIM(MB%storage)', TRIM(MB%storage) )
      MC = MA
!
!     Do the work.
!
      call String_Change_Case(MC%Data_Type,'u',matrixType)
      select case(TRIM(matrixType))
      case('REAL')
        MC%MatR = MC%MatR * MB%MatR
      case('INTEGER')
        MC%MatI = MC%MatI * MB%MatI
      case('COMPLEX')
        MC%MatC = MC%MatC * MB%MatC
      case default
        call mqc_error_A('UNKNOWN matrix type in MQC_MatrixMatrixProduct.', 6, &
             'TRIM(matrixType)', TRIM(matrixType) )
      end select
!
      return
      end function MQC_MatrixMatrixProduct
!
!
!     PROCEDURE MQC_MatrixMatrixDotProduct
      Function MQC_MatrixMatrixDotProduct(MA,MB) Result(MC)
!
!     This subroutine is used to calculate matrix multiplication for all types of MQC_matrix
!     MC = MA .dot. MB
!
!     X. Sheng   2017.2.
!
      Implicit None
      Type(MQC_Matrix),Intent(In)::MA,MB
      Type(MQC_Matrix)::MC
      Type(MQC_Matrix)::MAreal,MBreal,Tmp1,Tmp2
      integer(kind=int64)::i,j

      If (MA%NCol /= MB%NRow) call MQC_Error_I('The two matrices are not conformable for multiplication', 6, &
           'MA%NCol', MA%NCol, &
           'MB%NRow', MB%NRow )

      If (MQC_Matrix_HaveInteger(MA)) MAreal = MQC_Cast_Real(MA)
      If (MQC_Matrix_HaveInteger(MB)) MBreal = MQC_Cast_Real(MB)
      If (MQC_Matrix_HaveReal(MA) .and. MQC_Matrix_HaveReal(MB)) then

        If (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          MC%MatR = MatMul(MA%MatR,MB%MatR)
        ElseIf (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          mc = matmul(MA%MatR,Tmp1%MatR)
!          Call DSymm('r','L',MA%NRow,MB%NCol,1.0,Tmp1%MatR,MB%NRow,MA%MatR,MA%NRow,0, &
!            MC%MatR,MC%NRow)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          mc = matmul(Tmp1%MatR,MB%MatR)
!          Call DSymm('l','L',MA%NRow,MB%NCol,1.0,MA%MatR,MA%NRow,Tmp1%MatR,MB%NRow,0, &
!            MC%MatR,MC%NRow)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          mc = matmul(Tmp1%MatR,Tmp2%MatR)
!          Call DSymm('r','L',MA%nrow,MB%ncol,1.0,Tmp2%MatR,MB%NRow,Tmp1%MatR,MA%NRow,0, &
!            MC%MatR,MC%NRow)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          do i = 1, MB%NRow
            MC%matR(i,:) = MA%matR(i,1) * MB%matR(i,:)
          enddo
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          do i = 1, MB%NRow
            do j = 1, MB%NCol
              MC%matR(i,j) = MA%matR(i,1) * MB%matR(symIndexHash(i,j),1)
            enddo
          enddo
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorDiag')
          MC%matR = MA%matR * MB%matR
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorFull')
          do i = 1, MA%NCol
            MC%matR(:,i) = MA%matR(:,i) * MB%matR(i,1)
          enddo
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorFull')
          do i = 1, MA%NRow
            do j = 1, MA%NCol
              MC%matR(i,j) = MB%matR(i,1) * MA%matR(symIndexHash(i,j),1)
            enddo
          enddo
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorDiag') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorDiag')
          MC%matR = MB%matR * MA%matR
        EndIf
      ElseIf (MQC_Matrix_HaveReal(MA) .and. MQC_Matrix_HaveInteger(MB)) then
        If (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MA%NRow,MBreal%NCol,MC,'Real','StorFull')
          MC%MatR = MatMul(MA%MatR,MBreal%MatR)
        ElseIf (MA%Storage.eq.'StorFull' .and. MBreal%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MA%NRow,MBreal%NCol,MC,'Real','StorFull')
          Call MQC_Matrix_Symm2Full(MBreal)
          Call DSymm('r','L',MA%NRow,MBreal%NCol,1.0,MBreal%MatR,MBreal%NRow,MA%MatR, &
            MA%NRow,0,MC%MatR,MC%NRow)
          Call MQC_Matrix_Full2Symm(MBreal)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MBreal%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MA%NRow,MBreal%NCol,MC,'Real','StorFull')
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call DSymm('l','L',MA%NRow,MBreal%NCol,1.0,Tmp1%MatR,MA%NRow,MBreal%MatR, &
            MBreal%NRow,0,MC%MatR,MC%NRow)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MBreal%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MA%NRow,MBreal%NCol,MC,'Real','StorFull')
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(MBreal)
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call DSymm('r','L',MA%NRow,MBreal%NCol,1.0,MBreal%MatR,MBreal%NRow,Tmp1%MatR, &
            MA%NRow,0,MC%MatR,MC%NRow)
          Call MQC_Matrix_Full2Symm(MBreal)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          do i = 1, MB%NRow
            MC%matR(i,:) = MA%matR(i,1) * MB%matI(i,:)
          enddo
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          do i = 1, MB%NRow
            do j = 1, MB%NCol
              MC%matR(i,j) = MA%matR(i,1) * MB%matI(symIndexHash(i,j),1)
            enddo
          enddo
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorDiag')
          MC%matR = MA%matR * MB%matI

        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorFull')
          do i = 1, MA%NRow
            MC%matR(i,:) = MB%matR(i,1) * MA%matI(i,:)
          enddo
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorFull')
          do i = 1, MA%NRow
            do j = 1, MA%NCol
              MC%matR(i,j) = MB%matR(i,1) * MA%matI(symIndexHash(i,j),1)
            enddo
          enddo
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorDiag') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorDiag')
          MC%matR = MB%matR * MA%matI
        EndIf
      ElseIf (MQC_Matrix_HaveInteger(MA) .and. MQC_Matrix_HaveReal(MB)) then
        If (MAreal%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MAreal%nrow,MB%ncol,MC,'Real','StorFull')
          MC%MatR = MatMul(MAreal%MatR,MB%MatR)
        ElseIf (MAreal%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MAreal%nrow,MB%ncol,MC,'Real','StorFull')
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call DSymm('r','L',MAreal%NRow,MB%NCol,1.0,Tmp1%MatR,MB%NRow,MAreal%MatR, &
            MAreal%NRow,0,MC%MatR,MC%NRow)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MAreal%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MAreal%nrow,MB%ncol,MC,'Real','StorFull')
          Call MQC_Matrix_Symm2Full(MAreal)
          Call DSymm('l','L',MAreal%NRow,MB%NCol,1.0,MAreal%MatR,MAreal%NRow,MB%MatR, &
            MB%NRow,0,MC%MatR,MC%NRow)
          Call MQC_Matrix_Full2Symm(MAreal)
        ElseIf (MAreal%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MAreal%nrow,MB%ncol,MC,'Real','StorFull')
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(MAreal)
          Call DSymm('r','L',MAreal%nrow,MB%ncol,1.0,Tmp1%MatR,MB%NRow,MAreal%MatR, &
            MAreal%NRow,0,MC%MatR,MC%NRow)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Matrix_Full2Symm(MAreal)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          do i = 1, MB%NRow
            MC%matR(i,:) = MA%matI(i,1) * MB%matR(i,:)
          enddo
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorFull')
          do i = 1, MB%NRow
            do j = 1, MB%NCol
              MC%matR(i,j) = MA%matI(i,1) * MB%matR(symIndexHash(i,j),1)
            enddo
          enddo
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Real','StorDiag')
          MC%matR = MA%matI * MB%matR

        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorFull')
          do i = 1, MA%NRow
            MC%matR(i,:) = MB%matI(i,1) * MA%matR(i,:)
          enddo
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorFull')
          do i = 1, MA%NRow
            do j = 1, MA%NCol
              MC%matR(i,j) = MB%matI(i,1) * MA%matR(symIndexHash(i,j),1)
            enddo
          enddo
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorDiag') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Real','StorDiag')
          MC%matR = MB%matI * MA%matR
        EndIf
      ElseIf (MQC_Matrix_HaveInteger(MA) .and. MQC_Matrix_HaveInteger(MB)) then
        Call MQC_Allocate_Matrix(MA%nrow,MB%ncol,MC,'Integer','StorFull')
        If (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          MC%MatI = MatMul(MA%MatI,MB%MatI)
        ElseIf (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatI = MatMul(MA%MatI,Tmp1%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
          !call dsymm('r','L',MA%nrow,MB%ncol,1.0,MQC_Matrix_Symm2Full_func(MB)%mati,&
           !MB%nrow,MA%mati,MA%nrow,0,MC%mati,MC%nrow)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatI = MatMul(Tmp1%MatI,MB%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
          !call dsymm('l','L',MA%nrow,MB%ncol,1.0,MQC_Matrix_Symm2Full_func(MA),&
           !MA%nrow,MB%mati,MB%nrow,0,MC%mati,MC%nrow)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatI = MatMul(Tmp2%MatI,Tmp1%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
          !call dsymm('r','L',MA%nrow,MB%ncol,1.0,MQC_Matrix_Symm2Full_func(MB),&
           !MB%nrow,MQC_Matrix_Symm2Full_func(MA),MA%nrow,0,MC%mati,MC%nrow)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Integer','StorFull')
          do i = 1, MB%NRow
            MC%matI(i,:) = MA%matI(i,1) * MB%matI(i,:)
          enddo
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Integer','StorFull')
          do i = 1, MB%NRow
            do j = 1, MB%NCol
              MC%matI(i,j) = MA%matI(i,1) * MB%matI(symIndexHash(i,j),1)
            enddo
          enddo
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Call MQC_Allocate_Matrix(MA%NRow,MB%NCol,MC,'Integer','StorDiag')
          MC%matI = MA%matI * MB%matI

        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorFull') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Integer','StorFull')
          do i = 1, MA%NRow
            MC%matI(i,:) = MB%matI(i,1) * MA%matI(i,:)
          enddo
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorSymm') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Integer','StorFull')
          do i = 1, MA%NRow
            do j = 1, MA%NCol
              MC%matI(i,j) = MB%matI(i,1) * MA%matI(symIndexHash(i,j),1)
            enddo
          enddo
        ElseIf (MB%Storage.eq.'StorDiag' .and. MA%Storage.eq.'StorDiag') then
          Call MQC_Allocate_Matrix(MB%NRow,MA%NCol,MC,'Integer','StorDiag')
          MC%matI = MB%matI * MA%matI
        EndIf
      ElseIf (MQC_Matrix_HaveComplex(MA) .and. MQC_Matrix_HaveInteger(MB)) then
        Call MQC_Allocate_Matrix(MA%nrow,MB%ncol,MC,'Complex','StorFull')
        If (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          MC%MatC = MatMul(MA%MatC,MB%MatI)
        ElseIf (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(MA%MatC,Tmp1%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorFull'.and.MB%Storage.eq.'StorDiag') then
          Tmp1 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(MA%MatC,Tmp1%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatC,MB%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatC,MB%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        EndIf
      ElseIf (MQC_Matrix_HaveComplex(MA) .and. MQC_Matrix_HaveReal(MB)) then
        Call MQC_Allocate_Matrix(MA%nrow,MB%ncol,MC,'Complex','StorFull')
        If (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          MC%MatC = MatMul(MA%MatC,MB%MatR)
        ElseIf (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(MA%MatC,Tmp1%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorFull'.and.MB%Storage.eq.'StorDiag') then
          Tmp1 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(MA%MatC,Tmp1%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatC,MB%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatC,MB%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        EndIf
      ElseIf (MQC_Matrix_HaveInteger(MA) .and. MQC_Matrix_HaveComplex(MB)) then
        Call MQC_Allocate_Matrix(MA%nrow,MB%ncol,MC,'Complex','StorFull')
        If (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          MC%MatC = MatMul(MA%MatI,MB%MatC)
        ElseIf (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(MA%MatI,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorFull'.and.MB%Storage.eq.'StorDiag') then
          Tmp1 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(MA%MatI,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatI,MB%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatI,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatI,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatI,MB%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatI,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatI,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        EndIf
      ElseIf (MQC_Matrix_HaveReal(MA) .and. MQC_Matrix_HaveComplex(MB)) then
        Call MQC_Allocate_Matrix(MA%nrow,MB%ncol,MC,'Complex','StorFull')
        If (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          MC%MatC = MatMul(MA%MatR,MB%MatC)
        ElseIf (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(MA%MatR,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorFull'.and.MB%Storage.eq.'StorDiag') then
          Tmp1 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(MA%MatR,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatR,MB%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatR,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatR,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatR,MB%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatR,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatR,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        EndIf
      ElseIf (MQC_Matrix_HaveComplex(MA) .and. MQC_Matrix_HaveComplex(MB)) then
        Call MQC_Allocate_Matrix(MA%nrow,MB%ncol,MC,'Complex','StorFull')
        If (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorFull') then
          MC%MatC = MatMul(MA%MatC,MB%MatC)
        ElseIf (MA%Storage.eq.'StorFull' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(MA%MatC,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorFull'.and.MB%Storage.eq.'StorDiag') then
          Tmp1 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(MA%MatC,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatC,MB%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorSymm' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorFull') then
          Tmp1 = MA
          Call MQC_Matrix_Diag2Full(Tmp1)
          MC%MatC = MatMul(Tmp1%MatC,MB%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Symm2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        ElseIf (MA%Storage.eq.'StorDiag' .and. MB%Storage.eq.'StorDiag') then
          Tmp1 = MA
          Tmp2 = MB
          Call MQC_Matrix_Diag2Full(Tmp1)
          Call MQC_Matrix_Diag2Full(Tmp2)
          MC%MatC = MatMul(Tmp1%MatC,Tmp2%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
          Call MQC_Deallocate_Matrix(Tmp2)
        EndIf
      EndIf

      End Function MQC_MatrixMatrixDotProduct
!
!
!     PROCEDURE MQC_MatrixVectorDotProduct
      Function MQC_MatrixVectorDotProduct(MA,VB) Result(VC)
!
!     This subroutine is used to calculate matrix vector multiplication for
!     all types of MQC_matrix and MQC_vector
!                         VC = MA .dot. VB
!     Note that VB%Column has to be .true. and VC%column will be .true.
!
!     X. Sheng   2017.2.
!
      Implicit None
      Type(MQC_Matrix),Intent(In)::MA
      Type(MQC_Vector),Intent(In)::VB
      Type(MQC_Vector)::VC,VBreal
      Type(MQC_Matrix)::MAreal,Tmp1
      Integer(kind=int64),Allocatable::Result(:,:)

      If (.not. VB%Column .or. MA%NCol /= VB%Length) then
        Call MQC_Error_I('Matrix and vector are not conformable for multiplication', 6, &
             'MA%NCol', MA%NCol, &
             'VB%Length', VB%Length)
      EndIf

      If (MQC_Matrix_HaveInteger(MA)) MAreal = MQC_Cast_Real(MA)
      If (MQC_Vector_HaveInteger(VB)) VBreal = MQC_Cast_Real(VB)
      If (MQC_Matrix_HaveReal(MA) .and. MQC_Vector_HaveReal(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Real')
        If (MA%Storage.eq.'StorFull') then
          Call DGemV('N',MA%NRow,MA%NCol,1.0,MA%MatR,MA%NRow,VB%VecR,1,0.0,VC%VecR,1)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call DSymV('L',MA%NRow,1.0,Tmp1%MatR,MA%NRow,VB%VecR,1,0.0,VC%VecR,1)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecR = MA%matR(:,1) * VB%vecR
        EndIf
      ElseIf (MQC_Matrix_HaveReal(MA) .and. MQC_Vector_HaveInteger(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Real')
        If (MA%Storage.eq.'StorFull') then
          Call DGemV('N',MA%NRow,MA%NCol,1.0,MA%MatR,MA%NRow,VBreal%VecR,1,0.0,VC%VecR,1)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call DSymV('L',MA%NRow,1.0,Tmp1%MatR,MA%NRow,VBreal%VecR,1,0.0,VC%VecR,1)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecR = MA%matR(:,1) * VB%vecI
        EndIf
      ElseIf (MQC_Matrix_HaveInteger(MA) .and. MQC_Vector_HaveReal(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Real')
        If (MA%Storage.eq.'StorFull') then
          Call DGemV('N',MA%NRow,MA%NCol,1.0,MAreal%MatR,MA%NRow,VB%VecR,1,0.0,VC%VecR,1)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Call MQC_Matrix_Symm2Full(MAreal)
          Call DSymV('L',MA%NRow,1.0,MAreal%MatR,MA%NRow,VB%VecR,1,0.0,VC%VecR,1)
          Call MQC_Matrix_Full2Symm(MAreal)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecR = MA%matI(:,1) * VB%vecR
        EndIf
      ElseIf (MQC_Matrix_HaveInteger(MA) .and. MQC_Vector_HaveInteger(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Integer')
        Allocate(Result(MA%NRow,1))
        If (MA%Storage.eq.'StorFull') then
          Result = MatMul(MA%MatI,Reshape(VB%VecI,[MA%NCol,1]))
          VC%VecI = Result(:,1)
          !call dgemv('N',MA%nrow,MA%ncol,1.0,MA%mati,MA%nrow,VB%veci,1,0.0,VC%veci,1)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          Result = MatMul(Tmp1%MatI,Reshape(VB%VecI,[MA%NCol,1]))
          Call MQC_Deallocate_Matrix(Tmp1)
          VC%VecI = Result(:,1)
          !call dsymv('L',MA%nrow,1.0,MQC_Matrix_Symm2Full_func(MA), &
           !MA%nrow,VB%veci,1,0.0,VC%veci,1)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecI = MA%matI(:,1) * VB%vecI
        EndIf
      ElseIf (MQC_Matrix_HaveInteger(MA) .and. MQC_Vector_HaveComplex(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Complex')
        If (MA%Storage.eq.'StorFull') then
          VC%VecC = MatMul(MA%MatI,VB%VecC)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(Tmp1%MatI,VB%VecC)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecC = MA%matI(:,1) * VB%vecC
        EndIf
      ElseIf (MQC_Matrix_HaveReal(MA) .and. MQC_Vector_HaveComplex(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Complex')
        If (MA%Storage.eq.'StorFull') then
          VC%VecC = MatMul(MA%MatR,VB%VecC)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(Tmp1%MatR,VB%VecC)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecC = MA%matR(:,1) * VB%vecC
        EndIf
      ElseIf (MQC_Matrix_HaveComplex(MA) .and. MQC_Vector_HaveInteger(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Complex')
        If (MA%Storage.eq.'StorFull') then
          VC%VecC = MatMul(MA%MatC,VB%VecI)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(Tmp1%MatC,VB%VecI)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecC = MA%matC(:,1) * VB%vecI
        EndIf
      ElseIf (MQC_Matrix_HaveComplex(MA) .and. MQC_Vector_HaveReal(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Complex')
        If (MA%Storage.eq.'StorFull') then
          VC%VecC = MatMul(MA%MatC,VB%VecR)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(Tmp1%MatC,VB%VecR)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecC = MA%matC(:,1) * VB%vecR
        EndIf
      ElseIf (MQC_Matrix_HaveComplex(MA) .and. MQC_Vector_HaveComplex(VB)) then
        Call MQC_Allocate_Vector(MA%NRow,VC,'Complex')
        If (MA%Storage.eq.'StorFull') then
          VC%VecC = MatMul(MA%MatC,VB%VecC)
        ElseIf (MA%Storage.eq.'StorSymm') then
          Tmp1 = MA
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(Tmp1%MatC,VB%VecC)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MA%Storage.eq.'StorDiag') then
          VC%vecC = MA%matC(:,1) * VB%vecC
        EndIf
      EndIf

      End Function MQC_MatrixVectorDotProduct
!
!
!     PROCEDURE MQC_VectorMatrixDotProduct
      Function MQC_VectormatrixDotProduct(VA,MB) Result(VC)
!
!     This subroutine is used to calculate vector matrix multiplication for
!     all types of MQC_matrix and MQC_vector
!                         VC = VA .dot. MB
!     Note that VA%Column has to be .false. and VC%column will be .false.
!
!     X. Sheng   2017.2.
!
      Implicit None
      Type(MQC_Matrix),Intent(In)::MB
      Type(MQC_Vector),Intent(In)::VA
      Type(MQC_Vector)::VC,VAreal
      Type(MQC_Matrix)::MBreal,Tmp1
      Integer(kind=int64),Allocatable::Result(:,:)

      If (VA%Column .or. VA%Length /= MB%NRow) then
        Call MQC_Error_i('Matrix and vector are not conformable for multiplication', 6, &
             'VA%Length', VA%Length, &
             'MB%NRow', MB%NRow )
      EndIf

      If (.not. MQC_Matrix_HaveReal(MB)) MBreal = MQC_Cast_Real(MB)
      If (.not. MQC_Vector_HaveReal(VA)) VAreal = MQC_Cast_Real(VA)
      If (MQC_Matrix_HaveReal(MB) .and. MQC_Vector_HaveReal(VA)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Real')
        VC = MQC_Vector_Transpose(VC)
        If (MB%Storage.eq.'StorFull') then
          Call DGemV('T',MB%NRow,MB%NCol,1.0,MB%MatR,MB%NRow,VA%VecR,1,0.0,VC%VecR,1)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call DSymV('L',MB%NRow,1.0,Tmp1%MatR,MB%NRow,VA%VecR,1,0.0,VC%VecR,1)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MB%Storage.eq.'StorDiag') then
          VC%vecR = MB%matR(:,1) * VA%vecR
        EndIf
      ElseIf (MQC_Matrix_HaveInteger(MB) .and. MQC_Vector_HaveReal(VA)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Real')
        VC = MQC_Vector_Transpose(VC)
        If (MB%Storage.eq.'StorFull') then
          Call DGemV('T',MB%NRow,MB%NCol,1.0,MBreal%MatR,MB%NRow,VA%VecR,1,0.0,VC%VecR,1)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Call MQC_Matrix_Symm2Full(MBreal)
          Call DSymV('L',MB%NRow,1.0,MBreal%MatR,MB%NRow,VA%VecR,1,0.0,VC%VecR,1)
          Call MQC_Matrix_Full2Symm(MBreal)
        elseif (MB%Storage.eq.'StorDiag') then
          VC%vecR = MB%matI(:,1) * VA%vecR
        EndIf
      ElseIf (MQC_Matrix_HaveReal(MB) .and. MQC_Vector_HaveInteger(VA)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Real')
        VC = MQC_Vector_Transpose(VC)
        If (MB%Storage.eq.'StorFull') then
          Call DGemV('T',MB%nrow,MB%ncol,1.0,MB%matr,MB%nrow,VAreal%vecr,1,0.0,VC%vecr,1)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Call DSymV('L',MB%NRow,1.0,Tmp1%MatR,MB%NRow,VAreal%VecR,1,0.0,VC%VecR,1)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MB%Storage.eq.'StorDiag') then
          VC%vecR = MB%matR(:,1) * VA%vecI
        EndIf
      ElseIf (MQC_Matrix_HaveInteger(MB) .and. MQC_Vector_HaveInteger(VA)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Integer')
        VC = MQC_Vector_Transpose(VC)
        Allocate(Result(1,MB%NCol))
        If (MB%Storage.eq.'StorFull') then
          Result = MatMul(Reshape(VA%VecI,[1,MB%NRow]),MB%MatI)
          VC%VecI = Result(1,:)
          !call dgemv('T',MB%nrow,MB%ncol,1.0,MB%mati,MB%nrow,VA%veci,1,0.0,VC%veci,1)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          Result = MatMul(Reshape(VA%VecI,[1,MB%NRow]),Tmp1%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
          VC%VecI = Result(1,:)
          !call dsymv('L',MB%nrow,1.0,MQC_matrix_Symm2Full_func(MB), &
           !MB%nrow,VA%veci,1,0.0,VC%veci,1)
        elseif (MB%Storage.eq.'StorDiag') then
         VC%vecI = MB%matI(:,1) * VA%vecI
        EndIf
      ElseIf (MQC_Vector_HaveInteger(VA) .and. MQC_Matrix_HaveComplex(MB)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Complex')
        If (MB%Storage.eq.'StorFull') then
          VC%VecC = MatMul(VA%VecI,MB%MatC)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(VA%VecI,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MB%Storage.eq.'StorDiag') then
          VC%vecC = VA%vecI * MB%matC(:,1) 
        EndIf
      ElseIf (MQC_Vector_HaveReal(VA) .and. MQC_Matrix_HaveComplex(MB)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Complex')
        If (MB%Storage.eq.'StorFull') then
          VC%VecC = MatMul(VA%VecR,MB%MatC)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(VA%VecR,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MB%Storage.eq.'StorDiag') then
          VC%vecC = VA%VecR * MB%matC(:,1)
        EndIf
      ElseIf (MQC_Vector_HaveComplex(VA) .and. MQC_Matrix_HaveInteger(MB)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Complex')
        If (MB%Storage.eq.'StorFull') then
          VC%VecC = MatMul(VA%VecC,MB%MatI)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(VA%VecC,Tmp1%MatI)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MB%Storage.eq.'StorDiag') then
          VC%vecC = VA%vecC * MB%matI(:,1) 
        EndIf
      ElseIf (MQC_Vector_HaveComplex(VA) .and. MQC_Matrix_HaveReal(MB)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Complex')
        If (MB%Storage.eq.'StorFull') then
          VC%VecC = MatMul(VA%VecC,MB%MatR)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(VA%VecC,Tmp1%MatR)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MB%Storage.eq.'StorDiag') then
          VC%vecC = VA%VecC * MB%matR(:,1)
        EndIf
      ElseIf (MQC_Vector_HaveComplex(VA) .and. MQC_Matrix_HaveComplex(MB)) then
        Call MQC_Allocate_Vector(MB%NCol,VC,'Complex')
        If (MB%Storage.eq.'StorFull') then
          VC%VecC = MatMul(VA%VecC,MB%MatC)
        ElseIf (MB%Storage.eq.'StorSymm') then
          Tmp1 = MB
          Call MQC_Matrix_Symm2Full(Tmp1)
          VC%VecC = MatMul(VA%VecC,Tmp1%MatC)
          Call MQC_Deallocate_Matrix(Tmp1)
        elseif (MB%Storage.eq.'StorDiag') then
          VC%vecC = VA%VecC * MB%matC(:,1)
        EndIf
      EndIf

      End Function MQC_VectorMatrixDotProduct
!
!
!     PROCEDURE MQC_MatrixScalarProduct
      Function MQC_MatrixScalarProduct(Matrix,Scalar) Result(Matrix_Res)
!
!     This function calculates the product of an mqc_scalar and an mqc_vector
!
!     X. Sheng 2017.02
!
      Implicit None
      Type(MQC_Scalar),Intent(In)::Scalar
      Type(MQC_Matrix),Intent(In)::Matrix
      Type(MQC_Matrix)::Matrix_Res

      Call MQC_Deallocate_Matrix(Matrix_Res)
      If (Matrix%Data_Type.eq.'Integer' .and. Scalar%Data_Type.eq.'Integer') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Integer',Matrix%Storage)
        Matrix_Res%MatI = Scalar%ScaI * Matrix%MatI
      ElseIf (Matrix%Data_Type.eq.'Integer' .and. Scalar%Data_Type.eq.'Real') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Real',Matrix%Storage)
        Matrix_Res%MatR = Scalar%ScaR * Matrix%MatI
      Else If (Matrix%Data_Type.eq.'Real' .and. Scalar%Data_Type.eq.'Integer') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Real',Matrix%Storage)
        Matrix_Res%MatR = Scalar%ScaI * Matrix%MatR
      Else If (Matrix%Data_Type.eq.'Real' .and. Scalar%Data_Type.eq.'Real') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Real',Matrix%Storage)
        Matrix_Res%MatR = Scalar%ScaR * Matrix%MatR
      ElseIf (Matrix%Data_Type.eq.'Integer' .and. Scalar%Data_Type.eq.'Complex') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaC * Matrix%MatI
      Else If (Matrix%Data_Type.eq.'Real' .and. Scalar%Data_Type.eq.'Complex') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaC * Matrix%MatR
      ElseIf (Matrix%Data_Type.eq.'Complex' .and. Scalar%Data_Type.eq.'Integer') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaI * Matrix%MatC
      Else If (Matrix%Data_Type.eq.'Complex' .and. Scalar%Data_Type.eq.'Real') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaR * Matrix%MatC
      Else If (Matrix%Data_Type.eq.'Complex' .and. Scalar%Data_Type.eq.'Complex') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaC * Matrix%MatC
      Else
        Call MQC_Error_A('unrecognised data type in mqc_matrixscalarproduct', 6, &
             'Matrix%Data_Type', Matrix%Data_Type, &
             'Scalar%Data_Type', Scalar%Data_Type )
     EndIf
     

      End Function MQC_MatrixScalarProduct
!
!
!     PROCEDURE MQC_ScalarMatrixProduct
      Function MQC_ScalarMatrixProduct(Scalar,Matrix) Result(Matrix_Res)
!
!     This function calculates the product of an mqc_scalar and an mqc_matrix
!
!     X. Sheng 2017.02
!
      Implicit None
      Type(MQC_Scalar),Intent(In)::Scalar
      Type(MQC_Matrix),Intent(In)::Matrix
      Type(MQC_Matrix)::Matrix_Res

      Call MQC_Deallocate_Matrix(Matrix_Res)
      If (Matrix%Data_Type.eq.'Integer' .and. Scalar%Data_Type.eq.'Integer') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Integer',Matrix%Storage)
        Matrix_Res%MatI = Scalar%ScaI * Matrix%MatI
      ElseIf (Matrix%Data_Type.eq.'Integer' .and. Scalar%Data_Type.eq.'Real') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Real',Matrix%Storage)
        Matrix_Res%MatR = Scalar%ScaR * Matrix%MatI
      ElseIf (Matrix%Data_Type.eq.'Real' .and. Scalar%Data_Type.eq.'Integer') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Real',Matrix%Storage)
        Matrix_Res%MatR = Scalar%ScaI * Matrix%MatR
      Else If (Matrix%Data_Type.eq.'Real' .and. Scalar%Data_Type.eq.'Real') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Real',Matrix%Storage)
        Matrix_Res%MatR = Scalar%ScaR * Matrix%MatR
      ElseIf (Matrix%Data_Type.eq.'Integer' .and. Scalar%Data_Type.eq.'Complex') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaC * Matrix%MatI
      Else If (Matrix%Data_Type.eq.'Real' .and. Scalar%Data_Type.eq.'Complex') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaC * Matrix%MatR
      ElseIf (Matrix%Data_Type.eq.'Complex' .and. Scalar%Data_Type.eq.'Integer') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaI * Matrix%MatC
      Else If (Matrix%Data_Type.eq.'Complex' .and. Scalar%Data_Type.eq.'Real') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaR * Matrix%MatC
      Else If (Matrix%Data_Type.eq.'Complex' .and. Scalar%Data_Type.eq.'Complex') then
        Call MQC_Allocate_Matrix(Matrix%NRow,Matrix%NCol,Matrix_Res,'Complex',Matrix%Storage)
        Matrix_Res%MatC = Scalar%ScaC * Matrix%MatC
      Else
        Call MQC_Error_A('unrecognised data type in mqc_matrixscalarproduct', 6, &
             'Matrix%Data_Type', Matrix%Data_Type, &
             'Scalar%Data_Type', Scalar%Data_Type )
      EndIf

      End Function MQC_ScalarMatrixProduct
!
!
!     PROCEDURE MQC_Matrix_Matrix_Contraction
      Function MQC_Matrix_Matrix_Contraction(Matrix1,Matrix2) Result(Contraction)
!
!     This function returns the contration of two matrices, given as input
!     dummy arguments Matrix1 and Matrix2.
!
!     L. M. Thompson, 2016.
!
!
!     Variable Declarations.
!
      Implicit None
      Type(MQC_Matrix),Intent(In)::Matrix1,Matrix2
      Type(MQC_Scalar)::Contraction
      Integer(kind=int64)::N,M,i,j
!
      Contraction = 0.0
      N = Min(MQC_Matrix_Rows(Matrix1),MQC_Matrix_Rows(Matrix2))
      M = Min(MQC_Matrix_Columns(Matrix1),MQC_Matrix_Columns(Matrix2))
      Do i = 1,N
        Do j = 1,M
          Contraction = Contraction + Matrix1%at(i,j)*Matrix2%at(i,j)
        EndDo
      EndDo
!
      End Function MQC_Matrix_Matrix_Contraction
!
!
!     PROCEDURE MQC_Matrix_Scalar_Put
      Subroutine MQC_Matrix_Scalar_Put(Matrix,Scalar,I,J)
!
!     This subroutine updates the value of MQC_Matrix Matrix with the value
!     in MQC_Scalar Scalar at element I,J.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      Implicit None
      Integer(kind=int64),Intent(in)::I,J
      Class(MQC_Matrix),Intent(InOut)::Matrix
      Type(MQC_Scalar),Intent(In)::Scalar
      Integer(kind=int64)::IndI,IndJ

      IndI = I
      IndJ = J
      If (IndI.lt.0) IndI = Matrix%NRow + IndI + 1
      If (IndJ.lt.0) IndJ = Matrix%NRow + IndJ + 1
      If (IndI.eq.0.or.IndI.gt.Matrix%NRow) Call MQC_Error_I('Index I badly specified in mqc_matrix_scalar_put', 6, &
           'IndI', IndI, &
           'Matrix%NRow', Matrix%NRow )
      If (IndJ.eq.0.or.IndJ.gt.Matrix%NRow) Call MQC_Error_I('Index J badly specified in mqc_matrix_scalar_put', 6, &
           'IndJ', IndJ, &
           'Matrix%NRow', Matrix%NRow )
!     storage types only for full at the moment
      If (Matrix%Storage.ne.'StorFull') Call MQC_Error_A('MQC_Matrix_Scalar_Put only has StorFull implemented.', 6, &
           'Matrix%Storage', Matrix%Storage )
      If (Matrix%Data_Type.eq.'Integer') then
        If (Scalar%Data_Type.eq.'Integer') then
          Matrix%MatI(IndI,IndJ) = Scalar%ScaI
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Call MQC_Matrix_Copy_Int2Real(Matrix)
          Matrix%MatR(IndI,IndJ) = Scalar%ScaR
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Call MQC_Matrix_Copy_Int2Complex(Matrix)
          Matrix%MatC(IndI,IndJ) = Scalar%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Matrix_Scalar_Put', 6, &
              'Scalar%Data_Type', Scalar%Data_Type ) 
        EndIf
      ElseIf (Matrix%Data_Type.eq.'Real') then
        If (Scalar%Data_Type.eq.'Integer') then
          Matrix%MatR(IndI,IndJ) = Scalar%ScaI
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Matrix%MatR(IndI,IndJ) = Scalar%ScaR
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Call MQC_Matrix_Copy_Real2Complex(Matrix)
          Matrix%MatC(IndI,IndJ) = Scalar%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Matrix_Scalar_Put', 6, &
              'Scalar%Data_Type', Scalar%Data_Type ) 
        EndIf
      ElseIf (Matrix%Data_Type.eq.'Complex') then
        If (Scalar%Data_Type.eq.'Integer') then
          Matrix%MatC(IndI,IndJ) = Scalar%ScaI
        ElseIf (Scalar%Data_Type.eq.'Real') then
          Matrix%MatC(IndI,IndJ) = Scalar%ScaR
        ElseIf (Scalar%Data_Type.eq.'Complex') then
          Matrix%MatC(IndI,IndJ) = Scalar%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_Matrix_Scalar_Put', 6, &
               'Scalar%Data_Type', Scalar%Data_Type ) 
        EndIf
      Else
        Call MQC_Error_A('Matrix type not defined in MQC_Matrix_Scalar_Put', 6, &
             'Matrix%Data_Type', Matrix%Data_Type )
      EndIf

      End Subroutine MQC_Matrix_Scalar_Put
!
!
!     PROCEDURE MQC_Matrix_Initialize
      Subroutine MQC_Matrix_Initialize(Matrix,Rows,Columns,Scalar,Storage)
!
!     This subroutine initializes an MQC_Matrix Matrix of dimension Rows x
!     Columns with the value Scalar in all elements. <Storage> is a character
!     string that matches the meanings of argument <Storage> in routine
!     MQC_Allocate_Matrix.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      Implicit None
      Real(kind=real64)::Zero=0.0d0
      Integer(kind=int64),Intent(in)::Rows,Columns
      Class(*),Optional::Scalar
      Character(Len=*),Intent(in),OPTIONAL::Storage
      Class(MQC_Matrix),Intent(InOut)::Matrix
      Character(Len=64)::myStorage
!
      if(Present(Storage)) then
        myStorage = TRIM(Storage)
      else
        myStorage = 'StorFull'
      endIf
!
      If(Present(Scalar)) then
        Select Type (Scalar)
        Type is (Integer)
            Call MQC_Allocate_Matrix(Rows,Columns,Matrix,'Integer',TRIM(myStorage))
            Matrix%MatI = Scalar
        Type is (Real)
            Call MQC_Allocate_Matrix(Rows,Columns,Matrix,'Real',TRIM(myStorage))
            Matrix%MatR = Scalar
        Type is (Complex(Kind=real64))
            Call MQC_Allocate_Matrix(Rows,Columns,Matrix,'Complex',TRIM(myStorage))
            Matrix%MatC = Scalar
        Class Default
          Call MQC_Error_I('Scalar Type not defined in MQC_Matrix_Initialize', 6)
        End Select
      Else
        Call MQC_Allocate_Matrix(Rows,Columns,Matrix,'Real',TRIM(myStorage))
        Matrix%MatR = Zero
      EndIf
!
      Return
      End Subroutine MQC_Matrix_Initialize
!
!
!     PROCEDURE MQC_Matrix_Identity
      Subroutine MQC_Matrix_Identity(matrix,n,m)
!
!     This subroutine returns an identity matrix.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Matrix),intent(inOut)::matrix
      integer(kind=int64)::n,m,i
      type(MQC_Scalar)::one

      one = 1.0
      call matrix%initialize(n,m)
      Do i = 1,min(n,m)
        call matrix%put(one,i,i)
      EndDo

      End Subroutine MQC_Matrix_Identity
!
!
!     PROCEDURE MQC_Matrix_Set
      subroutine MQC_Matrix_Set(matrix,scalar,storage)
!
!     This subroutine sets a matrix that already exists to a value (default 0).
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Matrix),intent(inOut)::matrix
      class(*),optional::scalar
      character(len=*),intent(in),optional::storage
      integer(kind=int64)::rows,columns
      character(Len=64)::myStorage
      real(kind=real64)::zero=0.0d0
!
      if(Present(Storage)) then
        myStorage = TRIM(Storage)
      else
        myStorage = 'StorFull'
      endIf

      rows = mqc_matrix_rows(matrix) 
      columns = mqc_matrix_columns(matrix) 
      if(present(scalar)) then
        select type (scalar)
        type is (integer)
            call mqc_allocate_matrix(rows,columns,matrix,'Integer',TRIM(myStorage))
            Matrix%MatI = scalar
        type is (real)
            Call mqc_allocate_matrix(rows,columns,matrix,'Real',TRIM(myStorage))
            Matrix%MatR = scalar
        type is (complex(kind=real64))
            call mqc_allocate_matrix(rows,columns,matrix,'Complex',TRIM(myStorage))
            Matrix%MatC = scalar
        class default
          call mqc_error_i('scalar type not defined in MQC_Matrix_Set', 6)
        end select
      else
        if(Matrix%Data_Type.eq.'Integer') then 
          matrix%MatI = Zero
        elseIf(Matrix%Data_Type.eq.'Real') then 
          matrix%MatR = Zero
        elseIf(Matrix%Data_Type.eq.'Complex') then 
          matrix%MatC = Zero
        endIf
      endIf

      End Subroutine MQC_Matrix_Set
!
!
!     PROCEDURE MQC_Matrix_Norm
      Function MQC_Matrix_Norm(matrix,methodIn) result(norm)
!
!     This function returns the norm of an mqc matrix.
!     Different norms can be computed with optional argument methodIn.
!     'M' - max(abs(A(i,j)))
!     '1' - one norm
!     'I' - infinty norm
!     'F' - Frobenius norm (default)
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Matrix),intent(inOut)::matrix
      character(len=1),optional,intent(in)::methodIn
      character(len=1)::method
      character(len=64)::StorFlag,TypeFlag
      real(kind=real64),dimension(:),allocatable::work
      real(kind=real64)::dlange,zlange
      real(kind=real64),dimension(:,:),allocatable::temp
      Complex(Kind=real64),dimension(:,:),allocatable::tempC
      Type(MQC_Scalar)::norm

      if(Present(methodIn)) method = methodIn
      if(.not.Present(methodIn)) method = 'F'
      if(method.eq.'I') allocate(work(max(1,MQC_Matrix_Rows(matrix))))
      StorFlag = Matrix%Storage
      TypeFlag = Matrix%Data_Type
      if(Matrix%Storage.eq.'StorDiag') Call MQC_Matrix_Diag2Full(matrix)
      if(Matrix%Storage.eq.'StorSymm') Call MQC_Matrix_Diag2Symm(matrix)
      If(Matrix%Data_Type.eq.'Integer'.or.Matrix%Data_Type.eq.'Real') then 
        If(Matrix%Data_Type.eq.'Integer') Call MQC_Matrix_Copy_Int2Real(matrix)
        allocate(temp(MQC_Matrix_Rows(matrix),MQC_Matrix_Columns(matrix)))
        temp = matrix%matr
        norm = dlange(method,MQC_Matrix_Rows(matrix),MQC_Matrix_Columns(matrix), &
          temp,max(MQC_Matrix_Rows(matrix),1),work)
        if(allocated(work)) deallocate(work)
        if(allocated(temp)) deallocate(temp)
        if(TypeFlag.eq.'Integer') Call MQC_Matrix_Copy_Real2Int(matrix)
      ElseIf(matrix%Data_Type.eq.'Complex') then
        allocate(tempC(MQC_Matrix_Rows(matrix),MQC_Matrix_Columns(matrix)))
        tempC = matrix%matc
        norm = zlange(method,MQC_Matrix_Rows(matrix),MQC_Matrix_Columns(matrix), &
          tempC,max(MQC_Matrix_Rows(matrix),1),work)
        if(allocated(work)) deallocate(work)
        if(allocated(temp)) deallocate(temp)
      Else
        call mqc_error_a('unrecognised data type in mqc_matrix_norm', 6, &
             'matrix%Data_Type', matrix%Data_Type )
      EndIf
      if(StorFlag.eq.'StorDiag') Call MQC_Matrix_Full2Diag(matrix)
      if(StorFlag.eq.'StorSymm') Call MQC_Matrix_Full2Symm(matrix)

      End Function MQC_Matrix_Norm
!
!
!     PROCEDURE mqc_matrix_determinant
      function mqc_matrix_determinant(a) result(det)
!
!     This function returns the determinant of MQC matrix A.
!
!     Lee M. Thompson, 2016.
!
      implicit none
      class(mqc_matrix)::a
      type(mqc_scalar)::det
      character(len=64)::StorFlag,TypeFlag
      real(kind=real64),dimension(:,:),allocatable::temp
      Complex(Kind=real64),dimension(:,:),allocatable::tempC
      integer(kind=int64)::sgn,rows,cols,minrc,iError,i
      integer(kind=int32),dimension(:),allocatable::iPiv
!
!     Determinant of matrix A.
!
      rows = mqc_matrix_rows(a)
      cols = mqc_matrix_columns(a)
      minrc = min(rows,cols)
      storFlag = a%Storage
      typeFlag = a%Data_Type
      if(a%Storage.eq.'StorDiag') call mqc_matrix_diag2Full(a)
      if(a%Storage.eq.'StorSymm') call mqc_matrix_symm2Full(a)
      if(a%Data_Type.eq.'Integer'.or.a%Data_Type.eq.'Real') then
        If(a%Data_Type.eq.'Integer') call mqc_matrix_copy_int2Real(a)
        allocate(temp(rows,cols))
        temp = a%matr

        allocate(iPiv(minrc))
        iPiv = 0
        iError = 1
        call dgetrf(rows,cols,temp,rows,iPiv,iError)
        if(iError.lt.0) call mqc_error_i('DGETRF exited with error in mqc_matrix_determinant', 6, &
             'iError', iError )
        det = 1.0
        do i = 1,rows
          det = det%rval()*temp(i,i)
        enddo
        sgn = 1
        do i = 1,rows
          if(iPiv(i).ne.i) then
            sgn = -sgn
          endIf
        endDo
        det = sgn*det%rval()
        if(allocated(iPiv)) deallocate(iPiv)
        if(allocated(temp)) deallocate(temp)
        if(TypeFlag.eq.'Integer') Call MQC_Matrix_Copy_Real2Int(a)
      elseIf(a%Data_Type.eq.'Complex') then
        allocate(tempC(rows,cols))
        tempC = a%matc

        allocate(iPiv(minrc))
        iPiv = 0
        iError = 1
        call zgetrf(rows,cols,tempC,rows,iPiv,iError)
        if(iError.lt.0) call mqc_error_i('ZGETRF exited with error in mqc_matrix_determinant', 6, &
             'iError', iError )
        det = cmplx(1.0,0.0,kind=real64)
        do i = 1,rows
          det = det%cval()*tempC(i,i)
        enddo
        sgn = 1
        do i = 1,rows
          if(iPiv(i).ne.i) then
            sgn = -sgn
          endIf
        endDo
        det = sgn*det%cval()
        if(allocated(iPiv)) deallocate(iPiv)
        if(allocated(tempC)) deallocate(tempC)
      endIf
      if(StorFlag.eq.'StorDiag') Call MQC_Matrix_Full2Diag(a)
      if(StorFlag.eq.'StorSymm') Call MQC_Matrix_Full2Symm(a)

      end function mqc_matrix_determinant
!
!
!     PROCEDURE mqc_matrix_inverse
      function mqc_matrix_inverse(a) result(a_inv)
!
      implicit none
      class(mqc_matrix)::a
      type(mqc_matrix)::a_inv
      character(len=64)::StorFlag,TypeFlag
      real(kind=real64),dimension(:,:),allocatable::temp
      Complex(Kind=real64),dimension(:,:),allocatable::tempC
      Complex(Kind=real64),dimension(:),allocatable::workc
      integer(kind=int64)::iError,nDim1,nDim2
      integer(kind=int64),dimension(:),allocatable::iPiv,work
!
!     Inverse of matrix A.
!
      nDim1 = mqc_matrix_rows(a)
      nDim2 = mqc_matrix_columns(a)
      storFlag = a%Storage
      typeFlag = a%Data_Type
      if(a%Storage.eq.'StorDiag') call mqc_matrix_diag2Full(a)
      if(a%Storage.eq.'StorSymm') call mqc_matrix_diag2Symm(a)
      if(a%Data_Type.eq.'Integer'.or.a%Data_Type.eq.'Real') then
        If(a%Data_Type.eq.'Integer') call mqc_matrix_copy_int2Real(a)
        allocate(temp(nDim1,nDim2))
        temp = a%matr

        allocate(iPiv(min(nDim1,nDim2)),work(min(nDim1,nDim2)))
        iPiv = 0
        iError = 0
        call dgetrf(nDim1,nDim2,temp,nDim1,iPiv,iError)
        if (iError.gt.0) then 
          call mqc_error_i('Upper triangular matrix is singular', 6, 'iError', iError)
        else if (iError.lt.0) then
          call mqc_error_i('Illegal value in dgetrf matrix', 6, 'iError', iError)
        else
          call dgetri(NDim1,temp,NDim1,IPIV,Work,NDim1,IError)
          if(iError.gt.0) then
            call mqc_error_i('Matrix is not invertable', 6, 'iError', iError)
          else if (iError.lt.0) then
            call mqc_error_i('Illegal value in dgetri matrix', 6, 'iError', iError)
          end if
        end if
        a_inv = temp
        if(allocated(iPiv)) deallocate(iPiv)
        if(allocated(work)) deallocate(work)
        if(allocated(temp)) deallocate(temp)
        if(TypeFlag.eq.'Integer') Call MQC_Matrix_Copy_Real2Int(a)
      elseIf(a%Data_Type.eq.'Complex') then
        allocate(tempc(nDim1,nDim2))
        tempc = a%matc
        
        allocate(iPiv(min(nDim1,nDim2)),workc(min(nDim1,nDim2)))
        iPiv = 0
        iError = 0
        call zgetrf(nDim1,nDim2,tempc,nDim1,iPiv,iError)
        if (iError.gt.0) then
          call mqc_error_i('Upper triangular matrix is singular', 6, 'iError', iError)
        else if (iError.lt.0) then
          call mqc_error_i('Illegal value in zgetrf matrix', 6, 'iError', iError)
        else
          call zgetri(NDim1,tempC,NDim1,IPIV,WorkC,NDim1,IError)
          if(iError.gt.0) then
            call mqc_error_i('Matrix is not invertable', 6, 'iError', iError)
          else if(iError.lt.0) then
            call mqc_error_i('Illegal value in zgetri matrix', 6, 'iError', iError)
          end if
        end if
        a_inv = tempC
        if(allocated(iPiv)) deallocate(iPiv)
        if(allocated(workc)) deallocate(workc)
        if(allocated(tempC)) deallocate(tempC)
      endIf
      if(StorFlag.eq.'StorDiag') Call MQC_Matrix_Full2Diag(a)
      if(StorFlag.eq.'StorSymm') Call MQC_Matrix_Full2Symm(a)
!
      end function mqc_matrix_inverse
!
!
!     PROCEDURE mqc_matrix_trace
      function mqc_matrix_trace(matrix) result(trace)
!
!     This function computes the trace of a matrix
!
!     Variable Declarations...
!
      Implicit None
      class(mqc_matrix),intent(in)::matrix
      type(mqc_scalar)::trace
      integer(kind=int64)::i
!
      trace = 0.0
      do i = 1, min(mqc_matrix_rows(matrix),mqc_matrix_columns(matrix))
        trace = trace + matrix%at(i,i)
      enddo
!
      end function mqc_matrix_trace
!
!
!     PROCEDURE mqc_matrix_generalized_eigensystem
      subroutine mqc_matrix_generalized_eigensystem(a,b,eigenvals,reigenvecs,leigenvecs)
!
      implicit none
      class(mqc_matrix),intent(inOut)::a
      type(mqc_matrix),intent(inOut)::b
      type(mqc_vector),intent(out)::eigenvals
      type(mqc_matrix),intent(out),optional::reigenvecs,leigenvecs
      real(kind=real64),dimension(:,:),allocatable::tempA,tempB,rightvecs,leftvecs
      real(kind=real64),dimension(:),allocatable::eValsOut,work,eValsReal,eValsImag
      Complex(Kind=real64),dimension(:,:),allocatable::tempAC,tempBC,rightCvecs,leftCvecs
      Complex(Kind=real64),dimension(:),allocatable::eValsA,eValsB,workC
      integer(kind=int64)::iError=0,lWork,nDim1A,nDim2A,nDim1B,nDim2B,ilaenv,i,npdFlag
      character(len=64)::StorFlagA,StorFlagB,TypeFlagA,TypeFlagB
      type(mqc_vector)::tmp_vec
!
!     Generalized eigensystem solver.
!
      nDim1A = mqc_matrix_rows(a)
      nDim2A = mqc_matrix_columns(a)
      nDim1B = mqc_matrix_rows(b)
      nDim2B = mqc_matrix_columns(b)
      if(nDim1A.ne.nDim1B) call mqc_error_I('Order of matrices A and B is different in mqc_matrix_generalized_eigensystem', 6, &
           'nDim1A', nDim1A, &
           'nDim1B', nDim1B )
      if(nDim1A.ne.nDim2A) call mqc_error_i('Matrix A is not square in mqc_matrix_generalized_eigensystem', 6, &
           'nDim1A', nDim1A, &
           'nDim2A', nDim2A )
      if(nDim1B.ne.nDim2B) call mqc_error_i('Matrix B is not square in mqc_matrix_generalized_eigensystem', 6, &
           'nDim1B', nDim1B, &
           'nDim2B', nDim2B )
      storFlagA = a%Storage
      storFlagB = b%Storage
      typeFlagA = a%Data_Type
      typeFlagB = b%Data_Type
      if(a%Storage.eq.'StorDiag') call mqc_matrix_diag2Full(a)
      if(b%Storage.eq.'StorDiag') call mqc_matrix_diag2Full(b)
      if(a%Storage.eq.'StorSymm') call mqc_matrix_symm2Full(a)
      if(b%Storage.eq.'StorSymm') call mqc_matrix_symm2Full(b)
      If(a%Data_Type.eq.'Integer') call mqc_matrix_copy_int2Real(a)
      If(b%Data_Type.eq.'Integer') call mqc_matrix_copy_int2Real(b)
      If(a%Data_Type.eq.'Complex'.and.b%Data_Type.ne.'Complex') call mqc_matrix_copy_real2Complex(b)
      If(a%Data_Type.ne.'Complex'.and.b%Data_Type.eq.'Complex') call mqc_matrix_copy_real2Complex(a)
      npdFlag = 0

      if(MQC_Matrix_Test_Symmetric(a).and.MQC_Matrix_Test_Symmetric(b).and.a%Data_Type.eq.'Real' &
        .and.b%Data_Type.eq.'Real') then
        allocate(tempA(nDim1A,nDim2A),tempB(nDim1B,nDim2B))
        tempA = a%matR
        tempB = b%matR
        lWork = -1
        allocate(eValsOut(nDim2A),work(1))
        if(present(reigenvecs).or.present(leigenvecs)) then
          call dsygv(1,'V','L',nDim2A,tempA,nDim1A,tempB,nDim1B,eValsOut,work,lWork,iError)
          if(iError.ne.0) call mqc_error_i('Failure in SSYGV memory evaluation in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
          lWork = work(1)
          deallocate(work)
          allocate(work(lWork))
          call dsygv(1,'V','L',nDim2A,tempA,nDim1A,tempB,nDim1B,eValsOut,work,lWork,iError)
          if(iError.lt.0) call mqc_error_i('Failure in SSYGV routine in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
          if(iError.gt.0) then
            npdFlag = 1
          else
            if(present(reigenvecs)) reigenvecs = tempA
            if(present(leigenvecs)) leigenvecs = tempA
          endIf
        else
          call dsygv(1,'N','L',nDim2A,tempA,nDim1A,tempB,nDim1B,eValsOut,work,lWork,iError)
          if(iError.ne.0) call mqc_error_i('Failure in SSYGV memory evaluation in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
          lWork = work(1)
          deallocate(work)
          allocate(work(lWork))
          call dsygv(1,'N','L',nDim2A,tempA,nDim1A,tempB,nDim1B,eValsOut,work,lWork,iError)
          if(iError.lt.0) call mqc_error_i('Failure in SSYGV routine in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
          if(iError.gt.0) npdFlag = 1
        endIf
        if(iError.eq.0) eigenvals = eValsOut
        if(allocated(eValsOut)) deallocate(eValsOut)
        if(allocated(work)) deallocate(work)
        if(allocated(tempA)) deallocate(tempA)
        if(allocated(tempB)) deallocate(tempB)
      endIf
      
      if(.not.MQC_Matrix_Test_Symmetric(a).or..not.MQC_Matrix_Test_Symmetric(b).or.a%Data_Type.ne.'Real' &
        .or.b%Data_Type.ne.'Real'.or.npdFlag.eq.1) then
        if(a%Data_Type.eq.'Real'.and.b%Data_Type.eq.'Real') then
          allocate(tempA(nDim1A,nDim2A),tempB(nDim1B,nDim2B))
          tempA = a%matR
          tempB = b%matR
          lWork = -1
          allocate(eValsOut(nDim2A),work(1))
          allocate(eValsReal(nDim2A),eValsImag(NDim2A))
          if(present(reigenvecs).or.present(leigenvecs)) then
            allocate(rightvecs(nDim2A,nDim2A),leftvecs(nDim2A,nDim2A))
            call dggev('V','V',nDim2A,tempA,nDim1A,tempB,NDim1B,eValsReal,eValsImag,eValsOut,leftvecs,NDim2A, &
              rightvecs,nDim2A,work,lWork,iError)
            if(iError.ne.0) call mqc_error_i('Failure in DGGEV memory evaluation in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
            lWork = work(1)
            deallocate(work)
            allocate(work(lWork))
            call dggev('V','V',nDim2A,tempA,nDim1A,tempB,NDim1B,eValsReal,eValsImag,eValsOut,leftvecs,nDim2A, &
              rightvecs,nDim2A,work,lWork,iError)
            if(iError.ne.0) call mqc_error_i('Failure in DGGEV routine in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
            if(present(reigenvecs)) reigenvecs = rightvecs
            if(present(leigenvecs)) leigenvecs = leftvecs
            do i = 1, NDim2A
              if(eValsImag(i).ne.0.0) then
                if(eValsImag(i-1).ne.0.0) cycle
                if(present(reigenvecs)) then
                  tmp_vec = cmplx(rightvecs(:,i),rightvecs(:,i+1),kind=real64)
                  call reigenvecs%vput(tmp_vec,[0],[i])
                  tmp_vec = cmplx(rightvecs(:,i),-rightvecs(:,i+1),kind=real64)
                  call reigenvecs%vput(tmp_vec,[0],[i+1])
                endIf
                if(present(leigenvecs)) then
                  tmp_vec = cmplx(leftvecs(:,i),leftvecs(:,i+1),kind=real64)
                  call leigenvecs%vput(tmp_vec,[0],[i])
                  tmp_vec = cmplx(leftvecs(:,i),-leftvecs(:,i+1),kind=real64)
                  call leigenvecs%vput(tmp_vec,[0],[i+1])
                endIf
              endIf
            endDo
          else
            allocate(leftvecs(1,1),rightvecs(1,1))
            call dggev('N','N',nDim2A,tempA,nDim1A,tempB,NDim1B,eValsReal,eValsImag,eValsOut,leftvecs,1, &
              rightvecs,1,work,lWork,iError)
            if(iError.ne.0) call mqc_error_i('Failure in DGGEV memory evaluation in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
            lWork = work(1)
            deallocate(work)
            allocate(work(lwork))
            call dggev('N','N',nDim2A,tempA,nDim1A,tempB,NDim1B,eValsReal,eValsImag,eValsOut,leftvecs,1, &
              rightvecs,1,work,lwork,iError)
            if(iError.ne.0) call mqc_error_i('Failure in DGGEV routine in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
          endIf
          eigenvals = cmplx(eValsReal/eValsOut,eValsImag/eValsOut,kind=real64)
          if(allocated(eValsReal)) deallocate(eValsReal)
          if(allocated(eValsImag)) deallocate(eValsImag)
          if(allocated(rightvecs)) deallocate(rightvecs)
          if(allocated(leftvecs)) deallocate(leftvecs)
          if(allocated(eValsOut)) deallocate(eValsOut)
          if(allocated(work)) deallocate(work)
          if(allocated(tempA)) deallocate(tempA)
          if(allocated(tempB)) deallocate(tempB)
        elseIf(a%Data_Type.eq.'Complex'.and.b%Data_Type.eq.'Complex') then
          allocate(tempAC(nDim1A,nDim2A),tempBC(nDim1B,nDim2B))
          tempAC = a%matC
          tempBC = b%matC
          lWork  = -1
          allocate(eValsA(nDim2A),eValsB(NDim2A),workC(1),work(8*NDim2A))
          if(present(reigenvecs).or.present(leigenvecs)) then
            allocate(rightCvecs(nDim2A,nDim2A),leftCvecs(nDim2A,nDim2A))
            call zggev('V','V',nDim2A,tempAC,nDim1A,tempBC,NDim1B,eValsA,eValsB,leftCvecs,nDim2A, &
              rightCvecs,nDim2A,workC,lWork,work,iError)
            if(iError.ne.0) call mqc_error_i('Failure in ZGGEV routine in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
            lWork = workC(1)
            deallocate(workC)
            allocate(workC(lWork))
            call zggev('V','V',nDim2A,tempAC,nDim1A,tempBC,NDim1B,eValsA,eValsB,leftCvecs,nDim2A, &
              rightCvecs,nDim2A,workC,lWork,work,iError)
            if(iError.ne.0) call mqc_error_i('Failure in ZGGEV routine in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
            if(present(reigenvecs)) reigenvecs = rightCvecs
            if(present(leigenvecs)) leigenvecs = leftCvecs
          else
            allocate(rightCvecs(1,1),leftCvecs(1,1))
            call zggev('N','N',nDim2A,tempAC,nDim1A,tempBC,NDim1B,eValsA,eValsB,leftCvecs,1, &
              rightCvecs,1,workC,lWork,work,iError)
            if(iError.ne.0) call mqc_error_i('Failure in ZGGEV routine in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
            lWork = workC(1)
            deallocate(workC)
            allocate(workC(lWork))
            call zggev('N','N',nDim2A,tempAC,nDim1A,tempBC,NDim1B,eValsA,eValsB,leftCvecs,1, &
              rightCvecs,1,workC,lWork,work,iError)
            if(iError.ne.0) call mqc_error_i('Failure in ZGGEV routine in mqc_matrix_generalized_eigensystem', 6, &
                 'iError', iError)
          endIf
          eigenvals = eValsA/eValsB
          if(allocated(eValsA)) deallocate(eValsA)
          if(allocated(eValsB)) deallocate(eValsB)
          if(allocated(rightCvecs)) deallocate(rightCvecs)
          if(allocated(leftCvecs)) deallocate(leftCvecs)
          if(allocated(workC)) deallocate(workC)
          if(allocated(work)) deallocate(work)
          if(allocated(tempAC)) deallocate(tempAC)
          if(allocated(tempBC)) deallocate(tempBC)
        else
          call mqc_error_A('Incompatible data types in mqc_matrix_generalized_eigensystem',6,'a%Data_Type',a%Data_Type, &
            'b%Data_Type',b%Data_Type)
        endIf
      endIf

      if(typeFlagA.eq.'Integer') Call MQC_Matrix_Copy_Real2Int(a)
      if(typeFlagB.eq.'Integer') Call MQC_Matrix_Copy_Real2Int(b)
      if(typeFlagA.eq.'Real'.and.a%Data_Type.eq.'Complex') Call MQC_Matrix_Copy_Complex2Real(a)
      if(typeFlagB.eq.'Real'.and.b%Data_Type.eq.'Complex') Call MQC_Matrix_Copy_Complex2Real(b)
      if(storFlagA.eq.'StorDiag') Call MQC_Matrix_full2Diag(a)
      if(storFlagB.eq.'StorDiag') Call MQC_Matrix_full2Diag(b)
      if(storFlagA.eq.'StorSymm') Call MQC_Matrix_full2Symm(a)
      if(storFlagB.eq.'StorSymm') Call MQC_Matrix_full2Symm(b)
!
      end subroutine mqc_matrix_generalized_eigensystem
!
!
!     PROCEDURE MQC_Matrix_SVD
      Subroutine MQC_Matrix_SVD(A,EVals,EUVecs,EVVecs)
!
!     Do SVD decomposition of matrix A.
!   
!     L. M. Thompson, 2017
!
      Implicit None
      Class(MQC_Matrix),Intent(InOut)::A
      Type(MQC_Matrix),Optional,Intent(InOut)::EUVecs
      Type(MQC_Matrix),Optional,Intent(InOut)::EVVecs
      Type(MQC_Vector),Optional,Intent(InOut)::EVals
!
      Integer(kind=int64)::NDim1=0,NDim2=0,Work_Length=-1,IError=1
      Real(kind=real64),Dimension(:),Allocatable::A_SVals,Work
      Complex(Kind=real64),Dimension(:),Allocatable::WorkC
      Real(kind=real64),Dimension(:,:),Allocatable::A_Temp,A_UVecs,A_VVecs
      Complex(Kind=real64),Dimension(:,:),Allocatable::AC_Temp,AC_UVecs,AC_VVecs
      Character(len=64)::typeFlag,storFlag
!
      NDim1 = MQC_Matrix_Rows(A)
      NDim2 = MQC_Matrix_Columns(A)
      typeFlag = A%Data_Type
      storFlag = A%Storage
 1050 Format( A, I10 )
!
      if(storFlag.eq.'StorDiag') return
      if(storFlag.eq.'StorSymm') call mqc_matrix_symm2Full(A)
      if(storFlag.eq.'StorFull') then
        if(typeFlag.eq.'Real'.or.typeFlag.eq.'Integer') then
          if(typeFlag.eq.'Integer') call mqc_matrix_copy_int2Real(A)
          Allocate(A_Temp(NDim1,NDim2),A_SVals(Min(NDim1,NDim2)),A_UVecs(NDim1,NDim1), &
            A_VVecs(NDim2,NDim2),Work(1))
          A_Temp = A
          Call DGESVD('A','A',NDim1,NDim2,A_Temp,NDim1,A_SVals,A_UVecs,NDim1,A_VVecs, &
            NDim2,Work,Work_Length,IError)
          Work_Length = Int(Work(1))
          DeAllocate(Work)
          Allocate(Work(Work_Length))
          Call DGESVD('A','A',NDim1,NDim2,A_Temp,NDim1,A_SVals,A_UVecs,NDim1,A_VVecs, &
            NDim2,Work,Work_Length,IError)
          If(IError.ne.0) Write(*,1050)' SVD FAILED: IError = ',IError
          DeAllocate(Work)
          Work_Length = -1
        elseIf(typeFlag.eq.'Complex') then
          Allocate(AC_Temp(NDim1,NDim2),A_SVals(Min(NDim1,NDim2)),AC_UVecs(NDim1,NDim1), &
            AC_VVecs(NDim2,NDim2),WorkC(1),Work(5*min(NDim1,NDim2)))
          AC_Temp = A
          Call ZGESVD('A','A',NDim1,NDim2,AC_Temp,NDim1,A_SVals,AC_UVecs,NDim1,AC_VVecs, &
            NDim2,WorkC,Work_Length,Work,IError)
          Work_Length = Int(WorkC(1))
          DeAllocate(WorkC)
          Allocate(WorkC(Work_Length))
          Call ZGESVD('A','A',NDim1,NDim2,AC_Temp,NDim1,A_SVals,AC_UVecs,NDim1,AC_VVecs, &
            NDim2,WorkC,Work_Length,Work,IError)
          If(IError.ne.0) Write(*,1050)' SVD FAILED: IError = ',IError
          DeAllocate(WorkC,Work)
          Work_Length = -1
        else
          call mqc_error_a('typeFlag not recognised in mqc_matrix_svd', 6, &
               'typeFlag', typeFlag )
        endIf
!
        if(present(EVals)) EVals = A_SVals
        if(present(EUVecs)) then
          if(allocated(A_UVecs)) EUVecs = A_UVecs
          if(allocated(AC_UVecs)) EUVecs = AC_UVecs
        endIf
        if(present(EVVecs)) then
          if(allocated(A_VVecs)) EVVecs = A_VVecs
          if(allocated(AC_VVecs)) EVVecs = AC_VVecs
        endIf
        if(allocated(A_SVals)) deallocate(A_SVals)
        if(allocated(A_UVecs)) deallocate(A_UVecs)
        if(allocated(AC_UVecs)) deallocate(AC_UVecs)
        if(allocated(A_VVecs)) deallocate(A_VVecs)
        if(allocated(AC_VVecs)) deallocate(AC_VVecs)
        if(allocated(A_Temp)) deallocate(A_Temp)
        if(allocated(AC_Temp)) deallocate(AC_Temp)
!
      else
        call mqc_error_a('storFlag not recognised in mqc_matrix_svd', 6, &
             'storFlag', storFlag )
     endIf
!
      if(typeFlag.eq.'Integer') call mqc_matrix_copy_real2Int(A)
      if(storFlag.eq.'StorSymm') call mqc_matrix_full2Symm(A)
!
      Return
      End Subroutine MQC_Matrix_SVD
!
!
!     PROCEDURE MQC_Matrix_RMS_Max
      Subroutine MQC_Matrix_RMS_Max(A,rms_A,max_A)
!
!     This subroutine returns the RMS and Max Absolute values from the
!     INPUT matrix A.
!
      Implicit None
      Class(MQC_Matrix),Intent(InOut)::A
      Real(kind=real64),Dimension(:,:),Allocatable::A_Temp
      Complex(Kind=real64),Dimension(:,:),Allocatable::AC_Temp
      Type(MQC_Scalar),Intent(Out)::rms_A,max_A
!
      Integer(kind=int64)::N
      Character(len=64)::typeFlag,storFlag
!
!     Do the work.
!
      typeFlag = A%Data_Type
      storFlag = A%Storage
      N = MQC_Matrix_Rows(A)*MQC_Matrix_Columns(A)
      if(storFlag.eq.'StorDiag') call mqc_matrix_diag2Full(A) 
      if(storFlag.eq.'StorSymm') call mqc_matrix_symm2Full(A)
      if(storFlag.eq.'StorFull') then
        if(typeFlag.eq.'Real'.or.typeFlag.eq.'Integer') then
          if(typeFlag.eq.'Integer') call mqc_matrix_copy_int2Real(A)
          allocate(A_Temp(MQC_Matrix_Rows(A),MQC_Matrix_Columns(A)))
          rms_A = sqrt(dot_product(Reshape(A_Temp,[N]),Reshape(A_Temp,[N]))/N)
          max_A = maxval(abs(A_Temp))
          deallocate(A_Temp)
        elseIf(typeFlag.eq.'Complex') then
          allocate(AC_Temp(MQC_Matrix_Rows(A),MQC_Matrix_Columns(A)))
          rms_A = sqrt(dot_product(Reshape(AC_Temp,[N]),Reshape(AC_Temp,[N]))/N)
          max_A = real(maxval(abs(AC_Temp)))
          deallocate(AC_Temp)
        else
          call mqc_error_a('Unrecognised type in mqc_matrix_rms_max', 6, &
               'typeFlag', typeFlag )
        endIf
      else
        call mqc_error_a('storFlag not recognised in mqc_matrix_rms_max', 6, &
             'storFlag', storFlag )
      endIf
!
      if(typeFlag.eq.'Integer') call mqc_matrix_copy_real2Int(A)
      if(storFlag.eq.'StorSymm') call mqc_matrix_full2Symm(A)
      if(storFlag.eq.'StorDiag') call mqc_matrix_full2Diag(A)
!
      Return
      End Subroutine MQC_Matrix_RMS_Max
!
!
!     PROCEDURE MQC_Matrix_Sqrt
      subroutine mqc_matrix_sqrt(A,eVals,eVecs)
!
!     This subroutine returns the square root of INPUT matrix A.
!
      Implicit None
      class(mqc_matrix),intent(inOut)::A
      type(mqc_matrix),optional,intent(inOut)::eVecs
      type(mqc_vector),optional,intent(inOut)::eVals
      type(mqc_matrix)::tmpMat
      type(mqc_vector)::tmpVec
      integer(kind=int64)::i
!
!     Do the work.
!
      call A%diag(tmpVec,tmpMat)
      if(present(eVals)) eVals = tmpVec 
      if(present(eVecs)) eVecs = tmpMat 
      call A%set()
      do i=1,MQC_Length_Vector(tmpVec)
        call A%put(sqrt(tmpVec%at(i)),i,i)
      endDo
      A = matmul(matmul(tmpMat,A),transpose(tmpMat))
!
      end subroutine mqc_matrix_sqrt
!
!----------------------------------------------------------------
!                                                               |
!     R4TENSOR PROCEDURES                                       |
!                                                               |
!----------------------------------------------------------------
!
!
!     PROCEDURE MQC_Allocate_R4Tensor
      Subroutine MQC_Allocate_R4Tensor(I,J,K,L,Tensor,Data_Type,Storage)
!
!     This subroutine allocates an MQC_R4Tensor rank-4 tensor of dimension
!     I x J x K x L
!
      Implicit None
      Integer(kind=int64),Intent(In)::I,J,K,L
      Character(Len=*),Intent(In)::Data_Type,Storage
      Type(MQC_R4Tensor),Intent(InOut)::Tensor
      Integer(kind=int64)::ICur,JCur,KCur,LCur

      If (I.lt.0 .or. J.lt.0 .or. K.lt.0 .or. L.lt.0) then
        Call MQC_Error_I('Dimensions less than zero in MQC_Allocate_R4Tensor', 6, &
             'I', I, &
             'J', J, &
             'K', K, &
             'L', L )
      EndIf

      Call MQC_Deallocate_R4Tensor(Tensor)

      If(Storage.eq.'StorFull' .or. Storage.eq.'Full' .or. Storage.eq.'full' .or. Storage.eq.'storfull') then
        Tensor%Storage = 'StorFull'
        ICur = I
        JCur = J
        KCur = K
        LCur = L
      Else
        Call MQC_Error_A('Storage type not recognised in MQC_Allocate_R4Tensor', 6, &
             'Storage', Storage)
      EndIf

      If(Data_Type.eq.'Integer' .or. Data_Type.eq.'integer' .or. Data_Type.eq.'Int' .or. Data_Type.eq.'int') then
        Tensor%Data_Type = 'Integer'
        If(Allocated(Tensor%ITen)) Deallocate(Tensor%ITen)
        Allocate(Tensor%ITen(ICur,JCur,KCur,LCur))
      ElseIf(Data_Type.eq.'Real' .or. Data_Type.eq.'real') then
        Tensor%Data_Type = 'Real'
        If(Allocated(Tensor%RTen)) Deallocate(Tensor%RTen)
        Allocate(Tensor%RTen(ICur,JCur,KCur,LCur))
      ElseIf(Data_Type.eq.'Complex' .or. Data_Type.eq.'complex') then
        Tensor%Data_Type = 'Complex'
        If(Allocated(Tensor%CTen)) Deallocate(Tensor%CTen)
        Allocate(Tensor%CTen(ICur,JCur,KCur,LCur))
      EndIf

      Tensor%I = I
      Tensor%J = J
      Tensor%K = K
      Tensor%L = L

      End Subroutine mqc_allocate_r4tensor
!
!
!     PROCEDURE MQC_Deallocate_R4Tensor
      Subroutine MQC_Deallocate_R4Tensor(Tensor)
!
!     This subroutine deallocates an MQC_R4Tensor rank-4 tensor.
!
      Implicit None
      Type(MQC_R4Tensor),Intent(InOut)::Tensor


      If(Allocated(Tensor%ITen)) Deallocate(Tensor%ITen)
      If(Allocated(Tensor%RTen)) Deallocate(Tensor%RTen)
      If(Allocated(Tensor%CTen)) Deallocate(Tensor%CTen)
      Tensor%Data_Type = ''
      Tensor%Storage = ''
      Tensor%I = 0
      Tensor%J = 0
      Tensor%K = 0
      Tensor%L = 0

      End Subroutine MQC_Deallocate_R4Tensor
!
!
!     PROCEDURE MQC_R4Tensor_At
      Function MQC_R4Tensor_At(Tensor,I,J,K,L) Result(Element)
!
!     This Function returns the MQC_Scalar Scalar at I, J, K, L.
!
      Implicit None
      Integer(kind=int64),Intent(In)::I,J,K,L
      Class(MQC_R4Tensor),Intent(In)::Tensor
      Type(MQC_Scalar)::Element
      Integer(kind=int64)::IndI,IndJ,IndK,IndL

      IndI = I
      IndJ = J
      IndK = K
      IndL = L
      If (IndI.lt.0) IndI = Tensor%I + IndI + 1
      If (IndJ.lt.0) IndJ = Tensor%J + IndJ + 1
      If (IndK.lt.0) IndK = Tensor%K + IndI + 1
      If (IndL.lt.0) IndL = Tensor%L + IndJ + 1
      If (IndI.eq.0.or.IndI.gt.Tensor%I) Call MQC_Error_I('Index I badly specified in mqc_r4tensor_at', 6, &
           'IndI', IndI, &
           'Tensor%I', Tensor%I )
      If (IndJ.eq.0.or.IndJ.gt.Tensor%J) Call MQC_Error_I('Index J badly specified in mqc_r4tensor_at', 6, &
           'IndJ', IndJ, &
           'Tensor%J', Tensor%J )
      If (IndK.eq.0.or.IndK.gt.Tensor%K) Call MQC_Error_I('Index K badly specified in mqc_r4tensor_at', 6, &
           'IndK', IndK, &
           'Tensor%K', Tensor%K )
      If (IndL.eq.0.or.IndL.gt.Tensor%L) Call MQC_Error_I('Index L badly specified in mqc_r4tensor_at', 6, &
           'IndL', IndL, &
           'Tensor%L', Tensor%L )
      If(Tensor%Data_Type.eq.'Integer') then
        Element = (Tensor%ITen(IndI,IndJ,IndK,IndL))
      ElseIf(Tensor%Data_Type.eq.'Real') then
        Element = (Tensor%RTen(IndI,IndJ,IndK,IndL))
      ElseIf(Tensor%Data_Type.eq.'Complex') then
        Element = (Tensor%CTen(IndI,IndJ,IndK,IndL))
      Else
        Call MQC_Error_A('Tensor type not defined in MQC_R4Tensor_At', 6, &
             'Tensor%Data_Type', Tensor%Data_Type )
      EndIf

      End Function MQC_R4Tensor_At
!
!
!     PROCEDURE MQC_R4Tensor_Put
      Subroutine MQC_R4Tensor_Put(Tensor,Element,I,J,K,L)
!
!     This Subroutine places the MQC_Scalar Scalar at I, J, K, L.
!
      Implicit None
      Integer(kind=int64),Intent(In)::I,J,K,L
      Class(MQC_R4Tensor),Intent(InOut)::Tensor
      Type(MQC_Scalar),Intent(In)::Element
      Integer(kind=int64)::IndI,IndJ,IndK,IndL

      IndI = I
      IndJ = J
      IndK = K
      IndL = L
      If (IndI.lt.0) IndI = Tensor%I + IndI + 1
      If (IndJ.lt.0) IndJ = Tensor%J + IndJ + 1
      If (IndK.lt.0) IndK = Tensor%K + IndI + 1
      If (IndL.lt.0) IndL = Tensor%L + IndJ + 1
      If (IndI.eq.0.or.IndI.gt.Tensor%I) Call MQC_Error_I('Index I badly specified in mqc_r4tensor_put', 6, &
           'IndI', IndI, &
           'Tensor%I', Tensor%I)
      If (IndJ.eq.0.or.IndJ.gt.Tensor%J) Call MQC_Error_I('Index J badly specified in mqc_r4tensor_put', 6, &
           'IndJ', IndJ, &
           'Tensor%J', Tensor%J)
      If (IndK.eq.0.or.IndK.gt.Tensor%K) Call MQC_Error_I('Index K badly specified in mqc_r4tensor_put', 6, &
           'IndK', IndK, &
           'Tensor%K', Tensor%K)
      If (IndL.eq.0.or.IndL.gt.Tensor%L) Call MQC_Error_I('Index L badly specified in mqc_r4tensor_put', 6, &
           'IndL', IndL, &
           'Tensor%L', Tensor%L)
      !print*, 'tensor%ten(i,j,k,l)= ', tensor%ten(i,j,k,l)
      If (Tensor%Data_Type.eq.'Integer') then
        If (Element%Data_Type.eq.'Integer') then
          Tensor%ITen(IndI,IndJ,IndK,IndL) = Element%ScaI
        ElseIf (Element%Data_Type.eq.'Real') then
          Tensor%ITen(IndI,IndJ,IndK,IndL) = Element%ScaR
        ElseIf (Element%Data_Type.eq.'Complex') then
          Tensor%ITen(IndI,IndJ,IndK,IndL) = Element%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_R4Tensor_Put', 6, &
               'Element%Data_Type', Element%Data_Type, &
               'Tensor%Data_Type', Tensor%Data_Type )
        EndIf
      ElseIf (Tensor%Data_Type.eq.'Real') then
        If (Element%Data_Type.eq.'Integer') then
          Tensor%RTen(IndI,IndJ,IndK,IndL) = Element%ScaI
        ElseIf (Element%Data_Type.eq.'Real') then
          Tensor%RTen(IndI,IndJ,IndK,IndL) = Element%ScaR
        ElseIf (Element%Data_Type.eq.'Complex') then
          Tensor%RTen(IndI,IndJ,IndK,IndL) = Element%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_R4Tensor_Put', 6, &
               'Element%Data_Type', Element%Data_Type, &
               'Tensor%Data_Type', Tensor%Data_Type )
        EndIf
      ElseIf (Tensor%Data_Type.eq.'Complex') then
        If (Element%Data_Type.eq.'Integer') then
          Tensor%CTen(IndI,IndJ,IndK,IndL) = Element%ScaI
        ElseIf (Element%Data_Type.eq.'Real') then
          Tensor%CTen(IndI,IndJ,IndK,IndL) = Element%ScaR
        ElseIf (Element%Data_Type.eq.'Complex') then
          Tensor%CTen(IndI,IndJ,IndK,IndL) = Element%ScaC
        Else
          Call MQC_Error_A('Scalar type not defined in MQC_R4Tensor_Put', 6, &
               'Element%Data_Type', Element%Data_Type, &
               'Tensor%Data_Type', Tensor%Data_Type )
        EndIf
      Else
        Call MQC_Error_A('Tensor type not defined in MQC_R4Tensor_Put', 6, &
             'Tensor%Data_Type', Tensor%Data_Type )
      EndIf
      !call mqc_output_scalar(tensor%ten(i,j,k,l),element)
      !tensor%ten(i,j,k,l) = 60
      !print*, 'tensor%ten(i,j,k,l)= ', tensor%ten(i,j,k,l)

      End Subroutine MQC_R4Tensor_Put
!
!
!     PROCEDURE MQC_Print_R4Tensor_Algebra1
      Subroutine MQC_Print_R4Tensor_Algebra1(Tensor,IOut,Header,blank_at_top,blank_at_bottom)
!
!     This Subroutine prints the MQC_R4Tensor Tensor.
!
      Implicit None
      Integer(kind=int64),Intent(In)::IOut
      Class(MQC_R4Tensor),Intent(In)::Tensor
      Integer(kind=int64)::I,J,K,L
      Character(Len=*),intent(in),optional::Header
      Logical,Optional::blank_at_top,blank_at_bottom

 1000 Format(1x,A)
 1020 Format( " " )
 1100 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',I10)
 1200 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',F15.8)
 1300 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',F12.5,F11.5,"i")

      If(Present(blank_at_top)) then
        If(blank_at_top) Write(IOut,1020)
      EndIf
      If(Present(Header)) Write(IOut,1000) Trim(Header)
      Do I = 1, Tensor%I
        Do J = 1, Tensor%J
          Do K = 1, Tensor%K
            Do l = 1, Tensor%L
              If(Tensor%Data_Type.eq.'Integer') then
                Write(IOut,1100) I,J,K,L,Tensor%ITen(I,J,K,L)
              ElseIf(Tensor%Data_Type.eq.'Real') then
                Write(IOut,1200) I,J,K,L,Tensor%RTen(I,J,K,L)
              ElseIf(Tensor%Data_Type.eq.'Complex') then
                Write(IOut,1300) I,J,K,L,Tensor%CTen(I,J,K,L)
              EndIf
            EndDo
          EndDo
        EndDo
      EndDo
      If(present(blank_at_bottom)) then
        If(blank_at_bottom) Write(IOut,1020)
      EndIf

      End Subroutine MQC_Print_R4Tensor_Algebra1
!
!
!     PROCEDURE MQC_Set_Array2Tensor
      Subroutine MQC_Set_Array2Tensor(TensorOut,ArrayIn)
!
!     This subroutine sets an MQC_R4Tensor rank-4 tensor to an array.
!
      Implicit None
      Class(*),Intent(In)::ArrayIn(:,:,:,:)
      Type(MQC_R4Tensor),Intent(InOut)::TensorOut

      Call MQC_Deallocate_R4Tensor(TensorOut)
      Select Type (ArrayIn)
      Type Is (Integer)
        Call MQC_Allocate_R4Tensor(Size(ArrayIn,1),Size(ArrayIn,2),Size(ArrayIn,3), &
          Size(ArrayIn,4),TensorOut,'Integer','StorFull')
        TensorOut%ITen = Arrayin
      Type Is (Real)
        Call MQC_Allocate_R4Tensor(Size(ArrayIn,1),Size(ArrayIn,2),Size(ArrayIn,3), &
          Size(ArrayIn,4),TensorOut,'Real','StorFull')
        TensorOut%RTen = ArrayIn
      Type Is (Complex(Kind=real64))
        Call MQC_Allocate_R4Tensor(Size(ArrayIn,1),Size(ArrayIn,2),Size(ArrayIn,3), &
          Size(ArrayIn,4),TensorOut,'Complex','StorFull')
        TensorOut%CTen = ArrayIn
      Class Default
        Call MQC_Error_I('Array type not determined in MQC_Array2Tensor', 6)
      End Select

      End Subroutine MQC_Set_Array2Tensor
!
!
!     PROCEDURE MQC_R4Tensor_Initialize
      Subroutine MQC_R4Tensor_Initialize(R4Tensor,I,J,K,L,Scalar)
!
!     This subroutine initializes an MQC_R4Tensor rank-4 tensor of dimension
!     I x J x K x L with the value Scalar in all elements.
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
!
      Implicit None
      Real(kind=real64)::Zero=0.0d0
      Integer(kind=int64),Intent(In)::I,J,K,L
      Class(*),Optional::Scalar
      Class(MQC_R4Tensor),Intent(InOut)::R4Tensor

      If(Present(Scalar)) then
        Select Type (Scalar)
        Type is (Integer)
          Call MQC_Allocate_R4Tensor(I,J,K,L,R4Tensor,'Integer','StorFull')
          R4Tensor%ITen = Scalar
        Type is (Real)
          Call MQC_Allocate_R4Tensor(I,J,K,L,R4Tensor,'Real','StorFull')
          R4Tensor%RTen = Scalar
        Type is (Complex(Kind=real64))
          Call MQC_Allocate_R4Tensor(I,J,K,L,R4Tensor,'Complex','StorFull')
          R4Tensor%CTen = Scalar
        Class Default
          Call MQC_Error_I('Scalar Type not defined in MQC_R4Tensor_Initialize', 6)
        End Select
      Else
        Call MQC_Allocate_R4Tensor(I,J,K,L,R4Tensor,'Real','StorFull')
        R4Tensor%RTen = Zero
      EndIf

      End Subroutine MQC_R4Tensor_Initialize
!
!     PROCEDURE MQC_Matrix_SymmSymmR4Tensor_Put_Real
      Subroutine MQC_Matrix_SymmSymmR4Tensor_Put_Real(r4Tensor,symmSymmMatrixIn)
!
!     This subroutine loads a (MQC) rank-4 tensor <mat> from an input
!     (fortran intrinsic type) real packed-symmetric-symmetric matrix <symmSymmMatrixIn>,
!     which should be passed as a rank-1 array with size n(n+1)/2 where n is the
!     leading dimension of the matrix.
!
!     L. M. Thompson, 2017.
!
!
      implicit none
      class(MQC_R4Tensor),intent(inOut)::r4Tensor
      real(kind=real64),dimension(:),intent(in)::symmSymmMatrixIn
!
      real(kind=real64),dimension(:,:,:,:),allocatable::eRIs
      integer(kind=int64)::ii,jj,kk,ll,n,nSymm
      logical::DEBUG=.false.
!
 5000 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',F15.8)
!
!
!     Begin by figuring out the leading dimension of the matrix. Then, ensure
!     mat is allocated correctly and then use mat%mput routine.
!
      nSymm = SIZE(symmSymmMatrixIn)
      n = 0.5*(sqrt(4.0*sqrt(8.0*nSymm+1)-3)-1)
      allocate(eRIs(n,n,n,n))
      do ii = 1,n
        do jj = 1,ii
          do kk = 1,n
            do ll = 1,kk
              if( (((ii*(ii-1))/2)+jj) .lt. (((kk*(kk-1))/2)+ll) ) cycle 
              if(DEBUG) write(*,5000) ii,jj,kk,ll,symmSymmMatrixIn(symIndexHash(ii,jj,kk,ll))
              eRIs(ii,jj,kk,ll) = symmSymmMatrixIn(symIndexHash(ii,jj,kk,ll)) 
              eRIs(jj,ii,kk,ll) = eRIs(ii,jj,kk,ll) 
              eRIs(ii,jj,ll,kk) = eRIs(ii,jj,kk,ll) 
              eRIs(jj,ii,ll,kk) = eRIs(ii,jj,kk,ll) 
              eRIs(kk,ll,ii,jj) = eRIs(ii,jj,kk,ll) 
              eRIs(kk,ll,jj,ii) = eRIs(ii,jj,kk,ll) 
              eRIs(ll,kk,ii,jj) = eRIs(ii,jj,kk,ll) 
              eRIs(ll,Kk,jj,ii) = eRIs(ii,jj,kk,ll) 
            endDo
          endDo
        endDo
      endDo
      r4Tensor = eRIs
      deallocate(eRIs)
!
      return
      end subroutine MQC_Matrix_SymmSymmR4Tensor_Put_Real
!
!
!     PROCEDURE MQC_Matrix_SymmSymmR4Tensor_Put_Complex
      Subroutine MQC_Matrix_SymmSymmR4Tensor_Put_Complex(r4Tensor,symmSymmMatrixIn)
!
!     This subroutine loads a (MQC) rank-4 tensor <mat> from an input
!     (fortran intrinsic type) complex packed-symmetric-symmetric matrix <symmSymmMatrixIn>,
!     which should be passed as a rank-1 array with size n(n+1)/2 where n is the
!     leading dimension of the matrix.
!
!     L. M. Thompson, 2017.
!
!
      implicit none
      class(MQC_R4Tensor),intent(inOut)::r4Tensor
      complex(kind=real64),dimension(:),intent(in)::symmSymmMatrixIn
!
      complex(kind=real64),dimension(:,:,:,:),allocatable::eRIs
      integer(kind=int64)::ii,jj,kk,ll,n,nSymm
      logical::DEBUG=.false.
!
 5000 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',F12.5,F11.5,"i")
!
!
!     Begin by figuring out the leading dimension of the matrix. Then, ensure
!     mat is allocated correctly and then use mat%mput routine.
!
      nSymm = SIZE(symmSymmMatrixIn)
      n = 0.5*(sqrt(4.0*sqrt(8.0*nSymm+1)-3)-1)
      allocate(eRIs(n,n,n,n))
      do ii = 1,n
        do jj = 1,ii
          do kk = 1,n
            do ll = 1,kk
              if( (((ii*(ii-1))/2)+jj) .lt. (((kk*(kk-1))/2)+ll) ) cycle 
              if(DEBUG) write(*,5000) ii,jj,kk,ll,symmSymmMatrixIn(symIndexHash(ii,jj,kk,ll))
              eRIs(ii,jj,kk,ll) = symmSymmMatrixIn(symIndexHash(ii,jj,kk,ll)) 
              eRIs(jj,ii,kk,ll) = eRIs(ii,jj,kk,ll) 
              eRIs(ii,jj,ll,kk) = eRIs(ii,jj,kk,ll) 
              eRIs(jj,ii,ll,kk) = eRIs(ii,jj,kk,ll) 
              eRIs(kk,ll,ii,jj) = eRIs(ii,jj,kk,ll) 
              eRIs(kk,ll,jj,ii) = eRIs(ii,jj,kk,ll) 
              eRIs(ll,kk,ii,jj) = eRIs(ii,jj,kk,ll) 
              eRIs(ll,Kk,jj,ii) = eRIs(ii,jj,kk,ll) 
            endDo
          endDo
        endDo
      endDo
      r4Tensor = eRIs
      deallocate(eRIs)
!
      return
      end subroutine MQC_Matrix_SymmSymmR4Tensor_Put_Complex
!
!
      End Module MQC_Algebra
