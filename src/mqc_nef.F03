   Module MQC_NEF

      use omp_lib
      use mqc_algebra
      use mqc_est
      Use iso_fortran_env, only: int32, int64, real64

!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!
!>    \brief <b> </b>
      Type MQC_NEF_Vector
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Main_Vector
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha,Size_vector_beta
        Integer,Private::NumSpace=0
      Contains
        Procedure, Public::print => mqc_print_nef_vector
        Procedure, Public::init => mqc_nef_vector_initialize
        Procedure, public::getBlockVec => mqc_nef_vector_output_block_vector
        Procedure, Public::impv => mqc_nef_vector_import_vector
        Procedure, Public::expand => mqc_nef_vector_dimention_expand
      End Type
!
!>    \brief <b> </b>
      Type MQC_NEF_Matrix
        Type(MQC_SCF_Integral),Dimension(:),Allocatable,Private::Main_Matrix
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix_Col,Map_Matrix_Row
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_col,Size_vector_beta_col
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_row,Size_vector_beta_row
        Integer,Private::NumSpace_col=0,NumSpace_row=0
      Contains
        Procedure, Public::print => mqc_print_nef_matrix
        Procedure, Public::init => mqc_nef_matrix_initialize
        Procedure, public::collapseIndex => mqc_nef_matrix_collapse_index
        Procedure, public::collapseMap => mqc_nef_matrix_collapse_map
        Procedure, public::blockIndex => mqc_nef_matrix_index
        Procedure, public::at => mqc_nef_matrix_at
        Procedure, Public::imp => mqc_nef_matrix_import_integral
        Procedure, Public::export => mqc_nef_matrix_export_integral
        Procedure, Public::slice => mqc_nef_matrix_dimention_slice
        Procedure, Public::collapse => mqc_nef_matrix_dimention_collapse
!>    \brief <b> Outer Product of two variables into a two NEF R4</b>
        Procedure, Public::outer2R4 => mqc_NEFMatrix_NEFMatrix_outer_product
!>    \brief <b> Outer Product of two variables into a NEF R4</b>
        Procedure, Public::crossouter2R4 => mqc_NEFMatrix_NEFMatrix_crossouter_product
        Procedure, public::getBlockMat => mqc_nef_matrix_output_block_matrix
      End Type
!
!>    \brief <b> </b>
      Type MQC_NEF_R4
        Type(MQC_TwoERIs),Dimension(:),Allocatable,Private::Main_R4
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix_1,Map_Matrix_2,&
                                                                    Map_Matrix_3,Map_Matrix_4
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_1,Size_vector_beta_1
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_2,Size_vector_beta_2
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_3,Size_vector_beta_3
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_4,Size_vector_beta_4
        Integer,Private::NumSpace_1=0,NumSpace_2=0,NumSpace_3=0,NumSpace_4=0
      Contains
        Procedure, Public::print => mqc_print_nef_r4
        Procedure, public::at => mqc_nef_r4_at
        Procedure, public::blockSize => mqc_nef_r4_dimension
        Procedure, public::blockSizeVec => mqc_nef_r4_dimension_vector
        Procedure, public::blockMap => mqc_nef_r4_map
        Procedure, public::blockIndex => mqc_nef_r4_index
        Procedure, Public::collapse => mqc_nef_r4_dimention_collapse
        Procedure, Public::collapseIndex => mqc_nef_r4_collapse_index
        Procedure, public::collapseMap => mqc_nef_r4_collapse_map
      End Type
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!
!>    \brief <b> Multiplies two objects</b>
      Interface MatMul
        Module Procedure MQC_NEFMatrix_Integral_Multiply
        Module Procedure MQC_Integral_NEFMatrix_Multiply
        Module Procedure MQC_NEFMatrix_NEFMatrix_Multiply
      End Interface
!
!>    \brief <b> Returns the Hermitian conjugate</b>
      Interface Dagger
        Module Procedure MQC_NEFMatrix_Conjugate_Transpose
        Module Procedure MQC_NEFR4_Conjugate_Transpose
      End Interface
!
!>    \brief <b> Returns the conjugate</b>
      Interface Conjg
        Module Procedure MQC_NEFMatrix_Conjugate
      End Interface
!
!>    \brief <b> Returns the transpose</b>
      Interface Transpose
        Module Procedure MQC_NEFMatrix_Transpose
      End Interface
!
!>    \brief <b> Contracts two objects</b>
      Interface Contraction
        Module Procedure mqc_NEFMatrix_NEFMatrix_contraction
        Module Procedure mqc_NEF_R4_NEF_R4_contraction
        Module Procedure MQC_NEF_R4_NEF_Matrix_Contraction
      End Interface
!
!>    \brief <b> Contracts two objects</b>
      Interface PartialContraction
        Module Procedure mqc_NEF_R4_NEF_R4_partial_contraction
      End Interface
!
!>    \brief <b> Contracts two objects</b>
      Interface PartialContraction2R4
        Module Procedure mqc_NEF_R4_NEF_R4_partial_contraction_2NEF_R4
        Module Procedure mqc_NEF_R4_ERI_partial_contraction_2NEF_R4
        Module Procedure mqc_ERI_NEF_R4_partial_contraction_2NEF_R4
      End Interface
!
!>    \brief <b> Returns the block length along the specified dimension</b>
      Interface Size
        Module Procedure mqc_nef_r4_dimension
      End Interface
!

!
!
!----------------------------------------------------------------
!                                                               |
!     OPERATOR INTERFACES                                       |
!                                                               |
!----------------------------------------------------------------
!
!
!>    \brief <b> Sums two variables</b>
      Interface Operator (+)
        Module Procedure MQC_NEF_Vector_Sum 
        Module Procedure MQC_NEFMatrix_Sum
        Module Procedure MQC_NEF_R4_Sum 
      End Interface
!      
!>    \brief <b> Multiplies two variables</b>
      Interface Operator (*)
        Module Procedure MQC_Scalar_NEF_Vector_Multiply
        Module Procedure MQC_NEF_Vector_Scalar_Multiply
        Module Procedure MQC_Scalar_NEFMatrix_Multiply
        Module Procedure MQC_NEFMatrix_Scalar_Multiply
        Module Procedure MQC_Scalar_NEFR4_Multiply
        Module Procedure MQC_NEFR4_Scalar_Multiply
      End Interface

      !>    \brief <b> Subtracts two variables</b>
      Interface Operator (-)
        Module Procedure MQC_NEF_Vector_Difference
        Module Procedure MQC_NEFMatrix_Difference
        Module Procedure MQC_NEF_R4_Difference
      End Interface

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
!
!     PROCEDURE MQC_Print_NEF_Vector
!
!>    \brief <b> MQC_Print_NEF_Vector is a subroutine used to print an MQC
!>    NEF Vector object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_NEF_Vector is a subroutine used to print an MQC NEF Vector object. 
!>    Blank_At_Top and Blank_At_Bottom are optional logical arguments to print blank
!>    lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVec
!>    \verbatim
!>        Integral is Class(MQC_NEF_Vector)
!>        The MQC NEF Vector variable to be printed.
!>    \endverbatim
!>
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64) 
!>        The Fortran file number to print to. 
!>    \endverbatim
!>
!>    \param[in] Header
!>    \verbatim
!>        Header is Character(Len=*)
!>        The title to print along with Integral.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Top
!>    \verbatim
!>        Blank_At_Top is Logical,Optional
!>        = .True.:  print blank line above output
!>        = .False.: do not print blank line above output.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Bottom
!>    \verbatim
!>        Blank_At_Bottom is Logical,Optional
!>        = .True.:  print blank line below output
!>        = .False.: do not print blank line below output.
!>    \endverbatim
!>
!>    \param[in] FormatStr
!>    \verbatim
!>        FormatStr is Character(Len=*),Optional
!>        Formatting statement for printing.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M.M.F. Moraes
!>    \date 2026
!
      subroutine mqc_print_nef_vector(NEFVec,iOut,header, &
          blank_at_top,blank_at_bottom,formatStr)
!
        implicit none
        class(mqc_nef_vector)::NEFVec
        integer(kind=int64),intent(in)::iOut
        character(len=*),intent(in)::header
        character(len=64)::arrayType
        logical,intent(In),optional::blank_at_top,blank_at_bottom
        character(len=*),intent(in),optional::formatStr
        integer::i
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
        if(present(blank_at_top)) then
          if(blank_at_top) write(iout,1020)
        endif
        write(iout,1000) trim(header)
        do i = 1, NEFVec%NumSpace
          call mqc_print(NEFVec%Main_vector(i),iOut,'Sub-Space '//trim(num2char(i)))
        end do
        
      end subroutine mqc_print_nef_vector
!
!
!
      subroutine mqc_nef_vector_initialize(NEFVec,initial_map_matrix,scalar)
!
        implicit none
        class(mqc_nef_vector),intent(inOut)::NEFVec
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in)::initial_map_matrix
        Class(*),Optional::scalar
        integer::i
!       
        if(allocated(NEFVec%Main_vector)) deallocate(NEFVec%Main_vector)
        allocate(NEFVec%Main_vector(size(initial_map_matrix)))
        if(allocated(NEFVec%Size_vector_alpha)) deallocate(NEFVec%Size_vector_alpha)
        allocate(NEFVec%Size_vector_alpha(size(initial_map_matrix)))
        if(allocated(NEFVec%Size_vector_beta)) deallocate(NEFVec%Size_vector_beta)
        allocate(NEFVec%Size_vector_beta(size(initial_map_matrix)))
        NEFVec%NumSpace = size(initial_map_matrix)
        NEFVec%Map_Matrix = initial_map_matrix
        do i = 1, NEFVec%NumSpace 
          select case(initial_map_matrix(i)%type())
          case ('space')
            If(Present(Scalar)) then
              call NEFVec%Main_vector(i)%init(&
                     Size(initial_map_matrix(i),'alpha'),&
                     scalar=scalar)
            else
              call NEFVec%Main_vector(i)%init(&
                     Size(initial_map_matrix(i),'alpha'))
            end if
            NEFVec%Size_vector_alpha(i) = Size(initial_map_matrix(i),'alpha')
            NEFVec%Size_vector_beta(i) = Size(initial_map_matrix(i),'alpha')
          case ('spin')
            If(Present(Scalar)) then
              call NEFVec%Main_vector(i)%init(&
                     Size(initial_map_matrix(i),'alpha'),&
                     scalar=scalar,&
                     lengthBeta=Size(initial_map_matrix(i),'beta'))
            else
              call NEFVec%Main_vector(i)%init(&
                     Size(initial_map_matrix(i),'alpha'),&
                     lengthBeta=Size(initial_map_matrix(i),'beta'))
            end if
            NEFVec%Size_vector_alpha(i) = Size(initial_map_matrix(i),'alpha')
            NEFVec%Size_vector_beta(i) = Size(initial_map_matrix(i),'beta')
          case default
            call mqc_error_A('type not valid in mqc_nef_vector_initialize', 6, &
              'initial_map_matrix(i)%type()', initial_map_matrix(i)%type() )
          end select 
        end do
!
      end subroutine mqc_nef_vector_initialize
!
!     PROCEDURE MQC_Print_NEF_Matrix
!
!>    \brief <b> MQC_Print_NEF_Matrix is a subroutine used to print an MQC
!>    NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_NEF_Matrix is a subroutine used to print an MQC NEF Matrix object. 
!>    Blank_At_Top and Blank_At_Bottom are optional logical arguments to print blank
!>    lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        Integral is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix variable to be printed.
!>    \endverbatim
!>
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64) 
!>        The Fortran file number to print to. 
!>    \endverbatim
!>
!>    \param[in] Header
!>    \verbatim
!>        Header is Character(Len=*)
!>        The title to print along with Integral.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Top
!>    \verbatim
!>        Blank_At_Top is Logical,Optional
!>        = .True.:  print blank line above output
!>        = .False.: do not print blank line above output.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Bottom
!>    \verbatim
!>        Blank_At_Bottom is Logical,Optional
!>        = .True.:  print blank line below output
!>        = .False.: do not print blank line below output.
!>    \endverbatim
!>
!>    \param[in] FormatStr
!>    \verbatim
!>        FormatStr is Character(Len=*),Optional
!>        Formatting statement for printing.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M.M.F. Moraes
!>    \date 2026
!
      subroutine mqc_print_nef_matrix(NEFMat,iOut,header, &
          blank_at_top,blank_at_bottom,formatStr,blockPrint)
!
      implicit none
      class(mqc_nef_matrix)::NEFMat
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      logical,optional::blockPrint
      character(len=64)::arrayType
      type(mqc_scf_integral)::tmp_int
      logical,intent(In),optional::blank_at_top,blank_at_bottom
      character(len=*),intent(in),optional::formatStr
      integer::i,j,sub_i,sub_j,line_num,col_num
      logical::my_block_print
      Character(Len=256)::Line
!
      if (present(blockPrint)) then
        my_block_print = blockPrint
      else
        my_block_print = .true.
      end if
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
      if(present(blank_at_top)) then
        if(blank_at_top) write(iout,1020)
      endif
      write(iout,1000) trim(header)
      line_num = 0
      if (my_block_print) then
        write(6,*) 'Alpha-Alpha'
        do j = 1, NEFMat%NumSpace_row
          if (j.eq.1) then 
            Line = ' '
            do i = 1, NEFMat%NumSpace_col
              do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                Line = trim(Line)//'           '//num2char(NEFMat%Map_Matrix_Col(i)%at(sub_i,'alpha'))
              end do
            end do
            write(*,*) Line
          end if
          do sub_j = 1, NEFMat%Size_vector_alpha_row(j)
            line_num = line_num + 1
            col_num = 1 
            Line = num2char(NEFMat%Map_Matrix_Row(j)%at(sub_j,'alpha'))
            do i = 1, NEFMat%NumSpace_col
              tmp_int = NEFMat%Main_matrix(NEFMat%blockIndex([1,j],[2,i]))
              do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                if (tmp_int%at(sub_j,sub_i,'alpha').eq.0.0) then 
                  Line = trim(Line)//'     '//num2char(0.00000)
                else if (tmp_int%at(sub_j,sub_i,'alpha').lt.0.0) then 
                  Line = trim(Line)//'    '//num2char(tmp_int%at(sub_j,sub_i,'alpha'))
                else
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'alpha'))
                end if
                if (sub_i.ne.NEFMat%Size_vector_alpha_col(i)) Line = trim(Line)//' '
              end do
              if(i.ne.NEFMat%NumSpace_col) Line = trim(Line)//'|'
            end do
            write(6,*) Line 
          end do
          if(j.ne.NEFMat%NumSpace_row) &
             write(6,*) '----------------------------------------------------------------------&
                    &----------------------------------------------------------------------'    
        end do
        line_num = 0
        write(6,*) 'Beta-Beta'
        do j = 1, NEFMat%NumSpace_row
          do sub_j = 1, NEFMat%Size_vector_beta_row(j)
            line_num = line_num + 1
            col_num = 1
            Line = num2char(NEFMat%Map_Matrix_Row(j)%at(sub_j,'beta'))
            do i = 1, NEFMat%NumSpace_col
              tmp_int = NEFMat%Main_matrix(NEFMat%blockIndex([1,j],[2,i]))
              do sub_i = 1, NEFMat%Size_vector_beta_col(i)
                if (tmp_int%at(sub_j,sub_i,'beta').eq.0.0) then
                  Line = trim(Line)//'     '//num2char(0.00000)
                else if (tmp_int%at(sub_j,sub_i,'beta').lt.0.0) then
                  Line = trim(Line)//'    '//num2char(tmp_int%at(sub_j,sub_i,'beta'))
                else
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'beta'))
                end if
                if (sub_i.ne.NEFMat%Size_vector_beta_col(i)) Line = trim(Line)//' '
              end do
              if(i.ne.NEFMat%NumSpace_col) Line = trim(Line)//'|'
            end do
            write(6,*) trim(Line)
          end do
          if(j.ne.NEFMat%NumSpace_row) &
             write(6,*) '----------------------------------------------------------------------&
                    &----------------------------------------------------------------------'    
        end do
        line_num = 0
        write(6,*) 'Alpha-Beta'
        do j = 1, NEFMat%NumSpace_row
          do sub_j = 1, NEFMat%Size_vector_beta_row(j)
            line_num = line_num + 1
            col_num = 1
            Line = num2char(NEFMat%Map_Matrix_Row(j)%at(sub_j,'alpha'))
            do i = 1, NEFMat%NumSpace_col
              tmp_int = NEFMat%Main_matrix(NEFMat%blockIndex([1,j],[2,i]))
              do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                if (tmp_int%at(sub_j,sub_i,'alphabeta').eq.0.0) then
                  Line = trim(Line)//'     '//num2char(0.00000)
                else if (tmp_int%at(sub_j,sub_i,'alphabeta').lt.0.0) then
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'alphabeta'))
                else
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'alphabeta'))
                end if
                if (sub_i.ne.NEFMat%Size_vector_beta_col(i)) Line = trim(Line)//' '
              end do
              if(i.ne.NEFMat%NumSpace_col) Line = trim(Line)//'|'
            end do
            write(6,*) trim(Line)
          end do
          if(j.ne.NEFMat%NumSpace_row) &
             write(6,*) '----------------------------------------------------------------------&
                    &----------------------------------------------------------------------'    
        end do
        line_num = 0
        write(6,*) 'Beta-Alpha'
        do j = 1, NEFMat%NumSpace_row
          do sub_j = 1, NEFMat%Size_vector_alpha_row(j)
            line_num = line_num + 1
            col_num = 1
            Line = num2char(NEFMat%Map_Matrix_Row(j)%at(sub_j,'beta'))
            do i = 1, NEFMat%NumSpace_col
              tmp_int = NEFMat%Main_matrix(NEFMat%blockIndex([1,j],[2,i]))
              do sub_i = 1, NEFMat%Size_vector_beta_col(i)
                if (tmp_int%at(sub_j,sub_i,'betaalpha').eq.0.0) then
                  Line = trim(Line)//'     '//num2char(0.00000)
                else if (tmp_int%at(sub_j,sub_i,'betaalpha').lt.0.0) then
                  Line = trim(Line)//'    '//num2char(tmp_int%at(sub_j,sub_i,'betaalpha'))
                else
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'betaalpha'))
                end if
                if (sub_i.ne.NEFMat%Size_vector_alpha_col(i)) Line = trim(Line)//' '
              end do
              if(i.ne.NEFMat%NumSpace_col) Line = trim(Line)//'|'
            end do
            write(6,*) trim(Line)
          end do
          if(j.ne.NEFMat%NumSpace_row) &
             write(6,*) '----------------------------------------------------------------------&
                    &----------------------------------------------------------------------'
        end do
        write(6,*) ' '
      else
        do i = 1, NEFMat%NumSpace_col
          do j = 1, NEFMat%NumSpace_row
            call mqc_print(NEFMat%Main_matrix(NEFMat%blockIndex([1,j],[2,i])),iOut,&
                  'Sub-Space '//trim(num2char(j))//','//trim(num2char(i)))
          end do
        end do
      end if
! 
      end subroutine mqc_print_nef_matrix
!
!
!     PROCEDURE MQC_NEF_Matrix_Initialize
!
!>    \brief <b> MQC_NEF_Matrix_Initialize is a subroutine that initializes a
!>    NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_Matrix_Initialize is a subroutine that initializes a MQC_NEF_Matrix. Default
!>    element values are 0.0, or otherwise vector can be initialized with optional
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        Integral is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to initialize.
!>    \endverbatim
!>
!>    \param[in] initial_map_matrix_row
!>    \verbatim
!>        initial_map_matrix_row is a type(mqc_scf_eigenvalues),dimension(:),allocatable
!>        Maps each element from an original mqc_scf_integral basis into the nef basis.
!>    \endverbatim
!>
!>    \param[in] scalar
!>    \verbatim
!>        Scalar is Class(*),Optional
!>        Value to set each element of Vector. If not present,
!>        the value is set to 0.0. Can be of type integer, real,
!>        complex or MQC_Scalar.
!>    \endverbatim
!>
!>    \param[in] initial_map_matrix_col
!>    \verbatim
!>        initial_map_matrix_row is a type(mqc_scf_eigenvalues),dimension(:),allocatable,optional
!>        Maps each element from an original mqc_scf_integral basis into the nef basis on the 
!>        second axis (columns). If not present, the mapping is the same in both axes.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_matrix_initialize(NEFMat,initial_map_matrix_row,scalar,initial_map_matrix_col)
!
        implicit none
        class(mqc_nef_matrix),intent(inOut)::NEFMat
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in)::initial_map_matrix_row
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in),optional::initial_map_matrix_col
        type(mqc_scf_eigenvalues),dimension(:),allocatable::active_map_matrix_col
        Class(*),Optional::scalar
        logical::same_vec
        integer::i,j
       
        NEFMat%NumSpace_row = size(initial_map_matrix_row)
        NEFMat%Map_Matrix_row = initial_map_matrix_row
        if(present(initial_map_matrix_col)) then
          NEFMat%NumSpace_col = size(initial_map_matrix_col)
          active_map_matrix_col = initial_map_matrix_col
          NEFMat%Map_Matrix_col = initial_map_matrix_col
          same_vec = .false.
        else
          NEFMat%NumSpace_col = NEFMat%NumSpace_row
          active_map_matrix_col = initial_map_matrix_row
          NEFMat%Map_Matrix_col = initial_map_matrix_row
          same_vec = .true.
        end if
        if(allocated(NEFMat%Size_vector_alpha_row)) deallocate(NEFMat%Size_vector_alpha_row)
        if(allocated(NEFMat%Size_vector_alpha_col)) deallocate(NEFMat%Size_vector_alpha_col)
        if(allocated(NEFMat%Size_vector_beta_row)) deallocate(NEFMat%Size_vector_beta_row)
        if(allocated(NEFMat%Size_vector_beta_col)) deallocate(NEFMat%Size_vector_beta_col)
        if(allocated(NEFMat%Main_matrix)) deallocate(NEFMat%Main_matrix)
        allocate(NEFMat%Size_vector_alpha_row(NEFMat%NumSpace_row)) 
        allocate(NEFMat%Size_vector_beta_row(NEFMat%NumSpace_row)) 
        allocate(NEFMat%Size_vector_alpha_col(NEFMat%NumSpace_col)) 
        allocate(NEFMat%Size_vector_beta_col(NEFMat%NumSpace_col)) 
        
        allocate(NEFMat%Main_matrix(NEFMat%NumSpace_row*NEFMat%NumSpace_col))
        do i = 1, NEFMat%NumSpace_row
          NEFMat%Size_vector_alpha_row(i) = Size(initial_map_matrix_row(i),'alpha')
          NEFMat%Size_vector_beta_row(i) = Size(initial_map_matrix_row(i),'beta')
          do j = 1, NEFMat%NumSpace_col
            if(i.eq.1) then
              NEFMat%Size_vector_alpha_col(j) = Size(active_map_matrix_col(j),'alpha')
              NEFMat%Size_vector_beta_col(j) = Size(active_map_matrix_col(j),'beta')
            end if
            if (same_vec.and.i.eq.j) then
              If(Present(Scalar)) then
                call NEFMat%Main_matrix(NEFMat%blockIndex([1,i],[2,i]))%init(&
                       Size(initial_map_matrix_row(i),'alpha'),&
                       Size(initial_map_matrix_row(i),'beta'),&
                       scalar=scalar)
              else
                call NEFMat%Main_matrix(NEFMat%blockIndex([1,i],[2,i]))%init(&
                       Size(initial_map_matrix_row(i),'alpha'),&
                       Size(initial_map_matrix_row(i),'beta'))
              end if
            else
              If(Present(Scalar)) then
                call NEFMat%Main_matrix(NEFMat%blockIndex([1,i],[2,j]))%init(&
                       Size(initial_map_matrix_row(i),'alpha'),&
                       Size(initial_map_matrix_row(i),'beta'),&
                       scalar=scalar,&
                       nAlpha2=Size(active_map_matrix_col(j),'alpha'),&
                       nBeta2=Size(active_map_matrix_col(j),'beta'))
              else
                call NEFMat%Main_matrix(NEFMat%blockIndex([1,i],[2,j]))%init(&
                       Size(initial_map_matrix_row(i),'alpha'),&
                       Size(initial_map_matrix_row(i),'beta'),&
                       nAlpha2=Size(active_map_matrix_col(j),'alpha'),&
                       nBeta2=Size(active_map_matrix_col(j),'beta'))
              end if
            end if
          end do
        end do
      end subroutine mqc_nef_matrix_initialize
!
!
!     PROCEDURE MQC_NEF_Matrix_Import_Integral
!
!>    \brief <b> MQC_NEF_Matrix_Import_Integral is a subroutine that import the values 
!>    from a MQC SCF Integral object to an initiallized NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_Matrix_Import_Integral is a subroutine that import the values 
!>    from a MQC SCF Integral object to an initiallized NEF Matrix object
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to initialize.
!>    \endverbatim
!>
!>    \param[in] integral
!>    \verbatim
!>        integral is a type(mqc_scf_integral)
!>        A MQC SCF Integral object containing the NEF Matrix elements in the 
!>        original order (in general AO or MO)
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_matrix_import_integral(NEFMat,integral)
!
        implicit none
        class(mqc_nef_matrix),intent(inOut)::NEFMat
        type(mqc_scf_integral)::integral
        type(mqc_scf_integral)::tmp_integral
        type(mqc_vector)::tmp_alpha,tmp_beta
        integer::i,j,row,col
!
        do row = 1, NEFMat%NumSpace_row
          tmp_alpha = NEFMat%Map_Matrix_Row(row)%getBlock('alpha')
          tmp_beta = NEFMat%Map_Matrix_Row(row)%getBlock('beta')
          tmp_integral = &
             integral%orbitals('useStrings',&
              alphaOrbsIn = tmp_alpha%getArrayI(),&
               BetaOrbsIn = tmp_beta%getArrayI(),axis=1)
          do col = 1, NEFMat%NumSpace_col
            tmp_alpha = NEFMat%Map_Matrix_Col(col)%getBlock('alpha')
            tmp_beta = NEFMat%Map_Matrix_Col(col)%getBlock('beta')
            NEFMat%Main_matrix(NEFMat%blockIndex([1,row],[2,col]))=&
                tmp_integral%orbitals('useStrings',&
              alphaOrbsIn = tmp_alpha%getArrayI(),&
               BetaOrbsIn = tmp_beta%getArrayI(),axis=2)
          end do
        end do
!
      end subroutine mqc_nef_matrix_import_integral
!
!
!
!     PROCEDURE MQC_NEF_Matrix_Export_Integral
!
!>    \brief <b> MQC_NEF_Matrix_Export_Integral is a subroutine that export the values 
!>    to a MQC SCF Integral object from a NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_Matrix_Export_Integral is a subroutine that export the values 
!>    to a MQC SCF Integral object from a NEF Matrix object
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix.
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_export_integral(NEFMat) result(integral)
!
        implicit none
        class(mqc_nef_matrix),intent(inOut)::NEFMat
        type(mqc_scf_integral)::integral
        type(mqc_scf_integral)::tmp_integral
        type(mqc_vector)::tmp_alpha,tmp_beta
        integer::id,ia,ja,ib,jb,row,col,maxColAlpha,maxColBeta,maxRowAlpha,maxRowBeta
!
        maxColAlpha = 0
        maxColBeta  = 0
        maxRowAlpha = 0
        maxRowBeta  = 0
        do row = 1, NEFMat%NumSpace_row
          maxRowAlpha = maxRowAlpha + NEFMat%Size_vector_alpha_row(row)
          maxRowBeta = maxRowBeta + NEFMat%Size_vector_beta_row(row)
        end do
        do col = 1, NEFMat%NumSpace_col
          maxColAlpha = maxColAlpha + NEFMat%Size_vector_alpha_col(col)
          maxColBeta = maxColBeta + NEFMat%Size_vector_beta_col(col)
        end do
!
        call integral%init(maxRowAlpha,maxRowBeta,nAlpha2=maxColAlpha,nBeta2=maxColBeta)
        do row = 1, NEFMat%NumSpace_row
          do col = 1, NEFMat%NumSpace_col
            id = NEFMat%blockIndex([1,row],[2,col])
            do ia = 1, NEFMat%Size_vector_alpha_row(row)
              do ja = 1, NEFMat%Size_vector_alpha_col(col)
                call integral%put(NEFMat%Main_matrix(id)%at(ia,ja,'alpha'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'alpha'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'alpha'),'alpha')
              end do
              do ja = 1, NEFMat%Size_vector_beta_col(col)
                call integral%put(NEFMat%Main_matrix(id)%at(ia,ja,'alphabeta'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'alpha'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'beta'),'alphabeta')
              end do
            end do
            do ia = 1, NEFMat%Size_vector_beta_row(row)
              do ja = 1, NEFMat%Size_vector_alpha_col(col)
                call integral%put(NEFMat%Main_matrix(id)%at(ia,ja,'betaalpha'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'beta'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'alpha'),'betaalpha')
              end do
              do ja = 1, NEFMat%Size_vector_beta_col(col)
                call integral%put(NEFMat%Main_matrix(id)%at(ia,ja,'beta'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'beta'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'beta'),'beta')
              end do
            end do
          end do
        end do
!
      end function mqc_nef_matrix_export_integral
!
!
!     PROCEDURE MQC_NEF_Vector_Import_Vector
!
!>    \brief <b> MQC_NEF_Vector_Import_Vector is a subroutine that import the values 
!>    from a MQC Vector object to an initiallized NEF Vector object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    MQC_NEF_Vector_Import_Vector is a subroutine that import the values 
!>    from a MQC Vector object to an initiallized NEF Vector object
!>    \verbatim
!>       
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVec
!>    \verbatim
!>        NEFVec is Class(MQC_NEF_Vector)
!>        The MQC NEF Vector to initialize.
!>    \endverbatim
!>
!>    \param[in] vector
!>    \verbatim
!>        vector is a type(mqc_vector)
!>        A MQC Vector object containing the NEF Vector elements in the 
!>        original order (in general AO or MO)
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_vector_import_vector(NEFVec,vector)
!
        implicit none
        class(mqc_nef_vector),intent(inOut)::NEFVec
        type(mqc_vector)::vector
        type(mqc_vector)::tmp_alpha,tmp_beta
        integer::row,el
        type(mqc_scalar)::zero
!
        zero = 0.0
        do row = 1, NEFVec%NumSpace
          tmp_alpha = NEFVec%Map_Matrix(row)%getBlock('alpha')
          tmp_beta = NEFVec%Map_Matrix(row)%getBlock('beta')
          call NEFVec%Main_Vector(row)%init(size(tmp_alpha),zero,size(tmp_beta))
          do el = 1, size(tmp_alpha)
             call NEFVec%Main_Vector(row)%put(vector%at(tmp_alpha%at(el)),el,'alpha')
          end do
          do el = 1, size(tmp_beta)
             call NEFVec%Main_Vector(row)%put(vector%at(tmp_beta%at(el)),el,'beta')
          end do
        end do
!
      end subroutine mqc_nef_vector_import_vector
!     PROCEDURE MQC_NEFMatrix_NEFMatrix_Multiply
!       
!>    \brief <b> MQC_NEFMatrix_NEFMatrix_Multiply is used to multiply two MQC NEF Matrix
!>    type variables and returns an MQC NEF Matrix</b> 
!
!>    \par Purpose:
!     =============
!>      
!>    \verbatim
!>      
!>    MQC_NEFMatrix_NEFMatrix_Multiply is used to multiply two MQC NEF Matrix type 
!>    variables and returns an MQC NEF Matrix</b> 
!>
!>    \endverbatim
!       
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMAtA is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix type variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix)
!>        The second MQC NEF Matrix variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_nefmatrix_multiply(NEFMatrixA,NEFMatrixB) result(NEFMatrixOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatrixA,NEFMatrixB
      type(mqc_nef_matrix)::NEFMatrixOut
      integer(kind=int64)::rowA,colB,colA,final_index
!
      if(NEFMatrixA%NumSpace_col.ne.NEFMatrixB%NumSpace_row) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_nefmatrix_multiply',6,&
        'NEFMatA%NumSpace_col',NEFMatrixA%NumSpace_col,&
        'NEFMatB%NumSpace_row',NEFMatrixB%NumSpace_row)

      do colA = 1, NEFMatrixA%NumSpace_col
        if(NEFMatrixA%Size_vector_alpha_col(colA)&
           .ne.NEFMatrixB%Size_vector_alpha_row(colA).or. &
           NEFMatrixA%Size_vector_beta_col(colA)&
           .ne.NEFMatrixB%Size_vector_beta_row(colA)) &
          call mqc_error_i('NEF subspace dimensions are &
          &not conformable in mqc_nefmatrix_nefmatrix_multiply',6,&
          'Sub-space index (colA)', colA,&
          'NEFMatA%Size_vector_alpha_col(colA)',NEFMatrixA%Size_vector_alpha_col(colA),&
          'NEFMatB%Size_vector_alpha_row(colA)',NEFMatrixB%Size_vector_alpha_row(colA),&
          'NEFMatA%Size_vector_beta_col(colA)',NEFMatrixA%Size_vector_beta_col(colA),&
          'NEFMatB%Size_vector_beta_row(colA)',NEFMatrixB%Size_vector_beta_row(colA))
      end do

      call NEFMatrixOut%init(NEFMatrixA%Map_Matrix_Row,&
          initial_map_matrix_col=NEFMatrixB%Map_Matrix_Col)
      do rowA = 1, NEFMatrixA%NumSpace_row
        do colB = 1, NEFMatrixB%NumSpace_col
          final_index = NEFMatrixOut%blockIndex([1,rowA],[2,colB])
          do colA = 1, NEFMatrixA%NumSpace_col
            if (colA.ne.1) then
              NEFMatrixOut%Main_Matrix(final_index) = NEFMatrixOut%Main_Matrix(final_index) + &
                matmul(NEFMatrixA%Main_Matrix(NEFMatrixA%blockIndex([1,rowA],[2,colA])),&
                       NEFMatrixB%Main_Matrix(NEFMatrixB%blockIndex([1,colA],[2,colB])))
            else
              NEFMatrixOut%Main_Matrix(final_index) =  &
                matmul(NEFMatrixA%Main_Matrix(NEFMatrixA%blockIndex([1,rowA],[2,colA])),&
                       NEFMatrixB%Main_Matrix(NEFMatrixB%blockIndex([1,colA],[2,colB])))
            end if
          end do
        end do
      end do
!
      end function mqc_nefmatrix_nefmatrix_multiply
!       
!     PROCEDURE MQC_NEFMatrix_Integral_Multiply
!       
!>    \brief <b> MQC_NEFMatrix_Integral_Multiply is used to multiply a MQC NEF Matrix
!>    and MQC SCF Integral type variables and returns an MQC NEF Matrix</b> 
!
!>    \par Purpose:
!     =============
!>      
!>    \verbatim
!>      
!>    MQC_NEFMatrix_Integral_Multiply is used to multiply a MQC NEF Matrix
!>    and MQC SCF Integral type variables and returns an MQC NEF Matrix
!>
!>    \endverbatim
!       
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMAt is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix type variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] IntegralB
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The second MQC SCF Integral variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_integral_multiply(NEFMatA,IntegralB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_scf_integral),intent(in)::IntegralB
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::row
!
      if(NEFMatA%NumSpace_col.ne.1) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_integral_multiply',6,&
        'NEFMatA%NumSpace_col',NEFMatA%NumSpace_col)
     
      if(NEFMatA%Size_vector_alpha_col(1)&
         .ne.size(IntegralB,'alpha',1).or. &
         NEFMatA%Size_vector_beta_col(1)&
         .ne.size(IntegralB,'beta',1)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nefmatrix_integral_multiply',6,&
        'NEFMatA%Size_vector_alpha_col',NEFMatA%Size_vector_alpha_col(1),&
        'size(IntegralB,alpha,1)',size(IntegralB,'alpha',1),&
        'NEFMatA%Size_vector_beta_col',NEFMatA%Size_vector_beta_col(1),&
        'size(IntegralB,beta,1)',size(IntegralB,'beta',1))

      call NEFMatOut%init(NEFMatA%Map_Matrix_Row,&
          initial_map_matrix_col=NEFMatA%Map_Matrix_Col)
      do row = 1, NEFMatA%NumSpace_row
        NEFMatOut%Main_Matrix(row) =  matmul(NEFMatA%Main_Matrix(row),IntegralB)
      end do
!
      end function mqc_nefmatrix_integral_multiply
!       
!     PROCEDURE MQC_Integral_NEFMatrix_Multiply
!       
!>    \brief <b> MQC_Integral_NEFMatrix_Multiply is used to multiply a MQC NEF Matrix
!>    and MQC SCF Integral type variables and returns an MQC NEF Matrix</b> 
!
!>    \par Purpose:
!     =============
!>      
!>    \verbatim
!>      
!>    MQC_Integral_NEFMatrix_Multiply is used to multiply a MQC NEF Matrix
!>    and MQC SCF Integral type variables and returns an MQC NEF Matrix
!>
!>    \endverbatim
!       
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The second MQC SCF Integral variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMAt is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix type variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_integral_nefmatrix_multiply(IntegralA,NEFMatB) result(NEFMatOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::IntegralA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::col
!
      if(NEFMatB%NumSpace_row.ne.1) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_integral_nefmatrix_multiply',6,&
        'NEFMatA%NumSpace_row',NEFMatB%NumSpace_row)
 
      if(NEFMatB%Size_vector_alpha_row(1)&
         .ne.size(IntegralA,'alpha',2).or. & 
         NEFMatB%Size_vector_beta_row(1)&
         .ne.size(IntegralA,'beta',2)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_integral_nefmatrix_multiply',6,&
        'NEFMatB%Size_vector_alpha_row',NEFMatB%Size_vector_alpha_row(1),&
        'size(IntegralA,alpha,2)',size(IntegralA,'alpha',1),&
        'NEFMatB%Size_vector_beta_row',NEFMatB%Size_vector_beta_row(1),&
        'size(IntegralA,beta,2)',size(IntegralA,'beta',1))

      call NEFMatOut%init(NEFMatB%Map_Matrix_Row,&
          initial_map_matrix_col=NEFMatB%Map_Matrix_Col)
      do col = 1, NEFMatB%NumSpace_col
        NEFMatOut%Main_Matrix(col) = matmul(IntegralA,NEFMatB%Main_Matrix(col))
      end do
!
      end function mqc_integral_nefmatrix_multiply
!
!
!     PROCEDURE MQC_NEF_Vector_Sum
!
!>    \brief <b> MQC_NEF_Vector_Sum is used to sum two MQC NEF Vector type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Sum is used to sum two MQC NEF vector type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVecA
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The first MQC NEF Vector variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFVecB
!>    \verbatim               
!>        NEFVecB is type(mqc_nef_vector)
!>        The first MQC NEF Vector variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_sum(NEFVecA,NEFVecB) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA,NEFVecB
      type(mqc_nef_vector)::NEFVecOut
      integer::col
!
      if(NEFVecA%NumSpace.ne.NEFVecB%NumSpace) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_vector_sum',6,&
        'NEFVecA%NumSpace',NEFVecA%NumSpace,&
        'NEFVecB%NumSpace',NEFVecB%NumSpace)
!
      do col = 1, NEFVecA%NumSpace
        if (NEFVecA%Size_vector_alpha(col).ne.NEFVecB%Size_vector_alpha(col) .or. &
                NEFVecA%Size_vector_beta(col) .ne.NEFVecB%Size_vector_beta(col)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_vector_sum',6,&
          'Sub-space index (col)', col,&
          'NEFVecA%Size_vector_alpha(col)',NEFVecA%Size_vector_alpha(col),&
          'NEFVecB%Size_vector_alpha(col)',NEFVecB%Size_vector_alpha(col),&
          'NEFVecA%Size_vector_beta(col)',NEFVecA%Size_vector_beta(col),&
          'NEFVecB%Size_vector_beta(col)',NEFVecB%Size_vector_beta(col))
      end do
!
      call NEFVecOut%init(NEFVecA%Map_Matrix)
      do col = 1, NEFVecB%NumSpace
        NEFVecOut%Main_Vector(col) = NEFVecA%Main_Vector(col) + NEFVecB%Main_Vector(col)
      end do
!
      end function mqc_nef_vector_sum
!
!
!     PROCEDURE MQC_NEFMatrix_Sum
!
!>    \brief <b> MQC_NEFMatrix_Sum is used to sum two MQC NEF matrix type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Sum is used to sum two MQC NEF matrix type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim               
!>        NEFMatB is type(mqc_nef_matrix)
!>        The second MQC NEF Matrix variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_sum(NEFMatA,NEFMatB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA,NEFMatB
      type(mqc_nef_matrix)::NEFMatOut
      integer::col,row,id
!
      if(NEFMatA%NumSpace_col.ne.NEFMatB%NumSpace_col.or.&
         NEFMatA%NumSpace_row.ne.NEFMatB%NumSpace_row) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_sum',6,&
        'NEFMatA%NumSpace_col',NEFMatA%NumSpace_col,&
        'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col,&
        'NEFMatA%NumSpace_row',NEFMatA%NumSpace_row,&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)

      do col = 1, NEFMatA%NumSpace_col
        if(NEFMatA%Size_vector_alpha_col(col)&
           .ne.NEFMatB%Size_vector_alpha_col(col).or. &
               NEFMatA%Size_vector_beta_col(col)  &
           .ne.NEFMatB%Size_vector_beta_col(col)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_sum',6,&
          'Sub-space index (col)', col,&
          'NEFMatA%Size_vector_alpha_col(col)',NEFMatA%Size_vector_alpha_col(col),&
          'NEFMatB%Size_vector_alpha_col(col)',NEFMatB%Size_vector_alpha_col(col),&
          'NEFMatA%Size_vector_beta_col(col)',NEFMatA%Size_vector_beta_col(col),&
          'NEFMatB%Size_vector_beta_col(col)',NEFMatB%Size_vector_beta_col(col))
      end do
      do row = 1, NEFMatA%NumSpace_row
        if(NEFMatA%Size_vector_alpha_row(row)&
           .ne.NEFMatB%Size_vector_alpha_row(row).or. &
               NEFMatA%Size_vector_beta_row(row)  &
           .ne.NEFMatB%Size_vector_beta_row(row)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_sum',6,&
          'Sub-space index (row)', row,&
          'NEFMatA%Size_vector_alpha_row(row)',NEFMatA%Size_vector_alpha_row(row),&
          'NEFMatB%Size_vector_alpha_row(row)',NEFMatB%Size_vector_alpha_row(row),&
          'NEFMatA%Size_vector_beta_row(row)',NEFMatA%Size_vector_beta_row(row),&
          'NEFMatB%Size_vector_beta_row(row)',NEFMatB%Size_vector_beta_row(row))
      end do

      call NEFMatOut%init(NEFMatA%Map_Matrix_Row,&
          initial_map_matrix_col=NEFMatB%Map_Matrix_Col)
      do row = 1, NEFMatA%NumSpace_row
        do col = 1, NEFMatB%NumSpace_col
          id = NEFMatOut%blockIndex([1,row],[2,col])
          NEFMatOut%Main_Matrix(id) = NEFMatA%Main_Matrix(id)+NEFMatB%Main_Matrix(id)
        end do
      end do
!
      end function mqc_nefmatrix_sum
!
!
!     PROCEDURE MQC_NEF_R4_Sum
!
!>    \brief <b> MQC_NEF_R4_Sum is used to sum two MQC NEF R4 type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Sum is used to sum two MQC NEF R4 type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4)
!>        The first MQC NEF R4 variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4)
!>        The second MQC NEF R4 variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_sum(NEFR4A,NEFR4B) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      type(mqc_nef_r4)::NEFR4Out
      integer::id1,id2,id3,id4,id
!
      if(NEFR4A%NumSpace_1.ne.NEFR4B%NumSpace_1)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_sum',6,&
        'NEFR4A%NumSpace_1',NEFR4A%NumSpace_1,&
        'NEFR4B%NumSpace_1',NEFR4B%NumSpace_1)
      if(NEFR4A%NumSpace_2.ne.NEFR4B%NumSpace_2)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_sum',6,&
        'NEFR4A%NumSpace_2',NEFR4A%NumSpace_2,&
        'NEFR4B%NumSpace_2',NEFR4B%NumSpace_2)
      if(NEFR4A%NumSpace_3.ne.NEFR4B%NumSpace_3)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_sum',6,&
        'NEFR4A%NumSpace_3',NEFR4A%NumSpace_3,&
        'NEFR4B%NumSpace_3',NEFR4B%NumSpace_3)
      if(NEFR4A%NumSpace_4.ne.NEFR4B%NumSpace_4)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_sum',6,&
        'NEFR4A%NumSpace_4',NEFR4A%NumSpace_4,&
        'NEFR4B%NumSpace_4',NEFR4B%NumSpace_4)
     
      do id1 = 1, NEFR4A%NumSpace_1
        if(NEFR4A%Size_vector_alpha_1(id1)&
           .ne.NEFR4B%Size_vector_alpha_1(id1).or.&
           NEFR4A%Size_vector_beta_1(id1)&
           .ne.NEFR4B%Size_vector_beta_1(id1))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_sum',6,&
          'Sub-space index (id1)', id1,&
          'NEFR4A%Size_vector_alpha_1(id1)',NEFR4A%Size_vector_alpha_1(id1),&
          'NEFR4B%Size_vector_alpha_1(id1)',NEFR4B%Size_vector_alpha_1(id1),&
          'NEFR4A%Size_vector_beta_1(id1)',NEFR4A%Size_vector_beta_1(id1),&
          'NEFR4B%Size_vector_beta_1(id1)',NEFR4B%Size_vector_beta_1(id1))
      end do
      do id2 = 1, NEFR4A%NumSpace_2
        if(NEFR4A%Size_vector_alpha_2(id2)&
           .ne.NEFR4B%Size_vector_alpha_2(id2).or.&
           NEFR4A%Size_vector_beta_2(id2)&
           .ne.NEFR4B%Size_vector_beta_2(id2))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_sum',6,&
          'Sub-space index (id2)', id2,&
          'NEFR4A%Size_vector_alpha_2(id2)',NEFR4A%Size_vector_alpha_2(id2),&
          'NEFR4B%Size_vector_alpha_2(id2)',NEFR4B%Size_vector_alpha_2(id2),&
           'NEFR4A%Size_vector_beta_2(id2)', NEFR4A%Size_vector_beta_2(id2),&
           'NEFR4B%Size_vector_beta_2(id2)', NEFR4B%Size_vector_beta_2(id2))
      end do
      do id3 = 1, NEFR4A%NumSpace_3
        if(NEFR4A%Size_vector_alpha_3(id3)&
           .ne.NEFR4B%Size_vector_alpha_3(id3).or.&
           NEFR4A%Size_vector_beta_3(id3)&
           .ne.NEFR4B%Size_vector_beta_3(id3))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_sum',6,&
          'Sub-space index (id3)', id3,&
          'NEFR4A%Size_vector_alpha_3(id3)',NEFR4A%Size_vector_alpha_3(id3),&
          'NEFR4B%Size_vector_alpha_3(id3)',NEFR4B%Size_vector_alpha_3(id3),&
           'NEFR4A%Size_vector_beta_3(id3)', NEFR4A%Size_vector_beta_3(id3),&
           'NEFR4B%Size_vector_beta_3(id3)', NEFR4B%Size_vector_beta_3(id3))
      end do
      do id4 = 1, NEFR4A%NumSpace_4
        if(NEFR4A%Size_vector_alpha_4(id4)&
           .ne.NEFR4B%Size_vector_alpha_4(id4).or.&
           NEFR4A%Size_vector_beta_4(id4)&
           .ne.NEFR4B%Size_vector_beta_4(id4))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_sum',6,&
          'Sub-space index (id4)', id4,&
          'NEFR4A%Size_vector_alpha_4(id4)',NEFR4A%Size_vector_alpha_4(id4),&
          'NEFR4B%Size_vector_alpha_4(id4)',NEFR4B%Size_vector_alpha_4(id4),&
           'NEFR4A%Size_vector_beta_4(id4)', NEFR4A%Size_vector_beta_4(id4),&
           'NEFR4B%Size_vector_beta_4(id4)', NEFR4B%Size_vector_beta_4(id4))
      end do


      NEFR4Out = NEFR4A
      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              id = NEFR4Out%blockIndex([4,id4],[3,id3],[2,id2],[1,id1])
              NEFR4Out%Main_R4(id) = NEFR4A%Main_R4(id) + NEFR4B%Main_R4(id)
           end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_sum
!
!     PROCEDURE MQC_NEFMatrix_Conjugate
!
!>    \brief <b> MQC_NEFMatrix_Conjugate is used to compute the conjugate 
!>    of a MQC NEF Matrix type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Conjugate is used to compute the conjugate of 
!>    a MQC NEF Matrix type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix type variable to conjugate.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEFMatrix_conjugate(NEFMat) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::row,col
!
      call NEFMatOut%init(NEFMat%Map_Matrix_row,&
           initial_map_matrix_col=NEFMat%Map_Matrix_col)
      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          NEFMatOut%Main_Matrix(NEFMatOut%blockIndex([1,row],[2,col])) = &
            conjg(NEFMat%Main_Matrix(NEFMat%blockIndex([1,row],[2,col])))
        end do
      end do
!
      end function mqc_NEFMatrix_conjugate
!
!     PROCEDURE MQC_NEFMatrix_Conjugate_Transpose
!
!>    \brief <b> MQC_NEFMatrix_Conjugate_Transpose is used to compute the conjugate 
!>    transpose of a MQC NEF Matrix type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Conjugate_Transpose is used to compute the conjugate transpose of 
!>    a MQC NEF Matrix type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix type variable to conjugate transpose.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEFMatrix_conjugate_transpose(NEFMat) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::row,col
!
      call NEFMatOut%init(NEFMat%Map_Matrix_col,&
           initial_map_matrix_col=NEFMat%Map_Matrix_row)
      
      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          NEFMatOut%Main_Matrix(NEFMatOut%blockIndex([2,row],[1,col])) = &
            dagger(NEFMat%Main_Matrix(NEFMat%blockIndex([1,row],[2,col])))
        end do
      end do
!
      end function mqc_NEFMatrix_conjugate_transpose
!
!     PROCEDURE MQC_NEFR4_Conjugate_Transpose
!
!>    \brief <b> MQC_NEFR4_Conjugate_Transpose is used to compute the conjugate 
!>    transpose of a MQC NEF R4 type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFR4_Conjugate_Transpose is used to compute the conjugate transpose of 
!>    a MQC NEF R4 type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is type(mqc_nef_r4)
!>        The MQC NEF R4 type variable to conjugate transpose.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEFR4_conjugate_transpose(NEFR4) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4
      type(mqc_nef_r4)::NEFR4Out
      integer(kind=int64)::id1,id2,id3,id4,N1234
!
      NEFR4Out%NumSpace_1 = NEFR4%NumSpace_4
      NEFR4Out%NumSpace_2 = NEFR4%NumSpace_3
      NEFR4Out%NumSpace_3 = NEFR4%NumSpace_2
      NEFR4Out%NumSpace_4 = NEFR4%NumSpace_1
      N1234 = NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2*NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4
      allocate(NEFR4Out%Main_R4(N1234))
      NEFR4Out%Map_Matrix_1 = NEFR4%Map_Matrix_4
      NEFR4Out%Map_Matrix_2 = NEFR4%Map_Matrix_3
      NEFR4Out%Map_Matrix_3 = NEFR4%Map_Matrix_2
      NEFR4Out%Map_Matrix_4 = NEFR4%Map_Matrix_1
      NEFR4Out%Size_vector_alpha_1 = NEFR4%Size_vector_alpha_4
      NEFR4Out%Size_vector_alpha_2 = NEFR4%Size_vector_alpha_3
      NEFR4Out%Size_vector_alpha_3 = NEFR4%Size_vector_alpha_2
      NEFR4Out%Size_vector_alpha_4 = NEFR4%Size_vector_alpha_1
      NEFR4Out%Size_vector_beta_1 = NEFR4%Size_vector_beta_4
      NEFR4Out%Size_vector_beta_2 = NEFR4%Size_vector_beta_3
      NEFR4Out%Size_vector_beta_3 = NEFR4%Size_vector_beta_2
      NEFR4Out%Size_vector_beta_4 = NEFR4%Size_vector_beta_1

      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              NEFR4Out%Main_R4(NEFR4Out%blockIndex([4,id4],[3,id3],[2,id2],[1,id1])) = &
                dagger(NEFR4%Main_R4(NEFR4%blockIndex([1,id4],[2,id3],[3,id2],[4,id1])))
            end do
          end do
        end do
      end do
!
      end function mqc_NEFR4_conjugate_transpose
!
!
!     PROCEDURE MQC_NEFMatrix_Transpose
!
!>    \brief <b> MQC_NEFMatrix_Transpose is used to compute the 
!>    transpose of a MQC NEF Matrix type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Transpose is used to compute the transpose of 
!>    a MQC NEF Matrix type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix type variable to transpose.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEFMatrix_transpose(NEFMat) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::row,col
!
      call NEFMatOut%init(NEFMat%Map_Matrix_col,&
           initial_map_matrix_col=NEFMat%Map_Matrix_row)

      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          NEFMatOut%Main_Matrix(NEFMatOut%blockIndex([2,row],[1,col])) = &
            transpose(NEFMat%Main_Matrix(NEFMat%blockIndex([1,row],[2,col])))
        end do
      end do
!
      end function mqc_NEFMatrix_transpose
!
!     
!     PROCEDURE MQC_NEF_Scalar_Vector_Multiply
!       
!>    \brief <b> MQC_NEF_Scalar_Vector_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Vector type variable and returns an MQC NEF Vector </b> 
!       
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_NEF_Scalar_Vector_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Vector type variable and returns an MQC NEF Vector.
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] ScalarA
!>    \verbatim
!>        ScalarB class(*)
!>        The variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFVecB
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The MQC NEF Vector variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_scalar_nef_vector_multiply(scalarA,NEFVecB) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecB
      class(*),intent(in)::scalarA
      type(mqc_scalar)::scalar
      type(mqc_nef_vector)::NEFVecOut
      integer::row
!
      select type(scalarA)
      type is (integer)
        scalar = scalarA
      type is (real)
        scalar = scalarA
      type is (complex)
        scalar = scalarA
      type is (mqc_scalar)
        scalar = scalarA
      class default
        call mqc_error_I('scalarIn type not defined in mqc_NEF_Vector_scalar_multiply',6)
      end select
!
      NEFVecOut = NEFVecB
      do row = 1, NEFVecB%NumSpace
        NEFVecOut%Main_vector(row) = NEFVecB%Main_vector(row)*scalar
      end do
!
      end function mqc_scalar_nef_vector_multiply
!
!     
!     PROCEDURE MQC_NEF_Vector_Scalar_Multiply
!       
!>    \brief <b> MQC_NEF_Vector_Scalar_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Vector type variable and returns an MQC NEF Vector </b> 
!       
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_NEF_Vector_Scalar_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Vector type variable and returns an MQC NEF Vector.
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] NEFVecA
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The MQC NEF Vector variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] ScalarB
!>    \verbatim
!>        ScalarB class(*)
!>        The variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_scalar_multiply(NEFVecA,scalarB) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA
      class(*),intent(in)::scalarB
      type(mqc_scalar)::scalar
      type(mqc_nef_vector)::NEFVecOut
      integer::row
!
      select type(scalarB)
      type is (integer)
        scalar = scalarB
      type is (real)
        scalar = scalarB
      type is (complex)
        scalar = scalarB
      type is (mqc_scalar)
        scalar = scalarB
      class default
        call mqc_error_I('scalarIn type not defined in mqc_scalar_NEF_Vector_multiply',6)
      end select
!
      NEFVecOut = NEFVecA
      do row = 1, NEFVecA%NumSpace
        NEFVecOut%Main_vector(row) = scalar*NEFVecA%Main_vector(row)
      end do
!
      end function mqc_nef_vector_scalar_multiply
!
!     
!     PROCEDURE MQC_Scalar_NEFMatrix_Multiply
!       
!>    \brief <b> MQC_Scalar_NEFMatrix_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Matrix type variable and returns an MQC NEF Matrix</b> 
!       
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_Scalar_NEFMatrix_Multiply is used to multiply an MQC scalar with an MQC 
!>    NEF Matrix type variable and returns an MQC NEF Matrix
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] ScalarIn
!>    \verbatim
!>        ScalarIn class(*)
!>        The variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_scalar_nefmatrix_multiply(scalarIn,NEFMat) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      class(*),intent(in)::scalarIn
      type(mqc_scalar)::scalar
      type(mqc_nef_matrix)::NEFMatOut
      integer::col,row,id
!
      select type(scalarIn)
      type is (integer)
        scalar = scalarIn
      type is (real)
        scalar = scalarIn
      type is (complex)
        scalar = scalarIn
      type is (mqc_scalar)
        scalar = scalarIn
      class default
        call mqc_error_I('scalarIn type not defined in mqc_scalar_NEFMat_multiply',6)
      end select
!
      NEFMatOut = NEFMat
      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          id=NEFMat%blockIndex([1,row],[2,col])
          NEFMatOut%Main_matrix(id)=scalar*NEFMat%Main_matrix(id)
        end do
      end do
!
      end function mqc_scalar_nefmatrix_multiply
!
!     
!     PROCEDURE MQC_NEFMatrix_Scalar_Multiply
!       
!>    \brief <b> MQC_NEFMatrix_Scalar_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Matrix type variable and returns an MQC NEF Matrix</b> 
!       
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_NEFMatrix_Scalar_Multiply is used to multiply an MQC scalar with an MQC 
!>    NEF Matrix type variable and returns an MQC NEF Matrix
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] ScalarB
!>    \verbatim
!>        ScalarB class(*)
!>        The variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_scalar_multiply(NEFMat,scalarB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      class(*),intent(in)::scalarB
      type(mqc_scalar)::scalar
      type(mqc_nef_matrix)::NEFMatOut
      integer::col,row,id
!
      select type(scalarB)
      type is (integer)
        scalar = scalarB
      type is (real)
        scalar = scalarB
      type is (complex)
        scalar = scalarB
      type is (mqc_scalar)
        scalar = scalarB
      class default
        call mqc_error_I('scalarIn type not defined in mqc_scalar_NEFMat_multiply',6)
      end select
!
      NEFMatOut = NEFMat
      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          id = NEFMat%blockIndex([1,row],[2,col])
          NEFMatOut%Main_matrix(id)=NEFMat%Main_matrix(id)*scalar
        end do
      end do
!
      end function mqc_nefmatrix_scalar_multiply
!
!     
!     PROCEDURE MQC_Scalar_NEFR4_Multiply
!       
!>    \brief <b> MQC_Scalar_NEFR4_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF R4 type variable and returns an MQC NEF R4</b> 
!     
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_Scalar_NEFR4_Multiply is used to multiply an MQC scalar with an MQC 
!>    NEF R4 type variable and returns an MQC NEF R4
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] ScalarIn
!>    \verbatim
!>        ScalarIn class(*)
!>        The variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is type(mqc_nef_r4)
!>        The MQC NEF R4 variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_scalar_nefr4_multiply(scalarIn,NEFR4) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4
      class(*),intent(in)::scalarIn
      type(mqc_scalar)::scalar
      type(mqc_nef_r4)::NEFR4Out
      integer::id,id1,id2,id3,id4
!
      select type(scalarIn)
      type is (integer)
        scalar = scalarIn
      type is (real)
        scalar = scalarIn
      type is (complex)
        scalar = scalarIn
      type is (mqc_scalar)
        scalar = scalarIn
      class default
        call mqc_error_I('scalarIn type not defined in mqc_scalar_NEFR4_multiply',6)
      end select
!
      NEFR4Out = NEFR4
      do id1 = 1, NEFR4%NumSpace_1
        do id2 = 1, NEFR4%NumSpace_2
          do id3 = 1, NEFR4%NumSpace_3
            do id4 = 1, NEFR4%NumSpace_4
              id = NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              NEFR4Out%Main_R4(id) = scalar*NEFR4Out%Main_R4(id)
            end do
          end do
        end do
      end do
!
      end function mqc_scalar_nefr4_multiply
      function mqc_nefr4_scalar_multiply(NEFR4,scalarB) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4
      class(*),intent(in)::scalarB
      type(mqc_scalar)::scalar
      type(mqc_nef_r4)::NEFR4Out
      integer::id,id1,id2,id3,id4
!
      select type(scalarB)
      type is (integer)
        scalar = scalarB
      type is (real)
        scalar = scalarB
      type is (complex)
        scalar = scalarB
      type is (mqc_scalar)
        scalar = scalarB
      class default
        call mqc_error_I('scalarIn type not defined in mqc_NEFR4_scalar_multiply',6)
      end select
!
      NEFR4Out = NEFR4
      do id1 = 1, NEFR4%NumSpace_1
        do id2 = 1, NEFR4%NumSpace_2
          do id3 = 1, NEFR4%NumSpace_3
            do id4 = 1, NEFR4%NumSpace_4
              id = NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              NEFR4Out%Main_R4(id) = NEFR4Out%Main_R4(id)*scalar
            end do
          end do
        end do
      end do
!
      end function mqc_nefr4_scalar_multiply
!
!
!
!     PROCEDURE MQC_NEF_Vector_Difference
!
!>    \brief <b> MQC_NEF_Vector_Difference is used to subtract two MQC NEF Vector type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Difference is used to subtract two MQC NEF vector type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVecA
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The first MQC NEF Vector variable to subtract.
!>    \endverbatim
!>    
!>    \param[in] NEFVecB
!>    \verbatim               
!>        NEFVecB is type(mqc_nef_vector)
!>        The second MQC NEF Vector variable to subtract.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      function mqc_nef_vector_difference(NEFVecA,NEFVecB) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA,NEFVecB
      type(mqc_nef_vector)::NEFVecOut
      type(mqc_nef_vector)::tmpVectorAlpha,tmpVectorBeta
      integer::col
!
      if(NEFVecA%NumSpace.ne.NEFVecB%NumSpace) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_vector_difference',6,&
        'NEFVecA%NumSpace',NEFVecA%NumSpace,&
        'NEFVecB%NumSpace',NEFVecB%NumSpace)
!
        do col = 1, NEFVecA%NumSpace
          if (NEFVecA%Size_vector_alpha(col).ne.NEFVecB%Size_vector_alpha(col) .or. &
              NEFVecA%Size_vector_beta(col) .ne.NEFVecB%Size_vector_beta(col)) &
                call mqc_error_i('NEF subspace dimensions are not conformable &
                        &in mqc_nef_vector_difference',6,&
                        'Sub-space index (col)', col,&
                        'NEFVecA%Size_vector_alpha(col)',NEFVecA%Size_vector_alpha(col),&
                        'NEFVecB%Size_vector_alpha(col)',NEFVecB%Size_vector_alpha(col),&
                        'NEFVecA%Size_vector_beta(col)',NEFVecA%Size_vector_beta(col),&
                        'NEFVecB%Size_vector_beta(col)',NEFVecB%Size_vector_beta(col))      
        end do
!
      call NEFVecOut%init(NEFVecA%Map_Matrix)
      do col = 1, NEFVecB%NumSpace
        NEFVecOut%Main_Vector(col) = NEFVecA%Main_Vector(col) - NEFVecB%Main_Vector(col)
      end do
!
!
      end function mqc_nef_vector_difference
!

!
!     PROCEDURE MQC_NEFMatrix_Difference
!
!>    \brief <b> MQC_NEFMatrix_Difference is used to subtract two MQC NEF matrix type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Difference is used to subtract two MQC NEF matrix type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim               
!>        NEFMatB is type(mqc_nef_matrix)
!>        The second MQC NEF Matrix variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      function mqc_nefmatrix_difference(NEFMatA,NEFMatB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA,NEFMatB
      type(mqc_nef_matrix)::NEFMatOut
      integer::col,row,id
!
      if(NEFMatA%NumSpace_col.ne.NEFMatB%NumSpace_col.or.&
         NEFMatA%NumSpace_row.ne.NEFMatB%NumSpace_row) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_difference',6,&
        'NEFMatA%NumSpace_col',NEFMatA%NumSpace_col,&
        'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col,&
        'NEFMatA%NumSpace_row',NEFMatA%NumSpace_row,&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)

      do col = 1, NEFMatA%NumSpace_col
        if(NEFMatA%Size_vector_alpha_col(col)&
           .ne.NEFMatB%Size_vector_alpha_col(col).or. &
               NEFMatA%Size_vector_beta_col(col)  &
           .ne.NEFMatB%Size_vector_beta_col(col)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_difference',6,&
          'Sub-space index (col)', col,&
          'NEFMatA%Size_vector_alpha_col(col)',NEFMatA%Size_vector_alpha_col(col),&
          'NEFMatB%Size_vector_alpha_col(col)',NEFMatB%Size_vector_alpha_col(col),&
          'NEFMatA%Size_vector_beta_col(col)',NEFMatA%Size_vector_beta_col(col),&
          'NEFMatB%Size_vector_beta_col(col)',NEFMatB%Size_vector_beta_col(col))
      end do
      do row = 1, NEFMatA%NumSpace_row
        if(NEFMatA%Size_vector_alpha_row(row)&
           .ne.NEFMatB%Size_vector_alpha_row(row).or. &
               NEFMatA%Size_vector_beta_row(row)  &
           .ne.NEFMatB%Size_vector_beta_row(row)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_difference',6,&
          'Sub-space index (row)', row,&
          'NEFMatA%Size_vector_alpha_row(row)',NEFMatA%Size_vector_alpha_row(row),&
          'NEFMatB%Size_vector_alpha_row(row)',NEFMatB%Size_vector_alpha_row(row),&
          'NEFMatA%Size_vector_beta_row(row)',NEFMatA%Size_vector_beta_row(row),&
          'NEFMatB%Size_vector_beta_row(row)',NEFMatB%Size_vector_beta_row(row))
      end do

      call NEFMatOut%init(NEFMatA%Map_Matrix_Row,&
          initial_map_matrix_col=NEFMatB%Map_Matrix_Col)
      do row = 1, NEFMatA%NumSpace_row
        do col = 1, NEFMatB%NumSpace_col
          id = NEFMatOut%blockIndex([1,row],[2,col])
          NEFMatOut%Main_Matrix(id) = NEFMatA%Main_Matrix(id)-NEFMatB%Main_Matrix(id)
        end do
      end do
      end function mqc_nefmatrix_difference
!
!
!     PROCEDURE MQC_NEF_R4_Difference
!
!>    \brief <b> MQC_NEF_R4_Difference is used to subtract two MQC NEF R4 type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Difference is used to subtract two MQC NEF R4 type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4)
!>        The first MQC NEF R4 variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4)
!>        The second MQC NEF R4 variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      function mqc_nef_r4_difference(NEFR4A,NEFR4B) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      type(mqc_nef_r4)::NEFR4Out
      integer::id1,id2,id3,id4,id
!
      if(NEFR4A%NumSpace_1.ne.NEFR4B%NumSpace_1)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_difference',6,&
        'NEFR4A%NumSpace_1',NEFR4A%NumSpace_1,&
        'NEFR4B%NumSpace_1',NEFR4B%NumSpace_1)
      if(NEFR4A%NumSpace_2.ne.NEFR4B%NumSpace_2)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_difference',6,&
        'NEFR4A%NumSpace_2',NEFR4A%NumSpace_2,&
        'NEFR4B%NumSpace_2',NEFR4B%NumSpace_2)
      if(NEFR4A%NumSpace_3.ne.NEFR4B%NumSpace_3)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_difference',6,&
        'NEFR4A%NumSpace_3',NEFR4A%NumSpace_3,&
        'NEFR4B%NumSpace_3',NEFR4B%NumSpace_3)
      if(NEFR4A%NumSpace_4.ne.NEFR4B%NumSpace_4)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_difference',6,&
        'NEFR4A%NumSpace_4',NEFR4A%NumSpace_4,&
        'NEFR4B%NumSpace_4',NEFR4B%NumSpace_4)
     
      do id1 = 1, NEFR4A%NumSpace_1
        if(NEFR4A%Size_vector_alpha_1(id1)&
           .ne.NEFR4B%Size_vector_alpha_1(id1).or.&
           NEFR4A%Size_vector_beta_1(id1)&
           .ne.NEFR4B%Size_vector_beta_1(id1))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_difference',6,&
          'Sub-space index (id1)', id1,&
          'NEFR4A%Size_vector_alpha_1(id1)',NEFR4A%Size_vector_alpha_1(id1),&
          'NEFR4B%Size_vector_alpha_1(id1)',NEFR4B%Size_vector_alpha_1(id1),&
          'NEFR4A%Size_vector_beta_1(id1)',NEFR4A%Size_vector_beta_1(id1),&
          'NEFR4B%Size_vector_beta_1(id1)',NEFR4B%Size_vector_beta_1(id1))
      end do
      do id2 = 1, NEFR4A%NumSpace_2
        if(NEFR4A%Size_vector_alpha_2(id2)&
           .ne.NEFR4B%Size_vector_alpha_2(id2).or.&
           NEFR4A%Size_vector_beta_2(id2)&
           .ne.NEFR4B%Size_vector_beta_2(id2))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_difference',6,&
          'Sub-space index (id2)', id2,&
          'NEFR4A%Size_vector_alpha_2(id2)',NEFR4A%Size_vector_alpha_2(id2),&
          'NEFR4B%Size_vector_alpha_2(id2)',NEFR4B%Size_vector_alpha_2(id2),&
           'NEFR4A%Size_vector_beta_2(id2)', NEFR4A%Size_vector_beta_2(id2),&
           'NEFR4B%Size_vector_beta_2(id2)', NEFR4B%Size_vector_beta_2(id2))
      end do
      do id3 = 1, NEFR4A%NumSpace_3
        if(NEFR4A%Size_vector_alpha_3(id3)&
           .ne.NEFR4B%Size_vector_alpha_3(id3).or.&
           NEFR4A%Size_vector_beta_3(id3)&
           .ne.NEFR4B%Size_vector_beta_3(id3))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_difference',6,&
          'Sub-space index (id3)', id3,&
          'NEFR4A%Size_vector_alpha_3(id3)',NEFR4A%Size_vector_alpha_3(id3),&
          'NEFR4B%Size_vector_alpha_3(id3)',NEFR4B%Size_vector_alpha_3(id3),&
           'NEFR4A%Size_vector_beta_3(id3)', NEFR4A%Size_vector_beta_3(id3),&
           'NEFR4B%Size_vector_beta_3(id3)', NEFR4B%Size_vector_beta_3(id3))
      end do
      do id4 = 1, NEFR4A%NumSpace_4
        if(NEFR4A%Size_vector_alpha_4(id4)&
           .ne.NEFR4B%Size_vector_alpha_4(id4).or.&
           NEFR4A%Size_vector_beta_4(id4)&
           .ne.NEFR4B%Size_vector_beta_4(id4))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_difference',6,&
          'Sub-space index (id4)', id4,&
          'NEFR4A%Size_vector_alpha_4(id4)',NEFR4A%Size_vector_alpha_4(id4),&
          'NEFR4B%Size_vector_alpha_4(id4)',NEFR4B%Size_vector_alpha_4(id4),&
           'NEFR4A%Size_vector_beta_4(id4)', NEFR4A%Size_vector_beta_4(id4),&
           'NEFR4B%Size_vector_beta_4(id4)', NEFR4B%Size_vector_beta_4(id4))
      end do


      NEFR4Out = NEFR4A
      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              id = NEFR4Out%blockIndex([4,id4],[3,id3],[2,id2],[1,id1])
              NEFR4Out%Main_R4(id) = NEFR4A%Main_R4(id) - NEFR4B%Main_R4(id)
           end do
          end do
        end do
      end do
      end function mqc_nef_r4_difference
!
!
!     PROCEDURE MQC_Print_NEF_R4
!
!>    \brief <b> MQC_Print_NEF_R4 is a subroutine used to print an MQC
!>    NEF R4 Tensor object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_NEF_R4 is a subroutine used to print an MQC NEF R4 Tensor object. 
!>    Blank_At_Top and Blank_At_Bottom are optional logical arguments to print blank
!>    lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        Integral is Class(MQC_NEF_R4)
!>        The MQC NEF R4 tensor variable to be printed.
!>    \endverbatim
!>
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64) 
!>        The Fortran file number to print to. 
!>    \endverbatim
!>
!>    \param[in] Header
!>    \verbatim
!>        Header is Character(Len=*)
!>        The title to print along with Integral.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Top
!>    \verbatim
!>        Blank_At_Top is Logical,Optional
!>        = .True.:  print blank line above output
!>        = .False.: do not print blank line above output.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Bottom
!>    \verbatim
!>        Blank_At_Bottom is Logical,Optional
!>        = .True.:  print blank line below output
!>        = .False.: do not print blank line below output.
!>    \endverbatim
!>
!>    \param[in] FormatStr
!>    \verbatim
!>        FormatStr is Character(Len=*),Optional
!>        Formatting statement for printing.
!>    \endverbatim
!>
!>    \param[in] Short_Print_Threshold
!>    \verbatim
!>        Short_Print_Threshold is Real,Optional
!>        if present: print only elements with absolute value larger than Short_Print_Threshold
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M.M.F. Moraes
!>    \date 2026
!
      subroutine mqc_print_nef_r4(NEFR4,iOut,header, &
          blank_at_top,blank_at_bottom,formatStr,short_print_threshold)
!
      implicit none
      class(mqc_nef_r4)::NEFR4
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      character(len=64)::arrayType
      logical,intent(In),optional::blank_at_top,blank_at_bottom
      character(len=*),intent(in),optional::formatStr
      Real(kind=real64),intent(in),optional::short_print_threshold
      Real(kind=real64)::print_threshold
      integer::id,id1,id2,id3,id4
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
!  
      if(present(blank_at_top)) then
        if(blank_at_top) write(iout,1020)
      endif
      If(Present(short_print_threshold)) then
        print_threshold = short_print_threshold
      else
        print_threshold = -1.0
      EndIf

      write(iout,1000) trim(header)
      do id1 = 1, NEFR4%NumSpace_1
        do id2 = 1, NEFR4%NumSpace_2
          do id3 = 1, NEFR4%NumSpace_3
            do id4 = 1, NEFR4%NumSpace_4
              id = NEFR4%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              call mqc_print(NEFR4%Main_R4(id),iOut,&
   'Sub-Space ('//trim(num2char(id1))//trim(num2char(id2))//'|'//trim(num2char(id3))//trim(num2char(id4))//') '&
               &' '//trim(NEFR4%Main_R4(id)%storType())//' '//trim(NEFR4%Main_R4(id)%type()),&
              short_print_threshold=print_threshold)
            end do
          end do
        end do
      end do
        
      end subroutine mqc_print_nef_r4
!
!     
!     
!     PROCEDURE MQC_NEFMatrix_NEFMatrix_Outer_Product
!     
!>    \brief <b> MQC_NEFMatrix_NEFMatrix_Outer_Product is used to calculate the
!>    outer product between two MQC NEF Matrix type variabel and returns an MQC NEF R4</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_NEFMatrix_Outer_Product is used to calculate the outer product between
!>    two MQC NEF Matrix and returns an MQC NEF R4. 
!>
!>    \endverbatim
!     
!     Arguments: 
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix)
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix)
!>    \endverbatim
!
!         
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!     
      function mqc_nefmatrix_nefmatrix_outer_product(NEFMatA,NEFMatB) result(NEFR4Out)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_r4)::NEFR4Out

      integer::idA,idB,id,id1,id2,id3,id4,N1,N12,N123,N1234

      NEFR4Out%NumSpace_2 = NEFMatA%NumSpace_col
      NEFR4Out%NumSpace_1 = NEFMatA%NumSpace_row
      NEFR4Out%NumSpace_4 = NEFMatB%NumSpace_col
      NEFR4Out%NumSpace_3 = NEFMatB%NumSpace_row
      write(*,'(A,I3,I3,I3,I3)') 'Sizes: ',NEFR4Out%NumSpace_1,NEFR4Out%NumSpace_2,&
                                           NEFR4Out%NumSpace_3,NEFR4Out%NumSpace_4
      N1 = NEFR4Out%NumSpace_4
      N12 = N1*NEFR4Out%NumSpace_3
      N123 = N12*NEFR4Out%NumSpace_2
      N1234 = N123*NEFR4Out%NumSpace_1
      write(*,'(A,I5)') 'Total size: ',N1234

      allocate(NEFR4Out%Main_R4(N1234))
      NEFR4Out%Map_Matrix_2 = NEFMatA%Map_Matrix_Col
      NEFR4Out%Map_Matrix_1 = NEFMatA%Map_Matrix_Row
      NEFR4Out%Map_Matrix_4 = NEFMatB%Map_Matrix_Col
      NEFR4Out%Map_Matrix_3 = NEFMatB%Map_Matrix_Row
      NEFR4Out%Size_vector_alpha_2 = NEFMatA%Size_vector_alpha_Col
      NEFR4Out%Size_vector_alpha_1 = NEFMatA%Size_vector_alpha_Row
      NEFR4Out%Size_vector_alpha_4 = NEFMatB%Size_vector_alpha_Col
      NEFR4Out%Size_vector_alpha_3 = NEFMatB%Size_vector_alpha_Row
      NEFR4Out%Size_vector_beta_2 = NEFMatA%Size_vector_beta_Col
      NEFR4Out%Size_vector_beta_1 = NEFMatA%Size_vector_beta_Row
      NEFR4Out%Size_vector_beta_4 = NEFMatB%Size_vector_beta_Col
      NEFR4Out%Size_vector_beta_3 = NEFMatB%Size_vector_beta_Row

      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          idA = NEFMatA%blockIndex([1,id1],[2,id2])
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              idB = NEFMatB%blockIndex([1,id3],[2,id4])
              NEFR4Out%Main_R4(NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])) = &
                NEFMatA%Main_Matrix(idA)%outer2ERI(NEFMatB%Main_Matrix(idB))
            end do
          end do
        end do
      end do
!     
      end function mqc_nefmatrix_nefmatrix_outer_product
!
!     
!     
!     PROCEDURE MQC_NEFMatrix_NEFMatrix_Crossouter_Product
!     
!>    \brief <b> MQC_NEFMatrix_NEFMatrix_crossouter_Product is used to calculate the
!>    crossouter product between two MQC NEF Matrix type variabel and returns an MQC NEF R4</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_NEFMatrix_Crossouter_Product is used to calculate the crossouter product between
!>    two MQC NEF Matrix and returns an MQC NEF R4. 
!>
!>    \endverbatim
!     
!     Arguments: 
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix)
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix)
!>    \endverbatim
!
!         
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!     
      function mqc_nefmatrix_nefmatrix_crossouter_product(NEFMatA,NEFMatB) result(NEFR4Out)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_r4)::NEFR4Out

      integer::idA,idB,id,id1,id2,id3,id4,N1,N12,N123,N1234

      NEFR4Out%NumSpace_1 = NEFMatA%NumSpace_col
      NEFR4Out%NumSpace_4 = NEFMatA%NumSpace_row
      NEFR4Out%NumSpace_3 = NEFMatB%NumSpace_col
      NEFR4Out%NumSpace_2 = NEFMatB%NumSpace_row
      N1 = NEFR4Out%NumSpace_4
      N12 = N1*NEFR4Out%NumSpace_3
      N123 = N12*NEFR4Out%NumSpace_2
      N1234 = N123*NEFR4Out%NumSpace_1

      allocate(NEFR4Out%Main_R4(N1234))
      NEFR4Out%Map_Matrix_1 = NEFMatA%Map_Matrix_Col
      NEFR4Out%Map_Matrix_4 = NEFMatA%Map_Matrix_Row
      NEFR4Out%Map_Matrix_3 = NEFMatB%Map_Matrix_Col
      NEFR4Out%Map_Matrix_2 = NEFMatB%Map_Matrix_Row
      NEFR4Out%Size_vector_alpha_1 = NEFMatA%Size_vector_alpha_Col
      NEFR4Out%Size_vector_alpha_4 = NEFMatA%Size_vector_alpha_Row
      NEFR4Out%Size_vector_alpha_3 = NEFMatB%Size_vector_alpha_Col
      NEFR4Out%Size_vector_alpha_2 = NEFMatB%Size_vector_alpha_Row
      NEFR4Out%Size_vector_beta_1 = NEFMatA%Size_vector_beta_Col
      NEFR4Out%Size_vector_beta_4 = NEFMatA%Size_vector_beta_Row
      NEFR4Out%Size_vector_beta_3 = NEFMatB%Size_vector_beta_Col
      NEFR4Out%Size_vector_beta_2 = NEFMatB%Size_vector_beta_Row

      do id1 = 1, NEFR4Out%NumSpace_1
        do id4 = 1, NEFR4Out%NumSpace_4
          idA = NEFMatA%blockIndex([1,id1],[2,id4])
          do id3 = 1, NEFR4Out%NumSpace_3
            do id2 = 1, NEFR4Out%NumSpace_2
              idB = NEFMatB%blockIndex([1,id3],[2,id2])
              NEFR4Out%Main_R4(NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])) = &
                NEFMatA%Main_Matrix(idA)%crossouter2ERI(NEFMatB%Main_Matrix(idB))
            end do
          end do
        end do
      end do
!     
      end function mqc_nefmatrix_nefmatrix_crossouter_product
!
!
!     PROCEDURE MQC_NEFMatrix_NEFMatrix_Contraction
!
!
!>    \brief <b> MQC_NEFMatrix_NEFMatrix_Contraction is used to return the contraction of two 
!>    MQC NEF Matrix objects.</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim 
!>
!>    MQC_NEFMatrix_NEFMatrix_Contraction is used to return the contraction of two 
!>    MQC NEF Matrix objects.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix) 
!>        The first MQC NEF Matrix to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix) 
!>        The second MQC NEF Matrix to be contracted.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_nefmatrix_contraction(NEFMatA,NEFMatB) result(output)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA,NEFMatB
      type(mqc_scalar)::output
      integer::row,col,id
!
      if(NEFMatA%NumSpace_col.ne.NEFMatB%NumSpace_col.or.&
         NEFMatA%NumSpace_row.ne.NEFMatB%NumSpace_row) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_sum',6,&
        'NEFMatA%NumSpace_col',NEFMatA%NumSpace_col,&
        'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col,&
        'NEFMatA%NumSpace_row',NEFMatA%NumSpace_row,&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)

      do col = 1, NEFMatA%NumSpace_col
        if(NEFMatA%Size_vector_alpha_col(col)&
           .ne.NEFMatB%Size_vector_alpha_col(col).or. &
               NEFMatA%Size_vector_beta_col(col)  &
           .ne.NEFMatB%Size_vector_beta_col(col)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_sum',6,&
          'Sub-space index (col)', col,&
          'NEFMatA%Size_vector_alpha_col(col)',NEFMatA%Size_vector_alpha_col(col),&
          'NEFMatB%Size_vector_alpha_col(col)',NEFMatB%Size_vector_alpha_col(col),&
          'NEFMatA%Size_vector_beta_col(col)',NEFMatA%Size_vector_beta_col(col),&
          'NEFMatB%Size_vector_beta_col(col)',NEFMatB%Size_vector_beta_col(col))
      end do
      do row = 1, NEFMatA%NumSpace_row
        if(NEFMatA%Size_vector_alpha_row(row)&
           .ne.NEFMatB%Size_vector_alpha_row(row).or. &
               NEFMatA%Size_vector_beta_row(row)  &
           .ne.NEFMatB%Size_vector_beta_row(row)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_sum',6,&
          'Sub-space index (row)', row,&
          'NEFMatA%Size_vector_alpha_row(row)',NEFMatA%Size_vector_alpha_row(row),&
          'NEFMatB%Size_vector_alpha_row(row)',NEFMatB%Size_vector_alpha_row(row),&
          'NEFMatA%Size_vector_beta_row(row)',NEFMatA%Size_vector_beta_row(row),&
          'NEFMatB%Size_vector_beta_row(row)',NEFMatB%Size_vector_beta_row(row))
      end do
      output = 0.0
      do row = 1, NEFMatA%NumSpace_row
        do col = 1, NEFMatA%NumSpace_col
          id = NEFMatA%blockIndex([1,row],[2,col])
          output = output + Contraction(NEFMatA%Main_Matrix(id),NEFMatB%Main_Matrix(id))
        end do
      end do
!
      end function mqc_nefmatrix_nefmatrix_contraction
!
!
!     PROCEDURE MQC_NEF_R4_NEF_R4_Contraction
!
!
!>    \brief <b> MQC_NEF_R4_NEF_R4_Contraction is used to return the contraction of two 
!>    MQC NEF R4 objects.</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_NEF_R4_Contraction is used to return the contraction of two 
!>    MQC NEF R4 objects.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4) 
!>        The second MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is character(len=*),optional 
!>        = 'coulomb':   contracts with coulomb integrals.
!>        = 'exchange':  contracts with exchange integrals.
!>        = 'doublebar': contracts with double bar integrals
!>                       (default).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_nef_r4_contraction(NEFR4A,NEFR4B,label) result(output)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      character(len=*),optional,intent(in)::label
      character(len=64)::myLabel
      type(mqc_scalar)::output
      integer::id1,id2,id3,id4,id,idex
      type(mqc_twoeris):: debbug_ERI
!
      if(present(label)) then
        call string_change_case(label,'L',myLabel)
      else
        myLabel = 'doublebar'
      endIf
!
      if(myLabel.ne.'coulomb'.and.myLabel.ne.'exchange'.and.&
         myLabel.ne.'doublebar') &
        call mqc_error_a('Unrecognized 2ERI integral interaction specified&
                        &in mqc_nef_r4_nef_r4_contraction.',&
                          6,'myLabel',myLabel)
!
      if(NEFR4A%NumSpace_1.ne.NEFR4B%NumSpace_1)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
        'NEFR4A%NumSpace_1',NEFR4A%NumSpace_1,&
        'NEFR4B%NumSpace_1',NEFR4B%NumSpace_1)
      if(NEFR4A%NumSpace_2.ne.NEFR4B%NumSpace_2)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
        'NEFR4A%NumSpace_2',NEFR4A%NumSpace_2,&
        'NEFR4B%NumSpace_2',NEFR4B%NumSpace_2)
      if(NEFR4A%NumSpace_3.ne.NEFR4B%NumSpace_3)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
        'NEFR4A%NumSpace_3',NEFR4A%NumSpace_3,&
        'NEFR4B%NumSpace_3',NEFR4B%NumSpace_3)
      if(NEFR4A%NumSpace_4.ne.NEFR4B%NumSpace_4)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
        'NEFR4A%NumSpace_4',NEFR4A%NumSpace_4,&
        'NEFR4B%NumSpace_4',NEFR4B%NumSpace_4)
     
      do id1 = 1, NEFR4A%NumSpace_1
        if(NEFR4A%Size_vector_alpha_1(id1)&
           .ne.NEFR4B%Size_vector_alpha_1(id1).or.&
           NEFR4A%Size_vector_beta_1(id1)&
           .ne.NEFR4B%Size_vector_beta_1(id1))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
          'Sub-space index (id1)', id1,&
          'NEFR4A%Size_vector_alpha_1(id1)',NEFR4A%Size_vector_alpha_1(id1),&
          'NEFR4B%Size_vector_alpha_1(id1)',NEFR4B%Size_vector_alpha_1(id1),&
          'NEFR4A%Size_vector_beta_1(id1)',NEFR4A%Size_vector_beta_1(id1),&
          'NEFR4B%Size_vector_beta_1(id1)',NEFR4B%Size_vector_beta_1(id1))
      end do
      do id2 = 1, NEFR4A%NumSpace_2
        if(NEFR4A%Size_vector_alpha_2(id2)&
           .ne.NEFR4B%Size_vector_alpha_2(id2).or.&
           NEFR4A%Size_vector_beta_2(id2)&
           .ne.NEFR4B%Size_vector_beta_2(id2))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
          'Sub-space index (id2)', id2,&
          'NEFR4A%Size_vector_alpha_2(id2)',NEFR4A%Size_vector_alpha_2(id2),&
          'NEFR4B%Size_vector_alpha_2(id2)',NEFR4B%Size_vector_alpha_2(id2),&
           'NEFR4A%Size_vector_beta_2(id2)', NEFR4A%Size_vector_beta_2(id2),&
           'NEFR4B%Size_vector_beta_2(id2)', NEFR4B%Size_vector_beta_2(id2))
      end do
      do id3 = 1, NEFR4A%NumSpace_3
        if(NEFR4A%Size_vector_alpha_3(id3)&
           .ne.NEFR4B%Size_vector_alpha_3(id3).or.&
           NEFR4A%Size_vector_beta_3(id3)&
           .ne.NEFR4B%Size_vector_beta_3(id3))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
          'Sub-space index (id3)', id3,&
          'NEFR4A%Size_vector_alpha_3(id3)',NEFR4A%Size_vector_alpha_3(id3),&
          'NEFR4B%Size_vector_alpha_3(id3)',NEFR4B%Size_vector_alpha_3(id3),&
           'NEFR4A%Size_vector_beta_3(id3)', NEFR4A%Size_vector_beta_3(id3),&
           'NEFR4B%Size_vector_beta_3(id3)', NEFR4B%Size_vector_beta_3(id3))
      end do
      do id4 = 1, NEFR4A%NumSpace_4
        if(NEFR4A%Size_vector_alpha_4(id4)&
           .ne.NEFR4B%Size_vector_alpha_4(id4).or.&
           NEFR4A%Size_vector_beta_4(id4)&
           .ne.NEFR4B%Size_vector_beta_4(id4))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
          'Sub-space index (id4)', id4,&
          'NEFR4A%Size_vector_alpha_4(id4)',NEFR4A%Size_vector_alpha_4(id4),&
          'NEFR4B%Size_vector_alpha_4(id4)',NEFR4B%Size_vector_alpha_4(id4),&
           'NEFR4A%Size_vector_beta_4(id4)', NEFR4A%Size_vector_beta_4(id4),&
           'NEFR4B%Size_vector_beta_4(id4)', NEFR4B%Size_vector_beta_4(id4))
      end do
!
      output = 0.0
      do id1 = 1, NEFR4A%NumSpace_1
        do id2 = 1, NEFR4A%NumSpace_2
          do id3 = 1, NEFR4A%NumSpace_3
            do id4 = 1, NEFR4A%NumSpace_4
              id = NEFR4A%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
                output = output + contraction(NEFR4A%Main_R4(id),NEFR4B%Main_R4(id),'coulomb')
              end if
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
                 idex = NEFR4A%blockIndex([1,id1],[4,id2],[3,id3],[2,id4])
                 output = output - &
                   contraction(NEFR4A%Main_R4(idex),NEFR4B%Main_R4(id)%swap(2,4),'coulomb')
              end if
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_nef_r4_contraction
!
!
!     PROCEDURE MQC_NEF_R4_NEF_Matrix_Contraction
!
!
!>    \brief <b> MQC_NEF_R4_NEF_Matrix_Contraction is used to return the contraction of a
!>    MQC NEF R4 and MQC NEF Matrix objects.</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_NEF_Matrix_Contraction is used to return the contraction of a
!>    MQC NEF R4 and MQC NEF Matrix objects.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix) 
!>        The MQC NEF Matrix to be contracted.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is character(len=*),optional 
!>        = 'coulomb':   contracts with coulomb integrals.
!>        = 'exchange':  contracts with exchange integrals.
!>        = 'doublebar': contracts with double bar integrals
!>                       (default).
!>    \endverbatim
!>
!>    \param[in] IndexA
!>    \verbatim
!>        IndexA is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIsA to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_nef_matrix_contraction(NEFR4A,NEFMatB,label,indexA) result(output)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A
      type(mqc_nef_matrix),intent(in)::NEFMatB
      character(len=*),optional,intent(in)::label
      integer(kind=int64),dimension(4),intent(in),optional::indexA
      integer(kind=int64),dimension(4)::my_indexA
      character(len=64)::myLabel
      type(mqc_nef_matrix)::output
      integer::id1,id2,id3,id4,idA,idB,idex,idOut
      type(mqc_twoeris):: debbug_ERI
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la
!
      if(present(label)) then
        call string_change_case(label,'L',myLabel)
      else
        myLabel = 'doublebar'
      endIf
!
      if(myLabel.ne.'coulomb'.and.myLabel.ne.'exchange'.and.&
         myLabel.ne.'doublebar') &
        call mqc_error_a('Unrecognized 2ERI integral interaction specified&
                        &in mqc_nef_r4_nef_r4_contraction.',&
                          6,'myLabel',myLabel)
!
      if(present(indexA)) then
        ia=indexA(1)
        ja=indexA(2)
        ka=indexA(3)
        la=indexA(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ka', ka,'la',la)
        my_indexA = indexA
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_indexA = [1,2,3,4]
      endif
!
      if(NEFR4A%blockSize(ia).ne.NEFMatB%NumSpace_row)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_Matrix_contraction',6,&
        'NEFR4A%blockSize(ia)',NEFR4A%blockSize(ia),&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)
      if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
        if(NEFR4A%blockSize(ja).ne.NEFMatB%NumSpace_col)&
          call mqc_error_i('NEF number of subspaces are &
          &not conformable in mqc_nef_r4_nef_Matrix_contraction',6,&
          'NEFR4A%blockSize(ja)',NEFR4A%blockSize(ja),&
          'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col)
      end if
      if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
        if(NEFR4A%blockSize(la).ne.NEFMatB%NumSpace_col)&
          call mqc_error_i('NEF number of subspaces are &
          &not conformable for exchange in mqc_nef_r4_nef_Matrix_contraction',6,&
          'NEFR4A%blockSize(la)',NEFR4A%blockSize(la),&
          'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col)
      end if
     
      do id1 = 1, NEFR4A%blockSize(ia)
        if(NEFR4A%blockSize(ia,id1,'alpha')&
           .ne.NEFMatB%Size_vector_alpha_row(id1).or.&
           NEFR4A%blockSize(ia,id1,'beta')&
           .ne.NEFMatB%Size_vector_beta_row(id1))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_Matrix_contraction',6,&
          'Sub-space index (id1)', id1,&
          'NEFR4A%blockSize(ia,id1,alpha)',NEFR4A%blockSize(ia,id1,'alpha'),&
          'NEFMatB%Size_vector_alpha_row(id1)',NEFMatB%Size_vector_alpha_row(id1),&
          'NEFR4A%blockSize(ia,id1,beta)',NEFR4A%blockSize(ia,id1,'beta'),&
          'NEFMatB%Size_vector_beta_row(id1)',NEFMatB%Size_vector_beta_row(id1))
      end do

      if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
        do id1 = 1, NEFR4A%blockSize(ja)
          if(NEFR4A%blockSize(ja,id1,'alpha')&
             .ne.NEFMatB%Size_vector_alpha_col(id1).or.&
             NEFR4A%blockSize(ja,id1,'beta')&
             .ne.NEFMatB%Size_vector_beta_col(id1))&
            call mqc_error_i('NEF subspace dimensions are &
          &not conformable in mqc_nef_r4_nef_Matrix_contraction',6,&
            'Sub-space index (id1)', id1,&
            'NEFR4A%blockSize(ja,id1,alpha)',NEFR4A%blockSize(ja,id1,'alpha'),&
            'NEFMatB%Size_vector_alpha_col(id1)',NEFMatB%Size_vector_alpha_col(id1),&
            'NEFR4A%blockSize(ja,id1,beta)',NEFR4A%blockSize(ja,id1,'beta'),&
            'NEFMatB%Size_vector_beta_col(id1)',NEFMatB%Size_vector_beta_col(id1))
        end do
      end if
      if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
        do id1 = 1, NEFR4A%blockSize(la)
          if(NEFR4A%blockSize(la,id1,'alpha')&
             .ne.NEFMatB%Size_vector_alpha_col(id1).or.&
             NEFR4A%blockSize(la,id1,'beta')&
             .ne.NEFMatB%Size_vector_beta_col(id1))&
            call mqc_error_i('NEF subspace dimensions are &
          &not conformable for exchange in mqc_nef_r4_nef_Matrix_contraction',6,&
            'Sub-space index (id1)', id1,&
            'NEFR4A%blockSize(la,id1,alpha)',NEFR4A%blockSize(la,id1,'alpha'),&
            'NEFMatB%Size_vector_alpha_col(id1)',NEFMatB%Size_vector_alpha_col(id1),&
            'NEFR4A%blockSize(la,id1,beta)',NEFR4A%blockSize(la,id1,'beta'),&
            'NEFMatB%Size_vector_beta_col(id1)',NEFMatB%Size_vector_beta_col(id1))
        end do
      end if
!
      Map_Matrix_Col = NEFR4A%blockMap(la)
      Map_Matrix_Row = NEFR4A%blockMap(ka)
      call output%init(Map_Matrix_Row,initial_map_matrix_col=Map_Matrix_Col)
      do id4 = 1, NEFR4A%blockSize(la)
        do id3 = 1, NEFR4A%blockSize(ka)
          do id2 = 1, NEFR4A%blockSize(ja)
            do id1 = 1, NEFR4A%blockSize(ia)
              idA = NEFR4A%blockIndex([ia,id1],[ja,id2],[ka,id3],[la,id4])
              idB = NEFMatB%blockIndex([1,id1],[2,id2])
              idOut = NEFMatB%blockIndex([1,id3],[2,id4])
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
                if ((NEFR4A%blockSize(ia,id1,'alpha')+NEFR4A%blockSize(ia,id1,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ja,id2,'alpha')+NEFR4A%blockSize(ja,id2,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ka,id3,'alpha')+NEFR4A%blockSize(ka,id3,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(la,id4,'alpha')+NEFR4A%blockSize(la,id4,'beta')).eq.0) cycle
                output%Main_Matrix(idOut) = output%Main_Matrix(idOut) &
                    + contraction(NEFR4A%Main_R4(idA),&
                                  NEFMatB%Main_Matrix(idB),[ia,ja])
              end if
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
                if ((NEFR4A%blockSize(ia,id1,'alpha')+NEFR4A%blockSize(ia,id1,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(la,id2,'alpha')+NEFR4A%blockSize(la,id2,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ka,id3,'alpha')+NEFR4A%blockSize(ka,id3,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ja,id4,'alpha')+NEFR4A%blockSize(ja,id4,'beta')).eq.0) cycle
                idex = NEFR4A%blockIndex([ia,id1],[la,id2],[ka,id3],[ja,id4])
                idB = NEFMatB%blockIndex([1,id1],[2,id2])
                idOut = NEFMatB%blockIndex([1,id3],[2,id4])
                output%Main_Matrix(idOut) = output%Main_Matrix(idOut) &
                    - contraction(NEFR4A%Main_R4(idex),&
                                  NEFMatB%Main_Matrix(idB),[ia,la])
              end if
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_nef_matrix_contraction
!
!
!     PROCEDURE MQC_NEF_Matrix_Dimension_Slice
!
!
!>    \brief <b> MQC_NEF_Matrix_Dimension_Slice is used to return a slice of
!>    the input NEF Matrix object </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Dimension_Slice is used to return a slice of
!>    the input NEF Matrix object.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix) 
!>        NEF Matrix Object to be Sliced.
!>    \endverbatim
!>
!>    \param[in] rows
!>    \verbatim
!>        rows is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the remaining rows after the slicing.
!>        Default: all rows
!>    \endverbatim
!>
!>    \param[in] columns
!>    \verbatim
!>        columns is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the remaining columns after the slicing.
!>        Default: all columns
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_dimention_slice(NEFMat,rows,columns) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64),dimension(:),optional,intent(in)::rows,columns
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map_Col,New_Map_Row
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,id,idOld
      logical::full_row,full_col
!
      if(present(rows)) then
        my_row = rows
        full_row=.false.
        if(present(columns)) then
          my_col = columns
          full_col=.false.
        else
          full_col=.true.
        endIf
      else
        full_row=.true.
        if(present(columns)) then
          my_col = columns
          full_col=.false.
        else
          NEFMatOut = NEFMat
          return
        endIf
      endIf
!
      if(full_col) then
        New_Map_Col = NEFMat%Map_Matrix_Col
      else
        allocate(New_Map_Col(size(my_col)))
        do id = 1, size(my_col)
          if(my_col(id).gt.NEFMat%NumSpace_col.or.my_col(id).le.0) &
        call mqc_error_i('Column index out-of-bound in mqc_nef_matrix_dimention_slice',6,&
             'my_col(id)',my_col(id),'NEFMat%NumSpace_col',NEFMat%NumSpace_col)
          New_Map_Col(id) = NEFMat%Map_Matrix_Col(my_col(id))
        end do
      end if
!
      if(full_row) then
        New_Map_Row = NEFMat%Map_Matrix_Row
      else
        allocate(New_Map_Row(size(my_row)))
        do id = 1, size(my_row)
          if(my_row(id).gt.NEFMat%NumSpace_row.or.my_row(id).le.0) &
        call mqc_error_i('Row index out-of-bound in mqc_nef_matrix_dimention_slice',6,&
             'my_row(id)',my_row(id),'NEFMat%NumSpace_row',NEFMat%NumSpace_row)
          New_Map_Row(id) = NEFMat%Map_Matrix_Row(my_row(id))
        end do
      end if
!
      call NEFMatOut%init(New_Map_Row,&
          initial_map_matrix_col=New_Map_Col)
      if(full_row) then
        do col = 1, size(my_col)
          do row = 1, size(NEFMat%Map_Matrix_Row)
            id = NEFMatOut%blockIndex([1,row],[2,col])
            idOld = NEFMat%blockIndex([1,row],[2,my_col(col)])
            NEFMatOut%Main_matrix(id) = NEFMat%Main_matrix(idOld) 
          end do
        end do
      else if (full_col) then
        do row = 1, size(my_row)
          do col = 1, size(NEFMat%Map_Matrix_Col)
            id = NEFMatOut%blockIndex([1,row],[2,col])
            idOld = NEFMat%blockIndex([1,my_row(row)],[2,col])
            NEFMatOut%Main_matrix(id) = NEFMat%Main_matrix(idOld) 
          end do
        end do
      else
        do row = 1, size(my_row)
          do col = 1, size(my_col)
            id = NEFMatOut%blockIndex([1,row],[2,col])
            idOld = NEFMat%blockIndex([1,my_row(row)],[2,my_col(col)])
            NEFMatOut%Main_matrix(id) = NEFMat%Main_matrix(idOld) 
          end do
        end do
      end if
!
      end function mqc_nef_matrix_dimention_slice
!
!
!     PROCEDURE MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4
!
!>    \brief <b> MQC_NEF_R4_NEF_R4_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC NEF R4 </b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_NEF_R4_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC NEF R4 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4) 
!>        The second MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] Index1
!>    \verbatim
!>        Index1 is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIs1 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>        Ex: If Index1=[3,2,4,1], the third and second indexes of ERIs1 are contracted.
!>            The fourth and first are multiplied to form the first and second indices
!>            of the result matrix, respectivelly.
!>    \endverbatim
!>
!>    \param[in] Index2
!>    \verbatim
!>        Index2 is Integer(kind=int64),dimension(4),optional
!>        Contain the indexes of ERIs2 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>
!>        Ex: index1=[3,2,4,1] and index2=[2,1,4,3]
!>        sum_{i,j} R4A_{pjiq} R4B_{jisr} = R4Out_{pqrs}
!>
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_nef_r4_partial_contraction_2nef_r4(NEFR4A,NEFR4B,index1,index2) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
      logical::first
!
      if(present(index1)) then
        ia=index1(1)
        ja=index1(2)
        ka=index1(3)
        la=index1(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ka', ka,'la',la)
        my_index1 = index1
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_index1 = [1,2,3,4]
      endif
      if(present(index2)) then
        ib=index2(1)
        jb=index2(2)
        kb=index2(3)
        lb=index2(4)
        if (ib.eq.jb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'jb',jb)
        if (ib.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'kb',kb)
        if (ib.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'lb',lb)
        if (jb.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'jb', jb,'kb',kb)
        if (jb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'jb', jb,'lb',lb)
        if (kb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'kb', kb,'lb',lb)
        my_index2 = index2
      else
        ib=1
        jb=2
        kb=3
        lb=4
        my_index2 = [1,2,3,4]
      endif
      if (ia.gt.4.or.ia.le.0.or.ib.gt.4.or.ib.le.0.or.&
          ja.gt.4.or.ja.le.0.or.jb.gt.4.or.jb.le.0.or.&
          ka.gt.4.or.ka.le.0.or.kb.gt.4.or.kb.le.0.or.&
          la.gt.4.or.la.le.0.or.lb.gt.4.or.lb.le.0.) Call MQC_Error(&
           'Index out of range in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4', 6)

      if (nefr4A%blockSize(ia).ne.nefr4B%blockSize(ib).or.&
          nefr4A%blockSize(ja).ne.nefr4B%blockSize(jb)) &
          call mqc_error_i('NEFR4 not correctly dimensioned in mqc_nef_r4_nef_r4_contraction_2R4',6,&
           'dim1',nefr4A%blockSize(ia)-nefr4B%blockSize(ib),&
           'dim2',nefr4A%blockSize(ja)-nefr4B%blockSize(jb))

      NEFR4Out%NumSpace_1 = NEFR4A%blockSize(ka)
      NEFR4Out%NumSpace_2 = NEFR4A%blockSize(la)
      NEFR4Out%NumSpace_3 = NEFR4B%blockSize(kb)
      NEFR4Out%NumSpace_4 = NEFR4B%blockSize(lb)
      N1234 = NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2*NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4
      write(*,*) 'Numb of entries = ',N1234

      allocate(NEFR4Out%Main_R4(N1234))
      NEFR4Out%Map_Matrix_1 = NEFR4A%blockMap(ka) 
      NEFR4Out%Map_Matrix_2 = NEFR4A%blockMap(la) 
      NEFR4Out%Map_Matrix_3 = NEFR4B%blockMap(kb) 
      NEFR4Out%Map_Matrix_4 = NEFR4B%blockMap(lb) 
      NEFR4Out%Size_vector_alpha_1 = NEFR4A%blockSizeVec(ka,'alpha') 
      NEFR4Out%Size_vector_alpha_2 = NEFR4A%blockSizeVec(la,'alpha') 
      NEFR4Out%Size_vector_alpha_3 = NEFR4B%blockSizeVec(kb,'alpha') 
      NEFR4Out%Size_vector_alpha_4 = NEFR4B%blockSizeVec(lb,'alpha') 
      NEFR4Out%Size_vector_beta_1 = NEFR4A%blockSizeVec(ka,'beta')  
      NEFR4Out%Size_vector_beta_2 = NEFR4A%blockSizeVec(la,'beta')  
      NEFR4Out%Size_vector_beta_3 = NEFR4B%blockSizeVec(kb,'beta')  
      NEFR4Out%Size_vector_beta_4 = NEFR4B%blockSizeVec(lb,'beta')  

      do id1 = 1, nefr4Out%blockSize(1)
        do id2 = 1, nefr4Out%blockSize(2)
          do id3 = 1, nefr4Out%blockSize(3)
            do id4 = 1, nefr4Out%blockSize(4)
              out_index = NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              first=.true.
              do idiab = 1, nefr4B%blockSize(ib)
                do idjab = 1, nefr4B%blockSize(jb)
                  A_index = NEFR4A%blockIndex([ka,id1],[la,id2],[ia,idiab],[ja,idjab])
                  B_index = NEFR4B%blockIndex([kb,id3],[lb,id4],[ib,idiab],[jb,idjab])
                  if(first) then
                    NEFR4Out%Main_R4(out_index) =  &
                      mqc_eri_eri_partial_contraction_2eri(&
                      NEFR4A%Main_R4(A_index),NEFR4B%Main_R4(B_index),my_index1,my_index2)
                    first=.false.
                  else
                    NEFR4Out%Main_R4(out_index) =  NEFR4Out%Main_R4(out_index)+&
                      mqc_eri_eri_partial_contraction_2eri(&
                      NEFR4A%Main_R4(A_index),NEFR4B%Main_R4(B_index),my_index1,my_index2)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_nef_r4_partial_contraction_2nef_r4
!
!
!     PROCEDURE MQC_NEF_R4_ERI_Partial_Contraction_2R4
!
!>    \brief <b> MQC_NEF_R4_ERI_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC 2ERI </b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_ERI_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC ERI
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] ERIB
!>    \verbatim
!>        ERIB is type(mqc_twoeris) 
!>        The second MQC 2ERI to be contracted.
!>    \endverbatim
!>
!>    \param[in] Index1
!>    \verbatim
!>        Index1 is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIs1 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>        Ex: If Index1=[3,2,4,1], the third and second indexes of ERIs1 are contracted.
!>            The fourth and first are multiplied to form the first and second indices
!>            of the result matrix, respectivelly.
!>    \endverbatim
!>
!>    \param[in] Index2
!>    \verbatim
!>        Index2 is Integer(kind=int64),dimension(4),optional
!>        Contain the indexes of ERIs2 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>
!>        Ex: index1=[3,2,4,1] and index2=[2,1,4,3]
!>        sum_{i,j} R4A_{pjiq} ERIB_{jisr} = R4Out_{pqrs}
!>
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_eri_partial_contraction_2nef_r4(NEFR4A,ERIB,index1,index2) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A
      type(mqc_twoeris),intent(in)::ERIB
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      type(mqc_scalar)::val
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
      logical::first
!
      if(present(index1)) then
        ia=index1(1)
        ja=index1(2)
        ka=index1(3)
        la=index1(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ka', ka,'la',la)
        my_index1 = index1
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_index1 = [1,2,3,4]
      endif
      if(present(index2)) then
        ib=index2(1)
        jb=index2(2)
        kb=index2(3)
        lb=index2(4)
        if (ib.eq.jb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ib', ib,'jb',jb)
        if (ib.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ib', ib,'kb',kb)
        if (ib.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ib', ib,'lb',lb)
        if (jb.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'jb', jb,'kb',kb)
        if (jb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'jb', jb,'lb',lb)
        if (kb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'kb', kb,'lb',lb)
        my_index2 = index2
      else
        ib=1
        jb=2
        kb=3
        lb=4
        my_index2 = [1,2,3,4]
      endif
      if (ia.gt.4.or.ia.le.0.or.ib.gt.4.or.ib.le.0.or.&
          ja.gt.4.or.ja.le.0.or.jb.gt.4.or.jb.le.0.or.&
          ka.gt.4.or.ka.le.0.or.kb.gt.4.or.kb.le.0.or.&
          la.gt.4.or.la.le.0.or.lb.gt.4.or.lb.le.0.) Call MQC_Error(&
           'Index out of range in MQC_NEF_R4_ERI_Partial_Contraction_2R4', 6)

      if (nefr4A%blockSize(ia).ne.1.or.nefr4A%blockSize(ja).ne.1) &
          call mqc_error_i('NEFR4 not correctly dimensioned in mqc_nef_r4_eri_contraction_2R4',6,&
           'dim1',nefr4A%blockSize(ia),'dim2',nefr4A%blockSize(ja))
      if(NEFR4A%blockSize(ia,1,'alpha').ne.size(ERIB,'alpha',ib).or. &
         NEFR4A%blockSize(ia,1,'beta').ne.size(ERIB,'beta',ib)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nefNEF_R4_ERI_Partial_Contraction_2R4',6,&
        'NEFR4A%blockSize(ia,1,alpha)',NEFR4A%blockSize(ia,1,'alpha'),&
        'size(ERIB,alpha,ib)',size(ERIB,'alpha',ib),&
        'NEFR4A%blockSize(ia,1,beta)',NEFR4A%blockSize(ia,1,'beta'),&
        'size(ERIB,beta,ib)',size(ERIB,'beta',ib))
      if(NEFR4A%blockSize(ja,1,'alpha').ne.size(ERIB,'alpha',jb).or. &
         NEFR4A%blockSize(ja,1,'beta').ne.size(ERIB,'beta',jb)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nefNEF_R4_ERI_Partial_Contraction_2R4',6,&
        'NEFR4A%blockSize(ja,1,alpha)',NEFR4A%blockSize(ja,1,'alpha'),&
        'size(ERIB,alpha,jb)',size(ERIB,'alpha',jb),&
        'NEFR4A%blockSize(ja,1,beta)',NEFR4A%blockSize(ja,1,'beta'),&
        'size(ERIB,beta,jb)',size(ERIB,'beta',jb))

      NEFR4Out%NumSpace_1 = NEFR4A%blockSize(ka)
      NEFR4Out%NumSpace_2 = NEFR4A%blockSize(la)
      NEFR4Out%NumSpace_3 = 1 
      NEFR4Out%NumSpace_4 = 1
      N1234 = NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2

      allocate(NEFR4Out%Main_R4(N1234))
      NEFR4Out%Map_Matrix_1 = NEFR4A%blockMap(ka) 
      NEFR4Out%Map_Matrix_2 = NEFR4A%blockMap(la) 
      if(allocated(NEFR4Out%Map_Matrix_3)) deallocate(NEFR4Out%Map_Matrix_3)
      allocate(NEFR4Out%Map_Matrix_3(1))
      if(allocated(NEFR4Out%Map_Matrix_4)) deallocate(NEFR4Out%Map_Matrix_4)
      allocate(NEFR4Out%Map_Matrix_4(1))
      call NEFR4Out%Map_Matrix_3(1)%init(size(ERIB,'alpha',kb),scalar=1,lengthBeta=size(ERIB,'beta',kb))
      call NEFR4Out%Map_Matrix_4(1)%init(size(ERIB,'alpha',lb),scalar=1,lengthBeta=size(ERIB,'beta',lb))
      do id1 = 1, size(ERIB,'alpha',kb)
        val = id1
        call  NEFR4Out%Map_Matrix_3(1)%put(val,id1,'alpha')
      end do
      do id1 = 1, size(ERIB,'beta',kb)
        val = id1
        call  NEFR4Out%Map_Matrix_3(1)%put(val,id1,'beta')
      end do
      do id1 = 1, size(ERIB,'alpha',lb)
        val = id1
        call  NEFR4Out%Map_Matrix_4(1)%put(val,id1,'alpha')
      end do
      do id1 = 1, size(ERIB,'beta',lb)
        val = id1
        call  NEFR4Out%Map_Matrix_4(1)%put(val,id1,'beta')
      end do
    
      NEFR4Out%Size_vector_alpha_1 = NEFR4A%blockSizeVec(ka,'alpha') 
      NEFR4Out%Size_vector_alpha_2 = NEFR4A%blockSizeVec(la,'alpha') 
      NEFR4Out%Size_vector_beta_1 = NEFR4A%blockSizeVec(ka,'beta')  
      NEFR4Out%Size_vector_beta_2 = NEFR4A%blockSizeVec(la,'beta')  
      if(allocated(NEFR4Out%Size_vector_alpha_3)) deallocate(NEFR4Out%Size_vector_alpha_3)
      allocate(NEFR4Out%Size_vector_alpha_3(1))
      if(allocated(NEFR4Out%Size_vector_alpha_4)) deallocate(NEFR4Out%Size_vector_alpha_4)
      allocate(NEFR4Out%Size_vector_alpha_4(1))
      if(allocated(NEFR4Out%Size_vector_beta_3)) deallocate(NEFR4Out%Size_vector_beta_3)
      allocate(NEFR4Out%Size_vector_beta_3(1))
      if(allocated(NEFR4Out%Size_vector_beta_4)) deallocate(NEFR4Out%Size_vector_beta_4)
      allocate(NEFR4Out%Size_vector_beta_4(1))
      NEFR4Out%Size_vector_alpha_3(1) = size(ERIB,'alpha',kb)
      NEFR4Out%Size_vector_alpha_4(1) = size(ERIB,'alpha',lb)
      NEFR4Out%Size_vector_beta_3(1) = size(ERIB,'beta',kb)
      NEFR4Out%Size_vector_beta_4(1) = size(ERIB,'beta',lb)

      do id1 = 1, nefr4Out%blockSize(1)
        do id2 = 1, nefr4Out%blockSize(2)
          first=.true.
          out_index = NEFR4Out%blockIndex([1,id1],[2,id2],[3,1],[4,1])
          A_index = NEFR4A%blockIndex([ka,id1],[la,id2],[ia,1],[ja,1])
          if(first) then
            NEFR4Out%Main_R4(out_index) =  &
              mqc_eri_eri_partial_contraction_2eri(&
              NEFR4A%Main_R4(A_index),ERIB,my_index1,my_index2)
            first=.false.
          else
            NEFR4Out%Main_R4(out_index) = NEFR4Out%Main_R4(out_index)+&
              mqc_eri_eri_partial_contraction_2eri(&
              NEFR4A%Main_R4(A_index),ERIB,my_index1,my_index2)
          end if
        end do
      end do
!
      end function mqc_nef_r4_eri_partial_contraction_2nef_r4
!
!
!     PROCEDURE MQC_ERI_NEF_R4_Partial_Contraction_2R4
!
!>    \brief <b> MQC_ERI_NEF_R4_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC 2ERI </b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_ERI_NEF_R4_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC ERI
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] ERIA
!>    \verbatim
!>        ERIA is type(mqc_twoeris) 
!>        The second MQC 2ERI to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] Index1
!>    \verbatim
!>        Index1 is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIs1 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>        Ex: If Index1=[3,2,4,1], the third and second indexes of ERIs1 are contracted.
!>            The fourth and first are multiplied to form the first and second indices
!>            of the result matrix, respectivelly.
!>    \endverbatim
!>
!>    \param[in] Index2
!>    \verbatim
!>        Index2 is Integer(kind=int64),dimension(4),optional
!>        Contain the indexes of ERIs2 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>
!>        Ex: index1=[3,2,4,1] and index2=[2,1,4,3]
!>        sum_{i,j} R4A_{pjiq} ERIB_{jisr} = R4Out_{pqrs}
!>
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_eri_nef_r4_partial_contraction_2nef_r4(ERIA,NEFR4B,index1,index2) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4B
      type(mqc_twoeris),intent(in)::ERIA
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(mqc_scalar)::val
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
      logical::first
!
      if(present(index1)) then
        ia=index1(1)
        ja=index1(2)
        ka=index1(3)
        la=index1(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ka', ka,'la',la)
        my_index1 = index1
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_index1 = [1,2,3,4]
      endif
      if(present(index2)) then
        ib=index2(1)
        jb=index2(2)
        kb=index2(3)
        lb=index2(4)
        if (ib.eq.jb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'jb',jb)
        if (ib.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'kb',kb)
        if (ib.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'lb',lb)
        if (jb.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'jb', jb,'kb',kb)
        if (jb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'jb', jb,'lb',lb)
        if (kb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'kb', kb,'lb',lb)
        my_index2 = index2
      else
        ib=1
        jb=2
        kb=3
        lb=4
        my_index2 = [1,2,3,4]
      endif
      if (ia.gt.4.or.ia.le.0.or.ib.gt.4.or.ib.le.0.or.&
          ja.gt.4.or.ja.le.0.or.jb.gt.4.or.jb.le.0.or.&
          ka.gt.4.or.ka.le.0.or.kb.gt.4.or.kb.le.0.or.&
          la.gt.4.or.la.le.0.or.lb.gt.4.or.lb.le.0.) Call MQC_Error(&
           'Index out of range in MQC_ERI_NEF_R4_Partial_Contraction_2R4', 6)

      if (nefr4B%blockSize(ib).ne.1.or.nefr4B%blockSize(jb).ne.1) &
          call mqc_error_i('NEFR4 not correctly dimensioned in mqc_eri_nef_r4_contraction_2R4',6,&
           'dim1',nefr4B%blockSize(ib),'dim2',nefr4B%blockSize(jb))
      if(NEFR4B%blockSize(ib,1,'alpha').ne.size(ERIA,'alpha',ia).or. &
         NEFR4B%blockSize(ib,1,'beta').ne.size(ERIA,'beta',ia)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_ERI_NEF_R4_Partial_Contraction_2R4',6,&
        'NEFR4B%blockSize(ib,1,alpha)',NEFR4B%blockSize(ib,1,'alpha'),&
        'size(ERIA,alpha,ia)',size(ERIA,'alpha',ia),&
        'NEFR4B%blockSize(ib,1,beta)',NEFR4B%blockSize(ib,1,'beta'),&
        'size(ERIA,beta,ia)',size(ERIA,'beta',ia))
      if(NEFR4B%blockSize(jb,1,'alpha').ne.size(ERIA,'alpha',ja).or. &
         NEFR4B%blockSize(jb,1,'beta').ne.size(ERIA,'beta',ja)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_ERI_NEF_R4_Partial_Contraction_2R4',6,&
        'NEFR4B%blockSize(jb,1,alpha)',NEFR4B%blockSize(jb,1,'alpha'),&
        'size(ERIA,alpha,ja)',size(ERIA,'alpha',ja),&
        'NEFR4B%blockSize(jb,1,beta)',NEFR4B%blockSize(jb,1,'beta'),&
        'size(ERIA,beta,ja)',size(ERIA,'beta',ja))

      NEFR4Out%NumSpace_1 = 1
      NEFR4Out%NumSpace_2 = 1
      NEFR4Out%NumSpace_3 = NEFR4B%blockSize(kb) 
      NEFR4Out%NumSpace_4 = NEFR4B%blockSize(lb)
      N1234 = NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4

      allocate(NEFR4Out%Main_R4(N1234))
      NEFR4Out%Map_Matrix_3 = NEFR4B%blockMap(kb) 
      NEFR4Out%Map_Matrix_4 = NEFR4B%blockMap(lb) 
      if(allocated(NEFR4Out%Map_Matrix_1)) deallocate(NEFR4Out%Map_Matrix_1)
      allocate(NEFR4Out%Map_Matrix_1(1))
      if(allocated(NEFR4Out%Map_Matrix_2)) deallocate(NEFR4Out%Map_Matrix_2)
      allocate(NEFR4Out%Map_Matrix_2(1))
      call NEFR4Out%Map_Matrix_1(1)%init(size(ERIA,'alpha',ka),scalar=1,lengthBeta=size(ERIA,'beta',ka))
      call NEFR4Out%Map_Matrix_2(1)%init(size(ERIA,'alpha',la),scalar=1,lengthBeta=size(ERIA,'beta',la))
      do id1 = 1, size(ERIA,'alpha',ka)
        val = id1
        call  NEFR4Out%Map_Matrix_1(1)%put(val,id1,'alpha')
      end do
      do id1 = 1, size(ERIA,'beta',ka)
        val = id1
        call  NEFR4Out%Map_Matrix_1(1)%put(val,id1,'beta')
      end do
      do id1 = 1, size(ERIA,'alpha',la)
        val = id1
        call  NEFR4Out%Map_Matrix_2(1)%put(val,id1,'alpha')
      end do
      do id1 = 1, size(ERIA,'beta',la)
        val = id1
        call  NEFR4Out%Map_Matrix_2(1)%put(val,id1,'beta')
      end do
    
      NEFR4Out%Size_vector_alpha_3 = NEFR4B%blockSizeVec(kb,'alpha') 
      NEFR4Out%Size_vector_alpha_4 = NEFR4B%blockSizeVec(lb,'alpha') 
      NEFR4Out%Size_vector_beta_3 = NEFR4B%blockSizeVec(kb,'beta')  
      NEFR4Out%Size_vector_beta_4 = NEFR4B%blockSizeVec(lb,'beta')  
      if(allocated(NEFR4Out%Size_vector_alpha_1)) deallocate(NEFR4Out%Size_vector_alpha_1)
      allocate(NEFR4Out%Size_vector_alpha_1(1))
      if(allocated(NEFR4Out%Size_vector_alpha_2)) deallocate(NEFR4Out%Size_vector_alpha_2)
      allocate(NEFR4Out%Size_vector_alpha_2(1))
      if(allocated(NEFR4Out%Size_vector_beta_1)) deallocate(NEFR4Out%Size_vector_beta_1)
      allocate(NEFR4Out%Size_vector_beta_1(1))
      if(allocated(NEFR4Out%Size_vector_beta_2)) deallocate(NEFR4Out%Size_vector_beta_2)
      allocate(NEFR4Out%Size_vector_beta_2(1))
      NEFR4Out%Size_vector_alpha_1(1) = size(ERIA,'alpha',ka)
      NEFR4Out%Size_vector_alpha_2(1) = size(ERIA,'alpha',la)
      NEFR4Out%Size_vector_beta_1(1) = size(ERIA,'beta',ka)
      NEFR4Out%Size_vector_beta_2(1) = size(ERIA,'beta',la)

      do id1 = 1, nefr4Out%blockSize(3)
        do id2 = 1, nefr4Out%blockSize(4)
          out_index = NEFR4Out%blockIndex([1,1],[2,1],[3,id1],[4,id2])
          first=.true.
          A_index = NEFR4B%blockIndex([kb,id1],[lb,id2],[ib,1],[jb,1])
          if(first) then
            NEFR4Out%Main_R4(out_index) =  &
              mqc_eri_eri_partial_contraction_2eri(&
              ERIA,NEFR4B%Main_R4(A_index),my_index1,my_index2)
            first=.false.
          else
            NEFR4Out%Main_R4(out_index) = NEFR4Out%Main_R4(out_index)+&
              mqc_eri_eri_partial_contraction_2eri(&
              ERIA,NEFR4B%Main_R4(A_index),my_index1,my_index2)
          end if
        end do
      end do
!
      end function mqc_eri_nef_r4_partial_contraction_2nef_r4
!
!
!     PROCEDURE MQC_NEF_R4_NEF_R4_Partial_Contraction
!
!>    \brief <b> MQC_NEF_R4_NEF_R4_Contraction is used to return the contraction of 
!>    MQC NEF R4 with a MQC NEF R4 </b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_NEF_R4_Contraction is used to return the contraction of 
!>    MQC NEF R4 with a MQC NEF R4 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4) 
!>        The second MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] Index1
!>    \verbatim
!>        Index1 is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIs1 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>        Ex: If Index1=[3,2,4,1], the third and second indexes of ERIs1 are contracted.
!>            The fourth and first are multiplied to form the first and second indices
!>            of the result matrix, respectivelly.
!>    \endverbatim
!>
!>    \param[in] Index2
!>    \verbatim
!>        Index2 is Integer(kind=int64),dimension(4),optional
!>        Contain the indexes of ERIs2 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>
!>        Ex: index1=[3,2,4,1] and index2=[2,1,4,3]
!>        sum_{i,j} R4A_{ljik} R4B_{jilk} = MatOut_{kl}
!>
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_nef_r4_partial_contraction(NEFR4A,NEFR4B,index1,index2) result(MatOut)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_matrix)::MatOut
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4,mat_index
!
      if(present(index1)) then
        ia=index1(1)
        ja=index1(2)
        ka=index1(3)
        la=index1(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ka', ka,'la',la)
        my_index1 = index1
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_index1 = [1,2,3,4]
      endif
      if(present(index2)) then
        ib=index2(1)
        jb=index2(2)
        kb=index2(3)
        lb=index2(4)
        if (ib.eq.jb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ib', ib,'jb',jb)
        if (ib.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ib', ib,'kb',kb)
        if (ib.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ib', ib,'lb',lb)
        if (jb.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'jb', jb,'kb',kb)
        if (jb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'jb', jb,'lb',lb)
        if (kb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'kb', kb,'lb',lb)
        my_index2 = index2
      else
        ib=1
        jb=2
        kb=3
        lb=4
        my_index2 = [1,2,3,4]
      endif
      if (ia.gt.4.or.ia.le.0.or.ib.gt.4.or.ib.le.0.or.&
          ja.gt.4.or.ja.le.0.or.jb.gt.4.or.jb.le.0.or.&
          ka.gt.4.or.ka.le.0.or.kb.gt.4.or.kb.le.0.or.&
          la.gt.4.or.la.le.0.or.lb.gt.4.or.lb.le.0.) Call MQC_Error(&
           'Index out of range in MQC_NEF_R4_NEF_R4_Partial_Contraction', 6)

      if (nefr4A%blockSize(ia).ne.nefr4B%blockSize(ib).or.&
          nefr4A%blockSize(ja).ne.nefr4B%blockSize(jb).or.&
          nefr4A%blockSize(ka).ne.nefr4B%blockSize(kb).or.&
          nefr4A%blockSize(la).ne.nefr4B%blockSize(lb))&
          call mqc_error_i('NEFR4 not correctly dimensioned in mqc_nef_r4_nef_r4_contraction',6,&
           'dim1',nefr4A%blockSize(ia)-nefr4B%blockSize(ib),&
           'dim2',nefr4A%blockSize(ja)-nefr4B%blockSize(jb),&
           'dim3',nefr4A%blockSize(ka)-nefr4B%blockSize(kb),&
           'dim4',nefr4A%blockSize(la)-nefr4B%blockSize(lb))

      Map_Matrix_Col = NEFR4A%blockMap(la)
      Map_Matrix_Row = NEFR4A%blockMap(ka)
      call MatOut%init(Map_Matrix_Row,initial_map_matrix_col=Map_Matrix_Col)
      do id3 = 1, nefr4A%blockSize(ka)
        do id4 = 1, nefr4A%blockSize(la)
          mat_index = MatOut%blockIndex([1,id3],[2,id4])
          do id1 = 1, nefr4A%blockSize(ia)
            do id2 = 1, nefr4A%blockSize(ja)
              MatOut%Main_Matrix(mat_index) = MatOut%Main_Matrix(mat_index) +&
                partialcontraction(NEFR4A%Main_R4(NEFR4A%blockIndex([ka,id3],[la,id4],[ia,id1],[ja,id2])),&
                NEFR4B%Main_R4(NEFR4B%blockIndex([kb,id3],[lb,id4],[ib,id1],[jb,id2])),my_index1,my_index2)
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_nef_r4_partial_contraction
!
!
!     PROCEDURE MQC_NEF_R4_Dimension
!
!>    \brief <b> MQC_NEF_R4_Dimension is a function that returns the length of a spin
!>    block axis of an MQC NEF R4 object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Dimension is a function that returns the length of a spin block axis 
!>    of an MQC NEF R4 object. The required argument label specifies which block is queried. 
!>    The optional arguments axis, subaxis and spin specifies which dimension is returned. 
!>    The following options are available:
!>
!>    1. Axis = '1' selects the first dimension (default).
!>    2. Axis = '2' selects the second dimension.
!>    3. Axis = '3' selects the third dimension.
!>    4. Axis = '4' selects the fourth dimension.
!>
!>    1. .not.present(SubAxis) returns the dimension of sub-spaces in <Axis> (default).
!>    2. SubAxis = 'n' selects the n-th dimension of sub-spaces in <Axis>.
!>
!>    1. Spin = returns the 'alpha' dimension of the n-th dimension of  sub-spaces in <Axis> (default).
!>    2. Spin = returns the 'beta' dimension of the n-th dimension of  sub-spaces in <Axis>.
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The name of the MQC NEF R4 variable.
!>    \endverbatim
!>
!>    \param[in] Axis
!>    \verbatim
!>        Axis is Integer(kind=int64),Optional
!>        The spin block axis to return.
!>        = 1: number of first dimension (default).
!>        = 2: number of second dimension.
!>        = 3: number of third dimension.
!>        = 4: number of fourth dimension.
!>    \endverbatim
!>
!>    \param[in] SubAxis
!>    \verbatim
!>        SubAxis is Integer(kind=int64),Optional
!>       .not.present(SubAxis) returns the dimension of sub-spaces in <Axis> (default).
!>       SubAxis = 'n' selects the n-th dimension of sub-spaces in <Axis>.
!>    \endverbatim
!>
!>    \param[in] 
!>    \verbatim Spin
!>        Spin is Character(Len=*),optional
!>        The spin block to query.
!>        = 'alpha','a': alpha axis dimension.
!>        = 'beta','b':  beta axis dimension.
!>        Returns an error if present while subaxis is not present
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_dimension(nefr4,axis,subaxis,spin) result(dimBlock)
!
      implicit none
      class(mqc_nef_r4),intent(in)::nefr4
      integer(kind=int64)::dimBlock
      Character(Len=*),intent(in),optional::spin
      integer(kind=int64),intent(in),optional::axis,subaxis
      Character(Len=64)::my_spin
      integer(kind=int64)::my_axis,my_sub_axis
      logical::do_sub
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf
      if(my_axis.le.0.or.my_axis.gt.4) &
        call mqc_error_I('Unknown axis specified in MQC_NEF_R4_Dimension',6,'my_axis',my_axis)

      if(present(subaxis)) then
        my_sub_axis = subaxis
        do_sub=.true.
      else
        do_sub=.false.
      endIf
      if(present(spin)) then
        if (.not.do_sub) &
           call mqc_error('Spin specified without subaxis in MQC_NEF_R4_Dimension',6)
        call string_change_case(spin,'l',my_spin)
      end if
     
      if (do_sub) then
        select case(my_spin)
        case('alpha','a')
          if(my_axis.eq.1) then
            if (my_sub_axis.gt.NEFR4%NumSpace_1) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_1',NEFR4%NumSpace_1,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_alpha_1(my_sub_axis)
          else if(my_axis.eq.2) then
            if (my_sub_axis.gt.NEFR4%NumSpace_2) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_2',NEFR4%NumSpace_2,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_alpha_2(my_sub_axis)
          else if(my_axis.eq.3) then
            if (my_sub_axis.gt.NEFR4%NumSpace_3) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_3',NEFR4%NumSpace_3,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_alpha_3(my_sub_axis)
          else if(my_axis.eq.4) then
            if (my_sub_axis.gt.NEFR4%NumSpace_4) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_4',NEFR4%NumSpace_4,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_alpha_4(my_sub_axis)
          end if
        case('beta','b')
          if(my_axis.eq.1) then
            if (my_sub_axis.gt.NEFR4%NumSpace_1) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_1',NEFR4%NumSpace_1,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_beta_1(my_sub_axis)
          else if(my_axis.eq.2) then
            if (my_sub_axis.gt.NEFR4%NumSpace_2) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_2',NEFR4%NumSpace_2,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_beta_2(my_sub_axis)
          else if(my_axis.eq.3) then
            if (my_sub_axis.gt.NEFR4%NumSpace_3) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_3',NEFR4%NumSpace_3,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_beta_3(my_sub_axis)
          else if(my_axis.eq.4) then
            if (my_sub_axis.gt.NEFR4%NumSpace_4) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_4',NEFR4%NumSpace_4,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_beta_4(my_sub_axis)
          end if
        case default
          call mqc_error_A('Spin label not valid in mqc_nef_r4_dimension', 6, &
               'my_spin', my_spin )
        end select
      else
        if(my_axis.eq.1) then
          dimBlock = NEFR4%NumSpace_1
        else if(my_axis.eq.2) then
          dimBlock = NEFR4%NumSpace_2
        else if(my_axis.eq.3) then
          dimBlock = NEFR4%NumSpace_3
        else if(my_axis.eq.4) then
          dimBlock = NEFR4%NumSpace_4
        end if
      end if
!
      end function mqc_nef_r4_dimension
!
!
!     PROCEDURE MQC_NEF_R4_Dimension_vector
!
!>    \brief <b> MQC_NEF_R4_Dimension is a function that returns the size_vector_spin of a
!>    block axis of an MQC NEF R4 object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Dimension is a function that returns the size_vector_spin block axis 
!>    of an MQC NEF R4 object. The required argument label specifies which block is queried. 
!>    The optional arguments axis, subaxis and spin specifies which dimension is returned. 
!>    The following options are available:
!>
!>    1. Axis = '1' selects the first dimension (default).
!>    2. Axis = '2' selects the second dimension.
!>    3. Axis = '3' selects the third dimension.
!>    4. Axis = '4' selects the fourth dimension.
!>
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The name of the MQC NEF R4 variable.
!>    \endverbatim
!>
!>    \param[in] Axis
!>    \verbatim
!>        Axis is Integer(kind=int64),Optional
!>        The spin block axis to return.
!>        = 1: number of first dimension (default).
!>        = 2: number of second dimension.
!>        = 3: number of third dimension.
!>        = 4: number of fourth dimension.
!>    \endverbatim
!>
!>    \param[in] 
!>    \verbatim Spin
!>        Spin is Character(Len=*),optional
!>        The spin block to query.
!>        = 'alpha': alpha axis dimension.
!>        = 'beta':  beta axis dimension.
!>        Returns an error if present while subaxis is not present
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_dimension_vector(nefr4,axis,spin) result(dimVec)
!
      implicit none
      class(mqc_nef_r4),intent(in)::nefr4
      integer,dimension(:),allocatable::dimVec
      Character(Len=*),intent(in),optional::spin
      integer(kind=int64),intent(in),optional::axis
      Character(Len=64)::my_spin
      integer(kind=int64)::my_axis,my_sub_axis
      logical::do_sub=.false.
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf
      if(my_axis.le.0.or.my_axis.gt.4) &
        call mqc_error_I('Unknown axis specified in MQC_NEF_R4_Dimension',6,'my_axis',my_axis)
      if(present(spin)) then
        my_spin = spin
      else
        my_spin = 'alpha'
      endIf

      select case(my_spin)
      case('alpha')
        if(my_axis.eq.1) dimVec = NEFR4%Size_vector_alpha_1
        if(my_axis.eq.2) dimVec = NEFR4%Size_vector_alpha_2
        if(my_axis.eq.3) dimVec = NEFR4%Size_vector_alpha_3
        if(my_axis.eq.4) dimVec = NEFR4%Size_vector_alpha_4
      case('beta')
        if(my_axis.eq.1) dimVec = NEFR4%Size_vector_beta_1
        if(my_axis.eq.2) dimVec = NEFR4%Size_vector_beta_2
        if(my_axis.eq.3) dimVec = NEFR4%Size_vector_beta_3
        if(my_axis.eq.4) dimVec = NEFR4%Size_vector_beta_4
      case default
        call mqc_error_A('Unknown spin specified in MQC_NEF_R4_Dimension',6,'my_spin',my_spin)
      end select
!
      end function mqc_nef_r4_dimension_vector
!
!
!     PROCEDURE MQC_NEF_R4_Map
!
!>    \brief <b> MQC_NEF_R4_Map is a function that returns the Map_Matrix
!>    of an MQC NEF R4 object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Map is a function that returns the Map_Matrix
!>    of an MQC NEF R4 object
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The name of the MQC NEF R4 variable.
!>    \endverbatim
!>
!>    \param[in] Axis
!>    \verbatim
!>        Axis is Integer(kind=int64),Optional
!>        The spin block axis to return.
!>        = 1: number of first dimension (default).
!>        = 2: number of second dimension.
!>        = 3: number of third dimension.
!>        = 4: number of fourth dimension.
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_map(nefr4,axis) result(MapOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::nefr4
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::MapOut
      integer(kind=int64),intent(in),optional::axis
      integer(kind=int64)::my_axis
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf
      if(my_axis.le.0.or.my_axis.gt.4) &
        call mqc_error_I('Unknown axis specified in MQC_NEF_R4_Dimension',6,'my_axis',my_axis)

      if(my_axis.eq.1) then
        MapOut = NEFR4%Map_Matrix_1
      else if(my_axis.eq.2) then
        MapOut = NEFR4%Map_Matrix_2
      else if(my_axis.eq.3) then
        MapOut = NEFR4%Map_Matrix_3
      else if(my_axis.eq.4) then
        MapOut = NEFR4%Map_Matrix_4
      end if
!
      end function mqc_nef_r4_map
!
!
!     PROCEDURE MQC_NEF_R4_Index
!
!>    \brief <b> MQC_NEF_R4_Index is a function that returns the main index
!>    of an MQC NEF R4 object based on sub-indexes</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Index is a function that returns the main index
!>    of an MQC NEF R4 object based on sub-indexes. This function is used when
!>    the main index is required, but the order of the axis is variable.
!>
!>    Example: given the four sets:
!>      set1 = [1,4]; set2 = [3,1]; set3 = [4,10] and set4 = [2,3]
!>      This function returns the correspondent combined-index for ERI(4,3,1,10)
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The name of the MQC NEF R4 variable.
!>    \endverbatim
!>
!>    \param[in] set1
!>    \verbatim
!>        set1 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!>
!>    \param[in] set2
!>    \verbatim
!>        set2 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!>
!>    \param[in] set3
!>    \verbatim
!>        set3 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!>
!>    \param[in] set4
!>    \verbatim
!>        set4 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_index(nefr4,set1,set2,set3,set4) result(indexOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::nefr4
      integer(kind=int64),dimension(2),intent(in)::set1,set2,set3,set4
      integer(kind=int64)::indexOut
      integer(kind=int64)::id1,id2,id3,id4
      integer(kind=int64)::N1,N12,N123,N1234
!
      if (set1(1).le.0.or.set1(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,&
        'set1(1)',set1(1))
      if (set2(1).le.0.or.set2(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,&
        'set2(1)',set2(1))
      if (set3(1).le.0.or.set3(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,&
        'set3(1)',set3(1))
      if (set4(1).le.0.or.set4(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,&
        'set4(1)',set4(1))
    
      id1=-1
      id2=-1
      id3=-1
      id4=-1

      if (set1(1).eq.1) id1 = set1(2)
      if (set2(1).eq.1) then
        if (id1.gt.0) call mqc_error(&
          'Axis 1 index defined more than once in mqc_nef_r4_index',6)
        id1 = set2(2)
      endIf
      if (set3(1).eq.1) then
        if (id1.gt.0) call mqc_error(&
          'Axis 1 index defined more than once in mqc_nef_r4_index',6)
        id1 = set3(2)
      endIf
      if (set4(1).eq.1) then
        if (id1.gt.0) call mqc_error(&
          'Axis 1 index defined more than once in mqc_nef_r4_index',6)
        id1 = set4(2)
      endIf
      if (set1(1).eq.2) id2 = set1(2)
      if (set2(1).eq.2) then
        if (id2.gt.0) call mqc_error(&
          'Axis 2 index defined more than once in mqc_nef_r4_index',6)
        id2 = set2(2)
      endIf
      if (set3(1).eq.2) then
        if (id2.gt.0) call mqc_error(&
          'Axis 2 index defined more than once in mqc_nef_r4_index',6)
        id2 = set3(2)
      endIf
      if (set4(1).eq.2) then
        if (id2.gt.0) call mqc_error(&
          'Axis 2 index defined more than once in mqc_nef_r4_index',6)
        id2 = set4(2)
      endIf
      if (set1(1).eq.3) id3 = set1(2)
      if (set2(1).eq.3) then
        if (id3.gt.0) call mqc_error(&
          'Axis 3 index defined more than once in mqc_nef_r4_index',6)
        id3 = set2(2)
      endIf
      if (set3(1).eq.3) then
        if (id3.gt.0) call mqc_error(&
          'Axis 3 index defined more than once in mqc_nef_r4_index',6)
        id3 = set3(2)
      endIf
      if (set4(1).eq.3) then
        if (id3.gt.0) call mqc_error(&
          'Axis 3 index defined more than once in mqc_nef_r4_index',6)
        id3 = set4(2)
      endIf
      if (set1(1).eq.4) id4 = set1(2)
      if (set2(1).eq.4) then
        if (id4.gt.0) call mqc_error(&
          'Axis 4 index defined more than once in mqc_nef_r4_index',6)
        id4 = set2(2)
      endIf
      if (set3(1).eq.4) then
        if (id4.gt.0) call mqc_error(&
          'Axis 4 index defined more than once in mqc_nef_r4_index',6)
        id4 = set3(2)
      endIf
      if (set4(1).eq.4) then
        if (id4.gt.0) call mqc_error(&
          'Axis 4 index defined more than once in mqc_nef_r4_index',6)
        id4 = set4(2)
      endIf
      N1 = NEFR4%NumSpace_4
      N12 = N1*NEFR4%NumSpace_3
      N123 = N12*NEFR4%NumSpace_2
      N1234 = N123*NEFR4%NumSpace_1
      indexOut = (id1-1)*N123+(id2-1)*N12+(id3-1)*N1+id4
      !write(*,*) 'Sets = ', set1,set2,set3,set4
      !write(*,*) 'indexIn = ', id1,id2,id3,id4
      !write(*,*) 'indexOut = ', indexOut
!
      end function mqc_nef_r4_index
!
!
!     PROCEDURE MQC_NEF_Matrix_Index
!
!>    \brief <b> MQC_NEF_Matrix_Index is a function that returns the main index
!>    of an MQC NEF Matrix object based on sub-indexes</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Index is a function that returns the main index
!>    of an MQC NEF Matrix object based on sub-indexes. This function is used when
!>    the main index is required, but the order of the axis is variable.
!>
!>    Example: given the four sets:
!>      set1 = [2,4]; set2 = [1,2]
!>      This function returns the correspondent combined-index for ERI(2,4)
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The name of the MQC NEF Matrix variable.
!>    \endverbatim
!>
!>    \param[in] set1
!>    \verbatim
!>        set1 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!>
!>    \param[in] set2
!>    \verbatim
!>        set2 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_index(nefmat,set1,set2) result(indexOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::nefmat
      integer(kind=int64),dimension(2),intent(in)::set1,set2
      integer(kind=int64)::indexOut
      integer(kind=int64)::id1,id2
!
      if (set1(1).le.0.or.set1(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nefmatrix_index',6,&
        'set1(1)',set1(1))
      if (set2(1).le.0.or.set2(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nefmatrix_index',6,&
        'set2(1)',set2(1))
    
      id1=-1
      id2=-1

      if (set1(1).eq.1) id1 = set1(2)
      if (set2(1).eq.1) then
        if (id1.gt.0) call mqc_error(&
          'Axis 1 index defined more than once in mqc_nefmatrix_index',6)
        id1 = set2(2)
      endIf
      if (set1(1).eq.2) id2 = set1(2)
      if (set2(1).eq.2) then
        if (id2.gt.0) call mqc_error(&
          'Axis 2 index defined more than once in mqc_nefmatrix_index',6)
        id2 = set2(2)
      endIf
      indexOut = (id1-1)*NEFMat%NumSpace_col+id2
      !write(*,*) 'Sets = ', set1,set2
      !write(*,*) 'indexIn = ', id1,id2
      !write(*,*) 'indexOut = ', indexOut
!
      end function mqc_nef_matrix_index
!
!
!
!     PROCEDURE MQC_NEF_Matrix_At
!
!>    \brief <b> MQC_NEF_Matrix_At is a function that returns the value of an element 
!>    of a MQC NEF Matrix variable</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_Matrix_At is a function that returns the value of (I,J)th element of a 
!>    (p,q) block of a MQC NEF Matrix variable as an MQC scalar.  
!>    Optional argument SpinBlock specifies the spin block of MQC NEF Matrix.
!>    are available:
!>
!>    1. SpinBlock = 'alpha' uses the alpha spin block.
!>    2. SpinBlock = 'beta' uses the beta spin block.
!>    3. SpinBlock = 'alphaBeta' uses the alpha-beta spin block.
!>    4. SpinBlock = 'betaAlpha' uses the beta-alpha spin block.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] nefM
!>    \verbatim
!>        nefM is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to return the value of the (I,J)th for block (p,q)
!>        element.
!>    \endverbatim
!>    
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The first dimension of the element in MQC NEF Matrix
!>          If I>0 row count is from first index
!>    \endverbatim
!>    
!>    \param[in] J  
!>    \verbatim
!>        J is Integer(kind=int64)
!>        The second dimension of the element in MQC NEF Matrix
!>          If J>0 col count is from first index
!>    \endverbatim
!>    
!>    \param[in] p 
!>    \verbatim
!>        p is Integer(kind=int64)
!>        The first excitation dimension of the element in MQC NEF Matrix
!>          If p>0 row count is from first index
!>    \endverbatim
!>    
!>    \param[in] q 
!>    \verbatim
!>        q is Integer(kind=int64)
!>        The second excitation dimension of the element in MQC NEF Matrix
!>          If q>0 row count is from first index
!>    \endverbatim
!>
!>    \param[in] spinBlockIn
!>    \verbatim
!>        spinBlockIn is Character(Len=*),Optional 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_at(nefM,i,j,p,q,spinBlockIn) result(element)
!
      implicit none
      class(mqc_nef_matrix)::nefM
      class(*),intent(in)::i,j,p,q
      integer(kind=int64)::myI,myJ,myp,myq
      type(mqc_scalar)::element
      character(len=*),optional,intent(in)::spinBlockIn
      character(len=64)::spinBlock
      real(kind=real64)::zero=0.0d0
!
      element = zero
      select type (i)
      type is (integer)
        myI = i
      type is (mqc_scalar)
        myI = i
      class default
        call mqc_error('index type not defined in MQC_NEF_Matrix_at')
      end select
      select type (j)
      type is (integer)
        myJ = j
      type is (mqc_scalar)
        myJ = j
      class default
        call mqc_error('index type not defined in MQC_NEF_Matrix_at')
      end select
      select type (p)
      type is (integer)
        myp = p
      type is (mqc_scalar)
        myp = p
      class default
        call mqc_error('index type not defined in MQC_NEF_Matrix_at')
      end select
      select type (q)
      type is (integer)
        myq = q
      type is (mqc_scalar)
        myq = q
      class default
        call mqc_error('index type not defined in MQC_NEF_Matrix_at')
      end select
!
      if(myp.gt.nefM%NumSpace_row) &
        call mqc_error('Requested p dimension out of bounds in mqc_NEF_Matrix_at')
      if(myq.gt.nefM%NumSpace_col) &
        call mqc_error('Requested q dimension out of bounds in mqc_NEF_Matrix_at')
      if(myI.gt.(nefM%Size_vector_alpha_row(myp)+nefM%Size_vector_beta_row(myp))) &
        call mqc_error('Requested I dimension out of bounds in mqc_NEF_Matrix_at')
      if(myJ.gt.(nefM%Size_vector_alpha_col(myq)+nefM%Size_vector_beta_col(myq))) &
        call mqc_error('Requested J dimension out of bounds in mqc_NEF_Matrix_at')
!  
      if(present(spinBlockIn)) then
        call string_change_case(spinBlockIn,'L',spinBlock)
      else
        spinBlock = 'alpha'
      end if
      element = nefM%Main_Matrix(nefM%blockIndex([1,myp],[2,myq]))%at(myI,myJ,spinBlock)
!
      end function mqc_nef_matrix_at
!
!
!
!     PROCEDURE MQC_NEF_R4_At
!
!>    \brief <b> MQC_NEF_R4_At is a function that returns the value of an element 
!>    of a MQC NEF R4 variable</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_R4_At is a function that returns the value of (I,J,K,L)th element of a 
!>    (p,q,r,s) block of a MQC NEF R4 variable as an MQC scalar.  
!>    Optional argument SpinBlock specifies the spin block of MQC NEF R4.
!>    Each index is represented by a set = [dI,I,p,s], where s is 1 for alpha or -1 for beta
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] nefR4
!>    \verbatim
!>        nefR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to return the value 
!>        element.
!>    \endverbatim
!>    
!>    \param[in] Set1  
!>    \verbatim
!>        set1 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>    
!>    \param[in] Set2  
!>    \verbatim
!>        set2 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>
!>    \param[in] Set3  
!>    \verbatim
!>        set3 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>
!>    \param[in] Set4  
!>    \verbatim
!>        set4 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_at(nefR4,set1,set2,set3,set4) result(element)
!
      implicit none
      class(mqc_nef_r4)::nefR4
      Integer(kind=int64),dimension(4),intent(in)::set1,set2,set3,set4
      integer(kind=int64)::I,J,K,L,myp,myq,myr,mys,sp,sq,sr,ss,id1,id2,id3,id4
      integer(kind=int64)::dI,dJ,dK,dL,r,s,p,q
      type(mqc_scalar)::element
      character(len=1)::sI,sJ,sK,sL
      character(len=4)::spinBlock
      real(kind=real64)::zero=0.0d0
!
      element = zero
!
      dI = set1(1)
      dJ = set2(1)
      dK = set3(1)
      dL = set4(1)
      I = set1(2)
      J = set2(2)
      K = set3(2)
      L = set4(2)
      p = set1(3)
      q = set2(3)
      r = set3(3)
      s = set4(3)
      sp = set1(4)
      sq = set2(4)
      sr = set3(4)
      ss = set4(4)
!
      if (dI.le.0.or.dI.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,'dI',dI)
      if (dJ.le.0.or.dJ.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,'dJ',dJ)
      if (dK.le.0.or.dK.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,'dK',dK)
      if (dL.le.0.or.dL.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,'dL',dL)
!
      if (dI.eq.dJ) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dI',dI,'dJ',dJ)
      if (dI.eq.dK) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dI',dI,'dK',dK)
      if (dI.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dI',dI,'dL',dL)
      if (dJ.eq.dK) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dJ',dJ,'dK',dK)
      if (dJ.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dJ',dJ,'dL',dL)
      if (dK.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dK',dK,'dL',dL)
!
      if(I.gt.nefR4%blockSize(dI)) &
        call mqc_error('Requested I dimension out of bounds in mqc_NEF_R4_at')
      if(J.gt.nefR4%blockSize(dJ)) &
        call mqc_error('Requested J dimension out of bounds in mqc_NEF_R4_at')
      if(K.gt.nefR4%blockSize(dK)) &
        call mqc_error('Requested K dimension out of bounds in mqc_NEF_R4_at')
      if(L.gt.nefR4%blockSize(dL)) &
        call mqc_error('Requested L dimension out of bounds in mqc_NEF_R4_at')
!
      if (sp.eq.-1) then
        sI = 'b'
      else if (sp.eq.1) then
        sI = 'a'
      else
        call mqc_error_i('Requested I spin not recognized in mqc_NEF_R4_at',6,&
                       'sp',sp)
      end if
      if (sq.eq.-1) then
        sJ = 'b'
      else if (sq.eq.1) then
        sJ = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_at',6,&
                       'sq',sq)
      end if
      if (sr.eq.-1) then
        sK = 'b'
      else if (sr.eq.1) then
        sK = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_at',6,&
                       'sr',sr)
      end if
      if (ss.eq.-1) then
        sL = 'b'
      else if (ss.eq.1) then
        sL = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_at',6,&
                       'ss',ss)
      end if
!
      spinBlock=''
      if (dI.eq.1) then
        id1 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.1) then
        id1 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.1) then
        id1 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else 
        id1 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.2) then
        id2 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.2) then
        id2 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.2) then
        id2 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else 
        id2 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.3) then
        id3 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.3) then
        id3 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.3) then
        id3 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else 
        id3 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.4) then
        id4 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.4) then
        id4 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.4) then
        id4 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else 
        id4 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
!  
      element = nefR4%Main_R4(nefR4%blockIndex(&
                          [dI,I],[dJ,J],[dK,K],[dL,L]))%at(&
                           id1,id2,id3,id4,spinBlock)
!
      end function mqc_nef_r4_at
!
!
!     PROCEDURE MQC_NEF_Vector_Dimension_Expand
!
!
!>    \brief <b> MQC_NEF_Vector_Dimension_Expand is used to return an expanded
!>     NEF Matrix object based on a NEF Vector object. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Dimension_Expand is used to return an expanded
!>     NEF Matrix object based on a NEF Vector object.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFV
!>    \verbatim
!>        NEFV is type(mqc_nef_vector) 
!>        NEF Vector Object to be expanded.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        spinType is Character(Len=*),Optional 
!>        spinType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin 
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_dimention_expand(&
         NEFV,map_matrix_row,map_matrix_col,spinType) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_vector),intent(in)::NEFV
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in)::map_matrix_row,map_matrix_col
      character(len=*),intent(in),optional::spinType
      character(len=64)::mySpinType
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,i,j,id
!
      if(present(spinType)) then
        mySpinType = spinType 
      else
        mySpinType = 'all' 
      end if
!
      call NEFMatOut%init(map_matrix_row,initial_map_matrix_col=map_matrix_col)
      select case(mySpinType)
      case ('singlet')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            id = NEFMatOut%blockIndex([1,row],[2,col])
            do i = 1, NEFMatOut%Size_vector_alpha_row(row)
              do j = 1, NEFMatOut%Size_vector_alpha_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','a',row,col,'singlet')),'alpha'),&
                    i,j,'alpha')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','b',row,col,'singlet')),'beta'),&
                    i,j,'beta')
              end do
            end do
          end do
        end do
      case ('triplet')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            id = nefMatOut%blockIndex([1,row],[2,col])
            do i = 1, nefMatOut%Size_vector_alpha_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','b',row,col,'triplet')),'alpha'),&
                    i,j,'alphabeta')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','a',row,col,'triplet')),'beta'),&
                    i,j,'betaalpha')
              end do
            end do
          end do
        end do
      case ('all','full')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            id = nefMatOut%blockIndex([1,row],[2,col])
            do i = 1, nefMatOut%Size_vector_alpha_row(row)
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','a',row,col,'full')),'all'),&
                    i,j,'alpha')
              end do
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','b',row,col,'full')),'alpha'),&
                    i,j,'alphabeta')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','b',row,col,'full')),'beta'),&
                    i,j,'beta')
              end do
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','b',row,col,'full')),'beta'),&
                    i,j,'beta')
              end do
            end do
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_vector_dimention_excpanded: Unrecognized excitation type',6,&
          'spinType',spinType)
      end select  
!
      end function mqc_nef_vector_dimention_expand
!
!
!     PROCEDURE MQC_NEF_Matrix_Dimension_Collapse
!
!
!>    \brief <b> MQC_NEF_Matrix_Dimension_Collapse is used to return a collapse
!>    two indexes of the input NEF Matrix object into one. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Dimension_Collapse is used to return a collapse two indexes
!>    of the input NEF Matrix object into one.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFM
!>    \verbatim
!>        NEFM is type(mqc_nef_matrix) 
!>        NEF Matrix Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        spinType is Character(Len=*),Optional 
!>        spinType modiffies the excitations considered
!>          Options:
!>            'singlet' => return same spin 
!>            'triplet  => return opposite spin
!>            'all'     => return both types (default)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_dimention_collapse(NEFM,spinType) result(NEFVecOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFM
      character(len=*),intent(in),optional::spinType
      character(len=64)::mySpinType
      type(mqc_nef_vector)::NEFVecOut
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,i,j,id
!
      if(present(spinType)) then
        mySpinType = spinType 
      else
        mySpinType = 'all' 
      end if
!
      call NEFVecOut%init(NEFM%collapseMap(mySpinType))
      select case(mySpinType)
      case ('singlet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                call NEFVecOut%Main_vector(id)%put(&
                  NEFM%at(i,j,row,col,'alpha'),&
                  int(NEFM%collapseIndex(i,j,'a','a',row,col,'singlet')),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call NEFVecOut%Main_vector(id)%put(&
                  NEFM%at(i,j,row,col,'beta'),&
                  int(NEFM%collapseIndex(i,j,'b','b',row,col,'singlet')),&
                  'beta')
              end do
            end do
          end do
        end do
      case ('triplet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call NEFVecOut%Main_vector(id)%put(&
                  NEFM%at(i,j,row,col,'alphabeta'),&
                  int(NEFM%collapseIndex(i,j,'a','b',row,col,'triplet')),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                call NEFVecOut%Main_vector(id)%put(&
                  NEFM%at(i,j,row,col,'betaalpha'),&
                  int(NEFM%collapseIndex(i,j,'b','a',row,col,'triplet')),&
                  'beta')
              end do
            end do
          end do
        end do
      case ('all')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                call NEFVecOut%Main_vector(id)%put(&
                  NEFM%at(i,j,row,col,'alpha'),&
                  int(NEFM%collapseIndex(i,j,'a','a',row,col,'all')),&
                  'alpha')
              end do
              do j = 1, nefM%Size_vector_beta_col(col)
                call NEFVecOut%Main_vector(id)%put(&
                  NEFM%at(i,j,row,col,'alphabeta'),&
                  int(NEFM%collapseIndex(i,j,'a','b',row,col,'all')),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call NEFVecOut%Main_vector(id)%put(&
                  NEFM%at(i,j,row,col,'beta'),&
                  int(NEFM%collapseIndex(i,j,'b','b',row,col,'all')),&
                  'beta')
              end do
              do j = 1, nefM%Size_vector_alpha_col(col)
                call NEFVecOut%Main_vector(id)%put(&
                  NEFM%at(i,j,row,col,'betaalpha'),&
                  int(NEFM%collapseIndex(i,j,'b','a',row,col,'all')),&
                  'beta')
              end do
            end do
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_matrix_dimention_collapse: Unrecognized excitation type',6,&
          'spinType',spinType)
      end select  
!
      end function mqc_nef_matrix_dimention_collapse
!
!
!     PROCEDURE MQC_NEF_R4_Dimension_Collapse
!
!
!>    \brief <b> MQC_NEF_R4_Dimension_Collapse is used to return a collapse
!>    two indexes of the input NEF R4 object into one. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Dimension_Collapse is used to return a collapse two indexes
!>    of the input NEF R4 object into one.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFr4 is type(mqc_nef_r4) 
!>        NEF R4 Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        spinType is Character(Len=*),Optional 
!>        spinType modiffies the excitations considered
!>          Options:
!>            'singlet' => return same spin 
!>            'triplet  => return opposite spin
!>            'all'     => return both types (default)
!>    \endverbatim
!>
!>    \param[in] spinTypeKet
!>    \verbatim
!>        spinTypeKet is Character(Len=*),Optional 
!>        spinType modiffies the excitations considered over ket
!>          Options:
!>            'singlet' => return same spin 
!>            'triplet  => return opposite spin
!>            'all'     => return both types (default)
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is integer(kind=int64),dimension(4),intent(in),optional
!>        Order contains the indexes to be collapsed following:
!>        Order = [i,j,k,l] => i->j = p and  k->l = q
!>        (default [1,2,3,4])
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_dimention_collapse(NEFr4,spinType,spinTypeKet,Order) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFr4
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      integer(kind=int64),dimension(4)::my_order
      type(mqc_nef_matrix)::NEFMatOut
      character(len=64)::mySpinType,mySpinTypeKet
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,i,j,k,l,id,idij,idlk,fullid,finalid
      integer(kind=int64)::ia,ja,ka,la
      integer(kind=int64)::dimi,dimj,dimk,diml
!
      if(present(spinType)) then
        mySpinType = spinType 
        if(present(spinTypeKet)) then
          mySpinTypeKet = spinTypeKet 
        else
          mySpinTypeKet = spinType 
        end if
      else
        if(present(spinTypeKet)) then
          mySpinTypeKet = spinTypeKet 
          mySpinType = spinTypeKet
        else
          mySpinTypeKet = 'all'
          mySpinType = 'all' 
        end if
      end if
!
      if(present(Order)) then
        ia=Order(1)
        ja=Order(2)
        ka=Order(3)
        la=Order(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ka', ka,'la',la)
        my_order = Order
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_order = [1,2,3,4]
      endif
!
      call NEFr4%print(6,'initial R4')
      call NEFMatOut%init(NEFr4%collapseMap(mySpinType,[ia,ja]),&
         initial_map_matrix_col=NEFr4%collapseMap(mySpinTypeKet,[ka,la]))
      call NEFMatOut%print(6,'NEFMat Initial')
      select case(mySpinType)
      case ('singlet')
        do dimi = 1, NEFr4%blockSize(ia)
         do dimj = 1, NEFr4%blockSize(ja)
          idij = NEFr4%blockSize(ja)*(dimi-1) + dimj
          do dimk = 1, NEFr4%blockSize(ka)
           do diml = 1, NEFr4%blockSize(la)
            idlk = NEFr4%blockSize(la)*(dimk-1) + diml
            finalid = nefMatOut%blockIndex([1,idij],[2,idlk])
            write(*,*) idij, idlk, finalid
            do i = 1, NEFr4%blockSize(ia,dimi,'alpha')
             do j = 1, NEFr4%blockSize(ja,dimj,'alpha')
              write(*,*) 'a a a a'
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                write(*,*) dimi,dimj,dimk,diml,i,j,k,l
                call NEFR4%Main_R4(nefr4%blockIndex([ia,dimi],[ja,dimj],[la,diml],[ka,dimk]))%print(&
                    6,'MMFM: nefR4 out - aaaa')
                call NEFMatOut%Main_Matrix(finalid)%print(6,'MMFM: nef out - aa')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,1],[ka,dimk,k,1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'singlet')),&
                 int(NEFr4%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'singlet')),&
                 'alpha')
               end do
              end do
              write(*,*) 'a a b b'
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,1],[ka,dimk,k,-1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'singlet')),&
                 int(NEFr4%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'singlet')),&
                 'alphabeta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'beta')
             do j = 1, NEFr4%blockSize(ja,dimj,'beta')
              write(*,*) 'b b a a'
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                write(*,*) dimi,dimj,dimk,diml,i,j,k,l
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,-1],[ka,dimk,k,1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'singlet')),&
                 int(NEFr4%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'singlet')),&
                 'betaalpha')
               end do
              end do
              write(*,*) 'b b b b'
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                write(*,*) dimi,dimj,dimk,diml,i,j,k,l
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,-1],[ka,dimk,k,-1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'singlet')),&
                 int(NEFr4%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'singlet')),&
                 'beta')
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case ('triplet')
        do dimi = 1, NEFr4%blockSize(ia)
         do dimj = 1, NEFr4%blockSize(ja)
          idij = nefr4%blockIndex([1,1],[2,1],[3,dimi],[4,dimj])
          do dimk = 1, NEFr4%blockSize(ka)
           do diml = 1, NEFr4%blockSize(la)
            idlk = nefr4%blockIndex([1,1],[2,1],[3,diml],[4,dimk])
            finalid = nefMatOut%blockIndex([1,idij],[2,idlk])
            do i = 1, NEFr4%blockSize(ia,dimi,'alpha')
             do j = 1, NEFr4%blockSize(ja,dimj,'beta')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,-1],[ka,dimk,k,1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'triplet')),&
                 int(NEFr4%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'triplet')),&
                 'alpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,-1],[ka,dimk,k,-1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'triplet')),&
                 int(NEFr4%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'triplet')),&
                 'alphabeta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'beta')
             do j = 1, NEFr4%blockSize(ja,dimj,'alpha')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,1],[ka,dimk,k,1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'triplet')),&
                 int(NEFr4%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'triplet')),&
                 'betaalpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,1],[ka,dimk,k,-1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'triplet')),&
                 int(NEFr4%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'triplet')),&
                 'beta')
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case ('all')
        do dimi = 1, NEFr4%blockSize(ia)
         do dimj = 1, NEFr4%blockSize(ja)
          idij = nefr4%blockIndex([1,1],[2,1],[3,dimi],[4,dimj])
          do dimk = 1, NEFr4%blockSize(ka)
           do diml = 1, NEFr4%blockSize(la)
            idlk = nefr4%blockIndex([1,1],[2,1],[3,diml],[4,dimk])
            finalid = nefMatOut%blockIndex([1,idij],[2,idlk])
            do i = 1, NEFr4%blockSize(ia,dimi,'alpha')
             do j = 1, NEFr4%blockSize(ja,dimj,'alpha')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,1],[ka,dimk,k,1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'all')),&
                 'alpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,1],[ka,dimk,k,-1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'all')),&
                 'alphabeta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'beta')
             do j = 1, NEFr4%blockSize(ja,dimj,'beta')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,-1],[ka,dimk,k,1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'all')),&
                 'betaalpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,-1],[ka,dimk,k,-1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'all')),&
                 'beta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'alpha')
             do j = 1, NEFr4%blockSize(ja,dimj,'beta')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,-1],[ka,dimk,k,1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'all')),&
                 'alpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,-1],[ka,dimk,k,-1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'all')),&
                 'alphabeta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'beta')
             do j = 1, NEFr4%blockSize(ja,dimj,'alpha')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,1],[ka,dimk,k,1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'all')),&
                 'betaalpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,1],[ka,dimk,k,-1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'all')),&
                 'beta')
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case default
        call mqc_error_a('mqc_nef_r4_dimention_collapse: Unrecognized excitation type',6,&
          'spinType',spinType)
      end select  
      write(*,*) 'MMFM: <DONE mqc_nef_r4_dimention_collapse'
!
      end function mqc_nef_r4_dimention_collapse
!
!
!     PROCEDURE MQC_NEF_Matrix_Collapse_Index
!
!
!>    \brief <b> MQC_NEF_Matrix_Collapse_Index is used to return a collapsed
!>    two index excitation. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Collapse_Index is used to return a collapsed two index
!>    two index excitation.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFM
!>    \verbatim
!>        NEFM is type(mqc_nef_matrix) 
!>        NEF Matrix Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The first dimension of the element in MQC NEF Matrix
!>          If I>0 row count is from first index
!>    \endverbatim
!>    
!>    \param[in] J  
!>    \verbatim
!>        J is Integer(kind=int64)
!>        The second dimension of the element in MQC NEF Matrix
!>          If J>0 col count is from first index
!>    \endverbatim
!>    
!>    \param[in] i_spin 
!>    \verbatim
!>        i_spin is character(len=1)
!>        Spin of the annihilation operator natorb.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>        Options:
!>        a                 =>    alpha spin
!>        b                 =>    beta spin
!>    \endverbatim
!>
!>    \param[in] j_spin
!>    \verbatim
!>        j_spin is character(len=1)
!>        Spin of the creation operator natorb.
!>        If a_index or a_index is negative, this option does not 
!>        change the result.
!>        Options:
!>        a                 =>    alpha spin
!>        b                 =>    beta spin
!>    \endverbatim
!>
!>    \param[in] i_space
!>    \verbatim
!>        i_space is an integer
!>        Index of the annihilation operator natorb subspace.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>    \endverbatim
!>
!>    \param[in] j_space
!>    \verbatim
!>        j_space is an integer
!>        Index of the creation operator natorb subspace.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>    \endverbatim
!>
!>    \param[in] excType
!>    \verbatim
!>        excType is Character(Len=*),Optional 
!>        excType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin 
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_collapse_index(NEFM,i,j,i_spin,j_spin,i_space,j_space,excType) result(idOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFM
      integer::idOut
      integer(kind=int64)::i,j,i_space,j_space
      character(len=1)::i_spin,j_spin
      character(len=*),intent(in),optional::excType
      character(len=64)::selectedExcType
!
    if(i_spin.ne.'a'.and.i_spin.ne.'b') &
      call mqc_Error_a('mqc_nef_matrix_collapse_index: Unrecognized spin',&
           6,'i_spin',i_spin)
    if(j_spin.ne.'a'.and.j_spin.ne.'b') &
      call mqc_Error_a('mqc_nef_matrix_collapse_index: Unrecognized spin',&
            6,'j_spin',j_spin)
    if(i_space.gt.nefM%NumSpace_row.or.i_space.lt.1) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: Space index out-of-bound',&
            6,'i_space',i_space)
    if(j_space.gt.nefM%NumSpace_col.or.j_space.lt.1)&
      call mqc_Error_i('mqc_nef_matrix_collapse_index: Space index out-of-bound',&
            6,'j_space',j_space)
!
    idOut = 0
    if ((i_spin.eq.'b'.and.nefM%Size_vector_beta_row(i_space).eq.0).or.&
        (i_spin.eq.'a'.and.nefM%Size_vector_alpha_row(i_space).eq.0).or.&
        (j_spin.eq.'b'.and.nefM%Size_vector_beta_col(j_space).eq.0).or.& 
        (j_spin.eq.'a'.and.nefM%Size_vector_alpha_col(j_space).eq.0)) return
!
    if (i_spin.eq.'b'.and.i.gt.nefM%Size_vector_beta_row(i_space).or.i.lt.0) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: i beta index out-of-bound',&
            6,'i',i,'nefM%Size_vector_beta_row(i_space)',nefM%Size_vector_beta_row(i_space))
    if (i_spin.eq.'a'.and.i.gt.nefM%Size_vector_alpha_row(i_space).or.i.lt.0) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: i alpha index out-of-bound',&
            6,'i',i,'nefM%Size_vector_alpha_row(i_space)',nefM%Size_vector_alpha_row(i_space))
!
    if (j_spin.eq.'b'.and.j.gt.nefM%Size_vector_beta_col(j_space).or.j.lt.0) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: j beta index out-of-bound',&
            6,'j',j,'nefM%Size_vector_beta_col(j_space)',nefM%Size_vector_beta_col(j_space))
    if (j_spin.eq.'a'.and.j.gt.nefM%Size_vector_alpha_col(j_space).or.j.lt.0) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: j alpha index out-of-bound',&
            6,'j',j,'nefM%Size_vector_alpha_col(j_space)',nefM%Size_vector_alpha_col(j_space))
!
    if(present(excType)) then
      selectedExcType = trim(excType) !singlet, triplet, all
    else
      selectedExcType = 'singlet'
    end if
    if (selectedExcType.eq.'singlet'.and.i_spin.ne.j_spin) &
      call mqc_Error('inter_matrix_get_natural_single_excitation_index: wrong spin for triplet',6)
    if (selectedExcType.eq.'triplet'.and.i_spin.eq.j_spin) &
      call mqc_Error('inter_matrix_get_natural_single_excitation_index: wrong spin for triplet',6)
!
    idOut = 0
    select case(selectedExcType)
    case ('singlet')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%Size_vector_beta_col(j_space) + j
      else
        idOut = (i-1)*nefM%Size_vector_alpha_col(j_space) + j
      end if 
    case ('triplet')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%Size_vector_alpha_col(j_space) + j
      else
        idOut = (i-1)*nefM%Size_vector_beta_col(j_space) + j
      end if 
    case ('all','full')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%Size_vector_beta_col(j_space) + j
      else
        idOut = (i-1)*(nefM%Size_vector_beta_col(j_space)+nefM%Size_vector_alpha_col(j_space)) + &
              j + nefM%Size_vector_beta_col(j_space) 
      end if 
    case default
      call mqc_error_a('mqc_nef_matrix_collapse_index: Unrecognized',6,'selectedExcType',&
        selectedExcType)
    end select
!
      end function mqc_nef_matrix_collapse_index
!
!
!     PROCEDURE MQC_NEF_R4_Collapse_Index
!
!
!>    \brief <b> MQC_NEF_R4_Collapse_Index is used to return a collapsed
!>    two index excitation. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Collapse_Index is used to return a collapsed two index
!>    two index excitation.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is type(mqc_nef_r4) 
!>        NEF R4 Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The dimi of the element in MQC NEF R4
!>          If I>0 row count is from first index
!>    \endverbatim
!>    
!>    \param[in] J  
!>    \verbatim
!>        J is Integer(kind=int64)
!>        The dimj of the element in MQC NEF MR4
!>          If J>0 col count is from first index
!>    \endverbatim
!>    
!>    \param[in] i_spin 
!>    \verbatim
!>        i_spin is character(len=1)
!>        Spin of the annihilation operator natorb.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>        Options:
!>        a                 =>    alpha spin
!>        b                 =>    beta spin
!>    \endverbatim
!>
!>    \param[in] j_spin
!>    \verbatim
!>        j_spin is character(len=1)
!>        Spin of the creation operator natorb.
!>        If a_index or a_index is negative, this option does not 
!>        change the result.
!>        Options:
!>        a                 =>    alpha spin
!>        b                 =>    beta spin
!>    \endverbatim
!>
!>    \param[in] i_space
!>    \verbatim
!>        i_space is an integer
!>        Index of the annihilation operator natorb subspace.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>    \endverbatim
!>
!>    \param[in] j_space
!>    \verbatim
!>        j_space is an integer
!>        Index of the creation operator natorb subspace.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>    \endverbatim
!>
!>    \param[in] dimi
!>    \verbatim
!>        dimi is an integer
!>        Index of the creation operator natorb subspace.
!>    \endverbatim
!>
!>    \param[in] dimj
!>    \verbatim
!>        dimj is an integer
!>        Index of the creation operator natorb subspace.
!>    \endverbatim
!>
!>    \param[in] excType
!>    \verbatim
!>        excType is Character(Len=*),Optional 
!>        excType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin 
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_collapse_index(NEFM,i,j,i_spin,j_spin,&
                                          i_space,j_space,dimi,dimj,excType) result(idOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFM
      integer::idOut
      integer(kind=int64)::i,j,i_space,j_space,dimi,dimj
      character(len=1)::i_spin,j_spin
      character(len=*),intent(in),optional::excType
      character(len=64)::selectedExcType
!
    if(dimi.le.0.or.dimi.gt.4) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: dimention index out-of-bound',&
           6,'dimi',dimi)
    if(dimj.le.0.or.dimj.gt.4) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: dimention index out-of-bound',&
           6,'dimj',dimj)
!
    if(i_spin.ne.'a'.and.i_spin.ne.'b') &
      call mqc_Error_a('mqc_nef_r4_collapse_index: Unrecognized spin',&
           6,'i_spin',i_spin)
    if(j_spin.ne.'a'.and.j_spin.ne.'b') &
      call mqc_Error_a('mqc_nef_r4_collapse_index: Unrecognized spin',&
            6,'j_spin',j_spin)
    if(i_space.gt.nefM%blockSize(dimi).or.i_space.lt.1) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: Space index out-of-bound',&
            6,'i_space',i_space)
    if(j_space.gt.nefM%blockSize(dimj).or.j_space.lt.1)&
      call mqc_Error_i('mqc_nef_r4_collapse_index: Space index out-of-bound',&
            6,'j_space',j_space)
!
    idOut = 0
    if ((i_spin.eq.'b'.and.nefM%blockSize(dimi,i_space,'beta').eq.0).or.&
        (i_spin.eq.'a'.and.nefM%blockSize(dimi,i_space,'alpha').eq.0).or.&
        (j_spin.eq.'b'.and.nefM%blockSize(dimj,j_space,'beta').eq.0).or.& 
        (j_spin.eq.'a'.and.nefM%blockSize(dimj,j_space,'alpha').eq.0)) return
!
    if (i_spin.eq.'b'.and.i.gt.nefM%blockSize(dimi,i_space,'beta').or.i.lt.0) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: i beta index out-of-bound',&
            6,'i',i,'nefM%blockSize(dimi,i_space,beta)',nefM%blockSize(dimi,i_space,'beta'))
    if (i_spin.eq.'a'.and.i.gt.nefM%blockSize(dimi,i_space,'alpha').or.i.lt.0) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: i alpha index out-of-bound',&
            6,'i',i,'nefM%blockSize(dimi,i_space,alpha)',nefM%blockSize(dimi,i_space,'alpha'))
!
    if (j_spin.eq.'b'.and.j.gt.nefM%blockSize(dimj,j_space,'beta').or.j.lt.0) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: j beta index out-of-bound',&
            6,'j',j,'nefM%blockSize(dimj,j_space,beta)',nefM%blockSize(dimj,j_space,'beta'))
    if (j_spin.eq.'a'.and.j.gt.nefM%blockSize(dimj,j_space,'alpha').or.j.lt.0) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: j alpha index out-of-bound',&
            6,'j',j,'nefM%blockSize(dimj,j_space,alpha)',nefM%blockSize(dimj,j_space,'alpha'))
!
    if(present(excType)) then
      selectedExcType = trim(excType) !singlet, triplet, all
    else
      selectedExcType = 'singlet'
    end if
    if (selectedExcType.eq.'singlet'.and.i_spin.ne.j_spin) &
      call mqc_Error('inter_r4_get_natural_single_excitation_index: wrong spin for triplet',6)
    if (selectedExcType.eq.'triplet'.and.i_spin.eq.j_spin) &
      call mqc_Error('inter_r4_get_natural_single_excitation_index: wrong spin for triplet',6)
!
    idOut = 0
    select case(selectedExcType)
    case ('singlet')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'beta') + j
      else
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'alpha') + j
      end if 
    case ('triplet')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'alpha') + j
      else
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'beta') + j
      end if 
    case ('all','full')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'beta') + j
      else
        idOut = (i-1)*(nefM%blockSize(dimj,j_space,'beta') + &
                       nefM%blockSize(dimj,j_space,'alpha')) + &
                   j + nefM%blockSize(dimj,j_space,'beta')
      end if 
    case default
      call mqc_error_a('mqc_nef_r4_collapse_index: Unrecognized',6,'selectedExcType',&
        selectedExcType)
    end select
!
      end function mqc_nef_r4_collapse_index
!
!
!     PROCEDURE MQC_NEF_Matrix_Collapse_Map
!
!
!>    \brief <b> MQC_NEF_Matrix_Collapse_Map is used to return a collapsed
!>    two Map_Vectors. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Collapse_Map is used to return a collapsed two Map_Vectors.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFM
!>    \verbatim
!>        NEFM is type(mqc_nef_matrix) 
!>        NEF Matrix Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] excType
!>    \verbatim
!>        excType is Character(Len=*),Optional 
!>        excType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin 
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_collapse_map(NEFM,excType) result(MapOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFM
      character(len=*),intent(in),optional::excType
      character(len=64)::selectedExcType
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::MapOut
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map_Col,New_Map_Row
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,id,idOld,r0
      logical::full_row,full_col
!
      if(present(excType)) then
        selectedExcType = trim(excType) !singlet, triplet, all
      else
        selectedExcType = 'singlet'
      end if
!
      if(allocated(MapOut)) deallocate(MapOut)
      allocate(MapOut(NEFM%NumSpace_row*NEFM%NumSpace_col))
!
      select case(selectedExcType)
      case ('singlet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            call MapOut(id)%init(&
               nefM%Size_vector_alpha_col(col)*nefM%Size_vector_alpha_row(row),0,&
               nefM%Size_vector_beta_col(col)*nefM%Size_vector_beta_row(row))
          end do
        end do
      case ('triplet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            call MapOut(id)%init(&
               nefM%Size_vector_alpha_col(col)*nefM%Size_vector_beta_row(row),0,&
               nefM%Size_vector_beta_col(col)*nefM%Size_vector_alpha_row(row))
          end do
        end do
      case ('all','full')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            call MapOut(id)%init(&
               nefM%Size_vector_alpha_col(col)*nefM%Size_vector_alpha_row(row)+&
               nefM%Size_vector_alpha_col(col)*nefM%Size_vector_beta_row(row),0,&
               nefM%Size_vector_beta_col(col)*nefM%Size_vector_beta_row(row)+&
               nefM%Size_vector_beta_col(col)*nefM%Size_vector_alpha_row(row))
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_matrix_collapse_map: Unrecognized excitation type',6,&
          'selectedExcType',selectedExcType)
      end select  
!
      end function mqc_nef_matrix_collapse_map
!
!
!     PROCEDURE MQC_NEF_R4_Collapse_Map
!
!
!>    \brief <b> MQC_NEF_R4_Collapse_Map is used to return a collapsed
!>    two Map_Vectors. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Collapse_Map is used to return a collapsed two Map_Vectors.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFr4
!>    \verbatim
!>        NEFr4 is type(mqc_nef_r4) 
!>        NEF R4 Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] excType
!>    \verbatim
!>        excType is Character(Len=*),Optional 
!>        excType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin 
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is integer(kind=int64),dimension(2),intent(in),optional
!>        Order contains the indexes to be collapsed following:
!>        Order = [i,j] => i->j = p, where i and j are the dimentions
!>        (default [1,2])
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_collapse_map(NEFr4,excType,Order) result(MapOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFr4
      character(len=*),intent(in),optional::excType
      character(len=64)::selectedExcType
      integer(kind=int64),dimension(2),intent(in),optional::Order
      integer(kind=int64),dimension(2)::my_order
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::MapOut
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map_Col,New_Map_Row
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::ia,ja,ka,la,row,col,id
      logical::full_row,full_col
!
      if(present(excType)) then
        selectedExcType = trim(excType) !singlet, triplet, all
      else
        selectedExcType = 'singlet'
      end if
!
      if(present(Order)) then
        ia=Order(1)
        ja=Order(2)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ia', ia,'ja',ja)
        my_order = Order
      else
        ia=1
        ja=2
        my_order = [1,2]
      endif
!
      if(allocated(MapOut)) deallocate(MapOut)
      allocate(MapOut(int(NEFr4%blockSize(ia)*NEFr4%blockSize(ja))))
!
      select case(selectedExcType)
      case ('singlet')
        do row = 1, NEFr4%blockSize(ia)
          do col = 1, NEFr4%blockSize(ja)
            id = nefR4%blockIndex([1,1],[2,1],[3,row],[4,col])
            call MapOut(id)%init(&
               nefR4%blockSize(ia,row,'alpha')*nefR4%blockSize(ja,col,'alpha'),0,&
               nefR4%blockSize(ia,row,'beta')*nefR4%blockSize(ja,col,'beta'))
          end do
        end do
      case ('triplet')
        do row = 1, NEFr4%blockSize(ia)
          do col = 1, NEFr4%blockSize(ja)
            id = nefr4%blockIndex([1,1],[2,1],[3,row],[4,col])
            call MapOut(id)%init(&
               nefR4%blockSize(ia,row,'alpha')*nefR4%blockSize(ja,col,'beta'),0,&
               nefR4%blockSize(ia,row,'beta')*nefR4%blockSize(ja,col,'alpha'))
          end do
        end do
      case ('all','full')
        do row = 1, NEFr4%blockSize(ia)
          do col = 1, NEFr4%blockSize(ja)
            id = nefr4%blockIndex([1,1],[2,1],[3,row],[4,col])
            call MapOut(id)%init(&
               nefR4%blockSize(ia,row,'alpha')*(nefR4%blockSize(ja,col,'beta')+&
                                                nefR4%blockSize(ja,col,'alpha')),0,&
               nefR4%blockSize(ia,row,'beta')*(nefR4%blockSize(ja,col,'alpha')+&
                                               nefR4%blockSize(ja,col,'beta')))
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_r4_collapse_map: Unrecognized excitation type',6,&
          'selectedExcType',selectedExcType)
      end select  
!
      end function mqc_nef_r4_collapse_map
!
!
!     PROCEDURE mqc_nef_vector_output_block_vector
!
!>    \brief <b> mqc_nef_vector_output_block_vector is used to output a MQC Vector of an MQC 
!>    NEF Vector type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    mqc_nef_vector_output_block_vector is used to output a MQC Vector of an MQC 
!>    NEF Vector type variable
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] nefV
!>    \verbatim
!>        nefV is Class(MQC_NEF_Vector)
!>        The MQC_NEF_Vecotr variable to output.
!>    \endverbatim
!>    
!>    \param[in] indexes 
!>    \verbatim
!>        indexes is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the blocks to be returned.
!>        Default: all indexes
!>    \endverbatim
!>    
!>    \param[in] BlockName
!>    \verbatim
!>        BlockName is character(len=*),optional
!>        = 'full':        outputs the full MQC_NEF_Vector.(default) 
!>        = 'alpha','a':   outputs the alpha spin block.
!>        = 'beta','b':    outputs the beta spin block. 
!>    \endverbatim
!>    
!>    \param[in] BlockOrder
!>    \verbatim
!>        BlockOrder is character(len=*),optional
!>        = 'beta': orders the output as b->b, b->a, a->b, a->a (default)
!>        = 'alpha':orders the output as a->b, a->a, b->b, b->a
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_output_block_vector(nefV,indexes,blockName,BlockOrder) result(vectorOut)
!
      implicit none
      class(mqc_nef_vector)::nefV
      integer(kind=int64),dimension(:),optional::indexes
      character(len=*),optional::blockOrder
      character(len=*),optional::blockName
      type(mqc_vector)::vectorOut
      character(len=64)::myBlockName,myBlockOrder
      integer(kind=int64),dimension(:),allocatable::myId
      logical::full
      integer(kind=int64)::fsize,i,j
!
      if(present(blockOrder)) then
        call string_change_case(blockOrder,'l',myBlockOrder)
      else
        myBlockOrder = 'beta'
      endIf
      if(myBlockOrder.ne.'alpha'.and.myBlockOrder.ne.'beta') &
        call mqc_error_A('BlockOrder inrecognized in mqc_nef_vector_output_block_vector',6,&
           'myBlockOrder',myBlockOrder)
      if(present(blockName)) then
        call string_change_case(blockName,'l',myBlockName)
      else
        myBlockName = 'full'
      endIf
      if(myBlockName.ne.'alpha'.and.myBlockName.ne.'beta'.and.myBlockName.ne.'full') &
        call mqc_error_A('BlockName inrecognized in mqc_nef_vector_output_block_vector',6,&
           'myBlockName',myBlockName)
      if(present(indexes)) then
        myid = indexes
        full=.false.
      else
        full=.true.
      end if
!
      fsize = 0
      if (full) then
        do i = 1, nefV%NumSpace
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsize = fsize + nefV%Size_vector_alpha(i)
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsize = fsize + nefV%Size_vector_beta(i)
        end do
      else
        do i = 1, size(myid)
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsize = fsize + nefV%Size_vector_alpha(myid(i))
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsize = fsize + nefV%Size_vector_beta(myid(i))
        end do
      end if
      call vectorOut%init(fsize)
      if (full) then
        fsize = 1
        do i = 1, nefV%NumSpace
          if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then      
            if (nefV%Main_vector(i)%blockSize('full').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('full'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('full')
          else if (myBlockName.eq.'full') then
            if (nefV%Main_vector(i)%blockSize('beta').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('beta'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('beta')
            if (nefV%Main_vector(i)%blockSize('alpha').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('alpha'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('alpha')
          else if (myBlockName.eq.'alpha') then
            if (nefV%Main_vector(i)%blockSize('alpha').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('alpha'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('alpha')
          else
            if (nefV%Main_vector(i)%blockSize('beta').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('beta'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('beta')
          end if
        end do
      else
        do j = 1, size(myid)
          i = myid(j)
          if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then
            if (nefV%Main_vector(i)%blockSize('full').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('full'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('full')
          else if (myBlockName.eq.'full') then
            if (nefV%Main_vector(i)%blockSize('beta').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('beta'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('beta')
            if (nefV%Main_vector(i)%blockSize('alpha').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('alpha'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('alpha')
          else if (myBlockName.eq.'alpha') then
            if (nefV%Main_vector(i)%blockSize('alpha').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('alpha'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('alpha')
          else
            if (nefV%Main_vector(i)%blockSize('beta').gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('beta'),fsize)
            fsize = fsize + nefV%Main_vector(i)%blockSize('beta')
          end if
        end do
      end if
!
      end function mqc_nef_vector_output_block_vector
!
!
!     PROCEDURE mqc_nef_matrix_output_block_matrix
!
!>    \brief <b> mqc_nef_matrix_output_block_matrix is used to output a MQC Matrix of an MQC 
!>    NEF Matrix type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    mqc_nef_matrix_output_block_matrix is used to output a MQC Matrix of an MQC 
!>    NEF Matrix type variable
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] nefMat
!>    \verbatim
!>        nefMat is Class(MQC_NEF_Matrix)
!>        The MQC_NEF_Matrix variable to output.
!>    \endverbatim
!>    
!>    \param[in] indexes 
!>    \verbatim
!>        indexes is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the blocks to be returned.
!>        Default: all indexes
!>    \endverbatim
!>    
!>    \param[in] indexesKet
!>    \verbatim
!>        indexesKet is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the blocks to be returned.
!>        Default: all indexes
!>    \endverbatim
!>    
!>    \param[in] BlockName
!>    \verbatim
!>        BlockName is character(len=*),optional
!>        = 'full'          :   outputs the full MQC_NEF_Vector.(default) 
!>        = 'alpha','a','aa':   outputs the alpha spin block.
!>        = 'beta', 'b','bb':   outputs the beta spin block. 
!>        = 'alphabeta','ab':   outputs the alphabeta spin block.
!>        = 'betaalpha','ba':   outputs the betaalpha spin block. 
!>    \endverbatim
!>    
!>    \param[in] BlockOrder
!>    \verbatim
!>        BlockOrder is character(len=*),optional
!>        = 'beta': orders the output as b->b, b->a, a->b, a->a (default)
!>        = 'alpha':orders the output as a->b, a->a, b->b, b->a
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_output_block_matrix(nefMat,indexes,indexesKet,blockName,&
                                                  BlockOrder) result(MatrixOut)
!
      implicit none
      class(mqc_nef_matrix)::nefMat
      integer(kind=int64),dimension(:),optional::indexes,indexesKet
      character(len=*),optional::blockOrder
      character(len=*),optional::blockName
      type(mqc_Matrix)::MatrixOut
      character(len=64)::myBlockName,myBlockOrder
      integer(kind=int64),dimension(:),allocatable::myId,myIdKet
      logical::full,fullKet
      integer(kind=int64)::fsize,fsizeKet,i,j,k,p,size1,size2,size3,size4,id
!
      if(present(blockOrder)) then
        call string_change_case(blockOrder,'l',myBlockOrder)
      else
        myBlockOrder = 'beta'
      endIf
      if(myBlockOrder.ne.'alpha'.and.myBlockOrder.ne.'beta') &
        call mqc_error_A('BlockOrder inrecognized in mqc_nef_matrix_output_block_matrix',6,&
           'myBlockOrder',myBlockOrder)
      if(present(blockName)) then
        call string_change_case(blockName,'l',myBlockName)
      else
        myBlockName = 'full'
      endIf
      if(myBlockName.ne.'alpha'.and.myBlockName.ne.'beta'.and.myBlockName.ne.'full') &
        call mqc_error_A('BlockName inrecognized in mqc_nef_matrix_output_block_matrix',6,&
           'myBlockName',myBlockName)
      if(present(indexes)) then
        myid = indexes
        full = .false.
      else
        full = .true.
      end if
      if(present(indexesKet)) then
        myidKet = indexesKet
        fullKet = .false.
      else
        fullKet = .true.
      end if
!
      fsizeKet = 0
      if (full) then
        do i = 1, nefMat%NumSpace_col
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsizeKet = fsizeKet + nefMat%Size_vector_alpha_col(i)
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsizeKet = fsizeKet + nefMat%Size_vector_beta_col(i)
        end do
      else
        do i = 1, size(myidKet)
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsizeKet = fsizeKet + nefMat%Size_vector_alpha_col(myidKet(i))
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsizeKet = fsizeKet + nefMat%Size_vector_beta_col(myidKet(i))
        end do
      end if
      fsize = 0
      if (fullKet) then
        do i = 1, nefMat%NumSpace_row
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsize = fsize + nefMat%Size_vector_alpha_row(i)
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsize = fsize + nefMat%Size_vector_beta_row(i)
        end do
      else
        do i = 1, size(myid)
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsize = fsize + nefMat%Size_vector_alpha_row(myid(i))
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsize = fsize + nefMat%Size_vector_beta_row(myid(i))
        end do
      end if
      call matrixOut%init(fsize,fsizeKet)
      if (full) then
        fsize = 1
        do i = 1, nefMat%NumSpace_col
          if (fullKet) then
            fsizeKet = 1
            do j = 1, nefMat%NumSpace_row
              id = nefMat%blockIndex([1,i],[2,j])
              write(*,*) 'i = ',i,' j = ',j,' id = ',id
              if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then      
                size1 = nefMat%Main_matrix(id)%blockSize('full',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('full',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('full'),&
                                         [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'full') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                size3 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size4 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                call  mqc_print(nefMat%Main_matrix(id)%getBlock('beta'),6,'test beta')
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                call  mqc_print(nefMat%Main_matrix(id)%getBlock('betaalpha'),6,'test betaalpha')
                if (size3.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                               [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                call  mqc_print(nefMat%Main_matrix(id)%getBlock('alphabeta'),6,'test alphabeta')
                if (size1.gt.0.and.size4.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                               [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                call  mqc_print(nefMat%Main_matrix(id)%getBlock('alpha'),6,'test alpha')
                if (size3.gt.0.and.size4.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                   [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                fsizeKet = fsizeKet + size2 + size4
                size1 = size1 + size3
              else if (myBlockName.eq.'alpha') then
                size1 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'beta') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'alphabeta') then
                size1 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'betaalpha') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else
                call mqc_error_a('block name and order pair unrecognized in&
                                 & mqc_nef_matrix_output_block_matrix',6,&
                     'myBlockName',myBlockName,'myBlockOrder',myBlockOrder)
              end if
            end do
          else
            fsizeKet = 1
            do k = 1, size(myidKet)
              j = myidKet(k)
              id = nefMat%blockIndex([1,i],[2,j])
              if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then      
                size1 = nefMat%Main_matrix(id)%blockSize('full',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('full',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('full'),&
                                      [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'full') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                size3 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size4 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                         [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                if (size3.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                               [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                if (size1.gt.0.and.size4.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                               [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                if (size3.gt.0.and.size4.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                   [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                fsizeKet = fsizeKet + size2 + size4
                size1 = size1 + size3
              else if (myBlockName.eq.'alpha') then
                size1 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'beta') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'alphabeta') then
                size1 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'betaalpha') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else
                call mqc_error_a('block name and order pair unrecognized in&
                                 & mqc_nef_matrix_output_block_matrix',6,&
                     'myBlockName',myBlockName,'myBlockOrder',myBlockOrder)
              end if                                                                               
            end do
          end if  
          fsize = fsize + size1
        end do
      else
        fsize = 1
        do k = 1, size(myid)
          i = myid(k)
          if (fullKet) then
            fsizeKet = 1
            do j = 1, nefMat%NumSpace_row
              id = nefMat%blockIndex([1,i],[2,j])
              if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then      
                size1 = nefMat%Main_matrix(id)%blockSize('full',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('full',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('full'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'full') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                size3 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size4 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                         [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                if (size3.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                               [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                if (size1.gt.0.and.size4.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                               [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                if (size3.gt.0.and.size4.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                   [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                fsizeKet = fsizeKet + size2 + size4
                size1 = size1 + size3
              else if (myBlockName.eq.'alpha') then
                size1 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'beta') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'alphabeta') then
                size1 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'betaalpha') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else
                call mqc_error_a('block name and order pair unrecognized in&
                                 & mqc_nef_matrix_output_block_matrix',6,&
                     'myBlockName',myBlockName,'myBlockOrder',myBlockOrder)
              end if
            end do
          else
            fsizeKet = 1
            do p = 1, size(myidKet)
              j = myidKet(p)
              id = nefMat%blockIndex([1,i],[2,j])
              if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then      
                size1 = nefMat%Main_matrix(id)%blockSize('full',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('full',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('full'),&
                               [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'full') then
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('beta',2) 
                size3 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size4 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                if (size3.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                               [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                if (size1.gt.0.and.size4.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                               [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                if (size3.gt.0.and.size4.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                   [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                fsizeKet = fsizeKet + size2 + size4
                size1 = size1 + size3
              else if (myBlockName.eq.'alpha') then
                size1 = nefMat%Main_matrix(id)%blockSize('alpha',1) 
                size2 = nefMat%Main_matrix(id)%blockSize('alpha',2) 
                if (size1.gt.0.and.size2.gt.0) &
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'beta') then
                fsizeKet = fsizeKet + size2                                           
              else if (myBlockName.eq.'betaalpha') then                               
                size1 = nefMat%Main_matrix(id)%blockSize('beta',1)                    
                size2 = nefMat%Main_matrix(id)%blockSize('alpha',2)                   
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&  
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])     
                fsizeKet = fsizeKet + size2                                           
              else                                                                    
                call mqc_error_a('block name and order pair unrecognized in&
                                 & mqc_nef_matrix_output_block_matrix',6,&
                     'myBlockName',myBlockName,'myBlockOrder',myBlockOrder)
              end if                                                                               
            end do                                                                                 
          end if                                                                                   
          fsize = fsize + size1
        end do
      end if                                                                                   
      write(*,*) 'MMFM: <DONE mqc_nef_matrix_output_block_matrix'
!
      end function mqc_nef_matrix_output_block_matrix
      End Module MQC_NEF  
