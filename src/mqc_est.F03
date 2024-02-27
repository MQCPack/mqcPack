!>    \brief <b> MQC EST contains objects for the manipulation of intermediates used
!>    electronic structure theory </b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC EST contains objects that speed up and enhance the implementation of 
!>    electronic structure theories in Fortan. The objects allow for automatic 
!>    determination of spin structure and so enables code to be used without 
!>    modification whether a restricted, unrestricted or general formalism is used.
!>    The objects are able to store common intermediates in the calculation of 
!>    electronic structure, and carry their own procedures for manipulation of these 
!>    intermediates. The EST derived types defined in this package are:
!>
!>    1. MQC_SCF_Integral:       EST intermediate operator matrix representation
!>    2. MQC_SCF_Eigenvalues:    EST intermediate diagonal matrix variable
!>    3. MQC_Wavefunction:       EST intermediates of single determinant wavefunctions
!>    4. MQC_PSCF_Wavefunction:  EST intermediates of post-Hartree-Fock methods 
!>    5. MQC_Determinant:        Binary strings reresenting occupation number vectors
!>    6. MQC_TwoERIs:            Two-electron resonance integrals
!>
!>    This module is level 2 in the MQC hierarchy and so depends on level 0 and 1
!>    modules. 
!>
!>    \endverbatim
!
      Module MQC_EST  
!
!     **********************************************************************
!     **********************************************************************
!     **                                                                  **
!     **               The Merced Quantum Chemistry Package               **
!     **                            (MQCPack)                             **
!     **                                                                  **
!     ** Written By:                                                      **
!     **    Lee M. Thompson, Dave Mullaly, Xianghai Sheng, and Hrant P.   **
!     **    Hratchian                                                     **
!     **                                                                  **
!     **                         Version 24.02.1                          **
!     **                        Feburary 19, 2024                         **
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
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
      Use iso_fortran_env, only: int32, int64, real64
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     Single Reference Wavefunctions...
!
!>    \brief <b> EST intermediate operator matrix representation</b>
      Type MQC_SCF_Integral
!
!       Blocks are constructed in the order:
!
!              |' aa | ba '|
!              |-----|-----|
!              |. ab | bb .|
!
        Type(MQC_Matrix),Private::Alpha,Beta,AlphaBeta,BetaAlpha
        Character(Len=64),Private::Array_Name ! See below 
        Character(Len=64),Private::Array_Type ! Space,Spin,General
        Integer(kind=int64),Dimension(:),Allocatable,Private::Energy_List
      Contains 
!>      \brief <b> Print the MQC_SCF_Integral object</b>
        Procedure, Public::print => mqc_print_integral
        Procedure, Private::hasAlpha => mqc_integral_has_alpha
        Procedure, Private::hasBeta => mqc_integral_has_beta
        Procedure, Private::hasAlphaBeta => mqc_integral_has_alphabeta
        Procedure, Private::hasBetaAlpha => mqc_integral_has_betaalpha
        Procedure, Public::type => mqc_integral_array_type
        Procedure, Private::blockSize => mqc_integral_dimension 
!>      \brief <b> Return the label of the MQC_SCF_Integral object</b>
        Procedure, Public::getLabel => mqc_integral_array_name
!>      \brief <b> Set the label of the MQC_SCF_Integral object</b>
        Procedure, Public::addLabel => mqc_integral_add_name
!>      \brief <b> Return a spin block of MQC_SCF_Integral object as an MQC matrix</b>
        Procedure, Public::getBlock => mqc_integral_output_block
!>      \brief <b> Set an MQC_SCF_Integral object as the identity matrix</b>
        Procedure, Public::identity => mqc_integral_identity 
!>      \brief <b> Initializes an MQC_SCF_Integral object</b>
        Procedure, Public::init => mqc_integral_initialize 
!>      \brief <b> Diagonalizes an MQC_SCF_Integral object</b>
        Procedure, Public::diag => mqc_scf_integral_diagonalize
!>      \brief <b> Performs SVD on an MQC_SCF_Integral object</b>
        Procedure, Public::svd => mqc_scf_integral_svd
!>      \brief <b> Solves a generalized eigenvalue problem involving an 
!>      MQC_SCF_Integral object</b>
        Procedure, Public::eigensys => mqc_scf_integral_generalized_eigensystem
!>      \brief <b> Inverts an MQC_SCF_Integral object</b>
        Procedure, Public::inv => mqc_scf_integral_inverse
!>      \brief <b> Returns the trace of an MQC_SCF_Integral object</b>
        Procedure, Public::trace => mqc_scf_integral_trace
!>      \brief <b> Returns the determinant of an MQC_SCF_Integral object</b>
        Procedure, Public::det => mqc_scf_integral_determinant
!>      \brief <b> Returns the norm of an MQC_SCF_Integral object</b>
        Procedure, Public::norm => mqc_integral_norm
!>      \brief <b> Returns an MQC_SCF_Integral object raised to a power</b>
        Procedure, Public::power => mqc_scf_integral_power
        Procedure, Public::at => mqc_integral_at
        Procedure, Public::setEList => mqc_integral_set_energy_list
        Procedure, Public::getEList => mqc_integral_get_energy_list
        Procedure, Public::deleteEList => mqc_integral_delete_energy_list
!>      \brief <b> Returns a subset of an MQC_SCF_Integral object</b>
        Procedure, Public::orbitals => mqc_integral_output_orbitals 
!>      \brief <b> Swaps rows/columns of an MQC_SCF_Integral object</b>
        Procedure, Public::swap => mqc_integral_swap_orbitals
!>      \brief <b> Combines virtual and occupied MQC_SCF_Integral objects</b>
        Procedure, Public::combine => mqc_integral_combine_orbitals
!>      \brief <b> Swaps ab and ba blocks of an MQC_SCF_Integral objects</b>
        Procedure, Public::swapODB => mqc_integral_swap_offDiagBlocks
      End Type
!
!>    \brief <b> EST intermediate diagonal matrix variable</b>
!     MQC_SCF_Eigenvalues
      Type MQC_SCF_Eigenvalues
        Type(MQC_Vector),Private::Alpha,Beta
        Character(Len=64),Private::Array_Name ! See below 
        Character(Len=64),Private::Array_Type ! Space,Spin,General 
      Contains 
!>      \brief <b> Print the MQC_SCF_Eigenvalues object</b>
        Procedure, Public::print => mqc_print_eigenvalues
        Procedure, Private::hasAlpha => mqc_eigenvalues_has_alpha
        Procedure, Private::hasBeta => mqc_eigenvalues_has_beta
        Procedure, Private::type => mqc_eigenvalues_array_type
        Procedure, Private::blockSize => mqc_eigenvalues_dimension 
!>      \brief <b> Return the label of the MQC_SCF_Eigenvalues object</b>
        Procedure, Public::getLabel => mqc_eigenvalues_array_name
!>      \brief <b> Set the label of the MQC_SCF_Eigenvalues object</b>
        Procedure, Public::addLabel => mqc_eigenvalues_add_name
!>      \brief <b> Return a spin block of MQC_SCF_Eigenvalues object as an MQC matrix</b>
        Procedure, Public::getBlock => mqc_eigenvalues_output_block
!>      \brief <b> Returns values of an MQC_SCF_Eigenvalues object raised to a power</b>
        Procedure, Public::power => mqc_scf_eigenvalues_power
!>      \brief <b> Returns the value of a specified element in an MQC_SCF_Eigenvalues 
!>      object</b>
        Procedure, Public::at => mqc_eigenvalues_at
      End Type
!
!>    \brief <b> EST intermediates of single determinant wavefunctions</b>
!     MQC_Wavefunction
      Type MQC_Wavefunction
        Type(MQC_SCF_Integral)::MO_Coefficients
        Type(MQC_SCF_Eigenvalues)::MO_Energies
        Type(MQC_SCF_Eigenvalues)::MO_Symmetries
        Type(MQC_SCF_Integral)::Core_Hamiltonian
        Type(MQC_SCF_Integral)::Fock_Matrix
        Type(MQC_SCF_Integral)::Density_Matrix
        Type(MQC_SCF_Integral)::SCF_Density_Matrix
        Type(MQC_SCF_Integral)::Overlap_Matrix
        Type(MQC_Scalar)::NAlpha,NBeta,NElectrons,NBasis,Charge,Multiplicity
        Character(Len=256)::Basis,Symmetry,WF_Type
        Logical::WF_Complex
      Contains 
!>      \brief <b> Print the MQC_Wavefunction object</b>
        Procedure, Public::print => mqc_print_wavefunction
      End Type MQC_Wavefunction
!
!>    \brief <b> EST intermediates of one-electron basis</b>
!     MQC_Basis_Set
      Type MQC_Basis_Set
        Type(MQC_Matrix)::sh2AtMp
        Type(MQC_Matrix)::shlTyp
        Type(MQC_Matrix)::nPrmSh
        Type(MQC_Matrix)::prmExp
        Type(MQC_Matrix)::conCoef
        Type(MQC_Matrix)::conCoTwo
        Type(MQC_Matrix)::shCoor
      End Type MQC_Basis_Set
!
!
!     Post-SCF Wavefunctions...
!
!     Parent Type
!>    \brief <b> EST intermediates of post-Hartree-Fock methods</b>
      Type,Extends(MQC_Wavefunction)::MQC_PSCF_Wavefunction
        Integer(kind=int64)::NCore,NVal,NActive,NFrz
        Type(MQC_Matrix)::PSCF_Amplitudes
        Type(MQC_Vector)::PSCF_Energies
      End Type MQC_PSCF_Wavefunction
!
!    
!     Determinants....
!
!     MQC_Determinant_String
      Type,Private::MQC_Determinant_String
        Type(MQC_Matrix)::Alpha,Beta
      End Type MQC_Determinant_String
!
!     Parent Type
!>    \brief <b> Binary strings reresenting occupation number vectors</b>
      Type MQC_Determinant
        Type(MQC_Determinant_String)::Strings
        Character(Len=64)::Order !lexical,ci 
        Integer(kind=int64)::NDets,NAlpStr,NBetStr
        Type(MQC_Vector)::NSubsAlpha,NSubsBeta
      End Type MQC_Determinant
!
!    
!     Two electron integrals....
!
!     Parent Type
!>    \brief <b> Two-electron resonance integrals</b>
      type mqc_twoERIs
        type(mqc_r4tensor),private::alpha,beta,alphaBeta,betaAlpha,aaab,aaba,abaa,baaa,abab, &
          baab,baba,abba,bbba,bbab,babb,abbb
        character(Len=64),private::integralType !regular,raffenetti1,raffenetti2,raffenetti3,space,spin,general
        character(Len=64),private::storageType !full,symm
      contains 
!>      \brief <b> Print the MQC_TwoERIs object</b>
        procedure, public::print => mqc_print_twoERIs
        procedure, public::at => mqc_twoERIs_at
!>      \brief <b> Return the spin block dimention of the MQC_TwoERIs object</b>
        procedure, public::blockSize => mqc_eris_dimension
        procedure, private::hasSpinBlock => mqc_eris_has_spinBlock
!>      \brief <b> Return a spin block of MQC_TwoERIs object as an MQC rank-4 tensor</b>
        procedure, public:: getBlock => mqc_2eris_output_block
        procedure, public::type => mqc_eris_array_type
      end type mqc_twoERIs
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!>    \brief <b> Prints an object</b>
      interface mqc_print
        module procedure mqc_print_wavefunction
        module procedure mqc_print_integral
        module procedure mqc_print_eigenvalues
        module procedure mqc_print_twoERIS
      end interface
!
!>    \brief <b> Multiplies two objects</b>
      Interface MatMul
        Module Procedure MQC_Integral_Matrix_Multiply
        Module Procedure MQC_Matrix_Integral_Multiply
        Module Procedure MQC_Integral_Integral_Multiply
        Module Procedure MQC_Integral_Eigenvalues_Multiply
        Module Procedure MQC_Eigenvalues_Integral_Multiply
        Module Procedure MQC_Eigenvalues_Eigenvalues_Multiply
      End Interface
!
!>    \brief <b> Returns the dot product</b>
      Interface Dot_Product
        Module Procedure MQC_Eigenvalue_Eigenvalue_DotProduct
      End Interface
!
!>    \brief <b> Returns the transpose</b>
      Interface Transpose
        Module Procedure MQC_Integral_Transpose
      End Interface
!
!>    \brief <b> Returns the Hermitian conjugate</b>
      Interface Dagger
        Module Procedure MQC_Integral_Conjugate_Transpose
      End Interface
!
!>    \brief <b> Contracts two objects</b>
      Interface Contraction
        Module Procedure mqc_scf_integral_contraction
        Module Procedure mqc_eri_integral_contraction
        Module Procedure mqc_eri_r4tensor_contraction
      End Interface
!
!>    \brief <b> Converts object from alpha-beta order to energy order</b>
      Interface MQC_Matrix_UndoSpinBlockGHF
        Module Procedure MQC_Matrix_UndoSpinBlockGHF_Eigenvalues
        Module Procedure MQC_Matrix_UndoSpinBlockGHF_Integral
      End Interface
!
!>    \brief <b> Returns the complex conjugate</b>
      Interface Conjg
        Module Procedure mqc_scf_integral_conjg
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
!>    \brief <b> Assigns a variable to the value of another</b>
      Interface Assignment (=)
        Module Procedure MQC_Integral_Output_Array
        Module Procedure MQC_Eigenvalues_Output_Array
        Module Procedure mqc_2ERIs_output_array
      End Interface
!
!>    \brief <b> Sums two variables</b>
      Interface Operator (+)
        Module Procedure MQC_Integral_Sum
      End Interface
!
!>    \brief <b> Subtracts two variables</b>
      Interface Operator (-)
        Module Procedure MQC_Integral_Difference
        Module Procedure MQC_Matrix_Integral_Difference
      End Interface
!
!>    \brief <b> Multiplies two variables</b>
      Interface Operator (*)
        Module Procedure MQC_Scalar_Integral_Multiply
        Module Procedure MQC_Integral_Scalar_Multiply
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
!     PROCEDURE MQC_Print_Wavefunction
!
!>    \brief <b> MQC_Print_Wavefunction is a subroutine used to print an MQC
!>    EST wavefunction object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_Wavefunction is a subroutine used to print an MQC EST wavefunction 
!>    object. A optional argument can be added to specify the array that is printed
!>    from the wavefunction object. The following options are available:
!>
!>    1.  Label = 'overlap' prints the overlap matrix.
!>    2.  Label = 'core hamiltonian' prints the core hamiltonian matrix.
!>    3.  Label = 'orbital energies' prints the orbital energies vector.
!>    4.  Label = 'mo coefficients' prints the MO coefficient matrix.
!>    5.  Label = 'density' prints the final density matrix (SCF, MP2 etc.).
!>    6.  Label = 'scf density' prints the SCF density matrix.
!>    7.  Label = 'fock' prints the Fock matrix.
!>    8.  Label = 'nbasis' prints the number of basis functions.
!>    9.  Label = 'nalpha' prints the number of alpha electrons
!>    10. Label = 'nbeta' prints the number of beta electrons.
!>    11. Label = 'nelectrons' prints the total number of electrons.
!>    12. Label = 'charge' prints the molecular charge.
!>    13. Label = 'multiplicity' prints the electron multiplicity.
!>    14. Label = 'type' prints the wavefunction spin type (R,U,G).
!>    15. Label = 'complex' prints TRUE if the wavefunction is complex and FALSE if 
!>                it is real.
!>    16. Label = 'all' prints all EST wavefunction objects. This is the default.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Wavefunction
!>    \verbatim
!>        Wavefunction is Class(MQC_Wavefunction)
!>        The EST wavefunction variable to be printed.
!>    \endverbatim
!>
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64) 
!>        The Fortran file number to print to. 
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*)
!>        The specific object to print.
!>        = 'overlap':          overlap matrix.
!>        = 'core hamiltonian': core hamiltonian matrix.
!>        = 'orbital energies': orbital energies vector.
!>        = 'mo coefficients':  MO coefficient matrix.
!>        = 'density':          final density matrix.
!>        = 'scf density':      SCF density matrix.
!>        = 'fock':             Fock matrix.
!>        = 'nbasis':           number of basis functions.
!>        = 'nalpha':           number of alpha electrons.
!>        = 'nbeta':            number of beta electrons.
!>        = 'nelectrons':       total number of electrons.
!>        = 'charge':           molecular charge.
!>        = 'multiplicity':     electron multiplicity.
!>        = 'type':             wavefunction spin type (R,U,G).
!>        = 'complex':          wavefunction complex (T,F).
!>        = 'all':              all objects (default).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      subroutine mqc_print_wavefunction(wavefunction,iOut,label)
!
      implicit none
      class(mqc_wavefunction)::wavefunction
      integer(kind=int64),intent(in)::iOut
      character(len=*),optional,intent(in)::label
      character(len=64)::arrayType,myLabel
!
 1000 Format(1x,A)
 1050 Format( 2A )
 1060 Format( A,L10 )
!  
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = 'all'
      endIf
      
      select case (myLabel) 
      case('overlap')
        call wavefunction%overlap_matrix%print(iOut,'Overlap matrix',.true.,.true.)
      case('core hamiltonian')
        call wavefunction%core_hamiltonian%print(iOut,'Core Hamiltonian matrix',.true.,.true.)
      case('orbital energies')
        call wavefunction%mo_energies%print(iOut,'Orbital Energies',.true.,.true.)
      case('mo coefficients')
        call wavefunction%mo_coefficients%print(iOut,'Molecular orbital coefficients',.true.,.true.)
      case('density')
        call wavefunction%density_matrix%print(iOut,'Density matrix',.true.,.true.)
      case('scf density')
        call wavefunction%scf_density_matrix%print(iOut,'SCF Density matrix',.true.,.true.)
      case('fock')
        call wavefunction%fock_matrix%print(iOut,'Fock matrix',.true.,.true.)
      case('nbasis')
        call wavefunction%nBasis%print(iOut,'nBasis',.true.,.true.)
      case('nalpha')
        call wavefunction%nAlpha%print(iOut,'nAlpha',.true.,.true.)
      case('nbeta')
        call wavefunction%nBeta%print(iOut,'nBeta',.true.,.true.)
      case('nelectrons')
        call wavefunction%nElectrons%print(iOut,'nElectrons',.true.,.true.)
      case('charge')
        call wavefunction%charge%print(iOut,'Charge',.true.,.true.)
      case('multiplicity')
        call wavefunction%multiplicity%print(iOut,'Multiplicity',.true.,.true.)
      case('type')
        write(iOut,1050)' Wavefunction Type = ',Wavefunction%wF_type
      case('complex')
        write(iOut,1060)' Complex   = ',Wavefunction%wF_complex
      case('all')
        call wavefunction%nBasis%print(iOut,'nBasis',.true.,.false.)
        call wavefunction%nAlpha%print(iOut,'nAlpha',.true.,.false.)
        call wavefunction%nBeta%print(iOut,'nBeta',.true.,.false.)
        call wavefunction%nElectrons%print(iOut,'nElectrons',.true.,.false.)
        call wavefunction%charge%print(iOut,'Charge',.true.,.false.)
        call wavefunction%multiplicity%print(iOut,'Multiplicity',.true.,.true.)
        write(iOut,1050)'Wavefunction Type = ',Wavefunction%wF_type
        write(iOut,1060)'Complex   = ',Wavefunction%wF_complex

        if(mqc_integral_isAllocated(wavefunction%overlap_matrix)) &
          call wavefunction%overlap_matrix%print(iOut,'Overlap matrix',.true.,.false.)
        if(mqc_integral_isAllocated(wavefunction%core_hamiltonian)) &
          call wavefunction%core_hamiltonian%print(iOut,'Core Hamiltonian matrix',.true.,.false.)
        if(mqc_eigenvalues_isAllocated(wavefunction%mo_energies)) &
          call wavefunction%mo_energies%print(iOut,'Orbital Energies',.true.,.false.)
        if(mqc_integral_isAllocated(wavefunction%mo_coefficients)) &
          call wavefunction%mo_coefficients%print(iOut,'Molecular orbital coefficients',.true.,.false.)
        if(mqc_integral_isAllocated(wavefunction%density_matrix)) &
          call wavefunction%density_matrix%print(iOut,'Density matrix',.true.,.false.)
        if(mqc_integral_isAllocated(wavefunction%scf_density_matrix)) &
          call wavefunction%scf_density_matrix%print(iOut,'SCF Density matrix',.true.,.false.)
        if(mqc_integral_isAllocated(wavefunction%fock_matrix)) &
          call wavefunction%fock_matrix%print(iOut,'Fock matrix',.true.,.true.)
      case default
        call mqc_error_A('Invalid label sent to wavefunction print', 6, &
             'myLabel', myLabel )
      end select
!
      end subroutine mqc_print_wavefunction
!
!
!     PROCEDURE MQC_Print_Integral
!
!>    \brief <b> MQC_Print_Integral is a subroutine used to print an MQC
!>    EST SCF integral object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_Integral is a subroutine used to print an MQC EST SCF integral object. 
!>    Blank_At_Top and Blank_At_Bottom are optional logical arguments to print blank
!>    lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The EST SCF integral variable to be printed.
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
!>    \author L. M. Thompson
!>    \date 2016
!
      subroutine mqc_print_integral(integral,iOut,header, &
          blank_at_top,blank_at_bottom,formatStr)
!
      implicit none
      class(mqc_scf_integral)::integral
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      character(len=64)::arrayType
      logical,intent(In),optional::blank_at_top,blank_at_bottom
      character(len=*),intent(in),optional::formatStr
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
      if(present(blank_at_top)) then
        if(blank_at_top) write(iout,1020)
      endif
      write(iout,1000) trim(header)
      
      call string_change_case(integral%array_type,'L')
      if(integral%array_type.eq.'space') then
        if(present(formatStr)) then
          if(integral%hasAlpha()) call mqc_print(integral%alpha,iout,'',formatStr=formatStr)
        else
          if(integral%hasAlpha()) call mqc_print(integral%alpha,iout,'')
        endIf
      elseif(integral%array_type.eq.'spin') then
        if(present(formatStr)) then
          if(integral%hasAlpha()) call integral%alpha%print(iout,'Alpha Array',formatStr=formatStr)
          if(integral%hasBeta()) call integral%beta%print(iout,'Beta Array',formatStr=formatStr)
        else
          if(integral%hasAlpha()) call integral%alpha%print(iout,'Alpha Array')
          if(integral%hasBeta()) call integral%beta%print(iout,'Beta Array')
        endIf
      elseif(integral%array_type.eq.'general') then
        if(present(formatStr)) then
          if(integral%hasAlpha()) call integral%alpha%print(iout,'Alpha-Alpha Block',formatStr=formatStr)
          if(integral%hasBeta()) call integral%beta%print(iout,'Beta-Beta Block',formatStr=formatStr)
          if(integral%hasAlphaBeta()) call integral%alphabeta%print(iout,'Alpha-Beta Block',formatStr=formatStr)
          if(integral%hasBetaAlpha()) call integral%betaalpha%print(iout,'Beta-Alpha Block',formatStr=formatStr)
        else
          if(integral%hasAlpha()) call integral%alpha%print(iout,'Alpha-Alpha Block')
          if(integral%hasBeta()) call integral%beta%print(iout,'Beta-Beta Block')
          if(integral%hasAlphaBeta()) call integral%alphabeta%print(iout,'Alpha-Beta Block')
          if(integral%hasBetaAlpha()) call integral%betaalpha%print(iout,'Beta-Alpha Block')
        endIf
      else
        call mqc_error_A('Array type unrecogised in mqc_print_integral', 6, &
             'integral%array_type', integral%array_type )
      endIf

      if(present(blank_at_bottom)) then
        if(blank_at_bottom) write(iout,1020)
      endif
!
      end subroutine mqc_print_integral
!
!
!     PROCEDURE MQC_Print_Eigenvalues
!
!>    \brief <b> MQC_Print_Eigenvalues is a subroutine used to print an MQC EST
!>    eigenvalues object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_Eigenvalues is a subroutine used to print an MQC EST eigenvalues 
!>    object. Blank_At_Top and Blank_At_Bottom are optional logical arguments to 
!>    print blank lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_Eigenvalues)
!>        The EST eigenvalues variable to be printed.
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
!>        The title to print along with Eigenvalues.
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
!>    \author L. M. Thompson
!>    \date 2016
!
      subroutine mqc_print_eigenvalues(eigenvalues,iOut,header, &
          blank_at_top,blank_at_bottom,formatstr)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      character(len=64)::integralType
      logical,intent(In),optional::blank_at_top,blank_at_bottom
      character(len=*),intent(in),optional::formatStr
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
      if(present(blank_at_top)) then
        if(blank_at_top) write(iout,1020)
      endif
      write(iout,1000) trim(header)
      
      call string_change_case(eigenvalues%array_type,'L')
      if(eigenvalues%array_type.eq.'space') then
        if(present(formatStr)) then
          if(eigenvalues%hasAlpha()) call eigenvalues%alpha%print(iout,'',formatStr=formatStr)
        else
          if(eigenvalues%hasAlpha()) call eigenvalues%alpha%print(iout,'')
        endIf
      elseif(eigenvalues%array_type.eq.'spin'.or.eigenvalues%array_type.eq.'general') then
        if(present(formatStr)) then
          if(eigenvalues%hasAlpha()) call eigenvalues%alpha%print(iout,'Alpha Array',formatStr=formatStr)
          if(eigenvalues%hasBeta()) call eigenvalues%beta%print(iout,'Beta Array',formatStr=formatStr)
        else
          if(eigenvalues%hasAlpha()) call eigenvalues%alpha%print(iout,'Alpha Array')
          if(eigenvalues%hasBeta()) call eigenvalues%beta%print(iout,'Beta Array')
        endIf
      else
        call mqc_error_A('Array type unrecogised in mqc_print_eigenvalues', 6, &
             'eigenvalues%array_type', eigenvalues%array_type )
      endIf

      if(present(blank_at_bottom)) then
        if(blank_at_bottom) write(iout,1020)
      endif
!
      end subroutine mqc_print_eigenvalues
!
!
!     PROCEDURE MQC_Print_TwoERIs 
!
!>    \brief <b> MQC_Print_TwoERIs is a subroutine used to print an MQC EST 2ERIs
!>    object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_TwoERIs is a subroutine used to print an MQC EST 2ERIs object.
!>    Blank_At_Top and Blank_At_Bottom are optional logical arguments to print blank
!>    lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] TwoERIs
!>    \verbatim
!>        TwoERIs is Class(MQC_TwoERIs)
!>        The EST 2ERIs variable to be printed.
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
!>        The title to print along with TwoERIs.
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
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      subroutine mqc_print_twoERIs(twoERIs,iOut,header, &
          blank_at_top,blank_at_bottom)
!
      implicit none
      class(mqc_twoERIs)::twoERIs
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      character(len=64)::integralType
      logical,intent(In),optional::blank_at_top,blank_at_bottom
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
      if(present(blank_at_top)) then
        if(blank_at_top) write(iout,1020)
      endif
      write(iout,1000) trim(header)
      
      call string_change_case(twoERIs%integraltype,'L')
      if(twoERIs%integralType.eq.'regular') then
        call twoERIs%alpha%print(iout,'Regular 2ERIs')
      elseIf(twoERIs%integralType.eq.'raffenetti1') then
        call twoERIs%alpha%print(iout,'Raffenetti 1 2ERIs')
      elseIf(twoERIs%integralType.eq.'raffenetti2') then
        call twoERIs%alpha%print(iout,'Raffenetti 2 2ERIs')
      elseIf(twoERIs%integralType.eq.'raffenetti3') then
        call twoERIs%alpha%print(iout,'Raffenetti 3 2ERIs')
      elseIf(twoERIs%integralType.eq.'space') then
        call twoERIs%alpha%print(iout,'Restricted MO 2ERIs')
      elseIf(twoERIs%integralType.eq.'spin') then
        write(iout,'(A)') 'Unrestricted MO 2ERIs'
        call twoERIs%alpha%print(iout,'(aa|aa)')
        call twoERIs%beta%print(iout,'(bb|bb)')
        call twoERIs%alphaBeta%print(iout,'(aa|bb)')
        if(twoERIs%storageType.eq.'full') call twoERIs%betaAlpha%print(iout,'(bb|aa)')
      elseIf(twoERIs%integralType.eq.'general') then
        write(iout,'(A)') 'General MO 2ERIs'
        call twoERIs%alpha%print(iout,'(aa|aa)')
        call twoERIs%beta%print(iout,'(bb|bb)')
        call twoERIs%alphaBeta%print(iout,'(aa|bb)')
        if(twoERIs%storageType.eq.'full') call twoERIs%betaAlpha%print(iout,'(bb|aa)')
        call twoERIs%abab%print(iout,'(ab|ab)')
        call twoERIs%abba%print(iout,'(ab|ba)')
        if(twoERIs%storageType.eq.'full') call twoERIs%baab%print(iout,'(ba|ab)')
        if(twoERIs%storageType.eq.'full') call twoERIs%baba%print(iout,'(ba|ba)')
        call twoERIs%aaab%print(iout,'(aa|ab)')
        if(twoERIs%storageType.eq.'full') call twoERIs%aaba%print(iout,'(aa|ba)')
        if(twoERIs%storageType.eq.'full') call twoERIs%abaa%print(iout,'(ab|aa)')
        if(twoERIs%storageType.eq.'full') call twoERIs%baaa%print(iout,'(ba|aa)')
        if(twoERIs%storageType.eq.'full') call twoERIs%bbba%print(iout,'(bb|ba)')
        if(twoERIs%storageType.eq.'full') call twoERIs%bbab%print(iout,'(bb|ab)')
        if(twoERIs%storageType.eq.'full') call twoERIs%babb%print(iout,'(ba|bb)')
        call twoERIs%abbb%print(iout,'(ab|bb)')
      else
        call mqc_error_A('Integral type unrecogised in mqc_print_twoERIs', 6, &
             'twoERIs%integralType', twoERIs%integralType )
      endIf

      if(present(blank_at_bottom)) then
        if(blank_at_bottom) write(iout,1020)
      endif
!
      end subroutine mqc_print_twoERIs 
!
!
!     PROCEDURE MQC_Integral_isAllocated
!
!>    \brief <b> MQC_Integeral_isAllocated is a function that returns TRUE if an MQC
!>    EST integral is allocated and FALSE if it is not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_isAllocated is a function that returns TRUE if an MQC EST integral 
!>    is allocated and FALSE if it is not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      Function MQC_Integral_isAllocated(Integral) Result(isAllocated)
!
      Implicit None
      Class(MQC_SCF_Integral),Intent(InOut)::Integral
      Logical::isAllocated
!
      isAllocated = .False.
      If(MQC_Matrix_isAllocated(Integral%Alpha)) isAllocated = .True.
      If(MQC_Matrix_isAllocated(Integral%Beta)) isAllocated = .True.
      If(MQC_Matrix_isAllocated(Integral%AlphaBeta)) isAllocated = .True.
      If(MQC_Matrix_isAllocated(Integral%BetaAlpha)) isAllocated = .True.
!
      Return
      End Function MQC_Integral_isAllocated
!
!
!     PROCEDURE MQC_Eigenvalues_isAllocated
!
!>    \brief <b> MQC_Eigenvalues_isAllocated is a function that returns TRUE if an 
!>    MQC EST eigenvalues variable is allocated and FALSE if it is not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_isAllocated is a function that returns TRUE if an MQC 
!>    eigenvalues variable is allocated and FALSE if it is not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_Eigenvalues)
!>        The name of the MQC_Eigenvalues variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      Function MQC_Eigenvalues_isAllocated(Eigenvalues) Result(isAllocated)
!
      Implicit None
      Class(MQC_SCF_Eigenvalues),Intent(InOut)::Eigenvalues
      Logical::isAllocated
!
      isAllocated = .False.
      If(MQC_Vector_isAllocated(Eigenvalues%Alpha)) isAllocated = .True.
      If(MQC_Vector_isAllocated(Eigenvalues%Beta)) isAllocated = .True.
!
      Return
      End Function MQC_Eigenvalues_isAllocated
!
!
!     PROCEDURE MQC_Integral_Has_Alpha
!
!>    \brief <b> MQC_Integral_Has_Alpha is a function that returns TRUE if an MQC EST 
!>    SCF integral variable has an allocated alpha spin block and FALSE if it does not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Has_Alpha is a function that returns TRUE if an MQC EST SCF 
!>    integral variable has an allocated alpha spin block and FALSE if it does not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_has_alpha(integral) result(hasAlpha)
!
      implicit none
      class(mqc_scf_integral)::integral
      logical::hasAlpha
!
      hasAlpha = .false.
      if(MQC_Matrix_isAllocated(integral%alpha)) hasAlpha = .true. 
!
      end function mqc_integral_has_alpha
!
!
!     PROCEDURE MQC_Integral_Has_Beta 
!
!>    \brief <b> MQC_Integral_Has_Beta is a function that returns TRUE if an MQC EST 
!>    SCF integral variable has an allocated beta spin block and FALSE if it does not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Has_beta is a function that returns TRUE if an MQC EST SCF 
!>    integral variable has an allocated beta spin block and FALSE if it does not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_has_beta(integral) result(hasBeta)
!
      implicit none
      class(mqc_scf_integral)::integral
      logical::hasBeta
!
      hasBeta = .false.
      if(MQC_Matrix_isAllocated(integral%beta)) hasBeta = .true. 
!
      end function mqc_integral_has_beta
!
!
!     PROCEDURE MQC_Integral_Has_AlphaBeta
!
!>    \brief <b> MQC_Integral_Has_AlphaBeta is a function that returns TRUE if an MQC 
!>    EST SCF integral variable has an allocated alpha-beta spin block and FALSE if 
!>    it does not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Has_AlphaBeta is a function that returns TRUE if an MQC EST SCF 
!>    integral variable has an allocated alpha-beta spin block and FALSE if it does 
!>    not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_has_alphaBeta(integral) result(hasAlphaBeta)
!
      implicit none
      class(mqc_scf_integral)::integral
      logical::hasAlphaBeta
!
      hasAlphaBeta = .false.
      if(MQC_Matrix_isAllocated(integral%alphaBeta)) hasAlphaBeta = .true. 
!
      end function mqc_integral_has_alphaBeta
!
!
!     PROCEDURE MQC_Integral_Has_BetaAlpha 
!
!>    \brief <b> MQC_Integral_Has_BetaAlpha is a function that returns TRUE if an MQC 
!>    EST SCF integral variable has an allocated beta-alpha spin block and FALSE if 
!>    it does not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Has_BetaAlpha is a function that returns TRUE if an MQC EST SCF 
!>    integral variable has an allocated beta-alpha spin block and FALSE if it does 
!>    not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_has_betaAlpha(integral) result(hasBetaAlpha)
!
      implicit none
      class(mqc_scf_integral)::integral
      logical::hasBetaAlpha
!
      hasBetaAlpha = .false.
      if(MQC_Matrix_isAllocated(integral%betaAlpha)) hasBetaAlpha = .true. 
!
      end function mqc_integral_has_betaAlpha
!
!
!     PROCEDURE MQC_Eigenvalues_Has_Alpha
!
!>    \brief <b> MQC_Eigenvalues_Has_Alpha is a function that returns TRUE if an MQC 
!>    EST eigenvalues variable has an allocated alpha spin block and FALSE if it does 
!>    not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Has_Alpha is a function that returns TRUE if an MQC EST 
!>    eigenvalues integral variable has an allocated alpha spin block and FALSE if it 
!>    does not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_Eigenvalues)
!>        The name of the MQC_Eigenvalues variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_eigenvalues_has_alpha(eigenvalues) result(hasAlpha)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      logical::hasAlpha
!
      hasAlpha = .false.
      if(MQC_Vector_isAllocated(eigenvalues%alpha)) hasAlpha = .true. 
!
      end function mqc_eigenvalues_has_alpha
!
!
!     PROCEDURE MQC_Eigenvalues_Has_Beta 
!
!>    \brief <b> MQC_Eigenvalues_Has_Beta is a function that returns TRUE if an MQC 
!>    EST eigenvalues variable has an allocated beta spin block and FALSE if it does 
!>    not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Has_Beta is a function that returns TRUE if an MQC EST 
!>    eigenvalues integral variable has an allocated beta spin block and FALSE if it 
!>    does not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_Eigenvalues)
!>        The name of the MQC_Eigenvalues variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_eigenvalues_has_beta(eigenvalues) result(hasBeta)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      logical::hasBeta
!
      hasBeta = .false.
      if(MQC_Vector_isAllocated(eigenvalues%beta)) hasBeta = .true. 
!
      end function mqc_eigenvalues_has_beta
!
!
!     PROCEDURE MQC_ERIs_Has_SpinBlock
!
!>    \brief <b> MQC_ERIs_Has_SpinBlock is a function that returns TRUE if an MQC 2ERIs 
!>    variable has the specified spin block allocated and FALSE if it does not</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_ERIs_Has_SpinBlock is a function that returns TRUE if an MQC 2ERIs variable 
!>    has the specified spin block allocated and FALSE if it does not.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] ERIs
!>    \verbatim
!>        ERIs is Class(MQC_TwoERIs)
!>        The name of the MQC_TwoERIs variable.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is character(len=*)
!>        = 'alpha':       tests the alpha-alpha spin block.
!>        = 'alpha-alpha': tests the alpha-alpha spin block. 
!>        = 'aaaa':        tests the alpha-alpha spin block. 
!>        = 'beta':        tests the beta-beta spin block. 
!>        = 'beta-beta':   tests the beta-beta spin block. 
!>        = 'bbbb':        tests the beta-beta spin block. 
!>        = 'alpha-beta'   tests the alpha-beta spin block. 
!>        = 'aabb':        tests the alpha-beta spin block. 
!>        = 'beta-alpha'   tests the beta-alpha spin block. 
!>        = 'bbaa':        tests the beta-alpha spin block. 
!>        = 'abab':        tests the (ab|ab) spin block. 
!>        = 'abba':        tests the (ab|ba) spin block. 
!>        = 'baab':        tests the (ba|ab) spin block. 
!>        = 'baba':        tests the (ba|ba) spin block. 
!>        = 'aaab':        tests the (aa|ab) spin block. 
!>        = 'aaba':        tests the (aa|ba) spin block. 
!>        = 'abaa':        tests the (ab|aa) spin block. 
!>        = 'baaa':        tests the (ba|aa) spin block. 
!>        = 'bbba':        tests the (bb|ba) spin block. 
!>        = 'bbab':        tests the (bb|ab) spin block. 
!>        = 'babb':        tests the (ba|bb) spin block. 
!>        = 'abbb':        tests the (ab|bb) spin block. 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2019
!
      function mqc_eris_has_spinBlock(eris,label) result(hasBlock)
!
      implicit none
      class(mqc_twoeris),intent(in)::eris
      character(len=*),intent(in)::label
      character(len=12)::mylabel
      logical::hasBlock
!
      call string_change_case(label,'L',myLabel)
      
      hasBlock = .false.
      if(mylabel.eq.'aaaa'.or.mylabel.eq.'alpha'.or.mylabel.eq.'alpha-alpha') then
        if(MQC_R4Tensor_isAllocated(eris%alpha)) hasBlock = .true. 
      elseIf(mylabel.eq.'bbbb'.or.mylabel.eq.'beta'.or.mylabel.eq.'beta-beta') then
        if(MQC_R4Tensor_isAllocated(eris%beta)) hasBlock = .true. 
      elseIf(mylabel.eq.'aabb'.or.mylabel.eq.'alpha-beta') then
        if(MQC_R4Tensor_isAllocated(eris%alphaBeta)) hasBlock = .true. 
      elseIf(mylabel.eq.'bbaa'.or.mylabel.eq.'beta-alpha') then
        if(MQC_R4Tensor_isAllocated(eris%betaAlpha)) hasBlock = .true. 
      elseIf(mylabel.eq.'aaab') then
        if(MQC_R4Tensor_isAllocated(eris%aaab)) hasBlock = .true. 
      elseIf(mylabel.eq.'aaba') then
        if(MQC_R4Tensor_isAllocated(eris%aaba)) hasBlock = .true. 
      elseIf(mylabel.eq.'abaa') then
        if(MQC_R4Tensor_isAllocated(eris%abaa)) hasBlock = .true. 
      elseIf(mylabel.eq.'baaa') then
        if(MQC_R4Tensor_isAllocated(eris%baaa)) hasBlock = .true. 
      elseIf(mylabel.eq.'abab') then
        if(MQC_R4Tensor_isAllocated(eris%abab)) hasBlock = .true. 
      elseIf(mylabel.eq.'baab') then
        if(MQC_R4Tensor_isAllocated(eris%baab)) hasBlock = .true. 
      elseIf(mylabel.eq.'baba') then
        if(MQC_R4Tensor_isAllocated(eris%baba)) hasBlock = .true. 
      elseIf(mylabel.eq.'abba') then
        if(MQC_R4Tensor_isAllocated(eris%abba)) hasBlock = .true. 
      elseIf(mylabel.eq.'bbba') then
        if(MQC_R4Tensor_isAllocated(eris%bbba)) hasBlock = .true. 
      elseIf(mylabel.eq.'bbab') then
        if(MQC_R4Tensor_isAllocated(eris%bbab)) hasBlock = .true. 
      elseIf(mylabel.eq.'babb') then
        if(MQC_R4Tensor_isAllocated(eris%babb)) hasBlock = .true. 
      elseIf(mylabel.eq.'abbb') then
        if(MQC_R4Tensor_isAllocated(eris%abbb)) hasBlock = .true. 
      else
        call mqc_error_a('Unrecognized block requested in mqc_eris_has_spinblock',6,'myLabel',myLabel)
      endIf
!
      end function mqc_eris_has_spinBlock
!
!
!     PROCEDURE MQC_Integral_Array_Type
!
!>    \brief <b> MQC_Integral_Array_Type is a function that returns whether an MQC 
!>    EST SCF integral is of restricted (R), unrestricted (U) or general (G) type</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Array_Type is a function that returns whether an MQC EST SCF 
!>    integral is of restricted (R), unrestricted (U) or general (G) type.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_array_type(integral) result(arrayType)
!
      implicit none
      class(mqc_scf_integral)::integral
      Character(Len=64)::arrayType
!
      arrayType = integral%array_type 
!
      end function mqc_integral_array_type
!
!
!     PROCEDURE MQC_Eigenvalues_Array_Type
!
!>    \brief <b> MQC_Eigenvalues_Array_Type is a function that returns whether an MQC 
!>    EST eigenvalues is of restricted (R), unrestricted (U) or general (G) type</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Array_Type is a function that returns whether an MQC EST 
!>    eigenvalues is of restricted (R), unrestricted (U) or general (G) type.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_Eigenvalues)
!>        The name of the MQC_Eigenvalues variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_eigenvalues_array_type(eigenvalues) result(arrayType)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      Character(Len=64)::arrayType
!
      arrayType = eigenvalues%array_type 
!
      end function mqc_eigenvalues_array_type
!
!
!     PROCEDURE MQC_ERIs_Array_Type
!
!>    \brief <b> MQC_ERIs_Array_Type is a function that returns whether an MQC 2ERIs
!>    variable is of AO (regular, raffenettit1, raffenetti2, raffenetti3) or MO 
!>    (restricted, unrestricted, general) type</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    \brief <b> MQC_ERIs_Array_Type is a function that returns whether an MQC 2ERIs
!>    variable is of AO (regular, raffenettit1, raffenetti2, raffenetti3) or MO 
!>    (restricted, unrestricted, general) type</b>
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] ERIs
!>    \verbatim
!>        Integral is Class(MQC_TwoERIs)
!>        The name of the MQC_TwoERIs variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2019
!
      function mqc_eris_array_type(eris) result(arrayType)
!
      implicit none
      class(mqc_twoeris)::eris
      Character(Len=64)::arrayType
!
      arrayType = eris%integraltype 
!
      end function mqc_eris_array_type
!
!
!     PROCEDURE MQC_Integral_Array_Name    
!
!>    \brief <b> MQC_Integral_Array_Name is a function that returns the label of an 
!>    MQC EST SCF integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Array_Name is a function that returns the label of an MQC EST SCF 
!>    integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_array_name(integral) result(arrayName)
!
      implicit none
      class(mqc_scf_integral)::integral
      Character(Len=64)::arrayName
!
      arrayName = integral%array_name 
!
      end function mqc_integral_array_name
!
!
!     PROCEDURE MQC_Eigenvalues_Array_Name    
!
!>    \brief <b> MQC_Eigenvalues_Array_Name is a function that returns the label of 
!>    an MQC EST eigenvalues</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Array_Name is a function that returns the label of an MQC EST 
!>    eigenvalues.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_Eigenvalues)
!>        The name of the MQC_Eigenvalues variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_eigenvalues_array_name(eigenvalues) result(arrayName)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      Character(Len=64)::arrayName
!
      arrayName = eigenvalues%array_name 
!
      end function mqc_eigenvalues_array_name
!
!
!     PROCEDURE MQC_Integral_Add_Name    
!
!>    \brief <b> MQC_Integral_Add_Name is a function that sets the label of an MQC 
!>    EST SCF integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Add_Name is a function that sets the label of an MQC EST SCF 
!>    integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!>
!>    \param[in] ArrayName
!>    \verbatim
!>        ArrayName is Character(Len=*)
!>        The label of the MQC_SCF_Integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      subroutine mqc_integral_add_name(integral,arrayName) 
!
      implicit none
      class(mqc_scf_integral)::integral
      character(Len=*)::arrayName
!
      integral%array_name = arrayName
!
      end subroutine mqc_integral_add_name
!
!
!     PROCEDURE MQC_eigenvalues_Add_Name    
!
!>    \brief <b> MQC_Eigenvalues_Add_Name is a function that sets the label of an MQC 
!>    EST eigenvalues variable</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Add_Name is a function that sets the label of an MQC EST SCF 
!>    eigenvalues variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_Eigenvalues)
!>        The name of the MQC_SCF_Eigenvalues variable.
!>    \endverbatim
!>
!>    \param[in] ArrayName
!>    \verbatim
!>        ArrayName is Character(Len=*)
!>        The label of the MQC_SCF_Eigenvalues variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      subroutine mqc_eigenvalues_add_name(eigenvalues,arrayName) 
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      Character(Len=*)::arrayName
!
      eigenvalues%array_name = arrayName
!
      end subroutine mqc_eigenvalues_add_name
!
!
!     PROCEDURE MQC_Integral_Dimension
!
!>    \brief <b> MQC_Integral_Dimension is a function that returns the length of a
!>    spin block axis of an MQC_SCF_Integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Dimension is a function that returns the length of a spin block 
!>    of an MQC_SCF_Integral. The required argument label specifies which block is 
!>    queried. The optional argument axis specifies if rows or columns are returned.
!>    The following options are available:
!>
!>    1. Label = 'alpha' returns the alpha spin block.
!>    2. Label = 'beta' returns the beta spin block.
!>    3. Axis  = '1' returns the spin block rows (default).
!>    4. Axis  = '2' returns the spin block columns.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*)
!>        The spin block to query.
!>    \endverbatim
!>
!>    \param[in] Axis
!>    \verbatim
!>        Axis is Integer(kind=int64),Optional
!>        The spin block axis to return.
!>        = 1: number of rows (default).
!>        = 2: number of columns.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_dimension(integral,label,axis) result(dimBlock)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      integer(kind=int64)::dimBlock
      Character(Len=*),intent(in)::label
      integer(kind=int64),intent(in),optional::axis
      Character(Len=64)::myLabel
      integer(kind=int64)::my_axis
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf

      call string_change_case(label,'l',myLabel)
      select case (myLabel)
      case('alpha') 
        if(.not.integral%hasAlpha()) then
          dimBlock = 0
        elseIf(my_axis.eq.1) then
          dimBlock = mqc_matrix_rows(integral%alpha)
        elseIf(my_axis.eq.2) then
          dimBlock = mqc_matrix_columns(integral%alpha)
        else
          call mqc_error_I('Unknown axis specified in MQC_Integral_Dimension',6,'my_axis',my_axis)
        endIf
      case('beta')
        if(.not.integral%hasBeta()) then
          if(integral%type().eq.'space') then
            if(.not.integral%hasAlpha()) then
              dimBlock = 0
            elseIf(my_axis.eq.1) then
              dimBlock = mqc_matrix_rows(integral%alpha)
            elseIf(my_axis.eq.2) then
              dimBlock = mqc_matrix_columns(integral%alpha)
            else
              call mqc_error_I('Unknown axis specified in MQC_Integral_Dimension',6,'my_axis',my_axis)
            endIf
          else
            dimBlock = 0
          endIf
        elseIf(my_axis.eq.1) then
          dimBlock = mqc_matrix_rows(integral%beta)
        elseIf(my_axis.eq.2) then
          dimBlock = mqc_matrix_columns(integral%beta)
        else
          call mqc_error_I('Unknown axis specified in MQC_Integral_Dimension',6,'my_axis',my_axis)
        endIf
      case default
        call mqc_error_A('label not valid in mqc_integrals_dimension', 6, &
             'myLabel', myLabel )
      end select 
!
      end function mqc_integral_dimension 
!
!
!     PROCEDURE MQC_Eigenvalues_Dimension
!
!>    \brief <b> MQC_Eigenvalues_Dimension is a function that returns the length of a
!>    spin block axis of an MQC_Eigenvalues variable</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Dimension is a function that returns the length of a spin block 
!>    of an MQC_Eigenvalues variable. The required argument label specifies which block is 
!>    queried. The following options are available:
!>
!>    1. Label = 'alpha' returns the alpha spin block.
!>    2. Label = 'beta' returns the beta spin block.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_Eigenvalues)
!>        The name of the MQC_Eigenvalues variable.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*)
!>        The spin block to query.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_eigenvalues_dimension(eigenvalues,label) result(dimBlock)
!
      implicit none
      class(mqc_scf_eigenvalues),intent(in)::eigenvalues
      integer(kind=int64)::dimBlock
      Character(Len=*),intent(in)::label
      Character(Len=64)::myLabel
!
      call string_change_case(label,'l',myLabel)
      select case (myLabel)
      case('alpha')
        if(.not.eigenvalues%hasAlpha()) then
          dimBlock = 0
        else
          dimBlock = mqc_length_vector(eigenvalues%alpha)
        endIf
      case('beta')
        if(.not.eigenvalues%hasBeta()) then
          if(eigenvalues%type().eq.'space') then
            dimBlock = mqc_length_vector(eigenvalues%alpha)
          else
            dimBlock = 0
          endIf
        else
          dimBlock = mqc_length_vector(eigenvalues%beta)
        endIf
      case default
        call mqc_error_A('label not valid in mqc_eigenvalues_dimension', 6, &
             'myLabel', myLabel )
      end select 
!
      end function mqc_eigenvalues_dimension 
!
!
!     PROCEDURE MQC_ERIs_Dimension
!
!>    \brief <b> MQC_ERIs_Dimension is a function that returns the length of a spin
!>    block axis of an MQC_TwoERIs object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_ERIs_Dimension is a function that returns the length of a spin block of an
!>    MQC_ERIs_Integral. The required argument label specifies which block is queried. 
!>    The optional argument axis specifies which dimension is returned. The following 
!>    options are available:
!>
!>    1. Label = 'alpha' returns the alpha spin block dimension.
!>    2. Label = 'beta' returns the beta spin block dimension.
!>    3. Axis  = '1' returns the first spin block dimension (default).
!>    4. Axis  = '2' returns the second spin block dimension.
!>    5. Axis  = '3' returns the third spin block dimension.
!>    6. Axis  = '4' returns the fourth spin block dimension.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] ERIs
!>    \verbatim
!>        ERIs is Class(MQC_TwoERIs)
!>        The name of the MQC_TwoERIs variable.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*)
!>        The spin block to query.
!>        = 'alpha': alpha axis dimension.
!>        = 'beta':  beta axis dimension.
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
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2019
!
      function mqc_eris_dimension(eris,label,axis) result(dimBlock)
!
      implicit none
      class(mqc_twoERIs),intent(in)::eris
      integer(kind=int64)::dimBlock
      Character(Len=*),intent(in)::label
      integer(kind=int64),intent(in),optional::axis
      Character(Len=64)::myLabel
      integer(kind=int64)::my_axis
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf
      if(my_axis.le.0.or.my_axis.gt.4) &
        call mqc_error_I('Unknown axis specified in MQC_ERIs_Dimension',6,'my_axis',my_axis)

      call string_change_case(label,'l',myLabel)

      select case (myLabel)
      case('alpha') 
        if(.not.eris%hasSpinBlock('alpha')) then
          dimBlock = 0
        else
          dimBlock = size(eris%alpha,my_axis)
        endIf
      case('beta')
        if(.not.eris%hasSpinBlock('beta')) then
          if(eris%type().eq.'regular'.or.eris%type().eq.'raffenetti1'.or.eris%type().eq.'raffenetti2'.or.&
            eris%type().eq.'raffenetti3'.or.eris%type().eq.'space') then
            if(.not.eris%hasSpinBlock('alpha')) then
              dimBlock = 0
            else
              dimBlock = size(eris%alpha,my_axis)
            endIf
          else
            dimBlock = 0
          endIf
        else
          dimBlock = size(eris%beta,my_axis)
        endIf
      case default
        call mqc_error_A('label not valid in mqc_eris_dimension', 6, &
             'myLabel', myLabel )
      end select 
!
      end function mqc_eris_dimension 
!
!
!     PROCEDURE MQC_TwoERIs_Allocate
!
!>    \brief <b> MQC_TwoERIs_Allocate is used to allocate a two-electron integral 
!>    type variable of the MQC_twoERIs class</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_TwoERIs_Allocate is used to allocate a two-electron integral type variable 
!>    of the MQC_twoERIs class. The following options are available:
!>
!>    1. StorageType  = 'Full' declares the MQC_TwoERIs variable is unpacked.
!>    2. StorageType  = 'Symm' declares the MQC_TwoERIs variable is symmetric packed.
!>    3. IntegralType = 'Regular' declares the MQC_TwoERIs variable contains regular 
!>                      atomic orbital 2ERIs. 
!>    4. IntegralType = 'Raffenetti1' declares the MQC_TwoERIs variable contains 
!>                      raffenetti type 1 atomic orbital 2ERIs. 
!>    5. IntegralType = 'Raffenetti2 declares the MQC_TwoERIs variable contains 
!>                      raffenetti type 2 atomic orbital 2ERIs. 
!>    6. IntegralType = 'Raffenetti3' declares the MQC_TwoERIs variable contains 
!>                      raffenetti type 3 atomic orbital 2ERIs. 
!>    7. IntegralType = 'space' declares the MQC_TwoERIs variable contains restricted
!>                      molecular orbital 2ERIs.
!>    8. IntegralType = 'spin' declares the MQC_TwoERIs variable contains unrestricted
!>                      molecular orbital 2ERIs.
!>    9. IntegralType = 'general' declares the MQC_TwoERIs variable contains general
!>                      molecular orbital 2ERIs.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] TwoERIs
!>    \verbatim
!>        TwoERIs is class(mqc_twoERIs)
!>        The MQC_TwoERIs variable to allocate.
!>    \endverbatim
!>    
!>    \param[in,out] StorageType
!>    \verbatim
!>        StorageType is Character(len=*)
!>        = 'Full': the MQC_TwoERIs variable is upacked
!>        = 'Symm': if the MQC_TwoERIs variable is AO (regular, raffenetti) 
!>                  or space, then Symm=Full and the assignment is arbitrary,
!>                  whereas, if the MQC_TwoERIs variable spin or general, 
!>                  only unique spin blocks are stored (first given in 
!>                  equality): 
!> 
!>                    spin: [aa|aa]=[bb|bb]=[aa|bb]=[bb|aa]
!>                    space: [aa|aa], [bb|bb], [aa|bb]=[bb|aa]
!>                    general: [aa|aa], [bb|bb], [aa|bb]=[bb|aa], 
!>                      [abab]=[ba|ba]*, [ab|ba]=[ba|ab], 
!>                      [aa|ab]=[ab|aa]=[aa|ba]*=[ba|aa]*,
!>                      [ab|bb]=[bb|ab]=[bb|ba]*=[ba|bb]*
!>
!>    \endverbatim
!>    
!>    \param[in,out] IntegralType
!>    \verbatim
!>        IntegralType is Character(len=*)
!>        = 'Regular':     the MQC_TwoERIs variable contains 
!>                         regular atomic orbital 2ERIs.
!>                         R(i,j,k,l) = (ij|kl)
!>        = 'Raffenetti1': the MQC_TwoERIs variable contains 
!>                         raffenetti type 1 atomic orbital 2ERIs.
!>                         R(i,j,k,l) = (ij|kl) - 1/4[(ik|jl)+(il|jk)]
!>        = 'Raffenetti2': the MQC_TwoERIs variable contains 
!>                         raffenetti type 2 atomic orbital 2ERIs.
!>                         R(i,j,k,l) = (ik|jl) + (il|jk)  
!>        = 'Raffenetti3': the MQC_TwoERIs variable contains 
!>                         raffenetti type 3 atomic orbital 2ERIs.
!>                         R(i,j,k,l) = (ik|jl) - (il|jk)
!>        = 'space':       the MQC_TwoERIs variable contains 
!>                         restricted molecular orbital 2ERIs.
!>        = 'spin':        the MQC_TwoERIs variable contains 
!>                         unrestricted molecular orbital 2ERIs.
!>        = 'general':     the MQC_TwoERIs variable contains 
!>                         general molecular orbital 2ERIs.
!>    \endverbatim
!>
!>    \param[in] Alpha
!>    \verbatim
!>        Alpha is type(mqc_r4tensor),optional
!>        The alpha spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] Beta
!>    \verbatim
!>        Beta is type(mqc_r4tensor),optional
!>        The beta spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] AlphaBeta
!>    \verbatim
!>        AlphaBeta is type(mqc_r4tensor),optional
!>        The alpha-beta spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] BetaAlpha
!>    \verbatim
!>        AlphaBeta is type(mqc_r4tensor),optional
!>        The beta-alpha spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] ABAB 
!>    \verbatim
!>        ABAB is type(mqc_r4tensor),optional
!>        The abab spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] ABBA 
!>    \verbatim
!>        ABBA is type(mqc_r4tensor),optional
!>        The abba spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] BAAB 
!>    \verbatim
!>        BAAB is type(mqc_r4tensor),optional
!>        The baab spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] BABA 
!>    \verbatim
!>        BABA is type(mqc_r4tensor),optional
!>        The baba spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] AAAB 
!>    \verbatim
!>        AAAB is type(mqc_r4tensor),optional
!>        The aaab spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] AABA 
!>    \verbatim
!>        AABA is type(mqc_r4tensor),optional
!>        The aaba spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] ABAA 
!>    \verbatim
!>        ABAA is type(mqc_r4tensor),optional
!>        The abaa spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] BAAA 
!>    \verbatim
!>        BAAA is type(mqc_r4tensor),optional
!>        The baaa spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] BBBA 
!>    \verbatim
!>        BBBA is type(mqc_r4tensor),optional
!>        The bbba spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] BBAB 
!>    \verbatim
!>        BBAB is type(mqc_r4tensor),optional
!>        The bbab spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] BABB 
!>    \verbatim
!>        BABB is type(mqc_r4tensor),optional
!>        The babb spin block 2ERI integrals.
!>    \endverbatim
!>
!>    \param[in] ABBB 
!>    \verbatim
!>        ABBB is type(mqc_r4tensor),optional
!>        The abbb spin block 2ERI integrals.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      subroutine mqc_twoeris_allocate(twoERIs,storageType,integralType,alpha, &
          beta,alphaBeta,betaAlpha,abab,abba,baab,baba,aaab,aaba,abaa,baaa,bbba, &
          bbab,babb,abbb) 
!
      implicit none
      class(mqc_twoERIs)::twoERIs
      character(len=*)::storageType,integralType
      type(mqc_r4tensor),optional::alpha,beta,alphaBeta,betaAlpha,abab,abba,baab, &
        baba,aaab,aaba,abaa,baaa,bbba,bbab,babb,abbb
!
      call string_change_case(storageType,'L')
      call string_change_case(integralType,'L')
      twoERIs%storageType = storageType 
      twoERIs%integralType = integralType 
!
      if(present(alpha)) twoERIs%alpha = alpha
      if(present(beta))  twoERIs%beta = beta
      if(present(alphaBeta)) twoERIs%alphaBeta = alphaBeta
      if(present(betaAlpha)) twoERIs%betaAlpha = betaAlpha
      if(present(abab)) twoERIs%abab = abab
      if(present(abba)) twoERIs%abba = abba
      if(present(baab)) twoERIs%baab = baab
      if(present(baba)) twoERIs%baba = baba
      if(present(aaab)) twoERIs%aaab = aaab
      if(present(aaba)) twoERIs%aaba = aaba
      if(present(abaa)) twoERIs%abaa = abaa
      if(present(baaa)) twoERIs%baaa = baaa
      if(present(bbba)) twoERIs%bbba = bbba
      if(present(bbab)) twoERIs%bbab = bbab
      if(present(babb)) twoERIs%babb = babb
      if(present(abbb)) twoERIs%abbb = abbb
!
      end subroutine mqc_twoeris_allocate  
!
!
!     PROCEDURE MQC_Integral_Allocate
!
!>    \brief <b> MQC_Integral_Allocate is used to allocate a MQC SCF integral type
!>    variable of the MQC_SCF_Integral class</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Allocate is used to allocate a MQC SCF integral type variable of 
!>    the MQC_SCF_Integral class. The following options are available:
!>
!>    1. ArrayType = 'space' declares the MQC_SCF_Integral variable contains 
!>                   restricted spin blocks.
!>    2. ArrayType = 'spin' declares the MQC_SCF_Integral variable contains 
!>                   unrestricted spin blocks.
!>    3. ArrayType = 'general' declares the MQC_SCF_Integral variable contains 
!>                   general spin blocks.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is class(mqc_scf_integral)
!>        The MQC_SCF_Integral variable to allocate.
!>    \endverbatim
!>    
!>    \param[in,out] ArrayName
!>    \verbatim
!>        ArrayName is Character(len=*)
!>        The label of the MQC_SCF_Integral variable.
!>    \endverbatim
!>    
!>    \param[in,out] ArrayType
!>    \verbatim
!>        IntegralType is Character(len=*)
!>        = 'space':       the MQC_SCF_Integral variable contains 
!>                         restricted spin blocks.
!>        = 'spin':        the MQC_SCF_Integral variable contains
!>                         unrestricted spin blocks. 
!>        = 'general':     the MQC_SCF_Integral variable contains
!>                         general spin blocks.
!>    \endverbatim
!>
!>    \param[in] Alpha
!>    \verbatim
!>        Alpha is type(mqc_matrix),optional
!>        The alpha spin block integrals.
!>    \endverbatim
!>
!>    \param[in] Beta
!>    \verbatim
!>        Beta is type(mqc_matrix),optional
!>        The beta spin block integrals.
!>    \endverbatim
!>
!>    \param[in] AlphaBeta
!>    \verbatim
!>        AlphaBeta is type(mqc_matrix),optional
!>        The alpha-beta spin block integrals.
!>    \endverbatim
!>
!>    \param[in] BetaAlpha
!>    \verbatim
!>        AlphaBeta is type(mqc_matrix),optional
!>        The beta-alpha spin block integrals.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      subroutine mqc_integral_allocate(integral,arrayName,arrayType,alpha, &
          beta,alphaBeta,betaAlpha) 
!
      implicit none
      class(mqc_scf_integral)::integral
      character(len=*)::arrayName,arrayType
      type(mqc_matrix),optional::alpha,beta,alphaBeta,betaAlpha
!
      call string_change_case(arrayType,'L')
      call string_change_case(arrayName,'L')
      integral%array_type = arrayType 
      integral%array_name = arrayName 
!
      if(mqc_matrix_isallocated(integral%alpha)) call mqc_deallocate_matrix(integral%alpha)
      if(mqc_matrix_isallocated(integral%beta)) call mqc_deallocate_matrix(integral%beta)
      if(mqc_matrix_isallocated(integral%alphaBeta)) call mqc_deallocate_matrix(integral%alphaBeta)
      if(mqc_matrix_isallocated(integral%betaAlpha)) call mqc_deallocate_matrix(integral%betaAlpha)
!
      if(present(alpha)) integral%alpha = alpha
      if(arrayType.eq.'spin') then
        if(present(beta)) integral%beta = beta
      elseIf(arrayType.eq.'general') then
        if(present(beta)) integral%beta = beta
        ! This is completely general - we are not going to assume anything about block size
        if(present(alphaBeta)) integral%alphaBeta = alphaBeta
        if(present(betaAlpha)) integral%betaAlpha = betaAlpha
      endIf
!
      end subroutine mqc_integral_allocate
!
!
!     PROCEDURE MQC_Matrix_To_Integral
!
!>    \brief <b> MQC_Matrix_To_Integral is used to transform MQC_Matrix objects into
!>    MQC_SCF_Integral objects</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Matrix_To_Integral is used to transform MQC_Matrix objects into MQC_SCF_Integral
!>    objects.  The following options are available:
!>
!>    1.  integralType = 'space' declares the MQC_Matrix object to be of a MQC_SCF_Integral
!>                       object with an alpha-alpha component only (restricted spin).
!>    2.  integralType = 'spin' declares the MQC_Matrix object to be of a MQC_SCF_Integral
!>                       object with an alpha-alpha and beta-beta component (unrestricted spin).  
!>                       This option requires nBasis parameter to be set.
!>    3.  integralType = 'general' declares the MQC_Matrix object to be of a MQC_SCF_Integral
!>                       object with all four spin component blocks (general spin).  This option
!>                       requires nBasis parameter to be set.
!>
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Matrix
!>    \verbatim
!>      Matrix is type(mqc_matrix)
!>      The MQC_Matrix to be coverted to MQC_SCF_Integral
!>    \endverbatim
!>
!>    \param[in,out] Integral
!>    \verbatim
!>      Integral is class(mqc_scf_integral)
!>      The MQC_SCF_Integral that is to be allocated.
!>    \endverbatim
!>
!>    \param[in,out] IntegralName
!>    \verbatim
!>      IntegralName is Character(len=*)
!>      The label of the MQC_SCF_Integral object to be allocated.
!>    \endverbatim
!>
!>    \param[in,out] IntegralType
!>    \verbatim
!>      IntegralType is Character(len=*)
!>      = 'space':       the MQC_SCF_Integral variable contains
!>                       restricted spin blocks.
!>      = 'spin':        the MQC_SCF_Integral variable contains
!>                       unrestricted spin blocks.
!>      = 'general':     the MQC_SCF_Integral variable contains
!>                       general spin blocks.
!>    \endverbatim
!>
!>    \param[in,out] nBasis
!>    \verbatim
!>      nBasis is integer(kind=int64), optional
!>      The basis size used to partition Matrix into spin blocks.
!>      Required when IntegralType = 'spin' or 'general'.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author A. Mahler
!>    \date 2020
!
      subroutine mqc_matrix_to_integral(matrix,integral,integralName,integralType,nBasis)
!
        implicit none
        type(mqc_matrix)::matrix
        class(mqc_scf_integral)::integral
        character(len=*)::integralName,integralType
        integer(kind=int64),optional::nBasis
!
        call string_change_case(integralType,'L')
        integral%array_name = integralName
!
        if((integralType == 'space').or.(integralType == 'spin').or.&
          (integralType == 'general')) then         
          integral%array_type = integralType
        else
          call mqc_error_a('Bad integralType specification in MQC_Matrix_To_Integral',6,&
            'integralType = ', integralType)
        endif
!
        if((integralType /= 'space').and..not.(present(nBasis))) then
          call mqc_error('Integral types "spin" and "general" require nBasis argument '&
            //'in MQC_Matrix_To_Integral',6)
        end if
!
        if(integralType == 'space') then
          integral%alpha = matrix
        else if(integralType == 'spin') then
          integral%alpha = matrix%mat([1,nBasis],[1,nBasis])
          integral%beta = matrix%mat([nBasis+1,nBasis*2],[nBasis+1,nBasis*2])
        else
          integral%alpha = matrix%mat([1,nBasis],[1,nBasis])
          integral%beta = matrix%mat([nBasis+1,nBasis*2],[nBasis+1,nBasis*2])
          integral%alphaBeta = matrix%mat([nBasis+1,nBasis*2],[1,nBasis])
          integral%betaAlpha = matrix%mat([1,nBasis],[nBasis+1,nBasis*2])
        endif
!
        end subroutine mqc_matrix_to_integral
!
!
!     PROCEDURE MQC_Eigenvalues_Allocate
!
!>    \brief <b> MQC_Eigenvalues_Allocate is used to allocate a MQC eigenvalues type
!>    variable of the MQC_SCF_Eigenvalues class</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Allocate is used to allocate a MQC SCF eigenvalues type 
!>    variable of the MQC_SCF_Eigenvalues class. The following options are available:
!>
!>    1. ArrayType = 'space' declares the MQC_SCF_Eigenvalues variable contains 
!>                   restricted spin blocks.
!>    2. ArrayType = 'spin' or 'general' declares the MQC_SCF_Eigenvalues variable 
!>                   contains unrestricted spin blocks.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is class(mqc_scf_eigenvalues)
!>        The MQC_SCF_Eigenvalues variable to allocate.
!>    \endverbatim
!>    
!>    \param[in,out] ArrayName
!>    \verbatim
!>        ArrayName is Character(len=*)
!>        The label of the MQC_SCF_Eigenvalues variable.
!>    \endverbatim
!>    
!>    \param[in,out] ArrayType
!>    \verbatim
!>        ArrayType is Character(len=*)
!>        = 'space':       the MQC_SCF_Eigenvalues variable 
!>                         contains restricted spin blocks.
!>        = 'spin':        the MQC_SCF_Eigenvalues variable 
!>                         contains unrestricted spin blocks. 
!>        = 'general':     the MQC_SCF_Eigenvalues variable 
!>                         contains general spin blocks.
!>    \endverbatim
!>
!>    \param[in] Alpha
!>    \verbatim
!>        Alpha is type(mqc_vector),optional
!>        The alpha spin block vector.
!>    \endverbatim
!>
!>    \param[in] Beta
!>    \verbatim
!>        Beta is type(mqc_vector),optional
!>        The beta spin block vector.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      subroutine mqc_eigenvalues_allocate(eigenvalues,arrayName,arrayType,alpha, &
          beta) 
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      character(len=*)::arrayName,arrayType
      type(mqc_vector),optional::alpha,beta
!
      call string_change_case(arrayType,'L')
      call string_change_case(arrayName,'L')
      eigenvalues%array_type = arrayType 
      eigenvalues%array_name = arrayName 
!
      if(present(alpha)) eigenvalues%alpha = alpha
      if(arrayType.eq.'spin'.or.arrayType.eq.'general') then
        if(present(beta)) eigenvalues%beta = beta
      endIf
!
      end subroutine mqc_eigenvalues_allocate  
!
!
!     PROCEDURE MQC_TwoERIs_Transform
!
!>    \brief <b> MQC_TwoERIs_Transform is used to transform two-electron integral 
!>    type variables between different integral types</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_TwoERIs_Transform is used to transform two-electron integral type variables 
!>    between different integral types. nThe following options are available:
!>
!>    1. IntegralType = 'Regular' declares the MQC_TwoERIs variable contains regular 
!>                      atomic orbital 2ERIs. 
!>    2. IntegralType = 'Raffenetti1' declares the MQC_TwoERIs variable contains 
!>                      raffenetti type 1 atomic orbital 2ERIs. 
!>    3. IntegralType = 'Raffenetti2 declares the MQC_TwoERIs variable contains 
!>                      raffenetti type 2 atomic orbital 2ERIs. 
!>    4. IntegralType = 'Raffenetti3' declares the MQC_TwoERIs variable contains 
!>                      raffenetti type 3 atomic orbital 2ERIs. 
!>    5. IntegralType = 'Molecular' declares the MQC_TwoERIs variable 
!>                      contains molecular orbital basis 2ERIs.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] TwoERIs
!>    \verbatim
!>        TwoERIs is class(mqc_twoERIs)
!>        The MQC_TwoERIs variable to transform.
!>    \endverbatim
!>    
!>    \param[in] IntegralType
!>    \verbatim
!>        IntegralType is Character(len=*)
!>        = 'Regular':     the MQC_TwoERIs variable contains regular atomic orbital 
!>                         2ERIs.
!>                         R(i,j,k,l) = (ij|kl)
!>        = 'Raffenetti1': the MQC_TwoERIs variable contains raffenetti type 1 atomic 
!>                         orbital 2ERIs.
!>                         R(i,j,k,l) = (ij|kl) - 1/4[(ik|jl)+(il|jk)]
!>        = 'Raffenetti2': the MQC_TwoERIs variable contains raffenetti type 2 atomic 
!>                         orbital 2ERIs.
!>                         R(i,j,k,l) = (ik|jl) + (il|jk)  
!>        = 'Raffenetti3': the MQC_TwoERIs variable contains raffenetti type 3 atomic 
!>                         orbital 2ERIs.
!>                         R(i,j,k,l) = (ik|jl) - (il|jk)
!>        = 'Molecular':   the MQC_TwoERIs variable contains molecular orbital basis 
!>                         2ERIs.
!>    \endverbatim
!>    
!>    \param[in] mo_coefficients 
!>    \verbatim
!>        mo_coefficients is type(mqc_scf_integral),optional
!>        MO coefficients required for AO -> MO transformations. 
!>    \endverbatim
!>    
!>    \param[in] mo_coeff_2
!>    \verbatim
!>        mo_coeff_2 is type(mqc_scf_integral),optional
!>        Ket MO coefficients required for AO -> MO transformations
!>        if bra and ket orbitals are different. 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2019
!
!
      subroutine mqc_twoeris_transform(twoERIs,integralType,mo_coefficients,mo_coeff_2)
!
      implicit none
      type(mqc_twoERIs)::twoERIs
      character(len=*),intent(in)::integralType
      type(mqc_scf_integral),optional,intent(in)::mo_coefficients,mo_coeff_2
      character(len=64)::my_integralType
      type(mqc_r4tensor)::aaaa
      real(kind=real64),dimension(:),allocatable::aaaaLoc
      type(mqc_twoERIs)::temp
      integer::i,j,k,l,cnt,ip,kp,bs,bsntt,iOut=6,iPrint=4

      call string_change_case(integralType,'l',my_integralType)

      select case (twoERIs%type())
      case ('regular')
        if(my_integralType.eq.'regular') then
          return
        elseIf(my_integralType.eq.'raffenetti1') then
          bs = twoERIs%blockSize('alpha',1)
          bsntt = (bs*(bs+1))/2
          call aaaa%init(bs,bs,bs,bs,0.0)
          call mqc_r4Tensor_full2Symm(aaaa)
          allocate(aaaaLoc((bs**4+2*bs**3+3*bs**2+2*bs)/8))
          aaaaLoc = 0.0
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP REDUCTION(+:aaaaLoc), &
!$OMP SHARED(twoERIs,bs,bsntt),&
!$OMP PRIVATE(cnt,ip,kp,i,j,k,l)
          do cnt = 1, (bsntt*(bsntt+1))/2
#ifdef   PGI
            ip = floor((1+sqrt(dfloat(-7+8*cnt)))/2)
            kp = cnt - ip*(ip-1)/2
            i = floor((1+sqrt(dfloat(-7+8*ip)))/2)
            k = floor((1+sqrt(dfloat(-7+8*kp)))/2)
#else
            ip = floor((1+sqrt(float(-7+8*cnt)))/2)
            kp = cnt - ip*(ip-1)/2
            i = floor((1+sqrt(float(-7+8*ip)))/2)
            k = floor((1+sqrt(float(-7+8*kp)))/2)
#endif  
            j = ip - i*(i-1)/2
            l = kp - k*(k-1)/2

            aaaaLoc(symIndexHash(i,j,k,l)) = twoERIs%at(i,j,k,l) - &
              0.25*(twoERIs%at(i,l,k,j) + twoERIs%at(i,k,l,j))
            if(ip.eq.kp) aaaaLoc(symIndexHash(i,j,k,l)) = &
              0.5*aaaaLoc(symIndexHash(i,j,k,l))

          endDo
!$OMP END PARALLEL DO
          call MQC_Matrix_SymmSymmR4Tensor_Put_Real(aaaa,aaaaLoc)
          call mqc_twoeris_allocate(twoERIs,'symm','raffenetti1',aaaa)
          return
        elseIf(my_integralType.eq.'raffenetti2') then
          bs = twoERIs%blockSize('alpha',1)
          bsntt = (bs*(bs+1))/2
          call aaaa%init(bs,bs,bs,bs,0.0)
          call mqc_r4Tensor_full2Symm(aaaa)
          allocate(aaaaLoc((bs**4+2*bs**3+3*bs**2+2*bs)/8))
          aaaaLoc = 0.0
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP REDUCTION(+:aaaaLoc), &
!$OMP SHARED(twoERIs,bs,bsntt),&
!$OMP PRIVATE(cnt,ip,kp,i,j,k,l)
          do cnt = 1, (bsntt*(bsntt+1))/2
#ifdef   PGI
            ip = floor((1+sqrt(dfloat(-7+8*cnt)))/2)
            kp = cnt - ip*(ip-1)/2
            i = floor((1+sqrt(dfloat(-7+8*ip)))/2)
            k = floor((1+sqrt(dfloat(-7+8*kp)))/2)
#else
            ip = floor((1+sqrt(float(-7+8*cnt)))/2)
            kp = cnt - ip*(ip-1)/2
            i = floor((1+sqrt(float(-7+8*ip)))/2)
            k = floor((1+sqrt(float(-7+8*kp)))/2)
#endif  
            j = ip - i*(i-1)/2
            l = kp - k*(k-1)/2

            aaaaLoc(symIndexHash(i,j,k,l)) = twoERIs%at(i,l,k,j) + &
              twoERIs%at(i,k,l,j)
            if(ip.eq.kp) aaaaLoc(symIndexHash(i,j,k,l)) = &
              0.5*aaaaLoc(symIndexHash(i,j,k,l))

          endDo
!$OMP END PARALLEL DO
          call MQC_Matrix_SymmSymmR4Tensor_Put_Real(aaaa,aaaaLoc)
          call mqc_twoeris_allocate(twoERIs,'symm','raffenetti2',aaaa)
          return
        elseIf(my_integralType.eq.'raffenetti3') then
          bs = twoERIs%blockSize('alpha',1)
          bsntt = (bs*(bs+1))/2
          call aaaa%init(bs,bs,bs,bs,0.0)
          call mqc_r4Tensor_full2Symm(aaaa)
          allocate(aaaaLoc((bs**4+2*bs**3+3*bs**2+2*bs)/8))
          aaaaLoc = 0.0
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP REDUCTION(+:aaaaLoc), &
!$OMP SHARED(twoERIs,bs,bsntt),&
!$OMP PRIVATE(cnt,ip,kp,i,j,k,l)
          do cnt = 1, (bsntt*(bsntt+1))/2
#ifdef   PGI
            ip = floor((1+sqrt(dfloat(-7+8*cnt)))/2)
            kp = cnt - ip*(ip-1)/2
            i = floor((1+sqrt(dfloat(-7+8*ip)))/2)
            k = floor((1+sqrt(dfloat(-7+8*kp)))/2)
#else
            ip = floor((1+sqrt(float(-7+8*cnt)))/2)
            kp = cnt - ip*(ip-1)/2
            i = floor((1+sqrt(float(-7+8*ip)))/2)
            k = floor((1+sqrt(float(-7+8*kp)))/2)
#endif  
            j = ip - i*(i-1)/2
            l = kp - k*(k-1)/2

            aaaaLoc(symIndexHash(i,j,k,l)) = twoERIs%at(i,k,l,j) - &
              twoERIs%at(i,l,k,j)
            if(ip.eq.kp) aaaaLoc(symIndexHash(i,j,k,l)) = &
              0.5*aaaaLoc(symIndexHash(i,j,k,l))

          endDo
!$OMP END PARALLEL DO
          call MQC_Matrix_SymmSymmR4Tensor_Put_Real(aaaa,aaaaLoc)
          call mqc_twoeris_allocate(twoERIs,'symm','raffenetti3',aaaa)
          return
        elseIf(my_integralType.eq.'molecular') then
          if(.not.present(mo_coefficients)) call mqc_error_l('MO oefficients myst be present to&
            & transform AO 2ERIs to MO 2ERIs',6,'present(mo_coefficients)',present(mo_coefficients))
          if(present(mo_coeff_2)) then
            call twoERI_Trans(iOut,iPrint,mo_coefficients,twoERIs,temp,mo_coeff_2)
          else
            call twoERI_Trans(iOut,iPrint,mo_coefficients,twoERIs,temp)
          endIf
          twoERIs = temp
          return
        else
          call mqc_error_a('Unrecognized 2ERIs type requested in mqc_twoeris_transform',6, &
            'my_integralType',my_integralType)
        endIf
      case ('raffenetti1')
        call mqc_error_a('Transformation of Raffenetti 1 2ERIs NYI in mqc_twoeris_transform',6, &
          'twoERIs%type()',twoERIs%type())
      case ('raffenetti2')
        call mqc_error_a('Transformation of Raffenetti 2 2ERIs NYI in mqc_twoeris_transform',6, &
          'twoERIs%type()',twoERIs%type())
      case ('raffenetti3')
        call mqc_error_a('Transformation of Raffenetti 3 2ERIs NYI in mqc_twoeris_transform',6, &
          'twoERIs%type()',twoERIs%type())
      case ('space')
        call mqc_error_a('Transformation of spatial MO 2ERIs NYI in mqc_twoeris_transform',6, &
          'twoERIs%type()',twoERIs%type())
      case ('spin')
        call mqc_error_a('Transformation of spin MO 2ERIs NYI in mqc_twoeris_transform',6, &
          'twoERIs%type()',twoERIs%type())
      case ('general')
        call mqc_error_a('Transformation of general MO 2ERIs NYI in mqc_twoeris_transform',6, &
          'twoERIs%type()',twoERIs%type())
      case default
        call mqc_error_a('Unrecognized 2ERIs integral type in mqc_twoeris_transform',6, &
          'twoERIs%type()',twoERIs%type())
      end select

      end subroutine mqc_twoeris_transform
!
!
!     PROCEDURE MQC_Integral_Identity
!
!>    \brief <b> MQC_Integral_Identity is used to set an MQC integral type variable
!>    of the MQC_SCF_Integral class to the identity matrix</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Identity is used to set an MQC integral type variable of the 
!>    MQC_SCF_Integral class to the identity matrix. Optional argument label allows a
!>    description to be carried around with the variable. Optional arguments nAlpha2
!>    and nBeta2 allow the columns of spin blocks to be sized differently to the rows.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        integral is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable to set to the identiy 
!>        matrix.
!>    \endverbatim
!>    
!>    \param[in] NAlpha
!>    \verbatim
!>        NAlpha is integer
!>        The dimension of the alpha spin block.
!>    \endverbatim
!>    
!>    \param[in] NBeta
!>    \verbatim
!>        NBeta is integer
!>        The dimension of the beta spin block.
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(len=*),Optional
!>        Description of the MQC_SCF_Integral contents.
!>    \endverbatim
!>
!>    \param[in] NAlpha2
!>    \verbatim
!>        NAlpha2 is Integer,Optional
!>        The column number of the alpha spin block, if a 
!>        rectangular spin block is desired.
!>    \endverbatim
!>
!>    \param[in] NBeta2
!>    \verbatim
!>        NBeta2 is Integer,Optional
!>        The column number of the beta spin block, if a 
!>        rectangular spin block is desired.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      subroutine mqc_integral_identity(integral,nAlpha,nBeta,label,nAlpha2,nBeta2) 
!
      implicit none
      class(mqc_scf_integral),intent(inOut)::integral
      integer,intent(in)::nAlpha,nBeta
      integer,intent(in),optional::nAlpha2,nBeta2
      character(Len=*),optional,intent(in)::label
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta
      Character(Len=64)::myLabel
      integer::my_nAlpha2,my_nBeta2
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf
      if(.not.present(nAlpha2)) then
        my_nAlpha2 = nAlpha
      else
        my_nAlpha2 = nAlpha2
      endIf
      if(.not.present(nBeta2)) then
        my_nBeta2 = nBeta
      else
        my_nBeta2 = nBeta2
      endIf
!
      call tmpMatrixAlpha%identity(nAlpha,my_nAlpha2)

      if (nAlpha.eq.nBeta.and.my_nAlpha2.eq.my_nBeta2) then
        call mqc_integral_allocate(integral,myLabel,'space',tmpMatrixAlpha)
      else 
        call tmpMatrixBeta%identity(nBeta,my_nBeta2)
        call mqc_integral_allocate(integral,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
      endIf
!
      end subroutine mqc_integral_identity  
!
!
!     PROCEDURE MQC_Integral_Initialize
!
!>    \brief <b> MQC_Integral_Initialize is used to initialize an MQC integral type 
!>    variable with each element equal to a chosen value</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Initialize is used to initialize an MQC integral type variable 
!>    with each element equal to a chosen value. Optional argument scalar determines
!>    the value to set each element, with the default being zero. Optional argument 
!>    label allows a description to be carried around with the variable. Optional 
!>    arguments nAlpha2 and nBeta2 allow the columns of spin blocks to be sized 
!>    differently to the rows.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        integral is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable to initialize.
!>    \endverbatim
!>    
!>    \param[in] NAlpha
!>    \verbatim
!>        NAlpha is integer
!>        The dimension of the alpha spin block.
!>    \endverbatim
!>    
!>    \param[in] NBeta
!>    \verbatim
!>        NBeta is integer
!>        The dimension of the beta spin block.
!>    \endverbatim
!>
!>    \param[in] Scalar
!>    \verbatim
!>        Scalar is Class(*),Optional
!>        The value to set each element of Integral.
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(len=*),Optional
!>        Description of the MQC_SCF_Integral contents.
!>    \endverbatim
!>
!>    \param[in] NAlpha2
!>    \verbatim
!>        NAlpha2 is Integer,Optional
!>        The column number of the alpha spin block, if a 
!>        rectangular spin block is desired.
!>    \endverbatim
!>
!>    \param[in] NBeta2
!>    \verbatim
!>        NBeta2 is Integer,Optional
!>        The column number of the beta spin block, if a 
!>        rectangular spin block is desired.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      subroutine mqc_integral_initialize(integral,nAlpha,nBeta,scalar,label,nAlpha2,nBeta2) 
!
      implicit none
      class(mqc_scf_integral),intent(inOut)::integral
      integer,intent(in)::nAlpha,nBeta
      character(Len=*),optional,intent(in)::label
      class(*),optional,intent(in)::Scalar
      integer,intent(in),optional::nAlpha2,nBeta2
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta
      Character(Len=64)::myLabel
      integer::my_nAlpha2,my_nBeta2
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf
      if(.not.present(nAlpha2)) then
        my_nAlpha2 = nAlpha
      else
        my_nAlpha2 = nAlpha2
      endIf
      if(.not.present(nBeta2)) then
        my_nBeta2 = nBeta
      else
        my_nBeta2 = nBeta2
      endIf
      if(present(scalar)) then
        call tmpMatrixAlpha%init(nAlpha,my_nAlpha2,scalar)
        if (nAlpha.eq.nBeta.and.my_nAlpha2.eq.my_nBeta2) then
          call mqc_integral_allocate(integral,myLabel,'space',tmpMatrixAlpha)
        else 
          call tmpMatrixBeta%init(nBeta,my_nBeta2,scalar)
          call mqc_integral_allocate(integral,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
        endIf
      else
        call tmpMatrixAlpha%init(nAlpha,my_nAlpha2)
        if (nAlpha.eq.nBeta.and.my_nAlpha2.eq.my_nBeta2) then
          call mqc_integral_allocate(integral,myLabel,'space',tmpMatrixAlpha)
        else 
          call tmpMatrixBeta%init(nBeta,my_nBeta2)
          call mqc_integral_allocate(integral,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
        endIf
      endIf
!
!
      end subroutine mqc_integral_initialize
!
!
!     PROCEDURE MQC_Integral_Output_Block
!
!>    \brief <b> MQC_Integral_Output_Block is used to output a spin block of an MQC 
!>    integral type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Output_Block is used to output a spin block of an MQC integral 
!>    type variable. Optional argument blockName is used to specify which spin block
!>    should be output, with the default to output the all spin blocks. The following
!>    options are available:
!>
!>    1. BlockName = 'full' outputs the full MQC_SCF_Integral. 
!>    2. BlockName = 'alpha' or 'alpha-alpha' outputs the alpha-alpha spin block. 
!>    3. BlockName = 'beta' or 'beta-beta' outputs the beta-beta spin block. 
!>    4. BlockName = 'alpha-beta' outputs the alpha-beta spin block. 
!>    5. BlockName = 'beta-alpha' outputs the beta-alpha spin block. 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        integral is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable to output.
!>    \endverbatim
!>    
!>    \param[in,out] BlockName
!>    \verbatim
!>        BlockName is character(len=*),optional
!>        = 'full':        outputs the full MQC_SCF_Integral. 
!>        = 'alpha':       outputs the alpha-alpha spin block.
!>        = 'alpha-alpha': outputs the alpha-alpha spin block. 
!>        = 'beta':        outputs the beta-beta spin block. 
!>        = 'beta-beta':   outputs the beta-beta spin block. 
!>        = 'alpha-beta'   outputs the alpha-beta spin block. 
!>        = 'beta-alpha'   outputs the beta-alpha spin block. 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      function mqc_integral_output_block(integral,blockName) result(matrixOut)
!
      implicit none
      class(mqc_scf_integral)::integral
      character(len=*),optional::blockName
      character(len=64)::myBlockName
      type(mqc_matrix)::matrixOut
      integer(kind=int64)::nDimAlpha1=0,nDimBeta1=0,nDimTotal1=0,nDimAlpha2=0,nDimBeta2=0,nDimTotal2=0,zeroI
      real(kind=real64)::zeroR
      complex(kind=real64)::zeroC
      parameter(zeroI=0,zeroR=0.0,zeroC=(0.0,0.0))
!
      if(present(blockName)) then
        call string_change_case(blockName,'l',myBlockName)
      else
        myBlockName = 'full'
      endIf

      nDimAlpha1 = integral%blockSize('Alpha',1) 
      nDimBeta1 = integral%blockSize('Beta',1)
      nDimTotal1 = nDimAlpha1 + nDimBeta1
      nDimAlpha2 = integral%blockSize('Alpha',2) 
      nDimBeta2 = integral%blockSize('Beta',2)
      nDimTotal2 = nDimAlpha2 + nDimBeta2

      select case (myBlockName)
      case('full')
        if (integral%type().eq.'space') then
          if (integral%hasAlpha()) then
            if(mqc_matrix_haveInteger(integral%alpha)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroI)
            elseif(mqc_matrix_haveReal(integral%alpha)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroR)
            elseIf(mqc_matrix_haveComplex(integral%alpha)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroC)
            else
              call mqc_error('unrecognized integral data type in mqc_integral_output_block')
            endIf
            call matrixOut%mput(integral%alpha,[1,nDimAlpha1],[1,nDimAlpha2]) 
            call matrixOut%mput(integral%alpha,[nDimAlpha1+1,nDimTotal1],[nDimAlpha2+1,nDimTotal2]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlpha()', integral%hasAlpha() )
          endIf
        elseIf (integral%type().eq.'spin') then
          if (integral%hasAlpha().and.integral%hasBeta()) then
            if(mqc_matrix_haveComplex(integral%alpha).or.mqc_matrix_haveComplex(integral%beta)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroC)
            elseIf(mqc_matrix_haveReal(integral%alpha).or.mqc_matrix_haveReal(integral%beta)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroR)
            elseIf(mqc_matrix_haveInteger(integral%alpha).or.mqc_matrix_haveInteger(integral%beta)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroI)
            else
              call mqc_error('unrecognized integral data type in mqc_integral_output_block')
            endIf
            call matrixOut%mput(integral%alpha,[1,nDimAlpha1],[1,nDimAlpha2]) 
            call matrixOut%mput(integral%beta,[nDimAlpha1+1,nDimTotal1],[nDimAlpha2+1,nDimTotal2]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlpha()', integral%hasAlpha(), 'integral%hasBeta()', &
                 integral%hasBeta())
          endIf
        elseIf (integral%type().eq.'general') then
          if (integral%hasAlpha().and.integral%hasBeta().and.integral%hasAlphaBeta().and.integral%hasBetaAlpha()) then
            if(mqc_matrix_haveComplex(integral%alpha).or.mqc_matrix_haveComplex(integral%beta).or. &
              mqc_matrix_haveComplex(integral%alphaBeta).or.mqc_matrix_haveComplex(integral%betaAlpha)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroC)
            elseIf(mqc_matrix_haveReal(integral%alpha).or.mqc_matrix_haveReal(integral%beta).or. &
              mqc_matrix_haveReal(integral%alphaBeta).or.mqc_matrix_haveReal(integral%betaAlpha)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroR)
            elseIf(mqc_matrix_haveInteger(integral%alpha).or.mqc_matrix_haveInteger(integral%beta).or. &
              mqc_matrix_haveInteger(integral%alphaBeta).or.mqc_matrix_haveInteger(integral%betaAlpha)) then
              call matrixOut%init(nDimTotal1,nDimTotal2,zeroI)
            else
              call mqc_error('unrecognized integral data type in mqc_integral_output_block')
            endIf
            call matrixOut%mput(integral%alpha,[1,nDimAlpha1],[1,nDimAlpha2]) 
            call matrixOut%mput(integral%beta,[nDimAlpha1+1,nDimTotal1],[nDimAlpha2+1,nDimTotal2]) 
            call matrixOut%mput(integral%alphaBeta,[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2]) 
            call matrixOut%mput(integral%betaAlpha,[1,nDimAlpha1],[nDimAlpha2+1,nDimTotal2]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlpha()', integral%hasAlpha(), 'integral%hasBeta()', &
                 integral%hasBeta(), 'integral%hasAlphaBeta()', integral%hasAlphaBeta(), &
                 'integral%hasBetaAlpha()', integral%hasBetaAlpha() )
          endIf
        endIf
      case('alpha')
        if (integral%hasAlpha()) then
          matrixOut = integral%alpha
        else
          call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
               'integral%hasAlpha()', integral%hasAlpha() )
        endIf
      case('beta')
        if (integral%hasBeta()) then
          matrixOut = integral%beta
        elseIf (integral%type().eq.'space') then
          if (integral%hasAlpha()) then
            matrixOut = integral%alpha
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlpha()', integral%hasAlpha() )
          endIf
        else
          call mqc_error_A('block does not exist in mqc_integral_output_block', 6, &
               'integral%type()', integral%type() )
        endIf
      case('alpha-alpha')
        if (integral%hasAlpha()) then
          matrixOut = integral%alpha
        else
          call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
               'integral%hasAlpha()', integral%hasAlpha() )
        endIf
      case('alpha-beta')
        if (integral%hasAlphaBeta()) then
          matrixOut = integral%alphaBeta
        else
          if (integral%type().eq.'space') then
            if(mqc_matrix_haveInteger(integral%alpha)) then
              call matrixOut%init(NDimAlpha1,NDimAlpha2,zeroI)
            elseIf(mqc_matrix_haveReal(integral%alpha)) then
              call matrixOut%init(NDimAlpha1,NDimAlpha2,zeroR)
            elseIf(mqc_matrix_haveComplex(integral%alpha)) then
              call matrixOut%init(NDimAlpha1,NDimAlpha2,zeroC)
            else
              call matrixOut%init(NDimAlpha1,NDimAlpha2,zeroR)
            endIf
          elseIf (integral%type().eq.'spin') then
            call matrixOut%init(NDimAlpha1,NDimBeta2)
            if(mqc_matrix_haveComplex(integral%alpha).or.mqc_matrix_haveComplex(integral%beta)) then
              call matrixOut%init(NDimAlpha1,NDimBeta2,zeroC)
            elseIf(mqc_matrix_haveReal(integral%alpha).or.mqc_matrix_haveReal(integral%beta)) then
              call matrixOut%init(NDimAlpha1,NDimBeta2,zeroR)
            elseIf(mqc_matrix_haveInteger(integral%alpha).or.mqc_matrix_haveReal(integral%beta)) then
              call matrixOut%init(NDimAlpha1,NDimBeta2,zeroI)
            else
              call matrixOut%init(NDimAlpha1,NDimBeta2,zeroR)
            endIf
          endIf
        endIf
      case('beta-alpha')
        if (integral%hasBetaAlpha()) then
          matrixOut = integral%betaAlpha
        else
          if (integral%type().eq.'space') then
            if(mqc_matrix_haveInteger(integral%alpha)) then
              call matrixOut%init(NDimAlpha1,NDimAlpha2,zeroI)
            elseIf(mqc_matrix_haveReal(integral%alpha)) then
              call matrixOut%init(NDimAlpha1,NDimAlpha2,zeroR)
            elseIf(mqc_matrix_haveComplex(integral%alpha)) then
              call matrixOut%init(NDimAlpha1,NDimAlpha2,zeroC)
            else
              call matrixOut%init(NDimAlpha1,NDimAlpha2,zeroR)
            endIf
          elseIf (integral%type().eq.'spin') then
            if(mqc_matrix_haveComplex(integral%alpha).or.mqc_matrix_haveComplex(integral%beta)) then
              call matrixOut%init(NDimBeta1,NDimAlpha2,zeroC)
            elseIf(mqc_matrix_haveReal(integral%alpha).or.mqc_matrix_haveReal(integral%beta)) then
              call matrixOut%init(NDimBeta1,NDimAlpha2,zeroR)
            elseIf(mqc_matrix_haveInteger(integral%alpha).or.mqc_matrix_haveReal(integral%beta)) then
              call matrixOut%init(NDimBeta1,NDimAlpha2,zeroI)
            else
              call matrixOut%init(NDimBeta1,NDimAlpha2,zeroR)
            endIf
          endIf
        endIf
      case('beta-beta')
        if (integral%hasBeta()) then
          matrixOut = integral%beta
        elseIf (integral%type().eq.'space') then
          if (integral%hasAlpha()) then
            matrixOut = integral%alpha
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlpha()', integral%hasAlpha() )
          endIf
        else
          call mqc_error_A('block does not exist in mqc_integral_output_block', 6, &
               'integral%type()', integral%type() )
        endIf
      case default
        call mqc_error_A('unrecognised block name in mqc_integral_output_block', 6, &
             'myBlockName', myBlockName )
      end select
!
      end function mqc_integral_output_block  
!
!
!     PROCEDURE MQC_Integral_Output_Orbitals
!
!>    \brief <b> MQC_Integral_Output_Orbitals is used to output a vector subset from
!>    an MQC integral type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Output_Orbitals is used to output a vector subset from an MQC 
!>    integral type variable. Optional argument orbString enables automatic selection
!>    of common orbital subsets. Optional arguments alphaOrbsIn and betaOrbsIn can be
!>    used to manual specify which vectors should be selected, as well as their order.
!>    Optional argument axis specifies whether row or column vectors are selected. 
!>    The following options are available:
!>
!>    1. OrbString = 'useStrings' specifies that the vector subset will be specified
!>                   in alphaOrbsIn and betaOrbsIn is specified. This
!>                   is the default.
!>    2. OrbString = 'occupied' selects the set of occupied orbitals. Number of alpha 
!>                   electrons should be specified in a length 1 vector in 
!>                   alphaOrbsIn. Number of beta electrons should be specified in a 
!>                   length 1 vector in betaOrbsIn.
!>    3. OrbString = 'virtual' selects the set of virtual orbitals. Number of alpha 
!>                   electrons should be specified in a length 1 vector in 
!>                   alphaOrbsIn. Number of beta electrons should be specified in a 
!>                   length 1 vector in betaOrbsIn. The dimension of the fixed length 
!>                   axis is assumed to be the number of virtual orbitals.
!>    4. OrbString = 'real' selects the real part of input orbitals. A subset will be 
!>                   returned if specified in alphaOrbsIn and betaOrbsIn.
!>    5. OrbString = 'imaginary' selects the imaginary part of input orbitals. A subset 
!>                   will be returned if specified in alphaOrbsIn and betaOrbsIn.
!>    6. OrbString = 'complex' returns the input orbitals as complex, regardless of 
!>                   whether the imaginary part is zero. A subset will be returned if 
!>                   specified in alphaOrbsIn and betaOrbsIn.
!>    7. Axis =      1 specifies the row vector subset is selected. 
!>    8. Axis =      2 specifies the column vector subset is selected (default). 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        integral is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable to select orbitals.
!>    \endverbatim
!>    
!>    \param[in] OrbString 
!>    \verbatim
!>        OrbString is character(len=*),optional
!>        = 'useStrings': selects orbitals from strings.
!>        = 'occupied':   selects occupied orbitals.
!>        = 'virtual':    selects virtual orbitals.
!>        = 'real':       selects real part of orbitals.
!>        = 'imaginary':  selects imaginary part of orbitals.
!>        = 'complex':    returns complex orbitals.
!>    \endverbatim
!>    
!>    \param[in] AlphaOrbsIn
!>    \verbatim
!>        AlphaOrbsIn is integer(kind=int64),dimension(:),optional
!>        Vector of alpha vector subset.
!>    \endverbatim
!>    
!>    \param[in] BetaOrbsIn
!>    \verbatim
!>        BetaOrbsIn is integer(kind=int64),dimension(:),optional
!>        Vector of beta vector subset.
!>    \endverbatim
!>    
!>    \param[in] Axis
!>    \verbatim
!>        Axis is integer(kind=int64),optional
!>        = 1: specifies the row vector subset is selected. 
!>        = 2: specifies the column vector subset is selected (default). 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_integral_output_orbitals(integral,orbString,alphaOrbsIn,betaOrbsIn,axis) result(integralOut)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      character(len=*),optional,intent(in)::orbString
      character(len=64)::myOrbString
      integer(kind=int64),dimension(:),optional,intent(in)::alphaOrbsIn,betaOrbsIn
      integer(kind=int64),dimension(:),allocatable::alphaOrbs,betaOrbs
      integer(kind=int64),intent(in),optional::axis
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      type(mqc_matrix)::outMatrixAlpha,outMatrixBeta,outMatrixAlphaBeta,outMatrixBetaAlpha
      integer(kind=int64)::nDimAlpha1=0,nDimBeta1=0,nDimAlpha2=0,nDimBeta2=0,nAlpha=0,nBeta=0,i,my_axis

      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 2
      endIf

      if(present(orbString)) then
        call string_change_case(orbString,'l',myOrbString)
        if(myOrbString.eq.'occupied'.or.myOrbString.eq.'virtual') then
          if(.not.present(alphaOrbsIn).or..not.present(betaOrbsIn)) &
            call mqc_error_l('OrbString variable must be declared with &
            &number of alpha and/or beta electrons in string input',6,'present(orbString)', &
            present(orbString),'present(alphaOrbsIn)',present(alphaOrbsIn),'present(betaOrbsIn)', &
            present(betaOrbsIn))
          if(size(alphaOrbsIn).ne.1.or.size(betaOrbsIn).ne.1) &
            call mqc_error_i('Only one integer can be declared when inputting number of alpha and beta orbitals',&
            6,'size(alphaOrbsIn)',size(alphaOrbsIn),'size(betaOrbsIn)',size(betaOrbsIn))
          nAlpha = alphaOrbsIn(1)
          nBeta = betaOrbsIn(1) 
        endIf
      else
        myOrbString = 'usestrings'
      endIf

      nDimAlpha1 = integral%blockSize('Alpha',1) 
      nDimBeta1 = integral%blockSize('Beta',1)
      nDimAlpha2 = integral%blockSize('Alpha',2) 
      nDimBeta2 = integral%blockSize('Beta',2)

      select case (myOrbString)
      case ('occupied')
        allocate(alphaOrbs(nAlpha),betaOrbs(nBeta))
        alphaOrbs = [(i,i=1,nAlpha)]
        betaOrbs = [(i,i=1,nBeta)]
      case ('virtual')
        if(my_axis.eq.2) then
          allocate(alphaOrbs(nDimAlpha1-nAlpha),betaOrbs(nDimBeta1-nBeta))
          alphaOrbs = [(i,i=nAlpha+1,nDimAlpha1)]
          betaOrbs = [(i,i=nBeta+1,nDimBeta1)]
        elseIf(my_axis.eq.1) then
          allocate(alphaOrbs(nDimAlpha2-nAlpha),betaOrbs(nDimBeta2-nBeta))
          alphaOrbs = [(i,i=nAlpha,nDimAlpha2)]
          betaOrbs = [(i,i=nBeta,nDimBeta2)]
        else
          call mqc_error_I('Unknown array is in MQC_Integral_Output_Orbitals',6,'my_axis',my_axis)
        endIf
      case ('usestrings')
        if(present(alphaOrbsIn)) then
          alphaOrbs = alphaOrbsIn
        else
          alphaOrbs = [integer(kind=int64)::]
        endIf
        if(present(betaOrbsIn)) then
          betaOrbs = betaOrbsIn
        else
          betaOrbs = [integer(kind=int64)::]
        endIf
      case ('real','imaginary','complex')
        if(present(alphaOrbsIn)) then
          alphaOrbs = alphaOrbsIn
        elseif(my_axis.eq.2) then
          alphaOrbs = [(i,i=1,nDimAlpha1)]
        elseIf(my_axis.eq.1) then
          alphaOrbs = [(i,i=1,nDimAlpha2)]
        endIf
        if(present(betaOrbsIn)) then
          betaOrbs = betaOrbsIn
        elseif(my_axis.eq.2) then
          betaOrbs = [(i,i=1,nDimBeta1)]
        elseIf(my_axis.eq.1) then
          betaOrbs = [(i,i=1,nDimBeta2)]
        endIf
      case default
        call mqc_error_a('Unrecognized orbital string definition in MQC_Integral_Output_Orbitals',&
          6,'myOrbString',myOrbString)
      end select

      select case (integral%type())
      case('space')
        if(size(alphaOrbs).eq.size(betaOrbs).and.all(alphaOrbs.eq.betaOrbs)) then
          tmpMatrixAlpha = integral%getBlock('alpha')
          if(myOrbString.eq.'real') tmpMatrixAlpha = real(tmpMatrixAlpha)
          if(myOrbString.eq.'imaginary') tmpMatrixAlpha = aimag(tmpMatrixAlpha)
          if(myOrbString.eq.'complex') tmpMatrixAlpha = MQC_Matrix_Cmplx(real(tmpMatrixAlpha),aimag(tmpMatrixAlpha))
          if(my_axis.eq.2) then
            call outMatrixAlpha%init(nDimAlpha1,size(alphaOrbs))
            do i = 1,size(alphaOrbs)
              call outMatrixAlpha%vput(tmpMatrixAlpha%vat([0],[abs(alphaOrbs(i))]),[0],[i])
            endDo
          elseIf(my_axis.eq.1) then
            call outMatrixAlpha%init(size(alphaOrbs),nDimAlpha2)
            do i = 1,size(alphaOrbs)
              call outMatrixAlpha%vput(tmpMatrixAlpha%vat([abs(alphaOrbs(i))],[0]),[i],[0])
            endDo
          else
            call mqc_error_I('Unknown array is in MQC_Integral_Output_Orbitals',6,'my_axis',my_axis)
          endIf
          call mqc_integral_allocate(integralOut,'mo coefficients','space',outMatrixAlpha)
        else
          tmpMatrixAlpha = integral%getBlock('alpha')
          tmpMatrixBeta = integral%getBlock('beta')
          if(myOrbString.eq.'real') then
            tmpMatrixAlpha = real(tmpMatrixAlpha)
            tmpMatrixBeta = real(tmpMatrixBeta)
          endIf
          if(myOrbString.eq.'imaginary') then
            tmpMatrixAlpha = aimag(tmpMatrixAlpha)
            tmpMatrixBeta = aimag(tmpMatrixBeta)
          endIf
          if(myOrbString.eq.'complex') then
            tmpMatrixAlpha = MQC_Matrix_Cmplx(real(tmpMatrixAlpha),aimag(tmpMatrixAlpha))
            tmpMatrixBeta = MQC_Matrix_Cmplx(real(tmpMatrixBeta),aimag(tmpMatrixBeta))
          endIf
          if(my_axis.eq.2) then
            call outMatrixAlpha%init(nDimAlpha1,size(alphaOrbs))
            call outMatrixBeta%init(nDimBeta1,size(betaOrbs))
            do i = 1,size(alphaOrbs)
              call outMatrixAlpha%vput(tmpMatrixAlpha%vat([0],[abs(alphaOrbs(i))]),[0],[i])
            endDo
            do i = 1,size(betaOrbs)
              call outMatrixBeta%vput(tmpMatrixBeta%vat([0],[abs(betaOrbs(i))]),[0],[i])
            endDo
          elseIf(my_axis.eq.1) then
            call outMatrixAlpha%init(size(alphaOrbs),nDimAlpha2)
            call outMatrixBeta%init(size(betaOrbs),nDimBeta2)
            do i = 1,size(alphaOrbs)
              call outMatrixAlpha%vput(tmpMatrixAlpha%vat([abs(alphaOrbs(i))],[0]),[i],[0])
            endDo
            do i = 1,size(betaOrbs)
              call outMatrixBeta%vput(tmpMatrixBeta%vat([abs(betaOrbs(i))],[0]),[i],[0])
            endDo
          else
            call mqc_error_I('Unknown array is in MQC_Integral_Output_Orbitals',6,'my_axis',my_axis)
          endIf
          call mqc_integral_allocate(integralOut,'mo coefficients','spin',outMatrixAlpha, &
            outMatrixBeta)
        endIf
      case('spin')
        tmpMatrixAlpha = integral%getBlock('alpha')
        tmpMatrixBeta = integral%getBlock('beta')
        if(myOrbString.eq.'real') then
          tmpMatrixAlpha = real(tmpMatrixAlpha)
          tmpMatrixBeta = real(tmpMatrixBeta)
        endIf
        if(myOrbString.eq.'imaginary') then
          tmpMatrixAlpha = aimag(tmpMatrixAlpha)
          tmpMatrixBeta = aimag(tmpMatrixBeta)
        endIf
        if(myOrbString.eq.'complex') then
          tmpMatrixAlpha = MQC_Matrix_Cmplx(real(tmpMatrixAlpha),aimag(tmpMatrixAlpha))
          tmpMatrixBeta = MQC_Matrix_Cmplx(real(tmpMatrixBeta),aimag(tmpMatrixBeta))
        endIf
        if(my_axis.eq.2) then
          call outMatrixAlpha%init(nDimAlpha1,size(alphaOrbs))
          call outMatrixBeta%init(nDimBeta1,size(betaOrbs))
          do i = 1,size(alphaOrbs)
            if(alphaOrbs(i).ge.0) then
              call outMatrixAlpha%vput(tmpMatrixAlpha%vat([0],[alphaOrbs(i)]),[0],[i])
            else
              call outMatrixAlpha%vput(tmpMatrixBeta%vat([0],[abs(alphaOrbs(i))]),[0],[i])
            endIf
          endDo
          do i = 1,size(betaOrbs)
            if(betaOrbs(i).ge.0) then
              call outMatrixBeta%vput(tmpMatrixBeta%vat([0],[betaOrbs(i)]),[0],[i])
            else
              call outMatrixBeta%vput(tmpMatrixAlpha%vat([0],[abs(alphaOrbs(i))]),[0],[i])
            endIf
          endDo
        elseIf(my_axis.eq.1) then
          call outMatrixAlpha%init(size(alphaOrbs),nDimAlpha2)
          call outMatrixBeta%init(size(betaOrbs),nDimBeta2)
          do i = 1,size(alphaOrbs)
            if(alphaOrbs(i).ge.0) then
              call outMatrixAlpha%vput(tmpMatrixAlpha%vat([alphaOrbs(i)],[0]),[i],[0])
            else
              call outMatrixAlpha%vput(tmpMatrixBeta%vat([abs(alphaOrbs(i))],[0]),[i],[0])
            endIf
          endDo
          do i = 1,size(betaOrbs)
            if(betaOrbs(i).ge.0) then
              call outMatrixBeta%vput(tmpMatrixBeta%vat([betaOrbs(i)],[0]),[i],[0])
            else
              call outMatrixBeta%vput(tmpMatrixAlpha%vat([abs(alphaOrbs(i))],[0]),[i],[0])
            endIf
          endDo
        else
          call mqc_error_I('Unknown array is in MQC_Integral_Output_Orbitals',6,'my_axis',my_axis)
        endIf
        call mqc_integral_allocate(integralOut,'mo coefficients','spin',outMatrixAlpha, &
          outMatrixBeta)
      case('general')
        tmpMatrixAlpha = integral%getBlock('alpha')
        tmpMatrixBeta = integral%getBlock('beta')
        tmpMatrixAlphaBeta = integral%getBlock('alpha-beta')
        tmpMatrixBetaAlpha = integral%getBlock('beta-alpha')
        if(myOrbString.eq.'real') then
          tmpMatrixAlpha = real(tmpMatrixAlpha)
          tmpMatrixBeta = real(tmpMatrixBeta)
          tmpMatrixAlphaBeta = real(tmpMatrixAlphaBeta)
          tmpMatrixBetaAlpha = real(tmpMatrixBetaAlpha)
        endIf
        if(myOrbString.eq.'imaginary') then
          tmpMatrixAlpha = aimag(tmpMatrixAlpha)
          tmpMatrixBeta = aimag(tmpMatrixBeta)
          tmpMatrixAlphaBeta = aimag(tmpMatrixAlphaBeta)
          tmpMatrixBetaAlpha = aimag(tmpMatrixBetaAlpha)
        endIf
        if(myOrbString.eq.'complex') then
          tmpMatrixAlpha = MQC_Matrix_Cmplx(real(tmpMatrixAlpha),aimag(tmpMatrixAlpha))
          tmpMatrixBeta = MQC_Matrix_Cmplx(real(tmpMatrixBeta),aimag(tmpMatrixBeta))
          tmpMatrixAlphaBeta = MQC_Matrix_Cmplx(real(tmpMatrixAlphaBeta),aimag(tmpMatrixAlphaBeta))
          tmpMatrixBetaAlpha = MQC_Matrix_Cmplx(real(tmpMatrixBetaAlpha),aimag(tmpMatrixBetaAlpha))
        endIf
        if(my_axis.eq.2) then
          call outMatrixAlpha%init(nDimAlpha1,size(alphaOrbs))
          call outMatrixBeta%init(nDimBeta1,size(betaOrbs))
          call outMatrixBetaAlpha%init(nDimAlpha1,size(betaOrbs))
          call outMatrixAlphaBeta%init(nDimBeta1,size(alphaOrbs))
          do i = 1,size(alphaOrbs)
            if(alphaOrbs(i).ge.0) then
              call outMatrixAlpha%vput(tmpMatrixAlpha%vat([0],[alphaOrbs(i)]),[0],[i])
              call outMatrixAlphaBeta%vput(tmpMatrixAlphaBeta%vat([0],[alphaOrbs(i)]),[0],[i])
            else
              call outMatrixAlphaBeta%vput(tmpMatrixBeta%vat([0],[abs(alphaOrbs(i))]),[0],[i])
              call outMatrixAlpha%vput(tmpMatrixBetaAlpha%vat([0],[abs(alphaOrbs(i))]),[0],[i])
            endIf
          endDo
          do i = 1,size(betaOrbs)
            if(betaOrbs(i).ge.0) then
              call outMatrixBeta%vput(tmpMatrixBeta%vat([0],[betaOrbs(i)]),[0],[i])
              call outMatrixBetaAlpha%vput(tmpMatrixBetaAlpha%vat([0],[betaOrbs(i)]),[0],[i])
            else
              call outMatrixBetaAlpha%vput(tmpMatrixAlpha%vat([0],[abs(betaOrbs(i))]),[0],[i])
              call outMatrixBeta%vput(tmpMatrixAlphaBeta%vat([0],[abs(betaOrbs(i))]),[0],[i])
            endIf
          endDo
        elseIf(my_axis.eq.1) then
          call outMatrixAlpha%init(size(alphaOrbs),nDimAlpha2)
          call outMatrixBeta%init(size(betaOrbs),nDimBeta2)
          call outMatrixBetaAlpha%init(size(betaOrbs),nDimAlpha2)
          call outMatrixAlphaBeta%init(size(alphaOrbs),nDimBeta2)
          do i = 1,size(alphaOrbs)
            if(alphaOrbs(i).ge.0) then
              call outMatrixAlpha%vput(tmpMatrixAlpha%vat([alphaOrbs(i)],[0]),[i],[0])
              call outMatrixAlphaBeta%vput(tmpMatrixAlphaBeta%vat([alphaOrbs(i)],[0]),[i],[0])
            else
              call outMatrixAlphaBeta%vput(tmpMatrixBeta%vat([abs(alphaOrbs(i))],[0]),[i],[0])
              call outMatrixAlpha%vput(tmpMatrixBetaAlpha%vat([abs(alphaOrbs(i))],[0]),[i],[0])
            endIf
          endDo
          do i = 1,size(betaOrbs)
            if(betaOrbs(i).ge.0) then
              call outMatrixBeta%vput(tmpMatrixBeta%vat([betaOrbs(i)],[0]),[0],[i])
              call outMatrixBetaAlpha%vput(tmpMatrixBetaAlpha%vat([betaOrbs(i)],[0]),[0],[i])
            else
              call outMatrixBetaAlpha%vput(tmpMatrixAlpha%vat([abs(betaOrbs(i))],[0]),[i],[0])
              call outMatrixBeta%vput(tmpMatrixAlphaBeta%vat([abs(betaOrbs(i))],[0]),[i],[0])
            endIf
          endDo
        else
          call mqc_error_I('Unknown array is in MQC_Integral_Output_Orbitals',6,'my_axis',my_axis)
        endIf
        call mqc_integral_allocate(integralOut,'mo coefficients','general',outMatrixAlpha, &
          outMatrixBeta,outMatrixAlphaBeta,outMatrixBetaAlpha)
      case default
        call mqc_error_A('Unrecognised mqc scf integral type in mqc_integral_output_orbitals',6, &
          'integral%type()',integral%type())
      end select
!
      end function mqc_integral_output_orbitals
!
!
!     PROCEDURE MQC_Integral_Swap_Orbitals
!
!>    \brief <b> MQC_Integral_Swap_Orbitals is used to exchange a vector pair in an
!>    MQC integral type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Swap_Orbitals is used to exchange a vector pair in an MQC integral 
!>    type variable. Optional argument axis specifies whether row 
!>    or column vectors are selected. The following options are available:
!>
!>    1. Axis =      1 specifies the row vector subset is selected. 
!>    2. Axis =      2 specifies the column vector subset is selected (default). 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        integral is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable to select orbitals.
!>    \endverbatim
!>    
!>    \param[in] AlphaOrbsIn
!>    \verbatim
!>        AlphaOrbsIn is integer(kind=int64),dimension(2),optional
!>        Alpha/restricted vectors to swap.
!>    \endverbatim
!>    
!>    \param[in] BetaOrbsIn
!>    \verbatim
!>        BetaOrbsIn is integer(kind=int64),dimension(2),optional
!>        Beta vectors to swap.
!>    \endverbatim
!>    
!>    \param[in] Axis
!>    \verbatim
!>        Axis is integer(kind=int64),optional
!>        = 1: specifies the row vector subset is selected. 
!>        = 2: specifies the column vector subset is selected (default). 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_integral_swap_orbitals(integral,alphaOrbsIn,betaOrbsIn,axis) result(integralOut)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      integer(kind=int64),dimension(2),optional::alphaOrbsIn,betaOrbsIn
      integer(kind=int64),dimension(:),allocatable::alphaOrbs,betaOrbs
      integer(kind=int64),intent(in),optional::axis
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      type(mqc_matrix)::outMatrixAlpha,outMatrixBeta,outMatrixAlphaBeta,outMatrixBetaAlpha
      integer(kind=int64),dimension(:),allocatable::alphaVector,betaVector
      integer(kind=int64)::nDimAlpha=0,nDimBeta=0,nAlpha=0,nBeta=0,i,my_axis

      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 2
      endIf

      if(my_axis.eq.1) then
        nDimAlpha = integral%blockSize('Alpha',1) 
        nDimBeta = integral%blockSize('Beta',1)
      elseIf(my_axis.eq.2) then
        nDimAlpha = integral%blockSize('Alpha',2) 
        nDimBeta = integral%blockSize('Beta',2)
      else
        call mqc_error_I('Unknown array is in MQC_Integral_Swap_Orbitals',6,'my_axis',my_axis)
      endIf

      if(.not.present(alphaOrbsIn).and..not.present(betaOrbsIn)) then
        call mqc_error_l('Orbital variable must be declared in mqc_integral_swap_orbitals',&
          6,'present(alphaOrbsIn)',present(alphaOrbsIn),'present(betaOrbsIn)',present(betaOrbsIn))
      elseIf(.not.present(alphaOrbsIn)) then
        betaOrbs = betaOrbsIn
        if(betaOrbs(1).gt.betaOrbs(2)) betaOrbs = [betaOrbs(2),betaOrbs(1)]
        alphaVector = [(i,i=1,NDimAlpha)]
        if(betaOrbs(1).eq.betaOrbs(2)) then
          betaVector = [(i,i=1,NDimBeta)]
        else
          betaVector = [(i,i=1,betaOrbs(1)-1),betaOrbs(2),(i,i=betaOrbs(1)+1,betaOrbs(2)-1),&
            betaOrbs(1),(i,i=betaOrbs(2)+1,NDimBeta)]
        endIf
      elseIf(.not.present(betaOrbsIn)) then
        alphaOrbs = alphaOrbsIn
        if(alphaOrbs(1).gt.alphaOrbs(2)) alphaOrbs = [alphaOrbs(2),alphaOrbs(1)]
        if(alphaOrbs(1).eq.alphaOrbs(2)) then
          alphaVector = [(i,i=1,NDimAlpha)]
        else
          alphaVector = [(i,i=1,alphaOrbs(1)-1),alphaOrbs(2),(i,i=alphaOrbs(1)+1,alphaOrbs(2)-1),&
            alphaOrbs(1),(i,i=alphaOrbs(2)+1,NDimAlpha)]
        endIf
        betaVector = [(i,i=1,NDimBeta)]
      else
        alphaOrbs = alphaOrbsIn
        betaOrbs = betaOrbsIn
        if(alphaOrbs(1).gt.alphaOrbs(2)) alphaOrbs = [alphaOrbs(2),alphaOrbs(1)]
        if(betaOrbs(1).gt.betaOrbs(2)) betaOrbs = [betaOrbs(2),betaOrbs(1)]
        if(alphaOrbs(1).eq.alphaOrbs(2)) then
          alphaVector = [(i,i=1,NDimAlpha)]
        else
          alphaVector = [(i,i=1,alphaOrbs(1)-1),alphaOrbs(2),(i,i=alphaOrbs(1)+1,alphaOrbs(2)-1),&
            alphaOrbs(1),(i,i=alphaOrbs(2)+1,NDimAlpha)]
        endIf
        if(betaOrbs(1).eq.betaOrbs(2)) then
          betaVector = [(i,i=1,NDimBeta)]
        else
          betaVector = [(i,i=1,betaOrbs(1)-1),betaOrbs(2),(i,i=betaOrbs(1)+1,betaOrbs(2)-1),&
            betaOrbs(1),(i,i=betaOrbs(2)+1,NDimBeta)]
        endIf
      endIf

      integralOut = integral%orbitals(alphaOrbsIn=alphaVector,betaOrbsIn=betaVector,axis=my_axis)
!
      end function mqc_integral_swap_orbitals
!
!
!     PROCEDURE MQC_Integral_Combine_Orbitals
!
!>    \brief <b> MQC_Integral_Combine_Orbitals is used to combine two MQC integral 
!>    type variables containing occupied/virtual orbitals into a single MQC integral
!>    type</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Combine_Orbitals is used to combine two MQC integral type variables 
!>    containing occupied/virtual orbitals into a single MQC integral type.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable with combined orbital sets.
!>    \endverbatim
!>    
!>    \param[in] Occupieds
!>    \verbatim
!>        Occupieds is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable with the occupied orbital set.
!>    \endverbatim
!>    
!>    \param[in] Virtuals
!>    \verbatim
!>        Virtuals is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable with the virtual orbital set.
!>    \endverbatim
!>    
!>    \param[in] arrayNameIn
!>    \verbatim
!>        arrayNameIn is character(len=*),intent(in),optional 
!>        The array name label for the MQC_SCF_integrals variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2020
!
      subroutine mqc_integral_combine_orbitals(integral,occupieds,virtuals,arrayNameIn) 
!
      implicit none
      class(mqc_scf_integral),intent(inOut)::integral
      type(mqc_scf_integral),intent(in)::occupieds,virtuals
      character(len=*),intent(in),optional::arrayNameIn
      type(mqc_matrix)::tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,alpha,beta,alphaBeta,betaAlpha
      character(len=64)::arrayName

      if(present(arrayNameIn)) then
        arrayName = arrayNameIn
      else
        arrayName = ''
      endIf

      if(occupieds%type().eq.'space'.and.virtuals%type().eq.'space') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = virtuals%getBlock('alpha')
        if(size(tmp1,1).ne.size(tmp2,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp2,1)',&
          size(tmp2,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp2,2))
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp2,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp2,2)])
        call mqc_integral_allocate(integral,arrayName,'space',alpha)
      elseIf(occupieds%type().eq.'space'.and.virtuals%type().eq.'spin') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = virtuals%getBlock('alpha')
        tmp3 = virtuals%getBlock('beta')
        if(size(tmp1,1).ne.size(tmp2,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp2,1)',&
          size(tmp2,1))
        if(size(tmp1,1).ne.size(tmp3,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp3,1)',&
          size(tmp3,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp2,2))
        call beta%init(size(tmp1,1),size(tmp1,2)+size(tmp3,2))
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp2,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp2,2)])
        call beta%mput(tmp1,[0],[1,size(tmp1,2)])
        call beta%mput(tmp3,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp3,2)])
        call mqc_integral_allocate(integral,arrayName,'spin',alpha,beta)
      elseIf(occupieds%type().eq.'spin'.and.virtuals%type().eq.'space') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = occupieds%getBlock('beta')
        tmp3 = virtuals%getBlock('alpha')
        if(size(tmp1,1).ne.size(tmp3,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp3,1)',&
          size(tmp3,1))
        if(size(tmp2,1).ne.size(tmp3,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp2,1)',size(tmp2,1),'size(tmp3,1)',&
          size(tmp3,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp3,2))
        call beta%init(size(tmp2,1),size(tmp2,2)+size(tmp3,2))
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp3,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp3,2)])
        call beta%mput(tmp2,[0],[1,size(tmp2,2)])
        call beta%mput(tmp3,[0],[size(tmp2,2)+1,size(tmp2,2)+size(tmp3,2)])
        call mqc_integral_allocate(integral,arrayName,'spin',alpha,beta)
      elseIf(occupieds%type().eq.'spin'.and.virtuals%type().eq.'spin') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = occupieds%getBlock('beta')
        tmp3 = virtuals%getBlock('alpha')
        tmp4 = virtuals%getBlock('beta')
        if(size(tmp1,1).ne.size(tmp3,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp3,1)',&
          size(tmp3,1))
        if(size(tmp2,1).ne.size(tmp4,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp2,1)',size(tmp2,1),'size(tmp4,1)',&
          size(tmp4,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp3,2))
        call beta%init(size(tmp2,1),size(tmp2,2)+size(tmp4,2))
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp3,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp3,2)])
        call beta%mput(tmp2,[0],[1,size(tmp2,2)])
        call beta%mput(tmp4,[0],[size(tmp2,2)+1,size(tmp2,2)+size(tmp4,2)])
        call mqc_integral_allocate(integral,arrayName,'spin',alpha,beta)
      elseIf(occupieds%type().eq.'space'.and.virtuals%type().eq.'general') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = virtuals%getBlock('alpha')
        tmp3 = virtuals%getBlock('beta')
        tmp4 = virtuals%getBlock('alpha-beta')
        tmp5 = virtuals%getBlock('beta-alpha')
        if(size(tmp1,1).ne.size(tmp2,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp2,1)',&
          size(tmp2,1))
        if(size(tmp1,1).ne.size(tmp3,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp3,1)',&
          size(tmp3,1))
        if(size(tmp1,1).ne.size(tmp4,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp4,1)',&
          size(tmp4,1))
        if(size(tmp1,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp5,1)',&
          size(tmp5,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp2,2))
        call beta%init(size(tmp1,1),size(tmp1,2)+size(tmp3,2))
        call alphaBeta%init(size(tmp1,1),size(tmp1,2)+size(tmp4,2),0.0)
        call betaAlpha%init(size(tmp1,1),size(tmp1,2)+size(tmp5,2),0.0)
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp2,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp2,2)])
        call beta%mput(tmp1,[0],[1,size(tmp1,2)])
        call beta%mput(tmp3,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp3,2)])
        call alphaBeta%mput(tmp4,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp4,2)])
        call betaAlpha%mput(tmp5,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp5,2)])
        call mqc_integral_allocate(integral,arrayName,'general',alpha,beta,alphaBeta,betaAlpha)
      elseIf(occupieds%type().eq.'general'.and.virtuals%type().eq.'space') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = occupieds%getBlock('beta')
        tmp3 = occupieds%getBlock('alpha-beta')
        tmp4 = occupieds%getBlock('beta-alpha')
        tmp5 = virtuals%getBlock('alpha')
        if(size(tmp1,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp5,1)',&
          size(tmp5,1))
        if(size(tmp2,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp2,1)',size(tmp2,1),'size(tmp5,1)',&
          size(tmp5,1))
        if(size(tmp3,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp3,1)',size(tmp3,1),'size(tmp5,1)',&
          size(tmp5,1))
        if(size(tmp4,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp4,1)',size(tmp4,1),'size(tmp5,1)',&
          size(tmp5,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp5,2))
        call beta%init(size(tmp2,1),size(tmp2,2)+size(tmp5,2))
        call alphaBeta%init(size(tmp3,1),size(tmp3,2)+size(tmp5,2),0.0)
        call betaAlpha%init(size(tmp4,1),size(tmp4,2)+size(tmp5,2),0.0)
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp5,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp5,2)])
        call beta%mput(tmp2,[0],[1,size(tmp2,2)])
        call beta%mput(tmp5,[0],[size(tmp2,2)+1,size(tmp2,2)+size(tmp5,2)])
        call alphaBeta%mput(tmp3,[0],[1,size(tmp3,2)])
        call betaAlpha%mput(tmp4,[0],[1,size(tmp4,2)])
        call mqc_integral_allocate(integral,arrayName,'general',alpha,beta,alphaBeta,betaAlpha)
      elseIf(occupieds%type().eq.'spin'.and.virtuals%type().eq.'general') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = occupieds%getBlock('beta')
        tmp3 = virtuals%getBlock('alpha')
        tmp4 = virtuals%getBlock('beta')
        tmp5 = virtuals%getBlock('alpha-beta')
        tmp6 = virtuals%getBlock('beta-alpha')
        if(size(tmp1,1).ne.size(tmp3,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp3,1)',&
          size(tmp3,1))
        if(size(tmp2,1).ne.size(tmp4,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp2,1)',size(tmp2,1),'size(tmp4,1)',&
          size(tmp4,1))
        if(size(tmp1,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp5,1)',&
          size(tmp5,1))
        if(size(tmp2,1).ne.size(tmp6,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp2,1)',size(tmp2,1),'size(tmp6,1)',&
          size(tmp6,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp3,2))
        call beta%init(size(tmp2,1),size(tmp2,2)+size(tmp4,2))
        call alphaBeta%init(size(tmp1,1),size(tmp1,2)+size(tmp5,2),0.0)
        call betaAlpha%init(size(tmp2,1),size(tmp2,2)+size(tmp6,2),0.0)
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp3,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp3,2)])
        call beta%mput(tmp2,[0],[1,size(tmp2,2)])
        call beta%mput(tmp4,[0],[size(tmp2,2)+1,size(tmp2,2)+size(tmp4,2)])
        call alphaBeta%mput(tmp5,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp5,2)])
        call betaAlpha%mput(tmp6,[0],[size(tmp2,2)+1,size(tmp2,2)+size(tmp6,2)])
        call mqc_integral_allocate(integral,arrayName,'general',alpha,beta,alphaBeta,betaAlpha)
      elseIf(occupieds%type().eq.'general'.and.virtuals%type().eq.'spin') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = occupieds%getBlock('beta')
        tmp3 = occupieds%getBlock('alpha-beta')
        tmp4 = occupieds%getBlock('beta-alpha')
        tmp5 = virtuals%getBlock('alpha')
        tmp6 = virtuals%getBlock('beta')
        if(size(tmp1,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp5,1)',&
          size(tmp5,1))
        if(size(tmp2,1).ne.size(tmp6,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp2,1)',size(tmp2,1),'size(tmp6,1)',&
          size(tmp6,1))
        if(size(tmp3,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp3,1)',size(tmp3,1),'size(tmp5,1)',&
          size(tmp5,1))
        if(size(tmp4,1).ne.size(tmp6,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp4,1)',size(tmp4,1),'size(tmp6,1)',&
          size(tmp6,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp5,2))
        call beta%init(size(tmp2,1),size(tmp2,2)+size(tmp6,2))
        call alphaBeta%init(size(tmp3,1),size(tmp3,2)+size(tmp5,2),0.0)
        call betaAlpha%init(size(tmp4,1),size(tmp4,2)+size(tmp6,2),0.0)
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp5,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp5,2)])
        call beta%mput(tmp2,[0],[1,size(tmp2,2)])
        call beta%mput(tmp6,[0],[size(tmp2,2)+1,size(tmp2,2)+size(tmp6,2)])
        call alphaBeta%mput(tmp3,[0],[1,size(tmp3,2)])
        call betaAlpha%mput(tmp4,[0],[1,size(tmp4,2)])
        call mqc_integral_allocate(integral,arrayName,'general',alpha,beta,alphaBeta,betaAlpha)
      elseIf(occupieds%type().eq.'general'.and.virtuals%type().eq.'general') then
        tmp1 = occupieds%getBlock('alpha')
        tmp2 = occupieds%getBlock('beta')
        tmp3 = occupieds%getBlock('alpha-beta')
        tmp4 = occupieds%getBlock('beta-alpha')
        tmp5 = virtuals%getBlock('alpha')
        tmp6 = virtuals%getBlock('beta')
        tmp7 = virtuals%getBlock('alpha-beta')
        tmp8 = virtuals%getBlock('beta-alpha')
        if(size(tmp1,1).ne.size(tmp5,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp1,1)',size(tmp1,1),'size(tmp5,1)',&
          size(tmp5,1))
        if(size(tmp2,1).ne.size(tmp6,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp2,1)',size(tmp2,1),'size(tmp6,1)',&
          size(tmp6,1))
        if(size(tmp3,1).ne.size(tmp7,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp3,1)',size(tmp3,1),'size(tmp7,1)',&
          size(tmp7,1))
        if(size(tmp4,1).ne.size(tmp8,1)) call mqc_error_i('AO basis size does not match in &
          &mqc_integral_combine_orbitals',6,'size(tmp4,1)',size(tmp4,1),'size(tmp8,1)',&
          size(tmp8,1))
        call alpha%init(size(tmp1,1),size(tmp1,2)+size(tmp5,2))
        call beta%init(size(tmp2,1),size(tmp2,2)+size(tmp6,2))
        call alphaBeta%init(size(tmp3,1),size(tmp3,2)+size(tmp7,2),0.0)
        call betaAlpha%init(size(tmp4,1),size(tmp4,2)+size(tmp8,2),0.0)
        call alpha%mput(tmp1,[0],[1,size(tmp1,2)])
        call alpha%mput(tmp5,[0],[size(tmp1,2)+1,size(tmp1,2)+size(tmp5,2)])
        call beta%mput(tmp2,[0],[1,size(tmp2,2)])
        call beta%mput(tmp6,[0],[size(tmp2,2)+1,size(tmp2,2)+size(tmp6,2)])
        call alphaBeta%mput(tmp3,[0],[1,size(tmp3,2)])
        call alphaBeta%mput(tmp7,[0],[size(tmp3,2)+1,size(tmp3,2)+size(tmp7,2)])
        call betaAlpha%mput(tmp4,[0],[1,size(tmp4,2)])
        call betaAlpha%mput(tmp8,[0],[size(tmp4,2)+1,size(tmp4,2)+size(tmp8,2)])
        call mqc_integral_allocate(integral,arrayName,'general',alpha,beta,alphaBeta,betaAlpha)
      endIf
!
      end subroutine mqc_integral_combine_orbitals
!
!
!     PROCEDURE MQC_Integral_Swap_OffDiagBlocks
!
!>    \brief <b> MQC_Integral_Swap_OffDiagBlocks is used to exchange the alpha-beta and 
!>    beta-alpha blocks of an MQC integral type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Swap_OffDiagBlocks is used to exchange the alpha-beta and 
!>    beta-alpha blocks of an MQC integral type variable. This is required for performing
!>    the PG(P) contraction with non-symmetric GHF densities (most commonly with complex GHF
!>    or transition density matrices).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        integral is class(mqc_scf_integral)
!>        The MQC_SCF_integrals variable to swap blocks.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2023
!
      function mqc_integral_swap_offDiagBlocks(integral) result(integralOut)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmat1,tmat2,tmat3,tmat4

      if(integral%array_type.eq.'general') then
        tmat1 = integral%getBlock('alpha-alpha') 
        tmat2 = integral%getBlock('beta-beta')
        tmat3 = integral%getBlock('alpha-beta')
        tmat4 = integral%getBlock('beta-alpha')
        call mqc_integral_allocate(integralOut,integral%getLabel(),'general',tmat1,tmat2,tmat4,tmat3)
      else
        integralOut = integral
      endIf
!
      end function mqc_integral_swap_offDiagBlocks
!
!
!     PROCEDURE MQC_Bit2Num_String
!
!>    \brief <b> MQC_Bit2Num_String is used to convert a bit string implying occupied
!>    and virtual orbitals, and convert it to an integer vector with the numbers of 
!>    occupied orbitals</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Bit2Num_String is used to convert a bit string implying occupied
!>    and virtual orbitals, and convert it to an integer vector with the numbers of 
!>    occupied orbitals.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] BitString
!>    \verbatim
!>        BitString is Type(MQC_Vector)
!>        The bit string to be converted.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2022
!
      function mqc_bit2num_string(bitString) result(numString)

      implicit none
      type(mqc_vector),intent(in)::bitString
      integer(kind=int64),dimension(:),allocatable::numString

      integer::nBitInts,temp_int,i,j
      type(mqc_vector)::temp_vec

      nBitInts = size(bitString)
      do i = 1, nBitInts
        temp_int = bitString%at(i)
        do j = 0, Bit_Size(0)-1
          if(btest(temp_int,j)) then
            call temp_vec%push((i-1)*Bit_Size(0)+j+1)
          endIf
        endDo
      endDo
      numString = temp_vec

      end function mqc_bit2num_string
!
!
!     PROCEDURE MQC_DetString_Print
!
!>    \brief <b> MQC_DetString_Print takes as an argument an alpha bit string and a 
!>    beta bit string and returns a string for printing</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_DetString_Print takes as an argument an alpha bit string and a 
!>    beta bit string and returns a string for printing.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] AlphaString
!>    \verbatim
!>        AlphaString is Type(MQC_Vector)
!>        The alpha bit string to be printed.
!>    \endverbatim
!>
!>    \param[in] BetaString
!>    \verbatim
!>        BetaString is Type(MQC_Vector)
!>        The beta bit string to be printed.
!>    \endverbatim
!>
!>    \param[in] nBasisIn
!>    \verbatim
!>        nBasisIn is Integer, optional
!>        Number of basis functions.
!>    \endverbatim
!>
!>    \param[in] nCoreIn
!>    \verbatim
!>        nCoreIn is Integer, optional
!>        Number of low index orbitals to 
!>        exclude from printing.
!>    \endverbatim
!>
!>    \param[in] nVirtIn
!>    \verbatim
!>        nVirtIn is Integer, optional
!>        Number of high index orbitals to 
!>        exclude from printing.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2022
!
      function mqc_detString_print(alphaString,betaString,nBasisIn,nCoreIn,nVirtIn) result(returnString)

      implicit none
      type(mqc_vector),intent(in)::alphaString,betaString
      integer,optional,intent(in)::nBasisIn,nCoreIn,nVirtIn
      character(:),allocatable::returnString

      integer::nBitInts,temp_alpha_int,temp_beta_int,i,j,nBasis,nCore,nVirt

      if(size(alphaString).ne.size(betaString)) call mqc_error_i('Alpha and beta strings are different&
        & lengths in mqc_detString_print',6,'size(alphaString)',size(alphaString),'size(betaString)',&
        size(betaString))
      nBitInts = size(alphaString)

      if(present(nBasisIn)) then
        nBasis = nBasisIn
      else
        nBasis = nBitInts*Bit_Size(0)
      endIf
      if(present(nCoreIn)) then
        nCore = nCoreIn
      else
        nCore = 0
      endIf
      if(present(nVirtIn)) then
        nVirt = nVirtIn
      else
        nVirt = 0
      endIf

      returnString = '|'
      do i = 1, nBitInts
        temp_alpha_int = alphaString%at(i)
        temp_beta_int = betaString%at(i)
        do j = 0, Bit_Size(0)-1
          if((i-1)*Bit_Size(0)+(j+1).le.nCore.or.(i-1)*Bit_Size(0)+(j+1).gt.nBasis-nVirt) cycle
          if(btest(temp_alpha_int,j).and.btest(temp_beta_int,j)) then
            returnString = trim(returnString)//'1'
          elseIf(btest(temp_alpha_int,j)) then
            returnString = trim(returnString)//'a'
          elseIf(btest(temp_beta_int,j)) then
            returnString = trim(returnString)//'b'
          else
            returnString = trim(returnString)//'0'
          endIf
        endDo
      endDo
      returnString = trim(returnString)//'>'

      end function mqc_detstring_print
!
!     PROCEDURE MQC_Eigenvalues_Output_Block
!
!>    \brief <b> MQC_Eigenvalues_Output_Block is used to output a spin block of an MQC 
!>    eigenvalues type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Output_Block is used to output a spin block of an MQC 
!>    eigenvalues type variable. Optional argument blockName is used to specify which 
!>    spin block should be output, with the default to output the all spin blocks. 
!>    The following options are available:
!>
!>    1. BlockName = 'full' outputs the full MQC_SCF_Integral. 
!>    2. BlockName = 'alpha' outputs the alpha spin block. 
!>    3. BlockName = 'beta' outputs the beta spin block. 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_SCF_Eigenvalues)
!>        The MQC_SCF_Eigenvalues variable to output.
!>    \endverbatim
!>    
!>    \param[in] BlockName
!>    \verbatim
!>        BlockName is character(len=*),optional
!>        = 'full':        outputs the full MQC_SCF_Integral. 
!>        = 'alpha':       outputs the alpha spin block.
!>        = 'beta':        outputs the beta spin block. 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      function mqc_eigenvalues_output_block(eigenvalues,blockName) result(vectorOut)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      character(len=*),optional::blockName
      character(len=64)::myBlockName
      type(mqc_vector)::vectorOut
      integer(kind=int64)::nDimAlpha=0,nDimBeta=0,nDimTotal=0,zeroI
      real(kind=real64)::zeroR
      complex(kind=real64)::zeroC
      parameter(zeroI=0,zeroR=0.0,zeroC=(0.0,0.0))
!
      if(present(blockName)) then
        call string_change_case(blockName,'l',myBlockName)
      else
        myBlockName = 'full'
      endIf

      nDimAlpha = eigenvalues%blockSize('Alpha') 
      nDimBeta = eigenvalues%blockSize('Beta')
      nDimTotal = nDimAlpha + nDimBeta

      select case (myBlockName)
      case('full')
        if (eigenvalues%type().eq.'space') then
          if (eigenvalues%hasAlpha()) then
            if(mqc_vector_haveInteger(eigenvalues%alpha)) then
              call vectorOut%init(nDimAlpha,zeroI)
            elseif(mqc_vector_haveReal(eigenvalues%alpha)) then
              call vectorOut%init(nDimAlpha,zeroR)
            elseIf(mqc_vector_haveComplex(eigenvalues%alpha)) then
              call vectorOut%init(nDimAlpha,zeroC)
            else
              call mqc_error('unrecognized integral data type in mqc_integral_output_block')
            endIf
            call vectorOut%vput(eigenvalues%alpha,1) 
          else
            call mqc_error_L('block does not exist in mqc_eigenvalues_output_block', 6, &
                 'eigenvalues%hasAlpha()', eigenvalues%hasAlpha() )
          endIf
        elseIf ((eigenvalues%type().eq.'spin').or.(eigenvalues%type().eq.'general')) then
          if (eigenvalues%hasAlpha().and.eigenvalues%hasBeta()) then
            if(mqc_vector_haveInteger(eigenvalues%alpha).or.mqc_vector_haveInteger(eigenvalues%beta)) then
              call vectorOut%init(nDimTotal,zeroI)
            elseif(mqc_vector_haveReal(eigenvalues%alpha).or.mqc_vector_haveReal(eigenvalues%beta)) then
              call vectorOut%init(nDimTotal,zeroR)
            elseIf(mqc_vector_haveComplex(eigenvalues%alpha).or.mqc_vector_haveComplex(eigenvalues%beta)) then
              call vectorOut%init(nDimTotal,zeroC)
            else
              call mqc_error('unrecognized integral data type in mqc_integral_output_block')
            endIf
            call vectorOut%vput(eigenvalues%alpha,1) 
            call vectorOut%vput(eigenvalues%beta,nDimAlpha+1)
          else
            call mqc_error_L('block does not exist in mqc_eigenvalues_output_block', 6, &
                 'eigenvalues%hasAlpha()', eigenvalues%hasAlpha(),'eigenvalues%hasBeta()', & 
                 eigenvalues%hasBeta())
          endIf
        endIf
      case('alpha')
        if (eigenvalues%hasAlpha()) then
          vectorOut = eigenvalues%alpha
        else
          call mqc_error_L('block does not exist in mqc_eigenvalues_output_block', 6, &
               'eigenvalues%hasAlpha()', eigenvalues%hasAlpha() )
        endIf
      case('beta')
        if (eigenvalues%hasBeta()) then
          vectorOut = eigenvalues%beta
        elseIf (eigenvalues%type().eq.'space') then
          if (eigenvalues%hasAlpha()) then
            vectorOut = eigenvalues%alpha
          else
            call mqc_error_L('block does not exist in mqc_eigenvalues_output_block', 6, &
                 'eigenvalues%hasAlpha()', eigenvalues%hasAlpha() )
          endIf
        else
          call mqc_error_A('block does not exist in mqc_eigenvalues_output_block', 6, &
               'eigenvalues%type()', eigenvalues%type() )
        endIf
      case default
        call mqc_error_A('unrecognised block name in mqc_eigenvalues_output_block', 6, &
             'myBlockName', myBlockName )
      end select
!
      end function mqc_eigenvalues_output_block  
!
!
!     PROCEDURE MQC_2ERIs_Output_Block
!
!>    \brief <b> MQC_2ERIs_Output_Block is used to output a spin block of an MQC 
!>    2ERIs type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_2ERIs_Output_Block is used to output a spin block of an MQC 2ERIs type 
!>    variable. Optional argument blockName is used to specify which spin block
!>    should be output, with the default to output the all spin blocks. The following
!>    options are available:
!>
!>    1.  BlockName = 'full' outputs the full MQC_TwoERIs variable.
!>    2.  BlockName = 'alpha', 'alpha-alpha' or 'aaaa' outputs the alpha-alpha spin block. 
!>    3.  BlockName = 'beta', 'beta-beta' or 'bbbb' outputs the beta-beta spin block. 
!>    4.  BlockName = 'alpha-beta' or 'aabb' outputs the alpha-beta spin block. 
!>    5.  BlockName = 'beta-alpha' or 'bbaa' outputs the beta-alpha spin block. 
!>    6.  BlockName = 'abab' outputs the (ab|ab) spin block. 
!>    7.  BlockName = 'abba' outputs the (ab|ba) spin block. 
!>    8.  BlockName = 'baab' outputs the (ba|ab) spin block. 
!>    9.  BlockName = 'baba' outputs the (ba|ba) spin block. 
!>    10. BlockName = 'aaab' outputs the (aa|ab) spin block. 
!>    11. BlockName = 'aaba' outputs the (aa|ba) spin block. 
!>    12. BlockName = 'abaa' outputs the (ab|aa) spin block. 
!>    13. BlockName = 'baaa' outputs the (ba|aa) spin block. 
!>    14. BlockName = 'bbba' outputs the (bb|ba) spin block. 
!>    15. BlockName = 'bbab' outputs the (bb|ab) spin block. 
!>    16. BlockName = 'babb' outputs the (ba|bb) spin block. 
!>    17. BlockName = 'abbb' outputs the (ab|bb) spin block. 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] ERIs
!>    \verbatim
!>        eris is class(mqc_twoeris)
!>        The MQC_TwoERIs variable to output.
!>    \endverbatim
!>    
!>    \param[in,out] BlockName
!>    \verbatim
!>        BlockName is character(len=*),optional
!>        = 'full':        outputs the full MQC_SCF_Integral. 
!>        = 'alpha':       outputs the alpha-alpha spin block.
!>        = 'alpha-alpha': outputs the alpha-alpha spin block. 
!>        = 'aaaa':        outputs the alpha-alpha spin block. 
!>        = 'beta':        outputs the beta-beta spin block. 
!>        = 'beta-beta':   outputs the beta-beta spin block. 
!>        = 'bbbb':        outputs the beta-beta spin block. 
!>        = 'alpha-beta'   outputs the alpha-beta spin block. 
!>        = 'aabb':        outputs the alpha-beta spin block. 
!>        = 'beta-alpha'   outputs the beta-alpha spin block. 
!>        = 'bbaa':        outputs the beta-alpha spin block. 
!>        = 'abab':        outputs the (ab|ab) spin block. 
!>        = 'abba':        outputs the (ab|ba) spin block. 
!>        = 'baab':        outputs the (ba|ab) spin block. 
!>        = 'baba':        outputs the (ba|ba) spin block. 
!>        = 'aaab':        outputs the (aa|ab) spin block. 
!>        = 'aaba':        outputs the (aa|ba) spin block. 
!>        = 'abaa':        outputs the (ab|aa) spin block. 
!>        = 'baaa':        outputs the (ba|aa) spin block. 
!>        = 'bbba':        outputs the (bb|ba) spin block. 
!>        = 'bbab':        outputs the (bb|ab) spin block. 
!>        = 'babb':        outputs the (ba|bb) spin block. 
!>        = 'abbb':        outputs the (ab|bb) spin block. 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2019
!
      function mqc_2ERIs_output_block(eris,blockName) result(r4tOut)
!
      implicit none
      class(mqc_twoeris)::eris
      character(len=*),optional::blockName
      character(len=64)::myBlockName
      type(mqc_r4Tensor)::r4tOut
      integer(kind=int64)::nDimAlpha1=0,nDimBeta1=0,nDimTotal1=0,nDimAlpha2=0,nDimBeta2=0, &
        nDimTotal2=0,nDimAlpha3=0,nDimBeta3=0,nDimTotal3=0,nDimAlpha4=0,nDimBeta4=0, &
        nDimTotal4=0,zeroI
      real(kind=real64)::zeroR
      complex(kind=real64)::zeroC
      parameter(zeroI=0,zeroR=0.0,zeroC=(0.0,0.0))
!
      if(present(blockName)) then
        call string_change_case(blockName,'l',myBlockName)
      else
        myBlockName = 'full'
      endIf

      nDimAlpha1 = eris%blockSize('Alpha',1) 
      nDimBeta1 = eris%blockSize('Beta',1)
      nDimTotal1 = nDimAlpha1 + nDimBeta1
      nDimAlpha2 = eris%blockSize('Alpha',2) 
      nDimBeta2 = eris%blockSize('Beta',2)
      nDimTotal2 = nDimAlpha2 + nDimBeta2
      nDimAlpha3 = eris%blockSize('Alpha',3) 
      nDimBeta3 = eris%blockSize('Beta',3)
      nDimTotal3 = nDimAlpha3 + nDimBeta3
      nDimAlpha4 = eris%blockSize('Alpha',4) 
      nDimBeta4 = eris%blockSize('Beta',4)
      nDimTotal4 = nDimAlpha4 + nDimBeta4

      select case (myBlockName)
      case('full')
        if (eris%type().eq.'regular'.or.eris%type().eq.'raffenetti1'.or. &
          eris%type().eq.'raffenetti2'.or.eris%type().eq.'raffenetti3'.or. &
          eris%type().eq.'space') then
          call r4tOut%init(nDimAlpha1,nDimAlpha2,nDimAlpha3,nDimAlpha4)
          if(eris%hasSpinBlock('aaaa')) call r4tOut%tput(eris%alpha,[1,nDimAlpha1],[1,nDimAlpha2], &
            [1,nDimAlpha3],[1,nDimAlpha4]) 
        elseIf (eris%type().eq.'spin') then
          call r4tOut%init(nDimTotal1,nDimTotal2,nDimTotal3,nDimTotal4)
          if(eris%hasSpinBlock('aaaa')) call r4tOut%tput(eris%alpha,[1,nDimAlpha1],[1,nDimAlpha2], &
            [1,nDimAlpha3],[1,nDimAlpha4]) 
          if(eris%hasSpinBlock('bbbb')) call r4tOut%tput(eris%beta,[nDimAlpha1+1,nDimTotal1],&
            [nDimAlpha2+1,nDimTotal2],[nDimAlpha3+1,nDimTotal3],[nDimAlpha4+1,nDimTotal4]) 
          if(eris%hasSpinBlock('aabb')) call r4tOut%tput(eris%alphaBeta,[1,nDimAlpha1],[1,nDimAlpha2], &
            [nDimAlpha3+1,nDimTotal3],[nDimAlpha4+1,nDimTotal4]) 
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('bbaa')) call r4tOut%tput(eris%betaAlpha,[nDimAlpha1+1,nDimTotal1], &
              [nDimAlpha2+1,nDimTotal2],[1,nDimAlpha3],[1,nDimAlpha4]) 
          else
            if(eris%hasSpinBlock('aabb')) call r4tOut%tput(eris%alphaBeta,[nDimAlpha1+1,nDimTotal1], &
              [nDimAlpha2+1,nDimTotal2],[1,nDimAlpha3],[1,nDimAlpha4]) 
          endIf
        elseIf (eris%type().eq.'general') then
          call r4tOut%init(nDimTotal1,nDimTotal2,nDimTotal3,nDimTotal4)
          if(eris%hasSpinBlock('aaaa')) call r4tOut%tput(eris%alpha,[1,nDimAlpha1],[1,nDimAlpha2],&
            [1,nDimAlpha3],[1,nDimAlpha4]) 
          if(eris%hasSpinBlock('bbbb')) call r4tOut%tput(eris%beta,[nDimAlpha1+1,nDimTotal1],&
            [nDimAlpha2+1,nDimTotal2],[nDimAlpha3+1,nDimTotal3],[nDimAlpha4+1,nDimTotal4]) 
          if(eris%hasSpinBlock('aabb')) call r4tOut%tput(eris%alphaBeta,[1,nDimAlpha1],[1,nDimAlpha2], &
            [nDimAlpha3+1,nDimTotal3],[nDimAlpha4+1,nDimTotal4]) 
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('bbaa')) call r4tOut%tput(eris%betaAlpha,[nDimAlpha1+1,nDimTotal1], &
              [nDimAlpha2+1,nDimTotal2],[1,nDimAlpha3],[1,nDimAlpha4]) 
          else
            if(eris%hasSpinBlock('aabb')) call r4tOut%tput(eris%alphaBeta,[nDimAlpha1+1,nDimTotal1], &
              [nDimAlpha2+1,nDimTotal2],[1,nDimAlpha3],[1,nDimAlpha4]) 
          endIf
          if(eris%hasSpinBlock('abab')) call r4tOut%tput(eris%abab,[1,nDimAlpha1],[nDimAlpha2+1,nDimTotal2], &
              [1,nDimAlpha3], [nDimAlpha4+1,nDimTotal4])
          if(eris%hasSpinBlock('abba')) call r4tOut%tput(eris%abba,[1,nDimAlpha1],[nDimAlpha2+1,nDimTotal2], &
              [nDimAlpha3+1,nDimTotal3],[1,nDimAlpha4])
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('baab')) call r4tOut%tput(eris%baab,[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2], &
              [1,nDimAlpha3],[nDimAlpha4+1,nDimTotal4])
          else
            if(eris%hasSpinBlock('abba')) call r4tOut%tput(eris%abba,[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2], &
              [1,nDimAlpha3],[nDimAlpha4+1,nDimTotal4])
          endIf
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('baba')) call r4tOut%tput(eris%baba,[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2], &
              [nDimAlpha3+1,nDimTotal3],[1,nDimAlpha4])
          else
            if(eris%hasSpinBlock('abab')) call r4tOut%tput(conjg(eris%baba),[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2], &
              [nDimAlpha3+1,nDimTotal3],[1,nDimAlpha4])
          endIf

          if(eris%hasSpinBlock('aaab')) call r4tOut%tput(eris%aaab,[1,nDimAlpha1],[1,nDimAlpha2], &
              [1,nDimAlpha3],[nDimAlpha4+1,nDimTotal4])
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('aaba')) call r4tOut%tput(eris%aaba,[1,nDimAlpha1],[1,nDimAlpha2], &
              [nDimAlpha3+1,nDimTotal3],[1,nDimAlpha4])
          else
            if(eris%hasSpinBlock('aaab')) call r4tOut%tput(conjg(eris%aaba),[1,nDimAlpha1],[1,nDimAlpha2], &
              [nDimAlpha3+1,nDimTotal3],[1,nDimAlpha4])
          endIf
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('abaa')) call r4tOut%tput(eris%abaa,[1,nDimAlpha1],[nDimAlpha2+1,nDimTotal2], &
              [1,nDimAlpha3],[1,nDimAlpha4])
          else
            if(eris%hasSpinBlock('aaab')) call r4tOut%tput(eris%aaab,[1,nDimAlpha1],[nDimAlpha2+1,nDimTotal2], &
              [1,nDimAlpha3],[1,nDimAlpha4])
          endIf
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('baaa')) call r4tOut%tput(eris%baaa,[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2], &
              [1,nDimAlpha3],[1,nDimAlpha4])
          else
            if(eris%hasSpinBlock('aaab')) call r4tOut%tput(conjg(eris%aaab),[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2], &
              [1,nDimAlpha3],[1,nDimAlpha4])
          endIf
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('bbba')) call r4tOut%tput(eris%bbba,[nDimAlpha1+1,nDimTotal1], &
              [nDimAlpha2+1,nDimTotal2],[nDimAlpha3+1,nDimTotal3],[1,nDimAlpha4])
          else
            if(eris%hasSpinBlock('abbb')) call r4tOut%tput(conjg(eris%abbb),[nDimAlpha1+1,nDimTotal1], &
              [nDimAlpha2+1,nDimTotal2],[nDimAlpha3+1,nDimTotal3],[1,nDimAlpha4])
          endIf
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('bbab')) call r4tOut%tput(eris%bbab,[nDimAlpha1+1,nDimTotal1], &
              [nDimAlpha2+1,nDimTotal2],[1,nDimAlpha3],[nDimAlpha4+1,nDimTotal4])
          else
            if(eris%hasSpinBlock('abbb')) call r4tOut%tput(eris%abbb,[nDimAlpha1+1,nDimTotal1], &
              [nDimAlpha2+1,nDimTotal2],[1,nDimAlpha3],[nDimAlpha4+1,nDimTotal4])
          endIf
          if(eris%storageType.eq.'full') then
            if(eris%hasSpinBlock('babb')) call r4tOut%tput(eris%babb,[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2], &
              [nDimAlpha3+1,nDimTotal3],[nDimAlpha4+1,nDimTotal4])
          else
            if(eris%hasSpinBlock('abbb')) call r4tOut%tput(conjg(eris%abbb),[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2], &
              [nDimAlpha3+1,nDimTotal3],[nDimAlpha4+1,nDimTotal4])
          endIf
          if(eris%hasSpinBlock('abbb')) call r4tOut%tput(eris%abbb,[1,nDimAlpha1],[nDimAlpha2+1,nDimTotal2], &
            [nDimAlpha3+1,nDimTotal3],[nDimAlpha4+1,nDimTotal4])
        endIf
      case('alpha','aaaa','alpha-alpha')
        if (eris%hasSpinBlock('aaaa')) then
          r4tOut = eris%alpha
        else
          call r4tOut%init(NDimAlpha1,NDimAlpha2,NDimAlpha3,NDimAlpha4,zeroR)
        endIf
      case('beta','bbbb','beta-beta')
        if (eris%hasSpinBlock('bbbb')) then
          r4tOut = eris%beta
        elseIf (eris%type().eq.'regular'.or.eris%type().eq.'raffenetti1'.or. &
          eris%type().eq.'raffenetti2'.or.eris%type().eq.'raffenetti3'.or. &
          eris%type().eq.'space') then
          if (eris%hasSpinBlock('aaaa')) then
            r4tOut = eris%alpha
          else
            call r4tOut%init(NDimAlpha1,NDimAlpha2,NDimAlpha3,NDimAlpha4,zeroR)
          endIf
        else
          call r4tOut%init(NDimBeta1,NDimBeta2,NDimBeta3,NDimBeta4,zeroR)
        endIf
      case('alpha-beta','aabb')
        if (eris%hasSpinBlock('aabb')) then
          r4tOut = eris%alphaBeta
        elseIf (eris%type().eq.'regular'.or.eris%type().eq.'raffenetti1'.or. &
          eris%type().eq.'raffenetti2'.or.eris%type().eq.'raffenetti3'.or. &
          eris%type().eq.'space') then
          if (eris%hasSpinBlock('aaaa')) then
            r4tOut = eris%alpha
          else
            call r4tOut%init(NDimAlpha1,NDimAlpha2,NDimBeta3,NDimBeta4,zeroR)
          endIf
        else
          call r4tOut%init(NDimAlpha1,NDimAlpha2,NDimBeta3,NDimBeta4,zeroR)
        endIf
      case('beta-alpha','bbaa')
        if (eris%hasSpinBlock('bbaa')) then
          r4tOut = eris%betaAlpha
        elseif (eris%type().eq.'regular'.or.eris%type().eq.'raffenetti1'.or. &
          eris%type().eq.'raffenetti2'.or.eris%type().eq.'raffenetti3'.or. &
          eris%type().eq.'space') then
          if (eris%hasSpinBlock('aaaa')) then
            r4tOut = eris%alpha
          else
            call r4tOut%init(NDimBeta1,NDimBeta2,NDimAlpha3,NDimAlpha4,zeroR)
          endIf
        elseIf((eris%type().eq.'spin'.or.eris%type().eq.'general').and.eris%storageType.eq.'symm') then
          if (eris%hasSpinBlock('aabb')) then
            r4tOut = eris%alphaBeta
          else
            call r4tOut%init(NDimBeta1,NDimBeta2,NDimAlpha3,NDimAlpha4,zeroR)
          endIf
        else
          call r4tOut%init(NDimBeta1,NDimBeta2,NDimAlpha3,NDimAlpha4,zeroR)
        endIf
      case('abab')
        if (eris%hasSpinBlock('abab')) then
          r4tOut = eris%abab
        else
          call r4tOut%init(NDimAlpha1,NDimBeta2,NDimAlpha3,NDimBeta4,zeroR)
        endIf
      case('abba')
        if (eris%hasSpinBlock('abba')) then
          r4tOut = eris%abba
        else
          call r4tOut%init(NDimAlpha1,NDimBeta2,NDimBeta3,NDimAlpha4,zeroR)
        endIf
      case('baab')
        if (eris%hasSpinBlock('baab')) then
          r4tOut = eris%baab
        else
          if(eris%type().eq.'general'.and.eris%storageType.eq.'symm') then
            if (eris%hasSpinBlock('abba')) then
              r4tOut = eris%abba
            else
              call r4tOut%init(NDimbeta1,NDimAlpha2,NDimAlpha3,NDimBeta4,zeroR)
            endIf
          else
            call r4tOut%init(NDimbeta1,NDimAlpha2,NDimAlpha3,NDimBeta4,zeroR)
          endIf
        endIf
      case('baba')
        if (eris%hasSpinBlock('baba')) then
          r4tOut = eris%baba
        else
          if(eris%type().eq.'general'.and.eris%storageType.eq.'symm') then
            if (eris%hasSpinBlock('abab')) then
              r4tOut = conjg(eris%abab)
            else
              call r4tOut%init(NDimBeta1,NDimAlpha2,NDimBeta3,NDimAlpha4,zeroR)
            endIf
          else
            call r4tOut%init(NDimBeta1,NDimAlpha2,NDimBeta3,NDimAlpha4,zeroR)
          endIf
        endIf
      case('aaab')
        if (eris%hasSpinBlock('aaab')) then
          r4tOut = eris%aaab
        else
          call r4tOut%init(NDimAlpha1,NDimAlpha2,NDimAlpha3,NDimBeta4,zeroR)
        endIf
      case('aaba')
        if (eris%hasSpinBlock('aaba')) then
          r4tOut = eris%aaba
        else
          if(eris%type().eq.'general'.and.eris%storageType.eq.'symm') then
            if (eris%hasSpinBlock('aaab')) then
              r4tOut = conjg(eris%aaab)
            else
              call r4tOut%init(NDimAlpha1,NDimAlpha2,NDimBeta3,NDimAlpha4,zeroR)
            endIf
          else
            call r4tOut%init(NDimAlpha1,NDimAlpha2,NDimBeta3,NDimAlpha4,zeroR)
          endIf
        endIf
      case('abaa')
        if (eris%hasSpinBlock('abaa')) then
          r4tOut = eris%abaa
        else
          if(eris%type().eq.'general'.and.eris%storageType.eq.'symm') then
            if (eris%hasSpinBlock('aaab')) then
              r4tOut = eris%aaab
            else
              call r4tOut%init(NDimAlpha1,NDimBeta2,NDimAlpha3,NDimAlpha4,zeroR)
            endIf
          else
            call r4tOut%init(NDimAlpha1,NDimBeta2,NDimAlpha3,NDimAlpha4,zeroR)
          endIf
        endIf
      case('baaa')
        if (eris%hasSpinBlock('baaa')) then
          r4tOut = eris%baaa
        else
          if(eris%type().eq.'general'.and.eris%storageType.eq.'symm') then
            if (eris%hasSpinBlock('aaab')) then
              r4tOut = conjg(eris%aaab)
            else
              call r4tOut%init(NDimBeta1,NDimAlpha2,NDimAlpha3,NDimAlpha4,zeroR)
            endIf
          else
            call r4tOut%init(NDimBeta1,NDimAlpha2,NDimAlpha3,NDimAlpha4,zeroR)
          endIf
        endIf
      case('bbba')
        if (eris%hasSpinBlock('bbba')) then
          r4tOut = eris%bbba
        else
          if(eris%type().eq.'general'.and.eris%storageType.eq.'symm') then
            if (eris%hasSpinBlock('abbb')) then
              r4tOut = conjg(eris%abbb)
            else
              call r4tOut%init(NDimBeta1,NDimBeta2,NDimBeta3,NDimAlpha4,zeroR)
            endIf
          else
            call r4tOut%init(NDimBeta1,NDimBeta2,NDimBeta3,NDimAlpha4,zeroR)
          endIf
        endIf
      case('bbab')
        if (eris%hasSpinBlock('bbab')) then
          r4tOut = eris%bbab
        else
          if(eris%type().eq.'general'.and.eris%storageType.eq.'symm') then
            if (eris%hasSpinBlock('abbb')) then
              r4tOut = eris%abbb
            else
              call r4tOut%init(NDimBeta1,NDimBeta2,NDimAlpha3,NDimBeta4,zeroR)
            endIf
          else
            call r4tOut%init(NDimBeta1,NDimBeta2,NDimAlpha3,NDimBeta4,zeroR)
          endIf
        endIf
      case('babb')
        if (eris%hasSpinBlock('babb')) then
          r4tOut = eris%babb
        else
          if(eris%type().eq.'general'.and.eris%storageType.eq.'symm') then
            if (eris%hasSpinBlock('abbb')) then
              r4tOut = conjg(eris%abbb)
            else
              call r4tOut%init(NDimBeta1,NDimAlpha2,NDimBeta3,NDimBeta4,zeroR)
            endIf
          else
            call r4tOut%init(NDimBeta1,NDimAlpha2,NDimBeta3,NDimBeta4,zeroR)
          endIf
        endIf
      case('abbb')
        if (eris%hasSpinBlock('abbb')) then
          r4tOut = eris%abbb
        else
          call r4tOut%init(NDimAlpha1,NDimBeta2,NDimBeta3,NDimBeta4,zeroR)
        endIf
      case default
        call mqc_error_A('unrecognised block name in mqc_eris_output_block', 6, &
             'myBlockName', myBlockName )
      end select
!
      end function mqc_2ERIs_output_block  
!
!
!     PROCEDURE MQC_Integral_Output_Array
!
!>    \brief <b> MQC_Integral_Output_Array is used to output an MQC integral type
!>    variable to an MQC matrix</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Output_Array is used to output an MQC integral type variable to an 
!>    MQC matrix.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] MatrixOut 
!>    \verbatim
!>        MatrixOut is type(mqc_matrix)
!>        The MQC matrix variable to output.
!>    \endverbatim
!>    
!>    \param[in] IntegralIn
!>    \verbatim
!>        IntegralIn is Class(MQC_SCF_Integral)
!>        MQC integral variable to output as MQC matrix.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      subroutine mqc_integral_output_array(matrixOut,integralIn) 
!
      implicit none
      class(mqc_scf_integral),intent(in)::integralIn
      type(mqc_matrix),intent(inOut)::matrixOut
      integer(kind=int64)::nDimAlpha1=0,nDimBeta1=0,nDimTotal1=0,nDimAlpha2=0,nDimBeta2=0,nDimTotal2=0
!
      nDimAlpha1 = integralIn%blockSize('Alpha',1) 
      nDimBeta1 = integralIn%blockSize('Beta',1)
      nDimTotal1 = nDimAlpha1 + nDimBeta1
      nDimAlpha2 = integralIn%blockSize('Alpha',2) 
      nDimBeta2 = integralIn%blockSize('Beta',2)
      nDimTotal2 = nDimAlpha2 + nDimBeta2

      select case (integralIn%type())
      case('space')
        if(mqc_have_complex(integralIn%alpha)) then
          call matrixOut%init(nDimAlpha1,nDimAlpha2,cmplx(0.0,0.0))
        else
          call matrixOut%init(nDimAlpha1,nDimAlpha2)
        endIf
        if (integralIn%hasAlpha()) then
          call matrixOut%mput(integralIn%alpha,[1,nDimAlpha1],[1,nDimAlpha2]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasAlpha()', integralIn%hasAlpha() )
        endIf
      case('spin')
        if(mqc_have_complex(integralIn%alpha)) then
          call matrixOut%init(nDimTotal1,nDimTotal2,cmplx(0.0,0.0))
        else
          call matrixOut%init(nDimTotal1,nDimTotal2)
        endIf
        if (integralIn%hasAlpha()) then
          call matrixOut%mput(integralIn%alpha,[1,nDimAlpha1],[1,nDimAlpha2]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasAlpha()', integralIn%hasAlpha() )
        endIf
        if (integralIn%hasBeta()) then
          call matrixOut%mput(integralIn%beta,[nDimAlpha1+1,nDimTotal1],[nDimAlpha2+1,nDimTotal2]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasBeta()', integralIn%hasBeta() )
        endIf
      case ('general')
        if(mqc_have_complex(integralIn%alpha)) then
          call matrixOut%init(nDimTotal1,nDimTotal2,cmplx(0.0,0.0))
        else
          call matrixOut%init(nDimTotal1,nDimTotal2)
        endIf
        if (integralIn%hasAlpha()) then
          call matrixOut%mput(integralIn%alpha,[1,nDimAlpha1],[1,nDimAlpha2]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasAlpha()', integralIn%hasAlpha() )
        endIf
        if (integralIn%hasBeta()) then
          call matrixOut%mput(integralIn%beta,[nDimAlpha1+1,nDimTotal1],[nDimAlpha2+1,nDimTotal2]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasBeta()', integralIn%hasBeta() )
        endIf
        if (integralIn%hasAlphaBeta()) then
          call matrixOut%mput(integralIn%alphaBeta,[nDimAlpha1+1,nDimTotal1],[1,nDimAlpha2]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasAlphaBeta()', integralIn%hasAlphaBeta() )
        endIf
        if (integralIn%hasBetaAlpha()) then
          call matrixOut%mput(integralIn%betaAlpha,[1,nDimAlpha1],[nDimAlpha2+1,nDimTotal2]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasBetaAlpha()', integralIn%hasBetaAlpha() )
        endIf
      case default
        call mqc_error_A('unrecognised integer type in mqc_integral_output_array', 6, &
             'integralIn%type()', integralIn%type() )
      end select
!
      end subroutine mqc_integral_output_array  
!
!
!     PROCEDURE MQC_Eigenvalues_Output_Array
!
!>    \brief <b> MQC_Eigenvalues_Output_Array is used to output an MQC eigenvalues
!>    type variable to an MQC vector</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Output_Array is used to output an MQC eigenvalues type variable 
!>    to an MQC matrix.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] VectorOut 
!>    \verbatim
!>        VectorOut is type(mqc_vector)
!>        The MQC vector variable to output.
!>    \endverbatim
!>    
!>    \param[in] EigenvaluesIn
!>    \verbatim
!>        EigenvaluesIn is Class(MQC_SCF_Eigenvalues)
!>        MQC eigenvalues variable to output as MQC vector.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      subroutine mqc_eigenvalues_output_array(vectorOut,eigenvaluesIn) 
!
      implicit none
      class(mqc_scf_eigenvalues),intent(in)::eigenvaluesIn
      type(mqc_vector),intent(inOut)::vectorOut
      integer(kind=int64)::nDimAlpha=0,nDimBeta=0,nDimTotal=0
!
      nDimAlpha = eigenvaluesIn%blockSize('Alpha') 
      nDimBeta = eigenvaluesIn%blockSize('Beta')
      nDimTotal = nDimAlpha + nDimBeta

      select case (eigenvaluesIn%type())
      case('space')
        if(mqc_have_complex(eigenvaluesIn%alpha)) then
          call vectorOut%init(nDimAlpha,cmplx(0.0,0.0))
        else
          call vectorOut%init(nDimAlpha)
        endIf
        if (eigenvaluesIn%hasAlpha()) then
          call vectorOut%vput(eigenvaluesIn%alpha,1) 
        else
          call mqc_error_L('block does not exist in mqc_eigenvalues_output_array', 6, &
               'eigenvaluesIn%hasAlpha()', eigenvaluesIn%hasAlpha() )
        endIf
      case('spin','general')
        if(mqc_have_complex(eigenvaluesIn%alpha)) then
          call vectorOut%init(nDimTotal,cmplx(0.0,0.0))
        else
          call vectorOut%init(nDimTotal)
        endIf
        if (eigenvaluesIn%hasAlpha()) then
          call vectorOut%vput(eigenvaluesIn%alpha,1) 
        else
          call mqc_error_L('block does not exist in mqc_eigenvalues_output_array', 6, &
               'eigenvaluesIn%hasAlpha()', eigenvaluesIn%hasAlpha() )
        endIf
        if (eigenvaluesIn%hasBeta()) then
          call vectorOut%vput(eigenvaluesIn%beta,nDimAlpha+1) 
        else
          call mqc_error_L('block does not exist in mqc_eigenvalues_output_array', 6, &
               'eigenvaluesIn%hasBeta()', eigenvaluesIn%hasBeta() )
        endIf
      case default
        call mqc_error_A('unrecognised integer type in mqc_eigenvalues_output_array', 6, &
             'eigenvaluesIn%type()', eigenvaluesIn%type() )
      end select
!
      end subroutine mqc_eigenvalues_output_array  
!
!
!     PROCEDURE MQC_2ERIs_Output_Array
!
!>    \brief <b> MQC_2ERIs_Output_Array is used to output an MQC 2ERIs type variable
!>    to an MQC rank-4 tensor</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_2ERIs_Output_Array is used to output an MQC 2ERIs type variable to an MQC
!>    rank-4 tensor.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] TensorOut 
!>    \verbatim
!>        TensorOut is type(mqc_r4tensor)
!>        The MQC rank-4 tensor variable to output.
!>    \endverbatim
!>    
!>    \param[in] TwoERIsIn
!>    \verbatim
!>        TwoERIsIn is Class(MQC_TwoERIs)
!>        MQC 2ERIs variable to output as MQC rank-4 tensor.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2019
!
      subroutine mqc_2ERIs_output_array(tensorOut,twoERIsIn) 
!
      implicit none
      type(mqc_twoERIs),intent(in)::twoERIsIn
      type(mqc_r4tensor),intent(inOut)::tensorOut
      integer(kind=int64)::nDim=0
!
      tensorOut = twoERIsIn%getBlock('full')
!
      end subroutine mqc_2ERIs_output_array  
!
!
!     PROCEDURE MQC_Integral_Matrix_Multiply
!
!>    \brief <b> MQC_Integral_Matrix_Multiply is used to multiply an MQC integral
!>    type variable with a MQC matrix and returns an MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Matrix_Multiply is used to multiply an MQC integral type variable
!>    with a MQC matrix and returns an MQC integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        IntegralA is type(mqc_scf_integral)
!>        The MQC integral variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] MatrixB
!>    \verbatim
!>        MaxtrixB is type(mqc_matrix)
!>        MQC Matrix type variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        Description tag of output MQC integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_matrix_multiply(integralA,matrixB,label) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralA
      type(mqc_matrix),intent(in)::matrixB
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      integer(kind=int64)::rows,nBasisAlpha1,nBasisBeta1,nBasisTotal1,nBasisAlpha2,nBasisBeta2, &
        nBasisTotal2
      logical::doOffDiag
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      rows = mqc_matrix_rows(matrixB)
      if(integralA%hasAlpha()) then
        nBasisAlpha1 = integralA%blockSize('Alpha',1)
        nBasisBeta1 = integralA%blockSize('Beta',1)
        nBasisTotal1 = nBasisAlpha1 + nBasisBeta1
        nBasisAlpha2 = integralA%blockSize('Alpha',2)
        nBasisBeta2 = integralA%blockSize('Beta',2)
        nBasisTotal2 = nBasisAlpha2 + nBasisBeta2
        if((rows.eq.nBasisAlpha2).and.(nBasisAlpha2.eq.nBasisBeta2)) then
          doOffDiag = .False.
        elseIf(rows.eq.nBasisTotal2) then
          doOffDiag = .True.
        else
          call mqc_error_I('Integral and matrix are wrongly sized for multiplication', 6, &
               'rows', rows, &
               'nBasisAlpha2', nBasisAlpha2, &
               'nBasisBeta2', nBasisBeta2, &
               'nBasisTotal2', nBasisTotal2 )
        endIf
      else
        call mqc_error_L('Integral is unassigned in MQC_Integral_Matrix_Multiply', 6, &
             'integralA%hasAlpha()', integralA%hasAlpha() )
      endIf
     
      select case (integralA%array_type)
      case('space')
        if(.not.integralA%hasAlpha()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralA%hasAlpha()', integralA%hasAlpha())
        if(doOffDiag) then
          tmpMatrixAlpha = integralA%alpha.dot.matrixB%mat([1,nBasisAlpha2],[1,nBasisAlpha2])
          tmpMatrixBeta = integralA%alpha.dot.matrixB%mat([nBasisAlpha2+1,nBasisTotal2], &
            [nBasisAlpha2+1,nBasisTotal2])
          tmpMatrixAlphaBeta = integralA%alpha.dot.matrixB%mat([nBasisAlpha2+1,nBasisTotal2], &
            [1,nBasisAlpha2])
          tmpMatrixBetaAlpha = integralA%alpha.dot.matrixB%mat([1,nBasisAlpha2], &
            [nBasisAlpha2+1,nBasisTotal2])
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha, &
            tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha) 
        else
          tmpMatrixAlpha = integralA%alpha.dot.matrixB
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha)
        endIf
      case('spin')
        if(.not.integralA%hasAlpha()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralA%hasAlpha()', integralA%hasAlpha() )
        if(.not.integralA%hasBeta()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
'integralA%hasBeta()', integralA%hasBeta() )
        if(doOffDiag) then
          tmpMatrixAlpha = integralA%alpha.dot.matrixB%mat([1,nBasisAlpha2],[1,nBasisAlpha2])
          tmpMatrixBeta = integralA%beta.dot.matrixB%mat([nBasisAlpha2+1,nBasisTotal2], &
            [nBasisAlpha2+1,nBasisTotal2])
          tmpMatrixAlphaBeta = integralA%beta.dot.matrixB%mat([nBasisAlpha2+1,nBasisTotal2], &
            [1,nBasisAlpha2])
          tmpMatrixBetaAlpha = integralA%alpha.dot.matrixB%mat([1,nBasisAlpha2], &
            [nBasisAlpha2+1,nBasisTotal2])
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha, &
            tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha) 
        else
          tmpMatrixAlpha = integralA%alpha.dot.matrixB
          tmpMatrixBeta = integralA%beta.dot.matrixB
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
        endIf
      case('general')
        if(.not.integralA%hasAlpha()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralA%hasAlpha()', integralA%hasAlpha() )
        if(.not.integralA%hasBeta()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralA%hasBeta()', integralA%hasBeta() )
        if(.not.integralA%hasAlphaBeta()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralA%hasAlphaBeta()', integralA%hasAlphaBeta() )
        if(.not.integralA%hasBetaAlpha()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralA%hasBetaAlpha()', integralA%hasBetaAlpha() )
        if(doOffDiag) then
          tmpMatrixAlpha = (integralA%alpha.dot.matrixB%mat([1,nBasisAlpha2],[1,nBasisAlpha2])) + &
            (integralA%betaAlpha.dot.matrixB%mat([nBasisAlpha2+1,nBasisTotal2],[1,nBasisAlpha2]))
          tmpMatrixBeta = (integralA%beta.dot.matrixB%mat([nBasisAlpha2+1,nBasisTotal2],[nBasisAlpha2+1,nBasisTotal2])) + &
            (integralA%alphaBeta.dot.matrixB%mat([1,nBasisAlpha2],[nBasisAlpha2+1,nBasisTotal2]))
          tmpMatrixAlphaBeta = (integralA%beta.dot.matrixB%mat([nBasisAlpha2+1,nBasisTotal2],[1,nBasisAlpha2])) + &
            (integralA%alphaBeta.dot.matrixB%mat([1,nBasisAlpha2],[1,nBasisAlpha2]))
          tmpMatrixBetaAlpha = (integralA%alpha.dot.matrixB%mat([1,nBasisAlpha2],[nBasisAlpha2+1,nBasisTotal2])) + &
            (integralA%betaAlpha.dot.matrixB%mat([nBasisAlpha2+1,nBasisTotal2],[nBasisAlpha2+1,nBasisTotal2])) 
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha, &
            tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha) 
        else
          tmpMatrixAlpha = integralA%alpha.dot.matrixB
          tmpMatrixBeta = integralA%beta.dot.matrixB
          tmpMatrixAlphaBeta = integralA%alphaBeta.dot.matrixB
          tmpMatrixBetaAlpha = integralA%betaAlpha.dot.matrixB
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha, &
            tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        endIf
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_matrix_multiply', 6, &
             'integralA%array_type', integralA%array_type )
      end select
!
      end function mqc_integral_matrix_multiply
!
!
!     PROCEDURE MQC_Matrix_Integral_Multiply
!
!>    \brief <b> MQC_Matrix_Integral_Multiply is used to multiply an MQC matrix type
!>    variable with a MQC integral and returns an MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Matrix_Integral_Multiply is used to multiply an MQC matrix type variable
!>    with a MQC integral and returns an MQC integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] MatrixA
!>    \verbatim
!>        MaxtrixA is type(mqc_matrix)
!>        MQC Matrix type variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] IntegralB
!>    \verbatim
!>        IntegralB is type(mqc_scf_integral)
!>        The MQC integral variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        Description tag of output MQC integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_matrix_integral_multiply(matrixA,integralB,label) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralB
      type(mqc_matrix),intent(in)::matrixA
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      integer(kind=int64)::columns,nBasisAlpha1,nBasisBeta1,nBasisTotal1,nBasisAlpha2,nBasisBeta2, &
        nBasisTotal2
      logical::doOffDiag
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      columns = mqc_matrix_columns(matrixA)
      if(integralB%hasAlpha()) then
        nBasisAlpha1 = integralB%blockSize('Alpha',1)
        nBasisBeta1 = integralB%blockSize('Beta',1)
        nBasisTotal1 = nBasisAlpha1 + nBasisBeta1
        nBasisAlpha2 = integralB%blockSize('Alpha',2)
        nBasisBeta2 = integralB%blockSize('Beta',2)
        nBasisTotal2 = nBasisAlpha2 + nBasisBeta2
        if((columns.eq.nBasisAlpha1).and.(nBasisAlpha1.eq.nBasisBeta1)) then
          doOffDiag = .False.
        elseIf(columns.eq.nBasisTotal1) then
          doOffDiag = .True.
        else
          call mqc_error_I('Integral and matrix are wrongly sized for multiplication', 6, &
               'columns', columns, &
               'nBasisAlpha1', nBasisAlpha1, &
               'nBasisBeta1', nBasisBeta1, &
               'nBasisTotal1', nBasisTotal1 )
        endIf
      else
        call mqc_error_L('Integral is unassigned in MQC_Integral_Matrix_Multiply', 6, &
             'integralB%hasAlpha()', integralB%hasAlpha() )
      endIf
      select case (integralB%array_type)
      case('space')
        if(.not.integralB%hasAlpha()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralB%hasAlpha()', integralB%hasAlpha() )
        if(doOffDiag) then
          tmpMatrixAlpha = matrixA%mat([1,nBasisAlpha1],[1,nBasisAlpha1]).dot.integralB%alpha
          tmpMatrixBeta = matrixA%mat([nBasisAlpha1+1,nBasisTotal1],[nBasisAlpha1+1,nBasisTotal1]) &
            .dot.integralB%alpha
          tmpMatrixAlphaBeta = matrixA%mat([nBasisAlpha1+1,nBasisTotal1],[1,nBasisAlpha1]).dot. &
            integralB%alpha
          tmpMatrixBetaAlpha = matrixA%mat([1,nBasisAlpha1],[nBasisAlpha1+1,nBasisTotal1]).dot. &
            integralB%alpha
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha, &
            tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha) 
        else
          tmpMatrixAlpha = matrixA.dot.integralB%alpha
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha)
        endIf
      case('spin')
        if(.not.integralB%hasAlpha()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralB%hasAlpha()', integralB%hasAlpha() )
        if(.not.integralB%hasBeta()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralB%hasBeta()', integralB%hasBeta() )
        if(doOffDiag) then
          tmpMatrixAlpha = matrixA%mat([1,nBasisAlpha1],[1,nBasisAlpha1]).dot.integralB%alpha
          tmpMatrixBeta = matrixA%mat([nBasisAlpha1+1,nBasisTotal1],[nBasisAlpha1+1,nBasisTotal1]) &
            .dot.integralB%beta
          tmpMatrixAlphaBeta = matrixA%mat([nBasisAlpha1+1,nBasisTotal1],[1,nBasisAlpha1]).dot. &
            integralB%alpha
          tmpMatrixBetaAlpha = matrixA%mat([1,nBasisAlpha1],[nBasisAlpha1+1,nBasisTotal1]).dot. &
            integralB%beta
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha, &
            tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha) 
        else
          tmpMatrixAlpha = matrixA.dot.integralB%alpha
          tmpMatrixBeta = matrixA.dot.integralB%beta 
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
        endIf
      case('general')
        if(.not.integralB%hasAlpha()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralB%hasAlpha()', integralB%hasAlpha() )
        if(.not.integralB%hasBeta()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralB%hasBeta()', integralB%hasBeta() )
        if(.not.integralB%hasAlphaBeta()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralB%hasAlphaBeta()', integralB%hasAlphaBeta() )
        if(.not.integralB%hasBetaAlpha()) call mqc_error_L('Required integral element &
          & unassigned in MQC_Integral_Matrix_Multiply', 6, &
          'integralB%hasBetaAlpha()', integralB%hasBetaAlpha() )
        if(doOffDiag) then
          tmpMatrixAlpha = (matrixA%mat([1,nBasisAlpha1],[1,nBasisAlpha1]).dot.integralB%alpha) + &
            (matrixA%mat([1,nBasisAlpha1],[nBasisAlpha1+1,nBasisTotal1]).dot.integralB%alphaBeta)
          tmpMatrixBeta = (matrixA%mat([nBasisAlpha1+1,nBasisTotal1],[nBasisAlpha1+1,nBasisTotal1]).dot.integralB%beta) + &
            (matrixA%mat([nBasisAlpha1+1,nBasisTotal1],[1,nBasisAlpha1]).dot.integralB%betaAlpha)
          tmpMatrixAlphaBeta = (matrixA%mat([nBasisAlpha1+1,nBasisTotal1],[1,nBasisAlpha1]).dot.&
            integralB%alpha) + (matrixA%mat([nBasisAlpha1+1,nBasisTotal1],[nBasisAlpha1+1,nBasisTotal1]).dot.&
            integralB%alphaBeta)
          tmpMatrixBetaAlpha = (matrixA%mat([1,nBasisAlpha1],[1,nBasisAlpha1]).dot.integralB%betaAlpha) + &
            (matrixA%mat([1,nBasisAlpha1],[nBasisAlpha1+1,nBasisTotal1]).dot.integralB%beta) 
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha, &
            tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha) 
        else
          tmpMatrixAlpha = matrixA.dot.integralB%alpha
          tmpMatrixBeta = matrixA.dot.integralB%beta
          tmpMatrixAlphaBeta = matrixA.dot.integralB%alphaBeta
          tmpMatrixBetaAlpha = matrixA.dot.integralB%betaAlpha
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha, &
            tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        endIf
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_matrix_multiply', 6, &
          'integralB%array_type', integralB%array_type )
      end select
!
      end function mqc_matrix_integral_multiply
!
!
!     PROCEDURE MQC_Integral_Sum
!
!>    \brief <b> MQC_Integral_Sum is used to sum two MQC integral type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Sum is used to sum two MQC integral type variables.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        IntegralA is type(mqc_matrix)
!>        The first MQC integral variable to sum.
!>    \endverbatim
!>    
!>    \param[in] IntegralB
!>    \verbatim
!>        IntegralB is type(mqc_scf_integral)
!>        The second MQC integral variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_integral_sum(integralA,integralB) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralA,integralB
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
!
      select case(integralA%array_type)
      case('space')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          call mqc_integral_allocate(integralOut,'','space',tmpMatrixAlpha)
        case('spin')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          tmpMatrixBeta = integralA%alpha+integralB%beta
          call mqc_integral_allocate(integralOut,'','spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('general')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          tmpMatrixBeta = integralA%alpha+integralB%beta
          tmpMatrixAlphaBeta = integralB%alphaBeta
          tmpMatrixBetaAlpha = integralB%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_sum', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case('spin')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          tmpMatrixBeta = integralA%beta+integralB%alpha
          call mqc_integral_allocate(integralOut,'','spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('spin')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          tmpMatrixBeta = integralA%beta+integralB%beta
          call mqc_integral_allocate(integralOut,'','spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('general')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          tmpMatrixBeta = integralA%beta+integralB%beta
          tmpMatrixAlphaBeta = integralB%alphaBeta
          tmpMatrixBetaAlpha = integralB%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_sum', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case('general')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          tmpMatrixBeta = integralA%beta+integralB%alpha
          tmpMatrixAlphaBeta = integralA%alphaBeta
          tmpMatrixBetaAlpha = integralA%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case('spin')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          tmpMatrixBeta = integralA%beta+integralB%beta
          tmpMatrixAlphaBeta = integralA%alphaBeta
          tmpMatrixBetaAlpha = integralA%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case('general')
          tmpMatrixAlpha = integralA%alpha+integralB%alpha
          tmpMatrixBeta = integralA%beta+integralB%beta
          tmpMatrixAlphaBeta = integralA%alphaBeta+integralB%alphaBeta
          tmpMatrixBetaAlpha = integralA%betaAlpha+integralB%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_sum', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_sum', 6, &
             'integralA%array_type', integralA%array_type )
      end select

      end function mqc_integral_sum
!
!
!     PROCEDURE MQC_Integral_Difference
!
!>    \brief <b> MQC_Integral_Difference is used to subtract two MQC integral type 
!>    variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Difference is used to subtract two MQC integral type variables.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        IntegralA is type(mqc_matrix)
!>        The MQC integral variable to subtract from.
!>    \endverbatim
!>    
!>    \param[in] IntegralB
!>    \verbatim
!>        IntegralB is type(mqc_scf_integral)
!>        The MQC integral variable to subtract.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_integral_difference(integralA,integralB) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralA,integralB
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
!
      select case(integralA%array_type)
      case('space')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          call mqc_integral_allocate(integralOut,'','space',tmpMatrixAlpha)
        case('spin')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          tmpMatrixBeta = integralA%alpha-integralB%beta
          call mqc_integral_allocate(integralOut,'','spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('general')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          tmpMatrixBeta = integralA%alpha-integralB%beta
          tmpMatrixAlphaBeta = integralB%alphaBeta
          tmpMatrixBetaAlpha = integralB%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case('spin')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          tmpMatrixBeta = integralA%beta-integralB%alpha
          call mqc_integral_allocate(integralOut,'','spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('spin')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          tmpMatrixBeta = integralA%beta-integralB%beta
          call mqc_integral_allocate(integralOut,'','spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('general')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          tmpMatrixBeta = integralA%beta-integralB%beta
          tmpMatrixAlphaBeta = integralB%alphaBeta
          tmpMatrixBetaAlpha = integralB%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case('general')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          tmpMatrixBeta = integralA%beta-integralB%alpha
          tmpMatrixAlphaBeta = integralA%alphaBeta
          tmpMatrixBetaAlpha = integralA%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case('spin')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          tmpMatrixBeta = integralA%beta-integralB%beta
          tmpMatrixAlphaBeta = integralA%alphaBeta
          tmpMatrixBetaAlpha = integralA%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case('general')
          tmpMatrixAlpha = integralA%alpha-integralB%alpha
          tmpMatrixBeta = integralA%beta-integralB%beta
          tmpMatrixAlphaBeta = integralA%alphaBeta-integralB%alphaBeta
          tmpMatrixBetaAlpha = integralA%betaAlpha-integralB%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
             'integralA%array_type', integralA%array_type )
      end select

      end function mqc_integral_difference
!
!
!     PROCEDURE MQC_Matrix_Integral_Difference
!
!>    \brief <b> MQC_Matrix_Integral_Difference is used to subtract an MQC integral type 
!>    variable from an MQC matrix</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Matrix_Integral_Difference is used to subtract an MQC integral type variable 
!>    from an MQC variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] MatrixA
!>    \verbatim
!>        IntegralA is type(mqc_matrix)
!>        The MQC integral variable to subtract from.
!>    \endverbatim
!>    
!>    \param[in] IntegralB
!>    \verbatim
!>        IntegralB is type(mqc_scf_integral)
!>        The MQC integral variable to subtract.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2021
!
      function mqc_matrix_integral_difference(matrixA,integralB) result(integralOut)
!
      implicit none
      type(mqc_matrix),intent(in)::matrixA
      type(mqc_scf_integral),intent(in)::integralB
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
!
      select case (integralB%array_type)
      case('space')
        if(size(matrixA,1).eq.integralB%blockSize('alpha',1).and.&
          size(matrixA,2).eq.integralB%blockSize('alpha',2)) then
          tmpMatrixAlpha = matrixA-integralB%alpha
          call mqc_integral_allocate(integralOut,'','space',tmpMatrixAlpha)
        elseIf(size(matrixA,1).eq.2*integralB%blockSize('alpha',1).and.&
          size(matrixA,2).eq.2*integralB%blockSize('alpha',2)) then
          tmpMatrixAlpha = matrixA%mat([1,integralB%blockSize('alpha',1)],[1,integralB%blockSize('alpha',2)])&
            -integralB%alpha
          tmpMatrixBeta = matrixA%mat([integralB%blockSize('alpha',1)+1,-1],[integralB%blockSize('alpha',2)+1,-1])&
            -integralB%alpha
          tmpMatrixAlphaBeta = matrixA%mat([integralB%blockSize('alpha',1)+1,-1],[1,integralB%blockSize('alpha',2)])
          tmpMatrixBetaAlpha = matrixA%mat([1,integralB%blockSize('alpha',1)],[integralB%blockSize('alpha',2)+1,-1])
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        else
          call mqc_error('matrix A and integral B in mqc_matrix_integral_difference are not conformable')
        endIf
      case('spin')
        if(size(matrixA,1).eq.integralB%blockSize('alpha',1).and.&
          size(matrixA,2).eq.integralB%blockSize('alpha',2).and.&
          size(matrixA,1).eq.integralB%blockSize('beta',1).and.&
          size(matrixA,2).eq.integralB%blockSize('beta',2)) then
          tmpMatrixAlpha = matrixA-integralB%alpha
          tmpMatrixBeta = matrixA-integralB%beta
          call mqc_integral_allocate(integralOut,'','spin',tmpMatrixAlpha,tmpMatrixBeta)
        elseIf(size(matrixA,1).eq.integralB%blockSize('alpha',1)+integralB%blockSize('beta',1).and.&
          size(matrixA,2).eq.integralB%blockSize('alpha',2)+integralB%blockSize('beta',2)) then
          tmpMatrixAlpha = matrixA%mat([1,integralB%blockSize('alpha',1)],[1,integralB%blockSize('alpha',2)])&
            -integralB%alpha
          tmpMatrixBeta = matrixA%mat([integralB%blockSize('alpha',1)+1,-1],[integralB%blockSize('alpha',2)+1,-1])&
            -integralB%alpha
          tmpMatrixAlphaBeta = matrixA%mat([integralB%blockSize('alpha',1)+1,-1],[1,integralB%blockSize('alpha',2)])
          tmpMatrixBetaAlpha = matrixA%mat([1,integralB%blockSize('alpha',1)],[integralB%blockSize('alpha',2)+1,-1])
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        else
          call mqc_error('matrix A and integral B in mqc_matrix_integral_difference are not conformable')
        endIf
      case('general')
        if(size(matrixA,1).eq.integralB%blockSize('alpha',1).and.&
          size(matrixA,2).eq.integralB%blockSize('alpha',2).and.&
          size(matrixA,1).eq.integralB%blockSize('beta',1).and.&
          size(matrixA,2).eq.integralB%blockSize('beta',2)) then
          tmpMatrixAlpha = matrixA-integralB%alpha
          tmpMatrixBeta = matrixA-integralB%beta
          tmpMatrixAlphaBeta = (-1)*integralB%alphaBeta
          tmpMatrixBetaAlpha = (-1)*integralB%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        elseIf(size(matrixA,1).eq.integralB%blockSize('alpha',1)+integralB%blockSize('beta',1).and.&
          size(matrixA,2).eq.integralB%blockSize('alpha',2)+integralB%blockSize('beta',2)) then
          tmpMatrixAlpha = matrixA%mat([1,integralB%blockSize('alpha',1)],[1,integralB%blockSize('alpha',2)])&
            -integralB%alpha
          tmpMatrixBeta = matrixA%mat([integralB%blockSize('alpha',1)+1,-1],[integralB%blockSize('alpha',2)+1,-1])&
            -integralB%alpha
          tmpMatrixAlphaBeta = matrixA%mat([integralB%blockSize('alpha',1)+1,-1],[1,integralB%blockSize('alpha',2)])&
            -integralB%alphaBeta
          tmpMatrixBetaAlpha = matrixA%mat([1,integralB%blockSize('alpha',1)],[integralB%blockSize('alpha',2)+1,-1])&
            -integralB%betaAlpha
          call mqc_integral_allocate(integralOut,'','general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        else
          call mqc_error('matrix A and integral B in mqc_matrix_integral_difference are not conformable')
        endIf
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
             'integralB%array_type', integralB%array_type )
      end select

      end function mqc_matrix_integral_difference
!
!
!     PROCEDURE MQC_Integral_Integral_Multiply
!
!>    \brief <b> MQC_Integral_Integral_Multiply is used to multiply two MQC integral
!>    type variables and returns an MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Integral_Multiply is used to multiply two MQC integral type 
!>    variables and returns an MQC integral</b> 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        IntegralA is type(mqc_matrix)
!>        The first MQC integral type variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] IntegralB
!>    \verbatim
!>        IntegralB is type(mqc_scf_integral)
!>        The second MQC integral variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        Description tag of output MQC integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_integral_multiply(integralA,integralB,label) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralA,integralB
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      logical::doOffDiag
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      select case (integralA%array_type)
      case('space')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha)
        case('spin')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%alpha.dot.integralB%beta
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('general')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%alpha.dot.integralB%beta
          tmpMatrixAlphaBeta = integralA%alpha.dot.integralB%alphaBeta
          tmpMatrixBetaAlpha = integralA%alpha.dot.integralB%betaAlpha
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case('spin')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%beta.dot.integralB%alpha
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('spin')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%beta.dot.integralB%beta
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
        case('general')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%beta.dot.integralB%beta
          tmpMatrixAlphaBeta = integralA%beta.dot.integralB%alphaBeta
          tmpMatrixBetaAlpha = integralA%alpha.dot.integralB%betaAlpha
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case('general')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%beta.dot.integralB%alpha
          tmpMatrixAlphaBeta = integralA%alphaBeta.dot.integralB%alpha
          tmpMatrixBetaAlpha = integralA%betaAlpha.dot.integralB%alpha
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case('spin')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%beta.dot.integralB%beta
          tmpMatrixAlphaBeta = integralA%alphaBeta.dot.integralB%alpha
          tmpMatrixBetaAlpha = integralA%betaAlpha.dot.integralB%beta
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case('general')
          tmpMatrixAlpha = (integralA%alpha.dot.integralB%alpha) + (integralA%betaAlpha.dot.integralB%alphaBeta)
          tmpMatrixBeta = (integralA%alphaBeta.dot.integralB%betaAlpha) + (integralA%beta.dot.integralB%beta)
          tmpMatrixAlphaBeta = (integralA%alphaBeta.dot.integralB%alpha) + (integralA%beta.dot.integralB%alphaBeta)
          tmpMatrixBetaAlpha = (integralA%alpha.dot.integralB%betaAlpha) + (integralA%betaAlpha.dot.integralB%beta)
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_integral_multiply', 6, &
             'integralA%array_type', integralA%array_type )
      end select
!
      end function mqc_integral_integral_multiply
!
!
!     PROCEDURE MQC_Scalar_Integral_Multiply
!
!>    \brief <b> MQC_Scalar_Integral_Multiply is used to multiply an MQC scalar with
!>    an MQC integral type variable and returns an MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Scalar_Integral_Multiply is used to multiply an MQC scalar with an MQC 
!>    integral type variable and returns an MQC integral</b> 
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
!>    \param[in] Integral
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The MQC integral variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_scalar_integral_multiply(scalarIn,integral) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integral
      class(*),intent(in)::scalarIn
!      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      type(mqc_scalar)::scalar
      Character(Len=64)::myLabel
!
!      if(present(label)) then
!        call string_change_case(label,'l',myLabel)
!      else
        myLabel = ''
!      endIf

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
        call mqc_error_I('scalarIn type not defined in mqc_scalar_integral_multiply',6)
      end select

      select case (integral%array_type)
      case('space')
          tmpMatrixAlpha = scalar*integral%alpha
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha)
      case('spin')
          tmpMatrixAlpha = scalar*integral%alpha
          tmpMatrixBeta = scalar*integral%beta
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
      case('general')
          tmpMatrixAlpha = scalar*integral%alpha
          tmpMatrixBeta = scalar*integral%beta
          tmpMatrixAlphaBeta = scalar*integral%alphaBeta
          tmpMatrixBetaAlpha = scalar*integral%betaAlpha
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
      case default
        call mqc_error_A('Unknown integral type in mqc_scalar_integral_multiply', 6, &
             'integralA%array_type', integral%array_type )
      end select
!
      end function mqc_scalar_integral_multiply
!
!
!     PROCEDURE MQC_Integral_Scalar_Multiply
!
!>    \brief <b> MQC_Integral_Scalar_Multiply is used to multiply an MQC integral 
!>    type variable with an MQC scalar and returns an MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Scalar_Multiply is used to multiply an MQC integral type variable 
!>    with an MQC scalar and returns an MQC integral. 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The MQC integral variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] ScalarIn
!>    \verbatim
!>        ScalarIn is class(*)
!>        The MQC scalar type variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_integral_scalar_multiply(integral,scalarIn) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integral
      class(*),intent(in)::scalarIn
!      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      type(mqc_scalar)::scalar
      Character(Len=64)::myLabel
!
!      if(present(label)) then
!        call string_change_case(label,'l',myLabel)
!      else
        myLabel = ''
!      endIf

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
        call mqc_error_I('scalarIn type not defined in mqc_integral_scalar_multiply',6)
      end select

      select case (integral%array_type)
      case('space')
          tmpMatrixAlpha = scalar*integral%alpha
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha)
      case('spin')
          tmpMatrixAlpha = scalar*integral%alpha
          tmpMatrixBeta = scalar*integral%beta
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
      case('general')
          tmpMatrixAlpha = scalar*integral%alpha
          tmpMatrixBeta = scalar*integral%beta
          tmpMatrixAlphaBeta = scalar*integral%alphaBeta
          tmpMatrixBetaAlpha = scalar*integral%betaAlpha
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_scalar_multiply', 6, &
             'integralA%array_type', integral%array_type )
      end select
!
      end function mqc_integral_scalar_multiply
!
!
!     PROCEDURE MQC_Integral_Eigenvalues_Multiply
!
!>    \brief <b> MQC_Integral_Eigenvalues_Multiply is used to multiply an MQC 
!>    integral type variable with an MQC eigenvalues type variable (a diagonal 
!>    matrix) and returns an MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Eigenvalues_Multiply is used to multiply an MQC integral type 
!>    variable with an MQC eigenvalues type variable (a diagonal matrix) and returns 
!>    an MQC integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        IntegralA is type(mqc_scf_integral)
!>        The MQC integral variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] EigenvaluesB
!>    \verbatim
!>        EigenvaluesB is type(mqc_scf_eigenvalues)
!>        The MQC eigenvalues type variable (a diagonal matrix) 
!>        to multiply.
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        Description tag of output MQC integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_integral_eigenvalues_multiply(integralA,eigenvaluesB,label) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralA
      type(mqc_scf_eigenvalues),intent(in)::eigenvaluesB
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrix1,tmpMatrix2,tmpMatrix3,tmpMatrix4
      type(mqc_vector)::tmpVector1,tmpVector2
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      select case (integralA%array_type)
      case('space')
        select case (eigenvaluesB%array_type)
        case('space')
          tmpVector1 = eigenvaluesB%alpha
          call tmpVector1%diag(tmpMatrix1)
          tmpMatrix1 = matmul(integralA%alpha,tmpMatrix1)
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrix1)
        case('spin','general')
          tmpVector1 = eigenvaluesB%alpha
          tmpVector2 = eigenvaluesB%beta
          call tmpVector1%diag(tmpMatrix1)
          call tmpVector2%diag(tmpMatrix2)
          tmpMatrix1 = matmul(integralA%alpha,tmpMatrix1)
          tmpMatrix2 = matmul(integralA%alpha,tmpMatrix2)
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrix1,tmpMatrix2)
        case default
          call mqc_error_A('Unknown eigenvalues type in mqc_integral_eigenvalues_multiply', 6, &
               'eigenvaluesB%array_type', eigenvaluesB%array_type )
        end select
      case('spin')
        select case (eigenvaluesB%array_type)
        case('space')
          tmpVector1 = eigenvaluesB%alpha
          call tmpVector1%diag(tmpMatrix1)
          tmpMatrix2 = matmul(integralA%beta,tmpMatrix1)
          tmpMatrix1 = matmul(integralA%alpha,tmpMatrix1)
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrix1,tmpMatrix2)
        case('spin','general')
          tmpVector1 = eigenvaluesB%alpha
          tmpVector2 = eigenvaluesB%beta
          call tmpVector1%diag(tmpMatrix1)
          call tmpVector2%diag(tmpMatrix2)
          tmpMatrix1 = matmul(integralA%alpha,tmpMatrix1)
          tmpMatrix2 = matmul(integralA%beta,tmpMatrix2)
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrix1,tmpMatrix2)
        case default
          call mqc_error_A('Unknown eigenvalues type in mqc_integral_eigenvalues_multiply', 6, &
               'eigenvaluesB%array_type', eigenvaluesB%array_type )
        end select
      case('general')
        select case (eigenvaluesB%array_type)
        case('space')
          tmpVector1 = eigenvaluesB%alpha
          call tmpVector1%diag(tmpMatrix1)
          tmpMatrix2 = matmul(integralA%beta,tmpMatrix1)
          tmpMatrix3 = matmul(integralA%alphaBeta,tmpMatrix1)
          tmpMatrix4 = matmul(integralA%betaAlpha,tmpMatrix1)
          tmpMatrix1 = matmul(integralA%alpha,tmpMatrix1)
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrix1,tmpMatrix2,tmpMatrix3,tmpMatrix4)
        case('spin','general')
          tmpVector1 = eigenvaluesB%alpha
          tmpVector2 = eigenvaluesB%beta
          call tmpVector1%diag(tmpMatrix1)
          call tmpVector2%diag(tmpMatrix2)
          tmpMatrix3 = matmul(integralA%alphaBeta,tmpMatrix1)
          tmpMatrix4 = matmul(integralA%betaAlpha,tmpMatrix2)
          tmpMatrix1 = matmul(integralA%alpha,tmpMatrix1)
          tmpMatrix2 = matmul(integralA%beta,tmpMatrix2)
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrix1,tmpMatrix2,tmpMatrix3,tmpMatrix4)
        case default
          call mqc_error_A('Unknown eigenvalues type in mqc_integral_eigenvalues_multiply', 6, &
               'eigenvaluesB%array_type', eigenvaluesB%array_type )
        end select
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_eigenvalues_multiply', 6, &
             'integralA%array_type', integralA%array_type )
      end select
!
      end function mqc_integral_eigenvalues_multiply
!
!
!     PROCEDURE MQC_Eigenvalues_Integral_Multiply
!
!>    \brief <b> MQC_Integral_Eigenvalues_Multiply is used to multiply an MQC 
!>    integral type variable with an MQC eigenvalues type variable (a diagonal 
!>    matrix) and returns an MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Eigenvalues_Multiply is used to multiply an MQC integral type 
!>    variable with an MQC eigenvalues type variable (a diagonal matrix) and returns 
!>    an MQC integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] EigenvaluesA
!>    \verbatim
!>        EigenvaluesA is type(mqc_scf_eigenvalues)
!>        The MQC eigenvalues type variable (a diagonal matrix) 
!>        to multiply.
!>    \endverbatim
!>
!>    \param[in] IntegralB
!>    \verbatim
!>        IntegralB is type(mqc_scf_integral)
!>        The MQC integral variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        Description tag of output MQC integral variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_eigenvalues_integral_multiply(eigenvaluesA,integralB,label) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralB
      type(mqc_scf_eigenvalues),intent(in)::eigenvaluesA 
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMatrix1,tmpMatrix2,tmpMatrix3,tmpMatrix4
      type(mqc_vector)::tmpVector1,tmpVector2
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      select case (integralB%array_type)
      case('space')
        select case (eigenvaluesA%array_type)
        case('space')
          tmpVector1 = eigenvaluesA%alpha
          call tmpVector1%diag(tmpMatrix1)
          tmpMatrix1 = matmul(tmpMatrix1,integralB%alpha)
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrix1)
        case('spin','general')
          tmpVector1 = eigenvaluesA%alpha
          tmpVector2 = eigenvaluesA%beta
          call tmpVector1%diag(tmpMatrix1)
          call tmpVector2%diag(tmpMatrix2)
          tmpMatrix1 = matmul(tmpMatrix1,integralB%alpha)
          tmpMatrix2 = matmul(tmpMatrix2,integralB%alpha)
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrix1,tmpMatrix2)
        case default
          call mqc_error_A('Unknown eigenvalues type in mqc_integral_eigenvalues_multiply', 6, &
               'eigenvaluesA%array_type', eigenvaluesA%array_type )
        end select
      case('spin')
        select case (eigenvaluesA%array_type)
        case('space')
          tmpVector1 = eigenvaluesA%alpha
          call tmpVector1%diag(tmpMatrix1)
          tmpMatrix2 = matmul(tmpMatrix1,integralB%beta)
          tmpMatrix1 = matmul(tmpMatrix1,integralB%alpha)
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrix1,tmpMatrix2)
        case('spin','general')
          tmpVector1 = eigenvaluesA%alpha
          tmpVector2 = eigenvaluesA%beta
          call tmpVector1%diag(tmpMatrix1)
          call tmpVector2%diag(tmpMatrix2)
          tmpMatrix1 = matmul(tmpMatrix1,integralB%alpha)
          tmpMatrix2 = matmul(tmpMatrix2,integralB%beta)
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrix1,tmpMatrix2)
        case default
          call mqc_error_A('Unknown eigenvalues type in mqc_integral_eigenvalues_multiply', 6, &
               'eigenvaluesA%array_type', eigenvaluesA%array_type )
        end select
      case('general')
        select case (eigenvaluesA%array_type)
        case('space')
          tmpVector1 = eigenvaluesA%alpha
          call tmpVector1%diag(tmpMatrix1)
          tmpMatrix2 = matmul(tmpMatrix1,integralB%beta)
          tmpMatrix3 = matmul(tmpMatrix1,integralB%alphaBeta)
          tmpMatrix4 = matmul(tmpMatrix1,integralB%betaAlpha)
          tmpMatrix1 = matmul(tmpMatrix1,integralB%alpha)
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrix1,tmpMatrix2,tmpMatrix3,tmpMatrix4)
        case('spin','general')
          tmpVector1 = eigenvaluesA%alpha
          tmpVector2 = eigenvaluesA%beta
          call tmpVector1%diag(tmpMatrix1)
          call tmpVector2%diag(tmpMatrix2)
          tmpMatrix3 = matmul(tmpMatrix2,integralB%alphaBeta)
          tmpMatrix4 = matmul(tmpMatrix1,integralB%betaAlpha)
          tmpMatrix1 = matmul(tmpMatrix1,integralB%alpha)
          tmpMatrix2 = matmul(tmpMatrix2,integralB%beta)
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrix1,tmpMatrix2,tmpMatrix3,tmpMatrix4)
        case default
          call mqc_error_A('Unknown eigenvalues type in mqc_integral_eigenvalues_multiply', 6, &
               'eigenvaluesA%array_type', eigenvaluesA%array_type )
        end select
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_eigenvalues_multiply', 6, &
             'integralB%array_type', integralB%array_type )
      end select
!
      end function mqc_eigenvalues_integral_multiply
!
!
!     PROCEDURE MQC_Eigenvalues_Eigenvalues_Multiply
!
!>    \brief <b> MQC_Eigenvalues_Eigenvalues_Multiply is used to multiply two MQC 
!>    eigenvalues type variable (diagonal matrices) and returns an MQC eigenvalues
!>    type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalues_Eigenvalues_Multiply is used to multiply two MQC eigenvalues 
!>    type variable (diagonal matrices) and returns an MQC eigenvalues type variables.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] EigenvaluesA
!>    \verbatim
!>        EigenvaluesA is type(mqc_scf_eigenvalues)
!>        The first MQC eigenvalues type variable (a diagonal 
!>        matrix) to multiply.
!>    \endverbatim
!>
!>    \param[in] EigenvaluesB
!>    \verbatim
!>        EigenvaluesB is type(mqc_scf_eigenvalues)
!>        The second MQC eigenvalues type variable (a diagonal 
!>        matrix) to multiply.
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        Description tag of output MQC eigenvalues variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_eigenvalues_eigenvalues_multiply(eigenvaluesA,eigenvaluesB,label) result(eigenvaluesOut)
!
      implicit none
      type(mqc_scf_eigenvalues),intent(in)::eigenvaluesA,eigenvaluesB 
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_eigenvalues)::eigenvaluesOut
      type(mqc_vector)::tmpVector1,tmpVector2,tmpVector3,tmpVector4
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      select case (eigenvaluesB%array_type)
      case('space')
        select case (eigenvaluesA%array_type)
        case('space')
          tmpVector1 = eigenvaluesA%alpha
          tmpVector2 = eigenvaluesB%alpha
          tmpVector1 = tmpVector1.ewp.tmpVector2
          call mqc_eigenvalues_allocate(eigenvaluesOut,myLabel,'space',tmpVector1)
        case('spin','general')
          tmpVector1 = eigenvaluesA%alpha
          tmpVector2 = eigenvaluesA%beta
          tmpVector3 = eigenvaluesB%alpha
          tmpVector1 = tmpVector1.ewp.tmpVector3
          tmpVector2 = tmpVector2.ewp.tmpVector3
          call mqc_eigenvalues_allocate(eigenvaluesOut,myLabel,'spin',tmpVector1,tmpVector2)
        case default
          call mqc_error_A('Unknown eigenvalues type in mqc_eigenvalues_eigenvalues_multiply', 6, &
               'eigenvaluesA%array_type', eigenvaluesA%array_type )
        end select
      case('spin','general')
        select case (eigenvaluesA%array_type)
        case('space')
          tmpVector1 = eigenvaluesA%alpha
          tmpVector2 = eigenvaluesB%alpha
          tmpVector3 = eigenvaluesB%beta
          tmpVector2 = tmpVector1.ewp.tmpVector2
          tmpVector3 = tmpVector1.ewp.tmpVector3
          call mqc_eigenvalues_allocate(eigenvaluesOut,myLabel,'spin',tmpVector2,tmpVector3)
        case('spin','general')
          tmpVector1 = eigenvaluesA%alpha
          tmpVector2 = eigenvaluesA%beta
          tmpVector3 = eigenvaluesB%alpha
          tmpVector4 = eigenvaluesB%beta
          tmpVector1 = tmpVector1.ewp.tmpVector3
          tmpVector2 = tmpVector2.ewp.tmpVector4
          call mqc_eigenvalues_allocate(eigenvaluesOut,myLabel,'spin',tmpVector1,tmpVector2)
        case default
          call mqc_error_A('Unknown eigenvalues type in mqc_eigenvalues_eigenvalues_multiply', 6, &
               'eigenvaluesA%array_type', eigenvaluesA%array_type )
        end select
      case default
        call mqc_error_A('Unknown eigenvalues type in mqc_eigenvalues_eigenvalues_multiply', 6, &
             'eigenvaluesB%array_type', eigenvaluesB%array_type )
      end select
!
      end function mqc_eigenvalues_eigenvalues_multiply
!
!
!     PROCEDURE MQC_Eigenvalue_Eigenvalue_DotProduct
!
!>    \brief <b> MQC_Eigenvalue_Eigenvalue_DotProduct is used to compute the inner
!>    product of two MQC eigenvalues type variables and returns an MQC eigenvalues
!>    type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Eigenvalue_Eigenvalue_DotProduct is used to compute the inner product of 
!>    two MQC eigenvalues type variables and returns an MQC eigenvalues type variable. 
!>    This routine is consistent whether considering MQC eigenvalues as a vector or 
!>    diagonal matrix.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] EigenvalueA
!>    \verbatim
!>        EigenvaluesA is type(mqc_scf_eigenvalues)
!>        The first MQC eigenvalues type variable.
!>    \endverbatim
!>
!>    \param[in] EigenvalueB
!>    \verbatim
!>        EigenvaluesB is type(mqc_scf_eigenvalues)
!>        The second MQC eigenvalues type variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_eigenvalue_eigenvalue_dotProduct(eigenvalueA,eigenvalueB) result(scalarOut)
!
      implicit none
      type(mqc_scf_eigenvalues),intent(in)::eigenvalueA,eigenvalueB
      type(mqc_scalar)::scalarOut,two

      two = 2.0 

      select case (eigenvalueA%array_type)
      case('space')
        select case (eigenvalueB%array_type)
        case('space')
          scalarOut = two*(dot_product(eigenvalueA%alpha,eigenvalueB%alpha)) 
        case('spin')
          scalarOut = dot_product(eigenvalueA%alpha,eigenvalueB%alpha) + &
            dot_product(eigenvalueA%alpha,eigenvalueB%beta)
        case('general')
          scalarOut = dot_product(eigenvalueA%alpha,eigenvalueB%alpha) + &
            dot_product(eigenvalueA%alpha,eigenvalueB%beta)
        case default
          call mqc_error_A('Unknown integral type in mqc_eigenvalue_eigenvalue_dotProduct', 6, &
               'eigenvalueB%array_type', eigenvalueB%array_type )
        end select
      case('spin')
        select case (eigenvalueB%array_type)
        case('space')
          scalarOut = dot_product(eigenvalueA%alpha,eigenvalueB%alpha) + &
            dot_product(eigenvalueA%beta,eigenvalueB%alpha)
        case('spin')
          scalarOut = dot_product(eigenvalueA%alpha,eigenvalueB%alpha) + &
            dot_product(eigenvalueA%beta,eigenvalueB%beta)
        case('general')
          scalarOut = dot_product(eigenvalueA%alpha,eigenvalueB%alpha) + &
            dot_product(eigenvalueA%beta,eigenvalueB%beta)
        case default
          call mqc_error_A('Unknown integral type in mqc_eigenvalue_eigenvalue_dotProduct', 6, &
               'eigenvalueB%array_type', eigenvalueB%array_type )
        end select
      case('general')
        select case (eigenvalueB%array_type)
        case('space')
          scalarOut = dot_product(eigenvalueA%alpha,eigenvalueB%alpha) + &
            dot_product(eigenvalueA%beta,eigenvalueB%alpha)
        case('spin')
          scalarOut = dot_product(eigenvalueA%alpha,eigenvalueB%alpha) + &
            dot_product(eigenvalueA%beta,eigenvalueB%beta)
        case('general')
          scalarOut = dot_product(eigenvalueA%alpha,eigenvalueB%alpha) + &
            dot_product(eigenvalueA%beta,eigenvalueB%beta)
        case default
          call mqc_error_A('Unknown integral type in mqc_eigenvalue_eigenvalue_dotProduct', 6, &
               'eigenvalueB%array_type', eigenvalueB%array_type )
        end select
      case default
        call mqc_error_A('Unknown eigenvalue type in mqc_eigenvalue_eigenvalue_dotProduct', 6, &
             'eigenvalueA%array_type', eigenvalueA%array_type )
      end select
!
      end function mqc_eigenvalue_eigenvalue_dotProduct
!
!
!     PROCEDURE MQC_Integral_Transpose
!
!>    \brief <b> MQC_Integral_Transpose is used to compute the transpose of a MQC
!>    integral type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Transpose is used to compute the transpose of a MQC integral type 
!>    variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The MQC integral type variable to transpose.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        A description tag for the transposed integral.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_transpose(integral,label) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integral
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      logical::doOffDiag
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf
     
      select case (integral%array_type)
      case('space')
          tmpMatrixAlpha = transpose(integral%alpha)
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha)
      case('spin')
          tmpMatrixAlpha = transpose(integral%alpha)
          tmpMatrixBeta = transpose(integral%beta)
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
      case('general')
          tmpMatrixAlpha = transpose(integral%alpha)
          tmpMatrixBeta = transpose(integral%beta)
          tmpMatrixAlphaBeta = transpose(integral%betaAlpha)
          tmpMatrixBetaAlpha = transpose(integral%alphaBeta)
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_transpose', 6, &
             'integral%array_type', integral%array_type )
      end select
!
      end function mqc_integral_transpose
!
!
!     PROCEDURE MQC_Integral_Conjugate_Transpose
!
!>    \brief <b> MQC_Integral_Conjugate_Transpose is used to compute the conjugate 
!>    transpose of a MQC integral type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Conjugate_Transpose is used to compute the conjugate transpose of 
!>    a MQC integral type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The MQC integral type variable to conjugate transpose.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        A description tag for the conjugate transposed integral.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017
!
      function mqc_integral_conjugate_transpose(integral,label) result(integralOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integral
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      logical::doOffDiag
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      select case (integral%array_type)
      case('space')
          tmpMatrixAlpha = dagger(integral%alpha)
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha)
      case('spin')
          tmpMatrixAlpha = dagger(integral%alpha)
          tmpMatrixBeta = dagger(integral%beta)
          call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMatrixAlpha,tmpMatrixBeta)
      case('general')
          tmpMatrixAlpha = dagger(integral%alpha)
          tmpMatrixBeta = dagger(integral%beta)
          tmpMatrixAlphaBeta = dagger(integral%betaAlpha)
          tmpMatrixBetaAlpha = dagger(integral%alphaBeta)
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_conjugate_transpose', 6, &
             'integral%array_type', integral%array_type )
      end select
!
      end function mqc_integral_conjugate_transpose
!
!
!     PROCEDURE MQC_Integral_Norm
!
!>    \brief <b> MQC_Integral_Norm is used to compute the norm of a MQC integral type 
!>    variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Norm is used to compute the norm of a MQC integral type variable.
!>    The following options are available:
!>
!>    1. methodIn = 'M' uses the maximum absolute value max(abs(A(i,j))).
!>    2. methodIn = '1' uses the one norm.
!>    3. methodIn = 'I' uses the infinity norm.
!>    4. methodIn = 'F' uses the Frobenius norm (default).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The MQC integral type variable to determine the norm.
!>    \endverbatim
!>
!>    \param[in] MethodIn
!>    \verbatim
!>        MethodIn is Character(len=1)
!>        = 'M': max(abs(A(i,j)))
!>        = '1': one norm
!>        = 'I': infinity norm
!>        = 'F': Frobenius norm.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_integral_norm(integral,methodIn) result(norm)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      character(len=1),optional,intent(in)::methodIn
      character(len=1)::method
      type(mqc_scalar)::norm
      type(mqc_matrix)::ret_block
!
      if(Present(methodIn)) then
        method = methodIn
      else
        method = 'F'
      endIf
!      
       ret_block = integral%getBlock('full')
       norm = mqc_matrix_norm(ret_block,method)
!
      end function mqc_integral_norm
!
!
!
!     PROCEDURE MQC_Matrix_SpinBlockGHF
!
!>    \brief <b> MQC_Matrix_SpinBlockGHF is used to return a spin blocked MQC integral 
!>    type variable from a spin-interleaved array</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>     This subroutine takes a GHF spin-interleaved spin array (such as typically 
!>     generated by e.g. Gaussian) and returns it in spin blocked form i.e. the spinor 
!>     basis is reordered so that alpha coefficients are before beta coefficients.
!>
!>     If MOs are passed, the input matrix columns ae reordered to 'spin block' the 
!>     orbitals. The input MOs are typically ordered by energy. Spin blocking is 
!>     performed by determining alpha or beta dominance of each occupied orbital 
!>     vector through comparing the sum of squares of orbital coefficients. Orbitals 
!>     are places in their corresponding spin block from left to right in order of 
!>     spin dominance. If more alpha or beta dominant orbitals are identified than 
!>     possible through spin multiplicity, the least dominant orbitals are placed in 
!>     the unfilled spin block with order determined by the contribution of spin of 
!>     the block into which the orbitals are placed. The procedure is then repeated 
!>     for virtual orbitals.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Array
!>    \verbatim
!>        Array is class(*)
!>        The Array (MQC matrix or MQC vector) to be spin blocked.
!>    \endverbatim
!>
!>    \param[in] NElec
!>    \verbatim
!>        NElec is integer(kind=int64),optional
!>        The number of electrons. Only required when spin 
!>        blocking MOs.
!>    \endverbatim
!>     
!>    \param[in] Multi
!>    \verbatim
!>        Multi is Character(len=1)
!>        The multiplicity. Only required when spin blocking MOs.
!>    \endverbatim
!>     
!>    \param[in,out] EList
!>    \verbatim
!>        Elist is integer(kind=int64),dimension(:),allocatable,optional
!>        The energy list returns the energy order of MOs required
!>        for spin unblocking inverse procedure. When allocating 
!>        an MQC integral, EList should be used to fill the Elist 
!>        object. Only required when spin blocking MOs.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \author A. D. Mahler
!>    \date 2017, 2018
!
      subroutine mqc_matrix_spinBlockGHF(array,nelec,multi,elist)
!
      implicit none
      integer(kind=int64)::i,j,k,rows,columns,blockSize,acols,bcols,nalphas,nbetas, &
        athresh,bthresh
      class(*),intent(inOut)::array
      type(mqc_matrix)::tmpMatrix
      type(mqc_vector)::tmpVector
      type(mqc_scalar)::aBlock,bBlock
      character(len=1),dimension(:),allocatable::abarray
      integer(kind=int64),optional::nelec,multi
      integer(kind=int64),dimension(:),allocatable,optional::elist
!
!     Check that both optional arguments are defined or both
!     are absent
      if(present(nelec).neqv.present(multi)) then
        call mqc_error_I('NElectrons, Multiplicity, or energy list missing &
     &    in mqc_matrix_spinBlockGHF', 6)
      end if
!
!
!     Do the work...
!
      select type (array)
      type is (mqc_vector)
        rows = mqc_length_vector(array)
        call tmpVector%init(rows)
        j = 1
        k = rows/2+1
        do i = 1,rows
          if(mod(i,2).eq.1) then
            call tmpVector%put(array%at(i),j)
            j = j+1
          elseIf(mod(i,2).eq.0) then
            call tmpVector%put(array%at(i),k)
            k = k+1
          else
            call mqc_error_I('mqc_matrix_spinBlockGHF is confused.', 6, &
                 'mod(i,2)', mod(i,2) )
          endIf
        endDo
        array = tmpVector
      type is (mqc_matrix)
        if(present(nelec)) then
          rows = mqc_matrix_rows(array)
          columns = mqc_matrix_columns(array)
          blockSize = rows/2
          allocate(abarray(columns))
          if(present(elist)) allocate(elist(columns))
          call tmpMatrix%init(rows,columns)
          j = 1
          k = rows/2+1
          do i = 1,rows
            if(mod(i,2).eq.1) then
              call tmpMatrix%vput(array%vat([i],[0]),[j],[0])
              j = j+1
            elseIf(mod(i,2).eq.0) then
              call tmpMatrix%vput(array%vat([i],[0]),[k],[0])
              k = k+1
            else
              call mqc_error_i('mqc_matrix_spinBlockGHF is confused.', 6, &
                  'mod(i,2)', mod(i,2) )
            endIf
          endDo
!
!     Begin new column ordering procedure for columns
!     and intialize spin ordering array
!
!
!     Determine if each column is alpha- or beta-dominant
!
          DO I = 1, COLUMNS
            aBlock = 0.0
            bBlock = 0.0
            DO J = 1, blockSize
              aBlock = aBlock + (tmpMatrix%at(J,I) * tmpMatrix%at(J,I))
              bBlock = bBlock + (tmpMatrix%at(blockSize+J,I) * tmpMatrix%at(blockSize+J,I))
            END DO
            IF (aBlock%abs().ge.bBlock%abs()) THEN
              abarray(I) = 'a'
            ELSE
              abarray(I) = 'b'
            ENDIF
          END DO
!
!     Delace columns based on alpha/beta dominance
!
          acols = 1
          bcols = blocksize + 1
          athresh = (nelec+(multi-1))/2
          bthresh = nelec - athresh
          nalphas = 0
          nbetas = 0
!
!     Occupied orbitals first
!
          DO I = 1, NELEC
            !Check if a or b bucket is filled
            IF ((nalphas.lt.athresh).and.(nbetas.lt.bthresh)) THEN
              IF(abarray(I).eq.'a') THEN
                call array%vput(tmpMatrix%vat([0],[I]),[0],[acols])
                if(present(elist)) elist(acols) = I
                acols = acols + 1
                nalphas = nalphas + 1
              ELSE
                call array%vput(tmpMatrix%vat([0],[I]),[0],[bcols])
                if(present(elist)) elist(bcols) = I
                bcols = bcols + 1
                nbetas = nbetas + 1
              ENDIF
            !If a or b bucket is filled, check b bucket
            ELSE IF (nbetas.ge.bthresh) THEN
              call array%vput(tmpMatrix%vat([0],[I]),[0],[acols])
              if(present(elist)) elist(acols) = I
              acols = acols + 1
              nalphas = nalphas + 1
            !Otherwise, a bucket is unfilled
            ELSE
              call array%vput(tmpMatrix%vat([0],[I]),[0],[bcols])
              if(present(elist)) elist(bcols) = I
              bcols = bcols + 1
              nbetas = nbetas + 1
            ENDIF
          END DO
!
!     Virtual orbitals
!
          athresh = blocksize - athresh
          bthresh = blocksize - bthresh
          nalphas = 0
          nbetas = 0
          DO I = NELEC+1, COLUMNS
            !Check if a or b bucket is filled
            IF((nalphas.lt.athresh).and.(nbetas.lt.bthresh)) THEN
              IF(abarray(I).eq.'a') THEN
                call array%vput(tmpMatrix%vat([0],[I]),[0],[acols])
                if(present(elist)) elist(acols) = I
                acols = acols + 1
                nalphas = nalphas + 1
              ELSE
                call array%vput(tmpMatrix%vat([0],[I]),[0],[bcols])
                if(present(elist)) elist(bcols) = I
                bcols = bcols + 1
                nbetas = nbetas + 1
              ENDIF
            !If a or b bucket is filled, check b bucket
            ELSE IF(nbetas.ge.bthresh) THEN
              call array%vput(tmpMatrix%vat([0],[I]),[0],[acols])
              if(present(elist)) elist(acols) = I
              acols = acols + 1
              nalphas = nalphas + 1
            ELSE
              call array%vput(tmpMatrix%vat([0],[I]),[0],[bcols])
              if(present(elist)) elist(bcols) = I
              bcols = bcols + 1
              nbetas = nbetas + 1
            ENDIF
          END DO

          if(mqc_matrix_haveFull(array).and.mqc_matrix_test_symmetric(array)) call mqc_matrix_full2Symm(array)
        else
          rows = mqc_matrix_rows(array)
          columns = mqc_matrix_columns(array)
          call tmpMatrix%init(rows,columns)
          j = 1
          k = rows/2+1
          do i = 1,rows
            if(mod(i,2).eq.1) then
              call tmpMatrix%vput(array%vat([i],[0]),[j],[0])
              j = j+1
            elseIf(mod(i,2).eq.0) then
              call tmpMatrix%vput(array%vat([i],[0]),[k],[0])
              k = k+1
            else
              call mqc_error_i('mqc_matrix_spinBlockGHF is confused.', 6, &
                  'mod(i,2)', mod(i,2) )
            endIf
          endDo
          j = 1
          k = columns/2+1
          do i = 1,columns
            if(mod(i,2).eq.1) then
              call array%vput(tmpMatrix%vat([0],[i]),[0],[j])
              j = j+1
            elseIf(mod(i,2).eq.0) then
              call array%vput(tmpMatrix%vat([0],[i]),[0],[k])
              k = k+1
            else
              call mqc_error_I('mqc_matrix_spinBlockGHF is confused.', 6, &
                  'mod(i,2)', mod(i,2) )
            endIf
          endDo
          if(mqc_matrix_haveFull(array).and.mqc_matrix_test_symmetric(array)) call mqc_matrix_full2Symm(array)
        end if
      class default
        call mqc_error_I('unrecognised array type in mqc_matrix_spinBlockGHF', 6)
      end select 
!
      return
      end subroutine mqc_matrix_spinBlockGHF
!
!
!     PROCEDURE MQC_Matrix_UndoSpinBlockGHF_Eigenvalues
!
!>    \brief <b> MQC_Matrix_UndoSpinBlockGHF_Eigenvalues is used to return a spin 
!>    unblocked MQC vector from a spin blocked MQC eigenvalues variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Matrix_UndoSpinBlockGHF_Eigenvalues is used to return a spin unblocked MQC 
!>    vector from a spin blocked MQC eigenvalues variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] EigenvaluesIn
!>    \verbatim
!>        EigenvaluesIn is type(mqc_scf_eigenvalues) 
!>        The MQC eigenvalues variable to be spin unblocked.
!>    \endverbatim
!>
!>    \param[out] VectorOut
!>    \verbatim
!>        VectorOut is type(mqc_vector)
!>        The MQC vector output.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \author A. D. Mahler
!>    \date 2017, 2018
!
      subroutine mqc_matrix_undoSpinBlockGHF_EigenValues(eigenvaluesIn,vectorOut)
!
      implicit none
      integer(kind=int64)::i,j,nAlpha,nBeta,nElectrons
      type(mqc_scf_eigenvalues),intent(in)::eigenvaluesIn
      type(mqc_vector),intent(out)::vectorOut
      type(mqc_vector)::tmpVectorAlpha,tmpVectorBeta
!
!
!     Do the work...
!
      nAlpha = eigenvaluesIn%blockSize('alpha') 
      nBeta = eigenvaluesIn%blockSize('beta') 
      nElectrons = nAlpha+nBeta
      call vectorOut%init(nElectrons)
      tmpVectorAlpha = eigenvaluesIn%getBlock('alpha')
      tmpVectorBeta = eigenvaluesIn%getBlock('beta')
      j = 1
      do i = 1,nAlpha
        call vectorOut%put(tmpVectorAlpha%at(i),j)
        j = j+2
      endDo
      j = 2
      do i = 1,nBeta 
        call vectorOut%put(tmpVectorBeta%at(i),j)
        j = j+2
      endDo
!
      return
      end subroutine mqc_matrix_undoSpinBlockGHF_EigenValues
!
!
!     PROCEDURE MQC_Matrix_UndoSpinBlockGHF_Integral   
!
!>    \brief <b> MQC_Matrix_UndoSpinBlockGHF_Integral is used to return a spin 
!>    unblocked MQC matrix from a spin blocked MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>     MQC_Matrix_UndoSpinBlockGHF_Integral is used to return a spin unblocked MQC 
!>     matrix from a spin blocked MQC integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IntegralIn
!>    \verbatim
!>        IntegralIn is type(mqc_scf_integral) 
!>        The MQC integral to be spin unblocked.
!>    \endverbatim
!>
!>    \param[out] MatrixOut
!>    \verbatim
!>        MatrixOut is type(mqc_matrix)
!>        The MQC matrix output.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \author A. D. Mahler
!>    \date 2017, 2018
!
      subroutine mqc_matrix_undoSpinBlockGHF_integral(integralIn,matrixOut)
!
      implicit none
      integer(kind=int64)::i,j,k,iter,rows,columns,nAlpha1,nBeta1,nElectrons1, &
        nAlpha2,nBeta2,nElectrons2
      type(mqc_scf_integral),intent(in)::integralIn
      type(mqc_matrix),intent(out)::matrixOut
      type(mqc_matrix)::tmpMatrix
      type(mqc_vector)::tmpVector
      integer(kind=int64),dimension(:),allocatable::elist
!
      nAlpha1 = integralIn%blockSize('alpha',1) 
      nBeta1 = integralIn%blockSize('beta',1) 
      nElectrons1 = nAlpha1+nBeta1
      nAlpha2 = integralIn%blockSize('alpha',2) 
      nBeta2 = integralIn%blockSize('beta',2) 
      nElectrons2 = nAlpha2+nBeta2
      call matrixOut%init(nElectrons1,nElectrons2)
      tmpMatrix = integralIn%getBlock('full')
      j = 1
      do i = 1,nAlpha1
        call matrixOut%vput(tmpMatrix%vat([i],[0]),[j],[0])
        j = j+2
      endDo
      j = 2
      do i = nAlpha1+1,nElectrons1 
        call matrixOut%vput(tmpMatrix%vat([i],[0]),[j],[0])
        j = j+2
      endDo
      tmpMatrix = matrixOut

      if(allocated(integralIn%Energy_List)) then
        elist = integralIn%getEList()
        columns = size(elist)
        DO I = 1, columns
          DO ITER = 1, columns
            if(elist(iter).eq.I) then
              call matrixOut%vput(tmpMatrix%vat([0],[iter]),[0],[i])
              exit
            endif
          END DO
        END DO
      else
        if(integralIn%getLabel().eq.'mo coefficients') &
          Write(6,'(A)') 'WARNING: Outputting general MO coefficients but Energy_List is not set. &
            &Assuming alternating alpha and Beta orbitals.' 
        j = 1
        do i = 1,nAlpha2
          call matrixOut%vput(tmpMatrix%vat([0],[i]),[0],[j])
          j = j+2
        endDo
        j = 2
        do i = nAlpha2+1,nElectrons2
          call matrixOut%vput(tmpMatrix%vat([0],[i]),[0],[j])
          j = j+2
        endDo
      endIf
!
      return
      end subroutine mqc_matrix_undoSpinBlockGHF_integral    
!
!
!     PROCEDURE MQC_SCF_Integral_Contraction
!
!>    \brief <b> MQC_SCF_Integral_Contraction is used to return the contraction of 
!>    two MQC integrals</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Integral_Contraction is used to return the contraction of two MQC 
!>    integrals. 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral1
!>    \verbatim
!>        Integral1 is type(mqc_scf_integral) 
!>        The first MQC integral to be contracted.
!>    \endverbatim
!>
!>    \param[in] Integral2
!>    \verbatim
!>        Integral2 is type(mqc_scf_integral) 
!>        The second MQC integral to be contracted.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_scf_integral_contraction(integral1,integral2) result(output)
!
      implicit none
      type(mqc_scf_integral),intent(in)::integral1,integral2
      type(mqc_scalar)::output
!
      output = 0.0
      if(integral1%blockSize('alpha',1).ne.integral2%blockSize('alpha',1)) &
        call mqc_error_I('MQC_SCF_Integral_Contraction has different size inputs',6, &
        'integral 1 alpha 1 block size:',integral1%blockSize('alpha',1), &
        'integral 2 alpha 1 block size:',integral2%blockSize('alpha',1))
      if(integral1%blockSize('beta',1).ne.integral2%blockSize('beta',1)) &
        call mqc_error_I('MQC_SCF_Integral_Contraction has different size inputs',6, &
        'integral 1 beta 1 block size:',integral1%blockSize('beta',1), &
        'integral 2 beta 1 block size:',integral2%blockSize('beta',1))
      if(integral1%blockSize('alpha',2).ne.integral2%blockSize('alpha',2)) &
        call mqc_error_I('MQC_SCF_Integral_Contraction has different size inputs',6, &
        'integral 1 alpha 2 block size:',integral1%blockSize('alpha',2), &
        'integral 2 alpha 2 block size:',integral2%blockSize('alpha',2))
      if(integral1%blockSize('beta',2).ne.integral2%blockSize('beta',2)) &
        call mqc_error_I('MQC_SCF_Integral_Contraction has different size inputs',6, &
        'integral 1 beta 2 block size:',integral1%blockSize('beta',2), &
        'integral 2 beta 2 block size:',integral2%blockSize('beta',2))
      output = output + contraction(integral1%alpha,integral2%alpha)
      if(integral1%type().eq.'space') then
        if(integral2%type().eq.'space') then
          output = output*2.0
        elseIf(integral2%type().eq.'spin'.or.integral2%type().eq.'general') then
              output = output + contraction(integral1%alpha,integral2%beta)
        endIf
      elseIf(integral1%type().eq.'spin'.or.integral1%type().eq.'general') then
        if(integral2%type().eq.'space') then
          output = output + contraction(integral1%beta,integral2%alpha)
        elseIf(integral2%type().eq.'spin'.or.integral2%type().eq.'general') then
          output = output + contraction(integral1%beta,integral2%beta)
        endIf
      endIf
      if(integral1%type().eq.'general'.and.integral2%type().eq.'general') then
        output = output + contraction(integral1%alphaBeta,integral2%alphaBeta)
        output = output + contraction(integral1%betaAlpha,integral2%betaAlpha)
      endIf
!
      end function mqc_scf_integral_contraction
!
!
!     PROCEDURE MQC_ERI_Integral_Contraction
!
!>    \brief <b> MQC_ERI_Integral_Contraction is used to return the contraction of 
!>    two-ERIs with an MQC integral</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_ERI_Integral_Contraction is used to return the contraction of two-ERIs with 
!>    an MQC integral. Optional argument label can be used to return contractions 
!>    with single bar ERIs (coulomb or exchange), or the complete contraction with 
!>    double-bar integrals. The following options are available:
!>
!>    1. Label = 'coulomb' returns the contraction with single bar coulomb integrals.
!>    2. Label = 'exchange' returns the contraction with single bar exchange integrals.
!>    3. Label = 'doublebar' returns the contraction with double-bar integrals (default).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] ERIs
!>    \verbatim
!>        ERIs is type(mqc_twoERIs) 
!>        The 2ERIs to be contracted.
!>    \endverbatim
!>
!>    \param[in] Integral
!>    \verbatim
!>        Integral is type(mqc_scf_integral) 
!>        The MQC integral to be contracted.
!>    \endverbatim
!>
!>    \param[in] ZeroIn
!>    \verbatim
!>        ZeroIn is real(kind=real64),optional
!>        The threshold value for performing contraction between integrals and 
!>        density matrix elements. 
!>    \endverbatim
!>    
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        = 'coulomb':   contracts with coulomb integrals.
!>        = 'exchange':  contracts with exchange integrals.
!>        = 'doublebar': contracts with double-bar integrals 
!>                       (default).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_eri_integral_contraction(eris,integral,zeroIn,label) result(output)
!
      implicit none
      type(mqc_twoERIs),dimension(:),intent(in)::eris
      type(mqc_scf_integral),intent(in)::integral
      real(kind=real64),optional,intent(in)::zeroIn
      character(len=*),optional::label
      type(mqc_scf_integral)::output
      type(mqc_matrix)::alpha,beta,alphaBeta,betaAlpha
      integer(kind=int64)::cnt,kp,mp,k,l,m,n,nBasis,ntt,fact1,flag,i
      character(len=64)::mylabel,algorithm='regular'
      real(kind=real64)::zero
!     Pgfortran does not support declare reduce AND no way to restrict reduction to only 
!     defined arrays. For now define local copies by hand and then sum into global version
!     in critical section.
      type(mqc_matrix)::alphaLoc,betaLoc,alphaBetaLoc,betaAlphaLoc
      type(mqc_scf_integral)::hermden,ahermden
      type(mqc_scalar)::half
      real(kind=real64)::raf1,raf2,raf3
      logical::DEBUG=.false.,doHerm,doAhrm
      real(kind=real64),dimension(:),allocatable::aaHRden,bbHRden,aaHCden,bbHCden,abHHRden, &
        abAHRden,abHHCden,abAHCden,aaARden,bbARden,aaACden,bbACden,abHARden,abAARden, &
        abHACden,abAACden
      real(kind=real64),dimension(:),allocatable::aaHRfoc,bbHRfoc,aaHCfoc,bbHCfoc,abHHRfoc, &
        abAHRfoc,abHHCfoc,abAHCfoc,aaARfoc,bbARfoc,aaACfoc,bbACfoc,abHARfoc,abAARfoc, &
        abHACfoc,abAACfoc
      real(kind=real64),dimension(:,:),allocatable::temp
      type(mqc_matrix)::tmpMat,tmpMat2,tmpMat3
!
      half = 0.5
      if(algorithm.ne.'square'.and.algorithm.ne.'regular') &
        call mqc_error_A('Unknown algorithm option in mqc_eri_integral_contraction',6, &
        'algorithm',algorithm)
!
      if(present(zeroIn)) then
        zero = zeroIn
      else
        zero = 1.0d-10
      endIf
!
      if(present(label)) then
        call string_change_case(label,'L',myLabel)
        if(algorithm.eq.'regular'.and.myLabel.ne.'doublebar') &
          call mqc_error_a('Regular algorithm only returns double bar integrals',6, &
          'algorithm',algorithm,'myLabel',myLabel)
      else
        myLabel = 'doublebar'
      endIf

      if(integral%blockSize('alpha',1).ne.integral%blockSize('alpha',2)) &
        call mqc_error_I('Arrays are not sized as expected in MQC_ERI_Integral_Contraction', &
        6,'integral%blockSize(alpha,1)', integral%blockSize('alpha',1), &
        'integral%blockSize(alpha,2)',integral%blockSize('alpha',2))
      if(integral%blockSize('beta',1).ne.integral%blockSize('beta',2)) &
        call mqc_error_I('Arrays are not sized as expected in MQC_ERI_Integral_Contraction' &
        ,6,'integral%blockSize(beta,1)', integral%blockSize('beta',1), &
        'integral%blockSize(beta,2)',integral%blockSize('beta',2))
      if(integral%blockSize('alpha',1).ne.integral%blockSize('beta',2)) &
        call mqc_error_I('Arrays are not sized as expected in MQC_ERI_Integral_Contraction', &
        6,'integral%blockSize(alpha,1)', integral%blockSize('alpha',1), &
        'integral%blockSize(beta,2)',integral%blockSize('beta',2))
      if(integral%blockSize('beta',1).ne.integral%blockSize('alpha',2)) &
        call mqc_error_I('Arrays are not sized as expected in MQC_ERI_Integral_Contraction', &
        6,'integral%blockSize(beta,1)', integral%blockSize('beta',1), &
        'integral%blockSize(alpha,2)',integral%blockSize('alpha',2))
      if(myLabel.ne.'doublebar'.and.algorithm.eq.'regular') &
        call mqc_error('Regular Fock build cannot return separate interaction energies')

      flag = 0
      if(integral%array_type.eq.'spin') then
        flag = flag + 1
      elseIf(integral%array_type.eq.'general') then
        flag = flag + 2
      endIf
      if(mqc_have_complex(integral%getBlock('full'))) flag = flag + 10
      if(eris(1)%type().eq.'space'.or.eris(1)%type().eq.'spin'.or.eris(1)%type().eq.'general') &
        call mqc_error_a('Trying to contract matrix with MO basis 2ERIs',6,'eris(1)%type()',eris(1)%type())
      if(eris(1)%type().eq.'regular') then
        flag = flag + 100
      elseif(eris(1)%type().eq.'raffenetti1') then
        flag = flag + 200
      elseif(eris(1)%type().eq.'raffenetti2') then
        flag = flag + 300
      elseif(eris(1)%type().eq.'raffenetti3') then
        flag = flag + 400
      endIf
      if(size(eris).ge.2) then
        if(eris(2)%type().eq.'regular') then
          flag = flag + 1000
        elseif(eris(2)%type().eq.'raffenetti1') then
          flag = flag + 2000
        elseif(eris(2)%type().eq.'raffenetti2') then
          flag = flag + 3000
        elseif(eris(2)%type().eq.'raffenetti3') then
          flag = flag + 4000
        endIf
      endIf
      if(size(eris).ge.3) then
        if(eris(3)%type().eq.'regular') then
          flag = flag + 10000
        elseif(eris(3)%type().eq.'raffenetti1') then
          flag = flag + 20000
        elseif(eris(3)%type().eq.'raffenetti2') then
          flag = flag + 30000
        elseif(eris(3)%type().eq.'raffenetti3') then
          flag = flag + 40000
        endIf
      endIf
      
      if(debug) write(6,'(A,1x,I5)') 'Option flag:',flag

      nBasis = integral%blockSize('alpha')
      ntt = (nBasis*(nBasis+1))/2
      call alpha%init(nBasis,nBasis)
      if(integral%array_type.eq.'spin'.or.integral%array_type.eq.'general') &
        call beta%init(nBasis,nBasis)
      if(integral%array_type.eq.'general') then
        call alphaBeta%init(nBasis,nBasis)
        call betaAlpha%init(nBasis,nBasis)
      endIf

      if(algorithm.eq.'regular') then

        doHerm = .false.
        doAhrm = .false.
        
        hermden = half*(integral+dagger(integral))
        if(hermden%norm('M').gt.zero) doHerm = .true.
        if(doHerm) then
          tmpMat = real(hermden%getblock('alpha'))
          if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
          if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
          temp = tmpMat
          aaHRden = temp(:,1)
          if(mod(flag,100)/10.eq.1) then
            tmpMat = aimag(hermden%getblock('alpha'))
            if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
            if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
            temp = tmpMat
            aaHCden = temp(:,1)
            aaHCden = (-1)*aaHCden
          endIf
          if(mod(flag,10).gt.0) then
            tmpMat = real(hermden%getblock('beta'))
            if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
            if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
            temp = tmpMat
            bbHRden = temp(:,1)
            if(mod(flag,100)/10.eq.1) then
              tmpMat = aimag(hermden%getblock('beta'))
              if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
              if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
              temp = tmpMat
              bbHCden = temp(:,1)
              bbHCden = (-1)*bbHCden
            endIf
          endIf
          if(mod(flag,10).eq.2) then
            tmpMat = real(half*(hermden%getblock('alpha-beta')+ &
              dagger(hermden%getblock('alpha-beta'))))
            if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
            if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
            temp = tmpMat
            abHHRden = temp(:,1)
            tmpMat = real(half*(hermden%getblock('alpha-beta')- &
              dagger(hermden%getblock('alpha-beta'))))
            if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
            if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
            temp = tmpMat
            abAHRden = temp(:,1)
            if(mod(flag,100)/10.eq.1) then
              tmpMat = aimag(half*(hermden%getblock('alpha-beta')+ &
                dagger(hermden%getblock('alpha-beta'))))
              if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
              if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
              temp = tmpMat
              abHHCden = temp(:,1)
              abHHCden = (-1)*abHHCden
              tmpMat = aimag(half*(hermden%getblock('alpha-beta')- &
                dagger(hermden%getblock('alpha-beta'))))
              if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
              if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
              temp = tmpMat
              abAHCden = temp(:,1)
            endIf
          endIf
        endIf
        
        ahermden = half*(integral-dagger(integral))
        if(ahermden%norm('M').gt.zero) doAhrm = .true.
        if(doAhrm) then
          tmpMat = real(ahermden%getblock('alpha'))
          if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
          if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
          temp = tmpMat
          aaARden = temp(:,1)
          if(mod(flag,100)/10.eq.1) then
            tmpMat = aimag(ahermden%getblock('alpha'))
            if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
            if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
            temp = tmpMat
            aaACden = temp(:,1)
          endIf
          if((mod(flag,10).gt.0)) then
            tmpMat = real(ahermden%getblock('beta'))
            if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
            if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
            temp = tmpMat
            bbARden = temp(:,1)
            if(mod(flag,100)/10.eq.1) then
              tmpMat = aimag(ahermden%getblock('beta'))
              if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
              if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
              temp = tmpMat
              bbACden = temp(:,1)
            endIf
          endIf
          if((mod(flag,10).eq.2)) then
            tmpMat = real(half*(ahermden%getblock('alpha-beta')+ &
              dagger(ahermden%getblock('alpha-beta'))))
            if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
            if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
            temp = tmpMat
            abHARden = temp(:,1)
            tmpMat = real(half*(ahermden%getblock('alpha-beta')- &
              dagger(ahermden%getblock('alpha-beta'))))
            if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
            if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
            temp = tmpMat
            abAARden = temp(:,1)
            if(mod(flag,100)/10.eq.1) then
              tmpMat = aimag(half*(ahermden%getblock('alpha-beta')+ &
                dagger(ahermden%getblock('alpha-beta'))))
              if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
              if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
              temp = tmpMat
              abHACden = temp(:,1)
              tmpMat = aimag(half*(ahermden%getblock('alpha-beta')- &
                dagger(ahermden%getblock('alpha-beta'))))
              if(mqc_matrix_haveFull(tmpMat)) call mqc_matrix_full2symm(tmpMat)
              if(mqc_matrix_haveDiagonal(tmpMat)) call mqc_matrix_diag2symm(tmpMat)
              temp = tmpMat
              abAACden = temp(:,1)
            endIf
          endIf
        endIf

!       if integral is zero then output is also zero
        if(.not.doHerm.and..not.doAhrm) then
          output = integral
          return
        endIf

        if(DEBUG) then
          if(doHerm) then 
            call mqc_print(aaHRden,6,'Alpha-Alpha hermitian real density')
            if(mod(flag,100)/10.eq.1) &
              call mqc_print(aaHCden,6,'Alpha-Alpha hermitian complex density')
            if (mod(flag,10).gt.0) then
              call mqc_print(bbHRden,6,'Beta-Beta hermitian real density')
              if(mod(flag,100)/10.eq.1) &
                call mqc_print(bbHCden,6,'Beta-Beta hermitian complex density')
            endIf
            if (mod(flag,10).eq.2) then
              call mqc_print(abHHRden,6,'Alpha-Beta hermitian hermitian real density')
              call mqc_print(abAHRden,6,'Alpha-Beta hermitian antihermitian real density')
              if(mod(flag,100)/10.eq.1) then
                call mqc_print(abHHCden,6,'Alpha-Beta hermitian hermitian complex density')
                call mqc_print(abAHCden,6,'Alpha-Beta hermitian antihermitian complex density')
              endIf
            endIf
          endIf
          if(doAhrm) then 
            call mqc_print(aaARden,6,'Alpha-Alpha antihermitian real density')
            if(mod(flag,100)/10.eq.1) &
              call mqc_print(aaACden,6,'Alpha-Alpha antihermitian complex density')
            if (mod(flag,10).gt.0) then
              call mqc_print(bbARden,6,'Beta-Beta antihermitian real density')
              if(mod(flag,100)/10.eq.1) &
                call mqc_print(bbACden,6,'Beta-Beta antihermitian complex density')
            endIf
            if (mod(flag,10).eq.2) then
              call mqc_print(abHARden,6,'Alpha-Beta hermitian antihermitian real density')
              call mqc_print(abAARden,6,'Alpha-Beta antihermitian antihermitian real density')
              if(mod(flag,100)/10.eq.1) then
                call mqc_print(abHACden,6,'Alpha-Beta hermitian antihermitian complex density')
                call mqc_print(abAACden,6,'Alpha-Beta antihermitian antihermitian complex density')
              endIf
            endIf
          endIf
        endIf

        if(doHerm) then
          if (mod(flag,10).eq.0) then 
            aaHRden = 4*aaHRden 
            do i = 1,nBasis
              aaHRden(symIndexHash(i,i)) = aaHRden(symIndexHash(i,i))/2.0
            endDo
          else
            aaHRden = 2*(aaHRden + bbHRden) 
            bbHRden = 2*(aaHRden/2.0 - 2*bbHRden) 
            do i = 1,nBasis
              aaHRden(symIndexHash(i,i)) = aaHRden(symIndexHash(i,i))/2.0
              bbHRden(symIndexHash(i,i)) = bbHRden(symIndexHash(i,i))/2.0
            endDo
          endIf
          if (mod(flag,10).eq.2) then
            do i = 1,nBasis
              abHHRden(symIndexHash(i,i)) = abHHRden(symIndexHash(i,i))/2.0
              if(mod(flag,100)/10.eq.1) &
                abAHCden(symIndexHash(i,i)) = abAHCden(symIndexHash(i,i))/2.0
            endDo
          endIf
        endIf

        if(doAhrm) then
          aaACden = 2*(aaACden + bbACden) 
          bbACden = 2*(aaACden/2.0 - 2*bbACden) 
          do i = 1,nBasis
            if(mod(flag,10).eq.2) &
              abHARden(symIndexHash(i,i)) = abHARden(symIndexHash(i,i))/2.0
            if(mod(flag,100)/10.eq.1) then
              aaACden(symIndexHash(i,i)) = aaACden(symIndexHash(i,i))/2.0
              if(mod(flag,10).gt.0) &
                bbACden(symIndexHash(i,i)) = bbACden(symIndexHash(i,i))/2.0
              if(mod(flag,10).eq.2) &
                abAACden(symIndexHash(i,i)) = abAACden(symIndexHash(i,i))/2.0
            endIf
          endDo
        endIf

!       Logic here is commented out as cannot do reduction over unallocated arrays in
!       OMP section below. Therefore allocate and initialize all components of output 
!       matrix even though some will never be used.
!
!        if(doHerm) then
          allocate(aaHRfoc((nBasis*(nBasis+1))/2))
          aaHRfoc = 0.0
!          if(mod(flag,100)/10.eq.1) then
            allocate(aaHCfoc((nBasis*(nBasis+1))/2))
            aaHCfoc = 0.0
!          endIf
!          if(mod(flag,10).gt.0) then
            allocate(bbHRfoc((nBasis*(nBasis+1))/2))
            bbHRfoc = 0.0 
!            if(mod(flag,100)/10.eq.1) then
              allocate(bbHCfoc((nBasis*(nBasis+1))/2))
              bbHCfoc = 0.0 
!            endIf
!            if(mod(flag,10).eq.2) then
              allocate(abHHRfoc((nBasis*(nBasis+1))/2))
              abHHRfoc = 0.0
!              if(mod(flag,100)/10.eq.1) then
                allocate(abHHCfoc((nBasis*(nBasis+1))/2))
                abHHCfoc = 0.0
!              endIf
              allocate(abAHRfoc((nBasis*(nBasis+1))/2))
              abAHRfoc = 0.0
!              if(mod(flag,100)/10.eq.1) then
                allocate(abAHCfoc((nBasis*(nBasis+1))/2))
                abAHCfoc = 0.0
!              endIf
!            endIf
!          endIf
!        endIf

!        if(doAhrm) then
          allocate(aaARfoc((nBasis*(nBasis+1))/2))
          aaARfoc = 0.0
!          if(mod(flag,100)/10.eq.1) then
            allocate(aaACfoc((nBasis*(nBasis+1))/2))
            aaACfoc = 0.0
!          endIf
!          if(mod(flag,10).gt.0) then
            allocate(bbARfoc((nBasis*(nBasis+1))/2))
            bbARfoc = 0.0
!            if(mod(flag,100)/10.eq.1) then
              allocate(bbACfoc((nBasis*(nBasis+1))/2))
              bbACfoc = 0.0 
!            endIf
!            if(mod(flag,10).eq.2) then
              allocate(abHARfoc((nBasis*(nBasis+1))/2))
              abHARfoc = 0.0
!              if(mod(flag,100)/10.eq.1) then
                allocate(abHACfoc((nBasis*(nBasis+1))/2))
                abHACfoc = 0.0
!              endIf
              allocate(abAARfoc((nBasis*(nBasis+1))/2))
              abAARfoc = 0.0
!              if(mod(flag,100)/10.eq.1) then
                allocate(abAACfoc((nBasis*(nBasis+1))/2))
                abAACfoc = 0.0
!              endIf
!            endIf
!          endIf
!        endIf

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP REDUCTION(+:aaHRfoc,bbHRfoc,aaHCfoc,bbHCfoc,abHHRfoc,abAHRfoc,abHHCfoc,abAHCfoc, &
!$OMP aaARfoc,bbARfoc,aaACfoc,bbACfoc,abHARfoc,abAARfoc,abHACfoc,abAACfoc), &
!$OMP SHARED(eris,ntt,aaHRden,bbHRden,aaHCden,bbHCden,abHHRden,abAHRden,abHHCden,abAHCden, &
!$OMP aaARden,bbARden,aaACden,bbACden,abHARden,abAARden,abHACden,abAACden,flag,doHerm, &
!$OMP doAhrm,zero), &
!$OMP PRIVATE(cnt,kp,mp,k,l,m,n,raf1,raf2,raf3)
        do cnt = 1, (ntt*(ntt+1))/2
#ifdef   PGI
          kp = floor((1+sqrt(dfloat(-7+8*cnt)))/2)
          mp = cnt - kp*(kp-1)/2
          k = floor((1+sqrt(dfloat(-7+8*kp)))/2) 
          m = floor((1+sqrt(dfloat(-7+8*mp)))/2) 
#else
          kp = floor((1+sqrt(float(-7+8*cnt)))/2)
          mp = cnt - kp*(kp-1)/2
          k = floor((1+sqrt(float(-7+8*kp)))/2) 
          m = floor((1+sqrt(float(-7+8*mp)))/2) 
#endif  
          l = kp - k*(k-1)/2
          n = mp - m*(m-1)/2

          if(mod(flag,1000)/100.eq.1) then
            raf1 = eris(1)%at(k,l,m,n) - 0.25*(eris(1)%at(k,n,m,l) + eris(1)%at(k,m,n,l))
            if(kp.eq.mp) raf1 = 0.5*raf1 
          elseIf(mod(flag,1000)/100.eq.2) then
            raf1 = eris(1)%at(k,l,m,n)
          else
            call mqc_error_a('Unrecognized/incorrect integrals given in eris(1) in &
              & mqc_eri_integral_contraction',6,'eris(1)%type()',eris(1)%type())
          endIf
          if(mod(flag,10).gt.0) then
            if(mod(flag,1000)/100.eq.1) then
              raf2 = eris(1)%at(k,n,m,l) + eris(1)%at(k,m,n,l)
              if(kp.eq.mp) raf2 = 0.5*raf2 
            elseIf(mod(flag,10000)/1000.eq.3) then
              raf2 = eris(2)%at(k,l,m,n) 
            else
              call mqc_error_a('Unrecognized/incorrect integrals given in eris(2) in &
                & mqc_eri_integral_contraction',6,'eris(2)%type()',eris(2)%type())
            endIf
          else
            raf2 = 0.0 
          endIf
          if (m.ne.n) then
            if (doAhrm.or.mod(flag,10).eq.2.or.mod(flag,100)/10.eq.1) then
              if(mod(flag,1000)/100.eq.1) then
                raf3 = eris(1)%at(k,m,n,l) - eris(1)%at(k,n,m,l) 
                if(kp.eq.mp) raf3 = 0.5*raf3 
              elseIf(mod(flag,10000)/1000.eq.4) then
                raf3 = eris(2)%at(k,l,m,n) 
              elseIf(mod(flag,100000)/10000.eq.4) then
                raf3 = eris(3)%at(k,l,m,n) 
              else
                call mqc_error_a('Unrecognized/incorrect integrals given in eris(2) or & 
                  & eris(3) in mqc_eri_integral_contraction',6,'eris(2)%type()', &
                  eris(2)%type(),'eris(3)%type()',eris(3)%type())
              endIf
            endIf
          else
            raf3 = 0.0
          endIf

          if(doHerm) then
!         AA/BB Hermitian Fock
            if(abs(raf1).gt.zero.and.abs(aaHRden(symIndexHash(m,n))).gt.zero) &
              aaHRfoc(symIndexHash(k,l)) = aaHRfoc(symIndexHash(k,l)) + &
              raf1*aaHRden(symIndexHash(m,n))
            if(abs(raf1).gt.zero.and.abs(aaHRden(symIndexHash(k,l))).gt.zero) &
              aaHRfoc(symIndexHash(m,n)) = aaHRfoc(symIndexHash(m,n)) + &
              raf1*aaHRden(symIndexHash(k,l))
            if(mod(flag,10).gt.0) then
              if(abs(raf2).gt.zero.and.abs(bbHRden(symIndexHash(m,n))).gt.zero) &
                bbHRfoc(symIndexHash(k,l)) = bbHRfoc(symIndexHash(k,l)) + &
                raf2*bbHRden(symIndexHash(m,n))
              if(abs(raf2).gt.zero.and.abs(bbHRden(symIndexHash(k,l))).gt.zero) &
                bbHRfoc(symIndexHash(m,n)) = bbHRfoc(symIndexHash(m,n)) + &
                raf2*bbHRden(symIndexHash(k,l))
            endIf
            if(m.ne.n.and.mod(flag,100)/10.eq.1) then
              if(abs(raf3).gt.zero) then
                if(abs(aaHCden(symIndexHash(m,n))).gt.zero) &
                  aaHCfoc(symIndexHash(k,l)) = aaHCfoc(symIndexHash(k,l)) - &
                  raf3*aaHCden(symIndexHash(m,n))
                if(abs(aaHCden(symIndexHash(k,l))).gt.zero) &
                  aaHCfoc(symIndexHash(m,n)) = aaHCfoc(symIndexHash(m,n)) - &
                  raf3*aaHCden(symIndexHash(k,l))
                if(mod(flag,10).gt.0) then
                  if(abs(bbHCden(symIndexHash(m,n))).gt.zero) &
                    bbHCfoc(symIndexHash(k,l)) = bbHCfoc(symIndexHash(k,l)) - &
                    raf3*bbHCden(symIndexHash(m,n))
                  if(abs(bbHCden(symIndexHash(k,l))).gt.zero) &
                    bbHCfoc(symIndexHash(m,n)) = bbHCfoc(symIndexHash(m,n)) - &
                    raf3*bbHCden(symIndexHash(k,l))
                endIf
              endIf
            endIf
            if(mod(flag,10).eq.2) then
!         AB Hermitian Hermitian Fock
              if(abs(abHHRden(symIndexHash(m,n))).gt.zero) &
                abHHRfoc(symIndexHash(k,l)) = abHHRfoc(symIndexHash(k,l)) - &
                raf2*abHHRden(symIndexHash(m,n))
              if(abs(abHHRden(symIndexHash(k,l))).gt.zero) &
                abHHRfoc(symIndexHash(m,n)) = abHHRfoc(symIndexHash(m,n)) - &
                raf2*abHHRden(symIndexHash(k,l))
!              if(m.eq.n.and.abs(abAHRden(symIndexHash(k,l))).gt.zero) &
!                abHHRfoc(symIndexHash(m,n)) = abHHRfoc(symIndexHash(m,n)) - &
!                raf2*abAHRden(symIndexHash(k,l))
!              if(k.eq.l.and.abs(abAHRden(symIndexHash(m,n))).gt.zero) then
!                if(k.eq.m) then
!                  abHHRfoc(symIndexHash(k,l)) = abHHRfoc(symIndexHash(k,l)) + &
!                    raf2*abAHRden(symIndexHash(m,n))
!                else
!                  abHHRfoc(symIndexHash(k,l)) = abHHRfoc(symIndexHash(k,l)) - &
!                    raf2*abAHRden(symIndexHash(m,n))
!                endIf
!              endIf
              if(m.ne.n.and.mod(flag,100)/10.eq.1) then
                if(abs(abHHCden(symIndexHash(m,n))).gt.zero) &
                  abHHCfoc(symIndexHash(k,l)) = abHHCfoc(symIndexHash(k,l)) + &
                  raf3*abHHCden(symIndexHash(m,n))
                if(abs(abHHCden(symIndexHash(k,l))).gt.zero) &
                  abHHCfoc(symIndexHash(m,n)) = abHHCfoc(symIndexHash(m,n)) + &
                  raf3*abHHCden(symIndexHash(k,l))
              endIf
!         AB Antihermitian Hermitian Fock
              if(m.ne.n.and.abs(abAHRden(symIndexHash(m,n))).gt.zero) &
                abAHRfoc(symIndexHash(k,l)) = abAHRfoc(symIndexHash(k,l)) + &
                raf3*abAHRden(symIndexHash(m,n))
              if(m.ne.n.and.abs(abAHRden(symIndexHash(k,l))).gt.zero) &
                abAHRfoc(symIndexHash(m,n)) = abAHRfoc(symIndexHash(m,n)) + &
                raf3*abAHRden(symIndexHash(k,l))
              if (mod(flag,100)/10.eq.1) then
                if(abs(abAHCden(symIndexHash(m,n))).gt.zero) &
                  abAHCfoc(symIndexHash(k,l)) = abAHCfoc(symIndexHash(k,l)) - &
                  raf2*abAHCden(symIndexHash(m,n))
                if(abs(abAHCden(symIndexHash(k,l))).gt.zero) &
                  abAHCfoc(symIndexHash(m,n)) = abAHCfoc(symIndexHash(m,n)) - &
                  raf2*abAHCden(symIndexHash(k,l))
              endIf
            endIf
          endIf
!         AA Antihermitian Fock
          if(doAhrm) then
            if(m.ne.n.and.abs(aaARden(symIndexHash(m,n))).gt.zero) &
              aaARfoc(symIndexHash(k,l)) = aaARfoc(symIndexHash(k,l)) + &
              raf3*aaARden(symIndexHash(m,n))
            if(m.ne.n.and.abs(aaARden(symIndexHash(k,l))).gt.zero) &
              aaARfoc(symIndexHash(m,n)) = aaARfoc(symIndexHash(m,n)) + &
              raf3*aaARden(symIndexHash(k,l))
            if(mod(flag,100)/10.eq.1) then
              if(abs(aaACden(symIndexHash(m,n))).gt.zero) &
                aaACfoc(symIndexHash(k,l)) = aaACfoc(symIndexHash(k,l)) + &
                raf1*aaACden(symIndexHash(m,n))
              if(abs(aaACden(symIndexHash(k,l))).gt.zero) &
                aaACfoc(symIndexHash(m,n)) = aaACfoc(symIndexHash(m,n)) + &
                raf1*aaACden(symIndexHash(k,l))
            endIf
            if(mod(flag,10).gt.0) then
!         BB Antihermitian Fock
              if(m.ne.n.and.abs(bbARden(symIndexHash(m,n))).gt.zero) &
                bbARfoc(symIndexHash(k,l)) = bbARfoc(symIndexHash(k,l)) + &
                raf3*bbARden(symIndexHash(m,n))
              if(m.ne.n.and.abs(bbARden(symIndexHash(k,l))).gt.zero) &
                bbARfoc(symIndexHash(m,n)) = bbARfoc(symIndexHash(m,n)) + &
                raf3*bbARden(symIndexHash(k,l))
              if(mod(flag,100)/10.eq.1) then
                if(abs(bbACden(symIndexHash(m,n))).gt.zero) &
                  bbACfoc(symIndexHash(k,l)) = bbACfoc(symIndexHash(k,l)) + &
                  raf2*bbACden(symIndexHash(m,n))
                if(abs(bbACden(symIndexHash(k,l))).gt.zero) &
                  bbACfoc(symIndexHash(m,n)) = bbACfoc(symIndexHash(m,n)) + &
                  raf2*bbACden(symIndexHash(k,l))
              endIf
            endIf
!         AB Hermitian Antihermitian Fock
            if(mod(flag,10).eq.2) then
              if(abs(abHARden(symIndexHash(m,n))).gt.zero) &
                abHARfoc(symIndexHash(k,l)) = abHARfoc(symIndexHash(k,l)) - &
                raf2*abHARden(symIndexHash(m,n))
              if(abs(abHARden(symIndexHash(k,l))).gt.zero) &
                abHARfoc(symIndexHash(m,n)) = abHARfoc(symIndexHash(m,n)) - &
                raf2*abHARden(symIndexHash(k,l))
              if(mod(flag,100)/10.eq.1.and.m.ne.n) then
                if(abs(abHACden(symIndexHash(m,n))).gt.zero) &
                 abHACfoc(symIndexHash(k,l)) = abHACfoc(symIndexHash(k,l)) + &
                 raf3*abHACden(symIndexHash(m,n))
                if(abs(abHACden(symIndexHash(k,l))).gt.zero) &
                 abHACfoc(symIndexHash(m,n)) = abHACfoc(symIndexHash(m,n)) + &
                 raf3*abHACden(symIndexHash(k,l))
             endIf
!         AB Antihermitian Antihermitian Fock
              if(m.ne.n.and.abs(abAARden(symIndexHash(m,n))).gt.zero) &
                abAARfoc(symIndexHash(k,l)) = abAARfoc(symIndexHash(k,l)) + &
                raf3*abAARden(symIndexHash(m,n))
              if(m.ne.n.and.abs(abAARden(symIndexHash(k,l))).gt.zero) &
                abAARfoc(symIndexHash(m,n)) = abAARfoc(symIndexHash(m,n)) + &
                raf3*abAARden(symIndexHash(k,l))
              if(mod(flag,100)/10.eq.1) then
                if(abs(abAACden(symIndexHash(m,n))).gt.zero) &
                  abAACfoc(symIndexHash(k,l)) = abAACfoc(symIndexHash(k,l)) - &
                  raf2*abAACden(symIndexHash(m,n))
                if(abs(abAACden(symIndexHash(k,l))).gt.zero) &
                  abAACfoc(symIndexHash(m,n)) = abAACfoc(symIndexHash(m,n)) - &
                  raf2*abAACden(symIndexHash(k,l))
              endIf
            endIf
          endIf

        endDo
!$OMP END PARALLEL DO

        if(mod(flag,10).gt.0) then
          aaHRfoc = aaHRfoc - 0.25*bbHRfoc
          bbHRfoc = aaHRfoc + 0.5*bbHRfoc
          aaACfoc = aaACfoc - 0.25*bbACfoc
          bbACfoc = aaACfoc + 0.5*bbACfoc
        endIf

        if(DEBUG) then
          if(doHerm) then
            call mqc_print(aaHRfoc,6,'Alpha-Alpha hermitian real Fock')
            if(mod(flag,100)/10.eq.1) &
              call mqc_print((-1)*aaHCfoc,6,'Alpha-Alpha hermitian complex Fock')
            if(mod(flag,10).gt.0) then
              call mqc_print(bbHRfoc,6,'Beta-Beta hermitian real Fock')
              if(mod(flag,100)/10.eq.1) &
                call mqc_print((-1)*bbHCfoc,6,'Beta-Beta hermitian complex Fock')
              if(mod(flag,10).eq.2) then
                call mqc_print(abHHRfoc,6,'Alpha-Beta hermitian hermitian real Fock')
                if(mod(flag,100)/10.eq.1) &
                  call mqc_print(abHHCfoc,6,'Alpha-Beta hermitian hermitian complex Fock')
                call mqc_print(abAHRfoc,6,'Alpha-Beta hermitian antihermitian real Fock')
                if(mod(flag,100)/10.eq.1) &
                  call mqc_print(abAHCfoc,6,'Alpha-Beta hermitian antihermitian complex Fock')
              endIf
            endIf
          endIf
          if(doAhrm) then
            call mqc_print(aaARfoc,6,'Alpha-Alpha antihermitian real Fock')
            if(mod(flag,100)/10.eq.1) &
              call mqc_print(aaACfoc,6,'Alpha-Alpha antihermitian complex Fock')
            if(mod(flag,10).gt.0) then
              call mqc_print(bbARfoc,6,'Beta-Beta antihermitian real Fock')
              if(mod(flag,100)/10.eq.1) &
                call mqc_print(bbACfoc,6,'Beta-Beta antihermitian complex Fock')
              if(mod(flag,10).eq.2) then
                call mqc_print(abHARfoc,6,'Alpha-Beta antihermitian hermitian real Fock')
                if(mod(flag,100)/10.eq.1) &
                  call mqc_print(abHACfoc,6,'Alpha-Beta antihermitian hermitian complex Fock')
                call mqc_print(abAARfoc,6,'Alpha-Beta antihermitian antihermitian real Fock')
                if(mod(flag,100)/10.eq.1) &
                  call mqc_print(abAACfoc,6,'Alpha-Beta antihermitian antihermitian complex Fock')
              endIf
            endIf
          endIf
        endIf

        if(doHerm.and.doAhrm) then
          if(mod(flag,100)/10.eq.0) then
            call MQC_Matrix_SymmMatrix_Put(alpha,aaHRfoc,'symmetric')
            call MQC_Matrix_SymmMatrix_Put(tmpMat,aaARfoc,'antisymmetric')
            alpha = alpha + tmpMat
            if(mod(flag,10).gt.0) then
              call MQC_Matrix_SymmMatrix_Put(beta,bbHRfoc,'symmetric')
              call MQC_Matrix_SymmMatrix_Put(tmpMat,bbARfoc,'antisymmetric')
              beta = beta + tmpMat 
              if(mod(flag,10).eq.2) then
                call MQC_Matrix_SymmMatrix_Put(alphaBeta,abHHRfoc,'symmetric')
                call MQC_Matrix_SymmMatrix_Put(tmpMat,abAHRfoc,'antisymmetric')
                call MQC_Matrix_SymmMatrix_Put(tmpMat2,abHARfoc,'symmetric')
                call MQC_Matrix_SymmMatrix_Put(tmpMat3,abAARfoc,'antisymmetric')
                betaAlpha = alphaBeta + tmpMat - tmpMat2 + tmpMat3
                alphaBeta = alphaBeta - tmpMat + tmpMat2 + tmpMat3
              endIf
            endIf
          elseIf(mod(flag,100)/10.eq.1) then
            call MQC_Matrix_SymmMatrix_Put(alpha,cmplx(aaHRfoc,(-1)*aaHCfoc),'hermitian')
            call MQC_Matrix_SymmMatrix_Put(tmpMat,cmplx(aaARfoc,aaACfoc),'antihermitian')
            alpha = alpha + tmpMat
            if(mod(flag,10).gt.0) then
              call MQC_Matrix_SymmMatrix_Put(beta,cmplx(bbHRfoc,(-1)*bbHCfoc),'hermitian')
              call MQC_Matrix_SymmMatrix_Put(tmpMat,cmplx(bbARfoc,bbACfoc),'antihermitian')
              beta = beta + tmpMat
              if(mod(flag,10).eq.2) then
                call MQC_Matrix_SymmMatrix_Put(alphaBeta,cmplx(abHHRfoc,(-1)*abHHCfoc),'hermitian')
                call MQC_Matrix_SymmMatrix_Put(tmpMat,cmplx(abAHRfoc,(-1)*abAHCfoc),'antihermitian')
                call MQC_Matrix_SymmMatrix_Put(tmpMat2,cmplx(abHARfoc,abHACfoc),'hermitian')
                call MQC_Matrix_SymmMatrix_Put(tmpMat3,cmplx(abAARfoc,abAACfoc),'antihermitian')
                betaAlpha = alphaBeta + tmpMat - tmpMat2 + tmpMat3 
                alphaBeta = alphaBeta - tmpMat + tmpMat2 + tmpMat3 
              endIf
            endIf
          endIf
        elseIf(doHerm.and..not.doAhrm) then
          if(mod(flag,100)/10.eq.0) then
            call MQC_Matrix_SymmMatrix_Put(alpha,aaHRfoc,'symmetric')
            if(mod(flag,10).gt.0) then
              call MQC_Matrix_SymmMatrix_Put(beta,bbHRfoc,'symmetric')
              if(mod(flag,10).eq.2) then
                call MQC_Matrix_SymmMatrix_Put(alphaBeta,abHHRfoc,'symmetric')
                call MQC_Matrix_SymmMatrix_Put(tmpMat,abAHRfoc,'antisymmetric')
                betaAlpha = alphaBeta + tmpMat
                alphaBeta = alphaBeta - tmpMat
              endIf
            endIf
          elseIf(mod(flag,100)/10.eq.1) then
            call MQC_Matrix_SymmMatrix_Put(alpha,cmplx(aaHRfoc,aaHCfoc),'hermitian') 
            alpha = transpose(alpha)
            if(mod(flag,10).gt.0) then
               call MQC_Matrix_SymmMatrix_Put(beta,cmplx(bbHRfoc,bbHCfoc),'hermitian') 
               beta = transpose(beta)
              if(mod(flag,10).eq.2) then
                call MQC_Matrix_SymmMatrix_Put(alphaBeta,cmplx(abHHRfoc,abHHCfoc),'hermitian')
                call MQC_Matrix_SymmMatrix_Put(tmpMat,cmplx(abAHRfoc,(-1)*abAHCfoc),'antihermitian')
                betaAlpha = alphaBeta + tmpMat
                alphaBeta = alphaBeta - tmpMat
              endIf
            endIf
          endIf
        elseIf(.not.doHerm.and.doAhrm) then
          if(mod(flag,100)/10.eq.0) then
            call MQC_Matrix_SymmMatrix_Put(alpha,aaARfoc,'antisymmetric')
            if(mod(flag,10).gt.0) then
            call MQC_Matrix_SymmMatrix_Put(beta,bbARfoc,'antisymmetric')
              if(mod(flag,10).eq.2) then
                call MQC_Matrix_SymmMatrix_Put(alphaBeta,abHARfoc,'symmetric')
                call MQC_Matrix_SymmMatrix_Put(tmpMat,abAARfoc,'antisymmetric')
                betaAlpha = alphaBeta - tmpMat
                alphaBeta = alphaBeta + tmpMat
              endIf
            endIf
          elseIf(mod(flag,100)/10.eq.1) then
            call MQC_Matrix_SymmMatrix_Put(alpha,cmplx(aaARfoc,aaACfoc),'antihermitian') 
            if(mod(flag,10).gt.0) then
              call MQC_Matrix_SymmMatrix_Put(beta,cmplx(bbARfoc,bbACfoc),'antihermitian') 
              if(mod(flag,10).eq.2) then
                call MQC_Matrix_SymmMatrix_Put(alphaBeta,cmplx(abHARfoc,abHACfoc),'hermitian')
                call MQC_Matrix_SymmMatrix_Put(tmpMat,cmplx(abAARfoc,abAACfoc),'antihermitian')
                betaAlpha = tmpMat - alphaBeta 
                alphaBeta = alphaBeta + tmpMat
              endIf
            endIf
          endIf
        endIf
        
      elseIf(algorithm.eq.'square') then

        if(eris(1)%integralType.ne.'regular') call mqc_error_A('2ERI contraction only &
          & implemented with regular integrals in mqc_eri_integral_contraction',6, &
          'eris(1)%integralType:',eris(1)%integralType)

!        !$OMP PARALLEL REDUCTION(+:alpha,beta,alphaBeta,betaAlpha) DEFAULT(NONE), &
        !$OMP PARALLEL DEFAULT(NONE), &
        !$OMP SHARED(eris,integral,zero,nBasis,mylabel,alpha,beta,alphaBeta,betaAlpha), &
        !$OMP PRIVATE(k,l,m,n,alphaLoc,betaLoc,alphaBetaLoc,betaAlphaLoc)
        call alphaLoc%init(nBasis,nBasis)
        if(integral%array_type.eq.'spin'.or.integral%array_type.eq.'general') &
          call betaLoc%init(nBasis,nBasis)
        if(integral%array_type.eq.'general') then
          call alphaBetaLoc%init(nBasis,nBasis)
          call betaAlphaLoc%init(nBasis,nBasis)
        endIf
        !$OMP DO COLLAPSE(4) 
        do k = 1, nBasis
          do l = 1, nBasis
            do m = 1, nBasis
              do n = 1, nBasis
                if(abs(eris(1)%at(k,l,m,n)).le.zero.and.abs(eris(1)%at(k,n,m,l)).le.zero) cycle
                ! AA Block
                if(myLabel.eq.'coulomb'.or.myLabel.eq.'doublebar') then
                  call alphaLoc%put(alphaLoc%at(k,l)+eris(1)%at(k,l,m,n)*integral%at(m,n),k,l)
                  call alphaLoc%put(alphaLoc%at(k,l)+eris(1)%at(k,l,m,n)*integral%at(m+nBasis,n+nBasis),k,l)
                endIf
                if(myLabel.eq.'exchange'.or.myLabel.eq.'doublebar') then
                  call alphaLoc%put(alphaLoc%at(k,l)-eris(1)%at(k,n,m,l)*integral%at(m,n),k,l)
                endIf
                !BB Block
                if(integral%array_type.eq.'spin'.or.integral%array_type.eq.'general') then
                  if(myLabel.eq.'coulomb'.or.myLabel.eq.'doublebar') then
                    call betaLoc%put(betaLoc%at(k,l)+eris(1)%at(k,l,m,n)*integral%at(m+nBasis,n+nBasis),k,l)
                    call betaLoc%put(betaLoc%at(k,l)+eris(1)%at(k,l,m,n)*integral%at(m,n),k,l)
                  endIf
                  if(myLabel.eq.'exchange'.or.myLabel.eq.'doublebar') then
                    call betaLoc%put(betaLoc%at(k,l)-eris(1)%at(k,n,m,l)*integral%at(m+nBasis,n+nBasis),k,l)
                  endIf
                endIf
                !AB Block
                if(integral%array_type.eq.'general') then
                  if(myLabel.eq.'exchange'.or.myLabel.eq.'doublebar') then
                    call alphaBetaLoc%put(alphaBetaLoc%at(k,l)-eris(1)%at(k,n,m,l)*integral%at(m+nBasis,n),k,l)
                  endIf
                !BA Block
                  if(myLabel.eq.'exchange'.or.myLabel.eq.'doublebar') then
                    call betaAlphaLoc%put(betaAlphaLoc%at(k,l)-eris(1)%at(k,n,m,l)*integral%at(m,n+nBasis),k,l)
                  endIf
                endIf
              end do
            end do
          end do
        end do
        !$OMP END DO
        !$OMP CRITICAL
        alpha = alpha + alphaLoc
        if(integral%array_type.eq.'spin'.or.integral%array_type.eq.'general') &
          beta = beta + betaLoc
        if(integral%array_type.eq.'general') then
          alphaBeta = alphaBeta + alphaBetaLoc
          betaAlpha = betaAlpha + betaAlphaLoc
        endIf
        !$OMP END CRITICAL
        !$OMP END PARALLEL
      endIf
!
      if(integral%array_type.eq.'space') call mqc_integral_allocate(output,'contraction','space',alpha)
      if(integral%array_type.eq.'spin') call mqc_integral_allocate(output,'contraction','spin',alpha,beta)
      if(integral%array_type.eq.'general') call mqc_integral_allocate(output,'contraction','general',alpha,beta,alphaBeta,betaAlpha)
!
      if(DEBUG) call output%print(6,'Final Fock matrix')
      
      end function mqc_eri_integral_contraction
!
!
!     PROCEDURE MQC_ERI_R4Tensor_Contraction
!
!>    \brief <b> MQC_ERI_R4Tensor_Contraction is used to return the contraction of 
!>    two-ERIs with an MQC rank-4 tensor</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_ERI_R4Tensor_Contraction is used to return the contraction of two-ERIs with 
!>    an MQC rank-4 tensor. This routine takes care of whether the rank-4 tensor includes
!>    separate spin components. The following options are available: 
!>
!>    1. Label = 'coulomb' returns the contraction with single bar coulomb integrals.
!>    2. Label = 'exchange' returns the contraction with single bar exchange integrals.
!>    3. Label = 'doublebar' returns the contraction with double-bar integrals (default).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] ERIs
!>    \verbatim
!>        ERIs is type(mqc_twoERIs) 
!>        The 2ERIs to be contracted.
!>    \endverbatim
!>
!>    \param[in] R4Tensor
!>    \verbatim
!>        R4Tensor is type(mqc_R4Tensor) 
!>        The MQC rank-4 tensor to be contracted.
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
!>    \author L. M. Thompson
!>    \date 2019
!
      function mqc_eri_r4tensor_contraction(eris,r4tensor,label) result(output)
!
      implicit none
      type(mqc_twoERIs),intent(in)::eris
      type(mqc_r4Tensor),intent(in)::r4tensor
      character(len=*),optional,intent(in)::label
      type(mqc_scalar)::output
      character(len=64)::mylabel
      integer(kind=int64)::r4t1,r4t2,r4t3,r4t4,nB1,nB2,nB3,nB4
      type(mqc_r4tensor)::temp
!      integer(kind=int64)::k,l,m,n,nBasis
!      character(len=64)::mylabel
!      real(kind=real64)::zero=1.0d-10
!     
      if(present(label)) then
        call string_change_case(label,'L',myLabel)
      else
        myLabel = 'doublebar'
      endIf
!
      r4t1 = size(r4tensor,1)
      r4t2 = size(r4tensor,2)
      r4t3 = size(r4tensor,3)
      r4t4 = size(r4tensor,4)
      if(eris%type().eq.'regular') then
        nB1 = eris%blockSize('alpha',1)
        nB2 = eris%blockSize('alpha',2)
        nB3 = eris%blockSize('alpha',3)
        nB4 = eris%blockSize('alpha',4)
        if(r4t1.eq.nB1.and.r4t2.eq.nB2.and.r4t3.eq.nB3.and.r4t4.eq.nB4) then
          select case (myLabel)
          case ('coulomb')
            temp = eris
            output = contraction(temp,r4tensor)
          case ('exchange')
            temp = mqc_eri_interaction(eris,'exchange') 
            output = contraction(temp,r4tensor)
          case ('doublebar')
            temp = mqc_eri_interaction(eris,'doublebar') 
            output = contraction(temp,r4tensor)
          case default
            call mqc_error_a('Unrecognized 2ERI integral interaction specified',6,'myLabel',myLabel)
          end select
        elseIf(r4t1.eq.2*nB1.and.r4t2.eq.2*nB2.and.r4t3.eq.2*nB3.and.r4t4.eq.2*nB4) then
          select case (myLabel)
          case ('coulomb')
            temp = eris 
            output = contraction(temp,r4tensor%tat([1,nB1],[1,nB2],[1,nB3],[1,nB4]))
            output = output + contraction(temp,r4tensor%tat([nB1+1,nB1*2],[nB2+1,nB2*2],[nB3+1,nB3*2],[nB4+1,nB4*2]))
            output = output + contraction(temp,r4tensor%tat([nB1+1,nB1*2],[nB2+1,nB2*2],[1,nB3],[1,nB4]))
            output = output + contraction(temp,r4tensor%tat([1,nB1],[1,nB2],[nB3+1,nB3*2],[nB4+1,nB4*2]))
          case ('exchange')
            temp = mqc_eri_interaction(eris,'exchange') 
            output = contraction(temp,r4tensor%tat([1,nB1],[1,nB2],[1,nB3],[1,nB4]))
            output = output + contraction(temp,r4tensor%tat([nB1+1,nB1*2],[nB2+1,nB2*2],[nB3+1,nB3*2],[nB4+1,nB4*2]))
            output = output + contraction(temp,r4tensor%tat([nB1+1,nB1*2],[1,nB2],[1,nB3],[nB4+1,nB4*2]))
            output = output + contraction(temp,r4tensor%tat([1,nB1],[nB2+1,nB2*2],[nB3+1,nB3*2],[1,nB4]))
          case ('doublebar')
            temp = mqc_eri_interaction(eris,'doublebar') 
            output = contraction(temp,r4tensor%tat([1,nB1],[1,nB2],[1,nB3],[1,nB4]))
            output = output + contraction(temp,r4tensor%tat([nB1+1,nB1*2],[nB2+1,nB2*2],[nB3+1,nB3*2],[nB4+1,nB4*2]))
            temp = mqc_eri_interaction(eris,'coulomb') 
            output = output + contraction(temp,r4tensor%tat([nB1+1,nB1*2],[nB2+1,nB2*2],[1,nB3],[1,nB4]))
            output = output + contraction(temp,r4tensor%tat([1,nB1],[1,nB2],[nB3+1,nB3*2],[nB4+1,nB4*2]))
            temp = mqc_eri_interaction(eris,'exchange') 
            output = output - contraction(temp,r4tensor%tat([nB1+1,nB1*2],[1,nB2],[1,nB3],[nB4+1,nB4*2]))
            output = output - contraction(temp,r4tensor%tat([1,nB1],[nB2+1,nB2*2],[nB3+1,nB3*2],[1,nB4]))
          case default
            call mqc_error_a('Unrecognized 2ERI integral interaction specified',6,'myLabel',myLabel)
          end select
        else
          call mqc_error_i('2ERIs and rank-4 tensor not correctly dimensioned in &
            & mqc_eri_r4tensor_contraction',6,'r4t1-nB1',r4t1-nB1,'r4t2-nB2',r4t2-nB2,'r4t3-nB3',r4t3-nB3,'r4t4-nB4', &
            r4t4-nB4)
        endIf
      else
        call mqc_error_a('mqc_eri_r4tensor_contraction only implemented with regular 2ERIs', &
          6,'eris%type()',eris%type())
      endIf
!
      end function mqc_eri_r4tensor_contraction
!
!
!     PROCEDURE MQC_ERI_Interaction
!
!>    \brief <b> MQC_ERI_Interaction is used to return the two-electron integrals of
!>    the requested two-body interaction type</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_ERI_Interaction is used to return the two-electron integrals of the 
!>    requested two-body interaction type. The following options are available: 
!>
!>    1. Label = 'coulomb' returns the 2ERIs with single bar coulomb integrals (default).
!>    2. Label = 'exchange' returns the 2ERIs with single bar exchange integrals.
!>    3. Label = 'doublebar' returns the 2ERIs with double-bar integrals.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] ERIs
!>    \verbatim
!>        ERIs is type(mqc_twoERIs) 
!>        The 2ERIs from which to construct the interaction type.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is character(len=*),optional 
!>        = 'coulomb':   returns coulomb integrals.
!>        = 'exchange':  returns exchange integrals.
!>        = 'doublebar': returns double bar integrals (default).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2019
!
      function mqc_eri_interaction(eris,label) result(output)
!
      implicit none
      type(mqc_twoERIs),intent(in)::eris
      character(len=*),optional,intent(in)::label
      type(mqc_twoERIs)::output
      character(len=64)::mylabel
      integer(kind=int64)::nBasis,i,j
      type(mqc_r4tensor)::temp,temp2
!     
      if(present(label)) then
        call string_change_case(label,'L',myLabel)
      else
        myLabel = 'coulomb'
      endIf

      select case (eris%type())
      case ('regular')
        select case (myLabel)
        case('coulomb')
          output = eris
        case('exchange')
          temp = eris
          temp2 = eris
          nBasis = eris%blockSize('alpha')
          do i = 1, nBasis
            do j = 1, nBasis
              call temp%tput(temp2%tat([0],[i,i],[0],[j,j]),[0],[j,j],[0],[i,i])
            endDo
          endDo
          call mqc_twoeris_allocate(output,'full','regular',temp)
        case('doublebar')
          temp = eris
          temp2 = eris
          nBasis = eris%blockSize('alpha')
          do i = 1, nBasis
            do j = 1, nBasis
              call temp%tput(temp2%tat([0],[i,i],[0],[j,j]),[0],[j,j],[0],[i,i])
            endDo
          endDo
          temp = temp2 - temp
          call mqc_twoeris_allocate(output,'full','regular',temp)
        case default
          call mqc_error_a('Interaction type not recognized in mqc_eri_interaction',&
            6,'myLabel',myLabel) 
        end select
      case ('space')
        select case (myLabel)
        case('coulomb')
          output = eris
        case('exchange')
          temp = eris
          temp2 = eris
          nBasis = eris%blockSize('alpha')
          do i = 1, nBasis
            do j = 1, nBasis
              call temp%tput(temp2%tat([0],[i,i],[0],[j,j]),[0],[j,j],[0],[i,i])
            endDo
          endDo
          call mqc_twoeris_allocate(output,'full','space',temp)
        case('doublebar')
          temp = eris
          temp2 = eris
          nBasis = eris%blockSize('alpha')
          do i = 1, nBasis
            do j = 1, nBasis
              call temp%tput(temp2%tat([0],[i,i],[0],[j,j]),[0],[j,j],[0],[i,i])
            endDo
          endDo
          temp = temp2 - temp
          call mqc_twoeris_allocate(output,'full','space',temp)
        case default
          call mqc_error_a('Interaction type not recognized in mqc_eri_interaction',&
            6,'myLabel',myLabel) 
        end select
      case ('spin')
        select case (myLabel)
        case('coulomb')
          output = eris
        case('exchange')
!
!         Exchange spin blocks will require general 2ERI form
!
!         (aa|aa) -> (aa|aa)
!         (bb|bb) -> (bb|bb)
!         (aa|bb) -> 0
!         (bb|aa) -> 0       = (aa|bb) with symm
!         (ab|ba) -> (aa|bb)
!         (ba|ab) -> (bb|aa) = (ab|ba) with symm
!
!         Note that the result will have many zeros and so there are
!         probably better alternatives to this routine
!
          nBasis = eris%blockSize('alpha')
          call temp%init(2*nBasis,2*nBasis,2*nBasis,2*nBasis)
          temp2 = eris
          do i = 1, nBasis
            do j = 1, nBasis
              call temp%tput(temp2%tat([1,nBasis],[i,i],[1,nBasis],[j,j]),[1,nBasis],[j,j],[1,nBasis],[i,i])
              call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],[nBasis+1,2*nBasis],&
                [j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],&
                [i+nBasis,i+nBasis])
              call temp%tput(temp2%tat([1,nBasis],[i,i],[nBasis+1,2*nBasis],[j+nBasis,j+nBasis]),&
                [1,nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],[i,i])
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],&
                [1,nBasis],[j,j]),[nBasis+1,2*nBasis],[j,j],[1,nBasis],[i+nBasis,i+nBasis])
            endDo
          endDo
          if(eris%storageType.eq.'full') then
            call mqc_twoeris_allocate(output,'full','general',&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]))
          else
            call mqc_twoeris_allocate(output,'symm','general',&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              abab=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              abba=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              aaab=temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              abbb=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]))
          endIf
        case('doublebar')
!
!         Doublebar spin blocks will require general 2ERI form
!
!         (aa||aa) -> (aa|aa)-(aa|aa)
!         (bb||bb) -> (bb|bb)-(bb|bb) 
!         (aa||bb) -> (aa|bb) - 0
!         (bb||aa) -> (bb|aa) - 0     = (aa||bb) with symm
!         (ab||ba) -> 0 - (aa|bb)
!         (ba||ab) -> 0 - (bb|aa)     = (ab||ba) with symm
!
!         Note that the result will have many zeros and so there are
!         probably better alternatives to this routine
!
          nBasis = eris%blockSize('alpha')
          call temp%init(2*nBasis,2*nBasis,2*nBasis,2*nBasis)
          temp2 = eris
          do i = 1, nBasis
            do j = 1, nBasis
              call temp%tput(temp2%tat([1,nBasis],[i,i],[1,nBasis],[j,j]),[1,nBasis],[j,j],[1,nBasis],[i,i])
              call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],[nBasis+1,2*nBasis],&
                [j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],&
                [i+nBasis,i+nBasis])
              call temp%tput(temp2%tat([1,nBasis],[i,i],[nBasis+1,2*nBasis],[j+nBasis,j+nBasis]),&
                [1,nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],[i,i])
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],&
                [1,nBasis],[j,j]),[nBasis+1,2*nBasis],[j,j],[1,nBasis],[i+nBasis,i+nBasis])
            endDo
          endDo
          temp = temp2 - temp
          if(eris%storageType.eq.'full') then
            call mqc_twoeris_allocate(output,'full','general',&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]))
          else
            call mqc_twoeris_allocate(output,'symm','general',&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              abab=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              abba=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              aaab=temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              abbb=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]))
          endIf
        case default
          call mqc_error_a('Interaction type not recognized in mqc_eri_interaction',&
            6,'myLabel',myLabel) 
        end select
      case ('general')
        select case (myLabel)
        case('coulomb')
          output = eris
        case('exchange')
          nBasis = eris%blockSize('alpha')
          call temp%init(2*nBasis,2*nBasis,2*nBasis,2*nBasis)
          temp2 = eris
          do i = 1, nBasis
            do j = 1, nBasis
              !aaaa
              call temp%tput(temp2%tat([1,nBasis],[i,i],[1,nBasis],[j,j]),[1,nBasis],[j,j],[1,nBasis],[i,i])
              !bbbb
              call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],[nBasis+1,2*nBasis],&
                [j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],&
                [i+nBasis,i+nBasis])
              !aabb
              call temp%tput(temp2%tat([1,nBasis],[i+nBasis,i+nBasis],[nBasis+1,2*nBasis],[j,j]),&
                [1,nBasis],[j,j],[nBasis+1,2*nBasis],[i+nBasis,i+nBasis])
              !bbaa
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i,i],&
                [1,nBasis],[j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],[1,nBasis],[i,i])
              !abba
              call temp%tput(temp2%tat([1,nBasis],[i,i],[nBasis+1,2*nBasis],[j+nBasis,j+nBasis]),&
                [1,nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],[i,i])
              !baab
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],&
                [1,nBasis],[j,j]),[nBasis+1,2*nBasis],[j,j],[1,nBasis],[i+nBasis,i+nBasis])
              !abab
              call temp%tput(temp2%tat([1,nBasis],[i+nBasis,i+nBasis],[1,nBasis],[j+nBasis,j+nBasis]),&
                [1,nBasis],[j+nBasis,j+nBasis],[1,nBasis],[i+nBasis,i+nBasis])
              !baba
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i,i],&
                [nBasis+1,2*nBasis],[j,j]),[nBasis+1,2*nBasis],[j,j],[nBasis+1,2*nBasis],[i,i])
              !aaab
              call temp%tput(temp2%tat([1,nBasis],[i+nBasis,i+nBasis],[1,nBasis],[j,j]),[1,nBasis],[j,j],&
                [1,nBasis],[i+nBasis,i+nBasis])
              !aaba
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([1,nBasis],[i,i],[nBasis+1,2*nBasis],&
                [j,j]),[1,nBasis],[j,j],[nBasis+1,2*nBasis],[i,i])
              !abaa
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([1,nBasis],[i,i],[1,nBasis],&
                [j+nBasis,j+nBasis]),[1,nBasis],[j+nBasis,j+nBasis],[1,nBasis],[i,i])
              !baaa
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i,i],[1,nBasis],&
                [j,j]),[nBasis+1,2*nBasis],[j,j],[1,nBasis],[i,i])
              !bbba
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i,i],&
                [nBasis+1,2*nBasis],[j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],&
                [nBasis+1,2*nBasis],[i,i])
              !bbab
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],&
                [1,nBasis],[j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],[1,nBasis],&
                [i+nBasis,i+nBasis])
              !babb
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],&
                [nBasis+1,2*nBasis],[j,j]),[nBasis+1,2*nBasis],[j,j],[nBasis+1,2*nBasis],&
                [i+nBasis,i+nBasis])
              !abbb
              call temp%tput(temp2%tat([1,nBasis],[i+nBasis,i+nBasis],[nBasis+1,2*nBasis],[j+nBasis,j+nBasis]),&
                [1,nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],[i+nBasis,i+nBasis])
            endDo
          endDo
          if(eris%storageType.eq.'full') then
            call mqc_twoeris_allocate(output,'full','general',&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]))
          else
            call mqc_twoeris_allocate(output,'symm','general',&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              abab=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              abba=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              aaab=temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              abbb=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]))
          endIf
        case('doublebar')
          nBasis = eris%blockSize('alpha')
          call temp%init(2*nBasis,2*nBasis,2*nBasis,2*nBasis)
          temp2 = eris
          do i = 1, nBasis
            do j = 1, nBasis
              !aaaa
              call temp%tput(temp2%tat([1,nBasis],[i,i],[1,nBasis],[j,j]),[1,nBasis],[j,j],[1,nBasis],[i,i])
              !bbbb
              call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],[nBasis+1,2*nBasis],&
                [j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],&
                [i+nBasis,i+nBasis])
              !aabb
              call temp%tput(temp2%tat([1,nBasis],[i+nBasis,i+nBasis],[nBasis+1,2*nBasis],[j,j]),&
                [1,nBasis],[j,j],[nBasis+1,2*nBasis],[i+nBasis,i+nBasis])
              !bbaa
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i,i],&
                [1,nBasis],[j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],[1,nBasis],[i,i])
              !abba
              call temp%tput(temp2%tat([1,nBasis],[i,i],[nBasis+1,2*nBasis],[j+nBasis,j+nBasis]),&
                [1,nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],[i,i])
              !baab
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],&
                [1,nBasis],[j,j]),[nBasis+1,2*nBasis],[j,j],[1,nBasis],[i+nBasis,i+nBasis])
              !abab
              call temp%tput(temp2%tat([1,nBasis],[i+nBasis,i+nBasis],[1,nBasis],[j+nBasis,j+nBasis]),&
                [1,nBasis],[j+nBasis,j+nBasis],[1,nBasis],[i+nBasis,i+nBasis])
              !baba
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i,i],&
                [nBasis+1,2*nBasis],[j,j]),[nBasis+1,2*nBasis],[j,j],[nBasis+1,2*nBasis],[i,i])
              !aaab
              call temp%tput(temp2%tat([1,nBasis],[i+nBasis,i+nBasis],[1,nBasis],[j,j]),[1,nBasis],[j,j],&
                [1,nBasis],[i+nBasis,i+nBasis])
              !aaba
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([1,nBasis],[i,i],[nBasis+1,2*nBasis],&
                [j,j]),[1,nBasis],[j,j],[nBasis+1,2*nBasis],[i,i])
              !abaa
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([1,nBasis],[i,i],[1,nBasis],&
                [j+nBasis,j+nBasis]),[1,nBasis],[j+nBasis,j+nBasis],[1,nBasis],[i,i])
              !baaa
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i,i],[1,nBasis],&
                [j,j]),[nBasis+1,2*nBasis],[j,j],[1,nBasis],[i,i])
              !bbba
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i,i],&
                [nBasis+1,2*nBasis],[j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],&
                [nBasis+1,2*nBasis],[i,i])
              !bbab
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],&
                [1,nBasis],[j+nBasis,j+nBasis]),[nBasis+1,2*nBasis],[j+nBasis,j+nBasis],[1,nBasis],&
                [i+nBasis,i+nBasis])
              !babb
              if(eris%storageType.eq.'full') call temp%tput(temp2%tat([nBasis+1,2*nBasis],[i+nBasis,i+nBasis],&
                [nBasis+1,2*nBasis],[j,j]),[nBasis+1,2*nBasis],[j,j],[nBasis+1,2*nBasis],&
                [i+nBasis,i+nBasis])
              !abbb
              call temp%tput(temp2%tat([1,nBasis],[i+nBasis,i+nBasis],[nBasis+1,2*nBasis],[j+nBasis,j+nBasis]),&
                [1,nBasis],[j+nBasis,j+nBasis],[nBasis+1,2*nBasis],[i+nBasis,i+nBasis])
            endDo
          endDo
          temp = temp2 - temp
          if(eris%storageType.eq.'full') then
            call mqc_twoeris_allocate(output,'full','general',&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]))
          else
            call mqc_twoeris_allocate(output,'symm','general',&
              temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[1,nBasis]),&
              temp%tat([nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              temp%tat([1,nBasis],[1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]),&
              abab=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              abba=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[1,nBasis]),&
              aaab=temp%tat([1,nBasis],[1,nBasis],[1,nBasis],[nBasis+1,2*nBasis]),&
              abbb=temp%tat([1,nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis],[nBasis+1,2*nBasis]))
          endIf
        case default
          call mqc_error_a('Interaction type not recognized in mqc_eri_interaction',&
            6,'myLabel',myLabel) 
        end select
      case default
        call mqc_error_a('2ERI type either not recognized or not yet implemented in &
          & mqc_eri_interaction',6,'eris%type()',eris%type()) 
      end select

      end function mqc_eri_interaction
!
!     PROCEDURE MQC_SCF_Integral_Generalized_Eigensystem
!
!>    \brief <b> MQC_Integral_Generalized_Eigensystem is a subroutine that takes a 
!>    MQC integral and a metric and returns eigenvalues and eigenvectors</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Generalized_Eigensystem is a subroutine that takes a MQC integral 
!>    and an optional metric (also as an MQC integral) and optionally returns 
!>    eigenvalues to a MQC eigenvalues variable and/or left and/or right eigenvectors 
!>    to an MQC integral variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        IntegralA is Class(MQC_SCF_Integral)
!>        The MQC integral to diagonalize. 
!>    \endverbatim
!>
!>    \param[in] IntegralB
!>    \verbatim
!>        IntegralB is Type(MQC_SCF_Integral),Optional
!>        The metric. 
!>    \endverbatim
!>
!>    \param[out] EVals
!>    \verbatim
!>        EVals is Type(MQC_SCF_Eigenvalues),Optional
!>        Optional MQC eigenvalues variable containing the 
!>        eigenvalues.
!>    \endverbatim
!>
!>    \param[out] REVecs
!>    \verbatim
!>        REVecs is Type(MQC_SCF_Integral),Optional
!>        Optional MQC integral containing the right eigenvectors.
!>    \endverbatim
!>
!>    \param[out] LEVecs
!>    \verbatim
!>        LEVecs is Type(MQC_SCF_Integral),Optional
!>        Optional MQC integral containing the left eigenvectors.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      subroutine mqc_scf_integral_generalized_eigensystem(integralA,integralB,eVals,rEVecs,lEVecs)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integralA
      type(mqc_scf_integral),optional::integralB
      type(mqc_scf_eigenvalues),optional,intent(inOut)::eVals
      type(mqc_scf_integral),optional,intent(inOut)::rEVecs,lEVecs
      type(mqc_matrix)::tmpMat1,tmpMat2,tmpMat3,tmpMat4,tmpMat5,tmpMat6,tmpMat7,tmpMat8
      type(mqc_vector)::tmpVec1,tmpVec2
      integer(kind=int64)::nDimAlpha1=0,nDimAlpha2=0
!
      if(integralA%type().eq.'space') then
        if(present(integralB)) then
          if(integralB%type().eq.'space') then
            tmpMat1 = integralA%getblock('alpha')
            tmpMat2 = integralB%getblock('alpha')
            call tmpMat1%eigensys(tmpMat2,tmpVec1,tmpMat3,tmpMat4)
            If(present(eVals)) call mqc_eigenvalues_allocate(eVals,'eigenvalues','space',tmpVec1) 
            If(present(rEVecs)) call mqc_integral_allocate(rEVecs,'eigenvectors','space',tmpMat3) 
            If(present(lEVecs)) call mqc_integral_allocate(lEVecs,'eigenvectors','space',tmpMat4)
          elseIf(integralB%type().eq.'spin') then
            tmpMat1 = integralA%getblock('alpha')
            tmpMat2 = integralB%getblock('alpha')
            tmpMat3 = integralB%getblock('beta')
            call tmpMat1%eigensys(tmpMat2,tmpVec1,tmpMat4,tmpMat5)
            call tmpMat1%eigensys(tmpMat3,tmpVec2,tmpMat6,tmpMat7)
            If(present(eVals)) call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2) 
            If(present(rEVecs)) call mqc_integral_allocate(rEVecs,'eigenvectors','spin',tmpMat4,tmpMat6) 
            If(present(lEVecs)) call mqc_integral_allocate(lEVecs,'eigenvectors','spin',tmpMat5,tmpMat7)
          elseIf(integralB%type().eq.'general') then
            nDimAlpha1 = integralA%blockSize('alpha',1) 
            nDimAlpha2 = integralA%blockSize('alpha',2) 
            if (integralA%blockSize('alpha',1).ne.integralB%blockSize('alpha',1)) call mqc_error('Alpha block 1 sizes &
      &       must match in mqc_scf_integral_generalized_eigensystem') 
            if (integralB%blockSize('beta',1).ne.integralB%blockSize('beta',1)) call mqc_error('Beta block 1 sizes &
      &       must match in mqc_scf_integral_generalized_eigensystem') 
            if (integralA%blockSize('alpha',2).ne.integralB%blockSize('alpha',2)) call mqc_error('Alpha block 2 sizes &
      &       must match in mqc_scf_integral_generalized_eigensystem') 
            if (integralB%blockSize('beta',2).ne.integralB%blockSize('beta',2)) call mqc_error('Beta block 2 sizes &
      &       must match in mqc_scf_integral_generalized_eigensystem') 
            tmpMat1 = integralA%getblock('full')
            tmpMat2 = integralB%getblock('full')
            call tmpMat1%eigensys(tmpMat2,tmpVec1,tmpMat3,tmpMat4)
            if(present(rEVecs)) then
              tmpMat5 = tmpMat3%mat([1,NDimAlpha1],[1,NDimAlpha2])
              tmpMat6 = tmpMat3%mat([NDimAlpha1+1,-1],[NDimAlpha2+1,-1])
              tmpMat7 = tmpMat3%mat([NDimAlpha1+1,-1],[1,NDimAlpha2])
              tmpMat8 = tmpMat3%mat([1,NDimAlpha1],[NDimAlpha2+1,-1])
              call mqc_integral_allocate(rEVecs,'eigenvectors','general',tmpMat5,tmpMat6,tmpMat7,tmpMat8)
            endIf
            if(present(lEVecs)) then
              tmpMat5 = tmpMat4%mat([1,NDimAlpha1],[1,NDimAlpha2])
              tmpMat6 = tmpMat4%mat([NDimAlpha1+1,-1],[NDimAlpha2+1,-1])
              tmpMat7 = tmpMat4%mat([NDimAlpha1+1,-1],[1,NDimAlpha2])
              tmpMat8 = tmpMat4%mat([1,NDimAlpha1],[NDimAlpha2+1,-1])
              call mqc_integral_allocate(lEVecs,'eigenvectors','general',tmpMat5,tmpMat6,tmpMat7,tmpMat8)
            endIf
            if(present(eVals)) then
              tmpVec2 = tmpVec1%vat(NDimAlpha2+1,-1)
              tmpVec1 = tmpVec1%vat(1,NDimAlpha2)
              call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2)
            endIf
          endIf
        else
          tmpMat1 = integralA%getblock('alpha')
          call tmpMat1%eigensys(eigenvals=tmpVec1,reigenvecs=tmpMat3,leigenvecs=tmpMat4)
          If(present(eVals)) call mqc_eigenvalues_allocate(eVals,'eigenvalues','space',tmpVec1) 
          If(present(rEVecs)) call mqc_integral_allocate(rEVecs,'eigenvectors','space',tmpMat3) 
          If(present(lEVecs)) call mqc_integral_allocate(lEVecs,'eigenvectors','space',tmpMat4)
        endIf
      elseIf(integralA%type().eq.'spin') then
        if(present(integralB)) then
          if(integralB%type().eq.'space') then
            tmpMat1 = integralA%getblock('alpha')
            tmpMat2 = integralA%getblock('beta')
            tmpMat3 = integralB%getblock('alpha')
            call tmpMat1%eigensys(tmpMat3,tmpVec1,tmpMat4,tmpMat5)
            call tmpMat2%eigensys(tmpMat3,tmpVec2,tmpMat6,tmpMat7)
            If(present(eVals)) call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2) 
            If(present(rEVecs)) call mqc_integral_allocate(rEVecs,'eigenvectors','spin',tmpMat4,tmpMat6) 
            If(present(lEVecs)) call mqc_integral_allocate(lEVecs,'eigenvectors','spin',tmpMat5,tmpMat7)
          elseIf(integralB%type().eq.'spin') then
            tmpMat1 = integralA%getblock('alpha')
            tmpMat2 = integralB%getblock('alpha')
            tmpMat3 = integralA%getblock('beta')
            tmpMat4 = integralB%getblock('beta')
            call tmpMat1%eigensys(tmpMat2,tmpVec1,tmpMat5,tmpMat6)
            call tmpMat3%eigensys(tmpMat4,tmpVec2,tmpMat7,tmpMat8)
            If(present(eVals)) call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2) 
            If(present(rEVecs)) call mqc_integral_allocate(rEVecs,'eigenvectors','spin',tmpMat5,tmpMat7) 
            If(present(lEVecs)) call mqc_integral_allocate(lEVecs,'eigenvectors','spin',tmpMat6,tmpMat8)
          elseIf(integralB%type().eq.'general') then
            nDimAlpha1 = integralA%blockSize('alpha',1) 
            nDimAlpha2 = integralA%blockSize('alpha',2) 
            if (integralA%blockSize('alpha',1).ne.integralB%blockSize('alpha',1)) call mqc_error('Alpha block 1 sizes &
      &       must match in mqc_scf_integral_generalized_eigensystem') 
            if (integralB%blockSize('beta',1).ne.integralB%blockSize('beta',1)) call mqc_error('Beta block 1 sizes &
      &       must match in mqc_scf_integral_generalized_eigensystem') 
            if (integralA%blockSize('alpha',2).ne.integralB%blockSize('alpha',2)) call mqc_error('Alpha block 2 sizes &
      &       must match in mqc_scf_integral_generalized_eigensystem') 
            if (integralB%blockSize('beta',2).ne.integralB%blockSize('beta',2)) call mqc_error('Beta block 2 sizes &
      &       must match in mqc_scf_integral_generalized_eigensystem') 
            tmpMat1 = integralA%getblock('full')
            tmpMat2 = integralB%getblock('full')
            call tmpMat1%eigensys(tmpMat2,tmpVec1,tmpMat3,tmpMat4)
            if(present(rEVecs)) then
              tmpMat5 = tmpMat3%mat([1,NDimAlpha1],[1,NDimAlpha2])
              tmpMat6 = tmpMat3%mat([NDimAlpha1+1,-1],[NDimAlpha2+1,-1])
              tmpMat7 = tmpMat3%mat([NDimAlpha1+1,-1],[1,NDimAlpha2])
              tmpMat8 = tmpMat3%mat([1,NDimAlpha1],[NDimAlpha2+1,-1])
              call mqc_integral_allocate(rEVecs,'eigenvectors','general',tmpMat5,tmpMat6,tmpMat7,tmpMat8)
            endIf
            if(present(lEVecs)) then
              tmpMat5 = tmpMat4%mat([1,NDimAlpha1],[1,NDimAlpha2])
              tmpMat6 = tmpMat4%mat([NDimAlpha1+1,-1],[NDimAlpha2+1,-1])
              tmpMat7 = tmpMat4%mat([NDimAlpha1+1,-1],[1,NDimAlpha2])
              tmpMat8 = tmpMat4%mat([1,NDimAlpha1],[NDimAlpha2+1,-1])
              call mqc_integral_allocate(lEVecs,'eigenvectors','general',tmpMat5,tmpMat6,tmpMat7,tmpMat8)
            endIf
            if(present(eVals)) then
              tmpVec2 = tmpVec1%vat(NDimAlpha2+1,-1)
              tmpVec1 = tmpVec1%vat(1,NDimAlpha2)
              call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2)
            endIf
          endIf
        else
          tmpMat1 = integralA%getblock('alpha')
          tmpMat2 = integralA%getblock('beta')
          call tmpMat1%eigensys(eigenvals=tmpVec1,reigenvecs=tmpMat3,leigenvecs=tmpMat4)
          call tmpMat2%eigensys(eigenvals=tmpVec2,reigenvecs=tmpMat5,leigenvecs=tmpMat6)
          If(present(eVals)) call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2) 
          If(present(rEVecs)) call mqc_integral_allocate(rEVecs,'eigenvectors','space',tmpMat3,tmpMat4) 
          If(present(lEVecs)) call mqc_integral_allocate(lEVecs,'eigenvectors','space',tmpMat5,tmpMat6)
        endIf
      elseIf(integralA%type().eq.'general') then
        nDimAlpha1 = integralA%blockSize('alpha',1) 
        nDimAlpha2 = integralA%blockSize('alpha',2) 
        tmpMat1 = integralA%getblock('full')
        if(present(integralB)) then
          if (integralA%blockSize('alpha',1).ne.integralB%blockSize('alpha',1)) call mqc_error('Alpha block 1 sizes &
      &     must match in mqc_scf_integral_generalized_eigensystem') 
          if (integralB%blockSize('beta',1).ne.integralB%blockSize('beta',1)) call mqc_error('Beta block 1 sizes &
      &     must match in mqc_scf_integral_generalized_eigensystem') 
          if (integralA%blockSize('alpha',2).ne.integralB%blockSize('alpha',2)) call mqc_error('Alpha block 2 sizes &
      &     must match in mqc_scf_integral_generalized_eigensystem') 
          if (integralB%blockSize('beta',2).ne.integralB%blockSize('beta',2)) call mqc_error('Beta block 2 sizes &
      &     must match in mqc_scf_integral_generalized_eigensystem') 
          tmpMat2 = integralB%getblock('full')
          call tmpMat1%eigensys(tmpMat2,tmpVec1,tmpMat3,tmpMat4)
        else
          call tmpMat1%eigensys(eigenvals=tmpVec1,reigenvecs=tmpMat3,leigenvecs=tmpMat4)
        endIf
        if(present(rEVecs)) then
          tmpMat5 = tmpMat3%mat([1,NDimAlpha1],[1,NDimAlpha2])
          tmpMat6 = tmpMat3%mat([NDimAlpha1+1,-1],[NDimAlpha2+1,-1])
          tmpMat7 = tmpMat3%mat([NDimAlpha1+1,-1],[1,NDimAlpha2])
          tmpMat8 = tmpMat3%mat([1,NDimAlpha1],[NDimAlpha2+1,-1])
          call mqc_integral_allocate(rEVecs,'eigenvectors','general',tmpMat5,tmpMat6,tmpMat7,tmpMat8)
        endIf
        if(present(lEVecs)) then
          tmpMat5 = tmpMat4%mat([1,NDimAlpha1],[1,NDimAlpha2])
          tmpMat6 = tmpMat4%mat([NDimAlpha1+1,-1],[NDimAlpha2+1,-1])
          tmpMat7 = tmpMat4%mat([NDimAlpha1+1,-1],[1,NDimAlpha2])
          tmpMat8 = tmpMat4%mat([1,NDimAlpha1],[NDimAlpha2+1,-1])
          call mqc_integral_allocate(lEVecs,'eigenvectors','general',tmpMat5,tmpMat6,tmpMat7,tmpMat8)
        endIf
        if(present(eVals)) then
          tmpVec2 = tmpVec1%vat(NDimAlpha2+1,-1)
          tmpVec1 = tmpVec1%vat(1,NDimAlpha2)
          call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2)
        endIf
      endIf
!
      end subroutine mqc_scf_integral_generalized_eigensystem
!
!
!     PROCEDURE MQC_SCF_Integral_Diagonalize
!
!>    \brief <b> MQC_SCF_Integral_Diagonalize is a subroutine that takes a symmetric 
!>    or hermitian MQC integral and returns eigenvalues and eigenvectors</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Integral_Diagonalize is a subroutine that takes a symmetric or 
!>    hermitian MQC integral and optionally returns eigenvalues to a MQC 
!>    eigenvalues variable and/or eigenvectors to a MQC integral variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral to diagonalize. 
!>    \endverbatim
!>
!>    \param[in,out] EVals
!>    \verbatim
!>        EVals is Type(MQC_SCF_Eigenvalues),Optional
!>        Optional MQC eigenvalues variable containing the 
!>        eigenvalues.
!>    \endverbatim
!>
!>    \param[in,out] EVecs
!>    \verbatim
!>        EVecs is Type(MQC_SCF_Integral),Optional
!>        Optional MQC integral containing the eigenvectors.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      subroutine mqc_scf_integral_diagonalize(integral,eVals,eVecs)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      type(mqc_scf_eigenvalues),optional,intent(inOut)::eVals
      type(mqc_scf_integral),optional,intent(inOut)::eVecs
      type(mqc_matrix)::tmpMat1,tmpMat2,tmpMat3,tmpMat4
      type(mqc_vector)::tmpVec1,tmpVec2
      integer(kind=int64)::nDimAlpha1=0,nDimAlpha2=0
!
      if(integral%type().eq.'space') then
        tmpMat1 = integral%getblock('alpha')
        if(present(eVals).and.present(eVecs)) then
          call tmpMat1%diag(tmpVec1,tmpMat2)
          call mqc_integral_allocate(eVecs,'eigenvectors','space',tmpMat2)
          call mqc_eigenvalues_allocate(eVals,'eigenvalues','space',tmpVec1)
        elseIf(present(eVecs)) then
          call tmpMat1%diag(EVecs=tmpMat2)
          call mqc_integral_allocate(eVecs,'eigenvectors','space',tmpMat2)
        elseIf(present(eVals)) then
          call tmpMat1%diag(tmpVec1)
          call mqc_eigenvalues_allocate(eVals,'eigenvalues','space',tmpVec1)
        endIf
      elseIf(integral%type().eq.'spin') then
        tmpMat1 = integral%getblock('alpha')
        tmpMat2 = integral%getblock('beta')
        if(present(eVals).and.present(eVecs)) then
          call tmpMat1%diag(tmpVec1,tmpMat3)
          call tmpMat2%diag(tmpVec2,tmpMat4)
          call mqc_integral_allocate(eVecs,'eigenvectors','spin',tmpMat3,tmpMat4)
          call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2)
        elseIf(present(eVecs)) then
          call tmpMat1%diag(EVecs=tmpMat3)
          call tmpMat2%diag(EVecs=tmpMat4)
          call mqc_integral_allocate(eVecs,'eigenvectors','spin',tmpMat3,tmpMat4)
        elseIf(present(eVals)) then
          call tmpMat1%diag(tmpVec1)
          call tmpMat2%diag(tmpVec2)
          call mqc_eigenvalues_allocate(eVals,'eigenvalues','spin',tmpVec1,tmpVec2)
        endIf
      elseIf(integral%type().eq.'general') then
        nDimAlpha1 = integral%blockSize('alpha',1) 
        nDimAlpha2 = integral%blockSize('alpha',2) 
        tmpMat1 = integral%getblock('full')
        if(present(eVals).and.present(eVecs)) then
          call tmpMat1%diag(tmpVec1,tmpMat2)
          tmpMat1 = tmpMat2%mat([1,NDimAlpha1],[1,NDimAlpha2])
          tmpMat3 = tmpMat2%mat([NDimAlpha1+1,-1],[1,NDimAlpha2])
          tmpMat4 = tmpMat2%mat([1,NDimAlpha1],[NDimAlpha2+1,-1])
          tmpMat2 = tmpMat2%mat([NDimAlpha1+1,-1],[NDimAlpha2+1,-1])
          tmpVec2 = tmpVec1%vat(NDimAlpha2+1,-1)
          tmpVec1 = tmpVec1%vat(1,NDimAlpha2)
          call mqc_integral_allocate(eVecs,'eigenvectors','general',tmpMat1,tmpMat2,tmpMat3,tmpMat4)
          call mqc_eigenvalues_allocate(eVals,'eigenvalues','general',tmpVec1,tmpVec2)
        elseIf(present(eVecs)) then
          call tmpMat1%diag(EVecs=tmpMat2)
          tmpMat1 = tmpMat2%mat([1,NDimAlpha1],[1,NDimAlpha2])
          tmpMat3 = tmpMat2%mat([NDimAlpha1+1,-1],[1,NDimAlpha2])
          tmpMat4 = tmpMat2%mat([1,NDimAlpha1],[NDimAlpha2+1,-1])
          tmpMat2 = tmpMat2%mat([NDimAlpha1+1,-1],[NDimAlpha2+1,-1])
          call mqc_integral_allocate(eVecs,'eigenvectors','general',tmpMat1,tmpMat2,tmpMat3,tmpMat4)
        elseIf(present(eVals)) then
          call tmpMat1%diag(tmpVec1)
          tmpVec2 = tmpVec1%vat(NDimAlpha2+1,-1)
          tmpVec1 = tmpVec1%vat(1,NDimAlpha2)
          call mqc_eigenvalues_allocate(eVals,'eigenvalues','general',tmpVec1,tmpVec2)
        endIf
      endIf
!
      end subroutine mqc_scf_integral_diagonalize
!
!
!     PROCEDURE MQC_SCF_Integral_SVD
!
!>    \brief <b> MQC_SCF_Integral_SVD is a subroutine that performs singular value
!>    decomposition on a MQC integral and returns singular values and vectors</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Integral_SVD is a subroutine that performs singular value decomposition 
!>    on a MQC integral and returns singular values and vectors</b>
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral to perform SVD on. 
!>    \endverbatim
!>
!>    \param[in,out] EVals
!>    \verbatim
!>        EVals is Type(MQC_SCF_Eigenvalues),Optional
!>        Optional MQC eigenvalues variable containing the 
!>        singular values.
!>    \endverbatim
!>
!>    \param[in,out] EUVecs
!>    \verbatim
!>        EUVecs is Type(MQC_SCF_Integral),Optional
!>        Optional MQC integral containing the U vectors.
!>    \endverbatim
!>
!>    \param[in,out] EVVecs
!>    \verbatim
!>        EVVecs is Type(MQC_SCF_Integral),Optional
!>        Optional MQC integral containing the V vectors.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2021
!
      subroutine mqc_scf_integral_svd(integral,eVals,eUVecs,eVVecs)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      type(mqc_scf_eigenvalues),optional,intent(inOut)::eVals
      type(mqc_scf_integral),optional,intent(inOut)::eUVecs,eVVecs
      type(mqc_matrix)::tmpMat1,tmpMat2,tmpMat3,tmpMat4,tmpMat5,tmpMat6,tmpMat7,tmpMat8
      type(mqc_vector)::tmpVec1,tmpVec2
      integer(kind=int64)::nDimAlpha1=0,nDimAlpha2=0,nDimBeta1=0,nDimBeta2=0
!
      if(integral%type().eq.'space') then
        tmpMat1 = integral%getblock('alpha')
        if(present(eVals).and.present(eUVecs).and.present(eVVecs)) then
          call tmpMat1%svd(tmpVec1,tmpMat2,tmpMat3)
          call mqc_integral_allocate(eUVecs,'left singular vectors','space',tmpMat2)
          call mqc_integral_allocate(eVVecs,'right singular vectors','space',tmpMat3)
          call mqc_eigenvalues_allocate(eVals,'singular values','space',tmpVec1)
        elseIf(present(eVals).and.present(eUVecs)) then
          call tmpMat1%svd(tmpVec1,tmpMat2)
          call mqc_integral_allocate(eUVecs,'left singular vectors','space',tmpMat2)
          call mqc_eigenvalues_allocate(eVals,'singular values','space',tmpVec1)
        elseIf(present(eVals).and.present(eVVecs)) then
          call tmpMat1%svd(tmpVec1,eVVecs=tmpMat2)
          call mqc_integral_allocate(eVVecs,'right singular vectors','space',tmpMat2)
          call mqc_eigenvalues_allocate(eVals,'singular values','space',tmpVec1)
        elseIf(present(eUVecs).and.present(eVVecs)) then
          call tmpMat1%svd(eUVecs=tmpMat2,eVVecs=tmpMat3)
          call mqc_integral_allocate(eUVecs,'left singular vectors','space',tmpMat2)
          call mqc_integral_allocate(eVVecs,'right singular vectors','space',tmpMat3)
        elseIf(present(eUVecs)) then
          call tmpMat1%svd(EUVecs=tmpMat2)
          call mqc_integral_allocate(eUVecs,'left singular vectors','space',tmpMat2)
        elseIf(present(eVVecs)) then
          call tmpMat1%svd(EVVecs=tmpMat2)
          call mqc_integral_allocate(eVVecs,'right singular vectors','space',tmpMat2)
        elseIf(present(eVals)) then
          call tmpMat1%svd(tmpVec1)
          call mqc_eigenvalues_allocate(eVals,'singular values','space',tmpVec1)
        endIf
      elseIf(integral%type().eq.'spin') then
        tmpMat1 = integral%getblock('alpha')
        tmpMat2 = integral%getblock('beta')
        if(present(eVals).and.present(eUVecs).and.present(eVVecs)) then
          call tmpMat1%svd(tmpVec1,tmpMat3,tmpMat5)
          call tmpMat2%svd(tmpVec2,tmpMat4,tmpMat6)
          call mqc_integral_allocate(eUVecs,'left singular vectors','spin',tmpMat3,tmpMat4)
          call mqc_integral_allocate(eVVecs,'right singular vectors','spin',tmpMat5,tmpMat6)
          call mqc_eigenvalues_allocate(eVals,'singular values','spin',tmpVec1,tmpVec2)
        elseIf(present(eVals).and.present(eUVecs)) then
          call tmpMat1%svd(tmpVec1,tmpMat3)
          call tmpMat2%svd(tmpVec2,tmpMat4)
          call mqc_integral_allocate(eUVecs,'left singular vectors','spin',tmpMat3,tmpMat4)
          call mqc_eigenvalues_allocate(eVals,'singular values','spin',tmpVec1,tmpVec2)
        elseIf(present(eVals).and.present(eVVecs)) then
          call tmpMat1%svd(tmpVec1,eVVecs=tmpMat3)
          call tmpMat2%svd(tmpVec2,eVVecs=tmpMat4)
          call mqc_integral_allocate(eVVecs,'right singular vectors','spin',tmpMat3,tmpMat4)
          call mqc_eigenvalues_allocate(eVals,'singular values','spin',tmpVec1,tmpVec2)
        elseIf(present(eUVecs).and.present(eVVecs)) then
          call tmpMat1%svd(eUVecs=tmpMat3,eVVecs=tmpMat5)
          call tmpMat2%svd(eUVecs=tmpMat4,eVVecs=tmpMat6)
          call mqc_integral_allocate(eUVecs,'left singular vectors','spin',tmpMat3,tmpMat4)
          call mqc_integral_allocate(eVVecs,'right singular vectors','spin',tmpMat5,tmpMat6)
        elseIf(present(eUVecs)) then
          call tmpMat1%svd(EUVecs=tmpMat3)
          call tmpMat2%svd(EUVecs=tmpMat4)
          call mqc_integral_allocate(eUVecs,'left singular vectors','spin',tmpMat3,tmpMat4)
        elseIf(present(eVVecs)) then
          call tmpMat1%svd(EVVecs=tmpMat3)
          call tmpMat2%svd(EVVecs=tmpMat4)
          call mqc_integral_allocate(eVVecs,'right singular vectors','spin',tmpMat3,tmpMat4)
        elseIf(present(eVals)) then
          call tmpMat1%svd(tmpVec1)
          call tmpMat2%svd(tmpVec2)
          call mqc_eigenvalues_allocate(eVals,'singular values','spin',tmpVec1,tmpVec2)
        endIf
      elseIf(integral%type().eq.'general') then
        nDimAlpha1 = integral%blockSize('alpha',1) 
        nDimAlpha2 = integral%blockSize('alpha',2) 
        nDimBeta1 = integral%blockSize('beta',1) 
        nDimBeta2 = integral%blockSize('beta',2) 
        tmpMat1 = integral%getblock('full')
        if(present(eVals).and.present(eUVecs).and.present(eVVecs)) then
          call tmpMat1%svd(tmpVec1,tmpMat2,tmpMat3)
          tmpMat4 = tmpMat2%mat([1,NDimAlpha1],[1,NDimAlpha1])
          tmpMat5 = tmpMat2%mat([NDimAlpha1+1,-1],[1,NDimAlpha1])
          tmpMat6 = tmpMat2%mat([1,NDimAlpha1],[NDimAlpha1+1,-1])
          tmpMat2 = tmpMat2%mat([NDimAlpha1+1,-1],[NDimAlpha1+1,-1])
          tmpMat7 = tmpMat3%mat([1,NDimAlpha2],[1,NDimAlpha2])
          tmpMat8 = tmpMat3%mat([NDimAlpha2+1,-1],[1,NDimAlpha2])
          tmpMat1 = tmpMat3%mat([1,NDimAlpha2],[NDimAlpha2+1,-1])
          tmpMat3 = tmpMat3%mat([NDimAlpha2+1,-1],[NDimAlpha2+1,-1])
          tmpVec2 = tmpVec1%vat(1,min(NDimAlpha1,NDimAlpha2))
          tmpVec1 = tmpVec1%vat(min(NDimAlpha1,NDimAlpha2)+1,-1)
          call mqc_integral_allocate(eUVecs,'left singular vectors','general',tmpMat4,&
            tmpMat2,tmpMat5,tmpMat6)
          call mqc_integral_allocate(eVVecs,'right singular vectors','general',tmpMat7,&
            tmpMat3,tmpMat8,tmpMat1)
          call mqc_eigenvalues_allocate(eVals,'singular values','general',tmpVec2,tmpVec1)
        elseIf(present(eVals).and.present(eUVecs)) then
          call tmpMat1%svd(tmpVec1,tmpMat2)
          tmpMat3 = tmpMat2%mat([1,NDimAlpha1],[1,NDimAlpha1])
          tmpMat4 = tmpMat2%mat([NDimAlpha1+1,-1],[1,NDimAlpha1])
          tmpMat1 = tmpMat2%mat([1,NDimAlpha1],[NDimAlpha1+1,-1])
          tmpMat2 = tmpMat2%mat([NDimAlpha1+1,-1],[NDimAlpha1+1,-1])
          tmpVec2 = tmpVec1%vat(1,min(NDimAlpha1,NDimAlpha2))
          tmpVec1 = tmpVec1%vat(min(NDimAlpha1,NDimAlpha2)+1,-1)
          call mqc_integral_allocate(eUVecs,'left singular vectors','general',tmpMat3,&
            tmpMat2,tmpMat4,tmpMat1)
          call mqc_eigenvalues_allocate(eVals,'singular values','general',tmpVec2,tmpVec1)
        elseIf(present(eVals).and.present(eVVecs)) then
          call tmpMat1%svd(tmpVec1,eVVecs=tmpMat2)
          tmpMat3 = tmpMat2%mat([1,NDimAlpha2],[1,NDimAlpha2])
          tmpMat4 = tmpMat2%mat([NDimAlpha2+1,-1],[1,NDimAlpha2])
          tmpMat1 = tmpMat2%mat([1,NDimAlpha2],[NDimAlpha2+1,-1])
          tmpMat2 = tmpMat2%mat([NDimAlpha2+1,-1],[NDimAlpha2+1,-1])
          tmpVec2 = tmpVec1%vat(1,min(NDimAlpha1,NDimAlpha2))
          tmpVec1 = tmpVec1%vat(min(NDimAlpha1,NDimAlpha2)+1,-1)
          call mqc_integral_allocate(eUVecs,'right singular vectors','general',tmpMat3,&
            tmpMat2,tmpMat4,tmpMat1)
          call mqc_eigenvalues_allocate(eVals,'singular values','general',tmpVec2,tmpVec1)
        elseIf(present(eUVecs).and.present(eVVecs)) then
          call tmpMat1%svd(eUVecs=tmpMat2,eVVecs=tmpMat3)
          tmpMat4 = tmpMat2%mat([1,NDimAlpha1],[1,NDimAlpha1])
          tmpMat5 = tmpMat2%mat([NDimAlpha1+1,-1],[1,NDimAlpha1])
          tmpMat6 = tmpMat2%mat([1,NDimAlpha1],[NDimAlpha1+1,-1])
          tmpMat2 = tmpMat2%mat([NDimAlpha1+1,-1],[NDimAlpha1+1,-1])
          tmpMat7 = tmpMat3%mat([1,NDimAlpha2],[1,NDimAlpha2])
          tmpMat8 = tmpMat3%mat([NDimAlpha2+1,-1],[1,NDimAlpha2])
          tmpMat1 = tmpMat3%mat([1,NDimAlpha2],[NDimAlpha2+1,-1])
          tmpMat3 = tmpMat3%mat([NDimAlpha2+1,-1],[NDimAlpha2+1,-1])
          call mqc_integral_allocate(eUVecs,'left singular vectors','general',tmpMat4,&
            tmpMat2,tmpMat5,tmpMat6)
          call mqc_integral_allocate(eVVecs,'right singular vectors','general',tmpMat7,&
            tmpMat3,tmpMat8,tmpMat1)
        elseIf(present(eUVecs)) then
          call tmpMat1%svd(EUVecs=tmpMat2)
          tmpMat3 = tmpMat2%mat([1,NDimAlpha1],[1,NDimAlpha1])
          tmpMat4 = tmpMat2%mat([NDimAlpha1+1,-1],[1,NDimAlpha1])
          tmpMat1 = tmpMat2%mat([1,NDimAlpha1],[NDimAlpha1+1,-1])
          tmpMat2 = tmpMat2%mat([NDimAlpha1+1,-1],[NDimAlpha1+1,-1])
          call mqc_integral_allocate(eUVecs,'left singular vectors','general',tmpMat3,&
            tmpMat2,tmpMat4,tmpMat1)
        elseIf(present(eVVecs)) then
          call tmpMat1%svd(EVVecs=tmpMat2)
          tmpMat3 = tmpMat2%mat([1,NDimAlpha2],[1,NDimAlpha2])
          tmpMat4 = tmpMat2%mat([NDimAlpha2+1,-1],[1,NDimAlpha2])
          tmpMat1 = tmpMat2%mat([1,NDimAlpha2],[NDimAlpha2+1,-1])
          tmpMat2 = tmpMat2%mat([NDimAlpha2+1,-1],[NDimAlpha2+1,-1])
          call mqc_integral_allocate(eUVecs,'right singular vectors','general',tmpMat3,&
            tmpMat2,tmpMat4,tmpMat1)
        elseIf(present(eVals)) then
          call tmpMat1%svd(tmpVec1)
          tmpVec2 = tmpVec1%vat(1,min(NDimAlpha1,NDimAlpha2))
          tmpVec1 = tmpVec1%vat(min(NDimAlpha1,NDimAlpha2)+1,-1)
          call mqc_eigenvalues_allocate(eVals,'singular values','general',tmpVec2,tmpVec1)
        endIf
      endIf
!
      end subroutine mqc_scf_integral_svd
!
!
!     PROCEDURE MQC_SCF_Integral_Inverse    
!
!>    \brief <b> MQC_SCF_Integral_Inverse is a function that returns the inverse of 
!>    an MQC integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Integral_Inverse is a function that returns the inverse of an MQC 
!>    integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral which will be evaluated.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_scf_integral_inverse(integral) result(inverse)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      type(mqc_scf_integral)::inverse
      type(mqc_matrix)::tmpMat1,tmpMat2,tmpMat3,tmpMat4
      integer(kind=int64)::nDimAlpha=0
!
      if(integral%type().eq.'space') then
        tmpMat1 = integral%getblock('alpha')
        tmpMat1 = tmpMat1%inv()
        call mqc_integral_allocate(inverse,'integral inverse','space',tmpMat1)
      elseIf(integral%type().eq.'spin') then
        tmpMat1 = integral%getblock('alpha')
        tmpMat2 = integral%getblock('beta')
        tmpMat1 = tmpMat1%inv()
        tmpMat2 = tmpMat2%inv()
        call mqc_integral_allocate(inverse,'integral inverse','spin',tmpMat1,tmpMat2)
      elseIf(integral%type().eq.'general') then
        nDimAlpha = integral%blockSize('alpha') 
        tmpMat1 = integral%getblock('full')
        tmpMat1 = tmpMat1%inv()
        tmpMat2 = tmpMat1%mat([NDimAlpha+1,-1],[NDimAlpha+1,-1])
        tmpMat3 = tmpMat1%mat([NDimAlpha+1,-1],[1,NDimAlpha])
        tmpMat4 = tmpMat1%mat([1,NDimAlpha],[NDimAlpha+1,-1])
        tmpMat1 = tmpMat1%mat([1,NDimAlpha],[1,NDimAlpha])
        call mqc_integral_allocate(inverse,'integral inverse','general',tmpMat1,tmpMat2,tmpMat3,tmpMat4)
      endIf
!
      end function mqc_scf_integral_inverse
!
!
!     PROCEDURE MQC_SCF_Integral_Trace
!
!>    \brief <b> MQC_SCF_Integral_Trace is a function that returns the trace of an 
!>    MQC SCF integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Integral_Trace is a function that returns the trace of an  MQC SCF 
!>    integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral which will be evaluated.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_scf_integral_trace(integral) result(trace)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      type(mqc_scalar)::trace,two
      type(mqc_matrix)::tmpMat1,tmpMat2
!
      two = 2.0
      if(integral%type().eq.'space') then
        tmpMat1 = integral%getblock('alpha')
        trace = two*tmpMat1%trace()
      elseIf(integral%type().eq.'spin'.or.integral%type().eq.'general') then
        tmpMat1 = integral%getblock('alpha')
        tmpMat2 = integral%getblock('beta')
        trace = tmpMat1%trace() + tmpMat2%trace()
      endIf
!
      end function mqc_scf_integral_trace
!
!
!     PROCEDURE MQC_SCF_Integral_Determinant
!
!>    \brief <b> MQC_SCF_Integral_Determinant is a function that returns the 
!>    determinant of an MQC integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Integral_Determinant is a function that returns the determinant of an 
!>    MQC integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral which will be evaluated.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author A. D. Mahler
!>    \date 2018
!
      function mqc_scf_integral_determinant(integral) result(determinant)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      type(mqc_scalar)::determinant
      type(mqc_matrix)::tmpMat1
      integer(kind=int64)::nDimAlpha=0
!
      tmpMat1 = integral%getBlock('full')
      determinant = tmpMat1%det()
!
      end function mqc_scf_integral_determinant
!
!
!     PROCEDURE MQC_SCF_Integral_Conjg
!
!>    \brief <b> MQC_SCF_Integral_Conjg is a function that returns the complex conjugate 
!>    of an MQC SCF integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Integral_Conjg is a function that returns the complex conjugate of an  MQC SCF 
!>    integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral which will be evaluated.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        Label for the output SCF integral.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2022
!
      function mqc_scf_integral_conjg(integral,label) result(integralOut)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      type(mqc_matrix)::tmpMat1,tmpMat2,tmpMat3,tmpMat4
      Character(Len=64)::myLabel
!
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      if(integral%type().eq.'space') then
        tmpMat1 = conjg(integral%getblock('alpha'))
        call mqc_integral_allocate(integralOut,myLabel,'space',tmpMat1)
      elseIf(integral%type().eq.'spin') then
        tmpMat1 = conjg(integral%getblock('alpha'))
        tmpMat2 = conjg(integral%getblock('beta'))
        call mqc_integral_allocate(integralOut,myLabel,'spin',tmpMat1,tmpMat2)
      elseIf(integral%type().eq.'general') then
        tmpMat1 = conjg(integral%getblock('alpha'))
        tmpMat2 = conjg(integral%getblock('beta'))
        tmpMat3 = conjg(integral%getblock('alpha-beta'))
        tmpMat4 = conjg(integral%getblock('beta-alpha'))
        call mqc_integral_allocate(integralOut,myLabel,'general',tmpMat1,tmpMat2,tmpMat3,tmpMat4)
      endIf
!
      end function mqc_scf_integral_conjg
!
!
!     PROCEDURE MQC_Integral_Set_Energy_List
!
!>    \brief <b> MQC_Integral_Set_Energy_List is a subroutine that sets the energy
!>    list object of an MQC integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Set_Energy_List is a subroutine that sets the energy list object 
!>    of an MQC integral
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral for which EList will be set.
!>    \endverbatim
!>
!>    \param[in] EList
!>    \verbatim
!>        EList is integer(kind=int64),dimension(:),allocatable
!>        The energy order list to set EList in Integral.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author A. D. Mahler
!>    \date 2018
!
      subroutine mqc_integral_set_energy_list(integral,elist)
!
      implicit none
      class(mqc_scf_integral)::integral
      integer(kind=int64),dimension(:),allocatable::elist

      if(allocated(integral%Energy_List)) deallocate(integral%Energy_List)

      integral%Energy_List = elist
!
      end subroutine mqc_integral_set_energy_list
!
!
!     PROCEDURE MQC_Integral_Get_Energy_List
!
!>    \brief <b> MQC_Integral_Get_Energy_List is a subroutine that returns the energy
!>    list object of an MQC integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Get_Energy_List is a subroutine that returns the energy list object 
!>    of an MQC integral
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral from which EList will be returned.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author A. D. Mahler
!>    \date 2018
!
      function mqc_integral_get_energy_list(integral) result(elist)
!
      implicit none
      class(mqc_scf_integral)::integral
      integer(kind=int64),dimension(:),allocatable::elist

      if(.not.allocated(integral%Energy_List)) then
        call mqc_error('Energy_List is not set in MQC_Integral object')
      else
        elist = integral%Energy_List
      endif
!
      end function mqc_integral_get_energy_list
!
!
!     PROCEDURE MQC_Integral_Delete_Energy_List
!
!>    \brief <b> MQC_Integral_Delete_Energy_List is a subroutine that deletes the 
!>    energy list object of an MQC integral</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_Delete_Energy_List is a subroutine that deletes the energy list 
!>    object of an MQC integral
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral from which EList will be deleted.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author A. D. Mahler
!>    \date 2018
!
      subroutine mqc_integral_delete_energy_list(integral)
!
      implicit none
      class(mqc_scf_integral)::integral

      if(.not.allocated(integral%Energy_List)) then
        call mqc_error('Energy_List was not set in MQC_Integral object')
      else
        deallocate(integral%Energy_List)
      end if
!
      end subroutine mqc_integral_delete_energy_list
!
!
!     PROCEDURE MQC_SCF_Eigenvalues_Power
!
!>    \brief <b> MQC_SCF_Eigenvalues_Power is a function that returns the value of 
!>    all elements of an MQC eigenvalues variable raised to a power</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Eigenvalues_Power is a function that returns the value of all elements 
!>    of an MQC eigenvalues variable raised to a power.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_SCF_Eigenvalues)
!>        The name of the MQC_SCF_Eigenvalues variable.
!>    \endverbatim
!>
!>    \param[in] Power
!>    \verbatim
!>        Power is Class(*)
!>        The power to raise elements of the MQC eigenvalues.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2019
!
      subroutine mqc_scf_eigenvalues_power(eigenvalues,power)
!
      implicit none
      class(mqc_scf_eigenvalues),intent(inOut)::eigenvalues
      class(*)::power
      type(mqc_scalar)::scalar
      type(mqc_vector)::tmpVec1,tmpVec2
!
      select type(power)
      type is (integer)
        scalar = power
      type is (real)
        scalar = power
      type is (complex)
        scalar = power
      type is (mqc_scalar)
        scalar = power
      class default
        call mqc_error_I('power type not defined in MQC_SCF_Eigenvalues_Power',6)
      end select

      if(eigenvalues%type().eq.'space') then
        tmpVec1 = eigenvalues%getblock('alpha')
        call tmpVec1%power(scalar)
        call mqc_eigenvalues_allocate(eigenvalues,'eigenvalues exponentiated','space',tmpVec1)
      elseIf(eigenvalues%type().eq.'spin'.or.eigenvalues%type().eq.'general') then
        tmpVec1 = eigenvalues%getblock('alpha')
        tmpVec2 = eigenvalues%getblock('beta')
        call tmpVec1%power(scalar)
        call tmpVec2%power(scalar)
        call mqc_eigenvalues_allocate(eigenvalues,'eigenvalues exponentiated','spin',tmpVec1,tmpVec2)
      endIf
!
      end subroutine mqc_scf_eigenvalues_power   
!
!
!     PROCEDURE MQC_SCF_Integral_Power
!
!>    \brief <b> MQC_SCF_Integral_Power is a function that returns the value of 
!>    an MQC integral variable raised to a power</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_SCF_Integral_Power is a function that returns an MQC
!>    integral variable raised to a power.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The name of the MQC_SCF_Integral variable.
!>    \endverbatim
!>
!>    \param[in] Power
!>    \verbatim
!>        Power is Class(*)
!>        The power to raise the MQC integral.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is Character(Len=*),optional
!>        Label for the output SCF integral.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2022
!
      subroutine mqc_scf_integral_power(integral,power,label)
!
      implicit none
      class(mqc_scf_integral),intent(inOut)::integral
      class(*)::power
      type(mqc_scalar)::scalar
      type(mqc_matrix)::tmpMat1,tmpMat2,tmpMat3,tmpMat4
      integer(kind=int64)::nDimAlpha=0
      Character(Len=*),optional,intent(in)::label
      Character(Len=64)::myLabel
!
      select type(power)
      type is (integer)
        scalar = power
      type is (real)
        scalar = power
      type is (complex)
        scalar = power
      type is (mqc_scalar)
        scalar = power
      class default
        call mqc_error_I('power type not defined in MQC_SCF_Integral_Power',6)
      end select

      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = ''
      endIf

      if(integral%type().eq.'space') then
        tmpMat1 = integral%getblock('alpha')
        call tmpMat1%power(scalar)
        call mqc_integral_allocate(integral,myLabel,'space',tmpMat1)
      elseIf(integral%type().eq.'spin') then
        tmpMat1 = integral%getblock('alpha')
        tmpMat2 = integral%getblock('beta')
        call tmpMat1%power(scalar)
        call tmpMat2%power(scalar)
        call mqc_integral_allocate(integral,myLabel,'spin',tmpMat1,tmpMat2)
      elseIf(integral%type().eq.'general') then
        nDimAlpha = integral%blockSize('alpha') 
        tmpMat1 = integral%getblock('full')
        call tmpMat1%power(scalar)
        tmpMat2 = tmpMat1%mat([NDimAlpha+1,-1],[NDimAlpha+1,-1])
        tmpMat3 = tmpMat1%mat([NDimAlpha+1,-1],[1,NDimAlpha])
        tmpMat4 = tmpMat1%mat([1,NDimAlpha],[NDimAlpha+1,-1])
        tmpMat1 = tmpMat1%mat([1,NDimAlpha],[1,NDimAlpha])
        call mqc_integral_allocate(integral,myLabel,'general',tmpMat1,tmpMat2,tmpMat3,tmpMat4)
      endIf
!
      end subroutine mqc_scf_integral_power   
!
!
!     PROCEDURE MQC_twoERIs_At
!
!>    \brief <b> MQC_twoERIs_At is a function that returns the value of an element 
!>    of a MQC 2ERIs variable</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_twoERIs_At is a function that returns the value of (I,J,K,L)th element of 
!>    a MQC 2ERIs variable as an MQC scalar. If I, J, K or L is negative, the (N-I+1)th 
!>    index is selected. Optional argument SpinBlock specifies the spin block of MQC
!>    2ERIs if the contain MO integrals. If wavefunction is not restricted, the spin 
!>    block must be specified. The following options are available:
!>
!>    1.  SpinBlock = 'aaaa' or 'alpha' or 'alpha-alpha' uses the alpha spin block.
!>    2.  SpinBlock = 'bbbb' or 'beta' or 'beta-beta' uses the beta spin block.
!>    3.  SpinBlock = 'aabb' or 'alphaBeta' uses the alpha-beta spin block.
!>    4.  SpinBlock = 'bbaa' or 'betaAlpha' uses the beta-alpha spin block.
!>    5.  SpinBlock = 'abab' uses the abab spin block.
!>    6.  SpinBlock = 'abba' uses the abba spin block.
!>    7.  SpinBlock = 'baab' uses the baab spin block.
!>    8.  SpinBlock = 'baba' uses the baba spin block.
!>    9.  SpinBlock = 'aaab' uses the aaab spin block.
!>    10. SpinBlock = 'aaba' uses the aaba spin block.
!>    11. SpinBlock = 'abaa' uses the abaa spin block.
!>    12. SpinBlock = 'baaa' uses the baaa spin block.
!>    13. SpinBlock = 'bbba' uses the bbba spin block.
!>    14. SpinBlock = 'bbab' uses the bbab spin block.
!>    15. SpinBlock = 'babb' uses the babb spin block.
!>    16. SpinBlock = 'abbb' uses the abbb spin block.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] TwoERIs
!>    \verbatim
!>        TwoERIs is Class(MQC_TwoERIs)
!>        The MQC 2ERIs to return the value of the (I,J,K,L)th 
!>        element.
!>    \endverbatim
!>    
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The first dimension of the element in MQC 2ERIs
!>          If I>0 row count is from first index
!>          If I<0 row count is from last index.
!>    \endverbatim
!>    
!>    \param[in] J  
!>    \verbatim
!>        J is Integer(kind=int64)
!>        The second dimension of the element in MQC 2ERIs
!>          If J>0 row count is from first index
!>          If J<0 row count is from last index.
!>    \endverbatim
!>    
!>    \param[in] K  
!>    \verbatim
!>        K is Integer(kind=int64)
!>        The third dimension of the element in MQC 2ERIs
!>          If K>0 row count is from first index
!>          If K<0 row count is from last index.
!>    \endverbatim
!>    
!>    \param[in] L  
!>    \verbatim
!>        L is Integer(kind=int64)
!>        The fourth dimension of the element in MQC 2ERIs
!>          If L>0 row count is from first index
!>          If L<0 row count is from last index.
!>    \endverbatim
!>    
!>    \param[in] SpinBlockIn
!>    \verbatim
!>        SpinBlockIn is character(len=*),optional
!>        = 'aaaa' or 'alpha':     alpha spin block
!>        = 'bbbb' or 'beta':      beta spin block
!>        = 'aabb' or 'alphaBeta': alpha-beta spin block
!>        = 'bbaa' or 'betaAlpha': beta-alpha spin block.
!>        = 'abab':                abab spin block.
!>        = 'abba':                abba spin block.
!>        = 'baab':                baab spin block.
!>        = 'baba':                baba spin block.
!>        = 'aaab':                aaab spin block.
!>        = 'aaba':                aaba spin block.
!>        = 'abaa':                abaa spin block.
!>        = 'baaa':                baaa spin block.
!>        = 'bbba':                bbba spin block.
!>        = 'bbab':                bbab spin block.
!>        = 'babb':                babb spin block.
!>        = 'abbb':                abbb spin block.
!>    \endverbatim
!>    
!>    \param[in] InteractionIn
!>    \verbatim
!>        InteractionIn is character(len=*),optional
!>        = 'coulomb':   returns (ij|kl)
!>        = 'exchange':  returns (il|kj)
!>        = 'doublebar': returns (ij|kl) - (il|kj).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function MQC_twoERIs_At(twoERIs,i,j,k,l,spinBlockIn,interactionIn) result(element)
!
      implicit none
      class(mqc_twoERIs)::twoERIs
      class(*),intent(in)::i,j,k,l
      integer(kind=int64)::myI,myJ,myK,myL
      type(mqc_scalar)::element
      character(len=*),optional,intent(in)::spinBlockIn,interactionIn
      character(len=64)::spinBlock,interaction
      real(kind=real64)::zero=0.0d0

      element = zero
      select type (i)  
      type is (integer)
        myI = i
      type is (mqc_scalar)
        myI = i
      class default 
        call mqc_error('index type not defined in MQC_twoERIs_at')
      end select
      select type (j)  
      type is (integer)
        myJ = j
      type is (mqc_scalar)
        myJ = j
      class default 
        call mqc_error('index type not defined in MQC_twoERIs_at')
      end select
      select type (k)  
      type is (integer)
        myK = k
      type is (mqc_scalar)
        myK = k
      class default 
        call mqc_error('index type not defined in MQC_twoERIs_at')
      end select
      select type (l)  
      type is (integer)
        myL = l
      type is (mqc_scalar)
        myL = l
      class default 
        call mqc_error('index type not defined in MQC_twoERIs_at')
      end select

      if(myI.gt.twoERIs%blockSize('alpha',1)+twoERIs%blockSize('beta',1)) &
        call mqc_error('Requested I dimension out of bounds in MQC_twoERIs_At')
      if(myJ.gt.twoERIs%blockSize('alpha',2)+twoERIs%blockSize('beta',2)) &
        call mqc_error('Requested J dimension out of bounds in MQC_twoERIs_At')
      if(myK.gt.twoERIs%blockSize('alpha',3)+twoERIs%blockSize('beta',3)) &
        call mqc_error('Requested K dimension out of bounds in MQC_twoERIs_At')
      if(myL.gt.twoERIs%blockSize('alpha',4)+twoERIs%blockSize('beta',4)) &
        call mqc_error('Requested L dimension out of bounds in MQC_twoERIs_At')

      if(present(interactionIn)) then
        call string_change_case(interactionIn,'L',interaction)
      else
        interaction = 'coulomb'
      endIf
      if(interaction.ne.'coulomb'.and.interaction.ne.'exchange'.and.interaction.ne.'doublebar') &
        call mqc_error_a('Interaction type requested not recognized in mqc_twoERIs_at',6,&
        'interaction',interaction)

      if(present(spinBlockIn)) then
!       SPACE
        call string_change_case(spinBlockIn,'L',spinBlock)
        if(twoERIs%integralType.eq.'space') then
          if(spinBlock.eq.'alpha'.or.spinBlock.eq.'aaaa'.or.&
            spinBlock.eq.'beta'.or.spinBlock.eq.'bbbb') then
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL)
              elseIf(interaction.eq.'exchange') then
                element = twoERIs%alpha%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL) - twoERIs%alpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseIf(spinBlock.eq.'alphabeta'.or.spinBlock.eq.'aabb'.or.&
            spinBlock.eq.'betaalpha'.or.spinBlock.eq.'bbaa') then
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL)
              endIf
            endIf
          elseIf(spinBlock.eq.'abba'.or.spinBlock.eq.'baab') then
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'exchange') then
                element = twoERIs%alpha%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = (-1)*twoERIs%alpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          endIf
!       SPIN
        elseIf(twoERIs%integralType.eq.'spin') then
          if(spinBlock.eq.'alpha'.or.spinBlock.eq.'aaaa') then
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL)
              elseIf(interaction.eq.'exchange') then
                element = twoERIs%alpha%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL) - twoERIs%alpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinblock.eq.'beta'.or.spinBlock.eq.'bbbb') then
            if(twoERIs%hasSpinBlock('bbbb')) then
              if(interaction.eq.'coulomb') then
                element = twoeris%beta%at(myI,myJ,myK,myL)
              elseIf(interaction.eq.'exchange') then
                element = twoERIs%beta%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%beta%at(myI,myJ,myK,myL) - twoERIs%beta%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinblock.eq.'alphabeta'.or.spinBlock.eq.'aabb') then
            if(twoERIs%hasSpinBlock('aabb')) then
              if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
                element = twoeris%alphaBeta%at(myI,myJ,myK,myL)
              endIf
            endIf
          elseif(spinblock.eq.'betaalpha'.or.spinBlock.eq.'bbaa') then
            if(twoERIs%storageType.eq.'symm') then
              if(twoERIs%hasSpinBlock('aabb')) then
                if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
!                  element = twoeris%alphaBeta%at(myI,myJ,myK,myL)
                  element = twoeris%alphaBeta%at(myK,myL,myI,myJ)
                endIf
              endIf
            else
              if(twoERIs%hasSpinBlock('bbaa')) then
                if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
                  element = twoeris%betaAlpha%at(myI,myJ,myK,myL)
                endIf
              endIf
            endIf
          elseif(spinBlock.eq.'abba') then
            if(twoERIs%hasSpinBlock('aabb')) then
              if(interaction.eq.'exchange') then
                element = twoeris%alphaBeta%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                element = (-1)*twoeris%alphaBeta%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'baab') then
            if(twoERIs%storageType.eq.'symm') then
              if(twoERIs%hasSpinBlock('aabb')) then
                if(interaction.eq.'exchange') then
!                  element = twoeris%alphaBeta%at(myI,myL,myK,myJ)
                  element = twoeris%alphaBeta%at(myK,myJ,myI,myL)
                elseif(interaction.eq.'doublebar') then
!                  element = (-1)*twoeris%alphaBeta%at(myI,myL,myK,myJ)
                  element = (-1)*twoeris%alphaBeta%at(myK,myJ,myI,myL)
                endIf
              endIf
            else
              if(twoERIs%hasSpinBlock('bbaa')) then
                if(interaction.eq.'exchange') then
                  element = twoeris%betaAlpha%at(myI,myL,myK,myJ)
                elseif(interaction.eq.'doublebar') then
                  element = (-1)*twoeris%betaAlpha%at(myI,myL,myK,myJ)
                endIf
              endIf
            endIf
          endIf
!       GENERAL
        elseIf(twoERIs%integralType.eq.'general') then
          if(spinBlock.eq.'alpha'.or.spinBlock.eq.'aaaa') then
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL)
              elseIf(interaction.eq.'exchange') then
                element = twoERIs%alpha%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL) - twoERIs%alpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinblock.eq.'beta'.or.spinBlock.eq.'bbbb') then
            if(twoERIs%hasSpinBlock('bbbb')) then
              if(interaction.eq.'coulomb') then
                element = twoeris%beta%at(myI,myJ,myK,myL)
              elseIf(interaction.eq.'exchange') then
                element = twoERIs%beta%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%beta%at(myI,myJ,myK,myL) - twoERIs%beta%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinblock.eq.'alphabeta'.or.spinBlock.eq.'aabb') then
            if(interaction.eq.'coulomb') then
              if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myI,myJ,myK,myL)
            elseif(interaction.eq.'exchange') then
              if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myI,myL,myK,myJ)
            elseif(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myI,myJ,myK,myL)
              if(twoERIs%hasSpinBlock('abba')) element = element - twoeris%abba%at(myI,myL,myK,myJ)
            endIf
          elseif(spinblock.eq.'betaalpha'.or.spinBlock.eq.'bbaa') then
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
!                if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myK,myL,myI,myJ)
              elseif(interaction.eq.'exchange') then
!                if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myK,myJ,myI,myL)
              elseif(interaction.eq.'doublebar') then
!                if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myI,myJ,myK,myL)
!                if(twoERIs%hasSpinBlock('abba')) element = element - twoeris%abba%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myK,myL,myI,myJ)
                if(twoERIs%hasSpinBlock('abba')) element = element - twoeris%abba%at(myK,myJ,myI,myL)
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('bbaa')) element = twoeris%betaAlpha%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('baab')) element = twoeris%baab%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('bbaa')) element = twoeris%betaAlpha%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('baab')) element = element - twoeris%baab%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'abab') then
            if(twoERIs%hasSpinBlock('abab')) then
              if(interaction.eq.'coulomb') then
                element = twoeris%abab%at(myI,myJ,myK,myL)
              elseIf(interaction.eq.'exchange') then
                element = twoeris%abab%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%abab%at(myI,myJ,myK,myL) - twoERIs%abab%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'abba') then
            if(interaction.eq.'coulomb') then
              if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myI,myJ,myK,myL)
            elseif(interaction.eq.'exchange') then
              if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myI,myL,myK,myJ)
            elseif(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myI,myJ,myK,myL)
              if(twoERIs%hasSpinBlock('aabb')) element = element - twoeris%alphaBeta%at(myI,myL,myK,myJ)
            endIf
          elseif(spinBlock.eq.'baab') then
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
!                if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myK,myL,myI,myJ)
              elseif(interaction.eq.'exchange') then
!                if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myK,myJ,myI,myL)
              elseif(interaction.eq.'doublebar') then
!                if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myI,myJ,myK,myL)
!                if(twoERIs%hasSpinBlock('aabb')) element = element - twoeris%alphaBeta%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myK,myL,myI,myJ)
                if(twoERIs%hasSpinBlock('aabb')) element = element - twoeris%alphaBeta%at(myK,myJ,myI,myL)
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('baab')) element = twoeris%baab%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('bbaa')) element = twoeris%betaAlpha%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('baab')) element = twoeris%baab%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('bbaa')) element = element - twoeris%betaAlpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'baba') then
            if(twoERIs%storageType.eq.'symm') then
              if(twoERIs%hasSpinBlock('abab')) then
                if(interaction.eq.'coulomb') then
!                  element = conjg(twoeris%abab%at(myI,myJ,myK,myL))
                  element = conjg(twoeris%abab%at(myJ,myI,myL,myK))
                elseIf(interaction.eq.'exchange') then
!                  element = conjg(twoeris%abab%at(myI,myL,myK,myJ))
                  element = conjg(twoeris%abab%at(myL,myI,myJ,myK))
                elseIf(interaction.eq.'doublebar') then
!                  element = conjg(twoERIs%abab%at(myI,myJ,myK,myL)) - conjg(twoERIs%abab%at(myI,myL,myK,myJ))
                  element = conjg(twoERIs%abab%at(myJ,myI,myL,myK)) - conjg(twoERIs%abab%at(myL,myI,myJ,myK))
                endIf
              endIf
            else
              if(twoERIs%hasSpinBlock('baba')) then
                if(interaction.eq.'coulomb') then
                  element = twoeris%baba%at(myI,myJ,myK,myL)
                elseIf(interaction.eq.'exchange') then
                  element = twoeris%baba%at(myI,myL,myK,myJ)
                elseIf(interaction.eq.'doublebar') then
                  element = twoERIs%baba%at(myI,myJ,myK,myL) - twoERIs%baba%at(myI,myL,myK,myJ)
                endIf
              endIf
            endIf
          elseif(spinBlock.eq.'aaab') then
            if(interaction.eq.'coulomb') then
              if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myJ,myK,myL)
            elseif(interaction.eq.'exchange') then
              if(twoERIs%storageType.eq.'symm') then
!                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myK,myJ,myI,myL)
              else
                if(twoERIs%hasSpinBlock('abaa')) element = twoeris%abaa%at(myI,myL,myK,myJ)
              endIf
            elseif(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myJ,myK,myL)
              if(twoERIs%storageType.eq.'symm') then
!                if(twoERIs%hasSpinBlock('aaab')) element = element - twoeris%aaab%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('aaab')) element = element - twoeris%aaab%at(myK,myJ,myI,myL)
              else
                if(twoERIs%hasSpinBlock('abaa')) element = element - twoeris%abaa%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'aaba') then
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
!                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myJ,myK,myL))
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myJ,myL,myK))
              elseif(interaction.eq.'exchange') then
!                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myL,myK,myJ))
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myL,myJ,myK))
              elseif(interaction.eq.'doublebar') then
!                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myJ,myK,myL)) - &
!                  conjg(twoeris%aaab%at(myI,myL,myK,myJ))
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myJ,myL,myK)) - &
                  conjg(twoeris%aaab%at(myI,myL,myJ,myK))
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('aaba')) element = twoeris%aaba%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aaba')) element = twoeris%aaba%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('aaba')) element = twoeris%aaba%at(myI,myJ,myK,myL) - &
                  twoeris%aaab%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'abaa') then
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
!                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myK,myL,myI,myJ)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
!                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myJ,myK,myL) - &
!                  twoeris%aaab%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myK,myL,myI,myJ) - &
                  twoeris%aaab%at(myI,myL,myK,myJ)
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('abaa')) element = twoeris%abaa%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('abaa')) element = twoeris%abaa%at(myI,myJ,myK,myL) 
                if(twoERIs%hasSpinBlock('aaab')) element = element - twoeris%aaab%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'baaa') then
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
!                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myJ,myK,myL))
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myL,myK,myI,myJ))
              elseif(interaction.eq.'exchange') then
!                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myL,myK,myJ))
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myJ,myK,myI,myL))
              elseif(interaction.eq.'doublebar') then
!                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myI,myJ,myK,myL)) - &
!                  conjg(twoeris%aaab%at(myI,myL,myK,myJ))
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myL,myK,myI,myJ)) - &
                  conjg(twoeris%aaab%at(myJ,myK,myI,myL))
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('baaa')) element = twoeris%baaa%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('baaa')) element = twoeris%baaa%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('baaa')) element = twoeris%baaa%at(myI,myJ,myK,myL) - &
                  twoeris%baaa%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'bbba') then
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
!                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myI,myJ,myK,myL))
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myK,myL,myJ,myI))
              elseif(interaction.eq.'exchange') then
!                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myI,myL,myK,myJ))
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myL,myI,myK,myJ))
              elseif(interaction.eq.'doublebar') then
!                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myI,myJ,myK,myL)) - &
!                  conjg(twoeris%abbb%at(myI,myL,myK,myJ))
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myK,myL,myJ,myI)) - &
                  conjg(twoeris%abbb%at(myL,myI,myK,myJ))
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('bbba')) element = twoeris%bbba%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('babb')) element = twoeris%babb%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('bbba')) element = twoeris%bbba%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('babb')) element = element - twoeris%babb%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'bbab') then
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
!                if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myK,myL,myI,myJ)
              elseif(interaction.eq.'exchange') then
!                if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myK,myJ,myI,myL)
              elseif(interaction.eq.'doublebar') then
!                if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myI,myJ,myK,myL) - &
!                  twoeris%abbb%at(myI,myL,myK,myJ)
                if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myK,myL,myI,myJ) - &
                  twoeris%abbb%at(myK,myJ,myI,myL)
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('bbab')) element = twoeris%bbab%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('bbab')) element = twoeris%bbab%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('bbab')) element = twoeris%bbab%at(myI,myJ,myK,myL) - &
                  twoeris%bbab%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'babb') then
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
!                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myI,myJ,myK,myL))
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myJ,myI,myK,myL))
              elseif(interaction.eq.'exchange') then
!                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myI,myL,myK,myJ))
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myK,myJ,myL,myI))
              elseif(interaction.eq.'doublebar') then
!                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myI,myJ,myK,myL)) - &
!                  conjg(twoeris%abbb%at(myI,myL,myK,myJ))
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myJ,myI,myK,myL)) - &
                  conjg(twoeris%abbb%at(myK,myJ,myL,myI))
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('babb')) element = twoeris%babb%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('bbba')) element = twoeris%bbba%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('babb')) element = twoeris%babb%at(myI,myJ,myK,myL) 
                if(twoERIs%hasSpinBlock('bbba')) element = twoeris%bbba%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseif(spinBlock.eq.'abbb') then
            if(interaction.eq.'coulomb') then
              if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myI,myJ,myK,myL)
            elseif(interaction.eq.'exchange') then
              if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myI,myL,myK,myJ)
            elseif(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myI,myJ,myK,myL) - &
                twoeris%abbb%at(myI,myL,myK,myJ)
            endIf
          else
            call mqc_error_A('SpinBlock type not valid in MQC_TwoERIs_At',6,'SpinBlock',SpinBlock)
          endIf
        else
          call mqc_error_A('Integral type not valid in MQC_TwoERIs_At',6,'integralType',twoERIs%integralType)
        endIf
      else
        if(twoERIs%integralType.eq.'regular') then
          if(myI.gt.twoERIs%blockSize('alpha',1)) myI = myI - twoERIs%blockSize('alpha',1) 
          if(myJ.gt.twoERIs%blockSize('alpha',2)) myJ = myJ - twoERIs%blockSize('alpha',2) 
          if(myK.gt.twoERIs%blockSize('alpha',3)) myK = myK - twoERIs%blockSize('alpha',3) 
          if(myL.gt.twoERIs%blockSize('alpha',4)) myL = myL - twoERIs%blockSize('alpha',4) 
          if(twoERIs%hasSpinBlock('aaaa')) then
            if(interaction.eq.'coulomb') then
              element = twoERIs%alpha%at(myI,myJ,myK,myL)
            elseif(interaction.eq.'exchange') then
              element = twoERIs%alpha%at(myI,myL,myK,myJ)
            elseif(interaction.eq.'doublebar') then
              element = twoERIs%alpha%at(myI,myJ,myK,myL) - twoERIs%alpha%at(myI,myL,myK,myJ)
            endIf
          endIf
        elseIf(twoERIs%integralType.eq.'raffenetti1') then
          if(myI.gt.twoERIs%blockSize('alpha',1).or.myJ.gt.twoERIs%blockSize('alpha',2).or. &
            myK.gt.twoERIs%blockSize('alpha',3).or.myL.gt.twoERIs%blockSize('alpha',4)) &
            call mqc_error('Requested position is out of bounds in Raffenetti 1 ERIs')
          if(twoERIs%hasSpinBlock('aaaa')) element = twoERIs%alpha%at(myI,myJ,myK,myL)
        elseIf(twoERIs%integralType.eq.'raffenetti2') then
          if(myI.gt.twoERIs%blockSize('alpha',1).or.myJ.gt.twoERIs%blockSize('alpha',2).or. &
            myK.gt.twoERIs%blockSize('alpha',3).or.myL.gt.twoERIs%blockSize('alpha',4)) &
            call mqc_error('Requested position is out of bounds in Raffenetti 2 ERIs')
          if(twoERIs%hasSpinBlock('aaaa')) element = twoERIs%alpha%at(myI,myJ,myK,myL)
        elseIf(twoERIs%integralType.eq.'raffenetti3') then
          if(myI.gt.twoERIs%blockSize('alpha',1).or.myJ.gt.twoERIs%blockSize('alpha',2).or. &
            myK.gt.twoERIs%blockSize('alpha',3).or.myL.gt.twoERIs%blockSize('alpha',4)) &
            call mqc_error('Requested position is out of bounds in Raffenetti 3 ERIs')
          if(twoERIs%hasSpinBlock('aaaa')) element = twoERIs%alpha%at(myI,myJ,myK,myL)
        elseIf(twoERIs%integralType.eq.'space') then
          if(myI.le.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                element = twoERIs%alpha%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL) - twoERIs%alpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myK = myK - twoERIs%blockSize('alpha',3)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                element = twoERIs%alpha%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL) - twoERIs%alpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myK = myK - twoERIs%blockSize('alpha',3)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aaaa')) element = twoERIs%alpha%at(myI,myJ,myK,myL)
            endIf
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myJ = myJ - twoERIs%blockSize('alpha',2)
            if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aaaa')) element = twoERIs%alpha%at(myI,myJ,myK,myL)
            endIf
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myK = myK - twoERIs%blockSize('alpha',3)
            if(interaction.eq.'exchange') then
              if(twoERIs%hasSpinBlock('aaaa')) element = twoERIs%alpha%at(myI,myL,myK,myJ)
            elseIf(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aaaa')) element = (-1)*twoERIs%alpha%at(myI,myL,myK,myJ)
            endIf
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(interaction.eq.'exchange') then
              if(twoERIs%hasSpinBlock('aaaa')) element = twoERIs%alpha%at(myI,myL,myK,myJ)
            elseIf(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aaaa')) element = (-1)*twoERIs%alpha%at(myI,myL,myK,myJ)
            endIf
          endIf
        elseIf(twoERIs%integralType.eq.'spin') then
          !AAAA
          if(myI.le.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                element = twoERIs%alpha%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL) - twoERIs%alpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          !BBBB
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myK = myK - twoERIs%blockSize('alpha',3)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%hasSpinBlock('bbbb')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%beta%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                element = twoERIs%beta%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                element = twoERIs%beta%at(myI,myJ,myK,myL) - twoERIs%beta%at(myI,myL,myK,myJ)
              endIf
            endIf
          !AABB
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myK = myK - twoERIs%blockSize('alpha',3)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aabb')) element = twoERIs%alphaBeta%at(myI,myJ,myK,myL)
            endIf
          !BBAA
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myJ = myJ - twoERIs%blockSize('alpha',2)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('aabb')) element = twoERIs%alphaBeta%at(myK,myL,myI,myJ)
              endIf
            else
              if(interaction.eq.'coulomb'.or.interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('bbaa')) element = twoERIs%betaAlpha%at(myI,myJ,myK,myL)
              endIf
            endIf
          !ABBA
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myK = myK - twoERIs%blockSize('alpha',3)
            if(interaction.eq.'exchange') then
              if(twoERIs%hasSpinBlock('aabb')) element = twoERIs%alphaBeta%at(myI,myL,myK,myJ)
            elseif(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aabb')) element = (-1)*twoERIs%alphaBeta%at(myI,myL,myK,myJ)
            endIf
          !BAAB
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aabb')) element = twoERIs%alphaBeta%at(myK,myJ,myI,myL)
              elseIf(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('aabb')) element = (-1)*twoERIs%alphaBeta%at(myK,myJ,myI,myL)
              endIf
            else
              if(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('bbaa')) element = twoERIs%betaAlpha%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('bbaa')) element = (-1)*twoERIs%betaAlpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          endIf
        elseIf(twoERIs%integralType.eq.'general') then
          !AAAA
          if(myI.le.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            if(twoERIs%hasSpinBlock('aaaa')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                element = twoERIs%alpha%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%alpha%at(myI,myJ,myK,myL) - twoERIs%alpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          !BBBB
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myK = myK - twoERIs%blockSize('alpha',3)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%hasSpinBlock('bbbb')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%beta%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                element = twoERIs%beta%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%beta%at(myI,myJ,myK,myL) - twoERIs%beta%at(myI,myL,myK,myJ)
              endIf
            endIf
          !AABB
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myK = myK - twoERIs%blockSize('alpha',3)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(interaction.eq.'coulomb') then
              if(twoERIs%hasSpinBlock('aabb')) element = twoERIs%alphaBeta%at(myI,myJ,myK,myL)
            elseif(interaction.eq.'exchange') then
              if(twoERIs%hasSpinBlock('abba')) element = twoERIs%abba%at(myI,myL,myK,myJ)
            elseif(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myI,myJ,myK,myL)
              if(twoERIs%hasSpinBlock('abba')) element = element - twoeris%abba%at(myI,myL,myK,myJ)
            endIf
          !BBAA
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myJ = myJ - twoERIs%blockSize('alpha',2)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('aabb')) element = twoERIs%alphaBeta%at(myK,myL,myI,myJ)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myK,myJ,myI,myL)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myK,myL,myI,myJ)
                if(twoERIs%hasSpinBlock('abba')) element = element - twoeris%abba%at(myK,myJ,myI,myL)
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('bbaa')) element = twoERIs%betaAlpha%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('baab')) element = twoeris%baab%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('bbaa')) element = twoeris%betaAlpha%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('baab')) element = element - twoeris%baab%at(myI,myL,myK,myJ)
              endIf
            endIf
          !ABBA
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myK = myK - twoERIs%blockSize('alpha',3)
            if(interaction.eq.'coulomb') then
              if(twoERIs%hasSpinBlock('abba')) element = twoERIs%abba%at(myI,myJ,myK,myL)
            elseif(interaction.eq.'exchange') then
              if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myI,myL,myK,myJ)
            elseif(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myI,myJ,myK,myL)
              if(twoERIs%hasSpinBlock('aabb')) element = element - twoeris%alphaBeta%at(myI,myL,myK,myJ)
            endIf
          !BAAB
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('abba')) element = twoERIs%abba%at(myK,myL,myI,myJ)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aabb')) element = twoeris%alphaBeta%at(myK,myJ,myI,myL)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('abba')) element = twoeris%abba%at(myK,myL,myI,myJ)
                if(twoERIs%hasSpinBlock('aabb')) element = element - twoeris%alphaBeta%at(myK,myJ,myI,myL)
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('baab')) element = twoERIs%baab%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('bbaa')) element = twoeris%betaAlpha%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('baab')) element = twoeris%baab%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('bbaa')) element = element - twoeris%betaAlpha%at(myI,myL,myK,myJ)
              endIf
            endIf
          !ABAB
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%hasSpinBlock('abab')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%abab%at(myI,myJ,myK,myL)
              elseIf(interaction.eq.'exchange') then
                element = twoeris%abab%at(myI,myL,myK,myJ)
              elseIf(interaction.eq.'doublebar') then
                element = twoERIs%abab%at(myI,myJ,myK,myL) - twoERIs%abab%at(myI,myL,myK,myJ)
              endIf
            endIf
          !BABA
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myK = myK - twoERIs%blockSize('alpha',3)
            if(twoERIs%storageType.eq.'symm') then
              if(twoERIs%hasSpinBlock('abab')) then
                if(interaction.eq.'coulomb') then
                  element = conjg(twoERIs%abab%at(myJ,myI,myL,myK))
                elseIf(interaction.eq.'exchange') then
                  element = conjg(twoeris%abab%at(myL,myI,myJ,myK))
                elseIf(interaction.eq.'doublebar') then
                  element = conjg(twoERIs%abab%at(myJ,myI,myL,myK)) - conjg(twoERIs%abab%at(myL,myI,myJ,myK))
                endIf
              endIf
            else
              if(twoERIs%hasSpinBlock('baba')) then
                if(interaction.eq.'coulomb') then
                  element = twoERIs%baba%at(myI,myJ,myK,myL)
                elseIf(interaction.eq.'exchange') then
                 element = twoeris%baba%at(myI,myL,myK,myJ)
                elseIf(interaction.eq.'doublebar') then
                  element = twoERIs%baba%at(myI,myJ,myK,myL) - twoERIs%baba%at(myI,myL,myK,myJ)
                endIf
              endIf
            endIf
          !AAAB
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myL = myL - twoERIs%blockSize('alpha',4)
            if(interaction.eq.'coulomb') then
              if(twoERIs%hasSpinBlock('aaab')) element = twoERIs%aaab%at(myI,myJ,myK,myL)
            elseif(interaction.eq.'exchange') then
              if(twoERIs%storageType.eq.'symm') then
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myK,myJ,myI,myL)
              else
                if(twoERIs%hasSpinBlock('abaa')) element = twoeris%abaa%at(myI,myL,myK,myJ)
              endIf
            elseif(interaction.eq.'doublebar') then
              if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myJ,myK,myL)
              if(twoERIs%storageType.eq.'symm') then
                if(twoERIs%hasSpinBlock('aaab')) element = element - twoeris%aaab%at(myK,myJ,myI,myL)
              else
                if(twoERIs%hasSpinBlock('abaa')) element = element - twoeris%abaa%at(myI,myL,myK,myJ)
              endIf
            endIf
          !AABA
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myK = myK - twoERIs%blockSize('alpha',3)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoERIs%aaab%at(myJ,myI,myL,myK))
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myL,myI,myJ,myK))
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myJ,myI,myL,myK)) - &
                  conjg(twoeris%aaab%at(myL,myI,myJ,myK))
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('aaba')) element = twoERIs%aaba%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aaba')) element = twoeris%aaba%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('aaba')) element = twoeris%aaba%at(myI,myJ,myK,myL) - &
                  twoeris%aaab%at(myI,myL,myK,myJ)
              endIf
            endIf
          !ABAA
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myJ = myJ - twoERIs%blockSize('alpha',2)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('aaab')) element = twoERIs%aaab%at(myK,myL,myI,myJ)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myK,myL,myI,myJ) - &
                  twoeris%aaab%at(myI,myL,myK,myJ)
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('abaa')) element = twoERIs%abaa%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aaab')) element = twoeris%aaab%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('abaa')) element = twoeris%abaa%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('aaab')) element = element - twoeris%aaab%at(myI,myL,myK,myJ)
              endIf
            endIf
          !BAAA
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoERIs%aaab%at(myL,myK,myJ,myI))
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myJ,myK,myL,myI))
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('aaab')) element = conjg(twoeris%aaab%at(myL,myK,myJ,myI)) - &
                  conjg(twoeris%aaab%at(myJ,myK,myL,myI))
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('baaa')) element = twoERIs%baaa%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('baaa')) element = twoeris%baaa%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('baaa')) element = twoeris%baaa%at(myI,myJ,myK,myL) - &
                  twoeris%baaa%at(myI,myL,myK,myJ)
              endIf
            endIf
          !BBBA
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.le.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myK = myK - twoERIs%blockSize('alpha',3)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoERIs%abbb%at(myL,myK,myJ,myI))
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myL,myI,myJ,myK))
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myL,myK,myJ,myI)) - &
                  conjg(twoeris%abbb%at(myL,myI,myJ,myK))
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('bbba')) element = twoERIs%bbba%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('babb')) element = twoeris%babb%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('bbba')) element = twoeris%bbba%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('babb')) element = element - twoeris%babb%at(myI,myL,myK,myJ)
              endIf
            endIf
          !BBAB
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.le.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('abbb')) element = twoERIs%abbb%at(myK,myL,myI,myJ)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myK,myJ,myI,myL)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('abbb')) element = twoeris%abbb%at(myK,myL,myI,myJ) - &
                  twoeris%abbb%at(myK,myJ,myI,myL)
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('bbab')) element = twoERIs%bbab%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('bbab')) element = twoeris%bbab%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('bbab')) element = twoeris%bbab%at(myI,myJ,myK,myL) - &
                  twoeris%bbab%at(myI,myL,myK,myJ)
              endIf
            endIf
          !BABB
          elseIf(myI.gt.twoERIs%blockSize('alpha',1).and.myJ.le.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myI = myI - twoERIs%blockSize('alpha',1)
            myK = myK - twoERIs%blockSize('alpha',3)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%storageType.eq.'symm') then
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoERIs%abbb%at(myJ,myI,myL,myK))
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myK,myJ,myL,myI))
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('abbb')) element = conjg(twoeris%abbb%at(myJ,myI,myL,myK)) - &
                  conjg(twoeris%abbb%at(myJ,myK,myL,myI))
              endIf
            else
              if(interaction.eq.'coulomb') then
                if(twoERIs%hasSpinBlock('babb')) element = twoERIs%babb%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                if(twoERIs%hasSpinBlock('bbba')) element = twoeris%bbba%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                if(twoERIs%hasSpinBlock('babb')) element = twoeris%babb%at(myI,myJ,myK,myL)
                if(twoERIs%hasSpinBlock('bbba')) element = twoeris%bbba%at(myI,myL,myK,myJ)
              endIf
            endIf
          !ABBB
          elseIf(myI.le.twoERIs%blockSize('alpha',1).and.myJ.gt.twoERIs%blockSize('alpha',2).and. &
            myK.gt.twoERIs%blockSize('alpha',3).and.myL.gt.twoERIs%blockSize('alpha',4)) then
            myJ = myJ - twoERIs%blockSize('alpha',2)
            myK = myK - twoERIs%blockSize('alpha',3)
            myL = myL - twoERIs%blockSize('alpha',4)
            if(twoERIs%hasSpinBlock('abbb')) then
              if(interaction.eq.'coulomb') then
                element = twoERIs%abbb%at(myI,myJ,myK,myL)
              elseif(interaction.eq.'exchange') then
                element = twoeris%abbb%at(myI,myL,myK,myJ)
              elseif(interaction.eq.'doublebar') then
                element = twoERIs%abbb%at(myI,myJ,myK,myL) - twoeris%abbb%at(myI,myL,myK,myJ)
              endIf
            endIf
          endIf
        else
          call mqc_error_A('Integral type not valid in MQC_TwoERIs_At',6,'integralType',twoERIs%integralType)
        endIf
      endIf

      end function mqc_twoERIs_At
!
!
!     PROCEDURE MQC_Integral_At
!
!>    \brief <b> MQC_Integral_At is a function that returns the value of an element 
!>    of a MQC integral variable</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Integral_At is a function that returns the value of (I,J)th element of a 
!>    MQC integral variable as an MQC scalar. If I or J is negative, the (N-I+1)th 
!>    index is selected. Optional argument SpinBlock specifies the spin block of MQC
!>    integrals. If I or J is greater than the dimension of the alpha spin block, the
!>    remainder is used as the index of the beta spin block. The following options 
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
!>    \param[in] Integral
!>    \verbatim
!>        Integral is Class(MQC_SCF_Integral)
!>        The MQC integral to return the value of the (I,J)th 
!>        element.
!>    \endverbatim
!>    
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The first dimension of the element in MQC integral
!>          If I>0 row count is from first index
!>          If I<0 row count is from last index.
!>    \endverbatim
!>    
!>    \param[in] J  
!>    \verbatim
!>        J is Integer(kind=int64)
!>        The second dimension of the element in MQC integral
!>          If J>0 row count is from first index
!>          If J<0 row count is from last index.
!>    \endverbatim
!>    
!>    \param[in,out] SpinBlockIn
!>    \verbatim
!>        SpinBlock is character(len=64),optional
!>        = 'alpha':     alpha spin block
!>        = 'beta':      beta spin block
!>        = 'alphaBeta': alpha-beta spin block
!>        = 'betaAlpha': beta-alpha spin block.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_integral_at(integral,i,j,spinBlockIn) result(element)
!
      implicit none
      class(mqc_scf_integral)::integral
      class(*),intent(in)::i,j
      integer(kind=int64)::myI,myJ
      type(mqc_scalar)::element
      character(len=*),optional,intent(in)::spinBlockIn
      character(len=64)::spinBlock
      real(kind=real64)::zero=0.0d0

      element = zero
      select type (i)  
      type is (integer)
        myI = i
      type is (mqc_scalar)
        myI = i
      class default 
        call mqc_error('index type not defined in MQC_twoERIs_at')
      end select
      select type (j)  
      type is (integer)
        myJ = j
      type is (mqc_scalar)
        myJ = j
      class default 
        call mqc_error('index type not defined in MQC_twoERIs_at')
      end select

      if(myI.gt.integral%blockSize('alpha',1)+integral%blockSize('beta',1)) &
        call mqc_error('Requested I dimension out of bounds in mqc_integral_at')
      if(myJ.gt.integral%blockSize('alpha',2)+integral%blockSize('beta',2)) &
        call mqc_error('Requested J dimension out of bounds in mqc_integral_at')

      if(present(spinBlockIn)) then
        call string_change_case(spinBlockIn,'L',spinBlock)
        if(integral%array_type.eq.'space') then
          if(spinBlock.eq.'alpha'.or.spinBlock.eq.'beta') then
            if(integral%hasAlpha()) element = integral%alpha%at(myI,myJ)
          elseIf(spinBlock.ne.'alphabeta'.and.spinBlock.ne.'betaalpha') then
            call mqc_error_A('SpinBlock type not valid in MQC_Integral_At',6,'SpinBlock',SpinBlock)
          endIf
        elseIf(integral%array_type.eq.'spin') then
          if(spinBlock.eq.'alpha') then
            if(integral%hasAlpha()) element = integral%alpha%at(myI,myJ)
          elseif(spinblock.eq.'beta') then
            if(integral%hasBeta()) element = integral%beta%at(myI,myJ)
          elseIf(spinBlock.ne.'alphabeta'.and.spinBlock.ne.'betaalpha') then
            call mqc_error_A('SpinBlock type not valid in MQC_Integral_At',6,'SpinBlock',SpinBlock)
          endIf
        elseIf(integral%array_type.eq.'general') then
          if(spinBlock.eq.'alpha') then
            if(integral%hasAlpha()) element = integral%alpha%at(myI,myJ)
          elseif(spinblock.eq.'beta') then
            if(integral%hasBeta()) element = integral%beta%at(myI,myJ)
          elseif(spinblock.eq.'alphabeta') then
            if(integral%hasAlphaBeta()) element = integral%alphaBeta%at(myI,myJ)
          elseif(spinblock.eq.'betaalpha') then
            if(integral%hasBetaAlpha()) element = integral%betaAlpha%at(myI,myJ)
          else
            call mqc_error_A('SpinBlock type not valid in MQC_Integral_At',6,'SpinBlock',SpinBlock)
          endIf
        else
          call mqc_error_A('Array type not valid in MQC_Integral_At',6,'Array_Type',integral%array_type)
        endIf
      else 
        if(integral%array_type.eq.'space') then
          if(myI.le.integral%blockSize('alpha',1).and.myJ.le.integral%blockSize('alpha',2)) then
            if(integral%hasAlpha()) element = integral%alpha%at(myI,myJ)
          elseIf(myI.gt.integral%blockSize('alpha',1).and.myJ.gt.integral%blockSize('alpha',2)) then
            myI = myI - integral%blockSize('alpha',1)
            myJ = myJ - integral%blockSize('alpha',2)
            if(integral%hasAlpha()) element = integral%alpha%at(myI,myJ)
          endIf
        elseIf(integral%array_type.eq.'spin') then
          if(myI.le.integral%blockSize('alpha',1).and.myJ.le.integral%blockSize('alpha',2)) then
            if(integral%hasAlpha()) element = integral%alpha%at(myI,myJ)
          elseIf(myI.gt.integral%blockSize('alpha',1).and.myJ.gt.integral%blockSize('alpha',2)) then
            myI = myI - integral%blockSize('alpha',1)
            myJ = myJ - integral%blockSize('alpha',2)
            if(integral%hasBeta()) element = integral%beta%at(myI,myJ)
          endIf
        elseIf(integral%array_type.eq.'general') then
          if(myI.le.integral%blockSize('alpha',1).and.myJ.le.integral%blockSize('alpha',2)) then
            if(integral%hasAlpha()) element = integral%alpha%at(myI,myJ)
          elseIf(myI.le.integral%blockSize('alpha',1).and.myJ.gt.integral%blockSize('alpha',2)) then
            myJ = myJ - integral%blockSize('alpha',2)
            if(integral%hasAlphaBeta()) element = integral%alphaBeta%at(myI,myJ)
          elseIf(myI.gt.integral%blockSize('alpha',1).and.myJ.le.integral%blockSize('alpha',2)) then
            myI = myI - integral%blockSize('alpha',1)
            if(integral%hasBetaAlpha()) element = integral%betaAlpha%at(myI,myJ)
          elseIf(myI.gt.integral%blockSize('alpha',1).and.myJ.gt.integral%blockSize('alpha',2)) then
            myI = myI - integral%blockSize('alpha',1)
            myJ = myJ - integral%blockSize('alpha',2)
            if(integral%hasBeta()) element = integral%beta%at(myI,myJ)
          endIf
        else
          call mqc_error_A('Array type not valid in MQC_Integral_At',6,'Array_Type',integral%array_type)
        endIf
      endIf

      end function mqc_integral_at
!
!
!     PROCEDURE MQC_Eigenvalues_At
!
!>    \brief <b> MQC_Eigenvalues_At is a function that returns the value of an element 
!>    of a MQC eigenvalues variable</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Eigenvalues_At is a function that returns the value of Ith element of a MQC 
!>    eigenvalues variable as an MQC scalar. If I is negative, the (N-I+1)th index is 
!>    selected. Optional argument SpinBlock specifies the spin block of MQC 
!>    eigenvectors. If I is greater than the dimension of the alpha spin block, the
!>    remainder is used as the index of the beta spin block. The following options 
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
!>    \param[in] Eigenvalues
!>    \verbatim
!>        Eigenvalues is Class(MQC_SCF_Eigenvalues)
!>        The MQC eigenvalues to return the value of the Ith 
!>        element.
!>    \endverbatim
!>    
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The dimension of the element in Eigenvalues
!>          If I>0 row count is from first index
!>          If I<0 row count is from last index.
!>    \endverbatim
!>    
!>    \param[in,out] SpinBlock
!>    \verbatim
!>        SpinBlock is character(len=64),optional
!>        = 'alpha':     alpha spin block
!>        = 'beta':      beta spin block
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      function mqc_eigenvalues_at(eigenvalues,i,spinBlock) result(element)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      class(*),intent(in)::i
      integer(kind=int64)::myI
      type(mqc_scalar)::element
      character(len=64),optional::spinBlock
      real(kind=real64)::zero=0.0d0

      select type (i)  
      type is (integer)
        myI = i
      type is (mqc_scalar)
        myI = i
      class default 
        call mqc_error('index type not defined in MQC_twoERIs_at')
      end select

      if(present(spinBlock)) then
        call string_change_case(spinBlock,'L')
        if(eigenvalues%array_type.eq.'space') then
          if(spinBlock.eq.'alpha') then
            element = eigenvalues%alpha%at(myI)
          else
            call mqc_error_A('SpinBlock type not valid in MQC_Eigenvalues_At',6,'SpinBlock',SpinBlock)
          endIf
        elseIf(eigenvalues%array_type.eq.'spin'.or.eigenvalues%array_type.eq.'general') then
          if(spinBlock.eq.'alpha') then
            element = eigenvalues%alpha%at(myI)
          elseif(spinblock.eq.'beta') then
            element = eigenvalues%beta%at(myI)
          else
            call mqc_error_A('SpinBlock type not valid in MQC_Eigenvalues_At',6,'SpinBlock',SpinBlock)
          endIf
        else
          call mqc_error_A('Array type not valid in MQC_Eigenvalues_At',6,'Array_Type',eigenvalues%array_type)
        endIf
      else 
        if(eigenvalues%array_type.eq.'space') then
          if(myI.gt.eigenvalues%blockSize('alpha')) myI = myI - eigenvalues%blockSize('alpha')
          if(myI.lt.(-1)*eigenvalues%blockSize('alpha')) myI = myI + eigenvalues%blockSize('alpha')
          element = eigenvalues%alpha%at(myI)
        elseIf(eigenvalues%array_type.eq.'spin'.or.eigenvalues%array_type.eq.'general') then
          if(myI.gt.eigenvalues%blockSize('alpha')) then
            element = eigenvalues%beta%at(myI)
          elseIf(myI.ge.0) then
            element = eigenvalues%alpha%at(myI)
          elseIf(myI.lt.(-1)*eigenvalues%blockSize('alpha')) then
            myI = myI + eigenvalues%blockSize('alpha')
            element = eigenvalues%beta%at(myI)
          else
            element = eigenvalues%alpha%at(myI)
          endIf
        else
          call mqc_error_A('Array type not valid in MQC_Eigenvalues_At',6,'Array_Type',eigenvalues%array_type)
        endIf
      endIf

      end function mqc_eigenvalues_at
!
!
!     PROCEDURE MQC_SCF_Transformation_Matrix
!
!>    \brief <b> MQC_SCF_Transformation_Matrix is a subroutine that returns the
!>    atomic orbital to orthogonal atomic orbital transformation matrix and the
!>    number of linearly independent atomic orbital basis functions</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_SCF_Transformation_Matrix is a subroutine that returns the atomic orbital 
!>    to orthogonal atomic orbital transformation matrix and the number of linearly 
!>    independent atomic orbital basis functions.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Overlap
!>    \verbatim
!>        Overlap is type(mqc_scf_integral)
!>        The atomic orbital overlap matrix.
!>        element.
!>    \endverbatim
!>    
!>    \param[out] Transform_Matrix  
!>    \verbatim
!>        Transform_Matrix is type(mqc_scf_integral)
!>        The AO -> OAO transformation matrix.
!>    \endverbatim
!>    
!>    \param[out] NBasUse
!>    \verbatim
!>        NBasUse is integer(kind=int64),optional
!>        The number of linearly independent AO basis functions.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2018
!
      subroutine mqc_scf_transformation_matrix(overlap,transform_matrix,nBasUse)
!
      implicit none
      type(mqc_scf_integral),intent(in)::overlap
      type(mqc_scf_integral),intent(out)::transform_matrix
      integer(kind=int64),optional,intent(out)::nBasUse
      type(mqc_matrix)::transform_matrix_tmp,SVecs
      type(mqc_vector)::SVals
      integer(kind=int64)::i
      real(kind=real64)::thresh=1.0d-4
!
      transform_matrix_tmp = overlap%getBlock('alpha')
      call transform_matrix_tmp%diag(SVals,SVecs)
      if(present(nBasUse)) then
        nBasUse = 0
        do i = 1, SVals%size()
          if(MQC_Scalar_Get_Intrinsic_Real(SVals%at(i)).gt.thresh) nBasUse = nBasUse+1
        endDo
      endIf
      call SVals%power(-0.5)
!hph      call mqc_error('Hrant - STOP')
!hph      call SVals%diag(transform_matrix_tmp)
!hph      call SVals%hphdiag(transform_matrix_tmp)
      call MQC_Matrix_DiagMatrix_Put_Vector(SVals,transform_matrix_tmp)
      transform_matrix_tmp = matmul(SVecs,transform_matrix_tmp)
      call mqc_integral_allocate(transform_matrix,'transformation matrix','space',transform_matrix_tmp)
!
      end subroutine mqc_scf_transformation_matrix
!
!
!=====================================================================
!
!     POST-SCF ROUTINES
!
!=====================================================================
!     
!     PROCEDURE GEN_DET_STR 
!
!>    \brief <b> GEN_DET_STR is a subroutine that returns the list of binary strings 
!>    corresponding to occupation number vectors</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    GEN_DET_STR is a subroutine that returns the list of binary strings 
!>    corresponding to occupation number vectors.
!>
!>    The index to the determinant strings follows lexical order in which the 
!>    left most path has index zero (in the case where the right-most bit 
!>    corresponds to orbital 1 and the left-most bit corresponds to orbital 
!>    NBasis - the HF determinant is the smallest integer number). 
!>    Constructing the graph uses the following for the vertex and arc values:
!>    Yo(e,o)=0, Y1(e,o)=x(e,o+1), x(e,o)=x(e+1,o+1)+x(e,o+1) 
!>
!>    The strings are also ordered in reverse lexical order where the right most 
!>    path has index zero (in the case that the left-most bit corresponds to 
!>    orbital 1 and the right-most bit corresponds to orbital NBasis - the HF 
!>    determinant is the largest integer number)
!>    Constructing the graph uses the following for the vertex and arc values:
!>    Yo(e,o)=0, Y1(e,o)=x(e+1,o), x(e,o)=x(e,o-1)+x(e-1,o-1) 
!>
!>    The routine is set up to deal with arbitrary number of active orbitals, up to 
!>    the machine word size. It is unlikley you will need  more than this, as your 
!>    determinant expansion will be really large. 
!>
!>    Optional argument NCore allows inclusion of set of always occupied orbitals
!>    in determinant string. Note that there is no limit on the number of core +
!>    active orbitals. 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64)
!>        The FORTRAN file to print to. 
!>    \endverbatim
!>    
!>    \param[in] IPrint
!>    \verbatim
!>        IPrint is Integer(kind=int64)
!>        The print level for the subroutine. 
!>    \endverbatim
!>    
!>    \param[in] NBasisIn
!>    \verbatim
!>        NBasisIn is Class(*)
!>        The number of basis functions. 
!>    \endverbatim
!>    
!>    \param[in] NAlphaIn
!>    \verbatim
!>        NAlphaIn is Class(*)
!>        The number of alpha electrons.
!>    \endverbatim
!>    
!>    \param[in] NBetaIn
!>    \verbatim
!>        NBetaIn is Class(*)
!>        The number of beta electrons.
!>    \endverbatim
!>    
!>    \param[out] Determinants
!>    \verbatim
!>        Determinants is Type(MQC_Determinant)
!>        The output binary string occupation number vectors.
!>    \endverbatim
!>    
!>    \param[in] NCoreIn
!>    \verbatim
!>        NCoreIn is Class(*),optional
!>        The number of core orbitals to consider always occupied.
!>    \endverbatim
!>
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      Subroutine Gen_Det_Str(IOut,IPrint,NBasisIn,NAlphaIn,NBetaIn, &
          Determinants,NCoreIn)
!
      Implicit None

      Integer(kind=int64)::IOut,IPrint,NBasis,NAlpha,NBeta,NAlpha_Str,NBeta_Str,NBit_Ints,IOrb,&
        NElec,NEMax,NEMin,I,J,K,NCore,Pos
      Integer(kind=int64),Dimension(:,:),Allocatable::Alpha_Strings,Beta_Strings
      Integer(kind=int64),Dimension(:),Allocatable::Alpha_Core,Beta_Core,Temp_String
      Class(*)::NBasisIn,NAlphaIn,NBetaIn
      Class(*),Optional::NCoreIn
      Type(MQC_Determinant)::Determinants
      Type(MQC_LinkedList),Pointer::String_List_1,String_List_2,&
        String_Node_1,String_Node_2
      Integer(kind=int64)::New_Value,SubNum,Scratch
      Integer(kind=int64),Allocatable::Returned_Value
      Logical::Last
!      
 1001 Format(1x,I7,2x,*(B64))
 1050 Format( A )
!
      Select Type (NBasisIn)
      Type Is (Integer)
        NBasis = NBasisIn
      Type Is (MQC_Scalar)
        NBasis = NBasisIn
      Class Default
        Call MQC_Error('NBasisIn type is not integer or MQC scalar')
      End Select

      Select Type (NAlphaIn)
      Type Is (Integer)
        NAlpha = NAlphaIn
      Type Is (MQC_Scalar)
        NAlpha = NAlphaIn
      Class Default
        Call MQC_Error('NAlphaIn type is not integer or MQC scalar')
      End Select

      Select Type (NBetaIn)
      Type Is (Integer)
        NBeta = NBetaIn
      Type Is (MQC_Scalar)
        NBeta = NBetaIn
      Class Default
        Call MQC_Error('NBetaIn type is not integer or MQC scalar')
      End Select

      if(present(NCoreIn)) then
        Select Type (NCoreIn)
        Type Is (Integer)
          NCore = NCoreIn
        Type Is (MQC_Scalar)
          NCore = NCoreIn
        Class Default
          Call MQC_Error('NCoreIn type is not integer or MQC scalar')
        End Select
      EndIf
!
!     Initialize Arrays
!
      String_List_1=>Null()
      String_List_2=>Null()
      String_Node_1=>Null()
      String_Node_2=>Null()
      NAlpha_Str = Bin_Coeff(NBasis,NAlpha)
      NBeta_Str = Bin_Coeff(NBasis,NBeta)
      If(.not.Allocated(Alpha_Strings).and..not.Allocated(Beta_Strings)) then
        if(present(NCoreIn)) then
          Write(iOut,'(A,I4,A)') 'Adding ',NCore,' core orbitals to string'
          NBit_Ints = ((NBasis+NCore)/(Bit_Size(0)-1))+1
        else
          NBit_Ints = (NBasis/(Bit_Size(0)-1))+1
        endIf
        If(((NBasis/(Bit_Size(0)-1))+1).gt.1) then
          Call MQC_Error_I('Determinant generator limited to 1 integer', 6, &
               '(NBasis/(Bit_Size(0)-1))+1', (NBasis/(Bit_Size(0)-1))+1 )
        EndIf
        Allocate(Alpha_Strings(NAlpha_Str,NBit_Ints),Beta_Strings(NBeta_Str,NBit_Ints))
      EndIf
      Alpha_Strings = 0 
      Beta_Strings = 0
!
!    Forming alpha strings
      New_Value = 0
      Call LinkedList_Push(String_List_2,New_Value)
!    Still need to take care of cases where  number of orbitals is greater than one bit
      Do IOrb=1, NBasis
        Call LinkedList_Delete(String_List_1)
        String_Node_2 => String_List_2
        Allocate(String_List_1)
        If(Associated(String_Node_2)) then
          Last = .False.
          Do While(.not.Last)
            Call LinkedList_GetNext(String_Node_2,Last,.True.)
            Call LinkedList_Return_Value(String_Node_2,Returned_Value)
            Call LinkedList_Push(String_List_1,Returned_Value)
          EndDo
        EndIf
        Call LinkedList_Delete(String_List_2)
!      for each link in linked list 1
        String_Node_1 => String_List_1
        Allocate(String_List_2)
        If(Associated(String_Node_1)) then
          Last = .False.
          Do While(.not.Last)
            Call LinkedList_GetNext(String_Node_1,Last,.True.)
            Call LinkedList_Return_Value(String_Node_1,Returned_Value)
!         bit number zero
            New_Value = IBClr(Returned_Value,IOrb-1)
            NElec = PopCnt(New_Value)
            NEMax = NElec + (NBasis-IOrb)
            NEMin = NElec
            If(NEMax.ge.NAlpha.and.NEMin.le.NAlpha) then
              Call LinkedList_Push(String_List_2,New_Value)
            EndIf
!         bit number one
            New_Value = IBSet(Returned_Value,IOrb-1)
            NElec = PopCnt(New_Value)
            NEMax = NElec + (NBasis-IOrb)
            NEMin = NElec
            If(NEMax.ge.NAlpha.and.NEMin.le.NAlpha) then 
              Call LinkedList_Push(String_List_2,New_Value)
            EndIf
          EndDo
        EndIf
      EndDo 

      String_Node_2 => String_List_2
      If(Associated(String_Node_2)) then
        Last = .False.
        I = 1
        Do While(.not.Last)
          Call LinkedList_GetNext(String_Node_2,Last,.True.)
          Call LinkedList_Return_Value(String_Node_2,Returned_Value)
          Alpha_Strings(I,1) = Returned_Value
          I=I+1
        EndDo
      EndIf
!
!    Forming beta strings
      Call LinkedList_Delete(String_List_2)
      Allocate(String_List_2)
      New_Value = 0
      Call LinkedList_Push(String_List_2,New_Value)
!    Still need to take care of cases where  number of orbitals is greater than one bit
      Do IOrb=1, NBasis
        Call LinkedList_Delete(String_List_1)
        String_Node_2 => String_List_2
        Allocate(String_List_1)
        If(Associated(String_Node_2)) then
          Last = .False.
          Do While(.not.Last)
            Call LinkedList_GetNext(String_Node_2,Last,.True.)
            Call LinkedList_Return_Value(String_Node_2,Returned_Value)
            Call LinkedList_Push(String_List_1,Returned_Value)
          EndDo
        EndIf
        Call LinkedList_Delete(String_List_2)
!      for each link in linked list 1
        String_Node_1 => String_List_1
        Allocate(String_List_2)
        If(Associated(String_Node_1)) then
          Last = .False.
          Do While(.not.Last)
            Call LinkedList_GetNext(String_Node_1,Last,.True.)
            Call LinkedList_Return_Value(String_Node_1,Returned_Value)
!         bit number zero
            New_Value = IBClr(Returned_Value,IOrb-1)
            NElec = PopCnt(New_Value)
            NEMax = NElec + (NBasis-IOrb)
            NEMin = NElec
            If(NEMax.ge.NBeta.and.NEMin.le.NBeta) then
              Call LinkedList_Push(String_List_2,New_Value)
            EndIf
!         bit number one
            New_Value = IBSet(Returned_Value,IOrb-1)
            NElec = PopCnt(New_Value)
            NEMax = NElec + (NBasis-IOrb)
            NEMin = NElec
            If(NEMax.ge.NBeta.and.NEMin.le.NBeta) then 
              Call LinkedList_Push(String_List_2,New_Value)
            EndIf
          EndDo
        EndIf
      EndDo 

      String_Node_2 => String_List_2
      If(Associated(String_Node_2)) then
        Last = .False.
        I = 1
        Do While(.not.Last)
          Call LinkedList_GetNext(String_Node_2,Last,.True.)
          Call LinkedList_Return_Value(String_Node_2,Returned_Value)
          Beta_Strings(I,1) = Returned_Value
          I=I+1
        EndDo
      EndIf
!
!     Add core strings if NCore present
!
      if(present(NCoreIn)) then
!        NBit_Ints = ((NBasis+NCore)/(Bit_Size(0)-1))+1
        Allocate(Alpha_Core(NBit_Ints),Beta_Core(NBit_Ints),Temp_String(NBit_Ints))
        Alpha_Core = 0
        K = NBit_Ints + 1
        Do I = 0, NCore-1
          If(Mod(I,Bit_Size(0)).eq.0) K = K-1
          Alpha_Core(K) = IBSet(Alpha_Core(K),I)
        EndDo
        Do I = 1, NAlpha_Str
          Temp_String = 0
          Do J = 1, NBasis
            K = NBit_Ints - (NCore+J)/Bit_Size(0)
            Pos = Mod(NCore+J,Bit_Size(0)) - 1
            Call MvBits(Alpha_Strings(I,1),J-1,1,Temp_String(K),Pos)
          EndDo
          Do K = 1, NBit_Ints
            Alpha_Strings(I,K) = IOr(Temp_String(K),Alpha_Core(K))
          EndDo
        EndDo

        Beta_Core = 0
        K = NBit_Ints + 1
        Do I = 0, NCore-1
          If(Mod(I,Bit_Size(0)).eq.0) K = K-1
          Beta_Core(K) = IBSet(Beta_Core(K),I)
        EndDo
        Do I = 1, NBeta_Str
          Temp_String = 0
          Do J = 1, NBasis
            K = NBit_Ints - (NCore+J)/Bit_Size(0)
            Pos = Mod(NCore+J,Bit_Size(0)) - 1
            Call MvBits(Beta_Strings(I,1),J-1,1,Temp_String(K),Pos)
          EndDo
          Do K = 1, NBit_Ints
            Beta_Strings(I,K) = IOr(Temp_String(K),Beta_Core(K))
          EndDo
        EndDo
      endIf
!
      If(IPrint.ge.4) then
        Write(IOut,1050) "Alpha Strings"
        Do I = 1,NAlpha_Str
          Write(IOut,1001) I,(Alpha_Strings(I,J),J=1,NBit_Ints)
        EndDo
        Write(IOut,1050) "Beta Strings"
        Do I = 1,NBeta_Str
          Write(IOut,1001) I,(Beta_Strings(I,J),J=1,NBit_Ints)
        EndDo
      EndIf

      Determinants%Strings%Alpha = Alpha_Strings
      Determinants%Strings%Beta = Beta_Strings
      Determinants%Order = 'lexical'
      Determinants%NDets = NAlpha_Str*NBeta_Str 
      Determinants%NAlpStr = NAlpha_Str
      Determinants%NBetStr = NBeta_Str 
      subNum = 0
      call mqc_deallocate_vector(determinants%nSubsAlpha)
      do while (.true.)
        call determinants%nSubsAlpha%push(bin_coeff(nAlpha,subNum)*&
          bin_coeff(nBasis-nAlpha,subNum))
        if (determinants%nSubsAlpha%at(subNum+1).ne.0) then
          subNum = subNum + 1
        else
          scratch = determinants%nSubsAlpha%pop()
          exit
        endIf
      endDo
      subNum = 0
      call mqc_deallocate_vector(determinants%nSubsBeta)
      do while (.true.)
        call determinants%nSubsBeta%push(bin_coeff(nBeta,subNum)*&
          bin_coeff(nBasis-nBeta,subNum))
        if (determinants%nSubsBeta%at(subNum+1).ne.0) then
          subNum = subNum + 1
        else
          scratch = determinants%nSubsBeta%pop()
          exit
        endIf
      endDo

      Deallocate(Alpha_Strings,Beta_Strings)

      End Subroutine Gen_Det_Str
!
!=====================================================================
!     
!     PROCEDURE TRCI_DETS_STRING
!
!>    \brief <b> TRCI_DETS_STRING is a subroutine that takes an integer array 
!>    containing possible particle or hole indices for excitations from a reference 
!>    determinant and converts them to alpha and beta strings in a determinant
!>    object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    TRCI_DETS_STRING is a subroutine that takes an integer array containing 
!>    possible particle or hole indices for excitations from a reference determinant 
!>    and converts them to alpha and beta strings in a determinant object</b>
!>
!>    Optional argument NCore allows inclusion of set of always occupied (core) 
!>    orbitals in determinant string. Note that there is no limit on the number of 
!>    core, active, or core + active orbitals. 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64)
!>        The FORTRAN file to print to. 
!>    \endverbatim
!>    
!>    \param[in] IPrint
!>    \verbatim
!>        IPrint is Integer(kind=int64)
!>        The print level for the subroutine. 
!>    \endverbatim
!>    
!>    \param[in] NBasisIn
!>    \verbatim
!>        NBasisIn is Class(*)
!>        The number of basis functions. 
!>    \endverbatim
!>    
!>    \param[in] NAlphaIn
!>    \verbatim
!>        NAlphaIn is Class(*)
!>        The number of alpha electrons.
!>    \endverbatim
!>    
!>    \param[in] NBetaIn
!>    \verbatim
!>        NBetaIn is Class(*)
!>        The number of beta electrons.
!>    \endverbatim
!>    
!>    \param[in] Substitutions
!>    \verbatim
!>        Substitutions is Type(MQC_Vector)
!>        The substituion levels requested.
!>    \endverbatim
!>    
!>    \param[out] Determinants
!>    \verbatim
!>        Determinants is Type(MQC_Determinant)
!>        The output binary string occupation number vectors. 
!>        Note the reference is returned followed by 
!>        substitutions in increasing order.
!>    \endverbatim
!>    
!>    \param[in] NCoreIn
!>    \verbatim
!>        NCoreIn is Class(*),optional
!>        The number of core orbitals to consider always occupied.
!>    \endverbatim
!>
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2021
!
      subroutine TRCI_DETS_STRING(iOut,iPrint,nBasisIn,nAlphaIn,nBetaIn,&
          substitutions,determinants,nCoreIn)
!
      implicit none
      integer(kind=int64)::iOut,iPrint
      class(*),intent(in)::nBasisIn,nAlphaIn,nBetaIn
      class(*),Optional::NCoreIn
      integer(kind=int64)::NBasis,NAlpha,NBeta,NCore,NAlpha_Str,NBeta_Str,NBit_Ints,&
        arrStart,arrEnd,holeStart,partStart,holeInd,partInd
      integer(kind=int64),dimension(:,:),allocatable::Alpha_Strings,Beta_Strings
      integer(kind=int64),dimension(:),intent(in)::substitutions
      integer(kind=int64),dimension(:),allocatable::initstr
      type(mqc_vector)::subs
      type(mqc_matrix)::alpha_arr_hole,alpha_arr_part,beta_arr_hole,beta_arr_part
      type(mqc_vector)::nums_alpha_hole,nums_alpha_part,nums_beta_hole,&
        nums_beta_part,nah,nap,nbh,nbp
      integer::i,j,k,l,m,alpha_elems,beta_elems,arrayind
      type(mqc_determinant),intent(inOut)::determinants
!
 1001 Format(1x,I7,2x,*(B64))
 1050 Format( A )
!
      select type (nBasisIn)
      type is (integer)
        nBasis = nBasisIn
      type is (mqc_scalar)
        nBasis = nBasisIn
      class default
        call mqc_error('NBasisIn type is not integer or MQC scalar')
      end select

      select type (nAlphaIn)
      type is (integer)
        nAlpha = nAlphaIn
      type is (mqc_scalar)
        nAlpha = nAlphaIn
      class default
        call mqc_error('NAlphaIn type is not integer or MQC scalar')
      end select

      select type (nBetaIn)
      type is (integer)
        nBeta = nBetaIn
      type is (mqc_scalar)
        nBeta = nBetaIn
      class default
        call mqc_error('NBetaIn type is not integer or MQC scalar')
      end select

      if(present(nCoreIn)) then
        select type (nCoreIn)
        type is (integer)
          nCore = nCoreIn
        type is (mqc_scalar)
          nCore = nCoreIn
        class default
          call mqc_error('NCoreIn type is not integer or MQC scalar')
        end Select
        write(iOut,'(A,I4,A)') 'Adding ',nCore,' core orbitals to string.'
      else 
        nCore = 0
      endIf
      call trci_dets_arr(nBasis,nAlpha,nBeta,substitutions,alpha_arr_hole, &
        alpha_arr_part,beta_arr_hole,beta_arr_part,nah,nap,nbh,nbp)
!
!     Initialize Arrays
!
      nAlpha_Str = 1 
      nBeta_Str = 1
      do i = 1,size(nah)
        nAlpha_Str = nAlpha_Str + nah%at(i)*nap%at(i)
        nBeta_Str = nBeta_Str + nbh%at(i)*nbp%at(i)
      endDo
      if(allocated(alpha_strings)) deallocate(alpha_strings)
      if(allocated(beta_strings)) deallocate(beta_strings)
      nBit_Ints = ((nBasis+nCore)/(bit_size(0)-1))+1
      allocate(alpha_strings(nAlpha_Str,nBit_Ints),beta_strings(nBeta_str,nBit_ints))
      
      alpha_strings = 0
      beta_strings = 0
      allocate(initstr(nBit_Ints))
      initstr = 0
      j = 0
      do i = 0, nAlpha+nCore-1
        if(mod(i,bit_size(0)).eq.0) j = j+1
        initstr(j) = iBSet(initstr(j),mod(i,bit_size(0)))
      endDo
      do i = 1, nAlpha_Str
        alpha_strings(i,:) = initstr
      endDo
      if(nAlpha.ne.nBeta) then
        initstr = 0
        j = 0
        do i = 0, nBeta+nCore-1
          if(mod(i,bit_size(0)).eq.0) j = j+1
          initstr(j) = iBSet(initstr(j),mod(i,bit_size(0)))
        endDo
        do i = 1, nBeta_Str
          beta_strings(i,:) = initstr
        endDo
      else
        beta_strings = alpha_strings
      endIf

      arrstart = 2
      holestart = 0
      partstart = 0
      do i = 1,size(nah)
        if(i.ne.1) arrstart = arrstart + nah%at(i-1)*nap%at(i-1)
        arrend = arrstart + nah%at(i)*nap%at(i) - 1
        if(i.ne.1) holestart = sum(nah%vat(1,i-1))
        if(i.ne.1) partstart = sum(nap%vat(1,i-1))
        do j = arrstart, arrend 
          holeInd = int((j-arrstart)/nap%at(i))+holestart+1
          partInd = mod((j-arrstart),int(nap%at(i)))+partstart+1
          do k = 1,size(alpha_arr_hole,2)
            if(alpha_arr_hole%at(holeInd,k).ne.0) then
              l = mod(int(alpha_arr_hole%at(holeInd,k))+nCore-1,bit_size(0))+1
              m = (alpha_arr_hole%at(holeInd,k)+nCore)/(bit_size(0)+1)+1
              alpha_strings(j,m) = iBClr(alpha_strings(j,m),l-1)
            endIf
          endDo
          do k = 1,size(alpha_arr_part,2)
            if(alpha_arr_part%at(partInd,k).ne.0) then
              l = mod(int(alpha_arr_part%at(partInd,k))+nCore-1,bit_size(0))+1
              m = (alpha_arr_part%at(partInd,k)+nCore)/(bit_size(0)+1)+1
              alpha_strings(j,m) = iBSet(alpha_strings(j,m),l-1)
            endIf
          endDo
        endDo
      endDo

      arrstart = 2
      holestart = 0
      partstart = 0
      do i = 1,size(nbh)
        if(i.ne.1) arrstart = arrstart + nbh%at(i-1)*nbp%at(i-1)
        arrend = arrstart + nbh%at(i)*nbp%at(i) - 1
        if(i.ne.1) holestart = sum(nbh%vat(1,i-1))
        if(i.ne.1) partstart = sum(nbp%vat(1,i-1))
        do j = arrstart, arrend 
          holeInd = int((j-arrstart)/nbp%at(i))+holestart+1
          partInd = mod((j-arrstart),int(nbp%at(i)))+partstart+1
          do k = 1,size(beta_arr_hole,2)
            if(beta_arr_hole%at(holeInd,k).ne.0) then
              l = mod(int(beta_arr_hole%at(holeInd,k))+nCore-1,bit_size(0))+1
              m = (beta_arr_hole%at(holeInd,k)+nCore)/(bit_size(0)+1)+1
              beta_strings(j,m) = iBClr(beta_strings(j,m),l-1)
            endIf
          endDo
          do k = 1,size(beta_arr_part,2)
            if(beta_arr_part%at(partInd,k).ne.0) then
              l = mod(int(beta_arr_part%at(partInd,k))+nCore-1,bit_size(0))+1
              m = (beta_arr_part%at(partInd,k)+nCore)/(bit_size(0)+1)+1
              beta_strings(j,m) = iBSet(beta_strings(j,m),l-1)
            endIf
          endDo
        endDo
      endDo

      if(iPrint.ge.4) then
        write(iOut,1050) "Alpha Strings"
        do i = 1,nAlpha_Str
          write(iOut,1001) i,(alpha_strings(i,j),j=1,nBit_Ints)
        endDo
        write(iOut,1050) "Beta Strings"
        do i = 1,nBeta_Str
          write(iOut,1001) i,(beta_strings(i,j),j=1,nBit_Ints)
        endDo
      endIf

      determinants%strings%alpha = alpha_strings
      determinants%strings%beta = beta_strings
      determinants%order = 'ci'
      determinants%nDets = nAlpha_str*nBeta_str 
      determinants%nAlpStr = nAlpha_str
      determinants%nBetStr = nBeta_str 
      call determinants%nSubsAlpha%init(maxval(substitutions)+1,0)
      call determinants%nSubsAlpha%put(1,1)
      do i = 1, size(substitutions)
        call determinants%nSubsAlpha%put(nah%at(i)*nap%at(i),substitutions(i)+1)
      endDo
      call determinants%nSubsBeta%init(maxval(substitutions)+1,0)
      call determinants%nSubsBeta%put(1,1)
      do i = 1, size(substitutions)
        call determinants%nSubsBeta%put(nbh%at(i)*nbp%at(i),substitutions(i)+1)
      endDo

      deallocate(alpha_strings,beta_strings)
!
      end subroutine TRCI_DETS_STRING
!
!=====================================================================
!     
!     PROCEDURE TRCI_DETS_ARR 
!
!>    \brief <b> TRCI_DETS_ARR is a subroutine that returns the integer array 
!>    containing possible particle or hole indices for excitations from a reference 
!>    determinant</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    TRCI_DETS_ARR is a subroutine that returns the integer arrary containing 
!>    containing possible particle or hole labels for excitations from a reference 
!>    determinant. Arrays are returned for both alpha and beta electrons for all
!>    the substitution levels requested in order of substitution level (lowest to 
!>    highest). The product of alpha and beta substitutions gives the total 
!>    substitution level, which will be greater than the requested substitution 
!>    level. 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NBasisIn
!>    \verbatim
!>        NBasisIn is Class(*)
!>        The number of basis functions. 
!>    \endverbatim
!>    
!>    \param[in] NAlphaIn
!>    \verbatim
!>        NAlphaIn is Class(*)
!>        The number of alpha electrons.
!>    \endverbatim
!>    
!>    \param[in] NBetaIn
!>    \verbatim
!>        NBetaIn is Class(*)
!>        The number of beta electrons.
!>    \endverbatim
!>    
!>    \param[in] Substitutions
!>    \verbatim
!>        Substitutions is Type(MQC_Vector)
!>        The number of beta electrons.
!>    \endverbatim
!>    
!>    \param[out] Alpha_Arr_Hole
!>    \verbatim
!>        Alpha_Arr_Hole is Type(MQC_Matrix)
!>        The output alpha hole index array.
!>    \endverbatim
!>    
!>    \param[out] Alpha_Arr_Part
!>    \verbatim
!>        Alpha_Arr_Part is Type(MQC_Matrix)
!>        The output alpha particle index array.
!>    \endverbatim
!>    
!>    \param[out] Beta_Arr_Hole
!>    \verbatim
!>        Beta_Arr_Hole is Type(MQC_Matrix)
!>        The output beta hole index array.
!>    \endverbatim
!>    
!>    \param[out] Beta_Arr_Part
!>    \verbatim
!>        Beta_Arr_Part is Type(MQC_Matrix)
!>        The output beta particle index array.
!>    \endverbatim
!>    
!>    \param[out] Nums_Alpha_Hole_Out
!>    \verbatim
!>        Nums_Alpha_Hole_Out is Type(MQC_Vector),Optional
!>        Optional array containing the number of alpha hole 
!>        array elements at each substitution level.
!>    \endverbatim
!>    
!>    \param[out] Nums_Alpha_Part_Out
!>    \verbatim
!>        Nums_Alpha_Part_Out is Type(MQC_Vector),Optional
!>        Optional array containing the number of alpha particle 
!>        array elements at each substitution level.
!>    \endverbatim
!>    
!>    \param[out] Nums_Beta_Hole_Out
!>    \verbatim
!>        Nums_Beta_Hole_Out is Type(MQC_Vector),Optional
!>        Optional array containing the number of beta hole 
!>        array elements at each substitution level.
!>    \endverbatim
!>    
!>    \param[out] Nums_Beta_Part_Out
!>    \verbatim
!>        Nums_Beta_Part_Out is Type(MQC_Vector),Optional
!>        Optional array containing the number of beta particle
!>        array elements at each substitution level.
!>    \endverbatim
!>    
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2021
!
      subroutine TRCI_DETS_ARR(nBasisIn,nAlphaIn,nBetaIn,substitutions,alpha_arr_hole,&
          alpha_arr_part,beta_arr_hole,beta_arr_part,nums_alpha_hole_out,&
          nums_alpha_part_out,nums_beta_hole_out,nums_beta_part_out)
!
      implicit none
      class(*),intent(in)::nBasisIn,nAlphaIn,nBetaIn
      integer(kind=int64)::NBasis,NAlpha,NBeta
      integer(kind=int64),dimension(:),intent(in)::substitutions
      type(mqc_vector)::subs
      type(mqc_matrix),intent(inOut)::alpha_arr_hole,alpha_arr_part,beta_arr_hole,beta_arr_part
      integer(kind=int64),dimension(:,:),allocatable::alpha,beta
      type(mqc_vector)::nums_alpha_hole,nums_alpha_part,nums_beta_hole,nums_beta_part,tmpvec
      type(mqc_vector),intent(inOut),optional::nums_alpha_hole_out,&
        nums_alpha_part_out,nums_beta_hole_out,nums_beta_part_out
      integer::i,alpha_elems,beta_elems,arrayind
      select type (nBasisIn)
      type is (integer)
        nBasis = nBasisIn
      type is (mqc_scalar)
        nBasis = nBasisIn
      class default
        call mqc_error('NBasisIn type is not integer or MQC scalar')
      end select

      select type (nAlphaIn)
      type is (integer)
        nAlpha = nAlphaIn
      type is (mqc_scalar)
        nAlpha = nAlphaIn
      class default
        call mqc_error('NAlphaIn type is not integer or MQC scalar')
      end select

      select type (nBetaIn)
      type is (integer)
        nBeta = nBetaIn
      type is (mqc_scalar)
        nBeta = nBetaIn
      class default
        call mqc_error('NBetaIn type is not integer or MQC scalar')
      end select
      subs = substitutions
      call subs%sort()
      call nums_alpha_hole%init(size(subs),0)
      call nums_alpha_part%init(size(subs),0)
      if(nAlpha.ne.nBeta) then
        call nums_beta_hole%init(size(subs),0)
        call nums_beta_part%init(size(subs),0)
      endIf
      do i = 1, size(subs)
        alpha_elems = bin_coeff(nalpha,subs%at(i))
        call nums_alpha_hole%put(alpha_elems,i)
        if (nAlpha.ne.nBeta) then
          beta_elems = bin_coeff(nbeta,subs%at(i))
          call nums_beta_hole%put(beta_elems,i)
        else
          nums_beta_hole = nums_alpha_hole
        endIf
        alpha_elems = bin_coeff(nbasis-nalpha,subs%at(i))
        call nums_alpha_part%put(alpha_elems,i)
        if (nAlpha.ne.nBeta) then
          beta_elems = bin_coeff(nbasis-nbeta,subs%at(i))
          call nums_beta_part%put(beta_elems,i)
        else
          nums_beta_part = nums_alpha_part
        endIf
      endDo
      alpha_elems = int(sum(nums_alpha_hole))*int(sum(nums_alpha_part))
      beta_elems = int(sum(nums_beta_hole))*int(sum(nums_beta_part))
      call alpha_arr_hole%init(int(sum(nums_alpha_hole)),int(maxval(subs)),0)
      call beta_arr_hole%init(int(sum(nums_beta_hole)),int(maxval(subs)),0)
      call alpha_arr_part%init(int(sum(nums_alpha_part)),int(maxval(subs)),0)
      call beta_arr_part%init(int(sum(nums_beta_part)),int(maxval(subs)),0)
      arrayind=1
      call tmpvec%init(int(maxval(subs)),0)
      do i = 1,size(subs)
        call build_trci_ph_list(alpha_arr_hole,arrayind,1,nAlpha,1,int(subs%at(i)),tmpvec)
      endDo
      if(nAlpha.ne.nBeta) then
        arrayind=1
        call tmpvec%init(int(maxval(subs)),0)
        do i = 1,size(subs)
          call build_trci_ph_list(beta_arr_hole,arrayind,1,nBeta,1,int(subs%at(i)),tmpvec)
        endDo
      else
        beta_arr_hole = alpha_arr_hole
      endIf
      arrayind=1
      call tmpvec%init(int(maxval(subs)),0)
      do i = 1,size(subs)
        call build_trci_ph_list(alpha_arr_part,arrayind,NAlpha+1,nBasis,1,int(subs%at(i)),tmpvec)
      endDo
      if(nAlpha.ne.nBeta) then
        arrayind=1
        call tmpvec%init(int(maxval(subs)),0)
        do i = 1,size(subs)
          call build_trci_ph_list(beta_arr_part,arrayind,NBeta+1,nBasis,1,int(subs%at(i)),tmpvec)
        endDo
      else
        beta_arr_part = alpha_arr_part
      endIf
!
      if(present(nums_alpha_hole_out)) nums_alpha_hole_out = nums_alpha_hole
      if(present(nums_alpha_part_out)) nums_alpha_part_out = nums_alpha_part
      if(present(nums_beta_hole_out)) nums_beta_hole_out = nums_beta_hole
      if(present(nums_beta_part_out)) nums_beta_part_out = nums_beta_part
!
      end subroutine TRCI_DETS_ARR 
!
!=====================================================================
!     
!     PROCEDURE Build_TrCi_PH_List 
!
!>    \brief <b> Build_TrCi_PH_List is a recursive subroutine that builds an array 
!>    containing the particle-hole indices for a requested substitution level from a 
!>    reference determinant</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    Build_TrCi_PH_List is a recursive subroutine that builds an array containing 
!>    the particle-hole indices for a requested substitution level from a reference 
!>    determinant. If the looptotal is less than the length of store, the array is 
!>    right-padded with zeros.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] Array
!>    \verbatim
!>        Array is type(mqc_matrix)
!>        The array to update with particle/hole indices. 
!>    \endverbatim
!>    
!>    \param[in,out] ArrayInd 
!>    \verbatim
!>        ArrayInd is Integer
!>        The row of array which is the number of particle/home
!>        states.
!>    \endverbatim
!>    
!>    \param[in] IndStart
!>    \verbatim
!>        IndStart is Integer
!>        The initial index of the recursive loop.
!>    \endverbatim
!>    
!>    \param[in] IndEnd
!>    \verbatim
!>        IndEnd is Integer 
!>        The final index of the recursive loop.
!>    \endverbatim
!>    
!>    \param[in] LoopNum
!>    \verbatim
!>        LoopNum is Integer
!>        The current nested loop number.
!>    \endverbatim
!>    
!>    \param[in] LoopTotal
!>    \verbatim
!>        LoopTotal is Integer
!>        The total number of nested loops.
!>    \endverbatim
!>    
!>    \param[in,out] Store
!>    \verbatim
!>        Store is Type(MQC_Vector)
!>        Temporary vector that stores indices before writing to
!>        Array.
!>    \endverbatim
!>    
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2021
!
      recursive subroutine build_trci_ph_list(array,arrayind,indstart,indend,loopnum,looptotal,store)
!
      implicit none
      type(mqc_matrix),intent(inOut)::array
      integer,intent(in)::looptotal,loopnum,indstart,indend
      integer::i
      integer,intent(inOut)::arrayind
      type(mqc_vector)::store
!
      do i = indstart,indend
        if(arrayind.gt.size(array,1)) return
        call store%put(i,loopnum)
        call array%vput(store,[arrayind],[0])
        if(loopnum.lt.looptotal) then
          call build_trci_ph_list(array,arrayind,i+1,indend,loopnum+1,looptotal,store)
        else
          arrayind = arrayind + 1
        endIf
      endDo

      end subroutine build_trci_ph_list
!
!=====================================================================
!     
!     PROCEDURE MQC_DetString_to_OccArray
!
!>    \brief <b> MQC_DetString_to_OccArray takes as an argument a bit string
!>    and returns an integer array with locations of occupied orbitals</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_DetString_to_OccArray is a function that takes as an argument a 
!>    bit string array and returns an integer array with locations where 
!>    bits are set equal to one. The value of the routine is that it can 
!>    be used in conjuction with mqc_integral_output_orbitals to return
!>    an mqc_scf_integral object containing the set of occupied orbitals.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] String
!>    \verbatim
!>        String is Type(MQC_Vector)
!>        The bit string to be analyzed.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2024
!
      function mqc_detString_to_occArray(string) result(arrayOut)

      implicit none
      type(mqc_vector),intent(in)::string
      integer(kind=int64),dimension(:),allocatable::arrayOut

      integer::nBitInts,temp_int,i,j,nBasis
      type(mqc_vector)::tmpArray

      nBitInts = size(string)
      nBasis = nBitInts*Bit_Size(0)

      do i = 1, nBitInts
        temp_int = string%at(i)
        do j = 0, Bit_Size(0)-1
          if(btest(temp_int,j)) call tmpArray%push((i-1)*(Bit_Size(0)-1)+(j+1))
        endDo
      endDo
      allocate(arrayOut(size(tmpArray)))
      arrayOut = tmpArray

      end function mqc_detstring_to_occArray
!
!=====================================================================
!     
!     PROCEDURE SLATER_CONDON
!
!>    \brief <b> SLATER_CONDON is a function that returns an operator matrix 
!>    element using Slater-Condon rules</b> 
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    SLATER_CONDON is a function that returns an operator matrix element value
!>    using Slater-Condon rules for a given alpha and beta string combination.
!>    If UHF flag is set to true then Slater Condon rules are returned accounting for 
!>    different spatial occupation of alpha and beta orbitals. Only one set of 
!>    integrals need be present so that just the one or two particle terms, or both
!>    one and two particle terms, can be returned.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64)
!>        The FORTRAN file to print to. 
!>    \endverbatim
!>    
!>    \param[in] IPrint
!>    \verbatim
!>        IPrint is Integer(kind=int64)
!>        The print level for the subroutine. 
!>    \endverbatim
!>    
!>    \param[in] NBasisIn
!>    \verbatim
!>        NBasisIn is Type(MQC_Scalar)
!>        The number of basis functions. 
!>    \endverbatim
!>    
!>    \param[in] Alpha_String_1
!>    \verbatim
!>        Alpha_String_1 is Integer(kind=int64),Dimension(:),Allocatable
!>        The left alpha string.
!>    \endverbatim
!>    
!>    \param[in] Beta_String_1
!>    \verbatim
!>        Beta_String_1 is Integer(kind=int64),Dimension(:),Allocatable
!>        The left beta string.
!>    \endverbatim
!>    
!>    \param[in] Alpha_String_2
!>    \verbatim
!>        Alpha_String_2 is Integer(kind=int64),Dimension(:),Allocatable
!>        The right alpha string.
!>    \endverbatim
!>    
!>    \param[in] Beta_String_2
!>    \verbatim
!>        Beta_String_2 is Integer(kind=int64),Dimension(:),Allocatable
!>        The right beta string.
!>    \endverbatim
!>    
!>    \param[in] OnePartInts
!>    \verbatim
!>        OnePartInts is Type(MQC_SCF_Integral),Optional
!>        The MO basis one-particle integral matrix.
!>    \endverbatim
!>    
!>    \param[in] TwoPartInts
!>    \verbatim
!>        TwoPartInts is Type(MQC_TwoERIs),Optional 
!>        The MO basis two-particle integrals.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      Function Slater_Condon(IOut,IPrint,NBasisIn,Alpha_String_1,Beta_String_1, &
       Alpha_String_2,Beta_String_2,OnePartInts,TwoPartInts) &
       Result(MatEl)
!
!
!     Variable Declarations...
!
      Implicit None
      Integer(kind=int64),Intent(In)::IOut,IPrint
      Type(MQC_Scalar),Intent(In)::NBasisIn
      Type(MQC_SCF_Integral),Optional,Intent(In)::OnePartInts
      Integer(kind=int64),Dimension(:),Allocatable,Intent(In)::Alpha_String_1,Alpha_String_2,&
        Beta_String_1,Beta_String_2
      Type(MQC_TwoERIs),Optional,Intent(In)::TwoPartInts
      Type(MQC_Scalar)::MatEl,Sgn
      Integer(kind=int64)::NBasis,IPos,JPos,IDiff,Det_Diff,ISgn,NBit_Ints,I,J,K, &
        II,JJ,Alpha_Diff_Cnt,Beta_Diff_Cnt
      real(kind=real64)::ERI1,ERI2,Zero=0.0
      Integer(kind=int64),Dimension(4)::Orbs,Spin,Det
      Integer(kind=int64),Dimension(:),Allocatable::Alpha_Diff,Beta_Diff
      
      If(Present(OnePartInts)) then
        If(OnePartInts%type().ne.'space'.and.OnePartInts%type().ne.'spin') &
          call mqc_error_A('Slater_Condon only implemented for &
          & spin or space one-particle integrals',6,'OnePartInts%type()',OnePartInts%type())
      EndIf
      If(Present(TwoPartInts)) then
        If(TwoPartInts%type().ne.'space'.and.TwoPartInts%type().ne.'spin') call mqc_error_A('Slater_Condon only implemented for &
          & spin or space two-particle integrals',6,'TwoPartInts%type()',TwoPartInts%type())
      endIf
      If(.not.present(OnePartInts).and..not.present(TwoPartInts)) call mqc_error('No integrals given to Slater_Condon')

 1050 Format( A )
      NBasis = NBasisIn
      NBit_Ints = (NBasis/(Bit_Size(0)-1))+1 
!
!     Initialize Arrays
!
!      Write(IOut,*) '---alpha 1:---'
!      Write(IOut,'(B64)') Alpha_String_1
!      Write(IOut,*) '---alpha 2:---'
!      Write(IOut,'(B64)') Alpha_String_2
!      Write(IOut,*) '-------------'
!      Write(IOut,*) 'alpha XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Alpha_String_1,Alpha_String_2)   
!      Write(IOut,*) 'NDifa=', PopCnt(IEOR(Alpha_String_1,Alpha_String_2))
!      Write(IOut,*)
!      Write(IOut,*) '--beta 1:--'
!      Write(IOut,'(B64)') Beta_String_1
!      Write(IOut,*) '---beta 2:---'
!      Write(IOut,'(B64)') Beta_String_2
!      Write(IOut,*) '-------------'
!      Write(IOut,*) 'beta XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Beta_String_1,Beta_String_2)   
!      Write(IOut,*) 'NDifb=', PopCnt(IEOR(Beta_String_1,Beta_String_2))
!      Write(IOut,*)

      Allocate(Alpha_Diff(NBit_Ints),Beta_Diff(NBit_Ints))
      Det_Diff = 0
      Alpha_Diff_Cnt = 0
      Beta_Diff_Cnt = 0
      Do I = 1,NBit_Ints 
        Alpha_Diff(I) = IEOR(Alpha_String_1(I),Alpha_String_2(I))
!        Write(IOut,*) 'Alpha Diff',I,':'
!        Write(IOut,'(B64)') Alpha_Diff(I)
!        Write(IOut,*) '-------------'
        Alpha_Diff_Cnt = Alpha_Diff_Cnt + PopCnt(Alpha_Diff(I)) 
        Beta_Diff(I) = IEOR(Beta_String_1(I),Beta_String_2(I))
!        Write(IOut,*) 'Beta Diff',I,':'
!        Write(IOut,'(B64)') Beta_Diff(I)
!        Write(IOut,*) '-------------'
        Beta_Diff_Cnt = Beta_Diff_Cnt + PopCnt(Beta_Diff(I))
      EndDo
!      Write(IOut,*)'Alpha_Diff_Cnt:',Alpha_Diff_Cnt,'Beta_Diff_Cnt:',Beta_Diff_Cnt
      Det_Diff = Alpha_Diff_Cnt/2 + Beta_Diff_Cnt/2

      If(Mod(Alpha_Diff_Cnt,2).ne.0.or.Mod(Beta_Diff_Cnt,2).ne.0) & 
        Call MQC_Error_I('Slater_Condon has been handed spin non-conserving determinants', 6, &
        'Mod(Alpha_Diff_Cnt,2)', Mod(Alpha_Diff_Cnt,2), &
        'Mod(Beta_Diff_Cnt,2)', Mod(Beta_Diff_Cnt,2) )

!      Write(IOut,*) "Det_Diff:",Det_Diff
      Select Case (Det_Diff)
!
      Case(3:)
        MatEl = Zero 
        Return
!
      Case(2)
!     If we are just doing a one-particle operator, case 2 is zero
        If(.not.present(TwoPartInts)) then
          MatEl = Zero
          Return
        EndIf
!     I am going to comment the first logical block - the rest are different
!     permutations
!     First we need to determine the relevant orbital, spin and determinant
        IDiff = 1
        Do IPos = 0, NBasis-1
          I = IPos/Bit_Size(0)+1
          J = Mod(IPos,Bit_Size(0)) 
!          Write(IOut,*) 'Orbs:',IPos,' I:',I,' J:',J
          If(BTest(Alpha_Diff(I),J)) then
            Orbs(IDiff) = IPos+1
            Spin(IDiff) = 0
            If(BTest(Alpha_String_1(I),J)) then
              Det(IDiff) = 1
            Else 
              Det(IDiff) = 2
            EndIf
            IDiff = IDiff + 1
          EndIf
          If(BTest(Beta_Diff(I),J)) then
            Orbs(IDiff) = IPos+1
            Spin(IDiff) = 1
            If(BTest(Beta_String_1(I),J)) then
              Det(IDiff) = 1
            Else 
              Det(IDiff) = 2
            EndIf
            IDiff = IDiff + 1
          EndIf
        EndDo
!        Call MQC_Print(IOut,Orbs,'Orbs')
!        Call MQC_Print(IOut,Spin,'Spin')
!        Call MQC_Print(IOut,Det,'Det')

        If(Det(1).eq.Det(2).and.Det(3).eq.Det(4)) then
!          Write(IOut,*) 'Entering case 1'
!          
!         This section computes the Coulomb contribution to the matrix
!         element
          If(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
!            
!           This section detemines number of permutations of creation
!           operators and sets sign accordingly
!           It got quite a bit longer after generalizing max number of
!           orbitals
!           Probably want to make some functions for these then later
            ISgn = 0
            If(Spin(1).eq.0.and.Spin(2).eq.0) then
              If(Orbs(4)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(3)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
              If(Orbs(4)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(3)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
              If(Orbs(4)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(3)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
              If(Orbs(4)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(3)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            EndIf
!            
!           This section computes the value of the matrix element
!
            If(TwoPartInts%type().eq.'space'.or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(3),Orbs(2),Orbs(4),'alpha')
            ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(3),Orbs(2),Orbs(4),'alphaBeta')
            ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(3),Orbs(2),Orbs(4),'betaAlpha')
            ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(3),Orbs(2),Orbs(4),'beta')
            EndIf
!            
!           This section changes sign accounting for the ordering of alpha and
!           beta strings
            If(Spin(1).gt.Spin(2)) ERI1 = -ERI1
            If(Spin(3).gt.Spin(4)) ERI1 = -ERI1
            If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!            Write(IOut,*) 'Setting ERI1 to',ERI1
          Else
!            
!           Sets Coulomb contribution to zero if 2ERI is necessarily so by
!           spin order in the integral
!            Write(IOut,*) 'Setting ERI1 to zero'
            ERI1 = Zero
          EndIf
          If(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
            ISgn = 0
            If(Spin(1).eq.0.and.Spin(2).eq.0) then
              If(Orbs(3)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(4)-Orbs(1)-1.gt.0) then 
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
              If(Orbs(3)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(4)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
              If(Orbs(3)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(4)-Orbs(1)-1.gt.0) then 
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
              If(Orbs(3)-Orbs(2)-1.gt.0) then 
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(4)-Orbs(1)-1.gt.0) then 
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            EndIf
!
            If(TwoPartInts%type().eq.'space'.or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(4),Orbs(2),Orbs(3),'alpha')
            ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(4),Orbs(2),Orbs(3),'alphaBeta')
            ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(4),Orbs(2),Orbs(3),'betaAlpha')
            ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(4),Orbs(2),Orbs(3),'beta')
            EndIf
!
            If(Spin(1).gt.Spin(2)) ERI2 = -ERI2
            If(Spin(3).gt.Spin(4)) ERI2 = -ERI2
            If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!            Write(IOut,*) 'Setting ERI2 to',ERI2
          Else 
!            Write(IOut,*) 'Setting ERI2 to zero'
            ERI2 = Zero
          EndIf
!
        ElseIf(Det(1).eq.Det(3).and.Det(2).eq.Det(4)) then
!          Write(IOut,*) 'Entering case 2'
          If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
            ISgn = 0
            If(Spin(1).eq.0.and.Spin(3).eq.0) then
              If(Orbs(4)-Orbs(3)-1.gt.0) then 
                I  = (Orbs(3)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(2)-Orbs(1)-1.gt.0) then 
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(2)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
              If(Orbs(4)-Orbs(3)-1.gt.0) then
                I  = (Orbs(3)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(2)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(2)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
              If(Orbs(4)-Orbs(3)-1.gt.0) then
                I  = (Orbs(3)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(2)-Orbs(1)-1.gt.0) then 
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(2)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
              If(Orbs(4)-Orbs(3)-1.gt.0) then
                I  = (Orbs(3)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(2)-Orbs(1)-1.gt.0) then 
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(2)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            EndIf
            If(TwoPartInts%type().eq.'space'.or.(Spin(1).eq.0.and.Spin(3).eq.0)) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(2),Orbs(3),Orbs(4),'alpha')
            ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(2),Orbs(3),Orbs(4),'alphaBeta')
            ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(2),Orbs(3),Orbs(4),'betaAlpha')
            ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(2),Orbs(3),Orbs(4),'beta')
            EndIf
            If(Spin(1).gt.Spin(3)) ERI1 = -ERI1
            If(Spin(2).gt.Spin(4)) ERI1 = -ERI1
            If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!            Write(IOut,*) 'Setting ERI1 to',ERI1
          Else
!            Write(IOut,*) 'Setting ERI1 to zero'
            ERI1 = Zero
          EndIf
          If(Spin(1).eq.Spin(4).and.Spin(3).eq.Spin(2)) then
            ISgn = 0
            If(Spin(1).eq.0.and.Spin(3).eq.0) then
              If(Orbs(3)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(4)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
              If(Orbs(3)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(4)-Orbs(1)-1.gt.0) then 
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
              If(Orbs(3)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(4)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
              If(Orbs(3)-Orbs(2)-1.gt.0) then 
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(4)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            EndIf
            If(TwoPartInts%type().eq.'space'.or.(Spin(1).eq.0.and.Spin(3).eq.0)) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(4),Orbs(3),Orbs(2),'alpha')
            ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(4),Orbs(3),Orbs(2),'alphaBeta')
            ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(4),Orbs(3),Orbs(2),'betaAlpha')
            ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(4),Orbs(3),Orbs(2),'beta')
            EndIf
            If(Spin(1).gt.Spin(3)) ERI2 = -ERI2
            If(Spin(2).gt.Spin(4)) ERI2 = -ERI2
            If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!            Write(IOut,*) 'Setting ERI2 to',ERI2
          Else 
!            Write(IOut,*) 'Setting ERI2 to zero'
            ERI2 = Zero
          EndIf
!
        ElseIf(Det(1).eq.Det(4).and.Det(2).eq.Det(3)) then
!          Write(IOut,*) 'Entering case 3'
          If(Spin(1).eq.Spin(2).and.Spin(4).eq.Spin(3)) then
            ISgn = 0
            If(Spin(1).eq.0.and.Spin(4).eq.0) then
              If(Orbs(4)-Orbs(3)-1.gt.0) then
                I  = (Orbs(3)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(2)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(2)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
              If(Orbs(4)-Orbs(3)-1.gt.0) then
                I  = (Orbs(3)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(2)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(2)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
              If(Orbs(4)-Orbs(3)-1.gt.0) then
                I  = (Orbs(3)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(2)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(2)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
              If(Orbs(4)-Orbs(3)-1.gt.0) then
                I  = (Orbs(3)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(2)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(2)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
!                  Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            EndIf
            If(TwoPartInts%type().eq.'space'.or.(Spin(1).eq.0.and.Spin(4).eq.0)) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(2),Orbs(4),Orbs(3),'alpha')
            ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(2),Orbs(4),Orbs(3),'alphaBeta')
            ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(2),Orbs(4),Orbs(3),'betaAlpha')
            ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
              ERI1 = TwoPartInts%at(Orbs(1),Orbs(2),Orbs(4),Orbs(3),'beta')
            EndIf
            If(Spin(1).gt.Spin(4)) ERI1 = -ERI1
            If(Spin(2).gt.Spin(3)) ERI1 = -ERI1
            If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!            Write(IOut,*) 'Setting ERI1 to',ERI1
          Else
!            Write(IOut,*) 'Setting ERI1 to zero'
            ERI1 = Zero
          EndIf
          If(Spin(1).eq.Spin(3).and.Spin(4).eq.Spin(2)) then
            ISgn = 0
            If(Spin(1).eq.0.and.Spin(4).eq.0) then
              If(Orbs(4)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(3)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
              If(Orbs(4)-Orbs(2)-1.gt.0) then 
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(3)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
              If(Orbs(4)-Orbs(2)-1.gt.0) then 
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
                  ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(3)-Orbs(1)-1.gt.0) then 
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
              If(Orbs(4)-Orbs(2)-1.gt.0) then
                I  = (Orbs(2)-1)/Bit_Size(0)+1
                J  = (Orbs(4)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                  EndIf
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
              If(Orbs(3)-Orbs(1)-1.gt.0) then
                I  = (Orbs(1)-1)/Bit_Size(0)+1
                J  = (Orbs(3)-1)/Bit_Size(0)+1
                Do K = I, J, -1
                  If(I.ne.K) then 
                    II = 1
                  Else 
                    II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                  EndIf
                  If(J.ne.K) then 
                    JJ = Bit_Size(0)
                  Else
                    JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                  EndIf
                  ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                EndDo
              EndIf
            EndIf
            If(TwoPartInts%type().eq.'space'.or.(Spin(1).eq.0.and.Spin(4).eq.0)) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(3),Orbs(4),Orbs(2),'alpha')
            ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(3),Orbs(4),Orbs(2),'alphaBeta')
            ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(3),Orbs(4),Orbs(2),'betaAlpha')
            ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
              ERI2 = TwoPartInts%at(Orbs(1),Orbs(3),Orbs(4),Orbs(2),'beta')
            EndIf
            If(Spin(1).gt.Spin(4)) ERI2 = -ERI2
            If(Spin(2).gt.Spin(3)) ERI2 = -ERI2
            If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!            Write(IOut,*) 'Setting ERI2 to',ERI2
          Else 
!            Write(IOut,*) 'Setting ERI2 to zero'
            ERI2 = Zero
          EndIf
        Else
          Write(IOut,1050) 'Slater Condon has been handed confusing determinant info'
        EndIf
!
        MatEl = ERI1 - ERI2
!        Write(IOut,*) 'Slater Condon Final:',float(MatEl)
!        Deallocate(Orbs,Spin,Det)
        Return
!
      Case(1)
        IDiff = 1
        Do IPos = 0, NBasis-1
          I = IPos/Bit_Size(0)+1
          J = Mod(IPos,Bit_Size(0)) 
          If(BTest(Alpha_Diff(I),J)) then
            Orbs(IDiff) = IPos+1
            Spin(IDiff) = 0
            IDiff = IDiff + 1
          EndIf
          If(BTest(Beta_Diff(I),J)) then
            Orbs(IDiff) = IPos+1
            Spin(IDiff) = 1
            IDiff = IDiff + 1
          EndIf
        EndDo
!        Write(IOut,*)'Orb 1:',Orbs(1),' Orb 2:',Orbs(2)
!        Write(IOut,*)'Spin 1:',Spin(1),' Spin 2:',Spin(2)
!        Write(IOut,*)'Det 1:',Det(1),' Det 2:',Det(2)
!
        ISgn = 0
        If(Spin(1).eq.0.and.Spin(2).eq.0) then
          If(Orbs(2)-Orbs(1)-1.gt.0) then 
            I  = (Orbs(1)-1)/Bit_Size(0)+1
            J  = (Orbs(2)-1)/Bit_Size(0)+1
            Do K = I, J, -1
              If(I.ne.K) then 
                II = 1
              Else 
                II = Mod((Orbs(1)-1),Bit_Size(0)) + 1 
              EndIf
              If(J.ne.K) then 
                JJ = Bit_Size(0) 
              Else
                JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
              EndIf
              ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
            EndDo
          EndIf
        ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
          If(Orbs(2)-Orbs(1)-1.gt.0) then
            I  = (Orbs(1)-1)/Bit_Size(0)+1
            J  = (Orbs(2)-1)/Bit_Size(0)+1
            Do K = I, J, -1
              If(I.ne.K) then 
                II = 1
              Else 
                II = Mod((Orbs(1)-1),Bit_Size(0)) + 1 
              EndIf
              If(J.ne.K) then 
                JJ = Bit_Size(0) 
              Else
                JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
              EndIf
              ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
            EndDo
          EndIf
        EndIf
        Sgn = (-1)**ISgn
!        Write(IOut,*)'ISgn:',ISgn
!
        MatEl = Zero
        If(Present(TwoPartInts)) then
          Do IPos = 0, NBasis-1
            I = IPos/Bit_Size(0)+1
            J = Mod(IPos,Bit_Size(0)) 
!            Write(IOut,*) 'Alpha IPos:', IPos
            If(BTest(Alpha_String_1(I),J)) then
              If(TwoPartInts%type().eq.'space'.and.Spin(1).eq.Spin(2)) then
!                Call MQC_Print(Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'alpha'),iOut, &
!                  'ALPHA: UHF False and 1 and 2 same spin, adding:')
                MatEl = MatEl + Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'alpha')
                If(Spin(1).eq.0.and.Spin(2).eq.0) then
!                  Call MQC_Print(Sgn*TwoPartInts%at(Orbs(1),IPos+1,IPos+1,Orbs(2),'alpha'),iOut, &
!                    'ALPHA: UHF False and 1 and 2 both alpha, subtracting:')
                  MatEl = MatEl - Sgn*TwoPartInts%at(Orbs(1),IPos+1,IPos+1,Orbs(2),'alpha')
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.0) then
!                Call MQC_Print(Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'alpha') - &
!                  Sgn*TwoPartInts%at(Orbs(1),IPos+1,IPos+1,Orbs(2),'alpha'), iOut, &
!                  'ALPHA: UHF True and 1 and 2 both alpha, adding:')
                MatEl = MatEl + Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'alpha') - &
                  Sgn*TwoPartInts%at(Orbs(1),IPos+1,IPos+1,Orbs(2),'alpha')
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!                Call MQC_Print(Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'betaAlpha'),iOut, &
!                  'ALPHA: UHF True and 1 and 2 both beta, adding:')
                MatEl = MatEl + Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'betaAlpha')
              EndIf
            EndIf
          EndDo
!
          Do IPos = 0, NBasis-1
            I = IPos/Bit_Size(0)+1
            J = Mod(IPos,Bit_Size(0))
            If(BTest(Beta_String_1(I),J)) then
              If(TwoPartInts%type().eq.'space'.and.Spin(1).eq.Spin(2)) then
!                Call MQC_Print(Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'alpha'),iOut, &
!                  'BETA: UHF False and 1 and 2 same spin, adding:')
                MatEl = MatEl + Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'alpha')
                If(Spin(1).eq.1.and.Spin(2).eq.1) then
!                  Call MQC_Print(Sgn*TwoPartInts%at(Orbs(1),IPos+1,IPos+1,Orbs(2),'alpha'),iOut, &
!                    'BETA: UHF False and 1 and 2 both beta, subtracting:')
                  MatEl = MatEl - Sgn*TwoPartInts%at(Orbs(1),IPos+1,IPos+1,Orbs(2),'alpha')
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.0) then
!                Call MQC_Print(Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'alphaBeta'),iOut, &
!                  'BETA: UHF True and 1 and 2 both alpha, adding:')
                MatEl = MatEl + Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'alphaBeta') 
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!                Call MQC_Print(Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'beta') - &
!                  Sgn*TwoPartInts%at(Orbs(1),IPos+1,IPos+1,Orbs(2),'beta'),iOut, &
!                  'BETA: UHF True and 1 and 2 both beta, adding:')
                MatEl = MatEl + Sgn*TwoPartInts%at(Orbs(1),Orbs(2),IPos+1,IPos+1,'beta') - &
                  Sgn*TwoPartInts%at(Orbs(1),IPos+1,IPos+1,Orbs(2),'beta')
              EndIf
            EndIf
          EndDo
        EndIf
!
        If(Present(OnePartInts)) then
          If(OnePartInts%type().eq.'space'.or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
!            Call MQC_Print(Sgn*OnePartInts%Alpha%at(Orbs(1),Orbs(2)),iOut,'1 and 2 both alpha, adding core:')
            MatEl = MatEl + Sgn*OnePartInts%Alpha%at(Orbs(1),Orbs(2)) 
          ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!            Call MQC_Print(Sgn*OnePartInts%Beta%at(Orbs(1),Orbs(2)),iOut,'1 and 2 both beta, adding core:')
            MatEl = MatEl + Sgn*OnePartInts%Beta%at(Orbs(1),Orbs(2))
          EndIf
        EndIf
        
!        Call MQC_Print(MatEl,IOut,'CI Matrix Element')

        Return
!
      Case(0)
        MatEl = Zero
        If(Present(OnePartInts)) then
          Do IPos = 0, NBasis-1
            I = IPos/Bit_Size(0)+1
            J = Mod(IPos,Bit_Size(0))
            If(BTest(Alpha_String_1(I),J)) then
!              Write(IOut,*) 'Alpha:',IPos+1,' I:',I,' J:',J
!              Call MQC_Print(OnePartInts%Alpha%at(IPos+1,IPos+1),IOut,'(I|h|I)')
              MatEl = MatEl + OnePartInts%Alpha%at(IPos+1,IPos+1)
            EndIf
            If(BTest(Beta_String_1(I),J)) then
!              Write(IOut,*) 'Beta:',IPos+1,' I:',I,' J:',J
              If(TwoPartInts%type().eq.'spin') then
!                Call MQC_Print(OnePartInts%Beta%at(IPos+1,IPos+1),IOut,'(I|h|I)')
                MatEl = MatEl + OnePartInts%Beta%at(IPos+1,IPos+1)
              Else
!                Call MQC_Print(OnePartInts%Alpha%at(IPos+1,IPos+1),IOut,'(I|h|I)')
                MatEl = MatEl + OnePartInts%Alpha%at(IPos+1,IPos+1)
              EndIf
            EndIf
          EndDo
        EndIf
!
        If(Present(TwoPartInts)) then
          Do IPos = 0, NBasis-1
            I = IPos/Bit_Size(0)+1
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = IPos+1, NBasis-1
              II = JPos/Bit_Size(0)+1
              JJ = Mod(JPos,Bit_Size(0)) 
              If(BTest(Alpha_String_1(I),J)) then
                If(BTest(Alpha_String_1(II),JJ)) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
!                  Call MQC_Print(TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'alpha'),IOut,'(II|JJ)')
!                  Call MQC_Print(TwoPartInts%at(IPos+1,JPos+1,JPos+1,IPos+1,'alpha'),IOut,'(IJ|JI)') 
                  MatEl = MatEl + TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'alpha') - &
                    TwoPartInts%at(IPos+1,JPos+1,JPos+1,IPos+1,'alpha')
                EndIf
              EndIf
            EndDo
          EndDo
!
          Do IPos = 0, NBasis-1
            I = IPos/Bit_Size(0)+1
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = IPos+1, NBasis-1
              II = JPos/Bit_Size(0)+1
              JJ = Mod(JPos,Bit_Size(0))
              If(BTest(Beta_String_1(I),J)) then
                If(BTest(Beta_String_1(II),JJ)) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  If(TwoPartInts%type().eq.'spin') then
!                    Call MQC_Print(TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'beta'),IOut,'(II|JJ)')
!                    Call MQC_Print(TwoPartInts%at(IPos+1,JPos+1,JPos+1,IPos+1,'beta'),IOut,'(IJ|JI)') 
                    MatEl = MatEl + TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'beta') - & 
                      TwoPartInts%at(IPos+1,JPos+1,JPos+1,IPos+1,'beta')
                  Else
!                    Call MQC_Print(TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'alpha'),IOut,'(II|JJ)')
!                    Call MQC_Print(TwoPartInts%at(IPos+1,JPos+1,JPos+1,IPos+1,'alpha'),IOut,'(IJ|JI)') 
                    MatEl = MatEl + TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'alpha') - &
                      TwoPartInts%at(IPos+1,JPos+1,JPos+1,IPos+1,'alpha')
                  EndIf
                EndIf
              EndIf
            EndDo
          EndDo
!
          Do IPos = 0, NBasis-1
            I = IPos/Bit_Size(0)+1
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = 0, NBasis-1
              II = JPos/Bit_Size(0)+1
              JJ = Mod(JPos,Bit_Size(0))
              If(BTest(Alpha_String_1(I),J)) then
                If(BTest(Beta_String_1(II),JJ)) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  If(TwoPartInts%type().eq.'spin') then
!                    Call MQC_Print(TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'alphaBeta'),IOut,'(II|JJ)')
                    MatEl = MatEl + TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'alphaBeta')
                  Else
!                    Call MQC_Print(TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'alpha'),IOut,'(II|JJ)')
                    MatEl = MatEl + TwoPartInts%at(IPos+1,IPos+1,JPos+1,JPos+1,'alpha')
                  EndIf
                EndIf
              EndIf
            EndDo
          EndDo
        EndIf

!        Call MQC_Print(MatEl,IOut,'CI Matrix Element')

        Return
!      
      Case Default
        Call MQC_Error_I('Slater_Condon is confused about number of different orbitals', 6, &
             'Det_Diff', Det_Diff )

!
      End Select
!
      End Function Slater_Condon
!
!=====================================================================
!     
!     PROCEDURE TWOERI_TRANS
! 
      Subroutine TwoERI_Trans(IOut,IPrint,MO_Coeff,ERIs,MO_ERIs,Right_MOs)
!
!     This subroutine transforms two-electron integrals to the MO basis 
!     If the UHF flag is true it computes the set of MO integrals accounting
!     for different spacial occupation of alpha and beta electrons. The order
!     in the final MO ERI array is (aa|aa), (aa|bb), (bb|aa), (bb|bb). If the
!     optional argument Right_MOs is present, then MOs in MO_Coeff are the MOs 
!     of the bra determinant, and MOs in Right_MOs are the MOs of the ket 
!     determinant, otherwise bra and ket MOs are considered equal.
!
!             |' (aa|aa) | (bb|aa) '| 
!             |----------|----------|
!             |, (aa|bb) | (bb|bb) ,|
!
!     Variable Declarations...
!
      Implicit None
      Integer(kind=int64)::IOut,IPrint,NBasis,P,Q,R,S,Mu,Nu,Lambda,Sigma,IErr
      real(kind=real64)::Zero=0.0d0
      Logical::DoON5
      Type(MQC_Matrix)::X,Y
      Type(MQC_R4Tensor)::tmpR4TensorAlpha,tmpR4TensorAlphaBeta,tmpR4TensorBetaAlpha,tmpR4TensorBeta,&
        tmpR4TensorABAB,tmpR4TensorBABA,tmpR4TensorABBA,tmpR4TensorBAAB,tmpR4TensorAAAB,tmpR4TensorAABA,&
        tmpR4TensorABAA,tmpR4TensorBAAA,tmpR4TensorBBBA,tmpR4TensorBBAB,tmpR4TensorBABB,tmpR4TensorABBB
      Type(MQC_R4Tensor)::tmpR4TensorLoc,tmpR4TensorLoc2
      Type(MQC_R4Tensor)::tmpR4Tensor
      Type(MQC_TwoERIs),Intent(InOut)::ERIs
      Type(MQC_TwoERIs),Intent(InOut)::MO_ERIs
      Type(MQC_SCF_Integral),Intent(In)::MO_Coeff
      Type(MQC_SCF_Integral),Optional,Intent(In)::Right_MOs
      Type(MQC_SCF_Integral)::MO_Coeff_2
      Character(len=15)::storage,intType
!
 1000 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',F15.8)
 2000 Format(1x,'Doing O(N**8) integral transformation algorithm')
 3000 Format(1x,'Doing O(N**5) integral transformation algorithm')
!
!     Initialize Arrays
!
      DoON5 = .True.
!
      If(ERIs%type().ne.'regular') call mqc_error_A('TwoERI_trans only implemented for &
        & regular stored 2ERIs',6,'ERIs%type()',ERIs%type())
!
      If(Present(Right_MOs))then
        MO_Coeff_2 = Right_MOs
      Else
        MO_Coeff_2 = MO_Coeff
      EndIf
!
      If(MO_Coeff%type().eq.'general'.or.MO_Coeff_2%type().eq.'general') then
        intType = 'general'
      ElseIf(MO_Coeff%type().eq.'spin'.or.MO_Coeff_2%type().eq.'spin') then
        intType = 'spin'
      Else
        intType = 'space'
      EndIf
      If(mqc_integral_norm(MO_Coeff-MO_Coeff_2).lt.1.0e-14) then
        storage = 'symm'
      Else
        storage = 'full'
      EndIf
      If(intType.eq.'space') storage = 'full'

      NBasis = MQC_Matrix_Rows(MO_Coeff%Alpha) 

      call tmpR4TensorAlpha%init(NBasis,NBasis,NBasis,NBasis)
      If(intType.eq.'spin'.or.intType.eq.'general') then
        call tmpR4TensorBeta%init(NBasis,NBasis,NBasis,NBasis)
        call tmpR4TensorAlphaBeta%init(NBasis,NBasis,NBasis,NBasis)
        If(storage.eq.'full') call tmpR4TensorBetaAlpha%init(NBasis,NBasis,NBasis,NBasis)
      EndIf
      If(intType.eq.'general') then
        call tmpR4TensorABAB%init(NBasis,NBasis,NBasis,NBasis) 
        If(storage.eq.'full') call tmpR4TensorBABA%init(NBasis,NBasis,NBasis,NBasis) 
        call tmpR4TensorABBA%init(NBasis,NBasis,NBasis,NBasis) 
        If(storage.eq.'full') call tmpR4TensorBAAB%init(NBasis,NBasis,NBasis,NBasis) 
        call tmpR4TensorAAAB%init(NBasis,NBasis,NBasis,NBasis) 
        If(storage.eq.'full') call tmpR4TensorAABA%init(NBasis,NBasis,NBasis,NBasis) 
        If(storage.eq.'full') call tmpR4TensorABAA%init(NBasis,NBasis,NBasis,NBasis) 
        If(storage.eq.'full') call tmpR4TensorBAAA%init(NBasis,NBasis,NBasis,NBasis) 
        If(storage.eq.'full') call tmpR4TensorBBBA%init(NBasis,NBasis,NBasis,NBasis) 
        If(storage.eq.'full') call tmpR4TensorBBAB%init(NBasis,NBasis,NBasis,NBasis) 
        If(storage.eq.'full') call tmpR4TensorBABB%init(NBasis,NBasis,NBasis,NBasis) 
        call tmpR4TensorABBB%init(NBasis,NBasis,NBasis,NBasis) 
      EndIf
!
      If(.not.DoON5) then

        If(IPrint.ge.4) Write(IOut,2000) 
        Do P = 1, NBasis
          Do Q = 1, NBasis
            Do R = 1, NBasis
              Do S = 1, NBasis
                Do Mu = 1, NBasis
                  Do Nu = 1, NBasis
                    Do Lambda = 1, NBasis
                      Do Sigma = 1, NBasis
                        Call tmpR4TensorAlpha%put(tmpR4TensorAlpha%at(P,Q,R,S) + &
                          conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                          ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                          MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                        If(intType.eq.'spin'.or.intType.eq.'general') then
                          Call tmpR4TensorBeta%put(tmpR4TensorBeta%at(P,Q,R,S)+&
                            conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * &
                            conjg(MO_Coeff_2%Beta%at(Lambda,R)) * MO_Coeff_2%Beta%at(Sigma,S), &
                            P,Q,R,S)
                          Call tmpR4TensorAlphaBeta%put(tmpR4TensorAlphaBeta%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * &
                            conjg(MO_Coeff_2%Beta%at(Lambda,R)) * MO_Coeff_2%Beta%at(Sigma,S), &
                            P,Q,R,S)
                          If(storage.eq.'full') &
                            Call tmpR4TensorBetaAlpha%put(tmpR4TensorBetaAlpha%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * & 
                            conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * MO_Coeff_2%Alpha%at(Sigma,S), &
                            P,Q,R,S)
                        EndIf
                        If(intType.eq.'general') then
                          !AAAA
                          Call tmpR4TensorAlpha%put(tmpR4TensorAlpha%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%alphaBeta%at(Lambda,R)) * &
                            MO_Coeff_2%alphaBeta%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorAlpha%put(tmpR4TensorAlpha%at(P,Q,R,S) + &
                            conjg(MO_Coeff%alphaBeta%at(Mu,P)) * MO_Coeff%alphaBeta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                            MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorAlpha%put(tmpR4TensorAlpha%at(P,Q,R,S) + &
                            conjg(MO_Coeff%alphaBeta%at(Mu,P)) * MO_Coeff%alphaBeta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%alphaBeta%at(Lambda,R)) * &
                            MO_Coeff_2%alphaBeta%at(Sigma,S),P,Q,R,S)
                         !BBBB
                          Call tmpR4TensorBeta%put(tmpR4TensorBeta%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorBeta%put(tmpR4TensorBeta%at(P,Q,R,S) + &
                            conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                            MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorBeta%put(tmpR4TensorBeta%at(P,Q,R,S) + &
                            conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          !AABB
                          Call tmpR4TensorAlphaBeta%put(tmpR4TensorAlphaBeta%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorAlphaBeta%put(tmpR4TensorAlphaBeta%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                            MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorAlphaBeta%put(tmpR4TensorAlphaBeta%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          !BBAA
                          If(storage.eq.'full') then
                            Call tmpR4TensorBetaAlpha%put(tmpR4TensorBetaAlpha%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBetaAlpha%put(tmpR4TensorBetaAlpha%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBetaAlpha%put(tmpR4TensorBetaAlpha%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !ABAB
                          Call tmpR4TensorABAB%put(tmpR4TensorABAB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABAB%put(tmpR4TensorABAB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                            MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABAB%put(tmpR4TensorABAB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABAB%put(tmpR4TensorABAB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                            MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          !BABA
                          If(storage.eq.'full') then
                            Call tmpR4TensorBABA%put(tmpR4TensorBABA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBABA%put(tmpR4TensorBABA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBABA%put(tmpR4TensorBABA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBABA%put(tmpR4TensorBABA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !ABBA
                          Call tmpR4TensorABBA%put(tmpR4TensorABBA%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                            MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABBA%put(tmpR4TensorABBA%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                            MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABBA%put(tmpR4TensorABBA%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                            MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABBA%put(tmpR4TensorABBA%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                            MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                          !BAAB
                          If(storage.eq.'full') then
                            Call tmpR4TensorBAAB%put(tmpR4TensorBAAB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBAAB%put(tmpR4TensorBAAB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBAAB%put(tmpR4TensorBAAB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBAAB%put(tmpR4TensorBAAB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !AAAB
                          Call tmpR4TensorAAAB%put(tmpR4TensorAAAB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorAAAB%put(tmpR4TensorAAAB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                            MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorAAAB%put(tmpR4TensorAAAB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorAAAB%put(tmpR4TensorAAAB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                            MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          !AABA
                          If(storage.eq.'full') then
                            Call tmpR4TensorAABA%put(tmpR4TensorAABA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorAABA%put(tmpR4TensorAABA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorAABA%put(tmpR4TensorAABA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorAABA%put(tmpR4TensorAABA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !ABAA
                          If(storage.eq.'full') then
                            Call tmpR4TensorABAA%put(tmpR4TensorABAA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorABAA%put(tmpR4TensorABAA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorABAA%put(tmpR4TensorABAA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorABAA%put(tmpR4TensorABAA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !BAAA
                          If(storage.eq.'full') then
                            Call tmpR4TensorBAAA%put(tmpR4TensorBAAA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBAAA%put(tmpR4TensorBAAA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBAAA%put(tmpR4TensorBAAA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBAAA%put(tmpR4TensorBAAA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !BBBA
                          If(storage.eq.'full') then
                            Call tmpR4TensorBBBA%put(tmpR4TensorBBBA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBBBA%put(tmpR4TensorBBBA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBBBA%put(tmpR4TensorBBBA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                              MO_Coeff_2%Alpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBBBA%put(tmpR4TensorBBBA%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                              MO_Coeff_2%AlphaBeta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !BBAB
                          If(storage.eq.'full') then
                            Call tmpR4TensorBBAB%put(tmpR4TensorBBAB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBBAB%put(tmpR4TensorBBAB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBBAB%put(tmpR4TensorBBAB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Alpha%at(Lambda,R)) * &
                              MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBBAB%put(tmpR4TensorBBAB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%AlphaBeta%at(Lambda,R)) * &
                              MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !BABB
                          If(storage.eq.'full') then
                            Call tmpR4TensorBABB%put(tmpR4TensorBABB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                              MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBABB%put(tmpR4TensorBABB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%BetaAlpha%at(Mu,P)) * MO_Coeff%Alpha%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                              MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBABB%put(tmpR4TensorBABB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                              MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                            Call tmpR4TensorBABB%put(tmpR4TensorBABB%at(P,Q,R,S) + &
                              conjg(MO_Coeff%Beta%at(Mu,P)) * MO_Coeff%AlphaBeta%at(Nu,Q) * &
                              ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                              MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          endIf
                          !ABBB
                          Call tmpR4TensorABBB%put(tmpR4TensorABBB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABBB%put(tmpR4TensorABBB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%Alpha%at(Mu,P)) * MO_Coeff%BetaAlpha%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                            MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABBB%put(tmpR4TensorABBB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%BetaAlpha%at(Lambda,R)) * &
                            MO_Coeff_2%BetaAlpha%at(Sigma,S),P,Q,R,S)
                          Call tmpR4TensorABBB%put(tmpR4TensorABBB%at(P,Q,R,S) + &
                            conjg(MO_Coeff%AlphaBeta%at(Mu,P)) * MO_Coeff%Beta%at(Nu,Q) * &
                            ERIs%alpha%at(Mu,Nu,Lambda,Sigma) * conjg(MO_Coeff_2%Beta%at(Lambda,R)) * &
                            MO_Coeff_2%Beta%at(Sigma,S),P,Q,R,S)
                        EndIf
                      EndDo
                    EndDo
                  EndDo
                EndDo
              EndDo
            EndDo
          EndDo
        EndDo 

      ElseIf(DoON5) then

        If(IPrint.ge.4) Write(IOut,3000) 

        Call X%init(NBasis,NBasis)
        Call Y%init(NBasis,NBasis)

        call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
        !$OMP PARALLEL DEFAULT(NONE), &
        !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
        !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
        call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
        !$OMP DO COLLAPSE(2)
        Do P = 1, NBasis
          Do Q = 1, NBasis
            X = ERIs%alpha%mat([P],[Q],[0],[0]) 
            X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
            Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
          EndDo
        EndDo
        !$OMP END DO
        !$OMP CRITICAL
        tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
        !$OMP END CRITICAL
        !$OMP END PARALLEL

        !$OMP PARALLEL DEFAULT(NONE), &
        !$OMP SHARED(tmpR4Tensor,tmpR4TensorAlpha,MO_Coeff,MO_Coeff_2,nBasis), &
        !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
        call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
        !$OMP DO COLLAPSE(2)
        Do R = 1, NBasis
          Do S = 1, NBasis
            X = tmpR4Tensor%mat([0],[0],[R],[S])
            X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
            Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
          EndDo
        EndDo
        !$OMP END DO
        !$OMP CRITICAL
        tmpR4TensorAlpha = tmpR4TensorAlpha + tmpR4TensorLoc
        !$OMP END CRITICAL
        !$OMP END PARALLEL

        If(intType.eq.'spin'.or.intType.eq.'general') then

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0]) 
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorBeta,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorBeta = tmpR4TensorBeta + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL
         
          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAlphaBeta,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAlphaBeta = tmpR4TensorAlphaBeta + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          If(storage.eq.'full') then

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBetaAlpha,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBetaAlpha = tmpR4TensorBetaAlpha + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

          EndIf

        EndIf

        If(intType.eq.'general') then

!         AAAA
          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0]) 
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAlpha,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alpha).dot.X.dot.MO_Coeff%alpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAlpha = tmpR4TensorAlpha + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%alpha).dot.X.dot.MO_Coeff%alpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAlpha,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAlpha = tmpR4TensorAlpha + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAlpha,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAlpha = tmpR4TensorAlpha + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL
!
!         BBBB
          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorBeta,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorBeta = tmpR4TensorBeta + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorBeta,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%beta).dot.X.dot.MO_Coeff%beta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorBeta = tmpR4TensorBeta + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%betaAlpha).dot.X.dot.MO_Coeff%betaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorBeta,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%betaAlpha).dot.X.dot.MO_Coeff%betaAlpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorBeta = tmpR4TensorBeta + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

!         AABB
          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%betaAlpha).dot.X.dot.MO_Coeff%betaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAlphaBeta,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alpha).dot.X.dot.MO_Coeff%alpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAlphaBeta = tmpR4TensorAlphaBeta + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta 
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAlphaBeta,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAlphaBeta = tmpR4TensorAlphaBeta + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%betaAlpha).dot.X.dot.MO_Coeff%betaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAlphaBeta,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAlphaBeta = tmpR4TensorAlphaBeta + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

!         BBAA
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBetaAlpha,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%beta).dot.X.dot.MO_Coeff%beta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBetaAlpha = tmpR4TensorBetaAlpha + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%alpha).dot.X.dot.MO_Coeff%alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBetaAlpha,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%betaAlpha).dot.X.dot.MO_Coeff%betaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBetaAlpha = tmpR4TensorBetaAlpha + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBetaAlpha,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%betaAlpha).dot.X.dot.MO_Coeff%betaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBetaAlpha = tmpR4TensorBetaAlpha + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         ABAB
          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABAB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABAB = tmpR4TensorABAB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABAB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABAB = tmpR4TensorABAB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABAB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABAB = tmpR4TensorABAB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABAB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABAB = tmpR4TensorABAB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

!         BABA
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBABA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBABA = tmpR4TensorBABA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBABA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBABA = tmpR4TensorBABA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBABA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBABA = tmpR4TensorBABA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBABA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBABA = tmpR4TensorBABA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         ABBA
          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABBA,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABBA = tmpR4TensorABBA + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABBA,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABBA = tmpR4TensorABBA + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABBA,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABBA = tmpR4TensorABBA + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABBA,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABBA = tmpR4TensorABBA + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

!         BAAB
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBAAB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBAAB = tmpR4TensorBAAB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBAAB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBAAB = tmpR4TensorBAAB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBAAB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBAAB = tmpR4TensorBAAB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBAAB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBAAB = tmpR4TensorBAAB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         AAAB
          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAAAB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAAAB = tmpR4TensorAAAB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAAAB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAAAB = tmpR4TensorAAAB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAAAB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%AlphaBeta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAAAB = tmpR4TensorAAAB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorAAAB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorAAAB = tmpR4TensorAAAB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

!         AABA
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorAABA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorAABA = tmpR4TensorAABA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorAABA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorAABA = tmpR4TensorAABA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorAABA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorAABA = tmpR4TensorAABA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorAABA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%alphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorAABA = tmpR4TensorAABA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         ABAA
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorABAA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorABAA = tmpR4TensorABAA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorABAA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorABAA = tmpR4TensorABAA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorABAA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorABAA = tmpR4TensorABAA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorABAA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorABAA = tmpR4TensorABAA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         BAAA
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBAAA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBAAA = tmpR4TensorBAAA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBAAA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBAAA = tmpR4TensorBAAA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBAAA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBAAA = tmpR4TensorBAAA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBAAA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBAAA = tmpR4TensorBAAA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         BBBA
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBBBA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBBBA = tmpR4TensorBBBA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBBBA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBBBA = tmpR4TensorBBBA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBBBA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBBBA = tmpR4TensorBBBA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBBBA,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBBBA = tmpR4TensorBBBA + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         BBAB
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBBAB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBBAB = tmpR4TensorBBAB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBBAB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBBAB = tmpR4TensorBBAB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBBAB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBBAB = tmpR4TensorBBAB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBBAB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBBAB = tmpR4TensorBBAB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         BABB
          If(storage.eq.'full') then
            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBABB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBABB = tmpR4TensorBABB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBABB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%Alpha
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBABB = tmpR4TensorBABB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBABB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBABB = tmpR4TensorBABB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                X = ERIs%alpha%mat([P],[Q],[0],[0])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
                Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL

            !$OMP PARALLEL DEFAULT(NONE), &
            !$OMP SHARED(tmpR4Tensor,tmpR4TensorBABB,MO_Coeff,MO_Coeff_2,nBasis), &
            !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
            call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
            !$OMP DO COLLAPSE(2)
            Do R = 1, NBasis
              Do S = 1, NBasis
                X = tmpR4Tensor%mat([0],[0],[R],[S])
                X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%AlphaBeta
                Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
              EndDo
            EndDo
            !$OMP END DO
            !$OMP CRITICAL
            tmpR4TensorBABB = tmpR4TensorBABB + tmpR4TensorLoc
            !$OMP END CRITICAL
            !$OMP END PARALLEL
          EndIf

!         ABBB
          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABBB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABBB = tmpR4TensorABBB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABBB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%Alpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABBB = tmpR4TensorABBB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%BetaAlpha).dot.X.dot.MO_Coeff%BetaAlpha
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABBB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%AlphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABBB = tmpR4TensorABBB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          call tmpR4Tensor%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(ERIs,tmpR4Tensor,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,P,Q,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do P = 1, NBasis
            Do Q = 1, NBasis
              X = ERIs%alpha%mat([P],[Q],[0],[0])
              X = Dagger(MO_Coeff_2%Beta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[P],[Q],[0],[0])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4Tensor = tmpR4Tensor + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

          !$OMP PARALLEL DEFAULT(NONE), &
          !$OMP SHARED(tmpR4Tensor,tmpR4TensorABBB,MO_Coeff,MO_Coeff_2,nBasis), &
          !$OMP PRIVATE(X,R,S,tmpR4TensorLoc)
          call tmpR4TensorLoc%init(nBasis,nBasis,nBasis,nBasis)
          !$OMP DO COLLAPSE(2)
          Do R = 1, NBasis
            Do S = 1, NBasis
              X = tmpR4Tensor%mat([0],[0],[R],[S])
              X = Dagger(MO_Coeff_2%alphaBeta).dot.X.dot.MO_Coeff%Beta
              Call tmpR4TensorLoc%mput(X,[0],[0],[R],[S])
            EndDo
          EndDo
          !$OMP END DO
          !$OMP CRITICAL
          tmpR4TensorABBB = tmpR4TensorABBB + tmpR4TensorLoc
          !$OMP END CRITICAL
          !$OMP END PARALLEL

        EndIf
    
      EndIf
      
      If(intType.eq.'space') then
        call mqc_twoeris_allocate(mo_eris,'full','space',tmpR4TensorAlpha)
      ElseIf(intType.eq.'spin'.and.storage.eq.'full') then
        call mqc_twoeris_allocate(mo_eris,'full','spin',tmpR4TensorAlpha,tmpR4TensorBeta, &
          tmpR4TensorAlphaBeta,tmpR4TensorBetaAlpha)
      ElseIf(intType.eq.'spin'.and.storage.eq.'symm') then
        call mqc_twoeris_allocate(mo_eris,'symm','spin',tmpR4TensorAlpha,tmpR4TensorBeta, &
          tmpR4TensorAlphaBeta)
      ElseIf(intType.eq.'general'.and.storage.eq.'full') then
        call mqc_twoeris_allocate(mo_eris,'full','general',tmpR4TensorAlpha,tmpR4TensorBeta, &
          tmpR4TensorAlphaBeta,tmpR4TensorBetaAlpha,tmpR4TensorABAB,tmpR4TensorABBA,tmpR4TensorBAAB, &
          tmpR4TensorBABA,tmpR4TensorAAAB,tmpR4TensorAABA,tmpR4TensorABAA,tmpR4TensorBAAA,tmpR4TensorBBBA,&
          tmpR4TensorBBAB,tmpR4TensorBABB,tmpR4TensorABBB)
      ElseIf(intType.eq.'general'.and.storage.eq.'symm') then
        call mqc_twoeris_allocate(mo_eris,'symm','general',tmpR4TensorAlpha,tmpR4TensorBeta, &
          tmpR4TensorAlphaBeta,abab=tmpR4TensorABAB,abba=tmpR4TensorABBA,aaab=tmpR4TensorAAAB,&
          abbb=tmpR4TensorABBB)
      Else
        call mqc_error_a('TwoERI_Trans is confused about ERI storage or spin symmetry',6,'storage',&
          storage,'intType',intType)
      EndIf

      If(IPrint.ge.4) then
        Call MQC_Print(MO_ERIs%alpha,IOut,'Transformed MO 2ERIs (aa|aa)')
        If(intType.eq.'spin'.or.intType.eq.'general') then
          Call MQC_Print(MO_ERIs%beta,IOut,'Transformed MO 2ERIs (bb|bb)')
          Call MQC_Print(MO_ERIs%alphaBeta,IOut,'Transformed MO 2ERIs (aa|bb)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%BetaAlpha,IOut,'Transformed MO 2ERIs (bb|aa)') 
        EndIf
        If(intType.eq.'general') then
          Call MQC_Print(MO_ERIs%abab,IOut,'Transformed MO 2ERIs (ab|ab)') 
          Call MQC_Print(MO_ERIs%abba,IOut,'Transformed MO 2ERIs (ab|ba)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%baab,IOut,'Transformed MO 2ERIs (ba|ab)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%baba,IOut,'Transformed MO 2ERIs (ba|ba)')
          Call MQC_Print(MO_ERIs%aaab,IOut,'Transformed MO 2ERIs (aa|ab)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%aaba,IOut,'Transformed MO 2ERIs (aa|ba)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%abaa,IOut,'Transformed MO 2ERIs (ab|aa)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%baaa,IOut,'Transformed MO 2ERIs (ba|aa)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%bbba,IOut,'Transformed MO 2ERIs (bb|ba)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%bbab,IOut,'Transformed MO 2ERIs (bb|ab)')
          If(storage.eq.'full') Call MQC_Print(MO_ERIs%babb,IOut,'Transformed MO 2ERIs (ba|bb)')
          Call MQC_Print(MO_ERIs%abbb,IOut,'Transformed MO 2ERIs (ab|bb)')
        EndIf
      EndIf
!
      End Subroutine TwoERI_Trans


!
!=====================================================================
!
!     PROCEDURE S2_MAT_ELEM    
! 
!>    \brief <b> S2_Mat_Elem is a function returns the CI S**2 matrix elements</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    S2_Mat_Elem is a function returns the CI S**2 matrix elements used for computing 
!>    S**2 values of CI vectors for a given alpha and beta string combination. The MO 
!>    overlap matrix is required
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64)
!>        The FORTRAN file to print to. 
!>    \endverbatim
!>    
!>    \param[in] IPrint
!>    \verbatim
!>        IPrint is Integer(kind=int64)
!>        The print level for the subroutine. 
!>    \endverbatim
!>    
!>    \param[in] NBasisIn
!>    \verbatim
!>        NBasisIn is Type(MQC_Scalar)
!>        The number of basis functions. 
!>    \endverbatim
!>    
!>    \param[in] Alpha_String_1
!>    \verbatim
!>        Alpha_String_1 is Integer(kind=int64),Dimension(:),Allocatable
!>        The left alpha string.
!>    \endverbatim
!>    
!>    \param[in] Beta_String_1
!>    \verbatim
!>        Beta_String_1 is Integer(kind=int64),Dimension(:),Allocatable
!>        The left beta string.
!>    \endverbatim
!>    
!>    \param[in] Alpha_String_2
!>    \verbatim
!>        Alpha_String_2 is Integer(kind=int64),Dimension(:),Allocatable
!>        The right alpha string.
!>    \endverbatim
!>    
!>    \param[in] Beta_String_2
!>    \verbatim
!>        Beta_String_2 is Integer(kind=int64),Dimension(:),Allocatable
!>        The right beta string.
!>    \endverbatim
!>    
!>    \param[in] MO_Overlap
!>    \verbatim
!>        MO_Overlap is Type(MQC_SCF_Integral)
!>        The MO basis core hamiltonian matrix.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2016
!
      Function S2_Mat_Elem(IOut,IPrint,NBasisIn,Alpha_String_1,Beta_String_1, &
        Alpha_String_2,Beta_String_2,MO_Overlap)
!
!     Variable Declarations...
!
      Implicit None
      Integer(kind=int64),Intent(In)::IOut,IPrint
      Type(MQC_Scalar),Intent(In)::NBasisIn
!      Type(MQC_Determinant),Intent(In)::Determinants
      Integer(kind=int64),Dimension(4)::Orbs,Spin,Det
      Integer(kind=int64),Dimension(:),Allocatable,Intent(In)::Alpha_String_1,Alpha_String_2, &
        Beta_String_1,Beta_String_2
      Integer(kind=int64)::NBasis,IPos,JPos,IDiff,Det_Diff,NAlpha,NBeta, &
        IOcc,JOcc,KOcc,LOcc,Mat_Sign,Alpha_Diff_Cnt,Beta_Diff_Cnt,NBit_Ints, &
        I,J,II,JJ
      real(kind=real64)::Zero=0.0d0,Quarter=0.25d0,ABTerm,One=1.0d0
      type(mqc_scalar)::S2_Mat_Elem
      Type(MQC_SCF_Integral),Intent(In)::MO_Overlap
      Integer(kind=int64),Dimension(:),Allocatable::Alpha_Diff,Beta_Diff

      NBasis = NBasisIn
      NBit_Ints = (NBasis/Bit_Size(0))+1 
!      Allocate(Alpha_String_1(NBit_Ints),Alpha_String_2(NBit_Ints),Beta_String_1(NBit_Ints),Beta_String_2(NBit_Ints))
!      Alpha_String_1 = Determinants%Strings%Alpha%vat([L_A_String],[1,NBit_Ints])
!      Alpha_String_2 = Determinants%Strings%Alpha%vat([R_A_String],[1,NBit_Ints])
!      Beta_String_1 = Determinants%Strings%Beta%vat([L_B_String],[1,NBit_Ints])
!      Beta_String_2 = Determinants%Strings%Beta%vat([R_B_String],[1,NBit_Ints])
!
!      Write(IOut,*) 'alpha 1:'
!      Write(IOut,'(B64)') Alpha_String_1
!      Write(IOut,*) 'alpha 2:'
!      Write(IOut,'(B64)') Alpha_String_2
!      Write(IOut,*) 'alpha XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Alpha_String_1,Alpha_String_2)   
!      Write(IOut,*) 'NDifa=', PopCnt(IEOR(Alpha_String_1,Alpha_String_2))
!      Write(IOut,*)
!      Write(IOut,*) 'beta 1:'
!      Write(IOut,'(B64)') Beta_String_1
!      Write(IOut,*) 'beta 2:'
!      Write(IOut,'(B64)') Beta_String_2
!      Write(IOut,*) 'beta XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Beta_String_1,Beta_String_2)   
!      Write(IOut,*) 'NDifb=', PopCnt(IEOR(Beta_String_1,Beta_String_2))
!      Write(IOut,*)
!
      Allocate(Alpha_Diff(NBit_Ints),Beta_Diff(NBit_Ints))
      Det_Diff = 0
      Alpha_Diff_Cnt = 0
      Beta_Diff_Cnt = 0
      NAlpha = 0
      NBeta = 0
      Do I = 1,NBit_Ints
        Alpha_Diff(I) = IEOR(Alpha_String_1(I),Alpha_String_2(I))
!        Write(IOut,*) 'Alpha Diff',I,':'
!        Write(IOut,'(B64)') Alpha_Diff(I)
!        Write(IOut,*) '-------------'
        Alpha_Diff_Cnt = Alpha_Diff_Cnt + PopCnt(Alpha_Diff(I)) 
        NAlpha = NAlpha + PopCnt(Alpha_String_1(I))
        Beta_Diff(I) = IEOR(Beta_String_1(I),Beta_String_2(I))
!        Write(IOut,*) 'Beta Diff',I,':'
!        Write(IOut,'(B64)') Beta_Diff(I)
!        Write(IOut,*) '-------------'
        Beta_Diff_Cnt = Beta_Diff_Cnt + PopCnt(Beta_Diff(I))
        NBeta = NBeta + PopCnt(Beta_String_1(I))
      EndDo
!      Write(IOut,*)'Alpha_Diff_Cnt:',Alpha_Diff_Cnt,'Beta_Diff_Cnt:',Beta_Diff_Cnt
      Det_Diff = Alpha_Diff_Cnt/2 + Beta_Diff_Cnt/2

      If(Mod(Alpha_Diff_Cnt,2).ne.0.or.Mod(Beta_Diff_Cnt,2).ne.0) &
        Call MQC_Error_I('S2_Mat_Elem has been handed spin non-conserving &
        &determinants',6,'Mod(Alpha_Diff_Cnt,2)', Mod(Alpha_Diff_Cnt,2), &
        'Mod(Beta_Diff_Cnt,2)', Mod(Beta_Diff_Cnt,2) )
      
      Select Case (Det_Diff)
!
        Case(3:)
          S2_Mat_Elem = Zero 
          Return
!
        Case(2)
          IDiff = 1
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J)) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J)) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J)) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J)) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Call Print_Vector(IOut,Orbs,'Orbs')
!          Call Print_Vector(IOut,Spin,'Spin')
!          Call Print_Vector(IOut,Det,'Det')
!
          IOcc = 0
          Do IPos = Orbs(1)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(1).eq.0) then
              If(Det(1).eq.1) then
                If(BTest(Alpha_String_1(I),J)) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Alpha_String_2(I),J)) IOcc = IOcc + 1
              EndIf
            ElseIf(Spin(1).eq.1) then
              If(Det(1).eq.1) then
                If(BTest(Beta_String_1(I),J)) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Beta_String_2(I),J)) IOcc = IOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'IOcc:',IOcc
          JOcc = 0
          Do IPos = Orbs(2)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(2).eq.0) then
              If(Det(2).eq.1) then
                If(BTest(Alpha_String_1(I),J)) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Alpha_String_2(I),J)) JOcc = JOcc + 1
              EndIf
            ElseIf(Spin(2).eq.1) then
              If(Det(2).eq.1) then
                If(BTest(Beta_String_1(I),J)) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Beta_String_2(I),J)) JOcc = JOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'JOcc:',JOcc
          KOcc = 0
          Do IPos = Orbs(3)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(3).eq.0) then
              If(Det(3).eq.1) then
                If(BTest(Alpha_String_1(I),J)) KOcc = KOcc + 1
              ElseIf(Det(3).eq.2) then
                If(BTest(Alpha_String_2(I),J)) KOcc = KOcc + 1
              EndIf
            ElseIf(Spin(3).eq.1) then
              If(Det(3).eq.1) then
                If(BTest(Beta_String_1(I),J)) KOcc = KOcc + 1
              ElseIf(Det(3).eq.2) then
                If(BTest(Beta_String_2(I),J)) KOcc = KOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'KOcc:',KOcc
          LOcc = 0
          Do IPos = Orbs(4)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(4).eq.0) then
              If(Det(4).eq.1) then
                If(BTest(Alpha_String_1(I),J)) LOcc = LOcc + 1
              ElseIf(Det(4).eq.2) then
                If(BTest(Alpha_String_2(I),J)) LOcc = LOcc + 1
              EndIf
            ElseIf(Spin(4).eq.1) then
              If(Det(4).eq.1) then
                If(BTest(Beta_String_1(I),J)) LOcc = LOcc + 1
              ElseIf(Det(4).eq.2) then
                If(BTest(Beta_String_2(I),J)) LOcc = LOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'LOcc:',LOcc
!          Mat_Sign = -1
!          Mat_Sign = (-1)**(IOcc+JOcc+KOcc+LOcc-3)
!          Write(IOut,*) 'Permutations:',(IOcc+JOcc+KOcc+LOcc-3)
          Mat_Sign = (-1)**(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!          Write(IOut,*) 'Permutations:',(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!          Write(IOut,*) 'Mat_Sign:',Mat_Sign
!
          If(Det(1).eq.Det(2).and.Det(3).eq.Det(4)) then
            If(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(2).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(1),Orbs(4))*&
                  MO_Overlap%BetaAlpha%at(Orbs(2),Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(2),Orbs(3))*&
                  MO_Overlap%BetaAlpha%at(Orbs(1),Orbs(4))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              If(Spin(1).eq.0.and.Spin(2).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(1),Orbs(3))*&
                  MO_Overlap%BetaAlpha%at(Orbs(2),Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(2),Orbs(4))*&
                  MO_Overlap%BetaAlpha%at(Orbs(1),Orbs(3))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          ElseIf(Det(1).eq.Det(3).and.Det(2).eq.Det(4)) then
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(3).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(1),Orbs(4))*&
                  MO_Overlap%BetaAlpha%at(Orbs(3),Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(3),Orbs(2))*&
                  MO_Overlap%BetaAlpha%at(Orbs(1),Orbs(4))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              If(Spin(1).eq.0.and.Spin(3).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(1),Orbs(2))*&
                  MO_Overlap%BetaAlpha%at(Orbs(3),Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(3),Orbs(4))*&
                MO_Overlap%BetaAlpha%at(Orbs(1),Orbs(2))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          ElseIf(Det(1).eq.Det(4).and.Det(2).eq.Det(3)) then
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(4).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(1),Orbs(3))*&
                  MO_Overlap%BetaAlpha%at(Orbs(4),Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(4),Orbs(2))*&
                  MO_Overlap%BetaAlpha%at(Orbs(1),Orbs(3))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(4).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(1),Orbs(2))*&
                  MO_Overlap%BetaAlpha%at(Orbs(4),Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap%AlphaBeta%at(Orbs(4),Orbs(3))*&
                  MO_Overlap%BetaAlpha%at(Orbs(1),Orbs(2))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          EndIf
!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem

          Return

        Case(1)
          IDiff = 1
!          Allocate(Orbs(2),Spin(2))
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J)) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J)) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J)) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J)) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Write(IOut,*)'Orb 1:',Orbs(1),' Orb 2:',Orbs(2)
!          Write(IOut,*)'Spin 1:',Spin(1),' Spin 2:',Spin(2)
!          Write(IOut,*)'Det 1:',Det(1),' Det 2:',Det(2)
!
          S2_Mat_Elem = Zero 
          If(Spin(1).ne.Spin(2)) then
            S2_Mat_Elem = Zero
!
          ElseIf(Spin(1).eq.0) then
!
            IOcc = 0
            Do IPos = Orbs(1)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(1).eq.1) then
                If(BTest(Alpha_String_1(I),J)) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Alpha_String_2(I),J)) IOcc = IOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'IOcc:',IOcc
            JOcc = 0
            Do IPos = Orbs(2)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(2).eq.1) then
                If(BTest(Alpha_String_1(I),J)) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Alpha_String_2(I),J)) JOcc = JOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'JOcc:',JOcc
!
            Do IPos = 0, NBasis-1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(BTest(Beta_String_1(I),J)) then
                S2_Mat_Elem = S2_Mat_Elem + MO_Overlap%BetaAlpha%at(IPos+1,Orbs(1))*& 
                  MO_Overlap%AlphaBeta%at(Orbs(2),IPos+1)
              EndIf
            EndDo
!            S2_Mat_Elem = - S2_Mat_Elem
!            Write(IOut,*) 'Permutations:',(2*NAlpha+1-IOcc-JOcc)
!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NAlpha+1-IOcc-JOcc)
            S2_Mat_Elem = (-1)**(2*NAlpha+1-IOcc-JOcc) * S2_Mat_Elem
!            S2_Mat_Elem = (-1)**(IOcc+JOcc-1) * S2_Mat_Elem
!
          ElseIf(Spin(1).eq.1) then

            IOcc = 0
            Do IPos = Orbs(1)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(1).eq.1) then
                If(BTest(Beta_String_1(I),J)) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Beta_String_2(I),J)) IOcc = IOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'IOcc:',IOcc
            JOcc = 0
            Do IPos = Orbs(2)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(2).eq.1) then
                If(BTest(Beta_String_1(I),J)) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Beta_String_2(I),J)) JOcc = JOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'JOcc:',JOcc
!
            Do IPos = 0, NBasis-1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(BTest(Alpha_String_1(I),J)) then
                S2_Mat_Elem = S2_Mat_Elem + MO_Overlap%AlphaBeta%at(IPos+1,Orbs(1))*&
                  MO_Overlap%BetaAlpha%at(Orbs(2),IPos+1)
              EndIf
            EndDo
!            S2_Mat_Elem = - S2_Mat_Elem
!            Write(IOut,*) 'Permutations:',(2*NBeta+1-IOcc-JOcc)
!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NBeta+1-IOcc-JOcc)
            S2_Mat_Elem = (-1)**(2*NBeta+1-IOcc-JOcc) * S2_Mat_Elem
!            S2_Mat_Elem = (-1)**(IOcc+JOcc-1) * S2_Mat_Elem

          EndIf

!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem

          Return
!
        Case(0)
          ABTerm = Zero
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = 0, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0)) 
              If((BTest(Alpha_String_2(I),J)).and.(BTest(Beta_String_2(II),JJ))) then
                ABTerm = ABTerm + MO_Overlap%AlphaBeta%at(IPos+1,JPos+1)*&
                  MO_Overlap%BetaAlpha%at(JPos+1,IPos+1) 
              EndIf
            EndDo
          EndDo
          S2_Mat_Elem = Quarter*((NAlpha-NBeta)**2+2*(NAlpha+NBeta)) - ABTerm
!          Call S2_Mat_Elem%print(6,'S2_Mat_Elem')
!
          Return
!
      End Select
!
      End Function S2_Mat_Elem  
!
!
!=====================================================================
!     
!     PROCEDURE MQC_BUILD_CI_HAMILTONIAN
!
!>    \brief <b> MQC_BUILD_CI_HAMILTONIAN is a subroutine that builds the CI 
!>    operator matrix given MO integrals and determinant strings</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_BUILD_CI_HAMILTONIAN is a subroutine that builds the CI operator matrix 
!>    given MO integrals and determinant strings. The name of the routine reflects
!>    its original purpose of returning the CI Hamiltonian matrix, but in fact it
!>    can be used for any one, two, or one and two particle operator. The routine
!>    can also be used to build the S^2 operator matrix in the CI basis.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64)
!>        The FORTRAN file to print to. 
!>    \endverbatim
!>    
!>    \param[in] IPrint
!>    \verbatim
!>        IPrint is Integer(kind=int64)
!>        The print level for the subroutine. 
!>    \endverbatim
!>    
!>    \param[in] NBasis
!>    \verbatim
!>        NBasis is Type(MQC_Scalar)
!>        The number of basis functions. 
!>    \endverbatim
!>    
!>    \param[in] Determinants
!>    \verbatim
!>        Determinants is Type(MQC_Determinant)
!>        The binary string occupation number vectors. If 
!>        opional argument Dets2 is provided Determinants
!>        variable will contain the row basis strings.
!>    \endverbatim
!>    
!>    \param[in] MO_Core_Ham
!>    \verbatim
!>        MO_Core_Ham is Type(MQC_SCF_Integral),Optional
!>        The MO basis one-particle integrals.
!>    \endverbatim
!>    
!>    \param[in] MO_ERIs
!>    \verbatim
!>        MO_ERIs is Type(MQC_TwoERIs),Optional 
!>        The MO basis two-particle integrals.
!>    \endverbatim
!>    
!>    \param[out] CI_Hamiltonian
!>    \verbatim
!>        CI_Hamiltonian is Type(MQC_Matrix)
!>        The CI hamiltonian returned.
!>    \endverbatim
!>
!>    \param[in] Subs
!>    \verbatim
!>        Subs is Integer(kind=int64),Dimension(:),Optional
!>        Permitted substitutions to include in the Hamiltonian.
!>        If not present, all substitutions are included.
!>    \endverbatim
!>    
!>    \param[in] Dets2
!>    \verbatim
!>        Dets2 is Type(MQC_Determinant),Optional
!>        If provided, gives the binary string occupation number 
!>        vectors for the column basis, where Determinants variable
!>        will continue to provide the binary string occupation 
!>        number vectors for the row basis.
!>    \endverbatim
!>    
!>    \param[in] Subs2
!>    \verbatim
!>        Subs2 is Integer(kind=int64),Dimension(:),Optional
!>        If provided, gives the permitted substitutions to include 
!>        in the Hamiltonian for the column basis. If not present, 
!>        column basis substitutions are equal to row basis 
!>        substitutions provided by variable Subs.
!>    \endverbatim
!>    
!>    \param[in] doS2
!>    \verbatim
!>        DoS2 is Logical,Optional
!>        Logical flag to build the S^2 matrix instead of using 
!>        Slater Condon rules. The spatial molecular orbital basis
!>        overlap (i.e. do not multiply by spin overlap which 
!>        always leads to zero alpha-beta blocks) must be provided
!>        in MO_Core_Ham argument if this option is used) must be provided
!>        in MO_Core_Ham argument if this option is used..
!>    \endverbatim
!>    
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2017,2021,2024
!
      Subroutine MQC_Build_CI_Hamiltonian(IOut,IPrint,NBasis,Determinants, &
          MO_Core_Ham,MO_ERIs,CI_Hamiltonian,Subs,Dets2,Subs2,doS2)
!
      Implicit None
      Integer(kind=int64),Intent(In)::IOut,IPrint
      Type(MQC_Scalar),Intent(In)::NBasis
      Type(MQC_TwoERIs),Optional,Intent(In)::MO_ERIs
      Type(MQC_SCF_Integral),Optional,Intent(In)::MO_Core_Ham
      Type(MQC_Determinant),Intent(In)::Determinants
      Type(MQC_Determinant),Optional,Intent(In)::Dets2
      Logical,Optional,Intent(In)::doS2
      Type(MQC_Matrix),Intent(Out)::CI_Hamiltonian
      Integer(kind=int64)::I,J,K,NAlpha_Str1,NBeta_Str1,NDets1,L_A_String,L_B_String, &
        R_A_String,R_B_String,L_Index,R_Index,L_A_Start,L_B_Start,L_A_End,L_B_End, &
        R_A_Start,R_B_Start,R_A_End,R_B_End,L_A_Sub,L_B_Sub,R_A_Sub,R_B_Sub,Counter, &
        NBit_Ints,NAlpha_Str2,NBeta_Str2,NDets2,Counter2
      Integer(kind=int64),Dimension(:),Optional::Subs,Subs2
      Integer(kind=int64),Dimension(:),Allocatable::Alpha_String_1,Alpha_String_2,Beta_String_1, &
        Beta_String_2,rightSubs
      Type(MQC_Determinant)::Determinants2
      Logical::SymmFlag,S2Flag
      Character(len=12)::SymmStr
      Character(len=256)::detIndex1,detIndex2
      Character(len=44)::bar=" Build Progress: ???% |                    |"
!
      If(present(doS2)) then
        S2Flag = doS2
      else
        S2Flag = .false.
      EndIf
!
      If(Present(MO_Core_Ham)) then
        If(.not.S2Flag.and.MO_Core_Ham%type().ne.'space'.and.MO_Core_Ham%type().ne.'spin') &
          call mqc_error_A('Slater_Condon only implemented for &
          & spin or space one-particle integrals',6,'MO_Core_Ham%type()',MO_Core_Ham%type())
      EndIf
      If(Present(MO_ERIs)) then
        If(MO_ERIs%type().ne.'space'.and.MO_ERIs%type().ne.'spin') call mqc_error_A('Slater_Condon only implemented for &
          & spin or space two-particle integrals',6,'MO_ERIs%type()',MO_ERIs%type())
      endIf
      If(.not.present(MO_Core_Ham).and..not.present(MO_ERIs)) call mqc_error('No integrals given to Slater_Condon')
!
      If(present(Dets2)) then
        Determinants2 = Dets2
        SymmFlag = .False.
      Else
        Determinants2 = Determinants
        SymmFlag = .True.
      EndIf
!
!     Need to figure out dimensions of CI Hamiltonian, alpha strings and beta stings
      NAlpha_Str1 = MQC_Matrix_Rows(Determinants%Strings%Alpha)
      NBeta_Str1 = MQC_Matrix_Rows(Determinants%Strings%Beta)
      If(present(subs)) then
        NDets1 = 0
        Do I = 1, Size(Subs)
          Do J = 1, Size(Determinants%NSubsAlpha)
            Do K = 1, Size(Determinants%NSubsBeta)
              If((J-1)+(K-1).eq.Subs(i)) then
                NDets1 = NDets1 + Determinants%NSubsAlpha%at(J)*&
                  Determinants%NSubsBeta%at(K)
              EndIf
            EndDo
          EndDo
        EndDo
      Else
        NDets1 = NAlpha_Str1 * NBeta_Str1
      EndIf

      NAlpha_Str2 = MQC_Matrix_Rows(Determinants2%Strings%Alpha)
      NBeta_Str2 = MQC_Matrix_Rows(Determinants2%Strings%Beta)
      If(present(subs2)) then
        rightSubs = subs2
        SymmFlag = .False.
      elseIf(present(subs)) then
        rightSubs = subs
      endIf
      If(allocated(rightSubs)) then
        NDets2 = 0
        Do I = 1, Size(rightSubs)
          Do J = 1, Size(Determinants2%NSubsAlpha)
            Do K = 1, Size(Determinants2%NSubsBeta)
              If((J-1)+(K-1).eq.rightSubs(i)) then
                NDets2 = NDets2 + Determinants2%NSubsAlpha%at(J)*&
                  Determinants2%NSubsBeta%at(K)
              EndIf
            EndDo
          EndDo
        EndDo
      Else
        NDets2 = NAlpha_Str2 * NBeta_Str2 
      EndIf
!
      If(SymmFlag) then
        Call CI_Hamiltonian%init(NDets1,NDets2,Storage='StorHerm')
        SymmStr = 'hermitian'
      Else
        Call CI_Hamiltonian%init(NDets1,NDets2)
        SymmStr = 'element'
      EndIf
      NBit_Ints = (NBasis/(Bit_Size(0)-1))+1 
      Allocate(Alpha_String_1(NBit_Ints),Alpha_String_2(NBit_Ints), &
        Beta_String_1(NBit_Ints),Beta_String_2(NBit_Ints))
!
      If(Determinants%order.ne.Determinants2%order) &
        Call MQC_Error_A('Left and Right determinant sets have different order in &
        &MQC_Build_CI_Hamiltonian',6,'Determinants%order',Determinants%order,&
        'Determinants2%order',Determinants2%order)
      select case (Determinants%order)
      case ('lexical')
        Counter = 0
        Counter2 = 0
        Do L_A_String = 1, NAlpha_Str1  
          Do L_B_String = 1, NBeta_Str1  
            Do R_A_String = 1, NAlpha_Str2  
              Do R_B_String = 1, NBeta_Str2  
                L_Index = 1+(L_B_String-1)*NAlpha_Str1+(L_A_String-1) 
                R_Index = 1+(R_B_String-1)*NAlpha_Str2+(R_A_String-1) 
                if(SymmFlag.and.L_Index.lt.R_Index) cycle
                if(iPrint.ge.2.and.R_Index.eq.1) then
                  Counter = Counter + 1
                  if(L_Index.eq.1) then
                    if(SymmFlag) then
                      DetIndex1 = NEW_LINE('A')//' CI determinant index'//NEW_LINE('A')
                    else
                      DetIndex1 = NEW_LINE('A')//' Left CI determinant index'//NEW_LINE('A')
                    endIf
                  endIf
                  DetIndex1 = trim(DetIndex1)//NEW_LINE('A')//' '//trim(num2char(Counter,'I6'))//& 
                    ' = '//trim(mqc_detstring_print(Determinants%Strings%Alpha%vat([L_A_String],[0]),&
                    Determinants%Strings%Beta%vat([L_B_String],[0]),Int(NBasis)))
                endIf
                if(.not.SymmFlag.and.iPrint.ge.2.and.L_Index.eq.1) then
                  Counter2 = Counter2 + 1
                  if(R_Index.eq.1) DetIndex2 = NEW_LINE('A')//' Right CI determinant index'//NEW_LINE('A')
                  DetIndex2 = trim(DetIndex2)//NEW_LINE('A')//' '//trim(num2char(Counter2,'I6'))// & 
                    ' = '//trim(mqc_detstring_print(Determinants2%Strings%Alpha%vat([R_A_String],[0]),&
                    Determinants2%Strings%Beta%vat([R_B_String],[0]),Int(NBasis)))
                endIf
!                write(*,*) '-------------------------------------------'
!                write(*,*) ' L_Index: ',L_Index,'   R_Index:   ',R_Index 
!                write(*,*) '-------------------------------------------'
                Alpha_String_1 = Determinants%Strings%Alpha%vat([L_A_String],[1,NBit_Ints]) 
                Alpha_String_2 = Determinants2%Strings%Alpha%vat([R_A_String],[1,NBit_Ints]) 
                Beta_String_1 = Determinants%Strings%Beta%vat([L_B_String],[1,NBit_Ints]) 
                Beta_String_2 = Determinants2%Strings%Beta%vat([R_B_String],[1,NBit_Ints]) 
                If(S2Flag) then
                  Call CI_Hamiltonian%put(S2_Mat_Elem(IOut,IPrint,NBasis, &
                    Alpha_String_1,Beta_String_1,Alpha_String_2,Beta_String_2, &
                    MO_Core_Ham),L_Index,R_Index,SymmStr)
                Else
                  If(Present(MO_Core_Ham).and.Present(MO_ERIs)) then
                    Call CI_Hamiltonian%put(Slater_Condon(IOut,IPrint,NBasis, &
                      Alpha_String_1,Beta_String_1,Alpha_String_2,Beta_String_2, &
                      MO_Core_Ham,MO_ERIs),L_Index,R_Index,SymmStr)
                  ElseIf(Present(MO_Core_Ham).and..not.Present(MO_ERIs)) then
                    Call CI_Hamiltonian%put(Slater_Condon(IOut,IPrint,NBasis, &
                      Alpha_String_1,Beta_String_1,Alpha_String_2,Beta_String_2, &
                      MO_Core_Ham),L_Index,R_Index,SymmStr)
                  ElseIf(.not.Present(MO_Core_Ham).and.Present(MO_ERIs)) then
                    Call CI_Hamiltonian%put(Slater_Condon(IOut,IPrint,NBasis, &
                      Alpha_String_1,Beta_String_1,Alpha_String_2,Beta_String_2, &
                      TwoPartInts=MO_ERIs),L_Index,R_Index,SymmStr)
                  Else
                    Call MQC_Error('No integrals present in MQC_Build_CI_Hamiltonian')
                  EndIf
                EndIf

!                Counter = Counter + 1
!                If(iPrint.ge.1.and.Mod(100*Counter/(NDets*(NDets+1)/2),5).eq.0) then
!                  write(unit=bar(18:20),fmt="(i3)") 100*Counter/(NDets*(NDets+1)/2)
!                  do I = 1, 20*Counter/(NDets*(NDets+1)/2)
!                    bar(23+i:23+i) = '#'
!                  endDo
!                  write(unit=iOut,fmt='(a1,a1,a44)',advance='no') '+',char(13),bar
!                EndIf
              EndDo
            EndDo
          EndDo
        EndDo
        if(iPrint.ge.2) then
          write(6,'(A)') trim(detIndex1)
          if(.not.symmFlag) write(6,'(A)') trim(detIndex2)
        endIf
      case('ci') 
        L_Index = 0
        L_A_Start = 1
        Counter = 0
        Counter2 = 0
        Do L_A_Sub = 0, size(Determinants%NSubsAlpha)-1
          L_B_Start = 1
          Do L_B_Sub = 0, size(Determinants%NSubsBeta)-1
            If(.not.any(subs.eq.(L_A_Sub+L_B_Sub))) cycle 
            If (L_A_Sub.ne.0) L_A_Start = sum(Determinants%NSubsAlpha%vat(1,L_A_Sub))+1
            If (L_B_Sub.ne.0) L_B_Start = sum(Determinants%NSubsBeta%vat(1,L_B_Sub))+1
            If((L_A_Start.gt.NAlpha_Str1).or.(L_B_Start.gt.NBeta_Str1)) cycle
            L_A_End = L_A_Start + Determinants%NSubsAlpha%at(L_A_Sub+1) - 1
            L_B_End = L_B_Start + Determinants%NSubsBeta%at(L_B_Sub+1) - 1
            Do L_A_String = L_A_Start, L_A_End
              Do L_B_String = L_B_Start, L_B_End
                L_Index = L_Index + 1
                R_Index = 0
                R_A_Start = 1
                Do R_A_Sub = 0, size(Determinants2%NSubsAlpha)-1
                  R_B_Start = 1
                  Do R_B_Sub = 0, size(Determinants2%NSubsBeta)-1
                    If(.not.any(rightSubs.eq.(R_A_Sub+R_B_Sub))) cycle 
                    If (R_A_Sub.ne.0) R_A_Start = sum(Determinants2%NSubsAlpha%vat(1,R_A_Sub))+1
                    If (R_B_Sub.ne.0) R_B_Start = sum(Determinants2%NSubsBeta%vat(1,R_B_Sub))+1
                    If((R_A_Start.gt.NAlpha_Str2).or.(R_B_Start.gt.NBeta_Str2)) cycle
                    R_A_End = R_A_Start + Determinants2%NSubsAlpha%at(R_A_Sub+1) - 1
                    R_B_End = R_B_Start + Determinants2%NSubsBeta%at(R_B_Sub+1) - 1
                    Do R_A_String = R_A_Start, R_A_End
                      Do R_B_String = R_B_Start, R_B_End
                        R_Index = R_Index + 1
                        if(SymmFlag.and.L_Index.lt.R_Index) cycle
                        if(iPrint.ge.2.and.R_Index.eq.1) then
                          Counter = Counter + 1
                          if(L_Index.eq.1) then
                            if(SymmFlag) then
                              DetIndex1 = NEW_LINE('A')//' CI determinant index'//NEW_LINE('A')
                            else
                              DetIndex1 = NEW_LINE('A')//' Left CI determinant index'//NEW_LINE('A')
                            endIf
                          endIf
                          DetIndex1 = trim(DetIndex1)//NEW_LINE('A')//' '//trim(num2char(Counter,'I6'))//& 
                            ' = '//trim(mqc_detstring_print(Determinants%Strings%Alpha%vat([L_A_String],[0]),&
                            Determinants%Strings%Beta%vat([L_B_String],[0]),Int(NBasis)))
                        endIf
                        if(.not.SymmFlag.and.iPrint.ge.2.and.L_Index.eq.1) then
                          Counter2 = Counter2 + 1
                          if(R_Index.eq.1) DetIndex2 = NEW_LINE('A')//' Right CI determinant index'//NEW_LINE('A')
                          DetIndex2 = trim(DetIndex2)//NEW_LINE('A')//' '//trim(num2char(Counter2,'I6'))// & 
                            ' = '//trim(mqc_detstring_print(Determinants2%Strings%Alpha%vat([R_A_String],[0]),&
                            Determinants2%Strings%Beta%vat([R_B_String],[0]),Int(NBasis)))
                        endIf
!                        write(*,*) '-------------------------------------------'
!                        write(*,*) ' L_Index: ',L_Index,'   R_Index:   ',R_Index 
!                        write(*,*) '-------------------------------------------'
                        Alpha_String_1 = Determinants%Strings%Alpha%vat([L_A_String],[1,NBit_Ints]) 
                        Alpha_String_2 = Determinants2%Strings%Alpha%vat([R_A_String],[1,NBit_Ints]) 
                        Beta_String_1 = Determinants%Strings%Beta%vat([L_B_String],[1,NBit_Ints]) 
                        Beta_String_2 = Determinants2%Strings%Beta%vat([R_B_String],[1,NBit_Ints]) 
                        If(abs(R_A_Sub+R_B_Sub-L_A_Sub-L_B_Sub).gt.2) cycle 
                        If(S2Flag) then
                          Call CI_Hamiltonian%put(S2_Mat_Elem(IOut,IPrint,NBasis, &
                            Alpha_String_1,Beta_String_1,Alpha_String_2,Beta_String_2, &
                            MO_Core_Ham),L_Index,R_Index,SymmStr)
                        Else
                          If(Present(MO_Core_Ham).and.Present(MO_ERIs)) then
                            Call CI_Hamiltonian%put(Slater_Condon(IOut,IPrint,NBasis, &
                              Alpha_String_1,Beta_String_1,Alpha_String_2,Beta_String_2, &
                              MO_Core_Ham,MO_ERIs),L_Index,R_Index,SymmStr)
                          ElseIf(Present(MO_Core_Ham).and..not.Present(MO_ERIs)) then
                            Call CI_Hamiltonian%put(Slater_Condon(IOut,IPrint,NBasis, &
                              Alpha_String_1,Beta_String_1,Alpha_String_2,Beta_String_2, &
                              MO_Core_Ham),L_Index,R_Index,SymmStr)
                          ElseIf(.not.Present(MO_Core_Ham).and.Present(MO_ERIs)) then
                            Call CI_Hamiltonian%put(Slater_Condon(IOut,IPrint,NBasis, &
                              Alpha_String_1,Beta_String_1,Alpha_String_2,Beta_String_2, &
                              TwoPartInts=MO_ERIs),L_Index,R_Index,SymmStr)
                          Else
                            Call MQC_Error('No integrals present in MQC_Build_CI_Hamiltonian')
                          EndIf
                        EndIf
!                        Counter = Counter + 1
!                        If(iPrint.ge.1.and.Mod(100*Counter/(NDets*(NDets+1)/2),5).eq.0) then
!                          write(unit=bar(18:20),fmt="(i3)") 100*Counter/(NDets*(NDets+1)/2)
!                          do I = 1, 20*Counter/(NDets*(NDets+1)/2)
!                            bar(23+i:23+i) = '#'
!                          endDo
!                          write(unit=iOut,fmt='(a1,a1,a44)',advance='no') '+',char(13),bar
!                        EndIf
                      EndDo
                    EndDo
                  EndDo
                EndDo
              EndDo
            EndDo
          EndDo
        EndDo
        if(iPrint.ge.2) then
          write(6,'(A)') trim(detIndex1)
          if(.not.symmFlag) write(6,'(A)') trim(detIndex2)
        endIf
      case default
        call mqc_error_a('Unrecognized determinant storage type in mqc_build_ci_hamiltonian',&
          iOut,'Determinants%order',Determinants%order)
      end select
!
      If(iPrint.ge.2) write(iOut,*) NEW_LINE('a')

      End Subroutine MQC_Build_CI_Hamiltonian
!
!
!=====================================================================
!     
!     PROCEDURE MQC_GET_STRINGS_AT_INDEX
!
!>    \brief <b> MQC_GET_STRINGS_AT_INDEX is a subroutine that returns the strings
!>    corresponding to the index of matrices built by MQC_Build_CI_Hamiltonian</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_GET_STRINGS_AT_INDEX is a subroutine that returns the strings corresponding 
!>    to the index of matrices built by MQC_Build_CI_Hamiltonian.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64)
!>        The FORTRAN file to print to. 
!>    \endverbatim
!>    
!>    \param[in] IPrint
!>    \verbatim
!>        IPrint is Integer(kind=int64)
!>        The print level for the subroutine. 
!>    \endverbatim
!>    
!>    \param[in] Ind
!>    \verbatim
!>        Ind is Integer(kind=int64)
!>        The index from which determinant strings will be
!>        returned. If Ind is negative, the index will be
!>        counted from the end.
!>    \endverbatim
!>    
!>    \param[out] aString
!>    \verbatim
!>        aString is Integer(kind=int64),Dimension(:),Allocatable
!>        The alpha string at index Ind.
!>    \endverbatim
!>    
!>    \param[in] bString
!>    \verbatim
!>        bString is Integer(kind=int64),Dimension(:),Allocatable
!>        The beta string at index Ind.
!>    \endverbatim
!>    
!>    \param[in] NBasis
!>    \verbatim
!>        NBasis is Type(MQC_Scalar)
!>        The number of basis functions. 
!>    \endverbatim
!>    
!>    \param[in] Determinants
!>    \verbatim
!>        Determinants is Type(MQC_Determinant)
!>        The binary string occupation number vectors. If 
!>        opional argument Dets2 is provided Determinants
!>        variable will contain the row basis strings.
!>    \endverbatim
!>
!>    \param[in] Subs
!>    \verbatim
!>        Subs is Integer(kind=int64),Dimension(:),Optional
!>        Permitted substitutions to include in the Hamiltonian.
!>        If not present, all substitutions are included.
!>    \endverbatim
!>    
!     Authors:
!     ========
!>    \author L. M. Thompson
!>    \date 2024
!
      Subroutine MQC_Get_Strings_At_Index(IOut,IPrint,Ind,aString,bString,NBasis,Determinants, &
        Subs)
!
      Implicit None
      Integer(kind=int64),Intent(In)::IOut,IPrint,Ind
      Type(MQC_Scalar),Intent(In)::NBasis
      Type(MQC_Determinant),Intent(In)::Determinants
      Integer(kind=int64),Dimension(:),Intent(In),Optional::Subs
      Type(MQC_Vector),Intent(Out)::aString,bString

      Integer(kind=int64)::IndIn,I,J,K,NAlpha_Str1,NBeta_Str1,NDets1,L_A_String,L_B_String, &
        L_Index,L_A_Start,L_B_Start,L_A_End,L_B_End,L_A_Sub,L_B_Sub,NBit_Ints
!
!     Need to figure out dimensions of CI Hamiltonian, alpha strings and beta stings
      NAlpha_Str1 = MQC_Matrix_Rows(Determinants%Strings%Alpha)
      NBeta_Str1 = MQC_Matrix_Rows(Determinants%Strings%Beta)
      If(present(subs)) then
        NDets1 = 0
        Do I = 1, Size(Subs)
          Do J = 1, Size(Determinants%NSubsAlpha)
            Do K = 1, Size(Determinants%NSubsBeta)
              If((J-1)+(K-1).eq.Subs(i)) then
                NDets1 = NDets1 + Determinants%NSubsAlpha%at(J)*&
                  Determinants%NSubsBeta%at(K)
              EndIf
            EndDo
          EndDo
        EndDo
      Else
        NDets1 = NAlpha_Str1 * NBeta_Str1
      EndIf
      If(Ind.gt.NDets1) call mqc_error_i('Index requested is larger than dimension of matrix in &
        &MQC_Get_String_At_Index',6,'Ind',Ind,'NDets1',NDets1)
      If(Ind.eq.0) call mqc_error_i('Index zero requested in MQC_Get_String_At_Index',6,'Ind',Ind)
      If(Ind.lt.0) then
        IndIn = NDets1+Ind+1
      Else
        IndIn = Ind
      EndIf
!
      NBit_Ints = (NBasis/(Bit_Size(0)-1))+1 
!
      select case (Determinants%order)
      case ('lexical')
        Do L_A_String = 1, NAlpha_Str1  
          Do L_B_String = 1, NBeta_Str1  
            L_Index = 1+(L_B_String-1)*NAlpha_Str1+(L_A_String-1) 
            If(L_Index.eq.IndIn) then
              aString = Determinants%Strings%Alpha%vat([L_A_String],[1,NBit_Ints]) 
              bString = Determinants%Strings%Beta%vat([L_B_String],[1,NBit_Ints]) 
              If(iPrint.ge.3) then
                write(iOut,'(A)') ' Determinant at index '//trim(num2char(IndIn))//' = '//&
                  trim(mqc_detstring_print(Determinants%Strings%Alpha%vat([L_A_String],[0]),&
                  Determinants%Strings%Beta%vat([L_B_String],[0]),Int(NBasis)))
              EndIf
              Return
            EndIf
          EndDo
        EndDo
      case('ci') 
        L_Index = 0
        L_A_Start = 1
        Do L_A_Sub = 0, size(Determinants%NSubsAlpha)-1
          L_B_Start = 1
          Do L_B_Sub = 0, size(Determinants%NSubsBeta)-1
            If(.not.any(subs.eq.(L_A_Sub+L_B_Sub))) cycle 
            If (L_A_Sub.ne.0) L_A_Start = sum(Determinants%NSubsAlpha%vat(1,L_A_Sub))+1
            If (L_B_Sub.ne.0) L_B_Start = sum(Determinants%NSubsBeta%vat(1,L_B_Sub))+1
            If((L_A_Start.gt.NAlpha_Str1).or.(L_B_Start.gt.NBeta_Str1)) cycle
            L_A_End = L_A_Start + Determinants%NSubsAlpha%at(L_A_Sub+1) - 1
            L_B_End = L_B_Start + Determinants%NSubsBeta%at(L_B_Sub+1) - 1
            Do L_A_String = L_A_Start, L_A_End
              Do L_B_String = L_B_Start, L_B_End
                L_Index = L_Index + 1
                If(L_Index.eq.IndIn) then
                  aString = Determinants%Strings%Alpha%vat([L_A_String],[1,NBit_Ints]) 
                  bString = Determinants%Strings%Beta%vat([L_B_String],[1,NBit_Ints]) 
                  If(iPrint.ge.3) then
                    write(iOut,'(A)') ' Determinant at index '//trim(num2char(IndIn))//' = '//&
                      trim(mqc_detstring_print(Determinants%Strings%Alpha%vat([L_A_String],[0]),&
                      Determinants%Strings%Beta%vat([L_B_String],[0]),Int(NBasis)))
                  EndIf
                  Return
                EndIf
              EndDo
            EndDo
          EndDo
        EndDo
      case default
        call mqc_error_a('Unrecognized determinant storage type in mqc_get_strings_at_index',&
          iOut,'Determinants%order',Determinants%order)
      end select
!
      End Subroutine MQC_Get_Strings_At_Index 
!
!
!=====================================================================
!     
!     PROCEDURE Get_One_Gamma_Matrix
      function get_one_gamma_matrix(iOut,iPrint,nBasisIn,determinants,ci_amplitudes,UHF,nCoreIn,nOrbsIn,Subs) Result(onePDMint)
!
!     This function returns the CI one-particle density matrix 
!
!     Variable Declarations...
!
      Implicit None
      integer(kind=int64),intent(in)::iOut,iPrint
      Type(MQC_Scalar),Intent(In)::NBasisIn
      Type(MQC_Determinant),Intent(In)::Determinants
      type(mqc_vector),intent(in)::ci_amplitudes
      logical,intent(in)::UHF
      integer(kind=int64),optional,intent(in)::nCoreIn,nOrbsIn
      Integer(kind=int64),Dimension(:),Optional::Subs
      type(mqc_scf_integral)::onePDMint
      type(mqc_matrix)::onePDMalpha,onePDMbeta
      integer(kind=int64)::L_A_String,L_B_String,R_A_String,R_B_String
      Type(MQC_Scalar)::Sgn
      Integer(kind=int64)::NBasis,nCore,nOrbs,IPos,Det_Diff,ISgn,NBit_Ints,m,n,o,mm,nn,rdet,ldet, &
        Alpha_Diff_Cnt,Beta_Diff_Cnt,L_A_Start,L_B_Start,L_A_End,L_B_End,R_A_Start,R_B_Start,R_A_End,&
        R_B_End,L_A_Sub,L_B_Sub,R_A_Sub,R_B_Sub,L_A_Sub_End,L_B_Sub_End,R_A_Sub_End,R_B_Sub_End, &
        string_bits
      integer(kind=int64),dimension(2)::orb
      Integer(kind=int64),Dimension(:),Allocatable::Alpha_String_1,Alpha_String_2,Beta_String_1, &
        Beta_String_2,Alpha_Diff,Beta_Diff

      nBasis = nBasisIn
      if(present(nCoreIn)) then
        nCore = nCoreIn
      else
        nCore = 0
      endIf
      if(present(nOrbsIn)) then
        nOrbs = nOrbsIn
      else
        nOrbs = nBasis
      endIf
      nBit_Ints = (nBasis/(Bit_Size(0)-1))+1
      string_bits = ((nCore+nOrbs)/(Bit_Size(0)-1))+1
      Allocate(Alpha_String_1(string_bits),Alpha_String_2(string_bits),Beta_String_1(string_bits),Beta_String_2(string_bits))
      Allocate(Alpha_Diff(string_bits),Beta_Diff(string_bits))
      call onePDMalpha%init(nOrbs,nOrbs)
      call onePDMbeta%init(nOrbs,nOrbs)

      if(present(subs).and.determinants%order.ne.'ci') &
        call mqc_error('Truncated CI only implemented for CI order determinant strings')

      if(present(subs)) then
        L_A_Sub_End = size(Determinants%NSubsAlpha)-1
        L_B_Sub_End = size(Determinants%NSubsBeta)-1
        R_A_Sub_End = size(Determinants%NSubsAlpha)-1
        R_B_Sub_End = size(Determinants%NSubsBeta)-1
      else
        L_A_Sub_End = 0
        L_B_Sub_End = 0
        R_A_Sub_End = 0
        R_B_Sub_End = 0
      endIf
      ldet = 0
      L_A_Start = 1
      do L_A_Sub = 0, L_A_Sub_End
        L_B_Start =1
        do L_B_Sub = 0, L_B_Sub_End
          If(present(Subs)) then
            If(.not.any(subs.eq.(L_A_Sub+L_B_Sub))) cycle
            If (L_A_Sub.ne.0) L_A_Start = sum(Determinants%NSubsAlpha%vat(1,L_A_Sub))+1
            If (L_B_Sub.ne.0) L_B_Start = sum(Determinants%NSubsBeta%vat(1,L_B_Sub))+1
            If((L_A_Start.gt.determinants%NAlpStr).or.(L_B_Start.gt.determinants%NBetStr)) cycle
            L_A_End = L_A_Start + Determinants%NSubsAlpha%at(L_A_Sub+1) - 1
            L_B_End = L_B_Start + Determinants%NSubsBeta%at(L_B_Sub+1) - 1
          Else
            L_A_End = determinants%NAlpStr
            L_B_End = determinants%NBetStr
          EndIf
          Do L_A_String = L_A_Start, L_A_End 
            Do L_B_String = L_B_Start, L_B_End 
              If(present(subs)) then
                ldet = ldet + 1
              else
                ldet = 1+(L_B_String-1)*determinants%NAlpStr+(L_A_String-1)
              endIf
              rdet = 0
              R_A_Start = 1
              do R_A_Sub = 0, R_A_Sub_End
                R_B_Start =1
                do R_B_Sub = 0, R_B_Sub_End
                  If(present(Subs)) then
                    If(.not.any(subs.eq.(R_A_Sub+R_B_Sub))) cycle
                    If (R_A_Sub.ne.0) R_A_Start = sum(Determinants%NSubsAlpha%vat(1,R_A_Sub))+1
                    If (R_B_Sub.ne.0) R_B_Start = sum(Determinants%NSubsBeta%vat(1,R_B_Sub))+1
                    If((R_A_Start.gt.determinants%NAlpStr).or.(R_B_Start.gt.determinants%NBetStr)) cycle
                    R_A_End = R_A_Start + Determinants%NSubsAlpha%at(R_A_Sub+1) - 1
                    R_B_End = R_B_Start + Determinants%NSubsBeta%at(R_B_Sub+1) - 1
                  Else
                    R_A_End = determinants%NAlpStr
                    R_B_End = determinants%NBetStr
                  EndIf
                  Do R_A_String = R_A_Start, R_A_End 
                    Do R_B_String = R_B_Start, R_B_End 
                      If(present(subs)) then
                        rdet = rdet + 1
                      else
                        rdet = 1+(R_B_String-1)*determinants%NAlpStr+(R_A_String-1)
                      endIf
                      Alpha_String_1 = 0
                      Alpha_String_2 = 0
                      Beta_String_1 = 0
                      Beta_String_2 = 0
                      Alpha_Diff_Cnt = 0
                      Beta_Diff_Cnt = 0
                      Do m = 1,NBit_Ints
                        Alpha_String_1(m) = Determinants%Strings%Alpha%at(L_A_String,m)
                        Alpha_String_2(m) = Determinants%Strings%Alpha%at(R_A_String,m)
                        Beta_String_1(m) = Determinants%Strings%Beta%at(L_B_String,m)
                        Beta_String_2(m) = Determinants%Strings%Beta%at(R_B_String,m)
                        Alpha_Diff(m) = IEOR(Alpha_String_1(m),Alpha_String_2(m))
                        Alpha_Diff_Cnt = Alpha_Diff_Cnt + PopCnt(Alpha_Diff(m))
                        Beta_Diff(m) = IEOR(Beta_String_1(m),Beta_String_2(m))
                        Beta_Diff_Cnt = Beta_Diff_Cnt + PopCnt(Beta_Diff(m))
                      endDo
                      Det_Diff = Alpha_Diff_Cnt/2 + Beta_Diff_Cnt/2

                      select case (det_diff)
                      case(0)
                        Do IPos = nCore, nCore+nOrbs-1
                          m = IPos/Bit_Size(0)+1
                          n = Mod(IPos,Bit_Size(0))
                          If(BTest(Alpha_String_1(m),n).and.BTest(Beta_String_1(m),n)) then
                            call onePDMalpha%put(onePDMalpha%at(iPos+1-nCore,iPos+1-nCore)+ &
                              conjg(ci_amplitudes%at(ldet))*ci_amplitudes%at(rdet), &
                              iPos+1-nCore,iPos+1-nCore)
                            call onePDMbeta%put(onePDMbeta%at(iPos+1-nCore,iPos+1-nCore)+ &
                              conjg(ci_amplitudes%at(ldet))*ci_amplitudes%at(rdet), &
                              iPos+1-nCore,iPos+1-nCore)
                          elseIf(BTest(Alpha_String_1(m),n)) then
                            call onePDMalpha%put(onePDMalpha%at(iPos+1-nCore,iPos+1-nCore)+ &
                              conjg(ci_amplitudes%at(ldet))*ci_amplitudes%at(rdet), &
                              iPos+1-nCore,iPos+1-nCore)
                          elseIf(BTest(Beta_String_1(m),n)) then
                            call onePDMbeta%put(onePDMbeta%at(iPos+1-nCore,iPos+1-nCore)+ &
                              conjg(ci_amplitudes%at(ldet))*ci_amplitudes%at(rdet), &
                              iPos+1-nCore,iPos+1-nCore)
                          endIf
                        endDo
                      case(1)
                        Do IPos = nCore, nCore+nOrbs-1
                          m = IPos/Bit_Size(0)+1
                          n = Mod(IPos,Bit_Size(0))
                          If(Alpha_Diff_Cnt.eq.0) then
                            If(BTest(Beta_Diff(m),n)) then
                              If(BTest(Beta_String_1(m),n)) then
                                orb(1) = IPos+1-nCore
                              endif
                              If(BTest(Beta_String_2(m),n)) then
                                orb(2) = IPos+1-nCore
                              endif
                            endIf
                          elseIf(Beta_Diff_Cnt.eq.0) then
                            If(BTest(Alpha_Diff(m),n)) then
                              If(BTest(Alpha_String_1(m),n)) then
                                orb(1) = IPos+1-nCore
                              endif
                              If(BTest(Alpha_String_2(m),n)) then
                                orb(2) = IPos+1-nCore
                              endif
                            endIf
                          endIf
                        endDo
                        ISgn = 0
                        If(Orb(2)-Orb(1)-1.gt.0) then
                          m  = (Orb(1)-1+nCore)/Bit_Size(0)+1
                          n  = (Orb(2)-1+nCore)/Bit_Size(0)+1
                          Do o = m, n, -1
                            If(m.ne.o) then
                              mm = 1
                            Else
                              mm = Mod((Orb(1)-1+nCore),Bit_Size(0)) + 1
                            EndIf
                            If(n.ne.o) then
                              nn = Bit_Size(0)
                            Else
                              nn = Mod((Orb(2)-1+nCore),Bit_Size(0)) + 1
                            EndIf
                            If(Beta_Diff_Cnt.eq.0) then
                              ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(o),Alpha_String_2(o)),mm,nn-mm-1))
                            elseIf(Alpha_Diff_Cnt.eq.0) then
                              ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(o),Beta_String_2(o)),mm,nn-mm-1))
                            endIf
                          EndDo
                        elseIf(Orb(1)-Orb(2)-1.gt.0) then
                          m  = (Orb(2)-1+nCore)/Bit_Size(0)+1
                          n  = (Orb(1)-1+nCore)/Bit_Size(0)+1
                          Do o = m, n, -1
                            If(m.ne.o) then
                              mm = 1
                            Else
                              mm = Mod((Orb(2)-1+nCore),Bit_Size(0)) + 1
                            EndIf
                            If(n.ne.o) then
                              nn = Bit_Size(0)
                            Else
                              nn = Mod((Orb(1)-1+nCore),Bit_Size(0)) + 1
                            EndIf
                            If(Beta_Diff_Cnt.eq.0) then
                              ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(o),Alpha_String_2(o)),mm,nn-mm-1))
                            elseIf(Alpha_Diff_Cnt.eq.0) then
                              ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(o),Beta_String_2(o)),mm,nn-mm-1))
                            endIf
                          EndDo
                        EndIf
                        Sgn = (-1)**ISgn

                        if(Beta_Diff_Cnt.eq.0) then
                          call onePDMalpha%put(onePDMalpha%at(orb(1),orb(2))+sgn*&
                            conjg(ci_amplitudes%at(ldet))*ci_amplitudes%at(rdet), &
                            orb(1),orb(2))
                        elseIf(Alpha_Diff_Cnt.eq.0) then
                          call onePDMbeta%put(onePDMbeta%at(orb(1),orb(2))+sgn*&
                            conjg(ci_amplitudes%at(ldet))*ci_amplitudes%at(rdet), &
                            orb(1),orb(2))
                        endIf

                      case(2:)
                        cycle
                      case default
                        call mqc_error_i('Confused about string differences',iOut,'det_diff',det_diff)
                      end select
                    EndDo
                  EndDo
                EndDo
              EndDo
            EndDo
          EndDo
        EndDo
      EndDo

      If(.not.UHF) then
        call mqc_integral_allocate(onePDMint,'mo density','space',onePDMalpha)  
      else
        call mqc_integral_allocate(onePDMint,'mo density','spin',onePDMalpha,onePDMbeta)
      endIf

      end function get_one_gamma_matrix
!
      End Module MQC_EST  
