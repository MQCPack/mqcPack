      Module MQC_EST  
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
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     Single Reference Wavefunctions...
!
!     MQC_SCF_Integral
      Type MQC_SCF_Integral
!       we should make Alpha, Beta etc. private once MQC_Gaussian discards old 
!       matrix file reading routine
        Type(MQC_Matrix)::Alpha,Beta,AlphaBeta,BetaAlpha
        Character(Len=64),Private::Array_Name ! See below 
        Character(Len=64),Private::Array_Type ! Space,Spin,General 
      Contains 
        Procedure, Public::print => mqc_print_integral
        Procedure, Private::hasAlpha => mqc_integral_has_alpha
        Procedure, Private::hasBeta => mqc_integral_has_beta
        Procedure, Private::hasAlphaBeta => mqc_integral_has_alphabeta
        Procedure, Private::hasBetaAlpha => mqc_integral_has_betaalpha
        Procedure, Private::type => mqc_integral_array_type
        Procedure, Private::blockSize => mqc_integral_dimension 
        Procedure, Public::getLabel => mqc_integral_array_name
        Procedure, Public::addLabel => mqc_integral_add_name
        Procedure, Public::getBlock => mqc_integral_output_block
      End Type
!
!     MQC_SCF_Eigenvalues
      Type MQC_SCF_Eigenvalues
!       we should make Alpha, Beta private once MQC_Gaussian discards old 
!       matrix file reading routine
        Type(MQC_Vector)::Alpha,Beta
        Character(Len=64),Private::Array_Name ! See below 
        Character(Len=64),Private::Array_Type ! Space,Spin,General 
      Contains 
        Procedure, Public::print => mqc_print_eigenvalues
        Procedure, Private::hasAlpha => mqc_eigenvalues_has_alpha
        Procedure, Private::hasBeta => mqc_eigenvalues_has_beta
        Procedure, Private::type => mqc_eigenvalues_array_type
        Procedure, Private::blockSize => mqc_eigenvalues_dimension 
        Procedure, Public::getLabel => mqc_eigenvalues_array_name
        Procedure, Public::addLabel => mqc_eigenvalues_add_name
        Procedure, Public::getBlock => mqc_eigenvalues_output_block
      End Type
!
!     MQC_Wavefunction
      Type MQC_Wavefunction
        Type(MQC_SCF_Integral)::MO_Coefficients
        Type(MQC_SCF_Eigenvalues)::MO_Energies
        Type(MQC_SCF_Eigenvalues)::MO_Symmetries
        Type(MQC_SCF_Integral)::Core_Hamiltonian
        Type(MQC_SCF_Integral)::Fock_Matrix
        Type(MQC_SCF_Integral)::Density_Matrix
        Type(MQC_SCF_Integral)::Overlap_Matrix
        Type(MQC_Scalar)::NAlpha,NBeta,NElectrons,NBasis,Charge,Multiplicity
        Character(Len=256)::Basis,Symmetry,WF_Type
        Logical::WF_Complex
      Contains 
        Procedure, Public::print => mqc_print_wavefunction
      End Type MQC_Wavefunction
!
!
!     Post-SCF Wavefunctions...
!
!     Parent Type
      Type,Extends(MQC_Wavefunction)::MQC_PSCF_Wavefunction
        Integer::NCore,NVal,NActive
        Type(MQC_Matrix)::PSCF_Amplitudes
        Type(MQC_Vector)::PSCF_Energies
      End Type MQC_PSCF_Wavefunction
!
!    
!     Determinants....
!
!     MQC_Determinant_String
      Type MQC_Determinant_String
        Type(MQC_Matrix)::Alpha,Beta
      End Type MQC_Determinant_String
!
!     Parent Type
      Type MQC_Determinant
        Type(MQC_Determinant_String)::Strings
        Character(Len=64)::Order
      End Type MQC_Determinant
!
!    
!     Two electron integrals....
!
!     Parent Type
      Type MQC_TwoERIs
        Type(MQC_R4Tensor)::TwoERIs
        Logical::AO,UHF
        Character(Len=64)::Integral_Type,Storage_Type
      End Type MQC_TwoERIs
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
      Interface MQC_Print
        Module Procedure MQC_Print_Wavefunction
        Module Procedure MQC_Print_Integral
        Module Procedure MQC_Print_Eigenvalues
      End Interface
!
      Interface MatMul
        Module Procedure MQC_Integral_Matrix_Multiply
        Module Procedure MQC_Matrix_Integral_Multiply
        Module Procedure MQC_Integral_Integral_Multiply
      End Interface
!
      Interface Transpose
        Module Procedure MQC_Integral_Transpose
      End Interface
!
      Interface Dagger
        Module Procedure MQC_Integral_Conjugate_Transpose
      End Interface

!
      Interface MQC_Matrix_UndoSpinBlockGHF
        Module Procedure MQC_Matrix_UndoSpinBlockGHF_Eigenvalues
        Module Procedure MQC_Matrix_UndoSpinBlockGHF_Integral
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
      Interface Assignment (=)
        Module Procedure MQC_Integral_Output_Array
        Module Procedure MQC_Eigenvalues_Output_Array
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
      subroutine mqc_print_wavefunction(wavefunction,iOut,label)
!
      implicit none
      class(mqc_wavefunction)::wavefunction
      integer,intent(in)::iOut
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

        call wavefunction%overlap_matrix%print(iOut,'Overlap matrix',.true.,.false.)
        call wavefunction%core_hamiltonian%print(iOut,'Core Hamiltonian matrix',.true.,.false.)
        call wavefunction%mo_energies%print(iOut,'Orbital Energies',.true.,.false.)
        call wavefunction%mo_coefficients%print(iOut,'Molecular orbital coefficients',.true.,.false.)
        call wavefunction%density_matrix%print(iOut,'Density matrix',.true.,.false.)
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
      subroutine mqc_print_integral(integral,iOut,header, &
          blank_at_top,blank_at_bottom)
!
      implicit none
      class(mqc_scf_integral)::integral
      integer,intent(in)::iOut
      character(len=*),intent(in)::header
      character(len=64)::arrayType
      logical,intent(In),optional::blank_at_top,blank_at_bottom
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
        if(integral%hasAlpha()) call integral%alpha%print(iout,'')
      elseif(integral%array_type.eq.'spin') then
        if(integral%hasAlpha()) call integral%alpha%print(iout,'Alpha Array')
        if(integral%hasBeta()) call integral%beta%print(iout,'Beta Array')
      elseif(integral%array_type.eq.'general') then
        if(integral%hasAlpha()) call integral%alpha%print(iout,'Alpha-Alpha Block')
        if(integral%hasBeta()) call integral%beta%print(iout,'Beta-Beta Block')
        if(integral%hasAlphaBeta()) call integral%alphabeta%print(iout,'Alpha-Beta Block')
        if(integral%hasBetaAlpha()) call integral%betaalpha%print(iout,'Beta-Alpha Block')
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
      subroutine mqc_print_eigenvalues(eigenvalues,iOut,header, &
          blank_at_top,blank_at_bottom)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      integer,intent(in)::iOut
      character(len=*),intent(in)::header
      character(len=64)::arrayType
      logical,intent(In),optional::blank_at_top,blank_at_bottom
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
        if(eigenvalues%hasAlpha()) call eigenvalues%alpha%print(iout,'')
      elseif(eigenvalues%array_type.eq.'spin'.or.eigenvalues%array_type.eq.'general') then
        if(eigenvalues%hasAlpha()) call eigenvalues%alpha%print(iout,'Alpha Array')
        if(eigenvalues%hasBeta()) call eigenvalues%beta%print(iout,'Beta Array')
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
!     PROCEDURE MQC_Integral_Has_Alpha
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
!     PROCEDURE MQC_Integral_Array_Type
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
!     PROCEDURE MQC_Integral_Array_Name    
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
      function mqc_integral_dimension(integral,label) result(dimBlock)
!
      implicit none
      class(mqc_scf_integral),intent(in)::integral
      integer::dimBlock
      Character(Len=*),intent(in)::label
      Character(Len=64)::myLabel
!
      call string_change_case(label,'l',myLabel)
      select case (myLabel)
      case('alpha') 
        if(.not.integral%hasAlpha()) then
          dimBlock = 0
        else
          dimBlock = mqc_matrix_columns(integral%alpha)
        endIf
      case('beta')
        if(.not.integral%hasBeta()) then
          if(integral%type().eq.'space') then
            dimBlock = mqc_matrix_columns(integral%alpha)
          else
            dimBlock = 0
          endIf
        else
          dimBlock = mqc_matrix_columns(integral%beta)
        endIf
      case default
        call mqc_error_A('label not valid in mqc_integrals_dimension', 6, &
             'myLabel', myLabel )
      end select 
!
      end function mqc_integral_dimension 
!
!
!     PROCEDURE MQC_Eigenvaues_Dimension
      function mqc_eigenvalues_dimension(eigenvalues,label) result(dimBlock)
!
      implicit none
      class(mqc_scf_eigenvalues),intent(in)::eigenvalues
      integer::dimBlock
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
!     PROCEDURE MQC_Integral_Allocate
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
      if(present(alpha)) integral%alpha = alpha
      if(arrayType.eq.'spin') then
        if(present(beta)) integral%beta = beta
      elseIf(arrayType.eq.'general') then
        if(present(beta)) integral%beta = beta
!       add consistency checker for block sizes
        if(present(alphaBeta)) integral%alphaBeta = alphaBeta
        if(present(betaAlpha)) integral%betaAlpha = betaAlpha
      endIf
!
      end subroutine mqc_integral_allocate  
!
!
!     PROCEDURE MQC_Eigenvalues_Allocate
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
!     PROCEDURE MQC_Integral_Output_Block
      function mqc_integral_output_block(integral,blockName) result(matrixOut)
!
      implicit none
      class(mqc_scf_integral)::integral
      character(len=*),optional::blockName
      character(len=64)::myBlockName
      type(mqc_matrix)::matrixOut
      integer::nDimAlpha=0,nDimBeta=0,nDimTotal=0
!
      if(present(blockName)) then
        call string_change_case(blockName,'l',myBlockName)
      else
        myBlockName = 'full'
      endIf

      nDimAlpha = integral%blockSize('Alpha') 
      nDimBeta = integral%blockSize('Beta')
      nDimTotal = nDimAlpha + nDimBeta

      select case (myBlockName)
      case('full')
        if (integral%type().eq.'space') then
          call matrixOut%init(nDimTotal,nDimTotal)
          if (integral%hasAlpha()) then
            call matrixOut%mput(integral%alpha,[1,nDimAlpha],[1,nDimAlpha]) 
            call matrixOut%mput(integral%alpha,[nDimAlpha+1,nDimTotal],[nDimAlpha+1,nDimTotal]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlpha()', integral%hasAlpha() )
          endIf
        elseIf (integral%type().eq.'spin') then
          call matrixOut%init(nDimTotal,nDimTotal)
          if (integral%hasAlpha()) then
            call matrixOut%mput(integral%alpha,[1,nDimAlpha],[1,nDimAlpha]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlpha()', integral%hasAlpha() )
          endIf
          if (integral%hasBeta()) then
            call matrixOut%mput(integral%beta,[nDimAlpha+1,nDimTotal],[nDimAlpha+1,nDimTotal]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasBeta()', integral%hasBeta() )
          endIf
        elseIf (integral%type().eq.'general') then
          call matrixOut%init(nDimTotal,nDimTotal)
          if (integral%hasAlpha()) then
            call matrixOut%mput(integral%alpha,[1,nDimAlpha],[1,nDimAlpha]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlpha()', integral%hasAlpha() )
          endIf
          if (integral%hasBeta()) then
            call matrixOut%mput(integral%beta,[nDimAlpha+1,nDimTotal],[nDimAlpha+1,nDimTotal]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasBeta()', integral%hasBeta() )
          endIf
          if (integral%hasAlphaBeta()) then
            call matrixOut%mput(integral%alphaBeta,[nDimAlpha+1,nDimTotal],[1,nDimAlpha]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
                 'integral%hasAlphaBeta()', integral%hasAlphaBeta() )
          endIf
          if (integral%hasBetaAlpha()) then
            call matrixOut%mput(integral%betaAlpha,[1,nDimAlpha],[nDimAlpha+1,nDimTotal]) 
          else
            call mqc_error_L('block does not exist in mqc_integral_output_block', 6, &
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
            call matrixOut%init(NDimAlpha,NDimAlpha)
          elseIf (integral%type().eq.'spin') then
            call matrixOut%init(NDimAlpha,NDimBeta)
          endIf
        endIf
      case('beta-alpha')
        if (integral%hasBetaAlpha()) then
          matrixOut = integral%betaAlpha
        else
          if (integral%type().eq.'space') then
            call matrixOut%init(NDimAlpha,NDimAlpha)
          elseIf (integral%type().eq.'spin') then
            call matrixOut%init(NDimBeta,NDimAlpha)
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
!     PROCEDURE MQC_eigenvalues_Output_Block
      function mqc_eigenvalues_output_block(eigenvalues,blockName) result(vectorOut)
!
      implicit none
      class(mqc_scf_eigenvalues)::eigenvalues
      character(len=*),optional::blockName
      character(len=64)::myBlockName
      type(mqc_vector)::vectorOut
      integer::nDimAlpha=0,nDimBeta=0,nDimTotal=0
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
          call vectorOut%init(nDimAlpha)
          if (eigenvalues%hasAlpha()) then
            call vectorOut%vput(eigenvalues%alpha,1) 
          else
            call mqc_error_L('block does not exist in mqc_eigenvalues_output_block', 6, &
                 'eigenvalues%hasAlpha()', eigenvalues%hasAlpha() )
          endIf
        elseIf ((eigenvalues%type().eq.'spin').or.(eigenvalues%type().eq.'general')) then
          call vectorOut%init(nDimTotal)
          if (eigenvalues%hasAlpha()) then
            call vectorOut%vput(eigenvalues%alpha,1) 
          else
            call mqc_error_L('block does not exist in mqc_eigenvalues_output_block', 6, &
                 'eigenvalues%hasAlpha()', eigenvalues%hasAlpha() )
          endIf
          if (eigenvalues%hasBeta()) then
            call vectorOut%vput(eigenvalues%beta,nDimAlpha+1)
          else
            call mqc_error_L('block does not exist in mqc_eigenvalues_output_block', 6, &
                 'eigenvalues%hasBeta()', eigenvalues%hasBeta() )
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
!     PROCEDURE MQC_Integral_Output_Array
      subroutine mqc_integral_output_array(matrixOut,integralIn) 
!
      implicit none
      class(mqc_scf_integral),intent(in)::integralIn
      type(mqc_matrix),intent(inOut)::matrixOut
      integer::nDimAlpha=0,nDimBeta=0,nDimTotal=0
!
      nDimAlpha = integralIn%blockSize('Alpha') 
      nDimBeta = integralIn%blockSize('Beta')
      nDimTotal = nDimAlpha + nDimBeta

      select case (integralIn%type())
      case('space')
        call matrixOut%init(nDimAlpha,nDimAlpha)
        if (integralIn%hasAlpha()) then
          call matrixOut%mput(integralIn%alpha,[1,nDimAlpha],[1,nDimAlpha]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasAlpha()', integralIn%hasAlpha() )
        endIf
      case('spin')
        call matrixOut%init(nDimTotal,nDimTotal)
        if (integralIn%hasAlpha()) then
          call matrixOut%mput(integralIn%alpha,[1,nDimAlpha],[1,nDimAlpha]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasAlpha()', integralIn%hasAlpha() )
        endIf
        if (integralIn%hasBeta()) then
          call matrixOut%mput(integralIn%beta,[nDimAlpha+1,nDimTotal],[nDimAlpha+1,nDimTotal]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasBeta()', integralIn%hasBeta() )
        endIf
      case ('general')
        call matrixOut%init(nDimTotal,nDimTotal)
        if (integralIn%hasAlpha()) then
          call matrixOut%mput(integralIn%alpha,[1,nDimAlpha],[1,nDimAlpha]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasAlpha()', integralIn%hasAlpha() )
        endIf
        if (integralIn%hasBeta()) then
          call matrixOut%mput(integralIn%beta,[nDimAlpha+1,nDimTotal],[nDimAlpha+1,nDimTotal]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasBeta()', integralIn%hasBeta() )
        endIf
        if (integralIn%hasAlphaBeta()) then
          call matrixOut%mput(integralIn%alphaBeta,[nDimAlpha+1,nDimTotal],[1,nDimAlpha]) 
        else
          call mqc_error_L('block does not exist in mqc_integral_output_array', 6, &
               'integralIn%hasAlphaBeta()', integralIn%hasAlphaBeta() )
        endIf
        if (integralIn%hasBetaAlpha()) then
          call matrixOut%mput(integralIn%betaAlpha,[1,nDimAlpha],[nDimAlpha+1,nDimTotal]) 
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
      subroutine mqc_eigenvalues_output_array(vectorOut,eigenvaluesIn) 
!
      implicit none
      class(mqc_scf_eigenvalues),intent(in)::eigenvaluesIn
      type(mqc_vector),intent(inOut)::vectorOut
      integer::nDimAlpha=0,nDimBeta=0,nDimTotal=0
!
      nDimAlpha = eigenvaluesIn%blockSize('Alpha') 
      nDimBeta = eigenvaluesIn%blockSize('Beta')
      nDimTotal = nDimAlpha + nDimBeta

      select case (eigenvaluesIn%type())
      case('space')
        call vectorOut%init(nDimAlpha)
        if (eigenvaluesIn%hasAlpha()) then
          call vectorOut%vput(eigenvaluesIn%alpha,1) 
        else
          call mqc_error_L('block does not exist in mqc_eigenvalues_output_array', 6, &
               'eigenvaluesIn%hasAlpha()', eigenvaluesIn%hasAlpha() )
        endIf
      case('spin','general')
        call vectorOut%init(nDimTotal)
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
!     PROCEDURE MQC_Integral_Matrix_Multiply
      function mqc_integral_matrix_multiply(integralA,matrixB,label) result(integralOut)
!
!     This function multiplies a mqc integral type with an mqc matrix.
!
!     -L. M. Thompson, 2017.
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralA
      type(mqc_matrix),intent(in)::matrixB
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      integer::rows,nBasisAlpha,nBasisBeta,nBasisTotal
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
        nBasisAlpha = integralA%blockSize('Alpha')
        nBasisBeta = integralA%blockSize('Beta')
        nBasisTotal = nBasisAlpha + nBasisBeta
        if((rows.eq.nBasisAlpha).and.(nBasisAlpha.eq.nBasisBeta)) then
          doOffDiag = .False.
        elseIf(rows.eq.nBasisTotal) then
          doOffDiag = .True.
        else
          call mqc_error_I('Integral and matrix are wrongly sized for multiplication', 6, &
               'rows', rows, &
               'nBasisAlpha', nBasisAlpha, &
               'nBasisBeta', nBasisBeta, &
               'nBasisTotal', nBasisTotal )
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
          tmpMatrixAlpha = integralA%alpha.dot.matrixB%mat([1,nBasisAlpha],[1,nBasisAlpha])
          tmpMatrixBeta = integralA%alpha.dot.matrixB%mat([nBasisAlpha+1,nBasisTotal], &
            [nBasisAlpha+1,nBasisTotal])
          tmpMatrixAlphaBeta = integralA%alpha.dot.matrixB%mat([nBasisAlpha+1,nBasisTotal], &
            [1,nBasisAlpha])
          tmpMatrixBetaAlpha = integralA%alpha.dot.matrixB%mat([1,nBasisAlpha], &
            [nBasisAlpha+1,nBasisTotal])
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
          tmpMatrixAlpha = integralA%alpha.dot.matrixB%mat([1,nBasisAlpha],[1,nBasisAlpha])
          tmpMatrixBeta = integralA%beta.dot.matrixB%mat([nBasisAlpha+1,nBasisTotal], &
            [nBasisAlpha+1,nBasisTotal])
          tmpMatrixAlphaBeta = integralA%beta.dot.matrixB%mat([nBasisAlpha+1,nBasisTotal], &
            [1,nBasisAlpha])
          tmpMatrixBetaAlpha = integralA%alpha.dot.matrixB%mat([1,nBasisAlpha], &
            [nBasisAlpha+1,nBasisTotal])
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
          tmpMatrixAlpha = (integralA%alpha.dot.matrixB%mat([1,nBasisAlpha],[1,nBasisAlpha])) + &
            (integralA%betaAlpha.dot.matrixB%mat([nBasisAlpha+1,nBasisTotal],[1,nBasisAlpha]))
          tmpMatrixBeta = (integralA%beta.dot.matrixB%mat([nBasisAlpha+1,nBasisTotal],[nBasisAlpha+1,nBasisTotal])) + &
            (integralA%alphaBeta.dot.matrixB%mat([1,nBasisAlpha],[nBasisAlpha+1,nBasisTotal]))
          tmpMatrixAlphaBeta = (integralA%beta.dot.matrixB%mat([nBasisAlpha+1,nBasisTotal],[1,nBasisAlpha])) + &
            (integralA%alphaBeta.dot.matrixB%mat([1,nBasisAlpha],[1,nBasisAlpha]))
          tmpMatrixBetaAlpha = (integralA%alpha.dot.matrixB%mat([1,nBasisAlpha],[nBasisAlpha+1,nBasisTotal])) + &
            (integralA%betaAlpha.dot.matrixB%mat([nBasisAlpha+1,nBasisTotal],[nBasisAlpha+1,nBasisTotal])) 
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
      function mqc_matrix_integral_multiply(matrixA,integralB,label) result(integralOut)
!
!     This function multiplies a mqc integral type with an mqc matrix.
!
!     -L. M. Thompson, 2017.
!
      implicit none
      type(mqc_scf_integral),intent(in)::integralB
      type(mqc_matrix),intent(in)::matrixA
      Character(Len=*),optional,intent(in)::label
      type(mqc_scf_integral)::integralOut
      integer::columns,nBasisAlpha,nBasisBeta,nBasisTotal
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
        nBasisAlpha = integralB%blockSize('Alpha')
        nBasisBeta = integralB%blockSize('Beta')
        nBasisTotal = nBasisAlpha + nBasisBeta
        if((columns.eq.nBasisAlpha).and.(nBasisAlpha.eq.nBasisBeta)) then
          doOffDiag = .False.
        elseIf(columns.eq.nBasisTotal) then
          doOffDiag = .True.
        else
          call mqc_error_I('Integral and matrix are wrongly sized for multiplication', 6, &
               'columns', columns, &
               'nBasisAlpha', nBasisAlpha, &
               'nBasisBeta', nBasisBeta, &
               'nBasisTotal', nBasisTotal )
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
          tmpMatrixAlpha = matrixA%mat([1,nBasisAlpha],[1,nBasisAlpha]).dot.integralB%alpha
          tmpMatrixBeta = matrixA%mat([nBasisAlpha+1,nBasisTotal],[nBasisAlpha+1,nBasisTotal]) &
            .dot.integralB%alpha
          tmpMatrixAlphaBeta = matrixA%mat([nBasisAlpha+1,nBasisTotal],[1,nBasisAlpha]).dot. &
            integralB%alpha
          tmpMatrixBetaAlpha = matrixA%mat([1,nBasisAlpha],[nBasisAlpha+1,nBasisTotal]).dot. &
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
          tmpMatrixAlpha = matrixA%mat([1,nBasisAlpha],[1,nBasisAlpha]).dot.integralB%alpha
          tmpMatrixBeta = matrixA%mat([nBasisAlpha+1,nBasisTotal],[nBasisAlpha+1,nBasisTotal]) &
            .dot.integralB%beta
          tmpMatrixAlphaBeta = matrixA%mat([nBasisAlpha+1,nBasisTotal],[1,nBasisAlpha]).dot. &
            integralB%alpha
          tmpMatrixBetaAlpha = matrixA%mat([1,nBasisAlpha],[nBasisAlpha+1,nBasisTotal]).dot. &
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
          tmpMatrixAlpha = (matrixA%mat([1,nBasisAlpha],[1,nBasisAlpha]).dot.integralB%alpha) + &
            (matrixA%mat([1,nBasisAlpha],[nBasisAlpha+1,nBasisTotal]).dot.integralB%alphaBeta)
          tmpMatrixBeta = (matrixA%mat([nBasisAlpha+1,nBasisTotal],[nBasisAlpha+1,nBasisTotal]).dot.integralB%beta) + &
            (matrixA%mat([nBasisAlpha+1,nBasisTotal],[1,nBasisAlpha]).dot.integralB%betaAlpha)
          tmpMatrixAlphaBeta = (matrixA%mat([1,nBasisAlpha],[nBasisAlpha+1,nBasisTotal]).dot.integralB%alpha) + &
            (matrixA%mat([nBasisAlpha+1,nBasisTotal],[nBasisAlpha+1,nBasisTotal]).dot.integralB%alphaBeta)
          tmpMatrixBetaAlpha = (matrixA%mat([1,nBasisAlpha],[1,nBasisAlpha]).dot.integralB%betaAlpha) + &
            (matrixA%mat([1,nBasisAlpha],[nBasisAlpha+1,nBasisTotal]).dot.integralB%beta) 
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
!     PROCEDURE MQC_Integral_Integral_Multiply
      function mqc_integral_integral_multiply(integralA,integralB,label) result(integralOut)
!
!     This function multiplies a mqc integral type with an mqc integral.
!
!     -L. M. Thompson, 2017.
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
          call mqc_error_A('Unknown integral type in mqc_integral_matrix_multiply', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case('spin')
        select case (integralB%array_type)
        case('space')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%beta.dot.integralB%alpha
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha,tmpMatrixBeta)
        case('spin')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%beta.dot.integralB%beta
          call mqc_integral_allocate(integralOut,myLabel,'space',tmpMatrixAlpha,tmpMatrixBeta)
        case('general')
          tmpMatrixAlpha = integralA%alpha.dot.integralB%alpha
          tmpMatrixBeta = integralA%beta.dot.integralB%beta
          tmpMatrixAlphaBeta = integralA%beta.dot.integralB%alphaBeta
          tmpMatrixBetaAlpha = integralA%alpha.dot.integralB%betaAlpha
          call mqc_integral_allocate(integralOut,myLabel,'general',tmpMatrixAlpha,tmpMatrixBeta, &
            tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
        case default
          call mqc_error_A('Unknown integral type in mqc_integral_matrix_multiply', 6, &
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
          call mqc_error_A('Unknown integral type in mqc_integral_matrix_multiply', 6, &
               'integralB%array_type', integralB%array_type )
        end select
      case default
        call mqc_error_A('Unknown integral type in mqc_integral_matrix_multiply', 6, &
             'integralA%array_type', integralA%array_type )
      end select
!
      end function mqc_integral_integral_multiply
!
!
!     PROCEDURE MQC_Integral_Transpose
      function mqc_integral_transpose(integral,label) result(integralOut)
!
!     This function transposes a mqc integral type.
!
!     -L. M. Thompson, 2017.
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
      function mqc_integral_conjugate_transpose(integral,label) result(integralOut)
!
!     This function conjguate transposes a mqc integral type.
!
!     -L. M. Thompson, 2017.
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
!=====================================================================
!
!PROCEDURE MQC_Matrix_SpinBlockGHF
      subroutine mqc_matrix_spinBlockGHF(array)
!
!     This subroutine takes a GHF spin-interleaved spin array (such as typically 
!     generated by e.g. Gaussian) and returns it in spin blocked form i.e. the spinor 
!     basis is reordered so that alpha coefficients are before beta coefficients.
!
!     If MOs are passed, the input matrix columns ae reordered to 'spin block' the 
!     orbitals. As the orbitals can no longer be defined as being alpha or beta spin,
!     we use the definition that if a GHF calculation gives a UHF solution, the input
!     matrix would be reordered such that the off diagonal blocks are zero. 
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
!
      implicit none
      integer::i,j,k,rows,columns
      class(*),intent(inOut)::array
      type(mqc_matrix)::tmpMatrix
      type(mqc_vector)::tmpVector
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
        if(mqc_matrix_test_symmetric(array)) call mqc_matrix_full2Symm(array)
      class default
        call mqc_error_I('unrecognised array type in mqc_matrix_spinBlockGHF', 6)
      end select 
!
      return
      end subroutine mqc_matrix_spinBlockGHF


!=====================================================================
!
!PROCEDURE MQC_Matrix_UndoSpinBlockGHF_Eigenvalues
      subroutine mqc_matrix_undoSpinBlockGHF_EigenValues(eigenvaluesIn,vectorOut)
!
!     This subroutine takes a spin-blocked est eigenvalues object and 
!     returns it as a GHF spin-interleaved mqc vector (such as typically 
!     generated by e.g. Gaussian). 
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
!
      implicit none
      integer::i,j,nAlpha,nBeta,nElectrons
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


!=====================================================================
!
!PROCEDURE MQC_Matrix_UndoSpinBlockGHF_Integral   
      subroutine mqc_matrix_undoSpinBlockGHF_integral(integralIn,matrixOut)
!
!     This subroutine takes a spin-blocked est integral object and 
!     returns it as a GHF spin-interleaved mqc matrix (such as typically 
!     generated by e.g. Gaussian). 
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
!
      implicit none
      integer::i,j,k,rows,columns,nAlpha,nBeta,nElectrons
      type(mqc_scf_integral),intent(in)::integralIn
      type(mqc_matrix),intent(out)::matrixOut
      type(mqc_matrix)::tmpMatrix
!
!
!     Do the work...
!
      nAlpha = integralIn%blockSize('alpha') 
      nBeta = integralIn%blockSize('beta') 
      nElectrons = nAlpha+nBeta
      call matrixOut%init(nElectrons,nElectrons)
      tmpMatrix = integralIn%getBlock('full')
      j = 1
      do i = 1,nAlpha
        call matrixOut%vput(tmpMatrix%vat([i],[0]),[j],[0])
        j = j+2
      endDo
      j = 2
      do i = nAlpha+1,nElectrons 
        call matrixOut%vput(tmpMatrix%vat([i],[0]),[j],[0])
        j = j+2
      endDo
      tmpMatrix = matrixOut
      j = 1
      do i = 1,nAlpha
        call matrixOut%vput(tmpMatrix%vat([0],[i]),[0],[j])
        j = j+2
      endDo
      j = 2
      do i = nAlpha+1,nElectrons 
        call matrixOut%vput(tmpMatrix%vat([0],[i]),[0],[j])
        j = j+2
      endDo
!
      return
      end subroutine mqc_matrix_undoSpinBlockGHF_integral    


!=====================================================================
!
!     POST-SCF ROUTINES
!
!=====================================================================
!     
!     PROCEDURE GEN_DET_STR 
! 
      Subroutine Gen_Det_Str(IOut,IPrint,NBasisIn,NAlphaIn,NBetaIn, &
          Determinants)
!
!     This subroutine generates a list of alpha and beta strings in binary  
!     notation. 
!
!     The index to the determinant strings follows lexical order in which the 
!     left most path has index zero (in the case where the right-most bit 
!     corresponds to orbital 1 and the left-most bit corresponds to orbital 
!     NBasis - the HF determinant is the smallest integer number). 
!     Constructing the graph uses the following for the vertex and arc values:
!     Yo(e,o)=0, Y1(e,o)=x(e,o+1), x(e,o)=x(e+1,o+1)+x(e,o+1) 
!
!     The strings are also ordered in reverse lexical order where the right most 
!     path has index zero (in the case that the left-most bit corresponds to 
!     orbital 1 and the right-most bit corresponds to orbital NBasis - the HF 
!     determinant is the largest interger number)
!     Constructing the graph uses the following for the vertex and arc values:
!     Yo(e,o)=0, Y1(e,o)=x(e+1,o), x(e,o)=x(e,o-1)+x(e-1,o-1) 
!
!     While the routine is set up to deal with arbitrary number of orbitals, the 
!     maximum number of orbitals is currently limited to the machine word size.
!
!     Lee Thompson, 2016.
!
!     Variable Declarations...
!
      Implicit None

      Integer::IOut,IPrint,NBasis,NAlpha,NBeta,NAlpha_Str,NBeta_Str,NBit_Ints,IOrb,&
        NElec,NEMax,NEMin,I
      Integer,Dimension(:,:),Allocatable::Alpha_Strings,Beta_Strings
      Type(MQC_Scalar)::NAlphaIn,NBetaIn,NBasisIn
      Type(MQC_Determinant)::Determinants
      Type(MQC_LinkedList),Pointer::String_List_1=>Null(),String_List_2=>Null(),&
        String_Node_1=>Null(),String_Node_2=>Null()
      Integer::New_Value
      Integer,Allocatable::Returned_Value
      Logical::Last
!      
 1001 Format(1x,I7,2x,B64)
 1050 Format( A )
!
!     Initialize Arrays
!
      NAlpha = NAlphaIn 
      NBeta = NBetaIn
      NBasis = NBasisIn
      NAlpha_Str = Bin_Coeff(NBasis,NAlpha)
      NBeta_Str = Bin_Coeff(NBasis,NBeta)
      If(.not.Allocated(Alpha_Strings).and..not.Allocated(Beta_Strings)) then
        NBit_Ints = (NBasis/(Bit_Size(0)-1))+1
        If(NBit_Ints.gt.1) then
          Call MQC_Error_I('Determinant generator limited to 1 integer', 6, &
               'NBit_Ints', NBit_Ints )
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
      If(IPrint.ge.2) then
        Write(IOut,1050) "Alpha Strings"
        Do I = 1,NAlpha_Str
          Write(IOut,1001) I,Alpha_Strings(I,1)
        EndDo
        Write(IOut,1050) "Beta Strings"
        Do I = 1,NBeta_Str
          Write(IOut,1001) I,Beta_Strings(I,1)
        EndDo
      EndIf

      Determinants%Strings%Alpha = Alpha_Strings
      Determinants%Strings%Beta = Beta_Strings
      Determinants%Order = 'Lexical'

      Deallocate(Alpha_Strings,Beta_Strings)

      End Subroutine Gen_Det_Str
!
!=====================================================================
!     
!     PROCEDURE SLATER_CONDON
      Function Slater_Condon(IOut,IPrint,NBasisIn,Determinants,L_A_String, &
       L_B_String,R_A_String,R_B_String,Core_Hamiltonian,ERIs,UHF) &
       Result(MatEl)
!
!     This function returns the CI Hamiltonian matrix element value using 
!     Slater-Condon rules for a given alpha and beta string combination.
!     If UHF flag is set to true then Slater Condon rules are returned
!     accounting for different spacial occupation of alpha and beta orbitals but
!     the ERIs must be dimensioned to twice that of rhf case.
!
!     Variable Declarations...
!
      Implicit None
      Integer,Intent(In)::IOut,IPrint,L_A_String,L_B_String,R_A_String,R_B_String
      Logical,Intent(In)::UHF
      Type(MQC_Scalar),Intent(In)::NBasisIn
      Type(MQC_SCF_Integral),Intent(In)::Core_Hamiltonian
      Type(MQC_Determinant),Intent(In)::Determinants
      Type(MQC_TwoERIs),Intent(In)::ERIs
      Type(MQC_Scalar)::MatEl,Sgn
      Integer::NBasis,IPos,JPos,IDiff,Det_Diff,ISgn,NBit_Ints,I,J,K, &
        II,JJ,Alpha_Diff_Cnt,Beta_Diff_Cnt
      Real::ERI1,ERI2,Zero=0.0
      Integer,Dimension(4)::Orbs,Spin,Det
      Integer,Dimension(:),Allocatable::Alpha_String_1,Alpha_String_2,Beta_String_1, &
        Beta_String_2,Alpha_Diff,Beta_Diff
      
 1050 Format( A )
      NBasis = NBasisIn
      NBit_Ints = (NBasis/(Bit_Size(0)-1))+1 
      Allocate(Alpha_String_1(NBit_Ints),Alpha_String_2(NBit_Ints),Beta_String_1(NBit_Ints),Beta_String_2(NBit_Ints))
      Alpha_String_1 = Determinants%Strings%Alpha%vat([L_A_String],[1,NBit_Ints]) 
      Alpha_String_2 = Determinants%Strings%Alpha%vat([R_A_String],[1,NBit_Ints]) 
      Beta_String_1 = Determinants%Strings%Beta%vat([L_B_String],[1,NBit_Ints]) 
      Beta_String_2 = Determinants%Strings%Beta%vat([R_B_String],[1,NBit_Ints]) 
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
!       I am going to comment the first logical block - the rest are different
!       perutations
!       First we need to determine the relevant orbital, spin and determinant
          IDiff = 1
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
!            Write(IOut,*) 'Orbs:',IPos,' I:',I,' J:',J
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

          If(Det(1).eq.Det(2).and.Det(3).eq.Det(4)) then
!            Write(IOut,*) 'Entering case 1'
!            
!           This section computes the Coulomb contribution to the matrix
!           element
            If(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
!              
!             This section detemines number of permutations of creation
!             operators and sets sign accordingly
!             It got quite a bit longer after generalizing max number of
!             orbitals
!             Probably want to make some functions for these then later
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(2).eq.0) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
!              
!             This section computes the value of the matrix element
!
              If((.not.UHF).or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1),Orbs(3),Orbs(2),Orbs(4))
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1),Orbs(3),Orbs(2)+NBasis,Orbs(4)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(3)+NBasis,Orbs(2),Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(3)+NBasis,Orbs(2)+NBasis,Orbs(4)+NBasis)
              EndIf
!              
!             This section changes sign accounting for the ordering of alpha and
!             beta strings
              If(Spin(1).gt.Spin(2)) ERI1 = -ERI1
              If(Spin(3).gt.Spin(4)) ERI1 = -ERI1
              If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!              Write(IOut,*) 'Setting ERI1 to',ERI1
            Else
!              
!             Sets Coulomb contribution to zero if 2ERI is necessarily so by
!             spin order in the integral
!              Write(IOut,*) 'Setting ERI1 to zero'
              ERI1 = Zero
            EndIf
            If(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(2).eq.0) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
!
              If((.not.UHF).or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1),Orbs(4),Orbs(2),Orbs(3))
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1),Orbs(4),Orbs(2)+NBasis,Orbs(3)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(4)+NBasis,Orbs(2),Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(4)+NBasis,Orbs(2)+NBasis,Orbs(3)+NBasis)
              EndIf
!
              If(Spin(1).gt.Spin(2)) ERI2 = -ERI2
              If(Spin(3).gt.Spin(4)) ERI2 = -ERI2
              If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!              Write(IOut,*) 'Setting ERI2 to',ERI2
            Else 
!              Write(IOut,*) 'Setting ERI2 to zero'
              ERI2 = Zero
            EndIf
!
          ElseIf(Det(1).eq.Det(3).and.Det(2).eq.Det(4)) then
!            Write(IOut,*) 'Entering case 2'
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(3).eq.0) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
              If(.not.UHF.or.(Spin(1).eq.0.and.Spin(3).eq.0)) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1),Orbs(2),Orbs(3),Orbs(4))
              ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1),Orbs(2),Orbs(3)+NBasis,Orbs(4)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(2)+NBasis,Orbs(3),Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(2)+NBasis,Orbs(3)+NBasis,Orbs(4)+NBasis)
              EndIf
              If(Spin(1).gt.Spin(3)) ERI1 = -ERI1
              If(Spin(2).gt.Spin(4)) ERI1 = -ERI1
              If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!              Write(IOut,*) 'Setting ERI1 to',ERI1
            Else
!              Write(IOut,*) 'Setting ERI1 to zero'
              ERI1 = Zero
            EndIf
            If(Spin(1).eq.Spin(4).and.Spin(3).eq.Spin(2)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(3).eq.0) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
              If(.not.UHF.or.(Spin(1).eq.0.and.Spin(3).eq.0)) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1),Orbs(4),Orbs(3),Orbs(2))
              ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1),Orbs(4),Orbs(3)+NBasis,Orbs(2)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(4)+NBasis,Orbs(3),Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(4)+NBasis,Orbs(3)+NBasis,Orbs(2)+NBasis)
              EndIf
              If(Spin(1).gt.Spin(3)) ERI2 = -ERI2
              If(Spin(2).gt.Spin(4)) ERI2 = -ERI2
              If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!              Write(IOut,*) 'Setting ERI2 to',ERI2
            Else 
!              Write(IOut,*) 'Setting ERI2 to zero'
              ERI2 = Zero
            EndIf
!
          ElseIf(Det(1).eq.Det(4).and.Det(2).eq.Det(3)) then
!            Write(IOut,*) 'Entering case 3'
            If(Spin(1).eq.Spin(2).and.Spin(4).eq.Spin(3)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(4).eq.0) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
              If(.not.UHF.or.(Spin(1).eq.0.and.Spin(4).eq.0)) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1),Orbs(2),Orbs(4),Orbs(3))
              ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1),Orbs(2),Orbs(4)+NBasis,Orbs(3)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(2)+NBasis,Orbs(4),Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
                ERI1 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(2)+NBasis,Orbs(4)+NBasis,Orbs(3)+NBasis)
              EndIf
              If(Spin(1).gt.Spin(4)) ERI1 = -ERI1
              If(Spin(2).gt.Spin(3)) ERI1 = -ERI1
              If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!              Write(IOut,*) 'Setting ERI1 to',ERI1
            Else
!              Write(IOut,*) 'Setting ERI1 to zero'
              ERI1 = Zero
            EndIf
            If(Spin(1).eq.Spin(3).and.Spin(4).eq.Spin(2)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(4).eq.0) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
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
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
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
              If(.not.UHF.or.(Spin(1).eq.0.and.Spin(4).eq.0)) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1),Orbs(3),Orbs(4),Orbs(2))
              ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1),Orbs(3),Orbs(4)+NBasis,Orbs(2)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(3)+NBasis,Orbs(4),Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
                ERI2 = ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(3)+NBasis,Orbs(4)+NBasis,Orbs(2)+NBasis)
              EndIf
              If(Spin(1).gt.Spin(4)) ERI2 = -ERI2
              If(Spin(2).gt.Spin(3)) ERI2 = -ERI2
              If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!              Write(IOut,*) 'Setting ERI2 to',ERI2
            Else 
!              Write(IOut,*) 'Setting ERI2 to zero'
              ERI2 = Zero
            EndIf
          Else
            Write(IOut,1050) 'Slater Condon has been handed confusing determinant info'
          EndIf
!
          MatEl = ERI1 - ERI2
!          Write(IOut,*) 'Slater Condon Final:',MatEl
!          Deallocate(Orbs,Spin,Det)
          Return
!
        Case(1)
          IDiff = 1
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
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
!          Write(IOut,*)'Orb 1:',Orbs(1),' Orb 2:',Orbs(2)
!          Write(IOut,*)'Spin 1:',Spin(1),' Spin 2:',Spin(2)
!          Write(IOut,*)'Det 1:',Det(1),' Det 2:',Det(2)
!
          ISgn = 0
          If(Spin(1).eq.0.and.Spin(2).eq.0) then
            If(Orbs(2)-Orbs(1)-1.gt.0) then 
              I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
              J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
              I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
              J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
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
!          Write(IOut,*)'ISgn:',ISgn
!
          MatEl = Zero
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
!            Write(IOut,*) 'Alpha IPos:', IPos
            If(BTest(Alpha_String_1(I),J)) then
              If(.not.UHF.and.Spin(1).eq.Spin(2)) then
!                Write(IOut,*) 'ALPHA: UHF False and 1 and 2 same spin, adding:',ISgn*ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1)
                MatEl = MatEl + Sgn*ERIs%TwoERIs%at(Orbs(1),Orbs(2),IPos+1,IPos+1)
                If(Spin(1).eq.0.and.Spin(2).eq.0) then
!                  Write(IOut,*) 'ALPHA: UHF False and 1 and 2 both alpha, subtracting:',ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
                  MatEl = MatEl - Sgn*ERIs%TwoERIs%at(Orbs(1),IPos+1,IPos+1,Orbs(2))
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.0) then
!                Write(IOut,*) 'ALPHA: UHF True and 1 and 2 both alpha, adding:', &
!                  ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1) - ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
                MatEl = MatEl + Sgn*ERIs%TwoERIs%at(Orbs(1),Orbs(2),IPos+1,IPos+1) - &
                  Sgn*ERIs%TwoERIs%at(Orbs(1),IPos+1,IPos+1,Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!                Write(IOut,*) 'ALPHA: UHF True and 1 and 2 both beta, adding:', ISgn*ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,IPos+1,IPos+1)
                MatEl = MatEl + Sgn*ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(2)+NBasis,IPos+1,IPos+1)
              EndIf
            EndIf
          EndDo
!
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0))
            If(BTest(Beta_String_1(I),J)) then
              If(.not.UHF.and.Spin(1).eq.Spin(2)) then
!                Write(IOut,*) 'BETA: UHF False and 1 and 2 same spin, adding:',ISgn*ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1)
                MatEl = MatEl + Sgn*ERIs%TwoERIs%at(Orbs(1),Orbs(2),IPos+1,IPos+1)
                If(Spin(1).eq.1.and.Spin(2).eq.1) then
!                  Write(IOut,*) 'BETA: UHF False and 1 and 2 both beta, subtracting:',ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
                  MatEl = MatEl - Sgn*ERIs%TwoERIs%at(Orbs(1),IPos+1,IPos+1,Orbs(2))
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.0) then
!                Write(IOut,*) 'BETA: UHF True and 1 and 2 both alpha, adding:', ISgn*ERIs(Orbs(1),Orbs(2),IPos+1+NBasis,IPos+1+NBasis)
                MatEl = MatEl + Sgn*ERIs%TwoERIs%at(Orbs(1),Orbs(2),IPos+1+NBasis,IPos+1+NBasis) 
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!                Write(IOut,*) 'BETA: UHF True and 1 and 2 both beta, adding:', &
!                ISgn*ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,IPos+1+NBasis,IPos+1+NBasis) - &
!                  ISgn*ERIs(Orbs(1)+NBasis,IPos+1+NBasis,IPos+1+NBasis,Orbs(2)+NBasis)
                MatEl = MatEl + Sgn*ERIs%TwoERIs%at(Orbs(1)+NBasis,Orbs(2)+NBasis,IPos+1+NBasis,IPos+1+NBasis) - &
                  Sgn*ERIs%TwoERIs%at(Orbs(1)+NBasis,IPos+1+NBasis,IPos+1+NBasis,Orbs(2)+NBasis)
              EndIf
            EndIf
          EndDo
!
          If(.not.UHF.or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
!            Write(IOut,*) '1 and 2 both alpha, adding core:', ISgn*Core_Hamiltonian%Alpha(Orbs(1),Orbs(2))
            MatEl = MatEl + Sgn*Core_Hamiltonian%Alpha%at(Orbs(1),Orbs(2)) 
          ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!            Write(IOut,*) '1 and 2 both beta, adding core:', ISgn*Core_Hamiltonian%Beta(Orbs(1),Orbs(2))
            MatEl = MatEl + Sgn*Core_Hamiltonian%Beta%at(Orbs(1),Orbs(2))
          EndIf
          
!          Call MQC_Print(MatEl,IOut,'CI Matrix Element')

          Return
!
        Case(0)
          MatEl = Zero
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0))
            If(BTest(Alpha_String_1(I),J)) then
!              Write(IOut,*) 'Alpha:',IPos+1,' I:',I,' J:',J
!              Call MQC_Print(Core_Hamiltonian%Alpha%at(IPos+1,IPos+1),IOut,'(I|h|I)')
              MatEl = MatEl + Core_Hamiltonian%Alpha%at(IPos+1,IPos+1)
            EndIf
            If(BTest(Beta_String_1(I),J)) then
!              Write(IOut,*) 'Beta:',IPos+1,' I:',I,' J:',J
              If(UHF) then
!                Call MQC_Print(Core_Hamiltonian%Beta%at(IPos+1,IPos+1),IOut,'(I|h|I)')
                MatEl = MatEl + Core_Hamiltonian%Beta%at(IPos+1,IPos+1)
              Else
!                Call MQC_Print(Core_Hamiltonian%Alpha%at(IPos+1,IPos+1),IOut,'(I|h|I)')
                MatEl = MatEl + Core_Hamiltonian%Alpha%at(IPos+1,IPos+1)
              EndIf
            EndIf
          EndDo
!
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = IPos+1, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0)) 
              If(BTest(Alpha_String_1(I),J)) then
                If(BTest(Alpha_String_1(II),JJ)) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
!                  Call MQC_Print(ERIs%TwoERIs%at(IPos+1,IPos+1,JPos+1,JPos+1),IOut,'(II|JJ)')
!                  Call MQC_Print(ERIs%TwoERIs%at(IPos+1,JPos+1,JPos+1,IPos+1),IOut,'(IJ|JI)') 
                  MatEl = MatEl + ERIs%TwoERIs%at(IPos+1,IPos+1,JPos+1,JPos+1) - ERIs%TwoERIs%at(IPos+1,JPos+1,JPos+1,IPos+1)
                EndIf
              EndIf
            EndDo
          EndDo
!
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = IPos+1, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0))
              If(BTest(Beta_String_1(I),J)) then
                If(BTest(Beta_String_1(II),JJ)) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  If(UHF) then
!                    Call MQC_Print(ERIs%TwoERIs%at(IPos+1+NBasis,IPos+1+NBasis,JPos+1+NBasis,JPos+1+NBasis),IOut,'(II|JJ)')
!                    Call MQC_Print(ERIs%TwoERIs%at(IPos+1+NBasis,JPos+1+NBasis,JPos+1+NBasis,IPos+1+NBasis),IOut,'(IJ|JI)') 
                    MatEl = MatEl + ERIs%TwoERIs%at(IPos+1+NBasis,IPos+1+NBasis,JPos+1+NBasis,JPos+1+NBasis) - & 
                      ERIs%TwoERIs%at(IPos+1+NBasis,JPos+1+NBasis,JPos+1+NBasis,IPos+1+NBasis)
                  Else
!                    Call MQC_Print(ERIs%TwoERIs%at(IPos+1,IPos+1,JPos+1,JPos+1),IOut,'(II|JJ)')
!                    Call MQC_Print(ERIs%TwoERIs%at(IPos+1,JPos+1,JPos+1,IPos+1),IOut,'(IJ|JI)') 
                    MatEl = MatEl + ERIs%TwoERIs%at(IPos+1,IPos+1,JPos+1,JPos+1) - &
                      ERIs%TwoERIs%at(IPos+1,JPos+1,JPos+1,IPos+1)
                  EndIf
                EndIf
              EndIf
            EndDo
          EndDo
!
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = 0, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0))
              If(BTest(Alpha_String_1(I),J)) then
                If(BTest(Beta_String_1(II),JJ)) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
                  If(UHF) then
!                    Call MQC_Print(ERIs%TwoERIs%at(IPos+1,IPos+1,JPos+1+NBasis,JPos+1+NBasis),IOut,'(II|JJ)')
                    MatEl = MatEl + ERIs%TwoERIs%at(IPos+1,IPos+1,JPos+1+NBasis,JPos+1+NBasis)
                  Else
!                    Call MQC_Print(ERIs%TwoERIs%at(IPos+1,IPos+1,JPos+1,JPos+1),IOut,'(II|JJ)')
                    MatEl = MatEl + ERIs%TwoERIs%at(IPos+1,IPos+1,JPos+1,JPos+1)
                  EndIf
                EndIf
              EndIf
            EndDo
          EndDo

!          Call MQC_Print(MatEl,IOut,'CI Matrix Element')

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
      Subroutine TwoERI_Trans(IOut,IPrint,MO_Coeff,ERIs,MO_ERIs,UHF)
!
!     This subroutine transforms two-electron integrals to the MO basis 
!     If the UHF flag is true it computes the set of MO integrals accounting
!     for different spacial occupation of alpha and beta electrons. The order
!     in the final MO ERI array is (aa|aa), (aa|bb), (bb|aa), (bb|bb).
!
!             |' (aa|aa) | (bb|aa) '| 
!             |----------|----------|
!             |, (aa|bb) | (bb|bb) ,|
!
!     Variable Declarations...
!
      Implicit None
      Integer::IOut,IPrint,NBasis,P,Q,R,S,Mu,Nu,Lambda,Sigma,IErr
      Real::Zero=0.0d0
      Logical::DoON5,UHF
      Type(MQC_Matrix)::X,Y
      Type(MQC_TwoERIs),Intent(In)::ERIs
      Type(MQC_TwoERIs),Intent(Out)::MO_ERIs
      Type(MQC_SCF_Integral),Intent(In)::MO_Coeff
!
 1000 Format(1x,'(',I3,',',I3,'|',I3,',',I3,') = ',F15.8)
 2000 Format(1x,'Doing O(N**8) integral transformation algorithm')
 3000 Format(1x,'Doing O(N**5) integral transformation algorithm')
!
!     Initialize Arrays
!
      DoON5 = .True.

      NBasis = MQC_Matrix_Rows(MO_Coeff%Alpha) 
      If(.not.UHF) then
        Call MO_ERIs%TwoERIs%initialize(NBasis,NBasis,NBasis,NBasis)
      ElseIf(UHF) then
        Call MO_ERIs%TwoERIs%initialize(NBasis*2,NBasis*2,NBasis*2,NBasis*2)
      EndIf
      MO_ERIs%AO = .False.
      MO_ERIs%Storage_Type = 'Full'
      MO_ERIs%Integral_Type = 'Regular'

      If(.not.DoON5) then

        If(IPrint.ge.2) Write(IOut,2000) 
        Do P = 1, NBasis
          Do Q = 1, NBasis
            Do R = 1, NBasis
              Do S = 1, NBasis
                Do Mu = 1, NBasis
                  Do Nu = 1, NBasis
                    Do Lambda = 1, NBasis
                      Do Sigma = 1, NBasis
                        Call MO_ERIs%TwoERIs%put(MO_ERIs%TwoERIs%at(P,Q,R,S) + &
                          MO_Coeff%Alpha%at(Mu,P) * MO_Coeff%Alpha%at(Nu,Q) * &
                          ERIs%TwoERIs%at(Mu,Nu,Lambda,Sigma) * MO_Coeff%Alpha%at(Lambda,R) * &
                          MO_Coeff%Alpha%at(Sigma,S),P,Q,R,S)
                          If(UHF) then
                           Call MO_ERIs%TwoERIs%put(MO_ERIs%TwoERIs%at(P,Q,R+NBasis,S+NBasis) + &
                             MO_Coeff%Alpha%at(Mu,P) * MO_Coeff%Alpha%at(Nu,Q) * &
                             ERIs%TwoERIs%at(Mu,Nu,Lambda,Sigma) * &
                             MO_Coeff%Beta%at(Lambda,R) * MO_Coeff%Beta%at(Sigma,S), &
                             P,Q,R+NBasis,S+NBasis)
                           Call MO_ERIs%TwoERIs%put(MO_ERIs%TwoERIs%at(P+NBasis,Q+NBasis,R,S) + &
                             MO_Coeff%Beta%at(Mu,P) * MO_Coeff%Beta%at(Nu,Q) * &
                             ERIs%TwoERIs%at(Mu,Nu,Lambda,Sigma) * & 
                             MO_Coeff%Alpha%at(Lambda,R) * MO_Coeff%Alpha%at(Sigma,S), &
                             P+NBasis,Q+NBasis,R,S)
                           Call MO_ERIs%TwoERIs%put(MO_ERIs%TwoERIs%at(P+NBasis,Q+NBasis,R+NBasis,S+NBasis) + &
                             MO_Coeff%Beta%at(Mu,P) * MO_Coeff%Beta%at(Nu,Q) * &
                             ERIs%TwoERIs%at(Mu,Nu,Lambda,Sigma) * &
                             MO_Coeff%Beta%at(Lambda,R) * MO_Coeff%Beta%at(Sigma,S), &
                             P+NBasis,Q+NBasis,R+NBasis,S+NBasis)
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

        If(IPrint.ge.2) Write(IOut,3000) 

        Call X%initialize(NBasis,NBasis)
        Call Y%initialize(NBasis,NBasis)

        Do P = 1, NBasis
          Do Q = 1, NBasis
            Do R = 1, NBasis
              Do S = 1, NBasis
                Call X%put(ERIs%TwoERIs%at(P,Q,R,S),R,S) 
              EndDo
            EndDo
            Y = Transpose(MO_Coeff%Alpha).dot.X
            X = Y.dot.MO_Coeff%Alpha
            Do R = 1, NBasis
              Do S = 1, NBasis
                Call MO_ERIs%TwoERIs%put(X%at(R,S),P,Q,R,S)
              EndDo
            EndDo
          EndDo
        EndDo

        Do R = 1, NBasis
          Do S = 1, NBasis
            Do P = 1, NBasis
              Do Q = 1, NBasis 
                Call X%put(MO_ERIs%TwoERIs%at(P,Q,R,S),P,Q)
              EndDo
            EndDo
            Y = Transpose(MO_Coeff%Alpha).dot.X
            X = Y.dot.MO_Coeff%Alpha
            Do P = 1, NBasis
              Do Q = 1, NBasis
                Call MO_ERIs%TwoERIs%put(X%at(P,Q),P,Q,R,S)
              EndDo
            EndDo
          EndDo
        EndDo

        If(UHF) then
         
          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Call X%put(ERIs%TwoERIs%at(P,Q,R,S),R,S)
                EndDo
              EndDo
              Y = Transpose(MO_Coeff%Beta).dot.X
              X = Y.dot.MO_Coeff%Beta
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Call MO_ERIs%TwoERIs%put(X%at(R,S),P,Q,R+NBasis,S+NBasis)
                EndDo
              EndDo
            EndDo
          EndDo

          Do R = 1, NBasis
            Do S = 1, NBasis
              Do P = 1, NBasis
                Do Q = 1, NBasis 
                  Call X%put(MO_ERIs%TwoERIs%at(P,Q,R+NBasis,S+NBasis),P,Q)
                EndDo
              EndDo
              Y = Transpose(MO_Coeff%Alpha).dot.X
              X = Y.dot.MO_Coeff%Alpha
              Do P = 1, NBasis
                Do Q = 1, NBasis
                  Call MO_ERIs%TwoERIs%put(X%at(P,Q),P,Q,R+NBasis,S+NBasis)
                EndDo
              EndDo
            EndDo
          EndDo

          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Call X%put(ERIs%TwoERIs%at(P,Q,R,S),R,S) 
                EndDo
              EndDo
              Y = Transpose(MO_Coeff%Alpha).dot.X
              X = Y.dot.MO_Coeff%Alpha
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Call MO_ERIs%TwoERIs%put(X%at(R,S),P+NBasis,Q+NBasis,R,S)
                EndDo
              EndDo
            EndDo
          EndDo

          Do R = 1, NBasis
            Do S = 1, NBasis
              Do P = 1, NBasis
                Do Q = 1, NBasis 
                  Call X%put(MO_ERIs%TwoERIs%at(P+NBasis,Q+NBasis,R,S),P,Q)
                EndDo
              EndDo
              Y = Transpose(MO_Coeff%Beta).dot.X
              X = Y.dot.MO_Coeff%Beta
              Do P = 1, NBasis
                Do Q = 1, NBasis
                  Call MO_ERIs%TwoERIs%put(X%at(P,Q),P+NBasis,Q+NBasis,R,S)
                EndDo
              EndDo
            EndDo
          EndDo

          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Call X%put(ERIs%TwoERIs%at(P,Q,R,S),R,S) 
                EndDo
              EndDo
              Y = Transpose(MO_Coeff%Beta).dot.X
              X = Y.dot.MO_Coeff%Beta
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Call MO_ERIs%TwoERIs%put(X%at(R,S),P+NBasis,Q+NBasis,R+NBasis,S+NBasis)
                EndDo
              EndDo
            EndDo
          EndDo

          Do R = 1, NBasis
            Do S = 1, NBasis
              Do P = 1, NBasis
                Do Q = 1, NBasis 
                  Call X%put(MO_ERIs%TwoERIs%at(P+NBasis,Q+NBasis,R+NBasis,S+NBasis),P,Q)
                EndDo
              EndDo
              Y = Transpose(MO_Coeff%Beta).dot.X
              X = Y.dot.MO_Coeff%Beta
              Do P = 1, NBasis
                Do Q = 1, NBasis
                  Call MO_ERIs%TwoERIs%put(X%at(P,Q),P+NBasis,Q+NBasis,R+NBasis,S+NBasis)
                EndDo
              EndDo
            EndDo
          EndDo

        EndIf
    
      EndIf
      
      If(IPrint.ge.2) Call MQC_Print(MO_ERIs%TwoERIs,IOut,'Transformed MO 2ERIs')
!
      End Subroutine TwoERI_Trans
!
!=====================================================================
!!
!!     PROCEDURE S2_MAT_ELEM    
!! 
!      Function S2_Mat_Elem(IOut,IPrint,NBasis,Alpha_String_1,Beta_String_1, &
!      Alpha_String_2,Beta_String_2,MO_Overlap)
!!
!!     This function returns the CI S**2 matrix elements used for computing S**2
!!     values of CI vectors for a given alpha and beta string combination.
!!     The MO overlap matrix is required
!!
!!     Variable Declarations...
!!
!      Implicit None
!      Integer::IOut,IPrint,NBasis,IPos,JPos,IDiff,Det_Diff,NAlpha,NBeta, &
!        IOcc,JOcc,KOcc,LOcc,Mat_Sign,Alpha_Diff_Cnt,Beta_Diff_Cnt,NBit_Ints, &
!        I,J,II,JJ
!      Real::S2_Mat_Elem,Zero=0.0d0,Quarter=0.25d0,ABTerm,One=1.0d0
!      Integer,Dimension(4)::Orbs,Spin,Det
!      Integer,Dimension(:)::Alpha_String_1,Alpha_String_2,Beta_String_1, &
!        Beta_String_2
!      Integer,Dimension(:),Allocatable::Alpha_Diff,Beta_Diff
!      Real,Dimension(:,:),Allocatable::MO_Overlap
!!
!!      Write(IOut,*) 'alpha 1:'
!!      Write(IOut,'(B64)') Alpha_String_1
!!      Write(IOut,*) 'alpha 2:'
!!      Write(IOut,'(B64)') Alpha_String_2
!!      Write(IOut,*) 'alpha XOR results in slater condon'
!!      Write(IOut,'(B64)') IEOR(Alpha_String_1,Alpha_String_2)   
!!      Write(IOut,*) 'NDifa=', PopCnt(IEOR(Alpha_String_1,Alpha_String_2))
!!      Write(IOut,*)
!!      Write(IOut,*) 'beta 1:'
!!      Write(IOut,'(B64)') Beta_String_1
!!      Write(IOut,*) 'beta 2:'
!!      Write(IOut,'(B64)') Beta_String_2
!!      Write(IOut,*) 'beta XOR results in slater condon'
!!      Write(IOut,'(B64)') IEOR(Beta_String_1,Beta_String_2)   
!!      Write(IOut,*) 'NDifb=', PopCnt(IEOR(Beta_String_1,Beta_String_2))
!!      Write(IOut,*)
!!
!      NBit_Ints = (NBasis/Bit_Size(0))+1 
!      Allocate(Alpha_Diff(NBit_Ints),Beta_Diff(NBit_Ints))
!      Det_Diff = 0
!      Alpha_Diff_Cnt = 0
!      Beta_Diff_Cnt = 0
!      NAlpha = 0
!      NBeta = 0
!      Do I = 1,NBit_Ints
!        Alpha_Diff(I) = IEOR(Alpha_String_1(I),Alpha_String_2(I))
!!        Write(IOut,*) 'Alpha Diff',I,':'
!!        Write(IOut,'(B64)') Alpha_Diff(I)
!!        Write(IOut,*) '-------------'
!        Alpha_Diff_Cnt = Alpha_Diff_Cnt + PopCnt(Alpha_Diff(I)) 
!        NAlpha = NAlpha + PopCnt(Alpha_String_1(I))
!        Beta_Diff(I) = IEOR(Beta_String_1(I),Beta_String_2(I))
!!        Write(IOut,*) 'Beta Diff',I,':'
!!        Write(IOut,'(B64)') Beta_Diff(I)
!!        Write(IOut,*) '-------------'
!        Beta_Diff_Cnt = Beta_Diff_Cnt + PopCnt(Beta_Diff(I))
!        NBeta = NBeta + PopCnt(Beta_String_1(I))
!      EndDo
!!      Write(IOut,*)'Alpha_Diff_Cnt:',Alpha_Diff_Cnt,'Beta_Diff_Cnt:',Beta_Diff_Cnt
!      Det_Diff = Alpha_Diff_Cnt/2 + Beta_Diff_Cnt/2
!
!      If(Mod(Alpha_Diff_Cnt,2).ne.0.or.Mod(Beta_Diff_Cnt,2).ne.0) then
!        Write(IOut,*) "ERROR: S2_Mat_Elem has been handed spin non-conserving &
!        determinants"
!        Call Exit()
!      EndIf
!
!!      Write(IOut,*) 'Det_Diff:',Det_Diff
!      Select Case (Det_Diff)
!!
!        Case(3:)
!          S2_Mat_Elem = Zero 
!          Return
!!
!        Case(2)
!          IDiff = 1
!          Do IPos = 0, NBasis-1
!            I = NBit_Ints - IPos/Bit_Size(0)
!            J = Mod(IPos,Bit_Size(0)) 
!            If(BTest(Alpha_Diff(I),J).eq..True.) then
!              Orbs(IDiff) = IPos+1
!              Spin(IDiff) = 0
!              If(BTest(Alpha_String_1(I),J).eq..True.) then
!                Det(IDiff) = 1
!              Else 
!                Det(IDiff) = 2
!              EndIf
!              IDiff = IDiff + 1
!            EndIf
!            If(BTest(Beta_Diff(I),J).eq..True.) then
!              Orbs(IDiff) = IPos+1
!              Spin(IDiff) = 1
!              If(BTest(Beta_String_1(I),J).eq..True.) then
!                Det(IDiff) = 1
!              Else 
!                Det(IDiff) = 2
!              EndIf
!              IDiff = IDiff + 1
!            EndIf
!          EndDo
!!          Call Print_Vector(IOut,Orbs,'Orbs')
!!          Call Print_Vector(IOut,Spin,'Spin')
!!          Call Print_Vector(IOut,Det,'Det')
!!
!          IOcc = 0
!          Do IPos = Orbs(1)-1, 0, -1
!            I = NBit_Ints - IPos/Bit_Size(0)
!            J = Mod(IPos,Bit_Size(0)) 
!            If(Spin(1).eq.0) then
!              If(Det(1).eq.1) then
!                If(BTest(Alpha_String_1(I),J).eq..True.) IOcc = IOcc + 1
!              ElseIf(Det(1).eq.2) then
!                If(BTest(Alpha_String_2(I),J).eq..True.) IOcc = IOcc + 1
!              EndIf
!            ElseIf(Spin(1).eq.1) then
!              If(Det(1).eq.1) then
!                If(BTest(Beta_String_1(I),J).eq..True.) IOcc = IOcc + 1
!              ElseIf(Det(1).eq.2) then
!                If(BTest(Beta_String_2(I),J).eq..True.) IOcc = IOcc + 1
!              EndIf
!            EndIf
!          EndDo
!!          Write(IOut,*) 'IOcc:',IOcc
!          JOcc = 0
!          Do IPos = Orbs(2)-1, 0, -1
!            I = NBit_Ints - IPos/Bit_Size(0)
!            J = Mod(IPos,Bit_Size(0)) 
!            If(Spin(2).eq.0) then
!              If(Det(2).eq.1) then
!                If(BTest(Alpha_String_1(I),J).eq..True.) JOcc = JOcc + 1
!              ElseIf(Det(2).eq.2) then
!                If(BTest(Alpha_String_2(I),J).eq..True.) JOcc = JOcc + 1
!              EndIf
!            ElseIf(Spin(2).eq.1) then
!              If(Det(2).eq.1) then
!                If(BTest(Beta_String_1(I),J).eq..True.) JOcc = JOcc + 1
!              ElseIf(Det(2).eq.2) then
!                If(BTest(Beta_String_2(I),J).eq..True.) JOcc = JOcc + 1
!              EndIf
!            EndIf
!          EndDo
!!          Write(IOut,*) 'JOcc:',JOcc
!          KOcc = 0
!          Do IPos = Orbs(3)-1, 0, -1
!            I = NBit_Ints - IPos/Bit_Size(0)
!            J = Mod(IPos,Bit_Size(0)) 
!            If(Spin(3).eq.0) then
!              If(Det(3).eq.1) then
!                If(BTest(Alpha_String_1(I),J).eq..True.) KOcc = KOcc + 1
!              ElseIf(Det(3).eq.2) then
!                If(BTest(Alpha_String_2(I),J).eq..True.) KOcc = KOcc + 1
!              EndIf
!            ElseIf(Spin(3).eq.1) then
!              If(Det(3).eq.1) then
!                If(BTest(Beta_String_1(I),J).eq..True.) KOcc = KOcc + 1
!              ElseIf(Det(3).eq.2) then
!                If(BTest(Beta_String_2(I),J).eq..True.) KOcc = KOcc + 1
!              EndIf
!            EndIf
!          EndDo
!!          Write(IOut,*) 'KOcc:',KOcc
!          LOcc = 0
!          Do IPos = Orbs(4)-1, 0, -1
!            I = NBit_Ints - IPos/Bit_Size(0)
!            J = Mod(IPos,Bit_Size(0)) 
!            If(Spin(4).eq.0) then
!              If(Det(4).eq.1) then
!                If(BTest(Alpha_String_1(I),J).eq..True.) LOcc = LOcc + 1
!              ElseIf(Det(4).eq.2) then
!                If(BTest(Alpha_String_2(I),J).eq..True.) LOcc = LOcc + 1
!              EndIf
!            ElseIf(Spin(4).eq.1) then
!              If(Det(4).eq.1) then
!                If(BTest(Beta_String_1(I),J).eq..True.) LOcc = LOcc + 1
!              ElseIf(Det(4).eq.2) then
!                If(BTest(Beta_String_2(I),J).eq..True.) LOcc = LOcc + 1
!              EndIf
!            EndIf
!          EndDo
!!          Write(IOut,*) 'LOcc:',LOcc
!!          Mat_Sign = -1
!!          Mat_Sign = (-1)**(IOcc+JOcc+KOcc+LOcc-3)
!!          Write(IOut,*) 'Permutations:',(IOcc+JOcc+KOcc+LOcc-3)
!          Mat_Sign = (-1)**(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!!          Write(IOut,*) 'Permutations:',(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!!          Write(IOut,*) 'Mat_Sign:',Mat_Sign
!!
!          If(Det(1).eq.Det(2).and.Det(3).eq.Det(4)) then
!            If(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
!              If(Spin(1).eq.0.and.Spin(2).eq.1) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(4)+NBasis)*MO_Overlap(Orbs(2)+NBasis,Orbs(3))
!              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(2),Orbs(3)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(4))
!              Else
!!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
!                S2_Mat_Elem = Zero
!              EndIf
!            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
!              If(Spin(1).eq.0.and.Spin(2).eq.1) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(3)+NBasis)*MO_Overlap(Orbs(2)+NBasis,Orbs(4))
!              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(2),Orbs(4)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(3))
!              Else
!!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
!                S2_Mat_Elem = Zero
!              EndIf
!            Else
!!             This suggests that there are unbalanced spins between determinants 
!              S2_Mat_Elem = Zero
!            EndIf
!          ElseIf(Det(1).eq.Det(3).and.Det(2).eq.Det(4)) then
!            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
!              If(Spin(1).eq.0.and.Spin(3).eq.1) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(4)+NBasis)*MO_Overlap(Orbs(3)+NBasis,Orbs(2))
!              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(3),Orbs(2)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(4))
!              Else
!!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
!                S2_Mat_Elem = Zero
!              EndIf
!            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
!              If(Spin(1).eq.0.and.Spin(3).eq.1) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(2)+NBasis)*MO_Overlap(Orbs(3)+NBasis,Orbs(4))
!              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(3),Orbs(4)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(2))
!              Else
!!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
!                S2_Mat_Elem = Zero
!              EndIf
!            Else
!!             This suggests that there are unbalanced spins between determinants 
!              S2_Mat_Elem = Zero
!            EndIf
!          ElseIf(Det(1).eq.Det(4).and.Det(2).eq.Det(3)) then
!            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
!              If(Spin(1).eq.0.and.Spin(4).eq.1) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(3)+NBasis)*MO_Overlap(Orbs(4)+NBasis,Orbs(2))
!              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(4),Orbs(2)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(3))
!              Else
!!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
!                S2_Mat_Elem = Zero
!              EndIf
!            ElseIf(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
!              If(Spin(1).eq.0.and.Spin(4).eq.1) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(2)+NBasis)*MO_Overlap(Orbs(4)+NBasis,Orbs(3))
!              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
!                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(4),Orbs(3)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(2))
!              Else
!!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
!                S2_Mat_Elem = Zero
!              EndIf
!            Else
!!             This suggests that there are unbalanced spins between determinants 
!              S2_Mat_Elem = Zero
!            EndIf
!          EndIf
!!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem
!
!          Return
!
!        Case(1)
!          IDiff = 1
!!          Allocate(Orbs(2),Spin(2))
!          Do IPos = 0, NBasis-1
!            I = NBit_Ints - IPos/Bit_Size(0)
!            J = Mod(IPos,Bit_Size(0)) 
!            If(BTest(Alpha_Diff(I),J).eq..True.) then
!              Orbs(IDiff) = IPos+1
!              Spin(IDiff) = 0
!              If(BTest(Alpha_String_1(I),J).eq..True.) then
!                Det(IDiff) = 1
!              Else 
!                Det(IDiff) = 2
!              EndIf
!              IDiff = IDiff + 1
!            EndIf
!            If(BTest(Beta_Diff(I),J).eq..True.) then
!              Orbs(IDiff) = IPos+1
!              Spin(IDiff) = 1
!              If(BTest(Beta_String_1(I),J).eq..True.) then
!                Det(IDiff) = 1
!              Else 
!                Det(IDiff) = 2
!              EndIf
!              IDiff = IDiff + 1
!            EndIf
!          EndDo
!!          Write(IOut,*)'Orb 1:',Orbs(1),' Orb 2:',Orbs(2)
!!          Write(IOut,*)'Spin 1:',Spin(1),' Spin 2:',Spin(2)
!!          Write(IOut,*)'Det 1:',Det(1),' Det 2:',Det(2)
!!
!          S2_Mat_Elem = Zero 
!          If(Spin(1).ne.Spin(2)) then
!            S2_Mat_Elem = Zero
!!
!          ElseIf(Spin(1).eq.0) then
!!
!            IOcc = 0
!            Do IPos = Orbs(1)-1, 0, -1
!              I = NBit_Ints - IPos/Bit_Size(0)
!              J = Mod(IPos,Bit_Size(0)) 
!              If(Det(1).eq.1) then
!                If(BTest(Alpha_String_1(I),J).eq..True.) IOcc = IOcc + 1
!              ElseIf(Det(1).eq.2) then
!                If(BTest(Alpha_String_2(I),J).eq..True.) IOcc = IOcc + 1
!              EndIf
!            EndDo
!!            Write(IOut,*) 'IOcc:',IOcc
!            JOcc = 0
!            Do IPos = Orbs(2)-1, 0, -1
!              I = NBit_Ints - IPos/Bit_Size(0)
!              J = Mod(IPos,Bit_Size(0)) 
!              If(Det(2).eq.1) then
!                If(BTest(Alpha_String_1(I),J).eq..True.) JOcc = JOcc + 1
!              ElseIf(Det(2).eq.2) then
!                If(BTest(Alpha_String_2(I),J).eq..True.) JOcc = JOcc + 1
!              EndIf
!            EndDo
!!            Write(IOut,*) 'JOcc:',JOcc
!!
!            Do IPos = 0, NBasis-1
!              I = NBit_Ints - IPos/Bit_Size(0)
!              J = Mod(IPos,Bit_Size(0)) 
!              If(BTest(Beta_String_1(I),J).eq..True.) then
!                S2_Mat_Elem = S2_Mat_Elem + MO_Overlap(IPos+1+NBasis,Orbs(1)) * MO_Overlap(Orbs(2),IPos+1+NBasis)
!              EndIf
!            EndDo
!!            S2_Mat_Elem = - S2_Mat_Elem
!!            Write(IOut,*) 'Permutations:',(2*NAlpha+1-IOcc-JOcc)
!!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NAlpha+1-IOcc-JOcc)
!            S2_Mat_Elem = (-1)**(2*NAlpha+1-IOcc-JOcc) * S2_Mat_Elem
!!            S2_Mat_Elem = (-1)**(IOcc+JOcc-1) * S2_Mat_Elem
!!
!          ElseIf(Spin(1).eq.1) then
!
!            IOcc = 0
!            Do IPos = Orbs(1)-1, 0, -1
!              I = NBit_Ints - IPos/Bit_Size(0)
!              J = Mod(IPos,Bit_Size(0)) 
!              If(Det(1).eq.1) then
!                If(BTest(Beta_String_1(I),J).eq..True.) IOcc = IOcc + 1
!              ElseIf(Det(1).eq.2) then
!                If(BTest(Beta_String_2(I),J).eq..True.) IOcc = IOcc + 1
!              EndIf
!            EndDo
!!            Write(IOut,*) 'IOcc:',IOcc
!            JOcc = 0
!            Do IPos = Orbs(2)-1, 0, -1
!              I = NBit_Ints - IPos/Bit_Size(0)
!              J = Mod(IPos,Bit_Size(0)) 
!              If(Det(2).eq.1) then
!                If(BTest(Beta_String_1(I),J).eq..True.) JOcc = JOcc + 1
!              ElseIf(Det(2).eq.2) then
!                If(BTest(Beta_String_2(I),J).eq..True.) JOcc = JOcc + 1
!              EndIf
!            EndDo
!!            Write(IOut,*) 'JOcc:',JOcc
!!
!            Do IPos = 0, NBasis-1
!              I = NBit_Ints - IPos/Bit_Size(0)
!              J = Mod(IPos,Bit_Size(0)) 
!              If(BTest(Alpha_String_1(I),J).eq..True.) then
!                S2_Mat_Elem = S2_Mat_Elem + MO_Overlap(IPos+1,Orbs(1)+NBasis) * MO_Overlap(Orbs(2)+NBasis,IPos+1)
!              EndIf
!            EndDo
!!            S2_Mat_Elem = - S2_Mat_Elem
!!            Write(IOut,*) 'Permutations:',(2*NBeta+1-IOcc-JOcc)
!!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NBeta+1-IOcc-JOcc)
!            S2_Mat_Elem = (-1)**(2*NBeta+1-IOcc-JOcc) * S2_Mat_Elem
!!            S2_Mat_Elem = (-1)**(IOcc+JOcc-1) * S2_Mat_Elem
!
!          EndIf
!
!!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem
!
!          Return
!!
!        Case(0)
!          ABTerm = Zero
!          Do IPos = 0, NBasis-1
!            I = NBit_Ints - IPos/Bit_Size(0)
!            J = Mod(IPos,Bit_Size(0)) 
!            Do JPos = 0, NBasis-1
!              II = NBit_Ints - JPos/Bit_Size(0)
!              JJ = Mod(JPos,Bit_Size(0)) 
!              If((BTest(Alpha_String_2(I),J).eq..True.).and.(BTest(Beta_String_2(II),JJ).eq..True.)) then
!                ABTerm = ABTerm + MO_Overlap(IPos+1,JPos+1+NBasis)*MO_Overlap(JPos+1+NBasis,IPos+1) 
!              EndIf
!            EndDo
!          EndDo
!          S2_Mat_Elem = Quarter*((NAlpha-NBeta)**2+2*(NAlpha+NBeta)) - ABTerm
!!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem
!!
!          Return
!!
!      End Select
!!
!      End Function S2_Mat_Elem  
!!
!=====================================================================
!     
!     PROCEDURE MQC_BUILD_CI_HAMILTONIAN
      Subroutine MQC_Build_CI_Hamiltonian(IOut,IPrint,NBasis,Determinants, &
        MO_Core_Ham,MO_ERIs,UHF,CI_Hamiltonian)
!
!     This subroutine is used to build the CI Hamiltonian given MO
!     integrals and determinant strings.
!
!     Variable Declarations...
!
      Implicit None
      Integer,Intent(In)::IOut,IPrint
      Logical,Intent(In)::UHF
      Type(MQC_Scalar),Intent(In)::NBasis
      Type(MQC_TwoERIs),Intent(In)::MO_ERIs
      Type(MQC_SCF_Integral),Intent(In)::MO_Core_Ham
      Type(MQC_Determinant),Intent(In)::Determinants
      Type(MQC_Matrix),Intent(Out)::CI_Hamiltonian
      Integer::NAlpha_Str,NBeta_Str,NDets,L_A_String,L_B_String, &
        R_A_String,R_B_String,L_Index,R_Index
!
      NAlpha_Str = MQC_Matrix_Rows(Determinants%Strings%Alpha)
      NBeta_Str = MQC_Matrix_Rows(Determinants%Strings%Beta)
      NDets = NAlpha_Str * NBeta_Str 
!
      Call CI_Hamiltonian%initialize(NDets,NDets)
!
      Do L_A_String = 1, NAlpha_Str  
        Do L_B_String = 1, NBeta_Str  
          Do R_A_String = 1, NAlpha_Str  
            Do R_B_String = 1, NBeta_Str  
              L_Index = 1+(L_B_String-1)*NAlpha_Str+(L_A_String-1) 
              R_Index = 1+(R_B_String-1)*NAlpha_Str+(R_A_String-1) 
              Call CI_Hamiltonian%put(Slater_Condon(IOut,IPrint,NBasis, &
                Determinants,L_A_String,L_B_String,R_A_String,R_B_String, &
                MO_Core_Ham,MO_ERIs,UHF),L_Index,R_Index)
            EndDo
          EndDo
        EndDo
      EndDo
!
      If(IPrint.ge.2) Call MQC_Print(CI_Hamiltonian,IOut,'CI Hamiltonian')

      End Subroutine MQC_Build_CI_Hamiltonian
!
      End Module MQC_EST  
