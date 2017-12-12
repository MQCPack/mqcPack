! A generic linked list object
module MQC_GAUSSIAN_READ_WRITE
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
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
      Use MQC_EST
      use mqc_gaussian
!                                                                    
!----------------------------------------------------------------    
!                                                               |    
!     TYPE AND CLASS DEFINITIONS                                |    
!                                                               |    
!----------------------------------------------------------------    
  implicit none

  ! A public variable to use as a MOLD for transfer()
  integer, dimension(:), allocatable :: list_data


!       we should make Alpha, Beta etc. private once MQC_Gaussian discards old 
!       matrix file reading routine      
  Type MQC_ALL_Wavefunction
     Character(Len=256)::FileName,Basis,Symmetry,WF_Type
     Type(mqc_matrix)::Shell_To_Atom_Map
     Type(mqc_matrix)::Shell_Types
     Type(mqc_matrix)::Number_Of_Primitives_Per_Shell
     Type(mqc_matrix)::Primitive_Exponents
     Type(mqc_matrix)::Contraction_Coefficients
     Type(mqc_matrix)::P__S_P__Contraction_Coefficients
     Type(mqc_matrix)::COORDINATES_OF_EACH_Shell
     Type(MQC_SCF_Integral)::MO_Coefficients
     Type(MQC_SCF_Eigenvalues)::MO_Energies
     Type(MQC_SCF_Eigenvalues)::MO_Symmetries
     Type(MQC_SCF_Integral)::Core_Hamiltonian
     Type(MQC_SCF_Integral)::Fock_Matrix
     Type(MQC_SCF_Integral)::Density_Matrix
     Type(MQC_SCF_Integral)::Overlap_Matrix
     Type(MQC_Scalar)::NAlpha,NBeta,NElectrons,NBasis,Charge,Multiplicity
     Logical::WF_Complex
     Integer::NCore,NVal,NActive
     Type(MQC_Matrix)::PSCF_Amplitudes
     Type(MQC_Vector)::PSCF_Energies
  End Type MQC_ALL_Wavefunction

   ! Linked list node data type
  type mqc_link_list_node
     type(MQC_ALL_Wavefunction) :: current_wavefunction
     type(mqc_link_list_node),pointer :: next =>null()
  end type mqc_link_list_node

  type :: mqc_link_list
     private
     type(mqc_link_list_node),pointer,private :: first => null()
     type(mqc_link_list_node),pointer,private :: last => null()
     type(mqc_link_list_node),pointer,private :: current => null()
     integer ,private :: length = 0
   contains
     procedure :: start => mqc_link_list_start
     procedure :: add => mqc_link_list_add
     procedure :: next => mqc_link_list_next
  end type mqc_link_list

CONTAINS

      subroutine mqc_read_all_wavefunction( current_wavefunction)
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction

      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      integer::iOut=6,iPrint=2,nBasis,nAlpha,nBeta
      type(mqc_scf_integral)::mo_coefficients,overlap
!
!     The header can be loaded from the file using the load procedure.
!
      call fileInfo%load(current_wavefunction%FileName)
!
!     Using the getArray procedure, matrix file arrays can be read in exactly as found in the
!     matrix file.
!
      call fileInfo%getArray('SHELL TO ATOM MAP', &
           current_wavefunction%Shell_To_Atom_Map)
      call fileInfo%getArray('SHELL TYPES', &
           current_wavefunction%Shell_Types)
      call fileInfo%getArray('NUMBER OF PRIMITIVES PER SHELL', &
           current_wavefunction%Number_Of_Primitives_Per_Shell)
      call fileInfo%getArray('PRIMITIVE EXPONENTS', &
           current_wavefunction%Primitive_Exponents)
      call fileInfo%getArray('CONTRACTION COEFFICIENTS', &
           current_wavefunction%Contraction_Coefficients)
      call fileInfo%getArray('P(S=P) CONTRACTION COEFFICIENTS', &
           current_wavefunction%P__S_P__Contraction_Coefficients)
      call fileInfo%getArray('COORDINATES OF EACH SHELL', &
           current_wavefunction%COORDINATES_OF_EACH_Shell)

!
!     The getESTObj procedure allows certain matrix file arrays to be simultaneously read and
!     organized into spin blocks. Arrays are automatically be ordered regardless of whether the
!     matrix file contains RHF, UHF, GHF and/or complex integrals.!
     call fileInfo%getESTObj('overlap',est_integral=overlap)
      call fileInfo%getESTObj('mo coefficients',est_integral=mo_coefficients)


!      allocate(destination_wavefunction%MO_Coefficients, source=current_wavefunction%MO_Coefficients)
!      allocate(destination_wavefunction%MO_Energies, source=current_wavefunction%MO_Energies)
!      allocate(destination_wavefunction%MO_Symmetries, source=current_wavefunction%MO_Symmetries)
!      allocate(destination_wavefunction%Core_Hamiltonian, source=current_wavefunction%Core_Hamiltonian)
!      allocate(destination_wavefunction%Fock_Matrix, source=current_wavefunction%Fock_Matrix)
!      allocate(destination_wavefunction%Density_Matrix, source=current_wavefunction%Density_Matrix)
!      allocate(destination_wavefunction%Overlap_Matrix, source=current_wavefunction%Overlap_Matrix)
!      allocate(destination_wavefunction%NAlpha, source=current_wavefunction%NAlpha)
!      allocate(destination_wavefunction%NBeta, source=current_wavefunction%NBeta)
!      allocate(destination_wavefunction%NElectrons, source=current_wavefunction%NElectrons)
!      allocate(destination_wavefunction%NBasis, source=current_wavefunction%NBasis)
!      allocate(destination_wavefunction%Charge, source=current_wavefunction%Charge)
!      allocate(destination_wavefunction%Multiplicity, source=current_wavefunction%Multiplicity)
!      allocate(destination_wavefunction%Basis, source=current_wavefunction%Basis)
!      allocate(destination_wavefunction%Symmetry, source=current_wavefunction%Symmetry)
!      allocate(destination_wavefunction%WF_Type, source=current_wavefunction%WF_Type)
!      allocate(destination_wavefunction%PSCF_Amplitudes, source=current_wavefunction%PSCF_Amplitudes)
!      allocate(destination_wavefunction%PSCF_Energies, source=current_wavefunction%PSCF_Energies)
!
!      destination_wavefunction%WF_Complex = current_wavefunction%WF_Complex
!      destination_wavefunction%NCore = current_wavefunction%NCore
!      destination_wavefunction%NVal = current_wavefunction%NVal
!      destination_wavefunction%NActive = current_wavefunction%NActive
      call fileInfo%CloseFile()

      end subroutine mqc_read_all_wavefunction
!
!     PROCEDURE MQC_Print_all_Wavefunction     
      subroutine mqc_print_all_wavefunction(wavefunction,iOut,label)
!
      implicit none
      class(MQC_ALL_Wavefunction)::wavefunction
      integer,intent(in)::iOut
      character(len=*),optional,intent(in)::label
      character(len=64)::arrayType,myLabel
!                                                                               
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
      end subroutine mqc_print_all_wavefunction

      subroutine mqc_link_list_add(self,FileName)
        class(mqc_link_list),intent(inout) :: self
        character(len=256),intent(in) :: FileName

        if (.not.associated(self%last)) then
           allocate(self%first)
           self%last => self%first
        else
           allocate(self%last%next)
           self%last => self%last%next
        endif

        self%last%current_wavefunction%FileName=trim(FileName)

        call mqc_read_all_wavefunction( self%last%current_wavefunction)

        self%length = self%length + 1

      end subroutine mqc_link_list_add

      subroutine mqc_link_list_start(self)
        class(mqc_link_list),intent(inout) :: self

        self%current => self%first

      end subroutine mqc_link_list_start

      subroutine mqc_link_list_next(self,wavefunction)
        class(mqc_link_list),intent(inout) :: self
        class(MQC_ALL_Wavefunction),pointer,intent(out) :: wavefunction

        if (associated(self%current)) then
           wavefunction => self%current%current_wavefunction
           self%current => self%current%next
        else
           wavefunction => null()
        end if
      end subroutine mqc_link_list_next

      End Module MQC_GAUSSIAN_READ_WRITE

