! A generic linked list object
module MQC_FullWavefunction
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
      use iso_fortran_env
!                                                                    
!----------------------------------------------------------------    
!                                                               |    
!     TYPE AND CLASS DEFINITIONS                                |    
!                                                               |    
!----------------------------------------------------------------    
  implicit none

  ! A public variable to use as a MOLD for transfer()
  integer(kind=int64), dimension(:), allocatable :: list_data

!       we should make Alpha, Beta etc. private once MQC_Gaussian discards old 
!       matrix file reading routine      
  Type MQC_ALL_Wavefunction
     Character(Len=256)::FileName,Basis,Symmetry,WF_Type
     Character(Len=64)::LabFil
     Character(Len=64)::GVers
     Character(Len=64)::Title
     Integer(kind=int64)::UnitNumber,IVers,NLab,NAtoms,nbasisUse,Len12L,Len4L

     Integer(kind=int64)::IOpCl,ICGU,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB
     Integer(kind=int64)::NPrmDB,NBTot
     Integer(kind=int64)::NAlpha,NBeta,NElectrons,NBasis,iCharge,Multiplicity
     integer,dimension(:),allocatable::atomicNumbers,atomTypes, &
          basisFunction2Atom,IBasisFunctionType
     real(kind=real64),dimension(:),allocatable::atomicCharges,atomicWeights,cartesians

     Type(mqc_vector)::SHELL_TO_ATOM_MAP
     Type(mqc_vector)::SHELL_TYPES
     Type(mqc_vector)::NUMBER_OF_PRIMITIVES_PER_SHELL
     Type(mqc_vector)::PRIMITIVE_EXPONENTS
     Type(mqc_vector)::CONTRACTION_COEFFICIENTS
     Type(mqc_vector)::P__S_P__CONTRACTION_COEFFICIENTS
     Type(mqc_vector)::COORDINATES_OF_EACH_SHELL
     Type(mqc_vector)::INTEGER_ISO
     Type(mqc_vector)::INTEGER_SPIN
     Type(mqc_vector)::REAL_ZEFFECTIVE
     Type(mqc_vector)::REAL_GFACTOR
     Type(mqc_vector)::REAL_ZNUCLEAR
     Type(mqc_vector)::ELECTRIC_DIPOLE_MOMENT
     Type(mqc_vector)::ALPHA_ORBITAL_ENERGIES
     Type(mqc_vector)::BETA_ORBITAL_ENERGIES
     Type(mqc_vector)::BONDS_PER_ATOM
     Type(mqc_vector)::BONDED_ATOMS
     Type(mqc_vector)::BOND_TYPES
     Type(mqc_vector)::GAUSSIAN_SCALARS_I
     Type(mqc_vector)::GAUSSIAN_SCALARS_A

     Type(mqc_matrix)::NUCLEAR_GRADIENT
     Type(mqc_matrix)::OVERLAP
     Type(mqc_matrix)::CORE_HAMILTONIAN_ALPHA
     Type(mqc_matrix)::CORE_HAMILTONIAN_BETA
     Type(mqc_matrix)::KINETIC_ENERGY
     Type(mqc_matrix)::ORTHOGONAL_BASIS
     Type(mqc_matrix)::ALPHA_MO_COEFFICIENTS
     Type(mqc_matrix)::BETA_MO_COEFFICIENTS
     Type(mqc_matrix)::ALPHA_DENSITY_MATRIX
     Type(mqc_matrix)::BETA_DENSITY_MATRIX
     Type(mqc_matrix)::ALPHA_SCF_DENSITY_MATRIX
     Type(mqc_matrix)::BETA_SCF_DENSITY_MATRIX
     Type(mqc_matrix)::ALPHA_FOCK_MATRIX
     Type(mqc_matrix)::BETA_FOCK_MATRIX
     Type(MQC_R4Tensor)::DIPOLE_INTEGRALS
     Type(MQC_R4Tensor)::DIP_VEL_INTEGRALS
     Type(MQC_R4Tensor)::R_X_DEL_INTEGRALS
     Type(MQC_R4Tensor)::RAFFENETTI_2E_INTEGRALS
     Type(MQC_R4Tensor)::AA_MO_2E_INTEGRALS
     Type(MQC_R4Tensor)::AB_MO_2E_INTEGRALS
     Type(MQC_R4Tensor)::BA_MO_2E_INTEGRALS
     Type(MQC_R4Tensor)::BB_MO_2E_INTEGRALS
     Type(MQC_R4Tensor)::REGULAR_2E_INTEGRALS
! There is no logic for the types below
     Type(MQC_SCF_Integral)::MO_Coefficients
     Type(MQC_SCF_Eigenvalues)::MO_Energies
     Type(MQC_SCF_Eigenvalues)::MO_Symmetries
     Type(MQC_SCF_Integral)::Core_Hamiltonian
     Type(MQC_SCF_Integral)::Fock_Matrix
     Type(MQC_SCF_Integral)::Density_Matrix
     Type(MQC_SCF_Integral)::Overlap_Matrix
     Logical::WF_Complex
     Integer(kind=int64)::NCore,NVal,NActive
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
     integer(kind=int64),private :: length = 0
   contains
     procedure :: start => mqc_link_list_start
     procedure :: add => mqc_link_list_add
     procedure :: add_summary => mqc_summary_link_list_add
     procedure :: next => mqc_link_list_next
  end type mqc_link_list

CONTAINS
      subroutine mqc_read_all_wavefunction( current_wavefunction, Print_Options, iout)
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
      integer(kind=int64),intent(in)::Print_Options
      integer(kind=int64),intent(in)::iOut
      type(mqc_scf_integral)::mo_coefficients,overlap
      type(MQC_Matrix)::matrixOut
      type(MQC_Vector)::vectorOut
      Integer(kind=int64):: NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,NRI,LR
      logical ASym
!      Integer(kind=int64):: LStr64, rec_len, bufsize
      Character(Len=64)::CBuf
      logical EOF
      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
      integer(kind=int64),allocatable,dimension(:)::integerTmp
      real(kind=real64),allocatable,dimension(:)::arrayTmp
!
!     Format statements.
!
 1000 format( "Warning: MQC Pack cannot process file '",a,"' it may have been created with the raw attribute.")
 1010 format('CBuf  ',a,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
 1050 format( a,": implement proper handling of this data type")
 1060 format( "The data is identifyed as: ",a)
 1070 format('NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
 2010 Format( i15 )
 2020 Format ( "********************************************************", &
           /,10x,A,/, &
           "********************************************************",/, &
           "Gaussian Version:",A,/,&
           "Job Title:",A)
!
! Open the file and read the first 2 records
!
      call Open_Read(trim(current_wavefunction%FileName), &
           current_wavefunction%UnitNumber, &
           current_wavefunction%LabFil, &
           current_wavefunction%IVers, &
           current_wavefunction%NLab, &
           current_wavefunction%GVers, &
           current_wavefunction%Title, &
           current_wavefunction%NAtoms, &
           current_wavefunction%NBasis, &
           current_wavefunction%nbasisUse, &
           current_wavefunction%iCharge, &
           current_wavefunction%Multiplicity, &
           current_wavefunction%nelectrons, &
           current_wavefunction%Len12L, &
           current_wavefunction%Len4L, &
           current_wavefunction%IOpCl, &
           current_wavefunction%ICGU )

      if ( current_wavefunction%UnitNumber.eq.-2 ) then
         write( 6, 1000 ) trim(current_wavefunction%FileName)
         return
      endif
!
! Read the next 9 records
!
      allocate(current_wavefunction%atomicNumbers(current_wavefunction%natoms),  &
          current_wavefunction%atomTypes(current_wavefunction%natoms),  &
          current_wavefunction%atomicCharges(current_wavefunction%natoms),  &
          current_wavefunction%atomicWeights(current_wavefunction%natoms))

      allocate(current_wavefunction%cartesians(current_wavefunction%natoms*3))
      allocate(current_wavefunction%basisFunction2Atom(current_wavefunction%NBasis),  &
         current_wavefunction%IBasisFunctionType(current_wavefunction%NBasis))
      call Rd_Head(current_wavefunction%unitNumber, &
           current_wavefunction%NLab, &
           current_wavefunction%natoms, &
           current_wavefunction%nbasis, &
           current_wavefunction%atomicNumbers, &
           current_wavefunction%atomTypes, &
           current_wavefunction%atomicCharges, &
           current_wavefunction%cartesians, &
           current_wavefunction%basisFunction2Atom, &
           current_wavefunction%IBasisFunctionType, &
           current_wavefunction%atomicWeights, &
           current_wavefunction%NFC, &
           current_wavefunction%NFV, &
           current_wavefunction%ITran, &
           current_wavefunction%IDum9, &
           current_wavefunction%NShlAO, &
           current_wavefunction%NPrmAO, &
           current_wavefunction%NShlDB, &
           current_wavefunction%NPrmDB, &
           current_wavefunction%NBTot)
!
      if ( Print_Options == 1 ) then
         write( IOut, 2020) trim(current_wavefunction%FileName), &
              trim(current_wavefunction%GVers), &
              trim(current_wavefunction%Title)
         call MQC_PRINT_SCALAR_SUMMARY( iOut, current_wavefunction )
         call MQC_PRINT_DEFAULT_VECTOR_SUMMARY( iOut, current_wavefunction)
      endif

      EOF = .false.
      do while(.not.EOF)
         call Rd_Labl(current_wavefunction%unitNumber, &
              current_wavefunction%IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2, &
              N3,N4,N5,ASym,NRI,EOF)
         LR = LenArr(N1,N2,N3,N4,N5)
#ifdef DEBUG 
         write(IOut,1010) TRIM(cBuf),NI,NR,NRI,NTot,LenBuf,  &
              N1,N2,N3,N4,N5,ASym,LR
#endif
         call String_Change_Case(cBuf,'u')
         if ( Print_Options == 1 ) then
            if('END'.eq.trim(cBuf)) then
               Call Close_MatF(current_wavefunction%unitNumber)
               return
            else
               call MQC_PRINT_DATA_TYPE(Iout,cBuf, &
                 NI,NR,N1,N2,N3,N4,N5,NRI,NTot,LenBuf,LR,Asym)
               Call Rd_Skip(current_wavefunction%UnitNumber,NTot,LenBuf)
            endif
         else
            call String_Change_Case(cBuf,'u')
            select case(trim(CBuf))
            case('END')
               Call Close_MatF(current_wavefunction%unitNumber)
               return
            case('SHELL TO ATOM MAP')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%SHELL_TO_ATOM_MAP)
            case('SHELL TYPES')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%SHELL_TYPES)
            case('NUMBER OF PRIMITIVES PER SHELL')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%NUMBER_OF_PRIMITIVES_PER_SHELL)
            case('PRIMITIVE EXPONENTS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%PRIMITIVE_EXPONENTS)
            case('CONTRACTION COEFFICIENTS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%CONTRACTION_COEFFICIENTS)
            case('P(S=P) CONTRACTION COEFFICIENTS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%P__S_P__CONTRACTION_COEFFICIENTS)
            case('COORDINATES OF EACH SHELL')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%COORDINATES_OF_EACH_SHELL)
            case('GAUSSIAN SCALARS')
               allocate(integerTmp(LR))
               allocate(arrayTmp(LR))
!
! Need to make sure to read in all the bits without possible 32/64 bit
! conversion of the results.  So forcing a 64 bit integer read.
!
               call RdIBuf8(current_wavefunction%unitNumber,NTot,LenBuf,integerTmp)
               call mqc_copy_vec( integerTmp, integerTmp, arrayTmp, LR)
               current_wavefunction%GAUSSIAN_SCALARS_I = Reshape(integerTmp,[LR])
               current_wavefunction%GAUSSIAN_SCALARS_A = Reshape(arrayTmp,[LR])
               deallocate(integerTmp)
               deallocate(arrayTmp)
            case('INTEGER ISO')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%INTEGER_ISO)
            case('INTEGER SPIN')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%INTEGER_SPIN)
            case('REAL ZEFFECTIVE')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%REAL_ZEFFECTIVE)
            case('REAL GFACTOR')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%REAL_GFACTOR)
            case('REAL ZNUCLEAR')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%REAL_ZNUCLEAR)
            case('NUCLEAR GRADIENT')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%NUCLEAR_GRADIENT)
            case('ELECTRIC DIPOLE MOMENT')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%ELECTRIC_DIPOLE_MOMENT)
            case('OVERLAP')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%OVERLAP)
            case('CORE HAMILTONIAN ALPHA')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%CORE_HAMILTONIAN_ALPHA)
            case('CORE HAMILTONIAN BETA')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%CORE_HAMILTONIAN_BETA)
            case('KINETIC ENERGY')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%KINETIC_ENERGY)
            case('ORTHOGONAL BASIS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%ORTHOGONAL_BASIS)
            case('DIPOLE INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%DIPOLE_INTEGRALS)
            case('DIP VEL INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%DIP_VEL_INTEGRALS)
            case('R X DEL INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%R_X_DEL_INTEGRALS)
            case('RAFFENETTI 2E INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%RAFFENETTI_2E_INTEGRALS)
            case('AA MO 2E INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%AA_MO_2E_INTEGRALS)
            case('AB MO 2E INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%AB_MO_2E_INTEGRALS)
            case('BA MO 2E INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%BA_MO_2E_INTEGRALS)
            case('BB MO 2E INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%BB_MO_2E_INTEGRALS)
            case('ALPHA ORBITAL ENERGIES')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%ALPHA_ORBITAL_ENERGIES)
            case('BETA ORBITAL ENERGIES')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%BETA_ORBITAL_ENERGIES)
            case('ALPHA MO COEFFICIENTS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%ALPHA_MO_COEFFICIENTS)
            case('BETA MO COEFFICIENTS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%BETA_MO_COEFFICIENTS)
            case('ALPHA DENSITY MATRIX')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%ALPHA_DENSITY_MATRIX)
            case('BETA DENSITY MATRIX')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%BETA_DENSITY_MATRIX)
            case('ALPHA SCF DENSITY MATRIX')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%ALPHA_SCF_DENSITY_MATRIX)
            case('BETA SCF DENSITY MATRIX')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%BETA_SCF_DENSITY_MATRIX)
            case('ALPHA FOCK MATRIX')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%ALPHA_FOCK_MATRIX)
            case('BETA FOCK MATRIX')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    current_wavefunction%BETA_FOCK_MATRIX)
            case('REGULAR 2E INTEGRALS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut,vectorOut, &
                    current_wavefunction%REGULAR_2E_INTEGRALS)
            case('BONDS PER ATOM')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%BONDS_PER_ATOM)
            case('BONDED ATOMS')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%BONDED_ATOMS)
            case('BOND TYPES')
               call MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
                    NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf, &
                    matrixOut, current_wavefunction%BOND_TYPES)
            case default
               write( 6, 1050 ) trim(CBuf)
               write( 6, 1060 ) MQC_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI)
               write(IOut,1070) NI,NR,NRI,NTot,LenBuf,  &
                    N1,N2,N3,N4,N5,ASym,LR
            Call Rd_Skip(current_wavefunction%UnitNumber,NTot,LenBuf)
            end select
         endif
      enddo
      return
      end subroutine mqc_read_all_wavefunction
!
!     PROCEDURE MQC_Print_all_Wavefunction     
      subroutine mqc_print_all_wavefunction(wavefunction,iOut,label)
!
      implicit none
      class(MQC_ALL_Wavefunction)::wavefunction
      integer(kind=int64),intent(in)::iOut
      character(len=*),optional,intent(in)::label
      character(len=64)::arrayType,myLabel
!                                                                                    write( iout, 1050) "In mqc_print_all_wavefunction"
 1050 Format( 2A )
 1060 Format( A,L10 )
!                                                                              
      if(present(label)) then
        call string_change_case(label,'l',myLabel)
      else
        myLabel = 'all'
      endIf

      select case (myLabel)
      case('mo coefficients')
!        call wavefunction%fock_matrix%print(iOut,'Fock matrix',.true.,.true.)
      case default
!        call mqc_error_A('Invalid label sent to wavefunction print', 6, &
!             'myLabel', myLabel )
      end select
!
      write( iout, 1050) "Out mqc_print_all_wavefunction"
      end subroutine mqc_print_all_wavefunction

      subroutine mqc_summary_link_list_add(self,FileName,iout)
        class(mqc_link_list),intent(inout) :: self
        character(len=256),intent(in) :: FileName
        integer(kind=int64),intent(inout)::iOut
        integer(kind=int64)::Print_Options

        Print_Options = 1;
        if (.not.associated(self%last)) then
           allocate(self%first)
           self%last => self%first
        else
           allocate(self%last%next)
           self%last => self%last%next
        endif

        self%last%current_wavefunction%FileName=trim(FileName)

        call mqc_read_all_wavefunction( self%last%current_wavefunction, &
             Print_Options,iout)

        self%length = self%length + 1
      end subroutine mqc_summary_link_list_add


      subroutine mqc_link_list_add(self,FileName,iout)
        class(mqc_link_list),intent(inout) :: self
        character(len=256),intent(in) :: FileName
        integer(kind=int64),intent(inout)::iOut
        integer(kind=int64)::Print_Options

        Print_Options = 0;
        if (.not.associated(self%last)) then
           allocate(self%first)
           self%last => self%first
        else
           allocate(self%last%next)
           self%last => self%last%next
        endif

        self%last%current_wavefunction%FileName=trim(FileName)

        call mqc_read_all_wavefunction( self%last%current_wavefunction, &
             Print_Options,iout)

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

!=====================================================================
!
!PROCEDURE MQC_Unformatted_Matrix_Read_Array
      subroutine MQC_Unformatted_Matrix_Read_Array(current_wavefunction, &
       NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf,matrixOut,vectorOut,r4TensorOut)
!
!     This Routine is used to look-up a matrix in a unformatted matrix file load
!     that array into either (OPTIONAL) output dummy MQC_Matrix argument
!     <matrixOut>, (OPTIONAL) output dummy MQC_Vector argument <vectorOut>, or
!     (OPTIONAL) output dummy MQC_R4Tensor argument <r4TensorOut>. The character
!     label for the array of interest is sent to this routine in dummy argument
!     <label>.
!
!     Dummy argument <filename> is optional and is only used if fileinfo
!     hasn't already been defined using Routine
!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!     filename sent is different from the filename associated with object
!     fileinfo.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!     called before calling this routine. The expectation is that
!     MQC_Gaussian_Unformatted_Matrix_Read_Header is also called before this
!     routine. However, it is also OK to call this routine first. In that case,
!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     H. P. Hratchian, 2017.
!     L. M. Thompson, 2017
!
!     Variable Declarations.
!
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
      integer(kind=int64),intent(in)::NI,NR,N1,N2,N3,N4,N5,NRI,LR,NTot,LenBuf
      type(MQC_Matrix),intent(inout)::matrixOut
      type(MQC_Vector),intent(inout),OPTIONAL::vectorOut
      type(MQC_R4Tensor),intent(inout),OPTIONAL::r4TensorOut
!
! The last Out structure (matrixOut, vectorOut, r4TensorOut).  Just put in
! unused Out arguments before the desired structure.
!
      integer::iout=6
!
!     Temporary local variables used when calling the gauopen routines.
      character(len=64)::cBuffer,tmpLabel
!
!     Local temp variables.
      integer::i,nOutputArrays
      integer(kind=int64),allocatable,dimension(:)::integerTmp
      real(kind=real64),allocatable,dimension(:)::arrayTmp
      complex(kind=real64),allocatable,dimension(:)::complexTmp
      character(len=256)::my_filename,errorMsg
      logical::DEBUG=.true.,ok,found
      integer::Matrix=0;
      integer::Vector=1;
      integer::r4Tensor=2;
!
!
!     Format statements.
!
 1020 Format( " " )
!
 1040 Format( A, I15 )
 1050 Format( 2A )
 1060 Format( A )
!
!     Begin by seeing if a new file or filename has been sent by the calling
!     program unit. If so, then get the file declared before reading the
!     header information.
!
!
!     Ensure that one and only one output MQC-type array has been sent from the
!     calling program unit.
!
      nOutputArrays = 0
      if(Present(vectorOut)) nOutputArrays = nOutputArrays+1
      if(Present(r4TensorOut)) nOutputArrays = nOutputArrays+1
!
!     Look for the label sent by the calling program unit. If the label is
!     found, then load <matrixOut> with the data on the file.
!

      select case(MQC_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI))

      case('INTEGER-VECTOR')
         allocate(integerTmp(LR))
         call Rd_IBuf(current_wavefunction%unitNumber,NTot,LenBuf,integerTmp)
         if(Vector.eq.nOutputArrays ) then
            vectorOut = integerTmp
         elseif(Matrix.eq.nOutputArrays ) then
            call MQC_Matrix_DiagMatrix_Put(matrixOut,integerTmp)
         else
            call mqc_error_I('Reading vector from Gaussian matrix file, but NO VECTOR SENT to &
                 & procedure.', 6, &
                 'Matrix', Matrix, 'Vector', Vector, &
                 'nOutputArrays', nOutputArrays )
         endIf
         deallocate(integerTmp)
      case('INTEGER-MATRIX')
         if(Matrix.NE.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
              & procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(integerTmp(LR))
         call Rd_IBuf(current_wavefunction%unitNumber,NTot,LenBuf,integerTmp)
         matrixOut = Reshape(integerTmp,[N1,N2])
         deallocate(integerTmp)
      case('INTEGER-SYMMATRIX')
         if(Matrix.NE.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
              & procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(integerTmp(LR))
         call Rd_IBuf(current_wavefunction%unitNumber,NTot,LenBuf,integerTmp)
         call MQC_Matrix_SymmMatrix_Put(matrixOut,integerTmp)
         deallocate(integerTmp)
      case('INTEGER-MULTISYMMATRICES')
         if(Matrix.NE.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
              & procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(integerTmp(LR))
         call Rd_IBuf(current_wavefunction%unitNumber,NTot,LenBuf,integerTmp)
         matrixOut = Reshape(integerTmp,[((N2+1)*N2)/2,N3*N4])
         deallocate(integerTmp)
      case('REAL-VECTOR')
         allocate(arrayTmp(LR))
         call Rd_RBuf(current_wavefunction%unitNumber,NTot,LenBuf,arrayTmp)
         if(Matrix.eq.nOutputArrays ) then
            vectorOut = arrayTmp
         elseif(Vector.eq.nOutputArrays ) then
            call MQC_Matrix_DiagMatrix_Put(matrixOut,arrayTmp)
         else
            call mqc_error_i('Reading vector from Gaussian matrix file, but NO VECTOR SENT to &
                 & procedure.', 6, &
                 'Matrix', Matrix, 'Vector', Vector, &
                 'nOutputArrays', nOutputArrays )
         endIf
         deallocate(arrayTmp)
      case('REAL-MATRIX')
         if(Matrix.ne.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
              & procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(arrayTmp(LR))
         call Rd_RBuf(current_wavefunction%unitNumber,NTot,LenBuf,arrayTmp)
         matrixOut = Reshape(arrayTmp,[N1,N2])
         deallocate(arrayTmp)
      case('REAL-SYMMATRIX')
         if(Matrix.ne.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
              & procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(arrayTmp(LR))
         call Rd_RBuf(current_wavefunction%unitNumber,NTot,LenBuf,arrayTmp)
         call MQC_Matrix_SymmMatrix_Put(matrixOut,arrayTmp)
         deallocate(arrayTmp)
      case('REAL-MULTISYMMATRICES')
         if(Matrix.ne.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
              & procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(arrayTmp(LR))
         call Rd_RBuf(current_wavefunction%unitNumber,NTot,LenBuf,arrayTmp)
         matrixOut = Reshape(arrayTmp,[((N2+1)*N2)/2,N3*N4])
         deallocate(arrayTmp)
      case('COMPLEX-VECTOR')
         allocate(complexTmp(LR))
         call Rd_CBuf(current_wavefunction%unitNumber,NTot,LenBuf,complexTmp)
         if(Vector.eq.nOutputArrays ) then
            vectorOut = complexTmp
         elseif(Matrix.eq.nOutputArrays ) then
            call MQC_Matrix_DiagMatrix_Put(matrixOut,complexTmp)
         else
            call mqc_error_I('Reading vector from Gaussian matrix file, but NO VECTOR SENT &
                 & to procedure.', 6, &
                 'Matrix', Matrix, 'Vector', Vector, &
                 'nOutputArrays', nOutputArrays )
         endIf
         deallocate(complexTmp)
      case('COMPLEX-MATRIX')
         if(Matrix.ne.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix &
              & file, but NO MATRIX SENT to procedure.', 6, &
                 'Matrix', Matrix, 'Vector', Vector, &
                 'nOutputArrays', nOutputArrays )
         allocate(complexTmp(LR))
         call Rd_CBuf(current_wavefunction%unitNumber,NTot,LenBuf,complexTmp)
         matrixOut = Reshape(complexTmp,[N1,N2])
         deallocate(complexTmp)
      case('COMPLEX-SYMMATRIX')
         if(Matrix.ne.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix & 
              & file, but NO MATRIX SENT to procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(complexTmp(LR))
         call Rd_CBuf(current_wavefunction%unitNumber,NTot,LenBuf,complexTmp)
         call MQC_Matrix_SymmMatrix_Put(matrixOut,complexTmp)
         deallocate(complexTmp)
      case('COMPLEX-MULTISYMMATRICES')
         if(Matrix.ne.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix & 
              & file, but NO MATRIX SENT to procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(complexTmp(LR))
         call Rd_CBuf(current_wavefunction%unitNumber,NTot,LenBuf,complexTmp)
         matrixOut = Reshape(complexTmp,[((N2+1)*N2)/2,N3*N4])
         deallocate(complexTmp)
      case('MIXED')
         write(*,1020)
         write(*,1040)' Hrant - LR   = ',LR
         write(*,1040)' Hrant - NR   = ',NR
         write(*,1040)' Hrant - NI   = ',NI
         write(*,1040)' Hrant - NRI  = ',NRI
         write(*,1040)' Hrant - NTot = ',NTot
         write(*,1020)
         call mqc_error_a('No general way to load mixed types as of yet &
              &         We are doing it case-by-case at the moment and this does not match.', 6, &
              'MQC_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI)', &
              MQC_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI) )
      case('2ERIS-SYMSYMR4TENSOR')
         if(r4Tensor.ne.nOutputArrays ) call mqc_error_I('Reading r4 tensor from Gaussian matrix file, but NO R4TENSOR SENT to &
              & procedure.', 6, &
              'r4Tensor', r4Tensor, 'nOutputArrays', nOutputArrays )
         if(NRI.eq.1) then
            allocate(arrayTmp(LR))
            call Rd_2EN(current_wavefunction%unitNumber,NR,LR,NR*LR,NTot,LenBuf,arrayTmp)
            call MQC_Matrix_SymmSymmR4Tensor_Put_Real(r4TensorOut,arrayTmp)
            deallocate(arrayTmp)
         elseIf(NRI.eq.2) then
            allocate(arrayTmp(LR*2))
            !                allocate(complexTmp(LR))
            call Rd_2EN(current_wavefunction%unitNumber,NR,LR,NR*LR,2*NTot,2*LenBuf,arrayTmp)
!                call Rd_2EN(fileinfo%unitNumber,NR,LR,NR*LR,2*NTot,2*LenBuf,complexTmp)
            complexTmp = reshape(arrayTmp, shape(complexTmp))
            call MQC_Matrix_SymmSymmR4Tensor_Put_Complex(r4TensorOut,complexTmp)
            !                deallocate(complexTmp)
            deallocate(arrayTmp)
         endIf
      case('2ERIS-MULTISYMMATRICES')
         if(Matrix.ne.nOutputArrays ) call mqc_error_I('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
              & procedure.', 6, &
              'Matrix', Matrix, 'nOutputArrays', nOutputArrays )
         allocate(arrayTmp(LR))
         call Rd_RBuf(current_wavefunction%unitNumber,NTot,LenBuf,arrayTmp)
         matrixOut = Reshape(arrayTmp,[((N2+1)*N2)/2,N3*N4])
         deallocate(arrayTmp)
      case default
         write(*,1050)' Matrix type: ',Trim(MQC_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI))
         call mqc_error_A('Found strange matrix type in Gaussian matrix read routine.', 6, &
              'MQC_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI)', &
              MQC_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI) )
      end select
! if you want to skip:
!            Call Rd_Skip(current_wavefunction%UnitNumber,NTot,LenBuf)

!
      return
      end subroutine MQC_Unformatted_Matrix_Read_Array

!=====================================================================
!
!PROCEDURE MQC_PRINT_DATA_TYPE
      Subroutine MQC_PRINT_DATA_TYPE(Iout,DataName,NI,NR,N1,N2,N3,N4,N5,NRI,NTot,LenBuf,LR,Asym)
!
!     This function returns a character string indicating the type of array
!     found in a Gaussian matrix file. This is done using NI, NR, N1, N2, N3,
!     N4, N5 and NRIfrom a matrix header in a Gaussian unformatted matrix file to
!     determine the type of array the data corresponds to. The return value will
!     be prepended by "REAL-", "INTEGER-", or "COMPLEX-" and appended by one of
!     the following:
!
!           "VECTOR"           A vector.
!           "MATRIX"           A matrix that is allocated full (M x N).
!           "SYMMATRIX"        A symmetric matrix.
!           "MULTISYMMATRICES" A number of symmetric matrices.
!
!     If the input flags do not uniquely identify a known array type, then this
!     function returns "UNKNOWN".
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      integer(kind=int64),intent(in)::iOut
      character(len=64),intent(in)::DataName
      integer(kind=int64),intent(in)::NI,NR,N1,N2,N3,N4,N5,NRI
      integer(kind=int64),intent(in)::NTot,LenBuf,LR
      logical,intent(in)::Asym
      integer::i
      real(kind=real64)::a
      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
      character(len=8)::tab
!
!     Format statements.
!
 1010 Format( A )
 1020 Format( A,"  - ",A )
 1030 Format( A,"  - 2-D MATRIX allocated full (",a," x ",a,")" )
 1040 Format( A,"  - Total elements: ",a )
 1050 Format( A,"  - SYMMETRIC MATRIX ((",a,"+1)*",a,"/2)" )
 1060 Format( A,"  - VECTOR(",a,")" )
 1070 format('NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
 1080 Format( A,"  - MULTIPLE SYMMETIC MATRICES (((",a,"+1)x",a,"/2) x ",a,")" )
 1090 Format( A,"  - MULTIPLE SYMMETIC MATRICES (((",a,"+1)x",a,"/2) x (",a,"x",a,"))" )
 2010 Format( i15 )
!
!     Do the work...
!
      tab = char( 9 )
      write( iout, 1010) DataName
      if(NR.lt.0.or.NI.lt.0) then
         write( iout, 1020) tab, "UNKNOWN"
         write( OutMath1, 2010) NTot
         call mqc_StripSpaces(OutMath1)
         write( iout, 1040) tab, trim(OutMath1)
         return
      else if(NR.gt.0.and.NI.gt.0) then
        if(NR.eq.1.and.NI.eq.4) then
           write( iout, 1020) tab, "2ERIS"
           write( OutMath1, 2010) NTot
           call mqc_StripSpaces(OutMath1)
           write( iout, 1040) tab, trim(OutMath1)
        else
           write( iout, 1020) tab, "MIXED"
           write( OutMath1, 2010) NTot
           call mqc_StripSpaces(OutMath1)
           write( iout, 1040) tab, trim(OutMath1)
        endIf
      elseIf(NI.eq.1) then
         write( iout, 1020) tab, "INTEGER"
         write( OutMath1, 2010) NTot
         call mqc_StripSpaces(OutMath1)
         write( iout, 1040) tab, trim(OutMath1)
      elseIf(NR.eq.1) then
        if(NRI.eq.1) then
           write( iout, 1020) tab, "REAL"
           write( OutMath1, 2010) NTot
           call mqc_StripSpaces(OutMath1)
           write( iout, 1040) tab, trim(OutMath1)
        elseIf(NRI.eq.2) then
           write( iout, 1020) tab, "COMPLEX"
           write( OutMath1, 2010) NTot
           call mqc_StripSpaces(OutMath1)
           write( iout, 1040) tab, trim(OutMath1)
        endIf
      endIf
      if(N1.gt.1.and.N2.eq.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
         write( OutMath1, 2010) NTot
         call mqc_StripSpaces(OutMath1)
         write( iout, 1060) tab, trim(OutMath1)
      elseIf(N1.gt.1.and.N2.gt.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
         write( OutMath1, 2010) n1
         call mqc_StripSpaces(OutMath1)
         write( OutMath2, 2010) n2
         call mqc_StripSpaces(OutMath2)
         write( iout, 1030) tab, trim(OutMath1), trim(OutMath2)
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
         a = sqrt( 2.0*NTot )
         i = a
         write( OutMath1, 2010) i
         call mqc_StripSpaces(OutMath1)
         write( iout, 1050) tab, trim(OutMath1), trim(OutMath1)
      elseIf(N1.le.-1.and.N2.le.-1.and.N3.le.-1.and.N4.gt.1.and.N5.eq.1) then
         write( iout, 1020) tab, "SYM SYM R4 TENSOR"
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.gt.1.and.N4.eq.1.and.N5.eq.1) then
         if ( (abs(N1)+1)*N2/2*n3 .eq. ntot ) then
            write( OutMath1, 2010) abs(n1)
            call mqc_StripSpaces(OutMath1)
            write( OutMath2, 2010) n2
            call mqc_StripSpaces(OutMath2)
            write( OutMath3, 2010) n3
            call mqc_StripSpaces(OutMath3)
            write( iout, 1080) tab, trim(OutMath1), trim(OutMath2), trim(OutMath3)
         else
            write( iout, 1020) tab, "Bad assumption about size Unknown Data Form"
            write(IOut,1070) NI,NR,NRI,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,LR
         endif
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.gt.1.and.N4.gt.1.and.N5.eq.1) then
         if ( (abs(N1)+1)*N2/2*n3*n4 .eq. ntot ) then
            write( OutMath1, 2010) abs(n1)
            call mqc_StripSpaces(OutMath1)
            write( OutMath2, 2010) n2
            call mqc_StripSpaces(OutMath2)
            write( OutMath3, 2010) n3
            call mqc_StripSpaces(OutMath3)
            write( OutMath4, 2010) n4
            call mqc_StripSpaces(OutMath4)
            write( iout, 1090) tab, trim(OutMath1), trim(OutMath2), trim(OutMath3), trim(OutMath4)
         else
            write( iout, 1020) tab, "Bad assumption about size Unknown Data Form"
            write(IOut,1070) NI,NR,NRI,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,LR
         endif
      else
         write( iout, 1020) tab, "Unknown Data Form"
         write(IOut,1070) NI,NR,NRI,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,LR
      endIf
!
      return
      end subroutine MQC_PRINT_DATA_TYPE
      subroutine mqc_StripSpaces(string)
        character(len=*) :: string
        integer :: stringLen 
        integer :: last, actual

        stringLen = len (string)
        last = 1
        actual = 1

        do while (actual < stringLen)
           if (string(last:last) == ' ') then
              actual = actual + 1
              string(last:last) = string(actual:actual)
              string(actual:actual) = ' '
           else
              last = last + 1
              if (actual < last) &
                   actual = last
           endif
        end do
      end subroutine mqc_StripSpaces

      subroutine MQC_PRINT_SCALAR_SUMMARY( iOut, current_wavefunction)
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
      integer(kind=int64),intent(in)::iOut
      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
      character(len=8)::tab
!
!     Format statements.
!
 2010 Format( i15 )
 2020 Format ( "***************************",/," Scalar Data ",/, &
           "***************************")
 2030 Format( A,/,A,"  - Value: ",A,/,A,"  - ",A,/,A,"  - SCALAR" )
 2040 Format( A,/,A,"  - Value: ",A,/,A,"  - ",A,/,A,"  - SCALAR",/,A,"  - ",A )

      tab = char( 9 )
      write( IOut, 2020)
      write( OutMath1, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Atoms", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%nbasis
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Basis Functions", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%nbasis
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Basis Functions", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%nbasisUse
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Linearly Independent Basis Functions", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%iCharge
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Molecular Charge", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%Multiplicity
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Spin Multiplicity", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%nelectrons
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Electrons", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%IOpCl
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2040) 'Closed/open-shell flag', tab, trim(OutMath1), tab, "INTEGER", tab, tab, &
      'Set if the matrix element file is written after'//NEW_LINE('A')//tab//' &
         & an initial guess or the SCF has completed '//NEW_LINE('A')//tab//' &
         & (otherwise it is -1, meaning unspecified).'

      write( OutMath1, 2010 ) current_wavefunction%ICGU
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2040) "Encoded Indicator of whether the calculation is complex and/or GHF", &
           tab, trim(OutMath1), tab, "INTEGER", tab, tab, &
           'Its three-digit value is interpreted as klm,'//NEW_LINE('A')//tab//' &
           & where k is 1 for the spin-aligned case and'//NEW_LINE('A')//tab//' &
           & 2 for GHF; l is 1 for real and 2 for complex;'//NEW_LINE('A')//tab//' &
           & and m is 1 for RHF/GHF and 2 for UHF (i.e., '//NEW_LINE('A')//tab//' &
           & 1 vs. 2 spin blocks). When k=2, then NBasis is'//NEW_LINE('A')//tab//' &
           & the number of spatial basis functions, but the'//NEW_LINE('A')//tab//' &
           & operator matrices are over the spin orbital'//NEW_LINE('A')//tab//' &
           & basis and hence have dimension k*NBasis.'

      write( OutMath1, 2010 ) current_wavefunction%NFC
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Frozen Core Orbitals", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%NFV
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Frozen Virtual Orbitals", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%NShlAO
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2040) "Number of Contracted Shells of AO Basis Functions", &
           tab, trim(OutMath1), tab, "INTEGER", tab, tab, &
           "Needed if shell data is provided."

      write( OutMath1, 2010 ) current_wavefunction%NPrmAO
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Primitive AO Shells", tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%NShlDB
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2040) "Number of Contracted Shells of Density Fitting Functions", &
           tab, trim(OutMath1), tab, "INTEGER", tab, tab, &
           "Needed if fitting shell data is provided"

      write( OutMath1, 2010 ) current_wavefunction%NPrmDB
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Number of Primitive Density Fitting Shells", &
           tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%NBTot
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2030) "Total Number of Bonds In Connectivity Data, If Any", &
           tab, trim(OutMath1), tab, "INTEGER", tab

      write( OutMath1, 2010 ) current_wavefunction%ITran
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2040) "ITRAN", tab, trim(OutMath1), tab, &
           "INTEGER", tab, tab, &
           'If no MO integrals were stored, ITran=4 if'//NEW_LINE('A')//tab//' &
           & only MOs involving at least one occupied'//NEW_LINE('A')//tab//' &
           & orbital were stored, or ITran=5 if a full'//NEW_LINE('A')//tab//' &
           & transformation was done.'
      return
    end subroutine MQC_PRINT_SCALAR_SUMMARY

      subroutine MQC_PRINT_DEFAULT_VECTOR_SUMMARY( iOut, current_wavefunction)
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
      integer(kind=int64),intent(in)::iOut
      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
      character(len=8)::tab
!
!     Format statements.
!
 2010 Format( i15 )
 2050 Format ( "***************************",/," Non-Scalar Data ",/, &
              "***************************")
 2060 Format( A,/,A,"  - ",A,/,A,"  - Total elements: ",A,/,A,"  - VECTOR(",a,")" )
 2070 Format( A,/,A,"  - ",A,/,A,"  - Total elements: ",A,/,A,"  - VECTOR(",a,")",/,A,"  - ",A )
 2080 Format( A,/,A,"  - ",A,/,A,"  - Total elements: ",A,/,A,"  - 2-D MATRIX allocated full (",a," x ",a,")" )

      tab = char( 9 )
      write( IOut, 2050)
! current_wavefunction%atomicNumbers
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( IOut, 2060) "Atomic Numbers", tab, "INTEGER", tab, trim(OutMath2), tab, trim(OutMath2)

! current_wavefunction%atomTypes
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( IOut, 2060) "Atom Types", tab, "INTEGER", tab, trim(OutMath2), tab, trim(OutMath2)

! current_wavefunction%atomicCharges
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( IOut, 2070) "Atom Charges", tab, "REAL", tab, trim(OutMath2), tab, trim(OutMath2), tab, &
           'Nuclear charges may be different from atomic'//NEW_LINE('A')//tab//' &
           & numbers if ECPs were used.'

! current_wavefunction%cartesians
      write( OutMath1, 2010 ) current_wavefunction%natoms*3
      call mqc_StripSpaces(OutMath1)
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( OutMath3, 2010 ) 3
      call mqc_StripSpaces(OutMath3)
      write( IOut, 2080) "Cartesian Nuclear Coordinates In Bohr", tab, "REAL", &
           tab, trim(OutMath1), tab, trim(OutMath2), trim(OutMath3)

! First half of current_wavefunction%basisFunction2Atom
      write( OutMath1, 2010 ) current_wavefunction%nbasis
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2060) "Map From Basis Functions to Atoms", tab, &
           "INTEGER", tab, trim(OutMath1), tab, trim(OutMath1)

! Second half half of current_wavefunction%basisFunction2Atom
      write( OutMath1, 2010 ) current_wavefunction%nbasis
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2060) "Type Flag For Each Basis Function", tab, &
           "INTEGER", tab, trim(OutMath1), tab, trim(OutMath1)

!           current_wavefunction%atomicWeights
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( IOut, 2060) "Atomic Weights", tab, "REAL", tab, &
           trim(OutMath2), tab, trim(OutMath2)
      return
    end subroutine MQC_PRINT_DEFAULT_VECTOR_SUMMARY
!=====================================================================
!
!PROCEDURE MQC_Unformatted_Matrix_Array_Type
      Function MQC_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI)
!
!     This function returns a character string indicating the type of array
!     found in a Gaussian matrix file. This is done using NI, NR, N1, N2, N3,
!     N4, N5 and NRIfrom a matrix header in a Gaussian unformatted matrix file to
!     determine the type of array the data corresponds to. The return value will
!     be prepended by "REAL-", "INTEGER-", or "COMPLEX-" and appended by one of
!     the following:
!
!           "VECTOR"           A vector.
!           "MATRIX"           A matrix that is allocated full (M x N).
!           "SYMMATRIX"        A symmetric matrix.
!           "MULTISYMMATRICES" A number of symmetric matrices.
!
!     If the input flags do not uniquely identify a known array type, then this
!     function returns "UNKNOWN".
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      integer(kind=int64),intent(in)::NI,NR,N1,N2,N3,N4,N5,NRI
      character(len=64)::MQC_Unformatted_Matrix_Array_Type
!
!
!     Do the work...
!
      MQC_Unformatted_Matrix_Array_Type = "UNKNOWN"
      if(NR.lt.0.or.NI.lt.0) return
      if(NR.gt.0.and.NI.gt.0) then
        if((NR.eq.1.or.NR.eq.2.or.NR.eq.3).and.NI.eq.4) then
          MQC_Unformatted_Matrix_Array_Type = "2ERIS"
        else
           MQC_Unformatted_Matrix_Array_Type = "MIXED"
        endIf
      elseIf(NI.eq.1) then 
        MQC_Unformatted_Matrix_Array_Type = "INTEGER"
      elseIf(NR.eq.1) then
        if(NRI.eq.1) then
          MQC_Unformatted_Matrix_Array_Type = "REAL"
        elseIf(NRI.eq.2) then
          MQC_Unformatted_Matrix_Array_Type = "COMPLEX"
        endIf
      endIf
      if(N1.gt.1.and.N2.eq.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
        MQC_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Unformatted_Matrix_Array_Type)//"-VECTOR"
      elseIf(N1.gt.1.and.N2.gt.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
        MQC_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Unformatted_Matrix_Array_Type)//"-MATRIX"
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
        MQC_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Unformatted_Matrix_Array_Type)//"-SYMMATRIX"
      elseIf(N1.le.-1.and.N2.le.-1.and.N3.le.-1.and.N4.gt.1.and.N5.eq.1) then
        MQC_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Unformatted_Matrix_Array_Type)//"-SYMSYMR4TENSOR"
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.gt.1.and.N4.ge.1.and.N5.eq.1) then
        MQC_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Unformatted_Matrix_Array_Type)//"-MULTISYMMATRICES"
      else
        MQC_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Unformatted_Matrix_Array_Type)//"-UnknownDataForm"
      endIf
!
      return
      end function MQC_Unformatted_Matrix_Array_Type

    End module MQC_FullWavefunction


