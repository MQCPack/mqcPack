! A generic linked list object
module MQC_FullWavefunction
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
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
      Use MQC_EST
      use mqc_gaussian
      use iso_fortran_env, only: int64, real64

!hph+      
!!                                                                    
!!----------------------------------------------------------------    
!!                                                               |    
!!     TYPE AND CLASS DEFINITIONS                                |    
!!                                                               |    
!!----------------------------------------------------------------    
!  implicit none
!
!  ! Linked list array data type
!  type :: mqc_array_link_list_node
!     character(Len=64) :: Data_Name
!     integer(kind=int64) :: NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf
!     Logical Asym
!     character(Len=64):: Data_Type,Shape
!     character(Len=1024):: Formula
!     Real(Kind=real64),contiguous,Pointer::ArrR(:) =>null()
!     Integer(kind=int64),contiguous,Pointer::ArrI(:) =>null()
!     Complex(Kind=real64),contiguous,Pointer::ArrC(:) =>null()
!     type(mqc_array_link_list_node),pointer :: next =>null()
!  end type mqc_array_link_list_node
!
!  type :: mqc_array_link_list
!     type(mqc_array_link_list_node),pointer :: first => null()
!     type(mqc_array_link_list_node),pointer :: last => null()
!     type(mqc_array_link_list_node),pointer :: current => null()
!     integer(kind=int64),private :: length = 0
!   contains
!     procedure :: start => mqc_array_link_list_start
!     procedure :: next => mqc_array_link_list_next
!  end type mqc_array_link_list
!
!! MQC_SCF_Integral substitute
!  Type MQC_Matrix_ptr
!     Integer(kind=int64),Private::NCol=0,NRow=0
!     Character(Len=64),Private::Data_Type ! Real,Integer,Complex
!     Character(Len=64),Private::Storage ! StorFull,StorSymm,StorDiag,StorTriD
!     Real(Kind=real64),contiguous,Pointer::MatR(:) =>null()
!     Integer(kind=int64),contiguous,Pointer::MatI(:) =>null()
!     Complex(Kind=real64),contiguous,Pointer::MatC(:) =>null()
!  end type MQC_Matrix_ptr
!  Type MQC_SCF_Integral_ptr
!     Type(MQC_Matrix_ptr),Private::Alpha,Beta,AlphaBeta,BetaAlpha
!     Character(Len=64),Private::Array_Name ! See below
!     Character(Len=64),Private::Array_Type ! Space,Spin,General
!  end Type MQC_SCF_Integral_ptr
!
!  Type MQC_ALL_Wavefunction
!     Character(Len=256)::FileName,Basis,Symmetry,WF_Type
!     Character(Len=64)::LabFil
!     Character(Len=64)::GVers
!     Character(Len=64)::Title
!     Character(Len=64)::SpinType
!
!     Integer(kind=int64)::UnitNumber,IVers,NLab,NAtoms,nbasisUse,Len12L,Len4L
!     Integer(kind=int64)::IOpCl,ICGU,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB
!     Integer(kind=int64)::NPrmDB,NBTot
!     Integer(kind=int64)::NAlpha,NBeta,NElectrons,NBasis,iCharge,Multiplicity
!     Integer(kind=int64)::NCore,NVal,NActive
!     Logical Empty
!!
!! These Allocatable arrays are used when reading or writing the files.
!!
!     Integer(kind=int64),contiguous,Pointer::atomicNumbers(:) =>null()
!     Integer(kind=int64),contiguous,Pointer::atomTypes(:) =>null()
!     Integer(kind=int64),contiguous,Pointer::basisFunction2Atom(:) =>null()
!     INTEGER(kind=int64),contiguous,Pointer::IBasisFunctionType(:) =>null()
!     Real(kind=real64),contiguous,Pointer::atomicCharges(:) =>null()
!     Real(kind=real64),contiguous,Pointer::atomicWeights(:) =>null()
!     Real(kind=real64),contiguous,Pointer::cartesians(:) =>null()
!
!     Type(mqc_wavefunction)::est_wavefunction
!     Type(mqc_array_link_list)::Array_Data
!  End Type MQC_ALL_Wavefunction
!
!   ! Linked list node data type
!  type mqc_link_list_node
!     type(MQC_ALL_Wavefunction) :: current_wavefunction
!     type(mqc_link_list_node),pointer :: next =>null()
!  end type mqc_link_list_node
!
!  type :: mqc_link_list
!     type(mqc_link_list_node),pointer :: first => null()
!     type(mqc_link_list_node),pointer :: last => null()
!     type(mqc_link_list_node),pointer :: current => null()
!     integer(kind=int64),private :: length = 0
!   contains
!     procedure :: start => mqc_link_list_start
!     procedure :: add => mqc_link_list_add
!     procedure :: next => mqc_link_list_next
!  end type mqc_link_list
!
!  ! Linked list READONLY data type
!  type :: mqc_READONLY_list_node
!     character(Len=64) :: Data_Name
!     type(mqc_READONLY_list_node),pointer :: next =>null()
!  end type mqc_READONLY_list_node
!
!  type :: mqc_READONLY_list
!     type(mqc_READONLY_list_node),pointer :: first => null()
!     type(mqc_READONLY_list_node),pointer :: last => null()
!     type(mqc_READONLY_list_node),pointer :: current => null()
!     integer(kind=int64),private :: length = 0
!   contains
!     procedure :: start => mqc_READONLY_list_start
!     procedure :: next => mqc_READONLY_list_next
!  end type mqc_READONLY_list
!
!CONTAINS
!      subroutine mqc_read_all_wavefunction( current_wavefunction, &
!           my_READONLY_list, Print_Options, Storage_Options, iout)
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
!      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
!      integer(kind=int64),intent(in)::Print_Options
!      integer(kind=int64),intent(in)::Storage_Options
!      integer(kind=int64),intent(in)::iOut
!      Integer(kind=int64):: NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,NRI,LR
!      logical ASym
!      Character(Len=64)::CBuf
!      Character(Len=2048)::DataType,Shape,Formula
!      logical EOF
!      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
!      integer(kind=int32)::i
!      integer(kind=int64)::IERROR, natoms, NBasis
!!
!!     Format statements.
!!
! 1000 format( "Warning: MQC Pack cannot process file '",a,"' it may have been created with the raw attribute.")
! 1005 format( "Matrix Element File:  ",a)
! 1010 format('DataName  ',a,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
!        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
! 1050 format( 'ERROR:',a,': please implement proper handling of this data type')
! 1060 format( "The data is identifyed as: ",a)
! 1070 format('NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
!        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
! 2010 Format( i15 )
! 3000 format( "Finished reading file '",a,"'")
!!
!! Open the file and read the first 2 records
!!
!      call Open_Read(trim(current_wavefunction%FileName), &
!           current_wavefunction%UnitNumber, &
!           current_wavefunction%LabFil, &
!           current_wavefunction%IVers, &
!           current_wavefunction%NLab, &
!           current_wavefunction%GVers, &
!           current_wavefunction%Title, &
!           current_wavefunction%NAtoms, &
!           current_wavefunction%NBasis, &
!           current_wavefunction%nbasisUse, &
!           current_wavefunction%iCharge, &
!           current_wavefunction%Multiplicity, &
!           current_wavefunction%nelectrons, &
!           current_wavefunction%Len12L, &
!           current_wavefunction%Len4L, &
!           current_wavefunction%IOpCl, &
!           current_wavefunction%ICGU )
!
!      If ( current_wavefunction%UnitNumber.eq.-2 ) then
!         write( IOut, 1000 ) trim(current_wavefunction%FileName)
!!
!! Flag that this wavefunction is empty and should be ignored later
!!
!         current_wavefunction%Empty = .True.
!         return
!      EndIF
!
!      If (Storage_Options.eq. 0 ) then
!         write( IOut, 1005 ) trim(current_wavefunction%FileName)
!         current_wavefunction%Empty = .True.
!      Else
!         current_wavefunction%Empty = .False.
!      EndIF
!!
!! Read the next 9 records
!!
!      natoms = current_wavefunction%natoms
!      NBasis = current_wavefunction%NBasis
!      call MQC_Read_MatFile_9_Records( current_wavefunction, &
!           my_READONLY_list, Print_Options, Storage_Options, Natoms, &
!           NBASIS, IOut)
!
!      EOF = .false.
!      IERROR = 0
!      do while(.not.EOF)
!         call Rd_Labl(current_wavefunction%unitNumber, &
!              current_wavefunction%IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2, &
!              N3,N4,N5,ASym,NRI,EOF)
!
!!hph+
!         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,DataType,Shape,Formula)
!!         LR = LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,DataType,Shape,Formula)
!!hph-
!
!#ifdef DEBUG 
!         write(IOut,1010) TRIM(cBuf),NI,NR,NRI,NTot,LenBuf,  &
!              N1,N2,N3,N4,N5,ASym,LR
!#endif
!         call String_Change_Case(cBuf,'u')
!         If('END'.eq.trim(cBuf)) then
!            If (Storage_Options.eq. 0 ) then
!               call MQC_DeAlloc_Int_Ptr( Current_Wavefunction%atomicNumbers, &
!                    'ATOMIC NUMBERS', IOut)
!               call MQC_DeAlloc_Int_Ptr( Current_Wavefunction%atomTypes, &
!                    'ATOM TYPES', IOut)
!               call MQC_DeAlloc_Real_Ptr( Current_Wavefunction%atomicCharges, &
!                    'ATOM CHARGES', IOut)
!               call MQC_DeAlloc_Real_Ptr( Current_Wavefunction%cartesians, &
!                    'CARTESIAN NUCLEAR COORDINATES(BOHR)', IOut)
!               call MQC_DeAlloc_Int_Ptr( Current_Wavefunction%basisFunction2Atom,&
!                    'MAP FROM BASIS FUNCTIONS TO ATOMS', IOut)
!               call MQC_DeAlloc_Int_Ptr( Current_Wavefunction%IBasisFunctionType,&
!                    'TYPE FLAG FOR EACH BASIS FUNCTION', IOut)
!               call MQC_DeAlloc_Real_Ptr( Current_Wavefunction%atomicWeights, &
!                    'ATOM WEIGHTS', IOut)
!            Else
!               Call mqc_ALL_wavefunction_to_est_wavefunction( &
!                    Current_Wavefunction, IOut)
!            EndIf
!            Call Close_MatF(current_wavefunction%unitNumber)
!            write( IOut, 3000 ) trim(current_wavefunction%FileName)
!            return
!            If (IERROR .ne. 0 ) then
!               call mqc_error_I('mqc_read_all_wavefunction: Found errors reading Matrix Element File', IOut)
!            EndIf
!         Else
!            If ( Print_Options == 1 ) Then
!               call MQC_PRINT_DATA_TYPE(Iout,cBuf, &
!                 NI,NR,N1,N2,N3,N4,N5,NRI,NTot,LR, &
!                 DataType,Shape,Formula)
!               flush(IOut)
!            EndIf
!            call MQC_Unformatted_Read_Array(current_wavefunction, &
!                 my_READONLY_list, cBuf, &
!                 NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf,Asym,DataType, &
!                 Shape, Formula,IOut,IERROR,Print_Options,Storage_Options)
!         EndIf
!      EndDo
!      return
!      end subroutine mqc_read_all_wavefunction
!!
!!     PROCEDURE MQC_Read_MatFile_9_Records
!      Subroutine MQC_Read_MatFile_9_Records( wavefunction, &
!           my_READONLY_list, Print_Options, Storage_Options, &
!           Natoms, NBASIS, IOut)
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
!      integer(kind=int64),intent(in)::Print_Options
!      integer(kind=int64),intent(in)::Storage_Options
!      integer(kind=int64),intent(in)::NAtoms,NBasis
!      integer(kind=int64),intent(in)::iOut
!      INTEGER(kind=int64)::NR,NI,NRI,N1,N2,N3,N4,N5,LR,LenBuf
!      Logical Asym
!
!      Character(Len=2048)::Data_Type,Shape,Formula
!      Integer(kind=int32)::I,J,K
!!
!!     Format statements.
!!
! 2020 Format ( "********************************************************", &
!           /,10x,A,/, &
!           "********************************************************",/, &
!           "Gaussian Version:",A,/,&
!           "Job Title:",A)
!!
!! Allocate all the arrays required by rd_head.  Need these not to read.
!! May need them later to write the Matrix Element File
!!
!      call MQC_Alloc_Int_Ptr( wavefunction%atomicNumbers, NAtoms, &
!           'ATOMIC NUMBERS', IOut)
!      call MQC_Alloc_Int_Ptr( wavefunction%atomTypes, NAtoms, &
!           'ATOM TYPES', IOut)
!      call MQC_Alloc_Real_Ptr( wavefunction%atomicCharges, NAtoms, &
!           'ATOM CHARGES', IOut)
!      call MQC_Alloc_Real_Ptr( wavefunction%cartesians, 3*NAtoms, &
!           'CARTESIAN NUCLEAR COORDINATES(BOHR)', IOut)
!      call MQC_Alloc_Int_Ptr( wavefunction%basisFunction2Atom, NBasis, &
!           'MAP FROM BASIS FUNCTIONS TO ATOMS', IOut)
!      call MQC_Alloc_Int_Ptr( wavefunction%IBasisFunctionType, NBasis, &
!           'TYPE FLAG FOR EACH BASIS FUNCTION', IOut)
!      call MQC_Alloc_Real_Ptr( wavefunction%atomicWeights, NAtoms, &
!           'ATOM WEIGHTS', IOut)
!
!      call Rd_Head(wavefunction%unitNumber, wavefunction%NLab, &
!           natoms, nbasis, wavefunction%atomicNumbers, wavefunction%atomTypes, &
!           wavefunction%atomicCharges, wavefunction%cartesians, &
!           wavefunction%basisFunction2Atom, wavefunction%IBasisFunctionType, &
!           wavefunction%atomicWeights, wavefunction%NFC, wavefunction%NFV, &
!           wavefunction%ITran, wavefunction%IDum9, wavefunction%NShlAO, &
!           wavefunction%NPrmAO, wavefunction%NShlDB, wavefunction%NPrmDB, &
!           wavefunction%NBTot)
!!
!      if ( Print_Options .EQ. 1 ) then
!!        Write a Summary of the data.  No need to save it
!         write( IOut, 2020) trim(wavefunction%FileName), &
!              trim(wavefunction%GVers), trim(wavefunction%Title)
!         call MQC_PRINT_SCALAR_SUMMARY( iOut, wavefunction )
!         call MQC_PRINT_DEFAULT_VECTOR_SUMMARY( iOut, wavefunction)
!      Else If ( Print_Options == 2 ) then
!         call MQC_PRINT_SCALAR_HEAD( iOut, wavefunction )
!      EndIf
!      if ( Storage_Options .NE. 0 .OR. Print_Options .EQ. 2 ) then
!         NRI=1
!         N1=NAtoms
!         N2=1
!         N3=1
!         N4=1
!         N5=1
!         NR=0
!         NI=1
!
!!hph+
!         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!         LR = LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!hph-
!
!! LenBuf and ASM should never be used.
!         LenBuf = 5000
!         Asym = .False.
!         call mqc_array_link_list_add(wavefunction, "ATOMIC NUMBERS", &
!              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut,Data_Type,Shape,Formula)
!         wavefunction%Array_Data%last%ArrI(1:NAtoms) => &
!              wavefunction%atomicNumbers
!         If ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!         EndIf
!
!         call mqc_array_link_list_add(wavefunction, "ATOM TYPES", &
!              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut,Data_Type,Shape,Formula)
!         wavefunction%Array_Data%last%ArrI(1:NAtoms) => &
!              wavefunction%atomTypes
!         If ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!         EndIf
!
!         NR=1
!         NI=0
!
!!hph+
!         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!         LR = LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!hph-
!
!         call mqc_array_link_list_add(wavefunction,"ATOM CHARGES", &
!              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut,Data_Type,Shape,Formula)
!         wavefunction%Array_Data%last%ArrR(1:NAtoms) => &
!              wavefunction%atomicCharges
!         If ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!         EndIf
!
!         N1=3
!         N2=NAtoms
!
!!hph+
!         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!         LR = LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!hph-
!
!         call mqc_array_link_list_add(wavefunction, &
!              'CARTESIAN NUCLEAR COORDINATES(BOHR)', &
!              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut,Data_Type,Shape,Formula)
!         wavefunction%Array_Data%last%ArrR(1:3*NAtoms) => &
!              wavefunction%cartesians
!         If ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!         EndIf
!
!         N1=NBasis
!         N2=1
!         NR=0
!         NI=1
!
!!hph+
!         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!         LR = LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!hph-
!
!         call mqc_array_link_list_add(wavefunction, &
!              "MAP FROM BASIS FUNCTIONS TO ATOMS", &
!              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut,Data_Type,Shape,Formula)
!         wavefunction%Array_Data%last%ArrI(1:NBasis) => &
!              wavefunction%basisFunction2Atom
!         If ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!         EndIf
!
!         call mqc_array_link_list_add(wavefunction, &
!              "TYPE FLAG FOR EACH BASIS FUNCTION", &
!              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut,Data_Type,Shape,Formula)
!         wavefunction%Array_Data%last%ArrI(1:NBasis) => &
!              wavefunction%IBasisFunctionType
!         If ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!         EndIf
!
!         N1=NAtoms
!         NR=1
!         NI=0
!
!!hph+
!         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!         LR = LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
!!hph-
!
!         call mqc_array_link_list_add(wavefunction, "ATOM WEIGHTS", &
!              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut,Data_Type,Shape,Formula)
!         wavefunction%Array_Data%last%ArrR(1:NAtoms) => &
!              wavefunction%atomicWeights
!         If ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!         EndIf
!      endif
!
!      return
!      end Subroutine MQC_Read_MatFile_9_Records
!!
!!     PROCEDURE MQC_Print_all_Wavefunction     
!      subroutine mqc_print_all_wavefunction(wavefunction,iOut)
!!
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      integer(kind=int64),intent(in)::iOut
!      integer(kind=int64)::LR
!
! 1050 Format( 2A )
! 1060 Format( A,L10 )
!
!      call MQC_PRINT_SCALAR_HEAD( iOut, wavefunction )
!!                                                                              
!      wavefunction%Array_Data%last => wavefunction%Array_Data%first
!      Do While( associated(wavefunction%Array_Data%last) )
!         call mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!         wavefunction%Array_Data%last => wavefunction%Array_Data%last%next
!      enddo
!!
!      write( iout, 1050) "Finished mqc_print_all_wavefunction"
!      end subroutine mqc_print_all_wavefunction
!
!!
!!     PROCEDURE mqc_print_Last_Array_In_wavefunction
!      subroutine mqc_print_Last_Array_In_wavefunction(wavefunction,iOut)
!!
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      integer(kind=int64),intent(in)::iOut
!      integer(kind=int64)::LR
!
! 1050 Format( 2A )
! 1060 Format( A,L10 )
!
!      if (trim(wavefunction%Array_Data%last%Data_Type).eq."INTEGER" ) then
!         call MQC_Print_Array_Integer(IOut, wavefunction%Array_Data%last)
!      else if (trim(wavefunction%Array_Data%last%Data_Type).eq."REAL" ) then
!         call MQC_Print_Array_REAL(IOut, wavefunction%Array_Data%last)
!      else if (trim(wavefunction%Array_Data%last%Data_Type).eq."COMPLEX" ) then
!         call MQC_Print_Array_COMPLEX(IOut, wavefunction%Array_Data%last)
!      else if (trim(wavefunction%Array_Data%last%Data_Type).eq."MIXED" ) then
!! For some reason, LR and N1 are improperly much larger than NTot.  
!! For write purposes, set LR to NTot
!
!         LR = wavefunction%Array_Data%last%LR
!         wavefunction%Array_Data%last%LR = min( LR, &
!              wavefunction%Array_Data%last%NTot)
!         if (wavefunction%Array_Data%last%NRI.EQ.2 ) then
!            call MQC_Print_Array_CMIXED(IOut, &
!                 wavefunction%Array_Data%last)
!         else
!            call MQC_Print_Array_RMIXED(IOut, &
!                 wavefunction%Array_Data%last)
!         endif
!         wavefunction%Array_Data%last%LR=LR
!      else if (trim(wavefunction%Array_Data%last%Data_Type).eq."TENSOR" ) then
!         if (wavefunction%Array_Data%last%NR.lt.0 ) then
!            call MQC_Print_Array_COMPLEX(IOut, &
!                 wavefunction%Array_Data%last)
!         else
!            call MQC_Print_Array_REAL(IOut, &
!                 wavefunction%Array_Data%last)
!         endif
!      else
!         call MQC_PRINT_DATA_TYPE(IOut, &
!              wavefunction%Array_Data%last%Data_Name, &
!              wavefunction%Array_Data%last%NI, &
!              wavefunction%Array_Data%last%NR, &
!              wavefunction%Array_Data%last%N1, &
!              wavefunction%Array_Data%last%N2, &
!              wavefunction%Array_Data%last%N3, &
!              wavefunction%Array_Data%last%N4, &
!              wavefunction%Array_Data%last%N5, &
!              wavefunction%Array_Data%last%NRI, &
!              wavefunction%Array_Data%last%NTot, &
!              wavefunction%Array_Data%last%LR, &
!              wavefunction%Array_Data%last%Data_Type, &
!              wavefunction%Array_Data%last%Shape, &
!              wavefunction%Array_Data%last%Formula)
!      endif
!    End subroutine mqc_print_Last_Array_In_wavefunction
!
!!
!!     PROCEDURE MQC_WRITE_all_Wavefunction     
!      subroutine mqc_WRITE_all_wavefunction(wavefunction,iOut)
!!
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      integer(kind=int64),intent(in)::iOut
!      logical I_OPENED
!
! 1030 format( "File was Open. Close then open to write")
! 1050 Format( 2A )
! 1060 Format( A,L10 )
!           
!      INQUIRE (wavefunction%UnitNumber, OPENED=I_OPENED )
!      if ( I_OPENED ) then
!         write(IOut, 1030)
!         call Close_MatF(wavefunction%UnitNumber)
!      endif
!
!      call Open_Write(trim(wavefunction%FileName), &
!           wavefunction%UnitNumber, &
!           wavefunction%LabFil, &
!           wavefunction%GVers, &
!           wavefunction%Title, &
!           wavefunction%NAtoms, &
!           wavefunction%NBasis, &
!           wavefunction%nbasisUse, &
!           wavefunction%iCharge, &
!           wavefunction%Multiplicity, &
!           wavefunction%nelectrons, &
!           wavefunction%IOpCl, &
!           wavefunction%ICGU )
!      call Wr_Head(wavefunction%unitNumber, &
!           wavefunction%natoms,wavefunction%natoms*3, &
!           wavefunction%nbasis, wavefunction%atomicNumbers, &
!           wavefunction%atomTypes, &
!           wavefunction%atomicCharges, wavefunction%cartesians, &
!           wavefunction%basisFunction2Atom, wavefunction%IBasisFunctionType, &
!           wavefunction%atomicWeights, wavefunction%NFC, wavefunction%NFV, &
!           wavefunction%ITran, wavefunction%IDum9, wavefunction%NShlAO, &
!           wavefunction%NPrmAO, wavefunction%NShlDB, wavefunction%NPrmDB, &
!           wavefunction%NBTot)
!!                                                                              
!      wavefunction%Array_Data%last => wavefunction%Array_Data%first
!      Do While( associated(wavefunction%Array_Data%last) )
!! Do not rewrite arrays that have just been written in the first 9 records.
!         If ( wavefunction%Array_Data%last%Data_Name .NE. "ATOMIC NUMBERS" .AND. &
!            wavefunction%Array_Data%last%Data_Name .NE. "ATOM TYPES" .AND. &
!            wavefunction%Array_Data%last%Data_Name .NE. "ATOM CHARGES" .AND. &
!            wavefunction%Array_Data%last%Data_Name &
!            .NE. "CARTESIAN NUCLEAR COORDINATES(BOHR)" .AND. &
!            wavefunction%Array_Data%last%Data_Name &
!            .NE. "MAP FROM BASIS FUNCTIONS TO ATOMS" .AND. &
!            wavefunction%Array_Data%last%Data_Name &
!            .NE. "TYPE FLAG FOR EACH BASIS FUNCTION" .AND. &
!            wavefunction%Array_Data%last%Data_Name .NE. "ATOM WEIGHTS" ) then
!            call mqc_WRITE_Last_Array_In_wavefunction(wavefunction, IOut)
!         EndIf
!         wavefunction%Array_Data%last => wavefunction%Array_Data%last%next
!      enddo
!!
!      write( iout, 1050) "Finished mqc_print_all_wavefunction"
!      end subroutine mqc_WRITE_all_wavefunction
!
!!
!!     PROCEDURE mqc_WRITE_Last_Array_In_wavefunction
!      subroutine mqc_WRITE_Last_Array_In_wavefunction(wavefunction,iOut)
!!
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      integer(kind=int64),intent(in)::iOut
!
! 1050 Format("WARNING: Do not know how to write this data:")
!
!      if (trim(wavefunction%Array_Data%last%Data_Type).eq."INTEGER" ) then
!         call Wr_LIBuf(wavefunction%unitNumber, &
!              wavefunction%Array_Data%last%Data_Name, &
!              wavefunction%Array_Data%last%NI, &
!              wavefunction%Array_Data%last%LenBuf, &
!              wavefunction%Array_Data%last%N1, &
!              wavefunction%Array_Data%last%N2, wavefunction%Array_Data%last%N3, &
!              wavefunction%Array_Data%last%N4, &
!              wavefunction%Array_Data%last%N5, wavefunction%Array_Data%last%ASym, &
!              wavefunction%Array_Data%last%ArrI)
!      else if (trim(wavefunction%Array_Data%last%Data_Type).eq."REAL" ) then
!         call Wr_LRBuf(wavefunction%unitNumber, &
!              wavefunction%Array_Data%last%Data_Name, &
!              wavefunction%Array_Data%last%NR, &
!              wavefunction%Array_Data%last%LenBuf, &
!              wavefunction%Array_Data%last%N1, &
!              wavefunction%Array_Data%last%N2, wavefunction%Array_Data%last%N3, &
!              wavefunction%Array_Data%last%N4, &
!              wavefunction%Array_Data%last%N5, wavefunction%Array_Data%last%ASym, &
!              wavefunction%Array_Data%last%ArrR)
!      Else If (trim(wavefunction%Array_Data%last%Data_Type).eq."COMPLEX" ) then
!! Wr_LCBuf has a bug, so write the label and the the data seperately
!         call Wr_Labl(wavefunction%unitNumber, &
!              wavefunction%Array_Data%last%Data_Name, &
!              wavefunction%Array_Data%last%NI,&
!              -wavefunction%Array_Data%last%NR, &
!              wavefunction%Array_Data%last%NTot, &
!              wavefunction%Array_Data%last%LenBuf, &
!              wavefunction%Array_Data%last%N1, &
!              wavefunction%Array_Data%last%N2, &
!              wavefunction%Array_Data%last%N3, &
!              wavefunction%Array_Data%last%N4, &
!              wavefunction%Array_Data%last%N5, &
!              wavefunction%Array_Data%last%ASym)
!         call Wr_CBuf(wavefunction%unitNumber, &
!              wavefunction%Array_Data%last%NTOT* &
!              wavefunction%Array_Data%last%NR, &
!              wavefunction%Array_Data%last%LenBuf, &
!              wavefunction%Array_Data%last%ArrC)
!      Else If (trim(wavefunction%Array_Data%last%Data_Type).eq."TENSOR" ) then
!         call Wr_LAO2E(wavefunction%unitNumber, &
!              wavefunction%Array_Data%last%Data_Name, &
!              wavefunction%Array_Data%last%NR, &
!              wavefunction%Array_Data%last%LenBuf, &
!              abs(wavefunction%Array_Data%last%N1), &
!              wavefunction%Array_Data%last%ArrR)
!      Else  If (trim(wavefunction%Array_Data%last%Data_Type).eq."MIXED" ) then
!         If (wavefunction%Array_Data%last%NRI.eq.2 ) then
!            call Wr_Labl(wavefunction%unitNumber, &
!                 wavefunction%Array_Data%last%Data_Name, &
!                 wavefunction%Array_Data%last%NI,&
!                 -wavefunction%Array_Data%last%NR, &
!                 wavefunction%Array_Data%last%NTot, &
!                 wavefunction%Array_Data%last%LenBuf, &
!                 wavefunction%Array_Data%last%N1, &
!                 wavefunction%Array_Data%last%N2, &
!                 wavefunction%Array_Data%last%N3, &
!                 wavefunction%Array_Data%last%N4, &
!                 wavefunction%Array_Data%last%N5, &
!                 wavefunction%Array_Data%last%ASym)
!            call Wr_CMixed(wavefunction%unitNumber, &
!                 wavefunction%Array_Data%last%NTot,&
!                 wavefunction%Array_Data%last%LenBuf, &
!                 wavefunction%Array_Data%last%ArrI, &
!                 wavefunction%Array_Data%last%ArrC)
!         Else
!            call Wr_Labl(wavefunction%unitNumber, &
!                 wavefunction%Array_Data%last%Data_Name, &
!                 wavefunction%Array_Data%last%NI,&
!                 wavefunction%Array_Data%last%NR, &
!                 wavefunction%Array_Data%last%NTot, &
!                 wavefunction%Array_Data%last%LenBuf, &
!                 wavefunction%Array_Data%last%N1, &
!                 wavefunction%Array_Data%last%N2, &
!                 wavefunction%Array_Data%last%N3, &
!                 wavefunction%Array_Data%last%N4, &
!                 wavefunction%Array_Data%last%N5, &
!                 wavefunction%Array_Data%last%ASym)
!            Call Wr_RMixed(wavefunction%unitNumber, &
!                 wavefunction%Array_Data%last%NTot,&
!                 wavefunction%Array_Data%last%LenBuf, &
!                 wavefunction%Array_Data%last%ArrI, &
!                 wavefunction%Array_Data%last%ArrR)
!         Endif
!      Else
!         write( IOut, 1050)
!         call MQC_PRINT_DATA_TYPE(IOut, &
!              wavefunction%Array_Data%last%Data_Name, &
!              wavefunction%Array_Data%last%NI, &
!              wavefunction%Array_Data%last%NR, &
!              wavefunction%Array_Data%last%N1, &
!              wavefunction%Array_Data%last%N2, &
!              wavefunction%Array_Data%last%N3, &
!              wavefunction%Array_Data%last%N4, &
!              wavefunction%Array_Data%last%N5, &
!              wavefunction%Array_Data%last%NRI, &
!              wavefunction%Array_Data%last%NTot, &
!              wavefunction%Array_Data%last%LR, &
!              wavefunction%Array_Data%last%Data_Type, &
!              wavefunction%Array_Data%last%Shape, &
!              wavefunction%Array_Data%last%Formula)
!      EndIf
!    End subroutine mqc_WRITE_Last_Array_In_wavefunction
!
!!
!!     PROCEDURE MQC_Find_mqc_array_link_list_node
!      Integer(kind=int64) Function MQC_Find_mqc_array_link_list_node(wavefunction,Name,Found, &
!           iOut)
!!
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      character(len=*),intent(in) :: Name
!      integer(kind=int64),intent(in)::iOut
!      type(mqc_array_link_list_node),pointer,intent(out) ::Found
!      integer(kind=int64)::LR
! 1050 Format( 2A )
! 1060 Format( A,L10 )
!
!      wavefunction%Array_Data%last => wavefunction%Array_Data%first
!      Do While( associated(wavefunction%Array_Data%last) )
!         if (trim(wavefunction%Array_Data%last%Data_Name) .eq. trim(Name) ) then
!            Found => wavefunction%Array_Data%last
!            If ( associated(Found%ArrI) ) then
!               If ( associated(Found%ArrR) ) then
!                  MQC_Find_mqc_array_link_list_node = 4
!               Else If ( associated(Found%ArrC) ) then
!                  MQC_Find_mqc_array_link_list_node = 5
!               Else 
!                  MQC_Find_mqc_array_link_list_node = 1
!               EndIf
!            Else If ( associated(Found%ArrR) ) then
!               MQC_Find_mqc_array_link_list_node = 2
!            Else If ( associated(Found%ArrC) ) then
!               MQC_Find_mqc_array_link_list_node = 3
!            EndIf
!            return
!         endif
!         wavefunction%Array_Data%last => wavefunction%Array_Data%last%next
!      enddo
!!
!      Found => null()
!      MQC_Find_mqc_array_link_list_node = 0
!
!      end Function MQC_Find_mqc_array_link_list_node
!
!      Subroutine MQC_Print_Array_Integer(IOut, Array,Blank_At_Top, &
!        Blank_At_Bottom)
!      implicit none
!      integer(kind=int64),intent(in)::iOut
!      class(mqc_array_link_list_node),intent(in):: Array
!      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
!      Integer(kind=int64),Parameter::ColWidth=10 
!      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
!      Integer(kind=int64)::N1,N2,N3,N4,K,L,IBase,JBase,LMAX
!
!      Integer(Kind=int64),contiguous,Pointer::I2D(:,:) =>null()
!      Integer(Kind=int64),contiguous,Pointer::I3D(:,:,:) =>null()
!!
!!     Format statements.
!!
! 1000 Format(/,1x,A)
! 1001 Format(1x,I7,2x,I14)
! 1010 Format(1x,A,2x,I14)
! 1020 Format(" ")
! 2001 Format(/,5x,10(7x,I7))
! 2002 Format(/,5x,10(5x,I7,",",I1))
! 2003 Format(/,5x,10(4x,I7,",",I2))
! 2004 Format(/,5x,10(I7,",",I6))
! 2005 Format(/,5x,10(3x,I7,",",I1,",",I1))
! 2006 Format(/,5x,10(3x,I7,",",I7,",",I7))
! 2011 Format(1x,I7,10I14)
! 4010 format( "(",I5,",",I5,"|",I5,",",I5,") = ",I10 )
! 4020 format( "(",I5,",",I5,",",I5,") = ",I10 )
! 5010 Format("Warning: Do Not know how do deal with the shape of ",A)
!
!      If(PRESENT(Blank_At_Top)) then
!        If(Blank_At_Top) Write(IOut,1020)
!      EndIf
!      if ( Array%Shape.ne."SCALAR") then
!         Write(IOut,1000) TRIM(Array%Data_Name)
!      EndIf
!
!      if ( Array%Shape.eq."SCALAR") then
!         Write(IOut,1010) TRIM(Array%Data_Name), Array%ArrI(1)
!      else if ( Array%Shape.eq."VECTOR") then
!         Do I = 1, Array%LR
!            Write(IOut,1001) I,  Array%ArrI(I)
!         EndDo
!      else if ( Array%Shape.eq."MATRIX") then
!         NRows = Array%N1
!         NCols = Array%N2
!         I2D(1:NRows,1:NCols)=>Array%ArrI
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = 1,NRows
!               Write(IOut,2011) I, (I2D(I,J),J=IFirst,ILast)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMMATRIX") then
!         NRows = Array%N2
!         NCols = Array%N2
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = IFirst, ILAST
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Array%ArrI(IBase+J),J=IFirst,I)
!            EndDo
!            Do I = ILAST+1, NCols
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Array%ArrI(IBase+J),J=IFirst,ILast)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."MULTISYMMATRICES") then
!         NRows = Array%N2
!         NCols = Array%N2
!         N3 = Array%N3
!         N4 = Array%N4
!         I3D(1:(NCols+1)*NCols/2,1:N3,1:N4)=>Array%ArrI
!         Do L = 1,N4
!            Do K = 1,N3
!               IBase = 0
!               Do IFirst = 1,NCols,ColWidth
!                  ILast = Min(IFirst+ColWidth-1,NCols)
!                  if ( N4 .eq. 1 ) then
!                     if ( K .lt. 10 ) then
!                        Write(IOut,2002) (I,K,I=IFirst,ILast)
!                     else if ( K .lt. 100 ) then
!                        Write(IOut,2003) (I,K,I=IFirst,ILast)
!                     else
!                        Write(IOut,2004) (I,K,I=IFirst,ILast)
!                     endif
!                  else 
!                     If ( K .lt. 10 .and. L .lt. 10 ) then
!                        Write(IOut,2005) (I,K,L,I=IFirst,ILast)
!                     Else
!                        Write(IOut,2006) (I,K,L,I=IFirst,ILast)
!                     EndIf
!                  endif
!                  Do I = IFirst, ILAST
!                     JBase = IBase+(I*(I-1))/2
!                     Write(IOut,2011) I, (I3D(JBase+J,K,L),J=IFirst,I)
!                  EndDo
!                  Do I = ILAST+1, NCols
!                     JBase = IBase+(I*(I-1))/2
!                     Write(IOut,2011) I, (I3D(JBase+J,K,L),J=IFirst,ILast)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMSYM") then
!         N4 = Array%N4
!         IBase=0
!         Do I = 1,N4
!            Do J = 1,I
!               Do K = 1,I
!                  if ( I .eq. K ) then
!                     LMax = J
!                  else
!                     LMax = K
!                  endif
!                  Do L = 1,LMax
!                     IBase=IBase+1
!                     write(Iout,4010) I,J,K,L,Array%ArrI(IBase)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."TRISYMMATRICES") then
!         N3 = Array%N3
!         IBase=0
!         Do I = 1,N3
!            Do J = 1,I
!               Do K = 1,J
!                  IBase=IBase+1
!                  write(Iout,4020) I,J,K,Array%ArrI(IBase)
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."3DMATRIX") then
!         N1 = Array%N1
!         N2 = Array%N2
!         N3 = Array%N3
!         I3D(1:N1,1:N2,1:N3)=>Array%ArrI
!         IBase=0
!         Do K = 1,N3
!            Do J = 1,N2
!               Do I = 1,N1
!                  IBase=IBase+1
!                  write(Iout,4020) I,J,K,I3D(I,J,K)
!               EndDo
!            EndDo
!         EndDo
!      else
!         write(IOut,5010) trim(Array%Data_Name)
!         call MQC_PRINT_DATA_TYPE(IOut, Array%Data_Name, Array%NI, Array%NR, &
!              Array%N1, Array%N2, Array%N3, Array%N4, Array%N5, Array%NRI, &
!              Array%NTot, Array%LR, Array%Data_Type, Array%Shape, &
!              Array%Formula)
!      endif
!
!      If(PRESENT(Blank_At_Bottom)) then
!        If(Blank_At_Bottom) Write(IOut,1020)
!      EndIf
!      return
!      end Subroutine MQC_Print_Array_Integer
!
!      Subroutine MQC_Print_Array_REAL(IOut, Array,Blank_At_Top, &
!        Blank_At_Bottom)
!      implicit none
!      class(mqc_array_link_list_node),intent(in):: Array
!      integer(kind=int64),intent(in)::iOut
!      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
!      Integer(kind=int64),Parameter::ColWidth=10 
!      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
!      Integer(kind=int64)::N1,N2,N3,N4,K,L,IBase,JBase,LMAX
!      Real(Kind=real64),contiguous,Pointer::R2D(:,:) =>null()
!      Real(Kind=real64),Pointer::R3D(:,:,:) =>null()
!!
!!     Format statements.
!!
! 1000 Format(/,1x,A)
! 1001 Format(1x,I7,2x,F14.6)
! 1010 Format(1x,A,2x,F14.6)
! 1020 Format(" ")
! 2001 Format(/,5x,10(7x,I7))
! 2002 Format(/,5x,10(5x,I7,",",I1))
! 2003 Format(/,5x,10(4x,I7,",",I2))
! 2004 Format(/,5x,10(I7,",",I6))
! 2005 Format(/,5x,10(3x,I7,",",I1,",",I1))
! 2006 Format(/,5x,10(3x,I7,",",I7,",",I7))
! 2011 Format(1x,I7,10F14.6)
! 4010 format( "(",I5,",",I5,"|",I5,",",I5,") = ",F14.6 )
! 4020 format( "(",I5,",",I5,",",I5,") = ",F14.6 )
! 5010 Format("Warning: Do Not know how do deal with the shape of ",A)
!
!      If(PRESENT(Blank_At_Top)) then
!        If(Blank_At_Top) Write(IOut,1020)
!      EndIf
!      if ( Array%Shape.ne."SCALAR") then
!         Write(IOut,1000) TRIM(Array%Data_Name)
!      EndIf
!
!      if ( Array%Shape.eq."SCALAR") then
!         Write(IOut,1010) TRIM(Array%Data_Name), Array%ArrR(1)
!      else if ( Array%Shape.eq."VECTOR") then
!         Do I = 1, Array%LR
!            Write(IOut,1001) I,  Array%ArrR(I)
!         EndDo
!      else if ( Array%Shape.eq."MATRIX") then
!         NRows = Array%N1
!         NCols = Array%N2
!         R2D(1:NRows,1:NCols)=>Array%ArrR
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = 1,NRows
!               Write(IOut,2011) I, (R2D(I,J),J=IFirst,ILast)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMMATRIX") then
!         NRows = Array%N2
!         NCols = Array%N2
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = IFirst, ILAST
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Array%ArrR(IBase+J),J=IFirst,I)
!            EndDo
!            Do I = ILAST+1, NCols
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Array%ArrR(IBase+J),J=IFirst,ILast)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."MULTISYMMATRICES") then
!         NRows = Array%N2
!         NCols = Array%N2
!         N3 = Array%N3
!         N4 = Array%N4
!         R3D(1:(NCols+1)*NCols/2,1:N3,1:N4)=>Array%ArrR
!
!         Do L = 1,N4
!            Do K = 1,N3
!               Do IFirst = 1,NCols,ColWidth
!                  ILast = Min(IFirst+ColWidth-1,NCols)
!                  if ( N4 .eq. 1 ) then
!                     if ( K .lt. 10 ) then
!                        Write(IOut,2002) (I,K,I=IFirst,ILast)
!                     else if ( K .lt. 100 ) then
!                        Write(IOut,2003) (I,K,I=IFirst,ILast)
!                     else
!                        Write(IOut,2004) (I,K,I=IFirst,ILast)
!                     endif
!                  else 
!                     If ( K .lt. 10 .and. L .lt. 10 ) then
!                        Write(IOut,2005) (I,K,L,I=IFirst,ILast)
!                     Else
!                        Write(IOut,2006) (I,K,L,I=IFirst,ILast)
!                     EndIf
!                  endif
!                  Do I = IFirst, ILAST
!                     JBase = (I*(I-1))/2
!                     Write(IOut,2011) I, (R3D(JBase+J,K,L),J=IFirst,I)
!                  EndDo
!                  Do I = ILAST+1, NCols
!                     JBase = (I*(I-1))/2
!                     Write(IOut,2011) I, (R3D(JBase+J,K,L),J=IFirst,ILast)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMSYM") then
!         N4 = Array%N4
!         IBase=0
!         Do I = 1,N4
!            Do J = 1,I
!               Do K = 1,I
!                  if ( I .eq. K ) then
!                     LMax = J
!                  else
!                     LMax = K
!                  endif
!                  Do L = 1,LMax
!                     IBase=IBase+1
!                     write(Iout,4010) I,J,K,L,Array%ArrR(IBase)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!         If ( IBASE .ne. Array%LR ) then
!            call mqc_error_I('MQC_Print_Array_REAL: Calculated wrong size of Super Symmertric array.', IOut, &
!            'IBASE', IBASE, 'Array%LR', Array%LR )
!         EndIf
!      else if ( Array%Shape.eq."TRISYMMATRICES") then
!         N3 = Array%N3
!         IBase=0
!         Do I = 1,N3
!            Do J = 1,I
!               Do K = 1,J
!                  IBase=IBase+1
!                  write(Iout,4020) I,J,K,Array%ArrR(IBase)
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."3DMATRIX") then
!         N1 = Array%N1
!         N2 = Array%N2
!         N3 = Array%N3
!         R3D(1:N1,1:N2,1:N3)=>Array%ArrR
!         IBase=0
!         Do K = 1,N3
!            Do J = 1,N2
!               Do I = 1,N1
!                  IBase=IBase+1
!                  write(Iout,4020) I,J,K,R3D(I,J,K)
!               EndDo
!            EndDo
!         EndDo
!      else
!         write(IOut,5010) trim(Array%Data_Name)
!         call MQC_PRINT_DATA_TYPE(IOut, Array%Data_Name, Array%NI, Array%NR, &
!              Array%N1, Array%N2, Array%N3, Array%N4, Array%N5, Array%NRI, &
!              Array%NTot, Array%LR, Array%Data_Type, Array%Shape, &
!              Array%Formula)
!      EndIf
!
!      If(PRESENT(Blank_At_Bottom)) then
!        If(Blank_At_Bottom) Write(IOut,1020)
!      EndIf
!      return
!      end Subroutine MQC_Print_Array_REAL
!
!      Subroutine MQC_Print_Array_COMPLEX(IOut, Array,Blank_At_Top, &
!        Blank_At_Bottom)
!      implicit none
!      class(mqc_array_link_list_node),intent(in):: Array
!      integer(kind=int64),intent(in)::iOut
!      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
!      Integer(kind=int64),Parameter::ColWidth=5 
!      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
!      Integer(kind=int64)::N1,N2,N3,N4,K,L,IBase,JBase,LMAX
!      Complex(Kind=real64),contiguous,Pointer::C2D(:,:) =>null()
!      Complex(Kind=real64),contiguous,Pointer::C3D(:,:,:) =>null()
!      Character(len=8) :: LParenthesis
!!
!!     Format statements.
!!
! 1000 Format(/,1x,A)
! 1001 Format(1x,I7,2x,"(",F14.6,",",F14.6,"i)")
! 1010 Format(1x,A,2x,"(",F14.6,",",F14.6,"i)")
! 1020 Format(" ")
! 2001 Format(/,5x,5(20x,I7))
! 2002 Format(/,5x,5(18x,I7,",",I1))
! 2003 Format(/,5x,5(17x,I7,",",I2))
! 2004 Format(/,5x,5(13x,I7,",",I6))
! 2005 Format(/,5x,5(16x,I7,",",I1,",",I1))
! 2006 Format(/,5x,5(4x,I7,",",I7,",",I7))
! 2011 Format(1x,I7," ",A,F13.6,",",F13.6,"i)  ",A,F13.6,",",F13.6,"i)  ",A, &
!           F13.6,",",F13.6,"i)  ",A,F13.6,",",F13.6,"i)  ",A,F13.6,",", &
!           F13.6,"i)")
! 4010 format( "(",I5,",",I5,"|",I5,",",I5,") = (",F14.6,",",F14.6,"i)")
! 4020 format( "(",I5,",",I5,",",I5,") = (",F14.6,",",F14.6,"i)")
! 5010 Format("Warning: Do Not know how do deal with the shape of ",A)
!
!      If(PRESENT(Blank_At_Top)) then
!        If(Blank_At_Top) Write(IOut,1020)
!      EndIf
!      if ( Array%Shape.ne."SCALAR") then
!         Write(IOut,1000) TRIM(Array%Data_Name)
!      EndIf
!
!      if ( Array%Shape.eq."SCALAR") then
!         Write(IOut,1010) TRIM(Array%Data_Name), REAL(Array%ArrC(1)), &
!              AIMAG(Array%ArrC(1))
!      else if ( Array%Shape.eq."VECTOR") then
!         Do I = 1, Array%LR
!            Write(IOut,1001) I, REAL(Array%ArrC(I)), &
!                 AIMAG(Array%ArrC(I))
!         EndDo
!      else if ( Array%Shape.eq."MATRIX") then
!         NRows = Array%N1
!         NCols = Array%N2
!         LParenthesis = "("
!         C2D(1:NRows,1:NCols)=>Array%ArrC
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = 1,NRows
!               Write(IOut,2011) I, (Trim(LParenthesis), &
!                    REAL(C2D(I,J)), &
!                    AIMAG(C2D(I,J)),J=IFirst,ILast)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMMATRIX") then
!         NRows = Array%N2
!         NCols = Array%N2
!         LParenthesis = "("
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = IFirst, ILAST
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Trim(LParenthesis), &
!                    REAL(Array%ArrC(IBase+J)), &
!                    AIMAG(Array%ArrC(IBase+J)),J=IFirst,I)
!            EndDo
!            Do I = ILAST+1, NCols
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Trim(LParenthesis), &
!                    REAL(Array%ArrC(IBase+J)), &
!                    AIMAG(Array%ArrC(IBase+J)),J=IFirst,ILAST)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."MULTISYMMATRICES") then
!         NRows = Array%N2
!         NCols = Array%N2
!         N3 = Array%N3
!         N4 = Array%N4
!         LParenthesis = "("
!         C3D(1:(NCols+1)*NCols/2,1:N3,1:N4)=>Array%ArrC
!         Do L = 1,N4
!            Do K = 1,N3
!               Do IFirst = 1,NCols,ColWidth
!                  ILast = Min(IFirst+ColWidth-1,NCols)
!                  if ( N4 .eq. 1 ) then
!                     if ( K .lt. 10 ) then
!                        Write(IOut,2002) (I,K,I=IFirst,ILast)
!                     else if ( K .lt. 100 ) then
!                        Write(IOut,2003) (I,K,I=IFirst,ILast)
!                     else
!                        Write(IOut,2004) (I,K,I=IFirst,ILast)
!                     endif
!                  else 
!                     If ( K .lt. 10 .and. L .lt. 10 ) then
!                        Write(IOut,2005) (I,K,L,I=IFirst,ILast)
!                     Else
!                        Write(IOut,2006) (I,K,L,I=IFirst,ILast)
!                     EndIf
!                  endif
!                  Do I = IFirst, ILAST
!                     JBase = (I*(I-1))/2
!                     Write(IOut,2011) I, (Trim(LParenthesis), &
!                          REAL(C3D(JBase+J,K,L)), &
!                          AIMAG(C3D(JBase+J,K,L)),J=IFirst,I)
!                  EndDo
!                  Do I = ILAST+1, NCols
!                     JBase = (I*(I-1))/2
!                     Write(IOut,2011) I, (Trim(LParenthesis), &
!                          REAL(C3D(JBase+J,K,L)), &
!                          AIMAG(C3D(JBase+J,K,L)),J=IFirst,ILAST)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMSYM") then
!         N4 = Array%N4
!         IBase=0
!         Do I = 1,N4
!            Do J = 1,I
!               Do K = 1,I
!                  if ( I .eq. K ) then
!                     LMax = J
!                  else
!                     LMax = K
!                  endif
!                  Do L = 1,LMax
!                     IBase=IBase+1
!                     write(Iout,4010) I,J,K,L,REAL(Array%ArrC(IBase)), &
!                          AIMAG(Array%ArrC(IBase))
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."TRISYMMATRICES") then
!         N3 = Array%N3
!         IBase=0
!         Do I = 1,N3
!            Do J = 1,I
!               Do K = 1,J
!                  IBase=IBase+1
!                  write(Iout,4020) I,J,K,REAL(Array%ArrC(IBase)), &
!                       AIMAG(Array%ArrC(IBase))
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."3DMATRIX") then
!         N1 = Array%N1
!         N2 = Array%N2
!         N3 = Array%N3
!         C3D(1:N1,1:N2,1:N3)=>Array%ArrC
!         Do K = 1,N3
!            Do J = 1,N2
!               Do I = 1,N1
!                  write(Iout,4020) I,J,K,REAL(C3D(I,J,K)), &
!                       AIMAG(C3D(I,J,K))
!               EndDo
!            EndDo
!         EndDo
!      else
!         write(IOut,5010) trim(Array%Data_Name)
!         call MQC_PRINT_DATA_TYPE(IOut, Array%Data_Name, Array%NI, Array%NR, &
!              Array%N1, Array%N2, Array%N3, Array%N4, Array%N5, Array%NRI, &
!              Array%NTot, Array%LR, Array%Data_Type, Array%Shape, &
!              Array%Formula)
!      endif
!
!      If(PRESENT(Blank_At_Bottom)) then
!        If(Blank_At_Bottom) Write(IOut,1020)
!      EndIf
!      return
!      end Subroutine MQC_Print_Array_COMPLEX
!
!      Subroutine MQC_Print_Array_RMIXED(IOut, Array,Blank_At_Top, &
!        Blank_At_Bottom)
!      implicit none
!      class(mqc_array_link_list_node),intent(in):: Array
!      integer(kind=int64),intent(in)::iOut
!      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
!      Integer(kind=int64),Parameter::ColWidth=10 
!      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
!      Integer(kind=int64)::N1,N2,N3,N4,K,L,IBase,JBase,LMAX
!      Integer(Kind=int64),contiguous,Pointer::I2D(:,:) =>null()
!      Integer(Kind=int64),contiguous,Pointer::I3D(:,:,:) =>null()
!      Real(Kind=real64),contiguous,Pointer::R2D(:,:) =>null()
!      Real(Kind=real64),contiguous,Pointer::R3D(:,:,:) =>null()
!!
!!     Format statements.
!!
! 1000 Format(/,1x,A)
! 1001 Format(1x,I7,I10,2x,F14.6)
! 1010 Format(1x,A,I10,2x,F14.6)
! 1020 Format(" ")
! 2001 Format(/,5x,10(7x,I7))
! 2002 Format(/,5x,10(5x,I7,",",I1))
! 2003 Format(/,5x,10(4x,I7,",",I2))
! 2004 Format(/,5x,10(I7,",",I6))
! 2005 Format(/,5x,10(3x,I7,",",I1,",",I1))
! 2006 Format(/,5x,10(3x,I7,",",I7,",",I7))
! 2011 Format(1x,I7,10F14.6)
! 4010 format( "(",I5,",",I5,"|",I5,",",I5,") = ",F14.6 )
! 4020 format( "(",I5,",",I5,",",I5,") = ",F14.6 )
! 5010 Format("Warning: Do Not know how do deal with the shape of ",A)
! 7011 Format(1x,I7,10I14)
! 8010 format( "(",I5,",",I5,"|",I5,",",I5,") = ",I14 )
! 8020 format( "(",I5,",",I5,",",I5,") = ",I14 )
!
!      If(PRESENT(Blank_At_Top)) then
!        If(Blank_At_Top) Write(IOut,1020)
!      EndIf
!      if ( Array%Shape.ne."SCALAR") then
!         Write(IOut,1000) TRIM(Array%Data_Name)
!      EndIf
!
!      if ( Array%Shape.eq."SCALAR") then
!         Write(IOut,1010) TRIM(Array%Data_Name), Array%ArrI(1), &
!              Array%ArrR(1)
!      else if ( Array%Shape.eq."VECTOR") then
!         Do I = 1, Array%LR
!            Write(IOut,1001) I, Array%ArrI(I), &
!                 Array%ArrR(I)
!         EndDo
!      else if ( Array%Shape.eq."MATRIX") then
!         NRows = Array%N1
!         NCols = Array%N2
!         I2D(1:NRows,1:NCols)=>Array%ArrI
!         R2D(1:NRows,1:NCols)=>Array%ArrR
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = 1,NRows
!               Write(IOut,2011) I, (R2D(I,J),J=IFirst,ILast)
!               Write(IOut,7011) I, (I2D(I,J),J=IFirst,ILast)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMMATRIX") then
!         NRows = Array%N2
!         NCols = Array%N2
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = IFirst, ILAST
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Array%ArrR(IBase+J),J=IFirst,I)
!               Write(IOut,7011) I, (Array%ArrI(IBase+J),J=IFirst,I)
!            EndDo
!            Do I = ILAST+1, NCols
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Array%ArrR(IBase+J),J=IFirst,ILast)
!               Write(IOut,7011) I, (Array%ArrI(IBase+J),J=IFirst,ILast)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."MULTISYMMATRICES") then
!         NRows = Array%N2
!         NCols = Array%N2
!         N3 = Array%N3
!         N4 = Array%N4
!         I3D(1:(NCols+1)*NCols/2,1:N3,1:N4)=>Array%ArrI
!         R3D(1:(NCols+1)*NCols/2,1:N3,1:N4)=>Array%ArrR
!         Do L = 1,N4
!            Do K = 1,N3
!               Do IFirst = 1,NCols,ColWidth
!                  ILast = Min(IFirst+ColWidth-1,NCols)
!                  if ( N4 .eq. 1 ) then
!                     if ( K .lt. 10 ) then
!                        Write(IOut,2002) (I,K,I=IFirst,ILast)
!                     else if ( K .lt. 100 ) then
!                        Write(IOut,2003) (I,K,I=IFirst,ILast)
!                     else
!                        Write(IOut,2004) (I,K,I=IFirst,ILast)
!                     endif
!                  else 
!                     If ( K .lt. 10 .and. L .lt. 10 ) then
!                        Write(IOut,2005) (I,K,L,I=IFirst,ILast)
!                     Else
!                        Write(IOut,2006) (I,K,L,I=IFirst,ILast)
!                     EndIf
!                  endif
!                  Do I = IFirst, ILAST
!                     JBase = (I*(I-1))/2
!                     Write(IOut,2011) I, (R3D(JBase+J,K,L),J=IFirst,I)
!                     Write(IOut,7011) I, (I3D(JBase+J,K,L),J=IFirst,I)
!                  EndDo
!                  Do I = ILAST+1, NCols
!                     JBase = (I*(I-1))/2
!                     Write(IOut,2011) I, (R3D(JBase+J,K,L),J=IFirst,ILast)
!                     Write(IOut,7011) I, (I3D(JBase+J,K,L),J=IFirst,ILast)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMSYM") then
!         N4 = Array%N4
!         IBase=0
!         Do I = 1,N4
!            Do J = 1,I
!               Do K = 1,I
!                  if ( I .eq. K ) then
!                     LMax = J
!                  else
!                     LMax = K
!                  endif
!                  Do L = 1,LMax
!                     IBase=IBase+1
!                     write(Iout,4010) I,J,K,L,Array%ArrR(IBase)
!                     write(Iout,8010) I,J,K,L,Array%ArrI(IBase)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!         If ( IBASE .ne. Array%LR ) then
!            call mqc_error_I('MQC_Print_Array_REAL: Calculated wrong size of Super Symmertric array.', IOut, &
!            'IBASE', IBASE, 'Array%LR', Array%LR )
!         EndIf
!      else if ( Array%Shape.eq."TRISYMMATRICES") then
!         N3 = Array%N3
!         IBase=0
!         Do I = 1,N3
!            Do J = 1,I
!               Do K = 1,J
!                  IBase=IBase+1
!                  write(Iout,4020) I,J,K,Array%ArrR(IBase)
!                  write(Iout,8020) I,J,K,Array%ArrI(IBase)
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."3DMATRIX") then
!         N1 = Array%N1
!         N2 = Array%N2
!         N3 = Array%N3
!         I3D(1:N1,1:N2,1:N3)=>Array%ArrI
!         R3D(1:N1,1:N2,1:N3)=>Array%ArrR
!         Do K = 1,N3
!            Do J = 1,N2
!               Do I = 1,N1
!                  write(Iout,4020) I,J,K,R3D(I,J,K)
!                  write(Iout,8020) I,J,K,I3D(I,J,K)
!               EndDo
!            EndDo
!         EndDo
!      else
!         write(IOut,5010) trim(Array%Data_Name)
!         call MQC_PRINT_DATA_TYPE(IOut, Array%Data_Name, Array%NI, Array%NR, &
!              Array%N1, Array%N2, Array%N3, Array%N4, Array%N5, Array%NRI, &
!              Array%NTot, Array%LR, Array%Data_Type, Array%Shape, &
!              Array%Formula)
!      EndIf
!
!      return
!      end Subroutine MQC_Print_Array_RMIXED
!
!      Subroutine MQC_Print_Array_CMIXED(IOut, Array,Blank_At_Top, &
!        Blank_At_Bottom)
!      implicit none
!      class(mqc_array_link_list_node),intent(in):: Array
!      integer(kind=int64),intent(in)::iOut
!      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
!      Integer(kind=int64),Parameter::ColWidth=5 
!      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
!      Integer(kind=int64)::N1,N2,N3,N4,K,L,IBase,JBase,LMAX
!      Integer(Kind=int64),contiguous,Pointer::I2D(:,:) =>null()
!      Integer(Kind=int64),contiguous,Pointer::I3D(:,:,:) =>null()
!      Complex(Kind=real64),contiguous,Pointer::C2D(:,:) =>null()
!      Complex(Kind=real64),contiguous,Pointer::C3D(:,:,:) =>null()
!      Character(len=8) :: LParenthesis
!!
!!     Format statements.
!!
! 1000 Format(/,1x,A)
! 1001 Format(1x,I7,I10,2x,"(",F14.6,",",F14.6,"i)")
! 1010 Format(1x,A,I10,2x,"(",F14.6,",",F14.6,"i)")
! 1020 Format(" ")
! 2001 Format(/,5x,5(20x,I7))
! 2002 Format(/,5x,5(18x,I7,",",I1))
! 2003 Format(/,5x,5(17x,I7,",",I2))
! 2004 Format(/,5x,5(13x,I7,",",I6))
! 2005 Format(/,5x,5(16x,I7,",",I1,",",I1))
! 2006 Format(/,5x,5(4x,I7,",",I7,",",I7))
! 2011 Format(1x,I7," ",A,F13.6,",",F13.6,"i)  ",A,F13.6,",",F13.6,"i)  ",A, &
!           F13.6,",",F13.6,"i)  ",A,F13.6,",",F13.6,"i)  ",A,F13.6,",", &
!           F13.6,"i)")
! 4010 format( "(",I5,",",I5,"|",I5,",",I5,") = (",F14.6,",",F14.6,"i)")
! 4020 format( "(",I5,",",I5,",",I5,") = (",F14.6,",",F14.6,"i)")
! 5010 Format("Warning: Do Not know how do deal with the shape of ",A)
! 7011 Format(1x,I7,5I30)
! 8010 format( "(",I5,",",I5,"|",I5,",",I5,") = ",I14 )
! 8020 format( "(",I5,",",I5,",",I5,") = ",I14 )
!
!      If(PRESENT(Blank_At_Top)) then
!        If(Blank_At_Top) Write(IOut,1020)
!      EndIf
!      if ( Array%Shape.ne."SCALAR") then
!         Write(IOut,1000) TRIM(Array%Data_Name)
!      EndIf
!
!      if ( Array%Shape.eq."SCALAR") then
!         Write(IOut,1010) TRIM(Array%Data_Name), Array%ArrI(1), &
!              REAL(Array%ArrC(1)), &
!              AIMAG(Array%ArrC(1))
!      else if ( Array%Shape.eq."VECTOR") then
!         Do I = 1, Array%LR
!            Write(IOut,1001) I, Array%ArrI(I), &
!                 REAL(Array%ArrC(I)), &
!                 AIMAG(Array%ArrC(I))
!         EndDo
!      else if ( Array%Shape.eq."MATRIX") then
!         NRows = Array%N1
!         NCols = Array%N2
!         LParenthesis = "("
!         I2D(1:NRows,1:NCols)=>Array%ArrI
!         C2D(1:NRows,1:NCols)=>Array%ArrC
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = 1,NRows
!               Write(IOut,2011) I, (Trim(LParenthesis), &
!                    REAL(C2D(I,J)), &
!                    AIMAG(C2D(I,J)),J=IFirst,ILast)
!               Write(IOut,7011) I, (I2D(I,J),J=IFirst,ILast)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMMATRIX") then
!         NRows = Array%N2
!         NCols = Array%N2
!         LParenthesis = "("
!         Do IFirst = 1,NCols,ColWidth
!            ILast = Min(IFirst+ColWidth-1,NCols)
!            Write(IOut,2001) (I,I=IFirst,ILast)
!            Do I = IFirst, ILAST
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Trim(LParenthesis), &
!                    REAL(Array%ArrC(IBase+J)), &
!                    AIMAG(Array%ArrC(IBase+J)),J=IFirst,I)
!               Write(IOut,7011) I, (Array%ArrI(IBase+J),J=IFirst,I)
!            EndDo
!            Do I = ILAST+1, NCols
!               IBase = (I*(I-1))/2
!               Write(IOut,2011) I, (Trim(LParenthesis), &
!                    REAL(Array%ArrC(IBase+J)), &
!                    AIMAG(Array%ArrC(IBase+J)),J=IFirst,ILAST)
!               Write(IOut,7011) I, (Array%ArrI(IBase+J),J=IFirst,ILAST)
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."MULTISYMMATRICES") then
!         NRows = Array%N2
!         NCols = Array%N2
!         N3 = Array%N3
!         N4 = Array%N4
!         LParenthesis = "("
!         I3D(1:(NCols+1)*NCols/2,1:N3,1:N4)=>Array%ArrI
!         C3D(1:(NCols+1)*NCols/2,1:N3,1:N4)=>Array%ArrC
!         Do L = 1,N4
!            Do K = 1,N3
!               Do IFirst = 1,NCols,ColWidth
!                  ILast = Min(IFirst+ColWidth-1,NCols)
!                  if ( N4 .eq. 1 ) then
!                     if ( K .lt. 10 ) then
!                        Write(IOut,2002) (I,K,I=IFirst,ILast)
!                     else if ( K .lt. 100 ) then
!                        Write(IOut,2003) (I,K,I=IFirst,ILast)
!                     else
!                        Write(IOut,2004) (I,K,I=IFirst,ILast)
!                     endif
!                  else 
!                     If ( K .lt. 10 .and. L .lt. 10 ) then
!                        Write(IOut,2005) (I,K,L,I=IFirst,ILast)
!                     Else
!                        Write(IOut,2006) (I,K,L,I=IFirst,ILast)
!                     EndIf
!                  endif
!                  Do I = IFirst, ILAST
!                     JBase = (I*(I-1))/2
!                     Write(IOut,2011) I, (Trim(LParenthesis), &
!                          REAL(C3D(JBase+J,K,L)), &
!                          AIMAG(C3D(JBase+J,K,L)),J=IFirst,I)
!                     Write(IOut,7011) I, (I3D(JBase+J,K,L),J=IFirst,I)
!                  EndDo
!                  Do I = ILAST+1, NCols
!                     JBase = (I*(I-1))/2
!                     Write(IOut,2011) I, (Trim(LParenthesis), &
!                          REAL(C3D(JBase+J,K,L)), &
!                          AIMAG(C3D(JBase+J,K,L)),J=IFirst,ILAST)
!                     Write(IOut,7011) I, (I3D(JBase+J,K,L),J=IFirst,ILAST)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."SYMSYM") then
!         N4 = Array%N4
!         IBase=0
!         Do I = 1,N4
!            Do J = 1,I
!               Do K = 1,I
!                  if ( I .eq. K ) then
!                     LMax = J
!                  else
!                     LMax = K
!                  endif
!                  Do L = 1,LMax
!                     IBase=IBase+1
!                     write(Iout,4010) I,J,K,L,REAL(Array%ArrC(IBase)), &
!                          AIMAG(Array%ArrC(IBase))
!                     write(Iout,8010) I,J,K,L,Array%ArrI(IBase)
!                  EndDo
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."TRISYMMATRICES") then
!         N3 = Array%N3
!         IBase=0
!         Do I = 1,N3
!            Do J = 1,I
!               Do K = 1,J
!                  IBase=IBase+1
!                  write(Iout,4020) I,J,K,REAL(Array%ArrC(IBase)), &
!                       AIMAG(Array%ArrC(IBase))
!                  write(Iout,8020) I,J,K,Array%ArrI(IBase)
!               EndDo
!            EndDo
!         EndDo
!      else if ( Array%Shape.eq."3DMATRIX") then
!         N1 = Array%N1
!         N2 = Array%N2
!         N3 = Array%N3
!         I3D(1:N1,1:N2,1:N3)=>Array%ArrI
!         C3D(1:N1,1:N2,1:N3)=>Array%ArrC
!         Do K = 1,N3
!            Do J = 1,N2
!               Do I = 1,N1
!                  write(Iout,4020) I,J,K,REAL(C3D(I,J,K)), &
!                       AIMAG(C3D(I,J,K))
!                  write(Iout,8020) I,J,K,I3D(I,J,K)
!               EndDo
!            EndDo
!         EndDo
!      else
!         write(IOut,5010) trim(Array%Data_Name)
!         call MQC_PRINT_DATA_TYPE(IOut, Array%Data_Name, Array%NI, Array%NR, &
!              Array%N1, Array%N2, Array%N3, Array%N4, Array%N5, Array%NRI, &
!              Array%NTot, Array%LR, Array%Data_Type, Array%Shape, &
!              Array%Formula)
!      endif
!
!      If(PRESENT(Blank_At_Bottom)) then
!        If(Blank_At_Bottom) Write(IOut,1020)
!      EndIf
!      return
!      end Subroutine MQC_Print_Array_CMIXED
!
!
!!=====================================================================
!!
!! Functions for mqc_array_link_list class structures
!!
!!=====================================================================
!      subroutine mqc_array_link_list_add(wavefunction,DataName, &
!           NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut, &
!           Data_Type,Shape,Formula)
!        class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!        character(len=*),intent(in) :: DataName
!        integer(kind=int64),intent(in) :: NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf,IOut
!        logical,intent(in)::Asym
!        character(Len=*),intent(in):: Data_Type, Shape, Formula
!        Integer(kind=int32)::AllocateStatus
!
!        if (.not.associated(wavefunction%Array_Data%first)) then
!           allocate(wavefunction%Array_Data%first, STAT = AllocateStatus)
!           IF (AllocateStatus /= 0) call MQC_Error_A('Not Enough Memory to allocate data structure', IOut, 'DataName', DataName)
!           wavefunction%Array_Data%last => wavefunction%Array_Data%first
!        else
!           allocate(wavefunction%Array_Data%last%next, STAT = AllocateStatus)
!           IF (AllocateStatus /= 0) call MQC_Error_A('Not Enough Memory to allocate data structure', IOut, 'DataName', DataName)
!           wavefunction%Array_Data%last => wavefunction%Array_Data%last%next
!        endif
!        wavefunction%Array_Data%last%Data_Name=trim(DataName)
!        wavefunction%Array_Data%last%NR=NR
!        wavefunction%Array_Data%last%NI=NI
!        wavefunction%Array_Data%last%NRI=NRI
!        wavefunction%Array_Data%last%NTot=NTot
!        wavefunction%Array_Data%last%N1=N1
!        wavefunction%Array_Data%last%N2=N2
!        wavefunction%Array_Data%last%N3=N3
!        wavefunction%Array_Data%last%N4=N4
!        wavefunction%Array_Data%last%N5=N5
!        wavefunction%Array_Data%last%LR=LR
!        wavefunction%Array_Data%last%LenBuf=LenBuf
!        wavefunction%Array_Data%last%Asym=Asym
!        wavefunction%Array_Data%last%Data_Type=trim(Data_Type)
!        wavefunction%Array_Data%last%Shape=trim(Shape)
!        wavefunction%Array_Data%last%Formula=trim(Formula)
!      end subroutine mqc_array_link_list_add
!
!      subroutine mqc_array_link_list_start(self,wavefunction)
!        class(MQC_ALL_Wavefunction),intent(in) :: wavefunction
!        class(mqc_array_link_list),intent(inout) :: self
!
!        self%current => wavefunction%Array_Data%first
!      end subroutine mqc_array_link_list_start
!
!      subroutine mqc_array_link_list_next(self)
!        class(mqc_array_link_list),intent(inout) :: self
!
!        if (associated(self%current)) then
!           self%current => self%current%next
!        else
!           self%current => null()
!        end if
!      end subroutine mqc_array_link_list_next
!
!!=====================================================================
!!
!! Functions for mqc_link_list class structures
!!
!!=====================================================================
!      subroutine mqc_link_list_add(self,Print_Options,Storage_Options, &
!              my_READONLY_list,FileName,iout)
!        class(mqc_link_list),intent(inout) :: self
!        integer(kind=int64),intent(in)::Print_Options
!        integer(kind=int64),intent(in)::Storage_Options
!        class(mqc_READONLY_list),intent(inout):: my_READONLY_list
!        character(len=*),intent(in) :: FileName
!        integer(kind=int64),intent(in)::iOut
!
!        if (.not.associated(self%last)) then
!           allocate(self%first)
!           self%last => self%first
!        else
!           allocate(self%last%next)
!           self%last => self%last%next
!        endif
!
!        self%last%current_wavefunction%FileName=trim(FileName)
!
!        call mqc_read_all_wavefunction( self%last%current_wavefunction, &
!             my_READONLY_list, Print_Options,Storage_Options,iout)
!
!        self%length = self%length + 1
!      end subroutine mqc_link_list_add
!
!      subroutine mqc_link_list_start(self,Status,wavefunction)
!        class(mqc_link_list),intent(inout) :: self
!        integer(kind=int64),intent(out)::Status
!        class(MQC_ALL_Wavefunction),pointer,intent(out)::wavefunction
!
!        self%current => self%first
!        If (associated(self%current)) then
!           wavefunction => self%current%current_wavefunction
!           Status=1
!           if ( wavefunction%Empty ) then
!              call mqc_link_list_next(self,Status,wavefunction)
!           endif
!        Else
!           Status=0
!        Endif
!      end subroutine mqc_link_list_start
!
!      RECURSIVE subroutine mqc_link_list_next(self,Status,wavefunction)
!        class(mqc_link_list),intent(inout)::self
!        integer(kind=int64),intent(out)::Status
!        class(MQC_ALL_Wavefunction),pointer,intent(out)::wavefunction
!
!        if (associated(self%current)) then
!           self%current => self%current%next
!           If (associated(self%current)) then
!              wavefunction => self%current%current_wavefunction
!              Status=1
!              if ( wavefunction%Empty ) then
!!
!! if the wavefunction is empty, skip it and go to the next one
!!
!                 call mqc_link_list_next(self,Status,wavefunction)
!              EndIf
!           Else
!              Status=0
!           EndIf
!        else
!           wavefunction => null()
!           Status=0
!        end if
!      end subroutine mqc_link_list_next
!
!!=====================================================================
!!
!! Functions for mqc_link_list class structures
!!
!!=====================================================================
!      subroutine mqc_READONLY_list_add(self,Data_In)
!        class(mqc_READONLY_list),intent(inout) :: self
!        character(len=*),intent(in) :: Data_In
!
!        if (.not.associated(self%last)) then
!           allocate(self%first)
!           self%last => self%first
!        else
!           allocate(self%last%next)
!           self%last => self%last%next
!        endif
!
!        self%last%Data_Name=trim(Data_In)
!      end subroutine mqc_READONLY_list_add
!
!      subroutine mqc_READONLY_list_start(self,Status,Data_Name)
!        class(mqc_READONLY_list),intent(inout) :: self
!        integer(kind=int64),intent(out)::Status
!        character(len=64),intent(out) :: Data_Name
!
!        self%current => self%first
!        If (associated(self%current)) then
!           Data_Name = self%current%Data_Name
!           Status=1
!        Else
!           Status=0
!        Endif
!      end subroutine mqc_READONLY_list_start
!
!      subroutine mqc_READONLY_list_next(self,Status,Data_Name)
!        class(mqc_READONLY_list),intent(inout)::self
!        integer(kind=int64),intent(out)::Status
!        character(len=64),intent(out) :: Data_Name
!
!        if (associated(self%current)) then
!           If (associated(self%current%next)) then
!              self%current => self%current%next
!              Data_Name = self%current%Data_Name
!              Status=1
!           Else
!              Status=0
!           EndIf
!        else
!           Status=0
!        end if
!      end subroutine mqc_READONLY_list_next
!
!!=====================================================================
!!
!!PROCEDURE MQC_Unformatted_Read_Array
!      subroutine MQC_Unformatted_Read_Array(wavefunction, &
!           my_READONLY_list,Data_Name, &
!           NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf,Asym,Data_Type, &
!           Shape,Formula,IOut,IERROR,Print_Options,Storage_Options)
!!
!!     This Routine is used to look-up a matrix in a unformatted matrix file load
!!     that array into either (OPTIONAL) output dummy MQC_Matrix argument
!!     <matrixOut>, (OPTIONAL) output dummy MQC_Vector argument <vectorOut>, or
!!     (OPTIONAL) output dummy MQC_R4Tensor argument <r4TensorOut>. The character
!!     label for the array of interest is sent to this routine in dummy argument
!!     <label>.
!!
!!     Dummy argument <filename> is optional and is only used if fileinfo
!!     hasn't already been defined using Routine
!!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!!     filename sent is different from the filename associated with object
!!     fileinfo.
!!
!!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!!     called before calling this routine. The expectation is that
!!     MQC_Gaussian_Unformatted_Matrix_Read_Header is also called before this
!!     routine. However, it is also OK to call this routine first. In that case,
!!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!!
!!     H. P. Hratchian, 2017.
!!     L. M. Thompson, 2017
!!
!!     Variable Declarations.
!!
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
!      character(len=*),intent(in)::Data_Name
!      integer(kind=int64),intent(in)::NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf
!      Logical,intent(in)::Asym
!      character(len=*),intent(in)::Data_Type,Shape,Formula
!      integer(kind=int64),intent(in)::IOut
!      integer(kind=int64),intent(inout)::IERROR
!      integer(kind=int64),intent(in)::Print_Options
!      integer(kind=int64),intent(in)::Storage_Options
!!
!!     Temporary local variables used when calling the gauopen routines.
!!     Local temp variables.
!      character(len=256)::errorMsg
!!
!!     Format statements.
!!
! 1010 Format(A,/,' Data Name:',a,/,'  NI =',I9,/,'  NR =',I9,/,'  NRI =',I9,/,&
!          '  NTot =',I9,/, '  LenBuf =',I9,/, '  N1 =',I9,/, '  N2 =',I9,/ &
!          '  N3 =',I9,/, '  N4 =',I9,/, '  N5 =',I9,/, '  LR =',I9, &
!          '  Calculated Type: ',A )
! 2000 format("MQC_Unformatted_Read_Array: Unknown Shape: ",A,/,A,/, &
!           "N1-N5:",5I8)
!!
!! Set the Spin Type, if the information will be stored
!!
!      If ( Storage_Options .NE. 0 ) Then
!         select case(TRIM(Data_Type))
!         case('REAL')
!            If(trim(Data_Name).eq.'BETA MO COEFFICIENTS') then
!               wavefunction%SpinType='spin'
!            Else If(trim(Data_Name).eq.'ALPHA MO COEFFICIENTS') then
!               wavefunction%SpinType='space'
!            EndIf
!         case('COMPLEX')
!            If(trim(Data_Name).eq.'ALPHA MO COEFFICIENTS') then
!               wavefunction%SpinType='general'
!            EndIf
!         end select
!      EndIf
!
!      If ( (Print_Options .NE. 2 .AND. Storage_Options .EQ. 0 ) .OR. &
!           .not.mqc_read_this_MatrixFile_Element(my_READONLY_list,Data_Name )) then
!         Call Rd_Skip(wavefunction%UnitNumber,NTot,LenBuf)
!         return
!      EndIf
!
!      call mqc_array_link_list_add(wavefunction,Data_Name,NR,NI,NRI,NTot, & 
!           N1,N2,N3,N4,N5,LR,LenBuf,Asym,IOut,Data_Type,Shape,Formula)
!!
!      select case(TRIM(Data_Type))
!      case('INTEGER')
!         If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
!              Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
!              Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
!            If( Shape.eq."UNKNOWN") then
!               write(Iout,2000) TRIM(Data_Name),TRIM(Formula),N1,N2,N3,N4,N5
!            EndIf
!            call MQC_Alloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, LR, &
!                 Data_Name, IOut)
!         Else If ( Shape.eq."MATRIX") then
!            call MQC_Alloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, N1*N2, &
!                 Data_Name, IOut)
!         Else If ( Shape.eq."MULTISYMMATRICES") then
!            call MQC_Alloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, &
!                 ((N2+1)*N2*N3*N4)/2, Data_Name, IOut)
!         Else If ( Shape.eq."3DMATRIX") then
!            call MQC_Alloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, &
!                 N1*N2*N3, Data_Name, IOut)
!         Else
!            call MQC_Error_A("Do not know how to allocate array", Iout, &
!                 "Shape",Shape)
!         Endif
!         call Rd_IBuf(wavefunction%unitNumber,NTot,LenBuf, &
!              wavefunction%Array_Data%last%ArrI)
!         if ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction( wavefunction, IOut)
!            If (Storage_Options .EQ. 0 ) Then
!               call MQC_DeAlloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, &
!                    wavefunction%Array_Data%last%Data_Name, IOut)
!            Endif
!         endif
!      case('REAL')
!         If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
!              Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
!              Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
!            If( Shape.eq."UNKNOWN") then
!               write(Iout,2000) TRIM(Data_Name),TRIM(Formula),N1,N2,N3,N4,N5
!            EndIf
!            call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                 LR, Data_Name, IOut)
!         Else If ( Shape.eq."MATRIX") then
!            call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                 N1*N2, Data_Name, IOut)
!         Else If ( Shape.eq."MULTISYMMATRICES") then
!            call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                 ((N2+1)*N2*N3*N4)/2, Data_Name, IOut)
!         Else If ( Shape.eq."3DMATRIX") then
!            call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                 N1*N2*N3, Data_Name, IOut)
!         Else
!            call MQC_Error_A("Do not know how to allocate array", Iout, &
!                 "Shape",Shape)
!         Endif
!         call Rd_RBuf(wavefunction%unitNumber,NTot,LenBuf, &
!              wavefunction%Array_Data%last%ArrR)
!         if ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction( wavefunction, IOut)
!            If (Storage_Options .EQ. 0 ) Then
!               call MQC_DeAlloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                    wavefunction%Array_Data%last%Data_Name, IOut)
!            Endif
!         endif
!      case('COMPLEX')
!         If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
!              Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
!              Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
!            If( Shape.eq."UNKNOWN") then
!               write(Iout,2000) TRIM(Data_Name),TRIM(Formula),N1,N2,N3,N4,N5
!            EndIf
!            call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                 LR, Data_Name, IOut)
!         Else If ( Shape.eq."MATRIX") then
!            call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                 N1*N2, Data_Name, IOut)
!         Else If ( Shape.eq."MULTISYMMATRICES") then
!            call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                 ((N2+1)*N2*N3*N4)/2, Data_Name, IOut)
!         Else If ( Shape.eq."3DMATRIX") then
!            call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                 N1*N2*N3, Data_Name, IOut)
!         Else
!            call MQC_Error_A("Do not know how to allocate array", Iout, &
!                 "Shape",Shape)
!         Endif
!         call Rd_CBuf(wavefunction%unitNumber,NTot,LenBuf, &
!              wavefunction%Array_Data%last%ArrC)
!         if ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction( wavefunction, IOut)
!            If (Storage_Options .EQ. 0 ) Then
!               call MQC_DeAlloc_Complex_Ptr(wavefunction%Array_Data%last%ArrC,&
!                    wavefunction%Array_Data%last%Data_Name, IOut)
!            Endif
!         endif
!      case('MIXED')
!         If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
!              Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
!              Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
!            If( Shape.eq."UNKNOWN") then
!               write(Iout,2000) TRIM(Data_Name),TRIM(Formula),N1,N2,N3,N4,N5
!            EndIf
!            call MQC_Alloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, &
!                 LR, Data_Name, IOut)
!         Else If ( Shape.eq."MATRIX") then
!            call MQC_Alloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, &
!                 N1*N2, Data_Name, IOut)
!         Else If ( Shape.eq."MULTISYMMATRICES") then
!            call MQC_Alloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, &
!                 ((N2+1)*N2*N3*N4)/2, Data_Name, IOut)
!         Else If ( Shape.eq."3DMATRIX") then
!            call MQC_Alloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, &
!                 N1*N2*N3, Data_Name, IOut)
!         Else
!            call MQC_Error_A("Do not know how to allocate array", Iout, &
!                 "Shape",Shape)
!         Endif
!         if ( NR .ge. 0 ) then
!            If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
!                 Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
!                 Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
!               If( Shape.eq."UNKNOWN") then
!                  write(Iout,2000) TRIM(Data_Name),TRIM(Formula),N1,N2,N3,N4,N5
!               EndIf
!               call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                    LR, Data_Name, IOut)
!            Else If ( Shape.eq."MATRIX") then
!               call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                    N1*N2, Data_Name, IOut)
!            Else If ( Shape.eq."MULTISYMMATRICES") then
!               call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                    ((N2+1)*N2*N3*N4)/2, Data_Name, IOut)
!            Else If ( Shape.eq."3DMATRIX") then
!               call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                    N1*N2*N3, Data_Name, IOut)
!            Else
!               call MQC_Error_A("Do not know how to allocate array", Iout, &
!                    "Shape",Shape)
!            Endif
!            call Rd_RMixed(wavefunction%unitNumber,NTot,LenBuf, &
!                 wavefunction%Array_Data%last%ArrI, &
!                 wavefunction%Array_Data%last%ArrR)
!         else
!            If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
!                 Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
!                 Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
!               If( Shape.eq."UNKNOWN") then
!                  write(Iout,2000) TRIM(Data_Name),TRIM(Formula),N1,N2,N3,N4,N5
!               EndIf
!               call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                    LR, Data_Name, IOut)
!            Else If ( Shape.eq."MATRIX") then
!               call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                    N1*N2, Data_Name, IOut)
!            Else If ( Shape.eq."MULTISYMMATRICES") then
!               call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                    ((N2+1)*N2*N3*N4)/2, Data_Name, IOut)
!            Else If ( Shape.eq."3DMATRIX") then
!               call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                    N1*N2*N3, Data_Name, IOut)
!            Else
!               call MQC_Error_A("Do not know how to allocate array", Iout, &
!                    "Shape",Shape)
!            Endif
!            call Rd_CMixed(wavefunction%unitNumber,NTot,LenBuf, &
!                 wavefunction%Array_Data%last%ArrI, &
!                 wavefunction%Array_Data%last%ArrC)
!         endif
!         If ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction( wavefunction, IOut)
!            If (Storage_Options .EQ. 0 ) Then
!               call MQC_DeAlloc_Int_Ptr( wavefunction%Array_Data%last%ArrI, &
!                    wavefunction%Array_Data%last%Data_Name, IOut)
!               If ( NR .ge. 0 ) then
!                  call MQC_DeAlloc_Real_Ptr(wavefunction%Array_Data%last%ArrR,&
!                       wavefunction%Array_Data%last%Data_Name, IOut)
!               Else
!                  call MQC_DeAlloc_Complex_Ptr(wavefunction%Array_Data%last%ArrC,&
!                       wavefunction%Array_Data%last%Data_Name, IOut)
!               EndIf
!            EndIf
!         EndIf
!      case('TENSOR')
!         if(NRI.eq.1) then
!            call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                 LR*NR, Data_Name, IOut)
!            call Rd_2EN(wavefunction%unitNumber,NR,LR,NR*LR,NTot,LenBuf, &
!                 wavefunction%Array_Data%last%ArrR)
!         Else If(NRI.eq.2) then
!            call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                 LR*NR, Data_Name, IOut)
!            call Rd_2EN(wavefunction%unitNumber,NR,LR,NR*LR, &
!                 2*NTot,2*LenBuf,wavefunction%Array_Data%last%ArrR)
!!!                call Rd_2EN(fileinfo%unitNumber,NR,NTot,NR*NTot,2*NTot,2*LenBuf,complexTmp)
!!            complexTmp = reshape(arrayTmp, shape(complexTmp))
!!            call MQC_Matrix_SymmSymmR4Tensor_Put_Complex(r4TensorOut,complexTmp)
!         else
!            if ( NR .ge. 0 ) then
!               call MQC_Alloc_Real_Ptr( wavefunction%Array_Data%last%ArrR, &
!                    LR, Data_Name, IOut)
!               call Rd_RBuf(wavefunction%unitNumber,NTot,LenBuf, &
!                    wavefunction%Array_Data%last%ArrR)
!            else
!               call MQC_Alloc_Complex_Ptr( wavefunction%Array_Data%last%ArrC, &
!                    LR, Data_Name, IOut)
!               call Rd_CBuf(wavefunction%unitNumber,NTot,LenBuf, &
!                    wavefunction%Array_Data%last%ArrC)
!            endif
!         endIf
!         if ( Print_Options .EQ. 2 ) then
!            call mqc_print_Last_Array_In_wavefunction( wavefunction, IOut)
!            If (Storage_Options .EQ. 0 ) Then
!               if(NRI.eq.1) then
!                  call MQC_DeAlloc_Real_Ptr(wavefunction%Array_Data%last%ArrR,&
!                       wavefunction%Array_Data%last%Data_Name, IOut)
!               Else If(NRI.eq.2) then
!                  call MQC_DeAlloc_Complex_Ptr(wavefunction%Array_Data%last%ArrC,&
!                       wavefunction%Array_Data%last%Data_Name, IOut)
!               Else If ( NR .ge. 0 ) then
!                  call MQC_DeAlloc_Real_Ptr(wavefunction%Array_Data%last%ArrR,&
!                       wavefunction%Array_Data%last%Data_Name, IOut)
!               Else
!                  call MQC_DeAlloc_Complex_Ptr(wavefunction%Array_Data%last%ArrC,&
!                       wavefunction%Array_Data%last%Data_Name, IOut)
!               EndIf
!            EndIf
!         endif
!      case default
!!
!! UNKNOWN was treated as a vector, so should never reach here
!         write(IOut,1010) 'UnExpected data type, found:', &
!              TRIM(Data_Name),NI,NR,NRI,NTot,LenBuf,  &
!              N1,N2,N3,N4,N5,LR
!         Call Rd_Skip(wavefunction%UnitNumber,NTot,LenBuf)
!         IERROR = IERROR + 1
!      end select
!!
!      return
!    end subroutine MQC_Unformatted_Read_Array
!
!!=====================================================================
!!
!!PROCEDURE MQC_PRINT_DATA_TYPE
!      Subroutine MQC_PRINT_DATA_TYPE(Iout,DataName,NI,NR,N1,N2,N3,N4,N5,NRI,NTot,LR,DataType,Shape,Formula)
!!
!!     This function returns a character string indicating the type of array
!!     found in a Gaussian matrix file. This is done using NI, NR, N1, N2, N3,
!!     N4, N5 and NRIfrom a matrix header in a Gaussian unformatted matrix file to
!!     determine the type of array the data corresponds to. The return value will
!!     be prepended by "REAL-", "INTEGER-", or "COMPLEX-" and appended by one of
!!     the following:
!!
!!           "VECTOR"           A vector.
!!           "MATRIX"           A matrix that is allocated full (M x N).
!!           "SYMMATRIX"        A symmetric matrix.
!!           "MULTISYMMATRICES" A number of symmetric matrices.
!!
!!     If the input flags do not uniquely identify a known array type, then this
!!     function returns "UNKNOWN".
!!
!!
!!     H. P. Hratchian, 2017.
!!
!!
!!     Variable Declarations.
!!
!      implicit none
!      integer(kind=int64),intent(in)::iOut
!      character(len=*),intent(in)::DataName
!      integer(kind=int64),intent(in)::NI,NR,N1,N2,N3,N4,N5,NRI
!      integer(kind=int64),intent(in)::NTot,LR
!      character(len=*),intent(in)::DataType,Shape,Formula
!      integer(kind=int64)::i,size
!      real(kind=real64)::a
!      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
!      character(len=8)::tab
!!
!!     Format statements.
!!
! 1010 Format( A )
! 1020 Format( A,"  - ",A )
! 1030 Format( 2A )
! 1040 Format( A,"  - Total elements: ",a )
! 1070 format('NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
!        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
! 2010 Format( i15 )
!!
!!     Do the work...
!!
!      tab = char( 9 )
!      write( iout, 1010) DataName
!      write( iout, 1020) tab, trim(DataType)
!      write( OutMath1, 2010) NTot
!      call mqc_StripSpaces(OutMath1)
!      write( iout, 1040) tab, trim(OutMath1)
!      write( iout, 1030) tab, trim(formula)
!!
!      return
!      end subroutine MQC_PRINT_DATA_TYPE
!
!      subroutine MQC_PRINT_SCALAR_SUMMARY( iOut, current_wavefunction)
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
!      integer(kind=int64),intent(in)::iOut
!      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
!      character(len=8)::tab
!!
!!     Format statements.
!!
! 2010 Format( i15 )
! 2020 Format ( "***************************",/," Scalar Data ",/, &
!           "***************************")
! 2030 Format( A,/,A,"  - Value: ",A,/,A,"  - ",A,/,A,"  - SCALAR" )
! 2040 Format( A,/,A,"  - Value: ",A,/,A,"  - ",A,/,A,"  - SCALAR",/,A,"  - ",A )
!
!      tab = char( 9 )
!      write( IOut, 2020)
!      write( OutMath1, 2010 ) current_wavefunction%natoms
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Number of Atoms", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%nbasis
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Number of Basis Functions", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%nbasisUse
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Number of Linearly Independent Basis Functions", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%iCharge
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Molecular Charge", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%Multiplicity
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Spin Multiplicity", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%nelectrons
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Number of Electrons", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%IOpCl
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2040) 'Closed/open-shell flag', tab, trim(OutMath1), tab, "INTEGER", tab, tab, &
!      'Set if the matrix element file is written after'//NEW_LINE('A')//tab//' &
!         & an initial guess or the SCF has completed '//NEW_LINE('A')//tab//' &
!         & (otherwise it is -1, meaning unspecified).'
!
!      write( OutMath1, 2010 ) current_wavefunction%ICGU
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2040) "Encoded Indicator of whether the calculation is complex and/or GHF", &
!           tab, trim(OutMath1), tab, "INTEGER", tab, tab, &
!           'Its three-digit value is interpreted as klm,'//NEW_LINE('A')//tab//' &
!           & where k is 1 for the spin-aligned case and'//NEW_LINE('A')//tab//' &
!           & 2 for GHF; l is 1 for real and 2 for complex;'//NEW_LINE('A')//tab//' &
!           & and m is 1 for RHF/GHF and 2 for UHF (i.e., '//NEW_LINE('A')//tab//' &
!           & 1 vs. 2 spin blocks). When k=2, then NBasis is'//NEW_LINE('A')//tab//' &
!           & the number of spatial basis functions, but the'//NEW_LINE('A')//tab//' &
!           & operator matrices are over the spin orbital'//NEW_LINE('A')//tab//' &
!           & basis and hence have dimension k*NBasis.'
!
!      write( OutMath1, 2010 ) current_wavefunction%NFC
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Number of Frozen Core Orbitals", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%NFV
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Number of Frozen Virtual Orbitals", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%NShlAO
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2040) "Number of Contracted Shells of AO Basis Functions", &
!           tab, trim(OutMath1), tab, "INTEGER", tab, tab, &
!           "Needed if shell data is provided."
!
!      write( OutMath1, 2010 ) current_wavefunction%NPrmAO
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Number of Primitive AO Shells", tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%NShlDB
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2040) "Number of Contracted Shells of Density Fitting Functions", &
!           tab, trim(OutMath1), tab, "INTEGER", tab, tab, &
!           "Needed if fitting shell data is provided"
!
!      write( OutMath1, 2010 ) current_wavefunction%NPrmDB
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Number of Primitive Density Fitting Shells", &
!           tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%NBTot
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2030) "Total Number of Bonds In Connectivity Data, If Any", &
!           tab, trim(OutMath1), tab, "INTEGER", tab
!
!      write( OutMath1, 2010 ) current_wavefunction%ITran
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2040) "ITRAN", tab, trim(OutMath1), tab, &
!           "INTEGER", tab, tab, &
!           'If no MO integrals were stored, ITran=4 if'//NEW_LINE('A')//tab//' &
!           & only MOs involving at least one occupied'//NEW_LINE('A')//tab//' &
!           & orbital were stored, or ITran=5 if a full'//NEW_LINE('A')//tab//' &
!           & transformation was done.'
!      return
!    end subroutine MQC_PRINT_SCALAR_SUMMARY
!
!      subroutine MQC_PRINT_DEFAULT_VECTOR_SUMMARY( iOut, current_wavefunction)
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
!      integer(kind=int64),intent(in)::iOut
!      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
!      character(len=8)::tab
!!
!!     Format statements.
!!
! 2010 Format( i15 )
! 2050 Format ( "***************************",/," Non-Scalar Data ",/, &
!              "***************************")
! 2060 Format( A,/,A,"  - ",A,/,A,"  - Total elements: ",A,/,A,"  - VECTOR(",a,")" )
! 2070 Format( A,/,A,"  - ",A,/,A,"  - Total elements: ",A,/,A,"  - VECTOR(",a,")",/,A,"  - ",A )
! 2080 Format( A,/,A,"  - ",A,/,A,"  - Total elements: ",A,/,A,"  - 2-D MATRIX allocated full (",a," x ",a,")" )
!
!      tab = char( 9 )
!      write( IOut, 2050)
!! current_wavefunction%atomicNumbers
!      write( OutMath2, 2010 ) current_wavefunction%natoms
!      call mqc_StripSpaces(OutMath2)
!      write( IOut, 2060) "ATOMIC NUMBERS", tab, "INTEGER", tab, trim(OutMath2), tab, trim(OutMath2)
!
!! current_wavefunction%atomTypes
!      write( OutMath2, 2010 ) current_wavefunction%natoms
!      call mqc_StripSpaces(OutMath2)
!      write( IOut, 2060) "ATOM TYPES", tab, "INTEGER", tab, trim(OutMath2), tab, trim(OutMath2)
!
!! current_wavefunction%atomicCharges
!      write( OutMath2, 2010 ) current_wavefunction%natoms
!      call mqc_StripSpaces(OutMath2)
!      write( IOut, 2070) "ATOM CHARGES", tab, "REAL", tab, trim(OutMath2), tab, trim(OutMath2), tab, &
!           'Nuclear charges may be different from atomic'//NEW_LINE('A')//tab//' &
!           & numbers if ECPs were used.'
!
!! current_wavefunction%cartesians
!      write( OutMath1, 2010 ) current_wavefunction%natoms*3
!      call mqc_StripSpaces(OutMath1)
!      write( OutMath2, 2010 ) current_wavefunction%natoms
!      call mqc_StripSpaces(OutMath2)
!      write( OutMath3, 2010 ) 3
!      call mqc_StripSpaces(OutMath3)
!      write( IOut, 2080) "CARTESIAN NUCLEAR COORDINATES(BOHR)", tab, "REAL", &
!           tab, trim(OutMath1), tab, trim(OutMath2), trim(OutMath3)
!
!! First half of current_wavefunction%basisFunction2Atom
!      write( OutMath1, 2010 ) current_wavefunction%nbasis
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2060) "MAP FROM BASIS FUNCTIONS TO ATOMS", tab, &
!           "INTEGER", tab, trim(OutMath1), tab, trim(OutMath1)
!
!! Second half half of current_wavefunction%basisFunction2Atom
!      write( OutMath1, 2010 ) current_wavefunction%nbasis
!      call mqc_StripSpaces(OutMath1)
!      write( IOut, 2060) "TYPE FLAG FOR EACH BASIS FUNCTION", tab, &
!           "INTEGER", tab, trim(OutMath1), tab, trim(OutMath1)
!
!!           current_wavefunction%atomicWeights
!      write( OutMath2, 2010 ) current_wavefunction%natoms
!      call mqc_StripSpaces(OutMath2)
!      write( IOut, 2060) "ATOM WEIGHTS", tab, "REAL", tab, &
!           trim(OutMath2), tab, trim(OutMath2)
!      return
!    end subroutine MQC_PRINT_DEFAULT_VECTOR_SUMMARY
!!=====================================================================
!!
!
!      subroutine MQC_PRINT_SCALAR_HEAD( iOut, wavefunction)
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      integer(kind=int64),intent(in)::iOut
!!
!!     Format statements.
!!
! 1000 Format("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&",/, &
!           1x,A,/,"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
! 1010 Format(1x,A,2x,I14)
!
!      write( IOut, 1000) "Contents of Matrix Element File"
!      write( IOut, 1010) "NUMBER OF ATOMS", wavefunction%natoms
!      write( IOut, 1010) "NUMBER OF BASIS FUNCTIONS", wavefunction%nbasis
!      write( IOut, 1010) "NUMBER OF LINEARLY INDEPENDENT BASIS FUNCTIONS", &
!           wavefunction%nbasisUse
!      write( IOut, 1010) "MOLECULAR CHARGE", wavefunction%iCharge
!      write( IOut, 1010) "SPIN MULTIPLICITY", wavefunction%Multiplicity
!      write( IOut, 1010) "NUMBER OF ELECTRONS", wavefunction%nelectrons
!      write( IOut, 1010) 'CLOSED/OPEN-SHELL FLAG', wavefunction%IOpCl
!      write( IOut, 1010) "ENCODED INDICATOR OF WHETHER THE CALCULATION IS COMPLEX AND/OR GHF", &
!           wavefunction%ICGU
!      write( IOut, 1010) "NUMBER OF FROZEN CORE ORBITALS", wavefunction%NFC
!      write( IOut, 1010) "NUMBER OF FROZEN VIRTUAL ORBITALS", wavefunction%NFV
!      write( IOut, 1010) "NUMBER OF CONTRACTED SHELLS OF AO BASIS FUNCTIONS", &
!           wavefunction%NShlAO
!      write( IOut, 1010) "NUMBER OF PRIMITIVE AO SHELLS", wavefunction%NPrmAO
!      write( IOut, 1010) "NUMBER OF CONTRACTED SHELLS OF DENSITY FITTING FUNCTIONS", &
!           wavefunction%NShlDB
!      write( IOut, 1010) "NUMBER OF PRIMITIVE DENSITY FITTING SHELLS", &
!           wavefunction%NPrmDB
!      write( IOut, 1010) "TOTAL NUMBER OF BONDS IN CONNECTIVITY DATA", &
!           wavefunction%NBTot
!      write( IOut, 1010) "ITRAN", wavefunction%ITran
!      return
!    end subroutine MQC_PRINT_SCALAR_HEAD
!
!    Subroutine MQC_Create_Wavefunctions( Print_Options, Storage_Options, &
!         my_READONLY_list, mqc_input_file_list, Iout)
!      implicit none
!
!      integer(kind=int64),intent(in):: Print_Options
!      integer(kind=int64),intent(in):: Storage_Options
!      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
!      class(mqc_link_list),intent(inout) :: mqc_input_file_list
!      integer(kind=int64),intent(in)::iOut
!
!      character(len=1024)::program_name
!      character(len=1024)::no_program_name
!      character(len=1024)::FileName
!      integer(kind=int64) cmd_count, I
!      integer(kind=int64) status
!
!!
!!     Format statements.
!!
! 1010 Format( A )
!!
!!     Do the work...
!!
!      cmd_count = COMMAND_ARGUMENT_COUNT()
!      if ( cmd_count .eq. 0 ) then
!         write( iout, 1010 ) "Error in number of argments to this program."
!         write( iout, 1010 ) "Arguments should be:"
!         write( iout, 1010 ) "     - Full path to a Gaussian file"
!         write( iout, 1010 ) "        - This can either be FormChk File"
!         write( iout, 1010 ) "        - This can either be a Matrix Element File"
!         write( iout, 1010 ) "        - This can be a Gaussian input file"
!         write( iout, 1010 ) "     - Optional, only required if the first "
!         write( iout, 1010 ) "          argument is a Gaussian input file."
!         write( iout, 1010 ) "          The name of the binary for the "
!         write( iout, 1010 ) "          program that generates the"
!         write( iout, 1010 ) "          FormChk/MatrixElementFile, most likely"
!         write( iout, 1010 ) "          'g16'. The environment should be"
!         write( iout, 1010 ) "          setup to run the program."
!         call mqc_error_i( "No arguments to this program", iout)
!      endif
!!
!!     Create a list of Matrix Element Files.  Use Matrix Element Files supplied.
!!     If Gaussian input has been supplied:
!!     1. Modify the input include Matrix Element File options
!!     2. Run the job to generate Matrix Element Files
!!     If the name of the Gaussian binary is not given, and Gaussian is not
!!     in the path, this will fail.
!!
!      if ( cmd_count.eq. 1 ) then
!         no_program_name = 'Unknown'
!         call get_command_argument(1, FileName)
!         call mqc_create_File_Name_lists_F2C(FileName, no_program_name, iout)
!      else
!         call get_command_argument(cmd_count, program_name)
!         call mqc_File_Is_Executable_F2C( program_name, Status, iout )
!         if ( Status .eq. 1 ) then
!            do i = 1, cmd_count-1
!               call get_command_argument(i, FileName)
!               call mqc_create_File_Name_lists_F2C(FileName,Program_name,iout)
!            end do
!         else
!            no_program_name = 'Unknown'
!            do i = 1, cmd_count
!               call get_command_argument(i, FileName)
!               call mqc_create_File_Name_lists_F2C(FileName,no_program_name,iout)
!            end do
!         endif
!      endif
!#ifdef MATRIXFILE
!! 
!! Figure out is a Matrix Element File name was supplied.
!! Or if an input file has been supplied, and then run the job.
!! Whatever has happened, get a list of the Matrix Element Files
!!
!      call  mqc_get_MatrixFile_Name_F2C(FileName, iout)
!      if ( trim(FileName).eq.'Done with MatrixFiles' ) then
!         call mqc_error_I('No Matrix Element Files supplied.  (Remove this when FormCHK files are implimented', IOut)
!      endif
!!  
!      Do While( trim(FileName).ne.'Done with MatrixFiles' )
!!
!! Read the MatrixFile and put data into a wavefunction structure
!!
!         call mqc_link_list_add(mqc_input_file_list,Print_Options, &
!              Storage_Options, my_READONLY_list,FileName,iout)
!!
!! Get the next MatrixFile Name
!!
!         call mqc_get_MatrixFile_Name_F2C(FileName, iout)
!      EndDo
!#endif
!      return
!    end Subroutine MQC_Create_Wavefunctions
!
!!     PROCEDURE MQC_FULL_GET_NUCLEAR_REPULSION
!      Function MQC_FULL_Wave_To_Nuclear_Repulsion(wavefunction,IOut) Result(Vnn)
!!
!!     This function returns the nuclear repulsion energy given the
!!     molecular data object.
!!
!!     Lee M. Thompson, 2017.
!!
!!
!      Implicit None
!      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
!      real(kind=real64) Vnn
!      type(mqc_array_link_list_node),pointer::Coord
!      type(mqc_array_link_list_node),pointer::Charges
!      Real::Zero=0.0d0
!      Real::DX,DY,DZ
!      Integer::IOut,I,J,ICount,JCount
!      Integer(kind=int64)::NAtoms
!      Integer(kind=int64)::Found_Cart, Found_Charge
!!
!      Vnn = Zero
!      NAtoms = wavefunction%NAtoms
!      Found_Cart = MQC_Find_mqc_array_link_list_node(wavefunction,&
!           "CARTESIAN NUCLEAR COORDINATES(BOHR)", Coord, iout)
!      Found_Charge = MQC_Find_mqc_array_link_list_node(wavefunction,&
!           "ATOM CHARGES", Charges, iout)
!      If ( Found_Cart .ne. 2 .or. Found_Charge .ne. 2 ) then
!         call mqc_error_I('MQC_Wave_To_Nuclear_Repulsion: Did not find all compnents of Nuclear Replusion', IOut, &
!         'Found_Cart', Found_Cart, 'Found_Charge',Found_Charge )
!      EndIf
!      ICount = 0
!      Do I = 1,NAtoms-1
!         JCount = I*3
!         Do J = I+1,NAtoms
!            DX = Coord%ArrR(1+ICount) - Coord%ArrR(1+JCount)
!            DY = Coord%ArrR(2+ICount) - Coord%ArrR(2+JCount)
!            DZ = Coord%ArrR(3+ICount) - Coord%ArrR(3+JCount)
!            Vnn = Vnn + Charges%ArrR(I) * Charges%ArrR(J) / &
!                 SQRT(DX*DX+DY*DY+DZ*DZ)
!            JCount = JCount+3
!         EndDo
!         ICount = ICount+3
!      EndDo
!      End Function MQC_FULL_Wave_To_Nuclear_Repulsion
!
!      Function mqc_read_this_MatrixFile_Element(my_READONLY_list,Data_Name ) &
!           Result( ReadIt )
!      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
!      character(len=*),intent(in)::Data_Name
!      logical ReadIt
!      integer(kind=int64) Status
!      character(len=64)::List_Name
!
!      If (.not.associated(my_READONLY_list%first)) then
!         ReadIt = .true.
!         return
!      else
!         call mqc_READONLY_list_start(my_READONLY_list,Status,List_Name)
!         Do While( Status.eq.1 )
!            If ( trim(Data_Name) .eq. trim(List_Name) ) then
!               ReadIt = .true.
!               return
!            EndIf
!            call mqc_READONLY_list_next(my_READONLY_list,Status,List_Name)
!         EndDo
!      EndIf
!      ReadIt = .false.
!
!      Return
!    end Function mqc_read_this_MatrixFile_Element
!
!      subroutine mqc_ALL_wavefunction_to_est_wavefunction( &
!           current_wavefunction, IOut)
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
!      integer(kind=int64),intent(in)::iOut
!      type(mqc_array_link_list_node),pointer::Alpha
!      type(mqc_array_link_list_node),pointer::Beta
!      type(mqc_array_link_list_node),pointer::AlphaBeta
!      type(mqc_array_link_list_node),pointer::BetaAlpha
!
!      Integer(kind=int64):: Found_Alpha,Found_Beta,Found_AlphaBeta
!      Integer(kind=int64):: Found_BeatAlpha
!
!      Integer(kind=int64):: NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,NRI,LR
!      logical ASym
!      Character(Len=64)::CBuf
!      Character(Len=2048)::DataType,Shape,Formula
!      logical EOF
!      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
!      integer(kind=int32)::I
!      integer(kind=int64)::IERROR, NAtoms, NBasis
!      Type(MQC_Matrix)::Alpha_Matrix
!
!      current_wavefunction%est_wavefunction%nBasis = current_wavefunction%NBasis
!      current_wavefunction%est_wavefunction%nAlpha = current_wavefunction%NAlpha
!      current_wavefunction%est_wavefunction%nBeta = current_wavefunction%NBeta
!      current_wavefunction%est_wavefunction%nElectrons = current_wavefunction%NElectrons
!      current_wavefunction%est_wavefunction%charge = current_wavefunction%iCharge
!      current_wavefunction%est_wavefunction%multiplicity = current_wavefunction%Multiplicity
!
!!
!!     Integer(kind=int64)::NAtoms,nbasisUse,Len12L
!!     Integer(kind=int64)::ICGU,NFC,NFV,IDum9
!!
!
!
!!        if(fileinfo%isRestricted()) then
!!        elseIf(fileinfo%isGeneral()) then
!
!      call mqc_Array_to_Integral( current_wavefunction, &
!           "ALPHA MO COEFFICIENTS", "BETA MO COEFFICIENTS", &
!           'mo coefficients', &
!           current_wavefunction%est_wavefunction%MO_Coefficients, IOut )
!
!      call mqc_Array_to_Integral( current_wavefunction, &
!           "CORE HAMILTONIAN ALPHA", "CORE HAMILTONIAN BETA", &
!           'core hamiltonian', &
!           current_wavefunction%est_wavefunction%Core_Hamiltonian, IOut )
!
!      call mqc_Array_to_Integral( current_wavefunction, &
!           "ALPHA FOCK MATRIX", "BETA FOCK MATRIX", &
!           'fock', &
!           current_wavefunction%est_wavefunction%Fock_Matrix, IOut )
!
!      call mqc_Array_to_Integral( current_wavefunction, &
!           "ALPHA SCF DENSITY MATRIX", "BETA SCF DENSITY MATRIX", &
!           'density', &
!           current_wavefunction%est_wavefunction%Density_Matrix, IOut )
!
!      call mqc_Array_to_Integral( current_wavefunction, &
!           "OVERLAP", "There is no BETA", &
!           'overlap', &
!           current_wavefunction%est_wavefunction%overlap_Matrix, IOut )
!!
!!      Call mqc_Array_to_Eigenvalues( current_wavefunction, &
!!           'ALPHA ORBITAL ENERGIES', 'BETA ORBITAL ENERGIES', &
!!           'mo energies', &
!!           current_wavefunction%est_wavefunction%mo_energies, IOut )
!!
!
!
!!          call mqc_gaussian_ICGU(fileInfo%ICGU,est_wavefunction%wf_type,est_wavefunction%wf_complex)
!
!
!      return
!      end subroutine mqc_ALL_wavefunction_to_est_wavefunction
!
!      subroutine mqc_Array_to_Integral( current_wavefunction, &
!           Alpha_String, Beta_String, Integral_Type, Integral_Result, &
!           IOut)
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
!      character(len=*),intent(in)::Alpha_String
!      character(len=*),intent(in)::Beta_String
!      character(len=*),intent(in)::Integral_Type
!      type(MQC_SCF_Integral),intent(out):: Integral_Result
!      integer(kind=int64),intent(in)::IOut
!      type(mqc_array_link_list_node),pointer::Alpha_Array, Beta_Array
!      type(mqc_array_link_list_node),pointer::AlphaBeta_Array, BetaAlpha_Array
!      Type(MQC_Matrix)::Alpha_Matrix,Beta_Matrix
!      Type(MQC_Matrix)::AlphaBeta_Matrix,BetaAlpha_Matrix
!
!      Integer(kind=int64):: Found_Alpha,Found_Beta,Found_AlphaBeta
!      Integer(kind=int64):: Found_BeatAlpha
!      Integer(kind=int64):: N1,N2,NBasis
!
!      Found_Beta=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              Beta_String, Beta_Array, iout)
!      Found_Alpha=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              Alpha_String, Alpha_Array, iout)
!! If the data has not been saved, return
!      if ( Found_Alpha .eq. 0 ) return
!      N1 = current_wavefunction%Array_Data%last%N1
!      N2 = current_wavefunction%Array_Data%last%N2
!      If ( Found_Beta .eq. 0 ) then
!         If ( Found_Alpha .eq. 2 ) then
!            Call MQC_Copy_Array_To_Matrix(Alpha_Array, N1, N2, &
!                 Alpha_Matrix, IOut)
!            call mqc_integral_allocate(Integral_Result, &
!                 Integral_Type,current_wavefunction%SpinType,Alpha_Matrix)
!         Else
!            Call MQC_Copy_Array_To_Matrix(Alpha_Array, N1, N2, &
!                 Alpha_Matrix, IOut)
!            nBasis = current_wavefunction%nBasis
!            call mqc_matrix_spinBlockGHF(Alpha_Matrix)
!! %mat takes you to MQC_Matrix_Matrix_At
!            Beta_Matrix = Alpha_Matrix%mat([nBasis+1,-1],[nBasis+1,-1])
!            BetaAlpha_Matrix = Alpha_Matrix%mat([1,nBasis],[nBasis+1,-1])
!            AlphaBeta_Matrix = Alpha_Matrix%mat([nBasis+1,-1],[1,nBasis])
!            Alpha_Matrix = Alpha_Matrix%mat([1,nBasis],[1,nBasis])
!            call mqc_integral_allocate(Integral_Result, &
!                 Integral_Type,current_wavefunction%SpinType, &
!                 Alpha_Matrix,Beta_Matrix, &
!                 AlphaBeta_Matrix,BetaAlpha_Matrix)
!         EndIf
!      Else
!         Call MQC_Copy_Array_To_Matrix( Alpha_Array, N1, N2, Alpha_Matrix, &
!              IOut)
!         Call MQC_Copy_Array_To_Matrix( Beta_Array, N1, N2, Beta_Matrix, &
!              IOut)
!         call mqc_integral_allocate(Integral_Result, &
!              Integral_Type,current_wavefunction%SpinType, &
!              Alpha_Matrix,Beta_Matrix)
!      Endif
!      return
!      end subroutine mqc_Array_to_Integral
!
!      subroutine mqc_Array_to_Eigenvalues( current_wavefunction, &
!           Alpha_String, Beta_String, Eigenvalue_Type, Eigenvalue_Result, &
!           IOut)
!      implicit none
!      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
!      character(len=*),intent(in)::Alpha_String
!      character(len=*),intent(in)::Beta_String
!      character(len=*),intent(in)::Eigenvalue_Type
!      type(MQC_SCF_Eigenvalues),intent(out):: Eigenvalue_Result
!      integer(kind=int64),intent(in)::IOut
!      type(mqc_array_link_list_node),pointer::Alpha_Array, Beta_Array
!      Type(MQC_Vector)::Alpha_Vector,Beta_Vector
!
!      Integer(kind=int64):: Found_Alpha,Found_Beta
!      Integer(kind=int64):: N1,NBasis
!
!      Found_Beta=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              Beta_String, Beta_Array, iout)
!      Found_Alpha=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              Alpha_String, Alpha_Array, iout)
!! If the data has not been saved, return
!      if ( Found_Alpha .eq. 0 ) return
!#ifdef DEBUG 
!      write(6,3457) trim(Alpha_Array%Data_Name),N1
!3457  format( "Create E Vector ",A,2I5)
!#endif
!      N1 = current_wavefunction%Array_Data%last%N1
!      If ( trim(current_wavefunction%SpinType).eq.'space' ) then
!         Call MQC_Copy_Array_To_Vector(Alpha_Array, N1, Alpha_Vector, IOut)
!         call mqc_eigenvalues_allocate(Eigenvalue_Result, &
!              Eigenvalue_Type,current_wavefunction%SpinType,Alpha_Vector)
!      Else If ( trim(current_wavefunction%SpinType).eq.'spin' ) then
!         Call MQC_Copy_Array_To_Vector(Alpha_Array, N1, Alpha_Vector, IOut)
!         Call MQC_Copy_Array_To_Vector(Beta_Array, N1, Beta_Vector, IOut)
!         call mqc_eigenvalues_allocate(Eigenvalue_Result, &
!              Eigenvalue_Type,current_wavefunction%SpinType,Alpha_Vector, &
!              Beta_Vector)
!      Else If ( trim(current_wavefunction%SpinType).eq.'general' ) then
!         nBasis = current_wavefunction%nBasis
!         Call mqc_matrix_spinBlockGHF(Alpha_Vector)
!         Beta_Vector = Alpha_Vector%vat(nBasis+1,-1)
!         Alpha_Vector = Alpha_Vector%vat(1,nBasis)
!         Call mqc_eigenvalues_allocate(Eigenvalue_Result, &
!              Eigenvalue_Type,current_wavefunction%SpinType,Alpha_Vector, &
!              Beta_Vector)
!      Else
!         call MQC_Error_A( 'Unknown Spin type ', IOut, &
!              'current_wavefunction%SpinType', current_wavefunction%SpinType )
!      Endif
!      return
!      end subroutine mqc_Array_to_Eigenvalues
!
!      Subroutine MQC_Copy_Array_To_Vector( Alpha_Array, N1, Alpha_Vector, &
!           IOut )
!      type(mqc_array_link_list_node),pointer,intent(in)::Alpha_Array
!      Integer(kind=int64),intent(in):: N1
!      Type(MQC_Vector),intent(inout)::Alpha_Vector
!      integer(kind=int64),intent(in)::IOut
!      Integer(kind=int64)::ABSN1
!      Integer(kind=int64)::I
!      Type(MQC_Scalar)::Temp
!
!      If (associated(Alpha_Array%ArrC)) then
!         ABSN1 = Abs(N1)
!         Call MQC_Allocate_Vector(ABSN1,Alpha_Vector,'Complex')
!         Do I=1,ABSN1
!!            Alpha_Vector%VecC(I) = Alpha_Array%ArrC(I)
!            Temp = Alpha_Array%ArrC(I)
!            Call Alpha_Vector%put(Temp,I)
!         EndDo
!      Else If (associated(Alpha_Array%ArrR)) then
!         Call MQC_Allocate_Vector(N1,Alpha_Vector,'Real')
!         Do I=1,N1
!!            Alpha_Vector%VecR(I) = Alpha_Array%ArrR(I)
!            Temp = Alpha_Array%ArrR(I)
!            Call Alpha_Vector%put(Temp,I)
!         EndDo
!      Else If (associated(Alpha_Array%ArrI)) then
!         Call MQC_Allocate_Vector(N1,Alpha_Vector,'Integer')
!         Do I=1,N1
!!            Alpha_Vector%VecI(I) = Alpha_Array%ArrI(I)
!            Temp = Alpha_Array%ArrI(I)
!            Call Alpha_Vector%put(Temp,I)
!         EndDo
!      EndIf
!
!      Return
!      end Subroutine MQC_Copy_Array_To_Vector
!
!      Subroutine MQC_Copy_Array_To_Matrix( Alpha_Array, N1, N2, Alpha_Matrix, &
!           IOut )
!      type(mqc_array_link_list_node),pointer,intent(in)::Alpha_Array
!      Integer(kind=int64),intent(in):: N1,N2
!      Type(MQC_Matrix),intent(inout)::Alpha_Matrix
!      integer(kind=int64),intent(in)::IOut
!      Real(Kind=real64),contiguous,Pointer::R2D(:,:) =>null()
!      Complex(Kind=real64),contiguous,Pointer::C2D(:,:) =>null()
!      Integer(kind=int64)::I, J, IJ
!      Type(MQC_Scalar)::Temp
!
!      If ( N1.GT.0 .AND. N2.GT.0 ) Then
!         If (associated(Alpha_Array%ArrC)) then
!            C2D(1:N1,1:N2)=>Alpha_Array%ArrC
!            Call MQC_Allocate_Matrix(N1,N2,Alpha_Matrix,'Complex','StorFull')
!            Do J=1,N2
!               Do I=1,N1
!!                  Alpha_Matrix%MatC(I,J) = C2D(I,J)
!                  Temp = C2D(I,J)
!                  call Alpha_Matrix%put(Temp,I,J)
!               EndDo
!            EndDo
!         Else
!            R2D(1:N1,1:N2)=>Alpha_Array%ArrR
!            Call MQC_Allocate_Matrix(N1,N2,Alpha_Matrix,'Real','StorFull')
!            Do J=1,N2
!               Do I=1,N1
!!                  Alpha_Matrix%MatR(I,J) = R2D(I,J)
!                  Temp = R2D(I,J)
!                  call Alpha_Matrix%put(temp,I,J) 
!               EndDo
!            EndDo
!         EndIf
!      Else If ( N1.LT.0 .AND. N2.GT.0 .AND. -N1.EQ.N2 ) Then
!         If (associated(Alpha_Array%ArrC)) then
!            C2D(1:((N2+1)*N2)/2,1:1)=>Alpha_Array%ArrC
!            Call MQC_Allocate_Matrix(N2,N2,Alpha_Matrix,'Complex','StorFull')
!            IJ=0
!            Do J=1,N2
!               Do I=1,J
!                  IJ=IJ+1
!#ifdef DEBUG 
!      write(6,2345) trim(Alpha_Array%Data_Name),I,J,IJ,C2D(IJ,1)
! 2345 format( a, 3I5,2E20.10,"i" )
!#endif
!!                  Alpha_Matrix%MatC(I,J) = C2D(IJ,1)
!!                  Alpha_Matrix%MatC(J,I) = C2D(IJ,1)
!                   Temp = C2D(IJ,1)
!                   call Alpha_Matrix%put(temp,I,J)
!                   call Alpha_Matrix%put(temp,J,I)
!               EndDo
!            EndDo
!         Else
!            R2D(1:((N2+1)*N2)/2,1:1)=>Alpha_Array%ArrR
!            Call MQC_Allocate_Matrix(N2,N2,Alpha_Matrix,'Real','StorSymm')
!            IJ = 0 
!            Do J= 1,N2
!              Do I=1,J
!                IJ =IJ+1
!!               Alpha_Matrix%MatR(I,J) = R2D(I,J)
!                temp = R2D(IJ,1)
!                call Alpha_Matrix%put(temp,I,J)
!              EndDo
!            EndDo
!         EndIf
!      EndIf
!
!      Return
!      end Subroutine MQC_Copy_Array_To_Matrix
!
!      Subroutine MQC_Alloc_Int_Ptr( IntArray, Size, Data_Name, IOut)
!        Integer(kind=Int64),contiguous,Pointer,intent(Out)::IntArray(:)
!        Integer(kind=Int64),intent(In)::Size
!        character(len=*),intent(In) ::Data_Name
!        integer(kind=Int64),intent(in)::iOut
!        Integer(kind=Int32)::AllocateStatus
!!
!!     Format statements.
!!
! 1010   format( "Allocation of Size ", I12 )
!!
!!     Do the work...
!!
!        Allocate(IntArray(Size), STAT = AllocateStatus)
!        IF (AllocateStatus /= 0) Then
!           write( IOut, 1010 ) Size
!           call MQC_Error_A( 'MQC_Alloc_Int_Ptr: Not Enough Memory to allocate data structure', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        Return
!      End Subroutine MQC_Alloc_Int_Ptr
!      Subroutine MQC_Alloc_Real_Ptr( DRealArray, Size, Data_Name, IOut)
!        Real(kind=Real64),contiguous,Pointer,intent(Out)::DRealArray(:)
!        Integer(kind=Int64),intent(In)::Size
!        character(len=*),intent(In) ::Data_Name
!        integer(kind=Int64),intent(in)::iOut
!        Integer(kind=Int32)::AllocateStatus
!!
!!     Format statements.
!!
! 1010   format( "Allocation of Size ", I12 )
!!
!!     Do the work...
!!
!        Allocate(DRealArray(Size), STAT = AllocateStatus)
!        IF (AllocateStatus /= 0) Then
!           write( IOut, 1010 ) Size
!           call MQC_Error_A( 'MQC_Alloc_Real_Ptr: Not Enough Memory to allocate data structure', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        Return
!      End Subroutine MQC_Alloc_Real_Ptr
!      Subroutine MQC_Alloc_Complex_Ptr( CompArray, Size, Data_Name, IOut)
!        Complex(kind=Real64),contiguous,Pointer,intent(Out)::CompArray(:)
!        Integer(kind=Int64),intent(In)::Size
!        character(len=*),intent(In) ::Data_Name
!        integer(kind=Int64),intent(in)::iOut
!        Integer(kind=Int32)::AllocateStatus
!!
!!     Format statements.
!!
! 1010   format( "Complex Allocation of ",a," Size:", I12 )
!!
!!     Do the work...
!!
!#ifdef DEBUG 
!           write( IOut, 1010 ) trim(Data_Name),Size
!#endif
!
!        Allocate(CompArray(Size), STAT = AllocateStatus)
!        IF (AllocateStatus /= 0) Then
!           write( IOut, 1010 ) Size
!           call MQC_Error_A( 'MQC_Alloc_Complex_Ptr: Not Enough Memory to allocate data structure', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        Return
!      End Subroutine MQC_Alloc_Complex_Ptr
!
!      Subroutine MQC_DeAlloc_Int_Ptr( IntArray, Data_Name, IOut)
!        Integer(kind=Int64),contiguous,Pointer,intent(Out)::IntArray(:)
!        character(len=*),intent(In) ::Data_Name
!        integer(kind=Int64),intent(in)::iOut
!        Integer(kind=Int32)::AllocateStatus
!!
!!     Format statements.
!!
! 1900   Format( "MQC_DeAlloc_Int_Ptr: Deallocate ",A )
!!
!!     Do the work...
!!
!#ifdef DEBUG 
!        write( IOut, 1900) Data_Name
!        flush(IOut)
!#endif
!        If ( .not. Associated(IntArray) ) Then
!           call MQC_Error_A( 'MQC_DeAlloc_Int_Ptr: Tried to deallocate an array that is not there', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        DeAllocate(IntArray, STAT = AllocateStatus)
!        IF (AllocateStatus /= 0) Then
!           call MQC_Error_A( 'MQC_DeAlloc_Int_Ptr: Failure to deallocate data structure', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        IntArray=>null()
!        Return
!      End Subroutine MQC_DeAlloc_Int_Ptr
!      Subroutine MQC_DeAlloc_Real_Ptr( DRealArray, Data_Name, IOut)
!        Real(kind=Real64),contiguous,Pointer,intent(Out)::DRealArray(:)
!        character(len=*),intent(In) ::Data_Name
!        integer(kind=Int64),intent(in)::iOut
!        Integer(kind=Int32)::AllocateStatus
!!
!!     Format statements.
!!
! 1900   Format( "MQC_DeAlloc_Real_Ptr: Deallocate ",A )
!!
!!     Do the work...
!!
!#ifdef DEBUG 
!        write( IOut, 1900) Data_Name
!        flush(IOut)
!#endif
!        If ( .not. Associated(DRealArray) ) Then
!           call MQC_Error_A( 'MQC_DeAlloc_Real_Ptr: Tried to deallocate an array that is not there', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        DeAllocate(DRealArray, STAT = AllocateStatus)
!        IF (AllocateStatus /= 0) Then
!           call MQC_Error_A( 'MQC_DeAlloc_Real_Ptr: Failed to deallocate data structure', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        DRealArray=>null()
!        Return
!      End Subroutine MQC_DeAlloc_Real_Ptr
!      Subroutine MQC_DeAlloc_Complex_Ptr( CompArray, Data_Name, IOut)
!        Complex(kind=Real64),contiguous,Pointer,intent(Out)::CompArray(:)
!        Character(len=*),intent(In) ::Data_Name
!        Integer(kind=Int64),intent(in)::iOut
!        Integer(kind=Int32)::AllocateStatus
!!
!!     Format statements.
!!
! 1900   Format( "MQC_DeAlloc_Complex_Ptr: Deallocate ",A )
!!
!!     Do the work...
!!
!#ifdef DEBUG 
!        write( IOut, 1900) Data_Name
!        flush(IOut)
!#endif
!        If ( .not. Associated(CompArray) ) Then
!           call MQC_Error_A( 'MQC_DeAlloc_Complex_Ptr: Tried to deallocate an array that is not there', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        DeAllocate(CompArray, STAT = AllocateStatus)
!        IF (AllocateStatus /= 0) Then
!           call MQC_Error_A( 'MQC_DeAlloc_Complex_Ptr: Failure to deallocate data structure', &
!                IOut, 'DataName', Data_Name )
!        EndIf
!        CompArray=>null()
!        Return
!      End Subroutine MQC_DeAlloc_Complex_Ptr
    End module MQC_FullWavefunction
