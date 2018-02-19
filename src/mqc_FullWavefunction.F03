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

  ! Linked list array data type
  type :: mqc_array_link_list_node
     character(Len=64) :: Data_Name
     integer(kind=int64) :: NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR
     character(Len=64):: Data_Type,Shape
     character(Len=1024):: Formula
     Real(Kind=real64),Dimension(:,:,:,:,:),Allocatable::ArrR
     Integer(kind=int64),Dimension(:,:,:,:,:),Allocatable::ArrI
     Complex(Kind=real64),Dimension(:,:,:,:,:),Allocatable::ArrC
     type(mqc_array_link_list_node),pointer :: next =>null()
  end type mqc_array_link_list_node

  type :: mqc_array_link_list
     type(mqc_array_link_list_node),pointer :: first => null()
     type(mqc_array_link_list_node),pointer :: last => null()
     type(mqc_array_link_list_node),pointer :: current => null()
     integer(kind=int64),private :: length = 0
   contains
     procedure :: start => mqc_array_link_list_start
     procedure :: next => mqc_array_link_list_next
  end type mqc_array_link_list

  Type MQC_ALL_Wavefunction
     Character(Len=256)::FileName,Basis,Symmetry,WF_Type
     Character(Len=64)::LabFil
     Character(Len=64)::GVers
     Character(Len=64)::Title

     Integer(kind=int64)::UnitNumber,IVers,NLab,NAtoms,nbasisUse,Len12L,Len4L
     Integer(kind=int64)::IOpCl,ICGU,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB
     Integer(kind=int64)::NPrmDB,NBTot
     Integer(kind=int64)::NAlpha,NBeta,NElectrons,NBasis,iCharge,Multiplicity
     Integer(kind=int64)::NCore,NVal,NActive

     Logical::WF_Complex

     Type(mqc_array_link_list)::Array_Data

  End Type MQC_ALL_Wavefunction

   ! Linked list node data type
  type mqc_link_list_node
     type(MQC_ALL_Wavefunction) :: current_wavefunction
     type(mqc_link_list_node),pointer :: next =>null()
  end type mqc_link_list_node

  type :: mqc_link_list
     type(mqc_link_list_node),pointer :: first => null()
     type(mqc_link_list_node),pointer :: last => null()
     type(mqc_link_list_node),pointer :: current => null()
     integer(kind=int64),private :: length = 0
   contains
     procedure :: start => mqc_link_list_start
     procedure :: add => mqc_link_list_add
     procedure :: next => mqc_link_list_next
  end type mqc_link_list

  ! Linked list READONLY data type
  type :: mqc_READONLY_list_node
     character(Len=64) :: Data_Name
     type(mqc_READONLY_list_node),pointer :: next =>null()
  end type mqc_READONLY_list_node

  type :: mqc_READONLY_list
     type(mqc_READONLY_list_node),pointer :: first => null()
     type(mqc_READONLY_list_node),pointer :: last => null()
     type(mqc_READONLY_list_node),pointer :: current => null()
     integer(kind=int64),private :: length = 0
   contains
     procedure :: start => mqc_READONLY_list_start
     procedure :: next => mqc_READONLY_list_next
  end type mqc_READONLY_list

CONTAINS
      subroutine mqc_read_all_wavefunction( current_wavefunction, &
           my_READONLY_list, Print_Options, iout)
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
      integer(kind=int64),intent(in)::Print_Options
      integer(kind=int64),intent(in)::iOut
      Integer(kind=int64):: NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,NRI,LR
      logical ASym
      Character(Len=64)::CBuf
      Character(Len=2048)::DataType,Shape,Formula
      logical EOF
      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
      integer(kind=int32)::i
      integer(kind=int64)::IERROR, natoms, NBasis
!
!     Format statements.
!
 1000 format( "Warning: MQC Pack cannot process file '",a,"' it may have been created with the raw attribute.")
 1005 format( "MatrixFile:  ",a)
 1010 format('DataName  ',a,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
 1050 format( 'ERROR:',a,': please implement proper handling of this data type')
 1060 format( "The data is identifyed as: ",a)
 1070 format('NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
 2010 Format( i15 )
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

      If ( current_wavefunction%UnitNumber.eq.-2 ) then
         write( IOut, 1000 ) trim(current_wavefunction%FileName)
!
! Flag that this wavefunction is empty and should be ignored later
!
         current_wavefunction%NAtoms = -1
         return
      EndIF
      If ( Print_Options .eq. 1 .or. Print_Options .eq. 2 ) then
         write( IOut, 1005 ) trim(current_wavefunction%FileName)
      EndIF
!
! Read the next 9 records
!
      natoms = current_wavefunction%natoms
      NBasis = current_wavefunction%NBasis
      call MQC_Read_MatFile_9_Records( current_wavefunction, &
           my_READONLY_list, Print_Options, Natoms, NBASIS, IOut)

      EOF = .false.
      IERROR = 0
      do while(.not.EOF)
         call Rd_Labl(current_wavefunction%unitNumber, &
              current_wavefunction%IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2, &
              N3,N4,N5,ASym,NRI,EOF)
         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,DataType,Shape,Formula)
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
                 NI,NR,N1,N2,N3,N4,N5,NRI,NTot,LR, &
                 DataType,Shape,Formula)
               Call Rd_Skip(current_wavefunction%UnitNumber,NTot,LenBuf)
            endif
         else 
            if(trim(cBuf).eq.'END') then
               Call Close_MatF(current_wavefunction%unitNumber)
               if ( Print_Options == 2 ) then
                  call mqc_print_all_wavefunction(current_wavefunction, IOut)
               endif
               if (IERROR .ne. 0 ) then
                  call mqc_error_I('mqc_read_all_wavefunction: Found errors reading Matrix Element File', IOut)
               endif
               return
            endif

            call MQC_Unformatted_Read_Array(current_wavefunction, &
                 my_READONLY_list, cBuf, &
                 NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf,DataType,Shape,&
                 Formula,IOut,IERROR)
!
         endif
      enddo
      return
      end subroutine mqc_read_all_wavefunction
!
!     PROCEDURE MQC_Read_MatFile_9_Records
      Subroutine MQC_Read_MatFile_9_Records( wavefunction, &
           my_READONLY_list, Print_Options, Natoms, NBASIS, IOut)
      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
      integer(kind=int64),intent(in)::Print_Options
      integer(kind=int64),intent(in)::NAtoms,NBasis
      integer(kind=int64),intent(in)::iOut
      Integer(kind=int64)::atomicNumbers(natoms)
      Integer(kind=int64)::atomTypes(natoms)
      Integer(kind=int64)::basisFunction2Atom(NBasis)
      INTEGER(kind=int64)::IBasisFunctionType(NBasis)
      INTEGER(kind=int64)::NR,NI,NRI,N1,N2,N3,N4,N5,LR

      Real(kind=real64)::atomicCharges(natoms)
      Real(kind=real64)::atomicWeights(natoms)
      Real(kind=real64)::cartesians(3,natoms,1,1,1)
      Character(Len=2048)::Data_Type,Shape,Formula
      Integer(kind=int32)::I,J,K
      Integer(kind=int32)::AllocateStatus
!
!     Format statements.
!
 2020 Format ( "********************************************************", &
           /,10x,A,/, &
           "********************************************************",/, &
           "Gaussian Version:",A,/,&
           "Job Title:",A)

      call Rd_Head_G(wavefunction%unitNumber, wavefunction%NLab, &
           natoms, nbasis, atomicNumbers, atomTypes, atomicCharges, &
           cartesians, basisFunction2Atom, IBasisFunctionType, &
           atomicWeights, wavefunction%NFC, wavefunction%NFV, &
           wavefunction%ITran, wavefunction%IDum9, wavefunction%NShlAO, &
           wavefunction%NPrmAO, wavefunction%NShlDB, wavefunction%NPrmDB, &
           wavefunction%NBTot)
!
      if ( Print_Options == 1 ) then
!        Write a Summary of the data.  No need to save it
         write( IOut, 2020) trim(wavefunction%FileName), &
              trim(wavefunction%GVers), trim(wavefunction%Title)
         call MQC_PRINT_SCALAR_SUMMARY( iOut, wavefunction )
         call MQC_PRINT_DEFAULT_VECTOR_SUMMARY( iOut, wavefunction)
      else
         NRI=1
         N1=NAtoms
         N2=1
         N3=1
         N4=1
         N5=1
         NR=0
         NI=1
         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
         call mqc_array_link_list_add(wavefunction, "ATOMIC NUMBERS", &
              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
         allocate(wavefunction%Array_Data%last%ArrI(natoms,1,1,1,1), &
              STAT = AllocateStatus)
         IF (AllocateStatus /= 0) call MQC_Error_A( & 
              'Not Enough Memory to allocate data structure', Iout, &
              'DataName', 'ATOMIC NUMBERS')
         wavefunction%Array_Data%last%ArrI(1:NAtoms,1,1,1,1) = &
              atomicNumbers(1:NAtoms)

         call mqc_array_link_list_add(wavefunction, "ATOM TYPES", &
              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
         allocate(wavefunction%Array_Data%last%ArrI(NAtoms,1,1,1,1), &
              STAT = AllocateStatus)
         IF (AllocateStatus /= 0) call MQC_Error_A( & 
              'Not Enough Memory to allocate data structure', Iout, &
              'DataName', 'ATOM TYPES')
         wavefunction%Array_Data%last%ArrI(1:NAtoms,1,1,1,1) = &
              atomTypes(1:NAtoms)

         NR=1
         NI=0
         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
         call mqc_array_link_list_add(wavefunction,"ATOM CHARGES", &
              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
         allocate(wavefunction%Array_Data%last%ArrR(natoms,1,1,1,1), &
              STAT = AllocateStatus)
         IF (AllocateStatus /= 0) call MQC_Error_A( & 
              'Not Enough Memory to allocate data structure', Iout, &
              'DataName', 'ATOM CHARGES')
         wavefunction%Array_Data%last%ArrR(1:NAtoms,1,1,1,1) = &
              atomicCharges(1:NAtoms)
         N1=3
         N2=NAtoms
         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
         call mqc_array_link_list_add(wavefunction, &
              'CARTESIAN NUCLEAR COORDINATES(BOHR)', &
              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
         allocate(wavefunction%Array_Data%last%ArrR(3,natoms,1,1,1), &
              STAT = AllocateStatus)
         IF (AllocateStatus /= 0) call MQC_Error_A( & 
              'Not Enough Memory to allocate data structure', Iout, &
              'DataName', 'CARTESIAN NUCLEAR COORDINATES(BOHR)')
         wavefunction%Array_Data%last%ArrR(1:3,1:NAtoms,1,1,1) = &
              cartesians(1:3,1:NAtoms,1,1,1)
         N1=NBasis
         N2=1
         NR=0
         NI=1
         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
         call mqc_array_link_list_add(wavefunction, &
              "MAP FROM BASIS FUNCTIONS TO ATOMS", &
              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
         allocate(wavefunction%Array_Data%last%ArrI(NBasis,1,1,1,1), &
              STAT = AllocateStatus)
         IF (AllocateStatus /= 0) call MQC_Error_A( & 
              'Not Enough Memory to allocate data structure', Iout, &
              'DataName', 'MAP FROM BASIS FUNCTIONS TO ATOMS')
         wavefunction%Array_Data%last%ArrI(1:NBasis,1,1,1,1) = &
              basisFunction2Atom(1:NBasis)

         call mqc_array_link_list_add(wavefunction, &
              "TYPE FLAG FOR EACH BASIS FUNCTION", &
              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
         allocate(wavefunction%Array_Data%last%ArrI(NBasis,1,1,1,1), &
              STAT = AllocateStatus)
         IF (AllocateStatus /= 0) call MQC_Error_A( & 
              'Not Enough Memory to allocate data structure', Iout, &
              'DataName', 'TYPE FLAG FOR EACH BASIS FUNCTION')
         wavefunction%Array_Data%last%ArrI(1:NBasis,1,1,1,1) = &
              IBasisFunctionType(1:NBasis)

         N1=NAtoms
         NR=1
         NI=0
         LR = mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5,Data_Type,Shape,Formula)
         call mqc_array_link_list_add(wavefunction, "ATOM WEIGHTS", &
              NR,NI,NRI,LR,N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
         allocate(wavefunction%Array_Data%last%ArrR(NAtoms,1,1,1,1), &
              STAT = AllocateStatus)
         IF (AllocateStatus /= 0) call MQC_Error_A( & 
              'Not Enough Memory to allocate data structure', Iout, &
              'DataName', 'ATOM WEIGHTS')
         wavefunction%Array_Data%last%ArrR(1:NAtoms,1,1,1,1) = &
              atomicWeights(1:NAtoms)
      endif
      return
      end Subroutine MQC_Read_MatFile_9_Records
!
!     PROCEDURE MQC_Print_all_Wavefunction     
      subroutine mqc_print_all_wavefunction(wavefunction,iOut)
!
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
      integer(kind=int64),intent(in)::iOut
      integer(kind=int64)::LR

 1050 Format( 2A )
 1060 Format( A,L10 )

      call MQC_PRINT_SCALAR_HEAD( iOut, wavefunction )
!                                                                              
      wavefunction%Array_Data%last => wavefunction%Array_Data%first
      Do While( associated(wavefunction%Array_Data%last) )
         if (trim(wavefunction%Array_Data%last%Data_Type).eq."INTEGER" ) then
            call MQC_Print_Array_Integer(IOut, wavefunction%Array_Data%last)
         else if (trim(wavefunction%Array_Data%last%Data_Type).eq."REAL" ) then
            call MQC_Print_Array_REAL(IOut, wavefunction%Array_Data%last)
         else if (trim(wavefunction%Array_Data%last%Data_Type).eq."COMPLEX" ) then
            call MQC_Print_Array_COMPLEX(IOut, wavefunction%Array_Data%last)
         else if (trim(wavefunction%Array_Data%last%Data_Type).eq."MIXED" ) then
! For some reason, LR and N1 are improperly much larger than NTot.  
! For write purposes, set LR to NTot

            LR = wavefunction%Array_Data%last%LR
            wavefunction%Array_Data%last%LR = min( LR, &
                 wavefunction%Array_Data%last%NTot)
            call MQC_Print_Array_Integer(IOut, &
              wavefunction%Array_Data%last)
            if (wavefunction%Array_Data%last%NR.lt.0 ) then
               call MQC_Print_Array_COMPLEX(IOut, &
                    wavefunction%Array_Data%last)
            else
               call MQC_Print_Array_REAL(IOut, &
                    wavefunction%Array_Data%last)
            endif
            wavefunction%Array_Data%last%LR=LR
         else if (trim(wavefunction%Array_Data%last%Data_Type).eq."TENSOR" ) then
            if (wavefunction%Array_Data%last%NR.lt.0 ) then
               call MQC_Print_Array_COMPLEX(IOut, &
                    wavefunction%Array_Data%last)
            else
               call MQC_Print_Array_REAL(IOut, &
                    wavefunction%Array_Data%last)
            endif
         else
            call MQC_PRINT_DATA_TYPE(IOut, &
              wavefunction%Array_Data%last%Data_Name, &
              wavefunction%Array_Data%last%NI, &
              wavefunction%Array_Data%last%NR, &
              wavefunction%Array_Data%last%N1, &
              wavefunction%Array_Data%last%N2, &
              wavefunction%Array_Data%last%N3, &
              wavefunction%Array_Data%last%N4, &
              wavefunction%Array_Data%last%N5, &
              wavefunction%Array_Data%last%NRI, &
              wavefunction%Array_Data%last%NTot, &
              wavefunction%Array_Data%last%LR, &
              wavefunction%Array_Data%last%Data_Type, &
              wavefunction%Array_Data%last%Shape, &
              wavefunction%Array_Data%last%Formula)
         endif
         wavefunction%Array_Data%last => wavefunction%Array_Data%last%next
      enddo
!
      write( iout, 1050) "Finished mqc_print_all_wavefunction"
      end subroutine mqc_print_all_wavefunction

!
!     PROCEDURE MQC_Find_mqc_array_link_list_node
      Integer(kind=int64) Function MQC_Find_mqc_array_link_list_node(wavefunction,Name,Found, &
           iOut)
!
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
      character(len=*),intent(in) :: Name
      integer(kind=int64),intent(in)::iOut
      type(mqc_array_link_list_node),pointer,intent(out) ::Found
      integer(kind=int64)::LR
 1050 Format( 2A )
 1060 Format( A,L10 )

      wavefunction%Array_Data%last => wavefunction%Array_Data%first
      Do While( associated(wavefunction%Array_Data%last) )
         if (trim(wavefunction%Array_Data%last%Data_Name) .eq. trim(Name) ) then
            Found => wavefunction%Array_Data%last
            If ( allocated(Found%ArrI) ) then
               If ( allocated(Found%ArrR) ) then
                  MQC_Find_mqc_array_link_list_node = 4
               Else If ( allocated(Found%ArrC) ) then
                  MQC_Find_mqc_array_link_list_node = 5
               Else 
                  MQC_Find_mqc_array_link_list_node = 1
               EndIf
            Else If ( allocated(Found%ArrR) ) then
               MQC_Find_mqc_array_link_list_node = 2
            Else If ( allocated(Found%ArrC) ) then
               MQC_Find_mqc_array_link_list_node = 3
            EndIf
            return
         endif
         wavefunction%Array_Data%last => wavefunction%Array_Data%last%next
      enddo
!
      Found => null()
      MQC_Find_mqc_array_link_list_node = 0

      end Function MQC_Find_mqc_array_link_list_node

      Subroutine MQC_Print_Array_Integer(IOut, Array,Blank_At_Top, &
        Blank_At_Bottom)
      implicit none
      class(mqc_array_link_list_node),intent(in):: Array
      integer(kind=int64),intent(in)::iOut
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64),Parameter::ColWidth=10 
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64)::N1,N2,N3,N4,K,L,IBase,JBase,LMAX
!
!     Format statements.
!
 1000 Format(/,1x,A)
 1001 Format(1x,I7,2x,I14)
 1010 Format(1x,A,2x,I14)
 1020 Format(" ")
 2001 Format(5x,10(7x,I7))
 2002 Format(5x,10(5x,I7,",",I1))
 2003 Format(5x,10(4x,I7,",",I2))
 2004 Format(5x,10(I7,",",I6))
 2005 Format(5x,10(3x,I7,",",I1,",",I1))
 2006 Format(5x,10(3x,I7,",",I7,",",I7))
 2011 Format(1x,I7,10I14)
 4010 format( "(",I5,",",I5,"|",I5,",",I5,") = ",I10 )
 4020 format( "(",I5,",",I5,",",I5,") = ",I10 )
 5010 Format("Warning: Do Not know how do deal with the shape of ",A)

      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      if ( Array%Shape.ne."SCALAR") then
         Write(IOut,1000) TRIM(Array%Data_Name)
      EndIf

      if ( Array%Shape.eq."SCALAR") then
         Write(IOut,1010) TRIM(Array%Data_Name), Array%ArrI(1,1,1,1,1)
      else if ( Array%Shape.eq."VECTOR") then
         Do I = 1, Array%LR
            Write(IOut,1001) I,  Array%ArrI(I,1,1,1,1)
         EndDo
      else if ( Array%Shape.eq."MATRIX") then
         NRows = Array%N1
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = 1,NRows
               Write(IOut,2011) I, (Array%ArrI(I,J,1,1,1),J=IFirst,ILast)
            EndDo
         EndDo
      else if ( Array%Shape.eq."SYMMATRIX") then
         NRows = Array%N2
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = IFirst, ILAST
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (Array%ArrI(IBase+J,1,1,1,1),J=IFirst,I)
            EndDo
            Do I = ILAST+1, NCols
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (Array%ArrI(IBase+J,1,1,1,1),J=IFirst,ILast)
            EndDo
         EndDo
      else if ( Array%Shape.eq."MULTISYMMATRICES") then
         NRows = Array%N2
         NCols = Array%N2
         N3 = Array%N3
         N4 = Array%N4
         Do L = 1,N4
            Do K = 1,N3
               IBase = 0
               Do IFirst = 1,NCols,ColWidth
                  ILast = Min(IFirst+ColWidth-1,NCols)
                  if ( N4 .eq. 1 ) then
                     if ( K .lt. 10 ) then
                        Write(IOut,2002) (I,K,I=IFirst,ILast)
                     else if ( K .lt. 100 ) then
                        Write(IOut,2003) (I,K,I=IFirst,ILast)
                     else
                        Write(IOut,2004) (I,K,I=IFirst,ILast)
                     endif
                  else 
                     If ( K .lt. 10 .and. L .lt. 10 ) then
                        Write(IOut,2005) (I,K,L,I=IFirst,ILast)
                     Else
                        Write(IOut,2006) (I,K,L,I=IFirst,ILast)
                     EndIf
                  endif
                  Do I = IFirst, ILAST
                     JBase = IBase+(I*(I-1))/2
                     Write(IOut,2011) I, (Array%ArrI(JBase+J,K,L,1,1),J=IFirst,I)
                  EndDo
                  Do I = ILAST+1, NCols
                     Write(IOut,2011) I, (Array%ArrI(JBase+J,K,L,1,1),J=IFirst,ILast)
                  EndDo
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."SYMSYM") then
         N4 = Array%N4
         IBase=0
         Do I = 1,N4
            Do J = 1,I
               Do K = 1,I
                  if ( I .eq. K ) then
                     LMax = J
                  else
                     LMax = K
                  endif
                  Do L = 1,LMax
                     IBase=IBase+1
                     write(Iout,4010) I,J,K,L,Array%ArrI(IBase,1,1,1,1)
                  EndDo
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."TRISYMMATRICES") then
         N3 = Array%N3
         IBase=0
         Do I = 1,N3
            Do J = 1,I
               Do K = 1,J
                  IBase=IBase+1
                  write(Iout,4020) I,J,K,Array%ArrI(IBase,1,1,1,1)
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."3DMATRIX") then
         N1 = Array%N1
         N2 = Array%N2
         N3 = Array%N3
         IBase=0
         Do K = 1,N3
            Do J = 1,N2
               Do I = 1,N1
                  IBase=IBase+1
                  write(Iout,4020) I,J,K,Array%ArrI(I,J,K,1,1)
               EndDo
            EndDo
         EndDo
      else
         write(IOut,5010) trim(Array%Data_Name)
         call MQC_PRINT_DATA_TYPE(IOut, Array%Data_Name, Array%NI, Array%NR, &
              Array%N1, Array%N2, Array%N3, Array%N4, Array%N5, Array%NRI, &
              Array%NTot, Array%LR, Array%Data_Type, Array%Shape, &
              Array%Formula)
      endif

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
      return
      end Subroutine MQC_Print_Array_Integer

      Subroutine MQC_Print_Array_REAL(IOut, Array,Blank_At_Top, &
        Blank_At_Bottom)
      implicit none
      class(mqc_array_link_list_node),intent(in):: Array
      integer(kind=int64),intent(in)::iOut
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64),Parameter::ColWidth=10 
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64)::N1,N2,N3,N4,K,L,IBase,JBase,LMAX
!
!     Format statements.
!
 1000 Format(/,1x,A)
 1001 Format(1x,I7,2x,F14.6)
 1010 Format(1x,A,2x,F14.6)
 1020 Format(" ")
 2001 Format(5x,10(7x,I7))
 2002 Format(5x,10(5x,I7,",",I1))
 2003 Format(5x,10(4x,I7,",",I2))
 2004 Format(5x,10(I7,",",I6))
 2005 Format(5x,10(3x,I7,",",I1,",",I1))
 2006 Format(5x,10(3x,I7,",",I7,",",I7))
 2011 Format(1x,I7,10F14.6)
 4010 format( "(",I5,",",I5,"|",I5,",",I5,") = ",F14.6 )
 4020 format( "(",I5,",",I5,",",I5,") = ",F14.6 )
 5010 Format("Warning: Do Not know how do deal with the shape of ",A)

      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      if ( Array%Shape.ne."SCALAR") then
         Write(IOut,1000) TRIM(Array%Data_Name)
      EndIf

      if ( Array%Shape.eq."SCALAR") then
         Write(IOut,1010) TRIM(Array%Data_Name), Array%ArrR(1,1,1,1,1)
      else if ( Array%Shape.eq."VECTOR") then
         Do I = 1, Array%LR
            Write(IOut,1001) I,  Array%ArrR(I,1,1,1,1)
         EndDo
      else if ( Array%Shape.eq."MATRIX") then
         NRows = Array%N1
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = 1,NRows
               Write(IOut,2011) I, (Array%ArrR(I,J,1,1,1),J=IFirst,ILast)
            EndDo
         EndDo
      else if ( Array%Shape.eq."SYMMATRIX") then
         NRows = Array%N2
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = IFirst, ILAST
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (Array%ArrR(IBase+J,1,1,1,1),J=IFirst,I)
            EndDo
            Do I = ILAST+1, NCols
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (Array%ArrR(IBase+J,1,1,1,1),J=IFirst,ILast)
            EndDo
         EndDo
      else if ( Array%Shape.eq."MULTISYMMATRICES") then
         NRows = Array%N2
         NCols = Array%N2
         N3 = Array%N3
         N4 = Array%N4
         Do L = 1,N4
            Do K = 1,N3
               IBase = 0
               Do IFirst = 1,NCols,ColWidth
                  ILast = Min(IFirst+ColWidth-1,NCols)
                  if ( N4 .eq. 1 ) then
                     if ( K .lt. 10 ) then
                        Write(IOut,2002) (I,K,I=IFirst,ILast)
                     else if ( K .lt. 100 ) then
                        Write(IOut,2003) (I,K,I=IFirst,ILast)
                     else
                        Write(IOut,2004) (I,K,I=IFirst,ILast)
                     endif
                  else 
                     If ( K .lt. 10 .and. L .lt. 10 ) then
                        Write(IOut,2005) (I,K,L,I=IFirst,ILast)
                     Else
                        Write(IOut,2006) (I,K,L,I=IFirst,ILast)
                     EndIf
                  endif
                  Do I = IFirst, ILAST
                     JBase = IBase+(I*(I-1))/2
                     Write(IOut,2011) I, (Array%ArrR(JBase+J,K,L,1,1),J=IFirst,I)
                  EndDo
                  Do I = ILAST+1, NCols
                     Write(IOut,2011) I, (Array%ArrR(JBase+J,K,L,1,1),J=IFirst,ILast)
                  EndDo
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."SYMSYM") then
         N4 = Array%N4
         IBase=0
         Do I = 1,N4
            Do J = 1,I
               Do K = 1,I
                  if ( I .eq. K ) then
                     LMax = J
                  else
                     LMax = K
                  endif
                  Do L = 1,LMax
                     IBase=IBase+1
                     write(Iout,4010) I,J,K,L,Array%ArrR(IBase,1,1,1,1)
                  EndDo
               EndDo
            EndDo
         EndDo
         If ( IBASE .ne. Array%LR ) then
            call mqc_error_I('MQC_Print_Array_REAL: Calculated wrong size of Super Symmertric array.', IOut, &
            'IBASE', IBASE, 'Array%LR', Array%LR )
         EndIf
      else if ( Array%Shape.eq."TRISYMMATRICES") then
         N3 = Array%N3
         IBase=0
         Do I = 1,N3
            Do J = 1,I
               Do K = 1,J
                  IBase=IBase+1
                  write(Iout,4020) I,J,K,Array%ArrR(IBase,1,1,1,1)
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."3DMATRIX") then
         N1 = Array%N1
         N2 = Array%N2
         N3 = Array%N3
         IBase=0
         Do K = 1,N3
            Do J = 1,N2
               Do I = 1,N1
                  IBase=IBase+1
                  write(Iout,4020) I,J,K,Array%ArrR(I,J,K,1,1)
               EndDo
            EndDo
         EndDo
      else
         write(IOut,5010) trim(Array%Data_Name)
         call MQC_PRINT_DATA_TYPE(IOut, Array%Data_Name, Array%NI, Array%NR, &
              Array%N1, Array%N2, Array%N3, Array%N4, Array%N5, Array%NRI, &
              Array%NTot, Array%LR, Array%Data_Type, Array%Shape, &
              Array%Formula)
      EndIf

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
      return
      end Subroutine MQC_Print_Array_REAL

      Subroutine MQC_Print_Array_COMPLEX(IOut, Array,Blank_At_Top, &
        Blank_At_Bottom)
      implicit none
      class(mqc_array_link_list_node),intent(in):: Array
      integer(kind=int64),intent(in)::iOut
      Logical,Intent(In),Optional::Blank_At_Top,Blank_At_Bottom
      Integer(kind=int64),Parameter::ColWidth=5 
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
      Integer(kind=int64)::N1,N2,N3,N4,K,L,IBase,JBase,LMAX
!
!     Format statements.
!
 1000 Format(/,1x,A)
 1001 Format(1x,I7,2x,"(",F14.6,",",F14.6,"i)")
 1010 Format(1x,A,2x,"(",F14.6,",",F14.6,"i)")
 1020 Format(" ")
 2001 Format(5x,5(20x,I7))
 2002 Format(5x,5(18x,I7,",",I1))
 2003 Format(5x,5(17x,I7,",",I2))
 2004 Format(5x,5(13x,I7,",",I6))
 2005 Format(5x,5(16x,I7,",",I1,",",I1))
 2006 Format(5x,5(4x,I7,",",I7,",",I7))
 2011 Format(1x,I7," (",F13.6,",",F13.6,"i)  (",F13.6,",",F13.6,"i)  (", &
           F13.6,",",F13.6,"i)  (",F13.6,",",F13.6,"i)  (",F13.6,",",F13.6,"i)")
 4010 format( "(",I5,",",I5,"|",I5,",",I5,") = (",F14.6,",",F14.6,"i)")
 4020 format( "(",I5,",",I5,",",I5,") = (",F14.6,",",F14.6,"i)")
 5010 Format("Warning: Do Not know how do deal with the shape of ",A)

      If(PRESENT(Blank_At_Top)) then
        If(Blank_At_Top) Write(IOut,1020)
      EndIf
      if ( Array%Shape.ne."SCALAR") then
         Write(IOut,1000) TRIM(Array%Data_Name)
      EndIf

      if ( Array%Shape.eq."SCALAR") then
         Write(IOut,1010) TRIM(Array%Data_Name), REAL(Array%ArrC(1,1,1,1,1)), &
              AIMAG(Array%ArrC(1,1,1,1,1))
      else if ( Array%Shape.eq."VECTOR") then
         Do I = 1, Array%LR
            Write(IOut,1001) I, REAL(Array%ArrC(I,1,1,1,1)), &
                 AIMAG(Array%ArrC(I,1,1,1,1))
         EndDo
      else if ( Array%Shape.eq."MATRIX") then
         NRows = Array%N1
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = 1,NRows
               Write(IOut,2011) I, (REAL(Array%ArrC(I,J,1,1,1)), &
                    AIMAG(Array%ArrC(I,J,1,1,1)),J=IFirst,ILast)
            EndDo
         EndDo
      else if ( Array%Shape.eq."SYMMATRIX") then
         NRows = Array%N2
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = IFirst, ILAST
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (REAL(Array%ArrC(IBase+J,1,1,1,1)), &
                    AIMAG(Array%ArrC(IBase+J,1,1,1,1)),J=IFirst,I)
            EndDo
            Do I = ILAST+1, NCols
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (REAL(Array%ArrC(IBase+J,1,1,1,1)), &
                    AIMAG(Array%ArrC(IBase+J,1,1,1,1)),J=IFirst,ILAST)
            EndDo
         EndDo
      else if ( Array%Shape.eq."MULTISYMMATRICES") then
         NRows = Array%N2
         NCols = Array%N2
         N3 = Array%N3
         N4 = Array%N4
         Do L = 1,N4
            Do K = 1,N3
               IBase = 0
               Do IFirst = 1,NCols,ColWidth
                  ILast = Min(IFirst+ColWidth-1,NCols)
                  if ( N4 .eq. 1 ) then
                     if ( K .lt. 10 ) then
                        Write(IOut,2002) (I,K,I=IFirst,ILast)
                     else if ( K .lt. 100 ) then
                        Write(IOut,2003) (I,K,I=IFirst,ILast)
                     else
                        Write(IOut,2004) (I,K,I=IFirst,ILast)
                     endif
                  else 
                     If ( K .lt. 10 .and. L .lt. 10 ) then
                        Write(IOut,2005) (I,K,L,I=IFirst,ILast)
                     Else
                        Write(IOut,2006) (I,K,L,I=IFirst,ILast)
                     EndIf
                  endif
                  Do I = IFirst, ILAST
                     JBase = IBase+(I*(I-1))/2
                     Write(IOut,2011) I, (REAL(Array%ArrC(JBase+J,K,L,1,1)), &
                          AIMAG(Array%ArrC(JBase+J,K,L,1,1)),J=IFirst,I)
                  EndDo
                  Do I = ILAST+1, NCols
                     Write(IOut,2011) I, (REAL(Array%ArrC(JBase+J,K,L,1,1)), &
                          AIMAG(Array%ArrC(JBase+J,K,L,1,1)),J=IFirst,ILAST)
                  EndDo
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."SYMSYM") then
         N4 = Array%N4
         IBase=0
         Do I = 1,N4
            Do J = 1,I
               Do K = 1,I
                  if ( I .eq. K ) then
                     LMax = J
                  else
                     LMax = K
                  endif
                  Do L = 1,LMax
                     IBase=IBase+1
                     write(Iout,4010) I,J,K,L,REAL(Array%ArrC(IBase,1,1,1,1)), &
                          AIMAG(Array%ArrC(IBase,1,1,1,1))
                  EndDo
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."TRISYMMATRICES") then
         N3 = Array%N3
         IBase=0
         Do I = 1,N3
            Do J = 1,I
               Do K = 1,J
                  IBase=IBase+1
                  write(Iout,4020) I,J,K,REAL(Array%ArrC(IBase,1,1,1,1)), &
                       AIMAG(Array%ArrC(IBase,1,1,1,1))
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."3DMATRIX") then
         N1 = Array%N1
         N2 = Array%N2
         N3 = Array%N3
         IBase=0
         Do K = 1,N3
            Do J = 1,N2
               Do I = 1,N1
                  write(Iout,4020) I,J,K,REAL(Array%ArrC(I,J,K,1,1)), &
                       AIMAG(Array%ArrC(I,J,K,1,1))
               EndDo
            EndDo
         EndDo
      else
         write(IOut,5010) trim(Array%Data_Name)
         call MQC_PRINT_DATA_TYPE(IOut, Array%Data_Name, Array%NI, Array%NR, &
              Array%N1, Array%N2, Array%N3, Array%N4, Array%N5, Array%NRI, &
              Array%NTot, Array%LR, Array%Data_Type, Array%Shape, &
              Array%Formula)
      endif

      If(PRESENT(Blank_At_Bottom)) then
        If(Blank_At_Bottom) Write(IOut,1020)
      EndIf
      return
      end Subroutine MQC_Print_Array_COMPLEX

!=====================================================================
!
! Functions for mqc_array_link_list class structures
!
!=====================================================================
      subroutine mqc_array_link_list_add(wavefunction,DataName,NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,IOut, &
           Data_Type,Shape,Formula)
        class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
        character(len=*),intent(in) :: DataName
        integer(kind=int64),intent(in) :: NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,IOut
        character(Len=*),intent(in):: Data_Type, Shape, Formula
        Integer(kind=int32)::AllocateStatus

        if (.not.associated(wavefunction%Array_Data%first)) then
           allocate(wavefunction%Array_Data%first, STAT = AllocateStatus)
           IF (AllocateStatus /= 0) call MQC_Error_A('Not Enough Memory to allocate data structure', Iout, 'DataName', DataName)
           wavefunction%Array_Data%last => wavefunction%Array_Data%first
        else
           allocate(wavefunction%Array_Data%last%next, STAT = AllocateStatus)
           IF (AllocateStatus /= 0) call MQC_Error_A('Not Enough Memory to allocate data structure', Iout, 'DataName', DataName)
           wavefunction%Array_Data%last => wavefunction%Array_Data%last%next
        endif
        wavefunction%Array_Data%last%Data_Name=trim(DataName)
        wavefunction%Array_Data%last%NR=NR
        wavefunction%Array_Data%last%NI=NI
        wavefunction%Array_Data%last%NRI=NRI
        wavefunction%Array_Data%last%NTot=NTot
        wavefunction%Array_Data%last%N1=N1
        wavefunction%Array_Data%last%N2=N2
        wavefunction%Array_Data%last%N3=N3
        wavefunction%Array_Data%last%N4=N4
        wavefunction%Array_Data%last%N5=N5
        wavefunction%Array_Data%last%LR=LR
        wavefunction%Array_Data%last%Data_Type=trim(Data_Type)
        wavefunction%Array_Data%last%Shape=trim(Shape)
        wavefunction%Array_Data%last%Formula=trim(Formula)
      end subroutine mqc_array_link_list_add

      subroutine mqc_array_link_list_start(self,wavefunction)
        class(MQC_ALL_Wavefunction),intent(in) :: wavefunction
        class(mqc_array_link_list),intent(inout) :: self

        self%current => wavefunction%Array_Data%first
      end subroutine mqc_array_link_list_start

      subroutine mqc_array_link_list_next(self)
        class(mqc_array_link_list),intent(inout) :: self

        if (associated(self%current)) then
           self%current => self%current%next
        else
           self%current => null()
        end if
      end subroutine mqc_array_link_list_next

!=====================================================================
!
! Functions for mqc_link_list class structures
!
!=====================================================================
      subroutine mqc_link_list_add(self,Print_Options, &
              my_READONLY_list,FileName,iout)
        class(mqc_link_list),intent(inout) :: self
        integer(kind=int64),intent(in)::Print_Options
        class(mqc_READONLY_list),intent(inout):: my_READONLY_list
        character(len=256),intent(in) :: FileName
        integer(kind=int64),intent(in)::iOut

        if (.not.associated(self%last)) then
           allocate(self%first)
           self%last => self%first
        else
           allocate(self%last%next)
           self%last => self%last%next
        endif

        self%last%current_wavefunction%FileName=trim(FileName)

        call mqc_read_all_wavefunction( self%last%current_wavefunction, &
             my_READONLY_list, Print_Options,iout)

        self%length = self%length + 1
      end subroutine mqc_link_list_add

      subroutine mqc_link_list_start(self,Status,wavefunction)
        class(mqc_link_list),intent(inout) :: self
        integer(kind=int64),intent(out)::Status
        class(MQC_ALL_Wavefunction),pointer,intent(out)::wavefunction

        self%current => self%first
        If (associated(self%current)) then
           wavefunction => self%current%current_wavefunction
           Status=1
           if ( wavefunction%NAtoms .eq. -1 ) then
              call mqc_link_list_next(self,Status,wavefunction)
           endif
        Else
           Status=0
        Endif
      end subroutine mqc_link_list_start

      RECURSIVE subroutine mqc_link_list_next(self,Status,wavefunction)
        class(mqc_link_list),intent(inout)::self
        integer(kind=int64),intent(out)::Status
        class(MQC_ALL_Wavefunction),pointer,intent(out)::wavefunction

        if (associated(self%current)) then
           self%current => self%current%next
           If (associated(self%current)) then
              wavefunction => self%current%current_wavefunction
              Status=1
              if ( wavefunction%NAtoms .eq. -1 ) then
!
! if the wavefunction is empty, skip it and go to the next one
!
                 call mqc_link_list_next(self,Status,wavefunction)
              EndIf
           Else
              Status=0
           EndIf
        else
           wavefunction => null()
           Status=0
        end if
      end subroutine mqc_link_list_next

!=====================================================================
!
! Functions for mqc_link_list class structures
!
!=====================================================================
      subroutine mqc_READONLY_list_add(self,Data_In)
        class(mqc_READONLY_list),intent(inout) :: self
        character(len=*),intent(in) :: Data_In

        if (.not.associated(self%last)) then
           allocate(self%first)
           self%last => self%first
        else
           allocate(self%last%next)
           self%last => self%last%next
        endif

        self%last%Data_Name=trim(Data_In)
      end subroutine mqc_READONLY_list_add

      subroutine mqc_READONLY_list_start(self,Status,Data_Name)
        class(mqc_READONLY_list),intent(inout) :: self
        integer(kind=int64),intent(out)::Status
        character(len=64),intent(out) :: Data_Name

        self%current => self%first
        If (associated(self%current)) then
           Data_Name = self%current%Data_Name
           Status=1
        Else
           Status=0
        Endif
      end subroutine mqc_READONLY_list_start

      subroutine mqc_READONLY_list_next(self,Status,Data_Name)
        class(mqc_READONLY_list),intent(inout)::self
        integer(kind=int64),intent(out)::Status
        character(len=64),intent(out) :: Data_Name

        if (associated(self%current)) then
           If (associated(self%current%next)) then
              self%current => self%current%next
              Data_Name = self%current%Data_Name
              Status=1
           Else
              Status=0
           EndIf
        else
           Status=0
        end if
      end subroutine mqc_READONLY_list_next

!=====================================================================
!
!PROCEDURE MQC_Unformatted_Read_Array
      subroutine MQC_Unformatted_Read_Array(wavefunction, &
           my_READONLY_list,Data_Name, &
           NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf,Data_Type,Shape,Formula, &
           IOut,IERROR)
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
      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
      character(len=*),intent(in)::Data_Name
      integer(kind=int64),intent(in)::NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf
      character(len=*),intent(in)::Data_Type,Shape,Formula
      integer(kind=int64),intent(in)::IOut
      integer(kind=int64),intent(inout)::IERROR
!
!     Temporary local variables used when calling the gauopen routines.
!     Local temp variables.
      character(len=256)::errorMsg
      Integer(kind=int32)::AllocateStatus
!
!     Format statements.
!
 1010 Format(A,/,' Data Name:',a,/,'  NI =',I9,/,'  NR =',I9,/,'  NRI =',I9,/,&
          '  NTot =',I9,/, '  LenBuf =',I9,/, '  N1 =',I9,/, '  N2 =',I9,/ &
          '  N3 =',I9,/, '  N4 =',I9,/, '  N5 =',I9,/, '  LR =',I9, &
          '  Calculated Type: ',A )
 1100 Format( "Error Allocating Memory for ",A)
!
!     Begin by seeing if a new file or filename has been sent by the calling
!     program unit. If so, then get the file declared before reading the
!     header information.
!
!     Look for the label sent by the calling program unit. If the label is
!     found, then load <matrixOut> with the data on the file.
!
      If ( .not.mqc_read_this_MatrixFile_Element(my_READONLY_list,Data_Name )) then
         Call Rd_Skip(wavefunction%UnitNumber,NTot,LenBuf)
         return
      EndIf

      call mqc_array_link_list_add(wavefunction,Data_Name,NR,NI,NRI,NTot, & 
           N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
!
      select case(Data_Type)
      case('INTEGER')
         If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
              Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
              Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
            allocate(wavefunction%Array_Data%last%ArrI(LR,1,1,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."MATRIX") then
            allocate(wavefunction%Array_Data%last%ArrI(N1,N2,1,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."MULTISYMMATRICES") then
            allocate(wavefunction%Array_Data%last%ArrI(((N2+1)*N2)/2,N3,N4,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."3DMATRIX") then
            allocate(wavefunction%Array_Data%last%ArrI(N1,N2,N3,1,1), STAT = AllocateStatus)
         Else
            call MQC_Error_A("Do not know how to allocate array", Iout, &
                 "Shape",Shape)
         Endif
         IF (AllocateStatus /= 0) then
            write( errorMsg, 1100 ) trim(Data_Name)
            call MQC_Error_I(errorMsg, Iout, 'LR', LR)
         EndIf
         call Rd_IBuf_G(wavefunction%unitNumber,NTot,LenBuf, &
              wavefunction%Array_Data%last%ArrI)
      case('REAL')
         If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
              Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
              Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
            allocate(wavefunction%Array_Data%last%ArrR(LR,1,1,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."MATRIX") then
            allocate(wavefunction%Array_Data%last%ArrR(N1,N2,1,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."MULTISYMMATRICES") then
            allocate(wavefunction%Array_Data%last%ArrR(((N2+1)*N2)/2,N3,N4,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."3DMATRIX") then
            allocate(wavefunction%Array_Data%last%ArrR(N1,N2,N3,1,1), STAT = AllocateStatus)
         Else
            call MQC_Error_A("Do not know how to allocate array", Iout, &
                 "Shape",Shape)
         Endif
         IF (AllocateStatus /= 0) then
            write( errorMsg, 1100 ) trim(Data_Name)
            call MQC_Error_I(errorMsg, Iout, 'LR', LR)
         EndIf
         call Rd_RBuf_G(wavefunction%unitNumber,NTot,LenBuf, &
              wavefunction%Array_Data%last%ArrR)
      case('COMPLEX')
         If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
              Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
              Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
            allocate(wavefunction%Array_Data%last%ArrC(LR,1,1,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."MATRIX") then
            allocate(wavefunction%Array_Data%last%ArrC(N1,N2,1,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."MULTISYMMATRICES") then
            allocate(wavefunction%Array_Data%last%ArrC(((N2+1)*N2)/2,N3,N4,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."3DMATRIX") then
            allocate(wavefunction%Array_Data%last%ArrC(N1,N2,N3,1,1), STAT = AllocateStatus)
         Else
            call MQC_Error_A("Do not know how to allocate array", Iout, &
                 "Shape",Shape)
         Endif
         IF (AllocateStatus /= 0) then
            write( errorMsg, 1100 ) trim(Data_Name)
            call MQC_Error_I(errorMsg, Iout, 'LR', LR)
         EndIf
         call Rd_CBuf_G(wavefunction%unitNumber,NTot,LenBuf, &
              wavefunction%Array_Data%last%ArrC)
      case('MIXED')
         If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
              Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
              Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
            allocate(wavefunction%Array_Data%last%ArrI(LR,1,1,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."MATRIX") then
            allocate(wavefunction%Array_Data%last%ArrI(N1,N2,1,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."MULTISYMMATRICES") then
            allocate(wavefunction%Array_Data%last%ArrI(((N2+1)*N2)/2,N3,N4,1,1), STAT = AllocateStatus)
         Else If ( Shape.eq."3DMATRIX") then
            allocate(wavefunction%Array_Data%last%ArrI(N1,N2,N3,1,1), STAT = AllocateStatus)
         Else
            call MQC_Error_A("Do not know how to allocate array", Iout, &
                 "Shape",Shape)
         Endif
         IF (AllocateStatus /= 0) then
            write( errorMsg, 1100 ) trim(Data_Name)
            call MQC_Error_I(errorMsg, Iout, 'LR', LR)
         EndIf
         if ( NR .ge. 0 ) then
            If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
                 Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
                 Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
               allocate(wavefunction%Array_Data%last%ArrR(LR,1,1,1,1), STAT = AllocateStatus)
            Else If ( Shape.eq."MATRIX") then
               allocate(wavefunction%Array_Data%last%ArrR(N1,N2,1,1,1), STAT = AllocateStatus)
            Else If ( Shape.eq."MULTISYMMATRICES") then
               allocate(wavefunction%Array_Data%last%ArrR(((N2+1)*N2)/2,N3,N4,1,1), STAT = AllocateStatus)
            Else If ( Shape.eq."3DMATRIX") then
               allocate(wavefunction%Array_Data%last%ArrR(N1,N2,N3,1,1), STAT = AllocateStatus)
            Else
               call MQC_Error_A("Do not know how to allocate array", Iout, &
                    "Shape",Shape)
            Endif
            IF (AllocateStatus /= 0) then
               write( errorMsg, 1100 ) trim(Data_Name)
               call MQC_Error_I(errorMsg, Iout, 'LR', LR)
            EndIf
            call Rd_RMixed(wavefunction%unitNumber,NTot,LenBuf, &
                 wavefunction%Array_Data%last%ArrI, &
                 wavefunction%Array_Data%last%ArrR)
         else
            If ( Shape.eq."SCALAR".or.Shape.eq."VECTOR".OR. &
                 Shape.eq."SYMMATRIX".or.Shape.eq."SYMSYM".OR. &
                 Shape.eq."TRISYMMATRICES".or.Shape.eq."UNKNOWN") then
               allocate(wavefunction%Array_Data%last%ArrC(LR,1,1,1,1), STAT = AllocateStatus)
            Else If ( Shape.eq."MATRIX") then
               allocate(wavefunction%Array_Data%last%ArrC(N1,N2,1,1,1), STAT = AllocateStatus)
            Else If ( Shape.eq."MULTISYMMATRICES") then
               allocate(wavefunction%Array_Data%last%ArrC(((N2+1)*N2)/2,N3,N4,1,1), STAT = AllocateStatus)
            Else If ( Shape.eq."3DMATRIX") then
               allocate(wavefunction%Array_Data%last%ArrC(N1,N2,N3,1,1), STAT = AllocateStatus)
            Else
               call MQC_Error_A("Do not know how to allocate array", Iout, &
                    "Shape",Shape)
            Endif
            IF (AllocateStatus /= 0) then
               write( errorMsg, 1100 ) trim(Data_Name)
               call MQC_Error_I(errorMsg, Iout, 'LR', LR)
            EndIf
            call Rd_CMixed(wavefunction%unitNumber,NTot,LenBuf, &
                 wavefunction%Array_Data%last%ArrI, &
                 wavefunction%Array_Data%last%ArrC)
         endif
      case('TENSOR')
         if(NRI.eq.1) then
            allocate(wavefunction%Array_Data%last%ArrR(LR,NR,1,1,1), &
                 STAT=AllocateStatus)
            IF (AllocateStatus /= 0) then
               write( errorMsg, 1100 ) trim(Data_Name)
               call MQC_Error_I(errorMsg, Iout, 'LR', LR, 'NR', NR)
            EndIf
            call Rd_2EN_G(wavefunction%unitNumber,NR,LR,NR*LR,NTot,LenBuf, &
                 wavefunction%Array_Data%last%ArrR)
!            call MQC_Matrix_SymmSymmR4Tensor_Put_Real(r4TensorOut,arrayTmp)
!            deallocate(arrayTmp)
         elseIf(NRI.eq.2) then
            allocate(wavefunction%Array_Data%last%ArrR(LR,NR,1,1,1), &
                 STAT = AllocateStatus)
            IF (AllocateStatus /= 0) then
               write( errorMsg, 1100 ) trim(Data_Name)
               call MQC_Error_I(errorMsg, Iout, 'LR', LR, 'NR', NR)
            EndIf
            call Rd_2EN_G(wavefunction%unitNumber,NR,LR,NR*LR, &
                 2*NTot,2*LenBuf,wavefunction%Array_Data%last%ArrR)
!!                call Rd_2EN(fileinfo%unitNumber,NR,NTot,NR*NTot,2*NTot,2*LenBuf,complexTmp)
!            complexTmp = reshape(arrayTmp, shape(complexTmp))
!            call MQC_Matrix_SymmSymmR4Tensor_Put_Complex(r4TensorOut,complexTmp)
!            !                deallocate(complexTmp)
!            deallocate(arrayTmp)
         else
            if ( NR .ge. 0 ) then
               allocate(wavefunction%Array_Data%last%ArrR(LR,1,1,1,1), &
                    STAT = AllocateStatus)
               IF (AllocateStatus /= 0) then
                  write( errorMsg, 1100 ) trim(Data_Name)
                  call MQC_Error_I(errorMsg, Iout, 'LR', LR, 'NR', NR)
               EndIf
               call Rd_RBuf_G(wavefunction%unitNumber,NTot,LenBuf, &
                    wavefunction%Array_Data%last%ArrR)
            else
               allocate(wavefunction%Array_Data%last%ArrC(LR,1,1,1,1), &
                    STAT = AllocateStatus)
               IF (AllocateStatus /= 0) then
                  write( errorMsg, 1100 ) trim(Data_Name)
                  call MQC_Error_I(errorMsg, Iout, 'LR', LR, 'NR', NR)
               EndIf
               call Rd_CBuf_G(wavefunction%unitNumber,NTot,LenBuf, &
                    wavefunction%Array_Data%last%ArrC)
            endif
         endIf
      case default
         write(IOut,1010) 'UnExpected data type, found:', &
              TRIM(Data_Name),NI,NR,NRI,NTot,LenBuf,  &
              N1,N2,N3,N4,N5,LR
         Call Rd_Skip(wavefunction%UnitNumber,NTot,LenBuf)
         IERROR = IERROR + 1
      end select
!
      return
    end subroutine MQC_Unformatted_Read_Array

!=====================================================================
!
!PROCEDURE MQC_PRINT_DATA_TYPE
      Subroutine MQC_PRINT_DATA_TYPE(Iout,DataName,NI,NR,N1,N2,N3,N4,N5,NRI,NTot,LR,DataType,Shape,Formula)
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
      character(len=*),intent(in)::DataName
      integer(kind=int64),intent(in)::NI,NR,N1,N2,N3,N4,N5,NRI
      integer(kind=int64),intent(in)::NTot,LR
      character(len=*),intent(in)::DataType,Shape,Formula
      integer(kind=int64)::i,size
      real(kind=real64)::a
      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
      character(len=8)::tab
!
!     Format statements.
!
 1010 Format( A )
 1020 Format( A,"  - ",A )
 1030 Format( 2A )
 1040 Format( A,"  - Total elements: ",a )
 1070 format('NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
 1100 Format('Error: Bad shape for',a,': Calculated size: ',i9,' Actual size:' &
          ,i9,/,'  NI =',I9,/,'  NR =',I9,/,'  NRI =',I9,/,&
          '  LenBuf =',I9,/, '  N1 =',I9,/, '  N2 =',I9,/ &
          '  N3 =',I9,/, '  N4 =',I9,/, '  N5 =',I9 )
 2010 Format( i15 )
!
!     Do the work...
!
      tab = char( 9 )
      write( iout, 1010) DataName
      write( iout, 1020) tab, trim(DataType)
      write( OutMath1, 2010) NTot
      call mqc_StripSpaces(OutMath1)
      write( iout, 1040) tab, trim(OutMath1)
      write( iout, 1030) tab, trim(formula)
!
      return
      end subroutine MQC_PRINT_DATA_TYPE

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
      write( IOut, 2060) "ATOMIC NUMBERS", tab, "INTEGER", tab, trim(OutMath2), tab, trim(OutMath2)

! current_wavefunction%atomTypes
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( IOut, 2060) "ATOM TYPES", tab, "INTEGER", tab, trim(OutMath2), tab, trim(OutMath2)

! current_wavefunction%atomicCharges
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( IOut, 2070) "ATOM CHARGES", tab, "REAL", tab, trim(OutMath2), tab, trim(OutMath2), tab, &
           'Nuclear charges may be different from atomic'//NEW_LINE('A')//tab//' &
           & numbers if ECPs were used.'

! current_wavefunction%cartesians
      write( OutMath1, 2010 ) current_wavefunction%natoms*3
      call mqc_StripSpaces(OutMath1)
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( OutMath3, 2010 ) 3
      call mqc_StripSpaces(OutMath3)
      write( IOut, 2080) "CARTESIAN NUCLEAR COORDINATES(BOHR)", tab, "REAL", &
           tab, trim(OutMath1), tab, trim(OutMath2), trim(OutMath3)

! First half of current_wavefunction%basisFunction2Atom
      write( OutMath1, 2010 ) current_wavefunction%nbasis
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2060) "MAP FROM BASIS FUNCTIONS TO ATOMS", tab, &
           "INTEGER", tab, trim(OutMath1), tab, trim(OutMath1)

! Second half half of current_wavefunction%basisFunction2Atom
      write( OutMath1, 2010 ) current_wavefunction%nbasis
      call mqc_StripSpaces(OutMath1)
      write( IOut, 2060) "TYPE FLAG FOR EACH BASIS FUNCTION", tab, &
           "INTEGER", tab, trim(OutMath1), tab, trim(OutMath1)

!           current_wavefunction%atomicWeights
      write( OutMath2, 2010 ) current_wavefunction%natoms
      call mqc_StripSpaces(OutMath2)
      write( IOut, 2060) "ATOMIC WEIGHTS", tab, "REAL", tab, &
           trim(OutMath2), tab, trim(OutMath2)
      return
    end subroutine MQC_PRINT_DEFAULT_VECTOR_SUMMARY
!=====================================================================
!

      subroutine MQC_PRINT_SCALAR_HEAD( iOut, wavefunction)
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
      integer(kind=int64),intent(in)::iOut
!
!     Format statements.
!
 1000 Format("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&",/, &
           1x,A,/,"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
 1010 Format(1x,A,2x,I14)

      write( IOut, 1000) "Contents of Matrix File"
      write( IOut, 1010) "NUMBER OF ATOMS", wavefunction%natoms
      write( IOut, 1010) "NUMBER OF BASIS FUNCTIONS", wavefunction%nbasis
      write( IOut, 1010) "NUMBER OF LINEARLY INDEPENDENT BASIS FUNCTIONS", &
           wavefunction%nbasisUse
      write( IOut, 1010) "MOLECULAR CHARGE", wavefunction%iCharge
      write( IOut, 1010) "SPIN MULTIPLICITY", wavefunction%Multiplicity
      write( IOut, 1010) "NUMBER OF ELECTRONS", wavefunction%nelectrons
      write( IOut, 1010) 'CLOSED/OPEN-SHELL FLAG', wavefunction%IOpCl
      write( IOut, 1010) "ENCODED INDICATOR OF WHETHER THE CALCULATION IS COMPLEX AND/OR GHF", &
           wavefunction%ICGU
      write( IOut, 1010) "NUMBER OF FROZEN CORE ORBITALS", wavefunction%NFC
      write( IOut, 1010) "NUMBER OF FROZEN VIRTUAL ORBITALS", wavefunction%NFV
      write( IOut, 1010) "NUMBER OF CONTRACTED SHELLS OF AO BASIS FUNCTIONS", &
           wavefunction%NShlAO
      write( IOut, 1010) "NUMBER OF PRIMITIVE AO SHELLS", wavefunction%NPrmAO
      write( IOut, 1010) "NUMBER OF CONTRACTED SHELLS OF DENSITY FITTING FUNCTIONS", &
           wavefunction%NShlDB
      write( IOut, 1010) "NUMBER OF PRIMITIVE DENSITY FITTING SHELLS", &
           wavefunction%NPrmDB
      write( IOut, 1010) "TOTAL NUMBER OF BONDS IN CONNECTIVITY DATA", &
           wavefunction%NBTot
      write( IOut, 1010) "ITRAN", wavefunction%ITran
      return
    end subroutine MQC_PRINT_SCALAR_HEAD

    Subroutine MQC_Create_Wavefunctions( Print_Options, my_READONLY_list, &
         mqc_input_file_list, Iout)
      implicit none

      integer(kind=int64),intent(in):: Print_Options
      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
      class(mqc_link_list),intent(inout) :: mqc_input_file_list
      integer(kind=int64),intent(in)::iOut

      character(len=1024)::program_name
      character(len=1024)::no_program_name
      character(len=1024)::FileName
      integer(kind=int64) cmd_count, I
      integer(kind=int64) status

!
!     Format statements.
!
 1010 Format( A )
!
!     Do the work...
!
      cmd_count = COMMAND_ARGUMENT_COUNT()
      if ( cmd_count .eq. 0 ) then
         write( iout, 1010 ) "Error in number of argments to this program."
         write( iout, 1010 ) "Arguments should be:"
         write( iout, 1010 ) "     - Full path to a Gaussian file"
         write( iout, 1010 ) "        - This can either be FormChk File"
         write( iout, 1010 ) "        - This can either be a MatrixFile"
         write( iout, 1010 ) "        - This can be a Gaussian input file"
         write( iout, 1010 ) "     - Optional, only required if the first "
         write( iout, 1010 ) "          argument is a Gaussian input file."
         write( iout, 1010 ) "          The name of the binary for the "
         write( iout, 1010 ) "          program that generates the"
         write( iout, 1010 ) "          FormChk/MatrixFile, most likely"
         write( iout, 1010 ) "          'g16'. The environment should be"
         write( iout, 1010 ) "          setup to run the program."
         call mqc_error_i( "No arguments to this program", iout)
      endif
!
!     Create a list of MatrixFiles.  Use MatrixFiles supplied.
!     If Gaussian input has been supplied:
!     1. Modify the input include MatrixFile options
!     2. Run the job to generate MatrixFiles
!     If the name of the Gaussian binary is not given, and Gaussian is not
!     in the path, this will fail.
!
      if ( cmd_count.eq. 1 ) then
         no_program_name = 'Unknown'
         call get_command_argument(1, FileName)
         call mqc_create_File_Name_lists_F2C(FileName, no_program_name, iout)
      else
         call get_command_argument(cmd_count, program_name)
         call mqc_File_Is_Executable_F2C( program_name, Status, iout )
         if ( Status .eq. 1 ) then
            do i = 1, cmd_count-1
               call get_command_argument(i, FileName)
               call mqc_create_File_Name_lists_F2C(FileName,Program_name,iout)
            end do
         else
            no_program_name = 'Unknown'
            do i = 1, cmd_count
               call get_command_argument(i, FileName)
               call mqc_create_File_Name_lists_F2C(FileName,no_program_name,iout)
            end do
         endif
      endif
#ifdef MATRIXFILE
! 
! Figure out is a MatrixFile name was supplied.
! Or if an input file has been supplied, and then run the job.
! Whatever has happened, get a list of the MatrixFiles
!
      call  mqc_get_MatrixFile_Name_F2C(FileName, iout)
      if ( trim(FileName).eq.'Done with MatrixFiles' ) then
         call mqc_error_I('No MatrixFiles supplied.  (Remove this when FormCHK files are implimented', IOut)
      endif
!  
      Do While( trim(FileName).ne.'Done with MatrixFiles' )
!
! Read the MatrixFile and put data into a wavefunction structure
!
         call mqc_link_list_add(mqc_input_file_list,Print_Options, &
              my_READONLY_list,FileName,iout)
!
! Get the next MatrixFile Name
!
         call mqc_get_MatrixFile_Name_F2C(FileName, iout)
      EndDo
#endif
      return
    end Subroutine MQC_Create_Wavefunctions

!     PROCEDURE MQC_GET_NUCLEAR_REPULSION
      Function MQC_Wave_To_Nuclear_Repulsion(wavefunction,IOut) Result(Vnn)
!
!     This function returns the nuclear repulsion energy given the
!     molecular data object.
!
!     Lee M. Thompson, 2017.
!
!
      Implicit None
      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
      real(kind=real64) Vnn
      type(mqc_array_link_list_node),pointer::Coord
      type(mqc_array_link_list_node),pointer::Charges
      Real::Zero=0.0d0
      Real::DX,DY,DZ
      Integer::IOut,I,J
      Integer(kind=int64)::NAtoms
      Integer(kind=int64)::Found_Cart, Found_Charge
!
      Vnn = Zero
      NAtoms = wavefunction%NAtoms
      Found_Cart = MQC_Find_mqc_array_link_list_node(wavefunction,&
           "CARTESIAN NUCLEAR COORDINATES(BOHR)", Coord, iout)
      Found_Charge = MQC_Find_mqc_array_link_list_node(wavefunction,&
           "ATOM CHARGES", Charges, iout)
      If ( Found_Cart .ne. 2 .or. Found_Charge .ne. 2 ) then
         call mqc_error_I('MQC_Wave_To_Nuclear_Repulsion: Did not find all compnents of Nuclear Replusion', IOut, &
         'Found_Cart', Found_Cart, 'Found_Charge',Found_Charge )
      EndIf
      Do I = 1,NAtoms-1
        Do J = I+1,NAtoms
           DX = Coord%ArrR(1,I,1,1,1) - Coord%ArrR(1,J,1,1,1)
           DY = Coord%ArrR(2,I,1,1,1) - Coord%ArrR(2,J,1,1,1)
           DZ = Coord%ArrR(3,I,1,1,1) - Coord%ArrR(3,J,1,1,1)
           Vnn = Vnn + Charges%ArrR(I,1,1,1,1) * Charges%ArrR(J,1,1,1,1) / &
                SQRT(DX*DX+DY*DY+DZ*DZ)
        EndDo
      EndDo
      End Function MQC_Wave_To_Nuclear_Repulsion

      Function mqc_read_this_MatrixFile_Element(my_READONLY_list,Data_Name ) &
           Result( ReadIt )
      class(mqc_READONLY_list),intent(inout):: my_READONLY_list
      character(len=*),intent(in)::Data_Name
      logical ReadIt
      integer(kind=int64) Status
      character(len=64)::List_Name

      If (.not.associated(my_READONLY_list%first)) then
         ReadIt = .true.
         return
      else
         call mqc_READONLY_list_start(my_READONLY_list,Status,List_Name)
         Do While( Status.eq.1 )
            If ( trim(Data_Name) .eq. trim(List_Name) ) then
               ReadIt = .true.
               return
            EndIf
            call mqc_READONLY_list_next(my_READONLY_list,Status,List_Name)
         EndDo
      EndIf
      ReadIt = .false.

      Return
    end Function mqc_read_this_MatrixFile_Element

    End module MQC_FullWavefunction


