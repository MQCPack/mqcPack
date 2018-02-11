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
     private
     character(Len=64) :: Data_Name
     integer(kind=int64) :: NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR
     character(Len=64):: Data_Type,Shape
     character(Len=1024):: Formula
     Real(Kind=real64),Dimension(:),Allocatable::VecR
     Integer(kind=int64),Dimension(:),Allocatable::VecI
     Complex(Kind=real64),Dimension(:),Allocatable::VecC
     type(mqc_array_link_list_node),pointer :: next =>null()
  end type mqc_array_link_list_node

  type :: mqc_array_link_list
     private
     type(mqc_array_link_list_node),pointer,private :: first => null()
     type(mqc_array_link_list_node),pointer,private :: last => null()
     type(mqc_array_link_list_node),pointer,private :: current => null()
     integer(kind=int64),private :: length = 0
   contains
     procedure :: start => mqc_array_link_list_start
     procedure :: next => mqc_array_link_list_next
  end type mqc_array_link_list

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
     Logical::WF_Complex
     Integer(kind=int64)::NCore,NVal,NActive
     real(kind=real64),dimension(:),allocatable::atomicCharges,atomicWeights,cartesians

     Type(mqc_array_link_list)::Array_Data

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
     procedure :: next => mqc_link_list_next
  end type mqc_link_list

CONTAINS
      subroutine mqc_read_all_wavefunction( current_wavefunction, Print_Options, iout)
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::current_wavefunction
      integer(kind=int64),intent(in)::Print_Options
      integer(kind=int64),intent(in)::iOut
      Integer(kind=int64):: NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,NRI,LR
      logical ASym
      Character(Len=64)::CBuf
      Character(Len=2048)::DataType,Shape,Formula
      logical EOF
      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
      integer(kind=int32)::i
      integer(kind=int64)::IERROR
!
!     Format statements.
!
 1000 format( "Warning: MQC Pack cannot process file '",a,"' it may have been created with the raw attribute.")
 1010 format('DataName  ',a,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
 1050 format( 'ERROR:',a,': please implement proper handling of this data type')
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
            call String_Change_Case(cBuf,'u')
            if(trim(cBuf).eq.'END') then
               Call Close_MatF(current_wavefunction%unitNumber)
               call mqc_print_all_wavefunction(current_wavefunction, IOut)
               if (IERROR .ne. 0 ) then
                  call mqc_error_I('mqc_read_all_wavefunction: Found errors reading Matrix Element File', 6)
               endif
               return
            endif

            call MQC_Unformatted_Read_Array(current_wavefunction,cBuf, &
                 NR,NI,NRI,NTot,N1,N2,N3,N4,N5,LR,LenBuf,DataType,Shape,&
                 Formula,IOut,IERROR)
!
         endif
      enddo
      return
      end subroutine mqc_read_all_wavefunction
!
!     PROCEDURE MQC_Print_all_Wavefunction     
      subroutine mqc_print_all_wavefunction(wavefunction,iOut)
!
      implicit none
      class(MQC_ALL_Wavefunction)::wavefunction
      integer(kind=int64),intent(in)::iOut
      integer(kind=int64)::LR
!                                                                                    write( iout, 1050) "In mqc_print_all_wavefunction"
 1050 Format( 2A )
 1060 Format( A,L10 )

      call MQC_PRINT_SCALAR_HEAD( iOut, wavefunction )
      call MQC_PRINT_DEFAULT_VECTOR_HEAD( iOut, wavefunction)
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
 2011 Format(1x,I7,10I14)
 3010 Format("  Dimension 3:",I5 )
 3020 Format("  Dimension 3:",I5, "  Dimension 4:",I5 )
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
         Write(IOut,1010) TRIM(Array%Data_Name), Array%VecI(1)
      else if ( Array%Shape.eq."VECTOR") then
         Do I = 1, Array%LR
            Write(IOut,1001) I,  Array%VecI(I)
         EndDo
      else if ( Array%Shape.eq."MATRIX") then
         NRows = Array%N1
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = 1,NRows
               Write(IOut,2011) I, (Array%VecI(I+(J-1)*NRows),J=IFirst,ILast)
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
               Write(IOut,2011) I, (Array%VecI(IBase+J),J=IFirst,I)
            EndDo
            Do I = ILAST+1, NCols
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (Array%VecI(IBase+J),J=IFirst,ILast)
            EndDo
         EndDo
      else if ( Array%Shape.eq."MULTISYMMATRICES") then
         NRows = Array%N2
         NCols = Array%N2
         N3 = Array%N3
         N4 = Array%N4
         Do L = 1,N4
            Do K = 1,N3
               IBase = ((NCols+1)*NCols/2)*((L-1)*N3+K-1)
               Do IFirst = 1,NCols,ColWidth
                  ILast = Min(IFirst+ColWidth-1,NCols)
                  if ( N4 .eq. 1 ) then
                     write(Iout,3010) K
                     Write(IOut,2001) (I,I=IFirst,ILast)
                  else 
                     write(Iout,3020) K, L
                     Write(IOut,2001) (I,I=IFirst,ILast)
                  endif
                  Do I = IFirst, ILAST
                     JBase = IBase+(I*(I-1))/2
                     Write(IOut,2011) I, (Array%VecI(JBase+J),J=IFirst,I)
                  EndDo
                  Do I = ILAST+1, NCols
                     Write(IOut,2011) I, (Array%VecI(JBase+J),J=IFirst,ILast)
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
                     write(Iout,4010) I,J,K,L,Array%VecI(IBase)
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
                  write(Iout,4020) I,J,K,Array%VecI(IBase)
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."3DMATRIX") then
         N1 = Array%N1
         N2 = Array%N2
         N3 = Array%N3
         IBase=0
         Do I = 1,N1
            Do J = 1,N2
               Do K = 1,N3
                  IBase=IBase+1
                  write(Iout,4020) I,J,K,Array%VecI(IBase)
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
 2011 Format(1x,I7,10F14.6)
 3010 Format("  Dimension 3:",I5 )
 3020 Format("  Dimension 3:",I5, "  Dimension 4:",I5 )
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
         Write(IOut,1010) TRIM(Array%Data_Name), Array%VecR(1)
      else if ( Array%Shape.eq."VECTOR") then
         Do I = 1, Array%LR
            Write(IOut,1001) I,  Array%VecR(I)
         EndDo
      else if ( Array%Shape.eq."MATRIX") then
         NRows = Array%N1
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = 1,NRows
               Write(IOut,2011) I, (Array%VecR(I+(J-1)*NRows),J=IFirst,ILast)
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
               Write(IOut,2011) I, (Array%VecR(IBase+J),J=IFirst,I)
            EndDo
            Do I = ILAST+1, NCols
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (Array%VecR(IBase+J),J=IFirst,ILast)
            EndDo
         EndDo
      else if ( Array%Shape.eq."MULTISYMMATRICES") then
         NRows = Array%N2
         NCols = Array%N2
         N3 = Array%N3
         N4 = Array%N4
         Do L = 1,N4
            Do K = 1,N3
               IBase = ((NCols+1)*NCols/2)*((L-1)*N3+K-1)
               Do IFirst = 1,NCols,ColWidth
                  ILast = Min(IFirst+ColWidth-1,NCols)
                  if ( N4 .eq. 1 ) then
                     write(Iout,3010) K
                     Write(IOut,2001) (I,I=IFirst,ILast)
                  else 
                     write(Iout,3020) K, L
                     Write(IOut,2001) (I,I=IFirst,ILast)
                  endif
                  Do I = IFirst, ILAST
                     JBase = IBase+(I*(I-1))/2
                     Write(IOut,2011) I, (Array%VecR(JBase+J),J=IFirst,I)
                  EndDo
                  Do I = ILAST+1, NCols
                     Write(IOut,2011) I, (Array%VecR(JBase+J),J=IFirst,ILast)
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
                     write(Iout,4010) I,J,K,L,Array%VecR(IBase)
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
                  write(Iout,4020) I,J,K,Array%VecR(IBase)
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."3DMATRIX") then
         N1 = Array%N1
         N2 = Array%N2
         N3 = Array%N3
         IBase=0
         Do I = 1,N1
            Do J = 1,N2
               Do K = 1,N3
                  IBase=IBase+1
                  write(Iout,4020) I,J,K,Array%VecR(IBase)
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
 2011 Format(1x,I7," (",F13.6,",",F13.6,"i)  (",F13.6,",",F13.6,"i)  (", &
           F13.6,",",F13.6,"i)  (",F13.6,",",F13.6,"i)  (",F13.6,",",F13.6,"i)")
 3010 Format("  Dimension 3:",I5 )
 3020 Format("  Dimension 3:",I5, "  Dimension 4:",I5 )
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
         Write(IOut,1010) TRIM(Array%Data_Name), REAL(Array%VecC(I)), &
              AIMAG(Array%VecC(I))
      else if ( Array%Shape.eq."VECTOR") then
         Do I = 1, Array%LR
            Write(IOut,1001) I, REAL(Array%VecC(I)), AIMAG(Array%VecC(I))
         EndDo
      else if ( Array%Shape.eq."MATRIX") then
         NRows = Array%N1
         NCols = Array%N2
         Do IFirst = 1,NCols,ColWidth
            ILast = Min(IFirst+ColWidth-1,NCols)
            Write(IOut,2001) (I,I=IFirst,ILast)
            Do I = 1,NRows
               Write(IOut,2011) I, (REAL(Array%VecC(I+(J-1)*NRows)), &
                    AIMAG(Array%VecC(I+(J-1)*NRows)),J=IFirst,ILast)
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
               Write(IOut,2011) I, (REAL(Array%VecC(IBase+J)), &
                    AIMAG(Array%VecC(IBase+J)),J=IFirst,I)
            EndDo
            Do I = ILAST+1, NCols
               IBase = (I*(I-1))/2
               Write(IOut,2011) I, (REAL(Array%VecC(IBase+J)), &
                    AIMAG(Array%VecC(IBase+J)),J=IFirst,ILAST)
            EndDo
         EndDo
      else if ( Array%Shape.eq."MULTISYMMATRICES") then
         NRows = Array%N2
         NCols = Array%N2
         N3 = Array%N3
         N4 = Array%N4
         Do L = 1,N4
            Do K = 1,N3
               IBase = ((NCols+1)*NCols/2)*((L-1)*N3+K-1)
               Do IFirst = 1,NCols,ColWidth
                  ILast = Min(IFirst+ColWidth-1,NCols)
                  if ( N4 .eq. 1 ) then
                     write(Iout,3010) K
                     Write(IOut,2001) (I,I=IFirst,ILast)
                  else 
                     write(Iout,3020) K, L
                     Write(IOut,2001) (I,I=IFirst,ILast)
                  endif
                  Do I = IFirst, ILAST
                     JBase = IBase+(I*(I-1))/2
                     Write(IOut,2011) I, (REAL(Array%VecC(JBase+J)), &
                          AIMAG(Array%VecC(JBase+J)),J=IFirst,I)
                  EndDo
                  Do I = ILAST+1, NCols
                     Write(IOut,2011) I, (REAL(Array%VecC(JBase+J)), &
                          AIMAG(Array%VecC(JBase+J)),J=IFirst,ILAST)
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
                     write(Iout,4010) I,J,K,L,REAL(Array%VecC(IBase)), &
                          AIMAG(Array%VecC(IBase))
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
                  write(Iout,4020) I,J,K,REAL(Array%VecC(IBase)), &
                       AIMAG(Array%VecC(IBase))
               EndDo
            EndDo
         EndDo
      else if ( Array%Shape.eq."3DMATRIX") then
         N1 = Array%N1
         N2 = Array%N2
         N3 = Array%N3
         IBase=0
         Do I = 1,N1
            Do J = 1,N2
               Do K = 1,N3
                  IBase=IBase+1
                  write(Iout,4020) I,J,K,REAL(Array%VecC(IBase)), &
                       AIMAG(Array%VecC(IBase))
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
      subroutine mqc_link_list_add(self,Print_Options,FileName,iout)
        class(mqc_link_list),intent(inout) :: self
        integer(kind=int64),intent(in)::Print_Options
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
             Print_Options,iout)

        self%length = self%length + 1
      end subroutine mqc_link_list_add

      subroutine mqc_link_list_start(self)
        class(mqc_link_list),intent(inout) :: self

        self%current => self%first
      end subroutine mqc_link_list_start

      subroutine mqc_link_list_next(self,wavefunction)
        class(mqc_link_list),intent(inout)::self
        class(MQC_ALL_Wavefunction),pointer,intent(out)::wavefunction

        if (associated(self%current)) then
           wavefunction => self%current%current_wavefunction
           self%current => self%current%next
        else
           wavefunction => null()
        end if
      end subroutine mqc_link_list_next

!=====================================================================
!
!PROCEDURE MQC_Unformatted_Read_Array
      subroutine MQC_Unformatted_Read_Array(wavefunction,Data_Name, &
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
      call mqc_array_link_list_add(wavefunction,Data_Name,NR,NI,NRI,NTot, & 
           N1,N2,N3,N4,N5,LR,IOut,Data_Type,Shape,Formula)
!
      select case(Data_Type)
      case('INTEGER')
         allocate(wavefunction%Array_Data%last%VecI(LR), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) then
            write( errorMsg, 1100 ) trim(Data_Name)
            call MQC_Error_I(errorMsg, Iout, 'LR', LR)
         EndIf
         call Rd_IBuf(wavefunction%unitNumber,NTot,LenBuf, &
              wavefunction%Array_Data%last%VecI)
      case('REAL')
         allocate(wavefunction%Array_Data%last%VecR(LR), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) then
            write( errorMsg, 1100 ) trim(Data_Name)
            call MQC_Error_I(errorMsg, Iout, 'LR', LR)
         EndIf
         call Rd_RBuf(wavefunction%unitNumber,NTot,LenBuf, &
              wavefunction%Array_Data%last%VecR)
      case('COMPLEX')
         allocate(wavefunction%Array_Data%last%VecC(LR), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) then
            write( errorMsg, 1100 ) trim(Data_Name)
            call MQC_Error_I(errorMsg, Iout, 'LR', LR)
         EndIf
         call Rd_CBuf(wavefunction%unitNumber,NTot,LenBuf, &
              wavefunction%Array_Data%last%VecC)
      case('MIXED')
         allocate(wavefunction%Array_Data%last%VecI(LR), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) then
            write( errorMsg, 1100 ) trim(Data_Name)
            call MQC_Error_I(errorMsg, Iout, 'LR', LR)
         EndIf
         if ( NR .ge. 0 ) then
            allocate(wavefunction%Array_Data%last%VecR(LR), &
                 STAT = AllocateStatus)
            IF (AllocateStatus /= 0) then
               write( errorMsg, 1100 ) trim(Data_Name)
               call MQC_Error_I(errorMsg, Iout, 'LR', LR)
            EndIf
            call Rd_RMixed(wavefunction%unitNumber,NTot,LenBuf, &
                 wavefunction%Array_Data%last%VecI, &
                 wavefunction%Array_Data%last%VecR)
         else
            allocate(wavefunction%Array_Data%last%VecC(LR), &
                 STAT = AllocateStatus)
            IF (AllocateStatus /= 0) then
               write( errorMsg, 1100 ) trim(Data_Name)
               call MQC_Error_I(errorMsg, Iout, 'LR', LR)
            EndIf
            call Rd_CMixed(wavefunction%unitNumber,NTot,LenBuf, &
                 wavefunction%Array_Data%last%VecI, &
                 wavefunction%Array_Data%last%VecC)
         endif
      case('TENSOR')
         if(NRI.eq.1) then
            allocate(wavefunction%Array_Data%last%VecR(LR*NR), &
                 STAT=AllocateStatus)
            IF (AllocateStatus /= 0) then
               write( errorMsg, 1100 ) trim(Data_Name)
               call MQC_Error_I(errorMsg, Iout, 'LR', LR, 'NR', NR)
            EndIf
            call Rd_2EN(wavefunction%unitNumber,NR,LR,NR*LR,NTot,LenBuf, &
                 wavefunction%Array_Data%last%VecR)
!            call MQC_Matrix_SymmSymmR4Tensor_Put_Real(r4TensorOut,arrayTmp)
!            deallocate(arrayTmp)
         elseIf(NRI.eq.2) then
            allocate(wavefunction%Array_Data%last%VecR(LR*NR), &
                 STAT = AllocateStatus)
            IF (AllocateStatus /= 0) then
               write( errorMsg, 1100 ) trim(Data_Name)
               call MQC_Error_I(errorMsg, Iout, 'LR', LR, 'NR', NR)
            EndIf
            call Rd_2EN(wavefunction%unitNumber,NR,LR,NR*LR, &
                 2*NTot,2*LenBuf,wavefunction%Array_Data%last%VecR)
!!                call Rd_2EN(fileinfo%unitNumber,NR,NTot,NR*NTot,2*NTot,2*LenBuf,complexTmp)
!            complexTmp = reshape(arrayTmp, shape(complexTmp))
!            call MQC_Matrix_SymmSymmR4Tensor_Put_Complex(r4TensorOut,complexTmp)
!            !                deallocate(complexTmp)
!            deallocate(arrayTmp)
         else
            if ( NR .ge. 0 ) then
               allocate(wavefunction%Array_Data%last%VecR(LR), &
                    STAT = AllocateStatus)
               IF (AllocateStatus /= 0) then
                  write( errorMsg, 1100 ) trim(Data_Name)
                  call MQC_Error_I(errorMsg, Iout, 'LR', LR, 'NR', NR)
               EndIf
               call Rd_RBuf(wavefunction%unitNumber,NTot,LenBuf, &
                    wavefunction%Array_Data%last%VecR)
            else
               allocate(wavefunction%Array_Data%last%VecC(LR), &
                    STAT = AllocateStatus)
               IF (AllocateStatus /= 0) then
                  write( errorMsg, 1100 ) trim(Data_Name)
                  call MQC_Error_I(errorMsg, Iout, 'LR', LR, 'NR', NR)
               EndIf
               call Rd_CBuf(wavefunction%unitNumber,NTot,LenBuf, &
                    wavefunction%Array_Data%last%VecC)
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

      subroutine MQC_PRINT_DEFAULT_VECTOR_HEAD( iOut, wavefunction)
      implicit none
      class(MQC_ALL_Wavefunction),intent(inout)::wavefunction
      integer(kind=int64),intent(in)::iOut
      character(len=64)::OutMath1,OutMath2,OutMath3,OutMath4
      Integer(kind=int64),Parameter::ColWidth=10 
      Integer(kind=int64)::I,J,NCols,NRows,IFirst,ILast
!
!     Format statements.
!
 1000 Format(/,1x,A)
 1001 Format(1x,I7,2x,I14)
 2001 Format(1x,I7,2x,F14.6)
 2021 Format(5x,5(20x,I7))
 2011 Format(1x,I7,10F14.6)

! current_wavefunction%atomicNumbers
      write( IOut, 1000 ) "ATOMIC NUMBERS"
      Do I = 1, wavefunction%natoms
         Write(IOut,1001) I, wavefunction%atomicNumbers(I)
      EndDo

      write( IOut, 1000 ) "ATOMIC TYPES"
      Do I = 1, wavefunction%natoms
         Write(IOut,1001) I, wavefunction%atomTypes(I)
      EndDo

      write( IOut, 1000 ) "ATOM CHARGES"
      Do I = 1, wavefunction%natoms
         Write(IOut,2001) I, wavefunction%atomicCharges(I)
      EndDo

      write( IOut, 1000 ) "CARTESIAN NUCLEAR COORDINATES IN BOHR"
      NRows = wavefunction%natoms
      NCols = 3
      Do IFirst = 1,NCols,ColWidth
         ILast = Min(IFirst+ColWidth-1,NCols)
         Write(IOut,2021) (I,I=IFirst,ILast)
         Do I = 1,NRows
            Write(IOut,2011) I, (wavefunction%cartesians(I+(J-1)*NRows), &
                 J=IFirst,ILast)
         EndDo
      EndDo

      write( IOut, 1000 ) "MAP FROM BASIS FUNCTIONS TO ATOMS"
      Do I = 1, wavefunction%nbasis
         Write(IOut,1001) I, wavefunction%basisFunction2Atom(I)
      EndDo

      write( IOut, 1000 ) "TYPE FLAG FOR EACH BASIS FUNCTION"
      Do I = 1, wavefunction%nbasis
         Write(IOut,1001) I, wavefunction%IBasisFunctionType(I)
      EndDo

      write( IOut, 1000 ) "ATOMIC WEIGHTS"
      Do I = 1, wavefunction%natoms
         Write(IOut,2001) I, wavefunction%atomicWeights(I)
      EndDo
      return
    end subroutine MQC_PRINT_DEFAULT_VECTOR_HEAD

    End module MQC_FullWavefunction


