      module MQC_MatWrapper
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
!     **                          Version 26.8.1                           **
!     **                          August 31, 2026                          **
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
!     This module preserves the MQC_Gaussian compile-time MatrixFile boundary
!     in builds configured without GauOpen. Formatted-checkpoint operations do
!     not call these routines; every MatrixFile operation reports one clear
!     unsupported-feature error.
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env, only: int32, int64, real64
      use MQC_General, only: mqc_error
      implicit none
      private
      public::MQC_MatrixFile_Layout
      public::MQC_MatrixFile_Character_Record_Layout
      public::MQC_MatrixFile_Get_Character_Record_Layout
      public::MQC_MatrixFile_Get_Character_Write_Layout
      public::MQC_MatrixFile_Open_Read,MQC_MatrixFile_Open_Write
      public::MQC_MatrixFile_Close,MQC_MatrixFile_Len_Array
      public::MQC_MatrixFile_Read_Header,MQC_MatrixFile_Write_Header
      public::MQC_MatrixFile_Read_Label,MQC_MatrixFile_Write_Label
      public::MQC_MatrixFile_Read_Character_Buffer
      public::MQC_MatrixFile_Write_Character_Buffer
      public::MQC_MatrixFile_Read_Integer_Buffer
      public::MQC_MatrixFile_Read_Real_Buffer
      public::MQC_MatrixFile_Read_Complex_Buffer
      public::MQC_MatrixFile_Read_Two_Electron
      public::MQC_MatrixFile_Read_Indexed_Real,MQC_MatrixFile_Skip
      public::MQC_MatrixFile_Write_AO_Two_Electron
      public::MQC_MatrixFile_Write_Real_Array
      public::MQC_MatrixFile_Write_Integer_Array
      public::MQC_MatrixFile_Write_Complex_Array
!
      character(len=*),parameter::MQC_NO_GAUOPEN_MESSAGE =  &
        'MQCPack was configured without GauOpen support; MatrixFile '//  &
        'operations are unavailable.'
!
      type MQC_MatrixFile_Layout
        integer(kind=int32)::integerBytes=0
        integer(kind=int64)::fileVersion=0
        integer(kind=int64)::nHeaderRecords=0
        integer(kind=int64)::len12Bytes=0
        integer(kind=int64)::len4Bytes=0
        logical::rawFormat=.false.
      contains
        procedure,public::set => MQC_MatrixFile_Layout_Set
        procedure,public::setForWrite => MQC_MatrixFile_Layout_Set_For_Write
      end type MQC_MatrixFile_Layout
!
      type MQC_MatrixFile_Character_Record_Layout
        integer(kind=int64)::elementCount=0
        integer(kind=int64)::elementWidth=0
        integer(kind=int64)::payloadCharacters=0
        integer(kind=int64)::bufferCharacters=0
        integer(kind=int64)::payloadWords=0
        integer(kind=int64)::bufferWords=0
      end type MQC_MatrixFile_Character_Record_Layout
!
      contains
!
!PROCEDURE MQC_MatrixFile_Unsupported
      subroutine MQC_MatrixFile_Unsupported()
!
!     Report use of MatrixFile functionality from a FormChk-only build.
!
!     H. P. Hratchian, 2026.
!
      implicit none
!
      call mqc_error(MQC_NO_GAUOPEN_MESSAGE)
!
      return
      end subroutine MQC_MatrixFile_Unsupported
!
!PROCEDURE MQC_MatrixFile_Layout_Set
      subroutine MQC_MatrixFile_Layout_Set(layout,fileVersion,  &
        nHeaderRecords,len12Bytes,len4Bytes)
      implicit none
      class(MQC_MatrixFile_Layout),intent(inout)::layout
      integer(kind=int64),intent(in)::fileVersion,nHeaderRecords,  &
        len12Bytes,len4Bytes
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Layout_Set
!
!PROCEDURE MQC_MatrixFile_Layout_Set_For_Write
      subroutine MQC_MatrixFile_Layout_Set_For_Write(layout)
      implicit none
      class(MQC_MatrixFile_Layout),intent(inout)::layout
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Layout_Set_For_Write
!
!PROCEDURE MQC_MatrixFile_Get_Character_Record_Layout
      function MQC_MatrixFile_Get_Character_Record_Layout(layout,nTot,  &
        lenBuf,n1,n2,n3,n4,n5,typeA) result(recordLayout)
      implicit none
      type(MQC_MatrixFile_Layout),intent(in)::layout
      integer(kind=int64),intent(in)::nTot,lenBuf,n1,n2,n3,n4,n5,typeA
      type(MQC_MatrixFile_Character_Record_Layout)::recordLayout
!
      call MQC_MatrixFile_Unsupported()
      end function MQC_MatrixFile_Get_Character_Record_Layout
!
!PROCEDURE MQC_MatrixFile_Get_Character_Write_Layout
      function MQC_MatrixFile_Get_Character_Write_Layout(layout,  &
        elementCount,elementWidth,bufferWords) result(recordLayout)
      implicit none
      type(MQC_MatrixFile_Layout),intent(in)::layout
      integer(kind=int64),intent(in)::elementCount,elementWidth,bufferWords
      type(MQC_MatrixFile_Character_Record_Layout)::recordLayout
!
      call MQC_MatrixFile_Unsupported()
      end function MQC_MatrixFile_Get_Character_Write_Layout
!
!PROCEDURE MQC_MatrixFile_Open_Read
      subroutine MQC_MatrixFile_Open_Read(name,iUnit,labFil,iVers,nLab,  &
        gVers,title,nAtoms,  &
        nBasis,nBasisUse,iCharge,multiplicity,nElectrons,len12L,len4L,  &
        iOpCl,icgu)
      implicit none
      character(len=*),intent(in)::name
      character(len=64),intent(inout)::labFil,gVers,title
      integer(kind=int64),intent(out)::iUnit,iVers,nLab,nAtoms,nBasis,  &
        nBasisUse,iCharge,multiplicity,nElectrons,len12L,len4L,iOpCl,icgu
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Open_Read
!
!PROCEDURE MQC_MatrixFile_Open_Write
      subroutine MQC_MatrixFile_Open_Write(name,iUnit,labFil,gVers,title,  &
        nAtoms,nBasis,  &
        nBasisUse,iCharge,multiplicity,nElectrons,iOpCl,icgu)
      implicit none
      character(len=*),intent(in)::name
      character(len=64),intent(inout)::labFil,gVers,title
      integer(kind=int64),intent(out)::iUnit
      integer(kind=int64),intent(in)::nAtoms,nBasis,nBasisUse,iCharge,  &
        multiplicity,nElectrons,iOpCl,icgu
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Open_Write
!
!PROCEDURE MQC_MatrixFile_Close
      subroutine MQC_MatrixFile_Close(iUnit)
      implicit none
      integer(kind=int64),intent(in)::iUnit
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Close
!
!PROCEDURE MQC_MatrixFile_Len_Array
      function MQC_MatrixFile_Len_Array(n1,n2,n3,n4,n5)  &
        result(arrayLength)
      implicit none
      integer(kind=int64),intent(in)::n1,n2,n3,n4,n5
      integer(kind=int64)::arrayLength
!
      arrayLength = 0_int64
      call MQC_MatrixFile_Unsupported()
      end function MQC_MatrixFile_Len_Array
!
!PROCEDURE MQC_MatrixFile_Read_Header
      subroutine MQC_MatrixFile_Read_Header(iUnit,nLab,nAtoms,nBasis,  &
        atomicNumbers,atomTypes,  &
        atomicCharges,cartesians,basisFunction2Atom,basisFunctionTypes,  &
        atomicWeights,nfc,nfv,iTran,iDum9,nShellAO,nPrimitiveAO,nShellDB,  &
        nPrimitiveDB,nBasisTotal)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nLab,nAtoms,nBasis
      integer(kind=int64),dimension(nAtoms),intent(out)::atomicNumbers,  &
        atomTypes
      real(kind=real64),dimension(:),intent(out)::atomicCharges
      real(kind=real64),dimension(:),intent(inout)::cartesians,atomicWeights
      integer(kind=int64),dimension(nBasis),intent(out)::basisFunction2Atom,  &
        basisFunctionTypes
      integer(kind=int64),intent(out)::nfc,nfv,iTran,iDum9,nShellAO,  &
        nPrimitiveAO,nShellDB,nPrimitiveDB,nBasisTotal
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Read_Header
!
!PROCEDURE MQC_MatrixFile_Write_Header
      subroutine MQC_MatrixFile_Write_Header(iUnit,nAtoms,nAtomCoordinates,  &
        nBasis,atomicNumbers,  &
        atomTypes,atomicCharges,cartesians,basisFunction2Atom,  &
        basisFunctionTypes,atomicWeights,nfc,nfv,iTran,iDum9,nShellAO,  &
        nPrimitiveAO,nShellDB,nPrimitiveDB,nBasisTotal)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nAtoms,nAtomCoordinates,nBasis
      integer(kind=int64),dimension(nAtoms),intent(in)::atomicNumbers,  &
        atomTypes
      real(kind=real64),dimension(:),intent(in)::atomicCharges,cartesians,  &
        atomicWeights
      integer(kind=int64),dimension(nBasis),intent(in)::basisFunction2Atom,  &
        basisFunctionTypes
      integer(kind=int64),intent(in)::nfc,nfv,iTran,iDum9,nShellAO,  &
        nPrimitiveAO,nShellDB,nPrimitiveDB,nBasisTotal
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Write_Header
!
!PROCEDURE MQC_MatrixFile_Read_Label
      subroutine MQC_MatrixFile_Read_Label(iUnit,iVers,label,nInteger,  &
        nReal,nTotal,  &
        bufferLength,n1,n2,n3,n4,n5,typeA,nRealInteger,endOfFile)
      implicit none
      integer(kind=int64),intent(in)::iUnit,iVers
      character(len=64),intent(out)::label
      integer(kind=int64),intent(out)::nInteger,nReal,nTotal,bufferLength,  &
        n1,n2,n3,n4,n5,typeA,nRealInteger
      logical,intent(out)::endOfFile
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Read_Label
!
!PROCEDURE MQC_MatrixFile_Write_Label
      subroutine MQC_MatrixFile_Write_Label(iUnit,label,nInteger,nReal,  &
        nTotal,  &
        bufferLength,n1,n2,n3,n4,n5,typeA)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nInteger,nReal,nTotal,  &
        bufferLength,n1,n2,n3,n4,n5,typeA
      character(len=*),intent(in)::label
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Write_Label
!
!PROCEDURE MQC_MatrixFile_Read_Character_Buffer
      subroutine MQC_MatrixFile_Read_Character_Buffer(iUnit,nCharacters,  &
        bufferLength,characterArray)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nCharacters,bufferLength
      character(len=*),intent(inout)::characterArray
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Read_Character_Buffer
!
!PROCEDURE MQC_MatrixFile_Write_Character_Buffer
      subroutine MQC_MatrixFile_Write_Character_Buffer(iUnit,nCharacters,  &
        bufferLength,characterArray)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nCharacters,bufferLength
      character(len=*),intent(in)::characterArray
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Write_Character_Buffer
!
!PROCEDURE MQC_MatrixFile_Read_Integer_Buffer
      subroutine MQC_MatrixFile_Read_Integer_Buffer(iUnit,nTotal,  &
        bufferLength,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nTotal,bufferLength
      integer(kind=int64),dimension(nTotal),intent(out)::array
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Read_Integer_Buffer
!
!PROCEDURE MQC_MatrixFile_Read_Real_Buffer
      subroutine MQC_MatrixFile_Read_Real_Buffer(iUnit,nTotal,  &
        bufferLength,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nTotal,bufferLength
      real(kind=real64),dimension(:),intent(out)::array
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Read_Real_Buffer
!
!PROCEDURE MQC_MatrixFile_Read_Complex_Buffer
      subroutine MQC_MatrixFile_Read_Complex_Buffer(iUnit,nTotal,  &
        bufferLength,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nTotal,bufferLength
      complex(kind=real64),dimension(nTotal),intent(inout)::array
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Read_Complex_Buffer
!
!PROCEDURE MQC_MatrixFile_Read_Two_Electron
      subroutine MQC_MatrixFile_Read_Two_Electron(iUnit,nReal,recordLength,  &
        recordBlockLength,nTotal,  &
        bufferLength,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nReal,recordLength,  &
        recordBlockLength,nTotal,bufferLength
      real(kind=real64),dimension(:),intent(inout)::array
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Read_Two_Electron
!
!PROCEDURE MQC_MatrixFile_Read_Indexed_Real
      subroutine MQC_MatrixFile_Read_Indexed_Real(iUnit,nReal,recordLength,  &
        nTotal,bufferLength,  &
        nNonzero,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nReal,recordLength,nTotal,  &
        bufferLength
      integer(kind=int64),intent(out)::nNonzero
      real(kind=real64),dimension(:),intent(out)::array
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Read_Indexed_Real
!
!PROCEDURE MQC_MatrixFile_Skip
      subroutine MQC_MatrixFile_Skip(iUnit,nTotal,bufferLength)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nTotal,bufferLength
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Skip
!
!PROCEDURE MQC_MatrixFile_Write_AO_Two_Electron
      subroutine MQC_MatrixFile_Write_AO_Two_Electron(iUnit,label,nReal,  &
        bufferLength,nIntegrals,  &
        integrals)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::iUnit,nReal,bufferLength,nIntegrals
      real(kind=real64),dimension(:),intent(inout)::integrals
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Write_AO_Two_Electron
!
!PROCEDURE MQC_MatrixFile_Write_Real_Array
      subroutine MQC_MatrixFile_Write_Real_Array(iUnit,label,nReal,  &
        bufferLength,n1,n2,n3,n4,n5,  &
        asymmetric,array)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::iUnit,nReal,bufferLength,n1,n2,n3,  &
        n4,n5
      logical,intent(in)::asymmetric
      real(kind=real64),dimension(:),intent(inout)::array
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Write_Real_Array
!
!PROCEDURE MQC_MatrixFile_Write_Integer_Array
      subroutine MQC_MatrixFile_Write_Integer_Array(iUnit,label,nInteger,  &
        bufferLength,n1,n2,n3,n4,  &
        n5,asymmetric,array)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::iUnit,nInteger,bufferLength,n1,n2,  &
        n3,n4,n5
      logical,intent(in)::asymmetric
      integer(kind=int64),dimension(:),intent(in)::array
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Write_Integer_Array
!
!PROCEDURE MQC_MatrixFile_Write_Complex_Array
      subroutine MQC_MatrixFile_Write_Complex_Array(iUnit,label,nReal,  &
        bufferLength,n1,n2,n3,n4,n5,  &
        asymmetric,array)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::iUnit,nReal,bufferLength,n1,n2,n3,  &
        n4,n5
      logical,intent(in)::asymmetric
      complex(kind=real64),dimension(:),intent(inout)::array
!
      call MQC_MatrixFile_Unsupported()
      end subroutine MQC_MatrixFile_Write_Complex_Array
!
      end module MQC_MatWrapper
