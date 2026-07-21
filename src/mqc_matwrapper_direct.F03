#if defined(MQC_GAUOPEN_API_PUBLIC) && defined(MQC_GAUOPEN_API_FRONTIER)
#error "Only one GauOpen API family may be selected."
#elif !defined(MQC_GAUOPEN_API_PUBLIC) && !defined(MQC_GAUOPEN_API_FRONTIER)
#error "A GauOpen API family must be selected by configure."
#endif

      module MQC_MatWrapper
!
!     This module is the MQCPack-owned adapter for unmodified GauOpen
!     qcmatrix.F and qcmatrixio.F sources compiled with 8-byte integers. It
!     normalizes the public version-2 ASym API and the frontier TypeA API.
!     Raw GauOpen procedure names remain external and unmodified; prefixed
!     adapter names keep them separate from MQCPack's module symbols.
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env, only: int32, int64, real64
      use MQC_General, only: mqc_error,mqc_error_a,mqc_error_i
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
      character(len=*),parameter::publicCharacterError =  &
        'Character FAF I/O requires the GauOpen frontier TypeA API; '//  &
        'the configured public GauOpen API does not provide character '//  &
        'record support.'
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
!PROCEDURE MQC_MatrixFile_Layout_Set
      subroutine MQC_MatrixFile_Layout_Set(layout,fileVersion,  &
        nHeaderRecords,len12Bytes,len4Bytes)
!
!     Store layout metadata read from an 8-byte GauOpen MatrixFile header.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      class(MQC_MatrixFile_Layout),intent(inout)::layout
      integer(kind=int64),intent(in)::fileVersion,nHeaderRecords,  &
        len12Bytes,len4Bytes
!
      if(len12Bytes.ne.8_int64)  &
        call mqc_error_i('MatrixFile Len12L is incompatible with the '//  &
          'configured 8-byte GauOpen ABI.',6_int64,'len12Bytes',len12Bytes)
      if(len4Bytes.ne.8_int64)  &
        call mqc_error_i('MatrixFile Len4L is incompatible with the '//  &
          'configured 8-byte GauOpen ABI.',6_int64,'len4Bytes',len4Bytes)
#ifdef MQC_GAUOPEN_API_PUBLIC
      if(fileVersion.ge.3_int64)  &
        call mqc_error_i('The configured public GauOpen API cannot safely '//  &
          'identify character records in a frontier version-3 FAF.',  &
          6_int64,'fileVersion',fileVersion)
#endif
!
      layout%integerBytes = 8_int32
      layout%fileVersion = fileVersion
      layout%nHeaderRecords = nHeaderRecords
      layout%len12Bytes = len12Bytes
      layout%len4Bytes = len4Bytes
      layout%rawFormat = .false.
!
      return
      end subroutine MQC_MatrixFile_Layout_Set
!
!PROCEDURE MQC_MatrixFile_Layout_Set_For_Write
      subroutine MQC_MatrixFile_Layout_Set_For_Write(layout)
!
!     Initialize layout metadata for the configured 8-byte GauOpen ABI.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      class(MQC_MatrixFile_Layout),intent(inout)::layout
!
      layout%integerBytes = 8_int32
#ifdef MQC_GAUOPEN_API_PUBLIC
      layout%fileVersion = 2_int64
#else
      layout%fileVersion = 3_int64
#endif
      layout%nHeaderRecords = 11_int64
      layout%len12Bytes = 8_int64
      layout%len4Bytes = 8_int64
      layout%rawFormat = .false.
!
      return
      end subroutine MQC_MatrixFile_Layout_Set_For_Write
!
!PROCEDURE MQC_MatrixFile_Get_Character_Record_Layout
      function MQC_MatrixFile_Get_Character_Record_Layout(layout,nTot,  &
        lenBuf,n1,n2,n3,n4,n5,typeA) result(recordLayout)
!
!     Normalize GauOpen character-record words into fixed-width elements.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      type(MQC_MatrixFile_Layout),intent(in)::layout
      integer(kind=int64),intent(in)::nTot,lenBuf,n1,n2,n3,n4,n5,typeA
      type(MQC_MatrixFile_Character_Record_Layout)::recordLayout
      integer(kind=int64)::rawElementCount
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call mqc_error(publicCharacterError,6_int64)
      return
#else
      if(layout%len4Bytes.le.0)  &
        call mqc_error_i('Character FAF record has invalid Len4L.',  &
          6_int64,'layout%len4Bytes',layout%len4Bytes)
      if(nTot.le.0)  &
        call mqc_error_i('Character FAF record has invalid NTot.',  &
          6_int64,'nTot',nTot)
      if(lenBuf.le.0)  &
        call mqc_error_i('Character FAF record has invalid LenBuf.',  &
          6_int64,'lenBuf',lenBuf)
      if(typeA.le.0)  &
        call mqc_error_i('Character FAF record has invalid TypeA.',  &
          6_int64,'typeA',typeA)
!
      rawElementCount = MQC_MatrixFile_Len_Array(n1,n2,n3,n4,n5)
      if(rawElementCount.le.0)  &
        call mqc_error_i('Character FAF record has invalid raw dimensions.',  &
          6_int64,'rawElementCount',rawElementCount)
      recordLayout%payloadCharacters = layout%len4Bytes*nTot
      recordLayout%payloadWords = nTot
      if(typeA.eq.1) then
        recordLayout%elementWidth = recordLayout%payloadCharacters
      else
        recordLayout%elementWidth = typeA
      endIf
      if(mod(recordLayout%payloadCharacters,  &
        recordLayout%elementWidth).ne.0)  &
        call mqc_error('Character FAF payload length is not divisible by '//  &
          'the logical element width.',6_int64)
      recordLayout%elementCount = recordLayout%payloadCharacters/  &
        recordLayout%elementWidth
      recordLayout%bufferCharacters = layout%len4Bytes*lenBuf
      recordLayout%bufferWords = lenBuf
#endif
!
      return
      end function MQC_MatrixFile_Get_Character_Record_Layout
!
!PROCEDURE MQC_MatrixFile_Get_Character_Write_Layout
      function MQC_MatrixFile_Get_Character_Write_Layout(layout,  &
        elementCount,elementWidth,bufferWords) result(recordLayout)
!
!     Translate fixed-width character elements into GauOpen word counts.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      type(MQC_MatrixFile_Layout),intent(in)::layout
      integer(kind=int64),intent(in)::elementCount,elementWidth,bufferWords
      type(MQC_MatrixFile_Character_Record_Layout)::recordLayout
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call mqc_error(publicCharacterError,6_int64)
      return
#else
      if(layout%len4Bytes.le.0)  &
        call mqc_error_i('Character FAF write has invalid Len4L.',  &
          6_int64,'layout%len4Bytes',layout%len4Bytes)
      if(elementCount.le.0)  &
        call mqc_error_i('Character FAF write has no logical elements.',  &
          6_int64,'elementCount',elementCount)
      if(elementWidth.le.1)  &
        call mqc_error_i('Character FAF element width must be greater than '//  &
          'one; TypeA=1 is the legacy scalar marker.',6_int64,  &
          'elementWidth',elementWidth)
      if(bufferWords.le.0)  &
        call mqc_error_i('Character FAF write has invalid buffer words.',  &
          6_int64,'bufferWords',bufferWords)
!
      recordLayout%elementCount = elementCount
      recordLayout%elementWidth = elementWidth
      recordLayout%payloadCharacters = elementCount*elementWidth
      if(mod(recordLayout%payloadCharacters,layout%len4Bytes).ne.0)  &
        call mqc_error('Character FAF payload must occupy a whole number '//  &
          'of Len4L words.',6_int64)
      recordLayout%payloadWords = recordLayout%payloadCharacters/  &
        layout%len4Bytes
      recordLayout%bufferWords = bufferWords
      recordLayout%bufferCharacters = layout%len4Bytes*bufferWords
#endif
!
      return
      end function MQC_MatrixFile_Get_Character_Write_Layout
!
!PROCEDURE MQC_MatrixFile_Open_Read
      subroutine MQC_MatrixFile_Open_Read(name,iUnit,labFil,iVers,nLab,  &
        gVers,title,nAtoms,nBasis,nBasisUse,iCharge,multiplicity,  &
        nElectrons,len12L,len4L,iOpCl,icgu)
      implicit none
      character(len=*),intent(in)::name
      character(len=64),intent(inout)::labFil,gVers,title
      integer(kind=int64),intent(out)::iUnit,iVers,nLab,nAtoms,nBasis,  &
        nBasisUse,iCharge,multiplicity,nElectrons,len12L,len4L,iOpCl,icgu
!
      call Open_Read(name,iUnit,labFil,iVers,nLab,gVers,title,nAtoms,  &
        nBasis,nBasisUse,iCharge,multiplicity,nElectrons,len12L,len4L,  &
        iOpCl,icgu)
      if(iUnit.eq.-1_int64)  &
        call mqc_error_a('Unable to open MatrixFile for reading.',6_int64,  &
          'name',name)
!
      return
      end subroutine MQC_MatrixFile_Open_Read
!
!PROCEDURE MQC_MatrixFile_Open_Write
      subroutine MQC_MatrixFile_Open_Write(name,iUnit,labFil,gVers,title,  &
        nAtoms,nBasis,nBasisUse,iCharge,multiplicity,nElectrons,iOpCl,icgu)
      implicit none
      character(len=*),intent(in)::name
      character(len=64),intent(inout)::labFil,gVers,title
      integer(kind=int64),intent(out)::iUnit
      integer(kind=int64),intent(in)::nAtoms,nBasis,nBasisUse,iCharge,  &
        multiplicity,nElectrons,iOpCl,icgu
!
      call Open_Write(name,iUnit,labFil,gVers,title,nAtoms,nBasis,  &
        nBasisUse,iCharge,multiplicity,nElectrons,iOpCl,icgu)
      if(iUnit.eq.-1_int64)  &
        call mqc_error_a('Unable to open MatrixFile for writing.',6_int64,  &
          'name',name)
!
      return
      end subroutine MQC_MatrixFile_Open_Write
!
!PROCEDURE MQC_MatrixFile_Close
      subroutine MQC_MatrixFile_Close(iUnit)
      implicit none
      integer(kind=int64),intent(in)::iUnit
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call Close_MatF(iUnit)
#else
      call Close_FAF(.false.,iUnit)
#endif
      return
      end subroutine MQC_MatrixFile_Close
!
!PROCEDURE MQC_MatrixFile_Len_Array
      function MQC_MatrixFile_Len_Array(n1,n2,n3,n4,n5)  &
        result(arrayLength)
      implicit none
      integer(kind=int64),intent(in)::n1,n2,n3,n4,n5
      integer(kind=int64)::arrayLength
      integer(kind=int64),external::LenArr
!
      arrayLength = LenArr(n1,n2,n3,n4,n5)
      return
      end function MQC_MatrixFile_Len_Array
!
!PROCEDURE MQC_MatrixFile_Read_Header
      subroutine MQC_MatrixFile_Read_Header(iUnit,nLab,nAtoms,nBasis,  &
        atomicNumbers,atomTypes,atomicCharges,cartesians,basisFunction2Atom,  &
        basisFunctionTypes,atomicWeights,nfc,nfv,iTran,iDum9,nShellAO,  &
        nPrimitiveAO,nShellDB,nPrimitiveDB,nBasisTotal)
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
      call Rd_Head(iUnit,nLab,nAtoms,nBasis,atomicNumbers,atomTypes,  &
        atomicCharges,cartesians,basisFunction2Atom,basisFunctionTypes,  &
        atomicWeights,nfc,nfv,iTran,iDum9,nShellAO,nPrimitiveAO,nShellDB,  &
        nPrimitiveDB,nBasisTotal)
      return
      end subroutine MQC_MatrixFile_Read_Header
!
!PROCEDURE MQC_MatrixFile_Write_Header
      subroutine MQC_MatrixFile_Write_Header(iUnit,nAtoms,nAtomCoordinates,  &
        nBasis,atomicNumbers,atomTypes,atomicCharges,cartesians,  &
        basisFunction2Atom,basisFunctionTypes,atomicWeights,nfc,nfv,iTran,  &
        iDum9,nShellAO,nPrimitiveAO,nShellDB,nPrimitiveDB,nBasisTotal)
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
      call Wr_Head(iUnit,nAtoms,nAtomCoordinates,nBasis,atomicNumbers,  &
        atomTypes,atomicCharges,cartesians,basisFunction2Atom,  &
        basisFunctionTypes,atomicWeights,nfc,nfv,iTran,iDum9,nShellAO,  &
        nPrimitiveAO,nShellDB,nPrimitiveDB,nBasisTotal)
      return
      end subroutine MQC_MatrixFile_Write_Header
!
!PROCEDURE MQC_MatrixFile_Read_Label
      subroutine MQC_MatrixFile_Read_Label(iUnit,iVers,label,nInteger,  &
        nReal,nTotal,bufferLength,n1,n2,n3,n4,n5,typeA,nRealInteger,  &
        endOfFile)
      implicit none
      integer(kind=int64),intent(in)::iUnit,iVers
      character(len=64),intent(out)::label
      integer(kind=int64),intent(out)::nInteger,nReal,nTotal,bufferLength,  &
        n1,n2,n3,n4,n5,typeA,nRealInteger
      logical,intent(out)::endOfFile
#ifdef MQC_GAUOPEN_API_PUBLIC
      logical::asymmetric
#endif
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call Rd_Labl(iUnit,iVers,label,nInteger,nReal,nTotal,bufferLength,  &
        n1,n2,n3,n4,n5,asymmetric,nRealInteger,endOfFile)
      typeA = 0_int64
      if(asymmetric) typeA = -1_int64
#else
      call Rd_Labl(iUnit,iVers,label,nInteger,nReal,nTotal,bufferLength,  &
        n1,n2,n3,n4,n5,typeA,nRealInteger,endOfFile)
#endif
      return
      end subroutine MQC_MatrixFile_Read_Label
!
!PROCEDURE MQC_MatrixFile_Write_Label
      subroutine MQC_MatrixFile_Write_Label(iUnit,label,nInteger,nReal,  &
        nTotal,bufferLength,n1,n2,n3,n4,n5,typeA)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nInteger,nReal,nTotal,  &
        bufferLength,n1,n2,n3,n4,n5,typeA
      character(len=*),intent(in)::label
#ifdef MQC_GAUOPEN_API_PUBLIC
      logical::asymmetric
#endif
!
      if(nReal.lt.0)  &
        call mqc_error('MQC_MatrixFile_Write_Label: nReal must be '//  &
          'nonnegative.',6_int64)
#ifdef MQC_GAUOPEN_API_PUBLIC
      if(typeA.ne.0_int64.and.typeA.ne.-1_int64)  &
        call mqc_error(publicCharacterError,6_int64)
      asymmetric = typeA.eq.-1_int64
      call Wr_Labl(iUnit,label,nInteger,nReal,nTotal,bufferLength,  &
        n1,n2,n3,n4,n5,asymmetric)
#else
      call Wr_Labl(iUnit,label,nInteger,nReal,nTotal,bufferLength,  &
        n1,n2,n3,n4,n5,typeA)
#endif
      return
      end subroutine MQC_MatrixFile_Write_Label
!
!PROCEDURE MQC_MatrixFile_Read_Character_Buffer
      subroutine MQC_MatrixFile_Read_Character_Buffer(iUnit,nCharacters,  &
        bufferLength,characterArray)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nCharacters,bufferLength
      character(len=*),intent(inout)::characterArray
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call mqc_error(publicCharacterError,6_int64)
#else
      call Rd_ChBuf(iUnit,nCharacters,bufferLength,characterArray)
#endif
      return
      end subroutine MQC_MatrixFile_Read_Character_Buffer
!
!PROCEDURE MQC_MatrixFile_Write_Character_Buffer
      subroutine MQC_MatrixFile_Write_Character_Buffer(iUnit,nCharacters,  &
        bufferLength,characterArray)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nCharacters,bufferLength
      character(len=*),intent(in)::characterArray
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call mqc_error(publicCharacterError,6_int64)
#else
      call Wr_ChBuf(iUnit,nCharacters,bufferLength,characterArray)
#endif
      return
      end subroutine MQC_MatrixFile_Write_Character_Buffer
!
!PROCEDURE MQC_MatrixFile_Read_Integer_Buffer
      subroutine MQC_MatrixFile_Read_Integer_Buffer(iUnit,nTotal,  &
        bufferLength,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nTotal,bufferLength
      integer(kind=int64),dimension(nTotal),intent(out)::array
!
      call Rd_IBuf(iUnit,nTotal,bufferLength,array)
      return
      end subroutine MQC_MatrixFile_Read_Integer_Buffer
!
!PROCEDURE MQC_MatrixFile_Read_Real_Buffer
      subroutine MQC_MatrixFile_Read_Real_Buffer(iUnit,nTotal,  &
        bufferLength,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nTotal,bufferLength
      real(kind=real64),dimension(:),intent(out)::array
!
      call Rd_RBuf(iUnit,nTotal,bufferLength,array)
      return
      end subroutine MQC_MatrixFile_Read_Real_Buffer
!
!PROCEDURE MQC_MatrixFile_Read_Complex_Buffer
      subroutine MQC_MatrixFile_Read_Complex_Buffer(iUnit,nTotal,  &
        bufferLength,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nTotal,bufferLength
      complex(kind=real64),dimension(nTotal),intent(inout)::array
!
      call Rd_CBuf(iUnit,nTotal,bufferLength,array)
      return
      end subroutine MQC_MatrixFile_Read_Complex_Buffer
!
!PROCEDURE MQC_MatrixFile_Read_Two_Electron
      subroutine MQC_MatrixFile_Read_Two_Electron(iUnit,nReal,recordLength,  &
        recordBlockLength,nTotal,bufferLength,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nReal,recordLength,  &
        recordBlockLength,nTotal,bufferLength
      real(kind=real64),dimension(:),intent(inout)::array
!
      call Rd_2EN(iUnit,nReal,recordLength,recordBlockLength,nTotal,  &
        bufferLength,array)
      return
      end subroutine MQC_MatrixFile_Read_Two_Electron
!
!PROCEDURE MQC_MatrixFile_Read_Indexed_Real
      subroutine MQC_MatrixFile_Read_Indexed_Real(iUnit,nReal,recordLength,  &
        nTotal,bufferLength,nNonzero,array)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nReal,recordLength,nTotal,  &
        bufferLength
      integer(kind=int64),intent(out)::nNonzero
      real(kind=real64),dimension(:),intent(out)::array
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call Rd_RInd(iUnit,nReal,recordLength,nTotal,bufferLength,nNonzero,  &
        array)
#else
      call Rd_RInd(iUnit,nReal,recordLength,nReal*recordLength,nTotal,  &
        bufferLength,array)
      nNonzero = nTotal
#endif
      return
      end subroutine MQC_MatrixFile_Read_Indexed_Real
!
!PROCEDURE MQC_MatrixFile_Skip
      subroutine MQC_MatrixFile_Skip(iUnit,nTotal,bufferLength)
      implicit none
      integer(kind=int64),intent(in)::iUnit,nTotal,bufferLength
!
      call Rd_Skip(iUnit,nTotal,bufferLength)
      return
      end subroutine MQC_MatrixFile_Skip
!
!PROCEDURE MQC_MatrixFile_Write_AO_Two_Electron
      subroutine MQC_MatrixFile_Write_AO_Two_Electron(iUnit,label,nReal,  &
        bufferLength,nIntegrals,integrals)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::iUnit,nReal,bufferLength,nIntegrals
      real(kind=real64),dimension(:),intent(inout)::integrals
!
      call Wr_LAO2E(iUnit,label,nReal,bufferLength,nIntegrals,integrals)
      return
      end subroutine MQC_MatrixFile_Write_AO_Two_Electron
!
!PROCEDURE MQC_MatrixFile_Write_Real_Array
      subroutine MQC_MatrixFile_Write_Real_Array(iUnit,label,nReal,  &
        bufferLength,n1,n2,n3,n4,n5,asymmetric,array)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::iUnit,nReal,bufferLength,n1,n2,n3,  &
        n4,n5
      logical,intent(in)::asymmetric
      real(kind=real64),dimension(:),intent(inout)::array
#ifdef MQC_GAUOPEN_API_FRONTIER
      integer(kind=int64)::typeA
#endif
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call Wr_LRBuf(iUnit,label,nReal,bufferLength,n1,n2,n3,n4,n5,  &
        asymmetric,array)
#else
      typeA = 0_int64
      if(asymmetric) typeA = -1_int64
      call Wr_LRBuf(iUnit,label,nReal,bufferLength,n1,n2,n3,n4,n5,typeA,  &
        array)
#endif
      return
      end subroutine MQC_MatrixFile_Write_Real_Array
!
!PROCEDURE MQC_MatrixFile_Write_Integer_Array
      subroutine MQC_MatrixFile_Write_Integer_Array(iUnit,label,nInteger,  &
        bufferLength,n1,n2,n3,n4,n5,asymmetric,array)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::iUnit,nInteger,bufferLength,n1,n2,  &
        n3,n4,n5
      logical,intent(in)::asymmetric
      integer(kind=int64),dimension(:),intent(in)::array
#ifdef MQC_GAUOPEN_API_FRONTIER
      integer(kind=int64)::typeA
#endif
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call Wr_LIBuf(iUnit,label,nInteger,bufferLength,n1,n2,n3,n4,n5,  &
        asymmetric,array)
#else
      typeA = 0_int64
      if(asymmetric) typeA = -1_int64
      call Wr_LIBuf(iUnit,label,nInteger,bufferLength,n1,n2,n3,n4,n5,  &
        typeA,array)
#endif
      return
      end subroutine MQC_MatrixFile_Write_Integer_Array
!
!PROCEDURE MQC_MatrixFile_Write_Complex_Array
      subroutine MQC_MatrixFile_Write_Complex_Array(iUnit,label,nReal,  &
        bufferLength,n1,n2,n3,n4,n5,asymmetric,array)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::iUnit,nReal,bufferLength,n1,n2,n3,  &
        n4,n5
      logical,intent(in)::asymmetric
      complex(kind=real64),dimension(:),intent(inout)::array
#ifdef MQC_GAUOPEN_API_FRONTIER
      integer(kind=int64)::typeA
#endif
!
#ifdef MQC_GAUOPEN_API_PUBLIC
      call Wr_LCBuf(iUnit,label,nReal,bufferLength,n1,n2,n3,n4,n5,  &
        asymmetric,array)
#else
      typeA = 0_int64
      if(asymmetric) typeA = -1_int64
      call Wr_LCBuf(iUnit,label,ABS(nReal),bufferLength,n1,n2,n3,n4,n5,  &
        typeA,array)
#endif
      return
      end subroutine MQC_MatrixFile_Write_Complex_Array
!
      end module MQC_MatWrapper
