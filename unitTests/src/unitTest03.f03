      program unitTest03
!
!     Validate character scalar and fixed-width-vector FAF records through
!     the public MQCPack interface. The final raw record deliberately uses
!     Gaussian's legacy TypeA=1 scalar marker so the test also protects
!     normalization, character-aware skipping, rewind, rank, width, and value
!     behavior at the MatrixFile boundary.
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env, only: int64
      use mqc_general
      use mqc_algebra2
      use MQC_MatWrapper, only:  &
        Wr_Labl_TypeA => MQC_MatrixFile_Write_Label,  &
        Wr_ChBuf => MQC_MatrixFile_Write_Character_Buffer
      use mqc_gaussian
      implicit none
!
      character(len=*),parameter::fileName='unitTest03-character.faf'
      character(len=8),parameter::titleExpected='test'
      character(len=8),parameter::legacyExpected='legacy'
      character(len=72),parameter::routeExpected =  &
        '#p test geom=modela b3lyp sto-3g nosymm int=dkh output=faf '//  &
        'iop(3/33=2)'
      character(len=16),parameter::pointGroupExpected='C02V'
      character(len=8),dimension(4),parameter::irrepExpected =  &
        [character(len=8)::'A1','A2','B1','B2']
      integer(kind=int64),parameter::lenBufWords=4000_int64
!
      type(MQC_Gaussian_Unformatted_Matrix_File)::outputFile,inputFile
      type(MQC_Variable)::titleValue,routeValue,pointGroupValue,irrepValue
      type(MQC_Variable)::readValue
!
      call deleteFile(fileName)
!
      titleValue = titleExpected
      routeValue = routeExpected
      pointGroupValue = pointGroupExpected
      irrepValue = irrepExpected
!
      call outputFile%writeArray2('TITLE',titleValue,filename=fileName)
      call outputFile%writeArray2('ROUTE',routeValue)
      call outputFile%writeArray2('POINT GROUP',pointGroupValue)
      call outputFile%writeArray2('POINT GROUP IRREP NAMES',irrepValue)
!
!     Append one width-8 scalar using raw TypeA=1. NTot and LenBuf are file
!     words, while Wr_ChBuf receives the corresponding character counts.
!
      call Wr_Labl_TypeA(outputFile%unitNumber,'LEGACY TITLE',1_int64,  &
        0_int64,1_int64,lenBufWords,1_int64,1_int64,1_int64,1_int64,  &
        1_int64,1_int64)
      call Wr_ChBuf(outputFile%unitNumber,8_int64,  &
        outputFile%matrixFileLayout%len4Bytes*lenBufWords,legacyExpected)
      call outputFile%closeFile()
!
!     Reading the final record first forces the scanner to consume every
!     preceding character payload. Later requests require rewind and rescan.
!
      call inputFile%load(fileName)
      call inputFile%getArray('LEGACY TITLE',mqcVarOut=readValue)
      call validateScalar('LEGACY TITLE',readValue,legacyExpected,8_int64)
      call inputFile%getArray('POINT GROUP IRREP NAMES',mqcVarOut=readValue)
      call validateVector('POINT GROUP IRREP NAMES',readValue,irrepExpected)
      call inputFile%getArray('ROUTE',mqcVarOut=readValue)
      call validateScalar('ROUTE',readValue,routeExpected,72_int64)
      call inputFile%getArray('POINT GROUP',mqcVarOut=readValue)
      call validateScalar('POINT GROUP',readValue,pointGroupExpected,  &
        16_int64)
      call inputFile%getArray('TITLE',mqcVarOut=readValue)
      call validateScalar('TITLE',readValue,titleExpected,8_int64)
      call inputFile%closeFile()
!
      call deleteFile(fileName)
      write(*,'(1x,A)') 'unitTest03: PASS'
!
      contains


!PROCEDURE validateScalar
      subroutine validateScalar(label,value,expected,expectedWidth)
      implicit none
      character(len=*),intent(in)::label,expected
      type(MQC_Variable),intent(in)::value
      integer(kind=int64),intent(in)::expectedWidth
      character(len=:),allocatable::actual
!
      if(TRIM(value%getType()).ne.'CHARACTER')  &
        call mqc_error('Expected a character value for '//TRIM(label))
      if(RANK(value).ne.0)  &
        call mqc_error('Expected a scalar value for '//TRIM(label))
      if(SIZE(value).ne.1)  &
        call mqc_error('Expected one scalar element for '//TRIM(label))
      if(LEN(value).ne.expectedWidth)  &
        call mqc_error('Unexpected character width for '//TRIM(label))
      allocate(character(len=LEN(value))::actual)
      actual = value
      if(actual.ne.expected)  &
        call mqc_error('Unexpected character value for '//TRIM(label))
!
      return
      end subroutine validateScalar


!PROCEDURE validateVector
      subroutine validateVector(label,value,expected)
      implicit none
      character(len=*),intent(in)::label
      type(MQC_Variable),intent(in)::value
      character(len=*),dimension(:),intent(in)::expected
      character(len=:),dimension(:),allocatable::actual
!
      if(TRIM(value%getType()).ne.'CHARACTER')  &
        call mqc_error('Expected a character value for '//TRIM(label))
      if(RANK(value).ne.1)  &
        call mqc_error('Expected a vector value for '//TRIM(label))
      if(SIZE(value).ne.SIZE(expected))  &
        call mqc_error('Unexpected vector size for '//TRIM(label))
      if(LEN(value).ne.LEN(expected))  &
        call mqc_error('Unexpected vector width for '//TRIM(label))
      actual = value
      if(ANY(actual.ne.expected))  &
        call mqc_error('Unexpected vector values for '//TRIM(label))
!
      return
      end subroutine validateVector


!PROCEDURE deleteFile
      subroutine deleteFile(path)
      implicit none
      character(len=*),intent(in)::path
      integer(kind=int64)::deleteUnit,ioStatus
      logical::exists
!
      inquire(file=path,exist=exists)
      if(exists) then
        open(newunit=deleteUnit,file=path,status='old',form='unformatted',  &
          access='sequential',iostat=ioStatus)
        if(ioStatus.ne.0)  &
          call mqc_error('Unable to open temporary FAF for deletion.')
        close(deleteUnit,status='delete')
      endIf
!
      return
      end subroutine deleteFile


      end program unitTest03
