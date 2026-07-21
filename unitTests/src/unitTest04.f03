      program unitTest04
!
!     Validate disconnected deep-copy assignment for Gaussian MatrixFile
!     objects, assignment over an open destination, usable copied header
!     metadata, preservation of the source connection, and repeated access
!     through the Gaussian-scalar cache. The test creates both FAFs itself.
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env, only: int64,real64
      use mqc_general
      use mqc_algebra2
      use mqc_gaussian
      implicit none
!
      character(len=*),parameter::sourceFileName =  &
        'unitTest04-assignment-source.faf'
      character(len=*),parameter::outputFileName =  &
        'unitTest04-assignment-output.faf'
      character(len=8),parameter::titleExpected='test'
      character(len=72),parameter::routeExpected =  &
        '#p test geom=modela b3lyp sto-3g nosymm int=dkh output=faf '//  &
        'iop(3/33=2)'
      character(len=16),parameter::pointGroupExpected='C02V'
      character(len=8),dimension(4),parameter::irrepExpected =  &
        [character(len=8)::'A1','A2','B1','B2']
      real(kind=real64),parameter::scfEnergyExpected=-75.123456_real64
!
      type(MQC_Gaussian_Unformatted_Matrix_File)::createFile,sourceFile
      type(MQC_Gaussian_Unformatted_Matrix_File)::outputFile,verifyFile
      type(MQC_Variable)::characterValue,gaussianScalarsValue
      character(len=:),dimension(:),allocatable::irrepActual
      real(kind=real64),dimension(50)::gaussianScalarsExpected
      real(kind=real64)::scfEnergyFirst,scfEnergySecond,cachedEnergy
      integer(kind=int64)::sourceUnit,outputUnit
      logical::unitOpen
!
      call deleteFile(sourceFileName)
      call deleteFile(outputFileName)
!
!     Build a small source file containing both character records and enough
!     Gaussian scalar values for the public SCF-energy cache pathway.
!
      characterValue = titleExpected
      call createFile%writeArray2('TITLE',characterValue,  &
        filename=sourceFileName)
      characterValue = routeExpected
      call createFile%writeArray2('ROUTE',characterValue)
      characterValue = pointGroupExpected
      call createFile%writeArray2('POINT GROUP',characterValue)
      characterValue = irrepExpected
      call createFile%writeArray2('POINT GROUP IRREP NAMES',characterValue)
      gaussianScalarsExpected = 0.0_real64
      gaussianScalarsExpected(32) = scfEnergyExpected
      gaussianScalarsValue = gaussianScalarsExpected
      call createFile%writeArray2('GAUSSIAN SCALARS',gaussianScalarsValue)
      call createFile%closeFile()
!
      call sourceFile%load(sourceFileName)
      scfEnergyFirst = sourceFile%getValReal('scfenergy')
      scfEnergySecond = sourceFile%getValReal('scfenergy')
      call assertReal('first SCF energy',scfEnergyFirst,scfEnergyExpected)
      call assertReal('cached SCF energy',scfEnergySecond,scfEnergyExpected)
      call assertTrue(sourceFile%gaussianScalars_read,  &
        'Gaussian scalar retrieval did not set the cache flag.')
      sourceUnit = sourceFile%UnitNumber
      outputFile = sourceFile
!
!     Assignment must preserve the source connection and disconnect the copy.
!
      call assertTrue(sourceFile%isOpen(),  &
        'Assignment closed the source FAF.')
      call assertTrue(sourceFile%UnitNumber.eq.sourceUnit,  &
        'Assignment changed the source unit number.')
      call assertDisconnected(outputFile)
!
!     Check fixed metadata and allocation/value state.
!
      call assertTrue(outputFile%gaussianScalars_read.eqv.  &
        sourceFile%gaussianScalars_read,  &
        'gaussianScalars_read was not copied.')
      call assertTrue(outputFile%LabFil.eq.sourceFile%LabFil,  &
        'LabFil was not copied.')
      call assertTrue(outputFile%GVers.eq.sourceFile%GVers,  &
        'GVers was not copied.')
      call assertTrue(outputFile%Title.eq.sourceFile%Title,  &
        'Title was not copied.')
      call assertTrue(outputFile%matrixFileLayout%integerBytes.eq.  &
        sourceFile%matrixFileLayout%integerBytes,  &
        'Matrix-file integer width was not copied.')
      call assertTrue(outputFile%matrixFileLayout%fileVersion.eq.  &
        sourceFile%matrixFileLayout%fileVersion,  &
        'Matrix-file version was not copied.')
      call assertTrue(outputFile%matrixFileLayout%nHeaderRecords.eq.  &
        sourceFile%matrixFileLayout%nHeaderRecords,  &
        'Matrix-file header count was not copied.')
      call assertTrue(outputFile%matrixFileLayout%len12Bytes.eq.  &
        sourceFile%matrixFileLayout%len12Bytes,  &
        'Matrix-file Len12L size was not copied.')
      call assertTrue(outputFile%matrixFileLayout%len4Bytes.eq.  &
        sourceFile%matrixFileLayout%len4Bytes,  &
        'Matrix-file Len4L size was not copied.')
      call assertTrue(outputFile%matrixFileLayout%rawFormat.eqv.  &
        sourceFile%matrixFileLayout%rawFormat,  &
        'Matrix-file raw-format flag was not copied.')
!
      call assertIntegerScalar('natoms',outputFile%natoms,sourceFile%natoms)
      call assertIntegerScalar('nbasis',outputFile%nbasis,sourceFile%nbasis)
      call assertIntegerScalar('nbasisUse',outputFile%nbasisUse,  &
        sourceFile%nbasisUse)
      call assertIntegerScalar('icharge',outputFile%icharge,  &
        sourceFile%icharge)
      call assertIntegerScalar('multiplicity',outputFile%multiplicity,  &
        sourceFile%multiplicity)
      call assertIntegerScalar('nelectrons',outputFile%nelectrons,  &
        sourceFile%nelectrons)
      call assertIntegerScalar('icgu',outputFile%icgu,sourceFile%icgu)
      call assertIntegerScalar('NFC',outputFile%NFC,sourceFile%NFC)
      call assertIntegerScalar('NFV',outputFile%NFV,sourceFile%NFV)
      call assertIntegerScalar('ITran',outputFile%ITran,sourceFile%ITran)
      call assertIntegerScalar('IDum9',outputFile%IDum9,sourceFile%IDum9)
      call assertIntegerScalar('NShlAO',outputFile%NShlAO,  &
        sourceFile%NShlAO)
      call assertIntegerScalar('NPrmAO',outputFile%NPrmAO,  &
        sourceFile%NPrmAO)
      call assertIntegerScalar('NShlDB',outputFile%NShlDB,  &
        sourceFile%NShlDB)
      call assertIntegerScalar('NPrmDB',outputFile%NPrmDB,  &
        sourceFile%NPrmDB)
      call assertIntegerScalar('NBTot',outputFile%NBTot,sourceFile%NBTot)
      call assertIntegerVector('atomicNumbers',outputFile%atomicNumbers,  &
        sourceFile%atomicNumbers)
      call assertIntegerVector('atomTypes',outputFile%atomTypes,  &
        sourceFile%atomTypes)
      call assertIntegerVector('basisFunction2Atom',  &
        outputFile%basisFunction2Atom,sourceFile%basisFunction2Atom)
      call assertIntegerVector('IBasisFunctionType',  &
        outputFile%IBasisFunctionType,sourceFile%IBasisFunctionType)
      call assertIntegerVector('IgaussianScalars',  &
        outputFile%IgaussianScalars,sourceFile%IgaussianScalars)
      call assertRealVector('atomicCharges',outputFile%atomicCharges,  &
        sourceFile%atomicCharges)
      call assertRealVector('atomicWeights',outputFile%atomicWeights,  &
        sourceFile%atomicWeights)
      call assertRealVector('cartesians',outputFile%cartesians,  &
        sourceFile%cartesians)
      call assertRealVector('gaussianScalars',outputFile%gaussianScalars,  &
        sourceFile%gaussianScalars)
!
!     Allocatable assignment must produce independent cache storage.
!
      call assertTrue(allocated(outputFile%gaussianScalars),  &
        'The copied Gaussian scalar cache is not allocated.')
      cachedEnergy = sourceFile%gaussianScalars(32)
      outputFile%gaussianScalars(32) = cachedEnergy+1.0_real64
      call assertReal('independent Gaussian scalar cache',  &
        sourceFile%gaussianScalars(32),cachedEnergy)
      outputFile%gaussianScalars(32) = cachedEnergy
!
!     Copied header metadata must be sufficient to create a new FAF.
!
      call outputFile%create(outputFileName)
      call sourceFile%getArray('TITLE',mqcVarOut=characterValue)
      call outputFile%writeArray2('COPIED TITLE',characterValue)
      call sourceFile%getArray('ROUTE',mqcVarOut=characterValue)
      call outputFile%writeArray2('COPIED ROUTE',characterValue)
      call sourceFile%getArray('POINT GROUP',mqcVarOut=characterValue)
      call outputFile%writeArray2('COPIED POINT GROUP',characterValue)
      call sourceFile%getArray('POINT GROUP IRREP NAMES',  &
        mqcVarOut=characterValue)
      call outputFile%writeArray2('COPIED IRREPS',characterValue)
!
!     Replacing an open destination must close its old unit without changing
!     the source connection, then leave the replacement disconnected.
!
      outputUnit = outputFile%UnitNumber
      outputFile = sourceFile
      call assertDisconnected(outputFile)
      inquire(unit=outputUnit,opened=unitOpen)
      call assertTrue(.not.unitOpen,  &
        'Assignment did not close the previous output unit.')
      call assertTrue(sourceFile%isOpen(),  &
        'Replacing an open destination closed the source FAF.')
      call outputFile%closeFile()
      call sourceFile%closeFile()
!
      call verifyFile%load(outputFileName)
      call verifyFile%getArray('COPIED IRREPS',mqcVarOut=characterValue)
      call assertTrue(RANK(characterValue).eq.1,  &
        'The copied irrep record is not a vector.')
      call assertTrue(SIZE(characterValue).eq.SIZE(irrepExpected),  &
        'The copied irrep vector has the wrong size.')
      call assertTrue(LEN(characterValue).eq.LEN(irrepExpected),  &
        'The copied irrep vector has the wrong width.')
      irrepActual = characterValue
      call assertTrue(ALL(irrepActual.eq.irrepExpected),  &
        'The copied irrep values are incorrect.')
      call verifyFile%closeFile()
!
      call deleteFile(sourceFileName)
      call deleteFile(outputFileName)
      write(*,'(1x,A)') 'unitTest04: PASS'
!
      contains


!PROCEDURE assertTrue
      subroutine assertTrue(condition,message)
      implicit none
      logical,intent(in)::condition
      character(len=*),intent(in)::message
!
      if(.not.condition) call mqc_error(TRIM(message))
!
      return
      end subroutine assertTrue


!PROCEDURE assertReal
      subroutine assertReal(label,value,expected)
      implicit none
      character(len=*),intent(in)::label
      real(kind=real64),intent(in)::value,expected
!
      if(ABS(value-expected).gt.mqc_small)  &
        call mqc_error(TRIM(label)//' differs.')
!
      return
      end subroutine assertReal


!PROCEDURE assertDisconnected
      subroutine assertDisconnected(fileinfo)
      implicit none
      type(MQC_Gaussian_Unformatted_Matrix_File),intent(in)::fileinfo
!
      call assertTrue(.not.fileinfo%isOpen(),  &
        'The assigned FAF object is marked open.')
      call assertTrue(fileinfo%UnitNumber.eq.0,  &
        'The assigned FAF object has a nonzero unit number.')
      call assertTrue(LEN_TRIM(fileinfo%filename).eq.0,  &
        'The assigned FAF object retained a filename.')
      call assertTrue(.not.fileinfo%declared,  &
        'The assigned FAF object retained its declared flag.')
      call assertTrue(.not.fileinfo%header_read,  &
        'The assigned FAF object retained its header-read flag.')
      call assertTrue(.not.fileinfo%header_written,  &
        'The assigned FAF object retained its header-written flag.')
      call assertTrue(fileinfo%readWriteMode.eq.' ',  &
        'The assigned FAF object retained its read/write mode.')
!
      return
      end subroutine assertDisconnected


!PROCEDURE assertIntegerScalar
      subroutine assertIntegerScalar(label,value,expected)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),allocatable,intent(in)::value,expected
!
      call assertTrue(allocated(value).eqv.allocated(expected),  &
        TRIM(label)//' allocation state differs.')
      if(allocated(expected))  &
        call assertTrue(value.eq.expected,TRIM(label)//' value differs.')
!
      return
      end subroutine assertIntegerScalar


!PROCEDURE assertIntegerVector
      subroutine assertIntegerVector(label,value,expected)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),dimension(:),allocatable,intent(in)::value,expected
!
      call assertTrue(allocated(value).eqv.allocated(expected),  &
        TRIM(label)//' allocation state differs.')
      if(allocated(expected)) then
        call assertTrue(SIZE(value).eq.SIZE(expected),  &
          TRIM(label)//' size differs.')
        call assertTrue(ALL(value.eq.expected),TRIM(label)//' values differ.')
      endIf
!
      return
      end subroutine assertIntegerVector


!PROCEDURE assertRealVector
      subroutine assertRealVector(label,value,expected)
      implicit none
      character(len=*),intent(in)::label
      real(kind=real64),dimension(:),allocatable,intent(in)::value,expected
!
      call assertTrue(allocated(value).eqv.allocated(expected),  &
        TRIM(label)//' allocation state differs.')
      if(allocated(expected)) then
        call assertTrue(SIZE(value).eq.SIZE(expected),  &
          TRIM(label)//' size differs.')
        call assertTrue(ALL(ABS(value-expected).le.mqc_small),  &
          TRIM(label)//' values differ.')
      endIf
!
      return
      end subroutine assertRealVector


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


      end program unitTest04
