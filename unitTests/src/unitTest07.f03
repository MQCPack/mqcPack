      program unitTest07
!
!     Validate numerical MatrixFile writing and reading through the public
!     MQCPack interface shared by the public and frontier GauOpen profiles.
!     Reading the matrix first also exercises numerical record skipping and
!     rewind before the earlier vector is requested.
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env, only: int64,real64
      use mqc_general
      use mqc_algebra2
      use mqc_gaussian
      implicit none
!
      character(len=*),parameter::fileName='unitTest07-numerical.faf'
      real(kind=real64),dimension(4),parameter::vectorExpected =  &
        [1.25_real64,-2.5_real64,3.75_real64,-4.125_real64]
      real(kind=real64),dimension(2,3),parameter::matrixExpected =  &
        reshape([1.0_real64,2.0_real64,3.0_real64,4.0_real64,  &
          5.0_real64,6.0_real64],[2,3])
!
      type(MQC_Gaussian_Unformatted_Matrix_File)::outputFile,inputFile
      type(MQC_Variable)::writeValue,readValue
      real(kind=real64),dimension(:),allocatable::vectorActual
      real(kind=real64),dimension(:,:),allocatable::matrixActual
!
      call deleteFile(fileName)
!
      writeValue = vectorExpected
      call outputFile%writeArray2('MQC TEST REAL VECTOR',writeValue,  &
        filename=fileName)
      writeValue = matrixExpected
      call outputFile%writeArray2('MQC TEST REAL MATRIX',writeValue)
      call outputFile%closeFile()
!
!     Read the last record first, then rewind for the earlier record.
!
      call inputFile%load(fileName)
      call inputFile%getArray('MQC TEST REAL MATRIX',mqcVarOut=readValue)
      if(RANK(readValue).ne.2)  &
        call mqc_error('unitTest07: expected a rank-2 real matrix.')
      if(SIZE(readValue,1).ne.SIZE(matrixExpected,1).or.  &
        SIZE(readValue,2).ne.SIZE(matrixExpected,2))  &
        call mqc_error('unitTest07: unexpected real-matrix shape.')
      matrixActual = readValue
      if(ANY(ABS(matrixActual-matrixExpected).gt.mqc_small))  &
        call mqc_error('unitTest07: real-matrix round trip failed.')
!
      call inputFile%getArray('MQC TEST REAL VECTOR',mqcVarOut=readValue)
      if(RANK(readValue).ne.1)  &
        call mqc_error('unitTest07: expected a rank-1 real vector.')
      if(SIZE(readValue).ne.SIZE(vectorExpected))  &
        call mqc_error('unitTest07: unexpected real-vector size.')
      vectorActual = readValue
      if(ANY(ABS(vectorActual-vectorExpected).gt.mqc_small))  &
        call mqc_error('unitTest07: real-vector round trip failed.')
      call inputFile%closeFile()
!
      call deleteFile(fileName)
      write(*,'(1x,A)') 'unitTest07: PASS'
!
      contains


!PROCEDURE deleteFile
      subroutine deleteFile(path)
!
!     Delete a prior test MatrixFile without depending on shell cleanup.
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::path
      integer(kind=int64)::deleteUnit,ioStatus
      logical::exists
!
      inquire(file=path,exist=exists)
      if(exists) then
        open(newunit=deleteUnit,file=path,status='old',  &
          form='unformatted',access='sequential',iostat=ioStatus)
        if(ioStatus.ne.0)  &
          call mqc_error('unitTest07: unable to open old test file.')
        close(deleteUnit,status='delete',iostat=ioStatus)
        if(ioStatus.ne.0)  &
          call mqc_error('unitTest07: unable to delete old test file.')
      endIf
!
      return
      end subroutine deleteFile


      end program unitTest07
