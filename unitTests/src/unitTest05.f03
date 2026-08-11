      program unitTest05
!
!     Confirm that a FormChk-only MQCPack build links MQC_Gaussian and reads
!     the formatted-checkpoint title and route summary without GauOpen.
!
!     H. P. Hratchian, 2026.
!
      use MQC_Gaussian
      use MQC_General, only: mqc_error
      implicit none
!
      type(MQC_Gaussian_FChk_File)::fchkFile
      character(len=*),parameter::fileName='unitTest05.fchk'
      integer::fileUnit
      logical::ok
!
      open(newunit=fileUnit,file=fileName,status='replace',action='write')
      write(fileUnit,'(A)') 'FormChk-only MQC_Gaussian smoke test'
      write(fileUnit,'(A)') 'SP RHF STO-3G'
      close(fileUnit)
!
      call fchkFile%openFile(fileName,-1,ok)
      if(.not.ok)  &
        call mqc_error('unitTest05: unable to open the test FChk file.')
      if(trim(fchkFile%title).ne.'FormChk-only MQC_Gaussian smoke test')  &
        call mqc_error('unitTest05: incorrect FChk title.')
      if(trim(fchkFile%jobType).ne.'SP')  &
        call mqc_error('unitTest05: incorrect FChk job type.')
      if(trim(fchkFile%method).ne.'RHF')  &
        call mqc_error('unitTest05: incorrect FChk method.')
      if(trim(fchkFile%basisSet).ne.'STO-3G')  &
        call mqc_error('unitTest05: incorrect FChk basis set.')
!
      call fchkFile%closeFile()
      open(newunit=fileUnit,file=fileName,status='old')
      close(fileUnit,status='delete')
!
      write(*,'(A)') 'PASS: FormChk-only MQC_Gaussian link and read'
!
      end program unitTest05
