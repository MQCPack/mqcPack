      program unitTest06
!
!     A FormChk-only build must reject MatrixFile operations with the stable
!     unsupported-feature diagnostic checked by the surrounding shell test.
!
!     H. P. Hratchian, 2026.
!
      use MQC_Gaussian
      use MQC_General, only: mqc_error
      implicit none
!
      type(MQC_Gaussian_Unformatted_Matrix_File)::matrixFile
!
      call matrixFile%load('unitTest06.faf')
      call mqc_error('unitTest06: MatrixFile load unexpectedly returned.')
!
      end program unitTest06
