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
!
!
!     This include file provides a set of routines/functions that serve as
!     wrappers for LAPACK routines.
!
!
!
!----------------------------------------------------------------
!                                                               |
!     Eigensystem Solvers                                       |
!                                                               |
!----------------------------------------------------------------
!
!PROCEDURE mqc_dsyev_eigensystem_symmFull
      subroutine mqc_dsyev_eigensystem_symmFull(matrix,eigenValues,eigenVectors)
!
!     This subroutine wraps the LAPACK eigensystem routines for a symmetric
!     matrix in full/dense memory storage format.
!
!     H. P. Hratchian, 2020
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::matrix
      real(kind=real64),dimension(:),allocatable,intent(out)::eigenValues
      real(kind=real64),dimension(:,:),allocatable,OPTIONAL,intent(out)::eigenVectors
!
      integer::i,N,INFO,LWORK
      real(kind=real64),dimension(:,:),allocatable::matrixCopy
      real(kind=real64),dimension(:),allocatable::Work
      real(kind=real64),dimension(1)::WorkTmp
      character(len=1)::JOBZ
!
!
!     Do some initial checking of matrix to ensure it's square. Next, figure out
!     the dimension of <matrix>, <N>. Then, allocate space for <matrixCopy> and
!     copy matrix into that space so that the LAPACK routine DSYEV can use this
!     space as working memory. Finally, allocate space for <eigenValues> and, if
!     sent, <eigenVectors>.
!
      if(SIZE(matrix,1).ne.SIZE(matrix,2)) call mqc_error(  &
        'mqc_eigensystem_symmMatrixFull: Rectangular matrix sent.')
      N = SIZE(matrix,1)
      Allocate(matrixCopy(N,N))
      matrixCopy = matrix
      if(Allocated(eigenValues)) then
        if(SIZE(eigenValues).ne.N) then
          DeAllocate(eigenValues)
          Allocate(eigenValues(N))
        endIf
      else
        Allocate(eigenValues(N))
      endIf
      if(PRESENT(eigenVectors)) then
        JOBZ = 'V'
        if(Allocated(eigenVectors)) then
          if(SIZE(eigenVectors,1).ne.N.or.SIZE(eigenVectors,2).ne.N) then
            DeAllocate(eigenVectors)
            Allocate(eigenVectors(N,N))
          endIf
        else
          Allocate(eigenVectors(N,N))
        endIf
      else
        JOBZ = 'N'
      endIf
!
!     Call LAPCK routine DSYEV. 
!
      Info = 0
      call dsyev(JOBZ,'U',N,matrixCopy,N,eigenValues,WorkTmp,-1,INFO)    
      if(INFO.ne.0) call mqc_error(  &
        'mqc_eigensystem_symmMatrixFull: Failure in DSYEV getting WORK length.')
      LWork = WorkTmp(1)
      Allocate(Work(LWork))
      call dsyev(JOBZ,'U',N,matrixCopy,N,eigenValues,Work,LWork,INFO)    
      if(INFO.ne.0) call mqc_error(  &
        'mqc_eigensystem_symmMatrixFull: Failure in DSYEV diagonalizing matrix.')
      if(PRESENT(eigenVectors)) then
        eigenVectors = matrixCopy
        do i = 1,N
          call mqc_vectorPhase(eigenVectors(:,i),.true.)
        endDo
      endIf
      deAllocate(matrixCopy,Work)
!
      return
      end subroutine mqc_dsyev_eigensystem_symmFull


!
!PROCEDURE mqc_dgeev_eigensystem_Full
      subroutine mqc_dgeev_eigensystem_Full(matrix,eigenValues,  &
        leftEigenVectors,rightEigenVectors)
!
!     This subroutine wraps the LAPACK eigensystem routines for a general square
!     matrix. 
!
!     H. P. Hratchian, 2022
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::matrix
      real(kind=real64),dimension(:),allocatable,intent(out)::eigenValues
      real(kind=real64),dimension(:,:),allocatable,OPTIONAL,intent(out)::leftEigenVectors,  &
        rightEigenVectors
!
      integer::i,N,LDVL,LDVR,INFO,LWORK
      real(kind=real64),dimension(:,:),allocatable::matrixCopy
      real(kind=real64),dimension(:),allocatable::eigenValuesImaginary
      real(kind=real64),dimension(:),allocatable::Work
      real(kind=real64),dimension(1)::WorkTmp
      character(len=1)::JOBVL,JOBVR
!
!
!     Do some initial checking of matrix to ensure it's square. Next, figure out
!     the dimension of <matrix>, <N>. Then, allocate space for <matrixCopy> and
!     copy matrix into that space so that the LAPACK routine DSYEV can use this
!     space as working memory. Finally, allocate space for <eigenValues> and, if
!     sent, <eigenVectors>.
!
      if(SIZE(matrix,1).ne.SIZE(matrix,2)) call mqc_error(  &
        'mqc_eigensystem_matrixFull: Rectangular matrix sent.')
      N = SIZE(matrix,1)
      Allocate(matrixCopy(N,N))
      matrixCopy = matrix
      if(Allocated(eigenValues)) then
        if(SIZE(eigenValues).ne.N) then
          DeAllocate(eigenValues)
          Allocate(eigenValues(N))
        endIf
      else
        Allocate(eigenValues(N))
      endIf
      Allocate(eigenValuesImaginary(N))
      if(PRESENT(leftEigenVectors)) then
        JOBVL = 'V'
        LDVL = N
        if(Allocated(leftEigenVectors)) then
          if(SIZE(leftEigenVectors,1).ne.N.or.SIZE(leftEigenVectors,2).ne.N) then
            DeAllocate(leftEigenVectors)
            Allocate(leftEigenVectors(N,N))
          endIf
        else
          Allocate(leftEigenVectors(N,N))
        endIf
      else
        JOBVL = 'N'
        LDVL = 1
      endIf
      if(PRESENT(rightEigenVectors)) then
        JOBVR = 'V'
        LDVR = N
        if(Allocated(rightEigenVectors)) then
          if(SIZE(rightEigenVectors,1).ne.N.or.SIZE(rightEigenVectors,2).ne.N) then
            DeAllocate(rightEigenVectors)
            Allocate(rightEigenVectors(N,N))
          endIf
        else
          Allocate(rightEigenVectors(N,N))
        endIf
      else
        JOBVR = 'N'
        LDVR = 1
      endIf
!
!     Call LAPCK routine DSYEV. 
!
      Info = 0
      call dgeev(JOBVL,JOBVR,N,matrixCopy,N,eigenValues,  &
        eigenValuesImaginary,leftEigenVectors,LDVL,rightEigenVectors,LDVR,  &
        WorkTmp,-1,INFO)
      if(INFO.ne.0) call mqc_error(  &
        'mqc_eigensystem_MatrixFull: Failure in DGEEV getting WORK length.')
      LWork = WorkTmp(1)
      Allocate(Work(LWork))
      call dgeev(JOBVL,JOBVR,N,matrixCopy,N,eigenValues,  &
        eigenValuesImaginary,leftEigenVectors,LDVL,rightEigenVectors,LDVR,  &
        Work,LWork,INFO)
      if(INFO.ne.0) call mqc_error(  &
        'mqc_eigensystem_MatrixFull: Failure in DGEEV diagonalizing matrix.')
      if(PRESENT(leftEigenVectors)) then
        do i = 1,N
          call mqc_vectorPhase(leftEigenVectors(:,i),.true.)
        endDo
      endIf
      if(PRESENT(rightEigenVectors)) then
        do i = 1,N
          call mqc_vectorPhase(rightEigenVectors(:,i),.true.)
        endDo
      endIf
      deAllocate(matrixCopy,Work)
!
      return
      end subroutine mqc_dgeev_eigensystem_Full

!hph
!!
!!PROCEDURE mqc_dsyev_eigensystem_symm
!      subroutine mqc_dsyev_eigensystem_symm(matrix,eigenValues,eigenVectors)
!!
!!     This subroutine wraps the LAPACK eigensystem routines for a symmetric
!!     matrix in full/dense memory storage format.
!!
!!     H. P. Hratchian, 2020
!!
!      implicit none
!      real(kind=real64),dimension(:,:),intent(in)::matrix
!      real(kind=real64),dimension(:),allocatable,intent(out)::eigenValues
!      real(kind=real64),dimension(:,:),allocatable,OPTIONAL,intent(out)::eigenVectors
!!
!      integer::i,N,INFO,LWORK
!      real(kind=real64),dimension(:,:),allocatable::matrixCopy
!      real(kind=real64),dimension(:),allocatable::Work
!      real(kind=real64),dimension(1)::WorkTmp
!      character(len=1)::JOBZ
!!
!!
!!     Do some initial checking of matrix to ensure it's square. Next, figure out
!!     the dimension of <matrix>, <N>. Then, allocate space for <matrixCopy> and
!!     copy matrix into that space so that the LAPACK routine DSYEV can use this
!!     space as working memory. Finally, allocate space for <eigenValues> and, if
!!     sent, <eigenVectors>.
!!
!      if(SIZE(matrix,1).ne.SIZE(matrix,2)) call mqc_error(  &
!        'mqc_eigensystem_symmMatrixFull: Rectangular matrix sent.')
!      N = SIZE(matrix,1)
!      Allocate(matrixCopy(N,N))
!      matrixCopy = matrix
!      if(Allocated(eigenValues)) then
!        if(SIZE(eigenValues).ne.N) then
!          DeAllocate(eigenValues)
!          Allocate(eigenValues(N))
!        endIf
!      else
!        Allocate(eigenValues(N))
!      endIf
!      if(PRESENT(eigenVectors)) then
!        JOBZ = 'V'
!        if(Allocated(eigenVectors)) then
!          if(SIZE(eigenVectors,1).ne.N.or.SIZE(eigenVectors,2).ne.N) then
!            DeAllocate(eigenVectors)
!            Allocate(eigenVectors(N,N))
!          endIf
!        else
!          Allocate(eigenVectors(N,N))
!        endIf
!      else
!        JOBZ = 'N'
!      endIf
!!
!!    Call LAPCK routine DSYEV. 
!!
!      call dsyev(JOBZ,'U',N,matrixCopy,N,eigenValues,WorkTmp,-1,INFO)    
!      if(INFO.ne.0) call mqc_error(  &
!        'mqc_eigensystem_symmMatrixFull: Failure in DSYEV getting WORK length.')
!      LWork = WorkTmp(1)
!      Allocate(Work(LWork))
!      call dsyev(JOBZ,'U',N,matrixCopy,N,eigenValues,Work,LWork,INFO)    
!      if(INFO.ne.0) call mqc_error(  &
!        'mqc_eigensystem_symmMatrixFull: Failure in DSYEV diagonalizing matrix.')
!      if(PRESENT(eigenVectors)) then
!        eigenVectors = matrixCopy
!        do i = 1,N
!          call mqc_vectorPhase(eigenVectors(:,i),.true.)
!        endDo
!      endIf
!
!      deAllocate(matrixCopy,Work)
!!
!      return
!      end subroutine mqc_eigensystem_symmFull
!hph
