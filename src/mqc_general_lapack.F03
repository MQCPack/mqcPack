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
      write(*,*)' Hrant - INFO = ',INFO
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
