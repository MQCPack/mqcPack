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
!     **                          Version 25.11.2                          **
!     **                         November 11, 2025                         **
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
      integer::i,j,N,INFO,LWORK
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
      call dsyev(JOBZ,'L',N,matrixCopy,N,eigenValues,WorkTmp,-1,INFO)    
      if(INFO.ne.0) call mqc_error(  &
        'mqc_eigensystem_symmMatrixFull: Failure in DSYEV getting WORK length.')
      LWork = WorkTmp(1)
      Allocate(Work(LWork))
      call dsyev(JOBZ,'L',N,matrixCopy,N,eigenValues,Work,LWork,INFO)    
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
        eigenValuesImaginary,leftEigenVectors,rightEigenVectors)
!
!     This subroutine wraps the LAPACK eigensystem routines for a general square
!     matrix. 
!
!     H. P. Hratchian, 2022
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::matrix
      real(kind=real64),dimension(:),allocatable,intent(out)::eigenValues
      real(kind=real64),dimension(:),allocatable,OPTIONAL,intent(out)::eigenValuesImaginary
      real(kind=real64),dimension(:,:),allocatable,OPTIONAL,intent(out)::leftEigenVectors,  &
        rightEigenVectors
!
      integer::i,N,LDVL,LDVR,INFO,LWORK
      real(kind=real64)::minImaginary,maxImaginary
      real(kind=real64),dimension(:,:),allocatable::matrixCopy
      real(kind=real64),dimension(:),allocatable::myEigenValuesImaginary
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
      Allocate(myEigenValuesImaginary(N))
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
        myEigenValuesImaginary,leftEigenVectors,LDVL,rightEigenVectors,LDVR,  &
        WorkTmp,-1,INFO)
      if(INFO.ne.0) call mqc_error(  &
        'mqc_eigensystem_MatrixFull: Failure in DGEEV getting WORK length.')
      LWork = WorkTmp(1)
      Allocate(Work(LWork))
      call dgeev(JOBVL,JOBVR,N,matrixCopy,N,eigenValues,  &
        myEigenValuesImaginary,leftEigenVectors,LDVL,rightEigenVectors,LDVR,  &
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
!     Check to see if eigenvalues are complex.
!
      if(PRESENT(eigenValuesImaginary)) then
        if(Allocated(eigenValuesImaginary)) then
          if(SIZE(eigenValuesImaginary).ne.N) then
            DeAllocate(eigenValuesImaginary)
            Allocate(eigenValuesImaginary(N))
          endIf
        else
          Allocate(eigenValuesImaginary(N))
        endIf
        eigenValuesImaginary = myEigenValuesImaginary
      else
        minImaginary = ABS(MINVAL(myEigenValuesImaginary))
        maxImaginary = ABS(MAXVAL(myEigenValuesImaginary))
        if(minImaginary.gt.MQC_small.or.maxImaginary.gt.MQC_small)  &
          call mqc_error('mqc_dgeev_eigensystem_Full: Complex eigenvalues.')
      endIf
!
      return
      end subroutine mqc_dgeev_eigensystem_Full

!
!PROCEDURE mqc_dsyevr_eigensystem_symmFull
      subroutine mqc_dsyevr_eigensystem_symmFull(matrix,eigenValues,eigenVectors)
!
!     Wraps LAPACK DSYEVR to compute eigenvalues (and optionally eigenvectors)
!     of a real symmetric matrix in full/dense storage.
!
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::matrix
      real(kind=real64),dimension(:),allocatable,intent(out)::eigenValues
      real(kind=real64),dimension(:,:),allocatable,optional,intent(out)::eigenVectors
!
      integer::i,N,INFO,LDA,LDZ,LWORK,LIWORK,M
      character(len=1)::JOBZ
      real(kind=real64)::ABSTOL
      real(kind=real64),dimension(:,:),allocatable::matrixCopy,Zloc
      real(kind=real64),dimension(:),allocatable::Work,Wtmp
      integer,dimension(:),allocatable::IWork,IWtmp,ISUPPZ
!
!     Quick checks and setup
!
      if(SIZE(matrix,1).ne.SIZE(matrix,2)) call mqc_error( &
        'mqc_dsyevr_eigensystem_symmFull: Rectangular matrix sent.')
      N   = SIZE(matrix,1)
      LDA = N
!
!     Copy input so LAPACK can overwrite it.
!
      allocate(matrixCopy(N,N))
      matrixCopy = matrix
!
!     Prepare outputs...
!
      if (allocated(eigenValues)) then
         if (SIZE(eigenValues) .ne. N) then
            deallocate(eigenValues)
            allocate(eigenValues(N))
         end if
      else
         allocate(eigenValues(N))
      end if
      if(present(eigenVectors)) then
         JOBZ = 'V'
         LDZ  = N
         if (allocated(eigenVectors)) then
            if (SIZE(eigenVectors,1) .ne. N .or. SIZE(eigenVectors,2) .ne. N) then
               deallocate(eigenVectors)
               allocate(eigenVectors(N,N))
            end if
         else
            allocate(eigenVectors(N,N))
         end if
         allocate(Zloc(N,N))
      else
         JOBZ = 'N'
         LDZ  = 1
!        Minimal dummy Z to satisfy the interface (not referenced when JOBZ='N')
         allocate(Zloc(1,1))
      end if
!
!     ISUPPZ size: 2*max(1,M); we don't know M until after the call, but for
!     RANGE='A' we have M=N. Allocate for N.
!
      allocate(ISUPPZ(2*max(1,N)))
!
!     Workspace check...
!
      allocate(Wtmp(1), IWtmp(1))
      ABSTOL = 2.0_real64 * mqc_small
      INFO   = 0
      call dsyevr(JOBZ,'A','L',N,matrixCopy,LDA,0.0_real64,0.0_real64,0,0,&
        ABSTOL,M,eigenValues,Zloc,LDZ,ISUPPZ,Wtmp,-1,IWtmp,-1,INFO)
      if(INFO.ne.0) call mqc_error( &
     &  'mqc_dsyevr_eigensystem_symmFull: Workspace query failed in DSYEVR.')

      LWORK  = int(Wtmp(1))
      LIWORK = IWtmp(1)
      deallocate(Wtmp,IWtmp)
      allocate(Work(LWORK),IWork(LIWORK))
!
!     Actual computation (RANGE='A' ⇒ all eigenpairs; eigenvalues ascending)
!
      INFO = 0
      call dsyevr(JOBZ,'A','L',N,matrixCopy,LDA,0.0_real64,0.0_real64,0,0,&
        ABSTOL,M,eigenValues,Zloc,LDZ,ISUPPZ,Work,LWORK,IWork,LIWORK,INFO)
      if(INFO.ne.0) call mqc_error(  &
        'mqc_dsyevr_eigensystem_symmFull: Failure in DSYEVR.')
!
!     If vectors requested, return them (columns) and phase-normalize like DSYEV wrapper
!
      if(present(eigenVectors)) then
         eigenVectors = Zloc(:,1:M)   ! here M should equal N for RANGE='A'
         do i = 1,N
            call mqc_vectorPhase(eigenVectors(:,i),.true.)
         end do
      end if

!     Cleanup
      deallocate(matrixCopy,Zloc,ISUPPZ,Work,IWork)
!
      return
      end subroutine mqc_dsyevr_eigensystem_symmFull

! 
!PROCEDURE mqc_dgesvd_svd_Full
      subroutine mqc_dgesvd_svd_Full(matrix,singularValues,  &
        leftSingularVectors,rightSingularVectorsT)
!
!     This subroutine wraps the LAPACK SVD routine DGESVD for a general real
!     rectangular matrix A. If the optional outputs are PRESENT, it returns full
!     U (M×M) and full VT (N×N).
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(in)::matrix
      real(kind=real64),dimension(:),allocatable,intent(out)::singularValues
      real(kind=real64),dimension(:,:),allocatable,OPTIONAL,intent(out)::leftSingularVectors
      real(kind=real64),dimension(:,:),allocatable,OPTIONAL,intent(out)::rightSingularVectorsT
!
      integer::M,N,LDU,LDVT,INFO,LWORK
      real(kind=real64),dimension(:,:),allocatable::matrixCopy
      real(kind=real64),dimension(:),allocatable::Work
      real(kind=real64),dimension(1)::WorkTmp
      real(kind=real64),dimension(1,1)::Utmp,VTtmp
      character(len=1)::JOBU,JOBVT
!
!     Do initial checks and figure out dimensions for arrays needed below. Then
!     allocate space and figure out a few LAPack flags.
!
      M = SIZE(matrix,1)
      N = SIZE(matrix,2)
      Allocate(matrixCopy(M,N))
      matrixCopy = matrix
      if(Allocated(singularValues)) then
        if(SIZE(singularValues).ne.min(M,N)) then
          DeAllocate(singularValues)
          Allocate(singularValues(min(M,N)))
        endIf
      else
        Allocate(singularValues(min(M,N)))
      endIf
      if(PRESENT(leftSingularVectors)) then
        JOBU = 'A'
        LDU = M
        if(Allocated(leftSingularVectors)) then
          if(SIZE(leftSingularVectors,1).ne.M.or.SIZE(leftSingularVectors,2).ne.M) then
            DeAllocate(leftSingularVectors)
            Allocate(leftSingularVectors(M,M))
          endIf
        else
          Allocate(leftSingularVectors(M,M))
        endIf
      else
        JOBU = 'N'
        LDU = 1
      endIf
      if(PRESENT(rightSingularVectorsT)) then
        JOBVT = 'A'
        LDVT = N
        if(Allocated(rightSingularVectorsT)) then
          if(SIZE(rightSingularVectorsT,1).ne.N.or.SIZE(rightSingularVectorsT,2).ne.N) then
            DeAllocate(rightSingularVectorsT)
            Allocate(rightSingularVectorsT(N,N))
          endIf
        else
          Allocate(rightSingularVectorsT(N,N))
        endIf
      else
        JOBVT = 'N'
        LDVT = 1
      endIf
      INFO = 0
      if(JOBU.eq.'N') then
        Utmp(1,1) = mqc_float(0)
      endIf
      if(JOBVT.eq.'N') then
        VTtmp(1,1) = mqc_float(0)
      endIf
!
!     Call LAPack routine DGESVD.
!
      call dgesvd(JOBU,JOBVT,M,N,matrixCopy,M,singularValues,  &
        merge(leftSingularVectors,Utmp, PRESENT(leftSingularVectors)),LDU,  &
        merge(rightSingularVectorsT,VTtmp, PRESENT(rightSingularVectorsT)),LDVT,  &
        WorkTmp,-1,INFO)
      if(INFO.ne.0) call mqc_error(  &
        'mqc_dgesvd_svd_Full: Failure in DGESVD getting WORK length.')
      LWORK = int(WorkTmp(1))
      if(LWORK.lt.1) LWORK = 1
      Allocate(Work(LWORK))
      call dgesvd(JOBU,JOBVT,M,N,matrixCopy,M,singularValues,  &
        merge(leftSingularVectors,Utmp, PRESENT(leftSingularVectors)),LDU,  &
        merge(rightSingularVectorsT,VTtmp, PRESENT(rightSingularVectorsT)),LDVT,  &
        Work,LWORK,INFO)
      if(INFO.gt.0) then
        call mqc_error('mqc_dgesvd_svd_Full: DGESVD did not converge.')
      elseIf(INFO.lt.0) then
        call mqc_error('mqc_dgesvd_svd_Full: Illegal argument in DGESVD call.')
      endIf
      DeAllocate(matrixCopy,Work)
!
      return
      end subroutine mqc_dgesvd_svd_Full



!PROCEDURE mqc_dgesv_mat
      subroutine mqc_dgesv_mat(A,B,X,infoOut)
!
!     This routine is a wrapper for the LAPack routine DGESV for cases where the
!     right hand side of the linear system of equations is a matrix (or multiple
!     vectors).
!
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      real(kind=real64),intent(in)::A(:,:)
      real(kind=real64),intent(in)::B(:,:)
      real(kind=real64),allocatable,intent(out)::X(:,:)
      integer,intent(out),optional::infoOut
!
      integer::N,NRHS,LDA,LDB,INFO
      integer,dimension(:),allocatable::IPIV
      real(kind=real64),dimension(:,:),allocatable::Acopy,Bcopy
!
!     Carry out shape checks.
!
      if(SIZE(A,1).ne.SIZE(A,2)) call mqc_error('mqc_dgesv: A must be square.')
      N = SIZE(A,1)
      if (SIZE(B,1).ne.N) call mqc_error('mqc_dgesv: B has incompatible row dimension.')
      NRHS = SIZE(B,2)
      if (NRHS < 1) call mqc_error('mqc_dgesv: NRHS must be >= 1.')
!
!     Allocate and then fill working copies of A and B.
!
      allocate(Acopy(N,N),Bcopy(N,NRHS))
      Acopy = A
      Bcopy = B
!
!     Allocate output space.
!
      if(allocated(X)) then
        if(SIZE(X,1).ne.N .or. SIZE(X,2).ne.NRHS) then
          deallocate(X)
          allocate(X(N,NRHS))
        endIf
      else
        allocate(X(N,NRHS))
      endIf
!
!     Call DGESV...
!
      LDA = N
      LDB = N
      allocate(IPIV(N))
      INFO = 0
      call dgesv(N,NRHS,Acopy,LDA,IPIV,Bcopy,LDB,INFO)
!
!     Check for errors in the LAPack run.
!
      if(present(infoOut)) then
        infoOut = INFO
        if(INFO.ne.0) return
      endIf
      call mqc_error('mqc_dgesv_mat: DGESV failed (INFO /= 0).')
!
!     Move Bcopy into the output argument X. Then, DeAllocate the temporary
!     memory space this routine has allocated.
!
      X = Bcopy
      deallocate(Acopy,Bcopy,IPIV)
!
      return
      end subroutine mqc_dgesv_mat


!PROCEDURE mqc_dgesv_vec
      subroutine mqc_dgesv_vec(A,b,x,infoOut)
!
!     This routine is a wrapper for the LAPack routine DGESV for cases where the
!     right hand side of the linear system of equations is a vector.
!
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      real(kind=real64),intent(in)::A(:,:)
      real(kind=real64),intent(in)::b(:)
      real(kind=real64),allocatable,intent(out)::x(:)
      integer(kind=int64),intent(out),optional::infoOut
!
      integer(kind=int64)::N,NRHS,LDA,LDB,INFO
      integer(kind=int64),allocatable::IPIV(:)
      real(kind=real64),allocatable::Acopy(:,:),Bcopy(:,:)
!
!     Carry out shape checks.
!
      if(SIZE(A,1).ne.SIZE(A,2)) call mqc_error('mqc_dgesv: A must be square.')
      N = SIZE(A,1)
      if(SIZE(b).ne.N) call mqc_error('mqc_dgesv: b has incompatible length.')
!
!     Allocate and then fill working copies of A and B.
!
      NRHS = 1
      allocate(Acopy(N,N),Bcopy(N,NRHS))
      Acopy = A
      Bcopy(:,1) = b
!
!     Allocate output space.
!
      if(allocated(x)) then
        if(SIZE(x).ne.N) then
          deallocate(x)
          allocate(x(N))
        endIf
      else
        allocate(x(N))
      endIf
!
!     Call DGESV...
!
      LDA = N
      LDB = N
      allocate(IPIV(N))
      INFO = 0
      call dgesv(N,NRHS,Acopy,LDA,IPIV,Bcopy,LDB,INFO)
!
!     Check for errors in the LAPack run.
!
      if(present(infoOut)) then
        infoOut = INFO
        if(INFO.ne.0) return
      endIf
      if(INFO.ne.0)  &
        call mqc_error('mqc_dgesv_vec: DGESV failed (INFO /= 0).')
!
!     Move Bcopy into the output argument X. Then, DeAllocate the temporary
!     memory space this routine has allocated.
!
      x = Bcopy(:,1)
      deallocate(Acopy,Bcopy,IPIV)
!
      return
      end subroutine mqc_dgesv_vec
