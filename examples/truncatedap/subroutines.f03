      Module Subroutines 

      Use MQC_Algebra
      Use MQC_Gaussian

      Contains

!======================================================================================================
!     
!     SUBROUTINE BIORTHOG
! 
      Subroutine Biorthsub(IOut,IPrint,NBasis,NElectrons,Multiplicity,Overlap_Matrix, &
        MO_Coeff,Biorthog_Orbs,Biorthog_Overlap)
!
!     This subroutine performs biorthogonalization of 
!     unrestricted orbitals
!
!     Variable Declarations...
!
      Implicit None
      Integer::IOut,NBasis,NElectrons,Multiplicity,NAlpha,NBeta,NOccA, &
        NOccB,NVirtA,NVirtB,ND,NSingA,NSingB,I,Ind,NVirtP,IPrint
      Real::Thresh=0.1d0
      Type(MQC_MO_Coefficients)::MO_Coeff,Biorthog_Orbs
      Real,Dimension(:),Allocatable::A_SVals,Biorthog_Overlap
      Real,Dimension(:,:),Intent(In)::Overlap_Matrix
      Real,Dimension(:,:),Allocatable::A_UVecs,A_VVecs, &
        MO_Overlap

      NAlpha = (NElectrons + (Multiplicity-1))/2
      NBeta = NElectrons-NAlpha
      NOccA = NAlpha
      NOccB = NBeta
      NVirtA = NBasis-NAlpha
      NVirtB = NBasis-NBeta
!
!     Allocate arrays in which to collect results
!
      Allocate(Biorthog_Overlap(NBasis))
      Call MQC_Gaussian_Fill_MO_Coefficients(NBasis,MO_Coeff%Alpha, &
        MO_Coeff%Beta,Biorthog_Orbs)
!
!     Do SVD on occupied alpha and beta orbitals and then rotate MO coefficients
!
      Write(*,*) 'Occupied Alpha - Occupied Beta'
      Allocate(MO_Overlap(NOccA,NOccB))
      MO_Overlap = MatMul(MatMul(Transpose(Biorthog_Orbs%Alpha(:,1:NOccA)),Overlap_Matrix), &
        Biorthog_Orbs%Beta(:,1:NOccB))
      If(IPrint.ge.1) Call Print_Matrix(IOut,MO_Overlap,' MO Overlap')
      Call MQC_Matrix_SVD(MO_Overlap,A_UVecs,A_VVecs,A_SVals)
      If(IPrint.ge.1) then  
        Call Print_Matrix(IOut,A_UVecs,'U Vectors')
        Call Print_Matrix(IOut,A_VVecs,'V Vectors')
      EndIf
      Call Print_Vector(IOut,A_SVals,'Singular Values')
      Biorthog_Orbs%Alpha(:,1:NOccA) = MatMul(Biorthog_Orbs%Alpha(:,1:NOccA),A_UVecs) 
      Biorthog_Orbs%Beta(:,1:NOccB) = MatMul(Biorthog_Orbs%Beta(:,1:NOccB),Transpose(A_VVecs)) 
      Biorthog_Overlap(1:Size(A_SVals,1)) = A_SVals
      If(IPrint.ge.1) then  
        Call Print_Matrix(IOut,Biorthog_Orbs%Alpha(:,1:NOccA),'Biorthogonalized Alpha Orbitals')
        Call Print_Matrix(IOut,Biorthog_Orbs%Beta(:,1:NOccB),'Biorthogonalized Beta Orbitals')
      EndIf
      ND = 0
      Do I = 1, Min(NOccA,NOccB) 
        If(A_SVals(I).ge.Thresh) ND = ND + 1
      EndDo
!     Gaussian code has this always set so we never have singly occupied beta
      ND = NBeta
      NSingA = NAlpha - ND
      NSingB = NBeta - ND
!
!     Do SVD on left-over alpha occupieds and beta virtals and then rotate MO coefficients
!
      If(NSingA.gt.0) then
        Write(*,*) 'Singly Occupied Alpha - Virtual Beta'
        Deallocate(MO_Overlap,A_SVals,A_VVecs,A_UVecs)
        Allocate(MO_Overlap(NSingA,NVirtB))
        MO_Overlap = MatMul(MatMul(Transpose(Biorthog_Orbs%Alpha(:,ND+1:ND+NSingA)),Overlap_Matrix), &
          Biorthog_Orbs%Beta(:,NOccB+1:NBasis))
        If(IPrint.ge.1) Call Print_Matrix(IOut,MO_Overlap,' MO Overlap')
        Call MQC_Matrix_SVD(MO_Overlap,A_UVecs,A_VVecs,A_SVals)
        If(IPrint.ge.1) then  
          Call Print_Matrix(IOut,A_UVecs,'U Vectors')
          Call Print_Matrix(IOut,A_VVecs,'V Vectors')
        EndIf
        Call Print_Vector(IOut,A_SVals,'Singular Values')
        Biorthog_Orbs%Alpha(:,ND+1:ND+NSingA) = MatMul(Biorthog_Orbs%Alpha(:,ND+1:ND+NSingA),A_UVecs) 
        Biorthog_Orbs%Beta(:,NOccB+1:NOccB+NVirtB) = MatMul(Biorthog_Orbs%Beta(:,NOccB+1:NOccB+NVirtB),Transpose(A_VVecs)) 
        Biorthog_Overlap(ND+1:ND+NSingA) = A_SVals
        If(IPrint.ge.1) then  
          Call Print_Matrix(IOut,Biorthog_Orbs%Alpha(:,ND+1:ND+NSingA),'Biorthogonalized Alpha Orbitals')
          Call Print_Matrix(IOut,Biorthog_Orbs%Beta(:,NOccB+1:NOccB+NVirtB),'Biorthogonalized Beta Orbitals')
        EndIf
      EndIf
!
!     Do SVD on left-over beta occupieds and alpha virtals and then rotate MO coefficients
!
      If(NSingB.gt.0) then
        Write(*,*) 'Singly Occupied Beta - Virtual Alpha'
        Deallocate(MO_Overlap,A_SVals,A_VVecs,A_UVecs)
        Allocate(MO_Overlap(NVirtA,NSingB))
        MO_Overlap = MatMul(MatMul(Transpose(Biorthog_Orbs%Alpha(:,NOccA+1:NBasis)),Overlap_Matrix), &
          Biorthog_Orbs%Beta(:,ND+1:ND+NSingB))
        If(IPrint.ge.1) Call Print_Matrix(IOut,MO_Overlap,' MO Overlap')
        Call MQC_Matrix_SVD(MO_Overlap,A_UVecs,A_VVecs,A_SVals)
        If(IPrint.ge.1) then  
          Call Print_Matrix(IOut,A_UVecs,'U Vectors')
          Call Print_Matrix(IOut,A_VVecs,'V Vectors')
        EndIf
        Call Print_Vector(IOut,A_SVals,'Singular Values')
        Biorthog_Orbs%Alpha(:,NOccA+1:NOccA+NVirtA) = MatMul(Biorthog_Orbs%Alpha(:,NOccA+1:NOccA+NVirtA),A_UVecs) 
        Biorthog_Orbs%Beta(:,ND+1:ND+NSingB) = MatMul(Biorthog_Orbs%Beta(:,ND+1:ND+NSingB),Transpose(A_VVecs)) 
        Biorthog_Overlap(ND+NSingA+1:ND+NSingA+NSingB) = A_SVals
        If(IPrint.ge.1) then  
          Call Print_Matrix(IOut,Biorthog_Orbs%Alpha(:,ND+1:ND+NSingA),'Biorthogonalized Alpha Orbitals')
          Call Print_Matrix(IOut,Biorthog_Orbs%Beta(:,NOccB+1:NOccB+NVirtB),'Biorthogonalized Beta Orbitals')
        EndIf
      EndIf
!
!     Do SVD on alpha and beta virtal orbitals and then rotate MO coefficients
!
      Ind = ND + NSingA + NSingB
      NVirtP = NBasis - Ind
      If(NVirtP.gt.0) then
        Write(*,*) 'Virtual Alpha - Virtual Beta'
        Deallocate(MO_Overlap,A_SVals,A_VVecs,A_UVecs)
        Allocate(MO_Overlap(NVirtP,NVirtP))
        MO_Overlap = MatMul(MatMul(Transpose(Biorthog_Orbs%Alpha(:,Ind+1:Ind+NVirtP)),Overlap_Matrix), &
          Biorthog_Orbs%Beta(:,Ind+1:Ind+NVirtP))
        If(IPrint.ge.1) Call Print_Matrix(IOut,MO_Overlap,' MO Overlap')
        Call MQC_Matrix_SVD(MO_Overlap,A_UVecs,A_VVecs,A_SVals)
        If(IPrint.ge.1) then  
          Call Print_Matrix(IOut,A_UVecs,'U Vectors')
          Call Print_Matrix(IOut,A_VVecs,'V Vectors')
        EndIf
        Call Print_Vector(IOut,A_SVals,'Singular Values')
        Biorthog_Orbs%Alpha(:,NOccA+1:NOccA+NVirtA) = MatMul(Biorthog_Orbs%Alpha(:,Ind+1:Ind+NVirtP),A_UVecs) 
        Biorthog_Orbs%Beta(:,Ind+1:Ind+NVirtP) = MatMul(Biorthog_Orbs%Beta(:,Ind+1:Ind+NVirtP),Transpose(A_VVecs)) 
        Biorthog_Overlap(Ind+1:Ind+NVirtP) = A_SVals
        If(IPrint.ge.1) then  
          Call Print_Matrix(IOut,Biorthog_Orbs%Alpha(:,Ind+1:Ind+NVirtP),'Biorthogonalized Alpha Orbitals')
          Call Print_Matrix(IOut,Biorthog_Orbs%Beta(:,Ind+1:Ind+NVirtP),'Biorthogonalized Beta Orbitals')
        EndIf
      EndIf
!
      If(IPrint.ge.1) then
        Call Print_Matrix(IOut,Biorthog_Orbs%Alpha,'Biorthogonalized Alpha Orbitals')
        Call Print_Matrix(IOut,Biorthog_Orbs%Beta,'Biorthogonalized Beta Orbitals')
        Call Print_Vector(IOut,Biorthog_Overlap,'Biorthogonal Overlap')
      EndIf
!
 999  End Subroutine Biorthsub
!
!======================================================================================================
!     
!     SUBROUTINE JADE
! 
      Subroutine Jade(IOut,IPrint,NDets,AMat,BMat,Vecs,AVals,BVals)
!
!     This subroutine performs joint approiximate diagonalization of two
!     matrices based on Jacobi rotations method 
!
!     Variable Declarations...
!
      Implicit None
      Integer::IOut,IPrint,NDets,I,J
      Real::Thresh=1.0d-12,Zero=0.0d0,One=1.0d0,TOn,TOff,Theta,JFac
      Logical::Converged
      Real,Dimension(:),Allocatable::RowI,RowJ,Temp
      Real,Dimension(:,:),Intent(In)::AMat,BMat
      Real,Dimension(:,:),Allocatable::Vecs,AVals,BVals,Mp,Mq
      Real,Dimension(2,2)::GMat
!
      Allocate(Vecs(NDets,NDets),AVals(NDets,NDets),BVals(NDets,NDets), &
        Mp(NDets,2),Mq(NDets,2),RowI(NDets*2),RowJ(NDets*2),Temp(NDets))
!
      AVals = AMat
      BVals = BMat
      Vecs = Zero
      ForAll(I=1:NDets) Vecs(I,I) = One
      Converged = .False.
      Do While (.not.Converged.and.NDets.gt.1)
        Do I = 1,NDets-1
          Do J = I+1,NDets
            GMat(1,1) = BVals(I,I)-BVals(J,J)
            GMat(1,2) = AVals(I,I)-AVals(J,J)
            GMat(2,1) = BVals(I,J)+BVals(J,I)
            GMat(2,2) = AVals(I,J)+AVals(J,I)
            GMat = MatMul(GMat,Transpose(GMat))
            If(IPrint.ge.2) Call Print_Matrix(IOut,GMat,'G Matrix')
            TOn = GMat(1,1)-GMat(2,2)
            TOff = GMat(1,2)+GMat(2,1)
            Theta = 0.5*ATan2(TOff,TOn+Sqrt(TOn*TOn+TOff*TOff))
            If(IPrint.ge.2) Write(IOut,*)'TOn:',TOn,' TOff:',TOff,' Theta:',Theta
            If(Abs(Sin(Theta)).gt.Thresh) then
              Mp(:,1) = BVals(:,I) 
              Mp(:,2) = AVals(:,I) 
              Mq(:,1) = BVals(:,J) 
              Mq(:,2) = AVals(:,J) 
            If(IPrint.ge.2) Call Print_Matrix(IOut,Mp,'Mp Matrix')
            If(IPrint.ge.2) Call Print_Matrix(IOut,Mq,'Mq Matrix')
              BVals(:,I) = Cos(Theta)*Mp(:,1)+Sin(Theta)*Mq(:,1)
              BVals(:,J) = Cos(Theta)*Mq(:,1)-Sin(Theta)*Mp(:,1)
              AVals(:,I) = Cos(Theta)*Mp(:,2)+Sin(Theta)*Mq(:,2)
              AVals(:,J) = Cos(Theta)*Mq(:,2)-Sin(Theta)*Mp(:,2)
              RowI(1:NDets) = BVals(I,:) 
              RowI(NDets+1:NDets*2) = AVals(I,:) 
              RowJ(1:NDets) = BVals(J,:) 
              RowJ(NDets+1:NDets*2) = AVals(J,:) 
              BVals(I,:) = Cos(Theta)*RowI(1:NDets)+Sin(Theta)*RowJ(1:NDets)
              AVals(I,:) = Cos(Theta)*RowI(NDets+1:NDets*2)+Sin(Theta)*RowJ(NDets+1:NDets*2)
              BVals(J,:) = Cos(Theta)*RowJ(1:NDets)-Sin(Theta)*RowI(1:NDets)
              AVals(J,:) = Cos(Theta)*RowJ(NDets+1:NDets*2)-Sin(Theta)*RowI(NDets+1:NDets*2)
              If(IPrint.ge.2) then
                Call Print_Matrix(IOut,AVals,'A Vals')
                Call Print_Matrix(IOut,BVals,'B Vals')
              EndIf
              Temp = Vecs(:,I)
              Vecs(:,I) = Cos(Theta)*Vecs(:,I)+Sin(Theta)*Vecs(:,J)
              Vecs(:,J) = Cos(Theta)*Vecs(:,J)-Sin(Theta)*Temp
              If(IPrint.ge.2) Call Print_Matrix(IOut,Vecs,'Vecs')
            Else
              Converged = .True.
            EndIf
          EndDo
        EndDo
      EndDo
      JFac = Zero
      Do I = 1, NDets 
        Do J = 1, NDets
          If (I.ne.J) then
            JFac = JFac + AVals(I,J)**2 + BVals(I,J)**2
          EndIf
        EndDo
      EndDo
      Write(IOut,*) 'Joint Diagonalization Ended, Sim Factor:',0.5*Sqrt(JFac)

 999  End Subroutine Jade
!
!======================================================================================================
!     
      End Module Subroutines
