      Program Biorthog 
!
!     This program perfoms a borthogonalization between alpha and beta
!     orbitals obtained from an unrestricted Hartree-Fock or Density 
!     Functional Theory Calculation.
!
      Use MQC_Algebra
      Use MQC_Gaussian
!
!
!     Variable Declarations...
!
      Implicit None
      Character(Len=256)::FileName
      Integer::IOut=6,NBasis=0,NElectrons=0,Multiplicity=0,NAlpha=0, &
        NBeta=0,NOccA=0,NOccB=0,NVirtA=0,NVirtB=0,ND=0,NSingA=0,NSingB=0, &
        I=0,Ind=0,NVirtP=0,IPrint=1 
      Real::Thresh=0.1d0
      Type(MQC_Molecule_Data)::MoleculeInfo
      Type(MQC_MO_Coefficients)::MO_Coeff,Biorthog_Orbs
      Type(MQC_Fock_Matrix)::Fock_Matrix
      Type(MQC_Core_Hamiltonian)::Core_Hamiltonian
      Integer,Dimension(:),Allocatable::Basis2AtomMap
      Real,Dimension(:),Allocatable::A_SVals,Biorthog_Overlap
      Real,Dimension(:,:),Allocatable::Overlap_Matrix,A_UVecs,A_VVecs, &
        MO_Overlap
!
!
!     Get the user-defined filename from the command line and then call the
!     routine that reads the Gaussian matrix file.
!
      Call Get_Command_Argument(1,FileName)
      Call MQC_Gaussian_Read_Matrix_File(.True.,FileName,MoleculeInfo,NBasis, &
        NElectrons,Multiplicity,Basis2AtomMap,Overlap_Matrix,MO_Coeff,Core_Hamiltonian, &
        Fock_Matrix)

      If(IPrint.ge.1) then  
        Write(*,*)' Back from Read_Gaussian_Matrix_File'
        Write(*,*)' NAtoms       = ',MoleculeInfo%NAtoms
        Write(*,*)' NBasis       = ',NBasis
        Write(*,*)' NElectrons   = ',NElectrons
        Write(*,*)' Multiplicity = ',Multiplicity
      EndIf

      NAlpha = (NElectrons + (Multiplicity-1))/2
      NBeta = NElectrons-NAlpha
      NOccA = NAlpha
      NOccB = NBeta
      NVirtA = NBasis-NAlpha
      NVirtB = NBasis-NBeta

      If(IPrint.ge.1) then  
        Write(*,*)' NAlpha = ',NAlpha
        Write(*,*)' NBeta  = ',NBeta
        Write(*,*)' NOccA  = ',NOccA
        Write(*,*)' NOccB  = ',NOccB
        Write(*,*)' NVirtA = ',NVirtA
        Write(*,*)' NVirtB = ',NVirtB
        Call Print_Vector(IOut,MoleculeInfo%Atomic_Numbers,  &
          'Atomic Numbers')
        Call Print_Vector(IOut,MoleculeInfo%Atomic_Masses,  &
          'Atomic Masses')
        Call Print_Vector(IOut,MoleculeInfo%Atomic_Charges,  &
          'Atomic Charges')
        Call Print_Matrix(IOut,MoleculeInfo%Cartesian_Coordinates,  &
          'Cartesian Coordinates')
        Call Print_Vector(IOut,Basis2AtomMap,'Basis-Atom Mapping')
        Call Print_Matrix(IOut,Overlap_Matrix,'Overlap Matrix')
        Call Print_Matrix(IOut,Core_Hamiltonian%Alpha,'Core Hamiltonian Alpha')
        Call Print_Matrix(IOut,Core_Hamiltonian%Beta,'Core Hamiltonian Beta')
        Call Print_Matrix(IOut,Fock_Matrix%Alpha,'Alpha Fock Matrix')
        Call Print_Matrix(IOut,Fock_Matrix%Beta,'Beta Fock Matrix')
        Call Print_Matrix(IOut,MO_Coeff%Alpha,'Alpha MO Coefficient Matrix')
        Call Print_Matrix(IOut,MO_Coeff%Beta,'Beta MO Coefficient Matrix')
      EndIf
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
        Call Print_Matrix(IOut,Biorthog_Orbs%Alpha,'Biorthogonalized Alpha Orbitals')
        Call Print_Matrix(IOut,Biorthog_Orbs%Beta,'Biorthogonalized Beta Orbitals')
        Call Print_Vector(IOut,Biorthog_Overlap,'Biorthogonal Overlap')
!
 999  End Program Biorthog
