    Program hartreefock 
!
!     This program perfoms Hartree-Fock calculations in 
!     different ways. 
!
      Use MQC_Algebra
      Use MQC_Gaussian
!
!     Variable Declarations...
!
      Implicit None
      Character(Len=256)::FileName
      Integer::IOut=6,NBasis=0,NElectrons=0,Multiplicity=0,IPrint=1,I,J,K,L,M, &
        NAlpha=0,NBeta=0,NAlpha_SO
      Real::Zero=0.0d0,One=1.0d0,Vnn,UHF_Energy
      Type(MQC_Molecule_Data)::MoleculeInfo
      Type(MQC_MO_Coefficients)::MO_Coeff
      Type(MQC_Fock_Matrix)::Fock_Matrix
      Type(MQC_Core_Hamiltonian)::Core_Hamiltonian
      Integer,Dimension(:),Allocatable::Basis2AtomMap
      Real,Dimension(3)::Temp_3Vector
      Real,Dimension(:,:),Allocatable::Overlap_Matrix,Alpha_Density,Beta_Density, &
        Total_Density,Contract_T,Contract_A,Contract_B,AA_Density,BB_Density,AB_Density, &
        BA_Density,CouCon_AA,CouCon_BB,ExcCon_AA,ExcCon_BB,ExcCon_AB,ExcCon_BA, &
        Occupied_Orbs,Full_Density
      Real,Dimension(:,:,:,:),Allocatable::ERIs
!
!     Format Statements
!
 1000 Format(1x,'Nuclear Repulsion Energy = ',F20.10,' au.')
 1010 Format(1x,'Electronic Energy = ',F20.10,' au.')
 1020 Format(1x,'Total Energy = ',F20.10,' au.')
!
!     Get the user-defined filename from the command line and then call the
!     routine that reads the Gaussian matrix file.
!
      Call Get_Command_Argument(1,FileName)
      Call MQC_Gaussian_Read_Matrix_File(.True.,FileName,MoleculeInfo,NBasis, &
        NElectrons,Multiplicity,Basis2AtomMap,Overlap_Matrix,MO_Coeff,Core_Hamiltonian, &
        Fock_Matrix,ERIs)

      If(IPrint.ge.1) then  
        Write(*,*)' Back from Read_Gaussian_Matrix_File'
        Write(*,*)' NAtoms       = ',MoleculeInfo%NAtoms
        Write(*,*)' NBasis       = ',NBasis
        Write(*,*)' NElectrons   = ',NElectrons
        Write(*,*)' Multiplicity = ',Multiplicity
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
!     Compute the nuclear-nuclear repulsion energy.
!
      Vnn = Zero 
      Do I = 1,MoleculeInfo%NAtoms-1
        Do J = I+1,MoleculeInfo%NAtoms
          Temp_3Vector = MoleculeInfo%Cartesian_Coordinates(:,I) -  &
            MoleculeInfo%Cartesian_Coordinates(:,J)
          Vnn = Vnn +  &
            MoleculeInfo%Atomic_Charges(I) *  &
            MoleculeInfo%Atomic_Charges(J) /  &
            SQRT(dot_product(Temp_3Vector,Temp_3Vector))
        EndDo
      EndDo
      Write(IOut,1000) Vnn
      Write(IOut,*)
!
!     Compute Hartree-Fock Energy contracting with Fock Matrix 
!
      NAlpha = (NElectrons + (Multiplicity-1))/2
      NBeta = NElectrons-NAlpha
      NAlpha_SO = NAlpha - NBeta

      UHF_Energy = Zero
      UHF_Energy = UHF_Energy + 0.5*Contraction((MatMul(MO_Coeff%Alpha(:,1:NAlpha),Transpose(MO_Coeff%Alpha(:,1:NAlpha)))+MatMul(MO_Coeff%Beta(:,1:NBeta),Transpose(MO_Coeff%Beta(:,1:NBeta)))),Core_Hamiltonian%Alpha)
      UHF_Energy = UHF_Energy + 0.5*Contraction(MatMul(MO_Coeff%Alpha(:,1:NAlpha),Transpose(MO_Coeff%Alpha(:,1:NAlpha))), &
        Fock_Matrix%Alpha)
      UHF_Energy = UHF_Energy + 0.5*Contraction(MatMul(MO_Coeff%Beta(:,1:NBeta),Transpose(MO_Coeff%Beta(:,1:NBeta))), &
        Fock_Matrix%Beta)
      Write(IOut,1010) UHF_Energy
      Write(IOut,*)
!
      Write(IOut,1020) UHF_Energy + Vnn
      Write(IOut,*)
!
!     Compute Hartree-Fock Energy with contracting with 2ERIs explicitly
!
      Allocate(Alpha_Density(NBasis,NBasis),Beta_Density(NBasis,NBasis),Total_Density(NBasis,NBasis), &
        Contract_T(NBasis,NBasis),Contract_A(NBasis,NBasis),Contract_B(NBasis,NBasis))
      Alpha_Density = MatMul(MO_Coeff%Alpha(:,1:NAlpha),Transpose(MO_Coeff%Alpha(:,1:NAlpha)))
      Beta_Density = MatMul(MO_Coeff%Beta(:,1:NBeta),Transpose(MO_Coeff%Beta(:,1:NBeta)))
      Total_Density = Alpha_Density + Beta_Density
      UHF_Energy = Zero
      UHF_Energy = UHF_Energy + 0.5*Contraction(Total_Density,Core_Hamiltonian%Alpha)
      UHF_Energy = UHF_Energy + 0.5*Contraction(Alpha_Density,Core_Hamiltonian%Alpha)
      UHF_Energy = UHF_Energy + 0.5*Contraction(Beta_Density,Core_Hamiltonian%Alpha)
      Write(IOut,*) 'Core energy contraction:',UHF_Energy
      Contract_T = Zero
      Contract_A = Zero
      Contract_B = Zero
      Do J = 1,NBasis
        Do K = 1,NBasis
          Do L = 1,NBasis
            Do M = 1,NBasis
              Contract_T(J,K) = Contract_T(J,K) + ERIs(J,K,M,L)*Total_Density(L,M)
              Contract_A(J,K) = Contract_A(J,K) + ERIs(J,L,M,K)*Alpha_Density(L,M)
              Contract_B(J,K) = Contract_B(J,K) + ERIs(J,L,M,K)*Beta_Density(L,M)
            EndDo
          EndDo
        EndDo
      EndDo
      UHF_Energy = UHF_Energy + 0.5*Contraction(Contract_T,Alpha_Density)
      UHF_Energy = UHF_Energy + 0.5*Contraction(Contract_T,Beta_Density)
      UHF_Energy = UHF_Energy - 0.5*Contraction(Contract_A,Alpha_Density)
      UHF_Energy = UHF_Energy - 0.5*Contraction(Contract_B,Beta_Density)
      Write(IOut,1010) UHF_Energy
      Write(IOut,*)
!
      Write(IOut,1020) UHF_Energy + Vnn
      Write(IOut,*)
!
!     Construct a made up density and contract with 2ERIs               
!
      Alpha_Density = MatMul(MO_Coeff%Alpha(:,1:NAlpha),Transpose(MO_Coeff%Alpha(:,1:NAlpha)))
      Beta_Density = MatMul(MO_Coeff%Beta(:,1:NBeta),Transpose(MO_Coeff%Beta(:,1:NBeta)))
      Total_Density = Alpha_Density + Beta_Density
      Call Print_Matrix(IOut,Alpha_Density,'Alpha Density')
      Call Print_Matrix(IOut,Beta_Density,'Beta Density')
      Call Print_Matrix(IOut,Total_Density,'Total Density')
      UHF_Energy = Zero
      UHF_Energy = UHF_Energy + Contraction(Total_Density,Core_Hamiltonian%Alpha)
      Write(IOut,*) 'Core energy contraction:',UHF_Energy
      Contract_A = Zero
      Contract_B = Zero
      Do J = 1,NBasis
        Do K = 1,NBasis
          Do L = 1,NBasis
            Do M = 1,NBasis
              Contract_A(J,K) = Contract_A(J,K) + ERIs(J,K,M,L)*Total_Density(L,M)
              Contract_B(J,K) = Contract_B(J,K) + 0.5*ERIs(J,L,M,K)*Total_Density(L,M)
            EndDo
          EndDo
        EndDo
      EndDo
      UHF_Energy = UHF_Energy + 0.5*Contraction(Contract_A,Total_Density)
      UHF_Energy = UHF_Energy - 0.5*Contraction(Contract_B,Total_Density)
      Write(IOut,1010) UHF_Energy
      Write(IOut,*)
!
      Write(IOut,1020) UHF_Energy + Vnn
      Write(IOut,*)
!
!     General Hartree-Fock Calculation
!
      Allocate(Occupied_Orbs(2*NBasis,2*NBasis),Full_Density(2*NBasis,2*NBasis))
      Allocate(AA_Density(NBasis,NBasis),BB_Density(NBasis,NBasis),AB_Density(NBasis,NBasis),BA_Density(NBasis,NBasis))
      Allocate(CouCon_AA(NBasis,NBasis),CouCon_BB(NBasis,NBasis),ExcCon_AA(NBasis,NBasis),ExcCon_BB(NBasis,NBasis), &
        ExcCon_AB(NBasis,NBasis),ExcCon_BA(NBasis,NBasis))
      Occupied_Orbs = Zero
!      Occupied_Orbs(1:NBasis,1:NAlpha) = MO_Coeff%Alpha(1:NBasis,1:NAlpha)
!      Occupied_Orbs(1:NBasis,NElectrons+1:NBeta+NBasis) = MO_Coeff%Alpha(1:NBasis,NAlpha+1:NBasis)
!      Occupied_Orbs(NBasis+1:2*NBasis,NAlpha+1:NElectrons) = MO_Coeff%Beta(1:NBasis,1:NBeta)
!      Occupied_Orbs(NBasis+1:2*NBasis,NBeta+NBasis+1:2*NBasis) = MO_Coeff%Beta(1:NBasis,NBeta+1:NBasis)
!   1 1   H  1S   ar    -0.63175   0.65993   0.16056  -0.39247
!   1 1   H  1S   br     0.31127  -0.23406   0.70308  -0.60698
!   2 2   H  1S   ar     0.23406   0.31127   0.60698   0.70308
!   2 2   H  1S   br     0.65993   0.63175  -0.39247  -0.16056
      Occupied_Orbs(1,1) = -0.63175
      Occupied_Orbs(3,1) =  0.31127
      Occupied_Orbs(2,1) =  0.23406
      Occupied_Orbs(4,1) =  0.65993
      Occupied_Orbs(1,2) =  0.65993
      Occupied_Orbs(3,2) = -0.23406
      Occupied_Orbs(2,2) =  0.31127
      Occupied_Orbs(4,2) =  0.63175
      Occupied_Orbs(1,3) =  0.16056
      Occupied_Orbs(3,3) =  0.70308
      Occupied_Orbs(2,3) =  0.60698
      Occupied_Orbs(4,3) = -0.39247
      Occupied_Orbs(1,4) = -0.39247
      Occupied_Orbs(3,4) = -0.60698
      Occupied_Orbs(2,4) =  0.70308
      Occupied_Orbs(4,4) = -0.16056
      Call Print_Matrix(IOut,Occupied_Orbs,'Occupied Orbs')
      Full_Density = MatMul(Occupied_Orbs(:,1:NBasis),Transpose(Occupied_Orbs(:,1:NBasis)))
      Call Print_Matrix(IOut,Full_Density,'Density Matrix')
      AA_Density = Full_Density(1:NBasis,1:NBasis) 
      BB_Density = Full_Density(NBasis+1:2*NBasis,NBasis+1:2*NBasis) 
      AB_Density = Full_Density(1:NBasis,NBasis+1:2*NBasis)
      BA_Density = Full_Density(NBasis+1:2*NBasis,1:NBasis)
      Call Print_Matrix(IOut,AA_Density,'Alpha-Alpha Density')
      Call Print_Matrix(IOut,BB_Density,'Beta-Beta Density')
      Call Print_Matrix(IOut,AB_Density,'Alpha-Beta Density')
      Call Print_Matrix(IOut,BA_Density,'Beta-Alpha Density')
      UHF_Energy = Zero
      UHF_Energy = UHF_Energy + Contraction(AA_Density,Core_Hamiltonian%Alpha)
      UHF_Energy = UHF_Energy + Contraction(BB_Density,Core_Hamiltonian%Alpha)
      Write(IOut,*) 'Core Hamiltonian Energy', UHF_Energy
      CouCon_AA = Zero
      CouCon_BB = Zero
      ExcCon_AA = Zero
      ExcCon_BB = Zero
      ExcCon_AB = Zero
      ExcCon_BA = Zero
      Do J = 1,NBasis
        Do K = 1,NBasis
          Do L = 1,NBasis
            Do M = 1,NBasis
              CouCon_AA(J,K) = CouCon_AA(J,K) + ERIs(J,K,L,M)*AA_Density(L,M)
              CouCon_BB(J,K) = CouCon_BB(J,K) + ERIs(J,k,L,M)*BB_Density(L,M)
              ExcCon_AA(J,K) = ExcCon_AA(J,K) + ERIs(J,L,K,M)*AA_Density(L,M)
              ExcCon_BB(J,K) = ExcCon_BB(J,K) + ERIs(J,L,K,M)*BB_Density(L,M)
              ExcCon_AB(J,K) = ExcCon_AB(J,K) + ERIs(J,L,K,M)*AB_Density(L,M)
              ExcCon_BA(J,K) = ExcCon_BA(J,K) + ERIs(J,L,K,M)*BA_Density(L,M)
            EndDo
          EndDo
        EndDo
      EndDo
      Call Print_Matrix(IOut,CouCon_AA,'Alpha-Alpha Coulomb')
      Call Print_Matrix(IOut,CouCon_BB,'Beta-Beta Coulomb')
      Call Print_Matrix(IOut,ExcCon_AA,'Alpha-Alpha Exchange')
      Call Print_Matrix(IOut,ExcCon_BB,'Beta-Beta Exchange')
      Call Print_Matrix(IOut,ExcCon_AB,'Alpha-Beta Exchange')
      Call Print_Matrix(IOut,ExcCon_BA,'Beta-Alpha Exchange')
      UHF_Energy  = UHF_Energy + 0.5*Contraction(CouCon_AA,AA_Density)
      UHF_Energy  = UHF_Energy + 0.5*Contraction(CouCon_AA,BB_Density) 
      UHF_Energy  = UHF_Energy + 0.5*Contraction(CouCon_BB,AA_Density)
      UHF_Energy  = UHF_Energy + 0.5*Contraction(CouCon_BB,BB_Density)
      UHF_Energy  = UHF_Energy - 0.5*Contraction(ExcCon_AA,AA_Density)
      UHF_Energy  = UHF_Energy - 0.5*Contraction(ExcCon_BB,BB_Density)
      UHF_Energy  = UHF_Energy - 0.5*Contraction(ExcCon_BA,AB_Density)
      UHF_Energy  = UHF_Energy - 0.5*Contraction(ExcCon_AB,BA_Density)
      Write(IOut,1010) UHF_Energy
      Write(IOut,*)
!
      Write(IOut,1020) UHF_Energy + Vnn
      Write(IOut,*)
!
 999  End Program hartreefock
