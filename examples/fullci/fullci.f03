      Program FullCI     
!
!     This program perfoms a full CI calculation initially for a
!     restricted HF or DFT orbitals matrix file.
!
      Use MQC_Algebra
      Use MQC_Gaussian
      Use MQC_PSCF
!
!     Variable Declarations...
!
      Implicit None
      Character(Len=256)::FileName
      Integer::IOut=6,NBasis=0,NElectrons=0,Multiplicity=0,IPrint=0,NAlpha=0, &
        NBeta=0,NAlpha_Str,NBeta_Str,NDets,L_A_String,L_B_String,R_A_String, &
        R_B_String,L_Index,R_Index,I,J,K
      Real::One=1.0d0,Thresh=0.999999d0,Zero=0.0d0,Vnn
      Type(MQC_Molecule_Data)::MoleculeInfo
      Type(MQC_MO_Coefficients)::MO_Coeff
      Type(MQC_Fock_Matrix)::Fock_Matrix
      Type(MQC_Core_Hamiltonian)::Core_Hamiltonian,MO_Core_Ham
      Integer,Dimension(:),Allocatable::Basis2AtomMap
      Integer,Dimension(:,:),Allocatable::Alpha_Strings,Beta_Strings
      Real,Dimension(3)::Temp_3Vector
      Real,Dimension(:),Allocatable::CI_Values,Final_Energy,State_S2
      Real,Dimension(:,:),Allocatable::Overlap_Matrix,CI_Hamiltonian,CI_Vectors,Full_MO_Overlap,S2_Mat
      Real,Dimension(:,:,:,:),Allocatable::ERIs,MO_ERIs
      Logical::UHF
!
!     Format Statements
!
 1000 Format(1x,'Nuclear Repulsion Energy = ',F20.10,' au.')
!
!     Temporary logical for UHF while backend is sorted
!
      UHF = .True.
!      UHF = .False.
!
!     Get the user-defined filename from the command line and then call the
!     routine that reads the Gaussian matrix file.
!
      Call Get_Command_Argument(1,FileName)
      Call MQC_Gaussian_Read_Matrix_File(.False.,FileName,MoleculeInfo,NBasis, &
        NElectrons,Multiplicity,Basis2AtomMap,Overlap_Matrix,MO_Coeff,Core_Hamiltonian, &
        Fock_Matrix,ERIs)

      If(IPrint.ge.2) then  
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
!
!     Generate Slater Determinants       
!
      NAlpha = (NElectrons + (Multiplicity-1))/2
      NBeta = NElectrons-NAlpha
      Call Gen_Det_Str(IOut,IPrint,NBasis,NAlpha,NBeta,Alpha_Strings,Beta_Strings)

      NAlpha_Str = Size(Alpha_Strings,1) 
      NBeta_Str = Size(Beta_Strings,1) 
      NDets = NAlpha_Str * NBeta_Str 
      Write(IOut,*) "Number of configurations = ", NDets
      Write(IOut,*) "Alpha Strings"
      Write(IOut,'(B64)') Alpha_Strings
      Write(IOut,*) "Beta Strings"
      Write(IOut,'(B64)') Beta_Strings
!
!     Transform one and two-electron integrals to MO basis
!
      Allocate(MO_Core_Ham%Alpha(NBasis,NBasis),MO_Core_Ham%Beta(NBasis,NBasis))
      MO_Core_Ham%Alpha = Zero
      MO_Core_Ham%Beta = Zero
      MO_Core_Ham%Alpha = MatMul(MatMul(Transpose(MO_Coeff%Alpha),Core_Hamiltonian%Alpha), &
        MO_Coeff%Alpha)  
      If(UHF.eq..True.) MO_Core_Ham%Beta = MatMul(MatMul(Transpose(MO_Coeff%Beta), & 
        Core_Hamiltonian%Beta),MO_Coeff%Beta)  
      If(IPrint.ge.2) Call Print_Matrix(IOut,MO_Core_Ham%Alpha,' Alpha MO Core Hamiltonian')
      If(IPrint.ge.2.and.UHF.eq..True.) Call Print_Matrix(IOut,MO_Core_Ham%Beta, &
        ' Beta MO Core Hamiltonian')
      Call TwoERI_Trans(IOut,IPrint,MO_Coeff,ERIs,MO_ERIs,UHF)
!
!     Generate Hamiltonian Matrix
!
      Allocate(CI_Hamiltonian(NDets,NDets))
      CI_Hamiltonian = Zero
      If(IPrint.ge.2) Call Print_Matrix(IOut,CI_Hamiltonian,'CI Hamiltonian')
      Do L_A_String = 1, NAlpha_Str  
        Do L_B_String = 1, NBeta_Str  
          Do R_A_String = 1, NAlpha_Str  
            Do R_B_String = 1, NBeta_Str  
              L_Index = 1+(L_B_String-1)*NAlpha_Str+(L_A_String-1) 
              R_Index = 1+(R_B_String-1)*NAlpha_Str+(R_A_String-1) 
              CI_Hamiltonian(L_Index,R_Index) = Slater_Condon(IOut,IPrint,NBasis,Alpha_Strings(L_A_String,1), &
                Beta_Strings(L_B_String,1),Alpha_Strings(R_A_String,1),Beta_Strings(R_B_String,1),MO_Core_Ham,MO_ERIs,UHF)
            EndDo
          EndDo
        EndDO
      EndDo
      If(IPrint.ge.1) Call Print_Matrix(IOut,CI_Hamiltonian,'CI Hamiltonian')
!
!     Diagonalize Hamiltonian
!
      Write(IOut,*) 'Diagonalizing CI Hamiltonian'
      Allocate(CI_Vectors(NDets,NDets),CI_Values(NDets))
      Call MQC_Matrix_Diagonalize(CI_Hamiltonian,CI_Vectors,CI_Values)
      If(IPrint.ge.1) then 
        Call Print_Matrix(IOut,CI_Vectors,'CI Eigenvectors')
        Call Print_Vector(IOut,CI_Values,'CI Eigenvalues')
      EndIf
!
      Allocate(Final_Energy(NDets))
      Final_Energy = Vnn + CI_Values
      Call Print_Vector(IOut,Final_Energy,'Final Energy') 
!
 999  End Program FullCI     
