    Program ExtendedAP 
!
!     This program perfoms an extended AP calculation on an 
!     unrestricted Hartree-Fock or Density Functional Theory 
!     matrix file.
!
      Use Subroutines 
      Use MQC_PSCF
!
!
!     Variable Declarations...
!
      Implicit None
      Character(Len=256)::FileName
      Integer::IOut=6,NBasis=0,NElectrons=0,Multiplicity=0,IPrint=2,I=0,J=0,K=0,NBasRd=0, &
        NAlpha=0,NBeta=0,NAlpha_in_AS,NBeta_in_AS=0,NDets=0,NAlpha_Str,NBeta_Str, &
        NAlpha_in_Core,NBeta_in_Core,L_A_String,L_B_String,R_A_String,R_B_String, &
        L_Index,R_Index,Sz,NAlpha_SO,NBit_Ints,NBit_Ints_AS,Core_Population,Pos
      Real::Zero=0.0d0,One=1.0d0,Vnn
      Real::Thresh=0.99d0
      Type(MQC_Molecule_Data)::MoleculeInfo
      Type(MQC_MO_Coefficients)::MO_Coeff,Biorthog_Orbs
      Type(MQC_Fock_Matrix)::Fock_Matrix
      Type(MQC_Core_Hamiltonian)::Core_Hamiltonian,MO_Core_Ham
      Integer,Dimension(:),Allocatable::Basis2AtomMap,Magnetic_Orbitals,Alpha_Core_String, &
        Beta_Core_String,Temp_String
      Integer,Dimension(:,:),Allocatable::Alpha_Strings,Beta_Strings
      Real,Dimension(3)::Temp_3Vector
      Real,Dimension(:),Allocatable::Biorthog_Overlap,CI_Values,Final_Energy,State_S2, &
        Projected_UHF,Renorm_Val
      Real,Dimension(:,:),Allocatable::Overlap_Matrix,CI_Hamiltonian,CI_Vectors, &
        Full_Biorthog_Overlap,S2_Mat
      Real,Dimension(:,:,:,:),Allocatable::ERIs,MO_ERIs
!
!     Format Statements
!
 1000 Format(1x,'Nuclear Repulsion Energy = ',F20.10,' au.')
!
!     Get the user-defined filename from the command line and then call the
!     routine that reads the Gaussian matrix file.
!
      Call Get_Command_Argument(1,FileName)
      Call MQC_Gaussian_Read_Matrix_File(.False.,FileName,MoleculeInfo,NBasis, &
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
!     Compute Biorthogonalized Orbitals
!
      Call Biorthsub(IOut,IPrint,NBasis,NElectrons,Multiplicity, &
       Overlap_Matrix,MO_Coeff,Biorthog_Orbs,Biorthog_Overlap)           
      Call Print_Matrix(IOut,Biorthog_Orbs%Alpha,'Biorthogonalized Alpha Orbitals')
      Call Print_Matrix(IOut,Biorthog_Orbs%Beta,'Biorthogonalized Beta Orbitals')
      Call Print_Vector(IOut,Biorthog_Overlap,'Biorthogonal Overlap') 
!
      Allocate(Full_Biorthog_Overlap(NBasis*2,NBasis*2))
      Full_Biorthog_Overlap = Zero
      Full_Biorthog_Overlap(1:NBasis,1:NBasis) = MatMul(Transpose(Biorthog_Orbs%Alpha), &
        MatMul(Overlap_Matrix,Biorthog_Orbs%Alpha))
      Full_Biorthog_Overlap(1:NBasis,NBasis+1:NBasis*2) = MatMul(Transpose(Biorthog_Orbs%Alpha), &
        MatMul(Overlap_Matrix,Biorthog_Orbs%Beta))
      Full_Biorthog_Overlap(NBasis+1:NBasis*2,1:NBasis) = MatMul(Transpose(Biorthog_Orbs%Beta), &
        MatMul(Overlap_Matrix,Biorthog_Orbs%Alpha))
      Full_Biorthog_Overlap(NBasis+1:NBasis*2,NBasis+1:NBasis*2) = MatMul(Transpose(Biorthog_Orbs%Beta), &
        MatMul(Overlap_Matrix,Biorthog_Orbs%Beta))
      If(IPrint.ge.1) Call Print_Matrix(IOut,Full_Biorthog_Overlap,'Full Biorthogonalized Overlap')
      Deallocate(Overlap_Matrix)
!
!     Generate Slater Determinants       
!
      NAlpha = (NElectrons + (Multiplicity-1))/2
      NBeta = NElectrons-NAlpha
      NAlpha_SO = NAlpha - NBeta
      Write(IOut,*)'NAlpha:',NAlpha
      Write(IOut,*)'NBeta:',NBeta
      Write(IOut,*)'NAlpha_SO:',NAlpha_SO

      NAlpha_in_AS = 0
      NBeta_in_AS = 0
      NBasRd = 0
      Do I = 1, NBasis
        If(Biorthog_Overlap(I).le.Thresh) NBasRd = NBasRd + 1 
      EndDo
      If(NAlpha_SO.gt.0) NBasRd = NBasRd + NAlpha_SO
      Write(IOut,*) 'NBasRd', NBasRd
      Allocate(Magnetic_Orbitals(NBasRd))
      J = 1
      K = 1
      Do I = 1, NBasis
        If(Biorthog_Overlap(I).le.Thresh) then
          Magnetic_Orbitals(J) = I
          J = J + 1 
        ElseIf((NAlpha_SO.gt.0).and.(K.le.NAlpha_SO).and.((NAlpha-NAlpha_SO+K)).eq.I) then
          Magnetic_Orbitals(J) = I
          J = J + 1
          K = K + 1
        EndIf
      EndDo

      Deallocate(Biorthog_Overlap)

      Call Print_Vector(IOut,Magnetic_Orbitals,'Magnetic Orbital Pairs') 

      Do I = 1, NBasRd  
        If(Magnetic_Orbitals(I).le.NAlpha) NAlpha_in_AS = NAlpha_in_AS + 1
        If(Magnetic_Orbitals(I).le.NBeta) NBeta_in_AS = NBeta_in_AS + 1 
      EndDo
      NAlpha_in_Core = NAlpha - NAlpha_in_AS
      NBeta_in_Core = NBeta - NBeta_in_AS

      NAlpha_Str = Bin_Coeff(NBasRd,NAlpha_in_AS)
      NBeta_Str = Bin_Coeff(NBasRd,NBeta_in_AS)
      NDets = NAlpha_Str * NBeta_Str
      NBit_Ints = (NBasis/Bit_Size(0))+1
      NBit_Ints_AS = (NBasRd/Bit_Size(0))+1
      If(NBit_Ints_AS.gt.1) then
        Write(IOut,*) "ERROR: Determinant generator limited to ",Bit_Size(0)," orbitals"
        Call Exit()
      EndIf
      Allocate(Alpha_Strings(NAlpha_Str,NBit_Ints),Beta_Strings(NBeta_Str,NBit_Ints))
      Call Gen_Det_Str(IOut,IPrint,NBasRd,NAlpha_in_AS,NBeta_in_AS,Alpha_Strings,Beta_Strings)
      
      Allocate(Alpha_Core_String(NBit_Ints),Beta_Core_String(NBit_Ints),Temp_String(NBit_Ints))
      Alpha_Core_String = 0
      K = NBit_Ints
      Do I = 0, NAlpha-1
        If(Mod(I+1,Bit_Size(0)).eq.0) K = K-1
        Alpha_Core_String(K) = IBSet(Alpha_Core_String(K),I)
        Do J = 1, NBasRd
          If((I+1).eq.Magnetic_Orbitals(J)) Alpha_Core_String(K) = IBClr(Alpha_Core_String(K),I)
        EndDo
      EndDo
      Core_Population = 0 
      Do K = 1,NBit_Ints
        Core_Population = Core_Population + PopCnt(Alpha_Core_String(K))
      EndDo
!      Write(IOut,*) "Core Population", Core_Population
!      Write(IOut,'(B64)') Alpha_Core_String
      If(NAlpha_in_Core.ne.Core_Population) then
        Write(IOut,*) 'Number of bits set in Alpha_Core_String conflicts with &
          NAlpha_in_Core'
        Call Exit()
      EndIf
            
      Do I = 1, NAlpha_Str
        Temp_String = 0
        Do J = 1, NBasRd
          K = NBit_Ints - Magnetic_Orbitals(J)/Bit_Size(0)
          Pos = Mod(Magnetic_Orbitals(J),Bit_Size(0)) - 1
          Call MvBits(Alpha_Strings(I,1),J-1,1,Temp_String(K),Pos) 
        EndDo
        Do K = 1, NBit_Ints
          Alpha_Strings(I,K) = IOr(Temp_String(K),Alpha_Core_String(K))
        EndDo
      EndDo

      Beta_Core_String = 0
      K = NBit_Ints
      Do I = 0, NBeta-1
        If(Mod(I+1,Bit_Size(0)).eq.0) K = K-1
        Beta_Core_String(K) = IBSet(Beta_Core_String(K),I)
        Do J = 1, NBasRd
          If((I+1).eq.Magnetic_Orbitals(J)) Beta_Core_String(K) = IBClr(Beta_Core_String(K),I)
        EndDo
      EndDo
      Core_Population = 0 
      Do K = 1,NBit_Ints
        Core_Population = Core_Population + PopCnt(Beta_Core_String(K))
      EndDo
      If(NBeta_in_Core.ne.Core_Population) then
        Write(IOut,*) 'Number of bits set in Beta_Core_String conflicts with &
          NBeta_in_Core'
        Call Exit()
      EndIf
      
      Do I = 1, NBeta_Str
        Temp_String = 0
        Do J = 1, NBasRd
          K = NBit_Ints - Magnetic_Orbitals(J)/Bit_Size(0)
          Pos = Mod(Magnetic_Orbitals(J),Bit_Size(0)) - 1
          Call MvBits(Beta_Strings(I,1),J-1,1,Temp_String(K),Pos) 
        EndDo
        Do K = 1, NBit_Ints
          Beta_Strings(I,K) = IOr(Temp_String(K),Beta_Core_String(K))
        EndDo
      EndDo

      Write(IOut,*) "Number of configurations = ", NDets
      Write(IOut,*) "Alpha Strings"
      Write(IOut,'(B64)') Transpose(Alpha_Strings)
      Write(IOut,*) "Beta Strings"
      Write(IOut,'(B64)') Transpose(Beta_Strings)
!
!     Transform one and two-electron integrals to MO basis
!
      Allocate(MO_Core_Ham%Alpha(NBasis,NBasis),MO_Core_Ham%Beta(NBasis,NBasis))
      MO_Core_Ham%Alpha = Zero
      MO_Core_Ham%Beta = Zero
      MO_Core_Ham%Alpha = MatMul(MatMul(Transpose(Biorthog_Orbs%Alpha),Core_Hamiltonian%Alpha), &
        Biorthog_Orbs%Alpha)  
      MO_Core_Ham%Beta = MatMul(MatMul(Transpose(Biorthog_Orbs%Beta),Core_Hamiltonian%Beta), &
        Biorthog_Orbs%Beta)  
      If(IPrint.ge.2) then
        Call Print_Matrix(IOut,MO_Core_Ham%Alpha,' Alpha MO Core Hamiltonian')
        Call Print_Matrix(IOut,MO_Core_Ham%Beta,' Beta MO Core Hamiltonian')
      EndIf
      Call TwoERI_Trans(IOut,IPrint,Biorthog_Orbs,ERIs,MO_ERIs,.True.)
      Deallocate(ERIs)
!
!     Generate Hamiltonian and S**2 Matrix
!
      Allocate(CI_Hamiltonian(NDets,NDets),S2_Mat(NDets,NDets))
      CI_Hamiltonian = Zero
      S2_Mat = Zero
      If(IPrint.ge.2) Call Print_Matrix(IOut,CI_Hamiltonian,'CI Hamiltonian before')
      If(IPrint.ge.2) Call Print_Matrix(IOut,S2_Mat,'S^2 Matrix before')
      Do L_A_String = 1, NAlpha_Str  
        Do L_B_String = 1, NBeta_Str  
          Do R_A_String = 1, NAlpha_Str  
            Do R_B_String = 1, NBeta_Str  
              L_Index = 1+(L_B_String-1)*NAlpha_Str+(L_A_String-1) 
              R_Index = 1+(R_B_String-1)*NAlpha_Str+(R_A_String-1) 
!              Write(IOut,*)'L Index:',L_Index,' R Index:',R_Index
              CI_Hamiltonian(L_Index,R_Index) = Slater_Condon(IOut,IPrint,NBasis,Alpha_Strings(L_A_String,:), &
                Beta_Strings(L_B_String,:),Alpha_Strings(R_A_String,:),Beta_Strings(R_B_String,:),MO_Core_Ham,MO_ERIs,.True.)
              S2_Mat(L_Index,R_Index) = S2_Mat_Elem(IOut,IPrint,NBasis,Alpha_Strings(L_A_String,:), &
                Beta_Strings(L_B_String,:),Alpha_Strings(R_A_String,:),Beta_Strings(R_B_String,:),Full_Biorthog_Overlap)
            EndDo
          EndDo
        EndDo
      EndDo
      If(IPrint.ge.1) Call Print_Matrix(IOut,CI_Hamiltonian,'CI Hamiltonian')
      If(IPrint.ge.1) Call Print_Matrix(IOut,S2_Mat,'S^2 Matrix')
!
!     Diagonalize Hamiltonian and S**2 to get simultaneous eigenvectors
!
      Write(IOut,*) 'Diagonalizing CI Hamiltonian and S**2 Matrix'
      Do I = 1, NDets
        S2_Mat(I,I) = S2_Mat(I,I) + One
      EndDo
      Allocate(CI_Vectors(NDets,NDets),CI_Values(NDets))
      Call MQC_Matrix_Diagonalize(MatMul(CI_Hamiltonian,S2_Mat),CI_Vectors,CI_Values)
!      Call MQC_Matrix_Diagonalize(CI_Hamiltonian,CI_Vectors,CI_Values)
!      Call MQC_Matrix_Diagonalize(S2_Mat,CI_Vectors,CI_Values)
      If(IPrint.ge.1) then 
        Call Print_Matrix(IOut,CI_Vectors,'CI Eigenvectors')
        Call Print_Vector(IOut,CI_Values,'CI Eigenvalues')
      EndIf
      
      Allocate(State_S2(NDets))
      Do I = 1, NDets
        State_S2(I) = Zero
        Do J = 1, NDets
          Do K = 1, NDets 
            State_S2(I) = State_S2(I) + CI_Vectors(J,I)*CI_Vectors(K,I)*S2_Mat(J,K)
          EndDo
        EndDo
        State_S2(I) = State_S2(I) - One
      EndDo
      Call Print_Vector(IOut,State_S2,'State S**2 values') 
!
      Allocate(Final_Energy(NDets))
      Do I = 1, NDets
        CI_Values(I) = Zero
        Do J = 1, NDets
          Do K = 1, NDets 
            CI_Values(I) = CI_Values(I) + CI_Vectors(J,I)*CI_Vectors(K,I)*CI_Hamiltonian(J,K)
          EndDo
        EndDo
      EndDo
      Call Print_Vector(IOut,CI_Values,'Biorthogonal CI Energy') 
      Final_Energy = Vnn + CI_Values
      Call Print_Vector(IOut,Final_Energy,'Biorthogonal Energy including Vnn') 
!
!     Get UHF energy as function of pure spin states
!
      Allocate(Projected_UHF(NDets),Renorm_Val(NDets))
      Projected_UHF = MatMul(CI_Vectors**2,CI_Values) 
      Call Print_Vector(IOut,Projected_UHF,'Unrestricted Determinant Energies') 
!
!     Rearrange UHF equation and renormalize to get projected UHF energy
!
      Renorm_Val = Zero
      Do I = 1, NDets
        Do J = 1, NDets
!         Write(IOut,*)'State:',I,' S2:',State_S2(I),' Desired:',(Multiplicity**2)/4
!         Write(IOut,*)'Sz form S2:',(-1+SQRT(1+4*State_S2(I)))/2,' Sz from mult:',(Multiplicity-1)/2
!         Write(IOut,*)'na-nb form S2:',(-1+SQRT(1+4*State_S2(I))),' na-nb from mult:',(Multiplicity-1)
!         Write(IOut,*)'using this value to determine state:',Int((-1+SQRT(1+4*State_S2(I)))+0.1)
          If(Int((-1+SQRT(1+4*State_S2(J)))+0.1).eq.(Multiplicity-1)) then
            Renorm_Val(I) = Renorm_Val(I) + CI_Vectors(I,J)**2
          Else
            Projected_UHF(I) = Projected_UHF(I) -  CI_Vectors(I,J)**2 * CI_Values(J)
          EndIf
        EndDo
      EndDo
      Renorm_Val = 1/Renorm_Val
      If(IPrint.ge.2) Call Print_Vector(IOut,Renorm_Val,'Renormalization Constants')
      If(IPrint.ge.2) Call Print_Vector(IOut,Projected_UHF,'Raw Projected Value')
      Projected_UHF = Projected_UHF*Renorm_Val
      If(IPrint.ge.1) Call Print_Vector(IOut,Projected_UHF,'Full PUHF')
      Projected_UHF = Projected_UHF + Vnn
      Call Print_Vector(IOut,Projected_UHF,'Final Projected Energy')
      Write(IOut,*)'Projected Ground State Energy:',MinVal(Projected_UHF)
      Write(IOut,*)

 999  End Program ExtendedAP
