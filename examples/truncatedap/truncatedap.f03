    Program truncatedAP 
!
!     This program perfoms an truncated AP calculation on an 
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
      Integer::IOut=6,NBasis=0,NElectrons=0,Multiplicity=0,IPrint=1,I,J,K,P,Q,NBasRd=0, &
        NAlpha=0,NBeta=0,NAlpha_in_AS,NBeta_in_AS=0,NDets,NDMax,NAlpha_Str,NBeta_Str, &
        NAlpha_in_Core,NBeta_in_Core,Alpha_Core_String,Beta_Core_String,Temp_String, &
        L_A_String,L_B_String,R_A_String,R_B_String,L_Index,R_Index,Sz,NAlpha_SO,Max_SR, &
        SR_Trunc=1 
      Real::Zero=0.0d0,One=1.0d0,Vnn,UHF_Energy,Fib_Norm
      Real::Thresh=0.99d0
      Type(MQC_Molecule_Data)::MoleculeInfo
      Type(MQC_MO_Coefficients)::MO_Coeff,Biorthog_Orbs
      Type(MQC_Fock_Matrix)::Fock_Matrix
      Type(MQC_Core_Hamiltonian)::Core_Hamiltonian,MO_Core_Ham
      Integer,Dimension(:),Allocatable::Basis2AtomMap,Magnetic_Orbitals,SR_Index
      Integer,Dimension(:,:),Allocatable::Alpha_Strings,Beta_Strings
      Real,Dimension(3)::Temp_3Vector
!      Real,Dimension(:),Allocatable::Biorthog_Overlap,CI_Values,Final_Energy,State_S2, &
!        Det_Energy,Projected_UHF,Renorm_Val
      Real,Dimension(:),Allocatable::Biorthog_Overlap,Final_Energy,Det_Energy, &
        Final_S2,CI_Values
!      Real,Dimension(:,:),Allocatable::Overlap_Matrix,CI_Hamiltonian,CI_Vectors, &
!        Full_Biorthog_Overlap,S2_Mat
      Real,Dimension(:,:),Allocatable::Overlap_Matrix,CI_Hamiltonian,CI_Vectors, &
        Full_Biorthog_Overlap,S2_Mat,CI_Values2,State_S2,Fib_Mat,Renorm_Val,Projected_UHF, &
        State_Energies
      Real,Dimension(:,:,:,:),Allocatable::ERIs,MO_ERIs
!
!     Format Statements
!
 1000 Format(1x,'Nuclear Repulsion Energy = ',F20.10,' au.')
 1010 Format(1x,'[H,S2] Fibonacci Norm = ',F10.6)
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
!     Compute Hartree-Fock Energy 
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
      Write(IOut,*)'UHF_Energy:',UHF_Energy
      Write(IOut,*)
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
!
!     Generate Slater Determinants       
!
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
        ElseIf((NAlpha_SO.gt.0).and.((NAlpha-NAlpha_SO+K)).eq.I) then
          Magnetic_Orbitals(J) = I
          J = J + 1
          K = K + 1
        EndIf
      EndDo

      Call Print_Vector(IOut,Magnetic_Orbitals,'Magnetic Orbital Pairs') 

      Do I = 1, NBasRd  
        If(Magnetic_Orbitals(I).le.NAlpha) NAlpha_in_AS = NAlpha_in_AS + 1
        If(Magnetic_Orbitals(I).le.NBeta) NBeta_in_AS = NBeta_in_AS + 1 
      EndDo
      NAlpha_in_Core = NAlpha - NAlpha_in_AS
      NBeta_in_Core = NBeta - NBeta_in_AS

      Call Gen_Det_Str(IOut,IPrint,NBasRd,NAlpha_in_AS,NBeta_in_AS,Alpha_Strings,Beta_Strings)

      NAlpha_Str = Size(Alpha_Strings,1) 
      NBeta_Str = Size(Beta_Strings,1) 
      NDMax = NAlpha_Str * NBeta_Str 
      Max_SR = (NAlpha_in_AS + NBeta_in_AS - NAlpha_SO)/2
!      If(SR_Trunc.gt.Max_SR) then
!        Write(IOut,*) "ERROR: Requested seniority truncation greater than maximum possible"
!        Call Exit()
!      EndIf
      Allocate(SR_Index(SR_Trunc+1))
      NDets = 0
      Do I =1,SR_Trunc+1
        SR_Index(I) = Factorial(NBasRd)/(Factorial(Max_SR-I+1)*Factorial(NBasRd-Max_SR+I-1)) 
        SR_Index(I) = SR_Index(I) * (Factorial(NBasRd-Max_SR+I-1) / &
          (Factorial(NAlpha_in_AS-Max_SR+I-1)*Factorial((NBasRd-Max_SR+I-1)-(NAlpha_in_AS-Max_SR+I-1)))) 
        SR_Index(I) = SR_Index(I) * (Factorial((NBasRd-Max_SR+I-1)-(NAlpha_in_AS-Max_SR+I-1)) / & 
          (Factorial(NBeta_in_AS-Max_SR+I-1)*Factorial(((NBasRd-Max_SR+I-1)- &
          (NAlpha_in_AS-Max_SR+I-1))-(NBeta_in_AS-Max_SR+I-1)))) 
        NDets = NDets + SR_Index(I) 
      EndDo
      If(NDets.eq.0.or.NBasRd.eq.0) NDets = 1
      Write(IOut,*) "Truncation is removing ",NDMax-NDets," determinants"

      Alpha_Core_String = 0
      Do I = 0, NAlpha-1
        Alpha_Core_String = IBSet(Alpha_Core_String,I)
        Do J = 1, NBasRd
          If((I+1).eq.Magnetic_Orbitals(J)) Alpha_Core_String = IBClr(Alpha_Core_String,I)
        EndDo
      EndDo
      If(NAlpha_in_Core.ne.PopCnt(Alpha_Core_String)) then
        Write(IOut,*) 'Number of bits set in Alpha_Core_String conflicts with &
          NAlpha_in_Core'
        Call Exit()
      EndIf
      
      Do I = 1, NAlpha_Str
        Temp_String = 0
        Do J = 1, NBasRd
          Call MvBits(Alpha_Strings(I,1),J-1,1,Temp_String,Magnetic_Orbitals(J)-1) 
        EndDo
        Alpha_Strings(I,1) = IOr(Temp_String,Alpha_Core_String)
      EndDo

      Beta_Core_String = 0
      Do I = 0, NBeta-1
        Beta_Core_String = IBSet(Beta_Core_String,I)
        Do J = 1, NBasRd
          If((I+1).eq.Magnetic_Orbitals(J)) Beta_Core_String = IBClr(Beta_Core_String,I)
        EndDo
      EndDo
      If(NBeta_in_Core.ne.PopCnt(Beta_Core_String)) then
        Write(IOut,*) 'Number of bits set in Beta_Core_String conflicts with &
          NBeta_in_Core'
        Call Exit()
      EndIf
      
      Do I = 1, NBeta_Str
        Temp_String = 0
        Do J = 1, NBasRd
          Call MvBits(Beta_Strings(I,1),J-1,1,Temp_String,Magnetic_Orbitals(J)-1) 
        EndDo
        Beta_Strings(I,1) = IOr(Temp_String,Beta_Core_String)
      EndDo

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
      MO_Core_Ham%Alpha = MatMul(MatMul(Transpose(Biorthog_Orbs%Alpha),Core_Hamiltonian%Alpha), &
        Biorthog_Orbs%Alpha)  
      MO_Core_Ham%Beta = MatMul(MatMul(Transpose(Biorthog_Orbs%Beta),Core_Hamiltonian%Beta), &
        Biorthog_Orbs%Beta)  
      If(IPrint.ge.2) then
        Call Print_Matrix(IOut,MO_Core_Ham%Alpha,' Alpha MO Core Hamiltonian')
        Call Print_Matrix(IOut,MO_Core_Ham%Beta,' Beta MO Core Hamiltonian')
      EndIf
      Call TwoERI_Trans(IOut,IPrint,Biorthog_Orbs,ERIs,MO_ERIs,.True.)
!
!     Generate Hamiltonian and S**2 Matrix
!
      Allocate(CI_Hamiltonian(NDets,NDets),S2_Mat(NDets,NDets))
      CI_Hamiltonian = Zero
      S2_Mat = Zero
      If(IPrint.ge.2) Call Print_Matrix(IOut,CI_Hamiltonian,'CI Hamiltonian before')
      If(IPrint.ge.2) Call Print_Matrix(IOut,S2_Mat,'S^2 Matrix before')
      L_Index = 0
      R_Index = 0
      Do L_A_String = 1, NAlpha_Str  
        Do L_B_String = 1, NBeta_Str  
          If((PopCnt(IEOr(Alpha_Strings(L_A_String,1),Beta_Strings(L_B_String,1)))-NAlpha_SO)/2.le.SR_Trunc) then
            L_Index = L_Index + 1 
            R_Index = 0
            Do R_A_String = 1, NAlpha_Str  
              Do R_B_String = 1, NBeta_Str  
                If((PopCnt(IEOr(Alpha_Strings(R_A_String,1),Beta_Strings(R_B_String,1)))-NAlpha_SO)/2.le.SR_Trunc) then
                  R_Index = R_Index + 1 
!                  Write(IOut,*) 'L_Index:',L_Index,'R_Index:',R_Index
!                  Write(IOut,*) 'L String'
!                  Write(IOut,'(2B64)') Alpha_Strings(L_A_String,1),Beta_Strings(L_B_String,1)
!                  Write(IOut,*) 'R String'
!                  Write(IOut,'(2B64)') Alpha_Strings(R_A_String,1),Beta_Strings(R_B_String,1)
                  CI_Hamiltonian(L_Index,R_Index) = Slater_Condon(IOut,IPrint,NBasis,Alpha_Strings(L_A_String,1), &
                    Beta_Strings(L_B_String,1),Alpha_Strings(R_A_String,1),Beta_Strings(R_B_String,1),MO_Core_Ham,MO_ERIs,.True.)
                  S2_Mat(L_Index,R_Index) = S2_Mat_Elem(IOut,IPrint,NBasis,Alpha_Strings(L_A_String,1), &
                    Beta_Strings(L_B_String,1),Alpha_Strings(R_A_String,1),Beta_Strings(R_B_String,1),Full_Biorthog_Overlap)
                EndIf
              EndDo
            EndDo
          EndIf
        EndDo
      EndDo
      If(IPrint.ge.1) Call Print_Matrix(IOut,CI_Hamiltonian,'CI Hamiltonian')
      If(IPrint.ge.1) Call Print_Matrix(IOut,S2_Mat,'S^2 Matrix')
!
!     Diagonalize Hamiltonian and S**2 to get simultaneous eigenvectors
!
      Write(IOut,*) 'Diagonalizing CI Hamiltonian and S**2 Matrix'
      Write(IOut,*)
      Allocate(CI_Vectors(NDets,NDets),CI_Values(NDets),Fib_Mat(NDets,NDets))
!      Call Print_Matrix(IOut,MatMul(CI_Hamiltonian,S2_Mat)-MatMul(S2_Mat,CI_Hamiltonian),'[S2,H]')
      Fib_Mat = MatMul(Transpose(MatMul(CI_Hamiltonian,S2_Mat)-MatMul(S2_Mat,CI_Hamiltonian)),MatMul(CI_Hamiltonian,S2_Mat)-MatMul(S2_Mat,CI_Hamiltonian))
      Fib_Norm = Zero
      Do I=1,NDets
        Fib_Norm = Fib_Norm + Fib_Mat(I,I)
      EndDo
      Write(IOut,1010) Fib_Norm**0.5
      Deallocate(Fib_Mat)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      Call Jade(IOut,IPrint,NDets,CI_Hamiltonian,S2_Mat,CI_Vectors,CI_Values2,State_S22)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      Do I = 1, NDets
        S2_Mat(I,I) = S2_Mat(I,I) + One
      EndDo
      Call MQC_Matrix_Diagonalize(MatMul(CI_Hamiltonian,S2_Mat),CI_Vectors,CI_Values)
!      Call MQC_Matrix_Diagonalize(CI_Hamiltonian,CI_Vectors,CI_Values)
!      Call MQC_Matrix_Diagonalize(S2_Mat,CI_Vectors,CI_Values)
!      If(IPrint.ge.1) then 
        Call Print_Matrix(IOut,CI_Vectors,'CI Eigenvectors')
        Call Print_Vector(IOut,CI_Values,'CI Eigenvalues')
!      EndIf

      Allocate(State_S2(NDets,NDets),State_Energies(NDets,NDets))
      State_S2 = MatMul(Transpose(CI_Vectors),MatMul(S2_Mat,CI_Vectors)) 
      Do I = 1, NDets
        State_S2(I,I) = State_S2(I,I) - One
      EndDo
      State_Energies = MatMul(Transpose(CI_Vectors),MatMul(CI_Hamiltonian,CI_Vectors)) 
      Call Print_Matrix(IOut,State_Energies,'Rotated H')
      Call Print_Matrix(IOut,State_S2,'Rotated S**2')

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      Call Print_Matrix(IOut,CI_Vectors,'Joint Eigenvectors')
!      Call Print_Matrix(IOut,State_S22,'S2 Eigenvalues')
!      Call Print_Matrix(IOut,CI_Values2,'CI Eigenvalues')
!      
!      Allocate(Final_Energy(NDets),Final_S2(NDets))
!      Do I=1, NDets
!        Final_Energy(I) = CI_Values2(I,I)
!        Final_S2(I) = State_S22(I,I)
!      EndDo
!      Final_Energy = Vnn + Final_Energy
!      Call Print_Vector(IOut,Final_Energy,'Biorthogonal Energy including Vnn') 
!      Call Print_Vector(IOut,Final_S2,'State S**2 values') 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Rearrange UHF equation and renormalize to get projected UHF energy
!
      Allocate(Projected_UHF(NDets,NDets),Renorm_Val(NDets,NDets))
      Projected_UHF = CI_Hamiltonian 
      Renorm_Val = Zero
      Do I = 1, NDets
        Do J = 1, NDets
          Do P = 1, NDets
            Do Q = 1, NDets
!              Write(IOut,*)'State:',P,Q,' S2:',State_S2(P,Q),' Desired:',((Multiplicity**2)/4)+0.00001
              If(State_S2(P,Q).lt.((Multiplicity**2)/4)+0.00001) then
                If(I.eq.4.and.J.eq.4) then
!                  Write(IOut,*) 'State:',P,Q
!                  Write(IOut,*) 'CI_Vectors(I,P)',CI_Vectors(I,P)
!                  Write(IOut,*) 'CI_Vectors(J,Q)',CI_Vectors(J,Q)
                EndIf
                If(P.eq.Q) then
                  Renorm_Val(I,J) = Renorm_Val(I,J) + CI_Vectors(I,P)*CI_Vectors(J,Q)
                EndIf
              Else
                Projected_UHF(I,J) = Projected_UHF(I,J) - CI_Vectors(I,P)*CI_Vectors(J,Q)*State_Energies(P,Q)
              EndIf
            EndDo
          EndDo
        EndDo
      EndDo

      Renorm_Val = 1/Renorm_Val
      If(IPrint.ge.1) Call Print_Matrix(IOut,Renorm_Val,'Renormalization Constants')
      If(IPrint.ge.1) Call Print_Matrix(IOut,Projected_UHF,'Raw Projected Value')
      Projected_UHF = Projected_UHF*Renorm_Val
      If(IPrint.ge.1) Call Print_Matrix(IOut,Projected_UHF,'Full PUHF')
      Projected_UHF = Projected_UHF + Vnn
      Call Print_Matrix(IOut,Projected_UHF,'Final Projected Matrix')
      Allocate(Final_Energy(NDets))
      ForAll(I=1:NDets) Final_Energy(I) = Projected_UHF(I,I) 
      Call Print_Vector(IOut,Final_Energy,'Final Projected Energy')
      Write(IOut,*)'Projected Ground State Energy:',MinVal(Final_Energy)
      Write(IOut,*)
!
 999  End Program truncatedAP
