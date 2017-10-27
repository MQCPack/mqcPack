      Program FullCI
!
!     This program perfoms a full CI calculation.
!
!     L. M. Thompson, 2016
!
      Use MQC_General
      Use MQC_Algebra
      Use MQC_Gaussian
      Use MQC_EST 
!
!     Variable Declarations...
!
      Implicit None
      Character(Len=256)::FileName
      Integer::IOut=6,IPrint=0
      Logical::UHF
      Type(MQC_Scalar)::Vnn
      Type(MQC_Vector)::Final_Energy
      Type(MQC_Matrix)::CI_Hamiltonian
      Type(MQC_Molecule_Data)::MoleculeInfo
      Type(MQC_PSCF_Wavefunction)::Wavefunction_Data
      Type(MQC_Core_Hamiltonian)::MO_Core_Ham
      Type(MQC_Determinant)::Determinants
      Type(MQC_TwoERIs)::ERIs,MO_ERIs
!
      Write(IOut,*)
      Write(IOut,*) 'Full Configuration Interaction Energy Calculator'
      Write(IOut,*)
      Write(IOut,*) 'L. M. Thompson 2017'
      Write(IOut,*)
!
!     Get the user-defined filename from the command line and then call the
!     routine that reads the Gaussian matrix file.
!
      Call Get_Command_Argument(1,FileName)
      Call MQC_Gaussian_Read_Matrix_File(.False.,FileName,MoleculeInfo, &
        Wavefunction=Wavefunction_Data,ERIs_Full=ERIs)

      If(IPrint.ge.2) then  
        Write(*,*)' Back from Read_Gaussian_Matrix_File'
        Call MQC_Print(IOut,MoleculeInfo%NAtoms,'NAtoms')
        Call MQC_Print(IOut,Wavefunction_Data%NBasis,'NBasis')
        Call MQC_Print(IOut,Wavefunction_Data%NElectrons,'NElectrons')
        Write(*,*)' WF_Type      = ',Wavefunction_Data%WF_Type
        Write(*,*)' WF_Complex   = ',Wavefunction_Data%WF_Complex
        Call MQC_Print(IOut,Wavefunction_Data%Multiplicity,'Multiplicity')
        Call MQC_Print(IOut,Wavefunction_Data%Charge,'Charge')
        Call MQC_Print(IOut,MoleculeInfo%Atomic_Numbers,'Atomic Numbers')
        Call MQC_Print(IOut,MoleculeInfo%Atomic_Masses,'Atomic Masses')
        Call MQC_Print(IOut,MoleculeInfo%Nuclear_Charges,'Nuclear Charges')
        Call MQC_Print(IOut,MoleculeInfo%Cartesian_Coordinates,'Cartesian Coordinates')
        Call MQC_Print(IOut,Wavefunction_Data%Overlap_Matrix%Elements,'Overlap Matrix')
        Call MQC_Print(IOut,Wavefunction_Data%Core_Hamiltonian%Alpha,'Core Hamiltonian Alpha')
        Call MQC_Print(IOut,Wavefunction_Data%Core_Hamiltonian%Beta,'Core Hamiltonian Beta')
        Call MQC_Print(IOut,Wavefunction_Data%Fock_Matrix%Alpha,'Alpha Fock Matrix')
        Call MQC_Print(IOut,Wavefunction_Data%Fock_Matrix%Beta,'Beta Fock Matrix')
        Call MQC_Print(IOut,Wavefunction_Data%MO_Coefficients%Alpha,'Alpha MO Coefficient Matrix')
        Call MQC_Print(IOut,Wavefunction_Data%MO_Coefficients%Beta,'Beta MO Coefficient Matrix')
        Call MQC_Print(IOut,ERIs%TwoERIs,'ERIs')
      EndIf
!
      If(Wavefunction_Data%WF_Type.eq.'U') then
        UHF = .True.
      ElseIf(Wavefunction_Data%WF_Type.eq.'R') then
        UHF = .False.
      Else
        Call MQC_Error('Unsupported wavefunction type in fullci')
      EndIf 
!
      If (Wavefunction_Data%WF_Complex) Call MQC_Error('Complex wavefunctions unsupported in fullci')
!
!     Compute the nuclear-nuclear repulsion energy.
!
      Vnn = MQC_Get_Nuclear_Repulsion(IOut,MoleculeInfo)
!
!     Generate Slater Determinants       
!
      Call Gen_Det_Str(IOut,IPrint,Wavefunction_Data%NBasis,Wavefunction_Data%NAlpha, &
        Wavefunction_Data%NBeta,Determinants)
!
!     Transform one and two-electron integrals to MO basis
!
      If(IPrint.eq.1) Write(IOut,*) 'Transforming MO integrals'
      MO_Core_Ham%Alpha = Transpose(Wavefunction_Data%MO_Coefficients%Alpha) .dot. & 
        Wavefunction_Data%Core_Hamiltonian%Alpha .dot. Wavefunction_Data%MO_Coefficients%Alpha  
      If(UHF) MO_Core_Ham%Beta = Transpose(Wavefunction_Data%MO_Coefficients%Beta) .dot. & 
        Wavefunction_Data%Core_Hamiltonian%Beta .dot. Wavefunction_Data%MO_Coefficients%Beta  
      If(IPrint.ge.2) Call MQC_Print(IOut,MO_Core_Ham%Alpha,' Alpha MO Core Hamiltonian')
      If(IPrint.ge.2.and.UHF) Call MQC_Print(IOut,MO_Core_Ham%Beta, &
        ' Beta MO Core Hamiltonian')
      Call TwoERI_Trans(IOut,IPrint,Wavefunction_Data%MO_Coefficients,ERIs,MO_ERIs,UHF)
!
!     Generate Hamiltonian Matrix
!
      Call MQC_Build_CI_Hamiltonian(IOut,IPrint,Wavefunction_Data%NBasis,Determinants, &
        MO_Core_Ham,MO_ERIs,UHF,CI_Hamiltonian)
!
!     Diagonalize Hamiltonian
!
      If(IPrint.eq.1) Write(IOut,*) 'Diagonalizing CI Hamiltonian'
      Call MQC_Matrix_Diagonalize(CI_Hamiltonian,Wavefunction_Data%PSCF_Amplitudes,Wavefunction_Data%PSCF_Energies)
      If(IPrint.ge.1) then 
        Call MQC_Print(IOut,Wavefunction_Data%PSCF_Amplitudes,'CI Eigenvectors')
        Call MQC_Print(IOut,Wavefunction_Data%PSCF_Energies,'CI Eigenvalues')
      EndIf
!
      Final_Energy = Vnn + Wavefunction_Data%PSCF_Energies
      Call MQC_Print(IOut,Final_Energy,'Final Energy') 
!
 999  End Program FullCI     
