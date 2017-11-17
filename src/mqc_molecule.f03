      Module MQC_Molecule
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
      Use MQC_General
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
!
!     MQC_Molecule_Data
      Type MQC_Molecule_Data
        Type(MQC_Scalar)::NAtoms
        Type(MQC_Vector)::Atomic_Numbers,Atomic_Masses,Nuclear_Charges
        Type(MQC_Matrix)::Cartesian_Coordinates
      End Type MQC_Molecule_Data
!

      CONTAINS
!
!     PROCEDURE MQC_Molecule_Data_Allocate
      Subroutine MQC_Molecule_Data_Allocate(MQC_Molecule_Data_Object,  &
        NAtomsIn)
!
!     This subroutine is used to set the NAtoms element of an
!     MQC_Molecule_Data object and to allocate all of the object's array
!     elements.
!
!     -H. P. Hratchian, 2016
!
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(InOut)::MQC_Molecule_Data_Object
      Type(MQC_Scalar),Intent(In)::NAtomsIn
      Integer::NAtoms
!
      If(MQC_Molecule_Data_Object%NAtoms.ne.NAtomsIn) then
        Call MQC_Molecule_Data_DeAllocate(MQC_Molecule_Data_Object)
        MQC_Molecule_Data_Object%NAtoms = NAtomsIn
        NAtoms = NAtomsIn
        Call MQC_Molecule_Data_Object%Atomic_Numbers%initialize(NAtoms,0)
        Call MQC_Molecule_Data_Object%Atomic_Masses%initialize(NAtoms,0)
        Call MQC_Molecule_Data_Object%Nuclear_Charges%initialize(Natoms,0)
        Call MQC_Molecule_Data_Object%Cartesian_Coordinates%initialize(3,NAtoms,0)
      EndIf
!
      End Subroutine MQC_Molecule_Data_Allocate


!
!     PROCEDURE MQC_Molecule_Data_DeAllocate
      Subroutine MQC_Molecule_Data_DeAllocate(MQC_Molecule_Data_Object)
!
!     This subroutine is used to set the NAtoms element of an
!     MQC_Molecule_Data object and to deallocate all of the object's array
!     elements.
!
!     -H. P. Hratchian, 2016
!
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(InOut)::MQC_Molecule_Data_Object
!
      MQC_Molecule_Data_Object%NAtoms = 0
      Call MQC_Deallocate_Vector(MQC_Molecule_Data_Object%Atomic_Numbers)
      Call MQC_Deallocate_Vector(MQC_Molecule_Data_Object%Atomic_Masses)
      Call MQC_Deallocate_Vector(MQC_Molecule_Data_Object%Nuclear_Charges)
      Call MQC_Deallocate_Matrix(MQC_Molecule_Data_Object%Cartesian_Coordinates)
!
      End Subroutine MQC_Molecule_Data_DeAllocate


!
!     PROCEDURE MQC_Molecule_Data_GetNAtoms
      Function MQC_Molecule_Data_GetNAtoms(MQC_Molecule_Data_Object)
!
!     This function returns the number of atoms associated with the
!     MQC_Molecule_Data object MQC_Molecular_Data_Object.
!
!     -H. P. Hratchian, 2016
!
!
      Implicit None
      Type(MQC_Scalar)::MQC_Molecule_Data_GetNAtoms
      Class(MQC_Molecule_Data),Intent(InOut)::MQC_Molecule_Data_Object
!
      If(.not.MQC_Scalar_IsAllocated(MQC_Molecule_Data_Object%NAtoms)) & 
        MQC_Molecule_Data_Object%NAtoms = 0
      MQC_Molecule_Data_GetNAtoms = MQC_Molecule_Data_Object%NAtoms
      End Function MQC_Molecule_Data_GetNAtoms
!
!
!     PROCEDURE MQC_Molecule_Data_Fill
      Subroutine MQC_Molecule_Data_Fill(MQC_Molecule_Data_Object,  &
        NAtomsIn,Atomic_Numbers,Atomic_Masses,Nuclear_Charges,  &
        Cartesian_Coordinates)
!
!     This subroutine is used to put data into an MQC_Molecule_Data object.
!     All dummy arguments, after the first, are optional. Note that if the
!     number of atoms in the data changes wrt to previous data in the
!     object, the entire object is reset.
!
!     -H. P. Hratchian, 2016
!
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(InOut)::MQC_Molecule_Data_Object
      Integer,Optional,Intent(In)::NAtomsIn
      Integer,Dimension(:),Optional,Intent(In)::Atomic_Numbers
      Real,Dimension(:),Optional,Intent(In)::Atomic_Masses,  &
        Nuclear_Charges
      Real,Dimension(:,:),Optional,Intent(In)::Cartesian_Coordinates
!
      Type(MQC_Scalar)::NAtoms_Object,NAtoms_Sent,Zero
      Integer::NAtoms
!
!     Get NAtoms_Object and NAtoms_Sent set. Then allocate the object if
!     the two values do NOT agree or NAtoms_Object=0.
!

      NAtoms_Object =  &
        MQC_Molecule_Data_GetNAtoms(MQC_Molecule_Data_Object)
      If(Present(NAtomsIn)) NAtoms_Sent = NAtomsIn
      If(Present(Atomic_Numbers)) NAtoms_Sent = Size(Atomic_Numbers)
      If(Present(Atomic_Masses)) NAtoms_Sent = Size(Atomic_Masses)
      If(Present(Nuclear_Charges)) NAtoms_Sent = Size(Nuclear_Charges)
      If(Present(Cartesian_Coordinates)) NAtoms_Sent =  &
        Size(cartesian_coordinates,2)
      Zero = 0
      If(NAtoms_Object.le.Zero.or.NAtoms_Object.ne.NAtoms_Sent) & 
        Call MQC_Molecule_Data_Allocate(MQC_Molecule_Data_Object,  &
        NAtoms_Sent)
!
!     Now, fill the elements of the object corresponding to the data sent
!     to this routine. Along the way, we are careful to error check and
!     ensure the number of atoms matches that of the object.
!
      NAtoms = NAtoms_Sent
      If(Present(NAtomsIn)) then
        If(NAtomsIn.ne.NAtoms)  &
          Call MQC_Error_I('Error in MQC_Molecule_Data_Fill', 6, &
          'NAtomsIn', NAtomsIn, &
          'NAtoms', NAtoms )
      endIf
      If(Present(Atomic_Numbers)) then
        If(NAtoms.ne.Size(Atomic_Numbers))  &
          Call MQC_Error_I('Error in MQC_Molecule_Data_Fill', 6, &
          'NAtoms', NAtoms, &
          'Size(Atomic_Numbers)', Size(Atomic_Numbers))
        MQC_Molecule_Data_Object%Atomic_Numbers = Atomic_Numbers
      endIf
      If(Present(Atomic_Masses)) then
        If(NAtoms.ne.Size(Atomic_Masses))  &
          Call MQC_Error_I('Error in MQC_Molecule_Data_Fill', 6, &
          'NAtoms', NAtoms, &
          'Size(Atomic_Masses)', Size(Atomic_Masses) )
        MQC_Molecule_Data_Object%Atomic_Masses = Atomic_Masses
      endIf
      If(Present(Nuclear_Charges)) then
        If(NAtoms.ne.Size(Nuclear_Charges))  &
          Call MQC_Error_I('Error in MQC_Molecule_Data_Fill', 6, &
          'NAtoms', NAtoms, &
          'Size(Nuclear_Charges)', Size(Nuclear_Charges))
        MQC_Molecule_Data_Object%Nuclear_Charges = Nuclear_Charges
      endIf
      If(Present(Cartesian_Coordinates)) then
        If(NAtoms.ne.Size(Cartesian_Coordinates,2))  &
          Call MQC_Error_I('Error in MQC_Molecule_Data_Fill', 6, &
          'NAtoms', NAtoms, &
          'Size(Cartesian_Coordinates,2)', Size(Cartesian_Coordinates,2))
        MQC_Molecule_Data_Object%Cartesian_Coordinates =  &
          Cartesian_Coordinates
      endIf
!
      End Subroutine MQC_Molecule_Data_Fill
!
!
!     PROCEDURE MQC_GET_NUCLEAR_REPULSION 
      Function MQC_Get_Nuclear_Repulsion(IOut,Molecule_Info) Result(Vnn)  
!
!     This function returns the nuclear repulsion energy given the 
!     molecular data object.
!
!     Lee M. Thompson, 2017.
!
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(In)::Molecule_Info
      Type(MQC_Scalar)::Vnn
      Type(MQC_Vector)::Temp_3Vector
      Real::Zero=0.0d0
      Integer::IOut,I,J,NAtoms
!
      Vnn = Zero 
      NAtoms = Molecule_Info%NAtoms
      Do I = 1,NAtoms-1
        Do J = I+1,NAtoms
          Temp_3Vector = Molecule_Info%Cartesian_Coordinates%vat([1,3],[I]) -  &
            Molecule_Info%Cartesian_Coordinates%vat([1,3],[J])
          Vnn = Vnn +  &
            Molecule_Info%Nuclear_Charges%at(I) *  &
            Molecule_Info%Nuclear_Charges%at(J) /  &
            SQRT(Transpose(Temp_3Vector).dot.Temp_3Vector)
        EndDo
      EndDo
!
      Call MQC_Print(Vnn,IOut,'Nuclear Repulsion Energy (au)')
!
      End Function MQC_Get_Nuclear_Repulsion
!
      End Module MQC_Molecule
