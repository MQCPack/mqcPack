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
      Contains 
        Procedure, Public::print => MQC_Print_Nuclear_Geometry
        Procedure, Public::getNucRep => MQC_Get_Nuclear_Repulsion
      End Type MQC_Molecule_Data
!
      Interface MQC_Print
        Module Procedure MQC_Print_Nuclear_Geometry
      End Interface

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
      Function MQC_Get_Nuclear_Repulsion(Molecule_Info) Result(Vnn)  
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
      Integer::I,J,NAtoms
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
      End Function MQC_Get_Nuclear_Repulsion
!
!
!     PROCEDURE MQC_PRINT_NUCLEAR_GEOMETRY
      subroutine mqc_print_nuclear_geometry(molecule_info,iOut) 
!
!     This subroutine prints the nuclear geometry given the 
!     molecular data object.
!
!     Lee M. Thompson, 2019.
!
!
      implicit none
      class(mqc_molecule_data),intent(in)::molecule_info
      integer,intent(in)::iOut
      integer::i,nAtoms
      character(len=2)::elem
!
 1000 Format(1x,'          NUCLEAR CARTESIAN COORDINATES          ')
 1001 Format(1x,'-------------------------------------------------')
 1002 Format(2x,'Element                 Coordinates             ')
 1003 Format(17x,'X',13x,'Y',13x,'Z')
 1004 Format(5x,A,3F14.6)

      write(iOut,1000)
      write(iOut,1001)
      write(iOut,1002)
      write(iOut,1003)
      write(iOut,1001)

      nAtoms = molecule_info%nAtoms
      do i = 1,nAtoms
        elem = mqc_element_symbol(molecule_info%atomic_numbers%at(i)) 
        write(iOut,1004) elem, &
          MQC_Scalar_Get_Intrinsic_Real(molecule_info%cartesian_coordinates%at(1,i)), &
          MQC_Scalar_Get_Intrinsic_Real(molecule_info%cartesian_coordinates%at(2,i)), &
          MQC_Scalar_Get_Intrinsic_Real(molecule_info%cartesian_coordinates%at(3,i))
      endDo
!
      write(iOut,1001)
!
      end subroutine mqc_print_nuclear_geometry
!
!
!     PROCEDURE MQC_ELEMENT_SYMBOL
      function mqc_element_symbol(atomic_number) result(element)  
!
!     This function returns the nuclear element string given the 
!     atomic number.
!
!     Lee M. Thompson, 2019.
!
!
      implicit none
      class(*),intent(in)::atomic_number
      character(len=2)::element
      integer::temp
!
      select type(atomic_number)
      type is (integer)
        temp = atomic_number
      type is (mqc_scalar)
        temp = atomic_number
      class default
        call mqc_error_i('atomic number type unknown in mqc_element_symbol',6)
      end select
!
      select case (temp)
      case (1)
        element = 'H '
      case (2)
        element = 'He'
      case (3)
        element = 'Li'
      case (4)
        element = 'Be'
      case (5)
        element = 'B '
      case (6)
        element = 'C '
      case (7)
        element = 'N '
      case (8)
        element = 'O '
      case (9)
        element = 'F '
      case (10)
        element = 'Ne'
      case (11)
        element = 'Na'
      case (12)
        element = 'Mg'
      case (13)
        element = 'Al'
      case (14)
        element = 'Si'
      case (15)
        element = 'P '
      case (16)
        element = 'S '
      case (17)
        element = 'Cl'
      case (18)
        element = 'Ar'
      case (19)
        element = 'K '
      case (20)
        element = 'Ca'
      case (21)
        element = 'Sc'
      case (22)
        element = 'Ti'
      case (23)
        element = 'V '
      case (24)
        element = 'Cr'
      case (25)
        element = 'Mn'
      case (26)
        element = 'Fe'
      case (27)
        element = 'Co'
      case (28)
        element = 'Ni'
      case (29)
        element = 'Cu'
      case (30)
        element = 'Zn'
      case (31)
        element = 'Ga'
      case (32)
        element = 'Ge'
      case (33)
        element = 'As'
      case (34)
        element = 'Se'
      case (35)
        element = 'Br'
      case (36)
        element = 'Kr'
      case (37)
        element = 'Rb'
      case (38)
        element = 'Sr'
      case (39)
        element = 'Y '
      case (40)
        element = 'Zr'
      case (41)
        element = 'Nb'
      case (42)
        element = 'Mo'
      case (43)
        element = 'Tc'
      case (44)
        element = 'Ru'
      case (45)
        element = 'Rh'
      case (46)
        element = 'Pd'
      case (47)
        element = 'Ag'
      case (48)
        element = 'Cd'
      case (49)
        element = 'In'
      case (50)
        element = 'Sn'
      case (51)
        element = 'Sb'
      case (52)
        element = 'Te'
      case (53)
        element = 'I '
      case (54)
        element = 'Xe'
      case (55)
        element = 'Cs'
      case (56)
        element = 'Ba'
      case (57)
        element = 'La'
      case (58)
        element = 'Ce'
      case (59)
        element = 'Pr'
      case (60)
        element = 'Nd'
      case (61)
        element = 'Pm'
      case (62)
        element = 'Sm'
      case (63)
        element = 'Eu'
      case (64)
        element = 'Gd'
      case (65)
        element = 'Tb'
      case (66)
        element = 'Dy'
      case (67)
        element = 'Ho'
      case (68)
        element = 'Er'
      case (69)
        element = 'Tm'
      case (70)
        element = 'Yb'
      case (71)
        element = 'Lu'
      case (72)
        element = 'Hf'
      case (73)
        element = 'Ta'
      case (74)
        element = 'W '
      case (75)
        element = 'Re'
      case (76)
        element = 'Os'
      case (77)
        element = 'Ir'
      case (78)
        element = 'Pt'
      case (79)
        element = 'Au'
      case (80)
        element = 'Hg'
      case (81)
        element = 'Tl'
      case (82)
        element = 'Pb'
      case (83)
        element = 'Bi'
      case (84)
        element = 'Po'
      case (85)
        element = 'At'
      case (86)
        element = 'Rn'
      case (87)
        element = 'Fr'
      case (88)
        element = 'Ra'
      case (89)
        element = 'Ac'
      case (90)
        element = 'Th'
      case (91)
        element = 'Pa'
      case (92)
        element = 'U '
      case (93)
        element = 'Np'
      case (94)
        element = 'Pu'
      case (95)
        element = 'Am'
      case (96)
        element = 'Cm'
      case (97)
        element = 'Bk'
      case (98)
        element = 'Cf'
      case (99)
        element = 'Es'
      case (100)
        element = 'Fm'
      case (101)
        element = 'Md'
      case (102)
        element = 'No'
      case (103)
        element = 'Lr'
      case (104)
        element = 'Rf'
      case (105)
        element = 'Db'
      case (106)
        element = 'Sg'
      case (107)
        element = 'Bh'
      case (108)
        element = 'Hs'
      case (109)
        element = 'Mt'
      case (110)
        element = 'Ds'
      case (111)
        element = 'Rg'
      case default
        call mqc_error_i('Unknown atomic number in mqc_element_symbol',6,'temp',temp)
      end select
!
      end function mqc_element_symbol
!
      End Module MQC_Molecule
