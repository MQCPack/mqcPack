!>    \brief <b> MQC Molecule contains objects for assessing and manipulating  
!>    information related to the nuclei </b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC molecule contains objects that store and manipulate information related
!>    to nuclei and the nuclear structure of molecules in Fortran. The objects 
!>    enable straightforward retrieval of information such as the nucler-nuclear
!>    reulsion energy or printing of the nuclear coordinates. The MQC derived types
!>    defined in this package are
!>
!>    1. MQC_Molecule_Data:     EST object containing nuclear structure information
!>
!>    This module is level 2 in the MQC hierarchy and so depends on level 0 and 1
!>    modules. 
!>
!>    \endverbatim
!
      Module MQC_Molecule
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
!     **                          Version 25.8.0                           **
!     **                          August 1, 2025                           **
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
      Use MQC_General
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
!
!>    \brief <b> Nuclear structure information object</b>
      Type MQC_Molecule_Data
        Type(MQC_Scalar)::NAtoms
        Type(MQC_Vector)::Atomic_Numbers,Atomic_Masses,Nuclear_Charges
        Type(MQC_Matrix)::Cartesian_Coordinates
      Contains 
!>      \brief <b> Print the nuclear geometry in a pretty table</b>
        Procedure, Public::print => MQC_Print_Nuclear_Geometry
!>      \brief <b> Return the classical nuclear-nuclear repulsion energy </b>
        Procedure, Public::getNucRep => MQC_Get_Nuclear_Repulsion
!>      \brief <b> Return the number of atoms</b>
        Procedure, Public::getNumAtoms => MQC_Molecule_Data_GetNAtoms
!>      \brief <b> Update the information in the MQC_Molecule_Data object</b>
        Procedure, Public::updateMolData => MQC_Update_Mol_Data
      End Type MQC_Molecule_Data
!
!>    \brief <b> Print the nuclear geometry in a pretty table</b>
      Interface MQC_Print
        Module Procedure MQC_Print_Nuclear_Geometry
      End Interface

      CONTAINS
!
!     PROCEDURE MQC_Molecule_Data_Allocate
!
!>    \brief <b> This subroutine is used to set the NAtoms element of an
!>     MQC_Molecule_Data object and to allocate all of the object's array
!>     elements</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This subroutine is used to set the NAtoms element of an
!>    MQC_Molecule_Data object and to allocate all of the object's array
!>    elements.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[inOut] MQC_Molecule_Data_Object
!>    \verbatim
!>        MQC_Molecule_Data_Object is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data that is allocated.
!>    \endverbatim
!>
!>    \param[in] NAtomsIn
!>    \verbatim
!>        NAtomsIn is Type(MQC_Scalar) 
!>        The number of Atoms. 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author H. P. Hratchian
!>    \date 2016
!
      Subroutine MQC_Molecule_Data_Allocate(MQC_Molecule_Data_Object,  &
        NAtomsIn)
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
        Call MQC_Molecule_Data_Object%Atomic_Numbers%init(NAtoms,0)
        Call MQC_Molecule_Data_Object%Atomic_Masses%init(NAtoms,0)
        Call MQC_Molecule_Data_Object%Nuclear_Charges%init(Natoms,0)
        Call MQC_Molecule_Data_Object%Cartesian_Coordinates%init(3,NAtoms,0)
      EndIf
!
      End Subroutine MQC_Molecule_Data_Allocate
!
!
!     PROCEDURE MQC_Molecule_Data_DeAllocate
!
!>    \brief <b> This subroutine is used to deallocate all of the object's 
!>    array elements</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This subroutine is used to deallocate all of the object's array elements.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[inOut] MQC_Molecule_Data_Object
!>    \verbatim
!>        MQC_Molecule_Data_Object is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data that is deallocated.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author H. P. Hratchian
!>    \date 2016
!
      Subroutine MQC_Molecule_Data_DeAllocate(MQC_Molecule_Data_Object)
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
!
!     PROCEDURE MQC_Molecule_Data_GetNAtoms
!
!>    \brief <b> This function returns the number of atoms associated with the
!>    MQC_Molecule_Data object MQC_Molecular_Data_Object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This function returns the number of atoms associated with the
!>    MQC_Molecule_Data object MQC_Molecular_Data_Object.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[inOut] MQC_Molecule_Data_Object
!>    \verbatim
!>        MQC_Molecule_Data_Object is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data from which NAtoms is obtained.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author H. P. Hratchian
!>    \date 2016
!
      Function MQC_Molecule_Data_GetNAtoms(MQC_Molecule_Data_Object)
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
!
!>    \brief <b> This subroutine is used to put data into an MQC_Molecule_Data 
!>    object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This subroutine is used to put data into an MQC_Molecule_Data object.
!>    All dummy arguments, after the first, are optional. Note that if the
!>    number of atoms in the data changes wrt to previous data in the
!>    object, the entire object is reset.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[inOut] MQC_Molecule_Data_Object
!>    \verbatim
!>        MQC_Molecule_Data_Object is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data from which NAtoms is obtained.
!>    \endverbatim
!>
!>    \param[in] NAtomsIn
!>    \verbatim
!>        NAtomsIn is Integer,Optional
!>        The number of atoms.
!>    \endverbatim
!>
!>    \param[in] Atomic_Numbers
!>    \verbatim
!>        Atomic_Numbers is Integer,Dimension(:),Optional
!>        The list of atomic numbers.
!>    \endverbatim
!>
!>    \param[in] Atomic_Masses
!>    \verbatim
!>        Atomic_Masses is Real,Dimension(:),Optional
!>        The list of atomic masses.
!>    \endverbatim
!>
!>    \param[in] Nuclear_Charges 
!>    \verbatim
!>        Nuclear_Charges is Real,Dimension(:),Optional
!>        The list of nuclear charges.
!>    \endverbatim
!>
!>    \param[in] Cartesian_Coordinates
!>    \verbatim
!>        Cartesian_Coordinates is Real,Dimension(:,:),Optional
!>        The list of nuclear Cartesian coordinates.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author H. P. Hratchian
!>    \date 2016
!
      Subroutine MQC_Molecule_Data_Fill(MQC_Molecule_Data_Object,  &
        NAtomsIn,Atomic_Numbers,Atomic_Masses,Nuclear_Charges,  &
        Cartesian_Coordinates)
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
          Call MQC_Error_I('Error in MQC_Molecule_DataFill', 24, &
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
!     PROCEDURE MQC_Molecule_Data_Output
!
!>    \brief <b> This subroutine is used to output data from an 
!>    MQC_Molecule_Data object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This subroutine is used to output data from an MQC_Molecule_Data object.
!>    All dummy arguments, after the first, are optional.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] MQC_Molecule_Data_Object
!>    \verbatim
!>        MQC_Molecule_Data_Object is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data from which NAtoms is obtained.
!>    \endverbatim
!>
!>    \param[inOut] NAtomsIn
!>    \verbatim
!>        NAtomsIn is Integer,Optional
!>        The number of atoms.
!>    \endverbatim
!>
!>    \param[inOut] Atomic_Numbers
!>    \verbatim
!>        Atomic_Numbers is Integer,Dimension(:),Optional
!>        The list of atomic numbers.
!>    \endverbatim
!>
!>    \param[inOut] Atomic_Masses
!>    \verbatim
!>        Atomic_Masses is Real,Dimension(:),Optional
!>        The list of atomic masses.
!>    \endverbatim
!>
!>    \param[inOut] Nuclear_Charges 
!>    \verbatim
!>        Nuclear_Charges is Real,Dimension(:),Optional
!>        The list of nuclear charges.
!>    \endverbatim
!>
!>    \param[inOut] Cartesian_Coordinates
!>    \verbatim
!>        Cartesian_Coordinates is Real,Dimension(:,:),Optional
!>        The list of nuclear Cartesian coordinates.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson 
!>    \date 2020
!
      Subroutine MQC_Molecule_Data_Output(MQC_Molecule_Data_Object,  &
        NAtoms,Atomic_Numbers,Atomic_Masses,Nuclear_Charges,  &
        Cartesian_Coordinates)
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(In)::MQC_Molecule_Data_Object
      Integer(kind=int64),Optional,Allocatable,Intent(InOut)::NAtoms
      Integer(kind=int64),Dimension(:),Optional,Allocatable,Intent(InOut)::Atomic_Numbers
      Real(kind=int64),Dimension(:),Optional,Allocatable,Intent(InOut)::Atomic_Masses,  &
        Nuclear_Charges
      Real(kind=int64),Dimension(:,:),Optional,Allocatable,Intent(InOut)::Cartesian_Coordinates
!
      If(Present(NAtoms)) NAtoms = MQC_Molecule_Data_Object%NAtoms
      If(Present(Atomic_Numbers)) Atomic_Numbers = MQC_Molecule_Data_Object%Atomic_Numbers 
      If(Present(Atomic_Masses)) Atomic_Masses = MQC_Molecule_Data_Object%Atomic_Masses 
      If(Present(Nuclear_Charges)) Nuclear_Charges = MQC_Molecule_Data_Object%Nuclear_Charges
      If(Present(Cartesian_Coordinates)) Cartesian_Coordinates = &
        MQC_Molecule_Data_Object%Cartesian_Coordinates
!
      End Subroutine MQC_Molecule_Data_Output
!
!
!     PROCEDURE MQC_GET_NUCLEAR_REPULSION 
!
!>    \brief <b> This subroutine is used to output data from an 
!>    MQC_Molecule_Data object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This subroutine is used to output data from an MQC_Molecule_Data object.
!>    All dummy arguments, after the first, are optional.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Molecule_Info
!>    \verbatim
!>        Molecule_Info is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data from which nuclear repulsion is determined.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson 
!>    \date 2020
!
      Function MQC_Get_Nuclear_Repulsion(Molecule_Info) Result(Vnn)  
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
!     PROCEDURE MQC_GET_NUCLEAR_SYMBOLS
!
!>    \brief <b> This subroutine is used to output a list of nuclear symbols 
!>    for each element in the MQC_Molecule_Data object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This function is used to output a list of nuclear symbols for each element 
!>    in the MQC_Molecule_Data object.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Molecule_Info
!>    \verbatim
!>        Molecule_Info is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data from which nuclear symbols are obtained.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson 
!>    \date 2017
!
      Function MQC_Get_Nuclear_Symbols(Molecule_Info) Result(Symbols)  
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(In)::Molecule_Info
      Character(len=2),Dimension(:),allocatable::Symbols
      Integer::I
!
      Allocate(Symbols(int(Molecule_Info%NAtoms)))
      Do I = 1,int(Molecule_Info%NAtoms)
        Symbols(I) = mqc_element_symbol(molecule_info%atomic_numbers%at(i)) 
      EndDo
!
      End Function MQC_Get_Nuclear_Symbols
!
!
!     PROCEDURE MQC_PRINT_NUCLEAR_GEOMETRY
!
!>    \brief <b> This subroutine prints the nuclear geometry given the molecular 
!>    data object. Default units are Angstrom unless optional argument bohr is 
!>    passed as true</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This subroutine prints the nuclear geometry given the molecular data object. 
!>    Default units are Angstrom unless optional argument bohr is passed as true.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Molecule_Info
!>    \verbatim
!>        Molecule_Info is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data from which information is obtained.
!>    \endverbatim
!>
!>    \param[in] IOut
!>    \verbatim
!>        IOut is integer
!>        The output file number to which the the information is printed.
!>    \endverbatim
!>
!>    \param[in] Bohr
!>    \verbatim
!>        Bohr is logical,optional
!>        If true, the output coordinates are given in Bohr, otherwise, then
!>        are given in Angstrom (default).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson 
!>    \date 2019
!
      subroutine mqc_print_nuclear_geometry(molecule_info,iOut,bohr) 
!
      implicit none
      class(mqc_molecule_data),intent(in)::molecule_info
      integer,intent(in)::iOut
      logical,intent(in),optional::bohr
      integer::i,nAtoms
      character(len=2)::elem
      logical::myBohr
      real::sca,x,y,z
!
 1000 Format(1x,'          NUCLEAR CARTESIAN COORDINATES          ')
 1001 Format(1x,'-------------------------------------------------')
 1002 Format(2x,'Element              Coordinates (Ang)           ')
 1003 Format(17x,'X',13x,'Y',13x,'Z')
 1004 Format(5x,A,3F14.6)
 1005 Format(2x,'Element               Coordinates (au)           ')

      if(present(bohr)) then
        myBohr = bohr
      else
        myBohr = .false.
      endIf

      write(iOut,1000)
      write(iOut,1001)
      if(myBohr) then
        write(iOut,1005)
        sca = 1.0d0
      else
        sca = angPBohr
        write(iOut,1002)
      endIf
      write(iOut,1003)
      write(iOut,1001)

      nAtoms = molecule_info%nAtoms
      do i = 1,nAtoms
        elem = mqc_element_symbol(molecule_info%atomic_numbers%at(i)) 
        x = sca*molecule_info%cartesian_coordinates%at(1,i)
        y = sca*molecule_info%cartesian_coordinates%at(2,i)
        z = sca*molecule_info%cartesian_coordinates%at(3,i)
        write(iOut,1004) elem,x,y,z 
      endDo
!
      write(iOut,1001)
!
      end subroutine mqc_print_nuclear_geometry
!
!
!     PROCEDURE MQC_ELEMENT_SYMBOL
!
!>    \brief <b> This function returns the nuclear element string given the 
!>    atomic number</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This function returns the nuclear element string given the atomic number.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Atomic_Number
!>    \verbatim
!>        Atomic_Number is Class(*)
!>        The atomic number from which the element symbol is obtained.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson 
!>    \date 2019
!
      function mqc_element_symbol(atomic_number) result(element)  
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
      case (0)
        element = 'Bq'
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
!
!     PROCEDURE MQC_ATOMIC_NUMBER 
!
!>    \brief <b> This function returns the atomic number given the nuclear element 
!>    string</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This function returns the atomic number given the nuclear element string.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Element_Symbol
!>    \verbatim
!>        Element_Symbol is character(len=*)
!>        The element symbol from which the atomic number is obtained.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson 
!>    \date 2020
!
      function mqc_atomic_number(element_symbol) result(atomic_number)  
!
      implicit none
      character(len=*),intent(in)::element_symbol
      integer::atomic_number
      character(len=2)::temp
!
      temp = trim(element_symbol)
      call String_Change_Case(temp,'l')
      select case (temp)
      case ('bq')
        atomic_number = 0 
      case ('h')
        atomic_number = 1 
      case ('he')
        atomic_number = 2
      case ('li')
        atomic_number = 3
      case ('be')
        atomic_number = 4
      case ('b')
        atomic_number = 5
      case ('c')
        atomic_number = 6
      case ('n')
        atomic_number = 7
      case ('o')
        atomic_number = 8
      case ('f')
        atomic_number = 9
      case ('ne')
        atomic_number = 10
      case ('na')
        atomic_number = 11
      case ('mg')
        atomic_number = 12
      case ('al')
        atomic_number = 13
      case ('si')
        atomic_number = 14
      case ('p')
        atomic_number = 15
      case ('s')
        atomic_number = 16
      case ('cl')
        atomic_number = 17
      case ('ar')
        atomic_number = 18
      case ('k')
        atomic_number = 19
      case ('ca')
        atomic_number = 20
      case ('sc')
        atomic_number = 21
      case ('ti')
        atomic_number = 22
      case ('v')
        atomic_number = 23
      case ('cr')
        atomic_number = 24
      case ('mn')
        atomic_number = 25
      case ('fe')
        atomic_number = 26
      case ('co')
        atomic_number = 27
      case ('ni')
        atomic_number = 28
      case ('cu')
        atomic_number = 29
      case ('zn')
        atomic_number = 30
      case ('ga')
        atomic_number = 31
      case ('ge')
        atomic_number = 32
      case ('as')
        atomic_number = 33
      case ('se')
        atomic_number = 34
      case ('br')
        atomic_number = 35
      case ('kr')
        atomic_number = 36
      case ('rb')
        atomic_number = 37
      case ('sr')
        atomic_number = 38
      case ('y')
        atomic_number = 39
      case ('zr')
        atomic_number = 40
      case ('nb')
        atomic_number = 41
      case ('mo')
        atomic_number = 42
      case ('tc')
        atomic_number = 43
      case ('ru')
        atomic_number = 44
      case ('rh')
        atomic_number = 45
      case ('pd')
        atomic_number = 46
      case ('ag')
        atomic_number = 47
      case ('cd')
        atomic_number = 48
      case ('in')
        atomic_number = 49
      case ('sn')
        atomic_number = 50
      case ('sb')
        atomic_number = 51
      case ('te')
        atomic_number = 52
      case ('i')
        atomic_number = 53
      case ('xe')
        atomic_number = 54
      case ('cs')
        atomic_number = 55
      case ('ba')
        atomic_number = 56
      case ('la')
        atomic_number = 57
      case ('ce')
        atomic_number = 58
      case ('pr')
        atomic_number = 59
      case ('nd')
        atomic_number = 60
      case ('pm')
        atomic_number = 61
      case ('sm')
        atomic_number = 62
      case ('eu')
        atomic_number = 63
      case ('gd')
        atomic_number = 64
      case ('tb')
        atomic_number = 65
      case ('dy')
        atomic_number = 66
      case ('ho')
        atomic_number = 67
      case ('er')
        atomic_number = 68
      case ('tm')
        atomic_number = 69
      case ('yb')
        atomic_number = 70
      case ('lu')
        atomic_number = 71
      case ('hf')
        atomic_number = 72
      case ('ta')
        atomic_number = 73
      case ('w')
        atomic_number = 74
      case ('re')
        atomic_number = 75
      case ('os')
        atomic_number = 76
      case ('ir')
        atomic_number = 77
      case ('pt')
        atomic_number = 78
      case ('au')
        atomic_number = 79
      case ('hg')
        atomic_number = 80
      case ('tl')
        atomic_number = 81
      case ('pb')
        atomic_number = 82
      case ('bi')
        atomic_number = 83
      case ('po')
        atomic_number = 84
      case ('at')
        atomic_number = 85
      case ('rn')
        atomic_number = 86
      case ('fr')
        atomic_number = 87
      case ('ra')
        atomic_number = 88
      case ('ac')
        atomic_number = 89
      case ('th')
        atomic_number = 90
      case ('pa')
        atomic_number = 91
      case ('u')
        atomic_number = 92
      case ('np')
        atomic_number = 93
      case ('pu')
        atomic_number = 94
      case ('am')
        atomic_number = 95
      case ('cm')
        atomic_number = 96
      case ('bk')
        atomic_number = 97
      case ('cf')
        atomic_number = 98
      case ('es')
        atomic_number = 99
      case ('fm')
        atomic_number = 100
      case ('md')
        atomic_number = 101
      case ('no')
        atomic_number = 102
      case ('lr')
        atomic_number = 103
      case ('rf')
        atomic_number = 104
      case ('db')
        atomic_number = 105
      case ('sg')
        atomic_number = 106
      case ('bh')
        atomic_number = 107
      case ('hs')
        atomic_number = 108
      case ('mt')
        atomic_number = 109
      case ('ds')
        atomic_number = 110
      case ('rg')
        atomic_number = 111
      case default
        call mqc_error_a('Unknown element symbol in mqc_atomic_number',6,'temp',temp)
      end select
!
      end function mqc_atomic_number 
!
!
!     PROCEDURE MQC_UPDATE_MOL_DATA
!
!>    \brief <b> This subroutine updates the molecule data object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    This subroutine updates the molecule data object. Only one of the optional
!>    MQC_Scalar, MQC_Vector, or MQC_Matrix arrays should be provided. The 
!>    recognized labels and their meaning include:
!>    1. Label = 'natoms'          update the number of atoms.
!>    2. Label = 'atomic numbers'  update the atomic number vector.
!>    3. Label = 'atomic masses'   update the atomic mass vector.
!>    4. Label = 'nuclear charges' update the nuclear charge vector.
!>    5. Label = 'cartesians'      update the cartesian coordinate matrix.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[inOut] Molecule_Info
!>    \verbatim
!>        Molecule_Info is Class(MQC_Molecule_Data)
!>        The MQC_Molecule_Data that will be updated.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is character(len=*)
!>        The label indicating the data to update.
!>        = 'natoms'          update the number of atoms.
!>        = 'atomic numbers'  update the atomic number vector.
!>        = 'atomic masses'   update the atomic mass vector.
!>        = 'nuclear charges' update the nuclear charge vector.
!>        = 'cartesians'      update the cartesian coordinate matrix.
!>    \endverbatim
!>
!>    \param[in] Scalar
!>    \verbatim
!>        Scalar is Type(MQC_Scalar),Optional
!>        The name of the MQC_Scalar variable.
!>    \endverbatim
!>
!>    \param[in] Vector
!>    \verbatim
!>        Vector is Type(MQC_Vector),Optional
!>        The name of the MQC_Vector variable.
!>    \endverbatim
!>
!>    \param[in] Matrix
!>    \verbatim
!>        Matrix is Type(MQC_Matrix),Optional
!>        The name of the MQC_Matrix variable.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author L. M. Thompson 
!>    \date 2020
!
      Subroutine MQC_Update_Mol_Data(Molecule_Info,label,scalar,vector,matrix) 
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(InOut)::Molecule_Info
      character(len=*),intent(in)::label
      Type(MQC_Scalar),optional,intent(in)::scalar
      Type(MQC_Vector),optional,intent(in)::vector
      Type(MQC_Matrix),optional,intent(in)::matrix
      character(len=64)::myLabel
      Integer::nInputArrays
!
      nInputArrays = 0
      if(present(scalar)) nInputArrays = nInputArrays + 1
      if(present(vector)) nInputArrays = nInputArrays + 1
      if(present(matrix)) nInputArrays = nInputArrays + 1
      if(nInputArrays.ne.1) call mqc_error_i('Too many input values sent to mqc_update_mol_data',&
        6,'nInputArrays',nInputArrays)
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('natoms')
        if(.not.present(scalar)) call mqc_error('Update natoms requested in mqc_update_mol_data&
          & but input data is not scalar type')
        molecule_info%natoms = scalar
      case('atomic numbers')
        if(.not.present(vector)) call mqc_error('Update atomic numbers requested in mqc_update_mol_data&
          & but input data is not vector type')
        molecule_info%atomic_numbers = vector
      case('atomic masses')
        if(.not.present(vector)) call mqc_error('Update atomic masses requested in mqc_update_mol_data&
          & but input data is not vector type')
        molecule_info%atomic_masses = vector
      case('nuclear charges')
        if(.not.present(vector)) call mqc_error('Update nuclear charges requested in mqc_update_mol_data&
          & but input data is not vector type')
        molecule_info%nuclear_charges = vector
      case('cartesians')
        if(.not.present(matrix)) call mqc_error('Update Cartesian matrix requested in mqc_update_mol_data&
          & but input data is not matrix type')
        molecule_info%cartesian_coordinates = matrix
      case default
        call mqc_error_A('Invalid label sent to mqc_update_mol_data.', 6, &
             'mylabel', mylabel )
      end select
!
      End Subroutine MQC_Update_Mol_Data
!
      End Module MQC_Molecule
