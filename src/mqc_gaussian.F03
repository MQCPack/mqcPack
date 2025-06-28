      Module MQC_Gaussian
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
!     **                          Version 25.6.1                           **
!     **                           June 2, 2025                            **
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
!
!     Set-up USE Association with Key MQC Modules.
!
      USE MQC_General
      USE MQC_Files
      USE MQC_Algebra
      USE MQC_Algebra2
      USE MQC_EST
      USE MQC_molecule
      USE MQC_matwrapper
      USE iso_fortran_env
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     File types...
!
!     MQC_Gaussian_FChk_File
      Type,Extends(MQC_Text_FileInfo)::MQC_Gaussian_FChk_File
        Character(Len=72)::Title
        Character(Len=10)::JobType
        Character(Len=30)::Method,BasisSet
      Contains
        Procedure::OpenFile => MQC_Gaussian_FChk_Open
      End Type MQC_Gaussian_FChk_File
!
!
!     MQC_Gaussian_Unformatted_Matrix_File
!       This object extends and belongs to the Class MQC_FileInfo.
!       Specifically, it is intended for use with Gaussian unformatted
!       matrix files.
      Type,Extends(MQC_FileInfo)::MQC_Gaussian_Unformatted_Matrix_File
        logical::declared=.false.,header_read=.false.,header_written=.false.,  &
          gaussianScalars_read=.false.
        character(len=1)::readWriteMode=' '
        character(len=64)::LabFil=' ',GVers=' ',Title=' '
        integer,allocatable::natoms,nbasis,nbasisUse,icharge,multiplicity,nelectrons,icgu, &
          NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot
        integer,dimension(:),allocatable::atomicNumbers,atomTypes,basisFunction2Atom, &
          IBasisFunctionType,IgaussianScalars
        real,dimension(:),allocatable::atomicCharges,atomicWeights,cartesians,gaussianScalars
      Contains
        procedure,pass::OpenFile         => MQC_Gaussian_Unformatted_Matrix_Open
        procedure,pass::CloseFile        => MQC_Gaussian_Unformatted_Matrix_Close
        procedure,pass::load             => MQC_Gaussian_Unformatted_Matrix_Read_Header
        procedure,pass::create           => MQC_Gaussian_Unformatted_Matrix_Write_Header
        procedure,pass::isRestricted     => MQC_Gaussian_IsRestricted
        procedure,pass::isUnrestricted   => MQC_Gaussian_IsUnrestricted
        procedure,pass::isGeneral        => MQC_Gaussian_IsGeneral
        procedure,pass::isComplex        => MQC_Gaussian_IsComplex
        procedure,pass::getAtomicNumbers => MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Numbers
        procedure,pass::getAtomCarts     => MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Carts
        procedure,pass::getAtomWeights   => MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Weights
        procedure,pass::getVal           => MQC_Gaussian_Unformatted_Matrix_Get_Value_Integer
        procedure,pass::getValReal       => MQC_Gaussian_Unformatted_Matrix_Get_Value_Real
        procedure,pass::getArray         => MQC_Gaussian_Unformatted_Matrix_Read_Array
        procedure,pass::getAtomInfo      => MQC_Gaussian_Unformatted_Matrix_Get_Atom_Info
        procedure,pass::getBasisInfo     => MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Element
        procedure,pass::getBasisArray    => MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Array
        procedure,pass::getBasisData     => MQC_Gaussian_Unformatted_Matrix_Get_Basis_Set_Matrices
        procedure,pass::getMolData       => MQC_Gaussian_Unformatted_Matrix_Get_Molecule_Data
        procedure,pass::getESTObj        => MQC_Gaussian_Unformatted_Matrix_Get_EST_Object
        procedure,private::get2ERIObj    => MQC_Gaussian_Unformatted_Matrix_Get_twoERIs    
        procedure,private::get2ERISet    => MQC_Gaussian_Unformatted_Matrix_Get_twoERISet  
        generic,public::get2ERIs         => get2ERIObj, get2ERISet
        procedure,pass::writeArray2      => MQC_Gaussian_Unformatted_Matrix_Write_Array2
        procedure,pass::writeArray       => MQC_Gaussian_Unformatted_Matrix_Write_Array
        procedure,pass::writeBasisData   => MQC_Gaussian_Unformatted_Matrix_Write_Basis_Set_Matrices
        procedure,pass::writeESTObj      => MQC_Gaussian_Unformatted_Matrix_Write_EST_Object
        procedure,private::write2ERIObj  => MQC_Gaussian_Unformatted_Matrix_Write_twoERIs
        procedure,private::write2ERISet  => MQC_Gaussian_Unformatted_Matrix_Write_twoERISet
        generic,public::write2ERIs       => write2ERIObj, write2ERISet
        procedure,pass::updateHeader     => MQC_Gaussian_Unformatted_Matrix_Update_Header_Info
      End Type MQC_Gaussian_Unformatted_Matrix_File
!
!
!     Data containers...
!
!     MQC_Gaussian_Molecule_Data
      Type,Extends(MQC_Molecule_Data)::MQC_Gaussian_Molecule_Data
        Type(MQC_Scalar)::Charge,Multiplicity
      End Type MQC_Gaussian_Molecule_Data
!
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!
!     Interface BasisType2FunctionInfo...
      interface basisType2FunctionInfo
        module procedure MQC_Gaussian_Unformatted_Matrix_Basis_Type_To_Function_Info
      end interface
!
!
!----------------------------------------------------------------
!                                                               |
!     MODULE VARIABLES/CONSTANTS                                |
!                                                               |
!----------------------------------------------------------------
!
      Logical,Private::MQC_Gaussian_DEBUG=.False.
!
!
!     Subroutines/Functions...
!
      CONTAINS
!
!
!PROCEDURE MQC_Gaussian_SetDEBUG
      subroutine MQC_Gaussian_SetDEBUG(setDebugValue)
!
!     This subroutine is called to set a Module-Wide DEBUG flag to .TRUE. or
!     .FALSE.
!
!     H. P. Hratchian, 2020.
!
!
      implicit none
      logical,intent(IN)::setDebugValue
!
      MQC_Gaussian_DEBUG = setDebugValue
!
      return
      end subroutine MQC_Gaussian_SetDEBUG


!
!PROCEDURE MQC_Gaussian_ICGU
      subroutine MQC_Gaussian_ICGU(ICGU,wf_type,wf_complex)
!
!     This subroutine interprets the Gaussian ICGU flag (picked up from fchk and
!     matrix files). There are two wavefunction characteristics that can be
!     determined from ICGU: the spin type (Restricted, Unrestricted, General)
!     and the complex/real type (complex or real).
!
!     <ICGU> is an INPUT integer argument. <wf_type> is an OUTPUT character
!     argument that is filled with 'R', 'U', or 'G'. <wf_complex> is an OUTPUT
!     logical argument that is returned as TRUE is the wavefunction is complex.
!     Both <wf_type> and <wf_complex> are OPTIONAL arguments.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      integer,intent(IN)::ICGU
      character(len=*),intent(OUT),OPTIONAL::wf_type
      logical,intent(OUT),OPTIONAL::wf_complex
!
      if(PRESENT(wf_type)) then
        if (Mod(ICGU,1000)/100.eq.2) then
          wf_type = 'G'
        elseIf (Mod(ICGU,1000)/100.eq.1) then
          if (Mod(ICGU,10).eq.1) then
            wf_type = 'R'
          elseIf (Mod(ICGU,10).eq.2) then
            wf_type = 'U'
          else
            call MQC_Error_I('Unknown flag at ICGU 1st digit in MQC_Gaussian_ICGU ', 6, &
               'ICGU', ICGU )
          endIf
        else
          call MQC_Error_I('Unknown flag at ICGU 3rd digit in MQC_Gaussian_ICGU ', 6, &
               'ICGU', ICGU )
        endIf
      endIf
      if(PRESENT(wf_complex)) then
        if (Mod(ICGU,100)/10.eq.2) then
          wf_complex = .true.
        elseIf (Mod(ICGU,100)/10.eq.1) then
          wf_complex = .false.
        else
          Call MQC_Error_I('Unknown flag at ICGU 2nd digit in MQC_Gaussian_ICGU ', 6, &
               'ICGU', ICGU )
        endIf
      endIf
!
      return
      end subroutine MQC_Gaussian_ICGU

 
!
!PROCEDURE MQC_Gaussian_IsRestricted
      function MQC_Gaussian_IsRestricted(fileinfo)
!
!     This LOGICAL function returns TRUE if the job described by the Gaussian
!     matrix file defined by object <fileinfo> corresponds to a spin-restricted
!     wavefunction.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(in)::fileinfo
      logical::MQC_Gaussian_IsRestricted
      character(len=1)::wf_type
!
      call MQC_Gaussian_ICGU(fileinfo%icgu,wf_type=wf_type)
      MQC_Gaussian_IsRestricted = (wf_type=='R'.or.wf_type=='r')
!
      return
      end function MQC_Gaussian_IsRestricted
!
!
!PROCEDURE MQC_Gaussian_IsUnrestricted
      function MQC_Gaussian_IsUnrestricted(fileinfo)
!
!     This LOGICAL function returns TRUE if the job described by the Gaussian
!     matrix file defined by object <fileinfo> corresponds to a
!     spin-unrestricted wavefunction.
!
!     H. P. Hratchian, 2017.
!
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(in)::fileinfo
      logical::MQC_Gaussian_IsUnrestricted
      character(len=1)::wf_type
!
      call MQC_Gaussian_ICGU(fileinfo%icgu,wf_type=wf_type)
      MQC_Gaussian_IsUnrestricted = (wf_type=='U'.or.wf_type=='u')
!
      return
      end function MQC_Gaussian_IsUnrestricted
!
!
!PROCEDURE MQC_Gaussian_IsGeneral     
      function MQC_Gaussian_IsGeneral(fileinfo)
!
!     This LOGICAL function returns TRUE if the job described by the Gaussian
!     matrix file defined by object <fileinfo> corresponds to a general-spin 
!     wavefunction.
!
!     L. M. Thompson, 2017.
!
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(in)::fileinfo
      logical::MQC_Gaussian_IsGeneral
      character(len=1)::wf_type
!
      call MQC_Gaussian_ICGU(fileinfo%icgu,wf_type=wf_type)
      MQC_Gaussian_IsGeneral = (wf_type=='G'.or.wf_type=='g')
!
      return
      end function MQC_Gaussian_IsGeneral
!
!
!PROCEDURE MQC_Gaussian_IsComplex     
      function MQC_Gaussian_IsComplex(fileinfo)
!
!     This LOGICAL function returns TRUE if the job described by the Gaussian
!     matrix file defined by object <fileinfo> corresponds to a complex wavefunction.
!
!     L. M. Thompson, 2017.
!
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(in)::fileinfo
      logical::MQC_Gaussian_IsComplex,wf_complex
!
      call MQC_Gaussian_ICGU(fileinfo%icgu,wf_complex=wf_complex)
      MQC_Gaussian_IsComplex = wf_complex
!
      return
      end function MQC_Gaussian_IsComplex
!
!
!PROCEDURE MQC_Gaussian_FChk_Open
      Subroutine MQC_Gaussian_FChk_Open(FileInfo,FileName,UnitNumber, &
        OK,fileAction)
!
!     This Routine is used to connect a Gaussian formatted checkpoint file. Note
!     that optional dummy argument <fileAction> is ALWAYS IGNORED by this
!     routine.
!
!
!     Variable Declarations.
!
      Implicit None
      Class(MQC_Gaussian_FChk_File),Intent(InOut)::FileInfo
      Character(Len=*),Intent(In)::FileName
      Integer,Intent(In)::UnitNumber
      Logical,optional,Intent(Out)::OK
      character(len=*),intent(in),optional::fileAction
      Integer::myUnitNumber
!
      Character(Len=256)::Temp_Char
      Integer::IError
      Logical::EOF,Temp_Logical
!
!
!     Begin by opening the FChk file as an MQC text file and then set
!     Blank_Is_EOF to .False..
!
      if(UnitNumber.gt.0) then
        myUnitNumber = UnitNumber
      else
        myUnitNumber = MQC_getUnitNumber()
      endIf
      Call MQC_Open_Text_File(FileInfo,FileName,myUnitNumber,OK)
      If(.not.OK) Return
      call MQC_Files_Text_Parse_Controls(FileInfo,  &
        Do_Blank_Is_EOF=.false.)
!
!     Load the title, job type, method, and basis set from the first two
!     lines of the file.
!
      Call FileInfo%GetBuffer(FileInfo%Title)
      If (FileInfo%LoadBuffer()) Call MQC_Error_L('Failiure loading FChk buffer', 6, &
           'FileInfo%LoadBuffer()', FileInfo%LoadBuffer() )
      Call FileInfo%GetNextString(FileInfo%JobType,EOF,OK)
      Call FileInfo%GetNextString(FileInfo%Method,EOF,OK)
      Call FileInfo%GetNextString(FileInfo%BasisSet,EOF,OK)
!
      Return
      End Subroutine MQC_Gaussian_FChk_Open


!
!PROCEDURE Find_FChk_Entry
      Subroutine Find_FChk_Entry(FChkFile,EntryTag,FoundEntry,  &
        TypeFlag,NElements,Scalar_Integer,Scalar_Real,Scalar_Character, &
        Scalar_Logical,Vector_Integer,Vector_Real)
!
!     This subroutine searches a Gaussian formatted checkpoint file
!     (FChkFile) for an entry tag (EntryTag). It returns a logical
!     indicating if the entry was found (FoundEntry), what sort of entry
!     was found (TypeFlag), and the number of data elements of the entry
!     (NElements).
!
!     More specifically, FoundEntry is returned TRUE if the desired entry
!     is found. TypeFlag is a single character variable that will be
!     returned with 'I', 'R', or 'A' for integer, real, or character
!     (alpha) data type found. NElements will be returned with 0 if the
!     data entry is a scalar. If the data is an array, NElements will be
!     returned with the length of the array.
!
!     If a scalar value is found, the correct output variable of Scalar_*
!     is filled as long as it has been passed to this routine (as the
!     Scalar_* arguments are all OPTIONAL).
!
!     If an array is found, the correct output variable of Vector_* is
!     filled as long as it has been passed to this routine (as the Vector_*
!     arguments are all OPTIONAL).
!
!
!     H. P. Hratchian, 2016.
!
!
!     Variable Declarations
!
      Implicit None
      Character(Len=*),Intent(IN)::EntryTag
      Type(MQC_Gaussian_FChk_File),Intent(InOut)::FChkFile
      Logical,Intent(OUT)::FoundEntry
      Character(*),Intent(OUT)::TypeFlag
      Integer(kind=int64),Intent(OUT)::NElements
      Integer(kind=int64),Intent(OUT),OPTIONAL::Scalar_Integer
      Real(kind=real64),Intent(OUT),OPTIONAL::Scalar_Real
      Character(*),Intent(OUT),OPTIONAL::Scalar_Character
      Logical,Intent(OUT),OPTIONAL::Scalar_Logical
      Integer(kind=int64),Dimension(:),Intent(OUT),OPTIONAL::Vector_Integer
      Real(kind=real64),Dimension(:),Allocatable,OPTIONAL::Vector_Real
!
      Character(Len=40)::Lab
      Character(Len=2)::PlaceholderNEq
      Character(Len=22)::ValueHolder
      Character(Len=1024)::Temp_Char
      Logical::Fail,OK,EOF,IsArray
!
!
!     Format Statements
!
 1060 Format( A )
 2000 Format(A40,3X,A1,3X,A2,A22)
 8000 Format(1x,'Label: ',A40,' | Type=',A1,' | NEq=',A2,  &
        ' | Value=',A12)
!
!
!     Initialize FoundEntry, TypeFlag, and NElements. Then, ensure the fchk
!     is open.
!
      FoundEntry = .False.
      TypeFlag = ' '
      NElements = -1
      If(.not.FChkFile%IsOpen()) Return
!
!     Search the file for the entry tag. If it isn't found in the first
!     try, rewind the file and try again.
!
      Call MQC_Files_FF_Search(Fail,FChkFile,EntryTag,.True.)
      If(Fail) then
        Call FChkFile%Rewind(OK)
        If(.not.OK) Return
        Call MQC_Files_FF_Search(Fail,FChkFile,EntryTag,.True.)
      endIf
      If(Fail) Return
      FoundEntry = .True.
      Call FChkFile%GetBuffer(Temp_Char,.False.)
      Read(Temp_Char,2000) Lab,TypeFlag,PlaceholderNEq,ValueHolder
      Write(*,8000) Lab,TypeFlag,PlaceholderNEq,ValueHolder
      If(PlaceholderNEq.eq.'N=') then
        Read(ValueHolder,'(I22)') NElements
        fail = FChkFile%LoadBuffer()
        If(fail) Return
        Select Case (TypeFlag)
        Case('I')
          Write(*,1060)' Found an INTEGER ARRAY'
          If(Present(Vector_Integer)) then
            If(Size(Vector_Integer).ge.NElements) &
              Call MQC_Files_Text_File_Read_Int_Vec(FChkFile,  &
                Vector_Integer,EOF,OK)
          endIf
        Case('R')
          Write(*,1060)' Found a REAL ARRAY'
          If(Present(Vector_Real)) then
            if(.not.Allocated(Vector_Real)) then
              Allocate(Vector_Real(NElements))
            elseIf(Size(Vector_Real).ne.NElements) then
              DeAllocate(Vector_Real)
              Allocate(Vector_Real(NElements))
            endIf
            Call MQC_Files_Text_File_Read_Real_Vec(FChkFile,Vector_Real,EOF,OK)
          endIf
        Case('C')
          Write(*,1060)' Found a CHARACTER ARRAY'
        Case('L')
          Write(*,1060)' Found a LOGICAL ARRAY'
        Case Default
          Return
        End Select
      else
        Select Case (TypeFlag)
        Case('I')
          If(Present(Scalar_Integer))  &
            Read(ValueHolder,'(I22)') Scalar_Integer
        Case('R')
          If(Present(Scalar_Real))  &
            Read(ValueHolder,'(E22.15)') Scalar_Real
        Case('C')
          If(Present(Scalar_Character))  &
            Read(ValueHolder,'(A12)') Scalar_Character
        Case('L')
          If(Present(Scalar_Logical))  &
            Read(ValueHolder,'(L1)') Scalar_Logical
        Case Default
          Return
        End Select
      endIf
!
      Return
      End Subroutine Find_FChk_Entry


!
!PROCEDURE MQC_Gaussian_Fill_Molecule_Data_FChk
      Subroutine MQC_Gaussian_Fill_Molecule_Data_FChk(FChkFile,  &
        MoleculeData)
!
!     This subroutine is used to fill a MQC_Molecule_Data type variable
!     from a provided (and already open) Gaussian formatted checkpoint
!     file.
!
!     H. P. Hratchian, 2016.
!
!
      Implicit None
      Type(MQC_Gaussian_FChk_File),Intent(InOut)::FChkFile
      Class(MQC_Gaussian_Molecule_Data),Intent(InOut)::MoleculeData
!
      Integer::NElements,NAtoms,Charge,Multiplicity
      Integer,Dimension(:),Allocatable::AtomicNumbers
      Real,Dimension(:),Allocatable::AtomicMasses,NuclearCharges,  &
        Cartesians_1D
      Real,Dimension(:,:),Allocatable::Cartesians
      Character(Len=1024)::Temp_Char
      Logical::OK
!
!     Get the integer flags from the fchk file and then allocate the local
!     arrays appropriately.
!
      Call Find_FChk_Entry(FChkFile,                     &
        'Number of atoms                            I',  &
        OK,Temp_Char,NElements,Scalar_Integer=NAtoms)
      If(.not.OK)  &
        Call MQC_Error_L('MQC_Gaussian: FChk error loading NAtoms.', 6, &
        'OK', OK )
      Call Find_FChk_Entry(FChkFile,                     &
        'Charge                                     I',  &
        OK,Temp_Char,NElements,Scalar_Integer=Charge)
      If(.not.OK)  &
        Call MQC_Error_L('MQC_Gaussian: FChk error loading Charge.', 6, &
        'OK', OK )
      Call Find_FChk_Entry(FChkFile,                     &
        'Multiplicity                               I',  &
        OK,Temp_Char,NElements,Scalar_Integer=Multiplicity)
      If(.not.OK)  &
        Call MQC_Error_L('MQC_Gaussian: FChk error loading Multiplicity.', 6, &
        'OK', OK )
      Allocate(AtomicNumbers(NAtoms),AtomicMasses(NAtoms),  &
        NuclearCharges(NAtoms),Cartesians_1D(3*NAtoms),  &
        Cartesians(3,NAtoms))
!
!     Now, fill the arrays from the fchk file.
!
      Call Find_FChk_Entry(FChkFile,                     &
        'Atomic numbers                             I',  &
        OK,Temp_Char,NElements,Vector_Integer=AtomicNumbers)
      If(.not.OK)  &
        Call MQC_Error_L('MQC_Gaussian: FChk error - AtomicNumbers.', 6, &
        'OK', OK )
      Call Find_FChk_Entry(FChkFile,                     &
        'Nuclear charges                            R',  &
        OK,Temp_Char,NElements,Vector_Real=NuclearCharges)
      If(.not.OK)  &
        Call MQC_Error_L('MQC_Gaussian: FChk error - NuclearCharges.', 6, &
        'OK', OK )
      Call Find_FChk_Entry(FChkFile,                     &
        'Current cartesian coordinates              R',  &
        OK,Temp_Char,NElements,Vector_Real=Cartesians_1D)
      If(.not.OK)  &
        Call MQC_Error_L('MQC_Gaussian: FChk error - Cartesians.', 6, &
        'OK', OK )
      Cartesians = Reshape(Cartesians_1D,(/3,NAtoms/))
      Call Find_FChk_Entry(FChkFile,                     &
        'Real atomic weights                        R',  &
        OK,Temp_Char,NElements,Vector_Real=AtomicMasses)
      If(.not.OK)  &
        Call MQC_Error_L('MQC_Gaussian: FChk error - AtomicMasses.', 6, &
        'OK', OK )
!
!     Fill the molecular data object. Then de-allocate arrays.
!
      Call MQC_Gaussian_Fill_Molecule_Data(MoleculeData,NAtoms,  &
        AtomicNumbers,AtomicMasses,NuclearCharges,Cartesians,Charge,  &
        Multiplicity)
      DeAllocate(AtomicNumbers,AtomicMasses,NuclearCharges,  &
        Cartesians_1D,Cartesians)
!
      End Subroutine MQC_Gaussian_Fill_Molecule_Data_FChk


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Fill_Molecule_Data
      Subroutine MQC_Gaussian_Fill_Molecule_Data(MoleculeData,NAtoms,  &
        AtomicNumbers,AtomicMasses,NuclearCharges,Cartesians,  &
        TotalCharge,Multiplicity)
!
!     This subroutine is used to fill a MQC_Molecule_Data type variable
!     given INPUT dummy arguments for each of its constituents entries.
!     All the molecule objects need to be updated.
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(Out)::MoleculeData
      Integer,Intent(In)::NAtoms
      Integer,Dimension(NAtoms),Intent(In)::AtomicNumbers
      Real,Dimension(NAtoms),Intent(In)::AtomicMasses,NuclearCharges
      Real,Dimension(3,NAtoms),Intent(In)::Cartesians
      Integer,Intent(In)::TotalCharge,Multiplicity
!
!     Add data to the standard MQC_Molecule_Data elements.
!
      Call MQC_Molecule_Data_Fill(MoleculeData,NAtoms,  &
        AtomicNumbers,AtomicMasses,NuclearCharges,Cartesians)
      Select Type(MoleculeData)
      Type is(MQC_Gaussian_Molecule_Data)
        MoleculeData%Charge = TotalCharge
        MoleculeData%Multiplicity = Multiplicity
      endSelect
!
      End Subroutine MQC_Gaussian_Fill_Molecule_Data


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Output_Molecule_Data
      Subroutine MQC_Gaussian_Output_Molecule_Data(MoleculeData,NAtoms,  &
        AtomicNumbers,AtomicMasses,NuclearCharges,Cartesians,  &
        TotalCharge,Multiplicity)
!
!     This subroutine is used to output a MQC_Molecule_Data type variable.
!
      Implicit None
      Class(MQC_Molecule_Data),Intent(In)::MoleculeData
      Integer,Allocatable,Intent(InOut)::NAtoms
      Integer,Dimension(:),Allocatable,Intent(InOut)::AtomicNumbers
      Real,Dimension(:),Allocatable,Intent(InOut)::AtomicMasses,NuclearCharges,Cartesians
      Real,Dimension(:,:),Allocatable::TmpCart
      Integer,Optional,Intent(InOut)::TotalCharge,Multiplicity
!
      Call MQC_Molecule_Data_Output(MoleculeData,NAtoms,  &
        AtomicNumbers,AtomicMasses,NuclearCharges,TmpCart)
      Cartesians = reshape(TmpCart,[3*NAtoms])
      Select Type(MoleculeData)
      Type is(MQC_Gaussian_Molecule_Data)
        TotalCharge = MoleculeData%Charge
        Multiplicity = MoleculeData%Multiplicity
      endSelect
!
      End Subroutine MQC_Gaussian_Output_Molecule_Data


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Open
      subroutine MQC_Gaussian_Unformatted_Matrix_Open(fileinfo,filename,unitNumber,ok)
!
!     This Routine is used to set-up a Gaussian unformatted matrix file object
!     file. The dummy argument <FileName> is an input argument giving the the
!     name of the file. For now, the dummy argument <unitNumber> must be sent,
!     but it is ultimately ignored. The output dummy argument <ok> is returned
!     TRUE if everything proceeds without error; otherwise, <ok> is returned
!     FALSE.
!
!     NOTE: This routine does NOT actually open the file. Instead, this routine
!     serves as the required OpenFile deferred procedure binding for members of
!     the MQC_FileInfo class. Because the unformatted matrix file is a Fortran
!     unformatted file, the file will only be read OR write. So, there are two
!     different routines for working with the file -- one for READ mode and one
!     for WRITE mode. Also, note that the matrix file has header scalars and
!     arrays that must initially be read/written. Then arrays/matrices can be
!     read/written.
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::filename
      integer,intent(in)::unitNumber
      logical,intent(out)::ok
!
      integer::iout=6,unit_number
!
!
!     Begin by opening the file.
!
      if(fileinfo%IsOpen()) call MQC_Gaussian_Unformatted_Matrix_Close(fileinfo)
      ok = .true.
      fileinfo%filename       = TRIM(filename)
      fileinfo%unitNumber     = 0
      if(MQC_Gaussian_DEBUG) write(*,*)' In OpenFile, unitNumber = ',  &
        fileinfo%unitNumber
      fileinfo%declared       = .false.
      fileinfo%CurrentlyOpen  = .true.
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Open


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Close
      subroutine MQC_Gaussian_Unformatted_Matrix_Close(fileinfo)
!
!     This Routine is used to close a Gaussian unformatted matrix file.
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
!
!
!     Close the matrix file using the gauopen routines.
!
      if(fileinfo%isOpen()) call Close_MatF(fileinfo%UnitNumber)
      fileinfo%filename       = ' '
      fileinfo%CurrentlyOpen  = .false.
      fileinfo%UnitNumber     = 0
      fileinfo%declared       = .false.
      fileinfo%header_read    = .false.
      fileinfo%header_written = .false.
      fileinfo%readWriteMode  = ' '
      fileinfo%LabFil         = ' '
      fileinfo%GVers          = ' '
      fileinfo%Title          = ' '
      if(allocated(fileinfo%atomicNumbers)) deallocate(fileinfo%atomicNumbers)
      if(allocated(fileinfo%atomTypes)) deallocate(fileinfo%atomTypes)
      if(allocated(fileinfo%basisFunction2Atom)) deallocate(fileinfo%basisFunction2Atom)
      if(allocated(fileinfo%IBasisFunctionType)) deallocate(fileinfo%IBasisFunctionType)
      if(allocated(fileinfo%atomicCharges)) deallocate(fileinfo%atomicCharges)
      if(allocated(fileinfo%atomicWeights)) deallocate(fileinfo%atomicWeights)
      if(allocated(fileinfo%cartesians)) deallocate(fileinfo%cartesians)
      if(allocated(fileinfo%natoms)) deallocate(fileinfo%natoms)
      if(allocated(fileinfo%nbasis)) deallocate(fileinfo%nbasis)
      if(allocated(fileinfo%nbasisUse)) deallocate(fileinfo%nbasisUse)
      if(allocated(fileinfo%icharge)) deallocate(fileinfo%icharge)
      if(allocated(fileinfo%multiplicity)) deallocate(fileinfo%multiplicity)
      if(allocated(fileinfo%nelectrons)) deallocate(fileinfo%nelectrons)
      if(allocated(fileinfo%icgu)) deallocate(fileinfo%icgu)
      if(allocated(fileinfo%NFC)) deallocate(fileinfo%NFC)
      if(allocated(fileinfo%NFV)) deallocate(fileinfo%NFV)
      if(allocated(fileinfo%ITran)) deallocate(fileinfo%ITran)
      if(allocated(fileinfo%IDum9)) deallocate(fileinfo%IDum9)
      if(allocated(fileinfo%NShlAO)) deallocate(fileinfo%NShlAO)
      if(allocated(fileinfo%NPrmAO)) deallocate(fileinfo%NPrmAO)
      if(allocated(fileinfo%NShlDB)) deallocate(fileinfo%NShlDB)
      if(allocated(fileinfo%NPrmDB)) deallocate(fileinfo%NPrmDB)
      if(allocated(fileinfo%NBTot)) deallocate(fileinfo%NBTot)
      if(allocated(fileinfo%gaussianScalars)) deallocate(fileinfo%gaussianScalars)
      call MQC_setOff_unitNumbersOpened(fileinfo%UnitNumber)
!
      return
      end Subroutine MQC_Gaussian_Unformatted_Matrix_Close


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Read_Header
      subroutine MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
        filename)
!
!     This Routine is used to connect a Gaussian unformatted matrix file
!     and read the header records. If sent, this routine also loads an
!     MQC_Molecule_Data object. The dummy argument <FileName> is an input
!     argument giving the the name of the file.
!
!     Dummy argument <filename> is optional and is only used if fileinfo
!     hasn't already been defined using Routine
!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!     filename sent is different from the filename associated with object
!     fileinfo.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open may be called
!     before calling this routine. However, it is also OK to call this
!     routine first. In that case, this routine will first call
!     Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in),OPTIONAL::filename
!
      integer::iout=6
!
!     Temporary local variables used when calling the gauopen routines.
      integer::IVers,NLab,Len12L,Len4L,IOpCl
!      integer::NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot
      character(len=64)::cBuffer
!
!     Local temp variables.
      real,dimension(:),allocatable::tempArray
      character(len=256)::my_filename
      logical::DEBUG=.false.,ok
!
!
!     Format statements.
!
 1000 format(1x,'Reading data from file: ',A,/,1x,'Unit number: ',I3,/)
 1010 format(3x,' Label ',A,' IVers=',I2,' NLab=',I2,' Version=',A,  &
        /,3x,' Title ',A,  &
        /,3x,' NAtoms=',I6,' NBasis=',I6,' NBsUse=',I6,' ICharg=',I6,  &
        ' Multip=',I6,' NElec=',I6,' Len12L=',I1,' Len4L=',I1,' IOpCl=',I6,  &
        ' ICGU=',I3)
 1020 format(3x,'BasisFunction2Atom Map:')
 1022 format(3x,10(2x,i5))
!
!
!     Begin by seeing if a new file or filename has been sent by the calling
!     program unit. If so, then get the file declared before reading the
!     header information.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
             'ok', ok )
        else
          call MQC_Error_L('Error reading Gaussian matrix file header (1): Must include a filename.', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
             'ok', ok )
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'R' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call fileinfo%OPENFILE(my_filename,0,ok)
        if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
             'ok', ok )
      endIf
      if(fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call fileinfo%OPENFILE(my_filename,0,ok)
        if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
             'ok', ok )
      endIf
!
!     Set the readWriteMode flag in fileinfo to 'R' and then read the
!     header scalar flags.
!
      if(.not.fileinfo%header_read) then
        fileinfo%readWriteMode = 'R'
        allocate(fileinfo%natoms,fileinfo%nbasis,fileinfo%nbasisUse,fileinfo%icharge, &
          fileinfo%multiplicity,fileinfo%nelectrons,fileinfo%icgu,fileinfo%NFC, &
          fileinfo%NFV,fileinfo%ITran,fileinfo%IDum9,fileinfo%NShlAO,fileinfo%NPrmAO, &
          fileinfo%NShlDB,fileinfo%NPrmDB,fileinfo%NBTot)
        call Open_Read(TRIM(fileinfo%filename),fileinfo%UnitNumber,  &
          fileinfo%labfil,ivers,nlab,fileinfo%gvers,fileinfo%title,  &
          fileinfo%natoms,fileinfo%nbasis,fileinfo%nbasisUse,  &
          fileinfo%icharge,fileinfo%multiplicity,fileinfo%nelectrons,len12l,  &
          len4l,iopcl,fileinfo%icgu)
        call MQC_setOn_unitNumbersOpened(fileinfo%UnitNumber)
        allocate(fileinfo%atomicNumbers(fileinfo%natoms),  &
          fileinfo%atomTypes(fileinfo%natoms),  &
          fileinfo%atomicCharges(fileinfo%natoms),  &
          fileinfo%atomicWeights(fileinfo%natoms))
        allocate(fileinfo%cartesians(fileinfo%natoms*3))
        allocate(fileinfo%basisFunction2Atom(fileinfo%NBasis))
        allocate(fileinfo%IBasisFunctionType(fileinfo%NBasis))
        call Rd_Head(fileinfo%unitNumber,NLab,fileinfo%natoms,fileinfo%nbasis,  &
          fileinfo%atomicNumbers,fileinfo%atomTypes,fileinfo%atomicCharges,  &
          fileinfo%cartesians,fileinfo%basisFunction2Atom,fileinfo%IBasisFunctionType,  &
          fileinfo%atomicWeights,fileinfo%NFC,fileinfo%NFV,fileinfo%ITran, &
          fileinfo%IDum9,fileinfo%NShlAO,fileinfo%NPrmAO,fileinfo%NShlDB,  &
          fileinfo%NPrmDB,fileinfo%NBTot)
        fileinfo%CurrentlyOpen = .true.
        fileinfo%header_read   = .true.
      endIf
      if(DEBUG.or.MQC_Gaussian_DEBUG) then
        write(IOut,1010) TRIM(fileinfo%LabFil),IVers,NLab,  &
          TRIM(fileinfo%GVers),TRIM(fileinfo%Title),fileinfo%natoms,  &
          fileinfo%NBasis,fileinfo%nbasisUse,fileinfo%ICharge,  &
          fileinfo%Multiplicity,fileinfo%nelectrons,Len12L,Len4L,  &
          IOpCl,fileinfo%ICGU
        write(iOut,1020)
        write(iOut,1022) fileinfo%basisFunction2Atom
      endIf
!
!     Load the Gaussian scalars arrays from the Matrix File.
!

!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Read_Header


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Write_Header
      subroutine MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
        filename)
!
!     This Routine is used to connect a Gaussian unformatted matrix file
!     and write the header records stored in fileinfo. The dummy argument 
!     <FileName> is an input argument giving the the name of the file.
!
!     Dummy argument <filename> is optional. If filename is not specified,
!     the filename in fileinfo is used.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in),OPTIONAL::filename
!
      integer::iout=6
!
!     Temporary local variables used when calling the gauopen routines.
      integer::IOpCl=-1
      integer::NAt3,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot
!
!     Local temp variables.
      character(len=256)::my_filename
      logical::DEBUG=.true.,ok,openbool
!
!
!     Format statements.
!
 1010 format(3x,' Label ',A,' Version=',A,  &
        /,3x,' Title ',A,  &
        /,3x,' NAtoms=',I6,' NBasis=',I6,' NBsUse=',I6,' ICharg=',I6,  &
        ' Multip=',I6,' NElec=',I6,' IOpCl=',I6,  &
        ' ICGU=',I3)
!
!
!     Begin by seeing if a new file or filename has been sent by the calling
!     program unit. If so, then get the file declared before writing the
!     header information.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
               'ok', ok )
        else
          call MQC_Error_L('Error reading Gaussian matrix file header (2): Must include a filename.', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          fileinfo%filename=trim(filename)
        endIf
      endIf
      inquire(file=fileinfo%filename,number=fileinfo%unitNumber,opened=openbool)
      if(openbool) close(fileInfo%unitNumber)
!
!     Check if all the required information is in fileinfo, and if no then
!     initialize it.
!
      if(.not.allocated(fileinfo%natoms)) then
        allocate(fileinfo%natoms)
        fileinfo%natoms = 0
      endIf
      if(.not.allocated(fileinfo%nbasis)) then 
        allocate(fileinfo%nbasis)
        fileinfo%nbasis = 0
      endIf
      if(.not.allocated(fileinfo%nbasisUse)) then 
        allocate(fileinfo%nbasisUse)
        fileinfo%nbasisUse = 0
      endIf
      if(.not.allocated(fileinfo%icharge)) then 
        allocate(fileinfo%icharge)
        fileinfo%icharge = 0
      endIf
      if(.not.allocated(fileinfo%multiplicity)) then 
        allocate(fileinfo%multiplicity)
        fileinfo%multiplicity = 0
      endIf
      if(.not.allocated(fileinfo%nelectrons)) then 
        allocate(fileinfo%nelectrons)
        fileinfo%nelectrons = 0
      endIf
      if(.not.allocated(fileinfo%icgu)) then 
        allocate(fileinfo%icgu)
        fileinfo%icgu = 0
      endIf
      if(.not.allocated(fileinfo%NFC)) then 
        allocate(fileinfo%NFC)
        fileinfo%NFC = 0
      endIf
      if(.not.allocated(fileinfo%NFV)) then 
        allocate(fileinfo%NFV)
        fileinfo%NFV = 0
      endIf
      if(.not.allocated(fileinfo%ITran)) then 
        allocate(fileinfo%ITran)
        fileinfo%ITran = 0
      endIf
      if(.not.allocated(fileinfo%IDum9)) then 
        allocate(fileinfo%IDum9)
        fileinfo%IDum9 = 0
      endIf
      if(.not.allocated(fileinfo%NShlAO)) then 
        allocate(fileinfo%NShlAO)
        fileinfo%NShlAO = 0
      endIf
      if(.not.allocated(fileinfo%NPrmAO)) then 
        allocate(fileinfo%NPrmAO)
        fileinfo%NPrmAO = 0
      endIf
      if(.not.allocated(fileinfo%NShlDB)) then 
        allocate(fileinfo%NShlDB)
        fileinfo%NShlDB = 0
      endIf
      if(.not.allocated(fileinfo%NPrmDB)) then 
        allocate(fileinfo%NPrmDB)
        fileinfo%NPrmDB = 0
      endIf
      if(.not.allocated(fileinfo%NBTot)) then 
        allocate(fileinfo%NBTot)
        fileinfo%NBTot = 0
      endIf
      nAt3 = fileinfo%natoms*3
      if(.not.allocated(fileinfo%atomicNumbers)) then
        allocate(fileinfo%atomicNumbers(fileinfo%natoms))
        fileinfo%atomicNumbers = 0
      endIf
      if(.not.allocated(fileinfo%atomTypes)) then
        allocate(fileinfo%atomTypes(fileinfo%natoms))
        fileinfo%atomTypes = 0
      endIf
      if(.not.allocated(fileinfo%atomicCharges)) then
        allocate(fileinfo%atomicCharges(fileinfo%natoms))
        fileinfo%atomicCharges = 0.0
      endIf
      if(.not.allocated(fileinfo%atomicWeights)) then
        allocate(fileinfo%atomicWeights(fileinfo%natoms))
        fileinfo%atomicWeights = 0.0
      endIf
      if(.not.allocated(fileinfo%cartesians)) then
        allocate(fileinfo%cartesians(nAt3))
        fileinfo%cartesians = 0.0
      endIf
      if(.not.allocated(fileinfo%basisFunction2Atom)) then
        allocate(fileinfo%basisFunction2Atom(fileinfo%NBasis))
        fileinfo%basisFunction2Atom = 0
      endIf
     
      if(fileinfo%icgu.eq.111) then
        iopcl = 0 
      elseIf(fileinfo%icgu.eq.112) then
        iopcl = 1
      elseIf(fileinfo%icgu.eq.121) then
        iopcl = 2
      elseIf(fileinfo%icgu.eq.122) then
        iopcl = 3
      elseIf(fileinfo%icgu.eq.221) then
        iopcl = 6
      endIf
!
!     Set the readWriteMode flag in fileinfo to 'W' and then write the
!     header scalar flags.
!
      fileinfo%readWriteMode = 'W'
      call Open_Write(TRIM(fileinfo%filename),fileinfo%UnitNumber,  &
        fileinfo%labfil,fileinfo%gvers,fileinfo%title,  &
        fileinfo%natoms,fileinfo%nbasis,fileinfo%nbasisUse,  &
        fileinfo%icharge,fileinfo%multiplicity,fileinfo%nelectrons,iopcl, &
        fileinfo%icgu)
!
      call Wr_Head(fileinfo%unitNumber,fileinfo%natoms,nAt3,fileinfo%nbasis,  &
        fileinfo%atomicNumbers,fileinfo%atomTypes,fileinfo%atomicCharges,  &
        fileinfo%cartesians,fileinfo%basisFunction2Atom,fileinfo%IBasisFunctionType,  &
        fileinfo%atomicWeights,fileinfo%NFC,fileinfo%NFV,fileinfo%ITran,fileinfo%IDum9, &
        fileinfo%NShlAO,fileinfo%NPrmAO,fileinfo%NShlDB,fileinfo%NPrmDB,fileinfo%NBTot)
      fileinfo%CurrentlyOpen = .true.
      fileinfo%header_written   = .true.
!      
      if(DEBUG) write(IOut,1010) TRIM(fileinfo%LabFil),  &
        TRIM(fileinfo%GVers),TRIM(fileinfo%Title),fileinfo%natoms,  &
        fileinfo%NBasis,fileinfo%nbasisUse,fileinfo%ICharge,  &
        fileinfo%Multiplicity,fileinfo%nelectrons,  &
        IOpCl,fileinfo%ICGU
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Write_Header


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Read_Array
      subroutine MQC_Gaussian_Unformatted_Matrix_Read_Array(fileinfo,  &
       label,matrixOut,vectorOut,r4TensorOut,filename,mqcVarOut,foundOut, &
       arraynum)
!
!     This Routine is used to look-up a matrix in a unformatted matrix file load
!     that array into either (OPTIONAL) output dummy MQC_Matrix argument
!     <matrixOut>, (OPTIONAL) output dummy MQC_Vector argument <vectorOut>,
!     (OPTIONAL) output dummy MQC_R4Tensor argument <r4TensorOut>, or (OPTIONAL)
!     output dummy MQC_Variable argument <mqcVarOut>. The character label for
!     the array of interest is sent to this routine in dummy argument <label>.
!
!     Dummy argument <filename> is optional and is only used if fileinfo hasn't
!     already been defined using Routine MQC_Gaussian_Unformatted_Matrix_Open or
!     if it is determined that the filename sent is different from the filename
!     associated with object fileinfo.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!     called before calling this routine. The expectation is that
!     MQC_Gaussian_Unformatted_Matrix_Read_Header is also called before this
!     routine. However, it is also OK to call this routine first. In that case,
!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     Optional argument <arraynum> allows the user to specify a specific array
!     in the case that the matrix file has several arrays under a specific
!     label. Currently only implemented for the case '2ERIS-SYMSYMR4TENSOR' so
!     that we can access Raffenetti 2 and Raffenetti 3 integrals.
!
!     H. P. Hratchian, 2017, 2018.
!     L. M. Thompson, 2017, 2020
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      type(MQC_Matrix),intent(inout),OPTIONAL::matrixOut
      type(MQC_Vector),intent(inout),OPTIONAL::vectorOut
      type(MQC_R4Tensor),intent(inout),OPTIONAL::r4TensorOut
      character(len=*),intent(in),OPTIONAL::filename
      type(MQC_Variable),intent(inout),OPTIONAL::mqcVarOut
      logical,OPTIONAL::foundOut
      integer,intent(in),OPTIONAL::arrayNum
!
      integer::iout=6
!
!     Temporary local variables used when calling the gauopen routines.
      integer::IVers,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,NRI,LR
      character(len=64)::cBuffer,tmpLabel
      logical::EOF,ASym
!
!     Local temp variables.
      integer::i,nOutputArrays,LNZ,myArrayNum
!      integer,external::LenArr
      integer,allocatable,dimension(:)::integerTmp
      real,allocatable,dimension(:)::arrayTmp
      complex(kind=8),allocatable,dimension(:)::complexTmp
      character(len=256)::my_filename,errorMsg
      logical::DEBUG=.false.,ok,found
      integer::nBasis
      real,allocatable,dimension(:,:,:,:)::arrayTmp2
!
!
!     Format statements.
!
 1000 format(' MQC_Gaussian_Read is looking for array ',A,  &
        ' in matrix file.')
 1010 format(' Label ',A48,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
 1020 Format( " " )!
 1040 Format( A, I15 )
 1050 Format( 2A )
 1060 Format( A )
!
      if(MQC_Gaussian_DEBUG) DEBUG=.true.
!
!     Begin by seeing if a new file or filename has been sent by the calling
!     program unit. If so, then get the file declared before reading the
!     header information.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_l('Error reading Gaussian matrix file header (3): Must include a filename.', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'R' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Ensure that one and only one output MQC-type array has been sent from the
!     calling program unit.
!
      nOutputArrays = 0
      if(Present(matrixOut)) nOutputArrays = nOutputArrays+1
      if(Present(vectorOut)) nOutputArrays = nOutputArrays+1
      if(Present(r4TensorOut)) nOutputArrays = nOutputArrays+1
      if(Present(mqcVarOut)) nOutputArrays = nOutputArrays+1
      if(nOutputArrays.ne.1) call mqc_error_i(  &
           'Too many output arrays sent to Gaussian matrix file reading procedure.', 6, &
           'nOutputArrays', nOutputArrays )
!
      if(present(arrayNum)) then
        myArrayNum = arrayNum
      else
        myArrayNum = 1
      endIf
!
!     Look for the label sent by the calling program unit. If the label is
!     found, then load the appropriate output argument with the data on the
!     file.
!
      call String_Change_Case(label,'u',tmpLabel)
      if(DEBUG) write(IOut,1000) TRIM(tmpLabel)
      found = .false.
      outerLoop:do i = 1,2
        EOF = .false.
        Call Rd_Labl(fileinfo%UnitNumber,IVers,cBuffer,NI,NR,NTot,LenBuf,  &
          N1,N2,N3,N4,N5,ASym,NRI,EOF)
        LR = LenArr(N1,N2,N3,N4,N5)
        if(DEBUG) write(IOut,1010) TRIM(cBuffer),NI,NR,NRI,NTot,LenBuf,  &
          N1,N2,N3,N4,N5,ASym,LR
        do while(.not.EOF)
          call String_Change_Case(cBuffer,'u')
          if(TRIM(tmpLabel) == TRIM(cBuffer)) then
!           Dipole integral arrays are stored as a single rank 3 array, so 
!           to output them to an MQC algebra 1 object, we read them as 3
!           rank 2 arrays. This may cause unforeseen issues with outputting
!           MQC algebra 2 arrays.
            if((myArrayNum.gt.N3).and.(NR.ne.0.and.myArrayNum.gt.NR)) then
              if(present(foundOut)) then
                exit outerLoop
              else
                call mqc_error_I('Array number requested not found &
                  &under specified label',6,'NR',NR,'myArrayNum',myArrayNum)
              endIf
            endIf
!
!           This CASE block uses NI, NR, N1-N5, and NRI to determine the data
!           type (integer, real, etc.) and data structure (scalar, vector,
!           matrix, etc.).
!
            if(DEBUG) then
              write(IOut,*)' FOUND tmpLabel!!!'
              write(IOut,*)   &
                ' MQC_Gaussian_Unformatted_Matrix_Array_Type = ',  &
                MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,  &
                N2,N3,N4,N5,NRI,ASym)
              write(IOut,*)  &
                ' MQC_Gaussian_Unformatted_Matrix_Array_Type = ',  &
                MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI,ASym)
            endIf
            select case(MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI,ASym))
            case('INTEGER-VECTOR')
              allocate(integerTmp(LR))
              call Rd_IBuf(fileinfo%unitNumber,NTot,LenBuf,integerTmp)
              if(Present(vectorOut)) then
                vectorOut = integerTmp
              elseIf(Present(matrixOut)) then
                call MQC_Matrix_DiagMatrix_Put(matrixOut,integerTmp)
              elseIf(Present(mqcVarOut)) then
                mqcVarOut = integerTmp
              else
                call mqc_error_l('Reading vector from Gaussian matrix file, but NO VECTOR SENT to procedure.',  &
                  6,'Present(vectorOut)',Present(vectorOut),'Present(matrixOut)',Present(matrixOut) )
              endIf
              deallocate(integerTmp)
            case('INTEGER-MATRIX')
              if(.not.Present(matrixOut)) call mqc_error_l('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
                & procedure.', 6, &
                'Present(matrixOut)', Present(matrixOut) )
             allocate(integerTmp(LR))
              call Rd_IBuf(fileinfo%unitNumber,NTot,LenBuf,integerTmp)
              matrixOut = Reshape(integerTmp,[N1,N2])
              deallocate(integerTmp)
            case('INTEGER-SYMMATRIX')
              if(.not.Present(matrixOut)) call mqc_error_l('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
                & procedure.', 6, &
                'Present(matrixOut)', Present(matrixOut) )
              allocate(integerTmp(LR))
              call Rd_IBuf(fileinfo%unitNumber,NTot,LenBuf,integerTmp)
              call MQC_Matrix_SymmMatrix_Put(matrixOut,integerTmp)
              deallocate(integerTmp)
            case('INTEGER-ASYMMATRIX')
              if(.not.Present(matrixOut)) call mqc_error_l('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to &
                & procedure.', 6, &
                'Present(matrixOut)', Present(matrixOut) )
              allocate(integerTmp(LR))
              call Rd_IBuf(fileinfo%unitNumber,NTot,LenBuf,integerTmp)
              call MQC_Matrix_SymmMatrix_Put(matrixOut,integerTmp,'antisymmetric')
!             Matrix files have either symmetric/hermitian storage or antisymmetric/
!             anthermitian storage. MQC currently has a symmetric only storage for both real 
!             and complex parts so make nonsymmetric matrices square.
!              call mqc_matrix_symm2full(matrixOut,'antisymmetric')
!             Triangular matrices are stored in the order (A(J,I),J=1,I),I=1,N) on the matrix 
!             file, where first index is the row. Therefore, we need to transpose matrix file
!             storage to the MQC lower trangular matrix after reading for correct storage. 
!             This is only an issue for nonsymmetric matrices stored in symmetric form.
              matrixOut = transpose(matrixOut)
              deallocate(integerTmp)
            case('REAL-VECTOR')
              allocate(arrayTmp(LR))
              call Rd_RBuf(fileinfo%unitNumber,NTot,LenBuf,arrayTmp)
              if(Present(vectorOut)) then
                vectorOut = arrayTmp
              elseIf(Present(matrixOut)) then
                call MQC_Matrix_DiagMatrix_Put(matrixOut,arrayTmp)
              elseIf(Present(mqcVarOut)) then
                mqcVarOut = arrayTmp
              else
                call mqc_error_l('Reading vector from Gaussian matrix file, but NO VECTOR SENT to procedure.',  &
                  6,'Present(vectorOut)',Present(vectorOut),'Present(matrixOut)',Present(matrixOut) )
              endIf
              deallocate(arrayTmp)
            case('REAL-MATRIX')
              allocate(arrayTmp(LR))
              call Rd_RBuf(fileinfo%unitNumber,NTot,LenBuf,arrayTmp)
              if(Present(matrixOut)) then
                matrixOut = Reshape(arrayTmp,[N1,N2])
              elseIf(Present(mqcVarOut)) then
                mqcVarOut = Reshape(arrayTmp,[N1,N2])
              else
                call mqc_error_l('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to procedure.',  &
                  6,'Present(mqcVarOut)',Present(mqcVarOut),'Present(matrixOut)',Present(matrixOut))
              endIf
              deallocate(arrayTmp)
            case('REAL-SYMMATRIX')
              allocate(arrayTmp(LR))
              call Rd_RBuf(fileinfo%unitNumber,NTot,LenBuf,arrayTmp)
              if(Present(matrixOut)) then
                call MQC_Matrix_SymmMatrix_Put(matrixOut,arrayTmp)
              elseIf(Present(mqcVarOut)) then
                mqcVarOut = mqc_matrixSymm2Full(arrayTmp,'U')
              else
                call mqc_error_l('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to procedure.',  &
                  6,'Present(mqcVarOut)',Present(mqcVarOut),'Present(matrixOut)',Present(matrixOut))
              endIf
              deallocate(arrayTmp)
            case('REAL-ASYMMATRIX')
              allocate(arrayTmp(LR))
              call Rd_RBuf(fileinfo%unitNumber,NTot,LenBuf,arrayTmp)
!             Triangular matrices are stored in the order (A(J,I),J=1,I),I=1,N) on the matrix 
!             file, where first index is the row. Therefore, we need to transpose matrix file
!             storage to the MQC lower trangular matrix after reading for correct storage. 
!             This is only an issue for nonsymmetric matrices stored in symmetric form.
              if(Present(matrixOut)) then
                call MQC_Matrix_SymmMatrix_Put(matrixOut,arrayTmp,'antisymmetric')
!               Matrix files have either symmetric/hermitian storage or antisymmetric/
!               anthermitian storage. MQC currently has a symmetric only storage for both real 
!               and complex parts so make nonsymmetric matrices square.
!                call mqc_matrix_symm2full(matrixOut,'antisymmetric')
                matrixOut = transpose(matrixOut)
              elseIf(Present(mqcVarOut)) then
                mqcVarOut = mqc_matrixSymm2Full(arrayTmp,'U')
                mqcVarOut = transpose(mqcVarOut)
              else
                call mqc_error_l('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to procedure.',  &
                  6,'Present(mqcVarOut)',Present(mqcVarOut),'Present(matrixOut)',Present(matrixOut))
              endIf
              deallocate(arrayTmp)
            case('REAL-SYMMATRIXN')
              allocate(arrayTmp(LR))
              call Rd_RBuf(fileinfo%unitNumber,NTot,LenBuf,arrayTmp)
              if(Present(matrixOut)) then
                call MQC_Matrix_SymmMatrix_Put(matrixOut,arrayTmp((myArrayNum-1)*(LR/N3)+1:myArrayNum*(LR/N3)))
              elseIf(Present(mqcVarOut)) then
                mqcVarOut = mqc_matrixSymm2Full(arrayTmp((myArrayNum-1)*(LR/N3)+1:myArrayNum*(LR/N3)),'U')
              else
                call mqc_error_l('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to procedure.',  &
                  6,'Present(mqcVarOut)',Present(mqcVarOut),'Present(matrixOut)',Present(matrixOut))
              endIf
              deallocate(arrayTmp)
            case('REAL-ASYMMATRIXN')
              allocate(arrayTmp(LR))
              call Rd_RBuf(fileinfo%unitNumber,NTot,LenBuf,arrayTmp)
              if(Present(matrixOut)) then
!                call MQC_Matrix_SymmMatrix_Put(matrixOut,arrayTmp((myArrayNum-1)*(LR/N3)+1:myArrayNum*(LR/N3)),'antisymmetric')
                call MQC_Matrix_SymmMatrix_Put(matrixOut,arrayTmp((myArrayNum-1)*(LR/N3)+1:myArrayNum*(LR/N3)))
              elseIf(Present(mqcVarOut)) then
                mqcVarOut = mqc_matrixSymm2Full(arrayTmp((myArrayNum-1)*(LR/N3)+1:myArrayNum*(LR/N3)),'U')
              else
                call mqc_error_l('Reading matrix from Gaussian matrix file, but NO MATRIX SENT to procedure.',  &
                  6,'Present(mqcVarOut)',Present(mqcVarOut),'Present(matrixOut)',Present(matrixOut))
              endIf
              deallocate(arrayTmp)
            case('COMPLEX-VECTOR')
              allocate(complexTmp(LR))
              call Rd_CBuf(fileinfo%unitNumber,NTot,LenBuf,complexTmp)
              if(Present(vectorOut)) then
                vectorOut = complexTmp
              elseIf(Present(matrixOut)) then
                call MQC_Matrix_DiagMatrix_Put(matrixOut,complexTmp)
              else
                call mqc_error_l('Reading vector from Gaussian matrix file, but NO VECTOR SENT &
                  & to procedure.', 6, &
                  'Present(vectorOut)', Present(vectorOut), &
                  'Present(matrixOut)', Present(matrixOut) )
              endIf
              deallocate(complexTmp)
            case('COMPLEX-MATRIX')
              if(.not.Present(matrixOut)) call mqc_error_l('Reading matrix from Gaussian matrix &
                & file, but NO MATRIX SENT to procedure.', 6, &
                'Present(matrixOut)', Present(matrixOut) )
              allocate(complexTmp(LR))
 !             write(*,1060) 'reading matrix'
              call Rd_CBuf(fileinfo%unitNumber,NTot,LenBuf,complexTmp)
 !             write(*,1060) 'read matrix'
              matrixOut = Reshape(complexTmp,[N1,N2])
              deallocate(complexTmp)
            case('COMPLEX-SYMMATRIX')
              if(.not.Present(matrixOut)) call mqc_error_l('Reading matrix from Gaussian matrix &
                & file, but NO MATRIX SENT to procedure.', 6, &
                'Present(matrixOut)', Present(matrixOut) )
              allocate(complexTmp(LR))
              call Rd_CBuf(fileinfo%unitNumber,NTot,LenBuf,complexTmp)
              call MQC_Matrix_SymmMatrix_Put(matrixOut,complexTmp,'hermitian')
!             Matrix files have either symmetric/hermitian storage or antisymmetric/
!             anthermitian storage. MQC currently has a symmetric only storage for both real 
!             and complex parts so make nonsymmetric matrices square.
!              call mqc_matrix_symm2full(matrixOut,'hermitian')
!             Triangular matrices are stored in the order (A(J,I),J=1,I),I=1,N) on the matrix 
!             file, where first index is the row. Therefore, we need to transpose matrix file
!             storage to the MQC lower trangular matrix after reading for correct storage. 
!             This is only an issue for nonsymmetric matrices stored in symmetric form.
              matrixOut = transpose(matrixOut)
              deallocate(complexTmp)
            case('COMPLEX-ASYMMATRIX')
              if(.not.Present(matrixOut)) call mqc_error_l('Reading matrix from Gaussian matrix &
                & file, but NO MATRIX SENT to procedure.', 6, &
                'Present(matrixOut)', Present(matrixOut) )
              allocate(complexTmp(LR))
              call Rd_CBuf(fileinfo%unitNumber,NTot,LenBuf,complexTmp)
              call MQC_Matrix_SymmMatrix_Put(matrixOut,complexTmp,'antihermitian')
!             Triangular matrices are stored in the order (A(J,I),J=1,I),I=1,N) on the matrix 
!             file, where first index is the row. Therefore, we need to transpose matrix file
!             storage to the MQC lower trangular matrix after reading for correct storage. 
!             This is only an issue for nonsymmetric matrices stored in symmetric form.
              matrixOut = transpose(matrixOut)
              deallocate(complexTmp)
            case('COMPLEX-SYMMATRIXN')
              allocate(complexTmp(LR))
              call Rd_CBuf(fileinfo%unitNumber,NTot,LenBuf,complexTmp)
              call MQC_Matrix_SymmMatrix_Put(matrixOut,complexTmp((myArrayNum-1)*(LR/N3)+1:myArrayNum*(LR/N3)),'hermitian')
              matrixOut = transpose(matrixOut)
              deallocate(complexTmp)
            case('COMPLEX-ASYMMATRIXN')
              allocate(complexTmp(LR))
              call Rd_CBuf(fileinfo%unitNumber,NTot,LenBuf,complexTmp)
              call MQC_Matrix_SymmMatrix_Put(matrixOut,complexTmp((myArrayNum-1)*(LR/N3)+1:myArrayNum*(LR/N3)),'antihermitian')
              matrixOut = transpose(matrixOut)
              deallocate(complexTmp)
            case('MIXED')
              write(*,1020)
              write(*,1040)' LR   = ',LR
              write(*,1040)' NR   = ',NR
              write(*,1040)' NI   = ',NI
              write(*,1040)' NRI  = ',NRI
              write(*,1040)' NTot = ',NTot
              write(*,1020)
              call mqc_error_a('No general way to load mixed types as of yet &
      &         We are doing it case-by-case at the moment and this does not match.', 6, &
      'MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI,ASym)', &
      MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI,ASym) )
            case('2ERIS-SYMSYMR4TENSOR')
              if(Present(r4TensorOut)) then
                if(NRI.eq.1) then
                  allocate(arrayTmp(LR*NR))
                  call Rd_2EN(fileinfo%unitNumber,NR,LR,NR*LR,NTot,LenBuf,arrayTmp)
                  call MQC_Matrix_SymmSymmR4Tensor_Put_Real(r4TensorOut,&
                    arrayTmp((myArrayNum-1)*LR+1:myArrayNum*LR))
                  deallocate(arrayTmp)
                elseIf(NRI.eq.2) then
                  allocate(arrayTmp(NR*LR*2))
                  call Rd_2EN(fileinfo%unitNumber,NR,LR,NR*LR,2*NTot,2*LenBuf,arrayTmp)
                  complexTmp = reshape(arrayTmp,shape(complexTmp))
                  call MQC_Matrix_SymmSymmR4Tensor_Put_Complex(r4TensorOut, &
                    complexTmp((myArrayNum-1)*LR+1:myArrayNum*LR))
                  deallocate(arrayTmp,complexTmp)
                endIf
              elseIf(Present(mqcVarOut)) then
                if(NRI.eq.1) then
                  allocate(arrayTmp(LR*NR))
                  call Rd_2EN(fileinfo%unitNumber,NR,LR,NR*LR,NTot,LenBuf,arrayTmp)
                  nBasis = fileinfo%getVal('nbasis')
                  allocate(arraytmp2(nBasis,nBasis,nBasis,nBasis))
                  call mqc_packedSymmetricSymmetricR4Tensor2Full_real(  &
                    arrayTmp,arraytmp2)
                  mqcVarOut = arraytmp2
                  deAllocate(arraytmp2)
                  deallocate(arrayTmp)
                else
                  call mqc_error('Attempted to read complex rank-4 tensor into MQCVar type.')
                endIf
              else
                call mqc_error('Reading rank-4 tensor from Gaussian matrix file, but required &
                  & variable R4TENSOR or mqcVarOut not sent to procedure.')
              endIf
!
            case('REAL-SYMSYMR4TENSOR')
              if(Present(r4TensorOut)) then
                if(NRI.eq.1) then
                  allocate(arrayTmp(LR*NR))
                  call Rd_RBuf(fileinfo%unitNumber,NTot,LenBuf,arrayTmp)
                  call MQC_Matrix_SymmSymmR4Tensor_Put_Real(r4TensorOut,&
                    arrayTmp((myArrayNum-1)*LR+1:myArrayNum*LR))
                  deallocate(arrayTmp)
                else
                  call mqc_error('Attempted to read complex rank-4 tensor into MQC Rank4 type.')
                endIf
              elseIf(Present(mqcVarOut)) then
                if(NRI.eq.1) then
                  allocate(arrayTmp(LR*NR))
                  call Rd_RBuf(fileinfo%unitNumber,NTot,LenBuf,arrayTmp)
                  nBasis = ABS(N1)
                  allocate(arraytmp2(nBasis,nBasis,nBasis,nBasis))
                  call mqc_packedSymmetricSymmetricR4Tensor2Full_real(  &
                    arrayTmp,arraytmp2)
                  mqcVarOut = arraytmp2
                  deAllocate(arraytmp2)
                  deallocate(arrayTmp)
                else
                  call mqc_error('Attempted to read complex rank-4 tensor into MQCVar type.')
                endIf
              else
                call mqc_error('Reading rank-4 tensor from Gaussian matrix file, but required &
                  & variable R4TENSOR or mqcVarOut not sent to procedure.')
              endIf
!
            case('SCALARS-VECTOR') 
!             Just read Gaussian scalars into a real vector for now
              allocate(arrayTmp(NR*LR))
              call Rd_RInd(fileinfo%unitNumber,NR,LR,NTot,LenBuf,LNZ,arrayTmp)
              if(PRESENT(vectorOut)) then
                vectorOut = arrayTmp
              elseIf(PRESENT(mqcVarOut)) then
                mqcVarOut = arrayTmp
              else
                call mqc_error('SCALARS-VECTOR found in GMatrix array read without a vector to return it in.')
              endIf
              deallocate(arrayTmp)
!
            case default
              write(*,1050)' Matrix type: ',Trim(MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI,ASym))
              call mqc_error_A('Found strange matrix type in Gaussian matrix read routine.', 6, &
                   'MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI,ASym)', &
                   MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI,ASym) )
            end select
            found = .true.
            exit outerLoop
          elseIf(NTot.gt.0) then
            Call Rd_Skip(fileinfo%UnitNumber,NTot,LenBuf)
          endIf
          Call Rd_Labl(fileinfo%UnitNumber,IVers,cBuffer,NI,NR,NTot,LenBuf,  &
            N1,N2,N3,N4,N5,ASym,NRI,EOF)
          LR = LenArr(N1,N2,N3,N4,N5)
          EOF = EOF.or.cBuffer.eq.'END'
          if(DEBUG) write(IOut,1010) TRIM(cBuffer),NI,NR,NRI,NTot,LenBuf,  &
            N1,N2,N3,N4,N5,ASym,LR
        endDo
        if(i==1) then
          my_filename = TRIM(fileinfo%filename)
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            my_filename)
        endIf
      endDo outerLoop
      if(present(foundOut)) foundOut = found
      if(.not.found) then
        errorMsg = 'Could NOT find requested matrix file label "'//TRIM(label)//'".'
        if(present(foundOut)) then
          write(6,'(A)') errorMsg
        else
          call MQC_Error_L(errorMsg,6,'found',found)
        endIf
      endIf
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Read_Array


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Write_Array2
      subroutine MQC_Gaussian_Unformatted_Matrix_Write_Array2(fileinfo,  &
       label,mqcVariable,filename,storage)
!
!     This Routine is used to look-up a matrix in a unformatted matrix file and
!     write that array into either (OPTIONAL) output dummy MQC_Matrix argument
!     <matrixIn>, (OPTIONAL) output dummy MQC_Vector argument <vectorIn>, or
!     (OPTIONAL) output dummy MQC_R4Tensor argument <r4TensorIn>. The character
!     label for the array of interest is sent to this routine in dummy argument
!     <label>.
!
!     Dummy argument <filename> is optional and is only used if fileinfo
!     hasn't already been defined using Routine
!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!     filename sent is different from the filename associated with object
!     fileinfo.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!     called before calling this routine. The expectation is that
!     MQC_Gaussian_Unformatted_Matrix_Write_Header is also called before this
!     routine. However, it is also OK to call this routine first. In that case,
!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     L. M. Thompson, 2017.
!     H. P. Hratchian, 2022.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      type(MQC_Variable),intent(in)::mqcVariable
      character(len=*),intent(in),OPTIONAL::filename,storage
!
      integer::iout=6
!
!     Temporary local variables used when calling the gauopen routines.
      integer::LenBuf
      character(len=64)::tmpLabel
!
!     Local temp variables.
      integer::i,nInputArrays
      integer,parameter::IOne=1
      real,allocatable,dimension(:,:)::realMatrixTmp
      real,pointer,dimension(:)::realMatrixTmpV
      integer,allocatable,dimension(:,:)::intMatrixTmp
      complex(kind=real64),allocatable,dimension(:,:)::compMatrixTmp
      real,allocatable,dimension(:)::realVectorTmp
      integer,allocatable,dimension(:)::intVectorTmp
      complex(kind=real64),allocatable,dimension(:)::compVectorTmp
      type(MQC_Matrix)::matrixInUse 
      character(len=256)::my_filename,my_storage
      logical::DEBUG=.false.,ok
      Parameter(LenBuf=4000)
!
!     Format statements.
!
 1010 format(' Label ',A48,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
!
!
!     Begin by seeing if a new file or filename has been sent by the calling
!     program unit. If file has not been opened or there is no header data,
!     read the header data. If file is not in write mode, write the header
!     data and set to write mode
!
      if(present(storage)) then
        call String_Change_Case(storage,'l',my_storage)
      else
        my_storage = ''
      endIf
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
               'ok', ok )
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,filename)
        else
          call MQC_Error_L('Error reading Gaussian matrix file header (4): Must include a filename.', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
            'ok', ok )
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'W' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call fileinfo%OPENFILE(TRIM(my_filename),0,ok)
        if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
             'ok', ok )
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          my_filename)
      endIf
      if(.not.fileinfo%header_written) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call fileinfo%OPENFILE(TRIM(my_filename),0,ok)
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Now that the prelim checks are done, write the array sent by the calling
!     program unit to the matrix file.
!
      call String_Change_Case(label,'u',tmpLabel)
      select case(MQC_Variable_getTypeCode(mqcVariable))
!
!     Integer Array
      case(3)
        select case(RANK(mqcVariable))
!
!       vector...
        case(1)
          intVectorTmp = mqcVariable
          call wr_LIBuf(fileinfo%UnitNumber,TRIM(tmpLabel),IOne,LenBuf,  &
            SIZE(mqcVariable,1),0,0,0,0,.False.,intVectorTmp)
          DeAllocate(intVectorTmp)
!
!       matrix...
        case(2)
          Allocate(intVectorTmp(SIZE(mqcVariable,1)*SIZE(mqcVariable,2)))
          intMatrixTmp = mqcVariable
          intVectorTmp = RESHAPE(intMatrixTmp,[SIZE(intMatrixTmp)])
          call wr_LIBuf(fileinfo%UnitNumber,TRIM(tmpLabel),IOne,LenBuf,  &
            SIZE(mqcVariable,1),SIZE(mqcVariable,2),0,0,0,.False.,  &
            intVectorTmp)
          DeAllocate(intMatrixTmp)
          DeAllocate(intVectorTmp)
        case default
          call mqc_error('Invalid mqcVariable rank sent to '//  &
            'MQC_Gaussian_Unformatted_Matrix_Write_Array2.')
        end select
!
!     Real Array
      case(2)
        select case(RANK(mqcVariable))
!
!       vector...
        case(1)
          realVectorTmp = mqcVariable
          call wr_LRBuf(fileinfo%UnitNumber,TRIM(tmpLabel),IOne,LenBuf,  &
            SIZE(mqcVariable,1),0,0,0,0,.False.,realVectorTmp)
          DeAllocate(realVectorTmp)
!
!       matrix...
        case(2)
          Allocate(realVectorTmp(SIZE(mqcVariable,1)*SIZE(mqcVariable,2)))
          realMatrixTmp = mqcVariable
          realVectorTmp = RESHAPE(realMatrixTmp,[SIZE(realMatrixTmp)])
          call wr_LRBuf(fileinfo%UnitNumber,TRIM(tmpLabel),IOne,LenBuf,  &
            SIZE(mqcVariable,1),SIZE(mqcVariable,2),0,0,0,.False.,  &
            realVectorTmp)
          DeAllocate(realMatrixTmp)
          DeAllocate(realVectorTmp)
        case default
          call mqc_error('Invalid mqcVariable rank sent to '//  &
            'MQC_Gaussian_Unformatted_Matrix_Write_Array2.')
        end select
!
!     Unknown Array Type
      case default
        call mqc_error('Invalid mqcVariable type code found in '//  &
          'MQC_Gaussian_Unformatted_Matrix_Write_Array2.')
      end select
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Write_Array2


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Write_Array
      subroutine MQC_Gaussian_Unformatted_Matrix_Write_Array(fileinfo,  &
       label,matrixIn,vectorIn,r4TensorIn,filename,storage)
!
!     This Routine is used to look-up a matrix in a unformatted matrix file and
!     write that array into either (OPTIONAL) output dummy MQC_Matrix argument
!     <matrixIn>, (OPTIONAL) output dummy MQC_Vector argument <vectorIn>, or
!     (OPTIONAL) output dummy MQC_R4Tensor argument <r4TensorIn>. The character
!     label for the array of interest is sent to this routine in dummy argument
!     <label>.
!
!     Dummy argument <filename> is optional and is only used if fileinfo
!     hasn't already been defined using Routine
!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!     filename sent is different from the filename associated with object
!     fileinfo.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!     called before calling this routine. The expectation is that
!     MQC_Gaussian_Unformatted_Matrix_Write_Header is also called before this
!     routine. However, it is also OK to call this routine first. In that case,
!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      type(MQC_Matrix),intent(in),OPTIONAL::matrixIn 
      type(MQC_Vector),intent(in),OPTIONAL::vectorIn 
      type(MQC_R4Tensor),intent(in),OPTIONAL::r4TensorIn 
      character(len=*),intent(in),OPTIONAL::filename,storage
!
      integer::iout=6
!
!     Temporary local variables used when calling the gauopen routines.
      integer::LenBuf
      character(len=64)::tmpLabel
!
!     Local temp variables.
      integer::i,nInputArrays
      integer::Ione
      real,allocatable,dimension(:,:)::realMatrixTmp
      real,pointer,dimension(:)::realMatrixTmpV
      integer,allocatable,dimension(:,:)::intMatrixTmp
      complex(kind=real64),allocatable,dimension(:,:)::compMatrixTmp
      real,allocatable,dimension(:)::realVectorTmp
      integer,allocatable,dimension(:)::intVectorTmp
      complex(kind=real64),allocatable,dimension(:)::compVectorTmp
      type(MQC_Matrix)::matrixInUse 
      character(len=256)::my_filename,my_storage
      logical::DEBUG=.false.,ok
      Parameter(LenBuf=4000)
      Real(kind=real64)::Thresh=1e-10
!
!     Format statements.
!
 1010 format(' Label ',A48,' NI=',I2,' NR=',I2,' NRI=',I1,' NTot=',  &
        I8,' LenBuf=',I8,' N=',5I6,' ASym=',L1,' LR=',I5)
!
!
!     Begin by seeing if a new file or filename has been sent by the calling
!     program unit. If file has not been opened or there is no header data,
!     read the header data. If file is not in write mode, write the header
!     data and set to write mode
!
      if(present(storage)) then
        call String_Change_Case(storage,'l',my_storage)
      else
        my_storage = ''
      endIf
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
            'ok', ok )
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,filename)
        else
          call MQC_Error_L('Error reading Gaussian matrix file header (5): Must include a filename.', 6, &
            'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
            'ok', ok )
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'W' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call fileinfo%OPENFILE(TRIM(my_filename),0,ok)
        if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
             'ok', ok )
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          my_filename)
      endIf
      if(.not.fileinfo%header_written) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call fileinfo%OPENFILE(TRIM(my_filename),0,ok)
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Ensure that one and only one input MQC-type array has been sent from the
!     calling program unit.
!
      nInputArrays = 0
      if(Present(matrixIn)) nInputArrays = nInputArrays+1
      if(Present(vectorIn)) nInputArrays = nInputArrays+1
      if(Present(r4TensorIn)) nInputArrays = nInputArrays+1
      if(nInputArrays.ne.1) call mqc_error_i('Too many input arrays sent to Gaussian matrix file writing procedure.', 6, &
           'nInputArrays', nInputArrays )
!
!     Load the MQC variable into a regular array and get the required dimensions
!     Some routines will change when we upgrade to algebra2
!
      Ione = 1
      call String_Change_Case(label,'u',tmpLabel)
      If(present(matrixIn)) then
        matrixInUse = matrixIn
        if(mqc_matrix_haveReal(matrixInUse)) then 
          if((mqc_matrix_test_diagonal(matrixInUse).and.(my_storage.eq.'')).or.(my_storage.eq.'diag')) then
            if(.not.mqc_matrix_haveDiagonal(matrixInUse)) then
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Diag(matrixInUse)
              if(mqc_matrix_haveSymmetric(matrixInUse)) call mqc_matrix_symm2Diag(matrixInUse)
            endIf
            if(mqc_matrix_rows(matrixInUse).lt.mqc_matrix_columns(matrixInUse)) then
              allocate(realMatrixTmp(mqc_matrix_rows(matrixInUse),1))
            else
              allocate(realMatrixTmp(mqc_matrix_columns(matrixInUse),1))
            endIf
            allocate(realVectorTmp(size(realMatrixTmp,1)*size(realMatrixTmp,2)))
            realMatrixTmp = matrixInUse
            realVectorTmp = reshape(realMatrixTmp, shape(realVectorTmp))
            call wr_LRBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,size(realMatrixTmp,1), &
              0,0,0,0,.False.,realVectorTmp)
          elseIf((mqc_matrix_test_symmetric(matrixInUse).and.(my_storage.eq.'')).or.(my_storage.eq.'symm')) then
            if(.not.mqc_matrix_haveSymmetric(matrixInUse)) then
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Symm(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Symm(matrixInUse)
            endIf
            allocate(realMatrixTmp((mqc_matrix_rows(matrixInUse)*(mqc_matrix_rows(matrixInUse)+1))/2,1))
            allocate(realVectorTmp(size(realMatrixTmp,1)))
            realMatrixTmp = matrixInUse
            realVectorTmp = reshape(realMatrixTmp, shape(realVectorTmp))
            call wr_LRBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,-mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.False.,realVectorTmp)
          elseIf((mqc_matrix_test_symmetric(matrixInUse,'antisymmetric').and.(my_storage.eq.'')).or.(my_storage.eq.'asymm')) then
!           We store triangular matrices in the order (A(J,I),J=1,I),I=1,N) on the matrix file,
!           where first index is the row. Therefore, we need to transpose the MQC lower
!           trangular matrix before writing for correct matrix file storage. This is only an
!           issue for nonsymmetric matrices stored in LT form.
            realMatrixTmp = transpose(matrixInUse)
            if(.not.mqc_matrix_haveSymmetric(matrixInUse)) then
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Symm(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Symm(matrixInUse)
            endIf
            allocate(realMatrixTmp((mqc_matrix_rows(matrixInUse)*(mqc_matrix_rows(matrixInUse)+1))/2,1))
            allocate(realVectorTmp(size(realMatrixTmp,1)))
            realVectorTmp = reshape(realMatrixTmp, shape(realVectorTmp))
            call wr_LRBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,-mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.True.,realVectorTmp)
          elseIf((mqc_matrix_haveFull(matrixInUse).and.(my_storage.eq.'')).or.(my_storage.eq.'full')) then
            if(.not.mqc_matrix_haveFull(matrixInUse)) then
              if(mqc_matrix_haveSymmetric(matrixInUse)) call mqc_matrix_symm2Full(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Full(matrixInUse)
            endIf
            allocate(realMatrixTmp(mqc_matrix_rows(matrixInUse),mqc_matrix_columns(matrixInUse)))
            allocate(realVectorTmp(size(realMatrixTmp,1)*size(realMatrixTmp,2)))
            realMatrixTmp = matrixInUse
            realVectorTmp = reshape(realMatrixTmp, shape(realVectorTmp))
            call wr_LRBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.False.,realVectorTmp)
          else
            call mqc_error_l('type not recognised', 6, &
                 'mqc_matrix_test_diagonal(matrixInUse)', mqc_matrix_test_diagonal(matrixInUse), &
                 'mqc_matrix_test_symmetric(matrixInUse)', mqc_matrix_test_symmetric(matrixInUse), &
                 "mqc_matrix_test_symmetric(matrixInUse,'antisymmetric')", mqc_matrix_test_symmetric(matrixInUse,'antisymmetric'), &
                 'mqc_matrix_haveFull(matrixInUse)', mqc_matrix_haveFull(matrixInUse) )
          endIf
        elseIf(mqc_matrix_haveInteger(matrixInUse)) then 
          if((mqc_matrix_test_diagonal(matrixInUse).and.(my_storage.eq.'')).or.(my_storage.eq.'diag')) then
            if(.not.mqc_matrix_haveDiagonal(matrixInUse)) then
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Diag(matrixInUse)
              if(mqc_matrix_haveSymmetric(matrixInUse)) call mqc_matrix_symm2Diag(matrixInUse)
            endIf
            if(mqc_matrix_rows(matrixInUse).lt.mqc_matrix_columns(matrixInUse)) then
              allocate(intMatrixTmp(mqc_matrix_rows(matrixInUse),1))
            else
              allocate(intMatrixTmp(mqc_matrix_columns(matrixInUse),1))
            endIf
            intMatrixTmp = matrixInUse
             allocate(intVectorTmp(size(intMatrixTmp,1)*size(intMatrixTmp,2)))
            intVectorTmp = reshape(intMatrixTmp, shape(intVectorTmp))
            call wr_LIBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,size(intMatrixTmp,1), &
              0,0,0,0,.False.,intVectorTmp)
          elseIf((mqc_matrix_test_symmetric(matrixInUse).and.(my_storage.eq.'')).or.(my_storage.eq.'symm')) then
            if(.not.mqc_matrix_haveSymmetric(matrixInUse)) then
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Symm(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Symm(matrixInUse)
            endIf
            allocate(intMatrixTmp((mqc_matrix_rows(matrixInUse)*(mqc_matrix_rows(matrixInUse)+1))/2,1))
            allocate(intVectorTmp(size(intMatrixTmp,1)))
            intMatrixTmp = matrixInUse
            intVectorTmp = reshape(intMatrixTmp, shape(intVectorTmp))
            call wr_LIBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,-mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.False.,intVectorTmp)
          elseIf((mqc_matrix_test_symmetric(matrixInUse,'antisymmetric').and.(my_storage.eq.'')).or.(my_storage.eq.'asymm')) then
!           We store triangular matrices in the order (A(J,I),J=1,I),I=1,N) on the matrix file,
!           where first index is the row. Therefore, we need to transpose the MQC lower
!           trangular matrix before writing for correct matrix file storage. This is only an
!           issue for nonsymmetric matrices stored in LT form.
            intMatrixTmp = transpose(matrixInUse)
            if(.not.mqc_matrix_haveSymmetric(matrixInUse)) then
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Symm(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Symm(matrixInUse)
            endIf
            allocate(intMatrixTmp((mqc_matrix_rows(matrixInUse)*(mqc_matrix_rows(matrixInUse)+1))/2,1))
            allocate(intVectorTmp(size(intMatrixTmp,1)))
            intVectorTmp = reshape(intMatrixTmp, shape(intVectorTmp))
            call wr_LIBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,-mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.True.,intVectorTmp)
          elseIf((mqc_matrix_haveFull(matrixInUse).and.(my_storage.eq.'')).or.(my_storage.eq.'full')) then
            if(.not.mqc_matrix_haveFull(matrixInUse)) then
              if(mqc_matrix_haveSymmetric(matrixInUse)) call mqc_matrix_symm2Full(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Full(matrixInUse)
            endIf
            allocate(intMatrixTmp(mqc_matrix_rows(matrixInUse),mqc_matrix_columns(matrixInUse)))
            intMatrixTmp = matrixInUse
            allocate(intVectorTmp(size(intMatrixTmp,1)*size(intMatrixTmp,2)))
            intVectorTmp = reshape(intMatrixTmp, shape(intVectorTmp))
            call wr_LIBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.False.,intVectorTmp)
          else
            call mqc_error_l('type not recognised', 6, &
                 'mqc_matrix_test_diagonal(matrixInUse)', mqc_matrix_test_diagonal(matrixInUse), &
                 'mqc_matrix_test_symmetric(matrixInUse)', mqc_matrix_test_symmetric(matrixInUse), &
                 "mqc_matrix_test_symmetric(matrixInUse,'antisymmetric')", mqc_matrix_test_symmetric(matrixInUse,'antisymmetric'), &
                 'mqc_matrix_haveFull(matrixInUse)', mqc_matrix_haveFull(matrixInUse) )
          endIf
        elseIf(mqc_matrix_haveComplex(matrixInUse)) then 
!       There is a bug in Wr_LCBuf in qcmatrix.F of gauopen if you want to write complex.
!       NR should only be negative in Wr_Labl and positive everywhere else.
!       Please recomplile MQC after making these changes.
          Ione = -1
          if((mqc_matrix_test_diagonal(matrixInUse).and.(my_storage.eq.'')).or.(my_storage.eq.'diag')) then
            if(.not.mqc_matrix_haveDiagonal(matrixInUse)) then
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Diag(matrixInUse)
              if(mqc_matrix_haveSymmetric(matrixInUse)) call mqc_matrix_symm2Diag(matrixInUse)
            endIf
            if(mqc_matrix_rows(matrixInUse).lt.mqc_matrix_columns(matrixInUse)) then
              allocate(compMatrixTmp(mqc_matrix_rows(matrixInUse),1))
            else
              allocate(compMatrixTmp(mqc_matrix_columns(matrixInUse),1))
            endIf
            allocate(compVectorTmp(size(compMatrixTmp,1)*size(compMatrixTmp,2)))
            compMatrixTmp = matrixInUse
            compVectorTmp = reshape(compMatrixTmp, shape(compVectorTmp))
            call wr_LCBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,size(compMatrixTmp,1), &
              0,0,0,0,.False.,compVectorTmp)
          elseIf((mqc_matrix_test_symmetric(matrixInUse,'hermitian',.true.).and.(my_storage.eq.'')).or.&
              (my_storage.eq.'herm'.and.MQC_Vector_Norm(matrixInUse%diagonal()).gt.Thresh) &
              .or.(my_storage.eq.'symm')) then
!           We store triangular matrices in the order (A(J,I),J=1,I),I=1,N) on the matrix file,
!           where first index is the row. Therefore, we need to transpose the MQC lower
!           trangular matrix before writing for correct matrix file storage. This is only an
!           issue for nonsymmetric matrices stored in LT form.
            if(.not.mqc_matrix_haveSymmetric(matrixInUse)) then
              matrixInUse = transpose(matrixInUse)
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Symm(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Symm(matrixInUse)
            endIf
            allocate(compMatrixTmp((mqc_matrix_rows(matrixInUse)*(mqc_matrix_rows(matrixInUse)+1))/2,1))
            allocate(compVectorTmp(size(compMatrixTmp,1)))
            compMatrixTmp = matrixInUse
            compVectorTmp = reshape(compMatrixTmp, shape(compVectorTmp))
            call wr_LCBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,-mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.False.,compVectorTmp)
          elseIf((mqc_matrix_test_symmetric(matrixInUse,'antihermitian').and.(my_storage.eq.'')).or.(my_storage.eq.'aher') &
              .or.(my_storage.eq.'asym')) then
!           We store triangular matrices in the order (A(J,I),J=1,I),I=1,N) on the matrix file,
!           where first index is the row. Therefore, we need to transpose the MQC lower
!           trangular matrix before writing for correct matrix file storage. This is only an
!           issue for nonsymmetric matrices stored in LT form.
            if(.not.mqc_matrix_haveSymmetric(matrixInUse)) then
              matrixInUse = transpose(matrixInUse)
              if(mqc_matrix_haveFull(matrixInUse)) call mqc_matrix_full2Symm(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Symm(matrixInUse)
            endIf
            allocate(compMatrixTmp((mqc_matrix_rows(matrixInUse)*(mqc_matrix_rows(matrixInUse)+1))/2,1))
            allocate(compVectorTmp(size(compMatrixTmp,1)))
            compMatrixTmp = matrixInUse
            compVectorTmp = reshape(compMatrixTmp, shape(compVectorTmp))
            call wr_LCBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,-mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.True.,compVectorTmp)
          elseIf((mqc_matrix_haveFull(matrixInUse).and.(my_storage.eq.'')).or.(my_storage.eq.'full')) then
            if(.not.mqc_matrix_haveFull(matrixInUse)) then
              if(mqc_matrix_haveSymmetric(matrixInUse)) call mqc_matrix_symm2Full(matrixInUse)
              if(mqc_matrix_haveDiagonal(matrixInUse)) call mqc_matrix_diag2Full(matrixInUse)
            endIf
            allocate(compMatrixTmp(mqc_matrix_rows(matrixInUse),mqc_matrix_columns(matrixInUse)))
            allocate(compVectorTmp(size(compMatrixTmp,1)*size(compMatrixTmp,2)))
            compMatrixTmp = matrixInUse
            compVectorTmp = reshape(compMatrixTmp, shape(compVectorTmp))
            call wr_LCBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,mqc_matrix_rows(matrixInUse), &
              mqc_matrix_columns(matrixInUse),0,0,0,.False.,compVectorTmp)
          else
            call mqc_error_l('type not recognised', 6, &
                 'mqc_matrix_test_diagonal(matrixInUse)', mqc_matrix_test_diagonal(matrixInUse), &
                 'mqc_matrix_test_symmetric(matrixInUse)', mqc_matrix_test_symmetric(matrixInUse), &
                 "mqc_matrix_test_symmetric(matrixInUse,'hermitian')", mqc_matrix_test_symmetric(matrixInUse,'hermitian'), &
                 'mqc_matrix_haveFull(matrixInUse)', mqc_matrix_haveFull(matrixInUse) )
          endIf
        else
          call mqc_error_l('MatrixIn type not recognised in &
            &MQC_Gaussian_Unformatted_Matrix_Write_Array', 6, &
            'mqc_matrix_haveReal(matrixInUse)', mqc_matrix_haveReal(matrixInUse), &
            'mqc_matrix_haveInteger(matrixInUse)', mqc_matrix_haveInteger(matrixInUse), &
            'mqc_matrix_haveComplex(matrixInUse)', mqc_matrix_haveComplex(matrixInUse) )
        endIf
      elseIf(present(vectorIn)) then
        if(mqc_vector_haveReal(vectorIn)) then 
          allocate(realVectorTmp(mqc_length_vector(vectorIn)))
          realVectorTmp = vectorIn
          call wr_LRBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,mqc_length_vector(vectorIn), &
            0,0,0,0,.False.,realVectorTmp)
        elseIf(mqc_vector_haveInteger(vectorIn)) then 
          allocate(intVectorTmp(mqc_length_vector(vectorIn)))
          intVectorTmp = vectorIn
          call wr_LIBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,mqc_length_vector(vectorIn), &
            0,0,0,0,.False.,intVectorTmp)
        elseIf(mqc_vector_haveComplex(vectorIn)) then 
          IOne = -1
          allocate(compVectorTmp(mqc_length_vector(vectorIn)))
          compVectorTmp = vectorIn
          call wr_LCBuf(fileinfo%UnitNumber,tmpLabel,Ione,LenBuf,mqc_length_vector(vectorIn), &
            0,0,0,0,.False.,compVectorTmp)
        else
          call mqc_error_l('VectorIn type not recognised in &
            &MQC_Gaussian_Unformatted_Matrix_Write_Array', 6, &
            'mqc_vector_haveReal(vectorIn)', mqc_vector_haveReal(vectorIn), &
            'mqc_vector_haveInteger(vectorIn)', mqc_vector_haveInteger(vectorIn), &
            'mqc_vector_haveComplex(vectorIn)', mqc_vector_haveComplex(vectorIn) )
        endIf
      endIf
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Write_Array
!
!
!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Atom_Info
      Function MQC_Gaussian_Unformatted_Matrix_Get_Atom_Info(fileinfo,element,label)
!
!     This function is used to get info about specific atom information
!     associated with the Gaussian unformatted matrix file sent in object
!     fileinfo.
!
!     Input argument element refers to a specific atom in the molecule some
!     other element related to the info requested by input argument label.
!
!     The recognized labels and their meaning include:
!           'nuclearCharge'   return the atomic charge of atom number <element>.
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      integer::element
      character(len=*),intent(in)::label
      integer::MQC_Gaussian_Unformatted_Matrix_Get_Atom_Info
      integer::value_out=0
      character(len=64)::myLabel
      character(len=256)::my_filename
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen())  &
        call MQC_Error_L('Failed to retrieve atom info from Gaussian matrix file: File not open.', 6, &
        'fileinfo%isOpen()', fileinfo%isOpen() )
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('atomiccharge','nuclearcharge')
        if((element.le.0).or.(element.gt.fileinfo%natoms))  &
          call MQC_Error_I('element to %getAtomInfo is invalid.', 6, &
          'element', element, &
          'fileinfo%natoms', fileinfo%natoms )
        value_out = fileinfo%atomicCharges(element)
      case default
        call mqc_error_A('Invalid label sent to %getAtomInfo.', 6, &
             'mylabel', mylabel )
      endSelect
!
      MQC_Gaussian_Unformatted_Matrix_Get_Atom_Info = value_out
      return
      end Function MQC_Gaussian_Unformatted_Matrix_Get_Atom_Info


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Numbers
      Function MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Numbers(fileinfo)  &
        Result(arrayOut)
!
!     This function is used to get the array of atomic numbers from the Gaussian
!     matrix file corresponding to argument fileinfo.
!
!
!     H. P. Hratchian, 2020.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      integer(kind=int64),dimension(:),allocatable::arrayOut
      character(len=256)::my_filename
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen())  &
        call MQC_Error_L('Failed to retrieve atomic numbers array from Gaussian matrix file: File not open.', 6, &
        'fileinfo%isOpen()', fileinfo%isOpen() )
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      if(.not.allocated(fileinfo%atomicNumbers))  &
        call MQC_Error_L('Atomic numbers requested, but NOT available.', 6, &
        'allocated(fileinfo%atomicNumbers)', allocated(fileinfo%cartesians))
      allocate(arrayOut(fileinfo%natoms))
      arrayOut = fileinfo%atomicNumbers
!
      return
      end Function MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Numbers


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Carts
      Function MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Carts(fileinfo)  &
        Result(arrayOut)
!
!     This function is used to get the array of atomic Cartesian coordinates
!     from the Gaussian matrix file corresponding to argument fileinfo.
!
!
!     H. P. Hratchian, 2019.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      real(kind=real64),dimension(:),allocatable::arrayOut
      character(len=256)::my_filename
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen())  &
        call MQC_Error_L('Failed to atomic cartesian coordinates from Gaussian matrix file: File not open.', 6, &
        'fileinfo%isOpen()', fileinfo%isOpen() )
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      if(.not.allocated(fileinfo%cartesians))  &
        call MQC_Error_L('Atomic Cartesian coordinates requestion, but NOT available.', 6, &
        'allocated(fileinfo%cartesians)', allocated(fileinfo%cartesians))
      allocate(arrayOut(3*fileinfo%natoms))
      arrayOut = fileinfo%cartesians
!
      return
      end Function MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Carts


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Weights
      Function MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Weights(fileinfo)  &
        Result(arrayOut)
!
!     This function is used to get the array of atomic weights/masses from the
!     Gaussian matrix file corresponding to argument fileinfo.
!
!
!     H. P. Hratchian, 2019.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      real(kind=int64),dimension(:),allocatable::arrayOut
      character(len=256)::my_filename
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen())  &
        call MQC_Error_L('Failed to retrieve basis info from Gaussian matrix file: File not open.', 6, &
        'fileinfo%isOpen()', fileinfo%isOpen() )
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      if(.not.allocated(fileinfo%atomicWeights))  &
        call MQC_Error_L('Atomic weights requestion, but NOT available.', 6, &
        'allocated(fileinfo%atomicWeights)', allocated(fileinfo%atomicWeights))
      allocate(arrayOut(fileinfo%natoms))
      arrayOut = fileinfo%atomicWeights
!
      return
      end Function MQC_Gaussian_Unformatted_Matrix_Get_Atomic_Weights



!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Element
      Function MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Element(fileinfo,element,label)
!
!     This function is used to get info about specific basis functions
!     associated with the Gaussian unformatted matrix file sent in object
!     fileinfo.
!
!     Input argument element refers to a specific basis function by number or
!     some other element related to the info requested by input argument label.
!
!     The recognized labels and their meaning include:
!           'basis2Atom'      return the atomic center number on which basis
!                             function <element> is centered.
!           'basis type'      return the atomic orbital basis type of basis
!                             function <element> as numerical label. 
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      integer::element
      character(len=*),intent(in)::label
      integer::MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Element
      integer::value_out=0
      character(len=64)::myLabel
      character(len=256)::my_filename
!
 1000 format(1x,'getBasisInfo: element,value=',I5,2x,I5)
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen())  &
        call MQC_Error_L('Failed to retrieve basis info from Gaussian matrix file: File not open.', 6, &
        'fileinfo%isOpen()', fileinfo%isOpen() )
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('basis2atom')
        if(.not.allocated(fileinfo%basisFunction2Atom))  &
          call MQC_Error_L('Requested basis2Atom not possible.', 6, &
          'allocated(fileinfo%basisFunction2Atom)', allocated(fileinfo%basisFunction2Atom) )
        if((element.le.0).or.(element.gt.fileinfo%nbasis))  &
          call MQC_Error_I('element to %getBasisInfo is invalid.', 6, &
          'element', element, &
          'fileinfo%nbasis', fileinfo%nbasis )
        value_out = fileinfo%basisFunction2Atom(element)
        if(MQC_Gaussian_DEBUG) write(*,1000)element,value_out
      case('basis type')
        if(.not.allocated(fileinfo%IBasisFunctionType))  &
          call MQC_Error_l('Requested basis type not possible.', 6, &
          'allocated(fileinfo%IBasisFunctionType)', allocated(fileinfo%IBasisFunctionType) )
        if((element.le.0).or.(element.gt.fileinfo%nbasis))  &
          call MQC_Error_I('element to %getBasisInfo is invalid.', 6, &
          'element', element, &
          'fileinfo%nbasis', fileinfo%nbasis )
        value_out = fileinfo%IBasisFunctionType(element)
      case default
        call mqc_error_a('Invalid label sent to %getBasisInfo.', 6, &
             'mylabel', mylabel )
      endSelect
!
      MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Element = value_out
      return
      end Function MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Element


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Array
      Function MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Array(fileinfo,label)  &
        Result(arrayOut)
!
!     This function is used to get info about specific basis functions
!     associated with the Gaussian unformatted matrix file sent in object
!     fileinfo. This function returns the full array requested, not just a
!     single element.
!
!     Input argument element refers to a specific basis function by number or
!     some other element related to the info requested by input argument label.
!
!     The recognized labels and their meaning include:
!           'basis2Atom'      return an integer array giving the atomic center
!                             number on which each basis function is centered.
!           'basis type'      return an integer array giving the atomic orbital
!                             basis type of each basis function as a numerical
!                             label. 
!
!
!     H. P. Hratchian, 2018.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      integer,dimension(:),allocatable::arrayOut
      integer::value_out=0
      character(len=64)::myLabel
      character(len=256)::my_filename
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen())  &
        call MQC_Error_L('Failed to retrieve basis info from Gaussian matrix file: File not open.', 6, &
        'fileinfo%isOpen()', fileinfo%isOpen() )
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('basis2atom')
        if(.not.allocated(fileinfo%basisFunction2Atom))  &
          call MQC_Error_L('Requested basis2Atom not possible.', 6, &
          'allocated(fileinfo%basisFunction2Atom)', allocated(fileinfo%basisFunction2Atom) )
        allocate(arrayOut(fileinfo%NBasis))
        arrayOut = fileinfo%basisFunction2Atom
      case('basis type')
        if(.not.allocated(fileinfo%IBasisFunctionType))  &
          call MQC_Error_l('Requested basis type not possible.', 6, &
          'allocated(fileinfo%IBasisFunctionType)', allocated(fileinfo%IBasisFunctionType) )
        allocate(arrayOut(fileinfo%NBasis))
        arrayOut = fileinfo%IBasisFunctionType
      case default
        call mqc_error_a('Invalid label sent to %getBasisInfo.', 6, &
             'mylabel', mylabel )
      endSelect
!
      return
      end Function MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Array

      
!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Basis_Type_To_Function_Info
      Subroutine MQC_Gaussian_Unformatted_Matrix_Basis_Type_To_Function_Info(basisType,  &
        angularMomentum,component,isCartesian)
!
!     This routine is used to interpret an element of the basis function type
!     list. The basis type value is the only input argument (see functions
!     MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Element and
!     MQC_Gaussian_Unformatted_Matrix_Get_Basis_Info_Array to get the basis type
!     value for a basis function) The output of this routine is the angular
!     momentum of the basis function, the component (x, y, z for l=1, for
!     example), and a logical that is TRUE for Cartesian basis functions and
!     FALSE for pure basis functions.
!
!
!     H. P. Hratchian, 2024.
!
!
!     Variable Declarations.
!
      implicit none

      integer(kind=int64),intent(in)::basisType
      integer,intent(out),optional::angularMomentum,component
      logical,intent(out),optional::isCartesian
      integer::myBasisType,myAngularMomentum,myComponent
!
!     Decode the basisType number...
!
      myBasisType = ABS(basisType)
      myAngularMomentum = myBasisType/1000
      myComponent = mod(myBasisType,1000)
!
      if(PRESENT(angularMomentum)) angularMomentum = myAngularMomentum
      if(PRESENT(component)) component = myComponent
      if(PRESENT(isCartesian)) isCartesian = basisType.gt.0
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Basis_Type_To_Function_Info




!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Basis_Set_Matrices
      Subroutine MQC_Gaussian_Unformatted_Matrix_Get_Basis_Set_Matrices(fileinfo,basisData,fileName)
!
!     This subroutine is used to get the arrays containing the basis set coefficients, exponents
!     and coordinates etc. and load then into the basis set object from the Gaussian unformatted 
!     matrix file sent in object fileinfo.
!
!     L. M. Thompson, 2023.
!
!
!     Variable Declarations.
!
      implicit none
      class(mqc_gaussian_unformatted_matrix_file),intent(inout)::fileinfo
      class(mqc_basis_set),intent(out)::basisData
      character(len=*),intent(in),optional::filename
      character(len=256)::my_filename
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_L('Failed to retrieve geometry info from Gaussian matrix file: File not open.', 6, &
            'fileinfo%isOpen()', fileinfo%isOpen() )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'R' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      call fileInfo%getArray('SHELL TO ATOM MAP',basisData%sh2AtMp)
      call fileInfo%getArray('SHELL TYPES',basisData%shlTyp)
      call fileInfo%getArray('NUMBER OF PRIMITIVES PER SHELL',basisData%nPrmSh)
      call fileInfo%getArray('PRIMITIVE EXPONENTS',basisData%prmExp)
      call fileInfo%getArray('CONTRACTION COEFFICIENTS',basisData%conCoef)
      call fileInfo%getArray('P(S=P) CONTRACTION COEFFICIENTS',basisData%conCoTwo)
      call fileInfo%getArray('COORDINATES OF EACH SHELL',basisData%shCoor)
!
      return
      end Subroutine MQC_Gaussian_Unformatted_Matrix_Get_Basis_Set_Matrices
!
!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Molecule_Data
      Subroutine MQC_Gaussian_Unformatted_Matrix_Get_Molecule_Data(fileinfo,moleculeData,fileName)
!
!     This function is used to obtain the molecule information object
!     associated with the Gaussian unformatted matrix file sent in object
!     fileinfo.
!
!     L. M. Thompson, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(mqc_gaussian_unformatted_matrix_file),intent(inout)::fileinfo
      class(mqc_molecule_data),intent(out)::moleculeData
      character(len=*),intent(in),optional::filename
      character(len=256)::my_filename
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_L('Failed to retrieve geometry info from Gaussian matrix file: File not open.', 6, &
            'fileinfo%isOpen()', fileinfo%isOpen() )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'R' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      call mqc_gaussian_fill_molecule_data(moleculeData,fileInfo%nAtoms,fileInfo%atomicNumbers, &
        fileInfo%atomicWeights,fileInfo%atomicCharges,fileInfo%cartesians,fileInfo%iCharge, &
        fileInfo%multiplicity)
      return
      end Subroutine MQC_Gaussian_Unformatted_Matrix_Get_Molecule_Data


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Update_Header_Info
      Subroutine MQC_Gaussian_Unformatted_Matrix_Update_Header_Info(fileinfo,label,moleculeData,&
        integerIn,intVecIn,realVecIn)
!
!     This function is used to update data on an mqc_gaussian_unformatted matrix object. The 
!     recognized labels and their meaning include:
!           'molecule'          update the header with items stored in the molecule data object
!           'natoms'            update the number of atoms
!           'nbasis'            update the number of basis functions
!           'nbasis used'       update the number of used basis functions
!           'charge'            update the molecule charge
!           'multiplicity'      update the multiplicity
!           'nelectrons'        update the number of electrons
!           'icgu'              update the icgu flag (1XX = real, 2XX = complex, 1X = alpha/beta spin,
!                               2X = general spin, 1 = one spin block e.g. restricted or general, 
!                               2 = two spin blocks e.g. unrestricted)
!           'nfc'               update the NFC value
!           'nfv'               update the NFV value
!           'itran'             update the ITran value
!           'idum9'             update the IDum9 value
!           'nshlao'            update the NShlAO value
!           'nprmao'            update the NPrmAO value
!           'nshldb'            update the NShlDB value
!           'nprmdb'            update the NPrmDB value
!           'nbtot'             update the NBTot value
!           'atomic numbers'    update the atomic numbers vector
!           'atom types'        update the atom types vector
!           'basis2atom map'    update the vector that maps basis functions to atoms
!           'basis types'       update the vector of basis function types 
!           'gaussian scalars'  update the list of Gaussian scalars (containing e.g. energy etc.)
!           'atomic charges'    update the list of atomic charges
!           'atomic weights'    update the list of atomic weights
!           'cartesians'        update the list of cartesian coordinates
!
!     L. M. Thompson, 2020.
!
!
!     Variable Declarations.
!
      implicit none
      class(mqc_gaussian_unformatted_matrix_file),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      class(mqc_molecule_data),optional,intent(in)::moleculeData
      integer,optional,intent(in)::integerIn
      integer,dimension(:),allocatable,optional,intent(in)::intVecIn
      real,dimension(:),allocatable,optional,intent(in)::realVecIn
      character(len=64)::myLabel
      Integer::nInputArrays
!
      nInputArrays = 0
      if(present(moleculeData)) nInputArrays = nInputArrays + 1
      if(nInputArrays.ne.1) call mqc_error_i('Too many input values sent to mqc_update_header_info',&
        6,'nInputArrays',nInputArrays)
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('molecule')
        if(.not.present(moleculeData)) call mqc_error('Update molecule data requested in %updateHeader&
          & but input data is not molecule data type')
        call mqc_gaussian_output_molecule_data(moleculeData,fileInfo%nAtoms,fileInfo%atomicNumbers, &
          fileInfo%atomicWeights,fileInfo%atomicCharges,fileInfo%cartesians,fileInfo%iCharge, &
          fileInfo%multiplicity)
      case('natoms')
        if(.not.present(integerIn)) call mqc_error('Update atom number requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%nAtoms = integerIn
      case('nbasis')
        if(.not.present(integerIn)) call mqc_error('Update basis number requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%nBasis = integerIn
      case('nbasis used')
        if(.not.present(integerIn)) call mqc_error('Update used basis number requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%nBasisUse = integerIn
      case('charge')
        if(.not.present(integerIn)) call mqc_error('Update molecular charge requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%icharge = integerIn
      case('multiplicity')
        if(.not.present(integerIn)) call mqc_error('Update multiplicity requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%multiplicity = integerIn
      case('nelectrons')
        if(.not.present(integerIn)) call mqc_error('Update electron number requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%nElectrons = integerIn
      case('icgu')
        if(.not.present(integerIn)) call mqc_error('Update icgu flag requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%icgu = integerIn
      case('nfc')
        if(.not.present(integerIn)) call mqc_error('Update nfc parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%NFC = integerIn
      case('nfv')
        if(.not.present(integerIn)) call mqc_error('Update nfv parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%NFV = integerIn
      case('itran')
        if(.not.present(integerIn)) call mqc_error('Update itran parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%ITran = integerIn
      case('idum9')
        if(.not.present(integerIn)) call mqc_error('Update idum9 parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%IDum9 = integerIn
      case('nshlao')
        if(.not.present(integerIn)) call mqc_error('Update nshlao parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%NShlAO = integerIn
      case('nprmao')
        if(.not.present(integerIn)) call mqc_error('Update nprmao parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%NPrmAO = integerIn
      case('nshldb')
        if(.not.present(integerIn)) call mqc_error('Update nshldb parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%NShlDB = integerIn
      case('nprmdb')
        if(.not.present(integerIn)) call mqc_error('Update nprmdb parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%NPrmDB = integerIn
      case('nbtot')
        if(.not.present(integerIn)) call mqc_error('Update nbtot parameter requested in %updateHeader&
          & but input data is not integer type')
        fileInfo%NBTot = integerIn
      case('atomic numbers')
        if(.not.present(intVecIn)) call mqc_error('Update atomic numbers requested in %updateHeader&
          & but input data is not integer vector type')
        fileInfo%atomicNumbers = intVecIn
      case('atom types')
        if(.not.present(intVecIn)) call mqc_error('Update atom types requested in %updateHeader&
          & but input data is not integer vector type')
        fileInfo%atomTypes = intVecIn
      case('basis2atom map')
        if(.not.present(intVecIn)) call mqc_error('Update basis function to atom map requested in %updateHeader&
          & but input data is not integer vector type')
        fileInfo%basisFunction2Atom = intVecIn
      case('basis types')
        if(.not.present(intVecIn)) call mqc_error('Update basis types requested in %updateHeader&
          & but input data is not integer vector type')
        fileInfo%IbasisFunctionType = intVecIn
      case('gaussian scalars')
        if(.not.present(intVecIn)) call mqc_error('Update gaussian scalars requested in %updateHeader&
          & but input data is not integer vector type')
        fileInfo%IgaussianScalars = intVecIn
      case('atomic charges')
        if(.not.present(realVecIn)) call mqc_error('Update atomic charges requested in %updateHeader&
          & but input data is not real vector type')
        fileInfo%atomicCharges = realVecIn
      case('atomic weights')
        if(.not.present(realVecIn)) call mqc_error('Update atomic weights requested in %updateHeader&
          & but input data is not real vector type')
        fileInfo%atomicWeights = realVecIn
      case('cartesians')
        if(.not.present(realVecIn)) call mqc_error('Update Cartesian coordinates requested in %updateHeader&
          & but input data is not real vector type')
        fileInfo%cartesians = realVecIn
      case default
        call mqc_error_A('Invalid label sent to %updateHeader.',6,'mylabel',mylabel)
      end select

      return
      end Subroutine MQC_Gaussian_Unformatted_Matrix_Update_Header_Info


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Write_Basis_Set_Matrices
      Subroutine MQC_Gaussian_Unformatted_Matrix_Write_Basis_Set_Matrices(fileinfo,basisData,fileName)
!
!     This subroutine is used to write the arrays containing the basis set coefficients, exponents
!     and coordinates etc. from the basis set object to the Gaussian unformatted matrix file sent in 
!     object fileinfo.
!
!     L. M. Thompson, 2023.
!
!
!     Variable Declarations.
!
      implicit none
      class(mqc_gaussian_unformatted_matrix_file),intent(inout)::fileinfo
      class(mqc_basis_set),intent(out)::basisData
      character(len=*),intent(in),optional::filename
      character(len=256)::my_filename
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_L('Error writing basis set object to Gaussian matrix file', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'W' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_written) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      call fileInfo%writeArray('SHELL TO ATOM MAP',basisData%sh2AtMp)
      call fileInfo%writeArray('SHELL TYPES',basisData%shlTyp)
      call fileInfo%writeArray('NUMBER OF PRIMITIVES PER SHELL',basisData%nPrmSh)
      call fileInfo%writeArray('PRIMITIVE EXPONENTS',basisData%prmExp)
      call fileInfo%writeArray('CONTRACTION COEFFICIENTS',basisData%conCoef)
      call fileInfo%writeArray('P(S=P) CONTRACTION COEFFICIENTS',basisData%conCoTwo)
      call fileInfo%writeArray('COORDINATES OF EACH SHELL',basisData%shCoor)
!
      return
      end Subroutine MQC_Gaussian_Unformatted_Matrix_Write_Basis_Set_Matrices
!
!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Write_EST_Object
      subroutine mqc_gaussian_unformatted_matrix_write_EST_object(fileinfo,label, &
        est_wavefunction,est_integral,est_eigenvalues,filename,override,imagORide)
!
!     THIS SHOULD BE GAU_GET_EST_OBJ AND WE SHOULD HAVE A GENERAL ROUTINE IN EST OBJ
!     THAT CALLS THIS IF WE HAVE A GAUSSIAN FILE

!     This subroutine writes the desired MQC EST integral object as specified by
!     input argument <label> to a Gaussian unformatted matrix file sent in 
!     object <fileinfo>. The relevant information will be loaded from either 
!     (OPTIONAL) output dummy MQC_Wavefunction argument <est_wavefunction>, 
!     (OPTIONAL) output dummy MQC_SCF_Integral argument <est_integral>, 
!     (OPTIONAL) output dummy MQC_SCF_Eigenvalues argument <est_eigenvalues>.
!
!     Dummy argument <filename> is optional and is only used if fileinfo
!     hasn't already been defined using Routine
!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!     filename sent is different from the filename associated with object
!     fileinfo.
!
!     Dummy argument <override> is optional and can be used to write the EST
!     object as a particular wavefunction type (space, spin or general being 
!     the options). NOTE WE WILL WANT SOMETHING LIKE THIS FOR REAL AND COMPLEX 
!     WHEN MQC ALGEBRA HAS IT IMPLEMENTED.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!     called before calling this routine. The expectation is that
!     MQC_Gaussian_Unformatted_Matrix_Write_Header is also called before this
!     routine. However, it is also OK to call this routine first. In that case,
!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     The recognized labels and their meaning include:
!           'mo coefficients'    write the molecular orbital coefficients.
!           'mo energies'        write the molecular orbital energies.
!           'mo symmetries'      write the irreducible representation associated 
!                                  with each molecular orbital.*
!           'core hamiltonian'   write the core hamiltonian.
!           'fock'               write the fock matrix.
!           'density'            write the density matrix.
!           'scf density'        write the density matrix with the SCF label
!           'overlap'            write the overlap matrix.
!           'wavefunction'       export the wavefunction object.
!
!     * not yet implemented
!
!     Symmetric arrays are stored on the matrix file in the order (A(J,I),J=1,I),I=1,N)
!
!     L. M. Thompson, 2017.
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      type(mqc_wavefunction),optional::est_wavefunction
      type(mqc_scf_integral),optional::est_integral
      type(mqc_scf_eigenvalues),optional::est_eigenvalues
      character(len=*),intent(in),optional::filename,override,imagORide
      character(len=64)::myLabel
      character(len=256)::my_filename,my_override,my_integral_type,my_imagORide
      integer::nInputArrays,nBasis,nAlpha,nBeta
      type(mqc_matrix)::tmpMatrix
      type(mqc_vector)::tmpVector
      type(mqc_scalar)::tmpScalar
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_L('Error writing EST object to Gaussian matrix file', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'W' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_written) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Ensure that one and only one output MQC-type array has been sent from the
!     calling program unit.
!
      nInputArrays = 0
      if(Present(est_wavefunction)) nInputArrays = nInputArrays+1
      if(Present(est_integral)) nInputArrays = nInputArrays+1
      if(Present(est_eigenvalues)) nInputArrays = nInputArrays+1
      if(nInputArrays.ne.1) call mqc_error_i('Too many input arrays sent to Gaussian matrix file reading procedure.', 6, &
           'nInputArrays', nInputArrays )
!
!     Get the EST object in the desired wavefunction type format
!
      if(present(override)) then
        call String_Change_Case(override,'l',my_override)
        if(present(est_wavefunction)) call mqc_error_l('Overriding wavefunction types not implemented', 6, & 
             'present(est_wavefunction)', present(est_wavefunction) )
        if(my_override.eq.'space') then
          my_integral_type = 'space'
        elseIf(my_override.eq.'spin') then
          my_integral_type = 'spin'
        elseIf(my_override.eq.'general') then
          my_integral_type = 'general'
        else
          call mqc_error_A('Unrecognised override type in %writeESTObj', 6, &
               'my_override', my_override )
        endIf
      else
        if(present(est_integral)) then
          my_integral_type = mqc_integral_array_type(est_integral)
        elseIf(present(est_eigenvalues)) then
          my_integral_type = mqc_eigenvalues_array_type(est_eigenvalues)
        endIf
      endIf

      if(present(imagORide)) then
        call String_Change_Case(imagORide,'l',my_imagORide)
        if(present(est_wavefunction)) call mqc_error_l('Overriding wavefunction types not implemented', 6, &
          'present(est_wavefunction)', present(est_wavefunction) )
        if(present(est_eigenvalues)) call mqc_error_l('Overriding eigenvalue real/complex not implemented', 6, &
          'present(est_eigenvalues)', present(est_eigenvalues) )
        if(my_imagORide.eq.'real') then
          if(present(est_integral)) then
            est_integral = est_integral%orbitals('real') 
          endIf
        elseif(my_imagORide.eq.'complex') then
          if(present(est_integral)) then
            est_integral = est_integral%orbitals('complex') 
          endIf
        endIf
      endIf

!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('mo coefficients')
        if(.not.(Present(est_integral))) call mqc_error_l('wrong EST type in writeESTOBj', 6, &
             'Present(est_integral)', Present(est_integral) )
        if(my_integral_type.eq.'space') then
          call fileInfo%writeArray('ALPHA MO COEFFICIENTS', &
            matrixIn=est_integral%getBlock('alpha'),storage='full')
        elseIf(my_integral_type.eq.'spin') then
          call fileInfo%writeArray('ALPHA MO COEFFICIENTS', &
            matrixIn=est_integral%getBlock('alpha'),storage='full')
          call fileInfo%writeArray('BETA MO COEFFICIENTS', &
            matrixIn=est_integral%getBlock('beta'),storage='full')
        elseIf(my_integral_type.eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_integral,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('ALPHA MO COEFFICIENTS',matrixIn=tmpMatrix,storage='full')
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'my_integral_type', my_integral_type )
        endIf
      case('mo energies')
        if(.not.(Present(est_eigenvalues))) call mqc_error_l('wrong EST type in writeESTOBj', 6, &
             'Present(est_eigenvalues)', Present(est_eigenvalues) )
        if(my_integral_type.eq.'space') then
          call fileInfo%writeArray('ALPHA ORBITAL ENERGIES', &
            vectorIn=est_eigenvalues%getBlock('alpha'))
        elseIf(my_integral_type.eq.'spin') then
          call fileInfo%writeArray('ALPHA ORBITAL ENERGIES', &
            vectorIn=est_eigenvalues%getBlock('alpha'))
          call fileInfo%writeArray('BETA ORBITAL ENERGIES', &
            vectorIn=est_eigenvalues%getBlock('beta'))
        elseIf(my_integral_type.eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_eigenvalues,tmpVector)
          call fileInfo%writeArray('ALPHA ORBITAL ENERGIES',vectorIn=tmpVector)
        else
          call mqc_error_a('Unknown wavefunction type in getESTObj', 6, &
               'my_integral_type', my_integral_type )
        endIf
      case('core hamiltonian')
        if(.not.(Present(est_integral))) call mqc_error_l('wrong EST type in writeESTOBj', 6, &
             'Present(est_integral)', Present(est_integral) )
        if(my_integral_type.eq.'space') then
          call fileInfo%writeArray('CORE HAMILTONIAN ALPHA', &
            matrixIn=est_integral%getBlock('alpha'))
        elseIf(my_integral_type.eq.'spin') then
          call fileInfo%writeArray('CORE HAMILTONIAN ALPHA', &
            matrixIn=est_integral%getBlock('alpha'))
          call fileInfo%writeArray('CORE HAMILTONIAN BETA', &
            matrixIn=est_integral%getBlock('beta'))
        elseIf(my_integral_type.eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_integral,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('CORE HAMILTONIAN ALPHA',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'my_integral_type', my_integral_type )
        endIf
      case('fock')
        if(.not.(Present(est_integral))) call mqc_error_l('wrong EST type in writeESTOBj', 6, &
             'Present(est_integral)', Present(est_integral) )
        if(my_integral_type.eq.'space') then
          call fileInfo%writeArray('ALPHA FOCK MATRIX', &
            matrixIn=est_integral%getBlock('alpha'))
        elseIf(my_integral_type.eq.'spin') then
          call fileInfo%writeArray('ALPHA FOCK MATRIX', &
            matrixIn=est_integral%getBlock('alpha'))
          call fileInfo%writeArray('BETA FOCK MATRIX', &
            matrixIn=est_integral%getBlock('beta'))
        elseIf(my_integral_type.eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_integral,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('ALPHA FOCK MATRIX',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'my_integral_type', my_integral_type )
        endIf
      case('scf density')
        if(.not.(Present(est_integral))) call mqc_error_L('wrong EST type in writeESTOBj', 6, &
             'Present(est_integral)', Present(est_integral) )
        if(my_integral_type.eq.'space') then
          call fileInfo%writeArray('ALPHA SCF DENSITY MATRIX', &
            matrixIn=est_integral%getBlock('alpha'))
        elseIf(my_integral_type.eq.'spin') then
          call fileInfo%writeArray('ALPHA SCF DENSITY MATRIX', &
            matrixIn=est_integral%getBlock('alpha'))
          call fileInfo%writeArray('BETA SCF DENSITY MATRIX', &
            matrixIn=est_integral%getBlock('beta'))
        elseIf(my_integral_type.eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_integral,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('ALPHA SCF DENSITY MATRIX',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'my_integral_type', my_integral_type )
        endIf
      case('density')
         if(.not.(Present(est_integral))) call mqc_error_L('wrong EST type in writeESTOBj', 6, &
              'Present(est_integral)', Present(est_integral) )
         if(my_integral_type.eq.'space') then
           call fileInfo%writeArray('ALPHA DENSITY MATRIX', &
             matrixIn=est_integral%getBlock('alpha'))
         elseIf(my_integral_type.eq.'spin') then
           call fileInfo%writeArray('ALPHA DENSITY MATRIX', &
             matrixIn=est_integral%getBlock('alpha'))
           call fileInfo%writeArray('BETA DENSITY MATRIX', &
             matrixIn=est_integral%getBlock('beta'))
         elseIf(my_integral_type.eq.'general') then
           call mqc_matrix_undoSpinBlockGHF(est_integral,tmpMatrix)
           if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix)
           call fileInfo%writeArray('ALPHA DENSITY MATRIX',matrixIn=tmpMatrix)
         else
           call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
                'my_integral_type', my_integral_type )
         endIf
      case('overlap')
        if(.not.(Present(est_integral))) call mqc_error_L('wrong EST type in writeESTOBj', 6, &
             'Present(est_integral)', Present(est_integral) )
        if(my_integral_type.eq.'space') then
          call fileInfo%writeArray('OVERLAP', &
            matrixIn=est_integral%getBlock('alpha'))
        elseIf(my_integral_type.eq.'spin') then
          call fileInfo%writeArray('OVERLAP', &
            matrixIn=est_integral%getBlock('alpha'))
        elseIf(my_integral_type.eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_integral,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('OVERLAP',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'my_integral_type', my_integral_type )
        endIf
      case('wavefunction')
        if(.not.(Present(est_wavefunction))) call mqc_error_l('wrong EST type in writeESTOBj', 6, &
             'Present(est_wavefunction)', Present(est_wavefunction) )
        if(mqc_integral_array_type(est_wavefunction%overlap_matrix).eq.'space') then
          call fileInfo%writeArray('OVERLAP', &
            matrixIn=est_wavefunction%overlap_matrix%getBlock('alpha'))
        elseIf(mqc_integral_array_type(est_wavefunction%overlap_matrix).eq.'spin') then
          call fileInfo%writeArray('OVERLAP', &
            matrixIn=est_wavefunction%overlap_matrix%getBlock('alpha'))
        elseIf(mqc_integral_array_type(est_wavefunction%overlap_matrix).eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_wavefunction%overlap_matrix,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('OVERLAP',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'mqc_integral_array_type(est_wavefunction%overlap_matrix)', &
               mqc_integral_array_type(est_wavefunction%overlap_matrix) )
        endIf
        if(mqc_integral_array_type(est_wavefunction%core_hamiltonian).eq.'space') then
          call fileInfo%writeArray('CORE HAMILTONIAN ALPHA', &
            matrixIn=est_wavefunction%core_hamiltonian%getBlock('alpha'))
        elseIf(mqc_integral_array_type(est_wavefunction%core_hamiltonian).eq.'spin') then
          call fileInfo%writeArray('CORE HAMILTONIAN ALPHA', &
            matrixIn=est_wavefunction%core_hamiltonian%getBlock('alpha'))
          call fileInfo%writeArray('CORE HAMILTONIAN BETA', &
            matrixIn=est_wavefunction%core_hamiltonian%getBlock('beta'))
        elseIf(mqc_integral_array_type(est_wavefunction%core_hamiltonian).eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_wavefunction%core_hamiltonian,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('CORE HAMILTONIAN ALPHA',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'mqc_integral_array_type(est_wavefunction%core_hamiltonian)', &
               mqc_integral_array_type(est_wavefunction%core_hamiltonian) )
        endIf
        if(mqc_eigenvalues_array_type(est_wavefunction%mo_energies).eq.'space') then
          call fileInfo%writeArray('ALPHA ORBITAL ENERGIES', &
            vectorIn=est_wavefunction%mo_energies%getBlock('alpha'))
        elseIf(mqc_eigenvalues_array_type(est_wavefunction%mo_energies).eq.'spin') then
          call fileInfo%writeArray('ALPHA ORBITAL ENERGIES', &
            vectorIn=est_wavefunction%mo_energies%getBlock('alpha'))
          call fileInfo%writeArray('BETA ORBITAL ENERGIES', &
            vectorIn=est_wavefunction%mo_energies%getBlock('beta'))
        elseIf(mqc_eigenvalues_array_type(est_wavefunction%mo_energies).eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_wavefunction%mo_energies,tmpVector)
          call fileInfo%writeArray('ALPHA ORBITAL ENERGIES',vectorIn=tmpVector)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'mqc_eigenvalues_array_type(est_wavefunction%mo_energies)', &
               mqc_eigenvalues_array_type(est_wavefunction%mo_energies) )
        endIf
        if(mqc_integral_array_type(est_wavefunction%mo_coefficients).eq.'space') then
          call fileInfo%writeArray('ALPHA MO COEFFICIENTS', &
            matrixIn=est_wavefunction%mo_coefficients%getBlock('alpha'))
        elseIf(mqc_integral_array_type(est_wavefunction%mo_coefficients).eq.'spin') then
          call fileInfo%writeArray('ALPHA MO COEFFICIENTS', &
            matrixIn=est_wavefunction%mo_coefficients%getBlock('alpha'))
          call fileInfo%writeArray('BETA MO COEFFICIENTS', &
            matrixIn=est_wavefunction%mo_coefficients%getBlock('beta'))
        elseIf(mqc_integral_array_type(est_wavefunction%mo_coefficients).eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_wavefunction%mo_coefficients,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('ALPHA MO COEFFICIENTS',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
          'mqc_integral_array_type(est_wavefunction%mo_coefficients)', &
          mqc_integral_array_type(est_wavefunction%mo_coefficients) )
        endIf
        if(mqc_integral_array_type(est_wavefunction%density_matrix).eq.'space') then
          call fileInfo%writeArray('ALPHA DENSITY MATRIX', &
            matrixIn=est_wavefunction%density_matrix%getBlock('alpha'))
        elseIf(mqc_integral_array_type(est_wavefunction%density_matrix).eq.'spin') then
          call fileInfo%writeArray('ALPHA DENSITY MATRIX', &
            matrixIn=est_wavefunction%density_matrix%getBlock('alpha'))
          call fileInfo%writeArray('BETA DENSITY MATRIX', &
            matrixIn=est_wavefunction%density_matrix%getBlock('beta'))
        elseIf(mqc_integral_array_type(est_wavefunction%density_matrix).eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_wavefunction%density_matrix,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('ALPHA DENSITY MATRIX',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'mqc_integral_array_type(est_wavefunction%density_matrix)', &
               mqc_integral_array_type(est_wavefunction%density_matrix) )
        endIf
        if(mqc_integral_array_type(est_wavefunction%scf_density_matrix).eq.'space') then
          call fileInfo%writeArray('ALPHA SCF DENSITY MATRIX', &
            matrixIn=est_wavefunction%scf_density_matrix%getBlock('alpha'))
        elseIf(mqc_integral_array_type(est_wavefunction%scf_density_matrix).eq.'spin') then
          call fileInfo%writeArray('ALPHA SCF DENSITY MATRIX', &
            matrixIn=est_wavefunction%scf_density_matrix%getBlock('alpha'))
          call fileInfo%writeArray('BETA SCF DENSITY MATRIX', &
            matrixIn=est_wavefunction%scf_density_matrix%getBlock('beta'))
        elseIf(mqc_integral_array_type(est_wavefunction%scf_density_matrix).eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_wavefunction%scf_density_matrix,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('ALPHA SCF DENSITY MATRIX',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'mqc_integral_array_type(est_wavefunction%scf_density_matrix)', &
               mqc_integral_array_type(est_wavefunction%scf_density_matrix) )
        endIf
        if(mqc_integral_array_type(est_wavefunction%fock_matrix).eq.'space') then
          call fileInfo%writeArray('ALPHA FOCK MATRIX', &
            matrixIn=est_wavefunction%fock_matrix%getBlock('alpha'))
        elseIf(mqc_integral_array_type(est_wavefunction%fock_matrix).eq.'spin') then
          call fileInfo%writeArray('ALPHA FOCK MATRIX', &
            matrixIn=est_wavefunction%fock_matrix%getBlock('alpha'))
          call fileInfo%writeArray('BETA FOCK MATRIX', &
            matrixIn=est_wavefunction%fock_matrix%getBlock('beta'))
        elseIf(mqc_integral_array_type(est_wavefunction%fock_matrix).eq.'general') then
          call mqc_matrix_undoSpinBlockGHF(est_wavefunction%fock_matrix,tmpMatrix)
          if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
          call fileInfo%writeArray('ALPHA FOCK MATRIX',matrixIn=tmpMatrix)
        else
          call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
               'mqc_integral_array_type(est_wavefunction%fock_matrix)', &
               mqc_integral_array_type(est_wavefunction%fock_matrix) )
        endIf
      case default
        call String_Change_Case(label,'u',myLabel)
        If(present(est_eigenvalues)) then
          if(my_integral_type.eq.'space') then
            call fileInfo%writeArray('ALPHA '//trim(myLabel), &
              vectorIn=est_eigenvalues%getBlock('alpha'))
          elseIf(my_integral_type.eq.'spin') then
            call fileInfo%writeArray('ALPHA '//trim(myLabel), &
              vectorIn=est_eigenvalues%getBlock('alpha'))
            call fileInfo%writeArray('BETA '//trim(myLabel), &
              vectorIn=est_eigenvalues%getBlock('beta'))
          elseIf(my_integral_type.eq.'general') then
            call mqc_matrix_undoSpinBlockGHF(est_eigenvalues,tmpVector)
            call fileInfo%writeArray('ALPHA '//trim(myLabel),vectorIn=tmpVector)
          else
            call mqc_error_a('Unknown wavefunction type in getESTObj', 6, &
                 'my_integral_type', my_integral_type )
          endIf
        elseIf(present(est_integral)) then 
          if(my_integral_type.eq.'space') then
            call fileInfo%writeArray('ALPHA '//trim(myLabel), &
              matrixIn=est_integral%getBlock('alpha'),storage='full')
          elseIf(my_integral_type.eq.'spin') then
            call fileInfo%writeArray('ALPHA '//trim(myLabel), &
              matrixIn=est_integral%getBlock('alpha'),storage='full')
            call fileInfo%writeArray('BETA '//trim(myLabel), &
              matrixIn=est_integral%getBlock('beta'),storage='full')
          elseIf(my_integral_type.eq.'general') then
            call mqc_matrix_undoSpinBlockGHF(est_integral,tmpMatrix)
            if(.not.mqc_matrix_haveComplex(tmpMatrix)) call MQC_Matrix_Copy_Real2Complex(tmpMatrix) 
            call fileInfo%writeArray('ALPHA '//trim(myLabel),matrixIn=tmpMatrix,storage='full')
          else
            call mqc_error_a('Unknown wavefunction type in writeESTObj', 6, &
                 'my_integral_type', my_integral_type )
          endIf
        else
          call mqc_error('Unsuported open label EST object in writeESTObj')
        endIf
      end select
!
      return
!
      end subroutine MQC_Gaussian_Unformatted_Matrix_Write_EST_Object 
!
!
!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_EST_Object
      subroutine mqc_gaussian_unformatted_matrix_get_EST_object(fileinfo,label, &
        est_wavefunction,est_integral,est_eigenvalues,filename,foundObj)
!
!     IS IT POSSIBLE TO MAKE THIS GAU_GET_EST_OBJ AND MAKE A GENERAL ROUTINE IN 
!     EST OBJ THAT CALLS THIS IF WE HAVE A GAUSSIAN FILE? AS FAR AS I CAN TELL
!     WE CAN'T AS THIS REQUIRES MQC_EST TO CALL HIGHER LEVEL MODULES.
!
!     This subroutine loads the desired MQC EST integral object as specified by
!     input argument <label> from a Gaussian unformatted matrix file sent in 
!     object <fileinfo>. The relevant information will be loaded into either 
!     (OPTIONAL) output dummy MQC_Wavefunction argument <est_wavefunction>, 
!     (OPTIONAL) output dummy MQC_SCF_Integral argument <est_integral>, 
!     (OPTIONAL) output dummy MQC_SCF_Eigenvalues argument <est_eigenvalues>, or 
!     (OPTIONAL) output dummy MQC_Scalar argument <scalarOut>, or 
!     (OPTIONAL) output dummy character argument <characterOut>, or 
!     (OPTIONAL) output dummy logical argument <logicalOut>.
!
!     Dummy argument <filename> is optional and is only used if fileinfo
!     hasn't already been defined using Routine
!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!     filename sent is different from the filename associated with object
!     fileinfo.
!
!     If foundObj is present it is returned true if the EST object is sucessfully
!     loaded and false if the EST object is not. If foundObj is not present then
!     an error message is called when the EST object cannot be loaded. The 
!     exception is when the 'wavefunction' object is called, where foundObj is
!     returned true only if all ESTobjects are loaded, false if any EST object
!     cannot be loaded, and does not call an error if foundObj is not present.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!     called before calling this routine. The expectation is that
!     MQC_Gaussian_Unformatted_Matrix_Read_Header is also called before this
!     routine. However, it is also OK to call this routine first. In that case,
!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     The recognized labels and their meaning include:
!           'mo coefficients'    return the molecular orbital coefficients.
!           'mo energies'        return the molecular orbital energies.
!           'mo symmetries'      return the irreducible representation associated 
!                                with each molecular orbital.*
!           'core hamiltonian'   return the core hamiltonian.
!           'fock'               return the fock matrix.
!           'density'            return the density matrix.
!           'scf density'        return the SCF density matrix.
!           'overlap'            return the overlap matrix.
!           'dipole x'           return the length gauge x dipole AO integrals. 
!           'dipole y'           return the length gauge y dipole AO integrals. 
!           'dipole z'           return the length gauge z dipole AO integrals.
!           'vel dipole x'       return the velocity gauge x dipole AO integrals. 
!           'vel dipole y'       return the velocity gauge y dipole AO integrals. 
!           'vel dipole z'       return the velocity gauge z dipole AO integrals.
!           'wavefunction'       load the wavefunction object (does not contain 
!                                dipoles).
!
!     * not yet implemented
!
!     Symmetric arrays are stored on the matrix file in the order (A(J,I),J=1,I),I=1,N)
!
!     L. M. Thompson, 2017.
!
!     GHF routine has been updated to use modified mqc_matrix_spinBlockGHF
!     subroutine.  nAlpha electrons are passed as optional second dummy argument.
!     
!     -A. Mahler, 4/26/18
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      class(mqc_wavefunction),optional::est_wavefunction
      type(mqc_scf_integral),optional::est_integral
      type(mqc_scf_eigenvalues),optional::est_eigenvalues
      character(len=*),intent(in),optional::filename
      logical,optional::foundObj
      character(len=64)::myLabel
      character(len=256)::my_filename,errorMsg
      integer::nOutputArrays,nBasis,nElectrons,multiplicity
      integer(kind=int64),dimension(:),allocatable::elist
      type(mqc_matrix)::tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha
      type(mqc_vector)::tmpVectorAlpha,tmpVectorBeta
      type(mqc_scalar)::tmpScalar
      logical::found,OK
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_L('Error reading Gaussian matrix file header (6): Must include a filename.', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'R' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Ensure that one and only one output MQC-type array has been sent from the
!     calling program unit.
!
      nOutputArrays = 0
      if(Present(est_wavefunction)) nOutputArrays = nOutputArrays+1
      if(Present(est_integral)) nOutputArrays = nOutputArrays+1
      if(Present(est_eigenvalues)) nOutputArrays = nOutputArrays+1
      if(nOutputArrays.ne.1) call mqc_error_i('Too many output arrays sent to Gaussian matrix file reading procedure.', 6, &
           'nOutputArrays', nOutputArrays )
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('mo coefficients')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('ALPHA MO COEFFICIENTS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA MO COEFFICIENTS not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call mqc_integral_allocate(est_integral,'mo coefficients','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('ALPHA MO COEFFICIENTS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA MO COEFFICIENTS not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call fileInfo%getArray('BETA MO COEFFICIENTS',tmpMatrixBeta,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'BETA MO COEFFICIENTS not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
              call mqc_integral_allocate(est_integral,'mo coefficients','spin',tmpMatrixAlpha, &
                tmpMatrixBeta)
            endIf
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('ALPHA MO COEFFICIENTS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA MO COEFFICIENTS not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha,fileInfo%getVal('nElectrons'), &
              fileInfo%getVal('multiplicity'),elist)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'mo coefficients','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
            call est_integral%setEList(elist)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('mo energies')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('ALPHA ORBITAL ENERGIES',vectorOut=tmpVectorAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA ORBITAL ENERGIES not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call mqc_eigenvalues_allocate(est_eigenvalues,'mo energies','space',tmpVectorAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('ALPHA ORBITAL ENERGIES',vectorOut=tmpVectorAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA ORBITAL ENERGIES not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call fileInfo%getArray('BETA ORBITAL ENERGIES',vectorOut=tmpVectorBeta,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'BETA ORBITAL ENERGIES not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
              call mqc_eigenvalues_allocate(est_eigenvalues,'mo energies','spin',tmpVectorAlpha, &
                tmpVectorBeta)
            endIf
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('ALPHA ORBITAL ENERGIES',vectorOut=tmpVectorAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA ORBITAL ENERGIES not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpVectorAlpha)
            tmpVectorBeta = tmpVectorAlpha%vat(nBasis+1,-1)
            tmpVectorAlpha = tmpVectorAlpha%vat(1,nBasis)
            call mqc_eigenvalues_allocate(est_eigenvalues,'mo energies','general',tmpVectorAlpha, &
              tmpVectorBeta)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('core hamiltonian')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('CORE HAMILTONIAN ALPHA',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'CORE HAMILTONIAN ALPHA not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'core hamiltonian','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('CORE HAMILTONIAN ALPHA',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'CORE HAMILTONIAN ALPHA not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call fileInfo%getArray('CORE HAMILTONIAN BETA',tmpMatrixBeta,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'CORE HAMILTONIAN BETA not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
!              if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!                call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!                tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!              endIf
!              if(MQC_Matrix_HaveComplex(tmpMatrixBeta)) then
!                call mqc_matrix_symm2full(tmpMatrixBeta,'hermitian')
!                tmpMatrixBeta = transpose(tmpMatrixBeta)
!              endIf
              call mqc_integral_allocate(est_integral,'core hamiltonian','spin',tmpMatrixAlpha, &
                tmpMatrixBeta)
            endIf
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('CORE HAMILTONIAN ALPHA',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'CORE HAMILTONIAN ALPHA not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'core hamiltonian','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('fock')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('ALPHA FOCK MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA FOCK MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'fock','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('ALPHA FOCK MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA FOCK MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call fileInfo%getArray('BETA FOCK MATRIX',tmpMatrixBeta,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'BETA FOCK MATRIX not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
!              if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!                call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!                tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!              endIf
!              if(MQC_Matrix_HaveComplex(tmpMatrixBeta)) then
!                call mqc_matrix_symm2full(tmpMatrixBeta,'hermitian')
!                tmpMatrixBeta = transpose(tmpMatrixBeta)
!              endIf
              call mqc_integral_allocate(est_integral,'fock','spin',tmpMatrixAlpha, &
                tmpMatrixBeta)
            endIf
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('ALPHA FOCK MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA FOCK MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'fock','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('density')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('ALPHA DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA DENSITY MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'density','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('ALPHA DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA DENSITY MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call fileInfo%getArray('BETA DENSITY MATRIX',tmpMatrixBeta,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'BETA DENSITY MATRIX not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
!              if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!                call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!                tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!              endIf
!              if(MQC_Matrix_HaveComplex(tmpMatrixBeta)) then
!                call mqc_matrix_symm2full(tmpMatrixBeta,'hermitian')
!                tmpMatrixBeta = transpose(tmpMatrixBeta)
!              endIf
              call mqc_integral_allocate(est_integral,'density','spin',tmpMatrixAlpha, &
                tmpMatrixBeta)
            endIf
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('ALPHA DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA DENSITY MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'density','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('scf density')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('ALPHA SCF DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA SCF DENSITY MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'density','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('ALPHA SCF DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA SCF DENSITY MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call fileInfo%getArray('BETA SCF DENSITY MATRIX',tmpMatrixBeta,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'BETA SCF DENSITY MATRIX not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
!              if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!                call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!                tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!              endIf
!              if(MQC_Matrix_HaveComplex(tmpMatrixBeta)) then
!                call mqc_matrix_symm2full(tmpMatrixBeta,'hermitian')
!                tmpMatrixBeta = transpose(tmpMatrixBeta)
!              endIf
              call mqc_integral_allocate(est_integral,'density','spin',tmpMatrixAlpha, &
                tmpMatrixBeta)
            endIf
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('ALPHA SCF DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'ALPHA SCF DENSITY MATRIX not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'density','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('overlap')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('OVERLAP',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'OVERLAP not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'overlap','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('OVERLAP',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'OVERLAP not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'overlap','spin',tmpMatrixAlpha, &
              tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('OVERLAP',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'OVERLAP not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'overlap','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('dipole x')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (x) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'dipole x','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (x) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'dipole x','spin',tmpMatrixAlpha, &
              tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (x) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'dipole x','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('dipole y')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,arrayNum=2,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (y) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'dipole y','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,arrayNum=2,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (y) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'dipole y','spin',tmpMatrixAlpha, &
              tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,arrayNum=2,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (y) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'dipole y','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('dipole z')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,arrayNum=3,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (z) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'dipole z','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,arrayNum=3,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (z) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'dipole z','spin',tmpMatrixAlpha, &
              tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('DIPOLE INTEGRALS',tmpMatrixAlpha,arrayNum=3,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIPOLE INTEGRALS (z) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'dipole y','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('vel dipole x')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (x) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'vel dipole x','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (x) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'vel dipole x','spin',tmpMatrixAlpha, &
              tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (x) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'vel dipole x','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('vel dipole y')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,arrayNum=2,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (y) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'vel dipole y','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,arrayNum=2,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (y) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'vel dipole y','spin',tmpMatrixAlpha, &
              tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,arrayNum=2,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (y) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'vel dipole y','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('vel dipole z')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,arrayNum=3,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (z) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'vel dipole z','space',tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,arrayNum=3,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (z) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_integral,'vel dipole z','spin',tmpMatrixAlpha, &
              tmpMatrixAlpha)
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('DIP VEL INTEGRALS',tmpMatrixAlpha,arrayNum=3,foundOut=found)
          if(present(foundObj)) foundObj = found
          if(.not.found) then
            errorMsg = 'DIP VEL INTEGRALS (z) not present on file'
            if(present(foundObj)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            nBasis = fileInfo%getVal('nBasis')
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_integral,'vel dipole y','general',tmpMatrixAlpha, &
              tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case('wavefunction')
        if(present(foundObj)) foundObj = .true.
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('OVERLAP',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_wavefunction%overlap_matrix,'overlap','space', &
              tmpMatrixAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'OVERLAP not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          call fileInfo%getArray('CORE HAMILTONIAN ALPHA',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_wavefunction%core_hamiltonian,'core hamiltonian','space', &
              tmpMatrixAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'CORE HAMILTONIAN ALPHA not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          call fileInfo%getArray('ALPHA ORBITAL ENERGIES',vectorOut=tmpVectorAlpha,foundOut=found)
          if(found) then
            call mqc_eigenvalues_allocate(est_wavefunction%mo_energies,'mo energies','space', &
              tmpVectorAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA ORBITAL ENERGIES not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          call fileInfo%getArray('ALPHA MO COEFFICIENTS',tmpMatrixAlpha,foundOut=found)
          if(found) then
            call mqc_integral_allocate(est_wavefunction%mo_coefficients,'mo coefficients','space', &
              tmpMatrixAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA MO COEFFICIENTS not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          call fileInfo%getArray('ALPHA DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_wavefunction%density_matrix,'density','space', &
              tmpMatrixAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA DENSITY MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          call fileInfo%getArray('ALPHA SCF DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_wavefunction%scf_density_matrix,'density','space', &
              tmpMatrixAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA SCF DENSITY MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          call fileInfo%getArray('ALPHA FOCK MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_wavefunction%fock_matrix,'fock','space',tmpMatrixAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA FOCK MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          est_wavefunction%nBasis = fileInfo%getVal('nBasis')
          est_wavefunction%nAlpha = fileInfo%getVal('nAlpha')
          est_wavefunction%nBeta = fileInfo%getVal('nBeta')
          est_wavefunction%nElectrons = fileInfo%getVal('nElectrons')
          est_wavefunction%charge = fileInfo%getVal('charge')
          est_wavefunction%multiplicity = fileInfo%getVal('multiplicity')
          call mqc_gaussian_ICGU(fileInfo%ICGU,est_wavefunction%wf_type,est_wavefunction%wf_complex)
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('OVERLAP',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_integral_allocate(est_wavefunction%overlap_matrix,'overlap','spin', &
              tmpMatrixAlpha,tmpMatrixAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'OVERLAP not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          call fileInfo%getArray('CORE HAMILTONIAN ALPHA',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call fileInfo%getArray('CORE HAMILTONIAN BETA',tmpMatrixBeta,foundOut=found)
            if(found) then
!              if(MQC_Matrix_HaveComplex(tmpMatrixBeta)) then
!                call mqc_matrix_symm2full(tmpMatrixBeta,'hermitian')
!                tmpMatrixBeta = transpose(tmpMatrixBeta)
!              endIf
              call mqc_integral_allocate(est_wavefunction%core_hamiltonian,'core hamiltonian','spin', &
                tmpMatrixAlpha,tmpMatrixBeta)
            else
              if(present(foundObj)) foundObj = .false.
              write(6,'(A)') 'CORE HAMILTONIAN BETA not present on file - skipping'
              my_filename = TRIM(fileinfo%filename)
              call fileinfo%CLOSEFILE()
              call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
            endIf
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'CORE HAMILTONIAN ALPHA not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf
          call fileInfo%getArray('ALPHA ORBITAL ENERGIES',vectorOut=tmpVectorAlpha,foundOut=found)
          if(found) then
            call fileInfo%getArray('BETA ORBITAL ENERGIES',vectorOut=tmpVectorBeta,foundOut=found)
            if(found) then
              call mqc_eigenvalues_allocate(est_wavefunction%mo_energies,'mo energies','spin', &
                tmpVectorAlpha,tmpVectorBeta)
            else
              if(present(foundObj)) foundObj = .false.
              write(6,'(A)') 'BETA ORBITAL ENERGIES not present on file - skipping'
              my_filename = TRIM(fileinfo%filename)
              call fileinfo%CLOSEFILE()
              call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
            endIf
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA ORBITAL ENERGIES not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA MO COEFFICIENTS',tmpMatrixAlpha,foundOut=found)
          if(found) then
            call fileInfo%getArray('BETA MO COEFFICIENTS',tmpMatrixBeta,foundOut=found)
            if(found) then
              call mqc_integral_allocate(est_wavefunction%mo_coefficients,'mo coefficients','spin', &
                tmpMatrixAlpha,tmpMatrixBeta)
            else
              if(present(foundObj)) foundObj = .false.
              write(6,'(A)') 'BETA MO COEFFICIENTS not present on file - skipping'
              my_filename = TRIM(fileinfo%filename)
              call fileinfo%CLOSEFILE()
              call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
            endIf
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA MO COEFFICIENTS not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call fileInfo%getArray('BETA DENSITY MATRIX',tmpMatrixBeta,foundOut=found)
            if(found) then
!              if(MQC_Matrix_HaveComplex(tmpMatrixBeta)) then
!                call mqc_matrix_symm2full(tmpMatrixBeta,'hermitian')
!                tmpMatrixBeta = transpose(tmpMatrixBeta)
!              endIf
              call mqc_integral_allocate(est_wavefunction%density_matrix,'density','spin', &
                tmpMatrixAlpha,tmpMatrixBeta)
            else
              if(present(foundObj)) foundObj = .false.
              write(6,'(A)') 'BETA DENSITY MATRIX not present on file - skipping'
              my_filename = TRIM(fileinfo%filename)
              call fileinfo%CLOSEFILE()
              call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
            endIf
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA DENSITY MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA SCF DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call fileInfo%getArray('BETA SCF DENSITY MATRIX',tmpMatrixBeta,foundOut=found)
            if(found) then
!              if(MQC_Matrix_HaveComplex(tmpMatrixBeta)) then
!                call mqc_matrix_symm2full(tmpMatrixBeta,'hermitian')
!                tmpMatrixBeta = transpose(tmpMatrixBeta)
!              endIf
              call mqc_integral_allocate(est_wavefunction%scf_density_matrix,'density','spin', &
                tmpMatrixAlpha,tmpMatrixBeta)
            else
              if(present(foundObj)) foundObj = .false.
              write(6,'(A)') 'BETA SCF DENSITY MATRIX not present on file - skipping'
              my_filename = TRIM(fileinfo%filename)
              call fileinfo%CLOSEFILE()
              call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
            endIf
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA SCF DENSITY MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA FOCK MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call fileInfo%getArray('BETA FOCK MATRIX',tmpMatrixBeta,foundOut=found)
            if(found) then
!              if(MQC_Matrix_HaveComplex(tmpMatrixBeta)) then
!                call mqc_matrix_symm2full(tmpMatrixBeta,'hermitian')
!                tmpMatrixBeta = transpose(tmpMatrixBeta)
!              endIf
              call mqc_integral_allocate(est_wavefunction%fock_matrix,'fock','spin',tmpMatrixAlpha, &
                tmpMatrixBeta)
            else
              if(present(foundObj)) foundObj = .false.
              write(6,'(A)') 'BETA FOCK MATRIX not present on file - skipping'
              my_filename = TRIM(fileinfo%filename)
              call fileinfo%CLOSEFILE()
              call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
              endIf
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA FOCK MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          est_wavefunction%nBasis = fileInfo%getVal('nBasis')
          est_wavefunction%nAlpha = fileInfo%getVal('nAlpha')
          est_wavefunction%nBeta = fileInfo%getVal('nBeta')
          est_wavefunction%nElectrons = fileInfo%getVal('nElectrons')
          est_wavefunction%charge = fileInfo%getVal('charge')
          est_wavefunction%multiplicity = fileInfo%getVal('multiplicity')
          call mqc_gaussian_ICGU(fileInfo%ICGU,est_wavefunction%wf_type,est_wavefunction%wf_complex)
        elseIf(fileinfo%isGeneral()) then
          nBasis = fileInfo%getVal('nBasis')

          call fileInfo%getArray('OVERLAP',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_wavefunction%overlap_matrix,'overlap','general', &
              tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'OVERLAP not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('CORE HAMILTONIAN ALPHA',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_wavefunction%core_hamiltonian,'core hamiltonian','general', &
              tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'CORE HAMILTONIAN ALPHA not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA ORBITAL ENERGIES',vectorOut=tmpVectorAlpha,foundOut=found)
          if(found) then
            call mqc_matrix_spinBlockGHF(tmpVectorAlpha)
            tmpVectorBeta = tmpVectorAlpha%vat(nBasis+1,-1)
            tmpVectorAlpha = tmpVectorAlpha%vat(1,nBasis)
            call mqc_eigenvalues_allocate(est_wavefunction%mo_energies,'mo energies','general', &
              tmpVectorAlpha,tmpVectorBeta)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA ORBITAL ENERGIES not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA MO COEFFICIENTS',tmpMatrixAlpha,foundOut=found)
          if(found) then
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha,fileInfo%getVal('nElectrons'), &
              fileInfo%getVal('multiplicity'),elist) !MODIFIED
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_wavefunction%mo_coefficients,'mo_coefficients','general', &
              tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
            call est_wavefunction%mo_coefficients%setEList(elist)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA MO COEFFICIENTS not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              if(mqc_matrix_test_symmetric(tmpMatrixAlpha,'hermitian')) &
!                call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
!            tmpMatrixBetaAlpha = MQC_Matrix_Transpose(tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1]))
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_wavefunction%density_matrix,'density','general', &
              tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA DENSITY MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA SCF DENSITY MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              if(mqc_matrix_test_symmetric(tmpMatrixAlpha,'hermitian')) &
!                call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
!            tmpMatrixBetaAlpha = MQC_Matrix_Transpose(tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1]))
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_wavefunction%scf_density_matrix,'density','general', &
              tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA SCF DENSITY MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          call fileInfo%getArray('ALPHA FOCK MATRIX',tmpMatrixAlpha,foundOut=found)
          if(found) then
!            if(MQC_Matrix_HaveComplex(tmpMatrixAlpha)) then
!              call mqc_matrix_symm2full(tmpMatrixAlpha,'hermitian')
!              tmpMatrixAlpha = transpose(tmpMatrixAlpha)
!            endIf
            call mqc_matrix_spinBlockGHF(tmpMatrixAlpha)
            tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
            tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
            tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
            tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
            call mqc_integral_allocate(est_wavefunction%fock_matrix,'fock','general', &
              tmpMatrixAlpha,tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
          else
            if(present(foundObj)) foundObj = .false.
            write(6,'(A)') 'ALPHA FOCK MATRIX not present on file - skipping'
            my_filename = TRIM(fileinfo%filename)
            call fileinfo%CLOSEFILE()
            call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,my_filename)
          endIf

          est_wavefunction%nBasis = fileInfo%getVal('nBasis')
          est_wavefunction%nAlpha = fileInfo%getVal('nAlpha')
          est_wavefunction%nBeta = fileInfo%getVal('nBeta')
          est_wavefunction%nElectrons = fileInfo%getVal('nElectrons')
          est_wavefunction%charge = fileInfo%getVal('charge')
          est_wavefunction%multiplicity = fileInfo%getVal('multiplicity')
          call mqc_gaussian_ICGU(fileInfo%ICGU,est_wavefunction%wf_type,est_wavefunction%wf_complex)
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case default
        call String_Change_Case(label,'u',myLabel)
        If(present(est_eigenvalues)) then
          if(fileinfo%isRestricted()) then
            call fileInfo%getArray('ALPHA '//trim(myLabel),vectorOut=tmpVectorAlpha,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'ALPHA '//trim(myLabel)//' not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
              call mqc_eigenvalues_allocate(est_eigenvalues,'','space',tmpVectorAlpha)
            endIf
          elseIf(fileinfo%isUnrestricted()) then
            call fileInfo%getArray('ALPHA '//trim(myLabel),vectorOut=tmpVectorAlpha,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'ALPHA '//trim(myLabel)//' not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
              call fileInfo%getArray('BETA '//trim(myLabel),vectorOut=tmpVectorBeta,foundOut=found)
              if(present(foundObj)) foundObj = found
              if(.not.found) then
                errorMsg = 'BETA '//trim(myLabel)//' not present on file'
                if(present(foundObj)) then
                  write(6,'(A)') errorMsg
                  call fileinfo%load()
                else
                  call mqc_error_l(trim(errorMsg),6,'found',found)
                endIf
              else
                call mqc_eigenvalues_allocate(est_eigenvalues,'','spin',tmpVectorAlpha, &
                  tmpVectorBeta)
              endIf
            endIf
          elseIf(fileinfo%isGeneral()) then
            call fileInfo%getArray('ALPHA '//trim(myLabel),vectorOut=tmpVectorAlpha,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'ALPHA '//trim(myLabel)//' not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
              nBasis = fileInfo%getVal('nBasis')
              call mqc_matrix_spinBlockGHF(tmpVectorAlpha)
              tmpVectorBeta = tmpVectorAlpha%vat(nBasis+1,-1)
              tmpVectorAlpha = tmpVectorAlpha%vat(1,nBasis)
              call mqc_eigenvalues_allocate(est_eigenvalues,'','general',tmpVectorAlpha, &
                tmpVectorBeta)
            endIf
          else
            call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
                 'fileinfo%isRestricted()', fileinfo%isRestricted(), &
                 'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
                 'fileinfo%isGeneral()', fileinfo%isGeneral() )
          endIf
        elseIf(present(est_integral)) then 
          if(fileinfo%isRestricted()) then
            call fileInfo%getArray('ALPHA '//trim(myLabel),tmpMatrixAlpha,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'ALPHA '//trim(myLabel)//' not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
              call mqc_integral_allocate(est_integral,'','space',tmpMatrixAlpha)
            endIf
          elseIf(fileinfo%isUnrestricted()) then
            call fileInfo%getArray('ALPHA '//trim(myLabel),tmpMatrixAlpha,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'ALPHA '//trim(myLabel)//' not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
              call fileInfo%getArray('BETA '//trim(myLabel),tmpMatrixBeta,foundOut=found)
              if(present(foundObj)) foundObj = found
              if(.not.found) then
                errorMsg = 'BETA '//trim(myLabel)//' not present on file'
                if(present(foundObj)) then
                  write(6,'(A)') errorMsg
                  call fileinfo%load()
                else
                  call mqc_error_l(trim(errorMsg),6,'found',found)
                endIf
              else
                call mqc_integral_allocate(est_integral,'','spin',tmpMatrixAlpha, &
                  tmpMatrixBeta)
              endIf
            endIf
          elseIf(fileinfo%isGeneral()) then
            call fileInfo%getArray('ALPHA '//trim(myLabel),tmpMatrixAlpha,foundOut=found)
            if(present(foundObj)) foundObj = found
            if(.not.found) then
              errorMsg = 'ALPHA '//trim(myLabel)//' not present on file'
              if(present(foundObj)) then
                write(6,'(A)') errorMsg
                call fileinfo%load()
              else
                call mqc_error_l(trim(errorMsg),6,'found',found)
              endIf
            else
              nBasis = fileInfo%getVal('nBasis')
              call mqc_matrix_spinBlockGHF(tmpMatrixAlpha,fileInfo%getVal('nElectrons'), &
                fileInfo%getVal('multiplicity'),elist)
              tmpMatrixBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[nBasis+1,-1])
              tmpMatrixBetaAlpha = tmpMatrixAlpha%mat([1,nBasis],[nBasis+1,-1])
              tmpMatrixAlphaBeta = tmpMatrixAlpha%mat([nBasis+1,-1],[1,nBasis])
              tmpMatrixAlpha = tmpMatrixAlpha%mat([1,nBasis],[1,nBasis])
              call mqc_integral_allocate(est_integral,'','general',tmpMatrixAlpha, &
                tmpMatrixBeta,tmpMatrixAlphaBeta,tmpMatrixBetaAlpha)
              call est_integral%setEList(elist)
            endIf
          else
            call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
                 'fileinfo%isRestricted()', fileinfo%isRestricted(), &
                 'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
                 'fileinfo%isGeneral()', fileinfo%isGeneral() )
          endIf
        else
          call mqc_error('Unsuported open label EST object in getESTObj')
        endIf
      end select

      if(allocated(elist)) deallocate(elist)
!
      return

      end subroutine MQC_Gaussian_Unformatted_Matrix_Get_EST_Object 


!
!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_twoERIs   
      subroutine mqc_gaussian_unformatted_matrix_get_twoERIs(fileinfo,label, &
        est_twoeris,filename,foundERI)
!
!     This subroutine loads the two-electron resonance integrals from a 
!     Gaussian unformatted matrix file sent in object <fileinfo>. The 
!     relevant information will be loaded into output dummy mqc_twoERIs 
!     argument <est_twoeris>.
!
!     Dummy argument <filename> is optional and is only used if fileinfo
!     hasn't already been defined using Routine
!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!     filename sent is different from the filename associated with object
!     fileinfo.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!     called before calling this routine. The expectation is that
!     MQC_Gaussian_Unformatted_Matrix_Read_Header is also called before this
!     routine. However, it is also OK to call this routine first. In that case,
!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     The recognized labels and their meaning include:
!           'regular'            return the regular stored 2ERIs 
!                                  R(i,j,k,l) = (ij|kl)
!           'raffenetti1'        return the raffenetti stored 2ERIs for real
!                                restricted wavefunctions
!                                  R1(i,j,k,l) = (ij|kl) - 1/4[(ik|jl)+(il|jk)]
!           'raffenetti2'        return the raffenetti stored 2ERIs for real
!                                unrestricted wavefunctions
!                                  R2(i,j,k,l) = (ij|kl) + (il|jk)
!           'raffenetti3'        return the raffenetti stored 2ERIs for complex
!                                and/or general wavefunctions
!                                  R3(i,j,k,l) = (ik|jl) - (il|jk)
!           'molecular'          return the molecular orbital basis 2ERIs

!     L. M. Thompson, 2018.
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      logical,optional::foundERI
      type(mqc_twoERIs),optional::est_twoeris
      character(len=*),intent(in),optional::filename
      character(len=256)::my_filename,errorMsg
      character(len=64)::myLabel
      integer::nBasis,nAlpha,nBeta
      type(mqc_r4tensor)::tmpR4TensorAlpha,tmpR4TensorBeta,tmpR4TensorAlphaBeta,tmpR4TensorBetaAlpha
      logical::found
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_L('Error reading Gaussian matrix file header (7): Must include a filename.', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'R' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('regular')
        call fileInfo%getArray('REGULAR 2E INTEGRALS',r4TensorOut=tmpR4TensorAlpha,foundOut=found)
        if(present(foundERI)) foundERI = found
        if(.not.found) then
          errorMsg = 'REGULAR 2E INTEGRALS not present on file'
          if(present(foundERI)) then
            write(6,'(A)') errorMsg
            call fileinfo%load()
          else
            call mqc_error_l(trim(errorMsg),6,'found',found)
          endIf
        else
          call mqc_twoeris_allocate(est_twoeris,'symm','regular',tmpR4TensorAlpha)
        endIf
      case('raffenetti1')
        call fileInfo%getArray('RAFFENETTI 2E INTEGRALS',r4TensorOut=tmpR4TensorAlpha,foundOut=found)
        if(present(foundERI)) foundERI = found
        if(.not.found) then
          errorMsg = 'RAFFENETTI 1 2E INTEGRALS not present on file'
          if(present(foundERI)) then
            write(6,'(A)') errorMsg
            call fileinfo%load()
          else
            call mqc_error_l(trim(errorMsg),6,'found',found)
          endIf
        else
          call mqc_twoeris_allocate(est_twoeris,'symm','raffenetti1',tmpR4TensorAlpha)
        endIf
      case('raffenetti2')
        call fileInfo%getArray('RAFFENETTI 2E INTEGRALS',r4TensorOut=tmpR4TensorAlpha, &
          arraynum=2,foundOut=found)
        if(present(foundERI)) foundERI = found
        if(.not.found) then
          errorMsg = 'RAFFENETTI 2 2E INTEGRALS not present on file'
          if(present(foundERI)) then
            write(6,'(A)') errorMsg
            call fileinfo%load()
          else
            call mqc_error_l(trim(errorMsg),6,'found',found)
          endIf
        else
          call mqc_twoeris_allocate(est_twoeris,'symm','raffenetti2',tmpR4TensorAlpha)
        endIf
      case('raffenetti3')
        call fileInfo%getArray('RAFFENETTI 2E INTEGRALS',r4TensorOut=tmpR4TensorAlpha, &
          arraynum=3,foundOut=found)
        if(present(foundERI)) foundERI = found
        if(.not.found) then
          errorMsg = 'RAFFENETTI 3 2E INTEGRALS not present on file'
          if(present(foundERI)) then
            write(6,'(A)') errorMsg
            call fileinfo%load()
          else
            call mqc_error_l(trim(errorMsg),6,'found',found)
          endIf
        else
          call mqc_twoeris_allocate(est_twoeris,'symm','raffenetti3',tmpR4TensorAlpha)
        endIf
      case('molecular')
        if(fileinfo%isRestricted()) then
          call fileInfo%getArray('Write AA MO 2E INTEGRALS',r4TensorOut=tmpR4TensorAlpha,foundOut=found)
          if(present(foundERI)) foundERI = found
          if(.not.found) then
            errorMsg = 'AA MO 2E INTEGRALS not present on file'
            if(present(foundERI)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call mqc_twoeris_allocate(est_twoeris,'symm','space',tmpR4TensorAlpha)
          endIf
        elseIf(fileinfo%isUnrestricted()) then
          call fileInfo%getArray('Write AA MO 2E INTEGRALS',r4TensorOut=tmpR4TensorAlpha,foundOut=found)
          if(present(foundERI)) foundERI = found
          call fileInfo%getArray('Write BB MO 2E INTEGRALS',r4TensorOut=tmpR4TensorBeta,foundOut=found)
          if(present(foundERI)) found = foundERI.and.found
          call fileInfo%getArray('Write AB MO 2E INTEGRALS',r4TensorOut=tmpR4TensorAlphaBeta,foundOut=found)
          if(present(foundERI)) foundERI = foundERI.and.found
          call fileInfo%getArray('Write BA MO 2E INTEGRALS',r4TensorOut=tmpR4TensorBetaAlpha,foundOut=found)
          if(present(foundERI)) found = foundERI.and.found
          if(.not.found) then
            errorMsg = 'AA MO 2E INTEGRALS and/or BB MO 2E INTEGRALS and/or AB MO 2E INTEGRALS &
              &and/or BA MO 2E INTEGRALS not present on file'
            if(present(foundERI)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call mqc_twoeris_allocate(est_twoeris,'full','spin',tmpR4TensorAlpha,tmpR4TensorBeta, &
              tmpR4TensorAlphaBeta,tmpR4TensorBetaAlpha)
          endIf
        elseIf(fileinfo%isGeneral()) then
          call fileInfo%getArray('Write AA MO 2E INTEGRALS',r4TensorOut=tmpR4TensorAlpha,foundOut=found)
          if(present(foundERI)) foundERI = found
          call fileInfo%getArray('Write BB MO 2E INTEGRALS',r4TensorOut=tmpR4TensorBeta,foundOut=found)
          if(present(foundERI)) foundERI = foundERI.and.found
          call fileInfo%getArray('Write AB MO 2E INTEGRALS',r4TensorOut=tmpR4TensorAlphaBeta,foundOut=found)
          if(present(foundERI)) foundERI = foundERI.and.found
          call fileInfo%getArray('Write BA MO 2E INTEGRALS',r4TensorOut=tmpR4TensorBetaAlpha,foundOut=found)
          if(present(foundERI)) found = foundERI.and.found
          if(.not.found) then
            errorMsg = 'AA MO 2E INTEGRALS and/or BB MO 2E INTEGRALS and/or AB MO 2E INTEGRALS &
              &and/or BA MO 2E INTEGRALS not present on file'
            if(present(foundERI)) then
              write(6,'(A)') errorMsg
              call fileinfo%load()
            else
              call mqc_error_l(trim(errorMsg),6,'found',found)
            endIf
          else
            call mqc_twoeris_allocate(est_twoeris,'full','general',tmpR4TensorAlpha,tmpR4TensorBeta, &
              tmpR4TensorAlphaBeta,tmpR4TensorBetaAlpha)
          endIf
        else
          call mqc_error_L('Unknown wavefunction type in getESTObj', 6, &
               'fileinfo%isRestricted()', fileinfo%isRestricted(), &
               'fileinfo%isUnrestricted()', fileinfo%isUnrestricted(), &
               'fileinfo%isGeneral()', fileinfo%isGeneral() )
        endIf
      case default
        call mqc_error_A('Invalid label sent to %get2ERIs.', 6, &
             'mylabel', mylabel )
      end select
!
      return

      end subroutine MQC_Gaussian_Unformatted_Matrix_Get_twoERIs    


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Write_twoERIs   
      subroutine mqc_gaussian_unformatted_matrix_write_twoERIs(fileinfo,label, &
        ERIs,filename)
!
!     This subroutine writes two-electron resonance integrals to a 
!     Gaussian unformatted matrix file sent in object <fileinfo>. The 
!     relevant information will be loaded from output dummy mqc_twoERIs 
!     argument <est_twoeris>.
!
!     Dummy argument <filename> is optional and is only used if fileinfo
!     hasn't already been defined using Routine
!     MQC_Gaussian_Unformatted_Matrix_Open or if it is determined that the
!     filename sent is different from the filename associated with object
!     fileinfo.
!
!     NOTE: The routine MQC_Gaussian_Unformatted_Matrix_Open is meant to be
!     called before calling this routine. The expectation is that
!     MQC_Gaussian_Unformatted_Matrix_Read_Header is also called before this
!     routine. However, it is also OK to call this routine first. In that case,
!     this routine will first call Routine MQC_Gaussian_Unformatted_Matrix_Open.
!
!     The recognized labels and their meaning include:
!           'regular'            the array contains regular stored 2ERIs 
!           'raffenetti'         the array contains raffenetti stored 2ERIs
!           'molecular'          the array contains molecular orbital basis 2ERIs (NYI)

!     L. M. Thompson, 2020.
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      type(mqc_twoERIs),dimension(:),allocatable::ERIs
      character(len=*),intent(in),optional::filename
      character(len=256)::my_filename
      character(len=64)::myLabel
      type(mqc_r4tensor)::tmpR4Tensor
      real,dimension(:),allocatable::RvectorOut
      real(kind=real64),dimension(:,:,:,:),allocatable::Rr4TenOut
      integer,parameter::LenBuf=4000
      integer::i
      logical::ok
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
               'ok', ok )
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,filename)
        else
          call MQC_Error_L('Error reading Gaussian matrix file header (8): Must include a filename.', 6, &
               'PRESENT(filename)', PRESENT(filename) )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call fileinfo%OPENFILE(TRIM(filename),0,ok)
          if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
            'ok', ok )
          call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'W' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call fileinfo%OPENFILE(TRIM(my_filename),0,ok)
        if(.not.ok) Call MQC_Error_L('Error opening Gaussian matrix file.', 6, &
             'ok', ok )
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          my_filename)
      endIf
      if(.not.fileinfo%header_written) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call fileinfo%OPENFILE(TRIM(my_filename),0,ok)
        call MQC_Gaussian_Unformatted_Matrix_Write_Header(fileinfo,  &
          my_filename)
      endIf
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('regular')
        tmpR4Tensor = ERIs(1)
        if(.not.MQC_R4Tensor_HaveSymmetric(tmpR4Tensor)) call MQC_R4Tensor_Full2Symm(tmpR4Tensor)
        Rr4TenOut = tmpR4Tensor
          if(.not.allocated(RvectorOut)) allocate(RvectorOut(size(Rr4TenOut,1)))
        RvectorOut(:) = Rr4TenOut(:,1,1,1)
        Call Wr_LAO2E(fileinfo%UnitNumber,'REGULAR 2E INTEGRALS',1,LenBuf,size(tmpR4Tensor,1),&
          RvectorOut)
      case('raffenetti')
        do i = 1,size(ERIs)
          tmpR4Tensor = ERIs(i)
          if(.not.MQC_R4Tensor_HaveSymmetric(tmpR4Tensor)) call MQC_R4Tensor_Full2Symm(tmpR4Tensor)
          Rr4TenOut = tmpR4Tensor
          if(.not.allocated(RvectorOut)) allocate(RvectorOut(size(Rr4TenOut,1)*size(ERIs)))
          RvectorOut((i-1)*size(Rr4TenOut,1)+1:) = Rr4TenOut(:,1,1,1)
        endDo
        Call Wr_LAO2E(fileinfo%UnitNumber,'RAFFENETTI 2E INTEGRALS',size(ERIs),LenBuf,size(tmpR4Tensor,1),&
          RvectorOut)
      case('molecular')
        call mqc_error_A('Molecular integrals not yet implemented in %write2ERIs.', 6, &
             'mylabel', mylabel )
      case default
        call mqc_error_A('Invalid label sent to %write2ERIs.', 6, &
             'mylabel', mylabel )
      end select
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Write_twoERIs    



!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_twoERISet   
      subroutine mqc_gaussian_unformatted_matrix_get_twoERISet(fileinfo,label, &
        ERISet,filename,foundERI,replace)
!
!     This subroutine reads two-electron resonance integrals of the requested type
!     from input Gaussian unformatted matrix file sent in object <fileinfo> to
!     output dummy mqc_twoERISet argument <est_twoeriSet>. The recognized labels 
!     and their meaning include:
!           'regular'            return regular 2ERIs 
!           'raffenetti'         return the set of raffenetti 2ERIs present
!           'molecular'          return the molecular orbital basis 2ERIs

!     L. M. Thompson, 2025.
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inOut)::fileinfo
      character(len=*),intent(in)::label
      type(mqc_twoERISet),intent(inOut)::ERISet
      character(len=*),intent(in),optional::filename
      logical,optional,intent(out)::foundERI
      logical,optional,intent(in)::replace
!      
      type(mqc_twoERIs)::ERIs
      character(len=64)::myLabel
      logical::myReplace,found,raffound
!
!     Do the work...
!
      if(present(replace)) then
        myReplace = replace
      else
        myReplace = .true.
      endIf
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('regular')
        if(present(filename)) then
          call fileinfo%get2ERIs('regular',ERIs,filename,found)
        else
          call fileinfo%get2ERIs('regular',ERIs,foundERI=found)
        endIf
        if(.not.found) then
          if(present(foundERI)) then
            foundERI = .false.
            return
          else
            call mqc_error('EXITING DUE TO MISSING INTEGRALS IN &
              &mqc_gaussian_unformatted_matrix_get_twoERISet')
          endIf
        else
          if(present(foundERI)) foundERI = .true.
          call ERISet%addToSet(ERIs,myReplace)
        endIf
      case('raffenetti')
        raffound = .false.
        if(present(filename)) then
          call fileinfo%get2ERIs('raffenetti1',ERIs,filename,found)
        else
          call fileinfo%get2ERIs('raffenetti1',ERIs,foundERI=found)
        endIf
        if(found) then
           raffound = .true.
          if(present(foundERI)) foundERI = .true.
          call ERISet%addToSet(ERIs,myReplace)
        endIf
        if(present(filename)) then
          call fileinfo%get2ERIs('raffenetti2',ERIs,filename,found)
        else
          call fileinfo%get2ERIs('raffenetti2',ERIs,foundERI=found)
        endIf
        if(found) then
           raffound = .true.
          if(present(foundERI)) foundERI = .true.
          call ERISet%addToSet(ERIs,myReplace)
        endIf
        if(present(filename)) then
          call fileinfo%get2ERIs('raffenetti3',ERIs,filename,found)
        else
          call fileinfo%get2ERIs('raffenetti3',ERIs,foundERI=found)
        endIf
        if(found) then
           raffound = .true.
          if(present(foundERI)) foundERI = .true.
          call ERISet%addToSet(ERIs,myReplace)
        endIf
        if(.not.raffound) then
          if(present(foundERI)) then
            foundERI = .false.
            return
          else
            call mqc_error('EXITING DUE TO MISSING INTEGRALS IN &
              &mqc_gaussian_unformatted_matrix_get_twoERISet')
          endIf
        endIf
      case('molecular')
        if(present(filename)) then
          call fileinfo%get2ERIs('molecular',ERIs,filename,found)
        else
          call fileinfo%get2ERIs('molecular',ERIs,foundERI=found)
        endIf
        if(.not.found) then
          if(present(foundERI)) then
            foundERI = .false.
            return
          else
            call mqc_error('EXITING DUE TO MISSING INTEGRALS IN &
              &mqc_gaussian_unformatted_matrix_get_twoERISet')
          endIf
        else
          if(present(foundERI)) foundERI = .true.
          call ERISet%addToSet(ERIs,myReplace)
        endIf
      case default
        call mqc_error_A('Invalid label sent to MQC_Gaussian_Unformatted_Matrix_Get_twoERISet', 6, &
             'mylabel', mylabel )
      end select
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Get_twoERISet    



!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Write_twoERISet   
      subroutine mqc_gaussian_unformatted_matrix_write_twoERISet(fileinfo,label, &
        ERISet,filename)
!
!     This subroutine writes two-electron resonance integrals of the requested type
!     from input dummy mqc_twoERISet argument <est_twoeriSet> to a Gaussian unformatted 
!     matrix file sent in object <fileinfo>. The recognized labels and their meaning 
!     include:
!           'regular'            search for and write regular 2ERIs if present 
!           'raffenetti'         search for and write raffenetti 2ERIs if present
!           'molecular'          the array contains molecular orbital basis 2ERIs (NYI)

!     L. M. Thompson, 2025.
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      type(mqc_twoERISet),intent(in)::ERISet
      character(len=*),intent(in),optional::filename
!      
      type(mqc_twoERIs),dimension(:),allocatable::ERIs
      character(len=64)::myLabel
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('regular')
        if(ERISet%erinum('regular').eq.0) &
          call mqc_error('Attempting to write regular 2ERIs but regular integrals not &
          &present in ERISet in MQC_Gaussian_Unformatted_Matrix_Write_twoERISet')
        allocate(ERIs(1))
        ERIs(1) = ERISet%eris('regular')
        if(present(filename)) then
          call fileinfo%write2ERIs('regular',ERIs,fileName)
        else
          call fileinfo%write2ERIs('regular',ERIs)
        endIf
      case('raffenetti')
        if(ERISet%erinum('raffenetti1').eq.0) call mqc_error('Attempting to write Raffenetti 2ERIs &
          &but Raffenetti 1 integrals not present in ERISet in &
          &MQC_Gaussian_Unformatted_Matrix_Write_twoERISet')
        if(ERISet%erinum('raffenetti2').eq.0.and.ERISet%erinum('raffenetti3').ne.0) &
          call mqc_error('GauOpen expects Raffenetti 2 integrals to be present if Raffenetti 3 &
          &integrals are present in MQC_Gaussian_Unformatted_Matrix_Write_twoERISet')
        if(ERISet%erinum('raffenetti3').ne.0) then
          allocate(ERIs(3))
        elseIf(ERISet%erinum('raffenetti2').ne.0) then
          allocate(ERIs(2))
        else
          allocate(ERIs(1))
        endIf
        ERIs(1) = ERISet%eris('raffenetti1')
        if(size(ERIs).ge.2) ERIs(2) = ERISet%eris('raffenetti2')
        if(size(ERIs).eq.3) ERIs(3) = ERISet%eris('raffenetti3')
        if(present(filename)) then
          call fileinfo%write2ERIs('raffenetti',ERIs,fileName)
        else
          call fileinfo%write2ERIs('raffenetti',ERIs)
        endIf
      case('molecular')
        if(ERISet%erinum('molecular').eq.0) &
          call mqc_error('Attempting to write molecular 2ERIs but molecular integrals not &
          &present in ERISet in MQC_Gaussian_Unformatted_Matrix_Write_twoERISet')
        allocate(ERIs(1))
        ERIs(1) = ERISet%eris('molecular')
        if(present(filename)) then
          call fileinfo%write2ERIs('molecular',ERIs,fileName)
        else
          call fileinfo%write2ERIs('molecular',ERIs)
        endIf
      case default
        call mqc_error_A('Invalid label sent to MQC_Gaussian_Unformatted_Matrix_Write_twoERISet', 6, &
             'mylabel', mylabel )
      end select
!
      return
      end subroutine MQC_Gaussian_Unformatted_Matrix_Write_twoERISet    


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Value_Integer
      Function MQC_Gaussian_Unformatted_Matrix_Get_Value_Integer(fileinfo,label,fileName)
!
!     This function is used to get an integer scalar value that is stored in a
!     Gaussian matrix file object. Available label values that can be sent here
!     include:
!       natoms
!       nbasis
!       nbasisuse
!       charge
!       multiplicity
!       nelectrons
!       nalpha
!       nbeta
!       nshlao
!       nprmao
!
!
!     H. P. Hratchian, 2017.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      character(len=*),intent(in),optional::filename
      integer(kind=int64)::MQC_Gaussian_Unformatted_Matrix_Get_Value_Integer
      integer(kind=int64)::value_out=0
      character(len=64)::myLabel
      character(len=256)::my_filename
      type(MQC_Variable)::mqcTmp
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_L('Failed to retrieve value from Gaussian matrix file: File not open.', 6, &
            'fileinfo%isOpen()', fileinfo%isOpen() )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'R' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
!      if(.not.fileinfo%gaussianScalars_read) then
!        call MQC_Gaussian_Unformatted_Matrix_Read_Array(fileinfo,  &
!          'gaussian scalars',mqcVarOut=mqcTmp)
!        if(Allocated(fileinfo%gaussianScalars))  &
!          call MQC_Error('Logic error loading gaussianScalars.')
!        fileinfo%gaussianScalars = mqcTmp
!      endIf
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('natoms')
        value_out = fileinfo%natoms
      case('nbasis')
        value_out = fileinfo%nbasis
      case('nbasisuse')
        value_out = fileinfo%nbasisUse
      case('charge')
        value_out = fileinfo%icharge
      case('multiplicity')
        value_out = fileinfo%multiplicity
      case('nelectrons')
        value_out = fileinfo%nelectrons
      case('nalpha')
        value_out = (fileinfo%nelectrons + fileinfo%multiplicity - 1)/2 
      case('nbeta')
        value_out = (fileinfo%nelectrons - fileinfo%multiplicity + 1)/2 
      case('nshlao')
        value_out = fileinfo%NShlAO
      case('nprmao')
        value_out = fileinfo%NPrmAO
      case default
        call mqc_error_A('Failure finding requested integer value in Gaussian matrix file', 6, &
          'mylabel', mylabel )
      endSelect
!
      MQC_Gaussian_Unformatted_Matrix_Get_Value_Integer = value_out
      return
      end Function MQC_Gaussian_Unformatted_Matrix_Get_Value_Integer


!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Get_Value_Real
      Function MQC_Gaussian_Unformatted_Matrix_Get_Value_Real(fileinfo,label,fileName)
!
!     This function is used to get an integer scalar value that is stored in a
!     Gaussian matrix file object. Available label values that can be sent here
!     include:
!       energy
!       scfenergy
!       scfs2
!
!
!     H. P. Hratchian, 2023.
!
!
!     Variable Declarations.
!
      implicit none
      class(MQC_Gaussian_Unformatted_Matrix_File),intent(inout)::fileinfo
      character(len=*),intent(in)::label
      character(len=*),intent(in),optional::filename
      real(kind=real64)::MQC_Gaussian_Unformatted_Matrix_Get_Value_Real
      real(kind=real64)::value_out=0.0
      integer,dimension(:),allocatable::GScalars_IPositions
      real,dimension(:),allocatable::GScalars_Reals
      character(len=64)::myLabel
      character(len=256)::my_filename
      type(MQC_Variable)::mqcTmp
!
!
!     Ensure the matrix file has already been opened and the header read.
!
      if(.not.fileinfo%isOpen()) then
        if(PRESENT(filename)) then
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        else
          call MQC_Error_L('Failed to retrieve value from Gaussian matrix file: File not open.', 6, &
            'fileinfo%isOpen()', fileinfo%isOpen() )
        endIf
      endIf
      if(PRESENT(filename)) then
        if(TRIM(filename)/=TRIM(fileinfo%filename)) then
          call fileinfo%CLOSEFILE()
          call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
            filename)
        endIf
      endIf
      if(.not.(fileinfo%readWriteMode .eq. 'R' .or.  &
        fileinfo%readWriteMode .eq. ' ')) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          filename)
      endIf
      if(.not.fileinfo%header_read) then
        my_filename = TRIM(fileinfo%filename)
        call fileinfo%CLOSEFILE()
        call MQC_Gaussian_Unformatted_Matrix_Read_Header(fileinfo,  &
          my_filename)
      endIf
      if(.not.fileinfo%gaussianScalars_read) then
        call MQC_Gaussian_Unformatted_Matrix_Read_Array(fileinfo,  &
          'gaussian scalars',mqcVarOut=mqcTmp)
        if(Allocated(fileinfo%gaussianScalars))  &
          call MQC_Error('Logic error loading gaussianScalars.')
        fileinfo%gaussianScalars = mqcTmp
      endIf
!
!     Do the work...
!
      call String_Change_Case(label,'l',myLabel)
      select case (mylabel)
      case('energy')
        value_out = fileinfo%gaussianScalars(43)
      case('scfenergy')
        value_out = fileinfo%gaussianScalars(32)
      case('scfs2')
        value_out = fileinfo%gaussianScalars(44)
      case default
        call mqc_error_A('Failure finding requested real value in Gaussian matrix file', 6, &
          'mylabel', mylabel )
      endSelect
!
      MQC_Gaussian_Unformatted_Matrix_Get_Value_Real = value_out
      return
      end Function MQC_Gaussian_Unformatted_Matrix_Get_Value_Real

      
!=====================================================================
!
!PROCEDURE MQC_Gaussian_Unformatted_Matrix_Array_Type
      Function MQC_Gaussian_Unformatted_Matrix_Array_Type(NI,NR,N1,N2,N3,N4,N5,NRI,ASym)
!
!     This function returns a character string indicating the type of array
!     found in a Gaussian matrix file. This is done using NI, NR, N1, N2, N3,
!     N4, N5 and NRIfrom a matrix header in a Gaussian unformatted matrix file to
!     determine the type of array the data corresponds to. The return value will
!     be prepended by "REAL-", "INTEGER-", or "COMPLEX-" and appended by one of
!     the following:
!
!           "VECTOR"          A vector.
!           "MATRIX"          A matrix that is allocated full (M x N).
!           "SYMMATRIX"       A symmetric/hermitian matrix.
!           "ASYMMATRIX"      An antisymmetric/antihermitian matrix.
!           "SYMMATRIXN"      Multiple (N) symmetric/hermitian matrix.
!           "ASYMMATRIXN"     Multiple (N) asymmetric/antihermitian matrix.
!
!     If the input flags do not uniquely identify a known array type, then this
!     function returns "UNKNOWN".
!
!
!     H. P. Hratchian, 2017.
!     L. M. Thompson, 2018.
!
!
!     Variable Declarations.
!
      implicit none
      integer::NI,NR,N1,N2,N3,N4,N5,NRI
      logical::Asym
      character(len=64)::MQC_Gaussian_Unformatted_Matrix_Array_Type
!
!
!     Do the work...
!
      MQC_Gaussian_Unformatted_Matrix_Array_Type = "UNKNOWN"
      if(NR.lt.0.or.NI.lt.0) return
      if(NR.gt.0.and.NI.gt.0) then
        MQC_Gaussian_Unformatted_Matrix_Array_Type = "MIXED"
        if(NR.eq.1.and.NI.eq.1) then
          MQC_Gaussian_Unformatted_Matrix_Array_Type = "SCALARS"
        elseIf((NR.eq.1.or.NR.eq.2.or.NR.eq.3).and.NI.eq.4) then
          MQC_Gaussian_Unformatted_Matrix_Array_Type = "2ERIS"
        else
          return
        endIf
      elseIf(NI.eq.1) then 
        MQC_Gaussian_Unformatted_Matrix_Array_Type = "INTEGER"
      elseIf(NR.eq.1) then
        if(NRI.eq.1) then
          MQC_Gaussian_Unformatted_Matrix_Array_Type = "REAL"
        elseIf(NRI.eq.2) then
          MQC_Gaussian_Unformatted_Matrix_Array_Type = "COMPLEX"
        endIf
      endIf
      if(N1.gt.1.and.N2.eq.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
        MQC_Gaussian_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Gaussian_Unformatted_Matrix_Array_Type)//"-VECTOR"
      elseIf(N1.gt.1.and.N2.gt.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
        MQC_Gaussian_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Gaussian_Unformatted_Matrix_Array_Type)//"-MATRIX"
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1.and..not.ASym) then
        MQC_Gaussian_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Gaussian_Unformatted_Matrix_Array_Type)//"-SYMMATRIX"
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1.and.ASym) then
        MQC_Gaussian_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Gaussian_Unformatted_Matrix_Array_Type)//"-ASYMMATRIX"
      elseIf(N1.le.-1.and.N2.le.-1.and.N3.le.-1.and.N4.gt.1.and.N5.eq.1) then
        MQC_Gaussian_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Gaussian_Unformatted_Matrix_Array_Type)//"-SYMSYMR4TENSOR"
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.ge.1.and.N4.eq.1.and.N5.eq.1.and..not.ASym) then
        MQC_Gaussian_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Gaussian_Unformatted_Matrix_Array_Type)//"-SYMMATRIXN"
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.ge.1.and.N4.eq.1.and.N5.eq.1.and.ASym) then
        MQC_Gaussian_Unformatted_Matrix_Array_Type = &
          TRIM(MQC_Gaussian_Unformatted_Matrix_Array_Type)//"-ASYMMATRIXN"
      endIf
!
      return
      end function MQC_Gaussian_Unformatted_Matrix_Array_Type


!=====================================================================


      End Module MQC_Gaussian
