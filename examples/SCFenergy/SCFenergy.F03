! This program tests Reading everthing in a MartixFile
     program SCFenergy
!
!     This program contains examples of the matrix file read and write routines
!

!hph+
      Write(*,*)' Hrant - This example program is obsolescent.'
!!     Variable Declarations...       
!!
!      Use MQC_Algebra
!      Use MQC_DataStructures
!      Use MQC_Files
!      Use MQC_EST
!      use mqc_General
!      use mqc_FullWavefunction
!      use mqc_matwrapper
!      use iso_fortran_env
!      implicit none
!
!      real(kind=real64)::Vnn,EE,EN
!      real(kind=real64)::EE_AA,EN_AA
!      real(kind=real64)::EE_AB,EN_AB
!      real(kind=real64)::EE_BA,EN_BA
!      real(kind=real64)::EE_BB,EN_BB
!      real(kind=real64)::ENL
!      real(kind=real64)::SCFe
!      Complex(kind=real64)::EE_C,EN_C
!      integer(kind=int64) iout
!      integer(kind=int64) I, status
!      integer(kind=int64) Print_Options
!      integer(kind=int64) Storage_Options
!      integer(kind=int64) NBasis
!      integer(kind=int64) NTRI
!      integer(kind=int32) NBasis_32
!      integer(kind=int32) NTRI_32
!      integer(kind=int64) NAtoms
!      integer(kind=int32) J, K, L, IJ, N1, IN, IK, IL
!      integer(kind=int32) JN, KN, JK, LN, JL, KL
!      integer(kind=int32) IJKL, IKJL, ILJK
!
!integer(kind=int32) LMax
!
!      class(mqc_link_list), allocatable:: mqc_input_file_list
!      class(mqc_READONLY_list), allocatable:: my_READONLY_list
!      class(MQC_ALL_Wavefunction),pointer::current_wavefunction
!      type(mqc_array_link_list_node),pointer::DenA
!      type(mqc_array_link_list_node),pointer::DenB
!      type(mqc_array_link_list_node),pointer::HamA
!      type(mqc_array_link_list_node),pointer::HamB
!      type(mqc_array_link_list_node),pointer::FockA
!      type(mqc_array_link_list_node),pointer::FockB
!      type(mqc_array_link_list_node),pointer::Integrals
!
!      Real(Kind=Real64),Pointer,contiguous::Density_aaRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Density_abRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Density_aaIA(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Density_abRA(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Density_bbRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Density_abIS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Density_bbIA(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Density_abIA(:) =>null()
!
!      Real(Kind=Real64),Pointer,contiguous::Den2_aaRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Den2_abRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Den2_bbRS(:) =>null()
!
!      Real(Kind=Real64),Pointer,contiguous::Ham_aaRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ham_abRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ham_aaIA(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ham_abRA(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ham_bbRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ham_abIS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ham_bbIA(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ham_abIA(:) =>null()
!
!      Real(Kind=Real64),Pointer,contiguous::AlphaDensity(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::BetaDensity(:) =>null()
!
!      Real(Kind=Real64),Pointer,contiguous::Int2e(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Coulomb(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ex_aaRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ex_bbRS(:) =>null()
!      Real(Kind=Real64),Pointer,contiguous::Ex_abRS(:) =>null()
!
!      Integer(kind=int32)::NDIM,NMatI,NTT, ntt2
!
!      integer(kind=int64) Found_DenA,Found_DenB,Found_HamA,Found_HamB
!      integer(kind=int64) Found_FockA,Found_FockB,Found_2e_ints
!      character(len=8)::tab
!      logical dofock
!!
!!     Format statements.
!!
! 1010 Format( A )
! 1020 Format( "Matrix File:  ",A )
! 1030 Format( A,I6 )
! 1200 format( "   Electronic Energy:",a,a,E20.10,/, &
!           "   Electronic/Nuclear Energy:",a,E20.9,/, &
!           "   Nuclear Energy:",a,a,E20.10,/, &
!           "   Total Energy:",a,a,E20.10)
!!
!! To run this program plese invoke as:
!!
!!    SCFenergy {Matrix Element File}
!!            OR
!!    SCFenergy {Matrix Element File 1} ... {Matrix Element File 2}
!!            OR
!!    SCFenergy {Gaussian Input File} {Gaussian Binary}
!!            OR
!!    SCFenergy {Gaussian Input File} {Matrix Element File} {Gaussian Binary}
!!
!! When running in the form that accepts Gaussian input files:
!!    - You must set up the Gaussian environment before you run
!!    - Multiple Gaussian inputs or a mixture of Gaussian inputs and Matrix 
!!      Element Files  may be supplied.  The Gaussian binary must always be 
!!      the last argument.
!!
!!
!!     Now start the work:
!!
!! We want to write to standard out, unit 6
!!
!! The value needs to go into an 8 Byte Integer variable.  It will be passed 
!! to various subroutines.  The default size of Integers is 4 Bytes, passing
!! '6' in an argument list can cause problems.
!!
!      IOut = 6
!!
!! Matrix Element files do not contain information about functionals.  So
!! this program will not give proper energies for DFT calculations.  
!! Unfortunately, since Matrix Elements files contain nothing about functionals,
!! 
!!
!      write( IOut, 1010) "Warning, SCFenergy gives incorrect energies for DFT wavefunctions"
!!
!!     Create a list of data structures to try to read.
!!
!! The Options of this can be found by:
!!   In the output of the Gaussian Job that wrote the Matrix Element File.
!!
!!            OR
!!
!!   Run the example program DataSummary.  The executable can be found in 
!!   ${MQC_Home}/examples/DataSummary/DataSummary
!!   The arguments options for DataSummary are the as the options for this
!!   program.
!!
!      allocate(my_READONLY_list)
!! This data is needed for Nuclear Replusion
!      call mqc_READONLY_list_add(my_READONLY_list, &
!           "CARTESIAN NUCLEAR COORDINATES(BOHR)")
!      call mqc_READONLY_list_add(my_READONLY_list,"ATOM CHARGES")
!! This data is needed for components that include electronic interactions
!! There are always Alpha spin versions of these.
!      call mqc_READONLY_list_add(my_READONLY_list,"ALPHA SCF DENSITY MATRIX")
!      call mqc_READONLY_list_add(my_READONLY_list,"ALPHA DENSITY MATRIX")
!      call mqc_READONLY_list_add(my_READONLY_list,"CORE HAMILTONIAN ALPHA")
!      call mqc_READONLY_list_add(my_READONLY_list,"ALPHA FOCK MATRIX")
!! There are sometimes spin versions of the matrixies.  If they are in the
!! Matrix Element File, they will be saved.  These requests, not demands.  If
!! the data is not there, the request will be ignored.
!      call mqc_READONLY_list_add(my_READONLY_list,"BETA SCF DENSITY MATRIX")
!      call mqc_READONLY_list_add(my_READONLY_list,"BETA DENSITY MATRIX")
!      call mqc_READONLY_list_add(my_READONLY_list,"CORE HAMILTONIAN BETA")
!      call mqc_READONLY_list_add(my_READONLY_list,"BETA FOCK MATRIX")
!      call mqc_READONLY_list_add(my_READONLY_list,"REGULAR 2E INTEGRALS")
!!
!! Now start the heavy lifting.
!!
!! MQC_Create_Wavefunctions first deals with the input(s).
!!   - Matrix Element Files will be read.
!!   - Gaussian Input will be modified, so that Matrix Element files are 
!!     created.  And then the Gaussian job will be run. The resulting
!!     Matrix Element File will be read.
!!
!! The resulting structures that contain the information from the Matrix Element
!! files are found in a linked list.
!!
!
!!
!! Create the linked list that will contain the information saved from the
!! Matrix Element File.
!!
!      allocate(mqc_input_file_list)
!!
!! Select the option for printing.  The currently implemented options are:
!!   0 - Silent mode
!!   1 - Print a summary of the data in the Matrix Element File
!!   2 - Print the values of all data read
!!
!      Print_Options=0
!!
!! Select the option for storage.  The currently implemented options are:
!!   0 - Do not save the values read from the Matrix Element File
!!   1 - Save the values read from the Matrix Element File
!!
!      Storage_Options=1
!!
!! Read the Matrix Element Files or Generate Matrix Element Files and then
!! Read them.
!! Store the 
!      Call MQC_Create_Wavefunctions( Print_Options, Storage_Options, &
!           my_READONLY_list, mqc_input_file_list, Iout)
!!
!! Now the Matrix Element files have been read.
!! Process the information from each file that has been read.
!!
!#ifdef MATRIXFILE
!      call mqc_link_list_start(mqc_input_file_list,Status,current_wavefunction)
!!
!! If no information has been sotred, Status will be 0.
!! If something has been read, the information is stored in current_wavefunction
!!
!      Do While( Status.eq.1 )
!!
!! Print out the name of the Matrix Element file being processed.  It's nice
!! to have a reminder of what you are looking at.
!!
!         write(iout,1020) trim(current_wavefunction%FileName)
!!
!! Assign the 2 most commonly used integers used for dimensions of arrays:
!!   NBasis - the number of basis functions
!!   NAtoms - the number of atoms (Not used here, but often used)
!!
!         NBasis = current_wavefunction%nbasis
!         NAtoms = current_wavefunction%natoms
!!
!! Get Nuclear repulsion Energy
!!
!         Vnn = MQC_FULL_Wave_To_Nuclear_Repulsion(current_wavefunction,IOut)
!!
!! Get all possible components of Electronic part of energy
!!
!! Get the values read
!!
!! The data is stored with memory allocated with a pointer.  So there is a
!! single copy residing in memory.  The second to last argument is the pointer.
!!
!! The return argument gives the success of the getting the data as well as the
!! type of data.  The possible values are:
!!   0 - No Data found
!!   1 - Integer Data
!!   2 - Real Data
!!   3 - Complex Data
!!   4 - Mixed Data - both Integer and Real
!!   5 - Mixed Data - both Integer and Complex
!!
!         Found_DenA=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              "ALPHA SCF DENSITY MATRIX", DenA, IOut)
!! So, if data of the name "ALPHA SCF DENSITY MATRIX" is not found, try another
!! name
!         If ( Found_DenA .eq. 0 ) then
!            Found_DenA=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!                 "ALPHA DENSITY MATRIX", DenA, IOut)
!         EndIf
!         Found_DenB=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              "BETA SCF DENSITY MATRIX", DenB, IOut)
!         If ( Found_DenB .eq. 0 ) then
!            Found_DenB=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!                 "BETA DENSITY MATRIX", DenB, IOut)
!         EndIf
!         Found_HamA=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              "CORE HAMILTONIAN ALPHA", HamA, IOut)
!         Found_HamB=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              "CORE HAMILTONIAN BETA", HamB, IOut)
!         Found_FockA=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              "ALPHA FOCK MATRIX", FockA, IOut)
!         Found_FockB=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
!              "BETA FOCK MATRIX", FockB, IOut)
!!
!! Find the number of elements of these Triangular matricies
!!
!         NTRI = (NBasis+1)*NBasis/2
!         NTRI_32 = NTRI
!!
!! First Check if Real Data for the beta spin arrays is found.
!!
!         If ( Found_DenB .eq. 2 .and. Found_HamB .eq. 2 .and. &
!              Found_FockB .eq. 2 ) then
!            write(IOut, 1010) "UHF or ROHF"
!            flush(IOut)
!!
!! The easy way to do a Triangular Matrix/Matrix multiply is to double the 
!! off diagonal elements of 1 array and then do a dot product.  
!! It's computationally a bit faster to half the diagonal elements of the 
!! density Matrix and adjust result result as needed.
!!
!! Allocate Alpha and Beta Density arrays.  These are needed, because we are
!! going to double the off-Diagonal elements.  It's bad practice to modify
!! data pointed to with a pointer, unless the change will be required
!! by everything that follows it.
!!
!            call MQC_Alloc_Real_Ptr( AlphaDensity, &
!                 NTRI, 'Alpha Density', IOut)
!            call MQC_Alloc_Real_Ptr( BetaDensity, &
!                 NTRI, 'Beta Density', IOut)
!!
!            NBasis_32 = NBasis
!            Call Double_off_Diag( DenA%ArrR, AlphaDensity, NBasis_32 )
!            Call Double_off_Diag( DenB%ArrR, BetaDensity, NBasis_32 )
!
!            EE_AA = dot_product( AlphaDensity, FockA%ArrR )
!            EN_AA = dot_product( AlphaDensity, HamA%ArrR )
!            EE_BB = dot_product( BetaDensity, FockB%ArrR )
!            EN_BB = dot_product( BetaDensity, HamB%ArrR )
!
!            EE = (EE_AA+EE_BB)*0.5D0
!            EN = (EN_AA+EN_BB)*0.5D0
!            SCFe = EE + EN + Vnn
!         Else If ( Found_DenA .eq. 2 .and. Found_HamA .eq. 2 .and. &
!              Found_FockA .eq. 2 ) then
!!
!! Half Diagonal elements of Density
!!
!            write(IOut, 1010) "RHF"
!            flush(IOut)
!
!            call MQC_Alloc_Real_Ptr( AlphaDensity, &
!                 NTRI, 'Alpha Density', IOut)
!
!            NBasis_32 = NBasis
!            Call Double_off_Diag( DenA%ArrR, AlphaDensity, NBasis_32 )
!            EE = dot_product( AlphaDensity, FockA%ArrR )
!            EN = dot_product( AlphaDensity, HamA%ArrR )
!
!            SCFe = EE + EN + Vnn
!         Else If ( Found_DenB .eq. 3 .and. Found_HamB .eq. 3 .and. &
!              Found_FockB .eq. 3 ) then
!            call MQC_Error_A( 'ERROR: Found a Complex Beta Spin Data', IOut )
!
!! scalar = dot_product(vector1, vector2)
!
!         Else If ( Found_DenA .eq. 3 .and. Found_HamA .eq. 3 .and. &
!              Found_FockA .eq. 3 ) then
!!
!! Half Diagonal elements of Density
!!
!! Array sizes for this seem to be based on 2*NBasis, used sizes stored with arrays
!            write(IOut, 1010) "GHF"
!            flush(IOut)
!
!            Found_2e_ints=MQC_Find_mqc_array_link_list_node( &
!                 current_wavefunction, "REGULAR 2E INTEGRALS", Integrals, IOut)
!            If ( Found_2e_ints .ne. 2 ) Then
!               write(IOut, 1010) 'Warning SCFenergy: Did not find 2e &
!                    &  Integrals.  GHF energies require 2 Electron Integrals'
!               write(IOut, 1030) '  Found_2e_ints = ', Found_2e_ints
!               SCFe = 0.0
!               EE = 0.0
!               EN = 0.0
!               Vnn = 0.0
!            Else
!
!               call MQC_Alloc_Real_Ptr( Den2_aaRS, NTRI, 'Den2_aaRS', IOut)
!               call MQC_Alloc_Real_Ptr( Den2_abRS, NTRI, 'Den2_abRS', IOut)
!               call MQC_Alloc_Real_Ptr( Den2_bbRS, NTRI, 'Den2_bbRS', IOut)
!
!               call MQC_Alloc_Real_Ptr( Density_aaRS, NTRI, 'Density_aaRS', &
!                    IOut)
!               call MQC_Alloc_Real_Ptr( Density_abRS, NTRI, 'Density_abRS', &
!                    IOut)
!               call MQC_Alloc_Real_Ptr( Density_aaIA, NTRI, 'Density_aaIA', &
!                    IOut)
!               call MQC_Alloc_Real_Ptr( Density_abRA, NTRI, 'Density_abRA', &
!                    IOut)
!               call MQC_Alloc_Real_Ptr( Density_bbRS, NTRI, 'Density_bbRS', &
!                    IOut)
!               call MQC_Alloc_Real_Ptr( Density_abIS, NTRI, 'Density_abIS', &
!                    IOut)
!               call MQC_Alloc_Real_Ptr( Density_bbIA, NTRI, 'Density_bbIA', &
!                    IOut)
!               call MQC_Alloc_Real_Ptr( Density_abIA, NTRI, 'Density_abIA', &
!                    IOut)
!!     NDim is the dimension after
!!     NTT of the separated arrays,
!               NDIM = 1
!               NMatI = 1
!               NTT = (NBasis*(NBasis+1))/2
!               ntt2 = (NBasis*(NBasis*2+1))
!               NBasis_32 = NBasis
!!
!! Extract the 8 Density Triangular arrays from the Complex Matrix
!!
!               Call MQC_GHF_Complex_to_LT(DoFock,NDim,NMatI,NBasis_32, &
!                    NTT,NTT2,DenA%ArrC, &
!                    Density_aaRS,Density_abRS,Density_aaIA,Density_abRA, &
!                    Density_bbRS,Density_abIS,Density_bbIA,Density_abIA )
!!
!! Create the Density matricies with doubled off diagonal values
!!
!               Call Double_off_Diag( Density_aaRS, Den2_aaRS, NBasis_32 )
!               Call Double_off_Diag( Density_abRS, Den2_abRS, NBasis_32 )
!               Call Double_off_Diag( Density_bbRS, Den2_bbRS, NBasis_32 )
!               call MQC_Alloc_Real_Ptr( AlphaDensity, NTRI, 'Alpha Density', &
!                    IOut)
!! Sum the alpha and beta densities for use with e-e replusion
!               AlphaDensity = Den2_aaRS + Den2_bbRS
!!
!! Allocate the 8 Hamiltonian arrays
!!
!               call MQC_Alloc_Real_Ptr( Ham_aaRS, NTRI, 'Ham_aaRS', IOut)
!               call MQC_Alloc_Real_Ptr( Ham_abRS, NTRI, 'Ham_abRS', IOut)
!               call MQC_Alloc_Real_Ptr( Ham_aaIA, NTRI, 'Ham_aaIA', IOut)
!               call MQC_Alloc_Real_Ptr( Ham_abRA, NTRI, 'Ham_abRA', IOut)
!               call MQC_Alloc_Real_Ptr( Ham_bbRS, NTRI, 'Ham_bbRS', IOut)
!               call MQC_Alloc_Real_Ptr( Ham_abIS, NTRI, 'Ham_abIS', IOut)
!               call MQC_Alloc_Real_Ptr( Ham_bbIA, NTRI, 'Ham_bbIA', IOut)
!               call MQC_Alloc_Real_Ptr( Ham_abIA, NTRI, 'Ham_abIA', IOut)
!!
!! Extract the 8 Hamiltonian Triangular arrays from the Complex Matrix
!!
!               Call MQC_GHF_Complex_to_LT(DoFock,NDim,NMatI,NBasis_32,&
!                    NTT,NTT2,HamA%ArrC, &
!                    Ham_aaRS,Ham_abRS,Ham_aaIA,Ham_abRA, &
!                    Ham_bbRS,Ham_abIS,Ham_bbIA,Ham_abIA )
!!
!! Calculate the Electron-Nucleus Energy, this just requires 2 
!! Hamiltonian Matricies
!!
!               EN = dot_product( Ham_aaRS, Den2_aaRS ) + &
!                    dot_product( Ham_bbRS, Den2_bbRS )
!!
!! Deallocate the 8 Hamiltonian arrays
!!
!               call MQC_DeAlloc_Real_Ptr( Ham_aaRS, 'Ham_aaRS', IOut)
!               call MQC_DeAlloc_Real_Ptr( Ham_abRS, 'Ham_abRS', IOut)
!               call MQC_DeAlloc_Real_Ptr( Ham_aaIA, 'Ham_aaIA', IOut)
!               call MQC_DeAlloc_Real_Ptr( Ham_abRA, 'Ham_abRA', IOut)
!               call MQC_DeAlloc_Real_Ptr( Ham_bbRS, 'Ham_bbRS', IOut)
!               call MQC_DeAlloc_Real_Ptr( Ham_abIS, 'Ham_abIS', IOut)
!               call MQC_DeAlloc_Real_Ptr( Ham_bbIA, 'Ham_bbIA', IOut)
!               call MQC_DeAlloc_Real_Ptr( Ham_abIA, 'Ham_abIA', IOut)
!!
!! Allocate arrays for Coulomb and Exchange
!!
!               call MQC_Alloc_Real_Ptr( Coulomb, NTRI, 'Coulomb', IOut)
!               call MQC_Alloc_Real_Ptr( Ex_aaRS, NTRI, 'Ex_aaRS', IOut)
!               call MQC_Alloc_Real_Ptr( Ex_bbRS, NTRI, 'Ex_bbRS', IOut)
!               call MQC_Alloc_Real_Ptr( Ex_abRS, NTRI, 'Ex_abRS', IOut)
!! Initialize the Triangular Matricies
!               Coulomb = 0.0D0
!               Ex_aaRS = 0.0D0
!               Ex_bbRS = 0.0D0
!               Ex_abRS = 0.0D0
!!
!! The exchange loop needs work to simplify.  This works.
!!
!            IJ = 0
!            IN = 0
!            Do I = 1,NBasis
!               JN = 0
!               Do J = 1,I
!                  IJ = IJ + 1
!                  KN = 0
!                  Do K = 1,NBasis
!                     if ( K .gt. I ) THen
!                        IK = KN + I
!                     Else
!                        IK = IN + K
!                     EndIf
!                     LN = 0
!                     Do L = 1,NBASIS
!                        if ( L .gt. K ) THen
!                           KL = LN + K
!                        Else
!                           KL = KN + L
!                        EndIf
!                        if ( L .gt. J ) THen
!                           JL = LN + J
!                        Else
!                           JL = JN + L
!                        EndIf
!                        if ( IK .gt. JL ) Then
!                           IKJL = (IK*(IK-1))/2+JL
!                        Else
!                           IKJL = (JL*(JL-1))/2+IK
!                        Endif
!                        Ex_aaRS(IJ) = Ex_aaRS(IJ) + &
!                             Integrals%ArrR(IKJL)*Density_aaRS(KL)
!                        Ex_bbRS(IJ) = Ex_bbRS(IJ) + &
!                             Integrals%ArrR(IKJL)*Density_bbRS(KL)
!                        Ex_abRS(IJ) = Ex_abRS(IJ) + &
!                             Integrals%ArrR(IKJL)*Density_abRS(KL)
!                        LN = LN + L
!                     EndDo
!                     KN = KN + K
!                  EndDo
!                  JN = JN + J
!               EndDo
!               IN = IN + I
!            EndDo
!            IJ = 0
!            IJKL = 0
!            Do I = 1,NBasis
!               Do J = 1,I
!                  IJ = IJ + 1
!                  KL = 0
!                  Do K = 1,I-1
!                     Do L = 1,K
!                        KL = KL + 1
!                        IJKL = IJKL + 1
!                        Coulomb(IJ) = Coulomb(IJ) + &
!                             Integrals%ArrR(IJKL)*AlphaDensity(KL)
!                        Coulomb(KL) = Coulomb(KL) + &
!                             Integrals%ArrR(IJKL)*AlphaDensity(IJ)
!                     EndDo
!                  EndDo
!!                 K = I is a special case, because it includes IJ = KL
!                  Do L = 1,J - 1
!                     KL = KL + 1
!                     IJKL = IJKL + 1
!                     Coulomb(IJ) = Coulomb(IJ) + &
!                          Integrals%ArrR(IJKL)*AlphaDensity(KL)
!                     Coulomb(KL) = Coulomb(KL) + &
!                          Integrals%ArrR(IJKL)*AlphaDensity(IJ)
!                  EndDo
!!                 Integral symmemtry leads to 1 integral for IJ = KL
!                  IJKL = IJKL + 1
!                  Coulomb(IJ) = Coulomb(IJ) + &
!                       Integrals%ArrR(IJKL)*AlphaDensity(IJ)
!               EndDo
!            EndDo
!            EE = dot_product(Coulomb,AlphaDensity ) &
!                 - dot_product(Ex_aaRS,Den2_aaRS ) &
!                 - dot_product(Ex_bbRS,Den2_bbRS ) &
!                 - 2.0*dot_product(Ex_abRS,Den2_abRS )
!
!
!#ifdef DEBUG
!!
!! This data is for runs with test1132.  
!! To get the required 2 electron integrals in the Matrix Element File, 
!! add "scf=conven int=noraff" to the keywords for the last job in
!! examples/data/Gaussian_input/test1132.com 
!!
!            write( 6, 678)
!678         format( " EE should be 83.56814677766")
!            write(6,1432) EE
!1432        format( "Calculated EE",f20.10)
!            write( 6, 543 ) -149.610171147 - 83.56814677766 - Vnn
!543         format( "EN should be ",F20.10 )
!            write( 6,7654 ) EN
!7654        Format( 'Calculated EN',F20.10)
!#endif
!
!!
!! Deallocate arrays for 2e ints
!!
!            call MQC_DeAlloc_Real_Ptr( Coulomb, 'Coulomb', IOut)
!            call MQC_DeAlloc_Real_Ptr( Ex_aaRS, 'Ex_aaRS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Ex_bbRS, 'Ex_bbRS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Ex_abRS, 'Ex_abRS', IOut)
!!
!! Deallocate arrays for Density
!!
!            call MQC_DeAlloc_Real_Ptr( Den2_aaRS, 'Den2_aaRS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Den2_abRS, 'Den2_abRS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Den2_bbRS, 'Den2_bbRS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Density_aaRS, 'Density_aaRS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Density_abRS, 'Density_abRS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Density_aaIA, 'Density_aaIA', IOut)
!            call MQC_DeAlloc_Real_Ptr( Density_abRA, 'Density_abRA', IOut)
!            call MQC_DeAlloc_Real_Ptr( Density_bbRS, 'Density_bbRS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Density_abIS, 'Density_abIS', IOut)
!            call MQC_DeAlloc_Real_Ptr( Density_bbIA, 'Density_bbIA', IOut)
!            call MQC_DeAlloc_Real_Ptr( Density_abIA, 'Density_abIA', IOut)
!            call MQC_DeAlloc_Real_Ptr( AlphaDensity, 'AlphaDensity', IOut)
!
!            SCFe = EE + EN + Vnn
!            EndIf
!         Else
!            call mqc_error_I('SCFenergy: Do not have a case for these flags &
!                 & for components of Electrical components of Energy', &
!                 IOut,&
!                 'Found_DenA', Found_DenA, 'Found_DenB', Found_DenB, &
!                 'Found_HamA', Found_HamA, 'Found_HamB', Found_HamB, &
!                 'Found_FockA', Found_FockA, 'Found_FockB', Found_FockB )
!         EndIf
!         tab = char( 9 )
!         write(6, 1200)tab, tab, EE, tab, EN, tab, tab, Vnn, tab, tab, SCFe
!         call mqc_link_list_next(mqc_input_file_list,Status, &
!              current_wavefunction)
!      EndDo
!#endif


      end program SCFenergy


!    Subroutine DouDot( DResult, Density, Other, N )
!      use iso_fortran_env
!      implicit none
!
!      Real(Kind=real64),Intent(Out)::DResult
!      integer(kind=int32),Intent(In)::N
!      Real(Kind=real64),dimension(N),Intent(In)::Density
!      Real(Kind=real64),dimension(N),Intent(In)::Other
!
!      integer(kind=int32)::K
!      DResult = 0.0D0
!      DO K = 1, N
!         DResult = DResult + Density(K) * Other(K)
!      EndDo
!      Return
!    End Subroutine DouDot
!
!    Subroutine CompDot( DResult, Density, Other, NBasis )
!      use iso_fortran_env
!      implicit none
!
!      Real(Kind=real64),Intent(Out)::DResult
!      integer(kind=int32),Intent(In)::NBasis
!      Complex(Kind=real64),dimension(NBasis,NBasis),Intent(In)::Density
!      Complex(Kind=real64),dimension(NBasis,NBasis),Intent(In)::Other
!      integer(kind=int32)::I, J
!      DResult = 0.0D0
!      DO J = 1, NBasis
!         DO I = 1, NBasis
!            DResult = DResult + Real(Density(I,J)) * Real(Other(I,J))
!
!!         write(6, 1200) REAL(Density(I,J)), REAL(Other(I,J)), I,J, DResult
!! 1200   format( 2E20.8,2I5, 2E20.8)
!
!         EndDo
!      EndDo
!      Return
!    End Subroutine CompDot
!
!    Subroutine Double_off_Diag( Density_In, Density_Out, NBasis )
!      use iso_fortran_env
!      implicit none
!
!      integer(kind=int32),Intent(In)::NBasis
!      Real(Kind=real64),dimension(((NBasis+1)*NBasis)/2),Intent(In)::Density_In
!      Real(Kind=real64),dimension(((NBasis+1)*NBasis)/2),Intent(Out)::Density_Out
!      integer(kind=int32)::I, J, IJ
!!
!! Here's why we double the off diagonal elements of the density:
!!
!!  E = DotProduct(D,H) + DotProduct(D,F) + other stuff...
!!
!! Where D is the Density.  Note that there are 2 terms that use D.
!! The Dot Product of D and H is:
!!
!! Dprod = 0.
!! Do I=1->NBasis{ Do J=1->NBasis{ Dprod = Dprod + D(I,J)*H(I,J) } }
!!
!! Now seperate off diagonal from Diagonal terms:
!!
!! Dprod = 0.
!! Do I=1->NBasis{ Do J=1->I-1{ Dprod = Dprod + D(I,J)*H(I,J) } }
!! Do I=1->NBasis{ Dprod = Dprod + D(I,I)*H(I,I) }
!! Do I=1->NBasis{ Do J=I+1->NBasis{ Dprod = Dprod + D(I,J)*H(I,J) } }
!!
!! Now reorder the loops on the last term.
!!
!! Dprod = 0.
!! Do I=1->NBasis{ Do J=1->I-1{ Dprod = Dprod + D(I,J)*H(I,J) } }
!! Do I=1->NBasis{ Dprod = Dprod + D(I,I)*H(I,I) }
!! Do I=1->NBasis{ Do J=1->I-1{ Dprod = Dprod + D(J,I)*H(J,I) } }
!!
!! Modify last term using D(I,J) = D(J,I) and H(I,J) = H(J,I)
!!
!! Dprod = 0.
!! Do I=1->NBasis{ Do J=1->I-1{ Dprod = Dprod + D(I,J)*H(I,J) } }
!! Do I=1->NBasis{ Dprod = Dprod + D(I,I)*H(I,I) }
!! Do I=1->NBasis{ Do J=1->I-1{ Dprod = Dprod + D(I,J)*H(I,J) } }
!!
!! Since the first and last terms are the same:
!!
!! Dprod = 0.
!! Do I=1->NBasis{ Do J=1->I-1{ Dprod = Dprod + 2.0*D(I,J)*H(I,J) } }
!! Do I=1->NBasis{ Dprod = Dprod + D(I,I)*H(I,I) }
!!
!! This routine doubles the off diagonal elements of the density so, that 
!! about half the terms and half the multiplies are used to get each dot 
!! product
!
!      IJ = 0
!      DO I = 1, NBasis
!         DO J = 1, I-1
!            IJ = IJ+1
!            Density_Out(IJ) = 2.0D0*Density_In(IJ)
!         EndDo
!         IJ = IJ+1
!         Density_Out(IJ) = Density_In(IJ)
!      EndDo
!
!      Return
!    end Subroutine Double_off_Diag
!
!      Subroutine MQC_GHF_Complex_to_LT(DoFock,NDim,NMatI,N,NTT,NTT2,XA, &
!           VaaRS,VabRS,VaaIA,VabRA,VbbRS,VabIS,VbbIA,VabIA)
!      use iso_fortran_env
!      Implicit none
!!
!!     Convert a GHF matrix in XA from complex to lower-triangular real
!!     components in XA and XB.  Spin-orbital basis functions on input
!!     are assumed to be 1a, 1b, 2a, 2b, ... and XA is Hermetian and
!!     lower-triangular.
!!
!!     The order generated is:
!!     XA:  PaaRS, PabRS, PaaIA, PabRA.
!!     XB:  PbbRS, PabIS, PbbIA, PabIA.
!!
!      Logical DoFock
!      integer(kind=int32),intent(in)::NDim,NMatI,N,NTT,NTT2
!      complex(kind=real64),intent(in)::XA(NTT2,NMatI)
!      real(kind=real64),intent(out)::VaaRS(NTT,NDim)
!      real(kind=real64),intent(out)::VabRS(NTT,NDim)
!      real(kind=real64),intent(out)::VaaIA(NTT,NDim)
!      real(kind=real64),intent(out)::VabRA(NTT,NDim)
!      real(kind=real64),intent(out)::VbbRS(NTT,NDim)
!      real(kind=real64),intent(out)::VabIS(NTT,NDim)
!      real(kind=real64),intent(out)::VbbIA(NTT,NDim)
!      real(kind=real64),intent(out)::VabIA(NTT,NDim)
!
!      integer(kind=int32)::I,J,IM,JJ,JA,JB,JJA,JJB,JAIB,IJA,IJB
!      integer(kind=int32)::IJ,IA,IB,IAJA,IAJb
!      real(kind=real64)::FactI,FABI
!!
!!     Same-spin cases.
!!
!      If(DoFock) then
!         FactI = -1.0D0
!         FABI = -0.5D0
!      else
!         FactI = 1.0D0
!         FABI = 0.5D0
!      endIf
!      Do IM = 1, NMatI
!         Do J = 1, N
!            JA = 2*J - 1
!            JB = JA + 1
!            JJ = (J*(J-1))/2
!            JJA = (JA*(JA-1))/2
!            JJB = (JB*(JB-1))/2
!            Do I = 1, J
!               IA = 2*I - 1
!               IB = IA + 1
!               IJ = JJ + I
!               IJA = JJA + IA
!               IJB = JJB + IB
!               VaaRS(IJ,IM) = Real(XA(IJA,IM))
!               VaaIA(IJ,IM) = FactI*aImag(XA(IJA,IM))
!               VbbRS(IJ,IM) = Real(XA(IJB,IM))
!               VbbIA(IJ,IM) = FactI*aImag(XA(IJB,IM))
!            EndDo
!         EndDo
!!
!!       Mixed-spin.
!!
!         Do J = 1, N
!            JJ = (J*(J-1))/2
!            JA = 2*J - 1
!            JB = JA + 1
!            JJA = (JA*(JA-1))/2
!            JJB = (JB*(JB-1))/2
!            Do I = 1, (J-1)
!               IA = 2*I - 1
!               IB = IA + 1
!               IAJB = JJB + IA
!               JAIB = JJA + IB
!               IJ = JJ + I
!               VabRS(IJ,IM) = 0.5D0*(Real(XA(IAJB,IM))+Real(XA(JAIB,IM)))
!               VabRA(IJ,IM) = FABI*(Real(XA(IAJB,IM))-Real(XA(JAIB,IM)))
!               VabIS(IJ,IM) = 0.5D0*(aImag(XA(JAIB,IM))-aImag(XA(IAJB,IM)))
!               VabIA(IJ,IM) = FABI*(aImag(XA(JAIB,IM))+aImag(XA(IAJB,IM)))
!            EndDo
!            IJ = JJ + J
!            JAIB = JJB + JA
!            VabRS(IJ,IM) = Real(XA(JAIB,IM))
!            VabRA(IJ,IM) = 0.0D0
!            VabIS(IJ,IM) = -aImag(XA(JAIB,IM))
!            VabIA(IJ,IM) = 0.0D0
!         EndDo
!      EndDo
!      Return
!
!    End Subroutine MQC_GHF_Complex_to_LT
!hph-
