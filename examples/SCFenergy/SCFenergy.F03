! This program tests Reading everthing in a MartixFile
     program SCFenergy
!
!     This program contains examples of the matrix file read and write routines
!
!     Variable Declarations...       
!
      Use MQC_Algebra
      Use MQC_DataStructures
      Use MQC_Files
      Use MQC_EST
      use mqc_General
      use mqc_FullWavefunction
      use mqc_matwrapper
      use iso_fortran_env
      implicit none

      real(kind=real64)::Vnn,EE,EN
      real(kind=real64)::ENL
      real(kind=real64)::SCFe
      Complex(kind=real64)::EE_C,EN_C
      integer(kind=int64) iout
      integer(kind=int64) I, status
      integer(kind=int64) Print_Options
      integer(kind=int32) N, J, K
      class(mqc_link_list), allocatable:: mqc_input_file_list
      class(mqc_READONLY_list), allocatable:: my_READONLY_list
      class(MQC_ALL_Wavefunction),pointer::current_wavefunction
      type(mqc_array_link_list_node),pointer::DenA
      type(mqc_array_link_list_node),pointer::DenB
      type(mqc_array_link_list_node),pointer::HamA
      type(mqc_array_link_list_node),pointer::HamB
      type(mqc_array_link_list_node),pointer::FockA
      type(mqc_array_link_list_node),pointer::FockB
      integer(kind=int64) Found_DenA,Found_DenB,Found_HamA,Found_HamB
      integer(kind=int64) Found_FockA,Found_FockB
      character(len=8)::tab
!
!     Format statements.
!
 1010 Format( A )
 1020 Format( "Matrix File:  ",A )
 1200 format( "   Electronic Energy:",a,a,E20.10,/, &
           "   Electronic/Nuclear Energy:",a,E20.9,/, &
           "   Nuclear Energy:",a,a,E20.10,/, &
           "   Total Energy:",a,a,E20.10)
!
!     Do the work...
!
!
! This program has problems.  It does not work properly for DFT or GHF
!
      IOut = 6
      write( IOut, 1010) "Warning, SCFenergy gives incorrect energies for DFT and GHF wavefunctions"

!
!     Create a list of data structures to try to read.
!
      allocate(my_READONLY_list)
! This data is needed for Nuclear Replusion
      call mqc_READONLY_list_add(my_READONLY_list,"CARTESIAN NUCLEAR COORDINATES(BOHR)")
      call mqc_READONLY_list_add(my_READONLY_list,"ATOM CHARGES")
! This data is needed for components that include electronic interactions
      call mqc_READONLY_list_add(my_READONLY_list,"ALPHA SCF DENSITY MATRIX")
      call mqc_READONLY_list_add(my_READONLY_list,"ALPHA DENSITY MATRIX")
      call mqc_READONLY_list_add(my_READONLY_list,"BETA SCF DENSITY MATRIX")
      call mqc_READONLY_list_add(my_READONLY_list,"BETA DENSITY MATRIX")
      call mqc_READONLY_list_add(my_READONLY_list,"CORE HAMILTONIAN ALPHA")
      call mqc_READONLY_list_add(my_READONLY_list,"CORE HAMILTONIAN BETA")
      call mqc_READONLY_list_add(my_READONLY_list,"ALPHA FOCK MATRIX")
      call mqc_READONLY_list_add(my_READONLY_list,"BETA FOCK MATRIX")

      tab = char( 9 )
      allocate(mqc_input_file_list)
      Print_Options=0

      Call MQC_Create_Wavefunctions( Print_Options, my_READONLY_list, &
           mqc_input_file_list, Iout)
!
!     Create a list of MatrixFiles.  Use MatrixFiles supplied.
!     If Gaussian input has been supplied:
!     1. Modify the input include MatrixFile options
!     2. Run the job to generate MatrixFiles
!     If the name of the Gaussian binary is not given, and Gaussian is not
!     in the path, this will fail.
!
#ifdef MATRIXFILE
      call mqc_link_list_start(mqc_input_file_list,Status,current_wavefunction)
      Do While( Status.eq.1 )
         write(iout,1020) trim(current_wavefunction%FileName)
! Get Nuclear repulsion Energy
         Vnn = MQC_Wave_To_Nuclear_Repulsion(current_wavefunction,IOut)
! Get all possible components of Electronic part of energy
         Found_DenA=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
              "ALPHA SCF DENSITY MATRIX", DenA, iout)
         If ( Found_DenA .eq. 0 ) then
            Found_DenA=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
                 "ALPHA DENSITY MATRIX", DenA, iout)
         EndIf
         Found_DenB=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
              "BETA SCF DENSITY MATRIX", DenB, iout)
         If ( Found_DenB .eq. 0 ) then
            Found_DenB=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
                 "BETA DENSITY MATRIX", DenB, iout)
         EndIf
         Found_HamA=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
              "CORE HAMILTONIAN ALPHA", HamA, iout)
         Found_HamB=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
              "CORE HAMILTONIAN BETA", HamB, iout)
         Found_FockA=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
              "ALPHA FOCK MATRIX", FockA, iout)
         Found_FockB=MQC_Find_mqc_array_link_list_node(current_wavefunction,&
              "BETA FOCK MATRIX", FockB, iout)

         If ( Found_DenB .eq. 2 .and. Found_HamB .eq. 2 .and. &
              Found_FockB .eq. 2 ) then
!
! The easy way to do a Triangular Matrix/Matrix multiply is to double the 
! off diagonal elements of 1 array and then do a dot product.  
! It's computationally a bit faster to half the diagonal elements of the 
! density Matrix and adjust result result as needed.
!
            N = (current_wavefunction%nbasis+1)* &
                 current_wavefunction%nbasis/2
            J = 0
            DO K = 1, current_wavefunction%nbasis
               J = J+K
               DenA%ArrR(J,1,1,1,1) = DenA%ArrR(J,1,1,1,1)*0.5D0
               DenB%ArrR(J,1,1,1,1) = DenB%ArrR(J,1,1,1,1)*0.5D0
            EndDo
            EE = 0.0
            EN = 0.0
            DO K = 1, N
               EE = EE + DenA%ArrR(K,1,1,1,1) * HamA%ArrR(K,1,1,1,1) &
                    + DenB%ArrR(K,1,1,1,1) * HamB%ArrR(K,1,1,1,1)

               EN = EN + DenA%ArrR(K,1,1,1,1) * FockA%ArrR(K,1,1,1,1) &
                    + DenB%ArrR(K,1,1,1,1) * FockB%ArrR(K,1,1,1,1)
            enddo
            SCFe = EE + EN + Vnn
         Else If ( Found_DenA .eq. 2 .and. Found_HamA .eq. 2 .and. &
              Found_FockA .eq. 2 ) then
!
! Half Diagonal elements of Density
!
            N = (current_wavefunction%nbasis+1)* &
                 current_wavefunction%nbasis/2
            J = 0
            DO K = 1, current_wavefunction%nbasis
               J = J+K
               DenA%ArrR(J,1,1,1,1) = DenA%ArrR(J,1,1,1,1)*0.5D0
            EndDo
            EE = 0.0
            EN = 0.0
            DO K = 1, N
               EE = EE + DenA%ArrR(K,1,1,1,1) * HamA%ArrR(K,1,1,1,1)
               EN = EN + DenA%ArrR(K,1,1,1,1) * FockA%ArrR(K,1,1,1,1)
            enddo
            EE = EE*2.0D0
            EN = EN*2.0D0
            SCFe = EE + EN + Vnn
         Else If ( Found_DenB .eq. 3 .and. Found_HamB .eq. 3 .and. &
              Found_FockB .eq. 3 ) then
!
! Half Diagonal elements of Density
!
! Array sizes for this seem to be based on 2*NBasis, used sizes stored with arrays
            N = DenA%LR
            J = 0
            DO K = 1, DenA%N2
               J = J+K
               DenA%ArrC(J,1,1,1,1) = DenA%ArrC(J,1,1,1,1)*0.5D0
               DenB%ArrC(J,1,1,1,1) = DenB%ArrC(J,1,1,1,1)*0.5D0
            EndDo
            EE_C = 0.0
            EN_C = 0.0
            DO K = 1, N
               EE_C = EE_C + DenA%ArrC(K,1,1,1,1) * HamA%ArrC(K,1,1,1,1) &
                    + DenB%ArrC(K,1,1,1,1) * HamB%ArrC(K,1,1,1,1)
               EN_C = EN_C + DenA%ArrC(K,1,1,1,1) * FockA%ArrC(K,1,1,1,1) &
                    + DenB%ArrC(K,1,1,1,1) * FockB%ArrC(K,1,1,1,1)
            enddo
            EE = DBLE(EE_C)
            EN = DBLE(EN_C)
            SCFe = EE + EN + Vnn
         Else If ( Found_DenA .eq. 3 .and. Found_HamA .eq. 3 .and. &
              Found_FockA .eq. 3 ) then
!
! Half Diagonal elements of Density
!
! Array sizes for this seem to be based on 2*NBasis, used sizes stored with arrays
            N = DenA%LR
            J = 0
            DO K = 1, DenA%N2
               J = J+K
               DenA%ArrC(J,1,1,1,1) = DenA%ArrC(J,1,1,1,1)*0.5D0
            EndDo
            EE_C = 0.0
            EN_C = 0.0
            DO K = 1, N
               EE_C = EE_C + DenA%ArrC(K,1,1,1,1) * HamA%ArrC(K,1,1,1,1)
               EN_C = EN_C + DenA%ArrC(K,1,1,1,1) * FockA%ArrC(K,1,1,1,1)
            enddo
            EE = DBLE(EE_C)*2.0D0
            EN = DBLE(EN_C)*2.0D0
            SCFe = EE + EN + Vnn
         Else
            call mqc_error_I('SCFenergy: Do not have a case for these flags &
                 & for components of Electrical components of Energy', &
                 IOut,&
                 'Found_DenA', Found_DenA, 'Found_DenB', Found_DenB, &
                 'Found_HamA', Found_HamA, 'Found_HamB', Found_HamB, &
                 'Found_FockA', Found_FockA, 'Found_FockB', Found_FockB )
         EndIf
         write(6, 1200)tab, tab, EE, tab, EN, tab, tab, Vnn, tab, tab, SCFe
         call mqc_link_list_next(mqc_input_file_list,Status, &
              current_wavefunction)
      EndDo
#endif
    end program SCFenergy

