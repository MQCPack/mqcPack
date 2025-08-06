!hph #define UCMGAUOPEN
      Module MQC_MatWrapper
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
!     **                          Version 25.7.3                           **
!     **                          July 30, 2025                            **
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
        use mqc_general
        use iso_fortran_env, only: int32, int64, real64
!
!
!
!     Subroutines/Functions...
!
      CONTAINS
      Subroutine MQC_Read_Mat_Precision_Test(MatFileName, IntPrecision)
        Implicit None
!
!     IU receives the Fortran unit number or -1 if the open failed.
!
        Integer( kind=int32 ) LStr, IUUse, Len12D, Len4D
        Parameter (LStr=64,IUUse=57)

        Character(Len=*),Intent(In)::MatFileName
        Integer(kind=int32),Intent(Out)::IntPrecision

        Character(len=LStr)::LabFil, GVers, Title
        Integer(kind=int64) IU
        Integer(kind=int32) IVers,NLab,NAtoms,NBasis,NBsUse,ICharg,Multip,NE,Len12L,Len4L,IOpCl,ICGU

        Integer(kind=int64) IVers_8,NLab_8,NAtoms_8,NBasis_8,NBsUse_8,ICharg_8,Multip_8,NE_8,Len12L_8,Len4L_8,IOpCl_8,ICGU_8
        Integer( kind=int64 ) Len12D_8, Len4D_8
        Integer(kind=int32) IError
        Integer(kind=int64) IError_8
        Integer(kind=int64) IOut
        Integer(kind=int64) rec_len, bufsize
        Integer(kind=int64) LStr64
        Character(len=LStr)::DataType
        logical I_OPENED

!!!!!!!!!!!!!!!!!
        IntPrecision = 8
!!!!!!!!!!!!!!!!!
!
!        Len12D=4
!        Len4D=4
!        IU = IUUse
!!
!        Title = ' '
!        INQUIRE (IU, OPENED=I_OPENED )
!        if ( I_OPENED ) then
!           write(6, 1030) trim(MatFileName)
! 1030      format( "File was Open. Close then open ", a)
!           call Close_MatF(IU)
!        endif
!
!        LStr64 = LStr
!        call mqc_Matfile_Guide_F2C(MatFileName, DataType, LStr64, rec_len, bufsize)
!        if ( trim(DataType) .eq. "Long" ) then
!           IntPrecision = 8
!        else if ( trim(DataType) .eq. "Int" ) then
!           IntPrecision = 4
!        else if ( trim(DataType) .eq. "Raw" ) then
!           IntPrecision = 0
!        else if ( trim(DataType) .eq. "Could not open File" ) then
!           IntPrecision = -1
!        else
!           Iout = 6
!           call MQC_Error_A('Unknown DataType - should be Long, Int or Raw', Iout, 'trim(DataType)', trim(DataType))
!        endIf

        Return
      end subroutine MQC_Read_Mat_Precision_Test

      Subroutine MQC_Mat_Precision(MatFileName, IntPrecision, Call_Test)
        Implicit None
!
!     IU receives the Fortran unit number or -1 if the open failed.
!
        Character(Len=*),Intent(In)::MatFileName
        Integer(kind=int32),Intent(Out)::IntPrecision
        Integer(kind=int32),Intent(In)::Call_Test
        Integer(kind=int32),save::Save_IntPrecision
        Data Save_IntPrecision/0/

        If ( Call_Test .EQ. 0 ) then
           IntPrecision = Save_IntPrecision
        Else If ( Call_Test .EQ. 2 ) then
!
! Write the Matrix Element file.  Set Integer precision to what is was 
! for the Matrix Element file read.  If nothing has been read, set it 
! to 8 Bytes.
           If ( Save_IntPrecision .EQ. 0 ) then
              Save_IntPrecision = 8
           EndIf
           IntPrecision = Save_IntPrecision
        Else
           call MQC_Read_Mat_Precision_Test(MatFileName, IntPrecision)
           Save_IntPrecision = IntPrecision
        EndIf
        return
      end subroutine MQC_Mat_Precision

      Subroutine Wr_LRInd(IU,Label,NR,LenBuf,N1,N2,N3,N4,N5,ASym,X)
        Implicit None

        Character(Len=*),Intent(In)::Label
        Integer(kind=int64),Intent(In)::IU
        Logical,Intent(In)::ASym
        Integer(kind=int64),Intent(In)::NR,LenBuf,N1,N2,N3,N4,N5
        Real(kind=real64),Dimension(:),Intent(InOut)::X

        Integer(kind=int32)::NR_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4
        Integer(kind=int32)::IU_4
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name
        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NR_4 = NR
           LenBuf_4 = LenBuf
           N1_4 = N1
           N2_4 = N2
           N3_4 = N3
           N4_4 = N4
           N5_4 = N5
           call WrLRInd4(IU_4,Label,NR_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4,ASym,X)
        else
! File requires 64-bit integers, so no change
           call WrLRInd8(IU,Label,NR,LenBuf,N1,N2,N3,N4,N5,ASym,X)
        endif
      Return
      End Subroutine Wr_LRInd

      Subroutine Wr_LAO2E(IU,Label,NR,LenBuf,N,RInt)
        Implicit None

        Character(Len=*),Intent(In)::Label
        Integer(kind=int64),Intent(In)::IU
        Integer(kind=int64),Intent(In)::NR,LenBuf,N
        Real(kind=real64),Dimension(:),Intent(InOut)::RInt

        Integer(kind=int32)::NR_4,LenBuf_4,N_4
        Integer(kind=int32)::IU_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name
        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NR_4 = NR
           LenBuf_4 = LenBuf
           N_4 = N
           call WrLAO2E4(IU_4,Label,NR_4,LenBuf_4,N_4,RInt)
        else
! File requires 64-bit integers, so no change
           call WrLAO2E8(IU,Label,NR,LenBuf,N,RInt)
        endif
      Return
      End Subroutine Wr_LAO2E

      Subroutine Open_Read(Name,IU,LabFil,IVers,NLab,GVers,Title,NAtoms,NBasis,NBsUse,ICharg,Multip,NE,Len12L,Len4L,IOpCl,ICGU)
        Implicit None

        Character(Len=*),Intent(In)::Name
        Character(Len=64),Intent(InOut)::LabFil
        Character(Len=64),Intent(InOut)::GVers
        Character(Len=64),Intent(InOut)::Title
        Integer(kind=int64),Intent(Out)::IU,IVers,NLab,NAtoms,NBasis,NBsUse,ICharg,Multip,NE,Len12L,Len4L,IOpCl,ICGU

        Integer(kind=int32)::IU_4,IVers_4,NLab_4,NAtoms_4,NBasis_4,NBsUse_4, &
             ICharg_4,Multip_4,NE_4,Len12L_4,Len4L_4,IOpCl_4,ICGU_4
        Integer(kind=int32)::IntPrecision, Call_Test
        Integer(kind=int64)::IOut
        Call_Test = 1
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)

        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
          call mqc_error('Open_Read: Attempt to use single precision.')
          call OpenRead4(Name,IU_4,LabFil,IVers_4,NLab_4,GVers,Title, &
                NAtoms_4,NBasis_4,NBsUse_4,ICharg_4,Multip_4,NE_4,Len12L_4, &
                Len4L_4,IOpCl_4,ICGU_4)
          IU = IU_4
          IVers = IVers_4
          NLab = NLab_4
          NAtoms = NAtoms_4
          NBasis = NBasis_4
          NBsUse = NBsUse_4
          ICharg = ICharg_4
          Multip = Multip_4
          NE = NE_4
          Len12L = Len12L_4
          Len4L = Len4L_4
          IOpCl = IOpCl_4
          ICGU = ICGU_4
        else if ( IntPrecision .eq. 8 ) then
! File requires 64-bit integers, so no change
           call OpenRead8(Name,IU,LabFil,IVers,NLab,GVers,Title,NAtoms,NBasis,NBsUse,ICharg,Multip,NE,Len12L,Len4L,IOpCl,ICGU)
        else if ( IntPrecision .eq. 0 ) then
! Raw file
           IU = -2
        else if ( IntPrecision .eq. -1 ) then
! Failure to open file
           IU = -1
        else
           IOut = 6
           ICGU = IntPrecision
           call MQC_Error_I('Open_Read: IntPrecision should be 4 or 8', Iout, &
                'IntPrecision', ICGU )
        endif
        if ( IU .eq. -1 ) then
           IOut = 6
           call MQC_Error_A('Open_Read: Unable to Open File', Iout, 'Name', Name)
        endif
      Return
      End Subroutine Open_Read

      Subroutine Open_Write(Name,IU,LabFil,GVers,Title,NAtoms,NBasis,NBsUse,ICharg,Multip,NE,IOpCl,ICGU)
        Implicit None

        Character(Len=*),Intent(In)::Name
        Character(Len=64),Intent(InOut)::LabFil
        Character(Len=64),Intent(InOut)::GVers
        Character(Len=64),Intent(InOut)::Title
        Integer(kind=int64),Intent(Out)::IU
        Integer(kind=int64),Intent(In)::NAtoms,NBasis,NBsUse,ICharg,Multip,NE,IOpCl,ICGU
        Integer(kind=int32)::IU_4,NAtoms_4,NBasis_4,NBsUse_4,ICharg_4,Multip_4,NE_4,IOpCl_4,ICGU_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Integer(kind=int64)::IERR
        Integer(kind=int64)::IOut
        Call_Test = 2

        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           NAtoms_4 = NAtoms
           NBasis_4 = NBasis
           NBsUse_4 = NBsUse
           ICharg_4 = ICharg
           Multip_4 = Multip
           NE_4 = NE
           IOpCl_4 = IOpCl
           ICGU_4 = ICGU
           call OpenWrite4(Name,IU_4,LabFil,GVers,Title,NAtoms_4,NBasis_4,NBsUse_4,ICharg_4,Multip_4,NE_4,IOpCl_4,ICGU_4)
           IU = IU_4
        else if ( IntPrecision .eq. 8 ) then
! File requires 64-bit integers, so no change
           call OpenWrite8(Name,IU,LabFil,GVers,Title,NAtoms,NBasis,NBsUse,ICharg,Multip,NE,IOpCl,ICGU)
        else
           IOut = 6
           IERR = IntPrecision
           call MQC_Error_I('Open_Write: IntPrecision should be 4 or 8', Iout, &
                'IntPrecision', IERR)
        endif
      Return
      End Subroutine Open_Write

      Subroutine Close_MatF(IU)
        Implicit None
!
!     Close a matrix-element file.
!
        Integer(kind=int64),Intent(In)::IU

#ifdef UCMGAUOPEN
        Call CloseMatF8(.false.,IU)
#else
        Call CloseMatF8(IU)
#endif

      Return
      End Subroutine Close_MatF
!hph+
!
!      Function AOInts(CBuf) result(aoint)
!        Implicit None
!        Character(Len=*),Intent(InOut)::CBuf
!        Logical aoint
!        Logical, external::a4ints
!        AOInt = A4Ints(CBuf)
!      End Function AOInts
!
!hph-

      Function LenArr(N1,N2,N3,N4,N5) result(LenAr)
        Implicit None
        Integer(kind=int64),Intent(In)::N1,N2,N3,N4,N5
        Integer(kind=int64)::LenAr

        Integer(kind=int64), external::LenAr8

        LenAr = LenAr8(N1,N2,N3,N4,N5)
      End Function LenArr

!hph+

      Function LInd2C(Check,N1,N2,ASym,I,J,Sign) result(Lind8)
        Implicit None
        Integer(kind=int64),Intent(InOUT)::N1,N2,I,J,Sign
        Logical,Intent(InOut)::Check,ASym

        Integer(kind=int64), external::l8nd2c
        Integer(kind=int64)::Lind8

        Lind8 = L8nd2C(Check,N1,N2,ASym,I,J,Sign)
      End Function LInd2C

      Function Lind2(Check,N1,N2,ASym,I,J,Sign) result(Lind8)
        Implicit None
        Integer(kind=int64),Intent(InOUT)::N1,N2,I,J,Sign
        Logical,Intent(InOut)::Check,ASym

        Integer(kind=int64), external::L8nd2
        Integer(kind=int64)::Lind8

        Lind8 = L8nd2(Check,N1,N2,ASym,I,J,Sign)
      End Function Lind2

      Function Lind3C(Check,N1,N2,N3,ASym,I,J,K,Sign) result(Lind8)
        Implicit None

        Integer(kind=int64),Intent(InOUT)::N1,N2,N3,I,J,K,Sign
        Logical,Intent(InOut)::Check,ASym
        Integer(kind=int64), external::L8nd3C
        Integer(kind=int64)::Lind8

        Lind8 = L8nd3C(Check,N1,N2,N3,ASym,I,J,K,Sign)
      End Function Lind3C

      Function Lind3(Check,N1,N2,N3,ASym,I,J,K,Sign)  result(Lind8)
        Implicit None
        Integer(kind=int64),Intent(InOUT)::N1,N2,N3,I,J,K,Sign
        Logical,Intent(InOut)::Check,ASym

        Integer(kind=int64), external::L8nd3
        Integer(kind=int64)::Lind8

        Lind8 = L8nd3(Check,N1,N2,N3,ASym,I,J,K,Sign)
      End Function Lind3

      Function Lind4C(Check,N1,N2,N3,N4,ASym,I,J,K,L,Sign) result(Lind8)
        Implicit None
        Integer(kind=int64),Intent(InOUT)::N1,N2,N3,N4,I,J,K,L,Sign
        Logical,Intent(InOut)::Check,ASym

        Integer(kind=int64), external::L8nd4C
        Integer(kind=int64)::Lind8

        Lind8 = L8nd4C(Check,N1,N2,N3,N4,ASym,I,J,K,L,Sign)
      End Function Lind4C

      Function Lind4(Check,N1,N2,N3,N4,ASym,I,J,K,L,Sign) result(Lind8)
        Implicit None
        Integer(kind=int64),Intent(InOUT)::N1,N2,N3,N4,I,J,K,L,Sign
        Logical,Intent(InOut)::Check,ASym

        Integer(kind=int64), external::L8nd4
        Integer(kind=int64)::Lind8

        Lind8 = L8nd4(Check,N1,N2,N3,N4,ASym,I,J,K,L,Sign)
      End Function Lind4

      Function Lind5C(Check,N1,N2,N3,N4,N5,ASym,I,J,K,L,M,Sign) result(Lind8)
        Implicit None
        Integer(kind=int64),Intent(InOUT)::N1,N2,N3,N4,N5,I,J,K,L,M,Sign
        Logical,Intent(InOut)::Check,ASym

        Integer(kind=int64),external::L8nd5C
        Integer(kind=int64)::Lind8

        Lind8 = L8nd5C(Check,N1,N2,N3,N4,N5,ASym,I,J,K,L,M,Sign) 
      End Function Lind5C

      Function Lind5(Check,N1,N2,N3,N4,N5,ASym,I,J,K,L,M,Sign) result(Lind8)
        Implicit None
        Integer(kind=int64),Intent(InOUT)::N1,N2,N3,N4,N5,I,J,K,L,M,Sign
        Logical,Intent(InOut)::Check,ASym

        Integer(kind=int64), external::L8nd5
        Integer(kind=int64)::Lind8

        Lind8 = L8nd5(Check,N1,N2,N3,N4,N5,ASym,I,J,K,L,M,Sign)
      End Function Lind5

      Subroutine Rd_2E1(IU,LR,NTot,LenBuf,RArr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,LR,NTOT,LenBuf
        Real(kind=real64),Dimension(:),Intent(InOut)::RArr

        Integer(kind=int32)::IU_4,LR_4,NTOT_4,LenBuf_4
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name
        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           LR_4 = LR
           NTOT_4 = NTOT
           LenBuf_4 = LenBuf
           call Rd2E14(IU_4,LR_4,NTot_4,LenBuf_4,RArr)
        else
! File requires 64-bit integers, so no change
           call Rd2E18(IU,LR,NTot,LenBuf,RArr)
        endif
      Return
      End Subroutine Rd_2E1

      Subroutine Rd_IBuf(IU,NTot,LenBuf,Arr)
        Implicit None
        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Integer(kind=int64),Dimension(NTot),Intent(Out)::Arr

        Integer(kind=int32)::IU_4,NTot_4,LenBuf_4
        Integer(kind=int32),Dimension(NTot)::Arr_4

        Integer(kind=int64)::I
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NTot_4 = NTot
           LenBuf_4 = LenBuf

           Call RdIBuf4(IU_4,NTot_4,LenBuf_4,Arr_4)
           DO I = 1,NTot
              Arr(I) = Arr_4(I)
           ENDDO
        else
! File requires 64-bit integers, so no change
           Call RdIBuf8(IU,NTot,LenBuf,Arr)
        endif
      Return
      End Subroutine Rd_IBuf

      Subroutine Rd_Labl(IU,IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,NRI,EOF)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,IVers
        Character(Len=64),Intent(Out)::CBuf
        Integer(kind=int64),Intent(Out)::NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,NRI
        Logical,Intent(Out)::ASym,EOF

        Integer(kind=int32)::IU_4,IVers_4
        Integer(kind=int32)::NI_4,NR_4,NTot_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4,NRI_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           IVers_4 = IVers
           Call RdLabl4(IU_4,IVers_4,CBuf,NI_4,NR_4,NTot_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4,ASym,NRI_4,EOF)
           NI = NI_4
           NR = NR_4
           NTot = NTot_4
           LenBuf = LenBuf_4
           N1 = N1_4
           N2 = N2_4
           N3 = N3_4
           N4 = N4_4
           N5 = N5_4
           NRI = NRI_4
        else
! File requires 64-bit integers, so no change
           Call RdLabl8(IU,IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,NRI,EOF)
        endif
      Return
      End Subroutine Rd_Labl

      Subroutine Rd_RInd(IU,NR,LR,NTot,LenBuf,LNZ,RArr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NR,LR,NTot,LenBuf
        Integer(kind=int64),Intent(Out)::LNZ
        Real(kind=real64),Dimension(:),Intent(Out)::RArr

        Integer(kind=int32)::IU_4,NR_4,LR_4,NTot_4,LenBuf_4
        Integer(kind=int32)::LNZ_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name
        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NR_4 = NR
           LR_4 = LR
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           Call RdRInd4(IU_4,NR_4,LR_4,NTot_4,LenBuf_4,LNZ_4,RArr)
           LNZ = LNZ_4
        else
! File requires 64-bit integers, so no change
#ifdef UCMGAUOPEN
           Call RdRInd8(IU,NR,LR,NR*LR,NTot,LenBuf,RArr)
#else
           Call RdRInd8(IU,NR,LR,NTot,LenBuf,LNZ,RArr)
#endif
        endif
      Return
      End Subroutine Rd_RInd

      Subroutine Rd_Skip(IU,NTot,LenBuf)
        Implicit None
        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf

        Call RdSkip8(IU,NTot,LenBuf)
      Return
      End Subroutine Rd_Skip

      Subroutine Wr_2E(IU,NTot,NR,N,LR,LenBuf,RArr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,NR,N,LR,LenBuf
        Real(kind=real64),Dimension(:),Intent(Out)::RArr

        Integer(kind=int32)::IU_4,NTot_4,NR_4,N_4,LR_4,LenBuf_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           NTot_4 = NTot
           N_4 = N
           IU_4 = IU
           NR_4 = NR
           LR_4 = LR
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           CALL Wr2E4(IU_4,NTot_4,NR_4,N_4,LR_4,LenBuf_4,RArr)
        else
! File requires 64-bit integers, so no change
           CALL Wr2E8(IU,NTot,NR,N,LR,LenBuf,RArr)
        endif
      Return
      End Subroutine Wr_2E

      Subroutine Wr_IBuf(IU,NTot,LenBuf,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Integer(kind=int64),Dimension(NTot),Intent(In)::Arr

        Integer(kind=int32)::IU_4,NTot_4,LenBuf_4
        Integer(kind=int32),Dimension(NTot)::Arr_4
        Integer(kind=int32)::I

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           DO i = 1, NTot
              Arr_4(I) = Arr(I)
           Enddo
           Call WrIBuf4(IU_4,NTot_4,LenBuf_4,Arr_4)
        else
! File requires 64-bit integers, so no change
           Call WrIBuf8(IU,NTot,LenBuf,Arr)
        endif
      Return
      End Subroutine Wr_IBuf

      Subroutine Wr_Labl(IU,CBuf,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,ASym)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5
        Logical,Intent(In)::ASym
        Character(Len=*),Intent(In)::CBuf

        Integer(kind=int32)::IU_4,NI_4,NR_4,NTot_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NI_4 = NI
           NR_4 = NR
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           N1_4 = N1
           N2_4 = N2
           N3_4 = N3
           N4_4 = N4
           N5_4 = N5
           Call WrLabl4(IU_4,CBuf,NI_4,NR_4,NTot_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4,ASym)
        else
! File requires 64-bit integers, so no change
           Call WrLabl8(IU,CBuf,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,ASym)
        endif
      Return
      End Subroutine Wr_Labl

      Subroutine Wr_RBuf(IU,NTot,LenBuf,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Real(kind=real64),Dimension(:),Intent(In)::Arr

        Integer(kind=int32)::IU_4,NTot_4,LenBuf_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           Call WrRBuf4(IU_4,NTot_4,LenBuf_4,Arr)
        else
! File requires 64-bit integers, so no change
           Call WrRBuf8(IU,NTot,LenBuf,Arr)
        endif
      Return
      End Subroutine Wr_RBuf

      Subroutine Wr_RInd(IU,NR,LR,NTot,LenBuf,RArr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NR,LR,NTot,LenBuf
        Real(kind=real64),Dimension(:),Intent(In)::RArr

        Integer(kind=int32)::IU_4,NR_4,LR_4,NTot_4,LenBuf_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           NR_4 = NR
           LR_4 = LR
           Call WrRInd4(IU_4,NR_4,LR_4,NTot_4,LenBuf_4,RArr)
        else
! File requires 64-bit integers, so no change
           Call WrRInd8(IU,NR,LR,NTot,LenBuf,RArr)
        endif
      Return
      End Subroutine Wr_RInd

      Subroutine ExpAO1(N,LR,RI,RO)
        Implicit None

        Integer(kind=int64),Intent(In)::N, LR
        Real(kind=real64),Dimension(:),Intent(Out)::RI
        Real(kind=real64),Dimension(:),Intent(Out)::RO
! This routine does no I/O
        Call ExpA18(N,LR,RI,RO)
      Return
      End Subroutine ExpAO1

      Subroutine ExpAON(NE,N,LR,RI,RO)
        Implicit None

        Integer(kind=int64),Intent(In)::NE,N, LR
        Real(kind=real64),Dimension(:),Intent(Out)::RI
        Real(kind=real64),Dimension(:),Intent(Out)::RO
! This routine does no I/O
        Call ExpAN8(NE,N,LR,RI,RO)
      Return
      End Subroutine ExpAON

      Subroutine AClear(N,A)
        Implicit None

        Integer(kind=int64),Intent(In)::N
        Real(kind=real64),Dimension(:),Intent(Out)::A

        call AClea8(N,A)
      Return
      End Subroutine AClear

      Subroutine IClear(N,IA)
        Implicit None

        Integer(kind=int64),Intent(In)::N
        Integer(kind=int64),Dimension(N),Intent(Out)::IA

        Call IClea8(N,IA)
      Return
      End Subroutine IClear

!hph-

!hph+
      Integer(kind=int64) Function mqc_LenArr(NR,NI,NRI,N1,N2,N3,N4,N5, &
           s_datatype, s_shape, s_formula )
      Implicit None
!
! This is a modified Lenarr routine from GauOpen  with these changes:
!    1. Added code to capture the math, so that it's easier to figure out 
!       the dimensions of arrays.
!    2. Converted from FORTRAN 77 to FORTRAN 2003
!
      Integer(kind=int64),intent(in)::NR,NI,NRI
      Integer(kind=int64),intent(in)::N1,N2,N3,N4,N5
      character(LEN=*),intent(out)::s_datatype, s_shape, s_formula
      Integer(kind=int64) N1X,N2X,N3X,N4X,N5X,Sign
      Integer(kind=int64) size
      character(LEN=1024) s_N1X, s_N2X, s_N3X, s_N4X, s_N5X
      character(LEN=1024) s_absN1X,s_absN2X,s_absN3X,s_absN4X,s_absN5X
      character(LEN=1024) s_TMP, s_Lind5_result, s_Logic
      Logical::ASym
      character(len=8)::tab
!
!     Format statements.
!
 1030 Format( "  - 2-D MATRIX allocated full (",a," x ",a,")" )
 1050 Format( "  - SYMMETRIC MATRIX ((",a,"+1)*",a,"/2)" )
 1060 Format( "  - VECTOR(",a,")" )
 1080 Format( "  - MULTIPLE SYMMETIC MATRICES (((",a,"+1)x",a,"/2) x ",a,")" )
 1090 Format( "  - MULTIPLE SYMMETIC MATRICES (((",a,"+1)x",a,"/2) x (",a,"x",a,"))" )
 1110 Format( "  - SYMMETRIC SUPER MATRIX (((",a,"+1) x ",a,")/2+1) x ((",a,"+1) x ",a,")/4" )
 1120 Format( "  - TRISYMMETRIC MATRIX ((",a,"+2) x (",a,"+1) x ",a,"1)/6" )
 1130 Format( "  - 3-D MATRIX allocated full (",a," x ",a," x ",a,")" )
 2010 Format( i15 )
!
 1231 format( a)
 1232 format( 2a)
 1233 format( 3a)

      s_datatype= "UNKNOWN"
      if(NR.gt.0.and.NI.gt.0) then
        if((NR.eq.1.or.NR.eq.2.or.NR.eq.3).and.NI.eq.4) then
           s_datatype= "TENSOR"
        else
           s_datatype= "MIXED"
        endIf
      elseIf(NI.eq.1) then
         s_datatype= "INTEGER"
      elseIf(NR.eq.1) then
        if(NRI.eq.1) then
           s_datatype= "REAL"
        elseIf(NRI.eq.2) then
           s_datatype= "COMPLEX"
        endIf
      endIf
      if (N1.gt.1.and.N2.eq.1.and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
         s_shape="VECTOR"
         write( s_N1X,2010) N1
         call mqc_StripSpaces(s_N1X)
         write( s_formula, 1060) trim(s_N1X)
         mqc_LenArr = N1
         return
      else If(N1.gt.1.and.N2.gt.1.and.(N3.eq.1.or.N3.eq.0) &
           .and.(N4.eq.1.or.N4.eq.1).and.(N5.eq.1.or.N5.eq.0)) then
         s_shape="MATRIX"
         write( s_N1X,2010) N1
         call mqc_StripSpaces(s_N1X)
         write( s_N2X,2010) N2
         call mqc_StripSpaces(s_N2X)
         write( s_formula, 1030) trim(s_N1X),trim(s_N2X)
         mqc_LenArr = N1*N2
         return
      elseIf(N1.le.-1.and.N2.eq.(-N1).and.N3.eq.1.and.N4.eq.1.and.N5.eq.1) then
         s_shape="SYMMATRIX"
         write( s_N2X,2010) N2
         call mqc_StripSpaces(s_N2X)
         write( s_formula, 1050) trim(s_N2X),trim(s_N2X)
         mqc_LenArr = ((N2+1)*N2)/2
         return
      elseIf(N1.le.-1.and.N2.eq.N1.and.N3.eq.N1.and.N4.eq.(-N1).and.N5.eq.1) then
         s_shape="SYMSYM"
         Size = ((N4+1)*N4)/2
         mqc_LenArr = ((Size+1)*Size)/2
         write( s_N4X,2010) N4
         call mqc_StripSpaces(s_N4X)
         write( s_formula, 1110)trim(s_N4X),trim(s_N4X),trim(s_N4X),trim(s_N4X)
         return
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.gt.1.and.N4.eq.1.and.N5.eq.1) then
         s_shape="MULTISYMMATRICES"
         mqc_LenArr = ((N2+1)*N2)/2*N3
         write( s_N2X,2010) N2
         call mqc_StripSpaces(s_N2X)
         write( s_N3X,2010) N3
         call mqc_StripSpaces(s_N3X)
         write( s_formula, 1080)trim(s_N2X),trim(s_N2X),trim(s_N3X)
         return
      elseIf(N1.le.-1.and.N2.gt.1.and.N3.gt.1.and.N4.gt.1.and.N5.eq.1) then
         s_shape="MULTISYMMATRICES"
         mqc_LenArr = ((N2+1)*N2)/2*N3*N4
         write( s_N2X,2010) N2
         call mqc_StripSpaces(s_N2X)
         write( s_N3X,2010) N3
         call mqc_StripSpaces(s_N3X)
         write( s_N4X,2010) N4
         call mqc_StripSpaces(s_N4X)
         write( s_formula, 1090)trim(s_N2X),trim(s_N2X),trim(s_N3X),trim(s_N4X)
         return
      elseIf(N1.le.-1.and.N2.eq.N1.and.N3.eq.(-N1).and.N4.eq.1.and.N5.eq.1) then
         s_shape="TRISYMMATRICES"
         mqc_LenArr = ((N3+2)*(N3+1)*N3)/6
         write( s_N3X,2010) N3
         call mqc_StripSpaces(s_N3X)
         write( s_formula, 1120)trim(s_N3X),trim(s_N3X),trim(s_N3X)
         return
      else If((N1.eq.1.or.N1.eq.0).and.(N2.eq.1.or.N2.eq.0).and. &
           (N3.eq.1.or.N3.eq.0).and.(N4.eq.1.or.N4.eq.0).and.(N5.eq.1 &
           .or.N5.eq.0)) then
         s_shape="SCALAR"
         mqc_LenArr = 1
         s_formula = " 1 "
         return
      else If((N4.eq.1.or.N4.eq.0).and.(N5.eq.1.or.N5.eq.0).and.N1.gt.1 &
           .and.N2.gt.1.and.N3.gt.1) then
         s_shape="3DMATRIX"
         mqc_LenArr = N1*N2*N3
         write( s_N1X,2010) N1
         call mqc_StripSpaces(s_N1X)
         write( s_N2X,2010) N2
         call mqc_StripSpaces(s_N2X)
         write( s_N3X,2010) N3
         call mqc_StripSpaces(s_N3X)
         write( s_formula, 1130)trim(s_N1X),trim(s_N2X),trim(s_N3X)
         return
      endif
      s_shape="UNKNOWN"
      N1X = N1
      If(N1X.eq.0) N1X = 1
      N2X = N2
      If(N2X.eq.0) N2X = 1
      N3X = N3
      If(N3X.eq.0) N3X = 1
      N4X = N4
      If(N4X.eq.0) N4X = 1
      N5X = N5
      If(N5X.eq.0) N5X = 1

      s_Logic = " "
      if ( N1X .eq. 1 ) Then
         s_N1X = "1"
         s_absN1X = "1"
         s_Logic = "(NA.eq.1.or.NA.eq.0)"
      else if ( N1X .eq. -1 ) Then
         s_N1X = "-1"
         s_absN1X = "1"
         s_Logic = "NA.eq.-1"
      else
         s_N1X = "NA"
         if (  N1X .gt. 0 ) Then
            s_absN1X = s_N1X
         else
            s_absN1X = "(-NA)"
         endif
      endif
      if ( N2X .eq. 1 ) Then
         s_N2X = "1"
         s_absN2X = "1"
         write( s_TMP,1232) trim(s_Logic),".and.(NB.eq.1.or.NB.eq.0)"
         s_Logic=s_TMP
      else if ( N2X .eq. -1 ) Then
         s_N2X = "-1"
         s_absN2X = "1"
         write( s_TMP,1232) trim(s_Logic),".and.NB.eq.-1"
         s_Logic=s_TMP
      else if ( N1.eq.N2 ) then
         s_N2X = s_N1X
         s_absN2X = s_absN1X
         write( s_TMP,1232) trim(s_Logic),".and.NA.eq.NB"
         s_Logic=s_TMP
      else if ( (-N1).eq.N2 ) then
         s_N2X = "(-NA)"
         s_absN2X = s_absN1X
         write( s_TMP,1232) trim(s_Logic),".and.(-NA).eq.NB"
         s_Logic=s_TMP
      else
         s_N2X = "NB"
         if (  N2X .gt. 0 ) Then
            s_absN2X = s_N2X
         else
            s_absN2X = "(-NB)"
         endif
      endif
      if ( N3X .eq. 1 ) Then
         s_N3X = "1"
         s_absN3X = "1"
         write( s_TMP,1232) trim(s_Logic),".and.(NC.eq.1.or.NC.eq.0)"
         s_Logic=s_TMP
      else if ( N3X .eq. -1 ) Then
         s_N3X = "-1"
         s_absN3X = "1"
         write( s_TMP,1232) trim(s_Logic),".and.NC.eq.-1"
         s_Logic=s_TMP
      else if ( N1.eq.N3 ) then
         s_N3X = s_N1X
         s_absN3X = s_absN1X
         write( s_TMP,1232) trim(s_Logic),".and.NA.eq.NC"
         s_Logic=s_TMP
      else if ( (-N1).eq.N3 ) then
         s_N3X = "(-NA)"
         s_absN3X = s_absN1X
         write( s_TMP,1232) trim(s_Logic),".and.(-NA).eq.NC"
         s_Logic=s_TMP
      else if ( N2.eq.N3 ) then
         s_N3X = s_N2X
         s_absN3X = s_absN2X
         write( s_TMP,1232) trim(s_Logic),".and.NB.eq.NC"
         s_Logic=s_TMP
      else
         s_N3X = "NC"
         if (  N3X .gt. 0 ) Then
            s_absN3X = s_N3X
         else
            s_absN3X = "(-NC)"
         endif
      endif
      if ( N4X .eq. 1 ) Then
         s_N4X = "1"
         s_absN4X = "1"
         write( s_TMP,1232) trim(s_Logic),".and.(ND.eq.1.or.ND.eq.0)"
         s_Logic=s_TMP
      else if ( N4X .eq. -1 ) Then
         s_N4X = "-1"
         s_absN4X = "1"
         write( s_TMP,1232) trim(s_Logic),".and.ND.eq.-1"
         s_Logic=s_TMP
      else if ( N1.eq.N4 ) then
         s_N4X = s_N1X
         s_absN4X = s_absN1X
         write( s_TMP,1232) trim(s_Logic),".and.NA.eq.ND"
         s_Logic=s_TMP
      else if ( (-N1).eq.N4 ) then
         s_N4X = "(-NA)"
         s_absN4X = s_absN1X
         write( s_TMP,1232) trim(s_Logic),".and.(-NA).eq.ND"
         s_Logic=s_TMP
      else if ( N2.eq.N4 ) then
         s_N4X = s_N2X
         s_absN4X = s_absN2X
         write( s_TMP,1232) trim(s_Logic),".and.NB.eq.ND"
         s_Logic=s_TMP
      else if ( N3.eq.N4 ) then
         s_N4X = s_N3X
         s_absN4X = s_absN3X
         write( s_TMP,1232) trim(s_Logic),".and.NC.eq.ND"
         s_Logic=s_TMP
      else
         s_N4X = "ND"
         if (  N4X .gt. 0 ) Then
            s_absN4X = s_N4X
         else
            s_absN4X = "(-ND)"
         endif
      endif
      if ( N5X .eq. 1 ) Then
         s_N5X = "1"
         s_absN5X = "1"
         write( s_TMP,1232) trim(s_Logic),".and.(NE.eq.1.or.NE.eq.0)"
         s_Logic=s_TMP
      else if ( N5X .eq. -1 ) Then
         s_N5X = "-1"
         s_absN5X = "1"
         write( s_TMP,1232) trim(s_Logic),".and.NE.eq.-1"
         s_Logic=s_TMP
      else if ( N1.eq.N5 ) then
         s_N5X = s_N1X
         s_absN5X = s_absN1X
         write( s_TMP,1232) trim(s_Logic),".and.NA.eq.NE"
         s_Logic=s_TMP
      else if ( N2.eq.N5 ) then
         s_N5X = s_N2X
         s_absN5X = s_absN2X
         write( s_TMP,1232) trim(s_Logic),".and.NB.eq.NE"
         s_Logic=s_TMP
      else if ( N3.eq.N5 ) then
         s_N5X = s_N3X
         s_absN5X = s_absN3X
         write( s_TMP,1232) trim(s_Logic),".and.NC.eq.NE"
         s_Logic=s_TMP
      else if ( N4.eq.N5 ) then
         s_N5X = s_N4X
         s_absN5X = s_absN4X
         write( s_TMP,1232) trim(s_Logic),".and.ND.eq.NE"
         s_Logic=s_TMP
      else
         s_N5X = "NE"
         if ( N5X .gt. 0 ) Then
            s_absN5X = s_N5X
         else
            s_absN5X = "(-NE)"
         endif
      endif

      ASym=.False.
      mqc_LenArr = mqc_Lind5(.False.,N1X,N2X,N3X,N4X,N5X,ASym,Abs(N1X), &
           Abs(N2X),Abs(N3X),Abs(N4X),abs(N5X),Sign, &
           s_N1X,s_N2X,s_N3X,s_N4X,s_N5X, &
           s_absN1X,s_absN2X,s_absN3X,s_absN4X,s_absN5X, &
           s_Lind5_result,s_Logic)+1
      call mqc_simplify_formula_F2C(s_Lind5_result)
      write( s_formula, 1233) trim(s_Lind5_result),"+1  Here's the logic:",trim(s_Logic)

      Return
      End Function mqc_Lenarr
!hph-

!hph+

      Integer(kind=int64) Function mqc_Lind2(Check,N1,N2,ASym,I,J,Sign, &
           s_N1,s_N2,s_I,s_J,s_Lind2_result, s_Logic)

      Implicit None
!
!     Linear or square indexing, I,J are 1-based,
!     output is 0-based.  Sign is +/-1.
!
      Logical,intent(in)::Check,ASym
      Integer(kind=int64),intent(out)::Sign
      Integer(kind=int64),intent(in)::N1,N2,I,J
      character(LEN=*),intent(in)::s_N1,s_N2,s_I,s_J
      character(LEN=*),intent(out)::s_Lind2_result, s_Logic
      character(LEN=1024) s_TMP
         
      Sign = 1
      If(Check.and.(N2.le.0.or.N1.eq.0.or.(N1.lt.0.and.N1.ne.(-N2)).or. &
           I.lt.1.or.I.gt.Abs(N1).or.J.lt.1.or.J.gt.N2)) then
        mqc_Lind2 = -1
        s_Lind2_result = "-1"
        Return
        endIf
 1231   format(a)
 1232   format(2a)
 1233   format(3a)
 1234   format(4a)
 1235   format(5a)
 1236   format(6a)
 1237   format(7a)
      If(N1.lt.0) then
         write(s_TMP,1234) trim(s_Logic),".and.",trim(s_N1),".lt.0"
         s_Logic = s_TMP
        If(I.ge.J) then
           write(s_TMP,1235) trim(s_Logic),".and.",trim(s_I),".ge.",trim(s_J)
           s_Logic = s_TMP
           mqc_Lind2 = (I*(I-1))/2 + J - 1
           write( s_Lind2_result, 1237)"(",trim(s_I),"*(", &
                trim(s_I),"-1))/2 + ",trim(s_J), "-1"
        else
           write(s_TMP,1235) trim(s_Logic),".and.",trim(s_I),".lt.",trim(s_J)
           s_Logic = s_TMP
           mqc_Lind2 = (J*(J-1))/2 + I - 1
           write( s_Lind2_result, 1237) "(",trim(s_J),"*(", &
                trim(s_J),"-1))/2+",trim(s_I),"-1"
          If(ASym) Sign = -1
          endIf
      else
         If(N1.gt.1) then
            write(s_TMP,1234) trim(s_Logic),".and.",trim(s_N1),".gt.1"
         else
            write(s_TMP,1234) trim(s_Logic),".and.",trim(s_N1),".ge.0"
         endif
         s_Logic = s_TMP
         mqc_Lind2 = N1*(J-1) + I - 1
         write( s_Lind2_result, 1236)trim(s_N1),"*(", &
              trim(s_J),"-1)+",trim(s_I),"-1"
        endIf
      Return
      End Function mqc_Lind2

      Integer(kind=int64) Function mqc_Lind3(Check,N1,N2,N3,ASym,I,J,K,Sign, &
           s_N1,s_N2,s_N3,s_I,s_J,s_K,s_Lind3_result,s_Logic)
      Implicit None
!
!     Linear or square indexing, I,J,K are 1-based,
!     output is 0-based.  Sign is +/-1.
!
      Logical,intent(in)::Check
      Logical,intent(inout)::ASym
      Integer(kind=int64),intent(in)::N1,N2,N3,I,J,K
      Integer(kind=int64),intent(out)::Sign
      character(LEN=*),intent(in)::s_N1,s_N2,s_N3,s_I,s_J,s_K
      character(LEN=*),intent(out)::s_Lind3_result, s_Logic
      Integer(kind=int64) I1,J1,K1,N12,IJ
      character(LEN=1024) s_ij,s_i1,s_j1,s_k1,s_N12,s_Lind2_result
 1231   format(a)
 1232   format(2a)
 1233   format(3a)
 1234   format(4a)
 1235   format(5a)
 1236   format(6a)
 1237   format(7a)
 1238   format(8a)
 1239   format(9a)
 1240   format(10a)
 1242   format(12a)
      Sign = 1
      If(Check.and.(N3.le.0.or.(N1*N2).eq.0.or. &
           (N1.lt.0.and.N1.ne.(-Abs(N2))).or. &
           (N2.lt.0.and.N2.ne.(-N3)).or.I.lt.1.or.I.gt.Abs(N1).or. &
           J.lt.1.or.J.gt.Abs(N2).or.K.lt.1.or.K.gt.N3)) then
        mqc_Lind3 = -1
        s_Lind3_result = "-1"
        Return
        endIf
      I1 = I - 1
      if ( I1 .ne. 0 ) Then
         write( s_I1, 1233)"(",trim(s_I),"-1)"
      else
         s_I1 = "0"
      endif
      J1 = J - 1
      if ( J1 .ne. 0 ) Then
         write( s_J1, 1233)"(",trim(s_J),"-1)"
      else
         s_J1 = "0"
      endif

      K1 = K - 1
      if ( K1 .ne. 0 ) Then
         write( s_K1, 1233)"(",trim(s_K),"-1)"
      else
         s_K1 = "0"
      endif

      If(N1.gt.0) then
        If(N2.gt.0) then
          mqc_LInd3 = N1*(N2*K1+J1) + I1
          write( s_Lind3_result, 1239)trim(s_N1),"*(",trim(s_N2), &
               "*",trim(s_K1),"+",trim(s_J1),")+",trim(s_I1)
        else
          mqc_Lind3 = N1*mqc_Lind2(.False.,N2,N3,ASym,J,K,Sign, &
               s_N2,s_N3,s_J,s_K,s_Lind2_result,s_Logic)
          write( s_Lind3_result, 1234)trim(s_N1),"*(", &
               trim(s_Lind2_result),")"
          endIf
      else if(N2.gt.0) then
        N12 = (N2*(N2+1))/2
        write( s_N12, 1235)"(",trim(s_N2),"*(", &
             trim(s_N2),"+1))/2"
        If(I.ge.J) then
          IJ = (I*I1)/2 + J1
          write( s_IJ, 1236)"(",Trim(s_I),"*", &
               Trim(s_I1),")/2 + ",Trim(s_J1)
        else
          IJ = (J*J1)/2 + I1
          write( s_IJ, 1236)"(",Trim(s_J),"*", &
               Trim(s_J1),")/2+",Trim(s_I1)
          If(ASym) Sign = -1
          endIf
        mqc_Lind3 = (N12)*K1 + IJ
        write( s_Lind3_result, 1236)"(",trim(s_N12),")*", &
             trim(s_K1),"+",trim(s_IJ)
      else
        K1 = Max(I,J,K) - 1
        if ( I .eq. J ) then
           if ( I .eq. K ) then
              write( s_K1, 1233)"(",Trim(s_I),"-1)"
           Else
              write( s_K1, 1235)"(Max(",Trim(s_I),",", &
                   Trim(s_K),")-1)"
           EndIf
        Else
           write( s_K1, 1237)"(Max(",Trim(s_I),",", &
                Trim(s_J),",",Trim(s_K),")-1)"
        EndIf
        I1 = Min(I,J,K) - 1
        if ( I .eq. J ) then
           if ( I .eq. K ) then
              write( s_K1, 1233)"(",Trim(s_I),"-1)"
           Else
              write( s_I1, 1235)"(Min(",Trim(s_I),",", &
                   Trim(s_K),")-1)"
           EndIf
        Else
           write( s_I1, 1237)"(Min(",Trim(s_I),",", &
                Trim(s_J),",",Trim(s_K),")-1)"
        EndIf
        J1 = I + J + K - K1 - I1 - 3
        write( s_J1, 1240)Trim(s_I),"+",Trim(s_J),"+", &
             Trim(s_K),"-",Trim(s_K1),"-",Trim(s_I1),"-3"
        mqc_LInd3 = I1 + ((J1)*(J1+1))/2 + ((K1)*(K1+1)*(K1+2))/6
        write( s_Lind3_result, 1242)Trim(s_I1),"+((", &
             Trim(s_J1),")*(",Trim(s_J1),"+1))/2+((", &
             Trim(s_K1),")*(",Trim(s_K1),"+1)*(",Trim(s_K1),"+2))/6"

        endIf
      Return
      End Function mqc_Lind3

      Integer(kind=int64) Function mqc_Lind4(Check,N1,N2,N3,N4,ASym,I,J,K,L,Sign, &
           s_N1,s_N2,s_N3,s_N4,s_I,s_J,s_K,s_L, &
           s_Lind4_result,s_Logic)
      Implicit None
!
!     Linear or square indexing, I,J,K,L are 1-based,
!     output is 0-based.  Sign is +/-1.
!
      Logical,intent(in)::Check
      Logical,intent(inout)::ASym
      character(LEN=*),intent(in)::s_N1,s_N2,s_N3,s_N4,s_I,s_J,s_K,s_L
      character(LEN=*),intent(out)::s_Lind4_result,s_Logic
      Integer(kind=int64),intent(in)::N1,N2,N3,N4,I,J,K,L
      Integer(kind=int64),intent(inout)::Sign
      character(LEN=1024) s_I1,s_J1,s_K1,s_L1,s_ICase
      character(LEN=1024) s_TMP1, s_TMP2, s_JK, s_JKL, s_KL
      character(LEN=1024) s_ij, s_N12, s_n23, s_IJK, s_N123
      character(LEN=1024) s_Lind2_result, s_Lind3_result
      Integer(kind=int64) I1,J1,K1,L1,ICase,N23 
      Integer(kind=int64) N123,KL,IJ,IJK,JKL,JK,N12,SignIJ,SignKL
!
      If(Check.and.(N4.le.0.or.(N1*N2*N3).eq.0.or. &
           (N1.lt.0.and.N1.ne.(-Abs(N2))).or. &
           (N2.lt.0.and.N2.ne.(-Abs(N3))).or. &
           (N3.lt.0.and.N3.ne.(-N4)).or. &
           I.lt.1.or.I.gt.Abs(N1).or.J.lt.1.or.J.gt.Abs(N2).or. &
           K.lt.1.or.K.gt.Abs(N3).or.L.lt.1.or.L.gt.N4)) then
        mqc_Lind4 = -1
        Sign = 1
        s_Lind4_result = "-1"
        Return
        endIf
 1231   format(a)
 1232   format(2a)
 1233   format(3a)
 1235   format(5a)
 1236   format(6a)
 1237   format(7a)
 1238   format(8a)
 1239   format(9a)
 1243   format(13a)
 1244   format(14a)

      I1 = I - 1
      if ( I1.ne.0 ) then
         write(S_I1,1233)"(",trim(s_i),"-1)"
      else
         S_I1 = "0"
      endif
      J1 = J - 1
      if ( J1.ne.0 ) then
         write(S_J1,1233)"(",trim(s_j),"-1)"
      else
         S_J1 = "0"
      endif
      K1 = K - 1
      if ( K1.ne.0 ) then
         write(S_K1,1233)"(",trim(s_k),"-1)"
      else
         S_K1 = "0"
      endif
      L1 = L - 1
      if ( L1.ne.0 ) then
         write(S_L1,1233)"(",trim(s_L),"-1)"
      else
         S_L1 = "0"
      endif

      ICase = 0
      If(N1.lt.0 ) ICase = ICase + 1
      If(N2.lt.0) ICase = ICase + 2
      If(N3.lt.0) ICase = ICase + 4
      if ( ICase .eq. 0 ) then
!
!     No symmetries.
!
         If(N1.gt.1 ) then
            write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N1),".gt.1"
         else
            write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N1),".ge.0"
         endif
         s_Logic=s_TMP1
         If(N2.gt.1 ) then
            write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N2),".gt.1"
         else
            write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N2),".ge.0"
         endif
         s_Logic=s_TMP1
         If(N3.gt.1 ) then
            write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N3),".gt.1"
         else
            write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N3),".ge.0"
         endif
         s_Logic=s_TMP1

         mqc_Lind4 = N1*(N2*(N3*L1+K1)+J1) + I1
         write(s_Lind4_result,1243)trim(s_N1),"*(", &
              trim(s_N2),"*(", trim(s_N3),"*",trim(s_L1), &
              "+",trim(s_K1),")+",trim(s_J1),")+",trim(s_I1)
         Sign = 1
         Return
      else if ( ICase .eq. 1 ) then
!
!     I<=J
!
         write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N1),".lt.0.and.", &
              trim(s_N2),".ge.0.and.",trim(s_N3),".ge.0"
         s_Logic=s_TMP1

         IJ = mqc_Lind2(.False.,N1,N2,ASym,I,J,Sign, &
              s_N1,s_N2,s_I,s_J,s_Lind2_result,s_Logic)
         s_ij = s_Lind2_result
         N12 = (N2*(N2+1))/2
         write(s_n12,1237)"(",Trim(s_N2),"*(",trim(s_N2),"+1))/2"

         mqc_Lind4 = N12*(N3*L1+K1) + IJ
         write(s_Lind4_result,1239)trim(s_N12),"*(",Trim(s_N3),"*", &
              Trim(s_L1),"+",Trim(s_K1),")+",trim(s_IJ)

         Return
      else if ( ICase .eq. 2 ) then
!
!     I,J<=K,L
!
         write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N1),".ge.0.and.", &
              trim(s_N2),".lt.0.and.",trim(s_N3),".ge.0"
         s_Logic=s_TMP1

         JK = mqc_Lind2(.False.,N2,N3,ASym,J,K,Sign, &
              s_N2,s_N3,s_J,s_K,s_Lind2_result,s_Logic)
         s_JK = s_Lind2_result

         N23 = (N3*(N3+1))/2
         write(s_n23,1235)"(",Trim(s_N3),"*(", Trim(s_N3),"+1))/2"
         mqc_Lind4 = N1*(N23*L1+JK) + I1

         write(s_Lind4_result,1239)trim(s_N1),"*(",trim(s_N23), &
              "*",trim(s_L1),"+",trim(s_JK),")+",trim(s_I1)
         Return
      else if ( ICase .eq. 3 ) then
!
!     I<=J<=K,L
!
         write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N1),".lt.0.and.", &
              trim(s_N2),".lt.0.and.",trim(s_N3),".ge.0"

         s_Logic=s_TMP1

         IJK = mqc_Lind3(.False.,N1,N2,N3,ASym,I,J,K,Sign, &
              s_N1,s_N2,s_N3,s_I,s_J,s_K,s_Lind3_result,s_Logic)
         s_IJK = s_Lind3_result

         N123 = (N3*(N3+1)*(N3+2))/6
         write(s_N123,1237)"(",Trim(s_N3),"*(",Trim(s_N3),"+1)*(", &
              Trim(s_N3),"+2))/6"

         mqc_Lind4 = (N123)*L1 + IJK
         if ( L1 .eq. 0 .or. N123 .eq. 0 ) then
            s_Lind4_result = Trim(s_IJK)
         Else
            write(s_Lind4_result,1236)"(",Trim(s_N123),")*", &
                 Trim(s_L1),"+",Trim(s_IJK)
         EndIf
         Return
      else if ( ICase .eq. 4 ) then
!
!     I,J,K<=L
!
         write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N1),".ge.0.and.", &
              trim(s_N2),".ge.0.and.",trim(s_N3),".lt.0"
         s_Logic=s_TMP1

         KL = mqc_Lind2(.False.,N3,N4,ASym,K,L,Sign, &
              s_N3,s_N4,s_K,s_L,s_Lind2_result,s_Logic)
         s_KL = s_Lind2_result

         mqc_Lind4 = N1*(N2*(KL)+J1) + I1
         write(s_Lind4_result,1239)trim(s_N1),"*(",trim(s_N2), &
              "*(",trim(s_KL),")+",trim(s_J1),")+",trim(s_I1)
         Return
      else if ( ICase .eq. 5 ) then
!
!     I<=J,K<=L
!
         write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N1),".lt.0.and.", &
              trim(s_N2),".ge.0.and.",trim(s_N3),".lt.0"
         s_Logic=s_TMP1

         IJ = mqc_Lind2(.False.,N1,N2,ASym,I,J,SignIJ, &
              s_N1,s_N2,s_I,s_J,s_Lind2_result,s_Logic)
         s_IJ = s_Lind2_result

         KL = mqc_Lind2(.False.,N3,N4,ASym,K,L,SignKL, &
              s_N3,s_N4,s_K,s_L,s_Lind2_result,s_Logic)
         s_KL = s_Lind2_result

         N12 = (N2*(N2+1))/2
         write(s_N12,1235)"(",trim(s_N2),"*(",trim(s_N2),"+1))/2"

         mqc_Lind4 = N12*(KL) + IJ
         write(s_Lind4_result,1235)trim(s_N12),"*(",trim(s_KL),")+", &
           trim(s_IJ)
         Sign = SignIJ*SignKL
         Return
      else if ( ICase .eq. 6 ) then
!
!     I,J<=K<=L
!
         write(s_TMP1,1238)trim(s_Logic),".and.",trim(s_N1),".ge.0.and.", &
              trim(s_N2),".lt.0.and.",trim(s_N3),".lt.0"
         s_Logic=s_TMP1

         JKL = mqc_Lind3(.False.,N2,N3,N4,ASym,J,K,L,Sign, &
              s_N2,s_N3,s_N4,s_J,s_K,s_L,s_Lind3_result,s_Logic)
         s_JKL = s_Lind3_result

         mqc_Lind4 = N1*(JKL) + I1
         If ( N1.eq.0 .or. JKL.eq.0) then
            s_Lind4_result = s_I1
         Else
            write(s_Lind4_result,1235)trim(s_N1),"*(", &
                 trim(s_JKL),")+",trim(s_I1)
         EndIF
         Return
      else if ( ICase .eq. 7 ) then
!
!     I<=J<=K<=L
!
         write(s_TMP1,1238) trim(s_Logic),".and.",trim(s_N1),".lt.0.and.", &
              trim(s_N2),".lt.0.and.",trim(s_N3),".lt.0"
         s_Logic=s_TMP1

         IJ = mqc_Lind2(.False.,N1,N4,ASym,I,J,SignIJ, &
              s_N1,s_N4,s_I,s_J,s_Lind2_result,s_Logic)
         s_IJ = s_Lind2_result

         KL = mqc_Lind2(.False.,N3,N4,ASym,K,L,SignKL, &
              s_N3,s_N4,s_K,s_L,s_Lind2_result,s_Logic)
         s_kl = s_Lind2_result

         write(s_TMP1,1233)"(",trim(s_IJ),"+1)"
         write(s_TMP2,1233)"(",trim(s_KL),"+1)"
         mqc_Lind4 = mqc_Lind2(.False.,N1,N4,ASym,IJ+1,KL+1,Sign, &
              s_N1,s_N4,s_TMP1,s_TMP2,s_Lind2_result,s_Logic)
         s_Lind4_result = s_Lind2_result

         Sign = Sign*SignIJ*SignKL
      endif
      Return
      End Function mqc_Lind4

      Integer(kind=int64) Function mqc_Lind5(Check,N1,N2,N3,N4,N5,ASym,I,J,K,L, &
           M,Sign,s_N1,s_N2,s_N3,s_N4,s_N5,s_I,s_J,s_K,s_L,s_M, &
           s_Lind5_result,s_Logic)
      Implicit None
!
!     Linear or square indexing, I,J,K,L,M are 1-based, output is
!     0-based.  Sign is +/-1, LM indices can not be lower triangular.
!
      Logical,intent(in)::Check
      Logical,intent(inout)::ASym
      Integer(kind=int64),intent(in)::N1,N2,N3,N4,N5,I,J,K,L,M
      Integer(kind=int64),intent(inout)::Sign
      Integer(kind=int64) M1,N1A,N2A,N3A,N1234
      character(LEN=*),intent(in)::s_N1,s_N2,s_N3,s_N4,s_N5,s_I
      character(LEN=*),intent(in)::s_J,s_K,s_L,s_M
      character(LEN=*),intent(out)::s_Lind5_result,s_Logic
      character(LEN=1024) s_M1, s_N1A, s_N2A, s_N3A
      character(LEN=1024) s_N1234, s_Lind4_result
!
      If(Check.and.(N5.le.0.or.N4.le.0.or.(N1*N2*N3).eq.0.or. &
           (N1.lt.0.and.N1.ne.(-Abs(N2))).or. &
           (N2.lt.0.and.N2.ne.(-Abs(N3))).or. &
           (N3.lt.0.and.N2.ne.(-N4)).or. &
           I.lt.1.or.I.gt.Abs(N1).or.J.lt.1.or.J.gt.Abs(N2).or. &
           K.lt.1.or.K.gt.Abs(N3).or.L.lt.1.or.L.gt.N4)) then
        mqc_Lind5 = -1
        Sign = 1
        s_Lind5_result = "-1"
        Return
        endIf
 1231 format(a)
 1232 format(2a)
 1233 format(3a)
 1235 format(5a)
 1236 format(6a)
      M1 = M - 1
      if ( M1.ne.0 ) then
         write( s_M1, 1233) "(",trim(s_M),"-1)"
      else 
         s_M1 = "0"
      endif

      N1A = abs(N1)
      if ( N1.lt.0 ) then
         write( s_N1A, 1233) "(-(",trim(s_N1),"))"
      else 
         s_N1A = s_N1
      endif

      N2A = abs(N2)
      if ( N2.lt.0 ) then
         write( s_N2A, 1233) "(-(",trim(s_N2),"))"
      else
         s_N2A = s_N2
      endif

      N3A = abs(N3)
      if ( N3.lt.0 ) then
         write( s_N3A, 1233) "(-(",trim(s_N3),"))"
      else
         s_N3A = s_N3
      endif

      N1234 = mqc_Lind4(.False.,N1,N2,N3,N4,ASym,N1A,N2A,N3A,N4,Sign, &
           s_N1,s_N2,s_N3,s_N4,s_N1A,s_N2A,s_N3A, &
           s_N4,s_Lind4_result,s_Logic) + 1
      write( s_N1234, 1233) "(",trim(s_Lind4_result),"+1)"

      mqc_Lind5 = mqc_Lind4(.False.,N1,N2,N3,N4,ASym,I,J,K,L,Sign, &
           s_N1,s_N2,s_N3,s_N4,s_I,s_J,s_K,s_L, &
           s_Lind4_result,s_Logic ) + N1234*M1

      if ( M1.eq.0 .or. N1234.eq.0 ) then
         s_Lind5_result = s_Lind4_result
      Else
         write( s_Lind5_result, 1235) trim(s_Lind4_result),"+", &
              trim(s_N1234),"*",trim(s_M1)
      EndIF

      Return
      End Function mqc_Lind5

!hph-

      subroutine mqc_StripSpaces(string)
        character(len=*) :: string
        integer :: stringLen
        integer :: last, actual

        stringLen = len (string)
        last = 1
        actual = 1

        do while (actual < stringLen)
           if (string(last:last) == ' ') then
              actual = actual + 1
              string(last:last) = string(actual:actual)
              string(actual:actual) = ' '
           else
              last = last + 1
              if (actual < last) &
                   actual = last
           endif
        end do
      end subroutine mqc_StripSpaces

!hph+

      Subroutine Rd_RMixed(IU,NTot,LenBuf,IArr,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Integer(kind=int64),Dimension(:),Intent(Out)::IArr
        Real(kind=real64),Dimension(:),Intent(Out)::Arr

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           Call RdRMixed4(IU,NTot,LenBuf,IArr,Arr)
        else
! File requires 64-bit integers, so no change
           Call RdRMixed8(IU,NTot,LenBuf,IArr,Arr)
        endif
      Return
      End Subroutine Rd_RMixed

      subroutine RdRMixed8(IU,NTot,LenBuf,IArr,Arr)
      implicit none
      Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
      Integer(kind=int64),Dimension(:),Intent(Out)::IArr
      Real(kind=real64),Dimension(:),Intent(Out)::Arr

      Integer(kind=int64),Dimension(lenbuf)::iBuf
      Real(kind=real64),Dimension(lenbuf)::Buf

      integer ind,ndo,i
!
      do ind = 0, (ntot-1), lenbuf
         ndo = min(lenbuf,ntot-ind)
         read(iu) iBuf, buf
         do i = 1, ndo
            iarr(ind+i) = ibuf(i)
            arr(ind+i) = buf(i)
         enddo
      enddo
      return
      end subroutine RdRMixed8

      subroutine RdRMixed4(IU,NTot,LenBuf,IArr,Arr)
      implicit none
      Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
      Integer(kind=int64),Dimension(:),Intent(Out)::IArr
      Real(kind=real64),Dimension(:),Intent(Out)::Arr

      Integer(kind=int32),Dimension(lenbuf)::iBuf
      Real(kind=real64),Dimension(lenbuf)::Buf

      integer ind,ndo,i
!
      do ind = 0, (ntot-1), lenbuf
         ndo = min(lenbuf,ntot-ind)
         read(iu) iBuf, buf
         do i = 1, ndo
            iarr(ind+i) = ibuf(i)
            arr(ind+i) = buf(i)
         enddo
      enddo
      return
      end subroutine RdRMixed4

      Subroutine Wr_RMixed(IU,NTot,LenBuf,IArr,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Integer(kind=int64),Dimension(:),Intent(Out)::IArr
        Real(kind=real64),Dimension(:),Intent(Out)::Arr

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           Call WrRMixed4(IU,NTot,LenBuf,IArr,Arr)
        else
! File requires 64-bit integers, so no change
           Call WrRMixed8(IU,NTot,LenBuf,IArr,Arr)
        endif
      Return
      End Subroutine Wr_RMixed

      subroutine WrRMixed8(IU,NTot,LenBuf,IArr,Arr)
      implicit none
      Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
      Integer(kind=int64),Dimension(:),Intent(Out)::IArr
      Real(kind=real64),Dimension(:),Intent(Out)::Arr

      Integer(kind=int64),Dimension(lenbuf)::iBuf
      Real(kind=real64),Dimension(lenbuf)::Buf

      integer ind,ndo,i
!
      do ind = 0, (ntot-1), lenbuf
         ndo = min(lenbuf,ntot-ind)
         do i = 1, ndo
            ibuf(i) = iarr(ind+i)
            buf(i) = arr(ind+i)
         enddo
         write(iu) iBuf, buf
      enddo
      return
      end subroutine WrRMixed8

      subroutine WrRMixed4(IU,NTot,LenBuf,IArr,Arr)
      implicit none
      Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
      Integer(kind=int64),Dimension(:),Intent(Out)::IArr
      Real(kind=real64),Dimension(:),Intent(Out)::Arr

      Integer(kind=int32),Dimension(lenbuf)::iBuf
      Real(kind=real64),Dimension(lenbuf)::Buf

      integer ind,ndo,i
!
      do ind = 0, (ntot-1), lenbuf
         ndo = min(lenbuf,ntot-ind)
         do i = 1, ndo
            ibuf(i) = iarr(ind+i)
            buf(i) = arr(ind+i)
         enddo
         write(iu) iBuf, buf
      enddo
      return
      end subroutine WrRMixed4

      Subroutine Rd_CMixed(IU,NTot,LenBuf,IArr,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Integer(kind=int64),Dimension(:),Intent(Out)::IArr
        Complex(kind=real64),Dimension(:),Intent(Out)::Arr

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           Call RdCMixed4(IU,NTot,LenBuf,IArr,Arr)
        else
! File requires 64-bit integers, so no change
           Call RdCMixed8(IU,NTot,LenBuf,IArr,Arr)
        endif
      Return
      End Subroutine Rd_CMixed

      subroutine RdCMixed8(IU,NTot,LenBuf,IArr,Arr)
      implicit none
      Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
      Integer(kind=int64),Dimension(:),Intent(Out)::IArr
      Complex(kind=real64),Dimension(:),Intent(Out)::Arr

      Integer(kind=int64),Dimension(lenbuf)::iBuf
      Complex(kind=real64),Dimension(lenbuf)::Buf

      integer ind,ndo,i
!
      do ind = 0, (ntot-1), lenbuf
         ndo = min(lenbuf,ntot-ind)
         read(iu) iBuf, buf
         do i = 1, ndo
            iarr(ind+i) = ibuf(i)
            arr(ind+i) = buf(i)
         enddo
      enddo
      return
      end subroutine RdCMixed8

      subroutine RdCMixed4(IU,NTot,LenBuf,IArr,Arr)
      implicit none
      Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
      Integer(kind=int64),Dimension(NTot),Intent(Out)::IArr
      Complex(kind=real64),Dimension(NTot),Intent(Out)::Arr

      Integer(kind=int32),Dimension(lenbuf)::iBuf
      Complex(kind=real64),Dimension(lenbuf)::Buf

      integer ind,ndo,i
!
      do  ind = 0, (ntot-1), lenbuf
         ndo = min(lenbuf,ntot-ind)
         read(iu) iBuf, buf
         do  i = 1, ndo
            iarr(ind+i) = ibuf(i)
            arr(ind+i) = buf(i)
         enddo
      enddo
      return
      end subroutine RdCMixed4

      Subroutine Wr_CMixed(IU,NTot,LenBuf,IArr,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Integer(kind=int64),Dimension(:),Intent(Out)::IArr
        Complex(kind=real64),Dimension(:),Intent(Out)::Arr

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           Call WrCMixed4(IU,NTot,LenBuf,IArr,Arr)
        else
! File requires 64-bit integers, so no change
           Call WrCMixed8(IU,NTot,LenBuf,IArr,Arr)
        endif
      Return
      End Subroutine Wr_CMixed

      subroutine WrCMixed8(IU,NTot,LenBuf,IArr,Arr)
      implicit none
      Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
      Integer(kind=int64),Dimension(:),Intent(Out)::IArr
      Complex(kind=real64),Dimension(:),Intent(Out)::Arr

      Integer(kind=int64),Dimension(lenbuf)::iBuf
      Complex(kind=real64),Dimension(lenbuf)::Buf

      integer ind,ndo,i
!
      do ind = 0, (ntot-1), lenbuf
         ndo = min(lenbuf,ntot-ind)
         do i = 1, ndo
            ibuf(i) = iarr(ind+i)
            buf(i) = arr(ind+i)
         enddo
         write(iu) iBuf, buf
      enddo
      return
      end subroutine WrCMixed8

      subroutine WrCMixed4(IU,NTot,LenBuf,IArr,Arr)
      implicit none
      Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
      Integer(kind=int64),Dimension(NTot),Intent(Out)::IArr
      Complex(kind=real64),Dimension(NTot),Intent(Out)::Arr

      Integer(kind=int32),Dimension(lenbuf)::iBuf
      Complex(kind=real64),Dimension(lenbuf)::Buf

      integer ind,ndo,i
!
      do  ind = 0, (ntot-1), lenbuf
         ndo = min(lenbuf,ntot-ind)
         do  i = 1, ndo
            ibuf(i) = iarr(ind+i)
            buf(i) = arr(ind+i)
         enddo
         write(iu) iBuf, buf
      enddo
      return
      end subroutine WrCMixed4

      Subroutine Rd_2EN(IU,NR,LR,LRNR,NTot,LenBuf,RArr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NR,LR,LRNR,NTot,LenBuf
        Real(kind=real64),Dimension(:),Intent(InOut)::RArr

        Integer(kind=int32)::IU_4,NR_4,LR_4,LRNR_4,NTot_4,LenBuf_4
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NR_4 = NR
           LR_4 = LR
           LRNR_4 = LRNR
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           call Rd2EN4(IU,NR,LR,LRNR,NTot,LenBuf,RArr)
        else
! File requires 64-bit integers, so no change
           call Rd2EN8(IU,NR,LR,LRNR,NTot,LenBuf,RArr)
        endif
      Return
      End Subroutine Rd_2EN

      Subroutine Rd_RBuf(IU,NTot,LenBuf,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Real(kind=real64),Dimension(:),Intent(Out)::Arr

        Integer(kind=int32)::IU_4,NTot_4,LenBuf_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           Call RdRBuf4(IU_4,NTot_4,LenBuf_4,Arr)
        else
! File requires 64-bit integers, so no change
           Call RdRBuf8(IU,NTot,LenBuf,Arr)
        endif
      Return
      End Subroutine Rd_RBuf

      Subroutine Wr_LRBuf(IU,Label,NR,LenBuf,N1,N2,N3,N4,N5,ASym,X)
        Implicit None

        Character(Len=*),Intent(In)::Label
        Integer(kind=int64),Intent(In)::IU
        Logical,Intent(In)::ASym
        Integer(kind=int64),Intent(In)::NR,LenBuf,N1,N2,N3,N4,N5
        Real(kind=real64),Dimension(:),Intent(InOut)::X

        Integer(kind=int32)::NR_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4
        Integer(kind=int32)::IU_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name
        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NR_4 = NR
           LenBuf_4 = LenBuf
           N1_4 = N1
           N2_4 = N2
           N3_4 = N3
           N4_4 = N4
           N5_4 = N5
           call WrLRBuf4(IU_4,Label,NR_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4,ASym,X)
        else
! File requires 64-bit integers, so no change
           call WrLRBuf8(IU,Label,NR,LenBuf,N1,N2,N3,N4,N5,ASym,X)
        endif
      Return
      End Subroutine Wr_LRBuf

      Subroutine Wr_Head(IU,NAtoms,NAt3,NBasis,IAn,IAtTyp,AtmChg,C,IBfAtm, &
           IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NAtoms,NAt3,NBasis
        INTEGER(kind=int64),Dimension(NAtoms),Intent(In)::IAn
        INTEGER(kind=int64),Dimension(NAtoms),Intent(In)::IAtTyp
        REAL(kind=real64),Dimension(:),Intent(In)::AtmChg
        REAL(kind=real64),Dimension(:),Intent(In)::C
        INTEGER(kind=int64),Dimension(NBasis),Intent(In)::IBfTyp
        INTEGER(kind=int64),Dimension(NBasis),Intent(In)::IBfAtm
        REAL(kind=real64),Dimension(:),Intent(In)::AtmWgt
        Integer(kind=int64),Intent(In)::NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot

        Integer(kind=int32)::IU_4,NAtoms_4,NAt3_4,NBasis_4
        INTEGER(kind=int32),Dimension(NAtoms)::IAn_4
        INTEGER(kind=int32),Dimension(NAtoms)::IAtTyp_4
        INTEGER(kind=int32),Dimension(NBasis)::IBfTyp_4
        INTEGER(kind=int32),Dimension(NBasis)::IBfAtm_4
        Integer(kind=int32)::NFC_4,NFV_4,ITran_4,IDum9_4,NShlAO_4,NPrmAO_4,NShlDB_4,NPrmDB_4,NBTot_4

        Integer(kind=int64) I

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NAtoms_4 = NAtoms
           NAt3_4 = NAt3
           NBasis_4 = NBasis
           DO i = 1, NAtoms
              IAn_4(I) = IAn(I)
              IAtTyp_4(I) = IAtTyp(I)
           Enddo
           DO i = 1, NBasis
              IBfTyp_4(I) = IBfTyp(I)
              IBfAtm_4(I) = IBfAtm(I)
           Enddo
           NFC_4 = NFC
           NFV_4 = NFV
           ITran_4 = ITran
           IDum9_4 = IDum9
           NShlAO_4 = NShlAO
           NPrmAO_4 = NPrmAO
           NShlDB_4 = NShlDB
           NPrmDB_4 = NPrmDB
           NBTot_4 = NBTot
           Call WrHead4(IU_4,NAtoms_4,NAt3_4,NBasis_4,IAn_4,IAtTyp_4,AtmChg,C,IBfAtm_4, &
           IBfTyp_4,AtmWgt,NFC_4,NFV_4,ITran_4,IDum9_4,NShlAO_4,NPrmAO_4,NShlDB_4,NPrmDB_4,NBTot_4)

        else
! File requires 64-bit integers, so no change
           Call WrHead8(IU,NAtoms,NAt3,NBasis,IAn,IAtTyp,AtmChg,C,IBfAtm, &
                IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot)
        endif
      Return
      End Subroutine Wr_Head

      Subroutine Rd_Head(IU,NLab,NAtoms,NBasis,IAn,IAtTyp,AtmChg,C, &
       IBfAtm,IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB, &
       NPrmDB,NBTot)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NAtoms,NBasis,NLab
        Integer(kind=int64),Dimension(NAtoms),Intent(Out)::IAn
        Integer(kind=int64),Dimension(NAtoms),Intent(Out)::IAtTyp
        Real(kind=real64),Dimension(:),Intent(Out)::AtmChg
        Real(kind=real64),Dimension(:),Intent(InOut)::C
        Integer(kind=int64),Dimension(NBasis),Intent(Out)::IBfAtm
        Integer(kind=int64),Dimension(NBasis),Intent(Out)::IBfTyp
        Real(kind=real64),Dimension(:),Intent(InOut)::AtmWgt
        Integer(kind=int64),Intent(Out)::NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot

        Integer(kind=int32)::IU_4,NAtoms_4,NBasis_4,NLab_4
        Integer(kind=int32),Dimension(NAtoms)::IAn_4
        Integer(kind=int32),Dimension(NAtoms)::IAtTyp_4
        Integer(kind=int32),Dimension(NBasis)::IBfAtm_4
        Integer(kind=int32),Dimension(NBasis)::IBfTyp_4
        Integer(kind=int32)::NFC_4,NFV_4,ITran_4,IDum9_4,NShlAO_4,NPrmAO_4,NShlDB_4,NPrmDB_4,NBTot_4
        Integer(kind=int64)::I
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NAtoms_4 = NAtoms
           NBasis_4 = NBasis
           NLab_4 = NLab
           Call RdHead4(IU_4,NLab_4,NAtoms_4,NBasis_4,IAn_4,IAtTyp_4,AtmChg,C, &
                IBfAtm_4,IBfTyp_4,AtmWgt,NFC_4,NFV_4,ITran_4,IDum9_4,NShlAO_4,NPrmAO_4,NShlDB_4, &
                NPrmDB_4,NBTot_4)
           DO I = 1,NAtoms
              IAn(I) = IAn_4(I)
              IAtTyp(I) = IAtTyp_4(I)
           ENDDO
           DO I = 1,NBasis
              IBfAtm(I) = IBfAtm_4(I)
              IBfTyp(I) = IBfTyp_4(I)
           ENDDO
           NFC = NFC_4
           NFV = NFV_4
           ITran = ITran_4
           IDum9 = IDum9_4
           NShlAO = NShlAO_4
           NPrmAO = NPrmAO_4
           NShlDB = NShlDB_4
           NPrmDB = NPrmDB_4
           NBTot = NBTot_4
        else
! File requires 64-bit integers, so no change
           Call RdHead8(IU,NLab,NAtoms,NBasis,IAn,IAtTyp,AtmChg,C, &
                IBfAtm,IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB, &
                NPrmDB,NBTot)
        Endif
      Return
      End Subroutine Rd_Head

      Subroutine Wr_LIBuf(IU,Label,NI,LenBuf,N1,N2,N3,N4,N5,ASym,IX)
        Implicit None
!
        Character(Len=*),Intent(In)::Label
        Integer(kind=int64),Intent(In)::IU
        Logical,Intent(In)::ASym
        Integer(kind=int64),Intent(In)::NI,LenBuf,N1,N2,N3,N4,N5
        Integer(kind=int64),Dimension(:),Intent(In)::IX

        Integer(kind=int32) IU_4,NI_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4

        Integer(kind=int32),Dimension(:),Allocatable::IX_4
        Integer(kind=int32)::I
        Integer(kind=int64)::NTot
!
! This routine writes and integer array, so it needs to be copied.
!
        Integer(kind=int32)::AllocateStatus
        Integer(kind=int32)::DeAllocateStatus
        Integer(kind=int64), external::LenAr8
!
        Integer(kind=int64)::Iout
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name
        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NI_4 = NI
           LenBuf_4 = LenBuf
           N1_4 = N1
           N2_4 = N2
           N3_4 = N3
           N4_4 = N4
           N5_4 = N5
           NTot = LenAr8(N1,N2,N3,N4,N5)
           Iout = 6
           Allocate ( IX_4(NTot), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) call MQC_Error_I('Not Enough Memory', Iout, 'NTot', NTot)
           do i=1, NTot
              IX_4(i) = IX(i)
           enddo
           call WrLIBuf4(IU_4,Label,NI_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4,ASym,IX_4)
           DEALLOCATE (IX_4, STAT = DeAllocateStatus)
           IF (DeAllocateStatus /= 0) call MQC_Error_I('Not able to deallocate Memory', Iout, 'NTot', NTot)
        else
! File requires 64-bit integers, so no change
           call WrLIBuf8(IU,Label,NI,LenBuf,N1,N2,N3,N4,N5,ASym,IX)
        endif
      Return
      End Subroutine Wr_LIBuf

      Subroutine Wr_LCBuf(IU,Label,NR,LenBuf,N1,N2,N3,N4,N5,ASym,X)
        Implicit None

        Character(Len=*),Intent(In)::Label
        Integer(kind=int64),Intent(In)::IU
        Logical,Intent(In)::ASym
        Integer(kind=int64),Intent(In)::NR,LenBuf,N1,N2,N3,N4,N5
        Complex(kind=real64),Dimension(:),Intent(InOut)::X

        Integer(kind=int32)::NR_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4
        Integer(kind=int32)::IU_4
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name
        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NR_4 = NR
           LenBuf_4 = LenBuf
           N1_4 = N1
           N2_4 = N2
           N3_4 = N3
           N4_4 = N4
           N5_4 = N5
           call WrLCBuf4(IU_4,Label,NR_4,LenBuf_4,N1_4,N2_4,N3_4,N4_4,N5_4,ASym,X)
        else
! File requires 64-bit integers, so no change
           call WrLCBuf8(IU,Label,NR,LenBuf,N1,N2,N3,N4,N5,ASym,X)
        endif
      Return
      End Subroutine Wr_LCBuf

      Subroutine Rd_CBuf(IU,NTot,LenBuf,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Complex(kind=real64),Dimension(NTot),Intent(InOut)::Arr

        Integer(kind=int32)::IU_4,NTot_4,LenBuf_4
        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           call RdCBuf4(IU_4,NTot_4,LenBuf_4,Arr)
        else
! File requires 64-bit integers, so no change
           call RdCBuf8(IU,NTot,LenBuf,Arr)
        endif
      Return
      End Subroutine Rd_CBuf

      Subroutine Wr_CBuf(IU,NTot,LenBuf,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Complex(kind=real64),Dimension(NTot),Intent(In)::Arr

        Integer(kind=int32)::IU_4,NTot_4,LenBuf_4

        Integer(kind=int32)::IntPrecision, Call_Test
        Character(Len=64)::Name

        Call_Test = 0
        call MQC_Mat_Precision(Name, IntPrecision, Call_Test)
        if ( IntPrecision .eq. 4 ) then
! File requires 32-bit integers
           IU_4 = IU
           NTot_4 = NTot
           LenBuf_4 = LenBuf
           call WrCBuf4(IU_4,NTot_4,LenBuf_4,Arr)
        else
! File requires 64-bit integers, so no change
           call WrCBuf8(IU,NTot,LenBuf,Arr)
        endif
      Return
      End Subroutine Wr_CBuf
!hph-

      End Module MQC_MatWrapper
