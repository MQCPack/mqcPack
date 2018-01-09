      Module MQC_MatWrapper
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
!     **    Dave Mullally                                                 **
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
!                                                                           
!                                                                           
!                                                                           
!     This module includes procedures that provide a variety of general purpose
!     utilities. The subroutines and functions provided by this module are 
!     grouped into the following sections:  
!
! C to Fortran 2003 Interface
!
        use mqc_general
        use iso_fortran_env
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
        logical I_OPENED

        Len12D=4
        Len4D=4
        IU = IUUse
!
        Title = ' '
        INQUIRE (IU, OPENED=I_OPENED )
        if ( I_OPENED ) then
           write(6, 1030) trim(MatFileName)
 1030      format( "File was Open. Close then open ", a)
           call Close_MatF(IU)
        endif
! There can be trouble reading files that were just written.  The file system
! may need to be flushed

        Open (Unit=IU,File=TRIM(MatFileName),Form='Unformatted',Status='Old', IOStat=IError)
        If (IError.ne.0) then
           IOut = 6
           call MQC_Error_A('Error Opening File', Iout, 'MatFileName', MatFileName)
        endif

        Read(IU) LabFil(1:LStr), IVers, NLab, GVers(1:LStr)
        Read(IU,err=2000,end=2000) Title(1:LStr), NAtoms, NBasis, NBsUse, ICharg, Multip, NE, Len12L, Len4L
2000    continue
        If(Len4L.ne.Len4D.or.Len12L.ne.Len12D) then
           IntPrecision = 8
#ifdef DEBUG
           Write(6,1000) Len12D, Len4D, Len12L, Len4L, IntPrecision
1000       Format(' Read_Precision_Test compiled with Len12=',I1,' Len4=',I1, &
              ' but file has Len12=',I1,' Len4=', I1,' Precision=',I1,'.')
#endif
           call Close_MatF(IU)

           Len12D_8=8
           Len4D_8=8
           Open (Unit=IUUse,File=TRIM(MatFileName),Form='Unformatted',Status='Old', IOStat=IError)
           Read(IU) LabFil(1:LStr), IVers_8, NLab_8, GVers(1:LStr)
           Read(IU) Title(1:LStr), NAtoms_8, NBasis_8, NBsUse_8, ICharg_8, Multip_8, NE_8, Len12L_8, Len4L_8
           If(Len4L_8.ne.Len4D_8.or.Len12L_8.ne.Len12D_8) then
              Write(6,1080) trim(MatFileName)
1080          Format(' File:', A)
              IOut = 6
              call MQC_Error_I('Unable to make sense of File', Iout, 'Len4L_8', &
                   Len4L_8,'Len4D_8',Len4D_8,'Len12L_8',Len12L_8,'Len12D_8',Len12D_8)
           endif
        else
           IntPrecision = 4
        endIf
        call Close_MatF(IU)

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

        if ( Call_Test .eq. 0 ) then
           IntPrecision = Save_IntPrecision
        else
           call MQC_Read_Mat_Precision_Test(MatFileName, IntPrecision)
           Save_IntPrecision = IntPrecision
        endif
        return
      end subroutine MQC_Mat_Precision

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

      Subroutine Wr_LRInd(IU,Label,NR,LenBuf,N1,N2,N3,N4,N5,ASym,X)
        Implicit None

        Character(Len=*),Intent(In)::Label
        Integer(kind=int64),Intent(In)::IU
        Logical,Intent(In)::ASym
        Integer(kind=int64),Intent(In)::NR,LenBuf,N1,N2,N3,N4,N5
        Real(kind=real64),Dimension(:,:),Intent(InOut)::X

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
        else
! File requires 64-bit integers, so no change
           call OpenRead8(Name,IU,LabFil,IVers,NLab,GVers,Title,NAtoms,NBasis,NBsUse,ICharg,Multip,NE,Len12L,Len4L,IOpCl,ICGU)
        endif
        if ( IU .eq. -1 ) then
           IOut = 6
           call MQC_Error_A('Unable to Open File', Iout, 'Name', Name)
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
        Call_Test = 1
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
        else
! File requires 64-bit integers, so no change
           call OpenWrite8(Name,IU,LabFil,GVers,Title,NAtoms,NBasis,NBsUse,ICharg,Multip,NE,IOpCl,ICGU)
        endif
      Return
      End Subroutine Open_Write

      Subroutine Close_MatF(IU)
        Implicit None
!
!     Close a matrix-element file.
!
        Integer(kind=int64),Intent(In)::IU
        Call CloseMatF8(IU)
      Return
      End Subroutine Close_MatF

      Function AOInts(CBuf) result(aoint)
        Implicit None
        Character(Len=*),Intent(InOut)::CBuf
        Logical aoint
        Logical, external::a4ints
        AOInt = A4Ints(CBuf)
      End Function AOInts

      Function LenArr(N1,N2,N3,N4,N5) result(LenAr)
        Implicit None
        Integer(kind=int64),Intent(In)::N1,N2,N3,N4,N5
        Integer(kind=int64)::LenAr

        Integer(kind=int64), external::LenAr8

        LenAr = LenAr8(N1,N2,N3,N4,N5)
      End Function LenArr

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

      Function NumNZA(NR,NTot,X) result(Numz)
        Implicit None

        Integer(kind=int64),Intent(InOUT)::NR,NTot
        Real(kind=real64),Dimension(NTOT,NR),Intent(InOut)::X

        Integer(kind=int64), external::N8mNZA
        Integer(kind=int64)::Numz

        Numz = N8mNZA(NR,NTot,X)
      End Function NumNZA

      Function NumNZR(NR,NTot,X) result(Numz)
        Implicit None

        Integer(kind=int64),Intent(InOUT)::NR,NTot
        Real(kind=real64),Dimension(NR,NTOT),Intent(InOut)::X

        Integer(kind=int64), external::N8mNZr
        Integer(kind=int64)::Numz

        Numz = N8mNZR(NR,NTot,X)
      End Function NumNZR

      Subroutine Rd_2E1(IU,LR,NTot,LenBuf,RArr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,LR,NTOT,LenBuf
        Real(kind=real64),Dimension(LR),Intent(InOut)::RArr

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

      Subroutine Rd_2EN(IU,NR,LR,LRNR,NTot,LenBuf,RArr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NR,LR,LRNR,NTot,LenBuf
        Real(kind=real64),Dimension(LRNR),Intent(InOut)::RArr

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
           call RdCBuf8(IU_4,NTot_4,LenBuf_4,Arr)
        else
! File requires 64-bit integers, so no change
           call RdCBuf8(IU,NTot,LenBuf,Arr)
        endif
      Return
      End Subroutine Rd_CBuf

      Subroutine Rd_Head(IU,NLab,NAtoms,NBasis,IAn,IAtTyp,AtmChg,C, &
       IBfAtm,IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB, &
       NPrmDB,NBTot)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NAtoms,NBasis,NLab
        Integer(kind=int64),Dimension(NAtoms),Intent(Out)::IAn
        Integer(kind=int64),Dimension(NAtoms),Intent(Out)::IAtTyp
        Real(kind=real64),Dimension(NAtoms),Intent(Out)::AtmChg
        Real(kind=real64),Dimension(3*NAtoms),Intent(Out)::C
        Integer(kind=int64),Dimension(NBasis),Intent(Out)::IBfAtm
        Integer(kind=int64),Dimension(NBasis),Intent(Out)::IBfTyp
        Real(kind=real64),Dimension(NAtoms),Intent(Out)::AtmWgt
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

      write( 6, 1234) iu
1234  format( " DJM iu in rd_labl is ", I8 )

           Call RdLabl8(IU,IVers,CBuf,NI,NR,NTot,LenBuf,N1,N2,N3,N4,N5,ASym,NRI,EOF)
        endif
      Return
      End Subroutine Rd_Labl

      Subroutine Rd_RBuf(IU,NTot,LenBuf,Arr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NTot,LenBuf
        Real(kind=real64),Dimension(NTot),Intent(Out)::Arr

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

      Subroutine Rd_RInd(IU,NR,LR,NTot,LenBuf,LNZ,RArr)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NR,LR,NTot,LenBuf
        Integer(kind=int64),Intent(Out)::LNZ
        Real(kind=real64),Dimension(NR,LR),Intent(Out)::RArr

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
           Call RdRInd8(IU,NR,LR,NTot,LenBuf,LNZ,RArr)
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
        Real(kind=real64),Dimension(LR,NR),Intent(Out)::RArr

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

      Subroutine Wr_Head(IU,NAtoms,NAt3,NBasis,IAn,IAtTyp,AtmChg,C,IBfAtm, &
           IBfTyp,AtmWgt,NFC,NFV,ITran,IDum9,NShlAO,NPrmAO,NShlDB,NPrmDB,NBTot)
        Implicit None

        Integer(kind=int64),Intent(In)::IU,NAtoms,NAt3,NBasis
        INTEGER(kind=int64),Dimension(NAtoms),Intent(In)::IAn
        INTEGER(kind=int64),Dimension(NAtoms),Intent(In)::IAtTyp
        REAL(kind=real64),Dimension(NAtoms),Intent(In)::AtmChg
        REAL(kind=real64),Dimension(NAt3),Intent(In)::C
        INTEGER(kind=int64),Dimension(NBasis),Intent(In)::IBfTyp
        INTEGER(kind=int64),Dimension(NBasis),Intent(In)::IBfAtm
        REAL(kind=real64),Dimension(NAtoms),Intent(In)::AtmWgt
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
        Real(kind=real64),Dimension(NTot),Intent(In)::Arr

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
        Real(kind=real64),Dimension(NR,LR),Intent(In)::RArr

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
        Real(kind=real64),Dimension(LR),Intent(Out)::RI
        Real(kind=real64),Dimension(N,N,N,N),Intent(Out)::RO
! This routine does no I/O
        Call ExpA18(N,LR,RI,RO)
      Return
      End Subroutine ExpAO1

      Subroutine ExpAON(NE,N,LR,RI,RO)
        Implicit None

        Integer(kind=int64),Intent(In)::NE,N, LR
        Real(kind=real64),Dimension(NE,LR),Intent(Out)::RI
        Real(kind=real64),Dimension(NE,N,N,N,N),Intent(Out)::RO
! This routine does no I/O
        Call ExpAN8(NE,N,LR,RI,RO)
      Return
      End Subroutine

      Subroutine AClear(N,A)
        Implicit None

        Integer(kind=int64),Intent(In)::N
        Real(kind=real64),Dimension(N),Intent(Out)::A

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

      End Module MQC_MatWrapper
