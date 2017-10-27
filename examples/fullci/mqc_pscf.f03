      Module MQC_PSCF 
!
!     Call other Modules
!
      Use MQC_Algebra
      Use MQC_Gaussian
      Use MQC_DataStructures
!
!     Define Types and Classes.
!

!
!     Define Procedure Interfaces.
!

!
!     Define Operators.
!

!
!     Subroutines/Functions...
!
      Contains

!======================================================================================================
!     
!     SUBROUTINE GEN_DET_STR 
! 
      Subroutine Gen_Det_Str(IOut,IPrint,NBasis,NAlpha,NBeta, &
        Alpha_Strings,Beta_Strings)
!
!     This subroutine generates a list of alpha and beta strings in binary  
!     notation following lexical order 
!     I'm limiting this to one bit as for out purposes it doesn't make sense 
!     to go beyond this. If you need more then allocate the alpha and beta
!     string in the routine before this and then the strings will be built in
!     the first element
!
!     Variable Declarations...
!
      Implicit None

      Type Node
        Integer::i
        Real::r
        Character::a
      End Type Node

      Integer::IOut,IPrint,NBasis,NAlpha,NBeta,NAlpha_Str,NBeta_Str,NBit_Ints,IOrb,&
        NElec,NEMax,NEMin,I
      Integer,Dimension(:,:),Allocatable::Alpha_Strings,Beta_Strings
      Type(MQC_LinkedList),Pointer::String_List_1=>Null(),String_List_2=>Null(),&
        String_Node_1=>Null(),String_Node_2=>Null()
      Type(Node)::New_Value
      Type(Node),Allocatable::Returned_Value
      Logical::Last
!
!     Initialize Arrays
!
      NAlpha_Str = Bin_Coeff(NBasis,NAlpha)
      NBeta_Str = Bin_Coeff(NBasis,NBeta)
      If(.not.Allocated(Alpha_Strings).and..not.Allocated(Beta_Strings)) then
        NBit_Ints = (NBasis/Bit_Size(0))+1
        If(NBit_Ints.gt.1) then
          Write(IOut,*) "ERROR: Determinant generator limited to ",Bit_Size(0)," orbitals"
          Call Exit()
        EndIf
        Allocate(Alpha_Strings(NAlpha_Str,NBit_Ints),Beta_Strings(NBeta_Str,NBit_Ints))
      EndIf
      Alpha_Strings = 0 
      Beta_Strings = 0
!
!    Forming alpha strings
      New_Value%i = 0
      Call LinkedList_Push(String_List_2,New_Value)
!    Still need to take care of cases where  number of orbitals is greater than one bit
      Do IOrb=1, NBasis
        Call LinkedList_Delete(String_List_1)
        String_Node_2 => String_List_2
        Allocate(String_List_1)
        If(Associated(String_Node_2)) then
          Last = .False.
          Do While(.not.Last)
            Call LinkedList_GetNext(String_Node_2,Last,.True.)
            Call LinkedList_Return_Value(String_Node_2,Returned_Value)
            Call LinkedList_Push(String_List_1,Returned_Value)
          EndDo
        EndIf
        Call LinkedList_Delete(String_List_2)
!      for each link in linked list 1
        String_Node_1 => String_List_1
        Allocate(String_List_2)
        If(Associated(String_Node_1)) then
          Last = .False.
          Do While(.not.Last)
            Call LinkedList_GetNext(String_Node_1,Last,.True.)
            Call LinkedList_Return_Value(String_Node_1,Returned_Value)
!         bit number zero
            New_Value%i = IBClr(Returned_Value%i,IOrb-1)
            NElec = PopCnt(New_Value%i)
            NEMax = NElec + (NBasis-IOrb)
            NEMin = NElec
            If(NEMax.ge.NAlpha.and.NEMin.le.NAlpha) then
              Call LinkedList_Push(String_List_2,New_Value)
            EndIf
!         bit number one
            New_Value%i = IBSet(Returned_Value%i,IOrb-1)
            NElec = PopCnt(New_Value%i)
            NEMax = NElec + (NBasis-IOrb)
            NEMin = NElec
            If(NEMax.ge.NAlpha.and.NEMin.le.NAlpha) then 
              Call LinkedList_Push(String_List_2,New_Value)
            EndIf
          EndDo
        EndIf
      EndDo 

      String_Node_2 => String_List_2
      If(Associated(String_Node_2)) then
        Last = .False.
        I = 1
        Do While(.not.Last)
          Call LinkedList_GetNext(String_Node_2,Last,.True.)
          Call LinkedList_Return_Value(String_Node_2,Returned_Value)
          Alpha_Strings(I,1) = Returned_Value%i
          I=I+1
        EndDo
      EndIf
!
!    Forming beta strings
      Call LinkedList_Delete(String_List_2)
      Allocate(String_List_2)
      New_Value%i = 0
      Call LinkedList_Push(String_List_2,New_Value)
!    Still need to take care of cases where  number of orbitals is greater than one bit
      Do IOrb=1, NBasis
        Call LinkedList_Delete(String_List_1)
        String_Node_2 => String_List_2
        Allocate(String_List_1)
        If(Associated(String_Node_2)) then
          Last = .False.
          Do While(.not.Last)
            Call LinkedList_GetNext(String_Node_2,Last,.True.)
            Call LinkedList_Return_Value(String_Node_2,Returned_Value)
            Call LinkedList_Push(String_List_1,Returned_Value)
          EndDo
        EndIf
        Call LinkedList_Delete(String_List_2)
!      for each link in linked list 1
        String_Node_1 => String_List_1
        Allocate(String_List_2)
        If(Associated(String_Node_1)) then
          Last = .False.
          Do While(.not.Last)
            Call LinkedList_GetNext(String_Node_1,Last,.True.)
            Call LinkedList_Return_Value(String_Node_1,Returned_Value)
!         bit number zero
            New_Value%i = IBClr(Returned_Value%i,IOrb-1)
            NElec = PopCnt(New_Value%i)
            NEMax = NElec + (NBasis-IOrb)
            NEMin = NElec
            If(NEMax.ge.NBeta.and.NEMin.le.NBeta) then
              Call LinkedList_Push(String_List_2,New_Value)
            EndIf
!         bit number one
            New_Value%i = IBSet(Returned_Value%i,IOrb-1)
            NElec = PopCnt(New_Value%i)
            NEMax = NElec + (NBasis-IOrb)
            NEMin = NElec
            If(NEMax.ge.NBeta.and.NEMin.le.NBeta) then 
              Call LinkedList_Push(String_List_2,New_Value)
            EndIf
          EndDo
        EndIf
      EndDo 

      String_Node_2 => String_List_2
      If(Associated(String_Node_2)) then
        Last = .False.
        I = 1
        Do While(.not.Last)
          Call LinkedList_GetNext(String_Node_2,Last,.True.)
          Call LinkedList_Return_Value(String_Node_2,Returned_Value)
          Beta_Strings(I,1) = Returned_Value%i
          I=I+1
        EndDo
      EndIf
!
      If(IPrint.ge.2) then
        Write(IOut,*) "Alpha Strings"
        Write(IOut,'(B64)') Alpha_Strings(:,1)
        Write(IOut,*) "Beta Strings"
        Write(IOut,'(B64)') Beta_Strings(:,1)
      EndIf

      End Subroutine Gen_Det_Str
!
!======================================================================================================
!     
!     FUNCTION SLATER_CONDON
! 
      Function Slater_Condon(IOut,IPrint,NBasis,Alpha_String_1,Beta_String_1, &
      Alpha_String_2,Beta_String_2,Core_Hamiltonian,ERIs,UHF)
!
!     This function returns the CI Hamiltonian matrix element value using 
!     Slater-Condon rules for a given alpha and beta string combination.
!     If UHF flag is set to true then Slater Condon rules are returned
!     accounting for different spacial occupation of alpha and beta orbitals but
!     the ERIs must be dimensioned to twice that of rhf case.
!
!     Variable Declarations...
!
      Implicit None
      Integer::IOut,IPrint,NBasis,IPos,JPos,IDiff,Det_Diff,ISgn,NBit_Ints,I,J,K, &
        II,JJ,Alpha_Diff_Cnt,Beta_Diff_Cnt
      Real::Slater_Condon,ERI1,ERI2,Zero=0.0d0
      Logical::UHF
      Type(MQC_Core_Hamiltonian)::Core_Hamiltonian
      Integer,Dimension(4)::Orbs,Spin,Det
      Integer,Dimension(:)::Alpha_String_1,Alpha_String_2,Beta_String_1, &
        Beta_String_2
      Integer,Dimension(:),Allocatable::Alpha_Diff,Beta_Diff
      Real,Dimension(:,:,:,:),Allocatable::ERIs
!
!     Initialize Arrays
!
!      Write(IOut,*) '---alpha 1:---'
!      Write(IOut,'(B64)') Alpha_String_1
!      Write(IOut,*) '---alpha 2:---'
!      Write(IOut,'(B64)') Alpha_String_2
!      Write(IOut,*) '-------------'
!      Write(IOut,*) 'alpha XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Alpha_String_1,Alpha_String_2)   
!      Write(IOut,*) 'NDifa=', PopCnt(IEOR(Alpha_String_1,Alpha_String_2))
!      Write(IOut,*)
!      Write(IOut,*) '--beta 1:--'
!      Write(IOut,'(B64)') Beta_String_1
!      Write(IOut,*) '---beta 2:---'
!      Write(IOut,'(B64)') Beta_String_2
!      Write(IOut,*) '-------------'
!      Write(IOut,*) 'beta XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Beta_String_1,Beta_String_2)   
!      Write(IOut,*) 'NDifb=', PopCnt(IEOR(Beta_String_1,Beta_String_2))
!      Write(IOut,*)

      NBit_Ints = (NBasis/Bit_Size(0))+1 
      Allocate(Alpha_Diff(NBit_Ints),Beta_Diff(NBit_Ints))
      Det_Diff = 0
      Alpha_Diff_Cnt = 0
      Beta_Diff_Cnt = 0
      Do I = 1,NBit_Ints 
        Alpha_Diff(I) = IEOR(Alpha_String_1(I),Alpha_String_2(I))
!        Write(IOut,*) 'Alpha Diff',I,':'
!        Write(IOut,'(B64)') Alpha_Diff(I)
!        Write(IOut,*) '-------------'
        Alpha_Diff_Cnt = Alpha_Diff_Cnt + PopCnt(Alpha_Diff(I)) 
        Beta_Diff(I) = IEOR(Beta_String_1(I),Beta_String_2(I))
!        Write(IOut,*) 'Beta Diff',I,':'
!        Write(IOut,'(B64)') Beta_Diff(I)
!        Write(IOut,*) '-------------'
        Beta_Diff_Cnt = Beta_Diff_Cnt + PopCnt(Beta_Diff(I))
      EndDo
!      Write(IOut,*)'Alpha_Diff_Cnt:',Alpha_Diff_Cnt,'Beta_Diff_Cnt:',Beta_Diff_Cnt
      Det_Diff = Alpha_Diff_Cnt/2 + Beta_Diff_Cnt/2

      If(Mod(Alpha_Diff_Cnt,2).ne.0.or.Mod(Beta_Diff_Cnt,2).ne.0) then
        Write(IOut,*) "ERROR: Slater_Condon has been handed spin non-conserving &
        determinants"
        Call Exit()
      EndIf

!      Write(IOut,*) "Det_Diff:",Det_Diff
      Select Case (Det_Diff)

        Case(3:)
          Slater_Condon = Zero 
          Return
!
        Case(2)
!       I am going to comment the first logical block - the rest are different
!       perutations
!       First we need to determine the relevant orbital, spin and determinant
          IDiff = 1
!          Allocate(Orbs(4),Spin(4),Det(4))
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
!            Write(IOut,*) 'Orbs:',IPos,' I:',I,' J:',J
            If(BTest(Alpha_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Call Print_Vector(IOut,Orbs,'Orbs')
!          Call Print_Vector(IOut,Spin,'Spin')
!          Call Print_Vector(IOut,Det,'Det')

          If(Det(1).eq.Det(2).and.Det(3).eq.Det(4)) then
!            Write(IOut,*) 'Entering case 1'
!            
!           This section computes the Coulomb contribution to the matrix
!           element
            If(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
!              
!             This section detemines number of permutations of creation
!             operators and sets sign accordingly
!             It got quite a bit longer after generalizing max number of
!             orbitals
!             Probably want to make some functions for these then later
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(2).eq.0) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
!              
!             This section computes the value of the matrix element
!
              If((UHF.eq..False.).or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
                ERI1 = ERIs(Orbs(1),Orbs(3),Orbs(2),Orbs(4))
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
                ERI1 = ERIs(Orbs(1),Orbs(3),Orbs(2)+NBasis,Orbs(4)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                ERI1 = ERIs(Orbs(1)+NBasis,Orbs(3)+NBasis,Orbs(2),Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
                ERI1 = ERIs(Orbs(1)+NBasis,Orbs(3)+NBasis,Orbs(2)+NBasis,Orbs(4)+NBasis)
              EndIf
!              
!             This section changes sign accounting for the ordering of alpha and
!             beta strings
              If(Spin(1).gt.Spin(2)) ERI1 = -ERI1
              If(Spin(3).gt.Spin(4)) ERI1 = -ERI1
              If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!              Write(IOut,*) 'Setting ERI1 to',ERI1
            Else
!              
!             Sets Coulomb contribution to zero if 2ERI is necessarily so by
!             spin order in the integral
!              Write(IOut,*) 'Setting ERI1 to zero'
              ERI1 = Zero
            EndIf
            If(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(2).eq.0) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
!
              If((UHF.eq..False.).or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
                ERI2 = ERIs(Orbs(1),Orbs(4),Orbs(2),Orbs(3))
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.1) then
                ERI2 = ERIs(Orbs(1),Orbs(4),Orbs(2)+NBasis,Orbs(3)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                ERI2 = ERIs(Orbs(1)+NBasis,Orbs(4)+NBasis,Orbs(2),Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
                ERI2 = ERIs(Orbs(1)+NBasis,Orbs(4)+NBasis,Orbs(2)+NBasis,Orbs(3)+NBasis)
              EndIf
!
              If(Spin(1).gt.Spin(2)) ERI2 = -ERI2
              If(Spin(3).gt.Spin(4)) ERI2 = -ERI2
              If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!              Write(IOut,*) 'Setting ERI2 to',ERI2
            Else 
!              Write(IOut,*) 'Setting ERI2 to zero'
              ERI2 = Zero
            EndIf
!
          ElseIf(Det(1).eq.Det(3).and.Det(2).eq.Det(4)) then
!            Write(IOut,*) 'Entering case 2'
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(3).eq.0) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
              If(UHF.eq..False..or.(Spin(1).eq.0.and.Spin(3).eq.0)) then
                ERI1 = ERIs(Orbs(1),Orbs(2),Orbs(3),Orbs(4))
              ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
                ERI1 = ERIs(Orbs(1),Orbs(2),Orbs(3)+NBasis,Orbs(4)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                ERI1 = ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,Orbs(3),Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
                ERI1 = ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,Orbs(3)+NBasis,Orbs(4)+NBasis)
              EndIf
              If(Spin(1).gt.Spin(3)) ERI1 = -ERI1
              If(Spin(2).gt.Spin(4)) ERI1 = -ERI1
              If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!              Write(IOut,*) 'Setting ERI1 to',ERI1
            Else
!              Write(IOut,*) 'Setting ERI1 to zero'
              ERI1 = Zero
            EndIf
            If(Spin(1).eq.Spin(4).and.Spin(3).eq.Spin(2)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(3).eq.0) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
                If(Orbs(3)-Orbs(2)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(2),Orbs(3),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(4)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
              If(UHF.eq..False..or.(Spin(1).eq.0.and.Spin(3).eq.0)) then
                ERI2 = ERIs(Orbs(1),Orbs(4),Orbs(3),Orbs(2))
              ElseIf(Spin(1).eq.0.and.Spin(3).eq.1) then
                ERI2 = ERIs(Orbs(1),Orbs(4),Orbs(3)+NBasis,Orbs(2)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                ERI2 = ERIs(Orbs(1)+NBasis,Orbs(4)+NBasis,Orbs(3),Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.1) then
                ERI2 = ERIs(Orbs(1)+NBasis,Orbs(4)+NBasis,Orbs(3)+NBasis,Orbs(2)+NBasis)
              EndIf
              If(Spin(1).gt.Spin(3)) ERI2 = -ERI2
              If(Spin(2).gt.Spin(4)) ERI2 = -ERI2
              If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!              Write(IOut,*) 'Setting ERI2 to',ERI2
            Else 
!              Write(IOut,*) 'Setting ERI2 to zero'
              ERI2 = Zero
            EndIf
!
          ElseIf(Det(1).eq.Det(4).and.Det(2).eq.Det(3)) then
!            Write(IOut,*) 'Entering case 3'
            If(Spin(1).eq.Spin(2).and.Spin(4).eq.Spin(3)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(4).eq.0) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
                If(Orbs(4)-Orbs(3)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(3),Orbs(4),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(2)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
!                    Write(IOut,*) 'Orbs:',Orbs(1),Orbs(2),' I:',I,' J:',J,' II:',II,' JJ:',JJ
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
              If(UHF.eq..False..or.(Spin(1).eq.0.and.Spin(4).eq.0)) then
                ERI1 = ERIs(Orbs(1),Orbs(2),Orbs(4),Orbs(3))
              ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
                ERI1 = ERIs(Orbs(1),Orbs(2),Orbs(4)+NBasis,Orbs(3)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                ERI1 = ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,Orbs(4),Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
                ERI1 = ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,Orbs(4)+NBasis,Orbs(3)+NBasis)
              EndIf
              If(Spin(1).gt.Spin(4)) ERI1 = -ERI1
              If(Spin(2).gt.Spin(3)) ERI1 = -ERI1
              If(Mod(ISgn,2).eq.1) ERI1 = -ERI1
!              Write(IOut,*) 'Setting ERI1 to',ERI1
            Else
!              Write(IOut,*) 'Setting ERI1 to zero'
              ERI1 = Zero
            EndIf
            If(Spin(1).eq.Spin(3).and.Spin(4).eq.Spin(2)) then
              ISgn = 0
              If(Spin(1).eq.0.and.Spin(4).eq.0) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
                    ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then 
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
                If(Orbs(4)-Orbs(2)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(4)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(4)-1),Bit_Size(0)) + 1
                    EndIf
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
                If(Orbs(3)-Orbs(1)-1.gt.0) then
                  I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
                  J  = NBit_Ints - (Orbs(3)-1)/Bit_Size(0)
                  Do K = I, J, -1
                    If(I.ne.K) then 
                      II = 1
                    Else 
                      II = Mod((Orbs(1)-1),Bit_Size(0)) + 1
                    EndIf
                    If(J.ne.K) then 
                      JJ = Bit_Size(0)
                    Else
                      JJ = Mod((Orbs(3)-1),Bit_Size(0)) + 1
                    EndIf
                    ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
                  EndDo
                EndIf
              EndIf
              If(UHF.eq..False..or.(Spin(1).eq.0.and.Spin(4).eq.0)) then
                ERI2 = ERIs(Orbs(1),Orbs(3),Orbs(4),Orbs(2))
              ElseIf(Spin(1).eq.0.and.Spin(4).eq.1) then
                ERI2 = ERIs(Orbs(1),Orbs(3),Orbs(4)+NBasis,Orbs(2)+NBasis)
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                ERI2 = ERIs(Orbs(1)+NBasis,Orbs(3)+NBasis,Orbs(4),Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.1) then
                ERI2 = ERIs(Orbs(1)+NBasis,Orbs(3)+NBasis,Orbs(4)+NBasis,Orbs(2)+NBasis)
              EndIf
              If(Spin(1).gt.Spin(4)) ERI2 = -ERI2
              If(Spin(2).gt.Spin(3)) ERI2 = -ERI2
              If(Mod(ISgn,2).eq.1) ERI2 = -ERI2
!              Write(IOut,*) 'Setting ERI2 to',ERI2
            Else 
!              Write(IOut,*) 'Setting ERI2 to zero'
              ERI2 = Zero
            EndIf
          Else
            Write(IOut,*) 'Slater Condon has been handed confusing determinant info'
          EndIf
!
          Slater_Condon = ERI1 - ERI2
!          Write(IOut,*) 'Slater Condon Final:',Slater_Condon
!          Deallocate(Orbs,Spin,Det)
          Return
!
        Case(1)
          IDiff = 1
!          Allocate(Orbs(2),Spin(2))
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Write(IOut,*)'Orb 1:',Orbs(1),' Orb 2:',Orbs(2)
!          Write(IOut,*)'Spin 1:',Spin(1),' Spin 2:',Spin(2)
!          Write(IOut,*)'Det 1:',Det(1),' Det 2:',Det(2)
!
          ISgn = 0
          If(Spin(1).eq.0.and.Spin(2).eq.0) then
            If(Orbs(2)-Orbs(1)-1.gt.0) then 
              I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
              J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
              Do K = I, J, -1
                If(I.ne.K) then 
                  II = 1
                Else 
                  II = Mod((Orbs(1)-1),Bit_Size(0)) + 1 
                EndIf
                If(J.ne.K) then 
                  JJ = Bit_Size(0) 
                Else
                  JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                EndIf
                ISgn = ISgn + PopCnt(IBits(IAnd(Alpha_String_1(K),Alpha_String_2(K)),II,JJ-II-1))
              EndDo
            EndIf
          ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
            If(Orbs(2)-Orbs(1)-1.gt.0) then
              I  = NBit_Ints - (Orbs(1)-1)/Bit_Size(0)
              J  = NBit_Ints - (Orbs(2)-1)/Bit_Size(0)
              Do K = I, J, -1
                If(I.ne.K) then 
                  II = 1
                Else 
                  II = Mod((Orbs(1)-1),Bit_Size(0)) + 1 
                EndIf
                If(J.ne.K) then 
                  JJ = Bit_Size(0) 
                Else
                  JJ = Mod((Orbs(2)-1),Bit_Size(0)) + 1
                EndIf
                ISgn = ISgn + PopCnt(IBits(IAnd(Beta_String_1(K),Beta_String_2(K)),II,JJ-II-1))
              EndDo
            EndIf
          EndIf
          ISgn = (-1)**ISgn
!          Write(IOut,*)'ISgn:',ISgn
!
          Slater_Condon = Zero
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
!            Write(IOut,*) 'Alpha IPos:', IPos
            If(BTest(Alpha_String_1(I),J).eq..True.) then
              If(UHF.eq..False..and.Spin(1).eq.Spin(2)) then
!                Write(IOut,*) 'ALPHA: UHF False and 1 and 2 same spin, adding:',ISgn*ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1)
                Slater_Condon = Slater_Condon + ISgn*ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1)
                If(Spin(1).eq.0.and.Spin(2).eq.0) then
!                  Write(IOut,*) 'ALPHA: UHF False and 1 and 2 both alpha, subtracting:',ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
                  Slater_Condon = Slater_Condon - ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.0) then
!                Write(IOut,*) 'ALPHA: UHF True and 1 and 2 both alpha, adding:', &
!                  ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1) - ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
                Slater_Condon = Slater_Condon + ISgn*ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1) - ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!                Write(IOut,*) 'ALPHA: UHF True and 1 and 2 both beta, adding:', ISgn*ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,IPos+1,IPos+1)
                Slater_Condon = Slater_Condon + ISgn*ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,IPos+1,IPos+1)
              EndIf
            EndIf
          EndDo
!
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0))
            If(BTest(Beta_String_1(I),J).eq..True.) then
              If(UHF.eq..False..and.Spin(1).eq.Spin(2)) then
!                Write(IOut,*) 'BETA: UHF False and 1 and 2 same spin, adding:',ISgn*ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1)
                Slater_Condon = Slater_Condon + ISgn*ERIs(Orbs(1),Orbs(2),IPos+1,IPos+1)
                If(Spin(1).eq.1.and.Spin(2).eq.1) then
!                  Write(IOut,*) 'BETA: UHF False and 1 and 2 both beta, subtracting:',ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
                  Slater_Condon = Slater_Condon - ISgn*ERIs(Orbs(1),IPos+1,IPos+1,Orbs(2))
                EndIf
              ElseIf(Spin(1).eq.0.and.Spin(2).eq.0) then
!                Write(IOut,*) 'BETA: UHF True and 1 and 2 both alpha, adding:', ISgn*ERIs(Orbs(1),Orbs(2),IPos+1+NBasis,IPos+1+NBasis)
                Slater_Condon = Slater_Condon + ISgn*ERIs(Orbs(1),Orbs(2),IPos+1+NBasis,IPos+1+NBasis) 
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!                Write(IOut,*) 'BETA: UHF True and 1 and 2 both beta, adding:', &
!                ISgn*ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,IPos+1+NBasis,IPos+1+NBasis) - &
!                  ISgn*ERIs(Orbs(1)+NBasis,IPos+1+NBasis,IPos+1+NBasis,Orbs(2)+NBasis)
                Slater_Condon = Slater_Condon + ISgn*ERIs(Orbs(1)+NBasis,Orbs(2)+NBasis,IPos+1+NBasis,IPos+1+NBasis) - &
                  ISgn*ERIs(Orbs(1)+NBasis,IPos+1+NBasis,IPos+1+NBasis,Orbs(2)+NBasis)
              EndIf
            EndIf
          EndDo
!
          If(UHF.eq..False..or.(Spin(1).eq.0.and.Spin(2).eq.0)) then
!            Write(IOut,*) '1 and 2 both alpha, adding core:', ISgn*Core_Hamiltonian%Alpha(Orbs(1),Orbs(2))
            Slater_Condon = Slater_Condon + ISgn*Core_Hamiltonian%Alpha(Orbs(1),Orbs(2)) 
          ElseIf(Spin(1).eq.1.and.Spin(2).eq.1) then
!            Write(IOut,*) '1 and 2 both beta, adding core:', ISgn*Core_Hamiltonian%Beta(Orbs(1),Orbs(2))
            Slater_Condon = Slater_Condon + ISgn*Core_Hamiltonian%Beta(Orbs(1),Orbs(2))
          EndIf
          
!          Write(IOut,*) 'Final integral evaluated as:', Slater_Condon

!          Deallocate(Orbs,Spin)
          Return
!
        Case(0)
          Slater_Condon = Zero
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0))
            If(BTest(Alpha_String_1(I),J).eq..True.) then
!              Write(IOut,*) 'Alpha:',IPos+1,' I:',I,' J:',J
!              Write(IOut,*) '(I|h|I) =', Core_Hamiltonian%Alpha(IPos+1,IPos+1)
              Slater_Condon = Slater_Condon + Core_Hamiltonian%Alpha(IPos+1,IPos+1)
            EndIf
            If(BTest(Beta_String_1(I),J).eq..True.) then
!              Write(IOut,*) 'Beta:',IPos+1,' I:',I,' J:',J
!              Write(IOut,*) '(I|h|I) =', Core_Hamiltonian%Alpha(IPos+1,IPos+1)
              If(UHF.eq..True.) then
                Slater_Condon = Slater_Condon + Core_Hamiltonian%Beta(IPos+1,IPos+1)
              Else
                Slater_Condon = Slater_Condon + Core_Hamiltonian%Alpha(IPos+1,IPos+1)
              EndIf
            EndIf
          EndDo
!
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = IPos+1, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0)) 
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                If(BTest(Alpha_String_1(II),JJ).eq..True.) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
!                  Write(IOut,*)'(II|JJ) =', ERIs(IPos+1,IPos+1,JPos+1,JPos+1),'(IJ|JI) =',ERIs(IPos+1,JPos+1,JPos+1,IPos+1) 
                  Slater_Condon = Slater_Condon + ERIs(IPos+1,IPos+1,JPos+1,JPos+1) - ERIs(IPos+1,JPos+1,JPos+1,IPos+1)
                EndIf
              EndIf
            EndDo
          EndDo
!
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = IPos+1, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0))
              If(BTest(Beta_String_1(I),J).eq..True.) then
                If(BTest(Beta_String_1(II),JJ).eq..True.) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
!                  Write(IOut,*)'(II|JJ) =', ERIs(IPos+1,IPos+1,JPos+1,JPos+1),'(IJ|JI) =',ERIs(IPos+1,JPos+1,JPos+1,IPos+1) 
                  If(UHF.eq..True.) then
                    Slater_Condon = Slater_Condon + ERIs(IPos+1+NBasis,IPos+1+NBasis,JPos+1+NBasis,JPos+1+NBasis) - & 
                      ERIs(IPos+1+NBasis,JPos+1+NBasis,JPos+1+NBasis,IPos+1+NBasis)
                  Else
                    Slater_Condon = Slater_Condon + ERIs(IPos+1,IPos+1,JPos+1,JPos+1) - ERIs(IPos+1,JPos+1,JPos+1,IPos+1)
                  EndIf
                EndIf
              EndIf
            EndDo
          EndDo
!
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = 0, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0))
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                If(BTest(Beta_String_1(II),JJ).eq..True.) then 
!                  Write(IOut,*) 'IPos:',IPos+1,'JPos:',JPos+1,' I:',I,' J:',J,' II:',II,' JJ:',JJ
!                  Write(IOut,*)'(II|JJ) =', ERIs(IPos+1,IPos+1,JPos+1,JPos+1),'(IJ|JI) =',Zero
                  If(UHF.eq..True.) then
                    Slater_Condon = Slater_Condon + ERIs(IPos+1,IPos+1,JPos+1+NBasis,JPos+1+NBasis)
                  Else
                    Slater_Condon = Slater_Condon + ERIs(IPos+1,IPos+1,JPos+1,JPos+1)
                  EndIf
                EndIf
              EndIf
            EndDo
          EndDo

!          Write(IOut,*) 'Slater Condon:',Slater_Condon

          Return
!        
        Case Default
          Write(IOut,*) "ERROR: Slater_Condon is confused about number of &
            different orbitals"
          Call Exit()
!
      End Select
!
      End Function Slater_Condon
!
!======================================================================================================
!     
!     SUBROUTINE TWOERI_TRANS
! 
      Subroutine TwoERI_Trans(IOut,IPrint,MO_Coeff,ERIs,MO_ERIs,UHF)
!
!     This subroutine transforms two-electron integrals to the MO basis 
!     If the UHF flag is true it computes the set of MO integrals accounting
!     for different spacial occupation of alpha and beta electrons. The order
!     in the final MO ERI array is (aa|aa), (aa|bb), (bb|aa), (bb|bb).
!
!             |' (aa|aa) | (bb|aa) '| 
!             |----------|----------|
!             |, (aa|bb) | (bb|bb) ,|
!
!     Variable Declarations...
!
      Implicit None
      Integer::IOut,IPrint,NBasis,P,Q,R,S,Mu,Nu,Lambda,Sigma,IErr
      Real::Zero=0.0d0
      Logical::Noddy,Smart,UHF
      Real,Dimension(:,:),Allocatable::X,Y
      Real,Intent(In),Dimension(:,:,:,:),Allocatable::ERIs
      Real,Intent(Out),Dimension(:,:,:,:),Allocatable::MO_ERIs
      Type(MQC_MO_Coefficients)::MO_Coeff
!
!     Initialize Arrays
!
!      Noddy = .True.
!      Smart = .False.
      Noddy = .False.
      Smart = .True.

      NBasis = Size(ERIs,1)
      If(UHF.eq..False.) then
        Allocate(MO_ERIs(NBasis,NBasis,NBasis,NBasis))
      ElseIf(UHF.eq..True.) then
        IErr = 1
        Write(IOut,*) "IErr Before:",IErr
        Allocate(MO_ERIs(NBasis*2,NBasis*2,NBasis*2,NBasis*2),Stat=IErr)
        Write(IOut,*) "IErr After:",IErr
        If(IErr.ne.0) then
          Write(IOut,*) "Transformed arrays not allocated correctly"
          Call Exit()
        Else
          Write(IOut,*) "Transformed arrays allocated correctly"
        EndIf
      EndIf
      MO_ERIs(:,:,:,:) = Zero

      If(Noddy.eq..True.) then

        Write(IOut,*) 'Doing noddy integral transformation algorithm'
        Do P = 1, NBasis
          Do Q = 1, NBasis
            Do R = 1, NBasis
              Do S = 1, NBasis
                Do Mu = 1, NBasis
                  Do Nu = 1, NBasis
                    Do Lambda = 1, NBasis
                      Do Sigma = 1, NBasis
                        MO_ERIs(P,Q,R,S) =  MO_ERIs(P,Q,R,S) + MO_Coeff%Alpha(Mu,P)*MO_Coeff%Alpha(Nu,Q)* &
                          ERIs(Mu,Nu,Lambda,Sigma)*MO_Coeff%Alpha(Lambda,R)*MO_Coeff%Alpha(Sigma,S)
                      EndDo
                    EndDo
                  EndDo
                EndDo
                If(IPrint.ge.2) Write(IOut,*) '(',P,Q,'|',R,S,') =',MO_ERIs(P,Q,R,S)
              EndDo
            EndDo
          EndDo
        EndDo 

        If(UHF.eq..True.) then
          
          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Do Mu = 1, NBasis
                    Do Nu = 1, NBasis
                      Do Lambda = 1, NBasis
                        Do Sigma = 1, NBasis
                          MO_ERIs(P,Q,R+NBasis,S+NBasis) =  MO_ERIs(P,Q,R+NBasis,S+NBasis) + MO_Coeff%Alpha(Mu,P)*MO_Coeff%Alpha(Nu,Q)* &
                            ERIs(Mu,Nu,Lambda,Sigma)*MO_Coeff%Beta(Lambda,R)*MO_Coeff%Beta(Sigma,S)
                        EndDo
                      EndDo
                    EndDo
                  EndDo
                  If(IPrint.ge.2) Write(IOut,*) '(',P,Q,'|',R+NBasis,S+NBasis,') =',MO_ERIs(P,Q,R+NBasis,S+NBasis)
                EndDo
              EndDo
            EndDo
          EndDo 

          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Do Mu = 1, NBasis
                    Do Nu = 1, NBasis
                      Do Lambda = 1, NBasis
                        Do Sigma = 1, NBasis
                          MO_ERIs(P+NBasis,Q+NBasis,R,S) =  MO_ERIs(P+NBasis,Q+NBasis,R,S) + MO_Coeff%Beta(Mu,P)*MO_Coeff%Beta(Nu,Q)* &
                            ERIs(Mu,Nu,Lambda,Sigma)*MO_Coeff%Alpha(Lambda,R)*MO_Coeff%Alpha(Sigma,S)
                        EndDo
                      EndDo
                    EndDo
                  EndDo
                  If(IPrint.ge.2) Write(IOut,*) '(',P+NBasis,Q+NBasis,'|',R,S,') =',MO_ERIs(P+NBasis,Q+NBasis,R,S)
                EndDo
              EndDo
            EndDo
          EndDo 

          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  Do Mu = 1, NBasis
                    Do Nu = 1, NBasis
                      Do Lambda = 1, NBasis
                        Do Sigma = 1, NBasis
                          MO_ERIs(P+NBasis,Q+NBasis,R+NBasis,S+NBasis) =  MO_ERIs(P+NBasis,Q+NBasis,R+NBasis,S+NBasis) + MO_Coeff%Beta(Mu,P)*MO_Coeff%Beta(Nu,Q)* &
                            ERIs(Mu,Nu,Lambda,Sigma)*MO_Coeff%Beta(Lambda,R)*MO_Coeff%Beta(Sigma,S)
                        EndDo
                      EndDo
                    EndDo
                  EndDo
                  If(IPrint.ge.2) Write(IOut,*) '(',P+NBasis,Q+NBasis,'|',R+NBasis,S+NBasis,') =',MO_ERIs(P+NBasis,Q+NBasis,R+NBasis,S+NBasis)
                EndDo
              EndDo
            EndDo
          EndDo 

        EndIf

      ElseIf(Smart.eq..True.) then

        Write(IOut,*) 'Doing smart integral transformation algorithm'
        Allocate(X(NBasis,NBasis),Y(NBasis,NBasis))

        X = Zero 
        Y = Zero
        Do P = 1, NBasis
          Do Q = 1, NBasis
            Do R = 1, NBasis
              Do S = 1, NBasis
                X(R,S) = ERIs(P,Q,R,S) 
              EndDo
            EndDo
            Y = Zero
            Y = MatMul(Transpose(MO_Coeff%Alpha),X)
            X = Zero
            X = MatMul(Y,MO_Coeff%Alpha)
            Do R = 1, NBasis
              Do S = 1, NBasis
                MO_ERIs(P,Q,R,S) = X(R,S)
              EndDo
            EndDo
          EndDo
        EndDo

        Do R = 1, NBasis
          Do S = 1, NBasis
            X = Zero
            Y = Zero
            Do P = 1, NBasis
              Do Q = 1, NBasis 
                X(P,Q) = MO_ERIs(P,Q,R,S)
              EndDo
            EndDo
            Y = Zero
            Y = MatMul(Transpose(MO_Coeff%Alpha),X)
            X = Zero
            X = MatMul(Y,MO_Coeff%Alpha)
            Do P = 1, NBasis
              Do Q = 1, NBasis
                MO_ERIs(P,Q,R,S) = X(P,Q)
                If(IPrint.ge.2) Write(IOut,*) '(',P,Q,'|',R,S,') =',MO_ERIs(P,Q,R,S)
              EndDo
            EndDo
          EndDo
        EndDo

        If(UHF.eq..True.) then
         
          X = Zero 
          Y = Zero
          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  X(R,S) = ERIs(P,Q,R,S) 
                EndDo
              EndDo
              Y = Zero
              Y = MatMul(Transpose(MO_Coeff%Beta),X)
              X = Zero
              X = MatMul(Y,MO_Coeff%Beta)
              Do R = 1, NBasis
                Do S = 1, NBasis
                  MO_ERIs(P,Q,R+NBasis,S+NBasis) = X(R,S)
                EndDo
              EndDo
            EndDo
          EndDo

          Do R = 1, NBasis
            Do S = 1, NBasis
              X = Zero
              Y = Zero
              Do P = 1, NBasis
                Do Q = 1, NBasis 
                  X(P,Q) = MO_ERIs(P,Q,R+NBasis,S+NBasis)
                EndDo
              EndDo
              Y = Zero
              Y = MatMul(Transpose(MO_Coeff%Alpha),X)
              X = Zero
              X = MatMul(Y,MO_Coeff%Alpha)
              Do P = 1, NBasis
                Do Q = 1, NBasis
                  MO_ERIs(P,Q,R+NBasis,S+NBasis) = X(P,Q)
                  If(IPrint.ge.2) Write(IOut,*) '(',P,Q,'|',R+NBasis,S+NBasis,') =',MO_ERIs(P,Q,R+NBasis,S+NBasis)
                EndDo
              EndDo
            EndDo
          EndDo


          X = Zero 
          Y = Zero
          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  X(R,S) = ERIs(P,Q,R,S) 
                EndDo
              EndDo
              Y = Zero
              Y = MatMul(Transpose(MO_Coeff%Alpha),X)
              X = Zero
              X = MatMul(Y,MO_Coeff%Alpha)
              Do R = 1, NBasis
                Do S = 1, NBasis
                  MO_ERIs(P+NBasis,Q+NBasis,R,S) = X(R,S)
                EndDo
              EndDo
            EndDo
          EndDo

          Do R = 1, NBasis
            Do S = 1, NBasis
              X = Zero
              Y = Zero
              Do P = 1, NBasis
                Do Q = 1, NBasis 
                  X(P,Q) = MO_ERIs(P+NBasis,Q+NBasis,R,S)
                EndDo
              EndDo
              Y = Zero
              Y = MatMul(Transpose(MO_Coeff%Beta),X)
              X = Zero
              X = MatMul(Y,MO_Coeff%Beta)
              Do P = 1, NBasis
                Do Q = 1, NBasis
                  MO_ERIs(P+NBasis,Q+NBasis,R,S) = X(P,Q)
                  If(IPrint.ge.2) Write(IOut,*) '(',P+NBasis,Q+NBasis,'|',R,S,') =',MO_ERIs(P+NBasis,Q+NBasis,R,S)
                EndDo
              EndDo
            EndDo
          EndDo


          X = Zero 
          Y = Zero
          Do P = 1, NBasis
            Do Q = 1, NBasis
              Do R = 1, NBasis
                Do S = 1, NBasis
                  X(R,S) = ERIs(P,Q,R,S) 
                EndDo
              EndDo
              Y = Zero
              Y = MatMul(Transpose(MO_Coeff%Beta),X)
              X = Zero
              X = MatMul(Y,MO_Coeff%Beta)
              Do R = 1, NBasis
                Do S = 1, NBasis
                  MO_ERIs(P+NBasis,Q+NBasis,R+NBasis,S+NBasis) = X(R,S)
                EndDo
              EndDo
            EndDo
          EndDo

          Do R = 1, NBasis
            Do S = 1, NBasis
              X = Zero
              Y = Zero
              Do P = 1, NBasis
                Do Q = 1, NBasis 
                  X(P,Q) = MO_ERIs(P+NBasis,Q+NBasis,R+NBasis,S+NBasis)
                EndDo
              EndDo
              Y = Zero
              Y = MatMul(Transpose(MO_Coeff%Beta),X)
              X = Zero
              X = MatMul(Y,MO_Coeff%Beta)
              Do P = 1, NBasis
                Do Q = 1, NBasis
                  MO_ERIs(P+NBasis,Q+NBasis,R+NBasis,S+NBasis) = X(P,Q)
                  If(IPrint.ge.2) Write(IOut,*) '(',P+NBasis,Q+NBasis,'|',R+NBasis,S+NBasis,') =',MO_ERIs(P+NBasis,Q+NBasis,R+NBasis,S+NBasis)
                EndDo
              EndDo
            EndDo
          EndDo

        EndIf
        Deallocate(X,Y)
        
      EndIf
!
      End Subroutine TwoERI_Trans
!
!======================================================================================================
!
!     FUNCTION S2_MAT_ELEM    
! 
      Function S2_Mat_Elem(IOut,IPrint,NBasis,Alpha_String_1,Beta_String_1, &
      Alpha_String_2,Beta_String_2,MO_Overlap)
!
!     This function returns the CI S**2 matrix elements used for computing S**2
!     values of CI vectors for a given alpha and beta string combination.
!     The MO overlap matrix is required
!
!     Variable Declarations...
!
      Implicit None
      Integer::IOut,IPrint,NBasis,IPos,JPos,IDiff,Det_Diff,NAlpha,NBeta, &
        IOcc,JOcc,KOcc,LOcc,Mat_Sign,Alpha_Diff_Cnt,Beta_Diff_Cnt,NBit_Ints, &
        I,J,II,JJ
      Real::S2_Mat_Elem,Zero=0.0d0,Quarter=0.25d0,ABTerm,One=1.0d0
      Integer,Dimension(4)::Orbs,Spin,Det
      Integer,Dimension(:)::Alpha_String_1,Alpha_String_2,Beta_String_1, &
        Beta_String_2
      Integer,Dimension(:),Allocatable::Alpha_Diff,Beta_Diff
      Real,Dimension(:,:),Allocatable::MO_Overlap
!
!      Write(IOut,*) 'alpha 1:'
!      Write(IOut,'(B64)') Alpha_String_1
!      Write(IOut,*) 'alpha 2:'
!      Write(IOut,'(B64)') Alpha_String_2
!      Write(IOut,*) 'alpha XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Alpha_String_1,Alpha_String_2)   
!      Write(IOut,*) 'NDifa=', PopCnt(IEOR(Alpha_String_1,Alpha_String_2))
!      Write(IOut,*)
!      Write(IOut,*) 'beta 1:'
!      Write(IOut,'(B64)') Beta_String_1
!      Write(IOut,*) 'beta 2:'
!      Write(IOut,'(B64)') Beta_String_2
!      Write(IOut,*) 'beta XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Beta_String_1,Beta_String_2)   
!      Write(IOut,*) 'NDifb=', PopCnt(IEOR(Beta_String_1,Beta_String_2))
!      Write(IOut,*)
!
      NBit_Ints = (NBasis/Bit_Size(0))+1 
      Allocate(Alpha_Diff(NBit_Ints),Beta_Diff(NBit_Ints))
      Det_Diff = 0
      Alpha_Diff_Cnt = 0
      Beta_Diff_Cnt = 0
      NAlpha = 0
      NBeta = 0
      Do I = 1,NBit_Ints
        Alpha_Diff(I) = IEOR(Alpha_String_1(I),Alpha_String_2(I))
!        Write(IOut,*) 'Alpha Diff',I,':'
!        Write(IOut,'(B64)') Alpha_Diff(I)
!        Write(IOut,*) '-------------'
        Alpha_Diff_Cnt = Alpha_Diff_Cnt + PopCnt(Alpha_Diff(I)) 
        NAlpha = NAlpha + PopCnt(Alpha_String_1(I))
        Beta_Diff(I) = IEOR(Beta_String_1(I),Beta_String_2(I))
!        Write(IOut,*) 'Beta Diff',I,':'
!        Write(IOut,'(B64)') Beta_Diff(I)
!        Write(IOut,*) '-------------'
        Beta_Diff_Cnt = Beta_Diff_Cnt + PopCnt(Beta_Diff(I))
        NBeta = NBeta + PopCnt(Beta_String_1(I))
      EndDo
!      Write(IOut,*)'Alpha_Diff_Cnt:',Alpha_Diff_Cnt,'Beta_Diff_Cnt:',Beta_Diff_Cnt
      Det_Diff = Alpha_Diff_Cnt/2 + Beta_Diff_Cnt/2

      If(Mod(Alpha_Diff_Cnt,2).ne.0.or.Mod(Beta_Diff_Cnt,2).ne.0) then
        Write(IOut,*) "ERROR: S2_Mat_Elem has been handed spin non-conserving &
        determinants"
        Call Exit()
      EndIf

!      Write(IOut,*) 'Det_Diff:',Det_Diff
      Select Case (Det_Diff)
!
        Case(3:)
          S2_Mat_Elem = Zero 
          Return
!
        Case(2)
          IDiff = 1
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Call Print_Vector(IOut,Orbs,'Orbs')
!          Call Print_Vector(IOut,Spin,'Spin')
!          Call Print_Vector(IOut,Det,'Det')
!
          IOcc = 0
          Do IPos = Orbs(1)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(1).eq.0) then
              If(Det(1).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            ElseIf(Spin(1).eq.1) then
              If(Det(1).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'IOcc:',IOcc
          JOcc = 0
          Do IPos = Orbs(2)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(2).eq.0) then
              If(Det(2).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            ElseIf(Spin(2).eq.1) then
              If(Det(2).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'JOcc:',JOcc
          KOcc = 0
          Do IPos = Orbs(3)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(3).eq.0) then
              If(Det(3).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) KOcc = KOcc + 1
              ElseIf(Det(3).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) KOcc = KOcc + 1
              EndIf
            ElseIf(Spin(3).eq.1) then
              If(Det(3).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) KOcc = KOcc + 1
              ElseIf(Det(3).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) KOcc = KOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'KOcc:',KOcc
          LOcc = 0
          Do IPos = Orbs(4)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(4).eq.0) then
              If(Det(4).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) LOcc = LOcc + 1
              ElseIf(Det(4).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) LOcc = LOcc + 1
              EndIf
            ElseIf(Spin(4).eq.1) then
              If(Det(4).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) LOcc = LOcc + 1
              ElseIf(Det(4).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) LOcc = LOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'LOcc:',LOcc
!          Mat_Sign = -1
!          Mat_Sign = (-1)**(IOcc+JOcc+KOcc+LOcc-3)
!          Write(IOut,*) 'Permutations:',(IOcc+JOcc+KOcc+LOcc-3)
          Mat_Sign = (-1)**(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!          Write(IOut,*) 'Permutations:',(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!          Write(IOut,*) 'Mat_Sign:',Mat_Sign
!
          If(Det(1).eq.Det(2).and.Det(3).eq.Det(4)) then
            If(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(2).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(4)+NBasis)*MO_Overlap(Orbs(2)+NBasis,Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(2),Orbs(3)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(4))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              If(Spin(1).eq.0.and.Spin(2).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(3)+NBasis)*MO_Overlap(Orbs(2)+NBasis,Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(2),Orbs(4)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(3))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          ElseIf(Det(1).eq.Det(3).and.Det(2).eq.Det(4)) then
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(3).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(4)+NBasis)*MO_Overlap(Orbs(3)+NBasis,Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(3),Orbs(2)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(4))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              If(Spin(1).eq.0.and.Spin(3).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(2)+NBasis)*MO_Overlap(Orbs(3)+NBasis,Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(3),Orbs(4)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(2))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          ElseIf(Det(1).eq.Det(4).and.Det(2).eq.Det(3)) then
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(4).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(3)+NBasis)*MO_Overlap(Orbs(4)+NBasis,Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(4),Orbs(2)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(3))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(4).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(2)+NBasis)*MO_Overlap(Orbs(4)+NBasis,Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(4),Orbs(3)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(2))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          EndIf
!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem

          Return

        Case(1)
          IDiff = 1
!          Allocate(Orbs(2),Spin(2))
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Write(IOut,*)'Orb 1:',Orbs(1),' Orb 2:',Orbs(2)
!          Write(IOut,*)'Spin 1:',Spin(1),' Spin 2:',Spin(2)
!          Write(IOut,*)'Det 1:',Det(1),' Det 2:',Det(2)
!
          S2_Mat_Elem = Zero 
          If(Spin(1).ne.Spin(2)) then
            S2_Mat_Elem = Zero
!
          ElseIf(Spin(1).eq.0) then
!
            IOcc = 0
            Do IPos = Orbs(1)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(1).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'IOcc:',IOcc
            JOcc = 0
            Do IPos = Orbs(2)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(2).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'JOcc:',JOcc
!
            Do IPos = 0, NBasis-1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(BTest(Beta_String_1(I),J).eq..True.) then
                S2_Mat_Elem = S2_Mat_Elem + MO_Overlap(IPos+1+NBasis,Orbs(1)) * MO_Overlap(Orbs(2),IPos+1+NBasis)
              EndIf
            EndDo
!            S2_Mat_Elem = - S2_Mat_Elem
!            Write(IOut,*) 'Permutations:',(2*NAlpha+1-IOcc-JOcc)
!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NAlpha+1-IOcc-JOcc)
            S2_Mat_Elem = (-1)**(2*NAlpha+1-IOcc-JOcc) * S2_Mat_Elem
!            S2_Mat_Elem = (-1)**(IOcc+JOcc-1) * S2_Mat_Elem
!
          ElseIf(Spin(1).eq.1) then

            IOcc = 0
            Do IPos = Orbs(1)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(1).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'IOcc:',IOcc
            JOcc = 0
            Do IPos = Orbs(2)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(2).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'JOcc:',JOcc
!
            Do IPos = 0, NBasis-1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                S2_Mat_Elem = S2_Mat_Elem + MO_Overlap(IPos+1,Orbs(1)+NBasis) * MO_Overlap(Orbs(2)+NBasis,IPos+1)
              EndIf
            EndDo
!            S2_Mat_Elem = - S2_Mat_Elem
!            Write(IOut,*) 'Permutations:',(2*NBeta+1-IOcc-JOcc)
!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NBeta+1-IOcc-JOcc)
            S2_Mat_Elem = (-1)**(2*NBeta+1-IOcc-JOcc) * S2_Mat_Elem
!            S2_Mat_Elem = (-1)**(IOcc+JOcc-1) * S2_Mat_Elem

          EndIf

!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem

          Return
!
        Case(0)
          ABTerm = Zero
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = 0, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0)) 
              If((BTest(Alpha_String_2(I),J).eq..True.).and.(BTest(Beta_String_2(II),JJ).eq..True.)) then
                ABTerm = ABTerm + MO_Overlap(IPos+1,JPos+1+NBasis)*MO_Overlap(JPos+1+NBasis,IPos+1) 
              EndIf
            EndDo
          EndDo
          S2_Mat_Elem = Quarter*((NAlpha-NBeta)**2+2*(NAlpha+NBeta)) - ABTerm
!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem
!
          Return
!
      End Select
!
      End Function S2_Mat_Elem  
!
!======================================================================================================
!
      End Module MQC_PSCF 
