      Module MQC_DataStructures
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
      Use MQC_General
      Use MQC_Binary
      use iso_fortran_env, only: int64, real64
!
!
!     MQC_LinkedList Type
!
      Type MQC_LinkedList
        Private
          Class(*),Allocatable,Public::Val
          Type(MQC_LinkedList),Pointer::Previous=>null()
          Type(MQC_LinkedList),Pointer::Next=>null()
      End Type MQC_LinkedList

      Interface LinkedList_Return_Value
        Module Procedure LinkedList_Return_Integer_Value
        Module Procedure LinkedList_Return_Real_Value
        Module Procedure LinkedList_Return_Character_Value
        Module Procedure LinkedList_Return_Binary_Value
      End Interface


      CONTAINS
!
!     PROCEDURE LinkedList_Return_Integer_Value
      Subroutine LinkedList_Return_Integer_Value(LL_Cur,Return_Value)
!
!     This routine is used to return the value of a member of a linked
!     list. LL_Cur is an INPUT linked list member and Return_Value is an
!     OUTPUT argument that is allocatable of the same type as LL_Cur%Val.
!
!     -H. P. Hratchian, 2016
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Integer(kind=int64),Allocatable,Intent(InOut)::Return_Value
!
!     Allocate Return_Value, set it to LL_Cur%Val, and return to the
!     calling procedure.
!
      If(Associated(LL_Cur)) then
        If(Allocated(LL_Cur%Val)) then
          If(.not.Allocated(Return_Value)) Allocate(Return_Value)
          Select Type (Val => LL_Cur%Val)
          Type is (Integer)
            Return_Value = Val
          Class Default
            Call MQC_Error_I('Data types not compatible in LinkedList_Return_Integer_Value', 6)
          End Select
        endIf
      endIf
!
      Return
      End Subroutine LinkedList_Return_Integer_Value
!
!     PROCEDURE LinkedList_Return_Real_Value
      Subroutine LinkedList_Return_Real_Value(LL_Cur,Return_Value)
!
!     This routine is used to return the value of a member of a linked
!     list. LL_Cur is an INPUT linked list member and Return_Value is an
!     OUTPUT argument that is allocatable of the same type as LL_Cur%Val.
!
!     -H. P. Hratchian, 2016
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Real(kind=real64),Allocatable,Intent(InOut)::Return_Value
!
!     Allocate Return_Value, set it to LL_Cur%Val, and return to the
!     calling procedure.
!
      If(Associated(LL_Cur)) then
        If(Allocated(LL_Cur%Val)) then
          If(.not.Allocated(Return_Value)) Allocate(Return_Value)
          Select Type (Val =>LL_Cur%Val)
          Type is (Real)
            Return_Value = Val
          Class Default
            Call MQC_Error_I('Data types not compatible in LinkedList_Return_Integer_Value',6)
          End Select
        endIf
      endIf
!
      Return
      End Subroutine LinkedList_Return_Real_Value
!
!     PROCEDURE LinkedList_Return_Character_Value
      Subroutine LinkedList_Return_Character_Value(LL_Cur,Return_Value)
!
!     This routine is used to return the value of a member of a linked
!     list. LL_Cur is an INPUT linked list member and Return_Value is an
!     OUTPUT argument that is allocatable of the same type as LL_Cur%Val.
!
!     -H. P. Hratchian, 2016
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Character(Len=*),Allocatable,Intent(InOut)::Return_Value
!
!     Allocate Return_Value, set it to LL_Cur%Val, and return to the
!     calling procedure.
!
      If(Associated(LL_Cur)) then
        If(Allocated(LL_Cur%Val)) then
          If(.not.Allocated(Return_Value)) Allocate(Return_Value)
          Select Type (Val =>LL_Cur%Val)
          Type is (Character(Len=*))
            Return_Value = Val
          Class Default
            Call MQC_Error_I('Data types not compatible in LinkedList_Return_Integer_Value', 6)
          End Select
        endIf
      endIf
!
      Return
      End Subroutine LinkedList_Return_Character_Value
!
!
!     PROCEDURE LinkedList_Return_Binary_Value
      Subroutine LinkedList_Return_Binary_Value(LL_Cur,Return_Value)
!
!     This routine is used to return the value of a member of a linked
!     list. LL_Cur is an INPUT linked list member and Return_Value is an
!     OUTPUT argument that is allocatable of the same type as LL_Cur%Val.
!
!     L. M. Thompson, 2024
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Type(MQC_Bits),Intent(InOut)::Return_Value
!
!     Allocate Return_Value, set it to LL_Cur%Val, and return to the
!     calling procedure.
!
      If(Associated(LL_Cur)) then
        If(Allocated(LL_Cur%Val)) then
          Select Type (Val => LL_Cur%Val)
          Type is (MQC_Bits)
            Return_Value = Val
          Class Default
            Call MQC_Error('Data types not compatible in LinkedList_Return_Binary_Value')
          End Select
        endIf
      endIf
!
      Return
      End Subroutine LinkedList_Return_Binary_Value
!
!
!     PROCEDURE LinkedList_Print_Node
      Subroutine LinkedList_Print_Node(LL_Cur,IOut)
!
!     This routine prints the current node.
!
!     L. M. Thompson, 2018
!
      Implicit None
      Integer(kind=int64)::IOut
      Class(*),Allocatable::Output
      Type(MQC_LinkedList),Pointer::LL_Cur
      If(Associated(LL_Cur)) then
        If(Allocated(LL_Cur%Val)) Allocate(Output,Source=LL_Cur%Val)
        Select Type(Output)
        Type is (Integer)
          Write(IOut,'(I3)',advance="no") Output
        Type is (Character(*))
          Write(IOut,'(A10)',advance="no") Output
        Type is (Real)
          Write(IOut,'(D15.5)',advance="no") Output
        Type is (MQC_Bits)
          Call Output%Print(IOut)
        Class Default
          Stop 'Unexpected data type in link'
        End Select
        Deallocate(Output)
      EndIf
      write(IOut,*)
!
      Return
      End Subroutine LinkedList_Print_Node
!
!     PROCEDURE LinkedList_Print
      Subroutine LinkedList_Print(LL_Head,IOut)
!
!     This routine prints the whole list. 
!
!     L. M. Thompson, 2018
!
      Implicit None
      Integer(kind=int64)::IOut
      Class(*),Allocatable::Output
      Type(MQC_LinkedList),Pointer::LL_Head, LL_Cur
      If(Associated(LL_Head)) then
        LL_Cur => LL_Head
        If(Associated(LL_Cur%Previous)) then
          Do While(Associated(LL_Cur%Previous))
            LL_Cur => LL_Cur%Previous
          End Do
        EndIf
        Do While (Associated(LL_Cur))
          If(Allocated(LL_Cur%Val)) Allocate(Output,Source=LL_Cur%Val)
          Select Type(Output)
          Type is (Integer)
            Write(IOut,'(I3)',advance="no") Output
          Type is (Character(*))
            Write(IOut,'(A10)',advance="no") Output
          Type is (Real)
            Write(IOut,'(D15.5)',advance="no") Output
          Type is (MQC_Bits)
            Call Output%Print(IOut)
          Class Default
            Stop 'Unexpected data type in link'
          End Select
          Deallocate(Output)
          LL_Cur => LL_Cur%Next
          If(Associated(LL_Cur)) Write(IOut,'(A4)',advance="no") " <-> " 
        EndDo
      EndIf
      write(IOut,*)
!
      Return
      End Subroutine LinkedList_Print
!
!     PROCEDURE LinkedList_Insert_Node
      Subroutine LinkedList_Insert_Node(LL_Cur,New_Value,Direction)
!
!     This routine carries out a insertion of a particular node either after
!     (default) or before (specified by Direction) the current node LL_Cur.
!
!     L. M. Thompson
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Character(Len=*),Optional,Intent(In)::Direction
      Character(Len=32)::My_Direction
      Type(MQC_LinkedList),Pointer::Cur
      Class(*),Intent(In)::New_Value
!
      If(Present(Direction)) then
        Call String_Change_Case(Direction,'l',My_Direction)
      Else
        My_Direction = 'forward'
      EndIf
!
      If(Associated(LL_Cur)) then
        Allocate(Cur)
        Select Case (My_Direction)
        Case('forward')
          If(Associated(LL_Cur%Next)) then
            Cur%Next => LL_Cur%Next
            LL_Cur%Next%Previous => Cur
            LL_Cur%Next => Cur
            Cur%Previous => LL_Cur
          Else
            Cur%Next => Null()
            LL_Cur%Next => Cur
            Cur%Previous => LL_Cur
          EndIf
        Case('backward')
          If(Associated(LL_Cur%Previous)) then
            Cur%Previous => LL_Cur%Previous
            LL_Cur%Previous%Next => Cur
            LL_Cur%Previous => Cur
            Cur%Next => LL_Cur
          Else
            Cur%Previous => Null()
            LL_Cur%Previous => Cur
            Cur%Next => LL_Cur
          EndIf
        Case Default
          Call MQC_Error_A('Invalid direction sent to LinkedList_Insert_Node',6, &
            'Direction:',My_Direction)
        End Select
        Allocate(Cur%Val,Source=New_Value)
      EndIf
!
      Return
      End Subroutine LinkedList_Insert_Node  
!
!     PROCEDURE LinkedList_Push
      Subroutine LinkedList_Push(LL_Head,New_Value)
!
!     This routine carries out a push operation (append at end) on a linked
!     list.
!
!     -H. P. Hratchian, 2015
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Head
      Class(*),Intent(In)::New_Value
      Type(MQC_LinkedList),Pointer::Cur
!
!     Go through the linked list starting at LL_Head until we find the end,
!     then insert the new node at the end of the list.
!
      If(Associated(LL_Head)) then
        Cur => LL_Head
        Do While(Associated(Cur%Next))
          Cur => Cur%Next
        End Do
        Allocate(Cur%Next)
        Allocate(Cur%Next%Previous)
        Cur%Next%Previous => Cur
        Cur => Cur%Next
        Nullify(Cur%Next)
        Allocate(Cur%Val,Source=New_Value)
      else
        Allocate(LL_Head)
        Nullify(LL_Head%Next)
        Nullify(LL_Head%Previous)
        Allocate(LL_Head%Val,Source=New_Value)
      endIf
!
      Return
      End Subroutine LinkedList_Push
!
!     PROCEDURE LinkedList_Unshift
      Subroutine LinkedList_UnShift(LL_Head,New_Value)
!
!     This routine carries out an unshift operation (append at front) on a
!     linked list.
!
!     -H. P. Hratchian, 2015
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Head
      Class(*),Intent(In)::New_Value
      Type(MQC_LinkedList),Pointer::LL_Temp
!
      If(Associated(LL_Head)) then
        Do While(Associated(LL_Head%Previous))
          LL_Head => LL_Head%Previous
        End Do
        Allocate(LL_Temp)
        LL_Temp%Next => LL_Head
        LL_Head%Previous => LL_Temp
        Allocate(LL_Temp%Val,Source=New_Value)
        LL_Head => LL_Temp
      else
        Allocate(LL_Head)
        Nullify(LL_Head%Next)
        Nullify(LL_Head%Previous)
        Allocate(LL_Head%Val,Source=New_Value)
      endIf
!
      Return
      End Subroutine LinkedList_Unshift
!
!     PROCEDURE LinkedList_GoToBeginning
      Subroutine LinkedList_GoToBeginning(Current)
!
!     This routine moves the argument Current to the head of the linked
!     list.
!
!     L. M. Thompson, 2018
!
      Implicit None
      Type(MQC_LinkedList),Pointer,Intent(InOut)::Current
!
      Do While(Associated(Current%Previous))
          Current => Current%Previous
      End Do

      Return
      End Subroutine LinkedList_GoToBeginning
!
!     PROCEDURE LinkedList_GoToEnd
      Subroutine LinkedList_GoToEnd(Current)
!
!     This routine moves the argument Current to the end of the linked
!     list.
!
!     L. M. Thompson, 2018
!
      Implicit None
      Type(MQC_LinkedList),Pointer,Intent(InOut)::Current
!
      Do While(Associated(Current%Next))
          Current => Current%Next
      End Do

      Return
      End Subroutine LinkedList_GoToEnd
!
!     PROCEDURE LinkedList_GetNext
      Subroutine LinkedList_GetNext(Current,Last,Last_Looks_Ahead)
!
!     This routine moves the argument Current forward to the next member of
!     the link list to which Current belongs. If optional logical dummy
!     argument Last_Looks_Ahead, the logical argument Last is returned as
!     TRUE if Current is the last member of the list on return; otherwise,
!     Last is sent as TRUE if the requested next member of the linked list
!     does not exist because the end of the list has been reached.
!     Regardless of Last_Looks_Ahead, if the last member of the list is
!     sent as input Current then Last is always set to TRUE and this
!     routine returns to the calling procedure without updating Current.
!
!     If Last_Looks_Ahead is not sent, it is taken to be TRUE by default.
!
!     -H. P. Hratchian, 2015
!
      Implicit None
      Type(MQC_LinkedList),Pointer,Intent(InOut)::Current
      Logical,Intent(Out)::Last
      Logical,Optional,Intent(In)::Last_Looks_Ahead
      Logical::My_Last_Looks_Ahead
!
      If(Present(Last_Looks_Ahead)) then
        My_Last_Looks_Ahead = Last_Looks_Ahead
      else
        My_Last_Looks_Ahead = .True.
      endIf
      If(.not.Associated(Current%Next)) then
        Last = .True.
      else
        Current => Current%Next
        If(My_Last_Looks_Ahead) then
          Last = .not.Associated(Current%Next)
        else
          Last = .False.
        endIf
      endIf
      Return
      End Subroutine LinkedList_GetNext
!
!     PROCEDURE LinkedList_GetPrevious
      Subroutine LinkedList_GetPrevious(Current,First,First_Looks_Back)
!
!     This routine moves the argument Current backwards to the previous member of
!     the link list to which Current belongs. If optional logical dummy
!     argument First_Looks_Back, the logical argument First is returned as
!     TRUE if Current is the first member of the list on return; otherwise,
!     First is sent as TRUE if the requested previous member of the linked list
!     does not exist because the beginning of the list has been reached.
!     Regardless of First_Looks_Back, if the first member of the list is
!     sent as input Current then First is always set to TRUE and this
!     routine returns to the calling procedure without updating Current.
!
!     If First_Looks_Back is not sent, it is taken to be TRUE by default.
!
!     L. M. Thompson, 2018
!
      Implicit None
      Type(MQC_LinkedList),Pointer,Intent(InOut)::Current
      Logical,Intent(Out)::First
      Logical,Optional,Intent(In)::First_Looks_Back
      Logical::My_First_Looks_Back
!
      If(Present(First_Looks_Back)) then
        My_First_Looks_Back = First_Looks_Back
      else
        My_First_Looks_Back = .True.
      endIf
      If(.not.Associated(Current%Previous)) then
        First = .True.
      else
        Current => Current%Previous
        If(My_First_Looks_Back) then
          First = .not.Associated(Current%Previous)
        else
          First = .False.
        endIf
      endIf
      Return
      End Subroutine LinkedList_GetPrevious
!
!     PROCEDURE LinkedList_IsPreviousAssociated
      Function LinkedList_IsPreviousAssociated(Current)
!
!     This logical function returns TRUE is the LinkedList member Current
!     points to a previous list member; otherwise, FALSE is returned.
!
!     Lee M. Thompson, 2018
!
      Implicit None
      Logical::LinkedList_IsPreviousAssociated
      Type(MQC_LinkedList),Pointer,Intent(In)::Current
!
      LinkedList_IsPreviousAssociated = Associated(Current%Previous)
!
      Return
      End Function LinkedList_IsPreviousAssociated
!
!     PROCEDURE LinkedList_IsNextAssociated
      Function LinkedList_IsNextAssociated(Current)
!
!     This logical function returns TRUE is the LinkedList member Current
!     points to a next list member; otherwise, FALSE is returned.
!
!     -H. P. Hratchian, 2015
!
      Implicit None
      Logical::LinkedList_IsNextAssociated
      Type(MQC_LinkedList),Pointer,Intent(In)::Current
!
      LinkedList_IsNextAssociated = Associated(Current%Next)
!
      Return
      End Function LinkedList_IsNextAssociated
!
!     PROCEDURE LinkedList_Delete_Node
      Subroutine LinkedList_Delete_Node(LL_Cur)
!
!     This routine carries out a deletion of the particular node and links
!     previous and next nodes together. Puts us on the next node unless next
!     is null, in which case it puts us on the last node.
!
!     L. M. Thompson, 2018
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Type(MQC_LinkedList),Pointer::Cur
      Class(*),Allocatable::Output
!
!     Delete a linked list starting after LL_Cur until we reach the end.
!
      If(Associated(LL_Cur)) then
        Cur => LL_Cur
        If(Associated(Cur%Next).and.Associated(Cur%Previous)) then
          Cur%Next%Previous => Cur%Previous
          Cur%Previous%Next => Cur%Next
          LL_Cur => Cur%Next
          Nullify(Cur%Next)
          Nullify(Cur%Previous)
          Deallocate(Cur%Val)
          Nullify(Cur)
        ElseIf(Associated(Cur%Next).and..not.Associated(Cur%Previous)) then
          Cur%Next%Previous => Null()
          LL_Cur => Cur%Next
          Nullify(Cur%Next)
          Deallocate(Cur%Val)
          Nullify(Cur)
        ElseIf(.not.Associated(Cur%Next).and.Associated(Cur%Previous)) then
          Cur%Previous%Next => Null()
          LL_Cur => Cur%Previous
          Nullify(Cur%Previous)
          Deallocate(Cur%Val)
          Nullify(Cur)
        ElseIf(.not.Associated(Cur%Next).and..not.Associated(Cur%Previous)) then
          Nullify(Cur%Previous)
          Deallocate(Cur%Val)
          Nullify(Cur)
          Nullify(LL_Cur)
        EndIf
      EndIf
!
      Return
      End Subroutine LinkedList_Delete_Node  
!
!     PROCEDURE LinkedList_Delete_Branch
      Subroutine LinkedList_Delete_Branch(LL_Head,Direction)
!
!     This routine carries out a deletion of a links after (default) or before
!     the current node.
!
!     L. M. Thompson, 2018
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Head
      Type(MQC_LinkedList),Pointer::Cur
      Character(Len=*),Optional,Intent(In)::Direction
      Character(Len=32)::My_Direction
!
      If(Present(Direction)) then
        Call String_Change_Case(Direction,'l',My_Direction)
      Else
        My_Direction = 'forward'
      EndIf
!
      If(Associated(LL_Head)) then
        Cur => LL_Head
        Select Case (My_Direction)
        Case('forward')
          Do While(Associated(Cur%Next))
            Cur => Cur%Next
            LL_Head%Next => Cur%Next
            Nullify(Cur%Next)
            Nullify(Cur%Previous)
            Deallocate(Cur%Val)
            Nullify(Cur)
            Cur => LL_Head
          End Do
        Case('backward')
          Do While(Associated(Cur%Previous))
            Cur => Cur%Previous
            LL_Head%Previous => Cur%Previous
            Nullify(Cur%Next)
            Nullify(Cur%Previous)
            Deallocate(Cur%Val)
            Nullify(Cur)
            Cur => LL_Head
          End Do
        Case Default
          Call MQC_Error_A('Invalid direction sent to LinkedList_Delete_Branch',6, &
            'Direction:',My_Direction)
        End Select
      EndIf
!
      Return
      End Subroutine LinkedList_Delete_Branch
!
!     PROCEDURE LinkedList_Delete
      Subroutine LinkedList_Delete(LL_Head)
!
!     This routine carries out a deletion of a linked list.
!
!     L. M. Thompson, 2018
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Head,Cur,tmp
!
!     Delete a linked list starting at LL_Head until we reach the end.
!
      If(Associated(LL_Head)) then
        Cur => LL_Head
        Do While(Associated(Cur%Next))
          Cur => Cur%Next
        End Do
        Do While(Associated(Cur%Previous))
          Deallocate(Cur%Val)
          Cur => Cur%Previous
          Nullify(Cur%Next%Next)
          Nullify(Cur%Next%Previous)
        EndDo
        If(Allocated(Cur%Val)) Deallocate(Cur%Val)
        Nullify(Cur%Next)
        Nullify(Cur%Previous)
        Nullify(Cur)
        Nullify(LL_Head)
      EndIf
!
      Return
      End Subroutine LinkedList_Delete
!
      End Module MQC_DataStructures
