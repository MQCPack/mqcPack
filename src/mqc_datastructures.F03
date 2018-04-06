      Module MQC_DataStructures
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
      use iso_fortran_env
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
      End Interface


      CONTAINS
!
!
!PROCEDURE LinkedList_Return_Integer_Value
      Subroutine LinkedList_Return_Integer_Value(LL_Cur,Return_Value)
!
!     This routine is used to return the value of a member of a linked
!     list. LL_Cur is an INPUT linked list member and Return_Value is an
!     OUTPUT argument that is allocatable of the same type as LL_Cur%Val.
!
!     -H. P. Hratchian, 2016
!
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Integer(kind=int64),Allocatable,Intent(InOut)::Return_Value
!
!     Allocate Return_Value, set it to LL_Cur%Val, and return to the
!     calling procedure.
!
!
      If(Associated(LL_Cur)) then
        If(Allocated(LL_Cur%Val)) then
          If(.not.Allocated(Return_Value)) Allocate(Return_Value)
          Select Type (Val =>LL_Cur%Val)
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
!
!PROCEDURE LinkedList_Return_Real_Value
      Subroutine LinkedList_Return_Real_Value(LL_Cur,Return_Value)
!
!     This routine is used to return the value of a member of a linked
!     list. LL_Cur is an INPUT linked list member and Return_Value is an
!     OUTPUT argument that is allocatable of the same type as LL_Cur%Val.
!
!     -H. P. Hratchian, 2016
!
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Real(kind=real64),Allocatable,Intent(InOut)::Return_Value
!
!     Allocate Return_Value, set it to LL_Cur%Val, and return to the
!     calling procedure.
!
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
!
!PROCEDURE LinkedList_Return_Character_Value
      Subroutine LinkedList_Return_Character_Value(LL_Cur,Return_Value)
!
!     This routine is used to return the value of a member of a linked
!     list. LL_Cur is an INPUT linked list member and Return_Value is an
!     OUTPUT argument that is allocatable of the same type as LL_Cur%Val.
!
!     -H. P. Hratchian, 2016
!
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Character,Allocatable,Intent(InOut)::Return_Value
!
!     Allocate Return_Value, set it to LL_Cur%Val, and return to the
!     calling procedure.
!
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
!PROCEDURE LinkedList_Print
      Subroutine LinkedList_Print(LL_Head,IOut)
!
!     This routine prints the whole list. The head of the list is required.
!
      Implicit None
      Integer(kind=int64)::IOut
      Class(*),Allocatable::Val
      Type(MQC_LinkedList),Pointer::LL_Head, LL_Cur
      LL_Cur => LL_Head
      Do While (Associated(LL_Cur))
        If(Allocated(LL_Cur%Val)) Allocate(Val,Source=LL_Cur%Val)
        Select Type(Val)
        Type is (Integer)
          Write(IOut,'(I3,A4)',advance="no") Val," -> "
        Type is (Character(*))
          Write(IOut,'(A10,A4)',advance="no") Val," -> "
        Type is (Real)
          Write(IOut,'(D15.5,A4)',advance="no") Val," -> "
        Class Default
          Stop 'Unexpected data type in link'
        End Select
        LL_Cur => LL_Cur%Next
      EndDo
!
!
      Return
      End Subroutine LinkedList_Print
!
!
!PROCEDURE LinkedList_Push
      Subroutine LinkedList_Push(LL_Head,New_Value)
!
!     This routine carries out a push operation (append at end) on a linked
!     list.
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Head
      Class(*),Intent(In)::New_Value
      Type(MQC_LinkedList),Pointer::Cur
!
!     Go through the linked list starting at LL_Head until we find the end,
!     then insert the new node at the end of the list.
!
!
!     -H. P. Hratchian, 2015
!
!
      If(Associated(LL_Head)) then
        Cur => LL_Head
        Do While(Associated(Cur%Next))
          Cur => Cur%Next
        End Do
        Allocate(Cur%Next)
        Cur => Cur%Next
        Nullify(Cur%Next)
        Allocate(Cur%Val,Source=New_Value)
      else
        Allocate(LL_Head)
        Nullify(LL_Head%Next)
        Allocate(LL_Head%Val,Source=New_Value)
      endIf
!
      Return
      End Subroutine LinkedList_Push
!
!
!PROCEDURE LinkedList_Unshift
      Subroutine LinkedList_UnShift(LL_Head,New_Value)
!
!     This routine carries out an unshift operation (append at front) on a
!     linked list.
!
!
!     -H. P. Hratchian, 2015
!
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Head
      Class(*),Intent(In)::New_Value
      Type(MQC_LinkedList),Pointer::LL_Temp
!
      If(Associated(LL_Head)) then
        Allocate(LL_Temp)
        LL_Temp%Next => LL_Head
        Allocate(LL_Temp%Val,Source=New_Value)
        LL_Head => LL_Temp
      else
        Allocate(LL_Head)
        Nullify(LL_Head%Next)
        Allocate(LL_Head%Val,Source=New_Value)
      endIf
!
      Return
      End Subroutine LinkedList_Unshift
!
!
!PROCEDURE LinkedList_GetNext
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
!
!     -H. P. Hratchian, 2015
!
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
!
!PROCEDURE LinkedList_GetPrevious
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
!
!     L. M. Thompson, 2018
!
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
!
!PROCEDURE LinkedList_IsNextAssociated
      Function LinkedList_IsNextAssociated(Current)
!
!     This logical function returns TRUE is the LinkedList member Current
!     points to a next list member; otherwise, FALSE is returned.
!
!
!     -H. P. Hratchian, 2015
!
!
      Implicit None
      Logical::LinkedList_IsNextAssociated
      Type(MQC_LinkedList),Pointer,Intent(In)::Current
!
      LinkedList_IsNextAssociated = Associated(Current%Next)
!
      Return
      End Function LinkedList_IsNextAssociated


!hph+
!!
!!PROCEDURE LinkedList_Get_Value
!      Subroutine LinkedList_Get_Value(Val,CoordOut,OK)
!!
!!     This routine returns a coordinate type (CoordOut) given an input
!!     coordinate type given by Val.
!!
!!
!!     -H. P. Hratchian, 2015
!!
!!
!      Implicit None
!      Class(*),Intent(In)::Val
!      Class(*),Intent(Out)::CoordOut
!      Logical,Intent(Out)::OK
!!
!      OK = .False.
!      Select Type(Val)
!      Type is (Bond_Type)
!        Select Type(CoordOut)
!        Type is (Bond_Type)
!          CoordOut = Val
!          OK = .True.
!        End Select
!      Type is (Angle_Type)
!        Select Type(CoordOut)
!        Type is (Angle_Type)
!          CoordOut = Val
!          OK = .True.
!        End Select
!      Type is (Dihedral_Type)
!        Select Type(CoordOut)
!        Type is (Dihedral_Type)
!          CoordOut = Val
!          OK = .True.
!        End Select
!      Type is (Bond_Combination_Type)
!        Select Type(CoordOut)
!        Type is (Bond_Combination_Type)
!          CoordOut = Val
!          OK = .True.
!        End Select
!      Type is (Dihedral_Combination_Type)
!        Select Type(CoordOut)
!        Type is (Dihedral_Combination_Type)
!          CoordOut = Val
!          OK = .True.
!        End Select
!      End Select
!!
!      Return
!      End Subroutine LinkedList_Get_Value
!hph-
!
!
!!
!!PROCEDURE LinkedList_Print
!      Subroutine LinkedList_Print(IOut,Current_Link)
!!
!!     Quick way to print all values in a linked list
!!
!      Implicit None
!      Integer::IOut
!      Type(MQC_LinkedList)::Current_Link
!      Type(MQC_LinkedList),Pointer::Looped_Link
!
!      Call LinkedList_PrintValue(IOut,Current_Link)
!      Looped_Link => Current_Link%Next
!      Do While(Associated(Looped_Link))
!        Call LinkedList_PrintValue(IOut,Looped_Link)
!        Looped_Link => Looped_Link%Next
!      EndDo
!
!      End Subroutine LinkedList_Print
!
!
!!
!!PROCEDURE LinkedList_PrintValue
!      Subroutine LinkedList_PrintValue(IOut,Current_Link)
!!
!!     Print value at link in linked list
!!
!      Implicit None
!      End Type Node
!      Integer::IOut
!      Type(MQC_LinkedList)::Current_Link
!
!      Select Type(Value => Current_Link%Val)
!      Type is (Integer)
!        Write(IOut,*) Value
!      Type is (Character(*))
!        Write(IOut,*) Value(1:1)
!      Type is (Real)
!        Write(IOut,*) Value
!      Class Default
!        Stop 'Unexpected data type in link'
!      End Select
!
!      End Subroutine LinkedList_PrintValue

!
!
!PROCEDURE LinkedList_Delete_Branch
      Subroutine LinkedList_Delete_Branch(LL_Cur)
!
!     This routine carries out a deletion of a links after current node.
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Cur
      Type(MQC_LinkedList),Pointer::Cur
!
!     Delete a linked list starting after LL_Cur until we reach the end.
!
!
      If(Associated(LL_Cur)) then
        Cur => LL_Cur
        Do While(Associated(Cur%Next))
          Cur => Cur%Next
        End Do
        Do While(Cur.ne.LL_Cur)
          Deallocate(Cur%Val)
          Nullify(Cur%Next)
          Cur => Cur%Previous
        EndDo
      EndIf
!
      Return
      End Subroutine LinkedList_Delete_Branch

!
!PROCEDURE LinkedList_Delete
      Subroutine LinkedList_Delete(LL_Head)
!
!     This routine carries out a deletion of a linked list.
!
      Implicit None
      Type(MQC_LinkedList),Pointer::LL_Head
      Type(MQC_LinkedList),Pointer::Cur
!
!     Delete a linked list starting at LL_Head until we reach the end.
!
!
      If(Associated(LL_Head)) then
        Cur => LL_Head
        Do While(Associated(Cur%Next))
          Cur => Cur%Next
        End Do
        Do While(Associated(Cur%Previous))
          Deallocate(Cur%Val)
          Nullify(Cur%Next)
          Cur => Cur%Previous
        EndDo
        Deallocate(Cur%Val)
        Nullify(Cur%Next)
        Nullify(Cur)
        Nullify(LL_Head)
      EndIf
!
      Return
      End Subroutine LinkedList_Delete


      ! Type MQC_BinaryTree
      !   Private
      !     Class(*),Allocatable,Public::Val
      !     Type(MQC_BinaryTree),Pointer::Left
      !     Type(MQC_BinaryTree),Pointer::Right
      ! End Type MQC_BinaryTree
      !
      ! Subroutine BinaryTree_Insert(Node_Previous,New_Value,Left)
      !
      ! ! This subroutine inserts a node in a tree after a given node.
      ! ! Arguments:
      ! !     Node_Previous       the node to be inserted after.
      ! !     New_Value           the value to be inserted, can be any type.
      ! !     Left                True if new node is to be inserted as the
      ! !                         left child, vise versa.
      !
      !   Implicit None
      !   Type(MQC_BinaryTree),Pointer::Node_Previous,New_Node
      !   logical::Left
      !   Class(*),Intent(In)::New_Value
      !
      !   Allocate(New_Node%Val, source = New_Value)
      !   Nullify(New_Node%Left,New_Node%Right)
      !   if (.not. Associated(Node_Previous)) then
      !     Node_Previous => New_Node
      !   else
      !     if (Left) then
      !         New_Node%Left => Node_Previous%Left
      !         Node_Previous%Left => New_Node
      !     else
      !         New_Node%right => Node_Previous%right
      !         Node_Previous%right => New_Node
      !     endif
      !   endif
      ! End Subroutine BinaryTree_Insert

      ! Type MQC_Graph_AdjList
      !   Private
      !     Class(*)::Val
      !     Type(MQC_LinkedList)::adjacent
      ! End Type
      !
      ! Type MQC_ReactionMap
      !   Private
      !     Type(MQC_Graph_AdjList),dimension(:),allocatable::Map
      !     integer::capacity,size
      ! End Type
      !
      ! Subroutine ReactionMap_Insert(New_Node)
      !   Implicit None
      !   Type(MQC_Graph_AdjList)::New_Node
      !
      !
      !
      ! End subroutine



      End Module MQC_DataStructures
