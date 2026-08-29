      program unitTest02
!
!     This unit test validates character scalar and vector support in
!     MQC_Variable, including intrinsic-like inquiry, assignment, ADJUSTL,
!     ADJUSTR, LEN_TRIM, and scalar TRIM operations, printing, and
!     String_Change_Case delegation to MQC_General.
!
!     The internal assert* routines provide a lightweight, self-contained
!     validation framework. Successful comparisons return normally; failed
!     comparisons print the test label and available diagnostics, then use
!     ERROR STOP 1 to terminate the program with a nonzero status. Reaching the
!     final PASS message therefore means that every assertion succeeded.
!
!     H. P. Hratchian, 2026.
!
!
      use MQC_General
      use MQC_Algebra2
      use iso_fortran_env,only:int64
      implicit none
!
      type(MQC_Variable)::scalarVariable,vectorVariable,changedVariable
      type(MQC_Variable)::elementVariable,reshapedVariable,clearedVariable
      type(MQC_Variable)::adjustedVariable,lengthVariable,trimmedVariable
      type(MQC_Variable)::scalarKernelVariable
      character(len=8),dimension(3)::vectorInput
      character(len=8),dimension(3)::intrinsicCharacterVector
      character(len=:),dimension(:),allocatable::vectorOutput
      character(len=8),dimension(3)::vectorOutputFixed
      character(len=8),dimension(3)::expectedCharacterVector
      character(len=:),allocatable::trimmedOutput
      character(len=8)::intrinsicCharacterScalar
      character(len=8)::scalarOutput
      character(len=256)::line
      integer(kind=int64)::iOut
      integer(kind=int64)::integerScalar
      integer(kind=int64),dimension(:),allocatable::integerVector
!
!
!     Character scalar assignment and intrinsic-like inquiries.
!
      scalarVariable = 'Ab c9'
      call assertCharacter('scalar type',TRIM(scalarVariable%getType()),  &
        'CHARACTER')
      call assertInteger('scalar rank',RANK(scalarVariable),0_int64)
      call assertInteger('scalar size',SIZE(scalarVariable),1_int64)
      call assertInteger('scalar length',LEN(scalarVariable),5_int64)
      scalarOutput = scalarVariable
      call assertCharacter('scalar round trip',scalarOutput,'Ab c9')
      scalarOutput = 'MiXeD'
      call String_Change_Case(scalarOutput,'L')
      call assertCharacter('intrinsic string change case',scalarOutput,  &
        'mixed')
!
!     Both the generic and type-bound change-case entry points must delegate
!     to the intrinsic-character implementation without changing element width.
!
      call String_Change_Case(scalarVariable,'U',changedVariable)
      scalarOutput = scalarVariable
      call assertCharacter('change-case input retained',scalarOutput,'Ab c9')
      scalarOutput = changedVariable
      call assertCharacter('generic uppercase',scalarOutput,'AB C9')
      call changedVariable%change_case('L')
      scalarOutput = changedVariable
      call assertCharacter('type-bound lowercase',scalarOutput,'ab c9')
      call assertInteger('case preserves length',LEN(changedVariable),5_int64)
      scalarVariable = 7_int64
      scalarVariable = 'ReTyped'
      scalarOutput = scalarVariable
      call assertCharacter('numeric to character reassignment',scalarOutput,  &
        'ReTyped')
!
!     Intrinsic-character kernels in MQC_General and the corresponding
!     overloaded and type-bound MQC_Variable operations.
!
      intrinsicCharacterScalar = '  Ab'
      intrinsicCharacterVector(1) = '  One'
      intrinsicCharacterVector(2) = ' Two'
      intrinsicCharacterVector(3) = 'Three'
!
      scalarOutput = mqc_adjustl(intrinsicCharacterScalar)
      call assertCharacter('general scalar adjustl',scalarOutput,'Ab')
      scalarOutput = mqc_adjustr(intrinsicCharacterScalar)
      call assertCharacter('general scalar adjustr',scalarOutput,'      Ab')
      integerScalar = mqc_len_trim(intrinsicCharacterScalar)
      call assertInteger('general scalar len_trim',integerScalar,4_int64)
      trimmedOutput = mqc_trim(intrinsicCharacterScalar)
      call assertCharacter('general scalar trim',trimmedOutput,'  Ab')
      call assertInteger('general scalar trim length',  &
        LEN(trimmedOutput),4_int64)
!
      call testGeneralCharacterVectors(intrinsicCharacterVector)
!
      scalarKernelVariable = '  Ab    '
      adjustedVariable = ADJUSTL(scalarKernelVariable)
      scalarOutput = adjustedVariable
      call assertCharacter('overloaded scalar adjustl',scalarOutput,'Ab')
      call assertInteger('adjustl preserves scalar length',  &
        LEN(adjustedVariable),8_int64)
      adjustedVariable = scalarKernelVariable%adjustr()
      scalarOutput = adjustedVariable
      call assertCharacter('type-bound scalar adjustr',scalarOutput,  &
        '      Ab')
      lengthVariable = LEN_TRIM(scalarKernelVariable)
      integerScalar = lengthVariable
      call assertInteger('overloaded scalar len_trim',integerScalar,4_int64)
      trimmedVariable = scalarKernelVariable%trim()
      scalarOutput = trimmedVariable
      call assertCharacter('type-bound scalar trim',scalarOutput,'  Ab')
      call assertInteger('trim changes scalar length',  &
        LEN(trimmedVariable),4_int64)
!
      vectorVariable = intrinsicCharacterVector
      adjustedVariable = vectorVariable%adjustl()
      vectorOutput = adjustedVariable
      vectorOutputFixed = vectorOutput
      expectedCharacterVector(1) = 'One'
      expectedCharacterVector(2) = 'Two'
      expectedCharacterVector(3) = 'Three'
      call assertTrue('type-bound vector adjustl',  &
        ALL(vectorOutputFixed.eq.expectedCharacterVector))
      adjustedVariable = ADJUSTR(vectorVariable)
      vectorOutput = adjustedVariable
      vectorOutputFixed = vectorOutput
      expectedCharacterVector(1) = '     One'
      expectedCharacterVector(2) = '     Two'
      expectedCharacterVector(3) = '   Three'
      call assertTrue('overloaded vector adjustr',  &
        ALL(vectorOutputFixed.eq.expectedCharacterVector))
      lengthVariable = vectorVariable%len_trim()
      integerVector = [0_int64]
      integerVector = lengthVariable
      call assertIntegerVector('type-bound vector len_trim',integerVector,  &
        [5_int64,4_int64,5_int64])
!
!     Character vector assignment, extraction, insertion, reshape, and
!     allocatable intrinsic assignment.
!
      vectorInput(1) = 'Alpha'
      vectorInput(2) = 'bETA'
      vectorInput(3) = 'Gamma'
      vectorVariable = vectorInput
      call assertInteger('vector rank',RANK(vectorVariable),1_int64)
      call assertInteger('vector size',SIZE(vectorVariable),3_int64)
      call assertInteger('vector dimension',SIZE(vectorVariable,1),3_int64)
      call assertInteger('vector length',LEN(vectorVariable),8_int64)
      vectorOutput = vectorVariable
      vectorOutputFixed = vectorOutput
      call assertTrue('vector round trip',  &
        ALL(vectorOutputFixed.eq.vectorInput))
!
      call vectorVariable%put('Delta',[2_int64])
      elementVariable = vectorVariable%getVal([2_int64])
      scalarOutput = elementVariable
      call assertCharacter('put and getVal',scalarOutput,'Delta')
!
      reshapedVariable = RESHAPE(vectorVariable,[3_int64])
      vectorOutput = reshapedVariable
      vectorOutputFixed = vectorOutput
      expectedCharacterVector(1) = 'Alpha'
      expectedCharacterVector(2) = 'Delta'
      expectedCharacterVector(3) = 'Gamma'
      call assertTrue('reshape preserves values',  &
        ALL(vectorOutputFixed.eq.expectedCharacterVector))
      changedVariable = MQC_Variable_mqc2mqc(vectorVariable)
      vectorOutput = changedVariable
      vectorOutputFixed = vectorOutput
      call assertTrue('explicit MQC copy',  &
        ALL(vectorOutputFixed.eq.expectedCharacterVector))
!
!     Clear must retain the fixed element width across the new vector.
!
      call clearedVariable%clear('MiXeD',[3_int64])
      call assertInteger('clear length',LEN(clearedVariable),5_int64)
      call clearedVariable%change_case('U')
      vectorOutput = clearedVariable
      vectorOutputFixed = vectorOutput
      expectedCharacterVector = 'MIXED'
      call assertTrue('clear and change case',  &
        ALL(vectorOutputFixed.eq.expectedCharacterVector))
!
!     The scalar MQC conversion and the character printing wrappers are also
!     part of the public character pathway.
!
      scalarVariable = MQC('Print Me')
      scalarOutput = scalarVariable
      call assertCharacter('mqc character conversion',scalarOutput,'Print Me')
      open(newunit=iOut,status='scratch',action='readwrite',form='formatted')
      call mqc_print_scalar('General',iOut,header='scalar')
      call mqc_print_vector([ character(len=4)::'One','Two' ],iOut,  &
        header='vector')
      call mqc_print([.true.,.false.],iOut,header='logical mqc_print')
      call mqc_print_vector([.false.,.true.],iOut,  &
        header='logical mqc_print_vector')
      rewind(iOut)
      read(iOut,'(A)') line
      call assertContains('general scalar print',line,'General')
      read(iOut,'(A)') line
      call assertContains('general vector print header',line,'vector')
      read(iOut,'(A)') line
      call assertContains('general vector print element',line,'One')
      read(iOut,'(A)') line
      call assertContains('general vector print second element',line,'Two')
      read(iOut,'(A)') line
      call assertContains('logical mqc_print header',line,  &
        'logical mqc_print')
      read(iOut,'(A)') line
      call assertContains('logical mqc_print first element',line,'T')
      read(iOut,'(A)') line
      call assertContains('logical mqc_print second element',line,'F')
      read(iOut,'(A)') line
      call assertContains('logical mqc_print_vector header',line,  &
        'logical mqc_print_vector')
      read(iOut,'(A)') line
      call assertContains('logical mqc_print_vector first element',line,'F')
      read(iOut,'(A)') line
      call assertContains('logical mqc_print_vector second element',line,'T')
      close(iOut)
!
      open(newunit=iOut,status='scratch',action='readwrite',form='formatted')
      call scalarVariable%print(iOut=iOut,header='scalar')
      call vectorVariable%print(iOut=iOut,header='vector')
      rewind(iOut)
      read(iOut,'(A)') line
      call assertContains('scalar print',line,'Print Me')
      read(iOut,'(A)') line
      call assertContains('vector print header',line,'vector')
      read(iOut,'(A)') line
      call assertContains('vector print element',line,'Alpha')
      close(iOut)
!
      write(*,'(1x,A)') 'unitTest02: PASS'
!
      contains


!PROCEDURE testGeneralCharacterVectors
      subroutine testGeneralCharacterVectors(input)
      implicit none
      character(len=8),dimension(3),intent(in)::input
      character(len=8),dimension(3)::actual,expected
      integer(kind=int64),dimension(:),allocatable::lengths
!
      actual = mqc_adjustl(input)
      expected(1) = 'One'
      expected(2) = 'Two'
      expected(3) = 'Three'
      call assertTrue('general vector adjustl',ALL(actual.eq.expected))
!
      actual = mqc_adjustr(input)
      expected(1) = '     One'
      expected(2) = '     Two'
      expected(3) = '   Three'
      call assertTrue('general vector adjustr',ALL(actual.eq.expected))
!
      lengths = mqc_len_trim(input)
      call assertIntegerVector('general vector len_trim',lengths,  &
        [5_int64,4_int64,5_int64])
!
      return
      end subroutine testGeneralCharacterVectors

!
!PROCEDURE assertCharacter
      subroutine assertCharacter(label,actual,expected)
      implicit none
      character(len=*),intent(in)::label,actual,expected
!
      if(actual.ne.expected) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(label)
        write(*,'(3x,A,A,A)') 'actual   = <',actual,'>'
        write(*,'(3x,A,A,A)') 'expected = <',expected,'>'
        error stop 1
      endIf
!
      return
      end subroutine assertCharacter


!
!
!PROCEDURE assertContains
      subroutine assertContains(label,actual,expectedSubstring)
      implicit none
      character(len=*),intent(in)::label,actual,expectedSubstring
!
      if(INDEX(actual,expectedSubstring).eq.0) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(label)
        write(*,'(3x,A,A,A)') 'line = <',TRIM(actual),'>'
        error stop 1
      endIf
!
      return
      end subroutine assertContains


!
!PROCEDURE assertInteger
      subroutine assertInteger(label,actual,expected)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::actual,expected
!
      if(actual.ne.expected) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(label)
        write(*,'(3x,A,I0)') 'actual   = ',actual
        write(*,'(3x,A,I0)') 'expected = ',expected
        error stop 1
      endIf
!
      return
      end subroutine assertInteger


!
!PROCEDURE assertIntegerVector
      subroutine assertIntegerVector(label,actual,expected)
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),dimension(:),intent(in)::actual,expected
!
      if(SIZE(actual).ne.SIZE(expected)) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(label)//' size mismatch'
        error stop 1
      endIf
      if(ANY(actual.ne.expected)) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(label)
        write(*,'(3x,A,*(I0,1x))') 'actual   = ',actual
        write(*,'(3x,A,*(I0,1x))') 'expected = ',expected
        error stop 1
      endIf
!
      return
      end subroutine assertIntegerVector


!
!PROCEDURE assertTrue
      subroutine assertTrue(label,condition)
      implicit none
      character(len=*),intent(in)::label
      logical,intent(in)::condition
!
      if(.not.condition) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(label)
        error stop 1
      endIf
!
      return
      end subroutine assertTrue


      end program unitTest02
