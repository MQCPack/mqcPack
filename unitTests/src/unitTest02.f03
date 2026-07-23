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
      character(len=8),dimension(3)::vectorInput
      character(len=8),dimension(3)::intrinsicCharacterVector
      character(len=:),dimension(:),allocatable::vectorOutput
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
      intrinsicCharacterVector =  &
        [ character(len=8)::'  One',' Two','Three' ]
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
      vectorOutput = mqc_adjustl(intrinsicCharacterVector)
      call assertCharacterVector('general vector adjustl',vectorOutput,  &
        [ character(len=8)::'One','Two','Three' ])
      vectorOutput = mqc_adjustr(intrinsicCharacterVector)
      call assertCharacterVector('general vector adjustr',vectorOutput,  &
        [ character(len=8)::'     One','     Two','   Three' ])
      integerVector = mqc_len_trim(intrinsicCharacterVector)
      call assertIntegerVector('general vector len_trim',integerVector,  &
        [5_int64,4_int64,5_int64])
!
      scalarVariable = intrinsicCharacterScalar
      adjustedVariable = ADJUSTL(scalarVariable)
      scalarOutput = adjustedVariable
      call assertCharacter('overloaded scalar adjustl',scalarOutput,'Ab')
      call assertInteger('adjustl preserves scalar length',  &
        LEN(adjustedVariable),8_int64)
      adjustedVariable = scalarVariable%adjustr()
      scalarOutput = adjustedVariable
      call assertCharacter('type-bound scalar adjustr',scalarOutput,  &
        '      Ab')
      lengthVariable = LEN_TRIM(scalarVariable)
      integerScalar = lengthVariable
      call assertInteger('overloaded scalar len_trim',integerScalar,4_int64)
      trimmedVariable = scalarVariable%trim()
      scalarOutput = trimmedVariable
      call assertCharacter('type-bound scalar trim',scalarOutput,'  Ab')
      call assertInteger('trim changes scalar length',  &
        LEN(trimmedVariable),4_int64)
!
      vectorVariable = intrinsicCharacterVector
      adjustedVariable = vectorVariable%adjustl()
      vectorOutput = adjustedVariable
      call assertCharacterVector('type-bound vector adjustl',vectorOutput,  &
        [ character(len=8)::'One','Two','Three' ])
      adjustedVariable = ADJUSTR(vectorVariable)
      vectorOutput = adjustedVariable
      call assertCharacterVector('overloaded vector adjustr',vectorOutput,  &
        [ character(len=8)::'     One','     Two','   Three' ])
      lengthVariable = vectorVariable%len_trim()
      integerVector = [0_int64]
      integerVector = lengthVariable
      call assertIntegerVector('type-bound vector len_trim',integerVector,  &
        [5_int64,4_int64,5_int64])
!
!     Character vector assignment, extraction, insertion, reshape, and
!     allocatable intrinsic assignment.
!
      vectorInput = [ character(len=8)::'Alpha','bETA','Gamma' ]
      vectorVariable = vectorInput
      call assertInteger('vector rank',RANK(vectorVariable),1_int64)
      call assertInteger('vector size',SIZE(vectorVariable),3_int64)
      call assertInteger('vector dimension',SIZE(vectorVariable,1),3_int64)
      call assertInteger('vector length',LEN(vectorVariable),8_int64)
      vectorOutput = vectorVariable
      call assertCharacterVector('vector round trip',vectorOutput,vectorInput)
!
      call vectorVariable%put('Delta',[2_int64])
      elementVariable = vectorVariable%getVal([2_int64])
      scalarOutput = elementVariable
      call assertCharacter('put and getVal',scalarOutput,'Delta')
!
      reshapedVariable = RESHAPE(vectorVariable,[3_int64])
      vectorOutput = reshapedVariable
      call assertCharacterVector('reshape preserves values',vectorOutput,  &
        [ character(len=8)::'Alpha','Delta','Gamma' ])
      changedVariable = MQC_Variable_mqc2mqc(vectorVariable)
      vectorOutput = changedVariable
      call assertCharacterVector('explicit MQC copy',vectorOutput,  &
        [ character(len=8)::'Alpha','Delta','Gamma' ])
!
!     Clear must retain the fixed element width across the new vector.
!
      call clearedVariable%clear('MiXeD',[3_int64])
      call assertInteger('clear length',LEN(clearedVariable),5_int64)
      call clearedVariable%change_case('U')
      vectorOutput = clearedVariable
      call assertCharacterVector('clear and change case',vectorOutput,  &
        [ character(len=5)::'MIXED','MIXED','MIXED' ])
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
      rewind(iOut)
      read(iOut,'(A)') line
      call assertContains('general scalar print',line,'General')
      read(iOut,'(A)') line
      call assertContains('general vector print header',line,'vector')
      read(iOut,'(A)') line
      call assertContains('general vector print element',line,'One')
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
!PROCEDURE assertCharacterVector
      subroutine assertCharacterVector(label,actual,expected)
      implicit none
      character(len=*),intent(in)::label
      character(len=*),dimension(:),intent(in)::actual,expected
!
      if(SIZE(actual).ne.SIZE(expected)) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(label)//' size mismatch'
        error stop 1
      endIf
      if(ANY(actual.ne.expected)) then
        write(*,'(1x,A)') 'FAIL: '//TRIM(label)
        error stop 1
      endIf
!
      return
      end subroutine assertCharacterVector


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


      end program unitTest02
