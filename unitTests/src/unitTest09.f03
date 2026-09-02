      Program unitTest09
!
!     Validate element symbols and the Bragg-Slater and original-Becke
!     atomic-size radius APIs in MQC_Molecule.
!
!     H. P. Hratchian, 2026.
!
      Use iso_fortran_env, Only: int64,real64
      Use MQC_General, Only: angPBohr
      Use MQC_Molecule, Only: MQC_BRAGG_SLATER_MAX_ATOMIC_NUMBER,  &
        mqc_atomic_number,mqc_element_becke_1988_radius,  &
        mqc_element_bragg_slater_radius,  &
        mqc_element_has_bragg_slater_radius,mqc_element_symbol
      Implicit None
!
      Call assertTrue(mqc_element_has_bragg_slater_radius(1_int64),  &
        'Hydrogen Bragg-Slater radius is not reported as available.')
      Call assertTrue(  &
        mqc_element_has_bragg_slater_radius(  &
          MQC_BRAGG_SLATER_MAX_ATOMIC_NUMBER),  &
        'The last tabulated Bragg-Slater radius is not available.')
      Call assertTrue(.not.mqc_element_has_bragg_slater_radius(0_int64),  &
        'A ghost center was assigned a Bragg-Slater radius.')
      Call assertTrue(.not.mqc_element_has_bragg_slater_radius(87_int64),  &
        'An unsupported Bragg-Slater radius was reported as available.')
!
      Call assertNear(  &
        mqc_element_bragg_slater_radius(1_int64)*angPBohr,  &
        0.25_real64,1.0e-12_real64,'Slater hydrogen radius')
      Call assertNear(  &
        mqc_element_becke_1988_radius(1_int64)*angPBohr,  &
        0.35_real64,1.0e-12_real64,'Becke hydrogen radius')
      Call assertNear(  &
        mqc_element_bragg_slater_radius(6_int64)*angPBohr,  &
        0.70_real64,1.0e-12_real64,'Slater carbon radius')
      Call assertNear(  &
        mqc_element_becke_1988_radius(8_int64)*angPBohr,  &
        0.60_real64,1.0e-12_real64,'Becke oxygen radius')
      Call assertNear(  &
        mqc_element_bragg_slater_radius(86_int64)*angPBohr,  &
        2.10_real64,1.0e-12_real64,'Slater radon radius')
!
!     Preserve coverage of the existing complete symbol table while this
!     test is the molecule-module unit test.
!
      Call assertTrue(mqc_element_symbol(118_int64).eq.'Og',  &
        'Atomic number 118 did not map to Og.')
      Call assertTrue(mqc_atomic_number('Og').eq.118_int64,  &
        'Element symbol Og did not map to atomic number 118.')
!
      Write(*,'(1x,A)') 'unitTest09: PASS'
!
      Contains


!PROCEDURE assertNear
      Subroutine assertNear(value,target,tolerance,label)
      Implicit None
      Real(kind=real64),Intent(In)::value,target,tolerance
      Character(len=*),Intent(In)::label
!
      If(abs(value-target).gt.tolerance) then
        Write(*,'(1x,A)') 'FAIL: '//Trim(label)
        Write(*,'(1x,A,ES24.16)') 'value  = ',value
        Write(*,'(1x,A,ES24.16)') 'target = ',target
        Error Stop 1
      endIf
!
      Return
      End Subroutine assertNear


!PROCEDURE assertTrue
      Subroutine assertTrue(condition,message)
      Implicit None
      Logical,Intent(In)::condition
      Character(len=*),Intent(In)::message
!
      If(.not.condition) then
        Write(*,'(1x,A)') 'FAIL: '//Trim(message)
        Error Stop 1
      endIf
!
      Return
      End Subroutine assertTrue


      End Program unitTest09
