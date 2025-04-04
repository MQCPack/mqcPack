      module MQC_Binary
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
!     **                          Version 25.3.0                           **
!     **                          March 25, 2025                           **
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
!     NEED TO ADD COMMENTS
!
!
      USE MQC_General
      USE iso_fortran_env, only: int64, real64
!                                                                    
!----------------------------------------------------------------    
!                                                               |    
!     TYPE AND CLASS DEFINITIONS                                |    
!                                                               |    
!----------------------------------------------------------------    
      implicit none
!
!     mqc_bits
      type :: mqc_bits
        integer(kind=int64)::nBitsPerInteger=3
        integer(kind=int64)::nBits,nIntegers
        integer(kind=int64),dimension(:),allocatable::integerList
      contains
        procedure, public::print => mqc_bits_print
      end type mqc_bits
!
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!
!>    \brief <b> Print MQC_Bits</b>
      interface mqc_print
        module procedure MQC_Bits_Print
      end interface mqc_print
!
!>    \brief <b> Interface MQC_Bits to the custom constructor</b>
      interface mqc_bits
        module procedure MQC_Bits_Initialize
      end interface mqc_bits
!
!>    \brief <b> Tests if bit set</b>
      interface btest
        module procedure MQC_BTest
      end interface btest
!
!>    \brief <b> Set a bit</b>
      interface ibset
        module procedure MQC_IBitSet
      end interface ibset
!
!>    \brief <b> Clear a bit</b>
      interface ibclr
        module procedure MQC_IBitClr
      end interface ibclr
!
!>    \brief <b> Performs and test</b>
      interface iand
        module procedure MQC_IAnd
      end interface iand
!
!>    \brief <b> Performs exclusive or test</b>
      interface ieor
        module procedure MQC_IEOr
      end interface ieor
!
!>    \brief <b> Performs inclusive or test</b>
      interface ior
        module procedure MQC_IOr
      end interface ior
!
!>    \brief <b> Returns the number of set bits</b>
      interface popcnt
        module procedure MQC_PopCnt
      end interface popcnt
!
!>    \brief <b> Moves a set number of bits</b>
      interface mvBits
        module procedure MQC_mvBits
      end interface mvBits
!
!>    \brief <b> Extracts a set number of bits</b>
      interface iBits
        module procedure MQC_iBits
      end interface iBits
!
!>    \brief <b> Shifts by a set number of bits</b>
      interface iShft
        module procedure MQC_iShft
      end interface iShft
!
!
!----------------------------------------------------------------
!                                                               |
!     SUBROUTINES AND FUNCTIONS                                 |
!                                                               |
!----------------------------------------------------------------
!
!
      CONTAINS
!
!PROCEDURE MQC_Bits_Initialize
      function MQC_Bits_Initialize(nBits) result(mqc_bits_out)
!
!     This function is used to initialize an MQC_Bit object. Input dummy
!     argument <nBits> is the number of bits requested by the user. Keep in mind
!     that bit 0 counts as a bit in this count. As such, requesting 8 bits
!     provides bit numbers 0-7.
!
!     H. P. Hratchian, 2021.
!
!
      implicit none
      type(mqc_bits)::mqc_bits_out
      integer(kind=int64),intent(in),optional::nBits
!
      integer(kind=int64)::my_nBits,nIntegers
!
!
!     Start by setting up my_nBits.
!
      if(PRESENT(nBits)) then
        my_nBits = nBits
      else
        my_nBits = bit_size(nIntegers)
      endIf
!
!     Set nIntegers.
!
      nIntegers = my_nBits/mqc_bits_out%nBitsPerInteger
      if(MOD(my_nBits,mqc_bits_out%nBitsPerInteger).ne.0)  &
        nIntegers = nIntegers+1
!
!     Set mqc_bits_out%nBits. Then allocate the integer array in the bit object
!     and set all the values to 0.
!
      mqc_bits_out%nBits = my_nBits
      mqc_bits_out%nIntegers = nIntegers
      Allocate(mqc_bits_out%integerList(nIntegers))
      mqc_bits_out%integerList = 0
!
      return
      end function MQC_Bits_Initialize
!
!
!PROCEDURE MQC_Bits_Print
      subroutine MQC_Bits_Print(mqc_bits_in,iUnit,Header,Blank_At_Top,  &
        Blank_At_Bottom,verbose)
!
!     This subroutine is used to print out the bits in mqc_bits_in.
!
!     H. P. Hratchian, 2021.
!
!
      implicit none
      class(mqc_bits),intent(in)::mqc_bits_in
      integer(kind=int64),intent(in),optional::iUnit
      character(len=*),intent(in),optional::Header
      logical,intent(in),optional::Blank_At_Top,Blank_At_Bottom,verbose
!
      integer(kind=int64)::my_iUnit,i
      logical::my_verbose
      character(len=256)::formatText,charTemp,formatTextLast,charTempLast
!
!
!     Format Statements
!
 1000 format(1x,'Binary number: ')
 1100 format(1x,A,' ')
!
!
!     Start by setting up my_iUnit and my_verbose.
!
      if(PRESENT(iUnit)) then
        my_iUnit = iUnit
      else
        my_iUnit = 6
      endIf
      my_verbose = .false.
      if(PRESENT(verbose)) my_verbose = verbose
!
!     Print blank at top if requested.
!
      if(PRESENT(Blank_At_Top)) then
        if(Blank_At_Top) write(my_iUnit)
      endIf
!
!     Set up the format string.
!
      charTemp = num2char(mqc_bits_in%nBitsPerInteger)
      if(my_verbose) then
        if(PRESENT(Header)) then
          write(my_iUnit,1100) TRIM(Header)
        else
          write(my_iUnit,1000)
        endIf
        formatText = '(1x,I2,'':'',B'//TRIM(charTemp)//'.'//TRIM(charTemp)//')'
      else
        charTempLast = num2char(mod(mqc_bits_in%nBits,mqc_bits_in%nBitsPerInteger))
        if(PRESENT(Header)) then
          write(UNIT=my_iUnit,FMT=1100,ADVANCE="no") TRIM(Header)
        else
          write(UNIT=my_iUnit,FMT=1000,ADVANCE="no")
        endIf
        formatText = '(B'//TRIM(charTemp)//'.'//TRIM(charTemp)//')'
        formatTextLast = '(B'//TRIM(charTempLast)//'.'//TRIM(charTempLast)//')'
      endIf
!
!     Write out the binary number.
!
      if(my_verbose) then
        do i = 1,mqc_bits_in%nIntegers
          write(UNIT=my_iUnit,FMT=formatText,ADVANCE="yes") i,mqc_bits_in%integerList(i)
        endDo
      else
        do i = mqc_bits_in%nIntegers,2,-1
          if(i.eq.mqc_bits_in%nIntegers) then
            write(UNIT=my_iUnit,FMT=formatTextLast,ADVANCE="no") mqc_bits_in%integerList(i)
          else
            write(UNIT=my_iUnit,FMT=formatText,ADVANCE="no") mqc_bits_in%integerList(i)
          endIf
        endDo
        if(mqc_bits_in%nIntegers.eq.1) then
          write(UNIT=my_iUnit,FMT=formatTextLast,ADVANCE="yes") mqc_bits_in%integerList(1)
        else
          write(UNIT=my_iUnit,FMT=formatText,ADVANCE="yes") mqc_bits_in%integerList(1)
        endIf
      endIf
!
!     Print blank at bottom if requested.
!
      if(PRESENT(Blank_At_Bottom)) then
        if(Blank_At_Bottom) write(my_iUnit)
      endIf
!
      return
      end subroutine MQC_Bits_Print
!
!
!PROCEDURE MQC_BitPosition
      subroutine MQC_BitPosition(mqc_bits_in,bitnum,iInteger,iBitnum)
!
!     This subroutine is used to determine which integer <iInteger> and which
!     bit <iBitnum> within that integer in an MQC_Bits object correspond the
!     overall bit <bitnum>.
!
!     H. P. Hratchian, 2021.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_in
      integer(kind=int64),intent(in)::bitnum
      integer(kind=int64),intent(out)::iInteger,iBitnum
!
!
!     Do the work.
!
      iInteger = (bitnum+1)/mqc_bits_in%nBitsPerInteger
      if(MOD((bitnum+1),mqc_bits_in%nBitsPerInteger).ne.0)  &
        iInteger = iInteger+1
      iBitnum = MOD(bitnum,mqc_bits_in%nBitsPerInteger)
!
      return
      end subroutine MQC_BitPosition
!
!
!PROCEDURE MQC_BTest
      function MQC_BTest(mqc_bits_in,bitnum) result(bitValue)
!
!     This function tests the switch in bit number <bitNum> in <mqc_bits_in> and
!     returns TRUE or FALSE accordingly.
!
!     H. P. Hratchian, 2022.
!
!
      implicit none
      type(mqc_bits),intent(In)::mqc_bits_in
      integer(kind=int64),intent(in)::bitNum
      logical::bitValue
!
      integer(kind=int64)::iInteger,iBitNum
!
!
!     Do the work.
!
      if(bitnum.ge.mqc_bits_in%nBits) call mqc_error_i('Bit number requested &
        & greater than number of bits in MQC_BTest',6,'bitnum',bitnum,&
        'mqc_bits_in%nBits',mqc_bits_in%nBits)
      call MQC_BitPosition(mqc_bits_in,bitNum,iInteger,iBitNum)
      bitValue = BTest(mqc_bits_in%integerList(iInteger),iBitnum)
!
      return
      end function MQC_BTest
!
!
!PROCEDURE MQC_IBitSet
      function MQC_IBitSet(mqc_bits_in,bitnum) result(mqc_bits_out)
!
!     This subroutine is used to set a bit in mqc_bits_in.
!
!     H. P. Hratchian, 2021.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_in
      integer(kind=int64),intent(in)::bitnum
      type(mqc_bits)::mqc_bits_out
!
      integer(kind=int64)::iInteger,iBitnum
!
!
!     Do the work.
!
      mqc_bits_out = mqc_bits_in
      call MQC_BitPosition(mqc_bits_out,bitnum,iInteger,iBitnum)
      mqc_bits_out%integerList(iInteger) = IBSet(mqc_bits_out%integerList(iInteger),iBitnum)
!
      return
      end function MQC_IBitSet
!
!
!PROCEDURE MQC_IBitClr
      function MQC_IBitClr(mqc_bits_in,bitnum) result(mqc_bits_out)
!
!     This function is used to clear a bit in mqc_bits_in.
!
!     H. P. Hratchian, 2022.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_in
      integer(kind=int64),intent(in)::bitnum
      type(mqc_bits)::mqc_bits_out
!
      integer(kind=int64)::iInteger,iBitnum
!
!
!     Do the work.
!
      mqc_bits_out = mqc_bits_in
      call MQC_BitPosition(mqc_bits_out,bitnum,iInteger,iBitnum)
      mqc_bits_out%integerList(iInteger) = IBClr(mqc_bits_out%integerList(iInteger),iBitnum)
!
      return
      end function MQC_IBitClr
!
!
!PROCEDURE MQC_IAnd
      function MQC_IAnd(mqc_bits_in1,mqc_bits_in2) result(mqc_bits_out)
!
!     This function carries out the IAnd operation on two MQC_Bit objects. The
!     output is the resulting MQC_Bit object.
!
!     H. P. Hratchian, 2022.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_in1,mqc_bits_in2
      type(mqc_bits)::mqc_bits_out
!
      integer(kind=int64)::i
!
!
!     Do the work.
!
      if((mqc_bits_in1%nIntegers.eq.mqc_bits_in2%nIntegers).and.  &
        (mqc_bits_in1%nBits.eq.mqc_bits_in2%nBits)) then
        mqc_bits_out = MQC_Bits_Initialize(mqc_bits_in1%nBits)
        do i = 1,mqc_bits_in1%nIntegers
          mqc_bits_out%integerList(i) =   &
            IAND(mqc_bits_in1%integerList(i),mqc_bits_in2%integerList(i))
        endDo
      else
        call mqc_error('MQC_IAnd: Non-conformable bit objects.')
      endIf
!
      return
      end function MQC_IAnd
!
!
!PROCEDURE MQC_IEOR
      function MQC_IEOR(mqc_bits_in1,mqc_bits_in2) result(mqc_bits_out)
!
!     This function carries out the EXCLUSIVE OR operation on two MQC_Bit
!     objects. The output is the resulting MQC_Bit object.
!
!     H. P. Hratchian, 2022.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_in1,mqc_bits_in2
      type(mqc_bits)::mqc_bits_out
!
      integer(kind=int64)::i
!
!
!     Do the work.
!
      if((mqc_bits_in1%nIntegers.eq.mqc_bits_in2%nIntegers).and.  &
        (mqc_bits_in1%nBits.eq.mqc_bits_in2%nBits)) then
        mqc_bits_out = MQC_Bits_Initialize(mqc_bits_in1%nBits)
        do i = 1,mqc_bits_in1%nIntegers
          mqc_bits_out%integerList(i) =   &
            IEOR(mqc_bits_in1%integerList(i),mqc_bits_in2%integerList(i))
        endDo
      else
        call mqc_error('MQC_IEOR: Non-conformable bit objects.')
      endIf
!
      return
      end function MQC_IEOR
!
!
!PROCEDURE MQC_IOR
      function MQC_IOR(mqc_bits_in1,mqc_bits_in2) result(mqc_bits_out)
!
!     This function carries out the INCLUSIVE OR operation on two MQC_Bit
!     objects. The output is the resulting MQC_Bit object.
!
!     H. P. Hratchian, 2022.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_in1,mqc_bits_in2
      type(mqc_bits)::mqc_bits_out
!
      integer(kind=int64)::i
!
!
!     Do the work.
!
      if((mqc_bits_in1%nIntegers.eq.mqc_bits_in2%nIntegers).and.  &
        (mqc_bits_in1%nBits.eq.mqc_bits_in2%nBits)) then
        mqc_bits_out = MQC_Bits_Initialize(mqc_bits_in1%nBits)
        do i = 1,mqc_bits_in1%nIntegers
          mqc_bits_out%integerList(i) =   &
            IOR(mqc_bits_in1%integerList(i),mqc_bits_in2%integerList(i))
        endDo
      else
        call mqc_error('MQC_IOR: Non-conformable bit objects.')
      endIf
!
      return
      end function MQC_IOR
!
!
!PROCEDURE MQC_PopCnt
      function MQC_PopCnt(mqc_bits_in) result(bits_on)
!
!     This function returns the number of "ON" bits in input dummy argument
!     <mqc_bits_in>. The output of the function is an int64 integer.
!
!     H. P. Hratchian, 2022.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_in
      integer(kind=int64)::bits_on
!
      integer(kind=int64)::i
!
!
!     Do the work.
!
      bits_on = 0
      do i = 0,mqc_bits_in%nBits-1
        if(mqc_BTest(mqc_bits_in,i)) bits_on = bits_on+1
      endDo
!
      return
      end function MQC_PopCnt
!
!
!PROCEDURE MQC_MvBits
      subroutine MQC_MvBits(mqc_bits_src,fromPos,length,mqc_bits_des,toPos)
!
!     This subroutine moves length bits from position fromPos through
!     fromPos+length-1 bits of mqc_bits_src to positions toPos through
!     toPos+length-1 of mqc_bits_des.
!
!     L. M. Thompson, 2024.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_src
      type(mqc_bits),intent(inOut)::mqc_bits_des
      integer(kind=int64),intent(in)::fromPos,length,toPos
!
      integer(kind=int64)::i
      type(mqc_bits)::mqc_bits_temp
!
!
!     Do the work.
!
      mqc_bits_temp = mqc_bits(max(toPos+length,mqc_bits_des%nbits))
      do i = 0, toPos-1
        if(bTest(mqc_bits_des,i)) mqc_bits_temp = ibset(mqc_bits_temp,i)
      endDo
      do i = 0, length-1
        if(bTest(mqc_bits_src,fromPos+i)) mqc_bits_temp = ibset(mqc_bits_temp,toPos+i)
      endDo
      if(mqc_bits_des%nbits.gt.toPos+length) then
        do i = 0, mqc_bits_des%nbits-toPos+length
          if(bTest(mqc_bits_des,toPos+length+i)) &
            mqc_bits_temp = ibset(mqc_bits_temp,toPos+length+i)
        endDo
      endIf
      mqc_bits_des = mqc_bits_temp
!
      return
      end subroutine MQC_MvBits
!
!
!PROCEDURE MQC_iBits
      function MQC_iBits(mqc_bits_in,pos,length) result(mqc_bits_out)
!
!     This function extracts a field of length from mqc_bits_in, starting
!     from Pos and extending left for length bits. 
!
!     L. M. Thompson, 2024.
!
!
      implicit none
      type(mqc_bits),intent(in)::mqc_bits_in
      integer(kind=int64),intent(in)::pos,length
      type(mqc_bits)::mqc_bits_out
!
      integer(kind=int64)::i
      type(mqc_bits)::mqc_bits_temp
!
!
!     Do the work.
!
      mqc_bits_out = mqc_bits(length)
      do i = 0, length-1
        if(bTest(mqc_bits_in,pos+i)) mqc_bits_out = ibset(mqc_bits_out,i)
      endDo
!
      return
      end function MQC_iBits
!
!
!PROCEDURE MQC_iShft
      function MQC_iShft(bits,shift) result(bitsOut)
!
!     This function concatenates two mqc_bits.
!
!     L. M. Thompson, 2025.
!
!
      implicit none
      type(mqc_bits),intent(in)::bits
      integer(kind=int64),intent(in)::shift
      type(mqc_bits)::bitsOut
!
      integer(kind=int64)::i
!
!
!     Do the work.
!
      bitsOut = mqc_bits(bits%nBits+shift)
      do i = 0, bits%nBits-1
        if(i+shift.ge.0.or.i+shift.le.bitsOut%nBits) then
          if(bTest(bits,i)) bitsOut = ibset(bitsOut,i+shift)
        endIf
      endDo
!
      return
      end function MQC_iShft 
!
!
!PROCEDURE MQC_concatenate_bits
      function MQC_concatenate_bits(bits1,bits2) result(bitsOut)
!
!     This function concatenates two mqc_bits.
!
!     L. M. Thompson, 2025.
!
!
      implicit none
      type(mqc_bits),intent(in)::bits1,bits2
      type(mqc_bits)::bitsOut
!
      integer(kind=int64)::i
!
!
!     Do the work.
!
      bitsOut = ishft(bits1,bits2%nbits)
      do i = 0, bits2%nBits-1
        if(bTest(bits2,i)) bitsOut = ibset(bitsOut,i)
      endDo
!
      return
      end function MQC_concatenate_bits
!
      end module MQC_Binary
