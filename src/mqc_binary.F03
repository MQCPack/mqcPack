      module MQC_Binary
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
!     Interface MQC_Bits to the custom constructor.
      interface mqc_bits
        module procedure MQC_Bits_Initialize
      end interface mqc_bits
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
      type(mqc_bits),intent(out)::mqc_bits_out
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
      type(mqc_bits),intent(in)::mqc_bits_in
      integer(kind=int64),intent(in),optional::iUnit
      character(len=*),intent(in),optional::Header
      logical,intent(in),optional::Blank_At_Top,Blank_At_Bottom,verbose
!
      integer(kind=int64)::my_iUnit,i
      logical::my_verbose
      character(len=256)::formatText,charTemp
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
        if(PRESENT(Header)) then
          write(UNIT=my_iUnit,FMT=1100,ADVANCE="no") TRIM(Header)
        else
          write(UNIT=my_iUnit,FMT=1000,ADVANCE="no")
        endIf
        formatText = '(B'//TRIM(charTemp)//'.'//TRIM(charTemp)//')'
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
          write(UNIT=my_iUnit,FMT=formatText,ADVANCE="no") mqc_bits_in%integerList(i)
        endDo
        write(UNIT=my_iUnit,FMT=formatText,ADVANCE="yes") mqc_bits_in%integerList(1)
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
      logical,intent(out)::bitValue
!
      integer(kind=int64)::iInteger,iBitNum
!
!
!     Do the work.
!
      call MQC_BitPosition(mqc_bits_in,bitNum,iInteger,iBitNum)
      bitValue = BTest(mqc_bits_in%integerList(iInteger),iBitnum)
!
      return
      end function MQC_BTest
!
!
!PROCEDURE MQC_IBitSet
      subroutine MQC_IBitSet(mqc_bits_in,bitnum)
!
!     This subroutine is used to set a bit in mqc_bits_in.
!
!     H. P. Hratchian, 2021.
!
!
      implicit none
      type(mqc_bits),intent(inOut)::mqc_bits_in
      integer(kind=int64),intent(in)::bitnum
!
      integer(kind=int64)::iInteger,iBitnum
!
!
!     Do the work.
!
      call MQC_BitPosition(mqc_bits_in,bitnum,iInteger,iBitnum)
      mqc_bits_in%integerList(iInteger) = IBSet(mqc_bits_in%integerList(iInteger),iBitnum)
!
      return
      end subroutine MQC_IBitSet
!
!
!PROCEDURE MQC_IBitClr
      subroutine MQC_IBitClr(mqc_bits_in,bitnum)
!
!     This subroutine is used to clear a bit in mqc_bits_in.
!
!     H. P. Hratchian, 2022.
!
!
      implicit none
      type(mqc_bits),intent(inOut)::mqc_bits_in
      integer(kind=int64),intent(in)::bitnum
!
      integer(kind=int64)::iInteger,iBitnum
!
!
!     Do the work.
!
      call MQC_BitPosition(mqc_bits_in,bitnum,iInteger,iBitnum)
      mqc_bits_in%integerList(iInteger) = IBClr(mqc_bits_in%integerList(iInteger),iBitnum)
!
      return
      end subroutine MQC_IBitClr
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
      type(mqc_bits),intent(out)::mqc_bits_out
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
      type(mqc_bits),intent(out)::mqc_bits_out
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
      type(mqc_bits),intent(out)::mqc_bits_out
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
      integer(kind=int64),intent(out)::bits_on
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


      end module MQC_Binary
