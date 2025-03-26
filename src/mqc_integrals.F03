      module MQC_Integrals
!MD+
!#MQCPack: MQC_Integrals
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
!
!     This module contains objects and procedures for working with and
!     evaluating certain types of integrals over gaussian-type basis functions.
!MD-

!
      Use MQC_General
      use iso_fortran_env, only: int64, real64
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     MQC_CGTF
!
!     The MQC_CGTF class contains a contracted Gaussian-type function.
      type MQC_CGTF
        logical,private::initialized=.false.
        integer(kind=int64),private::angularMomentum
        integer(kind=int64),private::nPrimMax,nPrim
        real(kind=real64),dimension(3)::center
        real(kind=real64),dimension(:),allocatable::coefficient
        real(kind=real64),dimension(:),allocatable::alpha
        real(kind=real64),dimension(:),allocatable::primitiveNorm
      Contains
        procedure,public::init         => MQC_CGTF_init
        procedure,public::addPrimitive => MQC_CGTF_addPrimitive
        procedure,public::shell2nBasis => MQC_CGTF_nBasisFunctions
        procedure,public::selfOverlap  => MQC_CGTF_selfOverlap
        procedure,public::print        => MQC_CGTF_print
      end type MQC_CGTF
!
!
      CONTAINS
!
!PROCEDURE MQC_CGTG_init
      subroutine MQC_CGTF_init(myCGTF,angularMomentum,nPrimMax,center,  &
        coefficients,alphas)
!
!     This routine initializes a contracted Gaussian-type basis function object.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64),intent(in)::angularMomentum,nPrimMax
      real(kind=real64),dimension(3),intent(in)::center
      real(kind=real64),dimension(:),optional::coefficients,alphas
!
!     Begin by making sure the object's allocatable arrays aren't already
!     allocated.
!
      if(Allocated(myCGTF%coefficient)) DeAllocate(myCGTF%coefficient)
      if(Allocated(myCGTF%alpha)) DeAllocate(myCGTF%alpha)
      if(Allocated(myCGTF%primitiveNorm)) DeAllocate(myCGTF%primitiveNorm)
!
!     Fill scalars and center array. Also, allocate the allocatable arrays in
!     the object. Also, set the initialization flage to TRUE.
!
      myCGTF%angularMomentum = angularMomentum
      myCGTF%nPrimMax        = nPrimMax
      myCGTF%nPrim           = 0
      myCGTF%center          = center
      Allocate(myCGTF%coefficient(nPrimMax),myCGTF%alpha(nPrimMax))
!
!     If the primitive coefficients and alphas have been sent here, we load them
!     into the object.
!
      if(Present(coefficients).neqv.Present(alphas))  &
        call mqc_error('MQC_CGTF_init: Coefficients and exponents must both be sent.')
      if(Present(coefficients)) then
        if(Size(coefficients).ne.Size(alphas))  &
          call mqc_error('MQC_CGTF_init: Coefficients and exponents must be same size.')
        if(Size(coefficients).gt.nPrimMax)  &
          call mqc_error('MQC_CGTF_init: Too many coefficients and exponents sent.')
        myCGTF%nPrim                     = Size(coefficients)
        myCGTF%coefficient(1:myCGTF%nPrim) = coefficients
        myCGTF%alpha(1:myCGTF%nPrim)     = alphas
      endIf
      myCGTF%primitiveNorm = float(1)
!
!     Set the initialization flag to TRUE.
!
      myCGTF%initialized = .true.
!
      return
      end subroutine MQC_CGTF_init

!
!PROCEDURE MQC_CGTF_addPrimitive
      subroutine MQC_CGTF_addPrimitive(myCGTF,coefficient,alpha)
!
!     This routine adds one primitive Gaussian function to the contracted
!     Gaussian-type basis function object myCGTF.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      real(kind=real64),intent(in)::coefficient,alpha
!
!     Check that the myCGTF has been initialized.
!
      if(.not.myCGTF%initialized)  &
        call mqc_error('MQC_CGTF_addPrimitive: CGTF object has not been initiatilized.')
!
!     Add the sent primitive to myCGTF.
!
      if(myCGTF%nPrim.ge.myCGTF%nPrimMax)  &
        call mqc_error('MQC_CGTF_addPrimitive: nPrim > nPrimMax.')
      myCGTF%nPrim = myCGTF%nPrim + 1
      myCGTF%coefficient(myCGTF%nPrim) = coefficient
      myCGTF%alpha(myCGTF%nPrim)       = alpha
!
      return
      end subroutine MQC_CGTF_addPrimitive

!
!PROCEDURE MQC_CGTG_print
      subroutine MQC_CGTF_print(myCGTF,iOut)
!
!     This routine prints the data in a MQC_CGTF object.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64),intent(in)::iOut
      integer(kind=int64)::i
!
 1000 format(1x,'l=',i4,' shell at (',f8.4,',',f8.4,',',f8.4,'). NPrim=',i4)
 2000 format(3x,'Primitive Exponents and Coefficients.')
 2100 format(6x,f15.8,4x,f15.8)
!
      write(iOut,1000) myCGTF%angularMomentum,myCGTF%center,myCGTF%nPrim
      write(iOut,2000)
      do i = 1,myCGTF%nPrim
        write(iOut,2100) myCGTF%alpha(i),myCGTF%coefficient(i)
      endDo
!
      return
      end subroutine MQC_CGTF_print

!
!PROCEDURE MQC_CGTG_nBasisFunctions
      function MQC_CGTF_nBasisFunctions(myCGTF)  result(nBasis)
!
!     This function returns the number of basis functions associated with the
!     shell described by object myCGTF.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64)::nBasis
!
      nBasis = ((myCGTF%angularMomentum+1)*(myCGTF%angularMomentum+2))/2
!
      return
      end function MQC_CGTF_nBasisFunctions

!
!PROCEDURE MQC_CGTG_selfOverlap
      function MQC_CGTF_selfOverlap(myCGTF,i,lVector)  result(selfOverlap)
!
!     This routine evaluates and returns the self-overlap of the i-th primitive
!     basis function in the shell of myCGTF with angular momentum vector
!     lVector.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      integer(kind=int64),intent(in)::i
      integer(kind=int64),dimension(3),intent(in)::lVector
      real(kind=real64)::selfOverlap
      real(kind=real64)::fourAlpha,rootPiFourA3
!
!     Check that lVector is valid.
!
      if(Sum(lVector).ne.myCGTF%angularMomentum) call mqc_error('MQC_CGTF_selfOverlap: Invalid lVector sent.')
!
!     Evaluate the self-overlap.
!
      fourAlpha = 4*myCGTF%alpha(i)
      rootPiFourA3 = (SQRT(Pi/(2*myCGTF%alpha(i))))**3
      if(lVector(1).eq.0) then
        selfOverlap = float(1)
      else
        selfOverlap = mqc_doubleFactorial((2*lVector(1)-1))/(fourAlpha**lVector(1))
      endIf
      if(lVector(2).gt.0) selfOverlap = selfOverlap*(mqc_doubleFactorial((2*lVector(2)-1))/(fourAlpha**lVector(2)))
      if(lVector(3).gt.0) selfOverlap = selfOverlap*(mqc_doubleFactorial((2*lVector(3)-1))/(fourAlpha**lVector(3)))
      selfOverlap = selfOverlap*rootPiFourA3
!
      return
      end function MQC_CGTF_selfOverlap

!
!PROCEDURE MQC_Overlap_Primitive_XYZ_OS
      recursive subroutine MQC_Overlap_Primitive_XYZ_OS(ixyz,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
!
!     This recursive subroutine solves a 1-D Cartesian overlap between two
!     Cartesian Gaussian-type primitive functions using the Obara-Saika
!     recurrence relation.
!
      implicit none
      integer(kind=int64)::ixyz,jxyz
      real(kind=real64)::mu,p,XAB,XPA,XPB
      real(kind=real64),dimension(:,:)::Sij
      logical,dimension(:,:)::haveSij
!
!     Look at ixyz and jxyz to determine what we need to do here. There are 9
!     unique cases.
!
!     0,0
      if(ixyz.eq.0.and.jxyz.eq.0) then
        if(haveSij(1,1)) then
          write(*,*)' Hrant - asked for a value I already have (0,0).'
          return
        endIf
        Sij(1,1) = SQRT(pi/p)*exp(-mu*XAB**2)
        haveSij(1,1) = .true.
!
!     0,1
      elseIf(ixyz.eq.0.and.jxyz.eq.1) then
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(1,2) = XPB*Sij(1,1)
        haveSij(1,2) = .true.
!
!     0,>=2
      elseIf(ixyz.eq.0.and.jxyz.ge.2) then
        if(.not.haveSij(1,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,jxyz-1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-2,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(1,jxyz+1) = XPB*Sij(1,jxyz)+Sij(1,jxyz-1)*(jxyz-1)/(float(2)*p)
        haveSij(1,jxyz+1) = .true.
!
!     1,0
      elseIf(ixyz.eq.1.and.jxyz.eq.0) then
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,1) = XPA*Sij(1,1)
        haveSij(2,1) = .true.
!
!     >=2,0
      elseIf(ixyz.ge.2.and.jxyz.eq.0) then
        if(.not.haveSij(ixyz,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,1) = XPA*Sij(ixyz,1)+Sij(ixyz-1,1)*(ixyz-1)/(float(2)*p)
        haveSij(ixyz+1,1) = .true.
!
!     1,1  ==>  0,1
      elseIf(ixyz.eq.1.and.jxyz.eq.1) then
        if(.not.haveSij(1,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,2) = XPA*Sij(1,2)+Sij(1,1)/(float(2)*p)
        haveSij(2,2) = .true.
!
!     1,>=2
      elseIf(ixyz.eq.1.and.jxyz.ge.2) then
        if(.not.haveSij(1,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,jxyz+1) = XPA*Sij(1,jxyz+1)+jxyz*Sij(1,jxyz)/(float(2)*p)
        haveSij(2,jxyz+1) = .true.
!
!     >=2,1
      elseIf(ixyz.ge.2.and.jxyz.eq.1) then
        if(.not.haveSij(ixyz,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,2) = XPA*Sij(ixyz,2)+(ixyz*Sij(ixyz-1,2)+Sij(ixyz,1))/(float(2)*p)
        haveSij(ixyz+1,2) = .true.
!
!     >=2,>=2
      elseIf(ixyz.ge.2.and.jxyz.ge.2) then
        if(.not.haveSij(ixyz,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,jxyz+1) = XPA*Sij(ixyz,jxyz+1)+((ixyz)*Sij(ixyz-1,jxyz+1)+(jxyz+1)*Sij(ixyz,jxyz))/(float(2)*p)
        haveSij(ixyz+1,jxyz+1) = .true.
!
!     Uknown combination...fail out.
      else
        call mqc_error('MQC_Overlap_Primitive_XYZ_OS: Invalid i and j.')
      endIf
!
      return
      end subroutine MQC_Overlap_Primitive_XYZ_OS

!
!
      end module MQC_Integrals
