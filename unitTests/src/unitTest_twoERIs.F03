      program unitTest_twoERIs
!
!     unitTest_twoERIs: Tests mqc_eri_sum and mqc_eri_difference for
!     space/doublebar and space/ms integralType combinations.
!
!     -L. M. Thompson, 2026.
!
!
!     USE Connections
!
      use MQC_EST
      use iso_fortran_env, only: int64, real64
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64),parameter::n=2
      real(kind=real64),parameter::tol=1.0d-12
      integer::nfail
      type(mqc_r4tensor)::tA,tAB,tABBA
      type(mqc_r4tensor)::tB,tBB,tBBBA
      type(mqc_r4tensor)::blk
      type(mqc_scalar)::sc
      type(mqc_twoERIs)::eriA,eriB,eriC
      real(kind=real64)::got,expected
      real(kind=real64)::v2,v3,v4,v5,v6,v7
!
!     Format Statements
!
 1000 format(1x,'PASS  ',a)
 1010 format(1x,'FAIL  ',a,'  got=',f16.10,'  expected=',f16.10)
 9000 format(/,1x,'Tests completed: ',i0,' failure(s)')
!
!
      nfail = 0
      v2=2.0d0 ; v3=3.0d0 ; v4=4.0d0
      v5=5.0d0 ; v6=6.0d0 ; v7=7.0d0
!
!--------------------------------------------------------------------
!     Build input tensors.
!     Only (1,1,1,1) is set to a nonzero value; all others remain 0.
!
!     ERI A:  aaaa=2, aabb=4, abba=6
!     ERI B:  aaaa=3, aabb=5, abba=7
!
!     Sum expected:  aaaa=5,  aabb=9,  abba=13
!     Diff expected: aaaa=-1, aabb=-1, abba=-1
!--------------------------------------------------------------------
!
      call mqc_allocate_r4tensor(n,n,n,n,tA,   'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tAB,  'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tABBA,'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tB,   'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tBB,  'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tBBBA,'real','full')
!
      call tA%put(   v2,1,1,1,1)
      call tAB%put(  v4,1,1,1,1)
      call tABBA%put(v6,1,1,1,1)
      call tB%put(   v3,1,1,1,1)
      call tBB%put(  v5,1,1,1,1)
      call tBBBA%put(v7,1,1,1,1)
!
!--------------------------------------------------------------------
!     Test 1: space/doublebar + space/doublebar
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','space',tA, &
        alphaBeta=tAB,abba=tABBA)
      call mqc_twoeris_allocate(eriB,'full-doublebar','space',tB, &
        alphaBeta=tBB,abba=tBBBA)
      eriC = eriA + eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 5.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+doublebar: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+doublebar: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 9.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+doublebar: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+doublebar: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 13.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+doublebar: abba(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+doublebar: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 2: space/doublebar - space/doublebar
!--------------------------------------------------------------------
!
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -1.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar-doublebar: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar-doublebar: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -1.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar-doublebar: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar-doublebar: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -1.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar-doublebar: abba(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar-doublebar: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 3: space/ms + space/ms
!     ERI A: aaaa=2   ERI B: aaaa=3   Expected sum: aaaa=5
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      call mqc_twoeris_allocate(eriB,'full-ms','space',tB)
      eriC = eriA + eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 5.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+ms: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/ms+ms: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 4: space/ms - space/ms
!     ERI A: aaaa=2   ERI B: aaaa=3   Expected diff: aaaa=-1
!--------------------------------------------------------------------
!
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -1.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms-ms: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/ms-ms: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      end program unitTest_twoERIs
