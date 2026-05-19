      program unitTest_twoERIs
!
!     unitTest_twoERIs: Tests mqc_eri_sum and mqc_eri_difference for
!     space/space and space/spin integralType combinations.
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
      type(mqc_r4tensor)::tC,tD,tE,tF,tG,tZ
      type(mqc_r4tensor)::blk
      type(mqc_scalar)::sc
      type(mqc_twoERIs)::eriA,eriB,eriC
      real(kind=real64)::got,expected
      real(kind=real64)::v2,v3,v4,v5,v6,v7
      real(kind=real64)::v8,v9,v10,v11,v12
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
      v8=8.0d0 ; v9=9.0d0 ; v10=10.0d0 ; v11=11.0d0 ; v12=12.0d0
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
!     Spin ERI blocks: tC=8(aaaa), tD=9(bbbb), tE=10(aabb),
!                      tF=11(abba), tG=12(abab), tZ=0(zero)
!
      call mqc_allocate_r4tensor(n,n,n,n,tC,'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tD,'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tE,'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tF,'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tG,'real','full')
      call mqc_allocate_r4tensor(n,n,n,n,tZ,'real','full')
      call tC%put(v8, 1,1,1,1)
      call tD%put(v9, 1,1,1,1)
      call tE%put(v10,1,1,1,1)
      call tF%put(v11,1,1,1,1)
      call tG%put(v12,1,1,1,1)
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
!--------------------------------------------------------------------
!     Tests 5-8: space/coulomb and space/exchange cross-type operations.
!     Physical identities: C-E=DB, C-DB=E, DB-C=E, DB+E=C.
!
!     coulomb ERI: aaaa=2 (tA)
!     exchange ERI: aaaa=3 (tB)
!     doublebar ERI: aaaa=4 (tAB), aabb=5 (tBB), abba=6 (tABBA)
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tA)
      call mqc_twoeris_allocate(eriB,'full-exchange','space',tB)
      call mqc_twoeris_allocate(eriC,'full-doublebar','space',tAB, &
        alphaBeta=tBB,abba=tABBA)
!
!--------------------------------------------------------------------
!     Test 5: space/coulomb - space/exchange = space/doublebar
!     C-E: aaaa=2-3=-1, aabb=C.aaaa=2, abba=-E.aaaa=-3
!--------------------------------------------------------------------
!
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -1.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb-exchange: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb-exchange: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb-exchange: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb-exchange: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -3.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb-exchange: abba(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb-exchange: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 6: space/coulomb - space/doublebar = space/exchange
!     C-DB: aaaa=2-4=-2
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tA)
      call mqc_twoeris_allocate(eriB,'full-doublebar','space',tAB, &
        alphaBeta=tBB,abba=tABBA)
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb-doublebar: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb-doublebar: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 7: space/doublebar - space/coulomb = space/exchange
!     DB-C: aaaa=4-2=2
!--------------------------------------------------------------------
!
      eriC = eriB - eriA
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar-coulomb: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar-coulomb: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 8: space/doublebar + space/exchange = space/coulomb
!     DB+E: aaaa=4+3=7
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-exchange','space',tB)
      eriC = eriB + eriA
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+exchange: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+exchange: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 9-18: space/spin operations.
!     ERIA is always space type with aaaa=tA(2) unless noted.
!     ERIB is spin type: aaaa=tC(8), bbbb=tD(9), aabb=tE(10),
!                        abba=tF(11), abab=tG(12), zero=tZ(0).
!
!     For space/coulomb ERIA (aaaa=2):
!       getBlock: aaaa=2, bbbb=2, aabb=2, abba=0, abab=0
!     For space/exchange ERIA (aaaa=2):
!       getBlock: aaaa=2, bbbb=2, aabb=0, abba=2, abab=0
!     For space/inverted ERIA (aaaa=2):
!       getBlock: aaaa=2, bbbb=2, abba=0, abab=2, baba=2
!     For space/doublebar ERIA (aaaa=2,aabb=4,abba=6):
!       getBlock: aaaa=2, bbbb=2, aabb=4, abba=6
!     For space/ms ERIA (aaaa=2): all blocks=2
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!     Test 9: space/coulomb(aaaa=2) + spin/coulomb_symm(aaaa=8,bbbb=9,aabb=10)
!     Expected output spin/coulomb_symm: aaaa=10, bbbb=11, aabb=12
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tC,tD,tE)
      eriC = eriA + eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/coulomb: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/coulomb: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('bbbb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/coulomb: bbbb(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/coulomb: bbbb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/coulomb: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/coulomb: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 10: space/coulomb(aaaa=2) - spin/coulomb_symm(aaaa=8,bbbb=9,aabb=10)
!     Expected: aaaa=-6, bbbb=-7, aabb=-8
!--------------------------------------------------------------------
!
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb-spin/coulomb: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb-spin/coulomb: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('bbbb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb-spin/coulomb: bbbb(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb-spin/coulomb: bbbb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -8.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb-spin/coulomb: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb-spin/coulomb: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 11: space/exchange(aaaa=2) + spin/exchange_symm(aaaa=8,bbbb=9,abba=11)
!     Expected output spin/exchange_symm: aaaa=10, bbbb=11, abba=13
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-exchange','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-exchange','spin',tC,tD,tF)
      eriC = eriA + eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/exchange+spin/exchange: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/exchange+spin/exchange: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('bbbb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/exchange+spin/exchange: bbbb(1,1,1,1)'
      else
        write(*,1010) 'space/exchange+spin/exchange: bbbb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 13.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/exchange+spin/exchange: abba(1,1,1,1)'
      else
        write(*,1010) 'space/exchange+spin/exchange: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 12: space/exchange(aaaa=2) - spin/exchange_symm(aaaa=8,bbbb=9,abba=11)
!     Expected: aaaa=-6, abba=-9
!--------------------------------------------------------------------
!
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/exchange-spin/exchange: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/exchange-spin/exchange: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -9.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/exchange-spin/exchange: abba(1,1,1,1)'
      else
        write(*,1010) 'space/exchange-spin/exchange: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 13: space/inverted(aaaa=2) + spin/inverted_full(aaaa=8,bbbb=9,abab=12,baba=0)
!     Expected output spin/inverted_full: aaaa=10, abab=14
!     (verifies inverted+inverted variable name fix)
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-inverted','space',tA)
      call mqc_twoeris_allocate(eriB,'full-inverted','spin',tC,tD,tG,tZ)
      eriC = eriA + eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted+spin/inverted: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/inverted+spin/inverted: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abab')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted+spin/inverted: abab(1,1,1,1)'
      else
        write(*,1010) 'space/inverted+spin/inverted: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 14: space/inverted(aaaa=2) - spin/inverted_full(aaaa=8,bbbb=9,abab=12,baba=0)
!     Expected: aaaa=-6, abab=-10
!--------------------------------------------------------------------
!
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted-spin/inverted: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/inverted-spin/inverted: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abab')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted-spin/inverted: abab(1,1,1,1)'
      else
        write(*,1010) 'space/inverted-spin/inverted: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 15: space/doublebar(aaaa=2,aabb=4,abba=6)
!            + spin/doublebar_symm(aaaa=8,bbbb=9,aabb=10,abba=11)
!     Expected output spin/doublebar_symm: aaaa=10, bbbb=11, aabb=14, abba=17
!     (verifies allocate bug fix: abba now permitted for spin/doublebar)
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','space',tA, &
        alphaBeta=tAB,abba=tABBA)
      call mqc_twoeris_allocate(eriB,'symm-doublebar','spin',tC,tD,tE,abba=tF)
      eriC = eriA + eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+spin/doublebar: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+spin/doublebar: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('bbbb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+spin/doublebar: bbbb(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+spin/doublebar: bbbb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+spin/doublebar: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+spin/doublebar: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+spin/doublebar: abba(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+spin/doublebar: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 16: space/doublebar(aaaa=2,aabb=4,abba=6)
!            - spin/doublebar_symm(aaaa=8,bbbb=9,aabb=10,abba=11)
!     Expected: aaaa=-6, aabb=-6, abba=-5
!--------------------------------------------------------------------
!
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar-spin/doublebar: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar-spin/doublebar: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar-spin/doublebar: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar-spin/doublebar: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -5.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar-spin/doublebar: abba(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar-spin/doublebar: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 17: space/ms(aaaa=2)
!            + spin/ms_symm(aaaa=8,bbbb=9,aabb=10,abba=11,abab=12)
!     Expected output spin/ms_symm: aaaa=10, bbbb=11, aabb=12, abba=13
!     (verifies output type fix: ms+ms → ms, not doublebar)
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-ms','spin',tC,tD,tE,abba=tF,abab=tG)
      eriC = eriA + eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/ms: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/ms: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('bbbb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/ms: bbbb(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/ms: bbbb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/ms: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/ms: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 13.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/ms: abba(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/ms: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 18: space/ms(aaaa=2)
!            - spin/ms_symm(aaaa=8,bbbb=9,aabb=10,abba=11,abab=12)
!     Expected: aaaa=-6, aabb=-8, abba=-9
!--------------------------------------------------------------------
!
      eriC = eriA - eriB
!
      blk = eriC%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms-spin/ms: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/ms-spin/ms: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -8.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms-spin/ms: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/ms-spin/ms: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -9.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms-spin/ms: abba(1,1,1,1)'
      else
        write(*,1010) 'space/ms-spin/ms: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 19: Direct getBlock test for spin/ms_symm
!     Verifies that each stored block is retrieved correctly.
!     spin/ms_symm: aaaa=8, bbbb=9, aabb=10, abba=11, abab=12
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'symm-ms','spin',tC,tD,tE,abba=tF,abab=tG)
!
      blk = eriA%getBlock('aaaa')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 8.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms getBlock: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/ms getBlock: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriA%getBlock('bbbb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 9.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms getBlock: bbbb(1,1,1,1)'
      else
        write(*,1010) 'spin/ms getBlock: bbbb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriA%getBlock('aabb')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms getBlock: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/ms getBlock: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriA%getBlock('abba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms getBlock: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/ms getBlock: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      end program unitTest_twoERIs
