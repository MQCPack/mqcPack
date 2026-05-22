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
!     Expected output spin/ms_symm: aaaa=10, bbbb=11, aabb=12, abba=13, abab=14
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
      blk = eriC%getBlock('abab')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/ms: abab(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/ms: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('baba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/ms: baba(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/ms: baba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 18: space/ms(aaaa=2)
!            - spin/ms_symm(aaaa=8,bbbb=9,aabb=10,abba=11,abab=12)
!     Expected: aaaa=-6, aabb=-8, abba=-9, abab=-10
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
      blk = eriC%getBlock('abab')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms-spin/ms: abab(1,1,1,1)'
      else
        write(*,1010) 'space/ms-spin/ms: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriC%getBlock('baba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = -10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms-spin/ms: baba(1,1,1,1)'
      else
        write(*,1010) 'space/ms-spin/ms: baba(1,1,1,1)',got,expected
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
      blk = eriA%getBlock('abab')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms getBlock: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/ms getBlock: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      blk = eriA%getBlock('baba')
      sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms getBlock: baba(1,1,1,1)'
      else
        write(*,1010) 'spin/ms getBlock: baba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 20-24: spin+spin same-type, symm+symm → symm output.
!     ERI A (symm): aaaa=8, bbbb=9, aabb/abba/abab=10/11/12
!     ERI B (symm): aaaa=2, bbbb=3, aabb/abba/abab=4/6/7
!--------------------------------------------------------------------
!
!     Test 20: spin/coulomb_symm + spin/coulomb_symm → symm-coulomb
!     Expected: aaaa=10, bbbb=12, aabb=14
!
      call mqc_twoeris_allocate(eriA,'symm-coulomb','spin',tC,tD,tE)
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tA,tB,tAB)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+coulomb_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+coulomb_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('bbbb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+coulomb_s: bbbb(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+coulomb_s: bbbb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+coulomb_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+coulomb_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 21: spin/exchange_symm + spin/exchange_symm → symm-exchange
!     Expected: aaaa=10, bbbb=12, abba=17
!
      call mqc_twoeris_allocate(eriA,'symm-exchange','spin',tC,tD,tF)
      call mqc_twoeris_allocate(eriB,'symm-exchange','spin',tA,tB,tABBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+exchange_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+exchange_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+exchange_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+exchange_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 22: spin/inverted_symm + spin/inverted_symm → symm-inverted
!     Expected: aaaa=10, bbbb=12, abab=19
!
      call mqc_twoeris_allocate(eriA,'symm-inverted','spin',tC,tD,tG)
      call mqc_twoeris_allocate(eriB,'symm-inverted','spin',tA,tB,tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/inverted_s+inverted_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/inverted_s+inverted_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 19.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/inverted_s+inverted_s: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/inverted_s+inverted_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 23: spin/doublebar_symm + spin/doublebar_symm → symm-doublebar
!     Expected: aaaa=10, bbbb=12, aabb=14, abba=17
!
      call mqc_twoeris_allocate(eriA,'symm-doublebar','spin',tC,tD,tE,abba=tF)
      call mqc_twoeris_allocate(eriB,'symm-doublebar','spin',tA,tB,tAB,abba=tABBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+doublebar_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+doublebar_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+doublebar_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+doublebar_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+doublebar_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+doublebar_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 24: spin/ms_symm + spin/ms_symm → symm-ms
!     Expected: aaaa=10, bbbb=12, aabb=14, abba=17, abab=19
!
      call mqc_twoeris_allocate(eriA,'symm-ms','spin',tC,tD,tE,abba=tF,abab=tG)
      call mqc_twoeris_allocate(eriB,'symm-ms','spin',tA,tB,tAB,abba=tABBA,abab=tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_s+ms_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_s+ms_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_s+ms_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_s+ms_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_s+ms_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_s+ms_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 19.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_s+ms_s: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_s+ms_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 25-29: spin+spin same-type, full+full → full output.
!     ERI A (full): aaaa=8,bbbb=9; aabb/bbaa or abba/baab or abab/baba = 10/11/12
!     ERI B (full): aaaa=2,bbbb=3; aabb/bbaa or abba/baab or abab/baba = 4/6/7
!     For full-exchange/inverted, the 3rd/4th positional args store
!     abba/baab and abab/baba respectively (in alphaBeta/betaAlpha).
!--------------------------------------------------------------------
!
!     Test 25: spin/coulomb_full + spin/coulomb_full → full-coulomb
!     Expected: aaaa=10, bbbb=12, aabb=14, bbaa=14
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','spin',tC,tD,tE,tE)
      call mqc_twoeris_allocate(eriB,'full-coulomb','spin',tA,tB,tAB,tAB)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_f+coulomb_f: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_f+coulomb_f: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_f+coulomb_f: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_f+coulomb_f: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('bbaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_f+coulomb_f: bbaa(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_f+coulomb_f: bbaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 26: spin/exchange_full + spin/exchange_full → full-exchange
!     Positional 3/4 store abba/baab. Expected: aaaa=10, abba=17, baab=17
!
      call mqc_twoeris_allocate(eriA,'full-exchange','spin',tC,tD,tF,tF)
      call mqc_twoeris_allocate(eriB,'full-exchange','spin',tA,tB,tABBA,tABBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_f+exchange_f: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_f+exchange_f: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_f+exchange_f: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_f+exchange_f: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_f+exchange_f: baab(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_f+exchange_f: baab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 27: spin/inverted_full + spin/inverted_full → full-inverted
!     Positional 3/4 store abab/baba. Expected: aaaa=10, abab=19, baba=19
!
      call mqc_twoeris_allocate(eriA,'full-inverted','spin',tC,tD,tG,tG)
      call mqc_twoeris_allocate(eriB,'full-inverted','spin',tA,tB,tBBBA,tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/inverted_f+inverted_f: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/inverted_f+inverted_f: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 19.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/inverted_f+inverted_f: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/inverted_f+inverted_f: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 19.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/inverted_f+inverted_f: baba(1,1,1,1)'
      else
        write(*,1010) 'spin/inverted_f+inverted_f: baba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 28: spin/doublebar_full + spin/doublebar_full → full-doublebar
!     Expected: aaaa=10, bbbb=12, aabb=14, bbaa=14, abba=17, baab=17
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','spin',tC,tD,tE,tE,abba=tF,baab=tF)
      call mqc_twoeris_allocate(eriB,'full-doublebar','spin',tA,tB,tAB,tAB,abba=tABBA,baab=tABBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_f+doublebar_f: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_f+doublebar_f: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_f+doublebar_f: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_f+doublebar_f: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_f+doublebar_f: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_f+doublebar_f: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_f+doublebar_f: baab(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_f+doublebar_f: baab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 29: spin/ms_full + spin/ms_full → full-ms
!     Expected: aaaa=10, aabb=14, bbaa=14, abba=17, baab=17, abab=19, baba=19
!
      call mqc_twoeris_allocate(eriA,'full-ms','spin',tC,tD,tE,tE,abba=tF,baab=tF,abab=tG,baba=tG)
      call mqc_twoeris_allocate(eriB,'full-ms','spin',tA,tB,tAB,tAB,abba=tABBA,baab=tABBA,&
        abab=tBBBA,baba=tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_f+ms_f: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_f+ms_f: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_f+ms_f: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_f+ms_f: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_f+ms_f: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_f+ms_f: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 19.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_f+ms_f: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_f+ms_f: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 19.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/ms_f+ms_f: baba(1,1,1,1)'
      else
        write(*,1010) 'spin/ms_f+ms_f: baba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 30-36: spin+spin cross-type, symm inputs.
!     A: aaaa=8,bbbb=9,aabb/abba/abab=10/11/12
!     B: aaaa=2,bbbb=3,aabb/abba/abab=4/6/7
!--------------------------------------------------------------------
!
!     Test 30: spin/coulomb_symm(A) + spin/exchange_symm(B) → symm-doublebar
!     aaaa=10, bbbb=12, aabb=A.aabb=10, abba=B.abba=6
!
      call mqc_twoeris_allocate(eriA,'symm-coulomb','spin',tC,tD,tE)
      call mqc_twoeris_allocate(eriB,'symm-exchange','spin',tA,tB,tABBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+exchange_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+exchange_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+exchange_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+exchange_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+exchange_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+exchange_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 31: spin/coulomb_symm(A) + spin/inverted_symm(B) → full-ms
!     aaaa=10, bbbb=12, aabb=A.aabb=10, abab=B.abab=7
!
      call mqc_twoeris_allocate(eriA,'symm-coulomb','spin',tC,tD,tE)
      call mqc_twoeris_allocate(eriB,'symm-inverted','spin',tA,tB,tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+inverted_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+inverted_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+inverted_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+inverted_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+inverted_s: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+inverted_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 32: spin/coulomb_symm(A) + spin/ms_symm(B) → symm-ms
!     aaaa=10, aabb=A.aabb+B.aabb=14, abba=B.abba=6, abab=B.abab=7
!
      call mqc_twoeris_allocate(eriA,'symm-coulomb','spin',tC,tD,tE)
      call mqc_twoeris_allocate(eriB,'symm-ms','spin',tA,tB,tAB,abba=tABBA,abab=tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+ms_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+ms_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+ms_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+ms_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+ms_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+ms_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+ms_s: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+ms_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 33: spin/exchange_symm(A) + spin/inverted_symm(B) → full-ms
!     aaaa=10, abba=A.abba=11, abab=B.abab=7
!
      call mqc_twoeris_allocate(eriA,'symm-exchange','spin',tC,tD,tF)
      call mqc_twoeris_allocate(eriB,'symm-inverted','spin',tA,tB,tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+inverted_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+inverted_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+inverted_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+inverted_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+inverted_s: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+inverted_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 34: spin/exchange_symm(A) + spin/ms_symm(B) → symm-ms
!     aaaa=10, aabb=B.aabb=4, abba=A.abba+B.abba=17, abab=B.abab=7
!
      call mqc_twoeris_allocate(eriA,'symm-exchange','spin',tC,tD,tF)
      call mqc_twoeris_allocate(eriB,'symm-ms','spin',tA,tB,tAB,abba=tABBA,abab=tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 4.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+ms_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+ms_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+ms_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+ms_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+ms_s: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+ms_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 35: spin/doublebar_symm(A) + spin/inverted_symm(B) → full-ms
!     aaaa=10, aabb=A.aabb=10, abba=A.abba=11, abab=B.abab=7
!
      call mqc_twoeris_allocate(eriA,'symm-doublebar','spin',tC,tD,tE,abba=tF)
      call mqc_twoeris_allocate(eriB,'symm-inverted','spin',tA,tB,tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+inverted_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+inverted_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+inverted_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+inverted_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+inverted_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+inverted_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+inverted_s: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+inverted_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 36: spin/doublebar_symm(A) + spin/ms_symm(B) → symm-ms
!     aaaa=10, aabb=A.aabb+B.aabb=14, abba=A.abba+B.abba=17, abab=B.abab=7
!
      call mqc_twoeris_allocate(eriA,'symm-doublebar','spin',tC,tD,tE,abba=tF)
      call mqc_twoeris_allocate(eriB,'symm-ms','spin',tA,tB,tAB,abba=tABBA,abab=tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+ms_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+ms_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 14.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+ms_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+ms_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 17.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+ms_s: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+ms_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/doublebar_s+ms_s: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/doublebar_s+ms_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 37-43: space+spin cross-type combinations not yet covered.
!--------------------------------------------------------------------
!
!     Test 37: space/coulomb(aaaa=2) + spin/exchange_symm(8,9,abba=11) → full-doublebar
!     aaaa=10, bbbb=11, aabb=2 (space/coulomb), abba=11 (spin/exchange)
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-exchange','spin',tC,tD,tF)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/exchange_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/exchange_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/exchange_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/exchange_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/exchange_s: abba(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/exchange_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('bbaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/exchange_s: bbaa(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/exchange_s: bbaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/exchange_s: baab(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/exchange_s: baab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 38: space/exchange(aaaa=2) + spin/coulomb_symm(8,9,aabb=10) → symm-doublebar
!     aaaa=10, bbbb=11, aabb=spin/coulomb.aabb=10, abba=space/exchange.abba=2
!
      call mqc_twoeris_allocate(eriA,'full-exchange','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tC,tD,tE)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/exchange+spin/coulomb_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/exchange+spin/coulomb_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/exchange+spin/coulomb_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/exchange+spin/coulomb_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/exchange+spin/coulomb_s: abba(1,1,1,1)'
      else
        write(*,1010) 'space/exchange+spin/coulomb_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 39: space/doublebar(aaaa=2,aabb=4,abba=6) + spin/inverted_symm(8,9,abab=7)
!     → full-ms: aaaa=10, aabb=4, abba=6, abab=7
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','space',tA,alphaBeta=tAB,abba=tABBA)
      call mqc_twoeris_allocate(eriB,'symm-inverted','spin',tC,tD,tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+spin/inverted_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+spin/inverted_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 4.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+spin/inverted_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+spin/inverted_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+spin/inverted_s: abba(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+spin/inverted_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/doublebar+spin/inverted_s: abab(1,1,1,1)'
      else
        write(*,1010) 'space/doublebar+spin/inverted_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 40: space/ms(aaaa=2) + spin/coulomb_symm(8,9,aabb=10) → symm-ms
!     aaaa=10, bbbb=11, aabb=2+10=12, abba=2, abab=2
!
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tC,tD,tE)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/coulomb_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/coulomb_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/coulomb_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/coulomb_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+spin/coulomb_s: abab(1,1,1,1)'
      else
        write(*,1010) 'space/ms+spin/coulomb_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 41: space/inverted(aaaa=2) + spin/coulomb_symm(8,9,aabb=10) → symm-ms
!     aaaa=10, bbbb=11, aabb=spin/coulomb.aabb=10, abab=space/inverted.abab=2
!
      call mqc_twoeris_allocate(eriA,'full-inverted','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tC,tD,tE)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted+spin/coulomb_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/inverted+spin/coulomb_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted+spin/coulomb_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/inverted+spin/coulomb_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted+spin/coulomb_s: abab(1,1,1,1)'
      else
        write(*,1010) 'space/inverted+spin/coulomb_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 42-43: space+space cross-type with inverted and ms.
!--------------------------------------------------------------------
!
!     Test 42: space/inverted(aaaa=2) + space/coulomb(aaaa=3) → symm-ms spin
!     aaaa=5, aabb=0+3=3 (inverted.aabb=0, coulomb.aabb=3), abab=2+0=2
!
      call mqc_twoeris_allocate(eriA,'full-inverted','space',tA)
      call mqc_twoeris_allocate(eriB,'full-coulomb','space',tB)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 5.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted+space/coulomb: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/inverted+space/coulomb: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 3.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted+space/coulomb: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/inverted+space/coulomb: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/inverted+space/coulomb: abab(1,1,1,1)'
      else
        write(*,1010) 'space/inverted+space/coulomb: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 43: space/ms(aaaa=2) + space/inverted(aaaa=3) → symm-ms spin
!     aaaa=5, aabb=2+0=2, abba=2+0=2, abab=2+3=5
!
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      call mqc_twoeris_allocate(eriB,'full-inverted','space',tB)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 5.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+space/inverted: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/ms+space/inverted: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+space/inverted: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/ms+space/inverted: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 5.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/ms+space/inverted: abab(1,1,1,1)'
      else
        write(*,1010) 'space/ms+space/inverted: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!
!--------------------------------------------------------------------
!     Tests 44-45: space/coulomb + spin/inverted and spin/ms (symm).
!     These hit code paths that previously crashed due to betaAlpha
!     being passed as the 4th positional arg to symm-ms allocate.
!--------------------------------------------------------------------
!
!     Test 44: space/coulomb(aaaa=2) + spin/inverted_symm(8,9,abab=7) → symm-ms
!     aaaa=10, bbbb=11, aabb=2 (space/coulomb), abab=7 (spin/inverted)
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-inverted','spin',tC,tD,tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/inverted_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/inverted_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/inverted_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/inverted_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/inverted_s: abab(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/inverted_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 45: space/coulomb(aaaa=2) + spin/ms_symm(8,9,aabb=4,abba=6,abab=7) → symm-ms
!     aaaa=10, bbbb=11, aabb=2+4=6, abba=6 (spin/ms), abab=7 (spin/ms)
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-ms','spin',tC,tD,tAB,abba=tABBA,abab=tBBBA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/ms_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/ms_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/ms_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/ms_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/ms_s: abba(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/ms_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'space/coulomb+spin/ms_s: abab(1,1,1,1)'
      else
        write(*,1010) 'space/coulomb+spin/ms_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!
!--------------------------------------------------------------------
!     Tests 46-50: spin+space cross-type combinations (ERIA=spin, ERIB=space)
!     covering exchange/inverted for spin and coulomb/exchange/inverted for space.
!--------------------------------------------------------------------
!
!     Test 46: spin/coulomb_symm(8,9,10) + space/exchange(2) → symm-doublebar
!     aaaa=10, aabb=10, abba=2
!
      call mqc_twoeris_allocate(eriA,'symm-coulomb','spin',tC,tD,tE)
      call mqc_twoeris_allocate(eriB,'full-exchange','space',tA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+space/exchange: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+space/exchange: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+space/exchange: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+space/exchange: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+space/exchange: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+space/exchange: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 47: spin/coulomb_symm(8,9,10) + space/inverted(2) → symm-ms
!     aaaa=10, aabb=10, abab=2
!
      call mqc_twoeris_allocate(eriA,'symm-coulomb','spin',tC,tD,tE)
      call mqc_twoeris_allocate(eriB,'full-inverted','space',tA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+space/inverted: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+space/inverted: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+space/inverted: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+space/inverted: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/coulomb_s+space/inverted: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/coulomb_s+space/inverted: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 48: spin/exchange_symm(8,9,abba=11) + space/coulomb(2) → symm-doublebar
!     aaaa=10, aabb=2, abba=11
!
      call mqc_twoeris_allocate(eriA,'symm-exchange','spin',tC,tD,tF)
      call mqc_twoeris_allocate(eriB,'full-coulomb','space',tA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+space/coulomb: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+space/coulomb: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+space/coulomb: aabb(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+space/coulomb: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+space/coulomb: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+space/coulomb: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 49: spin/exchange_symm(8,9,abba=11) + space/exchange(2) → symm-exchange
!     aaaa=10, abba=13
!
      call mqc_twoeris_allocate(eriA,'symm-exchange','spin',tC,tD,tF)
      call mqc_twoeris_allocate(eriB,'full-exchange','space',tA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+space/exchange: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+space/exchange: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 13.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+space/exchange: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+space/exchange: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 50: spin/exchange_symm(8,9,abba=11) + space/inverted(2) → symm-ms
!     aaaa=10, abba=11, abab=2
!
      call mqc_twoeris_allocate(eriA,'symm-exchange','spin',tC,tD,tF)
      call mqc_twoeris_allocate(eriB,'full-inverted','space',tA)
      eriC = eriA + eriB
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+space/inverted: aaaa(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+space/inverted: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+space/inverted: abba(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+space/inverted: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'spin/exchange_s+space/inverted: abab(1,1,1,1)'
      else
        write(*,1010) 'spin/exchange_s+space/inverted: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!
!--------------------------------------------------------------------
!     Tests 51-58: Dagger (mqc_twoeris_conjugate_transpose) for spin
!     intTypes. Tests 52-56 cover the fixed bug: the old code read
!     getBlock('bbaa') for all intTypes, which returns 0 for exchange
!     and inverted, leaving all cross-spin blocks zero.
!--------------------------------------------------------------------
!
!     Test 51: spin/coulomb symm Dagger — regression: aabb preserved
!     Input: aaaa=8, bbbb=9, aabb=10; output aabb=dagger(bbaa)=dagger(alphaBeta)=10
!
      call mqc_twoeris_allocate(eriA,'symm-coulomb','spin',tC,tD,tE)
      eriC = Dagger(eriA)
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 8.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/coulomb_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/coulomb_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/coulomb_s: aabb(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/coulomb_s: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 52: spin/coulomb full Dagger — shows aabb<->bbaa swap
!     Input: aaaa=8, bbbb=9, aabb=10, bbaa=11; output aabb=11, bbaa=10
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','spin',tC,tD,tE,tF)
      eriC = Dagger(eriA)
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/coulomb_f: aabb(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/coulomb_f: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('bbaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/coulomb_f: bbaa(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/coulomb_f: bbaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 53: spin/exchange symm Dagger — was broken (returned 0 for abba)
!     Input: aaaa=8, bbbb=9, abba=11; abba is self-conjugate so output abba=11
!
      call mqc_twoeris_allocate(eriA,'symm-exchange','spin',tC,tD,tF)
      eriC = Dagger(eriA)
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 8.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/exchange_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/exchange_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/exchange_s: abba(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/exchange_s: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 54: spin/exchange full Dagger — was broken (returned 0 for abba/baab)
!     Input: aaaa=8, bbbb=9, abba=11, baab=12; both self-conjugate
!
      call mqc_twoeris_allocate(eriA,'full-exchange','spin',tC,tD,tF,tG)
      eriC = Dagger(eriA)
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/exchange_f: abba(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/exchange_f: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/exchange_f: baab(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/exchange_f: baab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 55: spin/inverted symm Dagger — was broken (returned 0 for abab)
!     Input: aaaa=8, bbbb=9, abab=12; for symm baba=abab so abab preserved
!
      call mqc_twoeris_allocate(eriA,'symm-inverted','spin',tC,tD,tG)
      eriC = Dagger(eriA)
      blk = eriC%getBlock('aaaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 8.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/inverted_s: aaaa(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/inverted_s: aaaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/inverted_s: abab(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/inverted_s: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 56: spin/inverted full Dagger — was broken; shows abab<->baba swap
!     Input: aaaa=8, bbbb=9, abab=12, baba=6; output abab=6, baba=12
!
      call mqc_twoeris_allocate(eriA,'full-inverted','spin',tC,tD,tG,tABBA)
      eriC = Dagger(eriA)
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/inverted_f: abab(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/inverted_f: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 12.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/inverted_f: baba(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/inverted_f: baba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 57: spin/doublebar full Dagger — aabb<->bbaa swap + abba/baab self
!     Input: aaaa=8, bbbb=9, aabb=10, bbaa=11, abba=6, baab=7
!     Output: aabb=11, bbaa=10, abba=6, baab=7
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','spin',tC,tD,tE,tF,&
                                abba=tABBA,baab=tBBBA)
      eriC = Dagger(eriA)
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/doublebar_f: aabb(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/doublebar_f: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('bbaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/doublebar_f: bbaa(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/doublebar_f: bbaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/doublebar_f: abba(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/doublebar_f: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/doublebar_f: baab(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/doublebar_f: baab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 58: spin/ms full Dagger — all 6 cross-spin blocks
!     Input: aaaa=8, bbbb=9, aabb=10, bbaa=11, abab=6, abba=7, baab=2, baba=3
!     Output: aabb=11, bbaa=10, abba=7(self), baab=2(self), abab=3(swap), baba=6(swap)
!
      call mqc_twoeris_allocate(eriA,'full-ms','spin',tC,tD,tE,tF,tABBA,tBBBA,tA,tB)
      eriC = Dagger(eriA)
      blk = eriC%getBlock('aabb') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 11.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/ms_f: aabb(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/ms_f: aabb(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('bbaa') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 10.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/ms_f: bbaa(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/ms_f: bbaa(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 7.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/ms_f: abba(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/ms_f: abba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 2.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/ms_f: baab(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/ms_f: baab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('abab') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 3.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/ms_f: abab(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/ms_f: abab(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
      blk = eriC%getBlock('baba') ; sc = blk%at(1,1,1,1)
      got = sc%rval() ; expected = 6.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'Dagger spin/ms_f: baba(1,1,1,1)'
      else
        write(*,1010) 'Dagger spin/ms_f: baba(1,1,1,1)',got,expected
        nfail = nfail+1
      endIf
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      end program unitTest_twoERIs
