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
      type(mqc_r4tensor)::tSO
      type(mqc_r4tensor)::blk
      type(mqc_r4tensor),dimension(16)::oracleBlocks
      type(mqc_scalar)::sc
      type(mqc_twoERIs)::eriA,eriB,eriC
      type(mqc_twoERIs)::eriAO,eriX,eriHalf,eriPartial,eriTrans,eriTransLR,eriFull,eriSymm
      type(mqc_matrix)::mA,mB,mC,mD,outBlk
      type(mqc_matrix)::moAA,moBB,moAB,moBA,moFull,moAA2,moBB2,moAB2,moBA2,moFull2
      type(mqc_scf_integral)::intSp,intSpin,intGen,intOut
      type(mqc_scf_integral)::generalMO,generalMO2
      real(kind=real64)::got,expected
      real(kind=real64)::partialError,transError,partialPairError,transPairError,fullError,atError
      real(kind=real64)::v1,v2,v3,v4,v5,v6,v7
      real(kind=real64)::v8,v9,v10,v11,v12
      complex(kind=real64)::directValue,partialValue,transValue
      complex(kind=real64)::coef1,coef2,coef3,coef4
      integer::p,q,r,s,mu,nu,lambda,sigma,spin1,spin2,iblock,jblock,iinteraction
      integer::colP,colQ,colR,colS
      character(len=4),dimension(16)::spinLabels
      character(len=9),dimension(3)::interactionLabels
      character(len=4)::spinLabel,exchangeLabel
!
!     Format Statements
!
 1000 format(1x,'PASS  ',a)
 1010 format(1x,'FAIL  ',a,'  got=',f16.10,'  expected=',f16.10)
 9000 format(/,1x,'Tests completed: ',i0,' failure(s)')
!
!
      nfail = 0
      v1=1.0d0 ; v2=2.0d0 ; v3=3.0d0 ; v4=4.0d0
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
!
!--------------------------------------------------------------------
!     Tests 59-64: mqc_eri_eri_contraction for newly implemented
!     space/doublebar and space/ms eris1 types.
!
!     Setup: eris1 is space type with aaaa=tA(2).
!            eris2 is spin type with aaaa=bbbb=tZ(0) so the
!            baseline aaaa*aaaa and aaaa*bbbb terms vanish,
!            isolating the cross-spin contributions from the new code.
!
!     Contraction sums: for n=2 with only (1,1,1,1) nonzero,
!     contraction(A,B) = A(1,1,1,1)*B(1,1,1,1).
!     For symm type, bXXX blocks equal the stored aXXX block.
!--------------------------------------------------------------------
!
!     Test 59: space/doublebar(aaaa=2) x spin/coulomb_symm(aabb=10),
!              label='coulomb'
!     New path: case('doublebar')/case('coulomb')
!     Expected: eris1.aaaa*(eris2.aabb+eris2.bbaa) = 2*(10+10) = 40
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tZ,tZ,tE)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 40.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction space/doublebar x spin/coulomb (J)'
      else
        write(*,1010) 'contraction space/doublebar x spin/coulomb (J)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 60: space/doublebar(aaaa=2) x spin/exchange_symm(abba=11),
!              label='coulomb'
!     New path: case('doublebar')/case('exchange')
!     Expected: eris1.aaaa*(eris2.abba+eris2.baab) = 2*(11+11) = 44
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-exchange','spin',tZ,tZ,tF)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 44.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction space/doublebar x spin/exchange (J)'
      else
        write(*,1010) 'contraction space/doublebar x spin/exchange (J)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 61: space/doublebar(aaaa=2) x spin/doublebar_symm(aabb=10,abba=11),
!              label='coulomb'
!     New path: case('doublebar')/case('doublebar')
!     Expected: eris1.aaaa*(aabb+bbaa+abba+baab) = 2*(10+10+11+11) = 84
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-doublebar','spin',tZ,tZ,tE,abba=tF)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 84.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction space/doublebar x spin/doublebar (J)'
      else
        write(*,1010) 'contraction space/doublebar x spin/doublebar (J)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 62: space/ms(aaaa=2) x spin/inverted_symm(abab=12),
!              label='coulomb'
!     New path: case('ms')/case('inverted') — previously would error
!     Expected: eris1.aaaa*(eris2.abab+eris2.baba) = 2*(12+12) = 48
!
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-inverted','spin',tZ,tZ,tG)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 48.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction space/ms x spin/inverted (J)'
      else
        write(*,1010) 'contraction space/ms x spin/inverted (J)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 63: space/ms(aaaa=2) x spin/ms_symm(aabb=10,abba=11,abab=12),
!              label='coulomb'
!     New path: case('ms')/case('ms') — all 6 cross-spin channels
!     Expected: eris1.aaaa*(aabb+bbaa+abba+baab+abab+baba)
!             = 2*(10+10+11+11+12+12) = 132
!
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-ms','spin',tZ,tZ,tE,abba=tF,abab=tG)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 132.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction space/ms x spin/ms (J)'
      else
        write(*,1010) 'contraction space/ms x spin/ms (J)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 64: space/ms(aaaa=2) x spin/inverted_symm(abab=12),
!              label='exchange'
!     New path: case('ms')/case('inverted') with exchange label
!     Expected: -(eris1.aaaa*(abab+baba)) = -(2*(12+12)) = -48
!
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      call mqc_twoeris_allocate(eriB,'symm-inverted','spin',tZ,tZ,tG)
      sc = mqc_eri_eri_contraction(eriA,eriB,'exchange')
      got = sc%rval() ; expected = -48.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction space/ms x spin/inverted (K)'
      else
        write(*,1010) 'contraction space/ms x spin/inverted (K)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 65: space/doublebar(aaaa=2) x general(aabb=10,bbaa=10,abba=11,baab=11),
!              label='coulomb'
!     New path: space(doublebar) x general, case('doublebar')
!     Expected: eris1.aaaa*(aabb+bbaa+abba+baab) = 2*(10+10+11+11) = 84
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','space',tA)
      call mqc_twoeris_allocate(eriB,'full-ms','general',tZ,tZ,tE,tE,abba=tF,baab=tF)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 84.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction space/doublebar x general/ms (J)'
      else
        write(*,1010) 'contraction space/doublebar x general/ms (J)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 66: space/ms(aaaa=2) x general(aabb=10,bbaa=10,abba=11,baab=11,abab=12,baba=12),
!              label='coulomb'
!     New path: space(ms) x general, case('ms')
!     Expected: eris1.aaaa*(aabb+bbaa+abba+baab+abab+baba) = 2*(10+10+11+11+12+12) = 132
!
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      call mqc_twoeris_allocate(eriB,'full-ms','general',tZ,tZ,tE,tE,abba=tF,baab=tF,  &
        abab=tG,baba=tG)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 132.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction space/ms x general/ms (J)'
      else
        write(*,1010) 'contraction space/ms x general/ms (J)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 67: general(aabb=10,bbaa=10,abba=11,baab=11) x space/doublebar(aaaa=2),
!              label='coulomb'
!     New path: general x space(doublebar), case('doublebar')
!     Expected: (aabb+bbaa+abba+baab)*eris2.aaaa = (10+10+11+11)*2 = 84
!
      call mqc_twoeris_allocate(eriA,'full-ms','general',tZ,tZ,tE,tE,abba=tF,baab=tF)
      call mqc_twoeris_allocate(eriB,'full-doublebar','space',tA)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 84.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction general/ms x space/doublebar (J)'
      else
        write(*,1010) 'contraction general/ms x space/doublebar (J)',got,expected
        nfail = nfail+1
      endIf
!
!     Test 68: general(aabb=10,bbaa=10,abba=11,baab=11,abab=12,baba=12) x space/ms(aaaa=2),
!              label='coulomb'
!     New path: general x space(ms), case('ms')
!     Expected: (aabb+bbaa+abba+baab+abab+baba)*eris2.aaaa = (10+10+11+11+12+12)*2 = 132
!
      call mqc_twoeris_allocate(eriA,'full-ms','general',tZ,tZ,tE,tE,abba=tF,baab=tF,  &
        abab=tG,baba=tG)
      call mqc_twoeris_allocate(eriB,'full-ms','space',tA)
      sc = mqc_eri_eri_contraction(eriA,eriB,'coulomb')
      got = sc%rval() ; expected = 132.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'contraction general/ms x space/ms (J)'
      else
        write(*,1010) 'contraction general/ms x space/ms (J)',got,expected
        nfail = nfail+1
      endIf
!
!
!--------------------------------------------------------------------
!     Tests 69-128: twoERIs%at element access.
!
!     blk: v5 at (1,1,1,2), v6 at (1,2,1,1).
!     At index (1,1,1,2): coulomb reads alpha(1,1,1,2)=v5;
!     exchange reads alpha(1,2,1,1)=v6.
!--------------------------------------------------------------------
!
      call mqc_allocate_r4tensor(n,n,n,n,blk,'real','full')
      call blk%put(v5,1,1,1,2)
      call blk%put(v6,1,2,1,1)
!
!--------------------------------------------------------------------
!     Tests 69-78: space integralType.
!     All space intTypes behave identically in %at (all read from
!     the alpha block).  Tests cover all 5 intTypes and both storageTypes.
!--------------------------------------------------------------------
!
!     Test 69: full-coulomb space, aaaa/J
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',blk)
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v5
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-coul-full/aaaa/J'
      else
        write(*,1010) 'at/space-coul-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 70: full-coulomb space, aaaa/K at (1,1,1,2) → alpha(1,2,1,1)=v6
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='exchange')
      got = sc%rval() ; expected = v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-coul-full/aaaa/K'
      else
        write(*,1010) 'at/space-coul-full/aaaa/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 71: full-exchange space, aaaa/J
      call mqc_twoeris_allocate(eriA,'full-exchange','space',blk)
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v5
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-exch-full/aaaa/J'
      else
        write(*,1010) 'at/space-exch-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 72: full-exchange space, aaaa/K
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='exchange')
      got = sc%rval() ; expected = v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-exch-full/aaaa/K'
      else
        write(*,1010) 'at/space-exch-full/aaaa/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 73: full-doublebar space, aaaa/JK at (1,1,1,2) → v5-v6
      call mqc_twoeris_allocate(eriA,'full-doublebar','space',blk)
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='doublebar')
      got = sc%rval() ; expected = v5-v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-dbar-full/aaaa/JK'
      else
        write(*,1010) 'at/space-dbar-full/aaaa/JK',got,expected ; nfail=nfail+1
      endIf
!
!     Test 74: full-ms space, aaaa/J at (1,1,1,1) → v2
      call mqc_twoeris_allocate(eriA,'full-ms','space',tA)
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v2
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-ms-full/aaaa/J'
      else
        write(*,1010) 'at/space-ms-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 75: full-inverted space, aaaa/J at (1,1,1,1) → v2
      call mqc_twoeris_allocate(eriA,'full-inverted','space',tA)
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v2
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-inv-full/aaaa/J'
      else
        write(*,1010) 'at/space-inv-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 76: symm-coulomb space, aaaa/J at (1,1,1,1) → v2
      call mqc_twoeris_allocate(eriA,'symm-coulomb','space',tA)
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v2
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-coul-symm/aaaa/J'
      else
        write(*,1010) 'at/space-coul-symm/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 77: symm-exchange space, aaaa/K at (1,1,1,2) → v6
      call mqc_twoeris_allocate(eriA,'symm-exchange','space',blk)
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='exchange')
      got = sc%rval() ; expected = v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-exch-symm/aaaa/K'
      else
        write(*,1010) 'at/space-exch-symm/aaaa/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 78: symm-ms space, aaaa/JK at (1,1,1,2) → v5-v6
      call mqc_twoeris_allocate(eriA,'symm-ms','space',blk)
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='doublebar')
      got = sc%rval() ; expected = v5-v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-ms-symm/aaaa/JK'
      else
        write(*,1010) 'at/space-ms-symm/aaaa/JK',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 79-98: spin integralType.
!     tC=v8(aaaa), tD=v9(bbbb), tE=v10(aabb), tF=v11(abba),
!     tG=v12(abab).  blk: v5@(1,1,1,2), v6@(1,2,1,1).
!--------------------------------------------------------------------
!
!     Tests 79-81: full-coulomb spin.
!     ERI: alpha=tC, beta=tD, alphaBeta=tE (aabb J integrals).
      call mqc_twoeris_allocate(eriA,'full-coulomb','spin',tC,tD,tE)
!
!     Test 79: aaaa/J at (1,1,1,1) → alpha(1,1,1,1)=v8
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v8
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-coul-full/aaaa/J'
      else
        write(*,1010) 'at/spin-coul-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 80: aabb/J at (1,1,1,1) → alphaBeta(1,1,1,1)=v10
      sc = eriA%at(1,1,1,1,spinBlockIn='aabb',interactionIn='coulomb')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-coul-full/aabb/J'
      else
        write(*,1010) 'at/spin-coul-full/aabb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 81: bbaa/J at (1,1,1,1) → betaAlpha absent (full, not symm) → 0
      sc = eriA%at(1,1,1,1,spinBlockIn='bbaa',interactionIn='coulomb')
      got = sc%rval() ; expected = 0.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-coul-full/bbaa/J=0'
      else
        write(*,1010) 'at/spin-coul-full/bbaa/J=0',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 82-83: symm-coulomb spin.
!     ERI: alpha=blk, beta=tD, alphaBeta=tE.
      call mqc_twoeris_allocate(eriA,'symm-coulomb','spin',blk,tD,alphaBeta=tE)
!
!     Test 82: aaaa/K at (1,1,1,2) → alpha(1,2,1,1)=v6
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='exchange')
      got = sc%rval() ; expected = v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-coul-symm/aaaa/K'
      else
        write(*,1010) 'at/spin-coul-symm/aaaa/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 83: bbaa/J at (1,1,1,1) → symm: alphaBeta(k,l,i,j)=(1,1,1,1)=v10
      sc = eriA%at(1,1,1,1,spinBlockIn='bbaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-coul-symm/bbaa/J'
      else
        write(*,1010) 'at/spin-coul-symm/bbaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 84-85: full-exchange spin.
!     ERI: alpha=tC, beta=tD, alphaBeta=tF (abba K integrals in alphaBeta slot).
      call mqc_twoeris_allocate(eriA,'full-exchange','spin',tC,tD,tF)
!
!     Test 84: aaaa/K at (1,1,1,1) → alpha(i,l,k,j)=alpha(1,1,1,1)=v8
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='exchange')
      got = sc%rval() ; expected = v8
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-exch-full/aaaa/K'
      else
        write(*,1010) 'at/spin-exch-full/aaaa/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 85: abba/K at (1,1,1,1) → alphaBeta(i,l,k,j)=(1,1,1,1)=v11
      sc = eriA%at(1,1,1,1,spinBlockIn='abba',interactionIn='exchange')
      got = sc%rval() ; expected = v11
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-exch-full/abba/K'
      else
        write(*,1010) 'at/spin-exch-full/abba/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 86: symm-exchange spin, baab/K at (1,1,1,1).
!     Symm derives baab from alphaBeta: alphaBeta(k,j,i,l)=(1,1,1,1)=v11.
      call mqc_twoeris_allocate(eriA,'symm-exchange','spin',tC,tD,tF)
      sc = eriA%at(1,1,1,1,spinBlockIn='baab',interactionIn='exchange')
      got = sc%rval() ; expected = v11
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-exch-symm/baab/K'
      else
        write(*,1010) 'at/spin-exch-symm/baab/K',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 87-89: full-doublebar spin.
!     ERI: alpha=blk, beta=tD, alphaBeta=tE (aabb), abba=tF.
      call mqc_twoeris_allocate(eriA,'full-doublebar','spin',blk,tD,  &
        alphaBeta=tE,abba=tF)
!
!     Test 87: aaaa/JK at (1,1,1,2) → alpha(1,1,1,2)-alpha(1,2,1,1)=v5-v6
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='doublebar')
      got = sc%rval() ; expected = v5-v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-dbar-full/aaaa/JK'
      else
        write(*,1010) 'at/spin-dbar-full/aaaa/JK',got,expected ; nfail=nfail+1
      endIf
!
!     Test 88: aabb/JK at (1,1,1,1) → alphaBeta(1,1,1,1)=v10 (no exchange term)
      sc = eriA%at(1,1,1,1,spinBlockIn='aabb',interactionIn='doublebar')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-dbar-full/aabb/JK'
      else
        write(*,1010) 'at/spin-dbar-full/aabb/JK',got,expected ; nfail=nfail+1
      endIf
!
!     Test 89: abba/JK at (1,1,1,1) → (-1)*alphaBeta(i,l,k,j)=-v10
      sc = eriA%at(1,1,1,1,spinBlockIn='abba',interactionIn='doublebar')
      got = sc%rval() ; expected = -v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-dbar-full/abba/JK'
      else
        write(*,1010) 'at/spin-dbar-full/abba/JK',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 90-93: full-ms and symm-ms spin.
!     ERI: alpha=tC, beta=tD, alphaBeta=tE (aabb), abba=tF, abab=tG.
      call mqc_twoeris_allocate(eriA,'full-ms','spin',tC,tD,alphaBeta=tE,  &
        abba=tF,abab=tG)
!
!     Test 90: aaaa/J at (1,1,1,1) → v8
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v8
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms-full/aaaa/J'
      else
        write(*,1010) 'at/spin-ms-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 91: bbbb/J at (1,1,1,1) → beta(1,1,1,1)=v9
      sc = eriA%at(1,1,1,1,spinBlockIn='bbbb',interactionIn='coulomb')
      got = sc%rval() ; expected = v9
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms-full/bbbb/J'
      else
        write(*,1010) 'at/spin-ms-full/bbbb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 92: aabb/J at (1,1,1,1) → alphaBeta(1,1,1,1)=v10
      sc = eriA%at(1,1,1,1,spinBlockIn='aabb',interactionIn='coulomb')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms-full/aabb/J'
      else
        write(*,1010) 'at/spin-ms-full/aabb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 93: symm-ms spin, bbaa/J at (1,1,1,1) → alphaBeta(k,l,i,j)=v10
      call mqc_twoeris_allocate(eriA,'symm-ms','spin',tC,tD,alphaBeta=tE,  &
        abba=tF,abab=tG)
      sc = eriA%at(1,1,1,1,spinBlockIn='bbaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms-symm/bbaa/J'
      else
        write(*,1010) 'at/spin-ms-symm/bbaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 94-95: full-inverted and symm-inverted spin.
!     Inverted ERI: alphaBeta stores abab data (tG=v12).
      call mqc_twoeris_allocate(eriA,'full-inverted','spin',tC,tD,tG)
!
!     Test 94: aaaa/J at (1,1,1,1) → alpha(1,1,1,1)=v8
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v8
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-inv-full/aaaa/J'
      else
        write(*,1010) 'at/spin-inv-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 95: aabb/J at (1,1,1,1) → alphaBeta(1,1,1,1)=v12 (abab data)
      sc = eriA%at(1,1,1,1,spinBlockIn='aabb',interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-inv-full/aabb/J'
      else
        write(*,1010) 'at/spin-inv-full/aabb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 96: symm-inverted spin, aabb/J → v12
      call mqc_twoeris_allocate(eriA,'symm-inverted','spin',tC,tD,tG)
      sc = eriA%at(1,1,1,1,spinBlockIn='aabb',interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-inv-symm/aabb/J'
      else
        write(*,1010) 'at/spin-inv-symm/aabb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 97-98: spin/ms abab BUG.
!     spin integralType has no abab/baba spinBlock cases; silently returns 0.
!     Both tests will FAIL until missing cases are added.
      call mqc_twoeris_allocate(eriA,'full-ms','spin',tC,tD,alphaBeta=tE,  &
        abba=tF,abab=tG)
!
!     Test 97: spinBlock='abab'/J → should be v12; BUG: returns 0
      sc = eriA%at(1,1,1,1,spinBlockIn='abab',interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms-full/abab/J'
      else
        write(*,1010) 'at/spin-ms-full/abab/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 98: no-spinBlock composite (1,3,1,3) → ABAB pattern → should be v12; BUG: 0
      sc = eriA%at(1,3,1,3,interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms-full/noSpinBlk-abab/J'
      else
        write(*,1010) 'at/spin-ms-full/noSpinBlk-abab/J',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 99-114: general integralType.
!--------------------------------------------------------------------
!
!     Tests 99-101: full-coulomb general and symm-coulomb general.
      call mqc_twoeris_allocate(eriA,'full-coulomb','general',tC,tD,tE)
!
!     Test 99: aaaa/J at (1,1,1,1) → alpha(1,1,1,1)=v8
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v8
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-coul-full/aaaa/J'
      else
        write(*,1010) 'at/gen-coul-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 100: aabb/J at (1,1,1,1) → alphaBeta(1,1,1,1)=v10
      sc = eriA%at(1,1,1,1,spinBlockIn='aabb',interactionIn='coulomb')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-coul-full/aabb/J'
      else
        write(*,1010) 'at/gen-coul-full/aabb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 101: symm-coulomb general, bbaa/J → alphaBeta(k,l,i,j)=v10
      call mqc_twoeris_allocate(eriA,'symm-coulomb','general',tC,tD,alphaBeta=tE)
      sc = eriA%at(1,1,1,1,spinBlockIn='bbaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-coul-symm/bbaa/J'
      else
        write(*,1010) 'at/gen-coul-symm/bbaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 102-104: full-exchange general and symm-exchange general.
      call mqc_twoeris_allocate(eriA,'full-exchange','general',tC,tD,abba=tF)
!
!     Test 102: aaaa/K at (1,1,1,1) → alpha(i,l,k,j)=(1,1,1,1)=v8
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='exchange')
      got = sc%rval() ; expected = v8
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-exch-full/aaaa/K'
      else
        write(*,1010) 'at/gen-exch-full/aaaa/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 103: abba/J at (1,1,1,1) → abba(1,1,1,1)=v11
      sc = eriA%at(1,1,1,1,spinBlockIn='abba',interactionIn='coulomb')
      got = sc%rval() ; expected = v11
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-exch-full/abba/J'
      else
        write(*,1010) 'at/gen-exch-full/abba/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 104: symm-exchange general, baab/J → abba(k,l,i,j)=v11
      call mqc_twoeris_allocate(eriA,'symm-exchange','general',tC,tD,abba=tF)
      sc = eriA%at(1,1,1,1,spinBlockIn='baab',interactionIn='coulomb')
      got = sc%rval() ; expected = v11
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-exch-symm/baab/J'
      else
        write(*,1010) 'at/gen-exch-symm/baab/J',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 105-107: full-doublebar and symm-doublebar general.
      call mqc_twoeris_allocate(eriA,'full-doublebar','general',blk,tD,  &
        alphaBeta=tE,abba=tF)
!
!     Test 105: aaaa/JK at (1,1,1,2) → alpha(1,1,1,2)-alpha(1,2,1,1)=v5-v6
      sc = eriA%at(1,1,1,2,spinBlockIn='aaaa',interactionIn='doublebar')
      got = sc%rval() ; expected = v5-v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-dbar-full/aaaa/JK'
      else
        write(*,1010) 'at/gen-dbar-full/aaaa/JK',got,expected ; nfail=nfail+1
      endIf
!
!     Test 106: aabb/JK at (1,1,1,1) → alphaBeta(1,1,1,1)-abba(i,l,k,j)=v10-v11
      sc = eriA%at(1,1,1,1,spinBlockIn='aabb',interactionIn='doublebar')
      got = sc%rval() ; expected = v10-v11
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-dbar-full/aabb/JK'
      else
        write(*,1010) 'at/gen-dbar-full/aabb/JK',got,expected ; nfail=nfail+1
      endIf
!
!     Test 107: symm-doublebar general, bbaa/JK → alphaBeta(k,l,i,j)-abba(k,j,i,l)=v10-v11
      call mqc_twoeris_allocate(eriA,'symm-doublebar','general',blk,tD,  &
        alphaBeta=tE,abba=tF)
      sc = eriA%at(1,1,1,1,spinBlockIn='bbaa',interactionIn='doublebar')
      got = sc%rval() ; expected = v10-v11
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-dbar-symm/bbaa/JK'
      else
        write(*,1010) 'at/gen-dbar-symm/bbaa/JK',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 108-111: full-ms and symm-ms general.
      call mqc_twoeris_allocate(eriA,'full-ms','general',tC,tD,tE,tE,  &
        abba=tF,baab=tF,abab=tG)
!
!     Test 108: abab/J at (1,1,1,1) → abab(1,1,1,1)=v12
      sc = eriA%at(1,1,1,1,spinBlockIn='abab',interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms-full/abab/J'
      else
        write(*,1010) 'at/gen-ms-full/abab/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 109: abab/K at (1,1,1,1) → abab(i,l,k,j)=(1,1,1,1)=v12
      sc = eriA%at(1,1,1,1,spinBlockIn='abab',interactionIn='exchange')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms-full/abab/K'
      else
        write(*,1010) 'at/gen-ms-full/abab/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 110: abab/JK at (1,1,1,1) → abab-abab(i,l,k,j)=v12-v12=0
      sc = eriA%at(1,1,1,1,spinBlockIn='abab',interactionIn='doublebar')
      got = sc%rval() ; expected = 0.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms-full/abab/JK=0'
      else
        write(*,1010) 'at/gen-ms-full/abab/JK=0',got,expected ; nfail=nfail+1
      endIf
!
!     Test 111: symm-ms general, baba/J → conjg(abab(j,i,l,k))=v12
      call mqc_twoeris_allocate(eriA,'symm-ms','general',tC,tD,alphaBeta=tE,  &
        abba=tF,abab=tG)
      sc = eriA%at(1,1,1,1,spinBlockIn='baba',interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms-symm/baba/J'
      else
        write(*,1010) 'at/gen-ms-symm/baba/J',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 112-114: full-inverted and symm-inverted general.
!     For general type, inverted uses a dedicated abab block (not alphaBeta).
      call mqc_twoeris_allocate(eriA,'full-inverted','general',tC,tD,abab=tG)
!
!     Test 112: aaaa/J at (1,1,1,1) → v8
      sc = eriA%at(1,1,1,1,spinBlockIn='aaaa',interactionIn='coulomb')
      got = sc%rval() ; expected = v8
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-inv-full/aaaa/J'
      else
        write(*,1010) 'at/gen-inv-full/aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 113: abab/J at (1,1,1,1) → abab(1,1,1,1)=v12
      sc = eriA%at(1,1,1,1,spinBlockIn='abab',interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-inv-full/abab/J'
      else
        write(*,1010) 'at/gen-inv-full/abab/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 114: symm-inverted general, baba/J → conjg(abab(j,i,l,k))=v12
      call mqc_twoeris_allocate(eriA,'symm-inverted','general',tC,tD,abab=tG)
      sc = eriA%at(1,1,1,1,spinBlockIn='baba',interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-inv-symm/baba/J'
      else
        write(*,1010) 'at/gen-inv-symm/baba/J',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 115-116: babb doublebar sign BUG (general, spinBlock and
!     no-spinBlock paths).
!     Expected: babb(i,j,k,l) - bbba(i,l,k,j) = v8-v9.
!     BUG: second assignment overwrites instead of subtracts.
!     Both tests FAIL until the minus sign is corrected.
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-ms','general',tZ,tZ,tZ,tZ,  &
        babb=tC,bbba=tD)
!
!     Test 115: spinBlock='babb', JK at (1,1,1,1) → v8-v9
      sc = eriA%at(1,1,1,1,spinBlockIn='babb',interactionIn='doublebar')
      got = sc%rval() ; expected = v8-v9
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms-full/spinBlk-babb/JK'
      else
        write(*,1010) 'at/gen-ms-full/spinBlk-babb/JK',got,expected ; nfail=nfail+1
      endIf
!
!     Test 116: no-spinBlock at composite (3,1,3,3) → BABB pattern → v8-v9
!     n=2: indices 1-2=alpha, 3-4=beta; (3,1,3,3) → (β,α,β,β)=babb.
      sc = eriA%at(3,1,3,3,interactionIn='doublebar')
      got = sc%rval() ; expected = v8-v9
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms-full/noSpinBlk-babb/JK'
      else
        write(*,1010) 'at/gen-ms-full/noSpinBlk-babb/JK',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 117-128: no-spinBlock path.
!     n=2: indices 1-2=alpha, 3-4=beta.
!--------------------------------------------------------------------
!
!     Tests 117-118: no-spinBlock, space integralType.
      call mqc_twoeris_allocate(eriA,'full-ms','space',blk)
!
!     Test 117: (1,1,1,2)/J → aaaa-pattern → alpha(1,1,1,2)=v5
      sc = eriA%at(1,1,1,2,interactionIn='coulomb')
      got = sc%rval() ; expected = v5
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-ms/noSpinBlk-aaaa/J'
      else
        write(*,1010) 'at/space-ms/noSpinBlk-aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 118: (1,1,1,2)/K → alpha(1,2,1,1)=v6
      sc = eriA%at(1,1,1,2,interactionIn='exchange')
      got = sc%rval() ; expected = v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/space-ms/noSpinBlk-aaaa/K'
      else
        write(*,1010) 'at/space-ms/noSpinBlk-aaaa/K',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 119-124: no-spinBlock, spin integralType (symm-ms).
      call mqc_twoeris_allocate(eriA,'symm-ms','spin',blk,tD,  &
        alphaBeta=tE,abba=tF,abab=tG)
!
!     Test 119: aaaa (1,1,1,2)/J → alpha(1,1,1,2)=v5
      sc = eriA%at(1,1,1,2,interactionIn='coulomb')
      got = sc%rval() ; expected = v5
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms/noSpinBlk-aaaa/J'
      else
        write(*,1010) 'at/spin-ms/noSpinBlk-aaaa/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 120: aaaa (1,1,1,2)/K → alpha(1,2,1,1)=v6
      sc = eriA%at(1,1,1,2,interactionIn='exchange')
      got = sc%rval() ; expected = v6
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms/noSpinBlk-aaaa/K'
      else
        write(*,1010) 'at/spin-ms/noSpinBlk-aaaa/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 121: bbbb (3,3,3,3)/J → beta(1,1,1,1)=v9
      sc = eriA%at(3,3,3,3,interactionIn='coulomb')
      got = sc%rval() ; expected = v9
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms/noSpinBlk-bbbb/J'
      else
        write(*,1010) 'at/spin-ms/noSpinBlk-bbbb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 122: aabb (1,1,3,3)/J → alphaBeta(1,1,1,1)=v10
      sc = eriA%at(1,1,3,3,interactionIn='coulomb')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms/noSpinBlk-aabb/J'
      else
        write(*,1010) 'at/spin-ms/noSpinBlk-aabb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 123: abba (1,3,3,1)/K → alphaBeta(i,l,k,j)=(1,1,1,1)=v10
      sc = eriA%at(1,3,3,1,interactionIn='exchange')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms/noSpinBlk-abba/K'
      else
        write(*,1010) 'at/spin-ms/noSpinBlk-abba/K',got,expected ; nfail=nfail+1
      endIf
!
!     Test 124: abab (1,3,1,3)/J → BUG: spin abab missing → returns 0 not v12
      sc = eriA%at(1,3,1,3,interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/spin-ms/noSpinBlk-abab/J'
      else
        write(*,1010) 'at/spin-ms/noSpinBlk-abab/J',got,expected ; nfail=nfail+1
      endIf
!
!     Tests 125-128: no-spinBlock, general integralType (full-ms).
      call mqc_twoeris_allocate(eriA,'full-ms','general',tC,tD,tE,tE,  &
        abba=tF,baab=tF,abab=tG)
!
!     Test 125: aabb (1,1,3,3)/J → alphaBeta(1,1,1,1)=v10
      sc = eriA%at(1,1,3,3,interactionIn='coulomb')
      got = sc%rval() ; expected = v10
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms/noSpinBlk-aabb/J'
      else
        write(*,1010) 'at/gen-ms/noSpinBlk-aabb/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 126: abba (1,3,3,1)/J → abba(1,1,1,1)=v11
      sc = eriA%at(1,3,3,1,interactionIn='coulomb')
      got = sc%rval() ; expected = v11
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms/noSpinBlk-abba/J'
      else
        write(*,1010) 'at/gen-ms/noSpinBlk-abba/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 127: abab (1,3,1,3)/J → abab(1,1,1,1)=v12
      sc = eriA%at(1,3,1,3,interactionIn='coulomb')
      got = sc%rval() ; expected = v12
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms/noSpinBlk-abab/J'
      else
        write(*,1010) 'at/gen-ms/noSpinBlk-abab/J',got,expected ; nfail=nfail+1
      endIf
!
!     Test 128: aabb (1,1,3,3)/JK → alphaBeta-abba(i,l,k,j)=v10-v11
      sc = eriA%at(1,1,3,3,interactionIn='doublebar')
      got = sc%rval() ; expected = v10-v11
      if(abs(got-expected).le.tol) then
        write(*,1000) 'at/gen-ms/noSpinBlk-aabb/JK'
      else
        write(*,1010) 'at/gen-ms/noSpinBlk-aabb/JK',got,expected ; nfail=nfail+1
      endIf
!
!
!
!--------------------------------------------------------------------
!     Tests 129-145: mqc_eri_integral_contraction
!
!     Setup:
!       mA=v2=2, mB=v3=3, mC=v4=4, mD=v5=5  (2x2 matrices, (1,1) only)
!       intSp   = space integral:   aa=mA
!       intSpin = spin  integral:   aa=mA, bb=mB
!       intGen  = general integral: aa=mA, bb=mB, ab=mC, ba=mD
!
!       eriA (space):     full-coulomb-space,    aaaa=tC(8)
!       eriB (spin symm): symm-coulomb-spin,     aaaa=8,bbbb=9,aabb=10(bbaa=10)
!       eriC (spin full): full-doublebar-spin,   aaaa=8,bbbb=9,aabb=10,bbaa=11
!       eriA (gen basic): full-ms-general,       aaaa=8,bbbb=9,aabb=bbaa=10,
!                                                abba=baab=11,abab=baba=12
!       eriB (gen+SF):    same + aaab=bbba=10
!
!     Contraction formula: for n=2 with only (1,1,1,1) nonzero in every
!     R4tensor and (1,1) only nonzero in every matrix, the contraction
!     of any R4tensor T with matrix M along any index pair gives:
!       result(1,1) = T(1,1,1,1) * M(1,1)
!--------------------------------------------------------------------
!
!     Build 2x2 integral matrices
!
      call mqc_allocate_matrix(n,n,mA,'real','StorFull')
      call mqc_allocate_matrix(n,n,mB,'real','StorFull')
      call mqc_allocate_matrix(n,n,mC,'real','StorFull')
      call mqc_allocate_matrix(n,n,mD,'real','StorFull')
      call mA%put(v2,1,1) ; call mB%put(v3,1,1)
      call mC%put(v4,1,1) ; call mD%put(v5,1,1)
!
!     Build SCF integrals
!
      call mqc_integral_allocate(intSp,  '','space',  mA)
      call mqc_integral_allocate(intSpin,'','spin',   mA,mB)
      call mqc_integral_allocate(intGen, '','general',mA,mB,mC,mD)
!
!--------------------------------------------------------------------
!     Tests 129-133: space ERI × {space, spin, general}
!--------------------------------------------------------------------
!
!     Test 129: space ERI × space integral, default ind (same=T, ia=1)
!     Code: aa = 2*R4(aaaa,M_aa) = 2*(8*2) = 32
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tC)
      intOut = Contraction(eriA,intSp)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 32.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x space/same=T: aa'
      else
        write(*,1010) 'IC/space x space/same=T: aa',got,expected ; nfail=nfail+1
      endIf
!
!     Test 130: space ERI × spin integral, same=T
!     Code: aa = R4(aaaa,M_aa) + R4(aaaa,M_bb) = 8*2+8*3 = 40
!
      intOut = Contraction(eriA,intSpin)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 40.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x spin/same=T: aa'
      else
        write(*,1010) 'IC/space x spin/same=T: aa',got,expected ; nfail=nfail+1
      endIf
!
!     Test 131: space ERI × spin integral, same=F ind=[1,3]
!     Code: aa=8*2=16, bb=8*3=24  (spin output)
!
      intOut = Contraction(eriA,intSpin,[1_int64,3_int64])
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 16.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x spin/same=F [1,3]: aa'
      else
        write(*,1010) 'IC/space x spin/same=F [1,3]: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 24.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x spin/same=F [1,3]: bb'
      else
        write(*,1010) 'IC/space x spin/same=F [1,3]: bb',got,expected ; nfail=nfail+1
      endIf
!
!     Test 132: space ERI × general integral, same=T
!     Code: aa = R4(aaaa,M_aa)+R4(aaaa,M_bb) = 8*2+8*3=40  (space output)
!
      intOut = Contraction(eriA,intGen)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 40.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x general/same=T: aa'
      else
        write(*,1010) 'IC/space x general/same=T: aa',got,expected ; nfail=nfail+1
      endIf
!
!     Test 133: space ERI × general integral, same=F ind=[1,3]
!     Code: aa=16, bb=24, ab=8*4=32, ba=8*5=40
!
      intOut = Contraction(eriA,intGen,[1_int64,3_int64])
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 16.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x general/same=F [1,3]: aa'
      else
        write(*,1010) 'IC/space x general/same=F [1,3]: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 24.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x general/same=F [1,3]: bb'
      else
        write(*,1010) 'IC/space x general/same=F [1,3]: bb',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ab') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 32.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x general/same=F [1,3]: ab'
      else
        write(*,1010) 'IC/space x general/same=F [1,3]: ab',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ba') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 40.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/space x general/same=F [1,3]: ba'
      else
        write(*,1010) 'IC/space x general/same=F [1,3]: ba',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 134-139: spin ERI × {space, spin, general}
!--------------------------------------------------------------------
!
!     Test 134: spin ERI (symm, aabb=bbaa=10) × space, same=T ind=[1,2] (ia=1 branch)
!     Code: aa = 8*2+10*2=36, bb = 9*2+10*2=38
!
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tC,tD,tE)
      intOut = Contraction(eriB,intSp)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 36.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin-symm x space/same=T [1,2]: aa'
      else
        write(*,1010) 'IC/spin-symm x space/same=T [1,2]: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 38.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin-symm x space/same=T [1,2]: bb'
      else
        write(*,1010) 'IC/spin-symm x space/same=T [1,2]: bb',got,expected ; nfail=nfail+1
      endIf
!
!     Test 135: spin ERI (full, aabb=10, bbaa=11) × space, same=T ind=[3,4] (ia=3 else branch)
!     Code: aa = 8*2+aabb*2=36, bb = 9*2+bbaa*2=9*2+11*2=40
!
      call mqc_twoeris_allocate(eriC,'full-doublebar','spin',tC,tD,tE,tF)
      intOut = Contraction(eriC,intSp,[3_int64,4_int64])
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 36.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin-full x space/same=T [3,4]: aa'
      else
        write(*,1010) 'IC/spin-full x space/same=T [3,4]: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 40.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin-full x space/same=T [3,4]: bb'
      else
        write(*,1010) 'IC/spin-full x space/same=T [3,4]: bb',got,expected ; nfail=nfail+1
      endIf
!
!     Test 136: spin ERI (symm, aabb=bbaa=10) × spin, same=T
!     Code: aa = 8*2+10*3=46, bb = 9*3+10*2=47
!
      intOut = Contraction(eriB,intSpin)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 46.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x spin/same=T: aa'
      else
        write(*,1010) 'IC/spin x spin/same=T: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 47.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x spin/same=T: bb'
      else
        write(*,1010) 'IC/spin x spin/same=T: bb',got,expected ; nfail=nfail+1
      endIf
!
!     Test 137: spin ERI × general integral, same=T
!     Code: same path as spin×spin, ab/ba blocks of intGen not used.
!     aa=46, bb=47 regardless of intGen%ab/ba values.
!
      intOut = Contraction(eriB,intGen)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 46.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x general/same=T: aa (ab/ba of intGen unused)'
      else
        write(*,1010) 'IC/spin x general/same=T: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 47.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x general/same=T: bb (ab/ba of intGen unused)'
      else
        write(*,1010) 'IC/spin x general/same=T: bb',got,expected ; nfail=nfail+1
      endIf
!
!     Test 138: spin ERI × general integral, same=F ind=[1,3] (ia<ib)
!     Bug-fix test: ba must use bbaa block, not aaaa.
!     Code: aa=8*2=16, bb=9*3=27, ab=aabb*M_ab=10*4=40, ba=bbaa*M_ba=10*5=50
!
      intOut = Contraction(eriB,intGen,[1_int64,3_int64])
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 16.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x general/same=F [1,3]: aa'
      else
        write(*,1010) 'IC/spin x general/same=F [1,3]: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 27.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x general/same=F [1,3]: bb'
      else
        write(*,1010) 'IC/spin x general/same=F [1,3]: bb',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ab') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 40.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x general/same=F [1,3]: ab'
      else
        write(*,1010) 'IC/spin x general/same=F [1,3]: ab',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ba') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 50.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x general/same=F [1,3]: ba (bbaa not aaaa)'
      else
        write(*,1010) 'IC/spin x general/same=F [1,3]: ba',got,expected ; nfail=nfail+1
      endIf
!
!     Test 139: spin ERI × general integral, same=F ind=[3,1] (ia>ib)
!     Bug-fix test: with swapped ind, ab/ba swap their integral block usage.
!     Code: aa=16, bb=27, ab=aabb*M_ba=10*5=50, ba=bbaa*M_ab=10*4=40
!
      intOut = Contraction(eriB,intGen,[3_int64,1_int64])
      outBlk = intOut%getBlock('ab') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 50.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x general/same=F [3,1]: ab (aabb*M_ba)'
      else
        write(*,1010) 'IC/spin x general/same=F [3,1]: ab',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ba') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 40.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/spin x general/same=F [3,1]: ba (bbaa*M_ab)'
      else
        write(*,1010) 'IC/spin x general/same=F [3,1]: ba',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 140-145: general ERI × {space, spin, general}
!
!     eriA = full-ms-general (no single-flip blocks):
!            aaaa=8,bbbb=9,aabb=bbaa=10,abba=baab=11,abab=baba=12
!     eriB = same + aaab=bbba=10 (for single-flip code paths)
!--------------------------------------------------------------------
!
      call mqc_twoeris_allocate(eriA,'full-ms','general',tC,tD,tE,tE, &
        abba=tF,baab=tF,abab=tG,baba=tG)
      call mqc_twoeris_allocate(eriB,'full-ms','general',tC,tD,tE,tE, &
        abba=tF,baab=tF,abab=tG,baba=tG,aaab=tE,bbba=tE)
!
!     Test 140: general ERI (with aaab=bbba=10) × space, same=T ind=[1,2] (ia=1)
!     Code: aa=8*2+10*2=36, bb=9*2+10*2=38
!           ab=R4(aaab,M_aa)+R4(bbab,M_aa)=10*2+0=20
!           ba=R4(aaba,M_aa)+R4(bbba,M_aa)=0+10*2=20
!
      intOut = Contraction(eriB,intSp)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 36.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen+SF x space/same=T [1,2]: aa'
      else
        write(*,1010) 'IC/gen+SF x space/same=T [1,2]: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 38.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen+SF x space/same=T [1,2]: bb'
      else
        write(*,1010) 'IC/gen+SF x space/same=T [1,2]: bb',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ab') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 20.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen+SF x space/same=T [1,2]: ab (aaab contribution)'
      else
        write(*,1010) 'IC/gen+SF x space/same=T [1,2]: ab',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ba') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 20.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen+SF x space/same=T [1,2]: ba (bbba contribution)'
      else
        write(*,1010) 'IC/gen+SF x space/same=T [1,2]: ba',got,expected ; nfail=nfail+1
      endIf
!
!     Test 141: general ERI (no single-flip) × spin, same=T ind=[1,2] (ia=1)
!     Code: aa=8*2+10*3=46, bb=9*3+10*2=47, ab=0, ba=0
!
      intOut = Contraction(eriA,intSpin)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 46.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x spin/same=T: aa'
      else
        write(*,1010) 'IC/gen x spin/same=T: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 47.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x spin/same=T: bb'
      else
        write(*,1010) 'IC/gen x spin/same=T: bb',got,expected ; nfail=nfail+1
      endIf
!
!     Test 142: general ERI (with bbba=10) × spin, same=F ind=[3,1] (ia=3, ib=1)
!     Bug-fix test: ba=R4(bbba,M_bb)=10*3=30 (old code: put into ab by mistake).
!     Code: aa=8*2+12*3=52, bb=9*3+12*2=51, ab=R4(aaab,M_aa)=10*2=20, ba=30
!
      intOut = Contraction(eriB,intSpin,[3_int64,1_int64])
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 52.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen+SF x spin/same=F [3,1]: aa'
      else
        write(*,1010) 'IC/gen+SF x spin/same=F [3,1]: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 51.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen+SF x spin/same=F [3,1]: bb'
      else
        write(*,1010) 'IC/gen+SF x spin/same=F [3,1]: bb',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ab') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 20.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen+SF x spin/same=F [3,1]: ab (aaab*M_aa, not bbba)'
      else
        write(*,1010) 'IC/gen+SF x spin/same=F [3,1]: ab',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ba') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 30.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen+SF x spin/same=F [3,1]: ba (bbba*M_bb)'
      else
        write(*,1010) 'IC/gen+SF x spin/same=F [3,1]: ba',got,expected ; nfail=nfail+1
      endIf
!
!     Test 143: general ERI × general integral, same=T ind=[1,2] (ia=1)
!     Code: aa=8*2+10*3=46, bb=9*3+10*2=47
!           ab=12*4+11*5=48+55=103, ba=12*5+11*4=60+44=104
!
      intOut = Contraction(eriA,intGen)
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 46.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=T: aa'
      else
        write(*,1010) 'IC/gen x general/same=T: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 47.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=T: bb'
      else
        write(*,1010) 'IC/gen x general/same=T: bb',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ab') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 103.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=T: ab'
      else
        write(*,1010) 'IC/gen x general/same=T: ab',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ba') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 104.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=T: ba'
      else
        write(*,1010) 'IC/gen x general/same=T: ba',got,expected ; nfail=nfail+1
      endIf
!
!     Test 144: general ERI × general integral, same=F ind=[1,3] (ia=1, ib=3)
!     Code: aa=8*2+12*3=52, bb=9*3+12*2=51
!           ab=10*4+11*5=40+55=95, ba=10*5+11*4=50+44=94
!
      intOut = Contraction(eriA,intGen,[1_int64,3_int64])
      outBlk = intOut%getBlock('aa') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 52.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=F [1,3]: aa'
      else
        write(*,1010) 'IC/gen x general/same=F [1,3]: aa',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('bb') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 51.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=F [1,3]: bb'
      else
        write(*,1010) 'IC/gen x general/same=F [1,3]: bb',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ab') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 95.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=F [1,3]: ab'
      else
        write(*,1010) 'IC/gen x general/same=F [1,3]: ab',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ba') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 94.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=F [1,3]: ba'
      else
        write(*,1010) 'IC/gen x general/same=F [1,3]: ba',got,expected ; nfail=nfail+1
      endIf
!
!     Test 145: general ERI × general integral, same=F ind=[3,1] (ia=3, ib=1)
!     Code: aa=52, bb=51, ab=10*5+11*4=50+44=94, ba=10*4+11*5=40+55=95
!
      intOut = Contraction(eriA,intGen,[3_int64,1_int64])
      outBlk = intOut%getBlock('ab') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 94.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=F [3,1]: ab'
      else
        write(*,1010) 'IC/gen x general/same=F [3,1]: ab',got,expected ; nfail=nfail+1
      endIf
      outBlk = intOut%getBlock('ba') ; sc = outBlk%at(1,1) ; got = sc%rval() ; expected = 95.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'IC/gen x general/same=F [3,1]: ba'
      else
        write(*,1010) 'IC/gen x general/same=F [3,1]: ba',got,expected ; nfail=nfail+1
      endIf
!
!
!--------------------------------------------------------------------
!     Tests 146-152: mqc_eri_r4tensor_contraction for space/spin/general
!
!     tC=8, tD=9, tE=10, tF=11, tG=12 at (1,1,1,1); n=2 basis functions.
!     tSO is a 4x4x4x4 spin-orbital tensor with:
!       1 at (1,1,1,1)=aaaa, (3,3,3,3)=bbbb,
!       (1,1,3,3)=aabb, (3,3,1,1)=bbaa
!--------------------------------------------------------------------
!
      call mqc_allocate_r4tensor(2*n,2*n,2*n,2*n,tSO,'real','full')
      call tSO%put(v1, 1,1,1,1)
      call tSO%put(v1, 3,3,3,3)
      call tSO%put(v1, 1,1,3,3)
      call tSO%put(v1, 3,3,1,1)
!
!     Test 146: space ERI (1x tensor, coulomb label)
!     aaaa alpha=8, tG(1,1,1,1)=12 → 8*12 = 96
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tC)
      sc = Contraction(eriA,tG,'coulomb')
      got = sc%rval() ; expected = 96.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/space/1x/coulomb'
      else
        write(*,1010) 'R4T/space/1x/coulomb',got,expected ; nfail=nfail+1
      endIf
!
!     Test 147: space ERI (1x tensor, exchange label)
!     K(1,1,1,1) = J(1,1,1,1) = 8 → 8*12 = 96
!
      sc = Contraction(eriA,tG,'exchange')
      got = sc%rval() ; expected = 96.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/space/1x/exchange'
      else
        write(*,1010) 'R4T/space/1x/exchange',got,expected ; nfail=nfail+1
      endIf
!
!     Test 148: space ERI (1x tensor, doublebar label)
!     (J-K)(1,1,1,1) = 8-8 = 0 → 0*12 = 0
!
      sc = Contraction(eriA,tG,'doublebar')
      got = sc%rval() ; expected = 0.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/space/1x/doublebar'
      else
        write(*,1010) 'R4T/space/1x/doublebar',got,expected ; nfail=nfail+1
      endIf
!
!     Test 149: space ERI (2x tensor, coulomb label)
!     spatial block (8) contracts with all 4 spin quadrants:
!     aaaa(1,1,1,1)=8*tSO(1,1,1,1)=8, bbbb(1,1,1,1)=8*tSO(3,3,3,3)=8,
!     aabb(1,1,1,1)=8*tSO(1,1,3,3)=8, bbaa(1,1,1,1)=8*tSO(3,3,1,1)=8 → 32
!
      sc = Contraction(eriA,tSO,'coulomb')
      got = sc%rval() ; expected = 32.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/space/2x/coulomb'
      else
        write(*,1010) 'R4T/space/2x/coulomb',got,expected ; nfail=nfail+1
      endIf
!
!     Test 150: symm-coulomb spin ERI (2x tensor, coulomb label, intType=label)
!     aaaa=8*tSO(1,1,1,1)=8, bbbb=9*1=9, aabb=10*1=10, bbaa=10*1=10 (symm)
!     expected = 37
!
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tC,tD,tE)
      sc = Contraction(eriB,tSO,'coulomb')
      got = sc%rval() ; expected = 37.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/spin/2x/coulomb (intType=label)'
      else
        write(*,1010) 'R4T/spin/2x/coulomb (intType=label)',got,expected ; nfail=nfail+1
      endIf
!
!     Test 151: full-doublebar spin ERI (2x tensor, doublebar label, intType=label)
!     aaaa=8*1=8, bbbb=9*1=9, aabb=10*1=10, bbaa=11*1=11 (full: stored separately)
!     expected = 38
!
      call mqc_twoeris_allocate(eriA,'full-doublebar','spin',tC,tD,tE,tF)
      sc = Contraction(eriA,tSO,'doublebar')
      got = sc%rval() ; expected = 38.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/spin/2x/doublebar (intType=label)'
      else
        write(*,1010) 'R4T/spin/2x/doublebar (intType=label)',got,expected ; nfail=nfail+1
      endIf
!
!     Test 152: full-coulomb general ERI (2x tensor, coulomb label, intType=label)
!     aaaa=8*1=8, bbbb=9*1=9, aabb=10*1=10, bbaa=10*1=10; no other blocks at tSO positions
!     expected = 37
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','general',tC,tD,tE,tE)
      sc = Contraction(eriA,tSO,'coulomb')
      got = sc%rval() ; expected = 37.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/general/2x/coulomb (intType=label)'
      else
        write(*,1010) 'R4T/general/2x/coulomb (intType=label)',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Tests 153-158: inverted and ms labels
!
!     Extend tSO with v2=2 at abab/baba/abba/baab positions:
!       tSO(1,3,1,3)=abab, tSO(3,1,3,1)=baba
!       tSO(1,3,3,1)=abba, tSO(3,1,1,3)=baab
!--------------------------------------------------------------------
!
      call tSO%put(v2, 1,3,1,3)
      call tSO%put(v2, 3,1,3,1)
      call tSO%put(v2, 1,3,3,1)
      call tSO%put(v2, 3,1,1,3)
!
!     Test 153: space ERI (1x tensor, inverted label)
!     K_inv(1,1,1,1)=J(1,1,1,1)=8 → 8*12 = 96
!
      call mqc_twoeris_allocate(eriA,'full-coulomb','space',tC)
      sc = Contraction(eriA,tG,'inverted')
      got = sc%rval() ; expected = 96.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/space/1x/inverted'
      else
        write(*,1010) 'R4T/space/1x/inverted',got,expected ; nfail=nfail+1
      endIf
!
!     Test 154: space ERI (1x tensor, ms label)
!     ms(1,1,1,1)=J(1,1,1,1)=8 → 8*12 = 96
!
      sc = Contraction(eriA,tG,'ms')
      got = sc%rval() ; expected = 96.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/space/1x/ms'
      else
        write(*,1010) 'R4T/space/1x/ms',got,expected ; nfail=nfail+1
      endIf
!
!     Test 155: space ERI (2x tensor, inverted label)
!     Inverted uses aaaa,bbbb,abab,baba (not aabb/bbaa):
!     aaaa: 8*tSO(1,1,1,1)=8, bbbb: 8*1=8,
!     abab: 8*tSO(1,3,1,3)=16, baba: 8*tSO(3,1,3,1)=16 → 48
!
      sc = Contraction(eriA,tSO,'inverted')
      got = sc%rval() ; expected = 48.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/space/2x/inverted'
      else
        write(*,1010) 'R4T/space/2x/inverted',got,expected ; nfail=nfail+1
      endIf
!
!     Test 156: space ERI (2x tensor, ms label)
!     All 8 quadrants: aaaa=8, bbbb=8, aabb=8, bbaa=8,
!     abba=8*2=16, baab=8*2=16, abab=16, baba=16 → 96
!
      sc = Contraction(eriA,tSO,'ms')
      got = sc%rval() ; expected = 96.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/space/2x/ms'
      else
        write(*,1010) 'R4T/space/2x/ms',got,expected ; nfail=nfail+1
      endIf
!
!     Test 157: symm-coulomb spin ERI, inverted label (calls mqc_eri_interaction)
!     K_inv: aaaa=8*1=8, bbbb=9*1=9, abab from aabb=10*2=20, baba=10*2=20
!     aabb/bbaa/abba/baab blocks of K_inv are zero → expected = 57
!
      call mqc_twoeris_allocate(eriB,'symm-coulomb','spin',tC,tD,tE)
      sc = Contraction(eriB,tSO,'inverted')
      got = sc%rval() ; expected = 57.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/spin/2x/coulomb→inverted'
      else
        write(*,1010) 'R4T/spin/2x/coulomb→inverted',got,expected ; nfail=nfail+1
      endIf
!
!     Test 158: symm-coulomb spin ERI, ms label (calls mqc_eri_interaction)
!     ms: aaaa=8*1=8, bbbb=9*1=9, aabb=10*1=10, bbaa=10*1=10,
!     abba(from aabb exchange)=10*2=20, baab=10*2=20,
!     abab(from aabb inverted)=10*2=20, baba=10*2=20 → 117
!
      sc = Contraction(eriB,tSO,'ms')
      got = sc%rval() ; expected = 117.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'R4T/spin/2x/coulomb→ms'
      else
        write(*,1010) 'R4T/spin/2x/coulomb→ms',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 159: complex-general AO->MO transformation.
!
!     The reference is an explicit spin-orbital four-index sum.  It
!     deliberately does not use MQC_TwoERIs%at, TwoERI_Trans, or the
!     partial-contraction spin dispatch.  Each of the 16 output spin
!     blocks is checked independently.
!--------------------------------------------------------------------
!
      call mqc_allocate_r4tensor(n,n,n,n,tA,'real','full')
      do mu=1,n
        do nu=1,n
          do lambda=1,n
            do sigma=1,n
              call tA%put(real((mu+nu)*(lambda+sigma)+mu*nu+lambda*sigma,&
                                  real64)/37.0_real64,&
                          mu,nu,lambda,sigma)
            end do
          end do
        end do
      end do
      call mqc_twoeris_allocate(eriAO,'full-coulomb','regular',tA)

      call moAA%init(n,n,(0.0_real64,0.0_real64))
      call moBB%init(n,n,(0.0_real64,0.0_real64))
      call moAB%init(n,n,(0.0_real64,0.0_real64))
      call moBA%init(n,n,(0.0_real64,0.0_real64))
      call moAA%put(( 0.81_real64, 0.07_real64),1,1)
      call moAA%put((-0.22_real64, 0.13_real64),1,2)
      call moAA%put(( 0.31_real64,-0.09_real64),2,1)
      call moAA%put(( 0.74_real64, 0.04_real64),2,2)
      call moBB%put(( 0.69_real64,-0.05_real64),1,1)
      call moBB%put(( 0.18_real64, 0.11_real64),1,2)
      call moBB%put((-0.27_real64, 0.08_real64),2,1)
      call moBB%put(( 0.77_real64,-0.06_real64),2,2)
      call moAB%put(( 0.17_real64, 0.14_real64),1,1)
      call moAB%put((-0.09_real64, 0.06_real64),1,2)
      call moAB%put(( 0.12_real64,-0.04_real64),2,1)
      call moAB%put(( 0.08_real64, 0.10_real64),2,2)
      call moBA%put((-0.11_real64, 0.09_real64),1,1)
      call moBA%put(( 0.07_real64,-0.12_real64),1,2)
      call moBA%put(( 0.05_real64, 0.03_real64),2,1)
      call moBA%put((-0.13_real64, 0.07_real64),2,2)
      call mqc_integral_allocate(generalMO,'complex general test MOs','general',&
                                 moAA,moBB,moAB,moBA)
      moFull=generalMO%getBlock('full')

      eriX=generalMO%outer2ERI(dagger(generalMO))
      eriHalf=mqc_eri_eri_partial_contraction_2ERI(eriAO,eriX,&
                                                   [4,3,1,2],[1,4,3,2])
      eriPartial=mqc_eri_eri_partial_contraction_2ERI(eriHalf,eriX,&
                                                      [1,2,3,4],[4,1,3,2])
      call TwoERI_Trans(6,0,generalMO,eriAO,eriTrans)

      spinLabels=[character(len=4)::'aaaa','aaab','aaba','aabb',&
        'abaa','abab','abba','abbb','baaa','baab','baba','babb',&
        'bbaa','bbab','bbba','bbbb']
      interactionLabels=[character(len=9)::'coulomb','exchange','doublebar']
      do iblock=1,size(spinLabels)
        spinLabel=spinLabels(iblock)
        tB=eriPartial%getBlock(spinLabel)
        tC=eriTrans%getBlock(spinLabel)
        partialError=0.0_real64
        transError=0.0_real64
        do p=1,n
          do q=1,n
            do r=1,n
              do s=1,n
                directValue=(0.0_real64,0.0_real64)
                colP=p+merge(0,n,spinLabel(1:1).eq.'a')
                colQ=q+merge(0,n,spinLabel(2:2).eq.'a')
                colR=r+merge(0,n,spinLabel(3:3).eq.'a')
                colS=s+merge(0,n,spinLabel(4:4).eq.'a')
                do spin1=0,1
                  do spin2=0,1
                    do mu=1,n
                      do nu=1,n
                        do lambda=1,n
                          do sigma=1,n
                            sc=moFull%at(mu+spin1*n,colP); coef1=sc%cval()
                            sc=moFull%at(nu+spin1*n,colQ); coef2=sc%cval()
                            sc=moFull%at(lambda+spin2*n,colR); coef3=sc%cval()
                            sc=moFull%at(sigma+spin2*n,colS); coef4=sc%cval()
                            sc=tA%at(mu,nu,lambda,sigma)
                            directValue=directValue+conjg(coef1)*coef2*&
                              conjg(coef3)*coef4*sc%cval()
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
                sc=tB%at(p,q,r,s); partialValue=sc%cval()
                sc=tC%at(p,q,r,s); transValue=sc%cval()
                partialError=max(partialError,abs(partialValue-directValue))
                transError=max(transError,abs(transValue-directValue))
              end do
            end do
          end do
        end do
        if(partialError.le.tol.and.transError.le.tol) then
          write(*,'(1X,A,1X,A4)') 'PASS  complex-general transform block',spinLabel
        else
          write(*,'(1X,A,1X,A4,2(A,ES14.6))') 'FAIL  complex-general transform block',&
            spinLabel,' partial=',partialError,' TwoERI_Trans=',transError
          nfail=nfail+1
        end if
      end do

!--------------------------------------------------------------------
!     Test 159b: different complex-general bra and ket transformations.
!
!     Use a genuinely different complex ket transformation so that this
!     distinguishes a bra/ket transformation from two independent
!     conventional transformations.
!--------------------------------------------------------------------
!
      moAA2=transpose(moAA)
      moBB2=transpose(moBB)
      moAB2=transpose(moAB)
      moBA2=transpose(moBA)
      call mqc_integral_allocate(generalMO2,'complex general ket MOs','general',&
                                 moAA2,moBB2,moAB2,moBA2)
      moFull2=generalMO2%getBlock('full')
      call TwoERI_Trans(6,0,generalMO,eriAO,eriTransLR,Right_MOs=generalMO2)
      tD=eriTransLR%getBlock('full')
      fullError=0.0_real64
      do p=1,2*n; do q=1,2*n; do r=1,2*n; do s=1,2*n
        directValue=(0.0_real64,0.0_real64)
        do spin1=0,1; do spin2=0,1
          do mu=1,n; do nu=1,n; do lambda=1,n; do sigma=1,n
            sc=moFull%at(mu+spin1*n,p); coef1=sc%cval()
            sc=moFull2%at(nu+spin1*n,q); coef2=sc%cval()
            sc=moFull%at(lambda+spin2*n,r); coef3=sc%cval()
            sc=moFull2%at(sigma+spin2*n,s); coef4=sc%cval()
            sc=tA%at(mu,nu,lambda,sigma)
            directValue=directValue+conjg(coef1)*coef2*&
              conjg(coef3)*coef4*sc%cval()
          end do; end do; end do; end do
        end do; end do
        sc=tD%at(p,q,r,s)
        fullError=max(fullError,abs(sc%cval()-directValue))
      end do; end do; end do; end do
      if(fullError.le.tol) then
        write(*,1000) 'different-bra/ket complex-general transformation'
      else
        write(*,'(1X,A,ES14.6)') 'FAIL  different-bra/ket transform error=',fullError
        nfail=nfail+1
      end if

!
!--------------------------------------------------------------------
!     Test 160: exhaustive complex-general MQC_TwoERIs%at oracle.
!
!     Construct full and symmetric storage from the already validated
!     transformed blocks.  For every element of every requested spin
!     block, compare Coulomb, exchange, and double-bar access against
!     the explicit block relationship
!       J(p,q,r,s)  = V[pqrs]
!       K(p,q,r,s)  = V[psrq]
!       JK(p,q,r,s) = V[pqrs] - V[psrq].
!--------------------------------------------------------------------
!
      do iblock=1,size(spinLabels)
        oracleBlocks(iblock)=eriTrans%getBlock(spinLabels(iblock))
      end do
      call mqc_twoeris_allocate(eriFull,'full-coulomb','general',&
        alpha=oracleBlocks(1),beta=oracleBlocks(16),&
        alphaBeta=oracleBlocks(4),betaAlpha=oracleBlocks(13),&
        abab=oracleBlocks(6),abba=oracleBlocks(7),baab=oracleBlocks(10),&
        baba=oracleBlocks(11),aaab=oracleBlocks(2),aaba=oracleBlocks(3),&
        abaa=oracleBlocks(5),baaa=oracleBlocks(9),bbba=oracleBlocks(15),&
        bbab=oracleBlocks(14),babb=oracleBlocks(12),abbb=oracleBlocks(8))
      call mqc_twoeris_allocate(eriSymm,'symm-coulomb','general',&
        alpha=oracleBlocks(1),beta=oracleBlocks(16),&
        alphaBeta=oracleBlocks(4),abab=oracleBlocks(6),&
        abba=oracleBlocks(7),aaab=oracleBlocks(2),abbb=oracleBlocks(8))

      do iinteraction=1,size(interactionLabels)
        do iblock=1,size(spinLabels)
          spinLabel=spinLabels(iblock)
          exchangeLabel=spinLabel(1:1)//spinLabel(4:4)//&
                        spinLabel(3:3)//spinLabel(2:2)
          jblock=0
          do p=1,size(spinLabels)
            if(spinLabels(p).eq.exchangeLabel) jblock=p
          end do
          if(jblock.eq.0) stop 'exchange spin block not found'
          do spin1=1,2
            atError=0.0_real64
            do p=1,n; do q=1,n; do r=1,n; do s=1,n
              sc=oracleBlocks(iblock)%at(p,q,r,s)
              directValue=sc%cval()
              sc=oracleBlocks(jblock)%at(p,s,r,q)
              if(interactionLabels(iinteraction).eq.'exchange') then
                directValue=sc%cval()
              elseif(interactionLabels(iinteraction).eq.'doublebar') then
                directValue=directValue-sc%cval()
              endIf
              if(spin1.eq.1) then
                sc=eriFull%at(p,q,r,s,spinLabel,interactionLabels(iinteraction))
              else
                sc=eriSymm%at(p,q,r,s,spinLabel,interactionLabels(iinteraction))
              endIf
              atError=max(atError,abs(sc%cval()-directValue))
            end do; end do; end do; end do
            if(atError.le.tol) then
              write(*,'(1X,A,1X,A,1X,A4,1X,A)') 'PASS  at oracle',&
                merge('full','symm',spin1.eq.1),spinLabel,trim(interactionLabels(iinteraction))
            else
              write(*,'(1X,A,1X,A,1X,A4,1X,A,A,ES14.6)') 'FAIL  at oracle',&
                merge('full','symm',spin1.eq.1),spinLabel,trim(interactionLabels(iinteraction)),&
                ' error=',atError
              nfail=nfail+1
            endIf
          end do
        end do
      end do
!
!     Exercise the independent no-spinBlock path using global
!     spin-orbital indices against the validated full tensor.
      tD=eriFull%getBlock('full')
      do iinteraction=1,size(interactionLabels)
        do spin1=1,2
          atError=0.0_real64
          do p=1,2*n; do q=1,2*n; do r=1,2*n; do s=1,2*n
            sc=tD%at(p,q,r,s)
            directValue=sc%cval()
            sc=tD%at(p,s,r,q)
            if(interactionLabels(iinteraction).eq.'exchange') then
              directValue=sc%cval()
            elseif(interactionLabels(iinteraction).eq.'doublebar') then
              directValue=directValue-sc%cval()
            endIf
            if(spin1.eq.1) then
              sc=eriFull%at(p,q,r,s,interactionIn=interactionLabels(iinteraction))
            else
              sc=eriSymm%at(p,q,r,s,interactionIn=interactionLabels(iinteraction))
            endIf
            atError=max(atError,abs(sc%cval()-directValue))
          end do; end do; end do; end do
          if(atError.le.tol) then
            write(*,'(1X,A,1X,A,1X,A)') 'PASS  at global oracle',&
              merge('full','symm',spin1.eq.1),trim(interactionLabels(iinteraction))
          else
            write(*,'(1X,A,1X,A,1X,A,A,ES14.6)') 'FAIL  at global oracle',&
              merge('full','symm',spin1.eq.1),trim(interactionLabels(iinteraction)),&
              ' error=',atError
            nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 161: real-general AO->MO transformation.
!
!     Retain nonzero alpha-beta and beta-alpha orbital mixing while
!     removing every imaginary coefficient.  This separates general-
!     spin block dispatch from complex conjugation errors.
!--------------------------------------------------------------------
!
      do p=1,n
        do q=1,n
          sc=moAA%at(p,q); call moAA%put(cmplx(real(sc%cval()),0.0_real64,real64),p,q)
          sc=moBB%at(p,q); call moBB%put(cmplx(real(sc%cval()),0.0_real64,real64),p,q)
          sc=moAB%at(p,q); call moAB%put(cmplx(real(sc%cval()),0.0_real64,real64),p,q)
          sc=moBA%at(p,q); call moBA%put(cmplx(real(sc%cval()),0.0_real64,real64),p,q)
        end do
      end do
      call mqc_integral_allocate(generalMO,'real general test MOs','general',&
                                 moAA,moBB,moAB,moBA)
      moFull=generalMO%getBlock('full')

      eriX=generalMO%outer2ERI(dagger(generalMO))
      eriHalf=mqc_eri_eri_partial_contraction_2ERI(eriAO,eriX,&
                                                   [4,3,1,2],[1,4,3,2])
      eriPartial=mqc_eri_eri_partial_contraction_2ERI(eriHalf,eriX,&
                                                      [1,2,3,4],[4,1,3,2])
      call TwoERI_Trans(6,0,generalMO,eriAO,eriTrans)

      do iblock=1,size(spinLabels)
        spinLabel=spinLabels(iblock)
        tB=eriPartial%getBlock(spinLabel)
        tC=eriTrans%getBlock(spinLabel)
        partialError=0.0_real64
        transError=0.0_real64
        do p=1,n
          do q=1,n
            do r=1,n
              do s=1,n
                directValue=(0.0_real64,0.0_real64)
                colP=p+merge(0,n,spinLabel(1:1).eq.'a')
                colQ=q+merge(0,n,spinLabel(2:2).eq.'a')
                colR=r+merge(0,n,spinLabel(3:3).eq.'a')
                colS=s+merge(0,n,spinLabel(4:4).eq.'a')
                do spin1=0,1
                  do spin2=0,1
                    do mu=1,n
                      do nu=1,n
                        do lambda=1,n
                          do sigma=1,n
                            sc=moFull%at(mu+spin1*n,colP); coef1=sc%cval()
                            sc=moFull%at(nu+spin1*n,colQ); coef2=sc%cval()
                            sc=moFull%at(lambda+spin2*n,colR); coef3=sc%cval()
                            sc=moFull%at(sigma+spin2*n,colS); coef4=sc%cval()
                            sc=tA%at(mu,nu,lambda,sigma)
                            directValue=directValue+conjg(coef1)*coef2*&
                              conjg(coef3)*coef4*sc%cval()
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
                sc=tB%at(p,q,r,s); partialValue=sc%cval()
                sc=tC%at(p,q,r,s); transValue=sc%cval()
                partialError=max(partialError,abs(partialValue-directValue))
                transError=max(transError,abs(transValue-directValue))
              end do
            end do
          end do
        end do
        if(partialError.le.tol.and.transError.le.tol) then
          write(*,'(1X,A,1X,A4)') 'PASS  real-general transform block',spinLabel
        else
          write(*,'(1X,A,1X,A4,2(A,ES14.6))') 'FAIL  real-general transform block',&
            spinLabel,' partial=',partialError,' TwoERI_Trans=',transError
          nfail=nfail+1
        end if
      end do

      call report_partner('aabb','bbaa',[3,4,1,2],.false.)
      call report_partner('aaab','aaba',[2,1,4,3],.true.)
      call report_partner('aaab','abaa',[3,4,1,2],.false.)
      call report_partner('aaab','baaa',[4,3,2,1],.true.)
      call report_partner('abba','baab',[3,4,1,2],.false.)
      call report_partner('abab','baba',[2,1,4,3],.true.)
      call report_partner('abbb','bbba',[4,3,2,1],.true.)
      call report_partner('abbb','bbab',[3,4,1,2],.false.)
      call report_partner('abbb','babb',[2,1,4,3],.true.)

      tD=eriTrans%getBlock('full')
      fullError=0.0_real64
      do p=1,2*n; do q=1,2*n; do r=1,2*n; do s=1,2*n
        directValue=(0.0_real64,0.0_real64)
        do spin1=0,1; do spin2=0,1
          do mu=1,n; do nu=1,n; do lambda=1,n; do sigma=1,n
            sc=moFull%at(mu+spin1*n,p); coef1=sc%cval()
            sc=moFull%at(nu+spin1*n,q); coef2=sc%cval()
            sc=moFull%at(lambda+spin2*n,r); coef3=sc%cval()
            sc=moFull%at(sigma+spin2*n,s); coef4=sc%cval()
            sc=tA%at(mu,nu,lambda,sigma)
            directValue=directValue+conjg(coef1)*coef2*&
              conjg(coef3)*coef4*sc%cval()
          end do; end do; end do; end do
        end do; end do
        sc=tD%at(p,q,r,s)
        fullError=max(fullError,abs(sc%cval()-directValue))
      end do; end do; end do; end do
      write(*,'(1X,A,ES14.6)') 'real-general full tensor oracle error=',fullError
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1

      contains

      subroutine report_partner(storedLabel,partnerLabel,order,doConjg)
        character(len=4),intent(in)::storedLabel,partnerLabel
        integer,dimension(4),intent(in)::order
        logical,intent(in)::doConjg
        type(mqc_r4tensor)::storedBlock,partnerBlock
        type(mqc_scalar)::storedScalar,partnerScalar
        integer::i1,i2,i3,i4
        integer,dimension(4)::idx
        complex(kind=real64)::storedValue,partnerValue

        storedBlock=eriPartial%getBlock(storedLabel)
        partnerBlock=eriPartial%getBlock(partnerLabel)
        partialPairError=0.0_real64
        do i1=1,size(partnerBlock,1); do i2=1,size(partnerBlock,2)
          do i3=1,size(partnerBlock,3); do i4=1,size(partnerBlock,4)
            idx=[i1,i2,i3,i4]
            storedScalar=storedBlock%at(idx(order(1)),idx(order(2)),&
                                        idx(order(3)),idx(order(4)))
            partnerScalar=partnerBlock%at(i1,i2,i3,i4)
            storedValue=storedScalar%cval()
            if(doConjg) storedValue=conjg(storedValue)
            partnerValue=partnerScalar%cval()
            partialPairError=max(partialPairError,abs(storedValue-partnerValue))
          end do; end do
        end do; end do

        storedBlock=eriTrans%getBlock(storedLabel)
        partnerBlock=eriTrans%getBlock(partnerLabel)
        transPairError=0.0_real64
        do i1=1,size(partnerBlock,1); do i2=1,size(partnerBlock,2)
          do i3=1,size(partnerBlock,3); do i4=1,size(partnerBlock,4)
            idx=[i1,i2,i3,i4]
            storedScalar=storedBlock%at(idx(order(1)),idx(order(2)),&
                                        idx(order(3)),idx(order(4)))
            partnerScalar=partnerBlock%at(i1,i2,i3,i4)
            storedValue=storedScalar%cval()
            if(doConjg) storedValue=conjg(storedValue)
            partnerValue=partnerScalar%cval()
            transPairError=max(transPairError,abs(storedValue-partnerValue))
          end do; end do
        end do; end do
        write(*,'(1X,A,1X,A4,A,A4,2(A,ES14.6))') 'real-general symmetry',&
          storedLabel,' -> ',partnerLabel,' partial=',partialPairError,&
          ' TwoERI_Trans=',transPairError
      end subroutine report_partner
!
      end program unitTest_twoERIs
