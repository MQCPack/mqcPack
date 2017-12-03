#p uhf/3-21g geom=modela test out=matrix symm=scf guess=alter pop=full

Gaussian Test Job 1132 (Part 1):
hof excited triplet saving unformatted matrix elements

0,3
o f h


7,10

test1132u.dat

--Link1--
#p uhf/3-21g geom=modela test pop=full

Gaussian Test Job 1132 (Part 2):
hof, lowest triplet

0,3
o f h

--Link1--
#p uhf/3-21g geom=modela test pop=full external("testunf test1132u.dat",postscf,iounf)

Gaussian Test Job 1132 (Part 3):
hof, lowest triplet but read in orbitals, etc. from unf from
higher triplet

0,3
o f h

--Link1--
#p uhf/3-21g geom=modela test out=raw symm=scf guess=alter pop=full

Gaussian Test Job 1132 (Part 4):
hof excited triplet saving raw matrix elements

0,3
o f h


7,10

test1132r.dat

--Link1--
#p uhf/3-21g geom=modela test pop=full external("testunf test1132r.dat",postscf,iounf,raw)

Gaussian Test Job 1132 (Part 5):
hof, lowest triplet but read in orbitals, etc. from raw matrix elements from
higher triplet

0,3
o f h

--Link1--
#P TEST GHF/6-31G* 5d Pop=Full SCF=NoVarAcc out=raw test iop(4/13=5)

Gaussian Test Job 1132 (Part 6):
O2 GHF starting as singlet, going down to triplet, saving matrix elements

 0 1
 O 
 O 1 R

R 1.220

test1132gr.dat

--Link1--
%chk=test1132
#p rhf/6-31g* 5d pop=full test

Gaussian Test Job 1132 (Part 7):
O2 singlet saving chk

0,1
 O 
 O 1 R

R 1.220

--Link1--
%chk=test1132
#p rhf/6-31g* 5d pop=full external("testunf test1132gr.dat",postscf,iounf,raw) test

Gaussian Test Job 1132 (Part 8):
O2 singlet reading ghf triplet raw file

0,1
 O 
 O 1 R

R 1.220

--Link1--
%chk=test1132
%nosave
#p ghf/6-31g* 5d pop=full geom=check guess=read test

Gaussian Test Job 1132 (Part 9):
O2 ghf reading from updated chk file

0,1
