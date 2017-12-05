#p external=("readmat.pl",iofchk,inunf,post,mo2el,raw) geom=modela test scf=(conven,qc) iop(3/11=3)

Gaussian Test Job 1130 (Part 1):
hof, test of raw matrix element file via perl rdmat run from external

1,2
o f h

--Link1--
#p geom=modela test scf=(conven,qc) iop(3/11=3) output=(raw,mo2ele)

Gaussian Test Job 1130 (Part 2):
hof, test of raw matrix element file save and perl rdmat run separately

1,2
o f h

test1130.dat

