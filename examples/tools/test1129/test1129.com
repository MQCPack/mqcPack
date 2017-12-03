#p external=("rdmat-ext",iofchk,inunf,post,mo2el) geom=modela test scf=(conven,qc) iop(3/11=3)

Gaussian Test Job 1129 (Part 1):
hof, test of matrix element file via fortran rdmat run from external

1,2
o f h

--Link1--
#p geom=modela test scf=(conven,qc) iop(3/11=3) output=(matrix,mo2ele)

Gaussian Test Job 1129 (Part 2):
hof, test of matrix element file save and fortran rdmat run separately

1,2
o f h

test1129.dat

