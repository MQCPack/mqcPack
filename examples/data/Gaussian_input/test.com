#p ROHF scf=conven output=mAt int=noraff nosymm geom=nocrowd

test

0 1
h
h 1 5.0

rhf_h2-sto3g.mat

