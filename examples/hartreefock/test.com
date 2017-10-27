#p rohf scf=conven output=mat int=noraff nosymm geom=nocrowd

test

0 1
h
h 1 5.0

rhf_h2-sto3g.mat

