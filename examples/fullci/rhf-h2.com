#P rohf scf=conven output=(mat,mo2elec) int=noraff nosymm

test

0 1
h 
h 1 2.0

rhf-h2.mat

