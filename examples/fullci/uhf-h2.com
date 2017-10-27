#P uhf guess=mix scf=conven output=(mat,MO2ElectronIntegrals) int=noraff nosymm

test

0 1
h 
h 1 2.0

uhf-h2.mat

