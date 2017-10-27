%chk=z.chk
%rwf=z.rwf
#P uhf guess=mix scf=conven output=mat int=noraff

test

0 1
h 
h 1 2.5

z.mat

