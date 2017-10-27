%chk=sxh.chk
%rwf=sxh.rwf
#P uhf guess=mix scf=conven output=mat

test

0 1
h 
h 1 0.7

sxh.mat

