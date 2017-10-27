%chk=y.chk
%rwf=y.rwf
#P uhf guess=mix scf=conven output=mat

test

1 4
h 
f 1 2.5

y.mat

