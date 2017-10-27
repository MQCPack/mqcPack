%chk=x.chk
%rwf=x.rwf
#P uhf guess=read scf=conven output=mat int=noraff

test

1 2
h 
f 1 2.5

x.mat

