%chk=rohf-hf.chk
%rwf=rohf-hf.rwf
#P rohf scf=conven output=mat int=noraff

test

1 2
h 
f 1 1.5

rohf-hf.mat

