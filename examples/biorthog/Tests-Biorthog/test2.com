%chk=test2.chk
%rwf=test2.rwf
#P uhf guess=mix scf=conven pop=(biorthog,savebio) IOp(6/33=3) 

test

1 2
h 
f 1 2.5


