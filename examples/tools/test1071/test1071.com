%chk=test1071
#p blyp/3-21g/auto geom(modela) opt(calcfc) freq scrf(saveq) test 

Gaussian Test Job 1071 (Part 1):
PCM pre-opt+freq saving the charges on the checkpoint file.

0 1
c o h h

--Link1--
%chk=test1071
#p b3lyp/6-311G(d,p) geom(check) guess(check) scrf(readq,writeq,read) test 

Gaussian Test Job 1071 (Part 2):
PCM energy with larger basis iterative reading and saving the PCM charges.

0 1

iterative qconv=tight

--Link1--
%chk=test1071
#p b3lyp/chkbasis geom(check) guess(check) opt(readfc) freq scrf(check,readq,writeq) test 

Gaussian Test Job 1071 (Part 3):
PCM opt+freq with larger basis iterative reading the PCM charges.

0 1

--Link1--
%chk=test1071
#p b3lyp/chkbasis geom(check) guess(check) scrf(external,1stvac,readq,writeq) test 

Gaussian Test Job 1071 (Part 4):
PCM external iterations reading and writing the PCM charges.

0 1

--Link1--
%chk=test1071
%nosave
#p b3lyp/chkbasis geom(check) guess(check) scrf(external,readq) test 

Gaussian Test Job 1071 (Part 5):
PCM external iterations reading the PCM charges.

0 1

