#! /bin/bash

cd readTest
./algebra_fun1 &> outfile
diff -b -B outfile OUTPUT/out
