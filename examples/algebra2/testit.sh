#! /bin/bash

cd algebra2
./algebra_fun &> outfile
diff -b -B outfile OUTPUT/out


