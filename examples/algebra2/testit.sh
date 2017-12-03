#! /bin/bash

./algebra_fun &> outfile
diff -b -B outfile OUTPUT/out


