#! /bin/bash

chmod 755 ../*/testit.sh
cd algebra2
./algebra_fun &> outfile
diff -b -B outfile OUTPUT/out


