#! /bin/bash

cd ../../unitTests/src
./unitTest_twoERIs &> outfile
diff -b -B outfile OUTPUT/out
