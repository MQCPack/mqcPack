#! /bin/bash

cd tools
rm outfile
cd test
echo "test-Intentional problems"&>> ../outfile
../CheckGauInput test.com >> ../outfile
../CheckGauInput test.com Gaussian_not_in_path >> ../outfile

if command -v g16 > /dev/null 2>&1; then
 echo "Found g16" >> ../outfile
else
 echo "g16 not found" >> ../outfile
 diff -b -B outfile OUTPUT/out
 exit
fi

echo "test-Should be no problems">> ../outfile
../CheckGauInput test.com g16 >> ../outfile
cd ..

cd test1071
echo "test1071-Should be no problems" >> ../outfile
../CheckGauInput test1071.com g16 >> ../outfile
cd ..

cd test1129
echo "test1129-Should be no problems" >> ../outfile
../CheckGauInput test1129.com g16 >> ../outfile
cd ..

cd test1130
echo "test1130-Should be no problems" >> ../outfile
../CheckGauInput test1130.com g16 >> ../outfile
cd ..

cd test1132
echo "test1132-Should be no problems" >> ../outfile
../CheckGauInput test1132.com g16 >> ../outfile
cd ..

 diff -b -B outfile OUTPUT/out

exit
# tests 1135 and 1136 require R.  So a pain to run.
cd test1135
echo "test1135" >> ../outfile
../CheckGauInput test1135.com g16 >> ../outfile
cd ..

cd test1136
echo "test1136" >> ../outfile
../CheckGauInput test1136.com g16 >> ../outfile
cd ..


