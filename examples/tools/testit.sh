#! /bin/bash -x

cd ../tools
rm outfile
cd test
# set up the Gaussian environment

echo "test-Intentional problems" >> ../outfile 2>&1
../CheckGauInput test.com >> ../outfile
../CheckGauInput test.com Gaussian_not_in_path >> ../outfile

export GauBINARY="${GAU_BINARY}"

echo "test-Should be no problems">> ../outfile
../CheckGauInput test.com ${GauBINARY} >> ../outfile
cd ..

cd test1071
pwd
echo "test1071-Should be no problems" >> ../outfile
../CheckGauInput test1071.com ${GauBINARY} >> ../outfile
cd ..

cd test1129
pwd
echo "test1129-Should be no problems" >> ../outfile
../CheckGauInput test1129.com ${GauBINARY} >> ../outfile
cd ..

cd test1130
pwd
echo "test1130-Should be no problems" >> ../outfile
../CheckGauInput test1130.com ${GauBINARY} >> ../outfile
cd ..

cd test1132
pwd
echo "test1132-Run Input should be no problems" >> ../outfile
../CheckGauInput test1132.com ${GauBINARY} >> ../outfile
echo "test1132-MatrixFile" >> ../outfile
../CheckGauInput test1132r.dat ${GauBINARY} >> ../outfile
echo "test1132-Run Input and MatrixFile" >> ../outfile
../CheckGauInput test1132.com test1132gr.dat test1132r.dat test1132u.dat ${GauBINARY} >> ../outfile
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


