#! /bin/bash -x

cd ../DataSummary
rm outfile
mkdir test
cd test
echo "-------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Only 1 input to program" >> ../outfile 2>&1
echo "-------------------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test.com . >> ../outfile
../DataSummary test.com >> ../outfile
echo "------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Program is not in path" >> ../outfile 2>&1
echo "------------------------------------" >> ../outfile 2>&1
../DataSummary test.com not_in_path >> ../outfile

# set up the Gaussian environment
export GauBINARY="${GAU_BINARY}"
cd ..
rm -r test

mkdir test2
cd test2
rm -f test.com
echo "---------------------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Test for when input file is not there" >> ../outfile 2>&1
echo "---------------------------------------------------" >> ../outfile 2>&1
../DataSummary test.com ${GauBINARY} >> ../outfile
echo "Finished with tests designed to fail">> ../outfile
echo "------------------------" >> ../outfile 2>&1
echo "Pass through MatFile test">> ../outfile
echo "------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test.com . >> ../outfile
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../DataSummary rhf_h2-sto3g.mat ${GauBINARY} >> ../outfile
echo "-------------------------------" >> ../outfile 2>&1
echo "test.com Create MatrixFile test">> ../outfile
echo "-------------------------------" >> ../outfile 2>&1
../DataSummary test.com ${GauBINARY} >> ../outfile
rm -f rhf_h2-sto3g.mat
rm -f test.com
cd ..
rm -r test2

mkdir test1071
cd test1071
echo "-------------------------------" >> ../outfile 2>&1
echo "test1071 Create MatrixFile test" >> ../outfile 2>&1
echo "-------------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test1071.com . >> ../outfile
../DataSummary test1071.com ${GauBINARY} >> ../outfile 2>&1
rm -f test1071.com
cd ..
rm -r  test1071

mkdir test1129
cd test1129
cp ../../data/Gaussian_input/test1129.com . >> ../outfile
echo "-------------------------------" >> ../outfile 2>&1
echo "test1129 Create MatrixFile test" >> ../outfile 2>&1
echo "-------------------------------" >> ../outfile 2>&1
../DataSummary test1129.com ${GauBINARY} >> ../outfile 2>&1
rm -f test1129.com
cd ..
rm -r test1129

mkdir test1132
cd test1132
cp ../../data/Gaussian_input/test1132.com . >> ../outfile
echo "--------------------------------" >> ../outfile 2>&1
echo "test1132 Create MatrixFiles test" >> ../outfile 2>&1
echo "--------------------------------" >> ../outfile 2>&1
../DataSummary test1132.com ${GauBINARY} >> ../outfile 2>&1
echo "-----------------------------------------" >> ../outfile 2>&1
echo "test1132 Pass one MatrixFile through test" >> ../outfile 2>&1
echo "-----------------------------------------" >> ../outfile 2>&1
../DataSummary test1132r.dat ${GauBINARY} >> ../outfile 2>&1
echo "------------------------------------------------------------------" >> ../outfile 2>&1
echo "test1132 Create One MartixFile and Pass through 3 MatrixFiles test" >> ../outfile 2>&1
echo "------------------------------------------------------------------" >> ../outfile 2>&1
../DataSummary test1132.com test1132gr.dat test1132r.dat test1132u.dat ${GauBINARY} >> ../outfile 2>&1
rm -f test1132.com
cd ..
rm -r test1132

grep -v "echo argv" outfile > tmpfile 2>&1
grep -v "cp test1132" tmpfile > outfile 2>&1
grep -v "Gaussian Version:" outfile > tmpfile 2>&1
mv tmpfile outfile

cd OUTPUT
rm -f out_runGau
var1=$(stat -c%s out_runGau_v2)
var2=$(stat -c%s ../outfile)
if [[ "$var1" -eq "$var2" ]]; then
    ln -s out_runGau_v2 out_runGau
else
    ln -s out_runGau_v1 out_runGau
fi
cd ..

diff -b -B outfile OUTPUT/out_runGau


exit
# tests 1130, 1135 and 1136 require R.  So a pain to run.
cd test1130
echo "test1130-Should be no problems" >> ../outfile 2>&1
../DataSummary test1130.com ${GauBINARY} >> ../outfile 2>&1
cd ..

cd test1135
echo "test1135" >> ../outfile
../DataSummary test1135.com g16 >> ../outfile
cd ..

cd test1136
echo "test1136" >> ../outfile
../DataSummary test1136.com g16 >> ../outfile
cd ..


