#! /bin/bash -x

cd ../CreateDataList
rm outfile
mkdir test
cd test
# set up the Gaussian environment
echo "---------------------------------" >> ../outfile 2>&1
echo "Failure Test: no input to program" >> ../outfile 2>&1
echo "---------------------------------" >> ../outfile 2>&1
../CheckInput >> ../outfile
echo "-------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Only 1 input to program" >> ../outfile 2>&1
echo "-------------------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test.com . >> ../outfile
../CheckInput test.com >> ../outfile
echo "------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Program is not in path" >> ../outfile 2>&1
echo "------------------------------------" >> ../outfile 2>&1
../CheckInput test.com not_in_path >> ../outfile
cd ..
rm -r test

mkdir test2
cd test2
rm -f test.com
echo "---------------------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Test for when input file is not there" >> ../outfile 2>&1
echo "---------------------------------------------------" >> ../outfile 2>&1
../CheckInput test.com not_in_path >> ../outfile
echo "Finished with tests designed to fail">> ../outfile
echo "-------------" >> ../outfile 2>&1
echo "MatFile test">> ../outfile
echo "-------------" >> ../outfile 2>&1
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../CheckInput rhf_h2-sto3g.mat not_in_path >> ../outfile

rm -f rhf_h2-sto3g.mat
rm -f test.com

cd ..
rm -r test2
grep -v "echo argv" outfile > tmpfile 2>&1
grep -v "cp test1132" tmpfile > outfile 2>&1
rm tmpfile

diff -b -B outfile OUTPUT/out_no_runGau

exit
