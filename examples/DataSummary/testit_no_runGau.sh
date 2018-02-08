#! /bin/bash -x

cd ../DataSummary
rm outfile
mkdir test
cd test
# set up the Gaussian environment
echo "-------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Only 1 input to program" >> ../outfile 2>&1
echo "-------------------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test.com . >> ../outfile
../DataSummary test.com >> ../outfile
echo "------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Program is not in path" >> ../outfile 2>&1
echo "------------------------------------" >> ../outfile 2>&1
../DataSummary test.com not_in_path >> ../outfile
cd ..
rm -r test

mkdir test2
cd test2
rm -f test.com
echo "---------------------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Test for when input file is not there" >> ../outfile 2>&1
echo "---------------------------------------------------" >> ../outfile 2>&1
../DataSummary test.com not_in_path >> ../outfile
echo "Finished with tests designed to fail">> ../outfile
echo "-------------" >> ../outfile 2>&1
echo "MatFile test">> ../outfile
echo "-------------" >> ../outfile 2>&1
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../DataSummary rhf_h2-sto3g.mat >> ../outfile

cd ..
rm -r test2

diff -b -B outfile OUTPUT/out_no_runGau

exit
