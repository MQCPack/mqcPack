#! /bin/bash -x

cd ../tools
rm outfile
cd test
echo "-------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Only 1 input to program" >> ../outfile 2>&1
echo "-------------------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test.com . >> ../outfile
../FullWavefunctionRW test.com >> ../outfile
echo "------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Program is not in path" >> ../outfile 2>&1
echo "------------------------------------" >> ../outfile 2>&1
../FullWavefunctionRW test.com not_in_path >> ../outfile

# set up the Gaussian environment
export GauBINARY="${GAU_BINARY}"
cd ../test2
rm -f test.com
echo "---------------------------------------------------" >> ../outfile 2>&1
echo "Failure Test: Test for when input file is not there" >> ../outfile 2>&1
echo "---------------------------------------------------" >> ../outfile 2>&1
../FullWavefunctionRW test.com ${GauBINARY} >> ../outfile
echo "Finished with tests designed to fail">> ../outfile
echo "------------------------" >> ../outfile 2>&1
echo "Pass through MatFile test">> ../outfile
echo "------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test.com . >> ../outfile
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../FullWavefunctionRW rhf_h2-sto3g.mat ${GauBINARY} >> ../outfile
echo "-------------------------------" >> ../outfile 2>&1
echo "test.com Create MatrixFile test">> ../outfile
echo "-------------------------------" >> ../outfile 2>&1
../FullWavefunctionRW test.com ${GauBINARY} >> ../outfile
rm -f rhf_h2-sto3g.mat
rm -f test.com
cd ..

cd test1071
echo "-------------------------------" >> ../outfile 2>&1
echo "test1071 Create MatrixFile test" >> ../outfile 2>&1
echo "-------------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test1071.com . >> ../outfile
../FullWavefunctionRW test1071.com ${GauBINARY} >> ../outfile 2>&1
rm -f test1071.com
cd ..

cd test1129
cp ../../data/Gaussian_input/test1129.com . >> ../outfile
echo "-------------------------------" >> ../outfile 2>&1
echo "test1129 Create MatrixFile test" >> ../outfile 2>&1
echo "-------------------------------" >> ../outfile 2>&1
../FullWavefunctionRW test1129.com ${GauBINARY} >> ../outfile 2>&1
rm -f test1129.com
cd ..

cd test1132
cp ../../data/Gaussian_input/test1132.com . >> ../outfile
echo "--------------------------------" >> ../outfile 2>&1
echo "test1132 Create MatrixFiles test" >> ../outfile 2>&1
echo "--------------------------------" >> ../outfile 2>&1
###../FullWavefunctionRW test1132.com ${GauBINARY} >> ../outfile 2>&1
echo "-----------------------------------------" >> ../outfile 2>&1
echo "test1132 Pass one MatrixFile through test" >> ../outfile 2>&1
echo "-----------------------------------------" >> ../outfile 2>&1
###../FullWavefunctionRW test1132r.dat ${GauBINARY} >> ../outfile 2>&1
echo "------------------------------------------------------------------" >> ../outfile 2>&1
echo "test1132 Create One MartixFile and Pass through 3 MatrixFiles test" >> ../outfile 2>&1
echo "------------------------------------------------------------------" >> ../outfile 2>&1
###../FullWavefunctionRW test1132.com test1132gr.dat test1132r.dat test1132u.dat ${GauBINARY} >> ../outfile 2>&1
rm -f test1132.com
cd ..

# Now, remove the process IDs from the names of Gaussian file names
../../.other_libs/build_fcns/remove_Gau_pid outfile
../../.other_libs/build_fcns/remove_Gau_pid outfile
../../.other_libs/build_fcns/remove_Gau_pid outfile
../../.other_libs/build_fcns/remove_Gau_pid outfile
../../.other_libs/build_fcns/remove_Gau_pid outfile
../../.other_libs/build_fcns/remove_Gau_pid outfile
# remove the location of the scratch directory from outfile
../../.other_libs/build_fcns/remove_SCRDIR outfile

 diff -b -B outfile OUTPUT/out_runGau

exit
# tests 1130, 1135 and 1136 require R.  So a pain to run.
cd test1130
echo "test1130-Should be no problems" >> ../outfile 2>&1
../FullWavefunctionRW test1130.com ${GauBINARY} >> ../outfile 2>&1
cd ..

cd test1135
echo "test1135" >> ../outfile
../FullWavefunctionRW test1135.com g16 >> ../outfile
cd ..

cd test1136
echo "test1136" >> ../outfile
../FullWavefunctionRW test1136.com g16 >> ../outfile
cd ..


