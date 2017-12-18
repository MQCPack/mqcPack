#! /bin/bash -x

cd ../tools
rm outfile
cd test
# set up the Gaussian environment

echo "test-Intentional problems" >> ../outfile 2>&1
../CheckGauInput test.com >> ../outfile
../CheckGauInput test.com Gaussian_not_in_path >> ../outfile

export GauBINARY="${GAU_BINARY}"
cd ../test2
echo "test-Should be no problems">> ../outfile
../CheckGauInput test.com ${GauBINARY} >> ../outfile 2>&1
cd ..

cd test1071
echo "test1071-Should be no problems" >> ../outfile 2>&1
../CheckGauInput test1071.com ${GauBINARY} >> ../outfile 2>&1
../CheckGauInput test.com ${GauBINARY} >> ../outfile 2>&1
cd ..

cd test1129
echo "test1129-Should be no problems" >> ../outfile 2>&1
../CheckGauInput test1129.com ${GauBINARY} >> ../outfile 2>&1
../CheckGauInput test.com ${GauBINARY} >> ../outfile 2>&1
cd ..

cd test1132
echo "test1132-Run Input should be no problems" >> ../outfile 2>&1
../CheckGauInput test1132.com ${GauBINARY} >> ../outfile 2>&1
echo "test1132-MatrixFile" >> ../outfile 2>&1
../CheckGauInput test1132r.dat ${GauBINARY} >> ../outfile 2>&1
echo "test1132-Run Input and MatrixFile" >> ../outfile 2>&1
../CheckGauInput test1132.com test1132gr.dat test1132r.dat test1132u.dat ${GauBINARY} >> ../outfile 2>&1
../CheckGauInput test.com ${GauBINARY} >> ../outfile 2>&1
cd ..

# Now, remove the process IDs from the names of Gaussian file names
cp outfile outfile_keep
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
../CheckGauInput test1130.com ${GauBINARY} >> ../outfile 2>&1
cd ..

cd test1135
echo "test1135" >> ../outfile
../CheckGauInput test1135.com g16 >> ../outfile
cd ..

cd test1136
echo "test1136" >> ../outfile
../CheckGauInput test1136.com g16 >> ../outfile
cd ..


