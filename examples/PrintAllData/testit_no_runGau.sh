#! /bin/bash -x

cd ../PrintAllData
rm outfile

mkdir test2
cd test2
rm -f test.com
echo "------------------------" >> ../outfile 2>&1
echo "Pass through MatFile test">> ../outfile
echo "------------------------" >> ../outfile 2>&1
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../DataPrint rhf_h2-sto3g.mat >> ../outfile
cd ..
rm -r test2

diff -b -B outfile OUTPUT/out_no_runGau

exit
