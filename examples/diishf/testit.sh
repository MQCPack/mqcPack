#! /bin/bash -x

cd ../diishf
rm outfile
# set up the Gaussian environment
mkdir workdir
cd workdir >> outfile 2>&1
echo "---------------------------------" >> ../outfile 2>&1
echo "Failure Test: no input to program" >> ../outfile 2>&1
echo "---------------------------------" >> ../outfile 2>&1
./diishf >> ../outfile
echo "Finished with tests designed to fail">> ../outfile
echo "-------------" >> ../outfile 2>&1
echo "MatFile test">> ../outfile
echo "-------------" >> ../outfile 2>&1
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../diishf -f rhf_h2-sto3g.mat >> ../outfile
cd ..
rm -r workdir

sed -e 'sZ-0.000000Z 0.000000Zg' < outfile > outfile_tmp
mv outfile_tmp outfile

diff -b -B outfile OUTPUT/out

exit
