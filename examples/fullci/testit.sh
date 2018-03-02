#! /bin/bash -x

cd ../fullci
rm outfile
# set up the Gaussian environment
mkdir workdir
cd workdir >> outfile 2>&1
echo "---------------------------------" >> ../outfile 2>&1
echo "Failure Test: no input to program" >> ../outfile 2>&1
echo "---------------------------------" >> ../outfile 2>&1
./fullci >> ../outfile
echo "Finished with tests designed to fail">> ../outfile
echo "-------------" >> ../outfile 2>&1
echo "MatFile test">> ../outfile
echo "-------------" >> ../outfile 2>&1
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../fullci rhf_h2-sto3g.mat >> ../outfile
cd ..
rm -r workdir

sed -e 'sZ-0.000000Z 0.000000Zg' < outfile > outfile_tmp
mv outfile_tmp outfile

diff -b -B outfile OUTPUT/out

exit
