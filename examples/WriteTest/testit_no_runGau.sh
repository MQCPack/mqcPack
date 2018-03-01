#! /bin/bash -x

cd ../WriteTest
rm outfile

mkdir test2
cd test2
echo "------------------------" >> ../outfile 2>&1
echo "Pass through MatFile test">> ../outfile
echo "------------------------" >> ../outfile 2>&1
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../../PrintAllData/DataPrint rhf_h2-sto3g.mat > Original
sed -e 'sZ-0.000000Z 0.000000Zg' < Original > ../Original
../ReadWrite rhf_h2-sto3g.mat
../../PrintAllData/DataPrint rhf_h2-sto3g.mat.copy > Copy
sed -e 'sZ-0.000000Z 0.000000Zg' < Copy > Copy2
sed -e 'sZ.copyZZg' < Copy2 > ../Copy
cd ..
rm -r test2

diff -b -B Original Copy

exit
