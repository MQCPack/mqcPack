#! /bin/bash -x

cd ../WriteTest
if [ -f Copy ]; then
    rm Copy
fi
if [ -f Original ]; then
    rm Original
fi

# set up the Gaussian environment
export GauBINARY="${GAU_BINARY}"

mkdir test2
cd test2
echo "------------------------"
echo "Pass through MatFile test"
echo "------------------------" 
cp ../../data/MatrixFile/rhf_h2-sto3g.mat .
../../PrintAllData/DataPrint rhf_h2-sto3g.mat > Original
sed -e 'sZ-0.000000Z 0.000000Zg' < Original > ../Original
../ReadWrite rhf_h2-sto3g.mat
../../PrintAllData/DataPrint rhf_h2-sto3g.mat.copy > Copy
sed -e 'sZ-0.000000Z 0.000000Zg' < Copy > Copy2
sed -e 'sZ.copyZZg' < Copy2 > ../Copy
cd ..
rm -r test2

mkdir test1071
cd test1071
echo "-------------------------------"
echo "test1071 Create MatrixFile test"
echo "-------------------------------"
!
! If mat file from DataSummary run is around, use it
!
if [ -d ../../DataSummary/test1071 ]; then
    cp ../../DataSummary/test1071/test1071.mat .
    ../../PrintAllData/DataPrint test1071.mat > Original
else
    cp ../../data/Gaussian_input/test1071.com .
    ../../PrintAllData/DataPrint test1071.com ${GauBINARY} > Original
    mkdir ../../DataSummary/test1071
    cp test1071.mat ../../DataSummary/test1071/test1071.mat
    rm -f test1071.com
fi
sed -e 'sZ-0.000000Z 0.000000Zg' < Original >> ../Original
../ReadWrite test1071.mat 
../../PrintAllData/DataPrint test1071.mat.copy > Copy
sed -e 'sZ-0.000000Z 0.000000Zg' < Copy > Copy2
sed -e 'sZ.copyZZg' < Copy2 >> ../Copy
cd ..
rm -r test1071

mkdir test1129
cd test1129
echo "-------------------------------"
echo "test1129 Create MatrixFile test"
echo "-------------------------------"
!
! If mat file from DataSummary run is around, use it
!
if [ -d ../../DataSummary/test1129 ]; then
    cp ../../DataSummary/test1129/test1129.dat .
    ../../PrintAllData/DataPrint test1129.dat > Original
else
    cp ../../data/Gaussian_input/test1129.com .
    ../../PrintAllData/DataPrint test1129.com ${GauBINARY} > Original
    mkdir ../../DataSummary/test1129
    cp test1129.dat ../../DataSummary/test1129/test1129.dat
    rm -f test1129.com
fi
sed -e 'sZ-0.000000Z 0.000000Zg' < Original >> ../Original
../ReadWrite test1129.dat 
../../PrintAllData/DataPrint test1129.dat.copy > Copy
sed -e 'sZ-0.000000Z 0.000000Zg' < Copy > Copy2
sed -e 'sZ.copyZZg' < Copy2 >> ../Copy
cd ..
rm -r test1129

mkdir test1132
cd test1132
echo "--------------------------------"
echo "test1132 Create MatrixFiles test"
echo "--------------------------------"
!
! If mat file from DataSummary run is around, use it
!
if [ -d ../../DataSummary/test1132 ]; then
    cp ../../DataSummary/test1132/test1132.mat .
    cp ../../DataSummary/test1132/test1132r.dat .
    cp ../../DataSummary/test1132/test1132gr.dat .
    cp ../../DataSummary/test1132/test1132u.dat .
    ../../PrintAllData/DataPrint test1132.mat > Original
else
    cp ../../data/Gaussian_input/test1132.com .
    ../../PrintAllData/DataPrint test1132.com ${GauBINARY} > Original
    mkdir ../../DataSummary/test1132
    cp test1132.mat ../../DataSummary/test1132/test1132.mat
    cp test1132r.dat ../../DataSummary/test1132/test1132r.dat
    cp test1132gr.dat ../../DataSummary/test1132/test1132gr.dat
    cp test1132u.dat ../../DataSummary/test1132/test1132u.dat
    rm -f test1132.com
fi
../../PrintAllData/DataPrint test1132u.dat >> Original
sed -e 'sZ-0.000000Z 0.000000Zg' < Original >> ../Original
../ReadWrite test1132.mat 
../../PrintAllData/DataPrint test1132.mat.copy > Copy
../ReadWrite test1132u.dat
../../PrintAllData/DataPrint test1132u.dat.copy >> Copy
sed -e 'sZ-0.000000Z 0.000000Zg' < Copy > Copy2
sed -e 'sZ.copyZZg' < Copy2 >> ../Copy
cd ..
rm -r test1132

diff -b -B Copy Original

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


