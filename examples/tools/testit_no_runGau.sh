#! /bin/bash -x

cd ../tools
rm outfile
cd test
# set up the Gaussian environment

echo "Intentional problems" >> ../outfile 2>&1
../CheckGauInput test.com >> ../outfile
../CheckGauInput test.com Gaussian_not_in_path >> ../outfile
cd ..

diff -b -B outfile OUTPUT/out_no_runGau

exit
