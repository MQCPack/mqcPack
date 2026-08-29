#! /bin/bash -x

: "${abs_top_srcdir:=$(cd ../.. && pwd)}"
: "${abs_top_builddir:=$(cd ../.. && pwd)}"

cd "${abs_top_builddir}/examples/hartreefock"
rm -f outfile
# set up the Gaussian environment
rm -rf workdir
mkdir workdir
cd workdir >> outfile 2>&1
echo "---------------------------------" >> ../outfile 2>&1
echo "Failure Test: no input to program" >> ../outfile 2>&1
echo "---------------------------------" >> ../outfile 2>&1
./hartreefock >> ../outfile
echo "Finished with tests designed to fail">> ../outfile
echo "-------------" >> ../outfile 2>&1
echo "MatFile test">> ../outfile
echo "-------------" >> ../outfile 2>&1
cp "${abs_top_srcdir}/examples/data/MatrixFile/rhf_h2-sto3g.mat" .  \
  >> ../outfile
../hartreefock -f rhf_h2-sto3g.mat >> ../outfile
cd ..
rm -r workdir

sed -e 'sZ-0.000000Z 0.000000Zg' < outfile > outfile_tmp
mv outfile_tmp outfile

diff -b -B outfile "${abs_top_srcdir}/examples/hartreefock/OUTPUT/out"

exit
