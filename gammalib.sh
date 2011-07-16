#! /bin/bash -f

# Set parameters
installdir=$(pwd)/build
makedoc="no"

echo "Checkout latest GammaLib"
rm -rf gammalib
cvs -Q -d /home/cvs export -N -r HEAD "gammalib"

# Extract version number
version=`cat gammalib/configure.ac | grep 'AC_INIT' | awk -F"[" '{print $3}' | sed 's/],//' | sed 's/\./ /g'`
version=`printf "%2.2d-%2.2d-%2.2d" $version`
echo "GammaLib version: "$version

# Set source directory
sourcedir=gammalib-$version
rm -rf $sourcedir
mv gammalib $sourcedir

#
echo "Create Makefile.in and configure scripts"
cd $sourcedir
./autogen.sh
rm -f gammalib.sh
cd ..

#
echo "Create gammalib_wrap.cpp and gammalib.py files for python binding (making swig obsolete)"
inc_inst="-I$sourcedir/inst/mwl/pyext -I$sourcedir/inst/lat/pyext -I$sourcedir/inst/cta/pyext"
opt_inst="-DWITH_INST_MWL -DWITH_INST_LAT -DWITH_INST_CTA"
swig -c++ -python -Wall -includeall -I$sourcedir/src $inc_inst $opt_inst -o $sourcedir/pyext/gammalib_wrap.cpp -outdir $sourcedir/pyext $sourcedir/pyext/gammalib.i

#
if [ "x$makedoc" == "xyes" ] ; then
  echo "Create doxygen documentation"
  cd $sourcedir/doc
  doxygen Doxyfile
  cd latex
  make
  make
  dvips -o refman.ps refman
  ps2pdf refman.ps 
  cd ../..
fi

#
echo "Create LaTeX documentation"
cd $sourcedir/doc
#
latex gammalib_um.tex
dvips gammalib_um.dvi -o gammalib_um.ps
ps2pdf gammalib_um.ps
rm -f gammalib_um.dvi gammalib_um.aux gammalib_um.log gammalib_um.tex gammalib_um.toc manual.sty gammalib_um.ps
#
cd dev/coding
latex gammalib_coding.tex
dvips gammalib_coding.dvi -o gammalib_coding.ps
ps2pdf gammalib_coding.ps
cd ../..
mv dev/coding/gammalib_coding.pdf .
#
cd dev/inst
latex gammalib_inst.tex
dvips gammalib_inst.dvi -o gammalib_inst.ps
ps2pdf gammalib_inst.ps
cd ../..
mv dev/inst/gammalib_inst.pdf .
#
cd dev/maths
latex gammalib_maths.tex
dvips gammalib_maths.dvi -o gammalib_maths.ps
ps2pdf gammalib_maths.ps
cd ../..
mv dev/maths/gammalib_maths.pdf .
#
cd dev/sdd 
latex gammalib_sdd.tex
dvips gammalib_sdd.dvi -o gammalib_sdd.ps
ps2pdf gammalib_sdd.ps
cd ../..
mv dev/sdd/gammalib_sdd.pdf .
#
cd dev/srs 
latex gammalib_srs.tex
dvips gammalib_srs.dvi -o gammalib_srs.ps
ps2pdf gammalib_srs.ps
cd ../..
mv dev/srs/gammalib_srs.pdf .
#
rm -rf dev
cd ../..

echo "Make files read/write"
chmod -R u+rw $sourcedir

echo "Create tarball"
tar cvf - $sourcedir > $sourcedir.tar
rm -f $sourcedir.tar.gz
gzip $sourcedir.tar

#
echo "Configure GammaLib"
cd $sourcedir
./configure --prefix=$installdir

#
echo "Compile GammaLib"
make -j10

#
echo "Check GammaLib"
make check

#
echo "Install GammaLib"
make install
