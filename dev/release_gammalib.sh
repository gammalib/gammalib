#! /bin/bash -f

# Set parameters
basedir=$(pwd)
installdir=$(pwd)/build
makedoc="no"
makelatex="no"

echo "Checkout latest GammaLib"
rm -rf gammalib
git clone -b release https://cta-gitlab.irap.omp.eu/gammalib/gammalib.git
rm -rf gammalib/.git

# Extract version number
version=`cat gammalib/configure.ac | grep 'AC_INIT' | awk -F"[" '{print $3}' | sed 's/],//' | sed 's/\./ /g'`
version=`printf "%d.%d.%d" $version`
echo "GammaLib version: "$version

# Set source directory
sourcedir=gammalib-$version
rm -rf $sourcedir
mv gammalib $sourcedir

# Create Python wrappers using SWIG. We need here a quite complicated kluge to get this done
# as the Makefile in pyext does not exist (and we don't want that it exists in the distro).
# We thus make a copy of gammalib, create the Makefile there, call 'make swig' in pyext,
# and copy over the wrappers into the gammalib distro
echo "Create Python wrappers using SWIG (making SWIG obsolete for installation)"
swig -version
tmpdir=$basedir/tmp_for_swig
rm -rf $tmpdir
cp -r $sourcedir $tmpdir
cd $tmpdir
./autogen.sh
./configure
cd pyext
make swig
cd $basedir
cp -r $tmpdir/pyext/gammalib $sourcedir/pyext
rm -rf $tmpdir

echo "Create Makefile.in and configure scripts"
cd $sourcedir
./autogen.sh
rm -f gammalib.sh
rm -rf dev
rm -rf autom4te.cache
cd $basedir

# Create Documentation
if [ "x$makedoc" == "xyes" ] ; then
  echo "Create Documentation"
  cd $sourcedir
  make doc
  cd $basedir
fi

# Create LaTeX documentation
if [ "x$makelatex" == "xyes" ] ; then
  echo "Create LaTeX documentation"
  cd $sourcedir/doc
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
  cd $basedir
fi

# Set permissions
echo "Make files read/write"
chmod -R u+rw $sourcedir

# Create tarball
echo "Create tarball"
tar cvf - $sourcedir > $sourcedir.tar
rm -f $sourcedir.tar.gz
gzip $sourcedir.tar
#exit 0

# Configure GammaLib
echo "Configure GammaLib"
cd $sourcedir
./configure --prefix=$installdir

# Build GammaLib
echo "Build GammaLib"
make -j10

# Check GammaLib
echo "Check GammaLib"
make check

# Install GammaLib
echo "Install GammaLib"
make install
