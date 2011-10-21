#!/bin/bash -f

# Set parameters
cfitsio=cfitsio3280
ncurses=ncurses-5.9
readline=readline-6.2
gammalib=gammalib-00-04-12
ctools=ctools-00-03-00
installdir=/usr/local/gamma
swig_path=/Users/jurgen/Software/swig/swig2.0.4/bin

# Set flags
do_delete=1
do_cfitsio=1
do_ncurses=1
do_readline=1
do_gammalib=1
do_ctools=1
do_python25=1
do_python26=1
do_python27=1

# Set configuration
export CFLAGS="-arch ppc -arch i386 -arch x86_64"
export CXXFLAGS="-arch ppc -arch i386 -arch x86_64"
export LDFLAGS="-arch ppc -arch i386 -arch x86_64"

# Remove old installation
if [ $do_delete -eq 1 ]; then
  sudo rm -rf $installdir
fi

# Dump header
rm install.log
echo > "Mac OS X packaging" > install.log

# Install cfitsio
if [ $do_cfitsio -eq 1 ]; then
  tar xfz $cfitsio.tar.gz
  cd cfitsio
  echo " " >> ../install.log
  echo "cfitsio" >> ../install.log
  ./configure --prefix=$installdir
  make shared
  sudo make install
  file $installdir/lib/libcfitsio.a >> ../install.log
  file $installdir/lib/libcfitsio.dylib >> ../install.log
  otool -L $installdir/lib/libcfitsio.dylib >> ../install.log
  cd ..
fi

# Install ncurses (not needed since available on Mac OS X)
if [ $do_ncurses -eq 1 ]; then
  tar xfz $ncurses.tar.gz
  cd $ncurses
  echo " " >> ../install.log
  echo "Ncurses" >> ../install.log
  ./configure --prefix=$installdir --with-shared
  make
  sudo make install
  file $installdir/lib/libncurses.dylib >> ../install.log
  otool -L $installdir/lib/libncurses.dylib >> ../install.log
  cd ..
fi

# Install readline (not needed since available on Mac OS X)
if [ $do_readline -eq 1 ]; then
  tar xfz $readline.tar.gz
  cd $readline
  echo " " >> ../install.log
  echo "Readline" >> ../install.log
  ./configure --prefix=$installdir
  make
  sudo make install
  file $installdir/lib/libreadline.dylib >> ../install.log
  otool -L $installdir/lib/libreadline.dylib >> ../install.log
  cd ..
fi

# Decide which Python versions to build
VERSIONS=
declare -a PYTHONS
let count=0
if [ $do_python25 -eq 1 ]; then
  VERSIONS=$VERSIONS"2.5 "
  PYTHONS[$count]="/System/Library/Frameworks/Python.framework/Versions/2.5/bin/python2.5"
  ((count++))
fi
if [ $do_python26 -eq 1 ]; then
  VERSIONS=$VERSIONS"2.6 "
  PYTHONS[$count]="/System/Library/Frameworks/Python.framework/Versions/2.6/bin/python2.6"
  ((count++))
fi
if [ $do_python27 -eq 1 ]; then
  VERSIONS=$VERSIONS"2.7 "
  PYTHONS[$count]="/Users/jurgen/Software/python/python27/10.5/3-way/bin/python2.7"
  ((count++))
fi  

# Install gammalib
if [ $do_gammalib -eq 1 ]; then

  # Secure PATH
  OLD_PATH=$PATH

  # Decompress tarball and step into GammaLib
  tar xfz $gammalib.tar.gz
  cd $gammalib
  
  # Make GammaLib builds
  let count=0
  for VERSION in $VERSIONS; do

    # Dump version
    echo " " >> ../install.log
    echo "GammaLib for Python $VERSION" >> install.log

    # Create symbolic link for this Python version
    rm -rf python$VERSION
    mkdir -p python$VERSION
    cd python$VERSION
    ln -s ${PYTHONS[$count]} python
    export PATH=$PWD:$swig_path:/usr/bin:/bin:/usr/sbin:/sbin
    echo $PATH >> ../../install.log
    cd ..
    
    # Dump Python version
    python -V >> ../install.log
    
    # Remove pyext/build
    rm -rf pyext/build
    
    # Build GammaLib
    ./configure --prefix=$installdir
    make -j8
    sudo make install
    file $installdir/lib/libgamma.la >> ../install.log
    file $installdir/lib/libgamma.dylib >> ../install.log
    file $installdir/lib/python$VERSION/site-packages/_gammalib.so >> ../install.log
    otool -L $installdir/lib/libgamma.dylib >> ../install.log
    otool -L $installdir/lib/python$VERSION/site-packages/_gammalib.so >> ../install.log

    # Increment version counter
    ((count++))

  done

  # GammaLib post install
  cd ..
  sudo cp gammalib-setup $installdir/bin/gammalib-setup

  # Recover PATH
  export PATH=$OLD_PATH

fi


# Install ctools
if [ $do_ctools -eq 1 ]; then

  # Secure PATH
  OLD_PATH=$PATH

  # Decompress tarball and step into ctools
  tar xfz $ctools.tar.gz
  cd $ctools
  
  # Make ctools builds
  let count=0
  for VERSION in $VERSIONS; do

    # Dump version
    echo " " >> ../install.log
    echo "ctools for Python $VERSION" >> install.log

    # Create symbolic link for this Python version
    rm -rf python$VERSION
    mkdir -p python$VERSION
    cd python$VERSION
    ln -s ${PYTHONS[$count]} python
    export PATH=$PWD:$swig_path:/usr/bin:/bin:/usr/sbin:/sbin
    echo $PATH >> ../../install.log
    cd ..
    
    # Dump Python version
    python -V >> ../install.log

    # Build ctools
    make distclean
    ./configure --prefix=$installdir
    make
    sudo make install
    file $installdir/lib/libctools.la >> ../install.log
    file $installdir/lib/libctools.dylib >> ../install.log
    file $installdir/lib/python$VERSION/site-packages/_ctools.so >> ../install.log
    otool -L $installdir/lib/libctools.dylib >> ../install.log
    otool -L $installdir/lib/python$VERSION/site-packages/_ctools.so >> ../install.log

    # Increment version counter
    ((count++))

  done

  # ctools post install
  cd ..

  # Recover PATH
  export PATH=$OLD_PATH

fi
