#!/bin/bash
# ==========================================================================
# gammalib CentOS package creation
#
# Copyright (C) 2017 Sylvie Brau-Nogué
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------
#
# This script installs the ctools and GammaLib packages in the location
# /usr/local/gamma and then builds a CentOS rpm package.
# Any existing content in /usr/local/gamma will be
# destroyed. The script needs the priviledges to create a folder in
# /usr/local.
#
# ==========================================================================

# ====================================== #
# Get GammaLib version from command line #
# ====================================== #
if [ $# -ne 1 ] ;  then
    echo "Please specify GammaLib in the form 'x.y.z'"
    exit 1
fi
VERSION=$1
PLATFORM=`uname -m`


# ====== #
# Header #
# ====== #
echo "================================================================="
echo "Build $PACKNAME rpm package, Version $VERSION, Platform $PLATFORM"
echo "================================================================="


# =============================== #
# Set software component versions #
# =============================== #
PACKNAME=gammalib
PACKAGE=$PACKNAME-$VERSION
RELEASE=1

# ============== #
# Set parameters #
# ============== #
WRKDIR=$PWD/rpm/rpmbuild
SRCDIR=$WRKDIR/SOURCES
PKGDIR=$WRKDIR/RPMS
SPKGDIR=$WRKDIR/SRPMS
PRODDIR=$WRKDIR/RPMS/$PLATFORM


# ======================= #
# Clean package directory #
# ======================= #
rm -rf $WRKDIR

# ============================= #
# Create package directory tree #
# ============================= #
mkdir -p $WRKDIR/{RPMS,SRPMS,BUILD,SOURCES,SPECS,tmp}


# ======================================== #
# Copy package configuration file (*.spec) #
# ======================================== #
cp $PACKNAME.spec $WRKDIR/SPECS/


# ======================================= #
# Create tarball and move it into $SRCDIR #
# ======================================= #
make dist
mv gammalib-*.tar.gz $SRCDIR/


# ================= #
# Post-process code #
# ================= #


# ==================== #
# Build CentOS package #
# ==================== #
rpmbuild --target $PLATFORM -ba --define "_topdir $WRKDIR" $PACKNAME.spec -v


# ============================== #
# Copy package to root directory #
# ============================== #
echo "===================================================="
echo "Packaging DONE"
echo " "
echo "New binary packages are located in $PRODDIR"
for package in $(find $PRODDIR -iname '*.rpm' ); do 
     echo "rpm file : $package"
     echo "===================================================="
     rpm -qip $package
     cp $package $PWD/
done
echo "===================================================="
echo " "
echo "New source packages are located in $SPKGDIR"
for package in $(find $SPKGDIR -iname '*.rpm' ); do
     echo "rpm file : $package"
     echo "===================================================="
     rpm -qip $package
     cp $package $PWD/
done
