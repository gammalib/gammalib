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
if [ $# -ne 2 ] ;  then
    echo "Please specify GammaLib in the form 'x.y.z x86_64/i386'"
    exit 1
fi
VERSION=$1
PLATFORM=$2

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
PRODDIR=$WRKDIR/RPMS/$PLATFORM
#RPMFILE=$PRODDIR/$PACKAGE-$RELEASE-*$PLATFORM.rpm
#LOGFILE=rpm/pkg_build/pkg_build.log
#LOGDEPILE=rpm/pkg_build/pkg_dependencies_rpm.log


# ==================================================== #
#  TODO : Create a secure RPM repository
# ==================================================== #
# Temporary solution, to be fixed
YUMREPO=$PWD/rpm/YUM-REPO

# ======================= #
# Clean package directory #
# ======================= #
rm -rf $WRKDIR

# ========================= #
# Clean temporary directory #
# ========================= #
#rm -rf $TMPDIR

# ============================= #
# Create package directory tree #
# ============================= #
mkdir -p $WRKDIR/{RPMS,SRPMS,BUILD,SOURCES,SPECS,tmp}

# =============== #
# Install ncurses #
# =============== #
# TODO : eventuellement fixer le numero de version du paquet
# rpm -qa | grep ncurses || sudo yum install ncurses
#
#if ! rpm -qa | grep -q ncurses; then
#    echo "ncurses need to be installed"
#    sudo yum install ncurses
#    if ! rpm -qa | grep -q ncurses; then
#        echo "================================================================================================="
#        echo "ERROR : ncurses need to be installed manually"
#        echo "   "
#        echo "================================================================================================="
#        exit
#    fi
#fi
# Check ncurses version
#v=$(rpm --qf '%{VERSION}\n' -q ncurses)
#echo "ncurses $v checked"

# ================ #
# Install readline #
# ================ #
# TODO : eventuellement fixer le numero de version du paquet
# rpm -qa | grep readline || sudo yum install readline

#if ! rpm -qa | grep -q readline; then
#    echo "readline need to be installed"
#    sudo yum install readline
#    if ! rpm -qa | grep -q readline; then
#        echo "================================================================================================="
#        echo "ERROR : readline need to be installed manually"
#        echo "   "
#        echo "================================================================================================="
#        exit
#    fi
#fi
# Check readline version
#v=$(rpm --qf '%{VERSION}\n' -q readline)
#echo "readline $v checked"

# =============== #
# Install cfitsio #
# =============== #
# TODO : eventuellement fixer le numero de version du paquet
# rpm -qa | grep cfitsio || sudo yum install cfitsio

#if ! rpm -qa | grep -q cfitsio; then
#    echo "cfitsio need to be installed"
#    sudo yum install cfitsio-devel
#    if ! rpm -qa | grep -q cfitsio; then
#        echo "================================================================================================="
#        echo " ERROR : cfitsio need to be installed"
#        echo "================================================================================================="
#        echo "You can download the latest version of CFITSIO on http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html"
#        echo "   "
#        echo "or "
#        echo "   "
#        echo "cfitsio rpm can be found in the epel-release repository ; Try : sudo yum install cfitsio-devel"
#        echo "   "
#        echo "================================================================================================="
#        exit
#    fi
#fi
# Check cfitsio version
#v=$(rpm --qf '%{VERSION}\n' -q cfitsio)
#echo "cfitsio $v checked"

# ======================================== #
# Copy package configuration file (*.spec) #
# ======================================== #
cp $PACKNAME.spec $WRKDIR/SPECS/


# ===================================== #
# Create tarball and move it in $SRCDIR #
# ===================================== #
make dist
mv gammalib-*.tar.gz $SRCDIR/
#cp $PACKAGE.tar.gz $SRCDIR/


# ================= #
# Post-process code #
# ================= #


# ==================== #
# Build CentOS package #
# ==================== #
#rpmbuild --target $PLATFORM -ba $WRKDIR/SPECS/$PACKNAME.spec -v
#rpmbuild --target $PLATFORM -ba --buildroot $WRKDIR $PACKNAME.spec -v
rpmbuild --target $PLATFORM -ba --define "_topdir $WRKDIR" $PACKNAME.spec -v

# ==================================================== #
# Sign package with an environment variable is set by the pass phrase during login
# for example in .bashrc file
# > GPG_PASSPHRASE="tombelaneige"
# > export GPG_PASSPHRASE
# Another possibility to create a specific file in .gnupg directory
# > at login script, extract pass phrase from this file
# ==================================================== #
(\
	echo set timeout -1;\
	echo spawn rpmsign --addsign $PKGDIR/*/*.rpm;\
	echo expect -re \"pass\";\
	echo send -- \"$GPG_PASSPHRASE\\r\";\
	echo expect eof;\
) | expect


echo " ===================================================="
echo " Packaging DONE :"
echo " New package is located at $PRODDIR"
for package in $(find $PRODDIR -iname *.rpm ); do 
     echo " rpm file : $package"
     echo " ===================================================="
     echo "  "     
     rpm -qip $package; 

	# ==================================================== #
	#  TODO : Move package in a secure RPM repository
	# ==================================================== #
	# Temporary solution, to be fixed
     cp $package $YUMREPO
done
exit
