#!/bin/bash
# ==========================================================================
# gammalib Mandriva package creation
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
# /usr/local/gamma and then builds a Mandriva rpm package.
# Any existing content in /usr/local/gamma will be
# destroyed. The script needs the priviledges to create a folder in
# /usr/local.
#
# ==========================================================================

# ============================================= #
# Get ctools/GammaLib version from command line #
# ============================================= #
if [ $# -ne 4 ] ;  then
echo "Please specify package name version platform in the form 'package x.y.z platform '"
exit 1
fi
PACKNAME=$1
VERSION=$2
PLATFORM=$3
GVERSION=$4

echo " ================================================================= "
echo " Build $PACKNAME rpm package, Version $VERSION, Platform $PLATFORM "
echo " ================================================================= "

# =============================== #
# Set software component versions #
# =============================== #
#CFITSIO=cfitsio3410
#NCURSES=lib64ncurses-devel
#READLINE=readline-6.3
CFITSIO=lib64cfitsio-devel
NCURSES=lib64ncurses-devel
READLINE=lib64readline-devel
#
PACKAGE=$PACKNAME-$VERSION
RELEASE=1

# ============== #
# Set parameters #
# ============== #
#
# System Configuration Hypothesis with file ~/.rpmmacros.
# This file must contain at least the follwing lines :
# %packager KNODLSEDER
# %_topdir %(echo $HOME)/rpmbuild
# %_tmppath %(echo $HOME)/rpmbuild/tmp
# %_smp_mflags  -j3
# %__arch_install_post   /usr/lib/rpm/check-rpaths   /usr/lib/rpm/check-buildroot
#
#
# SBN : Liste à fixer.
TMPDIR=$HOME/usr
WRKDIR=$HOME/rpmbuild
SRCDIR=$WRKDIR/SOURCES
PKGDIR=$WRKDIR/RPMS
PRODDIR=$WRKDIR/RPMS/$PLATFORM
RPMFILE=$PRODDIR/$PACKAGE-$RELEASE-*$PLATFORM.rpm
#
LOGFILE=$HOME/pkg_build/pkg_build.log
LOGDEPILE=$HOME/pkg_build/pkg_dependencies_rpm.log

# ==================================================== #
#  TODO : Create a secure RPM repository
# ==================================================== #
# Temporary solution, to be fixed
YUMREPO=$HOME/YUM-REPO

# ======================= #
# Clean package directory #
# ======================= #
rm -rf $WRKDIR

# ========================= #
# Clean temporary directory #
# ========================= #
# TODO : revoir ce point précis
# <=== ligne suivante commentée pour des problèmes de permissions de création de répertoire ensuite
# rm -rf $TMPDIR

# ============================= #
# Create package directory tree #
# ============================= #
mkdir -p $HOME/rpmbuild/{RPMS,SRPMS,BUILD,SOURCES,SPECS,tmp}

# =============== #
# Install ncurses #
# =============== #
# TODO : eventuellement fixer le numero de version du paquet
# rpm -qa | grep ncurses || sudo yum install ncurses

if ! rpm -qa | grep -q $NCURSES; then
    echo "$NCURSES need to be installed"
    sudo urpmi $NCURSES
    if ! rpm -qa | grep $NCURSES; then
        echo "================================================================================================="
        echo "ERROR : $NCURSES need to be installed manually"
        echo "   "
        echo "================================================================================================="
        exit
    fi
fi
# Check ncurses version
v=$(rpm --qf '%{VERSION}\n' -q $NCURSES)
echo "$NCURSES $v checked"


# ================ #
# Install readline #
# ================ #
# TODO : eventuellement fixer le numero de version du paquet
# rpm -qa | grep readline || sudo yum install readline

if ! rpm -qa | grep -q $READLINE; then
    echo "$READLINE need to be installed"
    sudo urpmi $READLINE
    if ! rpm -qa | grep $READLINE; then
        echo "================================================================================================="
        echo "ERROR : $READLINE need to be installed manually"
        echo "   "
        echo "================================================================================================="
        exit
    fi
fi
# Check readline version
v=$(rpm --qf '%{VERSION}\n' -q $READLINE)
echo "$READLINE $v checked"

# =============== #
# Install cfitsio #
# =============== #
# TODO : eventuellement fixer le numero de version du paquet
# rpm -qa | grep cfitsio || sudo yum install cfitsio

if ! rpm -qa | grep -q $CFITSIO; then
    echo "$CFITSIO need to be installed"
    sudo urpmi $CFITSIO
    if ! rpm -qa | grep $CFITSIO; then
        echo "================================================================================================="
        echo " ERROR : $CFITSIO need to be installed"
        echo "================================================================================================="
        echo "You can download the latest version of CFITSIO on http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html"
        echo "   "
        echo "================================================================================================="
        exit
    fi
fi
# Check cfitsio version
v=$(rpm --qf '%{VERSION}\n' -q $CFITSIO)
echo "$CFITSIO $v checked"

# ======================================== #
# Copy package configuration file (*.spec) #
# ======================================== #
cp $PACKNAME.spec $WRKDIR/SPECS/

# ================================================================== #
# Just copy dist in SOURCES dir                                      #
# ================================================================== #
cp $PACKAGE.tar.gz $SRCDIR/


# ================= #
# Post-process code #
# ================= #


# ====================== #
# Build centOS package #
# ====================== #

rpmbuild --target $PLATFORM -ba $WRKDIR/SPECS/$PACKNAME.spec -v

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
