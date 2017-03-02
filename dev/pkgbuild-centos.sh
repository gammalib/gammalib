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

# ============================================= #
# Get ctools/GammaLib version from command line #
# ============================================= #
if [ $# -ne 2 ] ;  then
echo "Please specify package name and version in the form 'package x.y.z'"
exit 1
fi
PACKNAME=$1
VERSION=$2


# =============================== #
# Set software component versions #
# =============================== #
#CFITSIO=cfitsio3410
#NCURSES=ncurses-5.9
#READLINE=readline-6.3
CFITSIO=cfitsio
NCURSES=ncurses
READLINE=readline
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
WRKDIR=$HOME/rpmbuild
SRCDIR=$WRKDIR/SOURCES
PKGDIR=$WRKDIR/RPMS
PRODDIR=$WRKDIR/RPMS/noarch
RPMFILE=$PRODDIR/$PACKAGE-$RELEASE.noarch.rpm
#
LOGFILE=$HOME/pkg_build/pkg_build.log
LOGDEPILE=$HOME/pkg_build/pkg_dependencies_rpm.log

# ======================= #
# Clean package directory #
# ======================= #
rm -rf $WRKDIR

# ============================= #
# Create package directory tree #
# ============================= #
mkdir -p $HOME/rpmbuild/{RPMS,SRPMS,BUILD,SOURCES,SPECS,tmp}

# =============== #
# Install ncurses #
# =============== #
# TODO : eventuellement fixer le numero de version du paquet
rpm -qa | grep ncurses || sudo yum install ncurses

# ================ #
# Install readline #
# ================ #
# TODO : eventuellement fixer le numero de version du paquet
rpm -qa | grep readline || sudo yum install readline

# =============== #
# Install cfitsio #
# =============== #
# TODO : eventuellement fixer le numero de version du paquet
rpm -qa | grep cfitsio || sudo yum install cfitsio


# ======================================== #
# Copy package configuration file (*.spec) #
# ======================================== #
cp $PACKNAME.spec $WRKDIR/SPECS/

# ================================================================== #
# Install GammaLib
# No need to install Gammalib since this script follows a 'make dist'#
# ================================================================== #
# Check gammalib version
rpm --qf '%{VERSION}\n' -q gammalib

# ================================================================== #
# Install ctools                                                     #
# No need to install ctools ; rpmbuild will package from 'scratch'   #
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

rpmbuild -ba $WRKDIR/SPECS/$PACKNAME.spec -v

echo " ===================================================="
echo " Packaging DONE :"
echo " New package is located at $PRODDIR"
echo " rpm file : $PACKAGE-$RELEASE.noarch.rpm"
echo " ===================================================="
echo "  "
rpm -qip $RPMFILE

exit
