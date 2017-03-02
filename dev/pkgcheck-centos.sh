#!/bin/bash
# ==========================================================================
# gammalib CentOS package validation
#
# Copyright (C) 2017 Sylvie Brau-Nogué
#
# This program is free software: you can redistribute it and/or modify
# it under the tesudo rms of the GNU General Public License as published by
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
# This script checks the ctools CentOS package by installing it into the
# /usr/local/gamma directory and running the Python tests.
#
# ==========================================================================

if [ $# -ne 2 ] ;  then
echo "Please specify package name and version in the form 'package x.y.z'"
exit 1
fi
PACKNAME=$1
VERSION=$2
# TODO : next step to fix RELEASE is not updated for now
RELEASE=1


# =============================== #
# Set software component versions #
# =============================== #
PACKAGE=$PACKNAME-$VERSION-$RELEASE


# ============== #
# Set parameters #
# ============== #
INSTALLDIR=/usr/local/gamma
WRKDIR=$HOME/pkg_check
LOCALINSTALL=$WRKDIR$INSTALLDIR
#
PKGDIR=$HOME/rpmbuild/RPMS/noarch
RPMFILE=$PKGDIR/$PACKAGE.noarch.rpm
#
LOGFILE=$WRKDIR/pkg_check.log
LOGRPMFILE=$WRKDIR/pkg_check_rpm.log


# ============================= #
# Secure installation directory #
# ============================= #
# SBN : pas besoin              #
# ============================= #

#if [ -d "$LOCALINSTALL/ctools" ]; then
#    mv $LOCALINSTALL/ctools $LOCALINSTALL/ctools.backup
#fi


# ======================= #
# Clean working directory #
# ======================= #
sudo rm -rf $WRKDIR

# ====================================== #
# Create and step into working directory #
# ====================================== #
mkdir -p $WRKDIR
cd $WRKDIR

# ========================================================================= #
# Install package in local directory                                        #
# ========================================================================= #
# 2 options with rpm                                                        #
#  --relocate <old>=<new>    relocate files from path <old> to <new>        #
#sudo rpm -U --relocate=$INSTALLDIR=$LOCALINSTALL $RPMFILE &>  $LOGRPMFILE  #
# or                                                                        #
#  --prefix=<dir>            relocate the package to <dir>, if relocatable  #
# ========================================================================= #

sudo rpm -U --prefix=$LOCALINSTALL $RPMFILE &>  $LOGRPMFILE


# ================= #
# Configure package #
# ================= #
export GAMMALIB=$LOCALINSTALL/gammalib
source $GAMMALIB/bin/gammalib-init.sh


# ============ #
# Test package #
# ============ #
python -c 'import gammalib; gammalib.test()' | tee -a $LOGFILE


# ======================= #
# Uninstall package       #
# ======================= #
sudo yum remove $PACKNAME

# ======================= #
# Clean package directory #
# ======================= #
sudo rm -rf $LOCALINSTALL

# ============================== #
# Recover installation directory #
# ============================== #
# SBN : pas besoin              #
# ============================= #

#if [ -d "$LOCALINSTALL/ctools.backup" ]; then
#    mv $LOCALINSTALL/ctools.backup $LOCALINSTALL/ctools
#fi
