#!/bin/bash
# ==========================================================================
# gammalib CentOS package validation
#
# Copyright (C) 2017 Sylvie Brau-Nogue
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


# =============================== #
# Set software component versions #
# =============================== #
PACKNAME=gammalib
PACKAGE=$PACKNAME-$VERSION
RELEASE=1


# ====== #
# Header #
# ====== #
echo "================================================================="
echo "Check $PACKNAME rpm package, Version $VERSION, Platform $PLATFORM"
echo "================================================================="


# ============== #
# Set parameters #
# ============== #
WRKDIR=$PWD/rpm_check
RPMDB=$WRKDIR/rpmdb
RPMFILE=$PACKAGE-$RELEASE.*$PLATFORM.rpm


# ========================= #
# Create local rpm database #
# ========================= #
rpm --initdb --root $RPMDB


# ======================== #
# Install GammaLib package #
# ======================== #
rpm --dbpath $RPMDB  --nodeps --prefix=$WRKDIR -i $RPMFILE


# ========================== #
# Configure GammaLib package #
# ========================== #
export LD_LIBRARY_PATH=$WRKDIR/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$WRKDIR/lib64/python2.7/site-packages:$PYTHONPATH


# ============ #
# Test package #
# ============ #
python -c 'import gammalib; gammalib.test()'
