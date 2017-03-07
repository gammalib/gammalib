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
# This script tests a CentOS package prior to release. It is assumed that
# the package has been installed in a local directory $LOCALINSTALL/gammalib.
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

if ! rpm -qa | grep -q ncurses; then
  echo "ncurses need to be installed"
  sudo yum install ncurses
  if ! rpm -qa | grep -q ncurses; then
     echo "================================================================================================="
     echo "ERROR : ncurses need to be installed manually"
     echo "================================================================================================="
     exit
  fi
fi
v=$(rpm --qf '%{VERSION}\n' -q ncurses)
echo "ncurses $v checked"

if ! rpm -qa | grep -q readline; then
  echo "readline need to be installed"
  sudo yum install readline
  if ! rpm -qa | grep -q readline; then
     echo "================================================================================================="
     echo "ERROR : readline need to be installed manually"
     echo "================================================================================================="
     exit
  fi
fi
v=$(rpm --qf '%{VERSION}\n' -q readline)
echo "readline $v checked"

if ! rpm -qa | grep -q cfitsio; then
   echo "cfitsio need to be installed"
   sudo yum install cfitsio-devel
   if ! rpm -qa | grep -q cfitsio; then
      echo "================================================================================================="
      echo " ERROR : cfitsio need to be installed"
      echo "================================================================================================="
      echo "You can download the latest version of CFITSIO on http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html"
      echo "   "
      echo "or "
      echo "   "
      echo "cfitsio rpm can be found in the epel-release repository ; Try : sudo yum install cfitsio-devel"
      echo "   "
      echo "================================================================================================="
      exit
   fi
fi
v=$(rpm --qf '%{VERSION}\n' -q cfitsio)
echo "cfitsio $v checked"

if ! sudo rpm -U --prefix=$LOCALINSTALL $RPMFILE &>  $LOGRPMFILE ; then
  echo "Cannot install gammalib package"
  cat  $LOGRPMFILE
  exit
fi

# ================= #
# Configure package #
# ================= #
export GAMMALIB=$LOCALINSTALL/gammalib
source $GAMMALIB/bin/gammalib-init.sh
ARCH=`lsb_release -irs`

# ============ #
# Test package #
# ============ #
# TODO : update following section with consistent set of tests
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

