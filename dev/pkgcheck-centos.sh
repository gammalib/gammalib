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
echo " Check $PACKNAME rpm package, Version $VERSION, Platform $PLATFORM "
echo " ================================================================= "

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
# TODO : next step to fix RELEASE is not updated for now
RELEASE=1

# ==================================================== #
#  TODO : Create a secure RPM repository
# Temporary solution, to be fixed
# ==================================================== #
YUMREPO=$HOME/YUM-REPO

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
WRKDIR=$HOME/pkg_check
SRCDIR=$WRKDIR/SOURCES
PKGDIR=$WRKDIR/RPMS

# ========================================================================= #
# Fix packages directory                                                    #
# ========================================================================= #
# If 'make pkg' produce multi-platform rpm, rpm files will be archived in $YUMREPO directory
#PRODDIR=$WRKDIR/RPMS/$PLATFORM
PRODIR=$YUMREPO
#
RPMFILE=$PRODIR/$PACKAGE-$RELEASE.*$PLATFORM.rpm
#
# ========================================================================= #
# Install package in local directory                                        #
# ========================================================================= #
# 2 options with rpm                                                        #
#  --relocate <old>=<new>    relocate files from path <old> to <new>        #
#sudo rpm -U --relocate=$INSTALLDIR=$LOCALINSTALL $RPMFILE &>  $LOGRPMFILE  #
# or                                                                        #
#  --prefix=<dir>            relocate the package to <dir>, if relocatable  #
# ========================================================================= #
INSTALLDIR=/usr/local/gamma
# LOCALINSTALL=$WRKDIR$INSTALLDIR  # <= another option, not selected for the moment
LOCALINSTALL=$HOME/$INSTALLDIR
#
LOGFILE=$WRKDIR/pkg_check.log
LOGRPMFILE=$WRKDIR/pkg_check_rpm.log

# ========================= #
# Clean temporary directory #
# ========================= #
rm -rf $TMPDIR

# ============================= #
# Secure installation directory #
# ============================= #
if [ -d "$LOCALINSTALL/$PACKNAME" ]; then
    mv $LOCALINSTALL/$PACKNAME $LOCALINSTALL/$PACKNAME.backup
fi

# ======================= #
# Clean working directory #
# ======================= #
sudo rm -rf $WRKDIR

# ====================================== #
# Create and step into working directory #
# ====================================== #
mkdir -p $WRKDIR
cd $WRKDIR

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

# ======================================================================= #
# Finally install package $HOME/$YUMREPO/$PACKAGE-$RELEASE-*$PLATFORM.rpm #
# ======================================================================= #
if ! sudo rpm -U --prefix=/$LOCALINSTALL $RPMFILE &>  $LOGRPMFILE ; then
   echo "================================================================================================="
   echo "ERROR : Cannot install $PACKAGE package"
   echo "================================================================================================="
   cat  $LOGRPMFILE
   echo "   "
   echo "================================================================================================="
   exit
fi

# ================= #
# Configure package #
# ================= #
export GAMMALIB=$LOCALINSTALL/$PACKNAME
source $GAMMALIB/bin/$PACKNAME-init.sh
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
sudo rm -rf $TMPDIR

# ============================== #
# Recover installation directory #
# ============================== #
if [ -d "$LOCALINSTALL/$PACKNAME.backup" ]; then
    mv $LOCALINSTALL/$PACKNAME.backup $LOCALINSTALL/$PACKNAME
fi

