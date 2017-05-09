#!/bin/bash
# ==========================================================================
# gammalib Ubuntu package validation
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
# This script tests a Ubuntu package prior to release. It is assumed that
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
# TODO : next step to fix RELEASE is not updated for now
RELEASE=1

# ==================================================== #
#  TODO : Create a secure RPM repository
# ==================================================== #
# Temporary solution, to be fixed
RPMREPO=$HOME/RPM-REPO
[ ! -d $RPMREPO ] && mkdir -p $RPMREPO

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
CHECKDIR=$HOME/usr
WRKDIR=$CHECKDIR/pkg_check
#SRCDIR=$WRKDIR/SOURCES
#PKGDIR=$WRKDIR/RPMS

# ========================================================================= #
# Fix packages directory                                                    #
# ========================================================================= #
# If 'make pkg' produce multi-platform rpm, rpm files will be archived in $RPMREPO directory
PRODDIR=$RPMREPO
#

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
LOGDEBFILE=$WRKDIR/pkg_check_deb.log

# ============================= #
# Secure installation directory #
# ============================= #
if [ -d "$LOCALINSTALL/$PACKNAME" ]; then
echo "Secure installation directory : $LOCALINSTALL/$PACKNAME"
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

# =============== #
# Install ncurses #
# =============== #
# TODO : eventuellement fixer le numero de version du paquet

echo "================================================================================================="
echo "Check ncurses install"
echo "================================================================================================="
if ! dpkg -l | grep '^ii' | grep libncurses | awk '{print $2}' ; then
    echo "ncurses need to be installed"
    sudo apt-get install libncurses-dev
    if ! dpkg -l | grep '^ii' | grep libncurses | awk '{print $2}'; then
        echo "================================================================================================="
        echo "ERROR : ncurses need to be installed manually"
        echo "   "
        echo "================================================================================================="
        exit
    fi
fi

# Check ncurses version
p=$(dpkg -l | grep '^ii' | grep libncurses | grep dev | awk '{print $2}')
if [ $(dpkg-query -W -f='${Status}' $p 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
  echo "$p is definitely not installed" ;
  echo "Next step should be : apt-get install $p";
  exit
fi
v=$(dpkg-query -W -f='${Version}\n' $p)
echo "Package $p (version : $v) checked"

# ================ #
# Install readline #
# ================ #
# TODO : eventuellement fixer le numero de version du paquet

echo "================================================================================================="
echo "Check readline install"
echo "================================================================================================="
if ! dpkg -l | grep '^ii' | grep libreadline | awk '{print $2}' ; then
    echo "readline need to be installed"
    sudo apt-get install libreadline6 libreadline6-dev
    if ! dpkg -l | grep '^ii' | grep libreadline | awk '{print $2}'; then
        echo "================================================================================================="
        echo "ERROR : readline need to be installed manually"
        echo "   "
        echo "================================================================================================="
        exit
    fi
fi
# Check readline version
p=$(dpkg -l | grep '^ii' | grep libreadline | grep dev | awk '{print $2}')

if [ $(dpkg-query -W -f='${Status}' $p 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
  echo "$p is definitely not installed" ;
  echo "Next step should be : apt-get install $p";
  exit
fi
v=$(dpkg-query -W -f='${Version}\n' $p)
echo "Package $p (version : $v) checked"

# =============== #
# Install cfitsio #
# =============== #
# TODO : eventuellement fixer le numero de version du paquet

echo "================================================================================================="
echo "Check cfitsio install"
echo "================================================================================================="
if ! dpkg -l | grep '^ii' | grep libcfitsio | awk '{print $2}' ; then
    echo "cfitsio need to be installed"
    apt-cache search cfitsio
    sudo apt-get install libcfitsio-dev
    if ! dpkg -l | grep '^ii' | grep libcfitsio | awk '{print $2}'; then
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
# Check cfitsio version
p=$(dpkg -l | grep '^ii' | grep libcfitsio | grep dev | awk '{print $2}')

if [ $(dpkg-query -W -f='${Status}' $p 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
  echo "$p is definitely not installed" ;
  echo "Next step should be : apt-get install $p";
  exit
fi
v=$(dpkg-query -W -f='${Version}\n' $p)
echo "Package $p (version : $v) checked"


# ======================================================================= #
# Finally install package $HOME/$RPMREPO/$DEBPKG-$RELEASE_$PLATFORM.deb  #
# ======================================================================= #
debformat=true # Check with Debian package


if $debformat ; then
   DEBPKG=$PACKNAME\_$VERSION
   DEBFILE=$PRODDIR/$DEBPKG-$RELEASE\_$PLATFORM.deb

   echo "Install and Check DEB Package : $DEBFILE"
   if ! dpkg-deb -x $DEBFILE $CHECKDIR &>  $LOGDEBFILE ; then
      echo "================================================================================================="
      echo "ERROR : Cannot install $DEBFILE package"
      echo "================================================================================================="
      cat  $LOGDEBFILE
      echo "   "
      echo "================================================================================================="
   fi
else
   # ======================================================================= #
   # Or install package $HOME/$RPMREPO/$RPMPKG-$RELEASE.$PLATFORM.rpm #
   # ======================================================================= #
   #if ! sudo rpm -U --prefix=/$LOCALINSTALL $RPMFILE &>  $LOGRPMFILE ; then
   RPMPKG=$PACKNAME-$VERSION
   RPMFILE=$PRODDIR/$RPMPKG-$RELEASE.$PLATFORM.rpm
   echo "Install and Check RPM Package : $RPMFILE"
   if ! alien -iv $RPMFILE &>  $LOGRPMFILE ; then
      echo "================================================================================================="
      echo "ERROR : Cannot install $RPMFILE package"
      echo "================================================================================================="
      cat  $LOGRPMFILE
      echo "   "
      echo "================================================================================================="
   fi
fi

# ================= #
# Configure package #
# ================= #
export GAMMALIB=$LOCALINSTALL/$PACKNAME
source $GAMMALIB/bin/$PACKNAME-init.sh

# ============ #
# Test package #
# ============ #
# TODO : update following section with consistent set of tests
python -c 'import gammalib; gammalib.test()' | tee -a $LOGFILE


# ======================= #
# Uninstall package       #
# ======================= #
sudo rm -rf $LOCALINSTALL/$PACKNAME

# ============================== #
# Recover installation directory #
# ============================== #
if [ -d "$LOCALINSTALL/$PACKNAME.backup" ]; then
    mv $LOCALINSTALL/$PACKNAME.backup $LOCALINSTALL/$PACKNAME
fi

