#!/bin/bash -f
#############################################################################
# Test Mac OS X package                                                     #
# ------------------------------------------------------------------------- #
# Copyright (C) 2011 Jurgen Knodlseder                                      #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
# ------------------------------------------------------------------------- #
# This script tests a Mac OS X package prior to release. It is assumed that #
# the package has been installed in ${GAMMALIB}.                            #
#############################################################################

# Set GAMMALIB environment variable
export GAMMALIB=$HOME/ctools
#GAMMALIB=$HOME/test/gamma
#GAMMALIB=/usr/local/gamma

# Set supported architectures (for command line tests)
ARCHS="ppc i386 x86_64"

# Set supported Python versions
declare -a PYTHONS
let count=0
PYTHONS[$count]=/System/Library/Frameworks/Python.framework/Versions/2.5/bin
((count++))
PYTHONS[$count]=/System/Library/Frameworks/Python.framework/Versions/2.6/bin
((count++))
PYTHONS[$count]=/Users/jurgen/Software/python/python27/10.5/3-way/bin
((count++))

# Loop over Python versions
let index=0
while [ $index -lt $count ]; do

  # Set PATH
  export PATH=${PYTHONS[$index]}:/opt/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
  
  # Reset PYTHONPATH
  export PYTHONPATH=
  
  # Source initialisation script
  source $GAMMALIB/bin/gammalib-init.sh

  # Determine Python version
  PYTHON_VERSION=`python -c "import sys, string; \
                  ver=string.split(sys.version)[0].split('.'); \
                  print ver[0]+'.'+ver[1]"`

  # Determine Python architectures
  PYTHON_BIN=`which python`
  PYTHON_ARCHS=
  if [ "x`file $PYTHON_BIN | grep ppc`" != x ]; then
    PYTHON_ARCHS=$PYTHON_ARCHS"ppc "
  fi
  if [ "x`file $PYTHON_BIN | grep i386`" != x ]; then
    PYTHON_ARCHS=$PYTHON_ARCHS"i386 "
  fi
  if [ "x`file $PYTHON_BIN | grep x86_64`" != x ]; then
    PYTHON_ARCHS=$PYTHON_ARCHS"x86_64 "
  fi
  
  # Loop over Python architectures
  let goods=0
  let bads=0
  ARCH_GOOD=
  ARCH_BAD=
  for ARCH in $PYTHON_ARCHS; do
  
    # Reset error counter
    let errors=0

    # Check if gammalib import works
    RES=`arch -$ARCH python -c "from gammalib import *" 2>&1`
    if [ "x$RES" != x ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): Unable to import gammalib ($RES)"
      ((errors++))
    fi

    # Check if ctools import works
    RES=`arch -$ARCH python -c "from ctools import *" 2>&1`
    if [ "x$RES" != x ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): Unable to import ctools ($RES)"
      ((errors++))
    fi

    # Check if gammalib and ctools import works
    RES=`arch -$ARCH python -c "from ctools import *; from gammalib import *" 2>&1`
    if [ "x$RES" != x ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): Unable to import ctools & gammalib ($RES)"
      ((errors++))
    fi
    
    # Check if GModels() works
    RES=`arch -$ARCH python -c "from gammalib import *; \
                                m=GModels(); \
                                print m" 2>&1`
    if [ "x`echo $RES | grep Traceback`" != x ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): $RES"
      ((errors++))
    fi
    
    # Remove existing result files
    rm -f events.fits
    rm -f cntmap.fits
    rm -f selected_events.fits
    rm -f results_unbinned.xml

    # Check if ctobssim works in Python
    RES=`arch -$ARCH python -c "from ctools import *; \
                                from gammalib import *; \
                                sim = ctobssim(); \
                                sim['infile'].filename('$GAMMALIB/share/models/crab.xml'); \
                                sim['outfile'].filename('events.fits'); \
                                sim['caldb'].string('$GAMMALIB/share/caldb/cta'); \
                                sim['irf'].string('kb_E_50h_v3'); \
                                sim['ra'].real(83.63); \
                                sim['dec'].real(22.01); \
                                sim['rad'].real(10.0); \
                                sim['tmin'].real(0.0); \
                                sim['tmax'].real(1800.0); \
                                sim['emin'].real(0.1); \
                                sim['emax'].real(100.0); \
                                sim.execute()" 2>&1`
    if [ "x`echo $RES | grep Traceback`" != x ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): $RES"
      ((errors++))
    fi
    if [ ! -s "events.fits" ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): ctobssim - No valid events.fits file."
      ((errors++))
    fi

    # Check if ctbin works in Python
    RES=`arch -$ARCH python -c "from ctools import *; \
                                from gammalib import *; \
                                bin = ctbin(); \
                                bin['evfile'].filename('events.fits'); \
                                bin['outfile'].filename('cntmap.fits'); \
                                bin['emin'].real(0.1); \
                                bin['emax'].real(100.0); \
                                bin['enumbins'].integer(20); \
                                bin['nxpix'].integer(200); \
                                bin['nypix'].integer(200); \
                                bin['binsz'].real(0.02); \
                                bin['coordsys'].string('CEL'); \
                                bin['xref'].real(83.63); \
                                bin['yref'].real(22.01); \
                                bin['proj'].string('CAR'); \
                                bin.execute()" 2>&1`
    if [ "x`echo $RES | grep Traceback`" != x ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): $RES"
      ((errors++))
    fi
    if [ ! -s "cntmap.fits" ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): ctbin - No valid cntmap.fits file."
      ((errors++))
    fi

    # Check if ctselect works in Python
    RES=`arch -$ARCH python -c "from ctools import *; \
                                from gammalib import *; \
                                select = ctselect(); \
                                select['infile'].filename('events.fits'); \
                                select['outfile'].filename('selected_events.fits'); \
                                select['ra'].real(83.63); \
                                select['dec'].real(22.01); \
                                select['rad'].real(4.0); \
                                select['emin'].real(0.1); \
                                select['emax'].real(100.0); \
                                select['tmin'].real(0.0); \
                                select['tmax'].real(1800); \
                                select.execute()" 2>&1`
    if [ "x`echo $RES | grep Traceback`" != x ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): $RES"
      ((errors++))
    fi
    if [ ! -s "selected_events.fits" ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): ctselect - No valid selected_events.fits file."
      ((errors++))
    fi

    # Check if ctlike works in Python
    RES=`arch -$ARCH python -c "from ctools import *; \
                                from gammalib import *; \
                                like = ctlike(); \
                                like['evfile'].filename('selected_events.fits'); \
                                like['srcmdl'].filename('$GAMMALIB/share/models/crab.xml'); \
                                like['outmdl'].filename('results_unbinned.xml'); \
                                like['caldb'].string('$GAMMALIB/share/caldb/cta'); \
                                like['irf'].string('kb_E_50h_v3'); \
                                like['method'].string('UNBINNED'); \
                                like.execute()" 2>&1`
    if [ "x`echo $RES | grep Traceback`" != x ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): $RES"
      ((errors++))
    fi
    if [ ! -s "results_unbinned.xml" ]; then
      echo "ERROR (Python$PYTHON_VERSION-$ARCH): ctlike - No valid results_unbinned.xml file."
      ((errors++))
    fi
    
    # Collect good and bad Python architectures
    if [ $errors == 0 ]; then
      ARCH_GOOD=$ARCH_GOOD" "$ARCH
      ((goods++))
    else
      ARCH_BAD=$ARCH_BAD" "$ARCH
      ((bads++))
    fi
  
  done
  
  # Signal success and failures
  if [ $goods -gt 0 ]; then
    echo "Python$PYTHON_VERSION support works for$ARCH_GOOD."
  fi
  if [ $bads -gt 0 ]; then
    echo "Python$PYTHON_VERSION support fails for$ARCH_BAD."
  fi
  
  # Increment index
  ((index++))
  
done

#
# Now test the command line interface
#====================================

# Source initialisation script
source $GAMMALIB/bin/gammalib-init.sh

# Create pfiles directory and set environment
rm -rf pfiles
mkdir -p pfiles
export PFILES=pfiles

# Loop over architectures
for ARCH in $ARCHS; do

  # Remove existing result files
  rm -f events.fits
  rm -f cntmap.fits
  rm -f selected_events.fits
  rm -f results_unbinned.xml

  # Check if ctobssim works from command line
  arch -$ARCH ctobssim \
         infile="$GAMMALIB/share/models/crab.xml" \
         outfile="events.fits" \
         caldb="$GAMMALIB/share/caldb/cta" \
         irf="kb_E_50h_v3" \
         ra=83.63 \
         dec=22.01 \
         rad=10.0 \
         tmin=0.0 \
         tmax=1800.0 \
         emin=0.1 \
         emax=100.0
  if [ -s "events.fits" ]; then
    echo "ctobssim($ARCH): ok."
  else
    echo "ERROR ($ARCH): ctobssim failed."
  fi

  # Check if ctbin works from command line
  arch -$ARCH ctbin \
       evfile="events.fits" \
       outfile="cntmap.fits" \
       emin=0.1 \
       emax=100.0 \
       enumbins=20 \
       nxpix=200 \
       nypix=200 \
       binsz=0.02 \
       coordsys="CEL" \
       xref=83.63 \
       yref=22.01 \
       proj="CAR"
  if [ -s "cntmap.fits" ]; then
    echo "ctbin($ARCH): ok."
  else
    echo "ERROR ($ARCH): ctbin failed."
  fi

  # Check if ctselect works from command line
  arch -$ARCH ctselect \
          infile="events.fits" \
          outfile="selected_events.fits" \
          ra=83.63 \
          dec=22.01 \
          rad=10.0 \
          tmin=0.0 \
          tmax=1800.0 \
          emin=0.1 \
          emax=100.0
  if [ -s "selected_events.fits" ]; then
    echo "ctselect($ARCH): ok."
  else
    echo "ERROR ($ARCH): ctselect failed."
  fi

  # Check if ctlike works from command line
  arch -$ARCH ctlike evfile="selected_events.fits" \
         srcmdl="$GAMMALIB/share/models/crab.xml" \
         outmdl="results_unbinned.xml" \
         method="UNBINNED" \
         caldb="$GAMMALIB/share/caldb/cta" \
         irf="kb_E_50h_v3"
  if [ -s "results_unbinned.xml" ]; then
    echo "ctlike($ARCH): ok."
  else
    echo "ERROR ($ARCH): ctlike failed."
  fi

done

