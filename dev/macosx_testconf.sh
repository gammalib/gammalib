#!/bin/bash -f
#############################################################################
# Test Mac OS X configuration                                               #
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
# This script tests GammaLib under various Mac OS X configurations. The     #
# main usage is for consolidation the ./configure script.                   #
#                                                                           #
# Here the various tests and the expected outcome:                          #
#                                                                           #
# Case A: Vanilla Mac OS X                                                  #
#         A1 - 1 of 15 tests failed                                         #
#         A2 - All 15 tests passed                                          #
#         A3 - All 15 tests passed                                          #
#         A4 - All 15 tests passed                                          #
#                                                                           #
# Case B: Pure MacPorts install                                             #
#         B1 - All 15 tests passed                                          #
#         B2 - All 15 tests passed                                          #
#                                                                           #
# Case C: MacPorts cfitsio and Mac OS X Python                              #
#         C1 - 1 of 15 tests failed                                         #
#         C2 - 1 of 15 tests failed                                         #
#         C3 - All 15 tests passed                                          #
#         C4 - All 15 tests passed                                          #
#                                                                           #
# Case D: MacPorts cfitsio and www.python.org                               #
#         D1 - 1 of 15 tests failed                                         #
#         D2 - 1 of 15 tests failed                                         #
#         D3 - All 15 tests passed                                          #
#         D4 - All 15 tests passed                                          #
#                                                                           #
# Case E: Using 3-way cfitsio                                               #
#         E1 - 1 of 15 tests failed                                         #
#         E2 - All 15 tests passed                                          #
#         E3 - All 15 tests passed                                          #
#         E4 - All 15 tests passed                                          #
#                                                                           #
#############################################################################

# Set working directory
wrkdir=/Users/jurgen/dev/gammalib

# Step into working directory
cd $wrkdir

# Set tests to be executed
do_A1="no"
do_A2="no"
do_A3="no"
do_A4="no"
do_B1="no"
do_B2="no"
do_C1="no"
do_C2="no"
do_C3="no"
do_C4="no"
#
do_D1="yes"
do_D2="yes"
do_D3="yes"
do_D4="yes"
do_E1="yes"
do_E2="yes"
do_E3="yes"
do_E4="yes"


#############################################################################
# Vanilla Mac OS X                                                          #
#############################################################################

# Case A.1: Python 2.5
if [ "$do_A1" = "yes" ]; then
    echo "Vanilla Build, Python 2.5 (./configure)"
    rm -rf tmp
    mkdir tmp
    cp -r /Users/jurgen/Software/cfitsio/cfitsio3280/intel/ tmp/
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.5/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --prefix=$wrkdir/tmp > CaseA1_conf.log 2>&1
    make -j10 > CaseA1_make.log 2>&1
    make check  > CaseA1_check.log 2>&1
fi

# Case A.2: Python 2.5
if [ "$do_A2" = "yes" ]; then
    echo "Vanilla Build, Python 2.5 (./configure --enable-universalsdk)"
    rm -rf tmp
    mkdir tmp
    cp -r /Users/jurgen/Software/cfitsio/cfitsio3280/intel/ tmp/
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.5/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --prefix=$wrkdir/tmp --enable-universalsdk > CaseA2_conf.log 2>&1
    make -j10 > CaseA2_make.log 2>&1
    make check  > CaseA2_check.log 2>&1
fi

# Case A.3: Python 2.6
if [ "$do_A3" = "yes" ]; then
    echo "Vanilla Build, Python 2.6 (./configure)"
    rm -rf tmp
    mkdir tmp
    cp -r /Users/jurgen/Software/cfitsio/cfitsio3280/intel/ tmp/
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.6/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --prefix=$wrkdir/tmp > CaseA3_conf.log 2>&1
    make -j10 > CaseA3_make.log 2>&1
    make check  > CaseA3_check.log 2>&1
fi

# Case A.4: Python 2.6
if [ "$do_A4" = "yes" ]; then
    echo "Vanilla Build, Python 2.6 (./configure --enable-universalsdk)"
    rm -rf tmp
    mkdir tmp
    cp -r /Users/jurgen/Software/cfitsio/cfitsio3280/intel/ tmp/
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.6/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --prefix=$wrkdir/tmp --enable-universalsdk > CaseA4_conf.log 2>&1
    make -j10 > CaseA4_make.log 2>&1
    make check  > CaseA4_check.log 2>&1
fi


#############################################################################
# Pure MacPorts                                                             #
#############################################################################

# Case B.1: Pure MacPorts
if [ "$do_B1" = "yes" ]; then
    echo "Pure MacPorts build (./configure)"
    export PATH=/opt/local/Library/Frameworks/Python.framework/Versions/2.6/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure > CaseB1_conf.log 2>&1
    make -j10 > CaseB1_make.log 2>&1
    make check  > CaseB1_check.log 2>&1
fi

# Case B.2: Pure MacPorts
if [ "$do_B2" = "yes" ]; then
    echo "Pure MacPorts build (./configure --enable-universalsdk)"
    export PATH=/opt/local/Library/Frameworks/Python.framework/Versions/2.6/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --enable-universalsdk > CaseB2_conf.log 2>&1
    make -j10 > CaseB2_make.log 2>&1
    make check  > CaseB2_check.log 2>&1
fi


#############################################################################
# MacPorts cfitsio & Mac OS X Python                                        #
#############################################################################

# Case C.1:
if [ "$do_C1" = "yes" ]; then
    echo "MacPorts cfitsio, Mac OS X Python 2.5 (./configure)"
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.5/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure > CaseC1_conf.log 2>&1
    make -j10 > CaseC1_make.log 2>&1
    make check  > CaseC1_check.log 2>&1
fi

# Case C.2:
if [ "$do_C2" = "yes" ]; then
    echo "MacPorts cfitsio, Mac OS X Python 2.5 (./configure --enable-universalsdk)"
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.5/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --enable-universalsdk > CaseC2_conf.log 2>&1
    make -j10 > CaseC2_make.log 2>&1
    make check  > CaseC2_check.log 2>&1
fi

# Case C.3:
if [ "$do_C3" = "yes" ]; then
    echo "MacPorts cfitsio, Mac OS X Python 2.6 (./configure)"
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.6/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure > CaseC3_conf.log 2>&1
    make -j10 > CaseC3_make.log 2>&1
    make check  > CaseC3_check.log 2>&1
fi

# Case C.4:
if [ "$do_C4" = "yes" ]; then
    echo "MacPorts cfitsio, Mac OS X Python 2.6 (./configure --enable-universalsdk)"
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.6/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --enable-universalsdk > CaseC4_conf.log 2>&1
    make -j10 > CaseC4_make.log 2>&1
    make check  > CaseC4_check.log 2>&1
fi


#############################################################################
# MacPorts cfitsio & www.python.org                                         #
#############################################################################

# Case D.1:
if [ "$do_D1" = "yes" ]; then
    echo "MacPorts cfitsio, www.python.org 32-bit (./configure)"
    export PATH=/Users/jurgen/Software/python/python27/10.4/32-bit/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure > CaseD1_conf.log 2>&1
    make -j10 > CaseD1_make.log 2>&1
    make check  > CaseD1_check.log 2>&1
fi

# Case D.2:
if [ "$do_D2" = "yes" ]; then
    echo "MacPorts cfitsio, www.python.org 32-bit (./configure --enable-universalsdk)"
    export PATH=/Users/jurgen/Software/python/python27/10.4/32-bit/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --enable-universalsdk > CaseD2_conf.log 2>&1
    make -j10 > CaseD2_make.log 2>&1
    make check  > CaseD2_check.log 2>&1
fi

# Case D.3:
if [ "$do_D3" = "yes" ]; then
    echo "MacPorts cfitsio, www.python.org intel (./configure)"
    export PATH=/Users/jurgen/Software/python/python27/10.6/intel/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure > CaseD3_conf.log 2>&1
    make -j10 > CaseD3_make.log 2>&1
    make check  > CaseD3_check.log 2>&1
fi

# Case D.4:
if [ "$do_D4" = "yes" ]; then
    echo "MacPorts cfitsio, www.python.org intel (./configure --enable-universalsdk)"
    export PATH=/Users/jurgen/Software/python/python27/10.6/intel/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --enable-universalsdk > CaseD4_conf.log 2>&1
    make -j10 > CaseD4_make.log 2>&1
    make check  > CaseD4_check.log 2>&1
fi


#############################################################################
# 3-way cfitsio                                                             #
#############################################################################

# Case E.1: Python 2.5
if [ "$do_E1" = "yes" ]; then
    echo "Vanilla Build, Python 2.5 (./configure)"
    rm -rf tmp
    mkdir tmp
    cp -r /Users/jurgen/Software/cfitsio/cfitsio3280/3-way/ tmp/
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.5/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --prefix=$wrkdir/tmp > CaseE1_conf.log 2>&1
    make -j10 > CaseE1_make.log 2>&1
    make check  > CaseE1_check.log 2>&1
fi

# Case E.2: Python 2.5
if [ "$do_E2" = "yes" ]; then
    echo "Vanilla Build, Python 2.5 (./configure --enable-universalsdk)"
    rm -rf tmp
    mkdir tmp
    cp -r /Users/jurgen/Software/cfitsio/cfitsio3280/3-way/ tmp/
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.5/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --prefix=$wrkdir/tmp --enable-universalsdk > CaseE2_conf.log 2>&1
    make -j10 > CaseE2_make.log 2>&1
    make check  > CaseE2_check.log 2>&1
fi

# Case E.3: Python 2.6
if [ "$do_E3" = "yes" ]; then
    echo "Vanilla Build, Python 2.6 (./configure)"
    rm -rf tmp
    mkdir tmp
    cp -r /Users/jurgen/Software/cfitsio/cfitsio3280/3-way/ tmp/
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.6/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --prefix=$wrkdir/tmp > CaseE3_conf.log 2>&1
    make -j10 > CaseE3_make.log 2>&1
    make check  > CaseE3_check.log 2>&1
fi

# Case E.4: Python 2.6
if [ "$do_E4" = "yes" ]; then
    echo "Vanilla Build, Python 2.6 (./configure --enable-universalsdk)"
    rm -rf tmp
    mkdir tmp
    cp -r /Users/jurgen/Software/cfitsio/cfitsio3280/3-way/ tmp/
    export PATH=/System/Library/Frameworks/Python.framework/Versions/2.6/bin:/usr/bin:/bin:/usr/sbin:/sbin
    make distclean > /dev/null 2>&1
    ./configure --prefix=$wrkdir/tmp --enable-universalsdk > CaseE4_conf.log 2>&1
    make -j10 > CaseE4_make.log 2>&1
    make check  > CaseE4_check.log 2>&1
fi
