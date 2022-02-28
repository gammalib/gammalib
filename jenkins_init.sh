#!/bin/sh
# =====================================================================
# Initialise for Jenkins Continuous Integration
#
# Copyright (C) 2022 Juergen Knoedlseder
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
# =====================================================================
#
# Set parameters
# ==============
lock="$HOME/.gamma/daemon.lock"
heartbeat="$HOME/.gamma/daemon.heartbeat"

#
# Kill daemon if it is running
# ============================
if [ -f "$lock" ]; then
    pid=`cat $lock`
    ps cax | grep gammalibd > /dev/null
    if [ $? -eq 0 ]; then
      kill -9 $pid
      rm -f $lock
      rm -f $heartbeat
    fi
fi

#
# Remove heartbeat file
# =====================
if [ -f "$heartbeat" ]; then
    rm -f $heartbeat
fi

