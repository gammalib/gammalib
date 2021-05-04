#!/usr/bin/python
# ==========================================================================
# Generate sloccount history (requires sloccount to be installed)
#
# Copyright (C) 2021 Juergen Knoedlseder
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
# ==========================================================================
import subprocess
import gammalib


# ========================== #
# Generate sloccount history #
# ========================== #
def generate_history():
    """
    Generate sloccount history
    """
    # Get commit history
    git        = subprocess.Popen(['git', 'log',
                                   '--pretty="%H %cd"',
                                   '--date=format:"%Y-%m-%dT%H:%M:%S"',
                                   '--reverse'], stdout=subprocess.PIPE)
    history, _ = git.communicate()

    # Open result file
    result = open('sloc.dat', 'w')

    # Loop over commit history
    for line in history.split('\n'):

        # Skip invalid lines
        if not line:
            continue

        # Get commit and date
        hash, datestr = line.split(' ', 1)
        hash          = hash.strip('"')
        datestr       = datestr.strip('"')

        # Checkout commit
        git = subprocess.Popen(['git', 'checkout', '-f', hash], stdout=subprocess.PIPE)
        commit, err = git.communicate()
        #print(commit, err)

        # Get SLOC count statistics
        git       = subprocess.Popen(['sloccount', '.'], stdout=subprocess.PIPE)
        sloc, err = git.communicate()

        # Initialise SLOC output scanning
        sloccount   = 0
        cppcount    = 0
        pythoncount = 0
        scan      = False

        # Loop over SLOC output
        for info in sloc.split('\n'):

            # Start scanning
            if 'Totals grouped by language (dominant language first):' in info:
                scan = True
                continue

            # Scan information
            if scan:

                # Stop scanning
                if info == '':
                    scan = False
                    continue

                # Extract information
                else:
                    language, loc, _ = gammalib.split(info, " ")
                    sloccount += int(loc)
                    if 'cpp' in language:
                        cppcount += int(loc)
                    if 'python' in language:
                        pythoncount += int(loc)

        # Show results
        print(datestr, sloccount, cppcount, pythoncount)

        # Write result to ASCII file
        result.write('%s %6d %6d %6d\n' % (datestr, sloccount, cppcount, pythoncount))

    # Close result file
    result.close()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Generate sloccount history
    generate_history()
