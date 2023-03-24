#! /usr/bin/env python
# ==========================================================================
# Test BGDLIXF rate fitting
# --------------------------------------------------------------------------
#
# Copyright (C) 2023 Juergen Knoedlseder
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
import gammalib
import glob
import math
import matplotlib.pyplot as plt


# =================== #
# Get Veto dome rates #
# =================== #
def get_veto_rates(xmlname):
    """
    Get Veto dome rates

    Parameters
    ----------
    xmlname : str
        Observation definition XML file name

    Returns
    -------
    rates : dict
        Veto rates
    """
    # Initialise rates dictionary
    rates = {'times': [], 'rates': []}

    # Load observations
    obs = gammalib.GObservations(xmlname)

    # Loop over observations
    for com in obs:

        # Get identifier
        rates['id'] = com.id()

        # Get 2nd Veto dome rates housekeeping data
        hkd = com.hkds()['SCV2M']

        # Extract data
        for i in range(hkd.size()):
            time = hkd.time(i).copy()
            rate = hkd.value(i) * 12.8 / 2.048
            if rate > 0.0:
                rates['times'].append(time)
                rates['rates'].append(rate)

    # Return rates
    return rates


# ====================================== #
# Check time ordering of Veto dome rates #
# ====================================== #
def check_veto_rates():
    """
    Check time ordering of Veto dome rates
    """
    # Get all XML files
    xmlnames = glob.glob('dbase/xml/vp*.xml')
    xmlnames.sort()

    # Loop over XML files
    for xmlname in xmlnames:

        # Load observations
        obs = gammalib.GObservations(xmlname)

        # Loop over observations
        for com in obs:

            # Get 2nd Veto dome rates housekeeping data
            hkd = com.hkds()['SCV2M']

            # Check times
            nbad = 0
            for i in range(hkd.size()):
                if i == 0:
                    last_time = hkd.time(i)
                else:
                    if hkd.time(i) <= last_time:
                        nbad += 1
                        print('%s: %.5f <= %.5f' % (com.id(), hkd.time(i).mjd(), last_time.mjd()))
            if nbad == 0:
                print('%s: ok.' % (com.id()))

    # Return
    return


# ==================== #
# Show Veto dome rates #
# ==================== #
def show_veto_rates(rates):
    """
    Show Veto dome rates
    """
    # Create figure
    fig = plt.figure(figsize=(14,4))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.93, bottom=0.12)
    plt.title('Veto dome rates for %s' % rates['id'])

    # Setup vectors
    x = [time.mjd() for time in rates['times']]
    y = [rate       for rate in rates['rates']]

    # Plot vectors
    plt.plot(x, y, 'r-')

    # Set attributes
    plt.xlabel(r'MJD (days)')
    plt.ylabel(r'Rate (triggers/sec)')

    # Show plot
    plt.show()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Dump header
    print('')
    print('*****************************')
    print('* Test BGDLIXF rate fitting *')
    print('*****************************')

    # Check time ordering of veto rates
    #check_veto_rates()

    # Set file names
    xmlname = 'dbase/xml/vp0001_0.xml'

    # Get Veto dome rates
    rates = get_veto_rates(xmlname)

    # Show Veto dome rates
    show_veto_rates(rates)
