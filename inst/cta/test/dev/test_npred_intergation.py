#! /usr/bin/env python
# ==========================================================================
# This Python script verifies the Npred computation for the CTA radial
# acceptance model.
#
# Copyright (C) 2012 Jurgen Knodlseder
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
import math


# ================= #
# Compute arclength #
# ================= #
def cta_roi_arclength(rad, dist, roi):
    """
    Returns length of circular arc within circular ROI

    Parameters:
     rad  - Circle radius in radians (<pi)
     dist - Circle centre distance to ROI centre (<pi)
     roi  - Radius of ROI in radians
    """
    # Initialise
    arclength = 0.0

    # Handle special case that circle centre matches ROI centre
    if dist == 0.0:
        if rad > roi:
            arclength = 0.0
        else:
            arclength = 2.0 * math.pi

    # ... otherwise circle and ROI centres are not identical
    else:

        # Handle special case that we evaluate exactly at the circle
        # centre. In this case we have in fact a point, and if this point
        # falls within the ROI it has a formal arclength of 2pi.
        if rad == 0.0:
            if dist > roi:
                arclength = 0.0
            else:
                arclength = 2.0 * math.pi

        # ... otherwise we have to handle the general case
        else:
            d = roi - dist
            if -rad >= d:
                arclength = 0.0
            elif rad <= d:
                arclength = 2.0 * math.pi
            else:
                cosrad = math.cos(rad)
                sinrad = math.sin(rad)
                cosroi = math.cos(roi)
                cosdist = math.cos(dist)
                sindist = math.sin(dist)
                cosang = (cosroi - cosdist * cosrad) / (sindist * sinrad)
                arclength = 2.0 * math.acos(cosang)

    # Return
    return arclength


# =================================== #
# Integration kernel for Npred method #
# =================================== #
def eval(offset, dist, roi):
    """
    Returns integration kernel

    Parameters:
     offset - Offset angle in radians
     dist   - Distance of ROI and pointing in radians
     roi    - ROI radius in radians
    """
    # Initialise Npred value
    value = 0.0

    # Continue only for positive offsets
    if offset > 0:

        # Get arclength
        phi = cta_roi_arclength(offset, dist, roi)

        # Get value if phi > 0
        if phi > 0.0:
            value = phi * math.sin(offset)

    # Return
    return value


# ======================= #
# Main script entry point #
# ======================= #
if __name__ == '__main__':
    """
    Analyse data using unbinned gtlike.
    """
    # Set parameters
    d2r = math.pi / 180.0
    dist = 1.0 * d2r
    roi = 1.0 * d2r

    # Perform numerical integration
    if dist > roi:
        rmin = dist - roi
    else:
        rmin = 0.0
    rmax = roi + dist
    nr = 1000
    dr = (rmax - rmin) / nr
    sum = 0.0
    r = rmin + dr
    while r <= rmax:
        sum += eval(r, dist, roi) * dr
        r += dr

    # Perform analytical integration
    ana = math.pi * roi * roi

    # Print result
    print "Numerical integration:  ", sum
    print "Analytical integration: ", ana
