/***************************************************************************
 *                     GCOMTools.cpp - COMPTEL tools                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCOMTools.hpp
 * @brief Implementation of COMPTEL tools
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTime.hpp"
#include "GCOMTools.hpp"

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/***********************************************************************//**
 * @brief Convert TJD and COMPTEL tics in GTime object
 *
 * @param[in] tjd Truncated Julian Days (days).
 * @param[in] tics COMPTEL tics (1/8 ms).
 * @return Time.
 *
 * Converts TJD and COMPTEL tics into a GTime object. COMPTEL times are
 * given in UTC, i.e. 8393:0 converts into 1991-05-17T00:00:00 UT
 * (see COM-RP-UNH-DRG-037).
 ***************************************************************************/
GTime gammalib::com_time(const int& tjd, const int& tics)
{
    // Compute MJD
    double mjd = double(tjd) + 40000.0;

    // Set time and retrieve result in native seconds
    GTime time;
    time.mjd(mjd, "UTC");
    double secs = time.secs();

    // Add tics and set time in native seconds
    secs += double(tics) * 0.000125;

    // Set time
    time.secs(secs);

    // Return time
    return time;
}


/***********************************************************************//**
 * @brief Convert GTime in COMPTEL TJD
 *
 * @param[in] time Time.
 * @return Truncated Julian Days (days).
 *
 * Converts GTime object in TJD.
 ***************************************************************************/
int gammalib::com_tjd(const GTime& time)
{
    // Compute TJD
    int tjd = int(time.mjd("UTC") - 40000.0);

    // Return TJD
    return tjd;
}


/***********************************************************************//**
 * @brief Convert GTime in COMPTEL tics
 *
 * @param[in] time Time.
 * @return COMPTEL tics (1/8 ms).
 *
 * Converts GTime object in COMPTEL tics (1/8 ms).
 ***************************************************************************/
int gammalib::com_tics(const GTime& time)
{
    // Compute COMPTEL time at 0 tics
    GTime tjd = com_time(com_tjd(time), 0);

    // Compute time difference in seconds
    int tics = int((time - tjd) * 8000.0 + 0.5); // rounding to nearest int

    // Return tics
    return tics;
}
