/***************************************************************************
 *                     GCOMTools.cpp - COMPTEL tools                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021-2022 by Juergen Knoedlseder                         *
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
#include "GTools.hpp"
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
 *
 * The method applies the CGRO clock correction. The CGRO clock was too
 * fast by 2.042144 seconds before 8798:28800000 (1992-06-25T01:00:00).
 * After that date the clock was corrected.
 ***************************************************************************/
GTime gammalib::com_time(const int& tjd, const int& tics)
{
    // Compute tics to days conversion factor
    const double tics2days     = 0.000125 * gammalib::sec2day;
    const double clockcor_days = 2.042144 * gammalib::sec2day;

    // Compute MJD
    double mjd = double(tjd) + 40000.0 + double(tics) * tics2days;

    // Apply CGRO clock correction as the clock was too fast by 2.042144 sec
    // before 8798:28800000 (1992-06-25T01:00:00).
    if ((tjd < 8798) || ((tjd == 8798) && (tics < 28800000))) {
        mjd -= clockcor_days;
    }

    // Set time in MJD in the UTC time system
    GTime time;
    time.mjd(mjd, "UTC");

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
 *
 * The method applies the CGRO clock correction. The CGRO clock was too
 * fast by 2.042144 seconds before 8798:28800000 (1992-06-25T01:00:00).
 * After that date the clock was corrected.
 ***************************************************************************/
int gammalib::com_tjd(const GTime& time)
{
    // Set MJD in UTC of clock correction
    const double mjd_of_clockcor = 48798.04166666666424134746;
    const double clockcor_days   = 2.042144 * gammalib::sec2day;

    // Compute MJD, applying the CGRO clock correction before 8798:28800000
    // (1992-06-25T01:00:00)
    double mjd = time.mjd("UTC");
    if (mjd < mjd_of_clockcor) {
        mjd += clockcor_days;
    }

    // Compute TJD
    int tjd = int(mjd - 40000.0);

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
 *
 * The method applies the CGRO clock correction. The CGRO clock was too
 * fast by 2.042144 seconds before 8798:28800000 (1992-06-25T01:00:00).
 * After that date the clock was corrected.
 ***************************************************************************/
int gammalib::com_tics(const GTime& time)
{
    // Set MJD in UTC of clock correction
    const double mjd_of_clockcor = 48798.04166666666424134746;
    const double clockcor_days   = 2.042144 * gammalib::sec2day;

    // Set number of tics in one day
    const double tics_in_day = 8000.0 * gammalib::sec_in_day;

    // Compute MJD, applying the CGRO clock correction before 8798:28800000
    // (1992-06-25T01:00:00)
    double mjd = time.mjd("UTC");
    if (mjd < mjd_of_clockcor) {
        mjd += clockcor_days;
    }

    // Compute number of tics, rounding to nearest int
    int tics = int((mjd - double(int(mjd))) * tics_in_day + 0.5);

    // Return tics
    return tics;
}
