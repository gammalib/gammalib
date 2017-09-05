/***************************************************************************
 *                GCOMSupport.cpp - COMPTEL support functions              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
 * @file GCOMSupport.hpp
 * @brief Implementation of support function used by COMPTEL classes
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GWcs.hpp"
#include "GWcsCAR.hpp"
#include "GTime.hpp"
#include "GSkyMap.hpp"
#include "GCOMSupport.hpp"

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/***********************************************************************//**
 * @brief Changes Mercator's projection to cartesian projection
 *
 * @param[in,out] map Skymap.
 *
 * Changes the World Coordinate System of a sky map from Mercartor's
 * projection to a cartesian projection. This transformation is needed to
 * correct the wrong WCS headers that are found in the COMPTEL files
 * distributed through the HEASARC web site.
 *
 * The method operates only on sky maps that use Mercator's projection.
 * Nothing is done for all other sky map. This method can thus be
 * transparently be applied to the sky maps read from COMPTEL data files.
 ***************************************************************************/
void com_wcs_mer2car(GSkyMap& map)
{
    // Get WCS poiunter
    const GWcs* wcs = dynamic_cast<const GWcs*>(map.projection());

    // Apply kluge only if the skymap is in Mercator's projection
    if (wcs != NULL && wcs->code() == "MER") {

        // Extract original WCS definition
        double crval1 = wcs->crval(0);
        double crval2 = wcs->crval(1);
        double crpix1 = wcs->crpix(0);
        double crpix2 = wcs->crpix(1);
        double cdelt1 = wcs->cdelt(0);
        double cdelt2 = wcs->cdelt(1);
        int    naxis1 = map.nx();
        int    naxis2 = map.ny();

        // Compute adjustment of WCS definition to map centre
        double adjust1 = (double(naxis1) + 1.0) / 2.0 - crpix1;
        double adjust2 = (double(naxis2) + 1.0) / 2.0 - crpix2;

        // Adjust WCS definition to map centre
        crval1 += adjust1;
        crpix1 += adjust1;
        crval2 += adjust2;
        crpix2 += adjust2;

        // Allocate WCS
        GWcsCAR car(wcs->coordsys(), crval1, crval2,
                                     crpix1, crpix2,
                                     cdelt1, cdelt2);

        // Set new WCS
        map.projection(car);

    } // endif: skymap was in Mercator projection

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return D1 energy deposit
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return D1 energy deposit in MeV.
 ***************************************************************************/
double com_energy1(const double& energy, const double& phigeo)
{
    // Compute D2 energy deposit
    double e2 = com_energy2(energy, phigeo);

    // Return D1 energy deposit
    return (energy - e2);
}


/***********************************************************************//**
 * @brief Return D2 energy deposit
 *
 * @param[in] energy Input photon energy (MeV).
 * @param[in] phigeo Geometrical scatter angle (deg).
 * @return D2 energy deposit in MeV.
 ***************************************************************************/
double com_energy2(const double& energy, const double& phigeo)
{
    // Compute 1-cos(phigeo)
    double one_minus_cos  = 1.0 - std::cos(phigeo * gammalib::deg2rad);

    // Compute D2 energy deposit
    double e2 = energy / (one_minus_cos * energy / gammalib::mec2 + 1.0);

    // Return D2 energy deposit
    return e2;
}


/***********************************************************************//**
 * @brief Convert TJD and COMPTEL ticks in GTime object
 *
 * @param[in] tjd Truncated Julian Days (days).
 * @param[in] tics COMPTEL ticks (1/8 ms).
 * @return Time.
 *
 * Converts TJD and COMPTEL ticks into a GTime object. COMPTEL times are
 * given in UTC, i.e. 8393:0 converts into 1991-05-17T00:00:00 UT
 * (see COM-RP-UNH-DRG-037).
 ***************************************************************************/
GTime com_time(const int& tjd, const int& tics)
{
    // Compute MJD
    double mjd = double(tjd) + 40000.0;

    // Set time and retrieve result in native seconds
    GTime time;
    time.mjd(mjd, "UTC");
    double secs = time.secs();

    // Add ticks and set time in native seconds
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
int com_tjd(const GTime& time)
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
 * @return COMPTEL ticks (1/8 ms).
 *
 * Converts GTime object in COMPTEL ticks (1/8 ms).
 ***************************************************************************/
int com_tics(const GTime& time)
{
    // Compute COMPTEL time at 0 tics
    GTime tjd = com_time(com_tjd(time), 0);

    // Compute time difference in seconds
    int tics = int((time - tjd) * 8000.0 + 0.5); // rounding to nearest int

    // Return tics
    return tics;
}
