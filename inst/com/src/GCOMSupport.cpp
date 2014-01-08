/***************************************************************************
 *                GCOMSupport.cpp - COMPTEL support functions              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
#include "GCOMSupport.hpp"
#include "GWcs.hpp"
#include "GWcsCAR.hpp"

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
void com_wcs_mer2car(GSkymap& map)
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
