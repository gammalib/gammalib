/***************************************************************************
 *                   GSPITools.cpp - INTEGRAL/SPI tools                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPITools.cpp
 * @brief Implementation of INTEGRAL/SPI tools
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSPITools.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/***********************************************************************//**
 * @brief Return FITS table
 *
 * @param[in] fits FITS file.
 * @param[in] extname Extension name.
 * @param[in] extver Extension version.
 * @return Pointer to FITS table.
 *
 * Returns the HDU with a specific extension name and version from the FITS
 * file. The method searched in the FITS file as well as grouping tables and
 * nested grouping tables.
 *
 * The method allocates a copy of the HDU, hence the client needs to
 * deallocate the HDU after usage.
 ***************************************************************************/
const GFitsTable* gammalib::spi_hdu(const GFits&       fits,
                                    const std::string& extname,
                                    const int&         extver)
{
    // Initialise HDU with NULL pointer
    const GFitsTable* hdu = NULL;

    // Loop over all extensions
    for (int extno = 0; extno < fits.size(); ++extno) {

        // Get FITS table extension
        const GFitsTable* ext = dynamic_cast<const GFitsTable*>(fits[extno]);

        // If extension is no table then skip extension
        if (ext == NULL) {
            continue;
        }

        // If HDU has the requested extension name and the requested extension
        // version then clone the HDU and break. If no extension version is
        // present in the file then only check on extension name.
        if (ext->extname() == extname) {
            if (ext->has_card("EXTVER")) {
                if (ext->integer("EXTVER") == extver) {
                    hdu = ext->clone();
                    break;
                }
            }
            else {
                hdu = ext->clone();
                break;
            }
        }

        // If HDU is a grouping table then search the requested extension
        // in the grouping table or any grouping tables that are found in
        // the grouping table
        else if (ext->extname() == "GROUPING") {

            // Loop over grouping table and search for a member with the
            // requested extension name and version. In case of success
            // get the URL, open the FITS file and search for the extension
            // with the requested name and version.
            for (int i = 0; i < ext->nrows(); ++i) {
                if (((*ext)["MEMBER_NAME"]->string(i)     == extname) &&
                    ((*ext)["MEMBER_VERSION"]->integer(i) == extver)) {
                    std::string filename = (*ext)["MEMBER_LOCATION"]->string(i);
                    if (!filename.empty()) {
                        std::string filepath = fits.filename().path();
                        GFits       fits_next(filepath+filename);
                        hdu = spi_hdu(fits_next, extname, extver);
                        if (hdu != NULL) {
                            break;
                        }
                    }
                }
            } // endfor: looped over grouping table

            // If we have still no HDU then search all grouping tables
            // that are found in the grouping table
            if (hdu == NULL) {
                for (int i = 0; i < ext->nrows(); ++i) {
                    if ((*ext)["MEMBER_NAME"]->string(i) == "GROUPING") {
                        std::string filename = (*ext)["MEMBER_LOCATION"]->string(i);
                        if (!filename.empty()) {
                            std::string filepath = fits.filename().path();
                            GFits       fits_next(filepath+filename);
                            hdu = spi_hdu(fits_next, extname, extver);
                            if (hdu != NULL) {
                                break;
                            }
                        }
                    }
                }
            } // endif: we had no HDU

            // If we have an HDU then break
            if (hdu != NULL) {
                break;
            }

        } // endelse: HDU was a grouping table

    } // endfor: looped over all extensions

    // Return pointer to FITS table
    return hdu;
}


/***********************************************************************//**
 * @brief Return number of HDU versions.
 *
 * @param[in] fits FITS file.
 * @param[in] extname Extension name.
 * @return Number of HDU versions.
 *
 * Returns the number of HDU versions with a specific extension name in the
 * FITS file and associated or nested grouping tables.
 ***************************************************************************/
int gammalib::spi_num_hdus(const GFits& fits, const std::string& extname)
{
    // Initialise number of versions
    int extvers = 0;

    // Loop over all extensions
    for (int extno = 0; extno < fits.size(); ++extno) {

        // Get HDU
        const GFitsHDU* ext = fits[extno];

        // If HDU has the requested extension name then increase version
        // counter
        if (ext->extname() == extname) {
            extvers++;
        }

        // If HDU is a grouping table then count the number of versions in
        // the grouping table
        else if (ext->extname() == "GROUPING") {

            // Cast extension to FITS table
            const GFitsTable* table = static_cast<const GFitsTable*>(ext);

            // Loop over grouping table and search for a member with the
            // requested extension name and version. In case of success
            // get the URL, open the FITS file and search for the extension
            // with the requested name and version.
            for (int i = 0; i < table->nrows(); ++i) {
                if ((*table)["MEMBER_NAME"]->string(i) == extname) {
                    extvers++;
                }
                else if ((*table)["MEMBER_NAME"]->string(i) == "GROUPING") {
                    std::string filename = (*table)["MEMBER_LOCATION"]->string(i);
                    if (!filename.empty()) {
                        std::string filepath = fits.filename().path();
                        GFits       fits_next(filepath+filename);
                        extvers += spi_num_hdus(fits_next, extname);
                    }
                }
            } // endfor: looped over grouping table

        } // endelse: HDU was a grouping table

    } // endfor: looped over all extensions

    // Return number of versions
    return extvers;
}


/***********************************************************************//**
 * @brief Convert IJD to GTime
 *
 * @param[in] ijd INTEGRAL Julian Days (days).
 * @return time.
 *
 * Converts time given in INTEGRAL Julian Days into a GTime object.
 ***************************************************************************/
GTime gammalib::spi_ijd2time(const double& ijd)
{
    // Convert IJD to MJD
    double mjd = ijd + 51544.0;

    // Set GTime
    GTime time;
    time.mjd(mjd);

    // Return time
    return time;
}
