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

            // If the extension has a GRPNAME then check if the requested
            // extension name if that grouping name
            if (ext->has_card("GRPNAME")) {
                if (ext->string("GRPNAME") == extname) {
                    hdu = ext->clone();
                    break;
                }
            }

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


/***********************************************************************//**
 * @brief Return start time of annealing operations
 *
 * @return Start time of annealing operations.
 *
 * Returns the start time of the SPI annealing operations.
 *
 * Source: https://www.cosmos.esa.int/web/integral/long-term-plan
 ***************************************************************************/
GTimes gammalib::spi_annealing_start_times(void)
{
    // Allocates times
    GTimes times;

    // Set annealing start times
    times.append(GTime("2003-02-06T09:00:00")); //  1.
    times.append(GTime("2003-07-15T12:00:00")); //  2.
    times.append(GTime("2003-11-11T12:00:00")); //  3.
    times.append(GTime("2004-06-17T09:00:00")); //  4.
    times.append(GTime("2005-01-19T09:00:00")); //  5.
    times.append(GTime("2005-06-14T00:00:00")); //  6.
    times.append(GTime("2006-01-09T00:00:00")); //  7.
    times.append(GTime("2006-06-08T00:00:00")); //  8.
    times.append(GTime("2006-12-04T10:00:00")); //  9.
    times.append(GTime("2007-05-29T12:00:00")); // 10.
    times.append(GTime("2008-01-12T04:00:00")); // 11.
    times.append(GTime("2008-08-17T04:00:00")); // 12.
    times.append(GTime("2009-04-20T08:49:33")); // 13. Rev.  796- 802
    times.append(GTime("2009-10-19T19:38:18")); // 14. Rev.  857- 862
    times.append(GTime("2010-03-30T09:52:05")); // 15. Rev.  911- 916
    times.append(GTime("2010-10-10T19:44:20")); // 16. Rev.  976- 981
    times.append(GTime("2011-04-26T06:35:08")); // 17. Rev. 1042-1047
    times.append(GTime("2011-11-21T15:36:48")); // 18. Rev. 1112-1118
    times.append(GTime("2012-06-06T01:30:06")); // 19. Rev. 1178-1183
    times.append(GTime("2013-01-04T10:32:36")); // 20. Rev. 1249-1253
    times.append(GTime("2013-08-01T21:00:00")); // 21. Rev. 1319-1324
    times.append(GTime("2014-01-07T08:32:20")); // 22. Rev. 1372-1377
    times.append(GTime("2014-08-28T16:05:05")); // 23. Rev. 1450-1453
    times.append(GTime("2015-02-15T12:54:04")); // 24. Rev. 1508-1512
    times.append(GTime("2015-09-08T09:17:16")); // 25. Rev. 1585-1590
    times.append(GTime("2016-03-04T11:23:26")); // 26. Rev. 1652-1656
    times.append(GTime("2016-07-23T11:19:51")); // 27. Rev. 1705-1710
    times.append(GTime("2017-01-14T23:26:20")); // 28. Rev. 1771-1776
    times.append(GTime("2017-07-25T09:56:01")); // 29. Rev. 1843-1848
    times.append(GTime("2018-01-27T14:29:29")); // 30. Rev. 1913-1918
    times.append(GTime("2018-07-22T00:28:13")); // 31. Rev. 1979-1984
    times.append(GTime("2019-01-21T13:31:40")); // 32. Rev. 2048-2053
    times.append(GTime("2019-09-23T01:41:56")); // 33. Rev. 2140-2145
    times.append(GTime("2020-03-05T21:56:30")); // 34. Rev. 2202-2207

    // Return times
    return times;
}


/***********************************************************************//**
 * @brief Return times of detector failures
 *
 * @return Times of detector failures.
 *
 * Returns the times of detector failures.
 ***************************************************************************/
GTimes gammalib::spi_gedfail_times(void)
{
    // Allocates times
    GTimes times;

    // Set detector failure times
    times.append(GTime("2003-12-06T09:58:30")); // Detector #2
    times.append(GTime("2004-07-17T11:00:00")); // Detector #17
    times.append(GTime("2009-02-19T12:00:00")); // Detector #5
    times.append(GTime("2010-05-27T16:00:00")); // Detector #1

    // Return times
    return times;
}
