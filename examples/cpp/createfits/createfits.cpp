/***************************************************************************
 *                   createfits.cpp - Create FITS file                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file createfits.cpp
 * @brief Create FITS file
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @brief Create FITS file
 *
 * This code illustrates the creation of a FITS file.
 ***************************************************************************/
int main(void) {

    // Allocate FITS object
    GFits fits;

    // Create FITS image
    GFitsImageDouble image(20,10);
    for (int x = 0; x < 20; ++x) {
        for (int y = 0; y < 10; ++y) {
            image(x,y) = x+100.0*y;
        }
    }

    // Append image to FITS object
    fits.append(image);

    // Create binary table
    GFitsBinTable       table;
    GFitsTableDoubleCol column("ENERGY", 10, 3);
    for (int row = 0; row < 10; ++row) {
        for (int index = 0; index < 3; ++index) {
            column(row, index) = row*100.0+index;
        }
    }
    table.append_column(column);

    // Append table to FITS object
    fits.append(table);

    // Save file to disk
    fits.saveto("my_fits_file.fits", true);

    // Close FITS object
    fits.close();

    // Exit
    return 0;
}

