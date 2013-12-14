/***************************************************************************
 *                    createskymap.cpp - Create sky map                    *
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
 * @file createskymap.cpp
 * @brief Create sky map
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @brief Create sky map
 *
 * This code illustrates the creation of a sky map.
 ***************************************************************************/
int main(void) {

    // Create a SNR shell model. The shell is centred on RA=0.3 deg and
    // Dec=0.1 deg. The shell has an inner radius of 0.5 deg and a width
    // of 0.1 deg.
    GSkyDir centre;
    centre.radec_deg(0.3, 0.1);
    GModelSpatialRadialShell model(centre, 0.5, 0.1);

    // Create an empty sky map in celestial coordinates using a cartesian
    // projection. The image is centred on RA=Dec=0, has a bin size of
    // 0.01 degrees and 300 pixels in RA and Dec
    double ra(0.0);
    double dec(0.0);
    double binsz(0.01);
    int    npix(300);
    GSkymap image("CAR", "CEL", ra, dec, -binsz, binsz, npix, npix, 1);

    // Fill the sky map with the model image
    GEnergy energy; // Dummy (not relevant as model is not energy dependent)
    GTime   time;   // Dummy (not relevant as model is not time dependent)
    for (int index = 0; index < image.npix(); ++index) {
        GSkyDir dir   = image.inx2dir(index); // sky coordinate of pixel
        double  theta = centre.dist(dir);     // distance in radians
        image(index)  = model.eval(theta, energy, time);
    }

    // Save the image to a FITS file
    image.save("my_sky_map.fits", true);

    // Exit
    return 0;
}

