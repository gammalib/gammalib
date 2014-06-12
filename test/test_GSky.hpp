/***************************************************************************
 *                  test_GSky.hpp - Test sky module                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Jean-Baptiste Cayrou                        *
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

#ifndef TEST_GSKY_HPP
#define TEST_GSKY_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"

/***********************************************************************//**
 * @class TestGSky
 *
 * @brief Test suite for sky module
 ***************************************************************************/
class TestGSky : public GTestSuite
{
public:
    // Constructors and destructors
    TestGSky(void) : GTestSuite() { }
    virtual ~TestGSky(void) { }

    // Methods
    virtual void      set(void);
    virtual TestGSky* clone(void) const;
    void              test_GWcs(void);
    void              test_GSkyPixel(void);
    void              test_GSkymap_healpix_construct(void);
    void              test_GSkymap_healpix_io(void);
    void              test_GSkymap_wcs_construct(void);
    void              test_GSkymap_wcs_io(void);
    void              test_GSkymap(void);
    void              test_GSkyRegions_io(void);
    void              test_GSkyRegionCircle_construct(void);
    void              test_GSkyRegionCircle_logic(void);
    void              test_GHorizDir(void);

// Private methods
private:
    double wcs_forth_back_pixel(GWcs* wcs, int nx, int ny, double& crpix1, double& crpix2);
    double wcs_copy(GWcs* wcs, int nx, int ny, double& crpix1, double& crpix2);
};

#endif /* TEST_GSKY_HPP */
