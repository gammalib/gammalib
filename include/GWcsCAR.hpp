/***************************************************************************
 *                 GWcsCAR.hpp  -  Healpix projection class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GWcsCAR.hpp
 * @brief GWcsCAR class definition.
 * @author J. Knodlseder
 */

#ifndef GWCSCAR_HPP
#define GWCSCAR_HPP

/* __ Includes ___________________________________________________________ */
#include "GWcs.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"


/***********************************************************************//**
 * @class GWcsCAR
 *
 * @brief GWcsCAR class interface defintion
 ***************************************************************************/
class GWcsCAR : public GWcs {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GWcsCAR& wcs);

public:
    // Constructors and destructors
    GWcsCAR(void);
    GWcsCAR(const std::string& coords,
            const double& crval1, const double& crval2,
            const double& crpix1, const double& crpix2,
            const double& cdelt1, const double& cdelt2,
            const GMatrix& cd, const GVector& pv2);
    GWcsCAR(const GFitsHDU* hdu);
    GWcsCAR(const GWcsCAR& wcs);
    virtual ~GWcsCAR(void);

    // Operators
    GWcsCAR& operator= (const GWcsCAR& wcs);

    // Implemented virtual methods
    std::string type(void) const;
    GSkyDir     pix2dir(const int& pix);
    GSkyDir     xy2dir(const GSkyPixel& pix);
    int         dir2pix(GSkyDir dir) const;
    GSkyPixel   dir2xy(GSkyDir dir) const;
    double      omega(const int& pix) const;
    double      omega(const GSkyPixel& pix) const;

private:
    // Private methods
    void     init_members(void);
    void     copy_members(const GWcsCAR& wcs);
    void     free_members(void);
    GWcsCAR* clone(void) const;

    // Implemented pure virtual private methods
    void wcsxy2sph(const double& x, const double& y, double* lon, double* lat) const;
    void wcssph2xy(const double& lon, const double& lat, double* x, double* y) const;

    // Private data area
};

#endif /* GWCSCAR_HPP */
