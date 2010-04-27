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
#include "GSkyDir.hpp"
#include "GWcs.hpp"


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
    GWcsCAR(GSkyDir& crval, const double& crpix1, const double& crpix2,
            const double& cdelt1, const double& cdelt2,
            const std::string& coords);
    GWcsCAR(const GFitsHDU* hdu);
    GWcsCAR(const GWcsCAR& wcs);
    virtual ~GWcsCAR(void);

    // Operators
    GWcsCAR& operator= (const GWcsCAR& wcs);

    // Implemented virtual methods
    std::string type(void) const;
    void        read(const GFitsHDU* hdu);
    void        write(GFitsHDU* hdu);
    GSkyDir     pix2dir(const int& ipix);
    int         dir2pix(GSkyDir dir) const;
    double      omega(const int& pix) const;

    // Class specific methods

private:
    // Private methods
    void     init_members(void);
    void     copy_members(const GWcsCAR& wcs);
    void     free_members(void);
    GWcsCAR* clone(void) const;

    // Private data area
};

#endif /* GWCSCAR_HPP */
