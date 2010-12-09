/***************************************************************************
 *                GWcsCAR.hpp  -  Cartesian projection class               *
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
#include "GFitsHDU.hpp"
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
    explicit GWcsCAR(const std::string& coords,
                     const double& crval1, const double& crval2,
                     const double& crpix1, const double& crpix2,
                     const double& cdelt1, const double& cdelt2,
                     const GMatrix& cd, const GVector& pv2);
    explicit GWcsCAR(const GFitsHDU* hdu);
    GWcsCAR(const GWcsCAR& wcs);
    virtual ~GWcsCAR(void);

    // Operators
    GWcsCAR& operator= (const GWcsCAR& wcs);

    // Implemented pure virtual base class methods
    void        clear(void);
    GWcsCAR*    clone(void) const;
    void        read(const GFitsHDU* hdu);
    void        write(GFitsHDU* hdu) const;
    std::string print(void) const;

    // Overloaded base class methods
    double omega(const GSkyPixel& pix) const;

private:
    // Private methods
    void     init_members(void);
    void     copy_members(const GWcsCAR& wcs);
    void     free_members(void);
    void     std2nat(GVector *coord) const;
    void     nat2std(GVector *coord) const;
};

#endif /* GWCSCAR_HPP */
