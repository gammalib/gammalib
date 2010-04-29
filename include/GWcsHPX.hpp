/***************************************************************************
 *                 GWcsHPX.hpp  -  Healpix projection class                *
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
 * @file GWcsHPX.hpp
 * @brief GWcsHPX class definition.
 * @author J. Knodlseder
 */

#ifndef GWCSHPX_HPP
#define GWCSHPX_HPP

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GWcs.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"


/***********************************************************************//**
 * @class GWcsHPX
 *
 * @brief GWcsHPX class interface defintion
 ***************************************************************************/
class GWcsHPX : public GWcs {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GWcsHPX& wcs);

public:
    // Constructors and destructors
    GWcsHPX(void);
    GWcsHPX(const int& nside, const std::string& ordering = "NESTED",
            const std::string& coordsys = "GAL");
    explicit GWcsHPX(const GFitsHDU* hdu);
    GWcsHPX(const GWcsHPX& wcs);
    virtual ~GWcsHPX(void);

    // Operators
    GWcsHPX& operator= (const GWcsHPX& wcs);

    // Implemented pure virtual methods
    std::string type(void) const;
    void        read(const GFitsHDU* hdu);
    void        write(GFitsHDU* hdu) const;

    // Overloaded base class methods
    double      omega(const int& pix) const;
    GSkyDir     pix2dir(const int& pix);
    int         dir2pix(GSkyDir dir) const;

    // Class specific methods
    int         npix(void) const;
    int         nside(void) const;
    std::string ordering(void) const;
    void        ordering(const std::string& ordering);

private:
    // Private methods
    void         init_members(void);
    void         copy_members(const GWcsHPX& wcs);
    void         free_members(void);
    GWcsHPX*     clone(void) const;
    int          nside2order(int nside);
    void         pix2xy(const int& ipix, int* x, int* y);
    int          xy2pix(int x, int y) const;
    void         pix2ang_ring(int ipix, double* theta, double* phi);
    void         pix2ang_nest(int ipix, double* theta, double* phi);
    int          ang2pix_z_phi_ring(double z, double phi) const;
    int          ang2pix_z_phi_nest(double z, double phi) const;
    unsigned int isqrt(unsigned int arg);

    // Private data area
    int      m_nside;        //!< Number of divisions of each base pixel (1-8192)
    int      m_npface;       //!<
    int      m_ncap;         //!<
    int      m_order;        //!< Order
    int      m_ordering;     //!< Pixel ordering (0=ring, 1=nested, -1=?)
    int      m_num_pixels;   //!< Number of pixels
    double   m_fact1;        //!<
    double   m_fact2;        //!<
    double   m_omega;        //!< Solid angle of pixel
};

#endif /* GWCSHPX_HPP */
