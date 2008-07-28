/***************************************************************************
 *             GHealpix.hpp  -  Healpix sky representation class           *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GHealpix.hpp
 * @brief GHealpix class definition.
 * @author J. Knodlseder
 */

#ifndef GHEALPIX_HPP
#define GHEALPIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableDblCol.hpp"
#include "GSkyDir.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GHealpix
 *
 * @brief GHealpix class interface defintion
 ***************************************************************************/
class GHealpix {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GHealpix& pixels);

public:
    // Constructors and destructors
    GHealpix(int nside, std::string scheme = "NESTED",
             std::string coordsys = "GAL", int dimension = 1);
    GHealpix(const GFitsHDU* hdu);
    GHealpix(const GHealpix& pixels);
    virtual ~GHealpix();

    // Pixel access operators
    double& operator() (int pixel, int element = 0);
    const double& operator() (int pixel, int element = 0) const;

    // Operators
    GHealpix& operator= (const GHealpix& pixels);

    // Methods
    void    read(const GFitsHDU* hdu);
    void    write(GFits* fits);
    int     nside(void) const;
    int     order(void) const;
    int     npix(void) const;
    double  omega(void) const;
    GSkyDir pix2ang(const int& ipix);
    int     ang2pix(GSkyDir dir) const;

private:
    // Private methods
    void      init_members(void);
    void      alloc_members(void);
    void      copy_members(const GHealpix& pixels);
    void      free_members(void);
    GHealpix* clone(void) const;
    int       nside2order(int nside);
    void      pix2xy(const int& ipix, int* x, int* y);
    int       xy2pix(int x, int y) const;
    void      pix2ang_ring(int ipix, double* theta, double* phi);
    void      pix2ang_nest(int ipix, double* theta, double* phi);
    int       ang2pix_z_phi_ring(double z, double phi) const;
    int       ang2pix_z_phi_nest(double z, double phi) const;

    // Private data area
    int      m_nside;        //!< Number of divisions of each base pixel (1-8192)
    int      m_npface;       //!<
    int      m_ncap;         //!<
    int      m_order;        //!< Order
    int      m_scheme;       //!< Ordering scheme (0=ring, 1=nested, -1=?)
    int      m_coordsys;     //!< Coordinate system (0=equatorial, 1=galactic, -1=?)
    int      m_num_pixels;   //!< Number of pixels
    int      m_size_pixels;  //!< Vector size of each pixel
    double   m_fact1;        //!<
    double   m_fact2;        //!<
    double   m_omega;        //!< Solid angle of pixel
    double*  m_pixels;       //!< Pixel array
    GSkyDir* m_dir;          //!< Sky direction for each pixel
};

#endif /* GHEALPIX_HPP */
