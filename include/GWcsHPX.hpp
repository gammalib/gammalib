/***************************************************************************
 *                 GWcsHPX.hpp  -  Healpix projection class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GWcsHPX.hpp
 * @brief GWcsHPX class definition.
 * @author J. Knodlseder
 */

#ifndef GWCSHPX_HPP
#define GWCSHPX_HPP

/* __ Includes ___________________________________________________________ */
#include "GWcs.hpp"
#include "GFitsHDU.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"


/***********************************************************************//**
 * @class GWcsHPX
 *
 * @brief GWcsHPX class interface defintion
 ***************************************************************************/
class GWcsHPX : public GWcs {

public:
    // Constructors and destructors
    GWcsHPX(void);
    explicit GWcsHPX(const int& nside, const std::string& ordering = "NESTED",
                     const std::string& coordsys = "GAL");
    explicit GWcsHPX(const GFitsHDU* hdu);
    GWcsHPX(const GWcsHPX& wcs);
    virtual ~GWcsHPX(void);

    // Operators
    GWcsHPX& operator= (const GWcsHPX& wcs);

    // Implemented pure virtual base class methods
    void        clear(void);
    GWcsHPX*    clone(void) const;
    void        read(const GFitsHDU* hdu);
    void        write(GFitsHDU* hdu) const;
    std::string print(void) const;

    // Overloaded base class methods
    double      omega(const int& pix) const;
    GSkyDir     pix2dir(const int& pix) const;
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
    void         std2nat(GVector *coord) const { return; }
    void         nat2std(GVector *coord) const { return; }
    int          nside2order(int nside);
    void         pix2xy(const int& ipix, int* x, int* y) const;
    int          xy2pix(int x, int y) const;
    void         pix2ang_ring(int ipix, double* theta, double* phi) const;
    void         pix2ang_nest(int ipix, double* theta, double* phi) const;
    int          ang2pix_z_phi_ring(double z, double phi) const;
    int          ang2pix_z_phi_nest(double z, double phi) const;
    unsigned int isqrt(unsigned int arg) const;

    // NEW VERSION
    void prj_set(void);
    int  prj_x2s(int nx, int ny, int sxy, int spt, 
                 const double* x, const double* y,
                 double* phi, double* theta, int* stat) const;
    int  prj_s2x(int nphi, int ntheta, int spt, int sxy,
                 const double* phi, const double* theta,
                 double* x, double* y, int* stat) const;

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
