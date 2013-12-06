/***************************************************************************
 *                  GHealpix.hpp - Healpix projection class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GHealpix.hpp
 * @brief HealPix projection class definition
 * @author Juergen Knoedlseder
 */

#ifndef GHEALPIX_HPP
#define GHEALPIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GSkyProjection.hpp"
#include "GFitsHDU.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"


/***********************************************************************//**
 * @class GHealpix
 *
 * @brief HealPix projection class interface defintion
 *
 * The HealPix projection class has been implemented by adapting code from
 * the HealPix library (version 2.1). For more information about HEALPix, see
 * http://healpix.jpl.nasa.gov
 ***************************************************************************/
class GHealpix : public GSkyProjection {

public:
    // Constructors and destructors
    GHealpix(void);
    explicit GHealpix(const int&         nside,
                      const std::string& ordering = "NESTED",
                      const std::string& coordsys = "GAL");
    explicit GHealpix(const GFitsHDU& hdu);
    GHealpix(const GHealpix& wcs);
    virtual ~GHealpix(void);

    // Operators
    GHealpix& operator=(const GHealpix& wcs);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GHealpix*   clone(void) const;
    virtual int         size(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
    virtual void        read(const GFitsHDU& hdu);
    virtual void        write(GFitsHDU& hdu) const;
    virtual double      omega(const GSkyPixel& pixel) const;
    virtual GSkyDir     pix2dir(const GSkyPixel& pixel) const;
    virtual GSkyPixel   dir2pix(const GSkyDir& dir) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const int&   npix(void) const;
    const int&   nside(void) const;
    std::string  ordering(void) const;
    void         ordering(const std::string& ordering);

private:
    // Private methods
    void         init_members(void);
    void         copy_members(const GHealpix& wcs);
    void         free_members(void);
    virtual bool compare(const GSkyProjection& proj) const;
    int          nside2order(int nside);
    void         pix2xy(const int& ipix, int* x, int* y) const;
    int          xy2pix(int x, int y) const;
    void         pix2ang_ring(int ipix, double* theta, double* phi) const;
    void         pix2ang_nest(int ipix, double* theta, double* phi) const;
    int          ang2pix_z_phi_ring(double z, double phi) const;
    int          ang2pix_z_phi_nest(double z, double phi) const;
    unsigned int isqrt(unsigned int arg) const;

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


/***********************************************************************//**
 * @brief Return dimension of projection
 *
 * @return Dimension of projection.
 *
 * Returns the dimension of the projection.
 ***************************************************************************/
inline
int GHealpix::size(void) const
{
    return 1;
}


/***********************************************************************//**
 * @brief Return projection code
 *
 * @return Projection code.
 *
 * Returns the projection code "HPX".
 ***************************************************************************/
inline
std::string GHealpix::code(void) const
{
    return "HPX";
}


/***********************************************************************//**
 * @brief Return projection name
 *
 * @return Projection name.
 *
 * Returns the projection name.
 ***************************************************************************/
inline
std::string GHealpix::name(void) const
{
    return "HealPix";
}


/***********************************************************************//**
 * @brief Returns number of pixels
 ***************************************************************************/
inline
const int& GHealpix::npix(void) const
{
    // Return number of pixels
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Returns number of divisions of the side of each base pixel.
 ***************************************************************************/
inline
const int& GHealpix::nside(void) const
{
    return m_nside;
}

#endif /* GHEALPIX_HPP */
