/***************************************************************************
 *                  GHealpix.hpp - Healpix projection class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2015 by Juergen Knoedlseder                         *
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
#include "GBilinear.hpp"


/***********************************************************************//**
 * @class GHealpix
 *
 * @brief HealPix projection class interface definition
 *
 * The HealPix projection class has been implemented by adapting code from
 * the HealPix library (version 2.1). For more information about HEALPix, see
 * http://healpix.jpl.nasa.gov
 ***************************************************************************/
class GHealpix : public GSkyProjection {

public:
    // Constructors and destructors
    GHealpix(void);
    GHealpix(const int&         nside,
             const std::string& ordering = "NESTED",
             const std::string& coordsys = "EQU");
    explicit GHealpix(const GFitsHDU& hdu);
    GHealpix(const GHealpix& wcs);
    virtual ~GHealpix(void);

    // Operators
    GHealpix& operator=(const GHealpix& wcs);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GHealpix*   clone(void) const;
    virtual std::string classname(void) const;
    virtual int         size(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
    virtual void        read(const GFitsHDU& hdu);
    virtual void        write(GFitsHDU& hdu) const;
    virtual double      solidangle(const GSkyPixel& pixel) const;
    virtual GSkyDir     pix2dir(const GSkyPixel& pixel) const;
    virtual GSkyPixel   dir2pix(const GSkyDir& dir) const;
    virtual GBilinear   interpolator(const GSkyDir& dir) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const int&           npix(void) const;
    const int&           nside(void) const;
    std::string          ordering(void) const;
    void                 ordering(const std::string& ordering);
    std::vector<int>     neighbours(const GSkyPixel& pixel) const;
    std::vector<GSkyDir> boundaries(const GSkyPixel& pixel, const int& step = 1) const;
    double               max_pixrad(void) const;

private:
    // Private methods
    void             init_members(void);
    void             copy_members(const GHealpix& wcs);
    void             free_members(void);
    virtual bool     compare(const GSkyProjection& proj) const;

    // Low-level HealPix methods
    int              compress_bits(const int& value) const;
    int              spread_bits(const int& value) const;
    void             pix2xyf(const int& pix, int* ix, int* iy, int* face) const;
    void             nest2xyf(const int& pix, int* ix, int* iy, int* face) const;
    void             ring2xyf(const int& pix, int* ix, int* iy, int* face) const;
    int              xyf2nest(const int& ix, const int& iy, const int& face) const;
    int              xyf2ring(const int& ix, const int& iy, const int& face) const;
    void             xyf2loc(const double& x, const double& y, const int& face,
                             double* z, double* phi) const;
    int              nside2order(const int& nside) const;
    int              nest2ring(const int& pix) const;
    int              ring2nest(const int& pix) const;
    void             pix2xy(const int& ipix, int* x, int* y) const;
    int              xy2pix(int x, int y) const;
    void             pix2ang_ring(int ipix, double* theta, double* phi) const;
    void             pix2ang_nest(int ipix, double* theta, double* phi) const;
    int              ang2pix_z_phi_ring(double z, double phi) const;
    int              ang2pix_z_phi_nest(double z, double phi) const;
    int              ring_above(const double& z) const;
    void             get_ring_info(const int& ring, int* startpix, int* ringpix,
                                   bool* shifted) const;
    void             get_ring_info(const int& ring, int* startpix, int* ringpix,
                                   double* theta, bool* shifted) const;
    GBilinear        interpolator(const double& theta, const double& phi) const;
    unsigned int     isqrt(unsigned int arg) const;
    GVector          set_z_phi(const double& z, const double& phi) const;
    GSkyDir          loc2dir(const double& z, const double& phi) const;

    // Private data area
    int      m_nside;        //!< Number of divisions of each base pixel (1-8192)
    int      m_npface;       //!<
    int      m_ncap;         //!< Number of pixels in polar cap
    int      m_order;        //!< Order
    int      m_ordering;     //!< Pixel ordering (0=ring, 1=nested, -1=?)
    int      m_num_pixels;   //!< Number of pixels in projection
    double   m_fact1;        //!<
    double   m_fact2;        //!<
    double   m_solidangle;   //!< Solid angle of pixel
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GHealpix").
 ***************************************************************************/
inline
std::string GHealpix::classname(void) const
{
    return ("GHealpix");
}


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


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pixel Sky pixel.
 * @return Solid angle of pixel.
 *
 * Returns the solid angle of the specified @p pixel. Note that HEALPix
 * pixels have all the same solid angle, hence the @p pixel argument is in
 * fact not used by the method.
 ***************************************************************************/
inline
double GHealpix::solidangle(const GSkyPixel& pixel) const
{
    return m_solidangle;
}

#endif /* GHEALPIX_HPP */
