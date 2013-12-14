/***************************************************************************
 *          GWcs.hpp - Abstract world coordinate system base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GWcs.hpp
 * @brief Abstract world coordinate system base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GWCS_HPP
#define GWCS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <iostream>
#include "GSkyProjection.hpp"
#include "GLog.hpp"
#include "GFitsHDU.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"


/***********************************************************************//**
 * @class GWcs
 *
 * @brief Abstract world coordinate system base class
 ***************************************************************************/
class GWcs : public GSkyProjection {

public:
    // Constructors and destructors
    GWcs(void);
    explicit GWcs(const std::string& coords,
                  const double& crval1, const double& crval2,
                  const double& crpix1, const double& crpix2,
                  const double& cdelt1, const double& cdelt2);
    explicit GWcs(const GFitsHDU& hdu);
    GWcs(const GWcs& wcs);
    virtual ~GWcs(void);

    // Operators
    virtual GWcs& operator= (const GWcs& wcs);

    // Pure virtual base class methods
    virtual void        clear(void) = 0;
    virtual GWcs*       clone(void) const = 0;
    virtual std::string code(void) const = 0;
    virtual std::string name(void) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;
    
    // Implemented virtual base class methods
    virtual int         size(void) const;
    virtual void        read(const GFitsHDU& hdu);
    virtual void        write(GFitsHDU& hdu) const;
    virtual double      solidangle(const GSkyPixel& pixel) const;
    virtual GSkyDir     pix2dir(const GSkyPixel& pixel) const;
    virtual GSkyPixel   dir2pix(const GSkyDir& dir) const;

    // Other methods
    void   set(const std::string& coords,
               const double& crval1, const double& crval2,
               const double& crpix1, const double& crpix2,
               const double& cdelt1, const double& cdelt2);
    double crval(const int& inx) const;
    double crpix(const int& inx) const;
    double cdelt(const int& inx) const;

private:
    // Static constants (set in GWcs.cpp)
    static const int    PVN = 32;
    static const double UNDEFINED;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GWcs& wcs);
    void         free_members(void);
    void         set_members(const std::string& coords,
                             const double& crval1, const double& crval2,
                             const double& crpix1, const double& crpix2,
                             const double& cdelt1, const double& cdelt2);
    virtual bool compare(const GSkyProjection& proj) const;
    bool         undefined(const double& value) const { return (value==UNDEFINED); }

    // Methods adapted from wcslib::wcs.c 
    void        wcs_ini(int naxis);
    void        wcs_set(void) const;
    void        wcs_set_ctype(void) const;
    void        wcs_p2s(int ncoord, int nelem, const double* pixcrd, double* imgcrd,
                        double* phi, double* theta, double* world, int* stat) const;
    void        wcs_s2p(int ncoord, int nelem, const double* world,
                        double* phi, double* theta,  double* imgcrd,
                        double* pixcrd, int* stat) const;
    std::string wcs_print(const GChatter& chatter = NORMAL) const;
    std::string wcs_print_value(const double& value) const;
    
    // Methods adapted from wcslib::cel.c
    void cel_ini(void) const;
    void cel_set(void) const;
    void cel_x2s(int nx, int ny, int sxy, int sll,
                 const double* x, const double* y,
                 double* phi, double* theta,
                 double* lng, double* lat, int* stat) const;
    void cel_s2x(int nlng, int nlat, int sll, int sxy,
                 const double* lng, const double* lat,
                 double* phi, double* theta,
                 double* x, double* y, int* stat) const;

    // Methods adapted from wcslib::sph.c
    void sph_x2s(int nphi, int ntheta, int spt, int sll,
                 const double* phi, const double* theta,
                 double* lng, double* lat) const;
    void sph_s2x(int nlng, int nlat, int sll, int spt,
                 const double* lng, const double* lat,
                 double* phi, double* theta) const;

    // Methods adapted from wcslib::spc.c
    void spc_ini(void);

    // Methods adapted from wcslib::lin.c
    void lin_ini(int naxis);
    void lin_set(void) const;
    void lin_p2x(int ncoord, int nelem, const double* pixcrd, double* imgcrd) const;
    void lin_x2p(int ncoord, int nelem, const double* imgcrd, double* pixcrd) const;
    void lin_matinv(const std::vector<double>& mat, std::vector<double>& inv) const;

    // Methods adapted from wcslib::prj.c
    void         prj_ini(void) const;
    void         prj_off(const double& phi0, const double& theta0) const;
    virtual void prj_set(void) const = 0;
    virtual void prj_x2s(int nx, int ny, int sxy, int spt, 
                         const double* x, const double* y,
                         double* phi, double* theta, int* stat) const = 0;
    virtual void prj_s2x(int nphi, int ntheta, int spt, int sxy,
                         const double* phi, const double* theta,
                         double* x, double* y, int* stat) const = 0;
    
    // World Coordinate System parameters
    mutable bool                     m_wcsset;  //!< WCS information is set
    int                              m_naxis;   //!< Number of axes
    std::vector<double>              m_crval;   //!< CRVALia keyvalues for each coord axis
    std::vector<std::string>         m_cunit;   //!< CUNITia keyvalues for each coord axis
    mutable std::vector<std::string> m_ctype;   //!< CTYPEia keyvalues for each coord axis
    mutable std::vector<std::string> m_ctype_c; //!< CTYPEia comments for each coord axis
    mutable double                   m_lonpole; //!< LONPOLEa keyvalue
    mutable double                   m_latpole; //!< LATPOLEa keyvalue
    double                           m_restfrq; //!< RESTFRQa keyvalue
    double                           m_restwav; //!< RESTWAVa keyvalue
    std::string                      m_radesys; //!< RADESYS keyvalue
    double                           m_equinox; //!< EQUINOX keyvalue
    std::vector<double>              m_cd;      //!< CDi_ja linear transformation matrix
    std::vector<double>              m_crota;   //!< CROTAia keyvalues for each coord axis
    int                              m_lng;     //!< Longitude axis
    int                              m_lat;     //!< Latitude axis
    int                              m_spec;    //!< Spectral axis
    
    // Linear transformation parameters
    mutable bool                     m_linset;  //!< Linear transformation is set
    mutable bool                     m_unity;   //!< Signals unity PC matrix
    std::vector<double>              m_crpix;   //!< CRPIXja keyvalues for each pixel axis
    std::vector<double>              m_pc;      //!< PCi_ja  linear transformation matrix
    std::vector<double>              m_cdelt;   //!< CDELTia keyvalues for each coord axis
    mutable std::vector<double>      m_piximg;  //!< Pixel to image transformation matrix
    mutable std::vector<double>      m_imgpix;  //!< Image to pixel transformation matrix

    // Celestial transformation parameters
    mutable bool                     m_celset;  //!< Celestial transformation is set
    mutable bool                     m_offset;  //!< Force (x,y) = (0,0) at (phi_0,theta_0)
    mutable double                   m_phi0;    //!< Native azimuth angle of fiducial point
    mutable double                   m_theta0;  //!< Native zenith angle of fiducial point
    mutable double                   m_ref[4];  //!< Celestial coordinates of fiducial
                                                //   point and native coordinates of celestial
                                                //   pole
    mutable double                   m_euler[5];//!< Euler angles and functions thereof
    mutable int                      m_latpreq; //!< LATPOLEa requirement
    mutable bool                     m_isolat;  //!< True if |latitude| is preserved
    
    // Projection parameters
    mutable bool                     m_prjset;  //!< Projection is set
    mutable double                   m_r0;      //!< Radius of the generating sphere
    mutable double                   m_pv[PVN]; //!< Projection parameters
    mutable bool                     m_bounds;  //!< Enable strict bounds checking
    mutable double                   m_x0;      //!< Fiducial x offset
    mutable double                   m_y0;      //!< Fiducial y offset
    mutable std::vector<double>      m_w;       //!< Intermediate values
    
    // Spectral transformation parameters
    //struct spcprm spc
};


/***********************************************************************//**
 * @brief Return dimension of projection
 *
 * @return Dimension of projection.
 *
 * Returns the dimension of the projection.
 ***************************************************************************/
inline
int GWcs::size(void) const
{
    return 2;
}

#endif /* GWCS_HPP */
