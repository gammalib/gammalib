/***************************************************************************
 *  GModelSpatialEllipticalGeneralGauss.hpp -                              *
 *  Elliptical generalised gauss source model class                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022- by Luigi Tibaldo                                   *
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
 * @file GModelSpatialEllipticalGeneralGauss.hpp
 * @brief Elliptical generalised gaussian model class interface definition
 * @author Luigi Tibaldo
 */

#ifndef GMODELSPATIALELLIPTICALGENERALGAUSS_HPP
#define GMODELSPATIALELLIPTICALGENERALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialElliptical.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GSkyRegionCircle.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelSpatialEllipticalGeneralGauss
 *
 * @brief Elliptical generalised gaussian source model class
 *
 * This class implements the spatial component of the factorised source
 * model for an elliptical generalised  gaussian source, i.e. surface brightness
 * according to an asymmetric generalised Gaussian.
 ***************************************************************************/
class GModelSpatialEllipticalGeneralGauss : public GModelSpatialElliptical {

public:
    // Constructors and destructors
    GModelSpatialEllipticalGeneralGauss(void);
    GModelSpatialEllipticalGeneralGauss(const bool& dummy, const std::string& type);
    GModelSpatialEllipticalGeneralGauss(const GSkyDir& dir,
					const double&  major,
					const double&  minor,
					const double&  posangle,
					const double& ridx);
    explicit GModelSpatialEllipticalGeneralGauss(const GXmlElement& xml);
    GModelSpatialEllipticalGeneralGauss(const GModelSpatialEllipticalGeneralGauss& model);
    virtual ~GModelSpatialEllipticalGeneralGauss(void);

    // Operators
    virtual GModelSpatialEllipticalGeneralGauss& operator=(const GModelSpatialEllipticalGeneralGauss& model);

    // Implemented pure virtual base class methods
    virtual void                          clear(void);
    virtual GModelSpatialEllipticalGeneralGauss* clone(void) const;
    virtual std::string                   classname(void) const;
    virtual double                        eval(const double&  theta,
                                               const double&  posangle,
                                               const GEnergy& energy,
                                               const GTime&   time,
                                               const bool&    gradients = false) const;
    virtual GSkyDir                       mc(const GEnergy& energy,
                                             const GTime& time,
                                             GRan& ran) const;
    virtual bool                          contains(const GSkyDir& dir,
                                                   const double&  margin = 0.0) const;
    virtual double                        theta_max(void) const;
    virtual void                          read(const GXmlElement& xml);
    virtual void                          write(GXmlElement& xml) const;
    virtual std::string                   print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  ridx(void) const;
    void    ridx(const double& ridx);


protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialEllipticalGeneralGauss& model);
    void         free_members(void);
    void         update(void) const;
    virtual void set_region(void) const;

    // Protected members
    GModelPar m_ridx;                   //!< reciprocal of exponent of the radial profile

    // Cached members used for pre-computations
    mutable double m_last_minor;        //!< Last semi-minor axis
    mutable double m_last_major;        //!< Last semi-major axis
    mutable double m_last_posangle;     //!< Last position angle
    mutable double m_cospos2;           //!< squared cosine of position angle
    mutable double m_sinpos2;           //!< squared sine of position angle
    mutable double m_sin2pos;           //!< sine of twice the position angle
    mutable double m_minor2;            //!< square of minor axis
    mutable double m_major2;            //!< square of major axis
    mutable double m_minor_rad;         //!< Minor axis in radians
    mutable double m_major_rad;         //!< Major axis in radians
    mutable double m_norm;              //!< Normalization
    mutable double m_term1;             //!< Help term 1
    mutable double m_term2;             //!< Help term 2
    mutable double m_term3;             //!< Help term 3
    mutable double m_last_ridx;        // Last reciprocal radial index
    mutable double m_inv_ridx;         // Spatial profile index

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialEllipticalGeneralGauss").
 ***************************************************************************/
inline
std::string GModelSpatialEllipticalGeneralGauss::classname(void) const
{
    return ("GModelSpatialEllipticalGeneralGauss");
}

/***********************************************************************//**
 * @brief Return ridx
 *
 * @return ridx
 *
 * Returns the reciprocal of the radial profile index.
 ***************************************************************************/
inline
double GModelSpatialEllipticalGeneralGauss::ridx(void) const
{
    return (m_ridx.value());
}

/***********************************************************************//**
 * @brief Set reciprocal index
 *
 * @param[in] ridx.
 *
 * Sets the reciprocal index of the radial profile.
 ***************************************************************************/
inline
void GModelSpatialEllipticalGeneralGauss::ridx(const double& ridx)
{
    m_ridx.value(ridx);
    return;
}


#endif /* GMODELSPATIALELLIPTICALGENERALGAUSS_HPP */
