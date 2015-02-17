/***************************************************************************
 *   GModelSpatialEllipticalGauss.hpp - Elliptical gauss source model class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Michael Mayer                               *
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
 * @file GModelSpatialEllipticalGauss.hpp
 * @brief Elliptical gauss model class interface definition
 * @author Michael Mayer
 */

#ifndef GMODELSPATIALELLIPTICALGAUSS_HPP
#define GMODELSPATIALELLIPTICALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialElliptical.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelSpatialEllipticalGauss
 *
 * @brief Elliptical gauss source model class
 *
 * This class implements the spatial component of the factorised source
 * model for an elliptical gauss source, i.e. constant surface brightness
 * according to an asymmetric Gaussian
 ***************************************************************************/
class GModelSpatialEllipticalGauss : public GModelSpatialElliptical {

public:
    // Constructors and destructors
    GModelSpatialEllipticalGauss(void);
    explicit GModelSpatialEllipticalGauss(const GSkyDir& dir,
                                         const double&  major,
                                         const double&  minor,
                                         const double&  posangle);
    explicit GModelSpatialEllipticalGauss(const GXmlElement& xml);
    GModelSpatialEllipticalGauss(const GModelSpatialEllipticalGauss& model);
    virtual ~GModelSpatialEllipticalGauss(void);

    // Operators
    virtual GModelSpatialEllipticalGauss& operator=(const GModelSpatialEllipticalGauss& model);

    // Implemented pure virtual base class methods
    virtual void                         clear(void);
    virtual GModelSpatialEllipticalGauss* clone(void) const;
    virtual std::string                  classname(void) const;
    virtual std::string                  type(void) const;
    virtual double                       eval(const double&  theta,
                                              const double&  posangle,
                                              const GEnergy& energy,
                                              const GTime&   time) const;
    virtual double                       eval_gradients(const double&  theta,
                                                        const double&  posangle,
                                                        const GEnergy& energy,
                                                        const GTime&   time) const;
    virtual GSkyDir                      mc(const GEnergy& energy,
                                            const GTime& time,
                                            GRan& ran) const;
    virtual double                       theta_max(void) const;
    virtual void                         read(const GXmlElement& xml);
    virtual void                         write(GXmlElement& xml) const;
    virtual std::string                  print(const GChatter& chatter = NORMAL) const;


protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialEllipticalGauss& model);
    void free_members(void);
    void update(void) const;

    // Cached members used for pre-computations
    mutable double m_last_minor;   //!< Last semi-minor axis
    mutable double m_last_major;   //!< Last semi-major axis
    mutable double m_last_posangle;   //!< Last position angle
    mutable double m_last_posangle_rad;   //!< Last position angle in radians
    mutable double m_cospos2; //!< squared cosine of position angle
    mutable double m_sinpos2; //!< squared sine of position angle
    mutable double m_sin2pos; //!< sine of twice the position angle
    mutable double m_minor2; //!< square of minor axis
    mutable double m_major2; //!< square of major axis
    mutable double m_minor_rad;    //!< Minor axis in radians
    mutable double m_major_rad;    //!< Major axis in radians
    mutable double m_norm;             //!< Normalization
    mutable double m_term1;    //!< Help term 1
    mutable double m_term2;    //!< Help term 2
    mutable double m_term3;    //!< Help term 3

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialEllipticalGauss").
 ***************************************************************************/
inline
std::string GModelSpatialEllipticalGauss::classname(void) const
{
    return ("GModelSpatialEllipticalGauss");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "EllipticalGauss".
 *
 * Returns the type of the elliptical gauss model.
 ***************************************************************************/
inline
std::string GModelSpatialEllipticalGauss::type(void) const
{
    return "EllipticalGauss";
}

#endif /* GMODELSPATIALELLIPTICALGAUSS_HPP */
