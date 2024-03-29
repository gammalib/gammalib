/***************************************************************************
 *    GModelSpatialRadialGauss.hpp - Radial Gaussian source model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2022 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialRadialGauss.hpp
 * @brief Radial Gaussian model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALRADIALGAUSS_HPP
#define GMODELSPATIALRADIALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadial.hpp"
#include "GModelPar.hpp"

/* __ Forward declarations _______________________________________________ */
class GEnergy;
class GTime;
class GPhoton;
class GRan;
class GSkyDir;
class GSkyRegion;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpatialRadialGauss
 *
 * @brief Radial Gaussian model class
 *
 * This class implements the spatial component of the factorised source
 * model for a Gaussian source.
 ***************************************************************************/
class GModelSpatialRadialGauss : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialGauss(void);
    GModelSpatialRadialGauss(const bool& dummy, const std::string& type);
    GModelSpatialRadialGauss(const GSkyDir&     dir,
                             const double&      sigma,
                             const std::string& coordsys = "CEL");
    explicit GModelSpatialRadialGauss(const GXmlElement& xml);
    GModelSpatialRadialGauss(const GModelSpatialRadialGauss& model);
    virtual ~GModelSpatialRadialGauss(void);

    // Operators
    virtual GModelSpatialRadialGauss& operator=(const GModelSpatialRadialGauss& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialRadialGauss* clone(void) const;
    virtual std::string               classname(void) const;
    virtual double                    eval(const double&  theta,
                                           const GEnergy& energy,
                                           const GTime&   time,
                                           const bool&    gradients = false) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual bool                      contains(const GSkyDir& dir,
                                               const double&  margin = 0.0) const;
    virtual double                    theta_max(void) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  sigma(void) const;
    void    sigma(const double& sigma);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialRadialGauss& model);
    void         free_members(void);
    void         update(const bool& gradients) const;
    virtual void set_region(void) const;

    // Protected members
    GModelPar m_sigma;                 //!< Gaussian width (deg)

    // Cached members used for pre-computations
    mutable double m_last_sigma;       //!< Last Gaussian sigma
    mutable double m_inv_sigma2_rad;   //!< sigma(rad)^-2
    mutable double m_value_norm;       //!< (2pi sigma(rad))^-2
    mutable double m_g_theta_norm;     //!< sigma(rad)^-2 deg2rad
    mutable double m_g_sigma_norm;     //!< sigma(rad)^-1 sigma_scale deg2rad
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialRadialGauss").
 ***************************************************************************/
inline
std::string GModelSpatialRadialGauss::classname(void) const
{
    return ("GModelSpatialRadialGauss");
}


/***********************************************************************//**
 * @brief Return Gaussian sigma
 *
 * @return Gaussian sigma (degrees).
 *
 * Returns the Gaussian sigma in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadialGauss::sigma(void) const
{
    return (m_sigma.value());
}


/***********************************************************************//**
 * @brief Set Gaussian sigma
 *
 * @param[in] sigma Gaussian sigma (degrees).
 *
 * Sets the Gaussian sigma in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadialGauss::sigma(const double& sigma)
{
    m_sigma.value(sigma);
    return;
}

#endif /* GMODELSPATIALRADIALGAUSS_HPP */
