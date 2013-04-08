/***************************************************************************
 *    GModelSpatialRadialGauss.hpp - Radial Gaussian source model class    *
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
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GXmlElement.hpp"


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
    explicit GModelSpatialRadialGauss(const GSkyDir& dir, const double& sigma);
    explicit GModelSpatialRadialGauss(const GXmlElement& xml);
    GModelSpatialRadialGauss(const GModelSpatialRadialGauss& model);
    virtual ~GModelSpatialRadialGauss(void);

    // Operators
    virtual GModelSpatialRadialGauss& operator=(const GModelSpatialRadialGauss& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialRadialGauss* clone(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const double&  theta,
                                           const GEnergy& energy,
                                           const GTime& time) const;
    virtual double                    eval_gradients(const double& theta,
                                                     const GEnergy& energy,
                                                     const GTime& time) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual double                    theta_max(void) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double  sigma(void) const;
    void    sigma(const double& sigma);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialRadialGauss& model);
    void free_members(void);

    // Protected members
    GModelPar m_sigma;      //!< Gaussian width (deg)
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "GaussFunction".
 *
 * Returns the type of the radial Gauss model.
 ***************************************************************************/
inline
std::string GModelSpatialRadialGauss::type(void) const
{
    return "GaussFunction";
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
