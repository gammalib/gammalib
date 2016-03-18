/***************************************************************************
 *   GModelSpatialRadialProfileGauss.hpp - Gaussian radial profile class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialRadialProfileGauss.hpp
 * @brief Radial Gaussian profile model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALRADIALPROFILEGAUSS_HPP
#define GMODELSPATIALRADIALPROFILEGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadialProfile.hpp"
#include "GModelPar.hpp"

/* __ Forward declaration ________________________________________________ */
class GXmlElement;


/**************************************************************************
 * @class GModelSpatialRadialProfileGauss
 *
 * @brief Radial Gaussian profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a Gaussian radial profile.
 ***************************************************************************/
class GModelSpatialRadialProfileGauss : public GModelSpatialRadialProfile {

public:
    // Constructors and destructors
    GModelSpatialRadialProfileGauss(void);
    explicit GModelSpatialRadialProfileGauss(const GXmlElement& xml);
    GModelSpatialRadialProfileGauss(const GModelSpatialRadialProfileGauss& model);
    virtual ~GModelSpatialRadialProfileGauss(void);

    // Operators
    virtual GModelSpatialRadialProfileGauss& operator=(const GModelSpatialRadialProfileGauss& model);

    // Implemented pure virtual base class methods
    virtual void                             clear(void);
    virtual GModelSpatialRadialProfileGauss* clone(void) const;
    virtual std::string                      classname(void) const;
    virtual std::string                      type(void) const;
    virtual double                           theta_min(void) const;
    virtual double                           theta_max(void) const;
    virtual void                             read(const GXmlElement& xml);
    virtual void                             write(GXmlElement& xml) const;
    virtual std::string                      print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GModelSpatialRadialProfileGauss& model);
    void           free_members(void);
    virtual double profile_value(const double& theta) const;

    // Protected members
    GModelPar m_sigma;      //!< Gaussian width (deg)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialRadialProfileGauss").
 ***************************************************************************/
inline
std::string GModelSpatialRadialProfileGauss::classname(void) const
{
    return ("GModelSpatialRadialProfileGauss");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "GaussianProfile".
 *
 * Returns the type of the radial profile model.
 ***************************************************************************/
inline
std::string GModelSpatialRadialProfileGauss::type(void) const
{
    return "GaussianProfile";
}

#endif /* GMODELSPATIALRADIALPROFILEGAUSS_HPP */
