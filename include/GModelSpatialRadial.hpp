/***************************************************************************
 *    GModelSpatialRadial.hpp - Abstract radial spatial model base class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialRadial.hpp
 * @brief Abstract radial spatial model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALRADIAL_HPP
#define GMODELSPATIALRADIAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPhoton.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GModelSpatialRadial
 *
 * @brief Abstract radial spatial model base class
 *
 * This class defines the interface for a radial model as spatial component
 * of the factorized source model. Typical examples of radial components are
 * axisymmetric Disk, Gaussian or Shell sources.
 ***************************************************************************/
class GModelSpatialRadial : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialRadial(void);
    GModelSpatialRadial(const GModelSpatialRadial& model);
    explicit GModelSpatialRadial(const GXmlElement& xml);
    virtual ~GModelSpatialRadial(void);

    // Operators
    virtual GModelSpatialRadial& operator=(const GModelSpatialRadial& model);

    // Pure virtual methods
    virtual void                 clear(void) = 0;
    virtual GModelSpatialRadial* clone(void) const = 0;
    virtual std::string          classname(void) const = 0;
    virtual std::string          type(void) const = 0;
    virtual double               eval(const double&  theta,
                                      const GEnergy& energy,
                                      const GTime& time) const = 0;
    virtual double               eval_gradients(const double& theta,
                                                const GEnergy& energy,
                                                const GTime& time) const = 0;
    virtual GSkyDir              mc(const GEnergy& energy,
                                    const GTime& time,
                                    GRan& ran) const = 0;
    virtual double               theta_max(void) const = 0;
    virtual std::string          print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented pure virtual base class methods
    virtual GClassCode code(void) const;
    virtual double     eval(const GPhoton& photon) const;
    virtual double     eval_gradients(const GPhoton& photon) const;
    virtual double     norm(const GSkyDir& dir, const double&  radius) const;
    virtual void       read(const GXmlElement& xml);
    virtual void       write(GXmlElement& xml) const;

    // Other methods
    double  ra(void) const;
    double  dec(void) const;
    void    ra(const double& ra);
    void    dec(const double& dec);
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialRadial& model);
    void free_members(void);

    // Proteced members
    GModelPar m_ra;    //!< Right Ascension (deg)
    GModelPar m_dec;   //!< Declination (deg)
};


/***********************************************************************//**
 * @brief Return class code
 *
 * @return GModelSpatialRadial.
 *
 * Returns the code GModelSpatialRadial of the class.
 ***************************************************************************/
inline
GClassCode GModelSpatialRadial::code(void) const
{
    return GMODEL_SPATIAL_RADIAL;
}


/***********************************************************************//**
 * @brief Return Right Ascencion of model centre
 *
 * @return Right Ascencion of model centre (degrees).
 *
 * Returns the Right Ascension of the model centre in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadial::ra(void) const
{
    return (m_ra.value());
}


/***********************************************************************//**
 * @brief Set Right Ascencion of model centre
 *
 * @param[in] ra Right Ascension (degrees).
 *
 * Sets the Right Ascension of the model centre in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadial::ra(const double& ra)
{
    m_ra.value(ra);
    return;
}


/***********************************************************************//**
 * @brief Return Declination of model centre
 *
 * @return Declination of model centre (degrees).
 *
 * Returns the Declination of the model centre in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadial::dec(void) const
{
    return (m_dec.value());
}


/***********************************************************************//**
 * @brief Set Declination of model centre
 *
 * @param[in] dec Declination (degrees).
 *
 * Sets the Declination of the model centre in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadial::dec(const double& dec)
{
    m_dec.value(dec);
    return;
}


/***********************************************************************//**
 * @brief Return normalization of radial source
 *
 * @return Normalization.
 *
 * Returns the normalization for a radial source within a circular region.
 * The normalization is 1 if the radial source falls within the circle, 0
 * otherwise
 ***************************************************************************/
inline
double GModelSpatialRadial::norm(const GSkyDir& dir,
                                 const double&  radius) const
{
    double norm = (dir.dist_deg(this->dir()) <= radius+theta_max()) ? 1.0 : 0.0;
    return (norm);
}

#endif /* GMODELSPATIALRADIAL_HPP */
