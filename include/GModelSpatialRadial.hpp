/***************************************************************************
 *    GModelSpatialRadial.hpp - Abstract radial spatial model base class   *
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
    virtual std::string          print(void) const = 0;

    // Implemented virtual base class methods
    virtual double eval(const GPhoton& photon) const;
    virtual double eval_gradients(const GPhoton& photon) const;
    virtual void   read(const GXmlElement& xml);
    virtual void   write(GXmlElement& xml) const;

    // Other methods
    double  ra(void) const;
    double  dec(void) const;
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

#endif /* GMODELSPATIALRADIAL_HPP */
