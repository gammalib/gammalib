/***************************************************************************
 *  GModelSpatialElliptical.hpp - Abstract elliptical spatial model class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialElliptical.hpp
 * @brief Abstract elliptical spatial model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALELLIPTICAL_HPP
#define GMODELSPATIALELLIPTICAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GModelSpatialElliptical
 *
 * @brief Abstract elliptical spatial model base class
 *
 * This class defines the interface for an elliptical model as spatial
 * component of the factorized source model. Typical examples of elliptical
 * components are elliptical Disk, Gaussian or Shell shaped sources.
 ***************************************************************************/
class GModelSpatialElliptical : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialElliptical(void);
    GModelSpatialElliptical(const GModelSpatialElliptical& model);
    explicit GModelSpatialElliptical(const GXmlElement& xml);
    virtual ~GModelSpatialElliptical(void);

    // Operators
    virtual GModelSpatialElliptical& operator=(const GModelSpatialElliptical& model);

    // Pure virtual methods
    virtual void                     clear(void) = 0;
    virtual GModelSpatialElliptical* clone(void) const = 0;
    virtual std::string              type(void) const = 0;
    virtual double                   eval(const double& theta,
                                          const double& posangle) const = 0;
    virtual double                   eval_gradients(const double& theta,
                                                    const double& posangle) const = 0;
    virtual GSkyDir                  mc(GRan& ran) const = 0;
    virtual double                   theta_max(void) const = 0;
    virtual std::string              print(void) const = 0;

    // Implemented virtual methods
    virtual double eval(const GSkyDir& srcDir) const;
    virtual double eval_gradients(const GSkyDir& srcDir) const;
    virtual void   read(const GXmlElement& xml);
    virtual void   write(GXmlElement& xml) const;

    // Other methods
    double  ra(void) const { return m_ra.value(); }
    double  dec(void) const { return m_dec.value(); }
    double  posangle(void) const { return m_posangle.value(); }
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
    void    posangle(const double& posangle) { m_posangle.value(posangle); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialElliptical& model);
    void free_members(void);

    // Proteced members
    GModelPar m_ra;       //!< Right Ascension (deg)
    GModelPar m_dec;      //!< Declination (deg)
    GModelPar m_posangle; //!< Position angle from North, counterclockwise (deg)
};

#endif /* GMODELSPATIALELLIPTICAL_HPP */
