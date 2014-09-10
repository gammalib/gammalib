/***************************************************************************
 *   GModelSpatialDiffuse.hpp - Abstract diffuse spatial model base class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuse.hpp
 * @brief Abstract diffuse spatial model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALDIFFUSE_HPP
#define GMODELSPATIALDIFFUSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GPhoton.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GModelSpatialDiffuse
 *
 * @brief Abstract diffuse spatial model base class
 *
 * This class defines the interface for a diffuse model as spatial component
 * of the factorized source model.
 ***************************************************************************/
class GModelSpatialDiffuse : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialDiffuse(void);
    GModelSpatialDiffuse(const GModelSpatialDiffuse& model);
    virtual ~GModelSpatialDiffuse(void);

    // Operators
    virtual GModelSpatialDiffuse& operator=(const GModelSpatialDiffuse& model);

    // Pure virtual base class methods
    virtual void                  clear(void) = 0;
    virtual GModelSpatialDiffuse* clone(void) const = 0;
    virtual std::string           classname(void) const = 0;
    virtual std::string           type(void) const = 0;
    virtual double                eval(const GPhoton& photon) const = 0;
    virtual double                eval_gradients(const GPhoton& photon) const = 0;
    virtual GSkyDir               mc(const GEnergy& energy, const GTime& time,
                                     GRan& ran) const = 0;
    virtual double                norm(const GSkyDir& dir,
                                       const double&  radius) const = 0;
    virtual void                  read(const GXmlElement& xml) = 0;
    virtual void                  write(GXmlElement& xml) const = 0;
    virtual std::string           print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented virtual base class methods
    virtual GClassCode code(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialDiffuse& model);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Return class code
 *
 * @return GModelSpatialDiffuse.
 *
 * Returns the code GModelSpatialDiffuse of the class.
 ***************************************************************************/
inline
GClassCode GModelSpatialDiffuse::code(void) const
{
    return GMODEL_SPATIAL_DIFFUSE;
}

#endif /* GMODELSPATIALDIFFUSE_HPP */
