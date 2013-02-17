/***************************************************************************
 *      GModelSpatialDiffuseConst.hpp - Spatial isotropic model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseConst.hpp
 * @brief Isotropic spatial model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALDIFFUSECONST_HPP
#define GMODELSPATIALDIFFUSECONST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialDiffuse.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialDiffuseConst
 *
 * @brief Isotropic spatial model
 *
 * This class implements the spatial component of the factorised source
 * model for an isotropic source.
 ***************************************************************************/
class GModelSpatialDiffuseConst : public GModelSpatialDiffuse {

public:
    // Constructors and destructors
    GModelSpatialDiffuseConst(void);
    explicit GModelSpatialDiffuseConst(const GXmlElement& xml);
    GModelSpatialDiffuseConst(const GModelSpatialDiffuseConst& model);
    virtual ~GModelSpatialDiffuseConst(void);

    // Operators
    GModelSpatialDiffuseConst& operator= (const GModelSpatialDiffuseConst& model);

    // Implemented pure virtual base class methods
    virtual void                       clear(void);
    virtual GModelSpatialDiffuseConst* clone(void) const;
    virtual std::string                type(void) const { return "ConstantValue"; }
    virtual double                     eval(const GSkyDir& srcDir) const;
    virtual double                     eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir                    mc(GRan& ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;
    virtual std::string                print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialDiffuseConst& model);
    void free_members(void);

    // Protected members
    GModelPar m_value;         //!< Value
};

#endif /* GMODELSPATIALDIFFUSECONST_HPP */
