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
 * model for an isotropic source. The model has a single parameter which
 * is the normalization factor of the model.
 ***************************************************************************/
class GModelSpatialDiffuseConst : public GModelSpatialDiffuse {

public:
    // Constructors and destructors
    GModelSpatialDiffuseConst(void);
    explicit GModelSpatialDiffuseConst(const GXmlElement& xml);
    explicit GModelSpatialDiffuseConst(const double& value);
    GModelSpatialDiffuseConst(const GModelSpatialDiffuseConst& model);
    virtual ~GModelSpatialDiffuseConst(void);

    // Operators
    GModelSpatialDiffuseConst& operator=(const GModelSpatialDiffuseConst& model);

    // Implemented pure virtual base class methods
    virtual void                       clear(void);
    virtual GModelSpatialDiffuseConst* clone(void) const;
    virtual std::string                type(void) const;
    virtual double                     eval(const GPhoton& photon) const;
    virtual double                     eval_gradients(const GPhoton& photon) const;
    virtual GSkyDir                    mc(const GEnergy& energy,
                                          const GTime& time,
                                          GRan& ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;
    virtual std::string                print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double value(void) const;
    void   value(const double& value);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialDiffuseConst& model);
    void free_members(void);

    // Protected members
    GModelPar m_value;    //!< Value
};


/***********************************************************************//**
 * @brief Return spatial model type
 *
 * @return "ConstantValue".
 *
 * Returns the type of the isotropic spatial model.
 ***************************************************************************/
inline
std::string GModelSpatialDiffuseConst::type(void) const
{
    return "ConstantValue";
}


/***********************************************************************//**
 * @brief Get model value
 *
 * @return Model value.
 *
 * Returns the value of the isotropic spatial model.
 ***************************************************************************/
inline
double GModelSpatialDiffuseConst::value(void) const
{
    return (m_value.value());
}


/***********************************************************************//**
 * @brief Set model value
 *
 * @param[in] value Model value.
 *
 * Set the value of the isotropic spatial model.
 ***************************************************************************/
inline
void GModelSpatialDiffuseConst::value(const double& value)
{
    m_value.value(value);
    return;
}

#endif /* GMODELSPATIALDIFFUSECONST_HPP */
