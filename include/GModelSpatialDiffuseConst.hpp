/***************************************************************************
 *      GModelSpatialDiffuseConst.hpp - Spatial isotropic model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
#include "GSkyRegionCircle.hpp"
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
    GModelSpatialDiffuseConst(const bool& dummy, const std::string& type);
    explicit GModelSpatialDiffuseConst(const GXmlElement& xml);
    explicit GModelSpatialDiffuseConst(const double& value);
    GModelSpatialDiffuseConst(const GModelSpatialDiffuseConst& model);
    virtual ~GModelSpatialDiffuseConst(void);

    // Operators
    GModelSpatialDiffuseConst& operator=(const GModelSpatialDiffuseConst& model);

    // Implemented pure virtual base class methods
    virtual void                       clear(void);
    virtual GModelSpatialDiffuseConst* clone(void) const;
    virtual std::string                classname(void) const;
    virtual std::string                type(void) const;
    virtual double                     eval(const GPhoton& photon,
                                            const bool& gradients = false) const;
    virtual GSkyDir                    mc(const GEnergy& energy,
                                          const GTime& time,
                                          GRan& ran) const;
    virtual double                     mc_norm(const GSkyDir& dir,
                                               const double&  radius) const;
    virtual bool                       contains(const GSkyDir& dir,
                                                const double&  margin = 0.0) const;
    virtual GSkyRegion*                region(void) const;
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
    void set_region(void) const;

    // Protected members
    std::string              m_type;          //!< Model type
    GModelPar                m_value;         //!< Value
    mutable GSkyRegionCircle m_region;        //!< Bounding circle
    mutable GSkyDir          m_mc_centre;     //!< Simulation cone centre
    mutable double           m_mc_cos_radius; //!< Cosine of sim. cone radius
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialDiffuseConst").
 ***************************************************************************/
inline
std::string GModelSpatialDiffuseConst::classname(void) const
{
    return ("GModelSpatialDiffuseConst");
}


/***********************************************************************//**
 * @brief Return spatial model type
 *
 * @return Model type.
 *
 * Returns the type of the isotropic spatial model.
 ***************************************************************************/
inline
std::string GModelSpatialDiffuseConst::type(void) const
{
    return (m_type);
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


/***********************************************************************//**
 * @brief Signals whether model contains sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] margin Margin to be added to sky direction (deg, default: 0.0).
 * @return True.
 *
 * Signals whether a sky direction falls within the bounding circle of
 * the diffuse model. As the constant model is defined on the entire sphere,
 * the method returns always true.
 ***************************************************************************/
inline
bool GModelSpatialDiffuseConst::contains(const GSkyDir& dir,
                                         const double&  margin) const
{
    return (true);
}


/***********************************************************************//**
 * @brief Return boundary sky region
 *
 * @return Boundary sky region.
 *
 * Returns a sky region that fully encloses the point source.
 ***************************************************************************/
inline
GSkyRegion* GModelSpatialDiffuseConst::region(void) const
{
    set_region();
    return (&m_region);
}

#endif /* GMODELSPATIALDIFFUSECONST_HPP */
