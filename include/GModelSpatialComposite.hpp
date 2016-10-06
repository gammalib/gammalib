/***************************************************************************
 *     GModelSpatialPointSource.hpp - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Domenico Tiziani                                 *
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
 * @file GModelSpatialComposite.hpp
 * @brief Spatial composite model class interface definition
 * @author Domenico Tiziani
 */

#ifndef GMODELSPATIALCOMPOSITE_HPP
#define GMODELSPATIALCOMPOSITE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialComposite
 *
 * @brief Spatial composite model
 *
 * This class implements the spatial model for a composition of spatial
 * models.
 *
 ***************************************************************************/
class GModelSpatialComposite : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialComposite(void);
    GModelSpatialComposite(const bool& dummy, const std::string& type);
    explicit GModelSpatialComposite(const GXmlElement& xml);
    GModelSpatialComposite(const GModelSpatialComposite& model);
    virtual ~GModelSpatialComposite(void);

    // Operators
    virtual GModelSpatialComposite& operator=(const GModelSpatialComposite& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialComposite* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual GClassCode                code(void) const;
    virtual double                    eval(const GPhoton& photon,
                                           const bool& gradients = false) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual double                    mc_norm(const GSkyDir& dir,
                                              const double&  radius) const;
    virtual bool                      contains(const GSkyDir& dir,
                                               const double&  margin = 0.0) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    int components(void) const;
    void append(const GModelSpatial& component, const std::string& name = "");
    const GModelSpatial* component(const int& index) const;
    const GModelSpatial* component(const std::string& name) const;


protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialComposite& model);
    void free_members(void);

    // Protected members
    std::string m_type;   //!< Model type
    std::vector<GModelSpatial*> m_components; //!< Components
    std::vector<std::string> m_component_names; //!< Component names
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialComposite").
 ***************************************************************************/
inline
std::string GModelSpatialComposite::classname(void) const
{
    return ("GModelSpatialComposite");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the spatial model.
 ***************************************************************************/
inline
std::string GModelSpatialComposite::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Return class code
 *
 * @return GModelSpatialComposite.
 *
 * Returns the code GModelSpatialComposite of the class.
 ***************************************************************************/
inline
GClassCode GModelSpatialComposite::code(void) const
{
    return GMODEL_SPATIAL_COMPOSITE;
}

/***********************************************************************//**
 * @brief Return normalization of point source for Monte Carlo simulations
 *
 * @param[in] dir Centre of simulation cone.
 * @param[in] radius Radius of simulation cone (degrees).
 * @return Normalization.
 *
 * Returns the normalization for a point source within a circular region.
 * The normalization is 1 if the point source falls within the circle
 * defined by @p dir and @p radius, 0 otherwise.
 ***************************************************************************/
inline
double GModelSpatialComposite::mc_norm(const GSkyDir& dir,
                                         const double&  radius) const
{
    //double norm = (dir.dist_deg(this->dir()) <= radius) ? 1.0 : 0.0;
    //return (norm);
    return 1.0;
}

#endif /* GMODELSPATIALCOMPOSITE_HPP */
