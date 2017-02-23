/***************************************************************************
 *          GModelSpatialDiffuseMap.hpp - Spatial map model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseMap.hpp
 * @brief Spatial map model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALDIFFUSEMAP_HPP
#define GMODELSPATIALDIFFUSEMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelSpatialDiffuse.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GSkyMap.hpp"
#include "GSkyRegionCircle.hpp"
#include "GXmlElement.hpp"
#include "GFilename.hpp"


/***********************************************************************//**
 * @class GModelSpatialDiffuseMap
 *
 * @brief Spatial map model
 *
 * This class implements the spatial component of the factorised source
 * model for a skymap.
 ***************************************************************************/
class GModelSpatialDiffuseMap : public GModelSpatialDiffuse {

public:
    // Constructors and destructors
    GModelSpatialDiffuseMap(void);
    GModelSpatialDiffuseMap(const bool& dummy, const std::string& type);
    explicit GModelSpatialDiffuseMap(const GXmlElement& xml);
    GModelSpatialDiffuseMap(const GFilename& filename,
                            const double&    value = 1.0,
                            const bool&      normalize = true);
    GModelSpatialDiffuseMap(const GSkyMap& map,
                            const double&  value = 1.0,
                            const bool&    normalize = true);
    GModelSpatialDiffuseMap(const GModelSpatialDiffuseMap& model);
    virtual ~GModelSpatialDiffuseMap(void);

    // Operators
    virtual GModelSpatialDiffuseMap& operator= (const GModelSpatialDiffuseMap& model);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GModelSpatialDiffuseMap* clone(void) const;
    virtual std::string              classname(void) const;
    virtual std::string              type(void) const;
    virtual double                   eval(const GPhoton& photon,
                                          const bool& gradients = false) const;
    virtual GSkyDir                  mc(const GEnergy& energy,
                                        const GTime& time,
                                        GRan& ran) const;
    virtual double                   mc_norm(const GSkyDir& dir,
                                             const double&  radius) const;
    virtual bool                     contains(const GSkyDir& dir,
                                              const double&  margin = 0.0) const;
    virtual GSkyRegion*              region(void) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double           value(void) const;
    void             value(const double& value);
    const GFilename& filename(void) const;
    void             load(const GFilename& filename);
    const GSkyMap&   map(void) const;
    void             map(const GSkyMap& map);
    bool             normalize(void) const;
    void             set_mc_cone(const GSkyDir& centre,
                                 const double&  radius) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialDiffuseMap& model);
    void free_members(void);
    void prepare_map(void);
    void set_region(void) const;

    // Protected members
    std::string              m_type;          //!< Model type
    GModelPar                m_value;         //!< Value
    GSkyMap                  m_map;           //!< Skymap
    GFilename                m_filename;      //!< Name of skymap
    GSkyDir                  m_centre;        //!< Centre of bounding circle
    double                   m_radius;        //!< Radius of bounding circle
    bool                     m_normalize;     //!< Normalize map (default: true)
    bool                     m_has_normalize; //!< XML has normalize attribute
    mutable GSkyRegionCircle m_region;        //!< Bounding circle

    // MC simulation cache
    mutable GSkyDir m_mc_centre;           //!< Centre of MC cone
    mutable double  m_mc_radius;           //!< Radius of MC cone
    mutable double  m_mc_one_minus_cosrad; //!< 1-cosine of radius
    mutable double  m_mc_norm;             //!< Map normalization
    mutable double  m_mc_max;              //!< Map maximum for MC
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialDiffuseMap").
 ***************************************************************************/
inline
std::string GModelSpatialDiffuseMap::classname(void) const
{
    return ("GModelSpatialDiffuseMap");
}


/***********************************************************************//**
 * @brief Return spatial model type
 *
 * @return Model type.
 *
 * Returns the type of the spatial map model.
 ***************************************************************************/
inline
std::string GModelSpatialDiffuseMap::type(void) const
{
    return (m_type);
}


/***********************************************************************//**
 * @brief Get model value
 *
 * @return Model value.
 *
 * Returns the value of the spatial map model.
 ***************************************************************************/
inline
double GModelSpatialDiffuseMap::value(void) const
{
    return (m_value.value());
}


/***********************************************************************//**
 * @brief Set model value
 *
 * @param[in] value Model value.
 *
 * Set the value of the spatial map model.
 ***************************************************************************/
inline
void GModelSpatialDiffuseMap::value(const double& value)
{
    m_value.value(value);
    return;
}


/***********************************************************************//**
 * @brief Get file name
 *
 * @return File name.
 *
 * Returns the file name of the spatial map model.
 ***************************************************************************/
inline
const GFilename& GModelSpatialDiffuseMap::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Get map
 *
 * @return Sky map.
 *
 * Returns the sky map.
 ***************************************************************************/
inline
const GSkyMap& GModelSpatialDiffuseMap::map(void) const
{
    return (m_map);
}


/***********************************************************************//**
 * @brief Set map cube
 *
 * @param[in] map Sky map.
 *
 * Set the sky map.
 ***************************************************************************/
inline
void GModelSpatialDiffuseMap::map(const GSkyMap& map)
{
    m_map = map;
    prepare_map();
    return;
}


/***********************************************************************//**
 * @brief Return normalization flag
 *
 * @return True if the map has been normalized, false otherwise.
 *
 * Signals whether a map has been normalized or not.
 ***************************************************************************/
inline
bool GModelSpatialDiffuseMap::normalize(void) const
{
    return (m_normalize);
}


/***********************************************************************//**
 * @brief Return boundary sky region
 *
 * @return Boundary sky region.
 *
 * Returns a sky region that fully encloses the point source.
 ***************************************************************************/
inline
GSkyRegion* GModelSpatialDiffuseMap::region(void) const
{
    set_region();
    return (&m_region);
}

#endif /* GMODELSPATIALDIFFUSEMAP_HPP */
