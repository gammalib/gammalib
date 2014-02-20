/***************************************************************************
 *          GModelSpatialDiffuseMap.hpp - Spatial map model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
#include "GSkymap.hpp"
#include "GXmlElement.hpp"


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
    explicit GModelSpatialDiffuseMap(const GXmlElement& xml);
    GModelSpatialDiffuseMap(const std::string& filename,
                            const double&      value = 1.0,
                            const bool&        normalize = true);
    GModelSpatialDiffuseMap(const GSkymap& map,
                            const double&  value = 1.0,
                            const bool&    normalize = true);
    GModelSpatialDiffuseMap(const GModelSpatialDiffuseMap& model);
    virtual ~GModelSpatialDiffuseMap(void);

    // Operators
    virtual GModelSpatialDiffuseMap& operator= (const GModelSpatialDiffuseMap& model);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GModelSpatialDiffuseMap* clone(void) const;
    virtual std::string              type(void) const;
    virtual double                   eval(const GPhoton& photon) const;
    virtual double                   eval_gradients(const GPhoton& photon) const;
    virtual GSkyDir                  mc(const GEnergy& energy,
                                        const GTime& time,
                                        GRan& ran) const;
    virtual double                   norm(const GSkyDir& dir,
                                          const double&  radius) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double             value(void) const;
    void               value(const double& value);
    const std::string& filename(void) const;
    void               load(const std::string& filename);
    const GSkymap&     map(void) const;
    void               map(const GSkymap& map);
    bool               normalize(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialDiffuseMap& model);
    void free_members(void);
    void prepare_map(void);

    // Protected members
    GModelPar           m_value;      //!< Value
    GSkymap             m_map;        //!< Skymap
    std::string         m_filename;   //!< Name of skymap
    std::vector<double> m_mc_cache;   //!< Monte Carlo cache
    bool                m_normalize;  //!< Normalize map (default: true)
    double              m_norm;       //!< Map normalization
};

/***********************************************************************//**
 * @brief Return spatial model type
 *
 * @return "SpatialMap".
 *
 * Returns the type of the spatial map model.
 ***************************************************************************/
inline
std::string GModelSpatialDiffuseMap::type(void) const
{
    return "SpatialMap";
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
const std::string& GModelSpatialDiffuseMap::filename(void) const
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
const GSkymap& GModelSpatialDiffuseMap::map(void) const
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
void GModelSpatialDiffuseMap::map(const GSkymap& map)
{
    m_map = map;
    prepare_map();
    return;
}


/***********************************************************************//**
 * @brief Return normalization flag
 *
 * @return True is the map has been normalized, false otherwise.
 *
 * Signals whether a map has been normalized or not.
 ***************************************************************************/
inline
bool GModelSpatialDiffuseMap::normalize(void) const
{
    return (m_normalize);
}


/***********************************************************************//**
 * @brief Return normalization of diffuse map
 *
 * @return Normalization.
 *
 * Returns the normalization of a diffuse map. The normalization is given
 * by the model value times the integrated flux in the sky map.
 ***************************************************************************/
inline
double GModelSpatialDiffuseMap::norm(const GSkyDir& dir,
                                     const double&  radius) const
{
    return (m_norm * value());
}

#endif /* GMODELSPATIALDIFFUSEMAP_HPP */
