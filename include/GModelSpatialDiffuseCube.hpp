/***************************************************************************
 *        GModelSpatialDiffuseCube.hpp - Spatial map cube model class      *
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
 * @file GModelSpatialDiffuseCube.hpp
 * @brief Spatial map cube model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALDIFFUSECUBE_HPP
#define GMODELSPATIALDIFFUSECUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialDiffuse.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GSkymap.hpp"
#include "GXmlElement.hpp"
#include "GEbounds.hpp"

/***********************************************************************//**
 * @class GModelSpatialDiffuseCube
 *
 * @brief Spatial map cube model
 *
 * This class implements the spatial component of the factorized source
 * model for a map cube. A map cube is a set of sky maps for different
 * energies.
 *
 * @todo Loading not yet implemented. Introduce a loading flag to load cube
 *       only when it is required. This avoid unnecessary memory allocation.
 * @todo eval() and eval_gradients() methods not yet implemented.
 ***************************************************************************/
class GModelSpatialDiffuseCube : public GModelSpatialDiffuse {

public:
    // Constructors and destructors
    GModelSpatialDiffuseCube(void);
    explicit GModelSpatialDiffuseCube(const GXmlElement& xml);
    explicit GModelSpatialDiffuseCube(const std::string& filename,
                                      const double&      value = 1.0);
    explicit GModelSpatialDiffuseCube(const GSkymap& map,
                                      const double&  value = 1.0);
    GModelSpatialDiffuseCube(const GModelSpatialDiffuseCube& model);
    virtual ~GModelSpatialDiffuseCube(void);

    // Operators
    virtual GModelSpatialDiffuseCube& operator=(const GModelSpatialDiffuseCube& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialDiffuseCube* clone(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GPhoton& photon) const;
    virtual double                    eval_gradients(const GPhoton& photon) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void               load(const std::string& filename);
    double             value(void) const;
    void               value(const double& value);
    const std::string& filename(void) const;
    void               filename(const std::string& filename);
    const GSkymap&     cube(void) const;
    void               cube(const GSkymap& map);
    const GEbounds& ebounds(void) const;
    void               ebounds(const GEbounds& bounds);
    bool               isloaded(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialDiffuseCube& model);
    void free_members(void);
    void prepare_cube(void);

    // Protected members
    GModelPar   m_value;      //!< Value
    std::string m_filename;   //!< Name of map cube
    GSkymap     m_cube;       //!< Map cube
    GEbounds   m_ebounds;  //!< Energy bounds of the maps
    bool        m_loaded;     //!< Signals that map cube has been loaded
    std::vector<double> m_mc_cache; //!< Monte Carlo cache
};


/***********************************************************************//**
 * @brief Return spatial model type
 *
 * @return "MapCubeFunction".
 *
 * Returns the type of the spatial map cube model.
 ***************************************************************************/
inline
std::string GModelSpatialDiffuseCube::type(void) const
{
    return "MapCubeFunction";
}


/***********************************************************************//**
 * @brief Get model value
 *
 * @return Model value.
 *
 * Returns the value of the spatial map cube model.
 ***************************************************************************/
inline
double GModelSpatialDiffuseCube::value(void) const
{
    return (m_value.value());
}


/***********************************************************************//**
 * @brief Set model value
 *
 * @param[in] value Model value.
 *
 * Set the value of the spatial map cube model.
 ***************************************************************************/
inline
void GModelSpatialDiffuseCube::value(const double& value)
{
    m_value.value(value);
    return;
}


/***********************************************************************//**
 * @brief Get file name
 *
 * @return File name.
 *
 * Returns the file name of the spatial map cube model.
 ***************************************************************************/
inline
const std::string& GModelSpatialDiffuseCube::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Set file name
 *
 * @param[in] filename File name.
 *
 * Set the file name of the spatial map cube model.
 ***************************************************************************/
inline
void GModelSpatialDiffuseCube::filename(const std::string& filename)
{
    m_filename = filename;
    return;
}


/***********************************************************************//**
 * @brief Get map cube
 *
 * @return Map cube.
 *
 * Returns the map cube.
 ***************************************************************************/
inline
const GSkymap& GModelSpatialDiffuseCube::cube(void) const
{
    return (m_cube);
}


/***********************************************************************//**
 * @brief Set map cube
 *
 * @param[in] map Sky map.
 *
 * Set the map cube of the spatial map cube model.
 ***************************************************************************/
inline
void GModelSpatialDiffuseCube::cube(const GSkymap& map)
{
    m_cube   = map;
    m_loaded = true;
    return;
}

/***********************************************************************//**
 * @brief Get map cube
 *
 * @return Map cube.
 *
 * Returns the map cube.
 ***************************************************************************/
inline
const GEbounds& GModelSpatialDiffuseCube::ebounds(void) const
{
    return (m_ebounds);
}


/***********************************************************************//**
 * @brief Set map cube
 *
 * @param[in] map Sky map.
 *
 * Set the map cube of the spatial map cube model.
 ***************************************************************************/
inline
void GModelSpatialDiffuseCube::ebounds(const GEbounds& bounds)
{
    m_ebounds   = bounds;
    return;
}



/***********************************************************************//**
 * @brief Signal if map cube has been loaded
 *
 * @return True if map cube has been loaded, false otherwise.
 *
 * Signals if a map cube is present (either by loading it from a file or by
 * assigning it from a GSkymap).
 ***************************************************************************/
inline
bool GModelSpatialDiffuseCube::isloaded(void) const
{
    return (m_loaded);
}

#endif /* GMODELSPATIALDIFFUSECUBE_HPP */
