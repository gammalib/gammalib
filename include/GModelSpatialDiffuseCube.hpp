/***************************************************************************
 *        GModelSpatialDiffuseCube.hpp - Spatial map cube model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
#include "GModelSpectral.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GSkymap.hpp"
#include "GNodeArray.hpp"
#include "GXmlElement.hpp"
#include "GEbounds.hpp"
#include "GEnergies.hpp"


/***********************************************************************//**
 * @class GModelSpatialDiffuseCube
 *
 * @brief Spatial map cube model
 *
 * This class implements the spatial component of the factorised source
 * model for a map cube. A map cube is a set of sky maps for different
 * energies. It is assumed that the pixels of the sky map are given in the
 * units ph/cm2/s/sr/MeV. 
 ***************************************************************************/
class GModelSpatialDiffuseCube : public GModelSpatialDiffuse {

public:
    // Constructors and destructors
    GModelSpatialDiffuseCube(void);
    explicit GModelSpatialDiffuseCube(const GXmlElement& xml);
    GModelSpatialDiffuseCube(const std::string& filename,
                             const double&      value = 1.0);
    GModelSpatialDiffuseCube(const GSkymap&   cube,
                             const GEnergies& energies,
                             const double&    value = 1.0);
    GModelSpatialDiffuseCube(const GModelSpatialDiffuseCube& model);
    virtual ~GModelSpatialDiffuseCube(void);

    // Operators
    virtual GModelSpatialDiffuseCube& operator=(const GModelSpatialDiffuseCube& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialDiffuseCube* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GPhoton& photon) const;
    virtual double                    eval_gradients(const GPhoton& photon) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual double                    norm(const GSkyDir& dir,
                                           const double&  radius) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    int                        maps(void) const;
    int                        pixels(void) const;
    void                       load(const std::string& filename);
    double                     value(void) const;
    void                       value(const double& value);
    const std::string&         filename(void) const;
    void                       filename(const std::string& filename);
    const GSkymap&             cube(void) const;
    void                       cube(const GSkymap& cube);
    GEnergies                  energies(void);
    void                       energies(const GEnergies& energies);
    const GModelSpectralNodes& spectrum(void) const;
    void                       set_mc_cone(const GSkyDir& centre,
                                           const double&  radius);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialDiffuseCube& model);
    void free_members(void);
    void fetch_cube(void) const;
    void set_energy_boundaries(void);
    void update_mc_cache(void);

    // Protected members
    GModelPar           m_value;       //!< Value
    std::string         m_filename;    //!< Name of map cube
    bool                m_loaded;      //!< Signals if map cube has been loaded
    GSkymap             m_cube;        //!< Map cube
    GNodeArray          m_logE;        //!< Log10(energy) values of the maps
    GEbounds            m_ebounds;     //!< Energy bounds of the maps
    std::vector<double> m_mc_cache;    //!< Monte Carlo cache
    GModelSpectralNodes m_mc_spectrum; //!< Map cube spectrum
    GSkyDir             m_mc_cone_dir; //!< Monte Carlo simulation cone centre
    double              m_mc_cone_rad; //!< Monte Carlo simulation cone radius
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialDiffuseCube").
 ***************************************************************************/
inline
std::string GModelSpatialDiffuseCube::classname(void) const
{
    return ("GModelSpatialDiffuseCube");
}


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
 * @brief Return number of maps in cube
 *
 * @return Number of maps in cube.
 *
 * Returns the number of maps in the cube.
 ***************************************************************************/
inline
int GModelSpatialDiffuseCube::maps(void) const
{
    return (m_cube.nmaps());
}


/***********************************************************************//**
 * @brief Return number of pixels in cube
 *
 * @return Number of pixels in cube.
 *
 * Returns the number of pixels in the cube.
 ***************************************************************************/
inline
int GModelSpatialDiffuseCube::pixels(void) const
{
    return (m_cube.npix());
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
 * @param[in] cube Sky map.
 *
 * Set the map cube of the spatial map cube model.
 ***************************************************************************/
inline
void GModelSpatialDiffuseCube::cube(const GSkymap& cube)
{
    m_filename.clear();
    m_loaded = false;
    m_cube   = cube;
    update_mc_cache();
    return;
}


/***********************************************************************//**
 * @brief Get map cube spectrum
 *
 * @return Map cube spectrum.
 *
 * Returns the map cube spectrum.
 ***************************************************************************/
inline
const GModelSpectralNodes& GModelSpatialDiffuseCube::spectrum(void) const
{
    return (m_mc_spectrum);
}

/***********************************************************************//**
 * @brief Return normalization of diffuse cube
 *
 * @return Normalization.
 *
 * Returns the normalization of a diffuse cube. The normalization is given
 * by the model value.
 ***************************************************************************/
inline
double GModelSpatialDiffuseCube::norm(const GSkyDir& dir,
                                      const double&  radius) const
{
    return (value());
}

#endif /* GMODELSPATIALDIFFUSECUBE_HPP */
