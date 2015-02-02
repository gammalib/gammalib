/***************************************************************************
 *      GCTACubeSourcePoint.hpp - CTA cube analysis point source class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Juergen Knoedlseder                         *
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
 * @file GCTACubeSourcePoint.hpp
 * @brief CTA cube analysis point source class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTACUBESOURCEPOINT_HPP
#define GCTACUBESOURCEPOINT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GSkyDir.hpp"
#include "GCTACubeSource.hpp"
#include "GNodeArray.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GModelSpatial;
class GObservation;


/***********************************************************************//**
 * @class GCTACubeSourcePoint
 *
 * @brief CTA source cube base class
 *
 * This class handles pre-computed response information for a point source
 * in a cube-style analysis. It derives from the abstract GCTACubeSource
 * class.
 ***************************************************************************/
class GCTACubeSourcePoint : public GCTACubeSource {

public:
    // Constructors and destructors
    GCTACubeSourcePoint(void);
    GCTACubeSourcePoint(const GCTACubeSourcePoint& source);
    virtual ~GCTACubeSourcePoint(void);

    // Operators
    GCTACubeSourcePoint& operator=(const GCTACubeSourcePoint & source);

    // Implemented pure virtual methods
    void                 clear(void);
    GCTACubeSourcePoint* clone(void) const;
    std::string          classname(void) const;
    GCTAClassCode        code(void) const;
    void                 set(const std::string&   name,
                             const GModelSpatial& model,
                             const GObservation&  obs);
    std::string          print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double         aeff(const int& index) const;
    double         delta(const int& index) const;
    double         psf(const int& ieng, const double& delta) const;
    const GSkyDir& dir(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTACubeSourcePoint& source);
    void free_members(void);

    // Data members
    GSkyDir             m_dir;       //!< Source direction
    std::vector<double> m_aeff;      //!< Deadtime corrected effective area (cm2)
    std::vector<double> m_delta_map; //!< Distance of bin from source (radians)
    std::vector<double> m_psf;       //!< Point spread function
    GNodeArray          m_deltas;    //!< Delta node array for PSF
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTACubeSourcePoint").
 ***************************************************************************/
inline
std::string GCTACubeSourcePoint::classname(void) const
{
    return ("GCTACubeSourcePoint");
}


/***********************************************************************//**
 * @brief Return class type code
 *
 * @return GCTA_SOURCE_CUBE_POINT_SOURCE.
 *
 * Returns the class type code GCTA_CUBE_SOURCE_POINT.
 ***************************************************************************/
inline
GCTAClassCode GCTACubeSourcePoint::code(void) const
{
    return (GCTA_CUBE_SOURCE_POINT);
}


/***********************************************************************//**
 * @brief Return deadtime corrected effective area in cm2
 *
 * @param[in] index Energy index [0,...,m_aeff.size()-1]
 * @return Deadtime corrected effective area (cm2).
 *
 * Returns the deadtime corrected effective area.
 ***************************************************************************/
inline
double GCTACubeSourcePoint::aeff(const int& index) const
{
    double aeff = (index >= 0 && index < (int)m_aeff.size()) ? m_aeff[index] : 0.0;
    return (aeff);
}


/***********************************************************************//**
 * @brief Return angular distance to source in radians
 *
 * @param[in] index Spatial index [0,...,m_delta_map.size()-1]
 * @return Angular distance to source (radians).
 *
 * Returns the angular seperation between a event cube pixel with @p index
 * and the source position. If the index is not within the valid range, 99.9
 * is returned.
 ***************************************************************************/
inline
double GCTACubeSourcePoint::delta(const int& index) const
{
    double delta = (index >= 0 && index < (int)m_delta_map.size())
                   ? m_delta_map[index] : 99.9;
    return (delta);
}


/***********************************************************************//**
 * @brief Return point source sky direction
 *
 * @return Point source sky direction.
 ***************************************************************************/
inline
const GSkyDir& GCTACubeSourcePoint::dir(void) const
{
    return (m_dir);
}

#endif /* GCTACUBESOURCEPOINT_HPP */
