/***************************************************************************
 *   GCTACubeSourceDiffuse.hpp - CTA cube analysis diffuse source class    *
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
 * @file GCTACubeSourceDiffuse.hpp
 * @brief CTA cube analysis diffuse source class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTACUBESOURCEDIFFUSE_HPP
#define GCTACUBESOURCEDIFFUSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCTACubeSource.hpp"
#include "GSkymap.hpp"
#include "GFunction.hpp"
#include "GCTAResponseCube.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GModelSpatial;
class GObservation;
class GModelSpatialDiffuse;
class GSkyDir;
class GEnergy;
class GTime;


/***********************************************************************//**
 * @class GCTACubeSourceDiffuse
 *
 * @brief CTA diffuse source cube class
 *
 * This class handles pre-computed response information for a diffuse source
 * in a stacked cube analysis. It derives from the abstract GCTACubeSource
 * class.
 ***************************************************************************/
class GCTACubeSourceDiffuse : public GCTACubeSource {

public:
    // Constructors and destructors
    GCTACubeSourceDiffuse(void);
    GCTACubeSourceDiffuse(const GCTACubeSourceDiffuse& source);
    virtual ~GCTACubeSourceDiffuse(void);

    // Operators
    GCTACubeSourceDiffuse& operator=(const GCTACubeSourceDiffuse & source);

    // Implemented pure virtual methods
    void                   clear(void);
    GCTACubeSourceDiffuse* clone(void) const;
    std::string            classname(void) const;
    GCTAClassCode          code(void) const;
    void                   set(const std::string&   name,
                               const GModelSpatial& model,
                               const GObservation&  obs);
    std::string            print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double par(const int& index) const;
    double irf(const int& pixel, const int& iebin) const;
    double psf(const GCTAResponseCube* rsp,
               const GModelSpatial*    model,
               const GSkyDir&          srcDir,
               const GEnergy&          srcEng,
               const GTime&            srcTime) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTACubeSourceDiffuse& source);
    void free_members(void);

    // Data members
    GSkymap m_cube;  //!< Diffuse map convolved with IRF
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTACubeSourceDiffuse").
 ***************************************************************************/
inline
std::string GCTACubeSourceDiffuse::classname(void) const
{
    return ("GCTACubeSourceDiffuse");
}


/***********************************************************************//**
 * @brief Return class type code
 *
 * @return GCTA_SOURCE_CUBE_DIFFUSE.
 *
 * Returns the class type code GCTA_CUBE_SOURCE_DIFFUSE.
 ***************************************************************************/
inline
GCTAClassCode GCTACubeSourceDiffuse::code(void) const
{
    return (GCTA_CUBE_SOURCE_DIFFUSE);
}


/***********************************************************************//**
 * @brief Return instrument response function
 *
 * @param[in] pixel Spatial pixel index [0,...,???]
 * @param[in] iebin Energy index [0,...,???]
 * @return Instrument response function.
 *
 * Returns the instrument response function.
 ***************************************************************************/
inline
double GCTACubeSourceDiffuse::irf(const int& pixel, const int& iebin) const
{
    return (m_cube(pixel, iebin));
}

#endif /* GCTACUBESOURCEDIFFUSE_HPP */
