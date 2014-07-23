/***************************************************************************
 *        GCTASourceCubeDiffuse.hpp - CTA diffuse source cube class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTASourceCubeDiffuse.hpp
 * @brief CTA diffuse source cube class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTASOURCECUBEDIFFUSE_HPP
#define GCTASOURCECUBEDIFFUSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCTASourceCube.hpp"
#include "GSkymap.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GModelSpatial;
class GObservation;


/***********************************************************************//**
 * @class GCTASourceCubeDiffuse
 *
 * @brief CTA diffuse source cube class
 *
 * This class handles pre-computed response information for a diffuse source
 * in a stacked cube analysis. It derives from the abstract GCTASourceCube
 * class.
 ***************************************************************************/
class GCTASourceCubeDiffuse : public GCTASourceCube {

public:
    // Constructors and destructors
    GCTASourceCubeDiffuse(void);
    GCTASourceCubeDiffuse(const GCTASourceCubeDiffuse& cube);
    virtual ~GCTASourceCubeDiffuse(void);

    // Operators
    GCTASourceCubeDiffuse& operator=(const GCTASourceCubeDiffuse & cube);

    // Implemented pure virtual methods
    void                   clear(void);
    GCTASourceCubeDiffuse* clone(void) const;
    virtual GCTAClassCode  code(void) const;
    void                   set(const std::string&   name,
                               const GModelSpatial& model,
                               const GObservation&  obs);
    std::string            print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double irf(const int& pixel, const int& iebin) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTASourceCubeDiffuse& cube);
    void free_members(void);

    // Data members
    GSkymap m_cube;  //!< Diffuse map convolved with IRF
};


/***********************************************************************//**
 * @brief Return class type code
 *
 * @return GCTA_SOURCE_CUBE_DIFFUSE.
 *
 * Returns the class type code GCTA_SOURCE_CUBE_DIFFUSE.
 ***************************************************************************/
inline
GCTAClassCode GCTASourceCubeDiffuse::code(void) const
{
    return (GCTA_SOURCE_CUBE_DIFFUSE);
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
double GCTASourceCubeDiffuse::irf(const int& pixel, const int& iebin) const
{
    return (m_cube(pixel, iebin));
}

#endif /* GCTASOURCECUBEDIFFUSE_HPP */
