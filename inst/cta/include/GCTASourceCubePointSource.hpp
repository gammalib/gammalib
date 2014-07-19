/***************************************************************************
 *       GCTASourceCubePointSource.hpp - CTA point source cube class       *
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
 * @file GCTASourceCubePointSource.hpp
 * @brief CTA point source cube class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTASOURCECUBEPOINTSOURCE_HPP
#define GCTASOURCECUBEPOINTSOURCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCTASourceCube.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GModelSpatial;
class GObservation;


/***********************************************************************//**
 * @class GCTASourceCubePointSource
 *
 * @brief CTA source cube base class
 *
 * This class handles pre-computed response information for a point source
 * in a cube-style analysis. It derives from the abstract GCTASourceCube
 * class.
 ***************************************************************************/
class GCTASourceCubePointSource : public GCTASourceCube {

public:
    // Constructors and destructors
    GCTASourceCubePointSource(void);
    GCTASourceCubePointSource(const GCTASourceCubePointSource& cube);
    virtual ~GCTASourceCubePointSource(void);

    // Operators
    GCTASourceCubePointSource& operator=(const GCTASourceCubePointSource & cube);

    // Implemented pure virtual methods
    void                       clear(void);
    GCTASourceCubePointSource* clone(void) const;
    void                       set(const GModelSpatial&  model,
                                   const GObservation&   obs);
    std::string                print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTASourceCubePointSource& cube);
    void free_members(void);

    // Data members
    std::vector<double> m_aeff;      //!< Effective area vector (exposure/ontime)
    std::vector<double> m_delta_map; //!< Distance of bin from source (radians)
};

#endif /* GCTASOURCECUBEPOINTSOURCE_HPP */
