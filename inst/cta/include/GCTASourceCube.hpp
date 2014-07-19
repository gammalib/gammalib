/***************************************************************************
 *        GCTASourceCube.hpp - Abstract CTA source cube base class         *
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
 * @file GCTASourceCube.hpp
 * @brief Abstract CTA source cube base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTASOURCECUBE_HPP
#define GCTASOURCECUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GModelSpatial;
class GObservation;


/***********************************************************************//**
 * @class GCTASourceCube
 *
 * @brief CTA source cube base class
 *
 * This class provides an abstract base class for handling of pre-computed
 * response information in a cube-style analysis.
 ***************************************************************************/
class GCTASourceCube : public GBase {

public:
    // Constructors and destructors
    GCTASourceCube(void);
    GCTASourceCube(const GCTASourceCube& cube);
    virtual ~GCTASourceCube(void);

    // Operators
    virtual GCTASourceCube& operator=(const GCTASourceCube & cube);

    // Pure virtual methods
    virtual void               clear(void) = 0;
    virtual GCTASourceCube*    clone(void) const = 0;
    virtual void               set(const std::string&   name,
                                   const GModelSpatial& model,
                                   const GObservation&  obs) = 0;
    virtual std::string        print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented methods
    virtual const std::string& name(void) const;
    virtual void               name(const std::string& name);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTASourceCube& cube);
    void free_members(void);

    // Data members
    std::string m_name;   //!< Unique source name
};



/***********************************************************************//**
 * @brief Return the source name
 *
 * @return Source name.
 ***************************************************************************/
inline
const std::string& GCTASourceCube::name(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Set the source name
 *
 * @param[in] name Source name.
 ***************************************************************************/
inline
void GCTASourceCube::name(const std::string& name)
{
    m_name = name;
    return;
}

#endif /* GCTASOURCECUBE_HPP */
