/***************************************************************************
 *                  GLATPointing.hpp - LAT pointing class                  *
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
 * @file GLATPointing.hpp
 * @brief GLATPointing class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GLATPOINTING_HPP
#define GLATPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include "GPointing.hpp"


/***********************************************************************//**
 * @class GLATPointing
 *
 * @brief Interface for the LAT pointing class.
 *
 * The LAT pointing class contains information for a specific LAT pointing.
 ***************************************************************************/
class GLATPointing : public GPointing {

public:
    // Constructors and destructors
    GLATPointing(void);
    GLATPointing(const GLATPointing& pnt);
    virtual ~GLATPointing(void);

    // Operators
    GLATPointing& operator= (const GLATPointing& pnt);

    // Implemented pure virtual methods
    void           clear(void);
    GLATPointing*  clone(void) const;
    const GSkyDir& dir(void) const { return m_dir; }
    std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void dir(const GSkyDir& dir) { m_dir=dir; }

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GLATPointing& pnt);
    void          free_members(void);

    // Protected members
    GSkyDir m_dir;  //!< Pointing direction
};

#endif /* GLATPOINTING_HPP */
