/***************************************************************************
 *                 GXXXPointing.hpp  -  XXX pointing class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GXXXPointing.hpp
 * @brief XXX pointing class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXXXPOINTING_HPP
#define GXXXPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GPointing.hpp"
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GXXXPointing
 *
 * @brief Interface for the XXX pointing class.
 ***************************************************************************/
class GXXXPointing : public GPointing {

public:
    // Constructors and destructors
    GXXXPointing(void);
    explicit GXXXPointing(const GSkyDir& dir);
    GXXXPointing(const GXXXPointing& pnt);
    virtual ~GXXXPointing(void);

    // Operators
    virtual GXXXPointing& operator=(const GXXXPointing& pnt);

    // Implemented pure virtual methods
    virtual void           clear(void);
    virtual GXXXPointing*  clone(void) const;
    virtual const GSkyDir& dir(void) const;
    virtual std::string    print(void) const;

    // Other methods
    void dir(const GSkyDir& dir);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXXXPointing& pnt);
    void free_members(void);

    // Protected members
    GSkyDir m_dir;    //!< Pointing direction in sky coordinates
};

#endif /* GXXXPOINTING_HPP */
