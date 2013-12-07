/***************************************************************************
 *                GCOMPointing.hpp - COMPTEL pointing class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GCOMPointing.hpp
 * @brief COMPTEL pointing class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMPOINTING_HPP
#define GCOMPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GCOMPointing
 *
 * @brief Interface for the COMPTEL pointing class.
 *
 * This class implements a COMPTEL pointing. For the time being it is assumed
 * that the pointing direction is time-independent.
 ***************************************************************************/
class GCOMPointing {

public:
    // Constructors and destructors
    GCOMPointing(void);
    explicit GCOMPointing(const GSkyDir& dir);
    GCOMPointing(const GCOMPointing& pnt);
    virtual ~GCOMPointing(void);

    // Operators
    virtual GCOMPointing& operator=(const GCOMPointing& pnt);

    // Implemented pure virtual methods
    virtual void           clear(void);
    virtual GCOMPointing*  clone(void) const;
    virtual const GSkyDir& dir(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void dir(const GSkyDir& dir);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMPointing& pnt);
    void free_members(void);

    // Protected members
    GSkyDir m_dir;    //!< Pointing direction in sky coordinates
};

#endif /* GCOMPOINTING_HPP */
