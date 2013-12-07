/***************************************************************************
 *                  GCTAPointing.hpp - CTA pointing class                  *
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
 * @file GCTAPointing.hpp
 * @brief CTA pointing class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAPOINTING_HPP
#define GCTAPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GSkyDir.hpp"
#include "GTime.hpp"
#include "GMatrix.hpp"


/***********************************************************************//**
 * @class GCTAPointing
 *
 * @brief Interface for the CTA pointing class.
 *
 * This class implements a CTA pointing. For the time being it is assumed
 * that the pointing direction is time-independent.
 ***************************************************************************/
class GCTAPointing {

public:
    // Constructors and destructors
    GCTAPointing(void);
    explicit GCTAPointing(const GSkyDir& dir);
    GCTAPointing(const GCTAPointing& pnt);
    virtual ~GCTAPointing(void);

    // Operators
    virtual GCTAPointing& operator=(const GCTAPointing& pnt);

    // Implemented pure virtual methods
    virtual void           clear(void);
    virtual GCTAPointing*  clone(void) const;
    virtual const GSkyDir& dir(void) const { return m_dir; }
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void   dir(const GSkyDir& dir);
    const  GMatrix& rot(void) const;
    double zenith(void) const { return m_zenith; }
    double azimuth(void) const { return m_azimuth; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAPointing& pnt);
    void free_members(void);
    void update(void) const;

    // Protected members
    GSkyDir         m_dir;        //!< Pointing direction in sky coordinates
    double          m_zenith;     //!< Pointing zenith angle
    double          m_azimuth;    //!< Pointing azimuth angle

    // Cached members
    mutable bool    m_has_cache;  //!< Has transformation cache
    mutable GMatrix m_Rback;      //!< Rotation matrix
};

#endif /* GCTAPOINTING_HPP */
