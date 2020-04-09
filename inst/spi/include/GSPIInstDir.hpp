/***************************************************************************
 *        GSPIInstDir.hpp - INTEGRAL/SPI instrument direction class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIInstDir.hpp
 * @brief INTEGRAL/SPI instrument direction class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIINSTDIR_HPP
#define GSPIINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"
#include "GSkyDir.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSPIInstDir
 *
 * @brief INTEGRAL/SPI instrument direction class
 *
 * The INTEGRAL/SPI instrument direction defines the spatial information
 * associated to an event.
 ***************************************************************************/
class GSPIInstDir : public GInstDir {

public:
    // Constructors and destructors
    GSPIInstDir(void);
    GSPIInstDir(const GSkyDir& dir, const int& detid);
    GSPIInstDir(const GSPIInstDir& dir);
    virtual ~GSPIInstDir(void);

    // Operators
    GSPIInstDir& operator=(const GSPIInstDir& dir);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GSPIInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual double       hash(void) const;
    virtual std::string  print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void           dir(const GSkyDir& dir);
    const GSkyDir& dir(void) const;
    void           detid(const int& detid);
    const int&     detid(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSPIInstDir& dir);
    void free_members(void);

    // Protected members
    GSkyDir m_dir;
    int     m_detid;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSPIInstDir").
 ***************************************************************************/
inline
std::string GSPIInstDir::classname(void) const
{
    return ("GSPIInstDir");
}


/***********************************************************************//**
 * @brief Return instrument direction hash value
 *
 * @return Hash value.
 *
 * Returns a hash value that can be used in the response cache.
 ***************************************************************************/
inline
double GSPIInstDir::hash(void) const
{
    double hash = m_dir.ra_deg() + m_dir.dec_deg() * 1.0e2 + m_detid * 1.0e4;
    return hash;
}


/***********************************************************************//**
 * @brief Set pointing direction
 *
 * @param[in] dir Pointing direction.
 *
 * Set the pointing direction.
 ***************************************************************************/
inline
void GSPIInstDir::dir(const GSkyDir& dir)
{
    m_dir = dir;
    return;
}


/***********************************************************************//**
 * @brief Return pointing direction
 *
 * @return Pointing direction.
 *
 * Returns the pointing direction.
 ***************************************************************************/
inline
const GSkyDir& GSPIInstDir::dir(void) const
{
    return (m_dir);
}


/***********************************************************************//**
 * @brief Set detector identifier
 *
 * @param[in] detid Detector identifier.
 *
 * Set the detector identifier.
 ***************************************************************************/
inline
void GSPIInstDir::detid(const int& detid)
{
    m_detid = detid;
    return;
}


/***********************************************************************//**
 * @brief Return detector identifier
 *
 * @return Detector identifier.
 *
 * Returns the detector identifier.
 ***************************************************************************/
inline
const int& GSPIInstDir::detid(void) const
{
    return (m_detid);
}

#endif /* GSPIINSTDIR_HPP */
