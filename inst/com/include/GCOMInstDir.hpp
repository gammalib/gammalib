/***************************************************************************
 *          GCOMInstDir.hpp - COMPTEL instrument direction class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
 * @file GCOMInstDir.hpp
 * @brief COMPTEL instrument direction class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMINSTDIR_HPP
#define GCOMINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GCOMInstDir
 *
 * @brief Interface for the COMPTEL instrument direction class
 ***************************************************************************/
class GCOMInstDir : public GInstDir {

public:
    // Constructors and destructors
    GCOMInstDir(void);
    GCOMInstDir(const GCOMInstDir& dir);
    GCOMInstDir(const GSkyDir& dir, const double& phibar);
    virtual ~GCOMInstDir(void);

    // Operators
    GCOMInstDir& operator=(const GCOMInstDir& dir);

    // Methods
    virtual void         clear(void);
    virtual GCOMInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual std::string  print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void           dir(const GSkyDir& dir);
    const GSkyDir& dir(void) const;
    void           phibar(const double& phibar);
    const double&  phibar(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMInstDir& dir);
    void free_members(void);

    // Protected members
    GSkyDir m_dir;     //!< Observed scatter direction of event
    double  m_phibar;  //!< Observed scatter angle of event
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMInstDir").
 ***************************************************************************/
inline
std::string GCOMInstDir::classname(void) const
{
    return ("GCOMInstDir");
}


/***********************************************************************//**
 * @brief Return event scatter direction
 *
 * @return Event scatter direction.
 *
 * Returns the event scatter direction.
 ***************************************************************************/
inline
const GSkyDir& GCOMInstDir::dir(void) const
{
    return m_dir;
}


/***********************************************************************//**
 * @brief Set event scatter direction
 *
 * @param[in] dir Event scatter direction.
 *
 * Set the event scatter direction.
 ***************************************************************************/
inline
void GCOMInstDir::dir(const GSkyDir& dir)
{
    m_dir = dir;
    return;
}


/***********************************************************************//**
 * @brief Return event scatter angle
 *
 * @return Event scatter angle.
 *
 * Returns the event scatter angle.
 ***************************************************************************/
inline
const double& GCOMInstDir::phibar(void) const
{
    return m_phibar;
}


/***********************************************************************//**
 * @brief Set event scatter angle
 *
 * @param[in] dir Event scatter angle.
 *
 * Set the event scatter angle.
 ***************************************************************************/
inline
void GCOMInstDir::phibar(const double& phibar)
{
    m_phibar = phibar;
    return;
}

#endif /* GCOMINSTDIR_HPP */
