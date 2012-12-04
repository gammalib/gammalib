/***************************************************************************
 *         GCOMInstDir.hpp  -  COMPTEL instrument direction class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
    virtual ~GCOMInstDir(void);

    // Operators
    GCOMInstDir& operator= (const GCOMInstDir& dir);

    // Methods
    virtual void         clear(void);
    virtual GCOMInstDir* clone(void) const;
    virtual std::string  print(void) const;

    // Other methods
    void           skydir(const GSkyDir& dir) { m_dir=dir; }
    void           phi(const double& phi) { m_phi=phi; }
    const GSkyDir& skydir(void) const { return m_dir; }
    const double&  phi(void) const { return m_phi; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMInstDir& dir);
    void free_members(void);

    // Protected members
    GSkyDir   m_dir;   //!< Observed scatter direction of event
    double    m_phi;   //!< Observed scatter angle of event
};

#endif /* GCOMINSTDIR_HPP */
