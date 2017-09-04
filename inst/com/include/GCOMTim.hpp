/***************************************************************************
 *             GCOMTim.hpp - COMPTEL Good Time Intervals class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knodlseder                               *
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
 * @file GCOMTim.hpp
 * @brief COMPTEL Good Time Intervals class definition
 * @author Juergen Knodlseder
 */

#ifndef GCOMTIM_HPP
#define GCOMTIM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GGti.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFitsTable;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMTim
 *
 * @brief COMPTEL Good Time Intervals class
 *
 * @todo Add class description.
 ***************************************************************************/
class GCOMTim : public GBase {

public:
    // Constructors and destructors
    GCOMTim(void);
    GCOMTim(const GCOMTim& tim);
    GCOMTim(const GFilename& filename, const std::string& usage = "",
                                       const std::string& mode  = "");
    virtual ~GCOMTim(void);

    // Operators
    GCOMTim& operator=(const GCOMTim& tim);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMTim*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void        read(const GFitsTable& table, const std::string& usage = "",
                                              const std::string& mode  = "");
    const GGti& gti(void) const;
    void        gti(const GGti& gti);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMTim& tim);
    void free_members(void);
    
    // Protected members
    GGti m_gti; //!< Good Time intervals
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMTim").
 ***************************************************************************/
inline
std::string GCOMTim::classname(void) const
{
    return ("GCOMTim");
}


/***********************************************************************//**
 * @brief Return Good Time Intervals
 *
 * @return Good Time Intervals.
 *
 * Returns the Good Time Intervals.
 ***************************************************************************/
inline
const GGti& GCOMTim::gti(void) const
{
    return m_gti;
}


/***********************************************************************//**
 * @brief Set Good Time Intervals
 *
 * @param[in] gti Good Time Intervals.
 *
 * Sets the Good Time Intervals.
 ***************************************************************************/
inline
void GCOMTim::gti(const GGti& gti)
{
    m_gti = gti;
    return;
}

#endif /* GCOMTIM_HPP */
