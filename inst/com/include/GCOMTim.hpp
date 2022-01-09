/***************************************************************************
 *             GCOMTim.hpp - COMPTEL Good Time Intervals class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2022 by Juergen Knodlseder                          *
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
class GFitsBinTable;

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
    explicit GCOMTim(const GGti& gti);
    GCOMTim(const GCOMTim& tim);
    GCOMTim(const GFilename& filename, const std::string& usage = "YES",
                                       const std::string& mode  = "NORMAL");
    virtual ~GCOMTim(void);

    // Operators
    GCOMTim& operator=(const GCOMTim& tim);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMTim*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    bool        contains(const GTime& time) const;
    void        reduce(const GGti& gti);
    const GGti& gti(void) const;
    void        gti(const GGti& gti);
    void        load(const GFilename& filename, const std::string& usage = "YES",
                                                const std::string& mode  = "NORMAL");
    void        save(const GFilename& filename, const bool&        clobber = false) const;
    void        read(const GFitsTable& table, const std::string& usage = "YES",
                                              const std::string& mode  = "NORMAL");
    void        write(GFitsBinTable& table) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMTim& tim);
    void free_members(void);

    // Protected members
    GGti m_gti;     //!< Good Time intervals
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
 * @brief Check if time is comprised in the Good Time Intervals
 *
 * @param[in] time Time.
 * @return True if time is within Good Time Intervals, false otherwise.
 *
 * Checks if a time is comprised in the Good Time Intervals.
 ***************************************************************************/
inline
bool GCOMTim::contains(const GTime& time) const
{
    return (m_gti.contains(time));
}


/***********************************************************************//**
 * @brief Reduces Good Time Intervals to intersection with intervals
 *
 * @param[in] gti Good Time Intervals.
 *
 * Reduces the Good Time Intervals to the intersection with the specified
 * list of Good Time Intervals.
 ***************************************************************************/
inline
void GCOMTim::contains(const GGti& gti)
{
    m_gti.reduce(gti);
    return;
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
