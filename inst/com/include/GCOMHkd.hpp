/***************************************************************************
 *        GCOMHkd.hpp - COMPTEL Housekeeping Data container class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2023 by Juergen Knodlseder                               *
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
 * @file GCOMHkd.hpp
 * @brief COMPTEL Housekeeping Data container class definition
 * @author Juergen Knodlseder
 */

#ifndef GCOMHKD_HPP
#define GCOMHKD_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GTimes.hpp"

/* __ Forward declarations _______________________________________________ */
class GTime;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMHkd
 *
 * @brief COMPTEL Housekeeping Data container class
 *
 * This class holds data for one housekeeping parameter.
 ***************************************************************************/
class GCOMHkd : public GContainer {

public:
    // Constructors and destructors
    GCOMHkd(void);
    explicit GCOMHkd(const std::string& name);
    GCOMHkd(const GCOMHkd& hkd);
    virtual ~GCOMHkd(void);

    // Operators
    GCOMHkd& operator=(const GCOMHkd& hkd);

    // Methods
    void               clear(void);
    GCOMHkd*           clone(void) const;
    std::string        classname(void) const;
    int                size(void) const;
    bool               is_empty(void) const;
    void               append(const GTime& time, const double& value);
    void               remove(const int& index);
    void               reserve(const int& num);
    void               extend(const GCOMHkd& hkd);
    const std::string& name(void) const;
    void               name(const std::string& name);
    const GTime&       time(const int& index) const;
    void               time(const int& index, const GTime& time);
    const double&      value(const int& index) const;
    void               value(const int& index, const double& value);
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMHkd& hkd);
    void free_members(void);

    // Protected members
    std::string         m_name;   //!< Name of housekeeping parameter
    GTimes              m_times;  //!< Times
    std::vector<double> m_values; //!< Values at times
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMHkd").
 ***************************************************************************/
inline
std::string GCOMHkd::classname(void) const
{
    return ("GCOMHkd");
}


/***********************************************************************//**
 * @brief Return number of Housekeeping Data in container
 *
 * @return Number of Housekeeping Data in container.
 *
 * Returns the number of Housekeeping Data in the container.
 ***************************************************************************/
inline
int GCOMHkd::size(void) const
{
    return (int)m_times.size();
}


/***********************************************************************//**
 * @brief Signals if there are no Housekeeping Data in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the Housekeeping Data container does not contain any
 * Housekeeping Data.
 ***************************************************************************/
inline
bool GCOMHkd::is_empty(void) const
{
    return (m_times.is_empty());
}


/***********************************************************************//**
 * @brief Reserves space for Housekeeping Data in container
 *
 * @param[in] num Number of Housekeeping Data.
 *
 * Reserves space for @p num Housekeeping Data in the container.
 ***************************************************************************/
inline
void GCOMHkd::reserve(const int& num)
{
    m_times.reserve(num);
    m_values.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Return Housekeeping Data name
 *
 * @return Housekeeping Data name.
 *
 * Returns the Housekeeping Data name.
 ***************************************************************************/
inline
const std::string& GCOMHkd::name(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Set Housekeeping Data name
 *
 * @param[in] name Housekeeping Data name.
 *
 * Sets the Housekeeping Data name.
 ***************************************************************************/
inline
void GCOMHkd::name(const std::string& name)
{
    m_name = name;
    return;
}

#endif /* GCOMHKD_HPP */
