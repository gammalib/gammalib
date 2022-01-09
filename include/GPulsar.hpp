/***************************************************************************
 *                        GPulsar.hpp - Pulsar class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file GPulsar.hpp
 * @brief Pulsar class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPULSAR_HPP
#define GPULSAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GPulsarEphemerides.hpp"

/* __ Forward declarations _______________________________________________ */
class GGti;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GPulsar
 *
 * @brief Pulsar class
 *
 * This class implements a pulsar, defined by a name and a list of pulsar
 * ephemerides.
 ***************************************************************************/
class GPulsar : public GBase {

public:
    // Constructors and destructors
    GPulsar(void);
    GPulsar(const GFilename& filename, const std::string& name = "");
    GPulsar(const GPulsar& pulsar);
    virtual ~GPulsar(void);

    // Operators
    GPulsar& operator=(const GPulsar& pulsar);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GPulsar*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    int  size(void) const;
    bool is_empty(void) const;
    GGti validity(void) const;
    void load(const GFilename& filename, const std::string& name = "");

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPulsar& pulsar);
    void free_members(void);
    void load_fits(const GFilename& filename, const std::string& name = "");
    void load_integral(const GFitsTable* table, const std::string& name = "");
    void load_fermi(const GFitsTable* table, const std::string& name = "");
    void load_psrtime(const GFilename& filename, const std::string& name = "");
    void load_parfile(const GFilename& filename);

    // Protected members
    std::string                     m_name;        //!< Pulsar name
    std::vector<GPulsarEphemerides> m_ephemerides; //!< Pulsar ephemerides
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GPulsar").
 ***************************************************************************/
inline
std::string GPulsar::classname(void) const
{
    return ("GPulsar");
}


/***********************************************************************//**
 * @brief Return number of ephemerides for pulsar
 *
 * @return Number of ephemerides for pulsar.
 *
 * Returns the number of ephemerides for pulsar.
 ***************************************************************************/
inline
int GPulsar::size(void) const
{
    return (int)m_ephemerides.size();
}


/***********************************************************************//**
 * @brief Signals if there are no ephemerides for pulsar
 *
 * @return True if there are no ephemerides for pulsar, false otherwise.
 *
 * Signals if there are no ephemerides for pulsar.
 ***************************************************************************/
inline
bool GPulsar::is_empty(void) const
{
    return (m_ephemerides.empty());
}

#endif /* GPULSAR_HPP */
