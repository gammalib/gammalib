/***************************************************************************
 *                   GEphemerides.hpp - Ephemerides class                  *
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
 * @file GEphemerides.hpp
 * @brief Ephemerides class definition
 * @author Juergen Knoedlseder
 */

#ifndef GEPHEMERIDES_HPP
#define GEPHEMERIDES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GFilename.hpp"
#include "GTime.hpp"
#include "GVector.hpp"

/* __ Forward declarations _______________________________________________ */
class GSkyDir;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GEphemerides
 *
 * @brief Ephemerides class
 *
 * This class implements the JPL ephemerides.
 ***************************************************************************/
class GEphemerides : public GBase {

public:
    // Constructors and destructors
    GEphemerides(void);
    GEphemerides(const GEphemerides& ephemerides);
    virtual ~GEphemerides(void);

    // Operators
    GEphemerides& operator=(const GEphemerides& ephemerides);

    // Implemented pure virtual base class methods
    virtual void          clear(void);
    virtual GEphemerides* clone(void) const;
    virtual std::string   classname(void) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

    // Other methods
    int                size(void) const;
    bool               is_empty(void) const;
    const std::string& name(void) const;
    void               name(const std::string& name);
    void               load(const GFilename& filename);
    void               ephemeris(const GTime& time,
                                 GVector*     rce,
                                 GVector*     rcs,
                                 GVector*     vce,
                                 double*      etut) const;
    double             geo2ssb(const GTime&   time,
                               const GSkyDir& srcdir) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEphemerides& ephemerides);
    void free_members(void);
    void fetch_data(void);
    
    // Protected members
    std::string          m_name;      //!< Ephemerides (e.g. DE200)
    GFilename            m_filename;  //!< Ephemerides filename
    GTime                m_tstart;    //!< Ephemerides validity start time
    GTime                m_tstop;     //!< Ephemerides validity stop time
    std::vector<GTime>   m_times;     //!< Times of vectors
    std::vector<GVector> m_earth;     //!< Earth vectors
    std::vector<GVector> m_earth_dt;  //!< First derivative of Earth vectors
    std::vector<GVector> m_earth_d2t; //!< Second derivative of Earth vectors
    std::vector<GVector> m_earth_d3t; //!< Third derivative of Earth vectors
    std::vector<GVector> m_sun;       //!< Sun vectors
    std::vector<double>  m_tdb2tt;    //!< TBD to TT conversion term
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GEphemerides").
 ***************************************************************************/
inline
std::string GEphemerides::classname(void) const
{
    return ("GEphemerides");
}


/***********************************************************************//**
 * @brief Return number of ephemerides
 *
 * @return Number of ephemerides.
 *
 * Returns the number of ephemerides.
 ***************************************************************************/
inline
int GEphemerides::size(void) const
{
    return (int)m_times.size();
}


/***********************************************************************//**
 * @brief Signals if there are no ephemerides
 *
 * @return True if there are no ephemerides, false otherwise.
 *
 * Signals if there are no ephemerides.
 ***************************************************************************/
inline
bool GEphemerides::is_empty(void) const
{
    return (m_times.empty());
}


/***********************************************************************//**
 * @brief Return ephemerides name
 *
 * @return Ephemerides name
 ***************************************************************************/
inline
const std::string& GEphemerides::name(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Set ephemerides name
 *
 * @param[in] name Ephemerides name.
 ***************************************************************************/
inline
void GEphemerides::name(const std::string& name)
{
    m_name = name;
    return;
}

#endif /* GEPHEMERIDES_HPP */
