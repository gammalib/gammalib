/***************************************************************************
 *                        GPhoton.hpp - Photon class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2021 by Juergen Knoedlseder                         *
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
 * @file GPhoton.hpp
 * @brief Photon class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPHOTON_HPP
#define GPHOTON_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPolarization.hpp"


/***********************************************************************//**
 * @class GPhoton
 *
 * @brief Class that handles photons.
 *
 * The GPhoton class stores the physical attributes of a photon such as the
 * photon arrival direction, its energy and its arrival time. This class is
 * mainly used for Monte Carlo simulations.
 ***************************************************************************/
class GPhoton : public GBase {

    // Operator friends
    friend bool operator==(const GPhoton &a, const GPhoton &b);
    friend bool operator!=(const GPhoton &a, const GPhoton &b);

public:
    // Constructors and destructors
    GPhoton(void);
    GPhoton(const GSkyDir&       dir,
            const GEnergy&       energy,
            const GTime&         time,
            const GPolarization& polarization,
            const int&           mc_id = -1);
    GPhoton(const GPhoton& photon);
    virtual ~GPhoton(void);
 
    // Operators
    GPhoton& operator=(const GPhoton& photon);

    // Methods
    void                 clear(void);
    GPhoton*             clone(void) const;
    std::string          classname(void) const;
    const GSkyDir&       dir(void) const;
    const GEnergy&       energy(void) const;
    const GTime&         time(void) const;
    const GPolarization& polarization(void) const;
    const int&           mc_id(void) const;
    void                 dir(const GSkyDir& dir);
    void                 energy(const GEnergy& energy);
    void                 time(const GTime& time);
    void                 polarization(const GPolarization& polarization);
    void                 mc_id(const int& mc_id);
    std::string          print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPhoton& photon);
    void free_members(void);

    // Protected data members
    GSkyDir       m_dir;           //!< Photon arrival direction
    GEnergy       m_energy;        //!< Photon energy
    GTime         m_time;          //!< Photon arrival time
    GPolarization m_polarization;  //!< Polarization
    int           m_mc_id;         //!< Monte Carlo simulation origin
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GPhoton").
 ***************************************************************************/
inline
std::string GPhoton::classname(void) const
{
    return ("GPhoton");
}


/***********************************************************************//**
 * @brief Return photon sky direction
 *
 * @return Sky direction of photon.
 ***************************************************************************/
inline
const GSkyDir& GPhoton::dir(void) const
{
    return m_dir;
}


/***********************************************************************//**
 * @brief Return photon energy
 *
 * @return Energy of photon.
 ***************************************************************************/
inline
const GEnergy& GPhoton::energy(void) const
{
    return m_energy;
}


/***********************************************************************//**
 * @brief Return photon time
 *
 * @return Time of photon.
 ***************************************************************************/
inline
const GTime& GPhoton::time(void) const
{
    return m_time;
}


/***********************************************************************//**
 * @brief Return photon polarization
 *
 * @return Polarization of photon.
 ***************************************************************************/
inline
const GPolarization& GPhoton::polarization(void) const
{
    return m_polarization;
}


/***********************************************************************//**
 * @brief Return photon Monte-Carlo identifier
 *
 * @return Photon Monte-Carlo identifier.
 ***************************************************************************/
inline
const int& GPhoton::mc_id(void) const
{
    return m_mc_id;
}


/***********************************************************************//**
 * @brief Set photon sky direction
 *
 * @param[in] dir Sky direction of photon.
 ***************************************************************************/
inline
void GPhoton::dir(const GSkyDir& dir)
{
    m_dir = dir;
    return;
}


/***********************************************************************//**
 * @brief Set photon energy
 *
 * @param[in] energy Photon energy.
 ***************************************************************************/
inline
void GPhoton::energy(const GEnergy& energy)
{
    m_energy = energy;
    return;
}


/***********************************************************************//**
 * @brief Set photon time
 *
 * @param[in] time Photon time.
 ***************************************************************************/
inline
void GPhoton::time(const GTime& time)
{
    m_time = time;
    return;
}


/***********************************************************************//**
 * @brief Set photon polarization
 *
 * @param[in] polarization Photon polarization.
 ***************************************************************************/
inline
void GPhoton::polarization(const GPolarization& polarization)
{
    m_polarization = polarization;
    return;
}


/***********************************************************************//**
 * @brief Set photon Monte-Carlo identifier
 *
 * @param[in] mc_id Photon Monte-Carlo identifier.
 ***************************************************************************/
inline
void GPhoton::mc_id(const int& mc_id)
{
    m_mc_id = mc_id;
    return;
}


/***********************************************************************//**
 * @brief Equality friend operator
 *
 * @param[in] a First photon.
 * @param[in] b Second photon.
 * @return True is first photon is identical to second photon
 ***************************************************************************/
inline
bool operator==(const GPhoton &a, const GPhoton &b)
{
    return (a.m_energy == b.m_energy && a.m_time == b.m_time &&
            a.m_dir.dist(b.m_dir) == 0.0);
}


/***********************************************************************//**
 * @brief Non-equality friend operator
 *
 * @param[in] a First photon.
 * @param[in] b Second photon.
 * @return True is first photon is not identical to second photon
 ***************************************************************************/
inline
bool operator!=(const GPhoton &a, const GPhoton &b)
{
    return (a.m_energy != b.m_energy ||  a.m_time != b.m_time ||
            a.m_dir.dist(b.m_dir) > 0.0);
}

#endif /* GPHOTON_HPP */
