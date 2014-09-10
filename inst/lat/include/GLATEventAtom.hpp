/***************************************************************************
 *              GLATEventAtom.hpp - Fermi/LAT event atom class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
 * @file GLATEventAtom.hpp
 * @brief Fermi/LAT event atom class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATEVENTATOM_HPP
#define GLATEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GEventAtom.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GLATInstDir.hpp"


/***********************************************************************//**
 * @class GLATEventAtom
 *
 * @brief Fermi/LAT event atom class
 ***************************************************************************/
class GLATEventAtom : public GEventAtom {

    // Friend classes
    friend class GLATEventList;

public:
    // Constructors and destructors
    GLATEventAtom(void);
    GLATEventAtom(const GLATEventAtom& atom);
    virtual ~GLATEventAtom(void);

    // Operators
    GLATEventAtom& operator=(const GLATEventAtom& atom);

    // Implemented pure virtual base class methods
    void               clear(void);
    GLATEventAtom*     clone(void) const;
    std::string        classname(void) const;
    const GLATInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATEventAtom& atom);
    void free_members(void);

    // Protected members
    GLATInstDir m_dir;                 //!< Event direction
    GEnergy     m_energy;              //!< Event energy
    GTime       m_time;                //!< Event time
    float       m_theta;               //!< Zenith angle in instrument system
    float       m_phi;                 //!< Azimuth angle in instrument system
    float       m_zenith_angle;        //!< Zenith angle in Earth system
    float       m_earth_azimuth_angle; //!< Azimuth angle in Earth system
    long        m_event_id;            //!< ID number of original event
    long        m_run_id;              //!< Run number of original event
    short       m_recon_version;       //!< Version of event reconstruction software
    short       m_calib_version[3];    //!< Version of calibration tables for ACD, CAL
    short       m_event_class;         //!< Event class: 0, 1, 2, ...
    short       m_conversion_type;     //!< Type of conversion: 0=Front, 1=Back
    double      m_livetime;            //!< Accumulated livetime since mission start
    double*     m_difrsp;              //!< Diffuse response components
    int         m_num_difrsp;          //!< Number of diffuse model components
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATEventAtom").
 ***************************************************************************/
inline
std::string GLATEventAtom::classname(void) const
{
    return ("GLATEventAtom");
}


/***********************************************************************//**
 * @brief Return event instrument direction
 *
 * @return Event instrument direction.
 *
 * Returns the reconstructed arrival direction of the photon on the sky.
 ***************************************************************************/
inline
const GLATInstDir& GLATEventAtom::dir(void) const
{
    return m_dir;
}


/***********************************************************************//**
 * @brief Return event energy
 *
 * @return Event energy.
 *
 * Returns the reconstructed energy of the photon on the sky.
 ***************************************************************************/
inline
const GEnergy& GLATEventAtom::energy(void) const
{
    return m_energy;
}


/***********************************************************************//**
 * @brief Return event time
 *
 * @return Event time.
 *
 * Returns the event triggering time.
 ***************************************************************************/
inline
const GTime& GLATEventAtom::time(void) const
{
    return m_time;
}

#endif /* GLATEVENTATOM_HPP */
