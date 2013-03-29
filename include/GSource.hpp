/***************************************************************************
 *                        GSource.hpp - Source class                       *
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
 * @file GSource.hpp
 * @brief Source class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSOURCE_HPP
#define GSOURCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GModelSpatial.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GSource
 *
 * @brief Class that handles gamma-ray sources
 *
 * The GSource class stores the physical attributes of a gamma-ray source,
 * which is the distribution of photon arrival directions (i.e. the spatial
 * source model), the photon energy, and the photon arrival time. The class
 * is very similar to the GPhoton class, yet instead of a given photon
 * arrival direction it contains a sky model member.
 *
 * Each source also has a name so that it can be identified.
 ***************************************************************************/
class GSource : public GBase {

public:
    // Constructors and destructors
    GSource(void);
    explicit GSource(const std::string& name,
                     GModelSpatial*     model,
                     const GEnergy&     energy,
                     const GTime&       time);
    GSource(const GSource& src);
    virtual ~GSource(void);
 
    // Operators
    GSource& operator= (const GSource& src);

    // Methods
    void                 clear(void);
    GSource*             clone(void) const;
    std::string          name(void) const { return m_name; }
    const GModelSpatial* model(void) const { return m_model; }
    const GEnergy&       energy(void) const { return m_energy; }
    const GTime&         time(void) const { return m_time; }
    void                 name(const std::string& name) { m_name=name; }
    void                 model(GModelSpatial* model) { m_model=model; }
    void                 energy(const GEnergy& energy) { m_energy=energy; }
    void                 time(const GTime& time) { m_time=time; }
    std::string          print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSource& src);
    void free_members(void);

    // Protected data members
    std::string    m_name;     //!< Source name
    GModelSpatial* m_model;    //!< Spatial model
    GEnergy        m_energy;   //!< Photon energy
    GTime          m_time;     //!< Photon arrival time
};

#endif /* GSOURCE_HPP */
