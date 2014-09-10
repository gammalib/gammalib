/***************************************************************************
 *                        GSource.hpp - Source class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * which is the distribution of photon arrival directions in form of a
 * spatial source model, the photon energy, and the photon arrival time.
 * The class is very similar to the GPhoton class, yet instead of a given
 * photon arrival direction it contains a spatial model pointer. Note that
 * the class does not allocate or deallocate the spatial model pointer,
 * it just is a carrier of a pointer which is handled outside the class.
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
    GSource& operator=(const GSource& src);

    // Methods
    void                 clear(void);
    GSource*             clone(void) const;
    std::string          classname(void) const;
    const std::string&   name(void) const;
    const GModelSpatial* model(void) const;
    const GEnergy&       energy(void) const;
    const GTime&         time(void) const;
    void                 name(const std::string& name);
    void                 model(GModelSpatial* model);
    void                 energy(const GEnergy& energy);
    void                 time(const GTime& time);
    std::string          print(const GChatter& chatter = NORMAL) const;

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


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSource").
 ***************************************************************************/
inline
std::string GSource::classname(void) const
{
    return ("GSource");
}


/***********************************************************************//**
 * @brief Return model name
 *
 * @return Model name.
 *
 * Returns the model name.
 ***************************************************************************/
inline
const std::string& GSource::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Return spatial model component
 *
 * @return Pointer to spatial model component.
 *
 * Returns spatial model component.
 ***************************************************************************/
inline
const GModelSpatial* GSource::model(void) const
{
    return m_model;
}


/***********************************************************************//**
 * @brief Return photon energy
 *
 * @return Photon energy.
 *
 * Returns photon energy.
 ***************************************************************************/
inline
const GEnergy& GSource::energy(void) const
{
    return m_energy;
}


/***********************************************************************//**
 * @brief Return photon arrival time
 *
 * @return Photon arrival time.
 *
 * Returns the photon arrival time.
 ***************************************************************************/
inline
const GTime& GSource::time(void) const
{
    return m_time;
}


/***********************************************************************//**
 * @brief Set model name
 *
 * @param[in] name Model name.
 *
 * Sets the model name.
 ***************************************************************************/
inline
void GSource::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Set spatial model component
 *
 * @param[in] model Spatial model component.
 *
 * Sets the spatial model component.
 ***************************************************************************/
inline
void GSource::model(GModelSpatial* model)
{
    m_model = model;
    return;
}


/***********************************************************************//**
 * @brief Set photon energy
 *
 * @param[in] energy Photon energy.
 *
 * Sets the photon energy.
 ***************************************************************************/
inline
void GSource::energy(const GEnergy& energy)
{
    m_energy = energy;
    return;
}


/***********************************************************************//**
 * @brief Set photon arrival time
 *
 * @param[in] time Photon arrival time.
 *
 * Sets the photon arrival time.
 ***************************************************************************/
inline
void GSource::time(const GTime& time)
{
    m_time = time;
    return;
}

#endif /* GSOURCE_HPP */
