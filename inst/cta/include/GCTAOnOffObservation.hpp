/***************************************************************************
 *          GCTAOnOffObservation.hpp - CTA on-off observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Chia-Chun Lu & Christoph Deil                    *
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
 * @file GCTAOnOffObservation.hpp
 * @brief CTA on-off observation class definition
 * @author Chia-Chun Lu & Christoph Deil
 */

#ifndef GCTAONOFFOBSERVATION_HPP
#define GCTAONOFFOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GPha.hpp"
#include "GArf.hpp"
#include "GRmf.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAObservation.hpp"
#include "GSkyRegions.hpp"


/***********************************************************************//**
 * @class GCTAOnOffObservation
 *
 * @brief CTA on-off observation class
 ***************************************************************************/
class GCTAOnOffObservation : public GBase {

public:
    // Constructors and destructors
    GCTAOnOffObservation(void);
    GCTAOnOffObservation(const GEbounds& ereco, const GSkyRegions& on,
                         const GSkyRegions& off);
    GCTAOnOffObservation(const GCTAOnOffObservation& obs);
    virtual ~GCTAOnOffObservation(void);
 
    // Operators
    GCTAOnOffObservation& operator=(const GCTAOnOffObservation& obs);

    // Methods
    void                  clear(void);
    GCTAOnOffObservation* clone(void) const;
    void                  name(const std::string& name);
    void                  instrument(const std::string& instrument);
    void                  id(const std::string& id);
    void                  on_regions(const GSkyRegions& regions);
    void                  off_regions(const GSkyRegions& regions);
    const std::string&    name(void) const;
    const std::string&    instrument(void) const;
    const std::string&    id(void) const;
    const GPha&           on_spec(void) const;
    const GPha&           off_spec(void) const;
    const GArf&           arf(void) const;
    const GRmf&           rmf(void) const;
    void                  fill(const GCTAObservation& obs);
    void                  compute_response(const GCTAObservation& obs,
                                           const GEbounds& etrue);
    void                  read(const GXmlElement& xml);
    void                  write(GXmlElement& xml) const;
    std::string           print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAOnOffObservation& obs);
    void free_members(void);
    void compute_arf(const GCTAObservation& obs);
    void compute_rmf(const GCTAObservation& obs, const GEbounds& etrue);

    // Protected data members
    std::string m_name;         //!< Name
    std::string m_instrument;   //!< Instrument name
    std::string m_id;           //!< Observation identifier
    GPha 		m_on_spec;
    GPha 		m_off_spec;
    GArf        m_arf;
    GRmf        m_rmf;
    GSkyRegions m_on_regions;
    GSkyRegions m_off_regions;
};


/***********************************************************************//**
 * @brief Set name of observation
 *
 * @param[in] name Observation name.
 ***************************************************************************/
inline
void GCTAOnOffObservation::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Set instrument
 *
 * @param[in] instrument Instrument.
 ***************************************************************************/
inline
void GCTAOnOffObservation::instrument(const std::string& instrument)
{
    m_instrument = instrument;
    return;
}


/***********************************************************************//**
 * @brief Set observation identifier
 *
 * @param[in] id Observation identifier.
 ***************************************************************************/
inline
void GCTAOnOffObservation::id(const std::string& id)
{
    m_id = id;
    return;
}


/***********************************************************************//**
 * @brief Set ON regions
 *
 * @param[in] regions ON regions.
 ***************************************************************************/
inline
void GCTAOnOffObservation::on_regions(const GSkyRegions& regions)
{
    m_on_regions = regions;
    return;
}


/***********************************************************************//**
 * @brief Set OFF regions
 *
 * @param[in] regions OFF regions.
 ***************************************************************************/
inline
void GCTAOnOffObservation::off_regions(const GSkyRegions& regions)
{
    m_off_regions = regions;
    return;
}


/***********************************************************************//**
 * @brief Return name of observation
 *
 * @return Observation name.
 ***************************************************************************/
inline
const std::string& GCTAOnOffObservation::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Return instrument
 *
 * @return Instrument.
 ***************************************************************************/
inline
const std::string& GCTAOnOffObservation::instrument(void) const
{
    return m_instrument;
}


/***********************************************************************//**
 * @brief Return observation identifier
 *
 * @return Observation identifier.
 ***************************************************************************/
inline
const std::string& GCTAOnOffObservation::id(void) const
{
    return m_id;
}


/***********************************************************************//**
 * @brief Return ON spectrum
 *
 * @return ON spectrum.
 ***************************************************************************/
inline
const GPha& GCTAOnOffObservation::on_spec(void) const
{
    return m_on_spec;
}


/***********************************************************************//**
 * @brief Return OFF spectrum
 *
 * @return OFF spectrum.
 ***************************************************************************/
inline
const GPha& GCTAOnOffObservation::off_spec(void) const
{
    return m_off_spec;
}


/***********************************************************************//**
 * @brief Return Auxiliary Response File
 *
 * @return Auxiliary Response File.
 ***************************************************************************/
inline
const GArf& GCTAOnOffObservation::arf(void) const
{
    return m_arf;
}


/***********************************************************************//**
 * @brief Return Redistribution Matrix File
 *
 * @return Redistribution Matrix File.
 ***************************************************************************/
inline
const GRmf& GCTAOnOffObservation::rmf(void) const
{
    return m_rmf;
}

#endif /* GCTAONOFFOBSERVATION_HPP */
