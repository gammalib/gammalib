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
    void                  name(const std::string& name) { m_name = name; }
    void                  instrument(const std::string& instrument) { m_instrument = instrument; }
    void                  id(const std::string& id) { m_id = id; }
    void                  off_regions(const GSkyRegions& regions) {m_off_regions = regions;}
    void                  on_regions(const GSkyRegions& regions) {m_on_regions = regions;}
    const std::string&    name(void) const { return m_name; }
    const std::string&    instrument(void) const { return m_instrument; }
    const std::string&    id(void) const { return m_id; }
    const GPha&           on_spec(void) const { return m_on_spec; }
    const GPha&           off_spec(void) const { return m_off_spec; }
    const GArf&           arf(void) { return m_arf; };
    const GRmf&           rmf(void) {return m_rmf;}
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

#endif /* GCTAONOFFOBSERVATION_HPP */
