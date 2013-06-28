/***************************************************************************
 *           GCTAOnOffObservation.i - CTA on-off observation class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Michael Mayer                                    *
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
 * @file GCTAOnOffObservation.i
 * @brief CTA on-off observation class definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAOnOffObservation.hpp"
%}


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
};


/***********************************************************************//**
 * @brief GCTOnOffAObservation class extension
 ***************************************************************************/
%extend GCTAOnOffObservation {
    GCTAOnOffObservation copy() {
        return (*self);
    }
}
