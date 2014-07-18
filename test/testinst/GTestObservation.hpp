/***************************************************************************
 *       GTestObservation.hpp  -  Test observation class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Jean-Baptiste Cayrou                        *
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


#ifndef GTestOBSERVATION_HPP
#define GTestOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GTestResponse.hpp"
#include "GTestEventList.hpp"
#include "GTools.hpp"

class GTestObservation : public GObservation {

public:
    // Constructors and destructors
    GTestObservation(void) : GObservation() { 
        init_members();
        return;
    }    
    GTestObservation(const GTestObservation& obs) : GObservation(obs) {
        copy_members(obs);
        return;
    }
    
    virtual ~GTestObservation(void){
        free_members();
        return;
    }

    // Operators
    virtual GTestObservation& operator=(const GTestObservation& obs) {
        if (this != &obs) {
            this->GObservation::operator=(obs);
            free_members();
            copy_members(obs);
        }
        return *this;
    }

    // Implement pure virtual methods    
    virtual void clear(void) {
        free_members();
        this->GObservation::free_members();
        this->GObservation::init_members();
        init_members();
        return;
    }
    virtual GTestObservation* clone(void) const {
        return new GTestObservation(*this);
    }
    virtual void response(const GResponse& rsp) {
        const GTestResponse* testrsp = dynamic_cast<const GTestResponse*>(&rsp);
        if (testrsp == NULL) {
            throw;
        }
        m_response = *(testrsp->clone());
        return;
    }
    virtual const GTestResponse* response(void) const { return &m_response;}
    virtual std::string          instrument(void) const { return m_instrument; }
    virtual double               ontime(void) const { return m_ontime; }
    virtual double               livetime(void) const { return m_ontime; }
    virtual double               deadc(const GTime& time) const { return 1.0; }
    virtual void                 read(const GXmlElement& xml) { return; }
    virtual void                 write(GXmlElement& xml) const { return; }
    virtual void                 ontime(const double& ontime) { m_ontime=ontime; }
    virtual std::string          print(const GChatter& chatter = NORMAL) const {
        std::string result;
        result.append("=== GTestObservation ===");
        result.append("\n"+gammalib::parformat("Instrument Name")+m_instrument);
        GTestEventList* list = dynamic_cast<GTestEventList*>(m_events);
        if(list != NULL) {
            result.append("\n"+gammalib::parformat("Number of events")+gammalib::str(list->number()));
            result.append("\n"+gammalib::parformat("GTI")+list->gti().print());
        }
        return result;
    }

protected:
    // Protected methods
    void init_members(void){
        m_instrument = "Test Instrument";
        m_ontime     = 0.0;
        return;
    }
    void copy_members(const GTestObservation& obs){
        m_instrument = obs.m_instrument;
        m_response   = obs.m_response;
        m_ontime     = obs.m_ontime;
        return;    
    }
    void free_members(void) {
        return;
    }

    // Protected members
    std::string   m_instrument;   //!< Instrument name
    GTestResponse m_response;     //!< Pointer to response functions
    double        m_ontime;       //!< ontime
};

#endif /* GTestOBSERVATION_HPP */
