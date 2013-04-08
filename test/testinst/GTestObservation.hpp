/***************************************************************************
 *       GTestObservation.hpp  -  Test observation class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jean-Baptiste Cayrou                             *
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
#include "GTestPointing.hpp"
#include "GTools.hpp"

class GTestObservation : public GObservation {

public:
    // Constructors and destructors
    GTestObservation(void) : GObservation(){ 
        init_members();
        return;
    }
    
    GTestObservation(const GTestObservation& obs) : GObservation(obs){
        init_members();
        copy_members(obs);

        return;
    }
    
    virtual ~GTestObservation(void){
        free_members();
        return;
    }

    // Operators
    virtual GTestObservation& operator= (const GTestObservation& obs){
        
        // Execute only if object is not identical
        if (this != &obs) {

            // Copy base class members
            this->GObservation::operator=(obs);

             // Free members
            free_members();

            // Initialise members
            init_members();

            // Copy members
            copy_members(obs);

        } // endif: object was not identical

        // Return this object
        return *this;
    }

    // Implement pure virtual methods
    
    virtual void clear(void){
        free_members();
        this->GObservation::free_members();

        // Initialise members
        this->GObservation::init_members();
        init_members();
        
        return;
    }
    virtual GTestObservation* clone(void) const{ return new GTestObservation(*this);}
    virtual void              response(const GResponse& rsp){
        // Get pointer on Test response
        const GTestResponse* testrsp = dynamic_cast<const GTestResponse*>(&rsp);
        if (testrsp == NULL) {
            throw;
        }

         // Delete old response function
        if (m_response != NULL) delete m_response;

        // Clone response function
        m_response = testrsp->clone();

        // Return
        return;
    }
    virtual GTestResponse*    response(void) const{return m_response;}
    virtual GTestPointing*    pointing(void) const{ return m_pointing;}
    virtual std::string      instrument(void) const { return m_instrument; }
    virtual double           ontime(void) const { return m_ontime; }
    virtual double           livetime(void) const { return m_ontime; }
    virtual double           deadc(const GTime& time) const { return 1.0; }
    virtual void             read(const GXmlElement& xml){ return; }
    virtual void             write(GXmlElement& xml) const{ return; }
    void ontime(const double& ontime) { m_ontime=ontime; }
    virtual std::string      print(const GChatter& chatter = NORMAL) const{
       // Initialise result string
        std::string result;
        
        // Append header
        result.append("=== GTestObservation ===");
        result.append("\n"+parformat("Instrument Name")+m_instrument);
        
        GTestEventList* list = dynamic_cast<GTestEventList*>(m_events);
        if(list!=NULL){
            result.append("\n"+parformat("Number of events")+str(list->number()));
            result.append("\n"+parformat("GTI")+list->gti().print());
        }
        
        return result;
    }


protected:
    // Protected methods
    void init_members(void){
        m_instrument="Test Instrument";
        m_response = new GTestResponse();
        m_pointing = new GTestPointing();
        m_ontime = 0.0;
        return;
    }
    void copy_members(const GTestObservation& obs){
        m_ontime=obs.m_ontime;
        m_response = new GTestResponse(*obs.m_response);
        m_pointing = new GTestPointing(*obs.m_pointing);
        m_instrument = obs.m_instrument;
        return;    
    }
    void free_members(void){ return; }

    // Protected members
    std::string   m_instrument;   //!< Instrument name
    GTestResponse* m_response;     //!< Pointer to response functions
    GTestPointing* m_pointing;     //!< Pointer to pointing direction
    double         m_ontime;       //!< ontime
};

#endif /* GTestOBSERVATION_HPP */
