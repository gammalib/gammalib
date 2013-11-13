/***************************************************************************
 *            GTestModelData.hpp - Test data model class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Jean-Baptiste Cayrou                        *
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
 * @file GTestModelData.hpp
 * @brief Test data model definition
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTMODELDATA_HPP
#define GTESTMODELDATA_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelData.hpp"
#include "GModelPar.hpp"
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GXmlElement.hpp"
#include "GTestEventList.hpp"
#include "GTestEventAtom.hpp"
#include "GTestEventCube.hpp"
#include "GTestEventBin.hpp"


/***********************************************************************//**
 * @class GTestModelData
 *
 * @brief Test data model
 ***************************************************************************/
class GTestModelData : public GModelData {

public:
    // Constructors and destructors
    GTestModelData(void) : GModelData() {
        init_members();
        set_pointers();
        return;
    }
    explicit GTestModelData(const GXmlElement& xml) : GModelData(xml) {
        init_members();
        set_pointers();
        return;
    }
    GTestModelData(const GTestModelData& model) : GModelData(model) { 
        copy_members(model);
        set_pointers();
        return;        
    }
    virtual ~GTestModelData(void) {
        free_members();
        return;
    }

    // Operators
    virtual GTestModelData& operator=(const GTestModelData& model) {
        if (this != &model) {
            this->GModelData::operator=(model);
            free_members();
            copy_members(model);
            set_pointers();
        }
        return *this;
    }

    // Implement Pure virtual methods
    virtual void clear(void){
        free_members();
        init_members();
        return;
    }
    virtual GTestModelData* clone(void) const { return new GTestModelData(*this); }
    virtual std::string type(void) const {return "=== GTestModelData ===";}
    virtual bool        isconstant(void) const {return true;}
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const { 
                                 double result = m_modelTps->eval(event.time());
                                 return result;
                             }
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const {
                                           double result = m_modelTps->eval_gradients(event.time());
                                           return result;
                                       }
    virtual double      npred(const GEnergy& obsEng, const GTime& obsTime,
                              const GObservation& obs) const {
                                  double result = m_pars[0]->value();
                                  return result;
                              }

   // Generate an EventList, rate is the number of event per second. Events have a time between tmin and tmax.
   virtual GTestEventList* generateList(const double &rate, const GTime &tmin, const GTime &tmax, GRan &ran)
   {
        // Create an event list
        GTestEventList * list = new GTestEventList();
        
        // Set min and max energy for ebounds
        // npred method integrate the model on time and energy.
        // In order to have a rate which not depend on energy we create an interval of 1 Mev.
        GEnergy engmin,engmax;
        engmin.MeV(1.0);
        engmax.MeV(2.0);
        
        // Instrument Direction
        GTestInstDir dir;
        
        // Generate an times list.
        GTimes times = m_modelTps->mc(rate,tmin, tmax,ran);
        
        for (int i = 0; i < times.size() ; ++i) {
            GTestEventAtom event;
            event.dir(dir);
            event.energy(engmin);
            event.time(times[i]);
            
            // Add the event to the list
            list->append(event);
        }
        
        // Create a time interval and add it to the list.
        GGti gti;
        gti.append(tmin,tmax);
        list->gti(gti);
        
        // Create an energy interval and add it to the list
        GEbounds ebounds;
        ebounds.append(engmin,engmax);
        list->ebounds(ebounds);
        
        return list;
    }
    
    // Generate an EventCube, rate is the number of event per second. Events have a time between tmin and tmax.
    virtual GTestEventCube* generateCube(const double &rate, const GTime &tmin, const GTime &tmax, GRan &ran)
    {
        
        // Create an event list
        GTestEventCube* cube = new GTestEventCube();
        
        // Set min and max energy for ebounds
        // npred method integrate the model on time and energy.
        // In order to have a rate which not depend on energy we create an interval of 1 Mev.
        GEnergy engmin,engmax;
        engmin.MeV(1.0);
        engmax.MeV(2.0);
        
        // Instrument Direction
        GTestInstDir dir;
        
        // Generate an times list.
        GTimes times = m_modelTps->mc(rate, tmin, tmax, ran);

        GTestEventBin bin;
        bin.time(times[0]);
        bin.energy(engmin);
        bin.ewidth(engmax-engmin);
        bin.dir(dir);
        bin.ontime(10); // 10 sec per bin

        for (int i = 0; i < times.size(); ++i) {
            if ((bin.time().secs() + bin.ontime()) < times[i].secs()) {

                // Add the event to the cube
                cube->append(bin);
                
                bin.counts(0.0);
                bin.time(times[i]);
                bin.energy(engmin);
                bin.ewidth(engmax-engmin);
                bin.dir(dir);
                bin.ontime(10); // 10 sec per bin
            }
  
            bin.counts(bin.counts()+1);

        }
        
        
        // Create a time interval and add it to the list.
        GGti gti;
        gti.append(tmin,tmax);
        cube->gti(gti);
        
        // Create an energy interval and add it to the list
        GEbounds ebounds;
        ebounds.append(engmin,engmax);
        cube->ebounds(ebounds);
        
        return cube;
    };
    
    
    virtual void        read(const GXmlElement& xml){ return; }
    virtual void        write(GXmlElement& xml) const { return; }
    virtual std::string print(const GChatter& chatter = NORMAL) const { 
        std::string result;
        
        result.append(gammalib::parformat("Number of temporal par's")+gammalib::str(m_modelTps->size()));
        for (int i = 0; i < m_modelTps->size(); ++i) {
            result.append("\n"+(((*m_modelTps)[i]).print()));
        }
            
        return result; }
    
    virtual GModelTemporalConst* temporal() const { return m_modelTps;}
    
protected:
    // Protected methods
    void init_members(void) {
        m_modelTps = new GModelTemporalConst();
    }
    void copy_members(const GTestModelData& model){
        m_modelTps = model.temporal()->clone();
    }
    void free_members(void) {
        delete m_modelTps;
    }
    void set_pointers(void){
        m_pars.clear();
        for(int i = 0; i < m_modelTps->size(); ++i) {
            (*m_modelTps)[i].free(); // Set free
            m_pars.push_back(&((*m_modelTps)[i]));
        }
    }

    // Protected members
    GModelTemporalConst* m_modelTps; //!< Temporal model
};

#endif /* GMODELDATA_HPP */
