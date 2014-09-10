/***************************************************************************
 *                 GTestEventBin.hpp  -  Test event bin class              *
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

#ifndef GTESTEVENTBIN_HPP
#define GTESTEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GTestInstDir.hpp"

class GTestEventBin : public GEventBin {

    // Friend classes
    friend class GTestEventCube;

public:
    // Constructors and destructors
    GTestEventBin(void) : GEventBin(){
        init_members();
        return;
    }
    GTestEventBin(const GTestEventBin& bin) : GEventBin(bin){
        init_members();
        copy_members(bin);
        return;
    }
    virtual ~GTestEventBin(void){
        free_members();
        return; 
    }

    // Operators
    virtual GTestEventBin& operator= (const GTestEventBin& bin){
            // Execute only if object is not identical
        if (this != &bin) {

            // Copy base class members
            this->GEventBin::operator=(bin);

            // Free members
            free_members();

            // Initialise private members for clean destruction
            init_members();

            // Copy members
            copy_members(bin);

        } // endif: object was not identical

        // Return this object
        return *this;
    }

    // Implemented pure virtual base class methods
    
    virtual void clear(void){
        // Free class members (base and derived classes, derived class first)
        free_members();
        this->GEventBin::free_members();
        this->GEvent::free_members();

        // Initialise members
        this->GEvent::init_members();
        this->GEventBin::init_members();
        init_members();

        // Return
        return;
    }
    virtual GTestEventBin* clone(void) const{ return new GTestEventBin(*this); }
    virtual std::string    classname(void) const { return "GTestEventBin"; }
    virtual double         size(void) const{
        // Compute bin size
        double size = ewidth().MeV() * ontime();

        // Return bin size
        return size;
    }
    virtual const GTestInstDir& dir(void) const { return m_dir; }
    virtual const GEnergy&     energy(void) const { return m_energy; }
    virtual const GTime&       time(void) const { return m_time; }
    virtual double             counts(void) const { return m_counts; }
    virtual double             error(void) const{
        // Compute uncertainty
        double error = sqrt(counts()+1.0e-50);

        // Return error
        return error;
    }
    virtual void               counts(const double& counts) { m_counts=counts; }
    virtual std::string        print(const GChatter& chatter = NORMAL) const{
         // Initialise result string
        std::string result;

        // Append number of counts
        result.append(gammalib::str(counts()));

        // Return result
        return result;
    }

    // Other methods
    const GEnergy ewidth(void) const { return m_ewidth; }
    const double& ontime(void) const { return m_ontime; }
    
    void energy(GEnergy energy){ m_energy=energy;}
    void dir(GTestInstDir dir){ m_dir=dir;}
    void time(GTime time){ m_time=time;}
    void ewidth(GEnergy ewidth){ m_ewidth=ewidth;}
    void ontime(double ontime){ m_ontime=ontime;}
protected:
    // Protected methods
    void init_members(void){
        // Initialise members
        m_counts=0;
        m_ontime=0;
        
        // Return
        return; 
    }
    void copy_members(const GTestEventBin& bin){
        // Copy members
        m_energy = bin.m_energy;
        m_dir    = bin.m_dir;
        m_time   = bin.m_time;
        m_counts = bin.m_counts;
        m_ewidth = bin.m_ewidth;
        m_ontime = bin.m_ontime;

        // Return
        return;
    }
    void free_members(void){ return; }

    // Protected members
    GEnergy      m_energy;      //!< Pointer to bin energy
    GTestInstDir m_dir;         //!< Pointer to bin direction
    GTime        m_time;        //!< Pointer to bin time
    double       m_counts;      //!< Pointer to number of counts
    GEnergy      m_ewidth;      //!< Pointer to energy width of bin
    double       m_ontime;      //!< Pointer to ontime of bin (seconds)
};

#endif /* GTestEVENTBIN_HPP */
