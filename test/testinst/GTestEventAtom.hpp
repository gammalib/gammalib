/***************************************************************************
 *                 GTestEventAtom.hpp - Test event atom class              *
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
/**
 * @file GTestEventAtom.hpp
 * @brief Test data atom definition
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTEVENTATOM_HPP
#define GTESTEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventAtom.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GTestInstDir.hpp"
#include "GTools.hpp"

class GTestEventAtom : public GEventAtom {

    // Friend classes
    friend class GTestEventList;

    public:
    // Constructors and destructors
        GTestEventAtom(void){
            init_members();
            return;
        }
        GTestEventAtom(const GTestEventAtom& atom){
            init_members();
            copy_members(atom);
            return;
        }
        
        virtual ~GTestEventAtom(void){
            free_members();
            return;
        }

    // Operators
        GTestEventAtom& operator= (const GTestEventAtom& atom){
            // Execute only if object is not identical
            if (this != &atom) {

             // Copy base class members
               this->GEventAtom::operator=(atom);

                // Free members
                free_members();

                // Initialise private members for clean destruction
                init_members();

                // Copy members
                copy_members(atom);

            } // endif: object was not identical

            // Return this object
            return *this;
        };

    // Implemented pure virtual base class methods
        void clear(void){
            // Free class members (base and derived classes, derived class first)
            free_members();
            this->GEventAtom::free_members();
            this->GEvent::free_members();

            // Initialise members
            this->GEvent::init_members();
            this->GEventAtom::init_members();
            init_members();
            
            return;
        }
        
        GTestEventAtom*     clone(void) const{ return new GTestEventAtom(*this);}
        const GTestInstDir& dir(void) const { return m_dir; }
        const GEnergy&      energy(void) const { return m_energy; }
        const GTime&        time(void) const { return m_time; }
        void                dir(const GTestInstDir& dir) { m_dir=dir; }
        void                energy(const GEnergy& energy) { m_energy=energy; }
        void                time(const GTime& time) { m_time=time; }
        std::string         print(const GChatter& chatter = NORMAL) const{
            std::string result("== GTestEventAtom == \n");
            result.append("Direction : "+m_dir.print()+"\n");
            result.append("Time : "+gammalib::str(m_time.secs())+"\n");
            result.append("Energy : "+m_energy.print()+"\n");
            return result;
        }

    protected:
    // Protected methods
        void init_members(void){
            m_dir.clear();
            m_energy.clear();
            m_time.clear();
            return;
        }
        
        void copy_members(const GTestEventAtom& atom){
            m_dir         = atom.dir();
            m_time        = atom.time();
            m_energy      = atom.energy();
            return;
        }
        
        void free_members(void){ return; }

    // Protected members
        GTestInstDir   m_dir;            //!< Event direction
        GEnergy       m_energy;         //!< Event energy
        GTime         m_time;           //!< Event time
};

#endif /* GTestEVENTATOM_HPP */
