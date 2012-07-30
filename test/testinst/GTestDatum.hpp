/***************************************************************************
 *         GTestDatum.hpp  -  Test spectral point class                    *
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

#ifndef GTestDATUM_HPP
#define GTestDATUM_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GTestInstDir.hpp"

class GTestDatum : public GEventBin {

    // Friend classes
    friend class GTestSpectrum;

public:
    // Constructors and destructors
    GTestDatum(void){
        init_members();
        return;
    }
    GTestDatum(const GTestDatum& datum){
        init_members();
        copy_members(datum);
        return;
    }
    virtual ~GTestDatum(void){
        free_members();
        return;
    }

    // Operators
    virtual GTestDatum& operator= (const GTestDatum& datum){
        // Execute only if object is not identical
        if (this != &datum) {

            // Copy base class members
            this->GEventBin::operator=(datum);

            // Free members
            free_members();

            // Initialise private members for clean destruction
            init_members();

            // Copy members
            copy_members(datum);

        } // endif: object was not identical

        // Return this object
        return *this;
    }

    // Implemented pure virtual base class methods
    virtual void clear(void){     // Free class members (base and derived classes, derived class first)
        free_members();
        this->GEventBin::free_members();
        this->GEvent::free_members();

    // Initialise members
        this->GEvent::init_members();
        this->GEventBin::init_members();
        init_members();

    // Return
        return;
    };
    virtual GTestDatum*        clone(void) const{ return new GTestDatum(*this); }
    virtual double             size(void) const { return 1.0; }
    virtual const GInstDir&    dir(void) const { return m_dir; }
    virtual const GEnergy&     energy(void) const { return m_eng; }
    virtual const GTime&       time(void) const { return m_time; }
    virtual double             counts(void) const { return 1; }
    virtual double             error(void) const { return 0; }
    virtual void               counts(const double& counts) { return; }
    virtual std::string        print(void) const { return "=== GTestDatum ===";}

protected:
    // Protected methods
    void init_members(void){
        m_dir.clear();
        m_time.clear();
        m_eng.clear();
        return;
    }
    void copy_members(const GTestDatum& datum){
        m_dir      = datum.m_dir;
        m_time     = datum.m_time;
        m_eng      = datum.m_eng;
        return;
    }
    void free_members(void){ return; }

    // Protected members
    GTestInstDir m_dir;       //!< Instrument direction of spectral point (not used)
    GTime       m_time;      //!< Time of spectral point (not used)
    GEnergy     m_eng;       //!< Energy of spectral point

};

#endif /* GTestDATUM_HPP */
