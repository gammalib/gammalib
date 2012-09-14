/***************************************************************************
 *            GTestEventCube.hpp  -  Test event bin container class        *
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
 
#ifndef GTestEVENTCUBE_HPP
#define GTestEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventCube.hpp"
#include "GTestEventBin.hpp"
#include "GSkymap.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GTestInstDir.hpp"
#include "GFitsTable.hpp"
#include "GFitsImage.hpp"

class GTestEventCube : public GEventCube {

public:
    // Constructors and destructors
    GTestEventCube(void){
         // Initialise members
        init_members();

         // Return
        return;
    }
    GTestEventCube(const GTestEventCube& cube) : GEventCube(cube){
        // Initialise members
        init_members();

        // Copy members
        copy_members(cube);

        // Return
        return;
    }
    virtual ~GTestEventCube(void){
        // Free members
        free_members();

        // Return
        return;
    }

    // Operators
    virtual GTestEventCube&      operator=(const GTestEventCube& cube){
        // Execute only if object is not identical
        if (this != &cube) {

            // Copy base class members
            this->GEventCube::operator=(cube);

            // Free members
            free_members();

            // Initialise members
            init_members();

            // Copy members
            copy_members(cube);

        } // endif: object was not identical

        // Return this object
        return *this;
    }
    
    virtual GTestEventBin* operator[](const int& index){
        #if defined(G_RANGE_CHECK)
            if (index < 0 || index >= size())
            throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
        #endif
        // Return pointer
        return &m_bins.at(index);
    }
    
    virtual const GTestEventBin* operator[](const int& index) const{
        #if defined(G_RANGE_CHECK)
            if (index < 0 || index >= size())
            throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
        #endif
        // Return pointer
        return &m_bins.at(index);
    }

    // Implemented pure virtual base class methods
    virtual void clear(void){
        // Free class members (base and derived classes, derived class first)
        free_members();
        this->GEventCube::free_members();
        this->GEvents::free_members();

        // Initialise members
        this->GEvents::init_members();
        this->GEventCube::init_members();
        init_members();

        // Return
        return;
    }
    virtual GTestEventCube* clone(void) const{ return new GTestEventCube(*this);}
    virtual int            size(void) const{ return m_bins.size(); }
    virtual int            dim(void) const{ return 1;}
    virtual int            naxis(int axis) const{ return 1;}
    virtual void           load(const std::string& filename){ return; }
    virtual void           save(const std::string& filename, bool clobber = false) const{ return; }
    virtual void           read(const GFits& file){ return; }
    virtual void           write(GFits& file) const{ return; }
    virtual int            number(void) const{ return m_counts; }
    virtual std::string    print(void) const{
        // Initialise result string
        std::string result;

        // Append header
        result.append("=== GTestEventCube ===");
        result.append("\n"+parformat("Number of events")+str(number()));
        result.append("\n"+parformat("Number of elements")+str(size()));

        // Return result
        return result;
    }

    // Other methods
    void append(GTestEventBin& bin){
        m_bins.push_back(bin);
        m_counts+=(int)bin.counts();
    }

protected:
    // Protected methods
    void         init_members(void){
        // Initialise members
        m_bins.clear();
        
        m_counts=0;

        // Return
        return;
    }
    
    void         copy_members(const GTestEventCube& cube){
        // Copy members
        m_bins      = cube.m_bins;
        m_counts = cube.m_counts;
        
        // Return
        return;
    }
    
    void         free_members(void){ return; }
    virtual void set_energies(void){}
    virtual void set_times(void){}
    
    // Protected members
    std::vector<GTestEventBin> m_bins;
    int                       m_counts; //!< Number of event.
};

#endif /* GTestEVENTCUBE_HPP */
