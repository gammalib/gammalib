/***************************************************************************
 *             GTestEventCube.hpp - Test event bin container class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Jean-Baptiste Cayrou                        *
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
 * @file GTestEventCube.hpp
 * @brief Event cube test class
 * @author Jean-Baptiste Cayrou
 */

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


/***********************************************************************//**
 * @class GTestEventCube
 *
 * @brief Event cube test class
 *
 * This class implements an event cube for unit testing.
 ***************************************************************************/
class GTestEventCube : public GEventCube {

public:

    // Constructors and destructors
    GTestEventCube(void){
        init_members();
        return;
    }
    GTestEventCube(const GTestEventCube& cube) : GEventCube(cube){
        init_members();
        copy_members(cube);
        return;
    }
    virtual ~GTestEventCube(void){
        free_members();
        return;
    }

    // Operators
    virtual GTestEventCube& operator=(const GTestEventCube& cube){
        if (this != &cube) {
            this->GEventCube::operator=(cube);
            free_members();
            init_members();
            copy_members(cube);
        } // endif: object was not identical
        return *this;
    }
    virtual GTestEventBin* operator[](const int& index){
        return &m_bins.at(index);
    }
    
    virtual const GTestEventBin* operator[](const int& index) const{
        return &m_bins.at(index);
    }

    // Implemented pure virtual base class methods
    virtual void clear(void){
        free_members();
        this->GEventCube::free_members();
        this->GEvents::free_members();
        this->GEvents::init_members();
        this->GEventCube::init_members();
        init_members();
        return;
    }
    virtual GTestEventCube* clone(void) const { return new GTestEventCube(*this);}
    virtual std::string     classname(void) const { return "GTestEventCube"; }
    virtual int            size(void) const { return m_bins.size(); }
    virtual int            dim(void) const { return 1; }
    virtual int            naxis(const int& axis) const { return 1;}
    virtual void           load(const std::string& filename) { }
    virtual void           save(const std::string& filename,
                                const bool& clobber = false) const { }
    virtual void           read(const GFits& file) { }
    virtual void           write(GFits& file) const { }
    virtual int            number(void) const { return m_counts; }
    virtual std::string    print(const GChatter& chatter = NORMAL) const {
        std::string result;
        result.append("=== GTestEventCube ===");
        result.append("\n"+gammalib::parformat("Number of events")+gammalib::str(number()));
        result.append("\n"+gammalib::parformat("Number of elements")+gammalib::str(size()));
        return result;
    }

    // Other methods
    void append(GTestEventBin& bin){
        m_bins.push_back(bin);
        m_counts+=(int)bin.counts();
    }

protected:
    // Protected methods
    void init_members(void){
        m_bins.clear();
        m_counts = 0;
        return;
    }
    
    void copy_members(const GTestEventCube& cube){
        m_bins   = cube.m_bins;
        m_counts = cube.m_counts;
        return;
    }
    
    void         free_members(void) {}
    virtual void set_energies(void) {}
    virtual void set_times(void) {}
    
    // Protected members
    std::vector<GTestEventBin> m_bins;
    int                        m_counts; //!< Number of counts
};

#endif /* GTestEVENTCUBE_HPP */
