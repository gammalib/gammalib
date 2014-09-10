/***************************************************************************
 *                GTestRoi.hpp - Test region of interest class             *
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

#ifndef GTestROI_HPP
#define GTestROI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GRoi.hpp"
#include "GTestInstDir.hpp"


class GTestRoi : public GRoi {

public:
    // Constructors and destructors
    GTestRoi(void){ 
        init_members();
        return;
    }
    GTestRoi(const GTestRoi& roi){
        init_members();
        copy_members(roi);
        return;
    }
    
    virtual ~GTestRoi(void){
        free_members();
        return;
    }

    // Operators
    GTestRoi& operator=(const GTestRoi& roi){
        if (this != &roi) {
            this->GRoi::operator=(roi);
            free_members();
            init_members();
            copy_members(roi);
        }
        return *this;
    }

    // Implemented pure virtual base class methods
    virtual void clear(void){
        free_members();
        this->GRoi::free_members();
        this->GRoi::init_members();
        init_members();
        return;
    }
    virtual GTestRoi* clone(void) const{
        return new GTestRoi(*this);
    }
    virtual std::string classname(void) const { return "GTestRoi"; }
    virtual bool contains(const GEvent& event) const {
        return true;
    }
    virtual std::string  print(const GChatter& chatter = NORMAL) const{
        return "=== GTestRoi ===";
    }

protected:
    // Protected methods
    void init_members(void){ return; }
    void copy_members(const GTestRoi& roi) { return; }
    void free_members(void){ return; }

};

#endif /* GTestROI_HPP */
