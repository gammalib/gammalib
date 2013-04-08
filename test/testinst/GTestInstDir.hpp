/***************************************************************************
 *     GTestInstDir.hpp  -  Test instrument direction class                *
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

#ifndef GTESTINSTDIR_HPP
#define GTESTINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"


class GTestInstDir : public GInstDir {

public:
    // Constructors and destructors
    GTestInstDir(void){
        init_members();
        return;
    }
        
    GTestInstDir(const GTestInstDir& dir) : GInstDir(dir){
        init_members();
        // Copy members
        copy_members(dir);
        return;
    }
    
    virtual ~GTestInstDir(void){
        free_members();
        return;
    }

    // Operators
    GTestInstDir& operator= (const GTestInstDir& dir){
         // Execute only if object is not identical
        if (this != &dir) {

            // Copy base class members
            this->GInstDir::operator=(dir);

            // Free members
            free_members();

            // Initialise private members
            init_members();

            // Copy members
            copy_members(dir);

        } // endif: object was not identical

        // Return this object
        return *this;   
    }

    // Methods
    void clear(void){
        // Free members
        free_members();
        this->GInstDir::free_members();

        // Initialise private members
        this->GInstDir::init_members();
        init_members();
        
        return;
    }
    
    GTestInstDir* clone(void) const{ return new GTestInstDir(*this); }
        
    std::string  print(const GChatter& chatter = NORMAL) const{ 
        std::string result = "=== GTestInstDir ===";
            return result;
        }
        

protected:
    // Protected methods
    void init_members(void){ return; }
    void copy_members(const GTestInstDir& dir){ return; }
    void free_members(void){ return; }
};

#endif /* GTesINSTDIR_HPP */
