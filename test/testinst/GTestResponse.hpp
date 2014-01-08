/***************************************************************************
 *          GTestResponse.hpp  -  Test response class                      *
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
 * @file GTestResponse.hpp
 * @brief Test response class interface definition
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTRESPONSE_HPP
#define GTESTRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GResponse.hpp"


/***********************************************************************//**
 * @class GTestResponse
 *
 * @brief Test response class
 *
 ***************************************************************************/
class GTestResponse : public GResponse {

public:
    // Constructors and destructors
    GTestResponse(void){
        init_members();
        return;
    }
    
    GTestResponse(const GTestResponse& rsp){
        init_members();
        copy_members(rsp);
        return;
    }
    
    virtual ~GTestResponse(void){
        free_members();
        return;
    }

    // Operators
    virtual GTestResponse& operator=(const GTestResponse& rsp)
    {
        // Execute only if object is not identical
        if (this != &rsp) {

            // Copy base class members
            this->GResponse::operator=(rsp);

            // Free members
            free_members();

            // Initialise private members
            init_members();

            // Copy members
            copy_members(rsp);

        } // endif: object was not identical

        // Return this object
        return *this;
    }

    // Implemented pure virtual methods
    
    virtual void clear(void){
        // Free members
        free_members();
        this->GResponse::free_members();

        // Initialise private members
        this->GResponse::init_members();
        init_members();
        
        return;
    }
    
    virtual GTestResponse* clone(void) const{
        return new GTestResponse(*this);
    }
    
    virtual bool          has_edisp(void) const { return false; }
    virtual bool          has_tdisp(void) const { return false; }
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const { return 1.0; }
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const { return 1.0; }
    virtual std::string   print(const GChatter& chatter = NORMAL) const{ return "=== GTestReponse ==="; }

protected:
    // Protected methods
    void init_members(void){ return; }
    void copy_members(const GTestResponse& pnt){ return; }
    void free_members(void){ return; }
};

#endif /* GTESTRESPONSE_HPP */
