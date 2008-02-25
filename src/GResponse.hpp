/***************************************************************************
 *              GResponse.hpp  -  Response abstract base class             *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GRESPONSE_HPP
#define GRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
//#include "GVector.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                          GResponse class definition                     *
 ***************************************************************************/
class GResponse {

// Public methods
public:
    // Constructors and destructors
    //GResponse();
    //GResponse(const GResponse& r);
    //virtual ~GResponse();

    // Operators

    // Methods
    virtual void load(std::string rspname) = 0;
    virtual void save(std::string rspname) = 0;
    virtual void set_caldb(std::string dbname) = 0;
  
// Methods and data that are available to derived classes
protected:
    // Protected methods

    // Protected data area
    std::string m_caldb;      // Calibration database name
    std::string m_rspname;    // Response name

// Methods that are available to the base class only
private:
};

#endif /* GRESPONSE_HPP */
