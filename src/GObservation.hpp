/***************************************************************************
 *           GObservation.hpp  -  Observation abstract base class          *
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

#ifndef GOBSERVATION_HPP
#define GOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                        GObservation class definition                    *
 ***************************************************************************/
class GObservation {
  // Friend classes
  friend class GResponse;
  friend class GData;

// Public methods
public:
    // Constructors and destructors
    //GObservation();
    //GObservation(const GObservation& o);
    //virtual ~GObservation();

    // Operators

    // Methods
  
// Methods and data that are available to derived classes
protected:
    // Protected methods

    // Protected data area
    std::string m_obsname;    // Name of observation
	GResponse*  m_response;   // Response associated with observation
	GData*      m_data;       // Data associated with observations

// Methods that are available to the base class only
private:
};

#endif /* GOBSERVATION_HPP */
