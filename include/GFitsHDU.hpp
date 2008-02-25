/***************************************************************************
 *                  GFitsHDU.hpp  - FITS HDU handling class                *
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

#ifndef GFITSHDU_HPP
#define GFITSHDU_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHeader.hpp"
#include "GFitsData.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                          GFitsHDU class definition                      *
 ***************************************************************************/
class GFitsHDU {

// Public methods
public:
    // Constructors and destructors
    //GFitsHDU();
    //GFitsHDU(const GFitsHDU& hdu);
    //virtual ~GFitsHDU();

    // Operators

    // Methods
    GFitsHeader* header(void) { return m_header; }
    GFitsData*   data(void) { return m_data; }
    
// Methods and data that are available to derived classes
protected:
    // Protected methods

    // Protected data area
    GFitsHeader* m_header;
    GFitsData*   m_data;

// Methods that are available to the base class only
private:
};

#endif /* GFITSHDU_HPP */
