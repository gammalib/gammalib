/***************************************************************************
 *               GFitsHeader.hpp  - FITS header handling class             *
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

#ifndef GFITSHEADER_HPP
#define GFITSHEADER_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>                    // for std::vector
#include "GFitsHeaderCard.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                       GFitsHeader class definition                      *
 ***************************************************************************/
class GFitsHeader {

// Public methods
public:
    // Constructors and destructors
    //GFitsHeader();
    //GFitsHeader(const GFitsHeader& hdr);
    //virtual ~GFitsHeader();

    // Operators

    // Methods
    //GFitsHeaderCard* card(const std::string keyword);
    //GFitsHeaderCard* card(const int cardno);
    
// Methods and data that are available to derived classes
protected:
    // Protected methods

    // Protected data area
    std::vector<GFitsHeaderCard> m_card;     // List of FITS header cards

// Methods that are available to the base class only
private:
};

#endif /* GFITSHEADER_HPP */
