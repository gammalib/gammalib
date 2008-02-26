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
#include "GFitsHeaderCard.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                       GFitsHeader class definition                      *
 ***************************************************************************/
class GFitsHeader {

// Public methods
public:
    // Constructors and destructors
    GFitsHeader();
    GFitsHeader(const GFitsHeader& header);
    ~GFitsHeader();

    // Operators
    GFitsHeader& operator= (const GFitsHeader& header);

    // Methods
    void             open(__fitsfile*  fptr);
    GFitsHeaderCard* card(const std::string keyname);
    GFitsHeaderCard* card(const int cardno);
    std::string      value(const std::string keyname);
    
private:
    // Private methods
    void             init_members(void);
    void             copy_members(const GFitsHeader& header);
    void             free_members(void);
    GFitsHeaderCard* GFitsHeader::card_ptr(const std::string keyname);

    // Private data area
    int              m_num_cards;
    GFitsHeaderCard* m_card;
};

#endif /* GFITSHEADER_HPP */
