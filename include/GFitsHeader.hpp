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
#include "GFitsCfitsio.hpp"
#include "GFitsHeaderCard.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                       GFitsHeader class definition                      *
 ***************************************************************************/
class GFitsHeader {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsHeader& header);

public:
    // Constructors and destructors
    GFitsHeader();
    GFitsHeader(const GFitsHeader& header);
    ~GFitsHeader();

    // Operators
    GFitsHeader& operator= (const GFitsHeader& header);

    // Methods
    void             open(__fitsfile* fptr);
    void             close(void);
    void             update(const GFitsHeaderCard& card);
    GFitsHeaderCard* card(const std::string& keyname);
    GFitsHeaderCard* card(const int& cardno);
    std::string      string(const std::string& keyname);
    std::string      string(const int& cardno);
    double           real(const std::string& keyname);
    double           real(const int& cardno);
    int              integer(const std::string& keyname);
    int              integer(const int& cardno);
    GFitsHeader*     clone(void) const;
    
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


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline 
GFitsHeaderCard* GFitsHeader::card(const std::string& keyname)
{
    return GFitsHeader::card_ptr(keyname);
}
inline
GFitsHeaderCard* GFitsHeader::card(const int& cardno)
{
    return (cardno >= 0 && cardno < m_num_cards) ? &(m_card[cardno]) : NULL;
}
inline 
GFitsHeader* GFitsHeader::clone(void) const 
{
    return new GFitsHeader(*this);
}

#endif /* GFITSHEADER_HPP */
