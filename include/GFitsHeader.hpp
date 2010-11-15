/***************************************************************************
 *               GFitsHeader.hpp  - FITS header handling class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsHeader.hpp
 * @brief GFitsHeader class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSHEADER_HPP
#define GFITSHEADER_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHeaderCard.hpp"


/***********************************************************************//**
 * @class GFitsHeader
 *
 * @brief Interface for FITS header class
 *
 * The FITS header class contains all cards that are found in the header of
 * a HDU. All cards will be hold in memory, so no link to a FITS file is
 * required. Cards may be read from one file (using the 'open' method) and
 * saved into another file (using the 'save' method). Cards are added or
 * changed using the 'update' method.
 ***************************************************************************/
class GFitsHeader {

    // Friend classes
    friend class GFitsHDU;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsHeader& header);

public:
    // Constructors and destructors
    GFitsHeader(void);
    GFitsHeader(const GFitsHeader& header);
    virtual ~GFitsHeader(void);

    // Operators
    GFitsHeader& operator= (const GFitsHeader& header);

    // Methods
    void             clear(void);
    int              size(void) const;
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
    int              num_cards(void) const;

private:
    // Private methods
    void             init_members(void);
    void             copy_members(const GFitsHeader& header);
    void             free_members(void);
    GFitsHeaderCard* card_ptr(const std::string& keyname);
    void             open(void* vptr);
    void             save(void* vptr);
    void             close(void);

    // Private data area
    int              m_num_cards;
    GFitsHeaderCard* m_card;
};

#endif /* GFITSHEADER_HPP */
