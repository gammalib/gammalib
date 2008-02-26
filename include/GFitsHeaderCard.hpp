/***************************************************************************
 *       GFitsHeaderCard.hpp  - FITS header card abstract base class       *
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

#ifndef GFITSHEADERCARD_HPP
#define GFITSHEADERCARD_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                     GFitsHeaderCard class definition                    *
 ***************************************************************************/
class GFitsHeaderCard {
public:
    // Constructors & Destructors
    GFitsHeaderCard();
    GFitsHeaderCard(const GFitsHeaderCard& card);
    ~GFitsHeaderCard();
    
    // Operators
    GFitsHeaderCard& operator= (const GFitsHeaderCard& card);
    
    // Methods
    std::string  keyname(void);
    std::string  value(void);
    std::string  unit(void);
    std::string  comment(void);
    void         set_keyname(std::string& keyname);
    void         set_value(std::string& value);
    void         set_unit(std::string& unit);
    void         set_comment(std::string& comment);
    void         read(__fitsfile* fptr, int keynum);
    void         read(__fitsfile* fptr, std::string keyname);
    void         write(__fitsfile* fptr);
    
private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsHeaderCard& card);
    void free_members(void);

    // Private data area
    std::string m_keyname;
    std::string m_value;
    std::string m_unit;
    std::string m_comment;
};

#endif /* GFITSHEADERCARD_HPP */
