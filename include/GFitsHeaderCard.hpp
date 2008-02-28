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
    int          value_type(void);
    std::string  string(void);
    double       real(void);
    int          integer(void);
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
    int  get_value_type(void);

    // Private data area
    std::string m_keyname;
    std::string m_value;
    int         m_value_type;
    std::string m_unit;
    std::string m_comment;
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline std::string GFitsHeaderCard::keyname(void) { return m_keyname; }
inline std::string GFitsHeaderCard::value(void) { return m_value; }
inline int         GFitsHeaderCard::value_type(void) { return m_value_type; }
inline std::string GFitsHeaderCard::unit(void) { return m_unit; }
inline std::string GFitsHeaderCard::comment(void) { return m_comment; }
inline void        GFitsHeaderCard::set_keyname(std::string& keyname) { m_keyname = keyname; }
inline void        GFitsHeaderCard::set_value(std::string& value) { m_value = value; }
inline void        GFitsHeaderCard::set_unit(std::string& unit) { m_unit = unit; }
inline void        GFitsHeaderCard::set_comment(std::string& comment) { m_comment = comment; }


#endif /* GFITSHEADERCARD_HPP */
