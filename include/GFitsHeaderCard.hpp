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

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsHeaderCard& card);

public:
    // Constructors & Destructors
    GFitsHeaderCard();
    GFitsHeaderCard(const std::string& keyname, const std::string& value, 
                    const std::string& comment);
    GFitsHeaderCard(const GFitsHeaderCard& card);
    ~GFitsHeaderCard();

    // Operators
    GFitsHeaderCard& operator= (const GFitsHeaderCard& card);

    // Methods
    std::string  keyname(void) const;
    std::string  value(void) const;
    int          value_type(void) const;
    std::string  unit(void) const;
    std::string  comment(void) const;
    std::string  string(void);
    double       real(void);
    int          integer(void);
    void         set_keyname(const std::string& keyname);
    void         set_value(const std::string& value);
    void         set_unit(const std::string& unit);
    void         set_comment(const std::string& comment);
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
inline std::string GFitsHeaderCard::keyname(void) const { return m_keyname; }
inline std::string GFitsHeaderCard::value(void) const { return m_value; }
inline int         GFitsHeaderCard::value_type(void) const { return m_value_type; }
inline std::string GFitsHeaderCard::unit(void) const { return m_unit; }
inline std::string GFitsHeaderCard::comment(void) const { return m_comment; }
inline void        GFitsHeaderCard::set_keyname(const std::string& keyname) { m_keyname = keyname; }
inline void        GFitsHeaderCard::set_value(const std::string& value) { m_value = value; }
inline void        GFitsHeaderCard::set_unit(const std::string& unit) { m_unit = unit; }
inline void        GFitsHeaderCard::set_comment(const std::string& comment) { m_comment = comment; }


#endif /* GFITSHEADERCARD_HPP */
