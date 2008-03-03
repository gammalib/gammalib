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
/**
 * @file GFitsHeaderCard.hpp
 * @brief GFitsHeaderCard class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSHEADERCARD_HPP
#define GFITSHEADERCARD_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsHeaderCard
 *
 * @brief Implements FITS header card interface
 *
 * This class implements FITS header card. Each card consists of a
 * keyname (string), a value (string, floating pointer, integer or logical)
 * and a comment (string). COMMENT or HISTORY cards do not have any value.
 ***************************************************************************/
class GFitsHeaderCard {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsHeaderCard& card);

public:
    // Constructors & Destructors
    GFitsHeaderCard();
    GFitsHeaderCard(const std::string& keyname, const std::string& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const double& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const int& value,
                    const std::string& comment);
    GFitsHeaderCard(const GFitsHeaderCard& card);
    ~GFitsHeaderCard();

    // Operators
    GFitsHeaderCard& operator= (const GFitsHeaderCard& card);

    // Methods to set card properties
    void         keyname(const std::string& keyname);
    void         value(const std::string& value);
    void         value(const double& value);
    void         value(const int& value);
    void         unit(const std::string& unit);
    void         comment(const std::string& comment);

    // Methods to get card properties
    std::string  keyname(void) const;
    std::string  value(void) const;
    int          value_type(void) const;
    std::string  unit(void) const;
    std::string  comment(void) const;
    std::string  string(void);
    double       real(void);
    int          integer(void);

    // Other methods
    void         read(__fitsfile* fptr, int keynum);
    void         read(__fitsfile* fptr, const std::string& keyname);
    void         write(__fitsfile* fptr);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsHeaderCard& card);
    void free_members(void);
    int  get_value_type(const std::string& value);

    // Private data area
    std::string m_keyname;      //!< Name of the card
    std::string m_value;        //!< Value of the card as read from file
    int         m_value_type;   //!< Type of the card value
    std::string m_unit;         //!< Unit of the card value
    std::string m_comment;      //!< Card comment
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GFITSHEADERCARD_HPP */
