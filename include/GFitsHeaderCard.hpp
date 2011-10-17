/***************************************************************************
 *               GFitsHeaderCard.hpp  - FITS header card class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsHeaderCard.hpp
 * @brief FITS header card class definition
 * @author J. Knodlseder
 */

#ifndef GFITSHEADERCARD_HPP
#define GFITSHEADERCARD_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"


/***********************************************************************//**
 * @class GFitsHeaderCard
 *
 * @brief Implements FITS header card interface
 *
 * This class implements FITS header card. Each card consists of a
 * keyname (string), a value (string, floating pointer, integer or logical)
 * and a comment (string). COMMENT or HISTORY cards do not have any value.
 *
 * @todo Many more datatypes may exist for a header card.
 ***************************************************************************/
class GFitsHeaderCard {

    // Friend classes
    friend class GFitsHeader;

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GFitsHeaderCard& card);
    friend GLog&         operator<<(GLog& log,        const GFitsHeaderCard& card);

public:
    // Constructors & Destructors
    GFitsHeaderCard(void);
    GFitsHeaderCard(const std::string& keyname, const std::string& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const double& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const int& value,
                    const std::string& comment);
    GFitsHeaderCard(const GFitsHeaderCard& card);
    virtual ~GFitsHeaderCard(void);

    // Operators
    GFitsHeaderCard& operator= (const GFitsHeaderCard& card);

    // Methods to set card properties
    void         keyname(const std::string& keyname);
    void         value(const std::string& value);
    void         value(const bool& value);
    void         value(const float& value);
    void         value(const double& value);
    void         value(const unsigned short& value);
    void         value(const short& value);
    void         value(const unsigned int& value);
    void         value(const int& value);
    void         value(const long& value);
    void         value(const unsigned long& value);
    void         value(const long long& value);
    void         unit(const std::string& unit);
    void         comment(const std::string& comment);

    // Methods to get card properties
    std::string  keyname(void) const;
    std::string  value(void) const;
    int          decimals(void) const;
    std::string  unit(void) const;
    std::string  comment(void) const;
    std::string  string(void);
    double       real(void);
    int          integer(void);
    std::string  print(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsHeaderCard& card);
    void free_members(void);
    void copy_dtype(const GFitsHeaderCard& card);
    void free_dtype(void);
    void set_dtype(const std::string& value);
    void read(void* vptr, int keynum);
    void read(void* fptr, const std::string& keyname);
    void write(void* fptr);

    // Private data area
    std::string m_keyname;         //!< Name of the card
    std::string m_value;           //!< Value of the card as read from file
    void*       m_value_dtype;     //!< Value in native data type
    int         m_dtype;           //!< Native data type
    int         m_value_decimals;  //!< Decimals of value (for float)
    std::string m_unit;            //!< Unit of the card value
    std::string m_comment;         //!< Card comment
    bool        m_comment_write;   //!< Signals that comment should be written
};

#endif /* GFITSHEADERCARD_HPP */
