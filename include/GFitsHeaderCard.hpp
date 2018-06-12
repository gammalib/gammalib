/***************************************************************************
 *               GFitsHeaderCard.hpp - FITS header card class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2018 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

#ifndef GFITSHEADERCARD_HPP
#define GFITSHEADERCARD_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GFitsHeaderCard
 *
 * @brief Implements FITS header card interface
 *
 * This class implements a FITS header card. A header card consists of a
 * keyname (string), a value (string, floating point, integer or logical)
 * and a comment (string). COMMENT or HISTORY cards do not have a value.
 ***************************************************************************/
class GFitsHeaderCard : public GBase {

    // Friend classes
    friend class GFitsHeader;

public:
    // Constructors & Destructors
    GFitsHeaderCard(void);
    GFitsHeaderCard(const std::string& keyname, const std::string& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const double& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const int& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const bool& value,
                    const std::string& comment);
    GFitsHeaderCard(const std::string& keyname, const char* value,
                    const std::string& comment);
    GFitsHeaderCard(const GFitsHeaderCard& card);
    virtual ~GFitsHeaderCard(void);

    // Operators
    GFitsHeaderCard& operator=(const GFitsHeaderCard& card);

    // Methods
    void               clear(void);
    GFitsHeaderCard*   clone(void) const;
    std::string        classname(void) const;
    void               keyname(const std::string& keyname);
    const std::string& keyname(void) const;
    void               value(const std::string& value);
    void               value(const bool& value);
    void               value(const float& value);
    void               value(const double& value);
    void               value(const unsigned short& value);
    void               value(const short& value);
    void               value(const unsigned int& value);
    void               value(const int& value);
    void               value(const long& value);
    void               value(const unsigned long& value);
    void               value(const long long& value);
    const std::string& value(void) const;
    const int&         decimals(void) const;
    void               unit(const std::string& unit);
    const std::string& unit(void) const;
    void               comment(const std::string& comment);
    const std::string& comment(void) const;
    std::string        string(void) const;
    double             real(void) const;
    int                integer(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsHeaderCard& card);
    void free_members(void);
    void copy_dtype(const GFitsHeaderCard& card);
    void free_dtype(void);
    void set_dtype(const std::string& value);
    void read(void* vptr, const int& keynum);
    void read(void* fptr, const std::string& keyname);
    void write(void* fptr) const;

    // Private data area
    std::string m_keyname;         //!< Name of the card
    std::string m_value;           //!< Value of the card as read from file
    std::string m_unit;            //!< Unit of the card value
    std::string m_comment;         //!< Card comment
    int         m_dtype;           //!< Native data type
    int         m_value_decimals;  //!< Decimals of value (for float)
    bool        m_comment_write;   //!< Signals that comment should be written
    void*       m_value_dtype;     //!< Value in native data type
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFitsHeaderCard").
 ***************************************************************************/
inline
std::string GFitsHeaderCard::classname(void) const
{
    return ("GFitsHeaderCard");
}


/***********************************************************************//**
 * @brief Set unit of header card value
 *
 * @param[in] unit Unit of header card.
 ***************************************************************************/
inline
void GFitsHeaderCard::unit(const std::string& unit)
{
    m_unit = unit;
    return;
}


/***********************************************************************//**
 * @brief Set comment of header card
 *
 * @param[in] comment Header card comment.
 ***************************************************************************/
inline
void GFitsHeaderCard::comment(const std::string& comment)
{
    m_comment = comment;
    return;
}


/***********************************************************************//**
 * @brief Return header card keyname
  ***************************************************************************/
inline
const std::string& GFitsHeaderCard::keyname(void) const
{
    return m_keyname;
}


/***********************************************************************//**
 * @brief Return header card value
  ***************************************************************************/
inline
const std::string& GFitsHeaderCard::value(void) const
{
    return m_value;
}


/***********************************************************************//**
 * @brief Return header card decimals
  ***************************************************************************/
inline
const int& GFitsHeaderCard::decimals(void) const
{
    return m_value_decimals;
}


/***********************************************************************//**
 * @brief Return header card value unit
  ***************************************************************************/
inline
const std::string& GFitsHeaderCard::unit(void) const
{
    return m_unit;
}


/***********************************************************************//**
 * @brief Return header card comment
  ***************************************************************************/
inline
const std::string& GFitsHeaderCard::comment(void) const
{
    return m_comment;
}

#endif /* GFITSHEADERCARD_HPP */
