/***************************************************************************
 *       GFitsHeaderCard.cpp  - FITS header card abstract base class       *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsHeaderCard.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ_NUM   "GFitsHeaderCard::read(fitsfile*, int)"
#define G_READ_STR   "GFitsHeaderCard::read(fitsfile*, std::string)"
#define G_WRITE      "GFitsHeaderCard::write(fitsfile*)"

/* __ Definitions ________________________________________________________ */
#define CT_UNKNOWN 0
#define CT_INVALID 1
#define CT_STRING  2
#define CT_INT     3
#define CT_FLOAT   4
#define CT_BOOL    5
#define CT_COMMENT 6
#define CT_HISTORY 7

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Enumerations _______________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                     GFitsHDU constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor for string cards
 *
 * @param[in] keyname Name of the header card
 * @param[in] value String value of the header card
 * @param[in] comment Comment of the header card
 *
 * This constructor builds a header card from the keyname, value and comment.
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(const std::string& keyname,
                                 const std::string& value,
                                 const std::string& comment)
{
    // Initialise class members for clean destruction
    init_members();

    // Set members
    this->keyname(keyname);
    this->value("'"+value+"'");
    this->comment(comment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor for floating point cards
 *
 * @param[in] keyname Name of the header card
 * @param[in] value Floating point value of the header card
 * @param[in] comment Comment of the header card
 *
 * This constructor builds a header card from the keyname, value and comment.
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(const std::string& keyname,
                                 const double& value,
                                 const std::string& comment)
{
    // Initialise class members for clean destruction
    init_members();

    // Set members
    this->keyname(keyname);
    this->value(value);
    this->comment(comment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor for integer cards
 *
 * @param[in] keyname Name of the header card
 * @param[in] value Integer value of the header card
 * @param[in] comment Comment of the header card
 *
 * This constructor builds a header card from the keyname, value and comment.
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(const std::string& keyname,
                                 const int& value,
                                 const std::string& comment)
{
    // Initialise class members for clean destruction
    init_members();

    // Set members
    this->keyname(keyname);
    this->value(value);
    this->comment(comment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param card[in] Header card that will be used to construct GFitsHeaderCard
 *                 instance
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(const GFitsHeaderCard& card)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(card);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFitsHeaderCard::~GFitsHeaderCard()
{
    // Free members
    free_members();

    // Return
    return;
}

/*==========================================================================
 =                                                                         =
 =                         GFitsHeaderCard operators                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param card[in] Header card that will be assigned to GFitsHeaderCard 
 *                 instance
 ***************************************************************************/
GFitsHeaderCard& GFitsHeaderCard::operator= (const GFitsHeaderCard& card)
{
    // Execute only if object is not identical
    if (this != &card) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(card);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                      GFitsHeaderCard public methods                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set name of header card
 *
 * @param[in] keyname Name of header card
 ***************************************************************************/
void GFitsHeaderCard::keyname(const std::string& keyname)
{
    // Set name of card
    m_keyname = keyname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set string value of header card
 *
 * @param[in] value String value of header card
 ***************************************************************************/
void GFitsHeaderCard::value(const std::string& value)
{
    // Set value and value type
    m_value      = value;
    m_value_type = get_value_type(value);

    // Attach hyphens if required
    if (m_value_type == CT_STRING) {
        if (m_value[0] != '\x27')
            m_value = "'" + m_value;
        if (m_value[m_value.length()-1] != '\x27')
            m_value += "'";
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set floating point value of header card
 *
 * @param[in] value Floating point value of header card
 ***************************************************************************/
void GFitsHeaderCard::value(const double& value)
{
    // Convert value to string
    std::ostringstream s_value;
    s_value << std::scientific << value;

    // Assign value
    m_value = s_value.str();

    // Set value type
    m_value_type = get_value_type(m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set integer value of header card
 *
 * @param[in] value Integer value of header card
 ***************************************************************************/
void GFitsHeaderCard::value(const int& value)
{
    // Convert value to string
    std::ostringstream s_value;
    s_value << value;

    // Assign value
    m_value = s_value.str();

    // Set value type
    m_value_type = get_value_type(m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set unit of header card value
 *
 * @param[in] value Unit of header card
 ***************************************************************************/
void GFitsHeaderCard::unit(const std::string& unit)
{
    // Set unit
    m_unit = unit;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set comment of header card
 *
 * @param[in] value Header card comment
 ***************************************************************************/
void GFitsHeaderCard::comment(const std::string& comment)
{
    // Set comment
    m_comment = comment;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return header card keyname
  ***************************************************************************/
std::string GFitsHeaderCard::keyname(void) const
{
    // Return
    return m_keyname;
}


/***********************************************************************//**
 * @brief Return header card value
  ***************************************************************************/
std::string GFitsHeaderCard::value(void) const
{
    // Return
    return m_value;
}


/***********************************************************************//**
 * @brief Return header card value type
  ***************************************************************************/
int GFitsHeaderCard::value_type(void) const
{
    // Return
    return m_value_type;
}


/***********************************************************************//**
 * @brief Return header card decimals
  ***************************************************************************/
int GFitsHeaderCard::decimals(void) const
{
    // Return
    return m_value_decimals;
}


/***********************************************************************//**
 * @brief Return header card value unit
  ***************************************************************************/
std::string GFitsHeaderCard::unit(void) const
{
    // Return
    return m_unit;
}


/***********************************************************************//**
 * @brief Return header card comment
  ***************************************************************************/
std::string GFitsHeaderCard::comment(void) const
{
    // Return
    return m_comment;
}


/***********************************************************************//**
 * @brief Return header card value as string
 *
 * Convert header card value into a string.
 * Any hyphens that may occur in the FITS card will be automatically stripped.
 ***************************************************************************/
std::string GFitsHeaderCard::string(void)
{
    // Initialize return value to empty string
    std::string result;

    // Type dependent conversion
    switch (m_value_type) {
    case CT_STRING:
        if (m_value.length() > 2)
            result = strip_whitespace(m_value.substr(1, m_value.length() - 2));
        break;
    case CT_INT:
    case CT_FLOAT:
    case CT_BOOL:
        result = m_value;
        break;
    default:
        break;
    }

    // Return string
    return result;

}


/***********************************************************************//**
 * @brief Return header card value as double precision
 *
 * Convert header card value into a double precision value. In case that the
 * card did not contain a numerical value, 0 will be returned by the
 * method.
 ***************************************************************************/
double GFitsHeaderCard::real(void)
{
    // Initialize return value to 0.0
    double result = 0.0;

    // Type dependent conversion
    switch (m_value_type) {
    case CT_STRING:
        if (m_value.length() > 2)
            result = atof(m_value.substr(1, m_value.length() - 2).c_str());
        break;
    case CT_INT:
    case CT_FLOAT:
        result = atof(m_value.c_str());
        break;
    case CT_BOOL:
        result = (m_value == "T") ? 1.0 : 0.0;
        break;
    default:
        break;
    }

    // Return string
    return result;

}


/***********************************************************************//**
 * @brief Return header card value as integer
*
 * Convert header card value into a integer value. In case that the
 * card did not contain a numerical value, 0 will be returned by the
 * method.
 ***************************************************************************/
int GFitsHeaderCard::integer(void)
{
    // Initialize return value to 0
    int result = 0;

    // Type dependent conversion
    switch (m_value_type) {
    case CT_STRING:
        if (m_value.length() > 2)
            result = atoi(m_value.substr(1, m_value.length() - 2).c_str());
        break;
    case CT_INT:
    case CT_FLOAT:
        result = atoi(m_value.c_str());
        break;
    case CT_BOOL:
        result = (m_value == "T") ? 1 : 0;
        break;
    default:
        break;
    }

    // Return string
    return result;

}


/*==========================================================================
 =                                                                         =
 =                      GFitsHeaderCard private methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsHeaderCard::init_members(void)
{
    // Initialise members
    m_keyname.clear();
    m_value.clear();
    m_value_type     = CT_UNKNOWN;
    m_value_decimals = 10;
    m_unit.clear();
    m_comment.clear();
    m_comment_write = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] card Header card to be copied
 ***************************************************************************/
void GFitsHeaderCard::copy_members(const GFitsHeaderCard& card)
{
    // Copy membres
    m_keyname        = card.m_keyname;
    m_value          = card.m_value;
    m_value_type     = card.m_value_type;
    m_value_decimals = card.m_value_decimals;
    m_unit           = card.m_unit;
    m_comment        = card.m_comment;
    m_comment_write  = card.m_comment_write;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsHeaderCard::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Determine card value type
 *
 * @param[in] value Value string for which the type should be determined
 *
 * Method that determines the type of the card value. The following types are
 * defined:
 * CT_UNKNOWN=0 (inknown type)
 * CT_INVALID=1 (invalid type)
 * CT_STRING=2 (string)
 * CT_INT=3 (integer)
 * CT_FLOAT=4 (floating point)
 * CT_BOOL=5 (boolean)
 * CT_COMMENT=6 (comment)
 * CT_HISTORY=7 (history)
 * The type is determined by analyzing the card value. It is stored in the
 * private variable 'm_value_type'.
 ***************************************************************************/
int GFitsHeaderCard::get_value_type(const std::string& value)
{
    // Get index of last string element
    int last = m_value.length() - 1;

    // If value is empty then we either have a COMMENT or HISTORY card or
    // we don't know the type
    if (last < 0) {
        if (m_keyname == "COMMENT")
            return CT_COMMENT;
        if (m_keyname == "HISTORY")
            return CT_HISTORY;
        return CT_UNKNOWN;
    }

    // If value is enclosed in parentheses we have a string
    if (m_value[0] == '\x27' && m_value[last] == '\x27')
        return CT_STRING;

    // If value has only one digit and is either 'F' or 'T' we have a
    // Boolean
    if (last == 0 && (m_value[0] == 'F' || m_value[0] == 'T'))
        return CT_BOOL;

    // If value is composed only of numeric or '+' or '-' characters we
    // have an integer. If in addition we have '.' or 'e' or 'E' we (may)
    // have a float.
    int found_int   = 0;
    int found_float = 0;
    for (int i = 0; i <= last; ++i) {
        if ((m_value[i] >= '0' && m_value[i] <= '9') || m_value[i] == '+' ||
             m_value[i] == '-')
            found_int = 1;
        else if (m_value[i] == '.' || m_value[i] == 'e' || m_value[i] == 'E')
            found_float = 1;
        else
            return CT_INVALID;
    }

    // If we are still alive we should have now either a float or an integer
    if (found_int) {
        if (found_float) {
            return CT_FLOAT;
        }
        else {
            return CT_INT;
        }
    }
    else
        return CT_INVALID;
}


/***********************************************************************//**
 * @brief Read header card from FITS file
 *
 * @param[in] fptr Pointer to FITS file from which the card is read
 * @param[in] keynum Number of the header card
 ***************************************************************************/
void GFitsHeaderCard::read(__fitsfile* fptr, int keynum)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(fptr, (fptr->HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // Read keyword
    char keyname[80];
    char value[80];
    char comment[80];
    status = __ffgkyn(fptr, keynum, keyname, value, comment, &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // Store result
    m_keyname.assign(keyname);
    m_value.assign(value);
    m_comment.assign(comment);

    // Determine card type
    m_value_type = get_value_type(m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read header card from FITS file
 *
 * @param[in] fptr Pointer to FITS file from which the card is read
 * @param[in] keyname Name of the header card
 ***************************************************************************/
void GFitsHeaderCard::read(__fitsfile* fptr, const std::string& keyname)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(fptr, (fptr->HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // Read keyword
    char value[80];
    char comment[80];
    status = __ffgkey(fptr, (char*)keyname.c_str(), value, comment, &status);

    // Catch error
    if (status == 202)       // Keyword not found
        throw GException::fits_key_not_found(G_READ_STR, keyname, status);
    else if (status != 0)    // Any other error
        throw GException::fits_error(G_READ_STR, status);

    // Store result
    m_keyname = keyname;
    m_value.assign(value);
    m_comment.assign(comment);

    // Determine card type
    m_value_type = get_value_type(m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write header card
 *
 * @param fptr[in] Pointer to FITS file into which the header card should be
 *                  written
 *
 * Writes any kind of header card to a FITS file.
 ***************************************************************************/
void GFitsHeaderCard::write(__fitsfile* fptr)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(fptr, (fptr->HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // Write keyword
    switch (m_value_type) {
    case CT_UNKNOWN:
        break;
    case CT_INVALID:
        break;
    case CT_STRING:
        status = __ffukys(fptr, (char*)m_keyname.c_str(),
                          (char*)string().c_str(), (char*)m_comment.c_str(), 
                          &status);
        break;
    case CT_INT:
        status = __ffukyj(fptr, (char*)m_keyname.c_str(), integer(),
                          (char*)m_comment.c_str(), &status);
        break;
    case CT_FLOAT:
        status = __ffukyd(fptr, (char*)m_keyname.c_str(), real(),
                          decimals(), (char*)m_comment.c_str(), &status);
        break;
    case CT_BOOL:
        status = __ffukyl(fptr, (char*)m_keyname.c_str(), integer(),
                          (char*)m_comment.c_str(), &status);
        break;
    case CT_COMMENT:
        if (m_comment_write)
            status = __ffpcom(fptr, (char*)m_comment.c_str(), &status);
        break;
    case CT_HISTORY:
        if (m_comment_write)
            status = __ffphis(fptr, (char*)m_comment.c_str(), &status);
        break;
    default:
        break;
    }
    if (status != 0)
        throw GException::fits_error(G_WRITE, status);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GFitsHeaderCard friends                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] card Card that should be put in output stream
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GFitsHeaderCard& card)
{
    // Set buffer
    char buffer[100];

    // Setup left justified keyname
    std::string keyname = card.m_keyname;
    while (keyname.length() < 8)
        keyname.append(" ");

    // Setup buffer
    if (card.m_value_type != CT_COMMENT && card.m_value_type != CT_HISTORY) {
      if (card.m_unit.length() > 0) {
        sprintf(buffer, "%8s =%21s / [%s] %s",
                keyname.c_str(),
                card.m_value.c_str(),
                card.m_unit.c_str(),
                card.m_comment.c_str());
      }
      else {
        sprintf(buffer, "%8s =%21s / %s",
                keyname.c_str(),
                card.m_value.c_str(),
                card.m_comment.c_str());
      }
    }
    else {
      sprintf(buffer, "%8s %s",
              keyname.c_str(),
              card.m_comment.c_str());
    }

    // Put buffer in stream
    os << buffer;

    // Attach card type
    switch (card.m_value_type) {
    case CT_UNKNOWN:
        os << " <type unknown>" << std::endl;
        break;
    case CT_INVALID:
        os << " <invalid type>" << std::endl;
        break;
    case CT_STRING:
        os << " <string>" << std::endl;
        break;
    case CT_INT:
        os << " <int>" << std::endl;
        break;
    case CT_FLOAT:
        os << " <float>" << std::endl;
        break;
    case CT_BOOL:
        os << " <boolean>" << std::endl;
        break;
    case CT_COMMENT:
        os << " <comment>" << std::endl;
        break;
    case CT_HISTORY:
        os << " <history>" << std::endl;
        break;
    default:
        os << std::endl;
        break;
    }

    // Return output stream
    return os;
}
