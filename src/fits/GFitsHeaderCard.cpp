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
#include "GException.hpp"
#include "GFitsHeaderCard.hpp"

/* __ Namespaces _________________________________________________________ */

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

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***************************************************************************
 *                              Copy constructor                           *
 * ----------------------------------------------------------------------- *
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


/***************************************************************************
 *                               Destructor                                *
 * ----------------------------------------------------------------------- *
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

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
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

/***************************************************************************
 *                          Return value as string                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsHeaderCard::string(void)
{
    // Initialize return value to empty string
    std::string result;

    // Type dependent conversion
    switch (m_value_type) {
    case CT_STRING:
        if (m_value.length() > 2)
            result = m_value.substr(1, m_value.length() - 2);
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


/***************************************************************************
 *                          Return value as double                         *
 * ----------------------------------------------------------------------- *
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


/***************************************************************************
 *                            Return value as int                          *
 * ----------------------------------------------------------------------- *
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


/***************************************************************************
 *                             Read header card                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::read(__fitsfile* fptr, int keynum)
{
    // Read keyword
    char keyname[80];
    char value[80];
    char comment[80];
    int  status = 0;
    status      = __ffgkyn(fptr, keynum, keyname, value, comment, &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // Store result
    m_keyname.assign(keyname);
    m_value.assign(value);
    m_comment.assign(comment);

    // Determine card type
    m_value_type = get_value_type();

    // Return
    return;
}


/***************************************************************************
 *                             Read header card                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::read(__fitsfile* fptr, std::string keyname)
{
    // Read keyword
    char value[80];
    char comment[80];
    int  status = 0;
    status      = __ffgkey(fptr, (char*)keyname.c_str(), value, comment, &status);

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
    m_value_type = get_value_type();

    // Return
    return;
}


/***************************************************************************
 *                            Write header card                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::write(__fitsfile* fptr)
{
    // Write keyword
    int status = 0;
    status     = __ffuky(fptr, __TSTRING, (char*)m_keyname.c_str(), 
                         (char*)m_value.c_str(), (char*)m_comment.c_str(), &status);
    if (status != 0)
        throw GException::fits_error(G_WRITE, status);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                      GFitsHeaderCard private methods                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::init_members(void)
{
    // Initialise members
    m_keyname.clear();
    m_value.clear();
    m_value_type = CT_UNKNOWN;
    m_unit.clear();
    m_comment.clear();

    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::copy_members(const GFitsHeaderCard& card)
{
    // Copy membres
    m_keyname    = card.m_keyname;
    m_value      = card.m_value;
    m_value_type = card.m_value_type;
    m_unit       = card.m_unit;
    m_comment    = card.m_comment;

    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::free_members(void)
{
    // Return
    return;
}


/***************************************************************************
 *                        Determine card value type                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
int GFitsHeaderCard::get_value_type(void)
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
        if ((m_value[i] >= '0' && m_value[i] <= '9') || m_value[i] == '+' || m_value[i] == '-')
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
