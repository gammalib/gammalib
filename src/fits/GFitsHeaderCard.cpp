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

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

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
 *                               Get keyname                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsHeaderCard::keyname(void)
{
    return m_keyname;
}


/***************************************************************************
 *                                Get value                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsHeaderCard::value(void)
{
    return m_value;
}


/***************************************************************************
 *                                 Get unit                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsHeaderCard::unit(void)
{
    return m_unit;
}


/***************************************************************************
 *                               Get comment                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsHeaderCard::comment(void)
{
    return m_comment;
}


/***************************************************************************
 *                               Set keyname                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::set_keyname(std::string& keyname)
{
    m_keyname = keyname;
    return;
}


/***************************************************************************
 *                                Set value                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::set_value(std::string& value)
{
    m_value = value;
    return;
}


/***************************************************************************
 *                                 Set unit                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::set_unit(std::string& unit)
{
    m_unit = unit;
    return;
}


/***************************************************************************
 *                               Set comment                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeaderCard::set_comment(std::string& comment)
{
    m_comment = comment;
    return;
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
    if (status != 0) {
        throw GException::fits_error("GFitsHeaderCard::read(fitsfile*, int)", status);
    }
    
    // Store result
    m_keyname.assign(keyname);
    m_value.assign(value);
    m_comment.assign(comment);
    
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
    if (status == 202) {       // Keyword not found
        throw GException::fits_key_not_found("GFitsHeaderCard::read(fitsfile*, std::string)", keyname, status);
    }
    else if (status != 0) {    // Any other error
        throw GException::fits_error("GFitsHeaderCard::read(fitsfile*, std::string)", status);
    }
    
    // Store result
    m_keyname = keyname;
    m_value.assign(value);
    m_comment.assign(comment);
    
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
    if (status != 0) {
        throw GException::fits_error("GFitsHeaderCard::write(fitsfile*", status);
    }
    
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
    m_keyname = card.m_keyname;
    m_value   = card.m_value;
    m_unit    = card.m_unit;
    m_comment = card.m_comment;
    
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
