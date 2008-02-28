/***************************************************************************
 *               GFitsHeader.cpp  - FITS header handling class             *
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
#include "GFitsHeader.hpp"
#include <iostream>                           // cout, cerr

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_OPEN "GFitsHeader::open(int)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                    GFitsHeader constructors/destructors                 =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHeader::GFitsHeader()
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
GFitsHeader::GFitsHeader(const GFitsHeader& header)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(header);

    // Return
    return;
}


/***************************************************************************
 *                               Destructor                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHeader::~GFitsHeader()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsHeader operators                          =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHeader& GFitsHeader::operator= (const GFitsHeader& header)
{
    // Execute only if object is not identical
    if (this != &header) {
  
        // Free members
        free_members();
  
        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(header);
	
    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                        GFitsHeader public methods                       =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Open Header                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeader::open(__fitsfile* fptr)
{
    // Determine number of cards in header
    int status = 0;
    status     = __ffghsp(fptr, &m_num_cards, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);
    
    // Drop any old cards
    if (m_card != NULL) delete [] m_card;

    // Allocate memory for new cards
    m_card = new GFitsHeaderCard[m_num_cards];
    
    // Read all cards
    for (int i = 0; i < m_num_cards; ++i) {
        m_card[i].read(fptr, i+1);
        /*
        cout << m_card[i].keyname() << " " 
             << m_card[i].value() << " " 
             << m_card[i].value_type() << " " 
             << m_card[i].string() << " " 
             << m_card[i].real() << " " 
             << m_card[i].integer() << " " 
             << m_card[i].comment() << endl;
        */
    }

    // Return
    return;
}


/***************************************************************************
 *                              Close Header                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeader::close(void)
{
    // Free members
    free_members();
  
    // Initialise members
    init_members();
    
    // Return
    return;
}


/***************************************************************************
 *                 Get specified header card value as string               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsHeader::string(const std::string keyname)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(keyname);
    
    // Get string value
    std::string value = (ptr != NULL) ? ptr->string() : "";
    
    // Return string
    return value;
}


/***************************************************************************
 *                 Get specified header card value as string               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsHeader::string(const int cardno)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(cardno);
    
    // Get string value
    std::string value = (ptr != NULL) ? ptr->string() : "";
    
    // Return string
    return value;
}


/***************************************************************************
 *                 Get specified header card value as double               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GFitsHeader::real(const std::string keyname)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(keyname);
    
    // Get string value
    double value = (ptr != NULL) ? ptr->real() : 0.0;
    
    // Return double
    return value;
}


/***************************************************************************
 *                 Get specified header card value as double               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GFitsHeader::real(const int cardno)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(cardno);
    
    // Get string value
    double value = (ptr != NULL) ? ptr->real() : 0.0;
    
    // Return double
    return value;
}


/***************************************************************************
 *                  Get specified header card value as int                 *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
int GFitsHeader::integer(const std::string keyname)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(keyname);
    
    // Get string value
    int value = (ptr != NULL) ? ptr->integer() : 0;
    
    // Return double
    return value;
}


/***************************************************************************
 *                  Get specified header card value as int                 *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
int GFitsHeader::integer(const int cardno)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(cardno);
    
    // Get string value
    int value = (ptr != NULL) ? ptr->integer() : 0;
    
    // Return double
    return value;
}


/*==========================================================================
 =                                                                         =
 =                         GFitsHeader private methods                     =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeader::init_members(void)
{
    // Initialise members
    m_num_cards = 0;
    m_card      = NULL;
  
    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeader::copy_members(const GFitsHeader& header)
{
    // Copy attributes
    
    // Copy cards
    if (header.m_card != NULL && header.m_num_cards > 0) {
        m_num_cards = header.m_num_cards;
        m_card      = new GFitsHeaderCard[m_num_cards];
        for (int i = 0; i < m_num_cards; ++i)
            m_card[i] = header.m_card[i];
    }
    
    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHeader::free_members(void)
{
    // Free memory
    if (m_card != NULL) delete [] m_card;
    
    // Properly mark as free
    m_card = NULL;

    // Return
    return;
}


/***************************************************************************
 *                             Get card pointer                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHeaderCard* GFitsHeader::card_ptr(const std::string keyname)
{

    // Set card pointer to NULL (default)
    GFitsHeaderCard* ptr = NULL;
    
    // Search keyname in list
    for (int i = 0; i < m_num_cards; ++i) {
        if (m_card[i].keyname() == keyname) {
            ptr = &(m_card[i]);
            break;
        }
    }
    
    // Return pointer
    return ptr;
}


/*==========================================================================
 =                                                                         =
 =                            GFitsHeader friends                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                    Other functions used by GFitsHeader                  =
 =                                                                         =
 ==========================================================================*/
