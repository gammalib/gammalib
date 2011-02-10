/***************************************************************************
 *               GFitsHeader.cpp  - FITS header handling class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsHeader.cpp
 * @brief FITS header class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsHeader.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CARD                                      "GFitsHeader::card(int&)"
#define G_CARD_PTR                      "GFitsHeader::card_ptr(std::string&)"
#define G_OPEN                                     "GFitsHeader::open(void*)"
#define G_SAVE                                     "GFitsHeader::save(void*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsHeader::GFitsHeader(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param header Header from which GFitsHeader instance should be built
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


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFitsHeader::~GFitsHeader(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param header Header which should be assigned to GFitsHeader instance
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
 *
 * Sets the FITS header in an initial proper state.
 ***************************************************************************/
void GFitsHeader::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns number of header cards
 ***************************************************************************/
int GFitsHeader::size(void) const
{
    // Return number of header cards
    return m_num_cards;
}


/***********************************************************************//**
 * @brief Open Header
 *
 * @param[in] vptr FITS file void pointer.
 *
 * @exception GException::fits_error
 *            FITS error occured.
 *
 * Loads all header cards into memory. Any header cards that existed before
 * will be dropped.
 ***************************************************************************/
void GFitsHeader::open(void* vptr)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(FPTR(vptr), (FPTR(vptr)->HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Determine number of cards in header
    status = __ffghsp(FPTR(vptr), &m_num_cards, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Drop any old cards
    if (m_card != NULL) delete [] m_card;

    // Allocate memory for new cards
    m_card = new GFitsHeaderCard[m_num_cards];

    // Read all cards
    for (int i = 0; i < m_num_cards; ++i)
        m_card[i].read(FPTR(vptr), i+1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Header to FITS file
 *
 * @param[in] vptr FITS file void pointer.
 *
 * @exception GException::fits_error
 *            FITS error occured.
 *
 * Saves all header cards into HDU. This method does not write the following
 * mandatory keywords to the HDU (those will be written by methods handling
 * the data of the HDU):
 * 'SIMPLE'
 * 'BITPIX'
 * 'NAXIS', 'NAXIS1', 'NAXIS2', etc.
 * 'EXTEND'
 * 'PCOUNT'
 * 'GCOUNT'
 ***************************************************************************/
void GFitsHeader::save(void* vptr)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(FPTR(vptr), (FPTR(vptr)->HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_SAVE, status);

    // Save all cards
    for (int i = 0; i < m_num_cards; ++i) {
        if (m_card[i].keyname() != "SIMPLE" &&
            m_card[i].keyname() != "XTENSION" &&
            m_card[i].keyname() != "BITPIX" &&
            m_card[i].keyname() != "EXTEND" &&
            m_card[i].keyname() != "PCOUNT" &&
            m_card[i].keyname() != "GCOUNT" &&
            m_card[i].keyname().find("NAXIS") == std::string::npos) {
            m_card[i].write(FPTR(vptr));
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close Header
 *
 * Drops all header cards. Note that closing does not save the cards to the
 * FITS file. Use the save method in case that cards should be save to the
 * FITS file.
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


/***********************************************************************//**
 * @brief Update card in header or append card to header
 *
 * @param card FITS header card that should be updated
 *
 * This method updates one header card. Updating means replacing any
 * existing card with the specified one or appending a new card to the
 * list of existing cards.
 ***************************************************************************/
void GFitsHeader::update(const GFitsHeaderCard& card)
{
    // If card keyname is not COMMENT or HISTORY, then check first if
    // card exists. If yes then update existing card
    if (card.keyname() != "COMMENT" && card.keyname() != "HISTORY") {
        try {
            GFitsHeaderCard* ptr = this->card(card.keyname());
            *ptr = card;
            return;
        }
        catch (GException::fits_key_not_found) {
            ;
        }
    }

    // Create memory to hold cards
    GFitsHeaderCard* tmp = new GFitsHeaderCard[m_num_cards+1];
    if (tmp != NULL) {

        // Copy over existing cards and remove old ones
        if (m_card != NULL) {
            for (int i = 0; i < m_num_cards; ++i)
                tmp[i] = m_card[i];
            delete [] m_card;
        }

        // Connect the new memory to the card pointer
        m_card = tmp;

        // Append new card to list
        m_card[m_num_cards] = card;

        // Increment number of cards
        m_num_cards++;

    } // endif: new memory was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pointer of header card
 *
 * @param[in] keyname Name of header card
 *
 * Returns a pointer on the header card.
 * If card was not found a 'fits_key_not_found' error will be thrown.
 ***************************************************************************/
GFitsHeaderCard* GFitsHeader::card(const std::string& keyname)
{
    return GFitsHeader::card_ptr(keyname);
}


/***********************************************************************//**
 * @brief Return pointer of header card
 *
 * @param[in] cardno Number of card in header
 *
 * Returns a pointer on the header card.
 * If card was not found a 'out_of_range' error will be thrown.
 ***************************************************************************/
GFitsHeaderCard* GFitsHeader::card(const int& cardno)
{
    // If card number is out of range then throw an exception
    if (cardno < 0 || cardno >= m_num_cards)
        throw GException::out_of_range(G_CARD, cardno, 0, m_num_cards-1);

    // Get card pointer
    GFitsHeaderCard* ptr = &(m_card[cardno]);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get specified header card value as string
 *
 * @param[in] keyname Name of header card.
 *
 * If the header card is not found an empty string is returned.
 ***************************************************************************/
std::string GFitsHeader::string(const std::string& keyname)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(keyname);

    // Get string value
    std::string value = (ptr != NULL) ? ptr->string() : "";

    // Return string
    return value;
}


/***********************************************************************//**
 * @brief Get specified header card value as string
 *
 * @param[in] cardno Header card number (starting from 0).
 *
 * If the header card number does not exist an empty string is returned.
 ***************************************************************************/
std::string GFitsHeader::string(const int& cardno)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(cardno);

    // Get string value
    std::string value = (ptr != NULL) ? ptr->string() : "";

    // Return string
    return value;
}


/***********************************************************************//**
 * @brief Get specified header card value as double precision value
 *
 * @param[in] keyname Name of header card.
 *
 * If the header card is not found 0 is returned.
 ***************************************************************************/
double GFitsHeader::real(const std::string& keyname)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(keyname);

    // Get string value
    double value = (ptr != NULL) ? ptr->real() : 0.0;

    // Return double
    return value;
}


/***********************************************************************//**
 * @brief Get specified header card value as double precision value
 *
 * @param[in] cardno Header card number (starting from 0).
 *
 * If the header card number does not exist 0 is returned.
 ***************************************************************************/
double GFitsHeader::real(const int& cardno)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(cardno);

    // Get string value
    double value = (ptr != NULL) ? ptr->real() : 0.0;

    // Return double
    return value;
}


/***********************************************************************//**
 * @brief Get specified header card value as integer value
 *
 * @param keyname Name of header card.
 *
 * If the header card is not found 0 is returned.
 ***************************************************************************/
int GFitsHeader::integer(const std::string& keyname)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(keyname);

    // Get string value
    int value = (ptr != NULL) ? ptr->integer() : 0;

    // Return double
    return value;
}


/***********************************************************************//**
 * @brief Get specified header card value as integer value
 *
 * @param[in] cardno Header card number (starting from 0).
 *
 * If the header card number does not exist 0 is returned.
 ***************************************************************************/
int GFitsHeader::integer(const int& cardno)
{
    // Get card pointer
    GFitsHeaderCard* ptr = GFitsHeader::card(cardno);

    // Get string value
    int value = (ptr != NULL) ? ptr->integer() : 0;

    // Return double
    return value;
}


/***********************************************************************//**
 * @brief Clone header object
 ***************************************************************************/
GFitsHeader* GFitsHeader::clone(void) const
{
    return new GFitsHeader(*this);
}


/***********************************************************************//**
 * @brief Returns number of header cards
 *
 * @todo Depreceated method. Use size() instead.
 ***************************************************************************/
int GFitsHeader::num_cards(void) const
{
    // Return number of header cards
    return m_num_cards;
}


/***********************************************************************//**
 * @brief Print FITS file information
 ***************************************************************************/
std::string GFitsHeader::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GFitsHeader ("+str(m_num_cards)+" cards) ===");

    // Append cards
    for (int i = 0; i < m_num_cards; ++i) {
        result.append("\n");
        result.append(m_card[i].print());
    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsHeader::init_members(void)
{
    // Initialise members
    m_num_cards = 0;
    m_card      = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] header Header which should be copied
 ***************************************************************************/
void GFitsHeader::copy_members(const GFitsHeader& header)
{
    // Copy attributes
    m_num_cards = header.m_num_cards;

    // Copy cards
    if (header.m_card != NULL && m_num_cards > 0) {
        m_card = new GFitsHeaderCard[m_num_cards];
        for (int i = 0; i < m_num_cards; ++i)
            m_card[i] = header.m_card[i];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
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


/***********************************************************************//**
 * @brief Get pointer on header card
 *
 * @param[in] keyname Name of the header card
 *
 * @exception GException::fits_key_not_found
 *            Requested header card was not found in header.
 *
 * Returns a pointer on the header card.
 * If card was not found a 'fits_key_not_found' error will be thrown.
 ***************************************************************************/
GFitsHeaderCard* GFitsHeader::card_ptr(const std::string& keyname)
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

    // If no card was found then throw an error
    if (ptr == NULL)
        throw GException::fits_key_not_found(G_CARD_PTR, keyname);

    // Return pointer
    return ptr;
}


/*==========================================================================
 =                                                                         =
 =                            GFitsHeader friends                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] header FITS header.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GFitsHeader& header)
{
     // Write header in output stream
    os << header.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] header FITS header.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GFitsHeader& header)
{
    // Write header in logger
    log << header.print();

    // Return logger
    return log;
}
