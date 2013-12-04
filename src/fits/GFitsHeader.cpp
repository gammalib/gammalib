/***************************************************************************
 *            GFitsHeader.hpp - FITS header cards container class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GFitsHeader.cpp
 * @brief FITS header cards container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsHeader.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT1                        "GFitsHeaderCard& GFitsHeader::at(int&)"
#define G_AT2                "GFitsHeaderCard& GFitsHeader::at(std::string&)"
#define G_STRING1                                 "GFitsHeader::string(int&)"
#define G_STRING2                         "GFitsHeader::string(std::string&)"
#define G_REAL1                                     "GFitsHeader::real(int&)"
#define G_REAL2                             "GFitsHeader::real(std::string&)"
#define G_INTEGER1                               "GFitsHeader::integer(int&)"
#define G_INTEGER2                       "GFitsHeader::integer(std::string&)"
#define G_INSERT1               "GFitsHeader::insert(int&, GFitsHeaderCard&)"
#define G_INSERT2       "GFitsHeader::insert(std::string&, GFitsHeaderCard&)"
#define G_REMOVE1                                 "GFitsHeader::remove(int&)"
#define G_REMOVE2                         "GFitsHeader::remove(std::string&)"
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
 * @param header FITS header.
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
 * @param header FITS header.
 * @return FITS header.
 ***************************************************************************/
GFitsHeader& GFitsHeader::operator=(const GFitsHeader& header)
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
 * @brief Clear header
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
 * @brief Clone header
 *
 * @return Pointer to deep copy of header.
 ***************************************************************************/
GFitsHeader* GFitsHeader::clone(void) const
{
    return new GFitsHeader(*this);
}


/***********************************************************************//**
 * @brief Return header card
 *
 * @param[in] cardno Number of card in header [0,...,size()-1]
 * @return Header card.
 *
 * @exception GException::out_of_range
 *            Card number out of range.
 ***************************************************************************/
GFitsHeaderCard& GFitsHeader::at(const int& cardno)
{
    // Compile option: raise an exception if cardno is out of range
    #if defined(G_RANGE_CHECK)
    if (!contains(cardno)) {
        throw GException::out_of_range(G_AT1, "Header card number", cardno, size());
    }
    #endif

    // Return card
    return (m_cards[cardno]);
}


/***********************************************************************//**
 * @brief Return header card (const version)
 *
 * @param[in] cardno Number of card in header [0,...,size()-1]
 * @return Header card.
 *
 * @exception GException::out_of_range
 *            Card number out of range.
 ***************************************************************************/
const GFitsHeaderCard& GFitsHeader::at(const int& cardno) const
{
    // Compile option: raise an exception if cardno is out of range
    #if defined(G_RANGE_CHECK)
    if (!contains(cardno)) {
        throw GException::out_of_range(G_AT1, "Header card number", cardno, size());
    }
    #endif

    // Return card
    return (m_cards[cardno]);
}


/***********************************************************************//**
 * @brief Return header card
 *
 * @param[in] keyname Name of header card
 * @return Header card.
 *
 * @exception GException::invalid_argument
 *            Key name not found.
 ***************************************************************************/
GFitsHeaderCard& GFitsHeader::at(const std::string& keyname)
{
    // Get card number
    int cardno = get_index(keyname);

    // Throw an exception if keyname is not found
    if (cardno == -1) {
        std::string msg = "Keyword \""+keyname+"\" not found in FITS header.";
        throw GException::invalid_argument(G_AT2, msg);
    }

    // Return card
    return (m_cards[cardno]);
}


/***********************************************************************//**
 * @brief Return header card (const version)
 *
 * @param[in] keyname Name of header card
 * @return Header card.
 *
 * @exception GException::invalid_argument
 *            Key name not found.
 ***************************************************************************/
const GFitsHeaderCard& GFitsHeader::at(const std::string& keyname) const
{
    // Get card number
    int cardno = get_index(keyname);

    // Throw an exception if keyname is not found
    if (cardno == -1) {
        std::string msg = "Keyword \""+keyname+"\" not found in FITS header.";
        throw GException::invalid_argument(G_AT2, msg);
    }

    // Return card
    return (m_cards[cardno]);
}


/***********************************************************************//**
 * @brief Return header card value as string value
 *
 * @param[in] cardno Header card number [0,...,size()-1].
 * @return Header card string value.
 *
 * @exception GException::out_of_range
 *            Card number out of range.
 ***************************************************************************/
std::string GFitsHeader::string(const int& cardno) const
{
    // Compile option: raise an exception if cardno is out of range
    #if defined(G_RANGE_CHECK)
    if (!contains(cardno)) {
        throw GException::out_of_range(G_STRING1, "Header card number", cardno, size());
    }
    #endif

    // Return card value
    return (m_cards[cardno].string());
}


/***********************************************************************//**
 * @brief Return header card value as string value
 *
 * @param[in] keyname Header card key name.
 * @return Header card string value.
 *
 * @exception GException::invalid_argument
 *            Key name not found.
 ***************************************************************************/
std::string GFitsHeader::string(const std::string& keyname) const
{
    // Get card number
    int cardno = get_index(keyname);

    // Throw an exception if keyname is not found
    if (cardno == -1) {
        std::string msg = "Keyword \""+keyname+"\" not found in FITS header.";
        throw GException::invalid_argument(G_STRING2, msg);
    }

    // Return card value
    return (m_cards[cardno].string());
}


/***********************************************************************//**
 * @brief Return header card value as double precision value
 *
 * @param[in] cardno Header card number [0,...,size()-1].
 * @return Header card double precision value.
 *
 * @exception GException::out_of_range
 *            Card number out of range.
 ***************************************************************************/
double GFitsHeader::real(const int& cardno) const
{
    // Compile option: raise an exception if cardno is out of range
    #if defined(G_RANGE_CHECK)
    if (!contains(cardno)) {
        throw GException::out_of_range(G_REAL1, "Header card number", cardno, size());
    }
    #endif

    // Return card value
    return (m_cards[cardno].real());
}


/***********************************************************************//**
 * @brief Return header card value as double precision value
 *
 * @param[in] keyname Header card key name.
 * @return Header card double precision value.
 *
 * @exception GException::invalid_argument
 *            Key name not found.
 ***************************************************************************/
double GFitsHeader::real(const std::string& keyname) const
{
    // Get card number
    int cardno = get_index(keyname);

    // Throw an exception if keyname is not found
    if (cardno == -1) {
        std::string msg = "Keyword \""+keyname+"\" not found in FITS header.";
        throw GException::invalid_argument(G_REAL2, msg);
    }

    // Return card value
    return (m_cards[cardno].real());
}


/***********************************************************************//**
 * @brief Return header card value as integer value
 *
 * @param[in] cardno Header card number [0,...,size()-1].
 * @return Header card integer value.
 *
 * @exception GException::out_of_range
 *            Card number out of range.
 ***************************************************************************/
int GFitsHeader::integer(const int& cardno) const
{
    // Compile option: raise an exception if cardno is out of range
    #if defined(G_RANGE_CHECK)
    if (!contains(cardno)) {
        throw GException::out_of_range(G_INTEGER1, "Header card number", cardno, size());
    }
    #endif

    // Return card value
    return (m_cards[cardno].integer());
}


/***********************************************************************//**
 * @brief Return header card value as integer value
 *
 * @param[in] keyname Header card key name.
 * @return Header card integer value.
 *
 * @exception GException::invalid_argument
 *            Key name not found.
 ***************************************************************************/
int GFitsHeader::integer(const std::string& keyname) const
{
    // Get card number
    int cardno = get_index(keyname);

    // Throw an exception if keyname is not found
    if (cardno == -1) {
        std::string msg = "Keyword \""+keyname+"\" not found in FITS header.";
        throw GException::invalid_argument(G_INTEGER2, msg);
    }

    // Return card value
    return (m_cards[cardno].integer());
}


/***********************************************************************//**
 * @brief Append or update header card
 *
 * @param[in] card Header card.
 * @return Reference to appended header card.
 *
 * If the keyname of the header @p card does not yet exist in the header
 * (or if the keyname is COMMENT or HISTORY) then append the header card
 * to the header. If the keyname exists already, the header card is updated.
 ***************************************************************************/
GFitsHeaderCard& GFitsHeader::append(const GFitsHeaderCard& card)
{
    // If card keyname is not COMMENT or HISTORY, then check first if
    // card exists. If yes then update existing card
    int cardno = -1;
    if (card.keyname() != "COMMENT" && card.keyname() != "HISTORY") {
        cardno = get_index(card.keyname());
        if (cardno != -1) {
            m_cards[cardno] = card;
        }
    }

    // If card has not yet been updated then append card to header
    if (cardno == -1) {
        m_cards.push_back(card);
        cardno = size()-1;
    }

    // Return reference
    return (m_cards[cardno]);
}


/***********************************************************************//**
 * @brief Insert card into header
 *
 * @param[in] cardno Header card number [0,...,size()-1].
 * @param[in] card Header card.
 * @return Reference to inserted header card.
 *
 * @exception GException::out_of_range
 *            Card number out of range.
 *
 * Inserts a @p card into the header before the card with the specified card
 * number @p cardno.
 ***************************************************************************/
GFitsHeaderCard& GFitsHeader::insert(const int&             cardno,
                                     const GFitsHeaderCard& card)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (isempty()) {
        if (cardno > 0) {
            throw GException::out_of_range(G_INSERT1, "Header card number",
                                           cardno, size());
        }
    }
    else {
        if (cardno < 0 || cardno >= size()) {
            throw GException::out_of_range(G_INSERT1, "Header card number",
                                           cardno, size());
        }
    }
    #endif

    // Inserts card
    m_cards.insert(m_cards.begin()+cardno, card);
    
    // Return reference
    return m_cards[cardno];
}


/***********************************************************************//**
 * @brief Insert card into header
 *
 * @param[in] keyname Header card key name.
 * @param[in] card Header card.
 * @return Reference to inserted header card.
 *
 * @exception GException::invalid_argument
 *            Key name not found.
 *
 * Inserts a @p card into the header before the card with the specified
 * @p keyname.
 ***************************************************************************/
GFitsHeaderCard& GFitsHeader::insert(const std::string&     keyname,
                                     const GFitsHeaderCard& card)
{
    // Get card number
    int cardno = get_index(keyname);

    // Throw an exception if keyname is not found
    if (cardno == -1) {
        std::string msg = "Keyword \""+keyname+"\" not found in FITS header.";
        throw GException::invalid_argument(G_INSERT2, msg);
    }

    // Inserts card
    m_cards.insert(m_cards.begin()+cardno, card);
    
    // Return reference
    return m_cards[cardno];
}


/***********************************************************************//**
 * @brief Remove card from header
 *
 * @param[in] cardno Header card number [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Card number out of range.
 *
 * Remove header card of specified card number from container.
 ***************************************************************************/
void GFitsHeader::remove(const int& cardno)
{
    // Compile option: raise an exception if cardno is out of range
    #if defined(G_RANGE_CHECK)
    if (!contains(cardno)) {
        throw GException::out_of_range(G_REMOVE1, "Header card number", cardno, size());
    }
    #endif

    // Erase card from header
    m_cards.erase(m_cards.begin() + cardno);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove card from header
 *
 * @param[in] keyname Header card key name.
 *
 * @exception GException::invalid_argument
 *            Key name not found.
 *
 * Remove card with specified @p keyname from header.
 ***************************************************************************/
void GFitsHeader::remove(const std::string& keyname)
{
    // Get card number
    int cardno = get_index(keyname);

    // Throw an exception if keyname is not found
    if (cardno == -1) {
        std::string msg = "Keyword \""+keyname+"\" not found in FITS header.";
        throw GException::invalid_argument(G_REMOVE2, msg);
    }

    // Erase card from header
    m_cards.erase(m_cards.begin() + cardno);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append header
 *
 * @param[in] header FITS header.
 *
 * Append all cards from a FITS header to the actual header.
 ***************************************************************************/
void GFitsHeader::extend(const GFitsHeader& header)
{
    // Do nothing if header is empty
    if (!header.isempty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = header.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all card and append them to the header
        for (int i = 0; i < num; ++i) {
            m_cards.push_back(header[i]);
        }

    } // endif: header was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load header from FITS file
 *
 * @param[in] vptr FITS file void pointer.
 *
 * @exception GException::fits_error
 *            FITS error occured.
 *
 * Loads all header cards into memory. Any header cards that existed before
 * will be dropped.
 ***************************************************************************/
void GFitsHeader::load(void* vptr)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(FPTR(vptr), (FPTR(vptr)->HDUposition)+1, NULL, &status);
    if (status != 0) {
        throw GException::fits_error(G_OPEN, status);
    }

    // Determine number of cards in header
    int num_cards = 0;
    status = __ffghsp(FPTR(vptr), &num_cards, NULL, &status);
    if (status != 0) {
        throw GException::fits_error(G_OPEN, status);
    }

    // Drop any old cards and reserve space for new cards
    m_cards.clear();
    reserve(num_cards);

    // Read all cards
    for (int i = 0; i < num_cards; ++i) {
        GFitsHeaderCard card;
        card.read(FPTR(vptr), i+1);
        m_cards.push_back(card);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save header to FITS file
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
void GFitsHeader::save(void* vptr) const
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(FPTR(vptr), (FPTR(vptr)->HDUposition)+1, NULL, &status);
    if (status != 0) {
        throw GException::fits_error(G_SAVE, status);
    }

    // Save all cards
    for (int i = 0; i < size(); ++i) {
        if (m_cards[i].keyname() != "SIMPLE" &&
            m_cards[i].keyname() != "XTENSION" &&
            m_cards[i].keyname() != "BITPIX" &&
            m_cards[i].keyname() != "EXTEND" &&
            m_cards[i].keyname() != "PCOUNT" &&
            m_cards[i].keyname() != "GCOUNT" &&
            m_cards[i].keyname().find("NAXIS") == std::string::npos) {
            m_cards[i].write(FPTR(vptr));
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print FITS header information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing FITS header information.
 ***************************************************************************/
std::string GFitsHeader::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GFitsHeader ("+gammalib::str(size())+" cards) ===");

        // NORMAL: Append cards
        if (chatter >= NORMAL) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n"+m_cards[i].print());
            }
        }

    } // endif: chatter was not silent

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
    m_cards.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] header FITS header.
 ***************************************************************************/
void GFitsHeader::copy_members(const GFitsHeader& header)
{
    // Copy members
    m_cards = header.m_cards;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsHeader::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get index of header card
 *
 * @param[in] keyname Header card key name.
 * @return Index of header card (-1 if @p keyname is not found)
 *
 * Returns index of header card based on the @p keyname. If no header card
 * is found, -1 is returned.
 ***************************************************************************/
int GFitsHeader::get_index(const std::string& keyname) const
{
    // Initialise index
    int index = -1;

    // Search keyname in list
    for (int i = 0; i < size(); ++i) {
        if (m_cards[i].keyname() == keyname) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}
