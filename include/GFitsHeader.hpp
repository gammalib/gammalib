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
 * @file GFitsHeader.hpp
 * @brief FITS header cards container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSHEADER_HPP
#define GFITSHEADER_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GFitsHeaderCard.hpp"


/***********************************************************************//**
 * @class GFitsHeader
 *
 * @brief Interface for FITS header class
 *
 * The FITS header class is a container class for header cards.
 *
 * All cards of a FITS file extension will be held in memory, so no link to
 * a FITS file is required. Cards are read from a file using the load()
 * method, and cards are saved into a file using the save() method.
 ***************************************************************************/
class GFitsHeader : public GContainer {

public:
    // Constructors and destructors
    GFitsHeader(void);
    GFitsHeader(const GFitsHeader& header);
    virtual ~GFitsHeader(void);

    // Operators
    GFitsHeader&           operator=(const GFitsHeader& header);
    GFitsHeaderCard&       operator[](const int& cardno);
    const GFitsHeaderCard& operator[](const int& cardno) const;
    GFitsHeaderCard&       operator[](const std::string& keyname);
    const GFitsHeaderCard& operator[](const std::string& keyname) const;

    // Methods
    void                   clear(void);
    GFitsHeader*           clone(void) const;
    GFitsHeaderCard&       at(const int& cardno);
    const GFitsHeaderCard& at(const int& cardno) const;
    GFitsHeaderCard&       at(const std::string& keyname);
    const GFitsHeaderCard& at(const std::string& keyname) const;
    std::string            string(const int& cardno) const;
    std::string            string(const std::string& keyname) const;
    double                 real(const int& cardno) const;
    double                 real(const std::string& keyname) const;
    int                    integer(const int& cardno) const;
    int                    integer(const std::string& keyname) const;
    int                    size(void) const;
    bool                   is_empty(void) const;
    GFitsHeaderCard&       append(const GFitsHeaderCard& card);
    GFitsHeaderCard&       insert(const int& cardno, const GFitsHeaderCard& card);
    GFitsHeaderCard&       insert(const std::string& keyname, const GFitsHeaderCard& card);
    void                   remove(const int& cardno);
    void                   remove(const std::string& keyname);
    void                   reserve(const int& num);
    void                   extend(const GFitsHeader& header);
    bool                   contains(const int& cardno) const;
    bool                   contains(const std::string& keyname) const;
    void                   load(void* vptr);
    void                   save(void* vptr) const;
    std::string            print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsHeader& header);
    void free_members(void);
    int  get_index(const std::string& keyname) const;

    // Private data area
    std::vector<GFitsHeaderCard> m_cards; //!< Header cards
};


/***********************************************************************//**
 * @brief Return header card
 *
 * @param[in] cardno Number of card in header [0,...,size()-1]
 * @return Header card.
 ***************************************************************************/
inline
GFitsHeaderCard& GFitsHeader::operator[](const int& cardno)
{
    return (m_cards[cardno]);
}


/***********************************************************************//**
 * @brief Return pointer to model (const version)
 *
 * @param[in] cardno Number of card in header [0,...,size()-1]
 * @return Header card.
 ***************************************************************************/
inline
const GFitsHeaderCard& GFitsHeader::operator[](const int& cardno) const
{
    return (m_cards[cardno]);
}


/***********************************************************************//**
 * @brief Return header card
 *
 * @param[in] keyname Name of header card
 * @return Header card.
 ***************************************************************************/
inline
GFitsHeaderCard& GFitsHeader::operator[](const std::string& keyname)
{
    return (at(keyname));
}


/***********************************************************************//**
 * @brief Return header card (const version)
 *
 * @param[in] keyname Name of header card
 * @return Header card.
 ***************************************************************************/
inline
const GFitsHeaderCard& GFitsHeader::operator[](const std::string& keyname) const
{
    return (at(keyname));
}


/***********************************************************************//**
 * @brief Return number of cards in header
 *
 * @return Number of cards in header.
 *
 * Returns the number of cards in the extension header.
 ***************************************************************************/
inline
int GFitsHeader::size(void) const
{
    return (m_cards.size());
}


/***********************************************************************//**
 * @brief Signals if there are no cards in the FITS header
 *
 * @return True if there are no cards in the FITS header, false otherwise
 *
 * Signals if there are no cards in the FITS header.
 ***************************************************************************/
inline
bool GFitsHeader::is_empty(void) const
{
    return (m_cards.empty());
}


/***********************************************************************//**
 * @brief Reserves space for cards in FITS header
 *
 * @param[in] num Number of cards
 *
 * Reserves space for @p num cards in the FITS header.
 ***************************************************************************/
inline
void GFitsHeader::reserve(const int& num)
{
    m_cards.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Check if card is present in header
 *
 * @param[in] cardno Number of card in header.
 * @return True of card exists, false otherwise.
 *
 * Signals whether a card with specified card number exists in header.
 ***************************************************************************/
inline
bool GFitsHeader::contains(const int& cardno) const
{
    return (cardno >= 0 && cardno < size());
}


/***********************************************************************//**
 * @brief Check if card is present in header
 *
 * @param[in] keyname Name of header card.
 * @return True of card exists, false otherwise.
 *
 * Signals whether a card with specified @p keyname exists in header.
 ***************************************************************************/
inline
bool GFitsHeader::contains(const std::string& keyname) const
{
    return (get_index(keyname) != -1);
}

#endif /* GFITSHEADER_HPP */
