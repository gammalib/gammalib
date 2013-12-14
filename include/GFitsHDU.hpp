/***************************************************************************
 *            GFitsHDU.hpp - Abstract FITS extension base class            *
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
 * @file GFitsHDU.hpp
 * @brief Abstract FITS extension base class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GFITSHDU_HPP
#define GFITSHDU_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFitsHeader.hpp"
#include "GFitsHeaderCard.hpp"


/***********************************************************************//**
 * @class GFitsHDU
 *
 * @brief Abstract FITS extension base class
 *
 * This class defines the abstract interface for a FITS extension, also known
 * as Header Data Unit (HDU). Each HDU consists of a header and a data area.
 * The header is composed of cards and is implemented by the GFitsHeader
 * class. The data are is either an image or a table which are implemented
 * as derived classes from GFitsHDU.
 *
 * @todo Implement GFitsHDU* select(const std::string& expr) that applies
 * to a table HDU for row selection.
 ***************************************************************************/
class GFitsHDU : public GBase {

    // Friend classes
    friend class GFits;

public:
    // Constructors and destructors
    GFitsHDU(void);
    GFitsHDU(const GFitsHDU& hdu);
    virtual ~GFitsHDU(void);

    // Operators
    GFitsHDU& operator=(const GFitsHDU& hdu);

    // Public enumerators
    enum HDUType {
        HT_IMAGE = 0,
        HT_ASCII_TABLE = 1,
        HT_BIN_TABLE = 2
    };

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GFitsHDU*   clone(void) const = 0;
    virtual HDUType     exttype(void) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented methods
    int                    size(void) const;
    const std::string&     extname(void) const;
    void                   extname(const std::string& extname);
    const int&             extno(void) const;
    void                   extno(const int& extno);
    GFitsHeader&           header(void);
    const GFitsHeader&     header(void) const;
    bool                   hascard(const int& cardno) const;
    bool                   hascard(const std::string& keyname) const;
    GFitsHeaderCard&       card(const int& cardno);
    const GFitsHeaderCard& card(const int& cardno) const;
    GFitsHeaderCard&       card(const std::string& keyname);
    const GFitsHeaderCard& card(const std::string& keyname) const;
    void                   card(const GFitsHeaderCard& card);
    void                   card(const std::string& keyname,
                                const std::string& value,
                                const std::string& comment);
    void                   card(const std::string& keyname,
                                const double& value,
                                const std::string& comment);
    void                   card(const std::string& keyname,
                                const int& value,
                                const std::string& comment);
    std::string            string(const std::string& keyname) const;
    double                 real(const std::string& keyname) const;
    int                    integer(const std::string& keyname) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GFitsHDU& hdu);
    void         free_members(void);
    void         connect(void* fptr);
    void         move_to_hdu(void);
    HDUType      get_hdu_type(void) const;
    void         open(void* vptr, int hdunum);
    void         save(void);
    std::string  print_hdu(const GChatter& chatter = NORMAL) const;
    std::string  typecode(int type) const;

    // Pure virtual protected methods
    virtual void data_open(void* vptr) = 0;
    virtual void data_save(void) = 0;
    virtual void data_close(void) = 0;
    virtual void data_connect(void* vptr) = 0;

    // Protected data area
    void*        m_fitsfile;    //!< FITS file pointer pointing on actual HDU
    int          m_hdunum;      //!< HDU number (starting from 0)
    std::string  m_name;        //!< HDU name (extname)
    GFitsHeader  m_header;      //!< HDU header
};


/***********************************************************************//**
 * @brief Return number of cards in HDU header
 *
 * @return Number of cards in HDU header.
 *
 * Returns the number of cards in the header of the FITS extension.
 ***************************************************************************/
inline
int GFitsHDU::size(void) const
{
    return (m_header.size());
}


/***********************************************************************//**
 * @brief Return extension name
 *
 * @return Extension name.
 *
 * Returns the name of the FITS extension.
 ***************************************************************************/
inline
const std::string& GFitsHDU::extname(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Return extension number
 *
 * @return Extension number (starting from 0 for the primary image).
 *
 * Returns the number of the FITS extension.
 ***************************************************************************/
inline
const int& GFitsHDU::extno(void) const
{
    return (m_hdunum);
}


/***********************************************************************//**
 * @brief Set extension number
 *
 * @param[in] extno Extension number (starting from 0 for the primary image).
 *
 * Set the number of the FITS extension.
 ***************************************************************************/
inline
void GFitsHDU::extno(const int& extno)
{
    m_hdunum = extno;
    return;
}


/***********************************************************************//**
 * @brief Return extension header
 *
 * @return Extension header.
 *
 * Returns the extension header.
 ***************************************************************************/
inline
GFitsHeader& GFitsHDU::header(void)
{
    return (m_header);
}


/***********************************************************************//**
 * @brief Return extension header (const version)
 *
 * @return Extension header.
 *
 * Returns the extension header.
 ***************************************************************************/
inline
const GFitsHeader& GFitsHDU::header(void) const
{
    return (m_header);
}


/***********************************************************************//**
 * @brief Check existence of header card
 *
 * @param[in] cardno Number of card in header.
 * @return True if card with specified card number @p cardno exists in header.
 ***************************************************************************/
inline
bool GFitsHDU::hascard(const int& cardno) const
{
    return m_header.contains(cardno);
}


/***********************************************************************//**
 * @brief Checks for presence of header card
 *
 * @param[in] keyname Name of header card.
 * @return True if card with specified @p keyname exists in header.
 ***************************************************************************/
inline
bool GFitsHDU::hascard(const std::string& keyname) const
{
    return (m_header.contains(keyname));
}


/***********************************************************************//**
 * @brief Return header card
 *
 * @param[in] cardno Number of card in header.
 * @return Header card.
 ***************************************************************************/
inline
GFitsHeaderCard& GFitsHDU::card(const int& cardno)
{
    return (m_header.at(cardno));
}


/***********************************************************************//**
 * @brief Return header card (const version)
 *
 * @param[in] cardno Number of card in header.
 * @return Header card.
 ***************************************************************************/
inline
const GFitsHeaderCard& GFitsHDU::card(const int& cardno) const
{
    return (m_header.at(cardno));
}


/***********************************************************************//**
 * @brief Return header card
 *
 * @param[in] keyname Name of header card.
 * @return Header card.
 ***************************************************************************/
inline
GFitsHeaderCard& GFitsHDU::card(const std::string& keyname)
{
    return (m_header.at(keyname));
}


/***********************************************************************//**
 * @brief Return header card (const version)
 *
 * @param[in] keyname Name of header card.
 * @return Header card.
 ***************************************************************************/
inline
const GFitsHeaderCard& GFitsHDU::card(const std::string& keyname) const
{
    return (m_header.at(keyname));
}


/***********************************************************************//**
 * @brief Append or update header card
 *
 * @param[in] card Header card.
 ***************************************************************************/
inline
void GFitsHDU::card(const GFitsHeaderCard& card)
{
    m_header.append(card);
    return;
}


/***********************************************************************//**
 * @brief Append or update string value header card
 *
 * @param[in] keyname Name of the header card.
 * @param[in] value String value of the header card.
 * @param[in] comment Comment of the header card.
 ***************************************************************************/
inline
void GFitsHDU::card(const std::string& keyname, const std::string& value,
                    const std::string& comment)
{
    m_header.append(GFitsHeaderCard(keyname, value, comment));
    return;
}


/***********************************************************************//**
 * @brief Append or update double precision value header card
 *
 * @param[in] keyname Name of the header card.
 * @param[in] value Double precision value of the header card.
 * @param[in] comment Comment of the header card.
 ***************************************************************************/
inline
void GFitsHDU::card(const std::string& keyname, const double& value,
                    const std::string& comment)
{
    m_header.append(GFitsHeaderCard(keyname, value, comment));
    return;
}


/***********************************************************************//**
 * @brief Append or update integer value header card
 *
 * @param[in] keyname Name of the header card.
 * @param[in] value Integer value of the header card.
 * @param[in] comment Comment of the header card.
 ***************************************************************************/
inline
void GFitsHDU::card(const std::string& keyname, const int& value,
                    const std::string& comment)
{
    m_header.append(GFitsHeaderCard(keyname,  value, comment));
    return;
}


/***********************************************************************//**
 * @brief Return card value as string
 *
 * @param[in] keyname Name of header card.
 * @return String card value.
 ***************************************************************************/
inline
std::string GFitsHDU::string(const std::string& keyname) const
{
    return (m_header.at(keyname).string());
}


/***********************************************************************//**
 * @brief Return card value as double precision
 *
 * @param[in] keyname Name of header card.
 * @return Double precision card value.
 ***************************************************************************/
inline
double GFitsHDU::real(const std::string& keyname) const
{
    return (m_header.at(keyname).real());
}


/***********************************************************************//**
 * @brief Return card value as integer
 *
 * @param[in] keyname Name of header card.
 * @return Integer card value.
 ***************************************************************************/
inline
int GFitsHDU::integer(const std::string& keyname) const
{
    return (m_header.at(keyname).integer());
}

#endif /* GFITSHDU_HPP */
