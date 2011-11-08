/***************************************************************************
 *                  GFitsHDU.hpp  - FITS HDU handling class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
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
 * @brief GFitsHDU class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSHDU_HPP
#define GFITSHDU_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFitsHeader.hpp"
#include "GFitsHeaderCard.hpp"


/***********************************************************************//**
 * @class GFitsHDU
 *
 * @brief Implements the FITS Header Data Unit (HDU) interface
 *
 * The HDU is the basic unit of a FITS file. Each HDU consists of a header
 * and a data area. The header is composed of cards and is implemented by
 * the GFitsHeader class. The data are is either an image or a table and
 * is implemented by the abstract GFitsData base class.
 * The class holds a copy of the original FITS
 *
 * @todo Implement GFitsHDU* select(const std::string& expr) that applies
 * to a table HDU for row selection.
 ***************************************************************************/
class GFitsHDU {

    // Friend classes
    friend class GFits;

public:
    // Constructors and destructors
    GFitsHDU(void);
    GFitsHDU(const GFitsHDU& hdu);
    virtual ~GFitsHDU(void);

    // Operators
    GFitsHDU& operator= (const GFitsHDU& hdu);

    // Public enumerators
    enum HDUType {
        HT_IMAGE = 0,
        HT_ASCII_TABLE = 1,
        HT_BIN_TABLE = 2
    };

    // Pure virtual methods
    virtual HDUType     exttype(void) const = 0;
    virtual GFitsHDU*   clone(void) const = 0;
    virtual std::string print(void) const = 0;

    // Implemented methods
    int              size(void) const { return m_header.size(); }
    std::string      extname(void) const { return m_name; }
    void             extname(const std::string& extname);
    int              extno(void) const { return m_hdunum; }
    void             extno(int num) { m_hdunum=num; }
    GFitsHeader*     header(void) { return &m_header; }
    bool             hascard(const std::string& keyname) const;
    bool             hascard(const int& cardno) const;
    GFitsHeaderCard* card(const std::string& keyname);
    GFitsHeaderCard* card(const int& cardno);
    std::string      string(const std::string& keyname) const;
    double           real(const std::string& keyname) const;
    int              integer(const std::string& keyname) const;
    void             card(const std::string& keyname, const std::string& value,
                          const std::string& comment);
    void             card(const std::string& keyname, const double& value,
                          const std::string& comment);
    void             card(const std::string& keyname, const int& value,
                          const std::string& comment);

protected:
    // Protected methods
    void         connect(void* fptr);
    void         move_to_hdu(void);
    HDUType      get_hdu_type(void) const;
    void         open(void* vptr, int hdunum);
    void         save(void);
    std::string  print_hdu(void) const;
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

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsHDU& hdu);
    void free_members(void);
};

#endif /* GFITSHDU_HPP */
