/***************************************************************************
 *                    GFits.hpp  - FITS file access class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Juergen Knoedlseder                         *
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
 * @file GFits.hpp
 * @brief FITS file access class interface definition
 * @author J. Knoedlseder
 */

#ifndef GFITS_HPP
#define GFITS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <iostream>
#include "GLog.hpp"
#include "GFitsHDU.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GFits
 *
 * @brief FITS file access class
 *
 * GFits is the basic FITS file interface. All FITS file handlings operate
 * via members of GFits. A FITS file is composed of Header Data Units (HDU)
 * which are implemented by the GFitsHDU class. Each HDU is composed of a
 * header (GFitsHeader) and some data (GFitsData).
 ***************************************************************************/
class GFits {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GFits& fits);
    friend GLog&         operator<<(GLog& log,        const GFits& fits);

public:
    // Constructors and destructors
    GFits(void);
    explicit GFits(const std::string& filename, bool create = false);
    GFits(const GFits& fits);
    virtual ~GFits(void);

    // Operators
    GFits& operator=(const GFits& fits);

    // Methods
    void        clear(void);
    int         size(void) const;
    void        open(const std::string& filename, bool create = false);
    void        save(bool clobber = false);
    void        saveto(const std::string& filename, bool clobber = false);
    void        close(void);
    void        append(const GFitsHDU& hdu);
    bool        hashdu(const std::string& extname) const;
    bool        hashdu(int extno) const;
    GFitsHDU*   hdu(const std::string& extname) const;
    GFitsHDU*   hdu(int extno) const;
    GFitsImage* image(const std::string& extname) const;
    GFitsImage* image(int extno) const;
    GFitsTable* table(const std::string& extname) const;
    GFitsTable* table(int extno) const;
    std::string name(void) const { return m_filename; }
    std::string print(void) const;

    // Complex single precision type
    typedef struct {
        float re;
        float im;
    } cfloat;

    // Complex double precision type
    typedef struct {
        float re;
        float im;
    } cdouble;

private:
    // Private methods
    void        init_members(void);
    void        copy_members(const GFits& fits);
    void        free_members(void);
    GFitsImage* new_image(void);
    GFitsImage* new_primary(void);

    // Private data area
    std::vector<GFitsHDU*> m_hdu;        //!< Pointers to HDUs
    std::string            m_filename;   //!< FITS file name
    void*                  m_fitsfile;   //!< FITS file pointer
    bool                   m_readwrite;  //!< FITS file is readwrite (true/false)
    bool                   m_created;    //!< FITS file has been created (true/false)
};

#endif /* GFITS_HPP */
