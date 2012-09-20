/***************************************************************************
 *                     GFits.i  - FITS file access class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @file GFits.i
 * @brief Fits file access class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsImage.hpp"
#include "GFitsImageByte.hpp"
#include "GFitsImageDouble.hpp"
#include "GFitsImageFloat.hpp"
#include "GFitsImageLong.hpp"
#include "GFitsImageLongLong.hpp"
#include "GFitsImageSByte.hpp"
#include "GFitsImageShort.hpp"
#include "GFitsImageULong.hpp"
#include "GFitsImageUShort.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsAsciiTable.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GFits
 *
 * @brief FITS file access class
 *
 * GFits is the basic FITS file interface. All FITS file handlings operate
 * via members of GFits. A FITS file is composed of Header Data Units (HDU)
 * which are implemented by the GFitsHDU class.
 ***************************************************************************/
class GFits {

public:
    // Constructors and destructors
    GFits(void);
    explicit GFits(const std::string& filename, bool create = false);
    GFits(const GFits& fits);
    virtual ~GFits(void);

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
    std::string name(void) const;
};


/***********************************************************************//**
 * @brief GFits class SWIG extension
 ***************************************************************************/
%extend GFits {
    char *__str__() {
        return tochar(self->print());
    }
    GFits copy() {
        return (*self);
    }
}
