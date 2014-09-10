/***************************************************************************
 *                       GFits.hpp - FITS file class                       *
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
 * @file GFits.i
 * @brief FITS file class interface definition
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


/***********************************************************************//**
 * @class GFits
 *
 * @brief FITS file class
 ***************************************************************************/
class GFits : public GContainer {
public:
    // Constructors and destructors
    GFits(void);
    explicit GFits(const std::string& filename, const bool& create = false);
    GFits(const GFits& fits);
    virtual ~GFits(void);

    // Methods
    void               clear(void);
    GFits*             clone(void) const;
    GFitsHDU*          at(const int& extno);
    GFitsHDU*          at(const std::string& extname);
    GFitsImage*        image(const int& extno);
    GFitsImage*        image(const std::string& extname);
    GFitsTable*        table(const int& extno);
    GFitsTable*        table(const std::string& extname);
    int                size(void) const;
    bool               is_empty(void) const;
    GFitsHDU*          set(const int& extno, const GFitsHDU& hdu);
    GFitsHDU*          set(const std::string& extname, const GFitsHDU& hdu);
    GFitsHDU*          append(const GFitsHDU& hdu);
    GFitsHDU*          insert(const int& extno, const GFitsHDU& hdu);
    GFitsHDU*          insert(const std::string& extname, const GFitsHDU& hdu);
    void               remove(const int& extno);
    void               remove(const std::string& extname);
    void               reserve(const int& num);
    void               extend(const GFits& fits);
    bool               contains(const int& extno) const;
    bool               contains(const std::string& extname) const;
    const std::string& filename(void) const;
    int                extno(const std::string& extname) const;
    void               open(const std::string& filename,
                            const bool&        create = false);
    void               save(const bool& clobber = false);
    void               saveto(const std::string& filename,
                              const bool&        clobber = false);
    void               close(void);
};


/***********************************************************************//**
 * @brief GFits class SWIG extension
 ***************************************************************************/
%extend GFits {
    GFitsHDU* __getitem__(const int& extno) {
        if (extno >= 0 && extno < self->size()) {
            return (*self)[extno];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Extension number",
                                           extno, self->size());
        }
    }
    GFitsHDU* __getitem__(const std::string& extname) {
        return (*self)[extname];
    }
    void __setitem__(const int& extno, const GFitsHDU& hdu) {
        if (extno >= 0 && extno < self->size()) {
            self->set(extno, hdu);
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Extension number",
                                           extno, self->size());
        }
    }
    void __setitem__(const std::string& extname, const GFitsHDU& hdu) {
        self->set(extname, hdu);
        return;
    }
    GFits copy() {
        return (*self);
    }
}
