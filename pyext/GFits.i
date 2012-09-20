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

/* __ Typemaps ____________________________________________________________*/
%typemap(out) GFitsHDU* {
    if (dynamic_cast<GFitsImage*>($1) != NULL) {
        if (dynamic_cast<GFitsImageByte*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageByte, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageDouble*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageDouble, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageFloat*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageFloat, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageLong*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageFloat, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageLongLong*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageLongLong, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageSByte*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageSByte, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageShort*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageShort, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageULong*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageULong, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageUShort*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageUShort, 0 |  0 );
        }
        else {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImage, 0 |  0 );
        }
    }
    else if (dynamic_cast<GFitsTable*>($1) != NULL) {
        if (dynamic_cast<GFitsBinTable*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsBinTable, 0 |  0 );
        }
        else if (dynamic_cast<GFitsAsciiTable*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsAsciiTable, 0 |  0 );
        }
        else {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsTable, 0 |  0 );
        }
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsHDU, 0 |  0 );
    }
}
%typemap(out) GFitsImage* {
    if (dynamic_cast<GFitsImageByte*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageByte, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageDouble*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageDouble, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageFloat*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageFloat, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageLong*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageFloat, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageLongLong*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageLongLong, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageSByte*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageSByte, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageShort*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageShort, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageULong*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageULong, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageUShort*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageUShort, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImage, 0 |  0 );
    }
}
%typemap(out) GFitsTable* {
    if (dynamic_cast<GFitsBinTable*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsBinTable, 0 |  0 );
    }
    else if (dynamic_cast<GFitsAsciiTable*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsAsciiTable, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsTable, 0 |  0 );
    }
}


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
