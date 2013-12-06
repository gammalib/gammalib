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

/* __ Includes ___________________________________________________________ */
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
    bool               isempty(void) const;
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
