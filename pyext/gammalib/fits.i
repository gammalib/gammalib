/***************************************************************************
 *                          fits.i - FITS module                           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall fits.i                                          *
 ***************************************************************************/
/**
 * @file fits.i
 * @brief FITS module
 * @author Juergen Knoedlseder
 */
%module fits
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
#include "GException.hpp"
#include "GTools.hpp"
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.base") "GContainer.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Typemaps ___________________________________________________________ */
%typemap(out) GFitsHDU* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GFitsHDU;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GFitsImage* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GFitsImage;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GFitsTable* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GFitsTable;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GFitsTableCol* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GFitsTableCol;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}

/* __ FITS module ________________________________________________________ */
%include "GFits.i"
%include "GFitsHDU.i"
%include "GFitsHeader.i"
%include "GFitsHeaderCard.i"
%include "GFitsImage.i"
%include "GFitsImageByte.i"
%include "GFitsImageSByte.i"
%include "GFitsImageUShort.i"
%include "GFitsImageShort.i"
%include "GFitsImageULong.i"
%include "GFitsImageLong.i"
%include "GFitsImageLongLong.i"
%include "GFitsImageFloat.i"
%include "GFitsImageDouble.i"
%include "GFitsTable.i"
%include "GFitsAsciiTable.i"
%include "GFitsBinTable.i"
%include "GFitsTableCol.i"
%include "GFitsTableBitCol.i"
%include "GFitsTableByteCol.i"
%include "GFitsTableBoolCol.i"
%include "GFitsTableStringCol.i"
%include "GFitsTableUShortCol.i"
%include "GFitsTableShortCol.i"
%include "GFitsTableULongCol.i"
%include "GFitsTableLongCol.i"
%include "GFitsTableLongLongCol.i"
%include "GFitsTableFloatCol.i"
%include "GFitsTableDoubleCol.i"
%include "GFitsTableCFloatCol.i"
%include "GFitsTableCDoubleCol.i"
