/***************************************************************************
 *                     fits module  -  Python bindings                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
%module fits
%feature("autodoc", "1");

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="support") "GException.i";

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
