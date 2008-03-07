/***************************************************************************
 *           GFitsHeader.i  - FITS header handling class SWIG file         *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsHeader.i
 * @brief GFitsHeader class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHeader.hpp"
%}

%include stl.i


/***********************************************************************//**
 * @class GFitsHeader
 *
 * @brief Implements FITS header class SWIG interface
 *
 * The FITS header class contains all cards that are found in the header of
 * a HDU. All cards will be hold in memory, so no link to a FITS file is
 * required. Cards may be read from one file (using the 'open' method) and
 * saved into another file (using the 'save' method). Cards are added or
 * changed using the 'update' method.
 ***************************************************************************/
class GFitsHeader {
public:
    // Constructors and destructors
    GFitsHeader();
    GFitsHeader(const GFitsHeader& header);
    ~GFitsHeader();

    // Methods
    void             update(const GFitsHeaderCard& card);
    GFitsHeaderCard* card(const std::string& keyname);
    GFitsHeaderCard* card(const int& cardno);
    std::string      string(const std::string& keyname);
    std::string      string(const int& cardno);
    double           real(const std::string& keyname);
    double           real(const int& cardno);
    int              integer(const std::string& keyname);
    int              integer(const int& cardno);
    GFitsHeader*     clone(void) const;
    int              num_cards(void) const;
};


/***********************************************************************//**
 * @brief GFitsHeader class SWIG extension
 ***************************************************************************/
%extend GFitsHeader {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
        str_buffer[100000] = '\0';
        return str_buffer;
    }
    GFitsHeader copy() {
        return (*self);
    }
}
