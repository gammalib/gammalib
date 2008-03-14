/***************************************************************************
 *     GHealpix.i  -  Healpix sky representation class SWIG definition     *
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
 * @file GHealpix.i
 * @brief GHealpix class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GHealpix.hpp"
%}


/***********************************************************************//**
 * @class GHealpix
 *
 * @brief GHealpix class interface defintion
 ***************************************************************************/
class GHealpix {
public:
    // Constructors and destructors
    GHealpix();
    GHealpix(const GFitsHDU* hdu);
    GHealpix(const GHealpix& pixels);
    virtual ~GHealpix();

    // Methods
    void    read(const GFitsHDU* hdu);
    int     nside(void) const;
    int     num_pixels(void) const;
    double  omega(void) const;
    GSkyDir pix2ang(const int& ipix);
};


/***********************************************************************//**
 * @brief GHealpix class extension
 ***************************************************************************/
%extend GHealpix {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
	    str_buffer[1000] = '\0';
	    return str_buffer;
    }
    GHealpix copy() {
        return (*self);
    }
};
