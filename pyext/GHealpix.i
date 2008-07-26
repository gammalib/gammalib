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
    GHealpix(int nside, int scheme = 1, int coordsys = 1, int dimension = 1);
    GHealpix(const GFitsHDU* hdu);
    GHealpix(const GHealpix& pixels);
    virtual ~GHealpix();

    // Methods
    void    read(const GFitsHDU* hdu);
    void    write(GFits* fits);
    int     nside(void) const;
    int     npix(void) const;
    double  omega(void) const;
    GSkyDir pix2ang(const int& ipix);
    int     ang2pix(GSkyDir dir) const;
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
    double __getitem__(int pixel) {
	    if (pixel >= 0 && pixel < (int)self->npix())
            return (*self)(pixel);
        else
	        throw GException::out_of_range("__getitem__(int)", pixel, (int)self->npix());
    }
    void __setitem__(int pixel, const double val) {
	    if (pixel >= 0 && pixel < (int)self->npix())
	        (*self)(pixel) = val;
        else
            throw GException::out_of_range("__setitem__(int)", pixel, (int)self->npix());
    }
    GHealpix copy() {
        return (*self);
    }
};
