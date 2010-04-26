/***************************************************************************
 *                GSkymap.i  -  Sky map class SWIG definition              *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2010 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GSkymap.i
 * @brief GSkymap class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkymap.hpp"
%}
%include stl.i


/***********************************************************************//**
 * @class GSkymap
 *
 * @brief GSkymap class interface defintion
 ***************************************************************************/
class GSkymap {
public:
    // Constructors and destructors
    GSkymap(void);
    GSkymap(const std::string& wcs, const std::string& coords,
            const int& nside, const std::string& order,
            const int nmaps = 1);
    GSkymap(const std::string& wcs, const std::string& coords,
            GSkyDir& dir, const int& nlon, const int& nlat,
            const double& dlon, const double& dlat, const int nmaps = 1);
    GSkymap(const GSkymap& map);
    ~GSkymap(void);

    // Methods
    void    load(const std::string& filename);
    void    save(const std::string& filename, int clobber = 0);
    void    read(const GFitsHDU* hdu);
    void    write(GFitsHDU* hdu);
    GSkyDir pix2dir(const int& ipix);
    int     dir2pix(GSkyDir dir) const;
    double  omega(const int& pix) const;
    int     npix(void) const;
};


/***********************************************************************//**
 * @brief GSkymap class extension
 ***************************************************************************/
%extend GSkymap {
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
    GSkymap copy() {
        return (*self);
    }
};
