/***************************************************************************
 *                 GFits.i  - FITS file access class SWIG file             *
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
 * @file GFits.i
 * @brief GFits class SWIG definition.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFits.hpp"
%}

%include stl.i

/***********************************************************************//**
 * @class GFits
 *
 * @brief Implements FITS file SWIG interface
 *
 * GFits is the basic FITS file interface. All FITS file handlings operate
 * via members of GFits. A FITS file is composed of Header Data Units (HDU)
 * which are implemented by the GFitsHDU class.
 ***************************************************************************/
class GFits {
public:
    // Constructors and destructors
    GFits();
    GFits(const std::string& filename);
    GFits(const GFits& fits);
    ~GFits();

    // Operators
    //GFits& operator= (const GFits& fits);

    // Methods
    void      open(const std::string& filename);
    void      append_hdu(const GFitsHDU& hdu);
    void      save(void);
    void      saveto(const std::string& filename, int clobber = 0);
    void      close(void);
    GFitsHDU* hdu(const std::string& extname);
    GFitsHDU* hdu(int extno);
    int       num_hdus(void) const;

};


/***********************************************************************//**
 * @brief GFits class SWIG extension
 ***************************************************************************/
%extend GFits {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
        str_buffer[100000] = '\0';
        return str_buffer;
    }
    GFits copy() {
        return (*self);
    }
}
