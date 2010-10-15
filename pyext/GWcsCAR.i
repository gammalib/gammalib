/***************************************************************************
 *             GWcsCAR.i  -  Cartesian projection class SWIG file          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GWcsCAR.i
 * @brief GWcsCAR class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcsCAR.hpp"
%}


/***********************************************************************//**
 * @class GWcsCAR
 *
 * @brief GWcsCAR class interface defintion
 ***************************************************************************/
class GWcsCAR : public GWcs {
public:
    // Constructors and destructors
    GWcsCAR(void);
    explicit GWcsCAR(const std::string& coords,
                     const double& crval1, const double& crval2,
                     const double& crpix1, const double& crpix2,
                     const double& cdelt1, const double& cdelt2,
                     const GMatrix& cd, const GVector& pv2);
    explicit GWcsCAR(const GFitsHDU* hdu);
    GWcsCAR(const GWcsCAR& wcs);
    virtual ~GWcsCAR(void);

    // Implemented pure virtual methods
    void clear(void);
    void read(const GFitsHDU* hdu);
    void write(GFitsHDU* hdu) const;

    // Overloaded base class methods
    double omega(const GSkyPixel& pix) const;
};


/***********************************************************************//**
 * @brief GWcsCAR class extension
 ***************************************************************************/
%extend GWcsCAR {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
        str_buffer[1000] = '\0';
        return str_buffer;
    }
    GWcsCAR copy() {
        return (*self);
    }
};
