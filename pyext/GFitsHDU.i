/***************************************************************************
 *             GFitsHDU.i  - FITS HDU handling class SWIG file             *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsHDU.i
 * @brief GFitsHDU class SWIG file
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHDU.hpp"
%}

%include stl.i


/***********************************************************************//**
 * @class GFitsHDU
 *
 * @brief Implements the FITS Header Data Unit (HDU) SWIG interface
 *
 * The HDU is the basic unit of a FITS file. Each HDU consists of a header
 * and a data area. The header is composed of cards and is implemented by
 * the GFitsHeader class. The data are is either an image or a table and
 * is implemented by the abstract GFitsData base class.
 ***************************************************************************/
class GFitsHDU {
public:
    // Constructors and destructors
    GFitsHDU(void);
    GFitsHDU(const GFitsImage& image);
    GFitsHDU(const GFitsAsciiTable& table);
    GFitsHDU(const GFitsBinTable& table);
    GFitsHDU(const GFitsHDU& hdu);
    virtual ~GFitsHDU(void);

    // Methods
    std::string      extname(void) const;
    void             extname(const std::string& extname);
    int              extno(void) const;
    int              exttype(void) const;
    GFitsHeader*     header(void) const;
    GFitsData*       data(void) const;
    GFitsHeaderCard* card(const std::string& keyname) const;
    GFitsHeaderCard* card(const int& cardno) const;
    void             card(const std::string& keyname, const std::string& value,
                          const std::string& comment);
    void             card(const std::string& keyname, const double& value,
                          const std::string& comment);
    void             card(const std::string& keyname, const int& value,
                          const std::string& comment);
    GFitsTableCol*   column(const std::string& colname) const;
    GFitsTableCol*   column(const int& colnum) const;
    void             primary(void);
};


/***********************************************************************//**
 * @brief GFitsHDU class SWIG extension
 ***************************************************************************/
%extend GFitsHDU {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
        str_buffer[100000] = '\0';
        return str_buffer;
    }
    GFitsHDU copy() {
        return (*self);
    }
}
