/***************************************************************************
 *                   GFitsHDU.i  - FITS HDU handling class                 *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2011 by Jurgen Knodlseder                         *
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
 * @brief FITS HDU class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHDU.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GFitsHDU
 *
 * @brief FITS Header Data Unit (HDU) class
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
    GFitsHDU(const GFitsHDU& hdu);
    virtual ~GFitsHDU(void);

    // Public enumerators
    enum HDUType {
        HT_IMAGE = 0,
        HT_ASCII_TABLE = 1,
        HT_BIN_TABLE = 2
    };

    // Pure virtual methods
    virtual HDUType     exttype(void) const = 0;
    virtual GFitsHDU*   clone(void) const = 0;

    // Implemented methods
    virtual std::string      extname(void) const;
    virtual void             extname(const std::string& extname);
    virtual int              extno(void) const;
    virtual void             extno(int num);
    virtual GFitsHeader*     header(void);
    virtual GFitsHeaderCard* card(const std::string& keyname);
    virtual GFitsHeaderCard* card(const int& cardno);
    virtual int              cards(void) const;
    virtual std::string      string(const std::string& keyname) const;
    virtual double           real(const std::string& keyname) const;
    virtual int              integer(const std::string& keyname) const;
    virtual void             card(const std::string& keyname,
                                  const std::string& value,
                                  const std::string& comment);
    virtual void             card(const std::string& keyname,
                                  const double&      value,
                                  const std::string& comment);
    virtual void             card(const std::string& keyname,
                                  const int&         value,
                                  const std::string& comment);
};


/***********************************************************************//**
 * @brief GFitsHDU class SWIG extension
 ***************************************************************************/
%extend GFitsHDU {
    char *__str__() {
        return tochar(self->print());
    }
};
