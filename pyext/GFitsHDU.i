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
//%feature("notabstract") GFitsHDU;


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
    GFitsHDU(const GFitsHDU& hdu);
    virtual ~GFitsHDU(void);

    // Public enumerators
    enum HDUType {
        HT_IMAGE = 0,
        HT_ASCII_TABLE = 1,
        HT_BIN_TABLE = 2
    };

    // Pure virtual methods
    virtual HDUType   exttype(void) const = 0;
    virtual GFitsHDU* clone(void) const = 0;

    // Implemented methods
    virtual std::string      extname(void) const { return m_name; }
    virtual void             extname(const std::string& extname);
    virtual int              extno(void) const { return m_hdunum; }
    virtual void             extno(int num) { m_hdunum=num; }
    virtual GFitsHeader*     header(void) { return &m_header; }
    virtual GFitsHeaderCard* card(const std::string& keyname);
    virtual GFitsHeaderCard* card(const int& cardno);
    virtual int              cards(void) const { return m_header.size(); }
    virtual std::string      string(const std::string& keyname) const;
    virtual double           real(const std::string& keyname) const;
    virtual int              integer(const std::string& keyname) const;
    virtual void             card(const std::string& keyname, const std::string& value,
                                  const std::string& comment);
    virtual void             card(const std::string& keyname, const double& value,
                                  const std::string& comment);
    virtual void             card(const std::string& keyname, const int& value,
                                  const std::string& comment);
    virtual GFitsHDU*        primary(void);
};


/***********************************************************************//**
 * @brief GFitsHDU class SWIG extension
 ***************************************************************************/
