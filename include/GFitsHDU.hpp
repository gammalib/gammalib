/***************************************************************************
 *                  GFitsHDU.hpp  - FITS HDU handling class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsHDU.hpp
 * @brief GFitsHDU class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSHDU_HPP
#define GFITSHDU_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHeader.hpp"
#include "GFitsHeaderCard.hpp"


/***********************************************************************//**
 * @class GFitsHDU
 *
 * @brief Implements the FITS Header Data Unit (HDU) interface
 *
 * The HDU is the basic unit of a FITS file. Each HDU consists of a header
 * and a data area. The header is composed of cards and is implemented by
 * the GFitsHeader class. The data are is either an image or a table and
 * is implemented by the abstract GFitsData base class.
 * The class holds a copy of the original FITS
 *
 * @todo Implement GFitsHDU* select(const std::string& expr) that applies
 * to a table HDU for row selection.
 ***************************************************************************/
class GFitsHDU {

    // Friend classes
    friend class GFits;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFitsHDU& hdu);

public:
    // Constructors and destructors
    GFitsHDU(void);
    GFitsHDU(const GFitsHDU& hdu);
    virtual ~GFitsHDU(void);

    // Operators
    GFitsHDU& operator= (const GFitsHDU& hdu);

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

protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const GFitsHDU& hdu);
    void    free_members(void);
    void    connect(void* fptr);
    void    move_to_hdu(void);
    HDUType get_hdu_type(void) const;
    void    open(void* vptr, int hdunum);
    void    save(void);

    // Pure virtual protected methods
    virtual void data_open(void* vptr) = 0;
    virtual void data_save(void) = 0;
    virtual void data_close(void) = 0;
    virtual void data_connect(void* vptr) = 0;

    // Protected data area
    void*        m_fitsfile;    //!< FITS file pointer pointing on actual HDU
    int          m_hdunum;      //!< HDU number (starting from 0)
    std::string  m_name;        //!< HDU name (extname)
    GFitsHeader  m_header;      //!< HDU header
};

#endif /* GFITSHDU_HPP */
