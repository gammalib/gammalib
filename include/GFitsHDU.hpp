/***************************************************************************
 *                  GFitsHDU.hpp  - FITS HDU handling class                *
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
 * @file GFitsHDU.hpp
 * @brief GFitsHDU class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSHDU_HPP
#define GFITSHDU_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsHeader.hpp"
#include "GFitsData.hpp"
#include "GFitsImage.hpp"
#include "GFitsAsciiTable.hpp"
#include "GFitsBinTable.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsHDU
 *
 * @brief Implements the FITS Header Data Unit (HDU) interface
 *
 * The HDU is the basic unit of a FITS file. Each HDU consists of a header
 * and a data area. The header is composed of cards and is implemented by
 * the GFitsHeader class. The data are is either an image or a table and
 * is implemented by the abstract GFitsData base class.
 ***************************************************************************/
class GFitsHDU {

    // Friend classes
    friend class GFits;

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsHDU& hdu);

public:
    // Constructors and destructors
    GFitsHDU();
    GFitsHDU(const GFitsImage& image);
    GFitsHDU(const GFitsAsciiTable& table);
    GFitsHDU(const GFitsBinTable& table);
    GFitsHDU(const GFitsHDU& hdu);
    ~GFitsHDU();

    // Operators
    GFitsHDU& operator= (const GFitsHDU& hdu);

    // Methods
    std::string    extname(void) const;
    void           extname(const std::string& extname);
    int            extno(void) const;
    int            exttype(void) const;
    GFitsHeader*   header(void) const;
    GFitsData*     data(void) const;
    GFitsTableCol* column(const std::string& colname) const;
    GFitsTableCol* column(const int& colnum) const;
    void           primary(void);

private:
    // Private methods
    void        init_members(void);
    void        copy_members(const GFitsHDU& hdu);
    void        free_members(void);
    void        connect(__fitsfile* fptr);
    GFitsImage* new_image(void);
    void        open(__fitsfile* fptr, int hdunum);
    void        save(void);

    // Private data area
    __fitsfile   m_fitsfile;    //!< FITS file pointer
    int          m_hdunum;      //!< HDU number (starting from 1)
    std::string  m_name;        //!< HDU name (extname)
    int          m_type;        //!< HDU type
    GFitsHeader* m_header;      //!< HDU header
    GFitsData*   m_data;        //!< HDU data
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline std::string  GFitsHDU::extname(void) const { return m_name; }
inline int          GFitsHDU::extno(void) const { return m_hdunum; }
inline int          GFitsHDU::exttype(void) const { return m_type; }
inline GFitsHeader* GFitsHDU::header(void) const { return m_header; }
inline GFitsData*   GFitsHDU::data(void) const { return m_data; }

#endif /* GFITSHDU_HPP */
