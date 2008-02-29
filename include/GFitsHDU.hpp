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

#ifndef GFITSHDU_HPP
#define GFITSHDU_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHeader.hpp"
#include "GFitsData.hpp"
#include "GFitsImage.hpp"
#include "GFitsAsciiTable.hpp"
#include "GFitsBinTable.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                          GFitsHDU class definition                      *
 ***************************************************************************/
class GFitsHDU {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsHDU& hdu);

public:
    // Constructors and destructors
    GFitsHDU();
    GFitsHDU(const GFitsHDU& hdu);
    ~GFitsHDU();

    // Operators
    GFitsHDU& operator= (const GFitsHDU& hdu);

    // Methods
    void           open(__fitsfile* fptr, int hdunum);
    void           save(void);
    std::string    extname(void) const;
    int            extno(void) const;
    int            exttype(void) const;
    GFitsHeader*   header(void) const;
    GFitsData*     data(void) const;
    GFitsTableCol* column(const std::string& colname) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsHDU& hdu);
    void free_members(void);

    // Private data area
    __fitsfile   m_fitsfile;    // FITS file pointer
    std::string  m_name;        // HDU name
    int          m_type;        // HDU type
    GFitsHeader* m_header;      // HDU header
    GFitsData*   m_data;        // HDU data
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline std::string  GFitsHDU::extname(void) const { return m_name; }
inline int          GFitsHDU::extno(void) const { return m_fitsfile.HDUposition; }
inline int          GFitsHDU::exttype(void) const { return m_type; }
inline GFitsHeader* GFitsHDU::header(void) const { return m_header; }
inline GFitsData*   GFitsHDU::data(void) const { return m_data; }

#endif /* GFITSHDU_HPP */
