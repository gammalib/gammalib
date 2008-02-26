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

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                          GFitsHDU class definition                      *
 ***************************************************************************/
class GFitsHDU {

// Public methods
public:
    // Constructors and destructors
    GFitsHDU();
    GFitsHDU(const GFitsHDU& hdu);
    ~GFitsHDU();

    // Operators
    GFitsHDU& operator= (const GFitsHDU& hdu);

    // Methods
    void         open(__fitsfile*  fptr, int hdunum);
    void         close(void);
    GFitsHeader* header(void);
    GFitsData*   data(void);
    
private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsHDU& hdu);
    void free_members(void);
    void move2hdu(void);

    // Private data area
    __fitsfile*  m_fitsfile;    // FITS file pointer
    std::string  m_name;        // HDU name
    int          m_num;         // HDU number (starting from 1)
    int          m_type;        // HDU type
    GFitsHeader  m_header;      // HDU header
    GFitsData    m_data;        // HDU data
};

#endif /* GFITSHDU_HPP */
