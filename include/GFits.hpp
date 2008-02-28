/***************************************************************************
 *                    GFits.hpp  - FITS file access class                  *
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

#ifndef GFITS_HPP
#define GFITS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include "GFitsCfitsio.hpp"
#include "GFitsHDU.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                            GFits class definition                       *
 ***************************************************************************/
class GFits {

public:
    // Constructors and destructors
    GFits();
    GFits(const GFits& fits);
    ~GFits();

    // Operators
    GFits& operator= (const GFits& fits);

    // Methods
    void      open(const std::string filename);
    void      save(void);
    void      saveto(const std::string filename, int clobber);
    void      close(void);
    GFitsHDU* hdu(std::string extname);
    GFitsHDU* hdu(int extno);
    
private:
    // Private methods
    void init_members(void);
    void copy_members(const GFits& fits);
    void free_members(void);

    // Private data area
    std::string  m_filename;    // FITS file name
    __fitsfile*  m_fitsfile;    // FITS file pointer
    int          m_num_hdu;     // Number of HDUs in file
    GFitsHDU*    m_hdu;         // Pointers to HDUs
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline 
GFitsHDU* GFits::hdu(int extno) 
{
    #if defined(G_RANGE_CHECK)
    if (extno < 1 || extno > m_num_hdu)
        throw GException::out_of_range("GFits::hdu(int)", extno, 1, m_num_hdu);
    #endif
    return &(m_hdu[extno+1]); 
}


#endif /* GFITS_HPP */
