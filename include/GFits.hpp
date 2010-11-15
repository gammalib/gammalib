/***************************************************************************
 *                    GFits.hpp  - FITS file access class                  *
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
 * @file GFits.hpp
 * @brief GFits class definition.
 * @author J. Knodlseder
 */

#ifndef GFITS_HPP
#define GFITS_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHDU.hpp"


/***********************************************************************//**
 * @class GFits
 *
 * @brief Implements FITS file interface
 *
 * GFits is the basic FITS file interface. All FITS file handlings operate
 * via members of GFits. A FITS file is composed of Header Data Units (HDU)
 * which are implemented by the GFitsHDU class. Each HDU is composed of a
 * header (GFitsHeader) and some data (GFitsData).
 ***************************************************************************/
class GFits {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GFits& fits);

public:
    // Constructors and destructors
    GFits(void);
    GFits(const std::string& filename);
    GFits(const GFits& fits);
    virtual ~GFits(void);

    // Operators
    GFits& operator= (const GFits& fits);

    // Methods
    void      clear(void);
    int       size(void) const;
    void      open(const std::string& filename);
    void      close(void);
    void      save(void);
    void      saveto(const std::string& filename, bool clobber = false);
    void      append(const GFitsHDU& hdu);
    GFitsHDU* hdu(const std::string& extname) const;
    GFitsHDU* hdu(int extno) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFits& fits);
    void free_members(void);

    // Private data area
    std::string  m_filename;    //!< FITS file name
    void*        m_fitsfile;    //!< FITS file pointer
    int          m_readwrite;   //!< FITS file is read/write (1=true, 0=false)
    int          m_num_hdu;     //!< Number of HDUs in file
    GFitsHDU*    m_hdu;         //!< Pointers to HDUs
};

#endif /* GFITS_HPP */
