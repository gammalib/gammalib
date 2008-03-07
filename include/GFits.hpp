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
/**
 * @file GFits.hpp
 * @brief GFits class definition.
 * @author J. Knodlseder
 */

#ifndef GFITS_HPP
#define GFITS_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsHDU.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GFits
 *
 * @brief Implements FITS file interface
 *
 * GFits is the basic FITS file interface. All FITS file handlings operate
 * via members of GFits. A FITS file is composed of Header Data Units (HDU)
 * which are implemented by the GFitsHDU class.
 ***************************************************************************/
class GFits {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFits& fits);

public:
    // Constructors and destructors
    GFits();
    GFits(const GFits& fits);
    ~GFits();

    // Operators
    GFits& operator= (const GFits& fits);

    // Methods
    void      open(const std::string& filename);
    void      append_hdu(const GFitsHDU& hdu);
    void      save(void);
    void      saveto(const std::string& filename, int clobber = 0);
    void      close(void);
    GFitsHDU* hdu(const std::string& extname);
    GFitsHDU* hdu(int extno);
    int       num_hdus(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFits& fits);
    void free_members(void);

    // Private data area
    std::string  m_filename;    //!< FITS file name
    __fitsfile*  m_fitsfile;    //!< FITS file pointer
    int          m_num_hdu;     //!< Number of HDUs in file
    GFitsHDU*    m_hdu;         //!< Pointers to HDUs
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GFITS_HPP */
