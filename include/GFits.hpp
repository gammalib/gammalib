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
#include <string>
#include <vector>
#include <iostream>
#include "GLog.hpp"
#include "GFitsHDU.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"


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
    friend GLog&         operator<< (GLog& log, const GFits& fits);

public:
    // Constructors and destructors
    GFits(void);
    GFits(const std::string& filename);
    GFits(const GFits& fits);
    virtual ~GFits(void);

    // Operators
    GFits& operator= (const GFits& fits);

    // Methods
    void        clear(void);
    int         size(void) const;
    void        open(const std::string& filename);
    void        save(bool clobber = false);
    void        saveto(const std::string& filename, bool clobber = false);
    void        close(void);
    void        append(GFitsHDU* hdu);
    GFitsHDU*   hdu(const std::string& extname) const;
    GFitsHDU*   hdu(int extno) const;
    GFitsImage* image(const std::string& extname) const;
    GFitsImage* image(int extno) const;
    GFitsTable* table(const std::string& extname) const;
    GFitsTable* table(int extno) const;
    std::string name(void) const { return m_filename; }
    std::string print(void) const;

    // Complex single precision type
    typedef struct {
        float re;
        float im;
    } cfloat;

    // Complex double precision type
    typedef struct {
        float re;
        float im;
    } cdouble;

private:
    // Private methods
    void        init_members(void);
    void        copy_members(const GFits& fits);
    void        free_members(void);
    GFitsImage* new_image(void);
    GFitsImage* new_primary(void);

    // Private data area
    std::vector<GFitsHDU*> m_hdu;        //!< Pointers to HDUs
    std::string            m_filename;   //!< FITS file name
    void*                  m_fitsfile;   //!< FITS file pointer
    bool                   m_readwrite;  //!< FITS file is readwrite (true/false)
    bool                   m_created;    //!< FITS file has been created (true/false)
};

#endif /* GFITS_HPP */
