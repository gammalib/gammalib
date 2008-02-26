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

// Public methods
public:
    // Constructors and destructors
    GFits();
    GFits(const GFits& fits);
    ~GFits();

    // Operators
    GFits& operator= (const GFits& fits);

    // Methods
    void open(std::string filename);
    void close(void);
    //GFitsHDU* hdu(std::string extname);
    //GFitsHDU* hdu(int extno);
    
// Methods and data that are available to derived classes
protected:
    // Protected methods

    // Protected data area
    std::string  m_filename;    // FITS file name
    __fitsfile*  m_fitsfile;    // FITS file pointer
    int          m_num_hdu;     // Number of HDUs in file
    GFitsHDU*    m_hdu;         // Pointers to HDUs

// Methods that are available to the base class only
private:
  void init_members(void);
  void copy_members(const GFits& fits);
  void free_members(void);
};

#endif /* GFITS_HPP */
