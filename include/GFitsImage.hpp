/***************************************************************************
 *                    GFitsImage.hpp  - FITS image class                   *
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
 * @file GFitsImage.hpp
 * @brief GFitsImage class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGE_HPP
#define GFITSIMAGE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsImage
 *
 * @brief Implements a FITS image
 *
 ***************************************************************************/
class GFitsImage : public GFitsData {

    // I/O friends
    friend ostream& operator<< (ostream& os, const GFitsImage& image);

public:
    // Constructors and destructors
    GFitsImage();
    GFitsImage(const GFitsImage& image);
    ~GFitsImage();

    // Operators
    GFitsImage& operator= (const GFitsImage& image);

    // Methods
    void        open(__fitsfile* fptr);
    void        save(void);
    void        close(void);
    GFitsImage* clone(void) const;
    
private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsImage& image);
    void free_members(void);
    void connect(__fitsfile* fptr);

    // Private data area
    __fitsfile  m_fitsfile;
    int         m_bitpix;
    int         m_naxis;
    long*       m_naxes;
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline 
GFitsImage* GFitsImage::clone(void) const 
{
    return new GFitsImage(*this);
}

#endif /* GFITSIMAGE_HPP */
