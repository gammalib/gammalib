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

#ifndef GFITSIMAGE_HPP
#define GFITSIMAGE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                        GFitsImage class definition                      *
 ***************************************************************************/
class GFitsImage : public GFitsData {

public:
    // Constructors and destructors
    GFitsImage();
    GFitsImage(const GFitsImage& image);
    ~GFitsImage();

    // Operators
    GFitsImage& operator= (const GFitsImage& image);

    // Methods
    void        open(__fitsfile*  fptr);
    void        close(void);
    GFitsImage* clone(void) const;
    
private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsImage& image);
    void free_members(void);

    // Private data area
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
