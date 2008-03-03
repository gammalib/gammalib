/***************************************************************************
 *         GFitsDblImage.hpp  - FITS double precision image class          *
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
 * @file GFitsDblImage.hpp
 * @brief GFitsDblImage class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSDBLIMAGE_HPP
#define GFITSDBLIMAGE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"
#include "GFitsImage.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GFitsDblImage
 *
 * @brief Implements a FITS double precision image
 *
 ***************************************************************************/
class GFitsDblImage : public GFitsImage {

public:
    // Constructors and destructors
    GFitsDblImage();
    GFitsDblImage(int naxis, const int* naxes);
    GFitsDblImage(int naxis, const int* naxes, const double* pixels);
    GFitsDblImage(const GFitsDblImage& image);
    ~GFitsDblImage();

    // Operators
    GFitsDblImage& operator= (const GFitsDblImage& image);
    double&        operator() (const int& ix);
    const double&  operator() (const int& ix) const;
    double&        operator() (const int& ix, const int& iy);
    const double&  operator() (const int& ix, const int& iy) const;
    double&        operator() (const int& ix, const int& iy, const int& iz);
    const double&  operator() (const int& ix, const int& iy, const int& iz) const;
    double&        operator() (const int& ix, const int& iy, const int& iz, const int& it);
    const double&  operator() (const int& ix, const int& iy, const int& iz, const int& it) const;


    // Methods
    void           open(__fitsfile* fptr);
    void           link(double* pixels);
    void           save(void);
    void           close(void);
    GFitsDblImage* clone(void) const;
    void           set_nullval(const double* value);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsDblImage& image);
    void free_members(void);
    void fetch_pixels(void);

    // Private data area
    int     m_linked;        // Pixels are linked (don't delete them!)
    double* m_pixels;        // Pixels
    double* m_nulval;        // NULL value
    int     m_anynul;        // Number of NULLs encountered
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GFITSDBLIMAGE_HPP */
