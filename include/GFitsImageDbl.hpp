/***************************************************************************
 *         GFitsImageDbl.hpp  - FITS double precision image class          *
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
 * @file GFitsImageDbl.hpp
 * @brief GFitsImageDbl class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSIMAGEDBL_HPP
#define GFITSIMAGEDBL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"
#include "GFitsData.hpp"
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GFitsImageDbl
 *
 * @brief Implements a FITS double precision image
 *
 ***************************************************************************/
class GFitsImageDbl : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageDbl();
    GFitsImageDbl(int naxis, const int* naxes);
    GFitsImageDbl(int naxis, const int* naxes, const double* pixels);
    GFitsImageDbl(const GFitsImageDbl& image);
    ~GFitsImageDbl();

    // Operators
    GFitsImageDbl& operator= (const GFitsImageDbl& image);
    double&        operator() (const int& ix);
    const double&  operator() (const int& ix) const;
    double&        operator() (const int& ix, const int& iy);
    const double&  operator() (const int& ix, const int& iy) const;
    double&        operator() (const int& ix, const int& iy, const int& iz);
    const double&  operator() (const int& ix, const int& iy, const int& iz)
                               const;
    double&        operator() (const int& ix, const int& iy, const int& iz,
                               const int& it);
    const double&  operator() (const int& ix, const int& iy, const int& iz,
                               const int& it) const;

    // Methods
    void    link(double* pixels);
    void    set_nullval(const double* value);
    double* pixels(void);

private:
    // Private methods
    void           init_members(void);
    void           copy_members(const GFitsImageDbl& image);
    void           free_members(void);
    void           fetch_pixels(void);
    void           open(__fitsfile* fptr);
    void           save(void);
    void           close(void);
    GFitsImageDbl* clone(void) const;

    // Private data area
    int     m_linked;        // Pixels are linked (don't delete them!)
    double* m_pixels;        // Pixels
    double* m_nulval;        // NULL value
    int     m_anynul;        // Number of NULLs encountered
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GFITSIMAGEDBL_HPP */
