/***************************************************************************
 *               GFitsCfitsio.hpp  - CFITSIO interface header              *
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

#ifndef GFITSCFITSIO_HPP
#define GFITSCFITSIO_HPP

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                           CFITSIO is available                          *
 ***************************************************************************/
#if defined(G_CFITSIO)
#include <fitsio.h>

/* __ Macros _____________________________________________________________ */
#define __ffgerr(A, B) ffgerr(A, B)
#define __ffopen(A, B, C, D) ffopen(A, B, C, D)
#define __ffinit(A, B, C) ffinit(A, B, C)
#define __ffclos(A, B) ffclos(A, B)
#define __ffdelt(A, B) ffdelt(A, B)
#define __ffthdu(A, B, C) ffthdu(A, B, C)
#define __ffmahd(A, B, C, D) ffmahd(A, B, C, D)
#define __ffghdt(A, B, C) ffghdt(A, B, C)
#define __ffghsp(A, B, C, D) ffghsp(A, B, C, D)
#define __ffgkey(A, B, C, D, E) ffgkey(A, B, C, D, E)
#define __ffgkyn(A, B, C, D, E, F) ffgkyn(A, B, C, D, E, F)
#define __ffuky(A, B, C, D, E, F) ffuky(A, B, C, D, E, F)
#define __ffgnrw(A, B, C) ffgnrw(A, B, C)
#define __ffgncl(A, B, C) ffgncl(A, B, C)
#define __ffgtcl(A, B, C, D, E, F) ffgtcl(A, B, C, D, E, F)
#define __ffgcv(A, B, C, D, E, F, G, H, I, J) ffgcv(A, B, C, D, E, F, G, H, I, J)
#define __ffgcvs(A, B, C, D, E, F, G, H, I) ffgcvs(A, B, C, D, E, F, G, H, I)
#define __TLOGICAL TLOGICAL
#define __TSTRING  TSTRING
#define __TSHORT   TSHORT
#define __TLONG    TLONG
#define __TFLOAT   TFLOAT
#define __TDOUBLE  TDOUBLE

/* __ Type definition ____________________________________________________ */
typedef fitsfile __fitsfile;

/***************************************************************************
 *                          CFITSIO is not available                       *
 ***************************************************************************/
#else

/* __ Macros _____________________________________________________________ */
#define __ffgerr(A, B) __error(A, B)
#define __ffopen(A, B, C, D) __dummy()
#define __ffinit(A, B, C) __dummy()
#define __ffclos(A, B) __dummy()
#define __ffdelt(A, B) __dummy()
#define __ffthdu(A, B, C) __dummy()
#define __ffmahd(A, B, C, D) __dummy()
#define __ffghdt(A, B, C) __dummy()
#define __ffghsp(A, B, C, D) __dummy()
#define __ffgkey(A, B, C, D, E) __dummy()
#define __ffgkyn(A, B, C, D, E, F) __dummy()
#define __ffuky(A, B, C, D, E, F) __dummy()
#define __ffgnrw(A, B, C) __dummy()
#define __ffgncl(A, B, C) __dummy()
#define __ffgtcl(A, B, C, D, E, F) __dummy()
#define __ffgcv(A, B, C, D, E, F, G, H, I, J) __dummy()
#define __ffgcvs(A, B, C, D, E, F, G, H, I) __dummy()
#define __TBIT          1
#define __TSHORT       11
#define __TLOGICAL     14
#define __TSTRING      16
#define __TSHORT       21
#define __TLONG        41
#define __TLONGLONG    81
#define __TFLOAT       42
#define __TDOUBLE      82
#define __TCOMPLEX     83
#define __TDBLCOMPLEX 163

/* __ Type definition ____________________________________________________ */
typedef struct {
    int  HDUposition;  // HDU position in file; 0 = first HDU
    int* Fptr;         // Pointer to FITS file structure
} __fitsfile;

/* __ Dummy function _____________________________________________________ */
inline
void __error(int status, char* err_text) {
    strcpy(err_text, "CFITSIO not available");
}
inline
int __dummy(void) {
    return 0;
}

#endif

#endif /* GFITSCFITSIO_HPP */
