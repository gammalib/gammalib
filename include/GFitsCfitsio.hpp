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
#define __ffclos(A, B) ffclos(A, B)
#define __ffthdu(A, B, C) ffthdu(A, B, C)
#define __ffmahd(A, B, C, D) ffmahd(A, B, C, D)
#define __ffghdt(A, B, C) ffghdt(A, B, C)
#define __ffghsp(A, B, C, D) ffghsp(A, B, C, D)
#define __ffgkey(A, B, C, D, E) ffgkey(A, B, C, D, E)
#define __ffgkyn(A, B, C, D, E, F) ffgkyn(A, B, C, D, E, F)
#define __ffuky(A, B, C, D, E, F) ffuky(A, B, C, D, E, F)
#define __TSTRING TSTRING

/* __ Type definition ____________________________________________________ */
typedef fitsfile __fitsfile;

/***************************************************************************
 *                          CFITSIO is not available                       *
 ***************************************************************************/
#else

/* __ Macros _____________________________________________________________ */
#define __ffgerr(A, B) __error(A, B)
#define __ffopen(A, B, C, D) __dummy()
#define __ffclos(A, B, C, D) __dummy()
#define __ffthdu(A, B, C) __dummy()
#define __ffmahd(A, B, C, D) __dummy()
#define __ffghdt(A, B, C) __dummy()
#define __ffghsp(A, B, C, D) __dummy()
#define __ffgkey(A, B, C, D, E) __dummy()
#define __ffgkyn(A, B, C, D, E, F) __dummy()
#define __ffuky(A, B, C, D, E, F) __dummy()
#define __TSTRING 0

/* __ Type definition ____________________________________________________ */
typedef int __fitsfile;

/* __ Dummy function _____________________________________________________ */
void __error(int status, char* err_text) {
    strcpy(err_text, "CFITSIO not available");
}
int __dummy(void) {
    return 0;
}

#endif

#endif /* GFITSCFITSIO_HPP */
