/***************************************************************************
 *               GFitsCfitsio.hpp  - CFITSIO interface header              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsCfitsio.hpp
 * @brief CFITSIO interface header
 * @author Juergen Knoedlseder
 */

#ifndef GFITSCFITSIO_HPP
#define GFITSCFITSIO_HPP

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstring>

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                           CFITSIO is available                          *
 ***************************************************************************/
#if defined(HAVE_LIBCFITSIO)

#if defined(HAVE_CFITSIO_FITSIO_H)
#include <cfitsio/fitsio.h>
#elif defined(HAVE_LIBCFITSIO0_FITSIO_H)
#include <libcfitsio0/fitsio.h>
#else
#include <fitsio.h>
#endif

/* __ Macros _____________________________________________________________ */
#define __ffclos(A, B) ffclos(A, B)
#define __ffcrim(A, B, C, D, E) ffcrim(A, B, C, D, E)
#define __ffcrtb(A, B, C, D, E, F, G, H, I) ffcrtb(A, B, C, D, E, F, G, H, I)
#define __ffdcol(A, B, C) ffdcol(A, B, C)
#define __ffdelt(A, B) ffdelt(A, B)
#define __ffdhdu(A, B, C) ffdhdu(A, B, C)
#define __ffdrow(A, B, C, D) ffdrow(A, B, C, D)
#define __ffgabc(A, B, C, D, E, F) ffgabc(A, B, C, D, E, F)
#define __ffgcv(A, B, C, D, E, F, G, H, I, J) ffgcv(A, B, C, D, E, F, G, H, I, J)
#define __ffgcvb(A, B, C, D, E, F, G, H, I) ffgcvb(A, B, C, D, E, F, G, H, I)
#define __ffgcvs(A, B, C, D, E, F, G, H, I) ffgcvs(A, B, C, D, E, F, G, H, I)
#define __ffgdes(A, B, C, D, E, F) ffgdes(A, B, C, D, E, F)
#define __ffgerr(A, B) ffgerr(A, B)
#define __ffghdt(A, B, C) ffghdt(A, B, C)
#define __ffghsp(A, B, C, D) ffghsp(A, B, C, D)
#define __ffgidm(A, B, C) ffgidm(A, B, C)
#define __ffgidt(A, B, C) ffgidt(A, B, C)
#define __ffgiet(A, B, C) ffgiet(A, B, C)
#define __ffgipr(A, B, C, D, E, F) ffgipr(A, B, C, D, E, F)
#define __ffgisz(A, B, C, D) ffgisz(A, B, C, D)
#define __ffgky(A, B, C, D, E, F) ffgky(A, B, C, D, E, F)
#define __ffgkey(A, B, C, D, E) ffgkey(A, B, C, D, E)
#define __ffgkyn(A, B, C, D, E, F) ffgkyn(A, B, C, D, E, F)
#define __ffgnrw(A, B, C) ffgnrw(A, B, C)
#define __ffgncl(A, B, C) ffgncl(A, B, C)
#define __ffgsv(A, B, C, D, E, F, G, H, I) ffgsv(A, B, C, D, E, F, G, H, I)
#define __ffgtcl(A, B, C, D, E, F) ffgtcl(A, B, C, D, E, F)
#define __ffibin(A, B, C, D, E, F, G, H, I) ffibin(A, B, C, D, E, F, G, H, I)
#define __fficol(A, B, C, D, E) fficol(A, B, C, D, E)
#define __ffiimg(A, B, C, D, E) ffiimg(A, B, C, D, E)
#define __ffiimgll(A, B, C, D, E) ffiimgll(A, B, C, D, E)
#define __ffinit(A, B, C) ffinit(A, B, C)
#define __ffirow(A, B, C, D) ffirow(A, B, C, D)
#define __ffitab(A, B, C, D, E, F, G, H, I, J) ffitab(A, B, C, D, E, F, G, H, I, J)
#define __ffopen(A, B, C, D) ffopen(A, B, C, D)
#define __ffmahd(A, B, C, D) ffmahd(A, B, C, D)
#define __ffpcn(A, B, C, D, E, F, G, H, I) ffpcn(A, B, C, D, E, F, G, H, I)
#define __ffpcom(A, B, C) ffpcom(A, B, C)
#define __ffphis(A, B, C) ffphis(A, B, C)
#define __ffpss(A, B, C, D, E, F) ffpss(A, B, C, D, E, F)
#define __ffprec(A, B, C) ffprec(A, B, C)
#define __ffsrow(A, B, C, D) ffsrow(A, B, C, D)
#define __ffthdu(A, B, C) ffthdu(A, B, C)
#define __ffuky(A, B, C, D, E, F) ffuky(A, B, C, D, E, F)
#define __ffukye(A, B, C, D, E, F) ffukye(A, B, C, D, E, F)
#define __ffukyd(A, B, C, D, E, F) ffukyd(A, B, C, D, E, F)
#define __ffukyj(A, B, C, D, E) ffukyj(A, B, C, D, E)
#define __ffukyl(A, B, C, D, E) ffukyl(A, B, C, D, E)
#define __ffukys(A, B, C, D, E) ffukys(A, B, C, D, E)
#define __ffukyu(A, B, C, D) ffukyu(A, B, C, D)
#define __TNULL       0
#define __TBIT        TBIT
#define __TBYTE       TBYTE
#define __TSBYTE      TSBYTE
#define __TLOGICAL    TLOGICAL
#define __TSTRING     TSTRING
#define __TUSHORT     TUSHORT
#define __TSHORT      TSHORT
#define __TUINT       TUINT
#define __TINT        TINT
#define __TULONG      TULONG
#define __TLONG       TLONG
#define __TFLOAT      TFLOAT
#define __TLONGLONG   TLONGLONG
#define __TDOUBLE     TDOUBLE
#define __TCOMPLEX    TCOMPLEX
#define __TDBLCOMPLEX TDBLCOMPLEX

/* __ Type definition ____________________________________________________ */
typedef fitsfile __fitsfile;

/***************************************************************************
 *                          CFITSIO is not available                       *
 ***************************************************************************/
#else

/* __ Macros _____________________________________________________________ */
#define __ffclos(A, B) __dummy()
#define __ffcrim(A, B, C, D, E) __dummy()
#define __ffcrtb(A, B, C, D, E, F, G, H, I) __dummy()
#define __ffdcol(A, B, C) __dummy()
#define __ffdelt(A, B) __dummy()
#define __ffdhdu(A, B, C) __dummy()
#define __ffdrow(A, B, C, D) __dummy()
#define __ffgabc(A, B, C, D, E, F) __dummy()
#define __ffgcv(A, B, C, D, E, F, G, H, I, J) __dummy()
#define __ffgcvb(A, B, C, D, E, F, G, H, I) __dummy()
#define __ffgcvs(A, B, C, D, E, F, G, H, I) __dummy()
#define __ffgdes(A, B, C, D, E, F) __dummy()
#define __ffgerr(A, B) __error(A, B)
#define __ffghdt(A, B, C) __dummy()
#define __ffghsp(A, B, C, D) __dummy()
#define __ffgidm(A, B, C) __dummy()
#define __ffgidt(A, B, C) __dummy()
#define __ffgiet(A, B, C) __dummy()
#define __ffgipr(A, B, C, D, E, F) __dummy()
#define __ffgisz(A, B, C, D) __dummy()
#define __ffgky(A, B, C, D, E, F) __dummy()
#define __ffgkey(A, B, C, D, E) __dummy()
#define __ffgkyn(A, B, C, D, E, F) __dummy()
#define __ffgnrw(A, B, C) __dummy()
#define __ffgncl(A, B, C) __dummy()
#define __ffgsv(A, B, C, D, E, F, G, H, I) __dummy()
#define __ffgtcl(A, B, C, D, E, F) __dummy()
#define __ffibin(A, B, C, D, E, F, G, H, I) __dummy()
#define __fficol(A, B, C, D, E) __dummy()
#define __ffiimg(A, B, C, D, E) __dummy()
#define __ffiimgll(A, B, C, D, E) __dummy()
#define __ffinit(A, B, C) __dummy()
#define __ffirow(A, B, C, D) __dummy()
#define __ffitab(A, B, C, D, E, F, G, H, I, J) __dummy()
#define __ffmahd(A, B, C, D) __dummy()
#define __ffpcn(A, B, C, D, E, F, G, H, I) __dummy()
#define __ffopen(A, B, C, D) __dummy()
#define __ffpcom(A, B, C) __dummy()
#define __ffphis(A, B, C) __dummy()
#define __ffpss(A, B, C, D, E, F) __dummy()
#define __ffprec(A, B, C) __dummy()
#define __ffsrow(A, B, C, D) __dummy()
#define __ffthdu(A, B, C) __dummy()
#define __ffuky(A, B, C, D, E, F) __dummy()
#define __ffukye(A, B, C, D, E, F) __dummy()
#define __ffukyd(A, B, C, D, E, F) __dummy()
#define __ffukyj(A, B, C, D, E) __dummy()
#define __ffukyl(A, B, C, D, E) __dummy()
#define __ffukys(A, B, C, D, E) __dummy()
#define __ffukyu(A, B, C, D) __dummy()
#define __TNULL         0
#define __TBIT          1
#define __TBYTE        11
#define __TSBYTE       12
#define __TLOGICAL     14
#define __TSTRING      16
#define __TUSHORT      20
#define __TSHORT       21
#define __TUINT        30
#define __TINT         31
#define __TULONG       40
#define __TLONG        41
#define __TFLOAT       42
#define __TLONGLONG    81
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
void __error(int status, char* err_text)
{
    std::strcpy(err_text, "CFITSIO not available");
}
inline
int __dummy(void)
{
    return 0;
}

#endif

/* __ Common macros ______________________________________________________ */
#define FPTR(A)         ((__fitsfile*)A)
#define FHANDLE(A)      ((__fitsfile**)&A)
#define FPTR_COPY(A, B) *((__fitsfile*)A) = *((__fitsfile*)B)

#endif /* GFITSCFITSIO_HPP */
