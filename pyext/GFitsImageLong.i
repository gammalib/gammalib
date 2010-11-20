/***************************************************************************
 *    GFitsImageLong.i  - FITS long integer image class SWIG interface     *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2010 by Jurgen Knodlseder                              *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsImageLong.i
 * @brief GFitsImageLong class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsImageLong.hpp"
%}

/***********************************************************************//**
 * @class GFitsImageLong
 *
 * @brief SWIG interface for the FITS long integer image class.
 ***************************************************************************/
class GFitsImageLong : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageLong(void);
    explicit GFitsImageLong(int nx, const long* pixels = NULL);
    explicit GFitsImageLong(int nx, int ny, const long* pixels = NULL);
    explicit GFitsImageLong(int nx, int ny, int nz, const long* pixels = NULL);
    explicit GFitsImageLong(int nx, int ny, int nz, int nt, const long* pixels = NULL);
    GFitsImageLong(const GFitsImageLong& image);
    virtual ~GFitsImageLong(void);

    // Methods
    double          pixel(const int& ix) const;
    double          pixel(const int& ix, const int& iy) const;
    double          pixel(const int& ix, const int& iy, const int& iz) const;
    double          pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*           pixels(void);
    int             type(void) const;
    GFitsImageLong* clone(void) const;
};


/***********************************************************************//**
 * @brief GFitsImageLong class extension
 ***************************************************************************/
%extend GFitsImageLong {
    long __getitem__(int GFitsImageInx[]) {
        if (GFitsImageInx[0] == 1)
            return self->at(GFitsImageInx[1]);
        else if (GFitsImageInx[0] == 2)
            return self->at(GFitsImageInx[1], GFitsImageInx[2]);
        else if (GFitsImageInx[0] == 3)
            return self->at(GFitsImageInx[1], GFitsImageInx[2], GFitsImageInx[3]);
        else if (GFitsImageInx[0] == 4)
            return self->at(GFitsImageInx[1], GFitsImageInx[2], GFitsImageInx[3], GFitsImageInx[4]);
        else
            throw GException::fits_wrong_image_operator("__getitem__(int)",
                                                        self->naxis(), GFitsImageInx[0]);
    }
    void __setitem__(int GFitsImageInx[], long value) {
        if (GFitsImageInx[0] == 1)
            self->at(GFitsImageInx[1]) = value;
        else if (GFitsImageInx[0] == 2)
            self->at(GFitsImageInx[1], GFitsImageInx[2]) = value;
        else if (GFitsImageInx[0] == 3)
            self->at(GFitsImageInx[1], GFitsImageInx[2], GFitsImageInx[3]) = value;
        else if (GFitsImageInx[0] == 4)
            self->at(GFitsImageInx[1], GFitsImageInx[2], GFitsImageInx[3], GFitsImageInx[4]) = value;
        else
            throw GException::fits_wrong_image_operator("__setitem__(int)",
                                                        self->naxis(), GFitsImageInx[0]);
    }
    GFitsImageLong copy() {
        return (*self);
    }
};
