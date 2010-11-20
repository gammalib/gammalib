/***************************************************************************
 *  GFitsImageUShort.i  - FITS unsigned short image class SWIG interface   *
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
 * @file GFitsImageUShort.i
 * @brief GFitsImageUShort class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsImageUShort.hpp"
%}

/***********************************************************************//**
 * @class GFitsImageUShort
 *
 * @brief SWIG interface for the FITS unsigned short image class.
 ***************************************************************************/
class GFitsImageUShort : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageUShort(void);
    explicit GFitsImageUShort(int nx, const unsigned short* pixels = NULL);
    explicit GFitsImageUShort(int nx, int ny, const unsigned short* pixels = NULL);
    explicit GFitsImageUShort(int nx, int ny, int nz, const unsigned short* pixels = NULL);
    explicit GFitsImageUShort(int nx, int ny, int nz, int nt, const unsigned short* pixels = NULL);
    GFitsImageUShort(const GFitsImageUShort& image);
    virtual ~GFitsImageUShort(void);

    // Methods
    double            pixel(const int& ix) const;
    double            pixel(const int& ix, const int& iy) const;
    double            pixel(const int& ix, const int& iy, const int& iz) const;
    double            pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*             pixels(void);
    int               type(void) const;
    GFitsImageUShort* clone(void) const;
};


/***********************************************************************//**
 * @brief GFitsImageUShort class extension
 ***************************************************************************/
%extend GFitsImageUShort {
    unsigned short __getitem__(int GFitsImageInx[]) {
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
    void __setitem__(int GFitsImageInx[], unsigned short value) {
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
    GFitsImageUShort copy() {
        return (*self);
    }
};
