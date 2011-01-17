/***************************************************************************
 *  GFitsImageDouble.i  - FITS double precision image class SWIG interface *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsImageDouble.i
 * @brief GFitsImageDouble class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsImageDouble.hpp"
%}


/***********************************************************************//**
 * @class GFitsImageDouble
 *
 * @brief SWIG interface for the FITS double precision image class.
 ***************************************************************************/
class GFitsImageDouble : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageDouble(void);
    explicit GFitsImageDouble(int nx, const double* pixels = NULL);
    explicit GFitsImageDouble(int nx, int ny, const double* pixels = NULL);
    explicit GFitsImageDouble(int nx, int ny, int nz, const double* pixels = NULL);
    explicit GFitsImageDouble(int nx, int ny, int nz, int nt, const double* pixels = NULL);
    GFitsImageDouble(const GFitsImageDouble& image);
    virtual ~GFitsImageDouble(void);

    // Methods
    double            pixel(const int& ix) const;
    double            pixel(const int& ix, const int& iy) const;
    double            pixel(const int& ix, const int& iy, const int& iz) const;
    double            pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*             pixels(void);
    int               type(void) const;
    GFitsImageDouble* clone(void) const;
};


/***********************************************************************//**
 * @brief GFitsImageDouble class extension
 ***************************************************************************/
%extend GFitsImageDouble {
    double __getitem__(int GFitsImageInx[]) {
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
    void __setitem__(int GFitsImageInx[], double value) {
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
    GFitsImageDouble copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GFitsImageDouble type casts
 ***************************************************************************/
%inline %{
    GFitsImageDouble* cast_GFitsImageDouble(GFitsImage* image) {
        if (image->type() != 82)
            throw GException::fits_invalid_type("cast_GFitsImageDouble(GFitsImage*)",
                                                "Expecting double precision image.");
        return dynamic_cast<GFitsImageDouble*>(image);
    }
%};
