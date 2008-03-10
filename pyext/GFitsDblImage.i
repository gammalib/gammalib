/***************************************************************************
 *  GFitsDblImage.hpp  - FITS double precision image class SWIG interface  *
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
 * @file GFitsDblImage.i
 * @brief GFitsDblImage class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsDblImage.hpp"
%}
%feature("notabstract") GFitsDblImage;

/***********************************************************************//**
 * @class GFitsDblImage
 *
 * @brief Implements FITS double precision image SWIG interface
 ***************************************************************************/
class GFitsDblImage : public GFitsImage {

public:
    // Constructors and destructors
    GFitsDblImage();
    GFitsDblImage(const GFitsDblImage& image);
    ~GFitsDblImage();

    // Methods
    void link(double* pixels);
    void set_nullval(const double* value);
};


/***********************************************************************//**
 * @brief GFitsDblImage class extension
 ***************************************************************************/
%extend GFitsDblImage {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
	    std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
	    str_buffer[100000] = '\0';
	    return str_buffer;
    }
    GFitsDblImage(int nx) {
        int naxis   = 1;
        int naxes[] = {nx};
        GFitsDblImage* image = new GFitsDblImage(naxis, naxes);
        return image;
    }
    GFitsDblImage(int nx, int ny) {
        int naxis   = 2;
        int naxes[] = {nx, ny};
        GFitsDblImage* image = new GFitsDblImage(naxis, naxes);
        return image;
    }
    GFitsDblImage(int nx, int ny, int nz) {
        int naxis   = 3;
        int naxes[] = {nx, ny, nz};
        GFitsDblImage* image = new GFitsDblImage(naxis, naxes);
        return image;
    }
    GFitsDblImage(int nx, int ny, int nz, int nt) {
        int naxis   = 4;
        int naxes[] = {nx, ny, nz, nt};
        GFitsDblImage* image = new GFitsDblImage(naxis, naxes);
        return image;
    }
    double get(const int& ix) {
        return (*self)(ix);
    }
    void set(const int& ix, const double& value) {
        (*self)(ix) = value;
    }
    double get(const int& ix, const int& iy) {
        return (*self)(ix, iy);
    }
    void set(const int& ix, const int& iy, const double& value) {
        (*self)(ix, iy) = value;
    }
    double get(const int& ix, const int& iy, const int& iz) {
        return (*self)(ix, iy, iz);
    }
    void set(const int& ix, const int& iy, const int& iz, const double& value) {
        (*self)(ix, iy, iz) = value;
    }
    double get(const int& ix, const int& iy, const int& iz, const int& it) {
        return (*self)(ix, iy, iz, it);
    }
    void set(const int& ix, const int& iy, const int& iz, const int& it,
             const double& value) {
        (*self)(ix, iy, iz, it) = value;
    }
    GFitsDblImage copy() {
        return (*self);
    }
};
