/***************************************************************************
 *    GFitsImageLongLong.i  - FITS long long image class SWIG interface    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GFitsImageLongLong.i
 * @brief GFitsImageLongLong class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsImageLongLong.hpp"
#include "GException.hpp"
%}


/***********************************************************************//**
 * @class GFitsImageLongLong
 *
 * @brief SWIG interface for the FITS long long image class.
 ***************************************************************************/
class GFitsImageLongLong : public GFitsImage {

public:
    // Constructors and destructors
    GFitsImageLongLong(void);
    explicit GFitsImageLongLong(int nx, const long long* pixels = NULL);
    explicit GFitsImageLongLong(int nx, int ny, const long long* pixels = NULL);
    explicit GFitsImageLongLong(int nx, int ny, int nz, const long long* pixels = NULL);
    explicit GFitsImageLongLong(int nx, int ny, int nz, int nt, const long long* pixels = NULL);
    GFitsImageLongLong(const GFitsImageLongLong& image);
    virtual ~GFitsImageLongLong(void);

    // Methods
    double              pixel(const int& ix) const;
    double              pixel(const int& ix, const int& iy) const;
    double              pixel(const int& ix, const int& iy, const int& iz) const;
    double              pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*               pixels(void);
    int                 type(void) const;
    GFitsImageLongLong* clone(void) const;
};


/***********************************************************************//**
 * @brief GFitsImageLongLong class extension
 ***************************************************************************/
%extend GFitsImageLongLong {
    long long __getitem__(int GFitsImageInx[]) {
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
    void __setitem__(int GFitsImageInx[], long long value) {
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
    GFitsImageLongLong copy() {
        return (*self);
    }
};
