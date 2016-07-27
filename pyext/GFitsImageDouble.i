/***************************************************************************
 *         GFitsImageDouble.i - Double precision FITS image class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GFitsImageDouble.i
 * @brief Double precision FITS image class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsImageDouble.hpp"
#include "GException.hpp"
%}


/***********************************************************************//**
 * @class GFitsImageDouble
 *
 * @brief Double precision FITS image class
 ***************************************************************************/
class GFitsImageDouble : public GFitsImage {
public:
    // Constructors and destructors
    GFitsImageDouble(void);
    GFitsImageDouble(const int& nx, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const int& nz, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const int& nz, const int& nt, const double* pixels = NULL);
    GFitsImageDouble(const std::vector<int>& naxes, const double* pixels = NULL);
    GFitsImageDouble(const GFitsImageDouble& image);
    virtual ~GFitsImageDouble(void);

    // Methods
    void              clear(void);
    GFitsImageDouble* clone(void) const;
    std::string       classname(void) const;
    double            pixel(const int& ix) const;
    double            pixel(const int& ix, const int& iy) const;
    double            pixel(const int& ix, const int& iy, const int& iz) const;
    double            pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    void*             pixels(void);
    int               type(void) const;
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
