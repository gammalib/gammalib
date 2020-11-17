/***************************************************************************
 *                         GSkyMap.i - Sky map class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
 * @file GSkyMap.i
 * @brief Sky map class SWIG file.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyMap.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GSkyMap
 *
 * @brief GSkyMap class interface definition
 ***************************************************************************/
class GSkyMap : public GBase {

public:
    // Constructors and destructors
    GSkyMap(void);
    explicit GSkyMap(const GFilename& filename);
    explicit GSkyMap(const GFitsHDU& hdu);
    GSkyMap(const std::string& coords,
            const int&         nside,
            const std::string& order,
            const int&         nmaps = 1);
    GSkyMap(const std::string& wcs,
            const std::string& coords,
            const double&      x,
            const double&      y,
            const double&      dx,
            const double&      dy,
            const int&         nx,
            const int&         ny,
            const int&         nmaps = 1);
    GSkyMap(const GSkyMap& map);
    virtual ~GSkyMap(void);

    // Operators
    GSkyMap&      operator+=(const GSkyMap& map);
    GSkyMap&      operator+=(const double& value);
    GSkyMap&      operator-=(const GSkyMap& map);
    GSkyMap&      operator-=(const double& value);
    GSkyMap&      operator*=(const GSkyMap& map);
    GSkyMap&      operator*=(const double& factor);
    GSkyMap       operator+(const GSkyMap& map) const;
    GSkyMap       operator-(const GSkyMap& map) const;
    GSkyMap       operator*(const GSkyMap& map) const;
    double        operator()(const int& index, const int& map = 0);
    double        operator()(const GSkyPixel& pixel, const int& map = 0);
    double        operator()(const GSkyDir& dir, const int& map = 0) const;

    // Methods
    void                    clear(void);
    GSkyMap*                clone(void) const;
    std::string             classname(void) const;
    bool                    is_empty(void) const;
    const int&              npix(void) const;
    const int&              nx(void) const;
    const int&              ny(void) const;
    const int&              nmaps(void) const;
    void                    nmaps(const int& nmaps);
    const std::vector<int>& shape(void) const;
    void                    shape(const int& s1);
    void                    shape(const int& s1, const int& s2);
    void                    shape(const int& s1, const int& s2, const int& s3);
    void                    shape(const std::vector<int>& shape);
    int                     ndim(void) const;
    GSkyPixel               inx2pix(const int& index) const;
    GSkyDir                 inx2dir(const int& index) const;
    GSkyDir                 pix2dir(const GSkyPixel& pixel) const;
    int                     pix2inx(const GSkyPixel& pixel) const;
    int                     dir2inx(const GSkyDir& dir) const;
    GSkyPixel               dir2pix(const GSkyDir& dir) const;
    GNdarray                counts(void);
    double                  flux(const int& index, const int& map = 0) const;
    double                  flux(const GSkyPixel& pixel, const int& map = 0) const;
    GNdarray                flux(void);
    double                  solidangle(const int& index) const;
    double                  solidangle(const GSkyPixel& pixel) const;
    bool                    contains(const GSkyDir& dir) const;
    bool                    contains(const GSkyPixel& pixel) const;
    bool                    overlaps(const GSkyRegion& region) const;
    void                    smooth(const std::string& kernel, const double& par);
    const GSkyProjection*   projection(void) const;
    void                    projection(const GSkyProjection& proj);
    const double*           pixels(void) const;
    GSkyMap                 extract(const int& map, const int& nmaps = 1) const;
    GSkyMap                 extract(const int& startx, const int& stopx,
                                    const int& starty, const int& stopy) const;
    GSkyMap                 extract(const GSkyRegions& inclusions) const;
    void                    stack_maps(void);
    GSkyRegionCircle        region_circle(void) const;
    void                    load(const GFilename& filename);
    void                    save(const GFilename& filename,
                                 const bool&      clobber = false) const;
    void                    read(const GFitsHDU& hdu);
    GFitsHDU*               write(GFits& file,
                                  const std::string& extname = "") const;
    void                    publish(const std::string& name = "") const;
};


/***********************************************************************//**
 * @brief GSkyMap class extension
 *
 * @todo Implement __getitem__ and __setitem__ methods for GSkyPixel and
 * GSkyDir
 ***************************************************************************/
%extend GSkyMap {
    double __getitem__(int GTuple1D2D[]) {
        if (GTuple1D2D[1] < 0 || GTuple1D2D[1] >= self->npix()) {
            throw GException::out_of_range("GSkyMap::__getitem__(int*)",
                                           "Sky map index",
                                           GTuple1D2D[1], self->npix());
        }
        if (GTuple1D2D[0] == 1) {
            return (*self)(GTuple1D2D[1]);
        }
        else {
            if (GTuple1D2D[2] >= 0 && GTuple1D2D[2] < self->nmaps()) {
                return (*self)(GTuple1D2D[1], GTuple1D2D[2]);
            }
            else {
                throw GException::out_of_range("GSkyMap::__getitem__(int*)",
                                               "Sky map number",
                                               GTuple1D2D[2], self->nmaps());
            }
        }
    }
    /*
    double __getitem__(const GSkyPixel& pixel) {
        return (*self)(pixel);
    }
    */
    void __setitem__(int GTuple1D2D[], double value) {
        if (GTuple1D2D[1] < 0 || GTuple1D2D[1] >= self->npix()) {
            throw GException::out_of_range("GSkyMap::__setitem__(int*,double&)",
                                           "Sky map index",
                                           GTuple1D2D[1], self->npix());
        }
        if (GTuple1D2D[0] == 1) {
            (*self)(GTuple1D2D[1]) = value;
        }
        else {
            if (GTuple1D2D[2] >= 0 && GTuple1D2D[2] < self->nmaps()) {
                (*self)(GTuple1D2D[1], GTuple1D2D[2]) = value;
            }
            else {
                throw GException::out_of_range("GSkyMap::__setitem__(int*,double&)",
                                               "Sky map number",
                                               GTuple1D2D[2], self->nmaps());
            }
        }
    }
    /*
    void __setitem__(const GSkyPixel& pixel, double value) {
        (*self)(pixel) = value;
    }
    */
    GSkyMap copy() {
        return (*self);
    }
    GSkyMap sqrt() {
        return sqrt(*self);
    }
    GSkyMap log() {
        return log(*self);
    }
    GSkyMap log10() {
        return log10(*self);
    }
    GSkyMap clip(const double& thresh) {
        return clip(*self,thresh);
    }
    GSkyMap abs() {
        return abs(*self);
    }
    GSkyMap sign() {
        return sign(*self);
    }
    // Python 2.x operator/=
    GSkyMap __div__(const GSkyMap& map) {
        return ((*self) / map);
    }
    // Python 3.x operator/=
    GSkyMap __truediv__(const GSkyMap& map) {
        return ((*self) / map);
    }
    // Python 2.x operator/=
    GSkyMap __idiv__(const GSkyMap& map) {
        self->operator/=(map);
        return (*self);
    }
    // Python 3.x operator/=
    GSkyMap __itruediv__(const GSkyMap& map) {
        self->operator/=(map);
        return (*self);
    }
    // Python 2.x operator/=
    GSkyMap __idiv__(const double& factor) {
        self->operator/=(factor);
        return (*self);
    }
    // Python 3.x operator/=
    GSkyMap __itruediv__(const double& factor) {
        self->operator/=(factor);
        return (*self);
    }
    // Add pixel access operator as Python array
    PyObject* array(const int& imap = 0) {
        if (imap < 0 || imap >= self->nmaps()) {
            throw GException::out_of_range("GSkyMap::array(int&)", "Map index",
                                           imap, self->nmaps());
        }
        PyObject* array  = PyList_New(self->ny());
        int       offset = imap * self->nx() * self->ny();
        for (int iy = 0; iy < self->ny(); ++iy, offset += self->nx()) {
            PyObject* row = PyList_New(self->nx());
            for (int ix = 0; ix < self->nx(); ++ix) {
                PyList_SetItem(row, ix, PyFloat_FromDouble((self->pixels())[ix+offset]));
            }
            PyList_SetItem(array, iy, row);
        }
        return array;
    }
// Pickeling
%pythoncode {
    def __getstate__(self):
        fits = gammalib.GFits()
        self.write(fits)
        state = (fits,)
        return state
    def __setstate__(self, state):
        self.__init__()
        for hdu in state[0]:
            if ((hdu.exttype() == 0 and
                 hdu.has_card('NAXIS') and
                 hdu.integer('NAXIS') >= 2) or
                (hdu.exttype() == 0 and
                 hdu.has_card('PIXTYPE') and
                 hdu.string('PIXTYPE') == 'HEALPIX')):
                self.read(hdu)
                break
}
};
