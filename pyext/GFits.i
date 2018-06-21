/***************************************************************************
 *                        GFits.i - FITS file class                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2018 by Juergen Knoedlseder                         *
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
 * @file GFits.i
 * @brief FITS file class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsImage.hpp"
#include "GFitsImageByte.hpp"
#include "GFitsImageDouble.hpp"
#include "GFitsImageFloat.hpp"
#include "GFitsImageLong.hpp"
#include "GFitsImageLongLong.hpp"
#include "GFitsImageSByte.hpp"
#include "GFitsImageShort.hpp"
#include "GFitsImageULong.hpp"
#include "GFitsImageUShort.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsAsciiTable.hpp"
%}


/***********************************************************************//**
 * @class GFits
 *
 * @brief FITS file class
 ***************************************************************************/
class GFits : public GContainer {
public:
    // Constructors and destructors
    GFits(void);
    GFits(const GFilename& filename, const bool& create = false);
    GFits(const GFits& fits);
    virtual ~GFits(void);

    // Methods
    void             clear(void);
    GFits*           clone(void) const;
    std::string      classname(void) const;
    GFitsImage*      image(const int& extno);
    GFitsImage*      image(const std::string& extname);
    GFitsTable*      table(const int& extno);
    GFitsTable*      table(const std::string& extname);
    int              size(void) const;
    bool             is_empty(void) const;
    GFitsHDU*        set(const int& extno, const GFitsHDU& hdu);
    GFitsHDU*        set(const std::string& extname, const GFitsHDU& hdu);
    GFitsHDU*        append(const GFitsHDU& hdu);
    GFitsHDU*        insert(const int& extno, const GFitsHDU& hdu);
    GFitsHDU*        insert(const std::string& extname, const GFitsHDU& hdu);
    void             remove(const int& extno);
    void             remove(const std::string& extname);
    void             reserve(const int& num);
    void             extend(const GFits& fits);
    bool             contains(const int& extno) const;
    bool             contains(const std::string& extname) const;
    const GFilename& filename(void) const;
    int              extno(const std::string& extname) const;
    void             open(const GFilename& filename,
                          const bool&      create = false);
    void             save(const bool& clobber = false);
    void             saveto(const GFilename& filename,
                            const bool&      clobber = false);
    void             close(void);
    void             publish(const int& extno,
                             const std::string& name = "") const;
    void             publish(const std::string& extname,
                             const std::string& name = "") const;
};


/***********************************************************************//**
 * @brief GFits class SWIG extension
 ***************************************************************************/
%extend GFits {
    GFitsHDU* __getitem__(const int& extno) {
        // Counting from start, e.g. [2]
        if (extno >= 0 && extno < self->size()) {
            return (*self)[extno];
        }
        // Counting from end, e.g. [-1]
        else if (extno < 0 && self->size()+extno >= 0) {
            return (*self)[self->size()+extno];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Extension number",
                                           extno, self->size());
        }
    }
    GFitsHDU* __getitem__(const std::string& extname) {
        return (*self)[extname];
    }
    GFits* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GFits* fits = new GFits;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        fits->append(*(*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        fits->append(*(*self)[i]);
                    }
                }
                return fits;
            }
            else {
                throw GException::invalid_argument("__getitem__(PyObject)",
                                                   "Invalid slice indices");
            }
        }
        else {
            throw GException::invalid_argument("__getitem__(PyObject)","");
        }
    }
    void __setitem__(const int& extno, const GFitsHDU& hdu) {
        // Counting from start, e.g. [2]
        if (extno >= 0 && extno < self->size()) {
            self->set(extno, hdu);
        }
        // Counting from end, e.g. [-1]
        else if (extno < 0 && self->size()+extno >= 0) {
            self->set(self->size()+extno, hdu);
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Extension number",
                                           extno, self->size());
        }
    }
    void __setitem__(const std::string& extname, const GFitsHDU& hdu) {
        self->set(extname, hdu);
        return;
    }
    GFits copy() {
        return (*self);
    }
}
