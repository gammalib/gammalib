/***************************************************************************
 *                    GNodeArray.i - Array of nodes class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2021 by Juergen Knoedlseder                         *
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
 * @file GNodeArray.i
 * @brief Node array class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GNodeArray.hpp"
%}
%include "std_vector.i"
namespace std {
   %template(DoubleVector) vector<double>;
}


/***********************************************************************//**
 * @class GNodeArray
 *
 * @brief Node array class
 ***************************************************************************/
class GNodeArray : public GContainer {

public:
    // Constructors and destructors
    GNodeArray(void);
    explicit GNodeArray(const GFilename& filename);
    explicit GNodeArray(const GVector& vector);
    explicit GNodeArray(const std::vector<double>& vector);
    GNodeArray(const int& num, const double* array);
    GNodeArray(const GNodeArray& array);
    virtual ~GNodeArray(void);

    // Methods
    void          clear(void);
    GNodeArray*   clone(void) const;
    std::string   classname(void) const;
    int           size(void) const;
    bool          is_empty(void) const;
    void          append(const double& node);
    void          insert(const int& index, const double& node);
    void          remove(const int& index);
    void          reserve(const int& num);
    void          extend(const GNodeArray& nodes);
    void          nodes(const int& num, const double* array);
    void          nodes(const GVector& vector);
    void          nodes(const std::vector<double>& vector);
    double        interpolate(const double& value,
                              const std::vector<double>& vector) const;
    void          set_value(const double& value) const;
    const int&    inx_left(void) const;
    const int&    inx_right(void) const;
    const double& wgt_left(void) const;
    const double& wgt_right(void) const;
    const double& wgt_grad_left(void) const;
    const double& wgt_grad_right(void) const;
    void          load(const GFilename& filename);
    void          save(const GFilename& filename,
                       const bool&      clobber = false) const;
    void          read(const GFitsTable& table);
    void          write(GFits& fits,
                        const std::string& extname = "NODES") const;
};


/***********************************************************************//**
 * @brief GNodeArray class extension
 ***************************************************************************/
%extend GNodeArray {
    double __getitem__(const int& index) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            return (*self)[self->size()+index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Node index",
                                           index, self->size());
        }
    }
    GNodeArray* __getitem__(PyObject *param) {
        if (PySlice_Check(param)) {
            Py_ssize_t start = 0;
            Py_ssize_t stop  = 0;
            Py_ssize_t step  = 0;
            Py_ssize_t len   = self->size();
            if (PythonSlice_GetIndices(param, len, &start, &stop, &step) == 0) {
                GNodeArray* nodes = new GNodeArray;
                if (step > 0) {
                    for (int i = (int)start; i < (int)stop; i += (int)step) {
                        nodes->append((*self)[i]);
                    }
                }
                else {
                    for (int i = (int)start; i > (int)stop; i += (int)step) {
                        nodes->append((*self)[i]);
                    }
                }
                return nodes;
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
    void __setitem__(const int& index, const double& val) {
        // Counting from start, e.g. [2]
        if (index >= 0 && index < self->size()) {
            (*self)[index] = val;
        }
        // Counting from end, e.g. [-1]
        else if (index < 0 && self->size()+index >= 0) {
            (*self)[self->size()+index] = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Node index",
                                           index, self->size());
        }
    }
    GNodeArray copy() {
        return (*self);
    }
};
