/***************************************************************************
 *                    GNodeArray.i - Array of nodes class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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


/***********************************************************************//**
 * @class GNodeArray
 *
 * @brief Node array class
 ***************************************************************************/
class GNodeArray : public GContainer {
public:
    // Constructors and destructors
    GNodeArray(void);
    explicit GNodeArray(const int& num, const double* array);
    explicit GNodeArray(const GVector& vector);
    explicit GNodeArray(const std::vector<double>& vector);
    GNodeArray(const GNodeArray& array);
    virtual ~GNodeArray(void);

    // Methods
    void          clear(void);
    GNodeArray*   clone(void) const;
    double&       at(const int& index);
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
    void          load(const std::string& filename,
                       const std::string& extname = "NODES");
    void          save(const std::string& filename, const bool& clobber,
                       const std::string& extname = "NODES") const;
    void          read(const GFitsTable& table);
    void          write(GFits& file,
                        const std::string& extname = "NODES") const;
};


/***********************************************************************//**
 * @brief GNodeArray class extension
 ***************************************************************************/
%extend GNodeArray {
    double __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const double& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GNodeArray copy() {
        return (*self);
    }
};
