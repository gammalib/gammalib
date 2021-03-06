/***************************************************************************
 *               GCsv.i - Comma-separated values table class               *
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
 * @file GCsv.hpp
 * @brief Comma-separated values table class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCsv.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCsv
 *
 * @brief Comma-separated values table class definition
 ***************************************************************************/
class GCsv : public GBase {
public:
    // Constructors and destructors
    GCsv(void);
    GCsv(const int& nrows, const int& ncols);
    GCsv(const GFilename& filename, const std::string& sep = " ");
    GCsv(const GCsv& csv);
    virtual ~GCsv(void);
 
    // Methods
    void        clear(void);
    GCsv*       clone(void) const;
    std::string classname(void) const;
    int         size(void) const;
    const int&  ncols(void) const;
    const int&  nrows(void) const;
    const int&  precision(void) const;
    void        precision(const int& precision);
    void        append(const std::vector<std::string>& list);
    std::string string(const int& row, const int& col) const;
    double      real(const int& row, const int& col) const;
    int         integer(const int& row, const int& col) const;
    void        string(const int& row, const int& col, const std::string& value);
    void        real(const int& row, const int& col, const double& value);
    void        integer(const int& row, const int& col, const int& value);
    void        load(const GFilename& filename, const std::string& sep = " ");
    void        save(const GFilename& filename, const std::string& sep = " ",
                     const bool& clobber = false) const;
};


/***********************************************************************//**
 * @brief GCsv class extension
 ***************************************************************************/
%extend GCsv {
    std::string __getitem__(int GTuple2D[]) {
        return (*self)(GTuple2D[0], GTuple2D[1]);
    }
    void __setitem__(int GTuple2D[], std::string value) {
        (*self)(GTuple2D[0], GTuple2D[1]) = value;
    }
    int __len__() {
        return (self->size());
    }
    GCsv copy() {
        return (*self);
    }
};
