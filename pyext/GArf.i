/***************************************************************************
 *               GArf.i - XSPEC Auxiliary Response File class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2018 by Juergen Knoedlseder                         *
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
 * @file GArf.i
 * @brief XSPEC Auxiliary Response File class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GArf.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GArf
 *
 * @brief Auxiliary Response File class
 ***************************************************************************/
class GArf : public GBase {
public:
    // Constructors and destructors
    GArf(void);
    explicit GArf(const GFilename& filename);
    explicit GArf(const GEbounds& ebds);
    GArf(const GArf& arf);
    virtual ~GArf(void);

    // Operators
    GArf&   operator+=(const GArf& arf);
    GArf&   operator-=(const GArf& arf);
    GArf&   operator*=(const double& scale);
    double& operator()(const int& index, const int& col);
    double  operator()(const std::string& colname,
                       const GEnergy&     energy) const;

    // Methods
    void             clear(void);
    GArf*            clone(void) const;
    std::string      classname(void) const;
    int              size(void) const;
    int              columns(void) const;
    void             append(const std::string&         name,
                            const std::vector<double>& column);
    const GEbounds&  ebounds(void) const;
    void             load(const GFilename& filename);
    void             save(const GFilename& filename,
                          const bool&      clobber = false) const;
    void             read(const GFitsTable& table);
    void             write(GFits& fits) const;
    const GFilename& filename(void) const;
};


/***********************************************************************//**
 * @brief GArf class extension
 ***************************************************************************/
%extend GArf {
    double __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    std::vector<double> __getitem__(const std::string& colname) {
        return (*self)[colname];
    }
    void __setitem__(const int& index, const double& value) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = value;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    int __len__() {
        return (self->size());
    }
    GArf __add__(const GArf& arf) const {
        return ((*self) + arf);
    }
    GArf __sub__(const GArf& arf) const {
        return ((*self) - arf);
    }
    GArf __mul__(const double& scale) const {
        return ((*self) * scale);
    }
    // Python 2.x
    GArf __div__(const double& scale) const {
        return ((*self) / scale);
    }
    GArf __idiv__(const double& scale) {
        self->operator/=(scale);
        return (*self);
    }
    // Python 3.x
    GArf __truediv__(const double& scale) const {
        return ((*self) / scale);
    }
    GArf __itruediv__(const double& scale) {
        self->operator/=(scale);
        return (*self);
    }
    GArf copy() {
        return (*self);
    }
};
