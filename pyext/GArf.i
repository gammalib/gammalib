/***************************************************************************
 *               GArf.i - XSPEC Auxiliary Response File class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
    explicit GArf(const std::string& filename);
    explicit GArf(const GEbounds& ebds);
    GArf(const GArf& arf);
    virtual ~GArf(void);

    // Methods
    void            clear(void);
    GArf*           clone(void) const;
    int             size(void) const;
    double&         at(const int& index);
    const GEbounds& ebounds(void) const;
    void            load(const std::string& filename);
    void            save(const std::string& filename,
                         const bool& clobber = false) const;
    void            read(const GFitsTable& table);
    void            write(GFits& fits) const;
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
    GArf copy() {
        return (*self);
    }
};
