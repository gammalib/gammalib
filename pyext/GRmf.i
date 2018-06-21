/***************************************************************************
 *             GRmf.i - XSPEC Redistribution Matrix File class             *
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
 * @file GRmf.i
 * @brief XSPEC Redistribution Matrix File class
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GRmf.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GRmf
 *
 * @brief Redistribution Matrix File class
 ***************************************************************************/
class GRmf : public GBase {
public:
    // Constructors and destructors
    GRmf(void);
    explicit GRmf(const GFilename& filename);
    GRmf(const GEbounds& etrue, const GEbounds& emeasured);
    GRmf(const GRmf& rmf);
    virtual ~GRmf(void);

    // Operators
    GRmf& operator+=(const GRmf& rmf);
    GRmf& operator-=(const GRmf& rmf);
    GRmf& operator*=(const double& scale);

    // Methods
    void                 clear(void);
    GRmf*                clone(void) const;
    std::string          classname(void) const;
    int                  size(void) const;
    int                  ntrue(void) const;
    int                  nmeasured(void) const;
    const GEbounds&      etrue(void) const;
    const GEbounds&      emeasured(void) const;
    GEbounds             etrue(const GEnergy& emeasured) const;
    GEbounds             emeasured(const GEnergy& etrue) const;
    const GMatrixSparse& matrix(void) const;
    void                 load(const GFilename& filename);
    void                 save(const GFilename&   filename,
                              const bool&        clobber = false,
                              const std::string& unit = "keV") const;
    void                 read(const GFitsTable& table);
    void                 write(GFits& fits,
                               const std::string& unit = "keV") const;
    const GFilename&     filename(void) const;
};


/***********************************************************************//**
 * @brief GRmf class extension
 ***************************************************************************/
%extend GRmf {
    double __getitem__(int GTuple[2]) {
        return self->at(GTuple[0], GTuple[1]);
    }
    void __setitem__(int GTuple[2], const double& value) {
        self->at(GTuple[0], GTuple[1]) = value;
    }
    int __len__() {
        return (self->size());
    }
    GRmf __add__(const GRmf& rmf) const {
        return ((*self) + rmf);
    }
    GRmf __sub__(const GRmf& rmf) const {
        return ((*self) - rmf);
    }
    GRmf __mul__(const double& scale) const {
        return ((*self) * scale);
    }
    // Python 2.x
    GRmf __div__(const double& scale) const {
        return ((*self) / scale);
    }
    GRmf __idiv__(const double& scale) {
        self->operator/=(scale);
        return (*self);
    }
    // Python 3.x
    GRmf __truediv__(const double& scale) const {
        return ((*self) / scale);
    }
    GRmf __itruediv__(const double& scale) {
        self->operator/=(scale);
        return (*self);
    }
    GRmf copy() {
        return (*self);
    }
};
