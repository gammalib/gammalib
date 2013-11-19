/***************************************************************************
 *                GPha.i - XSPEC Pulse Height Analyzer class               *
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
 * @file GPha.i
 * @brief XSPEC Pulse Height Analyzer class
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPha.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GPha
 *
 * @brief Pulse Height Analyzer class
 ***************************************************************************/
class GPha : public GBase {
public:
    // Constructors and destructors
    GPha(void);
    explicit GPha(const std::string& filename);
    explicit GPha(const GEbounds& ebds);
    explicit GPha(const int& bins);
    GPha(const GPha& pha);
    virtual ~GPha(void);

    // Methods
    void               clear(void);
    GPha*              clone(void) const;
    int                size(void) const;
    double&            at(const int& index);
    const GEbounds&    ebounds(void) const;
    double             counts(void) const;
    const double&      underflow(void) const;
    const double&      overflow(void) const;
    const double&      outflow(void) const;
    void               fill(const GEnergy& energy, const double& value = 1.0);
    void               load(const std::string& filename);
    void               save(const std::string& filename, const bool& clobber = false) const;
    void               read(const GFitsTable* hdu);
    void               write(GFits& fits) const;
    const std::string& filename(void) const;
};


/***********************************************************************//**
 * @brief GPha class extension
 ***************************************************************************/
%extend GPha {
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
    GPha copy() {
        return (*self);
    }
};
