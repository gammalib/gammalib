/***************************************************************************
 *          GCOMEventCube.i  -  COMPTEL event bin container class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCOMEventCube.hpp
 * @brief COMPTEL event bin container class interface definition
 * @author J. Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMEventCube.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventCube
 *
 * @brief CTA event bin container class Python interface
 ***************************************************************************/
class GCOMEventCube : public GEventCube {
public:
    // Constructors and destructors
    GCOMEventCube(void);
    explicit GCOMEventCube(const std::string& filename);
    explicit GCOMEventCube(const GSkymap& map, const GEbounds& ebds,
                           const GGti& gti, const double& phimin,
                           const double& dphi);
    GCOMEventCube(const GCOMEventCube& cube);
    virtual ~GCOMEventCube(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCOMEventCube* clone(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(const int& axis) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename,
                                const bool& clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;

    // Other methods
    void                   map(const GSkymap& map, const double& phimin,
                               const double& dphi);
    const GSkymap&         map(void) const;
    int                    nchi(void) const;
    int                    npsi(void) const;
    int                    nphi(void) const;
    int                    npix(void) const;
};


/***********************************************************************//**
 * @brief GCOMEventCube class extension
 ***************************************************************************/
%extend GCOMEventCube {
    GCOMEventCube copy() {
        return (*self);
    }
    GCOMEventBin* __getitem__(int index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(int index, const GCOMEventBin& val) {
        if (index>=0 && index < self->size()) {
            *((*self)[index]) = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
};
