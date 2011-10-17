/***************************************************************************
 *               GEvents.i  -  Abstract event container class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GEvents.i
 * @brief Abstract event container class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEvents.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GEvents
 *
 * @brief Abstract event container class Python interface
 ***************************************************************************/
class GEvents {
public:
    // Constructors and destructors
    GEvents(void);
    GEvents(const GEvents& events);
    virtual ~GEvents(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEvents*    clone(void) const = 0;
    virtual int         size(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual void        save(const std::string& filename, bool clobber = false) const = 0;
    virtual void        read(const GFits& file) = 0;
    virtual void        write(GFits& file) const = 0;
    virtual int         number(void) const = 0;

    // Implemented methods
    void                ebounds(const GEbounds& ebounds);
    void                gti(const GGti& gti);
    GTime               tstart(void) const;
    GTime               tstop(void) const;
    GEnergy             emin(void) const;
    GEnergy             emax(void) const;
    const GEbounds&     ebounds(void) const;
    const GGti&         gti(void) const;
};


/***********************************************************************//**
 * @brief GEvents class extension
 ***************************************************************************/
%extend GEvents {
    char *__str__() {
        return tochar(self->print());
    }
    GEvent* __getitem__(int index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GEvent& val) {
        if (index>=0 && index < self->size())
            *((*self)[index]) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
};
