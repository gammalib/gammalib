/***************************************************************************
 *                GEvents.i - Abstract event container class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @brief Abstract event container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEvents.hpp"
#include "GEvent.hpp"
#include "GEventBin.hpp"
#include "GEventAtom.hpp"
#include "GTools.hpp"
%}

/* __ Typemaps ___________________________________________________________ */
%typemap(out) GEvent* {
    if (dynamic_cast<GEventBin*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GEventBin, 0 |  0 );
    }
    else if (dynamic_cast<GEventAtom*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GEventAtom, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GEvent, 0 |  0 );
    }
}


/***********************************************************************//**
 * @class GEvents
 *
 * @brief Abstract event container class Python interface
 ***************************************************************************/
class GEvents : public GBase {
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
    virtual void        save(const std::string& filename,
                             const bool& clobber = false) const = 0;
    virtual void        read(const GFits& file) = 0;
    virtual void        write(GFits& file) const = 0;
    virtual int         number(void) const = 0;

    // Implemented methods
    void                ebounds(const GEbounds& ebounds);
    void                gti(const GGti& gti);
    const GEbounds&     ebounds(void) const;
    const GGti&         gti(void) const;
    GTime               tstart(void) const;
    GTime               tstop(void) const;
    GEnergy             emin(void) const;
    GEnergy             emax(void) const;
};


/***********************************************************************//**
 * @brief GEvents class extension
 ***************************************************************************/
%extend GEvents {
    GEvent* __getitem__(int index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", "Event index",
                                           index, self->size());
        }
    }
    void __setitem__(int index, const GEvent& event) {
        if (index>=0 && index < self->size()) {
            *((*self)[index]) = event;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", "Event index",
                                           index, self->size());
        }
    }
};
