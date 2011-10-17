/***************************************************************************
 *             GMWLSpectrum.i  -  Multi-wavelength spectrum class          *
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
 * @file GMWLSpectrum.i
 * @brief Multi-wavelength spectrum class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMWLSpectrum.hpp"
%}


/***********************************************************************//**
 * @class GMWLSpectrum
 *
 * @brief Multi-wavelength spectrum class Python interface
 ***************************************************************************/
class GMWLSpectrum : public GEventCube {
public:
    // Constructors and destructors
    GMWLSpectrum(void);
    explicit GMWLSpectrum(const std::string& filename);
    GMWLSpectrum(const GMWLSpectrum& spec);
    virtual ~GMWLSpectrum(void);

    // Implemented pure virtual methods
    virtual void          clear(void);
    virtual GMWLSpectrum* clone(void) const;
    virtual int           size(void) const;
    virtual int           dim(void) const;
    virtual int           naxis(int axis) const;
    virtual void          load(const std::string& filename);
    virtual void          save(const std::string& filename, bool clobber = false) const;
    virtual void          read(const GFits& file);
    virtual void          write(GFits& file) const;
    virtual int           number(void) const;

    // Other methods
    void                  load(const std::string& filename, const std::string& extname);
    void                  load(const std::string& filename, int extno);
    void                  read(const GFits& file, const std::string& extname);
    void                  read(const GFits& file, int extno);
    std::string           telescope(void) const;
    std::string           instrument(void) const;
};


/***********************************************************************//**
 * @brief GMWLSpectrum class extension
 ***************************************************************************/
%extend GMWLSpectrum {
    GMWLSpectrum copy() {
        return (*self);
    }
    GMWLDatum* __getitem__(int index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
};


/***********************************************************************//**
 * @brief GMWLSpectrum type casts
 ***************************************************************************/
%inline %{
    GMWLSpectrum* cast_GMWLSpectrum(GEvents* events) {
        GMWLSpectrum* cube = dynamic_cast<GMWLSpectrum*>(events);
        if (cube == NULL)
            throw GException::bad_type("cast_GMWLSpectrum(GEvents*)",
                                       "GEvents not of type GMWLSpectrum");            
        return cube;
    }
%}
