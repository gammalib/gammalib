/***************************************************************************
 *              GMWLSpectrum.i - Multi-wavelength spectrum class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @brief Multi-wavelength spectrum class interface definition
 * @author Juergen Knoedlseder
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
    virtual int           naxis(const int& axis) const;
    virtual void          load(const std::string& filename);
    virtual void          save(const std::string& filename,
                               const bool& clobber = false) const;
    virtual void          read(const GFits& file);
    virtual void          write(GFits& file) const;
    virtual int           number(void) const;

    // Other methods
    void               load(const std::string& filename, const std::string& extname);
    void               load(const std::string& filename, const int& extno);
    void               read(const GFits& file, const std::string& extname);
    void               read(const GFits& file, const int& extno);
    const std::string& telescope(void) const;
    const std::string& instrument(void) const;
};


/***********************************************************************//**
 * @brief GMWLSpectrum class extension
 ***************************************************************************/
%extend GMWLSpectrum {
    GMWLSpectrum copy() {
        return (*self);
    }
    GMWLDatum* __getitem__(int index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", "Spectral point",
                                           index, self->size());
        }
    }
};
