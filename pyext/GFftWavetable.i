/***************************************************************************
 *   GFftWavetable.i - Lookup table class for Fast Fourier transformation  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GFftWavetable.i
 * @brief Lookup table class interface definition for Fast Fourier transformation
 * @author Juergen Knoedlseder
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFftWavetable.hpp"
%}


/***********************************************************************//**
 * @class GFftWavetable
 *
 * @brief Lookup table class for Fast Fourier Transformation
 ***************************************************************************/
class GFftWavetable : public GBase {

public:
    // Constructors and destructors
    GFftWavetable(void);
    explicit GFftWavetable(const int& size);
    GFftWavetable(const GFftWavetable& wavetable);
    virtual ~GFftWavetable(void);

    // Methods
    void           clear(void);
    GFftWavetable* clone(void) const;
    std::string    classname(void) const;
    int            size(void) const;
    int            factors(void) const;
    int            factor(const int& index) const;
};


/***********************************************************************//**
 * @brief GFftWavetable class extension
 ***************************************************************************/
%extend GFftWavetable {
    GFftWavetable copy() {
        return (*self);
    }
};
