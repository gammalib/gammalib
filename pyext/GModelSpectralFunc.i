/***************************************************************************
 *         GModelSpectralFunc.i  -  Spectral function model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GModelSpectralFunc.i
 * @brief Spectral function model class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralFunc.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralFunc
 *
 * @brief Spectral function model class
 ***************************************************************************/
class GModelSpectralFunc : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralFunc(void);
    explicit GModelSpectralFunc(const std::string& filename);
    explicit GModelSpectralFunc(const GXmlElement& xml);
    GModelSpectralFunc(const GModelSpectralFunc& model);
    virtual ~GModelSpectralFunc(void);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralFunc* clone(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng) const;
    virtual double              eval_gradients(const GEnergy& srcEng) const;
    virtual double              flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double              eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy             mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    double norm(void) const;
};


/***********************************************************************//**
 * @brief GModelSpectralFunc class extension
 ***************************************************************************/
%extend GModelSpectralFunc {
    GModelSpectralFunc copy() {
        return (*self);
    }
};
