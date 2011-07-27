/***************************************************************************
 *        GModelSpectralPlaw2.i  -  Spectral power law model class         *
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
 * @file GModelSpectralPlaw2.i
 * @brief Flux normalized power law spectral model class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralPlaw2.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralPlaw2
 *
 * @brief Flux normalized power law spectral model class
 ***************************************************************************/
class GModelSpectralPlaw2 : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralPlaw2(void);
    explicit GModelSpectralPlaw2(const double& integral, const double& index);
    explicit GModelSpectralPlaw2(const GXmlElement& xml);
    GModelSpectralPlaw2(const GModelSpectralPlaw2& model);
    virtual ~GModelSpectralPlaw2(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralPlaw2* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng) const;
    virtual double               eval_gradients(const GEnergy& srcEng) const;
    virtual double               flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;

    // Other methods
    double integral(void) const;
    double index(void) const;
    double emin(void) const;
    double emax(void) const;
};


/***********************************************************************//**
 * @brief GModelSpectralPlaw2 class extension
 ***************************************************************************/
%extend GModelSpectralPlaw2 {
    GModelSpectralPlaw2 copy() {
        return (*self);
    }
};
