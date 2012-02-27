/***************************************************************************
 *          GModelSpectralNodes.i  -  Spectral nodes model class           *
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
 * @file GModelSpectralNodes.i
 * @brief Spectral nodes model class Python interface definition
 * @author J. Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralNodes.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralNodes
 *
 * @brief Spectral nodes model class
 ***************************************************************************/
class GModelSpectralNodes : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralNodes(void);
    explicit GModelSpectralNodes(const GXmlElement& xml);
    GModelSpectralNodes(const GModelSpectralNodes& model);
    virtual ~GModelSpectralNodes(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralNodes* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng) const;
    virtual double               eval_gradients(const GEnergy& srcEng) const;
    virtual double               flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual double               eflux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GModelSpectralNodes class extension
 ***************************************************************************/
%extend GModelSpectralNodes {
    GModelSpectralNodes copy() {
        return (*self);
    }
};
