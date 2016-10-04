/***************************************************************************
 *          GModelSpectralComposite.i - Spectral power law model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 Michael Mayer                                       *
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
 * @file GModelSpecComposite spectral model class interface definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralComposite.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralComposite
 *
 * @brief Composite spectral model class
 ***************************************************************************/
class GModelSpectralComposite : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralComposite(void);
    explicit GModelSpectralComposite(const GXmlElement& xml);
    GModelSpectralComposite(const GModelSpectralComposite& model);
    virtual ~GModelSpectralComposite(void);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpectralComposite* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GEnergy& srcEng,
                                     const GTime&   srcTime = GTime(),
                                     const bool&    gradients = false) const;
    virtual double              flux(const GEnergy& emin,
                                     const GEnergy& emax) const;
    virtual double              eflux(const GEnergy& emin,
                                      const GEnergy& emax) const;
    virtual GEnergy             mc(const GEnergy& emin,
                                   const GEnergy& emax,
                                   const GTime&   time,
                                   GRan&          ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    int  components(void) const;
};


/***********************************************************************//**
 * @brief GModelSpectralComposite class extension
 ***************************************************************************/
%extend GModelSpectralComposite {
    GModelSpectralComposite copy() {
        return (*self);
    }
};
