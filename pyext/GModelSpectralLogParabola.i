/***************************************************************************
 *     GModelSpectralLogParabola.i - Log parabola spectral model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Michael Mayer                               *
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
 * @file GModelSpectralLogParabola.i
 * @brief Log parabola spectral model class definition
 * @author Michael Mayer 
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralLogParabola.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralLogParabola
 *
 * @brief LogParabola spectral model class
 ***************************************************************************/
class GModelSpectralLogParabola : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralLogParabola(void);
    GModelSpectralLogParabola(const std::string& type,
                              const std::string& prefactor,
                              const std::string& index,
                              const std::string& pivot,
                              const std::string& curvature);
    GModelSpectralLogParabola(const double&  prefactor,
                              const double&  index,
                              const GEnergy& pivot,
                              const double&  curvature);
    explicit GModelSpectralLogParabola(const GXmlElement& xml);
    GModelSpectralLogParabola(const GModelSpectralLogParabola& model);
    virtual ~GModelSpectralLogParabola(void);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GModelSpectralLogParabola* clone(void) const;
    virtual std::string                classname(void) const;
    virtual std::string                type(void) const;
    virtual double                     eval(const GEnergy& srcEng,
                                            const GTime&   srcTime = GTime(),
                                            const bool&    gradients = false) const;
    virtual double                     flux(const GEnergy& emin,
                                            const GEnergy& emax) const;
    virtual double                     eflux(const GEnergy& emin,
                                             const GEnergy& emax) const;
    virtual GEnergy                    mc(const GEnergy& emin,
                                          const GEnergy& emax,
                                          const GTime&   time,
                                          GRan&          ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;

    // Other methods
    double  prefactor(void) const;
    void    prefactor(const double& prefactor);
    double  index(void) const;
    void    index(const double& index);
    GEnergy pivot(void) const;
    void    pivot(const GEnergy& pivot);
    double  curvature(void) const;
    void    curvature(const double& curvature);
};


/***********************************************************************//**
 * @brief GModelSpectralLogParabola class extension
 ***************************************************************************/
%extend GModelSpectralLogParabola {
    GModelSpectralLogParabola copy() {
        return (*self);
    }
};
