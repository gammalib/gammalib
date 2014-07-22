/***************************************************************************
 *             GCTABackground.i - CTA background model base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTABackground.i
 * @brief CTA background model base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTABackground.hpp"
%}


/***********************************************************************//**
 * @class GCTABackground
 *
 * @brief Abstract base class for the CTA background model
 ***************************************************************************/
class GCTABackground : public GBase {
public:
    // Constructors and destructors
    GCTABackground(void);
    GCTABackground(const GCTABackground& bgd);
    virtual ~GCTABackground(void);

    // Pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety,
                              const bool&   etrue = false) const = 0;

    // Pure virtual methods
    virtual void                       clear(void) = 0;
    virtual GCTABackground*            clone(void) const = 0;
    virtual void                       load(const std::string& filename) = 0;
    virtual std::string                filename(void) const = 0;
    virtual GCTAInstDir                mc(const GEnergy& energy,
                                          const GTime& time,
                                          GRan& ran) const = 0;
    virtual const GModelSpectralNodes& spectrum(void) const = 0;
};


/***********************************************************************//**
 * @brief GCTABackground class extension
 ***************************************************************************/
%extend GCTABackground {
};
