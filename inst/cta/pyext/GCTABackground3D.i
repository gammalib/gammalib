/***************************************************************************
 *               GCTABackground3D.i - CTA 3D background class              *
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
 * @file GCTABackground3D.i
 * @brief CTA background model base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTABackground3D.hpp"
%}


/***********************************************************************//**
 * @class GCTABackground3D
 *
 * @brief Abstract base class for the CTA background model
 ***************************************************************************/
class GCTABackground3D : public GCTABackground {
public:
    // Constructors and destructors
    GCTABackground3D(void);
    explicit GCTABackground3D(const std::string& filename);
    GCTABackground3D(const GCTABackground3D& bgd);
    virtual ~GCTABackground3D(void);

    // Implemented pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety,
                              const bool&   etrue = false) const;

    // Implemented pure virtual methods
    void                       clear(void);
    GCTABackground3D*          clone(void) const;
    void                       load(const std::string& filename);
    std::string                filename(void) const;
    GCTAInstDir                mc(const GEnergy& energy,
                                  const GTime& time,
                                  GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;

    // Methods
    void read(const GFits& file);
};


/***********************************************************************//**
 * @brief GCTABackground3D class extension
 ***************************************************************************/
%extend GCTABackground3D {
    GCTABackground3D copy() {
        return (*self);
    }
};
