/***************************************************************************
 *               GCTACubeBackground.i - CTA cube background class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Michael Mayer                                    *
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
 * @file GCTACubeBackground.i
 * @brief CTA cube background class definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTACubeBackground.hpp"
%}


/***********************************************************************//**
 * @class GCTACubeBackground
 *
 * @brief CTA cube background class
 ***************************************************************************/
class GCTACubeBackground : public GBase {

public:
    // Constructors and destructors
    GCTACubeBackground(void);
    explicit GCTACubeBackground(const std::string& filename);
    GCTACubeBackground(const GCTACubeBackground& bgd);
    virtual ~GCTACubeBackground(void);
    
    // Operators
    GCTACubeBackground& operator=(const GCTACubeBackground& bgd);
    double            operator()(const GSkyDir& dir, const GEnergy& energy) const;
    
    
    // Implemented pure virtual methods
    void                       clear(void);
    GCTACubeBackground*          clone(void) const;
    std::string                classname(void) const;
    void                       set(const GSkymap& cube, const GEbounds& ebounds);
    void                       load(const std::string& filename);
    void                       read(const GFits& fits);
    void                       write(GFits& fits) const;
    const std::string&                filename(void) const;
    GSkyDir               mc(const GEnergy& energy,
                             const GTime& time,
                             GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;
    std::string                print(const GChatter& chatter = NORMAL) const;
    
    // Methods
    const GSkymap&     cube(void) const;
    const GEbounds&    ebounds(void) const;
    const GNodeArray&  elogmeans(void) const;
    void                       save(const std::string& filename,
                                    const bool& clobber = false) const;
    
    void                       set_mc_cone(const GSkyDir& centre, const double&  radius);
    double                  integral(const double& logE) const;
};


/***********************************************************************//**
 * @brief GCTACubeBackground class extension
 ***************************************************************************/
%extend GCTACubeBackground {
    GCTACubeBackground copy() {
        return (*self);
    }
};
