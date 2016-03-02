/***************************************************************************
 *      GCTACubeEdisp.i - CTA cube analysis energy dispersion class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Michael Mayer                               *
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
 * @file GCTACubeEdisp.i
 * @brief CTA cube analysis energy dispersion class definition
 * @author Michael Mayer
 */

%{
/* __ Includes ___________________________________________________________ */
#include "GCTACubeEdisp.hpp"
%}


/***********************************************************************//**
 * @class GCTACubeEdisp
 *
 * @brief CTA energy dispersion for cube analysis
 *
 * This class implements a mean CTA energy dispersion which provides the
 * average energy dispersion for binned analysis as function of sky
 * position, log10 energy and delta  migra (fraction of true and reconstructed energy)
 *
 ***************************************************************************/
class GCTACubeEdisp : public GBase {

public:
   
    // Constructors and destructors
    GCTACubeEdisp(void);
    GCTACubeEdisp(const GCTACubeEdisp& cube);
    explicit GCTACubeEdisp(const GFilename& filename);
    GCTACubeEdisp(const GCTAEventCube& cube,
                const double& mmax, const int& nmbins);
    GCTACubeEdisp(const std::string&   wcs,
                const std::string&   coords,
                const double&        x,
                const double&        y,
                const double&        dx,
                const double&        dy,
                const int&           nx,
                const int&           ny,
                const GEbounds&      ebounds,
                const double&        mmax,
                const int&           nmbins);
    virtual ~GCTACubeEdisp(void);

    // Interpolation operator
    double operator()(const GSkyDir& dir, 
                      const double & migra,
				      const GEnergy& energy) const;

    // Methods
    void               clear(void);
    GCTACubeEdisp*       clone(void) const;
    std::string        classname(void) const;
    void               set(const GCTAObservation& obs);
    void               fill(const GObservations& obs, GLog* log = NULL);
    const GSkyMap&     map(void) const;
    const GEbounds&    ebounds(void) const;
    const GNodeArray&  migras(void) const;
    const GNodeArray&  elogmeans(void) const;
    double             migra_max(void) const;
    int                offset(const int& imigra, const int& iebin) const;
    GEbounds     ebounds_src(const GSkyDir& dir,
                             const GEnergy obsEng) const;
    void               read(const GFits& fits);
    void               write(GFits& file) const;
    void               load(const GFilename& filename);
    void               save(const GFilename& filename,
                            const bool&      clobber = false) const;
    const GFilename&   filename(void) const;
};


/***********************************************************************//**
 * @brief GCTACubeEdisp class extension
 ***************************************************************************/
%extend GCTACubeEdisp {
    GCTACubeEdisp copy() {
        return (*self);
    }
};
