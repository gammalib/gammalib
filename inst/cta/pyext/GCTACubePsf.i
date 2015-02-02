/***************************************************************************
 *      GCTACubePsf.i - CTA cube analysis point spread function class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Chia-Chun Lu                                *
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
 * @file GCTACubePsf.i
 * @brief CTA cube analysis point spread function class definition
 * @author Chia-Chun Lu
 */

%{
/* __ Includes ___________________________________________________________ */
#include "GCTACubePsf.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTACubePsf
 *
 * @brief CTA point spread function for cube analysis
 *
 * This class implements a mean CTA point spread function which provides the
 * average point spread function for binned analysis as function of sky
 * position, log10 energy and delta angle between true and measured photon
 * direction.
 ***************************************************************************/
class GCTACubePsf : public GBase {

public:
   
    // Constructors and destructors
    GCTACubePsf(void);
    GCTACubePsf(const GCTACubePsf& cube);
    explicit GCTACubePsf(const std::string& filename);
    GCTACubePsf(const GCTAEventCube& cube,
                const double& dmax, const int& ndbins);
    GCTACubePsf(const std::string&   wcs,
                const std::string&   coords,
                const double&        x,
                const double&        y,
                const double&        dx,
                const double&        dy,
                const int&           nx,
                const int&           ny,
                const GEbounds&      ebounds,
                const double&        dmax,
                const int&           ndbins);
    virtual ~GCTACubePsf(void);

    // Interpolation operator
    double operator()(const GSkyDir& dir, 
                      const double & delta,
				      const GEnergy& energy) const;

    // Methods
    void               clear(void);
    GCTACubePsf*       clone(void) const;
    std::string        classname(void) const;
    void               set(const GCTAObservation& obs);
    void               fill(const GObservations& obs);
    const GSkymap&     map(void) const;
    const GEbounds&    ebounds(void) const;
    const GNodeArray&  deltas(void) const;
    const GNodeArray&  elogmeans(void) const;
    double             delta_max(void) const;
    int                offset(const int& idelta, const int& iebin) const;
    void               read(const GFits& fits);
    void               write(GFits& file) const;
    void               load(const std::string& filename);
    void               save(const std::string& filename, const bool& clobber) const;
    const std::string& filename(void) const;
};


/***********************************************************************//**
 * @brief GCTACubePsf class extension
 ***************************************************************************/
%extend GCTACubePsf {
    GCTACubePsf copy() {
        return (*self);
    }
};
