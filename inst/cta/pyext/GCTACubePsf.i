/***************************************************************************
 *      GCTACubePsf.i - CTA cube analysis point spread function class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2018 by Chia-Chun Lu                                *
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
    explicit GCTACubePsf(const GFilename& filename);
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
                const GEnergies&     energies,
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
    void               fill(const GObservations& obs, GLog* log = NULL);
    const GSkyMap&     cube(void) const;
    const GEnergies&   energies(void) const;
    const GNodeArray&  deltas(void) const;
    double             delta_max(void) const;
    int                offset(const int& idelta, const int& iebin) const;
    void               read(const GFits& fits);
    void               write(GFits& file) const;
    void               load(const GFilename& filename);
    void               save(const GFilename& filename,
                            const bool&      clobber = false) const;
    const GFilename&   filename(void) const;
};


/***********************************************************************//**
 * @brief GCTACubePsf class extension
 ***************************************************************************/
%extend GCTACubePsf {
    GCTACubePsf copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        if self.filename().is_empty():
            fits = gammalib.GFits()
            self.write(fits)
            state = (self.filename(), fits)
        else:
            state = (self.filename(),)
        return state
    def __setstate__(self, state):
        if state[0].is_empty():
            self.__init__()
            self.read(state[1])
        else:
            self.__init__(state[0])
}
};
