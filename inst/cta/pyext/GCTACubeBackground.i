/***************************************************************************
 *               GCTACubeBackground.i - CTA cube background class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2018 by Michael Mayer                               *
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
    explicit GCTACubeBackground(const GFilename& filename);
    explicit GCTACubeBackground(const GCTAEventCube& cube);
    GCTACubeBackground(const GCTACubeBackground& bgd);
    GCTACubeBackground(const std::string&   wcs,
                       const std::string&   coords,
                       const double&        x,
                       const double&        y,
                       const double&        dx,
                       const double&        dy,
                       const int&           nx,
                       const int&           ny,
                       const GEnergies&     energies);
    virtual ~GCTACubeBackground(void);

    // Operators
    double              operator()(const GCTAInstDir& dir,
                                   const GEnergy&     energy) const;

    // Methods
    void                clear(void);
    GCTACubeBackground* clone(void) const;
    std::string         classname(void) const;
    void                fill(const GObservations& obs, GLog* log = NULL);
    double              integral(const double& logE) const;
    void                read(const GFits& fits);
    void                write(GFits& file) const;
    void                load(const GFilename& filename);
    void                save(const GFilename& filename,
                             const bool&      clobber = false) const;
    const GSkyMap&      cube(void) const;
    const GEnergies&    energies(void) const;
    const GFilename&    filename(void) const;
};


/***********************************************************************//**
 * @brief GCTACubeBackground class extension
 ***************************************************************************/
%extend GCTACubeBackground {
    GCTACubeBackground copy() {
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
