/***************************************************************************
 *           GCTACubeExposure.i - CTA cube analysis exposure class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Chia-Chun Lu                                *
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
 * @file GCTACubeExposure.i
 * @brief CTA cube analysis exposure class definition
 * @author Chia-Chun Lu
 */

%{
/* __ Includes ___________________________________________________________ */
#include "GCTACubeExposure.hpp"
%}


/***********************************************************************//**
 * @class GCTACubeExposure
 *
 * @brief CTA exposure cube class
 *
 * This class implements a CTA exposure cube which provides the average
 * exposure for binned analysis as function of sky position and energy.
 ***************************************************************************/
class GCTACubeExposure : public GBase {

public:   
    // Constructors and destructors
    GCTACubeExposure(void);
    GCTACubeExposure(const GCTACubeExposure& cube);
    explicit GCTACubeExposure(const GFilename& filename);
    explicit GCTACubeExposure(const GCTAEventCube& cube);
    GCTACubeExposure(const std::string&   wcs,
                     const std::string&   coords,
                     const double&        x,
                     const double&        y,
                     const double&        dx,
                     const double&        dy,
                     const int&           nx,
                     const int&           ny,
                     const GEnergies&     energies);
    virtual ~GCTACubeExposure(void);
    
    // Interpolation Operator
    double operator()(const GSkyDir& dir, const GEnergy& energy) const;

    // Methods
    void               clear(void);
    GCTACubeExposure*  clone(void) const;
    std::string        classname(void) const;
    void               set(const GCTAObservation& obs);
    void               fill(const GObservations& obs, GLog* log = NULL);
    const GSkyMap&     cube(void) const;
    const GEnergies&   energies(void) const;
    const GGti&        gti(void) const;
    const double&      livetime(void) const;
    const double&      ontime(void) const;
    double             deadc(void) const;
    void               read(const GFits& fits);
    void               write(GFits& file) const;
    void               load(const GFilename& filename);
    void               save(const GFilename& filename,
                            const bool&      clobber = false) const;
    const GFilename&   filename(void) const;
};

/***********************************************************************//**
 * @brief GCTACubeExposure class extension
 ***************************************************************************/
%extend GCTACubeExposure {
    GCTACubeExposure copy() {
        return (*self);
    }
};
