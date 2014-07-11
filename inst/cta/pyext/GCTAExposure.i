/***************************************************************************
 *                 GCTAExposure.i - CTA exposure cube class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @file GCTAExposure.hpp
 * @brief CTA exposure cube class definition
 * @author Chia-Chun Lu
 */

%{
/* __ Includes ___________________________________________________________ */
#include "GCTAExposure.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAExposure
 *
 * @brief CTA exposure cube class
 *
 * This class implements a CTA exposure cube which provides the average
 * exposure for binned analysis as function of sky position and log10
 * energy.
 ***************************************************************************/
class GCTAExposure : public GBase {
public:   
    // Constructors and destructors
    GCTAExposure(void);
    GCTAExposure(const GCTAExposure& cube);
    GCTAExposure(const std::string&   wcs,
                 const std::string&   coords,
                 const double&        x,
                 const double&        y,
                 const double&        dx,
                 const double&        dy,
                 const int&           nx,
                 const int&           ny,
                 const GEbounds&      ebounds);
    virtual ~GCTAExposure(void);
    
    // Interpolation Operator
    double          operator()(const GSkyDir& dir, const GEnergy& energy) const;
    // Methods
    void            clear(void);
    GCTAExposure*   clone(void) const;
    void            set(const GCTAObservation& obs);
    void            fill(const GObservations& obs);
    const GSkymap&  cube(void) const;
    const GEbounds& ebounds(void) const;
    void            write(GFits& file) const;
    void            load(const std::string& filename);
    void            save(const std::string& filename,
                         const bool& clobber = false) const;
};

/***********************************************************************//**
 * @brief GCTAExposure class extension
 ***************************************************************************/
%extend GCTAExposure {
};
