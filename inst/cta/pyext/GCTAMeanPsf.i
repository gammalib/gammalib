/***************************************************************************
 *         GCTAMeanPsf.i - CTA mean point spread function cube class       *
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
 * @file GCTAMeanPsf.hpp
 * @brief CCTA point spread function cube class definition
 * @author Chia-Chun Lu
 */

%{
/* __ Includes ___________________________________________________________ */
#include "GCTAMeanPsf.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAMeanPsf
 *
 * @brief class for the CTA point spread function
 *
 * This class implements a mean CTA point spread function which provides the
 * average point spread function for binned analysis as function of sky
 * position, log10 energy and delta angle between true and measured photon
 * direction.
 ***************************************************************************/
class GCTAMeanPsf : public GBase {

public:
   
    // Constructors and destructors
    GCTAMeanPsf(void);
    GCTAMeanPsf(const GCTAMeanPsf& cube);
    explicit GCTAMeanPsf(const std::string& filename);
    GCTAMeanPsf(const std::string&   wcs,
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
    virtual ~GCTAMeanPsf(void);

    // Interpolation operator
    double operator()(const GSkyDir& dir, 
                      const double & delta,
				      const GEnergy& energy) const;

    // Methods
    void               clear(void);
    GCTAMeanPsf*       clone(void) const;
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
 * @brief GCTAMeanPsf class extension
 ***************************************************************************/
%extend GCTAMeanPsf {
};
