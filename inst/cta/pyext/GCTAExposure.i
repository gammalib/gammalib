/***************************************************************************
 *            GCTAMeanPsf.hpp - CTA mean point spread function class       *
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
 * @brief CTA point spread function base class definition
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
 * @brief Abstract base class for the CTA point spread function
 *
 * This class implements the abstract base class for the CTA point spread
 * function.
 ***************************************************************************/
class GCTAMeanPsf : public GBase {

public:
   
    // Constructors and destructors
    GCTAMeanPsf(void);
    GCTAMeanPsf(const GCTAMeanPsf& psf);
    GCTAMeanPsf(const GObservations& obs, 
		const double& x, const double& y, 
		const double& dx, const double& dy,
		const int& nx, const int& ny,
		const double& emin, const double& emax, const int& nebins,
		const double& min, const double& max, 
		const int& nbins);
    virtual ~GCTAMeanPsf(void);

    // Operators
    GCTAMeanPsf& operator=(const GCTAMeanPsf& psf);

    // Methods
    void         clear(void);
    GCTAMeanPsf* clone(void) const;
    void         load(const std::string& filename);
    void         write(GFits& file) const;
    void         save(const std::string& filename,
		      const bool& clobber) const;
    std::string  print(const GChatter& chatter = NORMAL) const;   
    const GSkymap& map(void) const;
    const GEbounds& ebounds(void) const;
    const GNodeArray& deltas(void) const;
};

/***********************************************************************//**
 * @brief GCTAMeanPsf class extension
 ***************************************************************************/
%extend GCTAMeanPsf {
};
