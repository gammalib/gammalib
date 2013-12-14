/***************************************************************************
 *                   GLATMeanPsf.i - Fermi/LAT mean PSF class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GLATMeanPsf.hpp
 * @brief Fermi/LAT mean PSF class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATMeanPsf.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GLATMeanPsf
 *
 * @brief Fermi/LAT mean PSF class
 ***************************************************************************/
class GLATMeanPsf : public GBase {
public:
    // Constructors and destructors
    GLATMeanPsf(void);
    GLATMeanPsf(const GSkyDir& dir, const GLATObservation& obs);
    GLATMeanPsf(const GLATMeanPsf& cube);
    virtual ~GLATMeanPsf(void);

    // Operators
    double       operator()(const double& offset, const double& logE);

    // Methods
    void               clear(void);
    GLATMeanPsf*       clone(void) const;
    int                size(void) const;
    void               set(const GSkyDir& dir, const GLATObservation& obs);
    int                noffsets(void) const;
    int                nenergies(void) const;
    const double&      offset(const int& inx) const;
    const double&      energy(const int& inx) const;
    const GSkyDir&     dir(void) const;
    const std::string& name(void) const;
    void               name(const std::string& name);
    const double&      thetamax(void) const;
    void               thetamax(const double& value);
    double             psf(const double& offset, const double& logE);
    double             exposure(const double& logE);
};


/***********************************************************************//**
 * @brief GLATMeanPsf class extension
 ***************************************************************************/
%extend GLATMeanPsf {
    GLATMeanPsf copy() {
        return (*self);
    }
};
