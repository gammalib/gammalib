/***************************************************************************
 *               GLATLtCube.i - Fermi/LAT lifetime cube class              *
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
 * @file GLATLtCube.i
 * @brief Fermi/LAT lifetime cube interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATLtCube.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GLATLtCube
 *
 * @brief Interface for the Fermi LAT lifetime cube.
 ***************************************************************************/
class GLATLtCube : public GBase {

public:
    // Constructors and destructors
    GLATLtCube(void);
    explicit GLATLtCube(const std::string& filename);
    GLATLtCube(const GLATLtCube& cube);
    virtual ~GLATLtCube(void);

    // Operators
    double operator()(const GSkyDir& dir, const GEnergy& energy,
                      _ltcube_ctheta fct) const;
    double operator()(const GSkyDir& dir, const GEnergy& energy,
                      _ltcube_ctheta_phi fct) const;
    double operator()(const GSkyDir& dir, const GEnergy& energy,
                      const GLATAeff& aeff) const;
    double operator()(const GSkyDir& dir, const GEnergy& energy,
                      const double& offset, const GLATPsf& psf,
                      const GLATAeff& aeff) const;

    // Methods
    void        clear(void);
    GLATLtCube* clone(void) const;
    void        load(const std::string& filename);
    void        save(const std::string& filename, bool clobber=false) const;
};


/***********************************************************************//**
 * @brief GLATLtCube class extension
 ***************************************************************************/
%extend GLATLtCube {
    GLATLtCube copy() {
        return (*self);
    }
};
