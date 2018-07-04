/***************************************************************************
 *               GLATLtCube.i - Fermi/LAT livetime cube class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @brief Fermi/LAT livetime cube interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATLtCube.hpp"
%}

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_lat_exposure     = "EXPOSURE";
    const std::string extname_lat_wgtexposure  = "WEIGHTED_EXPOSURE";
    const std::string extname_lat_cthetabounds = "CTHETABOUNDS";
}


/***********************************************************************//**
 * @class GLATLtCube
 *
 * @brief Interface for the Fermi LAT livetime cube.
 ***************************************************************************/
class GLATLtCube : public GBase {

public:
    // Constructors and destructors
    GLATLtCube(void);
    explicit GLATLtCube(const GFilename& filename);
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
    std::string classname(void) const;
    void        load(const GFilename& filename);
    void        save(const GFilename& filename,
                     const bool&      clobber = false) const;
    void        read(const GFits& file);
    void        write(GFits& file) const;
};


/***********************************************************************//**
 * @brief GLATLtCube class extension
 ***************************************************************************/
%extend GLATLtCube {
    GLATLtCube copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        fits = gammalib.GFits()
        self.write(fits)
        state = (fits,)
        return state
    def __setstate__(self, state):
        self.__init__()
        if state[0].size() > 3:
            self.read(state[0])
}
};
