/***************************************************************************
 *             GLATPsf.i - Fermi LAT point spread function class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2018 by Juergen Knoedlseder                         *
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
 * @file GLATPsf.i
 * @brief Fermi LAT point spread function class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATPsf.hpp"
%}


/***********************************************************************//**
 * @class GLATPsf
 *
 * @brief Interface for the Fermi/LAT point spread function
 ***************************************************************************/
class GLATPsf : public GBase {

public:
    // Constructors and destructors
    GLATPsf(void);
    GLATPsf(const GFilename& filename, const std::string& evtype);
    GLATPsf(const GLATPsf& psf);
    virtual ~GLATPsf(void);

    // Operators
    double operator()(const double& offset, const double& logE,
                      const double& ctheta);

    // Methods
    void               clear(void);
    GLATPsf*           clone(void) const;
    std::string        classname(void) const;
    const std::string& evtype(void) const;
    void               load(const GFilename&   filename,
                            const std::string& evtype);
    void               save(const GFilename& filename,
                            const bool&      clobber = false);
    void               read(const GFits& file, const std::string& evtype);
    void               write(GFits& file) const;
    int                size(void) const;
    int                nenergies(void) const;
    int                ncostheta(void) const;
    double             costhetamin(void) const;
    void               costhetamin(const double& ctheta);
    bool               has_phi(void) const;
    int                version(void) const;
};


/***********************************************************************//**
 * @brief GLATPsf class extension
 ***************************************************************************/
%extend GLATPsf {
    GLATPsf copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        fits = gammalib.GFits()
        self.write(fits)
        state = (fits, self.evtype())
        return state
    def __setstate__(self, state):
        self.__init__()
        if state[0].size() > 0:
            self.read(state[0], state[1])
}
};
