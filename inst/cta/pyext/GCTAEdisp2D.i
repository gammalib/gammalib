/***************************************************************************
 *             GCTAEdisp2D.i - CTA 2D energy dispersion class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2018 by Florent Forest                              *
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
 * @file GCTAEdisp2D.i
 * @brief CTA 2D energy dispersion class definition
 * @author Florent Forest
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEdisp2D.hpp"
%}


/***********************************************************************//**
 * @class GCTAEdisp2D
 *
 * @brief CTA 2D energy dispersion class
 ***************************************************************************/
class GCTAEdisp2D : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdisp2D(void);
    explicit GCTAEdisp2D(const GFilename& filename);
    GCTAEdisp2D(const GCTAEdisp2D& edisp);
    virtual ~GCTAEdisp2D(void);

    // Operators
    double operator()(const GEnergy& ereco,
                      const GEnergy& etrue,
                      const double&  theta = 0.0,
                      const double&  phi = 0.0,
                      const double&  zenith = 0.0,
                      const double&  azimuth = 0.0) const;

    // Implemented methods
    void         clear(void);
    GCTAEdisp2D* clone(void) const;
    std::string  classname(void) const;
    void         load(const GFilename& filename);
    GFilename    filename(void) const;
    GEnergy      mc(GRan&         ran,
                    const double& logE,
                    const double& theta = 0.0,
                    const double& phi = 0.0,
                    const double& zenith = 0.0,
                    const double& azimuth = 0.0) const;
    GEbounds     ebounds_obs(const double& logEsrc,
                             const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    GEbounds     ebounds_src(const double& logEobs,
                             const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    double       prob_erecobin(const GEnergy& ereco_min,
                               const GEnergy& ereco_max,
                               const GEnergy& etrue,
                               const double&  theta) const;

    // Methods
    void                     fetch(void) const;
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFitsTable& table);
    void                     write(GFitsBinTable& table) const;
    void                     save(const GFilename& filename,
                                  const bool&      clobber = false) const;
};


/***********************************************************************//**
 * @brief GCTAEdisp2D class extension
 ***************************************************************************/
%extend GCTAEdisp2D {
    GCTAEdisp2D copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        if self.filename().is_empty():
            hdu = gammalib.GFitsBinTable()
            self.write(hdu)
            state = (self.filename(), hdu)
        else:
            state = (self.filename(),)
        return state
    def __setstate__(self, state):
        if state[0].is_empty():
            self.__init__()
            if state[1].nrows() > 0:
                self.read(state[1])
        else:
            self.__init__(state[0])
}
};
