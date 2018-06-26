/***************************************************************************
 *       GCTAPsfKing.i - King profile CTA point spread function class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2018 by Michael Mayer                               *
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
 * @file GCTAPsfKing.i
 * @brief King profile CTA point spread function class definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAPsfKing.hpp"
//#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAPsfKing
 *
 * @brief CTA point spread function class using a King Profile
 *
 * This class implements the CTA point spread function response as function
 * of energy as determined from a FITS table.
 ***************************************************************************/
class GCTAPsfKing : public GCTAPsf {
public:
    // Constructors and destructors
    GCTAPsfKing(void);
    explicit GCTAPsfKing(const GFilename& filename);
    GCTAPsfKing(const GCTAPsfKing& psf);
    virtual ~GCTAPsfKing(void);

    // Operators
    double operator()(const double& delta,
                      const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void         clear(void);
    GCTAPsfKing* clone(void) const;
    std::string  classname(void) const;
    void         load(const GFilename& filename);
    GFilename    filename(void) const;
    double       mc(GRan&         ran,
                    const double& logE, 
                    const double& theta = 0.0, 
                    const double& phi = 0.0,
                    const double& zenith = 0.0,
                    const double& azimuth = 0.0,
                    const bool&   etrue = true) const;
    double       delta_max(const double& logE, 
                           const double& theta = 0.0, 
                           const double& phi = 0.0,
                           const double& zenith = 0.0,
                           const double& azimuth = 0.0,
                           const bool&   etrue = true) const;
    double       containment_radius(const double& fraction,
                                    const double& logE, 
                                    const double& theta = 0.0, 
                                    const double& phi = 0.0,
                                    const double& zenith = 0.0,
                                    const double& azimuth = 0.0,
                                    const bool&   etrue = true) const;
    // Methods
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFitsTable& table);
    void                     write(GFitsBinTable& table) const;
    void                     save(const GFilename& filename,
                                  const bool&      clobber = false) const;
};


/***********************************************************************//**
 * @brief GCTAPsfKing class extension
 ***************************************************************************/
%extend GCTAPsfKing {
    GCTAPsfKing copy() {
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
