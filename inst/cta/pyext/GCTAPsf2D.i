/***************************************************************************
 *            GCTAPsf2D.i - CTA 2D point spread function class             *
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
 * @file GCTAPsf2D.i
 * @brief CTA 2D point spread function class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAPsf2D.hpp"
%}


/***********************************************************************//**
 * @class GCTAPsf2D
 *
 * @brief CTA 2D point spread function class
 ***************************************************************************/
class GCTAPsf2D : public GCTAPsf {

public:
    // Constructors and destructors
    GCTAPsf2D(void);
    explicit GCTAPsf2D(const GFilename& filename);
    GCTAPsf2D(const GCTAPsf2D& psf);
    virtual ~GCTAPsf2D(void);

    // Operators
    double operator()(const double& delta,
                      const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void        clear(void);
    GCTAPsf2D*  clone(void) const;
    std::string classname(void) const;
    void        load(const GFilename& filename);
    GFilename   filename(void) const;
    double      mc(GRan&         ran,
                   const double& logE, 
                   const double& theta = 0.0, 
                   const double& phi = 0.0,
                   const double& zenith = 0.0,
                   const double& azimuth = 0.0,
                   const bool&   etrue = true) const;
    double      delta_max(const double& logE, 
                          const double& theta = 0.0, 
                          const double& phi = 0.0,
                          const double& zenith = 0.0,
                          const double& azimuth = 0.0,
                          const bool&   etrue = true) const;
    double      containment_radius(const double& fraction,
                                   const double& logE, 
                                   const double& theta = 0.0, 
                                   const double& phi = 0.0,
                                   const double& zenith = 0.0,
                                   const double& azimuth = 0.0,
                                   const bool&   etrue = true) const;
    // Methods
    const GCTAResponseTable&   table(void) const;
    void                       table(const GCTAResponseTable& table);
    void                       read(const GFitsTable& table);
    void                       write(GFitsBinTable& table) const;
    void                       save(const GFilename& filename,
                                    const bool&      clobber = false) const;
    
};


/***********************************************************************//**
 * @brief GCTAPsf2D class extension
 ***************************************************************************/
%extend GCTAPsf2D {
    GCTAPsf2D copy() {
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
