/***************************************************************************
 *               GCTABackground3D.i - CTA 3D background class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2018 by Juergen Knoedlseder                         *
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
 * @file GCTABackground3D.i
 * @brief CTA 3D background class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTABackground3D.hpp"
%}


/***********************************************************************//**
 * @class GCTABackground3D
 *
 * @brief CTA 3D background class
 ***************************************************************************/
class GCTABackground3D : public GCTABackground {

public:
    // Constructors and destructors
    GCTABackground3D(void);
    explicit GCTABackground3D(const GFilename& filename);
    GCTABackground3D(const GCTABackground3D& bgd);
    virtual ~GCTABackground3D(void);

    // Implemented pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety) const;

    // Implemented pure virtual methods
    void                       clear(void);
    GCTABackground3D*          clone(void) const;
    std::string                classname(void) const;
    void                       load(const GFilename& filename);
    GFilename                  filename(void) const;
    GCTAInstDir                mc(const GEnergy& energy,
                                  const GTime& time,
                                  GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;
    double                     rate_ebin(const GCTAInstDir& dir,
                                         const GEnergy&     emin,
                                         const GEnergy&     emax) const;

    // Methods
    bool                       is_valid(void) const;
    const GCTAResponseTable&   table(void) const;
    void                       table(const GCTAResponseTable& table);
    void                       read(const GFitsTable& table);
    void                       write(GFitsBinTable& table) const;
    void                       save(const GFilename& filename,
                                    const bool&      clobber = false) const;
};


/***********************************************************************//**
 * @brief GCTABackground3D class extension
 ***************************************************************************/
%extend GCTABackground3D {
    GCTABackground3D copy() {
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
