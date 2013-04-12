/***************************************************************************
 *                 GCTAAeffArf.i - CTA ARF effective area class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCTAAeffArf.hpp
 * @brief CTA ARF effective area class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAAeffArf.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAAeffArf
 *
 * @brief CTA performance table effective area class
 ***************************************************************************/
class GCTAAeffArf : public GCTAAeff {

public:
    // Constructors and destructors
    GCTAAeffArf(void);
    GCTAAeffArf(const std::string& filename);
    GCTAAeffArf(const GCTAAeffArf& cta);
    virtual ~GCTAAeffArf(void);

    // Operators
    double operator()(const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void         clear(void);
    GCTAAeffArf* clone(void) const;
    void         load(const std::string& filename);
    std::string  filename(void) const;

    // Methods
    int           size(void) const;
    void          sigma(const double& sigma);
    const double& sigma(void) const;
    void          thetacut(const double& thetacut);
    const double& thetacut(void) const;
    void          scale(const double& scale);
    const double& scale(void) const;
    void          read_arf(const GFitsTable* hdu);
};


/***********************************************************************//**
 * @brief GCTAAeffArf class extension
 ***************************************************************************/
%extend GCTAAeffArf {
    GCTAAeffArf copy() {
        return (*self);
    }
};
