/***************************************************************************
 *                 GLATAeff.i - Fermi/LAT effective area class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GLATAeff.i
 * @brief Fermi/LAT effective area class definition
 * @author J. Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATAeff.hpp"
%}


/***********************************************************************//**
 * @class GLATAeff
 *
 * @brief Interface for the Fermi/LAT effective area
 ***************************************************************************/
class GLATAeff : public GBase {

public:
    // Constructors and destructors
    GLATAeff(void);
    explicit GLATAeff(const std::string& filename);
    GLATAeff(const GLATAeff& aeff);
    virtual ~GLATAeff(void);

    // Operators
    double operator()(const double& logE, const double& ctheta);
    double operator()(const double& logE, const double& ctheta,
                      const double& phi);

    // Methods
    void         clear(void);
    GLATAeff*    clone(void) const;
    void         load(const std::string filename);
    void         save(const std::string filename,
                      const bool& clobber = false);
    void         read(const GFits& file);
    void         write(GFits& file) const;
    int          size(void) const;
    int          nenergies(void) const;
    int          ncostheta(void) const;
    double       costhetamin(void) const;
    void         costhetamin(const double& ctheta);
    bool         has_phi(void) const;
    bool         has_efficiency(void) const;
    double       efficiency_factor1(const GEnergy& srcEng) const;
    double       efficiency_factor2(const GEnergy& srcEng) const;
};


/***********************************************************************//**
 * @brief GLATAeff class extension
 ***************************************************************************/
%extend GLATAeff {
    GLATAeff copy() {
        return (*self);
    }
};
