/***************************************************************************
 *          GCOMIaq.i - COMPTEL instrument response representation         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMIaq.i
 * @brief COMPTEL instrument response representation class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMIaq.hpp"
%}


/***********************************************************************//**
 * @class GCOMIaq
 *
 * @brief Interface for the COMPTEL instrument response representation class
 ***************************************************************************/
class GCOMIaq : public GBase {

public:
    // Constructors and destructors
    GCOMIaq(void);
    GCOMIaq(const GCOMIaq& iaq);
    GCOMIaq(const double&   phigeo_max, const double& phigeo_bin_size,
            const double&   phibar_max, const double& phibar_bin_size,
            const GEbounds& ebounds);
    ~GCOMIaq(void);

    // Methods
    void        clear(void);
    GCOMIaq*    clone(void) const;
    std::string classname(void) const;
    void        save(const GFilename& filename, const bool& clobber = false) const;
    void        set(const GEnergy& energy);
    void        set(const GModelSpectral& spectrum);
};


/***********************************************************************//**
 * @brief GCOMIaq class extension
 ***************************************************************************/
%extend GCOMIaq {
    GCOMIaq copy() {
        return (*self);
    }
};
