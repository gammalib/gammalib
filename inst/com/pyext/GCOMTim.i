/***************************************************************************
 *              GCOMTim.i - COMPTEL Good Time Intervals class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knodlseder                               *
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
 * @file GCOMTim.i
 * @brief COMPTEL Good Time Intervals class definition
 * @author Juergen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMTim.hpp"
%}


/***********************************************************************//**
 * @class GCOMTim
 *
 * @brief COMPTEL Good Time Intervals class
 ***************************************************************************/
class GCOMTim : public GBase {

public:
    // Constructors and destructors
    GCOMTim(void);
    GCOMTim(const GCOMTim& tim);
    GCOMTim(const GFilename& filename, const std::string& usage = "",
                                       const std::string& mode  = "");
    virtual ~GCOMTim(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMTim*    clone(void) const;
    virtual std::string classname(void) const;

    // Other methods
    bool        contains(const GTime& time) const;
    const GGti& gti(void) const;
    void        gti(const GGti& gti);
    void        read(const GFitsTable& table, const std::string& usage = "",
                                              const std::string& mode  = "");
};


/***********************************************************************//**
 * @brief GCOMTim class extension
 ***************************************************************************/
%extend GCOMTim {
    GCOMTim copy() {
        return (*self);
    }
};
