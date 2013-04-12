/***************************************************************************
 *                         GSource.i - Source class                        *
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
 * @file GSource.i
 * @brief Source class interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSource.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GSource
 *
 * @brief Class that handles gamma-ray sources
 ***************************************************************************/
class GSource : public GBase {
public:
    // Constructors and destructors
    GSource(void);
    explicit GSource(const std::string& name,
                     GModelSpatial*     model,
                     const GEnergy&     energy,
                     const GTime&       time);
    GSource(const GSource& src);
    virtual ~GSource(void);

    // Methods
    void                 clear(void);
    GSource*             clone(void) const;
    const std::string&   name(void) const;
    const GModelSpatial* model(void) const;
    const GEnergy&       energy(void) const;
    const GTime&         time(void) const;
    void                 name(const std::string& name);
    void                 model(GModelSpatial* model);
    void                 energy(const GEnergy& energy);
    void                 time(const GTime& time);
};


/***********************************************************************//**
 * @brief GSource class extension
 ***************************************************************************/
%extend GSource {
    GSource copy() {
        return (*self);
    }
};
