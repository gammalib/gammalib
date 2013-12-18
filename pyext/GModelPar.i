/***************************************************************************
 *                    GModelPar.i - Model parameter class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GModelPar.i
 * @brief Model parameter class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelPar.hpp"
%}


/***********************************************************************//**
 * @class GModelPar
 *
 * @brief Model parameter class
 ***************************************************************************/
class GModelPar : public GOptimizerPar {
public:
    // Constructors and destructors
    GModelPar(void);
    GModelPar(const std::string& name, const double& value);
    GModelPar(const std::string& name, const double& factor,
              const double& scale);
    GModelPar(const GModelPar& par);
    virtual ~GModelPar(void);

    // Methods
    GModelPar* clone(void) const;
    void       read(const GXmlElement& xml);
    void       write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GModelPar class extension
 ***************************************************************************/
%extend GModelPar {
    GModelPar copy() {
        return (*self);
    }
};
