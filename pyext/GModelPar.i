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
 * @brief Model parameter class Python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelPar.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelPar
 *
 * @brief Model parameter class
 ***************************************************************************/
class GModelPar : public GBase {

public:
    // Constructors and destructors
    GModelPar(void);
    GModelPar(const std::string& name, const double& value, 
              const double& scale = 1.0);
    GModelPar(const GModelPar& par);
    virtual ~GModelPar(void);

    // Methods
    void        clear(void);
    GModelPar*  clone(void) const;
    std::string name(void) const;
    std::string unit(void) const;
    double      real_value(void) const;
    double      real_error(void) const;
    double      real_gradient(void) const;
    double      real_min(void) const;
    double      real_max(void) const;
    double      value(void) const;
    double      error(void) const;
    double      gradient(void) const;
    double      min(void) const;
    double      max(void) const;
    double      scale(void) const;
    bool        isfree(void) const;
    bool        hasmin(void) const;
    bool        hasmax(void) const;
    bool        hasgrad(void) const;
    void        name(const std::string& name);
    void        unit(const std::string& unit);
    void        real_value(const double& value);
    void        real_error(const double& error);
    void        value(const double& value);
    void        error(const double& error);
    void        gradient(const double& gradient);
    void        min(const double& min);
    void        max(const double& max);
    void        scale(const double& scale);
    void        range(const double& min, const double& max);
    void        remove_min(void);
    void        remove_max(void);
    void        remove_range(void);
    void        free(void);
    void        fix(void);
    void        hasgrad(const bool& grad);
    void        autoscale(void);
    void        read(const GXmlElement& xml);
    void        write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GModelPar class extension
 ***************************************************************************/
%extend GModelPar {
    char *__str__() {
        return tochar(self->print());
    }
    GModelPar copy() {
        return (*self);
    }
};
