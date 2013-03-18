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
    explicit GModelPar(const std::string& name, const double& value, int i);
    explicit GModelPar(const std::string& name, const double& factor, const double& scale);
    virtual ~GModelPar(void);

    // Attribute methods 
    double      Value(void) const;
    double      Error(void) const { return m_factor_error*m_scale; }
    double      Gradient(void) const { return m_factor_gradient*m_scale; }
    double      Min(void) const { return m_factor_min*m_scale; }
    double      Max(void) const { return m_factor_max*m_scale; }
    void        Value(const double& value);
    void        Error(const double& error);
    void        Gradient(const double& min);
    void        Min(const double& min);
    void        Max(const double& max);
    void        Range(const double& min, const double& max);

    // Factorization methods
    const double& factor_value(void) const { return m_factor_value; }
    const double& factor_error(void) const { return m_factor_error; }
    const double& factor_gradient(void) const { return m_factor_gradient; }
    const double& factor_min(void) const { return m_factor_min; }
    const double& factor_max(void) const { return m_factor_max; }
    const double& scale(void) const { return m_scale; }
    void          factor_value(const double& value);
    void          factor_error(const double& error) { m_factor_error=error; }
    void          factor_gradient(const double& gradient) { m_factor_gradient=gradient; }
    void          factor_min(const double& min);
    void          factor_max(const double& max);
    void          factor_range(const double& min, const double& max);
    void          Scale(const double& scale);

    // Boundary methods
    bool        hasmin(void) const { return m_hasmin; }
    bool        hasmax(void) const { return m_hasmax; }
    bool        hasrange(void) const { return m_hasmin && m_hasmax; }
    void        remove_min(void) { m_hasmin=false; }
    void        remove_max(void) { m_hasmax=false; }
    void        remove_range(void) { m_hasmin=false; m_hasmax=false; }

    // Property methods
    bool        isfree(void) const { return m_free; }
    bool        isfixed(void) const { return !m_free; }
    bool        hasgrad(void) const { return m_hasgrad; }
    void        free(void) { m_free=true; }
    void        fix(void) { m_free=false; }
    void        hasgrad(const bool& grad) { m_hasgrad=grad; }

    // Other methods
    void        clear(void);
    GModelPar*  clone(void) const;
    std::string name(void) const { return m_name; }
    std::string unit(void) const { return m_unit; }
    void        name(const std::string& name) { m_name=name; }
    void        unit(const std::string& unit) { m_unit=unit; }
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
