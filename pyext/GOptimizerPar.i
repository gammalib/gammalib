/***************************************************************************
 *                GOptimizerPar.i - Optimizer parameter class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GOptimizerPar.i
 * @brief Optimizer parameter class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GOptimizerPar.hpp"
%}


/***********************************************************************//**
 * @class GOptimizerPar
 *
 * @brief Model parameter class
 ***************************************************************************/
class GOptimizerPar : public GBase {

public:
    // Constructors and destructors
    GOptimizerPar(void);
    GOptimizerPar(const std::string& name, const double& value);
    GOptimizerPar(const std::string& name, const double& factor,
                  const double& scale);
    virtual ~GOptimizerPar(void);

    // Attribute methods 
    double value(void) const;
    double error(void) const;
    double gradient(void) const;
    double min(void) const;
    double max(void) const;
    void   value(const double& value);
    void   error(const double& error);
    void   gradient(const double& gradient);
    void   min(const double& min);
    void   max(const double& max);
    void   range(const double& min, const double& max);

    // Factorization methods
    const double& factor_value(void) const;
    const double& factor_error(void) const;
    const double& factor_gradient(void) const;
    const double& factor_min(void) const;
    const double& factor_max(void) const;
    const double& scale(void) const;
    void          factor_value(const double& value);
    void          factor_error(const double& error);
    void          factor_gradient(const double& gradient);
    void          factor_min(const double& min);
    void          factor_max(const double& max);
    void          factor_range(const double& min, const double& max);
    void          scale(const double& scale);

    // Boundary methods
    bool has_min(void) const;
    bool has_max(void) const;
    bool has_range(void) const;
    void remove_min(void);
    void remove_max(void);
    void remove_range(void);

    // Property methods
    bool is_free(void) const;
    bool is_fixed(void) const;
    bool has_grad(void) const;
    void free(void);
    void fix(void);
    void has_grad(const bool& grad);

    // Other methods
    void               clear(void);
    GOptimizerPar*     clone(void) const;
    const std::string& name(void) const;
    const std::string& unit(void) const;
    void               name(const std::string& name);
    void               unit(const std::string& unit);
    void               autoscale(void);
};


/***********************************************************************//**
 * @brief GOptimizerPar class extension
 ***************************************************************************/
%extend GOptimizerPar {
    GOptimizerPar copy() {
        return (*self);
    }
};
