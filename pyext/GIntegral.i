/***************************************************************************
 *                    GIntegral.i - Integration class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2023 by Juergen Knoedlseder                         *
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
 * @file GIntegral.i
 * @brief Integration class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GIntegral.hpp"
%}


/***********************************************************************//**
 * @class GIntegral
 *
 * @brief Integration class Python interface definition.
 *
 * This class allows to perform integration using various methods. The
 * integrand is implemented by a derived class of GFunction.
 ***************************************************************************/
class GIntegral : public GBase {
public:

    // Constructors and destructors
    explicit GIntegral(void);
    explicit GIntegral(GFunction* kernel);
    GIntegral(const GIntegral& integral);
    virtual ~GIntegral(void);

    // Methods
    void               clear(void);
    GIntegral*         clone(void) const;
    std::string        classname(void) const;
    void               max_iter(const int& iter);
    const int&         max_iter(void) const;
    void               fixed_iter(const int& iter);
    const int&         fixed_iter(void) const;
    void               eps(const double& eps);
    const double&      eps(void) const;
    void               silent(const bool& silent);
    const bool&        silent(void) const;
    const int&         iter(void) const;
    const int&         calls(void) const;
    const bool&        is_valid(void) const;
    const std::string& message(void) const;
    void               kernel(GFunction* kernel);
    const GFunction*   kernel(void) const;
    double             romberg(std::vector<double> bounds,
                               const int& order = 5);
    double             romberg(const double& a, const double& b,
                               const int& order = 5);
    double             trapzd(const double& a, const double& b,
                              const int& n = 1, double result = 0.0);
    double             adaptive_simpson(const double& a, const double& b) const;
    double             adaptive_gauss_kronrod(const double& a, const double& b) const;
    double             gauss_kronrod(const double& a, const double& b) const;
};


/***********************************************************************//**
 * @brief GIntegral class extension
 ***************************************************************************/
%extend GIntegral {
    GIntegral copy() {
        return (*self);
    }
};
