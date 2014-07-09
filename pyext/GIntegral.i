/***************************************************************************
 *                    GIntegral.i - Integration class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GIntegral.hpp
 * @brief Integration class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GIntegral.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GIntegral
 *
 * @brief Integration class Python interface definition.
 *
 * This class allows to perform integration using various methods. The
 * integrand is implemented by a derived class of GIntegrand.
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
    void               max_iter(const int& max_iter);
    void               eps(const double& eps);
    void               silent(const bool& silent);
    const int&         iter(void) const;
    const int&         max_iter(void) const;
    const double&      eps(void) const;
    const bool&        silent(void) const;
    const bool&        isvalid(void) const;
    const std::string& message(void) const;
    void               kernel(GFunction* kernel);
    const GFunction*   kernel(void) const;
    double             romb(const double& a, const double& b, const int& k = 5);
    double             trapzd(const double& a, const double& b, const int& n = 1,
                              double result = 0.0);
};


/***********************************************************************//**
 * @brief GIntegral class extension
 ***************************************************************************/
%extend GIntegral {
    GIntegral copy() {
        return (*self);
    }
};
