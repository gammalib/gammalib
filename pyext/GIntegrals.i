/***************************************************************************
 *           GIntegrals.i - Integration class for set of functions         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GIntegrals.i
 * @brief Integration class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GIntegrals.hpp"
%}


/***********************************************************************//**
 * @class GIntegrals
 *
 * @brief Integration class for set of functions
 *
 * This class allows to perform integration of a set of functions. The
 * integrand is implemented by a derived class of GFunctions.
 ***************************************************************************/
class GIntegrals : public GBase {
public:

    // Constructors and destructors
    explicit GIntegrals(void);
    explicit GIntegrals(GFunctions* kernels);
    GIntegrals(const GIntegrals& integral);
    virtual ~GIntegrals(void);

    // Methods
    void               clear(void);
    GIntegrals*        clone(void) const;
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
    void               kernels(GFunctions* kernels);
    const GFunctions*  kernels(void) const;
    GNdarray           romberg(std::vector<double> bounds,
                               const int&          order = 5);
    GNdarray           romberg(const double& a,
                               const double& b,
                               const int&    order = 5);
    GNdarray           trapzd(const double& a,
                              const double& b,
                              const int&    n,
                              GNdarray      result);
};


/***********************************************************************//**
 * @brief GIntegrals class extension
 ***************************************************************************/
%extend GIntegrals {
    GIntegrals copy() {
        return (*self);
    }
};
