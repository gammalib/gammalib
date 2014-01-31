/***************************************************************************
 *                  GRan.i - Random number generator class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GRan.i
 * @brief Random number generator class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GRan.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GRan
 *
 * @brief Random number generator class
 *
 * The GRan class implements a random number generator that is inspired from
 * the Ran structure given in Numerical Recipes, Third Edition (p. 341ff).
 ***************************************************************************/
class GRan : public GBase {
public:
    // Constructors and destructors
    GRan(void);
    GRan(unsigned long long int seed);
    GRan(const GRan& ran);
    virtual ~GRan(void);
 
    // Methods
    void                   clear(void);
    GRan*                  clone(void) const;
    void                   seed(unsigned long long int seed);
    unsigned long long int seed(void) const;
    unsigned long int      int32(void);
    unsigned long long int int64(void);
    double                 uniform(void);
    double                 normal(void);
    double                 exp(const double& arg);
    double                 poisson(const double& arg);
    double                 chisq2(void);
};


/***********************************************************************//**
 * @brief GRan class extension
 ***************************************************************************/
%extend GRan {
    GRan copy() {
        return (*self);
    }
};
