/***************************************************************************
 *       GIntegrand.hpp  -  Integrand function abstract base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
 * @file GIntegrand.hpp
 * @brief Integrand abstract virtual base class Python interface definition.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GIntegrand.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GIntegrand
 *
 * @brief Integrand class Python interface defintion.
 *
 * This class implements the abstract interface for the integration kernel
 * function. It has no members. The only pure virtual method that has to
 * be implemented by the derived class is the eval() method that provides
 * function evaluation at a given value x, e.g. y=eval(x).
 ***************************************************************************/
class GIntegrand {
public:

    // Constructors and destructors
    GIntegrand(void);
    GIntegrand(const GIntegrand& arg);
    virtual ~GIntegrand(void);

    // Methods
    virtual double eval(double x) = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GIntegrand& arg);
    void free_members(void);
};


/***********************************************************************//**
 * @brief GIntegrand class extension
 ***************************************************************************/
//%extend GIntegrand {
//};
