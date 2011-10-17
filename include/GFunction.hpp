/***************************************************************************
 *     GFunction.hpp  -  Single parameter function abstract base class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GFunction.hpp
 * @brief GFunction abstract virtual base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GFUNCTION_HPP
#define GFUNCTION_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GFunction
 *
 * @brief GFunction class interface defintion.
 *
 * This class implements the abstract interface for a one parameter function.
 * This function is for example used for integration or numerical computation
 * of derivatives. This class has no members. The only pure virtual method
 * that needs to be implemented by the derived class is the eval() method
 * that provides function evaluation at a given value x, e.g. y=eval(x).
 ***************************************************************************/
class GFunction {

public:

    // Constructors and destructors
    GFunction(void);
    GFunction(const GFunction& func);
    virtual ~GFunction(void);

    // Operators
    GFunction& operator= (const GFunction& func);

    // Methods
    virtual double eval(double x) = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GFunction& func);
    void free_members(void);
};

#endif /* GFUNCTION_HPP */
