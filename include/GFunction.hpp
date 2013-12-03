/***************************************************************************
 *      GFunction.hpp - Single parameter function abstract base class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @brief Single parameter function abstract base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFUNCTION_HPP
#define GFUNCTION_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GFunction
 *
 * @brief Single parameter function abstract base class
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
    GFunction(const GFunction& function);
    virtual ~GFunction(void);

    // Operators
    GFunction& operator=(const GFunction& function);

    // Methods
    virtual double eval(const double& x) = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GFunction& function);
    void free_members(void);
};

#endif /* GFUNCTION_HPP */
