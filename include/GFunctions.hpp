/***************************************************************************
 *     GFunctions.hpp - Single parameter functions abstract base class     *
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
 * @file GFunctions.hpp
 * @brief Single parameter functions abstract base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFUNCTIONS_HPP
#define GFUNCTIONS_HPP

/* __ Includes ___________________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GVector;


/***********************************************************************//**
 * @class GFunctions
 *
 * @brief Single parameter functions abstract base class
 *
 * This class implements the abstract interface for a set of one parameter
 * functions.
 *
 * The set of functions is for example used for integration or numerical
 * computation of derivatives. This class has no members. The only pure
 * virtual method that needs to be implemented by the derived class is the
 * eval() method that evaluates the functions at a given value x and returns
 * the evaluation result as an GNdarray object.
 ***************************************************************************/
class GFunctions {

public:

    // Constructors and destructors
    GFunctions(void);
    GFunctions(const GFunctions& functions);
    virtual ~GFunctions(void);

    // Operators
    GFunctions& operator=(const GFunctions& functions);

    // Methods
    virtual int     size(void) const = 0;
    virtual GVector eval(const double& x) = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GFunctions& functions);
    void free_members(void);
};

#endif /* GFUNCTIONS_HPP */
