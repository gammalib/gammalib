/***************************************************************************
 *                GMatrixBase.i  -  Abstract matrix base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2012 by Juergen Knoedlseder                         *
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
 * @file GMatrixBase.i
 * @brief Abstract matrix base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMatrixBase.hpp"
%}


/***********************************************************************//**
 * @class GMatrixBase
 *
 * @brief Abstract matrix base class interface defintion
 *
 * This is an abstract base class for all matrix classes. It defines the
 * common interface of the matrix objects and provides some common services
 * to the derived classes.
 ***************************************************************************/
class GMatrixBase {
public:
    // Constructors and destructors
    GMatrixBase(void);
    GMatrixBase(const GMatrixBase& matrix);
    virtual ~GMatrixBase(void);

    // Pure virtual operators
    virtual GVector       operator*(const GVector& vector) const = 0;

    // Implemented base class operators
    virtual bool          operator==(const GMatrixBase& matrix) const;
    virtual bool          operator!=(const GMatrixBase& matrix) const;

    // Pure virtual methods
    virtual void        clear() = 0;
    virtual void        transpose() = 0;
    virtual void        invert() = 0;
    virtual void        add_col(const GVector& vector, const int& col) = 0;
    virtual void        insert_col(const GVector& vector, const int& col) = 0;
    virtual GVector     extract_row(const int& row) const = 0;
    virtual GVector     extract_col(const int& col) const = 0;
    virtual double      fill(void) const = 0;
    virtual double      min(void) const = 0;
    virtual double      max(void) const = 0;
    virtual double      sum(void) const = 0;

    // Implemented base class methods
    virtual int         cols(void) const;
    virtual int         rows(void) const;
};
