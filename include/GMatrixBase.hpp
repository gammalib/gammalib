/***************************************************************************
 *               GMatrixBase.hpp  -  Abstract matrix base class            *
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
 * @file GMatrixBase.hpp
 * @brief Abstract matrix base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMATRIXBASE_HPP
#define GMATRIXBASE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GVector.hpp"


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
    virtual double&       operator()(const int& row, const int& col) = 0;
    virtual const double& operator()(const int& row, const int& col) const = 0;
    virtual GVector       operator*(const GVector& vector) const = 0;

    // Implemented base class operators
    virtual GMatrixBase&  operator=(const GMatrixBase& matrix);
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
    virtual std::string print(void) const = 0;

    // Implemented base class methods
    virtual int         cols(void) const { return m_cols; }
    virtual int         rows(void) const { return m_rows; }

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GMatrixBase& matrix);
    void        free_members(void);
    void        select_non_zero(void);
    void        negation(void);
    void        addition(const GMatrixBase& matrix);
    void        subtraction(const GMatrixBase& matrix);
    void        multiplication(const double& scalar);
    void        set_all_elements(const double& value);
    double      get_min_element(void) const;
    double      get_max_element(void) const;
    double      get_element_sum(void) const;
    std::string print_elements(void) const;
    std::string print_row_compression(void) const;
    std::string print_col_compression(void) const;

    // Protected data area
    int     m_rows;       //!< Number of rows
    int     m_cols;       //!< Number of columns
    int     m_elements;   //!< Number of elements stored in matrix
    int     m_alloc;      //!< Size of allocated matrix memory
    int     m_num_rowsel; //!< Number of selected rows (for compressed decomposition)
    int     m_num_colsel; //!< Number of selected columns (for compressed decomposition)
    int*    m_colstart;   //!< Column start indices (m_cols+1)
    int*    m_rowsel;     //!< Row selection (for compressed decomposition)
    int*    m_colsel;     //!< Column selection (for compressed decomposition)
    double* m_data;       //!< Matrix data
};

#endif /* GMATRIXBASE_HPP */
