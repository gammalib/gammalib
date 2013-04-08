/***************************************************************************
 *                GMatrixBase.hpp - Abstract matrix base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2013 by Juergen Knoedlseder                         *
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
#include "GBase.hpp"
#include "GVector.hpp"


/***********************************************************************//**
 * @class GMatrixBase
 *
 * @brief Abstract matrix base class interface defintion
 *
 * This class is an abstract base class from which all matrix classes will
 * be derived. It defines the common interface for all matrix classes and
 * implements common methods and data members.
 ***************************************************************************/
class GMatrixBase : public GBase {

public:
    // Constructors and destructors
    GMatrixBase(void);
    GMatrixBase(const GMatrixBase& matrix);
    virtual ~GMatrixBase(void);

    // Pure virtual operators
    virtual double&       operator()(const int& row, const int& column) = 0;
    virtual const double& operator()(const int& row, const int& column) const = 0;
    virtual GVector       operator*(const GVector& vector) const = 0;

    // Base class operators
    virtual GMatrixBase&  operator=(const GMatrixBase& matrix);
    virtual bool          operator==(const GMatrixBase& matrix) const;
    virtual bool          operator!=(const GMatrixBase& matrix) const;

    // Pure virtual methods
    virtual void         clear(void) = 0;
    virtual GMatrixBase* clone(void) const = 0;
    virtual GVector      row(const int& row) const = 0;
    virtual void         row(const int& row, const GVector& vector) = 0;
    virtual GVector      column(const int& column) const = 0;
    virtual void         column(const int& column, const GVector& vector) = 0;
    virtual void         add_to_row(const int& row, const GVector& vector) = 0;
    virtual void         add_to_column(const int& column, const GVector& vector) = 0;
    virtual void         transpose(void) = 0;
    virtual void         invert(void) = 0;
    virtual void         negate(void) = 0;
    virtual void         abs(void) = 0;
    virtual double       fill(void) const = 0;
    virtual double       min(void) const = 0;
    virtual double       max(void) const = 0;
    virtual double       sum(void) const = 0;
    virtual std::string  print(const GChatter& chatter = NORMAL) const = 0;

    // Base class methods
    const int& size(void) const;
    const int& columns(void) const;
    const int& rows(void) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GMatrixBase& matrix);
    void        free_members(void);
    void        select_non_zero(void);
    void        scale_elements(const double& scalar);
    void        set_all_elements(const double& value);
    double      get_min_element(void) const;
    double      get_max_element(void) const;
    double      get_element_sum(void) const;
    std::string print_elements(const GChatter& chatter = NORMAL,
                               const int&      num = 10) const;
    std::string print_row_compression(const GChatter& chatter = NORMAL) const;
    std::string print_col_compression(const GChatter& chatter = NORMAL) const;

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


/***********************************************************************//**
 * @brief Return number of matrix elements
 *
 * @return Number of matrix elements.
 *
 * Returns the number of matrix elements that are stored in the matrix.
 ***************************************************************************/
inline
const int& GMatrixBase::size(void) const
{
    return (m_elements);
}


/***********************************************************************//**
 * @brief Return number of matrix columns
 *
 * @return Number of matrix columns.
 *
 * Returns the number of matrix columns.
 ***************************************************************************/
inline
const int& GMatrixBase::columns(void) const
{
    return (m_cols);
}


/***********************************************************************//**
 * @brief Return number of matrix rows
 *
 * @return Number of matrix rows.
 *
 * Returns the number of matrix rows.
 ***************************************************************************/
inline
const int& GMatrixBase::rows(void) const
{
    return (m_rows);
}

#endif /* GMATRIXBASE_HPP */
