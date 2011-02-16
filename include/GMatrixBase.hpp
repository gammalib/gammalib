/***************************************************************************
 *               GMatrixBase.hpp  -  matrix abstract base class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GMatrixBase.hpp
 * @brief GMatrixBase class definition.
 * @author J. Knodlseder
 */

#ifndef GMATRIXBASE_HPP
#define GMATRIXBASE_HPP

/* __ Includes ___________________________________________________________ */
#include "GVector.hpp"


/***********************************************************************//**
 * @class GMatrixBase
 *
 * @brief GMatrixBase class interface defintion
 *
 * GMatrixBase is an abstract base class for all matrix classes. It defines
 * the common interfaces of the matrix objects and provides some common
 * services to the derived classes.
 ***************************************************************************/
class GMatrixBase {

public:
    // Constructors and destructors
    GMatrixBase(void);
    GMatrixBase(const GMatrixBase& m);
    virtual ~GMatrixBase(void);

    // Operators
    virtual GMatrixBase&  operator= (const GMatrixBase& m);
    virtual int           operator== (const GMatrixBase& m) const;
    virtual int           operator!= (const GMatrixBase& m) const;
    virtual double&       operator() (int row, int col) = 0;
    virtual const double& operator() (int row, int col) const = 0;
    virtual GVector       operator* (const GVector& v) const = 0;

    // Methods
    virtual int     cols(void) const { return m_cols; }
    virtual int     rows(void) const { return m_rows; }
    virtual void    clear() = 0;
    virtual void    transpose() = 0;
//    virtual void    invert() = 0;
    virtual GVector extract_row(int row) const = 0;
    virtual GVector extract_col(int col) const = 0;
    virtual double  fill(void) const = 0;
    virtual double  min(void) const = 0;
    virtual double  max(void) const = 0;
    virtual double  sum(void) const = 0;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GMatrixBase& m);
    void   free_members(void);
    void   select_non_zero(void);
    void   negation(void);
    void   addition(const GMatrixBase& m);
    void   subtraction(const GMatrixBase& m);
    void   multiplication(const double& s);
    void   set_all_elements(const double& s);
    double get_min_element(void) const;
    double get_max_element(void) const;
    double get_element_sum(void) const;
    void   dump_elements(std::ostream& os) const;
    void   dump_row_comp(std::ostream& os) const;
    void   dump_col_comp(std::ostream& os) const;

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
