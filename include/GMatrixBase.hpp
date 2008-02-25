/***************************************************************************
 *               GMatrixBase.hpp  -  matrix abstract base class            *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2006 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GMATRIXBASE_HPP
#define GMATRIXBASE_HPP

/* __ Includes ___________________________________________________________ */
#include "GVector.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                        GMatrixBase class definition                     *
 * ----------------------------------------------------------------------- *
 * GMatrixBase is an abstract base class for all matrix classes. It        *
 * defines the common interfaces of the matrix objects and provides some   *
 * common services to the derived classes.                                 *
 ***************************************************************************/
class GMatrixBase {

// Methods that are available to everybody
public:
  // Constructors and destructors
  GMatrixBase();
  GMatrixBase(const GMatrixBase& m);
  virtual ~GMatrixBase();

  // Operators (first virtual, then pure virtual)
  virtual GMatrixBase&  operator= (const GMatrixBase& m);
  virtual int           operator== (const GMatrixBase& m) const;
  virtual int           operator!= (const GMatrixBase& m) const;
  //
  virtual double&       operator() (int row, int col) = 0;
  virtual const double& operator() (int row, int col) const = 0;
  //virtual GMatrixBase   operator+ (const GMatrixBase& m) const = 0;
  //virtual GMatrixBase   operator- (const GMatrixBase& m) const = 0;
  //virtual GMatrixBase   operator* (const GMatrixBase& m) const = 0;
  virtual GVector       operator* (const GVector& v) const = 0;
  //virtual GMatrixBase   operator- () const = 0;
  //virtual GMatrixBase&  operator+= (const GMatrixBase& m) = 0;
  //virtual GMatrixBase&  operator-= (const GMatrixBase& m) = 0;
  //virtual GMatrixBase&  operator*= (const GMatrixBase& m) = 0;
  //virtual GMatrixBase&  operator*= (const double& s) = 0;
  //virtual GMatrixBase&  operator/= (const double& s) = 0;

  // Methods (first virtual, then pure virtual)
  virtual int     cols(void) const { return m_cols; }
  virtual int     rows(void) const { return m_rows; }
  //
  //virtual void    add_col(const GVector& v, int col) = 0;
  //virtual void    cholesky_decompose(int compress = 1) = 0;
  //virtual GVector cholesky_solver(const GVector& v, int compress = 1) = 0;
  //virtual void    cholesky_invert(int compress = 1) = 0;
  virtual void    clear() = 0;
  virtual GVector extract_row(int row) const = 0;
  virtual GVector extract_col(int col) const = 0;
  //virtual GMatrix extract_lower_triangle() const = 0;
  //virtual GMatrix extract_upper_triangle() const = 0;
  //virtual void    insert_col(const GVector& v, int col) = 0;
  virtual double  fill() const = 0;
  virtual double  min() const = 0;
  virtual double  max() const = 0;
  virtual double  sum() const = 0;
  virtual void    transpose() = 0;
  
// Methods and data that are available to derived classes
protected:
  // Protected methods
  void   select_non_zero();
  void   negation();
  void   addition(const GMatrixBase& m);
  void   subtraction(const GMatrixBase& m);
  void   multiplication(const double& s);
  void   set_all_elements(const double& s);
  double get_min_element() const;
  double get_max_element() const;
  double get_element_sum() const;
  void   dump_elements(ostream& os) const;
  void   dump_row_comp(ostream& os) const;
  void   dump_col_comp(ostream& os) const;

  // Protected data area
  int     m_rows;         // Number of rows
  int     m_cols;         // Number of columns
  int     m_elements;     // Number of elements stored in matrix
  int     m_alloc;        // Size of allocated matrix memory
  int     m_num_rowsel;   // Number of selected rows (for compressed decomposition)
  int     m_num_colsel;   // Number of selected columns (for compressed decomposition)
  int*    m_colstart;     // Column start indices (m_cols+1)
  int*    m_rowsel;       // Row selection (for compressed decomposition)
  int*    m_colsel;       // Column selection (for compressed decomposition)
  double* m_data;         // Matrix data

// Methods that are available to the base class only
private:
  void init_members(void);
  void copy_members(const GMatrixBase& m);
  void free_members(void);
};

#endif /* GMATRIXBASE_HPP */
