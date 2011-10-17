/***************************************************************************
 *          GMatrixBase.hpp  -  Matrix abstract base class SWIG file       *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2006 by Jurgen Knodlseder                   *
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
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMatrixBase.hpp"
%}


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
  virtual int           operator== (const GMatrixBase& m) const;
  virtual int           operator!= (const GMatrixBase& m) const;
  //
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
};
