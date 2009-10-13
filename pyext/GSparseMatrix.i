/***************************************************************************
 *            GSparseMatrix.i  -  Sparse Matrix class SWIG file            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2009 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSparseMatrix.hpp"
%}

%feature("notabstract") GSparseMatrix;

/***************************************************************************
 *                      GSparseMatrix class definition                     *
 * ----------------------------------------------------------------------- *
 * GMatrix implements a sparse matrix storage class. It derives from the   *
 * abstract base class GMatrixBase.                                        *
 ***************************************************************************/
class GSparseMatrix : public GMatrixBase {
public:
    // Constructors and destructors
    GSparseMatrix(int rows, int cols, int elements = 0);
    GSparseMatrix(const GSparseMatrix& m);
    virtual ~GSparseMatrix();

    // Operators
    //GSparseMatrix& operator= (const GSparseMatrix& m);
    //double&        operator() (int row, int col);
    //const double&  operator() (int row, int col) const;
    GSparseMatrix  operator+ (const GSparseMatrix& m) const;
    GSparseMatrix  operator- (const GSparseMatrix& m) const;
    GSparseMatrix  operator* (const GSparseMatrix& m) const;
    GVector        operator* (const GVector& v) const;
    int            operator== (const GSparseMatrix& m) const;
    int            operator!= (const GSparseMatrix& m) const;
    GSparseMatrix  operator- () const;
    GSparseMatrix& operator+= (const GSparseMatrix& m);
    GSparseMatrix& operator-= (const GSparseMatrix& m);
    GSparseMatrix& operator*= (const GSparseMatrix& m);
    GSparseMatrix& operator*= (const double& d);
    GSparseMatrix& operator/= (const double& d);

    // Matrix methods
    void          add_col(const GVector& v, int col);
    void          add_col(const double* values, const int* rows, int number, int col);
    void          cholesky_decompose(int compress = 1);
    GVector       cholesky_solver(const GVector& v, int compress = 1);
    void          cholesky_invert(int compress = 1);
    void          clear();
    GVector       extract_row(int row) const;
    GVector       extract_col(int col) const;
    //TBD GSparseMatrix extract_lower_triangle() const;
    //TBD GSparseMatrix extract_upper_triangle() const;
    double        fill() const;
    void          insert_col(const GVector& v, int col);
    void          insert_col(const double* values, const int* rows, int number, int col);
    double        min() const;
    double        max() const;
    void          set_mem_block(int block);
    void          stack_init(int size = 0, int entries = 0);
    int           stack_push_column(const GVector& v, int col);
    int           stack_push_column(const double* values, const int* rows,
                                    int number, int col);
    void          stack_flush(void);
    void          stack_destroy(void);
    double        sum() const;
    void          transpose();
};


/***************************************************************************
 *                  GSparseMatrix class SWIG extension                     *
 ***************************************************************************/
%extend GSparseMatrix {
    char *__str__() {
        static char str_buffer[100001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 100001);
        str_buffer[100000] = '\0';
        return str_buffer;
    }
    double __getslice__(int row, int col) {
        if (row >=0 && row < (int)self->rows() && col >= 0 || col < (int)self->cols())
            return (*self)(row,col);
        else
            throw GException::out_of_range("__getitem__(int,int)", row, col,
                                           (int)self->rows(), (int)self->cols());
    }
    void __setslice__(int row, int col, const double val) {
        if (row >=0 && row < (int)self->rows() && col >= 0 || col < (int)self->cols())
            (*self)(row,col) = val;
        else
            throw GException::out_of_range("__setitem__(int,int)", row, col,
                                           (int)self->rows(), (int)self->cols());
    }
    GSparseMatrix __mul__(const double &a) {
        return (*self) * a;
    }
    GSparseMatrix __div__(const double &a) {
        return (*self) / a;
    }
    GSparseMatrix copy() {
        return (*self);
    }
};
