/***************************************************************************
 *      GSparseSymbolic.hpp - Sparse matrix symbolic analysis class        *
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
 * @file GSparseSymbolic.hpp
 * @brief Sparse matrix symbolic analysis class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPARSESYMBOLIC_HPP
#define GSPARSESYMBOLIC_HPP

/* __ Includes ___________________________________________________________ */

/* __ Definitions ________________________________________________________ */

/* __ Enumerators ________________________________________________________ */

/* __ Structures _________________________________________________________ */

/* __ Prototypes (implementions in GMatrixSparse) ________________________ */
GMatrixSparse cs_symperm(const GMatrixSparse& matrix, const int* pinv);
GMatrixSparse cs_transpose(const GMatrixSparse& matrix, int values);
double        cs_cumsum(int* p, int* c, int n);


/***********************************************************************//**
 * @class GSparseSymbolic
 *
 * @brief Sparse matrix symbolic analysis class
 *
 * This class implements the symbolic analysis of a sparse matrix.
 ***************************************************************************/
class GSparseSymbolic {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GSparseSymbolic& s);

    // Friend classes
    friend class GMatrixSparse;
    friend class GSparseNumeric;

public:
    // Constructors and destructors
    GSparseSymbolic(void);
    virtual ~GSparseSymbolic(void);

    // Assignment operator
    GSparseSymbolic& operator= (const GSparseSymbolic& s);

    // Methods
    void cholesky_symbolic_analysis(int order, const GMatrixSparse& m);

private:
    // Private methods
    int*        cs_amd (int order, const GMatrixSparse* A);
    int*        cs_counts(const GMatrixSparse* A, const int* parent, const int* post, int ata);
    int*        cs_etree(const GMatrixSparse* A, int ata);
    int         cs_fkeep(GMatrixSparse* A, int(*fkeep)(int, int, double, void*), void* other);
    int         cs_leaf(int i, int j, const int* first, int* maxfirst, int* prevleaf, int* ancestor, int* jleaf);
    int*        cs_pinv(int const* p, int n);
    int*        cs_post(const int* parent, int n);
    int         cs_tdfs(int j, int k, int* head, const int* next, int* post, int* stack);
    static void init_ata(const GMatrixSparse* AT, const int* post, int* wrk_int, int** head, int** next);
    static int  cs_diag(int i, int j, double aij, void* other);
    static int  cs_wclear(int mark, int lemax, int* w, int n);

    // Data
    int*   m_pinv;        //!< Inverse row permutation for QR, fill reduce permutation for Cholesky
    int*   m_q;           //!< Fill-reducing column permutation for LU and QR
    int*   m_parent;      //!< elimination tree for Cholesky and QR
    int*   m_cp;	      //!< column pointers for Cholesky, row counts for QR
    int*   m_leftmost;    //!< leftmost[i] = min(find(A(i,:))), for QR
    int    m_m2;	      //!< # of rows for QR, after adding fictitious rows
    double m_lnz;         //!< # entries in L for LU or Cholesky; in V for QR
    double m_unz;         //!< # entries in U for LU; in R for QR
    int    m_n_pinv;      //!< Number of elements in m_pinv
    int    m_n_q;         //!< Number of elements in m_q
    int    m_n_parent;    //!< Number of elements in m_parent
    int    m_n_cp;        //!< Number of elements in m_cp
    int    m_n_leftmost;  //!< Number of elements in m_leftmost
};

#endif /* GSPARSESYMBOLIC_HPP */
