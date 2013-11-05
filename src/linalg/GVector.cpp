/***************************************************************************
 *                        GVector.cpp - Vector class                       *
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
 * @file GVector.cpp
 * @brief Vector class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GVector.hpp"
#include "GTools.hpp"


/* __ Method name definitions ____________________________________________ */
#define G_OP_ADD                              "GVector::operator+=(GVector&)"
#define G_OP_SUB                              "GVector::operator-=(GVector&)"
#define G_AT                                              "GVector::at(int&)"
#define G_CROSS                                   "cross(GVector&, GVector&)"
#define G_SCALAR                              "operator*(GVector&, GVector&)"


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void vector constructor
 ***************************************************************************/
GVector::GVector(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Vector constructor
 *
 * @param[in] num Number of elements in vector.
 *
 * Initialises a vector with @p num elements. All vector elements will be
 * set to 0.
 ***************************************************************************/
GVector::GVector(const int& num)
{
    // Initialise class members
    init_members();

    // Store vector size
    m_num = num;

    // Allocate vector (filled with 0)
    alloc_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Single element vector constructor
 *
 * @param[in] a Vector element.
 *
 * Initialises a vector with a single element.
 ***************************************************************************/
GVector::GVector(const double& a)
{
    // Initialise class members
    init_members();

    // Store vector size
    m_num = 1;

    // Allocate vector
    alloc_members();

    // Set value
    m_data[0] = a;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Two elements vector constructor
 *
 * @param[in] a First vector element.
 * @param[in] b Second vector element.
 *
 * Initialises a vector with two elements.
 ***************************************************************************/
GVector::GVector(const double& a, const double& b)
{
    // Initialise class members
    init_members();

    // Store vector size
    m_num = 2;

    // Allocate vector
    alloc_members();

    // Set values
    m_data[0] = a;
    m_data[1] = b;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Three elements vector constructor
 *
 * @param[in] a First vector element.
 * @param[in] b Second vector element.
 * @param[in] c Third vector element.
 *
 * Initialises a vector with three elements.
 ***************************************************************************/
GVector::GVector(const double& a, const double& b, const double& c)
{
    // Initialise class members
    init_members();

    // Store vector size
    m_num = 3;

    // Allocate vector
    alloc_members();

    // Set values
    m_data[0] = a;
    m_data[1] = b;
    m_data[2] = c;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] vector Vector.
 ***************************************************************************/
GVector::GVector(const GVector& vector)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(vector);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GVector::~GVector(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] vector Vector.
 * @return Vector.
 ***************************************************************************/
GVector& GVector::operator=(const GVector& vector)
{
    // Execute only if object is not identical
    if (this != &vector) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(vector);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Equality operator
 *
 * @param[in] vector Vector.
 * @return True if vectors are identical.
 *
 * Returns true if vectors are identical. Vectors are considered identical
 * if they have the same size and if all their elements are identical.
 ***************************************************************************/
bool GVector::operator==(const GVector& vector) const
{
    // Initalise result depending on vector size identity
    bool result = (m_num == vector.m_num);

    // Test for difference. Break at first difference
    if (result) {
        for (int i = 0; i < m_num; ++i) {
            if (m_data[i] != vector.m_data[i]) {
                result = false;
                break;
            }
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Non-equality operator
 *
 * @param[in] vector Vector.
 * @return True if both vectors are different.
 ***************************************************************************/
bool GVector::operator!=(const GVector& vector) const
{
    // Get negated result of equality operation
    bool result = !(this->operator==(vector));
	
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Unary addition operator
 *
 * @param[in] vector Vector.
 * @return Vector.
 *
 * @exception GException::vector_mismatch
 *            Vectors have not the same size.
 *
 * Adds vector.
 ***************************************************************************/
GVector& GVector::operator+=(const GVector& vector)
{
    // Raise exception if vectors mismatch
    if (m_num != vector.m_num) {
        throw GException::vector_mismatch(G_OP_ADD, m_num, vector.m_num);
    }

    // Add elements
    for (int i = 0; i < m_num; ++i) {
        m_data[i] += vector.m_data[i];
    }

    // Return vector
    return *this;
}


/***********************************************************************//**
 * @brief Unary subtraction operator
 *
 * @param[in] vector Vector.
 * @return Vector.
 *
 * @exception GException::vector_mismatch
 *            Vectors have not the same size.
 *
 * Subtracts vector.
 ***************************************************************************/
GVector& GVector::operator-=(const GVector& vector)
{
    // Raise exception if vectors mismatch
    if (m_num != vector.m_num) {
        throw GException::vector_mismatch(G_OP_SUB, m_num, vector.m_num);
    }

    // Subtract elements
    for (int i = 0; i < m_num; ++i) {
        m_data[i] -= vector.m_data[i];
    }

    // Return vector
    return *this;
}


/***********************************************************************//**
 * @brief Scalar assignment operator
 *
 * @param[in] scalar Scalar.
 * @return Vector.
 *
 * Subtracts vector.
 ***************************************************************************/
GVector& GVector::operator=(const double& scalar)
{
    // Set elements
    for (int i = 0; i < m_num; ++i) {
        m_data[i] = scalar;
    }

    // Return vector
    return *this;
}


/***********************************************************************//**
 * @brief Scalar unary addition operator
 *
 * @param[in] scalar Scalar.
 * @return Vector.
 *
 * Adds scalar to all vector elements.
 ***************************************************************************/
GVector& GVector::operator+=(const double& scalar)
{
    // Add scalar to elements
    for (int i = 0; i < m_num; ++i) {
        m_data[i] += scalar;
    }

    // Return vector
    return *this;
}


/***********************************************************************//**
 * @brief Scalar unary subtraction operator
 *
 * @param[in] scalar Scalar.
 * @return Vector.
 *
 * Subtract scalar to all vector elements.
 ***************************************************************************/
GVector& GVector::operator-=(const double& scalar)
{
    // Subtract scalar from elements
    for (int i = 0; i < m_num; ++i) {
        m_data[i] -= scalar;
    }

    // Return vector
    return *this;
}


/***********************************************************************//**
 * @brief Scalar unary multiplication operator
 *
 * @param[in] scalar Scalar.
 * @return Vector.
 *
 * Multiply all vector elements by scalar.
 ***************************************************************************/
GVector& GVector::operator*=(const double& scalar)
{
    // Multiply all elements
    for (int i = 0; i < m_num; ++i) {
        m_data[i] *= scalar;
    }

    // Return vector
    return *this;
}


/***********************************************************************//**
 * @brief Scalar unary division operator
 *
 * @param[in] scalar Scalar.
 * @return Vector.
 *
 * Divide all vector elements by scalar.
 ***************************************************************************/
GVector& GVector::operator/=(const double& scalar)
{
    // Divide all elements
    for (int i = 0; i < m_num; ++i) {
        m_data[i] /= scalar;
    }

    // Return vector
    return *this;
}


/***********************************************************************//**
 * @brief Unary minus operator
 *
 * @return Vector.
 *
 * Negate all vector elements.
 ***************************************************************************/
GVector GVector::operator-(void) const
{
    // Copy vector
    GVector result = *this;
    
    // Negate all elements
    for (int i = 0; i < m_num; ++i) {
        result.m_data[i] = -result.m_data[i];
    }

    // Return vector
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear vector
 ***************************************************************************/
void GVector::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone vector
 *
 * @return Pointer to deep copy of vector.
 ***************************************************************************/
GVector* GVector::clone(void) const
{
    // Clone vector
    return new GVector(*this);
}


/***********************************************************************//**
 * @brief Vector element access with range checking
 *
 * @param[in] index Element index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Element index is out of range.
 ***************************************************************************/
double& GVector::at(const int& index)
{
    // Raise an exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, size()-1);
    }

    // Return vector element
    return m_data[index];
}


/***********************************************************************//**
 * @brief Vector element access with range checking (const variant)
 *
 * @param[in] index Element index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Element index is out of range.
 ***************************************************************************/
const double& GVector::at(const int& index) const
{
    // Raise an exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, size()-1);
    }

    // Return vector element
    return m_data[index];
}


/***********************************************************************//**
 * @brief Returns number of non-zero elements in vector
 *
 * @return Number of non-zero elements in vector.
 ***************************************************************************/
int GVector::non_zeros(void) const
{
    // Initialise number of non-zeros
    int non_zeros = 0;

    // Gather all non-zero elements
    for (int i = 0; i < m_num; ++i) {
        if (m_data[i] != 0.0) {
            non_zeros++;
        }
    }

    // Return number of non-zero elements
    return non_zeros;
}


/***********************************************************************//**
 * @brief Print vector information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing vector information.
 ***************************************************************************/
std::string GVector::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result = "(";

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Put all elements in stream
        for (int i = 0; i < m_num; ++i) {
            result += gammalib::str((*this)[i]);
            if (i != m_num-1) {
                result += ", ";
            }
        }

        // Append )
        result += ")";

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GVector::init_members(void)
{
    // Initialise members
    m_num  = 0;
    m_data = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate vector
 ***************************************************************************/
void GVector::alloc_members(void)
{
    // Continue only if vector has non-zero length
    if (m_num > 0) {

        // Allocate vector and initialize elements to 0
        m_data = new double[m_num];
        for (int i = 0; i < m_num; ++i) {
            m_data[i] = 0.0;
        }

    } // endif: vector had non-zero length

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] vector Vector from which members should be copied.
 ***************************************************************************/
void GVector::copy_members(const GVector& vector)
{
    // Copy attributes
    m_num = vector.m_num;

    // Copy elements
    if (m_num > 0) {
        alloc_members();
        for (int i = 0; i <  m_num; ++i) {
            m_data[i] = vector.m_data[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GVector::free_members(void)
{
    // Free memory
    if (m_data != NULL) delete[] m_data;

    // Signal free pointers
    m_data = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Vector cross product
 *
 * @param[in] a Vector.
 * @param[in] b Vector.
 * @return Vector cross product.
 *
 * @exception GException::vector_mismatch
 *            Mismatch between vector size.
 *
 * Computes the cross product between two 3-element vectors (note that the
 * cross product is only defined for 3-element vectors).
 ***************************************************************************/
GVector cross(const GVector& a, const GVector& b)
{
    // Verify that vectors have same dimensions
    if (a.m_num != b.m_num) {
        throw GException::vector_mismatch(G_CROSS, a.m_num, b.m_num);
    }

    // Verify that vectors have 3 elements
    if (a.m_num != 3) {
       throw GException::vector_bad_cross_dim(G_CROSS, a.m_num);
    }

    // Compute cross product
    GVector result(3);
    result.m_data[0] = a.m_data[1]*b.m_data[2] - a.m_data[2]*b.m_data[1];
    result.m_data[1] = a.m_data[2]*b.m_data[0] - a.m_data[0]*b.m_data[2];
    result.m_data[2] = a.m_data[0]*b.m_data[1] - a.m_data[1]*b.m_data[0];

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Vector scalar product
 *
 * @param[in] a Vector.
 * @param[in] b Vector.
 * @return Product between vector @p a and @p b.
 *
 * @exception GException::vector_mismatch
 *            Mismatch between vector size.
 *
 * Returns the scalar product between vector @p a and @p b.
 ***************************************************************************/
double operator*(const GVector& a, const GVector& b)
{
    // Verify that vectors have same dimensions
    if (a.m_num != b.m_num) {
        throw GException::vector_mismatch(G_SCALAR, a.m_num, b.m_num);
    }

    // Compute scalar product
    double result = 0.0;
    for (int i = 0; i < a.m_num; ++i) {
        result += (a.m_data[i] * b.m_data[i]);
    }

    // Return scalar product
    return result;
}


/***********************************************************************//**
 * @brief Computes vector norm
 *
 * @param[in] vector Vector.
 * @return Vector norm.
 ***************************************************************************/
double norm(const GVector& vector)
{
    // Initialises result
    double result = 0.0;

    // Computes norm
    for (int i = 0; i < vector.m_num; ++i) {
        result += (vector.m_data[i] * vector.m_data[i]);
    }
    result = (result > 0.0) ? std::sqrt(result) : 0.0;

    // Returns norm
    return result;
}


/***********************************************************************//**
 * @brief Computes minimum vector element
 *
 * @param[in] vector Vector.
 * @return Minimum vector element.
 ***************************************************************************/
double min(const GVector& vector)
{
    // Initialises result
    double result = 0.0;

    // Continue only if we have elements
    if (vector.m_num > 0) {
    
        // Search for minimum
        result = vector.m_data[0];
        for (int i = 1; i < vector.m_num; ++i) {
            if (vector.m_data[i] < result) {
                result = vector.m_data[i];
            }
        }

    } // endif: there were elements

    // Returns minimum
    return result;
}


/***********************************************************************//**
 * @brief Computes maximum vector element
 *
 * @param[in] vector Vector.
 * @return Maximum vector element.
 ***************************************************************************/
double max(const GVector& vector)
{
    // Initialises result
    double result = 0.0;

    // Continue only if we have elements
    if (vector.m_num > 0) {
    
        // Search for maximum
        result = vector.m_data[0];
        for (int i = 1; i < vector.m_num; ++i) {
            if (vector.m_data[i] > result) {
                result = vector.m_data[i];
            }
        }

    } // endif: there were elements

    // Returns maximum
    return result;
}


/***********************************************************************//**
 * @brief Computes vector sum
 *
 * @param[in] vector Vector.
 * @return Sum of vector elements.
 ***************************************************************************/
double sum(const GVector& vector)
{
    // Compute sum
    double result = 0.0;
    for (int i = 0; i < vector.m_num; ++i) {
        result += vector.m_data[i];
    }

    // Returns sum
    return result;
}


/***********************************************************************//**
 * @brief Computes vector permutation
 *
 * @param[in] vector Vector.
 * @param[in] p Permutation array.
 * @return Permuted vector.
 ***************************************************************************/
GVector perm(const GVector& vector, const int* p)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Compute permutations
    if (p != NULL) {
        for (int i = 0; i < vector.m_num; ++i) {
            result.m_data[i] = vector.m_data[p[i]];
        }
    }
    else {
        result = vector;
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes vector inverse permutation
 *
 * @param[in] vector Vector.
 * @param[in] p Permutation array.
 * @return Inversely permuted vector.
 ***************************************************************************/
GVector iperm(const GVector& vector, const int* p)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Compute permutations
    if (p != NULL) {
        for (int i = 0; i < vector.m_num; ++i) {
            result.m_data[p[i]] = vector.m_data[i];
        }
    }
    else {
        result = vector;
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes arccos of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the arccos of every element.
 ***************************************************************************/
GVector acos(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::acos(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes acosh of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the acosh of every element.
 ***************************************************************************/
GVector acosh(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = acosh(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes arcsin of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the arcsin of every element.
 ***************************************************************************/
GVector asin(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::asin(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes asinh of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the asinh of every element.
 ***************************************************************************/
GVector asinh(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = asinh(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes arctan of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the arctan of every element.
 ***************************************************************************/
GVector atan(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::atan(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes atanh of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the atanh of every element.
 ***************************************************************************/
GVector atanh(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = atanh(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes cosine of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the cosine of every element.
 ***************************************************************************/
GVector cos(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::cos(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes cosh of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the cosh of every element.
 ***************************************************************************/
GVector cosh(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::cosh(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes exponential of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the exponential of every element.
 ***************************************************************************/
GVector exp(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::exp(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes absolute of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the absolute of every element.
 ***************************************************************************/
GVector abs(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::abs(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes natural logarithm of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the natural logarithm of every element.
 ***************************************************************************/
GVector log(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::log(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes base10 logarithm of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the base10 logarithm of every element.
 ***************************************************************************/
GVector log10(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::log10(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes sine of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the sine of every element.
 ***************************************************************************/
GVector sin(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::sin(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes sinh of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the sinh of every element.
 ***************************************************************************/
GVector sinh(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::sinh(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes square root of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the square root of every element.
 ***************************************************************************/
GVector sqrt(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::sqrt(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes tangens of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the tangens of every element.
 ***************************************************************************/
GVector tan(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::tan(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes tanh of vector elements
 *
 * @param[in] vector Vector.
 * @return Vector containing the tanh of every element.
 ***************************************************************************/
GVector tanh(const GVector& vector)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::tanh(vector.m_data[i]);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Computes tanh of vector elements
 *
 * @param[in] vector Vector.
 * @param[in] power Power.
 * @return Vector containing the power of every element.
 ***************************************************************************/
GVector pow(const GVector& vector, const double& power)
{
    // Initialise result vector
    GVector result(vector.m_num);

    // Evaluate each vector element
    for (int i = 0; i < vector.m_num; ++i) {
        result.m_data[i] = std::pow(vector.m_data[i], power);
    }

    // Return vector
    return result;
}
