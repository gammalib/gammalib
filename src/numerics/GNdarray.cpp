/***************************************************************************
 *                 GNdarray.cpp - N-dimensional array class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2017 by Juergen Knoedlseder                         *
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
 * @file GNdarray.cpp
 * @brief N-dimensional array class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GNdarray.hpp"


/* __ Method name definitions ____________________________________________ */
#define G_OP_ADD                            "GNdarray::operator+=(GNdarray&)"
#define G_OP_SUB                            "GNdarray::operator-=(GNdarray&)"
#define G_SHAPE                          "GNdarray::shape(std::vector<int>&)"
#define G_AT1                                            "GNdarray::at(int&)"
#define G_AT2                                      "GNdarray::at(int&, int&)"
#define G_AT3                                "GNdarray::at(int&, int&, int&)"
#define G_ATN                               "GNdarray::at(std::vector<int>&)"


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GNdarray::GNdarray(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief 1-dimensional constructor
 *
 * @param[in] nx Array dimension.
 *
 * Constructs a 1-dimensional array.
 ***************************************************************************/
GNdarray::GNdarray(const int& nx)
{
    // Initialise class members
    init_members();

    // Set shape
    m_shape.push_back(nx);

    // Set strides
    m_strides.push_back(1);

    // Initialise data array with zeros
    m_data.assign(nx, 0.0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief 2-dimensional constructor
 *
 * @param[in] nx Array's first dimension.
 * @param[in] ny Array's second dimension.
 *
 * Constructs a 2-dimensional array.
 ***************************************************************************/
GNdarray::GNdarray(const int& nx, const int& ny)
{
    // Initialise class members
    init_members();

    // Compute size of array
    int size = nx * ny;

    // Set shape
    m_shape.push_back(nx);
    m_shape.push_back(ny);

    // Set strides
    m_strides.push_back(1);
    m_strides.push_back(nx);

    // Initialise data array with zeros
    m_data.assign(size, 0.0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief 3-dimensional constructor
 *
 * @param[in] nx Array's first dimension.
 * @param[in] ny Array's second dimension.
 * @param[in] nz Array's third dimension.
 *
 * Constructs a 3-dimensional array.
 ***************************************************************************/
GNdarray::GNdarray(const int& nx, const int& ny, const int& nz)
{
    // Initialise class members
    init_members();

    // Compute size of array
    int size = nx * ny * nz;

    // Set shape
    m_shape.push_back(nx);
    m_shape.push_back(ny);
    m_shape.push_back(nz);

    // Set strides
    m_strides.push_back(1);
    m_strides.push_back(nx);
    m_strides.push_back(nx*ny);

    // Initialise data array with zeros
    m_data.assign(size, 0.0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief n-dimensional constructor
 *
 * @param[in] n Array dimensions.
 *
 * Constructs a n-dimensional array. If the array dimension vector @p n is
 * empty an empty n-dimensional array is constructed (similar to the void
 * constructor).
 ***************************************************************************/
GNdarray::GNdarray(const std::vector<int>& n)
{
    // Initialise class members
    init_members();

    // Continue only if there are dimensions in vector
    if (n.size() > 0) {

        // Compute size of array
        int size = n[0];
        for (int i = 1; i < n.size(); ++i) {
            size *= n[i];
        }

        // Set shape
        m_shape = n;

        // Set strides
        int stride = 1;
        for (int i = 0; i < n.size(); ++i) {
            m_strides.push_back(stride);
            stride *= n[i];
        }

        // Initialise data array with zeros
        m_data.assign(size, 0.0);

    } // endif: there were dimensions in vector

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] array Array.
 ***************************************************************************/
GNdarray::GNdarray(const GNdarray& array)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(array);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GNdarray::~GNdarray(void)
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
 * @param[in] array Array.
 * @return Array.
 ***************************************************************************/
GNdarray& GNdarray::operator=(const GNdarray& array)
{
    // Execute only if object is not identical
    if (this != &array) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(array);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Equality operator
 *
 * @param[in] array Array.
 * @return True if arrays are identical.
 *
 * Returns true if the arrays have identical shape and content.
 ***************************************************************************/
bool GNdarray::operator==(const GNdarray& array) const
{
    // Check of arrays have the same shape
    bool identity = has_same_shape(array);

    // If the arrays have identical dimensions then check if all elements
    // are identical. Break as soon as one array element differs.
    if (identity) {
        for (int i = 0; i < m_data.size(); ++i) {
            if (m_data[i] != array.m_data[i]) {
                identity = false;
                break;
            }
        }
    }

    // Return result
    return identity;
}


/***********************************************************************//**
 * @brief Non-equality operator
 *
 * @param[in] array Array.
 * @return True if arrays are different.
 *
 * Returns true if the arrays have different shape or content.
 ***************************************************************************/
bool GNdarray::operator!=(const GNdarray& array) const
{
    // Get negated result of equality operation
    bool difference = !(this->operator==(array));
	
    // Return result
    return difference;
}


/***********************************************************************//**
 * @brief Unary addition operator
 *
 * @param[in] array Array.
 * @return Array.
 *
 * Adds an array to the current array.
 ***************************************************************************/
GNdarray& GNdarray::operator+=(const GNdarray& array)
{
    // Throw an exception if the arrays have not the same shape
    require_same_shape(G_OP_ADD, array);

    // Add elements
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i] += array.m_data[i];
    }

    // Return array
    return *this;
}


/***********************************************************************//**
 * @brief Unary subtraction operator
 *
 * @param[in] array Array.
 * @return Array.
 *
 * Subtracts an array from the current array.
 ***************************************************************************/
GNdarray& GNdarray::operator-=(const GNdarray& array)
{
    // Throw an exception if the arrays have not the same shape
    require_same_shape(G_OP_SUB, array);

    // Subtract elements
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i] -= array.m_data[i];
    }

    // Return array
    return *this;
}


/***********************************************************************//**
 * @brief Unary value addition operator
 *
 * @param[in] value Value.
 * @return Array.
 *
 * Adds a value to all array elements.
 ***************************************************************************/
GNdarray& GNdarray::operator+=(const double& value)
{
    // Add value to elements
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i] += value;
    }

    // Return array
    return *this;
}


/***********************************************************************//**
 * @brief Unary value subtraction operator
 *
 * @param[in] value Value.
 * @return Array.
 *
 * Subtracts a value from all array elements.
 ***************************************************************************/
GNdarray& GNdarray::operator-=(const double& value)
{
    // Subtract value from elements
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i] -= value;
    }

    // Return array
    return *this;
}


/***********************************************************************//**
 * @brief Unary multiplication operator
 *
 * @param[in] value Value.
 * @return Array.
 *
 * Multiply all array elements by a value.
 ***************************************************************************/
GNdarray& GNdarray::operator*=(const double& value)
{
    // Multiply all elements
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i] *= value;
    }

    // Return array
    return *this;
}


/***********************************************************************//**
 * @brief Unary division operator
 *
 * @param[in] value Value.
 * @return Array.
 *
 * Divide all array elements by a value.
 ***************************************************************************/
GNdarray& GNdarray::operator/=(const double& value)
{
    // Divide all elements
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i] /= value;
    }

    // Return array
    return *this;
}


/***********************************************************************//**
 * @brief Unary minus operator
 *
 * @return Array.
 *
 * Negate all array elements.
 ***************************************************************************/
GNdarray GNdarray::operator-(void) const
{
    // Copy array
    GNdarray result = *this;
    
    // Negate all elements
    for (int i = 0; i < m_data.size(); ++i) {
        result.m_data[i] = -result.m_data[i];
    }

    // Return array
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear array
 ***************************************************************************/
void GNdarray::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone array
 *
 * @return Pointer to deep copy of array.
 ***************************************************************************/
GNdarray* GNdarray::clone(void) const
{
    // Clone array
    return new GNdarray(*this);
}


/***********************************************************************//**
 * @brief Set shape of array
 *
 * @param[in] shape Shape vector.
 *
 * @exception GException::invalid_argument
 *            Invalid shape factorisation specified.
 *
 * Set the shape of the array. The shape specifies how the information is
 * arranged in a n-dimensional array.
 ***************************************************************************/
void GNdarray::shape(const std::vector<int>& shape)
{
    // Computes the number of array elements
    int nelements = 0;
    if (shape.size() > 0) {
        nelements = 1;
        for (int i = 0; i < shape.size(); ++i) {
            nelements *= shape[i];
        }
    }

    // Throw an exception if resulting number of elements is not equal to
    // the existing number of elements
    if (nelements != size()) {
        std::string msg = "Number of elements "+gammalib::str(nelements)+
                          " in specified shape is not identical to the "
                          "number of elements "+gammalib::str(size())+
                          " in the array.";
        throw GException::invalid_argument(G_SHAPE, msg);
    }

    // Set shape
    m_shape = shape;

    // Return
    return;
}


/***********************************************************************//**
 * @brief 1-dimensional array element access with range checking
 *        (const version)
 *
 * @param[in] ix Element index [0,...,shape(0)-1].
 * @return Const reference to array element.
 *
 * @exception GException::invalid_value
 *            Array is not 1-dimensional.
 * @exception GException::out_of_range
 *            Element index is out of range.
 ***************************************************************************/
const double& GNdarray::at(const int& ix) const
{
    // Throw an exception if the array is not 1-dimensional
    if (m_shape.size() != 1) {
        std::string msg = "Invalid access of "+
                          gammalib::str(m_shape.size())+"-dimensional array "
                          "with 1-dimensional access operator.";
        throw GException::invalid_value(G_AT1, msg);
    }

    // Throw an exception if the index is outside the valid range
    if (ix < 0 || ix >= m_shape[0]) {
        throw GException::out_of_range(G_AT1, "Array index", ix, m_shape[0]);
    }

    // Return array element
    return (*this)(ix);
}


/***********************************************************************//**
 * @brief 2-dimensional array element access with range checking
 *        (const version)
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @return Const reference to array element.
 *
 * @exception GException::invalid_value
 *            Array is not 2-dimensional.
 * @exception GException::out_of_range
 *            Element index is out of range.
 ***************************************************************************/
const double& GNdarray::at(const int& ix, const int& iy) const
{
    // Throw an exception if the array is not 2-dimensional
    if (m_shape.size() != 2) {
        std::string msg = "Invalid access of "+
                          gammalib::str(m_shape.size())+"-dimensional array "
                          "with 2-dimensional access operator.";
        throw GException::invalid_value(G_AT2, msg);
    }

    // Throw an exception if the indices are outside the valid range
    if (ix < 0 || ix >= m_shape[0]) {
        throw GException::out_of_range(G_AT2, "First array index", ix,
                                       m_shape[0]);
    }
    if (iy < 0 || iy >= m_shape[1]) {
        throw GException::out_of_range(G_AT2, "Second array index", iy,
                                       m_shape[1]);
    }

    // Return array element
    return (*this)(ix,iy);
}


/***********************************************************************//**
 * @brief 3-dimensional array element access with range checking
 *        (const version)
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @param[in] iz Index in third dimension [0,...,shape(2)-1].
 * @return Const reference to array element.
 *
 * @exception GException::invalid_value
 *            Array is not 3-dimensional.
 * @exception GException::out_of_range
 *            Element index is out of range.
 ***************************************************************************/
const double& GNdarray::at(const int& ix, const int& iy, const int& iz) const
{
    // Throw an exception if the array is not 2-dimensional
    if (m_shape.size() != 3) {
        std::string msg = "Invalid access of "+
                          gammalib::str(m_shape.size())+"-dimensional array "
                          "with 3-dimensional access operator.";
        throw GException::invalid_value(G_AT3, msg);
    }

    // Throw an exception if the indices are outside the valid range
    if (ix < 0 || ix >= m_shape[0]) {
        throw GException::out_of_range(G_AT3, "First array index", ix,
                                       m_shape[0]);
    }
    if (iy < 0 || iy >= m_shape[1]) {
        throw GException::out_of_range(G_AT3, "Second array index", iy,
                                       m_shape[1]);
    }
    if (iz < 0 || iz >= m_shape[2]) {
        throw GException::out_of_range(G_AT3, "Third array index", iz,
                                       m_shape[2]);
    }

    // Return array element
    return (*this)(ix,iy,iz);
}


/***********************************************************************//**
 * @brief n-dimensional array element access with range checking
 *        (const version)
 *
 * @param[in] i Index vector.
 * @return Const reference to array element.
 *
 * @exception GException::invalid_value
 *            Array is not of correct dimension.
 * @exception GException::out_of_range
 *            Element index is out of range.
 ***************************************************************************/
const double& GNdarray::at(const std::vector<int>& i) const
{
    // Throw an exception if the array is not 2-dimensional
    if (m_shape.size() != i.size()) {
        std::string msg = "Invalid access of "+
                          gammalib::str(m_shape.size())+"-dimensional array "
                          "with "+
                          gammalib::str(i.size())+"-dimensional access "
                          "operator.";
        throw GException::invalid_value(G_ATN, msg);
    }

    // Throw an exception if the indices are outside the valid range
    for (int k = 0; k < m_shape.size(); ++k) {
        if (i[k] < 0 || i[k] >= m_shape[k]) {
            throw GException::out_of_range(G_ATN, "Dimension "+
                                           gammalib::str(k)+" index",
                                           i[k], m_shape[k]);
        }
    }

    // Return array element
    return (*this)(i);
}


/***********************************************************************//**
 * @brief Print array information
 *
 * @param[in] chatter Chattiness.
 * @return String containing array information.
 ***************************************************************************/
std::string GNdarray::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result = "(";

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Put all elements in string
        for (int i = 0; i < m_data.size(); ++i) {
            if (i > 0) {
                result += ", ";
            }
            result += gammalib::str(m_data[i]);
        }

        // Append closing parentheses
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
void GNdarray::init_members(void)
{
    // Initialise members
    m_shape.clear();
    m_strides.clear();
    m_data.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] array Array.
 ***************************************************************************/
void GNdarray::copy_members(const GNdarray& array)
{
    // Copy members
    m_shape   = array.m_shape;
    m_strides = array.m_strides;
    m_data    = array.m_data;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GNdarray::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute array element index
 *
 * @param[in] i Index vector.
 * @return Array element index.
 *
 * Computes the array element index from an index vector. The method assumes
 * that the index vector has the same dimension as the array. If this
 * condition is not guaranteed, the condition needs to be checked before
 * calling the method.
 ***************************************************************************/
int GNdarray::index(const std::vector<int>& i) const
{
    // Initialise element index
    int index = 0;

    // Compute element index
    for (int k = 0; k < m_shape.size(); ++k) {
        index += m_strides[k] * i[k];
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Check if array has the same shape
 *
 * @param[in] array Array.
 * @return True if the array has the same shape.
 ***************************************************************************/
bool GNdarray::has_same_shape(const GNdarray& array) const
{
    // Check if the array has the same dimension
    bool identity = m_shape.size() == array.m_shape.size();

    // Check if the shape is identical. Break as soon as one shape value
    // differs.
    if (identity) {
        for (int i = 0; i < m_shape.size(); ++i) {
            if (m_shape[i] != array.m_shape[i]) {
                identity = false;
                break;
            }
        }
    }

    // Return result
    return identity;
}


/***********************************************************************//**
 * @brief Throw exception if array shapes differ
 *
 * @param[in] method Method that throws exception.
 * @param[in] array Array.
 ***************************************************************************/
void GNdarray::require_same_shape(const std::string& method,
                                  const GNdarray&    array) const
{
    // If the shape differs then throw an exception
    if (!has_same_shape(array)) {

        // Compose exception message
        std::string msg = "Incompatible array dimensions (";
        for (int i = 0; i < m_shape.size(); ++i) {
            if (i > 0) {
                msg += ", ";
            }
            msg += gammalib::str(m_shape[i]);
        }
        msg += ") <=> (";
        for (int i = 0; i < array.m_shape.size(); ++i) {
            if (i > 0) {
                msg += ", ";
            }
            msg += gammalib::str(array.m_shape[i]);
        }
        msg += ").";

        // Throw exception
        throw GException::invalid_argument(method, msg);

    } // endif: array shapes differed

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Computes minimum array element
 *
 * @param[in] array Array.
 * @return Minimum array element.
 ***************************************************************************/
double min(const GNdarray& array)
{
    // Initialises result
    double result = 0.0;

    // Continue only if we have elements
    if (array.m_data.size() > 0) {
    
        // Search for minimum
        result = array.m_data[0];
        for (int i = 1; i < array.m_data.size(); ++i) {
            if (array.m_data[i] < result) {
                result = array.m_data[i];
            }
        }

    } // endif: there were elements

    // Returns minimum
    return result;
}


/***********************************************************************//**
 * @brief Computes maximum array element
 *
 * @param[in] array Array.
 * @return Maximum array element.
 ***************************************************************************/
double max(const GNdarray& array)
{
    // Initialises result
    double result = 0.0;

    // Continue only if we have elements
    if (array.m_data.size() > 0) {
    
        // Search for maximum
        result = array.m_data[0];
        for (int i = 1; i < array.m_data.size(); ++i) {
            if (array.m_data[i] > result) {
                result = array.m_data[i];
            }
        }

    } // endif: there were elements

    // Returns maximum
    return result;
}


/***********************************************************************//**
 * @brief Computes array sum
 *
 * @param[in] array Array.
 * @return Sum of array elements.
 ***************************************************************************/
double sum(const GNdarray& array)
{
    // Compute sum
    double result = 0.0;
    for (int i = 0; i < array.m_data.size(); ++i) {
        result += array.m_data[i];
    }

    // Returns sum
    return result;
}


/***********************************************************************//**
 * @brief Computes arccos of array elements
 *
 * @param[in] array Array.
 * @return Array containing the arccos of every element.
 ***************************************************************************/
GNdarray acos(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::acos(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes acosh of array elements
 *
 * @param[in] array Array.
 * @return Array containing the acosh of every element.
 ***************************************************************************/
GNdarray acosh(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = acosh(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes arcsin of array elements
 *
 * @param[in] array Array.
 * @return Array containing the arcsin of every element.
 ***************************************************************************/
GNdarray asin(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::asin(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes asinh of array elements
 *
 * @param[in] array Array.
 * @return Array containing the asinh of every element.
 ***************************************************************************/
GNdarray asinh(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = asinh(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes arctan of array elements
 *
 * @param[in] array Array.
 * @return Array containing the arctan of every element.
 ***************************************************************************/
GNdarray atan(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::atan(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes atanh of array elements
 *
 * @param[in] array Array.
 * @return Array containing the atanh of every element.
 ***************************************************************************/
GNdarray atanh(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = atanh(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes cosine of array elements
 *
 * @param[in] array Array.
 * @return Array containing the cosine of every element.
 ***************************************************************************/
GNdarray cos(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::cos(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes cosh of array elements
 *
 * @param[in] array Array.
 * @return Array containing the cosh of every element.
 ***************************************************************************/
GNdarray cosh(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::cosh(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes exponential of array elements
 *
 * @param[in] array Array.
 * @return Array containing the exponential of every element.
 ***************************************************************************/
GNdarray exp(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::exp(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes absolute of array elements
 *
 * @param[in] array Array.
 * @return Array containing the absolute of every element.
 ***************************************************************************/
GNdarray abs(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::abs(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes natural logarithm of array elements
 *
 * @param[in] array Array.
 * @return Array containing the natural logarithm of every element.
 ***************************************************************************/
GNdarray log(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::log(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes base10 logarithm of array elements
 *
 * @param[in] array Array.
 * @return Array containing the base10 logarithm of every element.
 ***************************************************************************/
GNdarray log10(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::log10(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes sine of array elements
 *
 * @param[in] array Array.
 * @return Array containing the sine of every element.
 ***************************************************************************/
GNdarray sin(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::sin(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes sinh of array elements
 *
 * @param[in] array Array.
 * @return Array containing the sinh of every element.
 ***************************************************************************/
GNdarray sinh(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::sinh(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes square root of array elements
 *
 * @param[in] array Array.
 * @return Array containing the square root of every element.
 ***************************************************************************/
GNdarray sqrt(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::sqrt(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes tangens of array elements
 *
 * @param[in] array Array.
 * @return Array containing the tangens of every element.
 ***************************************************************************/
GNdarray tan(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::tan(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes tanh of array elements
 *
 * @param[in] array Array.
 * @return Array containing the tanh of every element.
 ***************************************************************************/
GNdarray tanh(const GNdarray& array)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::tanh(array.m_data[i]);
    }

    // Return array
    return result;
}


/***********************************************************************//**
 * @brief Computes tanh of array elements
 *
 * @param[in] array Array.
 * @param[in] power Power.
 * @return Array containing the power of every element.
 ***************************************************************************/
GNdarray pow(const GNdarray& array, const double& power)
{
    // Initialise result array
    GNdarray result = array;

    // Evaluate each array element
    for (int i = 0; i < array.m_data.size(); ++i) {
        result.m_data[i] = std::pow(array.m_data[i], power);
    }

    // Return array
    return result;
}
