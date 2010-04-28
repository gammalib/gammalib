/***************************************************************************
 *                       GVector.cpp  -  vector class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2010         by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GVector.cpp
 * @brief GVector class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GVector.hpp"


/* __ Method name definitions ____________________________________________ */
#define G_CROSS                                      "cross(GVector,GVector)"


/*==========================================================================
 =                                                                         =
 =                      GVector constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void vector constructor (contains no elements)
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
 * Initialises a vector with num elements (all values are set to 0).
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
 * @param[in] a Value of first and single vector element.
 *
 * Initialises 1-element vector.
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
 * @param[in] a Value of first vector element.
 * @param[in] b Value of second vector element.
 *
 * Initialises 2-elements vector.
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
 * @param[in] a Value of first vector element.
 * @param[in] b Value of second vector element.
 * @param[in] c Value of third vector element.
 *
 * Initialises 3-elements vector.
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
 * @param[in] v Vector from which class should be instantiated.
 ***************************************************************************/
GVector::GVector(const GVector& v)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(v);

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
 =                            GVector operators                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] v GVector instance to be assigned
 ***************************************************************************/
GVector& GVector::operator= (const GVector& v)
{
    // Execute only if object is not identical
    if (this != &v) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(v);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GVector public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set vector elements to 0
 ***************************************************************************/
void GVector::clear(void)
{
    // Set all elements to 0
    for (int i = 0; i < m_num; ++i)
        m_data[i] = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns number of non-zero elements in vector
 ***************************************************************************/
int GVector::non_zeros(void) const
{
    // Initialise number of non-zeros
    int non_zeros = 0;

    // Gather all non-zero elements
    for (int i = 0; i < m_num; ++i) {
        if (m_data[i] != 0.0)
            non_zeros++;
    }

    // Return number of non-zero elements
    return non_zeros;
}


/*==========================================================================
 =                                                                         =
 =                          GVector private methods                        =
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
        for (int i = 0; i < m_num; ++i)
            m_data[i] = 0.0;

    } // endif: vector had non-zero length

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] v Vector from which members should be copied.
 ***************************************************************************/
void GVector::copy_members(const GVector& v)
{
    // Copy attributes
    m_num = v.m_num;

    // Copy elements
    if (m_num > 0) {
        alloc_members();
        for (int i = 0; i <  m_num; ++i)
            m_data[i] = v.m_data[i];
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
    if (m_data != NULL) delete m_data;

    // Signal free pointers
    m_data = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              GVector friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] v Vector to put in output stream.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GVector& v)
{
    // Prepend (
    os << "(";

    // Put all elements in stream
    for (int i = 0; i < v.m_num; ++i) {
        os << v(i);
        if (i != v.m_num-1)
            os << ", ";
    }

    // Append )
    os << ")";

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief GVector cross product
 *
 * @param[in] a First vector for cross product.
 * @param[in] b Second vector for cross product.
 *
 * Computes the cross product between two 3-element vectors (note that the
 * cross product is only defined for 3-element vectors).
 ***************************************************************************/
GVector cross(const GVector &a, const GVector &b)
{
    // Verify that vectors have same dimensions
    if (a.m_num != b.m_num)
        throw GException::vector_mismatch(G_CROSS, a.m_num, b.m_num);

    // Verify that vectors have 3 elements
    if (a.m_num != 3)
       throw GException::vector_bad_cross_dim(G_CROSS, a.m_num);

    // Compute cross product
    GVector result(3);
    result(0) = a.m_data[1]*b.m_data[2] - a.m_data[2]*b.m_data[1];
    result(1) = a.m_data[2]*b.m_data[0] - a.m_data[0]*b.m_data[2];
    result(2) = a.m_data[0]*b.m_data[1] - a.m_data[1]*b.m_data[0];

    // Return result
    return result;
}
