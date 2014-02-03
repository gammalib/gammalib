/***************************************************************************
 *             GCTAResponseTable.cpp - CTA response table class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAResponseTable.cpp
 * @brief CTA response table class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GException.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GCTAException.hpp"
#include "GCTAResponseTable.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR1                  "GCTAResponseTable::operator()(double&)"
#define G_OPERATOR2         "GCTAResponseTable::operator()(double&, double&)"
#define G_OPERATOR3         "GCTAResponseTable::operator()(double&, double&,"\
                                                                  " double&)"
#define G_INX_OPERATOR1        "GCTAResponseTable::operator()(int&, double&)"
#define G_INX_OPERATOR2        "GCTAResponseTable::operator()(int&, double&,"\
                                                                  " double&)"
#define G_INX_OPERATOR3        "GCTAResponseTable::operator()(int&, double&,"\
                                                         " double&, double&)"
#define G_AXIS                                "GCTAResponseTable::axis(int&)"
#define G_AXIS_LO_UNIT                "GCTAResponseTable::axis_lo_unit(int&)"
#define G_AXIS_HI_UNIT                "GCTAResponseTable::axis_hi_unit(int&)"
#define G_UNIT                                "GCTAResponseTable::unit(int&)"
#define G_AXIS_LO                    "GCTAResponseTable::axis_lo(int&, int&)"
#define G_AXIS_HI                    "GCTAResponseTable::axis_hi(int&, int&)"
#define G_AXIS_LINEAR                  "GCTAResponseTable::axis_linear(int&)"
#define G_AXIS_LOG10                    "GCTAResponseTable::axis_log10(int&)"
#define G_AXIS_RADIANS                "GCTAResponseTable::axis_radians(int&)"
#define G_SCALE                     "GCTAResponseTable::scale(int&, double&)"
#define G_READ                         "GCTAResponseTable::read(GFitsTable*)"
#define G_READ_COLNAMES       "GCTAResponseTable::read_colnames(GFitsTable*)"
#define G_READ_AXES               "GCTAResponseTable::read_axes(GFitsTable*)"
#define G_READ_PARS               "GCTAResponseTable::read_pars(GFitsTable*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Construct an empty CTA response table object.
 ***************************************************************************/
GCTAResponseTable::GCTAResponseTable(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] table Response table.
 *
 * Construct a CTA response table object by copying information from an
 * existing object. A deep copy is performed for construction, so that the
 * original object may be destroyed after construction the new object.
 ***************************************************************************/
GCTAResponseTable::GCTAResponseTable(const GCTAResponseTable& table)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS table constructor
 *
 * @param[in] hdu FITS table HDU.
 *
 * Construct FITS table object from the information that is found in the FITS
 * HDU. This method uses the GCTAResponseTable::read method to read the
 * response table information from the FITS HDU. Please refer to this method
 * for more details about what information is read from the HDU.
 *
 * If the FITS table HDU is not valid (i.e. NULL), the method initialises
 * the object to a clean state (equivalent of GCTAResponseTable::clean
 * method).
 ***************************************************************************/
GCTAResponseTable::GCTAResponseTable(const GFitsTable& hdu)
{
    // Initialise class members
    init_members();

    // Read table information from HDU
    read(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys the CTA response table.
 ***************************************************************************/
GCTAResponseTable::~GCTAResponseTable(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] table Response table.
 * @return Response table.
 *
 * Assigns a CTA response table to another object. The method is performing
 * a deep copy of the response table information, so that the original
 * object can be destroyed after assignment without any loss of information.
 ***************************************************************************/
GCTAResponseTable& GCTAResponseTable::operator=(const GCTAResponseTable& table)
{
    // Execute only if object is not identical
    if (this != &table) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(table);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Linear interpolation operator for 1D tables
 *
 * @param[in] arg Value.
 * @return Linearly interpolated response parameter vector.
 *
 * @exception GCTAException::bad_rsp_table_dim
 *            Response table has less than one dimension.
 *
 * Evaluates all response parameters at a given value for a one-dimensional
 * parameter vector. The evaluation is performed by a linear interpolation
 * of the vector. If the specified value lies outside the range covered by
 * the vector, the parameter is linearily extrapolated from using the first
 * (or last) two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
std::vector<double> GCTAResponseTable::operator()(const double& arg) const
{
    // Get parameter vector size
    int num = m_colname_par.size();

    // Optionally check that we have at least one dimension
    #if defined(G_RANGE_CHECK)
    if (num < 1) {
        throw GCTAException::bad_rsp_table_dim(G_OPERATOR1, num, 1);
    }
    #endif
    
    // Initialise result vector
    std::vector<double> result(num);
    
    // Set indices and weighting factors for interpolation
    update(arg);

    // Perform 1D interpolation
    for (int i = 0; i < num; ++i) {
        result[i] = m_wgt_left  * m_pars[i][m_inx_left] +
                    m_wgt_right * m_pars[i][m_inx_right];
    }
    
    // Return result vector
    return result;
}


/***********************************************************************//**
 * @brief Bilinear interpolation operator for 2D tables
 *
 * @param[in] arg1 Value for first axis.
 * @param[in] arg2 Value for second axis.
 * @return Bilinearly interpolated response parameter vector.
 *
 * @exception GCTAException::bad_rsp_table_dim
 *            Response table has less than one dimension.
 *
 * Evaluates all response parameters at a given pair of values for a
 * two-dimensional parameter vector. The evaluation is performed by a linear
 * interpolation of the vector. If the specified value lies outside the range
 * covered by the vector, the parameter is linearily extrapolated from using
 * the first (or last) two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
std::vector<double> GCTAResponseTable::operator()(const double& arg1,
                                                  const double& arg2) const
{
    // Get parameter vector size
    int num = m_colname_par.size();

    // Optionally check that we have at least two dimensions
    #if defined(G_RANGE_CHECK)
    if (num < 2) {
        throw GCTAException::bad_rsp_table_dim(G_OPERATOR2, num, 2);
    }
    #endif
    
    // Initialise result vector
    std::vector<double> result(num);

    // Set indices and weighting factors for interpolation
    update(arg1, arg2);

    // Perform 2D interpolation
    for (int i = 0; i < num; ++i) {
        result[i] = m_wgt1 * m_pars[i][m_inx1] +
                    m_wgt2 * m_pars[i][m_inx2] +
                    m_wgt3 * m_pars[i][m_inx3] +
                    m_wgt4 * m_pars[i][m_inx4];
    }
    
    // Return result vector
    return result;
}
/***********************************************************************//**
 * @brief Trilinear interpolation operator for 3D tables
 *
 * @param[in] arg1 Value for first axis.
 * @param[in] arg2 Value for second axis.
 * @param[in] arg3 Value for third axis.
 * @return Trilinearly interpolated response parameter vector.
 *
 * @exception GCTAException::bad_rsp_table_dim
 *            Response table has less than one dimension.
 *
 * Evaluates all response parameters at a given triplet of values for a
 * three-dimensional parameter vector. The evaluation is performed by a
 * tri-linear interpolation of the vector. If the specified value lies
 * outside the range covered by the vector, the parameter is linearily
 * extrapolated from using the first (or last) two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
std::vector<double> GCTAResponseTable::operator()(const double& arg1,
                                                  const double& arg2,
                                                  const double& arg3) const
{
    // Get parameter vector size
    int num = m_colname_par.size();

    // Optionally check that we have at least three dimensions
    #if defined(G_RANGE_CHECK)
    if (num < 3) {
        throw GCTAException::bad_rsp_table_dim(G_OPERATOR3, num, 3);
    }
    #endif
    
    // Initialise result vector
    std::vector<double> result(num);

    // Set indices and weighting factors for interpolation
    update(arg1, arg2, arg3);

    // Perform 3D interpolation
    for (int i = 0; i < num; ++i) {
        result[i] = m_wgt1 * m_pars[i][m_inx1] +
                    m_wgt2 * m_pars[i][m_inx2] +
                    m_wgt3 * m_pars[i][m_inx3] +
                    m_wgt4 * m_pars[i][m_inx4] +
                    m_wgt5 * m_pars[i][m_inx5] +
                    m_wgt6 * m_pars[i][m_inx6] +
                    m_wgt7 * m_pars[i][m_inx7] +
                    m_wgt8 * m_pars[i][m_inx8] ;
    }

    // Return result vector
    return result;
}


/***********************************************************************//**
 * @brief Linear interpolation operator for 1D tables
 *
 * @param[in] index Table index [0,...,size()-1].
 * @param[in] arg Value.
 * @return Linearly interpolated response parameter.
 *
 * @exception GCTAException::bad_rsp_table_dim
 *            Response table has less than one dimension.
 *
 * Evaluates one response parameter at a given value for a one-dimensional
 * parameter vector. The evaluation is performed by a linear interpolation
 * of the vector. If the specified value lies outside the range covered by
 * the vector, the parameter is linearily extrapolated from using the first
 * (or last) two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
double GCTAResponseTable::operator()(const int& index, const double& arg) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_INX_OPERATOR1, index, size()-1);
    }
    #endif
    
    // Set indices and weighting factors for interpolation
    update(arg);

    // Perform 1D interpolation
    double result = m_wgt_left  * m_pars[index][m_inx_left] +
                    m_wgt_right * m_pars[index][m_inx_right];
    
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Bilinear interpolation operator for 2D tables
 *
 * @param[in] index Table index [0,...,size()-1].
 * @param[in] arg1 Value for first axis.
 * @param[in] arg2 Value for second axis.
 * @return Bilinearly interpolated response parameter vector.
 *
 * @exception GCTAException::bad_rsp_table_dim
 *            Response table has less than one dimension.
 *
 * Evaluates all response parameters at a given pair of values for a
 * two-dimensional parameter vector. The evaluation is performed by a linear
 * interpolation of the vector. If the specified value lies outside the range
 * covered by the vector, the parameter is linearily extrapolated from using
 * the first (or last) two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
double GCTAResponseTable::operator()(const int& index, const double& arg1,
                                     const double& arg2) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_INX_OPERATOR2, index, size()-1);
    }
    #endif

    // Set indices and weighting factors for interpolation
    update(arg1, arg2);

    // Perform 2D interpolation
    double result = m_wgt1 * m_pars[index][m_inx1] +
                    m_wgt2 * m_pars[index][m_inx2] +
                    m_wgt3 * m_pars[index][m_inx3] +
                    m_wgt4 * m_pars[index][m_inx4];
    
    // Return result
    return result;
}
/***********************************************************************//**
 * @brief Trilinear interpolation operator for 3D tables
 *
 * @param[in] index Table index [0,...,size()-1].
 * @param[in] arg1 Value for first axis.
 * @param[in] arg2 Value for second axis.
 * @param[in] arg3 Value for second axis.
 * @return Trilinearly interpolated response parameter vector.
 *
 * @exception GCTAException::bad_rsp_table_dim
 *            Response table has less than one dimension.
 *
 * Evaluates all response parameters at a given triplet of values for a
 * three-dimensional parameter vector. The evaluation is performed by a
 * trilinear interpolation of the vector. If the specified value lies outside
 * the range covered by the vector, the parameter is linearily extrapolated
 * from using the first (or last) two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
double GCTAResponseTable::operator()(const int&    index,
                                     const double& arg1, 
                                     const double& arg2,
                                     const double& arg3) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_INX_OPERATOR3, index, size()-1);
    }
    #endif

    // Set indices and weighting factors for interpolation
    update(arg1, arg2, arg3);

    // Perform 3D interpolation
    double result = m_wgt1 * m_pars[index][m_inx1] +
                    m_wgt2 * m_pars[index][m_inx2] +
                    m_wgt3 * m_pars[index][m_inx3] +
                    m_wgt4 * m_pars[index][m_inx4] +
                    m_wgt5 * m_pars[index][m_inx5] +
                    m_wgt6 * m_pars[index][m_inx6] +
                    m_wgt7 * m_pars[index][m_inx7] +
                    m_wgt8 * m_pars[index][m_inx8];

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * Resets the object to a proper initial state. Any information that resided
 * before in the object will be lost.
 ***************************************************************************/
void GCTAResponseTable::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of response table instance.
 *
 * Clones an instance of the CTA response table. The cloned object will be
 * a deep copy of the original object, hence after cloning the original
 * object can be destroyed without any loss of information.
 ***************************************************************************/
GCTAResponseTable* GCTAResponseTable::clone(void) const
{
    return new GCTAResponseTable(*this);
}


/***********************************************************************//**
 * @brief Return axis length
 *
 * @param[in] index Axis index [0,...,axes()-1].
 * @return Number of bins along the specified axis.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the number of bins along the specified axis.
 ***************************************************************************/
int GCTAResponseTable::axis(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= axes()) {
        throw GException::out_of_range(G_AXIS, index, axes()-1);
    }
    #endif

    // Return axis length
    return (m_axis_lo[index].size());
}


/***********************************************************************//**
 * @brief Return axis lower boundary unit
 *
 * @param[in] index Axis index [0,...,axes()-1].
 * @return Unit of lower boundary.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the unit of the lower boundary for the specified axis.
 ***************************************************************************/
std::string GCTAResponseTable::axis_lo_unit(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= axes()) {
        throw GException::out_of_range(G_AXIS_LO_UNIT, index, axes()-1);
    }
    #endif

    // Return units
    return (m_units_lo[index]);
}


/***********************************************************************//**
 * @brief Return axis upper boundary unit
 *
 * @param[in] index Axis index [0,...,axes()-1].
 * @return Unit of upper boundary.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the unit of the upper boundary for the specified axis.
 ***************************************************************************/
std::string GCTAResponseTable::axis_hi_unit(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= axes()) {
        throw GException::out_of_range(G_AXIS_HI_UNIT, index, axes()-1);
    }
    #endif

    // Return units
    return (m_units_hi[index]);
}


/***********************************************************************//**
 * @brief Return parameter unit
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Unit of parameter.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the unit of the parameter for the specified index.
 ***************************************************************************/
std::string GCTAResponseTable::unit(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_UNIT, index, size()-1);
    }
    #endif

    // Return units
    return (m_units_par[index]);
}


/***********************************************************************//**
 * @brief Return lower bin boundary for bin in axis
 *
 * @param[in] index Axis index [0,...,axes()-1].
 * @param[in] bin Bin index (starting from 0).
 * @return Lower bin boundary.
 *
 * @exception GException::out_of_range
 *            Axis or bin index out of range.
 *
 * Returns the lower boundary for a given bin of a given axis.
 ***************************************************************************/
double GCTAResponseTable::axis_lo(const int& index, const int& bin) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= axes() ||
        bin   < 0 || bin   >= m_axis_lo[index].size()) {
        throw GException::out_of_range(G_AXIS_LO, index, bin,
                                       axes()-1, m_axis_lo[index].size()-1);
    }
    #endif

    // Return bin boundary
    return (m_axis_lo[index][bin]);
}


/***********************************************************************//**
 * @brief Return upper bin boundary for bin in axis
 *
 * @param[in] index Axis index [0,...,axes()-1].
 * @param[in] bin Bin index (starting from 0).
 * @return Uppser bin boundary.
 *
 * @exception GException::out_of_range
 *            Axis or bin index out of range.
 *
 * Returns the upper boundary for a given bin of a given axis.
 ***************************************************************************/
double GCTAResponseTable::axis_hi(const int& index, const int& bin) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= axes() ||
        bin   < 0 || bin   >= m_axis_hi[index].size()) {
        throw GException::out_of_range(G_AXIS_HI, index, bin,
                                       axes()-1, m_axis_hi[index].size()-1);
    }
    #endif

    // Return bin boundary
    return (m_axis_hi[index][bin]);
}


/***********************************************************************//**
 * @brief Set nodes for a linear axis
 *
 * @param[in] index Axis index [0,...,axes()-1].
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Set axis nodes so that each node is the linear mean of the lower and upper
 * bin boundary, i.e.
 * \f[ n_i = 0.5 \times ({\rm LO}_i + {\rm HI}_i) \f]
 * where
 * \f$n_i\f$ is node \f$i\f$,
 * \f${\rm LO}_i\f$ is the lower bin boundary for bin \f$i\f$, and
 * \f${\rm HI}_i\f$ is the upper bin boundary for bin \f$i\f$.
 ***************************************************************************/
void GCTAResponseTable::axis_linear(const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= axes()) {
        throw GException::out_of_range(G_AXIS_LINEAR, index, axes()-1);
    }
    #endif

    // Get number of bins
    int bins = axis(index);

    // Allocate node values
    std::vector<double> axis_nodes(bins);

    // Compute nodes
    for (int i = 0; i < bins; ++i) {
        axis_nodes[i] = 0.5*(m_axis_lo[index][i] + m_axis_hi[index][i]);
    }

    // Set node array
    m_axis_nodes[index] = GNodeArray(axis_nodes);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for a logarithmic (base 10) axis
 *
 * @param[in] index Axis index [0,...,axes()-1].
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Set axis nodes so that each node is the logarithmic mean of the lower and
 * upper bin boundary, i.e.
 * \f[ n_i = \log \sqrt{{\rm LO}_i \times {\rm HI}_i} \f]
 * where
 * \f$n_i\f$ is node \f$i\f$,
 * \f${\rm LO}_i\f$ is the lower bin boundary for bin \f$i\f$, and
 * \f${\rm HI}_i\f$ is the upper bin boundary for bin \f$i\f$.
 *
 * @todo Check that none of the axis boundaries is non-positive.
 ***************************************************************************/
void GCTAResponseTable::axis_log10(const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= axes()) {
        throw GException::out_of_range(G_AXIS_LINEAR, index, axes()-1);
    }
    #endif

    // Get number of bins
    int bins = axis(index);

    // Allocate node values
    std::vector<double> axis_nodes(bins);

    // Compute nodes
    for (int i = 0; i < bins; ++i) {
        axis_nodes[i] = std::log10(std::sqrt(m_axis_lo[index][i] * 
                                             m_axis_hi[index][i]));
    }

    // Set node array
    m_axis_nodes[index] = GNodeArray(axis_nodes);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for a radians axis
 *
 * @param[in] index Axis index [0,...,axes()-1].
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Set axis nodes so that each node is the lo mean of the lower and upper
 * bin boundary in radians, i.e.
 * \f[ n_i = \frac{\pi}{360} \times ({\rm LO}_i + {\rm HI}_i) \f]
 * where
 * \f$n_i\f$ is node \f$i\f$,
 * \f${\rm LO}_i\f$ is the lower bin boundary for bin \f$i\f$, and
 * \f${\rm HI}_i\f$ is the upper bin boundary for bin \f$i\f$.
 ***************************************************************************/
void GCTAResponseTable::axis_radians(const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= axes()) {
        throw GException::out_of_range(G_AXIS_RADIANS, index, axes()-1);
    }
    #endif

    // Get number of bins
    int bins = axis(index);

    // Allocate node values
    std::vector<double> axis_nodes(bins);

    // Compute nodes
    for (int i = 0; i < bins; ++i) {
        axis_nodes[i] = 0.5*(m_axis_lo[index][i] + m_axis_hi[index][i]) *
                        gammalib::deg2rad;
    }

    // Set node array
    m_axis_nodes[index] = GNodeArray(axis_nodes);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for a radians axis
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @param[in] scale Scaling factor.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Multiplies all values in a parameter table by the specified scaling
 * factor.
 ***************************************************************************/
void GCTAResponseTable::scale(const int& index, const double& scale)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SCALE, index, size()-1);
    }
    #endif

    // Scale parameter values
    for (int i = 0; i < m_nelements; ++i) {
        m_pars[index][i] *= scale;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read response table from FITS table HDU
 *
 * @param[in] table Response table.
 *
 * Reads CTA response table information from a FITS table. The FITS table
 * is expected to have a single row, and axes and parameter information are
 * found in vector columns. Axes information is expected to be placed
 * before parameter information. 
 *
 * Each axis is defined by two vector columns of equal width, describing the
 * lower and upper limits for each bin. The column names for this boundary
 * information terminates by "_LO" and "HI", respectively (upper case). It
 * is expected that the "_LO" column preceeds the "_HI" column.
 *
 * Following the axes columns are parameter columns of equal width. The width
 * of each parameter column is given by the product of the lengths of all
 * axes. It is furthermore expected that the first axis is the most rapidely
 * varying index of the vector.
 *
 * In case that the HDU table pointer is not valid (i.e. NULL), this method
 * clears the objects and does nothing else.
 ***************************************************************************/
void GCTAResponseTable::read(const GFitsTable& table)
{
    // Clear instance
    clear();

    // Read column names
    read_colnames(table);

    // Read axes
    read_axes(table);

    // Read parameter cubes
    read_pars(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write response table into FITS table
 *
 * @param[in] table Response table.
 *
 * @todo Implement method
 ***************************************************************************/
void GCTAResponseTable::write(GFitsTable& table) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print response table information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing response table information.
 *
 * Puts CTA response table information into a std::string object for
 * printing.
 ***************************************************************************/
std::string GCTAResponseTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAResponseTable ===");

        // Append information
        result.append("\n"+gammalib::parformat("Dimension") +
                      gammalib::str(axes()));
    
        // Append axes information
        for (int i = 0; i < axes(); ++i) {
            result.append("\n"+gammalib::parformat("Axis " +
                          gammalib::str(i)));
            result.append(gammalib::str(axis(i)));
            result.append(" (");
            result.append(m_colname_lo[i]);
            if (!m_units_lo[i].empty()) {
                result.append(" ["+m_units_lo[i]+"]");
            }
            result.append(", ");
            result.append(m_colname_hi[i]);
            if (!m_units_hi[i].empty()) {
                result.append(" ["+m_units_hi[i]+"]");
            }
            result.append(")");
        }

        // Append parameter information
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+gammalib::parformat("Parameter " +
                          gammalib::str(i)));
            result.append(m_colname_par[i]);
            if (!m_units_par[i].empty()) {
                result.append(" ["+m_units_par[i]+"]");
            }
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 *
 * Initialises all members of the class.
 ***************************************************************************/
void GCTAResponseTable::init_members(void)
{
    // Initialise members
    m_naxes     = 0;
    m_npars     = 0;
    m_nelements = 0;
    m_colname_lo.clear();
    m_colname_hi.clear();
    m_colname_par.clear();
    m_axis_lo.clear();
    m_axis_hi.clear();
    m_units_lo.clear();
    m_units_hi.clear();
    m_units_par.clear();
    m_axis_nodes.clear();
    m_pars.clear();

    // Initialise cache
    m_inx_left  = 0;
    m_inx_right = 0;
    m_wgt_left  = 0.0;
    m_wgt_right = 0.0;
    m_inx1      = 0;
    m_inx2      = 0;
    m_inx3      = 0;
    m_inx4      = 0;
    m_inx5      = 0;
    m_inx6      = 0;
    m_inx7      = 0;
    m_inx8      = 0;
    m_wgt1      = 0.0;
    m_wgt2      = 0.0;
    m_wgt3      = 0.0;
    m_wgt4      = 0.0;
    m_wgt5      = 0.0;
    m_wgt6      = 0.0;
    m_wgt7      = 0.0;
    m_wgt8      = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] table Response table.
 *
 * Copies all members of the class.
 ***************************************************************************/
void GCTAResponseTable::copy_members(const GCTAResponseTable& table)
{
    // Copy number of bins
    m_naxes       = table.m_naxes;
    m_npars       = table.m_npars;
    m_nelements   = table.m_nelements;
    m_colname_lo  = table.m_colname_lo;
    m_colname_hi  = table.m_colname_hi;
    m_colname_par = table.m_colname_par;
    m_axis_lo     = table.m_axis_lo;
    m_axis_hi     = table.m_axis_hi;
    m_units_lo    = table.m_units_lo;
    m_units_hi    = table.m_units_hi;
    m_units_par   = table.m_units_par;
    m_axis_nodes  = table.m_axis_nodes;
    m_pars        = table.m_pars;

    // Copy cache
    m_inx_left  = table.m_inx_left;
    m_inx_right = table.m_inx_right;
    m_wgt_left  = table.m_wgt_left;
    m_wgt_right = table.m_wgt_right;
    m_inx1      = table.m_inx1;
    m_inx2      = table.m_inx2;
    m_inx3      = table.m_inx3;
    m_inx4      = table.m_inx4;
    m_inx5      = table.m_inx5;
    m_inx6      = table.m_inx6;
    m_inx7      = table.m_inx7;
    m_inx8      = table.m_inx8;
    m_wgt1      = table.m_wgt1;
    m_wgt2      = table.m_wgt2;
    m_wgt3      = table.m_wgt3;
    m_wgt4      = table.m_wgt4;
    m_wgt5      = table.m_wgt5;
    m_wgt6      = table.m_wgt6;
    m_wgt7      = table.m_wgt7;
    m_wgt8      = table.m_wgt8;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * Deallocates any memory that was allocated by the object.
 ***************************************************************************/
void GCTAResponseTable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read column names from FITS HDU
 *
 * @param[in] hdu Response table HDU.
 *
 * @exception GCTAException::bad_rsp_table_format
 *            Bad response table format encountered in FITS HDU.
 *
 * Read the response table column names from the HDU. Column names
 * terminating with "_LO" and "_HI" define the parameter space axes, while
 * all other column names define response parameter data cubes. It is
 * assumed that parameter space axes are given in subsequent order, with
 * the first column corresponding to the lower bin boundaries of the first
 * axis. Lower bin boundaries are indicated by the "_LO" termination. Upper
 * bin boundaries are assumed to follow immediately the lower bin boundaires
 * and are designated by the "_HI" termination.
 *
 * This method sets the following members:
 *    m_colname_lo - Column names of lower boundaries
 *    m_colname_hi - Column names of upper boundaries
 *    m_colname_par - Column names of parameters
 *    m_naxes - Number of axes
 *    m_npars - Number of parameters
 *
 * In case that the HDU pointer is not valid (NULL), this method clears the
 * column names and does nothing else.
 *
 * @todo Implement exceptions for invalid HDU format
 ***************************************************************************/
void GCTAResponseTable::read_colnames(const GFitsTable& hdu)
{
    // Clear column name arrays
    m_naxes = 0;
    m_npars = 0;
    m_colname_lo.clear();
    m_colname_hi.clear();
    m_colname_par.clear();
    
    // Initialise search mode. There are three search modes:
    // 0 - we're looking for the next axis by searching for a column
    //     terminating with "_LO"
    // 1 - we're looking for the upper boundary of an axis, terminating
    //     with "_HI"
    // 2 - we're looking for a parameter column
    int         mode = 0;
    std::string lo_column;
    std::string next_column;

    // Extract column names for all axes
    for (int i = 0; i < hdu.ncols(); ++i) {

        // Get column name and unit
        std::string colname = hdu[i]->name();

        // If we search for a "_LO" column, check if we have one. If one
        // is found, change the search mode to 1 and set the expected name
        // for the "_HI" column. If none is found, change the search
        // mode to 2 since from now on we should only have parameter
        // columns.
        if (mode == 0) {
            size_t pos = colname.rfind("_LO");
            if (pos != std::string::npos) {
                mode        = 1;
                lo_column   = colname;
                next_column = colname.substr(0, pos) + "_HI";
            }
            else {
                if (colname.rfind("_HI") != std::string::npos) {
                    std::string message = "Column '" +
                                          colname +
                                          "' encountered without"
                                          " preceeding '_LO' column.";
                    throw GCTAException::bad_rsp_table_format(G_READ_COLNAMES,
                                                              message);
                }
                else {
                    mode = 2;
                    m_colname_par.push_back(colname);
                }
            }
        }

        // If we search for a "_HI" column, check if we have the
        // expected column name. If this is the case, switch back to
        // search mode 0 to get the next "_LO" column. Otherwise
        // throw an exception.
        else if (mode == 1) {
            if (colname == next_column) {
                mode = 0;
                m_colname_lo.push_back(lo_column);
                m_colname_hi.push_back(next_column);
            }
            else {
                std::string message = "Expected column '" +
                                      next_column +
                                      "' not found. '_HI' columns have"
                                      " to be placed immediately after"
                                      " corresponding '_LO' columns.";
                throw GCTAException::bad_rsp_table_format(G_READ_COLNAMES,
                                                          message);
            }
        }

        // If we search for a parameter column, make sure that we have
        // neither a "_LO" nor a "_HI" column.
        else {
            if (colname.rfind("_LO") != std::string::npos) {
                std::string message = "Column '" +
                                      colname +
                                      "' found. All '_LO' columns have to"
                                      " be placed before the parameter"
                                      " columns.";
                throw GCTAException::bad_rsp_table_format(G_READ_COLNAMES,
                                                          message);
            }
            else if (colname.rfind("_HI") != std::string::npos) {
                std::string message = "Column '" +
                                      colname +
                                      "' found. All '_HI' columns have to"
                                      " be placed before the parameter"
                                      " columns.";
                throw GCTAException::bad_rsp_table_format(G_READ_COLNAMES,
                                                          message);
            }
            else {
                m_colname_par.push_back(colname);
            }
        }

    } // endfor: looped over all columns

    // Store number of axes
    m_naxes = m_colname_lo.size();

    // Store number of parameters
    m_npars = m_colname_par.size();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read axes definitions from FITS HDU
 *
 * @param[in] hdu Response table HDU.
 *
 * @exception GCTAException::bad_rsp_table_format
 *            Incompatible axes columns found.
 *
 * Reads the lower and upper boundaries for all response table axes from the
 * FITS HDU. The method verifies if the lower and upper boundary axes have
 * the same vector dimensions. 
 *
 * The method also allocates the nodes for each axis, assuming that axis will
 * be used linearily, where each node is given by the linear mean of the
 * lower and upper boundary. In case that the axis should be used
 * logarithmically, the method axis_log10() can be used to change to a
 * logarithmic scale. The axis_linear() can be used to switch back to a
 * linear scale.
 *
 * This method sets the following members:
 *    m_axis_lo - Axes lower boundaries
 *    m_axis_hi - Axes upper boundaries
 *    m_axis_nodes - Axes mean values
 *
 * In case that the HDU pointer is not valid (NULL), this method clears the
 * axes boundaries and does nothing else.
 ***************************************************************************/
void GCTAResponseTable::read_axes(const GFitsTable& hdu)
{
    // Clear axes arrays
    m_axis_lo.clear();
    m_axis_hi.clear();
    m_axis_nodes.clear();
    
    // Loop over all dimensions
    for (int i = 0; i < axes(); ++i) {

        // Get pointers to table columns
        const GFitsTableCol* col_lo = hdu[m_colname_lo[i]];
        const GFitsTableCol* col_hi = hdu[m_colname_hi[i]];

        // Extract number of bins. Make sure that both columns have the
        // same number of bins
        int num = col_lo->number();
        if (num != col_hi->number()) {
            std::string message = "Incompatible number of bins found"
                                  " for columns '"+m_colname_lo[i]+"'"+
                                  " ("+gammalib::str(num)+") and '"+
                                  m_colname_hi[i]+"' ("+
                                  gammalib::str(col_hi->number())+".";
            throw GCTAException::bad_rsp_table_format(G_READ_AXES,
                                                      message);
        }

        // Initialise axis and node arrays
        std::vector<double> axis_lo(num);
        std::vector<double> axis_hi(num);
        std::vector<double> axis_nodes(num);

        // Copy axis information into arrays
        for (int k = 0; k < num; ++k) {
            axis_lo[k]    = col_lo->real(0,k);
            axis_hi[k]    = col_hi->real(0,k);
            axis_nodes[k] = 0.5*(axis_lo[k] + axis_hi[k]);
        }

        // Push axis array on storage
        m_axis_lo.push_back(axis_lo);
        m_axis_hi.push_back(axis_hi);

        // Push units on storage
        m_units_lo.push_back(col_lo->unit());
        m_units_hi.push_back(col_hi->unit());

        // Create node array
        m_axis_nodes.push_back(GNodeArray(axis_nodes));

    } // endfor: looped over all dimensions

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read parameter cubes
 *
 * @param[in] hdu Response table HDU.
 *
 * @exception GCTAException::bad_rsp_table_format
 *            Parameter vector of bad size encountered.
 *
 * Reads the parameter cubes from the response table. The method also sets
 * the data members m_npars (number of parameters) and m_nelements (number
 * of elements per parameter).
 *
 * This method sets the following members:
 *    m_pars - Parameter values
 *    m_nelements - Number of elements per parameter
 *
 * In case that the HDU pointer is not valid (NULL), this method clears the
 * axes boundaries and does nothing else.
 ***************************************************************************/
void GCTAResponseTable::read_pars(const GFitsTable& hdu)
{
    // Clear parameter cubes
    m_pars.clear();
    
    // Compute expected cube size
    m_nelements = axis(0);
    for (int i = 1; i < axes(); ++i) {
        m_nelements *= axis(i);
    }
    
    // Loop over all parameter cubes
    for (int i = 0; i < size(); ++i) {

        // Get pointer to table column
        const GFitsTableCol* col = hdu[m_colname_par[i]];

        // Extract number of bins. Verify that the number of bins
        // corresponds to the expectation.
        int num = col->number();
        if (num != m_nelements) {
            std::string message = "Parameter vector '"+m_colname_par[i]+
                                  "' has wrong size "+gammalib::str(num)+
                                  " (expected"+
                                  gammalib::str(m_nelements)+").";
            throw GCTAException::bad_rsp_table_format(G_READ_PARS,
                                                      message);
        }

        // Initialise parameter cube
        std::vector<double> pars(num);

        // Copy parameter values
        for (int k = 0; k < num; ++k) {
            pars[k] = col->real(0,k);
        }

        // Push cube into storage
        m_pars.push_back(pars);

        // Push units on storage
        m_units_par.push_back(col->unit());

    } // endfor: looped over all parameter cubes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update 1D cache
 *
 * @param[in] arg Argument.
 *
 * Updates the 1D interpolation cache. The interpolation cache is composed
 * of two indices and weights that define 2 data values of the 2D table
 * that are used for linear interpolation.
 *
 * @todo Write down formula
 *
 * @todo Makes GNodeArray::set_value method const and use mutable members
 ***************************************************************************/
void GCTAResponseTable::update(const double& arg) const
{
    // Get pointer to node array (circumvent const correctness)
    GNodeArray* nodes = const_cast<GNodeArray*>(&(m_axis_nodes[0]));

    // Set value for node array
    nodes->set_value(arg);

    // Set indices and weighting factors for interpolation
    m_inx_left  = nodes->inx_left();
    m_inx_right = nodes->inx_right();
    m_wgt_left  = nodes->wgt_left();
    m_wgt_right = nodes->wgt_right();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update 2D cache
 *
 * @param[in] arg1 Argument for first axis.
 * @param[in] arg2 Argument for second axis.
 *
 * Updates the 2D interpolation cache. The interpolation cache is composed
 * of four indices and weights that define 4 data values of the 2D table
 * that are used for bilinear interpolation.
 *
 * @todo Write down formula
 *
 * @todo Makes GNodeArray::set_value method const and use mutable members
 ***************************************************************************/
void GCTAResponseTable::update(const double& arg1, const double& arg2) const
{
    // Get pointers to node arrays (circumvent const correctness)
    GNodeArray* nodes1 = const_cast<GNodeArray*>(&(m_axis_nodes[0]));
    GNodeArray* nodes2 = const_cast<GNodeArray*>(&(m_axis_nodes[1]));

    // Set values for node arrays
    nodes1->set_value(arg1);
    nodes2->set_value(arg2);

    // Compute offsets
    int size1        = axis(0);
    int offset_left  = nodes2->inx_left()  * size1;
    int offset_right = nodes2->inx_right() * size1;

    // Set indices for bi-linear interpolation
    m_inx1 = nodes1->inx_left()  + offset_left;
    m_inx2 = nodes1->inx_left()  + offset_right;
    m_inx3 = nodes1->inx_right() + offset_left;
    m_inx4 = nodes1->inx_right() + offset_right;

    // Set weighting factors for bi-linear interpolation
    m_wgt1 = nodes1->wgt_left()  * nodes2->wgt_left();
    m_wgt2 = nodes1->wgt_left()  * nodes2->wgt_right();
    m_wgt3 = nodes1->wgt_right() * nodes2->wgt_left();
    m_wgt4 = nodes1->wgt_right() * nodes2->wgt_right();
    
    // Return
    return;
}
/***********************************************************************//**
 * @brief Update 3D cache
 *
 * @param[in] arg1 Argument for first axis.
 * @param[in] arg2 Argument for second axis.
 * @param[in] arg3 Argument for second axis.
 *
 * Updates the 3D interpolation cache. The interpolation cache is composed
 * of four indices and weights that define 4 data values of the 2D table
 * that are used for bilinear interpolation.
 *
 * @todo Write down formula
 *
 * @todo Makes GNodeArray::set_value method const and use mutable members
 ***************************************************************************/
void GCTAResponseTable::update(const double& arg1, const double& arg2,
                               const double& arg3) const
{
    // Get pointers to node arrays (circumvent const correctness)
    GNodeArray* nodes1 = const_cast<GNodeArray*>(&(m_axis_nodes[0]));
    GNodeArray* nodes2 = const_cast<GNodeArray*>(&(m_axis_nodes[1]));
    GNodeArray* nodes3 = const_cast<GNodeArray*>(&(m_axis_nodes[2]));

    // Set values for node arrays
    nodes1->set_value(arg1);
    nodes2->set_value(arg2);
    nodes3->set_value(arg3);

    // Compute offsets
    int size1          = axis(0);
    int size2          = axis(1);
    int offset_left_2  = nodes2->inx_left()  * size1;
    int offset_right_2 = nodes2->inx_right() * size1;
    int offset_left_3  = nodes3->inx_left()  * size1 * size2;
    int offset_right_3 = nodes3->inx_right() * size1 * size2;

    // Set indices for tri-linear interpolation
    m_inx1 = nodes1->inx_left()  + offset_left_2  + offset_left_3 ;
    m_inx2 = nodes1->inx_left()  + offset_left_2  + offset_right_3;
    m_inx3 = nodes1->inx_left()  + offset_right_2 + offset_left_3 ;
    m_inx4 = nodes1->inx_left()  + offset_right_2 + offset_right_3;
    m_inx5 = nodes1->inx_right() + offset_left_2  + offset_left_3 ;
    m_inx6 = nodes1->inx_right() + offset_left_2  + offset_right_3;
    m_inx7 = nodes1->inx_right() + offset_right_2 + offset_left_3 ;
    m_inx8 = nodes1->inx_right() + offset_right_2 + offset_right_3;

    // Set weighting factors for tri-linear interpolation
    m_wgt1 = nodes1->wgt_left()  * nodes2->wgt_left()  *  nodes3->wgt_left();
    m_wgt2 = nodes1->wgt_left()  * nodes2->wgt_left()  *  nodes3->wgt_right();
    m_wgt3 = nodes1->wgt_left()  * nodes2->wgt_right() *  nodes3->wgt_left();
    m_wgt4 = nodes1->wgt_left()  * nodes2->wgt_right() *  nodes3->wgt_right();
    m_wgt5 = nodes1->wgt_right() * nodes2->wgt_left()  *  nodes3->wgt_left();
    m_wgt6 = nodes1->wgt_right() * nodes2->wgt_left()  *  nodes3->wgt_right();
    m_wgt7 = nodes1->wgt_right() * nodes2->wgt_right() *  nodes3->wgt_left();
    m_wgt8 = nodes1->wgt_right() * nodes2->wgt_right() *  nodes3->wgt_right();

    // Return
    return;
}
