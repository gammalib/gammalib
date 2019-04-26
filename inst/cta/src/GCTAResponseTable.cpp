/***************************************************************************
 *             GCTAResponseTable.cpp - CTA response table class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2018 by Juergen Knoedlseder                         *
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
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GCTAResponseTable.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR1                  "GCTAResponseTable::operator()(double&)"
#define G_OPERATOR2         "GCTAResponseTable::operator()(double&, double&)"
#define G_OPERATOR3         "GCTAResponseTable::operator()(double&, double&,"\
                                                                  " double&)"
#define G_ELEMENT_OPERATOR1             "GCTAResponseTable::operator()(int&)"
#define G_ELEMENT_OPERATOR2       "GCTAResponseTable::operator()(int&, int&)"
#define G_INX_OPERATOR1        "GCTAResponseTable::operator()(int&, double&)"
#define G_INX_OPERATOR2        "GCTAResponseTable::operator()(int&, double&,"\
                                                                  " double&)"
#define G_INX_OPERATOR3        "GCTAResponseTable::operator()(int&, double&,"\
                                                         " double&, double&)"
#define G_TABLE                      "GCTAResponseTable::table(std::string&)"
#define G_SCALE                     "GCTAResponseTable::scale(int&, double&)"
#define G_AXIS                        "GCTAResponseTable::axis(std::string&)"
#define G_AXIS_BINS                      "GCTAResponseTable::axis_bins(int&)"
#define G_AXIS_LO_NAME                "GCTAResponseTable::axis_lo_name(int&)"
#define G_AXIS_HI_NAME                "GCTAResponseTable::axis_hi_name(int&)"
#define G_AXIS_LO_UNIT                "GCTAResponseTable::axis_lo_unit(int&)"
#define G_AXIS_HI_UNIT                "GCTAResponseTable::axis_hi_unit(int&)"
#define G_UNIT                                "GCTAResponseTable::unit(int&)"
#define G_AXIS_LO                    "GCTAResponseTable::axis_lo(int&, int&)"
#define G_AXIS_HI                    "GCTAResponseTable::axis_hi(int&, int&)"
#define G_AXIS_LINEAR                  "GCTAResponseTable::axis_linear(int&)"
#define G_AXIS_LOG10                    "GCTAResponseTable::axis_log10(int&)"
#define G_AXIS_RADIANS                "GCTAResponseTable::axis_radians(int&)"
#define G_AXIS_NODES                    "GCTAResponseTable::axis_nodes(int&)"
#define G_APPEND_AXIS  "GCTAResponseTable::append_axis(std::vector<double>&,"\
                         " std::vector<double>&, std::string&, std::string&)"
#define G_APPEND_TABLE        "GCTAResponseTable::append_table(std::string&,"\
                                                             " std::string&)"
#define G_READ                         "GCTAResponseTable::read(GFitsTable&)"
#define G_READ_COLNAMES       "GCTAResponseTable::read_colnames(GFitsTable&)"
#define G_READ_AXES               "GCTAResponseTable::read_axes(GFitsTable&)"
#define G_READ_TABLES           "GCTAResponseTable::read_tables(GFitsTable&)"

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
 * Construct an empty CTA response table.
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
 * Construct a CTA response table by copying information from an existing
 * response table. A deep copy is performed for construction, so that the
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
 * Construct a CTA response table from the information in a FITS table HDU.
 * Refer to the read() method for details.
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
 * Assigns a CTA response table to another object. The method is assigning
 * a deep copy of the response table, so that the original object can be
 * destroyed after assignment without any loss of information.
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
 * @return Vector of linearly interpolated response tables.
 *
 * @exception GException::invalid_value
 *            Invalid response table dimension.
 *
 * Evaluates all response tables at a given value. This method applies to
 * one-dimensional tables only. If the tables have a different dimension
 * an exception is thrown.
 *
 * The evaluation is performed by a linear interpolation of the table. If
 * the specified value lies outside the range covered by the table, the
 * table is linearly extrapolated from using either the first or the last
 * two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
std::vector<double> GCTAResponseTable::operator()(const double& arg) const
{
    // Throw exception if table is not 1D
    if (axes() != 1) {
        std::string msg = "Invalid response table dimension "+
                          gammalib::str(axes())+" encountered. Response "
                          "table needs to be one-dimensional.";
        throw GException::invalid_value(G_OPERATOR1, msg);
    }

    // Initialise result vector
    std::vector<double> result(m_ntables);

    // Set indices and weighting factors for interpolation
    update(arg);

    // Perform 1D interpolation
    for (int i = 0; i < m_ntables; ++i) {
        result[i] = m_wgt_left  * m_tables[i][m_inx_left] +
                    m_wgt_right * m_tables[i][m_inx_right];
    }

    // Return result vector
    return result;
}


/***********************************************************************//**
 * @brief Bilinear interpolation operator for 2D tables
 *
 * @param[in] arg1 Value for first axis.
 * @param[in] arg2 Value for second axis.
 * @return Vector of bilinearly interpolated response tables.
 *
 * @exception GException::invalid_value
 *            Invalid response table dimension.
 *
 * Evaluates all response tables at a given value. This method applies to
 * two-dimensional tables only. If the tables have a different dimension
 * an exception is thrown.
 *
 * The evaluation is performed by a bilinear interpolation of the table. If
 * the specified value lies outside the range covered by the table, the
 * table is linearly extrapolated from using either the first or the last
 * two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
std::vector<double> GCTAResponseTable::operator()(const double& arg1,
                                                  const double& arg2) const
{
    // Throw exception if table is not 2D
    if (axes() != 2) {
        std::string msg = "Invalid response table dimension "+
                          gammalib::str(axes())+" encountered. Response "
                          "table needs to be two-dimensional.";
        throw GException::invalid_value(G_OPERATOR2, msg);
    }

    // Initialise result vector
    std::vector<double> result(m_ntables);

    // Set indices and weighting factors for interpolation
    update(arg1, arg2);

    // Perform 2D interpolation
    for (int i = 0; i < m_ntables; ++i) {
        result[i] = m_wgt1 * m_tables[i][m_inx1] +
                    m_wgt2 * m_tables[i][m_inx2] +
                    m_wgt3 * m_tables[i][m_inx3] +
                    m_wgt4 * m_tables[i][m_inx4];
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
 * @return Vector of trilinearly interpolated response tables.
 *
 * @exception GException::invalid_value
 *            Invalid response table dimension.
 *
 * Evaluates all response tables at a given value. This method applies to
 * three-dimensional tables only. If the tables have a different dimension
 * an exception is thrown.
 *
 * The evaluation is performed by a trilinear interpolation of the table. If
 * the specified value lies outside the range covered by the table, the
 * table is linearly extrapolated from using either the first or the last
 * two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
std::vector<double> GCTAResponseTable::operator()(const double& arg1,
                                                  const double& arg2,
                                                  const double& arg3) const
{
    // Throw exception if table is not 3D
    if (axes() != 3) {
        std::string msg = "Invalid response table dimension "+
                          gammalib::str(axes())+" encountered. Response "
                          "table needs to be tri-dimensional.";
        throw GException::invalid_value(G_OPERATOR3, msg);
    }

    // Initialise result vector
    std::vector<double> result(m_ntables);

    // Set indices and weighting factors for interpolation
    update(arg1, arg2, arg3);

    // Perform 3D interpolation
    for (int i = 0; i < m_ntables; ++i) {
        result[i] = m_wgt1 * m_tables[i][m_inx1] +
                    m_wgt2 * m_tables[i][m_inx2] +
                    m_wgt3 * m_tables[i][m_inx3] +
                    m_wgt4 * m_tables[i][m_inx4] +
                    m_wgt5 * m_tables[i][m_inx5] +
                    m_wgt6 * m_tables[i][m_inx6] +
                    m_wgt7 * m_tables[i][m_inx7] +
                    m_wgt8 * m_tables[i][m_inx8] ;
    }

    // Return result vector
    return result;
}


/***********************************************************************//**
 * @brief Table element access operator
 *
 * @param[in] element Table element index [0,...,elements()-1].
 * @return Table element value.
 *
 * @exception GException::out_of_range
 *            Invalid element index.
 *
 * Returns the value of a table element for the first response table.
 ***************************************************************************/
double& GCTAResponseTable::operator()(const int& element)
{
    // Optionally check if the index and element is valid
    #if defined(G_RANGE_CHECK)
    if (element < 0 || element >= elements()) {
        throw GException::out_of_range(G_ELEMENT_OPERATOR1, "Element index",
                                       element, elements());
    }
    #endif

    // Return elements
    return (m_tables[0][element]);
}


/***********************************************************************//**
 * @brief Table element access operator (const version)
 *
 * @param[in] element Table element index [0,...,elements()-1].
 * @return Table element value.
 *
 * @exception GException::out_of_range
 *            Invalid element index.
 *
 * Returns the value of a table element for the first response table.
 ***************************************************************************/
const double& GCTAResponseTable::operator()(const int& element) const
{
    // Optionally check if the index and element is valid
    #if defined(G_RANGE_CHECK)
    if (element < 0 || element >= elements()) {
        throw GException::out_of_range(G_ELEMENT_OPERATOR1, "Element index",
                                       element, elements());
    }
    #endif

    // Return elements
    return (m_tables[0][element]);
}


/***********************************************************************//**
 * @brief Table element access operator
 *
 * @param[in] table Table index [0,...,tables()-1].
 * @param[in] element Element index [0,...,elements()-1].
 * @return Response table element.
 *
 * @exception GException::out_of_range
 *            Invalid table or element index.
 *
 * Returns the value of a table element for the response table with index
 * @p table.
 ***************************************************************************/
double& GCTAResponseTable::operator()(const int& table, const int& element)
{
    // Optionally check if the index and element is valid
    #if defined(G_RANGE_CHECK)
    if (table < 0 || table >= tables()) {
        throw GException::out_of_range(G_ELEMENT_OPERATOR2, "Table index",
                                       table, tables());
    }
    if (element < 0 || element >= elements()) {
        throw GException::out_of_range(G_ELEMENT_OPERATOR2, "Element index",
                                       element, elements());
    }
    #endif

    // Return elements
    return (m_tables[table][element]);
}


/***********************************************************************//**
 * @brief Table element access operator (const version)
 *
 * @param[in] table Table index [0,...,tables()-1].
 * @param[in] element Element index [0,...,elements()-1].
 * @return Response table element.
 *
 * @exception GException::out_of_range
 *            Invalid table or element index.
 *
 * Returns the value of a table element for the response table with index
 * @p table.
 ***************************************************************************/
const double& GCTAResponseTable::operator()(const int& table, const int& element) const
{
    // Optionally check if the index and element is valid
    #if defined(G_RANGE_CHECK)
    if (table < 0 || table >= tables()) {
        throw GException::out_of_range(G_ELEMENT_OPERATOR2, "Table index",
                                       table, tables());
    }
    if (element < 0 || element >= elements()) {
        throw GException::out_of_range(G_ELEMENT_OPERATOR2, "Element index",
                                       element, elements());
    }
    #endif

    // Return elements
    return (m_tables[table][element]);
}


/***********************************************************************//**
 * @brief Linear interpolation operator for 1D tables
 *
 * @param[in] table Table index [0,...,tables()-1].
 * @param[in] arg Value.
 * @return Linearly interpolated response table.
 *
 * @exception GException::invalid_value
 *            Invalid response table dimension.
 * @exception GException::out_of_range
 *            Invalid table index.
 *
 * Evaluates one response table at a given value. This method applies to
 * one-dimensional response tables. If the table has a different dimension
 * an exception is thrown.
 *
 * The evaluation is performed by a linear interpolation of the table values.
 * If the specified value lies outside the range covered by the table, the
 * table is linearly extrapolated from using either the first or the last
 * two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
double GCTAResponseTable::operator()(const int& table, const double& arg) const
{
    // Throw exception if table is not 1D
    if (axes() != 1) {
        std::string msg = "Invalid response table dimension "+
                          gammalib::str(axes())+" encountered. Response "
                          "table needs to be one-dimensional.";
        throw GException::invalid_value(G_INX_OPERATOR1, msg);
    }

    // Throw exception if index is not valid
    #if defined(G_RANGE_CHECK)
    if (table < 0 || table >= tables()) {
        throw GException::out_of_range(G_INX_OPERATOR1, "Table index",
                                       table, tables());
    }
    #endif

    // Set indices and weighting factors for interpolation
    update(arg);

    // Perform 1D interpolation
    double result = m_wgt_left  * m_tables[table][m_inx_left] +
                    m_wgt_right * m_tables[table][m_inx_right];

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Bilinear interpolation operator for 2D tables
 *
 * @param[in] table Table index [0,...,tables()-1].
 * @param[in] arg1 Value for first axis.
 * @param[in] arg2 Value for second axis.
 * @return Bilinearly interpolated response table.
 *
 * @exception GException::invalid_value
 *            Invalid response table dimension.
 * @exception GException::out_of_range
 *            Invalid table index.
 *
 * Evaluates one response table at a given value. This method applies to
 * two-dimensional response tables. If the table has a different dimension
 * an exception is thrown.
 *
 * The evaluation is performed by a bilinear interpolation of the table
 * values. If the specified value lies outside the range covered by the
 * table, the table is linearly extrapolated from using either the first
 * or the last two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
double GCTAResponseTable::operator()(const int&    table,
                                     const double& arg1,
                                     const double& arg2) const
{
    // Throw exception if table is not 2D
    if (axes() != 2) {
        std::string msg = "Invalid response table dimension "+
                          gammalib::str(axes())+" encountered. Response "
                          "table needs to be two-dimensional.";
        throw GException::invalid_value(G_INX_OPERATOR2, msg);
    }

    // Throw exception if index is not valid
    #if defined(G_RANGE_CHECK)
    if (table < 0 || table >= tables()) {
        throw GException::out_of_range(G_INX_OPERATOR2, "Table index",
                                       table, tables());
    }
    #endif

    // Set indices and weighting factors for interpolation
    update(arg1, arg2);

    // Perform 2D interpolation
    double result = m_wgt1 * m_tables[table][m_inx1] +
                    m_wgt2 * m_tables[table][m_inx2] +
                    m_wgt3 * m_tables[table][m_inx3] +
                    m_wgt4 * m_tables[table][m_inx4];

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Trilinear interpolation operator for 3D tables
 *
 * @param[in] table Table index [0,...,tables()-1].
 * @param[in] arg1 Value for first axis.
 * @param[in] arg2 Value for second axis.
 * @param[in] arg3 Value for second axis.
 * @return Trilinearly interpolated response table.
 *
 * @exception GException::invalid_value
 *            Invalid response table dimension.
 * @exception GException::out_of_range
 *            Invalid table index.
 *
 * Evaluates one response table at a given value. This method applies to
 * three-dimensional response tables. If the table has a different dimension
 * an exception is thrown.
 *
 * The evaluation is performed by a trilinear interpolation of the table
 * values. If the specified value lies outside the range covered by the
 * table, the table is linearly extrapolated from using either the first
 * or the last two vector elements.
 *
 * @todo Write down formula.
 ***************************************************************************/
double GCTAResponseTable::operator()(const int&    table,
                                     const double& arg1, 
                                     const double& arg2,
                                     const double& arg3) const
{
    // Throw exception if table is not 3D
    if (axes() != 3) {
        std::string msg = "Invalid response table dimension "+
                          gammalib::str(axes())+" encountered. Response "
                          "table needs to be three-dimensional.";
        throw GException::invalid_value(G_INX_OPERATOR3, msg);
    }

    // Throw exception if index is not valid
    #if defined(G_RANGE_CHECK)
    if (table < 0 || table >= tables()) {
        throw GException::out_of_range(G_INX_OPERATOR3, "Table index",
                                       table, tables());
    }
    #endif

    // Set indices and weighting factors for interpolation
    update(arg1, arg2, arg3);

    // Perform 3D interpolation
    double result = m_wgt1 * m_tables[table][m_inx1] +
                    m_wgt2 * m_tables[table][m_inx2] +
                    m_wgt3 * m_tables[table][m_inx3] +
                    m_wgt4 * m_tables[table][m_inx4] +
                    m_wgt5 * m_tables[table][m_inx5] +
                    m_wgt6 * m_tables[table][m_inx6] +
                    m_wgt7 * m_tables[table][m_inx7] +
                    m_wgt8 * m_tables[table][m_inx8];

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear response table
 *
 * Clears the response table.
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
 * @brief Clone response table
 *
 * @return Deep copy of response table instance.
 *
 * Returns a pointer to a deep copy of the response table.
 ***************************************************************************/
GCTAResponseTable* GCTAResponseTable::clone(void) const
{
    return new GCTAResponseTable(*this);
}


/***********************************************************************//**
 * @brief Check whether a table exists
 *
 * @param[in] name Column name
 * @return True if table exists.
 *
 * Checks whether a table exists.
 ***************************************************************************/
bool GCTAResponseTable::has_table(const std::string& name) const
{
    // Initialise existence flag
    bool exists = false;

    // Loop over all tables to find index
    int table = 0;
    for (; table < tables(); ++table) {
        if (m_colname_table[table] == name) {
            exists = true;
            break;
        }
    }

    // Return existence flag
    return exists;
}


/***********************************************************************//**
 * @brief Check whether an axis exists
 *
 * @param[in] name Column name
 * @return True if axis exists.
 *
 * Check whether an axis exists.
 ***************************************************************************/
bool GCTAResponseTable::has_axis(const std::string& name) const
{
    // Initialise existence flag
    bool exists = false;

    // Build column names
    std::string col_lo = name + "_LO";
    std::string col_hi = name + "_HI";

    // Loop over all axes to find index
    int axis = 0;
    for (; axis < axes(); ++axis) {
        if ((axis_lo_name(axis) == col_lo) &&
            (axis_hi_name(axis) == col_hi)) {
            exists = true;
            break;
        }
    }

    // Return existence flag
    return exists;
}


/***********************************************************************//**
 * @brief Determine index of table
 *
 * @param[in] name Column name
 * @return Table index.
 *
 * @exception GException::invalid_value
 *            Table not found.
 *
 * Determines the index of a table @p name. An exception is thrown if the
 * table is not found.
 ***************************************************************************/
int GCTAResponseTable::table(const std::string& name) const
{
    // Loop over all tables to find index
    int table = 0;
    for (; table < tables(); ++table) {
        if (m_colname_table[table] == name) {
            break;
        }
    }

    // Throw an exception if table has not been found
    if (table >= tables()) {
        std::string msg = "Table \""+name+"\" not found. Please verify the "
                          "response table.";
        throw GException::invalid_value(G_TABLE, msg);
    }

    // Return table index
    return table;
}


/***********************************************************************//**
 * @brief Return table unit
 *
 * @param[in] table Table index [0,...,tables()-1].
 * @return Unit of table.
 *
 * @exception GException::out_of_range
 *            Table index out of range.
 *
 * Returns the unit of the table.
 ***************************************************************************/
const std::string& GCTAResponseTable::unit(const int& table) const
{
    // Optionally check if the table index is valid
    #if defined(G_RANGE_CHECK)
    if (table < 0 || table >= tables()) {
        throw GException::out_of_range(G_UNIT, "Table index", table, tables());
    }
    #endif

    // Return units
    return (m_units_table[table]);
}


/***********************************************************************//**
 * @brief Scale table
 *
 * @param[in] table Table index [0,...,tables()-1].
 * @param[in] scale Scaling factor.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Multiplies all values in a table by the scaling factor @p scale.
 ***************************************************************************/
void GCTAResponseTable::scale(const int& table, const double& scale)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (table < 0 || table >= tables()) {
        throw GException::out_of_range(G_SCALE, "Table index", table,
                                       tables());
    }
    #endif

    // Scale table values
    for (int i = 0; i < m_nelements; ++i) {
        m_tables[table][i] *= scale;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Determine index of an axis
 *
 * @param[in] name Column name
 * @return Axis index.
 *
 * @exception GException::invalid_value
 *            Axis definition columns not found.
 *
 * Determines the index of an axis @p name, where the lower and upper bin
 * edges of the axis are stored in the columns @p name_LO and @p name_HI.
 * An exception is thrown if the axis is not found.
 ***************************************************************************/
int GCTAResponseTable::axis(const std::string& name) const
{
    // Build column names
    std::string col_lo = name + "_LO";
    std::string col_hi = name + "_HI";

    // Loop over all axes to find index
    int axis = 0;
    for (; axis < axes(); ++axis) {
        if ((axis_lo_name(axis) == col_lo) &&
            (axis_hi_name(axis) == col_hi)) {
            break;
        }
    }

    // Throw an exception if axis has not been found
    if (axis >= axes()) {
        std::string msg = "Axis definition columns for axis \""+name+"\" not "
                          "found. Please verify the axis names in the "
                          "response table.";
        throw GException::invalid_value(G_AXIS, msg);
    }

    // Return axis index
    return axis;
}


/***********************************************************************//**
 * @brief Return number bins in an axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
 * @return Number of bins along the specified axis.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the number of bins in the specified @p axis.
 ***************************************************************************/
int GCTAResponseTable::axis_bins(const int& axis) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_BINS, "Axis index", axis, axes());
    }
    #endif

    // Return axis length
    return (m_axis_lo[axis].size());
}


/***********************************************************************//**
 * @brief Return lower bin boundary for bin in axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
 * @param[in] bin Bin index [0,...,axis_bins(axis)-1].
 * @return Lower bin boundary.
 *
 * @exception GException::out_of_range
 *            Axis or bin index out of range.
 *
 * Returns the lower boundary for a given @p bin of the specified @p axis.
 ***************************************************************************/
const double& GCTAResponseTable::axis_lo(const int& axis, const int& bin) const
{
    // Check if index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_LO, "Axis index", axis, axes());
    }
    if (bin < 0 || bin >= m_axis_lo[axis].size()) {
        throw GException::out_of_range(G_AXIS_LO, "Bin index", bin, 
                                       m_axis_lo[axis].size());
    }
    #endif

    // Return bin boundary
    return (m_axis_lo[axis][bin]);
}


/***********************************************************************//**
 * @brief Return upper bin boundary for bin in axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
 * @param[in] bin Bin index [0,...,axis_bins(axis)-1].
 * @return Upper bin boundary.
 *
 * @exception GException::out_of_range
 *            Axis or bin index out of range.
 *
 * Returns the upper boundary for a given @p bin of the specified @p axis.
 ***************************************************************************/
const double& GCTAResponseTable::axis_hi(const int& axis, const int& bin) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_HI, "Axis index", axis, axes());
    }
    if (bin < 0 || bin >= m_axis_hi[axis].size()) {
        throw GException::out_of_range(G_AXIS_HI, "Bin index", bin,
                                       m_axis_hi[axis].size());
    }
    #endif

    // Return bin boundary
    return (m_axis_hi[axis][bin]);
}


/***********************************************************************//**
 * @brief Return lower bin boundary FITS table column name for axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
 * @return FITS table column name of lower bin boundary.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the FITS table column name of the lower boundary for the specified
 * @p axis.
 ***************************************************************************/
const std::string& GCTAResponseTable::axis_lo_name(const int& axis) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_LO_NAME, "Axis index", axis, axes());
    }
    #endif

    // Return units
    return (m_colname_lo[axis]);
}


/***********************************************************************//**
 * @brief Return upper bin boundary FITS table column name for axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
 * @return FITS table column name of upper bin boundary.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the FITS table column name of the upper boundary for the specified
 * @p axis.
 ***************************************************************************/
const std::string& GCTAResponseTable::axis_hi_name(const int& axis) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_HI_NAME, "Axis index", axis, axes());
    }
    #endif

    // Return units
    return (m_colname_hi[axis]);
}


/***********************************************************************//**
 * @brief Return lower bin boundary unit for axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
 * @return Unit of lower bin boundary.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the unit of the lower boundary for the specified @p axis.
 ***************************************************************************/
const std::string& GCTAResponseTable::axis_lo_unit(const int& axis) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_LO_UNIT, "Axis index", axis, axes());
    }
    #endif

    // Return units
    return (m_units_lo[axis]);
}


/***********************************************************************//**
 * @brief Return upper bin boundary unit for axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
 * @return Unit of upper bin boundary.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the unit of the upper boundary for the specified @p axis.
 ***************************************************************************/
const std::string& GCTAResponseTable::axis_hi_unit(const int& axis) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_HI_UNIT, "Axis index", axis, axes());
    }
    #endif

    // Return units
    return (m_units_hi[axis]);
}


/***********************************************************************//**
 * @brief Set nodes for a linear axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
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
void GCTAResponseTable::axis_linear(const int& axis)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_LINEAR, "Axis index", axis,
                                       axes());
    }
    #endif

    // Get number of bins
    int bins = axis_bins(axis);

    // Allocate node values
    std::vector<double> axis_nodes(bins);

    // Compute nodes
    for (int i = 0; i < bins; ++i) {
        axis_nodes[i] = 0.5*(m_axis_lo[axis][i] + m_axis_hi[axis][i]);
    }

    // Set node array
    m_axis_nodes[axis] = GNodeArray(axis_nodes);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for a logarithmic (base 10) axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
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
void GCTAResponseTable::axis_log10(const int& axis)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_LOG10, "Axis index", axis,
                                       axes());
    }
    #endif

    // Get number of bins
    int bins = axis_bins(axis);

    // Allocate node values
    std::vector<double> axis_nodes(bins);

    // Compute nodes
    for (int i = 0; i < bins; ++i) {
        axis_nodes[i] = std::log10(std::sqrt(m_axis_lo[axis][i] * 
                                             m_axis_hi[axis][i]));
    }

    // Set node array
    m_axis_nodes[axis] = GNodeArray(axis_nodes);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for a radians axis
 *
 * @param[in] axis Axis index [0,...,axes()-1].
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
void GCTAResponseTable::axis_radians(const int& axis)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_RADIANS, "Axis index", axis,
                                       axes());
    }
    #endif

    // Get number of bins
    int bins = axis_bins(axis);

    // Allocate node values
    std::vector<double> axis_nodes(bins);

    // Compute nodes
    for (int i = 0; i < bins; ++i) {
        axis_nodes[i] = 0.5*(m_axis_lo[axis][i] + m_axis_hi[axis][i]) *
                        gammalib::deg2rad;
    }

    // Set node array
    m_axis_nodes[axis] = GNodeArray(axis_nodes);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append an axis to the response table
 * 
 * @param[in] axis_lo Lower axis boundaries.
 * @param[in] axis_hi Upper axis boundaries. 
 * @param[in] name Axis name. 
 * @param[in] unit Axis unit.
 *
 * @exception GException::invalid_argument
 *            Incompatible axis boundaries.
 *
 * Append an axis to the response table.
 ***************************************************************************/
void GCTAResponseTable::append_axis(const std::vector<double>& axis_lo, 
                                    const std::vector<double>& axis_hi,
                                    const std::string&         name,
                                    const std::string&         unit)
{
    // Throw an exception if the length of axis_lo and axis_hi are different
    if (axis_lo.size() != axis_hi.size()) {
        std::string msg = "Number of elements in lower axis boundaries "
                          "mismatches the number of elements in upper axis "
                          "boundaries ("+gammalib::str(axis_lo.size())+
                          " != "+gammalib::str(axis_hi.size())+").";
        throw GException::invalid_argument(G_APPEND_AXIS, msg);
    }

    // Set axis names
    std::string name_lo = name + "_LO";
    std::string name_hi = name + "_HI";

    // Append axis  
    m_colname_lo.push_back(name_lo);
    m_colname_hi.push_back(name_hi);
    m_axis_lo.push_back(axis_lo);
    m_axis_hi.push_back(axis_hi);
    m_units_lo.push_back(unit);
    m_units_hi.push_back(unit);


    // Set node array
    std::vector<double> axis_nodes(axis_lo.size());
    for (int k = 0; k < axis_lo.size(); ++k) {
        axis_nodes[k] = 0.5*(axis_lo[k] + axis_hi[k]);
    }
    m_axis_nodes.push_back(GNodeArray(axis_nodes));

    // Increment number of axes
    m_naxes++;

    // Compute the cube size
    m_nelements = axis_bins(0);
    for (int i = 1; i < axes(); ++i) {
        m_nelements *= axis_bins(i);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append table to response table
 * 
 * @param[in] name Table name. 
 * @param[in] unit Table unit.
 *
 * @exception GException::invalid_value
 *            No axes have been specified.
 *
 * Append a table to the response table. The number of elements in the table
 * is the product of the length of all axes. All elements are set to zero by
 * default.
 ***************************************************************************/
void GCTAResponseTable::append_table(const std::string& name,
                                     const std::string& unit)
{
    // Throw an exception message when m_nelement is zero
    if (m_nelements == 0) {
        std::string msg = "No axis columns have been defined. Please append "
                          "axis columns before appending table columns.";
        throw GException::invalid_value(G_APPEND_TABLE, msg);
    }

    // Append table column name and unit
    m_colname_table.push_back(name);
    m_units_table.push_back(unit);
    
    // Initialise empty table column
    std::vector<double> table(m_nelements, 0.0);
    
    // Append column
    m_tables.push_back(table);

    // Increment number of table columns
    m_ntables++;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return axis nodes
 *
 * @param[in] axis Axis index [0,...,axes()-1].
 * @return Node array for axis.
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the node array of the specified @p axis.
 ***************************************************************************/
const GNodeArray& GCTAResponseTable::axis_nodes(const int& axis) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= axes()) {
        throw GException::out_of_range(G_AXIS_NODES, "Axis index", axis,
                                       axes());
    }
    #endif

    // Return node array
    return (m_axis_nodes[axis]);
}


/***********************************************************************//**
 * @brief Read response table from FITS table HDU
 *
 * @param[in] table FITS table.
 *
 * Reads CTA response table information from a FITS table. The FITS table
 * is expected to have a single row, and axes and table information are
 * found in vector columns. Axes information is expected to be placed
 * before table information. 
 *
 * Each axis is defined by two vector columns of equal width, describing the
 * lower and upper limits for each bin. The column names for this boundary
 * information terminates by "_LO" and "HI", respectively (upper case). It
 * is expected that the "_LO" column preceeds the "_HI" column.
 *
 * Following the axes columns are table columns of equal width. The width
 * of each table column is given by the product of the lengths of all
 * axes. It is furthermore expected that the first axis is the most rapidely
 * varying index of the vector.
 ***************************************************************************/
void GCTAResponseTable::read(const GFitsTable& table)
{
    // Clear instance
    clear();

    // Read column names
    read_colnames(table);

    // Read axes
    read_axes(table);

    // Read tables
    read_tables(table);

    // Read TELESCOP keyword
    m_telescope = (table.has_card("TELESCOP")) ? table.string("TELESCOP") : "";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write response table into FITS table HDU
 *
 * @param[in] table FITS table.
 *
 * Writes the response table in a FITS table HDU.
 *
 * Axes columns are written in floating point columns with the suffixes "_LO"
 * and "_HI". After the axes columns the response tables are written in
 * floating point columns.
 ***************************************************************************/
void GCTAResponseTable::write(GFitsTable& table) const
{
    // Initialise dimension vector
    std::vector<int> dim;

    // Loop over all response table axes
    for (int iaxis = 0; iaxis < m_naxes; ++iaxis) {

        // Create axis columns
        GFitsTableFloatCol col_lo(m_colname_lo[iaxis], 1,
                                  m_axis_lo[iaxis].size());
        GFitsTableFloatCol col_hi(m_colname_hi[iaxis], 1,
                                  m_axis_hi[iaxis].size());

        // Loop through all elements in this axis column
        for (int i = 0; i < m_axis_lo[iaxis].size(); ++i) {
            col_lo(0,i) = m_axis_lo[iaxis][i];
            col_hi(0,i) = m_axis_hi[iaxis][i];
        }

        // Set column units
        col_lo.unit(m_units_lo[iaxis]);
        col_hi.unit(m_units_hi[iaxis]);

        // Append column to FITS table
        table.append(col_lo);
        table.append(col_hi);

        // Append dimension
        dim.push_back(axis_bins(iaxis));

    } // endif: looped over all axes in response table

    // Loop over all tables
    for (int itable = 0; itable < m_ntables; ++itable) {

        // Create table column
        GFitsTableFloatCol col_table(m_colname_table[itable], 1,
                                     m_tables[itable].size());
        
        // Loop through elements in this table column
        for (int i = 0; i < m_tables[itable].size() ; ++i) {
            col_table(0,i) = m_tables[itable][i];
        }

        // Set column unit
        col_table.unit(m_units_table[itable]);

        // Set column dimension
        col_table.dim(dim);

        // Append column to table
        table.append(col_table);

    } // endfor: looped over all axes in response table

    // Write TELESCOPE keyword
    table.card("TELESCOP", m_telescope, "Telescope");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print response table information
 *
 * @param[in] chatter Chattiness.
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
        result.append("\n"+gammalib::parformat("Telescope") +
                      m_telescope);
        result.append("\n"+gammalib::parformat("Dimension") +
                      gammalib::str(axes()));

        // Append axes information
        for (int i = 0; i < axes(); ++i) {
            result.append("\n"+gammalib::parformat("Axis " +
                          gammalib::str(i)));
            result.append(gammalib::str(axis_bins(i)));
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

        // Append table information
        for (int i = 0; i < tables(); ++i) {
            result.append("\n"+gammalib::parformat("Table " +
                          gammalib::str(i)));
            result.append(m_colname_table[i]);
            if (!m_units_table[i].empty()) {
                result.append(" ["+m_units_table[i]+"]");
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
 * @brief Initialise response table members
 *
 * Initialises all members of the response table.
 ***************************************************************************/
void GCTAResponseTable::init_members(void)
{
    // Initialise members
    m_naxes     = 0;
    m_ntables   = 0;
    m_nelements = 0;
    m_colname_lo.clear();
    m_colname_hi.clear();
    m_colname_table.clear();
    m_axis_lo.clear();
    m_axis_hi.clear();
    m_units_lo.clear();
    m_units_hi.clear();
    m_units_table.clear();
    m_axis_nodes.clear();
    m_tables.clear();
    m_telescope.clear();

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
 * @brief Copy response table members
 *
 * @param[in] table Response table.
 *
 * Copies all response table members.
 ***************************************************************************/
void GCTAResponseTable::copy_members(const GCTAResponseTable& table)
{
    // Copy number of bins
    m_naxes         = table.m_naxes;
    m_ntables       = table.m_ntables;
    m_nelements     = table.m_nelements;
    m_colname_lo    = table.m_colname_lo;
    m_colname_hi    = table.m_colname_hi;
    m_colname_table = table.m_colname_table;
    m_axis_lo       = table.m_axis_lo;
    m_axis_hi       = table.m_axis_hi;
    m_units_lo      = table.m_units_lo;
    m_units_hi      = table.m_units_hi;
    m_units_table   = table.m_units_table;
    m_axis_nodes    = table.m_axis_nodes;
    m_tables        = table.m_tables;
    m_telescope     = table.m_telescope;

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
 * @brief Delete response table members
 *
 * De-allocates any memory that was allocated by the response table.
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
 * @exception GException::invalid_value
 *            Invalid response table format encountered.
 *
 * Read the response table column names from the HDU. Column names
 * terminating with "_LO" and "_HI" define the table axes, while all other
 * column names define response tables. It is assumed that table axes are
 * given in subsequent order, with the first column corresponding to the
 * lower bin boundaries of the first axis. Lower bin boundaries are indicated
 * by the "_LO" termination. Upper bin boundaries are assumed to follow
 * immediately the lower bin boundaires and are designated by the "_HI"
 * termination.
 *
 * This method sets the following members:
 *    m_colname_lo - Column names of lower boundaries
 *    m_colname_hi - Column names of upper boundaries
 *    m_colname_table - Column names of tables
 *    m_naxes - Number of axes
 *    m_ntables - Number of tables
 *
 * @todo Implement exceptions for invalid HDU format
 ***************************************************************************/
void GCTAResponseTable::read_colnames(const GFitsTable& hdu)
{
    // Clear column name arrays
    m_naxes   = 0;
    m_ntables = 0;
    m_colname_lo.clear();
    m_colname_hi.clear();
    m_colname_table.clear();

    // Initialise search mode. There are three search modes:
    // 0 - we're looking for the next axis by searching for a column
    //     terminating with "_LO"
    // 1 - we're looking for the upper boundary of an axis, terminating
    //     with "_HI"
    // 2 - we're looking for a table column
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
        // mode to 2 since from now on we should only have table
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
                    std::string msg = "Column \""+colname+"\" encountered "
                                      "without preceeding \"_LO\" column. "
                                      "Please correct response table format.";
                    throw GException::invalid_value(G_READ_COLNAMES, msg);
                }
                else {
                    mode = 2;
                    m_colname_table.push_back(colname);
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
                std::string msg = "Expected column \""+next_column+"\" not "
                                  "found. The axis upper bin boundary "
                                  "column has to follow immediately the "
                                  "lower bin boundary column.";
                throw GException::invalid_value(G_READ_COLNAMES, msg);
            }
        }

        // If we search for a table column, make sure that we have
        // neither a "_LO" nor a "_HI" column.
        else {
            if (colname.rfind("_LO") != std::string::npos) {
                std::string msg = "Column \""+colname+"\" encountered while "
                                  "searching for table columns. All "
                                  "lower bin boundary columns have to be "
                                  "placed before the table columns.";
                throw GException::invalid_value(G_READ_COLNAMES, msg);
            }
            else if (colname.rfind("_HI") != std::string::npos) {
                std::string msg = "Column \""+colname+"\" encountered while "
                                  "searching for table columns. All "
                                  "upper bin boundary columns have to be "
                                  "placed before the table columns.";
                throw GException::invalid_value(G_READ_COLNAMES, msg);
            }
            else {
                m_colname_table.push_back(colname);
            }
        }

    } // endfor: looped over all columns

    // Store number of axes
    m_naxes = m_colname_lo.size();

    // Store number of tables
    m_ntables = m_colname_table.size();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read axes definitions from FITS HDU
 *
 * @param[in] hdu Response table HDU.
 *
 * @exception GException::invalid_value
 *            Incompatible axes columns encoutered.
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
 *     m_axis_lo - Axes lower boundaries
 *     m_axis_hi - Axes upper boundaries
 *     m_axis_nodes - Axes mean values
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
            std::string msg = "Axis lower bin boundary column \""+
                              m_colname_lo[i]+"\" contains "+
                              gammalib::str(num)+" elements while upper "
                              "bin boundary column \""+m_colname_hi[i]+
                              "\" contains "+gammalib::str(col_hi->number())+
                              " elements. Both columns need to contain the "
                              "same number of elements.";
            throw GException::invalid_value(G_READ_AXES, msg);
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
 * @brief Read tables
 *
 * @param[in] hdu Response table HDU.
 *
 * @exception GException::invalid_value
 *            Table vector of bad size encountered.
 *
 * Reads the tables from the FITS HDU.
 *
 * The method also sets m_nelements which gives the number of elements in a
 * response table.
 ***************************************************************************/
void GCTAResponseTable::read_tables(const GFitsTable& hdu)
{
    // Clear tables
    m_tables.clear();

    // Compute expected cube size
    m_nelements = (axes() > 0) ? axis_bins(0) : 0;
    for (int i = 1; i < axes(); ++i) {
        m_nelements *= axis_bins(i);
    }

    // Loop over all tables
    for (int i = 0; i < tables(); ++i) {

        // Get pointer to table column
        const GFitsTableCol* col = hdu[m_colname_table[i]];

        // Extract number of bins. Verify that the number of bins
        // corresponds to the expectation.
        int num = col->number();
        if (num != m_nelements) {
            std::string msg = "Table column \""+m_colname_table[i]+"\" "
                              "contains "+gammalib::str(num)+" elements "
                              "while "+gammalib::str(m_nelements)+" "
                              "elements were expected from the axis "
                              "definitions. Please check the response "
                              "table for consistency.";
            throw GException::invalid_value(G_READ_TABLES, msg);
        }

        // Check table dimension if there is dimension information
        if (col->dim().size() > 0) {

            // Check that the table dimension is identical to the number of
            // axis pairs
            if (axes() != col->dim().size()) {
                std::string msg = "Table column \""+m_colname_table[i]+"\" "
                                  "dimension "+
                                  gammalib::str(col->dim().size())+" is "
                                  "inconsistent with "+
                                  gammalib::str(axes())+" axis pairs. Please "
                                  "check the response table for consistency.";
                throw GException::invalid_value(G_READ_TABLES, msg);
            }

            // Check for all axes that the length of the axis columns is
            // identical to the size of the table in that axis direction
            for (int i = 0; i < axes(); ++i) {
                if (axis_bins(i) != col->dim()[i]) {
                    std::string msg = "Table column \""+m_colname_table[i]+"\" "
                                      "has "+gammalib::str(col->dim()[i])+
                                      " bins in axis "+gammalib::str(i)+
                                      " while corresponding axis columns \""+
                                      m_colname_lo[i]+"\" and \""+
                                      m_colname_hi[i]+"\" have a length of "+
                                      gammalib::str(axis_bins(i))+". Please "
                                      "check the response table for "
                                      "consistency.";
                    throw GException::invalid_value(G_READ_TABLES, msg);
                }
            }

        } // endif: table had dimension information


        // Initialise table
        std::vector<double> table(num);

        // Copy table values
        for (int k = 0; k < num; ++k) {
            table[k] = col->real(0,k);
        }

        // Push table into storage
        m_tables.push_back(table);

        // Push units on storage
        m_units_table.push_back(col->unit());

    } // endfor: looped over all tables

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
 ***************************************************************************/
void GCTAResponseTable::update(const double& arg) const
{
    // Get pointer to node array
    const GNodeArray* nodes = &(m_axis_nodes[0]);

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
 ***************************************************************************/
void GCTAResponseTable::update(const double& arg1, const double& arg2) const
{
    // Get pointers to node arrays
    const GNodeArray* nodes1 = &(m_axis_nodes[0]);
    const GNodeArray* nodes2 = &(m_axis_nodes[1]);

    // Set values for node arrays
    nodes1->set_value(arg1);
    nodes2->set_value(arg2);

    // Compute offsets
    int size1        = axis_bins(0);
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
 * @param[in] arg3 Argument for third axis.
 *
 * Updates the 3D interpolation cache. The interpolation cache is composed
 * of four indices and weights that define 4 data values of the 2D table
 * that are used for bilinear interpolation.
 *
 * @todo Write down formula
 ***************************************************************************/
void GCTAResponseTable::update(const double& arg1, const double& arg2,
                               const double& arg3) const
{

    // Get pointers to node arrays (circumvent const correctness)
    const GNodeArray* nodes1 = &(m_axis_nodes[0]);
    const GNodeArray* nodes2 = &(m_axis_nodes[1]);
    const GNodeArray* nodes3 = &(m_axis_nodes[2]);

    // Set values for node arrays
    nodes1->set_value(arg1);
    nodes2->set_value(arg2);
    nodes3->set_value(arg3);

    // Compute offsets
    int size1          = axis_bins(0);
    int size2          = axis_bins(1);
    int size12         = size1 * size2;
    int offset_left_2  = nodes2->inx_left()  * size1;
    int offset_right_2 = nodes2->inx_right() * size1;
    int offset_left_3  = nodes3->inx_left()  * size12;
    int offset_right_3 = nodes3->inx_right() * size12;

    // Pre-compute stuff
    int inx1l2l = nodes1->inx_left()  + offset_left_2;
    int inx1l2r = nodes1->inx_left()  + offset_right_2;
    int inx1r2l = nodes1->inx_right() + offset_left_2;
    int inx1r2r = nodes1->inx_right() + offset_right_2;

    // Set indices for tri-linear interpolation
    m_inx1 = inx1l2l + offset_left_3;
    m_inx2 = inx1l2l + offset_right_3;
    m_inx3 = inx1l2r + offset_left_3;
    m_inx4 = inx1l2r + offset_right_3;
    m_inx5 = inx1r2l + offset_left_3;
    m_inx6 = inx1r2l + offset_right_3;
    m_inx7 = inx1r2r + offset_left_3;
    m_inx8 = inx1r2r + offset_right_3;

    // Pre-compute stuff
    double wgt1l2l = nodes1->wgt_left()  * nodes2->wgt_left();
    double wgt1l2r = nodes1->wgt_left()  * nodes2->wgt_right();
    double wgt1r2l = nodes1->wgt_right() * nodes2->wgt_left();
    double wgt1r2r = nodes1->wgt_right() * nodes2->wgt_right();

    // Set weighting factors for tri-linear interpolation
    m_wgt1 = wgt1l2l * nodes3->wgt_left();
    m_wgt2 = wgt1l2l * nodes3->wgt_right();
    m_wgt3 = wgt1l2r * nodes3->wgt_left();
    m_wgt4 = wgt1l2r * nodes3->wgt_right();
    m_wgt5 = wgt1r2l * nodes3->wgt_left();
    m_wgt6 = wgt1r2l * nodes3->wgt_right();
    m_wgt7 = wgt1r2r * nodes3->wgt_left();
    m_wgt8 = wgt1r2r * nodes3->wgt_right();

    // Return
    return;
}
