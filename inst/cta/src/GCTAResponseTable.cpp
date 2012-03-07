/***************************************************************************
 *            GCTAResponseTable.cpp  -  CTA response table class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GException.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GCTAException.hpp"
#include "GCTAResponseTable.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR1                  "GCTAResponseTable::operator()(double&)"
#define G_OPERATOR2          "GCTAResponseTable::operator()(double&,double&)"
#define G_AXIS                                "GCTAResponseTable::axis(int&)"
#define G_AXIS_LO                     "GCTAResponseTable::axis_lo(int&,int&)"
#define G_AXIS_HI                     "GCTAResponseTable::axis_hi(int&,int&)"
#define G_AXIS_LINEAR                  "GCTAResponseTable::axis_linear(int&)"
#define G_AXIS_LOG10                    "GCTAResponseTable::axis_log10(int&)"
#define G_AXIS_RADIANS                "GCTAResponseTable::axis_radians(int&)"
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
GCTAResponseTable::GCTAResponseTable(const GFitsTable* hdu)
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
 *
 * Assigns a CTA response table to another object. The method is performing
 * a deep copy of the response table information, so that the original
 * object can be destroyed after assignment without any loss of information.
 ***************************************************************************/
GCTAResponseTable& GCTAResponseTable::operator= (const GCTAResponseTable& table)
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
 * @param[in] arg Argument.
 *
 * @exception GCTAException::bad_rsp_table_dim
 *            Response table has less than one dimension.
 *
 * Performs a linear interpolation in all parameters for one-dimensional
 * parameter tables (i.e. vectors).
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

    // Get pointer to node array (circumvent const correctness)
    GNodeArray* nodes = (GNodeArray*)&(m_axis_nodes[0]);

    // Set value for node array
    nodes->set_value(arg);

    // Get indices and weighting factors
    int    inx_left  = nodes->inx_left();
    int    inx_right = nodes->inx_right();
    double wgt_left  = nodes->wgt_left();
    double wgt_right = nodes->wgt_right();

    // Perform 1D interpolation
    for (int i = 0; i < num; ++i) {
        result[i] = wgt_left  * m_pars[i][inx_left] +
                    wgt_right * m_pars[i][inx_right];
    }
    
    // Return result vector
    return result;
}


/***********************************************************************//**
 * @brief Bilinear interpolation operator for 2D tables
 *
 * @param[in] arg1 Argument for first axis.
 * @param[in] arg2 Argument for second axis.
 *
 * @exception GCTAException::bad_rsp_table_dim
 *            Response table has less than one dimension.
 *
 * Performs a bilinear interpolation in all parameters for two-dimensional
 * parameter tables.
 ***************************************************************************/
std::vector<double> GCTAResponseTable::operator()(const double& arg1,
                                                  const double& arg2) const
{
    // Get parameter vector size
    int num = m_colname_par.size();

    // Optionally check that we have at least two dimensions
    #if defined(G_RANGE_CHECK)
    if (num < 1) {
        throw GCTAException::bad_rsp_table_dim(G_OPERATOR2, num, 2);
    }
    #endif
    
    // Initialise result vector
    std::vector<double> result(num);

    // Get pointers to node arrays (circumvent const correctness)
    GNodeArray* nodes1 = (GNodeArray*)&(m_axis_nodes[0]);
    GNodeArray* nodes2 = (GNodeArray*)&(m_axis_nodes[1]);

    // Set values for node arrays
    nodes1->set_value(arg1);
    nodes2->set_value(arg2);

    // Get array indices for bi-linear interpolation
    int size1        = axis(0);
    int offset_left  = nodes2->inx_left()  * size1;
    int offset_right = nodes2->inx_right() * size1;
    int inx1         = nodes1->inx_left()  + offset_left;
    int inx2         = nodes1->inx_left()  + offset_right;
    int inx3         = nodes1->inx_right() + offset_left;
    int inx4         = nodes1->inx_right() + offset_right;

    // Get weighting factors for bi-linear interpolation
    double wgt1 = nodes1->wgt_left()  * nodes2->wgt_left();
    double wgt2 = nodes1->wgt_left()  * nodes2->wgt_right();
    double wgt3 = nodes1->wgt_right() * nodes2->wgt_left();
    double wgt4 = nodes1->wgt_right() * nodes2->wgt_right();

    // Perform 2D interpolation
    for (int i = 0; i < num; ++i) {
        result[i] = wgt1 * m_pars[i][inx1] +
                    wgt2 * m_pars[i][inx2] +
                    wgt3 * m_pars[i][inx3] +
                    wgt4 * m_pars[i][inx4];
    }
    
    // Return result vector
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
 * @param[in] index Axis index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Axis index out of range.
 *
 * Returns the number of bins in given axis.
 ***************************************************************************/
int GCTAResponseTable::axis(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AXIS, index, size()-1);
    }
    #endif

    // Return axis length
    return (m_axis_lo[index].size());
}


/***********************************************************************//**
 * @brief Return lower bin boundary for bin in axis
 *
 * @param[in] index Axis index (starting from 0).
 * @param[in] bin Bin index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Axis or bin index out of range.
 *
 * Returns the lower boundary for a given bin in a given axis.
 ***************************************************************************/
double GCTAResponseTable::axis_lo(const int& index, const int& bin) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size() ||
        bin   < 0 || bin   >= m_axis_lo[index].size()) {
        throw GException::out_of_range(G_AXIS_LO, index, bin,
                                       size()-1, m_axis_lo[index].size()-1);
    }
    #endif

    // Return bin boundary
    return (m_axis_lo[index][bin]);
}


/***********************************************************************//**
 * @brief Return upper bin boundary for bin in axis
 *
 * @param[in] index Axis index (starting from 0).
 * @param[in] bin Bin index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Axis or bin index out of range.
 *
 * Returns the upper boundary for a given bin in a given axis.
 ***************************************************************************/
double GCTAResponseTable::axis_hi(const int& index, const int& bin) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size() ||
        bin   < 0 || bin   >= m_axis_hi[index].size()) {
        throw GException::out_of_range(G_AXIS_HI, index, bin,
                                       size()-1, m_axis_hi[index].size()-1);
    }
    #endif

    // Return bin boundary
    return (m_axis_hi[index][bin]);
}


/***********************************************************************//**
 * @brief Set nodes for a linear axis
 *
 * @param[in] index Axis index (starting from 0).
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
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AXIS_LINEAR, index, size()-1);
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
 * @param[in] index Axis index (starting from 0).
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
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AXIS_LINEAR, index, size()-1);
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
 * @param[in] index Axis index (starting from 0).
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
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AXIS_RADIANS, index, size()-1);
    }
    #endif

    // Get number of bins
    int bins = axis(index);

    // Allocate node values
    std::vector<double> axis_nodes(bins);

    // Compute nodes
    for (int i = 0; i < bins; ++i) {
        axis_nodes[i] = 0.5*(m_axis_lo[index][i] + m_axis_hi[index][i])*deg2rad;
    }

    // Set node array
    m_axis_nodes[index] = GNodeArray(axis_nodes);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read response table from FITS table HDU
 *
 * @param[in] hdu Response table HDU pointer.
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
void GCTAResponseTable::read(const GFitsTable* hdu)
{
    // Clear instance
    clear();

    // Continue only if HDU pointer is valid
    if (hdu != NULL) {

        // Read column names
        read_colnames(hdu);

        // Read axes
        read_axes(hdu);

        // Read parameter cubes
        read_pars(hdu);

    } // endif: HDU pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write response table into FITS table
 *
 * @param[in] hdu Fits table HDU.
 *
 * @todo Implement method
 ***************************************************************************/
void GCTAResponseTable::write(GFitsTable* hdu) const
{
/*
    // Allocate floating point vector columns
    GFitsTableFloatCol col_energy_lo = GFitsTableFloatCol("ENERG_LO",  1, m_energy_num);
    GFitsTableFloatCol col_energy_hi = GFitsTableFloatCol("ENERG_HI",  1, m_energy_num);
    GFitsTableFloatCol col_ctheta_lo = GFitsTableFloatCol("CTHETA_LO", 1, m_ctheta_num);
    GFitsTableFloatCol col_ctheta_hi = GFitsTableFloatCol("CTHETA_HI", 1, m_ctheta_num);

    // Set column values
    for (int i = 0; i < m_energy_num; ++i) {
        col_energy_lo(0,i) = m_energy_lo[i];
        col_energy_hi(0,i) = m_energy_hi[i];
    }
    for (int i = 0; i < m_ctheta_num; ++i) {
        col_ctheta_lo(0,i) = m_ctheta_lo[i];
        col_ctheta_hi(0,i) = m_ctheta_hi[i];
    }

    // Append columns to boundary table
    hdu->append_column(col_energy_lo);
    hdu->append_column(col_energy_hi);
    hdu->append_column(col_ctheta_lo);
    hdu->append_column(col_ctheta_hi);
*/
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print response table information
 *
 * Puts CTA response table information into a std::string object for
 * printing.
 ***************************************************************************/
std::string GCTAResponseTable::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAResponseTable ===");
    result.append("\n"+parformat("Dimension")+str(size()));
    
    // Append axes information
    for (int i = 0; i < size(); ++i) {
        result.append("\n"+parformat("Axis "+str(i)));
        result.append(str(axis(i)));
        result.append(" ["+m_colname_lo[i]);
        result.append(", ");
        result.append(m_colname_hi[i]+"]");
    }

    // Append parameter information
    for (int i = 0; i < m_colname_par.size(); ++i) {
        result.append("\n"+parformat("Parameter "+str(i)));
        result.append(m_colname_par[i]);
    }

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
    m_colname_lo.clear();
    m_colname_hi.clear();
    m_colname_par.clear();
    m_axis_lo.clear();
    m_axis_hi.clear();
    m_axis_nodes.clear();
    m_pars.clear();

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
    m_colname_lo  = table.m_colname_lo;
    m_colname_hi  = table.m_colname_hi;
    m_colname_par = table.m_colname_par;
    m_axis_lo     = table.m_axis_lo;
    m_axis_hi     = table.m_axis_hi;
    m_axis_nodes  = table.m_axis_nodes;
    m_pars        = table.m_pars;

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
 * @param[in] hdu Response table HDU pointer.
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
 * In case that the HDU pointer is not valid (NULL), this method clears the
 * column names and does nothing else.
 *
 * @todo Implement exceptions for invalid HDU format
 ***************************************************************************/
void GCTAResponseTable::read_colnames(const GFitsTable* hdu)
{
    // Clear column name arrays
    m_colname_lo.clear();
    m_colname_hi.clear();
    m_colname_par.clear();
    
    // Continue only if HDU is valid
    if (hdu != NULL) {

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
        for (int i = 0; i < hdu->ncols(); ++i) {

            // Get column name
            std::string colname = (*hdu)[i].name();

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

        } // endfor: looped over all column names

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read axes definitions from FITS HDU
 *
 * @param[in] hdu Response table HDU pointer.
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
 * In case that the HDU pointer is not valid (NULL), this method clears the
 * axes boundaries and does nothing else.
 ***************************************************************************/
void GCTAResponseTable::read_axes(const GFitsTable* hdu)
{
    // Clear axes arrays
    m_axis_lo.clear();
    m_axis_hi.clear();
    m_axis_nodes.clear();
    
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Loop over all dimensions
        for (int i = 0; i < size(); ++i) {

            // Get pointers to table columns
            const GFitsTableCol* col_lo = &(*hdu)[m_colname_lo[i]];
            const GFitsTableCol* col_hi = &(*hdu)[m_colname_hi[i]];

            // Extract number of bins. Make sure that both columns have the
            // same number of bins
            int num = col_lo->number();
            if (num != col_hi->number()) {
                std::string message = "Incompatible number of bins found"
                                      " for columns '"+m_colname_lo[i]+"'"+
                                      " ("+str(num)+") and '"+
                                      m_colname_hi[i]+"' ("+
                                      str(col_hi->number())+".";
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

            // Create node array
            m_axis_nodes.push_back(GNodeArray(axis_nodes));

        } // endfor: looped over all dimensions

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read parameter cubes
 *
 * @param[in] hdu Response table HDU pointer.
 *
 * @exception GCTAException::bad_rsp_table_format
 *            Parameter vector of bad size encountered.
 *
 * Reads the parameter cubes from the response table. 
 *
 * In case that the HDU pointer is not valid (NULL), this method clears the
 * axes boundaries and does nothing else.
 ***************************************************************************/
void GCTAResponseTable::read_pars(const GFitsTable* hdu)
{
    // Clear parameter cubes
    m_pars.clear();
    
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Compute expected cube size
        int expected_size = axis(0);
        for (int i = 1; i < size(); ++i) {
            expected_size *= axis(i);
        }
    
        // Loop over all parameter cubes
        for (int i = 0; i < m_colname_par.size(); ++i) {

            // Get pointer to table column
            const GFitsTableCol* col = &(*hdu)[m_colname_par[i]];

            // Extract number of bins. Verify that the number of bins
            // corresponds to the expectation.
            int num = col->number();
            if (num != expected_size) {
                std::string message = "Parameter vector '"+m_colname_par[i]+
                                      "' has wrong size "+str(num)+" (expected"+
                                      str(expected_size)+").";
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

        } // endfor: looped over all parameter cubes

    } // endif: HDU was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] table Response table.
 *
 * Put response table in output stream.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTAResponseTable& table)
{
     // Write response table in output stream
    os << table.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] table Response table.
 *
 * Put response table in logger.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GCTAResponseTable& table)
{
    // Write response table into logger
    log << table.print();

    // Return logger
    return log;
}
