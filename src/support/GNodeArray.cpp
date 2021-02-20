/***************************************************************************
 *                  GNodeArray.cpp - Array of nodes class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2021 by Juergen Knoedlseder                         *
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
 * @file GNodeArray.cpp
 * @brief Node array class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GException.hpp"
#include "GFilename.hpp"
#include "GNodeArray.hpp"
#include "GVector.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableCol.hpp"
#include "GFitsTableDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                           "GNodeArray::at(int&)"
#define G_INSERT                          "GNodeArray::insert(int&, double&)"
#define G_REMOVE                                   "GNodeArray::remove(int&)"
#define G_INTERPOLATE                      "GNodeArray::interpolate(double&,"\
                                                     " std::vector<double>&)"
#define G_SET_VALUE                          "GNodeArray::set_value(double&)"

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
 ***************************************************************************/
GNodeArray::GNodeArray(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS file constructor
 *
 * @param[in] filename FITS file name.
 *
 * Constructs node array from a FITS file.
 ***************************************************************************/
GNodeArray::GNodeArray(const GFilename& filename)
{
    // Initialise members
    init_members();

    // Load FITS file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Element constructor
 *
 * @param[in] num Number of elements.
 * @param[in] array Array of elements.
 *
 * Constructs a node array from an @p array of double precision values of
 * length @p num.
 ***************************************************************************/
GNodeArray::GNodeArray(const int& num, const double* array)
{
    // Initialise class members for clean destruction
    init_members();

    // Use elements to initialise nodes
    nodes(num, array);

    // Return
    return;
}


/***********************************************************************//**
 * @brief GVector constructor
 *
 * @param[in] vector Vector.
 *
 * Constructs a node array from the elements in a @p vector, represented by
 * a GVector object.
 ***************************************************************************/
GNodeArray::GNodeArray(const GVector& vector)
{
    // Initialise class members for clean destruction
    init_members();

    // Use vector to initialise nodes
    nodes(vector);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Vector constructor
 *
 * @param[in] vector Vector.
 *
 * Constructs a node array from the elements in a double precision @p vector.
 ***************************************************************************/
GNodeArray::GNodeArray(const std::vector<double>& vector)
{
    // Initialise class members for clean destruction
    init_members();

    // Use vector to initialise nodes
    nodes(vector);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] array Node array.
 ***************************************************************************/
GNodeArray::GNodeArray(const GNodeArray& array)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(array);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GNodeArray::~GNodeArray(void)
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
 * @param[in] array Node array.
 * @return Node array.
 ***************************************************************************/
GNodeArray& GNodeArray::operator=(const GNodeArray& array)
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


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear node array
 ***************************************************************************/
void GNodeArray::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone node array
 *
 * @return Pointer to deep copy of node array.
 ***************************************************************************/
GNodeArray* GNodeArray::clone(void) const
{
    return new GNodeArray(*this);
}


/***********************************************************************//**
 * @brief Node access operator
 *
 * @param[in] index Node index [0,...,size()-1].
 * @return Node value.
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 *
 * Returns a reference to the node with the specified @p index. The @p index
 * is checked on its validity.
 ***************************************************************************/
double& GNodeArray::at(const int& index)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Node index", index, size());
    }
    #endif

    // Signal that setup needs to be called
    m_need_setup = false;

    // Return node
    return m_node[index];
}


/***********************************************************************//**
 * @brief Node access operator (const version)
 *
 * @param[in] index Node index [0,...,size()-1].
 * @return Node value.
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 *
 * Returns a reference to the node with the specified @p index. The @p index
 * is checked on its validity.
 ***************************************************************************/
const double& GNodeArray::at(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Node index", index, size());
    }
    #endif

    // Signal that setup needs to be called
    m_need_setup = false;

    // Return node
    return m_node[index];
}


/***********************************************************************//**
 * @brief Set node array
 *
 * @param[in] num Number of nodes
 * @param[in] array Node values \f$x_i\f$.
 *
 * Setup node array from an array of double precision values.
 ***************************************************************************/
void GNodeArray::nodes(const int& num, const double* array)
{
    // Free any existing memory
    free_members();

    // Initialise members
    init_members();

    // Set node values
    for (int i = 0; i < num; ++i) {
        m_node.push_back(array[i]);
    }

    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append one node to array.
 *
 * @param[in] node Node.
 ***************************************************************************/
void GNodeArray::append(const double& node)
{
    // Add node
    m_node.push_back(node);

    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert one node into array.
 *
 * @param[in] index Node index [0,...,size()-1].
 * @param[in] node Node.
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 *
 * Inserts a @p node into the node array before the node with the specified
 * @p index.
 ***************************************************************************/
void GNodeArray::insert(const int& index, const double& node)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "Node index", index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "Node index", index, size());
        }
    }
    #endif

    // Inserts node
    m_node.insert(m_node.begin()+index, node);

    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove one node into array.
 *
 * @param[in] index Node index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 *
 * Remove node of specified @p index from node array.
 ***************************************************************************/
void GNodeArray::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Node index", index, size());
    }
    #endif

    // Remove node
    m_node.erase(m_node.begin() + index);

    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append node array
 *
 * @param[in] nodes Node array.
 *
 * Append a node array to the node array.
 ***************************************************************************/
void GNodeArray::extend(const GNodeArray& nodes)
{
    // Do nothing if node array is empty
    if (!nodes.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = nodes.size();

        // Reserve enough space
        reserve(size() + num);

        // Append all nodes 
        for (int i = 0; i < num; ++i) {
            m_node.push_back(nodes[i]);
        }

        // Setup node distances and linear array handling
        setup();

    } // endif: node array was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set node array from vector
 *
 * @param[in] vector Vector from which node array will be built.
 *
 * Setup node array from a vector of values.
 ***************************************************************************/
void GNodeArray::nodes(const GVector& vector)
{
    // Free any existing memory
    free_members();

    // Initialise members
    init_members();

    // Set node values
    for (int i = 0; i < vector.size(); ++i) {
        m_node.push_back(vector[i]);
    }

    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set node array from vector
 *
 * @param[in] vector Vector from which node array will be built.
 *
 * Setup node array from a vector of double precision values.
 ***************************************************************************/
void GNodeArray::nodes(const std::vector<double>& vector)
{
    // Free any existing memory
    free_members();

    // Initialise members
    init_members();

    // Set node values
    m_node = vector;

    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Interpolate value.
 *
 * @param[in] value Value \f$x\f$ at which interpolation should be done.
 * @param[in] vector Vector \f$y_i\f$ that should be interpolated.
 *
 * @exception GException::not_enough_nodes
 *            Not enough nodes for interpolation in node array.
 * @exception GException::vector_mismatch
 *            Size of node vector does not match the size of vector argument.
 *
 * This method performs a linear interpolation of values \f$y_i\f$. The
 * corresponding values \f$x_i\f$ are stored in the node array.
 ***************************************************************************/
double GNodeArray::interpolate(const double& value,
                               const std::vector<double>& vector) const
{
    // Throw exception if there are not enough nodes
    if (m_node.size() < 2) {
        throw GException::not_enough_nodes(G_INTERPOLATE, m_node.size());
    }

    // Throw exception if vectors have not the same size
    if (m_node.size() != vector.size()) {
        throw GException::vector_mismatch(G_INTERPOLATE, m_node.size(),
                                          vector.size());
    }

    // Set interpolation value
    set_value(value);

    // Interpolate
    double y = vector[inx_left()]  * wgt_left() +
               vector[inx_right()] * wgt_right();

    // Return
    return y;
}


/***********************************************************************//**
 * @brief Set indices and weighting factors for interpolation
 *
 * @param[in] value Value for which the interpolation should be done.
 *
 * @exception GException::invalid_value
 *            No nodes are available for interpolation.
 *
 * Set the indices that bound the specified value and the corresponding
 * weighting factors for linear interpolation. If the array has a linear
 * form (i.e. the nodes are equidistant), an analytic formula is used to
 * determine the boundary indices. If the nodes are not equidistant the
 * boundary indices are searched by bisection. If there is only a single
 * node, no interpolation is done and the index of this node is returned.
 *
 * Note that this method needs to be called after changing the node array
 * to setup the corrected interpolation indices and weights.
 ***************************************************************************/
void GNodeArray::set_value(const double& value) const
{
    // Get number of nodes
    int nodes = m_node.size();

    // Throw an exception if there are no nodes
    if (nodes < 1) {
        std::string msg = "Attempting to set interpolating value without "
                          "having any nodes. Interpolation can only be "
                          "done if nodes are available.";
        throw GException::invalid_value(G_SET_VALUE, msg);
    }

    // Initialize computation flag
    bool compute = true;

    // Update cache if required. If cache was updated, computation is
    // enforced. Otherwise we check if the value has changed. The value
    // check is only done when a last value has been recorded.
    if (m_need_setup) {
        setup();
    }
    else {
        if (m_has_last_value && value == m_last_value) {
            compute = false;
        }
    }

    // Continue only if computation is required
    if (compute) {

        // Handle special case of a single node
        if (nodes == 1) {
            m_inx_left       = 0;
            m_inx_right      = 0;
            m_wgt_left       = 1.0;
            m_wgt_right      = 0.0;
            m_wgt_grad_left  = 0.0;
            m_wgt_grad_right = 0.0;
        }

        // Handle all other cases
        else {

            // If array is linear then get left index from analytic formula
            if (m_is_linear) {

                // Set left index
                m_inx_left = int(m_linear_slope * value + m_linear_offset);

                // Keep index in valid range
                if (m_inx_left < 0) {
                    m_inx_left = 0;
                }
                else if (m_inx_left >= nodes-1) {
                    m_inx_left = nodes - 2;
                }

            } // endif: array is linear

            // ... otherwise search the relevant indices by bisection
            else {

                // Set left index if value is before first node
                if (value < m_node[0]) {
                    m_inx_left = 0;
                }

                // Set left index if value is after last node
                else if (value >  m_node[nodes-1]) {
                    m_inx_left = nodes - 2;
                }

                // Set left index by bisection
                else {
                    int low  = 0;
                    int high = nodes - 1;
                    while ((high - low) > 1) {
                        int mid = (low+high) / 2;
                        if (m_node[mid] > value) {
                            high = mid;
                        }
                        else {
                            low = mid;
                        }
                    }
                    m_inx_left = low;
                } // endelse: did bisection
            }

            // Set right index
            m_inx_right = m_inx_left + 1;

            // Set weighting factors and gradients
            m_wgt_grad_right = 1.0 / m_step[m_inx_left];
            m_wgt_grad_left  = -m_wgt_grad_right;
            m_wgt_right      = (value - m_node[m_inx_left]) * m_wgt_grad_right;
            m_wgt_left       = 1.0 - m_wgt_right;

        } // endelse: more than one node was present

    } // endif: computation was required

    // Store last value and signal availability
    m_last_value     = value;
    m_has_last_value = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load nodes from FITS file
 *
 * @param[in] filename FITS filename.
 *
 * Loads the node array from a FITS file.
 *
 * If no extension name is provided, the node array is loaded from the
 * "NODES" extension.
 ***************************************************************************/
void GNodeArray::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get nodes table
    const GFitsTable& table = *fits.table(filename.extname("NODES"));

    // Read nodes from table
    read(table);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save node array into FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite an existing node array extension?
 *
 * Saves node array into a FITS file. If a file with the given @p filename
 * does not yet exist it will be created, otherwise the method opens the
 * existing file. Node arrays can only be appended to an existing file if the
 * @p clobber flag is set to `true` (otherwise an exception is thrown).
 *
 * The method will append a binary FITS table containing the node array to
 * the FITS file. The extension name can be specified as part of the
 * @p filename. For example the @p filename
 *
 *      myfile.fits[NODE ARRAY]
 *
 * will save the node array in the `NODE ARRAY` extension of the
 * `myfile.fits` file. If the extension exists already in the file it will be
 * replaced, otherwise a new extension will be created. If no extension name
 * is provided, the method will use `NODES` as the default extension name
 * for the node array.
 ***************************************************************************/
void GNodeArray::save(const GFilename& filename, const bool& clobber) const
{
    // Open or create FITS file (without extension name since the requested
    // extension may not yet exist in the file)
    GFits fits(filename.url(), true);

    // Write node array to FITS object
    write(fits, filename.extname("NODES"));

    // Save to file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read nodes from FITS table
 *
 * @param[in] table FITS table.
 *
 * Reads the nodes from a FITS @p table.
 ***************************************************************************/
void GNodeArray::read(const GFitsTable& table)
{
    // Clear nodes
    clear();

    // Extract number of nodes in FITS file
    int num = table.nrows();

    // Continue only if there are nodes bins
    if (num > 0) {

        // Get the column with the name "Value"
        const GFitsTableCol* col_nodes = table["Value"];

        // Set nodes
        for (int i = 0; i < num; ++i) {
            append(col_nodes->real(i));
        }

        // Setup node distances and linear array handling
        setup();

    } // endif: there were nodes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write nodes into FITS object
 *
 * @param[in] fits FITS file.
 * @param[in] extname Nodes extension name (defaults to "NODES")
 *
 * Writes nodes into FITS object.
 ***************************************************************************/
void GNodeArray::write(GFits& fits, const std::string& extname) const
{
    // Set number of nodes
    int num = m_node.size();

    // Create node column
    GFitsTableDoubleCol col_nodes("Value", num);

    // Fill node column
    for (int i = 0; i < num; ++i) {
        col_nodes(i) = m_node[i];
    }

    // Create nodes table
    GFitsBinTable table(num);
    table.append(col_nodes);
    table.extname(extname);

    // If the FITS object contains already an extension with the same
    // name then remove now this extension
    if (fits.contains(extname)) {
        fits.remove(extname);
    }

    // Append NODES table to FITS file
    fits.append(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print nodes
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing nodes information.
 ***************************************************************************/
std::string GNodeArray::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GNodeArray ===");

        // Append array type
        result.append("\n"+gammalib::parformat("Number of nodes in array"));
        result.append(gammalib::str(size()));
        if (m_is_linear) {
            result.append("\n"+gammalib::parformat("Array type")+"linear");
            result.append("\n"+gammalib::parformat("Linear slope"));
            result.append(gammalib::str(m_linear_slope));
            result.append("\n"+gammalib::parformat("Linear offset"));
            result.append(gammalib::str(m_linear_offset));
        }
        else {
            result.append("\n"+gammalib::parformat("Array type")+"nonlinear");
        }

        // EXPLICIT: Append indices and weights
        if (chatter >= EXPLICIT) {
            result.append("\n"+gammalib::parformat("Indices and weights"));
            result.append("("+gammalib::str(m_inx_left)+",");
            result.append(gammalib::str(m_inx_right)+")=");
            result.append("("+gammalib::str(m_wgt_left)+",");
            result.append(gammalib::str(m_wgt_right)+")");
            result.append(" grad=");
            result.append("("+gammalib::str(m_wgt_grad_left)+",");
            result.append(gammalib::str(m_wgt_grad_right)+")");
        }

        // VERBOSE: Append nodes
        if (chatter >= VERBOSE) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n"+gammalib::parformat("Node "+gammalib::str(i)));
                result.append(gammalib::str(m_node[i]));
                if (i < m_step.size()) {
                    result.append(" (delta="+gammalib::str(m_step[i])+")");
                }
            }
        }

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
void GNodeArray::init_members(void)
{
    // Initialise members
    m_node.clear();
    m_step.clear();
    m_is_linear      = false;
    m_has_last_value = false;
    m_last_value     = 0.0;
    m_linear_slope   = 0.0;
    m_linear_offset  = 0.0;
    m_inx_left       = 0;
    m_inx_right      = 0;
    m_wgt_left       = 0.0;
    m_wgt_right      = 0.0;
    m_wgt_grad_left  = 0.0;
    m_wgt_grad_right = 0.0;
    m_need_setup     = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] array Node array.
 ***************************************************************************/
void GNodeArray::copy_members(const GNodeArray& array)
{
    // Copy number of bins
    m_node           = array.m_node;
    m_step           = array.m_step;
    m_is_linear      = array.m_is_linear;
    m_has_last_value = array.m_has_last_value;
    m_last_value     = array.m_last_value;
    m_linear_slope   = array.m_linear_slope;
    m_linear_offset  = array.m_linear_offset;
    m_inx_left       = array.m_inx_left;
    m_inx_right      = array.m_inx_right;
    m_wgt_left       = array.m_wgt_left;
    m_wgt_right      = array.m_wgt_right;
    m_wgt_grad_left  = array.m_wgt_grad_left;
    m_wgt_grad_right = array.m_wgt_grad_right;
    m_need_setup     = array.m_need_setup;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GNodeArray::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute distance array and linear slope/offset
 *
 * Precomputes values for fast interpolation. The precomputation requires
 * at least 2 nodes to be present in the node array. If less than two
 * nodes are present, the distance vector m_step will be empty and no
 * computation is done.
 ***************************************************************************/
void GNodeArray::setup(void) const
{
    // Reset distance vector
    m_step.clear();

    // Get number of nodes
    int nodes = m_node.size();

    // Continue only if there are at least 2 nodes
    if (nodes > 1) {

        // Setup distance array between subsequent nodes
        for (int i = 0; i < nodes-1; ++i) {
            m_step.push_back(m_node[i+1] - m_node[i]);
        }

        // Evaluate linear slope and offset
        m_linear_slope  = double(nodes-1) / (m_node[nodes-1] - m_node[0]);
        m_linear_offset = -m_linear_slope * m_node[0];

        // Check if nodes form a linear array
        m_is_linear = true;
        for (int i = 0; i < nodes-1; ++i) {
            double eps = m_linear_slope * m_node[i] + m_linear_offset - double(i);
            if (std::abs(eps) > 1.0e-6) {
                m_is_linear = false;
                break;
            }
        }

    } // endif: there were at least two nodes

    // Signal that setup has been called
    m_need_setup = false;

    // Return
    return;
}
