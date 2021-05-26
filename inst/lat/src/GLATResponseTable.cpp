/***************************************************************************
 *          GLATResponseTable.cpp - Fermi/LAT Response table class         *
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
 * @file GLATResponseTable.cpp
 * @brief Fermi/LAT response table class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATResponseTable.hpp"
#include "GTools.hpp"
#include "GException.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                         "GLATResponseTable::read(GFitsTable&)"
#define G_INDEX                        "GLATResponseTable::index(int&, int&)"
#define G_ENERGY                            "GLATResponseTable::energy(int&)"
#define G_ENERGY_LO                      "GLATResponseTable::energy_lo(int&)"
#define G_ENERGY_HI                      "GLATResponseTable::energy_hi(int&)"
#define G_COSTHETA_LO                  "GLATResponseTable::costheta_lo(int&)"
#define G_COSTHETA_HI                  "GLATResponseTable::costheta_hi(int&)"

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
GLATResponseTable::GLATResponseTable(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param table Response table.
 ***************************************************************************/
GLATResponseTable::GLATResponseTable(const GLATResponseTable& table)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATResponseTable::~GLATResponseTable(void)
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
 * @param table Response table.
 * @return Response table.
 ***************************************************************************/
GLATResponseTable& GLATResponseTable::operator=(const GLATResponseTable& table)
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


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATResponseTable::clear(void)
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
 ***************************************************************************/
GLATResponseTable* GLATResponseTable::clone(void) const
{
    return new GLATResponseTable(*this);
}


/***********************************************************************//**
 * @brief Read response table from FITS table HDU
 *
 * @param[in] hdu Response table HDU.
 *
 * The response table definition is assumed to be stored in 4 vector columns
 * with names
 * ENERG_LO (energy bins lower boundary)
 * ENERG_HI (energy bins upper boundary)
 * CTHETA_LO (cos theta bins lower boundary)
 * CTHETA_HI (cos theta bins upper boundary)
 ***************************************************************************/
void GLATResponseTable::read(const GFitsTable& hdu)
{
    // Clear instance
    clear();

    // Get pointers to table columns
    const GFitsTableCol* energy_lo = hdu["ENERG_LO"];
    const GFitsTableCol* energy_hi = hdu["ENERG_HI"];
    const GFitsTableCol* ctheta_lo = hdu["CTHETA_LO"];
    const GFitsTableCol* ctheta_hi = hdu["CTHETA_HI"];

    // Extract number of bins
    m_energy_num = energy_lo->number();
    m_ctheta_num = ctheta_lo->number();

    // Allocate memory
    if (m_energy_num > 0) {
        m_energy_lo = new double[m_energy_num];
        m_energy_hi = new double[m_energy_num];
    }
    if (m_ctheta_num > 0) {
        m_ctheta_lo = new double[m_ctheta_num];
        m_ctheta_hi = new double[m_ctheta_num];
    }

    // Transfer data
    for (int i = 0; i < m_energy_num; ++i) {
        m_energy_lo[i] = energy_lo->real(0,i);
        m_energy_hi[i] = energy_hi->real(0,i);
    }
    for (int i = 0; i < m_ctheta_num; ++i) {
        m_ctheta_lo[i] = ctheta_lo->real(0,i);
        m_ctheta_hi[i] = ctheta_hi->real(0,i);
    }

    // Set energy nodes (log10 of energy)
    if (m_energy_num > 0) {
        double *logE = new double[m_energy_num];
        for (int i = 0; i < m_energy_num; ++i) {
            logE[i] = 0.5 * (log10(m_energy_lo[i]) + log10(m_energy_hi[i]));
//ST method logE[i] = log10(sqrt(m_energy_lo[i]*m_energy_hi[i]));
            m_energy.push_back(std::pow(10.0, logE[i]));
        }
        m_logE.nodes(m_energy_num, logE);
        delete [] logE;
    }

    // Set cos theta nodes
    if (m_ctheta_num > 0) {
        double *ctheta = new double[m_ctheta_num];
        for (int i = 0; i < m_ctheta_num; ++i) {
            ctheta[i] = 0.5 * (m_ctheta_lo[i] + m_ctheta_hi[i]);
        }
        m_ctheta.nodes(m_ctheta_num, ctheta);
        delete [] ctheta;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write response table into FITS table
 *
 * @param[in] hdu Fits table HDU.
 ***************************************************************************/
void GLATResponseTable::write(GFitsTable& hdu) const
{
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
    hdu.append(col_energy_lo);
    hdu.append(col_energy_hi);
    hdu.append(col_ctheta_lo);
    hdu.append(col_ctheta_hi);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return table index
 *
 * @param[in] ie Energy bin [0,...,m_energy_num[
 * @param[in] ic cos theta bin [0,...,m_ctheta_num[
 *
 * @exception GException::out_of_range
 *            Energy or cos theta bin index out of range
 *
 * Computes table index from energy and cos theta bin indices.
 ***************************************************************************/
int GLATResponseTable::index(const int& ie, const int& ic) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (ie < 0 || ie >= m_energy_num) {
        throw GException::out_of_range(G_INDEX, "Energy bin index",
                                       ie, m_energy_num);
    }
    if (ic < 0 || ic >= m_ctheta_num) {
        throw GException::out_of_range(G_INDEX, "cos theta bin index",
                                       ic, m_ctheta_num);
    }
    #endif

    // Compute index
    int index = ic * m_energy_num + ie;

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Return mean energy of bin (units: MeV)
 *
 * @param[in] ie Index of energy bin [0,...,m_energy_num[
 *
 * @exception GException::out_of_range
 *            Energy bin index out of range
 *
 * The mean energy is actually computed from the logarithmic average of 
 * lower and upper boundaries:
 * \f$E=10^{0.5 \times (\log{E_{\rm min}} + \log{E_{\rm max}})}\f$
 * Note that this energy is stored in the m_energy array, hence to convert
 * to MeV we simply have to the it to the power of 10.
 *
 * @todo Store also linear energies to avoid conversion.
 ***************************************************************************/
double GLATResponseTable::energy(const int& ie) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (ie < 0 || ie >= m_energy_num) {
        throw GException::out_of_range(G_ENERGY, "Energy bin index",
                                       ie, m_energy_num);
    }
    #endif

    // Determine mean energy of bin
    //double mean   = 0.5 * (log10(m_energy_lo[ie]) + log10(m_energy_hi[ie]));
    //double mean   = m_logE[ie];
    //double energy = pow(10.0, mean);

    // Return mean energy
    return (m_energy[ie]);
}


/***********************************************************************//**
 * @brief Set indices and weighting for bi-linear interpolation of 2D array
 *
 * @param[in] logE Base 10 logarithm of the energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 *
 * Bi-linear interpolation is performed in log10 of energy and in cos theta.
 ***************************************************************************/
void GLATResponseTable::set(const double& logE, const double& ctheta)
{
    // Flag no change of values
    bool change = false;

    // Set energy interpolation
    if (logE != m_last_energy) {
        m_logE.set_value(logE);
        m_last_energy = logE;
        change        = true;
    }

    // Set cos(theta) interpolation
    if (ctheta != m_last_ctheta) {
        m_ctheta.set_value(ctheta);
        m_last_ctheta = ctheta;
        change        = true;
    }

    // If change occured then update interpolation indices and weighting
    // factors
    if (change) {

        // Set array indices for bi-linear interpolation
        int inx_ctheta_left  = m_ctheta.inx_left()  * m_energy_num;
        int inx_ctheta_right = m_ctheta.inx_right() * m_energy_num;
        m_inx1 = m_logE.inx_left()  + inx_ctheta_left;
        m_inx2 = m_logE.inx_left()  + inx_ctheta_right;
        m_inx3 = m_logE.inx_right() + inx_ctheta_left;
        m_inx4 = m_logE.inx_right() + inx_ctheta_right;

        // Set weighting factors for bi-linear interpolation
        m_wgt1 = m_logE.wgt_left()  * m_ctheta.wgt_left();
        m_wgt2 = m_logE.wgt_left()  * m_ctheta.wgt_right();
        m_wgt3 = m_logE.wgt_right() * m_ctheta.wgt_left();
        m_wgt4 = m_logE.wgt_right() * m_ctheta.wgt_right();

    } // endif: logE or ctheta changed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Perform bi-linear interpolation of 2D array
 *
 * @param[in] logE Base 10 logarithm of the energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 * @param[in] array Array to be interpolated.
 *
 * Bi-linear interpolation is performed in log10 of energy and in cos theta.
 * The array is stored in a std::vector object with the logE axis varying
 * more rapidely.
 ***************************************************************************/
double GLATResponseTable::interpolate(const double&              logE,
                                      const double&              ctheta,
                                      const std::vector<double>& array)
{
    // Set interpolation indices and weights
    set(logE, ctheta);

    // Perform bi-linear interpolation
    double value = m_wgt1 * array[m_inx1] +
                   m_wgt2 * array[m_inx2] +
                   m_wgt3 * array[m_inx3] +
                   m_wgt4 * array[m_inx4];

    // Return bi-linear interpolated value
    return value;
}


/***********************************************************************//**
 * @brief Perform bi-linear interpolation of 3D array
 *
 * @param[in] logE Base 10 logarithm of the energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 * @param[in] array Array to be interpolated.
 * @param[in] offset Offset if 3D array in 1st dimension.
 * @param[in] size Size of 3D array in 1st dimension.
 *
 * Bi-linear interpolation is performed in log10 of energy and in cos theta.
 * The array is stored in a std::vector object.
 ***************************************************************************/
double GLATResponseTable::interpolate(const double&              logE, 
                                      const double&              ctheta, 
                                      const std::vector<double>& array,
                                      const int&                 offset,
                                      const int&                 size)
{
    // Set interpolation indices and weights
    set(logE, ctheta);

    // Perform bi-linear interpolation
    double value = m_wgt1 * array[m_inx1*size+offset] +
                   m_wgt2 * array[m_inx2*size+offset] +
                   m_wgt3 * array[m_inx3*size+offset] +
                   m_wgt4 * array[m_inx4*size+offset];

    // Return bi-linear interpolated value
    return value;
}


/***********************************************************************//**
 * @brief Return lower bin energy (units: MeV)
 *
 * @param[in] inx Index of energy bin [0,...,m_energy_num[
 *
 * @exception GException::out_of_range
 *            Energy bin index out of range
 ***************************************************************************/
double GLATResponseTable::energy_lo(const int& inx) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= m_energy_num) {
        throw GException::out_of_range(G_ENERGY_LO, "Energy bin index",
                                       inx, m_energy_num);
    }
    #endif

    // Return mean energy
    return m_energy_lo[inx];
}


/***********************************************************************//**
 * @brief Return upper bin energy (units: MeV)
 *
 * @param[in] inx Index of energy bin [0,...,m_energy_num[
 *
 * @exception GException::out_of_range
 *            Energy bin index out of range
 ***************************************************************************/
double GLATResponseTable::energy_hi(const int& inx) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= m_energy_num) {
        throw GException::out_of_range(G_ENERGY_HI, "Energy bin index",
                                       inx, m_energy_num);
    }
    #endif

    // Return mean energy
    return m_energy_hi[inx];
}


/***********************************************************************//**
 * @brief Return lower bin cos theta
 *[
 * @param[in] inx Index of cos theta bin [0,...,m_ctheta_num[
 *
 * @exception GException::out_of_range
 *            Cos theta bin index out of range
 ***************************************************************************/
double GLATResponseTable::costheta_lo(const int& inx) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= m_ctheta_num) {
        throw GException::out_of_range(G_COSTHETA_LO, "cos thera bin index",
                                       inx, m_ctheta_num);
    }
    #endif

    // Return mean energy
    return m_ctheta_lo[inx];
}


/***********************************************************************//**
 * @brief Return upper bin cos theta
 *
 * @param[in] inx Index of cos theta bin (starting from 0)
 *
 * @exception GException::out_of_range
 *            Cos theta bin index out of range
 ***************************************************************************/
double GLATResponseTable::costheta_hi(const int& inx) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= m_ctheta_num) {
        throw GException::out_of_range(G_COSTHETA_HI, "cos thera bin index",
                                       inx, m_ctheta_num);
    }
    #endif

    // Return mean energy
    return m_ctheta_hi[inx];
}


/***********************************************************************//**
 * @brief Return indices of 4 corners used for interpolation
 ***************************************************************************/
std::vector<int> GLATResponseTable::indices(void) const
{
    // Define vector
    std::vector<int> incides;

    // Push indices on vector
    incides.push_back(m_inx1);
    incides.push_back(m_inx2);
    incides.push_back(m_inx3);
    incides.push_back(m_inx4);

    // Return vector
    return incides;
}


/***********************************************************************//**
 * @brief Return energies of 4 corners used for interpolation
 ***************************************************************************/
std::vector<double> GLATResponseTable::energies(void) const
{
    // Define vector
    std::vector<double> energies;

    // Push energies on vector
    if (m_energy_num > 0) {
        energies.push_back(m_energy[m_logE.inx_left()]);
        energies.push_back(m_energy[m_logE.inx_left()]);
        energies.push_back(m_energy[m_logE.inx_right()]);
        energies.push_back(m_energy[m_logE.inx_right()]);
    }

    // Return vector
    return energies;
}


/***********************************************************************//**
 * @brief Return weights of 4 corners used for interpolation
 ***************************************************************************/
std::vector<double> GLATResponseTable::weights(void) const
{
    // Define vector
    std::vector<double> weights;

    // Push indices on vector
    weights.push_back(m_wgt1);
    weights.push_back(m_wgt2);
    weights.push_back(m_wgt3);
    weights.push_back(m_wgt4);

    // Return vector
    return weights;
}


/***********************************************************************//**
 * @brief Print response table information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing response table information.
 ***************************************************************************/
std::string GLATResponseTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATResponseTable ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(nenergies()));
        result.append("\n"+gammalib::parformat("Number of cos theta bins") +
                      gammalib::str(ncostheta()));

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
 ***************************************************************************/
void GLATResponseTable::init_members(void)
{
    // Initialise members
    m_energy_num  = 0;
    m_ctheta_num  = 0;
    m_energy.clear();
    m_logE.clear();
    m_ctheta.clear();
    m_energy_lo   = NULL;
    m_energy_hi   = NULL;
    m_ctheta_lo   = NULL;
    m_ctheta_hi   = NULL;
    m_last_energy = -1.0;
    m_last_ctheta = -1.0;
    m_inx1        = 0;
    m_inx2        = 0;
    m_inx3        = 0;
    m_inx4        = 0;
    m_wgt1        = 0.0;
    m_wgt2        = 0.0;
    m_wgt3        = 0.0;
    m_wgt4        = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] table Response table.
 ***************************************************************************/
void GLATResponseTable::copy_members(const GLATResponseTable& table)
{
    // Copy number of bins
    m_energy_num  = table.m_energy_num;
    m_ctheta_num  = table.m_ctheta_num;
    m_energy      = table.m_energy;
    m_logE        = table.m_logE;
    m_ctheta      = table.m_ctheta;
    m_last_energy = table.m_last_energy;
    m_last_ctheta = table.m_last_ctheta;
    m_inx1        = table.m_inx1;
    m_inx2        = table.m_inx2;
    m_inx3        = table.m_inx3;
    m_inx4        = table.m_inx4;
    m_wgt1        = table.m_wgt1;
    m_wgt2        = table.m_wgt2;
    m_wgt3        = table.m_wgt3;
    m_wgt4        = table.m_wgt4;

    // Copy energy bins
    if (m_energy_num > 0) {
        m_energy_lo   = new double[m_energy_num];
        m_energy_hi   = new double[m_energy_num];
        for (int i=0; i < m_energy_num; ++i) {
            m_energy_lo[i] = table.m_energy_lo[i];
            m_energy_hi[i] = table.m_energy_hi[i];
        }
    }

    // Copy cos theta bins
    if (m_ctheta_num > 0) {
        m_ctheta_lo = new double[m_ctheta_num];
        m_ctheta_hi = new double[m_ctheta_num];
        for (int i=0; i < m_ctheta_num; ++i) {
            m_ctheta_lo[i] = table.m_ctheta_lo[i];
            m_ctheta_hi[i] = table.m_ctheta_hi[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATResponseTable::free_members(void)
{
    // Free memory
    if (m_energy_lo != NULL) delete [] m_energy_lo;
    if (m_energy_hi != NULL) delete [] m_energy_hi;
    if (m_ctheta_lo != NULL) delete [] m_ctheta_lo;
    if (m_ctheta_hi != NULL) delete [] m_ctheta_hi;

    // Signal that memory is free
    m_energy_lo = NULL;
    m_energy_hi = NULL;
    m_ctheta_lo = NULL;
    m_ctheta_hi = NULL;

    // Return
    return;
}
