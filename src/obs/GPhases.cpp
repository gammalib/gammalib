/***************************************************************************
 *                  GPhases.cpp - Phase intervals class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2022 by Juergen Knoedlseder                         *
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
 * @file GPhases.cpp
 * @brief Phase intervals class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GPhases.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_REMOVE                                      "GPhases::remove(int&)"
#define G_PMIN                                          "GPhases::pmin(int&)"
#define G_PMAX                                          "GPhases::pmax(int&)"
#define G_INSERT_INTERVAL  "GPhases::insert_interval(int&, double&, double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs empty phase intervals.
 ***************************************************************************/
GPhases::GPhases(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] phases Phase intervals.
 *
 * Constructs phase intervals by copying other phase intervals.
 ***************************************************************************/
GPhases::GPhases(const GPhases& phases)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(phases);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Single phase interval constructor
 *
 * @param[in] pmin Lower boundary of phase interval.
 * @param[in] pmax Upper boundary of phase interval.
 *
 * Constructs Good Time Intervals from a single phase interval, given by
 * [p@ pmin, @p pmax].
 ***************************************************************************/
GPhases::GPhases(const double& pmin, const double& pmax)
{
    // Initialise members
    init_members();

    // Append phase interval
    append(pmin, pmax);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPhases::~GPhases(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] phases Phase intervals.
 * @return Phase intervals.
 ***************************************************************************/
GPhases& GPhases::operator=(const GPhases& phases)
{
    // Execute only if object is not identical
    if (this != &phases) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(phases);

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
 * @brief Clear phase intervals
 ***************************************************************************/
void GPhases::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone phase intervals
 *
 * @return Pointer to deep copy of phase intervals.
 ***************************************************************************/
GPhases* GPhases::clone(void) const
{
    return new GPhases(*this);
}


/***********************************************************************//**
 * @brief Check whether phase is contained in phases
 *
 * @param[in] phase Phase.
 * @returns Phase is contained in phases.
 *
 * Checks whether a phase is contained in one of the phase intervals. The
 * lower phase bound is included while the upper phase bound is excluded
 * by the check.
 ***************************************************************************/
bool GPhases::contains(const double& phase) const
{
    // Initialise containment flag
    bool contained = false;

    // Check if phase is contained in one of the phase intervals
    for (int i = 0; i < size(); ++i) {
        if (phase >= m_pmin[i] && phase < m_pmax[i]) {
            contained = true;
            break;
        }
    }

    // Return containment flag
    return contained;
}


/***********************************************************************//**
 * @brief Append phase interval
 *
 * @param[in] pmin Lower boundary of phase interval.
 * @param[in] pmax Upper boundary of phase interval.
 *
 * Appends a phase interval at the end of the container.
 ***************************************************************************/
void GPhases::append(const double& pmin, const double& pmax)
{
    // Insert phase interval at the end of the list
    insert_interval(size(), pmin, pmax);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove phase interval
 *
 * @param[in] index Phase interval index (0,...,size()-1).
 *
 * Removes phase interval at @p index from the container. All intervals
 * after the specified @p index are moved forward by one position.
 ***************************************************************************/
void GPhases::remove(const int& index)
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Phase interval", index,
                                       size());
    }
    #endif

    // Remove lower and upper boundaries
    m_pmin.erase(m_pmin.begin() + index);
    m_pmax.erase(m_pmax.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reserve space for phase intervals
 *
 * @param[in] num Number of elements.
 ***************************************************************************/
void GPhases::reserve(const int& num)
{
    // Reserves memory
    m_pmin.reserve(num);
    m_pmax.reserve(num);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append phase intervals
 *
 * @param[in] phases Phase intervals.
 *
 * Append phase intervals to the container.
 ***************************************************************************/
void GPhases::extend(const GPhases& phases)
{
    // Continue only if phase intervals are not empty
    if (!phases.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = phases.size();

        // Resize the boundary vectors
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_pmin.push_back(phases.pmin(i));
            m_pmax.push_back(phases.pmax(i));
        }

    } // endif: Phase intervals were not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns lower boundary for a given phase interval
 *
 * @param[in] index Phase interval index (0,...,size()-1).
 * @return Lower boundary of phase interval.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
double GPhases::pmin(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_PMIN, "Phase interval", index,
                                       size());
    }
    #endif

    // Return
    return (m_pmin[index]);
}


/***********************************************************************//**
 * @brief Returns upper boundary for a given phase interval
 *
 * @param[in] index Phase interval index (0,...,size()-1).
 * @return Upper boundary of phase interval.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
double GPhases::pmax(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_PMAX, "Phase interval", index,
                                       size());
    }
    #endif

    // Return
    return (m_pmax[index]);
}


/***********************************************************************//**
 * @brief Returns total length of phase intervals
 *
 * @return Total length of phase intervals.
 ***************************************************************************/
double GPhases::length(void) const
{
    // Initialise length
    double length = 0.0;

    // Loop over all intervals and add their length
    for (int i = 0; i < size(); ++i) {
        length += (m_pmax[i] - m_pmin[i]);
    }

    // Return
    return length;
}


/***********************************************************************//**
 * @brief Print phase intervals
 *
 * @param[in] chatter Chattiness.
 * @return String containing phase interval information.
 ***************************************************************************/
std::string GPhases::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GPhases ===");

        // Append phase interval information
        result.append("\n"+gammalib::parformat("Number of phase intervals"));
        result.append(gammalib::str(size()));

        // Append phases
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+gammalib::parformat("Phase interval "+
                          gammalib::str(i+1)));
            result.append(gammalib::str(m_pmin[i]));
            result.append(" - ");
            result.append(gammalib::str(m_pmax[i]));
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
void GPhases::init_members(void)
{
    // Initialise members
    m_pmin.clear();
    m_pmax.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] phases Phase intervals.
 ***************************************************************************/
void GPhases::copy_members(const GPhases& phases)
{
    // Copy attributes
    m_pmin = phases.m_pmin;
    m_pmax = phases.m_pmax;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPhases::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert phase interval
 *
 * @param[in] index Index after which interval is inserted.
 * @param[in] pmin Lower boundary of interval.
 * @param[in] pmax Upper boundary of interval.
 *
 * @exception GException::invalid_argument
 *            Invalid phase interval boundaries specified
 *
 * Inserts a phase interval before the interval with the specified @p index.
 * If no interval with the specified index exists then append the interval
 * at the end of the existing phase intervals.
 ***************************************************************************/
void GPhases::insert_interval(const int&    index,
                              const double& pmin,
                              const double& pmax)
{
    // Throw an exception if phase minimum is outside the range [0,1]
    if (pmin < 0.0 || pmin > 1.0) {
        std::string msg = "Phase minimum "+gammalib::str(pmin)+
                          " outside the valid range [0,1]. Please "
                          "specify phase interval boundaries comprised "
                          "within 0 and 1.";
        throw GException::invalid_argument(G_INSERT_INTERVAL, msg);
    }

    // Throw an exception if phase minimum is outside the range [0,1]
    if (pmax < 0.0 || pmax > 1.0) {
        std::string msg = "Phase maximum "+gammalib::str(pmax)+
                          " outside the valid range [0,1]. Please "
                          "specify phase interval boundaries comprised "
                          "within 0 and 1.";
        throw GException::invalid_argument(G_INSERT_INTERVAL, msg);
    }

    // If phase minimum is smaller than phase maximum then append or insert
    // a single interval ...
    if (pmin < pmax) {
        if (index >= size()) {
            m_pmin.push_back(pmin);
            m_pmax.push_back(pmax);
        }
        else {
            m_pmin.insert(m_pmin.begin()+index, pmin);
            m_pmax.insert(m_pmax.begin()+index, pmax);
        }
    }

    // ... otherwise if the phase minimum is larger than the maximum then
    // consider this a wrap around interval and append or insert one interval
    // from phase minimum to 1 and another interval from 0 to phase maximum
    else if (pmin > pmax) {
        if (pmin < 1.0) {
            if (index >= size()) {
                m_pmin.push_back(pmin);
                m_pmax.push_back(1.0);
            }
            else {
                m_pmin.insert(m_pmin.begin()+index, pmin);
                m_pmax.insert(m_pmax.begin()+index, 1.0);
            }
        }
        if (pmax > 0.0) {
            if (index >= size()) {
                m_pmin.push_back(0.0);
                m_pmax.push_back(pmax);
            }
            else {
                m_pmin.insert(m_pmin.begin()+index, 0.0);
                m_pmax.insert(m_pmax.begin()+index, pmax);
            }
        }
    }

    // Return
    return;
}
