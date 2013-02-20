/***************************************************************************
 *                  GEbounds.cpp - Energy boundary class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GEbounds.cpp
 * @brief Energy boundary class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <iostream>
#include "GException.hpp"
#include "GTools.hpp"
#include "GEbounds.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_POP                                           "GEbounds::pop(int&)"
#define G_EMIN                                         "GEbounds::emin(int&)"
#define G_EMAX                                         "GEbounds::emax(int&)"
#define G_EMEAN                                       "GEbounds::emean(int&)"
#define G_ELOGMEAN                                 "GEbounds::elogmean(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GEbounds::GEbounds(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ebds Energy boundaries.
 ***************************************************************************/
GEbounds::GEbounds(const GEbounds& ebds)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(ebds);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Interval constructor
 *
 * @param[in] num Number of energy intervals.
 * @param[in] emin Minimum energy of first interval.
 * @param[in] emax Maximum energy of last interval.
 * @param[in] log Use logarithmic spacing? (defaults to true).
 *
 * Constructs energy boundaries by defining @p num successive energy
 * intervals between @p emin and @p emax. The @p log parameter controls
 * whether the energy spacing is logarihmic (default) or linear.
 ***************************************************************************/
GEbounds::GEbounds(const int& num, const GEnergy& emin, const GEnergy& emax,
                   const bool& log)
{
    // Initialise members
    init_members();

    // Set intervals
    if (log) {
        this->setlog(num, emin, emax);
    }
    else {
        this->setlin(num, emin, emax);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename FITS filename.
 * @param[in] extname FITS extension name (defaults to "EBOUNDS").
 ***************************************************************************/
GEbounds::GEbounds(const std::string& filename, const std::string& extname)
{
    // Initialise members
    init_members();

    // Load FITS file
    load(filename, extname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEbounds::~GEbounds(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              Operators                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] ebds Energy boundaries to be assigned.
 * @return Energy boundaries.
 ***************************************************************************/
GEbounds& GEbounds::operator=(const GEbounds& ebds)
{
    // Execute only if object is not identical
    if (this != &ebds) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(ebds);

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
 * @brief Clear energy boundaries
 ***************************************************************************/
void GEbounds::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone energy boundaries
 *
 * @return Pointer to deep copy of energy boundaries.
 ***************************************************************************/
GEbounds* GEbounds::clone(void) const
{
    return new GEbounds(*this);
}


/***********************************************************************//**
 * @brief Append energy interval
 *
 * @param[in] emin Minimum energy of interval.
 * @param[in] emax Maximum energy of interval.
 *
 * Appends an energy interval to the end of the energy boundaries
 * container.
 *
 * The method does nothing if the energy interval is not valid
 * (i.e. @p emin >= @p emax).
 ***************************************************************************/
void GEbounds::append(const GEnergy& emin, const GEnergy& emax)
{
    // Continue only if energy interval is valid
    if (emax > emin) {

        // Insert interval
        insert_eng(m_num, emin, emax);

    } // endif: Energy interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert energy interval
 *
 * @param[in] emin Minimum energy of interval.
 * @param[in] emax Maximum energy of interval.
 *
 * Inserts an energy interval into the energy boundaries after the first
 * boundary that has a minimum energy smaller than @p emin. The method
 * implicitely assumes that the intervals are ordered by increasing minimum
 * energy.
 *
 * The method does nothing if the energy interval is not valid
 * (i.e. @p emin >= @p emax).
 ***************************************************************************/
void GEbounds::insert(const GEnergy& emin, const GEnergy& emax)
{
    // Continue only if energy interval is valid
    if (emax > emin) {

        // Determine index at which interval should be inserted
        int inx = 0;
        for (; inx < m_num; ++inx) {
            if (emin < m_min[inx])
                break;
        }

        // Insert interval
        insert_eng(inx, emin, emax);

    } // endif: Energy interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Merge all overlapping or connecting successive energy intervals
 *
 * Merges all overlapping or connecting successive energy intervals. The
 * method implicitely assumes that the intervals are ordered by increasing
 * minimum energy.
 *
 * Note that the method does not actually reduce the memory size but just
 * updates the information on the number of elements in the array.
 ***************************************************************************/
void GEbounds::merge(void)
{
    // Find overlaps
    int i   = 0;
    int num = m_num;
    while (i < num-1) {

        // If current energy interval overlaps with successor then merge both
        // intervals, move all remaining intervals one position up, and
        // reduce the number of elements
        if (m_min[i+1] <= m_max[i]) {
            m_min[i] = (m_min[i] < m_min[i+1]) ? m_min[i] : m_min[i+1];
            m_max[i] = (m_max[i] > m_max[i+1]) ? m_max[i] : m_max[i+1];
            for (int k = i+2; k < num; ++k) {
                m_min[k-1] = m_min[k];
                m_max[k-1] = m_max[k];
            }
            num--;
        }

        // Otherwise increment interval index
        else {
            i++;
        }

    } // endwhile: there were still intervals to check

    // Update number of elements in object
    m_num = num;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Merge energy interval into energy boundaries
 *
 * @param[in] emin Minimum energy of interval.
 * @param[in] emax Maximum energy of interval.
 *
 * Inserts an energy interval into the energy boundaries after the first
 * boundary that has a minimum energy smaller than @p emin and then merges
 * any overlapping or connecting energy boundaries. The method implicitely
 * assumes that the intervals are ordered by increasing minimum energy.
 *
 * The method does nothing if the energy interval is not valid
 * (i.e. @p emin >= @p emax).
 ***************************************************************************/
void GEbounds::merge(const GEnergy& emin, const GEnergy& emax)
{
    // Continue only if energy interval is valid
    if (emax > emin) {

        // Determine index at which interval should be inserted
        int inx = 0;
        for (; inx < m_num; ++inx) {
            if (emin < m_min[inx])
                break;
        }

        // Insert interval
        insert_eng(inx, emin, emax);

        // Merge any overlapping energy intervals
        merge();

    } // endif: Energy interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove energy interval
 *
 * @param[in] index Energy interval index (0,...,size()-1).
 *
 * Removes energy interval at @p index from the energy boundaries container.
 * All intervals after the specified @p index are moved forward by one
 * position.
 *
 * Note that the method does not actually reduce the memory size but just
 * updates the information on the number of elements in the array.
 ***************************************************************************/
void GEbounds::pop(const int& index)
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_POP, index, 0, m_num-1);
    }
    #endif

    // Move all elements located after index forward
    for (int i = index+1; i < m_num; ++i) {
        m_min[i-1] = m_min[i];
        m_min[i-1] = m_max[i];
    }

    // Reduce number of elements by one
    m_num--;

    // Update attributes
    set_attributes();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reserve space for energy intervals
 *
 * @param[in] num Number of elements.
 *
 * This method does nothing (it is needed to satisfy the GContainer
 * interface).
 ***************************************************************************/
void GEbounds::reserve(const int& num)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append energy boundaries
 *
 * @param[in] ebds Energy boundaries.
 *
 * Append energy boundaries to the container.
 ***************************************************************************/
void GEbounds::extend(const GEbounds& ebds)
{
    // Do nothing if energy boundaries are empty
    if (!ebds.isempty()) {

        // Allocate new intervals
        int      num = m_num+ebds.size();
        GEnergy* min = new GEnergy[num];
        GEnergy* max = new GEnergy[num];

        // Initialise index
        int inx = 0;
        
        // Copy existing intervals
        for (; inx < m_num; ++inx) {
            min[inx] = m_min[inx];
            max[inx] = m_max[inx];
        }

        // Append intervals
        for (int i = 0; i < ebds.size(); ++i, ++inx) {
            min[inx] = ebds.m_min[i];
            max[inx] = ebds.m_max[i];
        }

        // Free memory
        if (m_min != NULL) delete [] m_min;
        if (m_max != NULL) delete [] m_max;

        // Set new memory
        m_min = min;
        m_max = max;

        // Set number of elements
        m_num = num;

        // Set attributes
        set_attributes();

    } // endif: energy boundaries were not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set linearly spaced energy intervals
 *
 * @param[in] num Number of energy intervals.
 * @param[in] emin Minimum energy of first interval.
 * @param[in] emax Maximum energy of last interval.
 *
 * Creates @p num linearly spaced energy boundaries running from @p emin to
 * @p emax.
 ***************************************************************************/
void GEbounds::setlin(const int& num, const GEnergy& emin, const GEnergy& emax)
{
    // Initialise members
    clear();

    // Compute bin width
    GEnergy ebin = (emax - emin)/double(num); 

    // Append boundaries
    GEnergy min = emin;
    GEnergy max = emin + ebin;
    for (int i = 0; i < num; ++i) {
        append(min, max);
        min += ebin;
        max += ebin;
    }

    // Set attributes
    set_attributes();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set logarithmically spaced energy intervals
 *
 * @param[in] num Number of energy intervals.
 * @param[in] emin Minimum energy of first interval.
 * @param[in] emax Maximum energy of last interval.
 *
 * Creates @p num logarithmically spaced energy boundaries running from
 * @p emin to @p emax.
 ***************************************************************************/
void GEbounds::setlog(const int& num, const GEnergy& emin, const GEnergy& emax)
{
    // Initialise members
    clear();

    // Compute bin width
    double elogmin = std::log10(emin.MeV());
    double elogmax = std::log10(emax.MeV());
    double elogbin = (elogmax - elogmin)/double(num);

    // Append boundaries
    GEnergy min;
    GEnergy max;
    for (int i = 0; i < num; ++i) {
        min.MeV(std::pow(10.0, double(i)*elogbin   + elogmin));
        max.MeV(std::pow(10.0, double(i+1)*elogbin + elogmin));
        append(min, max);
    }

    // Set attributes
    set_attributes();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy boundaries from FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] extname FITS extension name (defaults to "EBOUNDS").
 *
 * Loads the energy boundaries from FITS file.
 ***************************************************************************/
void GEbounds::load(const std::string& filename, const std::string& extname)
{
    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(filename);

    // Get energy boundary HDU
    GFitsTable* hdu = file.table(extname);

    // Read energy boundaries from HDU
    read(hdu);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy boundaries into FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite any existing file?
 * @param[in] extname Energy boundary extension name (defaults to "EBOUNDS").
 *
 * Saves the energy boundaries into extension @p extname of a FITS file.
 ***************************************************************************/
void GEbounds::save(const std::string& filename, bool clobber,
                    const std::string& extname) const
{
    // Allocate FITS file
    GFits file;

    // Write energy boundaries to FITS file
    write(&file, extname);

    // Save to file
    file.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy boundaries from FITS table
 *
 * @param[in] hdu Pointer to FITS table.
 *
 * Reads the energy boundaries from a FITS table. It is assumed that the
 * energies are stored in units of keV.
 *
 * @todo Needs to interpret energy units. We could also add an optional
 *       string parameter that allows external specification about how
 *       the energies should be interpreted.
 ***************************************************************************/
void GEbounds::read(GFitsTable* hdu)
{
    // Free members
    free_members();

    // Initialise attributes
    init_members();

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Extract energy boundary information from FITS file
        m_num = hdu->integer("NAXIS2");
        if (m_num > 0) {

            // Allocate memory
            m_min = new GEnergy[m_num];
            m_max = new GEnergy[m_num];

            // Copy information
            for (int i = 0; i < m_num; ++i) {
                m_min[i].keV((*hdu)["E_MIN"].real(i));
                m_max[i].keV((*hdu)["E_MAX"].real(i));
            }

        } // endif: there were channels to read

        // Set attributes
        set_attributes();

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energy boundaries into FITS object
 *
 * @param[in] file Pointer to FITS file.
 * @param[in] extname Energy boundary extension name (defaults to "EBOUNDS")
 *
 * Writes the energy boundaries into a FITS object. Energies are written in
 * units of keV.
 *
 * @todo Write header keywords.
 * @todo Should we allow for the possibility to specify the energy units?
 ***************************************************************************/
void GEbounds::write(GFits* file, const std::string& extname) const
{
    // Create energy boundary columns
    GFitsTableDoubleCol cemin = GFitsTableDoubleCol("E_MIN", m_num);
    GFitsTableDoubleCol cemax = GFitsTableDoubleCol("E_MAX", m_num);

    // Fill energy boundary columns
    for (int i = 0; i < m_num; ++i) {
        cemin(i) = m_min[i].keV();
        cemax(i) = m_max[i].keV();
    }

    // Set energy units to keV
    cemin.unit("keV");
    cemax.unit("keV");

    // Create binary table
    GFitsBinTable* table = new GFitsBinTable(m_num);
    table->append_column(cemin);
    table->append_column(cemax);
    table->extname(extname);

    // Write to FITS file
    file->append(*table);

    // Free table
    delete table;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns energy bin index for a given energy
 *
 * @param[in] eng Energy.
 * @return Bin index.
 *
 * Returns the energy boundary bin index for a given energy. If the energy
 * falls outside all boundaries, -1 is returned.
 ***************************************************************************/
int GEbounds::index(const GEnergy& eng) const
{
    // Initialise index with 'not found'
    int index = -1;

    // Search all energy boundaries for containment
    for (int i = 0; i < m_num; ++i) {
        if (eng >= m_min[i] && eng <= m_max[i]) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Returns minimum energy for a given energy interval
 *
 * @param[in] index Energy interval index (0,...,size()-1).
 * @return Minimum energy of interval.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GEnergy GEbounds::emin(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_EMIN, index, 0, m_num-1);
    }
    #endif

    // Return
    return (m_min[index]);
}


/***********************************************************************//**
 * @brief Returns maximum energy for a given energy interval
 *
 * @param[in] index Energy interval index (0,...,size()-1).
 * @return Maximum energy of interval.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GEnergy GEbounds::emax(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_EMAX, index, 0, m_num-1);
    }
    #endif

    // Return
    return (m_max[index]);
}


/***********************************************************************//**
 * @brief Returns mean energy for a given energy interval
 *
 * @param[in] index Energy interval index (0,...,size()-1).
 * @return Mean energy of interval.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 *
 * Computes the mean energy
 * \f$0.5 * (E_{\rm min} + E_{\rm max})\f$
 * for the energy interval @p index.
 ***************************************************************************/
GEnergy GEbounds::emean(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_EMEAN, index, 0, m_num-1);
    }
    #endif

    // Compute mean energy
    GEnergy emean = 0.5 * (m_min[index] + m_max[index]);

    // Return
    return emean;
}


/***********************************************************************//**
 * @brief Returns logarithmic mean energy for a given energy interval
 *
 * @param[in] index Energy interval index (0,...,size()-1).
 * @return Logarithmic mean energy of interval.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 *
 * Computes the logarithmic mean energy
 * \f$10^{0.5 * (\log E_{\rm min} + \log E_{\rm max})}\f$
 * for the energy interval @p index.
 ***************************************************************************/
GEnergy GEbounds::elogmean(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_ELOGMEAN, index, 0, m_num-1);
    }
    #endif

    // Compute logarithmic mean energy
    GEnergy elogmean;
    double  elogmin = std::log10(m_min[index].MeV());
    double  elogmax = std::log10(m_max[index].MeV());
    elogmean.MeV(std::pow(10.0, 0.5 * (elogmin + elogmax)));

    // Return
    return elogmean;
}


/***********************************************************************//**
 * @brief Checks whether energy boundaries contain energy
 *
 * @param[in] eng Energy to be checked.
 * @return True if energy falls in at least one interval, false otherwise.
 *
 * Checks if the energy @p eng falls in at least one of the energy intervals.
 * The method exits when the first matching interval has been found.
 ***************************************************************************/
bool GEbounds::contains(const GEnergy& eng) const
{
    // Initialise test
    bool found = false;

    // Test all energy boundaries
    for (int i = 0; i < m_num; ++i) {
        if (eng >= m_min[i] && eng <= m_max[i]) {
            found = true;
            break;
        }
    }

    // Return result
    return found;
}


/***********************************************************************//**
 * @brief Print energy boundaries
 *
 * @return String containing energy boundary information.
 ***************************************************************************/
std::string GEbounds::print(void) const
{
    // Initialise result string
    std::string result;

    // Append Header
    result.append("=== GEbounds ===\n");

    // Append information
    result.append(parformat("Number of intervals")+str(size()));
    result.append("\n");
    result.append(parformat("Energy range"));
    result.append(emin().print());
    result.append(" - ");
    result.append(emax().print());
    
    // If there are multiple energy bins then append them
    if (size() > 1) {
        for (int i = 0; i < size(); ++i) {
            result.append("\n");
            result.append(parformat("Energy interval "+str(i)));
            result.append(emin(i).print());
            result.append(" - ");
            result.append(emax(i).print());
        }
    }

    // Return
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
void GEbounds::init_members(void)
{
    // Initialise members
    m_num = 0;
    m_emin.clear();
    m_emax.clear();
    m_min = NULL;
    m_max = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ebds Energy boundaries.
 ***************************************************************************/
void GEbounds::copy_members(const GEbounds& ebds)
{
    // Copy attributes
    m_num  = ebds.m_num;
    m_emin = ebds.m_emin;
    m_emax = ebds.m_emax;

    // Copy arrays
    if (m_num > 0) {
        m_min = new GEnergy[m_num];
        m_max = new GEnergy[m_num];
        for (int i = 0; i < m_num; ++i) {
            m_min[i] = ebds.m_min[i];
            m_max[i] = ebds.m_max[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEbounds::free_members(void)
{
    // Free memory
    if (m_min != NULL) delete [] m_min;
    if (m_max != NULL) delete [] m_max;

    // Signal free pointers
    m_min = NULL;
    m_max = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set class attributes
 *
 * Determines the minimum and maximum energy from all intervals. If no
 * interval is present the minimum and maximum energies are cleared.
 ***************************************************************************/
void GEbounds::set_attributes(void)
{
    // If there are intervals then determine the minimum and maximum
    // energy from these intervals ...
    if (m_num > 0) {
        m_emin = m_min[0];
        m_emax = m_max[0];
        for (int i = 1; i < m_num; ++i) {
            if (m_min[i] < m_emin) m_emin = m_min[i];
            if (m_max[i] > m_emax) m_emax = m_max[i];
        }
    }

    // ... otherwise clear the minimum and maximum energy
    else {
        m_emin.clear();
        m_emax.clear();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert energy interval
 *
 * @param[in] index Index after with interval is inserted.
 * @param[in] emin Minimum energy of interval.
 * @param[in] emax Maximum energy of interval.
 *
 * Inserts an energy interval after the specified @p index in the energy
 * boundaries. The method does not reorder the intervals by energy, instead
 * the client needs to determine the approriate @p index.
 *
 * Invalid parameters do not produce any exception, but are handled
 * transparently. If the interval is invalid (i.e. @p emin >= @p emax) then
 * nothing is done. If the @p index is out of the valid range, the index
 * will be adjusted to either the first or the last element.
 ***************************************************************************/
void GEbounds::insert_eng(const int&     index,
                          const GEnergy& emin,
                          const GEnergy& emax)
{
    // Continue only if energy interval is valid
    if (emax > emin) {

        // Set index
        int inx = index;

        // If inx is out of range then adjust it
        if (inx < 0)     inx = 0;
        if (inx > m_num) inx = m_num;

        // Allocate new intervals
        int      num = m_num+1;
        GEnergy* min = new GEnergy[num];
        GEnergy* max = new GEnergy[num];

        // Copy intervals before index to be inserted
        for (int i = 0; i < inx; ++i) {
            min[i] = m_min[i];
            max[i] = m_max[i];
        }

        // Insert interval
        min[inx] = emin;
        max[inx] = emax;

        // Copy intervals after index to be inserted
        for (int i = inx+1; i < num; ++i) {
            min[i] = m_min[i-1];
            max[i] = m_max[i-1];
        }

        // Free memory
        if (m_min != NULL) delete [] m_min;
        if (m_max != NULL) delete [] m_max;

        // Set new memory
        m_min = min;
        m_max = max;

        // Set number of elements
        m_num = num;

        // Set attributes
        set_attributes();

    } // endif: Energy interval was valid

    // Return
    return;
}
