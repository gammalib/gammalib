/***************************************************************************
 *                  GEbounds.cpp - Energy boundary class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2020 by Juergen Knoedlseder                         *
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
#include "GException.hpp"
#include "GTools.hpp"
#include "GFilename.hpp"
#include "GEbounds.hpp"
#include "GEnergies.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GXmlElement.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SET        "GEbounds::set(int&, GEnergy&, GEnergy&, std::string&, "\
                                                                   "double&)"
#define G_READ_XML                             "GEbounds::read(GXmlElement&)"
#define G_WRITE_XML                           "GEbounds::write(GXmlElement&)"
#define G_REMOVE                                     "GEbounds::remove(int&)"
#define G_EMIN_SET                           "GEbounds::emin(int&, GEnergy&)"
#define G_EMAX_SET                           "GEbounds::emax(int&, GEnergy&)"
#define G_EMIN_GET                                     "GEbounds::emin(int&)"
#define G_EMAX_GET                                     "GEbounds::emax(int&)"
#define G_EMEAN                                       "GEbounds::emean(int&)"
#define G_ELOGMEAN                                 "GEbounds::elogmean(int&)"
#define G_EWIDTH                                     "GEbounds::ewidth(int&)"
#define G_INSERT_ENG         "GEbounds::insert_eng(int&, GEnergy&, GEnergy&)"
#define G_SET_LIN               "GEbounds::set_lin(int&, GEnergy&, GEnergy&)"
#define G_SET_LOG               "GEbounds::set_log(int&, GEnergy&, GEnergy&)"
#define G_SET_POW      "GEbounds::set_pow(int&, GEnergy&, GEnergy&, double&)"

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
 * @brief FITS file constructor
 *
 * @param[in] filename FITS file name.
 *
 * Constructs energy boundaries from a FITS file.
 ***************************************************************************/
GEbounds::GEbounds(const GFilename& filename)
{
    // Initialise members
    init_members();

    // Load FITS file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML element constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs energy boundaries from an XML element.
 ***************************************************************************/
GEbounds::GEbounds(const GXmlElement& xml)
{
    // Initialise members
    init_members();

    // Read energy boundaries from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy container constructor
 *
 * @param[in] energies Energy container.
 *
 * Constructs energy boundaries from an energy container.
 ***************************************************************************/
GEbounds::GEbounds(const GEnergies& energies)
{
    // Initialise members
    init_members();

    // Set energy boundaries from energy container
    set(energies);

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
 * @brief Single energy band constructor
 *
 * @param[in] emin Minimum energy of the interval.
 * @param[in] emax Maximum energy of the interval.
 *
 * Constructs energy boundaries for one (emin, emax) energy band.
 ***************************************************************************/
GEbounds::GEbounds(const GEnergy& emin, const GEnergy& emax)
{
    // Initialise members
    init_members();

    // Append energies
    append(emin, emax);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Interval constructor
 *
 * @param[in] num Number of energy intervals.
 * @param[in] emin Minimum energy of first interval.
 * @param[in] emax Maximum energy of last interval.
 * @param[in] method Energy spacing method (one of "LIN", "LOG" or "POW").
 * @param[in] gamma Power law index for @p POW method.
 *
 * Constructs energy boundaries by defining @p num successive energy
 * intervals between @p emin and @p emax. The @p method parameter controls
 * the energy spacing of the energy boundaries. See the set() method for
 * more information.
 ***************************************************************************/
GEbounds::GEbounds(const int&         num,
                   const GEnergy&     emin,
                   const GEnergy&     emax,
                   const std::string& method,
                   const double&      gamma)
{
    // Initialise members
    init_members();

    // Set intervals
    set(num, emin, emax, method, gamma);

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
 ***************************************************************************/
void GEbounds::append(const GEnergy& emin, const GEnergy& emax)
{
    // Append interval at the end
    insert_eng(m_num, emin, emax);

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
 ***************************************************************************/
void GEbounds::insert(const GEnergy& emin, const GEnergy& emax)
{
    // Determine index at which interval should be inserted
    int inx = 0;
    for (; inx < m_num; ++inx) {
        if (emin < m_min[inx])
            break;
    }

    // Insert interval
    insert_eng(inx, emin, emax);

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
 ***************************************************************************/
void GEbounds::merge(const GEnergy& emin, const GEnergy& emax)
{
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
void GEbounds::remove(const int& index)
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_REMOVE, index, 0, m_num-1);
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
    if (!ebds.is_empty()) {

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
 * @brief Set energy boundaries from energy container
 *
 * @param[in] energies Energy container.
 *
 * Sets the energy boundaries from an energy container. Each two subsequent
 * energies in the energy container will form an energy boundary. This
 * means that n energies will lead to n-1 energy boundaries with the
 * following mapping:
 *
 *      [energies[0], energies[1]]
 *      [energies[1], energies[2]]
 *      ...
 *      [energies[n-2], energies[n-1]]
 *
 * If there is only one energy in the container the following empty energy
 * boundary will be appended:
 *
 *      [energies[0], energies[0]]
 ***************************************************************************/
void GEbounds::set(const GEnergies& energies)
{
    // Initialise members
    clear();

    // Get number of energies in container
    int num = energies.size();

    // If there is only one energy in the container then append an empty
    // energy boundary
    if (num == 1) {
        append(energies[0], energies[0]);
    }

    // ... otherwise if there is more than one energy in the container
    // then append subsequent energies as boundaries
    else if (num > 1) {
        for (int i = 0; i < num-1; ++i) {
            append(energies[i], energies[i+1]);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set linearly spaced energy intervals
 *
 * @param[in] num Number of energy intervals.
 * @param[in] emin Minimum energy of first interval.
 * @param[in] emax Maximum energy of last interval.
 * @param[in] method Energy spacing method (one of "LIN", "LOG" or "POW").
 * @param[in] gamma Power law index for @p POW method.
 *
 * @exception GException::invalid_argument
 *            Minimum energy is larger than the maximum energy.
 *
 * Sets energy boundaries by defining @p num successive energy intervals
 * between @p emin and @p emax. The @p method parameter controls the energy
 * spacing of the energy boundaries. See the set_lin(), set_log() and
 * set_pow() methods for more information.
 *
 * If @p num is not a positive integer, no energy boundaries will be set.
 ***************************************************************************/
void GEbounds::set(const int&         num,
                   const GEnergy&     emin,
                   const GEnergy&     emax,
                   const std::string& method,
                   const double&      gamma)
{
    // Convert method to upper-case string
    std::string umethod = gammalib::toupper(method);

    // Dispatch to corresponding set methos
    if (umethod == "LIN") {
        set_lin(num, emin, emax);
    }
    else if (umethod == "LOG") {
        set_log(num, emin, emax);
    }
    else if (umethod == "POW") {
        set_pow(num, emin, emax, gamma);
    }
    else {
        std::string msg = "Invalid energy spacing method \""+umethod+"\" "
                          "specified. Please provide one of \"LIN\", \"LOG\""
                          ", or \"POW\".";
        throw GException::invalid_argument(G_SET, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy boundaries from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the energy boundaries from a FITS file.
 *
 * If no extension name is provided, the energy boundaries are loaded from
 * the `EBOUNDS` extension.
 ***************************************************************************/
void GEbounds::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get energy boundary table
    const GFitsTable& table = *fits.table(filename.extname(gammalib::extname_ebounds));

    // Read energy boundaries from table
    read(table);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy boundaries into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite an existing energy boundaries extension?
 * @param[in] unit Energy unit
 *
 * Saves energy boundaries into a FITS file. If a file with the given
 * @p filename does not yet exist it will be created, otherwise the method
 * opens the existing file. Energy boundaries can only be appended to an
 * existing file if the @p clobber flag is set to "true" (otherwise an
 * exception is thrown).
 *
 * The method will append a binary FITS table containing the energy
 * boundaries to the FITS file. The extension name can be specified as part
 * of the @p filename. For example the @p filename
 *
 *      myfile.fits[ENERGY BOUNDARIES]
 *
 * will save the energy boundaries in the `ENERGY BOUNDARIES` extension of
 * the "myfile.fits" file. If the extension exists already in the file it
 * will be replaced, otherwise a new extension will be created. If no
 * extension name is provided, the method will use `EBOUNDS` as the default
 * extension name for energy boundaries.
 ***************************************************************************/
void GEbounds::save(const GFilename&   filename,
                    const bool&        clobber,
                    const std::string& unit) const
{
    // Open or create FITS file (without extension name since the requested
    // extension may not yet exist in the file)
    GFits fits(filename.url(), true);

    // Write energy boundaries to FITS file
    write(fits, filename.extname(gammalib::extname_ebounds), unit);

    // Save to file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy boundaries from FITS table
 *
 * @param[in] table FITS table.
 *
 * Reads the energy boundaries from a FITS table. The method interprets the
 * energy units provide in the FITS header. If no energy units are found it
 * is assumed that the energies are stored in units of keV.
 ***************************************************************************/
void GEbounds::read(const GFitsTable& table)
{
    // Free members
    free_members();

    // Initialise attributes
    init_members();

    // Extract energy boundary information from FITS table
    m_num = table.integer("NAXIS2");
    if (m_num > 0) {

        // Allocate memory
        m_min = new GEnergy[m_num];
        m_max = new GEnergy[m_num];

        // Get units
        std::string emin_unit = table["E_MIN"]->unit();
        std::string emax_unit = table["E_MAX"]->unit();
        if (emin_unit.empty()) {
            emin_unit = "keV";
        }
        if (emax_unit.empty()) {
            emax_unit = "keV";
        }

        // Copy information
        for (int i = 0; i < m_num; ++i) {
            m_min[i](table["E_MIN"]->real(i), emin_unit);
            m_max[i](table["E_MAX"]->real(i), emax_unit);
        }

        // Set attributes
        set_attributes();

    } // endif: there were channels to read

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energy boundaries into FITS object
 *
 * @param[in] fits FITS file.
 * @param[in] extname Energy boundary extension name.
 * @param[in] unit Energy units.
 *
 * Writes the energy boundaries into a FITS object. The @p unit parameter
 * specifies in which unit the energies are written. By default, the energy
 * units are keV.
 *
 * @todo Write header keywords.
 ***************************************************************************/
void GEbounds::write(GFits&             fits,
                     const std::string& extname,
                     const std::string& unit) const
{
    // Create energy boundary columns
    GFitsTableDoubleCol cemin("E_MIN", m_num);
    GFitsTableDoubleCol cemax("E_MAX", m_num);

    // Fill energy boundary columns
    for (int i = 0; i < m_num; ++i) {
        cemin(i) = m_min[i](unit);
        cemax(i) = m_max[i](unit);
    }

    // Set energy units
    cemin.unit(unit);
    cemax.unit(unit);

    // Create binary table
    GFitsBinTable table(m_num);
    table.append(cemin);
    table.append(cemax);
    table.extname(extname);

    // If the FITS object contains already an extension with the same
    // name then remove now this extension
    if (fits.contains(extname)) {
        fits.remove(extname);
    }

    // Append energy boundary table to FITS file
    fits.append(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy boundaries from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Read energy boundaries from an XML element. The format of the energy
 * boundaries is
 *
 *     <parameter name="EnergyBoundaries" emin="0.1" emax="10.0"/>
 *
 * The units of the @a emin and @a emax parameters are MeV.
 ***************************************************************************/
void GEbounds::read(const GXmlElement& xml)
{
    // Clear energy boundaries
    clear();

    // Get energy boundaries parameter
    const GXmlElement* par = gammalib::xml_get_par(G_READ_XML, xml,
                                                   "EnergyBoundaries");

    // Extract position attributes
    if (par->has_attribute("emin") && par->has_attribute("emax")) {
        double emin = gammalib::todouble(par->attribute("emin"));
        double emax = gammalib::todouble(par->attribute("emax"));
        append(GEnergy(emin, "MeV"), GEnergy(emax, "MeV"));
    }
    else {
        std::string msg = "Attributes \"emin\" and/or \"emax\" not found"
                          " in XML parameter \"EnergyBoundaries\"."
                          " Please verify the XML format.";
        throw GException::invalid_value(G_READ_XML, msg);
    }

    // Set attribues
    set_attributes();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energy boundaries into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes energy boundaries into an XML element. The format of the energy
 * boundaries is
 *
 *     <parameter name="EnergyBoundaries" emin="0.1" emax="10.0"/>
 *
 * The units of the @a emin and @a emax parameters are MeV.
 *
 * This method does nothing if the energy boundaries are empty.
 ***************************************************************************/
void GEbounds::write(GXmlElement& xml) const
{
    // Continue only if there are energy boundaries
    if (!is_empty()) {

        // Get parameter
        GXmlElement* par = gammalib::xml_need_par(G_WRITE_XML, xml,
                                                  "EnergyBoundaries");

        // Write attributes
        par->attribute("emin", gammalib::str(emin().MeV()));
        par->attribute("emax", gammalib::str(emax().MeV()));

    } // endif: energy boundaries were not empty
        
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns energy bin index for a given energy
 *
 * @param[in] eng Energy.
 * @return Bin index.
 *
 * Returns the energy boundary bin index for a given energy. By convention,
 * the limits for an energy bin are defined as
 *
 *      min <= energy < max
 *
 * i.e. and energy equals to max falls above the largest energy.
 *
 * If the energy falls outside all boundaries, -1 is returned.
 ***************************************************************************/
int GEbounds::index(const GEnergy& eng) const
{
    // Initialise index with 'not found'
    int index = -1;

    // Search all energy boundaries for containment
    for (int i = 0; i < m_num; ++i) {
        if (eng >= m_min[i] && eng < m_max[i]) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Set minimum energy for a given energy interval
 *
 * @param[in] index Energy interval index (0,...,size()-1).
 * @param[in] energy Minimum energy of interval.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 *
 * Sets the minimum energy for the energy interval @p index.
 ***************************************************************************/
void GEbounds::emin(const int& index, const GEnergy& energy)
{
    #if defined(G_RANGE_CHECK)
    // Throw an exception if index is outside valid range
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_EMIN_SET,
                                       "Minimum energy of interval",
                                       index, m_num);
    }
    #endif

    // Set minimum energy
    m_min[index] = energy;

    // Set attributes
    set_attributes();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set maximum energy for a given energy interval
 *
 * @param[in] index Energy interval index (0,...,size()-1).
 * @param[in] energy Maximum energy of interval.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 *
 * Sets the maximum energy for the energy interval @p index.
 ***************************************************************************/
void GEbounds::emax(const int& index, const GEnergy& energy)
{
    #if defined(G_RANGE_CHECK)
    // Throw an exception if index is outside valid range
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_EMAX_SET,
                                       "Maximum energy of interval",
                                       index, m_num);
    }
    #endif

    // Set maximum energy
    m_max[index] = energy;

    // Set attributes
    set_attributes();

    // Return
    return;
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
    // Throw an exception if index is outside valid range
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_EMIN_GET,
                                       "Minimum energy of interval",
                                       index, m_num);
    }
    #endif

    // Return minimum energy
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
    // Throw an exception if index is outside valid range
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_EMAX_GET,
                                       "Maximum energy of interval",
                                       index, m_num);
    }
    #endif

    // Return maximum energy
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
    elogmean.MeV(std::sqrt(m_min[index].MeV() * m_max[index].MeV()));

    // Return
    return elogmean;
}


/***********************************************************************//**
 * @brief Returns energy interval width
 *
 * @param[in] index Energy interval index (0,...,size()-1).
 * @return Energy interval width.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GEnergy GEbounds::ewidth(const int& index) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_num) {
        throw GException::out_of_range(G_EWIDTH, index, 0, m_num-1);
    }
    #endif

    // Return
    return (m_max[index]-m_min[index]);
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
 * @brief Checks whether energy boundaries contain and energy bin
 *
 * @param[in] emin Minimum energy of bin to be checked.
 * @param[in] emax Maximum energy of bin to be checked.
 * @return True if energy bin [emin, emax] is fully contained inside energy
 * the boundaries, false otherwiese
 *
 * Checks if the energy interval [@p emin, @p emax ] falls is fully contained
 * by energy boundaries
 *
 * @todo This method is so far only correct for contiguous energy boundaries.
 ***************************************************************************/
bool GEbounds::contains(const GEnergy& emin, const GEnergy& emax) const
{
    // Initialise result
    bool contained = false;

    // Check if emin, emax is fully contained within this ebounds
    if (emin >= m_emin && emax <= m_emax) {
        contained = true;
    }

    // Return result
    return contained;
}


/***********************************************************************//**
 * @brief Print energy boundaries
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy boundary information.
 ***************************************************************************/
std::string GEbounds::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append Header
        result.append("=== GEbounds ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of intervals"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(emin().print());
        result.append(" - ");
        result.append(emax().print());
    
        // If there are multiple energy bins then append them
        if (chatter >= EXPLICIT) {
            if (size() > 1) {
                for (int i = 0; i < size(); ++i) {
                    result.append("\n");
                    result.append(gammalib::parformat("Energy interval "+
                                  gammalib::str(i)));
                    result.append(emin(i).print());
                    result.append(" - ");
                    result.append(emax(i).print());
                }
            }
        } // endif: chatter was less than explicit

    } // endif: chatter was not silent

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
 * @exception GException::invalid_argument
 *            Minimum energy larger than maximum energy
 *
 * Inserts an energy interval after the specified @p index in the energy
 * boundaries. The method does not reorder the intervals by energy, instead
 * the client needs to determine the approriate @p index.
 *
 * Invalid parameters do not produce any exception, but are handled
 * transparently. If the interval is invalid (i.e. @p emin > @p emax) an
 * exception is thrown. If the @p index is out of the valid range, the
 * index will be adjusted to either the first or the last element.
 ***************************************************************************/
void GEbounds::insert_eng(const int&     index,
                          const GEnergy& emin,
                          const GEnergy& emax)
{
    // Throw an exception if energy interval is invalid
    if (emin > emax) {
        std::string msg = "Invalid energy interval specified. Minimum"
                          " energy "+emin.print(NORMAL)+" can not be"
                          " larger than maximum energy "+
                          emax.print(NORMAL)+".";
        throw GException::invalid_argument(G_INSERT_ENG, msg);
    }

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
 * @exception GException::invalid_argument
 *            Minimum energy is larger than the maximum energy.
 *
 * Creates \f$N\f$ linearly spaced energy boundaries running from
 * \f$E_{\rm min}\f$ to \f$E_{\rm max}\f$ defined by
 *
 * \f[
 *    E_{\rm min}^{(i)} = E_{\rm min} + \frac{E_{\rm max} - E_{\rm min}}{N}
 *                        \times i
 * \f]
 *
 * and
 *
 * \f[
 *    E_{\rm max}^{(i)} = E_{\rm min} + \frac{E_{\rm max} - E_{\rm min}}{N}
 *                        \times (i + 1)
 * \f]
 *
 * If @p num is not positive the method will do nothing.
 ***************************************************************************/
void GEbounds::set_lin(const int& num, const GEnergy& emin, const GEnergy& emax)
{
    // Continue only if the requested number of energy boundaries is positive
    if (num > 0) {

        // Throw an exception if the maximum energy is smaller than the
        // minimum energy
        if (emin > emax) {
            std::string msg = "Minimum energy "+emin.print()+" is larger than "
                              "maximum energy "+emax.print()+". Please provide "
                              "a minimum energy that is not larger than the "
                              "maximum energy.";
            throw GException::invalid_argument(G_SET_LIN, msg);
        }

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

    } // endif: number of energy boundaries was positive

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
 * @exception GException::invalid_argument
 *            Minimum or maximum energy are not positive or minimum energy
 *            is larger than the maximum energy.
 *
 * Creates \f$N\f$ logarithmically spaced energy boundaries running from
 * \f$E_{\rm min}\f$ to \f$E_{\rm max}\f$ defined by
 *
 * \f[
 *    E_{\rm min}^{(i)} = 10^{\log_{10} E_{\rm min} +
 *                            \frac{\log_{10} E_{\rm max} -
 *                                  \log_{10} E_{\rm min}}{N}
 *                            \times i}
 * \f]
 *
 * and
 *
 * \f[
 *    E_{\rm max}^{(i)} = 10^{\log_{10} E_{\rm min} +
 *                            \frac{\log_{10} E_{\rm max} -
 *                                  \log_{10} E_{\rm min}}{N}
 *                            \times (i + 1)}
 * \f]
 *
 * If @p num is not positive the method will do nothing.
 ***************************************************************************/
void GEbounds::set_log(const int& num, const GEnergy& emin, const GEnergy& emax)
{
    // Initialise members
    clear();

    // Continue only if the requested number of energy boundaries is positive
    if (num > 0) {

        // Throw an exception if the minimum or maximum energy is not positive
        if (emin.MeV() <= 0.0) {
            std::string msg = "Non-positive minimum energy "+emin.print()+
                              " specified. Please provide a positive minimum "
                              "energy value.";
            throw GException::invalid_argument(G_SET_LOG, msg);
        }
        if (emax.MeV() <= 0.0) {
            std::string msg = "Non-positive maximum energy "+emax.print()+
                              " specified. Please provide a positive minimum "
                              "energy value.";
            throw GException::invalid_argument(G_SET_LOG, msg);
        }

        // Throw an exception if the maximum energy is smaller than the
        // minimum energy
        if (emin > emax) {
            std::string msg = "Minimum energy "+emin.print()+" is larger than "
                              "maximum energy "+emax.print()+". Please provide "
                              "a minimum energy that is not larger than the "
                              "maximum energy.";
            throw GException::invalid_argument(G_SET_LOG, msg);
        }

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

    } // endif: number of energy boundaries was positive

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set power-law spaced energy intervals
 *
 * @param[in] num Number of energy intervals.
 * @param[in] emin Minimum energy of first interval.
 * @param[in] emax Maximum energy of last interval.
 * @param[in] gamma Power law index.
 *
 * @exception GException::invalid_argument
 *            Minimum or maximum energy are not positive or minimum energy
 *            is larger than the maximum energy.
 *
 * Creates \f$N\f$ power-law spaced energy boundaries running from
 * \f$E_{\rm min}\f$ to \f$E_{\rm max}\f$.
 *
 * If @p num is not positive the method will do nothing.
 ***************************************************************************/
void GEbounds::set_pow(const int&     num,
                       const GEnergy& emin,
                       const GEnergy& emax,
                       const double&  gamma)
{
    // Initialise members
    clear();

    // Continue only if the requested number of energy boundaries is positive
    if (num > 0) {

        // Throw an exception if the minimum or maximum energy is not positive
        if (emin.MeV() <= 0.0) {
            std::string msg = "Non-positive minimum energy "+emin.print()+
                              " specified. Please provide a positive minimum "
                              "energy value.";
            throw GException::invalid_argument(G_SET_POW, msg);
        }
        if (emax.MeV() <= 0.0) {
            std::string msg = "Non-positive maximum energy "+emax.print()+
                              " specified. Please provide a positive minimum "
                              "energy value.";
            throw GException::invalid_argument(G_SET_POW, msg);
        }

        // Throw an exception if the maximum energy is smaller than the
        // minimum energy
        if (emin > emax) {
            std::string msg = "Minimum energy "+emin.print()+" is larger than "
                              "maximum energy "+emax.print()+". Please provide "
                              "a minimum energy that is not larger than the "
                              "maximum energy.";
            throw GException::invalid_argument(G_SET_POW, msg);
        }

        // Precomputation
        double a = 1.0 - gamma;
        double c = (a == 0.0) ? 1.0 / (std::log(emax.MeV())   - std::log(emin.MeV()))
                              : a   / (std::pow(emax.MeV(),a) - std::pow(emin.MeV(),a));
        //double b = double(c*(num-1.0));
        double b = double(c*num);

        // Append boundaries
        double e = emin.MeV();
        for (int i = 0; i < num; ++i) {

            // Compute upper bin boundary
            double log_e_next = (a == 0.0)
                                ? 1.0/b + std::log(e)
                                : std::log(a/b + std::pow(e,a)) / a;
            double e_next = std::exp(log_e_next);

            // Append energy bin
            append(GEnergy(e, "MeV"), GEnergy(e_next, "MeV"));

            // Store upper boundary as next lower boundary
            e = e_next;

        }

        // Set attributes
        set_attributes();

    } // endif: number of energy boundaries was positive

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Energy boundaries equality operator friend
 *
 * @param[in] a First energy boundaries.
 * @param[in] b Second energy boundaries.
 * @return True if both energy boundaries are identical.
 ***************************************************************************/
bool operator==(const GEbounds& a, const GEbounds& b)
{
    // Initialise identify flag
    bool identity = true;

    // Check that both energy boundaries have the same size
    if (a.size() != b.size()) {
        identity = false;
    }

    // Check all energy boundaries
    else {
        for (int i = 0; i < a.size(); ++i) {
            if ((a.emin(i) != b.emin(i)) || (a.emax(i) != b.emax(i))) {
                identity = false;
                break;
            }
        }
    }

    // Return identity flag
    return identity;
}
