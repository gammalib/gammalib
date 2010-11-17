/***************************************************************************
 *                 GEbounds.cpp  -  Energy boundary class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEbounds.cpp
 * @brief Energy boundary class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <math.h>
#include "GException.hpp"
#include "GEbounds.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_EMIN                                          "GEbounds::emin(int)"
#define G_EMAX                                          "GEbounds::emax(int)"
#define G_EMEAN                                        "GEbounds::emean(int)"
#define G_ELOGMEAN                                  "GEbounds::elogmean(int)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                          Constructors/destructors                       =
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
 * @param[in] ebds Energy boundaries from which the instance should be built.
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
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] ebds Energy boundaries to be assigned.
 ***************************************************************************/
GEbounds& GEbounds::operator= (const GEbounds& ebds)
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
 * @brief Append energy boundaries
 *
 * @param[in] emin Minimum energy of interval to be appended.
 * @param[in] emax Maximum energy of interval to be appended.
 *
 * This method appends a new energy interval at the end of the list of
 * existing intervals. No checking for energy ordering or overlap is done.
 * If the energy interval is not valid (emax <= emin), nothing is done.
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
 * @brief Insert energy boundaries at the approprita location
 *
 * @param[in] emin Minimum energy of interval to be appended.
 * @param[in] emax Maximum energy of interval to be appended.
 *
 * This method inserts a new energy interval at the appropriate location
 * into the object. Energy interval ordering is done based on the minimum
 * energy. No checking is done for energy interval overlap.
 * If the energy interval is not valid (emax <= emin), nothing is done.
 ***************************************************************************/
void GEbounds::insert(const GEnergy& emin, const GEnergy& emax)
{
    // Continue only if energy interval is valid
    if (emax > emin) {

        // Determine index at which interval should be inserted
        int inx = 0;
        for (int i = 0; i < m_num; ++i) {
            if (emin < m_min[i])
                break;
        }
        
        // Insert interval
        insert_eng(inx, emin, emax);

    } // endif: Energy interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set linearly spaced energy bins
 *
 * @param[in] emin Minimum energy of boundaries.
 * @param[in] emax Maximum energy of boundaries.
 * @param[in] num Number of energy bins.
 ***************************************************************************/
void GEbounds::setlin(const GEnergy& emin, const GEnergy& emax, const int& num)
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set logarithmically spaced energy bins
 *
 * @param[in] emin Minimum energy of boundaries.
 * @param[in] emax Maximum energy of boundaries.
 * @param[in] num Number of energy bins.
 ***************************************************************************/
void GEbounds::setlog(const GEnergy& emin, const GEnergy& emax, const int& num)
{
    // Initialise members
    clear();
    
    // Compute bin width
    double elogmin = log10(emin.MeV());
    double elogmax = log10(emax.MeV());
    double elogbin = (elogmax - elogmin)/double(num);

    // Append boundaries
    GEnergy min;
    GEnergy max;
    for (int i = 0; i < num; ++i) {
        min.MeV(pow(10.0, double(i)*elogbin   + elogmin));
        max.MeV(pow(10.0, double(i+1)*elogbin + elogmin));
        append(min, max);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy boundaries from file.
 *
 * @param[in] filename FITS filename from which GEbounds is to be loaded.
 * @param[in] extname FITS extension name of the energy boundaries.
 *
 * This method loads the energy boundary definitions from a FITS file.
 ***************************************************************************/
void GEbounds::load(const std::string& filename, const std::string& extname)
{
    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(filename);

    // Get energy boundary HDU
    GFitsBinTable* hdu = (GFitsBinTable*)file.hdu(extname);

    // Read energy boundaries from HDU
    read(hdu);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy boundaries to FITS file.
 *
 * @param[in] filename Name of file into which energy boundaries are to be saved.
 * @param[in] clobber Overwrite any existing file.
 * @param[in] extname Energy boundary extension name (default is "EBOUNDS")
 ***************************************************************************/
void GEbounds::save(const std::string& filename, bool clobber,
                    const std::string& extname)
{
    // Allocate FITS file
    GFits file;

    // Write GTI to FITS file
    write(&file, extname);

    // Save to file
    file.saveto(filename, clobber);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy boundaries from HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which GEbounds are loaded.
 *
 * This method loads the energy boundary definitions from a FITS HDU. It
 * assumes that the energy is stored in units of keV.
 *
 * @todo Needs to interpret energy units. We could also add an optional
 *       string parameter that allows external specification about how
 *       the energies should be interpreted.
 ***************************************************************************/
void GEbounds::read(GFitsBinTable* hdu)
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
                m_min[i].keV(hdu->column("E_MIN")->real(i));
                m_max[i].keV(hdu->column("E_MAX")->real(i));
            }

        } // endif: there were channels to read

        // Set attributes
        set_attributes();

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energy boundaries into FITS file.
 *
 * @param[in] file Pointer to FITS file.
 * @param[in] extname Energy boundary extension name (default is "EBOUNDS")
 *
 * This methods write energy boundaries in units of keV.
 *
 * @todo Write header keywords.
 ***************************************************************************/
void GEbounds::write(GFits* file, const std::string& extname)
{
    // Create energy boundary columns
    GFitsTableDoubleCol cemin = GFitsTableDoubleCol("E_MIN", m_num);
    GFitsTableDoubleCol cemax = GFitsTableDoubleCol("E_MAX", m_num);

    // Fill energy boundary columns
    for (int i = 0; i < m_num; ++i) {
        cemin(i) = m_min[i].keV();
        cemax(i) = m_max[i].keV();
    }

    // Create binary table
    GFitsBinTable* table = new GFitsBinTable(m_num);
    table->append_column(cemin);
    table->append_column(cemax);
    table->extname(extname);

    // Write to FITS file
    file->append(table);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns energy bin index for a given energy
 *
 * @param[in] eng Energy for which index is to be returned.
 *
 * Returns the energy boundary bin index for a given energy. If the energy
 * falls outside any of the boundaries, -1 is returned.
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
 * @brief Returns minimum energy for a given bin
 *
 * @param[in] inx Energy bin.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GEnergy GEbounds::emin(int inx) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (inx < 0 || inx >= m_num)
        throw GException::out_of_range(G_EMIN, inx, 0, m_num-1);
    #endif

    // Return
    return (m_min[inx]);
}


/***********************************************************************//**
 * @brief Returns maximum energy for a given bin
 *
 * @param[in] inx Energy bin.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GEnergy GEbounds::emax(int inx) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (inx < 0 || inx >= m_num)
        throw GException::out_of_range(G_EMAX, inx, 0, m_num-1);
    #endif

    // Return
    return (m_max[inx]);
}


/***********************************************************************//**
 * @brief Returns mean energy for a given bin
 *
 * @param[in] inx Energy bin.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GEnergy GEbounds::emean(int inx) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (inx < 0 || inx >= m_num)
        throw GException::out_of_range(G_EMEAN, inx, 0, m_num-1);
    #endif

    // Compute mean energy
    GEnergy emean = 0.5 * (m_min[inx] + m_max[inx]);

    // Return
    return emean;
}


/***********************************************************************//**
 * @brief Returns logarithmic mean energy for a given bin
 *
 * @param[in] inx Energy bin.
 *
 * @exception GException::out_of_range
 *            Specified index is out of range.
 ***************************************************************************/
GEnergy GEbounds::elogmean(int inx) const
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (inx < 0 || inx >= m_num)
        throw GException::out_of_range(G_ELOGMEAN, inx, 0, m_num-1);
    #endif

    // Compute logarithmic mean energy
    GEnergy elogmean;
    double  elogmin  = log10(m_min[inx].MeV());
    double  elogmax  = log10(m_max[inx].MeV());
    elogmean.MeV(pow(10.0, 0.5 * (elogmin + elogmax)));

    // Return
    return elogmean;
}


/***********************************************************************//**
 * @brief Tests whether energy is within energy boundaries intervals
 *
 * @param[in] eng Energy to be tested.
 ***************************************************************************/
bool GEbounds::isin(const GEnergy& eng) const
{
    // Initialise test
    bool found = false;

    // Test all GTIs
    for (int i = 0; i < m_num; ++i) {
        if (eng >= m_min[i] && eng <= m_max[i]) {
            found = true;
            break;
        }
    }

    // Return result
    return found;
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
 * @param[in] ebds GEbounds members which should be copied.
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
 ***************************************************************************/
void GEbounds::set_attributes(void)
{
    // Continue only if there are intervals
    if (m_num > 0) {
    
        // Set attributes
        m_emin = m_min[0];
        m_emax = m_max[m_num-1];

    } // endif: there were intervals
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone GEbounds
 *
 * Cloning provides a copy of the GEbounds. Cloning is used to allocate
 * derived classes into a base class pointer.
 ***************************************************************************/
GEbounds* GEbounds::clone(void) const
{
    // Clone this image
    return new GEbounds(*this);
}


/***********************************************************************//**
 * @brief Insert energy interval
 *
 * @param[in] inx Index at which energy interval is to be inserted.
 * @param[in] emin Minimum energy of interval to be inserted.
 * @param[in] emax Maximum energy of interval to be inserted.
 *
 * Inserts an energy interval at a given position in the GEbounds array. This
 * method does not assure the energy ordering of the intervals, this has to
 * be done by the client that determines the appropriate value of inx.
 * Invalid parameters do not produce any exception, but are handled
 * transparently. If the interval is invalid (emax <= emin), no interval
 * will be inserted. If the index is out of the valid range, the index
 * will be adjusted to either the first or the last element.
 ***************************************************************************/
void GEbounds::insert_eng(int inx, const GEnergy& emin, const GEnergy& emax)
{
    // Continue only if energy interval is valid
    if (emax > emin) {

        // If inx is out of range then adjust it
        if (inx < 0)     inx = 0;
        if (inx > m_num) inx = m_num;
    
        // Allocate new intervals
        int      num = m_num+1;
        GEnergy* min = new GEnergy[num];
        GEnergy* max = new GEnergy[num];
        
        // Copy intervals before GTI to be inserted
        for (int i = 0; i < inx; ++i) {
            min[i] = m_min[i];
            max[i] = m_max[i];
        }
        
        // Insert interval
        min[inx] = emin;
        max[inx] = emax;

        // Copy intervals after GTI to be inserted
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
        
        // Set attributes
        m_num = num;
        set_attributes();
    
    } // endif: Energy interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Merge all overlapping energy intervals
 *
 * Merges all overlapping energy intervals into a single interval. This
 * method assumes that the intervals are ordered correctly by minimum
 * energy time. Note that the method does not actually reduce the memory
 * size but just updates the information on the number of elements in the
 * array.
 ***************************************************************************/
void GEbounds::merge_engs(void)
{
    // Find overlaps
    int i   = 0;
    int num = m_num;
    while (i < num-1) {

        // If energy interval overlaps with following one then merge both
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
        else
            i++;

    } // endwhile: there were still intervals to check
    
    // Update number of elements in object
    m_num = num;
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put energy boundaries in output stream
 *
 * @param[in] os Output stream into which the GEbounds will be dumped
 * @param[in] ebds GEbounds to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GEbounds& ebds)
{
    // Put GEbounds in stream
    os.precision(3);
    os << "=== GEbounds ===" << std::endl;
    os << " Number of boundaries ......: " << ebds.m_num << std::endl;
    os << " Energy range ..............: " << std::fixed
       << ebds.emin().MeV() << " - "
       << ebds.emax().MeV() << " MeV";

    // Return output stream
    return os;
}
