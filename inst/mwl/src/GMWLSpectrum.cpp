/***************************************************************************
 *           GMWLSpectrum.cpp  -  Multi-wavelength spectrum class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GMWLSpectrum.cpp
 * @brief GMWLSpectrum class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GLog.hpp"
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GException.hpp"
#include "GMWLException.hpp"
#include "GMWLSpectrum.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_POINTER                                "GMWLSpectrum::pointer(int)"
#define G_LOAD_FITS              "GMWLSpectrum::load_fits(std::string&, int)"
#define G_READ_FITS                    "GMWLSpectrum::read_fits(GFitsTable*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates instance of an undefined spectrum.
 ***************************************************************************/
GMWLSpectrum::GMWLSpectrum(void) : GEventCube()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Filename.
 *
 * Creates instance from file.
 ***************************************************************************/
GMWLSpectrum::GMWLSpectrum(const std::string& filename) : GEventCube()
{
    // Initialise class members for clean destruction
    init_members();

    // Load spectrum
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] spec Spectrum.
 *
 * Creates instance by copying an existing spectrum.
 ***************************************************************************/
GMWLSpectrum::GMWLSpectrum(const GMWLSpectrum& spec) : GEventCube(spec)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(spec);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroy instance.
 ***************************************************************************/
GMWLSpectrum::~GMWLSpectrum(void)
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
 * @param[in] spec Spectrum.
 *
 * Copies spectrum into the instance.
 ***************************************************************************/
GMWLSpectrum& GMWLSpectrum::operator= (const GMWLSpectrum& spec)
{
    // Execute only if object is not identical
    if (this != &spec) {

        // Copy base class members
        this->GEventCube::operator=(spec);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(spec);

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
 * @brief Clear object
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GMWLSpectrum::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GEventCube::free_members();
    this->GEvents::free_members();

    // Initialise members
    this->GEvents::init_members();
    this->GEventCube::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GMWLSpectrum* GMWLSpectrum::clone(void) const
{
    return new GMWLSpectrum(*this);
}


/***********************************************************************//**
 * @brief Return number of points in spectrum
 ***************************************************************************/
int GMWLSpectrum::size(void) const
{
    // Return
    return (m_data.size());
}


/***********************************************************************//**
 * @brief Load spectrum
 *
 * @param[in] filename File name
 *
 * This method loads a spectrum from a variety of file types. The method
 * analyses the format of the file that is presented and choses then the
 * appropriate method to load the specific format. The following file
 * formats are supported: FITS, TBW ...
 *
 * @todo Not yet implemented.
 ***************************************************************************/
void GMWLSpectrum::load(const std::string& filename)
{
    // Load from FITS file
    load_fits(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load spectrum from FITS file
 *
 * @param[in] filename File name.
 * @param[in] extno Extension number of table (optional).
 *
 * @exception GMWLException::file_open_error
 *            No table found in file.
 *
 * Loads the spectrum from a FITS table found in the specified extension.
 * In no extension number if specified (or if extno=0) then the spectrum
 * is loaded from the first table extension that is found in the file.
 ***************************************************************************/
void GMWLSpectrum::load_fits(const std::string& filename, int extno)
{
    // Clear object
    clear();

    // Open FITS file
    GFits file(filename);

    // If the extension number is 0 then load first FITS table in file.
    if (extno == 0) {
        for (int i = 0; i < file.size(); ++i) {
            if (file.hdu(i)->exttype() == GFitsHDU::HT_ASCII_TABLE ||
                file.hdu(i)->exttype() == GFitsHDU::HT_BIN_TABLE) {
                extno = i;
                break;
            }
        }
    }

    // If we found no table then throw an exception
    if (extno == 0)
        throw GMWLException::file_open_error(G_LOAD_FITS, filename,
                                             "No table found in file.");

    // Get table pointer
    GFitsTable* table = file.table(extno);

    // Read spectrum from table
    read_fits(table);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load spectrum from FITS file
 *
 * @param[in] filename File name.
 * @param[in] extname Extension name.
 *
 * Loads the spectrum from a FITS table with a given extension name.
 ***************************************************************************/
void GMWLSpectrum::load_fits(const std::string& filename,
                             const std::string& extname)
{
    // Clear object
    clear();

    // Open FITS file
    GFits file(filename);

    // Get table pointer
    GFitsTable* table = file.table(extname);

    // Read spectrum from table
    read_fits(table);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to spectral point
 *
 * @param[in] index Spectral point index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Spectral point index not in valid range.
 ***************************************************************************/
GMWLDatum* GMWLSpectrum::pointer(int index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_data.size())
        throw GException::out_of_range(G_POINTER, index, 0, m_data.size()-1);
    #endif

    // Return pointer
    return &(m_data[index]);
}


/***********************************************************************//**
 * @brief Return number of point in spectrum
 *
 * @todo Obsolete, but still needed as it is pure virtual in the base class.
 * Better use size().
 ***************************************************************************/
int GMWLSpectrum::number(void) const
{
    // Return
    return (size());
}


/***********************************************************************//**
 * @brief Print spectrum
 ***************************************************************************/
std::string GMWLSpectrum::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GMWLObservation ===\n");
    result.append(parformat("Number of points")+str(size()));

    // Append spectral points
    for (int i = 0; i < size(); ++i) {

        // Build energy string
        std::string energy = m_data[i].m_eng.print();
        if (m_data[i].m_eng_err.MeV() > 0.0)
            energy += " +/- "+m_data[i].m_eng_err.print();

        // Build flux string
        std::string flux = str(m_data[i].m_flux);
        if (m_data[i].m_flux_err > 0.0)
            flux += " +/- "+str(m_data[i].m_flux_err);
        flux += " ph/cm2/s/MeV";

        // Append to string
        result.append("\n"+parformat(energy));
        result.append(flux);

    } // endfor: looped over spectral points

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
void GMWLSpectrum::init_members(void)
{
    // Initialise members
    m_data.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] spec Instance to be copied.
 ***************************************************************************/
void GMWLSpectrum::copy_members(const GMWLSpectrum& spec)
{
    // Copy members
    m_data = spec.m_data;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMWLSpectrum::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read spectrum from FITS file
 *
 * @param[in] table File table.
 *
 * @exception GMWLException::bad_file_format
 *            Table has invalid format
 *
 * Read spectrum from FITS table. The table is expected to be in one of the
 * three following formats:
 * 2 columns: energy, flux
 * 3 columns: energy, flux, e_flux
 * 4 columns or more:  energy, e_energy, flux, e_flux, ...
 *
 * @todo Investigate whether we can exploit UCDs for identifying the correct
 * columns or for determining the units.
 ***************************************************************************/
void GMWLSpectrum::read_fits(const GFitsTable* table)
{
    // Reset spectrum
    m_data.clear();

    // Initialise column pointers columns
    GFitsTableCol* c_energy     = NULL;
    GFitsTableCol* c_energy_err = NULL;
    GFitsTableCol* c_flux       = NULL;
    GFitsTableCol* c_flux_err   = NULL;

    // Extract column pointers
    if (table->ncols() == 2) {
        c_energy = ((GFitsTable*)table)->column(0);
        c_flux   = ((GFitsTable*)table)->column(1);
    }
    else if (table->ncols() == 3) {
        c_energy   = ((GFitsTable*)table)->column(0);
        c_flux     = ((GFitsTable*)table)->column(1);
        c_flux_err = ((GFitsTable*)table)->column(2);
    }
    else if (table->ncols() > 3) {
        c_energy     = ((GFitsTable*)table)->column(0);
        c_energy_err = ((GFitsTable*)table)->column(1);
        c_flux       = ((GFitsTable*)table)->column(2);
        c_flux_err   = ((GFitsTable*)table)->column(3);
    }
    else {
        throw GMWLException::bad_file_format(G_READ_FITS,
                             "At least 2 columns are expected is table \""+
                              table->extname()+"\".");
    }

    // Read spectral points and add to spectrum
    for (int i = 0; i < table->nrows(); ++i) {
        GMWLDatum datum;
        if (c_energy     != NULL) datum.m_eng.MeV(c_energy->real(i));
        if (c_energy_err != NULL) datum.m_eng_err.MeV(c_energy_err->real(i));
        if (c_flux       != NULL) datum.m_flux = c_flux->real(i);
        if (c_flux_err   != NULL) datum.m_flux_err = c_flux_err->real(i);
        m_data.push_back(datum);
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GMWLSpectrum friends                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] spec Spectrum.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GMWLSpectrum& spec)
{
     // Write spectrum in output stream
    os << spec.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] spec Spectrum.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GMWLSpectrum& spec)
{
    // Write spectrum into logger
    log << spec.print();

    // Return logger
    return log;
}
