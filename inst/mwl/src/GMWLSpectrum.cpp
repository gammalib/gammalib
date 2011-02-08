/***************************************************************************
 *           GMWLSpectrum.cpp  -  Multi-wavelength spectrum class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @brief Multi-wavelength spectrum class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GEnergy.hpp"
#include "GException.hpp"
#include "GMWLException.hpp"
#include "GMWLSpectrum.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR                           "GMWLSpectrum::operator[](int&)"
#define G_READ                             "GMWLSpectrum::read(GFits&, int&)"
#define G_READ_FITS                    "GMWLSpectrum::read_fits(GFitsTable*)"
#define G_CONV_ENERGY       "GMWLSpectrum::conv_energy(double&,std::string&)"
#define G_CONV_FLUX  "GMWLSpectrum::conv_flux(GEnergy&,double&,std::string&)"

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
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name.
 *
 * Creates instance from file.
 ***************************************************************************/
GMWLSpectrum::GMWLSpectrum(const std::string& filename) : GEventCube()
{
    // Initialise members
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
GMWLSpectrum& GMWLSpectrum::operator=(const GMWLSpectrum& spec)
{
    // Execute only if object is not identical
    if (this != &spec) {

        // Copy base class members
        this->GEventCube::operator=(spec);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(spec);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Spectral point access operator
 *
 * @param[in] index Spectral point index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Spectral point index outside valid range.
 ***************************************************************************/
GMWLDatum* GMWLSpectrum::operator[](const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
    #endif

    // Return pointer
    return &(m_data[index]);
}


/***********************************************************************//**
 * @brief Spectral point access operator (const version)
 *
 * @param[in] index Spectral point index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Spectral point index outside valid range.
 ***************************************************************************/
const GMWLDatum* GMWLSpectrum::operator[](const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
    #endif

    // Return pointer
    return &(m_data[index]);
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
 * @brief Load spectrum
 *
 * @param[in] filename File name.
 *
 * This method loads a spectrum from a variety of file types. The method
 * analyses the format of the file that is presented and choses then the
 * appropriate method to load the specific format. The following file
 * formats are supported: FITS, TBW ...
 *
 * @todo So far only FITS file support is implemented.
 ***************************************************************************/
void GMWLSpectrum::load(const std::string& filename)
{
    // Clear object
    clear();

    // Open FITS file
    GFits file(filename);

    // Read spectrum
    read(file);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load spectrum from FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] extname FITS extension name.
 *
 * This method loads a spectrum from the specified extension from a FITS
 * file.
 ***************************************************************************/
void GMWLSpectrum::load(const std::string& filename, const std::string& extname)
{
    // Clear object
    clear();

    // Open FITS file
    GFits file(filename);

    // Read spectrum from FITS file
    read(file, extname);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load spectrum from FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] extno FITS extension number.
 *
 * This method loads a spectrum from the specified extension number from a
 * FITS file.
 ***************************************************************************/
void GMWLSpectrum::load(const std::string& filename, int extno)
{
    // Clear object
    clear();

    // Open FITS file
    GFits file(filename);

    // Read spectrum from FITS file
    read(file, extno);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save spectrum
 *
 * @param[in] filename File name.
 * @param[in] clobber Overwrite existing file (default=false).
 *
 * @todo To be implemented.
 ***************************************************************************/
void GMWLSpectrum::save(const std::string& filename, bool clobber) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read spectrum from FITS file
 *
 * @param[in] file FITS file.
 *
 * Read spectrum from first extension in FITS file.
 ***************************************************************************/
void GMWLSpectrum::read(const GFits& file)
{
    // Clear object
    clear();

    // Read spectrum from first extension in FITS file
    read(file, 0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read spectrum from FITS file
 *
 * @param[in] file FITS file.
 * @param[in] extname FITS extension name.
 ***************************************************************************/
void GMWLSpectrum::read(const GFits& file, const std::string& extname)
{
    // Clear object
    clear();

    // Get table pointer
    GFitsTable* table = file.table(extname);

    // Read spectrum from table
    read_fits(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read spectrum from FITS file
 *
 * @param[in] file FITS file.
 * @param[in] extno Extension number of spectrum.
 *
 * @exception GMWLException::file_open_error
 *            No table found in file.
 *
 * Read the spectrum from a FITS table found in the specified extension.
 * In no extension number if specified (or if extno=0) then the spectrum
 * is loaded from the first table extension that is found in the file.
 ***************************************************************************/
void GMWLSpectrum::read(const GFits& file, int extno)
{
    // Clear object
    clear();

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
        throw GMWLException::file_open_error(G_READ, file.name(),
                                             "No table found in file.");

    // Get table pointer
    GFitsTable* table = file.table(extno);

    // Read spectrum from table
    read_fits(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write spectrum into FITS file
 *
 * @param[in] file FITS file.
 *
 * @todo To be implemented.
 ***************************************************************************/
void GMWLSpectrum::write(GFits& file) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of points in spectrum
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
    result.append("=== GMWLSpectrum ===");
    result.append("\n"+parformat("Telescope")+m_telescope);
    result.append("\n"+parformat("Instrument")+m_instrument);
    result.append("\n"+parformat("Number of points")+str(size()));
    result.append("\n"+parformat("Time interval"));
    if (gti().size() > 0)
        result.append(str(tstart().met())+" - "+str(tstop().met()));
    else
        result.append("not defined");
    result.append("\n"+parformat("Energy range"));
    if (ebounds().size() > 0)
        result.append(emin().print()+" - "+emax().print());
    else
        result.append("not defined");

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
    m_telescope.clear();
    m_instrument.clear();
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
    m_telescope  = spec.m_telescope;
    m_instrument = spec.m_instrument;
    m_data       = spec.m_data;

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
 * @brief Set energy boundaries
 ***************************************************************************/
void GMWLSpectrum::set_ebounds(void)
{
    // Clear energy boundaries
    m_ebounds.clear();

    // Continue only if we have data
    if (size() > 0) {

        // Extract energy boundaries from spectrum
        GEnergy emin = m_data[0].energy();
        GEnergy emax = m_data[0].energy();
        for (int i = 0; i < m_data.size(); ++i) {
            if (m_data[i].energy() < emin)
                emin = m_data[i].energy();
            if (m_data[i].energy() > emax)
                emax = m_data[i].energy();
        }

        // Set energy boundaries
        m_ebounds.append(emin, emax);

    } // endif: we had data

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
        if (c_energy     != NULL) 
            datum.m_eng = conv_energy(c_energy->real(i), c_energy->unit());
        if (c_energy_err != NULL)
            datum.m_eng_err = conv_energy(c_energy_err->real(i), c_energy->unit());
        if (c_flux       != NULL)
            datum.m_flux = conv_flux(datum.m_eng, c_flux->real(i), c_flux->unit());
        if (c_flux_err   != NULL)
            datum.m_flux_err = conv_flux(datum.m_eng, c_flux_err->real(i), c_flux_err->unit());
        m_data.push_back(datum);
    }

    // Get telescope name
    try {
        m_telescope = table->string("TELESCOP");
    }
    catch (GException::fits_key_not_found &e) {
        m_telescope = "unknown";
    }

    // Get instrument name
    try {
        m_instrument = table->string("INSTRUME");
    }
    catch (GException::fits_key_not_found &e) {
        m_instrument = "unknown";
    }

    // Set energy boundaries
    set_ebounds();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert value into energy
 *
 * @param[in] energy Energy value.
 * @param[in] unit Unit of value.
 *
 * @exception GMWLException::invalid_unit
 *            Invalid unit string encountered
 *
 * Converts an energy value into a GEnergy object based on the specified
 * units. The following units are supported (case insensitive):
 * erg, keV, MeV, GeV, and TeV.
 ***************************************************************************/
GEnergy GMWLSpectrum::conv_energy(const double& energy, const std::string& unit)
{
    // Initialise energy
    GEnergy result;

    // Convert unit string to upper base without any leading/trailing
    // whitespace
    std::string str_unit = strip_whitespace(toupper(unit));

    // High-energy units
    if (str_unit == "KEV")
        result.keV(energy);
    else if (str_unit == "MEV")
        result.MeV(energy);
    else if (str_unit == "GEV")
        result.GeV(energy);
    else if (str_unit == "TEV")
        result.TeV(energy);

    // Other units
    else if (str_unit == "ERG")
        result.erg(energy);
    
    // ... otherwise throw exception
    else
        throw GMWLException::invalid_unit(G_CONV_ENERGY, unit);

    // Return energy
    return result;
}


/***********************************************************************//**
 * @brief Convert value into flux
 *
 * @param[in] energy Energy at which flux is given.
 * @param[in] flux Flux value.
 * @param[in] unit Unit of value.
 *
 * @exception GMWLException::invalid_unit
 *            Invalid unit string encountered
 *
 * Converts a flux value into units of ph/cm2/s/MeV based on the specified
 * units. The following units are supported (case insensitive):
 * ph/cm2/s/MeV, ph/s/cm2/MeV, erg/cm2/s and erg/s/cm2.
 ***************************************************************************/
double GMWLSpectrum::conv_flux(const GEnergy& energy, const double& flux,
                               const std::string& unit)
{
    // Initialise energy
    double result;

    // Convert unit string to upper base without any leading/trailing
    // whitespace
    std::string str_unit = strip_whitespace(toupper(unit));

    // High-energy units
    if (str_unit == "PH/CM2/S/MEV" || str_unit == "PH/S/CM2/MEV")
        result = flux;
    else if (str_unit == "ERG/CM2/S" || str_unit == "ERG/S/CM2")
        result = (erg2MeV*flux) / (energy.MeV()*energy.MeV());

    // ... otherwise throw exception
    else
        throw GMWLException::invalid_unit(G_CONV_FLUX, unit);

    // Return energy
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
