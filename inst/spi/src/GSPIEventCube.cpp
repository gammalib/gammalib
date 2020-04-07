/***************************************************************************
 *        GSPIEventCube.cpp - INTEGRAL/SPI event bin container class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIEventCube.hpp
 * @brief INTEGRAL/SPI event bin container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSPIEventCube.hpp"
#include "GSPITools.hpp"
#include "GSPIInstDir.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAXIS                                   "GSPIEventCube::naxis(int)"
#define G_MODEL_COUNTS                    "GSPIEventCube::model_counts(int&)"
#define G_SET_ENERGIES                        "GSPIEventCube::set_energies()"
#define G_SET_TIMES                              "GSPIEventCube::set_times()"
#define G_SET_BIN                              "GSPIEventCube::set_bin(int&)"
#define G_READ_FITS                             "GSPIEventCube::read(GFits&)"
#define G_READ_EBDS                   "GSPIEventCube::read_ebds(GFitsTable*)"
#define G_READ_PNT        "GSPIEventCube::read_pnt(GFitsTable*, GFitsTable*)"
#define G_READ_MODELS                    "GSPIEventCube::read_models(GFits&)"
#define G_PTID                                    "GSPIEventCube::ptid(int&)"

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
 * Constructs an empty INTEGRAL/SPI event cube.
 ***************************************************************************/
GSPIEventCube::GSPIEventCube(void) : GEventCube()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observation Group constructor
 *
 * @param[in] filename Observation Group FITS file name.
 *
 * Construct an INTEGRAL/SPI event cube from an Observation Group.
 ***************************************************************************/
GSPIEventCube::GSPIEventCube(const GFilename& filename) : GEventCube()
{
    // Initialise members
    init_members();

    // Load event cube
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube INTEGRAL/SPI event cube.
 ***************************************************************************/
GSPIEventCube::GSPIEventCube(const GSPIEventCube& cube) : GEventCube(cube)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSPIEventCube::~GSPIEventCube(void)
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
 * @param[in] cube INTEGRAL/SPI event cube.
 * @return INTEGRAL/SPI event cube.
 ***************************************************************************/
GSPIEventCube& GSPIEventCube::operator=(const GSPIEventCube& cube)
{
    // Execute only if object is not identical
    if (this != &cube) {

        // Copy base class members
        this->GEventCube::operator=(cube);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(cube);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Event bin access operator
 *
 * @param[in] index Event index [0,...,size()-1].
 * @return Pointer to event bin.
 *
 * Returns pointer to an event bin. Note that the returned pointer is in
 * fact always the same, but the method sets the pointers within the
 * event bin so that they point to the appropriate information.
 ***************************************************************************/
GSPIEventBin* GSPIEventCube::operator[](const int& index)
{
    // Set event bin
    set_bin(index);

    // Return pointer
    return (&m_bin);
}


/***********************************************************************//**
 * @brief Event bin access operator (const version)
 *
 * @param[in] index Event index [0,...,size()-1].
 * @return Const pointer to event bin.
 *
 * Returns pointer to an event bin. Note that the returned pointer is in
 * fact always the same, but the method sets the pointers within the
 * event bin so that they point to the appropriate information.
 ***************************************************************************/
const GSPIEventBin* GSPIEventCube::operator[](const int& index) const
{
    // Set event bin (circumvent const correctness)
    ((GSPIEventCube*)this)->set_bin(index);

    // Return pointer
    return (&m_bin);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear INTEGRAL/SPI event cube
 *
 * Clears INTEGRAL/SPI event cube by resetting all class members to an
 * initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GSPIEventCube::clear(void)
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
 * @brief Clone event cube
 *
 * @return Pointer to deep copy of INTEGRAL/SPI event cube.
 ***************************************************************************/
GSPIEventCube* GSPIEventCube::clone(void) const
{
    return new GSPIEventCube(*this);
}


/***********************************************************************//**
 * @brief Return number of bins in axis
 *
 * @param[in] axis Axis [0,...,dim()-1].
 * @return Number of bins in axis.
 *
 * @exception GException::out_of_range
 *            Axis is out of range.
 *
 * Returns the number of bins along a given event cube @p axis.
 ***************************************************************************/
int GSPIEventCube::naxis(const int& axis) const
{
    // Optionally check if the axis is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= dim()) {
        throw GException::out_of_range(G_NAXIS, "INTEGRAL/SPI event cube axis",
                                       axis, dim());
    }
    #endif

    // Setup const pointer array that points to relevant axis size
    const int* ptr[3] = {&m_num_pt, &m_num_det, &m_num_ebin};

    // Set number of bins dependent on axis
    int naxis = *ptr[axis];

    // Return result
    return naxis;
}


/***********************************************************************//**
 * @brief Load INTEGRAL/SPI event cube from Observation Group
 *
 * @param[in] filename Observation Group FITS file name.
 *
 * Loads data from an Observation Group FITS file into an INTEGRAL/SPI
 * event cube.
 ***************************************************************************/
void GSPIEventCube::load(const GFilename& filename)
{
    #pragma omp critical(GSPIEventCube_load)
    {
        // Clear object
        clear();

        // Open FITS file
        GFits fits(filename);

        // Read event cube from FITS file
        read(fits);

        // Close FITS file
        fits.close();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save INTEGRAL/SPI event cube into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing FITS file?
 *
 * Saves the INTEGRAL/SPI event cube into a FITS file.
 ***************************************************************************/
void GSPIEventCube::save(const GFilename& filename, const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write event cube into FITS file
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read INTEGRAL/SPI event cube from Observation Group FITS file
 *
 * @param[in] fits Observation Group FITS file.
 *
 * @exception GException::invalid_value
 *            Observation Group FITS file invalid.
 *
 * Reads an INTEGRAL/SPI event cube from an Observation Group FITS file.
 * The following extension are mandatory
 *
 *     "SPI.-EBDS-SET"
 *     "SPI.-OBS.-PNT"
 *     "SPI.-OBS.-GTI"
 *     "SPI.-OBS.-DSP"
 *     "SPI.-OBS.-DTI"
 *
 * Optional extensions are
 *
 *     "SPI.-SDET-SPE"
 *     "SPI.-BMOD-DSP"
 *
 ***************************************************************************/
void GSPIEventCube::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get table pointers
    const GFitsTable* ebds = gammalib::spi_hdu(fits, "SPI.-EBDS-SET");
    const GFitsTable* pnt  = gammalib::spi_hdu(fits, "SPI.-OBS.-PNT");
    const GFitsTable* gti  = gammalib::spi_hdu(fits, "SPI.-OBS.-GTI");
    const GFitsTable* dsp  = gammalib::spi_hdu(fits, "SPI.-OBS.-DSP");
    const GFitsTable* dti  = gammalib::spi_hdu(fits, "SPI.-OBS.-DTI");

    // Throw an exception if one of the mandatory HDUs is missing
    if (ebds == NULL) {
        std::string msg = "Extension \"SPI.-EBDS-SET\" not found in "
                          "Observation Group FITS file \""+
                          fits.filename().url()+"\". Please specify a "
                          "valid Observation Group.";
        throw GException::invalid_value(G_READ_FITS, msg);
    }
    if (pnt == NULL) {
        std::string msg = "Extension \"SPI.-OBS.-PNT\" not found in "
                          "Observation Group FITS file \""+
                          fits.filename().url()+"\". Please specify a "
                          "valid Observation Group.";
        throw GException::invalid_value(G_READ_FITS, msg);
    }
    if (gti == NULL) {
        std::string msg = "Extension \"SPI.-OBS.-GTI\" not found in "
                          "Observation Group FITS file \""+
                          fits.filename().url()+"\". Please specify a "
                          "valid Observation Group.";
        throw GException::invalid_value(G_READ_FITS, msg);
    }
    if (dsp == NULL) {
        std::string msg = "Extension \"SPI.-OBS.-DSP\" not found in "
                          "Observation Group FITS file \""+
                          fits.filename().url()+"\". Please specify a "
                          "valid Observation Group.";
        throw GException::invalid_value(G_READ_FITS, msg);
    }
    if (dti == NULL) {
        std::string msg = "Extension \"SPI.-OBS.-DTI\" not found in "
                          "Observation Group FITS file \""+
                          fits.filename().url()+"\". Please specify a "
                          "valid Observation Group.";
        throw GException::invalid_value(G_READ_FITS, msg);
    }

    // Determine dataspace dimensions from FITS tables
    m_num_pt   = dsp->integer("PT_NUM");
    m_num_det  = dsp->integer("DET_NUM");
    m_num_ebin = dsp->integer("EBIN_NUM");

    // Get number of sky and background models
    m_num_sky = gammalib::spi_num_hdus(fits, "SPI.-SDET-SPE");
    m_num_bgm = gammalib::spi_num_hdus(fits, "SPI.-BMOD-DSP");

    // Allocate data
    alloc_data();

    // Read energy boundaries
    read_ebds(ebds);

    // Read pointing information
    read_pnt(pnt, gti);

    // Read Good Time Intervals
    read_gti(gti);

    // Read dead time information
    read_dti(dti);

    // Read detector spectra
    read_dsp(dsp);

    // Read models
    read_models(fits);

    // Free HDU pointers
    if (ebds != NULL) delete ebds;
    if (pnt  != NULL) delete pnt;
    if (gti  != NULL) delete gti;
    if (dsp  != NULL) delete dsp;
    if (dti  != NULL) delete dti;

    // Prepare event bin
    init_bin();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write INTEGRAL/SPI event cube into FITS file.
 *
 * @param[in] file FITS file.
 *
 * Writes the INTEGRAL/SPI event cube into a FITS file.
 ***************************************************************************/
void GSPIEventCube::write(GFits& file) const
{
    // TODO: Here you have to implement the interface between the class and
    // the event cube FITS file.

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of events in cube
 *
 * @return Number of events in event cube.
 *
 * This method returns the number of events in the event cube rounded to the
 * nearest integer.
 ***************************************************************************/
int GSPIEventCube::number(void) const
{
    // Initialise result
    double number = 0.0;

    // Compute sum of all events
    if (m_dsp_size > 0) {
        for (int i = 0; i < m_dsp_size; ++i) {
            number += m_counts[i];
        }
    }

    // Return
    return int(number+0.5);
}


/***********************************************************************//**
 * @brief Return total ontime
 *
 * @return Total ontime.
 *
 * Returns the total ontime in the event cube. The total ontime is the sum
 * of the ontime of all pointings and detectors, divided by the number of
 * detectors that were active (determine by positive ontime values).
 ***************************************************************************/
double GSPIEventCube::ontime(void) const
{
    // Initialise result
    double ontime = 0.0;

    // Loop over all pointings
    for (int ipt = 0; ipt < m_num_pt; ++ipt) {

        // Initialise mean ontime and number of active detectors for this
        // pointing
        double mean_ontime    = 0.0;
        double num_active_det = 0.0;

        // Loop over all detectors
        for (int idet = 0; idet < m_num_det; ++idet) {

            // Get row index in table
            int irow = ipt * m_num_det + idet;

            // If detector has positive ontime then count it for the mean
            // ontime
            if (m_ontime[irow] > 0.0) {
                mean_ontime    += m_ontime[irow];
                num_active_det += 1.0;
            }

        } // endfor: looped over all detectors

        // Compute mean ontime for this pointing
        if (num_active_det > 0.0) {
            mean_ontime /= num_active_det;
        }
        else {
            mean_ontime = 0.0;
        }

        // Add mean ontime to total sum
        ontime += mean_ontime;

    } // endfor: looped over all pointings

    // Return ontime
    return ontime;
}


/***********************************************************************//**
 * @brief Return total livetime
 *
 * @return Total livetime.
 *
 * Returns the total livetime in the event cube. The total livetime is the
 * sum of the livetime of all pointings and detectors, divided by the number
 * of detectors that were active (determine by positive livetime values).
 ***************************************************************************/
double GSPIEventCube::livetime(void) const
{
    // Initialise result
    double livetime = 0.0;

    // Loop over all pointings
    for (int ipt = 0; ipt < m_num_pt; ++ipt) {

        // Initialise mean livetime and number of active detectors for this
        // pointing
        double mean_livetime  = 0.0;
        double num_active_det = 0.0;

        // Loop over all detectors
        for (int idet = 0; idet < m_num_det; ++idet) {

            // Get row index in table
            int irow = ipt * m_num_det + idet;

            // If detector has positive livetime then count it for the mean
            // livetime
            if (m_livetime[irow] > 0.0) {
                mean_livetime  += m_livetime[irow];
                num_active_det += 1.0;
            }

        } // endfor: looped over all detectors

        // Compute mean livetime for this pointing
        if (num_active_det > 0.0) {
            mean_livetime /= num_active_det;
        }
        else {
            mean_livetime = 0.0;
        }

        // Add mean livetime to total sum
        livetime += mean_livetime;

    } // endfor: looped over all pointings

    // Return livetime
    return livetime;
}


/***********************************************************************//**
 * @brief Return number of events in model
 *
 * @param[in] index Model index.
 * @return Number of events in event cube.
 *
 * @exception GException::out_of_range
 *            Invalid model index
 *
 * Returns the total number of counts in model.
 ***************************************************************************/
double GSPIEventCube::model_counts(const int& index) const
{
    // Initialise result
    double counts = 0.0;

    // Compute total number of models
    int num_models = m_num_sky + m_num_bgm;

    // Optionally check if the model index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= num_models) {
        throw GException::out_of_range(G_MODEL_COUNTS, "Invalid model index",
                                       index, num_models);
    }
    #endif

    // Compute sum of all events in model
    if (m_dsp_size > 0) {
        for (int i = 0; i < m_dsp_size; ++i) {
            counts += m_models[i*num_models + index];
        }
    }

    // Return
    return counts;
}


/***********************************************************************//**
 * @brief Return pointing identifier
 *
 * @param[in] ipt Pointing index.
 * @return Pointing identifier.
 *
 * @exception GException::out_of_range
 *            Invalid pointing index
 *
 * Returns the pointing identifier for the specified index.
 ***************************************************************************/
const std::string& GSPIEventCube::ptid(const int& ipt) const
{
    // Optionally check if the pointing index is valid
    #if defined(G_RANGE_CHECK)
    if (ipt < 0 || ipt >= m_num_pt) {
        throw GException::out_of_range(G_PTID, "Invalid pointing index",
                                       ipt, m_num_pt);
    }
    #endif

    // Return
    return (m_ptid[ipt]);
}


/***********************************************************************//**
 * @brief Print INTEGRAL/SPI event cube information
 *
 * @param[in] chatter Chattiness.
 * @return String containing INTEGRAL/SPI event cube information.
 ***************************************************************************/
std::string GSPIEventCube::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GSPIEventCube ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of events"));
        result.append(gammalib::str(number()));
        result.append("\n"+gammalib::parformat("Number of elements"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Pointings"));
        result.append(gammalib::str(m_num_pt));
        result.append("\n"+gammalib::parformat("Detectors"));
        result.append(gammalib::str(m_num_det));
        result.append("\n"+gammalib::parformat("Energy bins"));
        result.append(gammalib::str(m_num_ebin));
        result.append("\n"+gammalib::parformat("Sky models"));
        result.append(gammalib::str(m_num_sky));
        for (int i = 0; i < m_num_sky; ++i) {
            result.append("\n"+gammalib::parformat(" Model name "+gammalib::str(i+1)));
            result.append(m_modnames[i]);
            result.append("\n"+gammalib::parformat(" Number of events"));
            result.append(gammalib::str(model_counts(i)));
        }
        result.append("\n"+gammalib::parformat("Background models"));
        result.append(gammalib::str(m_num_bgm));
        for (int i = 0; i < m_num_bgm; ++i) {
            result.append("\n"+gammalib::parformat(" Model name "+gammalib::str(i+1)));
            result.append(m_modnames[i+m_num_sky]);
            result.append("\n"+gammalib::parformat(" Number of events"));
            result.append(gammalib::str(model_counts(i+m_num_sky)));
        }
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(gammalib::str(emin().MeV())+" - ");
        result.append(gammalib::str(emax().MeV())+" MeV");
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(ontime())+" s");
        result.append("\n"+gammalib::parformat("Livetime"));
        result.append(gammalib::str(livetime())+" s");
        result.append("\n"+gammalib::parformat("Time interval"));
        result.append(tstart().utc()+" - ");
        result.append(tstop().utc());

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
void GSPIEventCube::init_members(void)
{
    // Initialise members
    m_bin.clear();
    m_num_pt     = 0;
    m_num_det    = 0;
    m_num_ebin   = 0;
    m_num_sky    = 0;
    m_num_bgm    = 0;
    m_gti_size   = 0;
    m_dsp_size   = 0;
    m_model_size = 0;
    m_ontime     = NULL;
    m_livetime   = NULL;
    m_counts     = NULL;
    m_stat_err   = NULL;
    m_models     = NULL;
    m_size       = NULL;
    m_dir        = NULL;
    m_time       = NULL;
    m_energy     = NULL;
    m_ewidth     = NULL;
    m_ptid       = NULL;

    // Prepare event bin
    init_bin();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube INTEGRAL/SPI event cube.
 *
 * This method copies the class members from another event cube in the actual
 * object. It also prepares the event bin member that will be returned in
 * case of an operator access to the class.
 ***************************************************************************/
void GSPIEventCube::copy_members(const GSPIEventCube& cube)
{
    // Copy members. Note that the event bin is not copied as it will
    // be initialised later. The event bin serves just as a container of
    // pointers, hence we do not want to copy over the pointers from the
    // original class.
    m_num_pt     = cube.m_num_pt;
    m_num_det    = cube.m_num_det;
    m_num_ebin   = cube.m_num_ebin;
    m_num_sky    = cube.m_num_sky;
    m_num_bgm    = cube.m_num_bgm;
    m_gti_size   = cube.m_gti_size;
    m_dsp_size   = cube.m_dsp_size;
    m_model_size = cube.m_model_size;
    m_modnames   = cube.m_modnames;

    // Copy data
    if (m_num_ebin > 0) {
        m_energy = new GEnergy[m_num_ebin];
        m_ewidth = new GEnergy[m_num_ebin];
        for (int i = 0; i < m_num_ebin; ++i) {
            m_energy[i] = cube.m_energy[i];
            m_ewidth[i] = cube.m_ewidth[i];
        }
    }
    if (m_num_pt > 0) {
        m_ptid = new std::string[m_num_pt];
        for (int i = 0; i < m_num_pt; ++i) {
            m_ptid[i] = cube.m_ptid[i];
        }
    }
    if (m_gti_size > 0) {
        m_ontime   = new double[m_gti_size];
        m_livetime = new double[m_gti_size];
        m_time     = new GTime[m_gti_size];
        m_dir      = new GSPIInstDir[m_gti_size];
        for (int i = 0; i < m_gti_size; ++i) {
            m_ontime[i]   = cube.m_ontime[i];
            m_livetime[i] = cube.m_livetime[i];
            m_time[i]     = cube.m_time[i];
            m_dir[i]      = cube.m_dir[i];
        }
    }
    if (m_dsp_size > 0) {
        m_counts   = new double[m_dsp_size];
        m_stat_err = new double[m_dsp_size];
        m_size     = new double[m_dsp_size];
        for (int i = 0; i < m_dsp_size; ++i) {
            m_counts[i]   = cube.m_counts[i];
            m_stat_err[i] = cube.m_stat_err[i];
            m_size[i]     = cube.m_size[i];
        }
    }
    if (m_model_size > 0) {
        m_models = new double[m_model_size];
        for (int i = 0; i < m_model_size; ++i) {
            m_models[i] = cube.m_models[i];
        }
    }

    // Prepare event bin
    init_bin();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSPIEventCube::free_members(void)
{
    // Delete memory
    if (m_ontime   != NULL) delete [] m_ontime;
    if (m_livetime != NULL) delete [] m_livetime;
    if (m_counts   != NULL) delete [] m_counts;
    if (m_stat_err != NULL) delete [] m_stat_err;
    if (m_models   != NULL) delete [] m_models;
    if (m_size     != NULL) delete [] m_size;
    if (m_dir      != NULL) delete [] m_dir;
    if (m_time     != NULL) delete [] m_time;
    if (m_energy   != NULL) delete [] m_energy;
    if (m_ewidth   != NULL) delete [] m_ewidth;
    if (m_ptid     != NULL) delete [] m_ptid;

    // Set pointers to free
    m_ontime   = NULL;
    m_livetime = NULL;
    m_counts   = NULL;
    m_stat_err = NULL;
    m_models   = NULL;
    m_size     = NULL;
    m_dir      = NULL;
    m_time     = NULL;
    m_energy   = NULL;
    m_ewidth   = NULL;
    m_ptid     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate data
 ***************************************************************************/
void GSPIEventCube::alloc_data(void)
{
    // Make sure that data is free
    free_members();

    // Compute array sizes
    m_gti_size   = m_num_pt * m_num_det;
    m_dsp_size   = m_gti_size * m_num_ebin;
    m_model_size = m_dsp_size * (m_num_sky + m_num_bgm);

    // Allocate and initialise EBDS data
    if (m_num_ebin > 0) {
        m_energy = new GEnergy[m_num_ebin];
        m_ewidth = new GEnergy[m_num_ebin];
        for (int i = 0; i < m_num_ebin; ++i) {
            m_energy[i].clear();
            m_ewidth[i].clear();
        }
    }

    // Allocate and initialise pointing data
    if (m_num_pt > 0) {
        m_ptid = new std::string[m_num_pt];
        for (int i = 0; i < m_num_pt; ++i) {
            m_ptid[i].clear();
        }
    }

    // Allocate and initialise GTI data
    if (m_gti_size > 0) {
        m_ontime   = new double[m_gti_size];
        m_livetime = new double[m_gti_size];
        m_time     = new GTime[m_gti_size];
        m_dir      = new GSPIInstDir[m_gti_size];
        for (int i = 0; i < m_gti_size; ++i) {
            m_ontime[i]   = 0.0;
            m_livetime[i] = 0.0;
            m_time[i].clear();
            m_dir[i].clear();
        }
    }

    // Allocate and initialise DSP data
    if (m_dsp_size > 0) {
        m_counts   = new double[m_dsp_size];
        m_stat_err = new double[m_dsp_size];
        m_size     = new double[m_dsp_size];
        for (int i = 0; i < m_dsp_size; ++i) {
            m_counts[i]   = 0.0;
            m_stat_err[i] = 0.0;
            m_size[i]     = 0.0;
        }
    }

    // Allocate and initialise model data
    if (m_model_size > 0) {
        m_models = new double[m_model_size];
        for (int i = 0; i < m_model_size; ++i) {
            m_models[i] = 0.0;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read data from INTEGRAL/SPI "SPI.-EBDS-SET" extension
 *
 * @param[in] ebds Energy boundaries FITS table.
 *
 * @exception GException::invalid_value
 *            Incompatible number of energy bins encountered
 *
 * Reads data from an INTEGRAL/SPI "SPI.-EBDS-SET" extension. The method
 * sets up the m_energy array and the m_ebounds member.
 ***************************************************************************/
void GSPIEventCube::read_ebds(const GFitsTable* ebds)
{
    // Read energy boundaries
    m_ebounds.read(*ebds);

    // Throw an exception if the number of energy boundaries is not
    // consistent with DSP keyword
    if (m_ebounds.size() != m_num_ebin) {
        std::string msg = "Number of energy bins "+
                          gammalib::str(m_ebounds.size())+" found in "
                          "\"SPI.-EBDS-SET\" extension differes from value "+
                          gammalib::str(m_num_det)+" of \"DET_NUM\" keyword "
                          "in \"SPI.-OBS.-DSP\" extension. Please specify a "
                          "valid Observation Group.";
        throw GException::invalid_value(G_READ_EBDS, msg);
    }

    // Loop over all energy bins
    for (int iebin = 0; iebin < m_num_ebin; ++iebin) {

        // Store linear mean energy and bin width
        m_energy[iebin] = m_ebounds.emean(iebin);
        m_ewidth[iebin] = m_ebounds.ewidth(iebin);

    } // endfor: looped over all energy bins


    // Return
    return;
}


/***********************************************************************//**
 * @brief Read pointing information
 *
 * @param[in] pnt Pointing FITS table.
 * @param[in] gti GTI FITS table.
 *
 * @exception GException::invalid_value
 *            Incompatible PTID_SPI pointing identifiers encountered
 *
 * Reads pointing information from "SPI.-OBS.-PNT" and "SPI.-OBS.-GTI"
 * extensions. The method sets up the m_dir array.
 ***************************************************************************/
void GSPIEventCube::read_pnt(const GFitsTable* pnt, const GFitsTable* gti)
{
    // Get relevant columns
    const GFitsTableCol* pnt_ptid = (*pnt)["PTID_SPI"];
    const GFitsTableCol* ra_spix  = (*pnt)["RA_SPIX"];
    const GFitsTableCol* dec_spix = (*pnt)["DEC_SPIX"];
    const GFitsTableCol* gti_ptid = (*gti)["PTID_SPI"];
    const GFitsTableCol* det_id   = (*gti)["DET_ID"];

    // Loop over all pointings
    for (int ipt = 0; ipt < m_num_pt; ++ipt) {

        // Store pointing identifier
        m_ptid[ipt] = pnt_ptid->string(ipt);

        // Set pointing direction
        GSkyDir pnt_dir;
        pnt_dir.radec_deg(ra_spix->real(ipt), dec_spix->real(ipt));

        // Loop over all detectors
        for (int idet = 0; idet < m_num_det; ++idet) {

            // Get row index in table
            int irow = ipt * m_num_det + idet;

            // Throw an exception if the pointing identifier in the pointing
            // and the GTI extension is not the same
            if (pnt_ptid->string(ipt) != gti_ptid->string(irow)) {
                std::string msg = "PITD_SPI \""+pnt_ptid->string(ipt)+"\" in "
                                  "\"SPI.-OBS.-PNT\" differs from \""+
                                  gti_ptid->string(irow)+"\" in "
                                  "\"SPI.-OBS.-GTI\" extension for detector "+
                                  gammalib::str(det_id->real(irow))+". Please "
                                  "specify a valid Observation Group.";
                throw GException::invalid_value(G_READ_PNT, msg);
            }

            // Store pointing direction and detector identifier
            m_dir[irow].dir(pnt_dir);
            m_dir[irow].detid(det_id->real(irow));

        } // endfor: looped over all detectors

    } // endfor: looped over all pointings

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read data from INTEGRAL/SPI "SPI.-OBS.-GTI" extension
 *
 * @param[in] gti GTI FITS table.
 *
 * Reads data from an INTEGRAL/SPI "SPI.-OBS.-GTI" extension. The method
 * sets up the m_ontime array and the m_gti member.
 ***************************************************************************/
void GSPIEventCube::read_gti(const GFitsTable* gti)
{
    // Get relevant columns
    const GFitsTableCol* ontime = (*gti)["ONTIME"];
    const GFitsTableCol* tstart = (*gti)["TSTART"];
    const GFitsTableCol* tstop  = (*gti)["TSTOP"];

    // Loop over all pointings
    for (int ipt = 0; ipt < m_num_pt; ++ipt) {

        // Initialise minimum TSTART and maximum TSTOP for GTI (they should
        // all be identical and are anyways not used, but we want to have a
        // reasonable GTI object
        double t_start = 0.0;
        double t_stop  = 0.0;

        // Loop over all detectors
        for (int idet = 0; idet < m_num_det; ++idet) {

            // Get row index in table
            int irow = ipt * m_num_det + idet;

            // Store ontime
            m_ontime[irow] = ontime->real(irow);

            // Compute and store mean time
            double ijd   = 0.5 * (tstart->real(irow) + tstop->real(irow));
            m_time[irow] = gammalib::spi_ijd2time(ijd);

            // Update TSTART and TSTOP (only if ontime is > 0.0 since TSTART
            // and TSTOP may be zero for zero ontime)
            if (ontime->real(irow) > 0.0) {
                if (t_start == 0.0 || (tstart->real(irow) < t_start)) {
                    t_start = tstart->real(irow);
                }
                if (t_stop == 0.0 || (tstop->real(irow) > t_stop)) {
                    t_stop = tstop->real(irow);
                }
            }

        } // endfor: looped over all detectors

        // Append GTI (only if TSTART and TSTOP are not zero)
        if (t_start != 0.0 && t_stop != 0.0) {
            m_gti.append(gammalib::spi_ijd2time(t_start),
                         gammalib::spi_ijd2time(t_stop));
        }

    } // endfor: looped over all pointings

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read data from INTEGRAL/SPI "SPI.-OBS.-DTI" extension
 *
 * @param[in] dti Dead time information FITS table.
 *
 * Reads data from an INTEGRAL/SPI "SPI.-OBS.-DTI" extension. The method
 * sets up the m_ontime array and the m_gti member.
 ***************************************************************************/
void GSPIEventCube::read_dti(const GFitsTable* dti)
{
    // Get relevant columns
    const GFitsTableCol* livetime = (*dti)["LIVETIME"];

    // Loop over all pointings
    for (int ipt = 0; ipt < m_num_pt; ++ipt) {

        // Loop over all detectors
        for (int idet = 0; idet < m_num_det; ++idet) {

            // Get row index in table
            int irow = ipt * m_num_det + idet;

            // Store livetime
            m_livetime[irow] = livetime->real(irow);

        } // endfor: looped over all detectors

    } // endfor: looped over all pointings

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read data from INTEGRAL/SPI "SPI.-OBS.-DSP" extension
 *
 * @param[in] dsp DSP FITS table.
 *
 * Reads data from an INTEGRAL/SPI "SPI.-OBS.-DSP" extension. The method
 * sets up the m_counts, m_stat_err and m_size arrays.
 *
 * Note that the computation of the m_size arrays needs ontime and energy
 * width information that has been previously setup in the read_ebds() and
 * read_gti() methods.
 ***************************************************************************/
void GSPIEventCube::read_dsp(const GFitsTable* dsp)
{
    // Get relevant columns
    const GFitsTableCol* counts   = (*dsp)["COUNTS"];
    const GFitsTableCol* stat_err = (*dsp)["STAT_ERR"];

    // Loop over all pointings
    for (int ipt = 0; ipt < m_num_pt; ++ipt) {

        // Loop over all detectors
        for (int idet = 0; idet < m_num_det; ++idet) {

            // Get row index in table
            int irow  = ipt * m_num_det + idet;

            // Set start index in destination arrays
            int index = irow * m_num_ebin;

            // Copy energy bins
            for (int iebin = 0; iebin < m_num_ebin; ++iebin, ++index) {
                m_counts[index]   = counts->real(irow, iebin);
                m_stat_err[index] = stat_err->real(irow, iebin);
                m_size[index]     = m_ontime[irow] * m_ewidth[iebin].MeV();
            }

        } // endfor: looped over all detectors

    } // endfor: looped over all pointings

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read models from INTEGRAL/SPI Observation Group
 *
 * @param[in] fits Observation Group FITS file.
 ***************************************************************************/
void GSPIEventCube::read_models(const GFits& fits)
{
    // Clear model names
    m_modnames.clear();

    // Initialise model index
    int imodel = 0;

    // Compute total number of models
    int num_models = m_num_sky + m_num_bgm;

    // Loop over all sky models
    for (int i = 0; i < m_num_sky; ++i, ++imodel) {

        // Get FITS table
        const GFitsTable* model = gammalib::spi_hdu(fits, "SPI.-SDET-SPE", i+1);

        // Throw an exception if the FITS table is missing
        if (model == NULL) {
            std::string msg = "Extension \"SPI.-SDET-SPE\" version "+
                              gammalib::str(i+1)+" not found in Observation "
                              "Group FITS file \""+fits.filename().url()+
                              "\". Please specify a valid Observation Group.";
            throw GException::invalid_value(G_READ_MODELS, msg);
        }

        // Get model name
        std::string name = "MODEL" + gammalib::str(imodel+1);
        if (model->has_card("SOURCEID")) {
            name = model->string("SOURCEID");
        }
        m_modnames.push_back(name);

        // Get relevant columns
        const GFitsTableCol* counts = (*model)["COUNTS"];

        // Loop over all pointings
        for (int ipt = 0; ipt < m_num_pt; ++ipt) {

            // Loop over all detectors
            for (int idet = 0; idet < m_num_det; ++idet) {

                // Get row index in model table
                int irow  = ipt * m_num_det + idet;

                // Set start index in destination array
                int index = irow * m_num_ebin;

                // Copy energy bins
                for (int iebin = 0; iebin < m_num_ebin; ++iebin, ++index) {
                    m_models[index*num_models + imodel] = counts->real(irow, iebin);
                }

            } // endfor: looped over all detectors

        } // endfor: looped over all pointings

    } // endfor: looped over all sky models

    // Loop over all background models
    for (int i = 0; i < m_num_bgm; ++i, ++imodel) {

        // Get FITS table
        const GFitsTable* model = gammalib::spi_hdu(fits, "SPI.-BMOD-DSP", i+1);

        // Throw an exception if the FITS table is missing
        if (model == NULL) {
            std::string msg = "Extension \"SPI.-BMOD-DSP\" version "+
                              gammalib::str(i+1)+" not found in Observation "
                              "Group FITS file \""+fits.filename().url()+
                              "\". Please specify a valid Observation Group.";
            throw GException::invalid_value(G_READ_MODELS, msg);
        }

        // Get model name
        std::string name = "MODEL" + gammalib::str(imodel+1);
        if (model->has_card("BKGNAME")) {
            name = model->string("BKGNAME");
        }
        m_modnames.push_back(name);

        // Get relevant columns
        const GFitsTableCol* counts = (*model)["COUNTS"];

        // Loop over all pointings
        for (int ipt = 0; ipt < m_num_pt; ++ipt) {

            // Loop over all detectors
            for (int idet = 0; idet < m_num_det; ++idet) {

                // Get row index in model table
                int irow  = ipt * m_num_det + idet;

                // Set start index in destination array
                int index = irow * m_num_ebin;

                // Copy energy bins
                for (int iebin = 0; iebin < m_num_ebin; ++iebin, ++index) {
                    m_models[index*num_models + imodel] = counts->real(irow, iebin);
                }

            } // endfor: looped over all detectors

        } // endfor: looped over all pointings

    } // endfor: looped over all background models

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise event bin
 *
 * Initialises the event bin. All fixed content is set here, the content
 * that depends on the bin index is set by the set_bin() method which is
 * called before any event bin access.
 ***************************************************************************/
void GSPIEventCube::init_bin(void)
{
    // Prepare event bin
    m_bin.free_members();
    m_bin.m_index      = 0;     //!< Set by set_bin method
    m_bin.m_idir       = 0;     //!< Set by set_bin method
    m_bin.m_iebin      = 0;     //!< Set by set_bin method
    m_bin.m_num_models = m_num_sky + m_num_bgm;
    m_bin.m_dir        = NULL;  //!< Set by set_bin method
    m_bin.m_time       = NULL;  //!< Set by set_bin method
    m_bin.m_energy     = NULL;  //!< Set by set_bin method
    m_bin.m_counts     = NULL;  //!< Set by set_bin method
    m_bin.m_ontime     = NULL;  //!< Set by set_bin method
    m_bin.m_size       = NULL;  //!< Set by set_bin method
    m_bin.m_models     = NULL;  //!< Set by set_bin method

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set event bin
 *
 * @param[in] index Event index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Event index is outside valid range.
 *
 * Sets the pointers of the event bin to the event cube cell that corresponds
 * to the specified @p index.
 *
 * The event bin is in fact physically stored in the event cube, and only a
 * single event bin is indeed allocated. This method sets up the pointers in
 * the event bin so that a client can easily access the information of
 * individual bins as if they were stored in an array.
 ***************************************************************************/
void GSPIEventCube::set_bin(const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET_BIN, index, size());
    }
    #endif

    // Set indices
    m_bin.m_index = index;
    m_bin.m_idir  = index / m_num_ebin;
    m_bin.m_iebin = index % m_num_ebin;

    // Set GSPIEventBin pointers
    m_bin.m_dir    = m_dir    + m_bin.m_idir;
    m_bin.m_time   = m_time   + m_bin.m_idir;
    m_bin.m_energy = m_energy + m_bin.m_iebin;
    m_bin.m_counts = m_counts + m_bin.m_index;
    //m_bin.m_ontime = m_ontime + m_bin.m_idir;
    m_bin.m_ontime = m_livetime + m_bin.m_idir; // Use livetime instead of ontime?
    m_bin.m_size   = m_size   + m_bin.m_index;
    m_bin.m_models = m_models + m_bin.m_index * m_bin.m_num_models;

    // Return
    return;
}
