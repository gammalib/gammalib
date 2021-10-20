/***************************************************************************
 *        GXXXEventCube.cpp - [INSTRUMENT] event bin container class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GXXXEventCube.hpp
 * @brief [INSTRUMENT] event bin container class implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GXXXEventCube.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAXIS                                   "GXXXEventCube::naxis(int)"
#define G_SET_ENERGIES                        "GXXXEventCube::set_energies()"
#define G_SET_TIMES                              "GXXXEventCube::set_times()"
#define G_SET_BIN                              "GXXXEventCube::set_bin(int&)"

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
 * Constructs an empty [INSTRUMENT] event cube.
 ***************************************************************************/
GXXXEventCube::GXXXEventCube(void) : GEventCube()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename [INSTRUMENT] event cube FITS filename.
 *
 * Construct a [INSTRUMENT] event cube by loading it from a FITS file.
 ***************************************************************************/
GXXXEventCube::GXXXEventCube(const GFilename& filename) : GEventCube()
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
 * @param[in] cube [INSTRUMENT] event cube.
 ***************************************************************************/
GXXXEventCube::GXXXEventCube(const GXXXEventCube& cube) : GEventCube(cube)
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
GXXXEventCube::~GXXXEventCube(void)
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
 * @param[in] cube [INSTRUMENT] event cube.
 * @return [INSTRUMENT] event cube.
 ***************************************************************************/
GXXXEventCube& GXXXEventCube::operator=(const GXXXEventCube& cube)
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
GXXXEventBin* GXXXEventCube::operator[](const int& index)
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
const GXXXEventBin* GXXXEventCube::operator[](const int& index) const
{
    // Set event bin (circumvent const correctness)
    const_cast<GXXXEventCube*>(this)->set_bin(index);

    // Return pointer
    return (&m_bin);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear [INSTRUMENT] event cube
 *
 * Clears [INSTRUMENT] event cube by resetting all class members to an
 * initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GXXXEventCube::clear(void)
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
 * @return Pointer to deep copy of [INSTRUMENT] event cube.
 ***************************************************************************/
GXXXEventCube* GXXXEventCube::clone(void) const
{
    return new GXXXEventCube(*this);
}


/***********************************************************************//**
 * @brief Return number of bins in event cube
 *
 * @return Number of bins in event cube.
 *
 * @todo Implement method.
 ***************************************************************************/
int GXXXEventCube::size(void) const
{
    // Compute number of bins
    // TODO: Implement computation
    int nbins = 0;

    // Return number of bins
    return nbins;
}


/***********************************************************************//**
 * @brief Return dimension of event cube
 *
 * @return Number of dimensions in event cube.
 *
 * @todo Implement method.
 ***************************************************************************/
int GXXXEventCube::dim(void) const
{
    // Compute dimension
    // TODO: Implement computation
    int dim = 0;

    // Return dimension
    return dim;
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
int GXXXEventCube::naxis(const int& axis) const
{
    // Optionally check if the axis is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= dim()) {
        throw GException::out_of_range(G_NAXIS, "[INSTRUMENT] event cube axis",
                                       axis, dim());
    }
    #endif

    // Set result
    // TODO: Implement computation
    int naxis = 0;

    // Return result
    return naxis;
}


/***********************************************************************//**
 * @brief Load [INSTRUMENT] event cube from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * The method clears the object before loading, thus any events residing in
 * the object before loading will be lost.
 ***************************************************************************/
void GXXXEventCube::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read event cube from FITS file
    read(fits);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save [INSTRUMENT] event cube into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing FITS file?
 *
 * Saves the [INSTRUMENT] event cube into a FITS file.
 ***************************************************************************/
void GXXXEventCube::save(const GFilename& filename, const bool& clobber) const
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
 * @brief Read [INSTRUMENT] event cube from FITS file
 *
 * @param[in] fits FITS file.
 *
 * Reads a [INSTRUMENT] event cube from a FITS file.
 *
 * @todo Implement method.
 ***************************************************************************/
void GXXXEventCube::read(const GFits& fits)
{
    // Clear object
    clear();

    // TODO: Here you have to implement the interface between the class and
    // the event cube FITS file.

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write [INSTRUMENT] event cube into FITS file.
 *
 * @param[in] file FITS file.
 *
 * Writes the [INSTRUMENT] event cube into a FITS file.
 ***************************************************************************/
void GXXXEventCube::write(GFits& file) const
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
 *
 * @todo Implement method.
 ***************************************************************************/
int GXXXEventCube::number(void) const
{
    // Initialise result
    double number = 0.0;

    // TODO: Here you have to implement the computation of the number of
    // events in the cube.

    // Return
    return int(number+0.5);
}


/***********************************************************************//**
 * @brief Print [INSTRUMENT] event cube information
 *
 * @param[in] chatter Chattiness.
 * @return String containing [INSTRUMENT] event cube information.
 ***************************************************************************/
std::string GXXXEventCube::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GXXXEventCube ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of events"));
        result.append(gammalib::str(number()));
        result.append("\n"+gammalib::parformat("Number of elements"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(gammalib::str(emin().MeV())+" - ");
        result.append(gammalib::str(emax().MeV())+" MeV");
        result.append("\n"+gammalib::parformat("Mean energy"));
        result.append(m_energy.print(chatter));
        result.append("\n"+gammalib::parformat("Energy bin width"));
        result.append(m_ewidth.print(chatter));
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(m_ontime)+" s");
        result.append("\n"+gammalib::parformat("Time interval"));
        result.append(gammalib::str(tstart().jd())+" - ");
        result.append(gammalib::str(tstop().jd())+" Julian days");
        result.append("\n"+gammalib::parformat("Mean time"));
        result.append(m_time.print(chatter));

        // Append additional information
        // TODO: Add code to append any additional information that might
        // be relevant.

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
void GXXXEventCube::init_members(void)
{
    // Initialise members
    m_bin.clear();
    m_energy.clear();
    m_time.clear();
    m_polarization.clear();
    m_ontime = 0.0;
    m_ewidth.clear();

    // Prepare event bin
    init_bin();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube [INSTRUMENT] event cube.
 *
 * This method copies the class members from another event cube in the actual
 * object. It also prepares the event bin member that will be returned in
 * case of an operator access to the class.
 ***************************************************************************/
void GXXXEventCube::copy_members(const GXXXEventCube& cube)
{
    // Copy members. Note that the event bin is not copied as it will
    // be initialised later. The event bin serves just as a container of
    // pointers, hence we do not want to copy over the pointers from the
    // original class.
    m_energy       = cube.m_energy;
    m_time         = cube.m_time;
    m_polarization = cube.m_polarization;
    m_ontime       = cube.m_ontime;
    m_ewidth       = cube.m_ewidth;

    // Prepare event bin
    init_bin();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXXXEventCube::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise event bin
 *
 * Initialises the event bin. The event bin is cleared and all fixed pointers
 * are set. Only the m_counts and the m_solidangle member of the event bin
 * will be set to NULL, but these will be set by the set_bin() method
 * which is called before any event bin access.
 ***************************************************************************/
void GXXXEventCube::init_bin(void)
{
    // Prepare event bin
    m_bin.free_members();
    m_bin.m_counts       = NULL;            //!< Will be set by set_bin method
    m_bin.m_dir          = &m_dir;          //!< Content will be set by set_bin method
    m_bin.m_time         = &m_time;         //!< Fixed content
    m_bin.m_polarization = &m_polarization; //!< Fixed content
    m_bin.m_energy       = &m_energy;       //!< Fixed content

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set event bin
 *
 * @param[in] index Event index [0,...,size()[.
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
void GXXXEventCube::set_bin(const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET_BIN, "Event index", index, size());
    }
    #endif

    // Set bin index
    m_bin.m_index = index;

    // TODO: Set here the pointers of the GXXXEventBin to the event cube
    // cell that corresponds to the specified index.

    // Return
    return;
}
