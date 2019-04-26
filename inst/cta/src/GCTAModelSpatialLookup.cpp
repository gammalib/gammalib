/***************************************************************************
 *         GCTAModelSpatialLookup.cpp - Spatial lookup table model         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialLookup.cpp
 * @brief Spatial lookup table model interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GEbounds.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GObservations.hpp"
#include "GCTAObservation.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAModelSpatialLookup.hpp"
#include "GCTAModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelSpatialLookup   g_cta_spatial_lookup_seed;
const GCTAModelSpatialRegistry g_cta_spatial_lookup_registry(&g_cta_spatial_lookup_seed);

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCTOR   "GCTAModelSpatialLookup(double&, double&, GEbounds&)"
#define G_READ_XML               "GCTAModelSpatialLookup::read(GXmlElement&)"
#define G_WRITE_XML             "GCTAModelSpatialLookup::write(GXmlElement&)"
#define G_FILL_1             "GCTAModelSpatialLookup::fill(GCTAObservation&)"
#define G_FILL_2               "GCTAModelSpatialLookup::fill(GObservations&)"
#define G_READ_TABLE              "GCTAModelSpatialLookup::read(GFitsTable&)"
#define G_LOAD                     "GCTAModelSpatialLookup::load(GFilename&)"
#define G_FILL_BUFFER "GCTAModelSpatialLookup::fill_buffer(GCTAObservation&,"\
                                                     " std::vector<double>&)"

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
 ***************************************************************************/
GCTAModelSpatialLookup::GCTAModelSpatialLookup(void) : GCTAModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Lookup table file constructor
 *
 * @param[in] filename Lookup table file name.
 *
 * Creates instance of a spatial lookup table model from a lookup table FITS
 * file.
 ***************************************************************************/
GCTAModelSpatialLookup::GCTAModelSpatialLookup(const GFilename& filename) :
                        GCTAModelSpatial()
{
    // Initialise members
    init_members();

    // Load lookup table
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response table constructor
 *
 * @param[in] table Response table.
 *
 * Creates instance of a spatial lookup table model from a response table.
 ***************************************************************************/
GCTAModelSpatialLookup::GCTAModelSpatialLookup(const GCTAResponseTable& table) :
                        GCTAModelSpatial()
{
    // Initialise members
    init_members();

    // Set lookup table
    this->table(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a spatial lookup table model from an XML element. See
 * the read() method for information about the expected structure of the XML
 * element.
 ***************************************************************************/
GCTAModelSpatialLookup::GCTAModelSpatialLookup(const GXmlElement& xml) :
                        GCTAModelSpatial()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Lookup table constructor
 *
 * @param[in] maxrad Maximum radius (deg).
 * @param[in] radbin Radial binsize (deg).
 * @param[in] ebds Energy boundaries.
 *
 * @exception GException::invalid_argument
 *            Non positive radial binsize specified.
 *
 * Creates empty instance of a spatial lookup table model using the maximum
 * lookup radius @p maxrad, the radial lookup binsize @p radbin, and some
 * energy boundaries @p ebds. Following construction, the instance may be
 * filled using the fill() methods.
 *
 * The method throws an exception if a non-positive radial bin size is
 * specified.
 ***************************************************************************/
GCTAModelSpatialLookup::GCTAModelSpatialLookup(const double&   maxrad,
                                               const double&   radbin,
                                               const GEbounds& ebds) :
                        GCTAModelSpatial()
{
    // Throw an exception if the radian binsize is not positive
    if (radbin <= 0.0) {
        std::string msg = "Non-positive radial bin size "+
                          gammalib::str(radbin)+" specified.";
        throw GException::invalid_argument(G_CONSTRUCTOR, msg);
    }

    // Initialise members
    init_members();

    // Setup energy axis vectors
    std::vector<double> eng_lo(ebds.size());
    std::vector<double> eng_hi(ebds.size());
    for (int i = 0; i < ebds.size(); ++i) {
        eng_lo[i] = ebds.emin(i).TeV();
        eng_hi[i] = ebds.emax(i).TeV();
    }

    // Setup radial axis vector
    int nrad = maxrad/radbin;
    std::vector<double> rad_lo(nrad);
    std::vector<double> rad_hi(nrad);
    for (int i = 0; i < nrad; ++i) {
        rad_lo[i] = double(i)*radbin;
        rad_hi[i] = double(i+1)*radbin;
    }

    // Setup response table
    m_lookup.append_axis(eng_lo, eng_hi, "ENERG", "TeV");
    m_lookup.append_axis(rad_lo, rad_hi, "THETA", "deg");
    m_lookup.append_table("BKG", "1/sr");

    // Prepare table
    prepare_table();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Lookup table spatial model.
 ***************************************************************************/
GCTAModelSpatialLookup::GCTAModelSpatialLookup(const GCTAModelSpatialLookup& model) :
                        GCTAModelSpatial(model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelSpatialLookup::~GCTAModelSpatialLookup(void)
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
 * @param[in] model Lookup table spatial model.
 ***************************************************************************/
GCTAModelSpatialLookup& GCTAModelSpatialLookup::operator=(const GCTAModelSpatialLookup& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GCTAModelSpatial::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GCTAModelSpatialLookup::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAModelSpatial::free_members();

    // Initialise members
    this->GCTAModelSpatial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GCTAModelSpatialLookup* GCTAModelSpatialLookup::clone(void) const
{
    return new GCTAModelSpatialLookup(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] dir Event direction.
 * @param[in] energy Event energy.
 * @param[in] time Event time.
 * @param[in] gradients Compute gradients?
 * @return Function value
 *
 * Evaluates the lookup table model for a given event direction and energy.
 * The evaluation is done using bi-linear interpolation in the logarithm
 * of the energy and the offset angle radius. The lookup table value is
 * multiplied with a normalisation parameter. The method always returns a
 * non-negative result.
 *
 * If the @p gradients flag is true the method will also compute the partial
 * derivatives of the normalisation parameter.
 ***************************************************************************/
double GCTAModelSpatialLookup::eval(const GCTAInstDir& dir,
                                    const GEnergy&     energy,
                                    const GTime&       time,
                                    const bool&        gradients) const
{
    // Compute offset angle in radians
    double offset = dir.theta();

    // Determine function value by bilinear interpolation
    double value = ((m_lookup.elements() > 0) &&
                    (m_lookup.tables()   > 0))
                   ? m_lookup(0, energy.log10TeV(), offset) : 0.0;

    // Make sure that function does not become negative
    if (value < 0.0) {
        value = 0.0;
    }

    // Optionally compute partial derivatives
    if (gradients) {

        // Compute partial derivative of the normalisation parameter
        double g_norm = (m_norm.is_free()) ? value * m_norm.scale() : 0.0;

        // Set gradient
        m_norm.factor_gradient(g_norm);

    } // endif: computed partial derivatives

    // Return intensity times normalization factor
    return (value * m_norm.value());
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Read the lookup table model information from an XML element. The XML
 * element needs to be of the following format:
 *
 *     <spatialModel type="LookupTable" file="lookuptable.fits">
 *       <parameter name="Normalization" ../>
 *     </spatialModel>
 *
 * The @p file attribute provides the filename of the lookup table FITS file.
 * The filename may be either an absolute filename (starting with '/') or a
 * relative filename. If no access path is given, the file is expected to
 * reside in the same location as the XML file.
 ***************************************************************************/
void GCTAModelSpatialLookup::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Load lookup table
    load(gammalib::xml_file_expand(xml, xml.attribute("file")));

    // Get parameter pointers
    const GXmlElement* norm  = gammalib::xml_get_par(G_READ_XML, xml, m_norm.name());

    // Read parameters
    m_norm.read(*norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Spatial model is not of valid type.
 *
 * Write the lookup table model information into an XML element. The XML
 * element will be of the following format:
 *
 *     <spatialModel type="LookupTable" file="lookuptable.fits">
 *       <parameter name="Normalization" ../>
 *     </spatialModel>
 *
 * The @p file attribute provides the filename of the lookup table FITS file.
 * The filename may be either an absolute filename (starting with '/') or a
 * relative filename. If no access path is given, the file is expected to
 * reside in the same location as the XML file.
 ***************************************************************************/
void GCTAModelSpatialLookup::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        std::string msg = "Spatial model \""+xml.attribute("type")+
                          "\" is not of type \""+type()+"\".";
        throw GException::invalid_value(G_WRITE_XML, msg);
    }

    // Set lookup table file name
    xml.attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Get XML parameters
    GXmlElement* norm  = gammalib::xml_need_par(G_WRITE_XML, xml, m_norm.name());

    // Write parameters
    m_norm.write(*norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill lookup table with events from one CTA observation
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No lookup table has been defined.
 *
 * Fills lookup table with events from one CTA observation. The method needs
 * that a lookup table is defined before. In case that no lookup table is
 * defined the method throws an exception.
 ***************************************************************************/
void GCTAModelSpatialLookup::fill(const GCTAObservation& obs)
{
    // Throw an exception if response table is not initialised
    if (m_lookup.elements() == 0) {
        std::string msg = "Attempting to fill a yet undefined lookup table. "
                          "Please allocate a valid lookup table before "
                          "calling this method.";
        throw GException::invalid_value(G_FILL_1, msg);
    }
    if (m_lookup.tables() == 0) {
        std::string msg = "Attempting to fill non-existing lookup table. "
                          "While the lookup table axes were defined, no table "
                          "exists. Please allocate a valid lookup table before "
                          "calling this method.";
        throw GException::invalid_value(G_FILL_1, msg);
    }

    // Allocate buffer for events
    std::vector<double> counts(m_lookup.elements());

    // Fill counts into buffer
    fill_buffer(obs, counts);

    // Set lookup table from buffer
    set_from_buffer(counts);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill lookup table with CTA events in observation container
 *
 * @param[in] obs Observation container.
 *
 * @exception GException::invalid_value
 *            No lookup table has been defined.
 *
 * Fills lookup table with CTA events in observation container. The method
 * needs that a lookup table is defined before. In case that no lookup table
 * is defined the method throws an exception.
 ***************************************************************************/
void GCTAModelSpatialLookup::fill(const GObservations& obs)
{
    // Throw an exception if response table is not initialised
    if (m_lookup.elements() == 0) {
        std::string msg = "Attempting to fill a yet undefined lookup table. "
                          "Please allocate a valid lookup table before "
                          "calling this method.";
        throw GException::invalid_value(G_FILL_2, msg);
    }
    if (m_lookup.tables() == 0) {
        std::string msg = "Attempting to fill non-existing lookup table. "
                          "While the lookup table axes were defined, no table "
                          "exists. Please allocate a valid lookup table before "
                          "calling this method.";
        throw GException::invalid_value(G_FILL_2, msg);
    }

    // Allocate buffer for events
    std::vector<double> counts(m_lookup.elements());

    // Loop over all observations in container
    for (int i = 0; i < obs.size(); ++i) {

        // Get observation and continue only if it is a CTA observation
        const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(obs[i]);

        // Skip observation if it's not CTA
        if (cta == NULL) {
            continue;
        }

        // Skip observation if we have a binned observation
        if (cta->eventtype() == "CountsCube") {
            continue;
        }

        // Fill counts into buffer
        fill_buffer(*cta, counts);

    }

    // Set lookup table from buffer
    set_from_buffer(counts);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Assign lookup table
 *
 * @param[in] table Lookup table.
 *
 * Assigns lookup table.
 ***************************************************************************/
void GCTAModelSpatialLookup::table(const GCTAResponseTable& table)
{
    // Assign lookup table
    m_lookup = table;
    
    // Prepare table
    prepare_table();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read lookup table from FITS table
 *
 * @param[in] table FITS table.
 *
 * @exception GException::invalid_value
 *            Response table is not two-dimensional.
 *
 * Reads the lookup table form the FITS @p table. The following column names
 * are mandatory:
 *
 *     ENERG_LO - Energy lower bin boundaries
 *     ENERG_HI - Energy upper bin boundaries
 *     THETA_LO - Offset angle lower bin boundaries
 *     THETA_HI - Offset angle upper bin boundaries
 *
 * After reading, the method assures that the lookup table is properly
 * normalised, i.e. that for each energy vector the radial profile has a
 * maximum value of one.
 ***************************************************************************/
void GCTAModelSpatialLookup::read(const GFitsTable& table)
{
    // Clear lookup table
    m_lookup.clear();

    // Read lookup table
    m_lookup.read(table);

    // Throw an exception if the table is not two-dimensional
    if (m_lookup.axes() != 2) {
        std::string msg = "Expected two-dimensional lookup table but found "+
                          gammalib::str(m_lookup.axes())+" dimensions. Please "
                          "specify a two-dimensional lookup table.";
        throw GException::invalid_value(G_READ_TABLE, msg);
    }

    // Prepare table
    prepare_table();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write lookup table into FITS binary table
 *
 * @param[in] table FITS binary table.
 *
 * Writes lookup table into a FITS binary @p table.
 ***************************************************************************/
void GCTAModelSpatialLookup::write(GFitsBinTable& table) const
{
    // Write lookup table
    m_lookup.write(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load lookup table
 *
 * @param[in] filename Load lookup table file name.
 *
 * Loads lookup table from FITS file.
 ***************************************************************************/
void GCTAModelSpatialLookup::load(const GFilename& filename)
{
    // Clear sigma spectrum
    clear();

    // Continue only if filename is not empty
    if (!filename.is_empty()) {

        // Open FITS file
        GFits fits(filename);

        // Open FITS table
        const GFitsTable* table = NULL;
        if (filename.has_extname()) {
            table = fits.table(filename.extname());
        }
        else if (filename.has_extno()) {
            table = fits.table(filename.extno());
        }
        else {
            table = fits.table(1);
        }

        // If no FITS table is found then throw an exception
        if (table == NULL) {
            std::string msg = "No lookup table found in file \""+
                              filename.url()+"\".";
            throw GException::file_error(G_LOAD, msg);
        }

        // Read lookup table
        read(*table);

        // Close FITS file
        fits.close();

        // Store filename
        m_filename = filename;

    } // endif: filename was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save lookup table into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file?
 *
 * Saves lookup table into a FITS file. If a file with the given @p filename
 * does not yet exist it will be created, otherwise the method opens the
 * existing file. The method will create a (or replace an existing)
 * effective area extension. The extension name can be specified as part
 * of the @p filename, or if no extension name is given, is assumed to be
 * `RADIAL BACKGROUND LOOKUP`.
 *
 * An existing file will only be modified if the @p clobber flag is set to
 * true.
 ***************************************************************************/
void GCTAModelSpatialLookup::save(const GFilename& filename, const bool& clobber) const
{
    // Get extension name
    std::string extname = filename.extname(gammalib::extname_cta_spatial_lookup);

    // Open or create FITS file (without extension name since the requested
    // extension may not yet exist in the file)
    GFits fits(filename.url(), true);

    // Remove extension if it exists already
    if (fits.contains(extname)) {
        fits.remove(extname);
    }

    // Create binary table
    GFitsBinTable table;

    // Write the lookup table
    write(table);

    // Set binary table extension name
    table.extname(extname);

    // Append table to FITS file
    fits.append(table);

    // Save to file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print lookup table information
 *
 * @param[in] chatter Chattiness.
 * @return String containing lookup table information.
 ***************************************************************************/
std::string GCTAModelSpatialLookup::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Initialise information
        int    nebins     = 0;
        int    nthetabins = 0;
        double emin       = 0.0;
        double emax       = 0.0;
        double omin       = 0.0;
        double omax       = 0.0;

        // Extract information if there are axes in the response table
        if (m_lookup.axes() > 0) {
            nebins     = m_lookup.axis_bins(m_inx_energy);
            nthetabins = m_lookup.axis_bins(m_inx_theta);
            emin       = m_lookup.axis_lo(m_inx_energy,0);
            emax       = m_lookup.axis_hi(m_inx_energy,nebins-1);
            omin       = m_lookup.axis_lo(m_inx_theta,0);
            omax       = m_lookup.axis_hi(m_inx_theta,nthetabins-1);
        }

        // Append header
        result.append("=== GCTAModelSpatialLookup ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(nebins));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(nthetabins));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");
        result.append("\n"+gammalib::parformat("Offset angle range"));
        result.append(gammalib::str(omin)+" - "+gammalib::str(omax)+" deg");

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
void GCTAModelSpatialLookup::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_lookup.clear();
    m_inx_energy = 0;
    m_inx_theta  = 1;

    // Initialise lookup table normalisation
    m_norm.clear();
    m_norm.name("Normalization");
    m_norm.value(1.0);
    m_norm.scale(1.0);
    m_norm.range(0.001, 1000.0);
    m_norm.gradient(0.0);
    m_norm.fix();
    m_norm.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Lookup table model.
 ***************************************************************************/
void GCTAModelSpatialLookup::copy_members(const GCTAModelSpatialLookup& model)
{
    // Copy members
    m_filename   = model.m_filename;
    m_lookup     = model.m_lookup;
    m_norm       = model.m_norm;
    m_inx_energy = model.m_inx_energy;
    m_inx_theta  = model.m_inx_theta;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelSpatialLookup::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Prepare lookup table indices
 *
 * Prepares a lookup table for usage.
 *
 * The method sets the data members m_inx_energy and m_inx_theta which
 * determine the axes order.
 *
 * The method furthermore sets the energy axis to logarithmic scale and the
 * offset angles to radians.
 *
 * Finally, the method assures that the maximum values for each radial
 * profile are one.
 ***************************************************************************/
void GCTAModelSpatialLookup::prepare_table(void)
{
    // Get mandatory indices (throw exception if not found)
    m_inx_energy = m_lookup.axis("ENERG");
    m_inx_theta  = m_lookup.axis("THETA");

    // Set energy axis to logarithmic scale
    m_lookup.axis_log10(m_inx_energy);

    // Set offset angle axis to radians
    m_lookup.axis_radians(m_inx_theta);

    // Normalise radial profiles
    normalise_table();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Normalise lookup table
 *
 * Normalises lookup table so that each radial profile has a maximum value
 * of one.
 ***************************************************************************/
void GCTAModelSpatialLookup::normalise_table(void)
{
    // Get axes dimensions
    int energy_size = m_lookup.axis_bins(m_inx_energy);
    int theta_size  = m_lookup.axis_bins(m_inx_theta);

    // Loop over energies
    for (int i_energy = 0; i_energy < energy_size; ++i_energy) {

        // Initialise maximum radial profile value
        double rad_max = 0.0;

        // Determine maximum radial profile value
        for (int i_theta = 0; i_theta < theta_size; ++i_theta) {
            int i = table_index(i_energy, i_theta);
            if (m_lookup(0,i) > rad_max) {
                rad_max = m_lookup(0,i);
            }
        }

        // If maximum value is positive then normalise profile
        if (rad_max > 0.0) {
            for (int i_theta = 0; i_theta < theta_size; ++i_theta) {
                int i          = table_index(i_energy, i_theta);
                m_lookup(0,i) /= rad_max;
            }
        }

    } // endfor: looped over all energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return index of lookup table element
 *
 * @param[in] ienergy Energy index.
 * @param[in] itheta Offset index.
 * @return Index of lookup table element.
 ***************************************************************************/
int GCTAModelSpatialLookup::table_index(const int& ienergy,
                                        const int& itheta) const
{
    // Set index vector
    int inx[2];
    inx[m_inx_energy] = ienergy;
    inx[m_inx_theta]  = itheta;

    // Compute index
    int index = inx[0] + inx[1] * m_lookup.axis_bins(0);

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Fill buffer from events in CTA observation
 *
 * @param[in] obs CTA observation.
 * @param[in,out] buffer Counts buffer.
 ***************************************************************************/
void GCTAModelSpatialLookup::fill_buffer(const GCTAObservation& obs,
                                         std::vector<double>&   buffer)
{
    // Make sure that the observation holds a CTA event list. If this is
    // not the case then throw an exception.
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs.events());
    if (events == NULL) {
        std::string msg = "CTA Observation does not contain an event "
                          "list. An event list is needed to fill the "
                          "lookup table.";
        throw GException::invalid_value(G_FILL_BUFFER, msg);
    }

    // Get axes dimensions
    int energy_size = m_lookup.axis_bins(m_inx_energy);
    int theta_size  = m_lookup.axis_bins(m_inx_theta);

    // Fill buffer
    for (int i = 0; i < events->size(); ++i) {

        // Get event
        const GCTAEventAtom* event = (*events)[i];

        // Determine offset angle in degrees
        GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
        double       theta = inst->theta() * gammalib::rad2deg;

        // Determine offset angle bin
        int i_theta = 0;
        for (; i_theta < theta_size; ++i_theta) {
            if ((theta >= m_lookup.axis_lo(m_inx_theta,i_theta)) &&
                (theta <= m_lookup.axis_hi(m_inx_theta,i_theta))) {
                break;
            }
        }

        // Skip event if no valid offset angle bin was found
        if (i_theta >= theta_size) {
            continue;
        }

        // Determine event energy in TeV
        double energy = event->energy().TeV();

        // Determine energy bin
        int i_energy = 0;
        for (; i_energy < energy_size; ++i_energy) {
            if ((energy >= m_lookup.axis_lo(m_inx_energy,i_energy)) &&
                (energy <= m_lookup.axis_hi(m_inx_energy,i_energy))) {
                break;
            }
        }

        // Skip event if no valid offset angle bin was found
        if (i_energy >= theta_size) {
            continue;
        }

        // Fill event into bugger
        int inx      = table_index(i_energy, i_theta);
        buffer[inx] += 1.0;

    } // endfor: looped over all events

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set lookup table from buffer
 *
 * @param[in] buffer Counts buffer.
 ***************************************************************************/
void GCTAModelSpatialLookup::set_from_buffer(const std::vector<double>& buffer)
{
    // Get axes dimensions
    int energy_size = m_lookup.axis_bins(m_inx_energy);
    int theta_size  = m_lookup.axis_bins(m_inx_theta);

    // Loop over offset angle
    for (int i_theta = 0; i_theta < theta_size; ++i_theta) {

        // Compute solid angle for offset angle bin
        double theta_min = m_lookup.axis_lo(m_inx_theta,i_theta) * gammalib::deg2rad;
        double theta_max = m_lookup.axis_hi(m_inx_theta,i_theta) * gammalib::deg2rad;
        double area      = gammalib::pi * (theta_max*theta_max - theta_min*theta_min);

        // Fill energy vector into lookup table
        for (int i_energy = 0; i_energy < energy_size; ++i_energy) {
            int inx         = table_index(i_energy, i_theta);
            m_lookup(0,inx) = buffer[inx] / area;
        }
        
    } // endfor: looped over offset angles
    
    // Normalise radial profiles
    normalise_table();

    // Return
    return;
}
