/***************************************************************************
 *               GCTAModelSkyCube.cpp - CTA sky cube model class           *
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
 * @file GCTAModelSkyCube.cpp
 * @brief CTA sky cube model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GModelRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"
#include "GObservation.hpp"
#include "GCTAModelSkyCube.hpp"
#include "GCTASupport.hpp"

/* __ Globals ____________________________________________________________ */
const GCTAModelSkyCube g_cta_inst_sky_seed;
const GModelRegistry   g_cta_inst_sky_registry(&g_cta_inst_sky_seed);

/* __ Method name definitions ____________________________________________ */
#define G_EVAL        "GCTAModelSkyCube::eval(GEvent&, GObservation&, bool&)"
#define G_NPRED    "GCTAModelSkyCube::npred(GEnergy&, GTime&, GObservation&)"
#define G_MC                     "GCTAModelSkyCube::mc(GObservation&, GRan&)"
#define G_READ_XML_SPATIAL "GCTAModelSkyCube::read_xml_spatial(GXmlElement&)"
#define G_WRITE_XML_SPATIAL            "GCTAModelSkyCube::write_xml_spatial("\
                                                              "GXmlElement&)"
#define G_XML_SPECTRAL         "GCTAModelSkyCube::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL         "GCTAModelSkyCube::xml_temporal(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_NPRED_CACHE

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_NPRED

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelSkyCube::GCTAModelSkyCube(void) : GModelData()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a CTA sky cube model from the information provided by an
 * XML element (see GCTAModelSkyCube::read method).
 ***************************************************************************/
GCTAModelSkyCube::GCTAModelSkyCube(const GXmlElement& xml) : GModelData(xml)
{
    // Initialise members
    init_members();

    // Read XML
    read(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from filename and spectral component
 *
 * @param[in] filename Cube file name.
 * @param[in] spectral Spectral model component.
 *
 * Constructs a CTA sky cube model from a cube @p filename and a @p spectral
 * model component. The temporal component is assumed to be constant.
 ***************************************************************************/
GCTAModelSkyCube::GCTAModelSkyCube(const GFilename&      filename,
                                   const GModelSpectral& spectral) :
                  GModelData()
{
    // Initialise members
    init_members();

    // Allocate temporal constant model
    GModelTemporalConst temporal;

    // Load filename
    load(filename);

    // Clone model components
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from model components
 *
 * @param[in] filename Cube file name.
 * @param[in] spectral Spectral model component.
 * @param[in] temporal Temporal model component.
 *
 * Constructs a CTA sky cube model from a cube @p filename, a @p spectral
 * model component and a @p temporal component.
 ***************************************************************************/
GCTAModelSkyCube::GCTAModelSkyCube(const GFilename&      filename,
                                   const GModelSpectral& spectral,
                                   const GModelTemporal& temporal) :
                  GModelData()
{
    // Initialise members
    init_members();

    // Load filename
    load(filename);
    
    // Clone model components
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] sky CTA sky cube model.
 ***************************************************************************/
GCTAModelSkyCube::GCTAModelSkyCube(const GCTAModelSkyCube& sky) :
                  GModelData(sky)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(sky);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelSkyCube::~GCTAModelSkyCube(void)
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
 * @param[in] sky CTA sky cube model.
 * @return CTA sky cube model.
 ***************************************************************************/
GCTAModelSkyCube& GCTAModelSkyCube::operator=(const GCTAModelSkyCube& sky)
{
    // Execute only if object is not identical
    if (this != &sky) {

        // Copy base class members
        this->GModelData::operator=(sky);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(sky);

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
 * @brief Clear CTA sky cube model
 *
 * This method properly resets the CTA sky cube model to an initial state.
 ***************************************************************************/
void GCTAModelSkyCube::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelData::free_members();

    // Initialise members
    this->GModelData::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone CTA sky cube model
 *
 * @return Pointer to deep copy of CTA sky cube model.
 ***************************************************************************/
GCTAModelSkyCube* GCTAModelSkyCube::clone(void) const
{
    return new GCTAModelSkyCube(*this);
}


/***********************************************************************//**
 * @brief Load sky cube
 *
 * @param[in] filename Sky cube filename.
 *
 * Load the sky cube from the primary image extension.
 ***************************************************************************/
void GCTAModelSkyCube::load(const GFilename& filename)
{
    // Clear object
    clear();

    // Put into OpenMP criticial zone
    #pragma omp critical(GCTAModelSkyCube_load)
    {
        // Open FITS file
        GFits fits(filename);

        // Get HDUs
        const GFitsImage& hdu_skycube = *fits.image("Primary");
        const GFitsTable& hdu_ebounds = *fits.table(gammalib::extname_ebounds);

        // Read sky cube
        m_cube.read(hdu_skycube);

        // Read energy boundaries
        m_ebounds.read(hdu_ebounds);

        // Set energy node array
        for (int i = 0; i < m_ebounds.size(); ++i) {
            m_elogmeans.append(m_ebounds.elogmean(i).log10TeV());
        }

        // Close FITS file
        fits.close();

    } // end of OpenMP critical zone

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set spectral model component
 *
 * @param[in] spectral Pointer to spectral model component.
 *
 * Sets the spectral model component of the model.
 ***************************************************************************/
void GCTAModelSkyCube::spectral(const GModelSpectral* spectral)
{
    // Free spectral model component only if it differs from current
    // component. This prevents unintential deallocation of the argument
    if ((m_spectral != NULL) && (m_spectral != spectral)) {
        delete m_spectral;
    }

    // Clone spectral model component if it exists, otherwise set pointer
    // to NULL
    m_spectral = (spectral != NULL) ? spectral->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set temporal model component
 *
 * @param[in] temporal Pointer to temporal model component.
 *
 * Sets the temporal model component of the model.
 ***************************************************************************/
void GCTAModelSkyCube::temporal(const GModelTemporal* temporal)
{
    // Free temporal model component only if it differs from current
    // component. This prevents unintential deallocation of the argument
    if ((m_temporal != NULL) && (m_temporal != temporal)) {
        delete m_temporal;
    }

    // Clone temporal model component if it exists, otherwise set pointer
    // to NULL
    m_temporal = (temporal != NULL) ? temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return event rate in units of events MeV\f$^{-1}\f$
 *        s\f$^{-1}\f$ sr\f$^{-1}\f$
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 * @param[in] gradients Compute gradients?
 * @return Event rate (events MeV\f$^{-1}\f$ s\f$^{-1}\f$ sr\f$^{-1}\f$).
 *
 * Evaluates the sky cube model. The method returns a measured rate, defined
 * as the number of counts per MeV, steradian and ontime.
 *
 * If the @p gradients flag is true the method will also set the parameter
 * gradients of the model parameters.
 ***************************************************************************/
double GCTAModelSkyCube::eval(const GEvent&       event,
                              const GObservation& obs,
                              const bool&         gradients) const
{
    // Set indices and weighting factors for energy interpolation
    m_elogmeans.set_value(event.energy().log10TeV());
    int    inx_left  = m_elogmeans.inx_left();
    int    inx_right = m_elogmeans.inx_right();
    double wgt_left  = m_elogmeans.wgt_left();
    double wgt_right = m_elogmeans.wgt_right();

    // Get reference on CTA instrument direction from event
    const GCTAInstDir& dir = gammalib::cta_dir(G_EVAL, event);

    // Initialise non-normalised spatial component
    double spat_raw = 0.0;

    // If left weight is close to 1, use left node
    if (std::abs(wgt_left-1.0) < 1.0e-6) {
        spat_raw = m_cube(dir.dir(), inx_left);
    }

    // ... otherwise, if right weight is close to 1, use right node
    else if (std::abs(wgt_right-1.0) < 1.0e-6) {
        spat_raw = m_cube(dir.dir(), inx_right);
    }

    // ... otherwise perform interpolation
    else {
        double spat_left  = m_cube(dir.dir(), inx_left);
        double spat_right = m_cube(dir.dir(), inx_right);
        if (spat_left > 0.0 && spat_right > 0.0) {
            spat_raw = std::exp(wgt_left  * std::log(spat_left) +
                                wgt_right * std::log(spat_right));
        }
    }

    // Multiply-in normalisation
    double spat = spat_raw * m_norm.value();

    // Make sure that spatial component does not become negative
    if (spat < 0.0) {
        spat = 0.0;
    }

    // Divide spatial model value by event size
    if (event.size() != 0.0) {
        spat /= event.size();
    }
    else {
        spat = 0.0;
    }

    // Evaluate spectral and temporal components
    double spec = (spectral() != NULL)
                  ? spectral()->eval(event.energy(), event.time(), gradients)
                  : 1.0;
    double temp = (temporal() != NULL)
                  ? temporal()->eval(event.time(), gradients) : 1.0;

    // Compute value
    double value = spat * spec * temp;

    // If gradients were requested then multiply factors to spectral and
    // temporal gradients
    if (gradients) {

        // Multiply factors to spatial gradient
        if (m_norm.is_free()) {
            double fact   = spec * temp;
            double g_norm = spat_raw * m_norm.scale();
            m_norm.factor_gradient(g_norm * fact);
        }

        // Multiply factors to spectral gradients
        if (spectral() != NULL) {
            double fact = spat * temp;
            if (fact != 1.0) {
                for (int i = 0; i < spectral()->size(); ++i)
                    (*spectral())[i].factor_gradient((*spectral())[i].factor_gradient() * fact);
            }
        }

        // Multiply factors to temporal gradients
        if (temporal() != NULL) {
            double fact = spat * spec;
            if (fact != 1.0) {
                for (int i = 0; i < temporal()->size(); ++i)
                    (*temporal())[i].factor_gradient((*temporal())[i].factor_gradient() * fact);
            }
        }

    } // endif: gradients were requested

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return spatially integrated event rate in units of
 *        events MeV\f$^{-1}\f$ s\f$^{-1}\f$
 *
 * @param[in] obsEng Measured event energy.
 * @param[in] obsTime Measured event time.
 * @param[in] obs Observation.
 * @return Spatially integrated event rate
 *         (events MeV\f$^{-1}\f$ s\f$^{-1}\f$)
 *
 * Spatially integrates the sky cube model for a given measured event energy
 * and event time.
 ***************************************************************************/
double GCTAModelSkyCube::npred(const GEnergy&      obsEng,
                               const GTime&        obsTime,
                               const GObservation& obs) const
{
    // Initialise result
    double npred     = 0.0;
    bool   has_npred = false;

    // Build unique identifier
    std::string id = obs.instrument() + ":" + obs.id();

    // Check if Npred value is already in cache
    #if defined(G_USE_NPRED_CACHE)
    if (!m_npred_names.empty()) {

        // Search for unique identifier, and if found, recover Npred value
        // and break
        for (int i = 0; i < m_npred_names.size(); ++i) {
            if (m_npred_names[i] == id && m_npred_energies[i] == obsEng) {
                npred     = m_npred_values[i];
                has_npred = true;
                #if defined(G_DEBUG_NPRED)
                std::cout << "GCTAModelSkyCube::npred:";
                std::cout << " cache=" << i;
                std::cout << " npred=" << npred << std::endl;
                #endif
                break;
            }
        }

    } // endif: there were values in the Npred cache
    #endif

    // Continue only if no Npred cache value has been found
    if (!has_npred) {

        // Evaluate only if model is valid
        if (valid_model()) {

            // Set indices and weighting factors for energy interpolation
            m_elogmeans.set_value(obsEng.log10TeV());
            int    inx_left  = m_elogmeans.inx_left();
            int    inx_right = m_elogmeans.inx_right();
            double wgt_left  = m_elogmeans.wgt_left();
            double wgt_right = m_elogmeans.wgt_right();

            // Loop over all map pixels
            for (int i = 0; i < m_cube.npix(); ++i) {

                // Get bin value
                double value = wgt_left  * m_cube(i, inx_left) +
                               wgt_right * m_cube(i, inx_right);

                // Sum bin contents
                npred += value * m_cube.solidangle(i);

            }

            // Store result in Npred cache
            #if defined(G_USE_NPRED_CACHE)
            m_npred_names.push_back(id);
            m_npred_energies.push_back(obsEng);
            m_npred_times.push_back(obsTime);
            m_npred_values.push_back(npred);
            #endif

            // Debug: Check for NaN
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
                std::string origin  = "GCTAModelSkyCube::npred";
                std::string message = " NaN/Inf encountered (npred=" +
                                      gammalib::str(npred) + ")";
                gammalib::warning(origin, message);
            }
            #endif

        } // endif: model was valid

    } // endif: Npred computation required

    // Multiply in normalisation value
    npred *= m_norm.value();

    // Multiply in spectral and temporal components
    npred *= spectral()->eval(obsEng, obsTime);
    npred *= temporal()->eval(obsTime);

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 * @return Pointer to list of simulated events (needs to be de-allocated by
 *         client)
 *
 * @exception GException::feature_not_implemented
 *            Class does not support the generation of event lists.
 *
 * The simulation of an event list from a sky cube model is not implemented,
 * hence the method will always throw an exception.
 ***************************************************************************/
GCTAEventList* GCTAModelSkyCube::mc(const GObservation& obs, GRan& ran) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_MC,
          "MC computation not implemented for binned analysis.");

    // Return NULL pointer
    return NULL;

}


/***********************************************************************//**
 * @brief Read CTA sky cube model from XML element
 *
 * @param[in] xml XML element.
 *
 * Set up CTA sky cube model from the information provided by an XML element.
 * The XML element is expected to have the following structure
 *
 *     <source name="CTA sky cube" type="CTASkyCube" instrument="CTA">
 *       <spatialModel type="ModelCube" file="model_cube.fits">
 *         <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="0"/>
 *       </spatialModel>
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *     </source>
 *
 * Optionally, a temporal model may be provided using the following
 * syntax
 *
 *     <source name="CTA sky cube" type="CTASkyCube" instrument="CTA">
 *       <spatialModel type="ModelCube" file="model_cube.fits">
 *         <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="0"/>
 *       </spatialModel>
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *       <temporal type="...">
 *         ...
 *       </temporal>
 *     </source>
 *
 * If no temporal component is found a constant model is assumed.
 ***************************************************************************/
void GCTAModelSkyCube::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Get pointers on spatial and spectral components
    const GXmlElement* spat = xml.element("spatialModel", 0);
    const GXmlElement* spec = xml.element("spectrum", 0);

    // Set spatial and spectral models
    read_xml_spatial(*spat);
    m_spectral = xml_spectral(*spec);

    // Handle optional temporal model
    if (xml.elements("temporal")) {
        const GXmlElement* temporal = xml.element("temporal", 0);
        m_temporal = xml_temporal(*temporal);
    }
    else {
        GModelTemporalConst constant;
        m_temporal = constant.clone();
    }

    // Read model attributes
    read_attributes(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA cube background model into XML element
 *
 * @param[in] xml XML element.
 *
 * Write CTA cube background model information into an XML element.
 * The XML element will have the following structure
 *
 *     <source name="CTA sky cube" type="CTASkyCube" instrument="CTA">
 *       <spatialModel type="ModelCube" file="model_cube.fits">
 *         <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="0"/>
 *       </spatialModel>
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *     </source>
 *
 * If the model contains a non-constant temporal model, the temporal
 * component will also be written following the syntax
 *
 *     <source name="CTA sky cube" type="CTASkyCube" instrument="CTA">
 *       <spatialModel type="ModelCube" file="model_cube.fits">
 *         <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="0"/>
 *       </spatialModel>
 *       <spectrum type="...">
 *         ...
 *       </spectrum>
 *       <temporal type="...">
 *         ...
 *       </temporal>
 *     </source>
 *
 * If no temporal component is found a constant model is assumed.
 ***************************************************************************/
void GCTAModelSkyCube::write(GXmlElement& xml) const
{
    // Initialise pointer on source
    GXmlElement* src = NULL;

    // Search corresponding source
    int n = xml.elements("source");
    for (int k = 0; k < n; ++k) {
        GXmlElement* element = xml.element("source", k);
        if (element->attribute("name") == name()) {
            src = element;
            break;
        }
    }

    // If we have a temporal model that is either not a constant, or a
    // constant with a normalization value that differs from 1.0 then
    // write the temporal component into the XML element. This logic
    // assures compatibility with the Fermi/LAT format as this format
    // does not handle temporal components.
    bool write_temporal = ((m_temporal != NULL) &&
                           (m_temporal->type() != "Constant" ||
                            (*m_temporal)[0].value() != 1.0));

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = xml.append("source");
        src->append(GXmlElement("spatialModel"));
        src->append(GXmlElement("spectrum"));
        if (write_temporal) {
            src->append(GXmlElement("temporal"));
        }
    }

    // Write spatial model
    GXmlElement* spat = src->element("spatialModel", 0);
    write_xml_spatial(*spat);

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Optionally write temporal model
    if (write_temporal) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = src->element("temporal", 0);
            temporal()->write(*temp);
        }
    }

    // Write model attributes
    write_attributes(*src);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA cube background model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing CTA cube background model information.
 ***************************************************************************/
std::string GCTAModelSkyCube::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelSkyCube ===");

        // Determine number of parameters per type
        int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
        int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;

        // Append attributes
        result.append("\n"+print_attributes());

        // Append model type
        result.append("\n"+gammalib::parformat("Model type"));
        if (n_spectral > 0) {
            result.append("\""+spectral()->type()+"\"");
            if (n_temporal > 0) {
                result.append(" * ");
            }
        }
        if (n_temporal > 0) {
            result.append("\""+temporal()->type()+"\"");
        }

        // Append parameters
        result.append("\n"+gammalib::parformat("Number of parameters") +
                      gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of spectral par's") +
                      gammalib::str(n_spectral));
        for (int i = 0; i < n_spectral; ++i) {
            result.append("\n"+(*spectral())[i].print());
        }
        result.append("\n"+gammalib::parformat("Number of temporal par's") +
                      gammalib::str(n_temporal));
        for (int i = 0; i < n_temporal; ++i) {
            result.append("\n"+(*temporal())[i].print());
        }

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
void GCTAModelSkyCube::init_members(void)
{
    // Initialise Value
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

    // Initialise members
    m_filename.clear();
    m_cube.clear();
    m_ebounds.clear();
    m_elogmeans.clear();
    m_spectral = NULL;
    m_temporal = NULL;

    // Initialise Npred cache
    m_npred_names.clear();
    m_npred_energies.clear();
    m_npred_times.clear();
    m_npred_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] sky CTA sky cube model.
 ***************************************************************************/
void GCTAModelSkyCube::copy_members(const GCTAModelSkyCube& sky)
{
    // Copy members
    m_filename  = sky.m_filename;
    m_cube      = sky.m_cube;
    m_ebounds   = sky.m_ebounds;
    m_elogmeans = sky.m_elogmeans;
    m_norm      = sky.m_norm;

    // Copy cache
    m_npred_names    = sky.m_npred_names;
    m_npred_energies = sky.m_npred_energies;
    m_npred_times    = sky.m_npred_times;
    m_npred_values   = sky.m_npred_values;

    // Clone spectral and temporal model components
    m_spectral = (sky.m_spectral != NULL) ? sky.m_spectral->clone() : NULL;
    m_temporal = (sky.m_temporal != NULL) ? sky.m_temporal->clone() : NULL;

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelSkyCube::free_members(void)
{
    // Free memory
    if (m_spectral != NULL) delete m_spectral;
    if (m_temporal != NULL) delete m_temporal;

    // Signal free pointers
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * Set pointers to all model parameters. The pointers are stored in a vector
 * that is member of the GModelData base class.
 ***************************************************************************/
void GCTAModelSkyCube::set_pointers(void)
{
    // Clear parameters
    m_pars.clear();

    // Push normalisation parameter on stack
    m_pars.push_back(&m_norm);

    // Determine the number of parameters
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;
    int n_pars     = n_spectral + n_temporal;

    // Continue only if there are parameters
    if (n_pars > 0) {

        // Gather spectral parameters
        for (int i = 0; i < n_spectral; ++i) {
            m_pars.push_back(&((*spectral())[i]));
        }

        // Gather temporal parameters
        for (int i = 0; i < n_temporal; ++i) {
            m_pars.push_back(&((*temporal())[i]));
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Verifies if model has all components
 *
 * Returns 'true' if models has a spectral and a temporal component.
 * Otherwise returns 'false'.
 ***************************************************************************/
bool GCTAModelSkyCube::valid_model(void) const
{
    // Set result
    bool result = ((spectral() != NULL) && (temporal() != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Read spatial model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Read sky cube information from the spatial XML element. The expected
 * format of the XML element is
 *
 *     <spatialModel type="ModelCube" file="model_cube.fits">
 *       <parameter name="Normalization" ../>
 *     </spatialModel>
 ***************************************************************************/
void GCTAModelSkyCube::read_xml_spatial(const GXmlElement& xml)
{
    // Verify model type
    if (xml.attribute("type") != "ModelCube") {
        std::string msg = "Spatial model type \""+xml.attribute("type")+
                          "\" is not of type \"ModelCube\". Please verify "
                          "the XML format.";
        throw GException::invalid_value(G_READ_XML_SPATIAL, msg);
    }

    // Get filename
    GFilename filename = gammalib::xml_file_expand(xml, xml.attribute("file"));

    // Get XML parameter
    const GXmlElement* norm = gammalib::xml_get_par(G_READ_XML_SPATIAL, xml,
                                                    m_norm.name());

    // Read parameter
    m_norm.read(*norm);

    // Load sky cube
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write spatial model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Write sky cube information into the spatial XML element. The format of
 * the XML element is
 *
 *     <spatialModel type="ModelCube" file="model_cube.fits">
 *       <parameter name="Normalization" ../>
 *     </spatialModel>
 ***************************************************************************/
void GCTAModelSkyCube::write_xml_spatial(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "ModelCube");
    }

    // Verify model type
    if (xml.attribute("type") != "ModelCube") {
        std::string msg = "Spatial model type \""+xml.attribute("type")+
                          "\" is not of type \"ModelCube\". Please verify "
                          "the XML format.";
        throw GException::invalid_value(G_WRITE_XML_SPATIAL, msg);
    }

    // Set filename
    xml.attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Get XML parameter
    GXmlElement* norm = gammalib::xml_need_par(G_WRITE_XML_SPATIAL, xml,
                                               m_norm.name());

    // Write parameter
    m_norm.write(*norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pointer to spectral model from XML element
 *
 * @param[in] spectral XML element.
 * @return Pointer to spectral model.
 *
 * Returns pointer to spectral model that is defined in an XML element.
 ***************************************************************************/
GModelSpectral* GCTAModelSkyCube::xml_spectral(const GXmlElement& spectral) const
{
    // Get spectral model
    GModelSpectralRegistry registry;
    GModelSpectral*        ptr = registry.alloc(spectral);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Return pointer to temporal model from XML element
 *
 * @param[in] temporal XML element.
 * @return Pointer to temporal model.
 *
 * Returns pointer to temporal model that is defined in an XML element.
 ***************************************************************************/
GModelTemporal* GCTAModelSkyCube::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(temporal);

    // Return pointer
    return ptr;
}
