/***************************************************************************
 *     GCTAModelRaccDummy.cpp  -  Radial acceptance dummy model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAModelRaccDummy.cpp
 * @brief GCTAModelRaccDummy class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelRegistry.hpp"
#include "GModelSpatialGauss.hpp"
#include "GCTAException.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAModelRaccDummy.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelRaccDummy g_cta_racc_dummy_seed;
const GModelRegistry     g_cta_racc_dummy_registry(&g_cta_racc_dummy_seed);

/* __ Method name definitions ____________________________________________ */
#define G_MC                   "GCTAModelRaccDummy::mc(GObservation&, GRan&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DUMP_MC 1                                 //!< Dump MC information


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelRaccDummy::GCTAModelRaccDummy(void) : GModelData(),
                                               GModelFactorized()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 ***************************************************************************/
GCTAModelRaccDummy::GCTAModelRaccDummy(const GXmlElement& xml) :
                                       GModelData(xml),
                                       GModelFactorized()
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
 * @brief Copy constructor
 *
 * @param[in] model Radial acceptance model.
 ***************************************************************************/
GCTAModelRaccDummy::GCTAModelRaccDummy(const GCTAModelRaccDummy& model) :
                                       GModelData(model),
                                       GModelFactorized(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelRaccDummy::~GCTAModelRaccDummy(void)
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
 * @param[in] model Radial acceptance model.
 ***************************************************************************/
GCTAModelRaccDummy& GCTAModelRaccDummy::operator= (const GCTAModelRaccDummy& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelData::operator=(model);
        this->GModelFactorized::operator=(model);

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
void GCTAModelRaccDummy::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelFactorized::free_members();
    this->GModelData::free_members();
    this->GModel::free_members();

    // Initialise members
    this->GModel::init_members();
    this->GModelData::init_members();
    this->GModelFactorized::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GCTAModelRaccDummy* GCTAModelRaccDummy::clone(void) const
{
    return new GCTAModelRaccDummy(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation (not used).
 *
 * @todo Verify that CTA instrument direction pointer is valid
 ***************************************************************************/
double GCTAModelRaccDummy::eval(const GEvent& event, const GObservation& obs)
{
    // Get instrument direction
    const GInstDir*    inst_dir = &(event.dir());
    const GCTAInstDir* cta_dir  = dynamic_cast<const GCTAInstDir*>(inst_dir);

    // Evaluate function
    double value = 1.0;
    if (spatial()  != NULL) value *= spatial()->eval(cta_dir->skydir());
    if (spectral() != NULL) value *= spectral()->eval(event.energy());
    if (temporal() != NULL) value *= temporal()->eval(event.time());

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation (not used).
 *
 * @todo Verify that CTA instrument direction pointer is valid
 ***************************************************************************/
double GCTAModelRaccDummy::eval_gradients(const GEvent& event,
                                          const GObservation& obs)
{
    // Get instrument direction
    const GInstDir*    inst_dir = &(event.dir());
    const GCTAInstDir* cta_dir  = dynamic_cast<const GCTAInstDir*>(inst_dir);

    // Evaluate function and gradients
    double value = 1.0;
    if (spatial()  != NULL) value *= spatial()->eval_gradients(cta_dir->skydir());
    if (spectral() != NULL) value *= spectral()->eval_gradients(event.energy());
    if (temporal() != NULL) value *= temporal()->eval_gradients(event.time());

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return simulated list of events
 *
 * @param[in] obs Observation.
 * @param[in] ran Random number generator.
 *
 * @exception GCTAException::no_pointing
 *            No pointing found
 *
 * This method requires that the energy boundaries, the good time interval,
 * the ROI of the observation and the pointing has been set up previously.
 ***************************************************************************/
GCTAEventList* GCTAModelRaccDummy::mc(const GObservation& obs, GRan& ran)
{
    // Initialise new event list
    GCTAEventList* list = new GCTAEventList;

    // Continue only if model is valid)
    if (valid_model()) {

        // Extract CTA pointing direction
        GTime time; // not used
        GCTAPointing* pnt = dynamic_cast<GCTAPointing*>(obs.pointing(time));
        if (pnt == NULL)
            throw GCTAException::no_pointing(G_MC);

        // Continue only if the spatial model is a Gaussian
        GModelSpatialGauss* src = dynamic_cast<GModelSpatialGauss*>(spatial());
        if (src != NULL) {

            // Loop over all energy boundaries
            for (int ieng = 0; ieng < obs.ebounds().size(); ++ieng) {

                // Compute the on-axis background rate in model within the
                // energy boundaries from spectral component (units: cts/s/sr)
                double flux = spectral()->flux(obs.ebounds().emin(ieng),
                                               obs.ebounds().emax(ieng));

                // Get the on-axis model value
                double onaxis = spatial()->eval(pnt->dir());

                // Compute solid angle used for normalization
                double area = 1.0 / onaxis;

                // Derive expecting rate (units: cts/s)
                double rate = flux * area;

                // Debug option: dump rate
                #if G_DUMP_MC
                std::cout << "GCTAModelRaccDummy::mc(\"" << name() << "\": ";
                std::cout << "flux=" << flux << " cts/s/sr, ";
                std::cout << "area=" << area << " sr, ";
                std::cout << "rate=" << rate << " cts/s)" << std::endl;
                #endif

                // Loop over all good time intervals
                for (int itime = 0; itime < obs.gti().size(); ++itime) {

                    // Get event arrival times from temporal model
                    GTimes times = m_temporal->mc(rate, 
                                                  obs.gti().tstart(itime),
                                                  obs.gti().tstop(itime), ran);

                    // Reserve space for events
                    if (times.size() > 0)
                        list->reserve(times.size());

                    // Loop over events
                    for (int i = 0; i < times.size(); ++i) {

                        // Simulate offset from photon arrival direction
                        double theta = src->sigma() * ran.chisq2();
                        double phi   = 360.0 * ran.uniform();

                        // Rotate pointing direction by offset
                        GSkyDir sky_dir = pnt->dir();
                        sky_dir.rotate(phi, theta);

                        // Set measured event arrival direction
                        GCTAInstDir inst_dir;
                        inst_dir.skydir(sky_dir);

                        // Set event energy
                        GEnergy energy = spectral()->mc(obs.ebounds().emin(ieng),
                                                        obs.ebounds().emax(ieng),
                                                        ran);

                        // Allocate event
                        GCTAEventAtom event;

                        // Set event attributes
                        event.dir(inst_dir);
                        event.energy(energy);
                        event.time(times[i]);

                        // Append event to list
                        list->append(event);

                    } // endfor: looped over all events

                } // endfor: looped over all GTIs

            } // endfor: looped over all energy boundaries

        } // endif: model was Gaussian

    } // endif: model was valid

    // Return
    return list;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GCTAModelRaccDummy::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Read model
    xml_read(xml);

    // Store name and instruments
    this->name(xml.attribute("name"));
    this->instruments(xml.attribute("instrument"));

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GCTAModelRaccDummy::write(GXmlElement& xml) const
{
    // Write model
    xml_write(xml, name(), instruments());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 ***************************************************************************/
std::string GCTAModelRaccDummy::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAModelRaccDummy ===");

    // Append name and instrument
    result.append("\n"+print_name_instrument());

    // Append model
    result.append("\n"+print_model());

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
void GCTAModelRaccDummy::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GCTAModelRaccDummy::copy_members(const GCTAModelRaccDummy& model)
{
    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelRaccDummy::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * @todo Not yet implemented
 ***************************************************************************/
void GCTAModelRaccDummy::set_pointers(void)
{
    // Set pointers
    m_pars = set_par_pointers();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
