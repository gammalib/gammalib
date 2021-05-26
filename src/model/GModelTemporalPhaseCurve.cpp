/***************************************************************************
 *     GModelTemporalPhaseCurve.cpp - Temporal phase curve model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
 * @file GModelTemporalPhaseCurve.cpp
 * @brief Temporal phase curve model class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelTemporalPhaseCurve.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelTemporalPhaseCurve g_temporal_phase_seed;
const GModelTemporalRegistry   g_temporal_phase_registry(&g_temporal_phase_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                 "GModelTemporalPhaseCurve::read(GXmlElement&)"
#define G_WRITE               "GModelTemporalPhaseCurve::write(GXmlElement&)"
#define G_LOAD_NODES       "GModelTemporalPhaseCurve::load_nodes(GFilename&)"
#define G_NORMALISE_NODES       "GModelTemporalPhaseCurve::normalize_nodes()"

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
GModelTemporalPhaseCurve::GModelTemporalPhaseCurve(void) : GModelTemporal()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name of phase curve nodes.
 * @param[in] mjd Reference time.
 * @param[in] phase Phase at reference time.
 * @param[in] f0 Frequency at reference time (Hz).
 * @param[in] f1 First frequency derivative at reference time.
 * @param[in] f2 Second frequency derivative at reference time.
 * @param[in] norm Normalization factor.
 * @param[in] normalize Normalize phase curve?
 *
 * Constructs phase curve from a list of nodes that is found in the specified
 * FITS file, a reference time and phase information at the reference time.
 * See the load_nodes() method for more information about the expected
 * structure of the file.
 ***************************************************************************/
GModelTemporalPhaseCurve::GModelTemporalPhaseCurve(const GFilename& filename,
                                                   const GTime&     mjd,
                                                   const double&    phase,
                                                   const double&    f0,
                                                   const double&    f1,
                                                   const double&    f2,
                                                   const double&    norm,
                                                   const bool&      normalize) :
                          GModelTemporal()
{
    // Initialise members
    init_members();

    // Set parameters
    this->mjd(mjd);
    this->phase(phase);
    this->f0(f0);
    this->f1(f1);
    this->f2(f2);
    this->norm(norm);

    // Set normalization flag
    m_normalize = normalize;

    // Load nodes
    load_nodes(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs phase curve by extracting information from an XML element. See
 * the read() method for more information about the expected structure of the
 * XML element.
 ***************************************************************************/
GModelTemporalPhaseCurve::GModelTemporalPhaseCurve(const GXmlElement& xml) :
                          GModelTemporal()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Phase curve.
 ***************************************************************************/
GModelTemporalPhaseCurve::GModelTemporalPhaseCurve(const GModelTemporalPhaseCurve& model) : 
                          GModelTemporal(model)
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
GModelTemporalPhaseCurve::~GModelTemporalPhaseCurve(void)
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
 * @param[in] model Phase curve.
 * @return Phase curve.
 ***************************************************************************/
GModelTemporalPhaseCurve& GModelTemporalPhaseCurve::operator=(const GModelTemporalPhaseCurve& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelTemporal::operator=(model);

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
 * @brief Clear phase curve
 ***************************************************************************/
void GModelTemporalPhaseCurve::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelTemporal::free_members();

    // Initialise members
    this->GModelTemporal::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone phase curve
 *
 * @return Pointer to deep copy of phase curve.
 ***************************************************************************/
GModelTemporalPhaseCurve* GModelTemporalPhaseCurve::clone(void) const
{
    // Clone constant temporal model
    return new GModelTemporalPhaseCurve(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcTime True photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Value of phase curve.
 *
 * Computes
 *
 * \f[
 *    S_{\rm t}(t) = r(\Phi(t)) \times {\tt m\_norm}
 * \f]
 *
 * where
 *
 * \f$r(\Phi(t))\f$ is the phase dependent rate, defined by linear
 * interpolation between the nodes in a FITS file, \f$\Phi(t)\f$ is the time
 * dependent phase, and \f${\tt m\_norm}\f$ is a normalisation constant.
 ***************************************************************************/
double GModelTemporalPhaseCurve::eval(const GTime& srcTime,
                                      const bool&  gradients) const
{
    // Get phase
    double phase = this->phase(srcTime);

    // Interpolate phase curve nodes
    double func = m_nodes.interpolate(phase, m_values);

    // Compute function value
    double value  = m_norm.value() * func;

    // Optionally compute gradients
    if (gradients) {

        // Compute partial derivatives of the parameter values
        double g_norm  = (m_norm.is_free())  ? m_norm.scale() * func : 0.0;

        // Set gradients
        m_norm.factor_gradient(g_norm);

    } // endif: gradient computation was requested

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Return vector of random event times
 *
 * @param[in] rate Mean event rate (events per second).
 * @param[in] tmin Minimum event time.
 * @param[in] tmax Maximum event time.
 * @param[in,out] ran Random number generator.
 *
 * Returns a vector of random event times between @p tmin and @p tmax for
 * the phase curve.
 ***************************************************************************/
GTimes GModelTemporalPhaseCurve::mc(const double& rate, const GTime&  tmin,
                                    const GTime&  tmax, GRan& ran) const
{
    // Allocates empty vector of times
    GTimes times;

    // Compute event rate (in events per seconds)
    double lambda = rate * norm() * m_scale;

    // Initialise start and stop times in seconds
    double time  = tmin.secs();
    double tstop = tmax.secs();

    // Generate events until maximum event time is exceeded
    while (time <= tstop) {

        // Simulate next event time
        time += ran.exp(lambda);

        // Break if time is beyond the stop time
        if (time > tstop) {
            break;
        }

        // Set time object
        GTime srcTime;
        srcTime.secs(time);

        // Get value for current time
        double value = eval(srcTime);

        // Get uniform random number
        double uniform = ran.uniform() * m_scale;

        // If the random number is not larger than the phase curve value
        // then accept the time
        if (uniform <= value) {
            times.append(srcTime);
        }

    } // endwhile: loop until stop time is reached

    // Return vector of times
    return times;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the temporal information from an XML element. The XML element should
 * have the format
 *
 *     <temporalModel type="PhaseCurve" file="phase.fits">
 *       <parameter name="Normalization" scale="1" value="1"       min="0.1" max="10"   free="0"/>
 *       <parameter name="MJD"           scale="1" value="51544.5" min="0.0" max="1e10" free="0"/>
 *       <parameter name="Phase"         scale="1" value="0.0"     min="0.0" max="1e10" free="1"/>
 *       <parameter name="F0"            scale="1" value="1.0"     min="0.0" max="1e10" free="1"/>
 *       <parameter name="F1"            scale="1" value="0.1"     min="0.0" max="1e10" free="1"/>
 *       <parameter name="F2"            scale="1" value="0.01"    min="0.0" max="1e10" free="1"/>
 *     </temporalModel>
 ***************************************************************************/
void GModelTemporalPhaseCurve::read(const GXmlElement& xml)
{
    // Verify number of model parameters
    gammalib::xml_check_parnum(G_READ, xml, 6);

    // Get parameter pointers
    const GXmlElement* norm  = gammalib::xml_get_par(G_READ, xml, m_norm.name());
    const GXmlElement* mjd   = gammalib::xml_get_par(G_READ, xml, m_mjd.name());
    const GXmlElement* phase = gammalib::xml_get_par(G_READ, xml, m_phase.name());
    const GXmlElement* f0    = gammalib::xml_get_par(G_READ, xml, m_f0.name());
    const GXmlElement* f1    = gammalib::xml_get_par(G_READ, xml, m_f1.name());
    const GXmlElement* f2    = gammalib::xml_get_par(G_READ, xml, m_f2.name());

    // Read parameters
    m_norm.read(*norm);
    m_mjd.read(*mjd);
    m_phase.read(*phase);
    m_f0.read(*f0);
    m_f1.read(*f1);
    m_f2.read(*f2);

    // Get optional normalization attribute
    m_normalize     = true;
    m_has_normalize = false;
    if (xml.has_attribute("normalize")) {
        m_has_normalize = true;
        std::string arg = xml.attribute("normalize");
        if (arg == "0" || gammalib::tolower(arg) == "false") {
            m_normalize = false;
        }
    }

    // Load nodes from file
    load_nodes(gammalib::xml_file_expand(xml, xml.attribute("file")));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes the temporal information into an XML element. The XML element will
 * have the format
 *
 *     <temporalModel type="PhaseCurve" file="phase.fits">
 *       <parameter name="Normalization" scale="1" value="1"       min="0.1" max="10"   free="0"/>
 *       <parameter name="MJD"           scale="1" value="51544.5" min="0.0" max="1e10" free="0"/>
 *       <parameter name="Phase"         scale="1" value="0.0"     min="0.0" max="1e10" free="1"/>
 *       <parameter name="F0"            scale="1" value="1.0"     min="0.0" max="1e10" free="1"/>
 *       <parameter name="F1"            scale="1" value="0.1"     min="0.0" max="1e10" free="1"/>
 *       <parameter name="F2"            scale="1" value="0.01"    min="0.0" max="1e10" free="1"/>
 *     </temporalModel>
 ***************************************************************************/
void GModelTemporalPhaseCurve::write(GXmlElement& xml) const
{
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Get XML parameters
    GXmlElement* norm  = gammalib::xml_need_par(G_WRITE, xml, m_norm.name());
    GXmlElement* mjd   = gammalib::xml_need_par(G_WRITE, xml, m_mjd.name());
    GXmlElement* phase = gammalib::xml_need_par(G_WRITE, xml, m_phase.name());
    GXmlElement* f0    = gammalib::xml_need_par(G_WRITE, xml, m_f0.name());
    GXmlElement* f1    = gammalib::xml_need_par(G_WRITE, xml, m_f1.name());
    GXmlElement* f2    = gammalib::xml_need_par(G_WRITE, xml, m_f2.name());

    // Write parameters
    m_norm.write(*norm);
    m_mjd.write(*mjd);
    m_phase.write(*phase);
    m_f0.write(*f0);
    m_f1.write(*f1);
    m_f2.write(*f2);

    // Set file attribute
    xml.attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Set optional normalization attribute
    if (m_has_normalize || !m_normalize) {
        if (m_normalize) {
            xml.attribute("normalize", "1");
        }
        else {
            xml.attribute("normalize", "0");
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes the phase for a given time
 *
 * @param[in] time Time.
 *
 * Computes
 *
 * \f[
 *    \Phi(t) = \Phi_0 + f(t-t_0) + \frac{1}{2}\dot{f} (t-t_0)^2 +
 *                                  \frac{1}{6}\ddot{f} (t-t_0)^3
 * \f]
 *
 * where
 * \f$t_0\f$ is a reference time,
 * \f$\Phi_0\f$ is the phase at the reference time,
 * \f$f\f$ is the variation frequency at the reference time,
 * \f$\dot{f}\f$ is the first derivative of the variation frequency at the
 * reference time, and
 * \f$\dot{f}\f$ is the second derivative of the variation frequency at the
 * reference time.
 *
 * The phase \f$\Phi(t)\f$ is in the interval [0,1].
 ***************************************************************************/
double GModelTemporalPhaseCurve::phase(const GTime& time) const
{
    // Set constants
    const double c1 = 0.5;
    const double c2 = 1.0 / 6.0;

    // Compute time since reference time in seconds
    double t = time - mjd();

    // Computes phase
    double phase = this->phase() + t * (f0() + t * (c1 * f1() + c2 * f2() * t));

    // Put phase into interval [0,1]
    phase -= floor(phase);

    // Return phase
    return phase;
}


/***********************************************************************//**
 * @brief Evaluate phase curve value for a given phase
 *
 * @param[in] phase Phase.
 * @return Value of phase curve.
 *
 * Computes
 *
 * \f[
 *    S_{\rm t}(\Phi) = r(\Phi) \times {\tt m\_norm}
 * \f]
 *
 * where
 *
 * \f$r(\Phi)\f$ is the phase dependent rate, defined by linear interpolation
 * between the nodes in a FITS file, \f$\Phi\f$ is the phase, and
 * \f${\tt m\_norm}\f$ is a normalisation constant.
 ***************************************************************************/
double GModelTemporalPhaseCurve::value(const double& phase) const
{
    // Compute function value
    double value  = m_norm.value() * m_nodes.interpolate(phase, m_values);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Print phase curve information
 *
 * @param[in] chatter Chattiness.
 * @return String containing phase curve information.
 ***************************************************************************/
std::string GModelTemporalPhaseCurve::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelTemporalPhaseCurve ===");

        // Append information
        result.append("\n"+gammalib::parformat("Function file"));
        result.append(m_filename.url());
        if (normalize()) {
            result.append(" [normalized]");
        }
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }
        result.append("\n"+gammalib::parformat("Number of phase nodes"));
        result.append(gammalib::str(m_nodes.size()));

        // EXPLICIT: Append node values
        if (chatter >= EXPLICIT) {
            for (int i = 0; i < m_nodes.size(); ++i) {
                result.append("\n");
                result.append(gammalib::parformat(" Phase "+gammalib::str(m_nodes[i])));
                result.append(gammalib::str(m_values[i]));
            }
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
void GModelTemporalPhaseCurve::init_members(void)
{
    // Initialise normalisation parameter
    m_norm.clear();
    m_norm.name("Normalization");
    m_norm.unit("(relative value)");
    m_norm.scale(1.0);
    m_norm.value(1.0);
    m_norm.range(0.0,1000.0);
    m_norm.fix();
    m_norm.gradient(0.0);
    m_norm.has_grad(true);

    // Initialise reference MJD parameter
    m_mjd.clear();
    m_mjd.name("MJD");
    m_mjd.unit("days");
    m_mjd.scale(1.0);
    m_mjd.value(51544.5);
    m_mjd.range(1.0e-10,1.0e+10);
    m_mjd.fix();
    m_mjd.gradient(0.0);
    m_mjd.has_grad(false);

    // Initialise phase at reference MJD parameter
    m_phase.clear();
    m_phase.name("Phase");
    m_phase.unit("");
    m_phase.scale(1.0);
    m_phase.value(0.0);
    m_phase.range(-1.0,1.0);
    m_phase.free();
    m_phase.gradient(0.0);
    m_phase.has_grad(false);

    // Initialise frequency at reference MJD parameter
    m_f0.clear();
    m_f0.name("F0");
    m_f0.unit("Hz");
    m_f0.scale(1.0);
    m_f0.value(0.0);
    m_f0.range(0.0,1000.0);
    m_f0.free();
    m_f0.gradient(0.0);
    m_f0.has_grad(false);

    // Initialise first frequency derivative at reference MJD parameter
    m_f1.clear();
    m_f1.name("F1");
    m_f1.unit("s^-2");
    m_f1.scale(1.0);
    m_f1.value(0.0);
    m_f1.range(-1000.0,1000.0);
    m_f1.free();
    m_f1.gradient(0.0);
    m_f1.has_grad(false);

    // Initialise second frequency derivative at reference MJD parameter
    m_f2.clear();
    m_f2.name("F2");
    m_f2.unit("s^-3");
    m_f2.scale(1.0);
    m_f2.value(0.0);
    m_f2.range(-1000.0,1000.0);
    m_f2.free();
    m_f2.gradient(0.0);
    m_f2.has_grad(false);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_mjd);
    m_pars.push_back(&m_phase);
    m_pars.push_back(&m_f0);
    m_pars.push_back(&m_f1);
    m_pars.push_back(&m_f2);

    // Initialise members
    m_nodes.clear();
    m_values.clear();
    m_filename.clear();
    m_scale         = 1.0;
    m_normalize     = true;
    m_has_normalize = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Phase curve.
 ***************************************************************************/
void GModelTemporalPhaseCurve::copy_members(const GModelTemporalPhaseCurve& model)
{
    // Copy members
    m_norm          = model.m_norm;
    m_mjd           = model.m_mjd;
    m_phase         = model.m_phase;
    m_f0            = model.m_f0;
    m_f1            = model.m_f1;
    m_f2            = model.m_f2;
    m_nodes         = model.m_nodes;
    m_values        = model.m_values;
    m_filename      = model.m_filename;
    m_scale         = model.m_scale;
    m_normalize     = model.m_normalize;
    m_has_normalize = model.m_has_normalize;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_mjd);
    m_pars.push_back(&m_phase);
    m_pars.push_back(&m_f0);
    m_pars.push_back(&m_f1);
    m_pars.push_back(&m_f2);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelTemporalPhaseCurve::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load nodes from file
 *
 * @param[in] filename File name.
 *
 * @exception GException::invalid_value
 *            Phase curve FITS file is invalid
 *
 * Load the phase curve nodes from a FITS file. The information is expected
 * in the first extension of the FITS file, containing two mandatory columns
 * with names "PHASE" and "NORM".
 ***************************************************************************/
void GModelTemporalPhaseCurve::load_nodes(const GFilename& filename)
{
    // Set maximum phase curve normalization value, including a small margin
    const double max_norm = 1.0 + 1.0e-8;

    // Clear nodes and values
    m_nodes.clear();
    m_values.clear();

    // Set filename
    m_filename = filename;

    // Continue only if filename is not empty
    if (!filename.is_empty()) {

        // Load FITS file
        GFits fits = GFits(filename);

        // Extract binary table (so far always load extension 1 as table)
        GFitsTable* table = fits.table(1);

        // Extract columns
        GFitsTableCol* phase_col = (*table)["PHASE"];
        GFitsTableCol* norm_col  = (*table)["NORM"];

        // Check that there are at least two nodes in table
        if (phase_col->nrows() < 2) {
            std::string msg = "\"PHASE\" column contains "+
                              gammalib::str(phase_col->nrows())+" rows but at "
                              "least two rows are required. Please specify a valid "
                              "phase curve file.";
            throw GException::invalid_value(G_LOAD_NODES, msg);
        }

        // Check that both columns are consistent
        if (phase_col->nrows() != norm_col->nrows()) {
            std::string msg = "\"PHASE\" and \"NORM\" columns have inconsistent "
                              "number of rows ("+
                              gammalib::str(phase_col->nrows())+", "+
                              gammalib::str(norm_col->nrows())+"). Please "
                              "specify a valid phase curve file.";
            throw GException::invalid_value(G_LOAD_NODES, msg);
        }

        // Set number of nodes
        int nodes = phase_col->nrows();

        // Check that phase values are in ascending order and comprised between
        // 0 and 1, and that no node is larger than 1
        double last_phase = -1.0;
        for (int i = 0; i < nodes; ++i) {

            // Check if phase has increased
            if (last_phase >= 0.0 && phase_col->real(i) <= last_phase) {
                std::string msg = "Phase value "+gammalib::str(phase_col->real(i))+
                                  " in row "+gammalib::str(i+1)+" of \"PHASE\" "
                                  "column is equal to or smaller than preceeding "
                                  "value. Please provide a phase curve file with "
                                  "monotonically increasing phase values.";
                throw GException::invalid_value(G_LOAD_NODES, msg);
            }

            // Check if phase is within [0,1]
            if (phase_col->real(i) < 0.0 || phase_col->real(i) > 1.0) {
                std::string msg = "Phase value "+gammalib::str(phase_col->real(i))+
                                  " outside range [0,1]. Please provide a phase "
                                  "curve file with phase values comprised in the "
                                  "interval [0,1].";
                throw GException::invalid_value(G_LOAD_NODES, msg);
            }

            // Check if value is smaller than maximum allowed normalisation
            if (norm_col->real(i) > max_norm) {
                std::string msg = "Value "+gammalib::str(norm_col->real(i))+" at "
                                  "phase "+gammalib::str(phase_col->real(i))+" is "
                                  "larger than 1. Please provide a phase curve file "
                                  "with normalizations not exceeding 1.";
                throw GException::invalid_value(G_LOAD_NODES, msg);
            }

        } // endfor: looped over nodes

        // If first phase is larger than zero then add last node with phase-1 as
        // first node. This makes sure that a valid normalization exists for
        // phase=0.
        if (phase_col->real(0) > 0.0) {
            m_nodes.append(phase_col->real(nodes-1)-1.0);
            m_values.push_back(norm_col->real(nodes-1));
        }

        // Extract nodes
        for (int i = 0; i < nodes; ++i) {
            m_nodes.append(phase_col->real(i));
            m_values.push_back(norm_col->real(i));
        }

        // If last node is smaller than one then add first node with phase+1 as
        // last node. This makes sure that a valid normalization exists for
        // phase=1.
        if (phase_col->real(nodes-1) < 1.0) {
            m_nodes.append(phase_col->real(0)+1.0);
            m_values.push_back(norm_col->real(0));
        }

        // Make sure that no node exceeds 1
        for (int i = 0; i < m_values.size(); ++i) {
            if (m_values[i] > 1.0) {
                m_values[i] = 1.0;
            }
        }

        // Close FITS file
        fits.close();

        // Optionally normalize nodes
        if (m_normalize) {
            normalize_nodes();
        }

    } // endif: filename was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Normalise nodes
 *
 * @exception GException::invalid_value
 *            Phase curve FITS file is invalid
 *
 * Normalise the node values so that the integral over the phase interval
 * [0,1] gives an average normalisation of 1.
 ***************************************************************************/
void GModelTemporalPhaseCurve::normalize_nodes(void)
{
    // Initialise node integral
    double sum = 0.0;

    // Loop over all nodes-1
    for (int i = 0; i < m_nodes.size()-1; ++i) {

        // Compute phase values
        double a = (m_nodes[i]   < 0.0) ? 0.0 : m_nodes[i];
        double b = (m_nodes[i+1] > 1.0) ? 1.0 : m_nodes[i+1];

        // Compute normalisation at phase values
        double fa = m_nodes.interpolate(a, m_values);
        double fb = m_nodes.interpolate(b, m_values);

        // Compute intergal for this interval using the Trapezoid rule
        sum += 0.5 * (b - a) * (fa + fb);

    }

    // Throw an exception if the integral is not positive
    if (sum <= 0.0) {
        std::string msg = "Integral over phase curve is not positive ("+
                          gammalib::str(sum)+"). Please provide a valid phase "
                          "curve file.";
        throw GException::invalid_value(G_NORMALISE_NODES, msg);
    }

    // Set scale factor
    m_scale = 1.0 / sum;

    // Normalise the node values
    for (int i = 0; i < m_values.size(); ++i) {
        m_values[i] *= m_scale;
    }

    // Return
    return;
}
