/***************************************************************************
 *     GModelTemporalLightCurve.cpp - Temporal light curve model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GModelTemporalLightCurve.cpp
 * @brief Light curve model class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelTemporalLightCurve.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GFits.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelTemporalLightCurve g_temporal_lightcurve_seed;
const GModelTemporalRegistry   g_temporal_lightcurve_registry(&g_temporal_lightcurve_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                 "GModelTemporalLightCurve::read(GXmlElement&)"
#define G_WRITE               "GModelTemporalLightCurve::write(GXmlElement&)"
#define G_LOAD_NODES       "GModelTemporalLightCurve::load_nodes(GFilename&)"

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
GModelTemporalLightCurve::GModelTemporalLightCurve(void) : GModelTemporal()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name of nodes.
 * @param[in] norm Normalization factor.
 *
 * Constructs light curve model from a list of nodes that is found in the
 * specified FITS file. See the load_nodes() method for more information
 * about the expected structure of the file.
 ***************************************************************************/
GModelTemporalLightCurve::GModelTemporalLightCurve(const GFilename& filename,
                                                   const double&    norm) :
                          GModelTemporal()
{
    // Initialise members
    init_members();

    // Load nodes
    load_nodes(filename);

    // Set normalization
    m_norm.value(norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs light curve model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelTemporalLightCurve::GModelTemporalLightCurve(const GXmlElement& xml) :
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
 * @param[in] model Light curve model.
 ***************************************************************************/
GModelTemporalLightCurve::GModelTemporalLightCurve(const GModelTemporalLightCurve& model) : 
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
GModelTemporalLightCurve::~GModelTemporalLightCurve(void)
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
 * @param[in] model Light curve model.
 * @return Light curve model.
 ***************************************************************************/
GModelTemporalLightCurve& GModelTemporalLightCurve::operator=(const GModelTemporalLightCurve& model)
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
 * @brief Clear light curve model
 ***************************************************************************/
void GModelTemporalLightCurve::clear(void)
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
 * @brief Clone light curve model
 *
 * @return Pointer to deep copy of light curve model.
 ***************************************************************************/
GModelTemporalLightCurve* GModelTemporalLightCurve::clone(void) const
{
    // Clone constant temporal model
    return new GModelTemporalLightCurve(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcTime True photon arrival time (not used).
 * @param[in] gradients Compute gradients?
 * @return Value of light curve model.
 *
 * Computes
 *
 * \f[
 *    S_{\rm t}(t) = r(t) \times {\tt m\_norm}
 * \f]
 *
 * where
 * \f$r(t)\f$ is the light curve, computed by linear interpolation between
 * the nodes in a FITS file, and \f${\tt m\_norm}\f$ is the normalization
 * constant.
 *
 * If the @p gradients flag is true the method will also evaluate the partial
 * derivatives of the model with respect to the normalization parameter using
 *
 * \f[
 *    \frac{\delta S_{\rm t}(t)}{\delta {\tt m\_norm}} = r(t)
 * \f]
 ***************************************************************************/
double GModelTemporalLightCurve::eval(const GTime& srcTime,
                                      const bool&  gradients) const
{
    // Initialise value
    double value = 0.0;

    // Optionally initialise gradient
    if (gradients) {
        m_norm.factor_gradient(0.0);
    }

    // Set value if the time is within the validity range
    if (srcTime >= m_tmin && srcTime <= m_tmax) {

        // Interpolate file function
        double func = m_nodes.interpolate(srcTime.convert(m_timeref), m_values);
        
        // Compute function value
        value  = m_norm.value() * func;
    
        // Optionally compute gradients
        if (gradients) {

            // Compute partial derivatives of the parameter values
            double g_norm  = (m_norm.is_free())  ? m_norm.scale() * func : 0.0;

            // Set gradients
            m_norm.factor_gradient(g_norm);

        } // endif: gradient computation was requested
    
    } // endif: time was within the validity range

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns vector of random event times
 *
 * @param[in] rate Mean event rate (events per second).
 * @param[in] tmin Minimum event time.
 * @param[in] tmax Maximum event time.
 * @param[in,out] ran Random number generator.
 *
 * This method returns a vector of random event times between @p tmin and
 * @p tmax assuming a light curve specified in a FITS file.
 ***************************************************************************/
GTimes GModelTemporalLightCurve::mc(const double& rate, const GTime&  tmin,
                                    const GTime&  tmax, GRan& ran) const
{
    // Allocates empty vector of times
    GTimes times;

    // Update Monte Carlo cache
    mc_update(tmin, tmax);

    // Continue only if effective duration is positive and cache is not empty
    if (m_mc_eff_duration > 0.0 && m_mc_cum.size() > 0) {

        // Compute mean number of times by multiplying the rate with the
        // effective duration and the normalization factor
        double lambda = rate * norm() * m_mc_eff_duration;

        // Compute number of times to be sampled
        int ntimes = int(ran.poisson(lambda)+0.5);

        // Loop over number of times
        for (int i = 0; i < ntimes; ++i) {

            // Determine in which bin we reside
            int inx = 0;
            if (m_mc_cum.size() > 1) {
                double u = ran.uniform();
                for (inx = m_mc_cum.size()-1; inx > 0; --inx) {
                    if (m_mc_cum[inx-1] <= u) {
                        break;
                    }
                }
            }

            // Get random time
            double seconds;
            if (m_mc_slope[inx] == 0.0) {
                seconds = m_mc_dt[inx] * ran.uniform() + m_mc_time[inx];
            }
            else {
                seconds = (std::sqrt(m_mc_offset[inx]*m_mc_offset[inx] +
                                     2.0 * m_mc_slope[inx] * ran.uniform()) -
                           m_mc_offset[inx]) / m_mc_slope[inx] + m_mc_time[inx];
            }

            // Append random time
            GTime time(seconds, m_timeref);
            times.append(time);

        } // endfor: looped over times

    } // endif: cache was valid

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
 *     <temporalModel type="LightCurve" file="..">
 *       <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="1"/>
 *     </temporalModel>
 ***************************************************************************/
void GModelTemporalLightCurve::read(const GXmlElement& xml)
{
    // Get parameter pointers
    const GXmlElement* norm  = gammalib::xml_get_par(G_READ, xml, m_norm.name());

    // Read parameters
    m_norm.read(*norm);

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
 * @exception GExpection::invalid_value
 *            Invalid XML format encountered.
 *
 * Writes the temporal information into an XML element in the format
 *
 *     <temporalModel type="LightCurve" file="..">
 *       <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="1"/>
 *     </temporalModel>
 ***************************************************************************/
void GModelTemporalLightCurve::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        std::string msg = "Temporal model of type "+xml.attribute("type")+
                          " encountered while \""+type()+"\" was expected.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Get XML parameters
    GXmlElement* norm  = gammalib::xml_need_par(G_WRITE, xml, m_norm.name());

    // Write parameters
    m_norm.write(*norm);

    // Set file attribute
    xml.attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print light curve information
 *
 * @param[in] chatter Chattiness.
 * @return String containing light curve information.
 ***************************************************************************/
std::string GModelTemporalLightCurve::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelTemporalLightCurve ===");

        // Append information
        result.append("\n"+gammalib::parformat("Light curve file"));
        result.append(m_filename.url());
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

        // Append reference MJD and time range
        result.append("\n"+gammalib::parformat("Reference MDJ"));
        result.append(gammalib::str(m_timeref.mjdref()));
        result.append("\n"+gammalib::parformat("Time range"));
        result.append(gammalib::str(m_tmin.convert(m_timeref)));
        result.append(" - ");
        result.append(gammalib::str(m_tmax.convert(m_timeref)));
        result.append(" "+m_timeref.timeunit());
        result.append(" ("+m_timeref.timesys()+")");

        // Append node information
        result.append("\n"+gammalib::parformat("Number of nodes"));
        result.append(gammalib::str(m_nodes.size()));
        for (int i = 0; i < m_nodes.size(); ++i) {
            result.append("\n");
            result.append(gammalib::parformat(" Time "+gammalib::str(m_nodes[i])+" s"));
            result.append(gammalib::str(m_values[i]));
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
void GModelTemporalLightCurve::init_members(void)
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

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Initialise members
    m_nodes.clear();
    m_values.clear();
    m_filename.clear();
    m_timeref.clear();
    m_tmin.clear();
    m_tmax.clear();
    
    // Initialise cache
    m_mc_tmin.clear();
    m_mc_tmax.clear();
    m_mc_eff_duration = 0.0;
    m_mc_cum.clear();
    m_mc_slope.clear();
    m_mc_offset.clear();
    m_mc_time.clear();
    m_mc_dt.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Light curve model.
 ***************************************************************************/
void GModelTemporalLightCurve::copy_members(const GModelTemporalLightCurve& model)
{
    // Copy members
    m_norm     = model.m_norm;
    m_nodes    = model.m_nodes;
    m_values   = model.m_values;
    m_filename = model.m_filename;
    m_timeref  = model.m_timeref;
    m_tmin     = model.m_tmin;
    m_tmax     = model.m_tmax;

    // Copy cache
    m_mc_tmin         = model.m_mc_tmin;
    m_mc_tmax         = model.m_mc_tmax;
    m_mc_eff_duration = model.m_mc_eff_duration;
    m_mc_cum          = model.m_mc_cum;
    m_mc_slope        = model.m_mc_slope;
    m_mc_offset       = model.m_mc_offset;
    m_mc_time         = model.m_mc_time;
    m_mc_dt           = model.m_mc_dt;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelTemporalLightCurve::free_members(void)
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
 *            File function FITS file is invalid
 *
 * Load the light curve nodes from a FITS file. The light curve nodes are
 * expected in the first extension of the FITS file, containing two mandatory
 * columns with names "TIME" and "NORM".
 ***************************************************************************/
void GModelTemporalLightCurve::load_nodes(const GFilename& filename)
{
    // Set maximum light curve normalization value, including a small margin
    const double max_norm = 1.0 + 1.0e-8;

    // Clear nodes and values
    m_nodes.clear();
    m_values.clear();

    // Set filename
    m_filename = filename;

    // Load FITS file
    GFits fits = GFits(filename);

    // Extract binary table (so far always load extension 1 as table)
    GFitsTable* table = fits.table(1);

    // Read time reference from binary table
    m_timeref.read(*table);

    // Extract columns
    GFitsTableCol* time_col = (*table)["TIME"];
    GFitsTableCol* norm_col = (*table)["NORM"];

    // Check that there are at least two nodes in table
    if (time_col->nrows() < 2) {
        std::string msg = "\"TIME\" column contains "+
                          gammalib::str(time_col->nrows())+" rows but at "
                          "least two rows are required. Please specify a valid "
                          "temporal file function.";
        throw GException::invalid_value(G_LOAD_NODES, msg);
    }

    // Check that both columns are consistent
    if (time_col->nrows() != norm_col->nrows()) {
        std::string msg = "\"TIME\" and \"NORM\" columns have inconsistent "
                          "number of rows ("+
                          gammalib::str(time_col->nrows())+", "+
                          gammalib::str(norm_col->nrows())+"). Please "
                          "specify a valid temporal file function.";
        throw GException::invalid_value(G_LOAD_NODES, msg);
    }

    // Set number of nodes
    int nodes = time_col->nrows();

    // Check that time values are in ascending order and that no node is
    // larger than 1
    double last_time = -1.0;
    for (int i = 0; i < nodes; ++i) {

        // Check if time has increased
        if (last_time >= 0.0 && time_col->real(i) <= last_time) {
            std::string msg = "Time "+gammalib::str(time_col->real(i))+
                              " in row "+gammalib::str(i+1)+" of \"TIME\" "
                              "column is equal to or smaller than preceeding "
                              "value. Please provide a light curve file with "
                              "monotonically increasing times.";
            throw GException::invalid_value(G_LOAD_NODES, msg);
        }

        // Check if value is smaller than maximum allowed normalisation
        if (norm_col->real(i) > max_norm) {
            std::string msg = "Value "+gammalib::str(norm_col->real(i))+" at "
                              "time "+gammalib::str(time_col->real(i))+" is "
                              "larger than 1. Please provide a light curve file "
                              "with normalizations not exceeding 1.";
            throw GException::invalid_value(G_LOAD_NODES, msg);
        }
        
    } // endfor: looped over nodes

    // Extract nodes
    for (int i = 0; i < nodes; ++i) {
        m_nodes.append(time_col->real(i));
        m_values.push_back(norm_col->real(i));
    }

    // Make sure that no node exceeds 1
    for (int i = 0; i < m_values.size(); ++i) {
        if (m_values[i] > 1.0) {
            m_values[i] = 1.0;
        }
    }

    // Set minimum and maximum times (assumes that times are ordered)
    m_tmin.set(time_col->real(0), m_timeref);
    m_tmax.set(time_col->real(time_col->nrows()-1), m_timeref);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set MC pre-computation cache
 *
 * @param[in] tmin Minimum time of time interval.
 * @param[in] tmax Maximum time of time interval.
 *
 * This method sets up an array of indices and the cumulative distribution
 * function needed for MC simulations.
 ***************************************************************************/
void GModelTemporalLightCurve::mc_update(const GTime& tmin,
                                   const GTime& tmax) const
{
    // Update the cache only if the time interval has changed
    if (tmin != m_mc_tmin || tmax != m_mc_tmax) {
    
        // Store new time interval
        m_mc_tmin = tmin;
        m_mc_tmax = tmax;
        
        // Initialise cache
        m_mc_eff_duration = 0.0;
        m_mc_cum.clear();
        m_mc_slope.clear();
        m_mc_offset.clear();
        m_mc_time.clear();
        m_mc_dt.clear();

        // Initialise duration
        double duration = 0.0;

        // Continue only if time interval overlaps with temporal file function
        // and if time interval is valid
        if (tmax > m_tmin && tmin < m_tmax && tmax > tmin) {

            // Loop over all intervals between the nodes
            for (int i = 0; i < m_nodes.size()-1; ++i) {

                // Get start and stop time of interval
                GTime tstart(m_nodes[i], m_timeref);
                GTime tstop(m_nodes[i+1], m_timeref);

                // Make sure that we only consider the interval that overlaps
                // with the requested time interval
                if (tmin > tstart) {
                    tstart = tmin;
                }
                if (tmax < tstop) {
                    tstop = tmax;
                }

                // If time interval is empty then skip this node
                if (tstart >= tstop) {
                    continue;
                }

                // Evaluate rate at start and stop time
                double rstart = eval(tstart);
                double rstop  = eval(tstop);

                // Compute mean normalization for this interval
                double norm = 0.5 * (rstart + rstop);

                // Compute integral for this interval
                double dt  = tstop - tstart;
                double cum = norm * dt;

                // Compute slope, offset and start time of interval
                double slope  = (rstop - rstart) / dt;
                double offset = rstart;
                double renorm = (0.5 * slope * dt + offset) * dt;
                slope  /= renorm;
                offset /= renorm;

                // Update effective duration and real duration
                m_mc_eff_duration += cum;
                duration          += dt;

                // Put mean normalization and integral on stack
                m_mc_cum.push_back(cum);
                m_mc_slope.push_back(slope);
                m_mc_offset.push_back(offset);
                m_mc_time.push_back(tstart.convert(m_timeref));
                m_mc_dt.push_back(dt);

            } // endfor: next interval

            // Build cumulative distribution
            for (int i = 1; i < m_mc_cum.size(); ++i) {
                m_mc_cum[i] += m_mc_cum[i-1];
            }
            double norm = m_mc_cum[m_mc_cum.size()-1];
            for (int i = 0; i < m_mc_cum.size(); ++i) {
                m_mc_cum[i] /= norm;
            }

        } // endif: time interval overlaps
        
    } // endif: Update was required

    // Return
    return;
}
