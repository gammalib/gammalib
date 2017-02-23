/***************************************************************************
 *         GModelTemporalConst.cpp - Temporal constant model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2016 by Juergen Knoedlseder                         *
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
 * @file GModelTemporalConst.cpp
 * @brief Constant temporal model class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelTemporalConst.hpp"
#include "GModelTemporalRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelTemporalConst    g_temporal_const_seed;
const GModelTemporalRegistry g_temporal_const_registry(&g_temporal_const_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                      "GModelTemporalConst::read(GXmlElement&)"
#define G_WRITE                    "GModelTemporalConst::write(GXmlElement&)"

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
GModelTemporalConst::GModelTemporalConst(void) : GModelTemporal()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs constant temporal model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelTemporalConst::GModelTemporalConst(const GXmlElement& xml) :
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
 * @brief Value constructor
 *
 * @param[in] norm Normalization factor.
 *
 * Constructs constant temporal model by setting the normalization factor.
 ***************************************************************************/
GModelTemporalConst::GModelTemporalConst(const double& norm) :
                     GModelTemporal()
{
    // Initialise members
    init_members();

    // Set normalization factor
    m_norm.value(norm);

    // Return
    return;
}




/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Constant temporal model
 ***************************************************************************/
GModelTemporalConst::GModelTemporalConst(const GModelTemporalConst& model) : 
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
GModelTemporalConst::~GModelTemporalConst(void)
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
 * @param[in] model Constant temporal model.
 * @return Constant temporal model.
 ***************************************************************************/
GModelTemporalConst& GModelTemporalConst::operator=(const GModelTemporalConst& model)
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
 * @brief Clear constant temporal model
 ***************************************************************************/
void GModelTemporalConst::clear(void)
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
 * @brief Clone constant temporal model
 *
 * @return Pointer to deep copy of constant temporal model.
 ***************************************************************************/
GModelTemporalConst* GModelTemporalConst::clone(void) const
{
    // Clone constant temporal model
    return new GModelTemporalConst(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcTime True photon arrival time (not used).
 * @param[in] gradients Compute gradients?
 * @return Value to temporal model.
 *
 * Computes
 *
 * \f[
 *    S_{\rm t}(t) = {\tt m\_norm}
 * \f]
 *
 * where
 * \f${\tt m\_norm}\f$ is the normalization constant.
 *
 * If the @p gradients flag is true the method will also evaluate the partial
 * derivatives of the model with respect to the normalization parameter using
 *
 * \f[
 *    \frac{\delta S_{\rm t}(t)}{\delta {\tt m\_norm}} = 1
 * \f]
 ***************************************************************************/
double GModelTemporalConst::eval(const GTime& srcTime,
                                 const bool&  gradients) const
{
    // Compute function value
    double value = m_norm.value();

    // Optionally compute partial derivatives
    if (gradients) {

        // Compute partial derivatives of the parameter values
        double g_norm = (m_norm.is_free()) ? m_norm.scale() : 0.0;

        // Set factor gradient
        m_norm.factor_gradient(g_norm);

    } // endif: computed partial derivatives

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
 * This method returns a vector of random event times assuming a constant
 * event rate that is specified by the rate parameter.
 ***************************************************************************/
GTimes GModelTemporalConst::mc(const double& rate, const GTime&  tmin,
                               const GTime&  tmax, GRan& ran) const
{
    // Allocates empty vector of times
    GTimes times;

    // Compute event rate (in events per seconds)
    double lambda = rate * norm();

    // Initialise start and stop times in seconds
    double time  = tmin.secs();
    double tstop = tmax.secs();

    // Generate events until maximum event time is exceeded
    while (time <= tstop) {

        // Simulate next event time
        time += ran.exp(lambda);

        // Add time if it is not beyond the stop time
        if (time <= tstop) {
            GTime event;
            event.secs(time);
            times.append(event);
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
 * @exception GExpection::invalid_value
 *            Invalid XML format encountered.
 *
 * Writes the temporal information from an XML element having the format
 *
 *     <temporal type="Constant">
 *       <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="1"/>
 *     </temporal>
 ***************************************************************************/
void GModelTemporalConst::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        std::string msg = "Constant temporal model requires exactly 1 "
                          "<parameter> tag but "+
                          gammalib::str(xml.elements("parameter"))+
                          " <parameter> tags were found.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Get parameter element
    const GXmlElement* par = xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization") {
        m_norm.read(*par);
    }
    else {
        std::string msg = "Constant temporal model parameter "
                          "\"Normalization\" not found in XML element.";
        throw GException::invalid_value(G_READ, msg);
    }

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
 *     <temporal type="Constant">
 *       <parameter name="Normalization" scale="1" value="1" min="0.1" max="10" free="1"/>
 *     </temporal>
 ***************************************************************************/
void GModelTemporalConst::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", "Constant");
    }

    // Verify model type
    if (xml.attribute("type") != "Constant") {
        std::string msg = "Temporal model of type "+xml.attribute("type")+
                          " encountered while \"Constant\" was expected.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // If XML element has 0 nodes then append parameter node.
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Normalization\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        std::string msg = "Constant temporal model requires exactly 1 "
                          "<parameter> tag but "+
                          gammalib::str(xml.elements("parameter"))+
                          " <parameter> tags were found.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Get parameter element
    GXmlElement* par = xml.element("parameter", 0);

    // Set parameyter
    if (par->attribute("name") == "Normalization") {
        m_norm.write(*par);
    }
    else {
        std::string msg = "Constant temporal model parameter "
                          "\"Normalization\" not found in XML element.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print constant information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelTemporalConst::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelTemporalConst ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
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
void GModelTemporalConst::init_members(void)
{
    // Initialise normalisation parameter
    m_norm.clear();
    m_norm.name("Normalization");
    m_norm.unit("(relative value)");
    m_norm.value(1.0);
    m_norm.fix();
    m_norm.gradient(0.0);
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
 * @param[in] model Constant temporal model
 ***************************************************************************/
void GModelTemporalConst::copy_members(const GModelTemporalConst& model)
{
    // Copy members
    m_norm = model.m_norm;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelTemporalConst::free_members(void)
{
    // Return
    return;
}
