/***************************************************************************
 *       GCTAModelRadialPolynom.cpp - Radial Polynom CTA model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialPolynom.cpp
 * @brief Radial Polynom model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GIntegral.hpp"
#include "GCTAModelRadialPolynom.hpp"
#include "GCTAModelRadialRegistry.hpp"

/* __ Constants __________________________________________________________ */
const double g_cta_radial_polynom_offset_max = 4.0;

/* __ Globals ____________________________________________________________ */
const GCTAModelRadialPolynom  g_cta_radial_polynom_seed;
const GCTAModelRadialRegistry g_cta_radial_polynom_registry(&g_cta_radial_polynom_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                   "GCTAModelRadialPolynom::read(GXmlElement&)"
#define G_WRITE                 "GCTAModelRadialPolynom::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_DEBUG_MC                                     //!< Debug MC method

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelRadialPolynom::GCTAModelRadialPolynom(void) : GCTAModelRadial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] coeffs Vector of polynomial coefficients.
 ***************************************************************************/
GCTAModelRadialPolynom::GCTAModelRadialPolynom(const std::vector<double>& coeffs)
                                                          : GCTAModelRadial()
{
    // Initialise members
    init_members();

    // Assign coefficients
    for (int i = 0; i < coeffs.size(); ++i) {

        // Allocate parameter
        GModelPar par;

        // Set value
        par.value(coeffs[i]);

        // Set other attributes
        std::string name = "coeff" + gammalib::str(i);
        par.name(name);
        par.unit("");
        par.free();
        par.scale(1.0);
        par.gradient(0.0);
        par.has_grad(true);

        // Push coefficient on list
        m_coeffs.push_back(par);

    } // endfor: looped over coefficients

    // Update parameter mapping
    update_pars();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a radial Polynom model by extracting information
 * from an XML element. See GCTAModelRadialPolynom::read() for more
 * information about the expected structure of the XML element.
 ***************************************************************************/
GCTAModelRadialPolynom::GCTAModelRadialPolynom(const GXmlElement& xml)
                                                          : GCTAModelRadial()
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
 * @param[in] model Radial Polynom model.
 ***************************************************************************/
GCTAModelRadialPolynom::GCTAModelRadialPolynom(const GCTAModelRadialPolynom& model)
                                                     : GCTAModelRadial(model)
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
GCTAModelRadialPolynom::~GCTAModelRadialPolynom(void)
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
 * @param[in] model Radial Polynom model.
 ***************************************************************************/
GCTAModelRadialPolynom& GCTAModelRadialPolynom::operator=(const GCTAModelRadialPolynom& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GCTAModelRadial::operator=(model);

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
void GCTAModelRadialPolynom::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAModelRadial::free_members();

    // Initialise members
    this->GCTAModelRadial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GCTAModelRadialPolynom* GCTAModelRadialPolynom::clone(void) const
{
    return new GCTAModelRadialPolynom(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] offset Offset angle [degrees].
 *
 * Evaluates the radial polynomial model for a given offset. The model is
 * defined as
 * \f[f(\theta) = \sum_{i=0}^m c_i \theta^i\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$c_i\f$ are the polynomial coefficients.
 ***************************************************************************/
double GCTAModelRadialPolynom::eval(const double& offset) const
{
    // Initialise result
    double value = 0.0;

    // Determine polynomial degree
    int ncoeffs = m_coeffs.size();
    
    // Compute value and gradients (only if there are coefficients)
    if (ncoeffs > 0) {

        // Compute value
        value = m_coeffs[ncoeffs-1].value();
        for (int i = ncoeffs-2; i >= 0; i--) {
            value = value * offset + m_coeffs[i].value();
        }

    } // endif: there were coefficients

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTAModelRadialPolynom::eval";
        std::cout << "(offset=" << offset << "): NaN/Inf encountered";
        std::cout << " (value=" << value;
        for (int i = 0; i < ncoeffs; ++i) {
            std::cout << ", c" << i << "=" << m_coeffs[i].value();
        }
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] offset Offset angle [degrees].
 *
 * Evaluates the Polynom model for a given offset. The Polynom model is
 * defined as
 * Evaluates the Polynom model for a given offset. The Polynomian model is
 * defined as
 * \f[f(\theta) = \sum_{i=0}^m c_i \theta^i\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$c_i\f$ are the polynomial coefficients.
 *
 * The partial derivative of the Polynom are given by
 * \f[\frac{df}{d{c_i}_v} = {c_i}_s \theta^i\f]
 * where 
 * \f${c_i}_v\f$ is the value part, 
 * \f${c_i}_s\f$ is the scaling part, and 
 * \f$c_i = {c_i}_v {c_i}_s\f$. 
 *
 * Note that this method implements a function which is unity for
 * \f$\theta=0\f$.
 ***************************************************************************/
double GCTAModelRadialPolynom::eval_gradients(const double& offset) const
{
    // Initialise result
    double value = 0.0;

    // Determine polynomial degree
    int ncoeffs = m_coeffs.size();
    
    // Compute value and gradients (only if there are coefficients)
    if (ncoeffs > 0) {

        // Compute value
        value = m_coeffs[ncoeffs-1].value();
        for (int i = ncoeffs-2; i >= 0; i--) {
            value = value * offset + m_coeffs[i].value();
        }

        // Initialise theta^i for the first coefficient
        double offset_power = 1.0;

        // Compute gradients for all coefficients
        for (int i = 0; i < ncoeffs; ++i) {

            // Compute gradient
            double grad = offset_power * m_coeffs[i].scale();

            // Store gradient
            const_cast<GCTAModelRadialPolynom*>(this)->m_coeffs[i].factor_gradient(grad);

            // Increase offset power for next coefficient
            offset_power *= offset;

        } // endfor: looped over all coefficients

    } // endif: there were coefficients

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTAModelRadialPolynom::eval_gradients";
        std::cout << "(offset=" << offset << "): NaN/Inf encountered";
        std::cout << " (value=" << value;
        for (int i = 0; i < ncoeffs; ++i) {
            std::cout << ", c_" << i << "=" << m_coeffs[i].value();
        }
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in] dir Pointing direction.
 * @param[in] ran Random number generator.
 *
 * Draws an arbitrary CTA instrument position from
 * \f[f(\theta) = \sin(\theta) \left( \sum_{i=0}^m c_i \theta^i \right)\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$c_i\f$ are the polynomial coefficients,
 * using the rejection method.
 *
 * Note that the maximum offset angle is fixed by the constant
 * g_cta_radial_polynom_offset_max.
 * This needs eventually adjusting on real data. The main reason for the
 * tight limit is to avoid divergence at large offset angles.
 *
 * @todo This method actually assumes that the polynom is always < 1, which
 *       may not be necessarily the case. Ideally, the method should
 *       determine the maximum of the polynomial to throw events. This is
 *       a severe limitation and should rapidly be corrected.
 ***************************************************************************/
GCTAInstDir GCTAModelRadialPolynom::mc(const GCTAInstDir& dir, GRan& ran) const
{
    // Set constants
    const double u_max = sin(g_cta_radial_polynom_offset_max *
                             gammalib::deg2rad);

    // Debug option: initialise number if samples thrown for one value
    #if defined(G_DEBUG_MC)
    int n_samples = 0;
    #endif
    
    // Simulate offset from photon arrival direction until we're not
    // rejected anymore
    double value  = 0.0;
    double u      = 1.0;
    double offset = 0.0;
    do {
        // Throw offset angle
        offset = ran.uniform() * g_cta_radial_polynom_offset_max;
        
        // Compute function value at this offset angle
        value  = sin(offset * gammalib::deg2rad) * eval(offset);
        
        // Throw value for rejection method
        u = ran.uniform() * u_max;

        // Debug option: update number of samples
        #if defined(G_DEBUG_MC)
        n_samples++;
        #endif
        
    } while (u > value);

    // Debug option: print number if samples thrown for one value
    #if defined(G_DEBUG_MC)
    std::cout << "#=" << n_samples << " ";
    #endif

    // Simulate azimuth angle
    double phi = 360.0 * ran.uniform();

    // Rotate pointing direction by offset and azimuth angle
    GCTAInstDir mc_dir = dir;
    mc_dir.dir().rotate_deg(phi, offset);

    // Return MC direction
    return mc_dir;
}


/***********************************************************************//**
 * @brief Returns integral over radial model (in steradians)
 *
 * Computes
 * \f[\Omega = 2 \pi \int_0^{\pi} \sin \theta f(\theta) d\theta\f]
 * where
 * \f[f(\theta) = \sum_{i=0}^m c_i \theta^i\f]
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$c_i\f$ are the polynomial coefficients.
 *
 * The integration is performed numerically, and the upper integration bound
 * \f$\pi\f$
 * is fixed to < g_cta_radial_polynom_offset_max.
 * This needs eventually adjusting on real data. The main reason for the
 * tight limit is to avoid divergence at large offset angles.
 ***************************************************************************/
double GCTAModelRadialPolynom::omega(void) const
{
    // Set constants
    const double offset_max_rad = g_cta_radial_polynom_offset_max *
                                  gammalib::deg2rad;

    // Allocate integrand
    GCTAModelRadialPolynom::integrand integrand(this);

    // Allocate intergal
    GIntegral integral(&integrand);

    // Perform numerical integration
    double omega = integral.romb(0.0, offset_max_rad) * gammalib::twopi;

    // Return integral
    return omega;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the radial Polynom model information from an XML element. The XML
 * element is required to have at least one parameter. Parameters are named
 * "CoeffX", where "X" starts from "0".  
 *
 * @todo Implement a test of the coefficient boundaries.
 ***************************************************************************/
void GCTAModelRadialPolynom::read(const GXmlElement& xml)
{
    // Free space for coefficients
    m_coeffs.clear();

    // Get maximum number of coefficients from XML file
    int max_coeffs = xml.elements("parameter");

    // Throw an error if no parameters were found
    if (max_coeffs < 1) {
        std::string message = "Radial polynomial model requires at least"
                              " one coefficient.";
        throw GException::model_invalid_parnum(G_READ, xml, message);
    }

    // Verify XML file consistency and determine number of coefficients
    int ncoeffs = 0;
    std::vector<int> npar;
    npar.assign(max_coeffs, 0);
    for (int i = 0; i < max_coeffs; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Verify that parameter is indeed a coefficient
        if (par->attribute("name").compare(0,5,"Coeff") != 0) {
            throw GException::model_invalid_parnames(G_READ, xml, 
                  par->attribute("name"));
        }
        
        // Verify that parameter coefficient has index in valid range
        int index = -1;
        size_t nchars = par->attribute("name").length() - 5;
        if (nchars > 0) {
            index = gammalib::toint(par->attribute("name").substr(5, nchars));
        }
        else {
            throw GException::model_invalid_parnames(G_READ, xml,
                  "Radial polynomial \"Coeff\" parameter has no index.");
        }
        if (index < 0) {
            throw GException::model_invalid_parnames(G_READ, xml,
                  "Radial polynomial \"Coeff\" parameter index < 0.");
        }
        if (index >= max_coeffs) {
            std::string message = "There are "+gammalib::str(max_coeffs)+" parameters,"
                                  " hence polynomial coefficients are expected"
                                  " to run from 0 to "+gammalib::str(max_coeffs-1)+", yet"
                                  " a coefficient with index "+gammalib::str(index)+" was"
                                  " encountered.";
            throw GException::model_invalid_parnames(G_READ, xml, message);
        }

        // Increment parameter counter
        npar[index]++;
        
        // Update number of coefficients
        if (index+1 > ncoeffs) {
            ncoeffs = index + 1;
        }

    } // endfor: verified XML file consistency
    
    // Verify that the number of coefficients is between 1 and max_coeffs
    if (ncoeffs < 0 || ncoeffs > max_coeffs) {
        std::string message = "Radial polynomial model requires at between"
                              " 1 and "+gammalib::str(max_coeffs)+" parameters.";
        throw GException::model_invalid_parnum(G_READ, xml, message);
    }

    // Verify that all parameters were found
    for (int i = 0; i < ncoeffs; ++i) {
        if (npar[i] == 0) {
            std::string message = "Parameter \"Coeff"+gammalib::str(i)+"\" required,"
                                  " but not found in XML file.";
            throw GException::model_invalid_parnames(G_READ, xml, message);
        }
        else if (npar[i] > 1) {
            std::string message = "Multiple parameters \"Coeff"+gammalib::str(i)+"\""
                                  " found in XML file.";
            throw GException::model_invalid_parnames(G_READ, xml, message);
        }
    }

    // Finally get all coefficients in order
    for (int i = 0; i < ncoeffs; ++i) {

        // Set parameter name
        std::string name = "Coeff"+gammalib::str(i);

        // Get corresponding parameter element
        const GXmlElement* par = NULL;
        for (int k = 0; k < ncoeffs; ++k) {
            const GXmlElement* element = xml.element("parameter", k);
            if (element->attribute("name") == name) {
                par = element;
                break;
            }
        }

        // Make sure that we really have one (just a double check, this should
        // never fail)
        if (par == NULL) {
            std::string message = "Required parameter \""+name+"\" not found.";
            throw GException::model_invalid_parnames(G_READ, xml, message);
        }

        // Now read that parameter ...
        GModelPar coeff;
        coeff.read(*par);

        // ... set other attributes ...
        coeff.name(name);
        coeff.unit("");

        //TODO: Check parameter

        // ... and push it on the list
        m_coeffs.push_back(coeff);

    } // endfor: looped over all parameters

    // Update parameter mapping
    update_pars();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type 'Polynom'
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the polynomial radial model information into an XML element. The XML
 * element will have 1-10 parameter leafs named "CoeffX", where X runs from
 * 0 to 9.
 ***************************************************************************/
void GCTAModelRadialPolynom::write(GXmlElement& xml) const
{
    // Determine number of coefficients
    int ncoeffs = m_coeffs.size();

    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Radial polynomial model is not of type \""+type()+"\".");
    }

    // If XML element has 0 nodes then append nodes for coefficients
    if (xml.elements() == 0) {
        for (int i = 0; i < ncoeffs; ++i) {
            std::string name = "parameter name=\"Coeff"+gammalib::str(i)+"\"";
            xml.append(GXmlElement(name));
        }
    }

    // Verify that XML element has one parameter per coefficient
    if (xml.elements() != ncoeffs || xml.elements("parameter") != ncoeffs) {
        std::string message = "Radial polynomial model requires exactly " +
                              gammalib::str(ncoeffs) + " parameters.";
        throw GException::model_invalid_parnum(G_WRITE, xml, message);
    }

    // Set or update model parameter attributes
    std::vector<int> npar;
    npar.assign(ncoeffs, 0);
    for (int i = 0; i < ncoeffs; ++i) {

        // Set parameter name
        std::string name = "Coeff"+gammalib::str(i);

        // Get corresponding parameter element
        GXmlElement* par = NULL;
        for (int k = 0; k < ncoeffs; ++k) {
            GXmlElement* element = xml.element("parameter", k);
            if (element->attribute("name") == name) {
                par = element;
                break;
            }
        }

        // Make sure that we really have one (just a double check, this should
        // never fail)
        if (par == NULL) {
            std::string message = "Required parameter \""+name+"\" not found.";
            throw GException::model_invalid_parnames(G_WRITE, xml, message);
        }

        // Write parameter
        m_coeffs[i].write(*par);

        // Increment parameter counter
        npar[i]++;

    } // endfor: looped over all parameters

    // Verify that all parameters were present
    for (int i = 0; i < ncoeffs; ++i) {
        if (npar[i] != 1) {
            std::string message = "Parameter \"Coeff"+gammalib::str(i)+"\" required,"
                                  " but not found in XML file.";
            throw GException::model_invalid_parnames(G_WRITE, xml, message);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point source information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing point source information.
 ***************************************************************************/
std::string GCTAModelRadialPolynom::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelRadialPolynom ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of parameters") +
                      gammalib::str(size()));
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
void GCTAModelRadialPolynom::init_members(void)
{
    // Initialise polynomial coefficients
    m_coeffs.clear();

    // Update parameter mapping
    update_pars();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial Polynomian model.
 ***************************************************************************/
void GCTAModelRadialPolynom::copy_members(const GCTAModelRadialPolynom& model)
{
    // Copy members
    m_coeffs = model.m_coeffs;

    // Update parameter mapping
    update_pars();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelRadialPolynom::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update parameter mapping
 ***************************************************************************/
void GCTAModelRadialPolynom::update_pars(void)
{
    // Clear parameter pointer(s)
    m_pars.clear();

    // Get number of coefficients
    int ncoeffs = m_coeffs.size();

    // Set parameter pointers for all coefficients
    for (int i = 0; i < ncoeffs; ++i) {

        // Signal that we have gradients
        m_coeffs[i].has_grad(true);

        // Set pointer
        m_pars.push_back(&(m_coeffs[i]));

    } // endfor: looped over coefficients

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
