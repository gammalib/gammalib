/***************************************************************************
 *    GModelSpectralExpPlaw.cpp  -  Exponential cut off power law model    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralExpPlaw.cpp
 * @brief GModelSpectralExpPlaw class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralExpPlaw.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PAR                               "GModelSpectralExpPlaw::par(int)"
#define G_READ                    "GModelSpectralExpPlaw::read(GXmlElement&)"
#define G_WRITE                  "GModelSpectralExpPlaw::write(GXmlElement&)"

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
GModelSpectralExpPlaw::GModelSpectralExpPlaw(void) : GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] norm Power law normalization.
 * @param[in] index Power law index.
 * @param[in] ecut Cut off energy.
 *
 * Construct an exponential cut off power law from a normalization value,
 * a spectral index and a cut off energy.
 ***************************************************************************/
GModelSpectralExpPlaw::GModelSpectralExpPlaw(const double& norm,
                                             const double& index,
                                             const double& ecut) : GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Set parameters
    m_norm.value(norm);
    m_index.value(index);
    m_ecut.value(ecut);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element containing position information.
 *
 * Construct a spectral power law from a XML element.
 ***************************************************************************/
GModelSpectralExpPlaw::GModelSpectralExpPlaw(const GXmlElement& xml) : GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Model from which the instance should be built.
 ***************************************************************************/
GModelSpectralExpPlaw::GModelSpectralExpPlaw(const GModelSpectralExpPlaw& model) :
                                                        GModelSpectral(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpectralExpPlaw::~GModelSpectralExpPlaw(void)
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
 * @param[in] model Model which should be assigned.
 ***************************************************************************/
GModelSpectralExpPlaw& GModelSpectralExpPlaw::operator= (const GModelSpectralExpPlaw& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpectral::operator=(model);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GModelSpectralExpPlaw::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpectral::free_members();

    // Initialise members
    this->GModelSpectral::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelSpectralExpPlaw* GModelSpectralExpPlaw::clone(void) const
{
    return new GModelSpectralExpPlaw(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True energy of photon.
 *
 * The power law function is defined as
 * \f[I(E)=norm (E/pivot)^{index} \exp(-E/ecut)\f]
 * where
 * \f$norm\f$ is the normalization or prefactor,
 * \f$pivot\f$ is the pivot energy, and
 * \f$index\f$ is the spectral index, and
 * \f$ecut\f$ is the cut off energy.
 *
 * @todo For the moment the pivot and the cut off energies are fixed to
 * units of MeV. This may not be ideal and should eventually be improved
 * in the futur. Furthermore, the method expects that pivot!=0 and ecut!=0.
 * Otherwise Inf or NaN may result.
 ***************************************************************************/
double GModelSpectralExpPlaw::eval(const GEnergy& srcEng)
{
    // Compute function value
    double energy = srcEng.MeV();
    double e_norm = energy / pivot();
    double e_cut  = energy / ecut();
    double power  = pow(e_norm, index()) * exp(-e_cut);
    double value  = norm() * power;

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True energy of photon.
 *
 * The power law function is defined as
 * \f[I(E)=norm (E/pivot)^{index} \exp(-E/ecut)\f]
 * where
 * \f$norm=n_s n_v\f$ is the normalization or prefactor,
 * \f$index=i_s i_v\f$ is the spectral index,
 * \f$ecut=c_s c_v\f$ is the cut off energy, and
 * \f$pivot=p_s p_v\f$ is the pivot energy.
 * Note that each parameter is factorised into a scaling factor and a value
 * and that the method is expected to return the gradient with respect to
 * the parameter value (i.e. n_v, p_v, i_v, and c_v in this case).
 *
 * The partial derivatives of the parameter values are given by
 * \f[dI/dn_v=n_s  (E/pivot)^{index} \exp(-E/ecut)\f]
 * \f[dI/di_v=norm (E/pivot)^{index} \exp(-E/ecut) i_s \ln(E/pivot)\f]
 * \f[dI/dc_v=norm (E/pivot)^{index} \exp(-E/ecut) (E/ecut^2) c_s \f]
 * \f[dI/dp_v=norm (E/pivot)^{index} \exp(-E/ecut) (-index) / p_v\f]
 *
 * @todo For the moment the pivot and the cut off energies are fixed to
 * units of MeV. This may not be ideal and should eventually be improved
 * in the futur. Furthermore, the method expects that pivot!=0 and ecut!=0.
 * Otherwise Inf or NaN may result.
 ***************************************************************************/
double GModelSpectralExpPlaw::eval_gradients(const GEnergy& srcEng)
{
    // Compute function value
    double energy = srcEng.MeV();
    double e_norm = energy / pivot();
    double e_cut  = energy / ecut();
    double power  = pow(e_norm, index()) * exp(-e_cut);
    double value  = norm() * power;

    // Compute partial derivatives of the parameter values
    double g_norm  = (m_norm.isfree())  ? m_norm.scale() * power : 0.0;
    double g_index = (m_index.isfree()) ? value * m_index.scale() * log(e_norm) : 0.0;
    double g_ecut  = (m_ecut.isfree())  ? value * e_cut/m_ecut.value() : 0.0;
    double g_pivot = (m_pivot.isfree()) ? -value * index() / m_pivot.value() : 0.0;

    // Set gradients
    m_norm.gradient(g_norm);
    m_index.gradient(g_index);
    m_ecut.gradient(g_ecut);
    m_pivot.gradient(g_pivot);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Autoscale normalization
 *
 * Set the internal scaling of the normalization parameter so that the
 * parameter value equals 1.
 ***************************************************************************/
void GModelSpectralExpPlaw::autoscale(void)
{
    // Autoscale normalization to a value of 1.0
    if (m_norm.value() != 0.0) {

        // Get inverse scaling factor
        double invscale = 1.0 / m_norm.value();

        // Set values, error, min and max
        m_norm.value(m_norm.value() * invscale);
        m_norm.error(m_norm.error() * invscale);
        if (m_norm.hasmin())
            m_norm.min(m_norm.min() * invscale);
        if (m_norm.hasmax())
            m_norm.max(m_norm.max() * invscale);

        // Set scale
        m_norm.scale(1.0 / invscale);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing power law model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the spectral power law information from an XML element. The XML
 * element is required to have 4 parameters with names "Prefactor", "Index",
 * "Cutoff" and "Scale".
 ***************************************************************************/
void GModelSpectralExpPlaw::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Power law model requires exactly 4 parameters.");

    // Extract model parameters
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor") {
            m_norm.read(*par);
            npar[0]++;
        }

        // Handle index
        else if (par->attribute("name") == "Index") {
            m_index.read(*par);
            npar[1]++;
        }

        // Handle cutoff
        else if (par->attribute("name") == "Cutoff") {
            m_ecut.read(*par);
            npar[2]++;
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.read(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1)
        throw GException::model_invalid_parnames(G_READ, xml,
              "Require \"Prefactor\", \"Index\", \"Cutoff\" and \"Scale\""
              " parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "ExpCutoff"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the spectral power law information into an XML element. The XML
 * element has to be of type "ExpCutoff" and will have 4 parameter leafs
 * named "Prefactor", "Index", "Cutoff" and "Scale".
 ***************************************************************************/
void GModelSpectralExpPlaw::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "ExpCutoff");

    // Verify model type
    if (xml.attribute("type") != "ExpCutoff")
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"ExpCutoff\".");

    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"Prefactor\""));
        xml.append(new GXmlElement("parameter name=\"Index\""));
        xml.append(new GXmlElement("parameter name=\"Cutoff\""));
        xml.append(new GXmlElement("parameter name=\"Scale\""));
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Power law model requires exactly 4 parameters.");

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor") {
            npar[0]++;
            m_norm.write(*par);
        }

        // Handle index
        else if (par->attribute("name") == "Index") {
            npar[1]++;
            m_index.write(*par);
        }

        // Handle index
        else if (par->attribute("name") == "Cutoff") {
            npar[2]++;
            m_ecut.write(*par);
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.write(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1)
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Require \"Prefactor\", \"Index\", \"Cutoff\" and \"Scale\""
              " parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print powerlaw information
 ***************************************************************************/
std::string GModelSpectralExpPlaw::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralExpPlaw ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i) {
        result.append("\n"+parformat("Parameter "+str(i+1)));
        result.append(m_par[i]->print());
    }

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
void GModelSpectralExpPlaw::init_members(void)
{
    // Initialise parameters
    m_npars  = 4;
    m_par[0] = &m_norm;
    m_par[1] = &m_index;
    m_par[2] = &m_ecut;
    m_par[3] = &m_pivot;

    // Initialise powerlaw normalisation
    m_norm = GModelPar();
    m_norm.name("Prefactor");
    m_norm.unit("ph/cm2/s/MeV");
    m_norm.value(1.0);
    m_norm.scale(1.0);
    m_norm.free();

    // Initialise powerlaw index
    m_index = GModelPar();
    m_index.name("Index");
    m_index.value(-2.0);
    m_index.scale(1.0);
    m_index.free();

    // Initialise cut off energy
    m_ecut = GModelPar();
    m_ecut.name("Cutoff");
    m_ecut.unit("MeV");
    m_ecut.value(1000.0);
    m_ecut.scale(1.0);
    m_ecut.free();

    // Initialise pivot energy
    m_pivot = GModelPar();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.value(100.0);
    m_pivot.scale(1.0);
    m_pivot.fix();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpectralExpPlaw members which should be copied.
 ***************************************************************************/
void GModelSpectralExpPlaw::copy_members(const GModelSpectralExpPlaw& model)
{
    // Copy model parameters (we do not need to copy the rest since it is
    // static)
    m_norm  = model.m_norm;
    m_index = model.m_index;
    m_ecut  = model.m_ecut;
    m_pivot = model.m_pivot;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralExpPlaw::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] model Model.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelSpectralExpPlaw& model)
{
     // Write spectrum in output stream
    os << model.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] model Model.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GModelSpectralExpPlaw& model)
{
    // Write spectrum into logger
    log << model.print();

    // Return logger
    return log;
}
