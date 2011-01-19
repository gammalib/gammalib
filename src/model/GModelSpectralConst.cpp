/***************************************************************************
 *        GModelSpectralConst.cpp  -  Spectral constant model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralConst.cpp
 * @brief GModelSpectralConst class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralConst.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MC             "GModelSpectralConst::mc(GEnergy&, GEnergy&, GRan&)"
#define G_READ                      "GModelSpectralConst::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralConst::write(GXmlElement&)"

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
GModelSpectralConst::GModelSpectralConst(void) : GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 *
 * Construct an arbitrary spectral function from a XML element.
 ***************************************************************************/
GModelSpectralConst::GModelSpectralConst(const GXmlElement& xml) :
                     GModelSpectral()
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
 * @param[in] model Spectral component.
 ***************************************************************************/
GModelSpectralConst::GModelSpectralConst(const GModelSpectralConst& model) :
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
GModelSpectralConst::~GModelSpectralConst(void)
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
 * @param[in] model Spectral component.
 ***************************************************************************/
GModelSpectralConst& GModelSpectralConst::operator= (const GModelSpectralConst& model)
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
void GModelSpectralConst::clear(void)
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
GModelSpectralConst* GModelSpectralConst::clone(void) const
{
    return new GModelSpectralConst(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True energy of photon.
 *
 * The spectral model is defined as
 * \f[I(E)=norm\f]
 * where
 * \f$norm\f$ is the normalization of the function.
 ***************************************************************************/
double GModelSpectralConst::eval(const GEnergy& srcEng)
{
    // Compute function value
    double value = norm();

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True energy of photon.
 *
 * The spectral model is defined as
 * \f[I(E)=norm\f]
 * where
 * \f$norm=n_s n_v\f$ is the normalization of the function.
 * Note that the normalization is factorised into a scaling factor and a
 * value and that the method is expected to return the gradient with respect
 * to the parameter value \f$n_v\f$.
 *
 * The partial derivative of the normalization value is given by
 * \f[dI/dn_v=n_s\f]
 ***************************************************************************/
double GModelSpectralConst::eval_gradients(const GEnergy& srcEng)
{
    // Compute function value
    double value = norm();

    // Compute partial derivatives of the parameter values
    double g_norm = (m_norm.isfree())  ? m_norm.scale() : 0.0;

    // Set gradients
    m_norm.gradient(g_norm);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model flux between [emin, emax] (units: ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * Computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} I(E) dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$I(E)\f$ is the spectral model (units: ph/cm2/s/MeV).
 ***************************************************************************/
double GModelSpectralConst::flux(const GEnergy& emin, const GEnergy& emax) const
{
    // Compute flux for a constant model
    double flux = norm() * (emax.MeV() - emin.MeV());

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns MC energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] ran Random number generator.
 *
 * @todo To be implemented
 ***************************************************************************/
GEnergy GModelSpectralConst::mc(const GEnergy& emin, const GEnergy& emax,
                                GRan& ran) const
{
    // Allocate energy
    GEnergy energy;

    // Dump warning that method is not yet implemented
    throw GException::feature_not_implemented(G_MC);

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing power law model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
 *
 * Read the function information from an XML element and load the nodes
 * from the associated file. The XML element is required to have a parameter
 * named "Normalization" or "Value".
 ***************************************************************************/
void GModelSpectralConst::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Spectral constant requires exactly 1 parameter.");

    // Get parameter element
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value")
        m_norm.read(*par);
    else
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Spectral constant requires either"
                          " \"Normalization\" or \"Value\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "ConstantValue"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the spectral constant information into an XML element. The XML
 * element has to be of type "ConstantValue" and will have 1 parameter leaf
 * named "Normalization" or "Value" (default).
 ***************************************************************************/
void GModelSpectralConst::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "ConstantValue");

    // Verify model type
    if (xml.attribute("type") != "ConstantValue")
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"ConstantValue\".");

    // If XML element has 0 nodes then append parameter node. The name
    // of the node is "Value" as this is the Fermi-LAT standard.
    // We thus assure that the XML files will be compatible with
    // Fermi-LAT.
    if (xml.elements() == 0)
        xml.append(new GXmlElement("parameter name=\"Value\""));

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Spectral constant requires exactly 1 parameter.");

    // Get parameter element
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Set parameyter
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value")
        m_norm.write(*par);
    else
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Spectral constant requires either"
                          " \"Normalization\" or \"Value\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print powerlaw information
 ***************************************************************************/
std::string GModelSpectralConst::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralConst ===\n");
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
void GModelSpectralConst::init_members(void)
{
    // Initialise parameters
    m_npars  = 1;
    m_par[0] = &m_norm;

    // Initialise powerlaw normalisation
    m_norm = GModelPar();
    m_norm.name("Value");
    m_norm.value(1.0);
    m_norm.scale(1.0);
    m_norm.range(0.0, 1000.0);
    m_norm.free();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral model component.
 ***************************************************************************/
void GModelSpectralConst::copy_members(const GModelSpectralConst& model)
{
    // Copy model parameters
    m_norm = model.m_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralConst::free_members(void)
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
std::ostream& operator<< (std::ostream& os, const GModelSpectralConst& model)
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
GLog& operator<< (GLog& log, const GModelSpectralConst& model)
{
    // Write spectrum into logger
    log << model.print();

    // Return logger
    return log;
}
