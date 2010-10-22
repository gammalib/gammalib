/***************************************************************************
 *        GModelSpectralPlaw.cpp  -  Spectral power law model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralPlaw.cpp
 * @brief GModelSpectralPlaw class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralPlaw.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PAR                                  "GModelSpectralPlaw::par(int)"
#define G_READ                   "GModelSpectralPlaw::read(GXmlElement& xml)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GModelSpectralPlaw::GModelSpectralPlaw(void) : GModelSpectral()
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
 ***************************************************************************/
GModelSpectralPlaw::GModelSpectralPlaw(const double& norm, const double& index) :
    GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Set parameters
    m_norm.value(norm);
    m_index.value(index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element containing position information.
 ***************************************************************************/
GModelSpectralPlaw::GModelSpectralPlaw(const GXmlElement& xml) : GModelSpectral()
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
GModelSpectralPlaw::GModelSpectralPlaw(const GModelSpectralPlaw& model) :
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
GModelSpectralPlaw::~GModelSpectralPlaw(void)
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
GModelSpectralPlaw& GModelSpectralPlaw::operator= (const GModelSpectralPlaw& model)
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
 * @brief Returns pointer to a model parameter
 *
 * @param[in] index Parameter index.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of valid range
 ***************************************************************************/
GModelPar* GModelSpectralPlaw::par(int index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_npars)
        throw GException::out_of_range(G_PAR, index, 0, m_npars-1);

    // Return parameter pointer
    return m_par[index];
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True energy of photon.
 *
 * @todo For the moment the pivot energy is fixed to units of MeV. This may
 * not be ideal and should eventually be improved in the futur.
 * Furthermore, the method expects that energy!=0. Otherwise Inf or NaN
 * may result.
 ***************************************************************************/
double GModelSpectralPlaw::eval(const GEnergy& srcEng)
{
    // Compute function value
    double e     = srcEng.MeV() / pivot();
    double p     = pow(e, m_index.value());
    double value = m_norm.value() * p;

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True energy of photon.
 *
 * @todo For the moment the pivot energy is fixed to units of MeV. This may
 * not be ideal and should eventually be improved in the futur.
 * Furthermore, the method expects that energy!=0. Otherwise Inf or NaN
 * may result.
 ***************************************************************************/
double GModelSpectralPlaw::eval_gradients(const GEnergy& srcEng)
{
    // Compute function value
    double e     = srcEng.MeV() / pivot();
    double p     = pow(e, m_index.value());
    double value = m_norm.value() * p;

    // Set gradient to 0 if all parameters are fixed ...
    if ( !m_norm.isfree() && !m_index.isfree() && !m_pivot.isfree() ) {
        m_norm.gradient(0.0);
        m_index.gradient(0.0);
        m_pivot.gradient(0.0);
    }

    // ... otherwise compute gradient
    else {

        // Compute normalization gradient
        if (m_norm.isfree())
            m_norm.gradient(p);
        else
            m_norm.gradient(0.0);

        // Compute index gradient
        if (m_index.isfree()) {
            double g = m_norm.value() * p * log(srcEng.MeV());
            m_index.gradient(g);
        }
        else
            m_index.gradient(0.0);

        // Compute pivot energy gradient
        if (m_pivot.isfree()) {
            double g = -m_norm.value() * m_index.value() * p / pivot();
            m_pivot.gradient(g);
        }
        else
            m_pivot.gradient(0.0);

    } // endelse: gradient computation required

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Autoscale normalization
 ***************************************************************************/
void GModelSpectralPlaw::autoscale(void)
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
 ***************************************************************************/
void GModelSpectralPlaw::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 3 parameters
    if (xml.elements("parameter") != 3)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Power law model requires exactly 3 parameters.");

    // Extract model parameters
    int npar[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor") {
            m_norm.value(todouble(par->attribute("value")));
            m_norm.scale(todouble(par->attribute("scale")));
            m_norm.min(todouble(par->attribute("min")));
            m_norm.max(todouble(par->attribute("max")));
            if (par->attribute("free") == "1" ||
                tolower(par->attribute("free")) == "true")
                m_norm.free();
            else
                m_norm.fix();
            npar[0]++;
        }

        // Handle index
        else if (par->attribute("name") == "Index") {
            m_index.value(todouble(par->attribute("value")));
            m_index.scale(todouble(par->attribute("scale")));
            m_index.min(todouble(par->attribute("min")));
            m_index.max(todouble(par->attribute("max")));
            if (par->attribute("free") == "1" ||
                tolower(par->attribute("free")) == "true")
                m_index.free();
            else
                m_index.fix();
            npar[1]++;
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.value(todouble(par->attribute("value")));
            m_pivot.scale(todouble(par->attribute("scale")));
            m_pivot.min(todouble(par->attribute("min")));
            m_pivot.max(todouble(par->attribute("max")));
            if (par->attribute("free") == "1" ||
                tolower(par->attribute("free")) == "true")
                m_pivot.free();
            else
                m_pivot.fix();
            npar[2]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1)
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Require \"Prefactor\", \"Index\" and \"Scale\".");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @todo Implement method
 ***************************************************************************/
void GModelSpectralPlaw::write(GXmlElement& xml) const
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralPlaw::init_members(void)
{
    // Initialise parameters
    m_npars  = 3;
    m_par[0] = &m_norm;
    m_par[1] = &m_index;
    m_par[2] = &m_pivot;

    // Initialise powerlaw normalisation
    m_norm = GModelPar();
    m_norm.name("Prefactor");
    m_norm.unit("ph/cm2/s/MeV");
    m_norm.value(1.0);

    // Initialise powerlaw index
    m_index = GModelPar();
    m_index.name("Index");
    m_index.value(-2.0);

    // Initialise pivot energy
    m_pivot = GModelPar();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.value(100.0);
    m_pivot.fix();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpectralPlaw members which should be copied.
 ***************************************************************************/
void GModelSpectralPlaw::copy_members(const GModelSpectralPlaw& model)
{
    // Copy model parameters (we do not need to copy the rest since it is
    // static)
    m_norm  = model.m_norm;
    m_index = model.m_index;
    m_pivot = model.m_pivot;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralPlaw::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GModelSpectralPlaw* GModelSpectralPlaw::clone(void) const
{
    return new GModelSpectralPlaw(*this);
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put model in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] model Model to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelSpectralPlaw& model)
{
    // Put observation in stream
    os << "=== GModelSpectralPlaw ===" << std::endl;
    os << " Number of parameters ......: " << model.m_npars << std::endl;
    for (int i = 0; i < model.m_npars; ++i) {
        if (i > 0)
            os << std::endl;
        os << " Parameter .................: " << *(model.m_par[i]);
    }

    // Return output stream
    return os;
}
