/***************************************************************************
 *            GException_model.cpp  -  Model exception handlers            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GException_model.cpp
 * @brief Implement exceptions for the model module
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GModelRegistry.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Invalid model type
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Model type that has been encountered.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_invalid::model_invalid(std::string origin,
                                         std::string type,
                                         std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid model type \""+type+"\" encountered. " + message;

    // Add list of valid models
    GModelRegistry registry;
    if (registry.size() > 0) {
        m_message += "The following models are registered: ";
        for (int i = 0; i < registry.size(); ++i) {
            if (i > 0)
                m_message += ", ";
            m_message += "\"" + registry.name(i) + "\"";
        }
        m_message += ".";
    }
    else
        m_message += "No models are registered.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid spatial model type
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Spatial model type that has been encountered.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_invalid_spatial::model_invalid_spatial(std::string origin,
                                                         std::string type,
                                                         std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid spatial model type \""+type+"\" encountered. " +
                message;

    // Add list of valid spatial models
    GModelSpatialRegistry registry;
    if (registry.size() > 0) {
        m_message += "The following models are registered: ";
        for (int i = 0; i < registry.size(); ++i) {
            if (i > 0)
                m_message += ", ";
            m_message += "\"" + registry.name(i) + "\"";
        }
        m_message += ".";
    }
    else
        m_message += "No models are registered.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid spectral model type
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Spectral model type that has been encountered.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_invalid_spectral::model_invalid_spectral(std::string origin,
                                                           std::string type,
                                                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid spectral model type \""+type+"\" encountered. " +
                message;

    // Add list of valid spectral models
    GModelSpectralRegistry registry;
    if (registry.size() > 0) {
        m_message += "The following models are registered: ";
        for (int i = 0; i < registry.size(); ++i) {
            if (i > 0)
                m_message += ", ";
            m_message += "\"" + registry.name(i) + "\"";
        }
        m_message += ".";
    }
    else
        m_message += "No models are registered.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid temporal model type
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Temporal model type that has been encountered.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_invalid_temporal::model_invalid_temporal(std::string origin,
                                                           std::string type,
                                                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid temporal model type \""+type+"\" encountered. " +
                message;

    // Add list of valid spectral models
    GModelTemporalRegistry registry;
    if (registry.size() > 0) {
        m_message += "The following models are registered: ";
        for (int i = 0; i < registry.size(); ++i) {
            if (i > 0)
                m_message += ", ";
            m_message += "\"" + registry.name(i) + "\"";
        }
        m_message += ".";
    }
    else
        m_message += "No models are registered.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid number of model parameters in XML element
 *
 * @param[in] origin Method that throws the error.
 * @param[in] xml XML element.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_invalid_parnum::model_invalid_parnum(std::string origin,
                                                       GXmlElement xml,
                                                       std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid number of model parameters found in XML element. " +
                message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid model parameter names in XML element
 *
 * @param[in] origin Method that throws the error.
 * @param[in] xml XML element.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_invalid_parnames::model_invalid_parnames(std::string origin,
                                                           GXmlElement xml,
                                                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid model parameter names found in XML element. " +
                message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model parameter has invalid scale
 *
 * @param[in] origin Method that throws the error.
 * @param[in] scale Scale.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_invalid_parscale::model_invalid_parscale(std::string origin,
                                                           double      scale,
                                                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Model parameter has invalid scale="+str(scale)+". " +
                message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model parameter has invalid scale
 *
 * @param[in] origin Method that throws the error.
 * @param[in] xml XML element.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_invalid_parscale::model_invalid_parscale(std::string origin,
                                                           GXmlElement xml,
                                                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Model parameter with invalid scale found in XML element. " +
                message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Not enough nodes in file function
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename File function filename.
 * @param[in] num Number of nodes.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::file_function_data::file_function_data(std::string origin,
                                                   std::string filename,
                                                   int         num,
                                                   std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "File function \""+filename+"\" contains "+str(num)+
                " energy nodes while at least 2 are required to describe a"
                " spectral shape.";
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Not enough columns in file function
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename File function filename.
 * @param[in] num Number of columns.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::file_function_columns::file_function_columns(std::string origin,
                                                         std::string filename,
                                                         int         num,
                                                         std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "File function \""+filename+"\" contains "+str(num)+
                " columns while at least 2 are required to define"
                " energy and intensity.";
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid value in file function
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename File function filename.
 * @param[in] value Number of columns.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::file_function_value::file_function_value(std::string origin,
                                                     std::string filename,
                                                     double      value,
                                                     std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid value \""+str(value)+"\" encountered in"
                " file function \""+filename+"\".";
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter not found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] name Parameter name.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::par_not_found::par_not_found(std::string origin,
                                         std::string name,
                                         std::string message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Parameter \""+name+"\" not found in container.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model not found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] name Model name.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::model_not_found::model_not_found(std::string origin,
                                             std::string name,
                                             std::string message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Model \""+name+"\" not found in container.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief No point source model component found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] name Model name.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_point_source::no_point_source(std::string origin,
                                             std::string name,
                                             std::string message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "No point source model component found in model \""+name+"\".";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief No extended source model component found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] name Model name.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_extended_source::no_extended_source(std::string origin,
                                                   std::string name,
                                                   std::string message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "No extended source model component found in model \""+name+"\".";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}
