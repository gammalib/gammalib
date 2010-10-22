/***************************************************************************
 *            GException_model.cpp  -  Model exception handlers            *
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

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
//#include "GTools.hpp"


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
