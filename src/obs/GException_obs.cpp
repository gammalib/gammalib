/***************************************************************************
 *              GException_obs.cpp  -  obs exception handlers              *
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
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Response error: invalid response type specified
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Specified response type.
 ***************************************************************************/
GException::rsp_invalid_type::rsp_invalid_type(std::string origin,
                                               std::string type)
{
    m_origin  = origin;
    m_message = "Invalid response type '"+type+"' specified";
}


/***********************************************************************//**
 * @brief Mismatch between gradient vector and number of model parameters
 *
 * @param[in] origin Method that throws the error.
 * @param[in] nsize Gradient vector size.
 * @param[in] npars Number of parameters in model.
 ***************************************************************************/
GException::gradient_par_mismatch::gradient_par_mismatch(std::string origin,
                                                         int         nsize,
                                                         int         npars)
{
    m_origin  = origin;
    m_message = "Gradient vector size "+str(nsize)+
                " mismatches number "+str(npars)+" of model parameters.";
}
