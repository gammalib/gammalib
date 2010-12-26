/***************************************************************************
 *          GException_obs.cpp  -  observations exception handlers         *
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
 * @file GException_obs.cpp
 * @brief Observation exception handler interface implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief No valid instrument response.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_response::no_response(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "Instrument response not defined."
                " Define instrument response before using this method."+
                message;
    return;
}


/***********************************************************************//**
 * @brief No valid region of interest.
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_roi::no_roi(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "Region of interest not defined."
                " Define region of interest before using this method."+
                message;
    return;
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


/***********************************************************************//**
 * @brief Calibration database directory not found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] caldb Name of calibration database.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::caldb_not_found::caldb_not_found(std::string origin,
                                             std::string caldb,
                                             std::string message)
{
    m_origin  = origin;
    m_message = "Calibration database \""+caldb+"\" not found. "+message;
    return;
}


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
    m_message = "Invalid response type '"+type+"' specified.";
}


/***********************************************************************//**
 * @brief Good Time Intervals are not valid.
 *
 * @param[in] origin Method that throws the error.
 * @param[in] start Good Time Interval start time.
 * @param[in] start Good Time Interval stop time.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::gti_invalid::gti_invalid(std::string origin, GTime tstart,
                                     GTime tstop, std::string message)
{
    std::ostringstream s_tstart;
    std::ostringstream s_tstop;
    s_tstart  << tstart;
    s_tstop   << tstop;
    m_origin  = origin;
    m_message = "Invalid Good Time Interval ("+s_tstart.str()+"-"+
                s_tstop.str()+") specified. "+message;
}


/***********************************************************************//**
 * @brief Energy range is not valid.
 *
 * @param[in] origin Method that throws the error.
 * @param[in] emin Minimum energy in MeV.
 * @param[in] emax Maximum energy in MeV.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::erange_invalid::erange_invalid(std::string origin, double emin, 
                                           double emax, std::string message)
{
    m_origin  = origin;
    m_message = "Invalid energy range (Emin="+str(emin)+" MeV, Emax="+
                str(emax)+" MeV. "+message;
}


/***********************************************************************//**
 * @brief Invalid optimization statistics
 *
 * @param[in] origin Method that throws the error.
 * @param[in] statistics Encountered statistics.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::invalid_statistics::invalid_statistics(std::string origin,
                                                   std::string statistics,
                                                   std::string message)
{
    m_origin  = origin;
    m_message = "Invalid optimization statistics \""+statistics+
                "\" encountered. "+message;
}
