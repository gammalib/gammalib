/***************************************************************************
 *                  GCTAResponse.cpp  -  CTA Response class                *
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
 * @file GCTAResponse.cpp
 * @brief GCTAResponse class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <vector>
#include <string>
#include <unistd.h>           // access() function
#include "GCTAResponse.hpp"
#include "GCTAException.hpp"
#include "GTools.hpp"
#include "GVector.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ           "GCTAResponse::read_performance_table(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GCTAResponse::GCTAResponse(void) : GResponse()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
GCTAResponse::GCTAResponse(const GCTAResponse& rsp) : GResponse(rsp)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAResponse::~GCTAResponse(void)
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
 * @param[in] rsp Response to be assigned
 ***************************************************************************/
GCTAResponse& GCTAResponse::operator= (const GCTAResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GResponse::operator=(rsp);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(rsp);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return value of instrument response function.
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] instPntDir Instrument pointing direction (e.g. z-axis).
 * @param[in] instPosAng Instrument position angle.
 * @param[in] time Photon arrival time.
 ***************************************************************************/
double GCTAResponse::irf(const GSkyDir& obsDir, const GEnergy& obsEng,
                         const GSkyDir& srcDir, const GEnergy& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const GTime& time)
{
    // Get IRF components
    double irf;
    irf  = psf(obsDir, obsEng, srcDir, srcEng, instPntDir, instPosAng, time);
    irf *= aeff(obsDir, obsEng, srcDir, srcEng, instPntDir, instPosAng, time);
    irf *= edisp(obsDir, obsEng, srcDir, srcEng, instPntDir, instPosAng, time);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return effective area un units of cm2.
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] instPntDir Instrument pointing direction (e.g. z-axis).
 * @param[in] instPosAng Instrument position angle.
 * @param[in] time Photon arrival time.
 *
 * The actual implementation of this method assumes an effective area that
 * depends only on the true photon energy. No dependence on the photon's
 * arrival direction, observatory pointing and arrival time is assumed.
 ***************************************************************************/
double GCTAResponse::aeff(const GSkyDir& obsDir, const GEnergy& obsEng,
                          const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GSkyDir& instPntDir, const double& instPosAng,
                          const GTime& time)
{
    // Get log(E)
    double logE = log10(srcEng.TeV());
    
    // Interpolate effective area using node array
    m_nodes.set_value(logE);
    double aeff = m_aeff.at(m_nodes.inx_left())  * m_nodes.wgt_left() +
                  m_aeff.at(m_nodes.inx_right()) * m_nodes.wgt_right();

    // Convert from m2 to cm2
    aeff *= 10000.0;

    // Return effective area
    return aeff;
}


/***********************************************************************//**
 * @brief Return point spread function
 *
 * @param[in] obsDir Observed photon direction
 * @param[in] obsEng Observed energy of photon
 * @param[in] srcDir True photon direction
 * @param[in] srcEng True energy of photon
 * @param[in] instPntDir Instrument pointing direction (e.g. z-axis)
 * @param[in] instPosAng Instrument position angle
 * @param[in] time Photon arrival time
 *
 * The Point Spread Function defines the probability that a photon coming
 * from direction 'srcDir' is measured towards direction 'obsDir'.
 *
 * @todo Implement method (just a dummy for the moment)
 ***************************************************************************/
double GCTAResponse::psf(const GSkyDir& obsDir, const GEnergy& obsEng,
                         const GSkyDir& srcDir, const GEnergy& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const GTime& time)
{
    // Return Psf value
    return 1.0;
}


/***********************************************************************//**
 * @brief Return energy dispersion.
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] instPntDir Instrument pointing direction (e.g. z-axis).
 * @param[in] instPosAng Instrument position angle.
 * @param[in] time Photon arrival time.
 *
 * The actual implementation of this method assumes no energy dispersion,
 * which is equivalent of having a Dirac type energy dispersion.
 ***************************************************************************/
double GCTAResponse::edisp(const GSkyDir& obsDir, const GEnergy& obsEng,
                           const GSkyDir& srcDir, const GEnergy& srcEng,
                           const GSkyDir& instPntDir, const double& instPosAng,
                           const GTime& time)
{
    // Dirac energy dispersion
    double edisp = (obsEng == srcEng) ? 1.0 : 0.0;
    
    // Return energy dispersion
    return edisp;
}


/***********************************************************************//**
 * @brief Set the path to the calibration database.
 *
 * @param[in] caldb Absolute path to calibration database
 *
 * @todo Implement checking 
 ***************************************************************************/
#define G_SET_CALDB                   "GCTAResponse::set_caldb(std::string&)"
void GCTAResponse::set_caldb(const std::string& caldb)
{
    // Check if calibration database directory is accessible
    if (access(caldb.c_str(), R_OK) != 0)
        throw GException::caldb_not_found(G_SET_CALDB, caldb);
    
    // Store the path to the calibration database
    m_caldb = caldb;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load CTA response.
 *
 * @param[in] irfname Name of CTA response (without any file extension).
 *
 * The actually dummy version of the CTA response loads a CTA performance
 * table given in ASCII format into memory.
 ***************************************************************************/
void GCTAResponse::load(const std::string& irfname)
{
    // Initialise response members
    free_members();
    init_members();
    
    // Build filename
    std::string filename = m_caldb + "/" + irfname + ".dat";

    // Read performance table
    read_performance_table(filename);

    // Store response name
    m_rspname = irfname;

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
void GCTAResponse::init_members(void)
{
    // Initialise members
    m_logE.clear();
    m_aeff.clear();
    m_r68.clear();
    m_r80.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
void GCTAResponse::copy_members(const GCTAResponse& rsp)
{
    // Copy attributes
    m_nodes = rsp.m_nodes;
    m_logE  = rsp.m_logE;
    m_aeff  = rsp.m_aeff;
    m_r68   = rsp.m_r68;
    m_r80   = rsp.m_r80;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAResponse::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone response class
***************************************************************************/
GCTAResponse* GCTAResponse::clone(void) const
{
    return new GCTAResponse(*this);
}


/***********************************************************************//**
 * @brief Read CTA performance table
 *
 * @param[in] filename Filename of CTA performance table.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method reads a CTA performance table given in the format that is
 * distributed within the CTA collaboration.
 ***************************************************************************/
void GCTAResponse::read_performance_table(const std::string& filename)
{
    // Allocate line buffer
    const int n = 1000; 
    char  line[n];

    // Open performance table readonly
    FILE* fptr = fopen(filename.c_str(), "r");
    if (fptr == NULL)
        throw GCTAException::file_open_error(G_READ, filename);

    // Read lines
    while (fgets(line, n, fptr) != NULL) {

        // Split line in elements
        std::vector<std::string> elements = split(line, " ");
        for (std::vector<std::string>::iterator it = elements.begin();
             it != elements.end(); ++it) {
            if (strip_whitespace(*it).length() == 0)
                elements.erase(it);
        }
        
        // Skip header
        if (elements[0].find("log(E)") != std::string::npos)
            continue;

        // Break loop if end of data table has been reached
        if (elements[0].find("----------") != std::string::npos)
            break;

        // Push elements in vectors
        m_logE.push_back(todouble(elements[0]));
        m_aeff.push_back(todouble(elements[1]));
        m_r68.push_back(todouble(elements[2]));
        m_r80.push_back(todouble(elements[3]));

    } // endwhile: looped over lines
    
    // If we have nodes then setup node array
    int num = m_logE.size();
    if (num > 0) {
        GVector logE(num);
        for (int i = 0; i < num; ++i)
            logE(i) = m_logE.at(i);
        m_nodes.nodes(logE);
    }

    // Close file
    fclose(fptr);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] GCTAResponse Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTAResponse& rsp)
{
    // Put object in stream
    os << "=== GCTAResponse ===" << std::endl;
    os << " Calibration database ......: " << rsp.m_caldb << std::endl;
    os << " Name of instrument response: " << rsp.m_rspname << std::endl;
    os << " Response definiton ........: ";
    for (int i = 0; i < rsp.m_logE.size(); ++i) {
        os << std::endl
           << "  logE=" << rsp.m_logE.at(i)
           << ": Aeff=" << rsp.m_aeff.at(i)
           << " m2, r68=" << rsp.m_r68.at(i)
           << " deg, r80=" << rsp.m_r80.at(i)
           << " deg";
    }

    // Return output stream
    return os;
}
