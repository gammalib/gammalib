/***************************************************************************
 *        GModelSpatialPtsrc.cpp  -  Spatial point source model class      *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GModelSpatialPtsrc.cpp
 * @brief GModelSpatialPtsrc class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GModelSpatialPtsrc.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PAR                            "GModelSpatialPtsrc::par(int) const"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                GModelSpatialPtsrc constructors/destructors              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(void) : GModelSpatial()
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] dir Position of the point source.
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(GSkyDir dir) : GModelSpatial()
{
    // Initialise private members for clean destruction
    init_members();
    
    // Assign Right Ascension and Declination
    m_ra.value(dir.ra_deg());
    m_dec.value(dir.dec_deg());
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Model from which the instance should be built.
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(const GModelSpatialPtsrc& model) : 
   GModelSpatial(model)
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
GModelSpatialPtsrc::~GModelSpatialPtsrc()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        GModelSpatialPtsrc operators                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Model which should be assigned.
 ***************************************************************************/
GModelSpatialPtsrc& GModelSpatialPtsrc::operator= (const GModelSpatialPtsrc& model)
{ 
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatial::operator=(model);

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
 =                    GModelSpatialPtsrc public methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Get pointer to model parameter
 *
 * @param[in] index Parameter index.
 ***************************************************************************/
GModelPar* GModelSpatialPtsrc::par(int index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_npars)
        throw GException::out_of_range(G_PAR, index, 0, m_npars-1);
    
    // Return parameter pointer
    return m_par[index];
}


/*==========================================================================
 =                                                                         =
 =                    GModelSpatialPtsrc private methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatialPtsrc::init_members(void)
{
    // Initialise parameters
    m_npars  = 2;
    m_par[0] = &m_ra;
    m_par[1] = &m_dec;

    // Initialise Right Ascension
    m_ra = GModelPar();
    m_ra.name("RA");
    m_ra.unit("deg");
    m_ra.fix();
    
    // Initialise Declination
    m_dec = GModelPar();
    m_dec.name("DEC");
    m_dec.unit("deg");
    m_dec.fix();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpatialPtsrc members which should be copied.
 ***************************************************************************/
void GModelSpatialPtsrc::copy_members(const GModelSpatialPtsrc& model)
{
    // Copy model parameters (we do not need to copy the rest since it is
    // static)
    m_ra  = model.m_ra;
    m_dec = model.m_dec;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialPtsrc::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GModelSpatialPtsrc* GModelSpatialPtsrc::clone(void) const
{
    return new GModelSpatialPtsrc(*this);
}


/*==========================================================================
 =                                                                         =
 =                        GModelSpatialPtsrc friends                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put model in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] model Model to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelSpatialPtsrc& model)
{
    // Put observation in stream
    os << "=== GModelSpatialPtsrc ===" << std::endl;
    os << " Number of parameters ......: " << model.m_npars << std::endl;
    for (int i = 0; i < model.m_npars; ++i) {
        if (i > 0)
            os << std::endl;
        os << " Parameter .................: " << *(model.m_par[i]);
    }
        
    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                Other functions used by GModelSpatialPtsrc               =
 =                                                                         =
 ==========================================================================*/
