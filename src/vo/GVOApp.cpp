/***************************************************************************
 *                     GVOApp.hpp - VO SAMP Hub class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Thierry Louge                                    *
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
 * @file GVOApp.cpp
 * @brief VO  Application connected to SAMP HUB class implementation
 * @author Thierry Louge
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#include <errno.h>
#endif
#include <cstdlib>         // std::getenv() function
#include <cstdio>          // std::fopen(), etc. functions
#include <cstring>         // std::memset() function
#include <unistd.h>        // close() function
#include "GVOApp.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_HANDLE_REQUEST                  "GVOApp::handle_request(socklen_t)"
#define G_START_HUB                                     "GVOApp::start_hub()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GVOApp::GVOApp(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] hub VO Hub.
 ***************************************************************************/
GVOApp::GVOApp(const GVOApp& hub)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(hub);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GVOApp::~GVOApp(void)
{
    
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
 * @param[in] hub VO hub.
 * @return VO hub.
 ***************************************************************************/
GVOApp& GVOApp::operator=(const GVOApp& hub)
{
    // Execute only if object is not identical
    if (this != &hub) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(hub);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object.
 *
 * Reset object to a clean initial state.
 ***************************************************************************/
void GVOApp::clear(void)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GVOApp* GVOApp::clone(void) const
{
    // Clone client
    return new GVOApp(*this);
}


/***********************************************************************//**
 * @brief Creates nex application container
 ***************************************************************************/
void GVOApp::start_app(void)
{
    // Creates Hub
    printf("Creating new application Container\n");
    init_app();

    // Return
    return;
}

