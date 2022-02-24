/***************************************************************************
 *               gammalibd.cpp - GammaLib daemon launcher                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file gammalibd.cpp
 * @brief GammaLib daemon launcher
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <fcntl.h>         // for file locking
#include <unistd.h>        // access() function
#include <cstdlib>         // exit() function
#include <cstdio>          // std::fopen(), etc. functions
#include <signal.h>        // signal() function
#include "GammaLib.hpp"    // make GammaLib classes available


/***********************************************************************//**
 * @brief GammaLib daemon launcher
 *
 * This executable launched the GammaLib daemon. The code is a copy of the
 * code implemented in the GApplication::start_daemon() method. Calling this
 * executable during the initialisation of GammaLib will launch a daemon
 * that is independent of any specific application, and that has the correct
 * process name on all platforms.
 ***************************************************************************/
int main(void) {

    // Initialise daemon
    GDaemon daemon;

    // If daemon is not alive then create instance
    if (!daemon.alive()) {

        // Create child process to start the daemon. Do nothing if child
        // process creation fails.
        int pid = fork();
        if (pid >= 0) {

            // If we have a PID of 0 we are in the child process. In this
            // case we create and start the daemon ...
            if (pid == 0) {

                // The child process becomes session leader
                setsid();

                // Ignore signals
                signal(SIGCHLD, SIG_IGN); // Ignore if child has stopped
                signal(SIGHUP,  SIG_IGN); // Ignore death of controlling process

                // Fork for the second time and let the first fork
                // process terminate
                pid = fork();
                if (pid == 0) {
                    GDaemon daemon;
                    daemon.start();
                    exit(EXIT_SUCCESS);
                }
                else {
                    exit(EXIT_SUCCESS);
                }

            }

        } // endif: child proces created

    } // endif: daemon

    // Always exit with success
    return EXIT_SUCCESS;
}
