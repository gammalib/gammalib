/***************************************************************************
 *                GSPIException.cpp  - SPI exception handler               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GSPIException.cpp
 * @brief SPI exception handler interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GSPIException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Exception message
 ***************************************************************************/
const char* GSPIExceptionHandler::what() const throw()
{
    // Set error message
    std::string message = "*** ERROR in " + m_origin + ": " + m_message;

    // Return message as C character array
    return tochar(message);
}
