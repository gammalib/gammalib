/***************************************************************************
 *                   GSPITools.i - INTEGRAL/SPI tools                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPITools.i
 * @brief Definition of INTEGRAL/SPI tools
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSPITools.hpp"
%}

/* __ Constants __________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    const GFitsTable* spi_hdu(const GFits&       fits,
                              const std::string& extname,
                              const int&         extver = 1);
    int               spi_num_hdus(const GFits&       fits,
                                   const std::string& extname);
}
