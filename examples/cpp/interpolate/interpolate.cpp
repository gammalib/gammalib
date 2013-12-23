/***************************************************************************
 *              interpolate.cpp - Illustrates interpolation                *
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
 * @file interpolate.cpp
 * @brief Illustrates interpolation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @brief Interpolation between nodes
 *
 * This code illustrates the interpolation between nodes.
 ***************************************************************************/
int main(void) {

    // Create some nodes
    double x_i[] = {1.0, 4.0, 6.0};
    double y_i[] = {8.0, 7.0, 2.0};
    
    // Create node array from x_i values
    GNodeArray nodes(3, x_i);

    // Get interpolated values
    for (double x = 0; x < 10.0; x += 0.5) {
        nodes.set_value(x);
        double y = y_i[nodes.inx_left()]  * nodes.wgt_left() +
                   y_i[nodes.inx_right()] * nodes.wgt_right();
        std::cout << "x=" << x << " : y=" << y << std::endl;
    }

    // Exit
    return 0;
}
