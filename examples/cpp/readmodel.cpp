/***************************************************************************
 *             readmodel.cpp - Illustrates how to read a model             *
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
 * @file readmodel.cpp
 * @brief Illustrates how to read a model
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GammaLib.hpp"


/***********************************************************************//**
 * @brief Dummy read model test
 ***************************************************************************/
int main(void) {
    GModels models("source.xml");
    std::cout << models << std::endl;
    return 0;
}

