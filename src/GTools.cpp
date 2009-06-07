/***************************************************************************
 *                          GTools.cpp  -  GammaLib tools                  *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <sstream>
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Strip leading and trailing whitespace from string
 *
 * @param[in] arg String which should be stripped of whitespace.
 ***************************************************************************/
std::string strip_whitespace(const std::string& arg)
{
    // Return empty string if argument is empty
    if (arg.empty())
        return arg;

    // Get start and stop
    std::string::size_type start = arg.find_first_not_of(' ');
    std::string::size_type stop  = arg.find_last_not_of(' ');

    // If there is only whitespace in string then return empty string.
    // Otherwise strip off whitespace
    std::string result;
    if (start <= stop)
        result = arg.substr(start, stop-start+1);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Convert integer value into string
 *
 * @param[in] value Value to be converted into a string.
 ***************************************************************************/
std::string str(const int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Returns the remainder of the division \a v1/v2.
 *
 * @param[in] v1 Argument 1.
 * @param[in] v2 Argument 2.
 *
 * Returns the remainder of the division \a v1/v2.
 * The result is non-negative.
 * \a v1 can be positive or negative; \a v2 must be positive.
 ***************************************************************************/
double modulo(double v1, double v2)
{
    // Return
    return (v1 >= 0) ? ((v1 < v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2);
}
