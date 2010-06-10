/***************************************************************************
 *                          GTools.cpp  -  GammaLib tools                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GTools.cpp
 * @brief Gammalib tools implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cctype>
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Strip leading and trailing whitespace from string
 *
 * @param[in] arg String from which whitespace should be stripped.
 ***************************************************************************/
std::string strip_whitespace(const std::string& arg)
{
    // Return result
    return (strip_chars(arg, " "));
}


/***********************************************************************//**
 * @brief Strip leading and trailing character from string
 *
 * @param[in] arg String from which character should be stripped.
 * @param[in] chars Character(s) to be stripped.
 ***************************************************************************/
std::string strip_chars(const std::string& arg, const std::string& chars)
{
    // Initialise empty result string
    std::string result;
    
    // Continue only if argument is not empty
    if (!arg.empty()) {

        // Get start and stop
        std::string::size_type start = arg.find_first_not_of(chars);
        std::string::size_type stop  = arg.find_last_not_of(chars);

        // Continue only is start and stop are valid
        if (start != std::string::npos && stop != std::string::npos) {

            // Continue only if stop is larger then or equal to  start
            if (start <= stop)
                result = arg.substr(start, stop-start+1);

        } // endif: start and stop were valid

    } // endif: argument was not empty

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Convert integer value into string
 *
 * @param[in] value Integer to be converted into string.
 ***************************************************************************/
std::string str(const int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Convert single precision value into string
 *
 * @param[in] value Single precision value to be converted into string.
 ***************************************************************************/
std::string str(const float& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Convert double precision value into string
 *
 * @param[in] value Double precision value to be converted into string.
 ***************************************************************************/
std::string str(const double& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Convert string into integer value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
int toint(const std::string& arg)
{
    std::istringstream iss(arg);
    int                result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into float value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
float tofloat(const std::string& arg)
{
    std::istringstream iss(arg);
    float              result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into double value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
double todouble(const std::string& arg)
{
    std::istringstream iss(arg);
    double             result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string to upper case
 *
 * @param[in] arg String to be converted to upper case.
 ***************************************************************************/
std::string toupper(const std::string& arg)
{
    std::string s = arg;
    std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) toupper);
    return s;
}


/***********************************************************************//**
 * @brief Convert string to lower case
 *
 * @param[in] arg String to be converted to upper case.
 ***************************************************************************/
std::string tolower(const std::string& arg)
{
    std::string s = arg;
    std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
    return s;
}


/***********************************************************************//**
 * @brief Split string
 *
 * @param[in] s String to be splitted.
 * @param[in] sep Separator(s).
 ***************************************************************************/
std::vector<std::string> split(const std::string& s, const std::string& sep)
{
    // Allocate result string vector
    std::vector<std::string> result;
    
    // Initialise counters
    size_t pos = 0;
    size_t len = s.size();
    
    // Loop over string
    while (pos < len && pos != std::string::npos) {
        size_t index = s.find_first_of(sep, pos);
        size_t n     = std::string::npos;
        if (index != std::string::npos)
            n = index-pos;
        result.push_back(s.substr(pos, n));
        pos = (index != std::string::npos) ? index + 1 : std::string::npos;
    } // endwhile: there were still characters in the string
    
    // Return result
    return result;
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
