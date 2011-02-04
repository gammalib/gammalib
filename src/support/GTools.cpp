/***************************************************************************
 *                          GTools.cpp  -  GammaLib tools                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>        // access() function
#include <sys/stat.h>
#include <cmath>
#include <cfloat>
#include <cctype>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "GTools.hpp"

/* __ Coding definitions _________________________________________________ */
#define G_PARFORMAT_LENGTH 29


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
 * @brief Convert unsigned short integer value into string
 *
 * @param[in] value Unsigned short integer to be converted into string.
 ***************************************************************************/
std::string str(const unsigned short int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Convert unsigned integer value into string
 *
 * @param[in] value Unsigned integer to be converted into string.
 ***************************************************************************/
std::string str(const unsigned int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Convert unsigned long integer value into string
 *
 * @param[in] value Unsigned long integer to be converted into string.
 ***************************************************************************/
std::string str(const unsigned long int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Convert unsigned long long integer value into string
 *
 * @param[in] value Unsigned long long integer to be converted into string.
 ***************************************************************************/
std::string str(const unsigned long long int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Convert short integer value into string
 *
 * @param[in] value Short integer to be converted into string.
 ***************************************************************************/
std::string str(const short int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
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
 * @brief Convert long integer value into string
 *
 * @param[in] value Long integer to be converted into string.
 ***************************************************************************/
std::string str(const long int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}


/***********************************************************************//**
 * @brief Convert long long integer value into string
 *
 * @param[in] value Long long integer to be converted into string.
 ***************************************************************************/
std::string str(const long long int& value)
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
 * @brief Convert string to C string
 *
 * @param[in] arg String to be converted.
 *
 * Allocates a C string with the content of a C++ string.
 ***************************************************************************/
char* tochar(const std::string& arg)
{
    // Allocate C string
    char* str = new char[arg.length()+1];

    // Copy characters
    for (std::size_t i = 0; i < arg.length(); ++i)
        str[i] = arg[i];

    // Set line end character
    str[arg.length()] = '\0';

    // Return C string
    return str;
}


/***********************************************************************//**
 * @brief Convert string into short value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
short toshort(const std::string& arg)
{
    std::istringstream iss(arg);
    short              result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into unsigned short value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
unsigned short toushort(const std::string& arg)
{
    std::istringstream iss(arg);
    unsigned short     result;
    iss >> std::dec >> result;
    return result;
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
 * @brief Convert string into unsigned integer value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
unsigned int touint(const std::string& arg)
{
    std::istringstream iss(arg);
    unsigned int       result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into long value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
long tolong(const std::string& arg)
{
    std::istringstream iss(arg);
    long               result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into unsigned long value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
unsigned long toulong(const std::string& arg)
{
    std::istringstream iss(arg);
    unsigned long      result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into long long value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
long long tolonglong(const std::string& arg)
{
    std::istringstream iss(arg);
    long long          result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into unsigned long long value
 *
 * @param[in] arg String to be converted.
 ***************************************************************************/
unsigned long long toulonglong(const std::string& arg)
{
    std::istringstream iss(arg);
    unsigned long long result;
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
    std::size_t pos = 0;
    std::size_t len = s.size();

    // Loop over string
    while (pos < len && pos != std::string::npos) {
        std::size_t index = s.find_first_of(sep, pos);
        std::size_t n     = std::string::npos;
        if (index != std::string::npos)
            n = index-pos;
        if (n > 0)
            result.push_back(s.substr(pos, n));
        pos = (index != std::string::npos) ? index + 1 : std::string::npos;
    } // endwhile: there were still characters in the string

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Fill string with n strings of same type
 *
 * @param[in] s String to be filled.
 * @param[in] n Number of fillings.
 *
 * Replicates a given string n time.
 ***************************************************************************/
std::string fill(const std::string& s, int n)
{
    // Initialise result
    std::string result = "";

    // Replicate string
    for (int i = 0; i < n; ++i)
        result.append(s);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Left justify string to achieve a length of n characters
 *
 * @param[in] s String to be centred.
 * @param[in] n Requested total width.
 * @param[in] c Fill character.
 *
 * Left justify string by adding whitespace to the right to achieve a length
 * of n characters.
 ***************************************************************************/
std::string left(const std::string& s, int n, char c)
{
    // Compute number of characters to fill right
    int n_right  = n - s.length();

    // Set result
    std::string result = s + fill(std::string(1,c), n_right);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Right justify string to achieve a length of n characters
 *
 * @param[in] s String to be centred.
 * @param[in] n Requested total width.
 *
 * Right justify string by adding whitespace to the left to achieve a length
 * of n characters.
 ***************************************************************************/
std::string right(const std::string& s, int n, char c)
{
    // Compute number of characters to fill right
    int n_left  = n - s.length();

    // Set result
    std::string result = fill(std::string(1,c), n_left) + s;

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Center string to achieve a length of n characters
 *
 * @param[in] s String to be centred.
 * @param[in] n Requested total width.
 *
 * Center string by adding whitespace to the left and the right to achieve a
 * length of n characters.
 ***************************************************************************/
std::string center(const std::string& s, int n, char c)
{
    // Compute number of characters to fill left and right
    int n_right = (n-s.length()) / 2;
    int n_left  = n - s.length() - n_right;

    // Set result
    std::string result = fill(std::string(1,c), n_left)+s+fill(" ", n_right);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Convert string in parameter format
 *
 * @param[in] s String to be converted.
 *
 * Converts and string into the parameter format of type "s ......: " with a
 * total length of G_PARFORMAT_LENGTH.
 ***************************************************************************/
std::string parformat(const std::string& s)
{
    // Compute number of characters to fill right. Do not clip the string if
    // it is too long since we do not want to loose information.
    int n_right  = G_PARFORMAT_LENGTH - s.length() - 3;

    // Set result
    std::string result = " " + s + " " + fill(".", n_right) + ": ";

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
    return (v1 >= 0) ? ((v1 < v2) ? v1 : std::fmod(v1,v2)) : (std::fmod(v1,v2)+v2);
}


/***********************************************************************//**
 * @brief Checks if a file exists.
 *
 * @param[in] filename File name.
 *
 * Checks if a file exists. If a directory with the same name if found, false
 * is returned.
 ***************************************************************************/
bool file_exists(const std::string& filename)
{
    // Initialise result
    bool result = false;

    // Allocate file information structure
    struct stat info;

    // Get file information structure
    int ret = stat(filename.c_str(), &info);

    // Check if file is a regular file
    if (ret == 0 && S_ISREG(info.st_mode))
        result = true;

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Checks if argument is infinite
 *
 * @param[in] x Argument.
 *
 * This function has been copied from gnulib.
 ***************************************************************************/
bool isinfinite(double x)
{
  return (x < -DBL_MAX || x > DBL_MAX);
}
