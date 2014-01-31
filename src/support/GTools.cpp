/***************************************************************************
 *                       GTools.cpp - GammaLib tools                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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
 * @file GTools.cpp
 * @brief Gammalib tools implementation
 * @author Juergen Knoedlseder
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
#include <cstdlib>         // std::getenv() function
#include <cstring>         // std::strlen() function
#include <iostream>
#include <sstream>
#include <algorithm>
#include "GTools.hpp"

/* __ Compile options ____________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_PARFORMAT_LENGTH 29


/***********************************************************************//**
 * @brief Strip leading and trailing whitespace from string
 *
 * @param[in] arg String from which whitespace should be stripped.
 ***************************************************************************/
std::string gammalib::strip_whitespace(const std::string& arg)
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
std::string gammalib::strip_chars(const std::string& arg,
                                  const std::string& chars)
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
            if (start <= stop) {
                result = arg.substr(start, stop-start+1);
            }

        } // endif: start and stop were valid

    } // endif: argument was not empty

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Expand environment variables in string
 *
 * @param[in] arg String
 *
 * Expands any environment variable that is found in a string. Valid
 * delimiters for environment variables are $ENV{\<name\>}, $ENV(\<name\>),
 * ${\<name\>}, $(\<name\>) and $\<name\> (in the last case the terminating
 * delimiter is either a / or a blank character or the end of the string).
 * Environment variables occuring within single quotes (') are ignored.
 *
 * This function has been inspired by the function ape_util_expand_env_var
 * from ape_util.c in the ape software developed at HEASARC.
 ***************************************************************************/
std::string gammalib::expand_env(const std::string& arg)
{
    // Set environment variable delimiters
    static const char* begin_delim[] = { "$ENV{", "$ENV(", "${", "$(", "$" };
    static const char* end_delim[]   = { "}", ")", "}", ")", "/" };
    static const int   num_delim     = 5;

    // Initialise result with argument
    std::string result = arg;

    // Initialise parse parameters
    size_t index    = 0;
    bool   in_quote = false;

    // Loop over string
    while (index < result.length()) {

        // If we have an escaped character then skip the current character
        // and the next one
        if (result[index] == '\\') {
            index += 2;
            continue;
        }

        // If we have a single quote then toggle the quote state
        if (result[index] == '\'') {
            in_quote = !in_quote;
            index++;
            continue;
        }

        // Don't expand environment variables inside single quotes. Note
        // that double quotes are ok.
        if (in_quote) {
            index++;
            continue;
        }

        // Find delimiter which indicates beginning of an environment variable
        size_t begin_length = 0;
        int    delim_idx    = 0;
        for (; delim_idx < num_delim; ++delim_idx) {
            size_t len = std::strlen(begin_delim[delim_idx]);
            if (result.compare(index, len, begin_delim[delim_idx]) == 0) {
                begin_length = len;
                break;
            }
        }

        // If we found a delimiter then process the environment variable
        if (begin_length > 0) {

            // Search for the termination delimiter of the environment
            // variable. There is a special case for delimiter 4:
            // It has always an end_length of zero as the / is not a real
            // delimiter, but just an indicator that the environment variable
            // ends. Another indicator is a blank. If the end of the
            // string has been reached this is also acceptable.
            size_t i_start    = index + begin_length;
            size_t i_end      = i_start;
            size_t end_length = 0;
            if (delim_idx == 4) {
                while (i_end < result.length() &&
                       result.compare(i_end, 1, "/") != 0 &&
                       result.compare(i_end, 1, " ") != 0) {
                    i_end++;
                }
            }
            else {
                end_length = std::strlen(end_delim[delim_idx]);
                while (i_end < result.length() &&
                       result.compare(i_end, end_length, end_delim[delim_idx]) != 0) {
                    i_end++;
                }
            }

            // If termination delimiter has been found then expand the
            // environment variable
            if (i_end < result.length() || delim_idx == 4) {

                // Extract environment variable name
                std::string name        = result.substr(i_start, i_end-i_start);
                size_t      name_length = name.length();

                // Erase delimiters and environment variable
                result.erase(index, begin_length+name_length+end_length);

                // Get the environment variable from the operating system
                const char* env = std::getenv(name.c_str());

                // If the environment variable has been found then replace
                // it by its value
                if (env != NULL) {

                    // Set replacement string and its length
                    std::string replace(env);
                    size_t      replace_length = replace.length();

                    // Insert replacement string
                    result.insert(index, replace);

                    // Advance pointer
                    index += replace_length;

                } // endif: environment variable has been found

            } // endif: termination delimiter found

            // If no environment variable has been found then set index=i_end+1
            else {
                index = i_end + 1;
            }

        } // endif: we found an environment variable delimiter

        // ... otherwise advance to next character
        else {
            index++;
        }

    } // endwhile: looped over string

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Build file path from path name and file name
 *
 * @param[in] pathname Path name.
 * @param[in] filename File name.
 * @return File path.
 *
 * Builds a file path by combining the @p pathname and the @p filename
 * following
 *
 *      filepath = pathname/filename
 *
 * If @p pathname is an empty string, the method simply returns the
 * @p filename.
 ***************************************************************************/
std::string gammalib::filepath(const std::string& pathname,
                               const std::string& filename)
{
    // Initialise filepath
    std::string filepath;

    // If path name is empty, simply return the file name
    if (pathname.empty()) {
        filepath = filename;
    }

    // ... otherwise combine both
    else {
        filepath = pathname + "/" + filename;
    }

    // Return the file path
    return filepath;
}
    

/***********************************************************************//**
 * @brief Convert unsigned short integer value into string
 *
 * @param[in] value Unsigned short integer to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const unsigned short int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert unsigned integer value into string
 *
 * @param[in] value Unsigned integer to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const unsigned int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert unsigned long integer value into string
 *
 * @param[in] value Unsigned long integer to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const unsigned long int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert unsigned long long integer value into string
 *
 * @param[in] value Unsigned long long integer to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const unsigned long long int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert short integer value into string
 *
 * @param[in] value Short integer to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const short int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert integer value into string
 *
 * @param[in] value Integer to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert long integer value into string
 *
 * @param[in] value Long integer to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const long int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert long long integer value into string
 *
 * @param[in] value Long long integer to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const long long int& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert single precision value into string
 *
 * @param[in] value Single precision value to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const float& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert double precision value into string
 *
 * @param[in] value Double precision value to be converted into string.
 ***************************************************************************/
std::string gammalib::str(const double& value)
{
    std::ostringstream s_value;
    s_value << value;
    return s_value.str();
}


/***********************************************************************//**
 * @brief Convert string to C string
 *
 * @param[in] arg String to be converted.
 *
 * Allocates a C string with the content of a C++ string.
 ***************************************************************************/
char* gammalib::tochar(const std::string& arg)
{
    // Allocate C string
    char* str = new char[arg.length()+1];

    // Copy characters
    for (std::size_t i = 0; i < arg.length(); ++i) {
        str[i] = arg[i];
    }

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
short gammalib::toshort(const std::string& arg)
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
unsigned short gammalib::toushort(const std::string& arg)
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
int gammalib::toint(const std::string& arg)
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
unsigned int gammalib::touint(const std::string& arg)
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
long gammalib::tolong(const std::string& arg)
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
unsigned long gammalib::toulong(const std::string& arg)
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
long long gammalib::tolonglong(const std::string& arg)
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
unsigned long long gammalib::toulonglong(const std::string& arg)
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
float gammalib::tofloat(const std::string& arg)
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
double gammalib::todouble(const std::string& arg)
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
std::string gammalib::toupper(const std::string& arg)
{
    std::string s = arg;
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}


/***********************************************************************//**
 * @brief Convert string to lower case
 *
 * @param[in] arg String to be converted to upper case.
 ***************************************************************************/
std::string gammalib::tolower(const std::string& arg)
{
    std::string s = arg;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}


/***********************************************************************//**
 * @brief Split string
 *
 * @param[in] s String to be splitted.
 * @param[in] sep Separator(s).
 ***************************************************************************/
std::vector<std::string> gammalib::split(const std::string& s,
                                         const std::string& sep)
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
        if (index != std::string::npos) {
            n = index-pos;
        }
        if (n > 0) {
            result.push_back(s.substr(pos, n));
        }
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
std::string gammalib::fill(const std::string& s, const int& n)
{
    // Initialise result
    std::string result = "";

    // Replicate string
    for (int i = 0; i < n; ++i) {
        result.append(s);
    }

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
std::string gammalib::left(const std::string& s, const int& n, const char& c)
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
 * @param[in] c Fill character.
 *
 * Right justify string by adding whitespace to the left to achieve a length
 * of n characters.
 ***************************************************************************/
std::string gammalib::right(const std::string& s, const int& n, const char& c)
{
    // Compute number of characters to fill right
    int n_left  = n - s.length();

    // Set result
    std::string result = fill(std::string(1,c), n_left) + s;

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Centre string to achieve a length of n characters
 *
 * @param[in] s String to be centred.
 * @param[in] n Requested total width.
 * @param[in] c Left fill character.
 *
 * Centre string by adding whitespace to the left and the right to achieve a
 * length of n characters.
 ***************************************************************************/
std::string gammalib::centre(const std::string& s, const int& n, const char& c)
{
    // Compute number of characters to fill left and right
    int n_right = (n-s.length()) / 2;
    int n_left  = n - s.length() - n_right;

    // Set result
    std::string result = fill(std::string(1,c), n_left) + s + 
                         fill(std::string(1,c), n_right);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Convert string in parameter format
 *
 * @param[in] s String to be converted.
 * @param[in] indent Indentation of parameter (default: 0).
 *
 * Converts and string into the parameter format of type "s ......: " with a
 * total length of G_PARFORMAT_LENGTH.
 ***************************************************************************/
std::string gammalib::parformat(const std::string& s, const int& indent)
{
    // Compute number of characters to fill right. Do not clip the string if
    // it is too long since we do not want to loose information.
    int n_right = G_PARFORMAT_LENGTH - s.length() - 3 - indent;

    // Set result
    std::string result = " " + s + " " + fill(".", n_right) + ": ";

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Compute photon flux between two energies for a power law
 *
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @param[in] epivot Pivot energy.
 * @param[in] gamma Spectral index.
 *
 * Analytically computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} 
 *                             \left( E/E_{\rm pivot} \right)^{\gamma} dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$E_{\rm pivot}\f$ is the pivot energy, and
 * \f$\gamma\f$ is the spectral index.
 ***************************************************************************/
double gammalib::plaw_photon_flux(const double& emin, const double& emax,
                                  const double& epivot, const double& gamma)
{
    // Initialise flux
    double flux = 0.0;

    // Continue only if emax > emin
    if (emax > emin) {

        // Compute photon flux. Computations dependend on the exponent. We
        // add here a kluge to assure numerical accuracy.
        flux = std::pow(epivot, -gamma);
        double exponent = gamma + 1.0;
        if (std::abs(exponent) > 1.0e-11) {
            flux *= (std::pow(emax, exponent) -
                     std::pow(emin, exponent)) / exponent;
        }
        else {
            flux *= (std::log(emax) - std::log(emin));
        }

    } // endif: emax > emin

    // Return result
    return flux;
}


/***********************************************************************//**
 * @brief Compute energy flux between two energies for a power law
 *
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @param[in] epivot Pivot energy.
 * @param[in] gamma Spectral index.
 *
 * Analytically computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} 
 *                           \left( E/E_{\rm pivot} \right)^{\gamma} E dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$E_{\rm pivot}\f$ is the pivot energy, and
 * \f$\gamma\f$ is the spectral index.
 ***************************************************************************/
double gammalib::plaw_energy_flux(const double& emin, const double& emax,
                                  const double& epivot, const double& gamma)
{
    // Initialise flux
    double flux = 0.0;

    // Continue only if emax > emin
    if (emax > emin) {

        // Compute energy flux. Computations dependend on the exponent. We
        // add here a kluge to assure numerical accuracy.
        flux = std::pow(epivot, -gamma);
        double exponent = gamma + 2.0;
        if (std::abs(exponent) > 1.0e-11) {
            flux *= (std::pow(emax, exponent) -
                     std::pow(emin, exponent)) / exponent;
        }
        else {
            flux *= (std::log(emax) - std::log(emin));
        }

    } // endif: emax > emin

    // Return result
    return flux;
}


/***********************************************************************//**
 * @brief Checks if file exists
 *
 * @param[in] filename File name.
 *
 * Checks if a file exists. If a directory with the same name if found, false
 * is returned. The function expands any environment variable prior to
 * checking.
 ***************************************************************************/
bool gammalib::file_exists(const std::string& filename)
{
    // Initialise result
    bool result = false;

    // Allocate file information structure
    struct stat info;

    // Get file information structure
    int ret = stat(gammalib::expand_env(filename).c_str(), &info);

    // Check if file is a regular file
    if (ret == 0 && S_ISREG(info.st_mode)) {
        result = true;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Checks if directory exists
 *
 * @param[in] dirname Directory name.
 *
 * Checks if a directory exists. The function expands any environment
 * variable prior to checking.
 ***************************************************************************/
bool gammalib::dir_exists(const std::string& dirname)
{
    // Initialise result
    bool result = false;

    // Allocate file information structure
    struct stat info;

    // Get file information structure
    int ret = stat(gammalib::expand_env(dirname).c_str(), &info);

    // Check if we have a directory
    if (ret == 0 && S_ISDIR(info.st_mode)) {
        result = true;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Checks if a substring is in a string
 *
 * @param[in] str string you want to search in.
 * @param[in] substring string you are looking for in str.
 *
 * Checks if substring is contained in str
 ***************************************************************************/
bool gammalib::contains(const std::string& str, const std::string& substring)
{
    // Initialise result
    bool result = false;

    // checks if substring is in str
    if (str.find(substring) != std::string::npos){
    	result = true;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Emits warning
 *
 * @param[in] origin Name of method that emits warning.
 * @param[in] message Warning message.
 *
 * Writes a warning to the console.
 ***************************************************************************/
void gammalib::warning(const std::string& origin,
                       const std::string& message)
{
    // Compile option: enable/disable warnings
    #if defined(G_WARNINGS)

    // Set warning message
    std::string warning = "+++ WARNING in " + origin + ": " + message;

    // Writes warning to the console
    std::cout << warning << std::endl;

    // End of compile option
    #endif
    
    // Return
    return;
}
