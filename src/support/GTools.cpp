/***************************************************************************
 *                       GTools.cpp - GammaLib tools                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
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
#include <unistd.h>        // access() functions
#include <fcntl.h>         // fcntl() functions
#include <sys/stat.h>
#include <sys/time.h>      // timeval
#include <sys/select.h>    // select() function
#include <sys/types.h>
#include <sys/socket.h>    // recv() function
#include <pwd.h>
#include <cmath>
#include <cfloat>
#include <cctype>
#include <cstdlib>         // std::getenv() function
#include <cstring>         // std::strlen() function
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdio>          // sprintf
#include "GTools.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GEnergy.hpp"
#include "GFilename.hpp"
#include "GMath.hpp"
#include "GXmlElement.hpp"
#include "GXmlNode.hpp"
#include "GXmlText.hpp"

/* __ Compile options ____________________________________________________ */
//#define G_CHECK_FOR_NAN

/* __ Function name definitions __________________________________________ */
#define G_XML2STRING                     "gammalib::xml2string(std::string&)"

/* __ Coding definitions _________________________________________________ */
#define G_PARFORMAT_LENGTH 29
#define G_CHAR_BUFFER      256


/***********************************************************************//**
 * @brief Strip leading and trailing whitespace from string
 *
 * @param[in] arg String from which whitespace should be stripped.
 * @return String with stripped whitespace.
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
 * @return String with stripped characters.
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

        // Continue only if start and stop are valid
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
 * @brief Strip trailing character from string
 *
 * @param[in] arg String from which character should be stripped.
 * @param[in] chars Character(s) to be stripped.
 * @return String with stripped characters.
 ***************************************************************************/
std::string gammalib::rstrip_chars(const std::string& arg,
                                   const std::string& chars)
{
    // Initialise empty result string
    std::string result;

    // Continue only if argument is not empty
    if (!arg.empty()) {

        // Get stop
        std::string::size_type stop = arg.find_last_not_of(chars);

        // Continue only if stop is valid
        if (stop != std::string::npos) {
            result = arg.substr(0, stop+1);
        }

    } // endif: argument was not empty

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Replace string segment in string
 *
 * @param[in] arg String in which character segements are to be replaced
 * @param[in] segment String segment to be replaced.
 * @param[in] replacement Replacement string.
 * @return String with replaced segements.
 *
 * Replaces string segments by a replacement string in a given string.
 *
 * If the input string @p arg is "Wonderful", the @p segment is "onder" and
 * the @p replacement is "ish" the method will return "Wishful".
 ***************************************************************************/
std::string gammalib::replace_segment(const std::string& arg,
                                      const std::string& segment,
                                      const std::string& replacement)
{
    // Initialise result string by argument
    std::string result = arg;

    // Initialise character pointer
    std::string::size_type pos   = 0;
    std::string::size_type start = 0;

    // Loop over string
    while (start != std::string::npos) {
        start = result.find(segment, pos);
        if (start != std::string::npos) {
            result = result.replace(start, segment.length(), replacement);
            pos    = start + replacement.length();
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Expand environment variables in string
 *
 * @param[in] arg String.
 * @return String with expected environment variables.
 *
 * Expands any environment variable that is found in a string. Valid
 * delimiters for environment variables are $ENV{\<name\>}, $ENV(\<name\>),
 * ${\<name\>}, $(\<name\>) and $\<name\> (in the last case the terminating
 * delimiter is either a / or a blank character or the end of the string).
 * Environment variables occuring within single quotes (') are ignored.
 * Environment variables that do not exist will be kept as specified.
 *
 * The method also replaces ~ or ~user by the user's home directory, ~+
 * by the value of the PWD environment variable and ~- by the value of the
 * OLDPWD variable. If the user or the PWD or OLDPWD variable are not found,
 * no replacement is done.
 *
 * This function has been inspired by the function ape_util_expand_env_var
 * from ape_util.c in the ape software developed at HEASARC.
 ***************************************************************************/
std::string gammalib::expand_env(const std::string& arg)
{
    // Set environment variable delimiters
    static const char* begin_delim[] = { "$ENV{", "$ENV(", "${", "$(", "$", "~" };
    static const char* end_delim[]   = { "}", ")", "}", ")", "/", "/" };
    static const int   num_delim     = 6;

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
            // variable. There is a special case for delimiters 4 and 5:
            // It has always an end_length of zero as the / is not a real
            // delimiter, but just an indicator that the environment variable
            // ends. Another indicator is a blank. If the end of the
            // string has been reached this is also acceptable.
            size_t i_start    = index + begin_length;
            size_t i_end      = i_start;
            size_t end_length = 0;
            if (delim_idx == 4 || delim_idx == 5) {
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
            if (i_end < result.length() || delim_idx == 4 || delim_idx == 5) {

                // Extract environment variable name
                std::string name        = result.substr(i_start, i_end-i_start);
                size_t      name_length = name.length();

                // Initialise pointer on environment variable
                const char* env = NULL;

                // Handle ~
                if (delim_idx == 5) {
                    if (name_length == 0) {
                        if ((env = std::getenv("HOME")) == NULL) {
                            struct passwd *pw = getpwuid(getuid());
                            if (pw != NULL) {
                                env = pw->pw_dir;
                            }
                        }
                    }
                    else {
                        if (name == "+") {
                            env = std::getenv("PWD");
                        }
                        else if (name == "-") {
                            env = std::getenv("OLDPWD");
                        }
                        else {
                            struct passwd *pw = getpwnam(name.c_str());
                            if (pw != NULL) {
                                env = pw->pw_dir;
                            }
                        }
                    }
                }

                // ... otherwise get the environment variable
                else {
                    env = std::getenv(name.c_str());
                }

                // If the environment variable has been found then replace
                // it by its value
                if (env != NULL) {

                    // Erase delimiters and environment variable
                    result.erase(index, begin_length+name_length+end_length);

                    // Set replacement string and its length
                    std::string replace(env);
                    size_t      replace_length = replace.length();

                    // Insert replacement string
                    result.insert(index, replace);

                    // Advance pointer
                    index += replace_length;

                } // endif: environment variable has been found

                // If no environment variable has been found then set
                // index=i_end+1
                else {
                    index = i_end + 1;
                }

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
 * @brief Return temporary file name
 *
 * @return Temporary file name.
 *
 * Returns a temporary file name.
 ***************************************************************************/
std::string gammalib::tmpnam(void)
{
    // Set default temporary directory
    char default_tmpdir[] = "/tmp";

    // Get pointer to name of temporary directory by searching various
    // possible environment variables
    char *tmpdir = NULL;
    if ((tmpdir = std::getenv("TEMP")) == NULL) {
        if ((tmpdir = std::getenv("TMP")) == NULL) {
            if ((tmpdir = std::getenv("TMPDIR")) == NULL) {
                tmpdir = default_tmpdir;
            }
        }
    }

    // Allocate empty filename
    char filename[256];
    filename[0] = '\0';

    // Combine temporary directory with temporary filename
    strcat(filename, tmpdir);
    strcat(filename, "/gammalibXXXXXX");

    // Create temporary file
    int         fd = mkstemp(filename);
    std::string tmpname(filename);
    close(fd);
    unlink(filename);

    // Return temporary file name
    return tmpname;
}


/***********************************************************************//**
 * @brief Convert unsigned short integer value into string
 *
 * @param[in] value Unsigned short integer to be converted into string.
 * @return String with unsigned short integer value.
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
 * @return String with unsigned integer value.
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
 * @return String with unsigned long integer value.
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
 * @return String with unsigned long long integer value.
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
 * @return String with short integer value.
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
 * @param[in] fmt Format string.
 * @return String with integer value.
 ***************************************************************************/
std::string gammalib::str(const int& value, const std::string& fmt)
{
    // Allocate character buffer
    char buffer[G_CHAR_BUFFER];

    // Put integer into buffer
    sprintf(buffer, fmt.c_str(), value);

    // Convert buffer into string
    std::string str_buffer(buffer);

    // Return string buffer
    return str_buffer;
}


/***********************************************************************//**
 * @brief Convert long integer value into string
 *
 * @param[in] value Long integer to be converted into string.
 * @return String with long integer value.
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
 * @return String with long long integer value.
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
 * @param[in] precision Floating point precision.
 * @return String with single precision value.
 *
 * Converts a single precision value into a string. Any positive
 * @p precision argument specifies the exact number of digits after the
 * comma.
 ***************************************************************************/
std::string gammalib::str(const float& value, const int& precision)
{
    // Allocate output stream
    std::ostringstream s_value;

    // If specified then set the requested fixed point precision. Otherwise
    // use a precision that should be sufficient for floating point values.
    if (precision > 0) {
        s_value.precision(precision);
        s_value.setf(std::ios::fixed, std::ios::floatfield);
    }
    else {
        s_value.precision(7);
    }

    // Put floating point value in stream
    s_value << value;

    // Convert to a string
    std::string result = s_value.str();

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Convert double precision value into string
 *
 * @param[in] value Double precision value to be converted into string.
 * @param[in] precision Floating point precision.
 * @return String with double precision value.
 *
 * Converts a double precision value into a string. Any positive
 * @p precision argument specifies the exact number of digits after the
 * comma.
 ***************************************************************************/
std::string gammalib::str(const double& value, const int& precision)
{
    // Allocate output stream
    std::ostringstream s_value;

    // If specified then set the requested fixed point precision. Otherwise
    // use a precision that should be sufficient for floating point values.
    if (precision > 0) {
        s_value.precision(precision);
        s_value.setf(std::ios::fixed, std::ios::floatfield);
    }
    else {
        s_value.precision(15);
    }

    // Put double precision floating point value in stream
    s_value << value;

    // Convert to a string
    std::string result = s_value.str();

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Convert complex value into string
 *
 * @param[in] value Complex value to be converted into string.
 * @param[in] precision Floating point precision.
 * @return String with complex value.
 *
 * Converts a complex value into a string. Any positive @p precision argument
 * specifies the exact number of digits after the comma.
 ***************************************************************************/
std::string gammalib::str(const std::complex<double>& value,
                          const int&                  precision)
{
    // Allocate output stream
    std::ostringstream s_value;

    // If specified then set the requested fixed point precision. Otherwise
    // use a precision that should be sufficient for floating point values.
    if (precision > 0) {
        s_value.precision(precision);
        s_value.setf(std::ios::fixed, std::ios::floatfield);
    }
    else {
        s_value.precision(15);
    }

    // Put double precision floating point value in stream
    s_value << value.real();
    if (value.imag() < 0.0) {
        s_value << "-";
    }
    else {
        s_value << "+";
    }
    s_value << std::abs(value.imag()) << "j";

    // Convert to a string
    std::string result = s_value.str();

    // Return result
    return result;
}



/***********************************************************************//**
 * @brief Convert string to C string
 *
 * @param[in] arg String to be converted.
 * @return C string.
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
 * @return Short value.
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
 * @return Unsigned short value.
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
 * @return Integer value.
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
 * @return Unsigned integer value.
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
 * @return Long value.
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
 * @return Unsigned long value.
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
 * @return Long long value.
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
 * @return Unsigned long long value.
 ***************************************************************************/
unsigned long long gammalib::toulonglong(const std::string& arg)
{
    std::istringstream iss(arg);
    unsigned long long result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into single precision value
 *
 * @param[in] arg String to be converted.
 * @return Single precision value.
 ***************************************************************************/
float gammalib::tofloat(const std::string& arg)
{
    std::istringstream iss(arg);
    float              result;
    iss >> std::dec >> result;
    return result;
}


/***********************************************************************//**
 * @brief Convert string into double precision value
 *
 * @param[in] arg String to be converted.
 * @return Double precision value.
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
 * @return String converted to upper case.
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
 * @return String converted to lower case.
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
 * @return Vector of split strings.
 *
 * Splits a string on the basis of one or multiple separator characters. The
 * separator characters are provided by the @p sep argument. Subsequent
 * separator characters that are not seperated by some other characters will
 * lead to an empty string element, except for a blank separator where
 * subsequent blanks are takens as a single separator. Below a few examples
 * that illustrate how the function will split a given string.
 *
 *     "Name;RA;DEC" => ["Name","RA","DEC"] (sep=";")
 *     "My house is    red" => ["My","house","is","red"] (sep=" ")
 *     "IRF::FRONT" => ["IRF","","FRONT"] (sep=":")
 *     "Fields;RA,DEC,Flux" => ["Fields","RA","DEC","Flux"] (sep=";,")
 *     "Last;Field;" => ["Last","Field",""] (sep=";")
 ***************************************************************************/
std::vector<std::string> gammalib::split(const std::string& s,
                                         const std::string& sep)
{
    // Allocate result string vector
    std::vector<std::string> result;

    // Initialise counters
    std::size_t pos = 0;
    std::size_t len = s.length();

    // Loop over string
    while (pos < len && pos != std::string::npos) {

        // Get index of first separator occurence and preset the length
        // of the substring to the end of the string
        std::size_t index = s.find_first_of(sep, pos);
        std::size_t n     = std::string::npos;

        // If we did not reach the end then compute now the length of the
        // substring
        if (index != std::string::npos) {
            n = index-pos;
        }

        // If we have no whitespace separator and the length of the
        // substring is zero then push back an empty string. If the
        // length of the substring is positive then push back the
        // substring.
        if (sep != " " && n == 0) {
            result.push_back("");
        }
        else if (n > 0) {
            result.push_back(s.substr(pos, n));
        }

        // Go to the string position after the last separator
        pos = (index != std::string::npos) ? index + 1 : std::string::npos;

        // If the position is pointing right beyong the last string
        // character we terminated with a separator, hence we need to push
        // back one more empty string before we leave
        if (sep != " " && pos == len) {
            result.push_back("");
        }

    } // endwhile: there were still characters in the string

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Fill string with n strings of same type
 *
 * @param[in] s String to be filled.
 * @param[in] n Number of fillings.
 * @return Filled strings.
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
 * @return Left justified string.
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
 * @return Right justified string.
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
 * @return Centred string.
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
 * @param[in] indent Indentation of parameter.
 * @return Parameter string.
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
 * @brief Convert singular noun into number noun
 *
 * @param[in] noun Singular noun.
 * @param[in] number Number of instance of noun.
 * @return Converted noun.
 *
 * Converts a singular noun into a number noun by appending a "s" to the
 * noun if the @p number of the instances of the noun is not one.
 ***************************************************************************/
std::string gammalib::number(const std::string& noun, const int& number)
{
    // Copy input noun
    std::string result(noun);

    // Append "s" if number if not one
    if (number != 1) {
        result += "s";
    }

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
 * @return Photon flux under power law.
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
        double exponent = gamma + 1.0;
        if (std::abs(exponent) > 1.0e-11) {
            double xmin = emin   / epivot;
            double xmax = emax   / epivot;
            flux        = epivot / exponent * (std::pow(xmax, exponent) -
                                               std::pow(xmin, exponent));
        }
        else {
            double ratio = emax / emin;
            flux         = epivot * std::log(ratio);
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
 * @return Energy flux under power law.
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
        double exponent = gamma + 2.0;
        double epivot2  = epivot * epivot;
        if (std::abs(exponent) > 1.0e-11) {
            double xmin = emin    / epivot;
            double xmax = emax    / epivot;
            flux        = epivot2 / exponent * (std::pow(xmax, exponent) -
                                                std::pow(xmin, exponent));
        }
        else {
            double ratio = emax / emin;
            flux         = epivot2 * (std::log(ratio));
        }

    } // endif: emax > emin

    // Return result
    return flux;
}


/***********************************************************************//**
 * @brief Computes log mean energy
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return Log mean energy.
 *
 * Computes the logarithmic mean energy
 * \f$10^{0.5 * (\log E_{\rm a} + \log E_{\rm b})}\f$
 * for two energies.
 ***************************************************************************/
GEnergy gammalib::elogmean(const GEnergy& a, const GEnergy& b)
{
    // Compute logarithmic mean energy
    GEnergy elogmean;
    double  eloga = a.log10MeV();
    double  elogb = b.log10MeV();
    elogmean.MeV(std::pow(10.0, 0.5 * (eloga + elogb)));

    // Return
    return elogmean;
}


/***********************************************************************//**
 * @brief Checks if directory exists
 *
 * @param[in] dirname Directory name.
 * @return True if directory exists, false otherwise.
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
 * @param[in] str String you want to search in.
 * @param[in] substring String you are looking for in @p str.
 * @return True if a string contains a sub string.
 *
 * Checks if substring is contained in str
 ***************************************************************************/
bool gammalib::contains(const std::string& str, const std::string& substring)
{
    // Initialise result
    bool result = false;

    // checks if substring is in str
    if (str.find(substring) != std::string::npos) {
        result = true;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Checks if a string is contained in a vector of strings
 *
 * @param[in] strings Vector of strings you want to search in.
 * @param[in] string string you are looking for in strings.
 * @return True if a string is contained a vector.
 *
 * Checks if a string is contained in a vector of strings
 ***************************************************************************/
bool gammalib::contains(const std::vector<std::string>& strings,
                        const std::string&              string)
{
    // Compute result
    bool result = std::find(strings.begin(), strings.end(), string) !=
                  strings.end();

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


/***********************************************************************//**
 * @brief Convert XML character references in string to characters
 *
 * @param[in] arg String containing XML character references.
 * @return String with character reference replaced by respective characters.
 *
 * Converts all character references found in a string in their respective
 * characters. For more information about XML character references read
 * http://en.wikipedia.org/wiki/List_of_XML_and_HTML_character_entity_references
 ***************************************************************************/
std::string gammalib::xml2str(const std::string& arg)
{
    // Initialise string
    std::string result;

    // Iitialise position, lenghts and flags
    size_t length = arg.length();
    size_t pos    = 0;
    size_t start  = 0;
    size_t stop   = 0;
    size_t len    = 0;      // Length of string preceeding char. reference
    bool   found  = false;

    // Loop over string
    while (pos < length) {

        // If we have not yet found a character reference then search for
        // the next one. If we do not find one we break here.
        if (!found) {
            start = arg.find("&", pos);
            if (start != std::string::npos) {
                len   = start - pos;
                pos   = start;
                found = true;
            }
            else {
                break;
            }
        }

        // ... otherwise search for end of the actual character reference.
        // Throw an exception if no end is found.
        else {
            stop = arg.find(";", pos);
            if (stop != std::string::npos) {

                // First attach leading string to the result
                if (len > 0) {
                    result += arg.substr(start-len, len);
                }

                // Next extract character reference
                std::string cref = arg.substr(start+1, stop-start-1);
                len              = cref.length();

                // Check for a numerical character reference
                //TODO: Check that there are only valid characters in
                //      numerical field
                if (len >= 2 && cref[0] == '#') {
                    int number = -1;
                    if (cref[1] == 'x') {
                        number = (int)std::strtol(cref.substr(2,len-2).c_str(), NULL, 16);
                    }
                    else {
                        number = toint(cref.substr(1,len-1));
                    }
                    if (number != -1) {
                        result.push_back((char)number);
                    }
                    else {
                        std::string msg = "Could not extract number from "
                                          "numerical character reference &"+
                                          cref+";";
                        throw GException::invalid_argument(G_XML2STRING, msg);
                    }
                }

                // ... otherwise check for a character entity reference
                // and push back the corresponding character to the result
                // string
                else {
                    if (cref == "quot") {
                        result.push_back((char)34);
                    }
                    else if (cref == "amp") {
                        result.push_back((char)38);
                    }
                    else if (cref == "apos") {
                        result.push_back((char)39);
                    }
                    else if (cref == "lt") {
                        result.push_back((char)60);
                    }
                    else if (cref == "gt") {
                        result.push_back((char)62);
                    }
                    else {
                        std::string msg = "Unknown character entity reference "
                                          "&"+cref+"; encountered in XML string \""+
                                          arg+"\".";
                        throw GException::invalid_argument(G_XML2STRING, msg);
                    }
                }

                // Signal that we're done and that we search for the
                // next character reference
                found = false;
                pos   = stop + 1;

            } // endif: end of character reference found

            // ... otherwise throw an exception
            else {
                std::string msg = "Missing ; character at end of character "
                                  "reference in XML string \""+arg+"\".";
                throw GException::invalid_argument(G_XML2STRING, msg);
            }

        } // endelse: search for end of character reference

    } // endwhile

    // Append any pending string to the result
    if (pos < length) {
        len     = length - pos;
        result += arg.substr(pos, len);
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Convert special characters in string to XML character references
 *
 * @param[in] arg String.
 * @return String with special characters replaced by character references.
 *
 * Converts all special characters found in a string into character
 * references. For more information about XML character references read
 * http://en.wikipedia.org/wiki/List_of_XML_and_HTML_character_entity_references
 ***************************************************************************/
std::string gammalib::str2xml(const std::string& arg)
{
    // Initialise string
    std::string result;

    // Loop over string
    for (int i = 0; i < arg.length(); ++i) {

        // Initialize string to add
        std::string character = arg.substr(i, 1);

        // Replace special characters
        if (character == "\"") {
            character = "&quot;";
        }
        else if (character == "&") {
            character = "&amp;";
        }
        else if (character == "'") {
            character = "&apos;";
        }
        else if (character == "<") {
            character = "&lt;";
        }
        else if (character == ">") {
            character = "&gt;";
        }

        // Append character
        result += character;

    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Checks if parameter with given name in XML element exists
 *
 * @param[in] xml XML element.
 * @param[in] name Parameter name.
 * @return True if parameter exists, false otherwise.
 *
 * Checks whether a parameter with given @p name exists in XML element.
 ***************************************************************************/
bool gammalib::xml_has_par(const GXmlElement& xml, const std::string& name)
{
    // Initialize flag
    bool found = false;

    // Get number of elements
    int n = xml.elements("parameter");
    
    // Search for parameter with given name
    for (int i = 0; i < n; ++i) {
        const GXmlElement* element = xml.element("parameter", i);
        if (element->attribute("name") == name) {
            found = true;
            break;
        }
    }

    // Return
    return found;
}


/***********************************************************************//**
 * @brief Return pointer to parameter with given name in XML element
 *
 * @param[in] origin Method requesting parameter.
 * @param[in] xml XML element.
 * @param[in] name Parameter name.
 * @return Pointer to parameter XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Returns pointer to parameter with given @p name in XML element. If the
 * @p name is not found, a parameter with the given @p name is added. In
 * that respect the function differs from xml_get_par which does not add a
 * parameter element.
 *
 * The function checks for multiple occurences of a parameter and throws an
 * exception in case that more than one parameter with a given name is found.
 ***************************************************************************/
GXmlElement* gammalib::xml_need_par(const std::string& origin,
                                    GXmlElement&       xml,
                                    const std::string& name)
{
    // Initialize XML element pointer
    GXmlElement* par = NULL;

    // Number of elements
    int number = 0;

    // Get number of elements in XML element
    int n = xml.elements("parameter");

    // Search for parameter with given name
    for (int i = 0; i < n; ++i) {
        GXmlElement* element = xml.element("parameter", i);
        if (element->attribute("name") == name) {
            par = element;
            number++;
        }
    }

    // Create parameter if none was found
    if (number == 0) {
        par = static_cast<GXmlElement*>(xml.append(GXmlElement("parameter name=\""+name+"\"")));
        number++;
    }

    // Check that there are no multiple parameters
    gammalib::xml_check_par(origin, name, number);

    // Return
    return par;
}


/***********************************************************************//**
 * @brief Return pointer to parameter with given name in XML element
 *
 * @param[in] origin Method requesting parameter.
 * @param[in] xml XML element.
 * @param[in] name Parameter name.
 * @return Pointer to parameter XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Returns pointer to parameter with given @p name in XML element. The
 * function checks whether the parameter has been found and throws an
 * exception if no parameter or multiple occurences of a parameter with given
 * @p name are found.
 ***************************************************************************/
const GXmlElement* gammalib::xml_get_par(const std::string& origin,
                                         const GXmlElement& xml,
                                         const std::string& name)
{
    // Initialize XML element pointer
    const GXmlElement* par = NULL;

    // Number of elements
    int number = 0;

    // Get number of elements in XML element
    int n = xml.elements("parameter");

    // Search for parameter with given name
    for (int i = 0; i < n; ++i) {
        const GXmlElement* element = xml.element("parameter", i);
        if (element->attribute("name") == name) {
            par = element;
            number++;
        }
    }

    // Check that there is at least one parameter and that there are no
    // multiple parameters
    gammalib::xml_check_par(origin, name, number);

    // Return
    return par;
}


/***********************************************************************//**
 * @brief Return attribute value for a given parameter in XML element
 *
 * @param[in] origin Method requesting parameter.
 * @param[in] xml XML element.
 * @param[in] name Parameter name.
 * @param[in] attribute Attribute name.
 * @return Value of attribute.
 *
 * @exception GException::invalid_value
 *            Attribute not found.
 *
 * Returns the value of @p attribute of parameter @p name in XML element.
 * The function checks whether the parameter has been found and throws an
 * exception if no parameter or multiple occurences of a parameter with given
 * @p name are found. The function furthermore checks whether the attribute
 * exists. 
 ***************************************************************************/
std::string gammalib::xml_get_attr(const std::string& origin,
                                   const GXmlElement& xml,
                                   const std::string& name,
                                   const std::string& attribute)
{
    // Initialise attribute value
    std::string value = "";

    // Get parameter
    const GXmlElement* par = gammalib::xml_get_par(origin, xml, name);

    // Throw an exception if a parameter has not the requested attribute
    if (!par->has_attribute(attribute)) {
        std::string msg = "Attribute \""+attribute+"\" not found in XML "
                          "parameter \""+name+"\". Please verify the XML "
                          "format.";
        throw GException::invalid_value(origin, msg);
    }

    // Extract attribute
    value = par->attribute(attribute);

    // Return attribute value
    return value;
}


/***********************************************************************//**
 * @brief Checks whether a parameter has occured once
 *
 * @param[in] origin Method performing the check.
 * @param[in] name Parameter name.
 * @param[in] number Number of occurences of parameter.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Throws an exception if a given parameter has not exactly occured once.
 * The exception text is adapted to the case that none or multiple parameters
 * have been found.
 ***************************************************************************/
void gammalib::xml_check_par(const std::string& origin,
                             const std::string& name,
                             const int&         number)
{
    // Throw case dependent exception
    if (number < 1) {
        std::string msg = "Parameter \""+name+"\" not found in XML element."
                          " Please verify the XML format.";
        throw GException::invalid_value(origin, msg);
    }
    else if (number > 1) {
        std::string msg = "Parameter \""+name+"\" found "+
                          gammalib::str(number)+" times in XML element."
                          " Please verify the XML format.";
        throw GException::invalid_value(origin, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Expand file name provided as XML attribute for loading
 *
 * @param[in] xml XML element.
 * @param[in] filename File name.
 * @return Expanded file name.
 *
 * Expands file name provided as XML attribute for loading. If the file name
 * is not empty and has no path it is assumed that the file is located in the
 * same directory as the XML file, and the XML file access path is prepended
 * to the file name.
 ***************************************************************************/
GFilename gammalib::xml_file_expand(const GXmlElement& xml,
                                    const std::string& filename)
{
    // Set file name
    GFilename fname(filename);

    // If the file name is not empty and has no path we assume that the file
    // name is a relative file name with respect to the XML file access path
    // and we therefore prepend the XML file access path to the file name
    if (!fname.is_empty() && fname.path().length() == 0) {
        fname = xml.filename().path() + fname;
    }

    // Return file name
    return fname;
}


/***********************************************************************//**
 * @brief Reduce file name provided for writing as XML attribute
 *
 * @param[in] xml XML element.
 * @param[in] filename File name.
 * @return Reduced file name.
 *
 * Reduces file name provided for writing as XML attribute. If the file name
 * is not empty and has the same access path as the XML file it is assumed
 * that both files are located in the same directory, and the access path is
 * stripped from the file name.
 ***************************************************************************/
GFilename gammalib::xml_file_reduce(const GXmlElement& xml,
                                    const std::string& filename)
{
    // Set file name
    GFilename fname(filename);

    // If the file name is not empty and has the same access path as the XML
    // file it is assumed that both files are located in the same directory,
    // and the access path is stripped from the file name.
    if (!fname.is_empty() && fname.path() == xml.filename().path()) {
        fname = fname.file();
    }

    // Return file name
    return fname;
}


/***********************************************************************//**
 * @brief Extract name / value pair from XML node
 *
 * @param[in] node Pointer to XML node.
 * @param[out] name Name string.
 * @param[out] value Value string.
 *
 * Extracts a name / value pair from a XML node. If the XML node pointer is
 * NULL, the name and value strings will be empty.
 ***************************************************************************/
void gammalib::xml_get_name_value_pair(const GXmlNode* node,
                                       std::string&    name,
                                       std::string&    value)
{
    // Clear name and value strings
    name.clear();
    value.clear();

    // Continue only if node is valid
    if (node != NULL) {

        // Get name node and extract text content
        const GXmlNode* ptr = node->element("name", 0);
        if (ptr != NULL) {
            const GXmlText* text = static_cast<const GXmlText*>((*ptr)[0]);
            if (text != NULL) {
                name = text->text();
            }
        }

        // Get value node and extract text content
        ptr = node->element("value", 0);
        if (ptr != NULL) {
            const GXmlText* text = static_cast<const GXmlText*>((*ptr)[0]);
            if (text != NULL) {
                value = text->text();
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks whether a parameter has occured once
 *
 * @param[in] fd Socket file descriptor.
 * @param[out] buffer Buffer to hold data.
 * @param[in] len Maximum number of bytes to recv().
 * @param[in] flags Flags (as the fourth param to recv() ).
 * @param[in] timeout Timeout in milliseconds.
 * @return recv() error code, -2 == timeout
 *
 * This function implements the recv() function with a timeout. The timeout
 * is specified in milliseconds.
 ***************************************************************************/
int gammalib::recv(int fd, char *buffer, int len, int flags, int timeout)
{
    // Initialise error code with time out
    int error = -2;

    // Initialize the set
    fd_set readset;
    FD_ZERO(&readset);
    FD_SET(fd, &readset);

    // Initialize time out struct
    struct timeval tv;
    if (timeout >= 1000) {
        tv.tv_sec  = timeout/1000;
        tv.tv_usec = 0;
    }
    else {
        tv.tv_sec  = 0;
        tv.tv_usec = timeout*1000;
    }

    // select()
    int result = select(fd+1, &readset, NULL, NULL, &tv);

    // Check status
    if (result < 0) {
        error = -1;
    }
    else if (result > 0 && FD_ISSET(fd, &readset)) {

        // Initialise flags
        int iof = -1;

        // Set non-blocking mode
        if ((iof = fcntl(fd, F_GETFL, 0)) != -1) {
            fcntl(fd, F_SETFL, iof | O_NONBLOCK);
        }

        // Receive data
        result = ::recv(fd, buffer, len, flags);

        // Set as before
        if (iof != -1) {
            fcntl(fd, F_SETFL, iof);
        }

        // Set error
        error = result;
    }

    // Return error
    return error;
}


/***********************************************************************//**
 * @brief Returns length of circular arc within circular ROI
 *
 * @param[in] rad Circle radius in radians (<pi).
 * @param[in] dist Circle centre distance to ROI centre (<pi).
 * @param[in] cosdist Cosine of circle centre distance to ROI centre.
 * @param[in] sindist Sinus of circle centre distance to ROI centre.
 * @param[in] roi Radius of ROI in radians.
 * @param[in] cosroi Cosine of ROI radius.
 *
 * This method returns the arclength in radians of a circle of radius 'rad'
 * with a centre that is offset by 'dist' from the ROI centre, where the ROI
 * radius is given by 'roi'. To speed-up computations, the cosines and sinus
 * of 'roi' and 'psf' should be calculated by the client and be passed to
 * the method.
 ***************************************************************************/
double gammalib::roi_arclength(const double& rad,     const double& dist,
                               const double& cosdist, const double& sindist,
                               const double& roi,     const double& cosroi)
{
    // Declare arclength
    double arclength;

    // Handle special case that circle centre matches ROI centre
    if (dist == 0.0) {
        if (rad > roi) arclength = 0.0;             // Circle outside ROI
        else           arclength = gammalib::twopi; // Circle inside ROI
    }

    // ... otherwise circle and ROI centres are not identical
    else {

        // Handle special case that we evaluate exactly at the circle
        // centre. In this case we have in fact a point, and if this point
        // falls within the ROI it has a formal arclength of 2pi.
        //
        if (rad == 0.0) {
            if (dist > roi) arclength = 0.0;             // Circle centre outside ROI
            else            arclength = gammalib::twopi; // Circle centre inside ROI
        }

        // ... otherwise we have to handle the general case
        else {
            double d = roi - dist;
            if (-rad >= d) {
                arclength = 0.0;
            }
            else if (rad <= d) {
                arclength = gammalib::twopi;
            }
            else {
                double cosrad = std::cos(rad);
                double sinrad = std::sin(rad);
                double cosang = (cosroi - cosdist*cosrad) / (sindist*sinrad);
                arclength     = 2.0 * gammalib::acos(cosang);
                #if defined(G_CHECK_FOR_NAN)
                if (gammalib::is_infinite(arclength) ||
                    gammalib::is_notanumber(arclength)) {
                    std::cout << "roi_arclength: NaN occured";
                    std::cout << " rad=" << rad;
                    std::cout << " sinrad=" << sinrad;
                    std::cout << " cosrad=" << cosrad;
                    std::cout << " sindist=" << sindist;
                    std::cout << " cosdist=" << cosdist;
                    std::cout << " cosang=" << cosang;
                    std::cout << std::endl;
                }
                #endif
            }
        }

    } // endelse: Circle and ROI centres were not identical

    // Return arclength
    return arclength;
}
