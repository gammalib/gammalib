/***************************************************************************
 *                       GTools.hpp - GammaLib tools                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GTools.hpp
 * @brief Gammalib tools definition
 * @author Juergen Knoedlseder
 */

#ifndef GTOOLS_HPP
#define GTOOLS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <cmath>
#include <cfloat>

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const double MeV2erg    =  1.6021765e-6;
    const double erg2MeV    =  624150.96;
    const double pc2cm      =  3.08568025e18;
    const double sec_in_day = 86400.0;
}

/* __ Prototypes ________________________________________________________ */
namespace gammalib {
    std::string              strip_whitespace(const std::string& arg);
    std::string              strip_chars(const std::string& arg,
                                         const std::string& chars);
    std::string              expand_env(const std::string& arg);
    std::string              str(const unsigned short int& value);
    std::string              str(const unsigned int& value);
    std::string              str(const unsigned long int& value);
    std::string              str(const unsigned long long int& value);
    std::string              str(const short int& value);
    std::string              str(const int& value);
    std::string              str(const long int& value);
    std::string              str(const long long int& value);
    std::string              str(const float& value);
    std::string              str(const double& value);
    char*                    tochar(const std::string& arg);
    short                    toshort(const std::string& arg);
    unsigned short           toushort(const std::string& arg);
    int                      toint(const std::string& arg);
    unsigned int             touint(const std::string& arg);
    long                     tolong(const std::string& arg);
    unsigned long            toulong(const std::string& arg);
    long long                tolonglong(const std::string& arg);
    unsigned long long       toulonglong(const std::string& arg);
    float                    tofloat(const std::string& arg);
    double                   todouble(const std::string& arg);
    std::string              toupper(const std::string& s);
    std::string              tolower(const std::string& s);
    std::vector<std::string> split(const std::string& s, const std::string& sep);
    std::string              fill(const std::string& s, int n);
    std::string              left(const std::string& s, int n, char c = ' ');
    std::string              right(const std::string& s, int n, char c = ' ');
    std::string              center(const std::string& s, int n, char c = ' ');
    std::string              parformat(const std::string& s, const int& indent = 0);
    double                   plaw_photon_flux(const double& emin,
                                              const double& emax,
                                              const double& epivot,
                                              const double& gamma);
    double                   plaw_energy_flux(const double& emin,
                                              const double& emax,
                                              const double& epivot,
                                              const double& gamma);
    bool                     file_exists(const std::string& filename);
    bool                     isinfinite(const double& x);
    bool                     isnotanumber(const double& x);
    bool 					 contains(std::string str1, std::string str2);
    void                     warning(const std::string& origin,
                                     const std::string& message);
}

/***********************************************************************//**
 * @brief Signal if argument is infinite
 *
 * @param[in] x Argument.
 * @return True if argument @p x is infinite, false otherwise.
 *
 * Signals if the argument @p x is infinite.
 *
 * This function has been copied from gnulib.
 ***************************************************************************/
inline
bool gammalib::isinfinite(const double& x)
{
  return (x < -DBL_MAX || x > DBL_MAX);
}


/***********************************************************************//**
 * @brief Signal if argument is not a number
 *
 * @param[in] x Argument.
 * @return True if argument @p x is not a number, false otherwise.
 *
 * Signals if the argument @p x is not a number.
 *
 * This function is a very simple kluge. It may not work on all systems.
 ***************************************************************************/
inline
bool gammalib::isnotanumber(const double& x)
{
  return (x != x);
}

#endif /* GTOOLS_HPP */
