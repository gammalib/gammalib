/***************************************************************************
 *                        GEnergy.hpp - Energy class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GEnergy.hpp
 * @brief Energy value class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GENERGY_HPP
#define GENERGY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GEnergy
 *
 * @brief Class that handles energies in a unit independent way.
 *
 * The GEnergy class stores an energy value in units of MeV and implements
 * methods that provide automatic conversion of the energy values in other
 * units. This makes instrument specific implementations more robust and
 * reduces the risk of unit errors. The method also holds on request the
 * log10 of the energy in MeV, which is only recomputed when the energy
 * value is changed. This considerably reduces the number of log10 operations
 * in many applications and speeds up the performance.
 ***************************************************************************/
class GEnergy : public GBase {

    // Operator friends
    friend GEnergy operator+(const GEnergy &a, const GEnergy &b);
    friend GEnergy operator-(const GEnergy &a, const GEnergy &b);
    friend GEnergy operator*(const double &a, const GEnergy &b);
    friend GEnergy operator*(const GEnergy &a, const double &b);
    friend GEnergy operator/(const GEnergy &a, const double &b);
    friend bool    operator==(const GEnergy &a, const GEnergy &b);
    friend bool    operator!=(const GEnergy &a, const GEnergy &b);
    friend bool    operator<(const GEnergy &a, const GEnergy &b);
    friend bool    operator<=(const GEnergy &a, const GEnergy &b);
    friend bool    operator>(const GEnergy &a, const GEnergy &b);
    friend bool    operator>=(const GEnergy &a, const GEnergy &b);

public:
    // Constructors and destructors
    GEnergy(void);
    GEnergy(const GEnergy& eng);
    GEnergy(const double& eng, const std::string& unit);
    virtual ~GEnergy(void);
 
    // Operators
    GEnergy& operator=(const GEnergy& eng);
    GEnergy& operator+=(const GEnergy& eng);
    GEnergy& operator-=(const GEnergy& eng);
    GEnergy& operator*=(const double& scale);
    GEnergy& operator/=(const double& scale);
    void     operator()(const double& eng, const std::string& unit);
    double   operator()(const std::string& unit) const;

    // Methods
    void        clear(void);
    GEnergy*    clone(void) const;
    std::string classname(void) const;
    double      erg(void) const;
    double      keV(void) const;
    double      MeV(void) const;
    double      GeV(void) const;
    double      TeV(void) const;
    double      Angstrom(void) const;
    double      log10erg(void) const;
    double      log10keV(void) const;
    double      log10MeV(void) const;
    double      log10GeV(void) const;
    double      log10TeV(void) const;
    double      log10(const std::string& unit) const;
    void        erg(const double& eng);
    void        keV(const double& eng);
    void        MeV(const double& eng);
    void        GeV(const double& eng);
    void        TeV(const double& eng);
    void        Angstrom(const double& wavelength);
    void        log10erg(const double& eng);
    void        log10keV(const double& eng);
    void        log10MeV(const double& eng);
    void        log10GeV(const double& eng);
    void        log10TeV(const double& eng);
    void        log10(const double& eng, const std::string& unit);
    std::string print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEnergy& eng);
    void free_members(void);

    // Protected data members
    double         m_energy;      //!< Energy in MeV
    mutable double m_elog10;      //!< log10 of energy in MeV
    mutable bool   m_has_log10;   //!< log10 of energy is valid
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GEnergy").
 ***************************************************************************/
inline
std::string GEnergy::classname(void) const
{
    return ("GEnergy");
}


/***********************************************************************//**
 * @brief Energy unary addition operator
 *
 * @param[in] eng Energy.
 * @return Sum of energies.
 ***************************************************************************/
inline
GEnergy& GEnergy::operator+=(const GEnergy& eng)
{
    m_energy    += eng.m_energy;
    m_has_log10  = false;
    return *this;
}


/***********************************************************************//**
 * @brief Energy unary differnce operator
 *
 * @param[in] eng Energy.
 * @return Difference of energies.
 ***************************************************************************/
inline
GEnergy& GEnergy::operator-=(const GEnergy& eng)
{
    m_energy    -= eng.m_energy;
    m_has_log10  = false;
    return *this;
}


/***********************************************************************//**
 * @brief Energy unary multiplication operator
 *
 * @param[in] scale Scale.
 * @return Energy multiplied with scale.
 ***************************************************************************/
inline
GEnergy& GEnergy::operator*=(const double& scale)
{
    m_energy    *= scale;
    m_has_log10  = false;
    return *this;
}


/***********************************************************************//**
 * @brief Energy unary division operator
 *
 * @param[in] scale Scale.
 * @return Energy divided by scale.
 ***************************************************************************/
inline
GEnergy& GEnergy::operator/=(const double& scale)
{
    m_energy    /= scale;
    m_has_log10  = false;
    return *this;
}


/***********************************************************************//**
 * @brief Energy addition operator friend
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return Sum of energies.
 ***************************************************************************/
inline
GEnergy operator+(const GEnergy& a, const GEnergy& b)
{
    GEnergy result;
    result.m_energy = a.m_energy + b.m_energy;
    return result;
}


/***********************************************************************//**
 * @brief Energy subtraction operator friend
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return Difference of energies.
 ***************************************************************************/
inline
GEnergy operator-(const GEnergy& a, const GEnergy& b)
{
    GEnergy result;
    result.m_energy = a.m_energy - b.m_energy;
    return result;
}


/***********************************************************************//**
 * @brief Energy multiplication operator friend
 *
 * @param[in] a Scale.
 * @param[in] b Energy.
 * @return Energy multiplied by scale.
 ***************************************************************************/
inline
GEnergy operator*(const double& a, const GEnergy& b)
{
    GEnergy result;
    result.m_energy = a * b.m_energy;
    return result;
}


/***********************************************************************//**
 * @brief Energy multiplication operator friend
 *
 * @param[in] a Energy.
 * @param[in] b Scale.
 * @return Energy multiplied by scale.
 ***************************************************************************/
inline
GEnergy operator*(const GEnergy& a, const double& b)
{
    GEnergy result;
    result.m_energy = b * a.m_energy;
    return result;
}


/***********************************************************************//**
 * @brief Energy division operator friend
 *
 * @param[in] a Energy.
 * @param[in] b Scale.
 * @return Energy divided by scale.
 ***************************************************************************/
inline
GEnergy operator/(const GEnergy& a, const double& b)
{
    GEnergy result;
    result.m_energy = a.m_energy / b;
    return result;
}


/***********************************************************************//**
 * @brief Energy equality operator friend
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return True if both energies are equal.
 ***************************************************************************/
inline
bool operator==(const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy == b.m_energy);
}


/***********************************************************************//**
 * @brief Energy non-equality operator friend
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return True if both energies are not equal.
 ***************************************************************************/
inline
bool operator!=(const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy != b.m_energy);
}


/***********************************************************************//**
 * @brief Energy smaller than operator friend
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return True if first energy is smaller than second energy.
 ***************************************************************************/
inline
bool operator<(const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy < b.m_energy);
}


/***********************************************************************//**
 * @brief Energy smaller than or equal to operator friend
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return True if first energy is smaller than or equal to second energy.
 ***************************************************************************/
inline
bool operator<=(const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy <= b.m_energy);
}


/***********************************************************************//**
 * @brief Energy larger than operator friend
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return True if first energy is larger than second energy.
 ***************************************************************************/
inline
bool operator>(const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy > b.m_energy);
}


/***********************************************************************//**
 * @brief Energy larger than or equal to operator friend
 *
 * @param[in] a First energy.
 * @param[in] b Second energy.
 * @return True if first energy is larger than or equal to second energy.
 ***************************************************************************/
inline
bool operator>=(const GEnergy &a, const GEnergy &b)
{
    return (a.m_energy >= b.m_energy);
}

#endif /* GENERGY_HPP */
