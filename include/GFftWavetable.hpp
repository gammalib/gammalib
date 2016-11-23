/***************************************************************************
 *  GFftWavetable.hpp - Lookup table class for Fast Fourier transformation *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GFftWavetable.hpp
 * @brief Lookup table class interface definition for Fast Fourier transformation
 * @author Juergen Knoedlseder
 */

#ifndef GFFTWAVETABLE_HPP
#define GFFTWAVETABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <complex>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GFftWavetable
 *
 * @brief Lookup table class for Fast Fourier Transformation
 ***************************************************************************/
class GFftWavetable : public GBase {

public:
    // Constructors and destructors
    GFftWavetable(void);
    explicit GFftWavetable(const int& size);
    GFftWavetable(const GFftWavetable& wavetable);
    virtual ~GFftWavetable(void);

    // Operators
    GFftWavetable&              operator=(const GFftWavetable& wavetable);
    std::complex<double>&       operator[](const int& index);
    const std::complex<double>& operator[](const int& index) const;

    // Methods
    void           clear(void);
    GFftWavetable* clone(void) const;
    std::string    classname(void) const;
    int            size(void) const;
    const int&     index(const int& factor) const;
    int            factors(void) const;
    int            factor(const int& index) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GFftWavetable& wavetable);
    void free_members(void);
    void set_members(const int& n);
    void set_factors(const int& n);

    // Protected members
    std::vector<int>                   m_factors; //!< Wavetable factors
    std::vector<int>                   m_twiddle; //!< Start index of factors
    std::vector<std::complex<double> > m_trig;    //!< Trigonometric coefficients
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFftWavetable").
 ***************************************************************************/
inline
std::string GFftWavetable::classname(void) const
{
    return ("GFftWavetable");
}


/***********************************************************************//**
 * @brief Return reference to trigonometric coefficient
 *
 * @param[in] index Trigonometric coefficient index [0,...,size()-1].
 *
 * Returns a reference to the trigonometric coefficient with the specified
 * @p index.
 ***************************************************************************/
inline
std::complex<double>& GFftWavetable::operator[](const int& index)
{
    return (m_trig[index]);
}


/***********************************************************************//**
 * @brief Return reference to trigonometric coefficient (const version)
 *
 * @param[in] index Trigonometric coefficient index [0,...,size()-1].
 *
 * Returns a const reference to the trigonometric coefficient with the
 * specified @p index.
 ***************************************************************************/
inline
const std::complex<double>& GFftWavetable::operator[](const int& index) const
{
    return (m_trig[index]);
}


/***********************************************************************//**
 * @brief Return start index for a given factor
 *
 * @param[in] factor Factor index [0,...,factors()-1].
 *
 * Returns the index of the first trigonometric coefficient for a given
 * @p factor.
 ***************************************************************************/
inline
const int& GFftWavetable::index(const int& factor) const
{
    return (m_twiddle[factor]);
}


/***********************************************************************//**
 * @brief Return number of trigonometric coefficients
 *
 * @return Number of trigonometric coefficients.
 ***************************************************************************/
inline
int GFftWavetable::size(void) const
{
    return (m_trig.size());
}


/***********************************************************************//**
 * @brief Return number of factorisation factors
 *
 * @return Number of factorisation factors.
 ***************************************************************************/
inline
int GFftWavetable::factors(void) const
{
    return (m_factors.size());
}

#endif /* GFFTWAVETABLE_HPP */
