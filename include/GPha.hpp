/***************************************************************************
 *                GPha.hpp - XSPEC Pulse Height Analyzer class             *
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
 * @file GPha.hpp
 * @brief XSPEC Pulse Height Analyzer class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPHA_HPP
#define GPHA_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GEnergy.hpp"
#include "GEbounds.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GPha
 *
 * @brief Pulse Height Analyzer class
 ***************************************************************************/
class GPha : public GBase {

public:
    // Constructors and destructors
    GPha(void);
    explicit GPha(const std::string& filename);
    explicit GPha(const GEbounds& ebds);
    GPha(const GPha& pha);
    virtual ~GPha(void);

    // Operators
    GPha&         operator=(const GPha& pha);
    double&       operator[](const int& index);
    const double& operator[](const int& index) const;

    // Methods
    void               clear(void);
    GPha*              clone(void) const;
    int                size(void) const;
    double&            at(const int& index);
    const double&      at(const int& index) const;
    const GEbounds&    ebounds(void) const;
    double             counts(void) const;
    const double&      underflow(void) const;
    const double&      overflow(void) const;
    const double&      outflow(void) const;
    void               fill(const GEnergy& energy, const double& value = 1.0);
    void               load(const std::string& filename);
    void               save(const std::string& filename, const bool& clobber = false) const;
    void               read(const GFitsTable* hdu);
    void               write(GFits& fits) const;
    const std::string& filename(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPha& pha);
    void free_members(void);
    
    // Protected members
    mutable std::string m_filename;   //!< Filename of origin
    GEbounds            m_ebounds;    //!< Energy boundaries
    std::vector<double> m_counts;     //!< Counts data
    double              m_underflow;  //!< Number of underflowing events
    double              m_overflow;   //!< Number of overflowing events
    double              m_outflow;    //!< Number of outflowing events
};


/***********************************************************************//**
 * @brief Return content of spectral bin
 *
 * @param[in] index Bin index [0,...,size()-1].
 *
 * Returns reference to content of spectral bin with specified @p index.
 ***************************************************************************/
inline
double& GPha::operator[](const int& index)
{
    return (m_counts[index]);
}


/***********************************************************************//**
 * @brief Return content of spectral bin (const version)
 *
 * @param[in] index Bin index [0,...,size()-1].
 *
 * Returns reference to content of spectral bin with specified @p index.
 ***************************************************************************/
inline
const double& GPha::operator[](const int& index) const
{
    return (m_counts[index]);
}


/***********************************************************************//**
 * @brief Return number of bins in spectrum
 *
 * @return Number of bins in spectrum.
 *
 * Returns the number of bins in the spectrum.
 ***************************************************************************/
inline
int GPha::size(void) const
{
    return (m_counts.size());
}


/***********************************************************************//**
 * @brief Return energy boundaries
 *
 * @return Energy boundaries for all spectral bins.
 *
 * Returns the energy boundaries for all spectral bins.
 ***************************************************************************/
inline
const GEbounds& GPha::ebounds(void) const
{
    return m_ebounds;
}


/***********************************************************************//**
 * @brief Return number of underflow counts
 *
 * @return Number of counts with energies lower than the first energy bin.
 *
 * Returns the number of counts with energies lower than the first energy
 * bin.
 ***************************************************************************/
inline
const double& GPha::underflow(void) const
{
    return m_underflow;
}


/***********************************************************************//**
 * @brief Return number of overflow counts
 *
 * @return Number of counts with energies larger than the last energy bin.
 *
 * Returns the number of counts with energies larger than the last energy
 * bin.
 ***************************************************************************/
inline
const double& GPha::overflow(void) const
{
    return m_overflow;
}


/***********************************************************************//**
 * @brief Return number of outflow counts
 *
 * @return Number of counts with energies between energy bins.
 *
 * Returns the number of counts with energies between energy bins.
 ***************************************************************************/
inline
const double& GPha::outflow(void) const
{
    return m_outflow;
}


/***********************************************************************//**
 * @brief Return file name
 *
 * @return File name from which the PHA information has been read or into
 *         which PHA information has been saved.
 *
 * Returns the file name from which the PHA information has been read or into
 * which PHA information has been saved. The returned string will be empty if
 * no load() or save() method has been called before.
 ***************************************************************************/
inline
const std::string& GPha::filename(void) const
{
    return (m_filename);
}

#endif /* GPHA_HPP */
