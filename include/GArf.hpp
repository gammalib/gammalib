/***************************************************************************
 *              GArf.hpp - XSPEC Auxiliary Response File class             *
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
 * @file GArf.hpp
 * @brief XSPEC Auxiliary Response File class definition
 * @author Juergen Knoedlseder
 */

#ifndef GARF_HPP
#define GARF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GEbounds.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GArf
 *
 * @brief Auxiliary Response File class
 ***************************************************************************/
class GArf : public GBase {

public:
    // Constructors and destructors
    GArf(void);
    explicit GArf(const std::string& filename);
    explicit GArf(const GEbounds& ebds);
    GArf(const GArf& arf);
    virtual ~GArf(void);

    // Operators
    GArf&         operator=(const GArf& arf);
    double&       operator[](const int& index);
    const double& operator[](const int& index) const;

    // Methods
    void               clear(void);
    GArf*              clone(void) const;
    int                size(void) const;
    double&            at(const int& index);
    const double&      at(const int& index) const;
    const GEbounds&    ebounds(void) const;
    void               load(const std::string& filename);
    void               save(const std::string& filename, const bool& clobber = false) const;
    void               read(const GFitsTable* hdu);
    void               write(GFits& fits) const;
    const std::string& filename(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GArf& pha);
    void free_members(void);
    
    // Protected members
    mutable std::string m_filename;   //!< Filename of origin
    GEbounds            m_ebounds;    //!< Energy boundaries
    std::vector<double> m_specresp;   //!< Spectral response
};


/***********************************************************************//**
 * @brief Return content of spectral bin
 *
 * @param[in] index Bin index [0,...,size()-1].
 *
 * Returns reference to content of spectral bin with specified @p index.
 ***************************************************************************/
inline
double& GArf::operator[](const int& index)
{
    return (m_specresp[index]);
}


/***********************************************************************//**
 * @brief Return content of spectral bin (const version)
 *
 * @param[in] index Bin index [0,...,size()-1].
 *
 * Returns reference to content of spectral bin with specified @p index.
 ***************************************************************************/
inline
const double& GArf::operator[](const int& index) const
{
    return (m_specresp[index]);
}


/***********************************************************************//**
 * @brief Return number of spectral bins
 *
 * @return Number of spectral bins.
 *
 * Returns the number of spectral bins.
 ***************************************************************************/
inline
int GArf::size(void) const
{
    return (m_specresp.size());
}


/***********************************************************************//**
 * @brief Return energy boundaries
 *
 * @return Energy boundaries for all spectral bins.
 *
 * Returns the energy boundaries for all spectral bins.
 ***************************************************************************/
inline
const GEbounds& GArf::ebounds(void) const
{
    return m_ebounds;
}


/***********************************************************************//**
 * @brief Return file name
 *
 * @return File name from which the ARF information has been read or into
 *         which ARF information has been saved.
 *
 * Returns the file name from which the ARF information has been read or into
 * which ARF information has been saved. The returned string will be empty if
 * no load() or save() method has been called before.
 ***************************************************************************/
inline
const std::string& GArf::filename(void) const
{
    return (m_filename);
}

#endif /* GARF_HPP */
