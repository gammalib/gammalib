/***************************************************************************
 *              GArf.hpp - XSPEC Auxiliary Response File class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2018 by Juergen Knoedlseder                         *
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
#include <vector>
#include "GBase.hpp"
#include "GEbounds.hpp"
#include "GFilename.hpp"
#include "GNodeArray.hpp"
#include "GFitsHeader.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFitsTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_arf = "SPECRESP";
}


/***********************************************************************//**
 * @class GArf
 *
 * @brief Auxiliary Response File class
 ***************************************************************************/
class GArf : public GBase {

    // Operator friends
    friend GArf operator+(const GArf& a,       const GArf& b);
    friend GArf operator-(const GArf& a,       const GArf& b);
    friend GArf operator*(const GArf& arf,     const double& scale);
    friend GArf operator*(const double& scale, const GArf& arf);
    friend GArf operator/(const GArf& arf,     const double& scale);

public:
    // Constructors and destructors
    GArf(void);
    explicit GArf(const GFilename& filename);
    explicit GArf(const GEbounds& ebds);
    GArf(const GArf& arf);
    virtual ~GArf(void);

    // Operators
    GArf&         operator=(const GArf& arf);
    GArf&         operator+=(const GArf& arf);
    GArf&         operator-=(const GArf& arf);
    GArf&         operator*=(const double& scale);
    GArf&         operator/=(const double& scale);
    double&       operator[](const int& index);
    const double& operator[](const int& index) const;
    double&       operator()(const int& index, const int& col);
    const double& operator()(const int& index, const int& col) const;

    // Additional column access operators
    std::vector<double>&       operator[](const std::string& colname);
    const std::vector<double>& operator[](const std::string& colname) const;
    double                     operator()(const std::string& colname,
                                          const GEnergy&     energy) const;

    // Methods
    void               clear(void);
    GArf*              clone(void) const;
    std::string        classname(void) const;
    int                size(void) const;
    int                columns(void) const;
    double&            at(const int& index);
    const double&      at(const int& index) const;
    double&            at(const int& index, const int& col);
    const double&      at(const int& index, const int& col) const;
    void               append(const std::string&         name,
                              const std::vector<double>& column);
    const GEbounds&    ebounds(void) const;
    void               load(const GFilename& filename);
    void               save(const GFilename& filename,
                            const bool&      clobber = false) const;
    void               read(const GFits& fits);
    void               read(const GFitsTable& table);
    void               write(GFits& fits) const;
    const GFilename&   filename(void) const;
    const GFitsHeader& header(void) const;
    void               header(const GFitsHeader& header);
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GArf& pha);
    void free_members(void);
    void set_logetrue(void);
    int  column_index(const std::string& colname) const;
    
    // Protected members
    mutable GFilename                 m_filename; //!< Filename of origin
    GEbounds                          m_ebounds;  //!< Energy boundaries
    GNodeArray                        m_logetrue; //!< Log10 energies in TeV
    std::vector<double>               m_specresp; //!< Spectral response
    std::vector<std::string>          m_colnames; //!< Additional column names
    std::vector<std::vector<double> > m_coldata;  //!< Additional column data
    GFitsHeader                       m_header;   //!< FITS header cards
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GArf").
 ***************************************************************************/
inline
std::string GArf::classname(void) const
{
    return ("GArf");
}


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
 * @brief Return content of additional columns
 *
 * @param[in] index Bin index [0,...,size()-1].
 * @param[in] col Columns index [0,...,columns()-1].
 *
 * Returns reference to content of additional columns.
 ***************************************************************************/
inline
double& GArf::operator()(const int& index, const int& col)
{
    return (m_coldata[col][index]);
}


/***********************************************************************//**
 * @brief Return content of additional columns (const version)
 *
 * @param[in] index Bin index [0,...,size()-1].
 * @param[in] col Columns index [0,...,columns()-1].
 *
 * Returns reference to content of additional columns.
 ***************************************************************************/
inline
const double& GArf::operator()(const int& index, const int& col) const
{
    return (m_coldata[col][index]);
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
    return (int)m_specresp.size();
}


/***********************************************************************//**
 * @brief Return number of additional columns
 *
 * @return Number of additional columns.
 *
 * Returns the number of additional columns.
 ***************************************************************************/
inline
int GArf::columns(void) const
{
    return (int)m_colnames.size();
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
const GFilename& GArf::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Return FITS header
 *
 * @return FITS header or ARF file.
 ***************************************************************************/
inline
const GFitsHeader& GArf::header(void) const
{
    return (m_header);
}


/***********************************************************************//**
 * @brief Set FITS header
 *
 * @param[in] header FITS header.
 ***************************************************************************/
inline
void GArf::header(const GFitsHeader& header)
{
    m_header = header;
    return;
}


/***********************************************************************//**
 * @brief Auxiliary Response File addition operator friend
 *
 * @param[in] a First Auxiliary Response File.
 * @param[in] b Second Auxiliary Response File.
 * @return Sum of Auxiliary Response Files.
 ***************************************************************************/
inline
GArf operator+(const GArf& a, const GArf& b)
{
    GArf result = a;
    result     += b;
    return result;
}


/***********************************************************************//**
 * @brief Auxiliary Response File subtraction operator friend
 *
 * @param[in] a First Auxiliary Response File.
 * @param[in] b Second Auxiliary Response File.
 * @return Difference of Auxiliary Response Files.
 ***************************************************************************/
inline
GArf operator-(const GArf& a, const GArf& b)
{
    GArf result = a;
    result     -= b;
    return result;
}


/***********************************************************************//**
 * @brief Auxiliary Response File scaling operator friend
 *
 * @param[in] arf Auxiliary Response File.
 * @param[in] scale Scale factor.
 * @return Scaled Auxiliary Response File.
 ***************************************************************************/
inline
GArf operator*(const GArf& arf, const double& scale)
{
    GArf result = arf;
    result     *= scale;
    return result;
}


/***********************************************************************//**
 * @brief Auxiliary Response File scaling operator friend
 *
 * @param[in] scale Scale factor.
 * @param[in] arf Auxiliary Response File.
 * @return Scaled Auxiliary Response File.
 ***************************************************************************/
inline
GArf operator*(const double& scale, const GArf& arf)
{
    GArf result = arf;
    result     *= scale;
    return result;
}


/***********************************************************************//**
 * @brief Auxiliary Response File vision operator friend
 *
 * @param[in] arf Auxiliary Response File.
 * @param[in] scale Division factor.
 * @return Divided Auxiliary Response File.
 ***************************************************************************/
inline
GArf operator/(const GArf& arf, const double& scale)
{
    GArf result = arf;
    result     /= scale;
    return result;
}

#endif /* GARF_HPP */
