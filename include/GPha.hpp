/***************************************************************************
 *                GPha.hpp - XSPEC Pulse Height Analyzer class             *
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
 * @file GPha.hpp
 * @brief XSPEC Pulse Height Analyzer class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPHA_HPP
#define GPHA_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GEbounds.hpp"
#include "GFilename.hpp"
#include "GNdarray.hpp"
#include "GFitsHeader.hpp"

/* __ Forward declarations _______________________________________________ */
class GEnergy;
class GFits;
class GFitsTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_pha = "SPECTRUM";
}


/***********************************************************************//**
 * @class GPha
 *
 * @brief Pulse Height Analyzer class
 *
 * This class implements a Pulse Height Analyzer (PHA) spectrum that is used
 * as data container for an XSPEC analysis. A PHA spectrum is a vector that
 * provides the number of measured counts as function of the channel number.
 *
 * As an extension to the PHA format, the GPha class may also store the
 * energy boundaries for all PHA channels. If defined, the energy boundaries
 * will be written as a EBOUNDS extension to the same file where the PHA
 * spectrum resides. Upon loading, GPha will also load the energy boundaries
 * from an EBOUNDS extension if they are present.
 ***************************************************************************/
class GPha : public GBase {

    // Operator friends
    friend GPha operator+(const GPha& a,       const GPha& b);
    friend GPha operator-(const GPha& a,       const GPha& b);
    friend GPha operator*(const GPha& pha,     const double& scale);
    friend GPha operator*(const double& scale, const GPha& pha);
    friend GPha operator/(const GPha& pha,     const double& scale);

public:
    // Constructors and destructors
    GPha(void);
    explicit GPha(const GFilename& filename);
    explicit GPha(const GEbounds& ebds);
    explicit GPha(const int& bins);
    GPha(const GPha& pha);
    virtual ~GPha(void);

    // Operators
    GPha&         operator=(const GPha& pha);
    GPha&         operator+=(const GPha& pha);
    GPha&         operator-=(const GPha& pha);
    GPha&         operator*=(const double& scale);
    GPha&         operator/=(const double& scale);
    double&       operator[](const int& index);
    const double& operator[](const int& index) const;
    double&       operator()(const int& index, const int& col);
    const double& operator()(const int& index, const int& col) const;

    // Additional column access operators
    std::vector<double>&       operator[](const std::string& colname);
    const std::vector<double>& operator[](const std::string& colname) const;

    // Methods
    void               clear(void);
    GPha*              clone(void) const;
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
    double             counts(void) const;
    GNdarray           counts_spectrum(void) const;
    void               areascal(const int& index, const double& areascal);
    const double&      areascal(const int& index) const;
    void               backscal(const int& index, const double& backscal);
    const double&      backscal(const int& index) const;
    GNdarray           backscal_spectrum(void) const;
    const double&      underflow(void) const;
    const double&      overflow(void) const;
    const double&      outflow(void) const;
    void               exposure(const double& exposure);
    const double&      exposure(void) const;
    void               emin_obs(const GEnergy& emin_obs);
    const GEnergy&     emin_obs(void) const;
    void               emax_obs(const GEnergy& emax_obs);
    const GEnergy&     emax_obs(void) const;
    void               backfile(const std::string& backfile);
    const std::string& backfile(void) const;
    void               corrfile(const std::string& corrfile);
    const std::string& corrfile(void) const;
    void               respfile(const std::string& respfile);
    const std::string& respfile(void) const;
    void               ancrfile(const std::string& ancrfile);
    const std::string& ancrfile(void) const;
    void               fill(const GEnergy& energy, const double& value = 1.0);
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
    void copy_members(const GPha& pha);
    void free_members(void);
    void alloc(const int& size);
    int  column_index(const std::string& colname) const;

    // Protected members
    mutable GFilename                 m_filename;  //!< Filename of origin
    std::vector<double>               m_counts;    //!< Counts data
    std::vector<double>               m_areascal;  //!< Area scaling
    std::vector<double>               m_backscal;  //!< Background scaling
    std::vector<std::string>          m_colnames;  //!< Additional column names
    std::vector<std::vector<double> > m_coldata;   //!< Additional column data
    GEnergy                           m_emin_obs;  //!< Minimum energy of observation
    GEnergy                           m_emax_obs;  //!< Minimum energy of observation
    double                            m_underflow; //!< Number of underflowing events
    double                            m_overflow;  //!< Number of overflowing events
    double                            m_outflow;   //!< Number of outflowing events
    double                            m_exposure;  //!< Deadtime corr. exp. time (sec)
    GEbounds                          m_ebounds;   //!< Energy boundaries
    std::string                       m_backfile;  //!< Background file
    std::string                       m_corrfile;  //!< Correction file
    std::string                       m_respfile;  //!< RMF file
    std::string                       m_ancrfile;  //!< ARF file
    GFitsHeader                       m_header;    //!< FITS header cards
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GPha").
 ***************************************************************************/
inline
std::string GPha::classname(void) const
{
    return ("GPha");
}


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
 * @brief Return content of additional columns
 *
 * @param[in] index Bin index [0,...,size()-1].
 * @param[in] col Columns index [0,...,columns()-1].
 *
 * Returns reference to content of additional columns.
 ***************************************************************************/
inline
double& GPha::operator()(const int& index, const int& col)
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
const double& GPha::operator()(const int& index, const int& col) const
{
    return (m_coldata[col][index]);
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
    return (int)m_counts.size();
}


/***********************************************************************//**
 * @brief Return number of additional columns
 *
 * @return Number of additional columns.
 *
 * Returns the number of additional columns.
 ***************************************************************************/
inline
int GPha::columns(void) const
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
const GEbounds& GPha::ebounds(void) const
{
    return (m_ebounds);
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
    return (m_underflow);
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
    return (m_overflow);
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
    return (m_outflow);
}


/***********************************************************************//**
 * @brief Set exposure time
 *
 * @param[in] exposure Exposure time (seconds).
 *
 * Set the exposure time in seconds.
 ***************************************************************************/
inline
void GPha::exposure(const double& exposure)
{
    m_exposure = exposure;
    return;
}


/***********************************************************************//**
 * @brief Return exposure time
 *
 * @return Exposure time (seconds).
 *
 * Returns the exposure time in seconds.
 ***************************************************************************/
inline
const double& GPha::exposure(void) const
{
    return (m_exposure);
}


/***********************************************************************//**
 * @brief Set minimum energy of observations
 *
 * @param[in] emin_obs Minimum energy of observation.
 *
 * Set the minimum energy of the observation.
 ***************************************************************************/
inline
void GPha::emin_obs(const GEnergy& emin_obs)
{
    m_emin_obs = emin_obs;
    return;
}


/***********************************************************************//**
 * @brief Return minimum energy of observations
 *
 * @return Minimum energy of observation.
 *
 * Returns the minimum energy of the observation.
 ***************************************************************************/
inline
const GEnergy& GPha::emin_obs(void) const
{
    return (m_emin_obs);
}


/***********************************************************************//**
 * @brief Set maximum energy of observations
 *
 * @param[in] emax_obs Maximum energy of observation.
 *
 * Set the maximum energy of the observation.
 ***************************************************************************/
inline
void GPha::emax_obs(const GEnergy& emax_obs)
{
    m_emax_obs = emax_obs;
    return;
}


/***********************************************************************//**
 * @brief Return maximum energy of observations
 *
 * @return Maximum energy of observation.
 *
 * Returns the maximum energy of the observation.
 ***************************************************************************/
inline
const GEnergy& GPha::emax_obs(void) const
{
    return (m_emax_obs);
}


/***********************************************************************//**
 * @brief Set background file name
 *
 * @param[in] backfile Background file name.
 *
 * Set the background file name.
 ***************************************************************************/
inline
void GPha::backfile(const std::string& backfile)
{
    m_backfile = backfile;
    return;
}


/***********************************************************************//**
 * @brief Return background file name
 *
 * @return Background file name.
 *
 * Returns the background file name.
 ***************************************************************************/
inline
const std::string& GPha::backfile(void) const
{
    return (m_backfile);
}


/***********************************************************************//**
 * @brief Set correction file name
 *
 * @param[in] corrfile Correction file name.
 *
 * Set the correction file name.
 ***************************************************************************/
inline
void GPha::corrfile(const std::string& corrfile)
{
    m_corrfile = corrfile;
    return;
}


/***********************************************************************//**
 * @brief Return correction file name
 *
 * @return Correction file name.
 *
 * Returns the correction file name.
 ***************************************************************************/
inline
const std::string& GPha::corrfile(void) const
{
    return (m_corrfile);
}


/***********************************************************************//**
 * @brief Set Redistribution Matrix File name
 *
 * @param[in] respfile Redistribution Matrix File name.
 *
 * Set the Redistribution Matrix File name.
 ***************************************************************************/
inline
void GPha::respfile(const std::string& respfile)
{
    m_respfile = respfile;
    return;
}


/***********************************************************************//**
 * @brief Return Redistribution Matrix File name
 *
 * @return Redistribution Matrix File name.
 *
 * Returns the Redistribution Matrix File name.
 ***************************************************************************/
inline
const std::string& GPha::respfile(void) const
{
    return (m_respfile);
}


/***********************************************************************//**
 * @brief Set Ancilliary Response File name
 *
 * @param[in] ancrfile Ancilliary Response File name.
 *
 * Set the Ancilliary Response File name.
 ***************************************************************************/
inline
void GPha::ancrfile(const std::string& ancrfile)
{
    m_ancrfile = ancrfile;
    return;
}


/***********************************************************************//**
 * @brief Return Ancilliary Response File name
 *
 * @return Ancilliary Response File name.
 *
 * Returns the Ancilliary Response File name.
 ***************************************************************************/
inline
const std::string& GPha::ancrfile(void) const
{
    return (m_ancrfile);
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
const GFilename& GPha::filename(void) const
{
    return (m_filename);
}


/***********************************************************************//**
 * @brief Return FITS header
 *
 * @return FITS header or PHA file.
 ***************************************************************************/
inline
const GFitsHeader& GPha::header(void) const
{
    return (m_header);
}


/***********************************************************************//**
 * @brief Set FITS header
 *
 * @param[in] header FITS header.
 ***************************************************************************/
inline
void GPha::header(const GFitsHeader& header)
{
    m_header = header;
    return;
}


/***********************************************************************//**
 * @brief Spectrum addition operator friend
 *
 * @param[in] a First Pulse Height Analyzer spectrum.
 * @param[in] b Second Pulse Height Analyzer spectrum.
 * @return Sum of Pulse Height Analyzer spectra.
 ***************************************************************************/
inline
GPha operator+(const GPha& a, const GPha& b)
{
    GPha result = a;
    result     += b;
    return result;
}


/***********************************************************************//**
 * @brief Spectrum subtraction operator friend
 *
 * @param[in] a First Pulse Height Analyzer spectrum.
 * @param[in] b Second Pulse Height Analyzer spectrum.
 * @return Difference of Pulse Height Analyzer spectra.
 ***************************************************************************/
inline
GPha operator-(const GPha& a, const GPha& b)
{
    GPha result = a;
    result     -= b;
    return result;
}


/***********************************************************************//**
 * @brief Spectrum scaling operator friend
 *
 * @param[in] pha Pulse Height Analyzer spectrum.
 * @param[in] scale Scale factor.
 * @return Scaled Pulse Height Analyzer spectrum.
 ***************************************************************************/
inline
GPha operator*(const GPha& pha, const double& scale)
{
    GPha result = pha;
    result     *= scale;
    return result;
}


/***********************************************************************//**
 * @brief Spectrum scaling operator friend
 *
 * @param[in] scale Scale factor.
 * @param[in] pha Pulse Height Analyzer spectrum.
 * @return Scaled Pulse Height Analyzer spectrum.
 ***************************************************************************/
inline
GPha operator*(const double& scale, const GPha& pha)
{
    GPha result = pha;
    result     *= scale;
    return result;
}


/***********************************************************************//**
 * @brief Spectrum division operator friend
 *
 * @param[in] pha Pulse Height Analyzer spectrum.
 * @param[in] scale Division factor.
 * @return Divided Pulse Height Analyzer spectrum.
 ***************************************************************************/
inline
GPha operator/(const GPha& pha, const double& scale)
{
    GPha result = pha;
    result     /= scale;
    return result;
}

#endif /* GPHA_HPP */
