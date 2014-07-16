/***************************************************************************
 *              GCTABackground3D.hpp - CTA 3D background class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTABackground3D.hpp
 * @brief CTA 3D background class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTABACKGROUND3D_HPP
#define GCTABACKGROUND3D_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GCTABackground.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFitsBinTable;


/***********************************************************************//**
 * @class GCTABackground3D
 *
 * @brief CTA 3D background class
 ***************************************************************************/
class GCTABackground3D : public GCTABackground {

public:
    // Constructors and destructors
    GCTABackground3D(void);
    explicit GCTABackground3D(const std::string& filename);
    GCTABackground3D(const GCTABackground3D& bgd);
    virtual ~GCTABackground3D(void);

    // Implemented pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety,
                              const bool&   etrue = false) const;

    // Operators
    GCTABackground3D& operator=(const GCTABackground3D& bgd);

    // Implemented pure virtual methods
    void                       clear(void);
    GCTABackground3D*          clone(void) const;
    void                       load(const std::string& filename);
    std::string                filename(void) const;
    GCTAInstDir                mc(const GEnergy& energy,
                                  const GTime& time,
                                  GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;
    std::string                print(const GChatter& chatter = NORMAL) const;
   
    // Methods
    const GCTAResponseTable&   table(void) const;
    void                       read(const GFits& file);
    void                       write(GFitsBinTable& hdu) const;
    void                       save(const std::string& filename,
                                    const bool& clobber = false) const;
    
private:
    // Methods
    void init_members(void);
    void copy_members(const GCTABackground3D& bgd);
    void free_members(void);
    void init_mc_cache(const bool& etrue = false) const;

    // Members
    std::string       m_filename;    //!< Name of background response file
    GCTAResponseTable m_background;  //!< Background response table
    double            m_mc_max_bin;  //!< Maximum spatial binsize for MC
    double            m_mc_max_logE; //!< Maximum log energy binsize for MC

    // Monte Carlo cache
    mutable std::vector<double> m_mc_cache;    //!< Monte Carlo cache
    mutable GModelSpectralNodes m_mc_spectrum; //!< Response cube spectrum
    mutable int                 m_mc_nx;       //!< DETX pixels for MC
    mutable int                 m_mc_ny;       //!< DETY pixels for MC
    mutable int                 m_mc_npix;     //!< DETX*DETY pixels for MC
    mutable int                 m_mc_nmaps;    //!< Number of maps for MC
    mutable double              m_mc_detx_min; //!< DETX minimum
    mutable double              m_mc_detx_max; //!< DETX maximum
    mutable double              m_mc_detx_bin; //!< DETX binsize for MC
    mutable double              m_mc_dety_min; //!< DETY minimum
    mutable double              m_mc_dety_max; //!< DETY maximum
    mutable double              m_mc_dety_bin; //!< DETY binsize for MC
    mutable double              m_mc_logE_min; //!< log10 energy minimum (TeV)
    mutable double              m_mc_logE_max; //!< log10 energy maximum (TeV)
    mutable double              m_mc_logE_bin; //!< log10 energy binsize (TeV)
};


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the background was loaded.
 ***************************************************************************/
inline
std::string GCTABackground3D::filename(void) const
{
    // Return filename
    return m_filename;
}


/***********************************************************************//**
 * @brief Get response cube spectrum
 *
 * @return Response cube spectrum.
 *
 * Returns the response cube spectrum.
 ***************************************************************************/
inline
const GModelSpectralNodes& GCTABackground3D::spectrum(void) const
{
    if (m_mc_spectrum.nodes() == 0) {
        init_mc_cache();
    }
    return (m_mc_spectrum);
}


/***********************************************************************//**
 * @brief Return response table
 *
 * @return Response table.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTABackground3D::table(void) const
{
    return m_background;
}

#endif /* GCTABACKGROUND3D_HPP */
