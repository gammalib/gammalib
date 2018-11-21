/***************************************************************************
 *              GCTABackground3D.hpp - CTA 3D background class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2018 by Juergen Knoedlseder                         *
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
#include "GEnergies.hpp"
#include "GFilename.hpp"
#include "GModelSpectralNodes.hpp"
#include "GCTABackground.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFitsBinTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_cta_background3d = "BACKGROUND";
}


/***********************************************************************//**
 * @class GCTABackground3D
 *
 * @brief CTA 3D background class
 ***************************************************************************/
class GCTABackground3D : public GCTABackground {

public:
    // Constructors and destructors
    GCTABackground3D(void);
    explicit GCTABackground3D(const GFilename& filename);
    GCTABackground3D(const GCTABackground3D& bgd);
    virtual ~GCTABackground3D(void);

    // Implemented pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety) const;

    // Operators
    GCTABackground3D& operator=(const GCTABackground3D& bgd);

    // Implemented pure virtual methods
    void                       clear(void);
    GCTABackground3D*          clone(void) const;
    std::string                classname(void) const;
    void                       load(const GFilename& filename);
    GFilename                  filename(void) const;
    GCTAInstDir                mc(const GEnergy& energy,
                                  const GTime& time,
                                  GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;
    double                     rate_ebin(const GCTAInstDir& dir,
                                         const GEnergy&     emin,
                                         const GEnergy&     emax) const;
    std::string                print(const GChatter& chatter = NORMAL) const;

    // Methods
    bool                       is_valid(void) const;
    const GCTAResponseTable&   table(void) const;
    void                       table(const GCTAResponseTable& table);
    void                       read(const GFitsTable& table);
    void                       write(GFitsBinTable& table) const;
    void                       save(const GFilename& filename,
                                    const bool&      clobber = false) const;

private:
    // Methods
    void   init_members(void);
    void   copy_members(const GCTABackground3D& bgd);
    void   free_members(void);
    void   set_members(void);
    int    index(const int& idetx, const int& idety, const int& iebin) const;
    void   init_mc_cache(void) const;
    void   init_mc_max_rate(void) const;
    double solid_angle(const double& detx1, const double& dety1,
                       const double& detx2, const double& dety2,
                       const double& detx3, const double& dety3) const;
    double rate(const int& iebin, const double& detx, const double& dety) const;

    // Members
    GFilename         m_filename;    //!< Name of background response file
    GCTAResponseTable m_background;  //!< Background response table
    GEnergies         m_energy;      //!< Vector of energies
    int               m_inx_detx;    //!< DETX index
    int               m_inx_dety;    //!< DETY index
    int               m_inx_energy;  //!< Energy index
    int               m_inx_bgd;     //!< Background index
    int               m_num_detx;    //!< Number of DETX bins
    int               m_num_dety;    //!< Number of DETY bins
    int               m_num_energy;  //!< Number of energy bins
    int               m_num[3];      //!< Array of number of bins
    double            m_detx_min;    //!< DETX minimum (radians)
    double            m_detx_max;    //!< DETX maximum (radians)
    double            m_dety_min;    //!< DETY minimum (radians)
    double            m_dety_max;    //!< DETY maximum (radians)
    double            m_logE_min;    //!< Log10(E/TeV) minimum
    double            m_logE_max;    //!< Log10(E/TeV) maximum

    // Monte Carlo cache
    mutable std::vector<double> m_mc_max;      //!< Maximum background rate
    mutable GModelSpectralNodes m_mc_spectrum; //!< Response cube spectrum
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTABackground3D").
 ***************************************************************************/
inline
std::string GCTABackground3D::classname(void) const
{
    return ("GCTABackground3D");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the background was loaded.
 ***************************************************************************/
inline
GFilename GCTABackground3D::filename(void) const
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
 * @brief Return validity of background model
 *
 * @return True if background model is valid.
 ***************************************************************************/
inline
bool GCTABackground3D::is_valid(void) const
{
    return (m_background.axes() == 3);
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


/***********************************************************************//**
 * @brief Assign response table
 *
 * @param[in] table Response table.
 ***************************************************************************/
inline
void GCTABackground3D::table(const GCTAResponseTable& table)
{
     m_background = table;
     set_members();
}

#endif /* GCTABACKGROUND3D_HPP */
