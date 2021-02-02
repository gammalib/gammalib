/***************************************************************************
 *              GCTABackground2D.hpp - CTA 2D background class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GCTABackground2D.hpp
 * @brief CTA 2D background class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTABACKGROUND2D_HPP
#define GCTABACKGROUND2D_HPP

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
    const std::string extname_cta_background2d = "BKG";
}


/***********************************************************************//**
 * @class GCTABackground2D
 *
 * @brief CTA 2D background class
 ***************************************************************************/
class GCTABackground2D : public GCTABackground {

public:
    // Constructors and destructors
    GCTABackground2D(void);
    explicit GCTABackground2D(const GFilename& filename);
    GCTABackground2D(const GCTABackground2D& bgd);
    virtual ~GCTABackground2D(void);

    // Implemented pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& detx, 
                              const double& dety) const;

    // Operators
    GCTABackground2D& operator=(const GCTABackground2D& bgd);

    // Implemented pure virtual methods
    void                       clear(void);
    GCTABackground2D*          clone(void) const;
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
    void   copy_members(const GCTABackground2D& bgd);
    void   free_members(void);
    void   set_members(void);
    int    index(const int& itheta, const int& iebin) const;
    void   init_mc_cache(void) const;
    void   init_mc_max_rate(void) const;
    double solid_angle(const double& detx1, const double& dety1,
                       const double& detx2, const double& dety2,
                       const double& detx3, const double& dety3) const;
    double rate(const int& iebin, const double& theta) const;

    // Members
    GFilename         m_filename;    //!< Name of background response file
    GCTAResponseTable m_background;  //!< Background response table
    GEnergies         m_energy;      //!< Vector of energies
    int               m_inx_theta;   //!< THETA axis index
    int               m_inx_energy;  //!< Energy axis index
    int               m_inx_bgd;     //!< Background index
    int               m_num_theta;   //!< Number of THETA bins
    int               m_num_energy;  //!< Number of energy bins
    int               m_num[2];      //!< Array of number of bins
    double            m_theta_min;   //!< THETA minimum (radians)
    double            m_theta_max;   //!< THETA maximum (radians)
    double            m_logE_min;    //!< Log10(E/TeV) minimum
    double            m_logE_max;    //!< Log10(E/TeV) maximum

    // Monte Carlo cache
    mutable std::vector<double> m_mc_max;      //!< Maximum background rate
    mutable GModelSpectralNodes m_mc_spectrum; //!< Response cube spectrum
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTABackground2D").
 ***************************************************************************/
inline
std::string GCTABackground2D::classname(void) const
{
    return ("GCTABackground2D");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the background was loaded.
 ***************************************************************************/
inline
GFilename GCTABackground2D::filename(void) const
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
const GModelSpectralNodes& GCTABackground2D::spectrum(void) const
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
bool GCTABackground2D::is_valid(void) const
{
    return (m_background.axes() == 2);
}


/***********************************************************************//**
 * @brief Return response table
 *
 * @return Response table.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTABackground2D::table(void) const
{
    return m_background;
}


/***********************************************************************//**
 * @brief Assign response table
 *
 * @param[in] table Response table.
 ***************************************************************************/
inline
void GCTABackground2D::table(const GCTAResponseTable& table)
{
     m_background = table;
     set_members();
}

#endif /* GCTABACKGROUND2D_HPP */
