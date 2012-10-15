/***************************************************************************
 *                 GLATAeff.hpp  -  Fermi/LAT effective area               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GLATAeff.hpp
 * @brief Fermi/LAT effective area class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATAEFF_HPP
#define GLATAEFF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GLATPointing.hpp"
#include "GLATResponseTable.hpp"
#include "GLATEfficiency.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GLATAeff
 *
 * @brief Interface for the Fermi/LAT effective area
 *
 * This class handles the effective area information for Fermi/LAT. It also
 * handles the IRF efficiency information that has been introduced for
 * Pass 7 data analysis.
 ***************************************************************************/
class GLATAeff : public GBase {

public:
    // Constructors and destructors
    GLATAeff(void);
    GLATAeff(const std::string& filename);
    GLATAeff(const GLATAeff& aeff);
    virtual ~GLATAeff(void);

    // Operators
    GLATAeff& operator= (const GLATAeff& aeff);
    double    operator() (const double& logE, const double& ctheta);
    double    operator() (const double& logE, const double& ctheta, const double& phi);
    double    operator() (const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime, const GLATPointing& pnt);

    // Methods
    void         clear(void);
    GLATAeff*    clone(void) const;
    void         load(const std::string& filename);
    void         save(const std::string& filename, bool clobber = false);
    void         read(const GFits* file);
    void         write(GFits& file) const;
    int          size(void) const { return nenergies()*ncostheta(); }
    int          nenergies(void) const { return m_aeff_bins.nenergies(); }
    int          ncostheta(void) const { return m_aeff_bins.ncostheta(); }
    double       costhetamin(void) const { return m_min_ctheta; }
    void         costhetamin(const double& ctheta);
    bool         hasphi(void) const { return false; }
    bool         hasefficiency(void) const;
    double       efficiency_factor1(const GEnergy& srcEng) const;
    double       efficiency_factor2(const GEnergy& srcEng) const;
    std::string  print(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATAeff& aeff);
    void free_members(void);
    void read_aeff(const GFitsTable* hdu);
    void read_efficiency(const GFitsTable* hdu);
    void write_aeff(GFits& file) const;
    void write_efficiency(GFits& file) const;
    
    // Protected members
    GLATResponseTable   m_aeff_bins;    //!< Aeff energy and cos theta binning
    std::vector<double> m_aeff;         //!< Aeff array
    double              m_min_ctheta;   //!< Minimum valid cos(theta)
    bool                m_front;        //!< Response is for front section
    bool                m_back;         //!< Response is for back section
    GLATEfficiency*     m_eff_func1;    //!< Efficiency functor 1
    GLATEfficiency*     m_eff_func2;    //!< Efficiency functor 2
};

#endif /* GLATAEFF_HPP */
