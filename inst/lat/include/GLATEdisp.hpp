/***************************************************************************
 *              GLATEdisp.hpp  -  Fermi LAT energy dispersion              *
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
 * @file GLATEdisp.hpp
 * @brief Fermi LAT energy dispersion class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GLATEDISP_HPP
#define GLATEDISP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <iostream>
#include "GBase.hpp"
#include "GLog.hpp"
#include "GLATPointing.hpp"
#include "GLATResponseTable.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
//#include "GSkyDir.hpp"
//#include "GEnergy.hpp"
//#include "GTime.hpp"


/***********************************************************************//**
 * @class GLATEdisp
 *
 * @brief Interface for the Fermi LAT energy dispersion.
 *
 * @todo Implement support for older response functions?
 ***************************************************************************/
class GLATEdisp : public GBase {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATEdisp& edisp);
    friend GLog&         operator<< (GLog& log, const GLATEdisp& edisp);

public:
    // Constructors and destructors
    GLATEdisp(void);
    GLATEdisp(const std::string& filename);
    GLATEdisp(const GLATEdisp& edisp);
    virtual ~GLATEdisp(void);

    // Operators
    GLATEdisp& operator= (const GLATEdisp& edisp);
    //double    operator() (const double& logE, const double& ctheta);
    //double    operator() (const GSkyDir& srcDir, const GEnergy& srcEng,
    //                      const GTime& srcTime, const GLATPointing& pnt);

    // Methods
    void         clear(void);
    GLATEdisp*   clone(void) const;
    void         load(const std::string& filename);
    void         save(const std::string& filename, bool clobber = false);
    void         read(const GFits& file);
    void         write(GFits& file) const;
    int          size(void) const { return nenergies()*ncostheta(); }
    int          nenergies(void) const { return m_edisp_bins.nenergies(); }
    int          ncostheta(void) const { return m_edisp_bins.ncostheta(); }
    //double       costhetamin(void) const { return m_min_ctheta; }
    //void         costhetamin(const double& ctheta);
    bool         hasphi(void) const { return false; }
    std::string  print(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATEdisp& edisp);
    void free_members(void);
    void read_edisp(const GFitsTable* hdu);
    void write_edisp(GFits& file) const;
    
    // Protected members
    GLATResponseTable   m_edisp_bins;   //!< Energy dispersion energy and cos theta binning
    std::vector<double> m_norm;         //!< Energy dispersion normalization
    std::vector<double> m_ls1;          //!< Energy dispersion ...
    std::vector<double> m_scale;        //!< Energy dispersion scaling parameters
};

#endif /* GLATEDISP_HPP */
