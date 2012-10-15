/***************************************************************************
 *         GLATResponseTable.hpp  -  Fermi/LAT Response table class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @file GLATResponseTable.hpp
 * @brief Fermi/LAT response table class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATRESPONSETABLE_HPP
#define GLATRESPONSETABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GNodeArray.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GLATResponseTable
 *
 * @brief Interface for the Fermi LAT Response table class.
 *
 * A response table contains the binning information in energy and cos theta
 * for the Fermi LAT response function. From this binning information,
 * response values can be obtained using a bilinear interpolation in
 * log10 of energy and in cos theta.
 ***************************************************************************/
class GLATResponseTable : public GBase {

public:
    // Constructors and destructors
    GLATResponseTable(void);
    GLATResponseTable(const GLATResponseTable& table);
    virtual ~GLATResponseTable(void);

    // Operators
    GLATResponseTable& operator= (const GLATResponseTable & table);

    // Methods
    void                clear(void);
    GLATResponseTable*  clone(void) const;
    void                read(const GFitsTable* hdu);
    void                write(GFitsTable* hdu) const;
    int                 index(const int& ie, const int& ic) const;
    double              energy(const int& ie) const;
    void                set(const double& logE, const double& ctheta);
    double              interpolate(const double& logE, const double& ctheta, 
                                    const std::vector<double>& array);
    double              interpolate(const double& logE, const double& ctheta, 
                                    const std::vector<double>& array,
                                    const int& offset, const int& size);
    int                 size(void) const { return m_energy_num*m_ctheta_num; }
    int                 nenergies(void) const { return m_energy_num; }
    int                 ncostheta(void) const { return m_ctheta_num; }
    double              energy_lo(const int& inx) const;
    double              energy_hi(const int& inx) const;
    double              costheta_lo(const int& inx) const;
    double              costheta_hi(const int& inx) const;
    std::vector<int>    indices(void) const;
    std::vector<double> energies(void) const;
    std::vector<double> weights(void) const;
    std::string         print(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATResponseTable& table);
    void free_members(void);
    
    // Table nodes
    int                 m_energy_num;   //!< Number of energy bins in table
    int                 m_ctheta_num;   //!< Number of cos theta bins in table
    double*             m_energy_lo;    //!< Energy bins lower boundary (MeV)
    double*             m_energy_hi;    //!< Energy bins upper boundary (MeV)
    double*             m_ctheta_lo;    //!< cos(theta) bins lower boundary
    double*             m_ctheta_hi;    //!< cos(theta) bins upper boundary
    std::vector<double> m_energy;       //!< Energy nodes (MeV)
    GNodeArray          m_logE;         //!< Energy nodes (log10 mean energy)
    GNodeArray          m_ctheta;       //!< cos(theta) nodes

    // Bi-linear interpolation data
    double     m_last_energy;  //!< Last requested energy for interpolation
    double     m_last_ctheta;  //!< Last requested cos(theta) for interpolation
    int        m_inx1;         //!< Index 1
    int        m_inx2;         //!< Index 2
    int        m_inx3;         //!< Index 3
    int        m_inx4;         //!< Index 4
    double     m_wgt1;         //!< Weighting factor 1
    double     m_wgt2;         //!< Weighting factor 2
    double     m_wgt3;         //!< Weighting factor 3
    double     m_wgt4;         //!< Weighting factor 4
};

#endif /* GLATRESPONSETABLE_HPP */
