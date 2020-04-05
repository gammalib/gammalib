/***************************************************************************
 *            GSPIEventCube.hpp - INTEGRAL/SPI event cube class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIEventCube.hpp
 * @brief INTEGRAL/SPI event bin container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIEVENTCUBE_HPP
#define GSPIEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventCube.hpp"
#include "GSPIEventBin.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GEbounds;
class GGti;
class GFits;
class GSkyDir;
class GEnergy;
class GTime;
class GSPIInstDir;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSPIEventCube
 *
 * @brief INTEGRAL/SPI event bin container class
 *
 * This class is a container class for INTEGRAL/SPI event bins.
 ***************************************************************************/
class GSPIEventCube : public GEventCube {

public:
    // Constructors and destructors
    GSPIEventCube(void);
    explicit GSPIEventCube(const GFilename& filename);
    GSPIEventCube(const GSPIEventCube& cube);
    virtual ~GSPIEventCube(void);

    // Operators
    virtual GSPIEventCube&      operator=(const GSPIEventCube& cube);
    virtual GSPIEventBin*       operator[](const int& index);
    virtual const GSPIEventBin* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GSPIEventCube* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(const int& axis) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double ontime(void) const;
    double livetime(void) const;
    double model_counts(const int& index) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GSPIEventCube& cube);
    void         free_members(void);
    void         alloc_data(void);
    void         read_ebds(const GFitsTable* ebds);
    void         read_pnt(const GFitsTable* pnt, const GFitsTable* gti);
    void         read_gti(const GFitsTable* gti);
    void         read_dti(const GFitsTable* dti);
    void         read_dsp(const GFitsTable* dsp);
    void         read_models(const GFits& file);
    virtual void set_energies(void);
    virtual void set_times(void);
    void         init_bin(void);
    void         set_bin(const int& index);

    // Protected members
    GSPIEventBin             m_bin;        //!< Actual event bin
    int                      m_num_pt;     //!< Number of pointings
    int                      m_num_det;    //!< Number of detectors
    int                      m_num_ebin;   //!< Number of energy bins
    int                      m_num_sky;    //!< Number of sky models
    int                      m_num_bgm;    //!< Number of background models
    int                      m_gti_size;   //!< Size of GTI arrays
    int                      m_dsp_size;   //!< Size of DSP arrays
    int                      m_model_size; //!< Size of model arrays
    double*                  m_ontime;     //!< Ontime array
    double*                  m_livetime;   //!< Livetime array
    double*                  m_counts;     //!< Counts array
    double*                  m_stat_err;   //!< Statistical error array
    double*                  m_models;     //!< Models array
    double*                  m_size;       //!< Event bin size array
    GSPIInstDir*             m_dir;        //!< Event direction array
    GTime*                   m_time;       //!< Time array
    GEnergy*                 m_energy;     //!< Energy array
    GEnergy*                 m_ewidth;     //!< Energy bin width array
    std::vector<std::string> m_modnames;   //!< Model names
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSPIEventCube").
 ***************************************************************************/
inline
std::string GSPIEventCube::classname(void) const
{
    return ("GSPIEventCube");
}


/***********************************************************************//**
 * @brief Return event cube size
 *
 * @return Number of bins in event cube.
 *
 * Returns number of bins in event cube.
 ***************************************************************************/
inline
int GSPIEventCube::size(void) const
{
    return (m_dsp_size);
}


/***********************************************************************//**
 * @brief Return event cube dimension
 *
 * @return Number of dimensions in event cube.
 *
 * Returns the dimension of the event cube which is always 3.
 ***************************************************************************/
inline
int GSPIEventCube::dim(void) const
{
    return (3);
}


/***********************************************************************//**
 * @brief Set energies
 *
 * Sets energies.
 ***************************************************************************/
inline
void GSPIEventCube::set_energies(void)
{
    return;
}


/***********************************************************************//**
 * @brief Set times
 *
 * Sets times.
 ***************************************************************************/
inline
void GSPIEventCube::set_times(void)
{
    return;
}

#endif /* GSPIEVENTCUBE_HPP */
