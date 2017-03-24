/***************************************************************************
 *            GRmf.hpp - XSPEC Redistribution Matrix File class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2017 by Juergen Knoedlseder                         *
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
 * @file GRmf.hpp
 * @brief XSPEC Redistribution Matrix File class definition
 * @author Juergen Knoedlseder
 */

#ifndef GRMF_HPP
#define GRMF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFilename.hpp"
#include "GEbounds.hpp"
#include "GMatrixSparse.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;
class GFitsTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_rmf = "MATRIX";
}


/***********************************************************************//**
 * @class GRmf
 *
 * @brief Redistribution Matrix File class
 *
 * Matrix rows are true energy, columns are channels
 ***************************************************************************/
class GRmf : public GBase {

public:
    // Constructors and destructors
    GRmf(void);
    explicit GRmf(const GFilename& filename);
    GRmf(const GEbounds& etrue, const GEbounds& emeasured);
    GRmf(const GRmf& rmf);
    virtual ~GRmf(void);

    // Operators
    GRmf&         operator=(const GRmf& rmf);
    double&       operator()(const int& itrue, const int& imeasured);
    const double& operator()(const int& itrue, const int& imeasured) const;

    // Methods
    void                 clear(void);
    GRmf*                clone(void) const;
    std::string          classname(void) const;
    int                  size(void) const;
    int                  ntrue(void) const;
    int                  nmeasured(void) const;
    int                  itruemax(void) const;
    int                  imeasmax(void) const;
    double&              at(const int& itrue, const int& imeasured);
    const double&        at(const int& itrue, const int& imeasured) const;
    const GEbounds&      etrue(void) const;
    const GEbounds&      emeasured(void) const;
    GEbounds             etrue(const GEnergy& emeasured) const;
    GEbounds             emeasured(const GEnergy& etrue) const;
    const GMatrixSparse& matrix(void) const;
    void                 load(const GFilename& filename);
    void                 save(const GFilename&   filename,
                              const bool&        clobber = false,
                              const std::string& unit = "keV") const;
    void                 read(const GFitsTable& table);
    void                 write(GFits& fits,
                               const std::string& unit = "keV") const;
    const GFilename&     filename(void) const;
    std::string          print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GRmf& rmf);
    void   free_members(void);
    
    // Protected members
    mutable GFilename m_filename;      //!< Filename of origin
    GEbounds          m_ebds_true;     //!< True energy boundaries
    GEbounds          m_ebds_measured; //!< Measured energy boundaries
    GMatrixSparse     m_matrix;        //!< Sparse redistribution matrix
    int               m_imeasmax;      //!< Index of measured maximum
    int               m_itruemax;      //!< Index of true maximum
    
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GRmf").
 ***************************************************************************/
inline
std::string GRmf::classname(void) const
{
    return ("GRmf");
}


/***********************************************************************//**
 * @brief Return content of redistribution matrix bin
 *
 * @param[in] itrue True energy index [0,...,ntrue()-1].
 * @param[in] imeasured Measured energy index [0,...,nmeasured()-1].
 *
 * Returns reference to content of redistribution matrix bin bin with true
 * energy index @p itrue and measured energy index @p imeasured.
 ***************************************************************************/
inline
double& GRmf::operator()(const int& itrue, const int& imeasured)
{
    return (m_matrix(itrue, imeasured));
}


/***********************************************************************//**
 * @brief Return content of redistribution matrix bin (const version)
 *
 * @param[in] itrue True energy index [0,...,ntrue()-1].
 * @param[in] imeasured Measured energy index [0,...,nmeasured()-1].
 *
 * Returns reference to content of redistribution matrix bin bin with true
 * energy index @p itrue and measured energy index @p imeasured.
 ***************************************************************************/
inline
const double& GRmf::operator()(const int& itrue, const int& imeasured) const
{
    return (m_matrix(itrue, imeasured));
}


/***********************************************************************//**
 * @brief Return number of redistribution matrix bins
 *
 * @return Number of redistribution matrix bins.
 *
 * Returns the number of redistribution matrix bins.
 ***************************************************************************/
inline
int GRmf::size(void) const
{
    return (m_matrix.rows()*m_matrix.columns());
}


/***********************************************************************//**
 * @brief Return number of true energy bins in redistribution matrix
 *
 * @return Number true energy bins in redistribution matrix.
 *
 * Returns the number of true energy bins in redistribution matrix.
 ***************************************************************************/
inline
int GRmf::ntrue(void) const
{
    return (m_matrix.rows());
}


/***********************************************************************//**
 * @brief Return number of measured energy bins in redistribution matrix
 *
 * @return Number measured energy bins in redistribution matrix.
 *
 * Returns the number of measured energy bins in redistribution matrix.
 ***************************************************************************/
inline
int GRmf::nmeasured(void) const
{
    return (m_matrix.columns());
}




/***********************************************************************//**
 * @brief Return true energy index of maximum value of the redistribution matrix
 *
 * @return True energy index of maximum value of the redistribution matrix.
 *
 * Returns the true energy index of maximum value of the redistribution matrix.
 ***************************************************************************/
inline
int GRmf::itruemax(void) const
{
    return m_itruemax;
}


/***********************************************************************//**
 * @brief Return measured energy index of maximum value of the redistribution matrix
 *
 * @return Measured energy index of maximum value of the redistribution matrix.
 *
 * Returns the measured energy index of maximum value of the redistribution matrix.
 ***************************************************************************/
inline
int GRmf::imeasmax(void) const
{
    return m_imeasmax;
}


/***********************************************************************//**
 * @brief Return true energy boundaries
 *
 * @return True energy boundaries for redistribution matrix.
 *
 * Returns the true energy boundaries for redistribution matrix.
 ***************************************************************************/
inline
const GEbounds& GRmf::etrue(void) const
{
    return m_ebds_true;
}


/***********************************************************************//**
 * @brief Return measured energy boundaries
 *
 * @return Measured energy boundaries for redistribution matrix.
 *
 * Returns the measured energy boundaries for redistribution matrix.
 ***************************************************************************/
inline
const GEbounds& GRmf::emeasured(void) const
{
    return m_ebds_measured;
}


/***********************************************************************//**
 * @brief Return redistribution matrix
 *
 * @return Redistribution matrix.
 ***************************************************************************/
inline
const GMatrixSparse& GRmf::matrix(void) const
{
    return m_matrix;
}


/***********************************************************************//**
 * @brief Return file name
 *
 * @return File name from which the RMF information has been read or into
 *         which RMF information has been saved.
 *
 * Returns the file name from which the RMF information has been read or into
 * which RMF information has been saved. The returned string will be empty if
 * no load() or save() method has been called before.
 ***************************************************************************/
inline
const GFilename& GRmf::filename(void) const
{
    return (m_filename);
}

#endif /* GRMF_HPP */
