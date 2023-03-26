/***************************************************************************
 *             GCOMDris.hpp - COMPTEL Data Space container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2023 by Juergen Knodlseder                               *
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
 * @file GCOMDris.hpp
 * @brief COMPTEL Data Space container class definition
 * @author Juergen Knodlseder
 */

#ifndef GCOMDRIS_HPP
#define GCOMDRIS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GNdarray.hpp"
#include "GOptimizerFunction.hpp"
#include "GCOMDri.hpp"

/* __ Forward declarations _______________________________________________ */
class GOptimizer;
class GCOMObservation;
class GCOMEventList;
class GCOMSelection;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMDris
 *
 * @brief COMPTEL Data Space container class
 *
 * The COMPTEL Data Space container class holds instances of the COMPTEL
 * Data Space. It allows for an efficient computation of data and response
 * information for multiple energy bins.
 ***************************************************************************/
class GCOMDris : public GContainer {

public:
    // Constructors and destructors
    GCOMDris(void);
    GCOMDris(const GCOMDris& dris);
    virtual ~GCOMDris(void);

    // Operators
    GCOMDris&      operator=(const GCOMDris& dris);
    GCOMDri&       operator[](const int& index);
    const GCOMDri& operator[](const int& index) const;

    // Methods
    void           clear(void);
    GCOMDris*      clone(void) const;
    std::string    classname(void) const;
    GCOMDri&       at(const int& index);
    const GCOMDri& at(const int& index) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GCOMDri&       append(const GCOMDri& dri);
    GCOMDri&       insert(const int& index, const GCOMDri& dri);
    void           remove(const int& index);
    void           reserve(const int& num);
    void           extend(const GCOMDris& dris);
    void           compute_drws(const GCOMObservation& obs,
                                const GCOMSelection&   select = GCOMSelection(),
                                const double&          zetamin = 5.0,
                                const double&          timebin = 300.0,
                                const std::string&     method = "phibar");
    std::string    print(const GChatter& chatter = NORMAL) const;

    // Likelihood function
    class likelihood : public GOptimizerFunction {
    public:
        // Constructors and destructors
        likelihood(GCOMDris *dris, const int& ieng, const double& norm);

        // Implemented pure virtual base class methods
        virtual void           eval(const GOptimizerPars& pars);
        virtual double         value(void) const;
        virtual GVector*       gradient(void);
        virtual GMatrixSparse* curvature(void);

    protected:
        int           m_ieng;          //!< DRW energy bin
        double        m_norm;          //!< Normalisation value
        int           m_nsp;           //!< Number of superpackets
        int           m_nphibar;       //!< Number of phibar layers
        GNdarray      m_vetorate;      //!< Vetorate array multiplied by EHA cut correction
        GNdarray      m_activrate;     //!< Activation rate array multiplied by EHA cut correction
        GNdarray      m_diffrate;      //!< Vetorate - activation rate array
        GNdarray      m_vetorate_sum;  //!< Time integrated vetorate array
        GNdarray      m_activrate_sum; //!< Time integrated activation rate array
        GNdarray      m_diffrate_sum;  //!< Time integrated difference rate array
        double        m_value;         //!< Function value
        GVector       m_gradient;      //!< Gradient vector
        GMatrixSparse m_curvature;     //!< Curvature matrix
        GCOMDris*     m_this;          //!< Pointer to GCOMDris object
    };

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMDris& dris);
    void free_members(void);
    void compute_drws_energy(const GCOMObservation& obs,
                             const GCOMEventList*   events,
                             const GCOMSelection&   select,
                             const double&          zetamin,
                             const double&          timebin);
    void compute_drws_phibar(const GCOMObservation& obs,
                             const GCOMEventList*   events,
                             const GCOMSelection&   select,
                             const double&          zetamin,
                             const double&          timebin);
    void compute_drws_vetorate(const GCOMObservation& obs,
                               const GCOMEventList*   events,
                               const GCOMSelection&   select,
                               const double&          zetamin);
    void vetorate_setup(const GCOMObservation& obs,
                        const GCOMEventList*   events,
                        const GCOMSelection&   select,
                        const double&          zetamin);
    void vetorate_fit(void);
    void vetorate_update_activ(void);
    void vetorate_generate(const GCOMObservation& obs,
                           const GCOMSelection&   select,
                           const double&          zetamin);
    void vetorate_finish(const GCOMObservation& obs);
    void vetorate_save(const GFilename& filename) const;
    void vetorate_load(const GFilename& filename);

    // Protected data members
    std::vector<GCOMDri> m_dris; //!< Data space instances

    // Working arrays for vetorate computation
    GNdarray          m_wrk_counts;     //!< 3D event cube array
    GNdarray          m_wrk_ehacutcorr; //!< 2D geometry response array
    GNdarray          m_wrk_vetorate;   //!< 1D vetorates array
    GNdarray          m_wrk_activrate;  //!< 2D activation rate array
    GNdarray          m_wrk_rate;       //!< 2D rate array
    std::vector<bool> m_wrk_use_sp;     //!< 1D superpacket usage array
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMDris").
 ***************************************************************************/
inline
std::string GCOMDris::classname(void) const
{
    return ("GCOMDris");
}


/***********************************************************************//**
 * @brief Return reference to Data space instance
 *
 * @param[in] index Data space index [0,...,size()-1].
 *
 * Returns a reference to the Data space instance with the specified @p index.
 ***************************************************************************/
inline
GCOMDri& GCOMDris::operator[](const int& index)
{
    return (m_dris[index]);
}


/***********************************************************************//**
 * @brief Return reference to Data space instance (const version)
 *
 * @param[in] index Data space index [0,...,size()-1].
 *
 * Returns a reference to the Data space instance with the specified @p index.
 ***************************************************************************/
inline
const GCOMDri& GCOMDris::operator[](const int& index) const
{
    return (m_dris[index]);
}


/***********************************************************************//**
 * @brief Return number of Data space instances in container
 *
 * @return Number of Data space instances in container.
 *
 * Returns the number of Data space instances in the container.
 ***************************************************************************/
inline
int GCOMDris::size(void) const
{
    return (int)m_dris.size();
}


/***********************************************************************//**
 * @brief Signals if there are no Data space instances in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the Data space instances container does not contain any Data
 * space instances.
 ***************************************************************************/
inline
bool GCOMDris::is_empty(void) const
{
    return (m_dris.empty());
}


/***********************************************************************//**
 * @brief Reserves space for Data space instances in container
 *
 * @param[in] num Number of Data space instances.
 *
 * Reserves space for @p num Data space instances in the container.
 ***************************************************************************/
inline
void GCOMDris::reserve(const int& num)
{
    m_dris.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Return log-likelihood value of optimizer function
 *
 * @return Log-likelihood value of optimizer function.
 *
 * Returns the log-likelihood value of optimizer function.
 ***************************************************************************/
inline
double GCOMDris::likelihood::value(void) const
{
    return m_value;
}


/***********************************************************************//**
 * @brief Return pointer to gradient vector
 *
 * @return Pointer to gradient vector.
 *
 * Returns a pointer to the parameter gradient vector.
 ***************************************************************************/
inline
GVector* GCOMDris::likelihood::gradient(void)
{
    return &(m_gradient);
}


/***********************************************************************//**
 * @brief Return pointer to curvature matrix
 *
 * @return Pointer to curvature matrix.
 *
 * Returns a pointer to the parameter curvature matrix.
 ***************************************************************************/
inline
GMatrixSparse* GCOMDris::likelihood::curvature(void)
{
    return &(m_curvature);
}

#endif /* GCOMDRIS_HPP */
