/***************************************************************************
 *            GCTAEdisp2D.hpp - CTA 2D energy dispersion class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2018 by Florent Forest                              *
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
 * @file GCTAEdisp2D.hpp
 * @brief CTA 2D energy dispersion class definition
 * @author Florent Forest
 */

#ifndef GCTAEDISP2D_HPP
#define GCTAEDISP2D_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFilename.hpp"
#include "GFunction.hpp"
#include "GEnergy.hpp"
#include "GCTAEdisp.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFft;
class GRan;
class GFits;
class GFitsBinTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_cta_edisp2d = "ENERGY DISPERSION";
}


/***********************************************************************//**
 * @class GCTAEdisp2D
 *
 * @brief CTA 2D energy dispersion class
 *
 * This class implements the energy dispersion for the CTA 2D response. The
 * CTA 2D energy dispersion is in fact a 3-dimensional function, where the
 * energy dispersion
 *
 * \f[
 *    E_{\rm disp}(E_{\rm reco} | E_{\rm true}, \theta)
 * \f]
 *
 * is given as function of true energy \f$E_{\rm true}\f$, reconstructed
 * energy \f$E_{\rm reco}\f$ and offset angle \f$\theta\f$. The
 * GCTAEdisp2D::operator() returns the value of the function
 * \f$E_{\rm disp}(E_{\rm reco} | E_{\rm true}, \theta)\f$ in units of
 * MeV\f$^{-1}\f$ with
 *
 * \f[
 *    \int_{E_{\rm reco}^{\rm min}}^{E_{\rm reco}^{\rm max}}
 *      E_{\rm disp}(E_{\rm reco} | E_{\rm true}, \theta) \, dE_{\rm reco} = 1
 * \f]
 *
 * The energy dispersion is stored internally as a 3-dimensional response
 * table, spanned by true energy \f$E_{\rm true}\f$, migration
 * \f$m=E_{\rm reco}/E_{\rm true}\f$ and offset angle \f$\theta\f$, with
 *
 * \f[
 *    E_{\rm disp}(E_{\rm reco} | E_{\rm true}, \theta) =
 *    \frac{E_{\rm disp}(m | E_{\rm true}, \theta)}{E_{\rm true}}
 * \f]
 *
 * and
 *
 * \f[
 *    \int_{m_{\rm min}}^{m_{\rm max}}
 *      E_{\rm disp}(m | E_{\rm true}, \theta) \, dm = 1
 * \f]
 ***************************************************************************/
class GCTAEdisp2D : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdisp2D(void);
    explicit GCTAEdisp2D(const GFilename& filename);
    GCTAEdisp2D(const GCTAEdisp2D& edisp);
    virtual ~GCTAEdisp2D(void);

    // Operators
    GCTAEdisp2D& operator=(const GCTAEdisp2D& edisp);
    double operator()(const GEnergy& ereco,
                      const GEnergy& etrue,
                      const double&  theta = 0.0,
                      const double&  phi = 0.0,
                      const double&  zenith = 0.0,
                      const double&  azimuth = 0.0) const;

    // Implemented methods
    void         clear(void);
    GCTAEdisp2D* clone(void) const;
    std::string  classname(void) const;
    void         load(const GFilename& filename);
    GFilename    filename(void) const;
    GEnergy      mc(GRan&          ran,
                    const GEnergy& etrue,
                    const double&  theta = 0.0,
                    const double&  phi = 0.0,
                    const double&  zenith = 0.0,
                    const double&  azimuth = 0.0) const;
    GEbounds     ereco_bounds(const GEnergy& etrue,
                              const double&  theta = 0.0,
                              const double&  phi = 0.0,
                              const double&  zenith = 0.0,
                              const double&  azimuth = 0.0) const;
    GEbounds     etrue_bounds(const GEnergy& ereco,
                              const double&  theta = 0.0,
                              const double&  phi = 0.0,
                              const double&  zenith = 0.0,
                              const double&  azimuth = 0.0) const;
    double       prob_erecobin(const GEnergy& ereco_min,
                               const GEnergy& ereco_max,
                               const GEnergy& etrue,
                               const double&  theta) const;
    std::string  print(const GChatter& chatter = NORMAL) const;

    // Methods
    void                     fetch(void) const;
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFitsTable& table);
    void                     write(GFitsBinTable& table) const;
    void                     save(const GFilename& filename,
                                  const bool&      clobber = false) const;

private:
    // Methods
    void    init_members(void);
    void    copy_members(const GCTAEdisp2D& edisp);
    void    free_members(void);
    GEnergy etrue(const int& ietrue) const;
    double  migra(const int& imigra) const;
    double  theta(const int& itheta) const;
    void    compute_ereco_bounds(const double& theta = 0.0,
                                 const double& phi = 0.0,
                                 const double& zenith = 0.0,
                                 const double& azimuth = 0.0) const;
    void    compute_etrue_bounds(const double& theta = 0.0,
                                 const double& phi = 0.0,
                                 const double& zenith = 0.0,
                                 const double& azimuth = 0.0) const;
    void    set_table(void);
    void    set_boundaries(void);
    void    set_max_edisp(void);
    double  get_max_edisp(const GEnergy& etrue, const double&  theta) const;
    void    normalize_table(void);
    int     table_index(const int& ietrue,
                        const int& imigra,
                        const int& itheta) const;
    int     table_stride(const int& axis) const;
    double  table_value(const int& base_ll,
                        const int& base_lr,
                        const int& base_rl,
                        const int& base_rr,
                        const double& wgt_el,
                        const double& wgt_er,
                        const double& wgt_tl,
                        const double& wgt_tr,
                        const int& offset) const;
    double  table_value(const int& base_ll,
                        const int& base_lr,
                        const int& base_rl,
                        const int& base_rr,
                        const double& wgt_el,
                        const double& wgt_er,
                        const double& wgt_tl,
                        const double& wgt_tr,
                        const double& migra) const;

    // Kludge
    void     smooth_table(void);
    GNdarray smooth_array(const GNdarray& array,
                          const int&      nstart,
                          const int&      nstop,
                          const double&   limit) const;
    double   smoothed_array_value(const int&      inx,
                                  const GNdarray& array,
                                  const int&      nodes,
                                  const double&   limit) const;
    void     get_moments(const int& itheta,
                         GNdarray*  mean,
                         GNdarray*  rms,
                         GNdarray*  total) const;
    void     get_mean_rms(const GNdarray& array, double* mean, double* rms) const;
    GNdarray gaussian_array(const double& mean,
                            const double& rms,
                            const double& total) const;

    // Protected classes
    class edisp_ereco_kern : public GFunction {
    public:
        edisp_ereco_kern(const GCTAEdisp2D* parent,
                         const GEnergy&     etrue,
                         const double&      theta) :
                         m_parent(parent),
                         m_etrue(etrue),
                         m_theta(theta) { }
        double eval(const double& log10Ereco);
    protected:
        const GCTAEdisp2D* m_parent;  //!< Pointer to parent class
        GEnergy            m_etrue;   //!< True photon energy
        double             m_theta;   //!< Offset angle
    };

    // Members
    mutable GFilename   m_filename;    //!< Name of Edisp response file
    GCTAResponseTable   m_edisp;       //!< Edisp response table
    mutable bool        m_fetched;     //!< Signals that Edisp has been fetched
    int                 m_inx_etrue;   //!< True energy index
    int                 m_inx_migra;   //!< Migration index
    int                 m_inx_theta;   //!< Theta index
    int                 m_inx_matrix;  //!< Matrix
    double              m_logEsrc_min; //!< Minimum logE (log10(E/TeV))
    double              m_logEsrc_max; //!< Maximum logE (log10(E/TeV))
    double              m_migra_min;   //!< Minimum migration
    double              m_migra_max;   //!< Maximum migration
    double              m_theta_min;   //!< Minimum theta (radians)
    double              m_theta_max;   //!< Maximum theta (radians)
    std::vector<double> m_max_edisp;   //!< Maximum Edisp

    // Computation cache
    mutable bool                  m_ereco_bounds_computed;
    mutable bool                  m_etrue_bounds_computed;
    mutable double                m_last_theta_ereco;
    mutable double                m_last_theta_etrue;
    mutable GEnergy               m_last_etrue;
    mutable GEnergy               m_last_ereco;
    mutable int                   m_index_ereco;
    mutable int                   m_index_etrue;
    mutable std::vector<GEbounds> m_ereco_bounds;
    mutable std::vector<GEbounds> m_etrue_bounds;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAEdisp2D").
 ***************************************************************************/
inline
std::string GCTAEdisp2D::classname(void) const
{
    return ("GCTAEdisp2D");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Name of FITS file from which energy dispersion was loaded.
 ***************************************************************************/
inline
GFilename GCTAEdisp2D::filename(void) const
{
    // Return filename
    return (m_filename);
}


/***********************************************************************//**
 * @brief Return response table
 *
 * @return Reference to response table.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTAEdisp2D::table(void) const
{
    fetch();
    return (m_edisp);
}

#endif /* GCTAEDISP2D_HPP */
