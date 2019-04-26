/***************************************************************************
 *         GCTAModelSpatialLookup.hpp - Spatial lookup table model         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialLookup.hpp
 * @brief Spatial lookup table model interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELSPATIALLOOKUP_HPP
#define GCTAMODELSPATIALLOOKUP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GCTAModelSpatial.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GXmlElement;
class GCTAInstDir;
class GEnergy;
class GTime;
class GEbounds;
class GFitsTable;
class GFitsBinTable;
class GObservations;
class GCTAObservation;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_cta_spatial_lookup = "RADIAL BACKGROUND LOOKUP";
}


/***********************************************************************//**
 * @class GCTAModelSpatialLookup
 *
 * @brief Spatial lookup table model class
 ***************************************************************************/
class GCTAModelSpatialLookup  : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelSpatialLookup(void);
    explicit GCTAModelSpatialLookup(const GFilename& filename);
    explicit GCTAModelSpatialLookup(const GCTAResponseTable& table);
    explicit GCTAModelSpatialLookup(const GXmlElement& xml);
    GCTAModelSpatialLookup(const double&   maxrad,
                           const double&   radbin,
                           const GEbounds& ebds);
    GCTAModelSpatialLookup(const GCTAModelSpatialLookup& model);
    virtual ~GCTAModelSpatialLookup(void);

    // Operators
    virtual GCTAModelSpatialLookup& operator=(const GCTAModelSpatialLookup& model);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCTAModelSpatialLookup* clone(void) const;
    virtual std::string             classname(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const GCTAInstDir& dir,
                                         const GEnergy&     energy,
                                         const GTime&       time,
                                         const bool&        gradients = false) const;
    virtual double                  mc_max_value(const GCTAObservation& obs) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;
    virtual std::string             print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void                     fill(const GCTAObservation& obs);
    void                     fill(const GObservations& obs);
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFitsTable& table);
    void                     write(GFitsBinTable& table) const;
    void                     load(const GFilename& filename);
    void                     save(const GFilename& filename,
                                  const bool&      clobber = false) const;
    double                   norm(void) const;
    void                     norm(const double& norm);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelSpatialLookup& model);
    void free_members(void);
    void prepare_table(void);
    void normalise_table(void);
    int  table_index(const int& ienergy, const int& itheta) const;
    void fill_buffer(const GCTAObservation& obs, std::vector<double>& buffer);
    void set_from_buffer(const std::vector<double>& buffer);

    // Protected members
    GFilename         m_filename;    //!< Name of lookup table
    GCTAResponseTable m_lookup;      //!< Lookup table
    GModelPar         m_norm;        //!< Normalization factor
    int               m_inx_energy;  //!< Energy index
    int               m_inx_theta;   //!< Theta index
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAModelSpatialLookup").
 ***************************************************************************/
inline
std::string GCTAModelSpatialLookup::classname(void) const
{
    return ("GCTAModelSpatialLookup");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type "LookupTable".
 ***************************************************************************/
inline
std::string GCTAModelSpatialLookup::type(void) const
{
    return ("LookupTable");
}


/***********************************************************************//**
 * @brief Return maximum function value for Monte Carlo simulations
 *
 * @param[in] obs CTA Observation.
 * @return Maximum function value for Monte Carlo simulations.
 *
 * This method always returns 1.
 ***************************************************************************/
inline
double GCTAModelSpatialLookup::mc_max_value(const GCTAObservation& obs) const
{
    return 1.0;
}


/***********************************************************************//**
 * @brief Return lookup table
 *
 * @return Lookup table.
 *
 * Returns the lookup table.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTAModelSpatialLookup::table(void) const
{
    return m_lookup;
}


/***********************************************************************//**
 * @brief Get lookup table model normalisation
 *
 * @return Lookup table model normalisation.
 *
 * Returns the normalisation of the lookup table model.
 ***************************************************************************/
inline
double GCTAModelSpatialLookup::norm(void) const
{
    return (m_norm.value());
}


/***********************************************************************//**
 * @brief Set lookup table model normalisation
 *
 * @param[in] norm Lookup table model normalisation.
 *
 * Set the normalisation of the lookup table model.
 ***************************************************************************/
inline
void GCTAModelSpatialLookup::norm(const double& norm)
{
    m_norm.value(norm);
    return;
}

#endif /* GCTAMODELSPATIALLOOKUP_HPP */
