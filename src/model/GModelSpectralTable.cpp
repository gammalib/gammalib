/***************************************************************************
 *          GModelSpectralTable.cpp - Spectral table model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2020 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralTable.cpp
 * @brief Spectral table model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GRan.hpp"
#include "GEnergy.hpp"
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableStringCol.hpp"
#include "GFitsTableLongCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GModelSpectralTable.hpp"
#include "GModelSpectralTablePar.hpp"
#include "GModelSpectralTablePars.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralTable    g_spectral_table_seed;
const GModelSpectralRegistry g_spectral_table_registry(&g_spectral_table_seed);

/* __ Method name definitions ____________________________________________ */
#define G_CONST   "GModelSpectralTable(GEbounds&, GModelSpectralTablePars&, "\
                                                                 "GNdarray&)"
#define G_FLUX                "GModelSpectralTable::flux(GEnergy&, GEnergy&)"
#define G_EFLUX              "GModelSpectralTable::eflux(GEnergy&, GEnergy&)"
#define G_MC     "GModelSpectralTable::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                      "GModelSpectralTable::read(GXmlElement&)"
#define G_WRITE                    "GModelSpectralTable::write(GXmlElement&)"
#define G_LOAD                        "GModelSpectralTable::load(GFilename&)"
#define G_TABLE_PAR                    "GModelSpectralTable::table_par(int&)"
#define PAR_INDEX              "GModelSpectralTable::par_index(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_LOAD_SPEC                     //!< Debug load_spec() method
#define G_DEBUG_UPDATE                           //!< Debug update() method


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(void) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name of table model.
 * @param[in] norm Normalization factor.
 *
 * Constructs a spectral table model from a FITS file. See the load() method
 * for more information about the expected structure of the FITS file.
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(const GFilename& filename,
                                         const double&    norm) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Load table
    load(filename);

    // Set normalization
    m_norm.value(norm);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Table model constructor
 *
 * @param[in] ebounds Energy boundaries.
 * @param[in] pars Table model parameters.
 * @param[in] spectra Spectra.
 *
 * Constructs a spectral table model from energy boundaries, table model
 * parameters, and spectra.
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(const GEbounds&                ebounds,
                                         const GModelSpectralTablePars& pars,
                                         const GNdarray&                spectra) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Check consistency of spectra dimension
    if (spectra.dim() != pars.size()+1) {
        std::string msg = "Spectra array has "+gammalib::str(spectra.dim())+
                          " dimensions but an array with "+
                          gammalib::str(pars.size()+1)+" dimensions is "
                          "expected. Please specify a spectra array with the "
                          "correct dimension.";
        throw GException::invalid_argument(G_CONST, msg);
    }
    
    // Check consistency of parameter dimensions
    for (int i = 0; i < pars.size(); ++i) {
        int npars = pars[i]->size();
        int nspec = spectra.shape()[i];
        if (npars != nspec) {
            std::string msg = "Parameter \""+pars[i]->par().name()+"\" has "+
                              gammalib::str(npars)+" values but there are "+
                              gammalib::str(nspec)+" spectra for table axis "+
                              gammalib::str(i)+". Please specify a spectra "
                              "array with the correct number of spectra.";
            throw GException::invalid_argument(G_CONST, msg);
        }
    }
 
    // Check consistency of energy bins
    int nebins = ebounds.size();
    int nspec  = spectra.shape()[pars.size()];
    if (nebins != nspec) {
        std::string msg = "Spectra array has "+gammalib::str(nspec)+" energy "
                          "bins but there are "+gammalib::str(nebins)+
                          " energy boundaries. Please specify a spectra "
                          "array with the correct number of energy bins.";
        throw GException::invalid_argument(G_CONST, msg);
    }

    // Set members
    m_ebounds    = ebounds;
    m_table_pars = pars;
    m_spectra    = spectra;

    // Set energy nodes
    set_energy_nodes();

    // Set parameter pointers
    set_par_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a spectral table model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(const GXmlElement& xml) :
                     GModelSpectral()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Table model.
 ***************************************************************************/
GModelSpectralTable::GModelSpectralTable(const GModelSpectralTable& model) :
                     GModelSpectral(model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpectralTable::~GModelSpectralTable(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Table model.
 * @return Table model.
 ***************************************************************************/
GModelSpectralTable& GModelSpectralTable::operator=(const GModelSpectralTable& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpectral::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear table model
***************************************************************************/
void GModelSpectralTable::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpectral::free_members();

    // Initialise members
    this->GModelSpectral::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone table model
***************************************************************************/
GModelSpectralTable* GModelSpectralTable::clone(void) const
{
    // Clone table model
    return new GModelSpectralTable(*this);
}


/***********************************************************************//**
 * @brief Evaluate table model
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Model value (ph/cm2/s/MeV).
 *
 * Evaluates ...
 ***************************************************************************/
double GModelSpectralTable::eval(const GEnergy& srcEng,
                                 const GTime&   srcTime,
                                 const bool&    gradients) const
{
    // Evaluate function
    double func = 1.0; // TODO

    // Compute function value
    double value  = m_norm.value() * func;

    // Optionally compute gradients
    if (gradients) {

        // Compute partial derivatives of the parameter values
        double g_norm  = (m_norm.is_free())  ? m_norm.scale() * func : 0.0;

        // Set gradients
        m_norm.factor_gradient(g_norm);

    } // endif: gradient computation was requested

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GModelSpectralTable::eval";
        std::cout << "(srcEng=" << srcEng;
        std::cout << ", srcTime=" << srcTime << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", norm=" << norm();
        std::cout << ", func=" << func;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (units: ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Photon flux (ph/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 ***************************************************************************/
double GModelSpectralTable::flux(const GEnergy& emin,
                                 const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        // TODO

        // Multiply flux by normalisation factor
        flux *= m_norm.value();

    } // endif: integration range was valid

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (units: erg/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Energy flux (erg/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) E \, dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 ***************************************************************************/
double GModelSpectralTable::eflux(const GEnergy& emin,
                                  const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        // TODO

        // Multiply flux by normalisation factor
        eflux *= norm();

    } // endif: integration range was valid

    // Return flux
    return eflux;
}


/***********************************************************************//**
 * @brief Returns MC energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] time True photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Energy.
 *
 * @exception GException::invalid_argument
 *            Energy range is invalid (emin < emax required).
 *
 * Returns Monte Carlo energy by randomly drawing from a broken power law
 * defined by the file function.
 ***************************************************************************/
GEnergy GModelSpectralTable::mc(const GEnergy& emin,
                                const GEnergy& emax,
                                const GTime&   time,
                                GRan&          ran) const
{
    // Throw an exception if energy range is invalid
    if (emin >= emax) {
        std::string msg = "Minimum energy "+emin.print()+" is equal or "
                          "larger than maximum energy "+emax.print()+". "
                          "Please provide a minimum energy that is smaller "
                          "than the maximum energy.";
        throw GException::invalid_argument(G_MC, msg);
    }

    // Allocate energy
    GEnergy energy;

    // TODO

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing power law model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
 *
 * Reads the spectral information from an XML element. The format of the XML
 * elements is
 *
 *     <spectrum type="TableModel" file="..">
 *       <parameter name="Normalization" scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralTable::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Spectral function requires exactly 1 parameter.");
    }

    // Load table model from file (do this first since method calls clear())
    load(gammalib::xml_file_expand(xml, xml.attribute("file")));

    // Get parameter element
    const GXmlElement* par = xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization") {
        m_norm.read(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Require \"Normalization\" parameter.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "FileFunction"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 *
 *     <spectrum type="FileFunction" file="..">
 *       <parameter name="Normalization" scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 * 
 * Note that the function nodes will not be written since they will not be
 * altered by any method.
 ***************************************************************************/
void GModelSpectralTable::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \""+type()+"\".");
    }

    // If XML element has 0 nodes then append 1 parameter node
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Normalization\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Spectral function requires exactly 1 parameter.");
    }

    // Get parameter element
    GXmlElement* par = xml.element("parameter", 0);

    // Set parameter
    if (par->attribute("name") == "Normalization") {
        m_norm.write(*par);
    }
    else {
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Require \"Normalization\" parameter.");
    }

    // Set file attribute
    xml.attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load table from file
 *
 * @param[in] filename File name.
 *
 * Loads table model from FITS file. The format of the FITS file complies
 * with https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html
 ***************************************************************************/
void GModelSpectralTable::load(const GFilename& filename)
{
    // Clear table model
    clear();

    // Set filename
    m_filename = filename;

    // Continue only if filename is not empty
    if (!filename.is_empty()) {

        // Open FITS file
        GFits fits(filename);

        // Load data from extensions
        load_par(fits);
        load_eng(fits);
        load_spec(fits);

        // Close FITS file
        fits.close();

        // Set energy nodes
        set_energy_nodes();

        // Set parameter pointers
        set_par_pointers();

    } // endif: filename was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table into file
 *
 * @param[in] filename File name.
 * @param[in] clobber Overwrite existing file?
 *
 * Save the table model into a FITS file. The FITS file will contain three
 * binary table extensions:
 *
 *       * PARAMETERS - Table model parameters
 *       * ENERGIES - Table model energies
 *       * SPECTRA - Table model spectra
 *
 * The format of the FITS file complies with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html
 ***************************************************************************/
void GModelSpectralTable::save(const GFilename& filename,
                               const bool&      clobber) const
{
    // Create FITS file
    GFits fits;

    // Create binary tables
    GFitsBinTable par_table  = create_par_table();
    GFitsBinTable eng_table  = create_eng_table();
    GFitsBinTable spec_table = create_spec_table();

    // Append binary tables
    fits.append(par_table);
    fits.append(eng_table);
    fits.append(spec_table);

    // Set keywords in primary header
    GFitsHDU* primary = fits[0];
    primary->card("CONTENT",  "MODEL", "Spectrum file");
    primary->card("FILENAME", filename.url(), "FITS file name");
    primary->card("ORIGIN", PACKAGE_NAME, "Origin of FITS file");
    primary->card("MODLNAME", "model", "Model name");
    primary->card("MODLUNIT", "photons/cm^2/s/MeV", "Model units");
    primary->card("REDSHIFT", false, "If true then redshift will be included as a par");
    primary->card("ADDMODEL", true, "If true then this is an additive table model");
    primary->card("HDUCLASS", "OGIP", "Format conforms to OGIP standard");
    primary->card("HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC");
    primary->card("HDUVERS", "1.0.0", "Version of format");

    // Save to FITS file
    fits.saveto(filename, clobber);

    // Set filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return reference to table parameter
 *
 * @param[in] index Table parameter index.
 * @return Reference to table parameter.
 ***************************************************************************/
GModelSpectralTablePar& GModelSpectralTable::table_par(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_TABLE_PAR, index, size());
    }

    // Return reference
    return *(m_table_pars[index]);
}


/***********************************************************************//**
 * @brief Return const reference to table parameter
 *
 * @param[in] index Table parameter index.
 * @return Const reference to table parameter.
 ***************************************************************************/
const GModelSpectralTablePar& GModelSpectralTable::table_par(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_TABLE_PAR, index, size());
    }

    // Return reference
    return *(m_table_pars[index]);
}


/***********************************************************************//**
 * @brief Return reference to table parameter
 *
 * @param[in] name Table parameter name.
 * @return Reference to table parameter.
 ***************************************************************************/
GModelSpectralTablePar& GModelSpectralTable::table_par(const std::string& name)
{
    // Get index from name
    int index = par_index(name);

    // Return reference
    return *(m_table_pars[index]);
}


/***********************************************************************//**
 * @brief Return const reference to table parameter
 *
 * @param[in] name Table parameter name.
 * @return Const reference to table parameter.
 ***************************************************************************/
const GModelSpectralTablePar& GModelSpectralTable::table_par(const std::string& name) const
{
    // Get index from name
    int index = par_index(name);

    // Return reference
    return *(m_table_pars[index]);
}


/***********************************************************************//**
 * @brief Print table model information
 *
 * @param[in] chatter Chattiness.
 * @return String containing table model information.
 ***************************************************************************/
std::string GModelSpectralTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralTable ===");

        // Append information
        result.append("\n"+gammalib::parformat("Table file"));
        result.append(m_filename.url());
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(this->size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

        // Append parameter value intervals
        for (int i = 0; i < m_table_pars.size(); ++i) {
            std::string name = m_table_pars[i]->par().name();
            int         size = m_table_pars[i]->values().size();
            result.append("\n"+gammalib::parformat(name+" values"));
            result.append(gammalib::str(size));
            if (size > 0) {
                double min = m_table_pars[i]->values()[0];
                double max = m_table_pars[i]->values()[0];
                for (int k = 0; k < size; ++k) {
                    double value = m_table_pars[i]->values()[k];
                    if (value < min) {
                        min = value;
                    }
                    if (value > max) {
                        max = value;
                    }
                }
                result.append(" ["+gammalib::str(min)+", "+gammalib::str(max)+"]");
            }
            if (m_table_pars[i]->method() == 0) {
                result.append("  (linear)");
            }
            else if (m_table_pars[i]->method() == 1) {
                result.append("  (logarithmic)");
            }
        }

        // Append energy boundaries
        result.append("\n"+gammalib::parformat("Energies"));
        result.append(gammalib::str(m_ebounds.size()));
        result.append(" [");
        result.append(m_ebounds.emin().print());
        result.append(", ");
        result.append(m_ebounds.emax().print());
        result.append("]");

        // Append shape of spectra array
        int nspectra = 0;
        int nebins   = 0;
        if (m_spectra.dim() > 0) {
            nspectra = 1;
            for (int i = 0; i < m_spectra.dim()-1; ++i) {
                nspectra *= m_spectra.shape()[i];
            }
            nebins = m_spectra.shape()[m_spectra.dim()-1];
        }
        result.append("\n"+gammalib::parformat("Spectra array dimension"));
        result.append(gammalib::str(m_spectra.dim()));
        result.append("\n"+gammalib::parformat("Number of spectra"));
        result.append(gammalib::str(nspectra));
        result.append("\n"+gammalib::parformat("Number of spectral bins"));
        result.append(gammalib::str(nebins));

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralTable::init_members(void)
{
    // Initialise normalisation
    m_norm.clear();
    m_norm.name("Normalization");
    m_norm.scale(1.0);
    m_norm.value(1.0);
    m_norm.range(0.0,1000.0);
    m_norm.free();
    m_norm.gradient(0.0);
    m_norm.has_grad(true);

    // Initialize other members
    m_table_pars.clear();
    m_spectra.clear();
    m_ebounds.clear();
    m_filename.clear();

    // Initialize cache
    m_last_values.clear();
    m_lin_nodes.clear();
    m_log_nodes.clear();
    m_lin_values.clear();
    m_log_values.clear();

    // Set parameter pointer(s)
    set_par_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Table model.
 ***************************************************************************/
void GModelSpectralTable::copy_members(const GModelSpectralTable& model)
{
    // Copy members
    m_norm        = model.m_norm;
    m_table_pars  = model.m_table_pars;
    m_spectra     = model.m_spectra;
    m_ebounds     = model.m_ebounds;
    m_filename    = model.m_filename;

    // Copy cache
    m_last_values = model.m_last_values;
    m_lin_nodes   = model.m_lin_nodes;
    m_log_nodes   = model.m_log_nodes;
    m_lin_values  = model.m_lin_values;
    m_log_values  = model.m_log_values;

    // Set energy nodes
    set_energy_nodes();

    // Set parameter pointer(s)
    set_par_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralTable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set parameter pointers
 ***************************************************************************/
void GModelSpectralTable::set_par_pointers(void)
{
    // Clear pointers
    m_pars.clear();

    // Put normalisation parameter on stack
    m_pars.push_back(&m_norm);

    // Put table model parameters on stack
    for (int i = 0; i < m_table_pars.size(); ++i) {
        m_pars.push_back(&(m_table_pars[i]->par()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy nodes from energy boundaries
 ***************************************************************************/
void GModelSpectralTable::set_energy_nodes(void)
{
    // Determine number of energy bins
    int nebins = m_ebounds.size();

    // Continue only if there are energy bins
    if (nebins > 0) {

        // Initialise vectors for values
        m_lin_nodes  = GNodeArray();
        m_log_nodes  = GNodeArray();
        m_lin_values = std::vector<double>(nebins, 0.0);
        m_log_values = std::vector<double>(nebins, 0.0);

        // Compute node values
        for (int i = 0; i < nebins; ++i) {
            double energy_MeV       = m_ebounds.elogmean(i).MeV();
            double log10_energy_MeV = std::log10(energy_MeV);
            m_lin_nodes.append(energy_MeV);
            m_log_nodes.append(log10_energy_MeV);
        }

    } // endif: there were energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create PARAMETERS FITS table
 *
 * @return Binary FITS table containing table model parameters
 *
 * The method creates a binary FITS table that is compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node6.html
 ***************************************************************************/
GFitsBinTable GModelSpectralTable::create_par_table(void) const
{
    // Determine number of parameters
    int nrows = m_table_pars.size();

    // Determine maximum number of parameters values
    int ncols = 0;
    for (int i = 0; i < nrows; ++i) {
        int n = m_table_pars[i]->size();
        if (n > ncols) {
            ncols = n;
        }
    }

    // Create table columns
    GFitsTableStringCol col_name("NAME", nrows, 12);
    GFitsTableLongCol   col_method("METHOD", nrows);
    GFitsTableFloatCol  col_initial("INITIAL", nrows);
    GFitsTableFloatCol  col_delta("DELTA", nrows);
    GFitsTableFloatCol  col_minimum("MINIMUM", nrows);
    GFitsTableFloatCol  col_bottom("BOTTOM", nrows);
    GFitsTableFloatCol  col_top("TOP", nrows);
    GFitsTableFloatCol  col_maximum("MAXIMUM", nrows);
    GFitsTableLongCol   col_numbvals("NUMBVALS", nrows);
    GFitsTableFloatCol  col_value("VALUE", nrows, ncols);

    // Fill columns
    for (int i = 0; i < nrows; ++i) {
    
        // Get model parameter
        const GModelPar& par = m_table_pars[i]->par();

        // Set parameter name and initial value
        col_name(i)    = par.name();
        col_method(i)  = m_table_pars[i]->method();
        col_initial(i) = (float)par.value();

        // Handle free/fixed parameter attribute
        if (par.is_fixed()) {
            col_delta(i) = -1.0;
        }
        else {
            col_delta(i) =  1.0;  // Dummy value
        }

        // Set number of parameters
        col_numbvals(i) = m_table_pars[i]->size();

        // Set parameter values
        double min_value = 0.0;
        double max_value = 0.0;
        if (m_table_pars[i]->size() > 0) {
            min_value = m_table_pars[i]->values()[0];
            max_value = m_table_pars[i]->values()[0];
            for (int k = 0; k < m_table_pars[i]->size(); ++k) {
                double value = m_table_pars[i]->values()[k];
                if (value < min_value) {
                    min_value = value;
                }
                if (value > max_value) {
                    max_value = value;
                }
                col_value(i,k) = value;
            }
        }

        // Handle parameter minimum
        if (par.has_min()) {
            if (par.min() < min_value) {
                col_minimum(i) = (float)par.min();  // Hard minimum limit
                col_bottom(i)  = min_value;         // Soft minimum limit
            }
            else {
                col_minimum(i) = (float)par.min();  // Hard minimum limit
                col_bottom(i)  = (float)par.min();  // Soft minimum limit
            }
        }
        else {
            col_minimum(i) = min_value;             // Hard minimum limit
            col_bottom(i)  = min_value;             // Soft minimum limit
        }

        // Handle parameter maximum
        if (par.has_max()) {
            if (par.max() > max_value) {
                col_top(i)     = max_value;         // Soft maximum limit
                col_maximum(i) = (float)par.max();  // Hard maximum limit
            }
            else {
                col_top(i)     = (float)par.max();  // Soft maximum limit
                col_maximum(i) = (float)par.max();  // Hard maximum limit
            }
        }
        else {
            col_top(i)     = max_value;             // Soft maximum limit
            col_maximum(i) = max_value;             // Soft maximum limit
        }

    } // endfor: looped over parameters

    // Create binary table
    GFitsBinTable table;

    // Append columns to FITS table
    table.append(col_name);
    table.append(col_method);
    table.append(col_initial);
    table.append(col_delta);
    table.append(col_minimum);
    table.append(col_bottom);
    table.append(col_top);
    table.append(col_maximum);
    table.append(col_numbvals);
    table.append(col_value);

    // Set table extension name
    table.extname("PARAMETERS");

    // Set table keywords
    table.card("HDUCLASS", "OGIP", "Format conforms to OGIP standard");
    table.card("HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC");
    table.card("HDUCLAS2", "PARAMETERS", "Extension containing parameter info");
    table.card("HDUVERS", "1.0.0", "Version of format");
    table.card("NINTPARM", nrows, "Number of interpolation parameters");
    table.card("NADDPARM", 0, "Number of additional parameters");

    // Return table
    return table;
}


/***********************************************************************//**
 * @brief Create ENERGIES FITS table
 *
 * @return Binary FITS table containing table model energies
 *
 * The method creates a binary FITS table that is compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node7.html
 ***************************************************************************/
GFitsBinTable GModelSpectralTable::create_eng_table(void) const
{
    // Create table columns
    GFitsTableFloatCol col_lo("ENERG_LO", m_ebounds.size());
    GFitsTableFloatCol col_hi("ENERG_HI", m_ebounds.size());

    // Fill energy boundary columns
    for (int i = 0; i < m_ebounds.size(); ++i) {
        col_lo(i) = m_ebounds.emin(i).keV();
        col_hi(i) = m_ebounds.emax(i).keV();
    }

    // Set energy units
    col_lo.unit("keV");
    col_hi.unit("keV");

    // Create binary table
    GFitsBinTable table;

    // Append columns to FITS table
    table.append(col_lo);
    table.append(col_hi);

    // Set table extension name
    table.extname("ENERGIES");

    // Set table keywords
    table.card("HDUCLASS", "OGIP", "Format conforms to OGIP standard");
    table.card("HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC");
    table.card("HDUCLAS2", "ENERGIES", "Extension containing energy bin info");
    table.card("HDUVERS", "1.0.0", "Version of format");

    // Return table
    return table;
}


/***********************************************************************//**
 * @brief Create SPECTRA FITS table
 *
 * @return Binary FITS table containing table model spectra
 *
 * The method creates a binary FITS table that is compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node8.html
 ***************************************************************************/
GFitsBinTable GModelSpectralTable::create_spec_table(void) const
{
    // Compute number of rows
    int nrows = 1;
    for (int i = 0; i < m_spectra.dim()-1; ++i) {
        nrows *= m_spectra.shape()[i];
    }

    // Compute number of parameters
    int npars = m_table_pars.size();

    // Compute number of energy bins
    int nebins = m_ebounds.size();

    // Create columns
    GFitsTableFloatCol col_pars("PARAMVAL", nrows, npars);
    GFitsTableFloatCol col_spec("INTPSPEC", nrows, nebins);

    // Fill columns
    std::vector<int> inx(npars+1,0);
    for (int i = 0; i < nrows; ++i) {

        // Set parameter values
        for (int k = 0; k < npars; ++k) {
            col_pars(i,k) = m_table_pars[k]->values()[inx[k]];
        }

        // Set current index
        std::vector<int> index = inx;

        // Loop over energy bins
        for (int k = 0; k < nebins; ++k, ++index[npars]) {
            col_spec(i,k) = m_spectra(index);
        }

        // Increment parameter index. Last parameter index is changing fastest
        int ipar = npars-1;
        do {
            inx[ipar] += 1;
            if (inx[ipar] < m_spectra.shape()[ipar]) {
                break;
            }
            else {
                inx[ipar] = 0;
                ipar--;
            }
        } while (ipar >= 0);

    } // endfor: looped over rows

    // Create binary table
    GFitsBinTable table;

    // Append columns to FITS table
    table.append(col_pars);
    table.append(col_spec);

    // Set table extension name
    table.extname("SPECTRA");

    // Set table keywords
    table.card("HDUCLASS", "OGIP",              "Format conforms to OGIP standard");
    table.card("HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC");
    table.card("HDUCLAS2", "MODEL SPECTRA",     "Extension containing model spectra");
    table.card("HDUVERS",  "1.0.0",             "Version of format");

    // Return table
    return table;
}


/***********************************************************************//**
 * @brief Load data from PARAMETERS extension
 *
 * @param[in] fits FITS file.
 *
 * The method loads data from the PARAMETERS binary table. The format of the
 * table needs to be compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node6.html
 ***************************************************************************/
void GModelSpectralTable::load_par(const GFits& fits)
{
    // Get PARAMETERS extension
    const GFitsTable& table = *fits.table("PARAMETERS");

    // Get number of parameters
    int npars = table.integer("NAXIS2");

    // Loop over all parameters
    for (int i = 0; i < npars; ++i) {

        // Set model parameter
        GModelPar par(table["NAME"]->string(i), table["INITIAL"]->real(i));

        // Apply hard minimum and maximum as parameter range. Note that the
        // minimum and maximum apply to the value factor, hence in case of
        // a negative scale factor the minimum becomes the maximum and vice
        // versa. This is actually a bug in GammaLib, see
        // https://cta-redmine.irap.omp.eu/issues/3072
        double min = table["MINIMUM"]->real(i) / par.scale();
        double max = table["MAXIMUM"]->real(i) / par.scale();
        if (min > max) {
            par.factor_range(max, min);
        }
        else {
            par.factor_range(min, max);
        }

        // Fix or free parameter according to DELTA value
        if (table["DELTA"]->real(i) < 0.0) {
            par.fix();
        }
        else {
            par.free();
        }

        // Set vector of parameter values
        std::vector<double> values;
        for (int k = 0; k < table["NUMBVALS"]->integer(i); ++k) {
            values.push_back(table["VALUE"]->real(i,k));
        }

        // Set table model parameter
        GModelSpectralTablePar table_model_par(par, values);

        // Set interpolation method
        table_model_par.method(table["METHOD"]->integer(i));

        // Append table model parameter
        m_table_pars.append(table_model_par);

    } // endfor: looped over all parameters

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data from ENERGIES extension
 *
 * @param[in] fits FITS file.
 *
 * The method loads data from the ENERGIES binary table. The format of the
 * table needs to be compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node7.html
 *
 * If energy units are provided in the ENERGIES extension the method decodes
 * the units and interprets the energy values correspondingly.
 ***************************************************************************/
void GModelSpectralTable::load_eng(const GFits& fits)
{
    // Get ENERGIES extension
    const GFitsTable& table = *fits.table("ENERGIES");

    // Get number of energy bins
    int nebins = table.integer("NAXIS2");

    // Extract energy boundary information from FITS table
    if (nebins > 0) {

        // Get units
        std::string emin_unit = table["ENERG_LO"]->unit();
        std::string emax_unit = table["ENERG_HI"]->unit();
        if (emin_unit.empty()) {
            emin_unit = "keV";
        }
        if (emax_unit.empty()) {
            emax_unit = "keV";
        }

        // Read energy boundaries and append bins
        for (int i = 0; i < nebins; ++i) {
            GEnergy emin(table["ENERG_LO"]->real(i), emin_unit);
            GEnergy emax(table["ENERG_HI"]->real(i), emax_unit);
            m_ebounds.append(emin, emax);
        }

    } // endif: there were energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data from SPECTRA extension
 *
 * @param[in] fits FITS file.
 *
 * The method loads data from the SPECTRA binary table. The format of the
 * table needs to be compliant with
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/node8.html
 ***************************************************************************/
void GModelSpectralTable::load_spec(const GFits& fits)
{
    // Get number of parameters and energy bins
    int npars  = fits.table("PARAMETERS")->integer("NAXIS2");
    int nebins = fits.table("ENERGIES")->integer("NAXIS2");

    // Setup dimension of spectra array
    std::vector<int> naxis(npars+1, 0);
    for (int i = 0; i < npars; ++i) {
        naxis[i] = (*fits.table("PARAMETERS"))["NUMBVALS"]->integer(i);
    }
    naxis[npars] = nebins;

    // Setup spectra array
    m_spectra = GNdarray(naxis);

    // Get SPECTRA extension
    const GFitsTable& table = *fits.table("SPECTRA");

    // Get number of rows
    int nrows = table.integer("NAXIS2");

    // Loop over rows
    for (int i = 0; i < nrows; ++i) {

        // Initialise spectra array index of first energy bin
        std::vector<int> index(npars+1, 0);

        // Setup spectra array index
        std::vector<double> parval(npars, 0.0);
        for (int k = 0; k < npars; ++k) {

            // Get reference to node array
            const GNodeArray& nodes = m_table_pars[k]->values();

            // Set interpolation value
            nodes.set_value(table["PARAMVAL"]->real(i,k));

            // Get best index
            int inx = (nodes.wgt_left() > nodes.wgt_right()) ? nodes.inx_left()
                                                             : nodes.inx_right();

            // Store index
            index[k] = inx;

        } // endfor: looped over parameters

        // Debug: dump vector array
        #if defined(G_DEBUG_LOAD_SPEC)
        std::cout << i << ": (";
        for (int k = 0; k < index.size(); ++k) {
            if (k > 0) {
                std::cout << ", ";
            }
            std::cout << index[k];
        }
        std::cout << ")" << std::endl;
        #endif

        // Loop over energy bins and store spectrum
        for (int k = 0; k < nebins; ++k, ++index[npars]) {
            m_spectra(index) = table["INTPSPEC"]->real(i,k);
        }

    } // endfor: looped over rows

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return index for parameter name
 *
 * @param[in] name Parameter name.
 * @return Parameter index.
 *
 * @exception GException::invalid_argument
 *            Parameter name not found in spectral table.
 ***************************************************************************/
int GModelSpectralTable::par_index(const std::string& name) const
{
    // Get parameter index
    int index = 0;
    for (; index < size(); ++index) {
        if (m_table_pars[index]->par().name() == name) {
            break;
        }
    }

    // Throw exception if parameter name was not found
    if (index >= size()) {
        throw GException::par_not_found(PAR_INDEX, name);
        std::string msg = "Parameter name \""+name+"\" not found in spectral "
                          "table. Please specify one of the following parameter "
                          "names:";
        for (int i = 0; i < size(); ++i) {
            if (i > 0) {
                msg += ",";
            }
            msg += " \""+m_table_pars[i]->par().name()+"\"";
        }
        throw GException::invalid_argument(PAR_INDEX, msg);
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Update
 ***************************************************************************/
void GModelSpectralTable::update(void)
{
    // Debug option: write header
    #if defined(G_DEBUG_UPDATE)
    std::cout << "GModelSpectralTable::update() entered" << std::endl;
    #endif

    // Get dimension of spectral table
    int dim = m_table_pars.size();

    // Initialise update flag
    bool need_update = false;

    // If dimension of last cached parameter values differ from dimension
    // of spectral table then reallocate cache and request update
    if (m_last_values.size() != dim) {
        m_last_values = std::vector<double>(dim, 0.0);
        need_update   = true;
    }

    // ... otherwise check if some parameter values have changed
    else {
        for (int i = 0; i < dim; ++i) {
            if (m_table_pars[i]->par().value() != m_last_values[i]) {
                need_update = true;
                break;
            }
        }
    }

    // Continue only if update is required
    if (need_update) {

        // Initialise vectors for weights and indices
        std::vector<double> m_weights(2*dim, 0.0);
        std::vector<int>    m_indices(2*dim, 0);

        // Loop over all parameters
        for (int i = 0; i < dim; ++i) {

            // Get pointers to node array and parameter
            const GNodeArray* nodes = &(m_table_pars[i]->values());
            const GModelPar*  par   = &(m_table_pars[i]->par());

            // Get parameter value
            double value = par->value();

            // Set values for node array
            nodes->set_value(value);

            // Cache parameter value
            m_last_values[i] = value;

            // Compute left and right indices
            int il = 2*i;
            int ir = il + 1;

            // Push back weigths and indices
            m_weights[il] = nodes->wgt_left();
            m_weights[ir] = nodes->wgt_right();
            m_indices[il] = nodes->inx_left();
            m_indices[ir] = nodes->inx_right();

            // Debug option: print weights and indices
            #if defined(G_DEBUG_UPDATE)
            std::cout << " wgt_l=" << m_weights[il];
            std::cout << " wgt_r=" << m_weights[ir];
            std::cout << " inx_l=" << m_indices[il];
            std::cout << " inx_r=" << m_indices[ir] << std::endl;
            #endif

        } // endfor: looped over all parameters

        // Compute number of combinations
        int combinations = 1 << dim;

        // Initialise vectors for values
        m_lin_values = std::vector<double>(ebounds().size(), 0.0);
        m_log_values = std::vector<double>(ebounds().size(), 0.0);

        // Debug option: initial sum of weights
        #if defined(G_DEBUG_UPDATE)
        double weight_sum = 0.0;
        #endif

        // Loop over combinations
        for (int i = 0; i < combinations; ++i) {

            // Debug option: start printing combination
            #if defined(G_DEBUG_UPDATE)
            std::cout << " " << i << ": ";
            #endif

            // Initialise weight
            double weight = 1.0;

            // Initialise index vector (including the energy dimension)
            std::vector<int> index_shape(dim+1,0);

            // Loop over dimensions
            for (int k = 0, div = 1; k < dim; ++k, div *= 2) {

                // Compute index for each dimension
                int index = i/div % 2 + k * 2;

                // Update weight
                weight *= m_weights[index];

                // Add index
                index_shape[k] = m_indices[index];

                // Debug option: print information for dimension
                #if defined(G_DEBUG_UPDATE)
                std::cout << index;
                std::cout << " (" << m_weights[index];
                std::cout << " @ " << m_indices[index] << ")";
                #endif

            } // endfor: looped over dimensions

            // Get index of spectrum
            int index_spectra = m_spectra.index(index_shape);

            // Add


            // Debug option: print total weight and index for combination
            #if defined(G_DEBUG_UPDATE)
            std::cout << ": wgt=" << weight;
            std::cout << " (" << index_spectra << ")" << std::endl;
            weight_sum += weight;
            #endif
    
        } // endfor: looped over combinations

        // Debug option: print sum of weights
        #if defined(G_DEBUG_UPDATE)
        std::cout << " sum(wgt)=" << weight_sum << std::endl;
        #endif

    } // endif: updated requested

    // Return
    return;
}
