/***************************************************************************
 *            GModelData.cpp - Abstract virtual data model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2020 by Juergen Knoedlseder                         *
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
 * @file GModelData.cpp
 * @brief GModelData class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelData.hpp"
#include "GVector.hpp"
#include "GMatrixSparse.hpp"
#include "GObservation.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_EVAL                  "GModel::eval(GObservation&, GMatrixSparse*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelData::GModelData(void) : GModel()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 ***************************************************************************/
GModelData::GModelData(const GXmlElement& xml) : GModel(xml)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Data model.
 ***************************************************************************/
GModelData::GModelData(const GModelData& model) : GModel(model)
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
GModelData::~GModelData(void)
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
 * @param[in] model Data model.
 * @return Data model.
 ***************************************************************************/
GModelData& GModelData::operator=(const GModelData& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModel::operator=(model);

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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return model values and gradients
 *
 * @param[in] obs Observation.
 * @param[out] gradients Pointer to matrix of gradients.
 * @return Model values.
 *
 * Evaluates the model values and parameter gradients for all events in an
 * observation. Gradients are only returned if the @p gradients pointer is
 * not NULL.
 *
 * The matrix of gradients is a sparse matrix where the number of rows
 * corresponds to the number of model parameters and the number of columns
 * corresponds to the number of events. An exception is thrown if the
 * dimension of the @p gradients matrix is not compatible with the model
 * and the observations.
 ***************************************************************************/
GVector GModelData::eval(const GObservation& obs,
                         GMatrixSparse*      gradients) const
{
    // Get number of model parameters and number of events
    int npars   = size();
    int nevents = obs.events()->size();

    // Initialise gradients flag
    bool grad = ((gradients != NULL) && (npars > 0));

    // Check matrix consistency
    if (grad) {
        if (gradients->columns() != npars) {
            std::string msg = "Number of "+gammalib::str(gradients->columns())+
                              " columns in gradient matrix differs from number "
                              "of "+gammalib::str(npars)+" parameters "
                              "in model. Please provide a compatible gradient "
                              "matrix.";
            throw GException::invalid_argument(G_EVAL, msg);
        }
        if (gradients->rows() != nevents) {
            std::string msg = "Number of "+gammalib::str(gradients->rows())+
                              " rows in gradient matrix differs from number "
                              "of "+gammalib::str(nevents)+" events in "
                              "observation. Please provide a compatible "
                              "gradient matrix.";
            throw GException::invalid_argument(G_EVAL, msg);
        }
    }

    // Initialise values
    GVector values(nevents);

    // Initialise temporary vectors to hold gradients
    GVector* tmp_gradients = NULL;
    if (grad) {

        // Allocate temporary vectors to hold the gradients
        tmp_gradients = new GVector[npars];
        for (int i = 0; i < npars; ++i) {
            tmp_gradients[i] = GVector(nevents);
        }

        // Signal all parameters that will have analytical gradients. These
        // are all parameters that are free and for which the model provides
        // analytical gradients.
        for (int i = 0; i < npars; ++i) {
            const GModelPar& par = (*this)[i];
            if (par.is_free() && par.has_grad()) {
                obs.computed_gradient(par);
            }
        }

    } // endif: gradients were requested

    // Loop over events
    for (int k = 0; k < nevents; ++k) {

        // Get reference to event
        const GEvent& event = *((*obs.events())[k]);

        // Evaluate probability
        values[k] = eval(event, obs, grad);

        // Set gradients
        if (grad) {

            // Extract relevant parameter gradients
            for (int i = 0; i < npars; ++i) {
                const GModelPar& par = (*this)[i];
                if (par.is_free() && par.has_grad()) {
                    tmp_gradients[i][k] = par.factor_gradient();
                }
            }

        } // endif: looped over gradients

    } // endfor: looped over events

    // Post-process gradients
    if (grad) {

        // Fill gradients into matrix
        for (int i = 0; i < npars; ++i) {
            gradients->column(i, tmp_gradients[i]);
        }

        // Delete temporal gradients
        delete [] tmp_gradients;

    } // endif: post-processed gradients

    // Return values
    return values;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelData::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Data model.
 ***************************************************************************/
void GModelData::copy_members(const GModelData& model)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelData::free_members(void)
{
    // Return
    return;
}
