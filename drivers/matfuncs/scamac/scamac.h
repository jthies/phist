#ifndef SCAMAC_H
#define SCAMAC_H

/** \mainpage ScaMaC: A Scalable Matrix Collection
 *
 * \section overview_sec Overview
 * The scalable matrix collection ScaMaC provides routines for the generation of matrices for test and benchmark applications.
 * \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 * \date   October 2017 (of last major changes)
 */

/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ScaMaC top-level header file
 */

/** \defgroup library ScaMaC library
 *  \brief    High-level functions to access the ScaMaC
 *  \details  Only these functions should be called by the user.
 */

/** \defgroup toolbox ScaMaC toolbox
 *  \brief    Toolbox functions of the ScaMaC
 *  \details  These functions provide additional functionality for the investigation of the ScaMaC matrices. They are used in the "scamact" ScaMaC Toolbox program.
 */

/** \defgroup matrix ScaMaC matrices
 *  \brief    Internal matrix generators of the ScaMaC
 *  \details  These generators should not be called directly by the user.
 */

/** \defgroup internal ScaMaC internal
 *  \brief    Internal functions of the ScaMaC
 *  \details  These functions should not be called directly by the user.
 */

/** \defgroup mwe ScaMaC min. work. ex.
 *  \brief    Minimal working examples for ScaMaC functions
 *  \details  These examples demonstrate the basic usage of the ScaMaC matrix generators.
 */

#include "scamac_include.h"
#include "scamac_generator.h"
#include "scamac_collection.h"

#endif /* SCAMAC_H */
