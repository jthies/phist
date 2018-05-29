/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  high-level access to the ScaMaC
 *  \ingroup library
 */

#ifndef SCAMAC_COLLECTION_H
#define SCAMAC_COLLECTION_H

#include "scamac_include.h"

ScamacErrorCode scamac_generator_obtain(const char * matname, ScamacGenerator ** gen);

ScamacErrorCode scamac_parse_argstr(const char * argstr, ScamacGenerator ** gen, char ** errdesc);

ScamacErrorCode scamac_generator_set_int     (ScamacGenerator * gen, const char * parname, int val);
ScamacErrorCode scamac_generator_set_double  (ScamacGenerator * gen, const char * parname, double val);
ScamacErrorCode scamac_generator_set_bool    (ScamacGenerator * gen, const char * parname, bool val);
ScamacErrorCode scamac_generator_set_rngseed (ScamacGenerator * gen, const char * parname, uint64_t seed);
ScamacErrorCode scamac_generator_set_rngseed_str(ScamacGenerator * gen, const char * parname, const char * seedstr);
ScamacErrorCode scamac_generator_set_option  (ScamacGenerator * gen, const char * parname, const char * option);

/** textual representation of generator parameters.
 * Returns a string that, when passed to scamac_parse_argstr, produces an identical generator.
 * Useful, e.g., to pass around a generator between different processes
 *
 * format="argstr","desc"
 */
ScamacErrorCode scamac_generator_parameter_desc(const ScamacGenerator * gen, const char * format, char ** desc);


ScamacErrorCode scamac_list_examples(char ** desc);
ScamacErrorCode scamac_example_parameters(const char * matname, char ** desc);
ScamacErrorCode scamac_example_desc(const char * matname, int * valtype, int * symmetry, char ** desc);



ScamacPar scamac_identify_parameter(const char * matname, const char * parname);


/** \brief Return error name.
 */
const char * scamac_error_desc(ScamacErrorCode err);
         int scamac_error_par (ScamacErrorCode err);
const char * scamac_error_dpar(ScamacErrorCode err);

#endif /* SCAMAC_COLLECTION_H */
