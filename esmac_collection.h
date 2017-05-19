#ifndef ESMAC_COLLECTION_H
#define ESMAC_COLLECTION_H

#include <stdbool.h>

#include "esmac_generator.h"

esmac_generator_t * esmac_collection_example(const char * name, int * info);
void esmac_collection_example_free(esmac_generator_t * gen);

int esmac_collection_set_int(esmac_generator_t * gen, const char * parname, int val);
int esmac_collection_set_double(esmac_generator_t * gen, const char * parname, double val);
int esmac_collection_set_bool(esmac_generator_t * gen, const char * parname, bool val);
int esmac_collection_set_option(esmac_generator_t * gen, const char * parname, char * option);

char * esmac_collection_list_examples();
char * esmac_collection_list_parameters(const char * name, int * info);

int esmac_collection_identify_parameter(const char * name, const char * parname);

char * esmac_collection_example_data(const esmac_generator_t * gen);

esmac_generator_t * esmac_collection_parse(const char * argstring, int * info);

#endif /* ESMAC_COLLECTION_H */
