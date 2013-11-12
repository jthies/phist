#ifndef PHIST_ENUMS_H
#define PHIST_ENUMS_H

typedef enum eigSort_t {

NONE=0,
LM=1,   // largest magnitude
SM=2,   // smallest magnitude
LR=3,   // largest real part
SR=4   // smallest real part

} eigSort_t;

#ifdef __cplusplus
extern "C" {
#endif
// defined in phist_tools.c
const char* eigSort2str(eigSort_t s);
#ifdef __cplusplus
}
#endif
#endif
