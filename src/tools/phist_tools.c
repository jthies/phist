#include "phist_tools.h"
#include "phist_macros.h"
#include "phist_enums.h"


const char* phist_retcode2str(int code)
  {
  if (code==PHIST_SUCCESS) return "success";
  else if (code==PHIST_FUNCTIONAL_ERROR) return "functional error";
  else if (code==PHIST_CAUGHT_EXCEPTION) return "caught exception";
  else if (code==PHIST_BAD_CAST) return "bad cast or NULL pointer";
  else if (code==PHIST_NOT_IMPLEMENTED) return "not implemented";
  return "unknown error";
  }

const char* eigSort2str(eigSort_t s)
  {
  return s==LM?"LM":s==SM?"SM":s==LR?"LR":s==SR?"SR":s==NONE?"none":"unknown eigsort_t";
  }
      
