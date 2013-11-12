#include "phist_tools.h"
#include "phist_macros.h"
#include "phist_enums.h"


const char* phist_retcode2str(int code)
  {
  if (code==_PHIST_SUCCESS_) return "success";
  else if (code==_PHIST_FUNCTIONAL_ERROR_) return "functional error";
  else if (code==_PHIST_CAUGHT_EXCEPTION_) return "caught exception";
  else if (code==_PHIST_BAD_CAST_) return "bad cast";
  else if (code==_PHIST_NOT_IMPLEMENTED_) return "not implemented";
  return "unknown error";
  }

const char* eigSort2str(eigSort_t s)
  {
  return s==LM?"LM":s==SM?"SM":s==LR?"LR":s==SR?"SR":s==NONE?"none":"unknown eigsort_t";
  }
      
