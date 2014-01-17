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

const char* phist_ghost_error2str(ghost_error_t code)
  {
  if (code==GHOST_SUCCESS) return "success";
  else if (code==GHOST_ERR_INVALID_ARG) return "invalid argument";
  else if (code==GHOST_ERR_MPI) return "MPI error";
  else if (code==GHOST_ERR_CUDA) return "CUDA error";
  else if (code==GHOST_ERR_UNKNOWN) return "unknown error";
  else if (code==GHOST_ERR_INTERNAL) return "internal error";
  else if (code==GHOST_ERR_NOT_IMPLEMENTED) return "not implemented";
  else if (code==GHOST_ERR_IO) return "I/O error";
  return "unknown error code";
  }

const char* eigSort2str(eigSort_t s)
  {
  return s==LM?"LM":s==SM?"SM":s==LR?"LR":s==SR?"SR":s==NONE?"none":"unknown eigsort_t";
  }
      
