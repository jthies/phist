#include "tools.h"
#include "macros.h"


static const char* retcode2str(int code)
  {
  if (code==_ESSEX_FUNCTIONAL_ERROR_) return "functional error";
  else if (code==_ESSEX_CAUGHT_EXCEPTION_) return "caught exception";
  else if (code==_ESSEX_BAD_CAST_) return "bad cast";
  else if (code==_ESSEX_NOT_IMPLEMENTED_) return "not implemented";
  return "unknown error";
  }

