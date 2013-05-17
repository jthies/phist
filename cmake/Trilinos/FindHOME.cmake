macro(FindHOME PKG)
  STRING(TOUPPER ${PKG} UPKG)
  set(PKG_HOME ${PKG}_HOME)
  set(UPKG_HOME ${UPKG}_HOME)
  set(ENV_PKG_HOME $ENV{${PKG_HOME}})
  set(ENV_UPKG_HOME $ENV{${UPKG_HOME}})

  message(STATUS "looking for " ${PKG_HOME})
  
  if (${PKG_HOME})
  elseif (ENV_PKG_HOME)
    set(${PKG_HOME} ${ENV_PKG_HOME})
  elseif (${UPKG_HOME})
    set(${PKG_HOME} ${${UPKG_HOME}})
  elseif (ENV_UPKG_HOME)
    set (${PKG_HOME} ${ENV_UPKG_HOME})
  endif()

  if (${PKG_HOME})
    message(STATUS "Found in " ${${PKG_HOME}})
  else()
    message(STATUS "not found.")
  endif()
  
  
endmacro()
