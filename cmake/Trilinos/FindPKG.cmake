macro(FindPKG PKG)

include(FindHOME)
FindHOME(${PKG})
if (NOT ${PKG}_HOME)
  FindHOME(Trilinos)
  set(${PKG}_HOME ${Trilinos_HOME})
endif()

if (${PKG}_HOME)
  find_package(${PKG} PATHS ${${PKG}_HOME}/lib/cmake/${PKG} ${${PKG}_HOME}/include )
endif()

endmacro()
