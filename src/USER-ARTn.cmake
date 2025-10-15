if(PKG_USER-ARTn)
  # Find MKL library
  find_package(MKL REQUIRED)
  
  set(USER-ARTn_SOURCES
    ${LAMMPS_SOURCE_DIR}/USER-ARTn/min_artn.cpp
  )

  set(USER-ARTn_HEADERS
    ${LAMMPS_SOURCE_DIR}/USER-ARTn/min_artn.h
  )

  get_property(LAMMPS_HEADERS GLOBAL PROPERTY HEADERS)
  list(APPEND LAMMPS_HEADERS ${USER-ARTn_HEADERS})
  set_property(GLOBAL PROPERTY HEADERS "${LAMMPS_HEADERS}")

  get_target_property(LAMMPS_SOURCES lammps SOURCES)
  list(APPEND LAMMPS_SOURCES ${USER-ARTn_SOURCES})
  set_property(TARGET lammps PROPERTY SOURCES "${LAMMPS_SOURCES}")
  
  # Link MKL to LAMMPS target
  target_link_libraries(lammps PRIVATE MKL::MKL)
endif()