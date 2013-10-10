if (CAPD_INCLUDE_DIR AND CAPD_LIBRARIES)
  # Already in cache, be silent
  set(CAPD_FIND_QUIETLY TRUE)
endif (CAPD_INCLUDE_DIR AND CAPD_LIBRARIES)

find_path(CAPD_INCLUDE_DIR NAMES capdlib.h REQUIRED)
find_path(CAPD_CONFIG_DIR NAMES capd-config REQUIRED)
find_library(CAPD_LIBRARIES NAMES capd capddynsys REQUIRED)
#find_library(CAPD_LIBRARIES NAMES casdapd clkjapddynsys REQUIRED)

IF ("${CAPD_INCLUDE_DIR}" MATCHES "CAPD_INCLUDE_DIR-NOTFOUND" OR
    "${CAPD_LIBRARIES}" MATCHES "CAPD_LIBRARIES-NOTFOUND")
  MESSAGE(FATAL_ERROR "Could not find CAPD-DynSys 3.0. Please visit http://capd.ii.uj.edu.pl/download.php to download")
ELSE ()
  include(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(CAPD DEFAULT_MSG CAPD_INCLUDE_DIR CAPD_LIBRARIES)
  mark_as_advanced(CAPD_INCLUDE_DIR CAPD_LIBRARIES)

  execute_process(COMMAND ${CAPD_CONFIG_DIR}/capd-config "--libs"
    OUTPUT_VARIABLE CAPD_LD_FLAGS)
  string(STRIP ${CAPD_LD_FLAGS} CAPD_LD_FLAGS)

  execute_process(COMMAND ${CAPD_CONFIG_DIR}/capd-config "--cflags"
    OUTPUT_VARIABLE CAPD_CFLAGS)
  string(STRIP ${CAPD_CFLAGS} CAPD_CFLAGS)
ENDIF ()
