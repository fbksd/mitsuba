# main programs
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

include_directories(${XERCES_INCLUDE_DIRS})

#set(CMAKE_AUTOMOC ON)
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
find_package(Qt5Core)
find_package(Qt5Network)

# Qt5::Core Qt5::Network
add_mts_exe(mtssrv  mtssrv.cpp ${MTS_WINDOWS_STUB}
  RES_DESCRIPTION "Mitsuba compute node application")
add_mts_exe(mitsuba mitsuba.cpp LINK_LIBRARIES ${MTS_WINDOWS_STUB}
  RES_DESCRIPTION "Mitsuba command line interface frontend")
add_mts_exe(mtsutil mtsutil.cpp ${MTS_DARWIN_STUB} ${MTS_WINDOWS_STUB}
  RES_DESCRIPTION "Mitsuba utility launcher")
