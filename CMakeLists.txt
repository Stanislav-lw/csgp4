cmake_minimum_required(VERSION 3.20)
project(csgp4 VERSION 1.0.0)
include(GNUInstallDirs)

add_library(csgp4 SHARED csgp4/afuncs.c csgp4/datetime.c csgp4/get_el.c csgp4/sgp4.c)
set_target_properties(csgp4 PROPERTIES PUBLIC_HEADER "csgp4/afuncs.h;csgp4/datetime.h;csgp4/norad.h")
target_include_directories(csgp4 PRIVATE csgp4)

install(TARGETS csgp4
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/csgp4)
