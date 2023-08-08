# Install script for directory: /home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/analysisdata

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/export/base/shisuan/swyaow/online/luhao/zx/APP/gmx511")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gromacs/analysisdata" TYPE FILE FILES
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/analysisdata/abstractdata.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/analysisdata/analysisdata.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/analysisdata/arraydata.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/analysisdata/dataframe.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/analysisdata/datamodule.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/src/gromacs/analysisdata/modules/cmake_install.cmake")
  include("/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/src/gromacs/analysisdata/tests/cmake_install.cmake")

endif()

