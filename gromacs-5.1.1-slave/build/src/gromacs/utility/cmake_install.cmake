# Install script for directory: /home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gromacs/utility" TYPE FILE FILES
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/arrayref.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/basedefinitions.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/classhelpers.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/cstringutil.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/datafilefinder.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/errorcodes.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/exceptions.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/fatalerror.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/file.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/flags.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/futil.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/gmxassert.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/init.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/programcontext.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/real.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/smalloc.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/utility/stringutil.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/src/gromacs/utility/tests/cmake_install.cmake")

endif()

