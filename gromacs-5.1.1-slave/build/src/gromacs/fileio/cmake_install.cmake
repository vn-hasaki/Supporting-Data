# Install script for directory: /home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gromacs/fileio" TYPE FILE FILES
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/confio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/enxio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/filenm.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/gmxfio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/matio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/mdoutf.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/mtxio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/pdbio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/tpxio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/trajectory_writing.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/trnio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/trx.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/trxio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/xdr_datatype.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/xtcio.h"
    "/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/src/gromacs/fileio/xvgr.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/export/base/shisuan/swyaow/online/luhao/zx/test/gromacs-5.1.1-test/build/src/gromacs/fileio/tests/cmake_install.cmake")

endif()

